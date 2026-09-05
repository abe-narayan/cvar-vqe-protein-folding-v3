"""Sequence-conditioned per-residue torsion distribution: a mixture of von Mises over
(phi, psi).

WHY THIS EXISTS. Sprint 5 established two things about the current system. The candidate
pool is not the problem -- a pool of real protein fragments retrieved by sequence
similarity reaches 1.356 A pool-best on held-out targets against Sprint 4's 2.355 A. And
selection is the entire problem: the learned distogram returns 3.3-3.5 A from those pools
regardless of whether the pool's best member is at 1.79 A or 1.36 A. It is insensitive to
the quality of what it is ranking.

The distogram is the only structural signal the system has, and it is a CA-CA distance
signal. The published state of the art for this length class (APPTEST, Brief Bioinform
2021) predicts distances AND backbone torsions, and reaches 1.57 A best-of-100 on 6-12mers
from a third of our training data. Torsions are the channel this project has never used.

They are a genuinely independent signal, not a reparameterisation of the same one:
distances are a global, pairwise, O(n^2) description that says nothing directly about
chirality or local backbone conformation, while (phi, psi) are local, O(n), and fix exactly
the thing distances are blind to. A candidate structure can match a predicted distance
matrix closely while sitting in the wrong region of Ramachandran space.

WHAT THIS MODEL IS. Per residue, a K-component mixture of independent von Mises densities
over phi and psi, conditioned on the sequence through a dilated convolution stack. A
mixture, not a point prediction, because backbone torsions are genuinely multimodal --
alpha and beta basins are both accessible to most residues and a unimodal fit lands between
them, in a region no residue occupies. That failure mode is the torsion-space analogue of
the lumpy-distogram problem already recorded in this project.

TWO USES, AND ONE OF THEM IS SAFER THAN THE OTHER.

  * SCORING (`score_torsions`): the log-likelihood of an existing structure's torsions.
    This is a drop-in selection signal over any candidate pool and is the use that
    addresses the measured bottleneck. It cannot produce invalid geometry because it never
    produces geometry.
  * SAMPLING (`sample`): draw torsions and build a structure. Per-residue independent draws
    ignore correlations along the chain, so sampled structures accumulate lever-arm error
    and can self-intersect. At length 9-15 the lever arm is short, but this is the riskier
    use and is evaluated separately rather than assumed to work.

Leakage: `train_fold(f)` trains on peptides whose identity fold is not `f`, plus
`distogram._fold_fragments(f)`, which is exactly the production discipline. Nothing here
reads a native at inference.
"""
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import torch
import torch.nn as nn

import distogram as dgm
import peptide_db as pdb
from s5 import esm32

AA = "ARNDCQEGHILKMFPSTWYV"
AA_IDX = {a: i for i, a in enumerate(AA)}
MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "torsion_models")
D_ESM = 32
N_COMP = 8


# --------------------------------------------------------------------------- features
def features(seq: str, esm: Optional[np.ndarray] = None) -> np.ndarray:
    """``(n, D)`` per-residue features. One-hot identity, ESM projection, position.

    Position is encoded as three quantities rather than one index because terminal
    residues behave differently from interior ones and the chain is short enough that
    "how far from each end" is a large fraction of the available context.
    """
    n = len(seq)
    hot = np.zeros((n, 20), dtype=np.float32)
    for i, c in enumerate(seq):
        j = AA_IDX.get(c)
        if j is not None:
            hot[i, j] = 1.0
    if esm is None:
        esm = esm32.load().get(seq)
        if esm is None:
            import esm_features as ef
            esm = ef.embed(seq).astype(np.float32)
    esm = np.asarray(esm, dtype=np.float32)[:n]
    if len(esm) < n:                                    # defensive: pad rather than crash
        esm = np.vstack([esm, np.zeros((n - len(esm), esm.shape[1]), np.float32)])
    idx = np.arange(n, dtype=np.float32)
    pos = np.stack([idx / max(n - 1, 1),
                    (n - 1 - idx) / max(n - 1, 1),
                    np.full(n, 1.0 / n, dtype=np.float32)], 1)
    return np.concatenate([hot, esm, pos], 1).astype(np.float32)


D_IN = 20 + D_ESM + 3


# ----------------------------------------------------------------------------- model
class TorsionNet(nn.Module):
    """Dilated 1-D convolutions -> per-residue mixture of von Mises over (phi, psi).

    Dilated convolutions rather than attention: at length 9-15 a dilation stack of
    1/2/4 with width 5 already has a receptive field covering the whole chain, and it has
    far fewer parameters to fit from ~82,000 supervised residues than a transformer does.
    """

    def __init__(self, d_in: int = D_IN, width: int = 128, n_comp: int = N_COMP,
                 dropout: float = 0.1):
        super().__init__()
        self.n_comp = n_comp
        self.stem = nn.Conv1d(d_in, width, 5, padding=2)
        # padding is 2*d, not d: a width-5 kernel at dilation d spans 4*d+1, so preserving
        # length needs (5-1)*d/2 on each side. Getting this wrong silently shortens the
        # sequence and the residual add then fails on a shape mismatch.
        self.blocks = nn.ModuleList([
            nn.Conv1d(width, width, 5, padding=2 * d, dilation=d) for d in (1, 2, 4)])
        self.norms = nn.ModuleList([nn.GroupNorm(8, width) for _ in range(3)])
        self.drop = nn.Dropout(dropout)
        # per component: mu_phi(cos,sin), mu_psi(cos,sin), log kappa_phi, log kappa_psi, logit
        self.head = nn.Conv1d(width, n_comp * 7, 1)

    def forward(self, x: torch.Tensor) -> Dict[str, torch.Tensor]:
        """`x` is (B, n, D_IN). Returns mixture parameters, each (B, n, K)."""
        h = torch.nn.functional.gelu(self.stem(x.transpose(1, 2)))
        for conv, norm in zip(self.blocks, self.norms):
            h = h + self.drop(torch.nn.functional.gelu(norm(conv(h))))
        o = self.head(h).transpose(1, 2)                       # (B, n, K*7)
        B, n, _ = o.shape
        o = o.reshape(B, n, self.n_comp, 7)
        return {
            "mu_phi": torch.atan2(o[..., 1], o[..., 0]),
            "mu_psi": torch.atan2(o[..., 3], o[..., 2]),
            # kappa in [~0.05, ~150]: unbounded concentration lets one component collapse
            # onto a single training example and take all the mixture weight.
            "kappa_phi": torch.exp(o[..., 4].clamp(-3.0, 5.0)),
            "kappa_psi": torch.exp(o[..., 5].clamp(-3.0, 5.0)),
            "logit": o[..., 6],
        }


# ------------------------------------------------------------------------ likelihood
def _log_i0(k: torch.Tensor) -> torch.Tensor:
    """log I0(k), numerically stable for large k via the exponentially scaled Bessel."""
    return torch.log(torch.special.i0e(k)) + k


def _vm_logpdf(x: torch.Tensor, mu: torch.Tensor, kappa: torch.Tensor) -> torch.Tensor:
    return kappa * torch.cos(x - mu) - math.log(2.0 * math.pi) - _log_i0(kappa)


def mixture_logpdf(p: Dict[str, torch.Tensor], phi: torch.Tensor,
                   psi: torch.Tensor) -> torch.Tensor:
    """``(B, n)`` log-density of each residue's (phi, psi) under its mixture."""
    lw = torch.log_softmax(p["logit"], dim=-1)
    lp = (_vm_logpdf(phi[..., None], p["mu_phi"], p["kappa_phi"])
          + _vm_logpdf(psi[..., None], p["mu_psi"], p["kappa_psi"]))
    return torch.logsumexp(lw + lp, dim=-1)


# -------------------------------------------------------------------------- sampling
def _sample_vm(mu: np.ndarray, kappa: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Best-Fisher rejection sampler for the von Mises distribution.

    Wrapped-normal approximation is tempting and wrong at low kappa, which is exactly
    where a torsion model is uncertain and where the sampled diversity has to come from.
    """
    mu = np.asarray(mu, float)
    kappa = np.maximum(np.asarray(kappa, float), 1e-8)
    out = np.empty(mu.shape, float)
    flat_mu, flat_k, flat_out = mu.ravel(), kappa.ravel(), out.ravel()
    a = 1.0 + np.sqrt(1.0 + 4.0 * flat_k ** 2)
    b = (a - np.sqrt(2.0 * a)) / (2.0 * flat_k)
    r = (1.0 + b ** 2) / (2.0 * b)
    todo = np.arange(flat_mu.size)
    for _ in range(100):
        if not todo.size:
            break
        u1, u2, u3 = rng.random(todo.size), rng.random(todo.size), rng.random(todo.size)
        z = np.cos(np.pi * u1)
        f = (1.0 + r[todo] * z) / (r[todo] + z)
        c = flat_k[todo] * (r[todo] - f)
        ok = (c * (2.0 - c) - u2 > 0) | (np.log(c / np.maximum(u2, 1e-300)) + 1.0 - c >= 0)
        acc = todo[ok]
        flat_out[acc] = flat_mu[acc] + np.sign(u3[ok] - 0.5) * np.arccos(
            np.clip(f[ok], -1.0, 1.0))
        todo = todo[~ok]
    if todo.size:                                    # give up gracefully, never hang
        flat_out[todo] = flat_mu[todo]
    return np.mod(out + np.pi, 2.0 * np.pi) - np.pi


# ------------------------------------------------------------------------- interface
class TorsionPrior:
    """Trained model for one sequence: score torsions, or sample them."""

    def __init__(self, seq: str, net: TorsionNet):
        self.seq = seq
        self.net = net.eval()
        with torch.no_grad():
            x = torch.from_numpy(features(seq))[None]
            self.p = {k: v[0] for k, v in net(x).items()}

    def score_torsions(self, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
        """Mean per-residue negative log-likelihood. Lower is better. Batched over rows.

        Mean rather than sum so the number is comparable across chain lengths, which
        matters the moment it is combined with any other term.
        """
        P = torch.as_tensor(np.atleast_2d(phi), dtype=torch.float32)
        S = torch.as_tensor(np.atleast_2d(psi), dtype=torch.float32)
        B = P.shape[0]
        p = {k: v[None].expand(B, *v.shape) for k, v in self.p.items()}
        with torch.no_grad():
            ll = mixture_logpdf(p, P, S)
        return (-ll.mean(1)).numpy()

    def sample(self, n: int, seed: int = 0, temperature: float = 1.0
               ) -> Tuple[np.ndarray, np.ndarray]:
        """`n` independent per-residue draws. Temperature scales kappa: >1 sharpens.

        Per-residue independence is the known weakness -- see the module docstring.
        """
        rng = np.random.default_rng(seed)
        w = torch.softmax(self.p["logit"], -1).numpy()            # (L, K)
        L, K = w.shape
        comp = np.stack([rng.choice(K, size=n, p=w[i] / w[i].sum()) for i in range(L)], 1)
        take = lambda t: np.asarray(t)[np.arange(L)[None, :], comp]
        kp = take(self.p["kappa_phi"]) * temperature
        ks = take(self.p["kappa_psi"]) * temperature
        phi = _sample_vm(take(self.p["mu_phi"]), kp, rng)
        psi = _sample_vm(take(self.p["mu_psi"]), ks, rng)
        return phi, psi


# -------------------------------------------------------------------------- training
def _dataset(entries: Sequence[pdb.Peptide]) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                                      np.ndarray]:
    """Right-pad to the longest chain; a mask keeps padding out of the loss."""
    L = max(len(e.seq) for e in entries)
    X = np.zeros((len(entries), L, D_IN), np.float32)
    PHI = np.zeros((len(entries), L), np.float32)
    PSI = np.zeros((len(entries), L), np.float32)
    M = np.zeros((len(entries), L), np.float32)
    cache = esm32.load()
    for i, e in enumerate(entries):
        n = len(e.seq)
        X[i, :n] = features(e.seq, cache.get(e.seq))
        PHI[i, :n] = e.phi[:n]
        PSI[i, :n] = e.psi[:n]
        M[i, :n] = 1.0
    ok = np.isfinite(PHI) & np.isfinite(PSI)
    return X, np.nan_to_num(PHI), np.nan_to_num(PSI), M * ok


def _model_path(fold: int, seed: int) -> str:
    return os.path.join(MODEL_DIR, f"torsion_fold{fold}_s{seed}.pt")


def train_fold(fold: int, n_folds: int = 5, epochs: int = 300, seed: int = 0,
               lr: float = 2e-3, batch: int = 256, width: int = 128,
               verbose: bool = True, force: bool = False) -> TorsionNet:
    """Train (or load) the model that excludes identity fold `fold`."""
    os.makedirs(MODEL_DIR, exist_ok=True)
    path = _model_path(fold, seed)
    net = TorsionNet(width=width)
    if os.path.exists(path) and not force:
        net.load_state_dict(torch.load(path, map_location="cpu"))
        return net.eval()
    folds = pdb.folds(n_folds)
    entries = [p for p in pdb.load() if folds[p.seq] != fold]
    entries = entries + list(dgm._fold_fragments(fold, n_folds))
    X, PHI, PSI, M = _dataset(entries)
    if verbose:
        print(f"  torsion fold {fold}: {len(entries)} chains, "
              f"{int(M.sum())} supervised residues, {sum(p.numel() for p in net.parameters())} params",
              flush=True)
    Xt = torch.from_numpy(X); Pt = torch.from_numpy(PHI)
    St = torch.from_numpy(PSI); Mt = torch.from_numpy(M)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    g = torch.Generator().manual_seed(seed)
    net.train()
    for ep in range(epochs):
        perm = torch.randperm(len(Xt), generator=g)
        tot, cnt = 0.0, 0
        for s in range(0, len(perm), batch):
            idx = perm[s:s + batch]
            p = net(Xt[idx])
            ll = mixture_logpdf(p, Pt[idx], St[idx])
            m = Mt[idx]
            loss = -(ll * m).sum() / m.sum().clamp(min=1)
            opt.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(net.parameters(), 5.0)
            opt.step()
            tot += float(loss) * float(m.sum()); cnt += float(m.sum())
        sched.step()
        if verbose and (ep % 25 == 0 or ep == epochs - 1):
            print(f"    epoch {ep:4d}  NLL/residue {tot/max(cnt,1):.4f}", flush=True)
    tmp = path + ".tmp"
    torch.save(net.state_dict(), tmp)
    os.replace(tmp, path)
    return net.eval()


_PRIORS: Dict[Tuple[str, int], TorsionPrior] = {}


def for_target(seq: str, n_folds: int = 5, seed: int = 0, **kw) -> TorsionPrior:
    """Leave-fold-out prior for `seq`, matching the production discipline."""
    fold = pdb.folds(n_folds).get(seq, 0)
    key = (seq, seed)
    if key not in _PRIORS:
        net = train_fold(fold, n_folds=n_folds, seed=seed, verbose=False, **kw)
        _PRIORS[key] = TorsionPrior(seq, net)
    return _PRIORS[key]
