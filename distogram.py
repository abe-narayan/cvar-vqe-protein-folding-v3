"""Sequence -> CA-CA distance distribution, and the structure score built on it.

Why this is the centre of the redesign. On decoy populations drawn from the search's own
representation, scoring a structure by its agreement with the *true* distance matrix gives
Spearman +0.85 to +0.96 against CA-RMSD and picks a top-20 averaging 0.4-2.6 A -- against
+0.10 for the 16-descriptor learned scorer, ~0.0 for the knowledge-based energy and -0.40
for Amber. A distance matrix is therefore a sufficient statistic for ranking here, and the
whole ranking problem reduces to predicting one from sequence.

Scoring form. The predicted distribution over bins is used through its L1 Bayes risk,

    s(structure) = mean_{ij} w_ij * sum_b P(b | i, j) * | d_ij - c_b | ,

not through the log-likelihood of the occupied bin. The two were measured side by side:
the log-likelihood form scores +0.02 to -0.42 (it is dominated by the entropy of the
predicted distribution and is blind to *how far* a wrong distance is wrong), while the
Bayes-risk form scores +0.16 to +0.82 on the same predictions. ``w_ij`` is the inverse
predicted spread, so pairs the model is sure about count more.

Model. A small residual MLP over pair features. Trees were tried first and are what the
sequence-only features get; they cannot use a 1280-dimensional ESM-2 embedding, and the
embedding is where the remaining signal is. Training is leave-fold-out over identity
clusters (`peptide_db.folds`), so a benchmark target is always scored by a model that saw
neither it nor any sequence above 0.6 identity to it.
"""
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import peptide_db as pdb
import priors

BASE = os.path.dirname(os.path.abspath(__file__))
import fragment_db as _fdb
MODEL_DIR = os.path.join(
    BASE, "distogram_models_large" if _fdb.USE_LARGE else "distogram_models")

BIN_EDGES = priors.BIN_EDGES
NBINS = priors.NBINS
CENTRES = np.concatenate([[BIN_EDGES[0] - 0.5],
                          0.5 * (BIN_EDGES[1:] + BIN_EDGES[:-1]),
                          [BIN_EDGES[-1] + 2.0]])


# --------------------------------------------------------------------- features
def _esm_available() -> bool:
    try:
        import esm_features
        return esm_features.available() and os.path.exists(esm_features._PCA_PATH)
    except Exception:
        return False


def features(seq: str, use_esm: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(X (npairs, D), i, j)`` for one sequence."""
    base, i, j = priors.pair_features(seq)
    if not use_esm or not _esm_available():
        return base, i, j
    import esm_features
    emb = esm_features.embed(seq)                     # (n, N_PCA)
    con = esm_features.contacts(seq)                  # (n, n)
    n = len(seq)
    # Row statistics of the contact map: how "central" each residue is predicted to be.
    deg = con.sum(1, keepdims=True) / max(1, n - 1)
    # Contact mass in the sequence-separation shells that decide fold class.
    sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    shells = np.stack([np.where((sep >= a) & (sep < b), con, 0.0).sum(1) / max(1, n)
                       for a, b in ((3, 5), (5, 9), (9, 100))], axis=1)
    extra = np.column_stack([
        con[i, j],
        # local neighbourhood of the predicted contact, which is what a real ladder or
        # helical turn looks like: contacts come in diagonal blocks, not isolated cells.
        np.stack([con[np.clip(i + a, 0, n - 1), np.clip(j + b, 0, n - 1)]
                  for a, b in ((-1, -1), (1, 1), (-1, 1), (1, -1))], axis=1),
        deg[i], deg[j], shells[i], shells[j],
        emb[i], emb[j], emb[i] * emb[j], np.abs(emb[i] - emb[j]),
    ])
    return np.hstack([base, extra]).astype(np.float32), i, j


def dataset(entries: Sequence[pdb.Peptide], use_esm: bool = True,
            with_sep: bool = False):
    X, Y, S = [], [], []
    for p in entries:
        f, i, j = features(p.seq, use_esm)
        if not len(f):
            continue
        d = np.linalg.norm(p.ca[i] - p.ca[j], axis=1)
        X.append(f)
        Y.append(np.digitize(d, BIN_EDGES))
        if with_sep:
            S.append((j - i).astype(np.float32))
    if with_sep:
        return np.vstack(X), np.concatenate(Y), np.concatenate(S)
    return np.vstack(X), np.concatenate(Y)


def separation_weights(sep: np.ndarray, mode: str = "lin", cap: float = 6.0
                       ) -> np.ndarray:
    """Training-loss weights as a function of sequence separation."""
    s = np.asarray(sep, np.float32)
    if mode == "none":
        return np.ones_like(s)
    if mode == "lin":
        return np.minimum(s / 3.0, cap)
    if mode == "sq":
        return np.minimum((s / 3.0) ** 2, cap)
    if mode == "step":
        return np.where(s >= 6, 4.0, 1.0).astype(np.float32)
    if mode == "longonly":
        return np.where(s >= 6, 1.0, 0.15).astype(np.float32)
    raise ValueError(mode)


# --------------------------------------------------------------------- model
class MLP:
    """Residual MLP over pair features, trained with a soft-binned cross-entropy.

    The target is smoothed across neighbouring bins (a 4.9 A pair is mostly bin 0 but
    partly bin 1). Hard one-hot targets over 17 narrow bins make the model confident about
    a boundary it cannot resolve, and the Bayes-risk score then reads that confidence as
    information.
    """

    def __init__(self, d_in: int, width: int = 384, depth: int = 3, seed: int = 0,
                 dropout: float = 0.2):
        import torch
        import torch.nn as nn
        torch.manual_seed(seed)
        self.torch = torch
        # The Dropout modules are omitted entirely at p = 0 so that checkpoints trained
        # before dropout existed still load: inserting an identity-behaving module would
        # still shift every layer index in the state dict.
        drop = ([nn.Dropout(dropout)] if dropout > 0 else [])
        layers = [nn.Linear(d_in, width), nn.GELU()] + list(drop)
        for _ in range(depth - 1):
            layers += [nn.Linear(width, width), nn.GELU()] + list(drop)
        layers += [nn.Linear(width, NBINS)]
        self.net = nn.Sequential(*layers)
        self.mu = np.zeros(d_in, np.float32)
        self.sd = np.ones(d_in, np.float32)

    def _norm(self, X):
        return (X - self.mu) / self.sd

    def fit(self, X: np.ndarray, Y: np.ndarray, epochs: int = 40, batch: int = 4096,
            lr: float = 2e-3, smooth: float = 0.6, sample_w: Optional[np.ndarray] = None,
            verbose: bool = False):
        t = self.torch
        self.mu = X.mean(0)
        self.sd = X.std(0) + 1e-3
        Xn = t.as_tensor(self._norm(X))
        # soft targets
        idx = np.arange(NBINS)[None, :]
        T = np.exp(-((idx - Y[:, None]) ** 2) / (2 * smooth ** 2)).astype(np.float32)
        T /= T.sum(1, keepdims=True)
        T = t.as_tensor(T)
        # Per-example loss weights. Short-range pairs are both more numerous and much
        # easier -- the diagnosed failure is entirely in the long-range band, where the
        # prior's MAE is 1.5-2.2 A against a near-native spread of 0.36-0.48 A -- so an
        # unweighted loss spends capacity where the model is already right.
        W = (t.ones(len(X)) if sample_w is None
             else t.as_tensor(np.asarray(sample_w, np.float32)))
        W = W / W.mean()
        opt = t.optim.AdamW(self.net.parameters(), lr=lr, weight_decay=1e-4)
        sched = t.optim.lr_scheduler.OneCycleLR(
            opt, max_lr=lr, total_steps=epochs * max(1, len(X) // batch + 1))
        n = len(X)
        g = t.Generator().manual_seed(0)
        self.net.train()
        for ep in range(epochs):
            perm = t.randperm(n, generator=g)
            tot = 0.0
            for s in range(0, n, batch):
                b = perm[s:s + batch]
                logit = self.net(Xn[b])
                loss = (-(T[b] * t.log_softmax(logit, 1)).sum(1) * W[b]).mean()
                opt.zero_grad()
                loss.backward()
                opt.step()
                sched.step()
                tot += float(loss.detach()) * len(b)
            if verbose and ep % 10 == 0:
                print(f"    epoch {ep} loss {tot / n:.4f}", flush=True)
        return self

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        t = self.torch
        self.net.eval()
        with t.no_grad():
            return t.softmax(self.net(t.as_tensor(self._norm(X))), 1).numpy()

    def save(self, path: str):
        self.torch.save({"sd": self.net.state_dict(), "mu": self.mu, "s": self.sd}, path)

    def load(self, path: str):
        z = self.torch.load(path, weights_only=False)
        self.net.load_state_dict(z["sd"])
        self.mu, self.sd = z["mu"], z["s"]
        return self

#: Regularisation defaults -- and the reason they are OFF is the most counter-intuitive
#: measurement in this module, so it is recorded in full.
#:
#: The model is firmly OVERFIT on a held-out slice of one fold's own training pairs
#: (`work/reg_sweep.py`): at width 384 x depth 3 with no dropout, train NLL 0.944 against
#: validation NLL 2.025, and 90 epochs instead of 40 makes validation *worse* (2.526).
#: Dropout 0.2 fixes that -- validation NLL 1.688, expected-distance MAE 1.82 -> 1.76 A --
#: and a 3-model ensemble improves it again to 1.645 / 1.72 A. On the benchmark the
#: regularised ensemble's distance matrix is decisively better: correlation with the true
#: distances 0.744 -> 0.805, MAE 2.04 -> 1.79 A.
#:
#: And its RANKING is worse: in-band Spearman 0.425 -> 0.363, top-20 2.75 -> 2.78 A.
#: The mechanism is the score's weighting. `Distogram` weights each pair by the inverse of
#: its predicted spread, and dropout plus averaging widen the predicted distributions --
#: they are better calibrated and less overconfident, which flattens the weights and blurs
#: the Bayes-risk curve. Sharpening the averaged distribution back (`work/sharpen.py`,
#: p -> p^T for T up to 4) recovers part of it, 0.363 -> 0.388, and never reaches 0.425.
#:
#: So ranking here wants *sharp and discriminative* per-pair predictions, not unbiased
#: ones, which is the same lesson as the distance-geometry projection in section 13 of the
#: report. Production keeps the single sharp model; the regularised path stays available
#: because a downstream use that needs a calibrated *distance* rather than a ranking should
#: prefer it.
DROPOUT = 0.0
N_ENSEMBLE = 1


def _model_path(fold: int, use_esm: bool, frag: bool, seed: int = 0) -> str:
    tag = ("esm" if use_esm else "seq") + ("_frag" if frag else "")
    suffix = "" if seed == 0 else f"_s{seed}"
    return os.path.join(MODEL_DIR, f"fold{fold}_{tag}{suffix}.pt")


def _fold_fragments(fold: int, n_folds: int, threshold: float = pdb.IDENTITY_THRESHOLD):
    """Fragments safe to train on for `fold`: none above the identity threshold to any
    peptide the fold holds out."""
    import fragment_db as fdb
    held = [p.seq for p in pdb.load() if pdb.folds(n_folds)[p.seq] == fold]
    hk = [pdb._kmers(s) for s in held]
    out = []
    for f in fdb.load():
        fk = pdb._kmers(f.seq)
        leak = False
        for s, k in zip(held, hk):
            if fk and k and not (fk & k):
                continue
            if pdb.identity(s, f.seq) >= threshold:
                leak = True
                break
        if not leak:
            out.append(f)
    return out


def train_fold(fold: int, use_esm: bool = True, n_folds: int = 5,
               epochs: int = 40, fragments: bool = False, dropout: float = 0.0,
               n_models: int = 1, verbose: bool = True) -> List[MLP]:
    """Train (or load) the ensemble that excludes identity fold `fold`."""
    os.makedirs(MODEL_DIR, exist_ok=True)
    entries = [p for p in pdb.load() if pdb.folds(n_folds)[p.seq] != fold]
    d_in = features(entries[0].seq, use_esm)[0].shape[1]
    tag = f"d{int(dropout * 100)}" if dropout else ""
    out, missing = [], []
    for s in range(n_models):
        path = _model_path(fold, use_esm, fragments, s) + (f".{tag}" if tag else "")
        m = MLP(d_in, seed=s, dropout=dropout)
        if os.path.exists(path):
            out.append(m.load(path))
        else:
            missing.append((m, path))
    if not missing:
        return out
    if fragments:
        entries = entries + _fold_fragments(fold, n_folds)
    X, Y = dataset(entries, use_esm)
    if verbose:
        print(f"  distogram fold {fold} ({'esm' if use_esm else 'seq'}"
              f"{'+frag' if fragments else ''}, dropout {dropout}, "
              f"{len(missing)} model(s)): {len(entries)} chains, {len(X)} pairs, "
              f"{d_in} features", flush=True)
    for m, path in missing:
        m.fit(X, Y, epochs=epochs, verbose=verbose)
        m.save(path)
        out.append(m)
    return out

# --------------------------------------------------------------------- prior
#: Per-separation-shell weights and the exponent on the predicted spread, fitted on the
#: 24-peptide dev set by maximising in-band Spearman (`work/fit_weights.py`). The default
#: -- every shell equal, exponent 1 -- was a guess; the fit says long-range pairs
#: (separation >= 14) are worth about 2.2x a short-range one and the middle shells about
#: half. RETIRED: that fit predated the fold pinning and the dev-set cluster-disjointness
#: fix, and does not reproduce under corrected conditions -- refitting now drives two of the
#: six parameters to their clip bounds and loses on the benchmark (top-1 2.76 -> 2.91 A)
#: while gaining 65% on the dev objective. `score_weights.json` is intentionally absent, so
#: `_score_weights()` returns unit weights; see work/benchmark_shift.md.
SHELLS = [(2, 3), (4, 5), (6, 8), (9, 13), (14, 40)]
_WPATH = os.path.join(BASE, "score_weights.json")


def _score_weights():
    if os.path.exists(_WPATH):
        import json
        z = json.load(open(_WPATH))
        return np.asarray(z["ws"], float), float(z["gamma"])
    return np.ones(len(SHELLS)), 1.0


class Distogram:
    """Per-target predicted distance distribution and the score built on it."""

    def __init__(self, sequence: str, prob: np.ndarray, i: np.ndarray, j: np.ndarray,
                 shell_weights: Optional[np.ndarray] = None,
                 gamma: Optional[float] = None):
        self.seq = sequence
        self.n = len(sequence)
        self.i, self.j = i, j
        self.prob = prob                                     # (npairs, NBINS)
        self.expected = (prob * CENTRES).sum(1)
        var = (prob * (CENTRES[None] - self.expected[:, None]) ** 2).sum(1)
        self.sd = np.sqrt(np.maximum(var, 1e-6))
        fws, fg = _score_weights()
        ws = fws if shell_weights is None else np.asarray(shell_weights, float)
        g = fg if gamma is None else float(gamma)
        sep = (j - i).astype(float)
        shell = np.zeros(len(sep))
        for k, (a, b) in enumerate(SHELLS):
            shell[(sep >= a) & (sep <= b)] = ws[k]
        self.w = shell / np.maximum((self.sd + 0.5) ** g, 1e-6)
        self.w = self.w / max(self.w.mean(), 1e-12)
        # Bayes-risk lookup on a fine distance grid, so scoring a batch is a gather.
        self.grid = np.arange(2.0, 40.0, 0.05)
        self._risk = (prob[:, None, :] * np.abs(self.grid[None, :, None]
                                                - CENTRES[None, None, :])).sum(2)
        self._risk = (self._risk * self.w[:, None]).astype(np.float32)

    @classmethod
    def for_target(cls, sequence: str, model=None, fold: Optional[int] = None,
                   use_esm: bool = True, n_folds: int = 5,
                   fragments: bool = True, dropout: float = DROPOUT,
                   n_models: int = N_ENSEMBLE) -> "Distogram":
        if model is None:
            if fold is None:
                fold = pdb.folds(n_folds).get(sequence, 0)
            model = train_fold(fold, use_esm, n_folds, fragments=fragments,
                               dropout=dropout, n_models=n_models, verbose=False)
        models = list(model) if isinstance(model, (list, tuple)) else [model]
        X, i, j = features(sequence, use_esm)
        prob = np.mean([m.predict_proba(X) for m in models], axis=0)
        return cls(sequence, prob, i, j)

    def score(self, ca: np.ndarray) -> np.ndarray:
        """Mean weighted Bayes-risk distance error for a batch ``(B, n, 3)``."""
        arr = np.asarray(ca, float)
        single = arr.ndim == 2
        if single:
            arr = arr[None]
        d = np.linalg.norm(arr[:, self.i, :] - arr[:, self.j, :], axis=-1)
        g = np.clip(((d - self.grid[0]) / 0.05).astype(np.int32), 0, len(self.grid) - 1)
        out = np.take_along_axis(self._risk.T, g, axis=0).mean(1) if False else \
            self._risk[np.arange(len(self.i))[None, :], g].mean(1)
        return float(out[0]) if single else out

    def bagged_rank(self, ca: np.ndarray, n_bags: int = 12, frac: float = 0.6,
                    seed: int = 0) -> np.ndarray:
        """Mean rank of each structure under `n_bags` random pair subsets.

        One scalar score over all pairs hides which pairs it is built from, and the model
        is wrong about some of them. Re-scoring on random 60% subsets and averaging the
        RANK (not the score) makes the ordering depend on agreement across subsets rather
        than on a few confidently-wrong pairs, which is the failure mode that puts a
        4 A structure at the top of the pool on the targets where the prior is off.
        """
        from scipy.stats import rankdata
        arr = np.asarray(ca, float)
        if arr.ndim == 2:
            arr = arr[None]
        d = np.linalg.norm(arr[:, self.i, :] - arr[:, self.j, :], axis=-1)
        g = np.clip(((d - self.grid[0]) / 0.05).astype(np.int32), 0,
                    len(self.grid) - 1)
        per_pair = self._risk[np.arange(len(self.i))[None, :], g]     # (B, npairs)
        rng = np.random.default_rng(seed)
        k = max(4, int(frac * per_pair.shape[1]))
        acc = np.zeros(len(arr))
        for _ in range(n_bags):
            cols = rng.choice(per_pair.shape[1], size=k, replace=False)
            acc += rankdata(per_pair[:, cols].mean(1))
        return acc / n_bags

    def realize(self, seed: int = 0, restarts: int = 4, iters: int = 400,
                lr: float = 0.06, chain: float = 3.81, k_chain: float = 6.0
                ) -> np.ndarray:
        """``(n, 3)`` CA coordinates fitted to the predicted distances.

        Weighted stress majorisation with a CA-CA chain restraint, best of `restarts`.
        This is the distogram's own answer to "what structure is this?", with no search and
        no energy model. It is NOT used for ranking -- projecting the matrix and scoring
        against the projection measurably hurts (see the report) -- it is used as a SEED:
        the torsion assignment closest to this trace is where the VQE starts looking.
        """
        n = self.n
        E, W = self.expected, self.w
        i, j = self.i, self.j
        adj = np.arange(n - 1)
        rng = np.random.default_rng(seed)
        best, best_s = None, np.inf
        for r in range(restarts):
            X = rng.normal(0, 4.0, (n, 3))
            X[:, 0] += np.arange(n) * 2.0
            for _ in range(iters):
                diff = X[:, None] - X[None]
                D = np.linalg.norm(diff, axis=-1) + 1e-9
                g = np.zeros_like(X)
                u = diff[i, j] / D[i, j][:, None]
                rij = ((D[i, j] - E) * W)[:, None] * u
                np.add.at(g, i, rij)
                np.add.at(g, j, -rij)
                dc = D[adj, adj + 1]
                uc = diff[adj, adj + 1] / dc[:, None]
                rc = (k_chain * (dc - chain))[:, None] * uc
                np.add.at(g, adj, rc)
                np.add.at(g, adj + 1, -rc)
                X = X - lr * g
            D = np.linalg.norm(X[:, None] - X[None], axis=-1)
            s = float((W * (D[i, j] - E) ** 2).sum()
                      + k_chain * ((D[adj, adj + 1] - chain) ** 2).sum())
            if s < best_s:
                best_s, best = s, X.copy()
        return best

    def matrix(self) -> np.ndarray:
        """Predicted expected CA-CA distance as a full (n, n) matrix."""
        M = np.zeros((self.n, self.n))
        M[self.i, self.j] = self.expected
        M[self.j, self.i] = self.expected
        return M


class CombinedDistogram(Distogram):
    """Mean of several distance models' scores, on a shared probability average.

    The per-pair MLP and the pair-tensor network make *different* mistakes -- on 6H7I the
    network's in-band Spearman is +0.79 against the MLP's +0.25, and on other targets the
    order reverses -- so averaging them beats either alone (in-band +0.425 / +0.412 -> 
    +0.446, top-20 2.75 / 2.76 -> 2.62 A). The averaged score is used for ranking; the
    averaged *distribution* is what `expected`, `realize` and the seed read, so the
    distance-geometry seed also benefits.
    """

    #: Weight on the members' disagreement, and a cautionary result. Penalising the spread
    #: between the two models -- so a structure scores well only if BOTH like it -- is a
    #: real improvement when the members are combined by RANK: in-band Spearman +0.462 ->
    #: +0.505 over 35 targets and top-20 2.58 -> 2.54 A (`work/disagree.py`). The same
    #: penalty applied to standardised SCORES, which is what the search loop needs (ranks
    #: depend on the batch and would make the objective non-stationary), does the opposite:
    #: +0.454 -> +0.422 at lambda 0.5 and +0.375 at 1.0. Score differences are unbounded and
    #: dominated by outlying pairs; rank differences are not. Default off.
    DISAGREE = 0.0

    def __init__(self, members: Sequence[Distogram],
                 ref_ca: Optional[np.ndarray] = None,
                 disagree: Optional[float] = None):
        m0 = members[0]
        prob = np.mean([m.prob for m in members], axis=0)
        super().__init__(m0.seq, prob, m0.i, m0.j)
        self.members = list(members)
        self.disagree = self.DISAGREE if disagree is None else float(disagree)
        # Members' scores live on different scales, so a plain average is dominated by
        # whichever model happens to have the larger spread. Standardising each member
        # against a fixed reference sample of structures makes the average an equal-weight
        # combination and keeps the objective stationary across VQE iterations (a
        # batch-relative normalisation would not).
        self.mu = np.zeros(len(self.members))
        self.sigma = np.ones(len(self.members))
        if ref_ca is not None and len(ref_ca) > 8:
            for k, m in enumerate(self.members):
                v = np.asarray(m.score(np.asarray(ref_ca, float)), float)
                self.mu[k] = v.mean()
                self.sigma[k] = max(v.std(), 1e-6)

    def score(self, ca: np.ndarray):
        arr = np.asarray(ca, float)
        single = arr.ndim == 2
        vals = np.stack([(np.asarray(m.score(arr), float) - self.mu[k]) / self.sigma[k]
                         for k, m in enumerate(self.members)])
        out = vals.mean(0)
        if self.disagree and len(vals) > 1:
            out = out + self.disagree * (vals.max(0) - vals.min(0))
        return float(out) if single else out
