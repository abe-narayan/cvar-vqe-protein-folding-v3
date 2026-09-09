"""Does the INPUT REPRESENTATION matter? Re-judged on SELECTION instead of distance MAE.

The error this corrects
`s5/dist_ablate.py` compared four inputs to the distance predictor and reported

    one-hot 2.216 | PCA-32 2.133 | PCA-128 2.212 | raw 1280-d + learned projection 2.109

all p >= 0.59. That table is the entire evidence for the stored conclusions
`esm-adds-nothing-for-short-peptides` and `sequence-signal-is-the-ceiling`.

**Every number in it is a distance MAE**, and this sprint has established that MAE does
not measure what the pipeline does:

  * `s7.audit` stage `noise`: error SHAPES injected into the oracle matrix at MATCHED
    MAE ~2.2 A select 2.09 A (iid), 2.68 A (coordinate-space) or 3.83 A (shrinkage).
    The same MAE buys a 1.7 A spread in the answer.
  * `s7.debias`: over 28 post-hoc scoring arms the best-MAE arm ranked 16th of 28 on
    selected RMSD, and the best in-band Spearman (pure NLL) had a WORSE selected RMSD.
  * `s7.bigdata`: 3x-10x more training data lowered training loss and RAISED dev MAE.
  * `s7.debias` stage `calib`: calibration slope +0.376 at correlation +0.283 -- the
    predicted deviations have ~75% of the right amplitude pointing the wrong way.

So the four-way ablation was judged by an instrument that does not measure the objective,
exactly the mistake already recorded as `better-matrix-worse-ranking`. Two representations
can tie on MAE and differ on ranking; that is now measured, not hypothesised. This module
re-runs the comparison with **selected CA-RMSD as the primary metric**.

ARMS -- they differ ONLY in the feature block
Every arm shares the architecture, the training corpus, the optimiser, the seed, the
number of epochs, the loss and the scoring code. Only the per-residue (and, for `escon`,
the per-pair) feature block changes.

    onehot   20-d one-hot. The control, and the arm `s5/dist_ablate` called `onehot`.
    phys     one-hot + the 5 Chou-Fasman/Kyte-Doolittle/charge/volume numbers already in
             `priors._PROPS`, z-scored. A cheap hand-crafted chemistry block: if ESM helps
             only by re-deriving hydrophobicity, this arm should capture most of it.
    pca32    one-hot + the SHIPPED 32-d ESM-2 PCA projection (`esm_features.embed`).
    pca128   one-hot + a 128-d ESM-2 PCA **refit per fold on training sequences only**,
             so no held-out sequence influences the projection.
    raw      one-hot + all 1280 ESM-2 dimensions, projected by the encoder's own first
             layer. This is the arm that can find informative low-variance directions a
             PCA discards.
    escon    pca32 plus PAIR features from ESM-2's contact head (the contact logit and
             its four diagonal neighbours, plus each residue's contact degree). This is
             the shipped distogram's input. It is reported SEPARATELY and flagged,
             because the contact head is the one part of ESM-2 that was supervised on
             structure (see LEAKAGE below); the one-hot vs ESM-embedding contrast is the
             clean one.

METRICS -- and which one decides
    selected CA-RMSD    PRIMARY. argmin of the score over the 500-candidate pool.
    in-band Spearman    candidates within 1.5 A of pool best. Necessary, NOT sufficient:
                        in `s7.debias` the best-in-band arm of 28 selected worse.
    calibration slope   regress (true - pool-mean) on (predicted - pool-mean) over pairs;
    and correlation     slope = corr x amplitude-ratio. The sharpest diagnostic of whether
                        a representation carries more real signal.
    MAE                 reported for completeness and explicitly NOT used to choose.

PROTOCOL
Tuning instrument: `s7.debias`'s 126 cluster-representative peptides, length 9-16,
identity-cluster-disjoint from BOTH the 24-target dev set and the 60-target benchmark,
pooled by `s7.audit.build_target`. Its base-arm SE is 0.147 A against dev's 0.354 A.
All arms, all seeds, all comparisons happen there. The SINGLE best arm is then run on the
24 dev targets ONCE, at the end. Seeds 0,1,2 per arm give a seed-variance estimate; with
6 arms and SE 0.147 A no single-seed difference is claimable without it.

DEPLOYABLE vs DIAGNOSTIC
DEPLOYABLE: `Encoder.residues`, `Encoder.pair_extra`, `ReprNet.forward`, `predict` and
`score_pool` -- their inputs are a SEQUENCE and the pool's own candidate distance
matrices. No native coordinate, distance, torsion or RMSD of the scored target reaches
them; `test_repr_select.py` asserts it by NaN-poisoning `Dnat` and `rr` in a real cached
record and requiring bit-identical scores.

DIAGNOSTIC: `target_stats` and `calibration` read the target's native `Dnat`/`rr`, but
only AFTER the score vector exists, and only to report. Nothing fitted on a native
crosses back into any arm -- there is no fitted scalar in this study at all.

LEAKAGE
* Training corpus for fold f = peptides with `peptide_db.folds(5)[seq] != f` plus
  `distogram._fold_fragments(f, 5)`, i.e. the production discipline: nothing above 0.6
  identity to a held-out peptide. `stage_leak` re-checks this by brute force.
* The pca128 projection is refit per fold on TRAINING sequences only.
* ESM-2 per-residue representations are a deterministic function of the amino-acid
  sequence: `esm_features.compute` feeds the model tokens and nothing else, and the cache
  is keyed by sequence. `stage_verify` checks (a) the cached arrays are (n, 1280) and
  (n, n) with no coordinate-shaped array anywhere, (b) distinct PDB entries that share a
  sequence but differ in conformation get IDENTICAL embeddings -- so the embedding cannot
  encode an individual structure -- and (c) optionally recomputes a few sequences with
  fair-esm and compares to the cache bit-for-bit.
* CAVEAT, stated because it is real: ESM-2's *contact head* is a logistic regression
  trained on PDB structures. Its weights are fixed and were not fitted on this project's
  held-out natives, but it is not a pure language-model channel. It appears only in the
  `escon` arm, which is reported apart from the primary contrast.

    python -m s7.repr_select verify        # ESM cache provenance checks (cheap)
    python -m s7.repr_select esmcon        # compact contact cache, throwaway subprocess
    python -m s7.repr_select train         # ARMS/SEEDS env-controlled; resumable
    python -m s7.repr_select tune          # score the 126-target instrument
    python -m s7.repr_select dev           # the ONE winning arm on 24 dev targets
    python -m s7.repr_select leak          # brute-force leakage audit
    python -m s7.repr_select report
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.stats import spearmanr

os.environ.setdefault("NT", "2")
import torch                                                              # noqa: E402
import torch.nn as nn                                                     # noqa: E402

torch.set_num_threads(int(os.environ.get("NT", "2")))

import distogram as dgm                                                   # noqa: E402
import peptide_db as db                                                   # noqa: E402
import priors                                                             # noqa: E402
from s7 import audit, debias                                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "repr_models")
RCACHE = os.path.join(HERE, "repr_cache")
ESMCON = os.path.join(RCACHE, "esmcon.npz")

K = debias.K                      #: 500, the pool every number in this sprint reports on
AA = "ARNDCQEGHILKMFPSTWYV"
AA_IDX = {a: i for i, a in enumerate(AA)}

#: Pair separation one-hot width. Peptides here are 9-16 residues, so |i-j| <= 15.
MAXSEP = 24

ARMS = ("onehot", "phys", "pca32", "pca128", "raw", "escon")
PRIMARY_ARMS = ("onehot", "phys", "pca32", "pca128", "raw")   #: no structure-supervised channel
SEEDS = (0, 1, 2)

#: Training hyperparameters. Fixed across every arm -- that is the whole point -- and
#: chosen to match the shipped predictor rather than to favour any arm: 40 epochs and
#: batch 4096 are `distogram.MLP.fit`'s defaults, the soft-binned cross-entropy with
#: smooth 0.6 is its loss, and the 17 bins are `priors.BIN_EDGES`.
EPOCHS = int(os.environ.get("EPOCHS", "40"))
BATCH = 4096
LR = 2e-3
SMOOTH = 0.6
D_RES = 64
D_PAIR = 128


# --------------------------------------------------------------------- residue blocks
def onehot(seq):
    h = np.zeros((len(seq), 20), np.float32)
    for k, c in enumerate(seq):
        t = AA_IDX.get(c)
        if t is not None:
            h[k, t] = 1.0
    return h


def _phys_table():
    """(20, 5) z-scored physicochemical properties, in `AA` order."""
    P = np.array([priors._PROPS[a] for a in AA], np.float32)
    return (P - P.mean(0)) / (P.std(0) + 1e-6)


_PHYS = _phys_table()


def phys(seq):
    idx = np.array([AA_IDX.get(c, 0) for c in seq], int)
    return _PHYS[idx]


# --------------------------------------------------------------------- ESM caches
def _esm32():
    from s5 import esm32
    return esm32.load()


def _esmraw():
    from s5 import esmraw
    return esmraw.load()


_CON = [None]


def esmcon():
    """dict sequence -> (n, n) float32 ESM-2 contact logits. Built by `stage_esmcon`."""
    if _CON[0] is None:
        z = np.load(ESMCON, allow_pickle=True)
        k, f, o = z["keys"], z["flat"], z["offs"]
        _CON[0] = {str(s): f[o[i]:o[i + 1]].reshape(len(str(s)), len(str(s)))
                   for i, s in enumerate(k)}
    return _CON[0]


def stage_esmcon():
    """Extract contact maps for every training sequence into a small cache, then exit.

    Same pattern as `s5.esm32`/`s5.esmraw`: the 1.5 GB `esm_cache.npz` is materialised in
    a process that does nothing else and dies immediately, so no long-lived job ever holds
    it. Contact maps for ~6,800 peptides of median length 13 are a few MB.
    """
    os.makedirs(RCACHE, exist_ok=True)
    if os.path.exists(ESMCON):
        print(f"{ESMCON} exists")
        return
    import esm_features as ef
    import fragment_db as fdb
    seqs = sorted({p.seq for p in db.load()} | {f.seq for f in fdb.load()})
    print(f"contact maps for {len(seqs)} sequences", flush=True)
    keys, vals, offs = [], [], [0]
    for k, s in enumerate(seqs):
        c = np.asarray(ef.contacts(s), np.float32)
        assert c.shape == (len(s), len(s)), (s, c.shape)
        keys.append(s)
        vals.append(c.reshape(-1))               # (n, n) maps vary in n: store flat
        offs.append(offs[-1] + c.size)
        if (k + 1) % 1000 == 0:
            print(f"  {k+1}/{len(seqs)}", flush=True)
    flat = np.concatenate(vals, 0)
    tmp = ESMCON + ".tmp.npz"
    np.savez_compressed(tmp, keys=np.array(keys, dtype=object), flat=flat,
                        offs=np.array(offs, np.int64))
    os.replace(tmp, ESMCON)
    print(f"wrote {ESMCON}: {flat.shape}, {flat.nbytes/1e6:.1f} MB", flush=True)


# --------------------------------------------------------------------- pca128
def pca128_path(fold):
    return os.path.join(RCACHE, f"pca128_fold{fold}.npz")


def fit_pca128(fold, n_folds=5, n_comp=128, verbose=True):
    """Refit the ESM projection on fold `fold`'s TRAINING sequences only.

    Using the shipped `esm_pca.npz` would let every held-out sequence influence the
    projection -- a weak leak but a real and trivially avoidable one. Randomised range
    finder rather than a full SVD of a (92k, 1280) matrix.
    """
    os.makedirs(RCACHE, exist_ok=True)
    if os.path.exists(pca128_path(fold)):
        return
    t0 = time.time()
    train = sorted({e.seq for e in train_entries(fold, n_folds)})
    raw = _esmraw()
    X = np.concatenate([raw[s] for s in train if s in raw], 0).astype(np.float32)
    mu = X.mean(0)
    g = np.random.default_rng(0).normal(size=(X.shape[1], n_comp + 32)).astype(np.float32)
    Y = (X - mu) @ g
    Q, _ = np.linalg.qr(Y)
    B = Q.T @ (X - mu)
    _, S, Vt = np.linalg.svd(B, full_matrices=False)
    np.savez(pca128_path(fold), mu=mu, W=Vt[:n_comp].T.astype(np.float32),
             scale=(S[:n_comp] / math.sqrt(len(X))).astype(np.float32))
    del X, Y, Q, B
    if verbose:
        print(f"  pca128 fold {fold}: {len(train)} training sequences "
              f"({time.time()-t0:.0f}s)", flush=True)


# --------------------------------------------------------------------- encoders
class Encoder:
    """DEPLOYABLE. Sequence -> per-residue features (and, for `escon`, pair features).

    Its only input is a string. It holds a fold index solely to pick the fold's own
    training-set PCA basis; it never sees a coordinate.
    """

    def __init__(self, arm, fold, n_folds=5):
        assert arm in ARMS, arm
        self.arm, self.fold = arm, fold
        self._pca = None
        if arm == "pca128":
            fit_pca128(fold, n_folds, verbose=False)
            z = np.load(pca128_path(fold))
            self._pca = (z["mu"], z["W"], z["scale"])

    @property
    def d_in(self):
        return {"onehot": 20, "phys": 25, "pca32": 52, "pca128": 148,
                "raw": 1300, "escon": 52}[self.arm]

    @property
    def d_extra(self):
        return 7 if self.arm == "escon" else 0

    def residues(self, seq):
        h = onehot(seq)
        a = self.arm
        if a == "onehot":
            return h
        if a == "phys":
            return np.hstack([h, phys(seq)])
        if a in ("pca32", "escon"):
            e = _esm32().get(seq)
        elif a == "raw":
            e = _esmraw().get(seq)
        else:
            mu, W, sc = self._pca
            r = _esmraw().get(seq)
            e = None if r is None else ((r - mu) @ W) / np.maximum(sc, 1e-6)
        if e is None:
            raise KeyError(f"{self.arm}: no cached ESM vector for {seq!r} -- rebuild "
                           f"s5/esm32.npz or s5/esmraw.npz in a throwaway subprocess")
        e = np.asarray(e, np.float32)[:len(seq)]
        assert len(e) == len(seq), (seq, e.shape)
        return np.hstack([h, e]).astype(np.float32)

    def pair_extra(self, seq, i, j):
        """(npairs, d_extra). ESM-2 contact-head features -- `escon` only."""
        if self.arm != "escon":
            return np.zeros((len(i), 0), np.float32)
        con = esmcon()[seq]
        n = len(seq)
        nb = np.stack([con[np.clip(i + a, 0, n - 1), np.clip(j + b, 0, n - 1)]
                       for a, b in ((-1, -1), (1, 1), (-1, 1), (1, -1))], 1)
        deg = con.sum(1) / max(1, n - 1)
        return np.column_stack([con[i, j], nb, deg[i], deg[j]]).astype(np.float32)


def sep_features(sep):
    """(npairs, MAXSEP + 2). Identical for every arm."""
    s1 = np.zeros((len(sep), MAXSEP), np.float32)
    s1[np.arange(len(sep)), np.clip(sep - 2, 0, MAXSEP - 1)] = 1.0
    sf = np.stack([sep / 20.0, np.log1p(sep) / 3.0], 1).astype(np.float32)
    return np.hstack([s1, sf])


# --------------------------------------------------------------------- model
class ReprNet(nn.Module):
    """Per-residue encoder, then a symmetric pair MLP over the 17 distance bins.

    Binned rather than the direct distance regression `s5/dist_ablate` used, because the
    quantity under test is now SELECTION and selection runs through the shipped
    Bayes-risk score, which needs a distribution and its spread-derived pair weights. The
    bins are `priors.BIN_EDGES`, the loss the same soft-binned cross-entropy as
    `distogram.MLP`.

    `hi + hj` and `|hi - hj|` because a distance is symmetric in (i, j) and the features
    must be too. The first layer IS the learned projection the `raw` arm needs.
    """

    def __init__(self, d_in, d_extra=0, d_res=D_RES, d_pair=D_PAIR):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(d_in, d_res), nn.GELU(),
                                 nn.Linear(d_res, d_res))
        self.pair = nn.Sequential(
            nn.Linear(2 * d_res + MAXSEP + 2 + d_extra, d_pair), nn.GELU(),
            nn.Linear(d_pair, d_pair), nn.GELU(),
            nn.Linear(d_pair, dgm.NBINS))

    def forward(self, Ri, Rj, extra):
        hi, hj = self.enc(Ri), self.enc(Rj)
        return self.pair(torch.cat([hi + hj, (hi - hj).abs(), extra], -1))


# --------------------------------------------------------------------- training data
def train_entries(fold, n_folds=5):
    """Out-of-fold peptides plus the fold's identity-filtered fragments. Production rule."""
    folds = db.folds(n_folds)
    return ([p for p in db.load() if folds[p.seq] != fold]
            + list(dgm._fold_fragments(fold, n_folds)))


def build_training_arrays(entries, enc, verbose=False):
    """(R float16 residue table, ii, jj, extra, ybin). Memory-lean by construction.

    The per-residue table is stored ONCE in float16 and gathered per minibatch; the
    alternative -- materialising (npairs, 2*d_in) features -- is 6.5 GB for the `raw` arm
    on a box that has repeatedly breached 90% RAM.
    """
    Rs, II, JJ, EX, YB, off = [], [], [], [], [], 0
    for e in entries:
        n = len(e.seq)
        if n < 4:
            continue
        i, j = np.triu_indices(n, k=2)
        d = np.linalg.norm(e.ca[i] - e.ca[j], axis=-1)
        Rs.append(enc.residues(e.seq).astype(np.float16))
        II.append((i + off).astype(np.int32))
        JJ.append((j + off).astype(np.int32))
        EX.append(np.hstack([sep_features(j - i), enc.pair_extra(e.seq, i, j)]))
        YB.append(np.digitize(d, dgm.BIN_EDGES).astype(np.int8))
        off += n
    R = np.concatenate(Rs, 0)
    out = (R, np.concatenate(II), np.concatenate(JJ),
           np.concatenate(EX, 0).astype(np.float32), np.concatenate(YB))
    if verbose:
        print(f"    {len(entries)} chains, {R.shape[0]} residues, {len(out[1])} pairs, "
              f"d_in {R.shape[1]}, table {R.nbytes/1e6:.0f} MB", flush=True)
    return out


def _chunk_moments(R, chunk=8192):
    """Column mean/sd of a float16 table without materialising a float32 copy of it.

    `R` is 471 MB in float32 for the `raw` arm; on a box that has breached 90% RAM three
    times, a whole second copy just to take a mean is not affordable.
    """
    n, d = R.shape
    s1 = np.zeros(d, np.float64)
    s2 = np.zeros(d, np.float64)
    for a in range(0, n, chunk):
        b = R[a:a + chunk].astype(np.float64)
        s1 += b.sum(0)
        s2 += (b * b).sum(0)
    mu = s1 / n
    var = np.maximum(s2 / n - mu * mu, 0.0)
    return mu.astype(np.float32), (np.sqrt(var) + 1e-3).astype(np.float32)


def soft_targets(y, smooth=SMOOTH):
    idx = np.arange(dgm.NBINS)[None, :]
    T = np.exp(-((idx - y[:, None].astype(np.float32)) ** 2) / (2 * smooth ** 2))
    return (T / T.sum(1, keepdims=True)).astype(np.float32)


def model_path(arm, seed, fold):
    return os.path.join(MODELS, f"{arm}_s{seed}_f{fold}.pt")


def train_one(arm, seed, fold, epochs=EPOCHS, verbose=True):
    """Train (or load) one (arm, seed, fold) model. Only `arm` changes the features."""
    os.makedirs(MODELS, exist_ok=True)
    path = model_path(arm, seed, fold)
    enc = Encoder(arm, fold)
    if os.path.exists(path):
        z = torch.load(path, weights_only=False)
        net = ReprNet(enc.d_in, enc.d_extra)
        net.load_state_dict(z["sd"])
        net.eval()
        return net, enc, z["mu"], z["sd_"]
    t0 = time.time()
    R, ii, jj, EX, yb = build_training_arrays(train_entries(fold), enc, verbose=verbose)
    mu, sd = _chunk_moments(R)
    T = soft_targets(yb)
    torch.manual_seed(seed)
    net = ReprNet(enc.d_in, enc.d_extra)
    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=1e-4)
    npair = len(ii)
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=LR, total_steps=epochs * max(1, npair // BATCH + 1))
    Rt = torch.from_numpy(R)
    mut, sdt = torch.from_numpy(mu), torch.from_numpy(sd)
    iit, jjt = torch.from_numpy(ii.astype(np.int64)), torch.from_numpy(jj.astype(np.int64))
    EXt, Tt = torch.from_numpy(EX), torch.from_numpy(T)
    g = torch.Generator().manual_seed(seed)
    net.train()
    for ep in range(epochs):
        perm = torch.randperm(npair, generator=g)
        tot = 0.0
        for s in range(0, npair, BATCH):
            b = perm[s:s + BATCH]
            Ri = (Rt[iit[b]].float() - mut) / sdt
            Rj = (Rt[jjt[b]].float() - mut) / sdt
            logit = net(Ri, Rj, EXt[b])
            loss = -(Tt[b] * torch.log_softmax(logit, 1)).sum(1).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
            sched.step()
            tot += float(loss.detach()) * len(b)
        if verbose and (ep % 10 == 0 or ep == epochs - 1):
            print(f"    epoch {ep:3d} loss {tot/npair:.4f}", flush=True)
    net.eval()
    tmp = path + ".tmp"
    torch.save({"sd": net.state_dict(), "mu": mu, "sd_": sd,
                "arm": arm, "seed": seed, "fold": fold, "epochs": epochs}, tmp)
    os.replace(tmp, path)
    if verbose:
        print(f"  {arm:7} seed {seed} fold {fold}: {time.time()-t0:.0f}s", flush=True)
    del R, Rt, T, Tt, EX, EXt
    return net, enc, mu, sd


# --------------------------------------------------------------------- inference
def predict(seq, net, enc, mu, sd):
    """DEPLOYABLE. Sequence -> `distogram.Distogram`. No structure of any kind enters."""
    n = len(seq)
    i, j = np.triu_indices(n, k=2)
    R = enc.residues(seq).astype(np.float32)
    Rn = torch.from_numpy((R - mu) / sd)
    EX = torch.from_numpy(np.hstack([sep_features(j - i), enc.pair_extra(seq, i, j)]))
    with torch.no_grad():
        logit = net(Rn[torch.from_numpy(i.astype(np.int64))],
                    Rn[torch.from_numpy(j.astype(np.int64))], EX)
        prob = torch.softmax(logit, 1).numpy().astype(np.float64)
    return dgm.Distogram(seq, prob, i, j)


def score_pool(dg, D):
    """DEPLOYABLE. The SHIPPED Bayes-risk score of a pool's candidate distance matrices."""
    return debias.score_risk(dg._risk, dg.grid, np.asarray(D, float))


# --------------------------------------------------------------------- statistics
def target_stats(dg, rec):
    """DIAGNOSTIC. Everything reported for one target. Natives enter only here."""
    D = rec["D"][:K].astype(float)
    rr = rec["rr"][:K].astype(float)
    Dn = rec["Dnat"].astype(float)
    st = debias.rank_stats(score_pool(dg, D), rr)
    Dbar = D.mean(0)
    x, y = dg.expected - Dbar, Dn - Dbar
    st.update({
        "pdb": rec["pdb"], "n": rec["n"], "fold": rec["fold"],
        "mae": float(np.abs(dg.expected - Dn).mean()),
        "slope": float(np.polyfit(x, y, 1)[0]) if len(x) > 2 else float("nan"),
        "corr": float(np.corrcoef(x, y)[0, 1]),
        "sd_pred_dev": float(x.std()), "sd_true_dev": float(y.std()),
        "npairs": int(len(x)),
    })
    return st, x, y


def paired(a, b):
    return debias.paired(a, b)


def aggregate(rows, arm, seed, X=None, Y=None):
    sel = np.array([r["sel"] for r in rows])
    rb = np.array([r["rho_band"] if r["rho_band"] is not None else np.nan for r in rows])
    nb = np.array([r["n_band"] for r in rows], float)
    ok = ~np.isnan(rb)
    out = {"arm": arm, "seed": seed, "n_targets": len(rows),
           "sel": float(sel.mean()),
           "sel_se": float(sel.std(ddof=1) / np.sqrt(len(sel))),
           "pool": float(np.mean([r["pool"] for r in rows])),
           "mae": float(np.mean([r["mae"] for r in rows])),
           "rho_global": float(np.mean([r["rho_global"] for r in rows])),
           "rho_band": float(np.nanmean(rb)) if ok.any() else None,
           "rho_band_sizeweighted": (float((rb[ok] * nb[ok]).sum() / nb[ok].sum())
                                     if ok.any() else None),
           "n_band_scored": int(ok.sum()),
           "slope_per_target_median": float(np.median([r["slope"] for r in rows])),
           "corr_per_target_median": float(np.median([r["corr"] for r in rows]))}
    if X is not None:
        Xc, Yc = np.concatenate(X), np.concatenate(Y)
        s, ic = np.polyfit(Xc, Yc, 1)
        out.update({"slope_pooled": float(s), "intercept_pooled": float(ic),
                    "corr_pooled": float(np.corrcoef(Xc, Yc)[0, 1]),
                    "sd_pred_dev": float(Xc.std()), "sd_true_dev": float(Yc.std()),
                    "n_pairs": int(len(Xc))})
    return out


# --------------------------------------------------------------------- stages
def _arms():
    return tuple(os.environ.get("ARMS", ",".join(ARMS)).split(","))


def _seeds():
    return tuple(int(s) for s in os.environ.get("SEEDS", "0,1,2").split(","))


def _mem_pct():
    """Percent of physical RAM in use, or None if it cannot be read."""
    try:
        import ctypes

        class MS(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        m = MS()
        m.dwLength = ctypes.sizeof(MS)
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return int(m.dwMemoryLoad)
    except Exception:                                                    # noqa: BLE001
        return None


def _mem_guard(limit=90):
    p = _mem_pct()
    if p is not None and p >= limit:
        print(f"RAM at {p}% (limit {limit}%) -- stopping rather than pushing through.",
              flush=True)
        return False
    return True


def stage_train():
    """Train every (arm, seed, fold) not already on disk. Resumable, one model at a time."""
    for seed in _seeds():
        for arm in _arms():
            for fold in range(5):
                if os.path.exists(model_path(arm, seed, fold)):
                    continue
                if not _mem_guard():
                    return
                train_one(arm, seed, fold)


def _run(recs, arm, seed, tag=""):
    rows, X, Y = [], [], []
    by_fold = {}
    for k, rec in enumerate(recs):
        f = rec["fold"]
        if f not in by_fold:
            by_fold = {f: train_one(arm, seed, f, verbose=False)}      # one model resident
        net, enc, mu, sd = by_fold[f]
        dg = predict(rec["seq"], net, enc, mu, sd)
        st, x, y = target_stats(dg, rec)
        rows.append(st)
        X.append(x)
        Y.append(y)
        if (k + 1) % 25 == 0 or k == len(recs) - 1:
            print(f"  [{k+1:3d}/{len(recs)}] {arm} s{seed}{tag} "
                  f"running sel {np.mean([r['sel'] for r in rows]):.3f}", flush=True)
    return rows, X, Y


def stage_tune():
    """Every arm x seed on the 126-target instrument. Incremental: appends per cell."""
    recs = [debias._load(p.pdb) for p in debias.tuning_targets() if debias._has(p.pdb)]
    recs.sort(key=lambda r: r["fold"])
    print(f"{len(recs)} tuning targets", flush=True)
    out = (audit.read_json("repr_tune.json")
           if os.path.exists(os.path.join(HERE, "repr_tune.json"))
           else {"n_targets": len(recs), "K": K, "epochs": EPOCHS, "cells": []})
    have = {(c["arm"], c["seed"]) for c in out["cells"]}
    for seed in _seeds():
        for arm in _arms():
            if (arm, seed) in have:
                continue
            if not _mem_guard():
                return out
            t0 = time.time()
            rows, X, Y = _run(recs, arm, seed)
            cell = aggregate(rows, arm, seed, X, Y)
            cell["wall"] = round(time.time() - t0, 1)
            cell["per_target"] = rows
            out["cells"].append(cell)
            audit.write_json("repr_tune.json", out)
            print(f"{arm:8} seed {seed}: sel {cell['sel']:.3f} +- {cell['sel_se']:.3f}  "
                  f"MAE {cell['mae']:.3f}  "
                  f"in-band {cell['rho_band'] if cell['rho_band'] is None else round(cell['rho_band'], 3)}  "
                  f"corr {cell['corr_pooled']:+.3f}  slope {cell['slope_pooled']:+.3f}",
                  flush=True)
    return out


def _cells(out, arm):
    return [c for c in out["cells"] if c["arm"] == arm]


def _sel_matrix(out, arm, order=None, key="sel"):
    """(n_seeds, n_targets) per-target values, target order pinned by pdb id.

    `order` is passed in for every arm but the reference one, so a paired difference can
    never silently compare target k of one arm against target k of another.
    """
    cs = sorted(_cells(out, arm), key=lambda c: c["seed"])
    if not cs:
        return None, []
    if order is None:
        order = [r["pdb"] for r in cs[0]["per_target"]]
    M = []
    for c in cs:
        d = {r["pdb"]: r[key] for r in c["per_target"]}
        assert all(p in d for p in order), (arm, c["seed"])
        M.append([d[p] for p in order])
    return np.array(M), order


def stage_dev():
    """ONE pass on the 24 dev targets with the arm that won the tuning instrument."""
    out = audit.read_json("repr_tune.json")
    ranked = sorted(
        [(float(np.mean([c["sel"] for c in _cells(out, a)])), a)
         for a in PRIMARY_ARMS if _cells(out, a)])
    best = ranked[0][1]
    base = "onehot"
    arms = [base] if best == base else [base, best]
    print(f"dev run, arms {arms} (winner chosen on the 126-target instrument)", flush=True)
    recs = [debias._load(p.pdb, tune=False) for p in debias.dev_targets()]
    recs.sort(key=lambda r: r["fold"])
    dev = {"n_targets": len(recs), "K": K, "winning_arm": best,
           "tuning_ranking": ranked, "cells": []}
    for arm in arms:
        per_seed = []
        for seed in _seeds():
            rows, X, Y = _run(recs, arm, seed, tag=" [dev]")
            c = aggregate(rows, arm, seed, X, Y)
            c["per_target"] = rows
            per_seed.append(c)
        dev["cells"].extend(per_seed)
    audit.write_json("repr_dev.json", dev)
    _print_dev(dev)
    return dev


def stage_oracle():
    """DIAGNOSTIC reference rows for the same 126 pools: the ceiling and the floor.

    `native` ranks by weighted L1 against the target's OWN distance matrix. It is not an
    arm and could never be deployed -- it is the best any distance-channel ranker could
    do, and the number every representation is being asked to approach. `pool` is the
    best candidate present. Neither influences any arm's score.
    """
    recs = [debias._load(p.pdb) for p in debias.tuning_targets() if debias._has(p.pdb)]
    rows = []
    for rec in recs:
        D = rec["D"][:K].astype(float)
        rr = rec["rr"][:K].astype(float)
        st = debias.rank_stats(audit.score_l1(D, rec["Dnat"].astype(float)), rr)
        st["pdb"] = rec["pdb"]
        rows.append(st)
    out = {"n_targets": len(rows),
           "native_sel": float(np.mean([r["sel"] for r in rows])),
           "native_sel_se": float(np.std([r["sel"] for r in rows], ddof=1)
                                  / np.sqrt(len(rows))),
           "native_rho_band": float(np.nanmean(
               [r["rho_band"] if r["rho_band"] is not None else np.nan for r in rows])),
           "native_rho_global": float(np.mean([r["rho_global"] for r in rows])),
           "pool": float(np.mean([r["pool"] for r in rows])),
           "per_target": rows}
    audit.write_json("repr_oracle.json", out)
    print(f"native-matrix ranking selects {out['native_sel']:.3f} +- "
          f"{out['native_sel_se']:.3f} A, in-band rho {out['native_rho_band']:+.3f}; "
          f"pool best {out['pool']:.3f} A")
    return out


def stage_leak():
    """Brute-force audit: no training chain is above the identity threshold to a target."""
    tg = debias.tuning_targets() + list(debias.dev_targets())
    folds = db.folds(5)
    bad = []
    for fold in range(5):
        ents = train_entries(fold)
        seqs = [e.seq for e in ents]
        ks = [db._kmers(s) for s in seqs]
        for t in tg:
            if folds[t.seq] != fold:
                continue
            tk = db._kmers(t.seq)
            for s, k in zip(seqs, ks):
                if tk and k and not (tk & k):
                    continue
                if db.identity(t.seq, s) >= db.IDENTITY_THRESHOLD:
                    bad.append((t.pdb, fold, s))
                    break
        print(f"  fold {fold}: {len(ents)} training chains checked", flush=True)
    res = {"n_targets": len(tg), "n_violations": len(bad), "violations": bad[:20]}
    audit.write_json("repr_leak.json", res)
    print(f"{len(tg)} targets, {len(bad)} identity violations")
    return res


def stage_verify():
    """The ESM cache carries no structure: shapes, and same sequence -> same embedding."""
    raw = _esmraw()
    e32 = _esm32()
    res = {"n_raw": len(raw), "n_pca32": len(e32)}
    dims = {v.shape[1] for v in raw.values()}
    res["raw_dims"] = sorted(dims)
    assert dims == {1280}, dims
    # (a) every cached array is (n, 1280) with n = len(sequence): nothing coordinate-shaped
    bad = [s for s, v in list(raw.items())[:2000] if v.shape != (len(s), 1280)]
    res["shape_violations"] = len(bad)
    assert not bad, bad[:3]
    # (b) the featuriser is fed tokens and nothing else. Source-level, because it is the
    #     only place a coordinate could possibly enter, and it is eight lines long.
    import inspect

    import esm_features as ef
    src = inspect.getsource(ef.compute)
    res["compute_mentions_coordinates"] = any(
        t in src for t in (".ca", "coords", "native", "pdb.", "distance"))
    assert not res["compute_mentions_coordinates"], src
    assert "bc([" in src and "toks" in src, "expected a token-only batch converter call"
    # (c) DECISIVE: recompute a few cached vectors from the SEQUENCE ALONE with fair-esm
    #     and compare. If the cache held anything structural this cannot reproduce it.
    #     Off by default -- the 650M model is ~2.6 GB and this box runs other jobs.
    res["recompute"] = None
    if os.environ.get("RECOMPUTE"):
        res["recompute"] = _recompute_check(raw)
    # (d) the projection is a pure function of the raw cache.
    z = np.load(ef._PCA_PATH)
    pmu, pW, psc = z["mu"], z["W"], z["scale"]
    ks = [s for s in list(raw)[:200] if s in e32][:20]
    res["n_pca32_checked"] = len(ks)
    res["pca32_matches_projection"] = bool(all(
        np.allclose(e32[s], ((raw[s] - pmu) @ pW) / np.maximum(psc, 1e-6), atol=2e-2)
        for s in ks))
    assert res["pca32_matches_projection"], "esm32 cache is not the PCA of the raw cache"
    audit.write_json("repr_verify.json", res)
    print(json.dumps(res, indent=1))
    return res


def _recompute_check(raw, n_seqs=3):
    """Run ESM-2 on the raw sequence strings and compare with the cache."""
    import esm
    import torch as t
    seqs = [s for s in list(raw)[:400] if len(s) <= 16][:n_seqs]
    model, alphabet = getattr(esm.pretrained, "esm2_t33_650M_UR50D")()
    model.eval()
    bc = alphabet.get_batch_converter()
    layer = model.num_layers
    out = []
    with t.no_grad():
        _, _, toks = bc([(f"p{k}", s) for k, s in enumerate(seqs)])
        r = model(toks, repr_layers=[layer])["representations"][layer].numpy()
    for k, s in enumerate(seqs):
        fresh = r[k, 1:len(s) + 1].astype(np.float16).astype(np.float32)
        out.append({"seq": s, "max_abs_diff": float(np.abs(fresh - raw[s]).max())})
    del model
    return out


# --------------------------------------------------------------------- report
def _print_arm_table(out, title, base="onehot"):
    print(f"\n=== {title} ===")
    print(f"{'arm':8} {'seeds':>5} {'sel':>7} {'+-se':>6} {'sd_seed':>8} "
          f"{'d(onehot)':>10} {'ci95':>16} {'W/L':>9} {'in-band':>8} "
          f"{'corr':>7} {'slope':>7} {'MAE':>6}")
    Mb, order = _sel_matrix(out, base)
    rows = []
    for arm in ARMS:
        cs = _cells(out, arm)
        if not cs:
            continue
        M, _ = _sel_matrix(out, arm, order if arm != base else None)
        sel = M.mean(0)                       # seed-averaged, per target
        d = paired(sel, Mb.mean(0)) if (Mb is not None and arm != base) else None
        rows.append((float(sel.mean()), arm, cs, M, d))
    rows.sort()
    for mean_sel, arm, cs, M, d in rows:
        se = float(M.mean(0).std(ddof=1) / np.sqrt(M.shape[1]))
        sdseed = float(M.mean(1).std(ddof=1)) if M.shape[0] > 1 else 0.0
        rb = float(np.mean([c["rho_band"] for c in cs]))
        co = float(np.mean([c["corr_pooled"] for c in cs]))
        sl = float(np.mean([c["slope_pooled"] for c in cs]))
        mae = float(np.mean([c["mae"] for c in cs]))
        star = "*" if arm not in PRIMARY_ARMS else " "
        line = (f"{arm:7}{star} {M.shape[0]:5d} {mean_sel:7.3f} {se:6.3f} {sdseed:8.3f} ")
        if d is not None:
            line += (f"{d['mean_diff']:+10.3f} [{d['ci95'][0]:+.3f},{d['ci95'][1]:+.3f}] "
                     f"{d['n_better']:3d}/{d['n_worse']:<5d} ")
        else:
            line += " " * 38
        line += f"{rb:+8.3f} {co:+7.3f} {sl:+7.3f} {mae:6.3f}"
        print(line)
    print("  * contains ESM-2's structure-supervised contact head; not a clean "
          "sequence-only arm")
    if rows:
        pool = float(np.mean([r["pool"] for r in rows[0][2][0]["per_target"]]))
        print(f"  pool best (a perfect ranker): {pool:.3f} A")
    op = os.path.join(HERE, "repr_oracle.json")
    if os.path.exists(op):
        o = audit.read_json("repr_oracle.json")
        print(f"  ranking by the target's OWN distance matrix (not deployable, the "
              f"ceiling of this channel): {o['native_sel']:.3f} A, "
              f"in-band {o['native_rho_band']:+.3f}")


def _print_dev(dev):
    print(f"\n=== DEV SET ({dev['n_targets']} targets, K={dev['K']}) -- one pass ===")
    print(f"  winner on the tuning instrument: {dev['winning_arm']}")
    print(f"{'arm':8} {'seeds':>5} {'sel':>7} {'+-se':>6} {'in-band':>8} "
          f"{'corr':>7} {'slope':>7} {'MAE':>6}")
    for arm in ARMS:
        cs = _cells(dev, arm)
        if not cs:
            continue
        M, _ = _sel_matrix(dev, arm)
        sel = M.mean(0)
        print(f"{arm:8} {M.shape[0]:5d} {sel.mean():7.3f} "
              f"{sel.std(ddof=1)/np.sqrt(len(sel)):6.3f} "
              f"{np.mean([c['rho_band'] for c in cs]):+8.3f} "
              f"{np.mean([c['corr_pooled'] for c in cs]):+7.3f} "
              f"{np.mean([c['slope_pooled'] for c in cs]):+7.3f} "
              f"{np.mean([c['mae'] for c in cs]):6.3f}")


def stage_report():
    p = os.path.join(HERE, "repr_tune.json")
    if os.path.exists(p):
        out = audit.read_json("repr_tune.json")
        _print_arm_table(out, f"TUNING INSTRUMENT ({out['n_targets']} targets, K={out['K']})")
        print("\n  seed-by-seed selected RMSD")
        for arm in ARMS:
            cs = sorted(_cells(out, arm), key=lambda c: c["seed"])
            if cs:
                print(f"  {arm:8} " + "  ".join(f"s{c['seed']}={c['sel']:.3f}" for c in cs))
        print("\n  the MAE column, for the record -- and it is NOT what chose anything")
        print("  s5/dist_ablate (dev MAE): onehot 2.216  pca32 2.133  pca128 2.212  "
              "raw 2.109")
    else:
        print("(tune not run)")
    p = os.path.join(HERE, "repr_dev.json")
    if os.path.exists(p):
        _print_dev(audit.read_json("repr_dev.json"))


STAGES = {"esmcon": stage_esmcon, "verify": stage_verify, "train": stage_train,
          "tune": stage_tune, "dev": stage_dev, "leak": stage_leak,
          "oracle": stage_oracle, "report": stage_report}


def main(argv):
    if len(argv) < 2 or argv[1] not in STAGES:
        print(__doc__)
        print("stages: " + ", ".join(STAGES))
        return 1
    STAGES[argv[1]]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
