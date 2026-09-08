"""s17/sel_lib.py -- shared machinery for the SELECT workstream.

THE ONE JOB OF THIS MODULE: make it impossible to compare two selectors on different
candidate sets, and impossible to read a tie as a signal.

Everything here is native-free except `rr`, which is carried explicitly under that name and
is used for EVALUATION ONLY.  Nothing in this module lets a caller put `rr` into a score.

WHAT A "PACK" IS.  One target, one candidate width K, everything a selector could want:

    idx        (K,)      indices into the universe, in BLOSUM-ranked order (the shipped order)
    rr         (K,)      ORACLE CA-RMSD of each candidate to the native.  EVALUATION ONLY.
    D          (K,P)     CA-CA distances of every scored pair, float32
    sep        (P,)      sequence separation j - i of each pair
    prob       (P,17)    the shipped leave-fold-out distogram's per-pair distribution
    expected   (P,)      its mean, sd (P,) its sd
    risk_raw   (P,G)     UNWEIGHTED L1 Bayes risk on the 2..40 A grid
    w_ship     (P,)      the shipped per-pair weight, 1/(sd+0.5) normalised to mean 1

`risk_raw` and `w_ship` are recomputed here from `prob` rather than read from the cached
weighted table, so a variant can change the weighting.  `shipped_score` is reproduced
bit-comparably (asserted in `selfcheck`) so no variant is being compared against a
slightly-different baseline.

THE TIE TRAP (S8, and it cost a whole table).  `np.argmin` returns the FIRST index of a tie,
and the candidate order here is BLOSUM rank -- informative.  A score that is constant, or
integral, or exactly zero on many candidates, will therefore appear to "select" whatever the
ranking put first.  `sel_of` averages the outcome over the WHOLE tied argmin set, which is
the expectation under random tie-breaking; a constant score then scores exactly the set mean,
which is what a zero-information signal must score.

THE MATCHED-RANDOM CONTROL.  Random selection of one candidate from the same set has
expectation exactly `rr.mean()`.  This module uses that EXACT expectation rather than a few
seeded draws, because a 3-draw control carries its own standard error of ~1 A per target and
that noise lands directly in the paired difference every selector is priced by.  A seeded
finite-draw version is provided as `rand_draws` for anyone who needs the sampling
distribution rather than its mean.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
CACHE = os.path.join(HERE, "cache")
RESULTS = os.path.join(HERE, "results")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402

_TG = None


def targets():
    global _TG
    if _TG is None:
        _TG = I.targets()
    return _TG


def folds():
    return np.array([int(t["fold"]) for t in targets()], int)


def names():
    return [t["pdb"] for t in targets()]


# --------------------------------------------------------------------------- packs
def _wship(sd):
    w = 1.0 / np.maximum(sd + 0.5, 1e-6)
    return w / max(w.mean(), 1e-12)


def pack(pdb, K=500, want=("D",)):
    """Candidate pack for one target at width K.  K = 0 means the whole universe."""
    t = {x["pdb"]: x for x in targets()}[pdb]
    n, fold, seq = int(t["n"]), int(t["fold"]), t["seq"]
    u = I.load_univ(pdb)
    order = np.asarray(u["order"], int)
    U = len(order)
    k = U if K in (0, None) else min(int(K), U)
    idx = order[:k]

    dg = I.distogram(pdb, seq, fold)
    i, j = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
    prob = np.asarray(dg["prob"], float)
    centres = np.asarray(dg["centres"], float)
    grid = np.asarray(dg["grid"], float)
    sd = np.asarray(dg["sd"], float)
    expected = np.asarray(dg["expected"], float)
    risk_raw = (prob[:, None, :] * np.abs(grid[None, :, None] - centres[None, None, :])).sum(2)

    out = dict(pdb=pdb, n=n, fold=fold, seq=seq, U=int(U), k=int(k), idx=idx,
               rr=np.asarray(u["rr"], float)[idx], i=i, j=j, sep=(j - i).astype(float),
               prob=prob, expected=expected, sd=sd, centres=centres, grid=grid,
               risk_raw=risk_raw, w_ship=_wship(sd))
    if "D" in want:
        W = np.asarray(u["W"], float)[idx]
        out["D"] = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1).astype(np.float32)
    if "W" in want:
        out["W"] = np.asarray(u["W"], float)[idx]
        out["nat_ca"] = np.asarray(u["nat_ca"], float)
    if "T" in want:
        out["PHI"] = np.asarray(u["PHI"], float)[idx]
        out["PSI"] = np.asarray(u["PSI"], float)[idx]
    if "org" in want:
        out["org"] = np.asarray(u["org"], bool)[idx]
    return out


def gbin(D, grid):
    """Grid index of observed distances, exactly as `instrument.shipped_score` does it."""
    return np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)


def shipped(p, D=None):
    """The shipped Bayes-risk score, recomputed from `prob`.  Lower is better."""
    D = p["D"] if D is None else D
    g = gbin(D, p["grid"])
    R = p["risk_raw"][np.arange(p["risk_raw"].shape[0])[None, :], g]
    return (R * p["w_ship"][None, :]).mean(1)


# --------------------------------------------------------------------------- selection
def sel_of(score, rr, tol=0.0):
    """Tie-safe realised RMSD of `argmin score`, averaged over the whole tied argmin set."""
    s = np.asarray(score, float)
    rr = np.asarray(rr, float)
    m = s.min()
    hit = np.where(s <= m + tol)[0]
    return float(rr[hit].mean())


def sel_idx(score, rr, tol=0.0):
    s = np.asarray(score, float)
    return np.where(s <= s.min() + tol)[0]


def rand_expect(rr):
    """Matched random selection of one: the EXACT expectation, noise-free."""
    return float(np.asarray(rr, float).mean())


def rand_draws(rr, pdb, tag, n=200):
    """Seeded finite-draw matched-random control, for its sampling distribution."""
    rng = SD.stable_rng(pdb, "s17sel", tag)
    rr = np.asarray(rr, float)
    return rr[rng.integers(0, len(rr), size=int(n))]


def topm_random(rr, m, pdb, tag, n=200):
    """Random selection of the BEST-of-m random draws -- the matched control for any
    selector that is allowed to look at m candidates.  Returns the mean over n repeats."""
    rng = SD.stable_rng(pdb, "s17selm", tag)
    rr = np.asarray(rr, float)
    m = min(int(m), len(rr))
    return float(np.mean([rr[rng.choice(len(rr), m, replace=False)].min() for _ in range(int(n))]))


# --------------------------------------------------------------------------- statistics
def boot_target(diff, B=8000, seed="s17target"):
    """Paired bootstrap over TARGETS.  `diff` is one number per target."""
    d = np.asarray(diff, float)
    n = len(d)
    rng = SD.stable_rng(seed, str(n))
    m = d[rng.integers(0, n, size=(int(B), n))].mean(1)
    return dict(mean=float(d.mean()), median=float(np.median(d)),
                lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)),
                W=int((d < 0).sum()), L=int((d > 0).sum()), T=int((d == 0).sum()),
                n=n, p_gain=float((m < 0).mean()))


def boot_fold(diff, fold, B=8000, seed="s17fold"):
    """Paired CLUSTER bootstrap: resample the 5 FOLDS with replacement.

    Five clusters is coarse and the interval is correspondingly wide.  That is the honest
    price of clustering at the level the leakage structure actually lives at, and it is
    reported beside the target-level interval rather than instead of it.
    """
    d = np.asarray(diff, float)
    f = np.asarray(fold, int)
    uf = np.unique(f)
    rng = SD.stable_rng(seed, str(len(d)))
    groups = [d[f == g] for g in uf]
    out = np.empty(int(B))
    for b in range(int(B)):
        pick = rng.integers(0, len(uf), size=len(uf))
        out[b] = np.concatenate([groups[q] for q in pick]).mean()
    return dict(mean=float(d.mean()), lo=float(np.percentile(out, 2.5)),
                hi=float(np.percentile(out, 97.5)), n_folds=int(len(uf)))


def drop_top(diff, k=10):
    """Mean after removing the k targets contributing the largest gains -- the
    concentration check.  A result carried by <10 targets is not a result."""
    d = np.sort(np.asarray(diff, float))
    return float(d[k:].mean()) if len(d) > k else float("nan")


def report_pair(name, a, b, fold, k_drop=10):
    """One canonical comparison line: a - b, negative = a better."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    d = a - b
    bt = boot_target(d)
    bf = boot_fold(d, fold)
    return dict(name=name, mean_a=float(a.mean()), mean_b=float(b.mean()),
                diff=bt["mean"], median=bt["median"], lo=bt["lo"], hi=bt["hi"],
                W=bt["W"], L=bt["L"], fold_lo=bf["lo"], fold_hi=bf["hi"],
                drop_top=drop_top(d, k_drop), n=int(len(d)))


def fmt_pair(r):
    return (f"{r['name']:<30} {r['mean_a']:>7.3f} vs {r['mean_b']:>7.3f}   "
            f"{r['diff']:+.3f} [{r['lo']:+.3f},{r['hi']:+.3f}]  "
            f"fold [{r['fold_lo']:+.3f},{r['fold_hi']:+.3f}]  "
            f"med {r['median']:+.3f}  {r['W']}W/{r['L']}L  dtop {r['drop_top']:+.3f}")


# --------------------------------------------------------------------------- LFO variants
def lfo_variant(M, fold, ref_col=0, lower_is_better=True):
    """Leave-fold-out selection AMONG variants.

    `M` is (n_targets, n_variants) of realised RMSD.  For each fold, pick the variant with
    the best mean on the OTHER folds, apply it to this fold.  Returns the held-out realised
    vector, the per-fold choices, and the in-sample best (the overfitting upper bound).

    This is the guard the `score-axis-does-not-transfer` post-mortem says was missing: a
    best-of-N over all targets is an in-sample number and must never be quoted as a result.
    """
    M = np.asarray(M, float)
    f = np.asarray(fold, int)
    uf = np.unique(f)
    out = np.empty(len(M))
    choice = {}
    for g in uf:
        tr = f != g
        v = int(np.argmin(M[tr].mean(0))) if lower_is_better else int(np.argmax(M[tr].mean(0)))
        out[f == g] = M[f == g, v]
        choice[int(g)] = v
    ins = int(np.argmin(M.mean(0)))
    return out, choice, ins


# --------------------------------------------------------------------------- selfcheck
def selfcheck(n=12, verbose=True):
    """Reproduce the pinned constants on a subset, and prove `shipped` == the instrument."""
    tg = targets()[:n]
    mx = 0.0
    best, sel = [], []
    dsel = 0.0
    for t in tg:
        p = pack(t["pdb"], 500)
        s0 = np.asarray(I.shipped_score(I.distogram(t["pdb"]), p["D"].astype(float)), float)
        s1 = shipped(p)
        mx = max(mx, float(np.abs(s0 - s1).max()))
        dsel = max(dsel, abs(sel_of(s0, p["rr"]) - sel_of(s1, p["rr"])))
        best.append(p["rr"].min())
        sel.append(sel_of(s1, p["rr"]))
    # tie trap: a constant score must return exactly the set mean
    p = pack(tg[0]["pdb"], 500)
    assert abs(sel_of(np.zeros(len(p["rr"])), p["rr"]) - p["rr"].mean()) < 1e-9
    if verbose:
        print(f"selfcheck n={n}: max |shipped - instrument| = {mx:.3e} (float32 cache), max selection deviation {dsel:.1e}")
        print(f"  pool best {np.mean(best):.4f}   shipped-selected {np.mean(sel):.4f}")
    assert mx < 1e-4 and dsel < 1e-9, "recomputed shipped score does not match the instrument"
    return dict(max_dev=mx, pool_best=float(np.mean(best)), sel=float(np.mean(sel)))


if __name__ == "__main__":
    selfcheck()
