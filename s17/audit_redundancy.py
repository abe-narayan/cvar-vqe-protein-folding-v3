"""s17/audit_redundancy.py -- IS THE DEEP UNIVERSE NEW STRUCTURE, OR MORE COPIES?

THE CLAIM UNDER ATTACK.  s17/BRIEF.md section 2 and s17/oracle_map.py present the full
window universe (13k-27k windows) as a wider CEILING than the K = 500 pool.  A ceiling
over 13,000 near-copies of 500 things is not a wider ceiling; it is the same ceiling with
a longer tail.  The universe is built by SLIDING every length-n frame along every library
chain (`s8/generate.py::_windows_all`), so consecutive windows share n-1 of n residues by
construction -- redundancy is the null hypothesis here, not an exotic possibility.

PRE-REGISTRATION.

  H0   The extra windows beyond K = 500 are structurally redundant with what is already
       there.  Specifically: (a) the number of DISTINCT conformations grows far more
       slowly than K, and (b) the full-universe oracle-best window is a near-duplicate of
       a window already in the K = 500 pool, so the ceiling "gain" is a refinement within
       an already-retrieved conformational family, not a new family.

  H1   The extra windows contain genuinely new conformations, and the full-universe best
       is structurally distinct from everything in the top 500.

  MEASUREMENTS
    (1) NEAR-DUPLICATE TEST, all 126 targets.  d* = CA-RMSD from the full-universe
        oracle-best window to its nearest neighbour inside the K = 500 pool.  Compare it
        with the ceiling gain it supposedly delivers.  H0 predicts d* small.
        Matched reference: the same nearest-neighbour distance for a RANDOM window drawn
        from outside the top 500 -- otherwise "small d*" could just mean the universe is
        dense everywhere.
    (2) DISTINCT-CONFORMATION COUNT.  Greedy leader clustering at 1.0 / 1.5 / 2.0 A on
        the BLOSUM-ordered universe: how many clusters exist by rank K.  If C(full) is a
        small multiple of C(500) while K grows 26x, the universe is redundant.
    (3) COVERAGE.  Fraction of full-universe clusters already represented at K = 500.

  SUCCESS CRITERION for H1: d* comparable to the ceiling gain AND cluster count growing
  roughly in proportion to K.
  FALSIFIER for H0: if d* is large (new family) this audit's redundancy story fails.

ORACLE.  rr and nat_ca are labels; every number here is a diagnostic.  benchmark60 is not
touched.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I     # noqa: E402
from s15 import seed as SD          # noqa: E402

UNIV = os.path.join(ROOT, "s8", "generate_univ")
THRESH = (1.5, 2.0)   # 1.0 A dropped: leader count ~9k makes the O(U x L) sweep unaffordable here
KS_C = (75, 500, 2000, 5000, 0)
N_CLUSTER_TARGETS = 6             # stratified subsample for the O(U x leaders) part
CHUNK = 256


def _cross_rmsd(A, B):
    """RMSD of every A[a] to every B[b] after optimal superposition -> (a, b)."""
    A = np.asarray(A, float); B = np.asarray(B, float)
    n = A.shape[1]
    Ac = A - A.mean(1, keepdims=True); Bc = B - B.mean(1, keepdims=True)
    na = (Ac ** 2).sum((1, 2)); nb = (Bc ** 2).sum((1, 2))
    H = np.einsum("ani,bnj->abij", Ac, Bc)
    S = np.linalg.svd(H, compute_uv=False)
    # reflection correction: if det(H) < 0 the smallest singular value flips sign
    det = np.linalg.det(H)
    S = S.copy()
    S[..., -1] = np.where(det < 0, -S[..., -1], S[..., -1])
    num = na[:, None] + nb[None, :] - 2.0 * S.sum(-1)
    return np.sqrt(np.maximum(num, 0.0) / n)


def _leaders(W, thr, order):
    """Greedy leader clustering in BLOSUM-rank order.  Returns leader ranks (ascending)."""
    lead_idx = []
    Wl = None
    U = len(order)
    for s in range(0, U, CHUNK):
        blk = order[s:s + CHUNK]
        Wb = W[blk]
        if Wl is None:
            lead_idx.append(0)
            Wl = Wb[0:1].copy()
            start = 1
        else:
            start = 0
        if start < len(Wb):
            D = _cross_rmsd(Wb[start:], Wl)
            unassigned = np.where(D.min(1) >= thr)[0]
        else:
            unassigned = np.array([], int)
        # unassigned members must still be checked against each other, sequentially
        for q in unassigned:
            w = Wb[start + q]
            if len(Wl) and _cross_rmsd(w[None], Wl).min() < thr:
                continue
            Wl = np.concatenate([Wl, w[None]], 0)
            lead_idx.append(s + start + q)
    return np.asarray(lead_idx, int)


def run_dupe(verbose=True):
    """(1) the near-duplicate test, all 126 targets."""
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        z = np.load(os.path.join(UNIV, f"{pdb}.npz"), allow_pickle=True)
        W = np.asarray(z["W"], float)
        rr = np.asarray(z["rr"], float)
        order = np.asarray(z["order"], int)
        U = len(rr)
        k = min(500, U)
        pool = order[:k]
        rest = order[k:]
        rng = SD.stable_rng(pdb, "s17auditdupe")

        i_full = int(np.argmin(rr))
        i_500 = int(pool[np.argmin(rr[pool])])
        gain = float(rr[i_500] - rr[i_full])

        dstar = float(_cross_rmsd(W[i_full][None], W[pool]).min())
        # matched reference: random windows from OUTSIDE the pool, same NN measurement
        if len(rest):
            pick = rng.choice(rest, size=min(50, len(rest)), replace=False)
            dref = _cross_rmsd(W[pick], W[pool]).min(1)
            dref_mean = float(dref.mean()); dref_med = float(np.median(dref))
        else:
            dref_mean = dref_med = float("nan")
        d_best_to_best = float(_cross_rmsd(W[i_full][None], W[i_500][None])[0, 0])

        rows.append({"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]), "U": int(U),
                     "rr_full": float(rr[i_full]), "rr_500": float(rr[i_500]),
                     "gain": gain, "rank_full": int(np.where(order == i_full)[0][0]),
                     "d_star": dstar, "d_best_to_best": d_best_to_best,
                     "d_ref_mean": dref_mean, "d_ref_med": dref_med})
        if verbose and (c + 1) % 25 == 0:
            print(f"  dupe {c+1}/{len(tg)}", flush=True)
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "audit_dupe.json"), "w"))
    return rows


def run_clusters(verbose=True):
    """(2)/(3) distinct-conformation counts on a stratified subsample."""
    tg = I.targets()
    ns = np.array([t["n"] for t in tg])
    rng = SD.stable_rng("s17audit", "clustersample")
    # stratify by chain length so U (which is ~1/n) is spanned
    idx = []
    for lo, hi in ((9, 11), (12, 13), (14, 16)):
        pool = np.where((ns >= lo) & (ns <= hi))[0]
        idx += list(rng.choice(pool, size=min(2, len(pool)), replace=False))
    idx = sorted(set(int(x) for x in idx))
    rows = []
    for c, k in enumerate(idx):
        t = tg[k]; pdb = t["pdb"]
        z = np.load(os.path.join(UNIV, f"{pdb}.npz"), allow_pickle=True)
        W = np.asarray(z["W"], float)
        order = np.asarray(z["order"], int)
        U = len(order)
        ent = {"pdb": pdb, "n": int(t["n"]), "U": int(U), "counts": {}}
        for thr in THRESH:
            t0 = time.time()
            L = _leaders(W, thr, order)
            ent["counts"][str(thr)] = {
                ("full" if K == 0 else str(K)): int((L < (U if K == 0 else min(K, U))).sum())
                for K in KS_C}
            if verbose:
                print(f"  {pdb} thr={thr} U={U} clusters={len(L)} ({time.time()-t0:.0f}s)",
                      flush=True)
        rows.append(ent)
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "audit_clusters.json"), "w"))
    return rows


# ------------------------------------------------------------------ reporting
def _boot(d, rng, B=4000, folds=None):
    d = np.asarray(d, float)
    if folds is None:
        m = d[rng.integers(0, len(d), size=(B, len(d)))].mean(1)
    else:
        folds = np.asarray(folds); uf = np.unique(folds)
        idx = [np.where(folds == f)[0] for f in uf]
        m = np.empty(B)
        for b in range(B):
            m[b] = d[np.concatenate([idx[q] for q in rng.integers(0, len(uf), len(uf))])].mean()
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report():
    rng = SD.stable_rng("s17audit", "redundancy-report")
    p = os.path.join(RESULTS, "audit_dupe.json")
    if os.path.exists(p):
        rows = json.load(open(p))["rows"]
        fold = np.array([r["fold"] for r in rows])
        ds = np.array([r["d_star"] for r in rows])
        db = np.array([r["d_best_to_best"] for r in rows])
        dr = np.array([r["d_ref_med"] for r in rows])
        gn = np.array([r["gain"] for r in rows])
        rk = np.array([r["rank_full"] for r in rows], float)
        print("\nA. IS THE FULL-UNIVERSE BEST A NEW STRUCTURE OR A NEAR-COPY?")
        m, lo, hi = _boot(ds, rng, folds=fold)
        print(f"  d* = RMSD of the full-universe best window to its NEAREST NEIGHBOUR in")
        print(f"       the K=500 pool          mean {m:.3f} [{lo:.3f},{hi:.3f}]  median {np.median(ds):.3f}")
        print(f"  RMSD of full-best to the K=500 BEST window   mean {db.mean():.3f}  median {np.median(db):.3f}")
        print(f"  ceiling gain it delivers                     mean {gn.mean():.3f}  median {np.median(gn):.3f}")
        m2, lo2, hi2 = _boot(dr, rng, folds=fold)
        print(f"  MATCHED REFERENCE: median NN-to-pool distance of RANDOM out-of-pool")
        print(f"       windows                 mean {m2:.3f} [{lo2:.3f},{hi2:.3f}]")
        md, dlo, dhi = _boot(ds - dr, rng, folds=fold)
        print(f"  d* - reference              {md:+.3f} [{dlo:+.3f},{dhi:+.3f}]   "
              f"W/L {int((ds<dr).sum())}/{int((ds>dr).sum())}")
        print(f"  BLOSUM rank of the full-universe best: median {np.median(rk):.0f}, "
              f"mean {rk.mean():.0f};  in top-500 on {int((rk<500).sum())}/{len(rk)} targets")
        print(f"  targets where the full-best is within 1.0 A of something already in the")
        print(f"       pool: {int((ds<1.0).sum())}/{len(ds)};  within 1.5 A: {int((ds<1.5).sum())}")
    q = os.path.join(RESULTS, "audit_clusters.json")
    if os.path.exists(q):
        rows = json.load(open(q))["rows"]
        print("\nB. DISTINCT CONFORMATIONS BY RANK  (greedy leader clustering, BLOSUM order)")
        for thr in THRESH:
            ks = [("full" if K == 0 else str(K)) for K in KS_C]
            cs = {k: np.array([r["counts"][str(thr)][k] for r in rows], float) for k in ks}
            Uv = np.array([r["U"] for r in rows], float)
            print(f"  threshold {thr} A   (n={len(rows)} targets)")
            print("      " + "".join(f"{k:>9}" for k in ks) + f"{'U':>9}")
            print("      " + "".join(f"{cs[k].mean():>9.1f}" for k in ks) + f"{Uv.mean():>9.0f}")
            print(f"      clusters per window at K=500 {cs['500'].mean()/500:.3f}   "
                  f"at full {(cs['full']/Uv).mean():.3f}   "
                  f"C(full)/C(500) = {(cs['full']/np.maximum(cs['500'],1)).mean():.2f}x "
                  f"while K grows {(Uv/500).mean():.1f}x")


if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else "all"
    if a in ("all", "dupe"):
        run_dupe()
    if a in ("all", "clusters"):
        run_clusters()
    report()
