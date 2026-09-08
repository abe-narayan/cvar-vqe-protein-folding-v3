"""s17/audit_frame.py -- THE FUNCTIONAL CHECK ON THE ORACLE MAP, AND AN EXACT RANDOM CONTROL.

TWO INDEPENDENT ATTACKS.

F.  THE FRAME QUESTION (BRIEF section 5, "free-superposition member error read as the
    common-frame term the identity consumes; attribution wrong by ~2x").

    s17/oracle_map.py table A prints, on the same row:
        "mean"      = mean of I.kabsch_rmsd_batch(W, nat)  -- FREE superposition, each
                      window separately optimally rotated onto the native
        "diversity" = RMS deviation of the members about the coordinate average, in the
                      medoid frame the averaging operator builds -- a COMMON-frame term
        "coordavg"  = RMSD of the average to the native
    The sprint's central identity is  readout^2 = mean member error^2 - diversity^2.
    That identity is the parallel-axis theorem and it holds ONLY when the member errors
    are measured in the SAME fixed frame as the readout.  If it is evaluated with the
    free-superposition column that sits next to it in the table, it is being fed the wrong
    functional -- the exact error Sprint 16 paid ~2x for.

    HYPOTHESIS: the identity is exact with common-frame member errors and violated with
    the free-superposition column the map prints.
    FALSIFIER: if the two columns agree to within noise, the objection is empty.
    CONTROL: report the identity residual for BOTH functionals, at two K, per target.

C.  AN EXACT RANDOM-SELECTION CONTROL (BRIEF section 4, "every selector priced against a
    matched random selection of the same count").

    oracle_map.py estimates that control with N_RAND = 3 draws per target.  But the
    expectation of a uniform random pick from the top-K is available IN CLOSED FORM: it is
    the mean of rr over the top-K.  This module computes it exactly, at zero noise, for
    every K and every target, and re-prices the map's table D against it -- including the
    question the map's "recovered %" does not control for: how much of the realized change
    from K=75 to K is simply the pool's composition drifting, which a random selector
    experiences too.

Both parts read rr / nat_ca, which are ORACLE labels.  benchmark60 untouched.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I     # noqa: E402
from s15 import seed as SD          # noqa: E402

UNIV = os.path.join(ROOT, "s8", "generate_univ")
KS = (25, 50, 75, 100, 150, 300, 500, 1000, 2000, 5000, 0)
N_FRAME_TARGETS = 30
K_FRAME = (75, 500)


def _kname(K):
    return "full" if K == 0 else str(K)


# --------------------------------------------------------------- part F
def frame_check(verbose=True):
    tg = I.targets()
    rng = SD.stable_rng("s17audit", "framesample")
    pick = sorted(rng.choice(len(tg), size=min(N_FRAME_TARGETS, len(tg)), replace=False))
    rows = []
    for c, q in enumerate(pick):
        t = tg[int(q)]
        z = np.load(os.path.join(UNIV, f"{t['pdb']}.npz"), allow_pickle=True)
        W = np.asarray(z["W"], float); order = np.asarray(z["order"], int)
        nat = np.asarray(z["nat_ca"], float)
        ent = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "K": {}}
        for K in K_FRAME:
            Ws = W[order[:K]]
            P = I.pairwise_rmsd(Ws)
            avg, b = I.coordinate_average(Ws, P)
            Wm = I.superpose_batch(Ws, Ws[b])          # the frame the operator builds in
            n = Ws.shape[1]
            div = float(np.sqrt(np.mean(((Wm - avg) ** 2).sum(axis=(1, 2)) / n)))
            readout = I.ca_rmsd(avg, nat)

            # COMMON FRAME: rotate the whole ensemble (avg + members) rigidly onto the
            # native using the transform that superposes the AVERAGE.  Members are then
            # NOT re-superposed -- that is what the identity consumes.
            AV = np.concatenate([avg[None], Wm], 0)
            AVm = AV - avg.mean(0)                       # centre on the average
            natc = nat - nat.mean(0)
            H = AVm[0].T @ natc
            Uu, S, Vt = np.linalg.svd(H)
            d = np.sign(np.linalg.det(Vt.T @ Uu.T))
            D = np.diag([1.0, 1.0, d])
            R = Vt.T @ D @ Uu.T
            AVr = np.einsum("ij,bnj->bni", R, AVm)
            err_common = np.sqrt(((AVr[1:] - natc) ** 2).sum((1, 2)) / n)
            err_free = I.kabsch_rmsd_batch(Ws, nat)      # the column the map prints

            ent["K"][str(K)] = {
                "readout": float(readout), "div": div,
                "mean_err_free": float(err_free.mean()),
                "rms_err_free": float(np.sqrt((err_free ** 2).mean())),
                "mean_err_common": float(err_common.mean()),
                "rms_err_common": float(np.sqrt((err_common ** 2).mean())),
                # identity residuals, in A^2
                "resid_common": float((err_common ** 2).mean() - div ** 2 - readout ** 2),
                "resid_free_rms": float((err_free ** 2).mean() - div ** 2 - readout ** 2),
                "resid_free_arith": float(err_free.mean() ** 2 - div ** 2 - readout ** 2),
            }
        rows.append(ent)
        if verbose:
            print(f"  frame {c+1}/{len(pick)} {t['pdb']}", flush=True)
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "audit_frame.json"), "w"))
    return rows


# --------------------------------------------------------------- part C
def exact_random(verbose=True):
    tg = I.targets()
    rows = []
    for t in tg:
        z = np.load(os.path.join(UNIV, f"{t['pdb']}.npz"), allow_pickle=True)
        rr = np.asarray(z["rr"], float); order = np.asarray(z["order"], int)
        v = rr[order]
        U = len(v)
        e = {}
        for K in KS:
            k = U if K == 0 else min(K, U)
            e[_kname(K)] = float(v[:k].mean())
        rows.append({"pdb": t["pdb"], "fold": int(t["fold"]), "exact_rand": e})
    json.dump({"rows": rows}, open(os.path.join(RESULTS, "audit_exactrand.json"), "w"))
    return rows


# --------------------------------------------------------------- reporting
def _boot(d, rng, B=4000, folds=None):
    d = np.asarray(d, float)
    if folds is None:
        m = d[rng.integers(0, len(d), size=(B, len(d)))].mean(1)
    else:
        folds = np.asarray(folds); uf = np.unique(folds)
        idx = [np.where(folds == f)[0] for f in uf]
        m = np.array([d[np.concatenate([idx[q] for q in rng.integers(0, len(uf), len(uf))])].mean()
                      for _ in range(B)])
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report():
    rng = SD.stable_rng("s17audit", "frame-report")
    p = os.path.join(RESULTS, "audit_frame.json")
    if os.path.exists(p):
        rows = json.load(open(p))["rows"]
        print(f"\nF. THE IDENTITY  readout^2 = <member error^2> - diversity^2 "
              f"(n = {len(rows)} targets)")
        print("   residual in A^2; 'common' uses the frame the readout is measured in,")
        print("   'free' uses the free-superposition column oracle_map.py table A prints.")
        print(f"  {'K':>5}{'readout':>9}{'div':>8}{'err_common':>12}{'err_free':>10}"
              f"{'resid common':>14}{'resid free':>12}")
        for K in K_FRAME:
            g = lambda f: np.array([r["K"][str(K)][f] for r in rows])  # noqa: E731
            print(f"  {K:>5}{g('readout').mean():>9.3f}{g('div').mean():>8.3f}"
                  f"{g('rms_err_common').mean():>12.3f}{g('rms_err_free').mean():>10.3f}"
                  f"{g('resid_common').mean():>14.2e}{g('resid_free_rms').mean():>12.3f}")
        for K in K_FRAME:
            g = lambda f: np.array([r["K"][str(K)][f] for r in rows])  # noqa: E731
            rc, rf = np.abs(g("resid_common")), np.abs(g("resid_free_rms"))
            ec, ef = g("rms_err_common"), g("rms_err_free")
            print(f"  K={K}: max |resid common| {rc.max():.2e} A^2 (EXACT identity); "
                  f"max |resid free| {rf.max():.3f} A^2")
            print(f"        the free column understates the common-frame member error by "
                  f"{(ec - ef).mean():.3f} A on average ({(ec/ef).mean():.2f}x)")
            fa = np.abs(g("resid_free_arith"))
            print(f"        using the ARITHMETIC mean of the free column instead of its "
                  f"quadratic mean: max |resid| {fa.max():.3f} A^2")

    q = os.path.join(RESULTS, "audit_exactrand.json")
    m = os.path.join(RESULTS, "oracle_map.json")
    if os.path.exists(q) and os.path.exists(m):
        er = {r["pdb"]: r for r in json.load(open(q))["rows"]}
        mp = json.load(open(m))["rows"]
        pdbs = [r["pdb"] for r in mp]
        fold = np.array([er[p]["fold"] for p in pdbs])
        ks = [_kname(K) for K in KS]
        print(f"\nC. EXACT vs 3-DRAW RANDOM CONTROL  (map rows available: {len(mp)}/126)")
        print(f"  {'K':>6}{'map sel_rand (3 draws)':>24}{'EXACT E[random]':>18}{'diff':>9}")
        for k in ks:
            a = np.array([r["cells"][k]["sel_rand"] for r in mp])
            b = np.array([er[p]["exact_rand"][k] for p in pdbs])
            print(f"  {k:>6}{a.mean():>24.3f}{b.mean():>18.3f}{a.mean()-b.mean():>+9.3f}")
        print("\n  Re-pricing the map's table D against the EXACT control.")
        print("  'realized' = sel_dist(K) - sel_dist(75).  'random drift' is the same")
        print("  quantity for a random selector -- pool composition, not conversion.")
        s75 = np.array([r["cells"]["75"]["sel_dist"] for r in mp])
        b75 = np.array([r["cells"]["75"]["best"] for r in mp])
        r75 = np.array([er[p]["exact_rand"]["75"] for p in pdbs])
        print(f"  {'K':>6}{'ceiling d':>11}{'realized d':>12}{'  [95% CI]':>22}"
              f"{'random drift':>14}{'net of drift':>14}")
        for k in ("150", "300", "500", "1000", "2000", "5000", "full"):
            s = np.array([r["cells"][k]["sel_dist"] for r in mp])
            b = np.array([r["cells"][k]["best"] for r in mp])
            rd = np.array([er[p]["exact_rand"][k] for p in pdbs])
            mu, lo, hi = _boot(s - s75, rng, folds=fold)
            dm, _, _ = _boot(rd - r75, rng, folds=fold)
            nm, nlo, nhi = _boot((s - s75) - (rd - r75), rng, folds=fold)
            print(f"  {k:>6}{b.mean()-b75.mean():>+11.3f}{mu:>+12.3f}"
                  f"   [{lo:+.3f},{hi:+.3f}]{dm:>+14.3f}"
                  f"   {nm:+.3f} [{nlo:+.3f},{nhi:+.3f}]")


if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else "all"
    if a in ("all", "frame"):
        frame_check()
    if a in ("all", "rand"):
        exact_random()
    report()
