"""s23/gscale.py -- ONE GLOBAL SCALE CONSTANT, FITTED BY NESTED CV.

WHY THIS ARM EXISTS.  The incumbent emits a COORDINATE AVERAGE of the top-75 pool members, and
averaging CONTRACTS: mean virtual Ca-Ca bond 3.12 A against a native 3.82 A.  Kabsch fits rotation
and translation only -- IT DOES NOT FIT SCALE -- so any systematic scale error propagates directly
into Ca-RMSD.  Coordinator scoping at n=40 measured:

    ORACLE per-target scale     -0.158 A [-0.255,-0.085]   38W/2L   mean s* = 0.975

and, counter-intuitively, that every physically-motivated estimator points the WRONG WAY:
matching the distogram-predicted Rg wants s = 1.084 (+0.115 A WORSE), a least-squares fit of the
pair distances to `dhat` wants s = 1.079 (+0.101 WORSE), and restoring the physical bond length
(s = 1.277) costs +0.82 A.

THE GEOMETRIC READING, WHICH IS THE POINT.  s* < 1 while the BOND is 18% short means the defect is
NOT a uniform shrink.  Averaging smooths: it shortens local bonds while leaving the global envelope
roughly right, and can leave it marginally too large.  A single global scale cannot fix that -- it
trades one error for the other -- which is precisely why the local-geometry route (torsion-space
rebuild, Ca-preserving repair) is a different lane and a better-motivated one.  This file exists to
PRICE the scale route honestly, not to rescue it.

WHY A GLOBAL CONSTANT AND NOT A PER-TARGET RULE.  A finite-sample bound derived last sprint from
this project's own sigma ~ 0.41 A shows that at n ~ 100 per training fold, even a ONE-GLOBAL-
THRESHOLD router carries a generalisation gap of 0.39 A, and richer classes need 280-1045 targets
against a benchmark of 126.  A single scalar is the ONLY adaptive object this sample size supports,
and four per-target routers failed last sprint exactly as the bound predicts.

OPERATOR FORKS, per BRIEF SS4 rule 0.  Directional hypothesis, so each fork names its alternative.

    functional     DECLARED the shipped Bayes-risk distogram score for the filter, as deployed.
                   NOT TAKEN the squared functional.
    basis          DECLARED the pool's own window coordinates, and the arm is a POINT CLOUD on both
                   sides -- the incumbent is a point cloud, so the comparison is like for like.
                   NOT TAKEN the rebuilt chain, which is a different basis and belongs to another lane.
    readout        DECLARED top-75 coordinate average then a global scale about the CENTROID.
                   NOT TAKEN scaling about the medoid frame (reported as an audit arm) or scaling
                   before the average rather than after (incoherent -- scaling each member first is
                   equivalent to scaling the mean, since the average is linear).
    normalisation  DECLARED scale as a dimensionless multiplier on centred coordinates.
                   NOT TAKEN matching a target statistic (Rg, bond length) -- both MEASURED and both
                   WORSE, reported here rather than quietly dropped.
    null           DECLARED s = 1.0, the incumbent, on identical data.  The honest null: it asks
                   whether scaling beats NOT scaling.  NOT TAKEN the per-target oracle, which is the
                   ceiling under test and would beg the question.
    THE LABEL      DECLARED continuous Ca-RMSD.  NOT TAKEN any binarised "did it win", whose
                   threshold would depend on target difficulty.

  Hypothesis   a single global s < 1 captures a small but real fraction of the 0.158 A oracle.
  Falsifier    if the nested-CV arm fails to beat s = 1.0 past its own MDE with a fold-clustered CI
               excluding zero, the global-scale route is refuted and the residual is per-target,
               i.e. unreachable at this n.
  Null         s = 1.0.
  Nesting      s is fitted on the 4 training folds and applied to the held-out fold.  A constant
               fitted on the folds it is scored on is leakage and is reported separately as
               `s_insample` so the optimism is visible rather than hidden.
  Promotion    none without an ablation and a geometry report -- an arm that wins by contracting
               further is buying Ca-RMSD with geometry and that trade must be surfaced.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402

TOPM = 75
GRID = np.linspace(0.90, 1.10, 81)        # 0.0025 resolution around the observed s* = 0.975
CA_TRANS = 3.8046


def _save(o, name="gscale.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, grid %d points" % (len(tg), len(GRID)), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        C, b_med = I.coordinate_average(W[o[:TOPM]])
        C = np.asarray(C, float)
        cen = C - C.mean(0)
        #: the whole RMSD curve in s -- one row per target, so every fit below is free
        curve = np.array([I.ca_rmsd(cen * s, nat) for s in GRID], float)
        rows.append({
            "pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
            "curve": curve.tolist(),
            "rmsd_1": float(I.ca_rmsd(cen, nat)),
            "s_star": float(GRID[int(np.argmin(curve))]),
            "rmsd_star": float(curve.min()),
            "vbond": float(np.linalg.norm(C[1:] - C[:-1], axis=1).mean()),
            "vbond_nat": float(np.linalg.norm(nat[1:] - nat[:-1], axis=1).mean()),
            "rg": float(np.sqrt((cen ** 2).sum() / len(cen))),
            "rg_nat": float(np.sqrt(((nat - nat.mean(0)) ** 2).sum() / len(nat))),
        })
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("curve", "rmsd_1", "s_star", "vbond", "rg")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "grid": GRID.tolist(), "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "gscale.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    Cv = np.array([r["curve"] for r in rows], float)          # (T, G)
    fold = g("fold"); base = g("rmsd_1")
    rng = SD.stable_rng("gscale", "rep")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        b = x[rng.integers(0, k, size=(4000, k))].mean(1)
        F = sorted(set(fold.astype(int)))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)),
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    #: --- NESTED: fit s on the 4 training folds, apply to the held-out fold
    held = np.empty(len(rows)); s_used = np.empty(len(rows))
    for f in sorted(set(fold.astype(int))):
        tr = fold != f; te = fold == f
        s_hat = GRID[int(np.argmin(Cv[tr].mean(0)))]
        held[te] = Cv[te][:, int(np.argmin(Cv[tr].mean(0)))]
        s_used[te] = s_hat
    #: --- IN-SAMPLE: the same fit scored on the folds that chose it (the optimism, shown not hidden)
    s_all = GRID[int(np.argmin(Cv.mean(0)))]
    insample = Cv[:, int(np.argmin(Cv.mean(0)))]

    print("\nn = %d.  All arms POINT CLOUDS on the pool's window basis.\n" % len(rows))
    print("  geometry of the incumbent average:")
    print("    mean virtual Ca-Ca bond  %.3f A   (native %.3f A)  -> %.0f%% short"
          % (g("vbond").mean(), g("vbond_nat").mean(),
             100 * (1 - g("vbond").mean() / g("vbond_nat").mean())))
    print("    radius of gyration       %.3f A   (native %.3f A)"
          % (g("rg").mean(), g("rg_nat").mean()))
    print("\n  %-42s%9s" % ("arm", "RMSD"))
    print("  %-42s%9.4f" % ("incumbent, s = 1.0", base.mean()))
    print("  %-42s%9.4f   s = %.4f" % ("GLOBAL s, nested CV (held out)", held.mean(), s_used.mean()))
    print("  %-42s%9.4f   s = %.4f  <- leaky, for reference"
          % ("global s, in-sample", insample.mean(), s_all))
    print("  %-42s%9.4f   <- ORACLE per-target, CEILING" % ("per-target s*", g("rmsd_star").mean()))
    print("\n  s* distribution across targets: mean %.4f  sd %.4f  10th %.4f  90th %.4f"
          % (g("s_star").mean(), g("s_star").std(), np.percentile(g("s_star"), 10),
             np.percentile(g("s_star"), 90)))

    print("\n  PRIMARY (paired vs the incumbent):")
    for lab, x in (("GLOBAL s, NESTED CV  <-- the deployable arm", held - base),
                   ("global s, in-sample (leaky)", insample - base),
                   ("per-target s* (ORACLE ceiling)", g("rmsd_star") - base)):
        m, se, mde, lo, hi, flo, fhi, w = st(x)
        verdict = "BEATS" if (fhi < 0 and abs(m) > mde) else "not measured"
        print("    %-42s%+.4f SE %.4f MDE %.3f iid[%+.4f,%+.4f] fold[%+.4f,%+.4f] %3dW/%3dL  %s"
              % (lab, m, se, mde, lo, hi, flo, fhi, w, len(rows) - w, verdict))

    worst = (held - base).max()
    print("\n  worst single-target degradation under the nested arm: %+.3f A" % worst)
    print("  Falsifier: the nested arm failing to beat s=1.0 past its own MDE with a fold CI"
          " excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
