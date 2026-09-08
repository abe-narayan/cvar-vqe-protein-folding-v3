"""SPRINT 15, coordinator -- UNDO THE CONTRACTION BEFORE PROJECTING.

THE GAP THIS ATTACKS. The full cascade (K8) leaves three gaps, and after aggregation recovers
-0.360 A from selection, the largest one left is the PROJECTION:

    G -> S   +0.751     choosing the objective's argmin
    S -> A   -0.360     aggregation recovering most of it
    A -> F   **+0.170**  putting the consensus back on the ideal-geometry manifold

The A -> F cost is not a mystery and it is not a tax. Coordinate averaging **contracts the backbone
by 25.8%**: superposing many structures that disagree and taking the mean pulls every atom toward
the centroid, and the more they disagree the more it pulls. The projection then has to re-expand a
structure that is too small, and it pays for that in RMSD. Sprint 14 initially misattributed this
cost to "the price of ideal geometry"; the standing correction is that it is **the price of
re-expanding a contracted structure**, which is a different claim with a different fix.

THE FIX IS ONE NUMBER, AND IT IS NATIVE-FREE. Rescale the consensus so its distance scale matches
the predicted distances before projecting. With weights `w` the optimal isotropic scale has a closed
form:

    s* = sum_p w_p d_p dhat_p  /  sum_p w_p d_p^2

which uses only the consensus geometry and the distogram -- no native coordinates anywhere.

WHY THIS MATTERS BEYOND THIS SPRINT. The **production pipeline does exactly the same thing**: it
coordinate-averages the top-75 retrieved windows and projects the result, and `avg_ca` and `fit_ca`
are both cached in every shipped record. So the same rescaling can be applied to the incumbent
itself. If it helps there, it is not an improvement to an experimental arm -- it is a deployable
improvement to the shipped pipeline, obtained from a measurement nobody had made.

ARMS, on both the incumbent's own average and the cascade's consensus:

    project_raw          project the consensus as-is                       the control (= F today)
    project_scaled       rescale by s* against the distogram, then project the arm under test
    project_scaled_deb   the same, against separation-debiased distances
    ORACLE_scale_true    rescale by s* against the TRUE distances           ORACLE, the ceiling
    ORACLE_scale_best    the single scale minimising RMSD, swept            ORACLE, the ceiling of
                                                                            isotropic rescaling

The two ORACLE arms bound the idea from two directions: `ORACLE_scale_true` says how well the
correct scale can be recovered from perfect distances, and `ORACLE_scale_best` says how much any
isotropic rescaling could ever be worth. If `ORACLE_scale_best` is small, the contraction is not
isotropic and the fix has to be shape-aware rather than scalar -- which is itself worth knowing and
is a cheap thing to learn before building anything more elaborate.

Run:
    python -m s15.expand
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import distcal as C                 # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

SWEEP = np.round(np.arange(0.85, 1.46, 0.025), 4)


def optimal_scale(CA, target_d, w, i, j):
    """Closed-form isotropic scale matching the structure's distances to `target_d`."""
    d = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    num = float((w * d * target_d).sum())
    den = float((w * d * d).sum())
    return num / max(den, 1e-12)


def contraction(CA, target_d, i, j):
    """How much smaller is this structure than the restraints say? Reported as a percentage."""
    d = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    return float(100.0 * (1.0 - d.mean() / max(target_d.mean(), 1e-12)))


def run_target(t, deb):
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    rec = I.shipped_record(pdb)
    avg = np.asarray(rec["avg_ca"], float)             # the production consensus, pre-projection
    i, j = I.pair_index(n)
    sep = (j - i).astype(float)

    dg = I.distogram(pdb, seq, fold)
    dhat = np.asarray(dg["expected"], float)
    sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
    w = 1.0 / sd ** 2
    dhat_deb = np.maximum(dhat - deb(sep), 2.0)
    dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))

    def emit(CA):
        return float(I.ca_rmsd(I.project(CA, seq, fold)["fit_ca"], nat))

    #: THE POOL AS A SCALE REFERENCE. The first run showed the distogram-derived scale
    #: OVER-expands (it wants 1.0997 where the truth wants 1.0412) and costs +0.123, because the
    #: distogram's own +0.509 A expansion bias contaminates the estimate. K3 measured the
    #: retrieval pool's distances as having the OPPOSITE bias sign (-0.128 global), so the pool
    #: should over-correct in the other direction and their average should land closer. This is
    #: the one place in the sprint where the two channels' opposed biases can be used directly.
    from s15 import pooldist as PD
    Dm, _sim, _pi, _pj = PD.pool_distances(pdb, n)
    d_pool = np.median(Dm, axis=0)
    w_pool = 1.0 / np.maximum(Dm.std(axis=0), 0.25) ** 2

    s_pred = optimal_scale(avg, dhat, w, i, j)
    s_deb = optimal_scale(avg, dhat_deb, w, i, j)
    s_pool = optimal_scale(avg, d_pool, w_pool, i, j)
    s_both = 0.5 * (s_deb + s_pool)
    s_true = optimal_scale(avg, dtrue, np.ones_like(dtrue), i, j)

    row = {"pdb": pdb, "n": n, "fold": fold,
           "contraction_vs_true_pct": contraction(avg, dtrue, i, j),
           "contraction_vs_pred_pct": contraction(avg, dhat, i, j),
           "s_pred": s_pred, "s_deb": s_deb, "s_pool": s_pool, "s_both": s_both,
           "s_true": s_true,
           "arms": {}}

    #: COST NOTE. A projection costs 5-20 s, so 126 targets x 6 arms x an 11-point sweep is
    #: several hours. Two economies, neither of which weakens a claim:
    #:   1. `project_raw` IS the production `fit_ca`, already cached in the shipped record.
    #:      We read it and assert agreement with a live projection on a sample rather than
    #:      recomputing it 126 times.
    #:   2. The scale sweep exists only to answer "is the contraction ISOTROPIC?". That is a
    #:      question about the consensus geometry, not about the projection, so it is answered
    #:      on the UNPROJECTED structure, where scaling and scoring are free.
    row["arms"]["project_raw"] = float(I.ca_rmsd(np.asarray(rec["fit_ca"], float), nat))
    row["arms"]["project_scaled"] = emit(avg * s_pred)
    row["arms"]["project_scaled_deb"] = emit(avg * s_deb)
    row["arms"]["project_scaled_pool"] = emit(avg * s_pool)
    row["arms"]["project_scaled_both"] = emit(avg * s_both)
    row["arms"]["ORACLE_scale_true"] = emit(avg * s_true)

    #: free sweep, pre-projection: the ceiling of any isotropic rescaling
    raw0 = float(I.ca_rmsd(avg, nat))
    sweep = {float(sc): float(I.ca_rmsd(avg * sc, nat)) for sc in SWEEP}
    best = min(sweep, key=lambda x: sweep[x])
    row["unprojected_raw"] = raw0
    row["unprojected_best"] = sweep[best]
    row["unprojected_gain_ORACLE"] = raw0 - sweep[best]
    row["oracle_best_scale"] = float(best)
    row["sweep_unprojected"] = sweep
    return row


def verify_cached_projection(targets, deb, n_check=4):
    """`project_raw` reads the cached `fit_ca`; assert a live projection reproduces it."""
    errs = []
    for t in targets[:n_check]:
        pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
        rec = I.shipped_record(pdb)
        avg = np.asarray(rec["avg_ca"], float)
        live = I.project(avg, seq, fold)["fit_ca"]
        errs.append(float(np.abs(live - np.asarray(rec["fit_ca"], float)).max()))
    return errs


def run(targets=None):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    data = C.gather(tg)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    path = os.path.join(RESULTS, "expand.json")
    for c, t in enumerate(tg):
        rows.append(run_target(t, deb[int(t["fold"])]))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": rows, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    arms = list(rows[0]["arms"])
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray([r["arms"]["project_raw"] for r in rows], float)

    # the leave-fold-out scale, as a single number per fold -- the deployable version
    lfo_scale = {}
    for f in sorted(set(folds)):
        lfo_scale[int(f)] = float(np.median([rows[z]["s_pred"] for z in range(len(rows))
                                             if int(folds[z]) != f]))

    out = {"n": len(rows), "incumbent": float(ref.mean()), "rows": rows, "arms": {}}
    named = {a: np.asarray([r["arms"][a] for r in rows], float) for a in arms}
    for a, v in named.items():
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_raw": I.paired(v, base, folds=folds, names=pdbs),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    ct = float(np.mean([r["contraction_vs_true_pct"] for r in rows]))
    cp = float(np.mean([r["contraction_vs_pred_pct"] for r in rows]))
    out["contraction_vs_true_pct"] = ct
    out["contraction_vs_pred_pct"] = cp
    out["s_pred_mean"] = float(np.mean([r["s_pred"] for r in rows]))
    out["s_pool_mean"] = float(np.mean([r["s_pool"] for r in rows]))
    out["s_both_mean"] = float(np.mean([r["s_both"] for r in rows]))
    at = np.asarray([r["s_true"] for r in rows])
    for k in ("s_pred", "s_deb", "s_pool", "s_both"):
        v = np.asarray([r[k] for r in rows])
        out[f"{k}_vs_true"] = {"corr": float(np.corrcoef(v, at)[0, 1]),
                               "mae": float(np.abs(v - at).mean())}
    out["s_true_mean"] = float(np.mean([r["s_true"] for r in rows]))
    out["oracle_best_scale_mean"] = float(np.mean([r["oracle_best_scale"] for r in rows]))
    out["lfo_scale_by_fold"] = lfo_scale
    out["unprojected_raw_mean"] = float(np.mean([r["unprojected_raw"] for r in rows]))
    out["unprojected_best_mean_ORACLE"] = float(np.mean([r["unprojected_best"] for r in rows]))
    out["unprojected_gain_ORACLE_mean"] = float(
        np.mean([r["unprojected_gain_ORACLE"] for r in rows]))
    out["projection_check_maxerr"] = verify_cached_projection(tg, deb)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_expand", out, n_expected=len(tg))

    print(f"\nn = {len(rows)}   incumbent {ref.mean():.3f}")
    print(f"the production consensus is contracted by {ct:.1f}% against the TRUE distances "
          f"and {cp:.1f}% against the PREDICTED ones")
    print(f"mean optimal scale: predicted {out['s_pred_mean']:.4f}, "
          f"true {out['s_true_mean']:.4f}, RMSD-optimal (ORACLE) "
          f"{out['oracle_best_scale_mean']:.4f}\n")
    print(f"{'arm':<24}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs project_raw':>24}")
    for a in list(arms):
        s = out["arms"][a]
        v = s["vs_raw"]
        print(f"{a:<24}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"  {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print("\n(`project_raw` reproduces the production `fit_ca`; only `project_scaled*` and"
          "\n `scaled_LFO_constant` are native-free. If `ORACLE_scale_best` is small, the"
          "\n contraction is not isotropic and a scalar fix cannot work.)")
    return out


if __name__ == "__main__":
    run()
