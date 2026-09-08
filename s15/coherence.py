"""SPRINT 15, coordinator -- REALIZABILITY, its null, and the fusion law, on the full instrument.

===============================================================================================
SPRINT 16 CORRECTION, 2026-09-06, RETRACT workstream.  READ BEFORE QUOTING ANY FUSION NUMBER.
===============================================================================================
Three defects in the fusion half of this module.  Nothing in the REALIZABILITY half is affected.

(1) THE LAW IS NOT OURS AND IS NOT NEW.  `d_avg = sqrt(r^2 - (s/2)^2)` is the two-member case of
    the KROGH-VEDELSBY (1995) AMBIGUITY DECOMPOSITION -- ensemble error = average member error
    minus diversity -- standard ensemble-learning material for 31 years, with a structural-biology
    twin in the ensemble-RMSD/RMSF identity.  `s16/lit_FINDINGS.md` B.5.  It must be CITED, not
    claimed.  The docstring below calls it "a parameter-free geometric identity"; that description
    is accurate but the word "identity" was doing no work, because it was never treated as one.

(2) IT IS COMPUTED WITH THE WRONG MEAN.  Line ~123 below uses `r_bar = 0.5*(rd + rp)`, the
    ARITHMETIC mean.  The identity requires the QUADRATIC mean, sqrt((rd^2 + rp^2)/2).  On these
    126 targets the two channels' RMSDs differ by mean 0.980 A (median 0.567, max 4.798), so this
    is not a small substitution.  Correcting it moves the reported prediction error from
    +0.162 A to +0.075 A -- **54% of the published "prediction error" was an arithmetic
    artefact, not physics**: -0.087 [-0.119, -0.061] i.i.d.-target, [-0.129, -0.056]
    fold-clustered, W/L 126/0.  `s16/retract_law.py`, `s16/results/retract_law.json`.

(3) THE CAVEATS BELOW ARE UNNECESSARY, AND THE REMAINING RESIDUAL IS A FRAME CONVENTION, NOT
    PHYSICS.  The identity is EXACT in any inner-product space -- no coplanarity, no equidistance.
    `s16/retract_exact.py` reproduces this module's Xd/Xp bit-identically (max |delta| 0.000e+00
    on all four quantities) and shows that with BOTH structures superposed onto the native and the
    quadratic mean, the identity holds to |err| mean 2.65e-15, max 2.52e-14, sign balanced 55/61.
    The residual +0.075 A that survives the mean fix is entirely a frame mismatch:
      +0.174 A  `I.coordinate_average` averages in the MEDOID frame, not the target's frame;
      -0.099 A  `s` is Kabsch-MINIMISED (3.138 A) where the identity wants the disagreement in
                that same common frame (3.506 A);
      = +0.075 A, exactly.
    The final Kabsch of the average onto the native contributes EXACTLY ZERO (mean = median =
    max = 0.000000 over 126 targets): the average of structures each optimally superposed onto T
    is itself optimally superposed onto T, because the cross-covariance is linear in the moving
    structure and a sum of symmetric PSD matrices is symmetric PSD.

    **CONSEQUENCE.** "The law holds through Kabsch superposition on 126 real targets" is a
    THEOREM, not an empirical finding, and this module never measured it -- it measured a
    different frame convention.  Full account: `s16/retract_FINDINGS.md`.

NOTE ON THE ARTEFACT.  `s15/results/coherence.json` on disk was written by an EARLIER version of
this module (it carries the key `coordavg_vs_best_single`, which the mislabelling correction below
renamed).  Its per-target `rows` are correct and are what Sprint 16 recomputed from; its aggregate
`law_prediction` / `prediction_error` fields are SUPERSEDED by `s16/results/retract_law.json`.
This module has NOT been re-run; the `law_prediction_quadmean` field added below will only appear
in a future run.
===============================================================================================

WHY THIS MODULE EXISTS. Three results were established on 24-target smoke reads and they carry more
of the sprint's interpretation than any other measurements:

  K9   about 56% of the distogram's error is geometrically REALIZABLE -- the predicted distance
       matrix is a good description of a WRONG structure rather than a noisy version of the right
       one, so the fit follows it instead of averaging it away;
  K9b  the null control: the same fit leaves 2.058 A on i.i.d. noise of the same RMS against 0.879 A
       on the real prediction, a ratio of 0.427, so the low residual is a property of the prediction
       and not of the fit's flexibility;
  K10  two channels can have raw errors correlated at +0.627 and still build wrong structures 2.980 A
       apart, and the gain from combining them is predicted to 0.135 A by a parameter-free geometric
       identity, `d_avg = sqrt(r^2 - (s/2)^2)`, in which **s is native-free**.

A hostile reviewer's sharpest criticism of this programme is that headline effects rest on n = 6,
n = 8 and n = 24, in a project whose own history includes an 8-target read that put an arm at 2.792 A
where the full instrument gave 3.644 A and **reversed the sign of the comparison**. That criticism is
correct, and this module is the answer to it: every one of the three results above, on all 126
targets, in one artefact, with paired confidence intervals.

WHAT IS MEASURED, PER TARGET

  the two fits          fit torsions to the distogram's distances and to the retrieval pool's median
                        distances, separately, by the same optimiser from the same starts
  realizability         |d_fit - d_target| at convergence -- how far each channel's distance set is
                        from being achievable by ANY ideal-geometry chain
  the null ladder       the same quantity for i.i.d. noise of matched RMS, for the real errors
                        permuted across pairs, and for the real magnitudes with random signs. These
                        are ORACLE DIAGNOSTICS: destroying a property of the error requires knowing
                        the error.
  disagreement          RMSD between the two fitted structures -- NATIVE-FREE, and the quantity the
                        fusion law consumes
  the law               predicted `sqrt(r^2 - (s/2)^2)` against the measured coordinate average and
                        against restraint-level fusion

WHAT WOULD FALSIFY EACH. Realizability fails if the real prediction's residual is not materially
below the i.i.d. null -- then the low residual was the fit's flexibility and K9 collapses. The
fusion law fails if its prediction error grows with n or is comparable to the effect it predicts;
it has **no free parameters**, so it cannot be rescued by refitting, which is exactly why it is worth
testing at scale.

Run:
    python -m s15.coherence
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
from s15 import distgeo as D                 # noqa: E402
from s15 import pooldist as P                # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)


def _fit(d_target, w, d, starts):
    """Objective-selected multi-start fit. Returns the CA trace and its own restraint residual."""
    i, j = d["i"], d["j"]
    best = None
    for phi0, psi0, _tag in starts:
        phi, psi, f, _ = D.fit_distances(np.maximum(d_target, 2.0), w, i, j, phi0, psi0)
        if best is None or f < best[0]:
            best = (f, phi, psi)
    CA = np.asarray(I.build_ca(best[1], best[2]), float)
    dfit = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    return CA, float(np.abs(dfit - d_target).mean())


def run(targets=None, n_start=4):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    data = C.gather(tg)

    rows = []
    path = os.path.join(RESULTS, "coherence.json")
    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, n = d["i"], d["j"], d["n"]
        dtrue, dhat, sd = d["dtrue"], d["dhat"], d["sd"]
        e = dhat - dtrue
        rms = float(np.sqrt((e * e).mean()))
        wd = 1.0 / sd ** 2
        starts = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p, "coherence"))
        rng = SD.stable_rng(p, "coherence_null")

        Dm, _s, _i, _j = P.pool_distances(p, n)
        d_pool = np.median(Dm, axis=0)
        wp = 1.0 / np.maximum(Dm.std(axis=0), 0.25) ** 2

        Xd, res_d = _fit(dhat, wd, d, starts)
        Xp, res_p = _fit(d_pool, wp, d, starts)
        _Xi, res_iid = _fit(dtrue + rng.normal(0, rms, size=len(e)), wd, d, starts)
        _Xm, res_perm = _fit(dtrue + e[rng.permutation(len(e))], wd, d, starts)
        _Xs, res_sign = _fit(dtrue + np.abs(e) * rng.choice([-1.0, 1.0], size=len(e)),
                             wd, d, starts)

        rd = float(I.ca_rmsd(Xd, d["nat"]))
        rp = float(I.ca_rmsd(Xp, d["nat"]))
        s_dis = float(I.ca_rmsd(Xd, Xp))
        Cm, _ = I.coordinate_average(np.stack([Xd, Xp]))
        r_avg = float(I.ca_rmsd(Cm, d["nat"]))
        wsum = wd + wp
        Xf, _ = _fit((wd * dhat + wp * d_pool) / wsum, wsum, d, starts)
        r_fuse = float(I.ca_rmsd(Xf, d["nat"]))

        #: SPRINT 16 CORRECTION (2026-09-06, RETRACT).  `r_bar` is the ARITHMETIC mean and is
        #: WRONG for this identity, which needs the QUADRATIC mean.  The original line is kept,
        #: computed, and still emitted as `law_prediction` so the published number stays
        #: reproducible and traceable; `law_prediction_quadmean` is the CORRECT one and is the
        #: only one that may be quoted.  See the correction block at the top of this file.
        r_bar = 0.5 * (rd + rp)                                     # SUPERSEDED -- arithmetic
        pred = float(np.sqrt(max(r_bar * r_bar - (s_dis / 2.0) ** 2, 0.0)))
        r_q2 = 0.5 * (rd * rd + rp * rp)                            # CORRECT -- quadratic
        pred_q = float(np.sqrt(max(r_q2 - (s_dis / 2.0) ** 2, 0.0)))

        rows.append({
            "pdb": p, "n": n, "fold": int(d["fold"]),
            "err_rms": rms, "err_mae": float(np.abs(e).mean()),
            "res_disto": res_d, "res_pool": res_p, "res_iid_null": res_iid,
            "res_perm": res_perm, "res_sign": res_sign,
            "rmsd_disto": rd, "rmsd_pool": rp, "disagreement": s_dis,
            "rmsd_coordavg": r_avg, "rmsd_restraintfusion": r_fuse,
            "law_prediction": pred,                 # SUPERSEDED, arithmetic mean of r
            "law_prediction_quadmean": pred_q,      # CORRECT, quadratic mean of r
            "pool_err_corr": float(np.corrcoef(e, d_pool - dtrue)[0, 1]),
        })
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": rows, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    g = lambda k: np.asarray([r[k] for r in rows], float)          # noqa: E731
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)

    out = {"n": len(rows), "incumbent": float(ref.mean()), "rows": rows}
    out["realizability"] = {
        "res_disto": float(g("res_disto").mean()),
        "res_pool": float(g("res_pool").mean()),
        "res_iid_null": float(g("res_iid_null").mean()),
        "res_perm": float(g("res_perm").mean()),
        "res_sign": float(g("res_sign").mean()),
        "err_mae": float(g("err_mae").mean()),
        "ratio_real_over_null": float(g("res_disto").mean() / g("res_iid_null").mean()),
        "disto_vs_null_paired": I.paired(g("res_disto"), g("res_iid_null"),
                                         folds=folds, names=pdbs),
        "perm_vs_null_paired": I.paired(g("res_perm"), g("res_iid_null"),
                                        folds=folds, names=pdbs),
    }
    out["fusion"] = {
        "r_disto": float(g("rmsd_disto").mean()), "r_pool": float(g("rmsd_pool").mean()),
        "disagreement": float(g("disagreement").mean()),
        "disagreement_median": float(np.median(g("disagreement"))),
        "measured_coordavg": float(g("rmsd_coordavg").mean()),
        # SUPERSEDED by the *_quadmean fields (Sprint 16 correction, 2026-09-06).  Retained so
        # the published +0.162 A stays reproducible from this module.
        "law_prediction": float(g("law_prediction").mean()),
        "prediction_error": float((g("rmsd_coordavg") - g("law_prediction")).mean()),
        "prediction_error_abs": float(np.abs(g("rmsd_coordavg")
                                             - g("law_prediction")).mean()),
        "law_prediction_quadmean": float(g("law_prediction_quadmean").mean()),
        "prediction_error_quadmean": float((g("rmsd_coordavg")
                                            - g("law_prediction_quadmean")).mean()),
        "prediction_error_abs_quadmean": float(np.abs(g("rmsd_coordavg")
                                                      - g("law_prediction_quadmean")).mean()),
        "restraint_fusion": float(g("rmsd_restraintfusion").mean()),
        "pool_err_corr": float(np.nanmean(g("pool_err_corr"))),
        #: MISLABELLING CORRECTED. `np.minimum` per target is a PER-TARGET ORACLE SELECTION --
        #: it requires knowing which channel won on each target. Calling it "the best single
        #: channel" and comparing fusion against it makes fusion look worse than it is, and on
        #: the full instrument the two comparisons have OPPOSITE SIGNS. Both are now reported,
        #: labelled for what they are.
        "coordavg_vs_better_channel_NATIVE_FREE": I.paired(
            g("rmsd_coordavg"), g("rmsd_disto"), folds=folds, names=pdbs),
        "restraintfusion_vs_better_channel_NATIVE_FREE": I.paired(
            g("rmsd_restraintfusion"), g("rmsd_disto"), folds=folds, names=pdbs),
        "coordavg_vs_ORACLE_per_target_pick": I.paired(
            g("rmsd_coordavg"), np.minimum(g("rmsd_disto"), g("rmsd_pool")),
            folds=folds, names=pdbs),
        "law_vs_measured": I.paired(g("law_prediction"), g("rmsd_coordavg"),
                                    folds=folds, names=pdbs),
    }
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_coherence", out, n_expected=len(tg))

    r = out["realizability"]; f = out["fusion"]
    print(f"\nn = {len(rows)} targets\n")
    print("REALIZABILITY -- mean |d_fit - d_target|; lower means the distance set is closer to")
    print("being achievable by an ideal-geometry chain. Every null arm is an ORACLE DIAGNOSTIC.")
    print(f"  the REAL predicted distances            {r['res_disto']:.3f} A")
    print(f"  the retrieval pool's distances          {r['res_pool']:.3f} A")
    print(f"  ORACLE i.i.d. Gaussian, matched RMS     {r['res_iid_null']:.3f} A   <- the null")
    print(f"  ORACLE errors permuted across pairs     {r['res_perm']:.3f} A")
    print(f"  ORACLE magnitudes kept, signs random    {r['res_sign']:.3f} A")
    print(f"  (mean |error| for reference             {r['err_mae']:.3f} A)")
    v = r["disto_vs_null_paired"]
    print(f"  ratio real/null {r['ratio_real_over_null']:.3f}   paired "
          f"{v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print(f"\nFUSION -- the two-member KROGH-VEDELSBY (1995) ambiguity decomposition,")
    print(f"  d_avg = sqrt((r1^2+r2^2)/2 - (s/2)^2).  Prior art; cite, do not claim.")
    print(f"  distogram fit {f['r_disto']:.3f} | pool fit {f['r_pool']:.3f} | "
          f"raw per-pair error correlation {f['pool_err_corr']:+.3f}")
    print(f"  structural disagreement s = {f['disagreement']:.3f} A "
          f"(median {f['disagreement_median']:.3f}) -- NATIVE-FREE")
    print(f"  law prediction, QUADRATIC mean (correct) {f['law_prediction_quadmean']:.3f}  vs  "
          f"measured coordinate average {f['measured_coordavg']:.3f}")
    print(f"  prediction error {f['prediction_error_quadmean']:+.3f} A "
          f"(mean absolute {f['prediction_error_abs_quadmean']:.3f})")
    print(f"  [SUPERSEDED, arithmetic mean of r: prediction {f['law_prediction']:.3f}, "
          f"error {f['prediction_error']:+.3f} -- kept for traceability only]")
    print(f"  restraint-level fusion {f['restraint_fusion']:.3f}")
    w = f["coordavg_vs_better_channel_NATIVE_FREE"]
    print(f"  coord-avg vs the better channel (chosen globally, NATIVE-FREE): "
          f"{w['mean_diff']:+.3f} [{w['ci95'][0]:+.3f},{w['ci95'][1]:+.3f}]")
    w = f["restraintfusion_vs_better_channel_NATIVE_FREE"]
    print(f"  restraint-fusion, same control:                                 "
          f"{w['mean_diff']:+.3f} [{w['ci95'][0]:+.3f},{w['ci95'][1]:+.3f}]")
    w = f["coordavg_vs_ORACLE_per_target_pick"]
    print(f"  coord-avg vs an ORACLE per-target pick of the better channel:   "
          f"{w['mean_diff']:+.3f} [{w['ci95'][0]:+.3f},{w['ci95'][1]:+.3f}]  "
          f"(NOT native-free; a ceiling, not a control)")
    return out


if __name__ == "__main__":
    run()
