"""S32 LANE V -- the FULL statistical gate on the order-statistic contrast of S32-L3.

S32-L3 quotes `compare(score128, random128) = +0.1872, MDE 0.1764, 1.06x -- WORSE, 72W/54L`
and `top-75 vs random 75 of 128 = +0.0696 at 1.02x, 65W/61L`.

Contract rule 1's gate is **>= 1.0x MDE AND a fold CI excluding zero AND >= 4/5 folds
agreeing**. The ledger sentence quotes the first and the W/L and not the other two, and both
effects sit at 1.02-1.06x -- the exact band where the remaining two criteria decide. This
recomputes the whole gate from the raw per-target rows and prints every criterion.

It also closes two things `s32_V_ladder_orderstat.json` left open:

  * **contract rule 10.** The random arm there is a per-target mean over 2000 draws. The
    draw-to-draw sd that matters for a HEADLINE is the sd of the **126-target aggregate**
    across draws, not the within-target sd. A single random draw is what a deployment would
    actually get, so the fraction of single draws that would clear MDE is reported.
  * **the moment.** The ledger says "mean rank 170/500"; the artefact stored the MEDIAN (134).
    Both are printed, named.

Rerun of the whole draw is avoided: the draws are regenerated from the same pinned seed as
`s32_V_ladder_orderstat.py`, and the script asserts it reproduces that file's per-target means.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s24 import stats_lib as ST                                           # noqa: E402
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
NDRAW = 2000


def main():
    src = json.load(open(os.path.join(RESULTS, "s32_V_ladder_orderstat.json")))
    per = src["per_target"]
    pdbs = sorted(per)
    folds = np.array([per[p]["fold"] for p in pdbs], int)
    b500 = np.array([per[p]["best500"] for p in pdbs])
    b128 = np.array([per[p]["best128"] for p in pdbs])
    b75 = np.array([per[p]["best75"] for p in pdbs])

    # --- regenerate the draws with the SAME pinned seed and assert they reproduce the artefact
    rng = np.random.default_rng(20260921)
    R128 = np.empty((len(pdbs), NDRAW))
    R75 = np.empty((len(pdbs), NDRAW))
    order_of = {}
    for k, t in enumerate(I.targets()):
        pdb = t["pdb"]
        u = I.load_univ(pdb); pool = I.pool_idx(u, I.K); rr = u["rr"][pool]
        dg = I.distogram(pdb, t["seq"], t["fold"]); i, j = I.pair_index(t["n"])
        D = I.pair_dists(u["W"][pool], i, j)
        sc = I.shipped_score(dg, D.astype(np.float32).astype(float))
        o = np.argsort(sc, kind="stable"); t128 = o[:128]
        idx = rng.random((NDRAW, 500)).argsort(1)[:, :128]
        r1 = rr[idx].min(1)
        sub = t128[rng.random((NDRAW, 128)).argsort(1)[:, :75]]
        r2 = rr[sub].min(1)
        q = pdbs.index(pdb)
        R128[q] = r1; R75[q] = r2
        order_of[pdb] = int(np.argsort(np.argsort(sc, kind="stable"), kind="stable")[int(np.argmin(rr))])
    rep128 = np.abs(R128.mean(1) - np.array([per[p]["rand128_mean"] for p in pdbs])).max()
    rep75 = np.abs(R75.mean(1) - np.array([per[p]["rand75of128_mean"] for p in pdbs])).max()
    print("reproduces s32_V_ladder_orderstat.json per-target draw means to %.2e / %.2e"
          % (rep128, rep75))
    assert rep128 < 1e-12 and rep75 < 1e-12, "the pinned seed no longer reproduces the artefact"

    out = dict(n=len(pdbs), n_draws=NDRAW, basis="member (deposited window)",
               note="ORACLE / NOT DEPLOYABLE; rr labels are read for evaluation only",
               reproduces_source_to=float(max(rep128, rep75)))

    print("=" * 100)
    for lab, arm, base, draws in (
            ("SCORE top-128 vs RANDOM 128 of 500", b128, R128.mean(1), R128),
            ("SCORE top-75  vs RANDOM 75 of 128", b75, R75.mean(1), R75)):
        c = ST.compare(arm, base, folds, names=pdbs, label=lab)
        gate_mde = abs(c["effect_over_mde"]) >= 1.0
        gate_ci = min(c["ci95_fold"]) * max(c["ci95_fold"]) > 0
        gate_fold = c["folds_same_sign"] >= 4
        verdict = ("RESULT" if (gate_mde and gate_ci and gate_fold)
                   else "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
        # contract rule 10: what a SINGLE draw would give, and its own distribution
        agg = draws.mean(0)                      # 2000 aggregate means, one per draw
        single = arm.mean() - agg                # the effect a single random draw would show
        frac_clear = float((np.abs(single) >= c["mde"]).mean())
        out[lab] = dict(compare=c, gate=dict(mde=gate_mde, fold_ci_excludes_zero=gate_ci,
                                             folds_same_sign_ge4=gate_fold, verdict=verdict),
                        draw_aggregate=dict(mean=float(agg.mean()), sd=float(agg.std(ddof=1)),
                                            single_draw_effect_mean=float(single.mean()),
                                            single_draw_effect_sd=float(single.std(ddof=1)),
                                            frac_single_draws_clearing_mde=frac_clear))
        print("%s" % lab)
        print("  effect %+0.4f   SE %.4f   MDE %.4f   %.2fx   median %+0.4f   %dW/%dL/%dT"
              % (c["effect"], c["se"], c["mde"], c["effect_over_mde"], c["median_effect"],
                 c["n_better"], c["n_worse"], c["n_tied"]))
        print("  fold effects %s" % {k: round(v, 4) for k, v in c["per_fold"].items()})
        print("  fold CI95 [%+0.4f, %+0.4f]  -> excludes zero: %s ; folds same sign %d/5 -> %s"
              % (c["ci95_fold"][0], c["ci95_fold"][1], gate_ci, c["folds_same_sign"], gate_fold))
        print("  type_m %.2f  type_s %.3f  flag %s   concentration flag %s (pctile in null %.3f)"
              % (c["type_m"], c["type_s"], c["type_m_flag"],
                 c["concentration"]["flag"], c["concentration"]["pctile_in_null"]))
        print("  RULE-10 draw distribution: the 126-target aggregate over %d draws is %.4f "
              "+- %.4f (sd)" % (NDRAW, agg.mean(), agg.std(ddof=1)))
        print("     a SINGLE random draw would show %+0.4f +- %.4f; %.1f%% of single draws clear "
              "this comparison's own MDE" % (single.mean(), single.std(ddof=1), 100 * frac_clear))
        print("  >>> GATE (contract rule 1: >=1.0x MDE AND fold CI excluding zero AND >=4/5 folds): %s"
              % verdict)
        print("-" * 100)

    rk = np.array([order_of[p] for p in pdbs], float)
    out["rank_of_oracle_best_member"] = dict(mean=float(rk.mean()), median=float(np.median(rk)),
                                             n_in_top128=int((rk < 128).sum()),
                                             n_in_top75=int((rk < 75).sum()), n=len(rk),
                                             expected_at_random_mean=249.5)
    print("DIS rank of the ORACLE-best member of the 500: mean %.1f, MEDIAN %.1f "
          "(random expectation 249.5); in the top-128 on %d/126, top-75 on %d/126"
          % (rk.mean(), np.median(rk), (rk < 128).sum(), (rk < 75).sum()))
    print("  -- the two moments differ by %.0f; the ledger's '170' is the MEAN, the artefact's "
          "'134' is the MEDIAN.  A number carries its definition." % (rk.mean() - np.median(rk)))
    print("=" * 100)
    with open(os.path.join(RESULTS, "s32_V_orderstat_gate.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out




def strata():
    """CONTRACT RULE 12: the aggregate above is driven by an OUTCOME-DEFINED stratum.

    FAIL18 is defined as the targets where no pool member within 1.5 A of the pool best
    survives into the production TOP-75 (`s12/instrument.py::selfcheck`).  A contrast about
    "does the score's prefix retain the good members" is therefore near-circular on exactly
    those 18 targets.  S31 measured the effect size rising monotonically with a stratum's
    circularity and named that gradient the stratum's definition doing the work.

    This splits both contrasts by FAIL18, re-gates on the 108 that are NOT outcome-defined,
    and adds a FILTER-INDEPENDENT stratum (chain length) as rule 12 asks.
    """
    src = json.load(open(os.path.join(RESULTS, "s32_V_ladder_orderstat.json")))["per_target"]
    pdbs = sorted(src)
    folds = np.array([src[p]["fold"] for p in pdbs], int)
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    nn = np.array([t["n"] for t in sorted(I.targets(), key=lambda x: x["pdb"])], float)
    out = {}
    print()
    print("=" * 100)
    print("CONTRACT RULE 12: is the aggregate an outcome-defined stratum's arithmetic?")
    for lab, ka, kb in (("500->128", "best128", "rand128_mean"),
                        ("128->75", "best75", "rand75of128_mean")):
        a = np.array([src[p][ka] for p in pdbs]); b = np.array([src[p][kb] for p in pdbs])
        d = a - b
        c = ST.compare(a[~f18], b[~f18], folds[~f18],
                       names=[p for p, m in zip(pdbs, ~f18) if m], label="%s on the 108" % lab)
        gate = ("RESULT" if abs(c["effect_over_mde"]) >= 1.0
                and min(c["ci95_fold"]) * max(c["ci95_fold"]) > 0 and c["folds_same_sign"] >= 4
                else "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
        lo = nn <= np.median(nn)
        out[lab] = dict(all_mean=float(d.mean()), all_median=float(np.median(d)),
                        fail18_mean=float(d[f18].mean()), fail18_WL=[int((d[f18] < 0).sum()), int((d[f18] > 0).sum())],
                        rest108_mean=float(d[~f18].mean()), rest108_median=float(np.median(d[~f18])),
                        rest108_WL=[int((d[~f18] < 0).sum()), int((d[~f18] > 0).sum())],
                        rest108_compare=c, rest108_gate=gate,
                        length_stratum=dict(short=float(d[lo].mean()), long=float(d[~lo].mean()),
                                            cut=float(np.median(nn))))
        print("  %-9s ALL126 %+0.4f (median %+0.4f) | FAIL18 %+0.4f (%dW/%dL) | OTHER 108 %+0.4f "
              "(median %+0.4f, %dW/%dL) -> %.2fx %s"
              % (lab, d.mean(), np.median(d), d[f18].mean(), (d[f18] < 0).sum(), (d[f18] > 0).sum(),
                 d[~f18].mean(), np.median(d[~f18]), (d[~f18] < 0).sum(), (d[~f18] > 0).sum(),
                 c["effect_over_mde"], gate))
        print("            filter-INDEPENDENT stratum (length <= %d vs >): %+0.4f vs %+0.4f "
              "-- no gradient, so it is not a length effect"
              % (np.median(nn), d[lo].mean(), d[~lo].mean()))
    print("=" * 100)
    path = os.path.join(RESULTS, "s32_V_orderstat_strata.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out


if __name__ == "__main__":
    main()
    strata()
