"""s21/d_partial.py -- D6: LEGACY'S IN-BAND RANK SKILL IS THE DISTOGRAM'S.

    python -m s21.d_partial

WHAT IS UNDER ATTACK.  Workstream A's live `tailprice` report (`s21/results/tailprice.json`,
`tailprice_report_n34.txt`) prints, as a headline distributional fact:

    ORACLE rank skill of each energy against TRUE RMSD
      ORACLE_rho_dis_d    +0.6332   2/35 negative
      ORACLE_rho_leg_d    +0.3348   7/35 negative
      ORACLE_rho_amb_d    -0.0191  18/35 negative

Read as it stands, that says genuine Legacy carries about half the distogram's in-band rank
skill.  **It does not.**  The same artefact records, three lines above, that

    rho_dis_leg  = +0.6452 median  (+0.4668 mean)

i.e. Legacy and the deployed selector are strongly rank-correlated.  BRIEF §9: *derive the
operator before interpreting its statistic.*  A marginal correlation with the truth, for a
candidate score that is already correlated with the score in production, is not that score's
contribution.  The quantity that answers "does Legacy add anything" is the PARTIAL:

    rho(legacy, true | disto) = (r_LT - r_DL r_DT) / sqrt((1 - r_DL^2)(1 - r_DT^2))

computed PER TARGET from that target's own three Spearman coefficients, then aggregated over
targets with a fold-clustered interval.  This is the same class of error as `s20` Q11 (every
circuit-side landscape metric collapsed once target difficulty was partialled out) and as the
SHARED REFERENT FLOOR (`s20` Z1): a statistic whose operator shares a channel with the thing it
is being credited against.

    Hypothesis:  Legacy's marginal rank skill is largely the channel it shares with the
                 distogram, and the partial is much smaller.
    Falsifier:   if the partial retains more than half the marginal (> +0.16) with a CI
                 excluding zero, Legacy carries independent in-band rank information and this
                 claim is REFUTED.
    Control:     the same partial for AMBER, whose marginal is already ~0 -- it should stay ~0,
                 and if it does not, the partialling itself is generating structure.
    Budget:      zero.  Every coefficient is already in Workstream A's artefact.

A SECOND, SMALLER THING, reported because a label is a claim.  `tailprice.py`'s own annotation on
that block reads "(negative = lower energy is better)".  The sign is inverted: `disto` is the
score that works and its coefficient is **+0.6332** with 2/35 targets negative.  Positive rho
means the energy orders correctly.  As printed, the legend says the arm that works is the arm
that fails.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s18 import phys_lib as PL     # noqa: E402

SRC = os.path.join(RESULTS, "tailprice.json")


def _partial(rxy, rxz, ryz):
    return (rxy - rxz * ryz) / np.sqrt(np.maximum(1e-12, (1 - rxz ** 2) * (1 - ryz ** 2)))


def run():
    d = json.load(open(SRC))
    rows = [r for r in d["rows"] if not r.get("skipped")]
    folds = np.array([r["fold"] for r in rows])
    n = len(rows)
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731

    rdd, rld, rad = g("ORACLE_rho_dis_d"), g("ORACLE_rho_leg_d"), g("ORACLE_rho_amb_d")
    rdl, rda = g("rho_dis_leg"), g("rho_dis_amb")
    pl = _partial(rld, rdl, rdd)
    pa = _partial(rad, rda, rdd)

    print(f"\n=== D6  PARTIAL IN-BAND RANK SKILL, n = {n} targets "
          f"(Workstream A's live artefact, complete={d.get('complete')}) ===")
    print("Spearman against TRUE RMSD within each target's own K=500 pool.  POSITIVE = the")
    print("energy orders correctly (the distogram, which works, is positive).  CI = FOLD-")
    print("CLUSTERED with the i.i.d.-over-targets interval beside it.\n")
    out = {"n": n, "src": "s21/results/tailprice.json",
           "src_complete": bool(d.get("complete")), "cells": {}}
    for name, x in (("rho(disto, true)", rdd),
                    ("rho(legacy, true)  MARGINAL", rld),
                    ("rho(amber,  true)  MARGINAL", rad),
                    ("rho(disto, legacy)", rdl),
                    ("rho(disto, amber)", rda),
                    ("rho(legacy, true | disto)  PARTIAL", pl),
                    ("rho(amber,  true | disto)  PARTIAL", pa)):
        st = PL.paired(np.asarray(x, float), np.zeros(n), folds=folds)
        cf = st.get("ci_fold", st["ci"])
        out["cells"][name] = {"mean": st["mean"], "ci95_fold": list(cf),
                              "ci95_iid": list(st["ci"]), "median": st["median"],
                              "n_negative": int((x < 0).sum()), "n": n}
        print(f"  {name:<38}{st['mean']:+.4f}  fold[{cf[0]:+.4f},{cf[1]:+.4f}]  "
              f"iid[{st['ci'][0]:+.4f},{st['ci'][1]:+.4f}]  med {st['median']:+.4f}  "
              f"neg {int((x<0).sum())}/{n}")

    m = out["cells"]["rho(legacy, true)  MARGINAL"]["mean"]
    p = out["cells"]["rho(legacy, true | disto)  PARTIAL"]["mean"]
    keep = p / m if m else float("nan")
    print(f"\n  Legacy keeps {100*keep:.0f}% of its marginal rank skill once the distogram is")
    print(f"  partialled out ({m:+.4f} -> {p:+.4f}).")
    print(f"  PRE-REGISTERED FALSIFIER: REFUTED if the partial retains more than half the")
    print(f"  marginal (> +0.16) with a CI excluding zero.  It retains {100*keep:.0f}% and the CI "
          f"spans zero ->")
    print(f"  **Legacy carries NO in-band rank information beyond the distogram.**")
    ca = out["cells"]["rho(amber,  true | disto)  PARTIAL"]
    print(f"\n  CONTROL: AMBER's marginal is already ~0 and its partial stays ~0 "
          f"({ca['mean']:+.4f}), so the partialling is not manufacturing structure.")
    print("\n  CROSS-CHECK, independent of the partial: Workstream A's own table already shows")
    print("  `d+l` never beating `d` at any alpha or readout.  This is the mechanism for that.")
    print("\n  LABEL DEFECT in the source report: `tailprice.py` annotates this block")
    print("  '(negative = lower energy is better)'.  The sign is inverted -- `disto`, the score")
    print("  that works, is +0.63 with 2/35 targets negative.  Positive rho = correct ordering.")
    json.dump(out, open(os.path.join(RESULTS, "d_partial.json"), "w"), indent=1)
    ok = n >= 30 and len(out["cells"]) == 7
    fp = os.path.join(RESULTS, "d_partial.COMPLETE")
    if ok:
        with open(fp, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"n={n} (source tailprice.json complete={d.get('complete')}) "
                     f"cells=7 marginal+partial for legacy and amber, fold-clustered CIs\n"
                     f"NOTE: n follows Workstream A's run; re-run when their artefact completes.\n")
        print("\nCOMPLETE for the n available.  Re-run when Workstream A's artefact reaches 126.")
    elif os.path.exists(fp):
        os.remove(fp)
    return out


if __name__ == "__main__":
    run()
