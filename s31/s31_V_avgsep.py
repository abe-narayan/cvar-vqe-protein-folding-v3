#!/usr/bin/env python
"""s31/s31_V_avgsep.py -- LANE V, V3: adversarial audit of AVG_SEP, the sprint's only live
deployable candidate.  Recomputes lane F's endpoint contrast from `s31_F_terminal_rows.jsonl`
independently, and checks it against every gate that applies.

Gates applied (none of them chosen after seeing the number):
  1. lane F's REGISTERED falsifier (PREREG_S31_F.md s10.3): mean_chain(AVG_SEP) - mean_chain(AVG)
     must reach <= -1.0 x MDE with a fold CI excluding zero; > -0.7 x MDE REFUTES.
  2. contract rule 31 / lane D standing rule 1: no built-chain claim below 0.0107 A (MEAN floor).
  3. lane P S31-L18: the PER-TARGET floor is 0.0134 mean / 0.0329 p90 / 0.2285 max.
  4. lane D standing rule 2: both sides projected in the SAME JOB from the SAME stored clouds.
  5. contract rule 1: state the basis.  Endpoint = BUILT CHAIN.
  6. the in-job AVG must reproduce the canonical production 3.2105 within the instrument spread.
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s24 import stats_lib as ST            # noqa: E402

ROWS = os.path.join(HERE, "results", "s31_F_terminal_rows.jsonl")
OUT = os.path.join(HERE, "results", "s31_V_avgsep.json")
SEED = 31098
CANON = 3.210533994943299        # s29/results/s29_O_chain_rows.jsonl :: item=prod
MEAN_FLOOR = 0.0107              # lane D standing rule 1
PT_FLOOR = dict(mean=0.0134, p90=0.0329, max=0.2285)   # lane P, S31-L18


def brief(o):
    return dict(label=o["label"], mean_a=o["mean_a"], mean_b=o["mean_b"], effect=o["effect"],
                median_effect=o["median_effect"], se=o["se"], mde=o["mde"],
                effect_over_mde=o["effect_over_mde"], ci95_fold=o["ci95_fold"],
                folds_same_sign=o["folds_same_sign"], W=o["n_better"], L=o["n_worse"], n=o["n"])


def gate(o):
    r = abs(o["effect_over_mde"]); ci = o["ci95_fold"]
    ex = bool(ci is not None and (ci[0] > 0 or ci[1] < 0))
    if r >= 1.0 and ex:
        return "MEASURED"
    if r >= 1.0:
        return "NOT MEASURED (>=1.0x MDE, fold CI spans zero)"
    if r >= 0.7:
        return "NOT MEASURED (0.7-1.0x MDE)"
    return "NOT A RESULT (<0.7x MDE)"


def main():
    R = [json.loads(l) for l in open(ROWS) if l.strip()]
    R.sort(key=lambda r: r["pdb"])
    assert len(R) == 126 and len({r["pdb"] for r in R}) == 126, \
        "INCOMPLETE %d rows -- contract rule 15, check the JOB not the file" % len(R)
    names = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    avg = np.array([r["chain_AVG"] for r in R], float)
    sep = np.array([r["chain_AVG_SEP"] for r in R], float)
    rec = np.array([r["prod_chain_record"] for r in R], float)

    out = dict(seed=SEED, n=len(R), BASIS="BUILT CHAIN (the endpoint, contract rule 1)",
               source=ROWS,
               same_job="both arms read from the SAME row of the SAME job, so lane D standing "
                        "rule 2 is satisfied by construction")

    # ---- gate 6: does the in-job AVG reproduce canonical production?
    out["reproduction"] = dict(
        in_job_AVG_mean=float(avg.mean()), canonical_prod=CANON,
        diff=float(avg.mean() - CANON),
        record_column_mean=float(rec.mean()),
        max_abs_per_target_vs_record=float(np.abs(avg - rec).max()),
        note="the in-job AVG is the correct baseline for the contrast; the canonical 3.2105 is "
             "quoted only to show the instrument is in its usual place")

    # ---- the registered contrast
    o = ST.compare(sep, avg, folds=folds, names=names,
                   label="V3.AVG_SEP - AVG (BUILT CHAIN, same job, native-free)",
                   seed_parts=("s31V", str(SEED)))
    eff, mde = o["effect"], o["mde"]
    out["V3_endpoint"] = dict(
        cmp=brief(o), gate=gate(o),
        mean_AVG_SEP=float(sep.mean()), mean_AVG=float(avg.mean()))

    # ---- the registered falsifier, applied verbatim
    out["registered_falsifier"] = dict(
        rule="PREREG_S31_F.md s10.3: CONFIRMED iff effect <= -1.0 x MDE AND fold CI excludes "
             "zero; REFUTED iff effect > -0.7 x MDE; otherwise NOT MEASURED",
        effect=float(eff), mde=float(mde),
        ci_excludes_zero=bool(o["ci95_fold"][0] > 0 or o["ci95_fold"][1] < 0),
        VERDICT=("CONFIRMED" if (eff <= -1.0 * mde and
                                 (o["ci95_fold"][0] > 0 or o["ci95_fold"][1] < 0))
                 else "REFUTED" if eff > -0.7 * mde else "NOT MEASURED"))

    # ---- the floors
    d = sep - avg
    out["floors"] = dict(
        mean_floor=MEAN_FLOOR, effect_abs=float(abs(eff)),
        clears_mean_floor=bool(abs(eff) > MEAN_FLOOR),
        ratio_to_mean_floor=float(abs(eff) / MEAN_FLOOR),
        per_target_floor=PT_FLOOR,
        frac_targets_inside_per_target_mean_floor=float(np.mean(np.abs(d) < PT_FLOOR["mean"])),
        frac_targets_inside_per_target_p90_floor=float(np.mean(np.abs(d) < PT_FLOOR["p90"])),
        note="the per-target floor does NOT invalidate a MEAN effect -- independent per-target "
             "noise averages down by sqrt(n). It bounds per-target STATEMENTS. Reported so the "
             "distinction is explicit rather than assumed.")

    # ---- concentration: is the mean carried by a few targets? (median-vs-mean warning)
    out["concentration"] = dict(
        mean=float(d.mean()), median=float(np.median(d)),
        sd=float(d.std(ddof=1)), W=int((d < 0).sum()), L=int((d > 0).sum()),
        tied=int((d == 0).sum()),
        top5_gain=[[names[i], float(d[i])] for i in np.argsort(d)[:5]],
        top5_loss=[[names[i], float(d[i])] for i in np.argsort(d)[-5:]],
        share_of_total_from_best_5=float(np.sort(d)[:5].sum() / d.sum()) if d.sum() else None,
        note="median-vs-mean is the free early warning; a raw drop-top threshold is NOT a valid "
             "test and is not run here")

    # ---- GEOMETRY: the registered secondary (PREREG s3) -- is it even instrumented?
    k0 = set(R[0].keys())
    out["geometry_instrumentation"] = {
        arm: dict(bond_cloud=("bond_cloud_%s" % arm) in k0,
                  bond_chain=("bond_chain_%s" % arm) in k0,
                  rg_cloud=("rg_cloud_%s" % arm) in k0,
                  rg_chain=("rg_chain_%s" % arm) in k0,
                  move=("move_%s" % arm) in k0,
                  P=("P_%s" % arm) in k0)
        for arm in ("AVG", "MED", "AVG_RG", "AVG_SEP")}
    out["geometry_note"] = (
        "PREREG_S31_F.md s3 registers, as a reported secondary, the projection penalty P and the "
        "Rg contraction for every arm. Any arm with False above has a registered diagnostic that "
        "was not recorded.")

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
