#!/usr/bin/env python
"""s32/s32_D5_signchain.py -- S32 lane D rung D5: WHAT IS ONE ORACLE BIT PER TARGET WORTH AT THE
BUILT-CHAIN ENDPOINT?

Registered in s32/PREREG_S32_D.md, commit 34973b1b, under the lane's DEPLOYMENT CONDITION: an
in-band rho result is a diagnostic and counts only when carried to the built chain, projected in
the SAME job as its paired production rows.

Charter section 41 asks whether the five-bit oracle result is a hard limit or a pricing clue.  D1-N,
D1-T and D4-S have reduced the missing in-band quantity to exactly ONE BIT PER TARGET -- the sign
of the scorer's in-band ordering -- shown to be real (1.88-3.10x a matched permutation
null), to TRANSFER across a split half of a target's own band against a CROSS-TARGET null
(s32_D1_signtransfer_v2.json, deduplicated: AMBER +0.1163, LEG_total +0.2180, RG +0.3228, and a
NOISE self-test at +0.0056 that correctly returns NOT A RESULT), and to be UNPREDICTABLE above its
own MARGINAL from native-free per-target features.  This file prices that bit in Angstroms on the endpoint.

ARMS, all projected in THIS process so the pairing is within-job (contract rule 3):
  PROD           uniform coordinate average of the shipped top-75  -> STAGE 3b chain.
  CONST_m25      the best 25 of the 75 by the RAW scorer (sign fixed to +1) -> average -> chain.
                 This is the matched zero-information control: same m, same scorer, no bit.
  ORACLE_m25     the best 25 by sign(rho_true) * scorer.            ORACLE / NOT DEPLOYABLE.
The ONLY difference between CONST and ORACLE is one bit per target, so their paired difference IS
the price of the bit, matched in the operator's own space (contract rule 6).

Geometry secondaries ship on row one (contract rule 15): virtual-bond mean and sd, and the
displacement of each arm from the production chain it claims to improve on.

    python s32/s32_D5_signchain.py --scorer LEG_total
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I          # noqa: E402
from s24 import stats_lib as ST          # noqa: E402
from s16 import energy_lib as EL         # noqa: E402
from s32.s32_D1_inband import spear      # noqa: E402

RES = os.path.join(HERE, "results")
M = 25


def vb(X):
    d = np.linalg.norm(X[1:] - X[:-1], axis=-1)
    return float(d.mean()), float(d.std())


def one(t, scorer):
    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    u = I.load_univ(pdb)
    p = I.pool_idx(u, k=500)
    sub = np.asarray(I.shipped_record(pdb)["sub"], int)
    W = np.asarray(u["W"][p], float)[sub]
    rr = np.asarray(u["rr"][p], float)[sub]               # ORACLE label
    nat = np.asarray(u["nat_ca"], float)
    if scorer == "AMBER":
        z = np.load(os.path.join(ROOT, "s24", "cache_amber", pdb + ".npz"))
        v = np.asarray(z["e_amber"], float)[sub]
    else:
        comp = EL.legacy_components_of_windows(
            seq, np.asarray(u["PHI"][p], float), np.asarray(u["PSI"][p], float))
        v = EL.legacy_total_from(comp)[sub]
    rho = spear(v, rr)
    s = np.sign(rho) if rho != 0 else 1.0                 # the ORACLE BIT

    row = {"pdb": pdb, "n": int(t["n"]), "fold": fold, "rho_inband_ORACLE": float(rho),
           "bit_ORACLE": float(s)}
    chains = {}
    for lab, sg in (("PROD", None), ("CONST", 1.0), ("ORACLE", s)):
        if sg is None:
            C, _ = I.coordinate_average(W)
        else:
            idx = np.argsort(sg * v, kind="mergesort")[:M]
            C, _ = I.coordinate_average(W[idx])
        pr = I.project(C, seq, fold)
        ca = np.asarray(pr["ca"], float)
        chains[lab] = ca
        row["%s_cloud" % lab] = float(I.ca_rmsd(C, nat))          # CLOUD basis, diagnostic
        row["%s_chain" % lab] = float(I.ca_rmsd(ca, nat))         # BUILT CHAIN, the endpoint
        row["%s_vb_mean" % lab], row["%s_vb_sd" % lab] = vb(ca)
    for lab in ("CONST", "ORACLE"):
        row["%s_move_from_prod" % lab] = float(I.ca_rmsd(chains[lab], chains["PROD"]))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scorer", default="LEG_total", choices=["LEG_total", "AMBER"])
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    jl = os.path.join(RES, "s32_D5_signchain_%s.jsonl" % a.scorer)
    done = set()
    if os.path.exists(jl):
        with open(jl, encoding="utf-8") as fh:
            for ln in fh:
                try:
                    done.add(json.loads(ln)["pdb"])
                except Exception:                                  # noqa: BLE001
                    pass
    t0 = time.time()
    for i, t in enumerate(I.targets()):
        if t["pdb"] in done:
            continue
        r = one(t, a.scorer)
        with open(jl, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
        if (i + 1) % 10 == 0:
            print("[D5 %d/126] %s  %.0fs" % (i + 1, t["pdb"], time.time() - t0), flush=True)

    rows = [json.loads(l) for l in open(jl, encoding="utf-8") if l.strip()]
    rows.sort(key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    g = lambda k: np.array([r[k] for r in rows], float)            # noqa: E731
    out = {"prereg_commit": "34973b1b", "scorer": a.scorer, "m": M, "n": len(rows),
           "note": ("BUILT CHAIN is the endpoint; the cloud column is a diagnostic and is never "
                    "differenced against the chain.  All three arms projected in THIS process, so "
                    "the pairing is within-job (contract rule 3).  ORACLE arms are ORACLE / NOT "
                    "DEPLOYABLE: the bit is sign(rho_true) and reads native RMSD."),
           "PROD_chain_mean": float(g("PROD_chain").mean()),
           "PROD_cloud_mean": float(g("PROD_cloud").mean()), "arms": {}}
    for lab in ("CONST", "ORACLE"):
        for basis in ("chain", "cloud"):
            c = ST.compare(g("%s_%s" % (lab, basis)), g("PROD_%s" % basis), folds=folds,
                           names=pdbs, label="%s vs PROD (%s)" % (lab, basis))
            out["arms"]["%s_vs_PROD_%s" % (lab, basis)] = {
                "mean": float(g("%s_%s" % (lab, basis)).mean()), "effect": c["effect"],
                "median_effect": c["median_effect"], "se": c["se"], "mde": c["mde"],
                "x_mde": abs(c["effect"]) / c["mde"] if c["mde"] else 0.0,
                "ci95_fold": c["ci95_fold"], "folds_same_sign": c["folds_same_sign"],
                "n_better": c["n_better"], "n_worse": c["n_worse"], "verdict": c["verdict"]}
    # THE PRICE OF THE BIT: ORACLE against its own matched zero-information control.
    for basis in ("chain", "cloud"):
        c = ST.compare(g("ORACLE_%s" % basis), g("CONST_%s" % basis), folds=folds, names=pdbs,
                       label="PRICE OF ONE ORACLE BIT (%s)" % basis)
        out["arms"]["PRICE_OF_THE_BIT_%s" % basis] = {
            "effect": c["effect"], "median_effect": c["median_effect"], "se": c["se"],
            "mde": c["mde"], "x_mde": abs(c["effect"]) / c["mde"] if c["mde"] else 0.0,
            "ci95_fold": c["ci95_fold"], "folds_same_sign": c["folds_same_sign"],
            "n_better": c["n_better"], "n_worse": c["n_worse"], "verdict": c["verdict"],
            "ORACLE": "ORACLE / NOT DEPLOYABLE"}
    out["geometry"] = {lab: {"vb_mean": float(g("%s_vb_mean" % lab).mean()),
                             "vb_sd": float(g("%s_vb_sd" % lab).mean())}
                       for lab in ("PROD", "CONST", "ORACLE")}
    out["geometry"]["move_from_prod"] = {
        lab: {"mean": float(g("%s_move_from_prod" % lab).mean()),
              "max": float(g("%s_move_from_prod" % lab).max())} for lab in ("CONST", "ORACLE")}
    out["bit_marginal_frac_positive"] = float((g("bit_ORACLE") > 0).mean())
    ST.save_atomic(os.path.join(RES, "s32_D5_signchain_%s.json" % a.scorer), out,
                   n_expected=126, module_file=__file__)
    print(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
