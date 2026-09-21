#!/usr/bin/env python
"""s32/s32_D1_signtransfer.py -- S32 lane D rung D1-T, RE-EMITTED AS A COMMITTED SCRIPT.

Registered in s32/PREREG_S32_D.md, commit 34973b1b.

## Why this file exists at all -- a reproducibility defect of mine, fixed

D1-T's four numbers were first produced by an INLINE `python -c` one-liner.  The resulting
`s32/results/s32_D1_signtransfer.json` carried no `provenance` block, no module, no git commit,
no source hash and no pinned seed, and NO SCRIPT IN THE REPOSITORY PRODUCED IT -- while two
downstream scripts quoted its numbers as a ceiling and the ledger quoted them as a headline.
Charter section 61 requires verifiable experiments, artefacts and seeds; that artefact failed all
three.  Project memory records this exact shape twice already
(`findings-prose-is-not-evidence-of-code`); this was the third.  **The old artefact is superseded
by this file's output and must not be quoted.**

## The question

D1-N established `Var(rho) > 0` with `E[rho] = 0` on every scorer: in-band ordering content exists
and the per-target SIGN is missing.  D1-T asks whether that sign is a property of the TARGET.

## The three nulls, and why the first one is not enough

  nullPERM    permute `rr` inside half B.  This destroys ALL structure, so it sits at ~0 for every
              scorer INCLUDING PURE NOISE.  It answers "is there any relation?", NOT "is the
              relation's SIGN a property of the target?".  It is reported because it was
              registered, and it is NOT the null the claim needs.
  nullXTGT    apply ANOTHER TARGET's sign to this target's held-out half B.  THIS is the null the
              per-target claim needs: it holds the sign's marginal distribution fixed and destroys
              only the pairing between sign and target.
  globalSGN   a single leave-one-FOLD-out GLOBAL sign (the majority sign of the training folds),
              applied to every test target.  If a scorer's transfer is matched by this, the
              "latent" is a global bias, not a per-target one.

## Two controls that decide whether the instrument works at all

  NOISE       a per-candidate Gaussian scorer with no relation to anything.  A transfer here would
              mean the statistic is broken.  **This is the self-test on the defect that motivates
              the nullXTGT column** (contract rule 5: a verification must be able to fail).
  RG          plain radius of gyration -- one line of numpy, native-free, not a Hamiltonian and not
              physics.  Included because if a pure compactness scalar carries a LARGER per-target
              latent than any Hamiltonian, then "four unrelated Hamiltonians respond to one bit"
              overstates their independence.  D-20 already showed LEG_total's in-band sign agrees
              with Rg's on 68.3%.

## DEDUP

7.7% of band members are EXACT coordinate duplicates (122/126 targets, mean 69.2 distinct of 75).
A duplicate landing in both halves contributes the same (score, rr) to each and buys sign agreement
for free.  Both the raw and the DEDUPLICATED numbers are emitted; **the deduplicated ones are the
ones to quote.**

ORACLE / NOT DEPLOYABLE throughout: the sign is estimated from `rr` on half A.

    python s32/s32_D1_signtransfer.py
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

from s12 import instrument as I          # noqa: E402
from s24 import stats_lib as ST          # noqa: E402
from s16 import energy_lib as EL         # noqa: E402
from s32.s32_D1_inband import spear, rg  # noqa: E402

OUT = os.path.join(HERE, "results", "s32_D1_signtransfer_v2.json")
SEED = 32001007                 # PINNED
NSPLIT = 16                     # random half-splits per target
NXTGT = 200                     # cross-target sign draws per target
SCORERS = ("AMBER", "DIS", "LEG_total", "LEG_torsion", "RG", "NOISE")


def collect(dedup):
    """Per target: the score vectors, the ORACLE rr, and the fold."""
    rng = np.random.default_rng(SEED + (1 if dedup else 0))
    out = []
    for t in I.targets():
        pdb, seq = t["pdb"], t["seq"]
        u = I.load_univ(pdb)
        p = I.pool_idx(u, k=500)
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        W = np.asarray(u["W"][p], float)[sub]
        if dedup:
            _, keep = np.unique(np.round(W.reshape(len(W), -1), 6), axis=0, return_index=True)
            keep = np.sort(keep)
        else:
            keep = np.arange(len(sub))
        idx = sub[keep]
        rr = np.asarray(u["rr"][p], float)[idx]
        z = np.load(os.path.join(ROOT, "s24", "cache_amber", pdb + ".npz"))
        comp = EL.legacy_components_of_windows(
            seq, np.asarray(u["PHI"][p], float), np.asarray(u["PSI"][p], float))
        S = {"AMBER": np.asarray(z["e_amber"], float)[idx],
             "DIS": np.asarray(z["score_dist"], float)[idx],
             "LEG_total": EL.legacy_total_from(comp)[idx],
             "LEG_torsion": np.asarray(comp["torsion"], float)[idx],
             "RG": rg(W[keep]),
             "NOISE": rng.standard_normal(len(idx))}
        out.append({"pdb": pdb, "fold": int(t["fold"]), "rr": rr, "S": S,
                    "k_raw": int(len(sub)), "k_kept": int(len(idx))})
    return out


def run(dedup):
    data = collect(dedup)
    rng = np.random.default_rng(SEED + 7 + (1 if dedup else 0))
    pdbs = [d["pdb"] for d in data]
    fold = np.array([d["fold"] for d in data], int)
    folds = ST.pinned_folds(pdbs)
    assert np.array_equal(fold, folds)
    res = {}
    for sc in SCORERS:
        tr, pm, rB, sA, rho_full = [], [], [], [], []
        for d in data:
            v, rr, m = d["S"][sc], d["rr"], len(d["rr"])
            a, b, ra_, rb_ = [], [], [], []
            for _ in range(NSPLIT):
                q = rng.permutation(m)
                A, B = q[:m // 2], q[m // 2:]
                x, y = spear(v[A], rr[A]), spear(v[B], rr[B])
                if not (np.isfinite(x) and np.isfinite(y)):
                    continue
                a.append(np.sign(x) if x != 0 else 1.0)
                b.append(y)
                ra_.append(np.sign(x) * y)
                rb_.append(np.sign(x) * spear(v[B], rng.permutation(rr[B])))
            tr.append(np.mean(ra_)); pm.append(np.mean(rb_))
            rB.append(np.mean(b)); sA.append(np.mean(a))
            rho_full.append(spear(v, rr))
        tr = np.array(tr); pm = np.array(pm); rB = np.array(rB)
        sA = np.array(sA); rho_full = np.array(rho_full)
        sgn_marg = np.sign(np.where(sA == 0, 1.0, sA))
        # --- nullXTGT: another target's sign applied to this target's half B
        xt = []
        for _ in range(NXTGT):
            perm = rng.permutation(len(data))
            while np.any(perm == np.arange(len(data))):          # no self-assignment
                bad = perm == np.arange(len(data))
                perm[bad] = rng.permutation(len(data))[bad]
            xt.append(sgn_marg[perm] * rB)
        xt = np.array(xt)
        # --- globalSGN: one leave-one-FOLD-out majority sign for all test targets
        gs = np.zeros(len(data))
        for f in np.unique(fold):
            trn, te = fold != f, fold == f
            s = 1.0 if np.mean(sgn_marg[trn] > 0) >= 0.5 else -1.0
            gs[te] = s * rB[te]
        c_xt = ST.compare(tr, xt.mean(0), folds=folds, names=pdbs, label=sc + " vs nullXTGT")
        c_gs = ST.compare(tr, gs, folds=folds, names=pdbs, label=sc + " vs globalSGN")
        c_pm = ST.compare(tr, pm, folds=folds, names=pdbs, label=sc + " vs nullPERM")
        xm = abs(c_xt["effect"]) / c_xt["mde"] if c_xt["mde"] else 0.0
        res[sc] = {
            "transfer": float(tr.mean()), "transfer_median": float(np.median(tr)),
            "nullPERM": float(pm.mean()),
            "nullXTGT_mean": float(xt.mean()), "nullXTGT_draw_sd": float(xt.mean(1).std()),
            "nullXTGT_draws": NXTGT,
            "globalSGN": float(gs.mean()),
            "mean_rho_full": float(rho_full.mean()),
            "mean_abs_rho_halfB": float(np.abs(rB).mean()),
            "sign_marginal": float(max((sgn_marg > 0).mean(), (sgn_marg < 0).mean())),
            "vs_nullXTGT": {"effect": c_xt["effect"], "se": c_xt["se"], "mde": c_xt["mde"],
                            "x_mde": xm, "folds_same_sign": c_xt["folds_same_sign"],
                            "ci95_fold": c_xt["ci95_fold"]},
            "vs_globalSGN": {"effect": c_gs["effect"], "se": c_gs["se"],
                             "x_mde": (abs(c_gs["effect"]) / c_gs["mde"]) if c_gs["mde"] else 0.0,
                             "folds_same_sign": c_gs["folds_same_sign"]},
            "vs_nullPERM": {"effect": c_pm["effect"], "se": c_pm["se"],
                            "x_mde": (abs(c_pm["effect"]) / c_pm["mde"]) if c_pm["mde"] else 0.0,
                            "folds_same_sign": c_pm["folds_same_sign"]},
            "VERDICT_per_target": ("PER-TARGET" if (xm >= 1.0 and c_xt["folds_same_sign"] >= 4)
                                   else ("NOT MEASURED" if xm >= 0.7 else "NOT A RESULT"))}
    return res, data


def main():
    out = {"prereg_commit": "34973b1b", "seed": SEED, "n_splits": NSPLIT, "n_xtgt_draws": NXTGT,
           "provenance": ST.provenance(__file__),
           "supersedes": "s32/results/s32_D1_signtransfer.json (no provenance, no script; "
                         "do not quote)",
           "oracle": "ORACLE / NOT DEPLOYABLE -- the sign is estimated from rr on half A",
           "basis": ("in-band Spearman against true CA-RMSD inside the shipped top-75 band; a "
                     "DIAGNOSTIC, never differenced against a built-chain RMSD"),
           "crossing_price_provenance": ("memory `in-band-ordering-is-per-target` prices 2.0 A at "
                                         "in-band rho = 0.638"),
           "arms": {}}
    for dedup in (False, True):
        res, data = run(dedup)
        key = "DEDUP" if dedup else "RAW"
        out["arms"][key] = res
        if dedup:
            out["dedup_stats"] = {
                "mean_distinct_of_75": float(np.mean([d["k_kept"] for d in data])),
                "frac_duplicate": float(1 - np.mean([d["k_kept"] / d["k_raw"] for d in data])),
                "n_targets_with_duplicates": int(sum(d["k_kept"] < d["k_raw"] for d in data))}
    for key in ("RAW", "DEDUP"):
        print("\n=== %s ===" % key)
        print("%-13s %9s %9s %9s %9s %7s %6s  %s" % (
            "scorer", "transfer", "nullPERM", "nullXTGT", "globalSGN", "xMDE", "folds", "verdict"))
        for sc in SCORERS:
            r = out["arms"][key][sc]
            print("%-13s %+9.4f %+9.4f %+9.4f %+9.4f %7.2f %6s  %s" % (
                sc, r["transfer"], r["nullPERM"], r["nullXTGT_mean"], r["globalSGN"],
                r["vs_nullXTGT"]["x_mde"], r["vs_nullXTGT"]["folds_same_sign"],
                r["VERDICT_per_target"]))
    print("\ndedup:", json.dumps(out["dedup_stats"]))
    print("\nTHE CEILING (memory `in-band-ordering-is-per-target`: 2.0 A needs in-band rho 0.638)")
    for sc in SCORERS:
        r = out["arms"]["DEDUP"][sc]
        print("  %-13s perfect free sign -> %+0.4f   = %.0f%% of the 0.638 crossing price" % (
            sc, r["transfer"], 100 * r["transfer"] / 0.638))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    ST.save_atomic(OUT, out, module_file=__file__)
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
