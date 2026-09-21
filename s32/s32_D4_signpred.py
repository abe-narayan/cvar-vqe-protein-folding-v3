#!/usr/bin/env python
"""s32/s32_D4_signpred.py -- S32 lane D rung D4-S: can the per-target in-band SIGN be supplied
NATIVE-FREE?  This is the deployable form of the lane's headline.

Registered in s32/PREREG_S32_D.md, commit 34973b1b, as the follow-on to D1-T.

D1-T established that the per-target in-band sign is a real latent that TRANSFERS across a split
half of a target's own band -- but estimating it there needs ORACLE labels on half the band.  The
deployable question is whether a model trained STRICTLY OUT OF FOLD on other targets' signs, using
only NATIVE-FREE per-target features, can supply it.

READOUT, chosen to be directly comparable to D1-T's numbers:
      oriented in-band rho  =  mean over targets of  sign(rho_hat_lfo) * rho_true
  rho_hat_lfo is a leave-one-FOLD-out ridge prediction of the target's in-band rho.
  D1-T's ORACLE-partial ceiling on the same readout, from the RE-EMITTED and DEDUPLICATED
  s32_D1_signtransfer_v2.json (the original s32_D1_signtransfer.json is SUPERSEDED -- it had no
  script, no provenance and no seed): AMBER +0.1163, LEG_total +0.2180, DIS +0.1901,
  LEG_torsion +0.1492.  Plain Rg, not a Hamiltonian, reaches +0.3228 and beats all of them.
  The do-nothing floor is the unoriented in-band rho: AMBER +0.0000, LEG_total +0.0376.

CONTROLS
  (i)  LABEL-SHUFFLED: the same ridge with the training signs permuted WITHIN each training fold
       -- matched in every way except the association.  Own distribution, 32 draws, mean and sd
       reported, never a best draw (contract rule 10).
  (ii) CONSTANT-SIGN: always predict +1.  This is the plausible zero-information control, not a
       uniform coin flip (memory `zero-information-control-must-be-plausible`); it returns exactly
       the unoriented in-band rho.

LEAKAGE.  No feature reads `rr` or `nat_ca`.  The training TARGETS are other folds' in-band rho,
which are native-derived -- that is the same discipline the shipped leave-fold-out distogram uses
and is deployable, but it is stated here rather than assumed.

    python s32/s32_D4_signpred.py
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
from s32.s32_D1_inband import spear, rg, ca_pseudo_torsion   # noqa: E402

OUT = os.path.join(HERE, "results", "s32_D4_signpred.json")
SCORERS = ("AMBER", "DIS", "LEG_total", "LEG_torsion")
HYDRO = set("AVLIMFWCY")


def features_and_rho():
    rows = []
    for t in I.targets():
        pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
        u = I.load_univ(pdb)
        p = I.pool_idx(u, k=500)
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        W = np.asarray(u["W"][p], float)
        PHI = np.asarray(u["PHI"][p], float)
        PSI = np.asarray(u["PSI"][p], float)
        rr = np.asarray(u["rr"][p], float)[sub]                 # ORACLE, label only
        z = np.load(os.path.join(ROOT, "s24", "cache_amber", pdb + ".npz"))
        comp = EL.legacy_components_of_windows(seq, PHI, PSI)
        S = {"AMBER": np.asarray(z["e_amber"], float)[sub],
             "DIS": np.asarray(z["score_dist"], float)[sub],
             "LEG_total": EL.legacy_total_from(comp)[sub],
             "LEG_torsion": np.asarray(comp["torsion"], float)[sub]}
        Wb, G = W[sub], rg(W)[sub]
        pw = I.pairwise_rmsd(Wb)[np.triu_indices(len(sub), 1)]
        tau = ca_pseudo_torsion(Wb)
        dg = I.distogram(pdb)
        exp = np.asarray(dg["expected"], float)
        sd = np.asarray(dg["sd"], float)
        ii, jj = I.pair_index(n, min_sep=2)
        Dband = np.linalg.norm(Wb[:, ii] - Wb[:, jj], axis=-1)
        # predicted vs realised Rg-like scale, the native-free disagreement channel
        rg_pred = float(np.sqrt((np.sum(exp ** 2) * 2 + (n - 1) * 2 * 3.8 ** 2) / (2.0 * n * n)))
        f = {"n": float(n),
             "band_spread": float(np.median(pw)),
             "band_spread_sd": float(pw.std()),
             "rg_mean": float(G.mean()), "rg_sd": float(G.std()),
             "rg_disagree": rg_pred - float(G.mean()),
             "chi_mean": float(np.sin(tau).mean()), "chi_sd": float(np.sin(tau).std(0).mean()),
             "disto_sd": float(sd.mean()), "disto_sd_sd": float(sd.std()),
             "band_d_sd": float(Dband.std(0).mean()),
             "exp_minus_band": float((exp - Dband.mean(0)).mean()),
             "frac_hydro": float(sum(c in HYDRO for c in seq) / n),
             "frac_gly": float(seq.count("G") / n), "frac_pro": float(seq.count("P") / n),
             "org_frac": float(np.asarray(u["org"][p], bool)[sub].mean()),
             "sim_mean": float(np.asarray(u["sim"][p], float)[sub].mean())}
        for k, v in S.items():
            f["scoresd_" + k] = float(np.std(v) / (abs(np.mean(v)) + 1e-9))
        r = {"pdb": pdb, "fold": int(t["fold"]), "feat": f}
        for k, v in S.items():
            r["rho_" + k] = spear(v, rr)                         # ORACLE label
        rows.append(r)
    return rows


def lfo_ridge(X, y, fold, alphas=(1e-1, 1e0, 1e1, 1e2, 1e3)):
    """Leave-one-fold-out ridge with the alpha chosen by an INNER leave-one-fold-out loop on the
    training folds only, so no test fold touches the hyperparameter."""
    yh = np.zeros(len(y))
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        best, ba = np.inf, alphas[0]
        for a in alphas:
            err = []
            for g in np.unique(fold[tr]):
                i1 = tr & (fold != g)
                i2 = tr & (fold == g)
                if i2.sum() == 0 or i1.sum() < 5:
                    continue
                m, s = X[i1].mean(0), np.where(X[i1].std(0) > 1e-12, X[i1].std(0), 1.0)
                A = (X[i1] - m) / s
                B = (X[i2] - m) / s
                w = np.linalg.solve(A.T @ A + a * np.eye(A.shape[1]), A.T @ (y[i1] - y[i1].mean()))
                err.append(float(np.mean((B @ w + y[i1].mean() - y[i2]) ** 2)))
            if err and np.mean(err) < best:
                best, ba = np.mean(err), a
        m, s = X[tr].mean(0), np.where(X[tr].std(0) > 1e-12, X[tr].std(0), 1.0)
        A = (X[tr] - m) / s
        B = (X[te] - m) / s
        w = np.linalg.solve(A.T @ A + ba * np.eye(A.shape[1]), A.T @ (y[tr] - y[tr].mean()))
        yh[te] = B @ w + y[tr].mean()
    return yh


def main():
    rows = features_and_rho()
    pdbs = [r["pdb"] for r in rows]
    fold = np.array([r["fold"] for r in rows], int)
    folds = ST.pinned_folds(pdbs)
    assert np.array_equal(fold, folds)
    names = sorted(rows[0]["feat"])
    X = np.array([[r["feat"][k] for k in names] for r in rows], float)
    rng = np.random.default_rng(3200444)
    out = {"prereg_commit": "34973b1b", "n": len(rows), "features": names,
           "note": ("Oriented in-band rho = mean over targets of sign(rho_hat_lfo) * rho_true.  "
                    "D1-T's ORACLE-partial ceiling on the same readout is in `oracle_ceiling`; "
                    "the plausible zero-information control CONSTANT_PLUS returns the unoriented "
                    "in-band rho.  No feature reads rr or nat_ca; the training TARGETS are other "
                    "folds' in-band rho and are native-derived, which is the shipped "
                    "leave-fold-out discipline."),
           "oracle_ceiling_from_D1T_v2_DEDUP": {"AMBER": 0.1163, "DIS": 0.1901,
                                                "LEG_total": 0.2180, "LEG_torsion": 0.1492,
                                                "RG": 0.3228, "NOISE_selftest": 0.0056},
           "oracle_ceiling_source": "s32/results/s32_D1_signtransfer_v2.json "
                                    "(s32_D1_signtransfer.json is SUPERSEDED: no script/provenance/seed)",
           "arms": {}}
    print("%-13s %10s %10s %10s %9s %6s %6s" % (
        "scorer", "oriented", "const(+1)", "shuffled", "excess", "xMDE", "folds"))
    for sc in SCORERS:
        y = np.array([r["rho_" + sc] for r in rows], float)
        yh = lfo_ridge(X, y, fold)
        oriented = np.sign(yh) * y
        const = y                                   # sign == +1 always
        sh = []
        for _ in range(32):
            ys = y.copy()
            for f in np.unique(fold):               # permute WITHIN each fold
                m = fold == f
                ys[m] = rng.permutation(ys[m])
            sh.append(np.sign(lfo_ridge(X, ys, fold)) * y)
        sh = np.array(sh)
        shm = sh.mean(0)
        c = ST.compare(oriented, shm, folds=folds, names=pdbs, label=sc + " oriented vs shuffled")
        c2 = ST.compare(oriented, const, folds=folds, names=pdbs, label=sc + " oriented vs const")
        out["arms"][sc] = {
            "oriented_mean": float(oriented.mean()), "oriented_se": float(c["se"]),
            "CONSTANT_PLUS_zero_info": float(const.mean()),
            "shuffled_draw_mean": float(sh.mean()), "shuffled_draw_sd": float(sh.mean(1).std()),
            "n_shuffle_draws": 32,
            "vs_shuffled": {"effect": c["effect"], "se": c["se"], "mde": c["mde"],
                            "x_mde": (abs(c["effect"]) / c["mde"]) if c["mde"] > 0 else 0.0,
                            "folds_same_sign": c["folds_same_sign"], "ci95_fold": c["ci95_fold"]},
            "vs_constant": {"effect": c2["effect"], "se": c2["se"], "mde": c2["mde"],
                            "x_mde": (abs(c2["effect"]) / c2["mde"]) if c2["mde"] > 0 else 0.0,
                            "folds_same_sign": c2["folds_same_sign"]},
            "sign_accuracy_vs_truth": float((np.sign(yh) == np.sign(y)).mean()),
            "r2_lfo": float(1 - np.sum((yh - y) ** 2) / np.sum((y - y.mean()) ** 2))}
        print("%-13s %+10.4f %+10.4f %+10.4f %+9.4f %6.2f %6s" % (
            sc, oriented.mean(), const.mean(), sh.mean(), c["effect"],
            (abs(c["effect"]) / c["mde"]) if c["mde"] > 0 else 0.0, c["folds_same_sign"]))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    ST.save_atomic(OUT, out, n_expected=len(rows), module_file=__file__)
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
