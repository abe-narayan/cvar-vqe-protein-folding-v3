"""SPRINT 22 / WORKSTREAM A -- READOUT-H vs TRAINING-H SEPARATION on the candidate encoding.

Requested directly by the coordinator as the second of three numbers. Mirrors s21 L16/B6
("changing the READOUT Hamiltonian is the one non-null lever in the selector, at zero extra
budget on the identical evaluated set") on the CANDIDATE-BASIS register instead of the continuous
torsion basin latent.

Design: train ONCE per (target, seed) on H_TRAIN (Legacy or DIST), T=0, alpha=0.15, canonical
label, 150 iterations (matching A1). Then read out the SAME trained distribution p_theta against
BOTH energies -- same-H (sanity, should reproduce A1's pattern) and cross-H (train on one, score
the tail/argmin with the other, no further training). Zero extra training cost: this is exactly
"at zero extra budget on the identical evaluated set" because p_theta is not touched.

No AMBER. No OpenMM.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import s22.qcand_lib as QC                      # noqa: E402
from s20 import qb2_lib as L20                  # noqa: E402

ALPHA = 0.15
ITERS = 150
SEEDS = [0, 1, 2, 3]


def readout_with(scores_other, enc_train, pool, p):
    """Score the trained distribution `p` (from `enc_train`'s training) against a DIFFERENT
    energy array `scores_other`, with no further training -- the cross-H readout."""
    enc_other = QC.Encoding(scores_other, label=enc_train.label)
    face = QC.exact_face(enc_other.E, p, ALPHA)
    _bits, cset = QC.argmin_readout(enc_other, p)
    argmin_rr, _ = QC.tail_average_rmsd(pool, cset)
    cands = QC.tail_candidates(enc_other, face["tail"])
    r_tailavg, _ = QC.tail_average_rmsd(pool, cands)
    return {"argmin_rr": float(argmin_rr), "r_tailavg": r_tailavg,
           "tail_size": face["tail_size"], "ess": face["ess"]}


def main():
    tg = L20.subset(n=16)
    rows = []
    t_start = time.time()
    for ti, t in enumerate(tg):
        pool = QC.build_pool(t["pdb"])
        enc_dist = QC.Encoding(pool["score_dist"])
        enc_leg = QC.Encoding(pool["score_leg"])
        exact_dist = QC.classical_exact_sort(pool, pool["score_dist"], ALPHA)
        exact_leg = QC.classical_exact_sort(pool, pool["score_leg"], ALPHA)

        seed_rows = []
        for seed in SEEDS:
            # train on DIST, readout DIST (same-H) and LEG (cross-H)
            out_d = QC.train(enc_dist, alpha_schedule=ALPHA, seed=seed, iters=ITERS,
                             entangler="cnot")
            same_d = readout_with(pool["score_dist"], enc_dist, pool, out_d["p_final"])
            cross_dl = readout_with(pool["score_leg"], enc_dist, pool, out_d["p_final"])

            # train on LEG, readout LEG (same-H) and DIST (cross-H)
            out_l = QC.train(enc_leg, alpha_schedule=ALPHA, seed=seed, iters=ITERS,
                             entangler="cnot")
            same_l = readout_with(pool["score_leg"], enc_leg, pool, out_l["p_final"])
            cross_ld = readout_with(pool["score_dist"], enc_leg, pool, out_l["p_final"])

            seed_rows.append({"seed": seed, "train_dist_readout_dist": same_d,
                              "train_dist_readout_leg": cross_dl,
                              "train_leg_readout_leg": same_l,
                              "train_leg_readout_dist": cross_ld})
        rows.append({"pdb": t["pdb"], "n": pool["n"], "n_qubits": enc_dist.n_qubits,
                    "R_pool": float(pool["rr"].min()),
                    "R_score_argmin_dist": exact_dist["r_argmin"],
                    "R_score_tailavg_dist": exact_dist["r_tailavg"],
                    "R_score_argmin_leg": exact_leg["r_argmin"],
                    "R_score_tailavg_leg": exact_leg["r_tailavg"],
                    "seeds": seed_rows})
        print(f"[readoutH {ti+1}/16] {t['pdb']} elapsed={time.time()-t_start:.0f}s", flush=True)

    def avg(path):
        vals = []
        for r in rows:
            vals.append(np.mean([_dig(s, path) for s in r["seeds"]]))
        return float(np.mean(vals)), vals

    def _dig(s, path):
        d = s
        for k in path:
            d = d[k]
        return d

    m_same_d_tail, tv_same_d = avg(["train_dist_readout_dist", "r_tailavg"])
    m_cross_dl_tail, tv_cross_dl = avg(["train_dist_readout_leg", "r_tailavg"])
    m_same_l_tail, tv_same_l = avg(["train_leg_readout_leg", "r_tailavg"])
    m_cross_ld_tail, tv_cross_ld = avg(["train_leg_readout_dist", "r_tailavg"])

    ci_train_leg_readout_effect = QC.paired_ci(tv_cross_ld, tv_same_l)     # readout DIST vs LEG, same training (LEG)
    ci_train_dist_readout_effect = QC.paired_ci(tv_cross_dl, tv_same_d)    # readout LEG vs DIST, same training (DIST)

    out = {
        "experiment": "A_readoutH_trainingH_separation", "alpha": ALPHA, "iters": ITERS,
        "seeds": SEEDS, "n_targets": len(rows), "rows": rows,
        "means_tailavg": {
            "train=DIST,readout=DIST (same-H)": m_same_d_tail,
            "train=DIST,readout=LEG (cross-H)": m_cross_dl_tail,
            "train=LEG,readout=LEG (same-H)": m_same_l_tail,
            "train=LEG,readout=DIST (cross-H)": m_cross_ld_tail,
        },
        "readout_effect_holding_training=LEG (DIST-readout minus LEG-readout)":
            ci_train_leg_readout_effect,
        "readout_effect_holding_training=DIST (LEG-readout minus DIST-readout)":
            ci_train_dist_readout_effect,
        "wall_time_s": time.time() - t_start,
    }
    req = ["pdb", "n_qubits", "R_pool", "seeds"]
    path = QC.write("a_readouth", out, required_keys=req)
    print("WROTE", path)
    print("means_tailavg:", out["means_tailavg"])
    print("readout effect (training=LEG fixed, switching readout DIST-vs-LEG):",
         ci_train_leg_readout_effect)
    print("readout effect (training=DIST fixed, switching readout LEG-vs-DIST):",
         ci_train_dist_readout_effect)


if __name__ == "__main__":
    main()
