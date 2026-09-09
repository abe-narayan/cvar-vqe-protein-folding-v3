"""s23/c1_probweight.py -- H_C1: PROBABILITY-WEIGHTED CVaR-VQE READOUT.

Pre-registered in `s23/PREREG_C.md` before this file read a single RMSD. Candidate-identity
register, `H|i> = E_i|i>` diagonal on the shipped distogram score (`s22/qcand_lib.build_pool` /
`Encoding`, reused as-is -- the exact Sprint 22 Workstream A setup, K<=500, n_qubits=9, exact
`StatevectorCircuit`). No AMBER.

Readout under test: X_out = sum_{i in tail} (mass_i / sum mass) * X_i, mass from
`core.quantum.cvar_from_probs`, against:
  (a) size-matched classical top-m average (uniform), m = the trained tail's realised size
  (b) unweighted average over the identical trained-tail SET (expected == (a), Gate 1)
  (c) the untrained circuit (`iters=0`, identical seed/architecture), weighted the same way,
      at ITS OWN realised tail size

alpha=0.15 fixed. T primary=0.5 (continues s22 D5's own sweep), secondary T in {0.1, 0.2}.
4 seeds/target, n=126 dev targets, point-cloud basis throughout.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s22 import qcand_lib as QC        # noqa: E402
from s12 import instrument as I        # noqa: E402
from core import quantum as Q          # noqa: E402
from core import pipeline as PL        # noqa: E402
from s23 import qc_lib as QL           # noqa: E402

ALPHA = 0.15
T_PRIMARY = 0.5
T_SECONDARY = (0.1, 0.2)
SEEDS = (0, 1, 2, 3)
LAYERS = 3
ITERS = 80
LR = 0.15
N_QUBITS = 9


def weighted_average(W: np.ndarray, weights: np.ndarray | None) -> np.ndarray:
    """Superpose on the (weighted) consensus-medoid frame, weighted mean -- exactly
    `core.pipeline.average_weighted`'s convention, reimplemented on `s12.instrument`'s
    primitives so both the weighted and uniform arms share one code path."""
    if len(W) == 1:
        return W[0]
    P = I.pairwise_rmsd(W)
    b = PL.consensus_medoid(P, weights)
    Sup = I.superpose_batch(W, W[b])
    if weights is None:
        w = np.full(len(W), 1.0 / len(W))
    else:
        w = np.asarray(weights, float)
        w = w / w.sum()
    return np.tensordot(w, Sup, axes=(0, 0))


def one_arm(enc: QC.Encoding, alpha: float, T: float, seed: int, iters: int, W_real: np.ndarray):
    """Train (or, at iters=0, just initialise) the circuit; return the tail's real-candidate
    indices, its normalised weights, the weighted-average RMSD-ready coordinate, and the
    matched classical (uniform, same SET) coordinate -- plus the gate check."""
    p, cv, H, circ = Q.run_cvar_vqe(enc.E, alpha, T, n=enc.n_qubits, layers=LAYERS,
                                    iters=iters, restarts=1, seed=seed, lr=LR)
    val, q, mass = Q.cvar_from_probs(enc.E, p, alpha)
    tail_bits = np.flatnonzero(mass > 0.0)
    real_mask = enc.cand_of_bit[tail_bits] >= 0
    n_pad_in_tail = int((~real_mask).sum())
    cand_idx = enc.cand_of_bit[tail_bits[real_mask]]
    w = mass[tail_bits[real_mask]]
    w = w / w.sum()
    m = int(len(cand_idx))
    # ---- GATE 1: tail SET == classical argsort(scores)[:m] SET, exact ----
    classical_top_m = set(np.argsort(enc.scores, kind="stable")[:m].tolist())
    gate1_pass = bool(set(cand_idx.tolist()) == classical_top_m)
    Wt = W_real[cand_idx]
    weighted_ca = weighted_average(Wt, w)
    classical_ca = weighted_average(Wt, None)     # (a)/(b): uniform over the SAME set
    return {
        "m": m, "n_pad_in_tail": n_pad_in_tail, "gate1_pass": gate1_pass,
        "entropy_bits": float(H) / np.log(2.0), "cvar": float(cv),
        "weighted_ca": weighted_ca, "classical_ca": classical_ca,
    }


def run_target(pdb: str, nat_ca: np.ndarray) -> dict:
    pool = QC.build_pool(pdb)
    enc = QC.Encoding(pool["score_dist"])
    W_real = pool["W"]
    rows = {}
    for T in (T_PRIMARY,) + T_SECONDARY:
        seed_rows = []
        for seed in SEEDS:
            tr = one_arm(enc, ALPHA, T, seed, ITERS, W_real)
            un = one_arm(enc, ALPHA, T, seed, 0, W_real)
            r_w_trained = I.ca_rmsd(tr["weighted_ca"], nat_ca)
            r_c_trained = I.ca_rmsd(tr["classical_ca"], nat_ca)
            r_w_untrained = I.ca_rmsd(un["weighted_ca"], nat_ca)
            r_c_untrained = I.ca_rmsd(un["classical_ca"], nat_ca)
            seed_rows.append({
                "seed": seed,
                "m_trained": tr["m"], "m_untrained": un["m"],
                "gate1_trained": tr["gate1_pass"], "gate1_untrained": un["gate1_pass"],
                "pad_trained": tr["n_pad_in_tail"], "pad_untrained": un["n_pad_in_tail"],
                "entropy_bits_trained": tr["entropy_bits"],
                "entropy_bits_untrained": un["entropy_bits"],
                "rmsd_w_trained": r_w_trained, "rmsd_c_trained": r_c_trained,
                "rmsd_w_untrained": r_w_untrained, "rmsd_c_untrained": r_c_untrained,
            })
        rows[str(T)] = seed_rows
    return rows


def run(limit: int | None = None):
    tg = sorted(I.targets(), key=lambda r: r["pdb"])
    if limit:
        tg = tg[:limit]
    all_rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        u = I.load_univ(t["pdb"])
        nat_ca = np.asarray(u["nat_ca"], float)
        rows = run_target(t["pdb"], nat_ca)
        all_rows.append({"pdb": t["pdb"], "n": t["n"], "fold": t["fold"], "by_T": rows})
        if (c + 1) % 10 == 0:
            el = time.time() - t0
            print(f"  {c + 1}/{len(tg)}  ({el:.1f}s, {el / (c + 1):.2f}s/target)", flush=True)
    complete = (len(all_rows) == 126)
    out = {"rows": all_rows, "complete": bool(complete), "n_expected": 126,
          "config": {"alpha": ALPHA, "T_primary": T_PRIMARY, "T_secondary": list(T_SECONDARY),
                     "seeds": list(SEEDS), "layers": LAYERS, "iters": ITERS, "lr": LR,
                     "n_qubits": N_QUBITS}}
    QL.save_json("c1_probweight_raw.json", out)
    print(f"done: {len(all_rows)}/126, complete={complete}, {time.time() - t0:.1f}s total")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    run(limit=a.limit)
