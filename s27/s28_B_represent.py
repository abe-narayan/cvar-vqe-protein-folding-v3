#!/usr/bin/env python
"""s27/s28_B_represent.py -- the REPRESENTABILITY FIT lane D asked for in S28-L26 (the missing
piece of the three-way split's middle leg, S28-L2 caveat (b)).

Question: at J = 3 the exact ground state of H = diag(E) - J A (PR 305, hopping 0.911, sign
coherent) has F = -7.31 against the circuit's reached -5.51 / -5.27 (`s28_B_split.json`). Is
the 1.8 gap an OPTIMISATION shortfall (the 27-parameter RY+CNOT circuit HOLDS a state close to
the ground state and Adam did not find it) or an EXPRESSIVITY limit (no theta comes close)?

Measurement (native-free, no RMSD): for each of S27's 12 trainability targets
(`P.targets()[::11][:12]`), the 9-qubit register of the endpoint run (E = zrank(DIS) plus 12
padding states, A the REAL Gaussian graph, unit spectral norm, zero-padded), the exact ground
state g at J = 3, and

    theta* = argmax_theta  <psi(theta)| g g^T |psi(theta)>  =  (g . psi(theta))^2

by Adam (lr 0.15, the loop of `run_hop_vqe`) from 16 starts theta ~ N(0, 0.6^2) (draw 1 is the
seed-0 VQE's own initial theta), at the registered budget (80 iterations, "the same loop") and
at 5x that budget (400 iterations, so the answer is not the optimiser's budget; the
one-target probe read 0.811 at 80 and 0.814 at 800). The gradient
is the exact shift rule for the projector observable g g^T (`grad_overlap2`, tested equal to
`grad_hop_paramshift` with A = g g^T and to finite differences). Reported per target: the best squared overlap over the
16 starts at each budget, the circuit's objective F at that fitted state (`F_of_state`, the
split's arithmetic), its hopping value and sign coherence; beside it the trained VQE state's
own squared overlap with g (seeds 0 and 1, re-run bit-identically) and the untrained best-of-16
overlap. Also the same fit for the SIGN-ALIGNED VQE state (|psi_vqe|, the state the same p
would give with Perron-aligned signs), so "expressivity" is read on a second target vector.

Reading rule (S28-L26): overlap near 1 makes the 1.8 gap an optimisation shortfall; overlap
well below 1 makes it expressivity. Both halves are quoted as numbers; nothing here is a
plateau or its absence (contract rule 9).

Usage:  python s27/s28_B_represent.py [--iters-long 400] [--starts 16]
Writes `s27/results/s28_B_represent_rows.jsonl` (per target, resumable) and
`s27/results/s28_B_represent.json`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q              # noqa: E402
from s22 import qcand_lib as QC            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402
from s27.s28_B_split import F_of_state     # noqa: E402

OUT = os.path.join(B.RESULTS, "s28_B_represent.json")
ROWS = os.path.join(B.RESULTS, "s28_B_represent_rows.jsonl")
J_FIT = 3.0
N_STARTS = 16
ITERS_SHORT = B.ITERS          # 80: the registered loop
ITERS_LONG = 400               # 5x: so the answer is not the budget (the probe saturated by 80)


def overlap2(psi: np.ndarray, g: np.ndarray) -> float:
    """(g . psi)^2 for real unit vectors."""
    return float(np.dot(g, psi) ** 2)


def grad_overlap2(circ: Q.StatevectorCircuit, theta: np.ndarray, g: np.ndarray) -> np.ndarray:
    """d (g . psi(theta))^2 / d theta by the exact two-term shift rule (pi/2, factor 1/2).

    Identical arithmetic to `grad_hop_paramshift(circ, theta, g g^T)` -- for the projector
    observable <psi|g g^T|psi> = (g . psi)^2 the shifted expectations are (g . psi_+/-)^2, so
    only the 2P inner products S @ g are needed, not S @ (g g^T). Tested equal to the
    projector path to 1e-12 in `tests/test_s28_B.py`.
    """
    S = circ.states_batch(circ._shift_grid(np.asarray(theta, float), np.pi / 2))   # (2P, dim)
    vals = (S @ g) ** 2
    return (vals[0::2] - vals[1::2]) / 2.0


def fit_overlap(circ: Q.StatevectorCircuit, g: np.ndarray, rng: np.random.Generator,
                n_starts: int, iters: int, lr: float = B.LR) -> Dict:
    """Adam on -(g . psi(theta))^2 from `n_starts` draws theta ~ N(0, 0.6^2).

    The observable is the projector g g^T; the exact gradient is `grad_overlap2` (the shift
    rule specialised to a projector, equal to `grad_hop_paramshift` with A = g g^T). Returns
    the best final overlap, its theta, the per-start finals and the untrained (iteration-0)
    best.
    """
    Pn = circ.n_params()
    finals, inits = [], []
    best_ov, best_th = -1.0, None
    for _ in range(n_starts):
        th = rng.normal(0.0, 0.6, Pn)
        inits.append(overlap2(circ.state(th), g))
        m = np.zeros_like(th)
        v = np.zeros_like(th)
        for t in range(1, iters + 1):
            gr = -grad_overlap2(circ, th, g)                    # minimise -overlap
            m = 0.9 * m + 0.1 * gr
            v = 0.999 * v + 0.001 * gr * gr
            th = th - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
        ov = overlap2(circ.state(th), g)
        finals.append(ov)
        if ov > best_ov:
            best_ov, best_th = ov, th.copy()
    return dict(best=float(best_ov), theta=best_th, finals=[float(x) for x in finals],
                untrained_best=float(max(inits)), untrained_mean=float(np.mean(inits)))


def endpoint_F_of(pdb: str) -> Dict[int, float]:
    """The committed endpoint rows' F at REAL, J = 3, per seed (the bit-identity anchor)."""
    out = {}
    if not os.path.exists(B.ROWS):
        return out
    with open(B.ROWS, encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            if (r.get("pdb") == pdb and r.get("source") == "vqe" and r.get("graph") == "REAL"
                    and float(r.get("J", -1)) == J_FIT):
                out[int(r["seed"])] = float(r["F"])
    return out


def target_row(pdb: str, n_starts: int, iters_long: int) -> Dict:
    endpoint_F = endpoint_F_of(pdb)
    cand, ch, _ = RP.channels_for(pdb)
    E = RP.zr(ch["DIS"])
    enc = QC.Encoding(E)
    D, G = B.build_graphs(cand.W, pdb)
    A = B.pad_graph(G["REAL"]["A"], enc.dim)
    circ = Q.StatevectorCircuit(enc.n_qubits, B.LAYERS)
    gs = B.ground_state(enc.E, A, J_FIT)
    g = gs["psi"]
    f_gs = F_of_state(g, enc.E, A, J_FIT)
    rec = dict(pdb=pdb, J=J_FIT, n_qubits=enc.n_qubits, dim=enc.dim, n_starts=n_starts,
               iters_short=ITERS_SHORT, iters_long=iters_long,
               gs_pr=float(1.0 / np.sum(gs["p"] ** 2)), gs_gap=gs["gap"], F_gs=f_gs["F"],
               hop_gs=f_gs["hop"], coh_gs=f_gs["sign_coh"])
    # the trained VQE state (bit-identical re-run of the endpoint arm), its overlap with g;
    # anchor: F must equal the committed endpoint row's F (`s28_B_rows.jsonl`) bit for bit
    psi_abs0 = None
    for s in B.SEEDS:
        p, cv, Hn, hv, F, th, _ = B.run_hop_vqe(enc.E, A, J_FIT, seed=s, n=enc.n_qubits)
        psi = circ.state(th)
        rec[f"F_vqe_s{s}"] = float(F)
        rec[f"F_vqe_s{s}_row"] = endpoint_F.get(s)
        rec[f"anchor_bit_identical_s{s}"] = (endpoint_F.get(s) is not None
                                             and float(F) == endpoint_F[s])
        rec[f"hop_vqe_s{s}"] = float(hv)
        rec[f"coh_vqe_s{s}"] = B.sign_coherence(psi)
        rec[f"ov_vqe_s{s}"] = overlap2(psi, g)
        rec[f"ov_absvqe_s{s}"] = overlap2(np.abs(psi), g)      # the same p, aligned signs
        f_abs = F_of_state(np.abs(psi), enc.E, A, J_FIT)
        rec[f"F_absvqe_s{s}"] = f_abs["F"]                     # F the same p would reach aligned
        if s == 0:
            psi_abs0 = np.abs(psi)
    assert psi_abs0 is not None
    # the fit: target g, two budgets, the same 16 draws (seed 0; draw 1 = the VQE's theta_0)
    for tag, iters in (("short", ITERS_SHORT), ("long", iters_long)):
        fit = fit_overlap(circ, g, np.random.default_rng(0), n_starts, iters)
        psi_f = circ.state(fit["theta"])
        f_f = F_of_state(psi_f, enc.E, A, J_FIT)
        rec[f"ov_fit_{tag}"] = fit["best"]
        rec[f"ov_fit_{tag}_finals"] = fit["finals"]
        rec[f"ov_fit_{tag}_median_start"] = float(np.median(fit["finals"]))
        rec[f"F_fit_{tag}"] = f_f["F"]
        rec[f"hop_fit_{tag}"] = f_f["hop"]
        rec[f"coh_fit_{tag}"] = f_f["sign_coh"]
        rec[f"pr_fit_{tag}"] = float(1.0 / np.sum((psi_f ** 2 / np.sum(psi_f ** 2)) ** 2))
        if tag == "short":
            rec["ov_untrained_best16"] = fit["untrained_best"]
            rec["ov_untrained_mean16"] = fit["untrained_mean"]
    # the second target vector: the sign-aligned seed-0 VQE state (same p as the VQE, Perron
    # signs); if THIS is representable and g is not, the barrier is the amplitude profile
    fit2 = fit_overlap(circ, psi_abs0, np.random.default_rng(0), n_starts, iters_long)
    rec["ov_fit_absvqe_long"] = fit2["best"]
    rec["F_fit_absvqe_long"] = F_of_state(circ.state(fit2["theta"]), enc.E, A, J_FIT)["F"]
    return rec


def summarise(rows) -> Dict:
    keys = [k for k in rows[0] if isinstance(rows[0][k], (int, float)) and k not in ("J",)]
    summ = dict(n=len(rows), targets=[r["pdb"] for r in rows])
    for k in keys:
        v = np.array([r[k] for r in rows], float)
        summ[k] = dict(mean=float(v.mean()), median=float(np.median(v)), min=float(v.min()),
                       max=float(v.max()))
    return summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters-long", type=int, default=ITERS_LONG)
    ap.add_argument("--starts", type=int, default=N_STARTS)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    pick = P.targets()[::11][:12]
    if a.limit:
        pick = pick[:a.limit]
    done = B._done_pdbs(ROWS)
    t0 = time.time()
    for pdb in pick:
        if pdb in done:
            continue
        rec = target_row(pdb, a.starts, a.iters_long)
        with open(ROWS, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"  {pdb}  ov_vqe s0/s1 {rec['ov_vqe_s0']:.3f}/{rec['ov_vqe_s1']:.3f}  "
              f"ov_fit 80 {rec['ov_fit_short']:.3f}  {rec['iters_long']} {rec['ov_fit_long']:.3f}  "
              f"F_fit_long {rec['F_fit_long']:.3f}  F_gs {rec['F_gs']:.3f}  F_vqe {rec['F_vqe_s0']:.3f}  "
              f"ov_fit|vqe| {rec['ov_fit_absvqe_long']:.3f}  ({(time.time()-t0)/60:.1f} min)",
              flush=True)
    rows = [r for r in B.load_rows(ROWS) if r["pdb"] in set(pick)]
    summ = summarise(rows)
    ST.save_atomic(OUT, dict(kind="representability fit of the J = 3 ground state (S28-L26)",
                             lane="S28B", J=J_FIT, n_starts=a.starts, iters_short=ITERS_SHORT,
                             iters_long=a.iters_long, rows=rows, summary=summ),
                   module_file=__file__)
    print("\nmedian over targets (n = %d):" % len(rows))
    for k in ("gs_pr", "ov_untrained_best16", "ov_vqe_s0", "ov_vqe_s1", "ov_absvqe_s0",
              "ov_fit_short", "ov_fit_long", "ov_fit_absvqe_long",
              "F_gs", "F_vqe_s0", "F_vqe_s1", "F_absvqe_s0", "F_fit_short", "F_fit_long",
              "F_fit_absvqe_long", "hop_gs", "hop_vqe_s0", "hop_fit_long", "coh_vqe_s0",
              "coh_fit_long", "pr_fit_long"):
        s = summ[k]
        print(f"  {k:22s} median {s['median']:9.4f}  mean {s['mean']:9.4f}  [{s['min']:.4f}, {s['max']:.4f}]")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
