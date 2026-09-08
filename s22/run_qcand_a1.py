"""SPRINT 22 / WORKSTREAM A -- A1: candidate-state encoding, build + permutation/gauge GATE.

Pre-registered in `s22/PREREG_A.md` A1. Two gates, both reported with a fired-count:

  GATE 1 (soundness, EXACT).  The classical exact global argmin SET and the trained circuit's
  order-based argmin readout SET must agree on EVERY (target, permutation, seed) cell. Expected
  fired-count: 0/N (a gate that never fires is not evidence on its own -- BRIEF S6 rule 7 -- so
  the denominator N is printed beside the numerator every time).

  TIE-BREAKING HAZARD, caught in a pilot run before this number was read for the record: target
  `5H1H` carries an EXACT numerical tie between its top two distogram scores (gap 0.0). A single-
  index argmin readout is then gauge-DEPENDENT purely from tie-breaking order (a relabelling can
  flip which tied candidate a stable sort visits first) -- the project's documented
  "tie-breaking leaks the pool order" trap, reproduced in a new substrate. Both `argmin_readout`
  and `classical_exact_sort` (`s22/qcand_lib.py`) now return the FULL tied candidate SET and
  average ORACLE RMSD over it; the set is a set of CANDIDATE indices and is therefore
  gauge-invariant by construction, unlike a single arbitrarily-broken tie.

  GATE 2 (achieved gauge robustness).  Across P=8 random relabellings x 2 seeds, does the
  ACHIEVED outcome (final exact CVaR value; argmin RMSD; tail-average RMSD) exceed the canonical
  labelling's OWN seed-noise spread (4 seeds, identity permutation)? Fired = |perm outcome -
  canonical mean| > 2 x canonical seed sd, for each of the three outcome variables separately.

Runs H_DIST (primary, native-free) in full on all 16 targets; H_Legacy (secondary, pillar
separability) on a fixed half-subset (8 of 16) at reduced seed/permutation count, cheap and
declared in PREREG_A as a scope choice, not a shortcut hiding a result.

No AMBER. No OpenMM. No shared-budget bookkeeping needed (E is a lookup table).
"""
from __future__ import annotations

import math
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

ITERS = 150
ALPHA = 0.15
N_PERM = 8
SEEDS_CANON = [0, 1, 2, 3]
SEEDS_PERM = [0, 1]


def run_one(pool, which_h, entangler, seed, label=None, iters=ITERS):
    scores = pool["score_dist"] if which_h == "dist" else pool["score_leg"]
    enc = QC.Encoding(scores, label=label)
    out = QC.train(enc, alpha_schedule=ALPHA, seed=seed, iters=iters, entangler=entangler)
    _bits, cset = QC.argmin_readout(enc, out["p_final"])
    argmin_rr, _ = QC.tail_average_rmsd(pool, cset)
    face = QC.exact_face(enc.E, out["p_final"], ALPHA)
    cands = QC.tail_candidates(enc, face["tail"])
    r_tailavg, C = QC.tail_average_rmsd(pool, cands)
    return {
        "argmin_set": sorted(int(x) for x in cset), "argmin_rr": float(argmin_rr),
        "final_cvar": float(out["hist_val"][-1]), "tail_size": int(face["tail_size"]),
        "ess": float(face["ess"]), "r_tailavg": r_tailavg,
        "n_tail_cands": int(len(cands)),
    }


def main():
    tg = L20.subset(n=16)
    rows = []
    t_start = time.time()
    for ti, t in enumerate(tg):
        pdb = t["pdb"]
        pool = QC.build_pool(pdb)
        exact = QC.classical_exact_sort(pool, pool["score_dist"], ALPHA)
        n_qubits = QC.Encoding(pool["score_dist"]).n_qubits
        # ---- canonical labelling, seed-noise baseline (identity label)
        canon = [run_one(pool, "dist", "cnot", s, label=None) for s in SEEDS_CANON]
        canon_argmin_rr = np.array([r["argmin_rr"] for r in canon])
        canon_cvar = np.array([r["final_cvar"] for r in canon])
        canon_tailavg = np.array([r["r_tailavg"] for r in canon])
        canon_sd = {"argmin_rr": float(canon_argmin_rr.std(ddof=0) or 1e-9),
                    "cvar": float(canon_cvar.std(ddof=0) or 1e-9),
                    "tailavg": float(canon_tailavg.std(ddof=0) or 1e-9)}
        canon_mean = {"argmin_rr": float(canon_argmin_rr.mean()),
                     "cvar": float(canon_cvar.mean()),
                     "tailavg": float(canon_tailavg.mean())}

        # ---- GATE 1: exactness, over canonical + every permutation cell. Compared as SETS
        # (gauge-invariant) -- see the tie-breaking hazard note in the module docstring.
        exact_set = set(int(x) for x in exact["argmin_set"])
        gate1_checks, gate1_fired = 0, 0

        def check_gate1(res):
            nonlocal gate1_checks, gate1_fired
            gate1_checks += 1
            if set(res["argmin_set"]) != exact_set:
                gate1_fired += 1

        for r in canon:
            check_gate1(r)

        # ---- P random permutations x 2 seeds each; GATE 2 vs the canonical spread
        perm_rows = []
        gate2_fired = {"argmin_rr": 0, "cvar": 0, "tailavg": 0}
        gate2_checks = 0
        for k in range(N_PERM):
            label = QC.random_label(1 << n_qubits, pdb, "a1", k)
            for s in SEEDS_PERM:
                r = run_one(pool, "dist", "cnot", s, label=label)
                check_gate1(r)
                gate2_checks += 1
                r_key_of = {"argmin_rr": "argmin_rr", "cvar": "final_cvar",
                           "tailavg": "r_tailavg"}
                for mkey, rk in r_key_of.items():
                    if abs(r[rk] - canon_mean[mkey]) > 2.0 * canon_sd[mkey]:
                        gate2_fired[mkey] += 1
                perm_rows.append({"perm": k, "seed": s, **r})

        # ---- H_Legacy, secondary pillar, reduced scope (declared in PREREG_A)
        leg_row = None
        if ti % 2 == 0:  # 8 of 16 targets
            leg_canon = [run_one(pool, "leg", "cnot", s, label=None) for s in [0, 1]]
            leg_perm = []
            for k in range(4):
                label = QC.random_label(1 << n_qubits, pdb, "a1leg", k)
                leg_perm.append(run_one(pool, "leg", "cnot", 0, label=label))
            leg_row = {"canon": leg_canon, "perm": leg_perm}

        # ---- matched controls
        rnd = QC.random_subset_control(pool, exact["m"], n_draws=32, seed=1)
        sa = QC.simulated_annealing_control(pool["score_dist"], seed=0)
        p_untrained, _ = QC.untrained_distribution(n_qubits, "cnot", seed=0)
        untrained_face = QC.exact_face(QC.Encoding(pool["score_dist"]).E, p_untrained, ALPHA)
        untrained_cands = QC.tail_candidates(QC.Encoding(pool["score_dist"]), untrained_face["tail"])
        r_untrained_tailavg, _ = QC.tail_average_rmsd(pool, untrained_cands)
        _b0, c0set = QC.argmin_readout(QC.Encoding(pool["score_dist"]), p_untrained)
        r_untrained_argmin, _ = QC.tail_average_rmsd(pool, c0set)

        rows.append({
            "pdb": pdb, "n": pool["n"], "fold": pool["fold"], "K": enc_k(pool),
            "n_qubits": n_qubits, "pool_oracle_best": float(pool["rr"].min()),
            "R_score_argmin": exact["r_argmin"], "R_score_tailavg": exact["r_tailavg"],
            "R_pool": float(pool["rr"].min()),
            "canon_mean": canon_mean, "canon_sd": canon_sd,
            "canon_rows": canon,
            "gate1_checks": gate1_checks, "gate1_fired": gate1_fired,
            "gate2_checks": gate2_checks, "gate2_fired": gate2_fired,
            "perm_rows": perm_rows,
            "legacy": leg_row,
            "random_control": rnd, "sa_control": sa,
            "untrained_argmin_rr": float(r_untrained_argmin),
            "untrained_tailavg_rr": r_untrained_tailavg,
            "untrained_tail_size": int(untrained_face["tail_size"]),
        })
        print(f"[{ti+1}/16] {pdb} n={pool['n']} nq={n_qubits} gate1 {gate1_fired}/{gate1_checks} "
              f"gate2 {gate2_fired} R_score_argmin={exact['r_argmin']:.3f} "
              f"R_score_tailavg={exact['r_tailavg']:.3f} elapsed={time.time()-t_start:.0f}s",
              flush=True)

    total_gate1_fired = sum(r["gate1_fired"] for r in rows)
    total_gate1_checks = sum(r["gate1_checks"] for r in rows)
    total_gate2_fired = {k: sum(r["gate2_fired"][k] for r in rows) for k in
                         ("argmin_rr", "cvar", "tailavg")}
    total_gate2_checks = sum(r["gate2_checks"] for r in rows)
    out = {
        "experiment": "A1_candidate_encoding_gauge",
        "alpha": ALPHA, "iters": ITERS, "n_perm": N_PERM,
        "seeds_canon": SEEDS_CANON, "seeds_perm": SEEDS_PERM,
        "n_targets": len(rows),
        "gate1_fired_total": total_gate1_fired, "gate1_checks_total": total_gate1_checks,
        "gate2_fired_total": total_gate2_fired, "gate2_checks_total": total_gate2_checks,
        "rows": rows,
        "wall_time_s": time.time() - t_start,
    }
    req = ["pdb", "n_qubits", "R_score_argmin", "R_score_tailavg", "canon_mean", "canon_sd",
          "gate1_checks", "gate1_fired", "gate2_checks", "gate2_fired", "perm_rows",
          "random_control", "sa_control", "untrained_argmin_rr", "untrained_tailavg_rr"]
    path = QC.write("a1_gauge", out, required_keys=req)
    print("WROTE", path)
    print(f"GATE 1 (exactness): {total_gate1_fired}/{total_gate1_checks} fired "
          f"(expect 0/{total_gate1_checks})")
    print(f"GATE 2 (achieved gauge): {total_gate2_fired} fired of {total_gate2_checks} checks each")


def enc_k(pool):
    return len(pool["score_dist"])


if __name__ == "__main__":
    main()
