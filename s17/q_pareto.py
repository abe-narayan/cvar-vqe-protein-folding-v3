"""SPRINT 17 / QUANTUM -- T1, THE PARETO TEST (sprint section 20), and T5 mixtures (19).

PRE-REGISTERED IN `s17/PREREG_quantum.md` BEFORE ANY ARM RAN.  P1 and P5 are the rules this
module fires.

HYPOTHESIS (P1).  No CVaR-VQE configuration reaches an (M, D) point that the classical
attainable set -- a fixed-temperature Metropolis thermostat ladder at MATCHED objective
budget, plus greedy 1-opt, plus annealing, plus the retrieval product law, plus uniform --
does not already dominate.

EXPECTED OUTCOME.  Dominated everywhere.  Sprint 16 section 6b showed a diversity-matched
classical thermostat reproduces alpha's whole ensemble effect (Pearson +0.93/+0.98) and moves
the readout 1.6-1.7x FURTHER; a thermostat that can move support at matched budget should be
strictly stronger than the `tilt_exact` used there.

STRONGEST CONTROL.  `search_metropolis` at matched budget, swept over a temperature ladder so
the classical arm traces a CURVE.  Sprint 16 never ran this: `tilt_samples` cannot leave its
support and `tilt_exact` costs 32-128x the budget.

SUCCESS CRITERION.  At least one VQE configuration whose (M, D) point is not dominated by any
classical point at matched budget, on >= 10 of 19 targets, with the paired per-target
dominance margin's CI excluding zero.

FALSIFIER.  If no VQE configuration reaches a frontier point the classical thermostat plus
1-opt do not already dominate at matched budget, the quantum pillar has no accuracy role in
this architecture and is retained as a studied object only.

WHAT IS ORACLE HERE.  `M`, `readout`, `set_best` and the top-0.1% novelty sets read the
native and are post-hoc scoring only.  `D`, distinct counts and every arm's parameters are
native-free.  No temperature, alpha or stopping rule is chosen with a native quantity.

RUN:  python -m s17.q_pareto {verify,run,report}
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from s12 import instrument as I
from s14 import vqe_lib as V
from s14 import vqe_run as R
from s15 import seed as SD
from s16 import qphase_lib as QP
from s17 import q_lib as L

SEEDS = (0, 1, 2)
MS = (5, 20, 75)
NREP = {5: 8, 20: 6, 75: 4}

#: temperature ladder in units of the objective's own sd -- NATIVE-FREE
TLADDER = (0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 12.8)
#: the retrieval product law's own temperature ladder, in units of the prior's sd
PLADDER = (0.1, 0.3, 1.0, 3.0, 30.0)

VQE_ARMS = (
    ("vqe_a0.02", dict(alpha=0.02)),
    ("vqe_a0.05", dict(alpha=0.05)),
    ("vqe_a0.10", dict(alpha=0.10)),
    ("vqe_a0.25", dict(alpha=0.25)),
    ("vqe_a0.50", dict(alpha=0.50)),
    ("vqe_a1.00", dict(alpha=1.00)),                       # = plain VQE (no CVaR tail)
    ("vqe_sched_hi_lo", dict(alpha=(1.00, 0.05), alpha_anneal=True)),
    ("vqe_sched_lo_hi", dict(alpha=(0.05, 1.00), alpha_anneal=True)),
    ("vqe_warm_prior", dict(alpha=0.25, init="prior")),    # retrieval-informed warm start
    ("vqe_deep_mps3f", dict(alpha=0.25, ansatz="mps3f")),
)

#: THE CLASSICAL ATTAINABLE SET.  Every member is at or below the matched budget and none
#: touches a quantum circuit.  `untrained`, `tilt_samples` and `tilt_exact` are excluded on
#: purpose: the first two draw from the circuit and the third reads the whole register.
CLASSICAL_MATCHED = tuple(
    [f"metro_T{t}" for t in TLADDER]
    + [f"metro1opt_T{t}" for t in (0.05, 0.2, 0.8)]
    + ["greedy", "anneal", "anneal_q", "random"]
    + [f"retrieval_T{t}" for t in PLADDER])
QUANTUM_ARMS = tuple(n for n, _ in VQE_ARMS)
#: reported, never counted as classical
CIRCUIT_DERIVED = ("untrained", "tilt_samples", "tilt_exact")


# ------------------------------------------------------------------ set readouts
def set_points(ins, seen, rng, E=None, ms=MS, nrep=NREP):
    """The arm's (M, D) points at each operational set size, in BOTH consumption modes.

    `M<m>` / `D<m>` -- the set is drawn UNIFORMLY from the arm's own VISITED MULTISET, so an
    arm that revisits configurations gets the lower diversity it earned.  This is the
    sampler's own geometry, with no ranker anywhere.

    `selM<m>` / `selD<m>` -- the `m` DISTINCT visited configurations with the lowest
    objective.  This is what the shipped pipeline actually consumes (a top-75 cut), it is
    deterministic given the visited set, and it is the readout the sprint's architectural
    pivot is about.  Reporting only the first mode would answer a question nobody asks.

    Every quantity is a quadratic mean in the common frame (`L.md_plane`).
    """
    seen = np.asarray(seen, np.int64)
    out = {}
    if E is not None:
        u = np.unique(seen)
        o = u[np.argsort(np.asarray(E, float)[u], kind="mergesort")]
        for m in ms:
            r = L.md_plane(ins.ca(o[:m]), ins.nat)
            for k in ("M", "D", "readout", "M_free", "set_best"):
                out[f"sel{k}{m}"] = float(r[k])
    for m in ms:
        acc = {k: [] for k in ("M", "D", "readout", "M_free", "set_best")}
        for _ in range(nrep[m]):
            pick = rng.choice(seen, size=m, replace=seen.size < m)
            r = L.md_plane(ins.ca(pick), ins.nat)
            for k in acc:
                acc[k].append(r[k])
        for k, v in acc.items():
            # M and D enter the identity QUADRATICALLY: average the squares, not the roots
            if k in ("M", "D", "readout"):
                out[f"{k}{m}"] = float(np.sqrt(np.mean(np.square(v))))
            else:
                out[f"{k}{m}"] = float(np.mean(v))
    return out


def arm_readout(ins, E, seen, rng, cost, e_pct=None):
    seen = np.asarray(seen, np.int64)
    u = np.unique(seen)
    Ee = np.asarray(E, float)
    es = Ee[seen]
    mn = es.min()
    tied = np.unique(seen[es == mn])
    out = {"n_distinct": int(u.size),
           "draw_entropy_bits": QP.draw_entropy_bits(seen, ins.N),
           "best_e": float(mn),
           "objective_gap": float(mn - Ee.min()),
           "rmsd_returned": float(ins.rmsd[tied].mean()),      # ORACLE, tie-averaged
           "argmin_ties": int(tied.size),
           "set_mean_rmsd": float(ins.rmsd[seen].mean()),      # ORACLE
           "set_best_rmsd": float(ins.rmsd[u].min()),          # ORACLE
           "cost": cost}
    out.update(set_points(ins, seen, rng, E=Ee))
    return out, u


# --------------------------------------------------------------------- the arms
def run_cell(pdb, seed, targets_top=0.001):
    ins = L.inst(pdb)
    E = ins.hamil()
    sd = float(E.std())
    rmsd = ins.rmsd
    thr = float(np.quantile(rmsd, targets_top))
    elite = set(np.flatnonzero(rmsd <= thr).tolist())        # ORACLE, for P5 novelty only
    rng = SD.stable_rng(pdb, seed, "readout", salt=L.SALT)
    res = {}
    supports = {}

    def add(name, seen, cost):
        r, u = arm_readout(ins, E, seen, rng, cost)
        r["elite_hits"] = int(len(elite.intersection(u.tolist())))
        res[name] = r
        supports[name] = u

    # ---- quantum arms -----------------------------------------------------
    for name, kw in VQE_ARMS:
        c = L.Cost()
        prior = None
        if kw.get("init") == "prior":
            f, _, _ = L.prior_factors(ins)
            z = -f / max(f.std(), 1e-9)
            z -= z.max(1, keepdims=True)
            prior = np.exp(z); prior /= prior.sum(1, keepdims=True)
        v = R.run(E, ins.n_qubits, kw.get("alpha", 0.25), L.BUDGET, shots=L.SHOTS,
                  ansatz=kw.get("ansatz", L.ANSATZ), seed=seed, rmsd=rmsd,
                  bits_per_res=ins.bits_per_res,
                  # `exact_dist` enumerates the whole register to report the trained
                  # distribution's mode and entropy.  It runs AFTER training, consumes no
                  # randomness and touches no arm's sampled set, so switching it off for the
                  # arms whose exact distribution is not used changes NO measured quantity --
                  # only the wall clock.  Kept on for the reference arm (alpha = 0.25).
                  exact_dist=(kw.get("alpha") == 0.25 and not kw.get("alpha_anneal")),
                  keep_seen=True,
                  lr=0.15, baseline="const", alpha_anneal=kw.get("alpha_anneal", False),
                  init=kw.get("init", "random"), prior=prior)
        seen = np.asarray(v["seen"], np.int64)
        cost = c.stop(v["evals"], circuit_samples=v["evals"], grad_passes=v["iters"],
                      distinct=int(np.unique(seen).size))
        add(name, seen, cost)
        for k in ("entropy_bits", "mode_rmsd", "mean_rmsd_under_p", "max_prob",
                  "init_rmsd_best", "param_displacement"):
            if k in v:
                res[name][k] = float(v[k])

    # ---- classical thermostat ladder, MATCHED BUDGET ----------------------
    for t in TLADDER:
        c = L.Cost()
        cc = L.search_metropolis(E, ins.n, ins.k, L.BUDGET, t * sd,
                                 SD.stable_rng(pdb, seed, f"metro{t}", salt=L.SALT))
        s = cc.all_seen()
        add(f"metro_T{t}", s, c.stop(cc.used, distinct=int(np.unique(s).size)))
    # ---- thermostat + 1-opt, the falsifier's named arm --------------------
    for t in (0.05, 0.2, 0.8):
        c = L.Cost()
        s, used = L.search_metropolis_1opt(
            E, ins.n, ins.k, L.BUDGET, t * sd,
            SD.stable_rng(pdb, seed, f"m1opt{t}", salt=L.SALT))
        add(f"metro1opt_T{t}", s, c.stop(used, distinct=int(np.unique(s).size)))

    # ---- the other classical searches, matched and quarter budget ---------
    for name, fn, b in (("greedy", V.search_greedy, L.BUDGET),
                        ("anneal", V.search_anneal, L.BUDGET),
                        ("anneal_q", V.search_anneal, L.BUDGET // 4),
                        ("random", V.search_random, L.BUDGET)):
        c = L.Cost()
        cc = fn(E, ins.n, ins.k, b, SD.stable_rng(pdb, seed, name, salt=L.SALT))
        s = cc.all_seen()
        add(name, s, c.stop(cc.used, distinct=int(np.unique(s).size)))

    # ---- the untrained circuit (weak control, kept for continuity) --------
    c = L.Cost()
    ud = QP.untrained_draws(ins.n_qubits, L.BUDGET, seed)
    add("untrained", ud, c.stop(L.BUDGET, circuit_samples=L.BUDGET,
                                distinct=int(np.unique(ud).size)))

    # ---- Sprint 16's two thermostats, for continuity ----------------------
    h = res["vqe_a0.25"]["draw_entropy_bits"]
    c = L.Cost()
    ti = QP.tilt_samples(ud, E, h, SD.stable_rng(pdb, seed, "tilt", salt=L.SALT))
    add("tilt_samples", ti, c.stop(L.BUDGET, distinct=int(np.unique(ti).size)))
    c = L.Cost()
    p0 = QP.untrained_probs(ins.n_qubits, seed)
    ti, _ = QP.tilt_exact(p0, E, res["vqe_a0.25"]["n_distinct"], L.BUDGET,
                          SD.stable_rng(pdb, seed, "tiltx", salt=L.SALT))
    add("tilt_exact", ti, c.stop(ins.N, distinct=int(np.unique(ti).size),
                                 register_reads=1))

    # ---- the retrieval distribution itself, ZERO objective evaluations ----
    for t in PLADDER:
        c = L.Cost()
        pr_sd = float(ins.prior.std())
        s = L.prior_draws(ins, t * pr_sd, L.BUDGET,
                          SD.stable_rng(pdb, seed, f"prior{t}", salt=L.SALT))
        add(f"retrieval_T{t}", s, c.stop(0, distinct=int(np.unique(s).size)))

    # ---- P5: novelty of the elite (ORACLE) ---------------------------------
    # `CLASSICAL_MATCHED` is the honest classical attainable set: no arm in it touches a
    # quantum circuit and none exceeds the matched budget.  `untrained`, `tilt_samples` and
    # `tilt_exact` are circuit-derived and/or over budget, so they are reported but never
    # counted as classical.
    cl_union = set()
    for k in supports:
        if k in CLASSICAL_MATCHED or k.startswith(("metro", "retrieval")):
            cl_union |= set(supports[k].tolist())
    sets = {k: set(v.tolist()) for k, v in supports.items()}
    for k in supports:
        own = sets[k]
        others = set()
        for j in supports:
            if j != k:
                others |= sets[j]
        # LEAVE-ONE-OUT novelty: symmetric across arms, so an arm is not credited or
        # penalised merely for how many sibling arms its family happens to have.
        res[k]["elite_unique_loo"] = int(len(elite & own - others))
        if k.startswith("vqe"):
            res[k]["elite_unique_vs_classical"] = int(len(elite & own - cl_union))
    q_union = set()
    for k in supports:
        if k.startswith("vqe"):
            q_union |= sets[k]
    res["_novelty"] = {
        "n_elite": int(len(elite)), "elite_rmsd_thr": thr,
        "elite_found_by_quantum_union": int(len(elite & q_union)),
        "elite_found_by_classical_union": int(len(elite & cl_union)),
        "elite_only_quantum": int(len(elite & q_union - cl_union)),
        "elite_only_classical": int(len(elite & cl_union - q_union)),
        "n_quantum_union": int(len(q_union)), "n_classical_union": int(len(cl_union)),
        "quantum_support_subset_of_classical": bool(len(q_union - cl_union) == 0),
        "frac_quantum_support_novel": float(len(q_union - cl_union) / max(1, len(q_union))),
    }
    res["_meta"] = {"pdb": pdb, "seed": seed, "n": ins.n, "k": ins.k, "N": ins.N,
                    "fold": ins.fold, "E_sd": sd,
                    "rmsd_best_in_space": float(rmsd.min()),
                    "rmsd_mean_in_space": float(rmsd.mean()),
                    "certified_argmin_rmsd": float(rmsd[E == E.min()].mean())}
    return res


def run(targets=L.TARGETS19, seeds=SEEDS):
    """SEEDS OUTER, TARGETS INNER.  On a contended box an interrupted run must leave a
    COMPLETE sweep at the statistical unit, not a partial one: after 19 cells there are 19
    targets at one seed (n = 19, the unit of analysis), not 6 targets at three seeds."""
    done = L.ck_load("pareto")
    for sd in seeds:
        for pdb in targets:
            key = f"{pdb}_{sd}"
            if key in done:
                print(f"  skip {key}", flush=True)
                continue
            L.gate(key)
            t0 = time.time()
            r = run_cell(pdb, sd)
            L.ck("pareto", key, r)
            done[key] = r
            print(f"  {key}  {time.time()-t0:.0f}s  "
                  f"vqe25 M/D {r['vqe_a0.25']['M75']:.3f}/{r['vqe_a0.25']['D75']:.3f}  "
                  f"metroT0.2 {r['metro_T0.2']['M75']:.3f}/{r['metro_T0.2']['D75']:.3f}",
                  flush=True)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "verify":
        v = L.verify_identity()
        L.ck("verify", "identity", v)
        print(json.dumps({k: v[k] for k in ("max_abs_identity_residual_A2", "PASS")},
                         indent=1))
        ins = L.inst("1CS9")
        f, c, resid = L.prior_factors(ins)
        print("prior additivity residual:", resid)
        # budget honesty of the new thermostat
        E = ins.hamil()
        cc = L.search_metropolis(E, ins.n, ins.k, 8192, 0.2 * E.std(),
                                 SD.stable_rng("1CS9", 0, "v", salt=L.SALT))
        print("metropolis used", cc.used, "of 8192; distinct",
              int(np.unique(cc.all_seen()).size))
        L.ck("verify", "metropolis_budget", {"used": int(cc.used), "budget": 8192})
    elif mode == "run":
        tg = sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS19
        run(tg)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
