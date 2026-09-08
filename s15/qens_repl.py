"""SPRINT 15 / QENS -- REPLICATION AT POWER of the project's only positive quantum result.

THE CLAIM UNDER TEST (`s15/qgeom_FINDINGS.md` E2, relayed as coord_FINDINGS Q3.1):

    Consumed as an ENSEMBLE with NO ranker anywhere, a CVaR-VQE beats best-of-N from its own
    UNTRAINED circuit by 0.36-0.57 A on a coordinate-average readout, CIs excluding zero, at
    alpha = 1 and alpha = 0.05.  n = 27 (9 targets x 3 seeds).  Flagged LOW POWER.

FOUR CELLS.

R0  BIT-LEVEL REPRODUCTION.  The original 27 cells re-run through this module's Struct and
    the ORIGINAL `qgeom_ens.unranked_readouts` function, compared to the recorded JSON.  This
    establishes that the replication machinery is the same machinery before it is changed.

R1  THE REPLICATION AT POWER.  19 enumerated targets x 8 seeds = 152 cells per alpha
    against the original 27 (9 targets x 3 seeds), at alpha in {1.0, 0.25, 0.05}, with SIX
    arms:

      vqe          CVaR-VQE final distribution
      untrained    the same circuit at its untrained theta          [the incumbent control]
      uniform      uniform random over the sub-register             [WEAK -- see the brief]
      anneal       classical SA, 2,048 objective evaluations        [budget class A]
      anneal_cost  classical SA, 819,200 objective evaluations      [budget class B]
      tilt         classical Boltzmann tilt of p_init on the same objective, temperature
                   solved so its entropy equals the VQE's           [budget class B, diagnostic]

    BUDGET ACCOUNTING, stated because the original convention flatters the VQE.  The
    exact-gradient VQE reads ALL 4,096 objective values on each of its 200 iterations, i.e.
    **819,200 objective evaluations**, while the original control was charged **2,048 draws**.
    That is a 400x advantage, and it is the convention the original result was measured under.
    Class A reproduces it (equal readout draws).  Class B charges the classical arms the
    VQE's real evaluation count.  Both are reported; if they disagree, the disagreement is
    the result.

    UNIT OF ANALYSIS.  The original treated 27 (target, seed) cells as 27 independent units.
    Seeds within a target are not independent replicates of a TARGET effect, so this module
    reports BOTH: the cell-level test (n = 152) and the target-level test (n = 19, seeds
    averaged first).  The target-level test is the primary one.

R2  THE MATCHED-DIVERSITY DECOMPOSITION.  `mdN_coordavg` subsamples every arm's DISTINCT
    configuration set to the common size `D = min(distinct(vqe), distinct(untrained))` and
    averages N members drawn uniformly from it.  Equal diversity, equal budget.  If the VQE
    wins there the effect is about WHICH configurations it favours; if not, the effect is
    diversity, which the untrained circuit supplies more cheaply.

R3  THE ALPHA MAP.  Ten alphas x 19 targets x 4 seeds, vqe and untrained only, to say whether
    the reported non-monotonicity (a win at alpha=1 and alpha=0.05, nothing at 0.25) is real
    or two noisy cells either side of nothing.

ORACLE.  Every RMSD, every `set_*`/`rand*`/`md*` readout and every diversity-vs-native ratio
is computed from native coordinates POST HOC and is ORACLE.  No native quantity enters any
objective, any circuit, any temperature solve or any arm's sampling.

    python -m s15.qens_repl            # R0, R1/R2, R3 in order; resumable
    python -m s15.qens_repl R1         # one cell
"""
from __future__ import annotations

import sys
import time

import numpy as np

from s15 import qens_lib as Q
from s15 import qgeom_lib as G
from s15 import seed as SD
from s15.qgeom_cvar import make_circ, run_vqe

TAG = "repl"
ALPHAS_MAIN = (1.0, 0.25, 0.05)
ALPHAS_MAP = (1.0, 0.75, 0.50, 0.35, 0.25, 0.15, 0.10, 0.05, 0.025, 0.01)
SEEDS_MAIN = tuple(range(8))
SEEDS_MAP = tuple(range(4))
BUDGET = 2048
ITERS = 200
ARMS = ("vqe", "untrained", "uniform", "anneal", "anneal_cost", "tilt")

# the VQE's true objective-evaluation count: an exact statevector gradient reads every one
# of the M configuration energies at every iteration.
def vqe_evals(M, iters=ITERS):
    return int(M * iters)


# ============================================================================== R0
def r0_reproduce(pdbs=Q.TARGETS9, alphas=(1.0, 0.25, 0.05, 0.01), seeds=(0, 1, 2)):
    """Bit-level reproduction of the original 27-cell E2 table through this module's Struct."""
    import json
    import os
    from s15.qgeom_ens import unranked_readouts
    print("=" * 118)
    print("R0. BIT-LEVEL REPRODUCTION of qgeom_ens.E2 (original machinery, this module's Struct)")
    print("=" * 118)
    ref = json.load(open(os.path.join(G.RESULTS, "qgeom_ens.json")))["E2_generator"]
    worst = {}
    rows = []
    for pdb in pdbs:
        st = Q.Struct(pdb)
        E = Q.hamil_sub(st)
        circ = make_circ(st.n_qubits)
        for a in alphas:
            for s in seeds:
                key = f"{pdb}|a{a}|s{s}"
                if key not in ref:
                    continue
                rng = np.random.default_rng(1000 * s + 7)
                r = run_vqe(circ, E, a, ITERS, seed=s, lowmem=True)
                ci = rng.choice(len(E), size=BUDGET, p=r["p_init"])
                vi = rng.choice(len(E), size=BUDGET, p=r["p_final"])
                got = {"control": unranked_readouts(st, ci, rng=np.random.default_rng(s)),
                       "vqe": unranked_readouts(st, vi, rng=np.random.default_rng(s))}
                for side in ("control", "vqe"):
                    for kk, vv in got[side].items():
                        if kk not in ref[key][side]:
                            continue
                        d = abs(float(vv) - float(ref[key][side][kk]))
                        worst[kk] = max(worst.get(kk, 0.0), d)
                rows.append(key)
        print(f"  {pdb} done ({len(rows)} cells)", flush=True)
    print()
    print(f"  {len(rows)} cells reproduced.  Worst absolute difference per readout:")
    for kk in sorted(worst):
        print(f"    {kk:28s} {worst[kk]:.3e}")
    mx = max(worst.values()) if worst else float("nan")
    print(f"\n  MAX over every readout and every cell: {mx:.3e}")
    Q.ck(TAG, "R0_reproduce", {"n_cells": len(rows), "worst_per_readout": worst,
                               "max_abs_diff": float(mx)})
    return worst


# ============================================================================== R1/R2
AI_ARMS = ("untrained", "uniform", "anneal", "anneal_cost")     # alpha-independent
AD_ARMS = ("vqe", "tilt")                                      # alpha-dependent
ANNEAL_COST_SEEDS = 4          # the equal-evaluation SA control runs on seeds 0..3 only


def _alpha_indep(st, E, s, circ):
    """The four arms that do not depend on alpha, drawn once per (target, seed).

    `untrained` is the same circuit at `random_theta(circ, seed)`, i.e. exactly the control
    the original used.  `anneal_cost` is charged the VQE's true objective-evaluation count
    and its ensemble is the LAST `BUDGET` evaluations of the annealing trajectory -- the
    converged annealer's sampling regime, which is the analogue of the VQE's FINAL
    distribution.  Charging it the whole trajectory would let its hot phase masquerade as
    ensemble diversity.
    """
    M = len(E)
    rng = SD.stable_rng("qens_ai", st.pdb, s)
    p0 = circ.probs(G.random_theta(circ, s))
    d = {"untrained": rng.choice(M, size=BUDGET, p=p0),
         "uniform": rng.integers(0, M, BUDGET)}
    n_res, k = len(st.res), st.e.k
    idx_a, used_a = Q.anneal_ensemble(E, n_res, k, BUDGET, rng)
    d["anneal"] = idx_a[-BUDGET:]
    ev = {"untrained": 0, "uniform": 0, "anneal": int(used_a)}
    # `anneal_cost` costs 819,200 single-move evaluations, ~60% of this module's whole
    # runtime, and it is a CONTROL arm rather than a headline.  It is therefore run on
    # SEEDS 0-3 ONLY (19 targets x 4 seeds = 76 paired cells).  Every cell that carries it
    # is still exactly paired with the VQE cell of the same (target, alpha, seed); cells
    # without it are simply absent from that one comparison and its n is printed.
    if int(s) < ANNEAL_COST_SEEDS:
        nev = vqe_evals(M)
        idx_b, used_b = Q.anneal_ensemble(E, n_res, k, nev, rng)
        d["anneal_cost"] = idx_b[-BUDGET:]
        ev["anneal_cost"] = int(used_b)
    meta = {"init_entropy_bits": float(-(p0[p0 > 0] * np.log2(p0[p0 > 0])).sum()),
            "Ep_rmsd_init": float(p0 @ st.rmsd),                 # ORACLE, post-hoc
            "evals": ev}
    return d, meta, p0


def _alpha_dep(st, E, a, s, circ, p0):
    M = len(E)
    rng = SD.stable_rng("qens_ad", st.pdb, s, a)
    r = run_vqe(circ, E, a, ITERS, seed=s, lowmem=True)
    d = {"vqe": rng.choice(M, size=BUDGET, p=r["p_final"])}
    pt = Q.boltzmann_tilt(p0, E, r["final"]["entropy_bits"])
    d["tilt"] = rng.choice(M, size=BUDGET, p=pt)
    meta = {"entropy_bits": r["final"]["entropy_bits"], "max_p": r["final"]["max_p"],
            "Ep_rmsd_vqe": float(r["p_final"] @ st.rmsd),        # ORACLE, post-hoc
            "tilt_entropy_bits": float(-(pt[pt > 0] * np.log2(pt[pt > 0])).sum()),
            "evals": {"vqe": vqe_evals(M), "tilt": int(M)}}
    return d, meta


def r1_power(pdbs=Q.TARGETS19, alphas=ALPHAS_MAIN, seeds=SEEDS_MAIN):
    print()
    print("=" * 118)
    print(f"R1/R2. REPLICATION AT POWER -- {len(pdbs)} targets x {len(seeds)} seeds "
          f"= {len(pdbs)*len(seeds)} cells per alpha (original: 9 x 3 = 27)")
    print("=" * 118)
    prev = Q.ck_load(TAG)
    out = prev.get("R1_cells", {})
    base = prev.get("R1_base", {})
    t00 = time.time()
    for pdb in pdbs:
        need = [(a, s) for a in alphas for s in seeds if f"{pdb}|a{a}|s{s}" not in out]
        if not need:
            print(f"  {pdb} cached", flush=True)
            continue
        G.wait_mem(0.6, "qens_repl")
        st = Q.Struct(pdb)
        E = Q.hamil_sub(st)
        circ = make_circ(st.n_qubits)
        cache = {}
        for a, s in need:
            t0 = time.time()
            if s not in cache:
                di, mi, p0 = _alpha_indep(st, E, s, circ)
                ro = {arm: Q.readouts(st, di[arm], SD.stable_rng("qens_ro", pdb, s))
                      for arm in AI_ARMS if arm in di}
                cache[s] = (di, mi, p0, ro)
                base[f"{pdb}|s{s}"] = {"meta": mi, **ro}
                Q.ck(TAG, "R1_base", base)
            di, mi, p0, ro_i = cache[s]
            dd, md = _alpha_dep(st, E, a, s, circ, p0)
            dm = min(int(np.unique(dd["vqe"]).size), int(np.unique(di["untrained"]).size))
            cell = {"meta": {**mi, **md, "evals": {**mi["evals"], **md["evals"]}},
                    "dmatch": dm}
            for arm in AD_ARMS:
                cell[arm] = Q.readouts(st, dd[arm], SD.stable_rng("qens_ro", pdb, s))
            for arm in AI_ARMS:
                if arm in ro_i:
                    cell[arm] = dict(ro_i[arm])
            for arm in AD_ARMS + AI_ARMS:
                if arm not in cell:
                    continue
                src = dd[arm] if arm in AD_ARMS else di[arm]
                cell[arm].update(Q.md_readouts(st, src,
                                               SD.stable_rng("qens_md", pdb, s, a), dm))
            cell["secs"] = time.time() - t0
            out[f"{pdb}|a{a}|s{s}"] = cell
            Q.ck(TAG, "R1_cells", out)
        print(f"  {pdb} done  [{len(out)} cells, {time.time()-t00:.0f}s]", flush=True)
    return out


# ============================================================================== R3
def r3_alpha_map(pdbs=Q.TARGETS19, alphas=ALPHAS_MAP, seeds=SEEDS_MAP):
    print()
    print("=" * 118)
    print(f"R3. THE ALPHA MAP -- {len(alphas)} alphas x {len(pdbs)} targets x "
          f"{len(seeds)} seeds, vqe vs untrained only")
    print("=" * 118)
    prev = Q.ck_load(TAG)
    out = prev.get("R3_cells", {})
    t00 = time.time()
    for pdb in pdbs:
        need = [(a, s) for a in alphas for s in seeds if f"{pdb}|a{a}|s{s}" not in out]
        if not need:
            print(f"  {pdb} cached", flush=True)
            continue
        G.wait_mem(0.6, "qens_repl")
        st = Q.Struct(pdb)
        E = Q.hamil_sub(st)
        circ = make_circ(st.n_qubits)
        M = len(E)
        cache = {}
        for a, s in need:
            if s not in cache:                     # the untrained arm is alpha-independent
                rng0 = SD.stable_rng("qens_map_ai", pdb, s)
                p0 = circ.probs(G.random_theta(circ, s))
                ci = rng0.choice(M, size=BUDGET, p=p0)
                cache[s] = (ci, Q.readouts(st, ci, SD.stable_rng("qens_ro", pdb, s),
                                           light=True),
                            float(p0 @ st.rmsd))   # ORACLE, post-hoc
            ci, ro_c, ep0 = cache[s]
            rng = SD.stable_rng("qens_map", pdb, s, a)
            r = run_vqe(circ, E, a, ITERS, seed=s, lowmem=True)
            vi = rng.choice(M, size=BUDGET, p=r["p_final"])
            dm = min(int(np.unique(vi).size), int(np.unique(ci).size))
            rv = Q.readouts(st, vi, SD.stable_rng("qens_ro", pdb, s), light=True)
            rv.update(Q.md_readouts(st, vi, SD.stable_rng("qens_md", pdb, s, a), dm))
            rc = dict(ro_c)
            rc.update(Q.md_readouts(st, ci, SD.stable_rng("qens_md", pdb, s, a), dm))
            out[f"{pdb}|a{a}|s{s}"] = {
                "vqe": rv, "untrained": rc,
                "entropy_bits": r["final"]["entropy_bits"],
                "Ep_rmsd_vqe": float(r["p_final"] @ st.rmsd),     # ORACLE, post-hoc
                "Ep_rmsd_init": ep0,                              # ORACLE, post-hoc
                "dmatch": dm}
            Q.ck(TAG, "R3_cells", out)
        print(f"  {pdb} done  [{len(out)} cells, {time.time()-t00:.0f}s]", flush=True)
    return out


# ============================================================================== R4
def r4_objective_axis(pdbs=Q.TARGETS19, alphas=ALPHAS_MAIN, seeds=SEEDS_MAIN):
    """THE OBJECTIVE AXIS, reported separately from the structural one as the brief demands.

    For every cell: the expected objective value under the VQE's final distribution and under
    the untrained one, the probability the final distribution puts on the objective's exact
    global minimiser, and the same quantities for the classical annealer at both budgets.
    The objective is uniformised to [0, 1] by rank, so `E_p[E]` is directly readable as the
    mean rank-percentile of what the distribution samples and 0 is the certified optimum.
    """
    print()
    print("=" * 118)
    print("R4. THE OBJECTIVE AXIS (rank-uniformised: 0 = the certified global optimum,")
    print("    0.5 = a uniform random draw).  Reported SEPARATELY from the structural axis.")
    print("=" * 118)
    prev = Q.ck_load(TAG)
    out = prev.get("R4_cells", {})
    for pdb in pdbs:
        need = [(a, s) for a in alphas for s in seeds if f"{pdb}|a{a}|s{s}" not in out]
        if not need:
            continue
        st = Q.Struct(pdb)
        E = Q.hamil_sub(st)
        circ = make_circ(st.n_qubits)
        M = len(E)
        xstar = int(np.argmin(E))
        cache = {}
        for a, s in need:
            if s not in cache:
                rng = SD.stable_rng("qens_ai", pdb, s)
                p0 = circ.probs(G.random_theta(circ, s))
                ia, ua = Q.anneal_ensemble(E, len(st.res), st.e.k, BUDGET, rng)
                cache[s] = (p0, float(E[ia[-BUDGET:]].mean()), float(E[ia].min()))
            p0, e_ann, e_ann_best = cache[s]
            r = run_vqe(circ, E, a, ITERS, seed=s, lowmem=True)
            p = r["p_final"]
            out[f"{pdb}|a{a}|s{s}"] = {
                "Ep_E_vqe": float(p @ E), "Ep_E_init": float(p0 @ E),
                "Ep_E_anneal_tail": e_ann, "min_E_anneal": e_ann_best,
                "p_argmin": float(p[xstar]),
                "mode_E": float(E[int(np.argmax(p))]),
                "certified_optimum": float(E[xstar]),
                "Ep_rmsd_vqe": float(p @ st.rmsd),        # ORACLE, post-hoc
                "mode_rmsd": float(st.rmsd[int(np.argmax(p))])}   # ORACLE, post-hoc
            Q.ck(TAG, "R4_cells", out)
        print(f"  {pdb} done [{len(out)} cells]", flush=True)
    return out


def main():
    G.wait_mem(0.8, "qens_repl")
    Q.ck_load(TAG)
    which = sys.argv[1:] or ["R0", "R1", "R3"]
    if "R4" in which:
        r4_objective_axis()
    if "R0" in which:
        r0_reproduce()
    if "R1" in which:
        r1_power()
    if "R3" in which:
        r3_alpha_map()
    print("\nwritten -> s15/results/qens_repl.json")


if __name__ == "__main__":
    main()
