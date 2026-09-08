"""SPRINT 17 / QUANTUM -- T3, THEORY (sprint section 22): is this problem class structurally
hostile to a variational sampler, and is there a reformulation in which it is not?

PRE-REGISTERED IN `s17/PREREG_quantum.md` (P6).  The two questions are asked SEPARATELY and
never conflated.

WHAT IS MEASURED, EXACTLY, NOT SWEPT.

1. THE WALSH (PAULI-Z) SPECTRUM OF THE DEPLOYED OBJECTIVE, resolved by Pauli weight, on the
   full enumerated register.  A diagonal Hamiltonian on n qubits expands as
   `E(x) = sum_S c_S prod_{q in S} (-1)^{x_q}`; `|S|` is the Pauli weight, and the fraction
   of the variance at each weight says whether the objective is a LOW-DEGREE (QAOA/VQE
   friendly) or a FULL-REGISTER function.  The programme's own torsion-space locality
   theorem says `d_ij` depends on exactly the `j-i-1` residues BETWEEN i and j, so the
   distogram term for a distant pair is a function of many residues and must carry high
   Walsh weight.  This measures how much.

   THE LEDGER'S WARNING IS OBEYED: "a Pauli spectrum of an UNCONDITIONED energy measures its
   worst clash" -- the top-10 configurations carried 99.6% of raw AMBER's Walsh variance.  So
   the PRIMARY spectrum here is of the RANK-UNIFORMISED objective (a monotone conditioning
   that cannot change any ranking and therefore cannot change any argmin), and the raw and
   99th-percentile-winsorised spectra are printed beside it.

2. DEGREE TRUNCATION -- the consequential form of (1).  Zero every coefficient above weight
   `d`, invert, and ask what the truncated objective RETAINS: its rank correlation with the
   full objective, its in-band correlation, and where ITS certified argmin sits.  If a
   weight-2 truncation preserves the ranking, then a pairwise classical model is sufficient
   and no entangling circuit is needed; if it does not, the objective is genuinely
   full-register and a bounded-correlation ansatz cannot represent its Boltzmann law.

3. THE SAME SPECTRUM FOR THE ORACLE TRUTH, labelled ORACLE.  If the true RMSD is high-degree
   in this register, no low-depth ansatz can represent the good set, whatever the objective.

4. THE REPRESENTABILITY THEOREM AND ITS CONSEQUENCE FOR THE PARETO QUESTION (EXACT).  The
   `mps2f` family contains every computational basis state exactly (all RY angles at 0 or pi
   give a product basis state, and the CNOT chain maps basis states to basis states).  So the
   unbudgeted variational optimum of `CVaR_alpha` for ANY alpha is a DELTA on the objective's
   certified argmin -- diversity exactly zero.  Every VQE point with `D > 0` is therefore a
   PARTIALLY CONVERGED INTERMEDIATE, not a designed frontier point, and its position on the
   (member error, diversity) plane is a coordinate on a collapse trajectory whose endpoint is
   fixed by the objective.  This is verified numerically, not asserted.

RUN:  python -m s17.q_theory run
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from core import quantum as Q
from s14 import vqe_lib as V
from s14 import vqe_run as R
from s16 import qphase_lib as QP
from s17 import q_lib as L


# ------------------------------------------------------------------ Walsh transform
def fwht(a):
    """In-place fast Walsh-Hadamard transform, normalised so `c = FWHT(f)/2^n` are the
    Fourier coefficients of `f = sum_S c_S chi_S` with `chi_S(x) = prod_{q in S} (-1)^{x_q}`.
    """
    a = np.array(a, dtype=np.float64, copy=True)
    n = a.size
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(-1)
        h *= 2
    return a / n


def ifwht(c):
    """Inverse: `f = sum_S c_S chi_S`. The unnormalised FWHT of `c`."""
    a = np.array(c, dtype=np.float64, copy=True)
    n = a.size
    h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(-1)
        h *= 2
    return a


def degree_spectrum(f, nq):
    """Fraction of the (non-constant) Walsh variance at each Pauli weight."""
    c = fwht(f)
    w = np.bincount(np.arange(1 << nq, dtype=np.int64) >> 0)  # placeholder, replaced below
    deg = _popcount(np.arange(1 << nq, dtype=np.int64))
    p = c ** 2
    tot = float(p.sum() - p[0])
    out = np.zeros(nq + 1)
    for d in range(1, nq + 1):
        out[d] = float(p[deg == d].sum())
    return out / max(tot, 1e-300), c, deg


def _popcount(a):
    a = np.asarray(a, np.int64)
    c = np.zeros_like(a)
    x = a.copy()
    while x.any():
        c += (x & 1)
        x >>= 1
    return c


def truncate(c, deg, d):
    cc = c.copy()
    cc[deg > d] = 0.0
    return ifwht(cc)


# ------------------------------------------------------------------ the experiment
def spectrum_cell(pdb, degrees=(1, 2, 3, 4, 6, 8)):
    ins = L.inst(pdb)
    nq = ins.n_qubits
    E = ins.hamil()
    out = {"pdb": pdb, "n": ins.n, "n_qubits": nq, "N": ins.N, "fold": ins.fold}
    fields = {
        # PRIMARY: monotone rank conditioning -- cannot change any ranking or argmin
        "hamil_uniformised": V.uniformise(E),
        "hamil_raw": E.astype(np.float64),
        "prior_uniformised": V.uniformise(ins.prior),
        "legacy_uniformised": V.uniformise(ins.legacy),
        "TRUTH_uniformised_ORACLE": ins.u_truth,
    }
    for name, f in fields.items():
        spec, c, deg = degree_spectrum(f, nq)
        row = {"var_frac_by_weight": [float(x) for x in spec],
               "cum_var_frac": [float(x) for x in np.cumsum(spec)],
               "mean_weight": float((spec * np.arange(nq + 1)).sum())}
        if name.startswith("hamil_uni") or name.startswith("TRUTH"):
            tr = {}
            ref = ins.rmsd
            for d in degrees:
                g = truncate(c, deg, d)
                tied = np.flatnonzero(g == g.min())
                order = np.argsort(f, kind="mergesort")[:max(64, ins.N // 100)]
                gorder = np.argsort(g, kind="mergesort")[:max(64, ins.N // 100)]
                tr[str(d)] = {
                    "rho_with_full": float(V.spearman(g, f)),
                    "rho_inband_with_full": float(V.spearman(g[order], f[order])),
                    "rho_with_rmsd_ORACLE": float(V.spearman(g, ref)),
                    "rho_ownband_with_rmsd_ORACLE": float(
                        V.spearman(g[gorder], ref[gorder])),
                    "ownband_mean_rmsd_ORACLE": float(ref[gorder].mean()),
                    "ownband_best_rmsd_ORACLE": float(ref[gorder].min()),
                    "argmin_rmsd_ORACLE": float(ref[tied].mean()),
                    "argmin_ties": int(tied.size),
                    "argmin_pct_ORACLE": float((ref < ref[tied].mean()).mean()),
                }
            row["truncation"] = tr
        out[name] = row
    # the objective's own certified argmin, for the truncation rows to be read against
    out["certified_argmin_rmsd_ORACLE"] = float(ins.rmsd[E == E.min()].mean())
    out["rmsd_best_in_space_ORACLE"] = float(ins.rmsd.min())
    out["rho_full_with_rmsd_ORACLE"] = float(V.spearman(E, ins.rmsd))
    ob = np.argsort(E, kind="mergesort")[:max(64, ins.N // 100)]
    out["full_ownband_mean_rmsd_ORACLE"] = float(ins.rmsd[ob].mean())
    out["full_ownband_best_rmsd_ORACLE"] = float(ins.rmsd[ob].min())
    return out


def weight2_range(pdb):
    """WHERE the weight-2 Walsh mass sits, as a function of qubit separation.

    This is the question an MPS chain ansatz actually cares about.  A depth-2 CNOT chain
    carries correlations over a bounded number of sites; if the objective's pairwise mass is
    SHORT-range the ansatz is matched to it, and if it is LONG-range no bounded-depth chain
    can represent the corresponding correlations and the entanglement in the circuit is
    decorative either way.  Reported for the rank-uniformised objective (the ledger's
    monotone conditioning) and for the ORACLE truth beside it.
    """
    ins = L.inst(pdb)
    nq = ins.n_qubits
    idx = np.arange(1 << nq, dtype=np.int64)
    deg = _popcount(idx)
    two = np.flatnonzero(deg == 2)
    # decode the two set bits of each weight-2 index; qubit 0 is the MSB
    b = idx[two]
    hi = nq - 1 - np.floor(np.log2(b)).astype(int)
    lo = nq - 1 - np.floor(np.log2(b - (1 << (nq - 1 - hi)))).astype(int)
    sep = np.abs(hi - lo)
    out = {"pdb": pdb, "n_qubits": nq}
    for name, f in (("hamil_uniformised", V.uniformise(ins.hamil())),
                    ("TRUTH_uniformised_ORACLE", ins.u_truth)):
        c = fwht(f)
        p2 = c[two] ** 2
        tot = float(p2.sum())
        by = np.zeros(nq)
        for s in range(1, nq):
            by[s] = float(p2[sep == s].sum())
        out[name] = {"weight2_mass_by_qubit_separation": (by / max(tot, 1e-300)).tolist(),
                     "mean_separation": float((by / max(tot, 1e-300) *
                                               np.arange(nq)).sum()),
                     "frac_within_2_qubits": float(by[:3].sum() / max(tot, 1e-300)),
                     "frac_within_1_residue": float(by[:2].sum() / max(tot, 1e-300))}
    return out


def budget_floor(pdb, budgets=(20, 36, 64, 128, 256, 512, 1024, 2048, 8192),
                 seeds=range(8)):
    """HOW CHEAP IS THIS SEARCH PROBLEM, REALLY?  (sprint section 21, taken to its end.)

    Sprint 16 compared CVaR-VQE against classical searches at 8,192 objective evaluations
    and found the VQE loses at every rung.  It never asked the prior question: how many
    evaluations does the classical arm actually NEED?  If the answer is tens, then 8,192 was
    never a budget, it was a formality, and every matched-budget comparison in this line of
    work has been comparing two arms on the far side of saturation.

    The Walsh spectrum predicts the answer: the objective is a separable per-residue field to
    93% of its variance, and a separable objective's exact optimum is found by ONE pass of
    coordinate descent.  This measures the fraction of seeds reaching the CERTIFIED global
    optimum as a function of budget, for greedy 1-opt and for a cold Metropolis chain.
    """
    ins = L.inst(pdb)
    E = ins.hamil()
    opt = float(E.min())
    out = {"pdb": pdb, "n": ins.n, "N": ins.N,
           "certified_argmin_rmsd_ORACLE": float(ins.rmsd[E == E.min()].mean()),
           "one_pass_cost_n_times_k": int(ins.n * (ins.k - 1) + 1)}
    from s15 import seed as SD
    for b in budgets:
        hits, rets, used = [], [], []
        for s in seeds:
            c = V.search_greedy(E, ins.n, ins.k, b,
                                SD.stable_rng(pdb, s, f"bg{b}", salt=L.SALT))
            hits.append(bool(c.best_e <= opt + 1e-12))
            rets.append(float(ins.rmsd[c.best_i]))
            used.append(int(c.used))
        hm, rm = [], []
        for s in seeds:
            c = L.search_metropolis(E, ins.n, ins.k, b, 0.05 * float(E.std()),
                                    SD.stable_rng(pdb, s, f"bm{b}", salt=L.SALT))
            hm.append(bool(c.best_e <= opt + 1e-12))
            rm.append(float(ins.rmsd[c.best_i]))
        out[str(b)] = {"greedy_frac_certified": float(np.mean(hits)),
                       "greedy_rmsd_ORACLE": float(np.mean(rets)),
                       "greedy_used": float(np.mean(used)),
                       "metro_frac_certified": float(np.mean(hm)),
                       "metro_rmsd_ORACLE": float(np.mean(rm))}
    return out


def representability(nq=18, seed=0):
    """(4) EXACT. The ansatz contains every computational basis state, so the unbudgeted
    variational optimum of CVaR at any alpha is a delta on the certified argmin (D = 0)."""
    an = R.make_ansatz("mps2f", nq)
    rng = np.random.default_rng(seed)
    rows = []
    P = an.n_params()
    nb = P // nq
    for _ in range(4):
        target = rng.integers(0, 2, nq)
        th = np.zeros(P).reshape(nb, nq)
        # every RY angle 0 except the FINAL layer, which flips the wires that must be 1;
        # the CNOT chain then maps this basis state to another basis state, so we read the
        # realised basis state off the distribution rather than predicting it.
        th[-1] = np.pi * target
        p = np.asarray(an.probs(th.reshape(-1)), float)
        rows.append({"max_prob": float(p.max()),
                     "entropy_bits": float(-(p[p > 0] * np.log2(p[p > 0])).sum()),
                     "argmax": int(np.argmax(p))})
    reach = set()
    for _ in range(64):
        target = rng.integers(0, 2, nq)
        th = np.zeros(P).reshape(nb, nq)
        th[-1] = np.pi * target
        p = np.asarray(an.probs(th.reshape(-1)), float)
        reach.add(int(np.argmax(p)))
    return {"delta_states": rows,
            "distinct_basis_states_reached_in_64_tries": len(reach),
            "min_max_prob": float(min(r["max_prob"] for r in rows)),
            "PASS_delta_representable": bool(
                min(r["max_prob"] for r in rows) > 1 - 1e-9)}


def gradient_scaling(pdbs=("1CS9", "1N9U"), seeds=(0, 1, 2, 3), shots=512, alpha=0.25):
    """Gradient magnitude and its sampling noise, at the two register sizes available.

    Two registers is not a scaling law and is not reported as one; what it IS good for is
    the SIGNAL-TO-NOISE of the estimator the deployable arm actually uses: the cosine between
    the 512-shot score-function gradient and the exact-expectation gradient at the same
    theta.  If that cosine is high, the optimiser is not noise-limited and any failure is
    representational, not statistical.
    """
    rows = []
    for pdb in pdbs:
        ins = L.inst(pdb)
        nq = ins.n_qubits
        E = ins.hamil()
        a = np.arange(1 << nq, dtype=np.int64)
        bits_all = ((a[:, None] >> np.arange(nq - 1, -1, -1)[None, :]) & 1).astype(np.int8)
        an = R.make_ansatz("mps2f", nq)
        for s in seeds:
            rng = np.random.default_rng(s)
            th = R.init_theta(an, rng, 0.8, "random")
            ge, _ = Q.cvar_gradient_exact(an, th, bits_all, E, alpha)
            b = an.sample(th, shots, rng)
            idx = R._bits_to_index(b)
            gs, _ = Q.cvar_gradient(an, th, b, E[idx], alpha, baseline="const")
            cos = float(ge @ gs / (np.linalg.norm(ge) * np.linalg.norm(gs) + 1e-30))
            rows.append({"pdb": pdb, "n_qubits": nq, "seed": s,
                         "norm_exact": float(np.linalg.norm(ge)),
                         "norm_sampled": float(np.linalg.norm(gs)),
                         "cos_sampled_vs_exact": cos,
                         "per_param_exact_sd": float(ge.std()),
                         "n_params": int(an.n_params())})
    return rows


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        done = L.ck_load("theory")
        if "representability" not in done:
            L.ck("theory", "representability", representability())
        tg = sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS19
        for pdb in tg:
            if f"spec_{pdb}" in done:
                print("  skip", pdb, flush=True)
                continue
            L.gate(pdb)
            t0 = time.time()
            r = spectrum_cell(pdb)
            L.ck("theory", f"spec_{pdb}", r)
            done[f"spec_{pdb}"] = r
            s = r["hamil_uniformised"]["var_frac_by_weight"]
            print(f"  {pdb} {time.time()-t0:.0f}s  mean Pauli weight "
                  f"{r['hamil_uniformised']['mean_weight']:.2f}  "
                  f"w<=2 {sum(s[:3]):.3f}  w<=4 {sum(s[:5]):.3f}", flush=True)
    elif mode == "range":
        done = L.ck_load("theory")
        tg = sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS9
        for pdb in tg:
            if f"range_{pdb}" in done:
                continue
            r = weight2_range(pdb)
            L.ck("theory", f"range_{pdb}", r)
            print(f"  {pdb} mean sep {r['hamil_uniformised']['mean_separation']:.2f} "
                  f"within2 {r['hamil_uniformised']['frac_within_2_qubits']:.3f}  "
                  f"TRUTH mean sep {r['TRUTH_uniformised_ORACLE']['mean_separation']:.2f} "
                  f"within2 {r['TRUTH_uniformised_ORACLE']['frac_within_2_qubits']:.3f}",
                  flush=True)
    elif mode == "budget":
        done = L.ck_load("theory")
        for pdb in (sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS19):
            if f"budget_{pdb}" in done:
                continue
            r = budget_floor(pdb)
            L.ck("theory", f"budget_{pdb}", r)
            print(f"  {pdb} one-pass cost {r['one_pass_cost_n_times_k']}  "
                  + "  ".join(f"b{b}:{r[b]['greedy_frac_certified']:.2f}"
                              for b in ("36", "64", "128", "512", "8192")), flush=True)
    elif mode == "grad":
        L.ck("theory", "gradient", gradient_scaling())
        print(json.dumps(L.ck_load("theory")["gradient"], indent=1))
    else:
        raise SystemExit(mode)


if __name__ == "__main__":
    main()
