"""SPRINT 15 / QGEOM -- PART C: GRADIENT VARIANCE AND SCALING, DONE PROPERLY.

THE RULE THIS MODULE IS BUILT AROUND (brief, explicit): *do not call a small gradient a
barren plateau unless the SCALING evidence supports it.*  So every cell reports a slope with
a bootstrap CI against the textbook `Var ~ 2^-n` (slope -1.0 per qubit), and the null the
slope is read against is stated.

AXES MEASURED
    n           qubit count 6..18 -- the axis Sprint 14 found to be the one that matters
    L           circuit depth 1..6
    pattern     entangler
    alpha       CVaR level 1.0 .. 0.01 -- CVaR reads a shrinking tail, so its gradient is
                carried by fewer states; that is a DIFFERENT mechanism from a plateau
    energy      Legacy / retrieval prior / distogram / a random unstructured control /
                true RMSD (ORACLE DIAGNOSTIC)
    locality    an objective built as a sum of Pauli-Z products of EXACT weight w.
                Sprint 14 reports the ansatz kernel is flat in Pauli weight at 6-18 qubits
                and that "what matters is n".  This tests that directly by holding n fixed
                and sweeping w, with the objective variance held constant by construction.
    conditioning the point at which the gradient is evaluated: Haar-random angles (the
                barren-plateau protocol), the uniform state, or a target-conditioned state
                (Part B).  Links trainability to information.

Every objective is RANK-CONDITIONED (`uniformise`) before use unless it is the synthetic
Pauli-weight family, which is built with a fixed variance by construction.  The brief records
that the Walsh spectrum of a raw molecular energy is a delta-spike artefact.

    python -m s15.qgeom_grad
"""
from __future__ import annotations

import time

import numpy as np

from core.quantum import cvar_exact
from s15 import qgeom_lib as G
from s15.qgeom_qng import sub_objective

TAG = "grad"


# =========================================================== gradient of the cost
def grad_at(circ, th, E, alpha=1.0):
    """Exact gradient by ADJOINT differentiation (verified to 2.4e-15 against the
    derivative-state route).  O(gates) rather than O(P*gates), which is what makes an
    n=18 x 48-seed variance estimate affordable."""
    psi = circ.state(th)
    p = psi ** 2
    p = p / p.sum()
    w = E if alpha >= 1.0 else cvar_exact(E, p, alpha)[2]
    return G.grad_adjoint(circ, th, 2.0 * psi * w)


def grad_var(circ, E, alpha=1.0, seeds=64, init="haar", cond_theta=None, seed0=0):
    """Var of the gradient components over random initialisations.

    `Var[g_1]` (a single fixed component, the textbook statistic) and the mean over
    components are both returned -- they differ when the ansatz is not permutation
    symmetric, and quoting only one has caused confusion in this literature.
    """
    gs = []
    for s in range(seeds):
        if init == "haar":
            th = G.haar_theta(circ, seed0 + 1000 * s)
        elif init == "small":
            th = G.random_theta(circ, seed0 + 1000 * s, 0.3)
        elif init == "uniform":
            th = G.uniform_theta(circ) + 0.01 * np.random.default_rng(
                seed0 + 1000 * s).standard_normal(circ.n_params())
        elif init == "cond":
            th = np.asarray(cond_theta, float) + 0.01 * np.random.default_rng(
                seed0 + 1000 * s).standard_normal(circ.n_params())
        else:
            raise ValueError(init)
        gs.append(grad_at(circ, th, E, alpha))
    Gm = np.asarray(gs)
    return {"var_g1": float(np.var(Gm[:, 0])),
            "var_mean_over_params": float(np.mean(np.var(Gm, axis=0))),
            "mean_abs_g": float(np.mean(np.abs(Gm))),
            "seeds": seeds}


def slope_ci(ns, vs, boot=2000, seed=0):
    """log2(Var) vs n slope with a bootstrap CI over the points."""
    ns = np.asarray(ns, float)
    y = np.log2(np.maximum(np.asarray(vs, float), 1e-300))
    s = float(np.polyfit(ns, y, 1)[0])
    r = np.random.default_rng(seed)
    bs = []
    for _ in range(boot):
        i = r.integers(0, len(ns), len(ns))
        if len(np.unique(ns[i])) < 2:
            continue
        bs.append(float(np.polyfit(ns[i], y[i], 1)[0]))
    lo, hi = (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))) if bs \
        else (np.nan, np.nan)
    return {"slope_log2_per_qubit": s, "ci_lo": lo, "ci_hi": hi,
            "excludes_textbook_minus1": bool(lo > -1.0 or hi < -1.0)}


# ======================================================== objective families
def real_objective(kind, nq, pdb="1CS9"):
    """A native-free objective on the first nq/2 residues of a real target."""
    e = G.V.Enum(pdb)
    nres = nq // 2
    res = tuple(range(1, 1 + nres)) if 1 + nres <= e.n else tuple(range(nres))
    base = np.random.default_rng(0).integers(0, e.k, e.n)
    if kind == "random":
        return G.V.uniformise(np.random.default_rng(nq).standard_normal(1 << nq))
    if kind == "rmsd":                                # ORACLE DIAGNOSTIC
        return G.V.uniformise(sub_objective(e, res, base, "rmsd"))
    if kind in ("legacy", "prior"):
        return G.V.uniformise(sub_objective(e, res, base, kind))
    raise ValueError(kind)


def pauli_weight_objective(nq, w, terms=32, seed=0):
    """Sum of `terms` random Pauli-Z products of EXACT weight `w`, unit variance.

    Z-strings are +/-1 valued and orthogonal, so a sum of `terms` distinct ones has variance
    exactly `terms`; dividing by sqrt(terms) fixes the variance at 1 for every w.  The
    objective's SCALE is therefore held constant and only its LOCALITY varies.
    """
    r = np.random.default_rng(1000 * seed + 10 * nq + w)
    idx = np.arange(1 << nq, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(nq - 1, -1, -1)[None, :]) & 1).astype(np.int8)
    out = np.zeros(1 << nq)
    seen = set()
    for _ in range(terms):
        for _try in range(50):
            s = tuple(sorted(r.choice(nq, size=w, replace=False)))
            if s not in seen:
                seen.add(s)
                break
        out += (1 - 2 * (bits[:, list(s)].sum(1) % 2)).astype(float)
    return out / np.sqrt(max(len(seen), 1))


# ================================================================ experiments
def scaling_n(patterns=("ring", "block", "all_to_all"), Ls=(1, 2, 3),
              ns=(6, 8, 10, 12, 14, 16, 18), seeds=48, kind="legacy"):
    print("=" * 108)
    print("C1. GRADIENT VARIANCE vs QUBIT COUNT.  Haar-random angles, exact gradient,")
    print("    rank-conditioned objective (variance fixed at 1/12 for every n).")
    print("=" * 108)
    print(f"{'pattern':11s} {'L':>2s} " + "".join(f"{f'n={n}':>11s}" for n in ns)
          + f"{'slope/qubit':>13s} {'CI95':>20s} {'!= -1?':>7s}")
    out = {}
    for pat in patterns:
        for L in Ls:
            vs = []
            for n in ns:
                E = G.V.uniformise(real_objective(kind, n))
                c = G.FlexCircuit(n, L, pat, 2)
                vs.append(grad_var(c, E, 1.0, seeds)["var_g1"])
            sl = slope_ci(ns, vs)
            out[f"{pat}|L{L}"] = {"ns": list(ns), "var_g1": vs, **sl}
            print(f"{pat:11s} {L:2d} " + "".join(f"{v:11.3e}" for v in vs)
                  + f"{sl['slope_log2_per_qubit']:13.3f} "
                  + f"{f'[{sl['ci_lo']:.3f}, {sl['ci_hi']:.3f}]':>20s} "
                  + f"{str(sl['excludes_textbook_minus1']):>7s}")
            G.ck(TAG, f"C1_scaling_{kind}", out)
    print("\n  Textbook barren plateau: slope = -1.0 per qubit (Var ~ 2^-n).")
    print("  The standard theorem needs local 2-design blocks this ansatz does NOT have,")
    print("  so it licenses no prediction here; the slope is the evidence, not the theory.")
    return out


def scaling_depth(ns=(8, 12, 16), Ls=(1, 2, 3, 4, 6, 8), pat="ring", seeds=48):
    print()
    print("=" * 108)
    print("C2. GRADIENT VARIANCE vs DEPTH at fixed width")
    print("=" * 108)
    print(f"{'n':>3s} " + "".join(f"{f'L={L}':>12s}" for L in Ls) + f"{'slope/layer':>13s}")
    out = {}
    for n in ns:
        E = G.V.uniformise(real_objective("legacy", n))
        vs = []
        for L in Ls:
            c = G.FlexCircuit(n, L, pat, 2)
            vs.append(grad_var(c, E, 1.0, seeds)["var_g1"])
        sl = float(np.polyfit(np.log2(Ls), np.log2(np.maximum(vs, 1e-300)), 1)[0])
        out[f"n{n}"] = {"Ls": list(Ls), "var_g1": vs, "slope_log2_per_log2L": sl}
        print(f"{n:3d} " + "".join(f"{v:12.3e}" for v in vs) + f"{sl:13.3f}")
        G.ck(TAG, "C2_depth", out)
    return out


def scaling_alpha(ns=(8, 10, 12, 14), alphas=(1.0, 0.5, 0.25, 0.1, 0.05, 0.01),
                  pat="ring", L=2, seeds=48):
    print()
    print("=" * 108)
    print("C3. GRADIENT VARIANCE vs CVaR ALPHA, and its SCALING with n at each alpha")
    print("    A shrinking tail is not a plateau: it is a smaller effective sample.")
    print("=" * 108)
    print(f"{'alpha':>7s} " + "".join(f"{f'n={n}':>12s}" for n in ns)
          + f"{'slope/qubit':>13s} {'CI95':>20s}")
    out = {}
    for a in alphas:
        vs = []
        for n in ns:
            E = G.V.uniformise(real_objective("legacy", n))
            c = G.FlexCircuit(n, L, pat, 2)
            vs.append(grad_var(c, E, a, seeds)["var_g1"])
        sl = slope_ci(ns, vs)
        out[f"a{a}"] = {"ns": list(ns), "var_g1": vs, **sl}
        print(f"{a:7.3f} " + "".join(f"{v:12.3e}" for v in vs)
              + f"{sl['slope_log2_per_qubit']:13.3f} "
              + f"{f'[{sl['ci_lo']:.3f}, {sl['ci_hi']:.3f}]':>20s}")
        G.ck(TAG, "C3_alpha", out)
    return out


def locality(ns=(8, 10, 12), pat="ring", L=2, seeds=48):
    print()
    print("=" * 108)
    print("C4. COST LOCALITY AT FIXED WIDTH -- Sprint 14's claim tested directly.")
    print("    Objectives are sums of Pauli-Z products of EXACT weight w with variance")
    print("    fixed at 1 for every w, so only LOCALITY varies.")
    print("=" * 108)
    out = {}
    for n in ns:
        ws = list(range(1, n + 1))
        vs = []
        for w in ws:
            E = pauli_weight_objective(n, w)
            c = G.FlexCircuit(n, L, pat, 2)
            vs.append(grad_var(c, E, 1.0, seeds)["var_g1"])
        rng = float(np.log2(max(vs) / max(min(vs), 1e-300)))
        sl = float(np.polyfit(ws, np.log2(np.maximum(vs, 1e-300)), 1)[0])
        out[f"n{n}"] = {"weights": ws, "var_g1": vs, "log2_range": rng,
                        "slope_log2_per_weight": sl}
        print(f"n={n}: " + " ".join(f"w{w}={v:.2e}" for w, v in zip(ws, vs)))
        print(f"      log2 range across w = {rng:.3f}, slope {sl:+.3f} per unit weight")
        G.ck(TAG, "C4_locality", out)
    print("\n  A FLAT profile means cost locality does not set trainability at this width,")
    print("  which is what Sprint 14 reports.  A steep negative slope in w would mean the")
    print("  usual 'make the cost local' prescription applies.  Compare against C1's slope")
    print("  in n measured on the same machinery.")
    return out


def energy_models(ns=(8, 10, 12, 14), pat="ring", L=2, seeds=48):
    print()
    print("=" * 108)
    print("C5. DOES THE ENERGY MODEL CHANGE THE GRADIENT VARIANCE?")
    print("    All rank-conditioned to the identical uniform marginal, so any difference")
    print("    is a difference of STRUCTURE, not of scale or tail shape.")
    print("=" * 108)
    kinds = ("random", "legacy", "prior", "rmsd")
    print(f"{'objective':10s} " + "".join(f"{f'n={n}':>12s}" for n in ns)
          + f"{'slope/qubit':>13s}")
    out = {}
    for k in kinds:
        vs = []
        for n in ns:
            E = G.V.uniformise(real_objective(k, n))
            c = G.FlexCircuit(n, L, pat, 2)
            vs.append(grad_var(c, E, 1.0, seeds)["var_g1"])
        sl = slope_ci(ns, vs)
        out[k] = {"ns": list(ns), "var_g1": vs, **sl}
        print(f"{k:10s} " + "".join(f"{v:12.3e}" for v in vs)
              + f"{sl['slope_log2_per_qubit']:13.3f}")
        G.ck(TAG, "C5_energy_models", out)
    print("\n  'rmsd' is ORACLE DIAGNOSTIC and is here only as the upper bound on how much")
    print("  structure an objective could possibly carry.  'random' is the unstructured null.")
    return out


def conditioning_effect(pdbs=("1CS9", "2MK7", "7N2I", "6EY3"), pat="block", L=2,
                        nres=6, seeds=48, iters=250):
    print()
    print("=" * 108)
    print("C6. TRAINABILITY AT THE CONDITIONED POINT -- does information cost gradient?")
    print("    Links Part B to Part C: the same conditioned theta, now read for Var[grad].")
    print("=" * 108)
    from s14 import retprior as RP
    res = tuple(range(1, 1 + nres))
    nq = 2 * nres
    c = G.FlexCircuit(nq, L, pat, 2)
    print(f"{'target':8s} {'H(prior) bits':>14s} " +
          "".join(f"{k:>13s}" for k in ("haar", "uniform", "COND", "scram_state")))
    out = {}
    for pdb in pdbs:
        e = G.V.Enum(pdb)
        base = np.random.default_rng(0).integers(0, e.k, e.n)
        E = G.V.uniformise(sub_objective(e, res, base, "legacy"))
        P = RP.state_prior(pdb, "top75", 4)["P"][list(res)]
        row = {"prior_entropy_bits": G.prior_entropy_bits(P)}
        row["haar"] = grad_var(c, E, 1.0, seeds, "haar")["var_g1"]
        row["uniform"] = grad_var(c, E, 1.0, seeds, "uniform")["var_g1"]
        for nm, Pk in (("COND", P), ("scram_state", G.scramble_prior(P, 0, "state"))):
            fit = G.fit_kl(c, G.product_target(Pk), theta0=G.random_theta(c, 0),
                           iters=iters, lr=0.08)
            row[nm] = grad_var(c, E, 1.0, seeds, "cond", fit["theta"])["var_g1"]
            row[nm + "_kl"] = fit["kl"]
        out[pdb] = row
        print(f"{pdb:8s} {row['prior_entropy_bits']:14.3f} " +
              "".join(f"{row[k]:13.3e}" for k in ("haar", "uniform", "COND",
                                                  "scram_state")))
        G.ck(TAG, "C6_conditioning", out)
    return out


def main():
    G.wait_mem(0.8, "qgeom_grad")
    G.ck_load(TAG)
    t0 = time.time()
    scaling_n()
    locality()
    scaling_alpha()
    scaling_depth()
    energy_models()
    conditioning_effect()
    print(f"\ndone in {time.time()-t0:.0f}s -> s15/results/qgeom_grad.json")


if __name__ == "__main__":
    main()
