#!/usr/bin/env python
"""s29/s29_T_reach.py -- LANE T, sections 5 and 6.

(5) THE REACHABLE SET. The dynamical Lie algebra of the deployed RY/CNOT(chain+ring) ansatz by
    exact nested commutators on Pauli strings (re-using `s26/q_dla.py`, whose closure is
    cross-checked against a dense SVD closure at n = 4, 5), for n = 3..6 and depth 1..4; and the
    MPS cut-rank of the states the ansatz can prepare, which is what actually bounds "a uniform
    superposition over an arbitrary subset".

(6) THE ENTROPY TERM. The exact simplex minimiser of F(p) = CVaR_alpha(E; p) - T H(p) on the
    deployed rank ladder, by the Rockafellar-Uryasev representation: for fixed t the minimiser is
    Gibbs in the CLIPPED energy, p*(x) ~ exp(-(E_x - t)_+/(alpha T)), so the optimum is uniform on
    the prefix {E <= t} with an exponential tail of scale alpha*T; t is then the alpha-quantile of
    p*. This file solves that fixed point and reports the realised support, entropy and F, so the
    DEPLOYED state's m (74 to 79) can be read against the objective's own optimum.

PROPERTY MEASUREMENT. No native, no RMSD, no pool: sections 5 and 6 are properties of the ansatz
and of the objective, not of the data.

    python s29/s29_T_reach.py --dla --entropy
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s29_T_reach.json")


# ------------------------------------------------------------------------------ (5) the DLA
def dla_sweep(ns=(3, 4, 5, 6), depths=(1, 2, 3, 4)):
    from s26.q_dla import fixed_generators, lie_closure, numeric_closure, dim_so, dim_su
    out = {}
    for n in ns:
        for L in depths:
            gens = fixed_generators(n, L, ring=True, convention="right")
            cl = lie_closure(gens)
            row = dict(n=n, L=L, dim=int(cl["dim"]), dim_so=int(dim_so(n)),
                       dim_su_half=int(dim_su(1 << (n - 1))),
                       two_dim_su_quarter=int(2 * dim_su(1 << (n - 2))) if n >= 3 else None,
                       frac_of_so=float(cl["dim"] / dim_so(n)), n_gens=len(gens))
            if n <= 5:
                row["dim_numeric"] = int(numeric_closure(gens)["dim"])
                row["numeric_agrees"] = bool(row["dim_numeric"] == row["dim"])
            out[f"n{n}_L{L}"] = row
            print(f"n={n} L={L} dim={row['dim']} / so={row['dim_so']}"
                  + (f"  numeric={row.get('dim_numeric')}" if "dim_numeric" in row else ""),
                  flush=True)
    return out


def mps_cut_rank_of_sets(n=9, trials=200, sizes=(8, 16, 32, 75, 128), seed=0):
    """The Schmidt rank, across every cut, of the UNIFORM SUPERPOSITION over a subset S.

    The amplitude vector of |S> = 1/sqrt|S| sum_{x in S} |x> reshaped at cut k is the 0/1 indicator
    matrix of S; its rank is the number of DISTINCT suffix-sets among the prefixes present. A
    depth-L RY/CNOT chain prepares states of bond dimension <= 2^L (<= 2^(L+1) with the ring
    closure), so a set whose indicator has larger rank at some cut is NOT representable.
    Reports, per size: the max-over-cuts rank of a RANDOM subset, and of the PREFIX subset
    {0..m-1} (the classical top-m, which the deployed encoding hands the circuit).
    """
    rng = np.random.default_rng(seed)
    D = 1 << n
    out = {}
    for m in sizes:
        if m > D:
            continue
        ranks = []
        for _ in range(trials):
            S = rng.choice(D, size=m, replace=False)
            v = np.zeros(D)
            v[S] = 1.0
            r = max(np.linalg.matrix_rank(v.reshape(1 << k, 1 << (n - k)))
                    for k in range(1, n))
            ranks.append(int(r))
        v = np.zeros(D)
        v[:m] = 1.0
        rp = max(np.linalg.matrix_rank(v.reshape(1 << k, 1 << (n - k))) for k in range(1, n))
        out[str(m)] = dict(random_max_cut_rank_median=float(np.median(ranks)),
                           random_max_cut_rank_min=int(min(ranks)),
                           prefix_max_cut_rank=int(rp),
                           bond_dim_depth3_ring=1 << 4, trials=trials)
        print(f"|S|={m}: random max-cut rank median {np.median(ranks)}, "
              f"prefix {rp}, depth-3 ring bond dim <= 16", flush=True)
    return out


# --------------------------------------------------------------------------- (6) the entropy
def deployed_E(dim):
    r = np.arange(1, dim + 1, dtype=float)
    return (r - r.mean()) / r.std()


def simplex_optimum(E, alpha, T, n_t=20001):
    """Exact minimiser of F(p) = CVaR_alpha(E; p) - T H(p) over the simplex.

    The project's CVaR is the LOWER tail, whose Rockafellar-Uryasev form is a MAXIMUM:
        CVaR_alpha(E; p) = max_t [ t - (1/alpha) sum_x p_x (t - E_x)_+ ].
    The bracket is linear in p and concave in t and the simplex is compact convex, so Sion's
    minimax theorem lets the order be swapped:
        F* = min_p max_t [...] = max_t [ t + min_p ( sum_x p_x d_x(t) + T sum_x p_x log p_x ) ],
        d_x(t) = -(t - E_x)_+/alpha,
    and the inner minimum is the Gibbs law p_t(x) ~ exp( (t - E_x)_+ / (alpha T) ) with value
    -T log Z(t):
        F* = max_t [ t - T log sum_x exp( (t - E_x)_+ / (alpha T) ) ].
    The outer problem is one-dimensional and concave, so a fine grid is exact to the grid.
    A sanity check the code satisfies: at alpha = 1 this returns the Gibbs free energy
    -T log sum_x exp(-E_x/T) and p* is the Boltzmann law (S25 section 6.3).
    """
    E = np.asarray(E, float)
    ts = np.linspace(E.min() - 0.5, E.max() + 3.0, n_t)
    best = None
    for t in ts:
        c = np.maximum(t - E, 0.0) / alpha
        mx = c.max()
        Z_shift = np.exp((c - mx) / T).sum() if T > 0 else 1.0
        val = t - (mx + T * np.log(Z_shift)) if T > 0 else t - mx
        if best is None or val > best[0]:
            best = (val, t)
    val, t = best
    c = np.maximum(t - E, 0.0) / alpha
    w = np.exp((c - c.max()) / T)
    p = w / w.sum()
    H = float(-(p[p > 0] * np.log(p[p > 0])).sum())
    o = np.argsort(E, kind="stable")
    cum = np.cumsum(p[o])
    m_tail = int(np.searchsorted(cum, alpha, side="left") + 1)
    return dict(alpha=float(alpha), T=float(T), F=float(val), t_star=float(t),
                prefix_size=int((E <= t).sum()), entropy_nats=H,
                entropy_bits=float(H / np.log(2)),
                participation_ratio=float(1.0 / (p ** 2).sum()),
                realised_tail_m=m_tail,
                p_max=float(p.max()), p_min=float(p.min()),
                suppression_scale_in_ranks=float(alpha * T / (E[1] - E[0])))


def entropy_block():
    out = {}
    for n in (7, 9):
        D = 1 << n
        E = deployed_E(D)
        for (a, T) in ((0.18, 0.5), (0.25, 0.3), (1.0, 0.3), (0.18, 0.25), (0.18, 1.0),
                       (0.10, 0.5), (0.40, 0.5)):
            r = simplex_optimum(E, a, T)
            r["n"] = n
            r["dim"] = D
            out[f"n{n}_a{a}_T{T}"] = r
            print(f"n={n} alpha={a} T={T}: F={r['F']:.4f} prefix={r['prefix_size']} "
                  f"tail_m={r['realised_tail_m']} PR={r['participation_ratio']:.1f} "
                  f"H={r['entropy_bits']:.2f} bits", flush=True)
    return out



# ----------------------------------------------- (Q1) the free-energy gap on the REAL registers
def fgap_block(limit=12):
    """Per target: the exact simplex optimum of F at the deployed (alpha, T) on that target's OWN
    E (padding included), against the uniform state and against the trained circuit's F recorded
    in `s27/results/s28_B_rows.jsonl` (source vqe, J = 0). Native-free."""
    import json as _json
    from s25 import phys_lib as P
    from s27 import run_pool as RP
    from s22 import qcand_lib as QC
    from core import quantum as Q
    rows = {}
    for line in open(os.path.join(ROOT, "s27", "results", "s28_B_rows.jsonl"), encoding="utf-8"):
        r = _json.loads(line)
        if r["source"] == "vqe" and r["J"] == 0.0:
            rows.setdefault(r["pdb"], {})[r["seed"]] = r
    pdbs = P.targets()[::11][:limit]
    out = {}
    for pdb in pdbs:
        cand, ch, _ = RP.channels_for(pdb)
        E = QC.Encoding(RP.zr(ch["DIS"])).E
        D = len(E)
        r = simplex_optimum(E, 0.18, 0.5)
        pu = np.full(D, 1.0 / D)
        cvu, _, _ = Q.cvar_exact(E, pu, 0.18)
        Hu = float(-(pu * np.log(pu)).sum())
        Fu = float(cvu - 0.5 * Hu)
        ou = np.argsort(E, kind="stable")
        mu = int(np.searchsorted(np.cumsum(pu[ou]), 0.18, side="left") + 1)
        tr = rows.get(pdb, {})
        Fc = {s: float(v["F"]) for s, v in tr.items()}
        mc = {s: int(v["m"]) for s, v in tr.items()}
        Hc = {s: float(v["entropy_bits"]) for s, v in tr.items()}
        PRc = {s: float(v["pr"]) for s, v in tr.items()}
        gap = {s: (Fc[s] - Fu) / (r["F"] - Fu) for s in Fc}
        out[pdb] = dict(dim=D, F_opt=r["F"], F_uniform=Fu, F_circuit=Fc,
                        gap_closed=gap, m_opt=r["realised_tail_m"], m_uniform=mu, m_circuit=mc,
                        H_opt_bits=r["entropy_bits"], H_circuit_bits=Hc,
                        PR_opt=r["participation_ratio"], PR_circuit=PRc,
                        prefix_opt=r["prefix_size"], t_star=r["t_star"])
        print(f"{pdb}: F opt {r['F']:.4f} unif {Fu:.4f} circ {Fc} | m opt {r['realised_tail_m']} "
              f"unif {mu} circ {mc} | gap closed {[round(v,3) for v in gap.values()]}", flush=True)
    ks = list(out)
    agg = dict(
        n_targets=len(ks),
        F_opt_mean=float(np.mean([out[k]["F_opt"] for k in ks])),
        F_uniform_mean=float(np.mean([out[k]["F_uniform"] for k in ks])),
        F_circuit_mean=float(np.mean([np.mean(list(out[k]["F_circuit"].values())) for k in ks])),
        gap_closed_mean=float(np.mean([np.mean(list(out[k]["gap_closed"].values())) for k in ks])),
        m_opt_mean=float(np.mean([out[k]["m_opt"] for k in ks])),
        m_uniform_mean=float(np.mean([out[k]["m_uniform"] for k in ks])),
        m_circuit_mean=float(np.mean([np.mean(list(out[k]["m_circuit"].values())) for k in ks])),
        H_opt_bits_mean=float(np.mean([out[k]["H_opt_bits"] for k in ks])),
        H_circuit_bits_mean=float(np.mean([np.mean(list(out[k]["H_circuit_bits"].values())) for k in ks])),
        PR_opt_mean=float(np.mean([out[k]["PR_opt"] for k in ks])),
        PR_circuit_mean=float(np.mean([np.mean(list(out[k]["PR_circuit"].values())) for k in ks])),
        prefix_opt_mean=float(np.mean([out[k]["prefix_opt"] for k in ks])))
    print(_json.dumps(agg, indent=1))
    return dict(per_target=out, summary=agg)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dla", action="store_true")
    ap.add_argument("--sets", action="store_true")
    ap.add_argument("--entropy", action="store_true")
    ap.add_argument("--fgap", action="store_true")
    a = ap.parse_args(argv)
    out = {}
    if os.path.exists(OUT):
        out = json.load(open(OUT, encoding="utf-8"))
    out.setdefault("kind", "s29_T_reach")
    out.setdefault("lane", "T")
    out.setdefault("sprint", 29)
    if a.dla:
        out["dla"] = dla_sweep()
    if a.sets:
        out["subset_cut_rank"] = mps_cut_rank_of_sets()
    if a.entropy:
        out["entropy"] = entropy_block()
    if a.fgap:
        out["fgap"] = fgap_block()
    from s24 import stats_lib as ST
    ST.save_atomic(OUT, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
