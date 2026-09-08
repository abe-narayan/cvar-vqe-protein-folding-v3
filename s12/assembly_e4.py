"""E4 -- QUBO / Ising form of the deployable assembly problem, one instance per target.

Variables: x_{p,i} in {0,1}, one-hot over the k candidate pieces at each of m positions
(pieces = top-k by the E2 `combo` key, as in E3).  Energy (lower = better):

    E(x) = sum_p sum_i h[p,i] x[p,i] + sum_{p<q} sum_{i,j} J[(p,i),(q,j)] x[p,i] x[q,j]
           + A * sum_p (sum_i x[p,i] - 1)^2

    h[p,i]         = sum of shipped-distogram Bayes risk over residue pairs INSIDE piece i at p
                     (exact: they depend only on the piece's own torsions)
    J[(p,i),(q,j)] = sum of risk over pairs (a in p, b in q) of the two-piece sub-assembly
                     + W_J * (-log P(state_first(q,j) | state_last(p,i)))   [adjacent q = p+1]
                     For NON-adjacent (p,q) (m=3 only) the cross distances depend on the
                     middle piece: J is the MEAN over the middle candidates (mean-field), so
                     the m=3 QUBO is an approximation of the true objective; the m=2 QUBO is
                     exact (E_true = E_qubo up to the constant).
    A              = 2 * (max|h| + max row-sum |J|), so every infeasible x costs more than any
                     feasible one.
Stored as Q (N x N, upper triangular, x^T Q x + const = E) with N = m*k, index p*k + i.
Also stored: the classical optimum by exhaustive enumeration over the k^m feasible assignments,
the TRUE-objective optimum (full build of every assignment), and (ORACLE, diagnostic) the
CA-RMSD to the native of both optima.  Ising form: x = (1 - s)/2.
"""
from __future__ import annotations
import os, sys, json, time, itertools
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_bank as AB
from s12 import assembly_common as AC
from s12 import assembly_e2 as E2
from s12 import assembly_e3 as E3

W_J = 1.0


def risk_of_pairs(dg, D, pair_rows):
    grid = dg["grid"]; risk = dg["risk"][pair_rows]
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].sum(1)


def instance(t, K, comp, k, u, dg):
    n, fold = t["n"], t["fold"]; iv = AC.intervals(comp); m = len(iv); N = m * k
    R = AB.residues(fold); i, j = I.pair_index(n)
    tops, tors, ids = [], [], []
    for q, (s, L) in enumerate(iv):
        top = E3.top_pieces(t, K, s, L, q == 0, q == m - 1, k)
        ph, ps, _ = AB.piece_torsions(fold, L, top)
        tops.append(top); tors.append((ph, ps)); ids.append(AB.pieces(fold, L)["start"][top])
    # unary: intra-piece risk
    h = np.zeros((m, k))
    for q, (s, L) in enumerate(iv):
        rows = np.where((i >= s) & (j < s + L))[0]
        if len(rows) == 0:
            continue
        ca = AC.build_many(tors[q][0], tors[q][1])
        D = I.pair_dists(ca, i[rows] - s, j[rows] - s)
        h[q] = risk_of_pairs(dg, D, rows)
    # pairwise
    J = np.zeros((m, k, m, k))
    T = R["trans"]
    for p in range(m):
        for q in range(p + 1, m):
            sp, Lp = iv[p]; sq, Lq = iv[q]
            rows = np.where((i >= sp) & (i < sp + Lp) & (j >= sq) & (j < sq + Lq))[0]
            if q == p + 1:
                # two-piece sub-assembly: exact
                PHI = np.zeros((k * k, sq + Lq - sp)); PSI = np.zeros_like(PHI)
                for a in range(k):
                    for b in range(k):
                        PHI[a * k + b, :Lp] = tors[p][0][a]; PSI[a * k + b, :Lp] = tors[p][1][a]
                        PHI[a * k + b, Lp:] = tors[q][0][b]; PSI[a * k + b, Lp:] = tors[q][1][b]
                ca = AC.build_many(PHI, PSI)
                D = I.pair_dists(ca, i[rows] - sp, j[rows] - sp)
                J[p, :, q, :] = risk_of_pairs(dg, D, rows).reshape(k, k)
                sa = R["state"][ids[p] + Lp - 1]; sb = R["state"][ids[q]]
                J[p, :, q, :] += W_J * (-np.log(T[sa[:, None], sb[None, :]]))
            else:
                # mean-field over the middle pieces (m == 3 only: p=0, q=2, middle=1)
                mid = p + 1; sm, Lm = iv[mid]
                acc = np.zeros((k, k))
                for c in range(k):
                    PHI = np.zeros((k * k, n)); PSI = np.zeros_like(PHI)
                    for a in range(k):
                        for b in range(k):
                            PHI[a * k + b, sp:sp + Lp] = tors[p][0][a]; PSI[a * k + b, sp:sp + Lp] = tors[p][1][a]
                            PHI[a * k + b, sm:sm + Lm] = tors[mid][0][c]; PSI[a * k + b, sm:sm + Lm] = tors[mid][1][c]
                            PHI[a * k + b, sq:sq + Lq] = tors[q][0][b]; PSI[a * k + b, sq:sq + Lq] = tors[q][1][b]
                    ca = AC.build_many(PHI, PSI)
                    D = I.pair_dists(ca, i[rows], j[rows])
                    acc += risk_of_pairs(dg, D, rows).reshape(k, k)
                J[p, :, q, :] = acc / k
    Jf = J.reshape(N, N)
    A = 2.0 * (np.abs(h).max() + np.abs(Jf).sum(1).max())
    Q = np.zeros((N, N))
    Q[np.arange(N), np.arange(N)] = h.ravel() - A
    for p in range(m):
        for a in range(k):
            for b in range(a + 1, k):
                Q[p * k + a, p * k + b] = 2.0 * A
    Q += np.triu(Jf, 1)
    const = A * m
    # exhaustive over feasible assignments: QUBO energy and TRUE objective (full build)
    grids = np.array(list(itertools.product(range(k), repeat=m)))
    X = np.zeros((len(grids), N))
    for q in range(m):
        X[np.arange(len(grids)), q * k + grids[:, q]] = 1.0
    Eq = np.einsum("bi,ij,bj->b", X, Q, X) + const
    PHI = np.zeros((len(grids), n)); PSI = np.zeros_like(PHI)
    for q, (s, L) in enumerate(iv):
        PHI[:, s:s + L] = tors[q][0][grids[:, q]]; PSI[:, s:s + L] = tors[q][1][grids[:, q]]
    ca = AC.build_many(PHI, PSI)
    D = I.pair_dists(ca, i, j)
    Etrue = risk_of_pairs(dg, D, np.arange(len(i)))
    for q in range(m - 1):
        sa = R["state"][ids[q][grids[:, q]] + iv[q][1] - 1]; sb = R["state"][ids[q + 1][grids[:, q + 1]]]
        Etrue += W_J * (-np.log(T[sa, sb]))
    rr = I.kabsch_rmsd_batch(ca, u["nat_ca"])              # ORACLE diagnostic
    bq = int(np.argmin(Eq)); bt = int(np.argmin(Etrue))
    # a brute-force check over ALL 2^N states when small, to confirm the penalty
    brute = None
    if N <= 16:
        S = ((np.arange(2 ** N)[:, None] >> np.arange(N)[None, :]) & 1).astype(float)
        Eall = np.einsum("bi,ij,bj->b", S, Q, S) + const
        ba = int(np.argmin(Eall)); feasible = all(S[ba, p * k:(p + 1) * k].sum() == 1 for p in range(m))
        brute = {"min_energy": float(Eall[ba]), "feasible": bool(feasible), "matches_feasible_opt": bool(abs(Eall[ba] - Eq[bq]) < 1e-6)}
    out = {"Q": Q, "const": const, "A": A, "k": k, "m": m, "comp": np.array(comp), "piece_start": np.array(ids),
           "piece_bank_idx": np.array(tops), "n": n, "fold": fold, "W_J": W_J,
           "opt_qubo_choice": grids[bq], "opt_qubo_energy": float(Eq[bq]), "opt_qubo_true": float(Etrue[bq]), "opt_qubo_rr": float(rr[bq]),
           "opt_true_choice": grids[bt], "opt_true_energy": float(Etrue[bt]), "opt_true_qubo": float(Eq[bq]), "opt_true_rr": float(rr[bt]),
           "best_rr_in_instance": float(rr.min()), "corr_qubo_true": float(np.corrcoef(Eq, Etrue)[0, 1]) if m > 2 else 1.0,
           "energies_qubo": Eq.astype(np.float32), "energies_true": Etrue.astype(np.float32), "rr": rr.astype(np.float32)}
    return out, brute


def run_target(t):
    pdb, n = t["pdb"], t["n"]
    u = I.load_univ(pdb); dg = I.distogram(pdb); K = E2.target_keys(t)
    summ = {"pdb": pdb, "n": n, "fold": t["fold"], "instances": {}}
    for m in (2, 3):
        comp = AC.equal_split(n, m)
        for k in (4, 8):
            inst, brute = instance(t, K, comp, k, u, dg)
            sfx = "" if E3.KEY == "combo" else f"_{E3.KEY}"
            path = os.path.join(I.CACHE, f"qubo_{pdb}_m{m}_k{k}{sfx}.npz")
            np.savez_compressed(path, **inst)
            summ["instances"][f"m{m}_k{k}"] = {"path": os.path.relpath(path, ROOT), "N": m * k, "comp": list(comp), "A": inst["A"],
                                               "opt_qubo_energy": inst["opt_qubo_energy"], "opt_qubo_rr": inst["opt_qubo_rr"],
                                               "opt_true_energy": inst["opt_true_energy"], "opt_true_rr": inst["opt_true_rr"],
                                               "qubo_opt_is_true_opt": bool(np.array_equal(inst["opt_qubo_choice"], inst["opt_true_choice"])),
                                               "corr_qubo_true": inst["corr_qubo_true"], "best_rr_in_instance": inst["best_rr_in_instance"], "brute": brute}
    return summ


if __name__ == "__main__":
    tg = I.targets(); rows = []
    t0 = time.time()
    for q, t in enumerate(tg):
        rows.append(run_target(t))
        if q % 10 == 0:
            print(f"[{q+1}/126] {t['pdb']} {time.time()-t0:.0f}s", flush=True)
    pdbs = [r["pdb"] for r in rows]
    agg = {"n_targets": len(rows), "W_J": W_J, "shortlist_key": E3.KEY, "instances": {}}
    for key in rows[0]["instances"]:
        d = [r["instances"][key] for r in rows]
        agg["instances"][key] = {"N": d[0]["N"], "opt_qubo_rr": AC.group_means([x["opt_qubo_rr"] for x in d], pdbs),
                                 "opt_true_rr": AC.group_means([x["opt_true_rr"] for x in d], pdbs),
                                 "best_rr_in_instance": AC.group_means([x["best_rr_in_instance"] for x in d], pdbs),
                                 "qubo_opt_is_true_opt_frac": float(np.mean([x["qubo_opt_is_true_opt"] for x in d])),
                                 "corr_qubo_true_mean": float(np.mean([x["corr_qubo_true"] for x in d])),
                                 "brute_penalty_ok_frac": float(np.mean([x["brute"]["feasible"] and x["brute"]["matches_feasible_opt"] for x in d if x["brute"]])) if any(x["brute"] for x in d) else None}
    agg["per_target"] = rows
    print(I.write("assembly_e4_qubo" if E3.KEY == "combo" else f"assembly_e4_qubo_{E3.KEY}", agg))
    for key, v in agg["instances"].items():
        print(key, {kk: (round(vv["all"], 3) if isinstance(vv, dict) else vv) for kk, vv in v.items()})
