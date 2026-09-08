"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 6 -- pricing the encoding's neighbourhood.

Section 1 measured that a single state flip moves an in-support CA-CA distance by a median
1.75-1.95 A, i.e. the landscape is NOT smooth in the encoded variables.  This file prices
that directly, on the complete enumerated 4^9 spaces, by comparing the two neighbourhoods
the encoding choice actually determines:

    RESIDUE move  ("one-hot"):  change one residue to ANY other state.  9 x 3 = 27 neighbours.
                  This is the move a one-hot register makes with two qubit flips, and it is
                  also the move a classical residue-wise local search makes.
    QUBIT flip    ("binary"):   flip ONE qubit of the binary index.  9 x 2 = 18 neighbours,
                  and each residue can only reach the 2 of its 3 alternatives whose index
                  differs in one bit.  `gray` reaches a different pair.

Measured on every objective, exactly, over the whole space:

    number of 1-move local minima, and the fraction of configurations that are one
    steepest-descent basin's worth away from the global optimum;
    hit rate of steepest descent from all 262,144 starts;
    the RMSD actually reached;
    random-walk autocorrelation of the objective (the standard ruggedness measure), from
    which the correlation length  l = -1/ln(rho_1)  follows.

Output: `s13/results/qarch_landscape.json`.

    python -m s13.qarch_landscape
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402


def neighbours(n, k, kind):
    """(k^n, m) array of neighbour indices under the given move set."""
    B = k ** n
    idx = np.arange(B)
    digits = np.stack(np.unravel_index(idx, (k,) * n))          # (n, B)
    out = []
    for i in range(n):
        place = k ** (n - 1 - i)
        cur = digits[i]
        if kind == "residue":
            alts = [a for a in range(1, k)]
            for a in alts:
                out.append(idx + ((cur + a) % k - cur) * place)
        elif kind in ("binary", "gray"):
            b = int(np.log2(k))
            code = cur if kind == "binary" else (cur ^ (cur >> 1))
            for bit in range(b):
                c2 = code ^ (1 << bit)
                new = c2 if kind == "binary" else _gray_inv(c2, b)
                out.append(idx + (new - cur) * place)
        else:
            raise ValueError(kind)
    return np.stack(out, axis=1).astype(np.int32)


def _gray_inv(g, b):
    g = np.asarray(g)
    s = g.copy()
    sh = 1
    while sh < b:
        s = s ^ (s >> sh)
        sh <<= 1
    return s


def descend(F, N):
    """Steepest descent from EVERY configuration at once. Returns the endpoint index."""
    cur = np.arange(len(F), dtype=np.int32)
    for _ in range(200):
        Fn = F[N[cur]]                                # (B, m)
        j = np.argmin(Fn, axis=1)
        best = N[cur, j]
        move = F[best] < F[cur] - 1e-12
        if not move.any():
            break
        cur = np.where(move, best, cur).astype(np.int32)
    return cur


def autocorr(F, N, steps=64, walkers=2000, seed=0):
    """Random-walk autocorrelation of F in the given neighbourhood."""
    rng = np.random.default_rng(seed)
    cur = rng.integers(0, len(F), walkers)
    traj = np.empty((steps + 1, walkers))
    traj[0] = F[cur]
    for t in range(steps):
        cur = N[cur, rng.integers(0, N.shape[1], walkers)]
        traj[t + 1] = F[cur]
    r = []
    m, s = traj.mean(), traj.std()
    for lag in range(1, 9):
        a = traj[:-lag].ravel(); b = traj[lag:].ravel()
        r.append(float(((a - m) * (b - m)).mean() / (s * s)))
    return r


def analyse(pdb_id):
    z = np.load(os.path.join(Q.RESULTS, f"qarch_enum_{pdb_id}.npz"))
    n, k = int(z["n"]), int(z["k"])
    rmsd = np.asarray(z["rmsd"], float)
    objs = {"legacy_total": np.asarray(z["legacy"], float),
            "prior_empirical": np.asarray(z["prior"], float),
            "ORACLE_rmsd": rmsd}
    rec = {"pdb": pdb_id, "n": n, "k": k, "configs": int(k ** n),
           "pool_mean_rmsd": float(rmsd.mean()), "best_rmsd": float(rmsd.min()),
           "moves": {}}
    for kind in ("residue", "binary", "gray"):
        N = neighbours(n, k, kind)
        rec["moves"][kind] = {"n_neighbours": int(N.shape[1]), "objectives": {}}
        for nm, F in objs.items():
            end = descend(F, N)
            Fn = F[N]
            is_min = F <= Fn.min(1) + 1e-12
            gstar = int(np.argmin(F))
            hit = float((F[end] <= F[gstar] + 1e-9).mean())
            rec["moves"][kind]["objectives"][nm] = {
                "n_local_minima": int(is_min.sum()),
                "frac_local_minima": float(is_min.mean()),
                "descent_hit_global_rate": hit,
                "mean_objective_gap": float((F[end] - F[gstar]).mean()),
                "mean_rmsd_reached": float(rmsd[end].mean()),
                "rmsd_of_global_argmin": float(rmsd[np.flatnonzero(F == F.min())].mean()),
                "autocorr_lag1_8": autocorr(F, N),
            }
            a1 = rec["moves"][kind]["objectives"][nm]["autocorr_lag1_8"][0]
            rec["moves"][kind]["objectives"][nm]["correlation_length"] = (
                float(-1.0 / np.log(max(a1, 1e-6))) if a1 > 0 else 0.0)
        del N
    return rec


def main():
    pdbs = sorted(f[len("qarch_enum_"):-4] for f in os.listdir(Q.RESULTS)
                  if f.startswith("qarch_enum_") and f.endswith(".npz"))
    rows = []
    for p in pdbs:
        rows.append(analyse(p))
        r = rows[-1]["moves"]
        print(f"  {p}: legacy local minima res {r['residual' if False else 'residue']['objectives']['legacy_total']['frac_local_minima']:.5f} "
              f"bin {r['binary']['objectives']['legacy_total']['frac_local_minima']:.5f} | "
              f"rmsd-descent reaches "
              f"res {r['residue']['objectives']['ORACLE_rmsd']['mean_rmsd_reached']:.3f} "
              f"bin {r['binary']['objectives']['ORACLE_rmsd']['mean_rmsd_reached']:.3f} "
              f"(best {rows[-1]['best_rmsd']:.3f})", flush=True)
        Q.write("qarch_landscape", {"what": "landscape ruggedness and the price of the "
                                            "encoding's move set", "rows": rows})
    # aggregate
    print(f"\n{'move set':<10s} {'objective':<17s} {'frac local min':>15s} {'descent->global':>16s} "
          f"{'mean rmsd reached':>18s} {'corr length':>12s}")
    for kind in ("residue", "binary", "gray"):
        for nm in ("legacy_total", "prior_empirical", "ORACLE_rmsd"):
            g = [r["moves"][kind]["objectives"][nm] for r in rows]
            print(f"{kind:<10s} {nm:<17s} {np.mean([x['frac_local_minima'] for x in g]):15.5f} "
                  f"{np.mean([x['descent_hit_global_rate'] for x in g]):16.3f} "
                  f"{np.mean([x['mean_rmsd_reached'] for x in g]):18.3f} "
                  f"{np.mean([x['correlation_length'] for x in g]):12.3f}")
    print(f"\npool mean rmsd {np.mean([r['pool_mean_rmsd'] for r in rows]):.3f}, "
          f"space best {np.mean([r['best_rmsd'] for r in rows]):.3f}")
    return rows


if __name__ == "__main__":
    main()
