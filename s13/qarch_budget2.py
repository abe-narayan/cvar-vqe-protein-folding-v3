"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 8b -- the budget curve on all 126 targets.

`qarch_budget` measured, on the nine completely enumerated n=9 spaces, that the CA-RMSD
reached by minimising Legacy is NON-MONOTONE in the number of objective evaluations: it
improves to a minimum around 300 evaluations and then gets steadily worse, ending 0.25 A
worse at the certified global optimum than at best-of-300.

This file asks whether that shape holds on the whole 126-target tuning instrument, where the
space cannot be enumerated: draw `M` configurations per target, then read the best-of-m curve
off the same sample for m = 10 ... M.  Same estimator as `qarch_budget`, no enumeration
needed, and directly comparable to the k=4 row of `qarch_scale`.

    python -m s13.qarch_budget2 [M] [workers]
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

M_DEFAULT = 60000
BUDGETS = (10, 30, 100, 300, 1000, 3000, 10000, 30000, 60000)
REPS = 24
K = 4


def one(args):
    pdb_id, seq, n, fold, M, seed = args
    from s13 import qarch_lib as Q
    sp = Q.Space(pdb_id, K, seq=seq, n=n, fold=fold)
    rng = np.random.default_rng(seed)
    S = sp.uniform(M, rng)
    r = sp.rmsd(S)
    comp = Q.legacy_components(sp, S, chunk=8192)
    tot = Q.legacy_total(comp)
    P = Q.empirical_prior(sp)
    pri = Q.prior_energy(P, S)
    rec = {"pdb": pdb_id, "n": n, "k": K, "M": int(M),
           "pool_mean": float(r.mean()), "sample_best": float(r.min()), "curves": {}}
    for nm, F in (("legacy_total", tot), ("prior_empirical", pri),
                  ("leg_steric", comp["steric"]), ("ORACLE_rmsd", r)):
        F = np.asarray(F, float)
        cur = {}
        for m in BUDGETS:
            if m > M:
                continue
            vals = []
            for _ in range(REPS if m < M else 1):
                idx = rng.choice(M, m, replace=False) if m < M else np.arange(M)
                E = F[idx]
                tie = np.flatnonzero(E == E.min())
                vals.append(float(r[idx][tie].mean()))
            cur[m] = float(np.mean(vals))
        rec["curves"][nm] = cur
    return rec


def main(M=M_DEFAULT, workers=6):
    from s12 import instrument as I
    from s13 import qarch_lib as Q
    tg = I.targets()
    args = [(t["pdb"], t["seq"], t["n"], t["fold"], M, 0) for t in tg]
    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, rec in enumerate(ex.map(one, args, chunksize=1)):
            rows.append(rec)
            if i % 15 == 0:
                print(f"  {i+1}/{len(args)} [{time.time()-t0:.0f}s]", flush=True)
                Q.write("qarch_budget2", {"what": "best-of-m CA-RMSD vs objective-evaluation "
                                                 "budget, 126 targets, k=4",
                                          "M": M, "reps": REPS, "rows": rows})
    Q.write("qarch_budget2", {"what": "best-of-m CA-RMSD vs objective-evaluation budget, "
                                      "126 targets, k=4", "M": M, "reps": REPS, "rows": rows})
    report(rows)
    return rows


def report(rows):
    names = list(rows[0]["curves"])
    print(f"\n{'budget':>8s} " + " ".join(f"{nm:>16s}" for nm in names))
    for m in BUDGETS:
        line = f"{m:8d} "
        for nm in names:
            v = [r["curves"][nm][m] for r in rows if m in r["curves"][nm]]
            line += f"{np.mean(v):16.3f} " if v else " " * 17
        print(line)
    print(f"{'pool mean':>8s} {np.mean([r['pool_mean'] for r in rows]):16.3f}")
    print(f"{'sampled best':>8s} {np.mean([r['sample_best'] for r in rows]):16.3f}")
    print("\nLegacy curve by n:")
    ns = sorted({r["n"] for r in rows})
    print(f"{'n':>3s} " + " ".join(f"{m:>8d}" for m in BUDGETS))
    for n in ns:
        R = [r for r in rows if r["n"] == n]
        print(f"{n:3d} " + " ".join(
            f"{np.mean([r['curves']['legacy_total'][m] for r in R]):8.3f}" for m in BUDGETS))


if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else M_DEFAULT
    w = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    main(M, w)
