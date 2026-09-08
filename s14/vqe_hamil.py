"""SPRINT 14 / VQE -- the coordinator's structural objective against a CERTIFIED optimum.

`s14/hamil.py` is a native-free STRUCTURAL Hamiltonian (retrieval-conditioned 1-local
torsion prior + leave-fold-out distogram Bayes risk).  The coordinator measured it on a
20,000-configuration sample per target and found the objective improves monotonically with
budget while the emitted structure does not move.

That result can always be answered with "you did not search far enough".  Here it cannot:
this module TABULATES the objective over ALL 262,144 configurations of the nine fully
enumerated targets, so its **certified global optimum is known exactly** and the budget
defence is closed.

Reported on all four axes the shared brief now requires -- global rho, in-decile rho, argmin
RMSD, decile-mean RMSD -- because a single in-decile number is ambiguous between objectives
differing several-fold in real quality (see `s14/vqe_decile.py`).

Then every Part C arm is run on it at matched budget, including the CVaR-VQE arms.

    python -m s14.vqe_hamil
"""
from __future__ import annotations

import os
import time

import numpy as np

from s12 import instrument as I
from s14 import vqe_lib as V
from s14 import vqe_signal as SG

WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)


def tabulate(pdb, chunk=8192, cache=True):
    """Full-enumeration tabulation of both hamil terms. Cached to `s14/cache/`."""
    path = os.path.join(V.CACHE, f"hamil_full_{pdb}.npz")
    if cache and os.path.exists(path):
        d = np.load(path)
        return d["prior"], d["disto"]
    from s14 import hamil as H
    e = V.Enum(pdb)
    t = H.Terms(pdb, e.seq, e.n, e.fold, k=e.k)
    pr = np.empty(e.N, np.float32)
    ds = np.empty(e.N, np.float32)
    t0 = time.time()
    for a in range(0, e.N, chunk):
        S = e.states(np.arange(a, min(a + chunk, e.N)))
        pr[a:a + len(S)] = t.e_prior(S)
        ds[a:a + len(S)] = t.e_disto(S)
        if a % (chunk * 8) == 0:
            print(f"    {pdb} {a}/{e.N}  {time.time()-t0:.0f}s", flush=True)
    np.savez_compressed(path, prior=pr, disto=ds)
    return pr, ds


def combine(pr, ds, w):
    """Per-target standardisation then a weighted sum -- hamil's own normalisation."""
    def z(x):
        x = np.asarray(x, float)
        s = x.std()
        return (x - x.mean()) / (s if s > 1e-12 else 1.0)
    return (1.0 - w) * z(pr) + w * z(ds)


def in_decile_rho(E, rmsd, frac=0.1):
    k = max(3, int(round(frac * len(E))))
    idx = np.argpartition(E, k - 1)[:k]
    return V.spearman(E[idx], rmsd[idx]), float(rmsd[idx].mean()), float(rmsd[idx].min())


def main():
    V.wait_for_memory(1.2, "vqe_hamil")
    rows = {}
    print("tabulating hamil over the FULL enumeration (262,144 configurations/target)")
    for pdb in V.ENUM_TARGETS:
        rows[pdb] = tabulate(pdb)
        print(f"  {pdb} done", flush=True)

    print()
    print("=" * 104)
    print("A. hamil ON THE CERTIFIED ENUMERATION -- all four axes, nine targets")
    print("=" * 104)
    print(f"  {'w':>5s} {'global rho':>11s} {'in-decile rho':>14s} "
          f"{'RMSD @ certified argmin':>24s} {'decile mean':>12s} {'decile best':>12s}")
    out = {"axes": {}}
    for w in WEIGHTS:
        g, dr, am, dm, db = [], [], [], [], []
        for pdb in V.ENUM_TARGETS:
            e = V.Enum(pdb)
            E = combine(*rows[pdb], w)
            g.append(V.spearman(E, e.rmsd))
            r, mn, bs = in_decile_rho(E, e.rmsd)
            dr.append(r); dm.append(mn); db.append(bs)
            am.append(float(e.rmsd[int(np.argmin(E))]))
        out["axes"][str(w)] = {"global_rho": float(np.mean(g)),
                               "in_decile_rho": float(np.mean(dr)),
                               "rmsd_at_argmin": float(np.mean(am)),
                               "decile_mean_rmsd": float(np.mean(dm)),
                               "decile_best_rmsd": float(np.mean(db)),
                               "per_target_argmin": am}
        print(f"  {w:5.2f} {np.mean(g):+11.3f} {np.mean(dr):+14.3f} "
              f"{np.mean(am):24.3f} {np.mean(dm):12.3f} {np.mean(db):12.3f}")
    e0 = [V.Enum(p) for p in V.ENUM_TARGETS]
    print(f"\n  reference: space best {np.mean([e.rmsd.min() for e in e0]):.3f}, "
          f"random draw {np.mean([e.rmsd.mean() for e in e0]):.3f}, "
          f"Legacy certified argmin 3.920, torsion-prior certified argmin 3.636")
    print("\n  THE SELECTION GAP AT INFINITE BUDGET is (RMSD @ certified argmin) minus")
    print("  (decile best): the objective's own optimum against the best structure inside")
    print("  its own lowest decile.  No amount of search can close it.")
    for w in WEIGHTS:
        a = out["axes"][str(w)]
        print(f"    w={w:4.2f}  {a['rmsd_at_argmin']:.3f} - {a['decile_best_rmsd']:.3f} "
              f"= {a['rmsd_at_argmin'] - a['decile_best_rmsd']:+.3f} A")

    # ---------------------------------------------------------- the arm comparison
    print()
    print("=" * 104)
    print("B. EVERY PART C ARM ON hamil (w=0.25), matched budgets, certified optimum known")
    print("=" * 104)
    res = {"config": {"objective": "s14.hamil", "space": "qarch_enum n=9 k=4"},
           "cells": {}}
    _p = os.path.join(V.RESULTS, "vqe_signal_hamil.json")
    if os.path.exists(_p):          # merge, so the w arms accumulate across invocations
        import json
        res["cells"] = json.load(open(_p))["cells"]
        print(f"resuming: {len(res['cells'])} hamil cells already present", flush=True)
    import sys
    buds = (30000,) if "--fast" in sys.argv else (3000, 30000)
    # `--w 0,1` runs the SIGN TEST the coordinator asked for: if the tail-concentration
    # mechanism were right, a small-alpha CVaR advantage should appear at w=0 (best argmin,
    # chance-level bulk) and be absent at w=1 (best bulk, anti-ranking tail).  My revised
    # prediction after `s14/vqe_tailacc.py` is that it appears at NEITHER, because no
    # objective in the family is tail-concentrated.
    ws = [float(x) for x in sys.argv[sys.argv.index("--w") + 1].split(",")] \
        if "--w" in sys.argv else [0.25]
    for w in ws:
      for bud in buds:
        for pdb in V.ENUM_TARGETS:
            e = V.Enum(pdb)
            E = combine(*rows[pdb], w)
            key = f"{pdb}|hamil{w}|0.0|{bud}"
            if key in res["cells"]:
                continue
            res["cells"][key] = SG.cell(
                e, E, bud, seeds=(0, 1) if "--fast" in sys.argv else (0, 1, 2))
            r = res["cells"][key]
            print(f"  {key:30s} rho={r['rho_obj_rmsd']:+.3f} "
                  f"exact={r['exact']['rmsd_returned']:.3f} "
                  f"rnd={r['random']['rmsd_returned']:.3f} "
                  f"sa={r['anneal']['rmsd_returned']:.3f} "
                  f"vqe1={r['vqe_a1.0']['rmsd_returned']:.3f} "
                  f"vqe.05={r['vqe_a0.05']['rmsd_returned']:.3f}", flush=True)
            V.write("vqe_signal_hamil", res)
    V.write("vqe_signal_hamil", res)
    V.write("vqe_hamil_axes", out)
    print("\nwritten -> s14/results/vqe_hamil_axes.json, vqe_signal_hamil.json")
    print("analyse with:  python -m s14.vqe_report --base hamil")


if __name__ == "__main__":
    main()
