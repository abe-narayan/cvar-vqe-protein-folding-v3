"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 7 -- does the objective rank at scale?

Sections 3a-3e are exact but confined to nine 9-residue targets at k=4, and the coordinator's
independent random-sampling measurement reports Legacy rho around +0.3 to +0.4 where I
measure +0.081.  A disagreement that large has to be resolved rather than averaged, and the
obvious candidate is target LENGTH: RMSD spread grows with n, so a correlation measured on
9-mers and on 16-mers are not the same quantity.

So: ALL 126 tuning targets, k in {4, 8}, `n_draw` uniform-random configurations each,
genuine Legacy (11 terms, DEFAULT_WEIGHTS) and the leakage-safe 1-local prior, scored against
true CA-RMSD.  Reported per target and binned by n, with the same statistics as section 3 --
Spearman, the argmin's RMSD against a random draw, the native-snap percentile -- plus the
`no_steric` variant, because the coordinator finds AMBER's rho recovers when its clash term
is removed and the same question has to be asked of Legacy.

    python -m s13.qarch_scale [n_draw] [workers]
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

N_DRAW = 12000
KS = (4, 8)


def one(args):
    pdb_id, seq, n, fold, n_draw, seed = args
    from s13 import qarch_lib as Q
    out = []
    for k in KS:
        sp = Q.Space(pdb_id, k, seq=seq, n=n, fold=fold)
        rng = np.random.default_rng(seed)
        S = sp.uniform(n_draw, rng)
        snap = sp.ORACLE_snap()
        S = np.vstack([S, snap[None, :]])
        r = sp.rmsd(S)
        comp = Q.legacy_components(sp, S, chunk=4096)
        tot = Q.legacy_total(comp)
        no_st = tot - float(Q.et.DEFAULT_WEIGHTS.get("steric", 0.0)) * comp["steric"]
        P = Q.empirical_prior(sp)
        pri = Q.prior_energy(P, S)
        rec = {"pdb": pdb_id, "n": n, "k": k, "fold": fold, "n_draw": int(len(S)),
               "rmsd_mean": float(r.mean()), "rmsd_min": float(r.min()),
               "rmsd_sd": float(r.std()), "snap_rmsd": float(r[-1]),
               "objectives": {}}
        for nm, E in (("legacy_total", tot), ("legacy_no_steric", no_st),
                      ("leg_steric", comp["steric"]), ("leg_torsion", comp["torsion"]),
                      ("prior_empirical", pri)):
            E = np.asarray(E, float)
            tie = np.flatnonzero(E == E.min())
            rec["objectives"][nm] = {
                "spearman": Q.spearman(E, r),
                "argmin_rmsd": float(r[tie].mean()),
                "delta_vs_random": float(r[tie].mean() - r.mean()),
                "native_snap_pct": Q.percentile_of(E[-1], E),
                "top10_rmsd": float(r[np.argsort(E)[:10]].mean()),
            }
        out.append(rec)
    return out


def main(n_draw=N_DRAW, workers=6):
    from s12 import instrument as I
    tg = I.targets()
    args = [(t["pdb"], t["seq"], t["n"], t["fold"], n_draw, 0) for t in tg]
    rows = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, res in enumerate(ex.map(one, args, chunksize=1)):
            rows += res
            if i % 10 == 0:
                print(f"  {i+1}/{len(args)} [{time.time()-t0:.0f}s]", flush=True)
                _write(rows, n_draw)
    _write(rows, n_draw)
    report(rows)
    return rows


def _write(rows, n_draw):
    from s13 import qarch_lib as Q
    Q.write("qarch_scale", {"what": "does Legacy/prior rank the native across all 126 "
                                    "tuning targets and k in {4,8}?",
                            "n_draw": n_draw, "rows": rows})


def report(rows):
    import numpy as np
    for k in KS:
        R = [r for r in rows if r["k"] == k]
        if not R:
            continue
        print(f"\n=== k={k}, {len(R)} targets, {R[0]['n_draw']} random configurations each ===")
        print(f"{'objective':<20s} {'mean rho':>9s} {'median':>8s} {'sd':>7s} {'argmin':>8s} "
              f"{'random':>8s} {'delta':>8s} {'W/L':>8s} {'nat pct':>8s}")
        for nm in ("legacy_total", "legacy_no_steric", "leg_steric", "leg_torsion",
                   "prior_empirical"):
            rho = np.array([r["objectives"][nm]["spearman"] for r in R])
            am = np.array([r["objectives"][nm]["argmin_rmsd"] for r in R])
            pm = np.array([r["rmsd_mean"] for r in R])
            npc = np.array([r["objectives"][nm]["native_snap_pct"] for r in R])
            d = am - pm
            print(f"{nm:<20s} {np.nanmean(rho):+9.3f} {np.nanmedian(rho):+8.3f} "
                  f"{np.nanstd(rho):7.3f} {am.mean():8.3f} {pm.mean():8.3f} {d.mean():+8.3f} "
                  f"{int((d<0).sum())}/{int((d>0).sum()):<4d} {npc.mean():8.3f}")
        print("  Legacy rho and delta binned by n:")
        n = np.array([r["n"] for r in R])
        rho = np.array([r["objectives"]["legacy_total"]["spearman"] for r in R])
        am = np.array([r["objectives"]["legacy_total"]["argmin_rmsd"] for r in R])
        pm = np.array([r["rmsd_mean"] for r in R])
        for v in sorted(set(n.tolist())):
            m = n == v
            print(f"    n={v:2d} ({m.sum():2d} targets)  rho={rho[m].mean():+.3f}  "
                  f"argmin={am[m].mean():.3f}  random={pm[m].mean():.3f}  "
                  f"delta={(am[m]-pm[m]).mean():+.3f}")


if __name__ == "__main__":
    nd = int(sys.argv[1]) if len(sys.argv) > 1 else N_DRAW
    w = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    main(nd, w)
