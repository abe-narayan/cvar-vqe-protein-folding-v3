"""SPRINT 14, ENER-5 (arms F, G, J) -- SEARCH and REFINEMENT as separate questions.

AMBER may be a genuine physical refinement component even though it is a poor ranking
objective.  That is a real and useful conclusion if true, and it has to be measured on its
own axis, with the credit attributed to refinement rather than to whatever generated the
structure.  So every row reports the SAME structure at four stages:

    raw          the ideal-geometry build of the discrete configuration
    cleanup      projection back onto the ideal-geometry manifold (`I.project`), which is
                 the pipeline's own stage-3b geometric cleanup, with no energy in it
    legacy_ref   coordinate descent on the Legacy total in torsion space from that start
                 -- exact, because the full enumeration is cached (arm F)
    amber_ref    genuine restrained ff14SB/GBn2 `refine_coords(k=10, steps=0)`, 9.1 s

and five different STARTS, so the refinement effect can be separated from the start:

    random | prior argmin | Legacy argmin | AMBER argmin | ORACLE pool best

Arms F and G are the search arms: F is Legacy's CERTIFIED global optimum over all 262,144
configurations (free -- it is a cached argmin, and it is tie-averaged); G is coordinate
descent on genuine AMBER (28 ms per distinct evaluation, so ~100 evaluations per start).

    python -m s14.ener_refine
"""
from __future__ import annotations

import time

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E
from s14.ener_matrix import ensemble

STARTS = ["random", "prior_argmin", "legacy_argmin", "amber_argmin", "ORACLE_pool_best"]


# ------------------------------------------------------------------ arm F: exact Legacy
def legacy_descent(z, start_idx, max_sweeps=12):
    """Coordinate descent on the Legacy total, exact via the cached enumeration."""
    s = z.states(np.array([start_idx]))[0].astype(int)
    cur = float(z.legacy[start_idx])
    n, k = z.n, z.k
    ev = 0
    for _ in range(max_sweeps):
        moved = False
        for i in range(n):
            cand = np.repeat(s[None], k, axis=0)
            cand[:, i] = np.arange(k)
            idx = z.index_of(cand)
            vals = z.legacy[idx]
            ev += k
            j = int(np.argmin(vals))
            if vals[j] < cur - 1e-12:
                cur = float(vals[j]); s = cand[j]; moved = True
        if not moved:
            break
    return int(z.index_of(s[None])[0]), cur, ev


def legacy_certified(z):
    """Arm F at its limit: the CERTIFIED global optimum over the whole space."""
    return E.argmin_rmsd(z.legacy, z.rmsd), int(np.argmin(z.legacy))


# ------------------------------------------------------------------ arm G: AMBER descent
def amber_descent(z, start_idx, max_sweeps=3, budget=400):
    """Coordinate descent on GENUINE AMBER single points. 28 ms per distinct evaluation."""
    from s13 import qarch_lib as Q
    sp = z.space()
    s = z.states(np.array([start_idx]))[0].astype(int)
    memo = {}

    def ev(state):
        key = tuple(int(x) for x in state)
        if key not in memo:
            memo[key] = float(Q.amber_energies(sp, np.array(state)[None])[0])
        return memo[key]

    cur = ev(s)
    n, k = z.n, z.k
    for _ in range(max_sweeps):
        moved = False
        for i in range(n):
            if len(memo) >= budget:
                break
            best_j, best_v = int(s[i]), cur
            for kk in range(k):
                if kk == s[i]:
                    continue
                t = s.copy(); t[i] = kk
                v = ev(t)
                if v < best_v:
                    best_v, best_j = v, kk
            if best_j != s[i]:
                s[i] = best_j; cur = best_v; moved = True
        if not moved:
            break
    return int(z.index_of(s[None])[0]), cur, len(memo)


# ------------------------------------------------------------------ arm J: refinement
def amber_refine_one(z, idx):
    """Genuine restrained ff14SB/GBn2 refinement of one configuration. ~9.1 s."""
    from core import amber as am
    sp = z.space()
    s = z.states(np.array([idx]))[0].astype(int)
    bits = sp.rep.bitstring_from_states(np.asarray(s, int))
    coords = sp.rep.build_coords(bits)
    E.wait_for_memory(1.5, "amber-refine")
    t0 = time.time()
    r = am.refine_coords(sp.seq, sp.rep, coords, k_restraint=10.0, steps=0,
                         threads=1, memo=True)
    ca = np.asarray(r["ca"], float)
    return dict(rmsd=float(I.ca_rmsd(ca, sp.nat_ca)),
                energy=float(r["energy"]), energy_initial=float(r["energy_initial"]),
                restraint_rmsd=float(r.get("restraint_rmsd", np.nan)),
                wall=round(time.time() - t0, 2))


def cleanup_one(z, idx):
    """Stage-3b geometric cleanup: projection back onto the ideal-geometry manifold."""
    sp = z.space()
    ca = z.ca(np.array([idx]))[0]
    p = I.project(ca, sp.seq, sp.fold)
    return dict(rmsd=float(I.ca_rmsd(p["ca"], sp.nat_ca)),
                fit_rmsd=float(I.ca_rmsd(p["fit_ca"], sp.nat_ca)))


# ------------------------------------------------------------------ main
def starts_for(pdb):
    d = ensemble(pdb)
    z = d["_z"]
    R = d["rmsd"]
    rng = np.random.default_rng(E.SEED)
    out = {}
    out["random"] = int(d["idx"][rng.integers(0, d["n"])])
    for nm, key in (("prior_argmin", "prior"), ("legacy_argmin", "legacy"),
                    ("amber_argmin", "amber")):
        v = np.asarray(d[key], float)
        tie = np.flatnonzero(v == np.nanmin(v))
        out[nm] = int(d["idx"][tie[rng.integers(0, len(tie))]])
    out["ORACLE_pool_best"] = int(d["idx"][int(np.argmin(R))])
    return z, out


def main():
    rows = []
    for pdb in E.ENUM_TARGETS:
        z, st = starts_for(pdb)
        cert_rmsd, cert_idx = legacy_certified(z)
        rec = dict(pdb=pdb, fold=int(z.fold), starts={},
                   F_legacy_certified_rmsd=cert_rmsd,
                   F_legacy_certified_idx=cert_idx,
                   pool_mean=float(z.rmsd.mean()), space_best=float(z.rmsd.min()))
        for nm, idx in st.items():
            raw = float(z.rmsd[idx])
            li, lv, lev = legacy_descent(z, idx)
            gi, gv, gev = amber_descent(z, idx)
            cl = cleanup_one(z, idx)
            ar = amber_refine_one(z, idx)
            rec["starts"][nm] = dict(
                start_idx=idx, raw_rmsd=raw,
                cleanup_rmsd=cl["rmsd"], cleanup_fit_rmsd=cl["fit_rmsd"],
                F_legacy_descent_rmsd=float(z.rmsd[li]), F_legacy_evals=lev,
                G_amber_descent_rmsd=float(z.rmsd[gi]), G_amber_evals=gev,
                J_amber_refine_rmsd=ar["rmsd"], J_wall=ar["wall"],
                J_energy=ar["energy"], J_energy_initial=ar["energy_initial"])
            print(f"{pdb} {nm:18s} raw {raw:.3f} -> cleanup {cl['rmsd']:.3f} "
                  f"legacyF {z.rmsd[li]:.3f} amberG {z.rmsd[gi]:.3f} "
                  f"amberJ {ar['rmsd']:.3f}  ({ar['wall']:.1f}s, "
                  f"{gev} AMBER evals)", flush=True)
        rows.append(rec)
        E.write("ener_refine", dict(
            what="arms F (Legacy search), G (AMBER search), J (AMBER refinement)",
            per_target=rows), n_expected=len(E.ENUM_TARGETS))
    report(rows)
    return rows


def report(rows):
    print("\n=== ARMS F / G / J, pooled over the completed targets ===")
    folds = [r["fold"] for r in rows]
    print(f"{'start':20s} {'raw':>7s} {'cleanup':>8s} {'F leg':>7s} {'G amb':>7s} "
          f"{'J refine':>9s} {'J - raw':>9s} {'CI':>20s} {'W/L':>6s}")
    for nm in STARTS:
        raw = np.array([r["starts"][nm]["raw_rmsd"] for r in rows])
        cl = np.array([r["starts"][nm]["cleanup_rmsd"] for r in rows])
        f = np.array([r["starts"][nm]["F_legacy_descent_rmsd"] for r in rows])
        g = np.array([r["starts"][nm]["G_amber_descent_rmsd"] for r in rows])
        j = np.array([r["starts"][nm]["J_amber_refine_rmsd"] for r in rows])
        p = I.paired(j, raw, folds=folds)
        print(f"{nm:20s} {raw.mean():7.3f} {cl.mean():8.3f} {f.mean():7.3f} "
              f"{g.mean():7.3f} {j.mean():9.3f} {p['mean_diff']:+9.3f} "
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}")
    print(f"\nArm F at its LIMIT (Legacy CERTIFIED global optimum over all 262,144): "
          f"{np.mean([r['F_legacy_certified_rmsd'] for r in rows]):.3f}   "
          f"vs whole-space mean {np.mean([r['pool_mean'] for r in rows]):.3f}   "
          f"vs space best {np.mean([r['space_best'] for r in rows]):.3f}")
    print("\n--- the separate axes: does refinement help INDEPENDENT of the start? ---")
    allraw = np.concatenate([[r["starts"][nm]["raw_rmsd"] for r in rows]
                             for nm in STARTS])
    allj = np.concatenate([[r["starts"][nm]["J_amber_refine_rmsd"] for r in rows]
                           for nm in STARTS])
    p = I.paired(allj, allraw)
    print(f"AMBER refinement over all {len(allraw)} (target, start) cells: "
          f"{p['mean_a']:.3f} vs {p['mean_b']:.3f}, {p['mean_diff']:+.3f} "
          f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}  "
          f"drop-top10 {p['drop_top10_mean_diff']}")
    allc = np.concatenate([[r["starts"][nm]["cleanup_rmsd"] for r in rows]
                           for nm in STARTS])
    p = I.paired(allc, allraw)
    print(f"geometric cleanup only:                         "
          f"{p['mean_a']:.3f} vs {p['mean_b']:.3f}, {p['mean_diff']:+.3f} "
          f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}")
    allf = np.concatenate([[r["starts"][nm]["F_legacy_descent_rmsd"] for r in rows]
                           for nm in STARTS])
    allg = np.concatenate([[r["starts"][nm]["G_amber_descent_rmsd"] for r in rows]
                           for nm in STARTS])
    for nm, v in (("Legacy descent (arm F)", allf), ("AMBER descent (arm G)", allg)):
        p = I.paired(v, allraw)
        print(f"{nm:32s} {p['mean_a']:.3f} vs {p['mean_b']:.3f}, "
              f"{p['mean_diff']:+.3f} [{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
              f"{p['n_better']}/{p['n_worse']}")


if __name__ == "__main__":
    main()
