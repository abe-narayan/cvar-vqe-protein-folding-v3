"""S25 / PHYSICS LANE / F0 -- THE CACHE-INTEGRITY GATE. Pre-registered in `s25/PREREG_PHYS.md` S4.

    python s25/phys_gate.py

WHAT THIS IS FOR, AND WHY IT IS A SEPARATE PROGRAM
==================================================
`s24/cache_amber/` holds 63,000 genuine ff14SB/GBn2 single points and the brief instructs this
lane to USE it rather than re-derive it at OpenMM cost. A cache is only worth what its
integrity check is worth, so the check runs first, alone, and nothing downstream is read until
it passes. Three assertions, in increasing cost:

  A. POOL IDENTITY, all 126 targets, free.  The cached `universe_idx` must reproduce the live
     `Candidates.from_universe(pdb, k=500)` bit-for-bit. If the retrieval pool has moved since
     Sprint 24, every AMBER number in the cache belongs to different structures.
  B. SCORE REPRODUCTION, all 126 targets, cheap.  The cached distogram score must reproduce
     bit-for-bit under today's code, which pins the distogram models and the shipped score.
  C. ENERGY BIT-EXACTNESS, 5 fold-stratified targets x 4 candidates, under LOCK_AMBER.  A
     FRESH `AmberSP` single point is computed and compared to the cached `e_amber`. This is
     the only OpenMM call this lane makes.

`amber_verify_max_rel` recorded in the cache is Sprint 24's own gate (AmberSP vs
`core.amber.refine_coords`) and is 0.0 on every target; it is reported, but it is NOT this
check -- a stored flag proves nothing about whether the stored ENERGIES are still reproducible.
`pauli-spectrum-delta-spike-artefact`'s lesson generalises: a gate that cannot fire is not
evidence, so this one recomputes rather than re-reads.

THE CONDITIONING CENSUS RIDES ALONG
===================================
The same pass reports the native-free distribution facts the pre-registration's normalisation
fork rests on -- what share of the pool exceeds 1e4 kcal/mol, what share lands inside
|z_raw| < 0.1, and the worst single point. These are diagnosed from the energy distribution
alone, with no outcome variable inspected, which is the order s24 D1 got right.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402

N_SPOT_TARGETS = 5
N_SPOT_CANDS = 4


def main() -> int:
    pdbs = P.targets()
    print(f"F0 CACHE-INTEGRITY GATE -- {len(pdbs)} targets\n")

    # ---------------------------------------------------- A + B, all 126, no OpenMM
    rows, t0 = [], time.time()
    worst_raw, frac_hi, frac_z, ks = 0.0, [], [], []
    for ii, pdb in enumerate(pdbs):
        z = np.load(os.path.join(P.CACHE_AMB, f"{pdb}.npz"), allow_pickle=True)
        cand = H.Candidates.from_universe(pdb, k=int(z["k"]))
        ui_ok = bool(np.array_equal(np.asarray(cand.meta["universe_idx"], int),
                                    np.asarray(z["universe_idx"], int)))
        sc = H.score_shipped(cand)
        sc_ok = bool(np.array_equal(sc, np.asarray(z["score_dist"], float)))
        e = np.asarray(z["e_amber"], float)
        zr = (e - e.mean()) / max(e.std(), 1e-12)
        rows.append(dict(pdb=pdb, n=int(z["n"]), fold=int(z["fold"]), k=int(z["k"]),
                         pool_identical=ui_ok, score_reproduces=sc_ok,
                         amber_verify_max_rel=float(z["amber_verify_max_rel"]),
                         e_min=float(e.min()), e_max=float(e.max()),
                         frac_gt_1e4=float((e > 1e4).mean()),
                         frac_absz_lt_0p1=float((np.abs(zr) < 0.1).mean()),
                         legacy_min=float("nan"), legacy_max=float("nan")))
        L = P.legacy_cached(cand)
        rows[-1]["legacy_min"] = float(L.min()); rows[-1]["legacy_max"] = float(L.max())
        worst_raw = max(worst_raw, float(e.max()))
        frac_hi.append(rows[-1]["frac_gt_1e4"]); frac_z.append(rows[-1]["frac_absz_lt_0p1"])
        ks.append(int(z["k"]))
        if (ii + 1) % 25 == 0:
            print(f"  A/B {ii+1}/{len(pdbs)}  {time.time()-t0:.0f}s", flush=True)

    a_ok = all(r["pool_identical"] for r in rows)
    b_ok = all(r["score_reproduces"] for r in rows)
    print(f"\n  A pool identity      : {sum(r['pool_identical'] for r in rows)}/{len(rows)}  "
          f"{'PASS' if a_ok else 'FAIL'}")
    print(f"  B score reproduction : {sum(r['score_reproduces'] for r in rows)}/{len(rows)}  "
          f"{'PASS' if b_ok else 'FAIL'}")
    print(f"  cached amber_verify_max_rel: max over 126 = "
          f"{max(r['amber_verify_max_rel'] for r in rows):.3e}")
    print(f"  k: {sorted(set(ks))}   total single points = {sum(ks)}")

    print(f"\n  CONDITIONING CENSUS (native-free, no outcome inspected)")
    print(f"    worst single AMBER point over the 63,000 : {worst_raw:.3e} kcal/mol")
    print(f"    mean share of pool with E > 1e4 kcal/mol : {np.mean(frac_hi):.4f}")
    print(f"    mean share of pool inside |z_raw| < 0.1  : {np.mean(frac_z):.4f}")
    print(f"    -> a MOMENT z-score of raw AMBER measures which candidate has the worst clash.")

    # ---------------------------------------------------- C, fresh single points
    by_fold = {}
    for r in rows:
        by_fold.setdefault(r["fold"], []).append(r["pdb"])
    spot = [sorted(v)[0] for _, v in sorted(by_fold.items())][:N_SPOT_TARGETS]
    print(f"\n  C fresh-single-point bit-exactness on {spot} "
          f"x {N_SPOT_CANDS} candidates, under LOCK_AMBER")
    spot_rows, c_ok = [], True
    from s13 import qarch_lib as QA
    from s20 import qb2_lib as QB2
    t1 = time.time()
    with P.AmberLock(tag=f"F0 spot check {len(spot)}x{N_SPOT_CANDS}"):
        for pdb in spot:
            z = np.load(os.path.join(P.CACHE_AMB, f"{pdb}.npz"), allow_pickle=True)
            cand = H.Candidates.from_universe(pdb, k=int(z["k"]))
            sp = QB2.AmberSP(cand.seq, QA.Space(pdb, 4).rep)
            ver, nver = sp.verify(cand.PHI[:N_SPOT_CANDS], cand.PSI[:N_SPOT_CANDS],
                                  m=N_SPOT_CANDS)
            fresh = sp.batch(cand.PHI[:N_SPOT_CANDS], cand.PSI[:N_SPOT_CANDS])
            cached = np.asarray(z["e_amber"], float)[:N_SPOT_CANDS]
            rel = np.abs(fresh - cached) / np.maximum(1.0, np.abs(cached))
            bit = bool(np.array_equal(fresh, cached))
            c_ok = c_ok and bit
            spot_rows.append(dict(pdb=pdb, n_cands=N_SPOT_CANDS, bit_identical=bit,
                                  max_rel=float(rel.max()),
                                  refine_coords_verify_max_rel=float(ver),
                                  n_verify=int(nver),
                                  fresh=fresh.tolist(), cached=cached.tolist()))
            print(f"    {pdb}: bit-identical={bit}  max_rel={rel.max():.3e}  "
                  f"AmberSP-vs-refine_coords={ver:.3e} over {nver} cmp", flush=True)
    print(f"  C: {'PASS' if c_ok else 'FAIL'}   ({time.time()-t1:.0f}s under lock)")

    verdict = "PASS" if (a_ok and b_ok and c_ok) else "FAIL"
    print(f"\n  F0 GATE: {verdict}")

    out = dict(
        gate="F0_cache_integrity", verdict=verdict,
        a_pool_identity=a_ok, b_score_reproduction=b_ok, c_energy_bit_exact=c_ok,
        n_single_points=int(sum(ks)),
        conditioning=dict(worst_single_point=float(worst_raw),
                          mean_frac_gt_1e4=float(np.mean(frac_hi)),
                          mean_frac_absz_lt_0p1=float(np.mean(frac_z))),
        spot=spot_rows, rows=rows,
        alpha=P.ALPHA, T=P.TEMP, layers=P.LAYERS, iters=P.ITERS, seed=P.SEED,
        m_prod=P.M_PROD, k=P.K)
    ST.save_atomic(os.path.join(P.RESULTS, "phys_gate.json"), out,
                   complete_keys=("pdb", "fold", "pool_identical", "score_reproduces",
                                  "frac_gt_1e4", "frac_absz_lt_0p1"),
                   rows=rows, n_expected=len(pdbs), module_file=__file__)
    print(f"  wrote s25/results/phys_gate.json")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
