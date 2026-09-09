"""S24 / LANE D / D2 -- AN AMBER-INFORMED SCORE. Pre-registered in `s24/PREREG_D.md` (D2).

    python s24/d2_amberscore.py --phase energies    # takes LOCK_AMBER, persists E_AMBER
    python s24/d2_amberscore.py --phase analyse     # no lock, no OpenMM, reads the cache

THE QUESTION, AND WHY IT IS THE LAST ONE IN THIS LANE
=====================================================
D1's mechanism finding: `rho(Legacy, distogram) = +0.3875` but `rho(AMBER, distogram) = -0.0270`.
**AMBER is the one available functional orthogonal to the thing that already selects**, and
orthogonality is the only structural reason anyone has found to expect a different functional to
help. D1 tested AMBER as a PRE-FILTER and as a PARTITION. Both failed, and both failed the same
way -- they let AMBER *choose candidates*, importing AMBER's quality problem wholesale
(`AMBER_PREFERS`: 4.13 A, q = 1.520).

**D2 asks whether the orthogonal DIRECTION can be inherited without the QUALITY being
inherited**, by letting AMBER contribute to the SCORE while the distogram keeps doing the
choosing:

    s(w) = zrank(s_distogram) + w * zrank(E_AMBER)          w = 0 IS the incumbent, exactly

`zrank` is strictly monotone, so at `w = 0` the ordering is bit-identical to the shipped score's
and the arm reproduces the incumbent rather than approximating it. Rank-conditioning BOTH sides
makes `w` dimensionless and comparable across targets, and it is mandatory on the AMBER side --
per D1-A, raw ff14SB/GBn2 single points on unrelaxed windows are clash-dominated and a moment
z-score puts 99.4% of every pool inside |z| < 0.1.

REGISTERED PRIOR: I EXPECT THIS TO FAIL
=======================================
`conf.py` fitted one global exponent under nested CV over the *same* functional and moved the
endpoint +0.0314 with 33W/93L. D1-C showed every AMBER-flavoured candidate set costs 0.26-0.55 A.
A blend is a monotone-preserving perturbation of an ordering whose top-75 is already known to be
composition-insensitive (L2(d), +0.0022). **Expected: a null at small |w|, harm at large |w|,
nested CV selecting w ~ 0. If nested CV picks w = 0, that is the result and it closes the lane.**

WHAT IS GENUINE
===============
`E_AMBER` is a genuine ff14SB/GBn2 **single point** via `s20.qb2_lib.AmberSP`, whose `verify()`
asserts bit-exactness against `core.amber.refine_coords(k_restraint=0, steps=-1)` before any
number is read. **No minimisation, no relaxation, no repair** -- s23 L11 measured 17 of 17
repair settings at or worse than the no-op. AMBER is an ENERGY MEASUREMENT here and nothing else.

PHASE 1 PERSISTS THE ENERGIES, AND THAT IS THE POINT
====================================================
126 x 500 genuine single points is an asset this project would otherwise re-derive at AMBER cost
every time anyone asks an AMBER question. D1 persisted per-partition SUMMARIES and not energies,
which is exactly why `agentD_FINDINGS.md` §3.3 cannot be strengthened against Lane E's generic
floor without repeating a 541 s hold. Phase 1 writes `E_AMBER` per candidate plus the pool
indices, so any partition or blend can be rebuilt later **without touching OpenMM**.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from s24 import d_harness as H            # noqa: E402
from s24 import stats_lib as ST           # noqa: E402

RES = H.RESULTS
CACHE = os.path.join(HERE, "cache_amber")
os.makedirs(CACHE, exist_ok=True)
SALT = "s24laneD_D2"
K = 500
M = 75

#: the w grid, symmetric so the sign-flip control lives on the same lattice, and dense near 0
#: because that is where the registered expectation puts the optimum.
WGRID = np.array([-1.0, -0.5, -0.25, -0.15, -0.10, -0.06, -0.03, -0.015,
                  0.0,
                  0.015, 0.03, 0.06, 0.10, 0.15, 0.25, 0.5, 1.0])


def zrank(x) -> np.ndarray:
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(float(r.std()), 1e-12)


# ================================================================ PHASE 1: the energies
def phase_energies(pdbs: List[str]) -> int:
    """Genuine AMBER single points for every candidate of every target. Persisted per target.

    Serialised under LOCK_AMBER, announced on open and release. Resumable: a target whose cache
    file already exists is skipped, so an interrupted hold does not cost the whole run again.
    """
    from s13 import qarch_lib as QA
    from s20 import qb2_lib as QB2
    todo = [p for p in pdbs if not os.path.exists(os.path.join(CACHE, f"{p}.npz"))]
    print(f"phase=energies  {len(todo)} of {len(pdbs)} targets to compute "
          f"({len(pdbs) - len(todo)} already cached)", flush=True)
    if not todo:
        return 0
    t0 = time.time()
    with H.AmberLock(tag=f"D2 energies {len(todo)} targets x k={K}"):
        for ii, pdb in enumerate(todo):
            try:
                cand = H.Candidates.from_universe(pdb, k=K)
                rep = QA.Space(pdb, 4).rep
                sp = QB2.AmberSP(cand.seq, rep)
                ver = sp.verify(cand.PHI[:4], cand.PSI[:4], m=4)
                e = sp.batch(cand.PHI, cand.PSI)
                sc = H.score_shipped(cand)
                # NOTE: `np.savez_compressed` APPENDS ".npz" when the path does not already end
                # in it, which silently defeats a tmp+os.replace rename. Write through an open
                # file handle so the name on disk is exactly the one we chose.
                tmp = os.path.join(CACHE, f"{pdb}.tmp{os.getpid()}.npz")
                with open(tmp, "wb") as fh:
                    np.savez_compressed(
                        fh, pdb=pdb, n=cand.n, fold=cand.fold, k=cand.k,
                        e_amber=np.asarray(e, float), score_dist=np.asarray(sc, float),
                        universe_idx=np.asarray(cand.meta["universe_idx"], int),
                        amber_verify_max_rel=float(ver[0] if isinstance(ver, (tuple, list))
                                                   else ver))
                os.replace(tmp, os.path.join(CACHE, f"{pdb}.npz"))
            except Exception as ex:
                print(f"  !! {pdb}: {type(ex).__name__}: {ex}", flush=True)
            el = time.time() - t0
            if (ii + 1) % 5 == 0 or ii == len(todo) - 1:
                print(f"  [{ii+1}/{len(todo)}] {pdb}  {el:.0f}s ({el/(ii+1):.1f}s/target)",
                      flush=True)
    return 0


# ================================================================ PHASE 2: the analysis
def _rmsd_at(cand, blended: np.ndarray, m: int = M) -> float:
    top = np.argsort(blended, kind="stable")[:m]
    C, _ = I.coordinate_average(cand.W[top])
    return float(I.ca_rmsd(C, cand.nat_ca))


def build_curves(pdbs: List[str], permute: bool = False, seed: int = 0) -> Dict:
    """RMSD at every w on the grid, per target. No AMBER, no lock -- reads the cache."""
    rows = []
    for ii, pdb in enumerate(pdbs):
        f = os.path.join(CACHE, f"{pdb}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f, allow_pickle=True)
        cand = H.Candidates.from_universe(pdb, k=int(z["k"]))
        sd, ea = np.asarray(z["score_dist"], float), np.asarray(z["e_amber"], float)
        assert np.allclose(H.score_shipped(cand), sd, rtol=0, atol=0), \
            f"{pdb}: cached distogram score does not reproduce -- cache is stale"
        if permute:
            # NULL 3: preserve AMBER's marginal exactly, destroy its correspondence to the
            # candidates. Any gain that survives this is not about AMBER's ordering.
            ea = ea[SD.stable_rng(pdb, "permnull", int(seed), salt=SALT).permutation(len(ea))]
        za, zd = zrank(ea), zrank(sd)
        rows.append(dict(pdb=pdb, fold=int(z["fold"]), n=int(z["n"]),
                         rmsd=[_rmsd_at(cand, zd + w * za) for w in WGRID]))
        if (ii + 1) % 20 == 0:
            print(f"    curves {ii+1}/{len(pdbs)}", flush=True)
    return dict(rows=rows, wgrid=WGRID.tolist())


def nested_cv(rows: List[Dict]) -> Dict:
    """Nested leave-one-fold-out: `w` is fitted on the OTHER folds and applied to the held-out
    one, so no target is ever scored under a `w` its own fold helped choose.

    The per-fold selected `w` is reported, not just its mean -- picking 0 on some folds and a
    small nonzero elsewhere is a different finding from picking the same small w everywhere,
    and a mean hides that.
    """
    R = np.array([r["rmsd"] for r in rows], float)            # (n_targets, n_w)
    folds = np.array([r["fold"] for r in rows], int)
    uf = np.unique(folds)
    out = np.zeros(len(rows))
    picked = {}
    for f in uf:
        te = folds == f
        tr = ~te
        j = int(np.argmin(R[tr].mean(0)))                     # fit on the OTHER folds only
        picked[int(f)] = float(WGRID[j])
        out[te] = R[te, j]
    j0 = int(np.flatnonzero(WGRID == 0.0)[0])
    return dict(nested=out, w_per_fold=picked, w_zero_col=R[:, j0], R=R, j0=j0,
                n_folds_picking_zero=int(sum(1 for v in picked.values() if v == 0.0)))


def phase_analyse(pdbs: List[str], seed: int = 0) -> int:
    print("phase=analyse  (no AMBER, no lock)", flush=True)
    real = build_curves(pdbs, permute=False)
    rows = real["rows"]
    n = len(rows)
    folds = np.array([r["fold"] for r in rows], int)
    R = np.array([r["rmsd"] for r in rows], float)
    j0 = int(np.flatnonzero(WGRID == 0.0)[0])
    base = R[:, j0]
    print(f"\nn = {n} targets. w=0 mean = {base.mean():.4f} A "
          f"(incumbent 3.0483 on the full 126)")

    print("\n--- the w curve, dev mean RMSD (descriptive; nested CV is the arm) -----------")
    print(f"{'w':>8}{'mean':>10}{'vs w=0':>10}{'W/L':>10}")
    for j, w in enumerate(WGRID):
        d = R[:, j] - base
        print(f"{w:>8.3f}{R[:, j].mean():>10.4f}{d.mean():>+10.4f}"
              f"{int((d<0).sum()):>6}/{int((d>0).sum()):<4}")

    cv = nested_cv(rows)
    print(f"\n--- NESTED LOFO CV -- THE PRIMARY -------------------------------------------")
    print(f"  selected w per fold: {cv['w_per_fold']}")
    print(f"  folds picking w = 0 exactly: {cv['n_folds_picking_zero']} of "
          f"{len(cv['w_per_fold'])}")
    s = H.paired_stats(cv["nested"], base, folds, seed=seed,
                       name_a="blended_nestedCV", name_b="w=0 (incumbent)")
    print(f"  blended - incumbent  {s['mean']:+.4f}  SE {s['se']:.4f}  MDE {s['mde']:.4f}  "
          f"{s['eff_over_mde']:.2f}x")
    print(f"  CI iid [{s['ci_iid'][0]:+.4f},{s['ci_iid'][1]:+.4f}]  "
          f"CI fold [{s['ci_fold'][0]:+.4f},{s['ci_fold'][1]:+.4f}]  "
          f"{s['W']}W/{s['L']}L  worst {s['worst_degradation']:+.4f}")
    print(f"  VERDICT: {s['verdict']}")

    # ---- NULL 2: sign-flipped w, adversarial. Applied per fold at -w_selected.
    flip = np.zeros(n)
    for f, w in cv["w_per_fold"].items():
        jj = int(np.argmin(np.abs(WGRID + w)))
        flip[folds == f] = R[folds == f, jj]
    sf = H.paired_stats(flip, base, folds, seed=seed, name_a="sign_flipped", name_b="w=0")
    print(f"\n--- NULL 2 sign-flipped w (adversarial) ---")
    print(f"  {sf['mean']:+.4f}  SE {sf['se']:.4f}  {sf['eff_over_mde']:.2f}x  "
          f"fold[{sf['ci_fold'][0]:+.4f},{sf['ci_fold'][1]:+.4f}]  {sf['W']}W/{sf['L']}L  "
          f"{sf['verdict']}")

    # ---- NULL 3: rank-permuted AMBER -- marginal preserved, correspondence destroyed.
    print(f"\n--- NULL 3 rank-permuted AMBER (marginal kept, correspondence destroyed) ---")
    perm = build_curves(pdbs, permute=True, seed=seed)
    cvp = nested_cv(perm["rows"])
    Rp = np.array([r["rmsd"] for r in perm["rows"]], float)
    sp = H.paired_stats(cvp["nested"], Rp[:, j0], folds, seed=seed,
                        name_a="permuted_nestedCV", name_b="w=0")
    print(f"  selected w per fold: {cvp['w_per_fold']}")
    print(f"  {sp['mean']:+.4f}  SE {sp['se']:.4f}  {sp['eff_over_mde']:.2f}x  "
          f"fold[{sp['ci_fold'][0]:+.4f},{sp['ci_fold'][1]:+.4f}]  {sp['W']}W/{sp['L']}L  "
          f"{sp['verdict']}")

    # ---- the leakage-optimism figure, reported separately and never as the result
    j_best = int(np.argmin(R.mean(0)))
    opt = float(R[:, j_best].mean() - base.mean())
    print(f"\n--- ORACLE-w optimism (leakage; reported, NEVER the result) ---")
    print(f"  best single w on the FULL dev set: w = {WGRID[j_best]:+.3f}, "
          f"mean {R[:, j_best].mean():.4f}, {opt:+.4f} vs w=0")
    print(f"  nested CV recovers {s['mean']:+.4f} of that -- the gap is the optimism.")

    out = dict(
        n=n, wgrid=WGRID.tolist(), mean_by_w=R.mean(0).tolist(),
        w_zero_mean=float(base.mean()),
        nested_cv=dict(w_per_fold=cv["w_per_fold"],
                       n_folds_picking_zero=cv["n_folds_picking_zero"],
                       mean=float(cv["nested"].mean()), stats=s),
        null_sign_flipped=sf, null_rank_permuted=dict(w_per_fold=cvp["w_per_fold"], stats=sp),
        oracle_w=dict(w=float(WGRID[j_best]), mean=float(R[:, j_best].mean()),
                      optimism_vs_nested=float(opt - s["mean"])),
        rows=[dict(pdb=r["pdb"], fold=r["fold"], n=r["n"], rmsd=r["rmsd"],
                   nested=float(cv["nested"][i]), w0=float(base[i]))
              for i, r in enumerate(rows)],
    )
    ST.save_atomic(os.path.join(RES, "d2_amberscore.json"), out,
                   complete_keys=("pdb", "fold", "rmsd", "nested", "w0"),
                   rows=out["rows"], n_expected=len(pdbs), module_file=__file__)
    print(f"\nwrote results/d2_amberscore.json")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=("energies", "analyse"), required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    pdbs = sorted(t["pdb"] for t in I.targets())
    if a.limit:
        pdbs = pdbs[:a.limit]
    return phase_energies(pdbs) if a.phase == "energies" else phase_analyse(pdbs, a.seed)


if __name__ == "__main__":
    sys.exit(main())
