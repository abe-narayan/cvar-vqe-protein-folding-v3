"""SPRINT 15, INFO, PART C -- the coverage x uncertainty phase diagram.

The question the brief needs answered before Phase 2 designs anything: **given a
target-specific torsion channel of accuracy sigma covering a fraction c of the residues,
what CA-RMSD is reachable?**  Contours at 3.0, 2.5 and 2.0 A, with the regime the project
actually occupies marked on the plane.

CONSTRUCTION (ORACLE DIAGNOSTIC -- this measures a ceiling, it is not an emitter).

    covered residues     native (phi_i, psi_i) + noise from the stated error model
    uncovered residues   filled from the NATIVE-FREE fallback: a random window of the
                         shipped top-75 retrieval pool supplies (phi_i, psi_i).  So
                         c = 0 is the retrieval channel alone, and c = 1 with sigma = 0 is
                         the ideal-geometry build floor.
    structure            `I.build_ca` -- the bit-exact builder the pipeline projects onto
    score                CA-RMSD to the native trace, averaged over 126 targets and R reps

Two facts from sprint 14 are built in rather than rediscovered: all real emitters sit in the
near-i.i.d. band (lag-1 autocorrelation -0.092..+0.066), and torsion-error cost is a
symmetric mid-chain hump, so uniform dropout is the HARDER missingness model and is used for
the primary diagram; contiguous and terminal dropout are run as separate surfaces.

ERROR MODELS
    iid          wrapped normal, sd sigma, independent per torsion
    ar1_0.5/0.8  AR(1) along the chain at the same marginal sd -- correlated error
    resclass     one offset per residue TYPE per chain (systematic, sd sigma) + 0.3 sigma iid
    ssbias       one offset per native SS class (H/E/C) per chain, sd sigma, + 0.3 sigma iid
    outlier      a fraction f of residues get uniform-random torsions, the rest are exact

    python -m s15.info_phase            # primary surface + all error models
    python -m s15.info_phase quick      # coarse grid, for a smoke test
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s15 import info_lib as L            # noqa: E402
import peptide_db as pdb                  # noqa: E402

SIGMAS = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 65.0, 80.0, 110.0]
COVS = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
REPS = 24
ALPHA_HELIX = (np.radians(-63.0), np.radians(-42.0))


# ------------------------------------------------------------------ per-target inputs
def target_pack(t):
    """Native torsions (ORACLE), the native-free top-75 fallback bank, and the native CA."""
    p = pdb.by_pdb(t["pdb"])
    phi = np.asarray(p.phi, float); psi = np.asarray(p.psi, float)
    u = I.load_univ(t["pdb"])
    pi = I.pool_idx(u)
    rec = I.shipped_record(t["pdb"])
    sub = np.asarray(rec["sub"], int)
    q = pi[sub]
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "seq": t["seq"],
            "phi": phi, "psi": psi, "nat": np.asarray(u["nat_ca"], float),
            "bphi": np.asarray(u["PHI"], float)[q], "bpsi": np.asarray(u["PSI"], float)[q],
            "ss": I.ss_of(phi, psi)}


# ------------------------------------------------------------------ error models
def noise(model, sigma_deg, n, R, rng, pack):
    """(R, n) phi and psi angular errors in radians."""
    s = np.radians(sigma_deg)
    if s == 0.0 and model != "outlier":
        return np.zeros((R, n)), np.zeros((R, n))
    if model == "iid":
        return rng.normal(0, s, (R, n)), rng.normal(0, s, (R, n))
    if model.startswith("ar1"):
        rho = float(model.split("_")[1])
        out = []
        for _ in range(2):
            e = np.empty((R, n))
            e[:, 0] = rng.normal(0, s, R)
            innov = s * np.sqrt(1 - rho ** 2)
            for i in range(1, n):
                e[:, i] = rho * e[:, i - 1] + rng.normal(0, innov, R)
            out.append(e)
        return out[0], out[1]
    if model == "resclass":
        codes = np.asarray([ord(c) for c in pack["seq"]])
        uc = {c: k for k, c in enumerate(sorted(set(codes.tolist())))}
        idx = np.asarray([uc[c] for c in codes.tolist()])
        bp = rng.normal(0, s, (R, len(uc)))[:, idx]
        bs = rng.normal(0, s, (R, len(uc)))[:, idx]
        return bp + rng.normal(0, 0.3 * s, (R, n)), bs + rng.normal(0, 0.3 * s, (R, n))
    if model == "ssbias":
        ss = pack["ss"]
        cls = {c: k for k, c in enumerate("HEC")}
        idx = np.asarray([cls.get(c, 2) for c in ss[:n]])
        bp = rng.normal(0, s, (R, 3))[:, idx]
        bs = rng.normal(0, s, (R, 3))[:, idx]
        return bp + rng.normal(0, 0.3 * s, (R, n)), bs + rng.normal(0, 0.3 * s, (R, n))
    raise ValueError(model)


def dropout_mask(kind, c, n, R, rng):
    """(R, n) True where the residue IS covered by the target-specific channel."""
    k = int(round(c * n))
    M = np.zeros((R, n), bool)
    if k <= 0:
        return M
    if k >= n:
        return ~M
    miss = n - k
    if kind == "uniform":
        for r in range(R):
            M[r, rng.choice(n, k, replace=False)] = True
    elif kind == "terminal":
        # the CHEAP missingness: the covered residues form one contiguous block, so the
        # gaps sit at the two chain ends.  (An earlier "contig" model was identical to
        # this one and has been removed rather than double-counted.)
        for r in range(R):
            a = rng.integers(0, miss + 1)
            M[r, a:a + k] = True
    elif kind == "contig_gap":
        # the missing residues form ONE contiguous interior block at a random position
        for r in range(R):
            a = rng.integers(0, k + 1)
            M[r] = True
            M[r, a:a + miss] = False
    elif kind == "mid_gap":
        # the missing residues form one contiguous block CENTRED on the chain -- the
        # worst case implied by the mid-chain error hump
        a = max(0, (n - miss) // 2)
        M[:] = True
        M[:, a:a + miss] = False
    else:
        raise ValueError(kind)
    return M


def fallback(kind, pack, n, R, rng):
    if kind == "pool":
        j = rng.integers(0, len(pack["bphi"]), (R, n))
        rows = np.arange(n)[None, :]
        return pack["bphi"][j, rows], pack["bpsi"][j, rows]
    if kind == "helix":
        return (np.full((R, n), ALPHA_HELIX[0]), np.full((R, n), ALPHA_HELIX[1]))
    if kind == "uniform":
        return (rng.uniform(-np.pi, np.pi, (R, n)), rng.uniform(-np.pi, np.pi, (R, n)))
    raise ValueError(kind)


def cell(pack, sigma, c, rng, model="iid", drop="uniform", fb="pool", R=REPS,
         outlier_frac=None, rng_cov=None):
    """COMMON RANDOM NUMBERS: the dropout mask and the fallback fill are drawn from
    `rng_cov`, which depends only on (target, coverage), so every sigma row of a column
    sees the SAME missing residues and the SAME fill.  Without this the coverage=0 row
    is not constant across sigma and the contours pick up pure Monte-Carlo noise."""
    n = pack["n"]
    rc = rng_cov if rng_cov is not None else rng
    M = dropout_mask(drop, c, n, R, rc)
    fphi, fpsi = fallback(fb, pack, n, R, rc)
    if outlier_frac is None:
        ep, es = noise(model, sigma, n, R, rng, pack)
        phi = pack["phi"][None, :] + ep
        psi = pack["psi"][None, :] + es
    else:
        phi = np.tile(pack["phi"], (R, 1)); psi = np.tile(pack["psi"], (R, 1))
        bad = rng.random((R, n)) < outlier_frac
        phi[bad] = rng.uniform(-np.pi, np.pi, int(bad.sum()))
        psi[bad] = rng.uniform(-np.pi, np.pi, int(bad.sum()))
    phi = np.where(M, phi, fphi); psi = np.where(M, psi, fpsi)
    ca = I.build_ca(phi, psi)
    return I.kabsch_rmsd_batch(ca, pack["nat"])


def surface(packs, sigmas, covs, model="iid", drop="uniform", fb="pool", R=REPS, seed=0,
            outlier_fracs=None, log=None):
    """(len(sigmas), len(covs)) mean RMSD over targets, plus the per-target cube."""
    xs = outlier_fracs if outlier_fracs is not None else sigmas
    cube = np.empty((len(packs), len(xs), len(covs)))
    t0 = time.time()
    for a, pk in enumerate(packs):
        for b, x in enumerate(xs):
            rng = np.random.default_rng(seed + 1009 * a + 7919 * b)
            for d, c in enumerate(covs):
                rc = np.random.default_rng(seed + 1009 * a + 104729 * d)
                if outlier_fracs is None:
                    v = cell(pk, x, c, rng, model, drop, fb, R, rng_cov=rc)
                else:
                    v = cell(pk, 0.0, c, rng, model, drop, fb, R, outlier_frac=x,
                             rng_cov=rc)
                cube[a, b, d] = v.mean()
        if log and (a % 25 == 0 or a == len(packs) - 1):
            print(f"    {log}: {a + 1}/{len(packs)} [{time.time() - t0:.0f}s]", flush=True)
    return cube.mean(0), cube


def contour_sigma(grid, sigmas, covs, level):
    """For each coverage column, the largest sigma whose mean RMSD is <= level
    (linear interpolation between grid rows).  NaN if unreachable at sigma=0."""
    out = []
    for d in range(len(covs)):
        col = grid[:, d]
        if col[0] > level:
            out.append(float("nan")); continue
        k = np.flatnonzero(col > level)
        if len(k) == 0:
            out.append(float(sigmas[-1])); continue
        k = k[0]
        s0, s1 = sigmas[k - 1], sigmas[k]
        v0, v1 = col[k - 1], col[k]
        out.append(float(s0 + (level - v0) * (s1 - s0) / (v1 - v0)))
    return out


def main(quick=False):
    tg = I.targets()
    sig = SIGMAS if not quick else [0.0, 15.0, 30.0, 60.0]
    cov = COVS if not quick else [0.0, 0.5, 1.0]
    reps = REPS if not quick else 6
    if quick:
        tg = tg[:12]
    print(f"building target packs (n={len(tg)}) ...", flush=True)
    packs = [target_pack(t) for t in tg]
    res = {"sigmas": sig, "covs": cov, "reps": reps, "n_targets": len(packs),
           "targets": [p["pdb"] for p in packs], "surfaces": {}}

    print("PRIMARY SURFACE: i.i.d. angular noise, uniform dropout, top-75 pool fallback",
          flush=True)
    g, cube = surface(packs, sig, cov, "iid", "uniform", "pool", reps, log="primary")
    res["surfaces"]["primary_iid_uniform_pool"] = {
        "grid": g.tolist(),
        "contours": {str(lv): contour_sigma(g, sig, cov, lv) for lv in (3.0, 2.5, 2.0)}}
    np.savez_compressed(os.path.join(L.CACHE, "phase_primary.npz"), cube=cube,
                        sigmas=sig, covs=cov, targets=[p["pdb"] for p in packs])

    variants = [("ar1_0.5", "ar1_0.5", "uniform", "pool"),
                ("ar1_0.8", "ar1_0.8", "uniform", "pool"),
                ("resclass", "resclass", "uniform", "pool"),
                ("ssbias", "ssbias", "uniform", "pool"),
                ("iid_terminal_gap", "iid", "terminal", "pool"),
                ("iid_contig_gap", "iid", "contig_gap", "pool"),
                ("iid_mid_gap", "iid", "mid_gap", "pool"),
                ("iid_uniform_helixfill", "iid", "uniform", "helix"),
                ("iid_uniform_randomfill", "iid", "uniform", "uniform")]
    for nm, md, dr, fb in variants:
        print(f"SURFACE: {nm}", flush=True)
        gg, _ = surface(packs, sig, cov, md, dr, fb, reps, seed=7, log=nm)
        res["surfaces"][nm] = {
            "grid": gg.tolist(),
            "contours": {str(lv): contour_sigma(gg, sig, cov, lv)
                         for lv in (3.0, 2.5, 2.0)}}

    fr = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6]
    print("SURFACE: sparse gross outliers", flush=True)
    go, _ = surface(packs, sig, cov, "iid", "uniform", "pool", reps, seed=11,
                    outlier_fracs=fr, log="outlier")
    res["surfaces"]["outlier"] = {"fracs": fr, "grid": go.tolist()}

    L.jwrite("info_phase", res)
    report(res)
    return res


def report(res):
    sig, cov = res["sigmas"], res["covs"]
    g = np.asarray(res["surfaces"]["primary_iid_uniform_pool"]["grid"])
    print("\nPART C -- COVERAGE x UNCERTAINTY -> RMSD  "
          f"(mean over {res['n_targets']} targets, i.i.d. noise, uniform dropout, "
          "top-75 pool fallback)")
    print("rows = per-torsion sigma (deg), cols = coverage\n")
    print("sigma\\cov " + "".join(f"{c:>7.2f}" for c in cov))
    for a, s in enumerate(sig):
        print(f"{s:>8.0f}  " + "".join(f"{g[a, d]:>7.3f}" for d in range(len(cov))))
    print("\nCONTOURS -- the largest sigma (deg) that still reaches the level, per coverage")
    print("level\\cov " + "".join(f"{c:>7.2f}" for c in cov))
    for lv in (3.0, 2.5, 2.0):
        row = res["surfaces"]["primary_iid_uniform_pool"]["contours"][str(lv)]
        print(f"{lv:>8.1f}  " + "".join(
            ("    n/a" if np.isnan(v) else f"{v:>7.1f}") for v in row))
    print("\nERROR-MODEL COMPARISON at coverage 1.0 (sigma required for each level)")
    print(f"{'surface':<26}{'3.0 A':>9}{'2.5 A':>9}{'2.0 A':>9}")
    for nm, s in res["surfaces"].items():
        if "contours" not in s:
            continue
        v = [s["contours"][str(lv)][-1] for lv in (3.0, 2.5, 2.0)]
        print(f"{nm:<26}" + "".join(("      n/a" if np.isnan(x) else f"{x:>9.1f}")
                                    for x in v))
    print("\nSPARSE GROSS OUTLIERS (fraction of residues with a random torsion, "
          "the rest EXACT)")
    o = res["surfaces"]["outlier"]
    print("frac\\cov  " + "".join(f"{c:>7.2f}" for c in cov))
    for a, f in enumerate(o["fracs"]):
        print(f"{f:>8.2f}  " + "".join(f"{np.asarray(o['grid'])[a, d]:>7.3f}"
                                       for d in range(len(cov))))


if __name__ == "__main__":
    main(quick=(len(sys.argv) > 1 and sys.argv[1] == "quick"))
