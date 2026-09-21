#!/usr/bin/env python
"""s32/s32_D1_inband.py -- S32 lane D, rungs D1-E / D1-L / D0-X.  ALL FREE OR NEARLY FREE.

Registered in s32/PREREG_S32_D.md, commit 34973b1b.

D1-E  AMBER ff14SB/GBn2 single-point energy: IN-BAND Spearman rho against true CA-RMSD,
      inside each target's shipped top-75 band, n=126.  Charter section 32 says "do not
      assume AMBER is a good ranker" -- this measures it IN BAND on the shipped instrument,
      which the record does not contain.  Free: 126 x 500 genuine single points are already
      on disk at s24/cache_amber/<PDB>.npz.

D1-L  Legacy, all 11 terms + the weighted total, same basis, same band.  Charter section 33.

D0-X  The chiral degree of freedom's IN-BAND VARIANCE SHARE.  A channel can only carry
      in-band information about a degree of freedom that VARIES in band.  If the chiral
      coordinate is near-constant across a target's top-75 members, no chiral functional can
      have in-band skill -- which would convert S30/S31's "empty at 9-16 residues" scope note
      into a mechanism.  Registered prediction: in-band share < 15%.

Every rho ships with its Rg-PARTIALLED twin in the same row (memory
`physics-ranks-real-geometry-not-lattice`: AMBER's nonbonded+solvation reward compactness, so
any physics-vs-RMSD correlation runs through Rg unless Rg is partialled).

ORACLE use: u["rr"] (per-candidate CA-RMSD to native) is read to SCORE the rankers.  That is
an evaluation label, not an input to any ranker -- no arm here reads it.

    python s32/s32_D1_inband.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I          # noqa: E402
from s24 import stats_lib as ST          # noqa: E402
from s16 import energy_lib as EL         # noqa: E402
from core import energy as CE            # noqa: E402

OUT = os.path.join(HERE, "results", "s32_D1_inband.json")
CACHE_AMBER = os.path.join(ROOT, "s24", "cache_amber")


# ----------------------------------------------------------------------------- rank tools
def _rank(a):
    """Average-rank transform (ties averaged).  Spearman = Pearson of these."""
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    r[order] = np.arange(len(a), dtype=float)
    # average ties
    s = a[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = np.mean(r[order[i:j + 1]])
        i = j + 1
    return r


def _pear(x, y):
    x = x - x.mean()
    y = y - y.mean()
    dn = np.sqrt((x * x).sum() * (y * y).sum())
    return float(x.dot(y) / dn) if dn > 1e-300 else np.nan


def spear(a, b):
    return _pear(_rank(a), _rank(b))


def spear_partial(a, b, c):
    """Spearman(a,b) with c partialled out, on ranks (residualise then correlate)."""
    ra, rb, rc = _rank(a), _rank(b), _rank(c)
    rc0 = rc - rc.mean()
    den = float(rc0.dot(rc0))
    if den < 1e-300:
        return spear(a, b)
    ea = ra - ra.mean() - rc0 * (rc0.dot(ra - ra.mean()) / den)
    eb = rb - rb.mean() - rc0 * (rc0.dot(rb - rb.mean()) / den)
    return _pear(ea, eb)


def rg(W):
    """Radius of gyration per structure, W (b,n,3)."""
    C = W - W.mean(1, keepdims=True)
    return np.sqrt((C ** 2).sum(-1).mean(-1))


def ca_pseudo_torsion(W):
    """(b, n-3) dihedral over 4 consecutive CA atoms, radians.  CHIRAL: a reflection of W
    sends every value to its negative.  This is the cheapest faithful chiral coordinate of a
    CA trace."""
    b0 = W[:, 1:-2] - W[:, 0:-3]
    b1 = W[:, 2:-1] - W[:, 1:-2]
    b2 = W[:, 3:] - W[:, 2:-1]
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.linalg.norm(b1, axis=-1, keepdims=True))
    x = (n1 * n2).sum(-1)
    y = (m1 * n2).sum(-1)
    return np.arctan2(y, x)


# ----------------------------------------------------------------------------- main
def main():
    tg = I.targets()
    assert len(tg) == 126, len(tg)
    rows = []
    legacy_terms = list(CE.DEFAULT_WEIGHTS)

    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        p = I.pool_idx(u, k=500)
        rec = I.shipped_record(pdb)
        sub = np.asarray(rec["sub"], int)              # indices INTO the 500-pool
        W = np.asarray(u["W"][p], float)               # (500, n, 3) CA point cloud
        PHI = np.asarray(u["PHI"][p], float)
        PSI = np.asarray(u["PSI"][p], float)
        rr = np.asarray(u["rr"][p], float)             # ORACLE label, evaluation only
        G = rg(W)

        # ---- AMBER, free from the cache
        z = np.load(os.path.join(CACHE_AMBER, pdb + ".npz"))
        assert int(z["k"]) == 500
        ui = np.asarray(z["universe_idx"], int)
        assert np.array_equal(ui, p), "cache_amber pool order != I.pool_idx on %s" % pdb
        e_amb = np.asarray(z["e_amber"], float)
        sc_dis = np.asarray(z["score_dist"], float)    # the shipped cost, for the S31 bar

        # ---- Legacy, batched ideal-geometry rebuild
        comp = EL.legacy_components_of_windows(t["seq"], PHI, PSI)
        leg_tot = EL.legacy_total_from(comp)

        # ---- chirality: the CA pseudo-torsion profile
        tau = ca_pseudo_torsion(W)                     # (500, n-3)
        chi_scalar = np.sin(tau).mean(1)               # a chiral summary; odd under reflection

        r = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"])}
        for band, idx in (("b500", np.arange(500)), ("b75", sub)):
            rrb, Gb = rr[idx], G[idx]
            ch = {}
            ch["AMBER"] = e_amb[idx]
            ch["DIS"] = sc_dis[idx]
            ch["LEG_total"] = leg_tot[idx]
            ch["RG"] = Gb
            for tm in legacy_terms:
                ch["LEG_" + tm] = np.asarray(comp[tm], float)[idx]
            for name, v in ch.items():
                if np.std(v) < 1e-12 or np.std(rrb) < 1e-12:
                    r["%s_%s_rho" % (band, name)] = np.nan
                    r["%s_%s_rho_rgpart" % (band, name)] = np.nan
                    continue
                r["%s_%s_rho" % (band, name)] = spear(v, rrb)
                r["%s_%s_rho_rgpart" % (band, name)] = (
                    np.nan if name == "RG" else spear_partial(v, rrb, Gb))
            # ---- D0-X: chiral variance share
            r["%s_chi_sd" % band] = float(np.std(chi_scalar[idx]))
            r["%s_chi_mean" % band] = float(np.mean(chi_scalar[idx]))
            r["%s_tau_sd_mean" % band] = float(np.std(tau[idx], axis=0).mean())
            r["%s_chi_rho" % band] = spear(chi_scalar[idx], rrb)
            r["%s_frac_pos_tau" % band] = float((tau[idx] > 0).mean())
        rows.append(r)
        if len(rows) % 20 == 0:
            print("  %3d/126  %s" % (len(rows), pdb), flush=True)

    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    zero = np.zeros(len(rows))

    out = {"note": ("IN-BAND Spearman rho of each scorer against true CA-RMSD, per target, "
                    "then aggregated over 126 targets.  POSITIVE rho = the scorer orders "
                    "candidates in the correct direction.  stats_lib.compare is LOWER IS "
                    "BETTER and its BETTER/WORSE verdict labels are therefore INVERTED for "
                    "this quantity and are NOT quoted; only effect, se, mde, x_mde, fold CI "
                    "and folds_same_sign are read.  ORACLE rr is an evaluation label only."),
           "prereg_commit": "34973b1b", "n": len(rows), "rows": rows, "agg": {}}

    names = ["AMBER", "DIS", "LEG_total", "RG"] + ["LEG_" + x for x in legacy_terms]
    for band in ("b500", "b75"):
        for name in names:
            for suf in ("rho", "rho_rgpart"):
                key = "%s_%s_%s" % (band, name, suf)
                v = np.array([r[key] for r in rows], float)
                if not np.isfinite(v).all():
                    continue
                c = ST.compare(v, zero, folds=folds, names=pdbs, label=key)
                out["agg"][key] = {"mean": c["effect"], "median": c["median_effect"],
                                   "se": c["se"], "mde": c["mde"],
                                   "x_mde": abs(c["effect"]) / c["mde"] if c["mde"] else np.nan,
                                   "ci95_fold": c["ci95_fold"],
                                   "folds_same_sign": c["folds_same_sign"],
                                   "n_pos": int((v > 0).sum()), "n_neg": int((v < 0).sum())}
        # D0-X aggregate
        chi_in = np.array([r["%s_chi_sd" % band] for r in rows], float)
        out["agg"]["%s_chi_sd_mean" % band] = float(chi_in.mean())
        out["agg"]["%s_chi_mean_mean" % band] = float(
            np.mean([r["%s_chi_mean" % band] for r in rows]))
        out["agg"]["%s_frac_pos_tau" % band] = float(
            np.mean([r["%s_frac_pos_tau" % band] for r in rows]))

    # D0-X: the variance decomposition that decides door 1.
    for band in ("b500", "b75"):
        within = np.array([r["%s_chi_sd" % band] for r in rows], float) ** 2
        means = np.array([r["%s_chi_mean" % band] for r in rows], float)
        between = float(np.var(means))
        out["agg"]["D0X_%s_within_var_mean" % band] = float(within.mean())
        out["agg"]["D0X_%s_between_var" % band] = between
        out["agg"]["D0X_%s_inband_share" % band] = float(
            within.mean() / (within.mean() + between)) if (within.mean() + between) > 0 else np.nan

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    ST.save_atomic(OUT, out, n_expected=126, module_file=__file__)
    print(json.dumps({k: v for k, v in out["agg"].items() if k.startswith("D0X")}, indent=2))
    print("\nIN-BAND (top-75 band), mean rho over 126 targets  [positive = correct direction]")
    print("%-22s %9s %9s %7s %6s %8s %s" % ("scorer", "rho", "rho|Rg", "se", "xMDE", "folds", "pos/neg"))
    for name in names:
        a = out["agg"].get("b75_%s_rho" % name)
        b = out["agg"].get("b75_%s_rho_rgpart" % name)
        if a is None:
            continue
        print("%-22s %+9.4f %+9.4f %7.4f %6.2f %8s %d/%d" % (
            name, a["mean"], (b["mean"] if b else float("nan")), a["se"],
            a["x_mde"], a["folds_same_sign"], a["n_pos"], a["n_neg"]))
    print("\nWHOLE POOL (K=500) for contrast")
    for name in names:
        a = out["agg"].get("b500_%s_rho" % name)
        b = out["agg"].get("b500_%s_rho_rgpart" % name)
        if a is None:
            continue
        print("%-22s %+9.4f %+9.4f %7.4f %6.2f %8s" % (
            name, a["mean"], (b["mean"] if b else float("nan")), a["se"], a["x_mde"],
            a["folds_same_sign"]))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
