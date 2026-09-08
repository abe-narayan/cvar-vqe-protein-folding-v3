"""SPRINT 14, OBJ, STEP 2 -- the discrimination floor of the existing objectives.

Measured on the fully enumerated k=4 targets, so every number is exact over the whole
state space, not a sample estimate.  For Legacy, its 11 components, the 1-local prior,
AMBER and its 5 components:

  * Spearman rho with true CA-RMSD, globally and inside the objective's own low-energy
    decile (the decile is where a search actually lives)
  * percentile of the ORACLE-snap configuration in the objective's ranking
  * top-k enrichment: mean RMSD of the k lowest-energy configurations vs a random k
  * PAIRWISE DISCRIMINATION AS A FUNCTION OF STRUCTURAL SEPARATION -- P(objective orders
    the pair correctly | the pair is dRMSD apart).  0.5 is chance.  This says at what
    structural distance each objective stops being able to tell two structures apart.
  * basin structure: mean RMSD in each energy decile.

RMSD is used here ONLY as a post-hoc evaluation label.  ORACLE DIAGNOSTIC throughout;
none of this is an inference-time claim.

    python -m s14.obj_floor            # every enumerated target
    python -m s14.obj_floor 1CS9 2MK7
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

from s13.qarch_lib import spearman            # noqa: E402
from s12 import instrument as I               # noqa: E402
from s14 import obj_enum as E                 # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

DBINS = [(0.25, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0),
         (2.0, 3.0), (3.0, 4.0), (4.0, 5.0), (5.0, 99.0)]
NPAIR = 400_000


def pair_discrimination(e, r, rng, npair=NPAIR, sub=None):
    """P(e orders a pair the same way rmsd does), binned by |dRMSD|.

    Ties in e count as 0.5.  Returns {bin_label: (rate, n_pairs)}.
    """
    idx = np.arange(len(e)) if sub is None else np.asarray(sub)
    if len(idx) < 100:
        return {}
    i = rng.choice(idx, npair)
    j = rng.choice(idx, npair)
    ok = i != j
    i, j = i[ok], j[ok]
    dr = r[i] - r[j]
    de = e[i] - e[j]
    out = {}
    ad = np.abs(dr)
    for lo, hi in DBINS:
        m = (ad >= lo) & (ad < hi)
        nm = int(m.sum())
        if nm < 200:
            out[f"{lo}-{hi}"] = (float("nan"), nm)
            continue
        agree = np.sign(dr[m]) == np.sign(de[m])
        tie = de[m] == 0
        rate = float((agree & ~tie).sum() + 0.5 * tie.sum()) / nm
        out[f"{lo}-{hi}"] = (rate, nm)
    return out


def analyse_column(e, r, snap_idx, rng, tag):
    """All Step-2 statistics for one energy column `e` over configurations with RMSD `r`."""
    B = len(e)
    fin = np.isfinite(e)
    if fin.sum() < 100:
        return None
    order = np.argsort(e, kind="mergesort")
    ndec = max(50, B // 10)
    dec = order[:ndec]                                   # low-energy decile
    res = dict(
        n=int(B),
        rho_global=spearman(e, r),
        rho_decile=spearman(e[dec], r[dec]),
        rmsd_argmin=float(r[order[0]]),
        rmsd_mean=float(r.mean()),
        rmsd_min=float(r.min()),
        rmsd_decile_mean=float(r[dec].mean()),
        rmsd_decile_min=float(r[dec].min()),
    )
    for kk in (1, 10, 100, 1000):
        if kk <= B:
            res[f"top{kk}_mean_rmsd"] = float(r[order[:kk]].mean())
            res[f"top{kk}_min_rmsd"] = float(r[order[:kk]].min())
    if snap_idx is not None and snap_idx < B:
        res["snap_pct"] = float((e < e[snap_idx]).mean())
        res["snap_rmsd"] = float(r[snap_idx])
    # percentile of the best-RMSD configuration in the energy ranking
    best = int(np.argmin(r))
    res["best_rmsd_pct"] = float((e < e[best]).mean())
    # basin structure: mean rmsd by energy decile
    res["decile_rmsd"] = [float(r[order[a * B // 10:(a + 1) * B // 10]].mean())
                          for a in range(10)]
    res["pair_global"] = pair_discrimination(e, r, rng)
    res["pair_inband"] = pair_discrimination(e, r, rng, sub=dec)
    return res


def analyse_target(pdb_id, seed=0):
    z = E.load_enum(pdb_id)
    r = np.asarray(z["rmsd"], np.float64)
    snap = int(z["snap_index"])
    rng = np.random.default_rng(seed)
    cols = {"legacy": np.asarray(z["legacy"], np.float64),
            "prior": np.asarray(z["prior"], np.float64)}
    for f in z.files:
        if f.startswith("leg_"):
            cols[f] = np.asarray(z[f], np.float64)
    out = dict(pdb=str(z["pdb"]), n=int(z["n"]), fold=int(z["fold"]),
               seq=str(z["seq"]), configs=int(len(r)),
               rmsd_mean=float(r.mean()), rmsd_min=float(r.min()),
               snap_rmsd=float(r[snap]), full={}, amber={})
    for name, e in cols.items():
        out["full"][name] = analyse_column(e, r, snap, rng, name)

    # --- AMBER lives on the labelled subsample only -----------------------
    ai = np.asarray(z["amber_idx"], np.int64)
    ra = r[ai]
    sub_snap = int(np.where(ai == snap)[0][0]) if (ai == snap).any() else None
    acols = {"amber_total": np.asarray(z["amber_total"], np.float64)}
    for f in z.files:
        if f.startswith("amb_"):
            acols[f] = np.asarray(z[f], np.float64)
    # Legacy and the prior re-measured on the SAME subsample, so AMBER is comparable
    acols["legacy@sub"] = np.asarray(z["legacy"], np.float64)[ai]
    acols["prior@sub"] = np.asarray(z["prior"], np.float64)[ai]
    for name, e in acols.items():
        out["amber"][name] = analyse_column(e, ra, sub_snap, rng, name)
    return out


def pool(rows, section, name, field):
    v = [r[section][name][field] for r in rows
         if r[section].get(name) and np.isfinite(r[section][name].get(field, np.nan))]
    return (float(np.mean(v)), len(v)) if v else (float("nan"), 0)


def main(args):
    pdbs = [a for a in args if not a.startswith("-")]
    if not pdbs:
        pdbs = [t["pdb"] for t in I.targets() if E.have(t["pdb"])]
    rows = []
    t0 = time.time()
    for p in pdbs:
        rows.append(analyse_target(p))
        print(f"{p}: legacy rho_g={rows[-1]['full']['legacy']['rho_global']:+.3f} "
              f"rho_dec={rows[-1]['full']['legacy']['rho_decile']:+.3f} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    with open(os.path.join(RESULTS, "obj_floor.json"), "w") as fh:
        json.dump({"what": "Step 2 discrimination floor on enumerated k=4 targets",
                   "n_expected": len(pdbs), "complete": len(rows) == len(pdbs),
                   "per_target": rows}, fh, indent=1, default=str)

    print("\n== POOLED (mean over targets) ==")
    print(f"{'objective':22s} {'rho_glob':>9s} {'rho_dec':>9s} {'snap_pct':>9s} "
          f"{'top1':>7s} {'top100':>7s} {'dec_mean':>8s} {'all_mean':>8s}")
    names = ["legacy", "prior"] + sorted(
        {k for r in rows for k in r["full"] if k.startswith("leg_")})
    for nm in names:
        rg = pool(rows, "full", nm, "rho_global")[0]
        rd = pool(rows, "full", nm, "rho_decile")[0]
        sp = pool(rows, "full", nm, "snap_pct")[0]
        t1 = pool(rows, "full", nm, "top1_mean_rmsd")[0]
        t100 = pool(rows, "full", nm, "top100_mean_rmsd")[0]
        dm = pool(rows, "full", nm, "rmsd_decile_mean")[0]
        am = float(np.mean([r["rmsd_mean"] for r in rows]))
        print(f"{nm:22s} {rg:+9.3f} {rd:+9.3f} {sp:9.3f} {t1:7.3f} {t100:7.3f} "
              f"{dm:8.3f} {am:8.3f}")
    print("\n== AMBER SUBSAMPLE ==")
    anames = ["amber_total", "legacy@sub", "prior@sub"] + sorted(
        {k for r in rows for k in r["amber"] if k.startswith("amb_")})
    for nm in anames:
        rg = pool(rows, "amber", nm, "rho_global")[0]
        rd = pool(rows, "amber", nm, "rho_decile")[0]
        sp = pool(rows, "amber", nm, "snap_pct")[0]
        t1 = pool(rows, "amber", nm, "top1_mean_rmsd")[0]
        t100 = pool(rows, "amber", nm, "top100_mean_rmsd")[0]
        print(f"{nm:22s} {rg:+9.3f} {rd:+9.3f} {sp:9.3f} {t1:7.3f} {t100:7.3f}")

    print("\n== PAIRWISE DISCRIMINATION vs STRUCTURAL SEPARATION (pooled) ==")
    print("P(objective orders the pair correctly); 0.5 = chance")
    hdr = "  ".join(f"{lo}-{hi}" for lo, hi in DBINS)
    for scope in ("pair_global", "pair_inband"):
        print(f"-- {scope}")
        print(f"{'objective':22s} {hdr}")
        for nm in ["legacy", "prior", "leg_steric", "leg_contact", "leg_compactness"]:
            vals = []
            for lo, hi in DBINS:
                key = f"{lo}-{hi}"
                v = [r["full"][nm][scope][key][0] for r in rows
                     if r["full"].get(nm) and key in r["full"][nm][scope]
                     and np.isfinite(r["full"][nm][scope][key][0])]
                vals.append(np.mean(v) if v else float("nan"))
            print(f"{nm:22s} " + "  ".join(f"{v:7.3f}" for v in vals))
        for nm in ["amber_total"]:
            vals = []
            for lo, hi in DBINS:
                key = f"{lo}-{hi}"
                v = [r["amber"][nm][scope][key][0] for r in rows
                     if r["amber"].get(nm) and key in r["amber"][nm][scope]
                     and np.isfinite(r["amber"][nm][scope][key][0])]
                vals.append(np.mean(v) if v else float("nan"))
            print(f"{nm:22s} " + "  ".join(f"{v:7.3f}" for v in vals))
    return rows


if __name__ == "__main__":
    main(sys.argv[1:])
