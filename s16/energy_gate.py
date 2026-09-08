"""SPRINT 16, ENERGY -- THE FOUR FIXES.  Nothing downstream is quotable until these run.

1. THE CONVERGENCE GATE.  Declared in `core.amber.CONVERGE_MAX_KCAL` / `convergence_flags`
   BEFORE it was applied, in the code itself:

       a restrained minimisation is CONVERGED iff its final potential energy, with the
       restraint switched off, is finite and <= 1000.0 kcal/mol.

   One threshold.  No target-specific tuning.  No native-derived quantity.  The same
   1000 kcal/mol the codebase already uses for its bond+angle strain gate.  It is
   REPORTED, never silently applied: `refine_coords` still returns the structure and now
   sets `converged` / `converge_reason`, and every consumer must print how many targets
   the gate excludes and the effect BOTH gated and ungated.  `energy_lib.gated_paired`
   makes it impossible to report one without the other.

2. THE STANDING FRAME-INVARIANCE NULL.  ff14SB, GBn2, the positional restraint and every
   RMSD in this project are rigid-invariant, so relaxing the SAME structure in a rotated
   lab frame is ZERO BY CONSTRUCTION.  `energy_lib.frame_null` is the permanent
   regression check with a pre-declared PASS band (|mean| <= 0.005 A, max <= 0.05 A on
   the converged subset); `tests/test_amber_frame_invariance.py` runs it on a fixed
   12-target slice.

3. MONOTONE CONDITIONING.  Only order-preserving maps may touch these energies before a
   spectral or moment-based claim.  Winsorising at the 99th percentile is not
   conditioning when the 99th percentile is itself 1e12-1e15, and this module measures
   exactly how far short it falls.

4. THE BINDING AMBER DATA RULE, audited across all 19 enumerated caches:
   `amber_kind == 0 AND amber_idx != snap_index`, per-target n printed.
   `s14/results/obj_floor.json` is RETRACTED and is not opened by this module or by any
   other module in `s16/`.

    python -m s16.energy_gate                # the cached-artefact analysis (no AMBER)
    python -m s16.energy_gate --frame N      # run the standing null live on N targets
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import argparse
import json
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s16 import energy_lib as L            # noqa: E402

S15 = os.path.join(ROOT, "s15", "results")
K = 30.0
FORBIDDEN = os.path.join(ROOT, "s14", "results", "obj_floor.json")


def _rows(name):
    p = os.path.join(S15, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)["per_target"]


def _arms(rows, k=K):
    key = f"amber_k{k}"
    rows = [r for r in rows if key in r and "rmsd" in r[key]]
    return (rows,
            np.array([r["projected_rmsd"] for r in rows]),
            np.array([r[key]["rmsd"] for r in rows]),
            np.array([r[key]["energy"] for r in rows]),
            np.array([r["fold"] for r in rows]),
            [r["pdb"] for r in rows])


# ------------------------------------------------------------------ 1 + 2
def gate_and_null():
    """The k = 30 accuracy effect and the exact frame null, gated and ungated."""
    assert not any("obj_floor" in f for f in os.listdir(os.path.join(ROOT, "s16"))), \
        "obj_floor.json is retracted"
    out = {"gate_rule": {
        "declared_in": "core.amber.CONVERGE_MAX_KCAL / convergence_flags",
        "rule": "final potential energy (restraint off) finite and <= 1000.0 kcal/mol",
        "applied_to": "every AMBER minimisation in s16; reported, never silent"}}

    canon = _rows("phys_repl_canon0.json")
    if canon is None:
        return out
    rows, A, B, E, fold, pdbs = _arms(canon)
    out["effect_k30"] = L.gated_paired(B, A, E, folds=fold, names=pdbs)

    # -- the EXACT NULL: same structure, rotated lab frame.  Zero by construction.
    nulls = {}
    for f in (1, 2, 3):
        fr = _rows(f"phys_repl_frame{f}.json")
        if fr is None:
            continue
        key = f"amber_k{K}"
        cm = {r["pdb"]: r for r in canon}
        pairs = [(r["pdb"], r[key]["rmsd"], cm[r["pdb"]][key]["rmsd"],
                  r[key]["energy"], cm[r["pdb"]][key]["energy"], r["fold"])
                 for r in fr if key in r and "rmsd" in r[key]
                 and r["pdb"] in cm and key in cm[r["pdb"]]]
        nm = [p[0] for p in pairs]
        b = np.array([p[1] for p in pairs]); a = np.array([p[2] for p in pairs])
        e = np.maximum(np.array([p[3] for p in pairs]), np.array([p[4] for p in pairs]))
        fo = np.array([p[5] for p in pairs])
        g = L.gated_paired(b, a, e, folds=fo, names=nm)
        d = b - a
        keep = L.gate(e)
        g["max_abs_ungated"] = float(np.abs(d).max())
        g["max_abs_gated"] = float(np.abs(d[keep]).max()) if keep.any() else None
        g["sd_ungated"] = float(d.std(ddof=1))
        g["sd_gated"] = float(d[keep].std(ddof=1)) if keep.sum() > 1 else None
        g["worst_targets"] = sorted(zip(np.abs(d).tolist(), nm), reverse=True)[:5]
        g["PASS"] = bool(g.get("gated") is not None
                         and abs(g["gated"]["mean_diff"]) <= L.FRAME_TOL_MEAN
                         and (g["max_abs_gated"] or 9e9) <= L.FRAME_TOL_MAX)
        nulls[f"frame{f}"] = g
    out["exact_frame_null"] = nulls

    # -- which targets the gate removes, and why
    bad = [(pdbs[i], float(E[i])) for i in np.where(~L.gate(E))[0]]
    out["non_converged"] = {"n": len(bad), "targets": bad,
                            "frac": float(len(bad) / max(len(E), 1))}
    return out


# ------------------------------------------------------------------ 3
def conditioning():
    """37.6% of cached AMBER single points exceed 1e6; the 99th percentile is 1e12-1e15."""
    allE, per = [], []
    for p in L.ENUM_FILES:
        z = np.load(p)
        pos, brep = L.binding_mask(z)
        e = np.asarray(z["amber_total"], float)[pos]         # BINDING-MASKED
        allE.append(e)
        per.append({"pdb": str(z["pdb"]), "n_binding": int(len(e)),
                    "frac_gt_1e6": float((e > 1e6).mean()),
                    "frac_gt_1e12": float((e > 1e12).mean()),
                    "max": float(e.max()), "p99": float(np.percentile(e, 99)),
                    "tail10_raw": L.tail_share(e),
                    "tail10_winsor99": L.tail_share(L.condition(e, "winsor99")),
                    "tail10_signedlog": L.tail_share(L.condition(e, "signed_log")),
                    "tail10_rank": L.tail_share(L.condition(e, "rank"))})
    E = np.concatenate(allE)
    pooled = {"n": int(len(E)), "frac_gt_1e6": float((E > 1e6).mean()),
              "frac_gt_1e12": float((E > 1e12).mean()), "max": float(E.max()),
              "p99": float(np.percentile(E, 99)), "p999": float(np.percentile(E, 99.9)),
              "median": float(np.median(E))}
    for how in ("raw", "winsor99", "signed_log", "rank"):
        x = E if how == "raw" else L.condition(E, how)
        pooled[f"tail10_{how}"] = L.tail_share(x)
        pooled[f"kurtosis_{how}"] = float(
            ((x - x.mean()) ** 4).mean() / max(((x - x.mean()) ** 2).mean() ** 2, 1e-300))
    # a monotone map preserves every ordering statement -- verified, not assumed
    sub = E[:20000]
    checks = {}
    for how in ("winsor99", "signed_log", "rank"):
        y = L.condition(sub, how)
        o1 = np.argsort(sub, kind="mergesort"); o2 = np.argsort(y, kind="mergesort")
        checks[how] = {"order_preserved": bool(np.array_equal(sub[o1], sub[o2])),
                       "spearman_with_raw": L.spearman(sub, y)}
    return {"pooled": pooled, "per_target": per, "monotonicity_check": checks}


# ------------------------------------------------------------------ 4
def binding_audit():
    per, tot = [], {"n_amber": 0, "n_kind0": 0, "n_binding": 0, "snap_in_kind0": 0}
    for p in L.ENUM_FILES:
        z = np.load(p)
        pos, rep = L.binding_mask(z)
        rmsd = np.asarray(z["rmsd"], float)
        idx = np.asarray(z["amber_idx"], int)
        k0 = idx[np.asarray(z["amber_kind"], int) == 0]
        snap = int(z["snap_index"])
        rr0 = rmsd[k0]
        rep.update(pdb=str(z["pdb"]),
                   snap_is_min_of_kind0=bool(len(rr0) and rmsd[snap] <= rr0.min()),
                   snap_pctile_in_kind0=float((rr0 < rmsd[snap]).mean()) if len(rr0) else None,
                   mean_rmsd_kind0=float(rr0.mean()) if len(rr0) else None,
                   mean_rmsd_binding=float(rmsd[idx[pos]].mean()))
        per.append(rep)
        for k in ("n_amber", "n_kind0", "n_binding"):
            tot[k] += rep[k]
        tot["snap_in_kind0"] += int(rep["snap_in_kind0"])
    tot["n_files"] = len(per)
    tot["obj_floor_read"] = False
    tot["obj_floor_path"] = FORBIDDEN + " -- RETRACTED, not opened"
    return {"totals": tot, "per_target": per}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", type=int, default=0,
                    help="run the STANDING frame null live on the first N targets")
    args = ap.parse_args()

    out = {"gate_and_null": gate_and_null(), "conditioning": conditioning(),
           "binding_rule": binding_audit()}
    if args.frame:
        rows, verdict = L.frame_null(I.targets()[:args.frame])
        out["standing_frame_null"] = {"rows": rows, "verdict": verdict}
    print("wrote", L.write("energy_gate", out))

    g = out["gate_and_null"]
    e = g.get("effect_k30", {})
    if e:
        print(f"\nk=30 accuracy effect (AMBER - projection):")
        print(f"  ungated n={e['ungated']['n']} {e['ungated']['mean_diff']:+.5f} "
              f"[{e['ungated']['ci95'][0]:+.5f},{e['ungated']['ci95'][1]:+.5f}] "
              f"W/L {e['ungated']['n_better']}/{e['ungated']['n_worse']}")
        print(f"  GATED   n={e['gated']['n']} {e['gated']['mean_diff']:+.5f} "
              f"[{e['gated']['ci95'][0]:+.5f},{e['gated']['ci95'][1]:+.5f}] "
              f"W/L {e['gated']['n_better']}/{e['gated']['n_worse']}  "
              f"excluded {e['n_excluded']}: {e['excluded']}")
    for f, v in g.get("exact_frame_null", {}).items():
        print(f"\nEXACT FRAME NULL {f} (zero by construction):")
        print(f"  ungated n={v['ungated']['n']} {v['ungated']['mean_diff']:+.5f} "
              f"[{v['ungated']['ci95'][0]:+.5f},{v['ungated']['ci95'][1]:+.5f}] "
              f"sd {v['sd_ungated']:.4f} max|d| {v['max_abs_ungated']:.4f}")
        if v.get("gated"):
            print(f"  GATED   n={v['gated']['n']} {v['gated']['mean_diff']:+.5f} "
                  f"[{v['gated']['ci95'][0]:+.5f},{v['gated']['ci95'][1]:+.5f}] "
                  f"sd {v['sd_gated']:.4f} max|d| {v['max_abs_gated']:.4f}  "
                  f"PASS={v['PASS']}")
    c = out["conditioning"]["pooled"]
    print(f"\nCONDITIONING, pooled n={c['n']:,} binding-masked AMBER single points")
    print(f"  frac>1e6 {c['frac_gt_1e6']:.3f}  frac>1e12 {c['frac_gt_1e12']:.4f}  "
          f"max {c['max']:.3g}  p99 {c['p99']:.3g}")
    for how in ("raw", "winsor99", "signed_log", "rank"):
        print(f"  top-10 variance share, {how:10s} {c['tail10_'+how]:.6f}   "
              f"kurtosis {c['kurtosis_'+how]:.4g}")
    b = out["binding_rule"]["totals"]
    print(f"\nBINDING RULE across {b['n_files']} caches: {b['n_amber']:,} AMBER points -> "
          f"{b['n_kind0']:,} kind==0 -> {b['n_binding']:,} binding; "
          f"snap inside kind==0 on {b['snap_in_kind0']}/{b['n_files']} files")
