"""SPRINT 16, ENERGY -- the causal ablation's read-out.

Reads `s16/results/energy_ablate_s*.json` and answers, on ONE fixed generated ensemble:

  * what does LEGACY contribute            (leg - ctrl, and leg - rnd, its matched control)
  * what does AMBER contribute             (amb - ctrl)
  * do they INTERACT                       ((legamb - amb) - (leg - ctrl))
  * does either improve ACCURACY           CA-RMSD, ORACLE label
  * does either improve STRUCTURAL VALIDITY the stereochemistry panel, separate axis
  * does any benefit survive a CLASSICAL CONTROL  (the matched-count random filter)

ACCURACY AND VALIDITY ARE NEVER SUBSTITUTED FOR ONE ANOTHER.  They are reported as two
tables, and every panel statistic is scored on the SAME structure whose RMSD is quoted.

Every AMBER-involving row is reported BOTH gated and ungated, with the excluded count
printed (`energy_lib.gated_paired` makes the pair mandatory).  Set DIVERSITY is printed
beside every aggregate, because a filter that concentrates the set improves the set-mean
proxy while degrading the output.

    python -m s16.energy_report
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s16 import energy_lib as L            # noqa: E402
from s16 import energy_ablate as A         # noqa: E402

#: panel axis -> +1 if a LARGER value is BETTER, -1 if a SMALLER value is better.
#: Without this the W/L columns read backwards on `rama_favoured` and `min_heavy`, which
#: is exactly the kind of sign slip that produced two Sprint 14 reporting defects.
PANEL_SIGN = {"rama_favoured": +1, "rama_allowed": +1, "rama_outlier": -1,
              "n_clash_2A": -1, "n_clash_2p6A": -1, "min_heavy": +1,
              "bond_strain": -1, "angle_strain": -1, "cis_frac": -1, "omega_dev": -1,
              "geom_rms_rel_dev": -1}
PANEL_KEYS = tuple(PANEL_SIGN)


def collect(rows):
    ok = [r for r in rows if "rmsd" in r.get("ctrlamber", {})
          and "rmsd" in r.get("legamber", {})
          and "rmsd" in r.get("rnd1amber", {}) and "rmsd" in r.get("rnd2amber", {})]
    names = [r["pdb"] for r in ok]
    fold = np.array([r["fold"] for r in ok])
    acc, ener, pan = {}, {}, {}
    for arm in ("ctrl", "leg", "rnd1", "rnd2",
                "ctrlamber", "legamber", "rnd1amber", "rnd2amber"):
        if arm.endswith("amber"):
            acc[arm] = np.array([r[arm]["rmsd"] for r in ok])
            ener[arm] = np.array([r[arm]["energy"] for r in ok])
            pan[arm] = {k: np.array([r[arm]["panel"].get(k, np.nan) for r in ok])
                        for k in PANEL_KEYS}
        else:
            acc[arm] = np.array([r[arm]["rmsd"] for r in ok])
            ener[arm] = np.full(len(ok), -1.0)          # no minimisation: always "converged"
            pan[arm] = {k: np.array([r[arm]["panel"].get(k, np.nan) for r in ok])
                        for k in PANEL_KEYS}
    # the mean of the two random-filter draws, as the matched classical control
    acc["rnd"] = (acc["rnd1"] + acc["rnd2"]) / 2.0
    acc["rndamber"] = (acc["rnd1amber"] + acc["rnd2amber"]) / 2.0
    ener["rnd"] = np.full(len(ok), -1.0)
    ener["rndamber"] = np.maximum(ener["rnd1amber"], ener["rnd2amber"])
    for k in PANEL_KEYS:
        pan.setdefault("rnd", {})[k] = (pan["rnd1"][k] + pan["rnd2"][k]) / 2.0
        pan.setdefault("rndamber", {})[k] = (pan["rnd1amber"][k] + pan["rnd2amber"][k]) / 2.0
    return ok, names, fold, acc, ener, pan


CONTRASTS = [
    ("AMBER main effect            amb - ctrl", "amb", "ctrl", "ctrlamber"),
    ("LEGACY main effect           leg - ctrl", "leg", "ctrl", None),
    ("LEGACY vs its RANDOM control leg - rnd", "leg", "rnd", None),
    ("RANDOM filter                rnd - ctrl", "rnd", "ctrl", None),
    ("LEGACY given AMBER on     legamb - amb", "legamb", "amb", "both"),
    ("AMBER given LEGACY on     legamb - leg", "legamb", "leg", "legamber"),
    ("BOTH vs neither          legamb - ctrl", "legamb", "ctrl", "legamber"),
    ("BOTH vs random+AMBER  legamb - rndamb", "legamb", "rndamber", "both"),
]
_ALIAS = {"amb": "ctrlamber", "legamb": "legamber", "rndamb": "rndamber"}


def _a(name):
    return _ALIAS.get(name, name)


def contrasts(acc, ener, names, fold, label=""):
    out = {}
    for tag, x, y, egate in CONTRASTS:
        ax, ay = acc[_a(x)], acc[_a(y)]
        if egate is None:
            e = np.full(len(ax), -1.0)
        elif egate == "both":
            e = np.maximum(ener[_a(x)], ener[_a(y)])
        else:
            e = ener[_a(egate)]
        out[tag] = L.gated_paired(ax, ay, e, folds=fold, names=names)
        out[tag]["mean_x"] = float(ax.mean()); out[tag]["mean_y"] = float(ay.mean())
    # -- the interaction, on the gated set
    e = np.maximum.reduce([ener["ctrlamber"], ener["legamber"]])
    keep = L.gate(e)
    inter = (acc["legamber"] - acc["ctrlamber"]) - (acc["leg"] - acc["ctrl"])
    z = np.zeros(len(inter))
    out["INTERACTION (legamb-amb)-(leg-ctrl)"] = {
        "ungated": I.paired(inter, z, folds=fold, names=names),
        "gated": I.paired(inter[keep], z[keep],
                          folds=fold[keep],
                          names=[names[i] for i in np.where(keep)[0]]),
        "n_excluded": int((~keep).sum()),
        "excluded": [names[i] for i in np.where(~keep)[0]]}
    return out


def run():
    rows = A.load_all()
    ok, names, fold, acc, ener, pan = collect(rows)
    print(f"n complete targets = {len(ok)} of {len(rows)} rows", flush=True)

    out = {"n_targets": len(ok), "pdbs": names,
           "arm_means": {k: float(v.mean()) for k, v in acc.items()},
           "arm_medians": {k: float(np.median(v)) for k, v in acc.items()},
           "accuracy": contrasts(acc, ener, names, fold)}

    # ---------------- VALIDITY, the separate axis
    val = {}
    for key in PANEL_KEYS:
        v = {}
        for tag, x, y, egate in CONTRASTS:
            ax, ay = pan[_a(x)][key], pan[_a(y)][key]
            m = np.isfinite(ax) & np.isfinite(ay)
            if m.sum() < 5:
                continue
            if egate is None:
                e = np.full(len(ax), -1.0)
            elif egate == "both":
                e = np.maximum(ener[_a(x)], ener[_a(y)])
            else:
                e = ener[_a(egate)]
            v[tag] = L.gated_paired(ax[m], ay[m], e[m], folds=fold[m],
                                    names=[names[i] for i in np.where(m)[0]])
            v[tag]["mean_x"] = float(ax[m].mean()); v[tag]["mean_y"] = float(ay[m].mean())
            # orientation-corrected win/loss: "improved" means BETTER on this axis
            g = v[tag].get("gated") or v[tag]["ungated"]
            s = PANEL_SIGN[key]
            g["n_improved"] = g["n_worse"] if s > 0 else g["n_better"]
            g["n_degraded"] = g["n_better"] if s > 0 else g["n_worse"]
            v[tag]["higher_is_better"] = bool(s > 0)
        val[key] = {"arm_means": {k: float(np.nanmean(pan[k][key])) for k in pan},
                    "higher_is_better": bool(PANEL_SIGN[key] > 0), "contrasts": v}
    out["validity"] = val

    # ---------------- what the filter did to the ENSEMBLE (set-mean law bookkeeping)
    ens = {}
    for tag in ("ctrl", "leg", "rnd1", "rnd2"):
        ens[tag] = {
            "set_diversity_meanpair_rmsd": float(np.mean([r[f"{tag}_setdiv"] for r in ok])),
            "ORACLE_set_mean_rmsd": float(np.mean([r[f"ORACLE_setmean_{tag}"] for r in ok])),
            "ORACLE_set_best_rmsd": float(np.mean([r[f"ORACLE_setbest_{tag}"] for r in ok])),
            "raw_average_rmsd": float(np.mean([r[f"{tag}_rawavg_rmsd"] for r in ok]))}
    out["ensemble"] = ens
    out["legacy_steric_gate"] = {
        "mean_frac_windows_flagged": float(np.mean([r["legacy_steric_nonzero"] for r in ok])),
        "max_frac": float(np.max([r["legacy_steric_nonzero"] for r in ok])),
        "n_targets_with_zero_flags": int(sum(r["legacy_steric_nonzero"] == 0 for r in ok)),
        "note": "why the pre-registered LEGACY-ON arm is a rank filter and not a pure "
                "clash gate: the clash gate is too sparse to move a 75-window average"}
    out["convergence"] = {
        arm: {"n_nonconverged": int((~L.gate(ener[arm])).sum()),
              "targets": [names[i] for i in np.where(~L.gate(ener[arm]))[0]]}
        for arm in ("ctrlamber", "legamber", "rnd1amber", "rnd2amber")}
    out["cost_s"] = {
        "projection_mean": float(np.mean([r["ctrl"]["wall"] for r in ok])),
        "amber_minimisation_mean": float(np.mean([r["ctrlamber"]["wall"] for r in ok])),
        "amber_minimisation_median": float(np.median([r["ctrlamber"]["wall"] for r in ok]))}
    # ---------------- the ledger's mandated concentration check, on every accuracy
    # contrast whose CI excludes zero with a near-even W/L.  A raw drop-top threshold is
    # NOT a test; `concentration_null` calibrates it against a UNIFORM-effect null.
    from s14.ener_avgrefine import concentration_null
    conc = {}
    for tag, x, y, _eg in CONTRASTS:
        d = out["accuracy"][tag]
        g = d.get("gated") or d.get("ungated")
        lo, hi = g["ci95"]
        wl = g["n_better"] / max(g["n_better"] + g["n_worse"], 1)
        if lo * hi <= 0 or not (0.35 <= wl <= 0.65):
            continue     # only fires on the ledger's early warning: near-even W/L + CI
        dd = acc[_a(x)] - acc[_a(y)]
        nl = concentration_null(dd)
        o = np.argsort(dd)
        obs10 = float(dd[o[10:]].mean())
        conc[tag] = {"observed_mean": float(dd.mean()),
                     "observed_median": float(np.median(dd)),
                     "observed_drop_top10": obs10,
                     "null_drop_top10_mean": nl["dt10"], "null_ci": nl["dt10_ci"],
                     "observed_percentile_in_null":
                         float((np.asarray(nl["dt10_draws"]) < obs10).mean()),
                     "verdict": ("PASS (not concentrated)"
                                 if nl["dt10_ci"][0] <= obs10 <= nl["dt10_ci"][1]
                                 else "CONCENTRATED")}
    out["concentration_null"] = conc
    return out


def _line(tag, d):
    g = d.get("gated") or d["ungated"]
    u = d["ungated"]
    return (f"  {tag:40s} {g['mean_diff']:+8.4f} [{g['ci95'][0]:+.4f},{g['ci95'][1]:+.4f}] "
            f"med {g['median_diff']:+8.4f} W/L {g['n_better']:3d}/{g['n_worse']:<3d} "
            f"n={g['n']:3d}  (ungated {u['mean_diff']:+.4f} n={u['n']}, "
            f"excl {d.get('n_excluded', 0)})")


if __name__ == "__main__":
    out = run()
    print("\nwrote", L.write("energy_report", out))
    print(f"\nARM MEANS (CA-RMSD, ORACLE):")
    for k, v in out["arm_means"].items():
        print(f"  {k:12s} {v:.4f}   median {out['arm_medians'][k]:.4f}")
    print(f"\nACCURACY (CA-RMSD, ORACLE label; negative = better)")
    for tag, d in out["accuracy"].items():
        print(_line(tag, d))
    print(f"\nVALIDITY (separate axis; each on the SAME structure whose RMSD is quoted)")
    for key, v in out["validity"].items():
        am = v["arm_means"]
        print(f"\n {key} ({'higher' if v['higher_is_better'] else 'lower'} is better):  "
              f"ctrl {am['ctrl']:.4f}  amb {am['ctrlamber']:.4f}  "
              f"leg {am['leg']:.4f}  legamb {am['legamber']:.4f}  rnd {am['rnd']:.4f}")
        for tag, d in v["contrasts"].items():
            g = d.get("gated") or d["ungated"]
            print(_line(tag, d) + f" [improved {g['n_improved']}/{g['n_degraded']}]")
    print("\nENSEMBLE (what the filter did, with DIVERSITY beside every aggregate):")
    for k, v in out["ensemble"].items():
        print(f"  {k:6s} div {v['set_diversity_meanpair_rmsd']:.3f}  "
              f"ORACLE setmean {v['ORACLE_set_mean_rmsd']:.3f}  "
              f"setbest {v['ORACLE_set_best_rmsd']:.3f}  "
              f"rawavg {v['raw_average_rmsd']:.3f}")
    print("\nconvergence:", {k: v["n_nonconverged"] for k, v in out["convergence"].items()})
    print("cost (s):", {k: round(v, 2) for k, v in out["cost_s"].items()})
    print("legacy steric gate:", {k: v for k, v in out["legacy_steric_gate"].items()
                                  if k != "note"})
