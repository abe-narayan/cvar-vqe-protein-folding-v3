"""SPRINT 16 / QPHASE -- every table, regenerated from the checkpoints.

Nothing in `s16/qphase_FINDINGS.md` is typed by hand: this module reads
`s16/results/qphase_*.json` and prints the tables, so a partial run is readable and a
finished one is reproducible.  The statistical unit is always the TARGET -- seeds are
averaged within a target before any interval is formed.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from s16 import qphase_lib as L
from s12 import instrument as I

ARMS = ("untrained", "tilt_samples", "tilt_exact", "anneal", "anneal_q", "greedy", "random")
READOUTS = ("rmsd_returned", "rand5_coordavg_rmsd", "rand20_coordavg_rmsd",
            "set_mean_rmsd", "set_best_rmsd", "objective_gap", "n_distinct")
FLOOR = 0.08


# ------------------------------------------------------------------ the long table
def long_table(family="blend"):
    """rows: one per (rung, target), seeds averaged WITHIN the target."""
    d = L.ck_load(f"ladder_{family}")
    rows = []
    for pdb, byrung in d.items():
        if pdb.startswith("_"):
            continue
        ins_fold = None
        for rung, blob in byrung.items():
            p = blob["profile"]
            r = {"pdb": pdb, "rung": float(rung), "family": family,
                 "rho_global": p["rho_global"], "rho_tail": p["rho_tail_0.01"],
                 "rho_tail_null_sd": p["rho_tail_0.01_null_sd"],
                 "argmin_pct": p["argmin_pct"], "argmin_rmsd": p["argmin_rmsd"],
                 "tail_gap": p["tail_0.01_gap"],
                 "rmsd_best_in_space": p["rmsd_best_in_space"],
                 "rmsd_mean_in_space": p["rmsd_mean_in_space"]}
            seeds = blob["seeds"]
            for arm in ("vqe",) + ARMS:
                vals = [s[arm] for s in seeds.values() if arm in s]
                if not vals:
                    continue
                for k in READOUTS + ("draw_entropy_bits", "cost_evals", "argmin_ties"):
                    if k in vals[0]:
                        r[f"{arm}.{k}"] = float(np.mean([v[k] for v in vals]))
            rows.append(r)
    return rows


def _folds(pdbs):
    return np.asarray([L.inst(p).fold for p in pdbs])


def rung_stats(rows, arm, readout, family):
    """Per-rung paired VQE-minus-control, unit = TARGET. Positive = VQE WORSE."""
    out = []
    rungs = sorted({r["rung"] for r in rows if r["family"] == family})
    for g in rungs:
        sub = [r for r in rows if r["family"] == family and r["rung"] == g
               and f"vqe.{readout}" in r and f"{arm}.{readout}" in r]
        if len(sub) < 6:
            continue
        pdbs = [r["pdb"] for r in sub]
        a = np.asarray([r[f"vqe.{readout}"] for r in sub])
        b = np.asarray([r[f"{arm}.{readout}"] for r in sub])
        st = I.paired(a, b, folds=_folds(pdbs), names=pdbs, seed=L.seed_of(f"{arm}{g}"))
        rec = {"rung": g, "n": st["n"],
               "rho_global": float(np.mean([r["rho_global"] for r in sub])),
               "rho_tail": float(np.mean([r["rho_tail"] for r in sub])),
               "argmin_pct": float(np.mean([r["argmin_pct"] for r in sub])),
               "vqe": st["mean_a"], "ctrl": st["mean_b"], "diff": st["mean_diff"],
               "median": st["median_diff"], "ci": st["ci95"],
               "W": st["n_better"], "L": st["n_worse"],
               "sig": bool(st["ci95"][1] < 0 or st["ci95"][0] > 0),
               "below_floor": bool(abs(st["mean_diff"]) <= FLOOR)}
        if rec["sig"] and st["n"] > 10:
            try:
                c = L.nullconc(a - b, seed=L.seed_of(f"nc{arm}{g}{readout}"))
                rec["conc"] = c["verdict"]
                rec["conc_p_share"] = c["p_share_vs_null"]
                rec["mean_over_sd"] = c["mean_over_sd"]
            except Exception as ex:
                rec["conc"] = f"error {ex}"
        out.append(rec)
    return out


def rho_star(stats, xkey="rho_global"):
    """The LOWEST-quality rung at which the VQE beats the control with an interval excluding
    zero AND an effect above the empirical false-positive floor.

    The scan is ordered by the QUALITY AXIS, not by the knob: the noise family's rungs are
    listed by `sigma`, which runs the opposite way to quality, and scanning in knob order
    would report the wrong end of the ladder.  `None` means no rung on the ladder crosses --
    which is itself the answer, and the answer for annealing and greedy.
    """
    for r in sorted(stats, key=lambda z: z[xkey]):
        if r["diff"] < 0 and r["sig"] and not r["below_floor"]:
            return {"rho_star": r[xkey], "rung": r["rung"], "diff": r["diff"],
                    "ci": r["ci"], "W": r["W"], "L": r["L"],
                    "conc": r.get("conc", "n/a")}
    return None


def crossing(rows, arm, readout, family, xkey="rho_global", nboot=2000, seed=0):
    """CONTINUOUS estimate of the boundary: cell-level OLS of the paired difference on the
    quality axis, with a CLUSTER bootstrap over TARGETS (never over cells -- cells inside a
    target share a target and are not independent).  Returns the x at which the fit crosses
    zero, and the CI of that crossing."""
    sub = [r for r in rows if r["family"] == family and f"vqe.{readout}" in r
           and f"{arm}.{readout}" in r]
    if len(sub) < 20:
        return None
    pdbs = sorted({r["pdb"] for r in sub})
    by = {p: [r for r in sub if r["pdb"] == p] for p in pdbs}

    def fit(rs):
        x = np.asarray([r[xkey] for r in rs])
        y = np.asarray([r[f"vqe.{readout}"] - r[f"{arm}.{readout}"] for r in rs])
        A = np.vstack([np.ones_like(x), x]).T
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        return c

    c0 = fit(sub)
    if abs(c0[1]) < 1e-12:
        return {"slope": float(c0[1]), "intercept": float(c0[0]), "cross": None}
    x0 = -c0[0] / c0[1]
    rng = np.random.default_rng(seed)
    xs = []
    for _ in range(nboot):
        pick = rng.choice(len(pdbs), len(pdbs))
        rs = [r for k in pick for r in by[pdbs[k]]]
        c = fit(rs)
        if abs(c[1]) > 1e-12:
            xs.append(-c[0] / c[1])
    xs = np.asarray(xs)
    xs = xs[np.isfinite(xs)]
    return {"slope": float(c0[1]), "intercept": float(c0[0]), "cross": float(x0),
            "cross_ci": [float(np.percentile(xs, 2.5)), float(np.percentile(xs, 97.5))]
            if xs.size > 50 else None,
            "frac_cross_in_unit": float(((xs >= -1) & (xs <= 1)).mean()) if xs.size else None,
            "n_cells": len(sub), "n_targets": len(pdbs)}


# ------------------------------------------------------------------------- printing
def fmt_stats(name, stats, xkey="rho_global"):
    lines = [f"  {'rung':>6} {'rho_g':>7} {'rho_t':>7} {'apct':>6} "
             f"{'VQE':>7} {'ctrl':>7} {'diff':>8} {'median':>8} "
             f"{'CI95':>20} {'W/L':>7} {'verdict':>10}"]
    for r in stats:
        v = ("SIG" if r["sig"] else "null")
        if r["below_floor"]:
            v += "<=FLOOR"
        lines.append(f"  {r['rung']:>6.3g} {r['rho_global']:>+7.3f} {r['rho_tail']:>+7.3f} "
                     f"{r['argmin_pct']:>6.3f} {r['vqe']:>7.3f} {r['ctrl']:>7.3f} "
                     f"{r['diff']:>+8.3f} {r['median']:>+8.3f} "
                     f"[{r['ci'][0]:+.3f},{r['ci'][1]:+.3f}]".rjust(20) +
                     f" {r['W']:>3}/{r['L']:<3} {v:>10}"
                     + (f" {r.get('conc','')}" if r["sig"] else ""))
    return f"\n{name}\n" + "\n".join(lines)


def boundary(readout="rmsd_returned", families=("blend", "noise")):
    out = {}
    for fam in families:
        rows = long_table(fam)
        if not rows:
            continue
        out[fam] = {}
        for arm in ARMS:
            st = rung_stats(rows, arm, readout, fam)
            if not st:
                continue
            out[fam][arm] = {
                "stats": st,
                "rho_star_global": rho_star(st, "rho_global"),
                "rho_star_tail": rho_star(st, "rho_tail"),
                "rho_star_argminpct": rho_star(st, "argmin_pct"),
                "crossing_rho_global": crossing(rows, arm, readout, fam, "rho_global"),
                "crossing_rho_tail": crossing(rows, arm, readout, fam, "rho_tail"),
                "crossing_argmin_pct": crossing(rows, arm, readout, fam, "argmin_pct"),
            }
    return out


# ------------------------------------------------- WHICH QUALITY AXIS GOVERNS THE BOUNDARY
def axis_test(arm="untrained", readout="rmsd_returned", nboot=2000, seed=0):
    """The blend and the noise families DISSOCIATE the candidate quality axes, and that is
    what makes this a test rather than a description.

    A blend at rho = 0.5 has an argmin that is still the TRUE objective's near-optimum (the
    real objective's structure is underneath and the truth's ordering dominates the extreme
    tail), while a noised truth at the SAME rho = 0.5 has an argmin chosen by the largest
    negative noise draw over 262,144 configurations, i.e. essentially at random.  So the two
    families reach identical rho with very different argmin locations.  Whichever axis the
    boundary is a function of, the two families must AGREE on it.

    Reported: cell-level R^2 of the paired difference on each candidate axis, pooled over
    both families and separately within each, with a TARGET-cluster bootstrap.
    """
    rows = long_table("blend") + long_table("noise")
    rows = [r for r in rows if f"vqe.{readout}" in r and f"{arm}.{readout}" in r]
    if not rows:
        return {}
    y = np.asarray([r[f"vqe.{readout}"] - r[f"{arm}.{readout}"] for r in rows])
    out = {"n_cells": len(rows), "arm": arm, "readout": readout}
    axes = {"rho_global": [r["rho_global"] for r in rows],
            "rho_tail": [r["rho_tail"] for r in rows],
            "argmin_pct": [r["argmin_pct"] for r in rows],
            "log10_argmin_pct": [np.log10(max(r["argmin_pct"], 1e-6)) for r in rows],
            "argmin_rmsd": [r["argmin_rmsd"] for r in rows]}
    pdbs = sorted({r["pdb"] for r in rows})
    idx_by = {p: np.flatnonzero(np.asarray([r["pdb"] for r in rows]) == p) for p in pdbs}
    rng = np.random.default_rng(seed)
    for nm, x in axes.items():
        x = np.asarray(x, float)
        A = np.vstack([np.ones_like(x), x]).T
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        r2 = 1.0 - ((y - A @ c) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-30)
        bs = []
        for _ in range(nboot):
            pick = np.concatenate([idx_by[pdbs[k]]
                                   for k in rng.choice(len(pdbs), len(pdbs))])
            xx, yy = x[pick], y[pick]
            AA = np.vstack([np.ones_like(xx), xx]).T
            cc, *_ = np.linalg.lstsq(AA, yy, rcond=None)
            bs.append(cc[1])
        bs = np.asarray(bs)
        out[nm] = {"slope": float(c[1]), "intercept": float(c[0]), "r2": float(r2),
                   "slope_ci": [float(np.percentile(bs, 2.5)),
                                float(np.percentile(bs, 97.5))],
                   "rho_spearman_with_diff": float(_sp(x, y))}
        for fam in ("blend", "noise"):
            m = np.asarray([r["family"] == fam for r in rows])
            if m.sum() > 20:
                xx, yy = x[m], y[m]
                AA = np.vstack([np.ones_like(xx), xx]).T
                cc, *_ = np.linalg.lstsq(AA, yy, rcond=None)
                r2f = 1.0 - ((yy - AA @ cc) ** 2).sum() / \
                    max(((yy - yy.mean()) ** 2).sum(), 1e-30)
                out[nm][f"{fam}_slope"] = float(cc[1])
                out[nm][f"{fam}_r2"] = float(r2f)
                out[nm][f"{fam}_cross"] = (float(-cc[0] / cc[1]) if abs(cc[1]) > 1e-12
                                           else None)
    return out


def _sp(a, b):
    from s14.vqe_lib import spearman
    return spearman(a, b)


# ------------------------------------------------------------------- T2: REAL OBJECTIVES
def realarms_table(readout="rmsd_returned"):
    d = L.ck_load("realarms")
    objs = {}
    for pdb, byobj in d.items():
        if pdb.startswith("_"):
            continue
        for nm, blob in byobj.items():
            r = {"pdb": pdb, "rho_global": blob["profile"]["rho_global"],
                 "rho_tail": blob["profile"]["rho_tail_0.01"],
                 "rho_tail_null_sd": blob["profile"]["rho_tail_0.01_null_sd"],
                 "argmin_pct": blob["profile"]["argmin_pct"],
                 "argmin_rmsd": blob["profile"]["argmin_rmsd"]}
            for arm in ("vqe",) + ARMS:
                vals = [s[arm] for s in blob["seeds"].values() if arm in s]
                if vals:
                    for k in READOUTS:
                        if k in vals[0]:
                            r[f"{arm}.{k}"] = float(np.mean([v[k] for v in vals]))
            objs.setdefault(nm, []).append(r)
    out = {}
    for nm, rs in objs.items():
        pdbs = [r["pdb"] for r in rs]
        rec = {"n": len(rs),
               "rho_global": float(np.mean([r["rho_global"] for r in rs])),
               "rho_global_sd": float(np.std([r["rho_global"] for r in rs])),
               "rho_tail": float(np.mean([r["rho_tail"] for r in rs])),
               "rho_tail_z": float(np.mean([r["rho_tail"] / max(r["rho_tail_null_sd"], 1e-9)
                                            for r in rs])),
               "argmin_pct": float(np.mean([r["argmin_pct"] for r in rs])),
               "argmin_rmsd": float(np.mean([r["argmin_rmsd"] for r in rs]))}
        for arm in ARMS:
            a = np.asarray([r[f"vqe.{readout}"] for r in rs if f"{arm}.{readout}" in r])
            b = np.asarray([r[f"{arm}.{readout}"] for r in rs if f"{arm}.{readout}" in r])
            if a.size < 6:
                continue
            st = I.paired(a, b, folds=_folds(pdbs), names=pdbs,
                          seed=L.seed_of(f"ra{nm}{arm}"))
            rec[f"vs_{arm}"] = {"diff": st["mean_diff"], "median": st["median_diff"],
                                "ci": st["ci95"], "W": st["n_better"], "L": st["n_worse"],
                                "sig": bool(st["ci95"][1] < 0 or st["ci95"][0] > 0),
                                "below_floor": bool(abs(st["mean_diff"]) <= FLOOR)}
        out[nm] = rec
    return out


# ---------------------------------------------------------------------- T3: THE LR LADDER
def lr_table():
    d = L.ck_load("lr")
    arms = None
    out = {}
    for pdb, byrung in d.items():
        if pdb.startswith("_"):
            continue
        for rung, blob in byrung.items():
            for s, cell in blob["seeds"].items():
                for arm, r in cell.items():
                    out.setdefault((rung, arm), {}).setdefault(pdb, []).append(r)
    agg = {}
    for (rung, arm), bypdb in out.items():
        rec = {"rung": rung, "arm": arm, "n": len(bypdb)}
        for k in ("rmsd_returned", "rand5_coordavg_rmsd", "n_distinct",
                  "draw_entropy_bits", "objective_gap", "grad_norm_first",
                  "grad_norm_mean", "param_displacement", "set_mean_rmsd"):
            v = [float(np.mean([x[k] for x in rs])) for rs in bypdb.values()
                 if k in rs[0]]
            if v:
                rec[k] = float(np.mean(v))
        agg[f"{rung}|{arm}"] = rec
    # paired VQE-minus-untrained per arm, unit = target
    for (rung, arm), bypdb in out.items():
        if arm == "untrained":
            continue
        ctrl = out.get((rung, "untrained"))
        if ctrl is None:
            continue
        pdbs = sorted(set(bypdb) & set(ctrl))
        for k in ("rmsd_returned", "rand5_coordavg_rmsd"):
            a = np.asarray([np.mean([x[k] for x in bypdb[p]]) for p in pdbs])
            b = np.asarray([np.mean([x[k] for x in ctrl[p]]) for p in pdbs])
            if a.size < 6:
                continue
            st = I.paired(a, b, folds=_folds(pdbs), names=pdbs,
                          seed=L.seed_of(f"lr{rung}{arm}{k}"))
            agg[f"{rung}|{arm}"][f"vs_untrained_{k}"] = {
                "diff": st["mean_diff"], "median": st["median_diff"], "ci": st["ci95"],
                "W": st["n_better"], "L": st["n_worse"],
                "sig": bool(st["ci95"][1] < 0 or st["ci95"][0] > 0),
                "below_floor": bool(abs(st["mean_diff"]) <= FLOOR)}
    return agg


# ------------------------------------------------------------------------ T4: WHAT ALPHA IS
def alpha_table():
    d = L.ck_load("alpha")
    out = {}
    for pdb, byrung in d.items():
        if pdb.startswith("_"):
            continue
        for rung, blob in byrung.items():
            for s, cell in blob["seeds"].items():
                for arm, r in cell.items():
                    out.setdefault((rung, arm), {}).setdefault(pdb, []).append(r)
    agg = {}
    for (rung, arm), bypdb in out.items():
        rec = {"rung": rung, "arm": arm, "n": len(bypdb)}
        for k in ("rmsd_returned", "rand5_coordavg_rmsd", "rand20_coordavg_rmsd",
                  "set_mean_rmsd", "set_best_rmsd", "n_distinct", "draw_entropy_bits",
                  "objective_gap"):
            v = [float(np.mean([x[k] for x in rs])) for rs in bypdb.values() if k in rs[0]]
            if v:
                rec[k] = float(np.mean(v))
        agg[f"{rung}|{arm}"] = rec
    # the decisive pairing: VQE at alpha a versus the entropy-matched classical tilt
    for (rung, arm), bypdb in out.items():
        if not arm.startswith("vqe_a"):
            continue
        a_lab = arm[len("vqe_a"):]
        for ctrl_name in (f"tiltmatch_a{a_lab}", f"tiltx_a{a_lab}", "untrained"):
            ctrl = out.get((rung, ctrl_name))
            if ctrl is None:
                continue
            pdbs = sorted(set(bypdb) & set(ctrl))
            for k in ("rmsd_returned", "rand5_coordavg_rmsd", "set_mean_rmsd"):
                a = np.asarray([np.mean([x[k] for x in bypdb[p]]) for p in pdbs])
                b = np.asarray([np.mean([x[k] for x in ctrl[p]]) for p in pdbs])
                if a.size < 6:
                    continue
                st = I.paired(a, b, folds=_folds(pdbs), names=pdbs,
                              seed=L.seed_of(f"al{rung}{arm}{ctrl_name}{k}"))
                agg[f"{rung}|{arm}"][f"vs_{ctrl_name}_{k}"] = {
                    "diff": st["mean_diff"], "median": st["median_diff"],
                    "ci": st["ci95"], "W": st["n_better"], "L": st["n_worse"],
                    "sig": bool(st["ci95"][1] < 0 or st["ci95"][0] > 0),
                    "below_floor": bool(abs(st["mean_diff"]) <= FLOOR)}
    return agg


# ------------------------------------------------------------- T5: WHAT THE VQE IS DOING
def doing_table():
    d = L.ck_load("doing")
    out = {}
    for pdb, byrung in d.items():
        if pdb.startswith("_"):
            continue
        for rung, blob in byrung.items():
            for s, cell in blob["seeds"].items():
                out.setdefault(rung, {}).setdefault(pdb, []).append(
                    dict(cell, rho_global=blob["profile"]["rho_global"],
                         argmin_pct=blob["profile"]["argmin_pct"]))
    agg = {}
    for rung, bypdb in out.items():
        rec = {"rung": rung, "n": len(bypdb)}
        keys = set()
        for rs in bypdb.values():
            keys |= set(rs[0])
        for k in sorted(keys):
            v = [float(np.mean([x[k] for x in rs])) for rs in bypdb.values() if k in rs[0]]
            if v:
                rec[k] = float(np.mean(v))
        agg[rung] = rec
    return agg


# ------------------------------------------------------------------- THE SCHEDULE CONTROL
SCHED_ARMS = ("untrained", "tilt_samples", "tilt_exact", "tilt_anneal", "anneal")


def sched_table():
    d = L.ck_load("sched")
    rows = {}
    for pdb, byr in d.items():
        if pdb.startswith("_"):
            continue
        for g, blob in byr.items():
            rec = {"pdb": pdb, "rho": blob["profile"]["rho_global"],
                   "apct": blob["profile"]["argmin_pct"]}
            for a in ("vqe",) + SCHED_ARMS:
                v = [s[a] for s in blob["seeds"].values() if a in s]
                for k in ("rmsd_returned", "rand5_coordavg_rmsd", "n_distinct"):
                    rec[f"{a}|{k}"] = float(np.mean([x[k] for x in v]))
            rows.setdefault(g, []).append(rec)
    out = {}
    for g, s in rows.items():
        pdbs = [r["pdb"] for r in s]
        rec = {"rung": g, "n": len(s), "rho": float(np.mean([r["rho"] for r in s]))}
        for k in ("rmsd_returned", "rand5_coordavg_rmsd", "n_distinct"):
            rec[f"vqe|{k}"] = float(np.mean([r[f"vqe|{k}"] for r in s]))
            for a in SCHED_ARMS:
                rec[f"{a}|{k}"] = float(np.mean([r[f"{a}|{k}"] for r in s]))
                if k == "n_distinct":
                    continue
                A = np.asarray([r[f"vqe|{k}"] for r in s])
                B = np.asarray([r[f"{a}|{k}"] for r in s])
                st = I.paired(A, B, folds=_folds(pdbs), names=pdbs,
                              seed=L.seed_of(f"sched{g}{a}{k}"))
                rec[f"vs_{a}_{k}"] = {"diff": st["mean_diff"], "ci": st["ci95"],
                                      "W": st["n_better"], "L": st["n_worse"],
                                      "sig": bool(st["ci95"][1] < 0 or st["ci95"][0] > 0),
                                      "below_floor": bool(abs(st["mean_diff"]) <= FLOOR)}
        out[g] = rec
    return out


def alpha_curves():
    """Does a one-parameter classical thermostat reproduce alpha's whole effect?

    The correlation, across the six alphas, between the VQE's readout curve and each
    classical tilt's readout curve, plus each curve's own range.  This is §6b of the
    findings, regenerated rather than typed.
    """
    t = alpha_table()
    ALS = (0.02, 0.05, 0.1, 0.25, 0.5, 1.0)
    out = {}
    for rung in ("0.0", "1.0"):
        if f"{rung}|vqe_a0.25" not in t:
            continue
        rec = {}
        for k in ("rand5_coordavg_rmsd", "set_mean_rmsd", "rmsd_returned",
                  "n_distinct", "draw_entropy_bits"):
            v = np.asarray([t[f"{rung}|vqe_a{a}"][k] for a in ALS])
            rec[k] = {"vqe_range": float(v.max() - v.min()),
                      "vqe_curve": v.tolist()}
            for ctrl in ("tiltx", "tiltmatch"):
                c = np.asarray([t[f"{rung}|{ctrl}_a{a}"][k] for a in ALS])
                rec[k][f"{ctrl}_pearson"] = float(np.corrcoef(v, c)[0, 1])
                rec[k][f"{ctrl}_range"] = float(c.max() - c.min())
                rec[k][f"{ctrl}_curve"] = c.tolist()
        out[rung] = rec
    return out


def real_axis_test():
    """§5b: which candidate quality axis predicts the VQE-minus-untrained difference on the
    REAL objectives.  95 (objective x target) cells, no oracle knob anywhere in the arms."""
    from s14.vqe_lib import spearman
    d = L.ck_load("realarms")
    rows = []
    for pdb, byobj in d.items():
        if pdb.startswith("_"):
            continue
        for nm, blob in byobj.items():
            p = blob["profile"]
            r = {"pdb": pdb, "obj": nm, "rho_global": p["rho_global"],
                 "rho_tail": p["rho_tail_0.01"], "argmin_pct": p["argmin_pct"],
                 "argmin_rmsd": p["argmin_rmsd"]}
            for a in ("vqe", "untrained"):
                v = [s[a] for s in blob["seeds"].values() if a in s]
                r[a] = float(np.mean([x["rmsd_returned"] for x in v]))
            rows.append(r)
    y = np.asarray([r["vqe"] - r["untrained"] for r in rows])
    out = {"n_cells": len(rows)}
    for nm in ("rho_global", "rho_tail", "argmin_pct", "argmin_rmsd"):
        x = np.asarray([r[nm] for r in rows])
        out[nm] = {"spearman": float(spearman(x, y)),
                   "pearson": float(np.corrcoef(x, y)[0, 1])}
    x = np.asarray([r["argmin_rmsd"] - r["untrained"] for r in rows])
    out["argmin_rmsd_minus_control"] = {"spearman": float(spearman(x, y)),
                                        "pearson": float(np.corrcoef(x, y)[0, 1])}
    return out


def main():
    what = sys.argv[2] if len(sys.argv) > 2 else "all"
    if what in ("all", "boundary"):
        for readout in ("rmsd_returned", "rand5_coordavg_rmsd", "objective_gap"):
            b = boundary(readout)
            print("=" * 120)
            print(f"READOUT: {readout}   (positive diff = VQE WORSE; unit = TARGET)")
            for fam, byarm in b.items():
                for arm, blob in byarm.items():
                    print(fmt_stats(f"[{fam}] vqe - {arm}", blob["stats"]))
                    rs = blob["rho_star_global"]
                    print(f"    rho* (global) = "
                          f"{('%.3f' % rs['rho_star']) if rs else 'NO CROSSING ON LADDER'}"
                          f"   continuous crossing = {blob['crossing_rho_global']}")
            L.ck("boundary", readout, b)
    if what in ("all", "real"):
        print("=" * 120)
        print("REAL OBJECTIVES")
        d = L.ck_load("real")
        agg = {}
        for pdb, row in d.items():
            if pdb.startswith("_"):
                continue
            for nm, prof in row.items():
                if not isinstance(prof, dict) or "rho_global" not in prof:
                    continue
                agg.setdefault(nm, []).append(prof)
        for nm, ps in sorted(agg.items()):
            print(f"  {nm:>14} n={len(ps):>2} "
                  f"rho_global {np.mean([p['rho_global'] for p in ps]):+.3f} "
                  f"rho_tail1% {np.mean([p['rho_tail_0.01'] for p in ps]):+.3f} "
                  f"argmin_pct {np.mean([p['argmin_pct'] for p in ps]):.3f} "
                  f"argmin_rmsd {np.mean([p['argmin_rmsd'] for p in ps]):.3f}")
    if what in ("all", "axis"):
        print("=" * 120)
        print("WHICH QUALITY AXIS GOVERNS THE BOUNDARY (blend and noise pooled)")
        for arm in ("untrained", "anneal"):
            for ro in ("rmsd_returned", "rand5_coordavg_rmsd"):
                a = axis_test(arm, ro)
                if a:
                    print(f"\n  arm={arm} readout={ro} cells={a['n_cells']}")
                    for k in ("rho_global", "rho_tail", "argmin_pct",
                              "log10_argmin_pct", "argmin_rmsd"):
                        v = a[k]
                        print(f"    {k:>17} R2 {v['r2']:.3f} slope {v['slope']:+.3f} "
                              f"[{v['slope_ci'][0]:+.3f},{v['slope_ci'][1]:+.3f}] "
                              f"blendR2 {v.get('blend_r2', float('nan')):.3f} "
                              f"noiseR2 {v.get('noise_r2', float('nan')):.3f} "
                              f"blendX {v.get('blend_cross')} noiseX {v.get('noise_cross')}")
                    L.ck("axis", f"{arm}|{ro}", a)
    if what in ("all", "realarms"):
        print("=" * 120)
        print("T2 DIRECT: the REAL objectives through the identical arm table")
        for ro in ("rmsd_returned", "rand5_coordavg_rmsd"):
            t = realarms_table(ro)
            print(f"\n  readout {ro}")
            for nm, rec in sorted(t.items(), key=lambda kv: -kv[1]["rho_global"]):
                s = (f"  {nm:>16} n={rec['n']:>2} rho {rec['rho_global']:+.3f} "
                     f"(sd {rec['rho_global_sd']:.3f}) tail {rec['rho_tail']:+.3f} "
                     f"z {rec['rho_tail_z']:+.1f} apct {rec['argmin_pct']:.3f} "
                     f"argminRMSD {rec['argmin_rmsd']:.2f} |")
                for arm in ARMS:
                    v = rec.get(f"vs_{arm}")
                    if v:
                        s += (f" {arm} {v['diff']:+.3f}"
                              f"{'*' if v['sig'] and not v['below_floor'] else ''}")
                print(s)
            L.ck("realarms_table", ro, t)
    if what in ("all", "lr"):
        print("=" * 120)
        print("T3 THE WEAK-OPTIMISER TEST")
        t = lr_table()
        for k in sorted(t):
            r = t[k]
            v = r.get("vs_untrained_rmsd_returned", {})
            v2 = r.get("vs_untrained_rand5_coordavg_rmsd", {})
            print(f"  {k:>22} |g|1 {r.get('grad_norm_first', float('nan')):>7.2f} "
                  f"|g|m {r.get('grad_norm_mean', float('nan')):>7.2f} "
                  f"dtheta {r.get('param_displacement', float('nan')):>6.2f} "
                  f"ndist {r.get('n_distinct', float('nan')):>7.0f} "
                  f"H {r.get('draw_entropy_bits', float('nan')):>5.2f} "
                  f"argmin {r.get('rmsd_returned', float('nan')):.3f} "
                  f"({v.get('diff', float('nan')):+.3f}"
                  f"{'*' if v.get('sig') and not v.get('below_floor') else ' '}) "
                  f"rand5 {r.get('rand5_coordavg_rmsd', float('nan')):.3f} "
                  f"({v2.get('diff', float('nan')):+.3f}"
                  f"{'*' if v2.get('sig') and not v2.get('below_floor') else ' '})")
        L.ck("lr_table", "all", t)
    if what in ("all", "alpha"):
        print("=" * 120)
        print("T4 WHAT ALPHA IS")
        t = alpha_table()
        for k in sorted(t):
            r = t[k]
            s = (f"  {k:>24} ndist {r.get('n_distinct', float('nan')):>7.0f} "
                 f"H {r.get('draw_entropy_bits', float('nan')):>5.2f} "
                 f"argmin {r.get('rmsd_returned', float('nan')):.3f} "
                 f"rand5 {r.get('rand5_coordavg_rmsd', float('nan')):.3f} "
                 f"setmean {r.get('set_mean_rmsd', float('nan')):.3f}")
            for c in [x for x in r if x.startswith("vs_")]:
                v = r[c]
                s += f" | {c}: {v['diff']:+.3f}{'*' if v['sig'] else ''}"
            print(s)
        L.ck("alpha_table", "all", t)
    if what in ("all", "sched"):
        print("=" * 120)
        print("THE SCHEDULE CONTROL (positive = VQE worse)")
        t = sched_table()
        for g in sorted(t, key=float):
            r = t[g]
            for k in ("rmsd_returned", "rand5_coordavg_rmsd"):
                s = f"  rung {g:>4} rho {r['rho']:+.3f} {k:>20} vqe {r[f'vqe|{k}']:.3f} |"
                for a in SCHED_ARMS:
                    v = r[f"vs_{a}_{k}"]
                    s += (f" {a} {r[f'{a}|{k}']:.3f}({v['diff']:+.3f}"
                          f"{'*' if v['sig'] and not v['below_floor'] else ''})")
                print(s)
            print("       distinct: " + " ".join(
                f"{a} {r[f'{a}|n_distinct']:.0f}" for a in ("vqe",) + SCHED_ARMS))
        L.ck("sched_table", "all", t)
    if what in ("all", "curves"):
        print("=" * 120)
        print("T4 ALPHA CURVES vs a one-parameter classical thermostat")
        c = alpha_curves()
        for rung, rec in c.items():
            print(f"  rung {rung}")
            for k, v in rec.items():
                print(f"    {k:>22} VQE range {v['vqe_range']:.3f} | "
                      f"tiltx pearson {v['tiltx_pearson']:+.3f} range {v['tiltx_range']:.3f}"
                      f" | tiltmatch pearson {v['tiltmatch_pearson']:+.3f} "
                      f"range {v['tiltmatch_range']:.3f}")
        L.ck("alpha_curves", "all", c)
        print("\nT2 which axis predicts on the REAL objectives (no oracle knob in any arm)")
        a = real_axis_test()
        for k, v in a.items():
            if isinstance(v, dict):
                print(f"    {k:>28} spearman {v['spearman']:+.3f} pearson {v['pearson']:+.3f}")
        print(f"    n_cells = {a['n_cells']}")
        L.ck("real_axis_test", "all", a)
    if what in ("all", "doing"):
        print("=" * 120)
        print("T5 WHAT THE VQE IS DOING")
        t = doing_table()
        for k in sorted(t, key=float):
            r = t[k]
            print(f"  rung {k:>4} rho {r.get('rho_global', float('nan')):+.3f} "
                  f"apct {r.get('argmin_pct', float('nan')):.3f} | "
                  f"ogap {r.get('objective_gap', float('nan')):.2e} "
                  f"p(argmin) {r.get('p_argmin', float('nan')):.4f} "
                  f"(p0 {r.get('p0_argmin', float('nan')):.2e}) "
                  f"E_p[Epct] {r.get('Ep_E_pct', float('nan')):.4f} "
                  f"(p0 {r.get('Ep0_E_pct', float('nan')):.4f}) "
                  f"H {r.get('entropy_bits', float('nan')):.2f} "
                  f"(H0 {r.get('entropy0_bits', float('nan')):.2f}) "
                  f"ndist {r.get('n_distinct', float('nan')):.0f} "
                  f"mode {r.get('mode_rmsd', float('nan')):.3f} "
                  f"(mode0 {r.get('mode0_rmsd', float('nan')):.3f}) "
                  f"E_p[RMSD] {r.get('mean_rmsd_under_p', float('nan')):.3f} "
                  f"(p0 {r.get('mean_rmsd_under_p0', float('nan')):.3f}) "
                  f"P(<2A) {r.get('pmass_below_2.0', float('nan')):.4f} "
                  f"(p0 {r.get('p0mass_below_2.0', float('nan')):.4f}) "
                  f"ret {r.get('rmsd_returned', float('nan')):.3f} "
                  f"bestseen {r.get('rmsd_best_seen', float('nan')):.3f} "
                  f"certopt {r.get('rmsd_of_certified_optimum', float('nan')):.3f}")
        L.ck("doing_table", "all", t)


if __name__ == "__main__":
    main()
