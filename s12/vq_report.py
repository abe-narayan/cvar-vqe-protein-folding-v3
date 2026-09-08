"""s12 QUANTUM-ROLE -- pooled tables and paired statistics over the vq_* result files."""
from __future__ import annotations
import os, sys, json, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

RES = os.path.join(ROOT, "s12", "results")


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


def split(rows, key):
    v = np.array([r[key] for r in rows], float)
    pdb = [r["pdb"] for r in rows]
    f18 = np.array([p in I.FAIL18 for p in pdb])
    return v, f18


def s1_table(fam):
    d = load(f"vq_stage1_{fam}.json")
    rows = d["rows"]
    pdb = [r["pdb"] for r in rows]
    f18 = np.array([p in I.FAIL18 for p in pdb])
    st = {k: np.array([r["structure"][k] for r in rows], float)
          for k in rows[0]["structure"]}
    out = {"family": fam, "n_targets": len(rows), "cfg": d["cfg"],
           "n_configs": float(st["n_configs"][0]),
           "interaction_share_mean": float(st["interaction_share"].mean()),
           "interaction_share_median": float(np.median(st["interaction_share"])),
           "interaction_share_max": float(st["interaction_share"].max()),
           "anova_residual_max": float(st["anova_residual"].max()),
           "n_local_minima_mean": float(st["n_local_minima"].mean()),
           "n_local_minima_max": float(st["n_local_minima"].max()),
           "frac_multimodal": float((st["n_local_minima"] > 1).mean())}
    names = sorted({k for r in rows for k in r["solvers"]})
    opt = np.array([r["solvers"]["exhaustive"]["value"] for r in rows], float)
    solv = {}
    for nm in names:
        vals, ev, wa = [], [], []
        for r in rows:
            s = r["solvers"].get(nm, {})
            if "value" not in s:
                vals.append(np.nan); ev.append(np.nan); wa.append(np.nan); continue
            vals.append(s["value"]); ev.append(s.get("evals", np.nan)); wa.append(s.get("wall", np.nan))
        vals = np.array(vals, float)
        gap = vals - opt
        solv[nm] = dict(mean_value=float(np.nanmean(vals)),
                        mean_gap=float(np.nanmean(gap)),
                        frac_optimal=float(np.nanmean(gap <= 1e-9)),
                        mean_evals=float(np.nanmean(ev)),
                        median_wall_ms=float(np.nanmedian(wa) * 1e3))
    out["solvers"] = solv
    rm = {}
    keys = sorted({k for r in rows for k in r["rmsd"]})
    for k in keys:
        v = np.array([r["rmsd"].get(k, np.nan) for r in rows], float)
        rm[k] = dict(mean=float(np.nanmean(v)),
                     FAIL18=float(np.nanmean(v[f18])),
                     other108=float(np.nanmean(v[~f18])),
                     frac_under_2=float(np.nanmean(v < 2.0)))
    out["rmsd"] = rm
    pct = np.array([r["rmsd_percentile_of_optimum"] for r in rows], float)
    cor = np.array([r["corr_E_rmsd"] for r in rows], float)
    out["rmsd_percentile_of_exact_optimum"] = dict(mean=float(pct.mean()),
                                                   FAIL18=float(pct[f18].mean()),
                                                   other108=float(pct[~f18].mean()))
    out["corr_energy_rmsd"] = dict(mean=float(np.nanmean(cor)),
                                   FAIL18=float(np.nanmean(cor[f18])),
                                   other108=float(np.nanmean(cor[~f18])))
    # paired: exact optimum vs anchor, and vs greedy_ls
    a = np.array([r["rmsd"]["exhaustive"] for r in rows], float)
    b = np.array([r["rmsd"]["anchor_fit_ca"] for r in rows], float)
    folds = np.array([r["fold"] for r in rows])
    out["paired_exactopt_vs_anchor"] = I.paired(a, b, folds=folds, names=pdb)
    g = np.array([r["rmsd"].get("greedy_ls", np.nan) for r in rows], float)
    out["paired_exactopt_vs_greedyls"] = I.paired(a, g, folds=folds, names=pdb)
    return out


def s3_table(fam):
    d = load(f"vq_stage3_{fam}.json")
    rows = d["rows"]
    pdb = [r["pdb"] for r in rows]
    f18 = np.array([p in I.FAIL18 for p in pdb])
    opt = np.array([r["opt"] for r in rows], float)
    out = {"family": fam, "n_targets": len(rows), "alpha": d["alpha"],
           "shots": d["shots"], "layers": d["layers"], "restarts": d["restarts"],
           "dim": rows[0]["dim"], "n_qubits": rows[0]["n_qubits"], "budgets": {}}
    for B in sorted(rows[0]["budgets"], key=int):
        e = {}
        for nm in ("vqe", "random", "anneal", "greedy_ls"):
            vals, ev, wa = [], [], []
            for r in rows:
                x = r["budgets"][B].get(nm, {})
                v = x.get("best_energy", x.get("best_seen_energy", np.nan))
                vals.append(v if v is not None else np.nan)
                ev.append(x.get("n_energy_evaluations", np.nan))
                wa.append(x.get("wall", np.nan))
            vals = np.array(vals, float)
            gap = vals - opt
            e[nm] = dict(mean_gap=float(np.nanmean(gap)),
                         median_gap=float(np.nanmedian(gap)),
                         frac_optimal=float(np.nanmean(gap <= 1e-9)),
                         mean_evals_used=float(np.nanmean(ev)),
                         median_wall_ms=float(np.nanmedian(wa) * 1e3))
            if nm == "vqe":
                e[nm]["mean_spsa_iters"] = float(np.nanmean(
                    [r["budgets"][B]["vqe"].get("spsa_iters", np.nan) for r in rows]))
                e[nm]["mean_entropy_bits"] = float(np.nanmean(
                    [r["budgets"][B]["vqe"].get("entropy_bits", np.nan) for r in rows]))
        # paired VQE vs each classical, on the energy gap
        gv = np.array([r["budgets"][B]["vqe"].get("best_energy", np.nan) for r in rows], float) - opt
        for nm in ("random", "anneal", "greedy_ls"):
            gc = np.array([r["budgets"][B][nm]["best_energy"] for r in rows], float) - opt
            ok = np.isfinite(gv) & np.isfinite(gc)
            e[f"paired_vqe_minus_{nm}"] = I.paired(gv[ok], gc[ok],
                                                   folds=np.array([r["fold"] for r in rows])[ok])
        out["budgets"][B] = e
    return out


def s2_table(fam):
    d = load(f"vq_stage2_{fam}.json")
    rows = d["rows"]
    combos = {}
    for r in rows:
        for a in r["arms"]:
            key = (a["alpha"], a["T"])
            combos.setdefault(key, []).append((r, a))
    out = {"family": fam, "n_targets": len(rows), "iters": d["iters"],
           "layers": d["layers"], "n_qubits": rows[0]["n_qubits"],
           "dim": rows[0]["dim"], "arms": []}
    for (a, T), lst in sorted(combos.items()):
        g = np.array([x["gap_argmax"] for _, x in lst], float)
        gs = np.array([x["gap_best_in_support"] for _, x in lst], float)
        H = np.array([x["entropy_bits"] for _, x in lst], float)
        ra = np.array([x["rmsd_argmax"] if x["rmsd_argmax"] is not None else np.nan for _, x in lst], float)
        rp = np.array([x["rmsd_p_weighted"] for _, x in lst], float)
        rs = np.array([x["rmsd_best_in_support"] if x["rmsd_best_in_support"] is not None else np.nan for _, x in lst], float)
        fm = np.array([x["feasible_mass"] for _, x in lst], float)
        rg = np.array([x.get("rmsd_gibbs_matched", np.nan) for _, x in lst], float)
        cv = np.array([x.get("cover8_vqe") if x.get("cover8_vqe") is not None else np.nan for _, x in lst], float)
        cg = np.array([x.get("cover8_gibbs") if x.get("cover8_gibbs") is not None else np.nan for _, x in lst], float)
        ce = np.array([x.get("cover8_lowest_energy") if x.get("cover8_lowest_energy") is not None else np.nan for _, x in lst], float)
        out["arms"].append(dict(alpha=a, T=T,
                                mean_gap_argmax=float(np.nanmean(g)),
                                frac_argmax_optimal=float(np.nanmean(g <= 1e-9)),
                                mean_gap_best_in_support=float(np.nanmean(gs)),
                                frac_support_optimal=float(np.nanmean(gs <= 1e-9)),
                                mean_entropy_bits=float(np.nanmean(H)),
                                max_entropy_bits=float(rows[0]["n_qubits"]),
                                mean_feasible_mass=float(np.nanmean(fm)),
                                rmsd_argmax=float(np.nanmean(ra)),
                                rmsd_best_in_support=float(np.nanmean(rs)),
                                rmsd_p_weighted=float(np.nanmean(rp)),
                                rmsd_gibbs_matched=float(np.nanmean(rg)),
                                cover8_vqe=float(np.nanmean(cv)),
                                cover8_gibbs=float(np.nanmean(cg)),
                                cover8_lowest_energy=float(np.nanmean(ce)),
                                paired_pw_minus_gibbs=I.paired(rp[np.isfinite(rp) & np.isfinite(rg)],
                                                               rg[np.isfinite(rp) & np.isfinite(rg)]),
                                mean_wall=float(np.nanmean([x["wall"] for _, x in lst]))))
    ref = dict(exact_opt_rmsd=float(np.mean([r["exact_opt_rmsd"] for r in rows])),
               oracle_best_rmsd=float(np.mean([r["oracle_best_rmsd"] for r in rows])),
               anchor_rmsd=float(np.mean([r["anchor_rmsd"] for r in rows])))
    out["reference_rmsd"] = ref
    out["grad_audit_mean"] = {k: float(np.mean([r["grad_audit"][k] for r in rows]))
                              for k in rows[0]["grad_audit"]}
    return out


def e2e_table(fam):
    d = load(f"vq_e2e_{fam}.json")
    rows = d["rows"]
    pdb = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    f18 = np.array([p in I.FAIL18 for p in pdb])
    names = sorted({k for r in rows for k in r["proj"]})
    out = {"family": fam, "budget": d["budget"], "n_targets": len(rows), "arms": {}}
    ref = np.array([r["proj"]["anchor"] for r in rows], float)
    for nm in names:
        raw = np.array([r["raw"].get(nm, np.nan) for r in rows], float)
        prj = np.array([r["proj"].get(nm, np.nan) for r in rows], float)
        e = dict(raw_mean=float(np.nanmean(raw)),
                 proj_mean=float(np.nanmean(prj)),
                 proj_FAIL18=float(np.nanmean(prj[f18])),
                 proj_other108=float(np.nanmean(prj[~f18])),
                 frac_under_2=float(np.nanmean(prj < 2.0)))
        ok = np.isfinite(prj) & np.isfinite(ref)
        if nm != "anchor":
            e["paired_vs_anchor"] = I.paired(prj[ok], ref[ok], folds=folds[ok],
                                             names=[pdb[i] for i in np.where(ok)[0]])
        out["arms"][nm] = e
    # vqe vs exact and vs greedy_ls
    for a, b in (("vqe", "exact"), ("vqe", "greedy_ls"), ("vqe_cvar_tail_avg", "vqe"),
                 ("vqe_cvar_tail_avg", "anchor")):
        x = np.array([r["proj"].get(a, np.nan) for r in rows], float)
        y = np.array([r["proj"].get(b, np.nan) for r in rows], float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() > 5:
            out[f"paired_{a}_vs_{b}"] = I.paired(x[ok], y[ok], folds=folds[ok])
    out["vqe_energy_gap"] = float(np.nanmean([r["vqe_best_energy"] - r["opt"] for r in rows]))
    out["vqe_entropy_bits"] = float(np.nanmean([r["vqe_entropy_bits"] for r in rows]))
    return out


def b_table(name="vq_stage3_B.json"):
    d = load(name)
    rows = d["rows"]
    pdb = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    f18 = np.array([p in I.FAIL18 for p in pdb])
    opt = np.array([r["opt"] for r in rows], float)
    out = {"problem": "B_subset", "n_targets": len(rows), "validation": d["validation"],
           "dim": rows[0]["dim"], "n_qubits": rows[0]["n_qubits"],
           "opt_rmsd": float(np.mean([r["opt_rmsd"] for r in rows])),
           "oracle_best_rmsd": float(np.mean([r["oracle_best_rmsd"] for r in rows])),
           "avg75_rmsd": float(np.mean([r["avg75_rmsd"] for r in rows])),
           "uniform16_rmsd": float(np.mean([r["uniform16_rmsd"] for r in rows])),
           "rmsd_percentile_of_optimum": float(np.mean([r["rmsd_percentile_of_optimum"] for r in rows])),
           "corr_E_rmsd": float(np.nanmean([r["corr_E_rmsd"] for r in rows])),
           "budgets": {}}
    for B in sorted(rows[0]["budgets"], key=int):
        e = {}
        for nm in ("vqe", "random", "anneal", "greedy_ls"):
            v = np.array([r["budgets"][B][nm].get("best_energy", np.nan) for r in rows], float)
            rm = np.array([r["budgets"][B][nm].get("rmsd") if r["budgets"][B][nm].get("rmsd") is not None else np.nan for r in rows], float)
            wa = np.array([r["budgets"][B][nm].get("wall", np.nan) for r in rows], float)
            e[nm] = dict(mean_gap=float(np.nanmean(v - opt)),
                         frac_optimal=float(np.nanmean((v - opt) <= 1e-9)),
                         rmsd=float(np.nanmean(rm)),
                         rmsd_FAIL18=float(np.nanmean(rm[f18])),
                         rmsd_other108=float(np.nanmean(rm[~f18])),
                         median_wall_s=float(np.nanmedian(wa)))
            if nm == "vqe":
                e[nm]["spsa_iters"] = float(np.nanmean([r["budgets"][B]["vqe"].get("spsa_iters", np.nan) for r in rows]))
                e[nm]["entropy_bits"] = float(np.nanmean([r["budgets"][B]["vqe"].get("entropy_bits", np.nan) for r in rows]))
        gv = np.array([r["budgets"][B]["vqe"].get("best_energy", np.nan) for r in rows], float) - opt
        for nm in ("random", "anneal", "greedy_ls"):
            gc = np.array([r["budgets"][B][nm]["best_energy"] for r in rows], float) - opt
            ok = np.isfinite(gv) & np.isfinite(gc)
            e[f"paired_vqe_minus_{nm}"] = I.paired(gv[ok], gc[ok], folds=folds[ok])
        out["budgets"][B] = e
    return out


if __name__ == "__main__":
    what = sys.argv[1]
    fams = sys.argv[2:]
    if what == "B":
        res = {"B": b_table()}
        print(json.dumps(res, indent=1)); I.write("vq_report_B", res); sys.exit(0)
    fn = {"s1": s1_table, "s2": s2_table, "s3": s3_table, "e2e": e2e_table}[what]
    res = {}
    for f in fams:
        try:
            res[f] = fn(f)
        except FileNotFoundError:
            print(f"  (no result file for {f})", file=sys.stderr)
    print(json.dumps(res, indent=1))
    I.write(f"vq_report_{what}", res)
