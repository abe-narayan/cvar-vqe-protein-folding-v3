"""SPRINT 20 / WORKSTREAM B -- the report.  Every number in `agentB_FINDINGS.md` is printed here.

    python -m s20.qb2_report land | opt | link | cvar | enc | q2 | all

`link` is the mandatory one (BRIEF section 4, PREREG section 2 endpoint P1-link): for EVERY
landscape metric, does it predict final Ca-RMSD?  A metric with no predictive relationship to
structural outcome is reported as such and not interpreted further.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from s20 import qb2_lib as L
from s20 import qb2_run as R

KINDS = R.KINDS


def _load(tag):
    p = os.path.join(L.RESULTS, f"qb2_{tag}.json")
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        return json.load(fh)


def _agg_land(d):
    """Per (target, objective) scalar summary of the landscape panel."""
    rows = {}
    for pdb, rec in d.items():
        for k, o in rec["obj"].items():
            sp = o["spec"]
            ln = o["lines"]
            rows.setdefault(k, {})[pdb] = {
                "fold": rec["fold"], "n": rec["n"],
                "gnorm": float(np.nanmean(o["gnorm_std"])),
                "gnorm_cv": float(np.nanstd(o["gnorm_std"]) / max(1e-30, abs(np.nanmean(o["gnorm_std"])))),
                "log10_range": rec["cal"][k]["log10_range"],
                "sd_over_iqr": rec["cal"][k]["sd_over_iqr"],
                "skew": rec["cal"][k]["skew_p99_over_p50"],
                "frac_neg": float(np.nanmean([s["frac_neg"] for s in sp])),
                "log10_cond": float(np.nanmean([np.log10(max(s["cond"], 1.0)) for s in sp])),
                "frac_nearzero": float(np.nanmean([s["frac_nearzero"] for s in sp])),
                "aniso": float(np.nanmean([s["aniso"] for s in sp])),
                "n_localmin": float(np.nanmean([x["n_localmin"] for x in ln])) if ln else np.nan,
                "acorr_len": float(np.nanmean([x["acorr_len_rad"] for x in ln])) if ln else np.nan,
                "tv_over_range": float(np.nanmean([x["tv_over_range"] for x in ln])) if ln else np.nan,
                "barrier": float(np.nanmean(o["barrier"])),
            }
    return rows


METRICS = ["gnorm", "gnorm_cv", "log10_range", "sd_over_iqr", "skew", "frac_neg",
           "log10_cond", "frac_nearzero", "aniso", "n_localmin", "acorr_len",
           "tv_over_range", "barrier"]


def rep_land():
    d = _load("land")
    if not d:
        print("no qb2_land.json"); return {}
    rows = _agg_land(d)
    pdbs = sorted(d)
    print(f"\n=== LANDSCAPE PANEL, n = {len(pdbs)} targets, robust scale (median, IQR/1.349) ===")
    print(f"{'metric':<16}" + "".join(f"{k:>14}" for k in KINDS))
    out = {"n": len(pdbs), "targets": pdbs, "metric": {}}
    for m in METRICS:
        line = f"{m:<16}"
        out["metric"][m] = {}
        for k in KINDS:
            v = np.array([rows[k][p][m] for p in pdbs], float)
            ci = L.boot_mean_ci(v)
            out["metric"][m][k] = ci
            line += f"{ci['mean']:>14.4g}"
        print(line)
    # the AMB vs AMBc contrast: EXACT -- a monotone transform, identical argmin and ranking
    print("\n--- AMB vs AMBc (a strictly monotone reparameterisation of the SAME energy) ---")
    out["AMB_vs_AMBc"] = {}
    for m in METRICS:
        a = np.array([rows["AMB"][p][m] for p in pdbs], float)
        b = np.array([rows["AMBc"][p][m] for p in pdbs], float)
        t = L.paired_ci(a, b, folds=[rows["AMB"][p]["fold"] for p in pdbs])
        out["AMB_vs_AMBc"][m] = t
        print(f"  {m:<16} {t.get('mean', float('nan')):>12.4g}  "
              f"CI [{t.get('ci95',[np.nan,np.nan])[0]:.4g}, {t.get('ci95',[np.nan,np.nan])[1]:.4g}]"
              f"  W/L {t.get('W','-')}/{t.get('L','-')}  sig={t.get('sig')}")
    print("\n--- AMBc vs LEG (is the CONDITIONED AMBER landscape Legacy-like?) ---")
    out["AMBc_vs_LEG"] = {}
    for m in METRICS:
        a = np.array([rows["AMBc"][p][m] for p in pdbs], float)
        b = np.array([rows["LEG"][p][m] for p in pdbs], float)
        t = L.paired_ci(a, b, folds=[rows["AMB"][p]["fold"] for p in pdbs])
        out["AMBc_vs_LEG"][m] = t
        print(f"  {m:<16} {t.get('mean', float('nan')):>12.4g}  "
              f"CI [{t.get('ci95',[np.nan,np.nan])[0]:.4g}, {t.get('ci95',[np.nan,np.nan])[1]:.4g}]"
              f"  W/L {t.get('W','-')}/{t.get('L','-')}  sig={t.get('sig')}")
    ex = [d[p]["amber_exact_maxrel"] for p in pdbs if d[p].get("amber_exact_maxrel") is not None]
    nc = [d[p]["amber_exact_ncmp"] for p in pdbs if d[p].get("amber_exact_ncmp") is not None]
    print(f"\nAMBER bit-exactness gate: max relative difference {max(ex):.3e} "
          f"over {int(np.sum(nc))} comparisons on {len(ex)} targets  (the gate FIRED {int(np.sum(nc))} times)")
    out["amber_exact_maxrel"] = float(max(ex)) if ex else None
    out["amber_exact_ncmp"] = int(np.sum(nc)) if nc else 0
    return out


def _opt_table(tag="opt"):
    d = _load(tag)
    rows = {}
    for pdb, rec in d.items():
        for key, r in rec["rows"].items():
            k, arm, sd_ = key.split("|")
            rows.setdefault((k, arm), {}).setdefault(pdb, []).append(r)
    return d, rows


def rep_opt(tag="opt"):
    d, rows = _opt_table(tag)
    if not d:
        print(f"no qb2_{tag}.json"); return {}
    pdbs = sorted(d)
    folds = [d[p]["fold"] for p in pdbs]
    arms = sorted({a for (_k, a) in rows})
    kinds = sorted({k for (k, _a) in rows})
    print(f"\n=== S1 OPTIMISER BATTERY ({tag}), n = {len(pdbs)} targets, "
          f"B = {d[pdbs[0]]['budget']} evaluations per (objective, arm, seed) ===")
    print("Ca-RMSD of the best-seen configuration (ORACLE, post-hoc), mean over targets and seeds")
    print(f"{'arm':<12}" + "".join(f"{k:>10}" for k in kinds))
    out = {"n": len(pdbs), "targets": pdbs, "budget": d[pdbs[0]]["budget"], "rmsd": {},
           "obj": {}, "vs_bestofN": {}}
    for a in arms:
        line = f"{a:<12}"
        for k in kinds:
            v = np.array([np.mean([x["rmsd_ORACLE"] for x in rows[(k, a)][p]]) for p in pdbs])
            out["rmsd"].setdefault(a, {})[k] = L.boot_mean_ci(v)
            line += f"{v.mean():>10.3f}"
        print(line)
    print("\nbest-seen OBJECTIVE reached (robust-standardised units; lower is better)")
    print(f"{'arm':<12}" + "".join(f"{k:>12}" for k in kinds))
    for a in arms:
        line = f"{a:<12}"
        for k in kinds:
            v = np.array([np.mean([x["best_std"] for x in rows[(k, a)][p]]) for p in pdbs])
            out["obj"].setdefault(a, {})[k] = L.boot_mean_ci(v)
            line += f"{v.mean():>12.4g}"
        print(line)
    if "best_of_N" in arms:
        print("\n--- vs best_of_N (MANDATORY control), paired, negative = the optimiser WINS ---")
        for k in kinds:
            print(f"  [{k}]")
            b = np.array([np.mean([x["rmsd_ORACLE"] for x in rows[(k, "best_of_N")][p]]) for p in pdbs])
            bo = np.array([np.mean([x["best_std"] for x in rows[(k, "best_of_N")][p]]) for p in pdbs])
            for a in arms:
                if a == "best_of_N":
                    continue
                v = np.array([np.mean([x["rmsd_ORACLE"] for x in rows[(k, a)][p]]) for p in pdbs])
                vo = np.array([np.mean([x["best_std"] for x in rows[(k, a)][p]]) for p in pdbs])
                t = L.paired_ci(v, b, folds=folds)
                to = L.paired_ci(vo, bo, folds=folds)
                out["vs_bestofN"].setdefault(k, {})[a] = {"rmsd": t, "obj": to}
                print(f"    {a:<11} RMSD {t['mean']:+7.3f} [{t['ci95'][0]:+.3f},{t['ci95'][1]:+.3f}]"
                      f" W/L {t['W']}/{t['L']} sig={int(t['sig'])} | "
                      f"OBJ {to['mean']:+10.4g} [{to['ci95'][0]:+.4g},{to['ci95'][1]:+.4g}]"
                      f" sig={int(to['sig'])}")
    # the ranking axis, and the reference RMSDs
    if "rank" in d[pdbs[0]]:
        print("\n--- DOES THE OBJECTIVE RANK?  (the same K=500 pool the pipeline uses) ---")
        out["rank"] = {}
        for k in kinds:
            sp = np.array([d[p]["rank"][k]["spearman_E_vs_rmsd_ORACLE"] for p in pdbs])
            am = np.array([d[p]["rank"][k]["argmin_rmsd_ORACLE"] for p in pdbs])
            t1 = np.array([d[p]["rank"][k]["top1pct_mean_rmsd_ORACLE"] for p in pdbs])
            out["rank"][k] = {"spearman": L.boot_mean_ci(sp), "argmin_rmsd": L.boot_mean_ci(am),
                              "top1pct_rmsd": L.boot_mean_ci(t1)}
            print(f"  {k:<5} rho(E, RMSD) {sp.mean():+.3f} [{out['rank'][k]['spearman']['ci95'][0]:+.3f},"
                  f"{out['rank'][k]['spearman']['ci95'][1]:+.3f}]   argmin RMSD {am.mean():.3f}"
                  f"   top-1% mean {t1.mean():.3f}")
        pb = np.array([d[p]["ref"]["pool_best_ORACLE"] for p in pdbs])
        pm = np.array([d[p]["ref"]["pool_mean_ORACLE"] for p in pdbs])
        st = np.array([d[p]["ref"]["start_rmsd_ORACLE"] for p in pdbs])
        out["ref"] = {"pool_best": L.boot_mean_ci(pb), "pool_mean": L.boot_mean_ci(pm),
                      "start": L.boot_mean_ci(st)}
        print(f"  REFERENCE  pool best {pb.mean():.3f}   pool mean {pm.mean():.3f}"
              f"   shared start {st.mean():.3f}")
    # SPSA gradient quality
    cq = {}
    for k in kinds:
        if (k, "spsa") in rows:
            c = [np.mean(x["spsa_cos_fd"]) for p in pdbs for x in rows[(k, "spsa")][p]
                 if x.get("spsa_cos_fd")]
            if c:
                cq[k] = L.boot_mean_ci(np.array(c))
    if cq:
        print("\n--- SPSA GRADIENT QUALITY: cosine(SPSA estimate, central-FD gradient at the SAME point) ---")
        for k, v in cq.items():
            print(f"  {k:<5} {v['mean']:+.4f}  CI [{v['ci95'][0]:+.4f}, {v['ci95'][1]:+.4f}]  n={v['n']}")
        out["spsa_cos"] = cq
    # the sentinel guard: how many times did it fire?
    nf = {k: int(sum(x["n_nonfinite"] for p in pdbs for a in arms for x in rows[(k, a)][p]))
          for k in kinds}
    print(f"\nnon-finite sentinel fired: {nf}   (a guard that never fires is not evidence)")
    out["n_nonfinite"] = nf
    return out


def rep_link():
    """P1-link, MANDATORY.  Does each landscape metric predict final Ca-RMSD?"""
    dl = _load("land")
    do, rows = _opt_table("opt")
    if not dl or not do:
        print("need both qb2_land.json and qb2_opt.json"); return {}
    pdbs = sorted(set(dl) & set(do))
    lrows = _agg_land(dl)
    arms = sorted({a for (_k, a) in rows})
    # THE DIFFICULTY CONTROL (BINDING, s20 LEDGER L4).  On this instrument a RAW rho against
    # Ca-RMSD is largely a target-difficulty measurement -- every circuit-side distributional
    # metric collapsed from |rho| ~ 0.25 to ~0.05 once pool difficulty was partialled out.  So
    # every landscape metric is reported BOTH ways and the PARTIAL is the one that counts.
    diff = np.array([do[p]["ref"]["pool_mean_ORACLE"] for p in pdbs])
    src = "pool_mean_ORACLE (K=500 retrieval pool mean Ca-RMSD)"
    q2 = _load("q2")
    if q2 and all(p in q2 and "pool500|0" in q2[p].get("fam", {}) for p in pdbs):
        diff = np.array([q2[p]["fam"]["pool500|0"]["rmsd_built_ORACLE"] for p in pdbs])
        src = "pool500 realised built-chain Ca-RMSD through the frozen terminal operator"
    print(f"\n=== P1-LINK: does any landscape metric predict final Ca-RMSD?  n = {len(pdbs)} ===")
    print("Spearman across targets, WITHIN each objective, of the metric against the mean")
    print("final Ca-RMSD over every optimiser arm of that objective.  Bootstrap CI on rho.")
    print(f"DIFFICULTY CONTROL partialled out: {src}")
    out = {"n": len(pdbs), "rho": {}, "difficulty_control": src}
    for k in KINDS:
        if not any(kk == k for (kk, _a) in rows):
            continue
        y = np.array([np.mean([np.mean([x["rmsd_ORACLE"] for x in rows[(k, a)][p]])
                               for a in arms if (k, a) in rows]) for p in pdbs])
        print(f"  [{k}]  mean final RMSD {y.mean():.3f}   "
              f"rho(difficulty, final RMSD) {L.spearman(diff, y):+.3f}")
        for m in METRICS:
            x = np.array([lrows[k][p][m] for p in pdbs], float)
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 5:
                continue
            rho = L.spearman(x[ok], y[ok])
            pr = L.partial_spearman(x, y, diff)
            rng = np.random.default_rng(0)
            bs, bp = [], []
            for _ in range(2000):
                idx = rng.integers(0, ok.sum(), ok.sum())
                bs.append(L.spearman(x[ok][idx], y[ok][idx]))
                bp.append(L.partial_spearman(x[ok][idx], y[ok][idx], diff[ok][idx]))
            lo, hi = np.nanpercentile(bs, [2.5, 97.5])
            plo, phh = np.nanpercentile(bp, [2.5, 97.5])
            sig = bool(lo * hi > 0)
            psig = bool(plo * phh > 0)
            out["rho"].setdefault(k, {})[m] = {
                "rho": float(rho), "ci95": [float(lo), float(hi)], "sig": sig,
                "partial": float(pr), "partial_ci95": [float(plo), float(phh)],
                "partial_sig": psig}
            print(f"     {m:<16} raw {rho:+.3f} [{lo:+.3f},{hi:+.3f}]   "
                  f"PARTIAL {pr:+.3f} [{plo:+.3f},{phh:+.3f}]  "
                  f"{'PREDICTS' if psig else 'no'}")
    # MULTIPLICITY, stated rather than hidden: 13 metrics x 4 objectives = 52 tests at 95%,
    # so ~2.6 "PREDICTS" flags are expected by chance alone.  A panel is only informative if
    # the flag COUNT exceeds that, and the surviving metrics must be the same ones across
    # objectives.  Reported both raw and partialled.
    nraw = sum(1 for k in out["rho"] for m in out["rho"][k] if out["rho"][k][m]["sig"])
    npar = sum(1 for k in out["rho"] for m in out["rho"][k] if out["rho"][k][m]["partial_sig"])
    ntot = sum(len(out["rho"][k]) for k in out["rho"])
    out["multiplicity"] = {"n_tests": ntot, "n_sig_raw": nraw, "n_sig_partial": npar,
                           "expected_by_chance_at_95pct": 0.05 * ntot}
    print("")
    print(f"  MULTIPLICITY: {ntot} tests at 95% -> {0.05*ntot:.1f} flags expected by chance. "
          f"Observed: raw {nraw}, PARTIALLED {npar}.")
    surv = {}
    for k in out["rho"]:
        for m in out["rho"][k]:
            if out["rho"][k][m]["partial_sig"]:
                surv.setdefault(m, []).append(k)
    out["multiplicity"]["partial_survivors"] = surv
    print(f"  metrics surviving the partial, and on which objectives: "
          + (", ".join(f"{m}:{'+'.join(v)}" for m, v in surv.items()) if surv else "NONE"))
    # ACROSS objectives: does the objective with the harder landscape do worse?
    print("\n  ACROSS objectives (the pooled question the brief asks):")
    allx, ally, lab = {m: [] for m in METRICS}, [], []
    for k in KINDS:
        if not any(kk == k for (kk, _a) in rows):
            continue
        for p in pdbs:
            ally.append(np.mean([np.mean([x["rmsd_ORACLE"] for x in rows[(k, a)][p]])
                                 for a in arms if (k, a) in rows]))
            lab.append(k)
            for m in METRICS:
                allx[m].append(lrows[k][p][m])
    ally = np.array(ally)
    out["pooled"] = {}
    for m in METRICS:
        x = np.array(allx[m], float)
        ok = np.isfinite(x) & np.isfinite(ally)
        rho = L.spearman(x[ok], ally[ok])
        out["pooled"][m] = float(rho)
        print(f"     {m:<16} rho {rho:+.3f}   (pooled over objectives, NOT target-paired)")
    return out


def rep_cvar():
    d = _load("cvar")
    if not d:
        print("no qb2_cvar.json"); return {}
    pdbs = sorted(d)
    folds = [d[p]["fold"] for p in pdbs]
    rows = {}
    for p in pdbs:
        for key, r in d[p]["rows"].items():
            k, arm, sd_ = key.split("|")
            rows.setdefault((k, arm), {}).setdefault(p, []).append(r)
    kinds = sorted({k for (k, _a) in rows})
    arms = sorted({a for (_k, a) in rows})
    print(f"\n=== S2 CVaR alpha PANEL, n = {len(pdbs)} targets, genuine CVaR-VQE ===")
    hdr = ["rmsd", "best_std", "ess_frac", "latent_entropy_bits", "cos_grad_vs_alpha1",
           "gnorm_mean", "latent_max_prob"]
    out = {"n": len(pdbs), "targets": pdbs, "rows": {}}
    for k in kinds:
        print(f"  [{k}]")
        print("    " + f"{'arm':<16}" + "".join(f"{h:>22}" for h in hdr))
        for a in arms:
            if (k, a) not in rows:
                continue
            vals = []
            for h in hdr:
                key = "rmsd_ORACLE" if h == "rmsd" else h
                v = np.array([np.mean([x.get(key, np.nan) for x in rows[(k, a)][p]]) for p in pdbs])
                vals.append(np.nanmean(v))
                out["rows"].setdefault(k, {}).setdefault(a, {})[h] = L.boot_mean_ci(v)
            print("    " + f"{a:<16}" + "".join(f"{v:>22.4g}" for v in vals))
        # alpha=1 vs alpha=0.05, and every alpha vs the MANDATORY untrained control
        if (k, "vqe_untrained") in rows:
            u = np.array([np.mean([x["rmsd_ORACLE"] for x in rows[(k, "vqe_untrained")][p]])
                          for p in pdbs])
            for a in arms:
                if a == "vqe_untrained":
                    continue
                v = np.array([np.mean([x["rmsd_ORACLE"] for x in rows[(k, a)][p]]) for p in pdbs])
                t = L.paired_ci(v, u, folds=folds)
                out.setdefault("vs_untrained", {}).setdefault(k, {})[a] = t
                print(f"      {a:<14} vs untrained: {t['mean']:+.3f} "
                      f"[{t['ci95'][0]:+.3f},{t['ci95'][1]:+.3f}] W/L {t['W']}/{t['L']} "
                      f"folds {t.get('folds_same_sign','-')}/5 sig={int(t['sig'])}")
    # ---- DERIVE THE OPERATOR BEFORE INTERPRETING ITS STATISTIC (BRIEF section 9).
    # The CVaR weight vector is w_i = -(q_alpha - e_i)/(alpha N) on the tail and 0 elsewhere,
    # so ESS = (sum|w|)^2 / sum w^2 is a functional of the TAIL'S ORDER STATISTICS ALONE.  It
    # contains no reference to the Hamiltonian.  Any LEG-vs-AMB difference in measured CVaR
    # concentration is therefore a difference in their ENERGY SPECTRA passed through one
    # SHARED, Hamiltonian-independent operator -- an identity, not a property of AMBER.
    # Demonstrated numerically: ESS is regressed on the objective's own spectrum shape.
    land = _load("land")
    if land:
        print("")
        print("--- CVaR concentration is a property of the RULE, not of the Hamiltonian ---")
        print("    (ESS depends on the batch's order statistics only; the spectrum is what differs)")
        out["operator"] = {}
        for a in arms:
            xs, ys, ls = [], [], []
            for k in kinds:
                if (k, a) not in rows:
                    continue
                for p in pdbs:
                    if p not in land or k not in land[p]["cal"]:
                        continue
                    sk = land[p]["cal"][k]["skew_p99_over_p50"]
                    ess = np.mean([x.get("ess_frac", np.nan) for x in rows[(k, a)][p]])
                    if np.isfinite(sk) and np.isfinite(ess):
                        xs.append(np.log10(max(sk, 1e-12))); ys.append(ess); ls.append(k)
            if len(xs) > 6:
                r = L.spearman(np.array(xs), np.array(ys))
                out["operator"][a] = {"rho_log10skew_vs_ess": float(r), "n": len(xs)}
                print(f"    {a:<16} rho(log10 spectrum skew, ESS fraction) = {r:+.3f}  n={len(xs)}")
    return out


def rep_enc():
    """S3: theta versus (cos theta, sin theta).  Physically equivalent, geometrically different."""
    do, ro = _opt_table("opt")
    de, re_ = _opt_table("enc")
    if not de:
        print("no qb2_enc.json"); return {}
    pdbs = sorted(set(do) & set(de))
    folds = [do[p]["fold"] for p in pdbs]
    arms = sorted({a for (_k, a) in re_})
    kinds = sorted({k for (k, _a) in re_})
    print(f"\n=== S3 ENCODING: theta vs (cos, sin), n = {len(pdbs)}  (negative = the EMBEDDING wins) ===")
    out = {"n": len(pdbs), "targets": pdbs, "delta": {}}
    for k in kinds:
        print(f"  [{k}]")
        for a in arms:
            if (k, a) not in ro or (k, a) not in re_:
                continue
            th = np.array([np.mean([x["rmsd_ORACLE"] for x in ro[(k, a)][p]]) for p in pdbs])
            em = np.array([np.mean([x["rmsd_ORACLE"] for x in re_[(k, a)][p]]) for p in pdbs])
            to = np.array([np.mean([x["best_std"] for x in ro[(k, a)][p]]) for p in pdbs])
            eo = np.array([np.mean([x["best_std"] for x in re_[(k, a)][p]]) for p in pdbs])
            t = L.paired_ci(em, th, folds=folds)
            t2 = L.paired_ci(eo, to, folds=folds)
            out["delta"].setdefault(k, {})[a] = {"rmsd": t, "obj": t2}
            print(f"    {a:<11} RMSD {t['mean']:+7.3f} [{t['ci95'][0]:+.3f},{t['ci95'][1]:+.3f}]"
                  f" W/L {t['W']}/{t['L']} sig={int(t['sig'])} | OBJ {t2['mean']:+10.4g}"
                  f" [{t2['ci95'][0]:+.4g},{t2['ci95'][1]:+.4g}] sig={int(t2['sig'])}")
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if mode in ("land", "all"):
        out["land"] = rep_land()
    if mode in ("opt", "all"):
        out["opt"] = rep_opt("opt")
    if mode in ("link", "all"):
        out["link"] = rep_link()
    if mode in ("cvar", "all"):
        out["cvar"] = rep_cvar()
    if mode in ("enc", "all"):
        out["enc"] = rep_enc()
    if out:
        L.write(f"qb2_report_{mode}", out, complete=True)
        print(f"\nwritten s20/results/qb2_report_{mode}.json")


if __name__ == "__main__":
    main()


def rep_decouple(tag="opt"):
    """THE DECISIVE TEST (BRIEF section 4, PREREG endpoint P1-link).

    Across every (objective, arm) cell, is BEATING `best_of_N` ON THE OBJECTIVE associated with
    beating it on Ca-RMSD?  Two forms are reported:

      (i)  the SIGN TABLE over the 36 (objective, arm) cells -- how many cells improve the
           objective and the RMSD together, versus how many trade one for the other, against
           the 50/50 null a coin would give;
      (ii) the WITHIN-CELL, TARGET-LEVEL correlation: for each cell, Spearman across targets of
           (objective gained vs best_of_N) against (RMSD gained vs best_of_N).  A positive rho
           means the targets where the arm optimised harder are the targets where it built a
           better structure.  This is the relationship the brief demands, measured directly, and
           it does not depend on any landscape metric at all.
    """
    d, rows = _opt_table(tag)
    if not d:
        print(f"no qb2_{tag}.json"); return {}
    pdbs = sorted(d)
    folds = [d[p]["fold"] for p in pdbs]
    arms = sorted({a for (_k, a) in rows} - {"best_of_N"})
    kinds = sorted({k for (k, _a) in rows})
    print(f"\n=== DOES OPTIMISING HARDER BUILD A BETTER STRUCTURE?  n = {len(pdbs)} targets, "
          f"{len(kinds)}x{len(arms)} = {len(kinds)*len(arms)} cells ===")
    out = {"n": len(pdbs), "cells": {}, "sign_table": {}}
    both, trade, neither, rhos = 0, 0, 0, []
    for k in kinds:
        for a in arms:
            if (k, a) not in rows or (k, "best_of_N") not in rows:
                continue
            dobj = np.array([np.mean([x["best_std"] for x in rows[(k, a)][p]])
                             - np.mean([x["best_std"] for x in rows[(k, "best_of_N")][p]])
                             for p in pdbs])
            drms = np.array([np.mean([x["rmsd_ORACLE"] for x in rows[(k, a)][p]])
                             - np.mean([x["rmsd_ORACLE"] for x in rows[(k, "best_of_N")][p]])
                             for p in pdbs])
            to = L.paired_ci(dobj, np.zeros_like(dobj), folds=folds)
            tr = L.paired_ci(drms, np.zeros_like(drms), folds=folds)
            rho = L.spearman(dobj, drms)
            rhos.append(rho)
            sig_o, sig_r = to["sig"], tr["sig"]
            if sig_o and sig_r:
                if np.sign(to["mean"]) == np.sign(tr["mean"]):
                    both += 1
                else:
                    trade += 1
            else:
                neither += 1
            out["cells"][f"{k}|{a}"] = {"d_obj": to, "d_rmsd": tr, "rho_within": float(rho)}
    print(f"  cells where objective and RMSD BOTH move significantly, SAME direction: {both}")
    print(f"  cells where they BOTH move significantly, OPPOSITE directions (a TRADE): {trade}")
    print(f"  cells where at most one moves significantly:                            {neither}")
    out["sign_table"] = {"both_same": both, "both_opposite": trade, "at_most_one": neither}
    r = np.array(rhos, float)
    ci = L.boot_mean_ci(r)
    out["rho_within_cell"] = ci
    out["rho_within_cell_positive"] = int((r > 0).sum())
    out["rho_within_cell_n"] = int(len(r))
    print(f"\n  WITHIN-CELL target-level rho(objective gained, RMSD gained), over "
          f"{len(r)} cells:")
    print(f"    mean {ci['mean']:+.4f}  median {ci['median']:+.4f}  "
          f"CI [{ci['ci95'][0]:+.4f}, {ci['ci95'][1]:+.4f}]   positive in "
          f"{int((r>0).sum())}/{len(r)} cells")
    print("    (positive = optimising harder on a target builds a better structure on that target)")
    L.write(f"qb2_report_decouple_{tag}", out, complete=True)
    return out
