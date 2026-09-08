"""SPRINT 14, ENER-3 -- component decomposition, complementarity (H7), and Pareto structure.

Three questions, in the order that makes the third one worth asking:

  1. DECOMPOSITION.  What fraction of each total's variance is one term?  The sprint brief
     records AMBER `nonbonded` at covariance share 1.000 and Legacy `steric` at 0.955 --
     "both models are a steric potential plus rounding error".  Confirmed or refuted here
     on the FULL 262,144-configuration enumeration for Legacy (no sampling error at all)
     and on the uniform AMBER stratum for AMBER.

  2. COMPLEMENTARITY (H7).  Raw correlation between two models is the wrong measure: two
     models that both merely track compactness look complementary if compactness is noisy.
     The right measure is the TRUTH-PARTIALLED error correlation -- regress each model's
     rank on the rank of the true RMSD and correlate the residuals.  Then price any
     proposed fusion against the project's established arithmetic: the gain from fusing a
     strong channel with a weak one goes as the SQUARE of the weaker channel's skill.

  3. PARETO.  Do the objectives actually CONFLICT?  If every pair of axes is positively
     rank-correlated there is no frontier worth having and multi-objective machinery is
     pure complexity.  Measure the conflict before proposing the structure.

Plus: a leave-one-target-out linear REWEIGHTING of the eleven Legacy terms, to test whether
the default weights destroy signal that the components carry.  In-sample is the ceiling;
LOTO is the honest number; both are reported and the gap is the overfit.

    python -m s14.ener_complement
"""
from __future__ import annotations

import itertools

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E
from s14.ener_matrix import ensemble


# ------------------------------------------------------------------ 1. decomposition
def decompose(pdb):
    """Covariance share of each term in its own total. Shares sum to 1 by construction."""
    from core import energy as et
    z = E.enum(pdb)
    out = {}
    # LEGACY -- full enumeration, exact
    tot = z.legacy
    v = tot.var()
    leg = {}
    for t in E.LEG_TERMS:
        w = float(et.DEFAULT_WEIGHTS.get(t, 0.0))
        leg[t] = dict(weight=w,
                      share=float(w * np.cov(z.leg[t], tot, ddof=0)[0, 1] / v),
                      sd=float(z.leg[t].std()), w_sd=float(abs(w) * z.leg[t].std()),
                      rho_with_total=E.spearman(z.leg[t], tot),
                      rho_with_rmsd=E.spearman(z.leg[t], z.rmsd),
                      frac_zero=float((z.leg[t] == 0).mean()))
    out["legacy"] = dict(n=int(z.B), terms=leg,
                         top_share=max(leg.items(), key=lambda kv: kv[1]["share"]))
    # AMBER -- uniform stratum only (unconditioned)
    m = z.uniform_mask
    at = z.amber_total[m]
    va = at.var()
    amb = {}
    for t in E.AMB_TERMS:
        x = z.amb[t][m]
        amb[t] = dict(share=float(np.cov(x, at, ddof=0)[0, 1] / va) if va > 0 else np.nan,
                      sd=float(x.std()), rho_with_total=E.spearman(x, at),
                      rho_with_rmsd=E.spearman(x, z.rmsd[z.amber_idx[m]]),
                      n_distinct=int(len(np.unique(x))))
    out["amber"] = dict(n=int(m.sum()), terms=amb,
                        top_share=max(amb.items(), key=lambda kv: kv[1]["share"]))
    # the same on a MONOTONE-CONDITIONED AMBER: the raw variance is a delta spike, so the
    # raw covariance share is an artefact of the worst clash. rank-normalise and redo.
    ar = E.rank_norm(at)
    amb_r = {}
    for t in E.AMB_TERMS:
        x = E.rank_norm(z.amb[t][m])
        amb_r[t] = dict(share=float(np.cov(x, ar, ddof=0)[0, 1] / ar.var()),
                        rho_with_total=E.spearman(x, ar))
    out["amber_rank"] = dict(terms=amb_r)
    return out


# ------------------------------------------------------------------ 2. complementarity
def _resid_on_truth(x, R):
    """Residual of rank(x) after removing rank(R). The truth-partialled error signal."""
    a = E.rank_norm(x); b = E.rank_norm(R)
    a = a - a.mean(); b = b - b.mean()
    if b.var() == 0:
        return a
    return a - (a @ b) / (b @ b) * b


def complementarity(pdb, channels=("legacy", "amber", "prior", "leg_steric",
                                   "amb_nonbonded", "legacy_nosteric", "leg_torsion",
                                   "rg")):
    d = ensemble(pdb)
    R = d["rmsd"]
    raw, part, skill = {}, {}, {}
    for c in channels:
        skill[c] = dict(rho=E.spearman(d[c], R),
                        pair=E.pair_accuracy(d[c], R, 200000,
                                             np.random.default_rng(E.SEED)))
    res = {c: _resid_on_truth(d[c], R) for c in channels}
    rk = {c: E.rank_norm(d[c]) for c in channels}
    for a, b in itertools.combinations(channels, 2):
        raw[f"{a}|{b}"] = float(np.corrcoef(rk[a], rk[b])[0, 1])
        part[f"{a}|{b}"] = float(np.corrcoef(res[a], res[b])[0, 1])
    return dict(pdb=pdb, skill=skill, raw_corr=raw, partialled_err_corr=part)


def price_fusion(rho_a, rho_b, raw_corr):
    """Optimal linear fusion of two channels -- the multiple correlation.

        R^2 = (r_a^2 + r_b^2 - 2 r_a r_b r_ab) / (1 - r_ab^2)

    r_ab is the RAW correlation between the two channels (this is the standard multiple-
    correlation identity; using the truth-partialled correlation here is an error and was
    one in the first version of this module).  All correlations are SIGNED and computed on
    ranks, so they are Spearman.

    The GAIN over the stronger channel alone is second order in the weaker channel's
    unique skill -- which is the project's established arithmetic:
        R - |r_a|  ~  (r_b - r_ab r_a)^2 / (2 |r_a| (1 - r_ab^2))   for small unique skill.
    """
    ra, rb, r = float(rho_a), float(rho_b), float(raw_corr)
    if abs(ra) < abs(rb):
        ra, rb = rb, ra
    if abs(r) >= 0.999 or not np.isfinite(r):
        return dict(rho_fused=abs(ra), gain=0.0, r_strong=ra, r_weak=rb, r_ab=r,
                    unique_skill_weak=0.0, second_order=0.0, note="channels collinear")
    den = 1 - r ** 2
    f = float(np.sqrt(max((ra ** 2 + rb ** 2 - 2 * ra * rb * r) / den, 0.0)))
    return dict(rho_fused=f, gain=float(f - abs(ra)), r_strong=ra, r_weak=rb, r_ab=r,
                unique_skill_weak=float(rb - r * ra),
                second_order=float((rb - r * ra) ** 2 / (2 * max(abs(ra), 1e-9) * den)))


# ------------------------------------------------------------------ 3. reweighting
def reweight(loto=True):
    """Leave-one-target-out linear reweighting of the eleven Legacy terms.

    Fits rank(RMSD) ~ sum_t beta_t rank(term_t) on eight targets, evaluates on the ninth.
    In-sample is the CEILING of what a linear reweighting could do; LOTO is the honest
    number.  This is a REWEIGHTING OF PHYSICAL TERMS, not a learned objective -- the
    learned-objective workstream (OBJ) owns that.
    """
    data = {p: ensemble(p) for p in E.ENUM_TARGETS}
    feats = ["leg_" + t for t in E.LEG_TERMS]
    rows = []
    for held in E.ENUM_TARGETS:
        tr = [p for p in E.ENUM_TARGETS if p != held]
        X = np.vstack([np.column_stack([E.rank_norm(data[p][f]) for f in feats])
                       for p in tr])
        y = np.concatenate([E.rank_norm(data[p]["rmsd"]) for p in tr])
        X = np.column_stack([X, np.ones(len(X))])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        d = data[held]
        Xh = np.column_stack([np.column_stack([E.rank_norm(d[f]) for f in feats]),
                              np.ones(d["n"])])
        pred = Xh @ beta
        rows.append(dict(pdb=held, fold=int(d["_z"].fold),
                         rho=E.spearman(pred, d["rmsd"]),
                         rho_decile=E.decile_rho(pred, d["rmsd"]),
                         sel_rmsd=E.argmin_rmsd(pred, d["rmsd"]),
                         pool_mean=float(d["rmsd"].mean()),
                         legacy_sel=E.argmin_rmsd(d["legacy"], d["rmsd"]),
                         beta={f: float(b) for f, b in zip(feats, beta[:-1])}))
    # in-sample ceiling
    X = np.vstack([np.column_stack([E.rank_norm(data[p][f]) for f in feats])
                   for p in E.ENUM_TARGETS])
    y = np.concatenate([E.rank_norm(data[p]["rmsd"]) for p in E.ENUM_TARGETS])
    X = np.column_stack([X, np.ones(len(X))])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    ins = []
    for p in E.ENUM_TARGETS:
        d = data[p]
        Xh = np.column_stack([np.column_stack([E.rank_norm(d[f]) for f in feats]),
                              np.ones(d["n"])])
        pred = Xh @ beta
        ins.append(dict(pdb=p, rho=E.spearman(pred, d["rmsd"]),
                        rho_decile=E.decile_rho(pred, d["rmsd"]),
                        sel_rmsd=E.argmin_rmsd(pred, d["rmsd"])))
    sel = np.array([r["sel_rmsd"] for r in rows])
    base = np.array([r["pool_mean"] for r in rows])
    leg = np.array([r["legacy_sel"] for r in rows])
    return dict(loto=rows, in_sample=ins,
                loto_vs_random=I.paired(sel, base, folds=[r["fold"] for r in rows]),
                loto_vs_legacy=I.paired(sel, leg, folds=[r["fold"] for r in rows]),
                loto_mean_rho=float(np.mean([r["rho"] for r in rows])),
                in_sample_mean_rho=float(np.mean([r["rho"] for r in ins])),
                loto_mean_decile=float(np.mean([r["rho_decile"] for r in rows])),
                in_sample_mean_decile=float(np.mean([r["rho_decile"] for r in ins])),
                in_sample_sel=float(np.mean([r["sel_rmsd"] for r in ins])))


# ------------------------------------------------------------------ 4. Pareto
PARETO_AXES = ["prior", "legacy", "amber", "rg", "hbond", "helix"]


def pareto(pdb):
    """Is there genuine conflict between the axes, and does the frontier hold anything?"""
    d = ensemble(pdb)
    z = d["_z"]
    R = d["rmsd"]
    ax = {"prior": d["prior"], "legacy": d["legacy"], "amber": d["amber"], "rg": d["rg"],
          "hbond": d["leg_hbond_local"] + d["leg_hbond_longrange"]}
    # secondary structure: helix fraction from the decoded torsions (native-free)
    S = z.states(d["idx"]).astype(int)
    rows = np.arange(z.n)
    phi = np.degrees(z.PHI[rows[None, :], S]); psi = np.degrees(z.PSI[rows[None, :], S])
    hel = ((phi > -100) & (phi < -30) & (psi > -80) & (psi < -5)).mean(1)
    ax["helix"] = -hel                       # minimise -> more helical
    names = list(ax)
    M = np.column_stack([E.rank_norm(ax[a]) for a in names])
    conflict = {}
    for a, b in itertools.combinations(range(len(names)), 2):
        conflict[f"{names[a]}|{names[b]}"] = float(E.spearman(M[:, a], M[:, b]))
    front = _pareto_mask(M)
    # conflict measured properly: an axis pair conflicts when rho < 0
    neg = [v for v in conflict.values() if v < -0.05]
    return dict(pdb=pdb, n=int(d["n"]), axes=names, conflict=conflict,
                n_negative_pairs=len(neg), mean_abs_conflict=float(np.mean(
                    [abs(v) for v in conflict.values()])),
                front_size=int(front.sum()), front_frac=float(front.mean()),
                front_best_rmsd=float(R[front].min()),
                front_mean_rmsd=float(R[front].mean()),
                pool_best_rmsd=float(R.min()), pool_mean_rmsd=float(R.mean()),
                legacy_argmin_rmsd=E.argmin_rmsd(d["legacy"], R),
                amber_argmin_rmsd=E.argmin_rmsd(d["amber"], R),
                # does the frontier RETAIN the good structures a scalar discards?
                frac_of_best1pct_on_front=float(
                    front[np.argsort(R)[:max(len(R) // 100, 1)]].mean()),
                # random-set control of the same size
                null_best_rmsd=float(np.mean([
                    R[np.random.default_rng(E.SEED + s).choice(
                        len(R), int(front.sum()), replace=False)].min()
                    for s in range(64)])))


def _pareto_mask(M):
    """Non-dominated rows of M (minimising every column)."""
    n = len(M)
    keep = np.ones(n, bool)
    order = np.argsort(M.sum(1))
    for i in order:
        if not keep[i]:
            continue
        dom = np.all(M <= M[i], axis=1) & np.any(M < M[i], axis=1)
        keep[dom] = False
    return keep


# ------------------------------------------------------------------ main
def main():
    dec = {p: decompose(p) for p in E.ENUM_TARGETS}
    print("=== 1. COVARIANCE SHARE OF EACH TERM IN ITS OWN TOTAL ===")
    print("LEGACY (full 262,144-config enumeration, exact)")
    print(f"{'term':22s} {'weight':>8s} {'share':>8s} {'|w|*sd':>9s} {'rho_tot':>8s} "
          f"{'rho_rmsd':>9s} {'frac=0':>7s}")
    for t in E.LEG_TERMS:
        r = [dec[p]["legacy"]["terms"][t] for p in E.ENUM_TARGETS]
        nfin = int(np.isfinite([x["rho_with_total"] for x in r]).sum())
        note = f"   (constant on {len(r)-nfin}/{len(r)} targets)" if nfin < len(r) else ""
        with np.errstate(all="ignore"):
            print(f"{t:22s} {r[0]['weight']:8.3f} "
                  f"{np.nanmean([x['share'] for x in r]):8.3f} "
                  f"{np.nanmean([x['w_sd'] for x in r]):9.3f} "
                  f"{np.nanmean([x['rho_with_total'] for x in r]):+8.3f} "
                  f"{np.nanmean([x['rho_with_rmsd'] for x in r]):+9.3f} "
                  f"{np.nanmean([x['frac_zero'] for x in r]):7.3f}{note}")
    print("\nAMBER (uniform stratum only)")
    print(f"{'term':22s} {'share':>8s} {'share_rank':>11s} {'rho_tot':>8s} "
          f"{'rho_rmsd':>9s} {'distinct':>9s}")
    for t in E.AMB_TERMS:
        r = [dec[p]["amber"]["terms"][t] for p in E.ENUM_TARGETS]
        rr = [dec[p]["amber_rank"]["terms"][t] for p in E.ENUM_TARGETS]
        print(f"{t:22s} {np.nanmean([x['share'] for x in r]):8.3f} "
              f"{np.nanmean([x['share'] for x in rr]):11.3f} "
              f"{np.nanmean([x['rho_with_total'] for x in r]):+8.3f} "
              f"{np.nanmean([x['rho_with_rmsd'] for x in r]):+9.3f} "
              f"{np.mean([x['n_distinct'] for x in r]):9.0f}")

    comp = [complementarity(p) for p in E.ENUM_TARGETS]
    print("\n=== 2. COMPLEMENTARITY: raw rank corr vs TRUTH-PARTIALLED error corr ===")
    keys = list(comp[0]["raw_corr"])
    print(f"{'pair':38s} {'raw':>8s} {'partialled':>11s}")
    for k in keys:
        print(f"{k:38s} {np.mean([c['raw_corr'][k] for c in comp]):+8.3f} "
              f"{np.mean([c['partialled_err_corr'][k] for c in comp]):+11.3f}")
    print("\nchannel skill (uniform ensemble):")
    for c in comp[0]["skill"]:
        print(f"  {c:20s} rho {np.mean([x['skill'][c]['rho'] for x in comp]):+.3f}"
              f"   pair {np.mean([x['skill'][c]['pair'] for x in comp]):.3f}")
    fus = {}
    for k in keys:
        a, b = k.split("|")
        # per target, then average the PRICE -- not the price of the averages
        per = [price_fusion(x["skill"][a]["rho"], x["skill"][b]["rho"], x["raw_corr"][k])
               for x in comp]
        fus[k] = {kk: float(np.mean([p[kk] for p in per]))
                  for kk in ("rho_fused", "gain", "unique_skill_weak", "second_order")}
        fus[k]["per_target_gain"] = [p["gain"] for p in per]
    print("\nPRICED FUSION (optimal linear multiple correlation, per target then averaged):")
    print(f"  {'pair':38s} {'rho_fused':>10s} {'gain':>8s} {'weak unique':>12s}")
    for k, v in sorted(fus.items(), key=lambda kv: -kv[1]["gain"]):
        print(f"  {k:38s} {v['rho_fused']:+10.4f} {v['gain']:+8.4f} "
              f"{v['unique_skill_weak']:+12.4f}")

    rw = reweight()
    print("\n=== 3. LEAVE-ONE-TARGET-OUT LINEAR REWEIGHTING OF THE 11 LEGACY TERMS ===")
    print(f"  in-sample (CEILING)  rho {rw['in_sample_mean_rho']:+.3f}  "
          f"decile {rw['in_sample_mean_decile']:+.3f}  sel {rw['in_sample_sel']:.3f}")
    print(f"  LOTO (honest)        rho {rw['loto_mean_rho']:+.3f}  "
          f"decile {rw['loto_mean_decile']:+.3f}  "
          f"sel {rw['loto_vs_random']['mean_a']:.3f}")
    p = rw["loto_vs_random"]
    print(f"  LOTO vs random draw  {p['mean_diff']:+.3f} "
          f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}  "
          f"drop-top10 {p['drop_top10_mean_diff']}")
    p = rw["loto_vs_legacy"]
    print(f"  LOTO vs Legacy total {p['mean_diff']:+.3f} "
          f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}")

    par = [pareto(p) for p in E.ENUM_TARGETS]
    print("\n=== 4. PARETO: is there conflict at all? ===")
    for k in par[0]["conflict"]:
        print(f"  rho({k:26s}) = {np.mean([x['conflict'][k] for x in par]):+.3f}")
    print(f"\n  mean |conflict| {np.mean([x['mean_abs_conflict'] for x in par]):.3f}   "
          f"negative pairs {np.mean([x['n_negative_pairs'] for x in par]):.1f} of "
          f"{len(par[0]['conflict'])}")
    print(f"  frontier size {np.mean([x['front_size'] for x in par]):.0f} "
          f"({np.mean([x['front_frac'] for x in par])*100:.1f}% of the pool)")
    print(f"  frontier best RMSD {np.mean([x['front_best_rmsd'] for x in par]):.3f}  "
          f"vs same-size RANDOM set {np.mean([x['null_best_rmsd'] for x in par]):.3f}  "
          f"vs pool best {np.mean([x['pool_best_rmsd'] for x in par]):.3f}")
    print(f"  fraction of the truly-best 1% retained on the frontier "
          f"{np.mean([x['frac_of_best1pct_on_front'] for x in par]):.3f}")

    E.write("ener_complement", dict(
        what="decomposition, truth-partialled complementarity, reweighting, Pareto",
        decomposition=dec, complementarity=comp, fusion_price=fus,
        reweighting=rw, pareto=par), n_expected=len(E.ENUM_TARGETS))
    return dec, comp, rw, par


if __name__ == "__main__":
    main()
