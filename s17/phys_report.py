"""s17/phys_report.py -- all analysis of the identical-candidate instrument.  No physics run
here; this module only reads `s17/results/phys_ident*.json` and reports.

E2  Legacy decomposed -- which components carry information, and what does Legacy contain
    that is NOT already in the distance model (nested ablation, not raw correlation).
E3  Legacy as a GATE, and the danger -- near-native recall through the gate, against matched
    random rejection of the same count.
E4  the decisive Legacy-vs-AMBER experiment, on identical candidates, in four roles.
E5  disagreement as information.

STANDING CONVENTIONS ENFORCED HERE.
* IN-BAND is the programme's only meaningful ranking metric: `d <= pool_best + 1.5 A`
  (`s12.instrument.BAND`).  Global correlations are reported beside it and are known to be
  inflated by garbage rejection.
* Every selection goes through `energy_lib.argmin_tied`, which averages over the tied argmin
  set.  `np.argmin` on a tied score reads the cache's ORACLE sort order.
* Every operation is priced against a MATCHED RANDOM operation of the same count, and against
  a ZERO-INFORMATION reference.
* AMBER statistics are RANK-based only.  Raw single points on unrelaxed structures put 99.9%
  of their variance in ten configurations.
* Legacy weights are DEFAULT_WEIGHTS and are never fitted.
* `d` is ORACLE and is used only as a label and a ceiling.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I                     # noqa: E402
from s15.seed import stable_rng                     # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s17 import phys_lib as P                       # noqa: E402

BAND = I.BAND           # 1.5 A above the pool best -- the standing in-band definition
NEAR = 2.0
GATE_F = (0.10, 0.25, 0.50)
N_RAND = 3

LEG = ["s_leg_" + t for t in EL.LEG_TERMS]
AMB = ["s_amb_" + t for t in ("bond", "angle", "torsion", "nonbonded", "solvation")]


def load(name="phys_ident_leg.json"):
    return json.load(open(os.path.join(RESULTS, name)))["rows"]


#: ------------------------------------------------------------------ zero-information
#: THE FLOOR CONTROLS.  A physics ranker's bar is NOT "beat the distance objective" -- the
#: distance objective has essentially no in-band skill (SELECT: -0.014 [-0.118, +0.094] at
#: n = 126 against random inside its own top-25).  The bar is a ZERO-INFORMATION constant
#: secondary-structure reference, and the incumbent native-free signal, pool typicality.
#: All four are added here and every physics score is quoted beside them.
HELIX = (-57.0, -47.0)
STRAND = (-139.0, 135.0)


def augment(rows, want_typicality=True, verbose=True):
    """Add the zero-information and incumbent native-free score columns to each row.

    s_helix  CA-RMSD of the candidate to a CONSTANT ideal alpha-helix of the same length
    s_strand CA-RMSD of the candidate to a CONSTANT ideal beta-strand of the same length
             -- both know nothing about the target beyond its length
    s_typic  mean pairwise CA-RMSD to the rest of the candidate set (pool typicality /
             consensus), the only native-free signal this programme has found with positive
             in-band skill
    s_rand   a uniform random score, the matched-random rank feature
    """
    for c, r in enumerate(rows):
        pdb, n = r["pdb"], int(r["n"])
        u = I.load_univ(pdb)
        order = np.asarray(u["order"], int)[:len(r["d"])]
        W = np.asarray(u["W"], float)[order]
        #: UNITS: `phys_lib.constant_backbone` converts to radians, which is what
        #: `core.geometry.build_backbone` actually takes.  Passing degrees silently builds
        #: a different constant conformation (see the note in phys_lib).
        for tag, (ph, ps) in (("helix", HELIX), ("strand", STRAND)):
            ref = P.constant_backbone(n, ph, ps)["CA"]
            r["s_" + tag] = I.kabsch_rmsd_batch(W, np.asarray(ref, float)).tolist()
        if want_typicality:
            Pm = I.pairwise_rmsd(W)
            r["s_typic"] = Pm.mean(1).tolist()
        rng = stable_rng("s17phys", "randscore", pdb)
        r["s_rand"] = rng.random(len(r["d"])).tolist()
        if verbose and (c + 1) % 40 == 0:
            print(f"    augment {c+1}/{len(rows)}", flush=True)
    return rows


def _rank(x):
    return EL._avg_rank(np.asarray(x, float))


def _partial_rho(a, b, c):
    """Spearman(a, b) with c partialled out, all rank-transformed first."""
    ra, rb, rc = _rank(a), _rank(b), _rank(c)
    rc = rc - rc.mean()
    n2 = float((rc ** 2).sum())
    if n2 <= 0:
        return float("nan")
    ea = (ra - ra.mean()) - ((ra - ra.mean()) * rc).sum() / n2 * rc
    eb = (rb - rb.mean()) - ((rb - rb.mean()) * rc).sum() / n2 * rc
    den = np.sqrt((ea ** 2).sum() * (eb ** 2).sum())
    return float((ea * eb).sum() / den) if den > 0 else float("nan")


def _agg(vals, folds, names=None, seed=0, n_boot=4000):
    """Mean, median, fold-clustered and i.i.d. bootstrap CI over targets."""
    v = np.asarray(vals, float)
    ok = np.isfinite(v)
    v = v[ok]; f = np.asarray(folds)[ok]
    if len(v) < 3:
        return {"n": int(len(v)), "mean": float("nan")}
    rng = np.random.default_rng(seed)
    bs = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(n_boot)])
    uf = np.unique(f); idx = {u: np.flatnonzero(f == u) for u in uf}
    cb = np.array([v[np.concatenate([idx[uf[p]] for p in
                                     rng.integers(0, len(uf), len(uf))])].mean()
                   for _ in range(n_boot)])
    return {"n": int(len(v)), "mean": float(v.mean()), "median": float(np.median(v)),
            "sd": float(v.std(ddof=1)),
            "ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "ci_fold": [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))],
            "pos": int((v > 0).sum()), "neg": int((v < 0).sum())}


def _ci_txt(a):
    if not np.isfinite(a.get("mean", np.nan)):
        return "     --"
    return (f"{a['mean']:+.3f} [{a['ci'][0]:+.3f},{a['ci'][1]:+.3f}] "
            f"f[{a['ci_fold'][0]:+.3f},{a['ci_fold'][1]:+.3f}]")


# =====================================================================  E2
def e2(rows, scores=None, verbose=True):
    """Legacy (and AMBER) decomposed: detection, ranking, and INCREMENTAL value over the
    distance model by nested ablation."""
    folds = np.array([r["fold"] for r in rows])
    zero = [c for c in ("s_helix", "s_strand", "s_typic", "s_rand") if c in rows[0]]
    cols = scores if scores is not None else (
        ["s_dist"] + zero + ["s_legacy"] + LEG
        + (["s_amber"] + AMB if "s_amber" in rows[0] else []))
    out = {}
    for col in cols:
        per = {k: [] for k in ("rho_glob", "rho_band", "auroc_band", "auroc_dec",
                               "sel", "sel_band", "pct_best", "fnr_near", "fpr",
                               "part_rho_band", "part_rho_glob")}
        for r in rows:
            if col not in r:
                for k in per:
                    per[k].append(np.nan)
                continue
            s = np.asarray(r[col], float)
            d = np.asarray(r["d"], float)
            dist = np.asarray(r["s_dist"], float)
            ok = np.isfinite(s)
            if ok.sum() < 20 or np.nanstd(s[ok]) == 0:
                for k in per:
                    per[k].append(np.nan)
                continue
            s, d, dist = s[ok], d[ok], dist[ok]
            band = d <= d.min() + BAND
            dec = d <= np.percentile(d, 10)
            per["rho_glob"].append(EL.spearman(s, d))
            per["rho_band"].append(EL.spearman(s[band], d[band]) if band.sum() >= 5
                                   else np.nan)
            per["auroc_band"].append(EL.auroc(-s, band))
            per["auroc_dec"].append(EL.auroc(-s, dec))
            per["sel"].append(EL.argmin_tied(s, d)[0])
            per["sel_band"].append(EL.argmin_tied(s[band], d[band])[0]
                                   if band.sum() >= 5 else np.nan)
            per["pct_best"].append(float((s < s[np.argmin(d)]).mean()))
            # FNR on near-native: near-native candidates landing in the worst quartile of s
            near = d < NEAR
            q = np.percentile(s, 75)
            per["fnr_near"].append(float((s[near] > q).mean()) if near.sum() else np.nan)
            per["fpr"].append(float((s[d > np.median(d)] <= q).mean()))
            per["part_rho_glob"].append(_partial_rho(s, d, dist))
            per["part_rho_band"].append(_partial_rho(s[band], d[band], dist[band])
                                        if band.sum() >= 8 else np.nan)
        out[col] = {k: _agg(v, folds) for k, v in per.items()}
        out[col]["_per_target"] = {k: list(map(float, v)) for k, v in per.items()}
    if verbose:
        print("\nE2.  LEGACY (and AMBER) DECOMPOSED, on the identical K=500 candidate set")
        print("     rho: Spearman(score, CA-RMSD), positive = the score orders correctly.")
        print("     IN-BAND is the only meaningful column; GLOBAL is inflated by garbage "
              "rejection.")
        print(f"     {'score':>18}{'rho glob':>10}{'rho BAND':>10}{'AUROC band':>12}"
              f"{'sel A':>8}{'pct of best':>12}{'FNR<2A':>9}")
        for col in cols:
            o = out[col]
            print(f"     {col:>18}{o['rho_glob']['mean']:>10.3f}"
                  f"{o['rho_band']['mean']:>10.3f}{o['auroc_band']['mean']:>12.3f}"
                  f"{o['sel']['mean']:>8.3f}{o['pct_best']['mean']:>12.3f}"
                  f"{o['fnr_near']['mean']:>9.3f}")
        print("\n     INCREMENTAL VALUE OVER THE DISTANCE MODEL -- partial Spearman with the")
        print("     shipped distance score partialled out, IN-BAND.  This is the nested")
        print("     ablation the sprint asks for; a raw correlation is not it.")
        print(f"     {'score':>18}   partial rho in-band [i.i.d. CI] [fold-clustered CI]")
        for col in cols:
            if col == "s_dist":
                continue
            print(f"     {col:>18}   {_ci_txt(out[col]['part_rho_band'])}")

        print("\n     THE OPERATIVE BAR -- IN-BAND SELECTION against MATCHED RANDOM "
              "selection from the same band.")
        print("     `s_rand` IS that matched random operation, so its row is the null by "
              "construction.")
        print("     `s_helix` / `s_strand` are ZERO-INFORMATION constant secondary-structure "
              "references.")
        print(f"     {'score':>18}{'sel in-band':>13}   vs matched-random in-band "
              f"[i.i.d.] [fold]")
        if "s_rand" in out:
            base = np.asarray(out["s_rand"]["_per_target"]["sel_band"], float)
            for col in cols:
                v = np.asarray(out[col]["_per_target"]["sel_band"], float)
                m = np.isfinite(v) & np.isfinite(base)
                if m.sum() < 5:
                    continue
                pr = P.paired(v[m], base[m], folds=folds[m])
                out[col]["sel_band_vs_rand"] = pr
                print(f"     {col:>18}{np.nanmean(v):>13.3f}   {pr['mean']:+.4f} "
                      f"[{pr['ci'][0]:+.4f},{pr['ci'][1]:+.4f}] "
                      f"f[{pr['ci_fold'][0]:+.4f},{pr['ci_fold'][1]:+.4f}] "
                      f"{pr['W']:>3}W/{pr['L']:<3}L")
    return out


# =====================================================================  E3
def e3(rows, gate_col="s_legacy", verbose=True):
    """Legacy as a GATE, and the danger: does it remove exactly the candidates that matter?"""
    folds = np.array([r["fold"] for r in rows])
    res = {}
    for f in GATE_F:
        keep = {k: [] for k in ("sel", "best", "near_surv", "best_surv", "band_surv")}
        rnd = {k: [] for k in keep}
        base = {k: [] for k in ("sel", "best", "near_surv", "best_surv", "band_surv")}
        for r in rows:
            s = np.asarray(r[gate_col], float)
            d = np.asarray(r["d"], float)
            dist = np.asarray(r["s_dist"], float)
            k = len(d)
            nkeep = int(round(k * (1.0 - f)))
            band = d <= d.min() + BAND
            near = d < NEAR
            ib = int(np.argmin(d))

            base["sel"].append(EL.argmin_tied(dist, d)[0])
            base["best"].append(float(d.min()))
            base["near_surv"].append(1.0)
            base["best_surv"].append(1.0)
            base["band_surv"].append(1.0)

            srt = np.argsort(s, kind="stable")[:nkeep]
            keep["sel"].append(EL.argmin_tied(dist[srt], d[srt])[0])
            keep["best"].append(float(d[srt].min()))
            keep["near_surv"].append(float(near[srt].sum() / max(1, near.sum()))
                                     if near.sum() else np.nan)
            keep["best_surv"].append(float(ib in set(srt.tolist())))
            keep["band_surv"].append(float(band[srt].sum() / max(1, band.sum())))

            rs, rb, rn, rbs, rband = [], [], [], [], []
            for t in range(N_RAND):
                rng = stable_rng("s17phys", "gate", gate_col, f, r["pdb"], t)
                pick = rng.permutation(k)[:nkeep]
                rs.append(EL.argmin_tied(dist[pick], d[pick])[0])
                rb.append(float(d[pick].min()))
                rn.append(float(near[pick].sum() / max(1, near.sum()))
                          if near.sum() else np.nan)
                rbs.append(float(ib in set(pick.tolist())))
                rband.append(float(band[pick].sum() / max(1, band.sum())))
            rnd["sel"].append(np.mean(rs)); rnd["best"].append(np.mean(rb))
            rnd["near_surv"].append(np.nanmean(rn) if np.isfinite(rn).any() else np.nan)
            rnd["best_surv"].append(np.mean(rbs)); rnd["band_surv"].append(np.mean(rband))

        res[f] = {
            "keep": {k: _agg(v, folds) for k, v in keep.items()},
            "rand": {k: _agg(v, folds) for k, v in rnd.items()},
            "base": {k: _agg(v, folds) for k, v in base.items()},
            "vs_rand_sel": P.paired(keep["sel"], rnd["sel"], folds=folds),
            "vs_base_sel": P.paired(keep["sel"], base["sel"], folds=folds),
            "vs_rand_near": P.paired(keep["near_surv"], rnd["near_surv"], folds=folds),
            "vs_rand_best": P.paired(keep["best"], rnd["best"], folds=folds),
        }
    if verbose:
        print(f"\nE3.  {gate_col.upper()} AS A GATE -- reject the worst f, then select by the "
              "distance score.")
        print("     The primary danger, tested directly: does the gate remove the candidates "
              "that matter?")
        print(f"     {'f':>6}{'sel A':>9}{'rand sel':>10}{'gate-rand':>26}"
              f"{'set best':>10}{'rand best':>11}{'near<2A kept':>14}{'rand kept':>11}"
              f"{'pool-best kept':>16}")
        for f in GATE_F:
            R = res[f]
            v = R["vs_rand_sel"]
            print(f"     {f:>6.2f}{R['keep']['sel']['mean']:>9.3f}"
                  f"{R['rand']['sel']['mean']:>10.3f}"
                  f"   {v['mean']:+.4f} [{v['ci'][0]:+.4f},{v['ci'][1]:+.4f}] {v['W']:>3}/{v['L']:<3}"
                  f"{R['keep']['best']['mean']:>10.3f}{R['rand']['best']['mean']:>11.3f}"
                  f"{R['keep']['near_surv']['mean']:>14.3f}"
                  f"{R['rand']['near_surv']['mean']:>11.3f}"
                  f"{R['keep']['best_surv']['mean']:>16.3f}")
        print(f"     no gate: sel {res[GATE_F[0]]['base']['sel']['mean']:.3f}   "
              f"pool best {res[GATE_F[0]]['base']['best']['mean']:.3f}")
        print("     THE DANGER, paired against the matched random gate of the SAME count:")
        print(f"     {'f':>6}{'d near<2A recall':>34}{'d pool-best survival':>34}"
              f"{'d set BEST (A)':>32}")
        for f in GATE_F:
            R = res[f]
            a = R["vs_rand_near"]; b = R["vs_rand_best"]
            bs = P.paired(R["keep"]["best_surv"]["mean"] * np.ones(1),
                          R["rand"]["best_surv"]["mean"] * np.ones(1))
            print(f"     {f:>6.2f}"
                  f"   {a['mean']:+.4f} [{a['ci'][0]:+.4f},{a['ci'][1]:+.4f}] "
                  f"{a['W']:>3}/{a['L']:<3}"
                  f"   {R['keep']['best_surv']['mean']:.3f} vs rand "
                  f"{R['rand']['best_surv']['mean']:.3f}      "
                  f"{b['mean']:+.4f} [{b['ci'][0]:+.4f},{b['ci'][1]:+.4f}] "
                  f"{b['W']:>3}/{b['L']:<3}")
    return res


# =====================================================================  E5
def e5(rows, verbose=True):
    """Disagreement as information: do the rankers' disagreements predict when the distance
    selector fails?"""
    folds = np.array([r["fold"] for r in rows])
    feats, resp = {}, {"err": [], "regret": []}
    have_amber = "s_amber" in rows[0]
    names = []
    for r in rows:
        d = np.asarray(r["d"], float)
        dist = np.asarray(r["s_dist"], float)
        leg = np.asarray(r["s_legacy"], float)
        band = d <= d.min() + BAND
        pick_d = int(np.argmin(dist)); pick_l = int(np.argmin(leg))
        f = {}
        f["rho_dl"] = EL.spearman(dist, leg)
        f["rho_dl_band"] = EL.spearman(dist[band], leg[band]) if band.sum() >= 5 else np.nan
        top = 25
        sd_ = set(np.argsort(dist, kind="stable")[:top].tolist())
        sl_ = set(np.argsort(leg, kind="stable")[:top].tolist())
        f["overlap_dl"] = len(sd_ & sl_) / top
        f["rank_of_distpick_in_leg"] = float((leg < leg[pick_d]).mean())
        f["rank_of_legpick_in_dist"] = float((dist < dist[pick_l]).mean())
        if have_amber:
            am = np.asarray(r["s_amber"], float)
            m = np.isfinite(am)
            f["rho_da"] = EL.spearman(dist[m], am[m])
            f["rho_la"] = EL.spearman(leg[m], am[m])
            sa = set(np.argsort(np.where(m, am, np.inf), kind="stable")[:top].tolist())
            f["overlap_da"] = len(sd_ & sa) / top
            f["overlap_la"] = len(sl_ & sa) / top
        #: NATIVE-FREE incumbent predictor: pool typicality of the distance pick.
        f["score_margin"] = float((np.sort(dist)[1] - dist[pick_d]) /
                                  (np.std(dist) + 1e-12))
        f["dist_score"] = float(dist[pick_d])
        for k, v in f.items():
            feats.setdefault(k, []).append(v)
        resp["err"].append(EL.argmin_tied(dist, d)[0])
        resp["regret"].append(EL.argmin_tied(dist, d)[0] - float(d.min()))
        names.append(r["pdb"])

    #: leave-fold-out Spearman of each single feature with the response
    out = {"features": {}, "response": {k: list(map(float, v)) for k, v in resp.items()}}
    for k, v in feats.items():
        v = np.asarray(v, float)
        row = {}
        for rk in ("err", "regret"):
            y = np.asarray(resp[rk], float)
            row[rk + "_rho"] = EL.spearman(v, y)
            #: leave-fold-out: sign fixed on 4 folds, skill scored on the 5th
            lo = []
            for fo in np.unique(folds):
                tr = folds != fo; te = folds == fo
                if te.sum() < 4 or np.nanstd(v[tr]) == 0:
                    continue
                sgn = np.sign(EL.spearman(v[tr], y[tr]))
                if sgn == 0:
                    sgn = 1.0
                lo.append(sgn * EL.spearman(v[te], y[te]))
            row[rk + "_lfo"] = float(np.mean(lo)) if lo else float("nan")
            row[rk + "_lfo_agg"] = _agg(lo, np.arange(len(lo))) if len(lo) >= 3 else None
        row["values"] = list(map(float, v))
        out["features"][k] = row
    #: INCREMENTAL over the trivial predictor.  `dist_score` -- how good the objective thinks
    #: its own pick is -- is a native-free target-level feature that needs no second ranker.
    #: A disagreement feature is only a result if it adds to THAT, so every feature is
    #: re-scored with `dist_score` partialled out, leave-fold-out.
    trivial = np.asarray(feats["dist_score"], float)
    for k, v in feats.items():
        if k == "dist_score":
            continue
        v = np.asarray(v, float)
        for rk in ("err", "regret"):
            y = np.asarray(resp[rk], float)
            lo = []
            for fo in np.unique(folds):
                tr = folds != fo; te = folds == fo
                if te.sum() < 4 or np.nanstd(v[tr]) == 0:
                    continue
                sgn = np.sign(_partial_rho(v[tr], y[tr], trivial[tr])) or 1.0
                lo.append(sgn * _partial_rho(v[te], y[te], trivial[te]))
            out["features"][k][rk + "_lfo_partial"] = (float(np.nanmean(lo)) if lo
                                                       else float("nan"))
    out["features"]["dist_score"]["err_lfo_partial"] = float("nan")
    out["features"]["dist_score"]["regret_lfo_partial"] = float("nan")

    if verbose:
        print("\nE5.  DISAGREEMENT AS INFORMATION -- does ranker disagreement predict when "
              "the distance selector fails?")
        print("     LFO = leave-fold-out: the sign of the feature-response relation is fixed "
              "on 4 folds and")
        print("     the Spearman is scored on the held-out fold.  A raw in-sample rho is not "
              "a result.")
        print("     The last two columns partial out `dist_score` -- how good the objective")
        print("     thinks its own pick is.  That is a native-free target-level feature that")
        print("     needs NO second ranker, so it is the baseline a disagreement feature must")
        print("     beat.  A feature is a result only if its PARTIAL column survives.")
        print(f"     {'feature':>26}{'rho(err)':>10}{'LFO(err)':>10}"
              f"{'rho(regret)':>13}{'LFOpart(err)':>14}{'LFOpart(reg)':>14}")
        for k, row in out["features"].items():
            print(f"     {k:>26}{row['err_rho']:>10.3f}{row['err_lfo']:>10.3f}"
                  f"{row['regret_rho']:>13.3f}"
                  f"{row.get('err_lfo_partial', float('nan')):>14.3f}"
                  f"{row.get('regret_lfo_partial', float('nan')):>14.3f}")
    return out


def main(name="phys_ident_leg.json"):
    rows = augment(load(name))
    print(f"IDENTICAL-CANDIDATE INSTRUMENT: n = {len(rows)} targets, K = "
          f"{len(rows[0]['d'])} candidates each.")
    try:
        om = json.load(open(os.path.join(RESULTS, "oracle_map.json")))["rows"]
        by = {r["pdb"]: r for r in om}
        b5 = np.array([by[r["pdb"]]["cells"]["500"]["best"] for r in rows if r["pdb"] in by])
        bf = np.array([by[r["pdb"]]["cells"]["full"]["best"] for r in rows if r["pdb"] in by])
        print(f"SHORTLIST COST, stated before any physics number: the K = 500 candidate set's")
        print(f"  ORACLE best is {b5.mean():.3f} A against {bf.mean():.3f} A over the full "
              f"universe -- the truncation costs {b5.mean()-bf.mean():+.3f} A of ceiling "
              f"on n = {len(b5)}.")
    except Exception as ex:
        print(f"  (oracle_map not readable for the shortlist cost: {ex})")
    o = {"e2": e2(rows), "e3": e3(rows, "s_legacy")}
    #: the gate's best chance: its single best component rather than the eleven-term total
    o["e3_torsion"] = e3(rows, "s_leg_torsion")
    o["e3_steric"] = e3(rows, "s_leg_steric")
    o["e5"] = e5(rows)
    json.dump(o, open(os.path.join(RESULTS, "phys_report.json"), "w"),
              default=lambda x: float(x) if isinstance(x, (np.floating,)) else str(x))
    return o


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phys_ident_leg.json")
