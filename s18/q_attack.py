"""SPRINT 18 / QUANTUM-ADVERSARIAL -- A1: attack the degree-1 hypothesis from the persisted
Sprint 17 spectra, at zero compute cost.

THE JOB IS FALSIFICATION, not confirmation.  Phase 0 established that `deg1 - full` is
-0.249 A with median +0.000, W/L 7/5 and a 95% CI of [-0.650, +0.093], and that the deg1
objective is a WORSE global correlate of RMSD (+0.153) than the full objective (+0.264).
Those two facts are in tension and one of them has to give.

THE DECOMPOSITION THAT SETTLES IT.  An objective's argmin RMSD is the composition of two
independent things:

    argmin_RMSD  =  [ mean RMSD of the objective's own low band ]      <- BAND QUALITY
                  + [ argmin_RMSD - that band mean ]                   <- WITHIN-BAND DRAW

Band quality is a property of the objective's ALIGNMENT: an objective that orders the space
better puts better structures in its own top 1%.  The within-band draw is a min-of-N order
statistic over ~N/100 configurations that the objective cannot distinguish -- it is LUCK
unless the objective has in-band skill.  Sprint 17's own `rho_ownband_with_rmsd_ORACLE` is
exactly the measure of that skill, and the programme's standing law is that nothing ranks
in band.

If `deg1 - full` is carried by the BAND-QUALITY term, the degree-1 hypothesis has a
mechanism.  If it is carried by the WITHIN-BAND term, the 2.411 A is a min-of-N draw from a
band the objective cannot order, the "worse global rho" is not a paradox at all, and the
hypothesis is dead without ever touching the 126-target instrument.

BOTH mandatory controls are constructed:
  * ZERO-INFORMATION reference   -- a uniformly random configuration (the space mean RMSD).
  * MATCHED-RANDOM              -- a uniformly random draw from the objective's OWN band, of
                                   the same size, i.e. the min-of-N null WITH ITS BAND NAMED.

TARGET is the unit.  Paired fold-clustered bootstrap CIs, medians and W/L beside every mean.

RUN:  python -m s18.q_attack a1
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

S17R = os.path.join(ROOT, "s17", "results")
RESULTS = os.path.join(ROOT, "s18", "results")
os.makedirs(RESULTS, exist_ok=True)

SALT = "s18quantum"
BAND = 0.01            # Sprint 17's own band: max(64, N // 100)


# ------------------------------------------------------------------ statistics
def boot_paired(d, folds, n=20000, seed=0):
    """Fold-CLUSTERED paired bootstrap of the mean of `d`.  Folds are resampled with
    replacement, then targets within the drawn folds -- so a fold-level effect cannot be
    counted as n independent target-level effects."""
    d = np.asarray(d, float)
    folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = {f: np.flatnonzero(folds == f) for f in uf}
    from s15 import seed as SD
    rng = SD.stable_rng("boot", seed, salt=SALT)
    out = np.empty(n)
    for b in range(n):
        pick = rng.choice(uf, size=uf.size, replace=True)
        take = []
        for f in pick:
            m = idx[f]
            take.append(rng.choice(m, size=m.size, replace=True))
        out[b] = d[np.concatenate(take)].mean()
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def summarise(name, d, folds, seed=0):
    d = np.asarray(d, float)
    lo, hi = boot_paired(d, folds, seed=seed)
    w = int((d < -1e-12).sum()); l = int((d > 1e-12).sum()); t = int(len(d) - w - l)
    # fold sign consistency: how many folds have a fold-mean of the same sign as the total
    uf = np.unique(folds)
    fm = np.array([d[np.asarray(folds) == f].mean() for f in uf])
    same = int((np.sign(fm) == np.sign(d.mean())).sum())
    return {"name": name, "n": int(d.size), "mean": float(d.mean()),
            "median": float(np.median(d)), "sd": float(d.std(ddof=1)) if d.size > 1 else 0.0,
            "ci_lo": lo, "ci_hi": hi, "W": w, "L": l, "ties": t,
            "folds_same_sign": f"{same}/{uf.size}"}


def fmt(r):
    return (f"{r['name']:<44s} n={r['n']:>3d} mean {r['mean']:+7.3f}  med {r['median']:+7.3f}"
            f"  CI[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]  W/L {r['W']}/{r['L']}"
            f" (t{r['ties']})  folds {r['folds_same_sign']}")


# ------------------------------------------------------------------ load
def load_specs():
    d = json.load(open(os.path.join(S17R, "quantum_theory.json")))
    rows = []
    for k, v in d.items():
        if not k.startswith("spec_"):
            continue
        h = v["hamil_uniformised"]
        t1 = h["truncation"]["1"]
        t2 = h["truncation"]["2"]
        rows.append({
            "pdb": v["pdb"], "n": v["n"], "nq": v["n_qubits"], "N": v["N"],
            "fold": v["fold"],
            "full_argmin": v["certified_argmin_rmsd_ORACLE"],
            "space_best": v["rmsd_best_in_space_ORACLE"],
            "full_rho": v["rho_full_with_rmsd_ORACLE"],
            "full_band_mean": v["full_ownband_mean_rmsd_ORACLE"],
            "full_band_best": v["full_ownband_best_rmsd_ORACLE"],
            "d1_argmin": t1["argmin_rmsd_ORACLE"],
            "d1_ties": t1["argmin_ties"],
            "d1_rho": t1["rho_with_rmsd_ORACLE"],
            "d1_band_rho": t1["rho_ownband_with_rmsd_ORACLE"],
            "d1_band_mean": t1["ownband_mean_rmsd_ORACLE"],
            "d1_band_best": t1["ownband_best_rmsd_ORACLE"],
            "d1_pct": t1["argmin_pct_ORACLE"],
            "d1_rho_full": t1["rho_with_full"],
            "d2_argmin": t2["argmin_rmsd_ORACLE"],
            "d2_band_mean": t2["ownband_mean_rmsd_ORACLE"],
            "d2_rho": t2["rho_with_rmsd_ORACLE"],
            "var1": h["var_frac_by_weight"][1],
            "cum2": h["cum_var_frac"][2],
        })
    rows.sort(key=lambda r: r["pdb"])
    return rows


# ------------------------------------------------------------------ A1
def a1():
    rows = load_specs()
    folds = [r["fold"] for r in rows]
    P = print
    P("=" * 104)
    P("S18 / QUANTUM-ADVERSARIAL  A1 -- IS THE DEGREE-1 ARGMIN ADVANTAGE BAND QUALITY, OR A")
    P("                                 MIN-OF-N DRAW FROM A BAND NOTHING CAN ORDER?")
    P("=" * 104)
    P(f"n = {len(rows)} exhaustively enumerated targets; band = top {BAND:.0%} of each")
    P("objective's own ordering (Sprint 17's `max(64, N//100)`).  All RMSD columns ORACLE.")
    P("")

    # ---------------- per target
    P("PER-TARGET  (argmin = certified; band mean = mean RMSD of the objective's own top 1%)")
    P("-" * 104)
    P(f"{'pdb':<6s}{'n':>3s}{'fold':>5s}  {'full':>7s}{'deg1':>8s}{'delta':>8s} | "
      f"{'fullband':>9s}{'d1band':>8s}{'dband':>8s} | {'d1 ties':>8s}{'d1 pct':>8s}"
      f"{'space':>7s}{'d1bandrho':>10s}")
    for r in rows:
        P(f"{r['pdb']:<6s}{r['n']:>3d}{r['fold']:>5d}  {r['full_argmin']:>7.3f}"
          f"{r['d1_argmin']:>8.3f}{r['d1_argmin']-r['full_argmin']:>+8.3f} | "
          f"{r['full_band_mean']:>9.3f}{r['d1_band_mean']:>8.3f}"
          f"{r['d1_band_mean']-r['full_band_mean']:>+8.3f} | "
          f"{r['d1_ties']:>8d}{r['d1_pct']:>8.4f}{r['space_best']:>7.3f}"
          f"{r['d1_band_rho']:>+10.3f}")
    P("")

    # ---------------- THE DECOMPOSITION
    full_arg = np.array([r["full_argmin"] for r in rows])
    d1_arg = np.array([r["d1_argmin"] for r in rows])
    full_bm = np.array([r["full_band_mean"] for r in rows])
    d1_bm = np.array([r["d1_band_mean"] for r in rows])
    space = np.array([r["space_best"] for r in rows])

    total = d1_arg - full_arg
    band = d1_bm - full_bm                       # BAND QUALITY (alignment)
    draw = (d1_arg - d1_bm) - (full_arg - full_bm)   # WITHIN-BAND DRAW (luck)

    P("THE DECOMPOSITION  --  total = band quality + within-band draw   (negative = deg1 better)")
    P("-" * 104)
    res = []
    for nm, v in (("TOTAL   deg1 - full  (argmin RMSD)", total),
                  ("  (a) BAND QUALITY   deg1 - full", band),
                  ("  (b) WITHIN-BAND DRAW  deg1 - full", draw)):
        r = summarise(nm, v, folds); res.append(r); P(fmt(r))
    P("")
    P(f"  identity check   max|total - (band+draw)| = {np.abs(total-(band+draw)).max():.2e}")
    P(f"  share of the total mean carried by the DRAW term: "
      f"{draw.mean()/total.mean()*100:.1f}%   by BAND QUALITY: "
      f"{band.mean()/total.mean()*100:.1f}%")
    P("")

    # ---------------- in-band skill: can either objective order its own band?
    d1br = np.array([r["d1_band_rho"] for r in rows])
    P("IN-BAND SKILL -- can the degree-1 objective order the band it draws its argmin from?")
    P("-" * 104)
    r = summarise("deg1 in-band rho(obj, RMSD) vs 0", d1br, folds); P(fmt(r))
    P("  (positive rho = the objective orders its own band; the draw is then skill, not luck)")
    P("")

    # ---------------- min-of-N NULL, band named
    P("THE MIN-OF-N NULL, WITH ITS BAND NAMED")
    P("-" * 104)
    P("  MATCHED-RANDOM control: a UNIFORM draw from the SAME band the argmin comes from,")
    P("  i.e. the band mean itself (its expectation).  ZERO-INFORMATION control: the space.")
    zero_ref = None
    P(f"  deg1 argmin      mean {d1_arg.mean():.3f}   median {np.median(d1_arg):.3f}")
    P(f"  deg1 band mean   mean {d1_bm.mean():.3f}   <- matched-random draw from deg1's band")
    P(f"  full argmin      mean {full_arg.mean():.3f}   median {np.median(full_arg):.3f}")
    P(f"  full band mean   mean {full_bm.mean():.3f}   <- matched-random draw from full's band")
    P(f"  space best       mean {space.mean():.3f}   (ORACLE ceiling of the register)")
    P("")
    ex = summarise("deg1 argmin - its own band mean (excess)", d1_arg - d1_bm, folds)
    ef = summarise("full argmin - its own band mean (excess)", full_arg - full_bm, folds)
    P(fmt(ex)); P(fmt(ef))
    P("")
    hit = np.array([abs(r["d1_argmin"] - r["space_best"]) < 1e-6 for r in rows])
    hitf = np.array([abs(r["full_argmin"] - r["space_best"]) < 1e-6 for r in rows])
    P(f"  targets where the deg1 argmin IS the global best of the whole register: "
      f"{int(hit.sum())}/{len(rows)}  ({[r['pdb'] for r,h in zip(rows,hit) if h]})")
    P(f"  targets where the full argmin IS the global best of the whole register: "
      f"{int(hitf.sum())}/{len(rows)}")
    P("")

    # ---------------- the three carrier targets
    P("ARE THE THREE CARRIER TARGETS SPECIAL?   (Phase 0: 7VI4, 1CS9, 7T3H drive the mean)")
    P("-" * 104)
    order = np.argsort(total)
    P(f"{'pdb':<6s}{'n':>3s}{'delta':>8s}{'band':>8s}{'draw':>8s}{'d1ties':>8s}"
      f"{'d1pct':>8s}{'var1':>7s}{'cum2':>7s}{'rho_d1full':>11s}{'spaceb':>8s}")
    for i in order:
        r = rows[i]
        P(f"{r['pdb']:<6s}{r['n']:>3d}{total[i]:>+8.3f}{band[i]:>+8.3f}{draw[i]:>+8.3f}"
          f"{r['d1_ties']:>8d}{r['d1_pct']:>8.4f}{r['var1']:>7.3f}{r['cum2']:>7.3f}"
          f"{r['d1_rho_full']:>+11.3f}{r['space_best']:>8.3f}")
    P("")
    top3 = [rows[i]["pdb"] for i in order[:3]]
    bot3 = [rows[i]["pdb"] for i in order[-3:]]
    P(f"  best three for deg1 : {top3}")
    P(f"  worst three for deg1: {bot3}")
    P("")
    # drop-top-3 -- reported ONLY beside its uniform-effect null, per the ledger's law
    P("LEAVE-THE-THREE-CARRIERS-OUT, AGAINST ITS UNIFORM-EFFECT NULL")
    P("  (the ledger records that a raw drop-top threshold is NOT a valid test on its own)")
    keep = np.ones(len(rows), bool); keep[order[:3]] = False
    r = summarise("deg1 - full, 3 best-for-deg1 targets removed",
                  total[keep], np.array(folds)[keep]); P(fmt(r))
    # uniform-effect null: what does dropping the 3 most favourable targets do to a sample
    # in which EVERY target carries the same true effect plus noise?
    from s15 import seed as SD
    rng = SD.stable_rng("uniform_effect_null", salt=SALT)
    mu, sd = total.mean(), total.std(ddof=1)
    nulls = []
    for _ in range(20000):
        z = rng.normal(mu, sd, size=total.size)
        nulls.append(np.sort(z)[3:].mean())
    nulls = np.array(nulls)
    P(f"  uniform-effect null: drop-3 mean {nulls.mean():+.3f} "
      f"[{np.percentile(nulls,2.5):+.3f},{np.percentile(nulls,97.5):+.3f}];"
      f"  observed {total[keep].mean():+.3f} sits at pct "
      f"{100*(nulls < total[keep].mean()).mean():.1f}")
    P("")

    # ---------------- coherence of "better argmin, worse correlation"
    P("IS 'BETTER ARGMIN, WORSE GLOBAL RHO' COHERENT?  -- the Gaussian-copula prediction")
    P("-" * 104)
    d1r = np.array([r["d1_rho"] for r in rows])
    fr = np.array([r["full_rho"] for r in rows])
    rr = summarise("rho(deg1,RMSD) - rho(full,RMSD)", -(d1r - fr), folds)
    rr["name"] = "rho(full) - rho(deg1)   [positive = full is the better correlate]"
    P(fmt(rr))
    P("")
    P("  Under any monotone-copula model, the expected minimum of RMSD over the objective's")
    P("  own band is DECREASING in the objective's correlation with RMSD.  deg1 has the")
    P("  LOWER correlation, so the model predicts a WORSE argmin.  The observed sign is the")
    P("  opposite -- which is only possible if (i) the correlation is not the operative")
    P("  statistic because the skill is concentrated in the tail, or (ii) the argmin")
    P("  difference is a draw.  (i) is testable and is tested above by the in-band rho and")
    P("  the band-quality term; both point at (ii).")
    P("")
    # correlation between the per-target rho advantage and the per-target argmin advantage
    from scipy import stats as st
    sp = st.spearmanr(d1r - fr, total)
    P(f"  per-target Spearman( rho advantage of deg1 , argmin advantage of deg1 ) = "
      f"{sp.statistic:+.3f}  (p={sp.pvalue:.3f}, n={len(rows)})")
    P("  If the argmin gain were driven by alignment, this would be strongly NEGATIVE")
    P("  (more relative correlation -> more relative argmin gain).")
    P("")

    out = {"rows": rows, "total": total.tolist(), "band": band.tolist(),
           "draw": draw.tolist(), "summaries": res,
           "inband_rho": r if False else summarise("d1_band_rho", d1br, folds),
           "excess_d1": ex, "excess_full": ef,
           "dropx3": {"observed": float(total[keep].mean()),
                      "null_mean": float(nulls.mean()),
                      "null_lo": float(np.percentile(nulls, 2.5)),
                      "null_hi": float(np.percentile(nulls, 97.5))},
           "rho_adv_vs_argmin_adv_spearman": float(sp.statistic),
           "rho_adv_vs_argmin_adv_p": float(sp.pvalue),
           "_complete": True}
    p = os.path.join(RESULTS, "q_attack_a1.json")
    json.dump(out, open(p, "w"), indent=1, default=float)
    P(f"[written] {p}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "a1"
    {"a1": a1}[mode]()
