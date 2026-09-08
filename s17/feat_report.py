"""s17/feat_report.py -- the FEATURES lane's verdict table.

READING ORDER, and it is deliberate:

  A. the band definitions and what a zero-information signal scores in each.  Three bands are
     reported because the bar MOVES between them: SELECT measured the constant alpha-helix at
     +0.053 in the shipped top-25 while PHYSICS measured it at +0.291 on pool_best + 1.5 A.
     Those are not contradictory, they are different objects, and a falsifier stated against
     "beat a constant alpha-helix" is meaningless until the band is named.
  B. in-band discrimination of every representation-derived feature, on ALL THREE AXES
     (Spearman measured; pairwise ordering accuracy COUNTED, not converted; copula rho derived
     through `selection_theory`'s named functions).
  C. THE BAR.  selected CA-RMSD against random-in-band, against the constant alpha-helix and
     against the constant beta-strand, paired at the target level with fold-clustered intervals.
  D. THE DECIDING CONTROLS.  esm vs onehot at matched architecture; PLL vs its background and
     vs BLOSUM.  A representation result that does not separate from its one-hot control is not
     a representation result.
  E. the min-of-N null beside any per-target best-of-V statistic (LEDGER L8).
  F. if anything clears C and D, the CEILING it implies through the band curve -- not a model.
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

from s17 import sel_lib as L               # noqa: E402
from s17 import feat_lib as F              # noqa: E402
from s17 import selection_theory as ST     # noqa: E402

BANDS = ("top25", "top75", "oracle1.5")
LABEL = {"top25": "NATIVE-FREE, the shipped objective's own top-25 (a deployable shortlist)",
         "top75": "NATIVE-FREE, the shipped objective's own top-75 (a deployable shortlist)",
         "oracle1.5": "ORACLE-DEFINED, pool best + 1.5 A -- DIAGNOSTIC, needs the native to draw"}

#: the representation-derived arms, in the order the pre-registration lists them
FAMILY = [
    ("F6  window embeddings", ["emb_cos", "emb_l2", "emb_pool"]),
    ("F1  ESM contact head", ["esmcon_bce", "esmcon_agree", "esmcon_corr", "esmcon_topd"]),
    ("F2  exogenous pair weights", ["risk_wcon", "risk_wconf", "risk_wshipcon", "risk_wone"]),
    ("F3  window-sequence PLL", ["pll_mask", "pll_nat", "pll_bg", "blosum_sim"]),
    ("F4  embedding -> local env", ["env_esm", "env_onehot", "env_const"]),
    ("F5  predicted vs realised SS", ["ss_esm", "ss_onehot", "ss_marg", "ss_const"]),
    ("ESM-FREE TWINS of F1 (compactness controls)", ["cf_uniform", "cf_topd_uniform", "rg_raw"]),
    ("REFERENCE (not a bar)", ["dist_shipped"]),
    ("ZERO-INFORMATION BAR", ["NULL_alpha", "NULL_beta", "NULL_chance"]),
]


def load(path=None):
    p = path or os.path.join(F.RESULTS, "feat_inband.json")
    return json.load(open(p))["rows"]


def _get(rows, b, nm, key):
    return np.array([r["band"][b]["sig"][nm][key] for r in rows], float)


def report(rows=None, out=None):
    rows = rows if rows is not None else load()
    buf = []

    def P(s=""):
        buf.append(s)
        print(s, flush=True)

    #: filter PER BAND, not globally.  A few targets have an ORACLE-defined band with fewer than
    #: 8 members; dropping them from the NATIVE-FREE top-25/top-75 arms too would silently change
    #: the target set between tables.
    allrows = rows
    n = len(allrows)
    have = sorted(allrows[0]["band"]["top25"]["sig"].keys())

    P("=" * 112)
    P(f"FEATURES WORKSTREAM -- DO LEARNED SEQUENCE REPRESENTATIONS DISCRIMINATE IN-BAND?"
      f"   n = {n}, K = 500")
    P("=" * 112)
    P("Pre-registration: s17/PREREG_features.md.  Sign of every feature fixed a priori,")
    P("LOWER = PREDICTED BETTER, so this battery is not a hidden best-of-2V screen.")
    P("")
    P("A. THE BANDS, and what nothing-at-all scores in each")
    P(f"   {'band':<12}{'size':>7}{'ORACLE best':>13}{'random-in-band':>16}"
      f"{'alpha-helix':>13}{'beta-strand':>13}{'dist argmin':>13}")
    for b in BANDS:
        rs = [r for r in allrows if r["band"][b].get("sig")]
        sz = np.array([r["band"][b]["size"] for r in rs], float)
        bb = np.array([r["band"][b]["best"] for r in rs], float)
        bm = np.array([r["band"][b]["mean"] for r in rs], float)
        P(f"   {b:<12}{sz.mean():>7.1f}{bb.mean():>13.3f}{bm.mean():>16.3f}"
          f"{_get(rs, b, 'NULL_alpha', 'sel').mean():>13.3f}"
          f"{_get(rs, b, 'NULL_beta', 'sel').mean():>13.3f}"
          f"{_get(rs, b, 'dist_shipped', 'sel').mean():>13.3f}   n = {len(rs)}")
    P("   random-in-band is the EXACT expectation of picking one member uniformly (= band mean),")
    P("   not a few seeded draws, so no sampling noise enters the paired difference.")

    #: GLOBAL vs IN-BAND is the axis the whole sprint turns on (L4): the shipped objective has
    #: global Spearman 0.568 and in-band 0.131, i.e. its skill is garbage-versus-plausible.  A new
    #: feature has to be placed on BOTH axes or there is no way to tell a second garbage filter
    #: from a genuine in-band discriminator.
    gk = [k for k in have if k in rows[0].get("global", {})]
    if gk:
        P("")
        P("A2. GLOBAL vs IN-BAND, the axis the sprint turns on (L4)")
        P(f"   {'feature':<18}{'GLOBAL rhoS':>13}{'in-band top25':>15}{'in-band top75':>15}"
          f"{'ratio in/glob':>15}")
        for nm in sorted(gk, key=lambda k: -np.mean([r["global"][k] for r in allrows])):
            g = float(np.mean([r["global"][nm] for r in allrows]))
            i25 = float(np.mean([r["band"]["top25"]["sig"][nm]["rho"]
                                 for r in allrows if nm in r["band"]["top25"]["sig"]]))
            i75 = float(np.mean([r["band"]["top75"]["sig"][nm]["rho"]
                                 for r in allrows if nm in r["band"]["top75"]["sig"]]))
            P(f"   {nm:<18}{g:>13.3f}{i25:>15.3f}{i75:>15.3f}"
              f"{(i75 / g if abs(g) > 0.05 else float('nan')):>15.2f}")
        P("   A large global rho with a small in-band rho is a GARBAGE FILTER, which is what the")
        P("   shipped objective is.  A feature whose in-band rho is comparable to its global rho")
        P("   is a different kind of object, whatever its magnitude.")

    for b in BANDS:
        rows = [r for r in allrows if r["band"][b].get("sig")]
        fold = np.array([r["fold"] for r in rows], int)
        bm = np.array([r["band"][b]["mean"] for r in rows], float)
        alpha = _get(rows, b, "NULL_alpha", "sel")
        beta = _get(rows, b, "NULL_beta", "sel")
        P("")
        P("-" * 112)
        P(f"BAND = {b}   {LABEL[b]}   n = {len(rows)} targets")
        P("-" * 112)
        P("B/C. in-band skill on all three axes, then the SELECTED RMSD against the three bars.")
        P(f"   {'feature':<16}{'rhoS':>8}{'95% CI':>17}{'acc*':>7}{'rho_c':>7}"
          f"{'sel':>8}{'vs RANDOM [95% CI]':>26}{'W/L':>8}{'vs a-helix':>12}{'vs b-strand':>13}")
        for fam, keys in FAMILY:
            ks = [k for k in keys if k in have]
            if not ks:
                continue
            P(f"  {fam}")
            for nm in ks:
                #: a band-restricted feature (F6) exists on a subset of targets; every arm in its
                #: row is then restricted to the SAME subset so the pairing stays exact.
                keep = [q for q, r in enumerate(rows) if nm in r["band"][b]["sig"]]
                if len(keep) < 8:
                    continue
                rs = [rows[q] for q in keep]
                bm_, alpha_, beta_ = bm[keep], alpha[keep], beta[keep]
                fold_ = fold[keep]
                rho = _get(rs, b, nm, "rho")
                acc = _get(rs, b, nm, "acc")
                sel = _get(rs, b, nm, "sel")
                st = L.boot_target(rho)
                pr = L.report_pair(nm, sel, bm_, fold_)
                pa = L.report_pair(nm, sel, alpha_, fold_)
                pb = L.report_pair(nm, sel, beta_, fold_)
                tag = "" if len(keep) == len(rows) else f"  n={len(keep)}"
                rc, _ = F.axes(rho.mean())
                star = "*" if (st["lo"] > 0 or st["hi"] < 0) else " "
                sig = "*" if pr["fold_lo"] > 0 or pr["fold_hi"] < 0 else " "
                P(f"   {star}{nm:<15}{rho.mean():+8.3f}"
                  f"  [{st['lo']:+.3f},{st['hi']:+.3f}]{acc.mean():>7.3f}{rc:>7.3f}"
                  f"{sel.mean():>8.3f}   {pr['diff']:+.3f} [{pr['lo']:+.3f},{pr['hi']:+.3f}]{sig}"
                  f"{pr['W']:>5}/{pr['L']}{pa['diff']:>+12.3f}{pb['diff']:>+13.3f}{tag}")
        P("   rhoS = MEASURED Spearman.  acc* = pairwise ordering accuracy COUNTED DIRECTLY")
        P("   (null 0.500), not converted.  rho_c = Gaussian-copula rho from rhoS via")
        P("   selection_theory.rho_from_spearman.  Three different axes (LEDGER L3).")
        P("   * on rhoS = target-level CI excludes zero.  * after the RANDOM interval = the")
        P("   FOLD-CLUSTERED interval on selected RMSD excludes zero.  Negative = better.")

        # ---------------------------------------------------------- D. deciding controls
        P("")
        P("D. THE CONTROLS THAT DECIDE THE LANE (paired, target-level and fold-clustered)")
        for a, c, why in (("env_esm", "env_onehot", "matched architecture, ESM block vs one-hot"),
                          ("env_esm", "env_const", "vs no sequence at all"),
                          ("ss_esm", "ss_onehot", "matched architecture, ESM block vs one-hot"),
                          ("pll_mask", "pll_bg", "sequence CONTEXT vs composition only"),
                          ("pll_mask", "blosum_sim", "learned vs substitution matrix"),
                          ("emb_cos", "blosum_sim", "embedding similarity vs substitution matrix"),
                          ("emb_cos", "NULL_alpha", "vs the zero-information bar"),
                          ("esmcon_bce", "dist_shipped", "ESM contact head vs the shipped score"),
                          ("risk_wcon", "risk_wone", "exogenous weight vs uniform")):
            keep = [q for q, r in enumerate(rows)
                    if a in r["band"][b]["sig"] and c in r["band"][b]["sig"]]
            if len(keep) < 8:
                continue
            rs = [rows[q] for q in keep]
            pr = L.report_pair(f"{a} - {c}", _get(rs, b, a, "sel"), _get(rs, b, c, "sel"),
                               fold[keep])
            P(f"   {L.fmt_pair(pr)}   [{why}]")

        # ------------------------------------------------- D2. rho contrasts (PAIRED)
        #: `core/predict.py` lines 144 and 612: the ESM contact map is an INPUT FEATURE of the
        #: shipped distogram.  So "the raw contact map orders the band better than the trained
        #: predictor built on it" is a claim about what training does to in-band content, and it
        #: needs a PAIRED test on per-target rho, not two means side by side.
        P("")
        P("D2. IN-BAND rho, paired per target (POSITIVE = the FIRST signal orders the band better;")
        P("    note the sign convention is opposite to the RMSD tables above, where lower wins)")
        for a, c, why in (("esmcon_agree", "cf_uniform", "*** THE DECIDING CONTROL: is F1 the "
                           "language model, or is it compactness? ***"),
                          ("esmcon_agree", "rg_raw", "vs raw radius of gyration"),
                          ("esmcon_topd", "cf_topd_uniform", "ESM's top pairs vs the distogram's"),
                          ("esmcon_agree", "dist_shipped", "raw ESM contact map vs the trained "
                           "distogram that CONSUMES it as an input feature"),
                          ("esmcon_agree", "NULL_alpha", "vs the zero-information bar"),
                          ("esmcon_agree", "NULL_chance", "vs a random score")):
            keep = [q for q, r in enumerate(rows)
                    if a in r["band"][b]["sig"] and c in r["band"][b]["sig"]]
            if len(keep) < 8:
                continue
            rs = [rows[q] for q in keep]
            d = _get(rs, b, a, "rho") - _get(rs, b, c, "rho")
            bt = L.boot_target(d)
            bf = L.boot_fold(d, fold[keep])
            P(f"   rho({a}) - rho({c})  {bt['mean']:+.3f} "
              f"[{bt['lo']:+.3f},{bt['hi']:+.3f}]  fold [{bf['lo']:+.3f},{bf['hi']:+.3f}]  "
              f"med {bt['median']:+.3f}  {bt['L']}W/{bt['W']}L   [{why}]")

        # ---------------------------------------------------------- E. min-of-N
        have = sorted(set.intersection(*[set(r["band"][b]["sig"]) for r in rows]))
        best = np.array([min(r["band"][b]["sig"][k]["sel"] for k in have) for r in rows])
        #: match V to the number of features actually compared; the artefact stores a ladder.
        ladder = sorted(int(v) for v in rows[0]["band"][b]["minN"])
        V = min(ladder, key=lambda v: abs(v - len(have)))
        nullV = np.array([r["band"][b]["minN"][str(V)] for r in rows])
        pr = L.report_pair("best-of-V real minus zero-info", best, nullV, fold)
        P("")
        P(f"E. MIN-OF-N CHECK (LEDGER L8), {len(have)} real features vs a V = {V} null")
        P(f"   ORACLE min over the {len(have)} REAL features   {best.mean():.3f}")
        P(f"   ZERO-INFORMATION min over {V} random picks from the same band"
          f"   {nullV.mean():.3f}")
        P(f"   {L.fmt_pair(pr)}")
        P("   Positive = the real per-target minimum LOSES to drawing the same number of")
        P("   candidates at random from the same band, i.e. there is no routing headroom.")

    # ------------------------------------------------------------------ G. stratification
    #: BRIEF s4: "a mean improvement carried by easy targets is not a result".  The same applies
    #: to a mean NULL -- it can hide a real effect on one stratum.  Reported for the lane's
    #: headline feature in the primary deployable band.
    P("")
    P("G. STRATIFICATION of the lane's headline feature (top-25 band), per BRIEF s4")
    rows = allrows
    fold = np.array([r["fold"] for r in rows], int)
    nres = np.array([r["n"] for r in rows], int)
    orc = np.array([r["oracle"] for r in rows], float)
    bm = np.array([r["band"]["top25"]["mean"] for r in rows], float)
    for nm in ("esmcon_agree", "pll_mask", "emb_cos"):
        if nm not in rows[0]["band"]["top25"]["sig"]:
            continue
        sel = np.array([r["band"]["top25"]["sig"][nm]["sel"]
                        if nm in r["band"]["top25"]["sig"] else np.nan for r in rows])
        rho = np.array([r["band"]["top25"]["sig"][nm]["rho"]
                        if nm in r["band"]["top25"]["sig"] else np.nan for r in rows])
        ok = ~np.isnan(sel)
        P(f"   {nm}   (n = {int(ok.sum())})")
        strata = [("length 9-11", nres <= 11), ("length 12-13", (nres >= 12) & (nres <= 13)),
                  ("length 14-16", nres >= 14),
                  ("pool best < 2.0 (easy)", orc < 2.0), ("pool best >= 2.0 (hard)", orc >= 2.0)]
        strata += [(f"fold {f}", fold == f) for f in sorted(set(fold.tolist()))]
        for lab, m in strata:
            m = m & ok
            if m.sum() < 5:
                continue
            d = sel[m] - bm[m]
            P(f"     {lab:<24}n={int(m.sum()):>4}  rho {np.nanmean(rho[m]):+.3f}   "
              f"sel-vs-random {d.mean():+.3f}   {int((d < 0).sum())}W/{int((d > 0).sum())}L")

    if out:
        open(out, "w").write("\n".join(buf))
    return buf


def ceiling(rho_s, band="25", K0="500"):
    """F. If a feature clears the bar, what would its in-band rho DELIVER?  The band curve,
    read at the copula rho implied by the measured Spearman.  A CEILING, not a model."""
    path = os.path.join(F.RESULTS, "sel_theory.json")
    if not os.path.exists(path):
        print("  band curve artefact absent (s17/results/sel_theory.json); ceiling not computed")
        return None
    rows = json.load(open(path))["rows"]
    rc = ST.rho_from_spearman(rho_s)
    grid = list(ST.RHOS)
    print(f"  measured in-band Spearman {rho_s:+.3f} -> copula rho {rc:+.3f} "
          f"-> pairwise accuracy {ST.accuracy_from_rho(rc):.3f}")
    for B in (5, 10, 25, 50, 100):
        c = [x["cells"][K0]["band"].get(str(B)) for x in rows]
        c = [y for y in c if y]
        if not c:
            continue
        lo = max([g for g in grid if g <= max(rc, 0.0)] or [0.0])
        hi = min([g for g in grid if g >= max(rc, 0.0)] or [1.0])
        vlo = np.mean([y["curve"][f"{lo:.2f}"] for y in c])
        vhi = np.mean([y["curve"][f"{hi:.2f}"] for y in c])
        w = 0.0 if hi == lo else (max(rc, 0.0) - lo) / (hi - lo)
        print(f"   shortlist B = {B:>3}   band best {np.mean([y['best'] for y in c]):.3f}   "
              f"delivered at this rho {vlo + w * (vhi - vlo):.3f}   "
              f"(band mean {np.mean([y['mean'] for y in c]):.3f})")


if __name__ == "__main__":
    report(out=os.path.join(F.RESULTS, "feat_inband.log"))
