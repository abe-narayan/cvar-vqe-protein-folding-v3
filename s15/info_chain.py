"""SPRINT 15, INFO, PART B.2 -- the DIRECT chain, proxy -> per-target sign -> emitted RMSD.

Sprint 14 reported chained correlations (proxy -> native Rg -> per-target skill) and
flagged the composition as an assumption.  This measures the whole chain end to end on
the quantity that actually matters -- the RMSD a selection operator emits.

The experiment, per enumerated target:

    objective   E_s(c) = s * Rg(c),  s in {-1, +1}         (the compactness axis itself,
                which s14/obj_resid showed the 4069-parameter learned model reduces to)
                and the same with Legacy in place of Rg.
    operator    top-100 by E, mean CA-RMSD  (the terminal operator sprint 14 priced at
                1.26 A against a decile)
    arms        random          no objective: the space mean
                fixed sign      the sign that wins on the OTHER targets (leave-one-out)
                proxy sign      the sign a NATIVE-FREE proxy picks (threshold also LOO)
                ORACLE sign     the better of the two signs for this target -- the CEILING
                random sign     coin flip, averaged exactly over both signs

Everything except the `ORACLE sign` arm and the RMSD labels is inference-legal.

    python -m s15.info_chain
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s15 import info_lib as L           # noqa: E402

TOPM = 100


def emit(score, rmsd, m=TOPM):
    """Mean RMSD of the m lowest-scoring configurations, ties averaged over the tied set.

    Tie handling matters: `np.argsort` on a tied signal reads the array's own order, which
    on an RMSD-sorted array leaks the oracle (project finding).  Rg is continuous so ties
    are rare, but Legacy is not, so the cut is taken at the m-th value and the tied block
    is averaged rather than truncated arbitrarily.
    """
    score = np.asarray(score, float); rmsd = np.asarray(rmsd, float)
    thr = np.partition(score, m - 1)[m - 1]
    strict = score < thr
    ns = int(strict.sum())
    tied = score == thr
    if ns >= m:
        return float(rmsd[strict].mean())
    w = np.zeros(len(score))
    w[strict] = 1.0
    w[tied] = (m - ns) / float(tied.sum())
    return float((w * rmsd).sum() / w.sum())


def proxies_native_free():
    """Per-target native-free scalars, 126 targets where available."""
    with open(os.path.join(ROOT, "s14", "results", "signpred.json")) as fh:
        sp = json.load(fh)
    return {r["pdb"]: r for r in sp["per_target"]}


def build():
    sp = proxies_native_free()
    rows = []
    for p in L.enum_targets():
        g = L.geom_cache(p)
        rg = np.asarray(g["u_rg"], float)
        leg = np.asarray(g["u_legacy"], float)
        rr = np.asarray(g["u_rmsd"], float)
        r = {"pdb": p, "n": int(g["n"]), "fold": int(g["fold"]),
             "space_mean": float(rr.mean()), "space_min": float(rr.min()),
             "emit_rg_plus": emit(rg, rr),      # prefer COMPACT (low Rg)
             "emit_rg_minus": emit(-rg, rr),    # prefer EXTENDED (high Rg)
             "emit_leg_plus": emit(leg, rr),    # Legacy as published
             "emit_leg_minus": emit(-leg, rr),
             # the same two axes RESTRICTED TO THE NEAR-NATIVE BAND (ORACLE DIAGNOSTIC:
             # band membership reads the native).  This is the region sprint 14 showed
             # every objective is at chance inside, and where the 2.0 A answer lives.
             "bemit_rg_plus": emit(np.asarray(g["b_rg"], float),
                                   np.asarray(g["b_rmsd"], float)),
             "bemit_rg_minus": emit(-np.asarray(g["b_rg"], float),
                                    np.asarray(g["b_rmsd"], float)),
             "bemit_leg_plus": emit(np.asarray(g["b_legacy"], float),
                                    np.asarray(g["b_rmsd"], float)),
             "bemit_leg_minus": emit(-np.asarray(g["b_legacy"], float),
                                     np.asarray(g["b_rmsd"], float)),
             "band_mean": float(np.asarray(g["b_rmsd"], float).mean()),
             "band_min": float(np.asarray(g["b_rmsd"], float).min()),
             "rg_space_mean": float(rg.mean()), "rg_space_sd": float(rg.std()),
             "rg_nat_ORACLE": float(g["rg_native_ORACLE"]),
             "z_space_ORACLE": float((g["rg_native_ORACLE"] - rg.mean()) / rg.std())}
        s = sp[p]
        # native-free proxies, expressed on the SAME z scale as the oracle variable:
        # the space Rg mean/sd are native-free (they are properties of the state space).
        for k in ("rg_disto", "rg_pool", "rg_emit"):
            r["z_" + k] = float((s[k] - rg.mean()) / rg.std())
            r[k] = float(s[k])
        rows.append(r)
    return rows


def loo_sign(vals, better_minus):
    """Leave-one-out fixed-sign rule: use the majority sign of the OTHER targets."""
    out = []
    for i in range(len(vals)):
        o = [better_minus[j] for j in range(len(vals)) if j != i]
        out.append(-1 if np.mean(o) > 0.5 else +1)
    return np.asarray(out)


def loo_proxy_sign(x, better_minus):
    """Leave-one-out proxy rule: threshold at the median of the OTHER targets, and take
    the orientation that wins on the OTHER targets.  Nothing about target i is used."""
    x = np.asarray(x, float)
    n = len(x)
    out = np.empty(n, int)
    for i in range(n):
        idx = [j for j in range(n) if j != i]
        thr = np.median(x[idx])
        hi = np.asarray([x[j] > thr for j in idx])
        bm = np.asarray([better_minus[j] for j in idx])
        # orientation: does "above threshold" mean the MINUS sign is better?
        agree = (hi == bm).mean()
        pred_hi = x[i] > thr
        out[i] = (-1 if pred_hi else +1) if agree >= 0.5 else (+1 if pred_hi else -1)
    return out


def run_axis(rows, plus_key, minus_key, name, base_key="space_mean"):
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows])
    ep = np.asarray([r[plus_key] for r in rows], float)
    em = np.asarray([r[minus_key] for r in rows], float)
    base = np.asarray([r[base_key] for r in rows], float)
    better_minus = em < ep
    orac = np.minimum(ep, em)
    coin = 0.5 * (ep + em)
    fixed = np.where(loo_sign(ep, better_minus) < 0, em, ep)
    out = {"axis": name, "n": len(rows),
           "space_mean": float(base.mean()), "oracle_sign": float(orac.mean()),
           "coin_sign": float(coin.mean()), "fixed_sign_LOO": float(fixed.mean()),
           "always_plus": float(ep.mean()), "always_minus": float(em.mean()),
           "frac_minus_better": float(better_minus.mean()), "arms": {}}
    out["oracle_vs_fixed"] = L.report(orac, fixed, folds, pdbs, "ORACLE sign - fixed sign")
    out["fixed_vs_coin"] = L.report(fixed, coin, folds, pdbs, "fixed - coin")
    for pk in ("z_rg_disto", "z_rg_pool", "z_rg_emit"):
        x = np.asarray([r[pk] for r in rows], float)
        sg = loo_proxy_sign(x, better_minus)
        ev = np.where(sg < 0, em, ep)
        out["arms"][pk] = {
            "emitted": float(ev.mean()),
            "sign_accuracy": float((np.where(sg < 0, True, False) == better_minus).mean()),
            "vs_fixed": L.report(ev, fixed, folds, pdbs, pk + " - fixed"),
            "vs_oracle": L.report(ev, orac, folds, pdbs, pk + " - ORACLE")}
    # ORACLE compactness itself, as the ceiling on any compactness-based sign rule
    x = np.asarray([r["z_space_ORACLE"] for r in rows], float)
    sg = loo_proxy_sign(x, better_minus)
    ev = np.where(sg < 0, em, ep)
    out["arms"]["z_space_ORACLE"] = {
        "emitted": float(ev.mean()),
        "sign_accuracy": float((np.where(sg < 0, True, False) == better_minus).mean()),
        "vs_fixed": L.report(ev, fixed, folds, pdbs, "ORACLE compactness - fixed"),
        "vs_oracle": L.report(ev, orac, folds, pdbs, "ORACLE compactness - ORACLE sign")}
    return out


def main():
    rows = build()
    res = {"per_target": rows, "axes": {}}
    for pk, mk, nm, bk in (
            ("emit_rg_plus", "emit_rg_minus", "compactness axis (+/- Rg)", "space_mean"),
            ("emit_leg_plus", "emit_leg_minus", "Legacy (+/- legacy total)", "space_mean"),
            ("bemit_rg_plus", "bemit_rg_minus",
             "IN-BAND compactness axis (+/- Rg, RMSD<=2.5)", "band_mean"),
            ("bemit_leg_plus", "bemit_leg_minus",
             "IN-BAND Legacy (+/- legacy, RMSD<=2.5)", "band_mean")):
        res["axes"][nm] = run_axis(rows, pk, mk, nm, bk)
    L.jwrite("info_chain", res)

    print("PART B.2 -- THE DIRECT CHAIN: native-free proxy -> per-target sign -> "
          f"emitted RMSD\n  operator = top-{TOPM} of a 60,000-config uniform sample, "
          f"n = {len(rows)} enumerated targets\n")
    for nm, a in res["axes"].items():
        print(f"--- {nm} ---")
        print(f"  random selection (region mean)    {a['space_mean']:.3f} A")
        print(f"  always +  (as published)          {a['always_plus']:.3f} A")
        print(f"  always -  (sign flipped)          {a['always_minus']:.3f} A")
        print(f"  coin-flip sign (exact average)    {a['coin_sign']:.3f} A")
        print(f"  FIXED sign, leave-one-out         {a['fixed_sign_LOO']:.3f} A")
        print(f"  ORACLE sign (the ceiling)         {a['oracle_sign']:.3f} A")
        print(f"  fraction of targets where '-' wins: {a['frac_minus_better']:.3f}")
        print(f"  ORACLE - fixed : {L.fmt_paired(a['oracle_vs_fixed'])}")
        print(f"  fixed  - coin  : {L.fmt_paired(a['fixed_vs_coin'])}")
        print(f"  {'proxy':<22}{'emitted':>9}{'signacc':>9}   vs fixed")
        for pk, v in a["arms"].items():
            print(f"  {pk:<22}{v['emitted']:>9.3f}{v['sign_accuracy']:>9.3f}   "
                  f"{L.fmt_paired(v['vs_fixed'])}")
        print()
    return res


if __name__ == "__main__":
    main()
