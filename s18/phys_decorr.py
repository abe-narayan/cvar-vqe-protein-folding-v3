"""s18/phys_decorr.py -- THE GO / NO-GO ON THE LAMBDA LADDER, run before its result is read.

The coordinator's n = 126 result (`s18/COORD_FINDING.md`, ledger L2, `s18/results/objceil.json`)
re-poses my hypothesis and sharpens it:

  * with PERFECT distances the identical functional form reaches **1.152 A** -- the form is sound,
    so "the objective is missing inter-residue structure" is the wrong diagnosis;
  * destroying only the DIRECTION of the distogram's residual, keeping its magnitude, buys
    **-1.001 [-1.226, -0.784]** -- the error is structured and the structure is actively harmful.

> **So the sharp question for `leg_contact` is not whether it carries information. It is whether
> the direction it pushes each pair is CORRELATED WITH THE DISTOGRAM'S ERROR ON THAT PAIR.**
>   * pushes AGAINST the error  -> it corrects, and the ladder should work;
>   * DECORRELATED               -> it adds independent information, and the ladder is worth its
>                                   compute but has no reason to be large;
>   * pushes WITH the error      -> it amplifies the harm, and the ladder must fail. Stating that
>                                   in advance is worth more than discovering it afterwards.

THIS MEASUREMENT IS AN **ORACLE DIAGNOSTIC**.  It reads `d_true` to form the residual.  It is
never a selector, never a predictor, and no parameter anywhere in this workstream is chosen with
it.  It is a go/no-go on how to *read* the lambda ladder, and it is labelled ORACLE in every line
of output.

THE TWO QUANTITIES, defined so the sign is unambiguous.

    r_p    = dhat_p - dtrue_p                      the distogram's SIGNED error on pair p.
                                                   r_p > 0: the objective is aiming too FAR.
    need_p = -r_p                                  the correction that pair needs.
    push_p = -dE_contact/dd_p                      the contact term's DESCENT direction on that
                                                   pair's distance.  `contact_term` is
                                                   sum_p MJ_p * switch(d_p), switch' <= 0, so
                                                   push_p = MJ_p * |switch'(d_p)|:
                                                   MJ_p < 0 (favourable) pulls the pair TOGETHER.
    c_p    = MJ_p * switch(d_p)                    the coordinator's literal "per-pair
                                                   contribution", reported beside push_p because
                                                   the two answer slightly different questions.

    corr(push, need) > 0  -> CORRECTIVE
    corr(push, need) ~ 0  -> INDEPENDENT
    corr(push, need) < 0  -> AMPLIFYING

THE CONTROL.  The same statistic with the MJ table built on a `stable_rng` permutation of the
residue labels -- identical functional form, identical magnitudes, no sequence information.  A
correlation that the shuffled table reproduces is a property of the geometry, not of Legacy.

WHERE IT IS EVALUATED.  At the structure the refinement actually starts from -- the ideal-geometry
projection of the shipped coordinate average -- because that is where the two gradients meet.  It
is repeated at the NATIVE, where the contact term's push is what a correct potential would want,
and over the K = 500 pool for a distribution rather than a point.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                  # noqa: E402
from s14 import avgspace as AV                   # noqa: E402
from s15 import distcal as C                     # noqa: E402
from s15 import seed as SD                       # noqa: E402
from s18 import phys_lib as PL                   # noqa: E402

CONTACT_D0, CONTACT_DC = PL.CONTACT_D0, PL.CONTACT_DC


def switch_deriv(d):
    """d/dd of `core.energy.switch(d, 4.5, 8.5)`.  Zero outside the shell, <= 0 inside."""
    d = np.asarray(d, float)
    g = np.zeros_like(d)
    mid = (d > CONTACT_D0) & (d < CONTACT_DC)
    if np.any(mid):
        g[mid] = -0.5 * math.pi / (CONTACT_DC - CONTACT_D0) * np.sin(
            math.pi * (d[mid] - CONTACT_D0) / (CONTACT_DC - CONTACT_D0))
    return g


def _pair_map(ct, i, j, n):
    """Index of each CONTACT pair inside the OBJECTIVE's pair list, and the shared mask.

    The two pair sets differ: the objective uses `I.pair_index(n)` (min separation 2, CA) and the
    contact term uses |i-j| >= 3 over CB.  Only their intersection can be correlated, and the
    fraction shared is reported rather than assumed.
    """
    key = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(i, j))}
    idx = np.array([key.get((int(a), int(b)), -1) for a, b in zip(ct.di, ct.dj)], int)
    ok = idx >= 0
    return idx, ok


def _stats(push, c, need, r, sw):
    """Every correlation this module reports, on one structure.  Positive push-need = CORRECTIVE."""
    live = np.abs(sw) > 1e-12          # pairs the switch actually touches
    out = {"n_pairs": int(len(push)), "frac_live": float(live.mean()),
           "corr_push_need": PL.spearman(push, need),
           "corr_push_need_pearson": (float(np.corrcoef(push, need)[0, 1])
                                      if push.std() > 1e-12 and need.std() > 1e-12 else float("nan")),
           "corr_c_r": PL.spearman(c, r),
           "corr_absc_absr": PL.spearman(np.abs(c), np.abs(r))}
    if live.sum() >= 4:
        out["corr_push_need_live"] = PL.spearman(push[live], need[live])
        out["corr_c_r_live"] = PL.spearman(c[live], r[live])
    else:
        out["corr_push_need_live"] = float("nan")
        out["corr_c_r_live"] = float("nan")
    return out


def run_target(t, data, deb):
    from core import geometry as geo
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    d = data[pdb]
    i, j, nat = d["i"], d["j"], np.asarray(d["nat"], float)
    dtrue = np.asarray(d["dtrue"], float)
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
    r = dhat - dtrue                                    # ORACLE: the distogram's signed error
    need = -r

    ct = PL.ContactTerm(seq)
    ct0 = PL.ContactTerm(seq, mj=PL.shuffled_mj(seq, SD.stable_rng(pdb, "s18mjnull")))
    idx, ok = _pair_map(ct, i, j, n)

    W, PH, PS, u = AV.top75_windows(pdb)
    W = np.asarray(W, float)
    P = I.pairwise_rmsd(W)
    avg, _b = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    phi0 = np.asarray(pr["phi"], float)
    psi0 = np.asarray(pr["psi"], float)

    rec = {"pdb": pdb, "n": n, "fold": fold,
           "frac_contact_pairs_shared": float(ok.mean()),
           "resid_rms": float(np.sqrt((r ** 2).mean())),
           "resid_mae": float(np.abs(r).mean())}

    def at(nm, CB):
        dcb = np.linalg.norm(CB[ct.di] - CB[ct.dj], axis=1)
        from core import energy as et
        sw = np.asarray(et.switch(dcb, CONTACT_D0, CONTACT_DC), float)
        for tag, mjv in (("", ct.mj_pair), ("_null", ct0.mj_pair)):
            push = mjv * np.abs(switch_deriv(dcb))
            c = mjv * sw
            rec[nm + tag] = _stats(push[ok], c[ok], need[idx[ok]], r[idx[ok]], sw[ok])

    at("start", np.asarray(geo.build_backbone(phi0, psi0)["CB"], float))
    #: the NATIVE structure: what a correct potential's push would look like where the truth is.
    #: The native has no torsions in this instrument, so CB is taken from the CA trace's ideal
    #: rebuild of the native's own projected torsions -- recorded as an approximation, not hidden.
    prn = I.project(nat, seq, fold)
    at("native", np.asarray(geo.build_backbone(np.asarray(prn["phi"], float),
                                               np.asarray(prn["psi"], float))["CB"], float))

    #: the POOL: the same statistic averaged over the first 50 retrieval candidates, so the
    #: point estimate above is read against a distribution rather than one structure.
    m = min(50, len(PH))
    CBs = geo.build_backbone_batch(np.asarray(PH, float)[:m], np.asarray(PS, float)[:m])["CB"]
    vals, vals0 = [], []
    from core import energy as et
    for b in range(m):
        dcb = np.linalg.norm(CBs[b][ct.di] - CBs[b][ct.dj], axis=1)
        g = np.abs(switch_deriv(dcb))
        vals.append(PL.spearman((ct.mj_pair * g)[ok], need[idx[ok]]))
        vals0.append(PL.spearman((ct0.mj_pair * g)[ok], need[idx[ok]]))
    rec["pool_corr_push_need_mean"] = float(np.nanmean(vals))
    rec["pool_corr_push_need_sd"] = float(np.nanstd(vals))
    rec["pool_corr_push_need_null_mean"] = float(np.nanmean(vals0))
    return rec


def run(targets=None, out="decorr.json", verbose=True):
    tg = targets if targets is not None else I.targets()
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    deb = {}
    for f in sorted({int(t["fold"]) for t in tg}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))
    rows, t0 = [], time.time()
    for t in tg:
        rows.append(run_target(t, data, deb))
        if verbose and len(rows) % 20 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
    PL.write(out, {"rows": rows, "label": "ORACLE DIAGNOSTIC -- reads d_true; never a selector"},
             n_expected=len(tg))
    report(rows)
    return rows


def report(rows=None, out="decorr.json"):
    import json
    if rows is None:
        rows = json.load(open(os.path.join(PL.RESULTS, out)))["rows"]
    folds = np.array([r["fold"] for r in rows], int)
    names = [r["pdb"] for r in rows]
    print(f"\nORACLE DIAGNOSTIC -- n = {len(rows)}.  Reads d_true to form the residual; it is a "
          f"go/no-go on how to READ the lambda ladder and is a selector for nothing.\n")
    print(f"  contact pairs also in the objective's pair list: "
          f"{np.mean([r['frac_contact_pairs_shared'] for r in rows]):.3f}")
    print(f"  distogram residual RMS {np.mean([r['resid_rms'] for r in rows]):.3f} A, "
          f"MAE {np.mean([r['resid_mae'] for r in rows]):.3f} A")
    print("\n  POSITIVE = the contact term pushes AGAINST the distogram's error (CORRECTIVE)")
    print("  NEGATIVE = it pushes WITH the error (AMPLIFYING);  ~0 = INDEPENDENT\n")
    #: NOTE ON THE COUNTS.  `PL.paired` counts W as (a - b < 0), which for an RMSD means "better".
    #: These rows are CORRELATIONS, where HIGHER is more corrective, so the columns are printed as
    #: `n_hi / n_lo` -- the number of targets on which the real MJ table is MORE / LESS corrective
    #: than its shuffled null.  Reading them as W/L would invert the sign of every conclusion.
    print(f"  {'where':<22}{'stat':<26}{'mean':>9}{'median':>9}"
          f"{'vs its MJ-shuffled null [95% fold CI]':>40}{'n_hi/n_lo':>12}")
    for where in ("start", "native"):
        for st in ("corr_push_need", "corr_push_need_live", "corr_c_r"):
            a = np.array([r[where].get(st, np.nan) for r in rows], float)
            b = np.array([r[where + "_null"].get(st, np.nan) for r in rows], float)
            p = PL.paired(a, b, folds=folds, names=names, seed=7)
            ci = p.get("ci_fold", p["ci"])
            print(f"  {where:<22}{st:<26}{np.nanmean(a):>9.3f}{np.nanmedian(a):>9.3f}"
                  f"   {p['mean']:+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}]{p['L']:>8}/{p['W']}")
    a = np.array([r["pool_corr_push_need_mean"] for r in rows], float)
    b = np.array([r["pool_corr_push_need_null_mean"] for r in rows], float)
    p = PL.paired(a, b, folds=folds, names=names, seed=7)
    ci = p.get("ci_fold", p["ci"])
    print(f"  {'pool (50 cands)':<22}{'corr_push_need':<26}{np.nanmean(a):>9.3f}"
          f"{np.nanmedian(a):>9.3f}   {p['mean']:+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}]"
          f"{p['L']:>8}/{p['W']}")
    print("\nREAD.  If `corr_push_need` is at zero and indistinguishable from its MJ-shuffled")
    print("null, `leg_contact` supplies NO directional correction to the distogram's error, and")
    print("the lambda ladder cannot recover the 1.001 A the coordinator's shuffled control buys.")
    print("That would not falsify the ladder -- an INDEPENDENT signal is still a signal -- but it")
    print("bounds what the ladder can be expected to do, and it is stated before the ladder is read.")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        tg = I.targets()
        if a.limit:
            tg = tg[:a.limit]
        run(tg)
