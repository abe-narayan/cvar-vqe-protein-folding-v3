"""SPRINT 14, ENER-7 -- NORMALISATION CONSTANTS, derived rather than searched.

Raw AMBER spans 16 orders of magnitude on this space.  Before anything is combined the
terms have to be on comparable scales, and the scale has to be derived from a property of
the objective rather than grid-searched.  Three candidates, in increasing relevance to an
optimiser:

    energy variance     sd of the objective over a native-free reference population.
                        INVALID for raw AMBER: 40-68% of ideal-geometry configurations
                        exceed 1e4 kcal/mol (Sprint 13), so the sample variance is set by
                        the worst clash and is not an estimate of anything.

    gradient variance   the variance of a SINGLE-FLIP difference,
                            G(f) = (1/n) sum_i E_{s_-i}[ Var_{s_i}( f ) ],
                        computed EXACTLY from the full enumeration.  This is the scale an
                        optimiser actually experiences: it is the typical energy change of
                        the only move the search can make.  Robust to a delta spike in a
                        way the global variance is not, because a single catastrophic
                        configuration contributes to one flip-neighbourhood, not to all.

    rank               scale-free by construction, finite on any distribution, and
                        MONOTONE -- so it cannot change any ranking metric.

The last point is the load-bearing one.  The brief records that a monotone log compression
leaves rho at -0.036 while collapsing the range from 16.1 decades to 1.8.  That is not a
curiosity: it PROVES the damage is in the ordering.  Any proposed fix that only changes the
scale is treating the wrong problem, and this module verifies the invariance rather than
assuming it, on every metric this sprint uses.

Constants are written to `s14/cache/ener_norm.json` for the other workstreams.

    python -m s14.ener_norm
"""
from __future__ import annotations

import json
import os

import numpy as np

from s14 import ener_lib as E


def grad_var(f, n, k):
    """Exact mean single-flip variance: the scale a local search experiences."""
    f = np.asarray(f, float)
    tot = 0.0
    for i in range(n):
        a = f.reshape(k ** i, k, k ** (n - 1 - i))
        tot += a.var(axis=1).mean()
    return float(tot / n)


def constants(pdb):
    z = E.enum(pdb)
    n, k = z.n, z.k
    out = {}
    fields = {"legacy": z.legacy, "prior": z.prior,
              "legacy_nosteric": z.obj("legacy_nosteric")}
    for t in E.LEG_TERMS:
        fields["leg_" + t] = z.leg[t]
    for nm, v in fields.items():
        gv = grad_var(v, n, k)
        out[nm] = dict(mean=float(v.mean()), sd=float(v.std()),
                       median=float(np.median(v)),
                       mad=float(np.median(np.abs(v - np.median(v)))),
                       grad_var=gv, grad_sd=float(np.sqrt(gv)),
                       range_decades=_decades(v), exact=True)
    m = z.uniform_mask
    at = z.amber_total[m]
    out["amber"] = dict(mean=float(at.mean()), sd=float(at.std()),
                        median=float(np.median(at)),
                        mad=float(np.median(np.abs(at - np.median(at)))),
                        range_decades=_decades(at), exact=False,
                        note="uniform stratum only; sd is not an estimate of anything")
    out["amber_softcore"] = dict(sd=float(E.softcore(at).std()),
                                 median=float(np.median(E.softcore(at))),
                                 mad=float(np.median(np.abs(E.softcore(at)
                                                            - np.median(E.softcore(at))))),
                                 range_decades=_decades(E.softcore(at)), exact=False)
    return out


def _decades(v):
    v = np.abs(np.asarray(v, float))
    v = v[np.isfinite(v) & (v > 0)]
    if len(v) < 2:
        return 0.0
    return float(np.log10(v.max()) - np.log10(v.min()))


def monotone_invariance(pdb):
    """VERIFY, do not assume: a monotone transform changes no ranking metric."""
    from s14.ener_matrix import ensemble
    d = ensemble(pdb)
    R = d["rmsd"]
    rows = {}
    for nm, v in (("amber", d["amber"]), ("legacy", d["legacy"])):
        variants = {
            "raw": v,
            "softcore": E.softcore(v),
            "rank": E.rank_norm(v),
            "robust_z": E.robust_z(v),
            "shift+scale": 3.7 * v - 11.0,
            "cap_p90": np.minimum(v, np.percentile(v, 90)),          # NOT monotone
            "cap_p99": np.minimum(v, np.percentile(v, 99)),          # NOT monotone
        }
        rows[nm] = {}
        for tn, tv in variants.items():
            rows[nm][tn] = dict(
                rho=E.spearman(tv, R), rho_decile=E.decile_rho(tv, R),
                sel_rmsd=E.argmin_rmsd(tv, R),
                pair=E.pair_accuracy(tv, R, 100000, np.random.default_rng(1)),
                decades=_decades(tv))
    return rows


def main():
    per = {p: constants(p) for p in E.ENUM_TARGETS}
    fields = list(per[E.ENUM_TARGETS[0]])
    print("=== NORMALISATION CONSTANTS (pooled over the nine enumerated targets) ===")
    print(f"{'objective':22s} {'sd':>12s} {'MAD':>12s} {'grad_sd':>12s} "
          f"{'sd/grad_sd':>11s} {'decades':>9s}")
    pooled = {}
    for nm in fields:
        r = [per[p][nm] for p in E.ENUM_TARGETS]
        sd = float(np.mean([x["sd"] for x in r if "sd" in x]))
        mad = float(np.mean([x["mad"] for x in r if "mad" in x]))
        gs = float(np.mean([x["grad_sd"] for x in r])) if "grad_sd" in r[0] else None
        dec = float(np.mean([x["range_decades"] for x in r]))
        pooled[nm] = dict(sd=sd, mad=mad, grad_sd=gs, decades=dec,
                          exact=bool(r[0].get("exact", False)))
        print(f"{nm:22s} {sd:12.4g} {mad:12.4g} "
              f"{(f'{gs:12.4g}' if gs else '         n/a')} "
              f"{(f'{sd/gs:11.3f}' if gs and gs > 0 else '        n/a')} {dec:9.2f}")

    print("\n=== MONOTONE INVARIANCE -- VERIFIED, NOT ASSUMED ===")
    inv = [monotone_invariance(p) for p in E.ENUM_TARGETS]
    for nm in ("amber", "legacy"):
        print(f"\n-- {nm} --")
        print(f"{'transform':14s} {'rho':>8s} {'decile':>8s} {'pair':>7s} "
              f"{'sel_rmsd':>9s} {'decades':>9s}  monotone?")
        for tn in inv[0][nm]:
            r = [x[nm][tn] for x in inv]
            mono = tn not in ("cap_p90", "cap_p99")
            print(f"{tn:14s} {np.nanmean([x['rho'] for x in r]):+8.4f} "
                  f"{np.nanmean([x['rho_decile'] for x in r]):+8.4f} "
                  f"{np.nanmean([x['pair'] for x in r]):7.4f} "
                  f"{np.mean([x['sel_rmsd'] for x in r]):9.3f} "
                  f"{np.mean([x['decades'] for x in r]):9.2f}  "
                  f"{'yes' if mono else 'NO'}")

    path = os.path.join(E.CACHE, "ener_norm.json")
    with open(path, "w") as fh:
        json.dump(dict(
            what="derived normalisation constants for the s14 objectives",
            rule=("Use `grad_sd` (exact mean single-flip sd) to put terms on the scale an "
                  "optimiser experiences.  Do NOT use the sd of raw AMBER -- it is set by "
                  "the worst clash and is not an estimate.  For AMBER use rank "
                  "normalisation, which is scale-free and monotone."),
            pooled=pooled, per_target=per,
            monotone_invariance=inv), fh, indent=1,
            default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"\nwrote {path}")
    E.write("ener_norm", dict(what="normalisation constants and monotone-invariance check",
                              pooled=pooled, monotone_invariance=inv))
    return pooled, inv


if __name__ == "__main__":
    main()
