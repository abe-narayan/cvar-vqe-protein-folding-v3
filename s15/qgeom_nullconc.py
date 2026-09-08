"""SPRINT 15 / QGEOM -- the NULL-CALIBRATED concentration check the brief mandates.

The brief's rule: *a raw drop-top threshold is NOT a valid concentration test. When mean/sd is
small, discarding the ten best removes much of the total EVEN IF every target carries an
identical effect. Compare every concentration statistic to a SIMULATED UNIFORM-EFFECT NULL,
and print mean/sd so a reader can see when the test has no power. Emit the verdict as one
PASS/FAIL, never as fields a reader can select from.*

This applies that rule to the only positive VQE result in this workstream -- the unranked
ensemble-generator comparison of `qgeom_ens.E2` -- and to the CVaR control test of
`qgeom_cvar.D2`.

    python -m s15.qgeom_nullconc
"""
from __future__ import annotations

import json
import os

import numpy as np

from s15 import qgeom_lib as G

TAG = "nullconc"


def null_calibrated(d, sims=4000, seed=0):
    """One PASS/FAIL verdict for whether a paired effect is CONCENTRATED.

    Simulates `d_i = mu + sigma z_i` with the OBSERVED mu and sigma -- a uniform effect plus
    the observed noise -- and asks where the observed drop-top statistics sit in that null.
    """
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    n = d.size
    if n < 6:
        return {"n": n, "verdict": "TOO SMALL"}
    mu, sd = float(d.mean()), float(d.std(ddof=1))
    k = max(1, n // 10)
    obs_share = float(np.sort(np.abs(d))[::-1][:k].sum() / max(np.abs(d).sum(), 1e-30))
    keep = np.sort(d)
    drop = float(np.mean(keep[: n - k]) if mu < 0 else np.mean(keep[k:]))
    r = np.random.default_rng(seed)
    z = r.standard_normal((sims, n))
    sim = mu + sd * z
    sim_share = np.sort(np.abs(sim), axis=1)[:, ::-1][:, :k].sum(1) / \
        np.maximum(np.abs(sim).sum(1), 1e-30)
    sim_sorted = np.sort(sim, axis=1)
    sim_drop = (sim_sorted[:, : n - k].mean(1) if mu < 0 else sim_sorted[:, k:].mean(1))
    p_share = float((sim_share >= obs_share).mean())
    p_drop = float((np.abs(sim_drop) <= np.abs(drop)).mean())
    concentrated = (p_share < 0.05) or (p_drop < 0.05)
    return {"n": n, "mean": mu, "sd": sd, "mean_over_sd": mu / sd if sd else np.inf,
            "top10pct_share": obs_share,
            "top10pct_share_null_median": float(np.median(sim_share)),
            "p_share_vs_null": p_share,
            "drop_top10pct_mean": drop,
            "drop_top10pct_null_median": float(np.median(sim_drop)),
            "p_drop_vs_null": p_drop,
            "power_warning": bool(abs(mu / sd) < 0.5) if sd else False,
            "verdict": "CONCENTRATED (FAIL)" if concentrated else "DIFFUSE (PASS)"}


def main():
    ens = json.load(open(os.path.join(G.RESULTS, "qgeom_ens.json")))["E2_generator"]
    cv = json.load(open(os.path.join(G.RESULTS, "qgeom_cvar.json")))["D2_control_sub12"]
    print("=" * 112)
    print("NULL-CALIBRATED CONCENTRATION CHECK -- one PASS/FAIL verdict per comparison")
    print("  PASS = the effect is spread across cells (a uniform-effect null explains the")
    print("  drop-top statistics).  FAIL = it is carried by a few cells and must not be")
    print("  reported as a general effect.  `mean/sd` below 0.5 means the test has little")
    print("  power and the PASS is weak evidence, not strong evidence.")
    print("=" * 112)
    out = {}
    jobs = []
    for ro in ("rand5_coordavg_rmsd", "rand20_coordavg_rmsd", "rand75_coordavg_rmsd",
               "set_coordavg_rmsd", "set_mean_rmsd", "set_best_rmsd"):
        for a in (1.0, 0.25, 0.05, 0.01):
            v = [x["vqe"][ro] - x["control"][ro] for k, x in ens.items()
                 if f"|a{a}|" in k and ro in x["vqe"]]
            if v:
                jobs.append((f"E2 {ro}|a{a}", v))
    for ro in ("argmin_rmsd", "top20_coordavg_rmsd", "top75_coordavg_rmsd"):
        for a in (1.0, 0.25, 0.05, 0.01):
            v = [x["vqe"][ro] - x["control"][ro] for k, x in cv.items()
                 if f"|a{a}|" in k and ro in x["vqe"]]
            if v:
                jobs.append((f"D2 {ro}|a{a}", v))
    print(f"{'comparison':38s} {'n':>3s} {'mean':>9s} {'mean/sd':>8s} "
          f"{'share':>7s} {'null':>7s} {'p':>6s} {'p_drop':>7s} {'verdict':>21s} pow")
    for name, v in jobs:
        r = null_calibrated(v)
        out[name] = r
        if r.get("verdict") == "TOO SMALL":
            continue
        print(f"{name:38s} {r['n']:3d} {r['mean']:+9.4f} {r['mean_over_sd']:+8.2f} "
              f"{r['top10pct_share']:7.3f} {r['top10pct_share_null_median']:7.3f} "
              f"{r['p_share_vs_null']:6.3f} {r['p_drop_vs_null']:7.3f} "
              f"{r['verdict']:>21s} {'LOW' if r['power_warning'] else ''}")
    G.ck(TAG, "null_calibrated", out)
    print()
    print("  `share` is the fraction of the total |effect| carried by the top 10% of cells;")
    print("  `null` is its median under a uniform-effect-plus-observed-noise simulation.")
    print("  Reading `share` alone against 1/n is the error the brief warns about.")
    return out


if __name__ == "__main__":
    main()
