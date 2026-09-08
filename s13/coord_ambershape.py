"""SPRINT 13, COORDINATOR EXPERIMENT 2 -- the shape of the AMBER objective in torsion space.

`s13/coord_objval.py` found, on one target, that AMBER ff14SB/GBn2 single-point energy over
random configurations of the discrete torsion space has a mean of ~1e16 kcal/mol and a rank
correlation with CA-RMSD of +0.011.  It separates the native from garbage trivially and then
carries no information about which non-native configuration is better.  Legacy, on the same
configurations, gives rho = +0.32 with a sane dynamic range.

If that generalises it is the single most important fact about the AMBER arm of this sprint,
and it has to be characterised properly rather than worked around, because the honest version
of the sprint's headline question is not "is AMBER worse" but:

    at which OPERATING POINT is an all-atom force field a usable variational objective,
    and what does moving to that operating point cost?

So this module measures the objective's SHAPE, per target, for five genuine variants:

    raw        ff14SB/GBn2 single-point on the ideal-geometry build.  Unmodified.
    capped     the same energy clipped at a percentile of its own distribution.  A MODIFIED
               POTENTIAL -- labelled as such everywhere, never called "AMBER".
    softcore   the same, with the steric blow-up tamed by a log transform of the positive
               part: sign-preserving, monotone, so it changes the SCALE not the ORDER.
               Also a modified potential.
    nonclash   electrostatics + solvation components only, i.e. AMBER minus the term that
               explodes.  The closest analogue to what Legacy actually models.
    minimised  energy after the production restrained minimisation (`refine_coords`,
               k=10, steps=0).  Genuine AMBER at the operating point the pipeline ships,
               and ~1000x the cost, so it runs on a small subset only.

Reported per variant: dynamic range, rho with CA-RMSD, the native's percentile, and -- the
quantity that actually decides usability -- rho computed over the LOW-ENERGY DECILE only,
because an optimiser never sees the bulk of the distribution.  An objective can look
informative on random configurations and be blind exactly where a search spends its time.

Native quantities are labels.  Every arm is DIAGNOSTIC.

    python -m s13.coord_ambershape [n_targets] [n_samples]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s13 import ceiling as CE              # noqa: E402
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)
K = 4
SEED = 20260905
N_MIN = 24            # how many configurations get the expensive minimised arm


def _wait_for_amber():
    import core.amber as AM
    for _ in range(90):
        if AM.memory_percent() < AM.MEMORY_LIMIT_PERCENT - 1.0:
            return True
        time.sleep(10.0)
    return False


def variants(seq, rep, S, want_min=True):
    """The five genuine variants for a batch of state assignments S (B, n)."""
    import core.amber as AM
    import core.energy as EN
    from core import geometry as geo
    B, n = S.shape
    idx = np.arange(n)
    PHI = np.asarray(rep._phi, float); PSI = np.asarray(rep._psi, float)
    phi = PHI[idx[None, :], S]; psi = PSI[idx[None, :], S]
    BB = geo.build_backbone_batch(phi, psi)

    lg = np.asarray(EN.totals_batch(EN.components_batch(seq, BB, phi=phi, psi=psi)), float)

    _wait_for_amber()
    raw = np.full(B, np.nan); nonclash = np.full(B, np.nan)
    for b in range(B):
        try:
            r = AM.single_point(seq, rep, S[b], components=True)
            raw[b] = float(r["energy"])
            c = r.get("components", {}) or {}
            #: electrostatics + implicit solvation only: AMBER minus the r^-12 term that
            #: blows up on an unrelaxed build.  Named for what it is.
            nonclash[b] = float(c.get("nonbonded_elec", c.get("electrostatic", 0.0))
                                + c.get("gb", c.get("solvation", 0.0)))
        except Exception:                                              # noqa: BLE001
            pass

    fin = np.isfinite(raw)
    capped = raw.copy()
    if fin.sum() > 4:
        hi = np.percentile(raw[fin], 90)
        capped = np.minimum(raw, hi)
    #: sign-preserving log compression of the positive part.  MONOTONE, so it cannot change
    #: the rank order -- it is a scale fix, and any rho difference against `raw` is purely a
    #: numerical-precision effect and should be ~0.  Reported precisely to prove that.
    softcore = np.where(raw > 0, np.log1p(np.maximum(raw, 0.0)), raw)

    mini = np.full(B, np.nan)
    if want_min:
        sel = np.arange(min(N_MIN, B))
        import torsion_lib2 as _tl
        for b in sel:
            try:
                cd = {k: v[b] for k, v in BB.items()}
                rr = AM.refine_coords(seq, rep, cd, k_restraint=10.0, steps=0,
                                      components=False)
                mini[b] = float(rr["energy"])
            except Exception:                                          # noqa: BLE001
                pass
    return {"legacy": lg, "raw": raw, "capped": capped, "softcore": softcore,
            "nonclash": nonclash, "minimised": mini}


def stats(E, rr, n_rand):
    from scipy.stats import spearmanr
    E = np.asarray(E, float); ok = np.isfinite(E)
    o = {"n_finite": int(ok[:n_rand].sum())}
    if ok[:n_rand].sum() < 8:
        return o
    e = E[:n_rand][ok[:n_rand]]; r = rr[:n_rand][ok[:n_rand]]
    s = spearmanr(e, r)
    o["rho_all"] = float(s.statistic); o["p"] = float(s.pvalue)
    o["dynamic_range"] = float(np.nanmax(e) - np.nanmin(e))
    o["log10_range"] = (float(np.log10(max(o["dynamic_range"], 1e-12))))
    o["median"] = float(np.median(e)); o["mad"] = float(np.median(np.abs(e - np.median(e))))
    #: THE DECIDING NUMBER: rho inside the low-energy decile, where a search actually lives.
    q = np.percentile(e, 10)
    m = e <= q
    if m.sum() >= 6:
        s2 = spearmanr(e[m], r[m])
        o["rho_low_decile"] = float(s2.statistic); o["n_low"] = int(m.sum())
        o["rmsd_low_decile_mean"] = float(r[m].mean())
    o["argmin_rmsd"] = float(r[int(np.argmin(e))])
    if np.isfinite(E[n_rand]):
        o["native_percentile"] = float(np.mean(e < E[n_rand]))
    if np.isfinite(E[n_rand + 1]):
        o["descent_percentile"] = float(np.mean(e < E[n_rand + 1]))
    return o


def run_target(t, n_samples=120, want_min=True):
    rng = np.random.default_rng(SEED + hash(t["pdb"]) % 99991)
    seq, n = t["seq"], t["n"]
    u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
    tab = tl2.library_for(seq, K, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    PHI = np.asarray(rep._phi, float); PSI = np.asarray(rep._psi, float)
    p = pdb.by_pdb(t["pdb"])
    s_snap = CE.snap(PHI, PSI, np.asarray(p.phi, float), np.asarray(p.psi, float))   # ORACLE
    s_desc, _ = CE.descent(PHI, PSI, nat, s_snap)                                    # ORACLE

    S = np.vstack([rng.integers(0, K, size=(n_samples, n)), s_snap[None, :], s_desc[None, :]])
    V = variants(seq, rep, S, want_min)
    W = I.build_ca(PHI[np.arange(n)[None, :], S], PSI[np.arange(n)[None, :], S])
    rr = I.kabsch_rmsd_batch(W, nat)                                                 # ORACLE

    return {"pdb": t["pdb"], "n": n, "fold": t["fold"],
            "rmsd_rand_mean": float(rr[:n_samples].mean()),
            "rmsd_rand_best": float(rr[:n_samples].min()),
            "rmsd_native_snap": float(rr[n_samples]), "rmsd_descent": float(rr[n_samples + 1]),
            "variants": {k: stats(v, rr, n_samples) for k, v in V.items()}}


def main(n_targets=24, n_samples=120):
    tg = I.targets()
    #: a length- and fold-spread subset, chosen by position not by performance
    pick = list(range(0, len(tg), max(1, len(tg) // n_targets)))[:n_targets]
    sel = [tg[i] for i in pick]
    path = os.path.join(RESULTS, "coord_amber_shape.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for i, t in enumerate(sel):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t, n_samples))
        json.dump({"what": "shape of the AMBER objective in discrete torsion space",
                   "k": K, "n_samples": n_samples, "per_target": rows},
                  open(path, "w"), indent=1)
        v = rows[-1]["variants"]
        print(f"  {i+1}/{len(sel)} {t['pdb']} rand={rows[-1]['rmsd_rand_mean']:.2f} | "
              + " ".join(f"{nm}:rho={v[nm].get('rho_all', float('nan')):+.2f}"
                         f"/lo={v[nm].get('rho_low_decile', float('nan')):+.2f}"
                         for nm in ("legacy", "raw", "nonclash", "minimised"))
              + f"  [{time.time()-t0:.0f}s]", flush=True)
    report(rows)


def report(rows):
    names = ("legacy", "raw", "capped", "softcore", "nonclash", "minimised")
    out = {"n_targets": len(rows), "variants": {}}
    print(f"\n{len(rows)} targets, k={K}. random-config mean RMSD "
          f"{np.mean([r['rmsd_rand_mean'] for r in rows]):.3f}, "
          f"native-snap {np.mean([r['rmsd_native_snap'] for r in rows]):.3f}\n")
    hdr = (f"{'variant':11s} {'rho_all':>8s} {'rho_low10':>10s} {'log10 range':>12s} "
           f"{'nat pct':>8s} {'desc pct':>9s} {'argmin RMSD':>12s}")
    print(hdr); print("-" * len(hdr))
    for nm in names:
        g = lambda key: np.array([r["variants"][nm].get(key, np.nan) for r in rows], float)  # noqa: E731
        o = {key: float(np.nanmean(g(key))) for key in
             ("rho_all", "rho_low_decile", "log10_range", "native_percentile",
              "descent_percentile", "argmin_rmsd", "rmsd_low_decile_mean")}
        o["rho_all_frac_positive"] = float(np.nanmean(g("rho_all") > 0))
        o["n_targets_with_data"] = int(np.isfinite(g("rho_all")).sum())
        out["variants"][nm] = o
        print(f"{nm:11s} {o['rho_all']:+8.3f} {o['rho_low_decile']:+10.3f} "
              f"{o['log10_range']:12.2f} {o['native_percentile']:8.3f} "
              f"{o['descent_percentile']:9.3f} {o['argmin_rmsd']:12.3f}")
    print("\nREAD: rho_low10 is rho INSIDE the low-energy decile -- where an optimiser "
          "actually lives. An objective with rho_all > 0 and rho_low10 ~ 0 is a garbage "
          "detector, not a ranker, and optimising it harder will not help.")
    with open(os.path.join(RESULTS, "coord_amber_shape_report.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    ns = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    main(nt, ns)
