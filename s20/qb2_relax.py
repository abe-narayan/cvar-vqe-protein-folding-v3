"""SPRINT 20 / WORKSTREAM B -- THE RELAXATION ARM.  `H_AMBER = E o Relax` is a composition,
and this module measures what the composition contributes so it is not silently attributed to
"AMBER physics".

    python -m s20.qb2_relax run [n]
    python -m s20.qb2_relax report

WHY THIS EXISTS.  The coordinator's Q1 as first briefed compared `H_Legacy` (a bare energy) with
`H_AMBER` as implemented in `amber_hamiltonian.AmberHamiltonian.energy`, which is
build -> restrain (k = 100 kcal/mol/A^2 on every heavy atom) -> `LocalEnergyMinimizer.minimize`
(50 iterations) -> report the UNRESTRAINED energy at the MINIMISED coordinates.  Every landscape
metric measured on that object is confounded by the relaxation operator, which both SMOOTHS the
field and INJECTS DISCONTINUITIES wherever the minimiser's basin assignment flips.

THIS LANE'S PRIMARY Q1 ARM IS NOT CONFOUNDED.  `qb2_lib.AmberSP` calls
`core.amber._run(..., steps < 0)`, which is a **true single point with no minimisation at all**,
asserted bit-exact against `core.amber.refine_coords(k_restraint=0, steps=-1)` on every target.
Source check (BRIEF section 10): `_run` never calls `_is_collapsed` -- the ONLY call site is
`AmberHamiltonian._evaluate` at core/amber.py:1168 -- so the bare arm has **no infinite-valued
region by construction**, not merely an unfired guard.

So this module adds the missing arms rather than relabelling the existing ones:

    AMB      bare single point, no relaxation             (the un-confounded Q1 arm)
    AMBr1    E o Relax_1   (k = 100, steps = 1)           the smallest honest relaxation --
                                                          OpenMM reads maxIterations = 0 as
                                                          UNBOUNDED, so 1 is the floor, not 0
    AMBr50   E o Relax_50  (k = 100, steps = 50)          the DEPLOYED VQE Hamiltonian setting

and reports, per target: the rank agreement between the three (the audit lane's F-D2 test),
the dynamic range each one presents to an optimiser, the ruggedness of each, whether each RANKS
the pool against Ca-RMSD, and `n_collapsed` as a MEASURED count on every arm.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s20 import qb2_lib as L
from s20 import qb2_run as R
from s15 import seed as SD

N_REF = 120          # reference-ensemble size for the relaxed arms (declared; compute-bound)
N_GRAD_STARTS = 3
LINE_PTS = 65


def _panel(e, scale_med, scale_iqr):
    """The scale-free ruggedness indices of a 1-D scan, on the robust-standardised field."""
    es = (np.asarray(e, float) - scale_med) / scale_iqr
    lm = int(((es[1:-1] < es[:-2]) & (es[1:-1] < es[2:])).sum())
    x = es - es.mean()
    if x.std() <= 0:
        acl = float("nan")
    else:
        ac = np.correlate(x, x, "full")[len(x) - 1:]
        ac = ac / ac[0]
        below = np.where(ac < 1.0 / np.e)[0]
        acl = float(below[0] * (2 * np.pi / (len(es) - 1))) if len(below) else 2 * np.pi
    rg = float(es.max() - es.min())
    return {"n_localmin": lm, "acorr_len_rad": acl,
            "tv_over_range": float(np.abs(np.diff(es)).sum() / rg) if rg > 0 else float("nan"),
            "range_std": rg}


def _robust(e):
    ev = np.asarray(e, float)
    ev = ev[np.isfinite(ev)]
    q1, q2, q3 = np.percentile(ev, [25, 50, 75])
    return float(q2), float(max((q3 - q1) / 1.349, 1e-12))


def run(pdbs):
    tag = "relax"
    out = R.ck_load(tag)
    for ti, pdb in enumerate(pdbs):
        if pdb in out:
            continue
        R.mem_hold(tag=tag)
        t0 = time.time()
        tgt = L.target(pdb)
        n = tgt["n"]
        Z = R.starts_for(tgt)
        sp = L.AmberSP(tgt["seq"], tgt["rep"])
        r1 = L.AmberRelax(tgt["seq"], tgt["rep"], k_restraint=100.0, steps=1)
        r50 = L.AmberRelax(tgt["seq"], tgt["rep"], k_restraint=100.0, steps=50)
        ver = {"AMB": sp.verify(tgt["PHI"], tgt["PSI"], 3),
               "AMBr1": r1.verify(tgt["PHI"], tgt["PSI"], 2),
               "AMBr50": r50.verify(tgt["PHI"], tgt["PSI"], 2)}
        PHI = tgt["PHI"][:N_REF]; PSI = tgt["PSI"][:N_REF]
        rr = L.rmsd_of(L.pack(PHI, PSI), tgt)                       # ORACLE, post-hoc
        E = {"AMB": sp.batch(PHI, PSI), "AMBr1": r1.batch(PHI, PSI),
             "AMBr50": r50.batch(PHI, PSI)}
        E["AMBc"] = np.sign(E["AMB"]) * np.log1p(np.abs(E["AMB"]))
        # Legacy on the SAME reference ensemble, for the scale comparison
        hl = L.Ham("LEG", tgt)
        E["LEG"] = hl.raw(PHI, PSI)
        rec = {"n": n, "fold": tgt["fold"], "n_ref": int(N_REF),
               "verify": {k: [float(v[0]), int(v[1])] for k, v in ver.items()},
               "n_collapsed": {"AMB": 0, "AMBr1": int(r1.n_collapsed),
                               "AMBr50": int(r50.n_collapsed)},
               "n_calls": {"AMB": int(sp.n_calls), "AMBr1": int(r1.n_calls),
                           "AMBr50": int(r50.n_calls)},
               "obj": {}, "rank_agreement": {}}
        for k, e in E.items():
            fin = np.isfinite(e)
            med, iqr = _robust(e[fin])
            rec["obj"][k] = {
                "log10_range": float(np.log10(np.ptp(e[fin]) + 1e-30)),
                "sd_over_iqr": float(np.std(e[fin]) / iqr),
                "skew_p99_over_p50": float(np.percentile(np.abs(e[fin]), 99)
                                           / max(1e-30, np.percentile(np.abs(e[fin]), 50))),
                "med": med, "iqr": iqr,
                "spearman_E_vs_rmsd_ORACLE": float(L.spearman(e[fin], rr[fin])),
                "argmin_rmsd_ORACLE": float(rr[fin][int(np.argmin(e[fin]))]),
                "top10pct_mean_rmsd_ORACLE": float(
                    rr[fin][np.argsort(e[fin])[:max(1, int(0.10 * fin.sum()))]].mean()),
                "n_nonfinite": int((~fin).sum()),
            }
        # --- F-D2: does the relaxation change the ORDER, or only the conditioning?
        for a, b in (("AMBr50", "AMBr1"), ("AMBr50", "AMB"), ("AMBr1", "AMB"),
                     ("AMBr50", "AMBc"), ("AMB", "LEG"), ("AMBr50", "LEG")):
            rec["rank_agreement"][f"{a}|{b}"] = float(L.spearman(E[a], E[b]))
        # --- gradient and ruggedness for the two relaxed arms and the bare one
        rng = SD.stable_rng(pdb, "relax", salt=L.SALT)
        u = rng.normal(size=2 * n); u /= np.linalg.norm(u)
        ts = np.linspace(-np.pi, np.pi, LINE_PTS)
        for k, src in (("AMB", sp), ("AMBr1", r1), ("AMBr50", r50)):
            med, iqr = rec["obj"][k]["med"], rec["obj"][k]["iqr"]
            gs = []
            for r in range(N_GRAD_STARTS):
                d = 2 * n
                Zp = np.repeat(Z[r][None], 2 * d, 0)
                for c in range(d):
                    Zp[2 * c, c] += L.FD_H; Zp[2 * c + 1, c] -= L.FD_H
                ev = src.batch(*L.unpack(Zp, n))
                if k == "AMBr50" or k == "AMBr1":
                    pass
                g = (ev[0::2] - ev[1::2]) / (2 * L.FD_H)
                gs.append(float(np.linalg.norm(g) / iqr) if np.isfinite(g).all() else np.nan)
            Zl = Z[0][None] + ts[:, None] * u[None, :]
            el = src.batch(*L.unpack(Zl, n))
            rec["obj"][k].update({"gnorm_std": [float(x) for x in gs],
                                  "gnorm_std_mean": float(np.nanmean(gs)),
                                  "line": _panel(el, med, iqr)})
        rec["n_collapsed"]["AMBr1"] = int(r1.n_collapsed)
        rec["n_collapsed"]["AMBr50"] = int(r50.n_collapsed)
        rec["n_calls"] = {"AMB": int(sp.n_calls), "AMBr1": int(r1.n_calls),
                          "AMBr50": int(r50.n_calls)}
        out[pdb] = rec
        R.ck_save(tag, out)
        R.log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={n} {time.time()-t0:.0f}s  "
              f"log10range AMB {rec['obj']['AMB']['log10_range']:.2f} "
              f"r1 {rec['obj']['AMBr1']['log10_range']:.2f} "
              f"r50 {rec['obj']['AMBr50']['log10_range']:.2f} "
              f"LEG {rec['obj']['LEG']['log10_range']:.2f} | "
              f"rho(r50,sp) {rec['rank_agreement']['AMBr50|AMB']:+.3f} "
              f"rho(r50,r1) {rec['rank_agreement']['AMBr50|AMBr1']:+.3f} | "
              f"collapsed r50 {rec['n_collapsed']['AMBr50']}/{rec['n_calls']['AMBr50']}")
    with open(os.path.join(L.RESULTS, "qb2_relax_COMPLETE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return out


ARMS = ["LEG", "AMB", "AMBc", "AMBr1", "AMBr50"]


def report():
    d = R.ck_load("relax")
    if not d:
        print("no qb2_relax.json"); return {}
    pdbs = sorted(d)
    print(f"\n=== THE RELAXATION ARM: H_AMBER = E o Relax, decomposed.  n = {len(pdbs)} "
          f"targets, reference ensemble {d[pdbs[0]]['n_ref']} pool windows ===")
    out = {"n": len(pdbs), "targets": pdbs, "obj": {}, "rank_agreement": {}}
    keys = ["log10_range", "sd_over_iqr", "skew_p99_over_p50", "gnorm_std_mean",
            "spearman_E_vs_rmsd_ORACLE", "argmin_rmsd_ORACLE", "top10pct_mean_rmsd_ORACLE"]
    print(f"{'arm':<8}" + "".join(f"{k[:16]:>18}" for k in keys))
    for a in ARMS:
        line = f"{a:<8}"
        for k in keys:
            v = np.array([d[p]["obj"][a].get(k, np.nan) for p in pdbs], float)
            ci = L.boot_mean_ci(v)
            out["obj"].setdefault(a, {})[k] = ci
            line += f"{ci['mean']:>18.4g}"
        print(line)
    print(f"\n{'arm':<8}{'n_localmin':>14}{'acorr_len':>14}{'tv/range':>14}{'range_std':>14}")
    for a in ("AMB", "AMBr1", "AMBr50"):
        row = []
        for k in ("n_localmin", "acorr_len_rad", "tv_over_range", "range_std"):
            v = np.array([d[p]["obj"][a]["line"][k] for p in pdbs], float)
            ci = L.boot_mean_ci(v)
            out["obj"].setdefault(a, {})[f"line_{k}"] = ci
            row.append(ci["mean"])
        print(f"{a:<8}" + "".join(f"{x:>14.4g}" for x in row))
    print("\n--- RANK AGREEMENT (F-D2: does the relaxation reorder, or only condition?) ---")
    for k in d[pdbs[0]]["rank_agreement"]:
        v = np.array([d[p]["rank_agreement"][k] for p in pdbs], float)
        ci = L.boot_mean_ci(v)
        out["rank_agreement"][k] = ci
        flag = "REORDERS" if ci["ci95"][1] < 0.95 else ("conditioner" if ci["ci95"][0] >= 0.95
                                                        else "spans 0.95")
        print(f"  Spearman({k.replace('|', ', ')}) = {ci['mean']:+.3f} "
              f"[{ci['ci95'][0]:+.3f}, {ci['ci95'][1]:+.3f}]   {flag}")
    print("\n--- the collapse sentinel, as a MEASURED count (a guard that never fires is not evidence) ---")
    for a in ("AMB", "AMBr1", "AMBr50"):
        nc = int(sum(d[p]["n_collapsed"][a] for p in pdbs))
        nn = int(sum(d[p]["n_calls"][a] for p in pdbs))
        out.setdefault("collapsed", {})[a] = {"n_collapsed": nc, "n_calls": nn}
        print(f"  {a:<8} {nc} / {nn} calls = {100.0*nc/max(1,nn):.2f}%"
              + ("   (the bare path never calls _is_collapsed -- no infinite region by construction)"
                 if a == "AMB" else ""))
    ve = {a: max(d[p]["verify"][a][0] for p in pdbs) for a in ("AMB", "AMBr1", "AMBr50")}
    vn = {a: int(sum(d[p]["verify"][a][1] for p in pdbs)) for a in ("AMB", "AMBr1", "AMBr50")}
    print("\nbit-exactness gate vs core.amber.refine_coords: "
          + "  ".join(f"{a} max rel {ve[a]:.2e} over {vn[a]} comparisons" for a in ve))
    out["verify"] = {a: [ve[a], vn[a]] for a in ve}
    L.write("qb2_report_relax", out, complete=True)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    n_t = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    if mode == "run":
        run([t["pdb"] for t in L.subset(L.N_SUBSET)[:n_t]])
    report()


if __name__ == "__main__":
    main()
