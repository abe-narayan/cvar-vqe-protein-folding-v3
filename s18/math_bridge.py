"""SPRINT 18 / MATH -- B1: the continuous machinery checked against the EXACT enumerated ANOVA.

`math_lattice.py` proves the ALGEBRA (order-1 ANOVA == Walsh weight-<=1, to 1e-15).  It says
nothing about whether the CONTINUOUS implementation -- frozen Monte-Carlo conditional
expectations on a 16x16 mesh, carried by trigonometric interpolation -- actually computes that
same object.  This module closes that gap, on the objective the PRODUCTION pipeline uses.

THE CONSTRUCTION.  Take one of the 19 exhaustively enumerated targets.  Its k = 4 torsion
library gives `k^n` configurations.  Evaluate the PRODUCTION objective
`sum_p ((d_p - dhat_p)/sd_p)^2` -- not the enumerated `hamil`, which is a different objective
-- on every one of them, and take the EXACT residue-additive ANOVA under the uniform measure
on those `k` atoms per residue.  Then build the continuous object with `mu = "lattice"`, the
SAME reference measure, and compare.

WHAT THE COMPARISON MEASURES.  Both error sources of the deliverable at once: the frozen
Monte-Carlo quadrature (S draws) and the mesh + trigonometric interpolation.  Reported as the
worst and RMS discrepancy over all lattice configurations, relative to `sd(E_le1_exact)` --
the scale on which the object is used.

A BONUS THAT IS WORTH MORE THAN THE CHECK.  Because the enumeration is exhaustive, this also
yields the CERTIFIED argmin of the production objective and of its residue-additive
truncation, on the same 19 targets Phase 0 used, with **no search error and no Monte Carlo**.
That is the closest exact analogue of the 2.411/2.661 comparison that can be run on the
objective the pipeline actually deploys.  Scored ORACLE, post hoc.

RUN:  python -m s18.math_bridge run
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s18 import math_lib as M
from s18 import math_anova as A
from s18 import math_iface as MI

from s12 import instrument as I                    # noqa: E402
from s16 import qphase_lib as QP                   # noqa: E402

CONFIG = {"module": "s18/math_bridge.py", "mu": "lattice", "grid": A.GRID,
          "nsamp": A.NSAMP, "salt": M.SALT}


def _field_at(tc, phi, psi):
    """`(n,)` -- each residue's continuous field `f_r` evaluated at its own `(phi_r, psi_r)`,
    by the same exact trigonometric interpolation the deliverable uses."""
    G = tc.grid
    kk = np.fft.fftfreq(G, 1.0 / G)
    u = (np.asarray(phi, float) + np.pi) / (2 * np.pi) * G - 0.5
    v = (np.asarray(psi, float) + np.pi) / (2 * np.pi) * G - 0.5
    Ep = np.exp(2j * np.pi * u[:, None] * kk[None, :] / G)
    Eq = np.exp(2j * np.pi * v[:, None] * kk[None, :] / G)
    return np.einsum("nu,nv,nuv->n", Ep, Eq, tc.F).real / (G * G)


def enumerate_objective(tt, ins, chunk=16384):
    """The production objective on every lattice configuration.  EXACT, no sampling."""
    n, k, N = ins.n, ins.k, ins.N
    S = M.states_of(N, n, k)
    ar = np.arange(n)
    out = np.empty(N, float)
    for a in range(0, N, chunk):
        b = min(a + chunk, N)
        CA = np.asarray(I.build_ca(ins.PHI[ar[None, :], S[a:b]],
                                   ins.PSI[ar[None, :], S[a:b]]), float)
        d = np.sqrt(((CA[:, tt.i] - CA[:, tt.j]) ** 2).sum(-1))
        r = d - tt.dhat[None, :]
        out[a:b] = (r * r * tt.w[None, :]).sum(1)
    return out, S


def cell(pdb, target):
    t0 = time.time()
    ins = QP.inst(pdb)
    n, k = ins.n, ins.k
    d = A.gather_one(target)
    assert d["n"] == n, (d["n"], n)
    tt = A.Target(pdb, d["seq"], n, d["fold"], d["dhat"], d["sd"], d["i"], d["j"])
    E, S = enumerate_objective(tt, ins)

    # ---- EXACT residue-additive ANOVA of the production objective on the lattice
    E0x, fx, Ex_le1 = M.anova1_discrete(E, n, k)
    rmsd = ins.rmsd                                       # ORACLE, post hoc only

    # ---- the CONTINUOUS object under the SAME reference measure
    path = os.path.join(M.CACHE, f"anova_{pdb}_lattice_{A.NSAMP}_{A.GRID}.npz")
    if os.path.exists(path):
        try:
            tc = A.Target.load(path, pdb)
        except IOError:
            os.remove(path)
            tc = None
    else:
        tc = None
    if tc is None:
        tc = A.Target(pdb, d["seq"], n, d["fold"], d["dhat"], d["sd"],
                      d["i"], d["j"]).fit(mu="lattice", S=A.NSAMP, grid=A.GRID)
        tc.save(path)
    o = MI.Obj(tc)

    # ---- the continuous field evaluated AT the k lattice atoms of each residue.
    #      Because both objects are additive, one (n, k) table reconstructs the continuous
    #      E_res on every one of the k^n lattice configurations exactly -- no loop.
    ar = np.arange(n)
    fc = np.empty((n, k))
    for s in range(k):
        _, _ = 0, 0
        fc[:, s] = _field_at(tc, ins.PHI[:, s], ins.PSI[:, s])
    gall = tc.E0 + fc[ar[None, :], S].sum(1)
    want = Ex_le1
    err = gall - want
    scale = float(Ex_le1.std() + 1e-300)
    got = gall

    tied_full = np.flatnonzero(E == E.min())
    tied_res = np.flatnonzero(Ex_le1 == Ex_le1.min())
    tied_c = np.flatnonzero(gall == gall.min())

    out = {
        "pdb": pdb, "n": n, "k": k, "N": int(ins.N), "fold": int(ins.fold),
        "bridge_max_abs_err": float(np.abs(err).max()),
        "bridge_rms_err": float(np.sqrt((err ** 2).mean())),
        "bridge_max_rel_err": float(np.abs(err).max() / scale),
        "bridge_rms_rel_err": float(np.sqrt((err ** 2).mean()) / scale),
        "corr_continuous_vs_exact": float(np.corrcoef(got, want)[0, 1]),
        "spearman_continuous_vs_exact": M.spearman(got, want),
        "sd_E_le1_exact": scale,
        "sd_E_full": float(E.std()),
        # ---- the certified comparison on the PRODUCTION objective (ORACLE scoring)
        "argmin_rmsd_full_ORACLE": float(rmsd[tied_full].mean()),
        "argmin_rmsd_res_exact_ORACLE": float(rmsd[tied_res].mean()),
        "argmin_rmsd_res_continuous_ORACLE": float(rmsd[tied_c].mean()),
        "ties_full": int(tied_full.size), "ties_res": int(tied_res.size),
        "rmsd_best_in_space_ORACLE": float(rmsd.min()),
        "rmsd_mean_in_space_ORACLE": float(rmsd.mean()),
        "rho_full_rmsd_ORACLE": M.spearman(E, rmsd),
        "rho_res_rmsd_ORACLE": M.spearman(Ex_le1, rmsd),
        "var_frac_residue_le1": float(Ex_le1.var() / max(E.var(), 1e-300)),
        "E0_exact": float(E0x), "E0_continuous": float(tc.E0),
        "seconds": time.time() - t0,
    }
    return out


def report():
    d = M.ck_load("bridge")
    ks = sorted(k for k in d if k.startswith("cell_"))
    if not ks:
        raise SystemExit("no cells")
    folds = np.asarray([d[k]["fold"] for k in ks], int)
    g = lambda a: np.array([d[k][a] for k in ks])          # noqa: E731
    print(f"n = {len(ks)} enumerated targets.  complete = {bool(d.get('_complete'))}\n")
    print("B1  DOES THE CONTINUOUS IMPLEMENTATION COMPUTE THE EXACT ANOVA?")
    print(f"  worst  |E_res_continuous - E_res_exact| / sd(E_res_exact) : "
          f"{g('bridge_max_rel_err').max():.4f}")
    print(f"  mean RMS relative error                                   : "
          f"{g('bridge_rms_rel_err').mean():.4f}")
    print(f"  mean Spearman(continuous, exact) over lattice configs     : "
          f"{g('spearman_continuous_vs_exact').mean():.4f}")
    print(f"  mean Pearson                                              : "
          f"{g('corr_continuous_vs_exact').mean():.4f}")
    print()
    print("B2  THE CERTIFIED COMPARISON ON THE PRODUCTION OBJECTIVE (ORACLE scoring)")
    full = g("argmin_rmsd_full_ORACLE")
    res = g("argmin_rmsd_res_exact_ORACLE")
    con = g("argmin_rmsd_res_continuous_ORACLE")
    print(f"  {'object':<34}{'argmin RMSD':>12}{'median':>9}"
          f"{'vs full [95% CI]':>28}{'W/L/T':>10}{'fold sign':>11}")
    print(f"  {'full production objective':<34}{full.mean():>12.3f}{np.median(full):>9.3f}"
          f"{'--':>28}{'--':>10}{'--':>11}")
    for nm, v in (("residue-additive, EXACT", res),
                  ("residue-additive, the deliverable", con)):
        s = M.paired(v, full, folds)
        print(f"  {nm:<34}{v.mean():>12.3f}{np.median(v):>9.3f}"
              f"   {s['mean_diff']:+.3f} [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]"
              f"{s['n_better']:>4}/{s['n_worse']}/{s['n_tie']}{s['folds_same_sign']:>11}")
    print(f"  {'zero-information (uniform draw)':<34}"
          f"{g('rmsd_mean_in_space_ORACLE').mean():>12.3f}")
    print(f"  {'ORACLE ceiling (space best)':<34}"
          f"{g('rmsd_best_in_space_ORACLE').mean():>12.3f}")
    print(f"\n  retained variance of the residue-additive projection: "
          f"{g('var_frac_residue_le1').mean():.3f}")
    print(f"  rho(objective, RMSD): full {g('rho_full_rmsd_ORACLE').mean():+.3f}   "
          f"residue-additive {g('rho_res_rmsd_ORACLE').mean():+.3f}")
    M.ck("bridgereport", "aggregate", {"n": len(ks)})
    return d


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "report":
        report()
    else:
        tg = {t["pdb"]: t for t in I.targets()}
        done = M.ck_load("bridge")
        names = sys.argv[2].split(",") if len(sys.argv) > 2 else list(QP.TARGETS19)
        for pdb in names:
            if f"cell_{pdb}" in done:
                print("  skip", pdb, flush=True)
                continue
            if pdb not in tg:
                print(f"  {pdb} not in the 126-target instrument -- skipped", flush=True)
                continue
            r = cell(pdb, tg[pdb])
            M.ck("bridge", f"cell_{pdb}", r)
            print(f"  {pdb} {r['seconds']:.0f}s  bridge rel err "
                  f"max {r['bridge_max_rel_err']:.4f} rms {r['bridge_rms_rel_err']:.4f}  "
                  f"argmin full {r['argmin_rmsd_full_ORACLE']:.3f} "
                  f"res {r['argmin_rmsd_res_exact_ORACLE']:.3f} "
                  f"cont {r['argmin_rmsd_res_continuous_ORACLE']:.3f}", flush=True)
        if all(f"cell_{p}" in M.ck_load("bridge") for p in QP.TARGETS19 if p in tg):
            print("  sealing:", M.seal("bridge", CONFIG))
