"""SPRINT 18 / MATH -- S1: HOW SENSITIVE IS EACH OBJECTIVE'S OPTIMUM TO THE DISTOGRAM?

The coordinator's `s18/COORD_FINDING.md` establishes two things at n = 126: the functional
form is SOUND (with perfect distances the identical objective reaches 1.152 A), and the
distogram's error is WORSE THAN NOISE OF THE SAME SIZE (shuffling the residual directions
buys -1.001 A).  That makes one hypothesis about degree-1 precise and testable WITHOUT any
RMSD:

    > if the harmful part of the distogram's error acts through the PAIRWISE structure the
    > objective induces, the additive projection -- which averages each residue's
    > contribution over the distribution of the others -- should have an optimum that MOVES
    > LESS when dhat is perturbed.

That is a statement about `d theta* / d dhat`, and it is exactly computable from the
construction.  Nothing here reads a native.

------------------------------------------------------------------------------ THE ALGEBRA

At a stationary point `theta*` of an objective `E(theta; dhat)` the implicit function theorem
gives

    d theta* / d dhat_p  =  - H^{-1} d^2 E / d theta d dhat_p ,   H = grad^2_theta E(theta*)

FULL objective.  `E = sum_p w_p (d_p - dhat_p)^2`, so `dE/dtheta = sum_p 2 w_p (d_p - dhat_p)
grad d_p` and `d^2E/dtheta ddhat_p = -2 w_p grad d_p`.  Hence

    J_full = 2 H^{-1} G^T diag(w) ,   G = (P, 2n) pair Jacobian = d d_p / d theta

`G` is `s15.align_lib.pair_jacobian`, the programme's own exact restraint Jacobian.

RESIDUE-ADDITIVE objective.  `E_res = E_0(dhat) + sum_r f_r(theta_r; dhat)` and, from the
moment tables,  `f_r(t) = sum_p w_p (m2_rp(t) - 2 dhat_p m1_rp(t) + dhat_p^2) - E_0`, so

    d^2 f_r / dt d dhat_p = -2 w_p grad_t m1_rp(t)

and `H` is BLOCK DIAGONAL with 2x2 blocks `grad^2_t f_r`.  Every derivative comes from the
DFT of the tabulated moments, so `J_res` is analytic.

FROM TORSIONS TO ANGSTROMS.  A displacement of `theta*` is only meaningful as a displacement
of the STRUCTURE, and Ca-RMSD is superposition-invariant, so both Jacobians are pushed
through `align_lib.sup_jacobian` (which projects out the six rigid modes and therefore also
kills the four inert torsions).  The reported statistic is

    S = || J_sup J ||_F / sqrt(n)        Angstroms of Ca RMS motion of the optimum
                                         per Angstrom of i.i.d. perturbation of dhat

reported both unweighted and with the perturbation scaled by the distogram's own `sd_p`,
which is the physically matched version.

THE ONE MODELLING CHOICE, STATED.  `H` has exact zero modes (rigid + inert torsions) and, for
the full objective, near-zero ones: a flat direction means the optimum genuinely is not
determined there and the first-order sensitivity genuinely is unbounded.  A pseudo-inverse
cutoff is therefore required and it is a CHOICE.  It is swept (1e-8, 1e-6, 1e-4 relative) and
every rung is reported, together with the Hessian condition number and the number of modes
dropped -- because the conditioning gap may itself be the whole mechanism.

CONTROLS.  (i) ZERO-INFORMATION -- the sensitivity of the SHUFFLED-distogram objective, the
same functional at the same start with `dhat` permuted across pairs.  (ii) MATCHED-RANDOM --
a finite perturbation of matched RMS applied to both arms, re-optimising each from the SAME
start, which also validates the linearisation.

RUN:  python -m s18.math_sens run          (analytic, all cached targets)
      python -m s18.math_sens finite        (finite-perturbation validation subset)
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s18 import math_lib as M
from s18 import math_anova as A
from s18 import math_iface as MI

from s12 import instrument as I                   # noqa: E402
from s14 import avgspace as AV                    # noqa: E402
from s15 import align_lib as AL                   # noqa: E402
from s15 import seed as SD                        # noqa: E402

CUTOFFS = (1e-8, 1e-6, 1e-4)
N_PERT = 8
CONFIG = {"module": "s18/math_sens.py", "cutoffs": list(CUTOFFS), "n_pert": N_PERT,
          "salt": M.SALT, "grid": A.GRID, "nsamp": A.NSAMP, "mu": A.MU_DEFAULT}


def _start(pdb, seq, fold):
    """The production start: the ideal-geometry projection of the shipped top-75 coordinate
    average.  Identical to `s17/refine.py`.  NATIVE-FREE."""
    W = np.asarray(AV.top75_windows(pdb)[0], float)
    P = I.pairwise_rmsd(W)
    avg, _ = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    return np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)


def _hessian_fd(fn, phi, psi, eps=1e-5):
    """Central-difference Hessian of an objective exposing `(f, grad)`.  `4n` gradients."""
    n = len(phi)
    H = np.zeros((2 * n, 2 * n))
    for c in range(2 * n):
        gs = []
        for s in (+1, -1):
            p, q = phi.copy(), psi.copy()
            (p if c < n else q)[c % n] += s * eps
            gs.append(fn(p, q)[1])
        H[:, c] = (gs[0] - gs[1]) / (2 * eps)
    return 0.5 * (H + H.T)


def _pinv_sens(H, B, Jsup, n, cutoffs=CUTOFFS, direction=None):
    """`|| Jsup H^+ B ||_F / sqrt(n)` at several relative pseudo-inverse cutoffs.

    `direction` (P,) optionally adds the DIRECTIONAL response `|| Jsup H^+ B u ||/sqrt(n)`
    along a named unit perturbation of `dhat` -- used only for the ORACLE-labelled
    true-correction direction, never for a native-free claim.
    """
    u, s, vt = np.linalg.svd(H)
    out = {}
    for c in cutoffs:
        keep = s > c * max(s.max(), 1e-300)
        si = np.where(keep, 1.0 / np.where(keep, s, 1.0), 0.0)
        Hp = (vt.T * si) @ u.T
        X = Jsup @ (Hp @ B)
        out[f"{c:.0e}"] = {"S": float(np.linalg.norm(X) / np.sqrt(n)),
                           "modes_kept": int(keep.sum())}
        if direction is not None:
            out[f"{c:.0e}"]["S_dir_ORACLE"] = float(
                np.linalg.norm(X @ direction) / np.sqrt(n))
    out["cond"] = float(s.max() / max(s[s > 1e-14 * s.max()].min(), 1e-300))
    out["eig_min_pos"] = float(s[s > 1e-14 * s.max()].min())
    out["eig_max"] = float(s.max())
    return out


def cell(pdb, target):
    t0 = time.time()
    o = MI.build(pdb, target=target)
    tt = o.t
    n = tt.n
    w, sd = tt.w, tt.sd
    phi0, psi0 = _start(pdb, tt.seq, tt.fold)
    out = {"pdb": pdb, "n": n, "fold": tt.fold, "n_pairs": int(len(tt.i))}
    #: ORACLE DIAGNOSTIC, post hoc only.  The unit vector pointing from the predicted
    #: distances to the true ones -- the coordinator's alpha direction.  Used ONLY for the
    #: `S_dir_ORACLE` rows; no native-free statistic depends on it.
    nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)
    dtrue = np.sqrt(((nat[tt.i] - nat[tt.j]) ** 2).sum(1))
    dirv = dtrue - tt.dhat
    dirv = dirv / max(float(np.linalg.norm(dirv)), 1e-12)
    out["resid_rms_ORACLE"] = float(np.sqrt(((dtrue - tt.dhat) ** 2).mean()))

    # ---------------------------------------------------------------- FULL objective
    from s15 import align_lib as _AL
    pf, qf, _ = _AL.fit(tt.dhat, sd, tt.i, tt.j, phi0, psi0)
    H = _hessian_fd(o.E_full, pf, qf)
    Jr, CA = AL.raw_jacobian(pf, qf)
    G, _d = AL.pair_jacobian(Jr, CA, tt.i, tt.j)              # (P, 2n)
    Jsup, _, _ = AL.sup_jacobian(pf, qf)
    B = 2.0 * (G.T * w[None, :])                              # (2n, P)  = -d^2E/dth ddhat
    out["full"] = _pinv_sens(H, B, Jsup, n, direction=dirv)
    out["full_sdscaled"] = _pinv_sens(H, B * sd[None, :], Jsup, n)
    out["full_objective_at_opt"] = float(o.E_full(pf, qf)[0])

    # ---------------------------------------------------------------- RESIDUE-ADDITIVE
    pr_, qr_ = o.argmin_res(phi0, psi0)
    pr_, qr_ = _descend_res(o, pr_, qr_)
    Hr = _hessian_fd(o.E_res, pr_, qr_)
    Jsup_r, _, _ = AL.sup_jacobian(pr_, qr_)
    Br = _dEres_ddhat(tt, pr_, qr_)                           # (2n, P)
    out["res"] = _pinv_sens(Hr, Br, Jsup_r, n, direction=dirv)
    out["res_sdscaled"] = _pinv_sens(Hr, Br * sd[None, :], Jsup_r, n)
    out["res_objective_at_opt"] = float(o.E_res(pr_, qr_)[0])

    # ---------------------------------------------------------------- ANGLE-ADDITIVE
    pa, qa = o.argmin_ang(phi0, psi0)
    pa, qa = _descend(o.E_ang, pa, qa)
    Ha = _hessian_fd(o.E_ang, pa, qa)
    Jsup_a, _, _ = AL.sup_jacobian(pa, qa)
    Ba = _dEang_ddhat(tt, pa, qa)
    out["ang"] = _pinv_sens(Ha, Ba, Jsup_a, n, direction=dirv)
    out["ang_sdscaled"] = _pinv_sens(Ha, Ba * sd[None, :], Jsup_a, n)

    # ---------------------------------------------------------------- CONTROL: shuffled
    rng = SD.stable_rng(pdb, "shuffle", salt=M.SALT)
    perm = rng.permutation(len(tt.dhat))
    ps, qs, _ = _AL.fit(tt.dhat[perm], sd[perm], tt.i, tt.j, phi0, psi0)
    osh = MI.Obj(tt)
    osh.t = tt
    Hs = _hessian_fd(lambda p, q: _full_g(tt, p, q, tt.dhat[perm], 1.0 / sd[perm] ** 2),
                     ps, qs)
    Jrs, CAs = AL.raw_jacobian(ps, qs)
    Gs, _ = AL.pair_jacobian(Jrs, CAs, tt.i, tt.j)
    Jss, _, _ = AL.sup_jacobian(ps, qs)
    Bs = 2.0 * (Gs.T * (1.0 / sd[perm] ** 2)[None, :])
    out["CONTROL_shuffled_full"] = _pinv_sens(Hs, Bs, Jss, n)

    out["seconds"] = time.time() - t0
    return out


def _full_g(tt, phi, psi, dhat, w):
    from core import project as pj
    phi = np.asarray(phi, float).ravel()
    psi = np.asarray(psi, float).ravel()
    Gm = pj.frames(phi[None], psi[None])[0]
    CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
    rv = CA[tt.i] - CA[tt.j]
    d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
    r = d - dhat
    f = float((w * r * r).sum())
    coef = ((2.0 * w * r) / d)[:, None] * rv
    gCA = np.zeros_like(CA)
    np.add.at(gCA, tt.i, coef)
    np.add.at(gCA, tt.j, -coef)
    return f, np.asarray(pj._torsion_grad(Gm, CA, gCA), float)


def _descend(fn, phi, psi, maxiter=200):
    from scipy.optimize import minimize
    n = len(phi)

    def fg(x):
        f, g = fn(x[:n], x[n:])
        return f, g

    r = minimize(fg, np.concatenate([phi, psi]), jac=True, method="L-BFGS-B",
                 options={"maxiter": maxiter, "ftol": 1e-14, "gtol": 1e-12})
    return r.x[:n], r.x[n:]


def _descend_res(o, phi, psi, maxiter=200):
    return _descend(o.E_res, phi, psi, maxiter)


def _grad_tab2(F, phi, psi):
    """`d/dphi` and `d/dpsi` of a `(n, G, G)` DFT at one torsion vector."""
    G = F.shape[1]                      # (n, G, G, P) -- the mesh axes are 1 and 2
    k = np.fft.fftfreq(G, 1.0 / G)
    u = (np.asarray(phi, float) + np.pi) / (2 * np.pi) * G - 0.5
    v = (np.asarray(psi, float) + np.pi) / (2 * np.pi) * G - 0.5
    Ep = np.exp(2j * np.pi * u[:, None] * k[None, :] / G)
    Eq = np.exp(2j * np.pi * v[:, None] * k[None, :] / G)
    dp = np.einsum("nu,nv,nuvp->np", Ep * (1j * k[None, :]), Eq, F).real / (G * G)
    dq = np.einsum("nu,nv,nuvp->np", Ep, Eq * (1j * k[None, :]), F).real / (G * G)
    return dp, dq


def _dEres_ddhat(tt, phi, psi):
    """`(2n, P)` = `-d^2 E_res / d theta d dhat` = `2 w_p grad_theta m1_rp`, analytic."""
    n, G = tt.n, tt.grid
    m1 = tt.m1.reshape(n, G, G, -1)
    F = np.fft.fft2(m1, axes=(1, 2))
    dp, dq = _grad_tab2(F, phi, psi)                # (n, P)
    out = np.zeros((2 * n, m1.shape[-1]))
    out[:n] = 2.0 * dp * tt.w[None, :]
    out[n:] = 2.0 * dq * tt.w[None, :]
    out[0] = out[n - 1] = out[n] = out[2 * n - 1] = 0.0
    return out


def _dEang_ddhat(tt, phi, psi):
    n, G = tt.n, tt.grid
    k = np.fft.fftfreq(G, 1.0 / G)
    Fa = np.fft.fft(tt.am1, axis=1)
    Fb = np.fft.fft(tt.bm1, axis=1)
    u = (phi + np.pi) / (2 * np.pi) * G - 0.5
    v = (psi + np.pi) / (2 * np.pi) * G - 0.5
    Ep = np.exp(2j * np.pi * u[:, None] * k[None, :] / G) * (1j * k[None, :])
    Eq = np.exp(2j * np.pi * v[:, None] * k[None, :] / G) * (1j * k[None, :])
    da = np.einsum("nu,nup->np", Ep, Fa).real / G
    db = np.einsum("nu,nup->np", Eq, Fb).real / G
    out = np.zeros((2 * n, tt.am1.shape[-1]))
    out[:n] = 2.0 * da * tt.w[None, :]
    out[n:] = 2.0 * db * tt.w[None, :]
    out[0] = out[n - 1] = out[n] = out[2 * n - 1] = 0.0
    return out


# ------------------------------------------------------------------ finite validation
def finite_cell(pdb, target, n_pert=N_PERT, sigma=0.5):
    """MATCHED-RANDOM control and linearisation check: perturb `dhat` by i.i.d. Gaussian
    noise of the SAME magnitude for both arms, re-optimise each from the SAME start, and
    measure how far the optimum MOVES (Ca-RMSD to its own unperturbed optimum).  NATIVE-FREE.
    """
    from s15 import align_lib as _AL
    o = MI.build(pdb, target=target)
    tt = o.t
    phi0, psi0 = _start(pdb, tt.seq, tt.fold)
    pf0, qf0, _ = _AL.fit(tt.dhat, tt.sd, tt.i, tt.j, phi0, psi0)
    pr0, qr0 = _descend_res(o, *o.argmin_res(phi0, psi0))
    CAf0 = np.asarray(I.build_ca(pf0, qf0), float)
    CAr0 = np.asarray(I.build_ca(pr0, qr0), float)
    df, dr = [], []
    for t in range(n_pert):
        rng = SD.stable_rng(pdb, t, "pert", salt=M.SALT)
        dd = tt.dhat + sigma * rng.standard_normal(len(tt.dhat))
        pf, qf, _ = _AL.fit(dd, tt.sd, tt.i, tt.j, phi0, psi0)
        df.append(float(I.ca_rmsd(np.asarray(I.build_ca(pf, qf), float), CAf0)))
        o2 = o.with_dhat(dd)
        pr, qr = _descend_res(o2, *o2.argmin_res(phi0, psi0))
        dr.append(float(I.ca_rmsd(np.asarray(I.build_ca(pr, qr), float), CAr0)))
    return {"pdb": pdb, "n": tt.n, "fold": tt.fold, "sigma": sigma,
            "full_argmin_move": df, "res_argmin_move": dr,
            "full_mean": float(np.mean(df)), "res_mean": float(np.mean(dr))}


# ------------------------------------------------------------------------- driver
def run(mode="run", targets=None):
    tg = targets if targets is not None else I.targets()
    tag = "sens" if mode == "run" else "sensfin"
    done = M.ck_load(tag)
    t0 = time.time()
    ok = 0
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        if f"cell_{pdb}" in done:
            ok += 1
            continue
        path = os.path.join(M.CACHE,
                            f"anova_{pdb}_{A.MU_DEFAULT}_{A.NSAMP}_{A.GRID}.npz")
        if not os.path.exists(path):
            continue                       # coefficients not built yet; skip, rerun later
        try:
            r = cell(pdb, t) if mode == "run" else finite_cell(pdb, t)
        except Exception as e:                                   # pragma: no cover
            print(f"  {pdb} FAILED {type(e).__name__}: {e}", flush=True)
            continue
        M.ck(tag, f"cell_{pdb}", r)
        ok += 1
        if mode == "run":
            print(f"  {pdb} n{r['n']}  S_full {r['full_sdscaled']['1e-06']['S']:.3f}  "
                  f"S_res {r['res_sdscaled']['1e-06']['S']:.3f}  "
                  f"S_ang {r['ang_sdscaled']['1e-06']['S']:.3f}  "
                  f"cond {r['full']['cond']:.1e}/{r['res']['cond']:.1e}  "
                  f"({time.time()-t0:.0f}s)", flush=True)
        else:
            print(f"  {pdb} full {r['full_mean']:.3f}  res {r['res_mean']:.3f}",
                  flush=True)
    print(f"  {ok} cells", flush=True)
    if ok == len(tg):
        M.seal(tag, CONFIG)
    return ok


def report(tag="sens"):
    """Aggregate `math_sens.json`.  TARGET is the unit; median and W/L beside every mean."""
    import json
    d = M.ck_load(tag)
    ks = sorted(k for k in d if k.startswith("cell_"))
    if not ks:
        raise SystemExit("no cells")
    complete = bool(d.get("_complete"))
    folds = np.asarray([d[k]["fold"] for k in ks], int)
    print(f"n = {len(ks)} targets.  artefact complete = {complete}"
          + ("" if complete else "   <-- PARTIAL: SMOKE ONLY, no directional claim"))
    print("\nS = Angstroms of Ca RMS motion of the OPTIMUM per Angstrom of perturbation of")
    print("dhat, superposition-invariant.  Native-free.  Lower = the optimum is more robust.\n")
    for scale in ("", "_sdscaled"):
        for cut in ("1e-08", "1e-06", "1e-04"):
            g = lambda a: np.array([d[k][a + scale][cut]["S"] for k in ks])  # noqa: E731
            fu, re_, an = g("full"), g("res"), g("ang")
            sh = np.array([d[k]["CONTROL_shuffled_full"][cut]["S"] for k in ks])
            s1 = M.paired(re_, fu, folds)
            s2 = M.paired(an, fu, folds)
            lab = "sd-scaled" if scale else "unweighted"
            print(f"  --- {lab}, pinv cutoff {cut} ---")
            print(f"  {'arm':<28}{'S':>8}{'median':>9}"
                  f"{'vs full [95% CI]':>28}{'W/L':>9}{'fold sign':>11}")
            print(f"  {'full objective':<28}{fu.mean():>8.3f}{np.median(fu):>9.3f}"
                  f"{'--':>28}{'--':>9}{'--':>11}")
            for nm, v, s in (("residue-additive", re_, s1), ("angle-additive", an, s2)):
                print(f"  {nm:<28}{v.mean():>8.3f}{np.median(v):>9.3f}"
                      f"   {s['mean_diff']:+.3f} [{s['ci95'][0]:+.3f},{s['ci95'][1]:+.3f}]"
                      f"{s['n_better']:>4}/{s['n_worse']}{s['folds_same_sign']:>11}")
            print(f"  {'CONTROL shuffled dhat (full)':<28}{sh.mean():>8.3f}"
                  f"{np.median(sh):>9.3f}")
            print()
    cf = np.array([d[k]["full"]["cond"] for k in ks])
    cr = np.array([d[k]["res"]["cond"] for k in ks])
    ef = np.array([d[k]["full"]["eig_min_pos"] for k in ks])
    er = np.array([d[k]["res"]["eig_min_pos"] for k in ks])
    print("  Hessian conditioning at the optimum (median over targets)")
    print(f"    full     cond {np.median(cf):.2e}   smallest positive eigenvalue "
          f"{np.median(ef):.3g}   largest {np.median([d[k]['full']['eig_max'] for k in ks]):.3g}")
    print(f"    residue  cond {np.median(cr):.2e}   smallest positive eigenvalue "
          f"{np.median(er):.3g}   largest {np.median([d[k]['res']['eig_max'] for k in ks]):.3g}")
    od = np.array([d[k]["full"]["1e-06"]["S_dir_ORACLE"] for k in ks])
    orr = np.array([d[k]["res"]["1e-06"]["S_dir_ORACLE"] for k in ks])
    print("\n  ORACLE DIAGNOSTIC -- motion of the optimum along the TRUE-correction direction")
    print(f"    full {od.mean():.3f}   residue-additive {orr.mean():.3f}   "
          f"(residual RMS {np.mean([d[k]['resid_rms_ORACLE'] for k in ks]):.3f} A)")
    out = {"n": len(ks), "complete": complete}
    M.ck("sensreport", "aggregate", out)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "report":
        report()
    else:
        run(mode)
