"""SPRINT 16, AUDIT -- adversarial re-test of the two Sprint 15 measurements the flagship
rests on, plus the ORACLE ceiling's construction and the quiet-subspace theory check.

Everything native-derived is named ORACLE_*.  The statistical unit is the TARGET.

PART A  the error-direction surrogate  `theta_fit - theta_pool`, recorded |cos| 0.390 vs a
        0.157 null and alignment-tracking +0.499 [+0.33,+0.65].
        Audited against: (a) MATCHED nulls that preserve the surrogate's own per-coordinate
        magnitude profile and per-target dimensionality; (b) a ZERO-INFORMATION reference
        direction (the constant alpha-helix) in place of the pool mean; (c) restriction to the
        coordinates the builder actually reads and to the RMSD-relevant complement of the exact
        null space; (d) a fold-clustered bootstrap; (e) partialling out |e|.

PART B  the quiet-subspace estimate, recorded cos^2 0.827 vs a 0.499 null.
        Audited against: the analytic k/m value of that null, the shrinking-subspace ladder
        (half / quarter / eighth / the exactly-null 5), a ZERO-INFORMATION Jacobian control
        (constant alpha-helix) and a random-pool-window control, and the deflated overlap with
        the four known-inert coordinate axes removed from both subspaces.

PART C  what `ORACLE_kill_loud_norm_0.5` actually constructs, and the control that says whether
        it carries error-DIRECTION information at all:
            ORACLE_quiet_random_norm  =  theta_native + |e| * (random unit vector in the
            bottom-half right-singular subspace of J at the fit).
        Same magnitude, same subspace, ZERO information about which way the fit is wrong.

PART D  the theory check.  Curvature of the restraint OBJECTIVE and of RMSD along quiet vs loud
        torsion directions, plus a finite-step version at a realistic 20 deg RMS.

    python -m s16.audit_align            # all parts, 126 targets
    python -m s16.audit_align 10         # smoke on 10
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ["" + _v] = "2"

import numpy.linalg as la                     # noqa: E402
from s12 import instrument as I               # noqa: E402
from s15 import align_lib as A                # noqa: E402
from s15 import distcal as C                  # noqa: E402
from s15 import distgeo as D                  # noqa: E402
from s15 import info_lib as L                 # noqa: E402
from s15 import seed as SD                    # noqa: E402
from s15.info_regime import ALPHA             # noqa: E402
from s14 import retprior as RP                # noqa: E402
import peptide_db as pdb                      # noqa: E402

OUT = os.path.join(ROOT, "s16", "results")
os.makedirs(OUT, exist_ok=True)
PATH = os.path.join(OUT, "audit_align.json")


# ------------------------------------------------------------------ statistics helpers
def cluster_boot_corr(x, y, folds, n_boot=4000, seed=0):
    """95% CI for pearson(x,y) resampling WHOLE FOLDS (the clustered unit), not targets."""
    x = np.asarray(x, float); y = np.asarray(y, float); folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = {f: np.where(folds == f)[0] for f in uf}
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(uf), len(uf))
        k = np.concatenate([idx[uf[q]] for q in pick])
        v.append(L.pearson(x[k], y[k]))
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def cluster_boot_mean(d, folds, n_boot=4000, seed=0):
    d = np.asarray(d, float); folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = {f: np.where(folds == f)[0] for f in uf}
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(uf), len(uf))
        k = np.concatenate([idx[uf[q]] for q in pick])
        v.append(d[k].mean())
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def partial_corr(x, y, *covs):
    return L.pearson(L.residualise(x, *covs), L.residualise(y, *covs))


def principal_angles(A1, A2):
    return la.svd(A1.T @ A2, compute_uv=False)


# ------------------------------------------------------------------ the per-target work
def per_target(t, data, deb, n_start=6, n_null=200):
    p = t["pdb"]
    d = data[p]
    n = d["n"]
    m = 2 * n
    i, j, sd = d["i"], d["j"], d["sd"]
    dhat = np.maximum(d["dhat"] - deb[d["fold"]](d["sep"]), 2.0)
    nat = d["nat"]
    nt = pdb.by_pdb(p)
    nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
    w0 = 1.0 / (sd * sd)
    rng = SD.stable_rng(p, "s16audit")

    # ---- the identical `squared` fit that align_free / align_fit use
    S = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p))
    sols = []
    for phi0, psi0, _tag in S:
        ph, ps, fv = A.fit(dhat, sd, i, j, phi0, psi0, wpair=w0)
        sols.append((fv, ph, ps, np.concatenate([np.asarray(phi0, float),
                                                 np.asarray(psi0, float)])))
    fv, phi, psi, th0 = min(sols, key=lambda z: z[0])
    th = np.concatenate([phi, psi])
    thn = np.concatenate([nphi, npsi])

    Jn, _, CAn = A.sup_jacobian(nphi, npsi)               # ORACLE J
    Jf, Jraw_f, CAf = A.sup_jacobian(phi, psi)            # NATIVE-FREE J at the fit
    e_true = A.wrap(th - thn)                             # ORACLE error vector

    # ---- the surrogates
    P7, S7, _ = RP.windows(p, "top75")
    cm = np.concatenate([RP.circ_mean(P7, None, 0), RP.circ_mean(S7, None, 0)])
    th_hel = np.concatenate([np.full(n, ALPHA[0]), np.full(n, ALPHA[1])])
    kref = int(SD.stable_rng(p, "ref").integers(0, P7.shape[0]))
    th_rw = np.concatenate([P7[kref], S7[kref]])
    sur = {"d_fit_pool": A.wrap(th - cm),
           "d_fit_start": A.wrap(th - th0),
           "d_fit_helix_ZEROINFO": A.wrap(th - th_hel),
           "d_fit_randwindow_ZEROINFO": A.wrap(th - th_rw)}
    #: unwrapped twin of the recorded surrogate, to price the wrapping choice
    sur_unwrapped = th - cm
    e_unwrapped = th - thn

    r = {"pdb": p, "n": n, "fold": int(d["fold"]),
         "rmsd_squared": float(I.ca_rmsd(I.build_ca(phi, psi), nat)),
         "tors_rms_deg": float(np.degrees(np.sqrt((e_true ** 2).mean()))),
         "norm_e_ORACLE": float(la.norm(e_true)),
         "align_true_ORACLE": A.alignment(Jn, e_true),
         "align_true_selfJ": A.alignment(Jf, e_true)}

    ne = float(la.norm(e_true))

    # ---- masks.  The builder never reads phi[0], phi[n-1], psi[n-1]; the fit therefore
    #      leaves those three coordinates exactly at their START value, so BOTH the surrogate
    #      and the true error contain the same +theta_start there.
    live = np.ones(m, bool)
    live[[0, n - 1, m - 1]] = False                       # phi0, phi_{n-1}, psi_{n-1}
    #: RMSD-relevant complement of the EXACT null space of J at the fit (5 directions)
    s_f, V_f = A.spectrum(Jf)
    nullmask = s_f <= 1e-8 * s_f[0]
    Vnull = V_f[:, nullmask]
    r["n_exact_null_fit"] = int(nullmask.sum())

    def cosabs(u, v):
        nu, nv = la.norm(u), la.norm(v)
        return float(abs(u @ v) / (nu * nv)) if nu > 1e-12 and nv > 1e-12 else np.nan

    def cossgn(u, v):
        nu, nv = la.norm(u), la.norm(v)
        return float((u @ v) / (nu * nv)) if nu > 1e-12 and nv > 1e-12 else np.nan

    for k, v in sur.items():
        r["cos_" + k] = cosabs(v, e_true)
        r["cossgn_" + k] = cossgn(v, e_true)
        r["align_" + k] = A.alignment(Jf, v)
        # live-coordinate restriction
        r["coslive_" + k] = cosabs(v[live], e_true[live])
        # RMSD-relevant complement of the exact null space
        vq = v - Vnull @ (Vnull.T @ v)
        eq = e_true - Vnull @ (Vnull.T @ e_true)
        r["cosnonnull_" + k] = cosabs(vq, eq)
    r["cos_unwrapped_d_fit_pool"] = cosabs(sur_unwrapped, e_unwrapped)
    r["cos_wrapped_sur_unwrapped_e"] = cosabs(sur["d_fit_pool"], e_unwrapped)

    # ---- MATCHED NULLS for the headline surrogate
    v = sur["d_fit_pool"]
    nv = la.norm(v)
    c_iso, c_flip, c_perm = [], [], []
    for _ in range(n_null):
        g = rng.normal(size=m)
        c_iso.append(cosabs(g, e_true))
        c_flip.append(cosabs(v * rng.choice([-1.0, 1.0], m), e_true))
        q = v.copy()
        q[:n] = rng.permutation(q[:n]); q[n:] = rng.permutation(q[n:])
        c_perm.append(cosabs(q, e_true))
    r["null_iso"] = float(np.mean(c_iso))
    r["null_signflip_MATCHED"] = float(np.mean(c_flip))
    r["null_perm_MATCHED"] = float(np.mean(c_perm))
    r["null_analytic_iso"] = float(np.sqrt(2.0 / (np.pi * m)))
    # per-coordinate magnitude agreement -- the magnitude confound, directly
    r["corr_absv_abse"] = L.pearson(np.abs(v), np.abs(e_true))

    # ================= PART B : subspaces =================
    #: the emitted structure align_jac used (incumbent projected torsions) and two controls
    rec = I.shipped_record(p)
    iphi, ipsi = np.asarray(rec["phi"], float), np.asarray(rec["psi"], float)
    Ji, _, _ = A.sup_jacobian(iphi, ipsi)
    Jh, _, _ = A.sup_jacobian(np.full(n, ALPHA[0]), np.full(n, ALPHA[1]))
    k_rand = int(rng.integers(0, P7.shape[0]))
    Jw, _, _ = A.sup_jacobian(P7[k_rand], S7[k_rand])

    s_n, V_n = A.spectrum(Jn)
    spec = {"incumbent": Ji, "fit": Jf, "helix_ZEROINFO": Jh, "poolwindow_ZEROINFO": Jw}
    inert = [0, n - 1, n, m - 1]                          # phi0, phi_{n-1}, psi0, psi_{n-1}
    E4 = np.zeros((m, 4))
    for c_, k_ in enumerate(inert):
        E4[k_, c_] = 1.0

    for tag, Jx in spec.items():
        s_x, V_x = A.spectrum(Jx)
        r[f"speccorr_{tag}"] = L.pearson(s_n, s_x)
        for frac, lab in ((0.5, "half"), (0.25, "quarter"), (0.125, "eighth")):
            k = max(1, int(round(frac * m)))
            ca = principal_angles(V_n[:, m - k:], V_x[:, m - k:])
            r[f"ov_{tag}_{lab}"] = float((ca ** 2).mean())
            r[f"ovnull_analytic_{lab}"] = float(k) / m
        # the exactly-null 5
        k5 = int((s_n <= 1e-8 * s_n[0]).sum())
        k5x = int((s_x <= 1e-8 * s_x[0]).sum())
        r["n_exact_null_native"] = k5
        r[f"n_exact_null_{tag}"] = k5x
        kk = min(k5, k5x)
        if kk > 0:
            ca = principal_angles(V_n[:, m - kk:], V_x[:, m - kk:])
            r[f"ov_{tag}_null{kk}"] = float((ca ** 2).mean())
            r[f"ovnull_analytic_null"] = float(kk) / m
        # DEFLATED: remove the 4 known-inert coordinate axes from BOTH bottom-half subspaces
        k = max(1, m // 2)
        for lab, kk2 in (("half", m // 2), ("quarter", max(1, m // 4))):
            Bn = V_n[:, m - kk2:]
            Bx = V_x[:, m - kk2:]
            Bn_d = Bn - E4 @ (E4.T @ Bn)
            Bx_d = Bx - E4 @ (E4.T @ Bx)
            Qn = la.qr(Bn_d)[0][:, :kk2 - 4] if kk2 > 4 else None
            Qx = la.qr(Bx_d)[0][:, :kk2 - 4] if kk2 > 4 else None
            if Qn is not None:
                # re-orthonormalise properly via svd of the deflated block
                Un, sn_, _ = la.svd(Bn_d, full_matrices=False)
                Ux, sx_, _ = la.svd(Bx_d, full_matrices=False)
                Un = Un[:, sn_ > 1e-8]
                Ux = Ux[:, sx_ > 1e-8]
                ca = principal_angles(Un, Ux)
                r[f"ovdefl_{tag}_{lab}"] = float((ca ** 2).mean())
                r[f"ovdefl_null_{lab}"] = float(min(Un.shape[1], Ux.shape[1])) / (m - 4)

    #: the FIFTH exactly-null direction, native J -- distributed whole-chain crank?
    Vn5 = V_n[:, s_n <= 1e-8 * s_n[0]]
    r["inert_in_null_native"] = [float((Vn5[k_] ** 2).sum()) for k_ in inert] \
        if Vn5.shape[1] else [0.0] * 4
    if Vn5.shape[1] >= 5:
        Rres = Vn5 - E4 @ (E4.T @ Vn5)
        U5, s5, _ = la.svd(Rres, full_matrices=False)
        r["fifth_residual_sv"] = [float(x) for x in s5]
        f5 = U5[:, 0]
        f5 = f5 / la.norm(f5)
        r["fifth_sum_dphi"] = float(f5[:n].sum())
        r["fifth_sum_dpsi"] = float(f5[n:].sum())
        r["fifth_crank_ratio"] = float(-f5[n:].sum() / f5[:n].sum()) \
            if abs(f5[:n].sum()) > 1e-12 else np.nan
        r["fifth_max_abs_load"] = float(np.abs(f5).max())
        r["fifth_participation"] = float(1.0 / (f5 ** 4).sum())
        r["fifth_damage"] = float(la.norm(Jn @ f5))
        r["fifth_s1"] = float(s_n[0])

    # ================= PART C : what the ORACLE rotation constructs =================
    e = e_true
    kloud = max(1, int(round(0.5 * m)))
    P = V_f[:, :kloud]                                    # loud half at the FIT's J
    e2 = e - P @ (P.T @ e)
    ne2 = float(la.norm(e2))
    th_kill = thn + e2
    r["C_kill_loud_0.5"] = float(I.ca_rmsd(I.build_ca(th_kill[:n], th_kill[n:]), nat))
    if ne2 > 1e-12:
        th_kn = thn + e2 * (ne / ne2)
        r["C_kill_loud_norm_0.5"] = float(I.ca_rmsd(I.build_ca(th_kn[:n], th_kn[n:]), nat))
    #: THE CONTROL -- same |e|, same quiet subspace, ZERO error-direction information
    Q = V_f[:, kloud:]
    vals = []
    for _ in range(20):
        g = rng.normal(size=Q.shape[1])
        u = Q @ (g / la.norm(g))
        u = u / la.norm(u)
        thr = thn + u * ne
        vals.append(float(I.ca_rmsd(I.build_ca(thr[:n], thr[n:]), nat)))
    r["C_quiet_RANDOM_norm_0.5"] = float(np.mean(vals))
    r["C_quiet_RANDOM_norm_sd"] = float(np.std(vals))
    #: the same with the loud half instead, as the opposite bound
    vals = []
    for _ in range(20):
        g = rng.normal(size=P.shape[1])
        u = P @ (g / la.norm(g))
        u = u / la.norm(u)
        thr = thn + u * ne
        vals.append(float(I.ca_rmsd(I.build_ca(thr[:n], thr[n:]), nat)))
    r["C_loud_RANDOM_norm_0.5"] = float(np.mean(vals))
    #: and an isotropic direction at the same |e|, the plain magnitude control
    vals = []
    for _ in range(20):
        g = rng.normal(size=m)
        u = g / la.norm(g)
        thr = thn + u * ne
        vals.append(float(I.ca_rmsd(I.build_ca(thr[:n], thr[n:]), nat)))
    r["C_iso_RANDOM_norm"] = float(np.mean(vals))
    #: the LIMIT: use J at the NATIVE (ORACLE subspace) -- should collapse toward 0
    Pn = V_n[:, :kloud]
    e2n = e - Pn @ (Pn.T @ e)
    if la.norm(e2n) > 1e-12:
        thx = thn + e2n * (ne / la.norm(e2n))
        r["C_kill_loud_norm_NATIVE_J"] = float(
            I.ca_rmsd(I.build_ca(thx[:n], thx[n:]), nat))

    # ================= PART D : objective vs RMSD sensitivity =================
    Gp, dfit = A.pair_jacobian(Jraw_f, CAf, i, j)          # (P, 2n) d d_p / d theta
    inv = 1.0 / np.asarray(sd, float)
    #: Gauss-Newton curvature of F = sum_p w_p ((d_p-dhat_p)/sd_p)^2 :  H = 2 sum w/sd^2 g g'
    Wg = (w0 * inv * inv)[:, None] * Gp
    H = 2.0 * (Gp.T @ Wg)
    #: RMSD^2 curvature along v is ||J v||^2 / n
    Kr = (Jf.T @ Jf) / n
    quiet_idx = np.arange(kloud, m)
    loud_idx = np.arange(0, kloud)
    obj_q = float(np.mean([V_f[:, q] @ H @ V_f[:, q] for q in quiet_idx]))
    obj_l = float(np.mean([V_f[:, q] @ H @ V_f[:, q] for q in loud_idx]))
    rms_q = float(np.mean([V_f[:, q] @ Kr @ V_f[:, q] for q in quiet_idx]))
    rms_l = float(np.mean([V_f[:, q] @ Kr @ V_f[:, q] for q in loud_idx]))
    r["D_obj_curv_quiet"] = obj_q
    r["D_obj_curv_loud"] = obj_l
    r["D_rmsd_curv_quiet"] = rms_q
    r["D_rmsd_curv_loud"] = rms_l
    r["D_obj_ratio_q_over_l"] = obj_q / obj_l if obj_l > 0 else np.nan
    r["D_rmsd_ratio_q_over_l"] = rms_q / rms_l if rms_l > 0 else np.nan

    #: FINITE STEP at a small (5 deg, linear regime) and a realistic (20 deg) RMS,
    #: exact rebuild, no linearisation.
    def fobj(x):
        ca = np.asarray(I.build_ca(x[:n], x[n:]), float)
        dd = np.sqrt(((ca[i] - ca[j]) ** 2).sum(1))
        return float((w0 * ((dd - dhat) * inv) ** 2).sum()), float(I.ca_rmsd(ca, nat))
    f0, rm0 = fobj(th)
    r["D_step_f0"] = f0
    r["D_step_rmsd0"] = rm0
    for deg in (5.0, 20.0):
        step = np.deg2rad(deg) * np.sqrt(m)                # ||delta|| for `deg` RMS
        dfq, drq, dfl, drl = [], [], [], []
        for _ in range(20):
            g = rng.normal(size=Q.shape[1]); u = Q @ (g / la.norm(g)); u /= la.norm(u)
            f1, rm1 = fobj(th + step * u)
            dfq.append(f1 - f0); drq.append(rm1 - rm0)
            g = rng.normal(size=P.shape[1]); u = P @ (g / la.norm(g)); u /= la.norm(u)
            f1, rm1 = fobj(th + step * u)
            dfl.append(f1 - f0); drl.append(rm1 - rm0)
        d_ = int(deg)
        r[f"D{d_}_dobj_quiet"] = float(np.mean(dfq))
        r[f"D{d_}_dobj_loud"] = float(np.mean(dfl))
        r[f"D{d_}_drmsd_quiet"] = float(np.mean(drq))
        r[f"D{d_}_drmsd_loud"] = float(np.mean(drl))
    return r


def main(limit=None):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)
    folds_all = sorted({data[p]["fold"] for p in pdbs})
    deb = {}
    for f in folds_all:
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        rows.append(per_target(t, data, deb))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(PATH, "w") as fh:
                json.dump({"rows": rows, "n_done": len(rows), "n_expected": len(tg)}, fh,
                          default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
            el = time.time() - t0
            print(f"  {c+1}/{len(tg)}  ({el/60:.1f} min, {el/(c+1):.1f} s/target)", flush=True)
    report(rows)
    return rows


def col(rows, k):
    return np.asarray([r.get(k, np.nan) for r in rows], float)


def report(rows):
    n = len(rows)
    folds = np.asarray([r["fold"] for r in rows], int)
    out = {"n": n, "rows": rows}
    P = print
    P(f"\n{'='*100}\nSPRINT 16 AUDIT -- n = {n} targets, statistical unit = TARGET\n{'='*100}")

    # ---------------- PART A
    P("\nPART A -- THE ERROR-DIRECTION SURROGATE\n")
    P(f"{'quantity':<44}{'mean':>9}{'median':>9}{'fold-clustered 95% CI':>28}")
    A_ = {}
    for k, lab in (("cos_d_fit_pool", "|cos| theta_fit-theta_pool (RECORDED 0.390)"),
                   ("cossgn_d_fit_pool", "  signed cos, same pair"),
                   ("cos_d_fit_start", "|cos| theta_fit-theta_start (rec 0.354)"),
                   ("cos_d_fit_helix_ZEROINFO", "|cos| theta_fit-alpha-helix ZERO-INFO"),
                   ("cos_d_fit_randwindow_ZEROINFO", "|cos| theta_fit-one random window ZERO-INFO"),
                   ("null_analytic_iso", "null: isotropic, analytic sqrt(2/pi*2n)"),
                   ("null_iso", "null: isotropic, sampled"),
                   ("null_signflip_MATCHED", "MATCHED null: random signs, same |v_k|"),
                   ("null_perm_MATCHED", "MATCHED null: permuted within phi/psi"),
                   ("coslive_d_fit_pool", "|cos| on the 2n-3 coords the builder READS"),
                   ("cosnonnull_d_fit_pool", "|cos| off the exact null space of J(fit)"),
                   ("cos_unwrapped_d_fit_pool", "|cos| with NO wrapping anywhere"),
                   ("cos_wrapped_sur_unwrapped_e", "|cos| wrapped surrogate, unwrapped e"),
                   ("corr_absv_abse", "per-coord corr(|v_k|,|e_k|)  MAGNITUDE"),
                   ):
        v = col(rows, k)
        ok = np.isfinite(v)
        ci = cluster_boot_mean(v[ok], folds[ok])
        A_[k] = {"mean": float(v[ok].mean()), "median": float(np.median(v[ok])),
                 "ci95_foldclustered": ci, "n": int(ok.sum())}
        P(f"{lab:<44}{v[ok].mean():>9.3f}{np.median(v[ok]):>9.3f}"
          f"    [{ci[0]:+.3f}, {ci[1]:+.3f}]")

    at = col(rows, "align_true_ORACLE")
    ne = col(rows, "norm_e_ORACLE")
    tr = col(rows, "tors_rms_deg")
    P(f"\n{'alignment tracking':<44}{'corr':>9}{'iid CI':>20}{'fold-clustered CI':>24}")
    for k, lab in (("align_d_fit_pool", "corr(true align, surrogate align) REC +0.499"),
                   ("align_d_fit_start", "  theta_fit-theta_start"),
                   ("align_d_fit_helix_ZEROINFO", "  alpha-helix ZERO-INFO reference"),
                   ("align_d_fit_randwindow_ZEROINFO", "  one random window ZERO-INFO ref")):
        a = col(rows, k)
        ok = np.isfinite(a) & np.isfinite(at)
        rr = L.pearson(at[ok], a[ok])
        ci_i = L.boot_ci(at[ok], a[ok])
        ci_f = cluster_boot_corr(at[ok], a[ok], folds[ok])
        pc = partial_corr(at[ok], a[ok], ne[ok])
        pc2 = partial_corr(at[ok], a[ok], ne[ok], tr[ok])
        A_[k] = {"corr": rr, "ci95_iid": ci_i, "ci95_foldclustered": ci_f,
                 "partial_given_norm_e": pc, "partial_given_norm_e_and_tors": pc2}
        P(f"{lab:<44}{rr:>+9.3f}  [{ci_i[0]:+.3f},{ci_i[1]:+.3f}]   "
          f"[{ci_f[0]:+.3f},{ci_f[1]:+.3f}]   partial|||e|| {pc:+.3f}")

    cp = col(rows, "cos_d_fit_pool")
    P(f"\n  corr(|cos|, ||e|| ORACLE) = {L.pearson(cp, ne):+.3f}   "
      f"corr(|cos|, torsion RMS) = {L.pearson(cp, tr):+.3f}   "
      f"corr(|cos|, RMSD) = {L.pearson(cp, col(rows,'rmsd_squared')):+.3f}")
    d = cp - col(rows, "null_signflip_MATCHED")
    P(f"  |cos| MINUS its magnitude-matched sign-flip null, paired: "
      f"{d.mean():+.4f}  fold-clustered CI {cluster_boot_mean(d, folds)}  "
      f"W/L {(d>0).sum()}/{(d<0).sum()}")
    A_["cos_minus_signflip"] = {"mean": float(d.mean()),
                                "ci95_foldclustered": cluster_boot_mean(d, folds),
                                "n_better": int((d > 0).sum()), "n_worse": int((d < 0).sum())}
    out["A"] = A_

    # ---------------- PART B
    P("\n\nPART B -- THE QUIET-SUBSPACE ESTIMATE\n")
    P(f"{'J compared against J(native)':<34}{'half':>9}{'quarter':>9}{'eighth':>9}"
      f"{'null-5':>9}{'spec corr':>11}")
    B_ = {}
    for tag, lab in (("incumbent", "incumbent emitted (RECORDED)"),
                     ("fit", "the restraint fit itself"),
                     ("helix_ZEROINFO", "constant alpha-helix ZERO-INFO"),
                     ("poolwindow_ZEROINFO", "a random pool window ZERO-INFO")):
        vals = []
        for lab2 in ("half", "quarter", "eighth"):
            vals.append(np.nanmean(col(rows, f"ov_{tag}_{lab2}")))
        nl = np.nanmean([r[k] for r in rows for k in r if k.startswith(f"ov_{tag}_null")])
        sc = np.nanmean(col(rows, f"speccorr_{tag}"))
        B_[tag] = {"half": float(vals[0]), "quarter": float(vals[1]),
                   "eighth": float(vals[2]), "null5": float(nl), "spec_corr": float(sc)}
        P(f"{lab:<34}{vals[0]:>9.3f}{vals[1]:>9.3f}{vals[2]:>9.3f}{nl:>9.3f}{sc:>+11.3f}")
    nh = np.nanmean(col(rows, "ovnull_analytic_half"))
    nq = np.nanmean(col(rows, "ovnull_analytic_quarter"))
    ne8 = np.nanmean(col(rows, "ovnull_analytic_eighth"))
    n5 = np.nanmean(col(rows, "ovnull_analytic_null"))
    P(f"{'ANALYTIC random-subspace null k/m':<34}{nh:>9.3f}{nq:>9.3f}{ne8:>9.3f}{n5:>9.3f}")
    B_["analytic_null"] = {"half": float(nh), "quarter": float(nq), "eighth": float(ne8),
                           "null5": float(n5)}
    P("\n  DEFLATED (the 4 known-inert coordinate axes removed from BOTH subspaces):")
    P(f"{'':<34}{'half':>9}{'quarter':>9}{'null(h)':>10}{'null(q)':>10}")
    for tag, lab in (("incumbent", "incumbent emitted"), ("fit", "the fit"),
                     ("helix_ZEROINFO", "alpha-helix ZERO-INFO"),
                     ("poolwindow_ZEROINFO", "random pool window")):
        h = np.nanmean(col(rows, f"ovdefl_{tag}_half"))
        q = np.nanmean(col(rows, f"ovdefl_{tag}_quarter"))
        B_[tag]["defl_half"] = float(h); B_[tag]["defl_quarter"] = float(q)
        P(f"{lab:<34}{h:>9.3f}{q:>9.3f}"
          f"{np.nanmean(col(rows,'ovdefl_null_half')):>10.3f}"
          f"{np.nanmean(col(rows,'ovdefl_null_quarter')):>10.3f}")
    nn = col(rows, "n_exact_null_native")
    P(f"\n  exactly-null directions of J(native): mean {np.nanmean(nn):.3f}  "
      f"all==5 on {int((nn==5).sum())}/{n} targets")
    iin = np.asarray([r["inert_in_null_native"] for r in rows], float)
    P(f"  the 4 coordinate inert torsions' projection into the null space: "
      f"min {iin.min():.6f} over 126x4")
    fr = col(rows, "fifth_crank_ratio")
    P(f"  the FIFTH direction after deflating the 4: sum d-psi / -sum d-phi = "
      f"{np.nanmean(fr):.3f} (median {np.nanmedian(fr):.3f}), "
      f"max |load| {np.nanmean(col(rows,'fifth_max_abs_load')):.3f}, "
      f"participation {np.nanmean(col(rows,'fifth_participation')):.2f} coords, "
      f"damage ||J f5|| {np.nanmean(col(rows,'fifth_damage')):.2e} vs s1 "
      f"{np.nanmean(col(rows,'fifth_s1')):.3f}")
    B_["n_null_is_5"] = int((nn == 5).sum())
    B_["fifth_crank_ratio_mean"] = float(np.nanmean(fr))
    B_["fifth_max_abs_load"] = float(np.nanmean(col(rows, "fifth_max_abs_load")))
    B_["fifth_participation"] = float(np.nanmean(col(rows, "fifth_participation")))
    out["B"] = B_

    # ---------------- PART C
    P("\n\nPART C -- WHAT THE ORACLE ROTATION CONSTRUCTS\n")
    P(f"{'arm  (all built as theta_NATIVE + delta)':<52}{'RMSD':>9}{'median':>9}"
      f"{'  vs squared fit':<30}")
    base = col(rows, "rmsd_squared")
    C_ = {}
    for k, lab in (("rmsd_squared", "squared fit itself (RECORDED 3.677)"),
                   ("C_kill_loud_0.5", "ORACLE_kill_loud_0.5 (rec 1.443)"),
                   ("C_kill_loud_norm_0.5", "ORACLE_kill_loud_norm_0.5 (rec 1.855)"),
                   ("C_quiet_RANDOM_norm_0.5", "CONTROL: |e| x RANDOM dir in the SAME quiet half"),
                   ("C_loud_RANDOM_norm_0.5", "CONTROL: |e| x RANDOM dir in the loud half"),
                   ("C_iso_RANDOM_norm", "CONTROL: |e| x isotropic random direction"),
                   ("C_kill_loud_norm_NATIVE_J", "the same but subspace from J(NATIVE) -- limit")):
        v = col(rows, k)
        ok = np.isfinite(v)
        pr = I.paired(v[ok], base[ok], folds=folds[ok])
        C_[k] = {"mean": float(v[ok].mean()), "median": float(np.median(v[ok])),
                 "vs_squared": pr}
        P(f"{lab:<52}{v[ok].mean():>9.3f}{np.median(v[ok]):>9.3f}"
          f"  {pr['mean_diff']:+.3f} [{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}]"
          f" W/L {pr['n_better']}/{pr['n_worse']}")
    a = col(rows, "C_kill_loud_norm_0.5"); b = col(rows, "C_quiet_RANDOM_norm_0.5")
    pr = I.paired(a, b, folds=folds)
    P(f"\n  ORACLE rotation MINUS its own zero-information control: "
      f"{pr['mean_diff']:+.3f} [{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] "
      f"W/L {pr['n_better']}/{pr['n_worse']}")
    C_["oracle_minus_random_quiet"] = pr
    out["C"] = C_

    # ---------------- PART D
    P("\n\nPART D -- OBJECTIVE vs RMSD SENSITIVITY ALONG QUIET AND LOUD DIRECTIONS\n")
    oq, ol = col(rows, "D_obj_curv_quiet"), col(rows, "D_obj_curv_loud")
    rq, rl = col(rows, "D_rmsd_curv_quiet"), col(rows, "D_rmsd_curv_loud")
    D_ = {}
    P(f"  Gauss-Newton curvature of the restraint OBJECTIVE, quiet/loud: "
      f"{np.mean(oq/ol):.5f}  (median {np.median(oq/ol):.5f})")
    P(f"  curvature of RMSD^2,                        quiet/loud: "
      f"{np.mean(rq/rl):.5f}  (median {np.median(rq/rl):.5f})")
    P(f"  RATIO OF THE TWO RATIOS (objective-quiet-share / RMSD-quiet-share): "
      f"{np.mean(oq/ol)/np.mean(rq/rl):.3f}")
    D_["obj_ratio"] = float(np.mean(oq / ol))
    D_["rmsd_ratio"] = float(np.mean(rq / rl))
    D_["ratio_of_ratios"] = float(np.mean(oq / ol) / np.mean(rq / rl))
    f0 = col(rows, "D_step_f0")
    P(f"\n  FINITE STEP, exact rebuild (no linearisation). Objective at the fit: "
      f"{np.mean(f0):.2f}  (RMSD at the fit {np.mean(col(rows,'D_step_rmsd0')):.3f} A)")
    for deg in (5, 20):
        fq, fl = col(rows, f"D{deg}_dobj_quiet"), col(rows, f"D{deg}_dobj_loud")
        sq, sl = col(rows, f"D{deg}_drmsd_quiet"), col(rows, f"D{deg}_drmsd_loud")
        ro = np.mean(fq) / np.mean(fl)
        rr = np.mean(sq) / np.mean(sl)
        P(f"   {deg:>3d} deg RMS step:")
        P(f"    d(objective)  quiet {np.mean(fq):>12.3f}   loud {np.mean(fl):>12.3f}"
          f"   quiet/loud {ro:.4f}   quiet step as a fraction of the objective itself "
          f"{np.mean(fq)/np.mean(f0):.3f}")
        P(f"    d(CA-RMSD)/A  quiet {np.mean(sq):>12.3f}   loud {np.mean(sl):>12.3f}"
          f"   quiet/loud {rr:.4f}")
        P(f"    RATIO OF THE TWO RATIOS (objective/RMSD): {ro/rr:.3f}")
        D_[f"step{deg}"] = {"dobj_quiet": float(np.mean(fq)), "dobj_loud": float(np.mean(fl)),
                            "drmsd_quiet": float(np.mean(sq)), "drmsd_loud": float(np.mean(sl)),
                            "obj_ratio": float(ro), "rmsd_ratio": float(rr),
                            "ratio_of_ratios": float(ro / rr),
                            "dobj_quiet_over_f0": float(np.mean(fq) / np.mean(f0))}
    out["D"] = D_

    with open(PATH, "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    P(f"\nwritten {PATH}")
    return out


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
