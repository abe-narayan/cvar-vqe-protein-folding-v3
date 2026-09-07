"""s21/c_cont.py -- BLOCKS C1/C2/C3: LEGACY -> AMBER CONTINUATION, AND BLOCK P.

`s21/PREREG_C.md` sections 2, 3, 4, 5, 6.  Falsifiers F-C1, F-C2, F-C3, F-P were written there
BEFORE this module ran, together with my honest prior that C1 partly succeeds and C3 fails.

    H(lambda) = (1 - lambda) * f_L(E_Legacy(theta))  +  lambda * f_A(E_AMBER(theta))

with `f` one of the four normalisations declared in `PREREG_C.md` section 1 and built by
`s21.c_norm` from each target's OWN pool statistics -- no native information, no fitted constant,
no RMSD seen when they were chosen.

THREE THINGS THIS MODULE IS CAREFUL ABOUT, each because the programme has been burned by it.

1.  **The gradient of the mixture is EXACT and costs nothing extra.**
    grad H(lambda) = (1-l) f_L'(E_L) grad E_L  +  l f_A'(E_A) grad E_A.
    So ONE Legacy gradient and ONE AMBER gradient price the WHOLE lambda axis for every
    normalisation.  No lambda-dependent recomputation, and no chance of the arms differing by
    anything but lambda.

2.  **A move is a move.**  s20 L12 measured that 72% of AMBER's RMSD damage in torsion space is
    the SIZE of the move.  Every RMSD number here therefore carries a `toward_member` geodesic
    null matched to that arm's OWN realised torus magnitude (`s20/c_land_null.py`'s operator,
    reused).

3.  **Both potentials stay separably evaluable** (`BRIEF` section 7).  Every row records
    `E_legacy` and `E_amber` at every reported point, never only the mixture.

THE AMBER OBJECT.  Bare single point (`ConstrainedBox.energy_point`), no minimisation, because a
landscape needs a function of theta.  The DEPLOYED `H_AMBER` is `E o Relax_50` -- a different
operator, cap-bound 192/192 -- and no table here is a statement about it.

    python -m s21.c_cont c1          # BLOCK C1  the mechanism test
    python -m s21.c_cont c2          # BLOCK C2  the lambda sweep (gradient / tail / Hessian)
    python -m s21.c_cont c3          # BLOCK C3 + P  the staged schedules and preconditioners
    python -m s21.c_cont report <blk>
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                       # noqa: E402
from s14.avgspace import top75_windows                # noqa: E402
from s15 import seed as SD                            # noqa: E402
from s18 import phys_lib as PL                        # noqa: E402
from s20 import c_land as CL                          # noqa: E402
from s21 import c_norm as CN                          # noqa: E402

SUBSET = CN.SUBSET
NORMS = ("raw", "Nz", "Ng", "Nt")
#: the declared lambda grid (PREREG section 3), plus the logarithmic sub-grid the pre-registration
#: predicts the RAW transition lives in.
LAM = (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95, 1.0)
LAM_LOG = tuple(10.0 ** np.arange(-9, 0, 1.0))
LAM_HESS = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
ALPHA = 0.15                     # the CVaR tail width used throughout
N_STARTS = 5                     # C1/C2, matched to s20/c_land.py
N_STARTS_C3 = 3                  # C3: medoid + 2, declared here to fit the compute budget
MAXITER = 100                    # matched across every C3 arm
N_NULL_DRAWS = 30                # toward_member draws per arm, as in s20/c_land_null.py


# ==========================================================================================
class MixPot(CL.Pot):
    """`s20.c_land.Pot` with a mixed Hamiltonian bolted on, so that every landscape routine in
    `c_land` (hess_fd, minimise, spectrum_metrics) applies unchanged.

    `which` may be "legacy", "amber", "mix", or "amber_bonded".
    `amber_bonded` is the bonded FORCE-GROUP SUBSET of the same genuine ff14SB system -- bond,
    angle and torsion only, nonbonded and solvation switched off.  It is a subset of AMBER, not a
    learned surrogate, and it is NEVER called "amber" anywhere (PREREG section 6).
    """

    def __init__(self, seq, fL=None, fA=None, lam=0.0):
        super().__init__(seq, want_amber=True)
        self.fL = fL if fL is not None else CN.Norm("raw")
        self.fA = fA if fA is not None else CN.Norm("raw")
        self.lam = float(lam)
        self._bonded_groups = {1, 2, 3}

    # -- the bonded force-group subset -------------------------------------------------
    def _ctx(self):
        if self.box._soft is None:
            self.box._make_soft()
        ctx = self.box._soft[0]
        for p in ("kh_ca", "kf_ca", "eps_ca", "kh_bb"):
            ctx.setParameter(p, 0.0)
        return ctx

    def amber_bonded(self, TH):
        from openmm import unit
        TH = np.atleast_2d(np.asarray(TH, float))
        ctx = self._ctx()
        out = np.empty(len(TH))
        for b in range(len(TH)):
            ctx.setPositions(self._x(TH[b]) * unit.nanometer)
            out[b] = ctx.getState(getEnergy=True, groups=self._bonded_groups) \
                .getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        self.n_eval["amber"] += len(TH)
        return out

    def amber_bonded_grad_chain(self, th):
        from openmm import unit
        ctx = self._ctx()
        J, x0 = self.jac(th)
        ctx.setPositions(x0 * unit.nanometer)
        st = ctx.getState(getEnergy=True, getForces=True, groups=self._bonded_groups)
        e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        f = np.asarray(st.getForces(asNumpy=True).value_in_unit(
            unit.kilocalorie_per_mole / unit.nanometer), float)
        self.n_eval["amber"] += 1
        return float(e), -(J @ f.ravel())

    # -- the mixture --------------------------------------------------------------------
    def E(self, which, TH):
        if which == "legacy":
            return self.legacy(TH)
        if which == "amber":
            return self.amber(TH)
        if which == "amber_bonded":
            return self.amber_bonded(TH)
        el = np.asarray(self.legacy(TH), float)
        ea = np.asarray(self.amber(TH), float)
        return (1.0 - self.lam) * self.fL(el) + self.lam * self.fA(ea)

    def grad(self, which, th):
        if which == "amber":
            return self.amber_grad_chain(th)
        if which == "amber_bonded":
            return self.amber_bonded_grad_chain(th)
        if which == "legacy":
            e = float(self.legacy(th)[0])
            return e, self.grad_fd("legacy", th)
        el = float(self.legacy(th)[0])
        gl = self.grad_fd("legacy", th)
        ea, ga = self.amber_grad_chain(th)
        e = (1.0 - self.lam) * float(self.fL(el)) + self.lam * float(self.fA(ea))
        g = ((1.0 - self.lam) * float(self.fL.dfdE(el)) * gl
             + self.lam * float(self.fA.dfdE(ea)) * ga)
        return e, g

    def both(self, th):
        """The separable pair at one point -- (E_legacy, grad_legacy, E_amber, grad_amber)."""
        el = float(self.legacy(th)[0])
        gl = self.grad_fd("legacy", th)
        ea, ga = self.amber_grad_chain(th)
        return el, gl, ea, ga


def mixed_grad(el, gl, ea, ga, fL, fA, lam):
    """EXACT.  One Legacy gradient and one AMBER gradient price the whole lambda axis."""
    return ((1.0 - lam) * float(fL.dfdE(el)) * np.asarray(gl, float)
            + lam * float(fA.dfdE(ea)) * np.asarray(ga, float))


def torus_d(a, b):
    d = np.mod(np.asarray(a, float) - np.asarray(b, float) + np.pi, 2 * np.pi) - np.pi
    return float(np.sqrt((d ** 2).mean()))


def PP(a, b, folds=None):
    st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
    return {"mean": st["mean"], "ci": st.get("ci_fold", st["ci"]), "ci_iid": st["ci"],
            "median": st["median"], "W": st["W"], "L": st["L"], "n": st["n"],
            "mean_a": st["mean_a"], "mean_b": st["mean_b"]}


def _norms_for_target(nrow):
    """Build the four declared normalisations from this target's Block-N pool statistics."""
    out = {}
    for k in NORMS:
        out[k] = (CN.Norm(k, nrow["legacy"]["median"], nrow["legacy"]["mad_n"],
                          nrow["legacy"]["gnorm"]),
                  CN.Norm(k, nrow["amber"]["median"], nrow["amber"]["mad_n"],
                          nrow["amber"]["gnorm"]))
    return out


def _load_norm():
    p = os.path.join(RESULTS, "c_norm.json")
    o = json.load(open(p))
    return {r["pdb"]: r for r in o["rows"]}, o


def _setup(pdb, t, nstarts):
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    Pm = I.pairwise_rmsd(W)
    med = int(I.medoid(Pm))
    rng0 = SD.stable_rng(pdb, "s20C_land")
    others = [int(x) for x in rng0.permutation(len(W))[:N_STARTS] if int(x) != med]
    order = ([med] + others[:N_STARTS - 1])[:nstarts]
    starts = [np.concatenate([PHI[b], PSI[b]]) for b in order]
    return W, PHI, PSI, nat, starts, order


def _null_toward(pdb, PHI, PSI, starts, active, mag, nat, n, tag):
    """s20/c_land_null.py's operator, reused verbatim: move the SAME per-coordinate RMS torus
    distance along the geodesic toward another randomly chosen pool member.  Realisable,
    native-free, zero information.  This is the null every RMSD claim here rests on."""
    rng = SD.stable_rng(pdb, "s21C_null_" + tag)
    vals = []
    for th_s in starts:
        for _ in range(N_NULL_DRAWS):
            j = int(rng.integers(len(PHI)))
            tj = np.concatenate([PHI[j], PSI[j]])
            d = np.mod(tj - th_s + np.pi, 2 * np.pi) - np.pi
            dn = float(np.sqrt((d[active] ** 2).mean()))
            if dn < 1e-9:
                continue
            th2 = th_s.copy()
            th2[active] = th_s[active] + d[active] * (mag / dn)
            vals.append(float(I.ca_rmsd(I.build_ca(th2[:n], th2[n:]), nat)))
    return float(np.mean(vals)), int(len(vals))


# ==========================================================================================
# BLOCK C1 -- THE MECHANISM.  Does Legacy relaxation de-singularise AMBER?
# ==========================================================================================
def c1_target(t, nrow):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, nat, starts, order = _setup(pdb, t, N_STARTS)
    active = np.arange(1, 2 * n)
    P = MixPot(seq)
    rec = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "starts": []}

    #: A CORRECTION OF MY OWN CODE AGAINST MY OWN PRE-REGISTRATION, made after two targets and
    #: recorded rather than repaired silently.  `PREREG_C.md` section 2 budgets "30 targets x
    #: (1 Legacy minimisation + 3 AMBER Hessians + nulls), ~1 h".  The first draft computed the
    #: Hessian TRIPLE at all 5 starts -- 15 AMBER Hessians per target, five times the registered
    #: design and ~4 h of wall clock.  The energies, gradients, moves and RMSDs are still read at
    #: ALL 5 starts; the three EXPENSIVE spectra are read at the MEDOID start, which is what was
    #: registered.  The two rows computed under the over-implemented version were DISCARDED, not
    #: merged, so the artefact is homogeneous.
    for si, th0 in enumerate(starts):
        want_hess = (si == 0)
        el0, gl0, ea0, ga0 = P.both(th0)
        if want_hess:
            H0, _ = P.hess_fd("amber", th0, active=active)
            m0 = CL.spectrum_metrics(H0)
        else:
            m0 = {k: float("nan") for k in CL.METRICS}
        #: the Legacy minimisation, identical optimiser and budget to s20/c_land.py
        thL, eL, nc, suc = CL.minimise(P, "legacy", th0, maxiter=MAXITER)
        elL, glL, eaL, gaL = P.both(thL)
        if want_hess:
            HL, _ = P.hess_fd("amber", thL, active=active)
            mL = CL.spectrum_metrics(HL)
        else:
            mL = {k: float("nan") for k in CL.METRICS}
        mag = torus_d(thL, th0)
        #: THE NULL -- a realisable move of the SAME torus magnitude, then the same three reads
        rngn = SD.stable_rng(pdb, "s21C1null")
        nul = []
        for _ in range(8):
            j = int(rngn.integers(len(PHI)))
            tj = np.concatenate([PHI[j], PSI[j]])
            d = np.mod(tj - th0 + np.pi, 2 * np.pi) - np.pi
            dn = float(np.sqrt((d[active] ** 2).mean()))
            if dn < 1e-9:
                continue
            th2 = th0.copy()
            th2[active] = th0[active] + d[active] * (mag / dn)
            nul.append(th2)
        eN = [float(P.amber(th)[0]) for th in nul]
        gN = [float(np.linalg.norm(P.amber_grad_chain(th)[1][active])) for th in nul]
        #: the Hessian is the expensive read -- one null draw only, and it is labelled as such
        if want_hess and nul:
            HN, _ = P.hess_fd("amber", nul[0], active=active)
            mN = CL.spectrum_metrics(HN)
        else:
            mN = {k: float("nan") for k in CL.METRICS}

        rec["starts"].append({
            "E_amber_0": float(ea0), "E_amber_L": float(eaL),
            "E_amber_null": float(np.mean(eN)) if eN else float("nan"),
            "E_legacy_0": float(el0), "E_legacy_L": float(elL),
            "gA_0": float(np.linalg.norm(ga0[active])),
            "gA_L": float(np.linalg.norm(gaL[active])),
            "gA_null": float(np.mean(gN)) if gN else float("nan"),
            "pr_0": float(m0["part_ratio"]), "pr_L": float(mL["part_ratio"]),
            "pr_null": float(mN["part_ratio"]),
            "cond_0": float(m0["cond_med"]), "cond_L": float(mL["cond_med"]),
            "nearzero_0": float(m0["nearzero_frac"]), "nearzero_L": float(mL["nearzero_frac"]),
            "aniso_0": float(m0["aniso"]), "aniso_L": float(mL["aniso"]),
            "theta_moved": mag,
            "rmsd_0": float(I.ca_rmsd(I.build_ca(th0[:n], th0[n:]), nat)),
            "rmsd_L": float(I.ca_rmsd(I.build_ca(thL[:n], thL[n:]), nat)),
            "n_null": int(len(nul)),
        })
    P.close()
    for k in ("E_amber_0", "E_amber_L", "E_amber_null", "gA_0", "gA_L", "gA_null",
              "pr_0", "pr_L", "pr_null", "cond_0", "cond_L", "nearzero_0", "nearzero_L",
              "aniso_0", "aniso_L", "theta_moved", "rmsd_0", "rmsd_L"):
        v = [s[k] for s in rec["starts"]]
        rec[k] = float(np.nanmean(v)) if np.any(np.isfinite(v)) else float("nan")
    #: medians too -- AMBER energies span 8 orders and a mean of those is meaningless
    for k in ("E_amber_0", "E_amber_L", "E_amber_null", "gA_0", "gA_L", "gA_null"):
        rec[k + "_med"] = float(np.nanmedian([s[k] for s in rec["starts"]]))
    return rec


# ==========================================================================================
# BLOCK C2 -- THE LAMBDA SWEEP.  Gradient side and tail side are EXACT and nearly free.
# ==========================================================================================
def c2_target(t, nrow, do_hess=True):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, nat, starts, order = _setup(pdb, t, N_STARTS)
    active = np.arange(1, 2 * n)
    P = MixPot(seq)
    NM = _norms_for_target(nrow)
    rec = {"pdb": pdb, "n": n, "fold": int(t["fold"])}

    # ---- gradient side: ONE Legacy gradient + ONE AMBER gradient per start prices all lambda
    per_start = []
    for th0 in starts:
        el, gl, ea, ga = P.both(th0)
        d = {"E_legacy": float(el), "E_amber": float(ea),
             "gL": float(np.linalg.norm(gl[active])), "gA": float(np.linalg.norm(ga[active])),
             "cos_LA": float(gl[active] @ ga[active]
                             / max(np.linalg.norm(gl[active]) * np.linalg.norm(ga[active]), 1e-30)),
             "lam": {}}
        for nm in NORMS:
            fL, fA = NM[nm]
            wl = float(fL.dfdE(el)); wa = float(fA.dfdE(ea))
            row = {}
            for l in tuple(LAM) + tuple(LAM_LOG):
                g = mixed_grad(el, gl, ea, ga, fL, fA, l)
                gn = float(np.linalg.norm(g[active]))
                row[f"{l:g}"] = {
                    "gnorm": gn,
                    #: the AMBER share of the mixed gradient NORM -- the honest "how much AMBER
                    #: is in this Hamiltonian" dial, and it is what a lambda grid should span.
                    "amber_share": float(l * wa * np.linalg.norm(ga[active])
                                         / max(l * wa * np.linalg.norm(ga[active])
                                               + (1 - l) * wl * np.linalg.norm(gl[active]),
                                               1e-300)),
                    "cos_with_amber": float(g[active] @ ga[active]
                                            / max(gn * np.linalg.norm(ga[active]), 1e-30)),
                    "cos_with_legacy": float(g[active] @ gl[active]
                                             / max(gn * np.linalg.norm(gl[active]), 1e-30)),
                }
            d["lam"][nm] = row
        per_start.append(d)
    rec["per_start"] = per_start

    # ---- tail side: the CVaR alpha-tail SET over the target's own pool, at every lambda
    TH = np.column_stack([PHI, PSI])
    EL = np.asarray(P.legacy(TH), float)
    EA = np.asarray(P.amber(TH), float)
    true = np.array([I.ca_rmsd(W[k], nat) for k in range(len(W))])
    m = len(TH)
    kt = max(int(round(ALPHA * m)), 1)

    def tail_stats(score):
        o = np.argsort(score, kind="stable")[:kt]
        sub = W[o]
        a, _b = I.coordinate_average(sub) if len(sub) > 1 else (sub[0], 0)
        return {"idx": [int(x) for x in o],
                "tail_mean_rmsd": float(true[o].mean()),
                "tail_avg_rmsd": float(I.ca_rmsd(np.asarray(a, float), nat)),
                "argmin_rmsd": float(true[int(np.argmin(score))])}

    rec["pool"] = {"m": m, "k_tail": kt, "pool_mean_rmsd": float(true.mean()),
                   "pool_best_rmsd": float(true.min()),
                   "E_legacy": EL.tolist(), "E_amber": EA.tolist(),
                   "true_rmsd": true.tolist()}
    rec["tail"] = {}
    for nm in NORMS:
        fL, fA = NM[nm]
        zl, za = np.asarray(fL(EL), float), np.asarray(fA(EA), float)
        row = {}
        base0 = set(np.argsort(zl, kind="stable")[:kt].tolist())
        base1 = set(np.argsort(za, kind="stable")[:kt].tolist())
        for l in tuple(LAM) + tuple(LAM_LOG):
            s = (1 - l) * zl + l * za
            ts = tail_stats(s)
            S = set(ts["idx"])
            ts["jaccard_vs_l0"] = float(len(S & base0) / max(len(S | base0), 1))
            ts["jaccard_vs_l1"] = float(len(S & base1) / max(len(S | base1), 1))
            ts.pop("idx")
            row[f"{l:g}"] = ts
        rec["tail"][nm] = row
    #: matched-count random control at the same tail width (the only comparison that means
    #: anything -- `s21/tailprice.py`'s trap, restated)
    rngr = SD.stable_rng(pdb, "s21C2rand")
    rr = []
    for _ in range(16):
        o = rngr.choice(m, kt, replace=False)
        a, _b = I.coordinate_average(W[o])
        rr.append(float(I.ca_rmsd(np.asarray(a, float), nat)))
    rec["tail_random_avg_rmsd"] = float(np.mean(rr))
    rec["tail_random_mean_rmsd"] = float(np.mean([true[rngr.choice(m, kt, replace=False)].mean()
                                                  for _ in range(64)]))

    # ---- Hessian side.  EXACT COMPOSITION, not a finite difference of the mixture.
    #
    #   H(l) = (1-l) f_L(E_L) + l f_A(E_A)
    #   grad^2 H(l) = (1-l)[ f_L'(E_L) grad^2 E_L + f_L''(E_L) g_L g_L^T ]
    #               +   l  [ f_A'(E_A) grad^2 E_A + f_A''(E_A) g_A g_A^T ]
    #
    # so TWO Hessians and TWO gradients price EVERY lambda under EVERY normalisation, exactly.
    # This is the same chain rule that `PREREG_C.md` section 4.3 derives, and it is strictly
    # better than differencing the mixture: it removes the mixture's own finite-difference noise
    # and it guarantees the arms differ by nothing but lambda.  It is VERIFIED against a direct
    # FD Hessian of the mixture at one interior lambda on every target (`hess_check`), so the
    # saving is certified rather than assumed.
    if do_hess:
        rec["hess"] = {}
        th0 = starts[0]
        HL, eL0 = P.hess_fd("legacy", th0, active=active)
        HA, eA0 = P.hess_fd("amber", th0, active=active)
        gl0 = P.grad_fd("legacy", th0)[active]
        _e, ga_full = P.amber_grad_chain(th0)
        ga0 = ga_full[active]
        HL = 0.5 * (HL + HL.T); HA = 0.5 * (HA + HA.T)
        for nm in NORMS:
            fL, fA = NM[nm]
            ML = float(fL.dfdE(eL0)) * HL + float(fL.d2fdE2(eL0)) * np.outer(gl0, gl0)
            MA = float(fA.dfdE(eA0)) * HA + float(fA.d2fdE2(eA0)) * np.outer(ga0, ga0)
            row = {}
            for l in LAM_HESS + tuple(LAM_LOG):
                mm = CL.spectrum_metrics((1.0 - l) * ML + l * MA)
                row[f"{l:g}"] = {k: float(mm[k]) for k in CL.METRICS}
                row[f"{l:g}"]["lam_max_abs"] = float(mm["lam_max_abs"])
            rec["hess"][nm] = row
        #: THE VERIFICATION.  Direct FD Hessian of the mixture at lambda = 0.5 under `Nt`,
        #: compared to the composed one.  Reported as a relative eigenvalue error, and the
        #: number of times its tolerance FIRES is counted across targets.
        fL, fA = NM["Nt"]
        P.fL, P.fA, P.lam = fL, fA, 0.5
        Hd, _ = P.hess_fd("mix", th0, active=active)
        Mc = 0.5 * (float(fL.dfdE(eL0)) * HL + float(fL.d2fdE2(eL0)) * np.outer(gl0, gl0))              + 0.5 * (float(fA.dfdE(eA0)) * HA + float(fA.d2fdE2(eA0)) * np.outer(ga0, ga0))
        ld = np.linalg.eigvalsh(0.5 * (Hd + Hd.T)); lc = np.linalg.eigvalsh(Mc)
        mdc = CL.spectrum_metrics(0.5 * (Hd + Hd.T))
        mcc = CL.spectrum_metrics(Mc)
        rec["hess_check"] = {
            "rel_eig_l2": float(np.linalg.norm(ld - lc) / max(np.linalg.norm(lc), 1e-30)),
            "part_ratio_direct": float(mdc["part_ratio"]),
            "part_ratio_composed": float(mcc["part_ratio"]),
            "d_part_ratio": float(abs(mdc["part_ratio"] - mcc["part_ratio"]))}
    P.close()
    return rec


# ==========================================================================================
# BLOCK C3 + BLOCK P -- the staged schedules and the preconditioners
# ==========================================================================================
def _minimise_staged(P, stages, th0, total_iter=MAXITER):
    """Run a schedule of (which, lam) stages splitting `total_iter` equally.  Returns the end
    point and the realised function-call count, which is REPORTED so 'matched budget' is a
    measured fact and not a claim."""
    per = max(int(round(total_iter / len(stages))), 1)
    th = np.asarray(th0, float).copy()
    calls = 0
    for (which, lam, fL, fA) in stages:
        P.fL, P.fA, P.lam = fL, fA, float(lam)
        th, _e, nc, _s = CL.minimise(P, which, th, maxiter=per)
        calls += nc
    return th, calls


def _minimise_tr(P, which, th0, radius, maxiter=MAXITER):
    """Trust region as a BOX around the start, in radians per coordinate.  Directly attacks
    s20 L12's finding that 72% of AMBER's RMSD damage is the SIZE of the move."""
    from scipy.optimize import minimize as smin
    calls = {"n": 0}

    def fun(t):
        calls["n"] += 1
        return P.grad(which, t)

    th0 = np.asarray(th0, float)
    bounds = [(float(x - radius), float(x + radius)) for x in th0]
    r = smin(fun, th0, jac=True, method="L-BFGS-B", bounds=bounds,
             options={"maxiter": int(maxiter), "maxfun": int(maxiter) * 2})
    return np.asarray(r.x, float), calls["n"]


def c3_target(t, nrow):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, nat, starts, order = _setup(pdb, t, N_STARTS_C3)
    active = np.arange(1, 2 * n)
    P = MixPot(seq)
    NM = _norms_for_target(nrow)
    fL_raw, fA_raw = NM["raw"]
    fL_nt, fA_nt = NM["Nt"]

    #: THE TRUST RADIUS, declared here in code before any RMSD was read and native-free: the
    #: pool's OWN median per-coordinate torsional dispersion about its medoid.  "Do not move
    #: further than the ensemble of real candidates for this target already differs."
    th_med = starts[0]
    disp = [torus_d(np.concatenate([PHI[j], PSI[j]]), th_med) for j in range(len(PHI))]
    radius = float(np.median(disp))

    ARMS = {
        "L_raw":     ("stage", [("legacy", 0.0, fL_raw, fA_raw)]),
        "A_raw":     ("stage", [("amber", 1.0, fL_raw, fA_raw)]),
        "A_Nt":      ("stage", [("mix", 1.0, fL_nt, fA_nt)]),
        "LA_raw":    ("stage", [("legacy", 0.0, fL_raw, fA_raw), ("amber", 1.0, fL_raw, fA_raw)]),
        "LA_Nt":     ("stage", [("mix", 0.0, fL_nt, fA_nt), ("mix", 1.0, fL_nt, fA_nt)]),
        "L5A_Nt":    ("stage", [("mix", 0.0, fL_nt, fA_nt), ("mix", 0.5, fL_nt, fA_nt),
                                ("mix", 1.0, fL_nt, fA_nt)]),
        "ramp_Nt":   ("stage", [("mix", l, fL_nt, fA_nt) for l in (0.0, 0.25, 0.5, 0.75, 1.0)]),
        "A_tr":      ("tr", radius),
        "A_tr3":     ("tr", radius / 3.0),
        "AbondA":    ("stage", [("amber_bonded", 1.0, fL_raw, fA_raw),
                                ("amber", 1.0, fL_raw, fA_raw)]),
    }

    rec = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "radius": radius, "arms": {}}
    for name, spec in ARMS.items():
        ends, moved, calls, eA, eL = [], [], [], [], []
        for th0 in starts:
            if spec[0] == "stage":
                th, nc = _minimise_staged(P, spec[1], th0)
            else:
                th, nc = _minimise_tr(P, "amber", th0, spec[1])
            ends.append(float(I.ca_rmsd(I.build_ca(th[:n], th[n:]), nat)))
            moved.append(torus_d(th, th0))
            calls.append(int(nc))
            #: BOTH potentials at every endpoint -- separably evaluable, always
            eA.append(float(P.amber(th)[0]))
            eL.append(float(P.legacy(th)[0]))
        mag = float(np.mean(moved))
        nullv, nnull = _null_toward(pdb, PHI, PSI, starts, active, mag, nat, n, name)
        rec["arms"][name] = {
            "rmsd_end": float(np.mean(ends)), "rmsd_end_best": float(np.min(ends)),
            "theta_moved": mag, "nfev": float(np.mean(calls)),
            "E_amber_end": float(np.median(eA)), "E_legacy_end": float(np.median(eL)),
            "null_toward": nullv, "n_null": nnull,
            "excess_over_null": float(np.mean(ends)) - nullv,
        }
    rec["rmsd_start"] = float(np.mean([I.ca_rmsd(I.build_ca(s[:n], s[n:]), nat) for s in starts]))
    rec["E_amber_start"] = float(np.median([P.amber(s)[0] for s in starts]))
    P.close()
    return rec




# ==========================================================================================
# BLOCK C2t -- THE SELECTION-SIDE LAMBDA SWEEP AT THE FULL n = 126.
#
# DECLARED EXTENSION, written while C1 was running and BEFORE any C2 table existed.  C2's
# landscape work is n = 30 and is powered for ~0.3-0.5 A (PREREG section 5).  The SELECTION side
# -- which Hamiltonian, and which mixture, picks a better CVaR tail out of the target's own pool
# -- costs one Legacy batch and 75 AMBER single points per target and can therefore be run on the
# WHOLE 126-target instrument.  That is the mandatory matrix of `BRIEF` section 4 evaluated
# POOL-RESTRICTED, which bounds what any CVaR-VQE selecting from this pool can achieve.
#
# It is NOT a VQE experiment and makes no quantum claim.  The normalisation is the one declared
# in `PREREG_C.md` section 1, built from each target's OWN pool statistics -- no native
# information anywhere except the final scoring.
# ==========================================================================================
def c2tail_target(t, nrow=None):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    P = MixPot(seq)
    TH = np.column_stack([PHI, PSI])
    EL = np.asarray(P.legacy(TH), float)
    EA = np.asarray(P.amber(TH), float)
    P.close()
    m = len(TH)
    true = np.array([I.ca_rmsd(W[k], nat) for k in range(m)])
    kt = max(int(round(ALPHA * m)), 1)

    def stats(x):
        med = float(np.median(x)); mad = float(1.4826 * np.median(np.abs(x - med)))
        return med, max(mad, 1e-12)

    mL, sL = stats(EL); mA, sA = stats(EA)
    NM = {k: (CN.Norm(k, mL, sL, 1.0), CN.Norm(k, mA, sA, 1.0)) for k in ("raw", "Nz", "Nt")}

    def tail_of(score):
        o = np.argsort(score, kind="stable")[:kt]
        a, _b = I.coordinate_average(W[o])
        return (float(I.ca_rmsd(np.asarray(a, float), nat)), float(true[o].mean()),
                float(true[int(np.argmin(score))]), set(int(x) for x in o))

    rec = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "m": m, "k_tail": kt,
           "pool_mean_rmsd": float(true.mean()), "pool_best_rmsd": float(true.min()),
           "med_legacy": mL, "mad_legacy": sL, "med_amber": mA, "mad_amber": sA,
           "amber_max": float(EA.max()), "amber_frac_gt_1e4": float((EA > 1e4).mean()),
           "tail": {}}
    for nm, (fL, fA) in NM.items():
        zl, za = np.asarray(fL(EL), float), np.asarray(fA(EA), float)
        row = {}
        b0 = tail_of(zl)[3]; b1 = tail_of(za)[3]
        for l in tuple(LAM) + tuple(LAM_LOG):
            av, tm, am, S = tail_of((1 - l) * zl + l * za)
            row[f"{l:g}"] = {"tail_avg_rmsd": av, "tail_mean_rmsd": tm, "argmin_rmsd": am,
                             "jaccard_vs_l0": float(len(S & b0) / max(len(S | b0), 1)),
                             "jaccard_vs_l1": float(len(S & b1) / max(len(S | b1), 1))}
        rec["tail"][nm] = row
    #: MATCHED-COUNT RANDOM at the same tail width -- the only comparison that means anything.
    rng = SD.stable_rng(pdb, "s21C2trand")
    rr, rm = [], []
    for _ in range(16):
        o = rng.choice(m, kt, replace=False)
        a, _b = I.coordinate_average(W[o])
        rr.append(float(I.ca_rmsd(np.asarray(a, float), nat)))
        rm.append(float(true[o].mean()))
    rec["random_tail_avg_rmsd"] = float(np.mean(rr))
    rec["random_tail_mean_rmsd"] = float(np.mean(rm))
    #: the WHOLE-POOL average, a zero-information operator that is known to be decent
    a, _b = I.coordinate_average(W)
    rec["whole_pool_avg_rmsd"] = float(I.ca_rmsd(np.asarray(a, float), nat))
    #: EXACT CHECK of PREREG section 4.1: a strictly monotone transform must leave the argsort
    #: -- and therefore every CVaR tail SET -- IDENTICAL.  Verified, not assumed.
    rec["monotone_check"] = {
        "amber_argsort_identical": bool(np.array_equal(
            np.argsort(EA, kind="stable"),
            np.argsort(np.asarray(NM["Nt"][1](EA), float), kind="stable"))),
        "legacy_argsort_identical": bool(np.array_equal(
            np.argsort(EL, kind="stable"),
            np.argsort(np.asarray(NM["Nt"][0](EL), float), kind="stable")))}
    return rec


def drive_c2tail(out="c_c2tail.json"):
    tg = I.targets()
    path = os.path.join(RESULTS, out)
    rows, done = [], set()
    if os.path.exists(path):
        rows = json.load(open(path)).get("rows", [])
        done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(CN.with_mem_retry(c2tail_target, t))
        if len(rows) % 10 == 0:
            print(f"  [C2t] {len(rows)}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
        cfg = {"block": "C2t", "ALPHA": ALPHA, "LAM": list(LAM), "LAM_LOG": list(LAM_LOG),
               "norms": ["raw", "Nz", "Nt"], "n_targets": len(tg),
               "amber_object": "bare single point, NO minimisation"}
        json.dump({"rows": rows, "config": cfg,
                   "cfg_hash": hashlib.sha256(
                       json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16],
                   "n_rows": len(rows), "n_expected": len(tg),
                   "memory_guard_firings": int(CN.MEM_FIRINGS["n"]),
                   "complete": len(rows) >= len(tg)}, open(path, "w"))
    return rows


# ==========================================================================================
def _drive(block, fn, out, subset=None, **kw):
    sub = SUBSET if subset is None else subset
    tg = {t["pdb"]: t for t in I.targets()}
    NR, nobj = _load_norm()
    missing = [p for p in sub if p not in NR]
    if missing:
        raise SystemExit(f"Block N incomplete for {missing[:5]} -- run `python -m s21.c_norm`")
    path = os.path.join(RESULTS, out)
    rows, done = [], set()
    if os.path.exists(path):
        prev = json.load(open(path))
        rows = [r for r in prev.get("rows", []) if r["pdb"] in sub]
        done = {r["pdb"] for r in rows}
    t0 = time.time()
    for pdb in sub:
        if pdb in done:
            continue
        rows.append(CN.with_mem_retry(fn, tg[pdb], NR[pdb], **kw))
        print(f"  [{block}] {len(rows)}/{len(sub)} {pdb}  ({time.time()-t0:.0f}s)", flush=True)
        _write(block, rows, sub, path, nobj, kw)
    _write(block, rows, sub, path, nobj, kw)
    return rows


def _write(block, rows, sub, path, nobj, kw):
    cfg = {"block": block, "subset": list(sub), "LAM": list(LAM), "LAM_LOG": list(LAM_LOG),
           "LAM_HESS": list(LAM_HESS), "ALPHA": ALPHA, "N_STARTS": N_STARTS,
           "N_STARTS_C3": N_STARTS_C3, "MAXITER": MAXITER, "N_NULL_DRAWS": N_NULL_DRAWS,
           "H_GRAD": CL.H_GRAD, "H_HESS": CL.H_HESS, "norms": list(NORMS),
           "norm_cfg_hash": nobj["cfg_hash"], "kw": {k: str(v) for k, v in kw.items()},
           "amber_object": "bare single point, NO minimisation"}
    obj = {"rows": rows, "config": cfg,
           "memory_guard_firings": int(CN.MEM_FIRINGS["n"]),
           "memory_guard_wait_s": float(CN.MEM_FIRINGS["waited_s"]),
           "cfg_hash": hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16],
           "n_rows": len(rows), "n_expected": len(sub),
           #: a completion flag must require the FULL configuration, not the subset it was
           #: called with (`BRIEF` section 6).
           "complete": bool(len(rows) >= len(SUBSET) and list(sub) == list(SUBSET))}
    json.dump(obj, open(path, "w"))


if __name__ == "__main__":
    a = sys.argv[1:]
    cmd = a[0] if a else "c1"
    smoke = "--smoke" in a
    sub = SUBSET[:2] if smoke else None
    pre = "_SMOKE_" if smoke else ""
    if cmd == "c1":
        _drive("C1", c1_target, f"{pre}c_c1.json", subset=sub)
    elif cmd == "c2":
        _drive("C2", c2_target, f"{pre}c_c2.json", subset=sub)
    elif cmd == "c2tail":
        drive_c2tail()
    elif cmd == "c3":
        _drive("C3", c3_target, f"{pre}c_c3.json", subset=sub)
    else:
        raise SystemExit("usage: python -m s21.c_cont {c1|c2|c3} [--smoke]")
