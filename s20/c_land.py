"""s20/c_land.py -- Q2: WHAT MAKES THE AMBER LANDSCAPE DIFFERENT?  THE PHYSICS HALF.

Pre-registration: `s20/PREREG_C.md` section 3.  Falsifiers F-C2 and F-C2b registered there,
first.  Workstream B owns the optimizer/CVaR half; this module owns the energetic
characterisation and nothing else.

THE SPACE.  CONTINUOUS torsion space `theta = (phi, psi) in R^{2n}`, radians.  No lattice, no
binary encoding (BRIEF section 5).  Both potentials are evaluated through the SAME ideal-geometry
rebuild `core.geometry.build_backbone`, so the parameterisation, the candidate set, the
initialisation, the optimiser and the evaluation budget are held identical and the ONLY thing
that changes between arms is `H_Legacy <-> H_AMBER`.

AN EXACT FACT ABOUT THE PARAMETERISATION, ESTABLISHED BEFORE ANY SPECTRUM IS READ.
`core.geometry.build_backbone`'s own docstring: *"phi[0] is never read"*.  So coordinate 0 of
theta is EXACTLY inert -- its gradient and its Hessian row are identically zero for BOTH
potentials -- and it would contribute a spurious exact null mode to every "near-zero mode count"
and an infinite condition number to every arm.  It is verified to be zero (`inert_max`) and then
REMOVED; the active dimension is `2n - 1`.  This is an implementation consequence, not a
discovery, and is labelled EXACT.

CROSS-POTENTIAL COMPARISONS ARE SCALE-INVARIANT ONLY (falsifier F-C2b).  AMBER is kcal/mol and
Legacy is an arbitrary weighted score; a raw gradient norm or eigenvalue comparison between them
is a unit error, not a result.  Only condition number, negative-curvature fraction, relative
spectral gap, participation ratio, anisotropy, signed spectral skew, relative barrier and the
COSINE between the two gradients are compared across potentials.  Everything with units is
reported WITHIN a potential.

    python -m s20.c_land --gate     # GC20b: the FD gradient, the h-plateau, and the inert coord
    python -m s20.c_land --smoke    # 2 targets
    python -m s20.c_land            # the declared 30-target subset
"""
from __future__ import annotations

import os
import sys
import json
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

from s12 import instrument as I                     # noqa: E402
from s14.avgspace import top75_windows              # noqa: E402
from s15 import seed as SD                          # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s18 import phys_lib as PL                      # noqa: E402

#: THE FINITE-DIFFERENCE STEPS, IN RADIANS, MEASURED ON THE PLATEAU AND THEN FROZEN.
#:
#: FIRST DEVIATION FROM `PREREG_C.md` section 3, and it is a correction of my own instrument that
#: my own gate forced.  The pre-registration declared ONE step and one plateau check, on the
#: GRADIENT.  That is not sufficient: a second difference divides the energy noise by h^2, and
#: `core.amber`'s single point is not double-clean.  MEASURED (`_DIAG_hessian_step.log`):
#:
#:      AMBER  E = 370.283 kcal/mol, granularity ~1e-6 .. 1e-4 kcal/mol over a 1e-7 rad move
#:      Legacy E = -16.006          granularity ~6e-7
#:
#: and the AMBER Hessian's NEGATIVE-CURVATURE FRACTION at the medoid start of 1A13 runs
#:
#:      h        0.03   0.02   0.01   0.006  0.003  0.002  0.001  0.0003  0.0001
#:      neg_frac 0.296  0.259  0.259  0.259  0.407  0.556  0.593  0.704   0.778
#:
#: **At the pre-registered h = 1e-4 the AMBER Hessian is round-off noise and would have been
#: published as "78% negative curvature".**  Legacy is flat to 5.6e-07 across the whole decade.
#: The plateau is h in [6e-3, 2e-2]; the step is set to its centre, THE SAME for both potentials
#: so the cross-potential comparison stays matched, and GC20b now reports the plateau table.
H_HESS = 1e-2
#: the gradient tolerates a smaller step (round-off enters as eps/h, not eps/h^2) and the
#: chain-rule agreement is >0.9999 over the whole decade; 1e-3 is on that plateau.
H_GRAD = 1e-3
H_FD = H_GRAD          # default for `grad_fd`/`jac`; the Hessian passes H_HESS explicitly
#: relative threshold below which an eigenvalue counts as a near-zero mode.  Relative, so the
#: count is scale-free and comparable across potentials (F-C2b).
NEARZERO_REL = 1e-4
#: two minima are DISTINCT if their torus distance exceeds this, in radians per coordinate.
DISTINCT_TOL = 0.10
#: identical for both potentials -- the "same evaluation budget" the brief requires.
MAXITER = 100
N_STARTS = 5
N_TARGETS = 30
#: straight-line interpolation points for the barrier estimate between two minima.
N_INTERP = 21

METRICS = ("neg_frac", "nearzero_frac", "cond_med", "gap_rel", "part_ratio",
           "aniso", "spec_skew")


# ==========================================================================
# the two potentials, as functions of theta, through the SAME rebuild
# ==========================================================================
class Pot:
    """Both potentials over theta, with a batched evaluator each."""

    def __init__(self, seq, want_amber=True):
        from core import energy as et
        from core import geometry as geo
        import torsion_lib2 as tl2
        self.seq = seq
        self.n = len(seq)
        self.geo = geo
        self.et = et
        self.w = EL.legacy_weight_vector()
        self.box = None
        self.n_eval = {"legacy": 0, "amber": 0}
        if want_amber:
            from s17 import phys_lib as P17
            tab = tl2.library_for(seq, 4, seq)
            rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
            self.box = P17.ConstrainedBox(seq, rep)

    def close(self):
        if self.box is not None:
            self.box.close()
            self.box = None

    def split(self, th):
        th = np.asarray(th, float)
        return th[:self.n], th[self.n:]

    def legacy(self, TH):
        """Batched.  `TH` is (B, 2n) or (2n,)."""
        TH = np.atleast_2d(np.asarray(TH, float))
        PHI, PSI = TH[:, :self.n], TH[:, self.n:]
        c = self.geo.build_backbone_batch(PHI, PSI)
        comp = self.et.components_batch(self.seq, c, PHI, PSI)
        M = np.column_stack([np.asarray(comp[t], float) for t in EL.LEG_TERMS])
        self.n_eval["legacy"] += len(TH)
        return M @ self.w

    def amber(self, TH):
        TH = np.atleast_2d(np.asarray(TH, float))
        out = np.empty(len(TH))
        for b in range(len(TH)):
            c = self.geo.build_backbone(TH[b, :self.n], TH[b, self.n:])
            out[b] = self.box.energy_point({k: np.asarray(v, float) for k, v in c.items()})[
                "energy"]
        self.n_eval["amber"] += len(TH)
        return out

    def E(self, which, TH):
        return self.legacy(TH) if which == "legacy" else self.amber(TH)

    # ---- the CHANNEL DECOMPOSITION of the Legacy gradient -------------------
    #: FOUND BY GATE GC20b, NOT ANTICIPATED.  `core.energy.components_batch(seq, c, PHI, PSI)`
    #: reads the torsions BOTH through the rebuilt coordinates `c` AND directly as arguments.
    #: AMBER reads only coordinates.  So E_Legacy = f(seq, x(theta), theta) has an
    #: EXPLICIT channel that AMBER structurally cannot have, and the two are separable:
    #:
    #:      g_coord[k] = d/dtheta_k f(seq, x(theta), theta_0)     (rebuild only)
    #:      g_expl[k]  = d/dtheta_k f(seq, x(theta_0), theta)     (arguments only)
    #:      g_total    = g_coord + g_expl                          (to first order)
    #:
    #: This is what makes `phi[0]` -- a coordinate that `build_backbone` never reads -- carry a
    #: NON-ZERO Legacy gradient and an EXACTLY ZERO AMBER gradient.
    def legacy_channels(self, th, h=None):
        h = H_FD if h is None else h
        th = np.asarray(th, float)
        p = len(th)
        PHI0, PSI0 = th[None, :self.n], th[None, self.n:]
        c0 = self.geo.build_backbone_batch(PHI0, PSI0)

        def _tot(TH):
            return self.legacy(TH)

        def _coord(TH):
            PHI, PSI = TH[:, :self.n], TH[:, self.n:]
            c = self.geo.build_backbone_batch(PHI, PSI)
            comp = self.et.components_batch(self.seq, c,
                                            np.repeat(PHI0, len(TH), 0),
                                            np.repeat(PSI0, len(TH), 0))
            return np.column_stack([np.asarray(comp[t], float)
                                    for t in EL.LEG_TERMS]) @ self.w

        def _expl(TH):
            cc = {k: np.repeat(np.asarray(v, float), len(TH), 0) for k, v in c0.items()}
            comp = self.et.components_batch(self.seq, cc, TH[:, :self.n], TH[:, self.n:])
            return np.column_stack([np.asarray(comp[t], float)
                                    for t in EL.LEG_TERMS]) @ self.w

        TH = np.repeat(th[None], 2 * p, axis=0)
        for k in range(p):
            TH[2 * k, k] += h
            TH[2 * k + 1, k] -= h
        out = {}
        for nm, fn in (("total", _tot), ("coord", _coord), ("expl", _expl)):
            v = fn(TH)
            out[nm] = (v[0::2] - v[1::2]) / (2.0 * h)
        return out

    # ---- the CHAIN-RULE gradient, used for the minimisation and as GC20b's reference ----
    def jac(self, th, h=H_FD):
        """dx/dtheta of the ASSEMBLED all-atom positions, by forward difference.

        Forward, not central: the map theta -> x is the deterministic rebuild and is smooth, so
        the O(h) truncation is far below the tolerance GC20b checks, and it halves the cost.
        """
        H = self.box.H
        th = np.asarray(th, float)
        x0 = self._x(th)
        J = np.empty((len(th), x0.size))
        for k in range(len(th)):
            tp = th.copy(); tp[k] += h
            J[k] = (self._x(tp) - x0).ravel() / h
        return J, x0

    def _x(self, th):
        c = self.geo.build_backbone(th[:self.n], th[self.n:])
        H = self.box.H
        return np.asarray(H._assemble(H._heavy_positions(
            {k: np.asarray(v, float) for k, v in c.items()}, chi1=None)), float)

    def amber_grad_chain(self, th):
        """-J^T f, with f the analytic OpenMM force.  One force call plus 2n cheap rebuilds."""
        from openmm import unit
        H = self.box.H
        if self.box._soft is None:
            self.box._make_soft()
        ctx = self.box._soft[0]
        J, x0 = self.jac(th)
        ctx.setPositions(x0 * unit.nanometer)
        for p in ("kh_ca", "kf_ca", "eps_ca", "kh_bb"):
            ctx.setParameter(p, 0.0)
        st = ctx.getState(getEnergy=True, getForces=True)
        e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        f = np.asarray(st.getForces(asNumpy=True).value_in_unit(
            unit.kilocalorie_per_mole / unit.nanometer), float)
        self.n_eval["amber"] += 1
        return float(e), -(J @ f.ravel())

    # ---- finite-difference gradient (the science instrument) ----
    def grad_fd(self, which, th, h=H_FD):
        th = np.asarray(th, float)
        p = len(th)
        TH = np.repeat(th[None], 2 * p, axis=0)
        for k in range(p):
            TH[2 * k, k] += h
            TH[2 * k + 1, k] -= h
        v = self.E(which, TH)
        return (v[0::2] - v[1::2]) / (2.0 * h)

    def grad(self, which, th):
        """The gradient used inside the minimiser: chain rule for AMBER (cheap), FD for
        Legacy (already cheap because Legacy is batched)."""
        if which == "amber":
            return self.amber_grad_chain(th)
        e = float(self.legacy(th)[0])
        return e, self.grad_fd("legacy", th)

    # ---- Hessian, by central difference of the ENERGY (no Jacobian anywhere) ----
    def hess_fd(self, which, th, h=H_HESS, active=None):
        th = np.asarray(th, float)
        idx = np.arange(len(th)) if active is None else np.asarray(active, int)
        p = len(idx)
        e0 = float(self.E(which, th)[0])
        #: diagonal
        TH = []
        for a in range(p):
            for s in (+1, -1):
                t = th.copy(); t[idx[a]] += s * h; TH.append(t)
        #: off-diagonal (++, +-, -+, --)
        pairs = [(a, b) for a in range(p) for b in range(a + 1, p)]
        for a, b in pairs:
            for sa, sb in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                t = th.copy(); t[idx[a]] += sa * h; t[idx[b]] += sb * h; TH.append(t)
        v = self.E(which, np.asarray(TH, float))
        Hm = np.zeros((p, p))
        for a in range(p):
            Hm[a, a] = (v[2 * a] - 2.0 * e0 + v[2 * a + 1]) / (h * h)
        off = 2 * p
        for q, (a, b) in enumerate(pairs):
            pp, pm, mp, mm = v[off + 4 * q: off + 4 * q + 4]
            Hm[a, b] = Hm[b, a] = (pp - pm - mp + mm) / (4.0 * h * h)
        return Hm, e0


# ==========================================================================
# scale-invariant spectrum metrics
# ==========================================================================
def spectrum_metrics(Hm):
    lam = np.linalg.eigvalsh(0.5 * (Hm + Hm.T))
    p = len(lam)
    a = np.abs(lam)
    amax = a.max() if p else np.nan
    if not np.isfinite(amax) or amax <= 0:
        return {m: float("nan") for m in METRICS} | {"lam": lam.tolist()}
    srt = np.sort(lam)
    out = {
        "neg_frac": float((lam < -NEARZERO_REL * amax).mean()),
        "nearzero_frac": float((a < NEARZERO_REL * amax).mean()),
        #: condition number over |lambda| with the MEDIAN as the denominator, so a single
        #: near-zero mode cannot send it to infinity.  Declared before the run.
        "cond_med": float(amax / max(np.median(a), 1e-30)),
        #: relative spectral gap at the bottom of the spectrum
        "gap_rel": float((srt[1] - srt[0]) / (srt[-1] - srt[0])) if p > 1 else float("nan"),
        #: participation ratio of the spectrum, in [1/p, 1] -- how many modes carry the curvature
        "part_ratio": float((a.sum() ** 2) / (p * (a ** 2).sum())),
        "aniso": float(amax / a.mean()),
        #: signed skew: +1 = purely convex, -1 = purely concave, 0 = balanced saddle
        "spec_skew": float(lam.sum() / a.sum()),
    }
    out["lam_max_abs"] = float(amax)          # HAS UNITS -- never compared across potentials
    out["lam"] = lam.tolist()
    return out


def torus_d(a, b):
    """per-coordinate RMS torus distance, radians."""
    d = np.mod(np.asarray(a, float) - np.asarray(b, float) + np.pi, 2 * np.pi) - np.pi
    return float(np.sqrt((d ** 2).mean()))


def minimise(P, which, th0, maxiter=MAXITER):
    from scipy.optimize import minimize as smin
    calls = {"n": 0}

    def fun(t):
        calls["n"] += 1
        e, g = P.grad(which, t)
        return e, g

    r = smin(fun, np.asarray(th0, float), jac=True, method="L-BFGS-B",
             options={"maxiter": int(maxiter), "maxfun": int(maxiter) * 2})
    return np.asarray(r.x, float), float(r.fun), int(calls["n"]), bool(r.success)


# ==========================================================================
def target_row(t, verbose=False):
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    P = Pot(seq)
    rec = {"pdb": pdb, "n": n, "fold": fold}

    #: STARTS -- identical for both potentials.  The pool medoid plus deterministic pool draws.
    Pm = I.pairwise_rmsd(W)
    med = int(I.medoid(Pm))
    rng = SD.stable_rng(pdb, "s20C_land")
    others = [int(x) for x in rng.permutation(len(W))[:N_STARTS] if int(x) != med]
    order = [med] + others[:N_STARTS - 1]
    starts = [np.concatenate([PHI[b], PSI[b]]) for b in order]

    #: the EXACT inert coordinate -- verified, then removed.  (phi[0] is never read.)
    th0 = starts[0]
    g_l = P.grad_fd("legacy", th0)
    g_a = P.grad_fd("amber", th0)
    rec["inert_max"] = {"legacy": float(abs(g_l[0])), "amber": float(abs(g_a[0]))}
    active = np.arange(1, 2 * n)

    #: THE CHANNEL DECOMPOSITION OF THE LEGACY GRADIENT (found by GC20b, not anticipated).
    #: How much of what Legacy "sees" is invisible to any structure-only potential?
    ch = P.legacy_channels(th0, h=H_GRAD)
    nt = max(np.linalg.norm(ch["total"]), 1e-30)
    rec["legacy_channel"] = {
        "additivity_rel": float(np.linalg.norm(ch["coord"] + ch["expl"] - ch["total"]) / nt),
        "expl_frac": float(np.linalg.norm(ch["expl"]) / nt),
        "coord_frac": float(np.linalg.norm(ch["coord"]) / nt),
        "cos_coord_expl": float(ch["coord"] @ ch["expl"]
                                / max(np.linalg.norm(ch["coord"])
                                      * np.linalg.norm(ch["expl"]), 1e-30)),
        #: and the cross-potential one: does the coordinate channel of Legacy agree with AMBER
        #: better than the total does?  (scale-free -- a cosine)
        "cos_coord_amber": float(ch["coord"][active] @ g_a[active]
                                 / max(np.linalg.norm(ch["coord"][active])
                                       * np.linalg.norm(g_a[active]), 1e-30)),
        "cos_total_amber": float(ch["total"][active] @ g_a[active]
                                 / max(np.linalg.norm(ch["total"][active])
                                       * np.linalg.norm(g_a[active]), 1e-30)),
    }

    #: LOCAL CURVATURE at the medoid start, on the active subspace.
    rec["start"] = {}
    for which, g in (("legacy", g_l), ("amber", g_a)):
        Hm, e0 = P.hess_fd(which, th0, active=active)
        m = spectrum_metrics(Hm)
        m["grad_norm"] = float(np.linalg.norm(g[active]))       # HAS UNITS
        m["E0"] = float(e0)                                     # HAS UNITS
        m.pop("lam")
        rec["start"][which] = m
    #: SCALE-FREE cross-potential quantity: do the two potentials pull the same way?
    ga, gl = g_a[active], g_l[active]
    den = np.linalg.norm(ga) * np.linalg.norm(gl)
    rec["grad_cos"] = float((ga @ gl) / den) if den > 0 else float("nan")

    #: MULTI-START MINIMISATION.  Identical starts, optimiser and budget; only H changes.
    rec["multi"] = {}
    for which in ("legacy", "amber"):
        mins, ends, drops, moved, rs, re_, ok, depths = [], [], [], [], [], [], [], []
        for th_s in starts:
            e_s = float(P.E(which, th_s)[0])
            th_e, e_e, nc, suc = minimise(P, which, th_s)
            mins.append(th_e); ends.append(e_e)
            depths.append(e_s - e_e)
            drops.append((e_s - e_e) / max(abs(e_s), 1e-9))
            moved.append(torus_d(th_e, th_s))
            rs.append(float(I.ca_rmsd(I.build_ca(th_s[:n], th_s[n:]), nat)))
            re_.append(float(I.ca_rmsd(I.build_ca(th_e[:n], th_e[n:]), nat)))
            ok.append(suc)
        #: distinct minima, and straight-line torus barriers between them
        dist = []
        for a in range(len(mins)):
            if all(torus_d(mins[a], mins[b]) > DISTINCT_TOL for b in dist):
                dist.append(a)
        bars, npair = [], 0
        for ii in range(len(dist)):
            for jj in range(ii + 1, len(dist)):
                a, b = mins[dist[ii]], mins[dist[jj]]
                d = np.mod(b - a + np.pi, 2 * np.pi) - np.pi
                TH = np.array([a + s * d for s in np.linspace(0, 1, N_INTERP)])
                v = P.E(which, TH)
                top = float(v.max())
                base = max(v[0], v[-1])
                #: SCALE-FREE normaliser, and it must be a DEPTH, not an energy LEVEL.  The
                #: first draft divided by |mean(E_min)|, which is an absolute energy: Legacy
                #: sits near -16 and AMBER near +370 in their own units, so that ratio would
                #: have compared two different things and would have been a units artefact --
                #: exactly falsifier F-C2b.  The denominator is the depth the SAME optimiser
                #: recovered from the SAME starts on the SAME potential, so the barrier is
                #: measured in units of "how far this potential's own minimisation fell".
                depth = max(float(np.mean(depths)), 1e-9)
                bars.append((top - base) / depth)
                npair += 1
        gend = [float(np.linalg.norm(P.grad_fd(which, m_)[active])) for m_ in mins[:2]]
        rec["multi"][which] = {
            "e_drop_frac": float(np.mean(drops)),
            "g_drop_frac": float(np.mean(gend) / max(rec["start"][which]["grad_norm"], 1e-30)),
            "n_distinct": float(len(dist)),
            "basin_width": float(np.mean(moved)),
            "theta_moved": float(np.mean(moved)),
            "rmsd_start": float(np.mean(rs)),
            "rmsd_end": float(np.mean(re_)),
            "rmsd_end_best": float(np.min(re_)),
            "n_pairs": float(npair),
            "barrier_rel": float(np.mean(bars)) if bars else float("nan"),
            "frac_barrier": float(np.mean([b > 1e-6 for b in bars])) if bars else float("nan"),
            "grad_norm": float(rec["start"][which]["grad_norm"]),
            "success_frac": float(np.mean(ok)),
        }
    rec["n_eval"] = dict(P.n_eval)
    P.close()
    return rec


# ==========================================================================
def gate_GC20b(n_targets=3, verbose=True, hess_plateau=True):
    """GC20b -- the instrument gate, and the record of what it refuted about my own design.

    SECOND DEVIATION FROM `PREREG_C.md` section 3, kept rather than repaired silently.  The
    pre-registered pass condition was *"phi[0] is EXACTLY inert for BOTH potentials"*.  It is
    FALSE for Legacy, and the reason is a real structural asymmetry between the two models
    rather than a bug:

        `core.energy.components_batch(seq, c, PHI, PSI)` reads the torsions BOTH through the
        rebuilt coordinates AND directly as arguments, so  E_Legacy = f(seq, x(theta), theta).
        `core.amber` reads coordinates ONLY, so  E_AMBER = f(x(theta)).

    `build_backbone` never reads `phi[0]`, so that coordinate moves NO atom -- and Legacy's
    energy changes anyway (|g| = 1.6e-03 / 1.4e-01 / 6.8e-02 on the three gate targets) while
    AMBER's is EXACTLY 0.0.  **Legacy has a direction in torsion space along which its energy
    varies and the structure does not.**  The pass condition is therefore restated as:

        (1) AMBER's `phi[0]` gradient is EXACTLY zero          -- EXACT, an implementation fact
        (2) Legacy's is NOT, and its size is MEASURED and reported, not gated on
        (3) the FD gradient reproduces the chain-rule gradient -J^T f  (cos > 0.999)
        (4) BOTH steps sit on their measured plateau, with the plateau table printed

    The original condition is left in this docstring so no reader has to take my word that the
    goalpost moved for a reason.  Per s19 Z6 the number of times each tolerance FIRED is
    reported; a pass with zero firings certifies only the chosen step, not the tolerance.
    """
    rows, fired_g, fired_h = [], 0, 0
    for t in I.targets()[:n_targets]:
        seq, n = t["seq"], int(t["n"])
        W, PHI, PSI, u = top75_windows(t["pdb"])
        th = np.concatenate([np.asarray(PHI[0], float), np.asarray(PSI[0], float)])
        P = Pot(seq)
        _e, gc = P.amber_grad_chain(th)
        r = {"pdb": t["pdb"], "n": n, "plateau_grad": {}, "plateau_hess": {}}
        for h in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 1e-5):
            gf = P.grad_fd("amber", th, h=h)
            a, b = gf[1:], gc[1:]
            den = np.linalg.norm(a) * np.linalg.norm(b)
            r["plateau_grad"][f"{h:g}"] = {
                "cos": float((a @ b) / den) if den > 0 else float("nan"),
                "rel_l2": float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-30))}
        gl = P.grad_fd("legacy", th, h=H_GRAD)
        ga = P.grad_fd("amber", th, h=H_GRAD)
        ch = P.legacy_channels(th, h=H_GRAD)
        act = np.arange(1, 2 * n)
        r["inert_phi0"] = {"legacy": float(abs(gl[0])), "amber": float(abs(ga[0])),
                           "chain": float(abs(gc[0]))}
        #: how much of Legacy's gradient is the coordinate-INVISIBLE channel
        r["legacy_channel"] = {
            "additivity_rel": float(np.linalg.norm(ch["coord"] + ch["expl"] - ch["total"])
                                    / max(np.linalg.norm(ch["total"]), 1e-30)),
            "expl_frac": float(np.linalg.norm(ch["expl"]) / max(np.linalg.norm(ch["total"]),
                                                                1e-30)),
            "coord_frac": float(np.linalg.norm(ch["coord"]) / max(np.linalg.norm(ch["total"]),
                                                                  1e-30)),
            "cos_coord_expl": float(ch["coord"] @ ch["expl"]
                                    / max(np.linalg.norm(ch["coord"])
                                          * np.linalg.norm(ch["expl"]), 1e-30))}
        #: the HESSIAN plateau, on the metric that is most noise-sensitive (neg_frac) and on the
        #: eigenvalue vector, for both potentials.
        if hess_plateau:
            for which in ("amber", "legacy"):
                prev = None
                for h in (2e-2, 1e-2, 6e-3, 3e-3, 1e-3):
                    Hm, _e0 = P.hess_fd(which, th, h=h, active=act)
                    lam = np.linalg.eigvalsh(0.5 * (Hm + Hm.T))
                    m = spectrum_metrics(Hm)
                    rel = (float(np.linalg.norm(lam - prev) / np.linalg.norm(prev))
                           if prev is not None else float("nan"))
                    prev = lam
                    r["plateau_hess"][f"{which}@{h:g}"] = {
                        "neg_frac": m["neg_frac"], "cond_med": m["cond_med"],
                        "part_ratio": m["part_ratio"], "spec_skew": m["spec_skew"],
                        "rel_d_lam_vs_prev_h": rel}
        okg = r["plateau_grad"][f"{H_GRAD:g}"]["cos"] > 0.999
        oki = (r["inert_phi0"]["amber"] == 0.0 and r["inert_phi0"]["chain"] == 0.0)
        okh = True
        if hess_plateau:
            for which in ("amber", "legacy"):
                a1 = r["plateau_hess"][f"{which}@{H_HESS:g}"]["neg_frac"]
                a2 = r["plateau_hess"][f"{which}@{6e-3:g}"]["neg_frac"]
                if abs(a1 - a2) > 0.05:
                    okh = False
        if not okg:
            fired_g += 1
        if not okh:
            fired_h += 1
        r["passed"] = bool(okg and oki and okh)
        rows.append(r)
        P.close()
        if verbose:
            pg = r["plateau_grad"][f"{H_GRAD:g}"]
            print(f"GC20b {t['pdb']}: cos(FD,chain)@h={H_GRAD:g} {pg['cos']:.6f} "
                  f"(rel_l2 {pg['rel_l2']:.2e}) | |g_phi0| AMBER {r['inert_phi0']['amber']:.1e} "
                  f"(EXACT 0) Legacy {r['inert_phi0']['legacy']:.2e} (NOT 0) | "
                  f"Legacy expl-channel {r['legacy_channel']['expl_frac']:.3f} of |g| | "
                  f"hess neg_frac A {r['plateau_hess'].get(f'amber@{H_HESS:g}',{}).get('neg_frac',float('nan')):.3f} "
                  f"L {r['plateau_hess'].get(f'legacy@{H_HESS:g}',{}).get('neg_frac',float('nan')):.3f}"
                  f" -> {'PASS' if r['passed'] else 'FAIL'}", flush=True)
    out = {"rows": rows, "h_grad": H_GRAD, "h_hess": H_HESS,
           "n_grad_tolerance_firings": int(fired_g),
           "n_hess_plateau_firings": int(fired_h),
           "prereg_condition_as_written": "phi[0] EXACTLY inert for BOTH potentials",
           "prereg_condition_verdict": "FALSE for Legacy -- see the docstring; the condition was "
                                       "mis-specified and is restated, not deleted",
           "passed": bool(all(r["passed"] for r in rows))}
    if verbose:
        print(f"GC20b overall {'PASS' if out['passed'] else 'FAIL'};  gradient tolerance FIRED "
              f"{fired_g}/{n_targets}, Hessian plateau tolerance FIRED {fired_h}/{n_targets}. "
              f"A pass with zero firings certifies the chosen step, not the tolerance (s19 Z6).")
    return out


def subset(n=N_TARGETS):
    """The declared subset, drawn deterministically -- NOT the first n in list order."""
    tg = I.targets()
    rng = SD.stable_rng("s20C_land", "subset")
    idx = np.sort(rng.permutation(len(tg))[:n])
    return [tg[int(i)] for i in idx]


def run(n=N_TARGETS, out="c_land.json", verbose=True):
    tg = subset(n)
    path = os.path.join(RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    #: the gate is expensive (five Hessians per potential per target) and deterministic, so a
    #: PASSED persisted result is reused rather than recomputed.  A missing or failed one is
    #: re-run; the science run refuses to start without a pass.
    gp = os.path.join(RESULTS, "c_land_gate.json")
    gate = None
    if os.path.exists(gp):
        try:
            g = json.load(open(gp))
            if g.get("passed") and g.get("GC20b", {}).get("h_hess") == H_HESS \
                    and g.get("GC20b", {}).get("h_grad") == H_GRAD:
                gate = g["GC20b"]
                print("GC20b: reusing the persisted PASS "
                      f"(h_grad={H_GRAD:g}, h_hess={H_HESS:g})", flush=True)
        except Exception:
            gate = None
    if gate is None:
        gate = gate_GC20b(verbose=verbose)
        with open(gp, "w") as fh:
            json.dump({"GC20b": gate, "complete": True, "passed": gate["passed"]}, fh, indent=1)
    if not gate["passed"]:
        raise SystemExit("GC20b FAILED -- the Hessian block does not run (PREREG_C.md 3)")
    cfg = {"H_FD": H_FD, "NEARZERO_REL": NEARZERO_REL, "DISTINCT_TOL": DISTINCT_TOL,
           "MAXITER": MAXITER, "N_STARTS": N_STARTS, "N_INTERP": N_INTERP,
           "METRICS": list(METRICS), "subset": [t["pdb"] for t in tg]}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.55)
        rows.append(target_row(t))
        if verbose:
            r = rows[-1]
            print(f"  {len(rows)}/{len(tg)} {t['pdb']} n={t['n']} "
                  f"gradcos={r['grad_cos']:+.3f} "
                  f"negfrac L={r['start']['legacy']['neg_frac']:.3f} "
                  f"A={r['start']['amber']['neg_frac']:.3f} "
                  f"rmsd L {r['multi']['legacy']['rmsd_start']:.3f}->"
                  f"{r['multi']['legacy']['rmsd_end']:.3f} "
                  f"A {r['multi']['amber']['rmsd_start']:.3f}->"
                  f"{r['multi']['amber']['rmsd_end']:.3f} ({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg), gate)
    _write(out, rows, cfg, len(tg), gate)
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)


def _write(out, rows, cfg, n_expected, gate):
    obj = {"rows": rows, "config": cfg, "gate": gate, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected))}
    p = os.path.join(RESULTS, out)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--n", type=int, default=N_TARGETS)
    a = ap.parse_args()
    if a.gate:
        g = gate_GC20b()
        with open(os.path.join(RESULTS, "c_land_gate.json"), "w") as fh:
            json.dump({"GC20b": g, "complete": True, "passed": g["passed"]}, fh, indent=1)
        sys.exit(0 if g["passed"] else 1)
    if a.smoke:
        run(n=2, out="_SMOKE_c_land.json")
    else:
        run(n=a.n)
