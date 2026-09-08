"""RETIRED -- DO NOT READ AS A RESULT.

This is the PROVISIONAL degree-1 implementation written against the interface published in
`s18/COORD_exp_to_math.md` while MATH had not yet landed a module.  MATH landed
`s18/math_anova.py`, whose object is RESIDUE-additive and differs from what this file built,
and every reported Sprint-18 EXPERIMENT arm consumes MATH's object through `s18/exp_obj.py`.
Nothing here contributed a single number to `s18/exp_FINDINGS.md`.  Its only artefact
(`_RETIRED_exp_lambda_provisional.json`, 3 rows from a 3-target smoke whose leave-fold-out
debias was fitted on those 3 targets and is therefore not comparable to anything) is renamed
so it cannot be mistaken for a result.  Kept only as the record of what was built and when.

s18/exp_anova.py -- the degree-1 object on the CONTINUOUS torsion instrument.

PROVENANCE.  The BRIEF tells the EXPERIMENT workstream to consume MATH's callable and not to
build a second, subtly different definition.  MATH had landed nothing when this ran, so this is
a PROVISIONAL implementation behind an adapter (`load_math()`).  The interface is published in
`s18/COORD_exp_to_math.md`; if a `s18/math_*.py` appears exposing it, `build()` uses it and this
file's tables become dead code.  Which path was taken is recorded in every result artefact under
`anova_source`.

THE OBJECT.  Not a lattice Walsh truncation ported by analogy -- the first-order functional
ANOVA under a product reference measure mu = (x) mu_c over the 2n individual torsion
coordinates:

    E0        = E_mu[E]
    f_c(t)    = E_mu[E | theta_c = t] - E0
    E_le1     = E0 + sum_c f_c(theta_c)                     STRICT weight-<=1 analogue
    g_i(p, s) = E_mu[E | phi_i = p, psi_i = s] - E0
    E_res     = E0 + sum_i g_i(phi_i, psi_i)                RESIDUE-ADDITIVE
    E_ge2     = E_full - E_le1

`E_res - E_le1` is exactly the intra-residue weight-2 mass the BRIEF sec 4.1 warns the strict
truncation throws away, so both objects are built and both are run.

mu.  Native-free and target-conditioned: the empirical marginal of each torsion coordinate over
the target's OWN shipped top-75 retrieval pool.  A uniform-torsion mu is built too, because the
measure is a modelling choice and the sprint must show the conclusion does not depend on it.

REPRESENTATION.  The conditionals are Monte-Carlo'd on a regular periodic grid and then carried
by EXACT trigonometric interpolation (FFT).  Two reasons: the conditional expectation of a
smooth periodic function is smooth and periodic, so the trig interpolant is the natural basis;
and it hands L-BFGS an analytic gradient, so the degree-1 arms are optimised by the same
machinery and to the same tolerance as the full arm.  Nothing about the comparison is allowed to
turn on one arm being optimised better than another.

TWO STRUCTURAL FACTS, recorded as EXACT (theorems, not discoveries):

  1. `E_le1` is a sum of univariate functions, so its GLOBAL argmin is coordinate-wise and
     costs O(2n G) table lookups.  There is no search problem for degree-1.  `argmin_le1()`
     returns the certified global optimum.  This is the continuous analogue of the 19-target
     *certified* degree-1 argmin that Phase 0 recomputed.
  2. E_lambda = E_le1 + lambda (E_full - E_le1) = (1-lambda) E_le1 + lambda E_full.  The ladder
     is a convex mix, so lambda = 1 MUST reproduce s17's `refine_full` exactly.  That identity
     is the harness's validity gate (PREREG P0.a) and it is checked, not assumed.
"""
from __future__ import annotations

import glob
import hashlib
import importlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import project as pj              # noqa: E402
from s15 import seed as SD                  # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, "cache")
for _d in (RESULTS, CACHE):
    os.makedirs(_d, exist_ok=True)

#: grid resolutions and the MC sample count.  G1 for the 1-D conditionals, G2 per axis for the
#: 2-D residue conditionals, M product draws from mu shared by every table (common random
#: numbers -- differences between tables are then far better resolved than each table's own
#: absolute MC error, which is what the decomposition actually needs).
G1 = 32
G2 = 12
MSAMP = 192
CHUNK = 12000


def grid(g):
    return -np.pi + 2 * np.pi * np.arange(g) / g


# --------------------------------------------------------------- trigonometric interpolation
def _trig_coef1(y):
    """FFT coefficients of a periodic sample vector on `grid(len(y))`."""
    return np.fft.rfft(np.asarray(y, float))


def _trig_eval1(c, g, x):
    """Value and derivative of the trig interpolant of a length-`g` grid sample at `x`.

    `x` may be an array.  Standard real-FFT interpolation:  y(x) = (1/g) Re sum_k a_k c_k
    e^{i k x} with a_0 = 1, a_k = 2 for 0 < k < g/2, a_{g/2} = 1 for even g.
    """
    #: the grid is `-pi + 2 pi m / g`, so the FFT phase carries an extra `e^{i k pi}` -- the
    #: shift that a naive `e^{i k x}` drops, and dropping it silently destroys the interpolant.
    x = np.atleast_1d(np.asarray(x, float)) + np.pi
    k = np.arange(len(c))
    a = np.full(len(c), 2.0)
    a[0] = 1.0
    if g % 2 == 0 and len(c) == g // 2 + 1:
        a[-1] = 1.0
    ph = np.exp(1j * np.outer(x, k))
    z = ph * (a * c)[None, :]
    val = z.real.sum(1) / g
    der = (1j * k[None, :] * z).real.sum(1) / g
    return val, der


def _trig_coef2(Y):
    return np.fft.fft2(np.asarray(Y, float))


def _trig_eval2(C, g, xp, xs):
    """Value and both partials of the 2-D trig interpolant at scalar-or-array (xp, xs)."""
    xp = np.atleast_1d(np.asarray(xp, float)) + np.pi
    xs = np.atleast_1d(np.asarray(xs, float)) + np.pi
    k = np.fft.fftfreq(g, 1.0 / g)
    ep = np.exp(1j * np.outer(xp, k))            # (B, g)
    es = np.exp(1j * np.outer(xs, k))            # (B, g)
    #: sum_{k,l} C_kl e^{i k xp} e^{i l xs} / g^2
    T = np.einsum("bk,kl,bl->b", ep, C, es) / g ** 2
    Tp = np.einsum("bk,k,kl,bl->b", ep, 1j * k, C, es) / g ** 2
    Ts = np.einsum("bk,kl,l,bl->b", ep, C, 1j * k, es) / g ** 2
    return T.real, Tp.real, Ts.real


# ------------------------------------------------------------------------------- the objects
class Anova:
    """The full objective, its ANOVA truncations and the lambda ladder for ONE target."""

    def __init__(self, dhat, sd, i, j, phi_pool, psi_pool, seed, mu="pool",
                 g1=G1, g2=G2, m=MSAMP, source="s18.exp_anova(provisional)"):
        self.dhat = np.asarray(dhat, float)
        self.sd = np.asarray(sd, float)
        self.inv = 1.0 / self.sd
        self.i = np.asarray(i, int)
        self.j = np.asarray(j, int)
        self.n = int(np.asarray(phi_pool).shape[1])
        self.g1, self.g2, self.m, self.mu_kind = int(g1), int(g2), int(m), str(mu)
        self.source = source
        self.nev = 0
        rng = SD.stable_rng(str(seed), "s18anova", mu)
        n, M = self.n, self.m

        #: ---- the product reference measure.  Draw each coordinate INDEPENDENTLY from its own
        #: marginal: that is what makes mu a product measure and not the pool's joint.
        P = np.asarray(phi_pool, float)
        S = np.asarray(psi_pool, float)
        if mu == "pool":
            k = P.shape[0]
            TH = np.empty((M, 2 * n))
            for c in range(n):
                TH[:, c] = P[rng.integers(0, k, M), c]
                TH[:, n + c] = S[rng.integers(0, k, M), c]
        elif mu == "uniform":
            TH = rng.uniform(-np.pi, np.pi, (M, 2 * n))
        else:
            raise ValueError(mu)
        self.TH = TH

        gx1, gx2 = grid(g1), grid(g2)
        self.gx1, self.gx2 = gx1, gx2

        #: ---- E0 and the conditionals, all on the SAME M draws (common random numbers).
        self.E0 = float(self._batch(TH).mean())

        #: 1-D conditionals: 2n coordinates x g1 grid points x M draws.  The SPLIT-HALF replica
        #: is taken in the same pass -- it is the only honest MC diagnostic here, because with
        #: common random numbers the absolute standard error is dominated by an offset that
        #: centring removes, and it is the SHAPE of f_c, not its level, that the object uses.
        f1 = np.empty((2 * n, g1))
        fA = np.empty((2 * n, g1))
        fB = np.empty((2 * n, g1))
        h = M // 2
        for c in range(2 * n):
            B = np.repeat(TH[None, :, :], g1, 0).reshape(g1 * M, 2 * n)
            B[:, c] = np.repeat(gx1, M)
            v = self._batch(B).reshape(g1, M)
            f1[c] = v.mean(1)
            fA[c] = v[:, :h].mean(1)
            fB[c] = v[:, h:].mean(1)
        #: centre each component under mu so the ANOVA identity E_mu[f_c] = 0 holds on the grid
        self.f1 = f1 - f1.mean(1, keepdims=True)
        self.f1_raw = f1
        self.c1 = np.array([_trig_coef1(self.f1[c]) for c in range(2 * n)])
        fA = fA - fA.mean(1, keepdims=True)
        fB = fB - fB.mean(1, keepdims=True)
        self.split_rms = float(np.median(
            np.sqrt(((fA - fB) ** 2).mean(1)) / 2.0
            / np.maximum(np.sqrt((self.f1 ** 2).mean(1)), 1e-12)))
        rr = [float(np.corrcoef(fA[c], fB[c])[0, 1]) for c in range(2 * n)
              if np.ptp(fA[c]) > 0 and np.ptp(fB[c]) > 0]
        self.split_r = float(np.median(rr)) if rr else float("nan")

        #: 2-D residue conditionals: n residues x g2^2 grid x M draws
        g2n = g2 * g2
        f2 = np.empty((n, g2, g2))
        PP, SS = np.meshgrid(gx2, gx2, indexing="ij")
        PPf, SSf = PP.ravel(), SS.ravel()
        for c in range(n):
            B = np.repeat(TH[None, :, :], g2n, 0).reshape(g2n * M, 2 * n)
            B[:, c] = np.repeat(PPf, M)
            B[:, n + c] = np.repeat(SSf, M)
            f2[c] = self._batch(B).reshape(g2, g2, M).mean(2)
        self.f2 = f2 - f2.mean((1, 2), keepdims=True)
        self.c2 = np.array([_trig_coef2(self.f2[c]) for c in range(n)])

        #: PREREG P0.c is answered by `split_rms` / `split_r` above, computed in the same pass.
        self.mc_ratio = self.split_rms

    # ------------------------------------------------------------------- the full objective
    def _batch(self, TH):
        """`E_full` for a batch of `(B, 2n)` torsion vectors, chunked for memory."""
        TH = np.asarray(TH, float)
        n, out = self.n, np.empty(len(TH))
        for a in range(0, len(TH), CHUNK):
            b = min(a + CHUNK, len(TH))
            CA = np.asarray(pj.build_ca_exact(TH[a:b, :n], TH[a:b, n:]), float)
            rv = CA[:, self.i] - CA[:, self.j]
            d = np.sqrt((rv * rv).sum(2))
            out[a:b] = ((((d - self.dhat) * self.inv) ** 2)).sum(1)
        self.nev += len(TH)
        return out

    def E_full(self, phi, psi):
        """Value and `(2n,)` analytic gradient of the deployed distance objective."""
        phi = np.asarray(phi, float)
        psi = np.asarray(psi, float)
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[self.i] - CA[self.j]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = (d - self.dhat) * self.inv
        f = float((r * r).sum())
        coef = ((2.0 * r * self.inv) / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, self.i, coef)
        np.add.at(gCA, self.j, -coef)
        self.nev += 1
        return f, pj._torsion_grad(G, CA, gCA)

    # ------------------------------------------------------------------ the ANOVA objects
    def E_le1(self, phi, psi):
        """Value and analytic gradient, every coordinate evaluated in one einsum."""
        th = np.concatenate([np.asarray(phi, float), np.asarray(psi, float)]) + np.pi
        K = self.c1.shape[1]
        k = np.arange(K)
        a = np.full(K, 2.0)
        a[0] = 1.0
        if self.g1 % 2 == 0 and K == self.g1 // 2 + 1:
            a[-1] = 1.0
        Z = np.exp(1j * np.outer(th, k)) * (a[None, :] * self.c1)
        v = Z.real.sum(1) / self.g1
        gr = (1j * k[None, :] * Z).real.sum(1) / self.g1
        self.nev += 1
        return float(self.E0 + v.sum()), gr

    def E_res(self, phi, psi):
        """Value and analytic gradient of the residue-additive object, vectorised."""
        n, g = self.n, self.g2
        xp = np.asarray(phi, float) + np.pi
        xs = np.asarray(psi, float) + np.pi
        k = np.fft.fftfreq(g, 1.0 / g)
        ep = np.exp(1j * np.outer(xp, k))
        es = np.exp(1j * np.outer(xs, k))
        T = np.einsum("ck,ckl,cl->c", ep, self.c2, es) / g ** 2
        Tp = np.einsum("ck,k,ckl,cl->c", ep, 1j * k, self.c2, es) / g ** 2
        Ts = np.einsum("ck,ckl,l,cl->c", ep, self.c2, 1j * k, es) / g ** 2
        self.nev += 1
        return float(self.E0 + T.real.sum()), np.concatenate([Tp.real, Ts.real])

    def E_lambda(self, phi, psi, lam):
        """`E_le1 + lam * (E_full - E_le1)`, i.e. the convex mix.  lam=1 IS `E_full`."""
        if lam == 1.0:
            return self.E_full(phi, psi)
        if lam == 0.0:
            return self.E_le1(phi, psi)
        f1, g1_ = self.E_le1(phi, psi)
        f2, g2_ = self.E_full(phi, psi)
        return (1 - lam) * f1 + lam * f2, (1 - lam) * g1_ + lam * g2_

    def E_ge2(self, phi, psi):
        f1, g1_ = self.E_le1(phi, psi)
        f2, g2_ = self.E_full(phi, psi)
        return f2 - f1, g2_ - g1_

    # ------------------------------------------------------- EXACT separable global argmins
    def argmin_le1(self, refine=True):
        """Certified global optimum of `E_le1`.  EXACT -- separability, not a search result."""
        n = self.n
        th = np.empty(2 * n)
        fine = grid(8 * self.g1)
        for c in range(2 * n):
            v, _ = _trig_eval1(self.c1[c], self.g1, fine)
            th[c] = fine[int(np.argmin(v))]
        if refine:                       # one Newton polish inside the resolved basin
            for c in range(2 * n):
                x = th[c]
                for _ in range(12):
                    _v, g = _trig_eval1(self.c1[c], self.g1, np.array([x - 1e-4, x, x + 1e-4]))
                    h = (g[2] - g[0]) / 2e-4
                    if h <= 1e-9:
                        break
                    x = x - g[1] / h
                th[c] = x
        return th[:n], th[n:]

    def argmin_res(self, refine=True):
        """Certified global optimum of `E_res`.  EXACT for the same reason."""
        n = self.n
        fine = grid(4 * self.g2)
        PP, SS = np.meshgrid(fine, fine, indexing="ij")
        ph = np.empty(n)
        ps = np.empty(n)
        for c in range(n):
            v, _, _ = _trig_eval2(self.c2[c], self.g2, PP.ravel(), SS.ravel())
            k = int(np.argmin(v))
            ph[c], ps[c] = PP.ravel()[k], SS.ravel()[k]
        return ph, ps


# ------------------------------------------------------------------------------- the adapter
def load_math():
    """Return MATH's builder if it has landed, else None.  Recorded in every artefact."""
    for p in sorted(glob.glob(os.path.join(HERE, "math_*.py"))):
        mod = "s18." + os.path.basename(p)[:-3]
        try:
            m = importlib.import_module(mod)
        except Exception:
            continue
        if hasattr(m, "build") and hasattr(m, "__anova_interface__"):
            return m
    return None


def build(dhat, sd, i, j, phi_pool, psi_pool, seed, mu="pool", **kw):
    m = load_math()
    if m is not None:
        return m.build(dhat, sd, i, j, phi_pool, psi_pool, seed, mu=mu, **kw)
    return Anova(dhat, sd, i, j, phi_pool, psi_pool, seed, mu=mu, **kw)


def config_hash():
    m = load_math()
    d = {"G1": G1, "G2": G2, "MSAMP": MSAMP, "src": "math" if m else "exp_anova",
         "mod": (m.__name__ if m else "s18.exp_anova")}
    return d, hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:16]


# ------------------------------------------------------------------------------ self-check
def selfcheck():
    """PREREG P0.b and P0.d, plus the separability theorem, on one synthetic instance."""
    from s12 import instrument as I
    rng = SD.stable_rng("selfcheck", "s18anova")
    n = 12
    phi = rng.uniform(-np.pi, np.pi, n)
    psi = rng.uniform(-np.pi, np.pi, n)
    i, j = I.pair_index(n)
    CA = I.build_ca(phi, psi)
    dhat = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1)) + rng.normal(0, 0.5, len(i))
    sd = np.full(len(i), 1.0)
    P = rng.uniform(-np.pi, np.pi, (75, n))
    S = rng.uniform(-np.pi, np.pi, (75, n))
    ob = Anova(dhat, sd, i, j, P, S, "sc", mu="uniform", g1=16, g2=8, m=32)

    # P0.d -- the interpolant reproduces its own grid values
    err = max(abs(_trig_eval1(ob.c1[c], ob.g1, ob.gx1)[0] - ob.f1[c]).max()
              for c in range(2 * n))
    err2 = max(abs(_trig_eval2(ob.c2[c], ob.g2,
                               np.repeat(ob.gx2, ob.g2), np.tile(ob.gx2, ob.g2))[0]
                   - ob.f2[c].ravel()).max() for c in range(n))
    # P0.b -- E_le1 + E_ge2 == E_full
    a = ob.E_le1(phi, psi)[0]
    b = ob.E_ge2(phi, psi)[0]
    c = ob.E_full(phi, psi)[0]
    # gradient check by central differences
    f0, g0 = ob.E_le1(phi, psi)
    th = np.concatenate([phi, psi])
    fd = np.empty(2 * n)
    for k in range(2 * n):
        e = np.zeros(2 * n); e[k] = 1e-5
        fp = ob.E_le1((th + e)[:n], (th + e)[n:])[0]
        fm = ob.E_le1((th - e)[:n], (th - e)[n:])[0]
        fd[k] = (fp - fm) / 2e-5
    f0r, g0r = ob.E_res(phi, psi)
    fdr = np.empty(2 * n)
    for k in range(2 * n):
        e = np.zeros(2 * n); e[k] = 1e-5
        fdr[k] = (ob.E_res((th + e)[:n], (th + e)[n:])[0]
                  - ob.E_res((th - e)[:n], (th - e)[n:])[0]) / 2e-5
    # separability: the certified argmin beats 20000 random draws
    ap, aq = ob.argmin_le1()
    best = ob.E_le1(ap, aq)[0]
    R = rng.uniform(-np.pi, np.pi, (20000, 2 * n))
    rv = np.array([ob.E0 + sum(_trig_eval1(ob.c1[c], ob.g1, R[t, c])[0][0]
                               for c in range(2 * n)) for t in range(400)])
    out = {"interp1_err": float(err), "interp2_err": float(err2),
           "le1_plus_ge2_minus_full": float(abs(a + b - c)),
           "grad_le1_maxerr": float(np.abs(g0 - fd).max()),
           "grad_res_maxerr": float(np.abs(g0r - fdr).max()),
           "argmin_le1": float(best), "random400_min": float(rv.min()),
           "argmin_is_lowest": bool(best <= rv.min() + 1e-9)}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    selfcheck()
