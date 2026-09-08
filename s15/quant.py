"""SPRINT 15, RESTRAINT agent -- what does the BINNED distogram output cost, and can a
continuous density buy it back?

THE FINDING THAT CREATED THIS MODULE.  Two oracles that both know the true distance matrix
perfectly disagree by ~1.06 A:

    ORACLE_true_distances (continuous weighted least squares, s15/distgeo.py)   0.697
    ORACLE_ml_true        (ML against a 0.6 A Gaussian SAMPLED at the 17 bin      1.756
                           centres, s15/distml.py)

Same optimiser, same starts, same selection-by-objective rule, same geometry.  The task was
to attribute that gap to the binned output representation.

WHAT THIS MODULE FOUND FIRST, BEFORE ANY EXPERIMENT (see restraint_FINDINGS.md, R0).
`s15/distml.LogPTable` assumes the 17 bin centres are UNIFORMLY spaced and indexes with
`w = centres[1] - centres[0] = 0.75 A`.  They are not uniform: they run
[4.0 4.75 5.25 ... 7.75 8.5 9.5 10.5 11.75 13.25 15.0 17.5 21.0 25.0].
A pair at 8.0 A is therefore scored with the log-probability of the 6.75 A bin, and a pair at
16.0 A with that of the 21.0 A bin.  Both `ml_full` and `ORACLE_ml_true` in `s15/distml.py`
carry this defect.  It is a mis-indexing bug, not a property of binning, and it had to be
separated from the quantisation question before either could be answered.

STRUCTURE OF THIS MODULE

  section 1   representations of a per-pair distance law, all with an ANALYTIC derivative and
              the common interface  table(d) -> (logP, dlogP/dd)
                LogPGrid          piecewise-linear log P on the TRUE non-uniform centres
                                  (`density=True` divides bin mass by bin width first)
                KDEDensity        sum of Gaussians on the centres, bandwidth = alpha * width
                PchipDensity      monotone PCHIP through the CDF at the bin EDGES,
                                  differentiated -- preserves every bin's mass EXACTLY
                GMMDensity        2-3 component Gaussian mixture by grouped-data EM
                GaussDensity      an analytic Gaussian; ML against it IS weighted least
                                  squares, which is the point of the `_cont` ladder arms
  section 2   derivative verification against central differences (cosine, as K0 did)
  section 3   multimodality census of the predicted histograms
  section 4   the ORACLE quantisation ladder (task 1)
  section 5   the predictive arms (task 2), with leave-fold-out hyperparameter choice

Run:
    python -m s15.quant check      derivative cosines           (seconds)
    python -m s15.quant mm         multimodality census         (seconds)
    python -m s15.quant ladder     ORACLE quantisation ladder   (~1 h, 126 targets)
    python -m s15.quant pred       predictive arms              (~1 h, 126 targets)
    python -m s15.quant all        all four, in order
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
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import distml as M                  # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

# ------------------------------------------------------------------ the binning, from source
#: `core.predict.BIN_EDGES`, copied here so this module states its own assumption explicitly.
BIN_EDGES = np.array([4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0,
                      11.0, 12.5, 14.0, 16.0, 19.0, 23.0])
CENTRES = np.concatenate([[BIN_EDGES[0] - 0.5],
                          0.5 * (BIN_EDGES[1:] + BIN_EDGES[:-1]),
                          [BIN_EDGES[-1] + 2.0]])
#: bin 0 is (-inf, 4.5) and bin 16 is [23, inf); their centres imply nominal outer edges at
#: 3.5 and 27.0, which is what a density has to use to be normalisable at all.
KNOTS = np.concatenate([[3.5], BIN_EDGES, [27.0]])          # 18 knots, 17 bins
WIDTHS = np.diff(KNOTS)                                     # (17,)

#: a small uniform background over [2, 40] keeps log P finite and differentiable far outside
#: the predicted support.  exp(-12) * 38 ~ 2.3e-4, chosen so the effective floor matches the
#: LOGP_FLOOR = -12 that `s15/distml.py` uses, so the arms are comparable.
BG_LO, BG_HI = 2.0, 40.0
BG_EPS = 2.3e-4
BG_DENS = BG_EPS / (BG_HI - BG_LO)

SEP_BANDS = [("short_2_5", 2, 5), ("mid_6_7", 6, 7), ("long_8plus", 8, 10 ** 9)]


# ============================================================ 1. representations of P(d)
class LogPGrid:
    """Piecewise-linear log P(d) interpolated on the TRUE non-uniform bin centres.

    This is the corrected form of `s15.distml.LogPTable`, which assumed a uniform 0.75 A
    grid.  With `density=True` the bin MASS is divided by the bin WIDTH first, which is the
    only self-consistent thing to do when the bins are 0.5 A at short range and 4.0 A at
    long range -- otherwise a wide bin is rewarded merely for being wide.
    """

    def __init__(self, prob, centres=CENTRES, floor=M.LOGP_FLOOR, density=False):
        p = np.asarray(prob, float)
        p = np.maximum(p, 1e-12)
        p = p / p.sum(1, keepdims=True)
        if density:
            p = p / WIDTHS[None, :]
        self.lp = np.maximum(np.log(p), floor)              # (P, B)
        self.c = np.asarray(centres, float)
        self.B = len(self.c)
        self.h = np.diff(self.c)                            # (B-1,) non-uniform
        self.slope = (self.lp[:, 1:] - self.lp[:, :-1]) / self.h[None, :]

    def __call__(self, d):
        d = np.asarray(d, float)
        k = np.clip(np.searchsorted(self.c, d, side="right") - 1, 0, self.B - 2)
        rows = np.arange(len(d))
        g = self.slope[rows, k]
        val = self.lp[rows, k] + g * (d - self.c[k])
        inside = (d >= self.c[0]) & (d <= self.c[-1])
        return np.where(inside, val, M.LOGP_FLOOR), np.where(inside, g, 0.0)


class LogPTableUniformBUG:
    """A FROZEN copy of `s15.distml.LogPTable` as it stood before 2026-09-05.

    Kept deliberately, and never to be fixed.  It indexes the 17 non-uniform bin centres as
    if they were a uniform 0.75 A grid, so a query at 8.0 A reads the 6.75 A bin and a query
    at 16.0 A reads the 21.0 A bin.  The coordinator has since corrected the live class; this
    copy exists so the ORACLE ladder can measure HOW MUCH of the apparent 1.06 A
    "quantisation ceiling" was the bug and how much is real binning.  Its analytic derivative
    is self-consistent (cosine 1.0) -- the defect is the index map, not the calculus, which is
    exactly why it survived a gradient check.
    """

    def __init__(self, prob, centres=CENTRES, floor=M.LOGP_FLOOR):
        p = np.asarray(prob, float)
        p = np.maximum(p, 1e-12)
        p = p / p.sum(1, keepdims=True)
        self.lp = np.maximum(np.log(p), floor)
        self.c = np.asarray(centres, float)
        self.B = len(self.c)
        self.w = float(self.c[1] - self.c[0])           # <-- the defect: assumes uniform
        self.dlp = np.zeros_like(self.lp)
        self.dlp[:, :-1] = (self.lp[:, 1:] - self.lp[:, :-1]) / self.w

    def __call__(self, d):
        x = (d - self.c[0]) / self.w
        k = np.clip(np.floor(x).astype(int), 0, self.B - 2)
        t = np.clip(x - k, 0.0, 1.0)
        rows = np.arange(len(d))
        lo = self.lp[rows, k]
        g = self.dlp[rows, k]
        val = lo + g * (t * self.w)
        inside = (d >= self.c[0]) & (d <= self.c[-1])
        return np.where(inside, val, M.LOGP_FLOOR), np.where(inside, g, 0.0)


class KDEDensity:
    """p(d) = sum_b m_b N(d; c_b, h_b),  h_b = alpha * width_b, plus a uniform background.

    The bandwidth is tied to the LOCAL bin width, so the short-range 0.5 A bins stay sharp
    while the 4 A long-range bins are smeared by the amount their own width justifies.
    Fully analytic in d.
    """

    def __init__(self, prob, alpha=1.0, centres=CENTRES, widths=WIDTHS):
        p = np.asarray(prob, float)
        p = np.maximum(p, 1e-12)
        self.m = p / p.sum(1, keepdims=True)                # (P, B)
        self.c = np.asarray(centres, float)[None, :]
        self.h = np.maximum(alpha * np.asarray(widths, float), 1e-3)[None, :]
        self.alpha = float(alpha)
        self._norm = 1.0 / (self.h * np.sqrt(2.0 * np.pi))

    def __call__(self, d):
        z = (np.asarray(d, float)[:, None] - self.c) / self.h        # (P, B)
        t = self.m * self._norm * np.exp(-0.5 * z * z)
        s = t.sum(1)
        ds = (t * (-z / self.h)).sum(1)
        dens = (1.0 - BG_EPS) * s + BG_DENS
        ddens = (1.0 - BG_EPS) * ds
        return np.log(dens), ddens / dens


def _pchip_slopes(x, y):
    """Fritsch-Carlson monotone slopes.  x (K,) shared, y (P, K).  Returns m (P, K)."""
    h = np.diff(x)[None, :]                                  # (1, K-1)
    dl = np.diff(y, axis=1) / h                              # (P, K-1)
    m = np.zeros_like(y)
    if y.shape[1] > 2:
        d0, d1 = dl[:, :-1], dl[:, 1:]
        h0, h1 = h[:, :-1], h[:, 1:]
        w1 = 2.0 * h1 + h0
        w2 = h1 + 2.0 * h0
        ok = (d0 * d1) > 0
        den = np.where(ok, w1 / np.where(d0 == 0, 1.0, d0) + w2 / np.where(d1 == 0, 1.0, d1),
                       1.0)
        m[:, 1:-1] = np.where(ok, (w1 + w2) / den, 0.0)
    #: endpoint slopes set to the adjacent secant.  0 <= delta <= 3*delta, so monotonicity
    #: is preserved, and the density at the outer knots equals the outer bins' own density.
    m[:, 0] = dl[:, 0]
    m[:, -1] = dl[:, -1]
    return m


class PchipDensity:
    """Monotone PCHIP through the CDF at the bin EDGES, differentiated to a density.

    Because the interpolant passes exactly through the cumulative mass at every edge, the
    integral of the density over bin b is EXACTLY p_b.  Monotone interpolation cannot emit
    negative density.  The density is piecewise quadratic and its derivative piecewise
    linear, both analytic.
    """

    def __init__(self, prob, knots=KNOTS):
        p = np.asarray(prob, float)
        p = np.maximum(p, 1e-12)
        p = p / p.sum(1, keepdims=True)
        self.x = np.asarray(knots, float)                    # (K,)
        F = np.concatenate([np.zeros((len(p), 1)), np.cumsum(p, axis=1)], axis=1)   # (P, K)
        m = _pchip_slopes(self.x, F)
        h = np.diff(self.x)[None, :]
        y0, y1 = F[:, :-1], F[:, 1:]
        m0, m1 = m[:, :-1], m[:, 1:]
        self.h = h
        self.b = h * m0
        self.cc = 3.0 * (y1 - y0) - h * (2.0 * m0 + m1)
        self.e = -2.0 * (y1 - y0) + h * (m0 + m1)

    def __call__(self, d):
        d = np.asarray(d, float)
        K = len(self.x)
        k = np.clip(np.searchsorted(self.x, d, side="right") - 1, 0, K - 2)
        rows = np.arange(len(d))
        h = self.h[0, k]
        t = np.clip((d - self.x[k]) / h, 0.0, 1.0)
        b, c, e = self.b[rows, k], self.cc[rows, k], self.e[rows, k]
        dens = (b + 2.0 * c * t + 3.0 * e * t * t) / h                 # F'
        ddens = (2.0 * c + 6.0 * e * t) / (h * h)                      # F''
        inside = (d > self.x[0]) & (d < self.x[-1])
        dens = np.where(inside, dens, 0.0)
        ddens = np.where(inside, ddens, 0.0)
        dens = np.maximum(dens, 0.0)
        tot = (1.0 - BG_EPS) * dens + BG_DENS
        return np.log(tot), (1.0 - BG_EPS) * ddens / tot


class GMMDensity:
    """A 2-3 component Gaussian mixture fitted per pair to the histogram by grouped-data EM.

    Each bin contributes mass p_b located at c_b with a within-bin variance width_b^2 / 12
    (the variance of a uniform on the bin), which is added in the M-step.  That is the
    standard grouped-data correction and it is what stops the mixture from collapsing onto
    the bin centres -- i.e. it is what restores sub-bin resolution rather than re-imposing
    the grid.  Explicitly preserves bimodality, which is the mechanism ML is supposed to
    exploit and least squares cannot.
    """

    def __init__(self, prob, K=2, iters=60, sig_min=0.30, centres=CENTRES, widths=WIDTHS):
        p = np.asarray(prob, float)
        p = np.maximum(p, 1e-12)
        p = p / p.sum(1, keepdims=True)
        P, B = p.shape
        c = np.asarray(centres, float)
        v = (np.asarray(widths, float) ** 2) / 12.0
        # init: split the bins into K contiguous groups by cumulative mass
        F = np.cumsum(p, axis=1)
        mu = np.zeros((P, K)); pi = np.full((P, K), 1.0 / K); sg = np.zeros((P, K))
        for a in range(K):
            q = (a + 0.5) / K
            idx = np.argmin(np.abs(F - q), axis=1)
            mu[:, a] = c[idx] + 0.01 * a
        sg[:] = np.sqrt(((p * (c[None, :] - (p * c[None, :]).sum(1, keepdims=True)) ** 2)
                         ).sum(1, keepdims=True)) + 0.5
        sg = np.maximum(sg, sig_min)
        for _ in range(iters):
            # E step over the B weighted bin atoms
            z = (c[None, :, None] - mu[:, None, :]) / sg[:, None, :]
            lg = -0.5 * z * z - np.log(sg[:, None, :]) + np.log(np.maximum(pi, 1e-12))[:, None, :]
            lg -= lg.max(2, keepdims=True)
            r = np.exp(lg); r /= r.sum(2, keepdims=True)                # (P, B, K)
            wgt = r * p[:, :, None]                                    # (P, B, K)
            nk = wgt.sum(1) + 1e-12                                    # (P, K)
            pi = nk / nk.sum(1, keepdims=True)
            mu = (wgt * c[None, :, None]).sum(1) / nk
            var = ((wgt * ((c[None, :, None] - mu[:, None, :]) ** 2 + v[None, :, None])
                    ).sum(1) / nk)
            sg = np.maximum(np.sqrt(var), sig_min)
        self.mu, self.sg, self.pi = mu, sg, pi
        self.K = K
        self._norm = pi / (sg * np.sqrt(2.0 * np.pi))

    def __call__(self, d):
        z = (np.asarray(d, float)[:, None] - self.mu) / self.sg
        t = self._norm * np.exp(-0.5 * z * z)
        s = t.sum(1)
        ds = (t * (-z / self.sg)).sum(1)
        dens = (1.0 - BG_EPS) * s + BG_DENS
        return np.log(dens), (1.0 - BG_EPS) * ds / dens


class GaussDensity:
    """An analytic Gaussian per pair.  ML against this is EXACTLY weighted least squares
    with w = 1/(2 sigma^2), so it is the continuous control for the sampled-Gaussian arms."""

    def __init__(self, mu, sigma):
        self.mu = np.asarray(mu, float)
        self.sg = np.maximum(np.asarray(sigma, float), 1e-6)

    def __call__(self, d):
        z = (np.asarray(d, float) - self.mu) / self.sg
        return -0.5 * z * z - np.log(self.sg * np.sqrt(2.0 * np.pi)), -z / self.sg


# ============================================================ 2. derivative verification
def check_derivatives(pdb=None, eps=1e-5):
    """Central-difference check of every density's dlogP/dd.  Reports cosine, as K0 did."""
    tg = I.targets()
    t = tg[0] if pdb is None else [x for x in tg if x["pdb"] == pdb][0]
    dg = I.distogram(t["pdb"], t["seq"], int(t["fold"]))
    prob = np.asarray(dg["prob"], float)
    rng = np.random.default_rng(0)
    d = rng.uniform(4.0, 24.0, prob.shape[0])
    tabs = {
        "LogPGrid_mass": LogPGrid(prob),
        "LogPGrid_density": LogPGrid(prob, density=True),
        "KDE_alpha0.5": KDEDensity(prob, 0.5),
        "KDE_alpha1.0": KDEDensity(prob, 1.0),
        "KDE_alpha2.0": KDEDensity(prob, 2.0),
        "Pchip": PchipDensity(prob),
        "GMM_K2": GMMDensity(prob, 2),
        "GMM_K3": GMMDensity(prob, 3),
        "Gauss": GaussDensity(np.asarray(dg["expected"], float),
                              np.maximum(np.asarray(dg["sd"], float), 1e-3)),
        "distml.LogPTable_FIXED": M.LogPTable(prob, np.asarray(dg["centres"], float)),
        "LogPTableUniformBUG_frozen": LogPTableUniformBUG(prob),
    }
    out = {}
    for name, tb in tabs.items():
        _, g = tb(d)
        fp, _ = tb(d + eps)
        fm, _ = tb(d - eps)
        fd = (fp - fm) / (2 * eps)
        ok = np.isfinite(g) & np.isfinite(fd)
        # piecewise-linear/quadratic pieces have kinks; drop points within eps of a knot
        near = np.zeros_like(d, bool)
        #: kinks of every piecewise representation present, INCLUDING the frozen buggy
        #: table's spurious uniform-0.75 grid -- a central difference straddling a kink is
        #: not a derivative error, and excluding them is what lets the check say the bug is
        #: in the index map and not in the calculus.
        for kn in np.concatenate([CENTRES, KNOTS, CENTRES[0] + 0.75 * np.arange(20)]):
            near |= np.abs(d - kn) < 10 * eps
        ok &= ~near
        a, b = g[ok], fd[ok]
        cos = float((a * b).sum() / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-300))
        out[name] = {"cosine": cos, "max_abs_err": float(np.max(np.abs(a - b))),
                     "n": int(ok.sum())}
        print(f"{name:<30} cosine {cos:.12f}  max|err| {out[name]['max_abs_err']:.3e}"
              f"  n={out[name]['n']}")
    #: INDEPENDENT REPLICATION.  `LogPGrid` here and `s15.distml.LogPTable` were written
    #: separately, by two agents, from the same specification.  They agree to machine
    #: precision, which is a genuine cross-check of the corrected non-uniform lookup rather
    #: than two copies of one piece of code.
    a1, g1 = tabs["LogPGrid_mass"](d)
    a2, g2 = tabs["distml.LogPTable_FIXED"](d)
    out["independent_replication_LogPGrid_vs_distml"] = {
        "max_abs_logp_diff": float(np.max(np.abs(a1 - a2))),
        "max_abs_grad_diff": float(np.max(np.abs(g1 - g2)))}
    print(f"{'LogPGrid vs distml.LogPTable':<30} max|dlogP| "
          f"{np.max(np.abs(a1 - a2)):.3e}  max|dgrad| {np.max(np.abs(g1 - g2)):.3e}")
    # and the full torsion-gradient path, exactly as K0 did, for one density
    n = int(t["n"])
    i, j = I.pair_index(n)
    tb = KDEDensity(prob, 1.0)
    from core import project as pj

    def fg(x):
        phi, psi = x[:n], x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        r = CA[i] - CA[j]
        dd = np.maximum(np.sqrt((r * r).sum(1)), 1e-9)
        lp, dlp = tb(dd)
        coef = (-dlp / dd)[:, None] * r
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        return -float(lp.sum()), pj._torsion_grad(G, CA, gCA)

    x0 = rng.uniform(-np.pi, np.pi, 2 * n)
    f0, g0 = fg(x0)
    h = 1e-6
    fd = np.zeros_like(x0)
    for k in range(len(x0)):
        xp = x0.copy(); xp[k] += h
        xm = x0.copy(); xm[k] -= h
        fd[k] = (fg(xp)[0] - fg(xm)[0]) / (2 * h)
    live = np.abs(fd) > 1e-8
    cos = float((g0[live] * fd[live]).sum()
                / (np.linalg.norm(g0[live]) * np.linalg.norm(fd[live])))
    out["torsion_gradient_KDE_alpha1.0"] = {"cosine": cos, "n": int(live.sum())}
    print(f"{'torsion gradient (KDE a=1.0)':<30} cosine {cos:.12f}  n={int(live.sum())}")
    with open(os.path.join(RESULTS, "quant_derivcheck.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ============================================================ 3. multimodality census
def multimodality(min_gap_bins=2, dip_ratio=0.5, min_minor_mass=0.15, basis="density"):
    """How often is a predicted per-pair histogram GENUINELY bimodal?

    Operational definition, stated so it can be argued with.  Work on the DENSITY
    h_b = p_b / width_b -- the raw mass is not comparable across bins that differ eightfold
    in width, and using mass would manufacture spurious modes at every wide bin.  A pair is
    called bimodal when it has two local density maxima at least `min_gap_bins` apart, the
    minimum density strictly between them is at most `dip_ratio` times the smaller peak
    (a real valley, not a shoulder), and the smaller mode's basin -- the bins on its side of
    the valley -- carries at least `min_minor_mass` of the total probability.

    A looser count (any two local maxima) is reported beside it so the reader can see how
    much of the answer the strictness is doing.
    """
    tg = I.targets()
    rows = []
    for t in tg:
        pdb, n, fold = t["pdb"], int(t["n"]), int(t["fold"])
        dg = I.distogram(pdb, t["seq"], fold)
        p = np.asarray(dg["prob"], float)
        p = p / p.sum(1, keepdims=True)
        h = p / WIDTHS[None, :] if basis == "density" else p
        i, j = I.pair_index(n)
        sep = (j - i).astype(int)
        B = h.shape[1]
        for r in range(h.shape[0]):
            v = h[r]
            loc = [b for b in range(B)
                   if (b == 0 or v[b] > v[b - 1]) and (b == B - 1 or v[b] > v[b + 1])]
            strict = False
            if len(loc) >= 2:
                for a in range(len(loc) - 1):
                    b1, b2 = loc[a], loc[a + 1]
                    if b2 - b1 < min_gap_bins:
                        continue
                    seg = v[b1 + 1:b2]
                    if len(seg) == 0:
                        continue
                    vmin = seg.min()
                    if vmin > dip_ratio * min(v[b1], v[b2]):
                        continue
                    cut = b1 + 1 + int(np.argmin(seg))
                    m_lo = p[r, :cut].sum(); m_hi = p[r, cut:].sum()
                    if min(m_lo, m_hi) >= min_minor_mass:
                        strict = True
                        break
            rows.append((sep[r], len(loc) >= 2, strict, len(loc)))
    rows = np.array(rows, float)
    out = {"definition": {"min_gap_bins": min_gap_bins, "dip_ratio": dip_ratio,
                          "min_minor_mass": min_minor_mass,
                          "basis": basis},
           "bands": {}}
    print(f"{'band':<14}{'n pairs':>9}{'>=2 local maxima':>19}{'STRICT bimodal':>17}"
          f"{'mean n modes':>14}")
    for name, lo, hi in SEP_BANDS + [("ALL", 2, 10 ** 9)]:
        m = (rows[:, 0] >= lo) & (rows[:, 0] <= hi)
        out["bands"][name] = {"n": int(m.sum()), "frac_two_maxima": float(rows[m, 1].mean()),
                             "frac_strict_bimodal": float(rows[m, 2].mean()),
                             "mean_n_modes": float(rows[m, 3].mean())}
        b = out["bands"][name]
        print(f"{name:<14}{b['n']:>9}{b['frac_two_maxima']:>19.3f}"
              f"{b['frac_strict_bimodal']:>17.3f}{b['mean_n_modes']:>14.2f}")
    with open(os.path.join(RESULTS, f"quant_multimodal_{basis}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ============================================ 3b. the analytic quantisation-noise budget
def noise_budget():
    """An arithmetic check on the whole hypothesis, independent of any fit.

    Rounding a distance to its bin adds an error whose variance is at most the variance of a
    uniform on that bin, width^2 / 12.  Compare that with the distogram's MEASURED error at
    the same separation (s15/distacc.py, reproduced in coord_FINDINGS K2).  If quantisation
    noise is a small share of the error variance the predictor already has, then no
    reconstruction of the histogram can buy much, and the hypothesis is arithmetically
    capped before a single structure is fitted.

    ORACLE only in the sense that the reference RMSEs were scored post hoc against native
    distances; nothing here is used at inference.
    """
    tg = I.targets()
    # which bins do the true distances of each separation band actually fall in?
    rows = []
    for t in tg:
        n = int(t["n"])
        u = I.load_univ(t["pdb"])
        nat = np.asarray(u["nat_ca"], float)
        i, j = I.pair_index(n)
        dt = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        b = np.digitize(dt, BIN_EDGES)
        rows.append(np.column_stack([(j - i).astype(float), WIDTHS[b]]))
    rows = np.vstack(rows)
    #: measured RMSE by separation, from coord_FINDINGS K2 (s15/distacc.py, 126 targets)
    ref_rmse = {"short_2_5": None, "mid_6_7": 4.001, "long_8plus": None}
    ref_by_sep = [(2, 3, 1.463), (4, 5, 2.963), (6, 7, 4.001), (8, 10, 5.216), (11, 999, 6.328)]
    out = {"bands": {}, "by_sep_bin": []}
    print(f"{'band':<14}{'n':>7}{'mean bin width':>16}{'quant sd (w/sqrt12)':>22}"
          f"{'distogram RMSE':>16}{'quant share of var':>20}")
    for name, lo, hi in SEP_BANDS + [("ALL", 2, 10 ** 9)]:
        m = (rows[:, 0] >= lo) & (rows[:, 0] <= hi)
        w = rows[m, 1]
        qsd = float(np.sqrt((w ** 2 / 12.0).mean()))
        # separation-weighted reference RMSE over the same pairs
        r = np.zeros(int(m.sum()))
        s = rows[m, 0]
        for a, b, v in ref_by_sep:
            r[(s >= a) & (s <= b)] = v
        rmse = float(np.sqrt((r ** 2).mean()))
        share = qsd ** 2 / rmse ** 2
        out["bands"][name] = {"n": int(m.sum()), "mean_bin_width": float(w.mean()),
                              "quant_sd": qsd, "distogram_rmse": rmse,
                              "quant_share_of_error_variance": float(share)}
        print(f"{name:<14}{int(m.sum()):>7}{w.mean():>16.2f}{qsd:>22.3f}"
              f"{rmse:>16.3f}{share:>20.3f}")
    out["ref_rmse_source"] = "s15/distacc.py, 126 targets, coord_FINDINGS K2"
    out["_unused"] = ref_rmse
    with open(os.path.join(RESULTS, "quant_noisebudget.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ============================================================ shared fitting harness
def _prep(t):
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    i, j = I.pair_index(n)
    dg = I.distogram(pdb, seq, fold)
    return {"pdb": pdb, "seq": seq, "n": n, "fold": fold, "nat": nat, "i": i, "j": j,
            "sep": (j - i).astype(float),
            "dtrue": np.sqrt(((nat[i] - nat[j]) ** 2).sum(1)),
            "prob": np.asarray(dg["prob"], float),
            "dhat": np.asarray(dg["expected"], float),
            "sd": np.maximum(np.asarray(dg["sd"], float), 1e-3)}


def _best_ls(d, target, w, S):
    best = None
    for phi0, psi0, _t in S:
        phi, psi, f, _ = D.fit_distances(target, w, d["i"], d["j"], phi0, psi0)
        if best is None or f < best[0]:
            best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), d["nat"])))
    return best[1]


def _best_ml(d, table, S):
    best = None
    for phi0, psi0, _t in S:
        phi, psi, f = M.fit_ml(table, d["i"], d["j"], phi0, psi0)
        if best is None or f < best[0]:
            best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), d["nat"])))
    return best[1]


def _report(res, arms, pdbs, folds, base_arm, name, extra=None):
    fail = np.isin(pdbs, I.FAIL18)
    from s14 import ladder as L
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res[base_arm], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "base_arm": base_arm,
           "arms": {}, "per_target": {a: dict(zip(pdbs, map(float, res[a]))) for a in arms}}
    if extra:
        out.update(extra)
    for a in arms:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_base": I.paired(v, base, folds=folds, names=pdbs),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    with open(os.path.join(RESULTS, name + ".json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_" + name, out, n_expected=len(pdbs))
    print(f"\nincumbent {ref.mean():.3f}    base arm = {base_arm}\n")
    print(f"{'arm':<34}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}{'FAIL18':>9}"
          f"{'vs base (paired, 95% CI)':>28}")
    for a in arms:
        s = out["arms"][a]
        vb = s["vs_base"]
        print(f"{a:<34}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"   {vb['mean_diff']:+.3f} [{vb['ci95'][0]:+.3f},{vb['ci95'][1]:+.3f}]"
              f" W/L {vb['n_better']}/{vb['n_worse']}")
    return out


# ============================================================ 4. the ORACLE ladder
def _snap_bin(d):
    """The quantisation the network's output layer actually performs."""
    return CENTRES[np.digitize(d, BIN_EDGES)]


def _snap_uniform(d, g):
    return np.round(np.asarray(d, float) / g) * g


def ladder(targets=None, n_start=6, gauss_w=(0.3, 0.6, 1.2)):
    """ORACLE DIAGNOSTIC.  Every arm knows the TRUE distance matrix; only the REPRESENTATION
    of that knowledge changes.  Nothing here can ever be a predictive headline."""
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)

    arms = ["ORACLE_continuous_ls", "ORACLE_snap_bin", "ORACLE_snap_uniform_0.5",
            "ORACLE_snap_uniform_0.25", "ORACLE_snap_bin_SHORT_2_5_only",
            "ORACLE_snap_bin_MID_6_7_only", "ORACLE_snap_bin_LONG_8plus_only"]
    arms += [f"ORACLE_ml_gauss_w{w}_binned" for w in gauss_w]
    arms += ["ORACLE_ml_gauss_w0.6_continuous", "ORACLE_ml_true_VOID_uniformgrid_BUG",
             "ORACLE_ml_gauss_w0.6_kde", "ORACLE_ml_gauss_w0.6_pchip"]
    res = {a: [] for a in arms}
    path = os.path.join(RESULTS, "quant_ladder.json")
    t0 = time.time()

    for c, t in enumerate(tg):
        d = _prep(t)
        rng = SD.stable_rng(d["pdb"], "ladder")
        S = D.starts(d["pdb"], d["seq"], d["n"], d["fold"], n_start, rng)
        dt, sep = d["dtrue"], d["sep"]
        one = np.ones_like(dt)
        sb = _snap_bin(dt)
        band = {n_: (sep >= lo) & (sep <= hi) for n_, lo, hi in SEP_BANDS}

        res["ORACLE_continuous_ls"].append(_best_ls(d, dt, one, S))
        res["ORACLE_snap_bin"].append(_best_ls(d, sb, one, S))
        res["ORACLE_snap_uniform_0.5"].append(_best_ls(d, _snap_uniform(dt, 0.5), one, S))
        res["ORACLE_snap_uniform_0.25"].append(_best_ls(d, _snap_uniform(dt, 0.25), one, S))
        for nm, key in (("short_2_5", "ORACLE_snap_bin_SHORT_2_5_only"),
                        ("mid_6_7", "ORACLE_snap_bin_MID_6_7_only"),
                        ("long_8plus", "ORACLE_snap_bin_LONG_8plus_only")):
            res[key].append(_best_ls(d, np.where(band[nm], sb, dt), one, S))

        for w in gauss_w:
            tp = np.exp(-0.5 * ((CENTRES[None, :] - dt[:, None]) / w) ** 2)
            res[f"ORACLE_ml_gauss_w{w}_binned"].append(
                _best_ml(d, LogPGrid(tp, density=False), S))
        res["ORACLE_ml_gauss_w0.6_continuous"].append(
            _best_ml(d, GaussDensity(dt, np.full_like(dt, 0.6)), S))
        tp06 = np.exp(-0.5 * ((CENTRES[None, :] - dt[:, None]) / 0.6) ** 2)
        #: the arm that reproduces the VOID K4 number, through the frozen pre-fix table
        res["ORACLE_ml_true_VOID_uniformgrid_BUG"].append(
            _best_ml(d, LogPTableUniformBUG(tp06), S))
        #: and the same binned law read back through two continuous reconstructions -- this
        #: is the "does a continuous density undo the binning?" question in ORACLE form,
        #: where the answer is not confounded by prediction error at all
        res["ORACLE_ml_gauss_w0.6_kde"].append(_best_ml(d, KDEDensity(tp06, 1.0), S))
        res["ORACLE_ml_gauss_w0.6_pchip"].append(_best_ml(d, PchipDensity(tp06), S))

        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1, "pdbs": pdbs[:c + 1]}, fh)
            print(f"  ladder {c+1}/{len(tg)}  {time.time()-t0:.0f}s", flush=True)

    out = _report(res, arms, pdbs, folds, "ORACLE_continuous_ls", "quant_ladder")
    # separation-band decomposition of the total quantisation cost
    tot = (np.asarray(res["ORACLE_snap_bin"]) - np.asarray(res["ORACLE_continuous_ls"])).mean()
    parts = {k: float((np.asarray(res[v]) - np.asarray(res["ORACLE_continuous_ls"])).mean())
             for k, v in (("short_2_5", "ORACLE_snap_bin_SHORT_2_5_only"),
                          ("mid_6_7", "ORACLE_snap_bin_MID_6_7_only"),
                          ("long_8plus", "ORACLE_snap_bin_LONG_8plus_only"))}
    print(f"\nquantisation cost, all pairs binned: {tot:+.3f} A")
    for k, v in parts.items():
        print(f"  {k:<12} binned alone: {v:+.3f} A")
    print(f"  sum of parts {sum(parts.values()):+.3f} vs joint {tot:+.3f}")
    out["band_decomposition"] = {"joint": float(tot), "parts": parts,
                                 "sum_of_parts": float(sum(parts.values()))}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ============================================================ 5. the predictive arms
def predictive(targets=None, n_start=6, alphas=(0.5, 1.0, 2.0), Ks=(2, 3)):
    """Native-free arms.  The only ORACLE quantity anywhere is the post-hoc RMSD."""
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)

    tuned = [f"ml_kde_a{a}" for a in alphas] + [f"ml_gmm_K{k}" for k in Ks]
    arms = (["ls_mean_sd", "ml_full_VOID_uniformgrid_BUG", "ml_full_fixedgrid_mass",
             "ml_full_fixedgrid_density"] + tuned + ["ml_pchip"])
    res = {a: [] for a in arms}
    path = os.path.join(RESULTS, "quant_pred.json")
    t0 = time.time()

    for c, t in enumerate(tg):
        d = _prep(t)
        rng = SD.stable_rng(d["pdb"], "pred")
        S = D.starts(d["pdb"], d["seq"], d["n"], d["fold"], n_start, rng)
        pr = d["prob"]
        res["ls_mean_sd"].append(_best_ls(d, d["dhat"], 1.0 / d["sd"] ** 2, S))
        res["ml_full_VOID_uniformgrid_BUG"].append(_best_ml(d, LogPTableUniformBUG(pr), S))
        res["ml_full_fixedgrid_mass"].append(_best_ml(d, LogPGrid(pr, density=False), S))
        res["ml_full_fixedgrid_density"].append(_best_ml(d, LogPGrid(pr, density=True), S))
        for a in alphas:
            res[f"ml_kde_a{a}"].append(_best_ml(d, KDEDensity(pr, a), S))
        for k in Ks:
            res[f"ml_gmm_K{k}"].append(_best_ml(d, GMMDensity(pr, k), S))
        res["ml_pchip"].append(_best_ml(d, PchipDensity(pr), S))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1, "pdbs": pdbs[:c + 1]}, fh)
            print(f"  pred {c+1}/{len(tg)}  {time.time()-t0:.0f}s", flush=True)

    # ---- leave-fold-out hyperparameter choice.  THE honest control: alpha and K are never
    #      chosen using any target they are then reported on.
    lfo = {}
    for gname, grid in (("ml_kde_LFO", [f"ml_kde_a{a}" for a in alphas]),
                        ("ml_gmm_LFO", [f"ml_gmm_K{k}" for k in Ks]),
                        ("ml_cont_LFO", [f"ml_kde_a{a}" for a in alphas]
                         + [f"ml_gmm_K{k}" for k in Ks] + ["ml_pchip"])):
        v = np.zeros(len(pdbs)); pick = {}
        for f in sorted(set(folds.tolist())):
            tr = folds != f
            means = {a: float(np.asarray(res[a])[tr].mean()) for a in grid}
            best = min(means, key=means.get)
            pick[int(f)] = best
            v[folds == f] = np.asarray(res[best])[folds == f]
        res[gname] = v.tolist()
        arms.append(gname)
        lfo[gname] = pick
        print(f"{gname}: per-fold choice {pick}")

    out = _report(res, arms, pdbs, folds, "ls_mean_sd", "quant_pred",
                  extra={"lfo_choice": lfo,
                         "tuned_arms": tuned + ["ml_pchip"],
                         "note": "arms named ml_kde_a*/ml_gmm_K* are EXPLORATORY (their "
                                 "hyperparameter is read on the same targets they report). "
                                 "Only the *_LFO arms are honest; they pick per held-out "
                                 "fold using the other four folds only."})
    out["elapsed_s"] = time.time() - t0
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


class GaussDensityNoConst:
    """The SAME Gaussian log-density with the normalisation constant dropped.

    `log N(d; mu, s) = -z^2/2 - log(s sqrt(2 pi))`.  The second term does not depend on `d`,
    so it changes neither the objective's argmin nor its gradient -- it is, mathematically, a
    no-op.  It is here because it turns out not to be a numerical no-op: see `tol_check`.
    """

    def __init__(self, mu, sigma):
        self.mu = np.asarray(mu, float)
        self.sg = np.maximum(np.asarray(sigma, float), 1e-6)

    def __call__(self, d):
        z = (np.asarray(d, float) - self.mu) / self.sg
        return -0.5 * z * z, -z / self.sg


def tol_check(targets=None, n_start=6, sigma=0.6):
    """ORACLE DIAGNOSTIC, and a methodological control the whole sprint needs.

    Maximum likelihood against a Gaussian of CONSTANT width IS unweighted least squares:
    `-sum log N(d_ij; dhat_ij, s) = (1/2s^2) * sum (d_ij - dhat_ij)^2 + P log(s sqrt(2 pi))`,
    an affine map of the least-squares objective with a positive slope.  Identical argmin,
    identical minimiser, identical everything -- mathematically.

    Three arms, all fitting the TRUE distance matrix from the SAME starts:

        ls_ref            sum (d - dtrue)^2                     (s15/distgeo.fit_distances)
        ml_gauss_const    -sum log N(d; dtrue, sigma)           (s15/distml.fit_ml)
        ml_gauss_noconst  the same MINUS the constant term      (s15/distml.fit_ml)

    Any difference among the three is pure optimiser behaviour.  It bounds how much of every
    "ML versus least squares" comparison in Sprint 15 is the objective and how much is
    L-BFGS-B's convergence test, whose `ftol` is RELATIVE and therefore sensitive to an
    additive constant that the mathematics says is irrelevant.
    """
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    arms = ["ls_ref", "ml_gauss_const", "ml_gauss_noconst"]
    res = {a: [] for a in arms}
    path = os.path.join(RESULTS, "quant_tolcheck.json")
    for c, t in enumerate(tg):
        d = _prep(t)
        S = D.starts(d["pdb"], d["seq"], d["n"], d["fold"], n_start,
                     SD.stable_rng(d["pdb"], "tolcheck"))
        dt = d["dtrue"]
        sg = np.full_like(dt, sigma)
        res["ls_ref"].append(_best_ls(d, dt, np.ones_like(dt), S))
        res["ml_gauss_const"].append(_best_ml(d, GaussDensity(dt, sg), S))
        res["ml_gauss_noconst"].append(_best_ml(d, GaussDensityNoConst(dt, sg), S))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1}, fh)
            print(f"  tolcheck {c+1}/{len(tg)}", flush=True)
    out = _report(res, arms, pdbs, folds, "ls_ref", "quant_tolcheck",
                  extra={"sigma": sigma,
                         "note": "all three arms have the SAME argmin by construction; any "
                                 "difference is L-BFGS-B convergence behaviour, not the "
                                 "objective."})
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ================================ 5b. how much noise does the unseeded start draw carry?
def start_noise(targets=None, n_start=6, seeds=(0, 1, 2, 3)):
    """`s15/distgeo.py`, `distml.py`, `distcal.py` and this module all seed the multi-start
    RNG with `hash(pdb)`.  Python's string hash is SALTED PER PROCESS unless PYTHONHASHSEED
    is set, and `s15/results/audit_env_lock.json` records it as null -- so every run draws a
    DIFFERENT set of pool starts.

    Within one process all arms share the same starts, so a paired within-run comparison is
    unaffected.  ACROSS processes it is not: any figure quoted from one module and compared
    with a figure from another is comparing different start sets.  This measures how large
    that is, on the cleanest arm in the sprint (ORACLE fit to the true distance matrix), so
    every cross-module difference in Sprint 15 can be read against a noise floor.
    """
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    per = {}
    for s in seeds:
        v = []
        for t in tg:
            d = _prep(t)
            rng = SD.stable_rng(d["pdb"], "startnoise", s)
            S = D.starts(d["pdb"], d["seq"], d["n"], d["fold"], n_start, rng)
            v.append(_best_ls(d, d["dtrue"], np.ones_like(d["dtrue"]), S))
        per[s] = v
        print(f"  seed {s}: ORACLE_continuous_ls mean {np.mean(v):.4f}  "
              f"median {np.median(v):.4f}", flush=True)
    A = np.array([per[s] for s in seeds], float)
    out = {"n": len(pdbs), "seeds": list(seeds),
           "arm_means": [float(x) for x in A.mean(1)],
           "sd_of_mean_across_seeds": float(A.mean(1).std(ddof=1)),
           "range_of_mean": float(A.mean(1).max() - A.mean(1).min()),
           "mean_per_target_sd": float(A.std(0, ddof=1).mean()),
           "max_per_target_range": float((A.max(0) - A.min(0)).max()),
           "per_seed": {str(s): dict(zip(pdbs, map(float, per[s]))) for s in seeds}}
    print(f"\nmean across seeds: sd {out['sd_of_mean_across_seeds']:.4f} A, "
          f"range {out['range_of_mean']:.4f} A")
    print(f"per-target: mean sd {out['mean_per_target_sd']:.4f} A, "
          f"worst range {out['max_per_target_range']:.4f} A")
    with open(os.path.join(RESULTS, "quant_startnoise.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


def ladder_repro(targets=None, n_start=6, seeds=(1, 2)):
    """INDEPENDENT VERIFICATION of the headline +0.381 A quantisation cost.

    `start_noise` shows the MEAN of an arm swings by ~0.25 A across start draws.  That does
    not by itself invalidate a PAIRED difference, because both arms in a run share one start
    set and the shared component cancels.  This re-runs only the two arms that carry the
    headline -- continuous least squares and least squares on bin-snapped distances -- at
    independent start draws, and reports the paired difference each time.  If the paired
    difference is stable while the means are not, the headline survives and the instability
    is confined to absolute levels.
    """
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    out = {"seeds": list(seeds), "runs": {}}
    for s in seeds:
        cont, snap = [], []
        for t in tg:
            d = _prep(t)
            S = D.starts(d["pdb"], d["seq"], d["n"], d["fold"], n_start,
                         SD.stable_rng(d["pdb"], "ladder_repro", s))
            dt = d["dtrue"]; one = np.ones_like(dt)
            cont.append(_best_ls(d, dt, one, S))
            snap.append(_best_ls(d, _snap_bin(dt), one, S))
        p = I.paired(np.asarray(snap), np.asarray(cont), folds=folds, names=pdbs)
        out["runs"][str(s)] = {"continuous_mean": float(np.mean(cont)),
                              "snap_mean": float(np.mean(snap)), "paired": p}
        print(f"  seed {s}: continuous {np.mean(cont):.3f}  snap_bin {np.mean(snap):.3f}  "
              f"paired {p['mean_diff']:+.3f} [{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
              f"W/L {p['n_better']}/{p['n_worse']}", flush=True)
        with open(os.path.join(RESULTS, "quant_ladder_repro.json"), "w") as fh:
            json.dump(out, fh, indent=1)
    return out


# ====================================== 6. does the ML gain track the stated mechanism?
def per_target_ipred():
    """Per-target `I_pred` in bits/pair, the INFO workstream's calibration statistic.

    Reimplemented here from `s15/info_mi.distogram_information`'s definition -- bits the
    distogram adds over a leave-one-target-out pooled distance marginal, on the native bin
    label -- because `s15/results/info_mi.json` stores only the aggregate.  INFO measures
    **-0.818 bits/pair [-1.265, -0.419], positive on only 52 of 126 targets**: the distogram
    is a good ranker and a BAD PROBABILITY.  Maximum likelihood consumes the distribution AS
    a probability; least squares consumes only its first two moments.  So this is the
    variable that should predict where ML wins, if the INFO account is right.
    """
    from s15.info_mi import smooth
    tg = I.targets()
    data = []
    for t in tg:
        dg = I.distogram(t["pdb"], t["seq"], int(t["fold"]))
        u = I.load_univ(t["pdb"])
        nat = np.asarray(u["nat_ca"], float)
        i, j = I.pair_index(int(t["n"]))
        data.append({"pdb": t["pdb"], "prob": np.asarray(dg["prob"], float),
                     "dtrue": np.sqrt(((nat[i] - nat[j]) ** 2).sum(1)),
                     "centres": np.asarray(dg["centres"], float)})
    cen = data[0]["centres"]
    edges = np.concatenate([[-np.inf], (cen[1:] + cen[:-1]) / 2, [np.inf]])
    bins_all = [np.clip(np.digitize(d["dtrue"], edges) - 1, 0, len(cen) - 1) for d in data]
    tot = np.zeros(len(cen))
    for b in bins_all:
        tot += np.bincount(b, minlength=len(cen))
    out = {}
    for k, d in enumerate(data):
        h0 = smooth(tot - np.bincount(bins_all[k], minlength=len(cen)), len(cen))
        q = np.maximum(d["prob"], 1e-9)
        q = q / q.sum(1, keepdims=True)
        b = bins_all[k]
        out[d["pdb"]] = float(np.mean(-np.log2(h0[b]) + np.log2(q[np.arange(len(b)), b])))
    return out


def mechanism(pred_file="quant_pred.json", arm="ml_full_fixedgrid_density"):
    """The falsification test for the ML-over-least-squares mechanism, and the INFO account.

    If ML beats least squares BECAUSE it can represent asymmetry and bimodality that a mean
    and an sd cannot, then the per-target advantage of the ML arm must be larger on targets
    whose predicted histograms are more often genuinely bimodal, and larger on targets whose
    histograms are more skewed.  If instead the INFO account is right -- ML loses because it
    faithfully fits a MISCALIBRATED density while least squares discards the miscalibrated
    part -- then the advantage must track per-target `I_pred`, and ML should win on the 52
    targets where `I_pred > 0`.  Three Spearman correlations with permutation nulls, plus the
    split-mean test on the sign of `I_pred`.  A null on the first two with a signal on the
    third converts a negative result into a mechanism.
    """
    from scipy.stats import spearmanr
    P = json.load(open(os.path.join(RESULTS, pred_file)))
    per = P["per_target"]
    pdbs = list(per["ls_mean_sd"].keys())
    diff = np.array([per[arm][p] - per["ls_mean_sd"][p] for p in pdbs])   # neg = ML better
    frac, skew = [], []
    for t in I.targets():
        if t["pdb"] not in per["ls_mean_sd"]:
            continue
        dg = I.distogram(t["pdb"], t["seq"], int(t["fold"]))
        p = np.asarray(dg["prob"], float)
        p = p / p.sum(1, keepdims=True)
        h = p / WIDTHS[None, :]
        B = h.shape[1]
        nb = 0
        for r in range(h.shape[0]):
            v = h[r]
            loc = [b for b in range(B)
                   if (b == 0 or v[b] > v[b - 1]) and (b == B - 1 or v[b] > v[b + 1])]
            hit = False
            for a in range(len(loc) - 1):
                b1, b2 = loc[a], loc[a + 1]
                if b2 - b1 < 2:
                    continue
                seg = v[b1 + 1:b2]
                if len(seg) and seg.min() <= 0.5 * min(v[b1], v[b2]):
                    cut = b1 + 1 + int(np.argmin(seg))
                    if min(p[r, :cut].sum(), p[r, cut:].sum()) >= 0.15:
                        hit = True
                        break
            nb += hit
        frac.append(nb / h.shape[0])
        mu = (p * CENTRES[None, :]).sum(1)
        sd = np.sqrt((p * (CENTRES[None, :] - mu[:, None]) ** 2).sum(1))
        m3 = (p * (CENTRES[None, :] - mu[:, None]) ** 3).sum(1)
        skew.append(float(np.abs(m3 / np.maximum(sd, 1e-6) ** 3).mean()))
    frac = np.asarray(frac); skew = np.asarray(skew)
    ip = per_target_ipred()
    ipred = np.array([ip[p] for p in pdbs], float)
    out = {"arm": arm, "n": len(pdbs),
           "ipred_mean_bits_per_pair": float(ipred.mean()),
           "ipred_median": float(np.median(ipred)),
           "n_ipred_positive": int((ipred > 0).sum())}
    print(f"per-target I_pred reproduced: mean {ipred.mean():+.3f} bits/pair, "
          f"median {np.median(ipred):+.3f}, positive on {int((ipred>0).sum())}/{len(ipred)} "
          f"(INFO reports -0.818 / -0.138 / 52)")
    # the INFO account's specific prediction: ML should win where the density is calibrated
    a, b = diff[ipred > 0], diff[ipred <= 0]
    out["split_on_ipred_sign"] = {
        "n_calibrated": int(len(a)), "n_miscalibrated": int(len(b)),
        "mean_diff_where_ipred_positive": float(a.mean()) if len(a) else None,
        "mean_diff_where_ipred_negative": float(b.mean()) if len(b) else None,
        "paired_within_calibrated": I.paired(
            np.array([per[arm][p] for p, q in zip(pdbs, ipred) if q > 0]),
            np.array([per["ls_mean_sd"][p] for p, q in zip(pdbs, ipred) if q > 0])),
        "paired_within_miscalibrated": I.paired(
            np.array([per[arm][p] for p, q in zip(pdbs, ipred) if q <= 0]),
            np.array([per["ls_mean_sd"][p] for p, q in zip(pdbs, ipred) if q <= 0]))}
    sc = out["split_on_ipred_sign"]
    print(f"{arm} - ls_mean_sd:  where I_pred>0 (n={len(a)}) "
          f"{sc['paired_within_calibrated']['mean_diff']:+.3f} "
          f"{sc['paired_within_calibrated']['ci95']}   "
          f"where I_pred<=0 (n={len(b)}) "
          f"{sc['paired_within_miscalibrated']['mean_diff']:+.3f} "
          f"{sc['paired_within_miscalibrated']['ci95']}")
    rng = np.random.default_rng(0)
    for nm, x in (("frac_strict_bimodal", frac), ("mean_abs_skewness", skew),
                  ("ipred_bits_per_pair", ipred)):
        rho, pv = spearmanr(x, diff)
        perm = np.array([spearmanr(x, rng.permutation(diff))[0] for _ in range(2000)])
        out[nm] = {"spearman_rho": float(rho), "p": float(pv),
                   "perm_null_2.5_97.5": [float(np.percentile(perm, 2.5)),
                                          float(np.percentile(perm, 97.5))],
                   "x_mean": float(x.mean()), "x_range": [float(x.min()), float(x.max())]}
        print(f"rho( {nm:<22}, {arm} - ls ) = {rho:+.3f}  p={pv:.3f}  "
              f"perm null [{np.percentile(perm,2.5):+.3f},{np.percentile(perm,97.5):+.3f}]")
    with open(os.path.join(RESULTS, f"quant_mechanism_{arm}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    tg = I.targets()[:lim] if lim else None
    if cmd in ("check", "all"):
        print("=== derivative verification ===");      check_derivatives()
    if cmd in ("mm", "all"):
        print("\n=== multimodality census (density basis) ===");   multimodality()
        print("\n=== multimodality census (raw MASS basis, control) ===")
        multimodality(basis="mass")
    if cmd in ("budget", "all"):
        print("\n=== analytic quantisation-noise budget ==="); noise_budget()
    if cmd in ("ladder", "all"):
        print("\n=== ORACLE quantisation ladder ==="); ladder(tg)
    if cmd in ("pred", "all"):
        print("\n=== predictive arms ===");            predictive(tg)
    if cmd in ("tolcheck",):
        print("\n=== least squares IS Gaussian ML: optimiser control ==="); tol_check(tg)
    if cmd in ("repro",):
        print("\n=== headline reproduction at independent starts ==="); ladder_repro(tg)
    if cmd in ("startnoise",):
        print("\n=== multi-start draw noise ===");     start_noise(tg)
    if cmd in ("mech", "all"):
        print("\n=== mechanism test ===");             mechanism()
