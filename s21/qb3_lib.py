"""SPRINT 21 / WORKSTREAM B -- CVaR-VQE MECHANICS.

Pre-registered in `s21/PREREG_B.md` before this file existed.

THREE BLOCKS
------------
E   the ENCODING, as a GAUGE question.  `theta` and `(cos theta, sin theta)` are two charts of
    the same torus, and each carries a one-parameter family of PHYSICALLY IDENTICAL
    reparameterisations -- `theta = s*x` for the angle chart, `u = r*(cos, sin)` for the circle
    chart.  `arctan2` is scale-invariant, so `r` is pure gauge; `s` is a linear chart change.
    Neither moves the objective, its argmin, or any level set.  **If the outcome moves along a
    gauge orbit, the Sprint-20 encoding result is a step-size result.**

A   the ANSATZ against its own EXACT limit.  `core.quantum.MPSAnsatz` sets `chi = 2**layers` in
    the source, so the deployed `layers=2` is bond dimension 4 and depth is the only knob.  Bond
    dimension is MEASURED here, not quoted: the amplitude tensor is contracted exactly and the
    Schmidt rank read at every bipartition.

Q   the OPTIMIZER, including a Fisher-preconditioned natural gradient.  For a SAMPLING objective
    optimised by a score-function estimator the natural metric is the classical Fisher
    information of `p_theta(b)`, estimable from the SAME shots at ZERO extra objective
    evaluations.  It is not the Fubini-Study metric of the state and is not called one.

WHAT IS HELD CONSTANT.  Hamiltonians, standardisation, sentinel penalty, starts, budget counting
and RMSD readout are Sprint 20's, imported unchanged from `s20/qb2_lib.py` and `s20/qb2_opt.py`.
One budget unit = one evaluation of the objective at one continuous configuration, finite
difference probes INCLUDED.  The readout is ALWAYS the argmin over everything evaluated
(`Field.best_z`), on the BUILT CHAIN basis -- never a point cloud, never an average.
"""
from __future__ import annotations

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

from s20 import qb2_lib as L        # noqa: E402  -- Hamiltonians, subset, stats, io
from s20 import qb2_opt as OP       # noqa: E402  -- Field, arms, gradient estimators
from s15 import seed as SD          # noqa: E402
from s19 import qb_lib as QB        # noqa: E402
from core import quantum as Q       # noqa: E402

SALT = "s21qb"
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

MDE = 0.084


# ==================================================================== BLOCK E: the charts
class ThetaChart:
    """`x` such that `theta = scale * x`.  A LINEAR chart change of the same torus.

    `scale` is pure gauge for the PROBLEM (the map is a bijection and the objective is
    evaluated at `theta`), and pure step size for the OPTIMISER: a fixed-size step in `x` is a
    step of `scale` radians.  `scale` is restricted to exact powers of two so that
    `(theta/scale)*scale == theta` bit-exactly and the harness gate (F-E5) is meaningful.

    `do_wrap=False` disables the post-step re-wrap.  Wrapping is the identity on the objective
    (it is periodic) but NOT on a stateful optimiser: Adam's moments and SPSA's iterate jump by
    2*pi/scale when it fires.  `nelder`/`lbfgs` never call it -- scipy owns the iterate -- so
    that cell is a wrap-free control by construction.
    """

    def __init__(self, F, scale=1.0, do_wrap=True):
        self.F = F
        self.scale = float(scale)
        self.do_wrap = bool(do_wrap)
        self.n = int(F.n)                 # arms use d = 2*n; angle chart dim = 2*n_res
        self.t = F.t
        self.n_wrap_fired = 0

    def to_z(self, X):
        return self.scale * np.atleast_2d(np.asarray(X, float))

    def from_z(self, Z):
        return np.atleast_2d(np.asarray(Z, float)) / self.scale

    def wrap(self, x):
        if not self.do_wrap:
            return x
        z = self.scale * np.asarray(x, float)
        zw = L.wrap(z)
        self.n_wrap_fired += int(np.sum(np.abs(zw - z) > 1e-12))
        return zw / self.scale

    def restart(self, rng, k=1):
        return self.from_z(self.F.restart(rng, k))

    @property
    def left(self):
        return self.F.left

    def __call__(self, X):
        return self.F(self.to_z(X))

    def info(self):
        return {"chart": "theta", "scale": self.scale, "wrap": self.do_wrap,
                "dim": 2 * self.n, "n_wrap_fired": self.n_wrap_fired}


class EmbChart:
    """`x = r * (cos theta, sin theta)`, retracted by `arctan2`.

    `arctan2(y, x)` is invariant under `(x, y) -> (c*x, c*y)` for `c > 0`, so **`r` is exact
    gauge**: every radius names the same physical configuration and the same objective value.
    What `r` changes is the OPTIMISER's geometry -- a step of size `d` in the tangent direction
    is `d/r` radians -- and, because Adam and SPSA steps leave the circle, how fast the radius
    drifts, which is an implicit step-size schedule (F-E4).

    `renorm=True` retracts each `(cos, sin)` pair back to radius `r` after every step: the same
    chart, the same dimension, the same smoothness, no drift.
    """

    def __init__(self, F, r=1.0, renorm=False):
        self.F = F
        self.r = float(r)
        self.renorm = bool(renorm)
        self.nz = int(F.n)                # residues
        self.n = 2 * int(F.n)             # arms use d = 2*n; circle chart dim = 4*n_res
        self.t = F.t
        self.rad_min = np.inf
        self.rad_max = 0.0
        self.rad_last = float("nan")

    def _d(self):
        return 2 * self.nz

    def to_z(self, X):
        X = np.atleast_2d(np.asarray(X, float))
        d = self._d()
        return np.arctan2(X[:, d:], X[:, :d])

    def from_z(self, Z):
        Z = np.atleast_2d(np.asarray(Z, float))
        return np.concatenate([self.r * np.cos(Z), self.r * np.sin(Z)], axis=1)

    def wrap(self, x):
        x = np.asarray(x, float)
        d = self._d()
        rad = np.hypot(x[..., :d], x[..., d:])
        self.rad_min = min(self.rad_min, float(np.min(rad)))
        self.rad_max = max(self.rad_max, float(np.max(rad)))
        self.rad_last = float(np.median(rad))
        if not self.renorm:
            return x
        s = self.r / np.maximum(rad, 1e-12)
        out = x.copy()
        out[..., :d] *= s
        out[..., d:] *= s
        return out

    def restart(self, rng, k=1):
        return self.from_z(self.F.restart(rng, k))

    @property
    def left(self):
        return self.F.left

    def __call__(self, X):
        return self.F(self.to_z(X))

    def info(self):
        return {"chart": "emb", "r": self.r, "renorm": self.renorm, "dim": 4 * self.nz,
                "radius_min": None if not np.isfinite(self.rad_min) else self.rad_min,
                "radius_max": self.rad_max, "radius_last_med": self.rad_last}


# All eight parameterisations of PREREG section 1.  Every one is physically identical.
CHARTS = {
    "th":        lambda F: ThetaChart(F, 1.0, True),
    "th_nowrap": lambda F: ThetaChart(F, 1.0, False),
    "th_s05":    lambda F: ThetaChart(F, 0.5, True),
    "th_s2":     lambda F: ThetaChart(F, 2.0, True),
    "emb_r1":    lambda F: EmbChart(F, 1.0, False),
    "emb_r05":   lambda F: EmbChart(F, 0.5, False),
    "emb_r2":    lambda F: EmbChart(F, 2.0, False),
    "emb_norm":  lambda F: EmbChart(F, 1.0, True),
}


# ==================================================================== BLOCK A: the ansatz zoo
class SharedMPS:
    """`MPSAnsatz` with ONE RY angle per block, shared across every wire.

    `nblocks` parameters instead of `n*nblocks`.  The chain rule for a shared parameter is the
    SUM of the per-wire partials, which is exactly what broadcasting the scalar into the full
    angle vector and summing the returned gradient gives.  Bond dimension is unchanged
    (`2**layers`): sharing removes parameters, not entanglement.
    """

    def __init__(self, n, layers=2, final_ry=True, entangler="cnot"):
        self.a = Q.MPSAnsatz(n, layers=layers, final_ry=final_ry, entangler=entangler)
        self.n = n
        self.nblocks = self.a.nblocks
        self.chi = self.a.chi

    def n_params(self):
        return self.nblocks

    def _full(self, theta):
        return np.repeat(np.asarray(theta, float), self.n)

    def sample(self, theta, shots, rng):
        return self.a.sample(self._full(theta), shots, rng)

    def logp(self, theta, bits):
        return self.a.logp(self._full(theta), bits)

    def grad_logp(self, theta, bits, weights):
        g = self.a.grad_logp(self._full(theta), bits, weights)
        return np.asarray(g, float).reshape(self.nblocks, self.n).sum(1)

    def probs(self, theta):
        return self.a.probs(self._full(theta))

    def amplitudes(self, theta):
        return mps_amplitudes(self.a, self._full(theta))


class SVAnsatz:
    """Exact real statevector, `layers x (RY on every wire, CNOT chain, RING closure)`.

    A topology the OPEN MPS chain cannot reach: the ring gate is long-range for an open MPS and
    raises the bond dimension across every cut.  Simulated exactly (`2**n` amplitudes,
    `n <= 16` here), so `p_theta` and `grad log p` are ANALYTIC, not sampled.

    `grad log p(b) = (1/p(b)) dp(b)/dtheta_j`, and `p(b) = <psi|b><b|psi>` is an expectation of
    a projector, so the two-term PARAMETER-SHIFT rule for an RY generator is exact:
    `dp/dtheta_j = [p(theta + (pi/2) e_j) - p(theta - (pi/2) e_j)] / 2`.  One `probs_batch` call
    on the `2P` shifted parameter vectors gives every partial for every bitstring at once.
    """

    MAXN = 16

    def __init__(self, n, layers=3, ring=True):
        if n > self.MAXN:
            raise ValueError(f"SVAnsatz refuses n={n} > {self.MAXN} (2**n statevector)")
        self.c = Q.StatevectorCircuit(n, layers=layers, ring=ring)
        self.n = n
        self.layers = layers
        self.ring = ring
        self.chi = None                   # measured, not analytic

    def n_params(self):
        return self.c.n_params()

    def probs(self, theta):
        return self.c.probs(np.asarray(theta, float))

    def amplitudes(self, theta):
        return self.c.state(np.asarray(theta, float))

    @staticmethod
    def _idx(bits):
        b = np.asarray(bits, np.int64)
        w = (1 << np.arange(b.shape[1] - 1, -1, -1)).astype(np.int64)
        return b @ w

    def sample(self, theta, shots, rng):
        p = self.probs(theta)
        k = rng.choice(len(p), size=shots, p=p / p.sum())
        n = self.n
        return ((k[:, None] >> np.arange(n - 1, -1, -1)[None, :]) & 1).astype(np.int64)

    def logp(self, theta, bits):
        p = self.probs(theta)
        return np.log(np.maximum(p[self._idx(bits)], 1e-300))

    def grad_logp(self, theta, bits, weights):
        th = np.asarray(theta, float)
        P = th.size
        TH = self.c._shift_grid(th, np.pi / 2.0)
        pr = self.c.probs_batch(TH)                     # (2P, dim)
        dp = 0.5 * (pr[0::2] - pr[1::2])                # (P, dim)
        k = self._idx(bits)
        p0 = np.maximum(self.probs(th)[k], 1e-300)      # (B,)
        gl = dp[:, k] / p0[None, :]                     # (P, B)
        return gl @ np.asarray(weights, float)


def make_ansatz(kind, n):
    """The zoo.  Returns (ansatz, static_info)."""
    if kind == "prod":
        a = Q.MPSAnsatz(n, layers=2, final_ry=True, entangler="none")
    elif kind == "mps_L1":
        a = Q.MPSAnsatz(n, layers=1, final_ry=True, entangler="cnot")
    elif kind == "mps_L2":
        a = Q.MPSAnsatz(n, layers=2, final_ry=True, entangler="cnot")
    elif kind == "mps_L3":
        a = Q.MPSAnsatz(n, layers=3, final_ry=True, entangler="cnot")
    elif kind == "mps_L4":
        a = Q.MPSAnsatz(n, layers=4, final_ry=True, entangler="cnot")
    elif kind == "mps_L2_nofinal":
        a = Q.MPSAnsatz(n, layers=2, final_ry=False, entangler="cnot")
    elif kind == "share_L2":
        a = SharedMPS(n, layers=2, final_ry=True, entangler="cnot")
    elif kind == "sv_ring_L3":
        a = SVAnsatz(n, layers=3, ring=True)
    else:
        raise ValueError(kind)
    info = {"kind": kind, "n_params": int(a.n_params()),
            "chi_analytic": (int(a.chi) if getattr(a, "chi", None) else None)}
    return a, info


def mps_amplitudes(ans, theta):
    """Exact amplitude vector of an `MPSAnsatz`, `amp(x) = e0^T prod_q A[q][:, x_q, :] e0`.

    Contracted left to right over the physical indices; `(2**n,)` for `n <= 16`.  Used ONLY to
    measure the Schmidt rank -- never inside a budgeted arm.
    """
    A = np.asarray(ans.build(np.asarray(theta, float)), float)   # (n, chi, 2, chi)
    n, chi = A.shape[0], A.shape[1]
    v = np.zeros((1, chi))
    v[0, 0] = 1.0
    for q in range(n):
        v = np.einsum("kc,csr->ksr", v, A[q]).reshape(-1, chi)
    return v[:, 0]


def schmidt_ranks(vec, n, tol=1e-10):
    """Schmidt rank of a `2**n` real amplitude (or probability) tensor at EVERY bipartition.

    `rank_k = #{ sigma : sigma > tol * sigma_max }` of the `2**k x 2**(n-k)` unfolding.  This is
    the bond dimension the state actually realises, as opposed to the one the circuit's
    topology permits.
    """
    v = np.asarray(vec, float).ravel()
    out = []
    for k in range(1, n):
        M = v.reshape(1 << k, 1 << (n - k))
        s = np.linalg.svd(M, compute_uv=False)
        smax = s.max() if s.size else 0.0
        out.append(int((s > tol * max(smax, 1e-300)).sum()))
    return out


def ansatz_geometry(ans, theta, n):
    """The measured expressivity panel: realised bond dimension, entropy, participation."""
    rec = {}
    try:
        if hasattr(ans, "amplitudes"):
            amp = ans.amplitudes(theta)
        else:
            amp = mps_amplitudes(ans, theta)
        amp = np.asarray(amp, float)
        nrm = np.linalg.norm(amp)
        if nrm > 0:
            amp = amp / nrm
        rec["chi_amp_max"] = int(max(schmidt_ranks(amp, n)))
        rec["chi_amp"] = schmidt_ranks(amp, n)
        p = amp ** 2
    except Exception as exc:                                   # pragma: no cover
        rec["chi_amp_max"] = None
        rec["chi_err"] = repr(exc)[:120]
        p = np.asarray(ans.probs(theta), float)
        p = p / p.sum()
    p = np.maximum(np.asarray(p, float), 0)
    p = p / max(p.sum(), 1e-300)
    try:
        rec["chi_prob_max"] = int(max(schmidt_ranks(p, n)))
    except Exception:
        rec["chi_prob_max"] = None
    pos = p[p > 0]
    rec["entropy_bits"] = float(-(pos * np.log2(pos)).sum())
    rec["max_prob"] = float(p.max())
    rec["latent_ess"] = float(1.0 / (p ** 2).sum())
    return rec


# ==================================================================== BLOCK Q: optimisers
class SGD:
    def __init__(self, n, lr=0.15):
        self.lr = float(lr)

    def step(self, x, g):
        return np.asarray(x, float) - self.lr * np.asarray(g, float)


class QNG:
    """Fisher-preconditioned natural gradient on the SAMPLER's own information metric.

    `F = E_{b~p}[ grad log p(b) grad log p(b)^T ]`, estimated on the SAME shots the CVaR
    gradient uses, so it costs ZERO extra objective evaluations.  The update is
    `x <- x - lr * (F + lam*I)^{-1} g`, with `lam` a fixed Levenberg damping.

    THIS IS THE CLASSICAL FISHER OF THE MEASUREMENT DISTRIBUTION, not the Fubini-Study metric
    of the state.  For an objective that depends on the state only through
    `p_theta(b) = |<b|psi>|^2` -- which is exactly our case, the Hamiltonian being diagonal in
    the computational basis -- the classical Fisher is the metric the KL divergence between
    successive sampling laws induces, and it is the one a natural-gradient step should use.
    The two coincide on the diagonal block and differ by the off-diagonal quantum geometric
    tensor; that difference is recorded, not glossed.

    `diag=True` keeps only the diagonal (a per-parameter rescale, `P` flops instead of a
    `P x P` solve) -- the cheap approximation, run as its own arm.
    """

    def __init__(self, n, lr=0.15, lam=1e-2, diag=False):
        self.lr = float(lr)
        self.lam = float(lam)
        self.diag = bool(diag)
        self.cond = []

    def precondition(self, g, S):
        """`S` is the (shots, P) matrix of per-sample scores."""
        S = np.asarray(S, float)
        P = S.shape[1]
        if self.diag:
            d = (S ** 2).mean(0)
            return g / (d + self.lam)
        F = (S.T @ S) / max(1, S.shape[0])
        w = np.linalg.eigvalsh((F + F.T) / 2.0)
        a = np.abs(w)
        self.cond.append(float(a.max() / max(a.min(), 1e-300)) if a.size else float("nan"))
        return np.linalg.solve(F + self.lam * np.eye(P), g)

    def step(self, x, g):
        return np.asarray(x, float) - self.lr * np.asarray(g, float)


def _scores(ans, theta, bits):
    """`(shots, P)` per-sample score vectors `grad log p(b_i)`, one autograd pass per sample.

    Exact for `SVAnsatz` (analytic parameter shift, one batched call).  For the MPS family the
    accumulator only returns weighted sums, so the per-sample rows are obtained by calling it
    with `shots` one-hot weight vectors -- correct, and the reason `qng` is run at a modest
    shot count.
    """
    if isinstance(ans, SVAnsatz):
        th = np.asarray(theta, float)
        TH = ans.c._shift_grid(th, np.pi / 2.0)
        pr = ans.c.probs_batch(TH)
        dp = 0.5 * (pr[0::2] - pr[1::2])
        k = ans._idx(bits)
        p0 = np.maximum(ans.probs(th)[k], 1e-300)
        return (dp[:, k] / p0[None, :]).T
    B = bits.shape[0]
    P = ans.n_params()
    S = np.empty((B, P))
    eye = np.eye(B)
    for i in range(B):
        S[i] = Q.grad_logp_weighted(ans, theta, bits, eye[i])
    return S


# ==================================================================== the VQE arm
def arm_vqe2(F, z0, rng, ansatz="mps_L2", opt="adam", alpha=0.25, shots=64, lr=0.15,
             train=True, lam=1e-2, need_scores=False, geom=True):
    """GENUINE CVaR-VQE over the CONTINUOUS torsion representation, configurable ansatz and
    optimiser.

    The bitstring selects a conformer basin per residue and the ANGLE is drawn CONTINUOUSLY
    from that basin's von Mises component, so the sampler's law is a continuous density with
    full support on the torus.  Gradient: `core.quantum.cvar_gradient(..., baseline="const")` --
    the CORRECTED estimator, never the recorded `tail` defect.

    READOUT: `Field.best_z`, the argmin over everything evaluated.  CVaR is the TRAINING
    objective only (BRIEF section 2, as corrected 2026-09-07).
    """
    t = F.t
    n = int(t["n"])
    ans, ainfo = make_ansatz(ansatz, n)
    P = ans.n_params()
    th = np.pi / 2 + rng.normal(0.0, 0.8, P)
    th0 = th.copy()
    if opt == "adam":
        O = Q.Adam(P, lr=lr)
    elif opt == "sgd":
        O = SGD(P, lr=lr)
    elif opt == "qng":
        O = QNG(P, lr=lr, lam=lam, diag=False)
    elif opt == "qng_diag":
        O = QNG(P, lr=lr, lam=lam, diag=True)
    elif opt == "spsa_ansatz":
        O = None
    else:
        raise ValueError(opt)
    want_scores = need_scores or opt in ("qng", "qng_diag")
    ess, gnorm, cvals, fcond = [], [], [], []
    it = 0
    while F.left >= shots:
        it += 1
        if opt == "spsa_ansatz":
            # SPSA on the ANSATZ parameters: one symmetric-Bernoulli perturbation, two
            # sampled CVaR evaluations, `shots//2` shots each so the budget matches exactly.
            half = max(8, shots // 2)
            if F.left < 2 * half:
                break
            d = rng.choice([-1.0, 1.0], size=P)
            ck = 0.15 / it ** 0.101
            ak = 0.30 / (it + 5.0) ** 0.602
            vals = []
            for sgn in (+1.0, -1.0):
                b = ans.sample(th + sgn * ck * d, half, rng)
                phi, psi = QB.draw_from_basins(b, t["mu"], t["kap"], rng)
                e = F(L.pack(phi, psi))
                if e.size < half:
                    vals = []
                    break
                vals.append(Q.cvar(e, alpha)[0])
            if len(vals) < 2:
                break
            g = (vals[0] - vals[1]) / (2 * ck) * (1.0 / d)
            cvals.append(float(0.5 * (vals[0] + vals[1])))
            gnorm.append(float(np.linalg.norm(g)))
            if train:
                th = th - ak * g
            continue
        bits = ans.sample(th, shots, rng)
        phi, psi = QB.draw_from_basins(bits, t["mu"], t["kap"], rng)
        e = F(L.pack(phi, psi))
        if e.size < shots:
            break
        val, q, tail = Q.cvar(e, alpha)
        w = np.zeros(len(e)); w[tail] = -(q - e[tail]) / (alpha * len(e))
        w = w - w.mean()
        aw = np.abs(w)
        ess.append(float(aw.sum() ** 2 / max(1e-300, (aw ** 2).sum())))
        cvals.append(float(val))
        g, _v = Q.cvar_gradient(ans, th, bits, e, alpha, baseline="const")
        g = np.asarray(g, float)
        gnorm.append(float(np.linalg.norm(g)))
        if want_scores:
            S = _scores(ans, th, bits)
            g = O.precondition(g, S)
            if getattr(O, "cond", None):
                fcond = O.cond
        if train:
            th = O.step(th, g)
    out = {"iters": it, "ansatz": ansatz, "opt": opt, "alpha": alpha, "shots": shots,
           "n_params": int(P), **{f"a_{k}": v for k, v in ainfo.items() if k != "kind"},
           "ess_mean": float(np.mean(ess)) if ess else float("nan"),
           "gnorm_mean": float(np.mean(gnorm)) if gnorm else float("nan"),
           "gnorm_sd": float(np.std(gnorm)) if gnorm else float("nan"),
           "gnorm_first": float(gnorm[0]) if gnorm else float("nan"),
           "gnorm_last": float(gnorm[-1]) if gnorm else float("nan"),
           "cvar_first": cvals[0] if cvals else float("nan"),
           "cvar_last": cvals[-1] if cvals else float("nan"),
           "param_disp": float(np.linalg.norm(th - th0)),
           "fisher_cond_med": float(np.median(fcond)) if fcond else float("nan")}
    if geom and n <= 16:
        try:
            out.update({f"g0_{k}": v for k, v in ansatz_geometry(ans, th0, n).items()
                        if not isinstance(v, list)})
            gt = ansatz_geometry(ans, th, n)
            out.update({f"g1_{k}": v for k, v in gt.items() if not isinstance(v, list)})
            out["chi_amp_profile"] = gt.get("chi_amp")
        except Exception as exc:                                # pragma: no cover
            out["geom_err"] = repr(exc)[:120]
    return out


# ==================================================================== stats helpers
paired_ci = L.paired_ci
boot_mean_ci = L.boot_mean_ci
write = L.write


def write21(name, obj, complete=True):
    import json
    import time
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, path)
            break
        except PermissionError:
            time.sleep(0.4)
    if complete:
        with open(path.replace(".json", "") + "_COMPLETE", "w") as fh:
            fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return path


def read21(name):
    import json
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    with open(path) as fh:
        return json.load(fh)


def sign_test_p(k, n):
    """Two-sided exact sign test, `k` successes out of `n`."""
    from math import comb
    k = min(k, n - k)
    return float(min(1.0, 2 * sum(comb(n, i) for i in range(k + 1)) / 2 ** n))
