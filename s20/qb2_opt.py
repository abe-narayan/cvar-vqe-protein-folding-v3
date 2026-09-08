"""SPRINT 20 / WORKSTREAM B -- the optimiser battery and the CVaR-VQE arms.

EVERY arm here sees the SAME standardised field  Ehat = (E - mu_ref)/sd_ref  (the free,
unbudgeted calibration on the target's own K=500 retrieval pool), the SAME starts, the SAME
seeds and the SAME evaluation budget.  One budget unit = one evaluation of the objective at one
continuous configuration, **finite-difference probes included**.  Nothing here is given a gain,
temperature or step size that depends on which Hamiltonian it is running against.

`best_of_N` is MANDATORY and is the only honest no-optimisation control: N = B i.i.d. draws from
the same von Mises basin mixtures every other arm is initialised from, at exactly the budget.
No arm in this lane is ever compared against an initialisation mean.

THE NON-FINITE SENTINEL.  `core.amber` returns +inf for a sterically collapsed structure.  Every
arm maps it to the SAME finite penalty `PEN` (10 standardised units above the reference maximum)
so that no optimiser gets an accidental advantage from how it happens to handle inf.  The number
of times that substitution fires is counted and reported per arm: a guard that never fires is not
evidence (BRIEF section 6 rule 3).
"""
from __future__ import annotations

import numpy as np

from s20 import qb2_lib as L
from s15 import seed as SD
from s19 import qb_lib as QB
from core import quantum as Q

SALT = L.SALT


class Field:
    """The standardised, sentinel-guarded, budget-counted view of one Hamiltonian."""

    def __init__(self, ham, tgt, pen_std=10.0):
        self.h = ham
        self.t = tgt
        self.n = int(tgt["n"])
        self.n_nonfinite = 0
        self.best = np.inf
        self.best_z = None
        self.PEN = None
        self.pen_std = pen_std

    def set_pen(self, ref_std_max):
        self.PEN = float(ref_std_max + self.pen_std)

    @staticmethod
    def wrap(z):
        """The angle is periodic, so wrapping is the identity on the objective.  The embedded
        field (S3) overrides it with a true identity, because its coordinates are not angles."""
        return L.wrap(z)

    def restart(self, rng, k=1):
        """A fresh start in the SPACE THIS FIELD IS PARAMETERISED IN.  Overridden by the
        embedded field (S3) so that a restart is matched in the space the operator actually
        works in (BRIEF section 6 rule 1)."""
        return L.basin_starts(self.t, k, rng)

    @property
    def left(self):
        return self.h.left

    def __call__(self, Z):
        Z = np.atleast_2d(np.asarray(Z, float))
        phi, psi = L.unpack(Z, self.n)
        e = self.h(phi, psi)
        if e.size == 0:
            return e
        s = self.h.std(e)
        bad = ~np.isfinite(s)
        if bad.any():
            self.n_nonfinite += int(bad.sum())
            s = np.where(bad, self.PEN, s)
        k = int(np.argmin(s))
        if s[k] < self.best:
            self.best = float(s[k]); self.best_z = Z[k].copy()
        return s


# =========================================================== gradient estimators
def _grad_fd(F, z, h=L.FD_H):
    """Central FD gradient.  Costs 2*2n budget units."""
    d = 2 * F.n
    if F.left < 2 * d:
        return None
    Z = np.repeat(np.atleast_2d(z), 2 * d, axis=0)
    for k in range(d):
        Z[2 * k, k] += h; Z[2 * k + 1, k] -= h
    e = F(Z)
    if len(e) < 2 * d:
        return None
    return (e[0::2] - e[1::2]) / (2 * h)


def _grad_spsa(F, z, c, rng):
    """SPSA: ONE symmetric-Bernoulli perturbation, 2 budget units."""
    if F.left < 2:
        return None, None
    d = 2 * F.n
    delta = rng.choice([-1.0, 1.0], size=d)
    e = F(np.stack([z + c * delta, z - c * delta]))
    if len(e) < 2:
        return None, None
    return (e[0] - e[1]) / (2 * c) * (1.0 / delta), delta


# =========================================================== the arms
def arm_spsa(F, z0, rng, a=0.20, c=L.FD_H, A_frac=0.1, track_quality=True):
    """SPSA with the standard decaying gains.  `track_quality` measures, at 3 checkpoints,
    the cosine between the SPSA estimate and the central-FD gradient at the SAME point -- the
    separation the brief asks for between gradient-estimation noise and landscape difficulty.
    The quality probes are BUDGETED like everything else."""
    z = np.asarray(z0, float).copy()
    it = 0
    max_it = max(1, F.left // 2)
    A = max(1.0, A_frac * max_it)
    cos, unorm = [], []
    checks = {int(max_it * f) for f in (0.05, 0.4, 0.8)}
    while F.left >= 2:
        it += 1
        ak = a / (it + A) ** 0.602
        ck = c / it ** 0.101
        g, _d = _grad_spsa(F, z, ck, rng)
        if g is None:
            break
        if track_quality and it in checks and F.left >= 4 * F.n + 2:
            gt = _grad_fd(F, z)
            if gt is not None and np.linalg.norm(gt) > 0 and np.linalg.norm(g) > 0:
                cos.append(float(g @ gt / (np.linalg.norm(g) * np.linalg.norm(gt))))
        step = ak * g
        unorm.append(float(np.linalg.norm(step)))
        z = F.wrap(z - step)
    return {"iters": it, "spsa_cos_fd": cos,
            "update_norm_med": float(np.median(unorm)) if unorm else float("nan"),
            "update_norm_p95": float(np.percentile(unorm, 95)) if unorm else float("nan")}


def arm_adam(F, z0, rng, lr=0.10):
    z = np.asarray(z0, float).copy()
    opt = Q.Adam(len(z), lr=lr)
    it, unorm = 0, []
    while F.left >= 4 * F.n:
        g = _grad_fd(F, z)
        if g is None:
            break
        it += 1
        zn = opt.step(z, g)
        unorm.append(float(np.linalg.norm(zn - z)))
        z = F.wrap(zn)
    return {"iters": it, "update_norm_med": float(np.median(unorm)) if unorm else float("nan"),
            "update_norm_p95": float(np.percentile(unorm, 95)) if unorm else float("nan")}


def arm_mom(F, z0, rng, lr=0.02, beta=0.9):
    z = np.asarray(z0, float).copy()
    v = np.zeros_like(z)
    it, unorm = 0, []
    while F.left >= 4 * F.n:
        g = _grad_fd(F, z)
        if g is None:
            break
        it += 1
        v = beta * v + g
        unorm.append(float(np.linalg.norm(lr * v)))
        z = F.wrap(z - lr * v)
    return {"iters": it, "update_norm_med": float(np.median(unorm)) if unorm else float("nan"),
            "update_norm_p95": float(np.percentile(unorm, 95)) if unorm else float("nan")}


def arm_lbfgs(F, z0, rng):
    """L-BFGS-B with the SAME central-FD gradient; every probe counted."""
    from scipy.optimize import minimize
    z = np.asarray(z0, float).copy()
    n_starts = 0

    def fun(x):
        if F.left < 1:
            raise _Stop()
        e = F(x[None])
        if e.size == 0:
            raise _Stop()
        return float(e[0])

    def jac(x):
        g = _grad_fd(F, x)
        if g is None:
            raise _Stop()
        return g

    while F.left > 4 * F.n + 4:
        n_starts += 1
        try:
            minimize(fun, z, jac=jac, method="L-BFGS-B",
                     options={"maxiter": 200, "maxfun": 10 ** 6})
        except _Stop:
            break
        # multi-start restart from a fresh basin draw, as the incumbent optimiser does
        z = F.restart(rng, 1)[0]
    return {"restarts": n_starts}


def arm_nelder(F, z0, rng):
    return _scipy_df(F, z0, rng, "Nelder-Mead")


def arm_powell(F, z0, rng):
    return _scipy_df(F, z0, rng, "Powell")


class _Stop(Exception):
    pass


def _scipy_df(F, z0, rng, method):
    from scipy.optimize import minimize
    z = np.asarray(z0, float).copy()
    n_starts = 0

    def fun(x):
        if F.left < 1:
            raise _Stop()
        e = F(x[None])
        if e.size == 0:
            raise _Stop()
        return float(e[0])

    while F.left > 4:
        n_starts += 1
        try:
            minimize(fun, z, method=method,
                     options={"maxiter": 10 ** 6, "maxfev": 10 ** 6})
        except _Stop:
            break
        z = F.restart(rng, 1)[0]
    return {"restarts": n_starts}


def arm_bestofn(F, z0, rng):
    """MANDATORY.  N = B i.i.d. draws from the same basin mixtures.  Zero optimisation."""
    n_batch = 0
    while F.left > 0:
        k = min(64, F.left)
        F(F.restart(rng, k))
        n_batch += 1
    return {"batches": n_batch}


def arm_metro(F, z0, rng, Tfrac=0.10, anneal=False, chains=8):
    """Continuous single-residue Metropolis / simulated annealing in torsion space.
    The move class is matched to every other arm's: one residue's (phi, psi) redrawn from its
    own von Mises basin mixture.  T is set from the FREE reference sd, which is 1.0 in
    standardised units -- so the temperature is IDENTICAL across Hamiltonians by construction."""
    t = F.t
    n = F.n
    Z = F.restart(rng, chains)
    cur = F(Z)
    if cur.size < chains:
        return {"chains": chains, "steps": 0}
    steps, acc_n = 0, 0
    ar = np.arange(chains)
    budget0 = F.h.budget
    while F.left >= chains:
        steps += 1
        frac = min(1.0, (budget0 - F.left) / max(1, budget0))
        T = Tfrac * (4.0 * (1 / 64.0) ** frac if anneal else 1.0)
        i = rng.integers(0, n, chains)
        c = (rng.random(chains) < t["wmarg"][i, 1]).astype(np.int64)
        Zn = Z.copy()
        Zn[ar, i] = L.wrap(rng.vonmises(t["mu"][i, c, 0], t["kap"][i, c, 0]))
        Zn[ar, n + i] = L.wrap(rng.vonmises(t["mu"][i, c, 1], t["kap"][i, c, 1]))
        ev = F(Zn)
        if ev.size < chains:
            break
        de = ev - cur
        ok = (de <= 0) | (rng.random(chains) < np.exp(np.clip(-de / max(T, 1e-12), -700, 0)))
        acc_n += int(ok.sum())
        Z[ok] = Zn[ok]; cur[ok] = ev[ok]
    return {"chains": chains, "steps": steps,
            "accept_rate": acc_n / max(1, steps * chains)}


def arm_anneal(F, z0, rng):
    return arm_metro(F, z0, rng, Tfrac=0.10, anneal=True)


def arm_greedy(F, z0, rng, chains=8):
    """Greedy coordinate descent: identical proposal class to `metro`, accept only strict
    improvements (T = 0).  The zero-temperature end of the same thermostat family."""
    return arm_metro(F, z0, rng, Tfrac=1e-12, anneal=False, chains=chains)


ARMS = {
    "best_of_N": arm_bestofn,      # MANDATORY no-optimisation control
    "spsa": arm_spsa,
    "adam_fd": arm_adam,
    "mom_fd": arm_mom,
    "lbfgs_fd": arm_lbfgs,
    "nelder": arm_nelder,
    "powell": arm_powell,
    "metro": arm_metro,
    "anneal": arm_anneal,
    "greedy": arm_greedy,
}


# =========================================================== the CVaR-VQE arm (S2)
def arm_vqe(F, z0, rng, alpha=0.25, shots=64, lr=0.15, entangler="cnot", train=True):
    """GENUINE CVaR-VQE over the CONTINUOUS torsion representation.

    `core.quantum.MPSAnsatz` (RY + CNOT chain + final RY), simulated exactly; its bitstring
    selects a conformer basin per residue and the ANGLE is drawn CONTINUOUSLY from that basin's
    von Mises component, so the sampler's law is a continuous density with full support on the
    torus.  Gradient: `core.quantum.cvar_gradient(..., baseline="const")` -- the CORRECTED
    estimator, never the recorded `tail` defect.

    Reports the CVaR-specific panel the brief asks for: effective sample size of the tail,
    latent entropy, gradient variance, and the bias of the alpha-tail gradient against the
    alpha = 1 full-mean gradient computed on the SAME batch (free, no extra evaluations).
    """
    t = F.t
    n = F.n
    an = Q.MPSAnsatz(n, layers=2, final_ry=True, entangler=entangler)
    th = np.pi / 2 + rng.normal(0.0, 0.8, an.n_params())
    th0 = th.copy()
    opt = Q.Adam(an.n_params(), lr=lr)
    ess, gvar, gbias, cvals = [], [], [], []
    it = 0
    while F.left >= shots:
        it += 1
        bits = an.sample(th, shots, rng)
        phi, psi = QB.draw_from_basins(bits, t["mu"], t["kap"], rng)
        e = F(L.pack(phi, psi))
        if e.size < shots:
            break
        # -- CVaR panel (all free: no objective evaluation)
        val, q, tail = Q.cvar(e, alpha)
        w = np.zeros(len(e)); w[tail] = -(q - e[tail]) / (alpha * len(e))
        w = w - w.mean()
        aw = np.abs(w)
        ess.append(float(aw.sum() ** 2 / max(1e-300, (aw ** 2).sum())))
        cvals.append(float(val))
        g, _v = Q.cvar_gradient(an, th, bits, e, alpha, baseline="const")
        g1, _v1 = Q.cvar_gradient(an, th, bits, e, 1.0, baseline="const")
        ng, n1 = np.linalg.norm(g), np.linalg.norm(g1)
        if ng > 0 and n1 > 0:
            gbias.append(float(g @ g1 / (ng * n1)))
        gvar.append(float(ng))
        if train:
            th = opt.step(th, g)
    out = {"iters": it, "ess_mean": float(np.mean(ess)) if ess else float("nan"),
           "ess_frac": float(np.mean(ess) / shots) if ess else float("nan"),
           "gnorm_mean": float(np.mean(gvar)) if gvar else float("nan"),
           "gnorm_sd": float(np.std(gvar)) if gvar else float("nan"),
           "cos_grad_vs_alpha1": float(np.mean(gbias)) if gbias else float("nan"),
           "cvar_first": cvals[0] if cvals else float("nan"),
           "cvar_last": cvals[-1] if cvals else float("nan"),
           "param_disp": float(np.linalg.norm(th - th0))}
    if n <= 18:
        allb = Q.all_bitstrings(n)
        p = np.exp(np.asarray(an.logp(th, allb), float))
        p = np.maximum(p, 0); p = p / p.sum()
        pos = p[p > 0]
        out["latent_entropy_bits"] = float(-(pos * np.log2(pos)).sum())
        out["latent_ess"] = float(1.0 / (p ** 2).sum())
        out["latent_max_prob"] = float(p.max())
    return out
