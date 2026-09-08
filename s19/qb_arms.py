"""SPRINT 19 / AGENT B -- the arms.  Continuous-torsion samplers at MATCHED budget.

Every arm emits a SET of continuous torsion configurations and is scored by the identical
`qb_lib.readout`: generation ceiling, the (M, D) plane, selection ceiling, realised -- the
four quantities `s19/BRIEF.md` section 6 forbids collapsing.

THE ACCOUNTING.  One objective evaluation = one read of `E` at one continuous configuration
(`qb_lib.Obj.__call__`).  Circuit simulation, latent sampling, grad-log-p and the von Mises
draws are free.  Arms that read `E` zero times are labelled with their zero.

THE ONE ASYMMETRY, DECLARED: `c_metro*` and `c_anneal` calibrate their temperature from an
UNBUDGETED marginal batch (`Obj.raw`).  That free calibration is given to the CLASSICAL arms
and withheld from the quantum ones, i.e. it is conservative against this lane's own
hypothesis.
"""
from __future__ import annotations

import numpy as np

from s19 import qb_lib as L
from s15 import seed as SD
from core import quantum as Q

SALT = L.SALT
ALPHAS = (0.05, 0.25, 1.00)


# =============================================================== latent samplers
def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))


class ChainLatent:
    """CLASSICAL learned proposal: a first-order Markov chain over the latent b.

    `p(b) = Bern(sig(a)) * prod_{i>=1} Bern(sig(u[i, b_{i-1}]))`.  This is the classical
    model with EXACTLY the nearest-neighbour correlation structure a CNOT chain supplies,
    trained by the IDENTICAL CVaR score-function gradient at the IDENTICAL budget.  It is
    the decisive control: if the circuit's advantage is "correlated basin moves", this has
    them and costs nothing.
    """

    def __init__(self, n, rng, scale=0.3):
        self.n = int(n)
        self.a = float(rng.normal(0, scale))
        self.u = rng.normal(0, scale, (n, 2))

    def params(self):
        return np.concatenate([[self.a], self.u.reshape(-1)])

    def set_params(self, p):
        self.a = float(p[0]); self.u = p[1:].reshape(self.n, 2)

    def n_params(self):
        return 1 + 2 * self.n

    def sample(self, shots, rng):
        b = np.zeros((shots, self.n), np.int64)
        b[:, 0] = (rng.random(shots) < _sigmoid(self.a)).astype(np.int64)
        for i in range(1, self.n):
            p = _sigmoid(self.u[i, b[:, i - 1]])
            b[:, i] = (rng.random(shots) < p).astype(np.int64)
        return b

    def grad_logp(self, bits, w):
        """sum_s w_s grad log p(b_s).  Bernoulli logit gradient: x - sigmoid(u)."""
        bits = np.asarray(bits, np.int64)
        w = np.asarray(w, float)
        g = np.zeros(self.n_params())
        g[0] = float((w * (bits[:, 0] - _sigmoid(self.a))).sum())
        gu = np.zeros((self.n, 2))
        for i in range(1, self.n):
            prev = bits[:, i - 1]
            r = w * (bits[:, i] - _sigmoid(self.u[i, prev]))
            np.add.at(gu, (i, prev), r)
        g[1:] = gu.reshape(-1)
        return g


def _cvar_weights(e, alpha, baseline="const"):
    """The same envelope-theorem weights `core.quantum.cvar_gradient` uses, exposed so the
    classical latent models are trained by the IDENTICAL estimator, not a lookalike."""
    e = np.asarray(e, float)
    val, q, tail = Q.cvar(e, alpha)
    N = len(e)
    w = np.zeros(N)
    w[tail] = -(q - e[tail]) / (alpha * N)
    if baseline == "const":
        w = w - w.mean()
    return w, float(val)


# =============================================================== the quantum arm
def run_circuit(tgt, obj, alpha, seed, ansatz="mps2f", train=True, anneal=None,
                shots=L.SHOTS, lr=L.LR):
    """GENUINE CVaR-VQE over CONTINUOUS torsions.

    The circuit is `core.quantum.MPSAnsatz` (RY + CNOT chain + final RY), simulated exactly.
    Its bitstring selects a conformer basin per residue; the ANGLE is then drawn continuously
    from that basin's von Mises component.  The gradient is
    `core.quantum.cvar_gradient(..., baseline="const")` -- the corrected estimator.

    `train=False` is the MANDATORY best-of-N control: the same circuit at theta_0, the same
    number of samples, ZERO gradient steps.  Never an initialisation mean.
    """
    n = tgt["n"]
    an = Q.MPSAnsatz(n, layers=2, final_ry=True,
                     entangler=("none" if ansatz.endswith("n") else "cnot"))
    # NOTE the seed key deliberately EXCLUDES `train`: the untrained best-of-N control
    # therefore starts from the IDENTICAL theta_0 and draws the IDENTICAL first batch as
    # its trained twin, so `trained - untrained` isolates the gradient steps and nothing
    # else.  It is never compared against an initialisation mean.
    rng = SD.stable_rng(tgt["pdb"], seed, ansatz, alpha, salt=SALT)
    th = np.pi / 2 + rng.normal(0.0, 0.8, an.n_params())
    th0 = th.copy()
    opt = Q.Adam(an.n_params(), lr=lr)
    iters = max(1, obj.budget // shots)
    cvals = []
    for t in range(iters):
        if obj.left < shots:
            break
        a = alpha if anneal is None else Q.alpha_schedule(t / max(1, iters - 1), *anneal)
        bits = an.sample(th, shots, rng)
        phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
        e = obj(phi, psi)
        if e.size < shots:
            break
        if train:
            g, val = Q.cvar_gradient(an, th, bits, e, a, baseline="const")
            th = opt.step(th, g)
            cvals.append(float(val))
    # ---- native-free diagnostics of the LEARNED LATENT LAW (exact, n <= 20)
    diag = {"iters": len(cvals), "cvar_traj": cvals,
            "param_disp": float(np.linalg.norm(th - th0))}
    if n <= 18:
        allb = Q.all_bitstrings(n)
        p = np.exp(np.asarray(an.logp(th, allb), float))
        p = np.maximum(p, 0); p = p / p.sum()
        pos = p[p > 0]
        diag["entropy_bits"] = float(-(pos * np.log2(pos)).sum())
        diag["max_prob"] = float(p.max())
        diag["ess"] = float(1.0 / (p ** 2).sum())
        # nearest-neighbour latent mutual information: what the CNOTs actually buy
        mi = []
        for i in range(n - 1):
            j2 = np.stack([allb[:, i], allb[:, i + 1]], 1)
            P = np.zeros((2, 2))
            np.add.at(P, (j2[:, 0], j2[:, 1]), p)
            P = np.maximum(P, 1e-300)
            m = float((P * np.log2(P / (P.sum(1, keepdims=True) * P.sum(0, keepdims=True))))
                      .sum())
            mi.append(m)
        diag["nn_mi_bits_mean"] = float(np.mean(mi)) if mi else 0.0
        diag["nn_mi_bits_max"] = float(np.max(mi)) if mi else 0.0
    return diag


# =============================================================== classical latents
def run_chain(tgt, obj, alpha, seed, shots=L.SHOTS, lr=L.LR):
    rng = SD.stable_rng(tgt["pdb"], seed, "chain", alpha, salt=SALT)
    mdl = ChainLatent(tgt["n"], rng)
    opt = Q.Adam(mdl.n_params(), lr=lr)
    iters = max(1, obj.budget // shots)
    for _ in range(iters):
        if obj.left < shots:
            break
        bits = mdl.sample(shots, rng)
        phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
        e = obj(phi, psi)
        if e.size < shots:
            break
        w, _v = _cvar_weights(e, alpha)
        mdl.set_params(opt.step(mdl.params(), mdl.grad_logp(bits, w)))
    return {"iters": iters}


def run_cem(tgt, obj, alpha, seed, shots=L.SHOTS, smooth=0.7):
    """Cross-entropy method on a product latent: the standard classical learned proposal."""
    rng = SD.stable_rng(tgt["pdb"], seed, "cem", alpha, salt=SALT)
    n = tgt["n"]
    p = np.full(n, 0.5)
    iters = max(1, obj.budget // shots)
    for _ in range(iters):
        if obj.left < shots:
            break
        bits = (rng.random((shots, n)) < p[None, :]).astype(np.int64)
        phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
        e = obj(phi, psi)
        if e.size < shots:
            break
        k = max(2, int(round(alpha * len(e))))
        elite = bits[np.argsort(e, kind="stable")[:k]]
        p = smooth * p + (1 - smooth) * np.clip(elite.mean(0), 0.02, 0.98)
    return {"iters": iters}


# =============================================================== thermostats
def _scale(tgt, obj, rng, m=256):
    """FREE (unbudgeted) calibration of the objective's scale, given only to classical arms."""
    bits = (rng.random((m, tgt["n"])) < tgt["wmarg"][:, 1][None, :]).astype(np.int64)
    phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
    e, _ = obj.raw(phi, psi)
    return float(np.std(e)), phi, psi


def run_metro(tgt, obj, seed, Tfrac=0.1, anneal=False, chains=16):
    """Continuous single-residue Metropolis in TORSION space, budget-capped.

    The move class is matched to the samplers': one residue's (phi, psi) is redrawn from its
    own von Mises basin mixture.  `anneal=True` makes it simulated annealing on a geometric
    schedule from 4*T down to T/16.  `chains` independent replicas are advanced together --
    a batched thermostat WITH restarts built in (each chain is its own restart), which is a
    stronger classical control than one long chain and costs the identical budget: one
    objective read per proposal per chain.
    """
    rng = SD.stable_rng(tgt["pdb"], seed, "metro", Tfrac, bool(anneal), salt=SALT)
    s, _p0, _q0 = _scale(tgt, obj, rng)
    T0 = max(Tfrac * s, 1e-9)
    n = tgt["n"]
    b = (rng.random((chains, n)) < tgt["wmarg"][:, 1][None, :]).astype(np.int64)
    phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
    cur = obj(phi, psi)
    if cur.size < chains:
        return {"T0": T0, "scale": s}
    steps = 0
    ar = np.arange(chains)
    while obj.left >= chains:
        steps += 1
        frac = min(1.0, steps * chains / max(1, obj.budget))
        T = T0 * 4.0 * (1 / 64.0) ** frac if anneal else T0
        i = rng.integers(0, n, chains)
        c = (rng.random(chains) < tgt["wmarg"][i, 1]).astype(np.int64)
        np_ = phi.copy(); ns = psi.copy()
        np_[ar, i] = L.wrap(rng.vonmises(tgt["mu"][i, c, 0], tgt["kap"][i, c, 0]))
        ns[ar, i] = L.wrap(rng.vonmises(tgt["mu"][i, c, 1], tgt["kap"][i, c, 1]))
        ev = obj(np_, ns)
        if ev.size < chains:
            break
        de = ev - cur
        acc = (de <= 0) | (rng.random(chains)
                           < np.exp(np.clip(-de / max(T, 1e-12), -700.0, 0.0)))
        phi[acc] = np_[acc]; psi[acc] = ns[acc]; cur[acc] = ev[acc]
    return {"T0": T0, "scale": s, "chains": chains}


def run_lbfgs(tgt, obj, seed, maxiter=60):
    """Multi-start L-BFGS on the DEPLOYED objective, budget-counted.

    Identical arithmetic to `s15/align_lib.fit` (kind="squared", wpair=None): same
    `pj.frames` / `pj.build_ca_exact` / `pj._torsion_grad` chain.  Verified against
    `A.fit` in `qb_verify`.  This is the incumbent optimiser, put at matched budget.
    """
    from scipy.optimize import minimize
    from core import project as pj
    rng = SD.stable_rng(tgt["pdb"], seed, "lbfgs", salt=SALT)
    n = tgt["n"]
    i, j = tgt["i"], tgt["j"]
    inv = 1.0 / np.asarray(tgt["sd"], float)
    dhat = tgt["dhat"]

    def fg(x):
        phi, psi = x[:n], x[n:]
        obj.used += 1                                   # one objective read
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[i] - CA[j]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = (d - dhat) * inv
        f = float((r * r).sum())
        coef = ((2.0 * r * inv) / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        if obj.keep:
            obj._phi.append(phi[None].copy()); obj._psi.append(psi[None].copy())
            obj._e.append(np.array([f]))
        return f, pj._torsion_grad(G, CA, gCA)

    while obj.left > 0:
        b = (rng.random(n) < tgt["wmarg"][:, 1]).astype(np.int64)
        p0, q0 = L.draw_from_basins(b[None], tgt["mu"], tgt["kap"], rng)
        x0 = np.concatenate([p0[0], q0[0]])
        minimize(fg, x0, jac=True, method="L-BFGS-B",
                 options={"maxiter": maxiter, "maxfun": min(obj.left, 400),
                          "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return {}


# =============================================================== zero-information nulls
def run_marg(tgt, obj, seed):
    """PLAUSIBLE ZERO-INFORMATION NULL: matched empirical torsion marginals.

    Basins drawn i.i.d. at the pool marginal rate, angles from the same von Mises
    components.  Reads the objective ZERO times for its decisions -- the evaluations exist
    only so the emitted set is the same size as every other arm's.
    """
    rng = SD.stable_rng(tgt["pdb"], seed, "marg", salt=SALT)
    n = tgt["n"]
    while obj.left > 0:
        m = min(L.SHOTS, obj.left)
        b = (rng.random((m, n)) < tgt["wmarg"][:, 1][None, :]).astype(np.int64)
        phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
        obj(phi, psi)
    return {}


HELIX_PHI, HELIX_PSI = np.deg2rad(-57.0), np.deg2rad(-47.0)


def run_helix(tgt, obj, seed):
    """THE ZERO-INFORMATION CONTROL THE BRIEF NAMES: a constant ideal alpha-helix plus a
    matched isotropic torsion perturbation (sigma = the pool's mean per-residue circular sd).
    """
    rng = SD.stable_rng(tgt["pdb"], seed, "helix", salt=SALT)
    n = tgt["n"]
    R = np.hypot(np.sin(tgt["PHI"]).mean(0), np.cos(tgt["PHI"]).mean(0))
    R2 = np.hypot(np.sin(tgt["PSI"]).mean(0), np.cos(tgt["PSI"]).mean(0))
    sig = float(np.mean(np.sqrt(-2.0 * np.log(np.clip(np.r_[R, R2], 1e-6, 0.999999)))))
    while obj.left > 0:
        m = min(L.SHOTS, obj.left)
        phi = L.wrap(HELIX_PHI + rng.normal(0, sig, (m, n)))
        psi = L.wrap(HELIX_PSI + rng.normal(0, sig, (m, n)))
        obj(phi, psi)
    return {"sigma_rad": sig}


def pool_set(tgt, k=500):
    """THE INCUMBENT GENERATION MECHANISM: the shipped BLOSUM pool's own torsions, 0 evals.

    Put through the identical ideal-geometry builder as every sampler, so the comparison is
    on one manifold.
    """
    return tgt["PHI"][:k], tgt["PSI"][:k]


# =============================================================== the arm table
def arms():
    """name -> (kind, kwargs).  Every budgeted arm gets exactly `L.BUDGET` evaluations."""
    A = {}
    for a in ALPHAS:
        A[f"q_a{a:.2f}"] = ("circuit", dict(alpha=a, ansatz="mps2f", train=True))
    A["q_anneal"] = ("circuit", dict(alpha=0.25, ansatz="mps2f", train=True,
                                     anneal=(1.0, 0.05)))
    A["q_untrained"] = ("circuit", dict(alpha=0.25, ansatz="mps2f", train=False))
    for a in (0.05, 0.25):
        A[f"c_prod{a:.2f}"] = ("circuit", dict(alpha=a, ansatz="mps2fn", train=True))
    A["c_chain0.05"] = ("chain", dict(alpha=0.05))
    A["c_chain0.25"] = ("chain", dict(alpha=0.25))
    A["c_cem0.05"] = ("cem", dict(alpha=0.05))
    A["c_metroL"] = ("metro", dict(Tfrac=0.05))
    A["c_metroH"] = ("metro", dict(Tfrac=0.5))
    A["c_anneal"] = ("metro", dict(Tfrac=0.2, anneal=True))
    A["c_lbfgs"] = ("lbfgs", dict())
    A["c_marg"] = ("marg", dict())
    A["c_helix"] = ("helix", dict())
    A["pool500"] = ("pool", dict())
    return A


QUANTUM = ("q_a0.05", "q_a0.25", "q_a1.00", "q_anneal")
CLASSICAL = ("c_prod0.05", "c_prod0.25", "c_chain0.05", "c_chain0.25", "c_cem0.05",
             "c_metroL", "c_metroH", "c_anneal", "c_lbfgs")
NULLS = ("c_marg", "c_helix")


def run_arm(tgt, name, kind, kw, seed=0, budget=L.BUDGET):
    obj = L.new_obj(tgt, budget=budget)
    if kind == "pool":
        phi, psi = pool_set(tgt)
        e, _ = obj.raw(phi, psi)
        extra = {"evals": 0}
        r = L.readout(phi, psi, e, tgt, tag=name)
        r.update(extra)
        return r
    if kind == "circuit":
        extra = run_circuit(tgt, obj, seed=seed, **kw)
    elif kind == "chain":
        extra = run_chain(tgt, obj, seed=seed, **kw)
    elif kind == "cem":
        extra = run_cem(tgt, obj, seed=seed, **kw)
    elif kind == "metro":
        extra = run_metro(tgt, obj, seed=seed, **kw)
    elif kind == "lbfgs":
        extra = run_lbfgs(tgt, obj, seed=seed, **kw)
    elif kind == "marg":
        extra = run_marg(tgt, obj, seed=seed, **kw)
    elif kind == "helix":
        extra = run_helix(tgt, obj, seed=seed, **kw)
    else:
        raise ValueError(kind)
    phi, psi, e = obj.seen()
    r = L.readout(phi, psi, e, tgt, tag=name)
    r.update({k: v for k, v in extra.items() if k != "cvar_traj"})
    r["evals"] = int(obj.used)
    return r
