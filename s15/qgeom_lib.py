"""SPRINT 15 / QGEOM -- shared machinery for the quantum-geometry workstream.

WHAT THIS ADDS OVER SPRINT 14.

`s14/vqe_ansatz.fubini_study` computes the Fubini-Study metric by CENTRAL DIFFERENCES on the
exact statevector.  Everything here uses an **exact analytic derivative** instead, which is a
mathematically independent route and is exact to machine precision rather than to O(h^2):

    RY(t) = cos(t/2) I - i sin(t/2) Y
    d/dt RY(t) = -(1/2) sin(t/2) I - (i/2) cos(t/2) Y
               = (1/2) [ cos((t+pi)/2) I - i sin((t+pi)/2) Y ]
               = (1/2) RY(t + pi)

so for a circuit in which each parameter appears in exactly one RY gate,

    d|psi>/dtheta_j  =  (1/2) |psi(theta with theta_j -> theta_j + pi)>            (EXACT)

No finite difference, no h.  `fs_metric` and `fs_metric_fd` are cross-checked against each
other in `qgeom_metric.verify()`.

THE SECOND NEW OBJECT is the CLASSICAL Fisher information of the measured distribution,

    F_ij = sum_x (1/p) (dp/dtheta_i)(dp/dtheta_j),      p(x) = |psi(x)|^2

which is the metric a *sampling* optimiser of a *diagonal* cost actually lives on.  For a
real-amplitude circuit (RY + CNOT has real amplitudes) `F = 4 g` EXACTLY -- see
`qgeom_metric` for the proof and the numerical confirmation -- so the QNG question and the
classical-natural-gradient question are the same question here, up to a learning rate.

NATIVE-FREE.  Nothing in this module reads a native coordinate, torsion or RMSD.  The
`Enum.rmsd` column is used only for post-hoc scoring by the callers, and the retrieval prior
(`s14.retprior`) is a BLOSUM retrieval output.
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

RESULTS = os.path.join(ROOT, "s15", "results")
CACHE = os.path.join(ROOT, "s15", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

from s14.vqe_ansatz import FlexCircuit, PATTERNS, fubini_study as fs_metric_fd  # noqa: E402
from s14 import vqe_lib as V                                                    # noqa: E402


# =============================================================== incremental checkpoint
_CK = {}


def ck(tag: str, key: str, value) -> None:
    """Write a partial result to `s15/results/qgeom_<tag>.json` IMMEDIATELY.

    Sprint 14 lost an entire ansatz study to a module that wrote its JSON only on
    completion.  Every cell in this workstream calls this the moment it has a number.

    MERGE-ON-WRITE.  The file on disk is re-read and merged before every write, so a
    process that never called `ck_load` cannot destroy keys another process wrote.  I lost
    a completed 27-cell Part B table to exactly that on 2026-09-05; the analysis re-runs in
    five minutes but the failure mode is silent and would not be in a rerun-free study.
    """
    path0 = os.path.join(RESULTS, f"qgeom_{tag}.json")
    if os.path.exists(path0):
        try:
            with open(path0) as fh:
                disk = json.load(fh)
            mem = _CK.setdefault(tag, {})
            for k, v in disk.items():
                if isinstance(v, dict) and isinstance(mem.get(k), dict):
                    for kk, vv in v.items():          # one-level deep merge; memory wins
                        mem[k].setdefault(kk, vv)
                else:
                    mem.setdefault(k, v)
        except Exception:
            pass
    d = _CK.setdefault(tag, {})
    d[key] = value
    d["_written"] = time.strftime("%Y-%m-%d %H:%M:%S")
    path = os.path.join(RESULTS, f"qgeom_{tag}.json")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    # Windows os.replace fails with WinError 5 if any reader holds the destination open.
    # A monitoring `json.load` in another shell killed a run this way; retry, then fall
    # back to a plain write rather than losing the cell.
    for _try in range(30):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.4)
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    try:
        os.remove(tmp)
    except OSError:
        pass


def ck_load(tag: str) -> dict:
    path = os.path.join(RESULTS, f"qgeom_{tag}.json")
    if os.path.exists(path):
        with open(path) as fh:
            d = json.load(fh)
        _CK[tag] = d
        return d
    return {}


def wait_mem(min_gb=1.0, tag=""):
    return V.wait_for_memory(min_gb, tag)


# ================================================================= exact derivatives
def deriv_states(circ: FlexCircuit, theta) -> np.ndarray:
    """(P, dim) exact d|psi>/dtheta_j via the RY shift-by-pi identity.  No finite difference."""
    th = np.asarray(theta, float)
    P = th.size
    TH = np.repeat(th[None, :], P, axis=0)
    TH[np.arange(P), np.arange(P)] += np.pi
    return 0.5 * circ.states_batch(TH)


def deriv_states_lowmem(circ: FlexCircuit, theta) -> np.ndarray:
    """Same, one parameter at a time -- for wide registers where the batch does not fit."""
    th = np.asarray(theta, float)
    P = th.size
    D = np.empty((P, circ.dim))
    for j in range(P):
        t = th.copy()
        t[j] += np.pi
        D[j] = 0.5 * circ.state(t)
    return D


def fs_metric(circ: FlexCircuit, theta, lowmem=False) -> np.ndarray:
    """Fubini-Study metric ``g_ij = Re<d_i psi|d_j psi> - <d_i psi|psi><psi|d_j psi>``.

    Exact analytic derivatives.  Amplitudes are real so the Berry (imaginary) part is
    identically zero and `g` is the full quantum geometric tensor.
    """
    psi = circ.state(theta)
    D = deriv_states_lowmem(circ, theta) if lowmem else deriv_states(circ, theta)
    ov = D @ psi
    return D @ D.T - np.outer(ov, ov)


def classical_fim(circ: FlexCircuit, theta, lowmem=False) -> np.ndarray:
    """Classical Fisher information of the MEASURED distribution p = psi^2.

    ``F_ij = sum_x (dp_i)(dp_j)/p``.  States with p = 0 contribute nothing (both
    numerator factors carry a psi), and are masked rather than divided.
    """
    psi = circ.state(theta)
    D = deriv_states_lowmem(circ, theta) if lowmem else deriv_states(circ, theta)
    p = psi ** 2
    dp = 2.0 * psi[None, :] * D
    m = p > 1e-300
    W = dp[:, m] / np.sqrt(p[m])[None, :]
    return W @ W.T


# ============================================================== spectrum statistics
def spec_stats(g: np.ndarray, tol=1e-9) -> dict:
    """Everything the brief asks for from one symmetric metric."""
    P = g.shape[0]
    ev = np.linalg.eigvalsh((g + g.T) / 2.0)
    ev = np.clip(ev, 0.0, None)
    pos = ev[ev > tol * max(ev.max(), 1e-30)]
    rank = int(pos.size)
    cond = float(pos.max() / pos.min()) if rank else float("nan")
    cond_full = float(ev.max() / ev.min()) if ev.min() > 1e-15 else float("inf")
    s = pos / pos.sum() if rank else np.array([1.0])
    eff_rank = float(np.exp(-(s * np.log(s)).sum()))          # spectral (entropy) rank
    logvol = float(np.log(pos).sum()) if rank else float("-inf")   # log pseudo-determinant
    off = g - np.diag(np.diag(g))
    return {
        "n_params": int(P),
        "mean_gii": float(np.mean(np.diag(g))),
        "min_gii": float(np.min(np.diag(g))),
        "max_gii": float(np.max(np.diag(g))),
        "max_offdiag": float(np.abs(off).max()) if P > 1 else 0.0,
        "rank": rank,
        "cond": cond,
        "cond_full": cond_full,
        "eff_rank": eff_rank,
        "eff_rank_frac": float(eff_rank / P),
        "log_pseudo_det": logvol,
        "eig_max": float(ev.max()),
        "eig_min_pos": float(pos.min()) if rank else float("nan"),
        "eig_min": float(ev.min()),
        "trace": float(np.trace(g)),
        "max_dev_from_I4": float(np.abs(g - np.eye(P) / 4.0).max()),
        "eigs": [float(x) for x in ev],
    }


# ==================================================================== reference thetas
def uniform_theta(circ: FlexCircuit) -> np.ndarray:
    """theta giving the exactly uniform distribution: RY(pi/2) on the first layer, 0 after.

    CNOTs are permutations of basis states and preserve a uniform |+...+>, and RY(0) = I,
    so p(x) = 2^-n exactly for every x.  Verified in `qgeom_metric.verify()`.
    """
    th = np.zeros(circ.n_params())
    th[: circ.n] = np.pi / 2
    return th


def random_theta(circ: FlexCircuit, seed=0, scale=0.9) -> np.ndarray:
    return np.random.default_rng(seed).normal(0.0, scale, circ.n_params())


def haar_theta(circ: FlexCircuit, seed=0) -> np.ndarray:
    return np.random.default_rng(seed).uniform(0.0, 2 * np.pi, circ.n_params())


# ================================================= target-conditioned distributions
def product_target(P: np.ndarray) -> np.ndarray:
    """(n, k) per-residue distribution -> full (k**n,) product distribution.

    Ordering is the project's big-endian base-k convention (residue 0 most significant),
    which is bit-identical to the MSB-first binary register at power-of-two k.
    """
    P = np.asarray(P, float)
    P = P / P.sum(1, keepdims=True)
    out = np.array([1.0])
    for i in range(P.shape[0]):
        out = np.outer(out, P[i]).ravel()
    return out


def prior_entropy_bits(P: np.ndarray) -> float:
    P = np.asarray(P, float)
    P = P / P.sum(1, keepdims=True)
    return float(-(P * np.log2(np.maximum(P, 1e-300))).sum())


def scramble_prior(P: np.ndarray, seed=0, mode="state") -> np.ndarray:
    """Entropy-MATCHED null: the same per-residue distributions, wrong assignment.

    mode="state"   permute the k state labels within each residue independently.
                   Per-residue entropies are IDENTICAL, the target information is destroyed.
    mode="residue" permute which residue gets which distribution.
                   The multiset of per-residue distributions is identical.
    """
    r = np.random.default_rng(seed)
    P = np.asarray(P, float).copy()
    if mode == "state":
        for i in range(P.shape[0]):
            P[i] = P[i][r.permutation(P.shape[1])]
        return P
    if mode == "residue":
        return P[r.permutation(P.shape[0])]
    raise ValueError(mode)


# ============================================================== KL fitting (conditioning)
def fit_kl(circ: FlexCircuit, target: np.ndarray, theta0=None, iters=400, lr=0.10,
           seed=0, lowmem=False, track=None, verbose=False):
    """Fit the ansatz to a target distribution by minimising KL(target || p_theta).

    Exact gradient:  d/dtheta_j KL = - sum_x target(x) * dp_j(x) / p(x),
    with dp_j = 2 psi D_j from the analytic derivative.  Adam.

    `track` -- optional list of iteration numbers at which to record the metric spectrum,
    which is how "how does the geometry evolve during optimisation" is measured.
    """
    from core.quantum import Adam
    t = np.asarray(target, float)
    t = t / t.sum()
    th = (random_theta(circ, seed) if theta0 is None else np.asarray(theta0, float).copy())
    opt = Adam(circ.n_params(), lr=lr)
    hist, snaps = [], {}
    trackset = set(track or [])
    for it in range(iters):
        psi = circ.state(th)
        p = psi ** 2
        p = p / p.sum()
        kl = float((t[t > 0] * np.log(t[t > 0] / np.maximum(p[t > 0], 1e-300))).sum())
        hist.append(kl)
        if it in trackset:
            snaps[it] = spec_stats(fs_metric(circ, th, lowmem))
            snaps[it]["kl"] = kl
        D = deriv_states_lowmem(circ, th) if lowmem else deriv_states(circ, th)
        dp = 2.0 * psi[None, :] * D
        g = -(dp @ (t / np.maximum(p, 1e-300)))
        th = opt.step(th, g)
        if verbose and it % 50 == 0:
            print(f"    it {it:4d}  KL {kl:.5f}", flush=True)
    psi = circ.state(th)
    p = psi ** 2
    p = p / p.sum()
    kl = float((t[t > 0] * np.log(t[t > 0] / np.maximum(p[t > 0], 1e-300))).sum())
    hist.append(kl)
    return {"theta": th, "kl": kl, "kl_history": hist, "snaps": snaps,
            "tv": float(0.5 * np.abs(p - t).sum())}


# ================================================================ natural gradient
def nat_solve(g: np.ndarray, grad: np.ndarray, lam=1e-3, mode="tikhonov") -> np.ndarray:
    """Solve for the (regularised) natural-gradient direction.

    mode="tikhonov"  (g + lam*tr(g)/P * I)^-1 grad          -- scale-aware ridge
    mode="pinv"      pseudo-inverse with rcond=lam
    mode="diag"      diagonal preconditioner only (the cheap approximation QNG papers use)
    """
    P = g.shape[0]
    if mode == "diag":
        d = np.maximum(np.diag(g), lam)
        return grad / d
    if mode == "pinv":
        return np.linalg.pinv(g, rcond=lam) @ grad
    sc = max(np.trace(g) / P, 1e-30)
    return np.linalg.solve(g + lam * sc * np.eye(P), grad)


# ================================================================ misc statistics
def boot_ci(x, n=4000, seed=0, lo=2.5, hi=97.5):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return (float("nan"), float("nan"))
    r = np.random.default_rng(seed)
    m = r.integers(0, x.size, (n, x.size))
    b = x[m].mean(1)
    return (float(np.percentile(b, lo)), float(np.percentile(b, hi)))


def paired(a, b, seed=0):
    """Paired comparison a - b with everything the brief demands in one dict."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    d = a[ok] - b[ok]
    if d.size == 0:
        return {"n": 0}
    lo, hi = boot_ci(d, seed=seed)
    w = int((d < 0).sum())
    l = int((d > 0).sum())
    sd = float(d.std(ddof=1)) if d.size > 1 else 0.0
    top = np.sort(np.abs(d))[::-1]
    share = float(top[: max(1, d.size // 10)].sum() / max(np.abs(d).sum(), 1e-30))
    return {"n": int(d.size), "mean": float(d.mean()), "median": float(np.median(d)),
            "sd": sd, "mean_over_sd": float(d.mean() / sd) if sd > 0 else float("inf"),
            "ci_lo": lo, "ci_hi": hi, "win": w, "loss": l,
            "top10pct_abs_share": share,
            "verdict": ("SIGNIFICANT" if (lo > 0) == (hi > 0) else "NULL")}


def spearman(a, b):
    return V.spearman(a, b)


# ============================================================ adjoint (reverse-mode) gradient
def _ry_apply(v, n, q, t):
    """Apply RY(t) on wire q (MSB-first) to a real statevector, in place on a copy."""
    c, s = np.cos(t / 2.0), np.sin(t / 2.0)
    w = v.reshape(1 << q, 2, 1 << (n - q - 1))
    a = w[:, 0, :].copy()
    b = w[:, 1, :]
    out = np.empty_like(w)
    out[:, 0, :] = c * a - s * b
    out[:, 1, :] = s * a + c * b
    return out.reshape(-1)


def grad_adjoint(circ, theta, u):
    """Exact gradient of ``C = u . psi`` in O(gates) full-state operations, not O(P*gates).

    The forward pass builds `psi`; the reverse pass undoes one gate at a time, reads the
    derivative through the SAME `d/dt RY(t) = (1/2) RY(t+pi)` identity used by
    `deriv_states`, and back-propagates the adjoint.  Both routes are exact; this one is
    ~P times faster and is what makes an 18-qubit CVaR trajectory affordable.

    For a cost ``C = sum_x w(x) p(x)`` with ``p = psi^2``, pass ``u = 2 * w * psi``.
    Verified against `deriv_states` to ~1e-14 in `qgeom_lib.verify_adjoint`.
    """
    n = circ.n
    nb = circ.layers + int(circ.final_ry)
    TH = np.asarray(theta, float).reshape(nb, n)
    psi = circ.state(theta).copy()
    lam = np.asarray(u, float).copy()
    if circ.pairs:
        inv = np.empty_like(circ._perm)
        inv[circ._perm] = np.arange(circ.dim)
    grad = np.zeros((nb, n))
    for L in range(nb - 1, -1, -1):
        if L < circ.layers and circ.pairs:
            psi = psi[inv]                       # undo the entangler
            lam = lam[inv]                       # permutation is orthogonal: P^T = P^-1
        for q in range(n - 1, -1, -1):
            t = TH[L, q]
            psi = _ry_apply(psi, n, q, -t)       # psi is now the gate's INPUT
            grad[L, q] = lam @ _ry_apply(psi, n, q, t + np.pi) * 0.5
            lam = _ry_apply(lam, n, q, -t)       # RY^T = RY(-t)
    return grad.ravel()


def verify_adjoint(n=10, patterns=("ring", "block", "chain", "all_to_all", "none"),
                   Ls=(1, 2, 3), seeds=3):
    """Independent check: adjoint gradient vs the O(P) derivative-state gradient."""
    worst = 0.0
    for pat in patterns:
        for L in Ls:
            c = FlexCircuit(n, L, pat, 2)
            for s in range(seeds):
                th = random_theta(c, 17 * s + L)
                E = np.random.default_rng(s).standard_normal(c.dim)
                psi = c.state(th)
                D = deriv_states(c, th)
                g1 = 2.0 * (D @ (psi * E))
                g2 = grad_adjoint(c, th, 2.0 * psi * E)
                worst = max(worst, float(np.abs(g1 - g2).max()
                                         / max(np.abs(g1).max(), 1e-30)))
    return worst
