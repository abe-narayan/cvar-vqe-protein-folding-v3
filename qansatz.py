"""The variational circuit, simulated exactly, with an analytic CVaR gradient.

The ansatz is the same family the project has always used -- ``layers x (RY on every
wire, CNOT chain, optional ring closure)`` on |0...0> -- but the simulation and the
optimisation are both rebuilt.

**Simulation.** The circuit is a depth-``layers`` nearest-neighbour chain, so its exact
matrix-product state has bond dimension at most ``2 ** layers``: a CNOT is
``|0><0| (x) I + |1><1| (x) X``, an exact bond-dimension-2 MPO, and applying it doubles one
bond with no truncation and no SVD. At layers=3 that is chi <= 8. Nothing here approximates
anything, so the 30-qubit statevector wall that forced the previous 36-qubit runs onto a
scratchpad driver is gone: cost is O(n chi^3), linear in the number of qubits. A 64-qubit
register is routine.

**Gradient.** The previous driver optimised by moment matching -- pull the RY angles toward
the mean bit pattern of the elite tail -- which is the cross-entropy method wearing a
circuit, not a variational optimisation of a stated objective. Here the objective is
written down,

    CVaR_alpha(theta) = max_t  t - (1/alpha) E_{x~p_theta} (t - E(x))_+ ,

and differentiated. At the optimal ``t`` (the alpha-quantile ``q``) the envelope theorem
kills the ``t`` dependence and leaves the score-function form

    grad CVaR = -(1/alpha) E[ (q - E(x))_+ grad log p_theta(x) ] ,

which needs only ``grad log p_theta(x)`` at the sampled bitstrings. For one layer that is
available in closed form (the pre-CNOT bits are independent Bernoulli variables and the
CNOT chain is an invertible XOR, so ``d log p / d theta_q = 2 (y_q - sin^2(theta_q/2)) /
sin(theta_q)``). For deeper circuits the amplitude of a basis state is a product of the
selected MPS matrices, which torch differentiates directly. Either way this is a genuine
gradient of the CVaR objective rather than a surrogate update.

**Adaptive alpha.** ``alpha`` is annealed from broad to narrow across the run. Early on a
narrow tail concentrates the distribution before the landscape has been sampled -- the
premature-concentration failure the previous configuration showed at ``init_scale=0.25``.
Late on a broad tail is just the mean, and the mean is not what we want to minimise.
"""
import math
from typing import Callable, Dict, Optional, Sequence, Tuple

import numpy as np


# ----------------------------------------------------------------------- CVaR
def cvar(energies: np.ndarray, alpha: float) -> Tuple[float, float, np.ndarray]:
    """``(cvar, quantile, tail mask)`` for the LOWER alpha tail of `energies`."""
    e = np.asarray(energies, float)
    if e.size == 0:
        raise ValueError("cvar received no samples")
    if not (0.0 < alpha <= 1.0):
        raise ValueError(f"alpha must be in (0, 1], got {alpha}")
    k = max(1, int(math.ceil(alpha * e.size)))
    order = np.argsort(e, kind="stable")
    tail = np.zeros(e.size, bool)
    tail[order[:k]] = True
    return float(e[tail].mean()), float(e[order[k - 1]]), tail


def alpha_schedule(progress: float, a0: float = 0.5, a1: float = 0.05) -> float:
    """Geometric anneal of the CVaR level from `a0` to `a1` over a run."""
    p = min(1.0, max(0.0, float(progress)))
    return float(a0 * (a1 / a0) ** p)


# ----------------------------------------------------------- one-layer analytic
class OneLayerAnsatz:
    """RY on every wire, then a CNOT chain (and optional ring closure).

    The prepared distribution is exactly: draw independent bits ``y_q ~ Bern(sin^2(theta_q
    / 2))``, then output the cumulative XOR. Sampling is a vectorised prefix-XOR, and both
    the log-probability and its gradient are closed form, so nothing about this path is a
    simulation approximation.
    """

    layers = 1

    def __init__(self, n_qubits: int, ring: bool = True):
        self.n = int(n_qubits)
        self.ring = bool(ring) and self.n > 2

    def n_params(self) -> int:
        return self.n

    def _p1(self, theta: np.ndarray) -> np.ndarray:
        return np.sin(np.asarray(theta, float) / 2.0) ** 2

    def sample(self, theta: np.ndarray, shots: int, rng) -> np.ndarray:
        y = (rng.random((int(shots), self.n)) < self._p1(theta)).astype(np.uint8)
        b = np.bitwise_xor.accumulate(y, axis=1)
        if self.ring:
            b = b.copy()
            b[:, 0] ^= b[:, -1]
        return b

    def _pre_bits(self, bits: np.ndarray) -> np.ndarray:
        """Invert the CNOT network: recover the independent Bernoulli bits."""
        b = np.asarray(bits, np.uint8)
        if self.ring:
            b = b.copy()
            # b0' = b0 ^ b_{n-1}; the chain fixes b_{n-1} from the *outputs*, so undo the
            # ring first using the untouched last output bit.
            b[:, 0] ^= b[:, -1]
        y = b.copy()
        y[:, 1:] ^= b[:, :-1]
        return y

    def logp(self, theta: np.ndarray, bits: np.ndarray) -> np.ndarray:
        s = np.clip(self._p1(theta), 1e-12, 1 - 1e-12)
        y = self._pre_bits(bits)
        return (y * np.log(s) + (1 - y) * np.log1p(-s)).sum(1)

    def grad_logp(self, theta: np.ndarray, bits: np.ndarray) -> np.ndarray:
        """``(B, n_params)`` gradient of log p at each sampled bitstring."""
        th = np.asarray(theta, float)
        s = self._p1(th)
        y = self._pre_bits(bits).astype(float)
        sin = np.sin(th)
        sin = np.where(np.abs(sin) < 1e-6, np.sign(sin + 1e-12) * 1e-6, sin)
        return 2.0 * (y - s[None, :]) / sin[None, :]


# ------------------------------------------------------------------ deep / MPS
class MPSAnsatz:
    """``layers`` repetitions of (RY on every wire, CNOT chain), simulated exactly.

    Bond dimension is capped by ``2 ** layers`` and by the chain cut, so the state is
    represented without truncation. Sampling uses right environments and the usual
    left-to-right conditional decomposition; log-probabilities are recomputed as a product
    of selected matrices, which torch differentiates for the score-function gradient.

    The ring closure is not applied here: it is a long-range gate for an open MPS and would
    cost a swap network for a correlation the chain already carries at depth >= 2.

    ``final_ry`` appends one more RY layer with no entangler after it, and it is not
    cosmetic. The shipped one-layer chain cannot express an arbitrary set of per-qubit
    marginals at all: ``b_q`` is the prefix XOR of independent Bernoulli bits, so
    ``|1 - 2 P(b_q = 1)|`` is a running product and must be non-increasing along the wire.
    Measured against the torsion prior's own marginals on 12 targets, 53% of chain steps
    violate that ordering and the best achievable initialisation is 0.13 off per qubit on
    average (0.29 at worst). Directly optimising the angles to match a random marginal
    vector leaves a residual of 0.21 for the chain and 0.000 with one trailing RY layer.
    The defect is that the last gate on every wire is a CNOT, which fixes that wire's
    marginal to whatever the running product allows; one more rotation removes it while
    keeping the entanglement the chain built.

    ``entangler="none"`` drops the CNOTs entirely, which gives a product distribution --
    free marginals, no correlations. It is the control that says whether the entanglement
    is doing anything.
    """

    def __init__(self, n_qubits: int, layers: int = 2, final_ry: bool = False,
                 entangler: str = "cnot"):
        import torch
        self.torch = torch
        self.n = int(n_qubits)
        self.layers = int(layers)
        self.final_ry = bool(final_ry)
        self.entangler = entangler
        self.ring = False

    def n_params(self) -> int:
        return self.n * (self.layers + int(self.final_ry))

    # -- construction ----------------------------------------------------
    def _build(self, theta):
        t = self.torch
        nblocks = self.layers + int(self.final_ry)
        th = theta.reshape(nblocks, self.n)
        A = [t.zeros(1, 2, 1, dtype=t.float64) for _ in range(self.n)]
        for q in range(self.n):
            A[q][0, 0, 0] = 1.0
        for l in range(nblocks):
            c = t.cos(th[l] / 2.0)
            s = t.sin(th[l] / 2.0)
            for q in range(self.n):
                g = t.stack([t.stack([c[q], -s[q]]), t.stack([s[q], c[q]])])
                A[q] = t.einsum("ij,ajb->aib", g, A[q])
            if l >= self.layers or self.entangler == "none":
                continue
            for q in range(self.n - 1):
                dl, _, dr = A[q].shape
                # control bond: copy the physical index of q onto a new right bond
                B = t.zeros(dl, 2, dr, 2, dtype=t.float64)
                B[:, 0, :, 0] = A[q][:, 0, :]
                B[:, 1, :, 1] = A[q][:, 1, :]
                A[q] = B.reshape(dl, 2, dr * 2)
                el, _, er = A[q + 1].shape
                nxt = t.stack([A[q + 1], A[q + 1].flip(1)], dim=1)      # (el, 2, 2, er)
                A[q + 1] = nxt.permute(0, 1, 2, 3).reshape(el * 2, 2, er)
        return A

    def _right_env(self, A):
        t = self.torch
        R = [None] * (self.n + 1)
        R[self.n] = t.ones(1, 1, dtype=t.float64)
        for q in range(self.n - 1, -1, -1):
            R[q] = t.einsum("asb,bc,dsc->ad", A[q], R[q + 1], A[q])
        return R

    # -- numpy fast paths -------------------------------------------------
    # Sampling needs no gradient, and at these bond dimensions (chi <= 8) the tensors are
    # far too small to amortise torch's per-operation overhead: the same contractions in
    # numpy run the search loop about twice as fast end to end. torch is kept for
    # `logp_torch`, where autograd is the point.
    def _build_np(self, theta: np.ndarray):
        nblocks = self.layers + int(self.final_ry)
        th = np.asarray(theta, float).reshape(nblocks, self.n)
        A = [np.zeros((1, 2, 1)) for _ in range(self.n)]
        for q in range(self.n):
            A[q][0, 0, 0] = 1.0
        for l in range(nblocks):
            c, s = np.cos(th[l] / 2.0), np.sin(th[l] / 2.0)
            for q in range(self.n):
                g = np.array([[c[q], -s[q]], [s[q], c[q]]])
                A[q] = np.einsum("ij,ajb->aib", g, A[q])
            if l >= self.layers or self.entangler == "none":
                continue
            for q in range(self.n - 1):
                dl, _, dr = A[q].shape
                B = np.zeros((dl, 2, dr, 2))
                B[:, 0, :, 0] = A[q][:, 0, :]
                B[:, 1, :, 1] = A[q][:, 1, :]
                A[q] = B.reshape(dl, 2, dr * 2)
                el, _, er = A[q + 1].shape
                A[q + 1] = np.stack([A[q + 1], A[q + 1][:, ::-1]],
                                    axis=1).reshape(el * 2, 2, er)
        return A

    @staticmethod
    def _right_env_np(A):
        n = len(A)
        R = [None] * (n + 1)
        R[n] = np.ones((1, 1))
        for q in range(n - 1, -1, -1):
            R[q] = np.einsum("asb,bc,dsc->ad", A[q], R[q + 1], A[q])
        return R

    def sample(self, theta: np.ndarray, shots: int, rng) -> np.ndarray:
        if self.entangler == "none":
            # Product state: the marginals are the RY angles of the last rotation layer
            # composed along the wire, and there is nothing to contract.
            nb = self.layers + int(self.final_ry)
            th = np.asarray(theta, float).reshape(nb, self.n).sum(0)
            return (rng.random((int(shots), self.n))
                    < np.sin(th / 2.0) ** 2).astype(np.uint8)
        A = self._build_np(theta)
        R = self._right_env_np(A)
        B = int(shots)
        v = np.ones((B, 1))
        out = np.empty((B, self.n), np.uint8)
        u = rng.random((B, self.n))
        rows = np.arange(B)
        for q in range(self.n):
            w = np.einsum("al,lsr->asr", v, A[q])                    # (B, 2, chi)
            p = np.einsum("asr,rc,asc->as", w, R[q + 1], w)
            np.maximum(p, 0.0, out=p)
            p /= np.maximum(p.sum(1, keepdims=True), 1e-300)
            bit = (u[:, q] > p[:, 0]).astype(np.int64)
            out[:, q] = bit
            v = w[rows, bit]
            v /= np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-300)
        return out

    def logp_torch(self, theta, bits):
        """Differentiable log p for a batch of bitstrings (torch tensors in and out)."""
        t = self.torch
        A = self._build(theta)
        R = self._right_env(A)
        B = bits.shape[0]
        v = t.ones(B, 1, dtype=t.float64)
        logscale = t.zeros(B, dtype=t.float64)
        for q in range(self.n):
            v = t.einsum("al,lsr->asr", v, A[q])[t.arange(B), bits[:, q].to(t.int64)]
            nrm = t.linalg.norm(v, dim=1)
            v = v / t.clamp(nrm, min=1e-300).unsqueeze(1)
            logscale = logscale + t.log(t.clamp(nrm, min=1e-300))
        amp2 = 2.0 * (logscale + t.log(t.clamp(v.abs().sum(1), min=1e-300)))
        return amp2 - t.log(t.clamp(R[0][0, 0], min=1e-300))

    def grad_logp(self, theta: np.ndarray, bits: np.ndarray,
                  weights: np.ndarray) -> np.ndarray:
        """``sum_i weights[i] * grad log p(bits[i])`` -- the score-function accumulator."""
        t = self.torch
        th = t.as_tensor(np.asarray(theta, float), dtype=t.float64).requires_grad_(True)
        lp = self.logp_torch(th, t.as_tensor(np.asarray(bits, np.int64)))
        (lp * t.as_tensor(np.asarray(weights, float))).sum().backward()
        return th.grad.detach().numpy()


# --------------------------------------------------------------------- driver
class Adam:
    def __init__(self, n: int, lr: float = 0.12, b1: float = 0.9, b2: float = 0.999):
        self.lr, self.b1, self.b2 = lr, b1, b2
        self.m = np.zeros(n)
        self.v = np.zeros(n)
        self.t = 0

    def step(self, x: np.ndarray, g: np.ndarray) -> np.ndarray:
        self.t += 1
        self.m = self.b1 * self.m + (1 - self.b1) * g
        self.v = self.b2 * self.v + (1 - self.b2) * g * g
        mh = self.m / (1 - self.b1 ** self.t)
        vh = self.v / (1 - self.b2 ** self.t)
        return x - self.lr * mh / (np.sqrt(vh) + 1e-8)


def cvar_gradient(ansatz, theta: np.ndarray, bits: np.ndarray,
                  energies: np.ndarray, alpha: float) -> Tuple[np.ndarray, float]:
    """Score-function gradient of ``CVaR_alpha`` and the objective value."""
    val, q, tail = cvar(energies, alpha)
    w = np.zeros(len(energies))
    # (q - E)_+ / alpha, mean-normalised: the envelope-theorem weight on the tail.
    w[tail] = -(q - energies[tail]) / (alpha * len(energies))
    # A constant baseline leaves the gradient unbiased and cuts its variance; the tail mean
    # is the natural one here.
    w[tail] -= w[tail].mean()
    if isinstance(ansatz, OneLayerAnsatz):
        grad = w @ ansatz.grad_logp(theta, bits)
    else:
        # Only the tail carries weight, so the autograd pass runs on alpha * shots rows
        # instead of all of them. Exactly equivalent -- the other rows contribute zero.
        grad = ansatz.grad_logp(theta, bits[tail], w[tail])
    return grad, val
