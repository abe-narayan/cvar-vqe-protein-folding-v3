"""Analytic-gradient continuous refinement of (phi, psi) -- the replacement for
`foldvqe.refine`.

Why
---
`foldvqe.refine` relaxes the selected structure off the discrete torsion grid with a
(1 + 24) evolution strategy: 40 generations of isotropic Gaussian proposals with a
geometrically shrinking step, twice. That is a derivative-free method applied to an
objective that is in fact smooth and cheap to differentiate exactly. Both terms
`FoldObjective.score_angles` evaluates -- the distogram's weighted L1 Bayes risk and the
CA/CB excluded-volume penalty -- are piecewise-linear in the interatomic distances, and
the distances are smooth in the torsions through the NeRF chain. The gradient in torsion
space therefore exists almost everywhere, and this module computes it exactly.

Two implementations, and the second is the one that ships
``TorchObjective``   the NeRF chain and the score rewritten in torch, so autograd supplies
                     the gradient. Correct, and it is the ORACLE the fast path is verified
                     against -- but it is slow here for a reason that has nothing to do
                     with the mathematics: an 11-residue chain is ~110 tensor ops on
                     (B, 3) arrays, and this box dispatches them at ~110 us each, so one
                     Adam step costs ~40 ms almost independently of the batch size.

``FastObjective``    the same energy and the same gradient in pure numpy, with the torsion
                     derivative in closed form instead of by autograd. Rotating torsion
                     ``k`` by ``dtheta`` rotates every atom downstream of its bond RIGIDLY
                     about that bond, so

                         dE/dtheta_k = e_k . sum_{a downstream} (x_a - o_k) x dE/dx_a

                     and the sums over "everything downstream" are one reverse cumulative
                     sum over residues. That turns the whole gradient into ~30 numpy
                     operations on ``(B, n, 3)`` arrays -- no per-atom Python, no autograd
                     tape -- and it is ~250x faster per structure than the torch path at
                     the batch sizes multi-start actually uses.

Both are checked, not asserted: `check_builder` (torch chain vs
`protein_geometry.build_backbone_batch`, ~2e-14), `check_gradient` (torch autograd vs
central differences), `check_analytic` (the numpy closed form vs torch autograd) and
`check_energy` (`FastObjective` vs `FoldObjective.score_angles`).

One deliberate numerical difference. `Distogram.score` does not evaluate the risk, it
GATHERS it from a 0.05 A lookup table, so as a function of distance it is piecewise
CONSTANT and its true gradient is zero everywhere. `FastObjective` reads the same table
but interpolates linearly between its nodes (and extrapolates past the ends, so a chain
stretched beyond 40 A still feels a restoring force). The value then differs from
`Distogram.score` by at most one cell of the table -- ~4e-4 in objective units against a
spread of ~1.5 -- and the derivative is the one the underlying continuous risk has.
Everything this module RETURNS as an energy is nevertheless recomputed with
`objective.score_angles`, so no reported number depends on that choice.

API
---
``refine(objective, states, seed=0, **kw) -> (phi, psi, energy, ca)``
    Drop-in for `foldvqe.refine`.
``refine_pool(objective, states_list, seed=0, **kw) -> (phi, psi, energy, ca)`` arrays
    Refine many reservoir members in ONE batched relaxation.
``refine_angles(objective, phi0, psi0, **kw)``
    The same, from angles rather than discrete states.

Keyword arguments (all of `refine`, `refine_pool`, `refine_angles`)
    ``method``    ``"adam"`` (default, numpy), ``"lbfgs"`` (scipy L-BFGS-B on the numpy
                  gradient), ``"torch-adam"``, ``"torch-lbfgs"`` (the reference path).
    ``n_starts``  perturbed restarts per structure, run as one batch. Row 0 is the
                  unperturbed incumbent, so the operator can never do worse than doing
                  nothing except through the objective's own noise.
    ``sigma0``    perturbation width for those restarts, DEGREES.
    ``steps``     optimiser iterations.
    ``lr``        Adam step, radians.
    ``trust``     cap on |angle - discrete start|, DEGREES, enforced by projection.
    ``hops``      basin-hopping rounds: kick ``hop_res`` residues by ``hop_sigma``
                  degrees, re-relax, keep the result only where it improved.
    ``fast_obj`` / ``torch_obj``  a prebuilt objective view, to amortise setup across
                  calls on the same target.
"""
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import protein_geometry as geo

try:
    import torch
    _HAVE_TORCH = True
    DTYPE = torch.float64
except Exception:                                    # pragma: no cover
    _HAVE_TORCH = False
    DTYPE = None

#: Threads the torch reference path may use. These problems are tiny -- a (64, 14)
#: parameter tensor and ~100 pair distances -- so intra-op threading costs more in
#: synchronisation than it saves, and other jobs share this box. Applied when a
#: `TorchObjective` is built rather than at import, so merely importing this module does
#: not reconfigure torch for whatever else is running in the same process.
TORCH_THREADS = 1


# ==================================================================== torch reference
def _place_atom_torch(a, b, c, length: float, angle: float, torsion):
    """Differentiable `protein_geometry._place_atom_batch`.

    The numpy version's degenerate-normal fallback (|n| < 1e-9 -> +z) becomes a clamped
    norm: the branch is unreachable for a real chain and a `where` on it would put a NaN
    into the backward pass.
    """
    bc = c - b
    bc = bc / bc.norm(dim=-1, keepdim=True).clamp_min(1e-9)
    ab = b - a
    nrm = torch.cross(ab, bc, dim=-1)
    nrm = nrm / nrm.norm(dim=-1, keepdim=True).clamp_min(1e-9)
    m = torch.cross(nrm, bc, dim=-1)
    d0 = -length * math.cos(angle)
    sa = length * math.sin(angle)
    return (c + d0 * bc
            + (sa * torch.cos(torsion)).unsqueeze(-1) * m
            + (sa * torch.sin(torsion)).unsqueeze(-1) * nrm)


def build_backbone_torch(phi, psi, omega: float = geo.OMEGA_TRANS,
                         with_o: bool = False) -> Dict[str, "torch.Tensor"]:
    """`geo.build_backbone_batch` in torch. ``phi``/``psi`` are ``(B, n)`` radians."""
    B, n = phi.shape
    dev, dt = phi.device, phi.dtype
    zero = torch.zeros(B, 3, dtype=dt, device=dev)
    ca0 = torch.tensor([geo.BOND_N_CA, 0.0, 0.0], dtype=dt, device=dev).expand(B, 3)
    c0 = torch.tensor([geo.BOND_N_CA - geo.BOND_CA_C * math.cos(geo.ANGLE_N_CA_C),
                       geo.BOND_CA_C * math.sin(geo.ANGLE_N_CA_C), 0.0],
                      dtype=dt, device=dev).expand(B, 3)
    om = torch.full((B,), float(omega), dtype=dt, device=dev)
    N, CA, C = [zero], [ca0], [c0]
    for i in range(n - 1):
        N.append(_place_atom_torch(N[i], CA[i], C[i],
                                   geo.BOND_C_N, geo.ANGLE_CA_C_N, psi[:, i]))
        CA.append(_place_atom_torch(CA[i], C[i], N[i + 1],
                                    geo.BOND_N_CA, geo.ANGLE_C_N_CA, om))
        C.append(_place_atom_torch(C[i], N[i + 1], CA[i + 1],
                                   geo.BOND_CA_C, geo.ANGLE_N_CA_C, phi[:, i + 1]))
    Nt, CAt, Ct = torch.stack(N, 1), torch.stack(CA, 1), torch.stack(C, 1)
    b, d = CAt - Nt, Ct - CAt
    CBt = -0.58273431 * torch.cross(b, d, dim=-1) + 0.56802827 * b - 0.54067466 * d + CAt
    out = {"N": Nt, "CA": CAt, "C": Ct, "CB": CBt}
    if with_o:
        out["O"] = torch.stack(
            [_place_atom_torch(Nt[:, i], CAt[:, i], Ct[:, i], geo.BOND_C_O,
                               geo.ANGLE_CA_C_O, psi[:, i] + math.pi) for i in range(n)], 1)
    return out


class _MemberRisk:
    """One `Distogram`'s weighted L1 Bayes risk as a differentiable function of CA."""

    def __init__(self, dist, device, dtype, exact: bool = False):
        import distogram as dgm
        self.i = torch.as_tensor(np.asarray(dist.i, np.int64), device=device)
        self.j = torch.as_tensor(np.asarray(dist.j, np.int64), device=device)
        coef = np.asarray(dist.prob, float) * np.asarray(dist.w, float)[:, None]
        self.coef = torch.as_tensor(coef, dtype=dtype, device=device)
        self.centres = torch.as_tensor(np.asarray(dgm.CENTRES, float),
                                       dtype=dtype, device=device)
        self.exact = bool(exact)
        self.g0 = float(dist.grid[0])
        self.dg = float(dist.grid[1] - dist.grid[0])
        self.gmax = float(len(dist.grid) - 1)

    def __call__(self, ca):
        d = (ca[:, self.i, :] - ca[:, self.j, :]).norm(dim=-1)
        if self.exact:
            # Snap to `Distogram`'s lookup grid in the forward pass, let the CONTINUOUS
            # risk's gradient through unchanged (straight-through), so the value matches
            # `score_angles` while the gradient stays non-zero.
            q = torch.trunc((d - self.g0) / self.dg).clamp(0.0, self.gmax)
            d = d + ((self.g0 + self.dg * q) - d).detach()
        return (self.coef * (d.unsqueeze(-1) - self.centres).abs()).sum(-1).mean(-1)


class TorchObjective:
    """Autograd view of `objective.FoldObjective.score_angles`, used as the oracle.

    Reproduces exactly the terms that function evaluates -- ``w_dist * dist.score(CA) +
    w_clash * clash(CA, CB)``, plus the optional consensus restraint -- and nothing else.
    The MRF term is a function of the DISCRETE states and is absent from `score_angles`,
    so it is absent here; ``trust`` is the stand-in for it.
    """

    def __init__(self, objective, device: str = "cpu", dtype=None, exact: bool = False):
        dtype = DTYPE if dtype is None else dtype
        if TORCH_THREADS:
            torch.set_num_threads(int(TORCH_THREADS))
        self.obj = objective
        self.n = objective.n
        self.device = torch.device(device)
        self.dtype = dtype
        self.exact = bool(exact)
        self.w_dist = float(objective.w_dist)
        self.w_clash = float(objective.w_clash)
        d = objective.dist
        self.members: List[_MemberRisk] = []
        self.mu = self.sigma = None
        self.disagree = 0.0
        if d is not None and self.w_dist:
            mem = getattr(d, "members", None)
            if mem:
                self.members = [_MemberRisk(m, self.device, dtype, exact) for m in mem]
                self.mu = torch.as_tensor(np.asarray(d.mu, float),
                                          dtype=dtype, device=self.device)
                self.sigma = torch.as_tensor(np.asarray(d.sigma, float),
                                             dtype=dtype, device=self.device)
                self.disagree = float(getattr(d, "disagree", 0.0))
            else:
                self.members = [_MemberRisk(d, self.device, dtype, exact)]
        iu = np.triu_indices(self.n, k=3)
        self.ci = torch.as_tensor(iu[0].astype(np.int64), device=self.device)
        self.cj = torch.as_tensor(iu[1].astype(np.int64), device=self.device)
        self.cut_ca = float(objective.clash_ca)
        self.cut_cb = float(objective.clash_cb)
        self.w_cons = float(getattr(objective, "w_cons", 0.0) or 0.0)
        cd = getattr(objective, "_cons_d", None)
        cw = getattr(objective, "_cons_w", None)
        self.cons_d = None if cd is None else torch.as_tensor(
            np.asarray(cd, float), dtype=dtype, device=self.device)
        self.cons_w = None if cw is None else torch.as_tensor(
            np.asarray(cw, float), dtype=dtype, device=self.device)
        self.cons_i = self.cons_j = None
        if self.cons_d is not None and d is not None:
            self.cons_i = torch.as_tensor(np.asarray(d.i, np.int64), device=self.device)
            self.cons_j = torch.as_tensor(np.asarray(d.j, np.int64), device=self.device)

    def set_exact(self, flag: bool) -> None:
        self.exact = bool(flag)
        for m in self.members:
            m.exact = bool(flag)

    def dist_score(self, ca):
        if not self.members:
            return torch.zeros(len(ca), dtype=self.dtype, device=self.device)
        if self.mu is None:
            return self.members[0](ca)
        v = torch.stack([m(ca) for m in self.members])
        v = (v - self.mu[:, None]) / self.sigma[:, None]
        out = v.mean(0)
        if self.disagree and len(v) > 1:
            out = out + self.disagree * (v.max(0).values - v.min(0).values)
        return out

    def clash(self, coords):
        out = torch.zeros(len(coords["CA"]), dtype=self.dtype, device=self.device)
        for key, cut in (("CA", self.cut_ca), ("CB", self.cut_cb)):
            d = (coords[key][:, self.ci, :] - coords[key][:, self.cj, :]).norm(dim=-1)
            out = out + (cut - d).clamp_min(0.0).sum(-1)
        return out

    def energy(self, phi, psi):
        coords = build_backbone_torch(phi, psi)
        e = self.w_dist * self.dist_score(coords["CA"])
        if self.w_clash:
            e = e + self.w_clash * self.clash(coords)
        if self.w_cons and self.cons_d is not None:
            ca = coords["CA"]
            dd = (ca[:, self.cons_i, :] - ca[:, self.cons_j, :]).norm(dim=-1)
            err = (dd - self.cons_d).abs()
            if self.cons_w is not None:
                err = err * self.cons_w
            e = e + self.w_cons * err.mean(-1)
        return e

    def _t(self, a):
        return torch.as_tensor(np.asarray(a, float), dtype=self.dtype, device=self.device)

    def energy_np(self, phi, psi) -> np.ndarray:
        with torch.no_grad():
            return self.energy(self._t(phi), self._t(psi)).cpu().numpy()


# ==================================================================== fast numpy path
def _incidence(n: int, i: np.ndarray, j: np.ndarray) -> np.ndarray:
    """``M[a, p] = [a == i_p] - [a == j_p]``, so a per-pair vector quantity scatters onto
    atoms with one ``einsum`` instead of an ``np.add.at``."""
    M = np.zeros((n, len(i)))
    M[i, np.arange(len(i))] += 1.0
    M[j, np.arange(len(j))] -= 1.0
    return M


class _FastMember:
    """A `Distogram` member's risk and its d(risk)/d(distance), by table interpolation."""

    def __init__(self, dist, n: int):
        self.i = np.asarray(dist.i, np.int64)
        self.j = np.asarray(dist.j, np.int64)
        self.table = np.asarray(dist._risk, float)          # (npairs, ngrid), w applied
        self.g0 = float(dist.grid[0])
        self.dg = float(dist.grid[1] - dist.grid[0])
        self.ngrid = self.table.shape[1]
        self.slope = np.diff(self.table, axis=1) / self.dg   # (npairs, ngrid - 1)
        self.npairs = len(self.i)
        self.M = _incidence(n, self.i, self.j)
        self.p = np.arange(self.npairs)

    def risk_and_dd(self, d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """``(risk (B,), d risk / d d (B, npairs))`` for pair distances ``d`` ``(B, npairs)``.

        Linear interpolation on `Distogram`'s own table, extrapolated past both ends by the
        end slopes so a collapsed or a hyper-extended chain still has a restoring gradient.
        """
        t = (d - self.g0) / self.dg
        g = np.clip(np.floor(t), 0, self.ngrid - 2).astype(np.int64)
        frac = t - g
        s = self.slope[self.p[None, :], g]
        val = self.table[self.p[None, :], g] + frac * (s * self.dg)
        return val.mean(1), s / self.npairs


class FastObjective:
    """`FoldObjective.score_angles` and its exact torsion gradient, in numpy.

    ``energy_and_grad(phi, psi)`` returns ``(E (B,), dE/dphi (B, n), dE/dpsi (B, n), CA)``.
    """

    def __init__(self, objective):
        self.obj = objective
        self.n = n = objective.n
        self.w_dist = float(objective.w_dist)
        self.w_clash = float(objective.w_clash)
        self.cut = {"CA": float(objective.clash_ca), "CB": float(objective.clash_cb)}
        d = objective.dist
        self.members: List[_FastMember] = []
        self.mu = self.sigma = None
        self.disagree = 0.0
        if d is not None and self.w_dist:
            mem = getattr(d, "members", None)
            if mem:
                self.members = [_FastMember(m, n) for m in mem]
                self.mu = np.asarray(d.mu, float)
                self.sigma = np.asarray(d.sigma, float)
                self.disagree = float(getattr(d, "disagree", 0.0))
                if self.disagree:
                    raise NotImplementedError(
                        "CombinedDistogram.disagree > 0 is not differentiable here")
            else:
                self.members = [_FastMember(d, n)]
        ci, cj = np.triu_indices(n, k=3)
        self.ci, self.cj = ci, cj
        self.cM = _incidence(n, ci, cj)
        self.w_cons = float(getattr(objective, "w_cons", 0.0) or 0.0)
        self.cons_d = getattr(objective, "_cons_d", None)
        self.cons_w = getattr(objective, "_cons_w", None)
        if self.w_cons and self.cons_d is not None:
            self.cons_i = np.asarray(d.i, np.int64)
            self.cons_j = np.asarray(d.j, np.int64)
            self.consM = _incidence(n, self.cons_i, self.cons_j)
        self.n_calls = 0

    # -- energy and coordinate gradient ---------------------------------
    def _pair(self, X, i, j):
        v = X[:, i, :] - X[:, j, :]
        d = np.linalg.norm(v, axis=-1)
        return v, np.maximum(d, 1e-9)

    def energy_and_grad(self, phi: np.ndarray, psi: np.ndarray):
        coords = geo.build_backbone_batch(phi, psi)
        CA, CB, N, C = coords["CA"], coords["CB"], coords["N"], coords["C"]
        B, n, _ = CA.shape
        self.n_calls += 1
        e = np.zeros(B)
        gCA = np.zeros_like(CA)
        gCB = np.zeros_like(CB)

        # distogram
        if self.members:
            vals = np.empty((len(self.members), B))
            for k, m in enumerate(self.members):
                v, d = self._pair(CA, m.i, m.j)
                risk, dd = m.risk_and_dd(d)
                vals[k] = risk
                fac = (self.w_dist if self.mu is None
                       else self.w_dist / (self.sigma[k] * len(self.members)))
                gCA += fac * np.einsum("ap,bpc->bac", m.M, (dd / d)[:, :, None] * v)
            if self.mu is None:
                e += self.w_dist * vals[0]
            else:
                e += self.w_dist * ((vals - self.mu[:, None])
                                    / self.sigma[:, None]).mean(0)

        # clash
        if self.w_clash:
            for key, g in (("CA", gCA), ("CB", gCB)):
                X = CA if key == "CA" else CB
                v, d = self._pair(X, self.ci, self.cj)
                pen = np.maximum(self.cut[key] - d, 0.0)
                e += self.w_clash * pen.sum(1)
                dd = -self.w_clash * (pen > 0.0)
                g += np.einsum("ap,bpc->bac", self.cM, (dd / d)[:, :, None] * v)

        # optional consensus restraint
        if self.w_cons and self.cons_d is not None:
            v, d = self._pair(CA, self.cons_i, self.cons_j)
            err = np.abs(d - self.cons_d[None, :])
            w = 1.0 if self.cons_w is None else self.cons_w[None, :]
            e += self.w_cons * (err * w).mean(1)
            dd = self.w_cons * w * np.sign(d - self.cons_d[None, :]) / len(self.cons_i)
            gCA += np.einsum("ap,bpc->bac", self.consM,
                             (dd / d)[:, :, None] * v)

        dphi, dpsi = self._torsion_grad(N, CA, C, CB, gCA, gCB)
        return e, dphi, dpsi, CA

    # -- coordinate gradient -> torsion gradient ------------------------
    @staticmethod
    def _torsion_grad(N, CA, C, CB, gCA, gCB):
        """Exact dE/d(phi, psi) from dE/d(CA, CB).

        Increasing a torsion is a right-handed rotation about the unit vector from the
        THIRD to the FOURTH atom of its defining frame -- ``CA[i] -> C[i]`` for ``psi[i]``,
        ``N[k] -> CA[k]`` for ``phi[k]`` (read off `protein_geometry._place_atom`, whose
        triad ``(m, nrm, bc)`` is right-handed with ``m x nrm = bc``). Everything the NeRF
        chain places after that bond moves rigidly with the rotation, so

            dE/dtheta = e . sum_{a downstream} (x_a - o) x g_a
                      = e . [ sum x_a x g_a  -  o x sum g_a ] ,

        and both sums over "downstream" are reverse cumulative sums over residues. The set
        downstream of ``psi[i]`` is residues > i; of ``phi[k]`` it is residues > k plus
        ``CB[k]`` (which is built from ``C[k]`` and so turns with it) -- ``CB[i]`` does NOT
        move with ``psi[i]``, which is why the two cases differ by exactly that one term.
        ``phi[0]`` and ``psi[n-1]`` are structurally unused and get zero.
        """
        B, n, _ = CA.shape
        Gres = gCA + gCB
        Ares = np.cross(CA, gCA) + np.cross(CB, gCB)
        Gsuf = np.zeros((B, n + 1, 3))
        Asuf = np.zeros((B, n + 1, 3))
        Gsuf[:, :n] = np.cumsum(Gres[:, ::-1], axis=1)[:, ::-1]
        Asuf[:, :n] = np.cumsum(Ares[:, ::-1], axis=1)[:, ::-1]

        dpsi = np.zeros((B, n))
        if n > 1:
            ax = C[:, :n - 1] - CA[:, :n - 1]
            ax = ax / np.maximum(np.linalg.norm(ax, axis=-1, keepdims=True), 1e-12)
            torque = Asuf[:, 1:n] - np.cross(CA[:, :n - 1], Gsuf[:, 1:n])
            dpsi[:, :n - 1] = (ax * torque).sum(-1)

        dphi = np.zeros((B, n))
        if n > 1:
            ax = CA[:, 1:] - N[:, 1:]
            ax = ax / np.maximum(np.linalg.norm(ax, axis=-1, keepdims=True), 1e-12)
            T = Asuf[:, 2:n + 1] + np.cross(CB[:, 1:], gCB[:, 1:])
            G = Gsuf[:, 2:n + 1] + gCB[:, 1:]
            torque = T - np.cross(CA[:, 1:], G)
            dphi[:, 1:] = (ax * torque).sum(-1)
        return dphi, dpsi

    def energy(self, phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
        return self.energy_and_grad(phi, psi)[0]


# ==================================================================== optimisers
def _adam_np(F: FastObjective, phi: np.ndarray, psi: np.ndarray,
             anchor_phi: np.ndarray, anchor_psi: np.ndarray,
             steps: int, lr: float, trust: Optional[float],
             b1: float = 0.9, b2: float = 0.999, eps: float = 1e-8):
    """Batched Adam. Every row is an independent structure and the loss is their sum, so
    Adam's per-parameter step makes the batch exactly equivalent to B separate runs --
    which is what makes multi-start nearly free."""
    phi, psi = phi.copy(), psi.copy()
    m = [np.zeros_like(phi), np.zeros_like(psi)]
    v = [np.zeros_like(phi), np.zeros_like(psi)]
    best_e = np.full(len(phi), np.inf)
    best_phi, best_psi = phi.copy(), psi.copy()
    for t in range(1, steps + 1):
        e, gp, gs, _ = F.energy_and_grad(phi, psi)
        imp = e < best_e
        if imp.any():
            best_e[imp] = e[imp]
            best_phi[imp] = phi[imp]
            best_psi[imp] = psi[imp]
        for k, (x, g) in enumerate(((phi, gp), (psi, gs))):
            m[k] = b1 * m[k] + (1 - b1) * g
            v[k] = b2 * v[k] + (1 - b2) * g * g
            x -= lr * (m[k] / (1 - b1 ** t)) / (np.sqrt(v[k] / (1 - b2 ** t)) + eps)
        if trust is not None:
            np.clip(phi, anchor_phi - trust, anchor_phi + trust, out=phi)
            np.clip(psi, anchor_psi - trust, anchor_psi + trust, out=psi)
    e, _, _, _ = F.energy_and_grad(phi, psi)
    imp = e < best_e
    if imp.any():
        best_e[imp] = e[imp]
        best_phi[imp] = phi[imp]
        best_psi[imp] = psi[imp]
    return best_phi, best_psi, best_e


def _lbfgs_np(F: FastObjective, phi: np.ndarray, psi: np.ndarray,
              anchor_phi: np.ndarray, anchor_psi: np.ndarray,
              steps: int, lr: float, trust: Optional[float]):
    """scipy L-BFGS-B on the summed batch energy, with the trust region as box bounds.

    The sum decomposes over rows, so its stationary points are exactly the per-row ones;
    only the line search and the curvature pairs are shared, which costs a little
    efficiency and buys a real quasi-Newton step.
    """
    from scipy.optimize import minimize
    B, n = phi.shape
    x0 = np.concatenate([phi.ravel(), psi.ravel()])

    def fun(x):
        p = x[:B * n].reshape(B, n)
        s = x[B * n:].reshape(B, n)
        e, gp, gs, _ = F.energy_and_grad(p, s)
        return float(e.sum()), np.concatenate([gp.ravel(), gs.ravel()])

    bounds = None
    if trust is not None:
        lo = np.concatenate([(anchor_phi - trust).ravel(), (anchor_psi - trust).ravel()])
        hi = np.concatenate([(anchor_phi + trust).ravel(), (anchor_psi + trust).ravel()])
        bounds = list(zip(lo, hi))
    r = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=bounds,
                 options={"maxiter": steps, "maxfun": 4 * steps, "ftol": 1e-14,
                          "gtol": 1e-12})
    p = r.x[:B * n].reshape(B, n)
    s = r.x[B * n:].reshape(B, n)
    e, _, _, _ = F.energy_and_grad(p, s)
    # Never return a row the optimiser made worse.
    e0, _, _, _ = F.energy_and_grad(phi, psi)
    bad = e0 < e
    p[bad], s[bad], e[bad] = phi[bad], psi[bad], e0[bad]
    return p, s, e


def _torch_relax(T: TorchObjective, phi, psi, anchor_phi, anchor_psi,
                 steps, lr, trust, kind: str):
    """The autograd reference path. Kept because it is the oracle `check_analytic` uses,
    and because a second independent implementation is what makes the fast one credible."""
    p = T._t(phi).requires_grad_(True)
    s = T._t(psi).requires_grad_(True)
    ap, asx = T._t(anchor_phi), T._t(anchor_psi)
    if kind == "lbfgs":
        opt = torch.optim.LBFGS([p, s], lr=lr, max_iter=steps, history_size=20,
                                line_search_fn="strong_wolfe")

        def closure():
            opt.zero_grad(set_to_none=True)
            e = T.energy(p, s).sum()
            e.backward()
            return e

        was = T.exact
        T.set_exact(False)          # a line search cannot work on a piecewise-constant value
        try:
            opt.step(closure)
        finally:
            T.set_exact(was)
    else:
        opt = torch.optim.Adam([p, s], lr=lr)
        for _ in range(steps):
            opt.zero_grad(set_to_none=True)
            e = T.energy(p, s)
            e.sum().backward()
            opt.step()
            if trust is not None:
                with torch.no_grad():
                    p.clamp_(ap - trust, ap + trust)
                    s.clamp_(asx - trust, asx + trust)
    with torch.no_grad():
        if trust is not None:
            p.clamp_(ap - trust, ap + trust)
            s.clamp_(asx - trust, asx + trust)
        e = T.energy(p, s).cpu().numpy()
    return p.detach().cpu().numpy(), s.detach().cpu().numpy(), e


# ==================================================================== the operator
def _start_angles(objective, states: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    rows = np.arange(objective.n)
    st = np.asarray(states, int)
    return objective.rep._phi[rows, st].copy(), objective.rep._psi[rows, st].copy()


def refine_angles(objective, phi0: np.ndarray, psi0: np.ndarray,
                  method: str = "adam", n_starts: int = 16, sigma0: float = 12.0,
                  steps: int = 200, lr: float = 0.03,
                  trust: Optional[float] = None, hops: int = 0,
                  hop_sigma: float = 30.0, hop_res: int = 2, seed: int = 0,
                  fast_obj: Optional[FastObjective] = None,
                  torch_obj: Optional[TorchObjective] = None,
                  return_all: bool = False):
    """Relax ``(phi0, psi0)`` -- ``(n,)`` or ``(B, n)`` -- and return the best per row.

    ``return_all=True`` instead returns every restart, shaped ``(B, n_starts, ...)``, which
    is what an oracle experiment needs: it separates "the operator could not find a better
    structure" from "the objective did not prefer the better structure it found".
    """
    use_torch = method.startswith("torch")
    if use_torch:
        T = torch_obj if torch_obj is not None else TorchObjective(objective)
        relax = lambda p, s, ap, asx: _torch_relax(
            T, p, s, ap, asx, steps, lr, tr, "lbfgs" if "lbfgs" in method else "adam")
    else:
        F = fast_obj if fast_obj is not None else FastObjective(objective)
        base = _lbfgs_np if method == "lbfgs" else _adam_np
        relax = lambda p, s, ap, asx: base(F, p, s, ap, asx, steps, lr, tr)

    P0 = np.atleast_2d(np.asarray(phi0, float))
    S0 = np.atleast_2d(np.asarray(psi0, float))
    B, n = P0.shape
    rng = np.random.default_rng(seed)
    k = max(1, int(n_starts))
    sig = math.radians(sigma0)
    dphi = rng.normal(0.0, sig, (B, k, n))
    dpsi = rng.normal(0.0, sig, (B, k, n))
    dphi[:, 0] = 0.0                       # row 0 is the untouched incumbent
    dpsi[:, 0] = 0.0
    phi = (P0[:, None, :] + dphi).reshape(B * k, n)
    psi = (S0[:, None, :] + dpsi).reshape(B * k, n)
    anchor_phi = np.repeat(P0, k, 0)
    anchor_psi = np.repeat(S0, k, 0)
    tr = None if trust is None else math.radians(trust)

    phi, psi, e = relax(phi, psi, anchor_phi, anchor_psi)
    for _ in range(max(0, hops)):
        idx = rng.integers(0, n, size=(B * k, max(1, hop_res)))
        rows = np.repeat(np.arange(B * k), idx.shape[1])
        kp, ks = np.zeros_like(phi), np.zeros_like(psi)
        amp = math.radians(hop_sigma)
        kp[rows, idx.ravel()] = rng.normal(0.0, amp, idx.size)
        ks[rows, idx.ravel()] = rng.normal(0.0, amp, idx.size)
        p2, s2, e2 = relax(phi + kp, psi + ks, anchor_phi, anchor_psi)
        better = e2 < e
        phi[better], psi[better], e[better] = p2[better], s2[better], e2[better]

    # Final ranking on the NUMPY objective, so nothing this returns can be an artefact of
    # the interpolation or of the torch path.
    en, ca = objective.score_angles(phi, psi)
    en = en.reshape(B, k)
    ca = ca.reshape(B, k, n, 3)
    ph = phi.reshape(B, k, n)
    ps = psi.reshape(B, k, n)
    if return_all:
        return ph, ps, en, ca
    pick = en.argmin(1)
    r = np.arange(B)
    out = ph[r, pick], ps[r, pick], en[r, pick], ca[r, pick]
    if np.ndim(phi0) == 1:
        return out[0][0], out[1][0], float(out[2][0]), out[3][0]
    return out


def refine(objective, states: np.ndarray, seed: int = 0, **kw):
    """Drop-in replacement for `foldvqe.refine`: one ``(n,)`` state assignment in,
    ``(phi, psi, energy, ca)`` out."""
    phi0, psi0 = _start_angles(objective, states)
    return refine_angles(objective, phi0, psi0, seed=seed, **kw)


def refine_pool(objective, states_list: Sequence[np.ndarray], seed: int = 0, **kw):
    """Refine many assignments in ONE batch -> ``(phi, psi, e, ca)`` stacked arrays."""
    rows = np.arange(objective.n)
    S = np.asarray(states_list, int)
    return refine_angles(objective, objective.rep._phi[rows[None, :], S],
                         objective.rep._psi[rows[None, :], S], seed=seed, **kw)


# ==================================================================== validation
def check_builder(n: int = 14, B: int = 5, seed: int = 0) -> Dict[str, float]:
    """Max abs difference between the torch and numpy backbone builders."""
    rng = np.random.default_rng(seed)
    phi = rng.uniform(-np.pi, np.pi, (B, n))
    psi = rng.uniform(-np.pi, np.pi, (B, n))
    ref = geo.build_backbone_batch(phi, psi)
    got = build_backbone_torch(torch.as_tensor(phi, dtype=DTYPE),
                               torch.as_tensor(psi, dtype=DTYPE), with_o=True)
    return {k: float(np.abs(got[k].numpy() - ref[k]).max()) for k in ref}


def _random_angles(objective, B, seed, jitter=0.2):
    rng = np.random.default_rng(seed)
    rows = np.arange(objective.n)
    S = rng.integers(0, objective.rep.n_states, (B, objective.n))
    return (objective.rep._phi[rows[None], S] + rng.normal(0, jitter, (B, objective.n)),
            objective.rep._psi[rows[None], S] + rng.normal(0, jitter, (B, objective.n)))


def check_objective(objective, B: int = 32, seed: int = 0) -> Dict[str, float]:
    """`TorchObjective` (grid-snapped) against `FoldObjective.score_angles`."""
    phi, psi = _random_angles(objective, B, seed)
    ref, _ = objective.score_angles(phi, psi)
    got = TorchObjective(objective, exact=True).energy_np(phi, psi)
    return {"max_abs": float(np.abs(got - ref).max()),
            "ref_spread": float(ref.std())}


def check_energy(objective, B: int = 64, seed: int = 0) -> Dict[str, float]:
    """`FastObjective` (interpolated) against `FoldObjective.score_angles`."""
    phi, psi = _random_angles(objective, B, seed)
    ref, _ = objective.score_angles(phi, psi)
    got = FastObjective(objective).energy(phi, psi)
    return {"max_abs": float(np.abs(got - ref).max()),
            "rms": float(np.sqrt(((got - ref) ** 2).mean())),
            "corr": float(np.corrcoef(got, ref)[0, 1]),
            "ref_spread": float(ref.std()),
            "grid_cell": float(objective.dist.grid[1] - objective.dist.grid[0])}


def check_gradient(objective, seed: int = 0, eps: float = 1e-5) -> Dict[str, float]:
    """Torch autograd against a central finite difference on the torch energy."""
    rng = np.random.default_rng(seed)
    n = objective.n
    T = TorchObjective(objective, exact=False)
    phi = torch.as_tensor(rng.uniform(-2, 2, (1, n)), dtype=DTYPE).requires_grad_(True)
    psi = torch.as_tensor(rng.uniform(-2, 2, (1, n)), dtype=DTYPE).requires_grad_(True)
    T.energy(phi, psi).sum().backward()
    g = np.concatenate([phi.grad.numpy().ravel(), psi.grad.numpy().ravel()])
    fd = np.zeros_like(g)
    with torch.no_grad():
        p, s = phi.detach().clone(), psi.detach().clone()
        for a in range(2 * n):
            arr, idx = (p, a) if a < n else (s, a - n)
            arr[0, idx] += eps
            hi = float(T.energy(p, s)[0])
            arr[0, idx] -= 2 * eps
            lo = float(T.energy(p, s)[0])
            arr[0, idx] += eps
            fd[a] = (hi - lo) / (2 * eps)
    return {"max_abs": float(np.abs(g - fd).max()),
            "rel": float(np.abs(g - fd).max() / max(np.abs(fd).max(), 1e-12))}


def check_analytic(objective, B: int = 8, seed: int = 0) -> Dict[str, float]:
    """The closed-form numpy torsion gradient against torch autograd.

    The two share no code: one differentiates the NeRF chain op by op, the other applies
    the rigid-rotation identity to a coordinate gradient computed by hand.
    """
    phi, psi = _random_angles(objective, B, seed, jitter=0.4)
    F = FastObjective(objective)
    _, gp, gs, _ = F.energy_and_grad(phi, psi)
    T = TorchObjective(objective, exact=False)
    p = T._t(phi).requires_grad_(True)
    s = T._t(psi).requires_grad_(True)
    T.energy(p, s).sum().backward()
    rp = p.grad.numpy()
    rs = s.grad.numpy()
    # phi[0] and psi[n-1] are structurally unused; autograd reports the same zeros.
    num = np.concatenate([gp.ravel(), gs.ravel()])
    ref = np.concatenate([rp.ravel(), rs.ravel()])
    scale = max(np.abs(ref).max(), 1e-12)
    return {"max_abs": float(np.abs(num - ref).max()),
            "rel": float(np.abs(num - ref).max() / scale),
            "cos": float((num @ ref) / (np.linalg.norm(num) * np.linalg.norm(ref))),
            "grad_scale": float(scale)}

if __name__ == "__main__":
    print("builder:", check_builder())
