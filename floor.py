"""Representation floor: the best CA-RMSD an encoding can express for a given native.

Nearest-state projection followed by batched coordinate descent, then a short annealed
polish. Every candidate for one sweep is built in a single `build_backbone_batch` call,
so a sweep costs one batched geometry build rather than n*k scalar ones.

This is a property of the *encoding*, measured without any search or energy model. If the
floor for a target is above 2 A, no optimiser can put that target under 2 A and the
representation is the binding constraint.
"""
from typing import Optional, Tuple

import numpy as np

import protein_geometry as geo


def _rmsd_batch(phi: np.ndarray, psi: np.ndarray, native_ca: np.ndarray) -> np.ndarray:
    ca = geo.build_backbone_batch(phi, psi)["CA"]
    return geo.ca_rmsd_batch(ca, native_ca)


def project(rep, native_phi: np.ndarray, native_psi: np.ndarray) -> np.ndarray:
    """Per-residue nearest state by angular distance, ignoring the undefined termini."""
    n = rep.n_residues
    d = np.zeros((n, rep.n_states))
    for i in range(n):
        if i > 0:
            d[i] += (np.angle(np.exp(1j * (rep._phi[i] - native_phi[i])))) ** 2
        if i < n - 1:
            d[i] += (np.angle(np.exp(1j * (rep._psi[i] - native_psi[i])))) ** 2
    return d.argmin(1)


def descend(rep, native_ca: np.ndarray, states: np.ndarray,
            sweeps: int = 12) -> Tuple[np.ndarray, float]:
    n, k = rep.n_residues, rep.n_states
    rows = np.arange(n)
    cur = np.array(states, int)
    best = float(_rmsd_batch(rep._phi[rows, cur][None], rep._psi[rows, cur][None],
                             native_ca)[0])
    for _ in range(sweeps):
        improved = False
        for i in range(n):
            cand = np.repeat(cur[None], k, axis=0)
            cand[:, i] = np.arange(k)
            r = _rmsd_batch(rep._phi[rows[None], cand], rep._psi[rows[None], cand],
                            native_ca)
            j = int(r.argmin())
            if r[j] < best - 1e-9:
                best, cur, improved = float(r[j]), cand[j], True
        if not improved:
            break
    return cur, best


def floor(rep, native_ca: np.ndarray, native_phi: np.ndarray, native_psi: np.ndarray,
          restarts: int = 6, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    n, k = rep.n_residues, rep.n_states
    s0 = project(rep, native_phi, native_psi)
    rows = np.arange(n)
    proj = float(_rmsd_batch(rep._phi[rows, s0][None], rep._psi[rows, s0][None],
                             native_ca)[0])
    best_s, best = descend(rep, native_ca, s0)
    dproj = best
    for r in range(restarts):
        start = np.where(rng.random(n) < 0.35, rng.integers(0, k, n), best_s)
        s, v = descend(rep, native_ca, start)
        if v < best:
            best, best_s = v, s
    return {"projection": proj, "descent_from_projection": dproj,
            "floor": best, "states": best_s}
