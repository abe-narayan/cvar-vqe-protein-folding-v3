"""SPRINT 14, coordinator -- the retrieval-conditioned torsion prior.

WHY THIS EXISTS.  Sprint 13's torsion prior (`s13.qarch_lib.empirical_prior`) is built from
the residue-CLASS back-off pool: every GENERAL residue in the held-out database gets the same
distribution.  That is deliberately conservative, and Sprint 13 then concluded that 88% of
what the library buys is generic Ramachandran and only 12% is sequence-related.

But the retrieval half of the project already holds something strictly stronger and nobody
has ever fed it to the torsion half: each target's universe cache stores `PHI`/`PSI` for
every retrieved window, so the K=500 BLOSUM pool is a POSITION-SPECIFIC, SEQUENCE-CONDITIONED
empirical distribution over (phi_i, psi_i) -- the exact object a VQE prior Hamiltonian wants,
and the natural Level 2 of the sprint brief's information ladder.

This module builds it, prices it as an emitter on the shared `s14.ladder` harness, and --
the number that actually decides its worth -- measures its per-torsion angular error, so it
lands on the sigma scale where the incumbent pipeline is sigma ~29 deg and sigma 12 deg
builds a 1.486 A structure.

NATIVE-FREE.  `order`/`sim` are BLOSUM retrieval outputs; `rr` (the oracle RMSD column) and
`nat_ca` are touched ONLY inside functions whose name says ORACLE, and only post hoc.  Every
emitter here is checked by `s14.ladder.audit_emitter`, which poisons the native and demands
bit-identical output.

Run:
    python -m s14.retprior
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s14 import ladder as L                # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def circ_mean(A, w=None, axis=0):
    """Circular mean of angles A along `axis`, optional weights."""
    A = np.asarray(A, float)
    if w is None:
        c, s = np.cos(A).mean(axis), np.sin(A).mean(axis)
    else:
        w = np.asarray(w, float)
        sh = [1] * A.ndim
        sh[axis] = -1
        w = w.reshape(sh)
        tot = w.sum(axis)
        c = (np.cos(A) * w).sum(axis) / tot
        s = (np.sin(A) * w).sum(axis) / tot
    return np.arctan2(s, c)


def circ_R(A, w=None, axis=0):
    """Mean resultant length in [0,1].  1 = all windows agree, 0 = uniform. A concentration."""
    A = np.asarray(A, float)
    if w is None:
        c, s = np.cos(A).mean(axis), np.sin(A).mean(axis)
    else:
        w = np.asarray(w, float)
        sh = [1] * A.ndim
        sh[axis] = -1
        w = w.reshape(sh)
        tot = w.sum(axis)
        c = (np.cos(A) * w).sum(axis) / tot
        s = (np.sin(A) * w).sum(axis) / tot
    return np.hypot(c, s)


# ------------------------------------------------------------------ the retrieved windows
def windows(pdb: str, arm: str = "pool"):
    """Native-free (m, n) phi/psi of the retrieved windows, plus BLOSUM weights.

    arm = "pool"   the shipped K=500 BLOSUM pool
          "top75"  the pool rows the shipped distogram filter keeps
          "top20"  the 20 highest-similarity pool rows
    """
    u = I.load_univ(pdb)
    p = I.pool_idx(u)
    if arm == "top75":
        rec = I.shipped_record(pdb)
        sub = np.asarray(rec["sub"], int)          # indices INTO THE POOL
        p = p[sub]
    elif arm == "top20":
        p = p[:20]
    elif arm != "pool":
        raise ValueError(arm)
    sim = np.asarray(u["sim"], float)[p]
    return np.asarray(u["PHI"], float)[p], np.asarray(u["PSI"], float)[p], sim


def _softmax_w(sim, tau):
    if tau is None:
        return None
    z = (sim - sim.max()) / max(tau, 1e-9)
    return np.exp(z)


# -------------------------------------------------------------------------------- emitters
def _make_circmean(arm, tau=None):
    def emit(pdb, seq, n, fold, rng):
        PHI, PSI, sim = windows(pdb, arm)
        w = _softmax_w(sim, tau)
        return circ_mean(PHI, w, axis=0), circ_mean(PSI, w, axis=0)
    return emit


def _make_state_argmax(arm, k=4, tau=None):
    """Assign every retrieved window's (phi_i, psi_i) to its nearest library state, then
    take the per-residue argmax.  This is the retrieval-conditioned analogue of
    `empirical_prior`, over the SAME k-state space the quantum encoding uses."""
    def emit(pdb, seq, n, fold, rng):
        sp = L.space(pdb, k)
        PHI, PSI, sim = windows(pdb, arm)
        w = _softmax_w(sim, tau)
        P = _state_counts(sp, PHI, PSI, w)
        S = np.argmax(P, axis=1)
        r = np.arange(sp.n)
        return sp.PHI[r, S], sp.PSI[r, S]
    return emit


def _state_counts(sp, PHI, PSI, w=None):
    """(n, k) Laplace-smoothed occupancy of each library state under the retrieved windows.

    Nearest state in the same 4-D (cos/sin phi, cos/sin psi) embedding the library was
    clustered in, which is the embedding `empirical_prior` uses -- so the two priors are
    directly comparable and differ only in WHICH observations they count.
    """
    n, k = sp.n, sp.k
    P = np.zeros((n, k))
    m = PHI.shape[0]
    ww = np.ones(m) if w is None else np.asarray(w, float)
    for i in range(n):
        X = np.column_stack([np.cos(PHI[:, i]), np.sin(PHI[:, i]),
                             np.cos(PSI[:, i]), np.sin(PSI[:, i])])
        C = np.column_stack([np.cos(sp.PHI[i]), np.sin(sp.PHI[i]),
                             np.cos(sp.PSI[i]), np.sin(sp.PSI[i])])
        lab = ((X[:, None, :] - C[None]) ** 2).sum(-1).argmin(1)
        cnt = np.bincount(lab, weights=ww, minlength=k).astype(float)
        P[i] = (cnt + 1.0) / (cnt.sum() + k)
    return P


def state_prior(pdb, arm="top75", k=4, tau=None):
    """The deliverable: a native-free (n, k) position-specific torsion prior.

    This is what an uncertainty-aware prior Hamiltonian encodes.  Returned with its
    per-residue concentration so a downstream penalty can be weighted by confidence.
    """
    sp = L.space(pdb, k)
    PHI, PSI, sim = windows(pdb, arm)
    w = _softmax_w(sim, tau)
    P = _state_counts(sp, PHI, PSI, w)
    return {"P": P, "R_phi": circ_R(PHI, w, axis=0), "R_psi": circ_R(PSI, w, axis=0),
            "entropy": -(P * np.log(P)).sum(1) / np.log(k), "m": int(PHI.shape[0])}


EMITTERS = {
    "L2a_pool500_circmean": (_make_circmean("pool"), False),
    "L2b_top75_circmean": (_make_circmean("top75"), False),
    "L2c_top20_circmean": (_make_circmean("top20"), False),
    "L2d_top75_simweighted": (_make_circmean("top75", tau=2.0), False),
    "L2e_top75_state_argmax": (_make_state_argmax("top75"), False),
    "L2f_pool500_state_argmax": (_make_state_argmax("pool"), False),
}


# ------------------------------------------------------------- ORACLE: place it on sigma
def ORACLE_angular_error(emit, targets=None):
    """POST-HOC ONLY.  Mean absolute phi/psi error of an emitter against native torsions.

    Places an emitter on the sigma scale the restraint surface is indexed by: the incumbent
    3.204 A pipeline corresponds to sigma ~29 deg, and sigma 12 deg builds 1.486 A.
    """
    import peptide_db as pdbm
    tg = targets if targets is not None else I.targets()
    ephi, epsi, rows = [], [], []
    for t in tg:
        p = pdbm.by_pdb(t["pdb"])
        nphi, npsi = np.asarray(p.phi, float), np.asarray(p.psi, float)
        rng = np.random.default_rng(0)
        phi, psi = emit(t["pdb"], t["seq"], int(t["n"]), int(t["fold"]), rng)
        # residues whose torsion is undefined at a terminus are excluded from the average
        mphi = np.abs(wrap(phi - nphi))[1:]
        mpsi = np.abs(wrap(psi - npsi))[:-1]
        ephi.append(np.rad2deg(mphi).mean()); epsi.append(np.rad2deg(mpsi).mean())
        rows.append({"pdb": t["pdb"], "phi": float(np.rad2deg(mphi).mean()),
                     "psi": float(np.rad2deg(mpsi).mean())})
    return {"mae_phi_deg": float(np.mean(ephi)), "mae_psi_deg": float(np.mean(epsi)),
            "mae_mean_deg": float(np.mean(ephi) / 2 + np.mean(epsi) / 2), "per_target": rows}


def run():
    L.EMITTERS.update(EMITTERS)
    out = L.run(names=list(EMITTERS))

    # concentration of the prior, and where it lands on the sigma scale
    extra = {}
    for name, (emit, _) in EMITTERS.items():
        extra[name] = ORACLE_angular_error(emit)
        print(f"  {name:<26} ORACLE angular error  phi {extra[name]['mae_phi_deg']:6.1f} deg"
              f"   psi {extra[name]['mae_psi_deg']:6.1f} deg")

    conc = []
    for t in I.targets():
        sp0 = state_prior(t["pdb"], "top75")
        conc.append({"pdb": t["pdb"], "R_phi": float(sp0["R_phi"].mean()),
                     "R_psi": float(sp0["R_psi"].mean()),
                     "entropy": float(sp0["entropy"].mean())})
    out["angular_error"] = extra
    out["top75_concentration"] = conc
    out["reference_scale"] = {
        "incumbent_equivalent_sigma_deg": 29.0,
        "sigma12_full_coverage_rmsd": 1.486,
        "best_sequence_only_sigma_deg": 67.7,
        "sequence_blind_marginal_phi_deg": 36.4,
        "sequence_blind_marginal_psi_deg": 72.8,
        "full_context_phi_deg": 36.1,
        "full_context_psi_deg": 62.4,
    }
    with open(os.path.join(RESULTS, "retprior.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    o = run()
    print()
    L._fmt(o)
    c = o["top75_concentration"]
    print(f"\ntop-75 prior concentration: mean R_phi {np.mean([x['R_phi'] for x in c]):.3f}"
          f"  R_psi {np.mean([x['R_psi'] for x in c]):.3f}"
          f"  normalised state entropy {np.mean([x['entropy'] for x in c]):.3f}")
