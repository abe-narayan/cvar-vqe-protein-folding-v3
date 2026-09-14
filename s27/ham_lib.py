#!/usr/bin/env python
"""s27/ham_lib.py -- the Sprint 27 library of alternative Hamiltonians (energy channels).

Every channel is a NATIVE-FREE function of one target's candidate pool that returns one
energy per candidate, LOWER IS BETTER, on the shipped K = 500 BLOSUM pool (`s24.d_harness
.Candidates.from_universe`).  Nothing here reads `nat_ca` or `oracle_rr`; the ORACLE
diagnostics live in `s27/run_singles.py`, in functions whose names start `oracle_`.

Where a channel is STATISTICAL it is fitted, per target, on that target's own leakage-safe
window universe (`s8/generate_univ/<pdb>.npz`: every peptide of the other four folds plus
the fold's identity-filtered fragments, the same rule the shipped distogram trains under), so
the fit never sees the target's own native or its fold.  Where a channel is PHYSICAL it uses
the ideal-geometry backbone rebuilt from the candidate's phi/psi (`core.geometry
.build_backbone_batch`), exactly as the S25 suite did for Legacy and AMBER.  Where a channel
is a CONSISTENCY term it uses only the pool itself.

Families and channels (the name in the results files is the key):

  geometric / compactness
    RG_LAW     |Rg - 2.2 n^0.38| / 2.2 n^0.38          compactness against the folded-protein
                                                       scaling law (Rg ~ N^0.34-0.38)
    RG_UNIV    |Rg - median Rg of the universe|          the same, data-driven target
    EXVOL      CA-CA excluded volume, |i-j| >= 3, d < 4.0 A, soft quadratic
    CAGEO      -log p(theta_i, tau_i): CA virtual-angle / virtual-torsion statistics of the
               universe (Levitt-style CA-trace potential)
  torsion
    RAMA       -sum_i log p_fold(aa_i, phi_i, psi_i), the per-fold leakage-safe table
               `s8/generate_rama.npz` (10-degree bins, pseudocount)
  statistical residue-pair (fitted on the universe, quasi-chemical reference state)
    CONTACT    Miyazawa-Jernigan-style contact potential e(a,b) = -log(N_obs/N_exp) on CA-CA
               contacts < 7.5 A, |i-j| >= 3
    DISTPOT    Sippl / DFIRE-style distance-dependent pair potential e(a,b,s,d) with 1 A bins
               3..14 A, separation classes {3, 4-5, >= 6}, reference = separation-only
    ENV        Rosetta-env-like one-body burial term: -log p(nbr bin | aa) / p(nbr bin),
               neighbours = CA within 8 A, |i-j| >= 2
    HP         hydrophobic burial, parameter-free: -sum_i KD(aa_i) * nbr_i (Kyte-Doolittle)
  backbone physics on the ideal rebuild
    DSSPHB     Kabsch-Sander backbone H-bond energy, sum of E < -0.5 kcal/mol, |i-j| >= 2,
               H placed 1.0 A from N opposite the previous C=O
    ELEC       Debye-Hueckel screened Coulomb between charged side-chain proxies at CB
               (K,R +1; D,E -1; H +0.5), eps 80, kappa 0.1 / A, |i-j| >= 2
  pool / ensemble consistency
    CONS       mean CA-RMSD to the other pool members (typicality; the S12 medoid criterion)
    DMAP_CONS  L1 distance of the candidate's distance map to the pool's median distance map
    TORS_CONS  summed circular distance of (phi, psi) to the pool's per-residue circular mean
    POOLGO     Go-like reward for the pool's consensus contacts: -sum_ij f_ij [d_ij < 8],
               f_ij = pool contact frequency, |i-j| >= 3
  sequence-structure compatibility
    SS_MATCH   mismatch between a CA-geometry secondary-structure call and the sequence's
               Chou-Fasman helix / strand propensity
  distogram re-readings (controls for "another functional of the same posterior")
    CONTACT_LL -log-likelihood of the candidate's 8 A contact map under the distogram's
               contact probabilities
    DIS_MEAN   L1 error to the posterior MEAN distance (the shipped score is the L1 Bayes
               risk, i.e. the median), unit weights
  reference channels (the S25 suite's three, read from their caches)
    DIS        the shipped Bayes-risk score      LEG   Legacy 11-term     AMB   AMBER single point
    LEG_<term> the eleven Legacy terms one at a time (`s16.energy_lib.legacy_components_of_windows`)
"""
from __future__ import annotations

import math
import os
import sys
from typing import Dict, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import data as cdata            # noqa: E402
from core import geometry as geo           # noqa: E402
from s12 import instrument as I            # noqa: E402

ALPHABET = cdata.ALPHABET                  # "ARNDCQEGHILKMFPSTWYV"
AIDX = {a: i for i, a in enumerate(ALPHABET)}

#: Kyte-Doolittle hydropathy, ALPHABET order
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
      "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
      "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
#: Chou-Fasman helix / strand propensities (Chou & Fasman 1978, Table)
CF_A = {"E": 1.51, "M": 1.45, "A": 1.42, "L": 1.21, "K": 1.16, "F": 1.13, "Q": 1.11, "W": 1.08,
        "I": 1.08, "V": 1.06, "D": 1.01, "H": 1.00, "R": 0.98, "T": 0.83, "S": 0.77, "C": 0.70,
        "Y": 0.69, "N": 0.67, "P": 0.57, "G": 0.57}
CF_B = {"V": 1.70, "I": 1.60, "Y": 1.47, "F": 1.38, "W": 1.37, "L": 1.30, "C": 1.19, "T": 1.19,
        "Q": 1.10, "M": 1.05, "R": 0.93, "N": 0.89, "H": 0.87, "A": 0.83, "S": 0.75, "G": 0.75,
        "K": 0.74, "P": 0.55, "D": 0.54, "E": 0.37}
CHARGE = {"K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0, "H": 0.5}

CONTACT_CUT = 7.5     # A, CA-CA, statistical contact potential
NBR_CUT = 8.0         # A, neighbour count for burial
DBINS = np.arange(3.0, 15.0, 1.0)     # distance-dependent potential bins, 3..14 A
PSEUDO = 1.0


# ============================================================================ helpers
def pair_index(n, min_sep=2):
    return I.pair_index(n, min_sep)


def pair_dists(W, i, j):
    return I.pair_dists(W, i, j)


def rg(W):
    """(k,) radius of gyration of CA clouds."""
    W = np.asarray(W, float)
    c = W.mean(1, keepdims=True)
    return np.sqrt(((W - c) ** 2).sum(-1).mean(-1))


def codes(seq):
    return np.array([AIDX.get(a, 0) for a in seq], int)


def dihedral(p0, p1, p2, p3):
    """Batched dihedral (..., 3) -> (...) radians."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-9)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.arctan2(y, x)


def vangle(a, b, c):
    """Angle at b, batched."""
    u = a - b
    v = c - b
    cos = (u * v).sum(-1) / np.maximum(np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1), 1e-9)
    return np.arccos(np.clip(cos, -1, 1))


def ca_geometry(W):
    """CA virtual bond angles theta (k, n-2) and virtual torsions tau (k, n-3)."""
    W = np.asarray(W, float)
    th = vangle(W[:, :-2], W[:, 1:-1], W[:, 2:])
    ta = dihedral(W[:, :-3], W[:, 1:-2], W[:, 2:-1], W[:, 3:])
    return th, ta


# ============================================================================ per-target context
class Context:
    """Everything a channel may read for one target: the pool, its ideal rebuild, the
    leakage-safe universe, the distogram, the fold's Ramachandran table.  Built once."""

    def __init__(self, cand, universe, dg, rama_cnt):
        self.cand = cand
        self.seq = cand.seq
        self.n = cand.n
        self.fold = cand.fold
        self.k = cand.k
        self.W = np.asarray(cand.W, float)                    # (k, n, 3)
        self.PHI = np.asarray(cand.PHI, float)
        self.PSI = np.asarray(cand.PSI, float)
        self.aa = codes(cand.seq)
        self.i, self.j = pair_index(self.n, 2)
        self.D = pair_dists(self.W, self.i, self.j)          # (k, npairs) pool CA distances
        self.sep = (self.j - self.i)
        # universe (float32 -> float64 once)
        self.UW = np.asarray(universe["W"], float)             # (U, n, 3)
        self.US = np.asarray(universe["S"], int)               # (U, n) residue codes
        self.UPHI = np.asarray(universe["PHI"], float)
        self.UPSI = np.asarray(universe["PSI"], float)
        self.dg = dg
        self.rama_cnt = rama_cnt                               # (20, 36, 36) counts, this fold
        self._bb = None
        self._UD = None

    @property
    def bb(self):
        """Ideal-geometry backbone of every candidate: dict of (k, n, 3) N, CA, C, O, CB."""
        if self._bb is None:
            self._bb = {k: np.asarray(v, float) for k, v in
                        geo.build_backbone_batch(self.PHI, self.PSI).items()}
        return self._bb

    @property
    def UD(self):
        """Universe pair distances (U, npairs), computed once, float32 to save memory."""
        if self._UD is None:
            self._UD = pair_dists(self.UW, self.i, self.j).astype(np.float32)
        return self._UD


# ============================================================================ the channels
def h_rg_law(cx: Context):
    r0 = 2.2 * cx.n ** 0.38
    return np.abs(rg(cx.W) - r0) / r0


def h_rg_univ(cx: Context):
    r0 = float(np.median(rg(cx.UW)))
    return np.abs(rg(cx.W) - r0) / r0


def h_exvol(cx: Context, cut=4.0):
    m = cx.sep >= 3
    d = cx.D[:, m]
    return (np.maximum(0.0, cut - d) ** 2).sum(1)


def h_cageo(cx: Context, nb_th=24, nb_ta=36):
    """-log p(theta, tau) from the universe's CA-trace statistics, summed along the chain."""
    uth, uta = ca_geometry(cx.UW)
    th, ta = ca_geometry(cx.W)
    # joint (theta_i, tau_i) for i = 0..n-4 uses theta of the same start index
    H = np.zeros((nb_th, nb_ta))
    eth = np.linspace(0, math.pi, nb_th + 1)
    eta = np.linspace(-math.pi, math.pi, nb_ta + 1)
    ti = np.clip(np.digitize(uth[:, :-1].ravel(), eth) - 1, 0, nb_th - 1)
    tj = np.clip(np.digitize(uta.ravel(), eta) - 1, 0, nb_ta - 1)
    np.add.at(H, (ti, tj), 1.0)
    P = (H + PSEUDO) / (H + PSEUDO).sum()
    logp = np.log(P)
    ci = np.clip(np.digitize(th[:, :-1], eth) - 1, 0, nb_th - 1)
    cj = np.clip(np.digitize(ta, eta) - 1, 0, nb_ta - 1)
    return -logp[ci, cj].sum(1)


def h_rama(cx: Context):
    """-sum_i log p_fold(aa_i, phi_i, psi_i) from the per-fold leakage-safe table."""
    cnt = cx.rama_cnt.astype(float) + PSEUDO                   # (20, 36, 36)
    P = cnt / cnt.sum(axis=(1, 2), keepdims=True)
    logp = np.log(P)
    b = 36
    pb = np.clip(((cx.PHI + math.pi) / (2 * math.pi) * b).astype(int), 0, b - 1)
    sb = np.clip(((cx.PSI + math.pi) / (2 * math.pi) * b).astype(int), 0, b - 1)
    aa = cx.aa[None, :].repeat(cx.k, 0)
    return -logp[aa, pb, sb].sum(1)


def _pair_types(cx: Context):
    """(npairs,) symmetric pair-type index for the pool's sequence, and the universe's."""
    return None


def h_contact(cx: Context, cut=CONTACT_CUT):
    """Quasi-chemical contact potential fitted on the universe: e(a,b) = -log N_obs/N_exp."""
    m = cx.sep >= 3
    ui, uj = cx.i[m], cx.j[m]
    UD = cx.UD[:, m]
    C = UD < cut                                               # (U, P)
    A = cx.US[:, ui]; B = cx.US[:, uj]
    N = np.zeros((20, 20))
    np.add.at(N, (A[C], B[C]), 1.0)
    N = N + N.T
    tot = N.sum() / 2.0
    f = N.sum(1) / (2.0 * tot)                                 # residue frequency among contacts
    Nexp = 2.0 * tot * np.outer(f, f)
    e = -np.log((N + PSEUDO) / (Nexp + PSEUDO))
    e = 0.5 * (e + e.T)
    a = cx.aa[ui]; b = cx.aa[uj]
    epair = e[a, b]                                            # (P,)
    Cc = cx.D[:, m] < cut
    return (Cc * epair[None, :]).sum(1)


def h_distpot(cx: Context):
    """Distance-dependent statistical pair potential (Sippl / DFIRE form, CA level)."""
    m = cx.sep >= 3
    ui, uj = cx.i[m], cx.j[m]
    sepc = np.where(cx.sep[m] == 3, 0, np.where(cx.sep[m] <= 5, 1, 2))   # separation class
    UD = cx.UD[:, m]
    db = np.clip(np.digitize(UD, DBINS) - 1, 0, len(DBINS) - 1)          # (U, P) bins
    A = cx.US[:, ui]; B = cx.US[:, uj]
    # pair type index symmetric: t = min*20 + max
    T = np.minimum(A, B) * 20 + np.maximum(A, B)                           # (U, P)
    nb = len(DBINS)
    H = np.zeros((400, 3, nb))
    Href = np.zeros((3, nb))
    S = np.broadcast_to(sepc[None, :], T.shape)
    np.add.at(H, (T.ravel(), S.ravel(), db.ravel()), 1.0)
    np.add.at(Href, (S.ravel(), db.ravel()), 1.0)
    P = (H + PSEUDO) / (H + PSEUDO).sum(-1, keepdims=True)
    Pref = (Href + PSEUDO) / (Href + PSEUDO).sum(-1, keepdims=True)
    e = -np.log(P / Pref[None])                                            # (400, 3, nb)
    a = cx.aa[ui]; b = cx.aa[uj]
    t = np.minimum(a, b) * 20 + np.maximum(a, b)
    cb = np.clip(np.digitize(cx.D[:, m], DBINS) - 1, 0, nb - 1)          # (k, P)
    return e[t[None, :], sepc[None, :], cb].sum(1)


def _nbr_count(W, i, j, sep, cut=NBR_CUT, min_sep=2):
    """(k, n) number of CA neighbours within cut, |i-j| >= min_sep."""
    D = pair_dists(W, i, j)
    m = (sep >= min_sep)
    C = (D[:, m] < cut).astype(float)
    n = W.shape[1]
    out = np.zeros((W.shape[0], n))
    np.add.at(out.T, i[m], C.T)
    np.add.at(out.T, j[m], C.T)
    return out


def h_env(cx: Context):
    """Rosetta-env-like burial statistics: -log p(bin | aa) / p(bin), fitted on the universe."""
    edges = np.array([0, 2, 4, 6, 99])
    un = _nbr_count(cx.UW, cx.i, cx.j, cx.sep)
    ub = np.clip(np.digitize(un, edges) - 1, 0, 3)
    H = np.zeros((20, 4))
    np.add.at(H, (cx.US.ravel(), ub.ravel()), 1.0)
    Pab = (H + PSEUDO) / (H + PSEUDO).sum(1, keepdims=True)
    Pb = (H.sum(0) + PSEUDO) / (H.sum() + 4 * PSEUDO)
    e = -np.log(Pab / Pb[None, :])
    cn = _nbr_count(cx.W, cx.i, cx.j, cx.sep)
    cb = np.clip(np.digitize(cn, edges) - 1, 0, 3)
    aa = cx.aa[None, :].repeat(cx.k, 0)
    return e[aa, cb].sum(1)


def h_hp(cx: Context):
    kd = np.array([KD[a] for a in cx.seq])
    cn = _nbr_count(cx.W, cx.i, cx.j, cx.sep)
    return -(kd[None, :] * cn).sum(1)


def h_dssphb(cx: Context, thresh=-0.5):
    """Kabsch-Sander backbone H-bond energy on the ideal rebuild, sum of bonds below -0.5."""
    bb = cx.bb
    N, C, O = bb["N"], bb["C"], bb["O"]
    k, n = N.shape[:2]
    # H_i for i >= 1: N_i + unit(C_{i-1} - O_{i-1}) * 1.0
    u = C[:, :-1] - O[:, :-1]
    u = u / np.maximum(np.linalg.norm(u, axis=-1, keepdims=True), 1e-9)
    Hh = np.full_like(N, np.nan)
    Hh[:, 1:] = N[:, 1:] + u
    E = np.zeros(k)
    q = 0.084 * 332.0
    for a in range(1, n):            # donor N-H of residue a
        for b in range(n):           # acceptor C=O of residue b
            if abs(a - b) < 2:
                continue
            rON = np.linalg.norm(O[:, b] - N[:, a], axis=-1)
            rCH = np.linalg.norm(C[:, b] - Hh[:, a], axis=-1)
            rOH = np.linalg.norm(O[:, b] - Hh[:, a], axis=-1)
            rCN = np.linalg.norm(C[:, b] - N[:, a], axis=-1)
            e = q * (1 / np.maximum(rON, 0.5) + 1 / np.maximum(rCH, 0.5)
                     - 1 / np.maximum(rOH, 0.5) - 1 / np.maximum(rCN, 0.5))
            E += np.where(e < thresh, e, 0.0)
    return E


def h_elec(cx: Context, eps=80.0, kappa=0.1):
    q = np.array([CHARGE.get(a, 0.0) for a in cx.seq])
    if np.count_nonzero(q) < 2:
        return np.zeros(cx.k)
    CB = cx.bb["CB"]
    m = cx.sep >= 2
    ii, jj = cx.i[m], cx.j[m]
    qq = q[ii] * q[jj]
    keep = qq != 0
    if not keep.any():
        return np.zeros(cx.k)
    ii, jj, qq = ii[keep], jj[keep], qq[keep]
    r = np.linalg.norm(CB[:, ii] - CB[:, jj], axis=-1)
    r = np.maximum(r, 3.0)
    return (332.0 * qq[None, :] * np.exp(-kappa * r) / (eps * r)).sum(1)


def h_cons(cx: Context):
    """Mean CA-RMSD to the other pool members (typicality)."""
    aud = __import__("core").backend("numerics")
    k = cx.k
    out = np.zeros(k)
    for a in range(k):
        r = aud.kabsch_rmsd_batch(cx.W, cx.W[a])
        out[a] = (r.sum() - r[a]) / (k - 1)
    return out


def h_dmap_cons(cx: Context):
    med = np.median(cx.D, axis=0)
    return np.abs(cx.D - med[None, :]).mean(1)


def h_tors_cons(cx: Context):
    def cmean(x):
        return np.arctan2(np.sin(x).mean(0), np.cos(x).mean(0))
    mphi, mpsi = cmean(cx.PHI), cmean(cx.PSI)
    dphi = np.abs(np.angle(np.exp(1j * (cx.PHI - mphi[None, :]))))
    dpsi = np.abs(np.angle(np.exp(1j * (cx.PSI - mpsi[None, :]))))
    return (dphi + dpsi).sum(1)


def h_poolgo(cx: Context, cut=8.0):
    m = cx.sep >= 3
    C = (cx.D[:, m] < cut).astype(float)
    f = C.mean(0)
    return -(C * f[None, :]).sum(1)


def h_ss_match(cx: Context):
    """CA-geometry secondary-structure call vs Chou-Fasman sequence propensity."""
    W = cx.W
    n = cx.n
    pa = np.array([CF_A[a] for a in cx.seq]); pb = np.array([CF_B[a] for a in cx.seq])
    # normalised propensities in [0, 1]
    pa_n = (pa - 0.5) / 1.1; pb_n = (pb - 0.3) / 1.5
    d3 = np.linalg.norm(W[:, 3:] - W[:, :-3], axis=-1)        # (k, n-3)
    d2 = np.linalg.norm(W[:, 2:] - W[:, :-2], axis=-1)        # (k, n-2)
    helix = np.zeros((cx.k, n)); strand = np.zeros((cx.k, n))
    h = (d3 > 4.5) & (d3 < 5.7)
    for s in range(4):
        helix[:, s:s + d3.shape[1]] += h
    helix = helix > 0
    st = d2 > 6.2
    for s in range(3):
        strand[:, s:s + d2.shape[1]] += st
    strand = (strand > 0) & ~helix
    return (helix * (1 - pa_n)[None, :] + strand * (1 - pb_n)[None, :]).sum(1)


def h_contact_ll(cx: Context, cut=8.0):
    dg = cx.dg
    centres = np.asarray(dg["centres"], float)
    pc = np.asarray(dg["prob"], float)[:, centres < cut].sum(1)
    pc = np.clip(pc, 1e-4, 1 - 1e-4)
    C = (cx.D < cut)
    return -(C * np.log(pc)[None, :] + (~C) * np.log(1 - pc)[None, :]).sum(1)


def h_dis_mean(cx: Context):
    mu = np.asarray(cx.dg["expected"], float)
    return np.abs(cx.D - mu[None, :]).mean(1)


CHANNELS = {
    "RG_LAW": h_rg_law, "RG_UNIV": h_rg_univ, "EXVOL": h_exvol, "CAGEO": h_cageo,
    "RAMA": h_rama, "CONTACT": h_contact, "DISTPOT": h_distpot, "ENV": h_env, "HP": h_hp,
    "DSSPHB": h_dssphb, "ELEC": h_elec, "CONS": h_cons, "DMAP_CONS": h_dmap_cons,
    "TORS_CONS": h_tors_cons, "POOLGO": h_poolgo, "SS_MATCH": h_ss_match,
    "CONTACT_LL": h_contact_ll, "DIS_MEAN": h_dis_mean,
}
FAMILY = {
    "RG_LAW": "geometric", "RG_UNIV": "geometric", "EXVOL": "geometric", "CAGEO": "geometric",
    "RAMA": "torsion", "CONTACT": "statistical", "DISTPOT": "statistical", "ENV": "statistical",
    "HP": "statistical", "DSSPHB": "physics", "ELEC": "physics", "CONS": "consistency",
    "DMAP_CONS": "consistency", "TORS_CONS": "consistency", "POOLGO": "consistency",
    "SS_MATCH": "compatibility", "CONTACT_LL": "distogram-rereading", "DIS_MEAN": "distogram-rereading",
    "DIS": "reference", "LEG": "reference", "AMB": "reference",
}
LEG_TERMS = ("steric", "contact", "hbond_local", "hbond_longrange", "coop_helix", "coop_sheet",
             "solvation", "electrostatic", "aromatic", "torsion", "compactness")


def legacy_terms(cx: Context) -> Dict[str, np.ndarray]:
    from s16 import energy_lib as EL
    comp = EL.legacy_components_of_windows(cx.seq, cx.PHI, cx.PSI)
    w = EL.legacy_weight_vector()
    out = {}
    for t, wt in zip(LEG_TERMS, w):
        out["LEG_" + t] = float(wt) * np.asarray(comp[t], float)
    return out


def zrank(x):
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(r.std(), 1e-12)


def asinh_std(x):
    """A monotone, BOUNDED-tail standardisation that keeps target-dependent gaps (unlike
    zrank) without the float64 collapse of the raw moment on 1e28 outliers."""
    x = np.asarray(x, float)
    med = np.median(x)
    mad = np.median(np.abs(x - med)) + 1e-12
    return np.arcsinh((x - med) / (1.4826 * mad))
