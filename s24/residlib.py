"""s24/residlib.py -- WORKSTREAM B shared machinery for the residual generator.

Everything the lane needs that is NOT a model: the fixed candidate set, the deployed readout, the
torsion substrate, the bias-cosine spec measurement (identical arithmetic to `s24/biasalign.py`),
the mode-collapse audit, and the geometric-validity audit.

BASIS.  Point cloud throughout.  Every arm in this lane is compared against `P0R` -- the SAME 75
retained windows rebuilt from their own torsions with no residual applied -- never against the
real-window incumbent, so no residual arm is ever charged for the rebuild operator.

    B-0, n=126, complete:  incumbent 3.0483   rebuild 3.0524   +0.0041 SE 0.0047 MDE 0.0131  NULL

WHAT IS EXACT BY CONSTRUCTION, and therefore not evidence.  Every structure this lane emits is
built by `core.geometry.build_backbone`, which fixes omega at OMEGA_TRANS and uses ideal bond
lengths and angles.  So cis fraction, omega deviation, bond-length deviation and bond-angle
deviation are 0 by construction and CA-CA spacing is 3.80 A by construction -- they are reported
because the BRIEF requires them with every generated pool, and they are labelled BY CONSTRUCTION so
they are never mistaken for a result.  The validity axes that carry information here are
RAMACHANDRAN and CLASHES, which the builder does NOT guarantee.
"""
from __future__ import annotations

import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I            # noqa: E402
from s12 import key_lib as KL              # noqa: E402

TOPM = 75


# --------------------------------------------------------------------------- angles
def wrap(a):
    """The single wrapping convention for this lane.  Radians, into (-pi, pi]."""
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def circ_mean(A, axis=0):
    """Genuine circular mean -- arctan2(mean sin, mean cos).  Same object as `s14.retprior`."""
    A = np.asarray(A, float)
    return np.arctan2(np.sin(A).mean(axis), np.cos(A).mean(axis))


def circ_R(A, axis=0):
    """Mean resultant length in [0,1].  1 = a delta, 0 = uniform on the circle."""
    A = np.asarray(A, float)
    return np.hypot(np.sin(A).mean(axis), np.cos(A).mean(axis))


# --------------------------------------------------------------------------- the fixed set
def retained(pdb, m=TOPM):
    """The lane's FIXED candidate set: the shipped top-m of the shipped K=500 pool.

    Selection happens ONCE, on the real windows, with the shipped functional, BEFORE any residual
    exists -- so no arm can move the answer by changing membership.  That is the `functional` fork.

    Returns dict: W (m,n,3) real coords, PHI/PSI (m,n), nat (n,3), nphi/npsi (n,) ORACLE, seq, fold.
    """
    u = I.load_univ(pdb)
    n = u["n"]
    idx = I.pool_idx(u)
    W = u["W"][idx]
    dg = I.distogram(pdb)
    i, j = I.pair_index(n)
    sc = I.shipped_score(dg, I.pair_dists(W, i, j))
    o = np.argsort(sc, kind="stable")[:m]
    nat = np.asarray(u["nat_ca"], float)
    nphi, npsi = KL.native_torsions()[pdb]                       # ORACLE, evaluation/design only
    out = dict(pdb=pdb, n=int(n), fold=int(u["fold"]), seq=u["seq"], nat=nat,
               W=np.asarray(W[o], float),
               PHI=wrap(u["PHI"][idx][o]), PSI=wrap(u["PSI"][idx][o]),
               nphi=wrap(nphi), npsi=wrap(npsi), dg=dg, ij=(i, j))
    del u
    return out


# --------------------------------------------------------------------------- the readout
def readout(members):
    """The deployed operator: superpose on the medoid, uniform mean.  Returns the point cloud."""
    members = np.asarray(members, float)
    P = I.pairwise_rmsd(members)
    return I.superpose_batch(members, members[I.medoid(P)]).mean(0)


def emit(phi, psi):
    """(m,n) torsions -> the emitted point cloud through the deployed readout."""
    return readout(I.build_ca(phi, psi))


# --------------------------------------------------------------------------- the spec numbers
def bias(C, nat):
    """Error vector of an emitted cloud in the NATIVE frame -- byte-for-byte `biasalign._bias`."""
    C = np.asarray(C, float); nat = np.asarray(nat, float)
    Pc = C - C.mean(0); Qc = nat - nat.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (C - C.mean(0)) @ R.T - (nat - nat.mean(0))


def cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


# --------------------------------------------------------------------------- mode-collapse audit
def collapse_audit(phi, psi, ca=None, dup_tol=0.10):
    """MANDATORY BEFORE ANY RMSD CLAIM from a generator (BRIEF SS4).

    `phi`,`psi` are (B, n) samples for ONE target; `ca` optional (B,n,3), built if absent.

      n_unique        distinct structures at `dup_tol` A pairwise Ca-RMSD
      dup_frac        1 - n_unique/B
      ess             effective sample size, exp(entropy of the cluster occupancy)
      tors_entropy    mean per-residue circular entropy proxy, 1 - R (0 = a delta, 1 = uniform)
      pair_rmsd       mean/median/min pairwise Ca-RMSD within the sample
      mode_occ_max    largest cluster's share -- 1.0 is total collapse
    """
    phi = np.asarray(phi, float); psi = np.asarray(psi, float)
    B = len(phi)
    if ca is None:
        ca = I.build_ca(phi, psi)
    P = I.pairwise_rmsd(ca)
    iu = np.triu_indices(B, 1)
    d = P[iu]
    # single-linkage clustering at dup_tol -> "distinct structures"
    lab = -np.ones(B, int); c = 0
    for a in range(B):
        if lab[a] >= 0:
            continue
        stack = [a]; lab[a] = c
        while stack:
            x = stack.pop()
            for y in np.where((P[x] <= dup_tol) & (lab < 0))[0]:
                lab[y] = c; stack.append(int(y))
        c += 1
    occ = np.bincount(lab, minlength=c) / B
    ent = -(occ[occ > 0] * np.log(occ[occ > 0])).sum()
    R = np.concatenate([circ_R(phi, 0), circ_R(psi, 0)])
    return dict(B=int(B), n_unique=int(c), dup_frac=float(1.0 - c / B),
                ess=float(np.exp(ent)), mode_occ_max=float(occ.max()),
                tors_entropy=float(1.0 - R.mean()),
                pair_rmsd_mean=float(d.mean()) if len(d) else 0.0,
                pair_rmsd_median=float(np.median(d)) if len(d) else 0.0,
                pair_rmsd_min=float(d.min()) if len(d) else 0.0)


# --------------------------------------------------------------------------- validity audit
# Ramachandran favoured/allowed boxes, general case, in degrees.  Coarse on purpose: this is a
# validity SCREEN, not a structure-validation package, and it is used identically across arms.
_FAV = (((-180, -30), (-80, 30)),          # alpha_R / bridge
        ((-180, -30), (60, 180)),          # beta / PPII
        ((30, 100), (-20, 90)))            # alpha_L


def _in_boxes(p, s, boxes):
    p = np.asarray(p, float); s = np.asarray(s, float)
    ok = np.zeros(p.shape, bool)
    for (a, b), (c, d) in boxes:
        ok |= (p >= a) & (p <= b) & (s >= c) & (s <= d)
    return ok


def validity_audit(phi, psi, clash_ca=4.0):
    """Geometric validity, reported with EVERY generated pool (BRIEF SS4).

    Everything marked BY CONSTRUCTION is fixed by `core.geometry.build_backbone` and is not
    evidence of anything; it is reported because the rule requires the axis to appear.
    """
    phi = np.rad2deg(wrap(phi)); psi = np.rad2deg(wrap(psi))
    fav = _in_boxes(phi[..., 1:], psi[..., 1:], _FAV)     # phi[0] is never read by the builder
    ca = I.build_ca(wrap(np.deg2rad(phi)), wrap(np.deg2rad(psi)))
    if ca.ndim == 2:
        ca = ca[None]
    n = ca.shape[1]
    i, j = np.triu_indices(n, 3)                          # |i-j|>=3 non-local Ca contacts
    D = np.linalg.norm(ca[:, i] - ca[:, j], axis=-1)
    seq = np.linalg.norm(np.diff(ca, axis=1), axis=-1)    # consecutive Ca-Ca
    return dict(rama_favoured=float(fav.mean()), rama_outlier=float(1.0 - fav.mean()),
                clash_frac=float((D < clash_ca).mean()),
                clash_per_struct=float((D < clash_ca).sum(1).mean()),
                ca_ca_mean=float(seq.mean()), ca_ca_sd=float(seq.std()),
                chain_break_frac=float((np.abs(seq - 3.80) > 0.5).mean()),
                omega_dev_deg=0.0, cis_frac=0.0,          # BY CONSTRUCTION (OMEGA_TRANS)
                bond_len_dev=0.0, bond_ang_dev=0.0)       # BY CONSTRUCTION (ideal builder)


# --------------------------------------------------------------------------- statistics
def stats(d, fold, seed=0):
    """Paired summary of a difference vector: mean, SE, MDE, iid and fold-clustered CI, W/L,
    worst degradation, and the drop-top concentration curve the project requires beside W/L."""
    d = np.asarray(d, float); n = len(d)
    se = d.std(ddof=1) / np.sqrt(n)
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(4000)])
    fold = np.asarray(fold); F = sorted(set(fold.tolist()))
    fs = np.array([np.concatenate([d[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
                   for _ in range(4000)])
    o = np.argsort(d)
    return dict(mean=float(d.mean()), median=float(np.median(d)), se=float(se), mde=float(2.8016 * se),
                eff_over_mde=float(abs(d.mean()) / (2.8016 * se)) if se > 0 else float("inf"),
                ci_iid=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                ci_fold=[float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5))],
                W=int((d < 0).sum()), L=int((d > 0).sum()), worst=float(d.max()),
                drop_top10=float(d[o[10:]].mean()), drop_top20=float(d[o[20:]].mean()))


def verdict(s):
    """BEATS / worse / ns, with the Type-M zone called out rather than swallowed."""
    lo, hi = s["ci_fold"]
    if s["eff_over_mde"] < 0.7:
        return "ns"
    tag = " TYPE-M" if s["eff_over_mde"] <= 1.3 else ""
    if hi < 0:
        return "BEATS" + tag
    if lo > 0:
        return "worse" + tag
    return "ns"
