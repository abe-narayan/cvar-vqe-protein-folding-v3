"""s17/feat_lib.py -- shared machinery for the FEATURES workstream.

THE ONE QUESTION THIS LANE ASKS.  Every native-free signal the programme has closed in-band is
GEOMETRIC or ENERGETIC -- computed from a candidate's own coordinates or from a physical
potential.  Nobody has asked whether a LEARNED SEQUENCE REPRESENTATION supplies in-band
discrimination.  This module supplies the representation-derived per-candidate features and the
in-band instrument that prices them.

WHAT IS AVAILABLE, established by `feat_inventory.py` and recorded here so a reader does not
have to re-derive it:

    esm_small.npz            per-residue ESM-2-650M representations for 360 short sequences,
                             covering all 126 tuning targets.  SEQUENCE-ONLY.
    esm_pca.npz              the shipped 1280 -> 32 whitening projection.
    s12/cache/esm_con_targets.npz
                             ESM-2 contact-head probabilities (n, n) for all 126 targets.
                             STRUCTURE-SUPERVISED (on other proteins, never on the target).
    esm_cache.npz            1.5 GB, ~22k library sequences.  Keyed by FULL sequence, not by
                             the length-n windows the universes hold.
    s8/generate_univ/<pdb>.npz `S` (nw, n) int8 -- the SOURCE SEQUENCE of every candidate window.

    PER-CANDIDATE ESM EMBEDDINGS DO NOT EXIST.  There are 13k-27k windows per target and the
    cache is keyed by parent sequence with no window->parent map.  Recomputing the full
    universe is out of budget; the BAND is not (see `feat_window.py`).

UNITS.  Spearman, pairwise ordering accuracy (null 0.500) and Gaussian-copula rho are three
different axes (LEDGER L3).  `inband_stats` measures Spearman AND counts pairwise accuracy
DIRECTLY, so neither is inferred from the other.  Conversions come from
`s17/selection_theory.py`'s named functions.

SIGN CONVENTION.  Every feature is returned with LOWER = PREDICTED BETTER, fixed a priori by the
feature's meaning.  This matters: a battery scored at both signs is a hidden best-of-2V screen
and would need the min-of-N null of LEDGER L8.  Where a sign is genuinely ambiguous it is
chosen LEAVE-FOLD-OUT and said so.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
CACHE = os.path.join(HERE, "cache")
RESULTS = os.path.join(HERE, "results")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s17 import sel_lib as L               # noqa: E402
from s17 import selection_theory as ST     # noqa: E402

K = 500
BANDS = ("top25", "top75", "oracle1.5")

#: ideal-geometry reference conformations, radians -- the ZERO-INFORMATION bar (BRIEF s4)
ALPHA = (np.deg2rad(-57.0), np.deg2rad(-47.0))
BETA = (np.deg2rad(-139.0), np.deg2rad(135.0))


# --------------------------------------------------------------------------- ESM caches
_ESM = None
_CON = None


def esm_reps():
    """seq -> (n, 1280) final-layer ESM-2 representation, for the short-sequence hot cache."""
    global _ESM
    if _ESM is None:
        z = np.load(os.path.join(ROOT, "esm_small.npz"), allow_pickle=True)
        _ESM = {str(s): np.asarray(v[0], np.float32) for s, v in zip(z["seqs"], z["vals"])}
    return _ESM


def esm_contacts():
    """seq -> (n, n) ESM-2 contact-head probabilities for the 126 tuning targets."""
    global _CON
    if _CON is None:
        z = np.load(os.path.join(ROOT, "s12", "cache", "esm_con_targets.npz"), allow_pickle=True)
        _CON = {str(s): np.asarray(c, np.float32) for s, c in zip(z["seqs"], z["con"])}
    return _CON


_PCA = None


def esm_pca32(seq):
    """(n, 32) whitened PCA projection of the per-residue representation -- the SHIPPED block."""
    global _PCA
    if _PCA is None:
        z = np.load(os.path.join(ROOT, "esm_pca.npz"))
        _PCA = (z["mu"], z["W"], z["scale"])
    mu, W, sc = _PCA
    return ((esm_reps()[seq] - mu) @ W) / np.maximum(sc, 1e-6)


def onehot(seq):
    """(n, 20) one-hot -- the MATCHED CONTROL for every ESM arm."""
    X = np.zeros((len(seq), 20), np.float32)
    for r, a in enumerate(seq):
        X[r, I.ALPHABET.index(a)] = 1.0
    return X


# --------------------------------------------------------------------------- CA descriptors
def env_desc(W):
    """Local CA environment of every residue: [d(r,r+2), d(r,r+3), d(r,r+4)].

    Defined for r = 0 .. n-5.  Computed identically for a native trace and for a candidate, so
    the two are commensurable -- which is the whole point of using CA-only descriptors when the
    native is available only as a CA trace.

    W: (n,3) or (B,n,3).  Returns (n-4, 3) or (B, n-4, 3).
    """
    W = np.asarray(W, float)
    single = W.ndim == 2
    if single:
        W = W[None]
    n = W.shape[1]
    m = n - 4
    out = np.empty((len(W), m, 3))
    for q, off in enumerate((2, 3, 4)):
        out[:, :, q] = np.linalg.norm(W[:, off:off + m] - W[:, :m], axis=-1)
    return out[0] if single else out


def ss_ca(W):
    """3-class secondary structure from CA coordinates ONLY.  0 = H, 1 = E, 2 = C.

    ONE FIXED RULE, applied identically to the native trace and to every candidate, so that a
    "predicted vs realised SS" comparison is between two commensurable objects.  The rule is the
    standard CA-only discriminator: an alpha helix has d(r,r+3) ~ 5.0-5.5 A and d(r,r+4) ~ 6.2 A;
    an extended strand has d(r,r+3) > 9 A.  Residues without a full r+4 window inherit the last
    assignable label.  This is a PROXY for DSSP and is labelled as such everywhere.

    W: (n,3) or (B,n,3).  Returns (n,) or (B,n) int8.
    """
    W = np.asarray(W, float)
    single = W.ndim == 2
    if single:
        W = W[None]
    n = W.shape[1]
    E = env_desc(W)                                  # (B, n-4, 3)
    d3, d4 = E[:, :, 1], E[:, :, 2]
    lab = np.full((len(W), n - 4), 2, np.int8)
    lab[(d3 < 6.4) & (d4 < 7.4)] = 0
    lab[d3 > 8.8] = 1
    out = np.full((len(W), n), 2, np.int8)
    if lab.shape[1]:
        out[:, :n - 4] = lab
        out[:, n - 4:] = lab[:, -1:]          # the trailing 4 residues inherit the last label
    return out[0] if single else out


def soft_contact(D, cut=8.0, tau=1.0):
    """Differentiable contact indicator, so a candidate's realised contact map is a probability
    on the same scale as ESM's contact head rather than a hard 0/1."""
    return 1.0 / (1.0 + np.exp((np.asarray(D, float) - cut) / tau))


# --------------------------------------------------------------------------- statistics
def spearman(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 3 or np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return 0.0
    ra = np.argsort(np.argsort(a, kind="stable"), kind="stable").astype(float)
    rb = np.argsort(np.argsort(b, kind="stable"), kind="stable").astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def pair_accuracy(score, truth):
    """DIRECTLY COUNTED pairwise ordering accuracy, null 0.500.

    Not converted from Spearman -- LEDGER L3 is exactly the cost of confusing the two axes.
    Ties in the score contribute 0.5, which is what a coin does.
    """
    s = np.asarray(score, float)
    y = np.asarray(truth, float)
    n = len(s)
    if n < 2:
        return 0.5
    ds = s[:, None] - s[None, :]
    dy = y[:, None] - y[None, :]
    m = np.triu(np.ones((n, n), bool), 1) & (dy != 0)
    if not m.any():
        return 0.5
    conc = (np.sign(ds[m]) == np.sign(dy[m])).sum()
    ties = (ds[m] == 0).sum()
    return float((conc + 0.5 * ties) / m.sum())


def inband_stats(score, rr):
    """Everything one feature is worth inside one band, on all three axes."""
    return dict(rho=spearman(score, rr), acc=pair_accuracy(score, rr), sel=L.sel_of(score, rr))


def axes(rho_s):
    """Spearman -> (copula rho, implied pairwise accuracy).  Named functions, per LEDGER L3."""
    r = ST.rho_from_spearman(rho_s)
    return r, ST.accuracy_from_rho(r)


def minN_null(rr_band, V, pdb, tag, nrep=400):
    """ZERO-INFORMATION minimum over V picks from the same band (LEDGER L8 / SELECT s6).

    Mandatory beside any 'best of V features chosen with labels' statistic.
    """
    rng = SD.stable_rng(pdb, "s17feat", tag)
    y = np.asarray(rr_band, float)
    pick = rng.integers(0, len(y), size=(int(nrep), int(V)))
    return float(np.mean(y[pick].min(axis=1)))


def bands_of(sc, rr):
    """The three band masks.  top25/top75 are NATIVE-FREE (the shipped shortlist); oracle1.5
    needs the native to draw and is a DIAGNOSTIC band, never deployable."""
    rk = np.argsort(sc, kind="stable")
    k = len(rr)
    m25 = np.zeros(k, bool); m25[rk[:25]] = True
    m75 = np.zeros(k, bool); m75[rk[:75]] = True
    return {"top25": m25, "top75": m75, "oracle1.5": rr <= rr.min() + 1.5}


def check_constants(n=13, verbose=True):
    """ASSERT that the zero-information references are real conformations, not a units bug.

    `core.geometry.build_backbone` takes RADIANS and passing degrees builds a DIFFERENT
    constant conformation silently -- one that still scores Ramachandran-favoured, so nothing
    raises.  The PHYSICS workstream hit this while building this exact control.  The falsifier of
    the whole FEATURES lane is stated against "beat a constant ideal alpha-helix", so the helix
    has to be an alpha-helix.  Measured, not asserted by inspection:

        ideal alpha : d(i,i+3) ~ 5.2, d(i,i+4) ~ 6.4, rise ~1.46 A/residue
        ideal beta  : d(i,i+3) ~ 10.5, rise ~3.46 A/residue
        the DEGREES bug gives d(i,i+3) = 7.49, which is neither.
    """
    out = {}
    for nm, (ph, ps) in (("alpha", ALPHA), ("beta", BETA)):
        ca = I.build_ca(np.full(n, ph), np.full(n, ps))
        d3 = float(np.linalg.norm(ca[3:] - ca[:-3], axis=1).mean())
        d4 = float(np.linalg.norm(ca[4:] - ca[:-4], axis=1).mean())
        rise = float(np.linalg.norm(ca[-1] - ca[0]) / (n - 1))
        out[nm] = (d3, d4, rise)
        if verbose:
            print(f"  ideal {nm}: d(i,i+3) {d3:.2f}  d(i,i+4) {d4:.2f}  rise {rise:.2f} A/res")
    assert 4.8 < out["alpha"][0] < 5.8 and 5.9 < out["alpha"][1] < 6.9, \
        "the 'ideal alpha helix' is not an alpha helix -- check radians vs degrees"
    assert 1.30 < out["alpha"][2] < 1.65, "alpha rise per residue is wrong"
    assert out["beta"][0] > 9.5 and out["beta"][2] > 3.0, \
        "the 'ideal beta strand' is not extended -- check radians vs degrees"
    #: and the CA-only SS rule must label them correctly, or `ss_ca` is not measuring SS
    assert ss_ca(I.build_ca(np.full(n, ALPHA[0]), np.full(n, ALPHA[1])))[0] == 0
    assert ss_ca(I.build_ca(np.full(n, BETA[0]), np.full(n, BETA[1])))[0] == 1
    return out


def null_signals(p, W, rng):
    """The ZERO-INFORMATION bar every arm must beat (BRIEF s4).  Lower = predicted better."""
    n = p["n"]
    out = {}
    for nm, (ph, ps) in (("NULL_alpha", ALPHA), ("NULL_beta", BETA)):
        ref = I.build_ca(np.full(n, ph), np.full(n, ps))
        out[nm] = I.kabsch_rmsd_batch(W, ref)
    out["NULL_chance"] = rng.standard_normal(len(W))
    return out
