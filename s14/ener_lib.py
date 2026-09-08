"""SPRINT 14, ENER -- shared machinery for the Legacy-vs-AMBER controlled comparison.

Everything the energy-model experiments share:

  * a cached accessor over `s13/results/qarch_enum_<PDB>.npz` (nine n=9 targets, k=4,
    262,144 configurations each, exact CA-RMSD + Legacy 11 terms + prior, plus genuine
    AMBER on a 2,955-config STRATIFIED subset),
  * the stratum labels for that subset (this matters: it is NOT uniform),
  * pairwise CA-RMSD between configurations, without materialising 262k structures,
  * controlled decoy sets at known structural distance, cached and importable,
  * scale-free / robust normalisation constants,
  * the discrimination metrics: in-decile rho, pairwise decision accuracy, top-k
    recovery, native percentile, calibration.

NOTHING here reads the native except through `rmsd` (a stored ORACLE label used for
post-hoc scoring) and functions named `ORACLE_*`.

    from s14 import ener_lib as E
    z = E.enum("1CS9")            # cached EnumTarget
    z.rmsd, z.legacy, z.prior     # (262144,)
    z.amber_total, z.amber_idx, z.amber_kind    # (2955,) 0=uniform 1=prior 2=ORACLE band
    z.uniform_mask                # the ONLY unconditioned AMBER stratum
"""
from __future__ import annotations

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
from s13 import qarch_lib as Q             # noqa: E402

S13_RESULTS = os.path.join(ROOT, "s13", "results")
RESULTS = os.path.join(ROOT, "s14", "results")
CACHE = os.path.join(ROOT, "s14", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

ENUM_TARGETS = ["1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5"]
K = 4
LEG_TERMS = ("steric", "contact", "hbond_local", "hbond_longrange", "coop_helix",
             "coop_sheet", "solvation", "electrostatic", "aromatic", "torsion",
             "compactness")
AMB_TERMS = ("bond", "angle", "torsion", "nonbonded", "solvation")
SEED = 20260905

spearman = Q.spearman
percentile_of = Q.percentile_of
wait_for_memory = Q.wait_for_memory


# ------------------------------------------------------------------ enumerated data
class EnumTarget:
    """Cached view of one fully-enumerated target."""

    def __init__(self, pdb_id):
        self.pdb = pdb_id
        z = np.load(os.path.join(S13_RESULTS, f"qarch_enum_{pdb_id}.npz"))
        self._z = {k: z[k] for k in z.files}
        self.n = int(self._z["n"])
        self.k = int(self._z["k"])
        self.seq = str(self._z["seq"])
        self.fold = int(self._z["fold"])
        self.rmsd = np.asarray(self._z["rmsd"], float)          # ORACLE label
        self.legacy = np.asarray(self._z["legacy"], float)
        self.prior = np.asarray(self._z["prior"], float)
        self.snap_index = int(self._z["snap_index"])
        self.snap_states = np.asarray(self._z["snap_states"], int)
        self.PHI = np.asarray(self._z["PHI"], float)
        self.PSI = np.asarray(self._z["PSI"], float)
        self.prior_table = np.asarray(self._z["prior_table"], float)
        self.amber_idx = np.asarray(self._z["amber_idx"], int)
        self.amber_kind = np.asarray(self._z["amber_kind"], int)
        self.amber_total = np.asarray(self._z["amber_total"], float)
        self.leg = {t: np.asarray(self._z["leg_" + t], float) for t in LEG_TERMS}
        self.amb = {t: np.asarray(self._z["amb_" + t], float) for t in AMB_TERMS}
        self.B = len(self.rmsd)

    # -- AMBER strata --------------------------------------------------------
    @property
    def uniform_mask(self):
        """The only stratum of the AMBER subset that is an unconditioned sample."""
        return self.amber_kind == 0

    @property
    def prior_mask(self):
        return self.amber_kind == 1

    @property
    def band_mask(self):
        """ORACLE near-native band: selected by TRUE RMSD. Conditioned on the label."""
        return self.amber_kind == 2

    # -- decode --------------------------------------------------------------
    def states(self, idx):
        """Config indices -> (m, n) state matrix (odometer order, base k)."""
        idx = np.atleast_1d(np.asarray(idx, np.int64))
        out = np.empty((len(idx), self.n), np.int8)
        rem = idx.copy()
        for j in range(self.n - 1, -1, -1):
            out[:, j] = rem % self.k
            rem //= self.k
        return out

    def index_of(self, S):
        S = np.atleast_2d(np.asarray(S, np.int64))
        w = self.k ** np.arange(self.n - 1, -1, -1, dtype=np.int64)
        return (S * w[None, :]).sum(1)

    def ca(self, idx, chunk=20000):
        """Config indices -> (m, n, 3) CA traces."""
        idx = np.atleast_1d(np.asarray(idx, np.int64))
        out = np.empty((len(idx), self.n, 3))
        rows = np.arange(self.n)
        for a in range(0, len(idx), chunk):
            S = self.states(idx[a:a + chunk]).astype(int)
            phi = self.PHI[rows[None, :], S]
            psi = self.PSI[rows[None, :], S]
            out[a:a + chunk] = I.build_ca(phi, psi)
        return out

    def space(self):
        if not hasattr(self, "_space"):
            self._space = Q.Space(self.pdb, self.k, seq=self.seq, n=self.n, fold=self.fold)
        return self._space

    # -- objective assembly --------------------------------------------------
    def obj(self, name):
        """Named objective over the FULL enumeration (AMBER excluded -- subset only)."""
        if name == "legacy":
            return self.legacy
        if name == "prior":
            return self.prior
        if name.startswith("leg_"):
            return self.leg[name[4:]]
        if name == "legacy_nosteric":
            from core import energy as et
            w = et.DEFAULT_WEIGHTS
            return sum(float(w.get(t, 0.0)) * self.leg[t]
                       for t in LEG_TERMS if t != "steric")
        raise KeyError(name)

    def amber_obj(self, name):
        """Named objective over the AMBER SUBSET, aligned to `amber_idx`."""
        if name in ("amber", "amber_total"):
            return self.amber_total
        if name.startswith("amb_"):
            return self.amb[name[4:]]
        if name == "amber_nonclash":
            # AMBER minus the term that explodes.  NOTE the cached decomposition lumps
            # vdW and Coulomb into `nonbonded`, so electrostatics cannot be separated
            # here; this is bond+angle+torsion+solvation.
            return (self.amb["bond"] + self.amb["angle"] + self.amb["torsion"]
                    + self.amb["solvation"])
        if name == "amber_soft":
            return softcore(self.amber_total)
        # fall back to a full-enumeration objective, restricted
        return self.obj(name)[self.amber_idx]


_ENUM = {}


def enum(pdb_id) -> EnumTarget:
    if pdb_id not in _ENUM:
        _ENUM[pdb_id] = EnumTarget(pdb_id)
    return _ENUM[pdb_id]


# ------------------------------------------------------------------ transforms
def softcore(x):
    """Sign-preserving monotone log compression. Changes SCALE, never ORDER."""
    x = np.asarray(x, float)
    return np.sign(x) * np.log1p(np.abs(x))


def rank_norm(x):
    """Rank -> [0, 1]. Scale-free, monotone, finite on any distribution."""
    x = np.asarray(x, float)
    r = Q._rank(x)
    return r / max(len(x) - 1, 1)


def robust_z(x):
    """Median / MAD standardisation -- finite even when the mean and sd are not."""
    x = np.asarray(x, float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    s = 1.4826 * mad
    if not np.isfinite(s) or s <= 0:
        s = 1.0
    return (x - med) / s


# ------------------------------------------------------------------ metrics
def decile_rho(E, R, frac=0.10):
    """Rank correlation with RMSD INSIDE the lowest-energy fraction of the population.

    This is where a search actually lives; the global rho is not the operative number.
    """
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    m = max(int(round(frac * len(E))), 10)
    if m > len(E):
        m = len(E)
    sel = np.argsort(E, kind="mergesort")[:m]
    return spearman(E[sel], R[sel])


def topk_recovery(E, R, k=10, good_frac=0.01):
    """Fraction of the energy top-k that lies in the truly-best `good_frac` of the pool."""
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    ng = max(int(round(good_frac * len(R))), 1)
    good = set(np.argsort(R, kind="mergesort")[:ng].tolist())
    sel = np.argsort(E, kind="mergesort")[:k]
    return float(np.mean([int(i in good) for i in sel]))


def pair_accuracy(E, R, n_pairs=200000, rng=None, min_gap=0.0):
    """P(sign(E_i - E_j) == sign(R_i - R_j)) over random pairs. 0.5 = no skill."""
    rng = rng or np.random.default_rng(0)
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    m = len(E)
    i = rng.integers(0, m, n_pairs); j = rng.integers(0, m, n_pairs)
    dE = E[i] - E[j]; dR = R[i] - R[j]
    keep = (dE != 0) & (np.abs(dR) > min_gap)
    if keep.sum() == 0:
        return float("nan")
    return float(((dE > 0) == (dR > 0))[keep].mean())


def tail_accuracy(E, R, q=0.01, n_pairs=200000, rng=None, gap_band=None, min_tail=30):
    """Pairwise decision accuracy INSIDE the objective's own lowest-q fraction.

    This is the regime an argmin actually samples, which neither global `pair_accuracy`
    (dominated by the 99% of the space a search never visits) nor `decile_rho` at a fixed
    10% measures.  Requested by the coordinator; exported so other objectives compose.

    TWO CONFOUNDS, both reported rather than hidden:

      * the tail is OBJECTIVE-DEFINED, so different objectives are scored on different
        populations.  `tail_mean_gap` and `tail_mean_rmsd` say what each population is.
      * accuracy is a function of the QUALITY GAP (finding E2a), and gaps inside a tail are
        small, so EVERY objective's tail accuracy falls toward 0.5 as q shrinks whether or
        not it has lost skill.  `gap_band` restricts to pairs with |dRMSD| in a fixed
        window so the comparison is like-for-like across objectives and across q.
    """
    rng = rng or np.random.default_rng(0)
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    m = max(int(round(q * len(E))), 2)
    sel = np.argsort(E, kind="mergesort")[:m]
    e, r = E[sel], R[sel]
    i = rng.integers(0, m, n_pairs); j = rng.integers(0, m, n_pairs)
    dE = e[i] - e[j]; dR = r[i] - r[j]
    gap = np.abs(dR)
    keep = dE != 0
    # A tail of a dozen configurations can still supply 200,000 sampled pairs, which is
    # pseudo-replication: the accuracy is then an average over a handful of DISTINCT
    # comparisons and reads 0.000 or 1.000.  Require enough distinct members, not enough
    # pairs.
    enough = m >= min_tail
    out = dict(q=float(q), n_tail=int(m), enough=bool(enough),
               tail_mean_rmsd=float(r.mean()), tail_best_rmsd=float(r.min()),
               tail_mean_gap=float(gap.mean()),
               tie_frac=float(1.0 - keep.mean()),
               acc=float(((dE > 0) == (dR > 0))[keep].mean())
               if (enough and keep.sum() >= 50) else float("nan"))
    if gap_band is not None:
        lo, hi = gap_band
        k2 = keep & (gap >= lo) & (gap < hi)
        n_distinct = len(np.unique(np.concatenate([i[k2], j[k2]]))) if k2.any() else 0
        out["acc_gap_matched"] = (float(((dE > 0) == (dR > 0))[k2].mean())
                                  if (enough and k2.sum() >= 200 and n_distinct >= min_tail)
                                  else float("nan"))
        out["n_gap_matched"] = int(k2.sum())
        out["n_distinct_gap_matched"] = int(n_distinct)
    return out


def selection_decomposition(E, R, q=0.01):
    """Exactly split an objective's selected RMSD into FILTERING and ORDERING.

        sel_rmsd = pool_mean + (tail_mean - pool_mean) + (sel_rmsd - tail_mean)
                              \\____ filtering ____/   \\____ ordering ____/

    Filtering is the value of being in the tail at all; ordering is the value of the
    objective's ranking WITHIN its own tail.  This is what explains an objective with
    excellent bulk discrimination and a poor argmin: all of its value is filtering, and its
    ordering term is zero or positive.
    """
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    m = max(int(round(q * len(E))), 2)
    sel = np.argsort(E, kind="mergesort")[:m]
    tail_mean = float(R[sel].mean())
    s = argmin_rmsd(E, R)
    return dict(q=float(q), n_tail=int(m), pool_mean=float(R.mean()),
                tail_mean=tail_mean, sel_rmsd=s,
                filtering=tail_mean - float(R.mean()),
                ordering=s - tail_mean)


def argmin_rmsd(E, R):
    """RMSD of the objective's argmin, averaging over the tied argmin set.

    Tie-averaging is mandatory: np.argmin on a tied signal reads the array order and
    (project finding) has invented a winner before.
    """
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    lo = E.min()
    tie = np.flatnonzero(E == lo)
    return float(R[tie].mean())


def calibration(E, R):
    """OLS slope / intercept / r of R on a rank-normalised E, plus Pearson r."""
    e = rank_norm(E); r = np.asarray(R, float)
    ok = np.isfinite(e) & np.isfinite(r)
    e, r = e[ok], r[ok]
    if len(e) < 3 or e.std() == 0:
        return dict(slope=float("nan"), intercept=float("nan"), pearson=float("nan"))
    b = np.cov(e, r, ddof=1)[0, 1] / e.var(ddof=1)
    a = r.mean() - b * e.mean()
    p = np.corrcoef(e, r)[0, 1]
    return dict(slope=float(b), intercept=float(a), pearson=float(p))


def energy_rmsd_curve(E, R, nbin=20):
    """Mean RMSD in each energy quantile bin -- the full energy-vs-RMSD curve."""
    E = np.asarray(E, float); R = np.asarray(R, float)
    ok = np.isfinite(E) & np.isfinite(R)
    E, R = E[ok], R[ok]
    order = np.argsort(E, kind="mergesort")
    edges = np.linspace(0, len(order), nbin + 1).astype(int)
    return [dict(bin=b, n=int(edges[b + 1] - edges[b]),
                 e_mean=float(E[order[edges[b]:edges[b + 1]]].mean()),
                 rmsd_mean=float(R[order[edges[b]:edges[b + 1]]].mean()),
                 rmsd_min=float(R[order[edges[b]:edges[b + 1]]].min()))
            for b in range(nbin) if edges[b + 1] > edges[b]]


# ------------------------------------------------------------------ pairwise structure
def pairwise_rmsd(z: EnumTarget, idx_a, idx_b=None, chunk=2000):
    """CA-RMSD between configurations (not to the native). Returns (len(a), len(b))."""
    A = z.ca(idx_a)
    Bc = A if idx_b is None else z.ca(idx_b)
    out = np.empty((len(A), len(Bc)))
    for a in range(0, len(A), chunk):
        for b in range(len(Bc)):
            out[a:a + chunk, b] = I.kabsch_rmsd_batch(A[a:a + chunk], Bc[b])
    return out


# ------------------------------------------------------------------ housekeeping
def write(name, obj, n_expected=None):
    import json
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    if isinstance(obj, dict) and n_expected is not None:
        rows = obj.get("per_target") or obj.get("rows") or []
        obj = dict(obj, n_rows=len(rows), n_expected=int(n_expected),
                   complete=bool(len(rows) == int(n_expected)))
    import time
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    # Windows raises PermissionError if the destination is momentarily held open (an
    # editor, an indexer, another reader).  Retry rather than lose a long run's output.
    for attempt in range(10):
        try:
            os.replace(tmp, path)
            return path
        except PermissionError:
            if attempt == 9:
                raise
            time.sleep(0.5)
    return path


def bootstrap_mean_ci(x, n_boot=4000, seed=0):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 2:
        return float(x.mean()) if len(x) else float("nan"), [float("nan")] * 2
    rng = np.random.default_rng(seed)
    bs = np.array([x[rng.integers(0, len(x), len(x))].mean() for _ in range(n_boot)])
    return float(x.mean()), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
