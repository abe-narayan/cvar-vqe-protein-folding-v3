"""SPRINT 13, QUANTUM-ARCHITECTURE -- shared machinery.

Everything the encoding / locality / validity / enumeration experiments share:

  * the discrete torsion state space  (`torsion_lib2.library_for`, target held out)
  * batched decode -> backbone -> Legacy energy  (`core.energy.components_batch`)
  * single-config AMBER ff14SB/GBn2 single points  (`core.amber.single_point`, 6 ms warm)
  * a leakage-safe 1-local torsion PRIOR over the k states, plus an ORACLE prior whose
    quality is an explicit parameter q (the predictor agent owns the real one).

NOTHING here reads the native except the functions whose name starts `ORACLE_` and the
evaluation helpers (`rmsd_batch`), which are post-hoc scoring only.
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
import peptide_db as pdb                   # noqa: E402
import representations as reps             # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402
from core import energy as et              # noqa: E402
from core import geometry as geo           # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# ------------------------------------------------------------------ state space
class Space:
    """The discrete torsion state space of one target at one k.

    `PHI`, `PSI` are (n, k) radian tables from the sequence-conditioned, leakage-safe
    library.  A configuration is an integer vector `s` in {0..k-1}^n.
    """

    def __init__(self, pdb_id: str, k: int, seq: str = None, n: int = None,
                 fold: int = None):
        self.pdb = pdb_id
        if seq is None:
            t = {x["pdb"]: x for x in I.targets()}[pdb_id]
            seq, n, fold = t["seq"], t["n"], t["fold"]
        self.seq, self.n, self.fold, self.k = seq, int(n), int(fold), int(k)
        tab = tl2.library_for(seq, k, seq)                 # target held out
        self.rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        self.PHI = np.ascontiguousarray(self.rep._phi, dtype=float)   # (n, k)
        self.PSI = np.ascontiguousarray(self.rep._psi, dtype=float)
        self.bits_per_res = int(np.log2(k))
        assert 2 ** self.bits_per_res == k, "k must be a power of two"
        self.n_qubits_binary = self.n * self.bits_per_res
        self._rows = np.arange(self.n)

    # -- oracle handles (post-hoc / diagnostics only) ------------------------
    @property
    def nat_ca(self):
        if not hasattr(self, "_nat"):
            self._nat = I.load_univ(self.pdb)["nat_ca"]
        return self._nat

    def ORACLE_native_torsions(self):
        p = pdb.by_pdb(self.pdb)
        return np.asarray(p.phi, float), np.asarray(p.psi, float)

    def ORACLE_snap(self):
        """Nearest library state per residue to the native torsions. Chain-blind."""
        phi0, psi0 = self.ORACLE_native_torsions()
        d = (np.abs(wrap(self.PHI - phi0[:, None]))
             + np.abs(wrap(self.PSI - psi0[:, None])))
        return np.argmin(d, axis=1)

    # -- decode --------------------------------------------------------------
    def angles(self, S):
        """(B, n) state indices -> phi, psi (B, n) radians."""
        S = np.atleast_2d(np.asarray(S, int))
        return self.PHI[self._rows[None, :], S], self.PSI[self._rows[None, :], S]

    def coords(self, S):
        phi, psi = self.angles(S)
        return geo.build_backbone_batch(phi, psi), phi, psi

    def ca(self, S):
        phi, psi = self.angles(S)
        return I.build_ca(phi, psi)

    def rmsd(self, S, chunk=20000):
        """ORACLE post-hoc CA-RMSD to the native trace."""
        S = np.atleast_2d(np.asarray(S, int))
        out = np.empty(len(S))
        for a in range(0, len(S), chunk):
            out[a:a + chunk] = I.kabsch_rmsd_batch(self.ca(S[a:a + chunk]), self.nat_ca)
        return out

    # -- random configurations ----------------------------------------------
    def uniform(self, B, rng):
        return rng.integers(0, self.k, size=(B, self.n))

    def sample_prior(self, P, B, rng):
        """B configurations drawn from a per-residue distribution P (n, k)."""
        cum = np.cumsum(P, axis=1)
        u = rng.random((B, self.n))
        return (u[:, :, None] > cum[None, :, :-1]).sum(2)


# ------------------------------------------------------------------ the prior
def empirical_prior(space: Space) -> np.ndarray:
    """(n, k) leakage-safe occupancy of each library state.

    Assign every observation in the residue's back-off pool (built from the database with
    the target and everything >0.6 identity to it removed) to its nearest library state,
    in the same 4-D (cos/sin phi, cos/sin psi) embedding the library was clustered in.
    Purely a property of the held-out database and the target sequence: no native.
    """
    entries = pdb.holdout(space.seq)
    lib = tl2.ContextLibrary(entries)
    cls = reps.residue_classes(space.seq, space.n)
    P = np.zeros((space.n, space.k))
    for i in range(space.n):
        pool = lib.pools.get((cls[i],)) or lib.pools[(reps.CLASS_GENERAL,)]
        Q = np.asarray(pool, float)
        X = np.column_stack([np.cos(Q[:, 0]), np.sin(Q[:, 0]),
                             np.cos(Q[:, 1]), np.sin(Q[:, 1])])
        C = np.column_stack([np.cos(space.PHI[i]), np.sin(space.PHI[i]),
                             np.cos(space.PSI[i]), np.sin(space.PSI[i])])
        lab = ((X[:, None, :] - C[None]) ** 2).sum(-1).argmin(1)
        cnt = np.bincount(lab, minlength=space.k).astype(float)
        P[i] = (cnt + 1.0) / (cnt.sum() + space.k)          # Laplace
    return P


def ORACLE_prior(space: Space, q: float) -> np.ndarray:
    """Predictor-quality-parameterised prior: mass q on the snapped native state.

    ORACLE DIAGNOSTIC.  Stands in for the predictor agent's per-residue distribution so
    that every result here is a function of a single, explicit quality knob q.  q = 1/k
    is uniform; q = 1 is a perfect predictor.
    """
    s = space.ORACLE_snap()
    q = min(float(q), 1.0 - 1e-9)          # keep -log P finite at q -> 1
    P = np.full((space.n, space.k), (1.0 - q) / (space.k - 1))
    P[np.arange(space.n), s] = q
    return P


def prior_energy(P: np.ndarray, S) -> np.ndarray:
    """1-local torsion-prior energy: -sum_i log P[i, s_i].  Exactly 1-local by construction."""
    S = np.atleast_2d(np.asarray(S, int))
    L = -np.log(P)
    return L[np.arange(S.shape[1])[None, :], S].sum(1)


# ------------------------------------------------------------------ Legacy energy
LEGACY_TERMS = tuple(et.TERM_NAMES)


def legacy_components(space: Space, S, chunk=4096):
    """Per-term Legacy energies for (B, n) state configurations. dict term -> (B,)."""
    S = np.atleast_2d(np.asarray(S, int))
    outs = {t: np.empty(len(S)) for t in LEGACY_TERMS}
    for a in range(0, len(S), chunk):
        c, phi, psi = space.coords(S[a:a + chunk])
        comp = et.components_batch(space.seq, c, phi, psi)
        for t in LEGACY_TERMS:
            outs[t][a:a + chunk] = comp[t]
    return outs


def legacy_total(comp, weights=None):
    w = et.DEFAULT_WEIGHTS if weights is None else weights
    return sum(float(w.get(t, 0.0)) * np.asarray(comp[t], float) for t in comp)


def legacy_energy(space: Space, S, chunk=4096, weights=None):
    return legacy_total(legacy_components(space, S, chunk), weights)


# ------------------------------------------------------------------ AMBER energy
_AMB = {}


def amber_energies(space: Space, S, components=False, progress=None):
    """Genuine ff14SB/GBn2 single points (no minimisation) for (B, n) configurations.

    ~6 ms per configuration warm.  Returns (B,) total in kcal/mol, or (dict, total).
    """
    from core import amber as am
    S = np.atleast_2d(np.asarray(S, int))
    tot = np.empty(len(S))
    comps = {t: np.empty(len(S)) for t in am.AMBER_TERMS} if components else None
    for b in range(len(S)):
        r = am.single_point(space.seq, space.rep, [int(x) for x in S[b]],
                            components=components, threads=1)
        tot[b] = float(r["energy"])
        if components:
            for t in am.AMBER_TERMS:
                comps[t][b] = float(r["components"][t])
        if progress and b % progress == 0:
            print(f"    amber {b}/{len(S)}", flush=True)
    return (comps, tot) if components else tot


# ------------------------------------------------------------------ housekeeping
def wait_for_memory(min_gb=1.5, tag=""):
    import time
    while True:
        g = I.free_gb()
        if g >= min_gb:
            return g
        print(f"  [{tag}] free {g:.2f} GB < {min_gb}; waiting", flush=True)
        time.sleep(20)


def write(name, obj):
    import json
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return path


def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3:
        return float("nan")
    ra = _rank(a); rb = _rank(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def _rank(x):
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    r[order] = np.arange(len(x), dtype=float)
    # average ties
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def percentile_of(value, population):
    """Fraction of the population strictly below `value` (0 = best)."""
    p = np.asarray(population, float)
    return float((p < value).mean())
