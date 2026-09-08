"""SPRINT 14 / VQE workstream -- shared machinery.

Everything the CVaR audit, the encoding study, the signal-tunable decisive experiment and
the ansatz study share.

THE CENTRAL OBJECT is `Enum`: one of the nine fully enumerated targets
(`s13/results/qarch_enum_<PDB>.npz`, n=9 residues, k=4, 4^9 = 262,144 configurations, every
one carrying a TRUE CA-RMSD).  Because the whole space is tabulated, an "objective" here is
just a length-262,144 float array, an objective evaluation is an array lookup, and the
certified global optimum of every objective is known exactly.  That is what makes a
budget-matched optimiser comparison possible at all.

INDEXING.  The enumeration is big-endian base-k: configuration index
``j = sum_i s_i k^(n-1-i)``, verified against the stored ``snap_index`` / ``snap_states``
pair.  At k=4 that means ``j``'s 18-bit binary expansion IS the binary-encoded register with
residue 0 in the most significant pair -- the same MSB-first convention
`core.quantum.all_bitstrings` and `StatevectorCircuit` use, so a statevector index and a
configuration index are the SAME integer with no relabelling.

LIVE QUBITS.  `core/project.py` leaves ``phi[0]``, ``psi[n-1]`` and ``phi[n-1]`` inert for
the CA trace.  At n=9 residue 8 is therefore entirely invisible to CA-RMSD, so of the 18
nominal qubits only **16 are live for CA-RMSD** (residue 0 keeps both of its qubits because
``psi[0]`` is live).  `Enum.live_qubits` measures this rather than asserting it.

NO NATIVE INFORMATION enters any objective built here unless the constructor is named
`ORACLE_*` or the `signal` parameter is non-zero -- the signal-tunable family in
`blend_objective` is an ORACLE DIAGNOSTIC by construction and is labelled as such
everywhere it is used.  Its purpose is to answer "at what objective quality does VQE become
useful", which cannot be asked without a quality knob.
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

RESULTS = os.path.join(ROOT, "s14", "results")
CACHE = os.path.join(ROOT, "s14", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

ENUM_TARGETS = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")
ENUM_DIR = os.path.join(ROOT, "s13", "results")


# ============================================================== the enumerated space
class Enum:
    """One fully enumerated target: every objective is an array over 4^n configurations."""

    def __init__(self, pdb: str):
        self.pdb = pdb
        d = np.load(os.path.join(ENUM_DIR, f"qarch_enum_{pdb}.npz"))
        self.n = int(d["n"])
        self.k = int(d["k"])
        self.seq = str(d["seq"])
        self.fold = int(d["fold"])
        self.bits_per_res = int(np.log2(self.k))
        self.n_qubits = self.n * self.bits_per_res
        self.N = self.k ** self.n
        self.rmsd = np.asarray(d["rmsd"], np.float64)          # ORACLE, post-hoc only
        self.legacy = np.asarray(d["legacy"], np.float64)
        self.prior = np.asarray(d["prior"], np.float64)
        self.PHI = np.asarray(d["PHI"], np.float64)
        self.PSI = np.asarray(d["PSI"], np.float64)
        self.snap_index = int(d["snap_index"])
        self.snap_states = np.asarray(d["snap_states"], int)
        self.leg = {c[4:]: np.asarray(d[c], np.float64)
                    for c in d.files if c.startswith("leg_")}
        self.amber_idx = np.asarray(d["amber_idx"], np.int64)
        self.amber_total = np.asarray(d["amber_total"], np.float64)
        assert self.rmsd.size == self.N
        # index <-> state vector, big-endian base-k (verified against snap_index)
        self._pw = self.k ** np.arange(self.n - 1, -1, -1)
        assert int(self.snap_states @ self._pw) == self.snap_index

    # -- index algebra ------------------------------------------------------
    def states(self, idx) -> np.ndarray:
        """(B,) configuration indices -> (B, n) state vectors."""
        j = np.atleast_1d(np.asarray(idx, np.int64))
        return (j[:, None] // self._pw[None, :]) % self.k

    def index(self, S) -> np.ndarray:
        return np.asarray(S, np.int64) @ self._pw

    def bits(self, idx) -> np.ndarray:
        """(B,) indices -> (B, n_qubits) MSB-first binary register."""
        j = np.atleast_1d(np.asarray(idx, np.int64))
        sh = np.arange(self.n_qubits - 1, -1, -1)
        return ((j[:, None] >> sh[None, :]) & 1).astype(np.uint8)

    # -- the live-qubit audit ----------------------------------------------
    def live_qubits(self) -> dict:
        """Which qubits actually move the CA trace? MEASURED, not asserted.

        A qubit is live iff flipping it changes the CA trace for at least one
        configuration.  Compares against `Enum.rmsd`, which is built from the same CA
        trace, so this measures the deployed pipeline's real dead-qubit defect.
        """
        rng = np.random.default_rng(0)
        base = rng.integers(0, self.N, 4000)
        live = []
        for q in range(self.n_qubits):
            flipped = base ^ (1 << (self.n_qubits - 1 - q))
            live.append(bool(np.any(self.rmsd[base] != self.rmsd[flipped])))
        return {"n_qubits_nominal": self.n_qubits,
                "live": [q for q in range(self.n_qubits) if live[q]],
                "dead": [q for q in range(self.n_qubits) if not live[q]],
                "n_live": int(sum(live))}


def load_all(targets=ENUM_TARGETS):
    return {p: Enum(p) for p in targets}


# ====================================================== rank / conditioning utilities
def ranks(x: np.ndarray) -> np.ndarray:
    """Average-tie ranks, 0-based."""
    x = np.asarray(x, float)
    o = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    r[o] = np.arange(len(x), dtype=float)
    xs = x[o]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def uniformise(x: np.ndarray) -> np.ndarray:
    """Monotone rank-preserving conditioning: map to uniform [0, 1] by rank.

    THIS IS THE MANDATORY PRE-STEP for any spectral comparison of objectives.  The brief
    records that the Walsh spectrum of a RAW molecular energy is a delta-spike artefact --
    the top-10 of 4,096 configurations carry a median 99.6% of raw AMBER's variance and the
    weight spectrum of a constant-plus-spike is exactly Binomial(m, 1/2).  Winsorisation at
    the 99th percentile is NOT enough because the spike survives it.  A rank transform is
    the strongest monotone conditioning available: it preserves every ordering the objective
    expresses (so every ranking statistic is untouched) while forcing all objectives onto
    the identical marginal, so a spectral difference between two conditioned objectives
    cannot be a difference of scale or of tail shape.
    """
    return ranks(x) / (len(x) - 1.0)


def spearman(a, b) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        return float("nan")
    ra, rb = ranks(a), ranks(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


# ============================================ the signal-tunable objective family
def blend_objective(base: np.ndarray, truth: np.ndarray, signal: float) -> np.ndarray:
    """ORACLE DIAGNOSTIC. Interpolate a realistic objective toward the true RMSD.

    Both inputs are put on the identical uniform-rank marginal first, so ``signal`` is a
    pure mixing weight in rank space and carries NO scale or tail-shape confound: at
    ``signal=0`` the result is a monotone image of `base` (identical ranking, identical
    argmin, identical optimiser behaviour up to ties), at ``signal=1`` a monotone image of
    `truth`.  The realised rank correlation with `truth` is not assumed -- every caller
    measures it with `spearman` and reports the measured value, because the map from mixing
    weight to realised rho is strongly non-linear.

    This family exists to answer the question "at what OBJECTIVE QUALITY does VQE start to
    beat classical search", which is unaskable without a quality knob.  It reads the native
    and is never a headline.
    """
    s = float(signal)
    return (1.0 - s) * uniformise(base) + s * uniformise(truth)


def noisy_truth(truth: np.ndarray, sigma: float, rng) -> np.ndarray:
    """ORACLE DIAGNOSTIC. Additive Gaussian noise on the uniformised truth.

    The second, independent parameterisation of objective quality.  A blend and a
    noise-corruption reach the same rho by different routes -- the blend keeps `base`'s own
    landscape structure underneath, the noise destroys structure isotropically -- so
    agreement between them is evidence the finding is about objective QUALITY and not about
    one particular family's shape.
    """
    u = uniformise(truth)
    return u + float(sigma) * rng.standard_normal(u.shape)


# ================================================================ classical searches
# Every search below counts objective evaluations the same way: one lookup into `E` is one
# evaluation, and the budget is a hard cap enforced by the caller-visible counter.
class Counter:
    """Hard evaluation budget. Never silently exceeded."""

    def __init__(self, E: np.ndarray, budget: int):
        self.E = np.asarray(E, float)
        self.budget = int(budget)
        self.used = 0
        self.best_i = -1
        self.best_e = np.inf
        self.seen = []

    @property
    def left(self) -> int:
        return self.budget - self.used

    def __call__(self, idx) -> np.ndarray:
        j = np.atleast_1d(np.asarray(idx, np.int64))
        if self.used + j.size > self.budget:
            j = j[: max(0, self.budget - self.used)]
            if j.size == 0:
                return np.empty(0)
        self.used += j.size
        e = self.E[j]
        m = int(np.argmin(e))
        if e[m] < self.best_e:
            self.best_e = float(e[m])
            self.best_i = int(j[m])
        self.seen.append(j)
        return e

    def all_seen(self) -> np.ndarray:
        return np.concatenate(self.seen) if self.seen else np.empty(0, np.int64)


def search_random(E, n, k, budget, rng):
    c = Counter(E, budget)
    N = k ** n
    while c.left > 0:
        c(rng.integers(0, N, min(c.left, 8192)))
    return c


def search_greedy(E, n, k, budget, rng, restarts=None):
    """1-opt coordinate descent with random restarts, budget-capped."""
    c = Counter(E, budget)
    pw = k ** np.arange(n - 1, -1, -1)
    while c.left > 0:
        s = rng.integers(0, k, n)
        cur = int(s @ pw)
        e = c([cur])
        if e.size == 0:
            break
        cur_e = float(e[0])
        improved = True
        while improved and c.left > 0:
            improved = False
            for i in rng.permutation(n):
                if c.left <= 0:
                    break
                cand = np.array([cur + (v - s[i]) * pw[i] for v in range(k)
                                 if v != s[i]], np.int64)
                ev = c(cand)
                if ev.size == 0:
                    break
                m = int(np.argmin(ev))
                if ev[m] < cur_e:
                    cur_e = float(ev[m])
                    cur = int(cand[m])
                    s = (cur // pw) % k
                    improved = True
    return c


def search_anneal(E, n, k, budget, rng, t0=None, t1=None):
    """Simulated annealing with single-residue moves; geometric temperature schedule."""
    c = Counter(E, budget)
    pw = k ** np.arange(n - 1, -1, -1)
    Es = np.asarray(E, float)
    sc = float(np.std(Es[rng.integers(0, len(Es), min(len(Es), 20000))]))
    t0 = 2.0 * sc if t0 is None else t0
    t1 = 0.02 * sc if t1 is None else t1
    t0 = max(t0, 1e-12); t1 = max(t1, 1e-15)
    s = rng.integers(0, k, n)
    cur = int(s @ pw)
    e0 = c([cur])
    if e0.size == 0:
        return c
    cur_e = float(e0[0])
    total = c.budget
    while c.left > 0:
        T = t0 * (t1 / t0) ** (1.0 - c.left / total)
        i = int(rng.integers(0, n))
        v = int(rng.integers(0, k - 1))
        v = v + 1 if v >= s[i] else v
        cand = cur + (v - s[i]) * pw[i]
        ev = c([cand])
        if ev.size == 0:
            break
        de = float(ev[0]) - cur_e
        if de <= 0 or rng.random() < np.exp(-de / T):
            cur_e = float(ev[0]); cur = int(cand); s[i] = v
    return c


def search_ga(E, n, k, budget, rng, pop=48, elite=12):
    """Generational GA: uniform crossover of the elite plus per-residue mutation."""
    c = Counter(E, budget)
    pw = k ** np.arange(n - 1, -1, -1)
    P = rng.integers(0, k, (pop, n))
    while c.left > 0:
        ev = c(P @ pw)
        if ev.size < len(P):
            break
        keep = P[np.argsort(ev)[:elite]]
        a = keep[rng.integers(0, elite, pop)]
        b = keep[rng.integers(0, elite, pop)]
        m = rng.random((pop, n)) < 0.5
        P = np.where(m, a, b)
        mut = rng.random((pop, n)) < (1.5 / n)
        P = np.where(mut, rng.integers(0, k, (pop, n)), P)
    return c


CLASSICAL = {"random": search_random, "greedy": search_greedy,
             "anneal": search_anneal, "ga": search_ga}


# ====================================================================== reporting
def summarise(c: Counter, rmsd: np.ndarray, E: np.ndarray, exact_i: int,
              near_thresh=(2.0, 2.5, 3.0)) -> dict:
    """Both axes, never substituted for one another.

    `objective_gap` is E_found - E_exact (the optimisation axis).
    `structural_gap` is RMSD_returned - RMSD_best_in_space (the accuracy axis).
    """
    seen = c.all_seen()
    out = {
        "evals": int(c.used),
        "best_e": float(c.best_e),
        "objective_gap": float(c.best_e - E[exact_i]),
        "rmsd_returned": float(rmsd[c.best_i]),
        "structural_gap": float(rmsd[c.best_i] - rmsd.min()),
        "rmsd_of_exact_optimum": float(rmsd[exact_i]),
    }
    if seen.size:
        r = rmsd[seen]
        out["rmsd_mean_seen"] = float(r.mean())
        out["rmsd_best_seen"] = float(r.min())
        out["n_distinct"] = int(np.unique(seen).size)
        for t in near_thresh:
            out[f"frac_below_{t}"] = float((r < t).mean())
    return out


def write(name, obj, results=RESULTS):
    import json
    path = os.path.join(results, name if name.endswith(".json") else name + ".json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)
    return path


def free_gb() -> float:
    from s12 import instrument as I
    return I.free_gb()


def wait_for_memory(min_gb=1.5, tag=""):
    import time
    while True:
        g = free_gb()
        if g >= min_gb:
            return g
        print(f"  [{tag}] free {g:.2f} GB < {min_gb}; waiting", flush=True)
        time.sleep(20)
