"""SPRINT 16 / QPHASE -- machinery for the OBJECTIVE-QUALITY PHASE BOUNDARY.

THE QUESTION.  Sprint 15 established that every apparent CVaR-VQE effect on this problem
dissolves under a classical Boltzmann reweighting of the SAME untrained circuit's samples,
and that annealing reaches the certified global optimum in 100% of cells at 1/400th of the
VQE's objective budget.  It also found the mechanism: the per-target VQE win correlates
+0.52 to +0.69 with WHERE THE OBJECTIVE'S OWN CERTIFIED ARGMIN SITS in the true RMSD
distribution.  Concentration pays exactly when the objective's optimum is the answer.

This module turns that observation into a measured phase boundary:

    at what objective quality rho does CVaR-VQE begin to beat each of its classical
    controls, and where does the real peptide-folding objective sit relative to that?

THE INSTRUMENT is the FULL enumerated register of 19 targets -- nine at n=9 (4^9 = 262,144
configurations, `s13/results/qarch_enum_*.npz`) and ten at n=10 (4^10 = 1,048,576,
`s14/cache/obj_enum_*.npz`).  Every configuration carries a true CA-RMSD, so the certified
global optimum of every objective is known exactly at every rung of the quality ladder.
The full register is used rather than Sprint 15's 12-qubit sub-register because on a
4,096-configuration space a budget large enough to train a VQE already exhausts the space,
which makes every argmin comparison degenerate.

TWO INDEPENDENT QUALITY KNOBS, both ORACLE DIAGNOSTICS by construction:
  * `V.blend_objective(base, truth, signal)` -- rank-space interpolation of a REAL objective
    toward the true RMSD.  Keeps the real objective's landscape structure underneath.
  * `V.noisy_truth(truth, sigma, rng)` -- additive Gaussian noise on the uniformised truth.
    Destroys structure isotropically instead.
Agreement between the two is the evidence that a boundary is about objective QUALITY and not
about one family's shape.  The realised Spearman rho is MEASURED at every rung and is the
x-axis of every result; the mixing weight is never used as a quality axis.

TIE DISCIPLINE.  The brief records that `np.argmin` on a tied signal reads the ORACLE sort
order and invented a 1.386 A winner.  Every argmin readout here averages the outcome over
the full tied argmin set (`argmin_readout`), and the tie multiplicity is returned.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I                      # noqa: E402
from s14 import vqe_lib as V                         # noqa: E402
from s14 import vqe_run as R                         # noqa: E402
from s14 import obj_enum as OE                       # noqa: E402
from s15 import seed as SD                           # noqa: E402

RESULTS = os.path.join(ROOT, "s16", "results")
CACHE = os.path.join(ROOT, "s16", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

TARGETS9 = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")
TARGETS10 = ("1N9U", "1TOR", "2BAO", "2MD2", "2MJQ", "5V5B", "6B9K", "7T3H", "7VI4", "8HVS")
TARGETS19 = TARGETS9 + TARGETS10

#: the sprint's primary budget convention, inherited from `s15/qrestraint.py`
BUDGET = 8192
SHOTS = 512
ANSATZ = "mps2f"
ALPHA = 0.25
SALT = "s16qphase"


# ===================================================================== resource gate
def cpu_pct():
    try:
        r = subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "(Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage "
             "-Average).Average"], capture_output=True, text=True, timeout=60)
        return float(r.stdout.strip())
    except Exception:
        return float("nan")


def wait_for_cpu(max_pct=97.0, tag="", poll=20, max_wait=180):
    t0 = time.time()
    while True:
        c = cpu_pct()
        if not np.isfinite(c) or c < max_pct:
            return c
        if time.time() - t0 > max_wait:
            print(f"  [{tag}] CPU still {c:.0f}% after {max_wait}s; PROCEEDING ON A "
                  f"CONTENDED BOX and recording it", flush=True)
            return c
        print(f"  [{tag}] CPU {c:.0f}% >= {max_pct}; waiting", flush=True)
        time.sleep(poll)


# ======================================================================= checkpoint
_CK = {}


def ck(tag, key, value):
    """Merge-on-write incremental checkpoint. Same contract as `s15.qens_lib.ck`."""
    path = os.path.join(RESULTS, f"qphase_{tag}.json")
    if os.path.exists(path):
        try:
            with open(path) as fh:
                disk = json.load(fh)
            mem = _CK.setdefault(tag, {})
            for k, v in disk.items():
                if isinstance(v, dict) and isinstance(mem.get(k), dict):
                    for kk, vv in v.items():
                        mem[k].setdefault(kk, vv)
                else:
                    mem.setdefault(k, v)
        except Exception:
            pass
    d = _CK.setdefault(tag, {})
    d[key] = value
    d["_written"] = time.strftime("%Y-%m-%d %H:%M:%S")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.4)
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def ck_load(tag):
    path = os.path.join(RESULTS, f"qphase_{tag}.json")
    if os.path.exists(path):
        with open(path) as fh:
            d = json.load(fh)
        _CK[tag] = d
        return d
    return {}


# ================================================================== the instrument
class Inst:
    """One fully enumerated target on its FULL register.

    Duck-types `s14.vqe_lib.Enum` and widens the file lookup to both enumeration
    directories via `s14.obj_enum.enum_path`, exactly as `s15.qrestraint.Enum` does.
    `rmsd`, `nat` and everything derived from them are ORACLE quantities used only for
    post-hoc scoring and for the deliberately-labelled quality knobs.
    """

    def __init__(self, pdb):
        d = np.load(OE.enum_path(pdb))
        self.pdb = pdb
        self.n = int(d["n"]); self.k = int(d["k"])
        self.seq = str(d["seq"]); self.fold = int(d["fold"])
        self.bits_per_res = int(np.log2(self.k))
        self.n_qubits = self.n * self.bits_per_res
        self.N = self.k ** self.n
        self.rmsd = np.asarray(d["rmsd"], np.float64)          # ORACLE
        self.legacy = np.asarray(d["legacy"], np.float64)
        self.prior = np.asarray(d["prior"], np.float64)
        self.PHI = np.asarray(d["PHI"], np.float64)
        self.PSI = np.asarray(d["PSI"], np.float64)
        self.snap_index = int(d["snap_index"])
        self.snap_states = np.asarray(d["snap_states"], int)
        self.amber_idx = np.asarray(d["amber_idx"], np.int64)
        self.amber_total = np.asarray(d["amber_total"], np.float64)
        self.amber_kind = np.asarray(d["amber_kind"], np.int64)
        self._pw = self.k ** np.arange(self.n - 1, -1, -1)
        assert self.rmsd.size == self.N
        assert int(self.snap_states @ self._pw) == self.snap_index
        self.nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)   # ORACLE
        self._u_truth = None
        self._hamil = None

    # ---------------------------------------------------------- index algebra
    def states(self, idx):
        j = np.atleast_1d(np.asarray(idx, np.int64))
        return (j[:, None] // self._pw[None, :]) % self.k

    def ca(self, idx):
        """(B, n, 3) CA traces, built with the same `I.build_ca` the labels used."""
        S = self.states(idx)
        ar = np.arange(self.n)
        return np.asarray(I.build_ca(self.PHI[ar[None, :], S],
                                     self.PSI[ar[None, :], S]), float)

    # ---------------------------------------------------------- the objectives
    def hamil(self, w=0.25):
        """The deployed native-free structural Hamiltonian (retrieval prior + distogram).

        Literally `s14.vqe_hamil.combine(*tabulate(pdb), w)` -- the objective the whole VQE
        line of this project optimises.  Cached full tabulations exist for all nineteen
        targets (`s14/cache/hamil_full_*.npz` and `s15/cache/qens_hamil_full_*.npz`), so
        this is a cache read.
        """
        if self._hamil is None:
            from s15.qens_lib import tabulate_any
            from s14.vqe_hamil import combine
            pr, ds = tabulate_any(self.pdb)
            self._hamil = combine(pr, ds, w).astype(np.float64)
        return self._hamil

    @property
    def u_truth(self):
        """ORACLE. `uniformise(rmsd)` -- cached to disk, it is a 1e6-element rank sort."""
        if self._u_truth is None:
            p = os.path.join(CACHE, f"utruth_{self.pdb}.npy")
            if os.path.exists(p):
                self._u_truth = np.load(p)
            else:
                self._u_truth = V.uniformise(self.rmsd)
                np.save(p, self._u_truth)
        return self._u_truth

    def restraint_objectives(self, batch=16384):
        """The two REAL objectives that consume the SAME distogram in opposite ways.

        `E_ls_pred` -- GENERATIVE: the weighted least-squares residual against the predicted
        distances, `sum_{|i-j|>=2} w_ij (d_ij(x) - dhat_ij)^2`, `w = 1/sd^2`.  NATIVE-FREE.
        `S14_disto_bayes` -- SELECTIVE: `I.shipped_score`, the distogram's Bayes-RISK score,
        the shipped pipeline's own consumption of the identical information.  NATIVE-FREE.

        Carrying both matters because Sprint 15 measured that they have very different
        in-tail ordering skill (in-tail rho +0.076 vs -0.075) and yielded the SAME VQE gain,
        which is the sharpest available test of whether a boundary drawn on rho predicts
        anything.  Streamed in batches so no (N, n_pairs) array is ever materialised.
        """
        p = os.path.join(CACHE, f"restraint_{self.pdb}.npz")
        if os.path.exists(p):
            d = np.load(p)
            return {k: np.asarray(d[k], np.float64) for k in d.files}
        i, j = I.pair_index(self.n)
        dg = I.distogram(self.pdb, self.seq, self.fold)
        assert np.array_equal(i, np.asarray(dg["i"])) and \
               np.array_equal(j, np.asarray(dg["j"])), "distogram pair order differs"
        dhat = np.asarray(dg["expected"], np.float32)
        w = (1.0 / np.maximum(np.asarray(dg["sd"], float), 1e-3) ** 2).astype(np.float32)
        els = np.empty(self.N, np.float64)
        bay = np.empty(self.N, np.float64)
        ar = np.arange(self.n)
        for s in range(0, self.N, batch):
            e = min(s + batch, self.N)
            S = self.states(np.arange(s, e))
            CA = np.asarray(I.build_ca(self.PHI[ar[None, :], S],
                                       self.PSI[ar[None, :], S]), np.float32)
            D = np.sqrt(((CA[:, i] - CA[:, j]) ** 2).sum(-1))
            r = D - dhat[None, :]
            els[s:e] = (r * r * w[None, :]).sum(1)
            bay[s:e] = I.shipped_score(dg, np.asarray(D, float))
        out = {"E_ls_pred": els, "S14_disto_bayes": bay}
        np.savez_compressed(p, **{k: v.astype(np.float32) for k, v in out.items()})
        return out

    def blend(self, base_u, signal):
        """ORACLE DIAGNOSTIC. `V.blend_objective` with the uniformisation hoisted.

        `blend_objective(base, truth, s)` is `(1-s)*uniformise(base) + s*uniformise(truth)`;
        both uniformisations are target-constant, so they are computed once and reused.
        Verified bit-identical to the library call in `qphase.py verify`.
        """
        return (1.0 - float(signal)) * base_u + float(signal) * self.u_truth

    def noisy(self, sigma, rng):
        """ORACLE DIAGNOSTIC. `V.noisy_truth(rmsd, sigma, rng)` with `u_truth` hoisted."""
        return self.u_truth + float(sigma) * rng.standard_normal(self.u_truth.shape)


_INST = {}


def inst(pdb):
    if pdb not in _INST:
        if len(_INST) > 2:                       # keep RAM bounded: 1e6 * ~9 arrays
            _INST.pop(next(iter(_INST)))
        _INST[pdb] = Inst(pdb)
    return _INST[pdb]


# ===================================================================== rank statistics
def ordinal_pct(x):
    """Tie-arbitrary percentile by argsort. DIAGNOSTIC only; claims use `V.ranks`."""
    o = np.argsort(np.asarray(x, float), kind="mergesort")
    r = np.empty(len(o), np.float64)
    r[o] = np.arange(len(o), dtype=np.float64)
    return r / max(1.0, len(o) - 1.0)


def rho_profile(E, rmsd, tails=(0.01, 0.001), n_null=40, seed=0):
    """The objective's realised quality: global rho, in-tail rho, and where its argmin sits.

    `rho_tail_*` is read against a MATCHED RANDOM-TAIL NULL (the sd of the identical
    statistic over tails of the same size drawn at random), never against zero -- the
    standing rule, and the reason Sprint 15's in-tail numbers were interpretable.
    `argmin_pct` is the mechanism variable Sprint 15 identified: the true-RMSD percentile of
    the objective's own certified global minimiser, averaged over the tied argmin set.
    """
    E = np.asarray(E, float); rmsd = np.asarray(rmsd, float)
    rng = np.random.default_rng(seed)
    out = {"rho_global": float(V.spearman(E, rmsd)), "n": int(E.size)}
    order = np.argsort(E, kind="mergesort")
    for t in tails:
        m = max(64, int(t * E.size))
        sub = order[:m]
        nulls = np.asarray([V.spearman(E[k], rmsd[k])
                            for k in rng.integers(0, E.size, (n_null, m))])
        out[f"rho_tail_{t}"] = float(V.spearman(E[sub], rmsd[sub]))
        out[f"rho_tail_{t}_null_sd"] = float(np.nanstd(nulls))
        out[f"rho_tail_{t}_null_mean"] = float(np.nanmean(nulls))
        out[f"tail_{t}_mean_rmsd"] = float(rmsd[sub].mean())
        out[f"tail_{t}_best_rmsd"] = float(rmsd[sub].min())
        out[f"tail_{t}_gap"] = float(rmsd[sub].mean() - rmsd[sub].min())
    tied = np.flatnonzero(E == E.min())
    out["argmin_rmsd"] = float(rmsd[tied].mean())
    out["argmin_ties"] = int(tied.size)
    out["argmin_pct"] = float((rmsd < rmsd[tied].mean()).mean())
    out["rmsd_best_in_space"] = float(rmsd.min())
    out["rmsd_mean_in_space"] = float(rmsd.mean())
    return out


# ========================================================================== the arms
_TH0 = {}


def theta0(n_qubits, seed, ansatz=ANSATZ, scale=0.8):
    """The untrained parameters. A function of (n_qubits, seed) only -- so the control arm
    is literally the arm the VQE starts from, and is shared across every objective."""
    key = (int(n_qubits), int(seed), ansatz, float(scale))
    if key not in _TH0:
        an = R.make_ansatz(ansatz, n_qubits)
        _TH0[key] = (an, R.init_theta(an, np.random.default_rng(seed), scale, "random"))
    return _TH0[key]


_UNTRAINED = {}


def untrained_draws(n_qubits, budget, seed, ansatz=ANSATZ, scale=0.8, shots=SHOTS):
    """THE control's draws: `budget` bitstrings from the circuit at its untrained theta.

    Identical to `s15.qrestraint._untrained_indices` (same construction, same streams), so
    this workstream's control arm is the same object Sprint 15's was.  Depends on E not at
    all, so it is drawn once per (n_qubits, budget, seed) and reused across objectives.
    """
    key = (int(n_qubits), int(budget), int(seed), ansatz, float(scale), int(shots))
    if key in _UNTRAINED:
        return _UNTRAINED[key]
    an, th0 = theta0(n_qubits, seed, ansatz, scale)
    rng = np.random.default_rng(seed)
    got, left = [], int(budget)
    while left > 0:
        bits = an.sample(th0, min(shots, left), rng)
        got.append(R._bits_to_index(bits))
        left -= got[-1].size
    idx = np.concatenate(got)[:budget]
    if len(_UNTRAINED) > 32:
        _UNTRAINED.pop(next(iter(_UNTRAINED)))
    _UNTRAINED[key] = idx
    return idx


_P0 = {}


def untrained_probs(n_qubits, seed, ansatz=ANSATZ, scale=0.8):
    """The EXACT untrained distribution over the whole register. Cached to disk: it is a
    2^18 / 2^20 enumeration costing 1.7 s / 8.6 s and depends only on (n_qubits, seed)."""
    key = (int(n_qubits), int(seed), ansatz, float(scale))
    if key in _P0:
        return _P0[key]
    path = os.path.join(CACHE, f"p0_{ansatz}_{n_qubits}_{seed}_{scale}.npy")
    if os.path.exists(path):
        p = np.load(path)
    else:
        an, th0 = theta0(n_qubits, seed, ansatz, scale)
        p = np.asarray(an.probs(th0), float)
        p = np.maximum(p, 0.0); p = p / p.sum()
        np.save(path, p)
    if len(_P0) > 3:
        _P0.pop(next(iter(_P0)))
    _P0[key] = p
    return p


def _entropy_bits(p):
    p = np.asarray(p, float)
    m = p > 0
    return float(-(p[m] * np.log2(p[m])).sum())


def draw_entropy_bits(idx, N):
    c = np.bincount(np.asarray(idx, np.int64), minlength=N).astype(float)
    return _entropy_bits(c / c.sum())


def tilt_samples(idx, E, target_entropy_bits, rng, iters=90):
    """CLASSICAL CONTROL, ZERO extra objective evaluations.

    `p ~ p_empirical * exp(-E/T)` over the multiset the UNTRAINED CIRCUIT ALREADY DREW, with
    `T` solved so the reweighted draw entropy equals the VQE's.  This is the brief's
    "classical Boltzmann reweighting of the same circuit's samples at matched entropy" in its
    strictest form: it reads no energy the control arm did not already pay for.

    NOTE, and it is a structural fact rather than a result: because it resamples from a set,
    its ARGMIN is identical to the untrained control's argmin by construction.  The tilt can
    only differ on ENSEMBLE readouts.  That identity is asserted in `qphase.py verify`.
    """
    idx = np.asarray(idx, np.int64)
    u, inv = np.unique(idx, return_inverse=True)
    w0 = np.bincount(inv).astype(float)
    w0 /= w0.sum()
    e = np.asarray(E, float)[u]
    e = e - e.min()

    def H_of(T):
        z = -e / max(T, 1e-30)
        q = w0 * np.exp(z - z.max())
        s = q.sum()
        if not np.isfinite(s) or s <= 0:
            return 0.0, w0
        q = q / s
        return _entropy_bits(q), q

    hi_h, hi_q = H_of(1e4)
    if target_entropy_bits >= hi_h:
        q = hi_q
    else:
        a, b = 1e-9, 1e4
        for _ in range(iters):
            m = np.sqrt(a * b)
            if H_of(m)[0] < target_entropy_bits:
                a = m
            else:
                b = m
        q = H_of(np.sqrt(a * b))[1]
    return u[rng.choice(len(u), size=idx.size, p=q)]


def expected_distinct(p, B):
    """E[# distinct configurations in B i.i.d. draws] = sum_i 1 - (1 - p_i)^B, closed form."""
    p = np.asarray(p, float)
    return float(np.sum(-np.expm1(B * np.log1p(-np.minimum(p, 1.0 - 1e-16)))))


def tilt_exact(p0, E, target_distinct, budget, rng, iters=45):
    """CLASSICAL CONTROL, charged the FULL register in objective evaluations.

    `p ~ p0 * exp(-E/T)` over the WHOLE space, with `T` solved so the EXPECTED number of
    distinct configurations in `budget` draws matches the VQE's observed distinct count --
    the brief's equal-diversity convention, in closed form so no matching noise enters.

    COST, stated wherever this arm is quoted: building it reads every entry of `E`, i.e.
    `N` = 262,144 or 1,048,576 objective evaluations against the VQE's 8,192.  On this
    instrument it is a 32x-to-128x BUDGET arm and is reported as an upper bound on what a
    reweighting of the untrained circuit can buy, never as a matched-cost control.  (On
    Sprint 15's 4,096-configuration sub-register the same arm was 200x CHEAPER than the
    exact-gradient VQE; the direction of that comparison is a property of the register size
    and must not be carried across.)
    """
    e = np.asarray(E, float)
    e = e - e.min()

    def D_of(T):
        z = -e / max(T, 1e-30)
        q = np.asarray(p0, float) * np.exp(z - z.max())
        s = q.sum()
        if not np.isfinite(s) or s <= 0:
            return 1.0, None
        q = q / s
        return expected_distinct(q, budget), q

    d_hi, q_hi = D_of(1e4)
    if target_distinct >= d_hi:
        q = q_hi
    else:
        a, b = 1e-9, 1e4
        for _ in range(iters):
            m = np.sqrt(a * b)
            if D_of(m)[0] < target_distinct:
                a = m
            else:
                b = m
        q = D_of(np.sqrt(a * b))[1]
    return rng.choice(len(q), size=budget, p=q), q


def tilt_anneal(p0, E, target_distinct, budget, rng, iters=16, iters_hi=1e3):
    """CLASSICAL CONTROL that isolates WHAT THE VQE'S SCHEDULE CONTRIBUTES.

    A VQE does not sample one distribution: it samples a SEQUENCE, broad early and narrow
    late, so its 8,192 draws carry an implicit annealing schedule that a single-temperature
    Boltzmann tilt does not have.  This arm gives the classical reweighting the same
    schedule and nothing else: a geometric temperature ladder from `iters_hi` (essentially
    `p0` itself) down to the temperature at which the expected distinct count matches the
    VQE's, drawing `budget/iters` samples at each step -- an ANNEALED INDEPENDENCE SAMPLER
    over exactly the VQE's own untrained support.

    It is charged the full register in objective evaluations, like `tilt_exact`.  Its
    purpose is diagnostic: if it reproduces the VQE, the VQE's advantage over a fixed tilt
    is its schedule and not its variational structure.
    """
    e = np.asarray(E, float)
    e = e - e.min()
    p0 = np.asarray(p0, float)

    def q_of(T):
        z = -e / max(T, 1e-30)
        q = p0 * np.exp(z - z.max())
        s = q.sum()
        return q / s if np.isfinite(s) and s > 0 else p0

    per = max(1, budget // iters)

    def D_of(T):
        return expected_distinct(q_of(T), budget)

    if target_distinct >= D_of(iters_hi):
        t_lo = iters_hi
    else:
        a, b = 1e-9, iters_hi
        for _ in range(40):
            m = np.sqrt(a * b)
            if D_of(m) < target_distinct:
                a = m
            else:
                b = m
        t_lo = np.sqrt(a * b)
    got = []
    left = budget
    for t in range(iters):
        T = iters_hi * (t_lo / iters_hi) ** (t / max(1, iters - 1))
        m = min(per, left)
        if m <= 0:
            break
        got.append(rng.choice(len(e), size=m, p=q_of(T)))
        left -= m
    if left > 0:
        got.append(rng.choice(len(e), size=left, p=q_of(t_lo)))
    return np.concatenate(got), float(t_lo)


# ======================================================================== the readout
def argmin_readout(E, rmsd, seen):
    """The argmin readout with the TIE RULE the brief mandates.

    `np.argmin` on a tied objective reads the array order, which on an enumerated space is
    the ORACLE sort order of nothing in particular but is still a leak of index structure.
    Every quantity here averages over the FULL tied argmin set and returns its size.
    """
    seen = np.asarray(seen, np.int64)
    es = np.asarray(E, float)[seen]
    m = es.min()
    tied = np.unique(seen[es == m])
    return {"best_e": float(m),
            "rmsd_returned": float(np.asarray(rmsd, float)[tied].mean()),
            "argmin_ties": int(tied.size)}


def readout(ins, E, seen, rng, e_pct=None, ms=(5, 20)):
    """Every axis, reported separately and never substituted for one another.

    OBJECTIVE axis: `objective_gap` (best E found minus the certified optimum) and
    `best_e_pct` (its rank percentile in the objective's own distribution).
    STRUCTURAL argmin axis: `rmsd_returned`, tie-averaged.
    STRUCTURAL ensemble axis: `rand{m}_coordavg_rmsd`, `set_mean_rmsd`, `set_best_rmsd`.
    DIVERSITY: `n_distinct`, `draw_entropy_bits` -- reported beside every aggregate,
    because the set-mean law is out of domain on a collapsed set.
    """
    seen = np.asarray(seen, np.int64)
    E = np.asarray(E, float)
    rmsd = ins.rmsd
    out = argmin_readout(E, rmsd, seen)
    tied_opt = np.flatnonzero(E == E.min())
    out["objective_gap"] = float(out["best_e"] - E.min())
    if e_pct is not None:
        out["best_e_pct"] = float(e_pct[seen[np.argmin(E[seen])]])
    out["rmsd_of_certified_optimum"] = float(rmsd[tied_opt].mean())
    out["structural_gap"] = float(out["rmsd_returned"] - rmsd.min())
    out["set_mean_rmsd"] = float(rmsd[seen].mean())
    out["set_best_rmsd"] = float(rmsd[seen].min())
    u = np.unique(seen)
    out["n_distinct"] = int(u.size)
    out["diversity"] = float(u.size / seen.size)
    out["draw_entropy_bits"] = draw_entropy_bits(seen, ins.N)
    for m in ms:
        pick = rng.choice(seen, size=m, replace=seen.size < m)
        W = ins.ca(np.unique(pick))
        if len(W) == 1:
            out[f"rand{m}_coordavg_rmsd"] = float(I.kabsch_rmsd_batch(W, ins.nat)[0])
        else:
            A = I.coordinate_average(W)
            A = A[0] if isinstance(A, tuple) else A
            out[f"rand{m}_coordavg_rmsd"] = float(
                I.kabsch_rmsd_batch(np.asarray(A, float)[None], ins.nat)[0])
        out[f"rand{m}_mean_rmsd"] = float(rmsd[pick].mean())
    return out


# ==================================================================== the arm bundle
def run_arms(ins, E, seed, budget=BUDGET, alpha=ALPHA, shots=SHOTS, lr=0.15,
             baseline="const", ansatz=ANSATZ, arms=None, e_pct=None, exact_dist=False):
    """Every arm on one (target, objective, seed), all at the SAME hard budget except the
    two arms whose different cost is stated in their name and in every table row."""
    arms = arms or ("vqe", "untrained", "tilt_samples", "tilt_exact",
                    "anneal", "anneal_q", "greedy", "random")
    rng = SD.stable_rng(ins.pdb, seed, "readout", salt=SALT)
    out = {}
    vq = None
    if "vqe" in arms:
        vq = R.run(E, ins.n_qubits, alpha, budget, shots=shots, ansatz=ansatz, seed=seed,
                   rmsd=ins.rmsd, bits_per_res=ins.bits_per_res, exact_dist=exact_dist,
                   baseline=baseline, lr=lr, keep_seen=True)
        r = readout(ins, E, vq["seen"], rng, e_pct)
        r.update({"iters": int(vq["iters"]), "evals": int(vq["evals"]),
                  "grad_norm_first": (float(vq["grad_norms"][0]) if vq["grad_norms"]
                                      else float("nan")),
                  "grad_norm_mean": (float(np.mean(vq["grad_norms"]))
                                     if vq["grad_norms"] else float("nan")),
                  "param_displacement": float(vq["param_displacement"]),
                  "init_rmsd_best": float(vq["init_rmsd_best"]),
                  "cost_evals": int(vq["evals"])})
        if exact_dist:
            r.update({k: vq[k] for k in ("entropy_bits", "max_prob", "mode_rmsd",
                                         "mean_rmsd_under_p", "eff_support")
                      if k in vq})
        out["vqe"] = r
    ud = untrained_draws(ins.n_qubits, budget, seed)
    if "untrained" in arms:
        out["untrained"] = readout(ins, E, ud, rng, e_pct)
        out["untrained"]["cost_evals"] = int(budget)
    if "tilt_samples" in arms and vq is not None:
        h = out["vqe"]["draw_entropy_bits"]
        ti = tilt_samples(ud, E, h, SD.stable_rng(ins.pdb, seed, "tilt", salt=SALT))
        out["tilt_samples"] = readout(ins, E, ti, rng, e_pct)
        out["tilt_samples"].update({"cost_evals": int(budget), "target_entropy": h})
    if "tilt_exact" in arms and vq is not None:
        p0 = untrained_probs(ins.n_qubits, seed)
        d = out["vqe"]["n_distinct"]
        ti, q = tilt_exact(p0, E, d, budget,
                           SD.stable_rng(ins.pdb, seed, "tiltx", salt=SALT))
        out["tilt_exact"] = readout(ins, E, ti, rng, e_pct)
        out["tilt_exact"].update({"cost_evals": int(ins.N), "target_distinct": d})
    if "tilt_anneal" in arms and vq is not None:
        p0 = untrained_probs(ins.n_qubits, seed)
        d = out["vqe"]["n_distinct"]
        ti, t_lo = tilt_anneal(p0, E, d, budget,
                              SD.stable_rng(ins.pdb, seed, "tilta", salt=SALT))
        out["tilt_anneal"] = readout(ins, E, ti, rng, e_pct)
        out["tilt_anneal"].update({"cost_evals": int(ins.N), "target_distinct": d,
                                   "T_final": t_lo})
    if "anneal" in arms:
        c = V.search_anneal(E, ins.n, ins.k, budget,
                            SD.stable_rng(ins.pdb, seed, "anneal", salt=SALT))
        out["anneal"] = readout(ins, E, c.all_seen(), rng, e_pct)
        out["anneal"]["cost_evals"] = int(c.used)
    if "anneal_q" in arms:
        b = max(64, budget // 4)
        c = V.search_anneal(E, ins.n, ins.k, b,
                            SD.stable_rng(ins.pdb, seed, "annealq", salt=SALT))
        out["anneal_q"] = readout(ins, E, c.all_seen(), rng, e_pct)
        out["anneal_q"]["cost_evals"] = int(c.used)
    if "greedy" in arms:
        c = V.search_greedy(E, ins.n, ins.k, budget,
                            SD.stable_rng(ins.pdb, seed, "greedy", salt=SALT))
        out["greedy"] = readout(ins, E, c.all_seen(), rng, e_pct)
        out["greedy"]["cost_evals"] = int(c.used)
    if "random" in arms:
        c = V.search_random(E, ins.n, ins.k, budget,
                            SD.stable_rng(ins.pdb, seed, "random", salt=SALT))
        out["random"] = readout(ins, E, c.all_seen(), rng, e_pct)
        out["random"]["cost_evals"] = int(c.used)
    return out


# ========================================================================= statistics
def paired(a, b, folds=None, names=None, seed=0):
    """Paired fold-aware bootstrap. THE UNIT IS THE TARGET -- seeds are averaged within a
    target before this is called, because treating (target x seed) cells as independent
    overstated n threefold in Sprint 15 and turned a null into a significant result."""
    return I.paired(np.asarray(a, float), np.asarray(b, float),
                    folds=folds, names=names, seed=seed)


def nullconc(d, sims=4000, seed=0):
    from s15.qgeom_nullconc import null_calibrated
    return null_calibrated(np.asarray(d, float), sims=sims, seed=seed)


def verdict(d, label="", floor=0.08):
    """One PASS/FAIL verdict line: the bootstrap CI, the median beside the mean, the W/L,
    the false-positive floor and the null-calibrated concentration check together."""
    d = np.asarray(d, float)
    p = I.paired(d, np.zeros_like(d))
    sig = (p["ci95"][0] > 0) or (p["ci95"][1] < 0)
    out = {"label": label, "n": p["n"], "mean": p["mean_diff"],
           "median": p["median_diff"], "ci95": p["ci95"],
           "W": p["n_better"], "L": p["n_worse"],
           "sig": bool(sig), "below_floor": bool(abs(p["mean_diff"]) <= floor),
           "mean_over_sd": float(p["mean_diff"] / (d.std(ddof=1) + 1e-12))}
    if sig and p["n"] > 10:
        try:
            c = nullconc(d, seed=seed_of(label))
            out["conc"] = {k: c[k] for k in c if not isinstance(c[k], (list, dict))}
        except Exception as ex:
            out["conc_error"] = str(ex)
    return out


def seed_of(s):
    return SD.stable_seed(s, salt=SALT) % (2 ** 31)
