"""SPRINT 15 / QENS -- machinery for the replication of the project's ONE positive VQE result.

WHAT IS UNDER TEST.  `s15/qgeom_ens.py` E2 measured that a CVaR-VQE, consumed as an ENSEMBLE
with NO ranker anywhere, beats best-of-N from its own untrained circuit by 0.36-0.57 A with
CIs excluding zero, at alpha=1 and alpha=0.05.  n = 27 paired cells (9 targets x 3 seeds).
The QGEOM agent flagged it LOW POWER (mean/sd = -0.40) and wrote "replicate before building
on it".  This module replicates it at power and decomposes it.

WHAT THIS ADDS OVER `qgeom_ens`.

1. **19 enumerated targets, not 9.**  `s13/results/qarch_enum_*.npz` (9 at n=9) plus
   `s14/cache/obj_enum_*.npz` (10 at n=10).  `GEnum` is a loader that duck-types
   `s14.vqe_lib.Enum` over both directories.  Every target keeps the same 12-qubit
   sub-register (residues 1..6), so the instance size is IDENTICAL across the 19 and the
   only thing that changes is which target.

2. **A precomputed CA table.**  The whole 4,096-configuration sub-register is built once per
   target as a (4096, n, 3) array, so a readout is an array index rather than 4,096 calls to
   `build_ca`.  Verified against the stored `Enum.rmsd` column to 1e-6 A.

3. **Arms the original did not have** (all at matched budget, all consumed with NO ranker):
      vqe            the CVaR-VQE final distribution
      untrained      the same circuit at its untrained theta -- the incumbent control
      tilt           a CLASSICAL Boltzmann tilt of p_init on the SAME objective, temperature
                     solved so its entropy equals the VQE's.  Native-free.  This is the
                     "did the circuit do anything a one-line classical reweighting cannot"
                     control, and it is the sharp test of whether the effect is the
                     OBJECTIVE's tilt rather than the variational optimisation.
      anneal         classical simulated annealing at EQUAL DRAWS (2,048 evaluations)
      anneal_cost    classical simulated annealing at EQUAL OBJECTIVE EVALUATIONS -- the
                     exact-gradient VQE reads all 4,096 objective values on every one of its
                     200 iterations, i.e. 819,200 evaluations, so convention B charges the
                     classical arm the same.
      uniform        uniform random draws.  A WEAK control: the brief records that a
                     zero-information constant alpha-helix beats uniform random by 0.457 A,
                     so beating it proves nothing and it is reported only for calibration.

4. **A matched-DIVERSITY readout.**  `mdN_coordavg` subsamples every arm's DISTINCT
   configuration set down to a common size `D` (the minimum over the arms in the cell) and
   averages `m` members drawn uniformly from that.  Equal diversity, equal budget: if the VQE
   still wins there the effect is about WHICH configurations it favours (location/shape); if
   it does not, the effect is diversity, which the untrained circuit supplies more cheaply.

NATIVE-FREE.  Every objective is the native-free structural Hamiltonian (retrieval torsion
prior + leave-fold-out distogram).  `Enum.rmsd`, `nat_ca` and every RMSD are ORACLE
quantities used only for post-hoc scoring, and are labelled ORACLE wherever they are printed.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

RESULTS = os.path.join(ROOT, "s15", "results")
CACHE = os.path.join(ROOT, "s15", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

from s12 import instrument as I                          # noqa: E402
from s14 import vqe_lib as V                             # noqa: E402
from s15 import qgeom_lib as G                           # noqa: E402
from s15 import seed as SD                               # noqa: E402

# the nine n=9 enumerations (s13) and the ten n=10 enumerations (s14 cache)
TARGETS9 = ("1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5")
TARGETS10 = ("1N9U", "1TOR", "2BAO", "2MD2", "2MJQ", "5V5B", "6B9K", "7T3H", "7VI4", "8HVS")
TARGETS19 = TARGETS9 + TARGETS10
SUB = (1, 2, 3, 4, 5, 6)                 # the 12-qubit sub-register, as in qgeom_ens


# ====================================================================== checkpoint
_CK = {}


def ck(tag, key, value):
    """Merge-on-write incremental checkpoint, same contract as `qgeom_lib.ck`."""
    path = os.path.join(RESULTS, f"qens_{tag}.json")
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
    path = os.path.join(RESULTS, f"qens_{tag}.json")
    if os.path.exists(path):
        with open(path) as fh:
            d = json.load(fh)
        _CK[tag] = d
        return d
    return {}


# ==================================================================== the enumeration
class GEnum:
    """One fully enumerated target from EITHER enumeration directory.

    Duck-types `s14.vqe_lib.Enum` for everything this workstream needs.  The n=10 caches use
    the identical schema (`s14/obj_enum.py` says so and the key list confirms it), so the only
    difference is the path and the dtype of the big columns.
    """

    def __init__(self, pdb):
        from s14 import obj_enum as OE
        d = np.load(OE.enum_path(pdb))
        self.pdb = pdb
        self.n = int(d["n"])
        self.k = int(d["k"])
        self.seq = str(d["seq"])
        self.fold = int(d["fold"])
        self.N = self.k ** self.n
        self.rmsd = np.asarray(d["rmsd"], np.float64)     # ORACLE, post-hoc scoring only
        self.PHI = np.asarray(d["PHI"], np.float64)
        self.PSI = np.asarray(d["PSI"], np.float64)
        self.snap_index = int(d["snap_index"])
        self.snap_states = np.asarray(d["snap_states"], int)
        self._pw = self.k ** np.arange(self.n - 1, -1, -1)
        assert self.rmsd.size == self.N
        assert int(self.snap_states @ self._pw) == self.snap_index

    def states(self, idx):
        j = np.atleast_1d(np.asarray(idx, np.int64))
        return (j[:, None] // self._pw[None, :]) % self.k

    def index(self, S):
        return np.asarray(S, np.int64) @ self._pw


class Struct:
    """A 12-qubit sub-register of one enumerated target, with every CA trace precomputed.

    Construction is bit-identical to `s15.qgeom_cvar.Struct`: the frozen residues take the
    same `np.random.default_rng(0).integers(0, k, n)` base configuration and the sub-register
    is the same big-endian odometer over `residues`.  The only change is that the (4096, n, 3)
    CA table is built once, in a single batched `build_ca` call, instead of per readout.
    """

    def __init__(self, pdb, residues=SUB):
        import itertools
        self.e = GEnum(pdb)
        e = self.e
        self.pdb = pdb
        self.nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)   # ORACLE, post-hoc
        self.res = tuple(residues)
        self.n_qubits = 2 * len(self.res)
        base = np.random.default_rng(0).integers(0, e.k, e.n)
        cfg = np.array(list(itertools.product(range(e.k), repeat=len(self.res))), np.int64)
        S = np.tile(base, (len(cfg), 1))
        S[:, list(self.res)] = cfg
        self.full = e.index(S)
        self.rmsd = e.rmsd[self.full]                    # ORACLE, post-hoc scoring only
        self.M = len(self.full)
        r = np.arange(e.n)
        phi = e.PHI[r[None, :], S]
        psi = e.PSI[r[None, :], S]
        self.CA = np.asarray(I.build_ca(phi, psi), float)          # (M, n, 3)
        # independent check that the rebuilt trace is the one the stored RMSD scored
        chk = I.kabsch_rmsd_batch(self.CA[:256], self.nat)
        self.rebuild_maxerr = float(np.abs(np.asarray(chk) - self.rmsd[:256]).max())

    def build(self, idx):
        """(B, n, 3) CA traces -- the same object `qgeom_cvar.Struct.build` returns."""
        return self.CA[np.asarray(idx, np.int64)]

    # ------------------------------------------------------------- readout primitives
    def coord_avg_rmsd(self, idx):
        """ORACLE post-hoc: CA-RMSD of the COORDINATE AVERAGE of a set of configurations."""
        idx = np.unique(np.asarray(idx, np.int64))
        W = self.CA[idx]
        if len(W) == 1:
            return float(I.kabsch_rmsd_batch(W, self.nat)[0])
        A = I.coordinate_average(W)
        A = A[0] if isinstance(A, tuple) else A
        return float(I.kabsch_rmsd_batch(np.asarray(A, float)[None], self.nat)[0])

    def mean_pairwise(self, idx, cap=40, rng=None):
        """ORACLE-free diversity: mean pairwise CA-RMSD inside a subsample of the set."""
        idx = np.unique(np.asarray(idx, np.int64))
        if idx.size < 2:
            return 0.0
        rng = rng or np.random.default_rng(0)
        if idx.size > cap:
            idx = rng.choice(idx, size=cap, replace=False)
        W = self.CA[idx]
        d = [float(x) for i in range(len(W) - 1)
             for x in I.kabsch_rmsd_batch(W[i + 1:], W[i])]
        return float(np.mean(d))


# ======================================================== the native-free objective
def tabulate_any(pdb, chunk=65536):
    """Full-enumeration tabulation of both hamil terms, for EITHER enumeration set.

    Returns the identical `(prior, disto)` float32 pair `s14.vqe_hamil.tabulate` returns, and
    REUSES that module's cache verbatim when it exists (the nine n=9 targets), so those nine
    objectives are the same bytes `qgeom_ens` consumed.  The ten n=10 targets are tabulated
    here on the same dtype and cached under `s15/cache/`.
    """
    p9 = os.path.join(V.CACHE, f"hamil_full_{pdb}.npz")
    if os.path.exists(p9):
        d = np.load(p9)
        return d["prior"], d["disto"]
    p = os.path.join(CACHE, f"qens_hamil_full_{pdb}.npz")
    if os.path.exists(p):
        d = np.load(p)
        return d["prior"], d["disto"]
    from s14 import hamil as H
    e = GEnum(pdb)
    t = H.Terms(pdb, e.seq, e.n, e.fold, k=e.k)
    pr = np.empty(e.N, np.float32)
    ds = np.empty(e.N, np.float32)
    for a in range(0, e.N, chunk):
        b = min(a + chunk, e.N)
        S = e.states(np.arange(a, b))
        pr[a:b] = t.e_prior(S)
        ds[a:b] = t.e_disto(S)
    np.savez_compressed(p, prior=pr, disto=ds)
    return pr, ds


def hamil_sub(st, w=0.25, cache=True):
    """The native-free structural objective on the sub-register, uniformised.

    Literally `s14.vqe_hamil.combine(*tabulate(pdb), w)[st.full]` put through
    `vqe_lib.uniformise` -- the same two lines `qgeom_ens` runs -- with `tabulate` widened to
    both enumeration sets.  Cached per target.
    """
    path = os.path.join(CACHE, f"qens_E_{st.pdb}_w{w}_{'-'.join(map(str, st.res))}.npy")
    if cache and os.path.exists(path):
        return np.load(path)
    from s14.vqe_hamil import combine
    pr, ds = tabulate_any(st.pdb)
    E = V.uniformise(combine(pr, ds, w)[st.full])
    if cache:
        np.save(path, E)
    return E


# ============================================================== the arms' distributions
def boltzmann_tilt(p0, E, target_entropy_bits, lo=1e-9, hi=1e4, iters=120):
    """CLASSICAL control: p ~ p0 * exp(-E/T), T solved so H(p) == the VQE's final entropy.

    Native-free: uses only the objective and the untrained circuit's own distribution.  This
    is the cheapest possible way to move mass toward low objective values while holding
    diversity fixed, and it is the arm that decides whether the variational optimisation
    contributed anything beyond "tilt the distribution toward the objective".
    """
    p0 = np.asarray(p0, float)
    E = np.asarray(E, float)

    def H_of(T):
        z = -(E - E.min()) / max(T, 1e-30)
        q = p0 * np.exp(z - z.max())
        s = q.sum()
        if not np.isfinite(s) or s <= 0:
            return 0.0, p0 / p0.sum()
        q = q / s
        m = q > 0
        return float(-(q[m] * np.log2(q[m])).sum()), q

    h_hi = H_of(hi)[0]
    if target_entropy_bits >= h_hi:          # cannot be broader than p0 itself
        return H_of(hi)[1]
    a, b = lo, hi
    for _ in range(iters):
        m = np.sqrt(a * b)
        if H_of(m)[0] < target_entropy_bits:
            a = m
        else:
            b = m
    return H_of(np.sqrt(a * b))[1]


def anneal_ensemble(E, n, k, budget, rng, sub_full=None):
    """CLASSICAL simulated annealing, consumed as an UNRANKED ensemble of visited configs.

    `search_anneal` works in the (n, k) product space of the SUB-register, so `n` here is
    `len(st.res)` and the returned indices are sub-register indices in the same big-endian
    order as `Struct.full`.
    """
    c = V.search_anneal(E, n, k, budget, rng)
    return c.all_seen().astype(np.int64), c.used


# ===================================================================== the readouts
MS = (5, 20, 75)


def readouts(st, idx, rng, ms=MS, dmatch=None, setcap=200, light=False):
    """Everything a downstream operator could consume from a drawn set WITHOUT a ranker.

    `dmatch` -- if given, additionally emit `mdN_coordavg_rmsd`: the coordinate average of
    `N` members drawn uniformly from a `dmatch`-sized random subsample of the set's DISTINCT
    configurations.  That is the equal-diversity, equal-budget readout.
    """
    idx = np.asarray(idx, np.int64)
    R = st.rmsd                                     # ORACLE, post-hoc scoring only
    u = np.unique(idx)
    out = {"set_mean_rmsd": float(R[idx].mean()),
           "set_best_rmsd": float(R[idx].min()),
           "n_drawn": int(idx.size),
           "n_distinct": int(u.size)}
    cnt = np.bincount(idx, minlength=st.M).astype(float)
    q = cnt / cnt.sum()
    m = q > 0
    out["draw_entropy_bits"] = float(-(q[m] * np.log2(q[m])).sum())
    # `set_coordavg_rmsd` keeps the ORIGINAL definition from `qgeom_ens.unranked_readouts`
    # so the reproduction is exact.  DEFECT, found here and recorded: `np.unique` returns
    # configuration indices SORTED, so `u[:cap]` is the `cap` LOWEST-INDEXED configurations,
    # not a random subsample -- a systematic slice of the register (residue 1 in its low
    # states).  It bites differently on arms with different distinct counts, so it is not a
    # valid cross-arm comparison.  `setrand_coordavg_rmsd` is the corrected, unbiased version
    # and is the one this workstream reads.
    # `light` skips the two `setcap`-sized coordinate averages, which cost about half of a
    # readout call.  Used by the alpha MAP (R3), where the `rand-m` family already covers the
    # same ground; never used by the primary comparison (R1).
    if not light:
        out["set_coordavg_rmsd"] = st.coord_avg_rmsd(u[:setcap] if u.size > setcap else u)
        out["setrand_coordavg_rmsd"] = st.coord_avg_rmsd(
            rng.choice(u, size=setcap, replace=False) if u.size > setcap else u)
    for mm in ms:
        pick = rng.choice(idx, size=mm, replace=idx.size < mm)
        out[f"rand{mm}_coordavg_rmsd"] = st.coord_avg_rmsd(pick)
        out[f"rand{mm}_mean_rmsd"] = float(R[pick].mean())
    out["mean_pairwise_rmsd"] = st.mean_pairwise(idx, rng=rng)
    if dmatch is not None:
        out.update(md_readouts(st, idx, rng, dmatch, ms=ms))
    return out


def md_readouts(st, idx, rng, dmatch, ms=MS):
    """The MATCHED-DIVERSITY readout: equal distinct-configuration count, equal budget.

    Subsample the set's DISTINCT configurations to a common size `dmatch` (set by the caller
    to the minimum over the arms being compared) and coordinate-average `m` members drawn
    uniformly from that subsample.  This removes both the diversity difference and the
    multiplicity (collapse) difference between arms, leaving only WHICH configurations each
    distribution favours.
    """
    R = st.rmsd                                     # ORACLE, post-hoc scoring only
    u = np.unique(np.asarray(idx, np.int64))
    out = {}
    if not u.size:
        return out
    D = min(int(dmatch), int(u.size))
    pool = rng.choice(u, size=D, replace=False)
    for mm in ms:
        pick = rng.choice(pool, size=mm, replace=D < mm)
        out[f"md{mm}_coordavg_rmsd"] = st.coord_avg_rmsd(pick)
        out[f"md{mm}_mean_rmsd"] = float(R[pick].mean())
    out["dmatch"] = D
    return out


# ===================================================================== statistics
def paired(a, b, seed=0):
    return G.paired(a, b, seed=seed)


def nullconc(d, sims=4000, seed=0):
    from s15.qgeom_nullconc import null_calibrated
    return null_calibrated(d, sims=sims, seed=seed)


def verdict_line(name, v, c, seed=0):
    """One row of the mandated report: paired mean, CI, W/L, median, mean/sd, concentration."""
    pr = paired(v, c, seed=seed)
    if pr.get("n", 0) < 6:
        return {**pr, "conc_verdict": "TOO SMALL"}
    d = np.asarray(v, float) - np.asarray(c, float)
    nc = nullconc(d, seed=seed)
    return {**pr, "vqe_mean": float(np.mean(v)), "ctrl_mean": float(np.mean(c)),
            "conc_verdict": nc["verdict"], "conc_p_share": nc["p_share_vs_null"],
            "conc_p_drop": nc["p_drop_vs_null"], "conc_power_low": nc["power_warning"]}


def fmt(name, r, width=34):
    return (f"  {name:<{width}s} {r['n']:4d} {r['mean']:+8.4f} "
            f"[{r['ci_lo']:+7.4f},{r['ci_hi']:+7.4f}] {r['median']:+8.4f} "
            f"{r['win']:4d}/{r['loss']:<4d} {r['mean_over_sd']:+6.2f} "
            f"{r['verdict']:>11s} {r.get('conc_verdict', ''):>19s}")


HDR = (f"  {'comparison':<34s} {'n':>4s} {'mean':>8s} {'CI95':>18s} {'median':>8s} "
       f"{'W/L':>9s} {'m/sd':>6s} {'verdict':>11s} {'concentration':>19s}")
