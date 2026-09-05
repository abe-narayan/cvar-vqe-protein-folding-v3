"""`core.pipeline` -- the four-stage production path, consolidated, parallel and resumable.

WHAT THIS IS
============
`s9/final.py` is the authoritative reference implementation of the shipped predictor.  It
is also a sprint script: one process, one target at a time, three separate CLI stages each
walking the whole manifest, and a pool cache that stores a 500x500 matrix of which 5625
entries are ever read.  This module is the same science with none of that.

    STAGE 1  RETRIEVE     out-of-fold peptides + `_fold_fragments(fold, 5)`, windowed at
                          the target length; K = 500 highest BLOSUM62 sums, taken with a
                          STABLE argsort (the similarity has large tie sets and a
                          different tie-break gives a different pool -- this is pinned).
    STAGE 2  FILTER       the shipped distogram Bayes-risk score's top M = 75.
    STAGE 3  SYNTHESISE   coordinate-average the 75 after superposing on their medoid,
                          then project back onto the manifold of ideal-geometry chains,
                          MULTI-START at every rung, hinged Ramachandran prior at
                          lambda = 0.3.
    STAGE 4  RELAX        one restrained ff14SB/GBn2 minimisation, k = 10 kcal/mol/A^2,
                          steps = 0 (converged).  A VALIDITY stage: it costs accuracy and
                          earns its place by removing builder strain.

Every constant above is inherited from the pre-registration in `s9/final.py` and NOTHING
here may be tuned against a reported number.

WHAT IS DIFFERENT, AND WHY EACH DIFFERENCE IS EXACT
===================================================
1. THE PAIRWISE RMSD MATRIX IS BUILT ON THE FILTERED SUBSET ONLY.
   `s9.final.build_pool` computes all K^2 = 250,000 superpositions; `synthesise` then
   reads `P[np.ix_(sub, sub)]` and nothing else ever reads P on the deployable path.  So
   45x of that stage is dead work.  `kabsch_rmsd_batch` is elementwise over its batch axis
   -- one 3x3 SVD, one determinant and one per-structure reduction per member, with no
   cross-batch reduction -- so the sub-block computed directly is BIT-IDENTICAL to the
   sub-block sliced out of the full matrix.  `tests/test_pipeline.py` asserts that at 0.0
   on real pools rather than arguing it.

2. THE LIBRARY AND ITS WINDOWS ARE MEMOISED PER (FOLD, LENGTH).
   `library_members(seq, fold)` filters `db.load()` by `folds[q.seq] != fold and q.seq !=
   target_seq`.  The second clause is unreachable: `folds` is keyed BY SEQUENCE, so any
   member sharing the target's sequence is in the target's own fold and the first clause
   already dropped it.  The library is therefore a function of the fold alone, and its
   windows a function of (fold, length).  `_Library` memoises on exactly that key and
   re-verifies the unreachability per target before it serves a cached entry, so a future
   change to `folds()` degrades to recomputation instead of to a wrong pool.

3. TARGETS RUN IN PARALLEL.  Targets are independent; this is the single largest
   structural win available and it changes no number.  See THREADING below.

4. ONE PASS INSTEAD OF THREE.  `s9.final` walks the manifest once per CLI stage and
   round-trips a 500x500 float32 matrix through `savez_compressed` between them.  Here a
   target is retrieved, filtered, synthesised and relaxed inside one function, and only
   the result is serialised.

THREADING -- THE OVERSUBSCRIPTION TRAP
======================================
Every worker imports numpy (BLAS), torch (the distogram predictor) and OpenMM.  With W
workers each spawning T threads the box runs W*T runnable threads on 8 cores and the
measured throughput FALLS.  Workers therefore pin `OMP_NUM_THREADS`, `MKL_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `NUMEXPR_NUM_THREADS` and `torch.set_num_threads` to
`Config.threads_per_worker` (default 1) BEFORE numpy is imported in the child.  OpenMM is
already pinned to one thread by `amber_refine` for determinism, and that is load-bearing
for reproducibility, not just for speed.

NATIVES
=======
`nat_ca`, `rr` and `Dnat` are REPORTING LABELS.  The deployable path sees a
`deployable_view` and nothing else; `label()` is the only function that opens a native and
it runs after the structure is final.  `tests/test_pipeline.py` NaN-poisons every native
quantity and asserts every emitted coordinate is bit-identical.

RESUMABILITY
============
Every target writes one atomic JSON checkpoint under
``bench_results/cache/<config-key>/<pdb>.json``.  The config key is a SHA-1 over every
parameter that can change a number -- K, M, penalty, lambda, the AMBER schedule, the fold
count, the pair separation, the optimiser iteration cap, the tie-break rule and the set of
live backends.  A run that dies at 99% resumes at 99%; a run with a changed parameter gets
a different key and cannot silently inherit the old answer.

    python -m core.pipeline run   --manifest tuning126 --workers 6
    python -m core.pipeline stats --manifest tuning126
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from collections import OrderedDict
from dataclasses import dataclass, asdict, field, replace

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np                                                       # noqa: E402

import core                                                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(_ROOT, "bench_results")
CACHE_ROOT = os.path.join(RESULTS, "cache")

#: The per-stage clock keys, in report order.  `startup` and `imports` are measured by the
#: harness (they are properties of the process, not of a target); the rest accumulate here.
STAGES = ("retrieval", "windows", "model_load", "distogram", "filter", "pairwise",
          "average", "projection", "quantum", "legacy", "amber", "amber_wait",
          "reference_arms", "reporting", "io", "serialise")

#: The leave-fold-out (alpha, T) for the CVaR selector, PINNED from `s8/integrate_vqe.json`
#: (`vqe_LFO`, 3.3135 A against the shipped 3.4540 on this same instrument).  It is a
#: LEAVE-FOLD-OUT table -- fold f's cell was chosen on the other four folds -- so reading
#: it here is not tuning on the instrument, and picking the grid's best cell instead
#: (a0.25_T0.3, 3.2825) would be.
VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}


# ============================================================ configuration
@dataclass(frozen=True)
class Config:
    """Every parameter that can change a number.  Frozen, hashed into the cache key.

    The defaults ARE the pre-registration.  A development mode may change them but must
    set `dev_mode` so the label travels with the result and the harness refuses to call it
    a headline.
    """
    k: int = 500                     #: BLOSUM62 retrieval depth
    m: int = 75                      #: distogram-score filter size
    penalty: str = "ramah"           #: the hinged Ramachandran prior
    lam: float = 0.3                 #: its weight
    amber_k: float = 10.0            #: restraint force constant, kcal/mol/A^2
    amber_steps: int = 0             #: 0 = minimise to convergence (OpenMM convention)
    amber: bool = True               #: run stage 4 at all
    n_folds: int = 5
    min_sep: int = 2                 #: CA-CA pair separation for the distogram score
    maxiter: int = 300               #: L-BFGS-B cap inside the projection
    multi_start: bool = True         #: multi-start at every rung -- load-bearing (S9-2)
    #: WHICH PROJECTION OBJECTIVE STAGE 3b MINIMISES, and it is a science parameter.
    #: `exact` is the reference: the bit-exact builder and the reference's own
    #: finite-difference gradient, reproducing `s8.project` bit-for-bit.  `analytic` is
    #: 4.6x faster and lands on DIFFERENT STRUCTURES -- 126/126 targets move, median
    #: 0.042 A, worst 1.90 A -- because the projection is degenerate and a 1e-13 A
    #: difference in the forward map routes some L-BFGS-B trajectories into the other
    #: torsion branch.  It is therefore a different pipeline and is hashed as one; a run
    #: under one setting can never be served out of the other's cache.
    project_grad: str = "exact"
    tie_break: str = "stable"        #: `np.argsort(kind=...)`; PINNED, changes the pool
    torsion_window: int = 8          #: `torsion_lib2.library_for` window
    #: THE REFERENCE PATH'S STORAGE PRECISION, and it is not cosmetic.
    #: `s9.final` writes its pool to a float32 `npz` between the CLI stages and
    #: `load_pool` reads it back as float64, so everything downstream of retrieval --
    #: the filter order, the medoid, the coordinate average, the projection -- sees
    #: coordinates and scores that have been ROUND-TRIPPED THROUGH FLOAT32.  This module
    #: has no such round trip, so reproducing the reference bit-for-bit means quantising
    #: at the same boundary.  Setting this False keeps full float64 through the pipeline;
    #: it is a numerical-precision CHANGE and `tests/test_pipeline.py` measures what it is
    #: worth before anyone may claim it as an improvement.
    reference_precision: bool = True
    #: Also emit the SINGLE-START projection from the medoid's own torsions.
    #: The project quotes two different numbers both called "the incumbent synthesis":
    #: `s9/synth_probe.json`'s `fit` = 3.2005, which is `s9.synth.fit_w` from the medoid's
    #: torsions, and `s9/synth.py`'s own `EXPECT["fit"]` = 3.204, which is the MULTI-START
    #: `lam_path[0.0]` that `s9/final.py` actually emits.  They are different
    #: constructions -- all 126 targets differ, by up to 0.339 A -- and reproducing only
    #: one of them would leave the other looking like a discrepancy.  This arm costs one
    #: extra L-BFGS solve per target and makes both reproducible.
    #: DEFAULT OFF, and that default is the point.  It is work `s9/final.py` does not do,
    #: so a run with it on is not like-for-like: it charged the optimised arm 70.8 CPU-s
    #: against the reference's 0.0 and made the measured speedup an understatement of an
    #: experiment nobody ran.  Turn it on for a declared VERIFICATION pass -- which is how
    #: `bench_results/optimised_tuning126_w4.json` reproduced 3.2007 -- never for a
    #: headline timing run.
    report_single_start_fit: bool = False

    # ---- the four mandated components.  `s9/final.py` is retrieve -> filter ->
    # synthesise -> AMBER: AMBER participates, Legacy is only imported for scoring, and
    # VQE/CVaR do not participate at all.  These two switches put the other two into the
    # roles the project's OWN measurements give them, and every arm is reported as an
    # ablation -- including where the contribution is zero, which for several of these
    # is the established scientific result and not a reason to omit them.
    quantum: bool = False            #: CVaR-VQE over the discrete hypothesis set
    legacy: bool = False             #: Legacy as a LATE refiner on an already-tight set
    vqe_qubits: int = 7              #: 2^7 = 128 hypotheses (s8/integrate.py's size)
    vqe_layers: int = 3
    vqe_iters: int = 50
    vqe_seed: int = 0
    legacy_top: int = 5              #: `refine_LegacyAll_top5`, the best-measured cell

    threads_per_worker: int = 1
    #: Pin each worker to its own core.  OFF by default, and the default is a measurement,
    #: not a guess: this box is an 8-core Intel Core Ultra 7 256V with NO hyperthreading
    #: and NON-UNIFORM cores.  Normalised per-core throughput on a fixed cache-resident
    #: workload measures 0.79 / 1.00 / 0.69 / 0.58 / 0.75 / 0.79 / 0.82 / 0.81 -- 6.22
    #: fastest-core-equivalents, not 8 -- so a fixed pin can strand a worker on a core
    #: worth 0.58 of another for the whole run, which is worse than letting the scheduler
    #: migrate it.  The switch exists so the hypothesis can be MEASURED; `affinity_cores`
    #: takes an explicit core list once someone has measured which ones are fast.
    affinity: bool = False
    affinity_stride: int = 1
    affinity_cores: tuple = ()
    #: How many window banks and fold libraries a worker may hold resident.  NOT science:
    #: an evicted entry is RECOMPUTED, never approximated.  See `_Library` for why the
    #: caps are small -- unbounded, the window banks alone were 435 MB per worker.
    window_cache: int = 2
    member_cache: int = 2
    #: How many AMBER minimisations may run AT ONCE, across all workers.  0 = no cap.
    #: The stages have very different memory profiles: GBn2 is ~98% of the AMBER call and
    #: is bandwidth-bound, while the projection is cache-resident numpy.  Eight workers
    #: all landing in AMBER together saturate memory bandwidth -- its CPU-work is +30% at
    #: w4 and +42% at w8 for identical science -- so admitting fewer at a time, and
    #: letting the rest do projection meanwhile, can lower TOTAL cpu-work by more than the
    #: lost concurrency costs.  This changes only WHEN work runs, never what is computed,
    #: so it cannot move a number and is not part of the cache key.
    amber_concurrency: int = 0
    dev_mode: str = ""               #: non-empty => not a headline number

    #: not part of the science, so not part of the key
    stop_pct: int = 92
    start_pct: int = 92

    #: Fields that cannot change a number, and so must NOT change the cache key: thread
    #: counts, core affinity, memo sizes and the memory gate.  An evicted memo entry is
    #: recomputed identically, and a pinned worker computes what an unpinned one does --
    #: `tests/test_pipeline.py` asserts both rather than taking them on trust.  Everything
    #: else in `Config` is in the key, including `report_single_start_fit`, because that
    #: one changes what the record CONTAINS.
    NOT_SCIENCE = ("threads_per_worker", "stop_pct", "start_pct", "affinity",
                   "affinity_stride", "affinity_cores", "window_cache", "member_cache",
                   "amber_concurrency")

    def science(self) -> dict:
        d = asdict(self)
        for k in self.NOT_SCIENCE:
            d.pop(k, None)
        return d

    def key(self, backends: dict | None = None) -> str:
        payload = {"science": self.science(),
                   "backends": backends if backends is not None else core.backend_report()}
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha1(blob.encode()).hexdigest()[:16]


PROD = Config()
#: A clearly-labelled development configuration.  EXCLUDED from every headline.
SMOKE = Config(dev_mode="smoke")


# ============================================================ resource discipline
def mem_pct() -> int:
    """Percent of physical RAM in use, from one `ctypes` GlobalMemoryStatusEx read.

    Not `Get-CimInstance`: on a loaded box that costs minutes per call, which turns a
    resource gate into the thing consuming the resource.
    """
    import ctypes

    class _MS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    st = _MS()
    st.dwLength = ctypes.sizeof(_MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return int(st.dwMemoryLoad)


def mem_bytes() -> tuple:
    import ctypes

    class _MS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    st = _MS()
    st.dwLength = ctypes.sizeof(_MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return int(st.ullTotalPhys), int(st.ullAvailPhys), int(st.dwMemoryLoad)


def wait_for_memory(stop_pct=92, timeout=300.0, poll=5.0, where=""):
    """Pause until the box is below the gate.  Returns False if it never gets there.

    A hard abort at the gate is the wrong behaviour on a shared box: the run dies, and the
    next invocation pays the process tax again to make the same three targets of progress.
    Waiting is free, the gate still holds, and every stage is checkpointed anyway -- so the
    worst case is that a run is slow, not that it is lost.  The caller still gets False
    after `timeout` so a genuinely wedged box does not hang a benchmark forever.
    """
    t0 = time.perf_counter()
    m = mem_pct()
    if m < stop_pct:
        return True
    drop_caches()
    while time.perf_counter() - t0 < timeout:
        m = mem_pct()
        if m < stop_pct:
            if where:
                print(f"  [mem {m}% -- resuming {where} after "
                      f"{time.perf_counter() - t0:.0f}s]", flush=True)
            return True
        time.sleep(poll)
    return False


def proc_rss():
    """This process's current and PEAK working set, in MB, via GetProcessMemoryInfo.

    Per-worker footprint is the quantity that decides how many workers fit on this box,
    so it is measured in every worker and carried in every per-target record rather than
    inferred from the system-wide RAM trace -- which cannot tell a worker's private pages
    from the page cache the workers share.
    """
    import ctypes

    class PMC(ctypes.Structure):
        _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t)]
    k32 = ctypes.WinDLL("kernel32")
    ps = ctypes.WinDLL("psapi")
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    ps.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(PMC),
                                        ctypes.c_ulong]
    c = PMC()
    c.cb = ctypes.sizeof(PMC)
    if not ps.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb):
        return None
    return {"rss_mb": c.WorkingSetSize / 2 ** 20,
            "peak_rss_mb": c.PeakWorkingSetSize / 2 ** 20}


def drop_caches():
    import gc
    try:
        core.backend("amber").clear_cache()
    except Exception:                                                    # noqa: BLE001
        pass
    _LIB.clear()
    gc.collect()


# ============================================================ atomic io
def _atomic_write_json(path, obj):
    tmp = path + ".tmp"
    for attempt in range(5):
        try:
            with open(tmp, "w") as fh:
                json.dump(obj, fh)
            os.replace(tmp, path)
            return
        except OSError:                                                  # noqa: PERF203
            if attempt == 4:
                raise
            time.sleep(0.5 + attempt)


def _read_json(path, default=None):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:                                                    # noqa: BLE001
        return default


# ============================================================ manifests
def manifest(name: str):
    """A DECLARED target list.  Never filtered on anything a result depends on.

    tuning126   the 126-target tuning instrument (`s7.debias.tuning_targets`).  This is
                where all optimisation work happens.
    dev24       the 24-target dev set.
    benchmark60 the held-out benchmark.  ITS SINGLE PRE-REGISTERED PASS IS SPENT -- the
                harness refuses to run it without `--i-am-spending-the-benchmark`, and it
                is never a tuning instrument.
    smoke8/24   DECLARED development subsets: the first 8 / 24 of tuning126 in pinned pdb
                order.  Labelled `dev_mode` everywhere; excluded from every headline.
    """
    db = core.backend("data")
    from s7 import debias
    if name == "tuning126":
        return list(debias.tuning_targets())
    if name == "dev24":
        return list(db.dev_set(24))
    if name == "benchmark60":
        return list(db.benchmark())
    if name.startswith("smoke"):
        n = int(name[5:] or 8)
        return list(debias.tuning_targets())[:n]
    raise KeyError(f"unknown manifest {name!r}")


MANIFESTS = ("tuning126", "dev24", "benchmark60", "smoke8", "smoke24")


# ============================================================ the library, memoised
class _Library:
    """Per-fold library members and per-(fold, length) windows, under a HARD memory cap.

    `s9.final` rebuilds both for every target.  They depend only on the fold (see the
    module docstring, point 2), so a run over 126 targets pays for a handful of libraries
    and window banks instead of 126 of each.  `check(seq, fold)` re-derives the claim per
    target rather than trusting it.

    WHY THE CAPS ARE SMALL, AND WHY THAT IS THE WHOLE POINT
    ------------------------------------------------------
    Cross-target parallelism on this box is MEMORY-bound, not core-bound: at 8 workers the
    run drops to 40% CPU at 99% RAM, which is swapping, and it is 2.1x SLOWER than 4
    workers.  So what each worker HOLDS decides how many workers fit, and an unbounded
    memo is a throughput bug wearing a speedup's clothes.  Measured, per worker:

        python + core.pipeline import                 31 MB
        the heavy closure (torch, the database)      ~230 MB   unavoidable, per process
        library members, all five folds               ~89 MB   17.8 MB of fragments each
        window banks, all 40 (fold, length) entries   435 MB   <-- the whole problem
        distogram fold models, all five                31 MB
        one target's AMBER context                   193 MB
                                                     -------
                                                     ~967 MB

    Targets are dispatched in (fold, length) order, so a worker that keeps the two most
    recent banks hits the cache on essentially every target and holds ~25 MB instead of
    435.  Capping is exact -- an evicted bank is RECOMPUTED, never approximated -- so it
    trades a little arithmetic for a lot of resident memory, which is the correct trade
    when memory is what is capping the worker count.
    """

    def __init__(self, n_folds=5, cap=2, cap_members=2):
        self.n_folds = n_folds
        self.cap = cap
        self.cap_members = cap_members
        self._members = OrderedDict()
        self._windows = OrderedDict()
        self._seqs = {}
        self.stats = {"window_hits": 0, "window_misses": 0, "window_evictions": 0,
                      "member_hits": 0, "member_misses": 0, "member_evictions": 0}

    def clear(self):
        self._members.clear()
        self._windows.clear()
        self._seqs.clear()

    def members(self, fold):
        if fold in self._members:
            self._members.move_to_end(fold)
            self.stats["member_hits"] += 1
            return self._members[fold]
        self.stats["member_misses"] += 1
        db = core.backend("data")
        dgm = core.backend("predict")
        folds = db.folds(self.n_folds)
        peps = [q for q in db.load() if folds[q.seq] != fold]
        frags = list(dgm._fold_fragments(fold, self.n_folds))
        while len(self._members) >= self.cap_members:
            old, _ = self._members.popitem(last=False)
            self._seqs.pop(old, None)
            self.stats["member_evictions"] += 1
        self._members[fold] = (peps, frags)
        self._seqs[fold] = {q.seq for q in peps}
        return self._members[fold]

    def check(self, target_seq, fold):
        """The memo is only valid if the target's own sequence is absent from the library.

        `s9.final.library_members` drops it explicitly; `folds` being keyed by sequence
        makes that drop unreachable.  If it ever becomes reachable this returns False and
        the caller recomputes for that target, so the memo can never make a wrong pool.
        """
        self.members(fold)
        return target_seq not in self._seqs[fold]

    def windows(self, fold, n):
        key = (fold, int(n))
        if key in self._windows:
            self._windows.move_to_end(key)
            self.stats["window_hits"] += 1
            return self._windows[key]
        self.stats["window_misses"] += 1
        peps, frags = self.members(fold)
        bank = windows_all(peps + frags, int(n))
        while len(self._windows) >= self.cap:
            self._windows.popitem(last=False)
            self.stats["window_evictions"] += 1
        self._windows[key] = bank
        return bank

    def nbytes(self):
        """Resident bytes in the two memos -- what this worker is holding, measured."""
        w = sum(sum(a.nbytes for a in v if hasattr(a, "nbytes"))
                for v in self._windows.values())
        return {"window_banks_mb": w / 2 ** 20, "n_window_banks": len(self._windows),
                "n_member_folds": len(self._members)}


_LIB = _Library()


def windows_all(pool, n):
    """Every length-n window of every library member: CA, phi, psi, encoded seq, origin.

    Byte-for-byte the iteration order of `s7.audit.windows_of`, `s8.generate._windows_all`
    and `s9.final.windows_all`.  That order IS the tie-break: the BLOSUM similarity has
    large tie sets and `argsort(kind="stable")` resolves them by position, so a different
    iteration order is a different pool.
    """
    aud = core.backend("numerics")
    cas, phis, psis, seqs, src = [], [], [], [], []
    for q in pool:
        m = len(q.seq)
        if m < n:
            continue
        e = aud.encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n])
            phis.append(np.asarray(q.phi, float)[s:s + n])
            psis.append(np.asarray(q.psi, float)[s:s + n])
            seqs.append(e[s:s + n])
            src.append(q.pdb)
    return (np.stack(cas), np.stack(phis), np.stack(psis), np.stack(seqs),
            np.array(src, dtype=object))


# ============================================================ the ESM guard
_GUARDED = [False]


def guard_esm(real_seqs=()):
    """Refuse to materialise the 1.5 GB `esm_cache.npz`.  Same guard as `s7.poolsize`.

    `train_fold` probes the feature width with ONE arbitrary TRAINING sequence which is
    not in the 21 MB hot cache; `esm_features.raw` would then load the full cache and take
    the box past its memory gate.  Every fold model is already on disk, so nothing trains
    and the probe's VALUES are never used -- only its shape, which fixes `d_in`.  The probe
    is served zeros of the right shape; every real (cached) sequence goes through
    unchanged; a MANIFEST sequence missing from the hot cache is a hard error, never a
    silent zero.
    """
    if _GUARDED[0]:
        return
    dgm = core.backend("predict")
    import esm_features as ef
    for f in range(_LIB.n_folds):
        path = dgm._model_path(f, True, True, 0)
        if not os.path.exists(path):
            raise RuntimeError(f"missing fold model {path}; this path must not train")
    hot = ef._load_small()
    real = set(real_seqs)
    _orig = ef.raw

    def raw(sequence):
        if sequence in hot:
            return _orig(sequence)
        if sequence in real:
            raise RuntimeError(f"manifest sequence {sequence!r} absent from the ESM hot "
                               f"cache; refusing to fall back to the 1.5 GB cache")
        return (np.zeros((len(sequence), 1280), np.float32),
                np.zeros((len(sequence), len(sequence)), np.float32))

    ef.raw = raw
    _guard_core_data(real)
    _GUARDED[0] = True


def _guard_core_data(real):
    """The same guard for the consolidated ESM provider, if that is the one in use.

    `core.predict` reads its features through `core.data.esm_raw`, not through
    `esm_features.raw`, so patching only the latter would leave the 1.5 GB bank one probe
    away.  `feature_width` probes with the literal sequence "ACDEFGHIKL", which is in
    neither the hot cache nor the bank, and `core.data.esm_raw` answers that by extracting
    from the bank in a subprocess and then running ESM-2 -- minutes of work and gigabytes
    of memory to fix `d_in`, a number that only depends on the SHAPE.  So the probe is
    served zeros of the right shape here too, and a manifest sequence missing from the hot
    cache stays a hard error rather than a silent zero vector.
    """
    if core.backend_name("predict") != "core.predict":
        return
    try:
        from core import data as cdata
    except ImportError:
        return
    if getattr(cdata, "_pipeline_guarded", False):
        return
    _orig = cdata.esm_raw

    def esm_raw(sequence):
        if sequence in cdata._hot():
            return _orig(sequence)
        if sequence in real:
            raise RuntimeError(f"manifest sequence {sequence!r} absent from the ESM hot "
                               f"cache; refusing to fall back to the 1.5 GB bank")
        d = _probe_width(cdata)
        return (np.zeros((len(sequence), d), np.float32),
                np.zeros((len(sequence), len(sequence)), np.float32))

    cdata.esm_raw = esm_raw
    cdata._pipeline_guarded = True


def _probe_width(cdata):
    """The raw ESM width, taken from the hot cache rather than assumed to be 1280."""
    hot = cdata._hot()
    for v in hot.values():
        return int(np.asarray(v[0]).shape[1])
    return 1280


_MODELS = {}


def fold_model(fold, n_folds=5):
    """The fold model, loaded off disk once per (process, fold).  Trains nothing.

    `train_fold` LOADS the model (all five exist on disk); the guard above makes certain
    the width probe cannot touch the 1.5 GB ESM bank on the way.
    """
    dgm = core.backend("predict")
    guard_esm()
    if fold not in _MODELS:
        _MODELS[fold] = dgm.train_fold(fold, True, n_folds, fragments=True, verbose=False)
    return _MODELS[fold]


def distogram_risk(seq, fold, n_folds=5, model=None):
    """The fold model that never saw this target, applied to its sequence.

    Only `_risk` is needed -- the shipped Bayes-risk lookup that `score_risk` gathers into.
    """
    dgm = core.backend("predict")
    d = dgm.Distogram.for_target(seq, model=model or fold_model(fold, n_folds))
    return np.asarray(d._risk, float)


# ============================================================ the stages
class Clock:
    """A per-target stage accumulator.  `with clk('projection'): ...`"""

    def __init__(self):
        self.t = {k: 0.0 for k in STAGES}
        self._k = None
        self._t0 = 0.0

    def __call__(self, key):
        self._k = key
        return self

    def __enter__(self):
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.t[self._k] = self.t.get(self._k, 0.0) + (time.perf_counter() - self._t0)
        return False

    def total(self):
        return float(sum(self.t.values()))


def _q(a, cfg: Config):
    """The reference path's float32 storage round trip (see `Config.reference_precision`)."""
    a = np.asarray(a, float)
    return a.astype(np.float32).astype(float) if cfg.reference_precision else a


def retrieve(target, fold, cfg: Config, clk: Clock):
    """STAGE 1.  The K highest BLOSUM62 sums over every window of the leakage-safe library.

    Returns the deployable pool: coordinates, torsions, encoded sequences, provenance and
    the similarity itself.  Nothing native is read or returned.
    """
    aud = core.backend("numerics")
    #: WINDOW EXTRACTION is clocked apart from the rest of retrieval on purpose.  It is
    #: the only part whose cost scales with the LIBRARY -- 7k to 39k windows per target --
    #: rather than with the pool, and it is exactly the part the per-(fold, length) memo
    #: removes.  Folding it into `retrieval` would hide both the cost and the saving.
    with clk("windows"):
        if _LIB.check(target.seq, fold):
            W, PHI, PSI, S, src = _LIB.windows(fold, target.n)
        else:                                     # the memo's premise failed: recompute
            db = core.backend("data")
            dgm = core.backend("predict")
            folds = db.folds(cfg.n_folds)
            peps = [q for q in db.load()
                    if folds[q.seq] != fold and q.seq != target.seq]
            frags = list(dgm._fold_fragments(fold, cfg.n_folds))
            W, PHI, PSI, S, src = windows_all(peps + frags, target.n)
    with clk("retrieval"):
        sim = aud.B62[S, aud.encode(target.seq)[None, :]].sum(1)
        #: STABLE, ALWAYS.  Stable reproduces all 126 deposited pools with zero
        #: mismatches; an unstable sort moves up to 47 of the 500 members, because the
        #: BLOSUM sum has large tie sets and the iteration order is the tie-break.
        order = np.argsort(-sim, kind=cfg.tie_break)
        idx = order[:cfg.k]
        W64 = np.asarray(W[idx], float)
        pool = {"W64": W64, "W": _q(W64, cfg),
                "PHI": _q(np.asarray(PHI[idx], float), cfg),
                "PSI": _q(np.asarray(PSI[idx], float), cfg),
                "S": S[idx].astype(int), "src": src[idx],
                "sim": _q(np.asarray(sim[idx], float), cfg),
                "n_windows": int(len(W))}
    return pool


def score(pool, seq, fold, n, cfg: Config, clk: Clock):
    """STAGE 2a.  The shipped distogram Bayes-risk score of every pool member.

    The score is computed from the FULL-PRECISION windows and quantised afterwards, which
    is exactly the order `s9.final.build_pool` does it in: it scores the float64 windows
    it just built and only then writes float32 to the cache.
    """
    aud = core.backend("numerics")
    from s7 import debias
    from s7.poolsize import GRID
    #: MODEL LOAD is clocked apart from SCORING.  The fold model is loaded off disk once
    #: per (process, fold) and the score is then a gather; charging the load to the first
    #: target of each fold would make one target in twenty-five look expensive and hide a
    #: fixed per-process cost that parallelism MULTIPLIES rather than divides.
    with clk("model_load"):
        model = fold_model(fold, cfg.n_folds)
    with clk("distogram"):
        risk = distogram_risk(seq, fold, cfg.n_folds, model=model)
        i, j = aud.pair_index(n, cfg.min_sep)
        D = aud.pair_dists(pool["W64"], i, j)
        sc = debias.score_risk(risk, GRID, D)
    pool["D"] = _q(D, cfg)
    pool["sc"] = _q(sc, cfg)
    pool["n_candidates"] = int(len(pool["W64"]))
    pool["n_pairs"] = int(len(i))
    return pool


def filter_pool(pool, cfg: Config, clk: Clock):
    """STAGE 2b.  The score's top M, and the pairwise RMSD matrix the medoid needs.

    The matrix is built on the FILTERED PREFIX, which is the whole of what any downstream
    stage reads (see the module docstring, point 1).  `s9.final` builds all K^2 of it.
    When the quantum stage is on, the prefix is widened to its 2^n hypothesis set, and the
    M x M block the classical synthesis needs is SLICED OUT of it -- one construction, two
    consumers, and the slice is bit-identical to computing it alone (asserted at 0.0 in
    `tests/test_pipeline.py`).

    DUPLICATE CANDIDATES.  Real pools contain byte-identical windows -- the same fragment
    retrieved through two parents -- and each duplicate costs a full row of Kabsch
    superpositions.  The row is computed once per DISTINCT structure and broadcast, which
    changes no number (identical inputs, identical outputs) and is reported as
    `n_distinct` so the saving is visible rather than assumed.
    """
    aud = core.backend("numerics")
    with clk("filter"):
        order = np.argsort(np.asarray(pool["sc"], float), kind=cfg.tie_break)
        want = max(cfg.m, (1 << cfg.vqe_qubits) if cfg.quantum else 0)
        top = order[:min(want, len(order))]
        sub = order[:cfg.m]
    with clk("pairwise"):
        #: the matrix is built from the FULL-PRECISION windows and quantised, because
        #: that is the order the reference builds and then stores it in
        W64top = pool["W64"][top]
        b = len(top)
        keys = {}
        rep_of = np.empty(b, int)
        for a in range(b):
            k = W64top[a].tobytes()
            rep_of[a] = keys.setdefault(k, a)
        Pt = np.zeros((b, b), float)
        for a in range(b):
            if rep_of[a] == a:
                Pt[a] = aud.kabsch_rmsd_batch(W64top, W64top[a])
            else:
                Pt[a] = Pt[rep_of[a]]
        Pt = _q(Pt, cfg)
        pool["n_distinct"] = int(len(keys))
        pool["n_top"] = int(b)
    where = {int(t): a for a, t in enumerate(top)}
    loc = np.array([where[int(t)] for t in sub], int)
    Ps = Pt[np.ix_(loc, loc)]
    Wsub = pool["W"][sub]
    return sub, Wsub, Ps, top, Pt


# ============================================================ the mandated components
def _zrank(x):
    """Standardised rank -- `s8.integrate.zrank`, the currency every channel is in."""
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(r.std(), 1e-12)


def consensus_medoid(Dblock, w=None):
    """The member closest to all the others, optionally under an ensemble weight.

    `s8.integrate.consensus_medoid`, restricted to an already-sliced block.  With `w=None`
    it is exactly S8-8's medoid, which is what makes the weighted arm comparable to it.
    """
    Dblock = np.asarray(Dblock, float)
    v = Dblock.mean(1) if w is None else Dblock @ (w / max(float(np.sum(w)), 1e-12))
    return int(np.argmin(v))


def quantum_stage(pool, top, Pt, fold, cfg: Config, clk: Clock):
    """CVaR-VQE as the SELECTOR over the discrete hypothesis set.  Deployable.

    WHY THIS IS THE HONEST PLACEMENT, and not a stub
    ------------------------------------------------
    The near-native band holds ~2-3 populated structural clusters per target, so the
    quantity a consensus readout needs is a DISTRIBUTION over hypotheses, not a point.
    `s8/integrate.py` established exactly this formulation and measured it: the energy is
    the standardised rank of the shipped distogram score over the top 2^n filtered
    candidates, the state is a 3-layer RY/CNOT ansatz, and the readout is the
    p_theta-weighted consensus medoid.

    WHY CVaR, MEASURED
    ------------------
    Minimising the mean energy is degenerate here: at T = 0.1 the alpha = 1 arm collapses
    to 0.076 bits of state entropy -- i.e. back to the argmin, 3.4540 A, the shipped
    selector, contributing exactly nothing.  At the same temperature alpha = 0.1 holds
    6.36 bits and returns 3.3414 A.  So the CVaR tail is worth +0.113 A **by preventing
    the collapse**, and that is the component's measured role.  (alpha, T) is read from a
    LEAVE-FOLD-OUT table, not chosen here.

    THE GRADIENT.  `run_cvar_vqe` optimises the EXACT parameter-shift gradient, which has
    no baseline at all, so the historical `baseline="tail"` CVaR defect -- cosine +0.6556
    at 0.758x norm against the exact gradient -- cannot arise on this path.  Where a
    SAMPLED estimator is used it must be `baseline="const"`; that is asserted below rather
    than left to a convention.
    """
    qm = core.backend("quantum")
    with clk("quantum"):
        dim = min(1 << cfg.vqe_qubits, len(top))
        if dim < (1 << cfg.vqe_qubits):
            raise ValueError(f"the hypothesis set needs 2^{cfg.vqe_qubits} = "
                             f"{1 << cfg.vqe_qubits} candidates, the pool has {len(top)}")
        o = np.asarray(top[:dim], int)
        E = _zrank(np.asarray(pool["sc"], float)[o])
        alpha, T = VQE_LFO[int(fold) % len(VQE_LFO)]
        p, cvar, H, _circ = qm.run_cvar_vqe(E, alpha, T, n=cfg.vqe_qubits,
                                            layers=cfg.vqe_layers, iters=cfg.vqe_iters,
                                            seed=cfg.vqe_seed)
        block = np.asarray(Pt, float)[:dim, :dim]
        local = consensus_medoid(block, p)
        local_unweighted = consensus_medoid(block)
    return {"o": o, "p": np.asarray(p, float), "block": block,
            "sel": int(o[local]), "sel_local": int(local),
            "sel_uniform": int(o[local_unweighted]),
            "cvar": float(cvar), "entropy_bits": float(H) / math.log(2.0),
            "alpha": float(alpha), "T": float(T),
            "collapsed": bool(float(H) / math.log(2.0) < 0.5)}


def average_weighted(Wo, block, w, clk: Clock):
    """The p_theta-weighted coordinate average -- how the quantum state ENTERS the answer.

    The unweighted form superposes on the medoid and takes a plain mean, which is the
    uniform special case of this.  Weighting by the state's own distribution is what makes
    the component materially participate in the emitted structure rather than only in a
    selection, and the uniform arm is kept beside it as the ablation.
    """
    from s8 import consensus2 as cc
    with clk("average"):
        b = consensus_medoid(block, w)
        Sup = cc.superpose_batch(np.asarray(Wo, float), np.asarray(Wo, float)[b])
        ww = np.asarray(w, float)
        ww = ww / max(float(ww.sum()), 1e-12)
        C = np.tensordot(ww, Sup, axes=(0, 0))
    return C


def legacy_stage(pool, seq, qs, cfg: Config, clk: Clock):
    """Legacy as a LATE refiner on an already-tight set.  Deployable.

    WHY LATE, and why not anywhere else
    -----------------------------------
    Legacy's 11-term field is measured to be HARMFUL as a global objective (-0.106 A when
    it enters generation) and mediocre as a global ranker (4.103 A on the normal pool
    against the shipped 3.454).  The one regime where it has measured skill is a pool that
    is already tight: 2.228 A against 2.406 on an oracle-tight pool.  The deployable
    analogue of a tight pool is the NEIGHBOURHOOD OF THE CONSENSUS MEDOID -- built with no
    native -- and `s8/integrate.py` priced exactly that: `refine_LegacyAll_top5` reaches
    3.3243 A against the shipped 3.4540, the best of the Legacy roles it tried.

    So this stage ranks the `legacy_top` pool members nearest the consensus centre by
    `lg_all`, the shipped fitted-weight combination, and returns that argmin.  It is
    reported as its own arm and as an ablation, including if the contribution is zero.
    """
    lf = core.backend("legacy")
    geo = core.backend("geometry")
    with clk("legacy"):
        o = qs["o"]
        centre = qs["sel_local"]
        k = min(cfg.legacy_top, len(o))
        nb_local = np.argsort(qs["block"][centre], kind=cfg.tie_break)[:k]
        nb = np.asarray(o, int)[nb_local]
        PHI = np.asarray(pool["PHI"], float)[nb]
        PSI = np.asarray(pool["PSI"], float)[nb]
        BB = geo.build_backbone_batch(PHI, PSI)
        bl = lf.BatchLegacy(seq, _RepShim(PHI, PSI))
        T = bl.terms_from_coords(BB, phi=PHI, psi=PSI)
        wv = np.array([lf.FITTED_WEIGHTS.get(t, 0.0) for t in lf.TERMS], float)
        lg_all = np.asarray(T, float) @ wv
        pick = int(nb[int(np.argmin(lg_all))])
    return {"neighbourhood": [int(x) for x in nb], "lg_all": [float(x) for x in lg_all],
            "sel": pick, "k": int(k)}


class _RepShim:
    """The minimal `representation` surface `BatchLegacy` reads: `_phi`, `_psi`, `n_states`.

    `BatchLegacy.__init__` builds a per-residue x per-state Ramachandran table from it
    which the coords path never reads (phi/psi are passed explicitly), so one dummy state
    is enough.  Lifted from `s8.inband._RepShim` rather than reinvented, so the Legacy
    columns this stage produces are the ones that study measured.
    """

    def __init__(self, PHI, PSI):
        self._phi = np.asarray(PHI, float)[:1].T.copy()
        self._psi = np.asarray(PSI, float)[:1].T.copy()
        self.n_states = 1
        self.n_residues = self._phi.shape[0]


def average(Wsub, Ps, clk: Clock):
    """STAGE 3a.  Superpose the subset on its medoid and take the coordinate average.

    Identical to `s8.consensus2._avg_struct`: the medoid is the member with the smallest
    mean pairwise RMSD to the rest, and the average is taken AFTER superposition onto it.
    The result has a mean CA-CA bond of ~2.96 A and is not a peptide; stage 3b is the
    inverse that makes it one.
    """
    from s8 import consensus2 as cc
    with clk("average"):
        b = int(np.argmin(cc._medoid(np.asarray(Ps, float))))
        C = cc.superpose_batch(np.asarray(Wsub, float),
                               np.asarray(Wsub, float)[b]).mean(0)
    return C, b


def project(C, seq, fold, cfg: Config, clk: Clock):
    """STAGE 3b.  The nearest ideal-geometry chain, multi-start, with the torsion prior.

    Returns (ca, phi, psi, fit_ca).  `fit` is the same projection with the prior switched
    off -- the S8-11 incumbent, kept as the internal reference the lambda = 0 rung is
    asserted to reproduce.

    Multi-start is load-bearing and is NOT a tuning knob: the projection is degenerate --
    a CA trace admits two ideal-geometry torsion solutions at near-equal distance, one
    Ramachandran-plausible and one not -- and warm-starting cannot cross between them
    (S9-2).  Switching it off would be weakening the science, not optimising it.
    """
    pj = core.backend("project")
    with clk("projection"):
        pen = pj.make_penalty(cfg.penalty, seq, int(fold))
        #: passed PER CALL, not set as a module global, so the mode cannot leak between
        #: targets inside a worker; the legacy backend has no such parameter and is not
        #: offered one.
        kw = ({"grad": cfg.project_grad}
              if getattr(pj, "SELECTABLE_GRADIENT", False) else {})
        path = pj.lam_path(C, pen, (0.0, cfg.lam), maxiter=cfg.maxiter,
                           multi=cfg.multi_start, **kw)
        fit, arm = path[0.0], path[cfg.lam]
    return arm[0], arm[1], arm[2], fit[0]


def project_single(C, pool, gi, cfg: Config, clk: Clock):
    """The SINGLE-START projection from the medoid's own torsions -- `s9.synth.fit_w`.

    This is the construction behind `s9/synth_probe.json`'s 3.2005, and it is NOT what
    `s9/final.py` emits (that is the multi-start `lam_path[0.0]`, 3.204).  Carried so both
    published numbers are reproducible from one code path instead of one of them looking
    like a discrepancy.  It is a REPORTING arm; nothing downstream consumes it.
    """
    pj = core.backend("project")
    #: clocked as a REFERENCE ARM, never as `projection`.  It is not on the deployable
    #: path -- nothing downstream consumes it -- and charging it to the projection stage
    #: would make the optimised arm look slower than the reference at the one stage where
    #: the two are supposed to be doing identical work.
    with clk("reference_arms"):
        ph0 = np.asarray(pool["PHI"], float)[gi].copy()
        ps0 = np.asarray(pool["PSI"], float)[gi].copy()
        kw = ({"grad": cfg.project_grad}
              if getattr(pj, "SELECTABLE_GRADIENT", False) else {})
        got = pj.fit_prior(C, ph0, ps0, pen=None, lam=0.0, w=None,
                           maxiter=cfg.maxiter, **kw)
    return np.asarray(got[0], float)


def relax(seq, phi, psi, cfg: Config, clk: Clock):
    """STAGE 4.  One restrained ff14SB/GBn2 minimisation of the projected arm.

    A VALIDITY stage.  S8-12 priced it exactly: k = 10-100 removes ~10^4 kcal/mol of
    builder strain for +0.011 to +0.026 A of CA-RMSD.  It is applied to what the pipeline
    emits, and its cost in accuracy is reported, never hidden.
    """
    geo = core.backend("geometry")
    ar = core.backend("amber")
    import torsion_lib2 as tl2
    gate = _WSTATE.get("amber_gate")
    if gate is not None:
        #: ADMISSION CONTROL, not a change to the computation.  Time spent QUEUEING is
        #: clocked as `amber_wait` so it can never be mistaken for time spent minimising,
        #: and the gate is released in a `finally` so a failed target cannot wedge the run.
        with clk("amber_wait"):
            gate.acquire()
    try:
        with clk("amber"):
            return _relax_inner(seq, phi, psi, cfg, geo, ar, tl2)
    finally:
        if gate is not None:
            gate.release()


def _relax_inner(seq, phi, psi, cfg, geo, ar, tl2):
    BB = geo.build_backbone_batch(np.asarray(phi, float)[None],
                                  np.asarray(psi, float)[None])
    cd = {k: v[0] for k, v in BB.items()}
    tab = tl2.library_for(seq, cfg.torsion_window, seq)   # holds the target's own seq out
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    try:
        r = ar.refine_coords(seq, rep, cd, k_restraint=cfg.amber_k,
                             steps=cfg.amber_steps, components=True)
    except Exception as exc:                                             # noqa: BLE001
        ar.clear_cache()
        return {"err": str(exc)[:300]}
    ca1 = np.asarray(r["ca"], float)
    bb = r.get("backbone", {}) or {}
    cm = r.get("components", {}) or {}
    out = {"ca": ca1, "e0": float(r["energy_initial"]), "e1": float(r["energy"]),
           "moved": float(r["restraint_rmsd"]),
           "strain_after": (float(cm.get("bond", 0.0) + cm.get("angle", 0.0))
                            if cm else None)}
    if all(k in bb for k in ("N", "CA", "C")):
        ph1, ps1 = geo.extract_torsions(bb["N"], bb["CA"], bb["C"])
        out["phi"] = np.asarray(ph1, float)
        out["psi"] = np.asarray(ps1, float)
    ar.clear_cache()
    return out


# ============================================================ one target, end to end
def deployable_view(pool, target, fold):
    """The ONLY object stages 2-4 are allowed to see.  No native quantity is in it."""
    return {"pdb": target.pdb, "n": int(target.n), "fold": int(fold), "seq": target.seq,
            "W": pool["W"], "PHI": pool["PHI"], "PSI": pool["PSI"], "S": pool["S"],
            "sim": pool["sim"], "sc": pool["sc"]}


def run_target(target, fold, cfg: Config = PROD, clk: Clock | None = None):
    """Stages 1-4 for one target.  Returns the deployable record and the stage clock.

    No native is opened anywhere in here.  `label()` does that, afterwards.
    """
    clk = clk or Clock()
    pool = retrieve(target, fold, cfg, clk)
    pool = score(pool, target.seq, fold, target.n, cfg, clk)
    v = deployable_view(pool, target, fold)
    sub, Wsub, Ps, top, Pt = filter_pool(pool, cfg, clk)
    C, b = average(Wsub, Ps, clk)
    ca, phi, psi, fit = project(C, v["seq"], fold, cfg, clk)
    fit1 = None
    if cfg.report_single_start_fit:
        fit1 = project_single(C, pool, int(sub[b]), cfg, clk)
    rec = {"pdb": target.pdb, "n": int(target.n), "fold": int(fold), "seq": target.seq,
           "fit1_ca": fit1,
           "n_windows": pool["n_windows"], "n_candidates": pool["n_candidates"],
           "n_distinct": pool.get("n_distinct"), "n_top": pool.get("n_top"),
           "sub": np.asarray(sub, int),
           "avg_ca": np.asarray(C, float), "fit_ca": np.asarray(fit, float),
           "ca": np.asarray(ca, float), "phi": np.asarray(phi, float),
           "psi": np.asarray(psi, float)}
    if cfg.amber:
        rec["amber"] = relax(target.seq, phi, psi, cfg, clk)

    # ---- the two components `s9/final.py` does not contain, each as its own arm ----
    if cfg.quantum:
        qs = quantum_stage(pool, top, Pt, fold, cfg, clk)
        rec["quantum"] = {k: qs[k] for k in ("sel", "sel_uniform", "cvar",
                                             "entropy_bits", "alpha", "T", "collapsed")}
        #: the QUANTUM SYNTHESIS: the state's own distribution weights the coordinate
        #: average, so the component participates in the emitted structure and not only
        #: in a selection.  The uniform-weight run beside it IS the ablation.
        Wo = pool["W"][qs["o"]]
        Cq = average_weighted(Wo, qs["block"], qs["p"], clk)
        qa, qp, qz, _qfit = project(Cq, v["seq"], fold, cfg, clk)
        rec["q_avg_ca"] = np.asarray(Cq, float)
        rec["q_ca"] = np.asarray(qa, float)
        rec["q_phi"] = np.asarray(qp, float)
        rec["q_psi"] = np.asarray(qz, float)
        Cu = average_weighted(Wo, qs["block"], np.ones(len(qs["o"])), clk)
        ua, _up, _uz, _ufit = project(Cu, v["seq"], fold, cfg, clk)
        rec["u_ca"] = np.asarray(ua, float)
        if cfg.amber:
            rec["q_amber"] = relax(target.seq, qp, qz, cfg, clk)
        if cfg.legacy:
            ls = legacy_stage(pool, target.seq, qs, cfg, clk)
            rec["legacy"] = ls
    elif cfg.legacy:
        #: Legacy still needs a tight set to refine INSIDE.  Without the quantum stage
        #: that set is the classical consensus medoid's neighbourhood over the same
        #: filtered prefix -- the same construction, uniform weights.
        block = np.asarray(Pt, float)
        qs = {"o": np.asarray(top, int), "block": block,
              "sel_local": consensus_medoid(block)}
        rec["legacy"] = legacy_stage(pool, target.seq, qs, cfg, clk)
    return rec, pool, clk


def label(rec, pool, target, clk: Clock, cfg: Config = PROD):
    """THE ONLY PLACE A NATIVE IS OPENED.  Runs after every structure is final.

    Emits the four reporting arms the project quotes: the shipped distogram argmin (the
    historical pipeline's own answer), the unfiltered pool's best member (the achievability
    CEILING -- not deployable), the synthesis, and the full system after relaxation.
    """
    aud = core.backend("numerics")
    from s8.inband import sel_of
    with clk("reporting"):
        #: `rr` is computed at full precision and stored float32 by the reference, and
        #: the native itself is read back from a float32 cache wherever an RMSD is scored.
        rr = _q(aud.kabsch_rmsd_batch(pool["W64"], np.asarray(target.ca, float)), cfg)
        nat = _q(np.asarray(target.ca, float), cfg)
        sub = np.asarray(rec["sub"], int)
        out = {
            "shipped": float(sel_of(pool["sc"], rr)),
            "pool_best": float(rr.min()),
            "pool_mean": float(rr.mean()),
            "top_m_best": float(rr[sub].min()),
            "top_m_mean": float(rr[sub].mean()),
            "rmsd_avg": float(aud.kabsch_rmsd_batch(rec["avg_ca"][None], nat)[0]),
            "rmsd_fit": float(aud.kabsch_rmsd_batch(rec["fit_ca"][None], nat)[0]),
            "rmsd_arm": float(aud.kabsch_rmsd_batch(rec["ca"][None], nat)[0]),
        }
        if rec.get("fit1_ca") is not None:
            out["rmsd_fit_single"] = float(aud.kabsch_rmsd_batch(
                np.asarray(rec["fit1_ca"], float)[None], nat)[0])
        a = rec.get("amber") or {}
        if "ca" in a:
            out["rmsd_full"] = float(aud.kabsch_rmsd_batch(a["ca"][None], nat)[0])
            out["d_rmsd"] = out["rmsd_full"] - out["rmsd_arm"]
        else:
            out["rmsd_full"] = None
            out["d_rmsd"] = None
        # ---- the mandated components, each as its own measured arm ----
        q = rec.get("quantum")
        if q:
            out["rmsd_vqe_sel"] = float(rr[q["sel"]])
            out["rmsd_vqe_sel_uniform"] = float(rr[q["sel_uniform"]])
            out["vqe_entropy_bits"] = q["entropy_bits"]
            out["vqe_alpha"] = q["alpha"]
            out["vqe_T"] = q["T"]
            out["vqe_collapsed"] = q["collapsed"]
            out["rmsd_q_avg"] = float(aud.kabsch_rmsd_batch(
                rec["q_avg_ca"][None], nat)[0])
            out["rmsd_q_synth"] = float(aud.kabsch_rmsd_batch(rec["q_ca"][None], nat)[0])
            out["rmsd_u_synth"] = float(aud.kabsch_rmsd_batch(rec["u_ca"][None], nat)[0])
            qa = rec.get("q_amber") or {}
            out["rmsd_system"] = (float(aud.kabsch_rmsd_batch(qa["ca"][None], nat)[0])
                                  if "ca" in qa else None)
        lg = rec.get("legacy")
        if lg:
            out["rmsd_legacy_sel"] = float(rr[lg["sel"]])
            out["legacy_k"] = lg["k"]
            #: the ablation partner: the SAME neighbourhood ranked by nothing, i.e. its
            #: own centre.  Legacy's contribution is the difference between these two.
            out["rmsd_legacy_centre"] = float(rr[lg["neighbourhood"][0]])
    return out


def _jsonable(rec, lab, clk, cfg, key):
    a = rec.get("amber") or {}
    out = {"pdb": rec["pdb"], "n": rec["n"], "fold": rec["fold"], "seq": rec["seq"],
           "n_windows": rec["n_windows"],
           "sub": np.asarray(rec["sub"], int).tolist(),
           "ca": rec["ca"].tolist(), "phi": rec["phi"].tolist(),
           "psi": rec["psi"].tolist(), "fit_ca": rec["fit_ca"].tolist(),
           "avg_ca": rec["avg_ca"].tolist(),
           "amber_ca": a["ca"].tolist() if "ca" in a else None,
           "amber_phi": a["phi"].tolist() if "phi" in a else None,
           "amber_psi": a["psi"].tolist() if "psi" in a else None,
           "amber_e0": a.get("e0"), "amber_e1": a.get("e1"),
           "amber_moved": a.get("moved"), "amber_strain_after": a.get("strain_after"),
           "amber_err": a.get("err"),
           "n_candidates": rec.get("n_candidates"), "n_distinct": rec.get("n_distinct"),
           "n_top": rec.get("n_top"),
           "q_ca": rec["q_ca"].tolist() if "q_ca" in rec else None,
           "q_amber_ca": (rec["q_amber"]["ca"].tolist()
                          if (rec.get("q_amber") or {}).get("ca") is not None else None),
           "quantum": rec.get("quantum"),
           "legacy": rec.get("legacy"),
           "timings": {k: round(v, 4) for k, v in clk.t.items()},
           "wall": round(clk.total(), 4),
           "cfg_key": key, "dev_mode": cfg.dev_mode}
    out.update(lab)
    return out


# ============================================================ worker plumbing
THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NT",
               "OPENMM_CPU_THREADS")


def _pin_threads(t: int):
    for var in THREAD_VARS:
        os.environ[var] = str(t)
    try:
        import torch
        torch.set_num_threads(t)
    except Exception:                                                    # noqa: BLE001
        pass


class _pinned_env:
    """Pin the thread-count environment for the lifetime of a worker pool, then restore."""

    def __init__(self, t):
        self.t = int(t)
        self.saved = {}

    def __enter__(self):
        self.saved = {v: os.environ.get(v) for v in THREAD_VARS}
        for v in THREAD_VARS:
            os.environ[v] = str(self.t)
        return self

    def __exit__(self, *exc):
        for v, old in self.saved.items():
            if old is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = old
        return False


_WSTATE = {}


def pin_affinity(core: int) -> bool:
    """Pin this process to one logical core, via `SetProcessAffinityMask`.

    Used only when `Config.affinity` is on.  Windows has no `sched_setaffinity`; this is
    the equivalent, and it is a PROCESS mask, so it must be set in the worker itself.
    """
    import ctypes
    k32 = ctypes.WinDLL("kernel32")
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    k32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    return bool(k32.SetProcessAffinityMask(k32.GetCurrentProcess(),
                                           1 << (core % (os.cpu_count() or 1))))


def _worker_init(cfg_dict, backends_env, seqs, slot=None, amber_gate=None):
    """Runs once per worker process, BEFORE any target.

    Pins the thread counts (the oversubscription trap: W workers x T threads on 8 cores),
    optionally pins this worker to one core, bounds the library memos so W workers fit in
    RAM, forces the requested backend set, and arms the ESM guard with the manifest's own
    sequences so a missing hot-cache entry is an error rather than a silent zero vector.
    """
    if backends_env is not None:
        os.environ["CORE_BACKENDS"] = backends_env
    cfg = Config(**cfg_dict)
    _pin_threads(cfg.threads_per_worker)
    if cfg.affinity and slot is not None:
        with slot.get_lock():
            k = int(slot.value)
            slot.value += 1
        core_id = (cfg.affinity_cores[k % len(cfg.affinity_cores)]
                   if cfg.affinity_cores else k * cfg.affinity_stride)
        _WSTATE["core"] = int(core_id)
        _WSTATE["affinity_ok"] = pin_affinity(int(core_id))
    _LIB.n_folds = cfg.n_folds
    _LIB.cap = cfg.window_cache
    _LIB.cap_members = cfg.member_cache
    _WSTATE["amber_gate"] = amber_gate
    _WSTATE["cfg"] = cfg
    _WSTATE["seqs"] = set(seqs)
    _WSTATE["t0"] = time.perf_counter()


def _worker_run(args):
    """One target inside a worker.  Returns (pdb, record | None, error | None)."""
    pdb, cache_dir, key = args
    cfg = _WSTATE["cfg"]
    db = core.backend("data")
    guard_esm(_WSTATE["seqs"])
    path = os.path.join(cache_dir, f"{pdb}.json")
    clk = Clock()
    with clk("io"):
        done = _read_json(path)
    if done is not None and done.get("cfg_key") == key:
        done["resumed"] = True
        return pdb, done, None
    if not wait_for_memory(cfg.stop_pct, where=pdb):
        return pdb, None, f"memory gate: {mem_pct()}% >= {cfg.stop_pct}% for 300s"
    if "targets" not in _WSTATE:
        _WSTATE["targets"] = {p.pdb: p for p in db.load()}
    target = _WSTATE["targets"][pdb]
    fold = db.folds(cfg.n_folds)[target.seq]
    try:
        rec, pool, clk = run_target(target, fold, cfg, clk)
        lab = label(rec, pool, target, clk, cfg)
    except Exception as exc:                                             # noqa: BLE001
        import traceback
        return pdb, None, traceback.format_exc()[-1500:]
    with clk("serialise"):
        blob = _jsonable(rec, lab, clk, cfg, key)
    with clk("io"):
        _atomic_write_json(path, blob)
    blob["timings"] = {k: round(v, 4) for k, v in clk.t.items()}
    blob["wall"] = round(clk.total(), 4)
    blob["resumed"] = False
    blob["worker"] = _worker_footprint()
    return pdb, blob, None


def _worker_footprint():
    """What this worker is holding right now -- measured, per target, in every run."""
    r = proc_rss() or {}
    out = {"pid": os.getpid(), "core": _WSTATE.get("core"),
           "affinity_ok": _WSTATE.get("affinity_ok"),
           **{k: round(v, 1) for k, v in r.items()}}
    out.update({k: (round(v, 1) if isinstance(v, float) else v)
                for k, v in _LIB.nbytes().items()})
    out["lib_cache"] = dict(_LIB.stats)
    return out


def dispatch_units(tg, folds):
    """The work units, and their order.  A SCHEDULING decision: it cannot move a number.

    THE TWO CONSTRAINTS PULL AGAINST EACH OTHER
    -------------------------------------------
    The library memo is keyed on (fold, length), so a worker wants CONSECUTIVE targets
    that share that key.  Load balance wants the opposite: targets sorted by (fold,
    length) are also sorted by cost, because cost climbs steeply with length -- measured,
    10.5 s at n=9 against 24.8 s at n=16 -- so any contiguous equal-COUNT chunking hands
    one worker all the short targets and another all the long ones.

    That is exactly what the previous chunked dispatch did, and it cost 16.5% of the
    machine.  Measured on the 8-worker run: workers were busy 351.4 s at most and 232.3 s
    at least, a 119.1 s spread, against a 367.1 s wall.  Occupancy 6.68 of 8 cores.

    THE FIX SATISFIES BOTH.  The unit is the (fold, length) GROUP itself -- every target
    in a unit shares the memo key, so locality is perfect by construction rather than
    approximate -- and the units are dispatched DYNAMICALLY, longest first.  Longest-first
    on a dynamic queue is the classic LPT schedule; with 39 units over 8 workers it lands
    within 2% of perfect balance.

    Simulated against the measured per-target costs of the 8-worker run:

        schedule                                 makespan   spread   memo hits
        contiguous equal-count chunks (before)     351.4 s   119.1 s     55%
        (fold, n) groups, fold order, no LPT       362.9 s   108.2 s     69%
        (fold, n) groups, LPT by measured cost      308.2 s     6.3 s     69%
        (fold, n) groups, LPT by total length       311.2 s    10.0 s     69%
        perfect balance at the same CPU-work        306.4 s        --      --

    Note the second row: grouping ALONE makes it worse.  The ordering is what does the
    work, and grouping is what keeps the memo.  Note also the fourth: the WEIGHT barely
    matters -- even ordering by target count reaches 312.2 s -- so the weight here is
    total sequence length, a structural property of the manifest.  Weighting by measured
    cost would buy 3 s and would make the scheduler depend on a results artefact, which
    is a far worse trade than it looks.

    THE ONE PLACE THIS SCHEDULE COSTS SOMETHING -- READ BEFORE RAISING THE WORKER COUNT
    -----------------------------------------------------------------------------------
    Longest-first orders by sequence length, and length is also what sets an AMBER
    system's atom count.  So the first W units dispatched are the W BIGGEST systems, and
    they are resident simultaneously.  Peak memory is therefore front-loaded by
    construction, and it is a property of the SCHEDULE COMBINED WITH THE ARM rather than
    of either alone:

        like-for-like (one minimisation per target)   w8 fits: 92% peak, 0 declines
        four-component (TWO per target)               w8 does NOT: AMBER declined on 7
                                                      targets and the full system on 26

    A decline is recorded per target, not as a run error, so that second row still
    reports 126/126 while `rmsd_full` is a mean over 119 and `rmsd_system` over 100 --
    an arm's mean moving without any single number being wrong.  `core.bench` reports
    `n_amber_declined` and refuses such a run as a headline, but the fix here is simply
    to run the heavier arm at w6.  Raising the worker count for an arm that relaxes more
    than once per target needs this checked, not assumed.
    """
    groups = {}
    for p in tg:
        groups.setdefault((folds[p.seq], int(p.n)), []).append(p)
    return sorted(groups.values(), key=lambda g: -sum(int(q.n) for q in g))


def _worker_run_group(tasks):
    """One dispatch unit: every target in it shares a (fold, length) memo key."""
    return [_worker_run(t) for t in tasks]


def _worker_targets(targets):
    _WSTATE["targets"] = {p.pdb: p for p in targets}


# ============================================================ the run
def cache_dir_for(cfg: Config, backends=None):
    d = os.path.join(CACHE_ROOT, cfg.key(backends))
    os.makedirs(d, exist_ok=True)
    return d


def run(manifest_name="tuning126", cfg: Config = PROD, workers: int = 1,
        verbose=True, backends_env=None, limit=None, order="fold"):
    """The whole manifest.  Resumable, parallel across targets, atomic per target.

    `order="fold"` groups targets by fold, so a worker that picks up consecutive tasks
    reuses one library, one window set and one distogram model instead of five.
    """
    tg = manifest(manifest_name)
    if limit:
        tg = tg[:limit]
    db = core.backend("data")
    folds = db.folds(cfg.n_folds)
    if order == "fold":
        tg = sorted(tg, key=lambda p: (folds[p.seq], p.n, p.pdb))
    else:
        tg = sorted(tg, key=lambda p: p.pdb)
    #: the cache key must describe the backends the WORKERS will use, not the parent's.
    backends = (core.backend_report() if backends_env is None
                else {"forced": backends_env})
    key = cfg.key(backends)
    cdir = cache_dir_for(cfg, backends)
    seqs = [p.seq for p in tg]

    t_start = time.perf_counter()
    rows, errs = [], []
    if workers <= 1:
        _worker_init(asdict(cfg), backends_env, seqs, None, None)
        _worker_targets(tg)
        for k, p in enumerate(tg):
            pdb, blob, err = _worker_run((p.pdb, cdir, key))
            if err:
                errs.append((pdb, err))
                if verbose:
                    print(f"  {pdb}: {err.splitlines()[-1][:120]}", flush=True)
                if "memory gate" in err:
                    break
                continue
            rows.append(blob)
            if verbose:
                _progress(k + 1, len(tg), blob)
    else:
        import concurrent.futures as cf
        ctx = _mp_context()
        #: THE OVERSUBSCRIPTION TRAP.  `_worker_init` runs AFTER the child has already
        #: imported numpy to unpickle the task, and OpenBLAS reads its thread count once,
        #: at load.  So the pin has to be in the environment the child STARTS with, which
        #: means setting it here, in the parent, before the pool exists -- `spawn` copies
        #: `os.environ` into the child.  Setting it only in the initializer looks correct
        #: and silently leaves W x 8 threads on an 8-core box.
        slot = ctx.Value("i", 0) if cfg.affinity else None
        gate = (ctx.Semaphore(cfg.amber_concurrency)
                if 0 < cfg.amber_concurrency < workers else None)
        with _pinned_env(cfg.threads_per_worker), \
                cf.ProcessPoolExecutor(max_workers=workers, mp_context=ctx,
                                       initializer=_worker_init,
                                       initargs=(asdict(cfg), backends_env, seqs,
                                                 slot, gate)) as ex:
            #: DISPATCH WHOLE (fold, length) GROUPS, LONGEST FIRST.  See `dispatch_units`.
            units = dispatch_units(tg, folds)
            futs = [ex.submit(_worker_run_group,
                              [(p.pdb, cdir, key) for p in u]) for u in units]
            k = 0
            for fu in cf.as_completed(futs):
                for pdb, blob, err in fu.result():
                    k += 1
                    if err:
                        errs.append((pdb, err))
                        if verbose:
                            print(f"  {pdb}: {err.splitlines()[-1][:120]}", flush=True)
                        continue
                    rows.append(blob)
                    if verbose:
                        _progress(k, len(tg), blob)
    wall = time.perf_counter() - t_start
    return {"manifest": manifest_name, "n_requested": len(tg), "n_done": len(rows),
            "wall": wall, "workers": workers,
            "threads_per_worker": cfg.threads_per_worker,
            "cfg": cfg.science(), "cfg_key": key, "dev_mode": cfg.dev_mode,
            "backends": backends,
            "cache_dir": cdir, "rows": rows,
            "errors": [{"pdb": p, "error": e} for p, e in errs],
            "stage_totals": stage_totals(rows), "summary": summarise(rows)}


def _mp_context():
    import multiprocessing as mp
    return mp.get_context("spawn")


def _progress(k, n, blob):
    full = blob.get("rmsd_full")
    print(f"[{k:3d}/{n}] {blob['pdb']:6} n={blob['n']:2d} "
          f"synth={blob['rmsd_arm']:.3f} "
          f"full={'  n/a' if full is None else format(full, '.3f')} "
          f"pool={blob['pool_best']:.3f} shipped={blob['shipped']:.3f} "
          f"({blob['wall']:.1f}s{' cached' if blob.get('resumed') else ''}) "
          f"[mem {mem_pct()}%]", flush=True)


def stage_totals(rows):
    """Summed per-stage seconds over the manifest, and separately over what ran now.

    A per-stage number is a PER-TARGET measurement, taken when that target was computed,
    and it stays valid in the checkpoint.  A wall clock is a property of one invocation
    and does not.  So a run that had to be resumed -- and on a shared box, most of them
    are -- can still report an honest per-stage table over all 126 targets while being
    disqualified, by `cold`, from reporting a wall-clock speedup.  Both totals are
    emitted so the two can never be confused:

        `<stage>`  summed over EVERY target in the manifest (measured when computed)
        `_fresh_*` summed over only the targets this invocation actually computed
    """
    fresh = [r for r in rows if not r.get("resumed")]
    out = {k: float(sum(r.get("timings", {}).get(k, 0.0) for r in rows)) for k in STAGES}
    out["_sum"] = float(sum(out[k] for k in STAGES))
    for k in STAGES:
        out[f"_fresh_{k}"] = float(sum(r["timings"].get(k, 0.0) for r in fresh))
    out["_fresh_sum"] = float(sum(out[f"_fresh_{k}"] for k in STAGES))
    out["_n_fresh"] = len(fresh)
    out["_n_resumed"] = len(rows) - len(fresh)
    return out


def summarise(rows):
    """The reported arms.  Means over the targets that completed, nothing dropped.

    REDUCED IN CANONICAL PDB ORDER, always.  Every per-target number is bit-identical
    across worker counts -- 17 fields x 126 targets, all exact -- but a mean accumulated
    in COMPLETION order is not: four workers finish in a different order every run, and
    floating-point addition is not associative, so seven aggregate statistics moved by up
    to 3.3e-16 relative (1-2 ulps) between a 1-worker and a 4-worker run.  That is a
    reduction artefact and it read as a scientific difference.  Sorting here removes it at
    the source, so `parallelism changes nothing` is exact rather than exact-to-a-tolerance.
    """
    if not rows:
        return {}
    rows = sorted(rows, key=lambda r: r["pdb"])
    def col(k):
        v = [r[k] for r in rows if r.get(k) is not None]
        return np.asarray(v, float)
    out = {"n": len(rows)}
    for k in ARMS:
        v = col(k)
        if len(v):
            out[k] = {"n": int(len(v)), "mean": float(v.mean()),
                      "median": float(np.median(v)),
                      "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                      "min": float(v.min()), "max": float(v.max()),
                      "frac_under_2.0": float((v < 2.0).mean())}
    # the paired gain over the shipped pipeline -- the property of the ARCHITECTURE
    paired = {}
    for arm in ("rmsd_arm", "rmsd_full", "rmsd_vqe_sel", "rmsd_q_synth", "rmsd_system",
                "rmsd_legacy_sel"):
        pair = [(r[arm], r["shipped"]) for r in rows if r.get(arm) is not None]
        if len(pair) > 1:
            paired[f"{arm}_vs_shipped"] = _paired([a for a, _ in pair],
                                                  [b for _, b in pair])
    out["paired"] = paired
    out["ablations"] = ablations(rows)
    ent = col("vqe_entropy_bits")
    if len(ent):
        out["vqe_entropy_bits"] = {"mean": float(ent.mean()), "min": float(ent.min()),
                                   "n_collapsed": int(sum(bool(r.get("vqe_collapsed"))
                                                          for r in rows))}
    dup = [(r["n_top"], r["n_distinct"]) for r in rows
           if r.get("n_top") and r.get("n_distinct")]
    if dup:
        out["duplicate_candidates"] = {
            "n_top_mean": float(np.mean([a for a, _ in dup])),
            "n_distinct_mean": float(np.mean([b for _, b in dup])),
            "frac_duplicate": float(1.0 - np.mean([b / a for a, b in dup]))}
    return out


#: Every arm the harness reports.  The first block is the like-for-like `s9/final.py`
#: pipeline; the second is the two mandated components that file does not contain.
ARMS = ("shipped", "pool_best", "top_m_best", "rmsd_avg", "rmsd_fit",
        "rmsd_fit_single", "rmsd_arm", "rmsd_full",
        "rmsd_vqe_sel", "rmsd_vqe_sel_uniform", "rmsd_q_avg", "rmsd_q_synth",
        "rmsd_u_synth", "rmsd_system", "rmsd_legacy_sel", "rmsd_legacy_centre")


def _paired(a, b):
    """Paired difference a - b with a 95% CI and win/loss counts."""
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[np.isfinite(d)]
    if len(d) < 2:
        return None
    se = float(d.std(ddof=1) / math.sqrt(len(d)))
    return {"mean_diff": float(d.mean()), "se": se,
            "ci95": [float(d.mean() - 1.96 * se), float(d.mean() + 1.96 * se)],
            "n_better": int((d < -1e-9).sum()), "n_worse": int((d > 1e-9).sum()),
            "n": int(len(d))}


#: THE COMPONENT ABLATIONS.  Each entry is (name, arm-with, arm-without, what it prices).
#: A zero here is a RESULT, not a failure: the project's own findings say VQE ties uniform
#: sampling at 4^n configurations and Legacy does not earn its place on accuracy.  The
#: mandate is that the components genuinely participate and are honestly priced, and an
#: ablation that comes out at zero does exactly that.
ABLATIONS = (
    ("cvar_vs_collapse", "rmsd_vqe_sel", "shipped",
     "CVaR-VQE selection against the shipped argmin it collapses to at alpha=1, T=0.1"),
    ("vqe_weights_vs_uniform", "rmsd_vqe_sel", "rmsd_vqe_sel_uniform",
     "the state's distribution against a uniform ensemble over the same hypotheses"),
    ("quantum_synthesis_vs_uniform", "rmsd_q_synth", "rmsd_u_synth",
     "the quantum-weighted coordinate average against the uniform one, same 2^n set"),
    ("quantum_synthesis_vs_classical", "rmsd_q_synth", "rmsd_arm",
     "the quantum synthesis against the s9 classical synthesis"),
    ("amber_validity_cost", "rmsd_full", "rmsd_arm",
     "what the ff14SB/GBn2 validity stage costs in accuracy (expected POSITIVE)"),
    ("legacy_refinement", "rmsd_legacy_sel", "rmsd_legacy_centre",
     "Legacy ranking inside the consensus neighbourhood against that neighbourhood's "
     "own centre"),
    ("full_system_vs_s9", "rmsd_system", "rmsd_full",
     "the four-component system against the like-for-like s9/final.py pipeline"),
)


def ablations(rows):
    out = {}
    rows = sorted(rows, key=lambda r: r["pdb"])     # canonical reduction order, as above
    for name, with_arm, without_arm, what in ABLATIONS:
        pair = [(r[with_arm], r[without_arm]) for r in rows
                if r.get(with_arm) is not None and r.get(without_arm) is not None]
        if len(pair) < 2:
            continue
        p = _paired([a for a, _ in pair], [b for _, b in pair])
        p["arms"] = [with_arm, without_arm]
        p["prices"] = what
        p["mean_with"] = float(np.mean([a for a, _ in pair]))
        p["mean_without"] = float(np.mean([b for _, b in pair]))
        p["contribution_is_zero_within_ci"] = bool(p["ci95"][0] <= 0.0 <= p["ci95"][1])
        out[name] = p
    return out


# ============================================================ reading a finished cache
def load_cached(cfg: Config = PROD, backends=None, manifest_name=None):
    """Every per-target checkpoint under this config's key."""
    d = os.path.join(CACHE_ROOT, cfg.key(backends))
    if not os.path.isdir(d):
        return []
    want = None
    if manifest_name:
        want = {p.pdb for p in manifest(manifest_name)}
    out = []
    for f in sorted(os.listdir(d)):
        if not f.endswith(".json"):
            continue
        if want is not None and f[:-5] not in want:
            continue
        r = _read_json(os.path.join(d, f))
        if r:
            out.append(r)
    return out


# ============================================================ cli
def main(argv=None):
    ap = argparse.ArgumentParser(prog="core.pipeline", description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=("run", "stats", "manifests"))
    ap.add_argument("--manifest", default="tuning126", choices=MANIFESTS)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--no-amber", action="store_true")
    ap.add_argument("--components", action="store_true",
                    help="run the full four-component system: the s9 pipeline PLUS the "
                         "CVaR-VQE hypothesis selector and the Legacy late refiner.  "
                         "Adds work the s9 path does not do, so the like-for-like "
                         "speedup must be measured WITHOUT it.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--backends", default=None,
                    help="CORE_BACKENDS override for the workers ('legacy' forces the "
                         "root modules)")
    ap.add_argument("--i-am-spending-the-benchmark", action="store_true")
    a = ap.parse_args(argv)

    if a.cmd == "manifests":
        for m in MANIFESTS:
            tg = manifest(m)
            print(f"{m:12} n={len(tg):3d}  lengths {min(p.n for p in tg)}-"
                  f"{max(p.n for p in tg)}")
        return 0
    if a.manifest == "benchmark60" and not a.i_am_spending_the_benchmark:
        print("REFUSED: the 60-target benchmark's single pre-registered pass is spent. "
              "Use tuning126 for all optimisation work.  Pass "
              "--i-am-spending-the-benchmark only for a real, pre-registered evaluation.")
        return 2
    cfg = replace(PROD, threads_per_worker=a.threads, amber=not a.no_amber,
                  quantum=a.components, legacy=a.components,
                  dev_mode=("no-amber" if a.no_amber else PROD.dev_mode))
    if a.cmd == "stats":
        rows = load_cached(cfg, manifest_name=a.manifest)
        print(json.dumps({"n": len(rows), "summary": summarise(rows)}, indent=2))
        return 0
    out = run(a.manifest, cfg, workers=a.workers, backends_env=a.backends,
              limit=a.limit)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=2,
                     default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
