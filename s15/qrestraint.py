"""SPRINT 15 -- THE QUANTUM EXPERIMENT: does Sprint 14's decisive negative survive a change
of OBJECTIVE CLASS, from selective to generative?

WHY THIS MODULE IS THE CENTRE OF THE QUANTUM STORY.

Sprint 14's headline negative was measured on **selective** objectives -- Legacy and AMBER
energies, which score a candidate for nativeness. On those, running CVaR-VQE was WORSE than
not running it: 0 wins out of 12 against best-of-N drawn from the *untrained* circuit, at
+0.65 to +1.32 A. The mechanism was traced precisely: optimisation concentrates the
distribution onto a tail in which every objective is at chance, so the argmin sits at the
41st-48th percentile of its own tail, and `tail mean - tail best` (1.93 A) accounts for the
entire loss.

That mechanism is a statement about **the objective**, not about VQE. It says: if the
objective cannot order structures within the region the optimiser concentrates on, then
concentrating is strictly harmful and any sampler that stays spread out will win. The obvious
and never-run follow-up is therefore:

    **Give the optimiser an objective whose ordering skill has been independently measured to
    be real, and re-run the identical experiment.**

`s15/distgeo.py` supplies exactly that. Restraint satisfaction is a **generative** objective:
it does not ask "is this structure native-like?" but "does this structure reproduce the
predicted distances?" -- a question with a checkable, non-chance answer. And K1-CORRECTED
measured its ceiling: from true distances the continuous fit reaches **0.611 A**, so the
objective demonstrably orders structures by accuracy when its restraints are good.

THE EXPERIMENT. One variable changes. Everything else is Sprint 14's apparatus, unmodified:
the same `Enum` index algebra, the same MPS ansatz, the same `Counter` hard budget, the same
searches, the same two-axis reporting. Only `E` changes -- from a physical energy to a
restraint residual.

    E_ls_pred     sum_{|i-j|>=2} w_ij (d_ij(x) - dhat_ij)^2 , w = 1/sd^2   NATIVE-FREE
    E_ml_pred     - sum_{|i-j|>=2} log P_ij(d_ij(x))                       NATIVE-FREE
    E_ls_pool     the same least squares against the retrieval pool's median distances
    E_combined    E_ml_pred plus the pool's log-likelihood                 NATIVE-FREE
    E_ORACLE_true the same least squares against the TRUE distance matrix  ORACLE, ceiling

The ORACLE arm is the control that separates "the optimiser cannot exploit this objective
class" from "the objective class is fine but these particular restraints are too noisy" --
a distinction Sprint 14 showed decides the interpretation of every negative result.

ARMS AT MATCHED BUDGET. Every arm gets an identical hard budget of objective evaluations,
enforced by `Counter`, which is the control Sprint 14 established as mandatory:

    untrained_bestofN   sample the ansatz at theta0 and keep the best -- THE control
    vqe_cvar_a{...}     CVaR-VQE, the arm under test
    random              uniform over the k^n space
    anneal, greedy      classical search, so a quantum win must beat classical too

WHAT WOULD CONSTITUTE A POSITIVE RESULT, declared before running. VQE must beat
`untrained_bestofN` on returned RMSD, paired across targets, with a confidence interval
excluding zero, at matched budget. Beating `random` proves nothing here: Sprint 13 measured
that a zero-information constant alpha-helix beats the random control on this instrument.

WHAT WOULD CONSTITUTE A NEGATIVE RESULT, equally declared. If VQE again fails to beat its own
untrained initialisation on a generative objective with measured ordering skill, then the
Sprint 14 negative is **not** an artefact of scoring-function pathology, and it generalises to
the objective class that this sprint was built to introduce. That is a stronger and more
publishable claim than the original, and it is the one the paper would lead with.

THE MECHANISTIC READOUT, which matters more than the arm table. For every objective we report
rho(E, RMSD) globally and **restricted to the tail the optimiser actually concentrates on**
(best 1% and best 0.1% by E). Sprint 14's negative was caused by that in-tail rho being at
chance. If the restraint objective's in-tail rho is genuinely positive and VQE still loses,
the explanation must be something other than tail-blindness, and we will have falsified our
own mechanism -- which is the outcome most worth looking for.

-----------------------------------------------------------------------------------------
QRESTRAINT AGENT'S ADDITIONS, 2026-09-05.  What the design above did not carry, and why.

1.  NINETEEN TARGETS, NOT NINE.  `s14/cache/obj_enum_<PDB>.npz` holds ten more fully
    enumerated n=10 targets on the identical schema (4^10 = 1,048,576 configurations each),
    so the experiment runs on 19 targets and 1.28e7 exactly-labelled structures.  The brief
    records that n <= 4 drawn from the enumerated nine has reversed a conclusion four times;
    nineteen paired targets is the single cheapest improvement to the result's credibility.
    `Enum` below is `s14.vqe_lib.Enum` with the file lookup widened -- nothing else changes.

2.  A BUDGET LADDER, NOT A BUDGET.  A single budget is one point on a curve and Sprint 14
    showed budget choice can flip a comparison.  Three budgets spanning 16x.

3.  THE COMPARATORS ARE COMPUTED HERE, not quoted.  Sprint 14's in-tail numbers for Legacy
    and AMBER were measured on a different tail definition.  `skill_table` computes the
    identical `ordering_skill` for `legacy`, `prior` and `amber_total` on the same targets
    and the same tails, so the "is the generative objective's in-tail skill better?" question
    is answered on ONE definition rather than across two.  The AMBER rows obey the Phase 0
    binding rule (`amber_kind == 0 AND amber_idx != snap_index`) and print their per-target n.

4.  A RANDOM-TAIL NULL.  The brief records that in-tail rank statistics must be read against
    a matched null, never against 0.  `ordering_skill` returns `rho_tail_*_null_sd`, the
    standard deviation of the same statistic over tails of the same size drawn at random,
    so an in-tail rho can be called at-chance or not.

5.  CONCENTRATION IS MEASURED, NOT ASSUMED.  The whole Sprint 14 mechanism presumes the
    optimiser concentrates onto the objective's low-E tail.  Every VQE run now records the
    E-percentile of what it sampled in its first and last quarter (`conc_*`), so "did the
    premise even hold" is a measurement rather than an inference.

6.  THE CVaR DEFECTS ARE ADDRESSED EXPLICITLY.  `R.run` defaults to `baseline="const"`, the
    CORRECTED estimator, so the recorded gradient-baseline defect is NOT active in the main
    result.  `alpha_sweep` additionally runs the defective `baseline="tail"` arm side by side
    at matched everything, which is stronger evidence than merely asserting the fix was used.

Run:
    python -m s15.qrestraint                  # main ladder, 19 targets, 3 budgets
    python -m s15.qrestraint skill            # the mechanistic readout only (cheap)
    python -m s15.qrestraint alpha            # the CVaR alpha sweep + baseline defect arm
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

from s12 import instrument as I              # noqa: E402
from s14 import vqe_lib as V                 # noqa: E402
from s14 import vqe_run as R                 # noqa: E402
from s14 import obj_enum as OE               # noqa: E402
from s15 import distml as M                  # noqa: E402
from s15 import pooldist as P                # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
CACHE = os.path.join(ROOT, "s15", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

OBJECTIVES = ("E_ls_pred", "E_ml_pred", "E_ls_pool", "E_combined", "E_ORACLE_true")
#: Sprint 14's SELECTIVE objectives, recomputed here on the identical targets, the identical
#: tails and the identical statistics.  `S14_disto_bayes` is the sharpest of the four: it is
#: the SAME distogram as the generative arms, consumed selectively.
COMPARATORS = ("S14_disto_bayes", "legacy", "prior", "amber_total")
BATCH = 16384

ENUM9 = tuple(V.ENUM_TARGETS)                                    # n = 9,  4^9  configs
ENUM10 = ("1N9U", "1TOR", "2BAO", "2MD2", "2MJQ",
          "5V5B", "6B9K", "7T3H", "7VI4", "8HVS")                # n = 10, 4^10 configs
ENUM19 = ENUM9 + ENUM10


# ======================================================= the resource gate
def cpu_pct():
    """Machine-wide CPU load, via the same WMI counter the coordinator's hold is stated in."""
    import subprocess
    try:
        r = subprocess.run(
            ["powershell", "-NoProfile", "-Command",
             "(Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage "
             "-Average).Average"],
            capture_output=True, text=True, timeout=60)
        return float(r.stdout.strip())
    except Exception:
        return float("nan")


def wait_for_cpu(max_pct=90.0, tag="", poll=30, max_wait=7200):
    """Block until the box has headroom. The coordinator's hold, enforced in code.

    `s14.vqe_lib.wait_for_memory` gates RAM; nothing gated CPU, and the Sprint 14 record is
    that agent processes were measured at 3-5% of one core when three heavy jobs overlapped.
    Returns the observed load; gives up after `max_wait` seconds and says so rather than
    blocking a whole session, because a permanently saturated box is a condition to report,
    not to wait out.
    """
    t0 = time.time()
    while True:
        c = cpu_pct()
        if not np.isfinite(c) or c < max_pct:
            return c
        if time.time() - t0 > max_wait:
            print(f"  [{tag}] CPU still {c:.0f}% after {max_wait}s; proceeding and "
                  f"RECORDING that this run was made on a contended box", flush=True)
            return c
        print(f"  [{tag}] CPU {c:.0f}% >= {max_pct}; waiting", flush=True)
        time.sleep(poll)


# ======================================================= the enumerated space, widened
class Enum(V.Enum):
    """`s14.vqe_lib.Enum` with the file lookup widened to `s14/cache/obj_enum_*.npz`.

    The schema is identical (`s14/obj_enum.py` says so and the audit below checks it), so
    nothing about the index algebra, the labels or the live-qubit accounting changes.  Only
    the ten n=10 targets become reachable.
    """

    def __init__(self, pdb: str):
        path = OE.enum_path(pdb)
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        d = np.load(path)
        self.pdb = pdb
        self.path = path
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
        self.amber_kind = np.asarray(d["amber_kind"], np.int64)
        assert self.rmsd.size == self.N
        self._pw = self.k ** np.arange(self.n - 1, -1, -1)
        assert int(self.snap_states @ self._pw) == self.snap_index


# ============================================================ the enumerated CA traces
def all_ca(z, batch=BATCH):
    """(N, n, 3) CA traces for every configuration of an enumerated target.

    Built with `I.build_ca`, the bit-exact builder the projection uses, so these live on the
    same manifold as `Enum.rmsd`. `_check_consistency` asserts that they do.
    """
    N, n = z.N, z.n
    out = np.empty((N, n, 3), np.float32)
    ar = np.arange(n)
    for s in range(0, N, batch):
        e = min(s + batch, N)
        S = z.states(np.arange(s, e))                       # (B, n)
        phi = z.PHI[ar[None, :], S]
        psi = z.PSI[ar[None, :], S]
        out[s:e] = I.build_ca(phi, psi).astype(np.float32)
    return out


def _check_consistency(z, CA, nat, n_check=512, tol=1e-6):
    """Our rebuilt traces must reproduce `Enum.rmsd` exactly, or E and rmsd index different
    structures -- the exact defect that cost Sprint 12 a whole table.

    MEASURED PROVENANCE OF THE RESIDUAL, 2026-09-05.  The check passes at 1e-6 with a
    residual of about 4e-7, and the residual is NOT ours: the cached `rmsd` label is stored
    as **float32** (`s14/obj_enum.py` writes `np.empty(B, np.float32)`), whose spacing at
    3 A is 2.4e-7.  Rebuilding the same 512 configurations entirely in float64 gives
    3.60e-07 against the cache and storing our own CA in float32 raises that only to
    3.86e-07, so 93% of the residual is the label's own dtype and 7% is our storage.  Both
    are returned so the tolerance never has to be argued about again -- and neither is
    weakened to make anything pass.
    """
    rng = np.random.default_rng(0)
    idx = rng.integers(0, z.N, n_check)
    got = np.asarray([I.ca_rmsd(np.asarray(CA[i], float), nat) for i in idx])
    err = float(np.abs(got - z.rmsd[idx]).max())
    # the float64 reference: same configurations, no float32 storage anywhere
    ar = np.arange(z.n)
    S = z.states(idx)
    got64 = I.kabsch_rmsd_batch(I.build_ca(z.PHI[ar[None, :], S], z.PSI[ar[None, :], S]), nat)
    err64 = float(np.abs(got64 - z.rmsd[idx]).max())
    assert err < tol, f"{z.pdb}: rebuilt CA disagrees with Enum.rmsd by {err:.3g} A"
    return {"err_f32_store": err, "err_f64_rebuild": err64,
            "label_dtype_spacing_at_3A": float(np.spacing(np.float32(3.0)))}


# ==================================================================== the objectives
def _pair_d(CA, i, j, batch=BATCH):
    """(N, npairs) pair distances, in float32, computed in batches."""
    N = CA.shape[0]
    out = np.empty((N, i.size), np.float32)
    for s in range(0, N, batch):
        e = min(s + batch, N)
        c = np.asarray(CA[s:e], np.float32)
        out[s:e] = np.sqrt(((c[:, i] - c[:, j]) ** 2).sum(-1))
    return out


def build_objectives(z, CA):
    """Every objective as a length-N vector. Only `E_ORACLE_true` reads native coordinates."""
    u = I.load_univ(z.pdb)
    nat = np.asarray(u["nat_ca"], float)
    n = z.n
    i, j = I.pair_index(n)
    D = _pair_d(CA, i, j)

    dg = I.distogram(z.pdb, z.seq, z.fold)
    assert np.array_equal(i, np.asarray(dg["i"])) and np.array_equal(j, np.asarray(dg["j"])), \
        f"{z.pdb}: distogram pair order differs from I.pair_index"
    dhat = np.asarray(dg["expected"], float)
    sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
    centres = np.asarray(dg["centres"], float)
    w = (1.0 / sd ** 2).astype(np.float32)

    Dm, _sim, _i, _j = P.pool_distances(z.pdb, n)
    d_pool = np.median(Dm, axis=0)
    w_pool = (1.0 / np.maximum(Dm.std(axis=0), 0.25) ** 2).astype(np.float32)

    tab_dg = M.LogPTable(np.asarray(dg["prob"], float), centres)
    tab_pool = P.hist_table(Dm, centres)
    tab_comb = P.SumLogP(tab_dg, tab_pool)

    dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))

    #: SPRINT 14's OWN DISTANCE OBJECTIVE, recomputed on this space.  `s14/hamil.py`'s
    #: `e_disto` -- the `w=1` column whose in-tail pairwise accuracy was **0.390**, the worst
    #: anti-ranking Sprint 14 measured -- is `I.shipped_score`, the distogram's Bayes-RISK
    #: score.  It reads the SAME distogram as `E_ls_pred` and `E_ml_pred` and differs only in
    #: how it CONSUMES it: a per-pair expected-loss lookup that scores a candidate for
    #: nativeness (SELECTIVE) rather than a residual against the predicted distances
    #: (GENERATIVE).  Carrying it as a comparator is what makes "generative vs selective" an
    #: experiment rather than a relabelling: same information, same targets, same tails.
    def bayes():
        out = np.empty(D.shape[0])
        for s in range(0, D.shape[0], BATCH):
            e = min(s + BATCH, D.shape[0])
            out[s:e] = I.shipped_score(dg, np.asarray(D[s:e], float))
        return out

    def ls(target, weight):
        r = D - np.asarray(target, np.float32)[None, :]
        return (r * r * weight[None, :]).sum(1).astype(np.float64)

    def nll(tab):
        #: `rowsum` is the vectorised per-configuration log-likelihood; negate for an energy.
        out = np.empty(D.shape[0])
        for s in range(0, D.shape[0], BATCH):
            e = min(s + BATCH, D.shape[0])
            out[s:e] = -tab.rowsum(np.asarray(D[s:e], float))
        return out

    return {
        "E_ls_pred": ls(dhat, w),
        "E_ml_pred": nll(tab_dg),
        "E_ls_pool": ls(d_pool, w_pool),
        "E_combined": nll(tab_comb),
        "E_ORACLE_true": ls(dtrue, np.ones(i.size, np.float32)),
    }, {"S14_disto_bayes": bayes()}


# ============================================================ the mechanistic readout
def ordering_skill(E, rmsd, tails=(0.01, 0.001), n_null=60, seed=0):
    """rho(E, RMSD) globally and inside the tail the optimiser concentrates on.

    Sprint 14's whole negative was that in-tail rho sits at chance. This is the number that
    says whether a generative objective is different in kind.

    THE NULL, added by QRESTRAINT.  The brief's standing rule is that an in-tail rank
    statistic is read against a matched null, never against zero.  `rho_tail_*_null_sd` is
    the sd of the identical statistic over `n_null` tails of the SAME SIZE drawn uniformly
    at random from the same space, so `rho / null_sd` is a z-score a reader can act on.
    (The random-tail null for rho is centred on 0 by construction -- it is a permutation of
    which rows enter, not of the labels within them -- and the mean is returned so that is
    checkable rather than asserted.)
    """
    E = np.asarray(E, float)
    rmsd = np.asarray(rmsd, float)
    rng = np.random.default_rng(seed)
    out = {"rho_global": float(V.spearman(E, rmsd)), "n": int(E.size)}
    order = np.argsort(E)
    for t in tails:
        m = max(32, int(t * E.size))
        sub = order[:m]
        nulls = np.asarray([V.spearman(E[k], rmsd[k])
                            for k in (rng.integers(0, E.size, (n_null, m)))])
        out[f"rho_tail_{t}"] = float(V.spearman(E[sub], rmsd[sub]))
        out[f"rho_tail_{t}_null_mean"] = float(np.nanmean(nulls))
        out[f"rho_tail_{t}_null_sd"] = float(np.nanstd(nulls))
        out[f"tail_{t}_n"] = int(m)
        out[f"tail_{t}_mean_rmsd"] = float(rmsd[sub].mean())
        out[f"tail_{t}_best_rmsd"] = float(rmsd[sub].min())
        #: the Sprint 14 selection gap: what concentrating costs if ordering is at chance
        out[f"tail_{t}_gap"] = float(rmsd[sub].mean() - rmsd[sub].min())
        #: where the tail's argmin sits INSIDE its own tail -- Sprint 14's 41st-48th pct
        b = int(np.argmin(E[sub]))
        out[f"tail_{t}_argmin_pct_in_tail"] = float((rmsd[sub] < rmsd[sub][b]).mean())
    k = int(np.argmin(E))
    out["rmsd_of_certified_optimum"] = float(rmsd[k])
    out["rmsd_best_in_space"] = float(rmsd.min())
    out["rmsd_mean_in_space"] = float(rmsd.mean())
    out["percentile_of_native_best"] = float((E < E[int(np.argmin(rmsd))]).mean())
    return out


#: SPRINT 14's OWN STATISTIC, so the answer to question (c) does not rest on a change of
#: metric.  Sprint 14 measured TAIL-RESTRICTED PAIRWISE ACCURACY at gap > 0.25 A and read it
#: against a RANDOM-TAIL NULL of 0.524-0.527, never against 0.500 (`s14/vqe_FINDINGS.md` 8a,
#: `s14/ener_lib.tail_accuracy`).  Recomputing it here means the generative objectives and
#: the selective energies are compared on TWO independent statistics -- Spearman rho over the
#: tail, and Sprint 14's pairwise accuracy -- rather than on one.
S14_FRACS = (1.0, 0.1, 0.01, 0.001)
S14_GAP = 0.25


def s14_tail_accuracy(E, rmsd, fracs=S14_FRACS, n_null=20, seed=0):
    """`s14.ener_lib.tail_accuracy` at Sprint 14's own settings, plus its random-tail null.

    The null is the SAME statistic on a tail of the same size drawn uniformly at random,
    which Sprint 14 measured at 0.524-0.527 rather than 0.500 because gaps inside a tight
    tail are small by construction.  It is recomputed per objective rather than quoted,
    because the null depends on the RMSD distribution of the space, which differs by target.
    """
    from s14 import ener_lib as EL
    E = np.asarray(E, float); rmsd = np.asarray(rmsd, float)
    rng = np.random.default_rng(seed)
    out = {}
    for q in fracs:
        r = EL.tail_accuracy(E, rmsd, q=q, rng=np.random.default_rng(seed),
                             gap_band=(S14_GAP, 1e9))
        nulls = [r["acc_gap_matched"]]        #: at q=1.0 the random tail IS the whole space
        m = max(int(round(q * E.size)), 2)
        if q < 1.0:
            nulls = []
            for _ in range(n_null):
                sel = rng.choice(E.size, size=m, replace=False)
                nulls.append(EL.tail_accuracy(E[sel], rmsd[sel], q=1.0,
                                              rng=np.random.default_rng(seed),
                                              gap_band=(S14_GAP, 1e9))["acc_gap_matched"])
        out[str(q)] = {"acc": r["acc"], "acc_gap_matched": r["acc_gap_matched"],
                       "n_tail": r["n_tail"], "enough": r["enough"],
                       "tail_mean_rmsd": r["tail_mean_rmsd"],
                       "tail_best_rmsd": r["tail_best_rmsd"],
                       "null_mean": float(np.nanmean(nulls)),
                       "null_sd": float(np.nanstd(nulls))}
    b, t = out[str(fracs[0])]["acc_gap_matched"], out[str(fracs[-1])]["acc_gap_matched"]
    out["tail_minus_bulk"] = float(t - b)
    out["selection_decomposition_q0.01"] = EL.selection_decomposition(E, rmsd, q=0.01)
    return out


def comparator_skill(z, extras=None, tails=(0.01, 0.001)):
    """The SELECTIVE energies, on the identical targets and the identical tail definition.

    Sprint 14's Legacy and AMBER in-tail numbers were computed on a different instrument
    (in-band pair accuracy, a top-100 band).  Quoting them beside a Spearman rho over the
    best 1% would be comparing two statistics, so both are recomputed here.

    AMBER obeys the Phase 0 BINDING RULE: `amber_kind == 0 AND amber_idx != snap_index`.
    That leaves only the unbiased stratum minus the force-included ORACLE snap, and the
    surviving n is returned with every row because the tail of a ~1,100-row subsample is
    ~11 structures and a reader must be able to see that.
    """
    out = {}
    for nm, E in (("legacy", z.legacy), ("prior", z.prior),
                  *(sorted((extras or {}).items()))):
        out[nm] = ordering_skill(E, z.rmsd, tails)
        out[nm]["stratum"] = "full space"
        out[nm]["s14_tail_accuracy"] = s14_tail_accuracy(E, z.rmsd)
    m = (z.amber_kind == 0) & (z.amber_idx != z.snap_index)
    idx = z.amber_idx[m]
    if idx.size >= 64:
        s = ordering_skill(z.amber_total[m], z.rmsd[idx], tails)
        s["stratum"] = "amber_kind==0 AND amber_idx!=snap_index"
        s["n_rows"] = int(idx.size)
        s["n_rows_before_mask"] = int(z.amber_idx.size)
        #: the tail of a ~1,100-row subsample is ~11 structures, below `tail_accuracy`'s own
        #: `min_tail=30` guard, so the 0.001 row comes back NaN by design rather than by
        #: accident.  It is left in so the reader sees WHY AMBER has no 0.1% row here.
        s["s14_tail_accuracy"] = s14_tail_accuracy(z.amber_total[m], z.rmsd[idx])
        out["amber_total"] = s
    return out


# ==================================================================== the arms
_BON_CACHE: dict = {}


def _untrained_indices(n_qubits, budget, seed, ansatz="mps2f", scale=0.8, shots=512):
    """The indices the UNTRAINED circuit draws. A function of (n_qubits, budget, seed) only.

    Cached because it does not depend on E, so the identical control arm is reused across
    the five objectives instead of being redrawn five times.  This is not an approximation:
    the draw is deterministic given the seed, so the cached array IS what the loop produced.
    """
    key = (int(n_qubits), int(budget), int(seed), ansatz, float(scale), int(shots))
    if key in _BON_CACHE:
        return _BON_CACHE[key]
    rng = np.random.default_rng(seed)
    an = R.make_ansatz(ansatz, n_qubits)
    th0 = R.init_theta(an, np.random.default_rng(seed), scale, "random")
    got, left = [], int(budget)
    while left > 0:
        bits = an.sample(th0, min(shots, left), rng)
        j = R._bits_to_index(bits)
        got.append(j)
        left -= j.size
    idx = np.concatenate(got)[:budget]
    _BON_CACHE[key] = idx
    if len(_BON_CACHE) > 64:
        _BON_CACHE.pop(next(iter(_BON_CACHE)))
    return idx


def untrained_best_of_n(E, n_qubits, budget, ansatz="mps2f", seed=0, scale=0.8):
    """THE control. Sample the ansatz at its UNTRAINED parameters, keep the best by E.

    This is the arm Sprint 14 found no VQE run could beat, and the literature agent found no
    precedent for anywhere -- the nearest published control is uniform random sampling, which
    is strictly weaker because it does not share the ansatz's support.
    """
    c = V.Counter(E, budget)
    c(_untrained_indices(n_qubits, budget, seed, ansatz, scale))
    return c


def ordinal_pct(x):
    """Ordinal (tie-arbitrary) percentile of every entry, by argsort only.

    `V.ranks` resolves ties with a PYTHON loop over the whole array, which is 1e6 iterations
    on an n=10 target and is called once per arm.  For a percentile DIAGNOSTIC the tie rule
    is immaterial (ties are a vanishing fraction of a continuous energy), so this is the
    argsort-only form, computed once per objective and reused.  Every RANKING STATISTIC that
    a claim rests on still goes through the audited `V.spearman` / `V.ranks`.
    """
    o = np.argsort(np.asarray(x, float), kind="mergesort")
    r = np.empty(len(o), np.float64)
    r[o] = np.arange(len(o), dtype=np.float64)
    return r / max(1.0, len(o) - 1.0)


def _concentration(seen, e_pct, r_pct):
    """Did the optimiser actually concentrate onto the objective's low-E tail?

    The whole Sprint 14 mechanism presumes it does.  Reported as the mean global E-percentile
    of the first and last quarter of what was sampled, plus the same for RMSD, plus the
    distinct-configuration fraction.  A run whose last quarter is not below its first has not
    concentrated, and no argument from tail-blindness applies to it.
    """
    seen = np.asarray(seen, np.int64)
    if seen.size < 8:
        return {}
    r, rr = e_pct, r_pct
    q = max(1, seen.size // 4)
    return {"conc_e_pct_first_q": float(r[seen[:q]].mean()),
            "conc_e_pct_last_q": float(r[seen[-q:]].mean()),
            "conc_rmsd_pct_first_q": float(rr[seen[:q]].mean()),
            "conc_rmsd_pct_last_q": float(rr[seen[-q:]].mean()),
            "conc_distinct_frac": float(np.unique(seen).size / seen.size)}


def _vqe_row(E, z, alpha, budget, seed, pc, baseline="const", shots=512, ansatz="mps2f"):
    r = R.run(E, z.n_qubits, alpha, budget, shots=shots, ansatz=ansatz, seed=seed,
              rmsd=z.rmsd, bits_per_res=z.bits_per_res, exact_dist=False,
              baseline=baseline, keep_seen=True)
    exact_i = int(np.argmin(E))
    row = {"evals": int(r["evals"]), "iters": int(r["iters"]),
           "best_e": float(r["best_e"]),
           "objective_gap": float(r["best_e"] - E[exact_i]),
           "rmsd_returned": float(r["rmsd_returned"]),
           "structural_gap": float(r["rmsd_returned"] - z.rmsd.min()),
           "rmsd_of_exact_optimum": float(z.rmsd[exact_i]),
           "rmsd_best_seen": float(r["rmsd_best_seen"]),
           "rmsd_mean_seen": float(r["rmsd_mean_seen"]),
           "init_rmsd_best": float(r["init_rmsd_best"]),
           "init_rmsd_mean": float(r["init_rmsd_mean"]),
           "param_displacement": float(r["param_displacement"]),
           "diversity": float(r["diversity"]),
           "baseline": baseline, "alpha": float(alpha),
           "grad_norm_first": (float(r["grad_norms"][0]) if r["grad_norms"] else float("nan")),
           "dead_start": bool(r["grad_norms"] and r["grad_norms"][0] == 0.0)}
    row.update(_concentration(r["seen"], *pc))
    return row


def _arm_row(c, z, E, exact_i, pc):
    s = V.summarise(c, z.rmsd, E, exact_i)
    s.update(_concentration(c.all_seen(), *pc))
    return s


def _pcache(name, E, z, store):
    """(E-percentile, RMSD-percentile) for one objective, computed once and reused."""
    if "_r" not in store:
        store["_r"] = ordinal_pct(z.rmsd)
    if name not in store:
        store[name] = (ordinal_pct(E), store["_r"])
    return store[name]


def run_target(z, Es, budget=8192, seeds=(0, 1, 2), alphas=(0.25,), store=None,
               only=None):
    """Every arm, every objective, one target, at one budget.

    `only` restricts which objectives run, which is how `S14_disto_bayes` gets a full arm
    table at the middle budget without paying for it at all three.  That control matters
    more than a third point on its curve: it is the SAME distogram consumed SELECTIVELY, so
    it separates "the objective CLASS changed the answer" from "the information changed".
    """
    store = {} if store is None else store
    rows = {}
    for name, E in Es.items():
        if only is not None and name not in only:
            continue
        exact_i = int(np.argmin(E))
        pc = _pcache(name, E, z, store)
        arms = {}
        for s in seeds:
            arms.setdefault("untrained_bestofN", []).append(
                _arm_row(untrained_best_of_n(E, z.n_qubits, budget, seed=s), z, E,
                         exact_i, pc))
            arms.setdefault("random", []).append(_arm_row(
                V.search_random(E, z.n, z.k, budget, np.random.default_rng(1000 + s)),
                z, E, exact_i, pc))
            arms.setdefault("anneal", []).append(_arm_row(
                V.search_anneal(E, z.n, z.k, budget, np.random.default_rng(2000 + s)),
                z, E, exact_i, pc))
            arms.setdefault("greedy", []).append(_arm_row(
                V.search_greedy(E, z.n, z.k, budget, np.random.default_rng(3000 + s)),
                z, E, exact_i, pc))
            for a in alphas:
                arms.setdefault(f"vqe_cvar_a{a}", []).append(
                    _vqe_row(E, z, a, budget, s, pc))
        rows[name] = {"arms": {k: _agg(v) for k, v in arms.items()}}
    return rows


_AGG_FIELDS = ("rmsd_returned", "objective_gap", "structural_gap", "rmsd_best_seen",
               "rmsd_mean_seen", "evals", "conc_e_pct_first_q", "conc_e_pct_last_q",
               "conc_rmsd_pct_first_q", "conc_rmsd_pct_last_q", "conc_distinct_frac",
               "init_rmsd_best", "param_displacement", "diversity",
               "grad_norm_first", "dead_start")


def _agg(v):
    out = {"n_seeds": len(v), "per_seed_rmsd": [float(x["rmsd_returned"]) for x in v]}
    for f in _AGG_FIELDS:
        vals = [x[f] for x in v if f in x and x[f] is not None]
        if vals:
            out[f] = float(np.mean(vals))
    return out


# ==================================================================== the mechanistic pass
def skill_table(targets=ENUM19, path=None):
    """The cheap, decisive pass: ordering skill of every objective and every comparator.

    This is the readout that answers question (c) -- is the generative objective's in-tail
    ordering skill materially better than the selective energies' -- and it costs one CA
    build per target, not a single optimiser run.
    """
    path = path or os.path.join(RESULTS, "qrestraint_skill.json")
    out = {"targets": {}, "consistency": {}, "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    for pdb in targets:
        V.wait_for_memory(1.6, tag=f"qrestraint skill {pdb}")
        t0 = time.time()
        z = Enum(pdb)
        CA = all_ca(z)
        nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)
        out["consistency"][pdb] = _check_consistency(z, CA, nat)
        Es, extras = build_objectives(z, CA)
        del CA
        row = {}
        for name, E in Es.items():
            row[name] = ordering_skill(E, z.rmsd)
            row[name]["s14_tail_accuracy"] = s14_tail_accuracy(E, z.rmsd)
        row.update(comparator_skill(z, extras))
        row["_meta"] = {"n": z.n, "k": z.k, "N": z.N, "n_qubits": z.n_qubits,
                        "fold": z.fold, "seq": z.seq}
        out["targets"][pdb] = row
        del Es, z
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
        print(f"  skill {pdb} done in {time.time()-t0:.1f}s", flush=True)
    print_skill(out)
    return out


def _cls(name):
    if name == "E_ORACLE_true":
        return "ORACLE"
    return "GENERATIVE" if name.startswith("E_") else "SELECTIVE"


def print_skill(out):
    pdbs = list(out["targets"])
    names = list(OBJECTIVES) + list(COMPARATORS)
    print(f"\nTABLE 1  ORDERING SKILL over {len(pdbs)} enumerated targets -- Spearman "
          f"rho(E, RMSD).\n         The tail is where the optimiser concentrates. z is "
          f"against a matched RANDOM-TAIL null.")
    print(f"{'objective':<17}{'class':>11}{'rho glob':>10}{'rho 1%':>9}{'z':>7}"
          f"{'rho 0.1%':>10}{'z':>7}{'gap 1%':>9}{'argminRMSD':>12}{'best':>8}{'n tgt':>7}")
    for name in names:
        rows = [out["targets"][p][name] for p in pdbs if name in out["targets"][p]]
        if not rows:
            continue

        def m(k):
            return float(np.nanmean([r[k] for r in rows]))
        z1 = m("rho_tail_0.01") / max(1e-9, m("rho_tail_0.01_null_sd"))
        z2 = m("rho_tail_0.001") / max(1e-9, m("rho_tail_0.001_null_sd"))
        print(f"{name:<17}{_cls(name):>11}{m('rho_global'):>10.3f}{m('rho_tail_0.01'):>9.3f}"
              f"{z1:>7.1f}{m('rho_tail_0.001'):>10.3f}{z2:>7.1f}"
              f"{m('tail_0.01_gap'):>9.3f}{m('rmsd_of_certified_optimum'):>12.3f}"
              f"{m('rmsd_best_in_space'):>8.3f}{len(rows):>7d}")

    if not any("s14_tail_accuracy" in out["targets"][p].get(names[0], {}) for p in pdbs):
        return
    print(f"\nTABLE 2  SPRINT 14's OWN STATISTIC -- tail-restricted pairwise accuracy at "
          f"gap > {S14_GAP} A.\n         Read against the RANDOM-TAIL NULL in the last "
          f"columns, NEVER against 0.500 (Sprint 14 measured 0.524-0.527).")
    print(f"{'objective':<17}{'class':>11}" + "".join(f"{'q=' + str(q):>10}"
                                                      for q in S14_FRACS)
          + f"{'tail-bulk':>11}{'null 1%':>9}{'null 0.1%':>11}")
    for name in names:
        rows = [out["targets"][p][name]["s14_tail_accuracy"] for p in pdbs
                if name in out["targets"][p] and "s14_tail_accuracy" in out["targets"][p][name]]
        if not rows:
            continue
        vals = [float(np.nanmean([r[str(q)]["acc_gap_matched"] for r in rows]))
                for q in S14_FRACS]
        tb = float(np.nanmean([r["tail_minus_bulk"] for r in rows]))
        n1 = float(np.nanmean([r["0.01"]["null_mean"] for r in rows]))
        n2 = float(np.nanmean([r["0.001"]["null_mean"] for r in rows]))
        print(f"{name:<17}{_cls(name):>11}"
              + "".join(f"{v:>10.3f}" for v in vals)
              + f"{tb:>11.3f}{n1:>9.3f}{n2:>11.3f}")

    print(f"\nTABLE 3  SELECTION DECOMPOSITION at q=0.01 -- sel = pool + FILTERING + ORDERING."
          f"\n         Sprint 14: 'the ORDERING term's sign IS optimise-harder-get-worse'.")
    print(f"{'objective':<17}{'class':>11}{'pool mean':>11}{'tail mean':>11}"
          f"{'sel(argmin)':>13}{'FILTERING':>11}{'ORDERING':>10}")
    for name in names:
        rows = [out["targets"][p][name]["s14_tail_accuracy"]["selection_decomposition_q0.01"]
                for p in pdbs
                if name in out["targets"][p] and "s14_tail_accuracy" in out["targets"][p][name]]
        if not rows:
            continue

        def mm(k):
            return float(np.nanmean([r[k] for r in rows]))
        print(f"{name:<17}{_cls(name):>11}{mm('pool_mean'):>11.3f}{mm('tail_mean'):>11.3f}"
              f"{mm('sel_rmsd'):>13.3f}{mm('filtering'):>11.3f}{mm('ordering'):>10.3f}")


# ==================================================================== the main ladder
def run(targets=ENUM19, budgets=(2048, 8192, 32768), seeds=(0, 1, 2), alphas=(0.25,),
        path=None):
    path = path or os.path.join(RESULTS, "qrestraint.json")
    out = {"budgets": list(budgets), "seeds": list(seeds), "alphas": list(alphas),
           "targets": {}, "consistency": {}, "shots": 512, "ansatz": "mps2f",
           "baseline": "const (the CORRECTED estimator; the recorded tail-baseline defect "
                       "is NOT active in this table)",
           "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    out["cpu_at_start"] = wait_for_cpu(90.0, tag="qrestraint main")
    for pdb in targets:
        V.wait_for_memory(1.6, tag=f"qrestraint {pdb}")
        t0 = time.time()
        z = Enum(pdb)
        CA = all_ca(z)
        nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)
        out["consistency"][pdb] = _check_consistency(z, CA, nat)
        Es, extras = build_objectives(z, CA)
        del CA
        cell = {"_skill": {}, "_meta": {"n": z.n, "N": z.N, "n_qubits": z.n_qubits,
                                        "fold": z.fold}}
        for name, E in Es.items():
            cell["_skill"][name] = ordering_skill(E, z.rmsd)
            cell["_skill"][name]["s14_tail_accuracy"] = s14_tail_accuracy(E, z.rmsd)
        cell["_skill"].update(comparator_skill(z, extras))
        store = {}
        mid = sorted(budgets)[len(budgets) // 2]
        Es_all = dict(Es, **extras)
        del extras
        for b in budgets:
            cell[str(b)] = run_target(z, Es_all if b == mid else Es, budget=b, seeds=seeds,
                                      alphas=alphas, store=store)
            with open(path, "w") as fh:
                out["targets"][pdb] = cell
                json.dump(out, fh, indent=1, default=float)
            print(f"  {pdb} budget {b} done ({time.time()-t0:.0f}s cumulative)", flush=True)
        out["targets"][pdb] = cell
        del Es, z
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
    summarise_run(out, path)
    return out


def summarise_run(out, path=None):
    """The two tables that matter, plus the paired test the positive result was declared on."""
    path = path or os.path.join(RESULTS, "qrestraint.json")
    budgets = out["budgets"]
    alphas = out["alphas"]
    pdbs = [p for p in out["targets"]
            if all(str(b) in out["targets"][p] for b in budgets)]
    if not pdbs:
        print("no complete target yet")
        return out
    print_skill({"targets": {p: out["targets"][p]["_skill"] for p in pdbs}})

    arm_names = ["untrained_bestofN", "random", "anneal", "greedy"] + \
                [f"vqe_cvar_a{a}" for a in alphas]
    out["summary"] = {}
    mid = sorted(budgets)[len(budgets) // 2]
    #: the budget is only interpretable as a FRACTION OF THE SPACE, and the ladder spans two
    #: register sizes: 2048/4^9 = 0.78% but 2048/4^10 = 0.20%.  An n=10 target at 32768 sits
    #: at the same coverage as an n=9 target at 8192, which is why the two are not pooled
    #: blindly and why the per-n split is printed.
    Ns = sorted({int(out["targets"][p]["_meta"]["N"]) for p in pdbs})
    print("\nBUDGET AS A FRACTION OF THE ENUMERATED SPACE")
    print("        " + "".join(f"{'N=' + str(N):>14}" for N in Ns))
    for b in budgets:
        print(f"{b:>7} " + "".join(f"{100.0 * b / N:>13.2f}%" for N in Ns))
    print(f"\nARM TABLE -- n={len(pdbs)} targets, matched hard budget enforced by Counter."
          f"\nBOTH axes are printed: RMSD (accuracy) and objective_gap (optimisation)."
          f"\n`S14_disto_bayes` is the SELECTIVE control on the SAME distogram, run at the "
          f"middle budget ({mid}) only.")
    for name in list(OBJECTIVES) + ["S14_disto_bayes"]:
        avail = [b for b in budgets
                 if name in out["targets"][pdbs[0]].get(str(b), {})]
        if not avail:
            continue
        print(f"\n  {name}   [{_cls(name)}]")
        print(f"    {'budget':>7} {'arm':<20}{'RMSD':>8}{'objgap':>10}"
              f"{'vs untrained_bestofN (paired)':>34}{'W/L':>9}")
        for b in avail:
            for a in arm_names:
                v = np.asarray([out["targets"][p][str(b)][name]["arms"][a]["rmsd_returned"]
                                for p in pdbs], float)
                g = np.asarray([out["targets"][p][str(b)][name]["arms"][a]["objective_gap"]
                                for p in pdbs], float)
                ctrl = np.asarray([out["targets"][p][str(b)][name]["arms"]
                                   ["untrained_bestofN"]["rmsd_returned"] for p in pdbs],
                                  float)
                st = I.paired(v, ctrl, names=pdbs)
                cv = concentration_verdict(v - ctrl)
                out["summary"][f"{name}|{b}|{a}"] = {
                    "mean_rmsd": float(v.mean()), "median_rmsd": float(np.median(v)),
                    "mean_objective_gap": float(g.mean()), "vs_control": st,
                    "concentration": cv}
                tag = "" if a == "untrained_bestofN" else \
                    f"  {st['mean_diff']:+.3f} [{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]"
                wl = "" if a == "untrained_bestofN" else \
                    f"  {st['n_better']}/{len(pdbs)}"
                print(f"    {b:>7} {a:<20}{v.mean():>8.3f}{g.mean():>10.3f}{tag:>34}{wl:>9}")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    I.write("s15_qrestraint", dict(out, per_target=out["targets"]), n_expected=len(pdbs))
    print("\n(a positive result requires beating untrained_bestofN with a CI excluding zero;"
          "\n beating `random` proves nothing -- a constant alpha-helix does that.)")
    return out


# ==================================================================== the CVaR dimension
ALPHAS = (0.02, 0.05, 0.1, 0.25, 0.5, 1.0)


def alpha_sweep(targets=ENUM19, objs=("E_ls_pred", "E_combined", "E_ORACLE_true"),
                budget=8192, seeds=(0, 1, 2), alphas=ALPHAS, path=None):
    """The CVaR dimension, with the recorded defects addressed rather than asserted away.

    Sprint 14 recorded three CVaR defects.  Their status in THIS experiment:

      1. `baseline="tail"` gradient bias.  NOT ACTIVE in the main table -- `R.run` defaults
         to `baseline="const"`, the corrected control variate.  Rather than merely say so,
         the defective estimator is run here as its own arm at matched everything, so the
         size of the defect is measured on this objective class instead of quoted from a
         different one.
      2. Sampled-CVaR upward bias at non-integer `alpha*N`.  ACTIVE wherever
         `alpha * shots` is not an integer.  At shots=512 the swept alphas give
         alpha*N = 10.24, 25.6, 51.2, 128, 256, 512 -- so it is active at 0.02, 0.05 and 0.1
         and inactive at 0.25, 0.5, 1.0.  Flagged per row.
      3. `dCVaR/dp` identically zero iff `p(argmin E) >= alpha`.  Measured per run rather
         than assumed: `dead_start` counts runs whose first-iteration gradient norm is 0.
    """
    path = path or os.path.join(RESULTS, "qrestraint_alpha.json")
    out = {"budget": budget, "seeds": list(seeds), "alphas": list(alphas),
           "objectives": list(objs), "shots": 512, "targets": {},
           "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    out["cpu_at_start"] = wait_for_cpu(90.0, tag="qrestraint alpha")
    for pdb in targets:
        V.wait_for_memory(1.6, tag=f"qrestraint alpha {pdb}")
        t0 = time.time()
        z = Enum(pdb)
        CA = all_ca(z)
        Es, _extras = build_objectives(z, CA)
        del CA, _extras
        cell = {}
        store = {}
        for name in objs:
            E = Es[name]
            exact_i = int(np.argmin(E))
            pc = _pcache(name, E, z, store)
            arms = {}
            for s in seeds:
                arms.setdefault("untrained_bestofN", []).append(
                    _arm_row(untrained_best_of_n(E, z.n_qubits, budget, seed=s),
                             z, E, exact_i, pc))
                for a in alphas:
                    arms.setdefault(f"vqe_const_a{a}", []).append(
                        _vqe_row(E, z, a, budget, s, pc, baseline="const"))
                for a in (0.1, 0.25):
                    arms.setdefault(f"vqe_TAILDEFECT_a{a}", []).append(
                        _vqe_row(E, z, a, budget, s, pc, baseline="tail"))
            cell[name] = {"arms": {k: _agg(v) for k, v in arms.items()}}
        out["targets"][pdb] = cell
        del Es, z
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
        print(f"  alpha {pdb} done ({time.time()-t0:.0f}s)", flush=True)
    summarise_alpha(out, path)
    return out


def summarise_alpha(out, path=None):
    path = path or os.path.join(RESULTS, "qrestraint_alpha.json")
    pdbs = list(out["targets"])
    shots = out["shots"]
    arm_names = ["untrained_bestofN"] + [f"vqe_const_a{a}" for a in out["alphas"]] + \
                [f"vqe_TAILDEFECT_a{a}" for a in (0.1, 0.25)]
    out["summary"] = {}
    for name in out["objectives"]:
        print(f"\n  CVaR alpha sweep -- {name}  (n={len(pdbs)}, budget {out['budget']})")
        print(f"    {'arm':<24}{'a*shots':>9}{'RMSD':>8}{'objgap':>10}"
              f"{'vs untrained (paired)':>30}{'W/L':>8}")
        ctrl = np.asarray([out["targets"][p][name]["arms"]["untrained_bestofN"]
                           ["rmsd_returned"] for p in pdbs], float)
        for a in arm_names:
            rows = [out["targets"][p][name]["arms"].get(a) for p in pdbs]
            if any(r is None for r in rows):
                continue
            v = np.asarray([r["rmsd_returned"] for r in rows], float)
            g = np.asarray([r["objective_gap"] for r in rows], float)
            st = I.paired(v, ctrl, names=pdbs)
            out["summary"][f"{name}|{a}"] = {"mean_rmsd": float(v.mean()),
                                             "mean_objective_gap": float(g.mean()),
                                             "vs_control": st}
            al = a.split("_a")[-1]
            an_ = f"{float(al)*shots:.2f}" if a != "untrained_bestofN" else "-"
            tag = "" if a == "untrained_bestofN" else \
                f"  {st['mean_diff']:+.3f} [{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]"
            wl = "" if a == "untrained_bestofN" else f"  {st['n_better']}/{len(pdbs)}"
            print(f"    {a:<24}{an_:>9}{v.mean():>8.3f}{g.mean():>10.3f}{tag:>30}{wl:>8}")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    return out


# ============================================ the null-calibrated concentration check
def concentration_verdict(d, k=None, n_sim=4000, seed=0):
    """ONE PASS/FAIL verdict on whether a paired effect is carried by a few targets.

    The brief is explicit that a raw drop-top threshold is NOT a valid test: when mean/sd is
    small, discarding the k best removes much of the total EVEN IF every target carries an
    identical effect.  So the observed drop-top statistic is compared to a **simulated
    uniform-effect null** -- every target carrying the SAME effect `mean(d)` plus i.i.d.
    noise at the observed `sd(d)` -- and the verdict is its percentile in that null.  It is
    returned as one field, not as separate numbers a reader can select from, because the
    brief records a claim published on the reassuring half of exactly this dict.

    `mean_over_sd` is printed with it: when it is small the test HAS NO POWER and the
    verdict says so rather than reading PASS.
    """
    d = np.asarray(d, float)
    n = d.size
    if n < 6:
        return {"verdict": "NO POWER (n < 6)", "n": int(n)}
    k = k or max(1, n // 4)
    rng = np.random.default_rng(seed)
    mu, sd = float(d.mean()), float(d.std(ddof=1))

    def stat(x):
        o = np.argsort(x)                       # most negative (largest gains) first
        return float(x[o[k:]].mean())
    obs = stat(d)
    sim = np.asarray([stat(mu + sd * rng.standard_normal(n)) for _ in range(n_sim)])
    pct = float((sim < obs).mean())
    mos = mu / sd if sd > 0 else float("inf")
    if abs(mos) < 0.3:
        v = "NO POWER (|mean/sd| < 0.3)"
    elif 0.025 <= pct <= 0.975:
        v = "PASS (effect is uniform; drop-top is what a uniform effect would give)"
    else:
        v = "FAIL (CONCENTRATED; the effect is carried by a few targets)"
    return {"verdict": v, "n": int(n), "k_dropped": int(k), "mean": mu, "sd": sd,
            "mean_over_sd": float(mos), "median": float(np.median(d)),
            "drop_top_k_mean": obs, "null_pctile_of_drop_top_k": pct,
            "null_mean": float(sim.mean()), "null_p2.5": float(np.percentile(sim, 2.5)),
            "null_p97.5": float(np.percentile(sim, 97.5))}


# ================================= did the VARIATIONAL STATE learn anything? (n=9 only)
def mode_probe(targets=ENUM9, objs=OBJECTIVES, budget=8192, seeds=(0, 1, 2), alpha=0.25,
               path=None):
    """The `mode` readout, which `exact_dist=False` cannot give.

    `s14/vqe_run.py` states the distinction and it is the sharpest one available: `best_seen`
    is a best-of-N SELECTOR, so it makes a VQE "a (bad) random search dressed in a circuit
    unless the distribution actually concentrates", while `mode` -- the argmax of the FINAL
    distribution -- "is the only readout for which 'the VQE did something' is a meaningful
    claim".  The main ladder runs with `exact_dist=False` because an n=10 register would need
    2**20 amplitudes per readout; this probe re-runs the identical configuration on the nine
    n=9 targets with `exact_dist=True` and reports mode, entropy, effective support and the
    probability mass the trained state puts below 2.0 / 2.5 / 3.0 A.

    It is a DIAGNOSTIC, not an arm: it changes no budget and enters no headline comparison.
    """
    path = path or os.path.join(RESULTS, "qrestraint_mode.json")
    out = {"budget": budget, "alpha": alpha, "seeds": list(seeds), "targets": {},
           "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    out["cpu_at_start"] = wait_for_cpu(90.0, tag="qrestraint mode")
    keep = ("mode_i", "mode_rmsd", "entropy_bits", "eff_support", "max_prob",
            "mean_rmsd_under_p", "pmass_below_2.0", "pmass_below_2.5", "pmass_below_3.0",
            "rmsd_returned", "rmsd_best_seen", "init_rmsd_best", "init_rmsd_mean",
            "param_displacement", "diversity")
    for pdb in targets:
        V.wait_for_memory(1.6, tag=f"qrestraint mode {pdb}")
        z = Enum(pdb)
        CA = all_ca(z)
        Es, _x = build_objectives(z, CA)
        del CA, _x
        cell = {}
        for name in objs:
            rows = []
            for s in seeds:
                r = R.run(Es[name], z.n_qubits, alpha, budget, shots=512, ansatz="mps2f",
                          seed=s, rmsd=z.rmsd, bits_per_res=z.bits_per_res,
                          exact_dist=True, baseline="const")
                rows.append({k: float(r[k]) for k in keep if k in r})
            cell[name] = {k: float(np.mean([x[k] for x in rows if k in x]))
                          for k in keep if any(k in x for x in rows)}
            #: the uniform-distribution references the entropy must be read against
            cell[name]["entropy_bits_uniform"] = float(z.n_qubits)
            cell[name]["rmsd_best_in_space"] = float(z.rmsd.min())
            cell[name]["rmsd_mean_in_space"] = float(z.rmsd.mean())
        out["targets"][pdb] = cell
        del Es, z
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
        print(f"  mode {pdb} done", flush=True)
    print_mode(out)
    return out


def untrained_enrichment(targets=ENUM9, seeds=(0, 1, 2), thresholds=(2.0, 2.5, 3.0),
                         path=None):
    """THE CONTROL `mode_probe` NEEDS: the UNTRAINED circuit's own near-native mass.

    `mode_probe` reports P(RMSD < t) under the TRAINED distribution.  Quoting that against a
    UNIFORM base rate would credit the optimiser with whatever bias the ansatz already carries
    at `theta0` -- and `theta0` is precisely the control this whole workstream rests on.  So
    the identical quantity is computed at `theta0`, on the same targets and the same seeds,
    by exact enumeration of the register (2**18 amplitudes, no sampling).

    Measured: the untrained circuit is slightly WORSE than uniform (0.85x at 2.0 A), so the
    enrichment measured against it is larger than against uniform, not smaller.
    """
    path = path or os.path.join(RESULTS, "qrestraint_untrained_enrich.json")
    from core import quantum as Qq
    out = {}
    for pdb in targets:
        V.wait_for_memory(1.6, tag=f"qrestraint enrich {pdb}")
        z = Enum(pdb)
        bits = Qq.all_bitstrings(z.n_qubits)
        an = R.make_ansatz("mps2f", z.n_qubits)
        rows = []
        for s in seeds:
            th0 = R.init_theta(an, np.random.default_rng(s), 0.8, "random")
            pr = np.exp(np.asarray(an.logp(th0, bits), float))
            pr = np.maximum(pr, 0.0)
            pr = pr / pr.sum()
            ent = float(-(pr[pr > 0] * np.log2(pr[pr > 0])).sum())
            rows.append({"entropy_bits": ent, "eff_support": float(2.0 ** ent),
                         "mode_rmsd": float(z.rmsd[int(np.argmax(pr))]),
                         "mean_rmsd_under_p": float(pr @ z.rmsd),
                         **{f"pmass_below_{t}": float(pr[z.rmsd < t].sum())
                            for t in thresholds}})
        out[pdb] = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
        #: the uniform base rate of the same space, so a reader can see BOTH references
        out[pdb]["uniform"] = {f"pmass_below_{t}": float((z.rmsd < t).mean())
                               for t in thresholds}
        print(f"  enrich {pdb} done", flush=True)
        del z, bits, an
        with open(path, "w") as fh:
            json.dump(out, fh, indent=1, default=float)
    return out


def print_mode(out):
    pdbs = list(out["targets"])
    print(f"\nTABLE 6  THE `mode` READOUT, n={len(pdbs)} n=9 targets, budget {out['budget']},"
          f" alpha {out['alpha']}.\n         `mode` is the argmax of the FINAL distribution "
          f"-- the only readout for which 'the VQE did\n         something' is meaningful. "
          f"Entropy is in bits against a uniform 18.")
    print(f"{'objective':<17}{'mode RMSD':>11}{'best_seen':>11}{'init best':>11}"
          f"{'space best':>12}{'space mean':>12}{'entropy':>9}{'eff supp':>10}"
          f"{'P(<2.5A)':>10}")
    for name in OBJECTIVES:
        rows = [out["targets"][p][name] for p in pdbs if name in out["targets"][p]]
        if not rows:
            continue

        def m(k):
            return float(np.nanmean([r[k] for r in rows if k in r]))
        print(f"{name:<17}{m('mode_rmsd'):>11.3f}{m('rmsd_best_seen'):>11.3f}"
              f"{m('init_rmsd_best'):>11.3f}{m('rmsd_best_in_space'):>12.3f}"
              f"{m('rmsd_mean_in_space'):>12.3f}{m('entropy_bits'):>9.2f}"
              f"{m('eff_support'):>10.0f}{m('pmass_below_2.5'):>10.4f}")


# ======================================== the PRE-DECLARED primary comparison, question (a)
def primary_table(out, alpha=0.25):
    """The one comparison the positive result was declared on, with everything a reader needs.

    Declared in the module docstring BEFORE any run: VQE must beat `untrained_bestofN` on
    returned RMSD, paired across targets, with a CI excluding zero, at matched budget.  This
    prints exactly that, plus the median (the brief's free early warning when it disagrees
    with the mean), the null-calibrated concentration verdict, and an n=9 / n=10 split --
    because the two register sizes are different instruments and pooling them without showing
    the split would hide a disagreement if one existed.
    """
    arm = f"vqe_cvar_a{alpha}"
    print(f"\nTABLE 0  THE PRE-DECLARED PRIMARY COMPARISON: {arm} vs untrained_bestofN."
          f"\n         Negative = VQE better. A positive result needs a CI excluding zero.")
    print(f"{'objective':<17}{'budget':>7}{'n':>4}{'mean diff':>11}{'95% CI':>20}"
          f"{'median':>9}{'W/L':>8}{'n=9':>8}{'n=10':>8}  concentration")
    res = {}
    for name in list(OBJECTIVES) + ["S14_disto_bayes"]:
        for b in out["budgets"]:
            pdbs = [p for p in out["targets"]
                    if str(b) in out["targets"][p] and name in out["targets"][p][str(b)]]
            if len(pdbs) < 3:
                continue

            def col(a_):
                return np.asarray([out["targets"][p][str(b)][name]["arms"][a_]
                                   ["rmsd_returned"] for p in pdbs], float)
            d = col(arm) - col("untrained_bestofN")
            nn = np.asarray([int(out["targets"][p]["_meta"]["n"]) for p in pdbs])
            st = I.paired(col(arm), col("untrained_bestofN"), names=pdbs)
            cv = concentration_verdict(d)
            s9 = float(d[nn == 9].mean()) if (nn == 9).any() else float("nan")
            s10 = float(d[nn == 10].mean()) if (nn == 10).any() else float("nan")
            res[f"{name}|{b}"] = {"paired": st, "concentration": cv,
                                  "mean_diff_n9": s9, "mean_diff_n10": s10}
            ci = f"[{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]"
            excl = "*" if (st["ci95"][0] > 0) or (st["ci95"][1] < 0) else " "
            print(f"{name:<17}{b:>7}{st['n']:>4}{st['mean_diff']:>+11.3f}{ci:>20}"
                  f"{st['median_diff']:>+9.3f}{st['n_better']:>4}/{st['n']:<3}"
                  f"{s9:>+8.3f}{s10:>+8.3f}  {excl}{cv['verdict'].split('(')[0].strip()}")
    print("  * = CI excludes zero")
    out["primary"] = res
    return res


# =============================================== the control that a quantum claim needs
def classical_table(out, alpha=0.25):
    """VQE against CLASSICAL search at the same matched budget.

    `untrained_bestofN` is the control that decides whether OPTIMISING beat NOT optimising --
    the Sprint 14 question.  It is not the control that decides whether anything QUANTUM
    happened.  For that, the arm has to beat 1-opt greedy and simulated annealing at the same
    hard budget, and the brief's rule 'a classical search is not VQE' cuts both ways: if
    greedy matches or beats the VQE, the correct claim is that SEARCH works on this
    objective, not that the circuit does.  This table is printed with equal prominence for
    exactly that reason.
    """
    arm = f"vqe_cvar_a{alpha}"
    print(f"\nTABLE 5  VQE vs CLASSICAL SEARCH at matched budget. Negative = VQE better."
          f"\n         A quantum claim needs this column, not only the untrained control.")
    print(f"{'objective':<17}{'budget':>8}{'vs greedy':>28}{'W/L':>7}"
          f"{'vs anneal':>28}{'W/L':>7}")
    res = {}
    for name in list(OBJECTIVES) + ["S14_disto_bayes"]:
        for b in out["budgets"]:
            pdbs = [p for p in out["targets"]
                    if str(b) in out["targets"][p] and name in out["targets"][p][str(b)]]
            if len(pdbs) < 2:
                continue

            def col(a_):
                return np.asarray([out["targets"][p][str(b)][name]["arms"][a_]
                                   ["rmsd_returned"] for p in pdbs], float)
            v = col(arm)
            line = f"{name:<17}{b:>8}"
            for c in ("greedy", "anneal"):
                st = I.paired(v, col(c), names=pdbs)
                res[f"{name}|{b}|vs_{c}"] = st
                line += (f"  {st['mean_diff']:+.3f} [{st['ci95'][0]:+.3f},"
                         f"{st['ci95'][1]:+.3f}]{st['n_better']:>4}/{len(pdbs)}")
            print(line)
    out["vs_classical"] = res
    return res


# ============================ is the budget trap a property of the OBJECTIVE? (the law)
def budget_trap_table(out, arm="greedy"):
    """Test the brief's standing law on a new objective class, per (objective, target) cell.

    The brief states: *"the budget trap is a property of BAD OBJECTIVES. Searching harder
    makes structures worse on a bad objective, is neutral on a mediocre one, and HELPS
    monotonically on a good one"*, and requires the objective-quality condition to be stated
    with any budget claim.  That law was measured on Sprint 14's signal-tunable ORACLE family.
    Here it is tested on real, native-free objectives.

    For every (objective, target) cell:
        x = in-tail ordering skill  = tail-restricted pairwise accuracy at q=0.01
                                      MINUS that objective's own matched random-tail null
        y = the value of searching harder = RMSD(arm @ max budget) - RMSD(arm @ min budget)

    `greedy` is the default arm because it is the one that actually reaches the certified
    optimum (objective gap 0.000 at the top budget on most cells), so `y` is as close as this
    instrument gets to "the cost of moving from a broad sample to the argmin".
    If the law holds, rho(x, y) is NEGATIVE: more in-tail skill, more benefit from searching.
    """
    lo, hi = min(out["budgets"]), max(out["budgets"])
    xs, ys, tags = [], [], []
    for name in list(OBJECTIVES) + ["S14_disto_bayes"]:
        for p in out["targets"]:
            c = out["targets"][p]
            if str(lo) not in c or str(hi) not in c:
                continue
            if name not in c[str(lo)] or name not in c[str(hi)]:
                continue
            sk = c["_skill"].get(name, {}).get("s14_tail_accuracy")
            if not sk:
                continue
            a, nul = sk["0.01"]["acc_gap_matched"], sk["0.01"]["null_mean"]
            if not (np.isfinite(a) and np.isfinite(nul)):
                continue
            xs.append(a - nul)
            ys.append(c[str(hi)][name]["arms"][arm]["rmsd_returned"]
                      - c[str(lo)][name]["arms"][arm]["rmsd_returned"])
            tags.append((name, p))
    xs, ys = np.asarray(xs), np.asarray(ys)
    if xs.size < 10:
        return {}
    r = float(V.spearman(xs, ys))
    #: THE CONFOUND, reported rather than hidden.  `x` is accuracy-minus-null, and the brief
    #: records that an in-tail rank statistic is NOT monotone in objective quality because a
    #: better objective has a NARROWER tail.  `E_ORACLE_true` is the extreme case: its bulk
    #: accuracy is ~0.87, so its tail sits far below its own null and its `x` is the most
    #: negative in the table while being by far the best objective.  The correlation is
    #: therefore also reported with the ORACLE cells removed, and both are quoted together.
    keep = np.asarray([t[0] != "E_ORACLE_true" for t in tags])
    r_no = float(V.spearman(xs[keep], ys[keep])) if keep.sum() >= 10 else float("nan")
    print(f"\nTABLE 7  IS THE BUDGET TRAP A PROPERTY OF THE OBJECTIVE?  arm={arm}, "
          f"budget {lo} -> {hi}, n={xs.size} (objective, target) cells."
          f"\n         x = in-tail accuracy at q=0.01 minus its own matched null."
          f"\n         y = RMSD({hi}) - RMSD({lo});  NEGATIVE y = searching harder HELPED."
          f"\n         The law predicts rho(x, y) < 0.")
    print(f"  rho(in-tail skill, value of searching harder) = {r:+.3f}   "
          f"(ORACLE cells removed: {r_no:+.3f}, n={int(keep.sum())})   "
          f"mean y = {ys.mean():+.3f} A   frac(y<0) = {(ys < 0).mean():.2f}")
    per = {}
    print(f"\n  {'objective':<17}{'class':>11}{'mean x':>9}{'mean y':>9}{'frac y<0':>10}"
          f"{'n':>5}")
    for name in list(OBJECTIVES) + ["S14_disto_bayes"]:
        m = np.asarray([t[0] == name for t in tags])
        if not m.any():
            continue
        per[name] = {"mean_x": float(xs[m].mean()), "mean_y": float(ys[m].mean()),
                     "frac_helped": float((ys[m] < 0).mean()), "n": int(m.sum())}
        print(f"  {name:<17}{_cls(name):>11}{xs[m].mean():>+9.3f}{ys[m].mean():>+9.3f}"
              f"{(ys[m] < 0).mean():>10.2f}{int(m.sum()):>5}")
    out.setdefault("budget_trap", {})[arm] = {
        "arm": arm, "rho": r, "rho_without_oracle": r_no, "n_cells": int(xs.size),
        "mean_y": float(ys.mean()), "frac_helped": float((ys < 0).mean()),
        "per_objective": per}
    out["budget_trap"] = out["budget_trap"]
    return out["budget_trap"][arm]


# =========================================================== the mechanism, question (d)
def mechanism_table(out, budget=None, alpha=0.25):
    """If VQE loses DESPITE better in-tail ordering, what is carrying the loss?

    Sprint 14's mechanism was tail-blindness: concentrate onto a tail, and if the objective
    cannot order inside it the argmin is a random draw.  That explanation is only available
    where in-tail ordering is at chance.  Where it is NOT, the loss has to come from
    somewhere else, and there are only a few places it can come from.  Each is a column:

      `conc`        did the optimiser concentrate at all?  E-percentile of the last quarter
                    of samples minus the first quarter.  Negative = concentrated onto low E.
                    If this is ~0 the Sprint 14 story does not even apply: nothing happened.
      `dist_frac`   fraction of the budget spent on DISTINCT configurations.  The readout of
                    both arms is best-of-N, and best-of-N is a MAXIMUM ORDER STATISTIC: it
                    improves with the number of independent draws.  Concentration reduces
                    that number.  This is a loss channel that is completely independent of
                    ordering skill, and it is the one that survives when tail-blindness does
                    not.
      `d_objgap`    VQE's objective gap minus the control's.  NEGATIVE means the VQE really
                    is the better OPTIMISER.  Sprint 14's whole point is that this and the
                    RMSD column come apart, so they are never substituted for one another.
      `d_rmsd`      VQE's returned RMSD minus the control's.  POSITIVE = VQE is worse.
      `selgap`      returned RMSD minus the BEST RMSD the arm actually sampled -- Sprint 14's
                    selection gap, the cost of the objective picking the wrong member of what
                    it saw.  It is printed for BOTH arms because it is a property of the
                    OBJECTIVE, not of the optimiser, and reading the VQE's alone would invite
                    exactly the misattribution this table exists to prevent.
      `erase`       returned RMSD minus the initialisation's own best.  Sprint 14 measured
                    that VQE "erases its own initialisation".  Read it against `selgap_c`:
                    the control has the same arithmetic working on it.
    """
    budgets = out["budgets"]
    b = str(budget or sorted(budgets)[len(budgets) // 2])
    pdbs = [p for p in out["targets"] if b in out["targets"][p]]
    arm = f"vqe_cvar_a{alpha}"
    print(f"\nTABLE 4  MECHANISM at budget {b}, n={len(pdbs)} targets, {arm} vs "
          f"untrained_bestofN.\n         `conc` < 0 means the optimiser concentrated; "
          f"`dist_frac` is the diversity a best-of-N readout consumes.")
    print(f"{'objective':<17}{'class':>11}{'conc':>8}{'dist_frac':>11}{'ctrl_dist':>11}"
          f"{'d_objgap':>10}{'d_rmsd':>9}{'selgap_v':>10}{'selgap_c':>10}{'erase':>8}")
    rows = {}
    for name in list(OBJECTIVES) + ["S14_disto_bayes"]:
        cells = [out["targets"][p][b][name]["arms"] for p in pdbs
                 if name in out["targets"][p][b]]
        if not cells:
            continue

        def m(a_, k):
            return float(np.nanmean([c[a_][k] for c in cells if k in c[a_]]))
        conc = m(arm, "conc_e_pct_last_q") - m(arm, "conc_e_pct_first_q")
        r = {"conc_shift": conc,
             "vqe_distinct_frac": m(arm, "conc_distinct_frac"),
             "control_distinct_frac": m("untrained_bestofN", "conc_distinct_frac"),
             "d_objective_gap": m(arm, "objective_gap") - m("untrained_bestofN",
                                                            "objective_gap"),
             "d_rmsd": m(arm, "rmsd_returned") - m("untrained_bestofN", "rmsd_returned"),
             "erase": m(arm, "rmsd_returned") - m(arm, "init_rmsd_best"),
             "selgap_vqe": m(arm, "rmsd_returned") - m(arm, "rmsd_best_seen"),
             "selgap_control": (m("untrained_bestofN", "rmsd_returned")
                                - m("untrained_bestofN", "rmsd_best_seen")),
             "n": len(cells)}
        rows[name] = r
        print(f"{name:<17}{_cls(name):>11}{r['conc_shift']:>8.3f}"
              f"{r['vqe_distinct_frac']:>11.3f}{r['control_distinct_frac']:>11.3f}"
              f"{r['d_objective_gap']:>10.3f}{r['d_rmsd']:>9.3f}"
              f"{r['selgap_vqe']:>10.3f}{r['selgap_control']:>10.3f}{r['erase']:>8.3f}")
    out.setdefault("mechanism", {})[b] = rows
    return rows


# ==================================================================== reporting on demand
def report(path=None):
    """Summarise whatever the checkpoint holds, complete or not.

    The brief records that Sprint 14 lost an ansatz study by writing its JSON only on
    completion.  `run` checkpoints after every (target, budget); this reads that file and
    prints the tables for the targets that are complete, so a long run is readable at any
    point and an interrupted one is not lost.  The n it prints is the n it actually used.
    """
    path = path or os.path.join(RESULTS, "qrestraint.json")
    with open(path) as fh:
        out = json.load(fh)
    done = [p for p in out["targets"]
            if all(str(b) in out["targets"][p] for b in out["budgets"])]
    print(f"checkpoint holds {len(out['targets'])} targets, {len(done)} complete "
          f"across all budgets {out['budgets']}")
    summarise_run(out, path)
    primary_table(out)
    classical_table(out)
    for b in out["budgets"]:
        mechanism_table(out, budget=b)
    budget_trap_table(out)
    budget_trap_table(out, arm="vqe_cvar_a0.25")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    return out


def report_alpha(path=None):
    path = path or os.path.join(RESULTS, "qrestraint_alpha.json")
    with open(path) as fh:
        out = json.load(fh)
    return summarise_alpha(out, path)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "main"
    if what == "skill":
        skill_table()
    elif what == "alpha":
        alpha_sweep()
    elif what == "mode":
        mode_probe()
    elif what == "enrich":
        untrained_enrichment()
    elif what == "main":
        run()
    elif what == "report":
        report()
    elif what == "report_alpha":
        report_alpha()
    else:
        raise SystemExit(f"unknown mode {what!r}")
