"""SPRINT 17 / QUANTUM -- machinery for the (member error, diversity) PARETO test.

THE ONLY LIVE QUESTION.  Sprint 16 closed the search framing (CLAIMS I1-I12): CVaR-VQE is
significantly worse than greedy 1-opt at 10/10 objective-quality rungs, indistinguishable
from annealing at a QUARTER of its budget at 9/10, and reaches the certified global optimum
in 0-32% of cells against greedy's 68-100%.  What survives is that the VQE *samples* and
*represents* -- so the only framing in which the pillar is still live is:

    can a CVaR-VQE sampler occupy a point on the (member error, diversity) Pareto frontier
    that is inaccessible to classical samplers at MATCHED computational budget?

THE PLANE IS FORCED BY AN EXACT IDENTITY, not chosen.  For a candidate set W = {X_1..X_m}
placed in a COMMON FRAME with mean A, and any rigid map taking A onto the native T,

    RMSD(A, T)^2  =  (1/m) sum_i RMSD_fixedframe(X_i, T)^2  -  (1/m) sum_i RMSD(X_i, A)^2
    readout^2     =  M^2                                    -  D^2

so a candidate-set readout is a function of EXACTLY those two coordinates.  Reporting only
the readout is what hid this mechanism for two sprints.

THE TWO TRAPS THE BRIEF NAMES, and how this module avoids them:
  * `M` is the COMMON-FRAME member error, NOT the free-superposition RMSD.  Sprint 16 read
    one as the other and its attribution was wrong by ~2x.  `md_plane` never re-superposes
    a member.
  * the mean is QUADRATIC.  Sprint 16 took an arithmetic mean of two RMSDs where the
    identity requires the quadratic mean and lost 54% of a "prediction error".  Every `M`
    here is `sqrt(mean(err^2))`.
`verify_identity` asserts the residual is < 1e-9 A^2 before any claim is made.

THE CLASSICAL CONTROL SPRINT 16 NEVER RAN.  Its thermostats were `tilt_samples` (a reweight
of a drawn multiset -- cannot leave its support) and `tilt_exact` (reads the whole register,
32-128x the budget).  Neither is a matched-budget thermostat that MOVES support.  This module
adds `search_metropolis`: a fixed-temperature single-residue Metropolis chain with restarts,
budget-capped by the same `Counter`, swept over a temperature ladder so the classical arm
traces a CURVE on the (M, D) plane rather than a point.  That curve is the frontier the
quantum arm has to beat.

ORACLE LABELLING.  `M`, `readout`, `set_best` and every RMSD are ORACLE quantities: they read
the native and are used for post-hoc scoring only.  `D`, `n_distinct` and the draw entropy
are NATIVE-FREE.  No arm's parameters, stopping rule or temperature is chosen with a native
quantity.
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

from s12 import instrument as I                      # noqa: E402
from s14 import vqe_lib as V                         # noqa: E402
from s14 import vqe_run as R                         # noqa: E402
from s15 import seed as SD                           # noqa: E402
from s16 import qphase_lib as QP                     # noqa: E402

RESULTS = os.path.join(ROOT, "s17", "results")
CACHE = os.path.join(ROOT, "s17", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

SALT = "s17quantum"
BUDGET = QP.BUDGET          # 8192 objective evaluations, the sprint-16 convention
SHOTS = QP.SHOTS            # 512
ANSATZ = QP.ANSATZ          # mps2f
TARGETS9 = QP.TARGETS9
TARGETS10 = QP.TARGETS10
TARGETS19 = QP.TARGETS19

inst = QP.inst
Inst = QP.Inst
paired = QP.paired
cpu_pct = QP.cpu_pct
wait_for_cpu = QP.wait_for_cpu

# ======================================================================= checkpoint
_CK = {}


def ck(tag, key, value):
    """Merge-on-write incremental checkpoint into `s17/results/quantum_<tag>.json`.

    Same contract as `s16.qphase_lib.ck` but with this sprint's own path, so an interrupted
    run is readable and a rerun skips completed cells.
    """
    path = os.path.join(RESULTS, f"quantum_{tag}.json")
    if os.path.exists(path):
        try:
            with open(path) as fh:
                disk = json.load(fh)
            mem = _CK.setdefault(tag, {})
            for k, v in disk.items():
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
    path = os.path.join(RESULTS, f"quantum_{tag}.json")
    if os.path.exists(path):
        with open(path) as fh:
            d = json.load(fh)
        _CK[tag] = d
        return d
    return _CK.setdefault(tag, {})

_GATE = {"n": 0, "samples": []}


def gate(tag, every=8):
    """A CHEAP contention record in place of a blocking CPU gate, and why.

    `s16.qphase_lib.wait_for_cpu` spawns a PowerShell `Get-CimInstance` per poll.  On this
    box that call itself takes tens of seconds under load, and with three sibling
    workstreams holding the machine at a steady 100% the gate NEVER OPENS -- so blocking on
    it costs more wall clock than the work and buys nothing, because the load it is
    protecting against is not mine to release.  Sprint 16 hit the same wall and recorded
    "PROCEEDING ON A CONTENDED BOX and recording it".

    This samples the CPU every `every` cells, records it, and never blocks.  Every wall-clock
    number in this workstream must therefore be read as CONTENDED, and the samples are
    written into the artefact so a reader can see how contended.
    """
    _GATE["n"] += 1
    if (_GATE["n"] - 1) % every == 0:
        c = cpu_pct()
        _GATE["samples"].append({"tag": tag, "cpu_pct": c, "free_gb": I.free_gb()})
        print(f"  [{tag}] contention sample: CPU {c:.0f}%, free {I.free_gb():.1f} GB",
              flush=True)
    return _GATE["samples"]


#: the programme's empirical false-positive floor, printed on every row
FLOOR = 0.08
#: the Pareto dominance tolerance, pre-registered at a quarter of the floor
TAU = 0.02


# ============================================================ the (M, D) plane
def md_plane(W, nat):
    """THE EXACT DECOMPOSITION.  `W` (m, n, 3) candidate traces, `nat` (n, 3) ORACLE native.

    Returns member error `M` (quadratic mean, COMMON FRAME), diversity `D`, the coordinate
    -average readout, and the residual of the identity.  Also the free-superposition member
    error `M_free` -- reported beside `M` precisely because reading one as the other is the
    recorded Sprint 16 error, and the gap between them is printed rather than assumed small.
    """
    W = np.asarray(W, float)
    nat = np.asarray(nat, float)
    m, n, _ = W.shape
    if m == 1:
        r = float(I.kabsch_rmsd_batch(W, nat)[0])
        return {"M": r, "D": 0.0, "readout": r, "identity_resid": 0.0,
                "M_free": r, "set_best": r, "set_worst": r, "m": 1}
    # 1. common frame: the shipped operator's frame -- superpose on the medoid
    P = I.pairwise_rmsd(W)
    b = I.medoid(P)
    Wc = I.superpose_batch(W, W[b])
    A = Wc.mean(0)
    # 2. the rigid map taking A onto the native, applied to EVERY member unchanged
    Ac = A - A.mean(0)
    Tc = nat - nat.mean(0)
    H = Ac.T @ Tc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Dm = np.diag([1.0, 1.0, d])
    Rot = Vt.T @ Dm @ U.T
    Y = np.einsum("ij,mnj->mni", Rot, Wc - A.mean(0)[None, None, :]) + nat.mean(0)
    Ay = Y.mean(0)
    # 3. the three terms, all quadratic means, no member re-superposed
    err2 = ((Y - nat[None]) ** 2).sum(-1).mean(-1)        # per-member MSD, fixed frame
    dev2 = ((Y - Ay[None]) ** 2).sum(-1).mean(-1)         # per-member MSD about the mean
    M2 = float(err2.mean())
    D2 = float(dev2.mean())
    rd2 = float(((Ay - nat) ** 2).sum(-1).mean())
    return {"M": float(np.sqrt(M2)), "D": float(np.sqrt(D2)),
            "readout": float(np.sqrt(max(rd2, 0.0))),
            "identity_resid": float(rd2 - (M2 - D2)),
            "M_free": float(np.sqrt((I.kabsch_rmsd_batch(W, nat) ** 2).mean())),
            "set_best": float(I.kabsch_rmsd_batch(W, nat).min()),
            "set_worst": float(I.kabsch_rmsd_batch(W, nat).max()),
            "m": int(m)}


def verify_identity(seed=0, trials=6):
    """P0. The identity must hold to machine precision before any claim is made."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    rows = []
    for pdb in ("1CS9", "1N9U"):
        ins = inst(pdb)
        for _ in range(trials):
            m = int(rng.integers(3, 80))
            idx = rng.integers(0, ins.N, m)
            r = md_plane(ins.ca(idx), ins.nat)
            worst = max(worst, abs(r["identity_resid"]))
            rows.append({"pdb": pdb, "m": m, **{k: r[k] for k in
                                                ("M", "D", "readout", "M_free",
                                                 "identity_resid")}})
    return {"max_abs_identity_residual_A2": worst, "rows": rows,
            "PASS": bool(worst < 1e-9)}


# ================================================== the matched-budget thermostat
def search_metropolis(E, n, k, budget, T, rng, restart_every=None):
    """THE CONTROL SPRINT 16 NEVER RAN: a fixed-temperature classical thermostat that MOVES
    SUPPORT, at MATCHED objective budget.

    Single-residue Metropolis on the same k^n register, temperature held CONSTANT (this is a
    thermostat, not an annealer -- `V.search_anneal` is the annealer and is run beside it),
    with periodic random restarts so the chain is not trapped in one basin for the whole
    budget.  Budget is enforced by the same `V.Counter` every other arm uses, so the
    comparison is at nominal parity AND at true parity: one objective read per proposal.

    Sweeping `T` traces the classical arm's whole curve on the (M, D) plane, which is what
    the Pareto question needs: a single temperature is a point, and a point cannot be a
    frontier.
    """
    c = V.Counter(E, budget)
    pw = k ** np.arange(n - 1, -1, -1)
    restart_every = restart_every or max(64, budget // 8)
    s = rng.integers(0, k, n)
    cur = int(s @ pw)
    e0 = c([cur])
    if e0.size == 0:
        return c
    cur_e = float(e0[0])
    step = 0
    T = max(float(T), 1e-12)
    while c.left > 0:
        step += 1
        if step % restart_every == 0:
            s = rng.integers(0, k, n)
            cur = int(s @ pw)
            ev = c([cur])
            if ev.size == 0:
                break
            cur_e = float(ev[0])
            continue
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


def search_metropolis_1opt(E, n, k, budget, T, rng, frac=0.5):
    """The falsifier's named classical arm: THE THERMOSTAT PLUS 1-OPT.

    `frac` of the budget on the fixed-T thermostat, the remainder on greedy 1-opt with
    restarts.  Both halves' visited sets are returned as one multiset, so the arm's (M, D)
    point is the point a consumer of the union would see.
    """
    b1 = int(round(frac * budget))
    c1 = search_metropolis(E, n, k, b1, T, rng)
    c2 = V.search_greedy(E, n, k, budget - b1, rng)
    seen = np.concatenate([c1.all_seen(), c2.all_seen()])
    used = int(c1.used + c2.used)
    return seen, used


# ============================================ the retrieval product law (0 evals)
def prior_factors(ins):
    """EXACT.  `ins.prior` is additive over residues (verified: residual 9.7e-07, the float32
    spacing).  Recover the per-residue table `f[i, a]` so its Boltzmann law is an exact
    PRODUCT distribution -- samplable in closed form at ZERO objective evaluations.
    """
    p = os.path.join(CACHE, f"priorfac_{ins.pdb}.npz")
    if os.path.exists(p):
        d = np.load(p)
        return np.asarray(d["f"], float), float(d["c"]), float(d["resid"])
    S = ins.states(np.arange(ins.N))
    pr = ins.prior
    f = np.empty((ins.n, ins.k))
    for i in range(ins.n):
        for a in range(ins.k):
            f[i, a] = pr[S[:, i] == a].mean()
    c = -(ins.n - 1) * float(pr.mean())
    resid = float(np.abs(f[np.arange(ins.n)[None, :], S].sum(1) + c - pr).max())
    np.savez(p, f=f, c=c, resid=resid)
    return f, c, resid


def prior_draws(ins, T, budget, rng):
    """Draw from `p(x) ~ exp(-prior(x)/T)`, exactly, as a product over residues.

    NATIVE-FREE and costs ZERO reads of the deployed objective -- the retrieval distribution
    itself, which the brief lists as an arm.  `T -> inf` is uniform random.
    """
    f, _, _ = prior_factors(ins)
    z = -f / max(float(T), 1e-12)
    z = z - z.max(1, keepdims=True)
    q = np.exp(z)
    q /= q.sum(1, keepdims=True)
    pw = ins.k ** np.arange(ins.n - 1, -1, -1)
    S = np.empty((budget, ins.n), np.int64)
    for i in range(ins.n):
        S[:, i] = rng.choice(ins.k, size=budget, p=q[i])
    return S @ pw


# ================================================== mean-field Boltzmann of E (P4)
def meanfield_boltzmann(ins, E, T, iters=300, tol=1e-11):
    """THE CORRECT CLASSICAL ANALOGUE OF A LOW-CORRELATION VARIATIONAL FAMILY.

    Naive mean-field: the fully factorised `q(x) = prod_i q_i(x_i)` that is a fixed point of
    `KL(q || exp(-E/T)/Z)`, by exact coordinate updates.  The register index is
    `sum_i s_i k^(n-1-i)`, so reshaping `E` to a `(k,)*n` tensor makes residue `i` axis `i`
    and every conditional marginal is one `logsumexp` over the other axes -- exact, no
    sampling, `n` passes over the register per sweep.

    A depth-2 MPS circuit is a BOUNDED-correlation variational family; this fully factorised
    law is its classical zero-correlation limit and is the control P4 requires.  If
    `q_theta` adds nothing over this, the variational family contributes factorised smoothing
    and nothing else.
    """
    from scipy.special import logsumexp
    n, k = ins.n, ins.k
    shape = (k,) * n
    negE = (-np.asarray(E, float) / max(float(T), 1e-12)).reshape(shape)
    negE = negE - negE.max()
    logq = np.zeros((n, k))
    for _ in range(iters):
        Rt = negE.copy()
        for i in range(n):
            Rt = Rt + logq[i].reshape([k if j == i else 1 for j in range(n)])
        new = np.empty((n, k))
        for i in range(n):
            Ri = Rt - logq[i].reshape([k if j == i else 1 for j in range(n)])
            li = logsumexp(Ri, axis=tuple(j for j in range(n) if j != i))
            new[i] = li - logsumexp(li)
        d = float(np.abs(np.exp(new) - np.exp(logq)).max())
        logq = new
        if d < tol:
            break
    q = np.exp(logq)
    q /= q.sum(1, keepdims=True)
    lp = np.zeros(shape)
    for i in range(n):
        lp = lp + np.log(np.maximum(q[i], 1e-300)).reshape(
            [k if j == i else 1 for j in range(n)])
    return q, lp.reshape(-1)


def meanfield_boltzmann_slow(ins, E, T, iters=200, tol=1e-10, seed=0):
    """THE CORRECT CLASSICAL ANALOGUE OF A LOW-CORRELATION VARIATIONAL FAMILY.

    Naive mean-field: the fully factorised `q(x) = prod_i q_i(x_i)` minimising
    `KL(q || exp(-E/T)/Z)`, by fixed-point iteration on the conditional means.  A
    depth-2 MPS circuit is a *bounded-correlation* variational family; the fully factorised
    Boltzmann law is its classical zero-correlation limit and is the control P4 requires --
    if `q_theta` adds nothing over this, the variational family is contributing factorised
    smoothing and nothing else.

    Exact on this instrument: the conditional expectations are computed by enumerating the
    register, so no sampling noise enters.
    """
    n, k, N = ins.n, ins.k, ins.N
    S = ins.states(np.arange(N))
    e = np.asarray(E, float) / max(float(T), 1e-12)
    e = e - e.min()
    q = np.full((n, k), 1.0 / k)
    logq = None
    for _ in range(iters):
        # log w(x) = sum_i log q_i(x_i); update each factor from the tilted marginal
        lw = np.log(np.maximum(q, 1e-300))[np.arange(n)[None, :], S].sum(1)
        newq = np.empty_like(q)
        for i in range(n):
            li = lw - np.log(np.maximum(q[i], 1e-300))[S[:, i]] - e
            for a in range(k):
                m = S[:, i] == a
                x = li[m]
                newq[i, a] = np.logaddexp.reduce(x)
            newq[i] -= newq[i].max()
            newq[i] = np.exp(newq[i]); newq[i] /= newq[i].sum()
        d = float(np.abs(newq - q).max())
        q = newq
        if d < tol:
            break
    logq = np.log(np.maximum(q, 1e-300))[np.arange(n)[None, :], S].sum(1)
    return q, logq


# =================================================================== accounting
class Cost:
    """P2. Every axis, always. Never a nominal-parity comparison without its true cost."""

    def __init__(self):
        self.t0 = time.time()
        self.c0 = time.process_time()

    def stop(self, obj_evals, circuit_samples=0, grad_passes=0, distinct=0,
             register_reads=0):
        return {"objective_evals": int(obj_evals),
                "circuit_samples": int(circuit_samples),
                "grad_passes": int(grad_passes),
                "distinct_configs": int(distinct),
                "full_register_reads": int(register_reads),
                "wall_s": float(time.time() - self.t0),
                "cpu_s": float(time.process_time() - self.c0)}


# =================================================================== statistics
def fold_of(pdb):
    return int(inst(pdb).fold)


def cluster_paired(d, folds, n_boot=4000, seed=0):
    """A FOLD-CLUSTERED bootstrap beside the target bootstrap.

    `I.paired` resamples TARGETS i.i.d. and only *reports* per-fold means; with five folds a
    true cluster bootstrap has very little power, so both are printed and neither is
    substituted for the other.  The per-fold sign vector is the honest small-n statement.
    """
    d = np.asarray(d, float)
    folds = np.asarray(folds)
    u = np.unique(folds)
    rng = np.random.default_rng(seed)
    groups = [d[folds == f] for f in u]
    bs = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, len(u), len(u))
        bs[b] = np.concatenate([groups[p] for p in pick]).mean()
    return {"n_folds": int(len(u)),
            "cluster_ci95": [float(np.percentile(bs, 2.5)),
                             float(np.percentile(bs, 97.5))],
            "per_fold_mean": {int(f): float(d[folds == f].mean()) for f in u},
            "folds_same_sign": int(max(
                sum(1 for f in u if d[folds == f].mean() < 0),
                sum(1 for f in u if d[folds == f].mean() > 0)))}


def verdict(d, folds=None, label="", floor=FLOOR, seed=0):
    """One line: mean, median, CI, W/L, the floor, and the fold-clustered interval."""
    d = np.asarray(d, float)
    p = I.paired(d, np.zeros_like(d), seed=seed)
    sig = (p["ci95"][0] > 0) or (p["ci95"][1] < 0)
    out = {"label": label, "n": p["n"], "mean": p["mean_diff"],
           "median": p["median_diff"], "ci95": p["ci95"],
           "W": p["n_better"], "L": p["n_worse"], "sig": bool(sig),
           "below_floor": bool(abs(p["mean_diff"]) <= floor)}
    if folds is not None:
        out.update(cluster_paired(d, folds, seed=seed))
    return out


# ===================================================================== dominance
def dominates(a, b, tau=TAU):
    """Pre-registered rule. Lower member error M is better; higher diversity D is better.

    `a` dominates `b` iff a is no worse on either coordinate and strictly better on one by
    more than `tau`.  `tau = 0.02 A` is a quarter of the programme's 0.08 A floor.
    """
    ok = (a["M"] <= b["M"] + tau) and (a["D"] >= b["D"] - tau)
    strict = (a["M"] < b["M"] - tau) or (a["D"] > b["D"] + tau)
    return bool(ok and strict)


def pareto_front(points, tau=TAU):
    """Indices of the non-dominated points under the pre-registered rule."""
    keep = []
    for i, p in enumerate(points):
        if not any(dominates(q, p, tau) for j, q in enumerate(points) if j != i):
            keep.append(i)
    return keep


def eps_dominance(p, classical):
    """THE PRIMARY PARETO STATISTIC: the additive epsilon-dominance indicator, in Angstroms.

    Write the plane as a two-objective MINIMISATION over `(M, -D)`.  A classical point `c`
    dominates `p` iff it is no worse on both.  The standard additive epsilon indicator is

        eps(p) = min_c max( M_c - M_p ,  D_p - D_c )

    `eps <= 0` means some classical point weakly dominates `p`.  `eps > 0` means NO classical
    point dominates `p`, and `eps` is how far the whole classical set would have to move to
    reach it -- so POSITIVE IS A QUANTUM WIN, in Angstroms, on a per-target basis.

    WHY THIS AND NOT `dominance_margin`.  The first statistic used here was
    `M_p - min{M_c : D_c >= D_p - tau}`, which is UNDEFINED when the VQE reaches a diversity
    no classical arm reaches -- and that case is not missing data, it is exactly the case
    where the VQE point cannot be dominated.  Reporting a mean over the defined cells would
    have silently deleted the only cells that could have been quantum wins.  The epsilon
    indicator is defined on every cell and has the right sign there by construction.
    """
    if not classical:
        return float("nan")
    return float(min(max(c["M"] - p["M"], p["D"] - c["D"]) for c in classical))


def dominance_margin(p, classical, tau=TAU):
    """How far OUTSIDE the classical attainable set a point sits, in Angstroms.

    For each classical point `c` that is at least as diverse as `p` (`D_c >= D_p - tau`), the
    margin is `M_p - min_c M_c`: how much worse `p` is than the BEST such classical point.
    Negative means `p` has lower member error than EVERY classical point of matched-or-better
    diversity, i.e. it is genuinely outside the classical frontier.

    NOTE, because getting this backwards is exactly the failure mode the brief names: the
    quantity is `M_p - min_c M_c`, NOT `min_c (M_p - M_c)`.  The latter is `M_p - max_c M_c`
    and measures the WORST classical point, which would score a dominated quantum arm as a
    win.  The first draft of this function had that bug and it was caught by hand-checking a
    printed table against the rule.

    This is the one scalar that answers the sprint's Pareto question per target, and it is
    signed so that a paired CI on it is meaningful.
    """
    cands = [c for c in classical if c["D"] >= p["D"] - tau]
    if not cands:
        return float("nan")
    return float(p["M"] - min(c["M"] for c in cands))
