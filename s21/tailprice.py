"""s21/tailprice.py -- PRICE THE MANDATORY MATRIX BEFORE ANY VQE RUNS.

WHY THIS EXISTS, AND WHY IT IS FIRST.  Sprint 20 measured **argmin-by-energy** over the shipped
K = 500 pool and found it worse than random for both physics energies (AMBER 4.990, Legacy 5.487,
pool mean 4.739, distogram 3.676).  The directive asks for **CVaR-VQE as the SELECTOR**.

    argmin is not CVaR.

CVaR selects an **alpha-tail** and — EXACT, Sprint 20 — its minimiser is a **face**: it constrains
the tail and is indifferent to everything outside it.  A tail MEAN and a MINIMUM are different
operators on the same energy, and on an energy that RANKS BADLY they can differ a lot: a minimum
chases the single worst-ranked outlier, a tail mean averages over many and regresses toward the
pool.

    THE QUESTION: for each Hamiltonian in the mandatory matrix, what does a TAIL selector achieve
    over the real candidate pool, and how does that compare to argmin, to random, and to the
    ORACLE ceiling?

Answering this costs minutes and **prices the whole matrix**: it bounds what any CVaR-VQE can
achieve by selection over this pool, per Hamiltonian, before a single circuit is built.  A VQE
that selects from this pool cannot beat the pool's own selector ceiling, and if a Hamiltonian's
tail is worse than random then no amount of variational skill on that Hamiltonian helps.

WHAT THIS IS NOT.  It is **not** a VQE experiment and makes no quantum claim.  It is the
selection-side ceiling that the quantum lanes must be measured against.  A real CVaR-VQE explores
a continuous torsional distribution and is not restricted to pool members; this bounds the
POOL-RESTRICTED case only, which is exactly the case the mandatory matrix runs.

ARMS, per Hamiltonian.  All native-free selection; the native is read only to score.

    argmin        the single lowest-energy candidate         Sprint 20's operator
    tail{a}       the coordinate average of the lowest-a fraction   the CVaR-tail readout
    random{m}     matched-count random subset, averaged      MATCHED CONTROL at each tail size
    ORACLE_best   the pool's best member                     the generation ceiling
    ORACLE_tail{a} the best a-fraction by TRUE RMSD          the selector ceiling at that width

HAMILTONIANS.  Every one separably evaluable, as section 7 of the brief requires:

    disto     the shipped Bayes-risk distogram score          the one that works
    legacy    genuine Legacy at DEFAULT_WEIGHTS               never fitted
    amber     genuine AMBER single point                      via the bare-single-point path
    d+l, d+a, l+a, d+l+a    z-scored sums, per target

NORMALISATION, declared before use.  Legacy and AMBER live on different scales and AMBER has a
heavy upper tail (median +16,062, max 5.5e23 on real rebuilds).  Combining raw sums would be
dominated by AMBER's outliers -- a normalisation artefact, not physics.  So every component is
**rank-transformed to a standard normal within the target's own pool** before summing.  This is
scale-free, outlier-robust, requires no fitted constant, and is declared here rather than chosen
after seeing RMSD.  Raw-sum arms are computed alongside so the choice is auditable.

PRE-REGISTRATION.

  Hypothesis:    a tail selector is materially better than argmin on the physics Hamiltonians,
                 because averaging a tail regresses toward the pool while a minimum chases the
                 worst-ranked outlier.
  Prediction:    tail(0.15) beats argmin by >0.3 A on legacy and amber; the gap shrinks as alpha
                 -> 0 and vanishes on disto, which already ranks well.
  Primary:       tail(0.15) - argmin per Hamiltonian, paired over targets, bootstrap CI.
  Falsifier:     if the tail arms do NOT beat their matched-count RANDOM subsets, tail selection
                 carries no information on that Hamiltonian and a CVaR-VQE over it is bounded by
                 chance no matter how well it optimises.
  Null:          random{m} at every matched tail size -- the operative control.
  Budget:        n = 126, no optimisation, pool-restricted.
  Promotion:     none.  This is a bound, not a method.

  THE TRAP.  A tail arm that beats argmin has not thereby beaten anything useful -- as alpha -> 1
  the tail becomes the whole-pool average, which is a zero-information operator that is known to
  be decent.  **The matched random control at the same count is the only comparison that means
  anything**, and it is reported beside every tail.

--------------------------------------------------------------------------------------------
WORKSTREAM A, 2026-09-07.  The pre-registration above is the coordinator's and is UNEDITED.
Two API calls in the body were wrong and are fixed here; the fixes changed no arm, no null and
no endpoint.  Recorded so the diff is auditable:

  1. `EL.legacy_components_of_windows(pdb, W, seq, fold)` -> `(seq, PHI, PSI)`.  The function
     takes torsions, not coordinates.  Pool torsions come from the same universe rows as W
     (`u["PHI"][pool_idx]`), which is exactly what `s14.retprior.windows(pdb,"pool")` returns.

  2. `QA.amber_energies(pdb, W, seq, fold)` -> the COORDINATE path.  `s13.qarch_lib.amber_energies`
     takes a discrete `Space` and integer states -- the lattice path -- and cannot score an
     arbitrary retrieved window.  Replaced with the genuine ff14SB/GBn2 single point used by
     `s18/phys_down.py` and certified by Sprint 20's gate GC20a:
         tab = tl2.library_for(seq,4,seq); rep = tl2.PerResidueTorsion(seq,tab,chi_bits=False)
         P17.ConstrainedBox(seq, rep).energy_point(build_backbone(phi,psi))
     No minimisation of any kind: `H_AMBER = E`, not `E o Relax`.  Sprint 20 L7c warns the
     relaxation is CONSTITUTIVE and that AMBER is not finite on 42% of the *lattice* register
     without it -- on these REAL rebuilds it is finite everywhere (L5c), and the bare point is
     the arm that isolates the energy from the repair operator.  Both facts are reported.

COORDINATOR-ISSUED AMENDMENT, 2026-09-07, logged before the table was read.  The
pre-registration above is left UNEDITED; this is the auditable diff.

  A. THE PRIMARY ENDPOINT WAS CONFOUNDED AND IS RE-PRIMARIED.  `argmin` returns ONE member;
     `tail{a}` returns the COORDINATE AVERAGE of a*m members.  They differ by the AVERAGING
     operator, which Sprint 19 priced at ~1.0 A.  `tail(0.15) - argmin` would therefore have
     "confirmed" the >0.3 A prediction on EVERY Hamiltonian including a random one, for
     entirely the wrong reason.
         NEW PRIMARY:  tail{a} - random{a} at MATCHED COUNT -- already in the file as the
         operative control, and it should have been the endpoint from the start.
     The old primary is retained and printed, relabelled as what it actually measures.

  B. THE TAIL OPERATOR AND THE AVERAGING OPERATOR ARE NOW SEPARATED.  Beside every averaged
     tail, two SINGLE-STRUCTURE readouts of the same tail:
         medoid(tail)      one real structure, no averaging
         member(tail)      a uniformly random member of the tail, selection-within-tail removed
     so that   medoid(tail) - argmin   prices the TAIL with averaging held out, and
               tail_avg   - medoid(tail)   prices the AVERAGING with the tail held fixed.
     Both have their own matched-count random controls (medoid and member of a RANDOM subset).
     EVERY AVERAGED ARM IS LABELLED AVERAGED.

  C. WHY THE SINGLE-STRUCTURE ARMS ARE THE ONES THAT CORRESPOND TO THE MATRIX.  The deployed
     pillar emits `vqe_bitstring` (argmin over the final distribution's samples),
     `vqe_modal_bitstring` (the mode) and `best_seen_bitstring` (argmin over everything seen).
     NONE is an average of a tail.  CVaR is the TRAINING objective; the READOUT is an argmin.
     The averaged arms are a different (Sprint-19) lever and are labelled as such.

  D. CROSS-READOUT CELLS ADDED, because a near-theorem makes them the only informative ones.
     CVaR's minimiser is a face supported on the alpha-tail, and the alpha-tail of H CONTAINS
     H's global pool minimum -- so for a POOL-RESTRICTED selector whose TRAINING and READOUT
     energy are the SAME H, CVaR-VQE and argmin have the same optimal answer.  The cells that
     break the scope conditions are (i) an AVERAGED readout and (ii) a READOUT H DIFFERENT
     FROM THE TRAINING H.  All 7x7 (train, read) pairs are therefore emitted at every alpha:
     take the alpha-tail of H_train, then argmin H_read INSIDE it.  Cost is nil -- no new
     geometry -- and it prices the genuinely new cell before any circuit is built.

BASIS, declared before the run and VERIFIED, not assumed:
  * the shipped top-75 filter is the distogram score on the REAL WINDOWS W (75/75 index
    agreement on 1A13; the ideal-geometry rebuild gives only 54/75).  So `disto` scores W.
  * Legacy and AMBER exist only on the ideal-geometry rebuild of each member's torsions.
    Scoring them through a geometry they never saw would be a confound, so EVERY arm is emitted
    and scored TWICE: once as the coordinate average of the selected real windows (`W` basis,
    the deployed operator) and once as the coordinate average of the selected rebuilds
    (`RB` basis).  Both are point clouds.  Neither is projected here.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402
from s18 import phys_lib as PL              # noqa: E402

ALPHAS = (0.01, 0.05, 0.15, 0.30, 0.50, 1.00)
HAMS = ("disto", "legacy", "amber", "d+l", "d+a", "l+a", "d+l+a")
COMBOS = ("d+l", "d+a", "l+a", "d+l+a")
N_RANDOM_DRAWS = 8


def _rank_normal(x):
    """Rank-transform to a standard normal within this target's pool.  Scale-free, outlier-robust,
    no fitted constant -- declared in the docstring before any RMSD was seen."""
    x = np.asarray(x, float)
    n = len(x)
    r = np.empty(n)
    r[np.argsort(x, kind="stable")] = np.arange(n)
    from scipy.special import ndtri
    return ndtri((r + 0.5) / n)


def _zscore(x):
    """Per-target z-score.  AUDIT arm only -- the declared normalisation is `_rank_normal`."""
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / s if s > 1e-12 else np.zeros_like(x)


def _stats(x):
    """The distributional facts §6 of the brief requires BEFORE combining."""
    x = np.asarray(x, float)
    f = np.isfinite(x)
    y = x[f]
    if len(y) < 4:
        return {"n": int(len(y)), "finite_frac": float(f.mean())}
    m, s = float(y.mean()), float(y.std())
    z = (y - m) / s if s > 1e-12 else np.zeros_like(y)
    q = np.percentile(y, [0, 1, 25, 50, 75, 99, 100])
    return {"n": int(len(y)), "finite_frac": float(f.mean()),
            "mean": m, "sd": s, "skew": float((z ** 3).mean()), "kurt": float((z ** 4).mean()),
            "q0": float(q[0]), "q1": float(q[1]), "q25": float(q[2]), "q50": float(q[3]),
            "q75": float(q[4]), "q99": float(q[5]), "q100": float(q[6]),
            "frac_gt_1e4": float((y > 1e4).mean()), "range_decades":
                float(np.log10(max(abs(q[6] - q[0]), 1e-30)))}


def components(pdb, seq, fold, n, W, PHI, PSI, box=None):
    """Per-candidate energies for every separable component.  Native-free.

    `disto` is scored on the REAL windows W -- verified to reproduce the shipped top-75 filter
    exactly.  `legacy` and `amber` are scored on the ideal-geometry rebuild of the same member's
    torsions, which is the only geometry on which either energy is defined; `disto_rb` is the
    distogram score on that same rebuild, kept so the basis asymmetry is auditable.
    """
    from core import geometry as geo
    from s16 import energy_lib as EL

    out = {}
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    out["disto"] = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)

    bb = geo.build_backbone_batch(PHI, PSI)
    CAb = np.asarray(bb["CA"], float)
    ii, jj = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
    D = np.linalg.norm(CAb[:, ii, :] - CAb[:, jj, :], axis=-1)
    out["disto_rb"] = np.asarray(I.shipped_score(dg, D), float)

    comp = EL.legacy_components_of_windows(seq, PHI, PSI)
    out["legacy"] = np.asarray(EL.legacy_total_from(comp), float)

    #: GENUINE ff14SB/GBn2 single point.  No minimisation: H_AMBER = E, not E o Relax.
    B = len(PHI)
    E = np.full(B, np.nan)
    for b in range(B):
        c = geo.build_backbone(PHI[b], PSI[b])
        E[b] = box.energy_point({a: np.asarray(v, float) for a, v in c.items()})["energy"]
    out["amber"] = E
    out["_ca_rebuild"] = CAb
    return out


def _make_box(seq, tries=40):
    """Build the OpenMM context, yielding to a contended box first.

    MEASURED, not assumed: this run was KILLED at target 13/126 by
    `core.amber.memory_guard` ("physical memory at 96% exceeds the 92% ceiling") because
    sibling workstreams filled the box.  The guard is correct; the caller was not resilient.
    Hold for memory, then retry with backoff, and never proceed by disabling the guard.
    """
    import time
    import torsion_lib2 as tl2
    from s16 import energy_lib as EL
    from s17 import phys_lib as P17
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    last = None
    for a in range(tries):
        EL.mem_hold(min_gb=1.45, tag="s21A/tailprice", max_wait=420.0)
        try:
            return P17.ConstrainedBox(seq, rep)
        except MemoryError as e:                                  # pragma: no cover
            last = e
            print(f"  [s21A] memory guard fired (attempt {a+1}/{tries}); backing off",
                  flush=True)
            time.sleep(30)
    raise last


def target_row(t):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    p = I.pool_idx(u)
    W = np.asarray(u["W"], float)[p]
    PHI = np.asarray(u["PHI"], float)[p]
    PSI = np.asarray(u["PSI"], float)[p]
    rng = SD.stable_rng(pdb, "s21tailprice")
    m = len(W)

    box = _make_box(seq)
    try:
        cmp_ = components(pdb, seq, fold, n, W, PHI, PSI, box=box)
    finally:
        box.close()
    RB = cmp_.pop("_ca_rebuild")

    raw = {k: v for k, v in cmp_.items() if isinstance(v, np.ndarray)}
    z = {k: _rank_normal(v) for k, v in raw.items()}
    zz = {k: _zscore(v) for k, v in raw.items()}

    def _mix(d):
        return {"disto": d["disto"], "legacy": d["legacy"], "amber": d["amber"],
                "d+l": d["disto"] + d["legacy"], "d+a": d["disto"] + d["amber"],
                "l+a": d["legacy"] + d["amber"],
                "d+l+a": d["disto"] + d["legacy"] + d["amber"]}

    E = _mix(z)                      # DECLARED normalisation: within-target rank-to-normal
    E_raw = _mix(raw)                # AUDIT: raw sums
    E_z = _mix(zz)                   # AUDIT: per-target z-score

    #: per-member TRUE RMSD on BOTH bases.  ORACLE -- scoring and ceilings only.
    true = np.asarray(I.kabsch_rmsd_batch(W, nat), float)
    true_rb = np.asarray(I.kabsch_rmsd_batch(RB, nat), float)

    e = {"m": int(m),
         "pool_mean": float(true.mean()), "ORACLE_best": float(true.min()),
         "pool_mean_RB": float(true_rb.mean()), "ORACLE_best_RB": float(true_rb.min()),
         "rebuild_shift": float(np.mean(true_rb - true)),
         "rho_true_trueRB": PL.spearman(true, true_rb),
         "rho_leg_amb": PL.spearman(raw["legacy"], raw["amber"]),
         "rho_dis_leg": PL.spearman(raw["disto"], raw["legacy"]),
         "rho_dis_amb": PL.spearman(raw["disto"], raw["amber"]),
         "ORACLE_rho_dis_d": PL.spearman(raw["disto"], true),
         "ORACLE_rho_leg_d": PL.spearman(raw["legacy"], true_rb),
         "ORACLE_rho_amb_d": PL.spearman(raw["amber"], true_rb),
         "stats": {k: _stats(v) for k, v in raw.items()}}

    memo = {}

    def readouts(idx):
        """Every readout of ONE candidate subset, on BOTH bases.  Memoised on the index set so
        the 16 whole-pool calls per target cost one.

            avg     AVERAGED -- the deployed coordinate average (a point cloud, NOT a structure)
            med     medoid of the subset -- ONE REAL STRUCTURE, no averaging
        `member` (a uniformly random element) is drawn by the caller, which needs the indices.
        """
        idx = np.sort(np.asarray(idx, int))
        key = idx.tobytes()
        if key not in memo:
            if len(idx) == 1:
                memo[key] = {"avg": float(true[idx[0]]), "avg_RB": float(true_rb[idx[0]]),
                             "med": float(true[idx[0]]), "med_RB": float(true_rb[idx[0]])}
            else:
                a, bW = I.coordinate_average(W[idx])
                a2, bR = I.coordinate_average(RB[idx])
                memo[key] = {"avg": float(I.ca_rmsd(np.asarray(a, float), nat)),
                             "avg_RB": float(I.ca_rmsd(np.asarray(a2, float), nat)),
                             "med": float(true[idx[bW]]), "med_RB": float(true_rb[idx[bR]])}
        return memo[key]

    orders = {h: np.argsort(E[h], kind="stable") for h in HAMS}

    for h in HAMS:
        order = orders[h]
        #: SINGLE STRUCTURE.  This is what the deployed pillar actually emits.
        e[f"argmin|{h}"] = float(true[order[0]])
        e[f"argmin_RB|{h}"] = float(true_rb[order[0]])
        for a in ALPHAS:
            k = max(int(round(a * m)), 1)
            T = order[:k]
            r = readouts(T)
            e[f"tail{a:g}|{h}"] = r["avg"]                      # AVERAGED
            e[f"tail{a:g}_RB|{h}"] = r["avg_RB"]                # AVERAGED
            e[f"tailmed{a:g}|{h}"] = r["med"]                   # single structure
            e[f"tailmed{a:g}_RB|{h}"] = r["med_RB"]             # single structure
            #: uniformly random member of the tail -- selection WITHIN the tail removed
            pick = rng.choice(k, size=min(N_RANDOM_DRAWS, k), replace=(k < N_RANDOM_DRAWS))
            e[f"tailmem{a:g}|{h}"] = float(np.mean(true[T[pick]]))
            e[f"tailmem{a:g}_RB|{h}"] = float(np.mean(true_rb[T[pick]]))
            #: CROSS-READOUT.  Train on h, read out by argmin of h2 INSIDE h's alpha-tail.
            #: h2 == h reproduces argmin exactly whenever the tail contains h's minimum,
            #: which it does by construction -- the near-theorem's scope condition, made visible.
            for h2 in HAMS:
                w = T[int(np.argmin(E[h2][T]))]
                e[f"x{a:g}|{h}>{h2}"] = float(true[w])
                e[f"x{a:g}_RB|{h}>{h2}"] = float(true_rb[w])

    #: MATCHED-COUNT RANDOM control at every tail size, for EVERY readout, and the ORACLE
    #: selector ceiling.  A control must be matched in the space the OPERATOR works in.
    for a in ALPHAS:
        k = max(int(round(a * m)), 1)
        acc = {"avg": [], "avg_RB": [], "med": [], "med_RB": [], "mem": [], "mem_RB": []}
        for _ in range(N_RANDOM_DRAWS):
            S = rng.choice(m, k, replace=False)
            r = readouts(S)
            for q in ("avg", "avg_RB", "med", "med_RB"):
                acc[q].append(r[q])
            j = S[rng.integers(0, k)]
            acc["mem"].append(float(true[j])); acc["mem_RB"].append(float(true_rb[j]))
        e[f"random{a:g}"] = float(np.mean(acc["avg"]))
        e[f"random{a:g}_RB"] = float(np.mean(acc["avg_RB"]))
        e[f"random{a:g}_sd"] = float(np.std(acc["avg"]))
        e[f"randommed{a:g}"] = float(np.mean(acc["med"]))
        e[f"randommed{a:g}_RB"] = float(np.mean(acc["med_RB"]))
        e[f"randommem{a:g}"] = float(np.mean(acc["mem"]))
        e[f"randommem{a:g}_RB"] = float(np.mean(acc["mem_RB"]))
        oW = readouts(np.argsort(true, kind="stable")[:k])
        oR = readouts(np.argsort(true_rb, kind="stable")[:k])
        e[f"ORACLEtail{a:g}"] = oW["avg"]
        e[f"ORACLEtail{a:g}_RB"] = oR["avg_RB"]
        e[f"ORACLEtailmed{a:g}"] = oW["med"]
        e[f"ORACLEtailmed{a:g}_RB"] = oR["med_RB"]

    #: AUDIT normalisations at the primary widths only.  Reported whatever they say.
    for tag, EE in (("raw", E_raw), ("zs", E_z)):
        for h in COMBOS:
            order = np.argsort(EE[h], kind="stable")
            e[f"AUD{tag}|argmin|{h}"] = float(true[order[0]])
            for a in (0.05, 0.15, 0.30):
                k = max(int(round(a * m)), 1)
                rr = readouts(order[:k])
                e[f"AUD{tag}|tail{a:g}|{h}"] = rr["avg"]
                e[f"AUD{tag}|tailmed{a:g}|{h}"] = rr["med"]

    return {"pdb": pdb, "n": n, "fold": fold, **e}


#: A row is COMPLETE only if every Hamiltonian x alpha x readout cell it promises is present.
#: The flag must require the FULL configuration, not the subset that happened to be called.
REQUIRED_KEYS = (["pool_mean", "pool_mean_RB", "ORACLE_best", "ORACLE_best_RB", "m", "stats",
                  "rho_leg_amb", "rebuild_shift"]
                 + [f"argmin{s}|{h}" for h in HAMS for s in ("", "_RB")]
                 + [f"{p}{a:g}{s}|{h}" for h in HAMS for a in ALPHAS for s in ("", "_RB")
                    for p in ("tail", "tailmed", "tailmem")]
                 + [f"x{a:g}{s}|{h}>{h2}" for h in HAMS for h2 in HAMS for a in ALPHAS
                    for s in ("", "_RB")]
                 + [f"{p}{a:g}{s}" for a in ALPHAS for s in ("", "_RB")
                    for p in ("random", "randommed", "randommem", "ORACLEtail", "ORACLEtailmed")]
                 + [f"AUD{tag}|tail{a:g}|{h}" for tag in ("raw", "zs") for h in COMBOS
                    for a in (0.05, 0.15, 0.30)])


def row_complete(r):
    if r.get("skipped"):
        return False
    return all(k in r for k in REQUIRED_KEYS)


PATH = os.path.join(RESULTS, "tailprice.json")
CFG = {"ALPHAS": list(ALPHAS), "HAMS": list(HAMS), "N_RANDOM_DRAWS": N_RANDOM_DRAWS,
       "pool": "shipped K=500 BLOSUM pool", "norm": "within-target rank-to-normal",
       "amber": "genuine ff14SB/GBn2 single point, NO minimisation (k=0, steps=-1 equivalent)",
       "basis": "W = coordinate average of real retrieved windows; RB = of ideal-geometry rebuilds"}


def _write(rows, n_expected):
    #: HAZARD, named in BRIEF section 6 and hit once already in this sprint: a completion flag
    #: that counts SKIPPED rows will certify a run of total failures as complete.  This one
    #: requires n_expected rows, ZERO skips, and every Hamiltonian x alpha x readout cell.
    n_skip = sum(1 for r in rows if r.get("skipped"))
    n_ok = sum(1 for r in rows if row_complete(r))
    obj = {"rows": rows, "config": CFG, "cfg_hash": PL.cfg_hash(CFG),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "n_skipped": n_skip, "n_complete_rows": n_ok,
           "n_required_keys": len(REQUIRED_KEYS),
           "complete": bool(n_skip == 0 and n_ok == int(n_expected)
                            and len(rows) == int(n_expected))}
    tmp = PATH + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, PATH)


def run(targets=None, out=None):
    global PATH
    if out:
        PATH = os.path.join(RESULTS, out)
    tg = targets if targets is not None else I.targets()
    rows = []
    if os.path.exists(PATH):
        try:
            rows = json.load(open(PATH)).get("rows", [])
        except Exception:
            rows = []
    #: RESUME ONLY FROM ROWS THAT ARE THEMSELVES COMPLETE.  A stale artefact from an earlier,
    #: broken version of this file wrote 3 SKIPPED rows; resuming naively would have silently
    #: dropped those 3 targets from n=126.  Rows failing `row_complete` are discarded and recomputed.
    n_before = len(rows)
    rows = [r for r in rows if row_complete(r)]
    if len(rows) != n_before:
        print(f"  RESUME: discarded {n_before - len(rows)} incomplete/skipped rows", flush=True)
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for c, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(target_row(t))
        _write(rows, len(tg))
        if len(rows) % 5 == 0:
            print(f"  {len(rows)}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    _write(rows, len(tg))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    report_all()
    return rows


def _pair(a, b, folds, names=None):
    p = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds, names=names)
    #: MDE IS PER COMPARISON, not a project constant.  The 0.084 A figure is the instrument's
    #: pooled value and is wrong by up to 84x in BOTH directions on an individual contrast, so
    #: every line carries its own SE and the MDE that its own paired sd implies at 80% power.
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[np.isfinite(d)]
    p["se"] = float(d.std(ddof=1) / max(np.sqrt(len(d)), 1e-12)) if len(d) > 1 else float("nan")
    p["mde"] = 2.8016 * p["se"]
    return p


def _line(p):
    ci, cf = p["ci"], p.get("ci_fold", p["ci"])
    return (f"{p['mean']:+.3f} med {p['median']:+.3f} se {p['se']:.3f} "
            f"(MDE {p['mde']:.3f})  iid[{ci[0]:+.3f},{ci[1]:+.3f}]  "
            f"fold[{cf[0]:+.3f},{cf[1]:+.3f}]  {p['W']}W/{p['L']}L")


def _folds_sign(p):
    if "per_fold" not in p:
        return ""
    v = list(p["per_fold"].values())
    neg = sum(1 for x in v if x < 0)
    return f"  folds {neg}/{len(v)} neg"


def report(rows=None, suffix="W"):
    if rows is None:
        rows = json.load(open(PATH))["rows"]
    rows = [r for r in rows if not r.get("skipped")]
    folds = np.array([r["fold"] for r in rows])
    names = [r["pdb"] for r in rows]
    sfx = "" if suffix == "W" else "_RB"
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731

    print(f"\n{'='*94}\nTAIL PRICING -- basis {suffix}  "
          f"({'coordinate average of REAL retrieved windows' if suffix=='W' else 'coordinate average of ideal-geometry REBUILDS'})")
    print(f"n = {len(rows)}.  POOL-RESTRICTED selection over the shipped K=500 pool.  "
          f"POINT-CLOUD basis, unprojected.\n{'='*94}")
    pm = g("pool_mean" + sfx).mean()
    print(f"  pool mean {pm:.3f}   ORACLE pool best {g('ORACLE_best'+sfx).mean():.3f}"
          f"   (median pool mean {np.median(g('pool_mean'+sfx)):.3f})\n")
    print(f"  {'alpha':<8}{'k':>5}{'random':>9}{'ORACLE':>9}   " + "".join(f"{h:>10}" for h in HAMS))
    kk = np.median(g("m"))
    print(f"  {'argmin':<8}{1:>5}{'--':>9}{'--':>9}   " +
          "".join(f"{g(f'argmin{sfx}|'+h).mean():>10.3f}" for h in HAMS))
    for a in ALPHAS:
        rk = f"random{a:g}" + ("" if sfx == "" else "_RB")
        ok = f"ORACLEtail{a:g}" + sfx
        print(f"  {a:<8g}{int(round(a*kk)):>5}{g(rk).mean():>9.3f}{g(ok).mean():>9.3f}   " +
              "".join(f"{g(f'tail{a:g}{sfx}|'+h).mean():>10.3f}" for h in HAMS))

    print("\n  MEDIANS over targets (same table):")
    print(f"  {'alpha':<8}{'random':>9}{'ORACLE':>9}   " + "".join(f"{h:>10}" for h in HAMS))
    print(f"  {'argmin':<8}{'--':>9}{'--':>9}   " +
          "".join(f"{np.median(g(f'argmin{sfx}|'+h)):>10.3f}" for h in HAMS))
    for a in ALPHAS:
        rk = f"random{a:g}" + ("" if sfx == "" else "_RB")
        print(f"  {a:<8g}{np.median(g(rk)):>9.3f}{np.median(g(f'ORACLEtail{a:g}'+sfx)):>9.3f}   " +
              "".join(f"{np.median(g(f'tail{a:g}{sfx}|'+h)):>10.3f}" for h in HAMS))

    print("\n  SINGLE-STRUCTURE readouts -- MEDOID of the alpha-tail (no averaging).")
    print("  This is the family the deployed pillar actually emits.  Means over targets.")
    print(f"  {'alpha':<8}{'randmed':>9}{'ORACLE':>9}   " + "".join(f"{h:>10}" for h in HAMS))
    for a in ALPHAS:
        print(f"  {a:<8g}{g(f'randommed{a:g}'+sfx).mean():>9.3f}"
              f"{g(f'ORACLEtailmed{a:g}'+sfx).mean():>9.3f}   " +
              "".join(f"{g(f'tailmed{a:g}{sfx}|'+h).mean():>10.3f}" for h in HAMS))

    print("\n  SINGLE-STRUCTURE readouts -- RANDOM MEMBER of the alpha-tail "
          "(selection WITHIN the tail removed).")
    print(f"  {'alpha':<8}{'randmem':>9}{'':>9}   " + "".join(f"{h:>10}" for h in HAMS))
    for a in ALPHAS:
        print(f"  {a:<8g}{g(f'randommem{a:g}'+sfx).mean():>9.3f}{'':>9}   " +
              "".join(f"{g(f'tailmem{a:g}{sfx}|'+h).mean():>10.3f}" for h in HAMS))

    print(f"\n{'-'*94}")
    print("  ***  PRIMARY (amended) -- ARM minus its MATCHED-COUNT RANDOM control.  ***")
    print("  Negative = the Hamiltonian's alpha-tail carries selection information over chance.")
    print("  Each readout is compared to the SAME readout of a random subset of the same size,")
    print("  because a control must be matched in the space the operator works in.")
    print(f"{'-'*94}")
    for lab, pre, rnd in (("AVERAGED tail (Sprint-19 lever, NOT the pillar's readout)", "tail", "random"),
                          ("MEDOID of tail   (single structure -- the pillar's family)", "tailmed", "randommed"),
                          ("MEMBER of tail   (single structure, no within-tail choice)", "tailmem", "randommem")):
        print(f"\n  == {lab}")
        for a in (0.01, 0.05, 0.15, 0.30, 0.50):
            print(f"    alpha={a:g}  (k={int(round(a*kk))})")
            for h in HAMS:
                p = _pair(g(f"{pre}{a:g}{sfx}|{h}"), g(f"{rnd}{a:g}{sfx}"), folds, names)
                lo, hi = p["ci"]
                flag = ("" if lo < 0 < hi else
                        ("  <-- BEATS random" if p["mean"] < 0 else "  <-- WORSE than random"))
                print(f"      {h:<8}{_line(p)}{_folds_sign(p)}{flag}")

    print(f"\n{'-'*94}")
    print("  OPERATOR DECOMPOSITION -- what the two operators are each worth, held apart.")
    print(f"{'-'*94}")
    print("  (i) THE TAIL, averaging held out:  medoid(tail) - argmin")
    for a in (0.05, 0.15, 0.30):
        print(f"    alpha={a:g}")
        for h in HAMS:
            p = _pair(g(f"tailmed{a:g}{sfx}|{h}"), g(f"argmin{sfx}|{h}"), folds, names)
            print(f"      {h:<8}{_line(p)}{_folds_sign(p)}")
    print("\n  (ii) THE AVERAGING, tail held fixed:  tail_avg - medoid(tail)")
    for a in (0.05, 0.15, 0.30):
        print(f"    alpha={a:g}")
        for h in HAMS:
            p = _pair(g(f"tail{a:g}{sfx}|{h}"), g(f"tailmed{a:g}{sfx}|{h}"), folds, names)
            print(f"      {h:<8}{_line(p)}{_folds_sign(p)}")
    print("\n  (iii) the CONFOUNDED original primary, retained and relabelled:")
    print("        tail_avg(0.15) - argmin  ==  (i) + (ii).  It is a sum of two operators.")
    for h in HAMS:
        p = _pair(g(f"tail0.15{sfx}|{h}"), g(f"argmin{sfx}|{h}"), folds, names)
        print(f"      {h:<8}{_line(p)}{_folds_sign(p)}")

    print("\n  vs the STRUCTURAL-ONLY reference (disto at the same alpha) -- does physics add?")
    for pre in ("tail", "tailmed"):
        print(f"    readout = {pre}")
        for a in (0.05, 0.15):
            print(f"      alpha={a:g}")
            for h in HAMS:
                if h == "disto":
                    continue
                p = _pair(g(f"{pre}{a:g}{sfx}|{h}"), g(f"{pre}{a:g}{sfx}|disto"), folds, names)
                print(f"        {h:<8}{_line(p)}{_folds_sign(p)}")


def report_cross(rows=None, suffix="W"):
    """The cells that break the near-theorem's scope: READOUT H != TRAINING H.

    Take the alpha-tail of H_train, then argmin H_read inside it.  Single structure throughout,
    so the averaging operator is entirely absent.  The diagonal must reproduce argmin exactly
    whenever the tail contains the training energy's own minimum -- which it does by
    construction, and that identity is printed as a check, never as a discovery.
    """
    if rows is None:
        rows = json.load(open(PATH))["rows"]
    rows = [r for r in rows if not r.get("skipped")]
    folds = np.array([r["fold"] for r in rows])
    sfx = "" if suffix == "W" else "_RB"
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    print(f"\n{'='*94}\nCROSS-READOUT MATRIX -- basis {suffix}.  "
          f"rows = TRAINING H (defines the tail); cols = READOUT H (argmin inside it).")
    print("Single structure, no averaging anywhere.  Mean Ca-RMSD over targets.")
    print(f"{'='*94}")
    for a in (0.05, 0.15, 0.30, 0.50):
        print(f"\n  alpha = {a:g}")
        print(f"    {'train\\read':<12}" + "".join(f"{h:>10}" for h in HAMS) + f"{'argminTr':>11}")
        for h in HAMS:
            row = "".join(f"{g(f'x{a:g}{sfx}|{h}>{h2}').mean():>10.3f}" for h2 in HAMS)
            print(f"    {h:<12}{row}{g(f'argmin{sfx}|{h}').mean():>11.3f}")
        d = max(abs(g(f"x{a:g}{sfx}|{h}>{h}").mean() - g(f"argmin{sfx}|{h}").mean())
                for h in HAMS)
        print(f"    IDENTITY CHECK (not a discovery): max |diag - argmin(train)| = {d:.2e}")
    print("\n  Best OFF-DIAGONAL cell vs the best DIAGONAL cell, paired over targets:")
    for a in (0.05, 0.15, 0.30):
        best_d = min(HAMS, key=lambda h: g(f"x{a:g}{sfx}|{h}>{h}").mean())
        offs = [(h, h2) for h in HAMS for h2 in HAMS if h != h2]
        best_o = min(offs, key=lambda p: g(f"x{a:g}{sfx}|{p[0]}>{p[1]}").mean())
        p = _pair(g(f"x{a:g}{sfx}|{best_o[0]}>{best_o[1]}"),
                  g(f"x{a:g}{sfx}|{best_d}>{best_d}"), folds)
        print(f"    alpha={a:g}  best off-diag {best_o[0]}>{best_o[1]} vs best diag "
              f"{best_d}>{best_d}: {_line(p)}")
        print("      (both chosen ON this instrument -- a max over 42 and 7 cells, "
              "reported as a CEILING, not a validated method)")


def report_aux(rows=None):
    if rows is None:
        rows = json.load(open(PATH))["rows"]
    rows = [r for r in rows if not r.get("skipped")]
    folds = np.array([r["fold"] for r in rows])
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731

    print(f"\n{'='*94}\nSTAGE 1 -- DISTRIBUTIONAL FACTS, measured BEFORE any combination")
    print(f"{'='*94}")
    keys = ("disto", "legacy", "amber")
    fields = ("mean", "sd", "skew", "kurt", "q0", "q25", "q50", "q75", "q99", "q100",
              "frac_gt_1e4", "finite_frac")
    print(f"  {'':<9}" + "".join(f"{f:>12}" for f in fields))
    for k in keys:
        v = {f: np.array([r["stats"][k][f] for r in rows], float) for f in fields}
        print(f"  {k:<9}" + "".join(
            (f"{np.median(v[f]):>12.4g}") for f in fields))
    print("  (per-target medians of the per-target statistic.  AMBER's q100 is the tail §6 warns about.)")

    print("\n  RANK AGREEMENT within each pool (per-target Spearman, median over targets):")
    for k in ("rho_leg_amb", "rho_dis_leg", "rho_dis_amb", "rho_true_trueRB"):
        v = g(k)
        print(f"    {k:<18}{np.median(v):+.4f}   mean {np.mean(v):+.4f}")
    #: LEGEND CORRECTED 2026-09-07 (coordinator).  It read "negative = lower energy is better",
    #: which is BACKWARDS: `disto`, the score that works, is +0.633 with only 2/35 targets
    #: negative.  As printed, the legend said the arm that works is the arm that fails.
    print("\n  ORACLE MARGINAL rank skill of each energy against TRUE RMSD")
    print("  (POSITIVE = lower energy goes with lower RMSD = CORRECT ordering):")
    for k in ("ORACLE_rho_dis_d", "ORACLE_rho_leg_d", "ORACLE_rho_amb_d"):
        v = g(k)
        p = PL.paired(v, np.zeros_like(v), folds=folds)
        cf = p.get("ci_fold", p["ci"])
        print(f"    {k:<20}{np.mean(v):+.4f} med {np.median(v):+.4f} "
              f"iid[{p['ci'][0]:+.4f},{p['ci'][1]:+.4f}] fold[{cf[0]:+.4f},{cf[1]:+.4f}]  "
              f"{int((v<0).sum())}/{len(v)} negative")

    #: ================== THE PARTIAL, which is the number that counts =====================
    #: A MARGINAL correlation with the truth, for a score already correlated with the score in
    #: PRODUCTION, is not that score's CONTRIBUTION.  Legacy's marginal (+0.33) sits beside
    #: rho(disto, legacy) = +0.48-0.65 in this same artefact.  Partialling the deployed selector
    #: out costs nothing -- all three coefficients are already recorded per target.
    #: Same class of error as Sprint 20 Q11 (every circuit-side metric collapsed once target
    #: difficulty was partialled out) and as the shared-referent floor.
    print("\n  ORACLE PARTIAL rank skill, GIVEN THE DEPLOYED DISTOGRAM -- the contribution:")
    rdt, rlt, rat = g("ORACLE_rho_dis_d"), g("ORACLE_rho_leg_d"), g("ORACLE_rho_amb_d")
    rdl, rda = g("rho_dis_leg"), g("rho_dis_amb")
    for nm, rxt, rdx in (("legacy | disto", rlt, rdl), ("amber  | disto", rat, rda)):
        den = np.sqrt(np.maximum((1 - rdx ** 2) * (1 - rdt ** 2), 1e-12))
        part = (rxt - rdx * rdt) / den
        p = PL.paired(part, np.zeros_like(part), folds=folds)
        cf = p.get("ci_fold", p["ci"])
        print(f"    {nm:<16}{np.mean(part):+.4f} med {np.median(part):+.4f} "
              f"iid[{p['ci'][0]:+.4f},{p['ci'][1]:+.4f}] fold[{cf[0]:+.4f},{cf[1]:+.4f}]  "
              f"{int((part<0).sum())}/{len(part)} negative")
    print("    (the AMBER row is the CONTROL: its marginal is already ~0, so a partial that")
    print("     also lands ~0 shows the partialling is not manufacturing structure.)")
    print("    CONSEQUENCE, already visible in the tables below: d+l never beats d, at any")
    print("    alpha or any readout.  This is the mechanism for that.")
    print(f"\n  BASIS SHIFT rebuild - window, per-member mean: {g('rebuild_shift').mean():+.4f} A")

    print(f"\n{'='*94}\nAUDIT -- the two normalisations that were NOT declared (raw sum, per-target z)")
    print(f"{'='*94}")
    print(f"  {'ham':<8}{'declared rank':>16}{'raw sum':>12}{'z-score':>12}   (tail 0.15, W basis, mean A)")
    for h in COMBOS:
        print(f"  {h:<8}{g(f'tail0.15|{h}').mean():>16.3f}"
              f"{g(f'AUDraw|tail0.15|{h}').mean():>12.3f}"
              f"{g(f'AUDzs|tail0.15|{h}').mean():>12.3f}")
    print("  If a non-declared arm wins, the declared choice is recorded as WRONG, not swapped.")


def read_verdict():
    print("\nREAD.  A tail beating argmin proves nothing on its own -- as alpha -> 1 the tail")
    print("becomes the whole-pool average, a zero-information operator that is already decent.")
    print("The matched-count RANDOM control is the only comparison that means anything, and any")
    print("Hamiltonian whose tail does not beat it bounds a CVaR-VQE over it at chance.")


def report_all():
    o = json.load(open(PATH))
    rows = o["rows"]
    print(f"ARTEFACT {PATH}\n  complete={o.get('complete')}  n_rows={o.get('n_rows')}"
          f"  n_expected={o.get('n_expected')}  n_skipped={o.get('n_skipped')}"
          f"  n_complete_rows={o.get('n_complete_rows')}"
          f"  required_keys/row={o.get('n_required_keys')}  cfg={o.get('cfg_hash')}")
    report_aux(rows)
    report(rows, "W")
    report_cross(rows, "W")
    report(rows, "RB")
    report_cross(rows, "RB")
    read_verdict()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report_all()
    elif len(sys.argv) > 1 and sys.argv[1] == "smoke":
        PATH = os.path.join(RESULTS, "_SMOKE_tailprice.json")
        run(targets=I.targets()[:3], out="_SMOKE_tailprice.json")
    else:
        run()
