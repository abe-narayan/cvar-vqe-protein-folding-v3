"""s21/d_cvarop.py -- D2: WHAT OPERATOR IS THE "CVaR TAIL" ARM ACTUALLY MEASURING?

    python -m s21.d_cvarop run | report

WHAT IS UNDER ATTACK.  `s21/BRIEF.md` section 2 and `s21/tailprice.py`: *"argmin is not CVaR ...
a tail MEAN and a MINIMUM are different operators on the same energy ... whether that difference
helps is the sprint's first measurable question."*  The pre-registered PRIMARY endpoint there is
`tail(0.15) - argmin` per Hamiltonian.

DERIVE THE OPERATOR BEFORE INTERPRETING ITS STATISTIC (BRIEF section 9).

(1) THE READOUT.  From source, `core/quantum.py:1600-1626`, a VQE run emits `vqe_bitstring`
    (argmin over the final distribution's samples), `vqe_modal_bitstring` (the mode) and
    `best_seen_bitstring` (argmin over everything seen).  **None of them is a tail mean.**  CVaR
    is the TRAINING OBJECTIVE; the READOUT is an argmin or a mode.  `tailprice.py`'s `tail{a}`
    arm is the COORDINATE AVERAGE of the alpha-tail's members, which is neither.

(2) THE POOL-RESTRICTED IDENTITY.  Q6 is EXACT: every CVaR-optimal law is supported inside the
    alpha-tail `T_a(H)`.  On a finite pool `P`, `argmin_P H` is a member of `T_a(H)` for every
    `a > 0`.  So if the readout is `argmin` over the law's support and the readout energy is the
    SAME `H` the law was trained on, the best attainable readout under a CVaR-trained law is
    exactly `argmin_P H`.  **A pool-restricted CVaR-VQE with an order-based readout and a
    single Hamiltonian cannot select anything argmin cannot.**  Scope, stated because the claim
    is worthless without it: it needs (i) pool restriction, (ii) the same H for training and
    readout, (iii) an order-based readout.  It FAILS if the readout is an average, if the
    readout energy differs from the training energy, or if the law is not converged -- and the
    UNCONVERGED case is the interesting one, which is why the measurement below prices a
    uniformly-random member of the tail, i.e. the worst case for a converged law.

(3) THE CONFOUND IN THE PRIMARY ENDPOINT.  `argmin` returns ONE structure; `tail{a}` returns the
    COORDINATE AVERAGE of `a*m` structures.  Those differ by the AVERAGING operator, which
    project memory `averaging-space-beats-the-objective` prices at ~1.0 A against torsion
    averaging while "any objective contributes only 0.171 A".  So `tail - argmin` is expected to
    be large and negative for ANY energy, including a meaningless one.

THE MEASUREMENT.  Four readouts on the same alpha-tail, so that "tail vs min" and "average vs
single structure" are separated instead of confounded, each beside its matched-count control:

    argmin(H)         single structure, lowest H            Sprint 20's operator
    tail_member(H)    a UNIFORMLY RANDOM member of the tail what a CVaR law + one-shot readout
                                                            gives in the worst (uniform) case
    tail_medoid(H)    the tail's medoid                     single structure, set-informed
    tail_avg(H)       the tail's coordinate average         `tailprice.py`'s arm

    rand_member       a uniformly random POOL member        control for tail_member and argmin
    rand_medoid{a}    medoid of a random size-k subset      control for tail_medoid
    rand_avg{a}       coordinate average of a random subset control for tail_avg

HAMILTONIANS.  `disto` (the shipped Bayes-risk score, the one that works), `legacy` (genuine
Legacy at DEFAULT_WEIGHTS, never fitted, on the ideal-geometry rebuild of each window's own
torsions -- the deployed construction), and **`rand`, a per-target random permutation: a
ZERO-INFORMATION energy.**  `rand` is the point of the experiment.  Genuine AMBER is NOT run
here; the operator decomposition does not depend on which H is used, and this lane is not
spending 63,000 continuous single points to make a point about an operator.

BASIS, stated (BRIEF section 1).  Every structural number is a **BUILT CHAIN** basis in the
weak sense that pool members are REAL retrieved windows with valid geometry; the `*_avg` arms
are **POINT CLOUDS** (a coordinate average is not a buildable backbone -- s20 X1) and are
labelled as such in the report.  Comparing an `_avg` number to a single-member number is a
basis mismatch and the report says so on the row.

PRE-REGISTRATION: `s21/PREREG_D.md` D2.
  Prediction: the zero-information `rand` energy reproduces most of the physics Hamiltonians'
              `tail(0.15) - argmin` gap.
  Falsifier:  if `rand`'s gap is LESS THAN HALF the physics gaps, the endpoint is not dominated
              by averaging and D2b is REFUTED.
  n = 126, no optimisation, pool-restricted, native read only to score.
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

from s12 import instrument as I      # noqa: E402
from s15 import seed as SD           # noqa: E402

ALPHAS = (0.05, 0.15, 0.30)
HAMS = ("disto", "legacy", "rand")
R_DRAW = 16          # repeats for every stochastic readout / control
CFG = {"alphas": list(ALPHAS), "hams": list(HAMS), "R_DRAW": R_DRAW, "K": I.K}


def energies(u, pool, seq, pdb, fold, n, rng):
    """Native-free per-candidate energies, plus one ZERO-INFORMATION energy."""
    W = np.asarray(u["W"], float)[pool]
    out = {}
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    out["disto"] = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
    from s16 import energy_lib as EL
    PHI = np.asarray(u["PHI"], float)[pool]
    PSI = np.asarray(u["PSI"], float)[pool]
    comp = EL.legacy_components_of_windows(seq, PHI, PSI)
    out["legacy"] = np.asarray(EL.legacy_total_from(comp), float)
    #: the zero-information control: a random permutation.  Not "uniform noise added to a real
    #: score" -- a permutation, so its marginal distribution is irrelevant and only its ORDER,
    #: which carries nothing, enters any of the four readouts.
    out["rand"] = rng.permutation(len(W)).astype(float)
    return out, W


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        pool = I.pool_idx(u)
        rng = SD.stable_rng(pdb, "s21_d_cvarop")
        E, W = energies(u, pool, seq, pdb, fold, n, rng)
        m = len(W)
        P = I.pairwise_rmsd(W)                       # ONE 500x500, reused by every readout
        true = I.kabsch_rmsd_batch(W, nat)           # ORACLE, scoring only

        def avg_of(idx):
            idx = np.asarray(idx, int)
            if len(idx) == 1:
                return float(true[idx[0]])
            a, _b = I.coordinate_average(W[idx], P[np.ix_(idx, idx)])
            return float(I.ca_rmsd(np.asarray(a, float), nat))

        def medoid_of(idx):
            idx = np.asarray(idx, int)
            return float(true[idx[int(I.medoid(P[np.ix_(idx, idx)]))]])

        e = {"pdb": pdb, "n": n, "fold": fold, "m": m,
             "pool_mean": float(true.mean()), "pool_median": float(np.median(true)),
             "ORACLE_pool_best": float(true.min())}
        def readouts(order):
            """Every readout for one energy ORDER.  Only the order enters, by construction."""
            r = {"argmin": float(true[order[0]])}
            for a in ALPHAS:
                k = max(int(round(a * m)), 1)
                T = order[:k]
                r[f"tail_avg{a:g}"] = avg_of(T)
                r[f"tail_medoid{a:g}"] = medoid_of(T)
                r[f"tail_member{a:g}"] = float(np.mean(
                    [true[T[rng.integers(0, k)]] for _ in range(R_DRAW)]))
                r[f"tail_min{a:g}"] = float(true[T[0]])      # == argmin, an IDENTITY check
            return r

        for h in HAMS:
            if h == "rand":
                #: R_DRAW INDEPENDENT permutations, averaged -- a single permutation would give
                #: a one-draw `argmin` whose variance swamps the comparison it is the control for
                reps = [readouts(rng.permutation(m)) for _ in range(R_DRAW)]
                rr = {k: float(np.mean([x[k] for x in reps])) for k in reps[0]}
            else:
                rr = readouts(np.argsort(E[h], kind="stable"))
            for k_, v_ in rr.items():
                e[f"{k_}|{h}"] = v_
        #: matched-count random controls, one set shared by every Hamiltonian
        e["rand_member"] = float(np.mean([true[rng.integers(0, m)] for _ in range(R_DRAW)]))
        for a in ALPHAS:
            k = max(int(round(a * m)), 1)
            S = [rng.choice(m, k, replace=False) for _ in range(R_DRAW)]
            e[f"rand_avg{a:g}"] = float(np.mean([avg_of(s) for s in S]))
            e[f"rand_medoid{a:g}"] = float(np.mean([medoid_of(s) for s in S]))
            #: the ORACLE selector ceiling at that width, labelled
            o = np.argsort(true, kind="stable")[:k]
            e[f"ORACLE_tail_avg{a:g}"] = avg_of(o)
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            _save(rows, len(tg))
    _save(rows, len(tg))
    return rows


def _save(rows, n_expected):
    json.dump({"rows": rows, "cfg": CFG, "n_expected": int(n_expected)},
              open(os.path.join(RESULTS, "d_cvarop.json"), "w"))
    #: COMPLETION FLAG requiring the FULL configuration -- every target, every Hamiltonian,
    #: every alpha, every readout, all finite.  Not the subset the run was called with.
    need = ([f"argmin|{h}" for h in HAMS]
            + [f"{r}{a:g}|{h}" for h in HAMS for a in ALPHAS
               for r in ("tail_avg", "tail_medoid", "tail_member", "tail_min")]
            + [f"{r}{a:g}" for a in ALPHAS for r in ("rand_avg", "rand_medoid",
                                                     "ORACLE_tail_avg")]
            + ["rand_member", "pool_mean"])
    ok = (len(rows) == n_expected == 126
          and all(np.isfinite(r.get(k, np.nan)) for r in rows for k in need))
    p = os.path.join(RESULTS, "d_cvarop.COMPLETE")
    if ok:
        with open(p, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"n=126 hams={list(HAMS)} alphas={list(ALPHAS)} R_DRAW={R_DRAW} "
                     f"readouts=argmin,tail_min,tail_member,tail_medoid,tail_avg "
                     f"controls=rand_member,rand_medoid,rand_avg keys={len(need)}\n")
    elif os.path.exists(p):
        os.remove(p)


def _boot(d, rng, B=20000):
    d = np.asarray(d, float); k = len(d)
    mm = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(mm, 2.5)), float(np.percentile(mm, 97.5))


def _fold_boot(d, folds, rng, B=20000):
    """FOLD-CLUSTERED interval.  `s12.instrument.paired` is i.i.d. over targets across 293 call
    sites; where a disposition could depend on the construction, quote this one too."""
    d = np.asarray(d, float); folds = np.asarray(folds)
    fs = np.unique(folds)
    groups = [d[folds == f] for f in fs]
    out = np.empty(B)
    for b in range(B):
        pick = [groups[i] for i in rng.integers(0, len(groups), len(groups))]
        out[b] = np.concatenate(pick).mean()
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def report():
    d = json.load(open(os.path.join(RESULTS, "d_cvarop.json")))
    rows = d["rows"]
    rng = np.random.default_rng(20210907)
    g = lambda k: np.array([r[k] for r in rows])          # noqa: E731
    folds = g("fold")
    n = len(rows)
    print(f"\n=== D2  THE READOUT OPERATOR, n = {n} targets, pool-restricted (K = {I.K}) ===")
    print("Native read only to score.  ORACLE rows labelled.  CI = paired target-level")
    print("bootstrap, i.i.d. over targets, with a FOLD-CLUSTERED interval beside the primary.\n")
    print(f"  pool mean {g('pool_mean').mean():.3f}   pool median {g('pool_median').mean():.3f}"
          f"   ORACLE pool best {g('ORACLE_pool_best').mean():.3f}"
          f"   random single member {g('rand_member').mean():.3f}\n")

    print("IDENTITY CHECK FIRST, so it is not mistaken for a finding: the MINIMUM of the")
    print("alpha-tail is the pool argmin, for every alpha and every H.  Q6 says a CVaR-optimal")
    print("law lives on the tail; an order-based readout over it therefore cannot beat argmin.")
    bad = 0
    for h in HAMS:
        for a in ALPHAS:
            bad += int(np.abs(g(f"tail_min{a:g}|{h}") - g(f"argmin|{h}")).max() > 0)
    print(f"  max |tail_min - argmin| over {len(HAMS)*len(ALPHAS)} (H, alpha) cells: "
          f"{bad} cells differ  -> {'IDENTITY HOLDS' if bad == 0 else 'CHECK FAILED'}\n")

    print("THE TABLE.  Rows are READOUTS; note the basis column -- an `_avg` row is a POINT")
    print("CLOUD (not a buildable backbone, s20 X1) and every other row is a real window.\n")
    print(f"  {'readout':<16}{'basis':<12}{'control':>9}   " + "".join(f"{h:>10}" for h in HAMS))
    print(f"  {'argmin':<16}{'window':<12}{g('rand_member').mean():>9.3f}   "
          + "".join(f"{g('argmin|'+h).mean():>10.3f}" for h in HAMS))
    for a in ALPHAS:
        print(f"  -- alpha = {a:g}  (k = {int(round(a*I.K))} of {I.K})")
        for ro, ctl, basis in (("tail_member", "rand_member", "window"),
                               ("tail_medoid", f"rand_medoid{a:g}", "window"),
                               ("tail_avg", f"rand_avg{a:g}", "POINT CLOUD")):
            print(f"  {ro+str(a):<16}{basis:<12}{g(ctl).mean():>9.3f}   "
                  + "".join(f"{g(f'{ro}{a:g}|'+h).mean():>10.3f}" for h in HAMS))
        print(f"  {'ORACLE_tail_avg':<16}{'POINT CLOUD':<12}{'--':>9}   "
              f"{g(f'ORACLE_tail_avg{a:g}').mean():>10.3f}  <- ORACLE selector ceiling")

    print("\n--- PRIMARY: `tailprice.py`'s endpoint, tail_avg(0.15) - argmin, per H ---")
    prim = {}
    for h in HAMS:
        dd = g("tail_avg0.15|" + h) - g("argmin|" + h)
        mu, lo, hi = _boot(dd, rng)
        flo, fhi = _fold_boot(dd, folds, rng)
        prim[h] = mu
        print(f"  {h:<8}{mu:+7.3f} [{lo:+.3f},{hi:+.3f}] iid  [{flo:+.3f},{fhi:+.3f}] "
              f"fold-clustered   W/L {(dd<0).sum()}/{(dd>0).sum()}")
    ratio = prim["rand"] / min(prim["disto"], prim["legacy"])
    print(f"\n  ZERO-INFORMATION `rand` reproduces {100*prim['rand']/np.mean([prim['disto'],prim['legacy']]):.0f}% "
          f"of the physics Hamiltonians' mean gap.")
    print(f"  PRE-REGISTERED FALSIFIER: D2b is REFUTED if rand's gap is LESS THAN HALF the "
          f"physics gaps.")
    print(f"  rand / weakest physics = {ratio:.2f}  -> "
          f"{'D2b SUPPORTED (endpoint dominated by averaging)' if ratio > 0.5 else 'D2b REFUTED'}")

    print("\n--- THE OPERATIVE COMPARISONS: each readout against its MATCHED-COUNT control ---")
    out = {"n": n, "cfg": CFG, "primary": {}, "matched": {}}
    for h in HAMS:
        out["primary"][h] = prim[h]
    for a in ALPHAS:
        print(f"  alpha = {a:g}")
        for ro, ctl in (("tail_member", "rand_member"),
                        ("tail_medoid", f"rand_medoid{a:g}"),
                        ("tail_avg", f"rand_avg{a:g}")):
            for h in HAMS:
                dd = g(f"{ro}{a:g}|{h}") - g(ctl)
                mu, lo, hi = _boot(dd, rng)
                flo, fhi = _fold_boot(dd, folds, rng)
                sig = "" if lo < 0 < hi else ("  beats its control" if mu < 0
                                              else "  WORSE than its control")
                out["matched"][f"{ro}|{a:g}|{h}"] = {"mean": mu, "ci95": [lo, hi],
                                                     "ci95_fold": [flo, fhi],
                                                     "W": int((dd < 0).sum()),
                                                     "L": int((dd > 0).sum())}
                print(f"    {ro:<12}{h:<8}{mu:+7.3f} [{lo:+.3f},{hi:+.3f}] "
                      f"W/L {(dd<0).sum():>3}/{(dd>0).sum():<3}{sig}")

    print("\n--- THE DECOMPOSITION: what each half of `tail_avg - argmin` is worth ---")
    print("  (tail_avg - tail_medoid) = the AVERAGING operator, tail held fixed")
    print("  (tail_medoid - argmin)   = the SET vs MINIMUM operator, readout kind held fixed")
    out["decomp"] = {}
    for a in ALPHAS:
        for h in HAMS:
            d1 = g(f"tail_avg{a:g}|{h}") - g(f"tail_medoid{a:g}|{h}")
            d2 = g(f"tail_medoid{a:g}|{h}") - g(f"argmin|{h}")
            m1, l1, h1 = _boot(d1, rng); m2, l2, h2 = _boot(d2, rng)
            out["decomp"][f"{a:g}|{h}"] = {"averaging": [m1, l1, h1], "set_vs_min": [m2, l2, h2]}
            print(f"    a={a:g} {h:<8} averaging {m1:+7.3f} [{l1:+.3f},{h1:+.3f}] | "
                  f"set-vs-min {m2:+7.3f} [{l2:+.3f},{h2:+.3f}]")

    print("\nMDE at 80% power = 0.084 A.  A zero-spanning CI without power is NOT MEASURED.")
    json.dump(out, open(os.path.join(RESULTS, "d_cvarop_report.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
