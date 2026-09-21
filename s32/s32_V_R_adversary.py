"""S32 LANE V -- the adversarial audit of LANE R's projection-branch arms.

Run INDEPENDENTLY of lane R, from lane R's raw rows (`s32/results/s32_R_branches*.jsonl`),
which carry every branch's ORACLE `rmsd_nat` and every native-free criterion side by side.

THE ATTACK, in the order the brief names it:

 1. **ORDER STATISTIC.**  `min(rmsd_nat)` over B ~ 300 branches is a best-of-B and will sit
    far below 3.2105 for arithmetic reasons alone.  It is priced here against
    `s24.stats_lib.best_of_k_within` (deviations resampled ACROSS targets) and, more
    importantly, against `split_half_transfer`, which nulls itself.  A per-target minimum over
    K variants is mostly an order statistic.

 2. **THE ONLY QUESTION THAT MATTERS IS THE DEPLOYABLE ONE.**  For every native-free criterion
    lane R records, the deployable arm is `argmin(criterion)` per target, scored on the built
    chain, paired against the SAME JOB's `prod_chain`.  Reported with SE, MDE, x MDE, fold
    structure, W/L and median -- contract rule 1.

 3. **ZERO-INFORMATION CONTROL MATCHED TO THE OPERATOR'S SPACE** (contract rule 6, and the
    project's most repeated error).  Not "uniform on the torus": a uniformly random branch
    from the SAME branch set, with its OWN draw distribution (>= 200 draws, mean and
    draw-to-draw sd -- contract rule 10).  Plus `d_to_C`, "take the branch closest to the
    cloud", which is the natural zero-information rule in this operator's own space.

 4. **COMPUTE MISMATCH.**  Production searches 4 starts.  An arm that searches 300 and wins is
    not a fairer selector, it is a bigger search.  Every criterion is therefore also scored on
    the GEN4-only sub-branch set, which is compute-matched to production.

 5. **TIE HANDLING.**  `argmin` on a tied criterion reads the storage order and can invent a
    winner (project memory: a tied argmin invented a 1.386 A winner).  Every argmin here is
    averaged over the tied argmin set via `s24.stats_lib.argmin_tied`.

 6. **BRANCH DISTINCTNESS.**  If the "branches" are one basin reached with different round-off
    they are arithmetic noise, not structural alternatives.  Reported: the pairwise-RMSD
    cluster count at lane R's own 1e-3 threshold, and the median spread of `rmsd_nat` within a
    cluster (which must be ~0 if the clustering means what it says).

Everything that reads `rmsd_nat` for SELECTION is labelled ORACLE / NOT DEPLOYABLE.
"""
from __future__ import annotations
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s24 import stats_lib as ST                                          # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
ROWS = os.path.join(RESULTS, "s32_R_branches*.jsonl")
N_DRAWS = 300

#: criteria lane R records per branch.  Sign convention: all are minimised unless listed in
#: MAXIMISE.  `rmsd_nat` and `d_to_prod` are ORACLE / diagnostic and never a deployable arm.
ORACLE_KEYS = {"rmsd_nat", "d_to_prod"}
MAXIMISE = {"typicality"}


def load_rows():
    out = {}
    for f in sorted(glob.glob(ROWS)):
        for line in open(f):
            line = line.strip()
            if line:
                r = json.loads(line)
                out[r["pdb"]] = r
    return out


def criteria_keys(r):
    B = int(r["n_branches"])
    return sorted(k for k, v in r.items()
                  if isinstance(v, list) and len(v) == B and k not in ORACLE_KEYS
                  and k not in ("fam", "tag", "cluster")
                  and all(isinstance(x, (int, float)) for x in v[:5]))


def pick(score, rmsd, maximise=False):
    """The outcome of argmin(score), AVERAGED over the tied argmin set (project memory:
    np.argmin on a tied signal reads the storage order and invents a winner)."""
    s = -np.asarray(score, float) if maximise else np.asarray(score, float)
    if not np.isfinite(s).any():
        return float("nan"), 0
    s = np.where(np.isfinite(s), s, np.inf)
    m = s.min()
    tied = np.flatnonzero(np.isclose(s, m, rtol=1e-9, atol=1e-12))
    return float(np.asarray(rmsd, float)[tied].mean()), len(tied)


def main():
    rows = load_rows()
    pdbs = sorted(rows)
    if len(pdbs) < 20:
        print("only %d targets present -- lane R is still running; rerun when complete" % len(pdbs))
    folds = np.array([rows[p]["fold"] for p in pdbs], int)
    prod = np.array([rows[p]["prod_chain"] for p in pdbs], float)
    out = dict(n=len(pdbs), prod_mean=float(prod.mean()), basis="chain (built chain, same job)",
               note="production chain recomputed IN LANE R'S OWN JOB; every arm is paired to it")

    print("=" * 104)
    print("S32-V adversarial audit of lane R.  n=%d targets, production chain (same job) %.4f"
          % (len(pdbs), prod.mean()))

    # ---------------------------------------------------------------- 6. branch distinctness
    B = np.array([rows[p]["n_branches"] for p in pdbs], float)
    D = np.array([rows[p]["n_distinct"] for p in pdbs], float)
    within = []
    for p in pdbs:
        lab = np.asarray(rows[p]["cluster"], int); rn = np.asarray(rows[p]["rmsd_nat"], float)
        for c in np.unique(lab):
            m = lab == c
            if m.sum() > 1:
                within.append(rn[m].max() - rn[m].min())
    out["branches"] = dict(mean_n_branches=float(B.mean()), mean_n_distinct=float(D.mean()),
                           frac_distinct=float((D / B).mean()),
                           within_cluster_rmsd_spread_median=float(np.median(within)) if within else None,
                           within_cluster_rmsd_spread_p95=float(np.percentile(within, 95)) if within else None)
    print("  branches/target %.1f, distinct at 1e-3 A %.1f (%.0f%%); within-cluster rmsd_nat "
          "spread median %.2e p95 %.2e"
          % (B.mean(), D.mean(), 100 * (D / B).mean(),
             out["branches"]["within_cluster_rmsd_spread_median"] or 0,
             out["branches"]["within_cluster_rmsd_spread_p95"] or 0))

    # ---------------------------------------------------------- 1. the ORACLE order statistic
    M = None
    kk = min(int(b) for b in B)
    M = np.array([np.asarray(rows[p]["rmsd_nat"], float)[:kk] for p in pdbs])   # (n, k) grid
    orc = np.array([np.asarray(rows[p]["rmsd_nat"], float).min() for p in pdbs])
    bok = ST.best_of_k_within(M)
    out["oracle_best_branch"] = dict(
        mean=float(orc.mean()), gain_vs_prod=float(orc.mean() - prod.mean()),
        oracle="ORACLE / NOT DEPLOYABLE",
        best_of_k_within=bok, k_used=int(kk))
    print()
    print("  1. ORACLE best branch                 %.4f  (%+.4f vs production)  ORACLE / NOT DEPLOYABLE"
          % (orc.mean(), orc.mean() - prod.mean()))
    print("     best-of-K null (deviations resampled ACROSS targets): observed %+.4f, null %+.4f, "
          "%.0f%% accounted" % (bok["observed_gain"], bok["null_across_targets"],
                                100 * bok["share_accounted"]))
    print("     SPLIT-HALF TRANSFER (the number to quote): %+.4f  = %.1f%% of the oracle  -> %s"
          % (bok["split_half"], 100 * bok["split_half_frac"], bok["verdict"]))

    # --------------------------------- 3. the zero-information control, with its own distribution
    rng = np.random.default_rng(32_0032)
    draws = np.empty((N_DRAWS, len(pdbs)))
    for d in range(N_DRAWS):
        for i, p in enumerate(pdbs):
            rn = np.asarray(rows[p]["rmsd_nat"], float)
            draws[d, i] = rn[rng.integers(0, len(rn))]
    dm = draws.mean(1)
    out["random_branch_control"] = dict(draw_mean=float(dm.mean()), draw_sd=float(dm.std(ddof=1)),
                                        n_draws=N_DRAWS, vs_prod=float(dm.mean() - prod.mean()),
                                        matched_to="the SAME branch set the criteria choose from")
    print()
    print("  3. ZERO-INFORMATION CONTROL (uniform random branch, %d draws): %.4f +- %.4f "
          "(draw-to-draw sd)  %+.4f vs production"
          % (N_DRAWS, dm.mean(), dm.std(ddof=1), dm.mean() - prod.mean()))
    print("     matched to the operator's own space: the same branch set every criterion picks from")

    # ------------------------------------------ 2/4/5. every native-free criterion, deployable
    keys = criteria_keys(rows[pdbs[0]])
    print()
    print("  2. DEPLOYABLE ARMS -- argmin of each native-free criterion, chain basis, paired to "
          "the same job's production")
    print("     %-16s %8s %9s %8s %7s %7s %6s %8s  %s"
          % ("criterion", "mean", "effect", "SE", "MDE", "xMDE", "folds", "W/L", "verdict"))
    out["arms"] = {}
    for sub, lab in (("ALL", "all branches"), ("GEN4", "GEN4 only (compute-matched to production)")):
        for k in keys:
            v, nt = [], []
            for p in pdbs:
                rn = np.asarray(rows[p]["rmsd_nat"], float)
                sc = np.asarray(rows[p][k], float)
                if sub == "GEN4":
                    m = np.array([f == "GEN4" for f in rows[p]["fam"]])
                    if m.sum() == 0:
                        v.append(np.nan); nt.append(0); continue
                    rn, sc = rn[m], sc[m]
                a, t = pick(sc, rn, maximise=(k in MAXIMISE))
                v.append(a); nt.append(t)
            v = np.array(v, float)
            if not np.isfinite(v).all():
                continue
            c = ST.compare(v, prod, folds, names=pdbs,
                           label="branch argmin(%s), %s vs production (BUILT CHAIN)" % (k, lab))
            out["arms"]["%s|%s" % (sub, k)] = dict(compare=c, mean_tied=float(np.mean(nt)), subset=lab)
            gate = ("RESULT" if abs(c["effect_over_mde"]) >= 1.0 and c["folds_same_sign"] >= 4
                    and max(c["ci95_fold"]) * min(c["ci95_fold"]) > 0
                    else "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "not a result")
            print("     %-16s %8.4f %+9.4f %8.4f %7.4f %+7.2f %4d/5 %4d/%-3d  %s%s"
                  % (k[:16], v.mean(), c["effect"], c["se"], c["mde"], c["effect_over_mde"],
                     c["folds_same_sign"], c["n_better"], c["n_worse"], gate,
                     "" if np.mean(nt) < 1.05 else "  (mean %.1f tied)" % np.mean(nt)))
        print("     --- %s above ---" % lab)

    # split-half transfer for the best deployable criterion, treating criteria as the K grid
    names = [k for k in out["arms"] if k.startswith("ALL|")]
    if names:
        G = np.array([[out["arms"][k]["compare"]["mean_a"] for k in names]])  # placeholder shape
        best = min(names, key=lambda k: out["arms"][k]["compare"]["effect"])
        out["best_deployable_arm"] = best
        print()
        print("  4/5. best deployable arm over %d criteria x 2 subsets = %d comparisons: %s"
              % (len(keys), 2 * len(keys), best))
        print("       THE SEARCH ITSELF IS A MULTIPLICITY COST -- the best of %d comparisons needs "
              "the search accounted (charter 45).  Its nominal MDE is not its achieved MDE."
              % (2 * len(keys)))
        out["n_comparisons_emitted"] = 2 * len(keys)

    with open(os.path.join(RESULTS, "s32_V_R_adversary.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("=" * 104)
    return out


if __name__ == "__main__":
    main()
