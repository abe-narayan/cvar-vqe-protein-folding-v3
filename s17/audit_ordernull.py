"""s17/audit_ordernull.py -- IS THE ORACLE MAP'S CEILING CURVE A RETRIEVAL RESULT,
OR THE ORDER STATISTIC OF MORE DRAWS?

THE CLAIM UNDER ATTACK (s17/BRIEF.md section 2, s17/oracle_map.py docstring):

    "The programme has been quoting a 'K = 500 pool oracle' of 1.711 A as the retrieval
     ceiling.  That is not the retrieval ceiling -- it is the ceiling of the first 500
     entries of a BLOSUM62 ranking over a universe of 13,000-27,000 windows per target."

The sentence is true as arithmetic.  The QUESTION it is read as answering -- "there is a
better structure out there that a wider retrieval would find" -- is a different question,
and this module measures that one.

PRE-REGISTRATION (written before the run).

  H0 (the null this audit believes)   The BLOSUM62 ordering carries no usable structural
      skill.  best(K) under the BLOSUM order is then, in distribution, the same as
      best(K) over K windows drawn uniformly at random from the same universe, and the
      whole "ceiling as a function of K" curve of oracle_map.py section A is the order
      statistic min(K draws) of a fixed distribution -- i.e. widening K buys nothing
      that is not already present at K = 500 in expectation, it just takes more samples
      of the same law.

  H1 (the coordinator's implicit reading)   BLOSUM retrieval concentrates good structures
      early, and the deep universe contains genuinely DIFFERENT and better conformations
      that truncation hides.

  MEASUREMENTS
    (1) ORDER SKILL.  For each target, best-of-K under the shipped BLOSUM order vs the
        mean of best-of-K over R uniform random permutations of the SAME universe.
        Skill(K) = random_best(K) - blosum_best(K).  Positive => retrieval concentrates
        quality.  ~0 or negative => the ordering is structurally uninformative and every
        ceiling number in the map is an order statistic.
    (2) THE MATCHED NULL FOR THE HEADLINE.  The map's headline is
        best(full) - best(500).  The matched null is the SAME statistic computed under a
        random ordering: random_best(U) - random_best(500).  If the observed gain is not
        larger than the null gain, "widening K raises the ceiling" is a property of
        min(N) and not of retrieval.
    (3) THE TIE FLOOR.  K = 500 cuts inside a BLOSUM tie group.  Measure how many of the
        500 shipped pool members are inside the boundary tie group, i.e. selected by
        corpus iteration order rather than by sequence similarity.  Report the spread of
        best-of-500 over random tie-breaks -- that is the irreducible arbitrariness of
        the "1.711 A pool oracle" itself.

  SUCCESS CRITERION for H1: skill(K) must be positive with a target-level CI excluding
  zero at the K values the sprint plans to use, AND the observed full-vs-500 gain must
  exceed the matched random-order null by an interval excluding zero.

  FALSIFIER for this audit's own H0: if skill(K) is clearly positive, the BLOSUM order is
  doing real work and the ceiling map is a retrieval result, not an order statistic.  This
  module then reports that and the attack fails.

All quantities read u["rr"], which is an ORACLE label; every number here is a CEILING or a
DIAGNOSTIC and none of it is predictive.  No benchmark60 target is touched.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I     # noqa: E402
from s15 import seed as SD          # noqa: E402

UNIV = os.path.join(ROOT, "s8", "generate_univ")
KS = (25, 50, 75, 100, 150, 300, 500, 1000, 2000, 5000, 0)
R_PERM = 200      # random permutations per target


def _kname(K):
    return "full" if K == 0 else str(K)


def _load_light(pdb):
    """Only the small arrays.  W is NOT loaded -- rr is precomputed in the npz."""
    z = np.load(os.path.join(UNIV, f"{pdb}.npz"), allow_pickle=True)
    return {
        "pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
        "seq": str(z["seq"]),
        "rr": np.asarray(z["rr"], float),
        "sim": np.asarray(z["sim"], float),
        "order": np.asarray(z["order"], int),
        "org": np.asarray(z["org"], bool),
    }


def run(verbose=True):
    tg = I.targets()
    rows = []
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = _load_light(pdb)
        rr, sim, order, org = u["rr"], u["sim"], u["order"], u["org"]
        U = len(rr)
        rng = SD.stable_rng(pdb, "s17auditorder")

        # ---- (1)/(2) order skill via random permutations -------------------
        rb = rr[order]
        cum = np.minimum.accumulate(rb)              # best-of-K under BLOSUM order
        blos = {}
        for K in KS:
            k = U if K == 0 else min(K, U)
            blos[_kname(K)] = float(cum[k - 1])

        # random permutations: best-of-K for every K in one pass
        ksz = [U if K == 0 else min(K, U) for K in KS]
        acc = np.zeros((len(ksz), R_PERM))
        for r in range(R_PERM):
            p = rng.permutation(U)
            cm = np.minimum.accumulate(rr[p])
            acc[:, r] = cm[np.asarray(ksz) - 1]
        rand_mean = acc.mean(1)
        rand_lo = np.percentile(acc, 2.5, axis=1)
        rand_hi = np.percentile(acc, 97.5, axis=1)

        # ---- (3) tie floor at the shipped K = 500 --------------------------
        s = sim[order]
        boundary = s[min(500, U) - 1]
        n_above = int((sim > boundary).sum())
        n_tied = int((sim == boundary).sum())
        n_from_tie = min(500, U) - n_above           # taken out of the tie group
        tie_idx = np.where(sim == boundary)[0]
        strict = np.where(sim > boundary)[0]
        tie_best = []
        if n_from_tie > 0 and n_tied > n_from_tie:
            base = rr[strict].min() if len(strict) else np.inf
            for _ in range(R_PERM):
                pick = rng.choice(tie_idx, size=n_from_tie, replace=False)
                tie_best.append(min(base, float(rr[pick].min())))
        tie_best = np.asarray(tie_best, float)

        rows.append({
            "pdb": pdb, "n": u["n"], "fold": u["fold"], "U": int(U),
            "pep_frac": float(org.mean()),
            "pep_frac_500": float(org[order[:min(500, U)]].mean()),
            "blosum": blos,
            "rand_mean": {_kname(K): float(v) for K, v in zip(KS, rand_mean)},
            "rand_lo": {_kname(K): float(v) for K, v in zip(KS, rand_lo)},
            "rand_hi": {_kname(K): float(v) for K, v in zip(KS, rand_hi)},
            "n_unique_sim": int(len(np.unique(sim))),
            "boundary_sim": float(boundary), "n_above": n_above, "n_tied": n_tied,
            "n_from_tie": int(n_from_tie),
            "tie_best_mean": float(tie_best.mean()) if len(tie_best) else None,
            "tie_best_lo": float(np.percentile(tie_best, 2.5)) if len(tie_best) else None,
            "tie_best_hi": float(np.percentile(tie_best, 97.5)) if len(tie_best) else None,
            "tie_best_sd": float(tie_best.std(ddof=1)) if len(tie_best) else None,
        })
        if verbose and (c + 1) % 25 == 0:
            print(f"  {c+1}/{len(tg)}", flush=True)

    json.dump({"rows": rows, "R_PERM": R_PERM},
              open(os.path.join(RESULTS, "audit_ordernull.json"), "w"))
    report(rows)
    return rows


# ------------------------------------------------------------------ reporting
def _boot(d, rng, B=4000, folds=None):
    d = np.asarray(d, float)
    if folds is None:
        m = d[rng.integers(0, len(d), size=(B, len(d)))].mean(1)
    else:                                    # fold-clustered
        folds = np.asarray(folds)
        uf = np.unique(folds)
        idx = [np.where(folds == f)[0] for f in uf]
        m = np.empty(B)
        for b in range(B):
            pick = [idx[k] for k in rng.integers(0, len(uf), len(uf))]
            m[b] = d[np.concatenate(pick)].mean()
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "audit_ordernull.json")))["rows"]
    rng = SD.stable_rng("s17audit", "ordernull-report")
    ks = [_kname(K) for K in KS]
    fold = np.array([r["fold"] for r in rows])
    B = {k: np.array([r["blosum"][k] for r in rows]) for k in ks}
    Rm = {k: np.array([r["rand_mean"][k] for r in rows]) for k in ks}

    print(f"\nn = {len(rows)} targets.  R = {R_PERM} random permutations per target.\n")

    print("A. DOES THE BLOSUM62 ORDER CONCENTRATE GOOD STRUCTURES?")
    print("   skill(K) = best-of-K under a RANDOM order  -  best-of-K under BLOSUM order.")
    print("   Positive = retrieval works.  ~0 = the ceiling curve is min(K draws).")
    print(f"  {'K':>6}{'BLOSUM':>9}{'random':>9}{'skill':>9}{'   [95% CI fold-clustered]':>28}{'  W/L':>9}")
    for k in ks:
        d = Rm[k] - B[k]
        mu, lo, hi = _boot(d, rng, folds=fold)
        w = int((B[k] < Rm[k]).sum()); l = int((B[k] > Rm[k]).sum())
        print(f"  {k:>6}{B[k].mean():>9.3f}{Rm[k].mean():>9.3f}{mu:>9.3f}"
              f"   [{lo:+.3f},{hi:+.3f}]{w:>7}/{l}")

    print("\nB. THE HEADLINE, AND ITS MATCHED NULL")
    obs = B["500"] - B["full"]
    nul = Rm["500"] - Rm["full"]
    mo, lo, hi = _boot(obs, rng, folds=fold)
    mn, nlo, nhi = _boot(nul, rng, folds=fold)
    md, dlo, dhi = _boot(obs - nul, rng, folds=fold)
    print(f"  observed  best(500) - best(full)          {mo:>7.3f}  [{lo:.3f},{hi:.3f}]"
          f"   median {np.median(obs):.3f}")
    print(f"  MATCHED NULL, same statistic random order  {mn:>7.3f}  [{nlo:.3f},{nhi:.3f}]"
          f"   median {np.median(nul):.3f}")
    print(f"  observed - null                           {md:>+7.3f}  [{dlo:+.3f},{dhi:+.3f}]"
          f"   W/L {int((obs>nul).sum())}/{int((obs<nul).sum())}")
    print("  If 'observed - null' does not exclude zero, the extra ceiling at large K is")
    print("  the order statistic of more draws from the same distribution.")

    print("\nC. THE K = 500 CUT LANDS INSIDE A BLOSUM TIE GROUP")
    na = np.array([r["n_above"] for r in rows], float)
    nt = np.array([r["n_tied"] for r in rows], float)
    nf = np.array([r["n_from_tie"] for r in rows], float)
    nu = np.array([r["n_unique_sim"] for r in rows], float)
    print(f"  distinct BLOSUM scores in a universe of ~10^4 windows: mean {nu.mean():.0f}"
          f"  (min {nu.min():.0f}, max {nu.max():.0f})")
    print(f"  windows strictly above the rank-500 score:  mean {na.mean():.0f}")
    print(f"  size of the tie group AT the rank-500 score: mean {nt.mean():.0f}")
    print(f"  pool members taken FROM that tie group by corpus order alone: "
          f"mean {nf.mean():.0f}  ({(nf/500).mean():.1%} of the shipped pool)")
    tb = np.array([r["tie_best_mean"] for r in rows if r["tie_best_mean"] is not None])
    tsd = np.array([r["tie_best_sd"] for r in rows if r["tie_best_sd"] is not None])
    sel = [r for r in rows if r["tie_best_mean"] is not None]
    b500 = np.array([r["blosum"]["500"] for r in sel])
    print(f"  shipped best-of-500                                  {b500.mean():.3f}")
    print(f"  best-of-500 averaged over RANDOM tie-breaks (n={len(sel)})  {tb.mean():.3f}"
          f"   (mean per-target sd {tsd.mean():.3f})")
    print("  The 1.711 A 'pool oracle' is one draw from this tie-break distribution.")

    print("\nD. WHAT WIDENING K DOES TO THE CORPUS MIX (a confound, not a ceiling)")
    pf = np.array([r["pep_frac"] for r in rows])
    pf5 = np.array([r["pep_frac_500"] for r in rows])
    print(f"  peptide-database share of the K=500 pool  {pf5.mean():.1%}")
    print(f"  peptide-database share of the universe    {pf.mean():.1%}")
    print("  Widening K therefore also DILUTES the corpus the retrieval favours; any")
    print("  K-effect is confounded with that mix change.")

    print("\nE. UNIVERSE SIZE -- is U a structural quantity or a bookkeeping one?")
    Uv = np.array([r["U"] for r in rows], float)
    nres = np.array([r["n"] for r in rows], float)
    print(f"  corr(U, chain length n) = {np.corrcoef(Uv, nres)[0,1]:+.3f}")
    print(f"  corr(U, best(full))     = {np.corrcoef(Uv, B['full'])[0,1]:+.3f}")
    print(f"  corr(U, best(500))      = {np.corrcoef(Uv, B['500'])[0,1]:+.3f}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
