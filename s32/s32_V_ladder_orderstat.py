"""S32 LANE V -- adversarial audit of the S29 O-ladder's NARROWING increments.

The coordinator's ladder writes

    best single member, K=500                     1.7078
      + retrieval filter K=500 -> 128            +0.4357  ->  2.1435
      + prefix 128 -> 75                         +0.1620  ->  2.3055

Each rung is a PER-TARGET MINIMUM over a set, and each increment is the difference of two
such minima over NESTED sets (s29_O_p128.json: overlap_top75_in_128 == 1.0 for all 126).
**A per-target minimum over K variants is mostly an order statistic** (project memory,
"Grid oracles are order statistics"; S32 contract rule 9).  Therefore an increment labelled
"the cost of the retrieval filter" is only that if a SIZE-MATCHED RANDOM subset of the same
parent set costs the same or more.  This script prices that null.

Decomposition, on the member (= deposited-structure) basis, where the projection price is
-0.0007 .. -0.0030 A, i.e. ~free (S32 contract rule 16):

    size effect  = E[min over a uniformly random 128 of the 500] - min over all 500
    order effect = min over the DIS top-128                      - E[random 128]

and identically for 128 -> 75.  A positive order effect means the DIS ordering is WORSE than
chance at retaining the oracle-best member; ~0 means the ordering contributes nothing and the
whole increment is the order statistic; negative means the ordering genuinely helps.

Everything here is ORACLE / NOT DEPLOYABLE (it reads `rr`, the native-CA RMSD labels).
"""
from __future__ import annotations
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)
NDRAW = 2000          # draws per target; the draw-to-draw sd is reported (contract rule 10)


def main():
    tg = I.targets()
    assert len(tg) == 126
    rng = np.random.default_rng(20260921)
    rows = {}
    for t in tg:
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb)
        pool = I.pool_idx(u, I.K)
        rr = u["rr"][pool]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        D = I.pair_dists(u["W"][pool], i, j)
        sc = I.shipped_score(dg, D.astype(np.float32).astype(float))
        order = np.argsort(sc, kind="stable")
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        assert set(order[:75].tolist()) == set(sub.tolist()), "top-75 does not reproduce for %s" % pdb

        t128, t75 = order[:128], order[:75]
        assert set(t75.tolist()) <= set(t128.tolist())

        # --- size-matched random nulls, each with its OWN draw distribution
        idx = rng.random((NDRAW, 500)).argsort(1)[:, :128]          # random 128 of 500
        r500_128 = rr[idx].min(1)
        sub128 = t128[rng.random((NDRAW, 128)).argsort(1)[:, :75]]   # random 75 of the DIS top-128
        r128_75 = rr[sub128].min(1)

        rows[pdb] = dict(
            fold=int(fold),
            best500=float(rr.min()), best128=float(rr[t128].min()), best75=float(rr[t75].min()),
            rand128_mean=float(r500_128.mean()), rand128_sd=float(r500_128.std(ddof=1)),
            rand75of128_mean=float(r128_75.mean()), rand75of128_sd=float(r128_75.std(ddof=1)),
            rank_of_best500=int(np.argsort(np.argsort(sc, kind="stable"), kind="stable")[int(np.argmin(rr))]),
            best500_in_top128=bool(int(np.argmin(rr)) in set(t128.tolist())),
            best500_in_top75=bool(int(np.argmin(rr)) in set(t75.tolist())),
        )
        print("  %-6s best500 %.3f | top128 %.3f rand128 %.3f | top75 %.3f rand75 %.3f | rank(best) %4d"
              % (pdb, rows[pdb]["best500"], rows[pdb]["best128"], rows[pdb]["rand128_mean"],
                 rows[pdb]["best75"], rows[pdb]["rand75of128_mean"], rows[pdb]["rank_of_best500"]), flush=True)

    pdbs = sorted(rows)
    G = lambda k: np.array([rows[p][k] for p in pdbs], float)
    b500, b128, b75 = G("best500"), G("best128"), G("best75")
    r128, r75 = G("rand128_mean"), G("rand75of128_mean")

    out = dict(n=len(pdbs), n_draws=NDRAW, basis="member (deposited window); projection price "
               "-0.0007..-0.0030 A for this object class (S32 contract rule 16)",
               oracle=True, deployable=False, per_target=rows)
    out["means"] = dict(best500=float(b500.mean()), best128=float(b128.mean()), best75=float(b75.mean()),
                        rand128=float(r128.mean()), rand75of128=float(r75.mean()))
    out["increments"] = {
        "K=500 -> 128 TOTAL": float(b128.mean() - b500.mean()),
        "  size effect (random 128)": float(r128.mean() - b500.mean()),
        "  ordering effect (DIS 128 - random 128)": float(b128.mean() - r128.mean()),
        "128 -> 75 TOTAL": float(b75.mean() - b128.mean()),
        "  size effect (random 75 of 128)": float(r75.mean() - b128.mean()),
        "  ordering effect (DIS 75 - random 75)": float(b75.mean() - r75.mean()),
    }
    out["draw_sd_mean"] = dict(rand128=float(G("rand128_sd").mean()), rand75of128=float(G("rand75of128_sd").mean()))
    out["recall"] = dict(best500_in_top128=float(G("best500_in_top128").mean()),
                         best500_in_top75=float(G("best500_in_top75").mean()),
                         median_rank_of_best500=float(np.median(G("rank_of_best500"))))
    print()
    print("=" * 92)
    print("ORACLE / NOT DEPLOYABLE.  member basis, n=%d, %d draws/target" % (len(pdbs), NDRAW))
    for k, v in out["means"].items():
        print("  %-14s %.4f" % (k, v))
    print()
    for k, v in out["increments"].items():
        print("  %-44s %+.4f" % (k, v))
    print()
    print("  best-of-500 member survives into the DIS top-128 on %.1f%% of targets, top-75 %.1f%%; "
          "median DIS rank of the oracle-best member = %.0f of 500"
          % (100 * out["recall"]["best500_in_top128"], 100 * out["recall"]["best500_in_top75"],
             out["recall"]["median_rank_of_best500"]))
    print("=" * 92)
    with open(os.path.join(RESULTS, "s32_V_ladder_orderstat.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    main()
