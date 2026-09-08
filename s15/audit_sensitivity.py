"""S15 AUDIT / Part B -- is a constant LOAD-BEARING?  Measured, not asserted.

For each pinned constant of the instrument, sweep a reasonable alternative and report how
far the published number moves.  A constant whose reasonable neighbourhood moves a
published conclusion is load-bearing and must carry its provenance in the paper.

    python -m s15.audit_sensitivity
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)

KS = (100, 250, 500, 1000, 2000)
MS = (25, 50, 75, 100, 150)
BANDS = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0)
MIN_SEPS = (1, 2, 3)


def main():
    from s12 import instrument as I
    tg = I.targets()
    per = []
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        rr_all = np.asarray(u["rr"], np.float64)
        order = np.asarray(u["order"], int)
        rec = I.shipped_record(pdb)
        sub = np.asarray(rec["sub"], int)
        dg = I.distogram(pdb, t["seq"], t["fold"])
        row = {"pdb": pdb, "n": t["n"], "fold": t["fold"]}
        # --- K sweep: pool_best and the shipped selection at min_sep = 2
        for K in KS:
            p = order[:K]
            row[f"pool_best_K{K}"] = float(rr_all[p].min())
            i, j = I.pair_index(t["n"], 2)
            D = I.pair_dists(u["W"][p], i, j).astype(np.float32).astype(float)
            sc = I.shipped_score(dg, D)
            row[f"shipped_K{K}"] = float(rr_all[p][int(np.argmin(sc))])
        # --- m sweep at K = 500, using the shipped score's own order
        p = order[:500]
        i, j = I.pair_index(t["n"], 2)
        D = I.pair_dists(u["W"][p], i, j).astype(np.float32).astype(float)
        sc = I.shipped_score(dg, D)
        o = np.argsort(sc, kind="stable")
        rrp = rr_all[p]
        for m in MS:
            row[f"topm_best_m{m}"] = float(rrp[o[:m]].min())
        row["sub_matches_score_top75"] = bool(np.array_equal(np.sort(o[:75]), np.sort(sub)))
        # --- BAND sweep: zero-recall count uses the PRODUCTION sub
        for b in BANDS:
            band = np.where(rrp <= rrp.min() + b)[0]
            row[f"zero_recall_band{b}"] = bool(not np.isin(band, sub).any())
        # --- min_sep sweep for the shipped score
        for ms in MIN_SEPS:
            ii, jj = I.pair_index(t["n"], ms)
            # the cached distogram stores its own i,j; only min_sep=2 is directly scorable
            if len(ii) == len(dg["i"]) and np.array_equal(ii, dg["i"]) and np.array_equal(jj, dg["j"]):
                Dm = I.pair_dists(u["W"][p], ii, jj).astype(np.float32).astype(float)
                row[f"shipped_minsep{ms}"] = float(rrp[int(np.argmin(I.shipped_score(dg, Dm)))])
            else:
                row[f"shipped_minsep{ms}"] = None
        per.append(row)
        if len(per) % 25 == 0:
            print(f"[{len(per)}/126]", flush=True)
            with open(os.path.join(OUT, "audit_sensitivity.json"), "w") as fh:
                json.dump({"per_target": per, "n_expected": 126,
                           "complete": False}, fh, indent=1)

    def mean(k):
        v = [r[k] for r in per if r.get(k) is not None]
        return float(np.mean(v)) if v else None

    agg = {
        "n": len(per),
        "pool_best_vs_K": {K: mean(f"pool_best_K{K}") for K in KS},
        "shipped_vs_K": {K: mean(f"shipped_K{K}") for K in KS},
        "topm_best_vs_m": {m: mean(f"topm_best_m{m}") for m in MS},
        "n_zero_recall_vs_BAND": {b: int(sum(r[f"zero_recall_band{b}"] for r in per)) for b in BANDS},
        "shipped_vs_min_sep": {ms: mean(f"shipped_minsep{ms}") for ms in MIN_SEPS},
        "n_targets_where_score_top75_equals_production_sub":
            int(sum(r["sub_matches_score_top75"] for r in per)),
    }
    print(json.dumps(agg, indent=1))
    with open(os.path.join(OUT, "audit_sensitivity.json"), "w") as fh:
        json.dump({"aggregate": agg, "per_target": per, "n_expected": 126,
                   "complete": len(per) == 126}, fh, indent=1)
    print("wrote", os.path.join(OUT, "audit_sensitivity.json"))


if __name__ == "__main__":
    main()
