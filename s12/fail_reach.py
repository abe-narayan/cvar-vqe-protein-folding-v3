"""FAIL18 forensics, step 12 (coordinator request): are the lasso and fibril classes
reachable IN PRINCIPLE by a linear-window retrieval architecture?

For every target: how many windows of the whole leakage-safe universe come within 2.0 /
2.5 / 3.0 A of the native CA trace, where they sit in the BLOSUM ranking (i.e. whether a
K=500 sequence-similarity query could ever have found them), and what the best single
window can do.  Reported per structural class.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I

RES = os.path.join(ROOT, "s12", "results")


def main():
    hdr = json.load(open(os.path.join(RES, "fail_headers.json")))["per_target"]
    doss = {d["pdb"]: d for d in json.load(open(os.path.join(RES, "fail_dossier.json")))}
    rows = []
    for t in I.targets():
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb)
        rr = u["rr"]; order = np.asarray(u["order"], int)
        rank = np.empty(len(rr), int); rank[order] = np.arange(len(rr))
        org = np.asarray(u["org"], bool)
        h = hdr[pdb]
        out = dict(pdb=pdb, n=n, fold=t["fold"], fail18=pdb in I.FAIL18,
                   lasso=bool(h["lasso"]), fibril=bool(h["fibril"]),
                   nw=int(len(rr)), uni_best=float(rr.min()),
                   pool_best=doss[pdb]["o_pool_best"],
                   emitted=doss[pdb]["rmsd_emitted"])
        for th in (2.0, 2.5, 3.0):
            m = rr <= th
            out[f"n_le_{th}"] = int(m.sum())
            out[f"frac_le_{th}"] = float(m.mean())
            out[f"medrank_le_{th}"] = float(np.median(rank[m])) if m.any() else None
            out[f"in500_le_{th}"] = int((rank[m] < I.K).sum()) if m.any() else 0
            out[f"org_le_{th}"] = float(org[m].mean()) if m.any() else None
        # what an ORACLE query (best-K by true rmsd rank vs BLOSUM rank) would need
        best_rank = int(rank[int(np.argmin(rr))])
        out["blosum_rank_of_universe_best"] = best_rank
        rows.append(out)
    I.write("fail_reach", rows)

    def grp(sel, name):
        g = [r for r in rows if sel(r)]
        if not g:
            return
        print(f"{name:26s} n={len(g):3d}  uni_best={np.mean([r['uni_best'] for r in g]):5.2f} "
              f"pool_best={np.mean([r['pool_best'] for r in g]):5.2f} "
              f"emitted={np.mean([r['emitted'] for r in g]):5.2f} | "
              f"windows<=2.0A: {np.mean([r['n_le_2.0'] for r in g]):8.1f} "
              f"({np.mean([r['frac_le_2.0'] for r in g])*100:5.2f}% of universe), "
              f"of which in K=500: {np.mean([r['in500_le_2.0'] for r in g]):6.1f} | "
              f"BLOSUM rank of universe best = {np.median([r['blosum_rank_of_universe_best'] for r in g]):8.0f}")

    print("== reachability by class ==")
    grp(lambda r: True, "ALL 126")
    grp(lambda r: not r["fail18"], "other-108")
    grp(lambda r: r["fail18"], "FAIL18")
    grp(lambda r: r["lasso"], "lasso (all 6)")
    grp(lambda r: r["fibril"] and r["fail18"], "fibril in FAIL18 (6)")
    grp(lambda r: r["fibril"] and not r["fail18"], "fibril not in FAIL18 (4)")
    grp(lambda r: r["fail18"] and not (r["lasso"] or r["fibril"]), "FAIL18 other (8)")

    print("\n== per-target, lasso + fibril ==")
    print(f"{'pdb':6s}{'n':>3s}{'cls':>8s}{'uni':>6s}{'pool':>6s}{'emit':>6s}{'<=2A':>7s}"
          f"{'in500':>7s}{'%uni':>7s}{'bestrank':>9s}")
    for r in rows:
        if not (r["lasso"] or r["fibril"]):
            continue
        cls = "lasso" if r["lasso"] else "fibril"
        print(f"{r['pdb']:6s}{r['n']:3d}{cls:>8s}{r['uni_best']:6.2f}{r['pool_best']:6.2f}"
              f"{r['emitted']:6.2f}{r['n_le_2.0']:7d}{r['in500_le_2.0']:7d}"
              f"{r['frac_le_2.0']*100:7.2f}{r['blosum_rank_of_universe_best']:9d}")


if __name__ == "__main__":
    main()
