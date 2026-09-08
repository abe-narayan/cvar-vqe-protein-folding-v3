"""ADVERSARIAL AUDIT 1c(ii) -- the measurement `s12/asm_deploy.py` did not make.

asm reports that deployable assembly improves the candidate set's ORACLE BEST from
1.711 -> 1.175 A while the emitted answer gets +0.218 A WORSE, and reads that as
"the objective cannot rank within the improved set".

`s12/adv_operator.py` shows the terminal operator's output is a function of the set MEAN
(pooled OLS slope 1.16, R^2 0.89) with a coefficient of 0.039 on the set BEST.  So the
asm result is fully explained without any claim about ranking -- UNLESS the assembly bank's
best-possible 75-member SET also beats the pool's.  That is the missing number:

    ORACLE FILTER on the assembly bank: the 75 lowest-rr assemblies, coordinate-averaged
    and projected -- versus the same thing on the K=500 pool (1.963 avg / 2.102 projected).

If the assembly bank's perfect-filter emission BEATS the pool's, assembly is a real
generation gain that a better objective could cash, and asm's reading survives.
If it does NOT, then assembly never had anything the operator could use and the
"objective cannot rank" reading is unsupported by this experiment.

ORACLE/DIAGNOSTIC: `rr` selects the members.  Not deployable; an attribution measurement.
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A
from s12.asm_deploy import assemble_candidates


def main(limit=24, N=200):
    tg = sorted(I.targets(), key=lambda t: (t["fold"], t["pdb"]))[:limit]
    res = {}; t0 = time.time(); cur = None
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); cur = t["fold"]
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]; p = I.pool_idx(u)
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(n)
        Wp = np.asarray(u["W"][p], float); rr_p = u["rr"][p]
        sc_p = I.shipped_score(dg, I.pair_dists(Wp, i, j))
        del u

        row = dict(n=n, fold=t["fold"], fail=int(pdb in I.FAIL18))
        # ---------- POOL reference arms
        sub = np.argsort(sc_p, kind="stable")[:75]
        orc = np.argsort(rr_p, kind="stable")[:75]
        for tag, idx in (("pool_score75", sub), ("pool_oracle75", orc)):
            C, _ = I.coordinate_average(Wp[idx])
            row[f"{tag}_avg"] = float(I.ca_rmsd(C, nat))
            row[f"{tag}_fit"] = float(I.ca_rmsd(I.project(C, t["seq"], t["fold"])["fit_ca"], nat))
            row[f"{tag}_setmean"] = float(rr_p[idx].mean())
            row[f"{tag}_setbest"] = float(rr_p[idx].min())
        row["pool_best"] = float(rr_p.min())

        # ---------- ASSEMBLY bank
        CA, meta = assemble_candidates(t, N=N)
        sc_a = np.empty(len(CA)); rr_a = np.empty(len(CA))
        for a in range(0, len(CA), 40000):
            blk = np.asarray(CA[a:a + 40000], float)
            sc_a[a:a + 40000] = I.shipped_score(dg, I.pair_dists(blk, i, j))
            rr_a[a:a + 40000] = I.kabsch_rmsd_batch(blk, nat)
        row["n_asm"] = int(len(CA))
        row["asm_best"] = float(rr_a.min())
        suba = np.argsort(sc_a, kind="stable")[:75]
        orca = np.argsort(rr_a, kind="stable")[:75]
        for tag, idx in (("asm_score75", suba), ("asm_oracle75", orca)):
            C, _ = I.coordinate_average(np.asarray(CA[idx], float))
            row[f"{tag}_avg"] = float(I.ca_rmsd(C, nat))
            row[f"{tag}_fit"] = float(I.ca_rmsd(I.project(C, t["seq"], t["fold"])["fit_ca"], nat))
            row[f"{tag}_setmean"] = float(rr_a[idx].mean())
            row[f"{tag}_setbest"] = float(rr_a[idx].min())
        # oracle-25 on both, since the operator likes small good sets (record C1)
        for tag, rr_x, CA_x in (("pool", rr_p, Wp), ("asm", rr_a, CA)):
            idx = np.argsort(rr_x, kind="stable")[:25]
            C, _ = I.coordinate_average(np.asarray(CA_x[idx], float))
            row[f"{tag}_oracle25_avg"] = float(I.ca_rmsd(C, nat))
        res[pdb] = row
        del CA, sc_a, rr_a
        print(f"{q+1:3d}/{len(tg)} {pdb} n={n} m={row['n_asm']} | "
              f"pool o75 {row['pool_oracle75_avg']:.3f} asm o75 {row['asm_oracle75_avg']:.3f} | "
              f"pool s75 {row['pool_score75_avg']:.3f} asm s75 {row['asm_score75_avg']:.3f} | "
              f"setmean {row['pool_oracle75_setmean']:.2f}/{row['asm_oracle75_setmean']:.2f} "
              f"{time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)

    names = list(res)
    fail = np.array([res[k]["fail"] for k in names], bool)
    agg = {"n": len(names), "n_fail": int(fail.sum())}
    for k in ("pool_score75_avg", "pool_score75_fit", "pool_oracle75_avg", "pool_oracle75_fit",
              "asm_score75_avg", "asm_score75_fit", "asm_oracle75_avg", "asm_oracle75_fit",
              "pool_score75_setmean", "pool_oracle75_setmean", "asm_score75_setmean",
              "asm_oracle75_setmean", "pool_score75_setbest", "pool_oracle75_setbest",
              "asm_score75_setbest", "asm_oracle75_setbest", "pool_best", "asm_best",
              "pool_oracle25_avg", "asm_oracle25_avg"):
        v = np.array([res[x][k] for x in names])
        agg[k] = dict(mean=float(v.mean()), fail=float(v[fail].mean()) if fail.any() else None,
                      other=float(v[~fail].mean()))
    for a, b in (("asm_oracle75_avg", "pool_oracle75_avg"),
                 ("asm_oracle75_fit", "pool_oracle75_fit"),
                 ("asm_oracle25_avg", "pool_oracle25_avg"),
                 ("asm_score75_fit", "pool_score75_fit")):
        agg[f"paired_{a}_vs_{b}"] = I.paired(np.array([res[x][a] for x in names]),
                                             np.array([res[x][b] for x in names]),
                                             names=names)
    print(json.dumps({k: v for k, v in agg.items() if not k.startswith("paired")}, indent=1))
    for k in agg:
        if k.startswith("paired"):
            v = agg[k]
            print(f"{k}: d={v['mean_diff']:+.3f} CI{v['ci95']} W/L {v['n_better']}/{v['n_worse']}")
    I.write("adv_asmoracle", dict(agg=agg, rows=res))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 24)
