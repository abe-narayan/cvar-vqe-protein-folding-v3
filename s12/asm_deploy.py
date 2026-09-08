"""ASM-5: DEPLOYABLE two-piece assembly — no native anywhere in the decision.

Candidate generation (deployable): for every cut c, the top-N pieces of each side by
BLOSUM62 against the target's SUB-sequence, from the fold's legal library.  Every
N x N torsion concatenation is rebuilt with `build_ca_exact` and scored with the SHIPPED
leave-fold-out distogram Bayes-risk score (`I.shipped_score`) — the same objective the
production filter uses.

Arms emitted per target:
  pool_argmin   the shipped baseline: distogram argmin over the K=500 real-coordinate pool
  pool1_argmin  same pool, ideal-rebuilt from each window's own torsions (isolates rebuild)
  asm_argmin    distogram argmin over the 2-piece assembly space
  mix_argmin    distogram argmin over (pool1 U assemblies)
  asm_avg75     coordinate-average of the 75 best-scoring assemblies -> I.project
  mix_avg75     coordinate-average of the 75 best-scoring of the union -> I.project
plus ORACLE diagnostics: the best member of each candidate set, and the score percentile
of that best member.

Usage: python -m s12.asm_deploy [--N 200] [--limit L] [--noproj]
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402

LMIN = 4


def assemble_candidates(t, N=200, chunk=40000):
    """(CA (m,n,3), meta) for every deployable 2-piece assembly considered."""
    n, fold = t["n"], t["fold"]
    codes = A.encode(t["seq"])
    B62 = A.B62()
    CAs, meta = [], []
    for c in range(LMIN, n - LMIN + 1):
        b1 = A.bank(fold, c); b2 = A.bank(fold, n - c)
        s1 = B62[b1["S"], codes[:c][None, :]].sum(1)
        s2 = B62[b2["S"], codes[c:][None, :]].sum(1)
        k1 = np.argsort(-s1, kind="stable")[:N]
        k2 = np.argsort(-s2, kind="stable")[:N]
        ii, jj = np.meshgrid(np.arange(len(k1)), np.arange(len(k2)), indexing="ij")
        ii = ii.ravel(); jj = jj.ravel()
        for a in range(0, len(ii), chunk):
            p1 = k1[ii[a:a + chunk]]; p2 = k2[jj[a:a + chunk]]
            phi = np.concatenate([b1["PHI"][p1], b2["PHI"][p2]], 1)
            psi = np.concatenate([b1["PSI"][p1], b2["PSI"][p2]], 1)
            CAs.append(A.build_ca(phi, psi).astype(np.float32))
            meta.append(np.stack([np.full(len(p1), c), p1, p2], 1).astype(np.int32))
    return np.concatenate(CAs, 0), np.concatenate(meta, 0)


def main():
    av = sys.argv
    N = int(av[av.index("--N") + 1]) if "--N" in av else 200
    limit = int(av[av.index("--limit") + 1]) if "--limit" in av else None
    noproj = "--noproj" in av
    tg = sorted(I.targets(), key=lambda t: (t["fold"], t["pdb"]))
    if limit:
        tg = tg[:limit]
    res = {}; t0 = time.time(); cur = None
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); cur = t["fold"]
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]; p = I.pool_idx(u)
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(n)

        Wp = u["W"][p]
        sc_p = I.shipped_score(dg, I.pair_dists(Wp, i, j))
        rr_p = u["rr"][p]
        P1 = A.build_ca(np.asarray(u["PHI"][p], float), np.asarray(u["PSI"][p], float))
        sc_p1 = I.shipped_score(dg, I.pair_dists(P1, i, j))
        rr_p1 = I.kabsch_rmsd_batch(P1, nat)
        del u

        CA, meta = assemble_candidates(t, N=N)
        sc_a = np.empty(len(CA)); rr_a = np.empty(len(CA))
        for a in range(0, len(CA), 40000):
            blk = np.asarray(CA[a:a + 40000], float)
            sc_a[a:a + 40000] = I.shipped_score(dg, I.pair_dists(blk, i, j))
            rr_a[a:a + 40000] = I.kabsch_rmsd_batch(blk, nat)

        row = {"n": n, "fold": t["fold"], "n_asm": int(len(CA)),
               "pool_argmin": float(rr_p[int(np.argmin(sc_p))]),
               "pool1_argmin": float(rr_p1[int(np.argmin(sc_p1))]),
               "asm_argmin": float(rr_a[int(np.argmin(sc_a))]),
               "pool_best": float(rr_p.min()), "pool1_best": float(rr_p1.min()),
               "asm_best": float(rr_a.min()),
               "asm_best_pctile": float((sc_a < sc_a[int(np.argmin(rr_a))]).mean() * 100)}
        sc_m = np.concatenate([sc_p1, sc_a]); rr_m = np.concatenate([rr_p1, rr_a])
        CA_m = np.concatenate([P1.astype(np.float32), CA], 0)
        row["mix_argmin"] = float(rr_m[int(np.argmin(sc_m))])
        row["mix_best"] = float(rr_m.min())
        if not noproj:
            for tag, sc_x, CA_x in (("asm", sc_a, CA), ("mix", sc_m, CA_m)):
                top = np.argsort(sc_x)[:75]
                C, _ = I.coordinate_average(np.asarray(CA_x[top], float))
                out = I.project(C, t["seq"], t["fold"])
                row[f"{tag}_avg75"] = float(I.ca_rmsd(out["fit_ca"], nat))
                row[f"{tag}_top75_best"] = float(rr_a[top].min() if tag == "asm" else rr_m[top].min())
        res[pdb] = row
        del CA, CA_m, sc_a, rr_a, sc_m, rr_m
        print(f"{q+1:3d}/{len(tg)} {pdb} n={n} m={row['n_asm']} pool {row['pool_argmin']:.2f} "
              f"pool1 {row['pool1_argmin']:.2f} asm {row['asm_argmin']:.2f} "
              f"mix {row['mix_argmin']:.2f} asmavg {row.get('asm_avg75', float('nan')):.2f} "
              f"mixavg {row.get('mix_avg75', float('nan')):.2f} | best asm {row['asm_best']:.2f} "
              f"pct {row['asm_best_pctile']:.1f}  {time.time()-t0:.0f}s free {I.free_gb():.2f}",
              flush=True)
    I.write(f"asm_deploy_N{N}", res)
    keys = [k for k in ("pool_argmin", "pool1_argmin", "asm_argmin", "mix_argmin",
                        "asm_avg75", "mix_avg75", "pool_best", "pool1_best", "asm_best",
                        "mix_best", "asm_best_pctile") if k in res[tg[0]["pdb"]]]
    f18 = np.array([t["pdb"] in I.FAIL18 for t in tg])
    print()
    for k in keys:
        v = np.array([res[t["pdb"]][k] for t in tg])
        print(f"{k:16s} {v.mean():7.3f}   FAIL18 {v[f18].mean():7.3f}  other {v[~f18].mean():7.3f}")
    folds = np.array([t["fold"] for t in tg]); names = [t["pdb"] for t in tg]
    for k in ("asm_argmin", "mix_argmin", "asm_avg75", "mix_avg75"):
        if k not in keys:
            continue
        a = np.array([res[t["pdb"]][k] for t in tg])
        b = np.array([res[t["pdb"]]["pool_argmin"] for t in tg])
        st = I.paired(a, b, folds=folds, names=names)
        print(f"\n{k} vs shipped pool argmin: d={st['mean_diff']:+.3f} "
              f"CI{st['ci95']} W/L {st['n_better']}/{st['n_worse']} "
              f"drop10 {st['drop_top10_mean_diff']:+.3f}")
    return res


if __name__ == "__main__":
    main()
