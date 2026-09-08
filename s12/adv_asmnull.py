"""ADVERSARIAL AUDIT 1c(iii) -- the CAPACITY NULL for the assembly oracle-filter result.

`s12/adv_asmoracle.py` finds that the assembly bank's ORACLE top-75 emits 1.054 A against
the K=500 pool's 2.073 A.  Before that can mean anything, it must survive the control the
record already used to kill the assembly ladder (`asm_FINDINGS.md` s2b, and record C4):
choosing the best 75 of 2.3e5 candidates is 460x more selection freedom than choosing the
best 75 of 500, and an ORACLE selection converts raw cardinality into apparent accuracy.

The null: m i.i.d. chains drawn from the LIBRARY'S OWN Ramachandran marginal -- same count,
same builder, same ideal geometry, ZERO sequence information and zero fragment structure --
put through the identical oracle top-75 -> coordinate average -> project path.

If the null reaches ~1.0 A too, the assembly "generation gain" is a capacity artefact and
the sprint must not treat it as headroom.
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A


def rama_marginal(fold, cmax=8, cap=40000, seed=0):
    """(phi, psi) pairs pooled from the fold's own piece banks -- the library marginal."""
    rng = np.random.default_rng(seed)
    P = []
    for c in range(4, cmax + 1):
        try:
            b = A.bank(fold, c)
        except Exception:                                              # noqa: BLE001
            continue
        ph = np.asarray(b["PHI"], float).ravel(); ps = np.asarray(b["PSI"], float).ravel()
        k = rng.choice(len(ph), min(cap, len(ph)), replace=False)
        P.append(np.stack([ph[k], ps[k]], 1))
    return np.concatenate(P, 0)


def main(limit=24, seed=0):
    with open(os.path.join(I.RESULTS, "adv_asmoracle.json")) as fh:
        prev = json.load(fh)["rows"]
    tg = [t for t in sorted(I.targets(), key=lambda t: (t["fold"], t["pdb"])) if t["pdb"] in prev][:limit]
    rng = np.random.default_rng(seed)
    res = {}; t0 = time.time(); cur = None; MARG = None
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); cur = t["fold"]; MARG = rama_marginal(cur, seed=seed)
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]; del u
        m = int(prev[pdb]["n_asm"])
        rr = np.empty(m); CAall = np.empty((m, n, 3), np.float32)
        for a in range(0, m, 40000):
            b = min(40000, m - a)
            k = rng.integers(0, len(MARG), (b, n))
            ca = I.build_ca(MARG[k, 0], MARG[k, 1])
            CAall[a:a + b] = ca.astype(np.float32)
            rr[a:a + b] = I.kabsch_rmsd_batch(ca, nat)
        row = dict(n=n, fold=t["fold"], m=m, null_best=float(rr.min()),
                   asm_best=float(prev[pdb]["asm_best"]),
                   pool_best=float(prev[pdb]["pool_best"]))
        for kk, tag in ((75, "o75"), (25, "o25")):
            idx = np.argsort(rr, kind="stable")[:kk]
            C, _ = I.coordinate_average(np.asarray(CAall[idx], float))
            row[f"null_{tag}_avg"] = float(I.ca_rmsd(C, nat))
            row[f"null_{tag}_setmean"] = float(rr[idx].mean())
            if kk == 75:
                row["null_o75_fit"] = float(
                    I.ca_rmsd(I.project(C, t["seq"], t["fold"])["fit_ca"], nat))
        row["asm_o75_avg"] = float(prev[pdb]["asm_oracle75_avg"])
        row["asm_o75_fit"] = float(prev[pdb]["asm_oracle75_fit"])
        row["asm_o75_setmean"] = float(prev[pdb]["asm_oracle75_setmean"])
        row["pool_o75_avg"] = float(prev[pdb]["pool_oracle75_avg"])
        row["pool_o75_fit"] = float(prev[pdb]["pool_oracle75_fit"])
        row["asm_o25_avg"] = float(prev[pdb]["asm_oracle25_avg"])
        row["pool_o25_avg"] = float(prev[pdb]["pool_oracle25_avg"])
        res[pdb] = row
        del CAall, rr
        print(f"{q+1:3d}/{len(tg)} {pdb} m={m} null_o75 {row['null_o75_avg']:.3f} "
              f"asm_o75 {row['asm_o75_avg']:.3f} pool_o75 {row['pool_o75_avg']:.3f} "
              f"| best null {row['null_best']:.3f} asm {row['asm_best']:.3f} "
              f"{time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)

    names = list(res)
    agg = {"n": len(names)}
    for k in ("null_best", "asm_best", "pool_best", "null_o75_avg", "asm_o75_avg",
              "pool_o75_avg", "null_o75_fit", "asm_o75_fit", "pool_o75_fit",
              "null_o25_avg", "asm_o25_avg", "pool_o25_avg", "null_o75_setmean",
              "asm_o75_setmean"):
        agg[k] = float(np.mean([res[x][k] for x in names]))
    for a, b in (("asm_o75_avg", "null_o75_avg"), ("asm_o75_fit", "null_o75_fit"),
                 ("null_o75_avg", "pool_o75_avg"), ("asm_o25_avg", "null_o25_avg")):
        agg[f"paired_{a}_vs_{b}"] = I.paired(np.array([res[x][a] for x in names]),
                                             np.array([res[x][b] for x in names]), names=names)
    print(json.dumps({k: v for k, v in agg.items() if not k.startswith("paired")}, indent=1))
    for k, v in agg.items():
        if k.startswith("paired"):
            print(f"{k}: d={v['mean_diff']:+.3f} CI{v['ci95']} W/L {v['n_better']}/{v['n_worse']}")
    I.write("adv_asmnull", dict(agg=agg, rows=res))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 24)
