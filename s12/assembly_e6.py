"""E6 (ORACLE/DIAGNOSTIC) -- is the better assembled pool worth anything if selection were SOLVED?

For the `random`-shortlist assembled space and for the retrieval K=500 pool, take the ORACLE top-25
and top-75 by true CA-RMSD, coordinate-average and project, and report the emitted structure. This is
the C1 comparison (oracle pool averaging 1.925 / oracle top-25 averaging 1.644 on retrieval) recomputed
on the assembled space. Nothing here is deployable; it prices the CEILING that a solved selector
would reach on each candidate space.
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_common as AC
from s12 import assembly_e2 as E2
from s12 import assembly_e3 as E3

KEY = os.environ.get("ASM_KEY", "random")
KTOP = 20


def run_target(t):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    out = os.path.join(I.CACHE, f"asm_e6_{KEY}_{pdb}.json")
    if os.path.exists(out):
        return json.load(open(out))
    u = I.load_univ(pdb); nat = u["nat_ca"]; p = I.pool_idx(u)
    E3.KEY = KEY; E3._TOP = {}
    K = E2.target_keys(t)
    allPHI, allPSI = [], []
    for comp in E3.structures(n):
        PHI, PSI, grids, tops, junc = E3.enumerate_structure(t, K, comp, KTOP)
        allPHI.append(PHI); allPSI.append(PSI)
    PHI = np.concatenate(allPHI); PSI = np.concatenate(allPSI)
    CA = AC.build_many(PHI, PSI); rr = I.kabsch_rmsd_batch(CA, nat)
    res = {"pdb": pdb, "n": n, "fold": fold, "asm_best": float(rr.min()), "retr_best": float(u["rr"][p].min())}
    for tag, W, r in (("asm", CA, rr), ("retr", u["W"][p], u["rr"][p])):
        for m in (25, 75):
            sub = np.argsort(r, kind="stable")[:m]
            C, _ = I.coordinate_average(W[sub])
            res[f"{tag}_oracle{m}_avg"] = I.ca_rmsd(C, nat)
            pr = I.project(C, seq, fold)
            res[f"{tag}_oracle{m}_fit"] = I.ca_rmsd(pr["fit_ca"], nat)
            res[f"{tag}_oracle{m}_ca"] = I.ca_rmsd(pr["ca"], nat)
    json.dump(res, open(out, "w"), default=float)
    return json.load(open(out))


def _init():
    os.environ["OMP_NUM_THREADS"] = "1"
    try:
        import torch; torch.set_num_threads(1)
    except Exception:
        pass


if __name__ == "__main__":
    import multiprocessing as mp
    tg = I.targets(); rows = []
    print("key", KEY, "free GB", I.free_gb(), flush=True)
    with mp.Pool(2, initializer=_init) as pool:
        for c, r in enumerate(pool.imap_unordered(run_target, tg)):
            rows.append(r)
            if c % 15 == 0:
                print(f"[{c+1}/126] {r['pdb']} free={I.free_gb():.1f}", flush=True)
    rows.sort(key=lambda r: r["pdb"]); pdbs = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    agg = {"key": KEY, "n": len(rows)}
    keys = [k for k in rows[0] if k not in ("pdb", "n", "fold")]
    for k in keys:
        agg[k] = AC.group_means([r[k] for r in rows], pdbs)
    for m in (25, 75):
        for w in ("avg", "fit"):
            agg[f"paired_asm_vs_retr_oracle{m}_{w}"] = I.paired([r[f"asm_oracle{m}_{w}"] for r in rows],
                                                                [r[f"retr_oracle{m}_{w}"] for r in rows], folds=folds, names=pdbs)
    agg["per_target"] = rows
    print(I.write(f"assembly_e6_oracle_{KEY}", agg))
    for k in keys:
        print(k, {a: round(b, 3) for a, b in agg[k].items()})
    for k, v in agg.items():
        if k.startswith("paired"):
            print(k, round(v["mean_diff"], 3), [round(x, 3) for x in v["ci95"]], f"{v['n_better']}W/{v['n_worse']}L")
