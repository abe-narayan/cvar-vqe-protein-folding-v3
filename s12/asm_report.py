"""ASM-8: pull the four result files together into the tables the findings need."""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I                        # noqa: E402
R = os.path.join(ROOT, "s12", "results")


def load(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    tg = I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.array([t["fold"] for t in tg])
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    rig = load("asm_rigid_ladder.json")
    ch = load("asm_chain2_M150.json")
    nul = load("asm_null_rigid.json")
    dep = load("asm_deploy_N200.json")
    out = {}

    def col(d, key, default=np.nan):
        return np.array([d[p].get(key, default) if p in d else default for p in pdbs], float)

    if ch:
        print("=== CHAIN-CONSISTENT assembly (ideal-geometry torsion concatenation), M=150 ===")
        print(f"{'arm':10s} {'k=1':>7s} {'k=2':>7s} {'k1 F18':>7s} {'k2 F18':>7s} {'k2 o108':>8s}")
        rows = {}
        for arm in ("real", "rama", "blos500", "rand500"):
            a1 = col(ch, f"k1_{arm}"); a2 = col(ch, f"k2_{arm}")
            rows[arm] = (a1, a2)
            print(f"{arm:10s} {np.nanmean(a1):7.3f} {np.nanmean(a2):7.3f} "
                  f"{np.nanmean(a1[f18]):7.3f} {np.nanmean(a2[f18]):7.3f} {np.nanmean(a2[~f18]):8.3f}")
        out["chain"] = {a: [float(np.nanmean(v[0])), float(np.nanmean(v[1]))] for a, v in rows.items()}
        print("\npaired real vs rama (capacity null), chain-consistent k=2:")
        st = I.paired(rows["real"][1], rows["rama"][1], folds=folds, names=pdbs)
        print(f"  d={st['mean_diff']:+.3f} CI{[round(x,3) for x in st['ci95']]} "
              f"W/L {st['n_better']}/{st['n_worse']} drop10 {st['drop_top10_mean_diff']:+.3f} "
              f"per-fold {[round(v,3) for v in st['per_fold'].values()]}")
        out["chain_real_vs_rama_k2"] = st
        print("paired blos500 vs rand500 (retrieval null), chain-consistent k=2:")
        st2 = I.paired(rows["blos500"][1], rows["rand500"][1], folds=folds, names=pdbs)
        print(f"  d={st2['mean_diff']:+.3f} CI{[round(x,3) for x in st2['ci95']]} "
              f"W/L {st2['n_better']}/{st2['n_worse']}")
        out["chain_blos_vs_rand_k2"] = st2
        if rig:
            rk2 = np.array([rig["arms"]["full_lmin4"]["2"]["per_target"][p] for p in pdbs], float)
            rk1 = np.array([rig["arms"]["full_lmin4"]["1"]["per_target"][p] for p in pdbs], float)
            print(f"\nPLACEMENT COST (chain-consistent minus rigid, full library):")
            print(f"  k=1  rigid {np.nanmean(rk1):.3f}  chain {np.nanmean(rows['real'][0]):.3f}"
                  f"   cost {np.nanmean(rows['real'][0]-rk1):+.3f}   (= ideal-rebuild penalty)")
            print(f"  k=2  rigid {np.nanmean(rk2):.3f}  chain {np.nanmean(rows['real'][1]):.3f}"
                  f"   cost {np.nanmean(rows['real'][1]-rk2):+.3f}")
            out["placement_cost"] = {"k1": float(np.nanmean(rows['real'][0] - rk1)),
                                     "k2": float(np.nanmean(rows['real'][1] - rk2))}

    if nul:
        print("\n=== RIGID ladder on the capacity-null bank (matched size) ===")
        for k in (1, 2, 3, 4):
            v = col(nul, f"k{k}")
            r = np.array([rig["arms"]["full_lmin4"][str(k)]["per_target"][p] for p in pdbs], float)
            print(f"  k={k}  null {np.nanmean(v):.3f}   real {np.nanmean(r):.3f}   "
                  f"real-minus-null {np.nanmean(r-v):+.3f}")
        out["null_rigid"] = {k: float(np.nanmean(col(nul, f"k{k}"))) for k in (1, 2, 3, 4)}

    if dep:
        print("\n=== DEPLOYABLE two-piece assembly (no native in the decision) ===")
        keys = ["pool_argmin", "pool1_argmin", "asm_argmin", "mix_argmin", "asm_avg75",
                "mix_avg75", "pool_best", "pool1_best", "asm_best", "mix_best",
                "asm_best_pctile", "asm_top75_best", "mix_top75_best"]
        for k in keys:
            v = col(dep, k)
            if np.all(np.isnan(v)):
                continue
            print(f"  {k:16s} {np.nanmean(v):7.3f}   FAIL18 {np.nanmean(v[f18]):7.3f}  "
                  f"other {np.nanmean(v[~f18]):7.3f}")
        base = col(dep, "pool_argmin")
        out["deploy"] = {}
        for k in ("asm_argmin", "mix_argmin", "asm_avg75", "mix_avg75"):
            a = col(dep, k)
            if np.all(np.isnan(a)):
                continue
            st = I.paired(a, base, folds=folds, names=pdbs)
            out["deploy"][k] = st
            print(f"\n  {k} vs shipped pool argmin (3.454): d={st['mean_diff']:+.3f} "
                  f"CI{[round(x,3) for x in st['ci95']]} W/L {st['n_better']}/{st['n_worse']} "
                  f"drop10 {st['drop_top10_mean_diff']:+.3f} "
                  f"per-fold {[round(v,3) for v in st['per_fold'].values()]}")
            d = a - base
            print(f"     FAIL18 {d[f18].mean():+.3f}   other108 {d[~f18].mean():+.3f}")
    I.write("asm_report", out)
    return out


if __name__ == "__main__":
    main()
