"""D4b. The deployable arms across a magnitude ladder, on the coordinate-average proxy.

Cheap companion to dir_emit.py: shows whether the head's effect is a knife-edge in the
correction magnitude Δ, and gives every arm at every Δ (the projected run only affords
one Δ).  Same nulls.  obj_FINDINGS 0b validates avg_rmsd as a proxy (r = 0.995).
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)

from s12 import instrument as I
from s12 import obj_common as OC
from s12 import dir_common as DC

PRED = os.path.join(ROOT, "s12", "cache", "dir_head_pred.npz")
LADDER = (1.0, 2.0, 3.0, 4.0, 6.0)


def main():
    hd = json.load(open(os.path.join(ROOT, "s12", "results", "dir_head.json")))
    acc = hd["per_target"]
    Z = np.load(PRED)
    rows, t0 = [], time.time()
    for k, t in enumerate(I.targets()):
        pdb = t["pdb"]
        while I.free_gb() < 1.5:
            print("  waiting on memory", flush=True); time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue = d["exp"], d["dtrue"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        p = {a: np.asarray(Z[f"{pdb}/{a}"], float) for a in
             ("full", "sep", "shellmaj", "const", "shuf")}
        sg = {a: np.where(p[a] > 0.5, 1.0, -1.0) for a in p}
        cells = {"bayes": OC.emit(d, OC.score_bayes(d["D"], d["risk"], d["grid"]), lam=None)["avg_rmsd"],
                 "pt": DC.emit_target(d, exp, lam=None)["avg_rmsd"]}
        a_head = float(acc[pdb]["full"])
        for D in LADDER:
            cells[f"o_sign@{D}"] = DC.emit_target(d, DC.sign_target(exp, s_true, D), lam=None)["avg_rmsd"]
            cells[f"head@{D}"] = DC.emit_target(d, DC.sign_target(exp, sg["full"], D), lam=None)["avg_rmsd"]
            cells[f"head_soft@{D}"] = DC.emit_target(d, DC.sign_target(exp, 2 * p["full"] - 1, D), lam=None)["avg_rmsd"]
            for a in ("sep", "shellmaj", "const", "shuf"):
                cells[f"{a}@{D}"] = DC.emit_target(d, DC.sign_target(exp, sg[a], D), lam=None)["avg_rmsd"]
            tm = 1.0 if s_true.mean() > 0 else -1.0
            cells[f"o_tmaj@{D}"] = DC.emit_target(d, DC.sign_target(exp, np.full(exp.shape, tm), D), lam=None)["avg_rmsd"]
            hm = 1.0 if sg["full"].mean() > 0 else -1.0
            cells[f"head_tmaj@{D}"] = DC.emit_target(d, DC.sign_target(exp, np.full(exp.shape, hm), D), lam=None)["avg_rmsd"]
            v = []
            for s in range(3):
                rng = np.random.default_rng(7000 + 13 * k + s)
                v.append(DC.emit_target(d, DC.sign_target(exp, DC.corrupt_sign(s_true, 1 - a_head, rng), D),
                                        lam=None)["avg_rmsd"])
            cells[f"rand_acc@{D}"] = float(np.mean(v))
        rows.append(dict(pdb=pdb, n=d["n"], fold=d["fold"], fail18=pdb in I.FAIL18,
                         acc_head=a_head, cells=cells))
        if (k + 1) % 20 == 0 or k == 0:
            print(f"  {k+1}/126 [{time.time()-t0:.0f}s free={I.free_gb():.1f}]", flush=True)
            I.write("dir_emitavg", rows)
    I.write("dir_emitavg", rows)

    arms = ["bayes"] + [f"{a}@{D}" for D in LADDER for a in
                        ("o_sign", "head", "head_soft", "rand_acc", "sep", "shellmaj",
                         "const", "shuf", "o_tmaj", "head_tmaj")]
    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    summ = {}
    for g, rs in grp.items():
        folds = np.array([r["fold"] for r in rs])
        pt = np.array([r["cells"]["pt"] for r in rs])
        by = np.array([r["cells"]["bayes"] for r in rs])
        tab = {"pt": float(pt.mean()), "bayes": float(by.mean())}
        for a in arms:
            x = np.array([r["cells"][a] for r in rs])
            s1 = I.paired(x, pt, folds=folds); s2 = I.paired(x, by, folds=folds)
            tab[a] = dict(mean=float(x.mean()), d_pt=s1["mean_diff"], ci_pt=s1["ci95"],
                          wl_pt=[s1["n_better"], s1["n_worse"]],
                          d_bayes=s2["mean_diff"], ci_bayes=s2["ci95"],
                          wl_bayes=[s2["n_better"], s2["n_worse"]],
                          drop10=s2["drop_top10_mean_diff"], drop20=s2["drop_top20_mean_diff"],
                          per_fold=s2.get("per_fold"))
        summ[g] = tab
    I.write("dir_emitavg_summary", summ)

    for g in grp:
        t = summ[g]
        print(f"\n== {g} == pt={t['pt']:.3f} bayes={t['bayes']:.3f}  (Δ vs BAYES, coordinate average)")
        hdr = ("o_sign", "head", "head_soft", "rand_acc", "sep", "shellmaj", "const", "shuf",
               "o_tmaj", "head_tmaj")
        print(f"{'Δ':>5s}" + "".join(f"{h:>11s}" for h in hdr))
        for D in LADDER:
            print(f"{D:5.1f}" + "".join(f"{t[f'{h}@{D}']['d_bayes']:11.3f}" for h in hdr))
        print(f"{'W/L':>5s}" + "".join(
            f"{t[f'{h}@4.0']['wl_bayes'][0]:5d}/{t[f'{h}@4.0']['wl_bayes'][1]:<5d}" for h in hdr))


if __name__ == "__main__":
    main()
