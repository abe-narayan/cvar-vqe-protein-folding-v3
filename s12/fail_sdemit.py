"""FAIL18 forensics, step 16: emitted-structure test of the sd-weighted objective.

E15 found that weighting the shipped Bayes-risk score by sd^-2 (the distogram's OWN
per-pair uncertainty -- fully deployable, no new model, no native information) lifts
FAIL18 band recall 0/18 -> 9/18 and top-75 best 4.677 -> 3.868 A at ZERO cost to the
other 108 (1.911 -> 1.912).  Does it survive the coordinate average + projection?

Baseline emitted values are re-used from s12/results/fail_emit.json (arm shipped75),
which is the identical selection, so only the new arm is projected here.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)
from s12 import instrument as I


def main():
    prev = {r["pdb"]: r for r in json.load(open(os.path.join(ROOT, "s12", "results", "fail_emit.json")))}
    tg = [t for t in I.targets() if t["pdb"] in prev]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]; nat = u["nat_ca"]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        grid = np.asarray(dg["grid"], float); risk = np.asarray(dg["risk"], float)
        sd = np.asarray(dg["sd"], float)
        D = I.pair_dists(Wc, i, j).astype(np.float32).astype(float)
        g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
        R = risk[np.arange(risk.shape[0])[None, :], g]
        w = 1.0 / np.maximum(sd, 1e-3) ** 2
        sc = (R * w[None, :]).sum(1) / w.sum()
        sub = np.argsort(sc, kind="stable")[:I.M]
        C, _ = I.coordinate_average(Wc[sub])
        o = I.project(C, seq, fold)
        rows.append(dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                         invsd2_ca=float(I.ca_rmsd(o["ca"], nat)),
                         invsd2_fit=float(I.ca_rmsd(o["fit_ca"], nat)),
                         invsd2_top75_best=float(rr[sub].min()),
                         invsd2_top75_mean=float(rr[sub].mean()),
                         base_ca=prev[pdb]["arms"]["shipped75"]["ca"],
                         base_top75_mean=prev[pdb]["arms"]["shipped75"]["set_mean"]))
        print(f"  {k+1}/{len(tg)} {pdb} F={int(rows[-1]['fail18'])} "
              f"base={rows[-1]['base_ca']:.2f} invsd2={rows[-1]['invsd2_ca']:.2f} "
              f"free={I.free_gb():.1f}", flush=True)
        I.write("fail_sdemit", rows)

    f = [r for r in rows if r["fail18"]]; c = [r for r in rows if not r["fail18"]]
    print(f"\n{'group':10s}{'base':>8s}{'invsd2':>9s}{'d':>8s}{'W/L':>8s}{'ci95':>20s}")
    summ = {}
    for tag, g in (("FAIL18", f), ("MATCH18", c), ("all36", rows)):
        a = np.array([r["invsd2_ca"] for r in g]); b = np.array([r["base_ca"] for r in g])
        st = I.paired(a, b, names=[r["pdb"] for r in g])
        summ[tag] = st
        print(f"{tag:10s}{b.mean():8.3f}{a.mean():9.3f}{st['mean_diff']:8.3f}"
              f"{st['n_better']:4d}/{st['n_worse']:<3d} [{st['ci95'][0]:7.3f},{st['ci95'][1]:7.3f}]")
    print("\nper-target FAIL18:")
    for r in sorted(f, key=lambda r: r["pdb"]):
        print(f"  {r['pdb']:6s} base={r['base_ca']:5.2f} invsd2={r['invsd2_ca']:5.2f} "
              f"d={r['invsd2_ca']-r['base_ca']:+5.2f}  set_mean {r['base_top75_mean']:5.2f}"
              f" -> {r['invsd2_top75_mean']:5.2f}")
    I.write("fail_sdemit_summary", summ)


if __name__ == "__main__":
    main()
