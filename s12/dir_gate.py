"""D2. The PRE-REGISTERED accuracy gate, set on the ORACLE sign before any model exists.

Two sweeps on all 126 targets, both ORACLE/DIAGNOSTIC:
  * correction magnitude  delta in {0.25 .. 3.0}     (how far to move a distance)
  * sign accuracy         rate of flipped signs in {0 .. 0.5}   (0.5 == random == the null)

Reported on the COORDINATE AVERAGE, which obj_FINDINGS 0b validated as a near-exact proxy
for the emitted structure (r = 0.995, a near-constant +0.16 A from projection).  The
headline cells are re-run through the real projection by dir_gateproj.py.

The reference for every cell is `pt` (delta = 0), i.e. the SAME point-estimate objective
with no sign correction -- so the number reported is the value of the sign channel alone,
not the value of replacing the Bayes risk with an L1 point estimate.
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

DELTAS = (0.25, 0.5, 1.0, 1.5, 2.0, 3.0)
RATES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
SEEDS = 3


def main():
    tg = I.targets()
    rows = []
    t0 = time.time()
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        if I.free_gb() < 1.5:
            print("  waiting on memory", flush=True); time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue, grid = d["exp"], d["dtrue"], d["grid"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        cells = {}
        cells["bayes"] = OC.emit(d, OC.score_bayes(d["D"], d["risk"], grid), lam=None)["avg_rmsd"]
        cells["pt"] = DC.emit_target(d, exp, lam=None)["avg_rmsd"]
        cells["o_true"] = DC.emit_target(d, dtrue, lam=None)["avg_rmsd"]
        for D in DELTAS:
            for r in RATES:
                vals = []
                for s in range(1 if r == 0.0 else SEEDS):
                    rng = np.random.default_rng(1000 * k + 17 * s + int(100 * r))
                    sc = DC.corrupt_sign(s_true, r, rng)
                    vals.append(DC.emit_target(d, DC.sign_target(exp, sc, D), lam=None)["avg_rmsd"])
                cells[f"d{D}_r{r}"] = float(np.mean(vals))
        rows.append(dict(pdb=pdb, n=d["n"], fold=d["fold"], fail18=pdb in I.FAIL18,
                         pool_best=float(d["rr"].min()), cells=cells))
        if (k + 1) % 5 == 0 or k == 0:
            print(f"  {k+1}/126 {pdb} pt={cells['pt']:.2f} d2r0={cells['d2.0_r0.0']:.2f} "
                  f"d2r0.5={cells['d2.0_r0.5']:.2f} true={cells['o_true']:.2f} "
                  f"[{time.time()-t0:.0f}s free={I.free_gb():.1f}]", flush=True)
            I.write("dir_gate", rows)
    I.write("dir_gate", rows)

    # ------------------------------------------------------------------ report
    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    summ = {}
    for g, rs in grp.items():
        pt = np.array([r["cells"]["pt"] for r in rs])
        tab = {"bayes": float(np.mean([r["cells"]["bayes"] for r in rs])),
               "pt": float(pt.mean()),
               "o_true": float(np.mean([r["cells"]["o_true"] for r in rs]))}
        for D in DELTAS:
            for r_ in RATES:
                x = np.array([r["cells"][f"d{D}_r{r_}"] for r in rs])
                st = I.paired(x, pt)
                tab[f"d{D}_r{r_}"] = dict(mean=float(x.mean()), d=st["mean_diff"],
                                          ci=st["ci95"], wl=[st["n_better"], st["n_worse"]],
                                          drop10=st["drop_top10_mean_diff"])
        summ[g] = tab
    I.write("dir_gate_summary", summ)

    for g in grp:
        print(f"\n== {g} ==  pt={summ[g]['pt']:.3f}  bayes={summ[g]['bayes']:.3f}  "
              f"o_true={summ[g]['o_true']:.3f}   (delta vs pt, coordinate average)")
        print("delta \\ flip-rate " + "".join(f"{r:>9.2f}" for r in RATES))
        for D in DELTAS:
            print(f"  {D:<16.2f}" + "".join(f"{summ[g][f'd{D}_r{r}']['d']:9.3f}" for r in RATES))


if __name__ == "__main__":
    main()
