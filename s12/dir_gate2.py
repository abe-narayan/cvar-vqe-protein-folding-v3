"""D2b. Extension of the magnitude sweep, and the LINEAR LIMIT of the sign objective.

dir_gate.py found the ORACLE-sign gain still rising at delta = 3.0, so the sweep is
extended to 4 and 6, and the delta -> infinity limit is scored directly.  As delta grows,

    |D_p - (E[d]_p + delta * s_p)|  ->  delta - s_p (D_p - E[d]_p)

so the objective becomes the pure sign-matching linear score

    score(cand) = - mean_p s_p (D_p - E[d]_p)                        ("lin")

which uses ONLY the direction of each pair error and the candidate's signed deviation from
the predicted distance.  A magnitude-free variant that also throws away the candidate's
deviation size is

    score(cand) = - mean_p s_p * sign(D_p - E[d]_p)                  ("linsgn")

Same corruption rates, same coordinate-average proxy, all 126.
"""
from __future__ import annotations
import os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)

from s12 import instrument as I
from s12 import obj_common as OC
from s12 import dir_common as DC

DELTAS = (4.0, 6.0)
RATES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
SEEDS = 3


def lin_score(D, exp, s, hard=False):
    dev = np.asarray(D, float) - np.asarray(exp, float)[None, :]
    if hard:
        dev = np.sign(dev)
    return -(dev * np.asarray(s, float)[None, :]).mean(1)


def main():
    tg = I.targets()
    rows, t0 = [], time.time()
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        while I.free_gb() < 1.5:
            print("  waiting on memory", flush=True); time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue = d["exp"], d["dtrue"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        cells = {"pt": DC.emit_target(d, exp, lam=None)["avg_rmsd"]}
        for D in DELTAS:
            for r in RATES:
                v = []
                for s in range(1 if r == 0.0 else SEEDS):
                    rng = np.random.default_rng(1000 * k + 17 * s + int(100 * r))
                    v.append(DC.emit_target(d, DC.sign_target(exp, DC.corrupt_sign(s_true, r, rng), D),
                                            lam=None)["avg_rmsd"])
                cells[f"d{D}_r{r}"] = float(np.mean(v))
        for nm, hard in (("lin", False), ("linsgn", True)):
            for r in RATES:
                v = []
                for s in range(1 if r == 0.0 else SEEDS):
                    rng = np.random.default_rng(1000 * k + 17 * s + int(100 * r))
                    sc = DC.corrupt_sign(s_true, r, rng)
                    v.append(OC.emit(d, lin_score(d["D"], exp, sc, hard), lam=None)["avg_rmsd"])
                cells[f"{nm}_r{r}"] = float(np.mean(v))
        rows.append(dict(pdb=pdb, n=d["n"], fold=d["fold"], fail18=pdb in I.FAIL18, cells=cells))
        if (k + 1) % 10 == 0 or k == 0:
            print(f"  {k+1}/126 {pdb} pt={cells['pt']:.2f} d4r0={cells['d4.0_r0.0']:.2f} "
                  f"lin={cells['lin_r0.0']:.2f} [{time.time()-t0:.0f}s free={I.free_gb():.1f}]",
                  flush=True)
            I.write("dir_gate2", rows)
    I.write("dir_gate2", rows)

    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    summ = {}
    arms = [f"d{D}" for D in DELTAS] + ["lin", "linsgn"]
    for g, rs in grp.items():
        pt = np.array([r["cells"]["pt"] for r in rs])
        tab = {"pt": float(pt.mean())}
        for a in arms:
            for r_ in RATES:
                x = np.array([r["cells"][f"{a}_r{r_}"] for r in rs])
                st = I.paired(x, pt)
                tab[f"{a}_r{r_}"] = dict(mean=float(x.mean()), d=st["mean_diff"], ci=st["ci95"],
                                         wl=[st["n_better"], st["n_worse"]],
                                         drop10=st["drop_top10_mean_diff"])
        summ[g] = tab
        print(f"\n== {g} ==  pt={tab['pt']:.3f}  (delta vs pt, coordinate average)")
        print("arm \\ flip-rate  " + "".join(f"{r:>9.2f}" for r in RATES))
        for a in arms:
            print(f"  {a:<14s}" + "".join(f"{tab[f'{a}_r{r}']['d']:9.3f}" for r in RATES))
    I.write("dir_gate2_summary", summ)


if __name__ == "__main__":
    main()
