"""D5c. Last chance for a positive: apply the head's sign only where the head is CONFIDENT.

The head's per-pair accuracy rises steeply with |p - 0.5|; the ladder result may simply be
that the low-confidence 2/3 of the pairs inject a coherent wrong field.  Arms, on the
coordinate-average proxy, at delta in {3, 4}:

  head_q      shift only the top q fraction of pairs by head confidence |p - 0.5|
  rand_q      the SAME gate set, but ORACLE signs corrupted to the head's CONDITIONAL
              accuracy on that set (3 seeds) -- null (c) matched inside the gate
  conf_acc    the head's realised accuracy inside each gate (reported, not an arm)
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")

from s12 import instrument as I
from s12 import obj_common as OC
from s12 import dir_common as DC

QS = (0.10, 0.25, 0.50, 1.00)
DS = (3.0, 4.0)


def main():
    Z = np.load(os.path.join(ROOT, "s12", "cache", "dir_head_pred.npz"))
    rows, t0 = [], time.time()
    for k, t in enumerate(I.targets()):
        pdb = t["pdb"]
        while I.free_gb() < 1.5:
            print("  waiting on memory", flush=True); time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue = d["exp"], d["dtrue"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        p = np.asarray(Z[f"{pdb}/full"], float)
        sh = np.where(p > 0.5, 1.0, -1.0)
        conf = np.abs(p - 0.5)
        order = np.argsort(-conf, kind="stable")
        cells = {"pt": DC.emit_target(d, exp, lam=None)["avg_rmsd"],
                 "bayes": OC.emit(d, OC.score_bayes(d["D"], d["risk"], d["grid"]), lam=None)["avg_rmsd"]}
        accs = {}
        for q in QS:
            m = np.zeros(len(p), bool); m[order[:max(1, int(round(q * len(p))))]] = True
            a = float((sh[m] == s_true[m]).mean())
            accs[f"q{q}"] = a
            for D in DS:
                s = np.where(m, sh, 0.0)
                cells[f"head_q{q}@{D}"] = DC.emit_target(d, DC.sign_target(exp, s, D), lam=None)["avg_rmsd"]
                v = []
                for sd_ in range(3):
                    rng = np.random.default_rng(9000 + 31 * k + sd_)
                    sr = s_true.copy()
                    idx = np.where(m)[0]
                    nf = int(round((1 - a) * len(idx)))
                    if nf:
                        sr[rng.choice(idx, nf, replace=False)] *= -1
                    v.append(DC.emit_target(d, DC.sign_target(exp, np.where(m, sr, 0.0), D),
                                            lam=None)["avg_rmsd"])
                cells[f"rand_q{q}@{D}"] = float(np.mean(v))
        rows.append(dict(pdb=pdb, fold=t["fold"], fail18=pdb in I.FAIL18, accs=accs, cells=cells))
        if (k + 1) % 25 == 0 or k == 0:
            print(f"  {k+1}/126 [{time.time()-t0:.0f}s free={I.free_gb():.1f}]", flush=True)
            I.write("dir_conf", rows)
    I.write("dir_conf", rows)

    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    out = {}
    for g, rs in grp.items():
        by = np.array([r["cells"]["bayes"] for r in rs])
        folds = np.array([r["fold"] for r in rs])
        tab = {"_bayes": float(by.mean()), "_pt_d": float(np.mean([r["cells"]["pt"] for r in rs]) - by.mean()),
               "n": len(rs)}
        for q in QS:
            tab[f"acc_q{q}"] = float(np.mean([r["accs"][f"q{q}"] for r in rs]))
            for D in DS:
                for a in ("head", "rand"):
                    x = np.array([r["cells"][f"{a}_q{q}@{D}"] for r in rs])
                    s = I.paired(x, by, folds=folds)
                    tab[f"{a}_q{q}@{D}"] = dict(mean=float(x.mean()), d=s["mean_diff"],
                                                ci=s["ci95"], wl=[s["n_better"], s["n_worse"]],
                                                drop10=s["drop_top10_mean_diff"])
        out[g] = tab
        print(f"\n== {g} (n={tab['n']}) bayes={tab['_bayes']:.3f}  pt d={tab['_pt_d']:+.3f} "
              f"(delta vs BAYES, coordinate average) ==")
        print(f"{'q':>6s}{'acc':>8s}" + "".join(f"{f'head@{D}':>11s}{f'rand@{D}':>11s}" for D in DS))
        for q in QS:
            print(f"{q:6.2f}{tab[f'acc_q{q}']:8.4f}" + "".join(
                f"{tab[f'head_q{q}@{D}']['d']:11.3f}{tab[f'rand_q{q}@{D}']['d']:11.3f}" for D in DS))
    I.write("dir_conf_summary", out)


if __name__ == "__main__":
    main()
