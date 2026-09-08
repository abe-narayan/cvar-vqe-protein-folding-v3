"""D2c. Merge the two ORACLE magnitude sweeps into the PRE-REGISTERED accuracy gate.

For each correction magnitude delta, linearly interpolate the emitted-RMSD gain as a
function of sign accuracy (= 1 - flip rate) and report the accuracy needed to reach
0 / 0.05 / 0.10 / 0.20 / 0.30 A of gain over the `pt` reference.  The magnitude used
downstream is chosen HERE, on the ORACLE sign, at the flip rate matching the head's
measured accuracy -- before any trained model's emitted RMSD is looked at.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I

RATES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
TARGETS_GAIN = (0.0, -0.05, -0.10, -0.20, -0.30, -0.50)


def main():
    g1 = json.load(open(os.path.join(ROOT, "s12", "results", "dir_gate_summary.json")))
    g2 = json.load(open(os.path.join(ROOT, "s12", "results", "dir_gate2_summary.json")))
    hd = json.load(open(os.path.join(ROOT, "s12", "results", "dir_head.json")))
    acc_head = {g: hd["groups"][g]["acc_full"] for g in ("tuning126", "fail18", "other108")}

    out = {"acc_head": acc_head, "curves": {}, "gate": {}}
    for g in ("all126", "fail18", "other108"):
        cur = {}
        for D in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0):
            cur[str(D)] = [g1[g][f"d{D}_r{r}"]["d"] for r in RATES]
        for D in (4.0, 6.0):
            cur[str(D)] = [g2[g][f"d{D}_r{r}"]["d"] for r in RATES]
        for nm in ("lin", "linsgn"):
            cur[nm] = [g2[g][f"{nm}_r{r}"]["d"] for r in RATES]
        out["curves"][g] = {"rates": list(RATES), "acc": [1 - r for r in RATES], "gain": cur}
        acc = np.array([1 - r for r in RATES])[::-1]        # ascending accuracy 0.5..1.0
        gate = {}
        for D, v in cur.items():
            y = np.array(v)[::-1]                            # gain at ascending accuracy
            req = {}
            for t in TARGETS_GAIN:
                k = np.where(y <= t)[0]
                if not len(k):
                    req[f"{t:+.2f}"] = None
                    continue
                k0 = k[0]
                if k0 == 0:
                    req[f"{t:+.2f}"] = float(acc[0])
                else:
                    a0, a1, y0, y1 = acc[k0 - 1], acc[k0], y[k0 - 1], y[k0]
                    req[f"{t:+.2f}"] = float(a0 + (t - y0) * (a1 - a0) / (y1 - y0))
            gate[D] = req
        out["gate"][g] = gate

    # ------------------------------------------------------------- pick delta
    a = acc_head["tuning126"]
    accs = np.array([1 - r for r in RATES])[::-1]
    best, bestv = None, 1e9
    for D, v in out["curves"]["all126"]["gain"].items():
        if D in ("lin", "linsgn"):
            continue
        y = np.array(v)[::-1]
        val = float(np.interp(a, accs, y))
        out["curves"]["all126"].setdefault("pred_at_head_acc", {})[D] = val
        if val < bestv:
            best, bestv = D, val
    pick = dict(delta=float(best), predicted_gain_all126=bestv, acc_head=a,
                rule="argmin over the ORACLE gate curve at the head's measured accuracy, "
                     "chosen before any trained model's emitted RMSD was inspected")
    out["pick"] = pick
    I.write("dir_gate_pick", pick)
    I.write("dir_gate_full", out)

    for g in out["gate"]:
        print(f"\n== {g}: sign ACCURACY required for a given gain over `pt` "
              f"(ORACLE sign, random flips, coordinate average) ==")
        print(f"{'delta':>8s}" + "".join(f"{t:>+9.2f}" for t in TARGETS_GAIN))
        for D in ("0.25", "0.5", "1.0", "1.5", "2.0", "3.0", "4.0", "6.0", "lin", "linsgn"):
            r = out["gate"][g][D]
            print(f"{D:>8s}" + "".join(
                ("   n/a   " if r[f'{t:+.2f}'] is None else f"{r[f'{t:+.2f}']:9.3f}")
                for t in TARGETS_GAIN))
    print(f"\nhead accuracy: {acc_head}")
    print(f"predicted gain at that accuracy, all126: "
          f"{ {k: round(v,3) for k,v in out['curves']['all126']['pred_at_head_acc'].items()} }")
    print(f"PICK delta = {pick['delta']}  (predicted {pick['predicted_gain_all126']:+.3f} A)")


if __name__ == "__main__":
    main()
