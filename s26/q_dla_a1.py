"""s26/q_dla_a1.py -- the DLA dimension of the ADAPT-grown circuit at every growth step, on the
126 A1 records (PREREG_A2 addendum 2).  Reads only the operator lists in
s26/results/a1/<pdb>.json (key adapt.*.ops); no native, no RMSD, no score.

    python s26/q_dla_a1.py            -> s26/results/q_dla_a1.json, s26/figures/a2_dla_grown_ladder.png
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s26.q_adapt import Pauli                                              # noqa: E402
from s26.q_dla import lie_closure, dim_so                                  # noqa: E402

A1 = os.path.join(HERE, "results", "a1")
OUT = os.path.join(HERE, "results", "q_dla_a1.json")
FIG = os.path.join(HERE, "figures", "a2_dla_grown_ladder.png")
RUNS = ("V_adam_best_zrank", "L2_adam_best_zrank", "V_lbfgs_zrank", "L2_lbfgs_zrank")


def ladder(ops):
    P = [Pauli.from_word(w) for w in ops]
    n = P[0].n
    out = []
    for k in range(n, len(P) + 1):
        c = lie_closure(P[:k])
        out.append(dict(P=k, dim=int(c["dim"]), exceeded=bool(c["exceeded_cap"])))
    return out


def main():
    from s24 import stats_lib as ST
    files = sorted(f for f in os.listdir(A1) if f.endswith(".json"))
    R = {"records": {}, "n_records": len(files)}
    if os.path.exists(OUT):
        try:
            R = json.load(open(OUT))["results"]
        except Exception:                                                # noqa: BLE001
            pass
    t0 = time.perf_counter()
    for i, f in enumerate(files):
        pdb = f[:-5]
        if pdb in R["records"]:
            continue
        r = json.load(open(os.path.join(A1, f)))
        rec = {"fold": r["fold"], "alpha": r["alpha"], "T": r["T"]}
        for key in RUNS:
            ad = r["adapt"][key]
            rec[key] = dict(ladder=ladder(ad["ops"]), n_appended=len(ad["sequence"]),
                            n_distinct=ad["n_distinct_ops"], stopped=ad["stopped"],
                            n_multiqubit=sum(1 for w in ad["sequence"] if Pauli.from_word(w).weight() >= 2))
        R["records"][pdb] = rec
        if (i + 1) % 10 == 0 or i + 1 == len(files):
            ST.save_atomic(OUT, dict(kind="property, DLA of ADAPT-grown sets on the A1 records, no RMSD",
                                     lane="Q", sprint=26, prereg="s26/PREREG_A2.md addendum 2",
                                     results=R), module_file=__file__)
            print(f"  [{i + 1}/{len(files)}] {pdb}  ({time.perf_counter() - t0:.0f} s)", flush=True)
    # ---- summary per (cell, run) at P = 7, 14, 21 (or final) --------------------------
    summ = {}
    for key in RUNS:
        for a in (1.0, 0.25):
            recs = [v for v in R["records"].values() if v["alpha"] == a]
            dims_final = np.array([v[key]["ladder"][-1]["dim"] for v in recs])
            d7 = np.array([v[key]["ladder"][0]["dim"] for v in recs])
            at = {}
            for P in (7, 14, 21):
                vals = [next((x["dim"] for x in v[key]["ladder"] if x["P"] == P), None) for v in recs]
                vals = [x for x in vals if x is not None]
                if vals:
                    at[P] = dict(n=len(vals), median=float(np.median(vals)), min=int(min(vals)), max=int(max(vals)))
            summ[f"{key}:alpha{a}"] = dict(
                n=len(recs), dim_at_P=at,
                final=dict(median=float(np.median(dims_final)), min=int(dims_final.min()), max=int(dims_final.max())),
                n_dim7_final=int((dims_final == 7).sum()),
                n_nothing_appended=int(sum(v[key]["n_appended"] == 0 for v in recs)),
                n_with_multiqubit=int(sum(v[key]["n_multiqubit"] > 0 for v in recs)),
                n_reach_so128=int((dims_final == dim_so(7)).sum()),
                dim7_iff_nothing_appended=bool(all((v[key]["ladder"][-1]["dim"] == 7) == (v[key]["n_multiqubit"] == 0) for v in recs)))
            print(f"  {key:22s} alpha={a}: n={len(recs)}  final dim median {np.median(dims_final):.0f} "
                  f"[{dims_final.min()}, {dims_final.max()}]  dim7 {int((dims_final == 7).sum())}  "
                  f"nothing appended {summ[f'{key}:alpha{a}']['n_nothing_appended']}  "
                  f"multi-qubit appended {summ[f'{key}:alpha{a}']['n_with_multiqubit']}  reach so(128) {int((dims_final == 8128).sum())}",
                  flush=True)
    R["summary"] = summ
    ST.save_atomic(OUT, dict(kind="property, DLA of ADAPT-grown sets on the A1 records, no RMSD",
                             lane="Q", sprint=26, prereg="s26/PREREG_A2.md addendum 2", results=R),
                   module_file=__file__)
    print("wrote", OUT)
    figure(R)
    return 0


def figure(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), facecolor="white", sharey=True)
    colors = {"V_adam_best_zrank": "#ff7f0e", "L2_adam_best_zrank": "#17becf",
              "V_lbfgs_zrank": "#d62728", "L2_lbfgs_zrank": "#9467bd"}
    for ax, a in zip(axes, (1.0, 0.25)):
        ax.set_facecolor("white")
        recs = [v for v in R["records"].values() if v["alpha"] == a]
        for key in RUNS:
            Ps = sorted({x["P"] for v in recs for x in v[key]["ladder"]})
            med, lo, hi, cnt = [], [], [], []
            for P in Ps:
                vals = [x["dim"] for v in recs for x in v[key]["ladder"] if x["P"] == P]
                med.append(np.median(vals)); lo.append(min(vals)); hi.append(max(vals)); cnt.append(len(vals))
            ax.plot(Ps, med, "-o", ms=3.5, color=colors[key], label=f"{key.replace('_zrank', '')}: median (min..max shaded)")
            ax.fill_between(Ps, lo, hi, color=colors[key], alpha=0.12)
        ax.axhline(dim_so(7), color="0.35", lw=1.2, ls="--", label="dim so(128) = 8128 (fixed ansatz from depth 2)")
        ax.axhline(7, color="0.6", lw=0.9, ls=":", label="dim 7 (abelian, the RY layer)")
        ax.set_yscale("log")
        ax.set_xlabel("parameters P (7 = the RY layer, then one string per growth step)")
        ax.set_title(f"alpha = {a}, T = 0.3  (n = {len(recs)} targets)", fontsize=9)
        ax.grid(True, which="major", color="0.9")
        if a == 1.0:
            ax.set_ylabel("dim of the Lie closure of the grown set")
        ax.legend(fontsize=6.6, loc="upper left")
    fig.suptitle("A2 per growth step, on the 126 A1 records: the algebra of the SET grows while the angles stay inert at alpha = 1 (L75)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG, dpi=190, facecolor="white")
    plt.close(fig)
    print("wrote", FIG)


if __name__ == "__main__":
    sys.exit(main())
