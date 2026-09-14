"""s26/q_figures.py -- the two slide figures for lane Q, from the results JSON only.

    python s26/q_figures.py            -> s26/figures/a2_dla_dimension.png, s26/figures/a4_variance_slopes.png

Both at 190 dpi on a white background.  A2 reads s26/results/q_dla.json; A4 reads
s26/results/q_var.json (a partial file is plotted with what it has and the title says so).
No number is typed in by hand.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

import matplotlib                                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                            # noqa: E402

DPI = 190
CELL_LABEL = {"linear_alpha1_T0": "alpha=1, T=0 (linear)", "cvar_alpha025_T0": "alpha=0.25, T=0",
              "cvar_alpha01_T0": "alpha=0.10, T=0", "deployed_a1_T03": "alpha=1, T=0.3 (deployed)",
              "deployed_a025_T03": "alpha=0.25, T=0.3 (deployed)"}


def fig_a2(path_json=os.path.join(RES, "q_dla.json"), out=os.path.join(FIG, "a2_dla_dimension.png")):
    R = json.load(open(path_json))["results"]
    fixed = R["fixed"]
    ns = sorted({v["n"] for v in fixed.values()})
    fig, ax = plt.subplots(figsize=(7.2, 4.6), facecolor="white")
    ax.set_facecolor("white")
    so = [2 ** (n - 1) * (2 ** n - 1) for n in ns]
    su = [4 ** n - 1 for n in ns]
    ax.plot(ns, su, "--", color="0.55", lw=1.2, label="dim su(2^n) = 4^n - 1")
    ax.plot(ns, so, "-", color="0.35", lw=1.6, label="dim so(2^n) = 2^(n-1)(2^n - 1)")
    styles = {1: dict(marker="v", color="#1f77b4"), 2: dict(marker="o", color="#d62728"),
              3: dict(marker="s", color="#2ca02c"), 4: dict(marker="^", color="#9467bd")}
    for L, st in styles.items():
        xs, ys = [], []
        for n in ns:
            k = f"n{n}_L{L}"
            if k in fixed and not fixed[k]["exceeded_cap"]:
                xs.append(n)
                ys.append(fixed[k]["dim"])
        if xs:
            ax.plot(xs, ys, linestyle="none", ms=7, mfc="white", mew=1.6,
                    label=f"fixed RY/CNOT ansatz, depth {L}", **st)
    pools = R.get("pools", {})
    for pn, st in (("V", dict(marker="D", color="#ff7f0e")), ("L2", dict(marker="x", color="#17becf"))):
        xs = sorted(int(k.split("n")[1]) for k in pools if k.startswith(pn + "_n"))
        ys = [pools[f"{pn}_n{n}"]["dim"] for n in xs]
        if xs:
            lab = ("Tang pool V (2n-2 strings), whole pool" if pn == "V"
                   else "2-local odd-Y pool, whole pool")
            ax.plot(xs, ys, linestyle=":", ms=7, mew=1.8, label=lab, **st)
    ad = R.get("adapt_sets", {})
    for key, v in ad.items():
        if v["pool"] == "L2" and v["alpha"] == 0.25 and v["optimiser"] == "adam_best":
            d = v["ladder"][-1]["dim"]
            ax.plot([v["n"]], [d], marker="*", ms=13, color="#8c564b", linestyle="none",
                    label=f"ADAPT-selected 21 strings (L2, alpha=0.25): {d}")
        if v["pool"] == "L2" and v["alpha"] == 1.0 and v["optimiser"] == "adam_best":
            d = v["ladder"][-1]["dim"]
            ax.plot([v["n"]], [d], marker="*", ms=13, color="#e377c2", linestyle="none",
                    label=f"ADAPT-selected 21 strings (L2, alpha=1): {d}")
    n7 = fixed.get("n7_L3", {}).get("dim")
    if n7:
        ax.annotate(f"deployed: n=7, depth 3\ndim = {n7}", xy=(7, n7), xytext=(4.1, 2.0e5),
                    fontsize=8.5, arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9))
    ax.set_yscale("log")
    ax.set_xlabel("qubits n")
    ax.set_ylabel("dim of the dynamical Lie algebra")
    ax.set_xticks(ns)
    ax.grid(True, which="major", color="0.9")
    ax.legend(fontsize=7.4, loc="lower right", frameon=True)
    ax.set_title("A2: exact Lie closure on Pauli strings (== dense SVD rank, rtol 1e-10, at n = 4, 5)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    return out


def _slope(rows, matched_only=True):
    rows = [r for r in rows if (r.get("P_matched", True) or not matched_only)]
    if len(rows) < 3:
        return float("nan"), rows
    ns = np.array([r["n"] for r in rows], float)
    v = np.array([r["var_g0"] for r in rows], float)
    ok = v > 0
    if ok.sum() < 3:
        return float("nan"), rows
    return float(np.polyfit(ns[ok], np.log2(v[ok]), 1)[0]), rows


def fig_a4(path_json=os.path.join(RES, "q_var.json"), out=os.path.join(FIG, "a4_variance_slopes.png")):
    J = json.load(open(path_json))
    R = J["results"]
    cells = [c for c in CELL_LABEL if c in R.get("fixed", {})]
    grown_keys = sorted(R.get("grown", {}).keys())
    pools = sorted({k.split(":")[0] for k in grown_keys})
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), facecolor="white",
                             gridspec_kw=dict(width_ratios=[1.35, 1.0]))
    # left: slopes per cell
    ax = axes[0]
    ax.set_facecolor("white")
    width = 0.8 / (1 + len(pools))
    x = np.arange(len(cells))
    partial = False
    sl_fixed = []
    for c in cells:
        s, rows = _slope(R["fixed"][c]["rows"], matched_only=False)
        sl_fixed.append(s)
        if len(rows) < 7:
            partial = True
    boot = None
    bpath = os.path.join(RES, "q_var_boot.json")
    if os.path.exists(bpath):
        try:
            boot = json.load(open(bpath))["results"].get("ci")
        except Exception:                                                # noqa: BLE001
            boot = None

    def err(cell, arm):
        if not boot or cell not in boot or arm not in boot[cell] or not boot[cell][arm].get("ci95"):
            return None
        pt, (lo, hi) = boot[cell][arm]["point"], boot[cell][arm]["ci95"]
        return [[pt - lo], [hi - pt]]

    ef = [err(c, "fixed") for c in cells]
    ax.bar(x, sl_fixed, width, color="0.35", label="fixed RY/CNOT ansatz, depth 3 (P = 3n)")
    for xi, v, e in zip(x, sl_fixed, ef):
        if e:
            ax.errorbar([xi], [v], yerr=e, fmt="none", ecolor="black", capsize=3, lw=1)
    colors = {"V": "#ff7f0e", "L2": "#17becf"}
    for j, pn in enumerate(pools):
        sl = []
        for c in cells:
            key = f"{pn}:adam_best:{c}"
            if key in R["grown"]:
                s, rows = _slope(R["grown"][key]["rows"], matched_only=True)
                n_un = sum(1 for r in R["grown"][key]["rows"] if not r.get("P_matched", True))
                if len(rows) + n_un < 7:
                    partial = True
                sl.append(s)
            else:
                sl.append(float("nan"))
                partial = True
        ax.bar(x + (j + 1) * width, sl, width, color=colors.get(pn, "0.6"),
               label=f"ADAPT-grown, pool {pn} (P = 3n, matched rows only)")
        for xi, v, c in zip(x + (j + 1) * width, sl, cells):
            e = err(c, f"grown_{pn}")
            if e and np.isfinite(v):
                ax.errorbar([xi], [v], yerr=e, fmt="none", ecolor="black", capsize=3, lw=1)
    if boot:
        ax.text(0.99, 0.985, "error bars: 95% percentile bootstrap over the theta draws\n(2,000 resamples, s26/results/q_var_boot.json)",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.8, color="0.3")
    ax.axhline(-1.0, color="#d62728", ls="--", lw=1.2, label="2-design rate (S13 depth 8: base 0.504)")
    ax.axhline(0.0, color="0.7", lw=0.8)
    # the degenerate cells: every grown row stopped at P = n (the collapse), no slope exists
    for i, c in enumerate(cells):
        degenerate = all(
            not any(r.get("P_matched", True) for r in R["grown"][f"{pn}:adam_best:{c}"]["rows"])
            for pn in pools if f"{pn}:adam_best:{c}" in R["grown"])
        if degenerate and pools:
            ax.text(x[i] + width * (len(pools) + 1) / 2, -0.06, "grown:\ndegenerate\n(stops at P = n)",
                    ha="center", va="top", fontsize=6.8, color="0.3")
    ax.set_xticks(x + width * len(pools) / 2)
    ax.set_xticklabels([CELL_LABEL[c] for c in cells], rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("fitted log2 Var[dF/dtheta_0] per qubit, n = 4..13")
    ax.legend(fontsize=7.2, loc="lower center", ncol=2, frameon=True)
    ax.grid(True, axis="y", color="0.9")
    # right: the deployed cells, Var vs n
    ax = axes[1]
    ax.set_facecolor("white")
    for c, ls in (("deployed_a1_T03", "-"), ("deployed_a025_T03", "--")):
        if c not in R["fixed"]:
            continue
        rows = R["fixed"][c]["rows"]
        ax.plot([r["n"] for r in rows], [r["var_g0"] for r in rows], ls, color="0.35", marker="o",
                ms=4, label=f"fixed, {CELL_LABEL[c]}")
        for pn in pools:
            key = f"{pn}:adam_best:{c}"
            if key in R["grown"]:
                rr = [r for r in R["grown"][key]["rows"] if r.get("P_matched", True)]
                ax.plot([r["n"] for r in rr], [r["var_g0"] for r in rr], ls, color=colors.get(pn, "0.6"),
                        marker="s", ms=4, label=f"grown {pn}, {CELL_LABEL[c]}")
    ax.set_yscale("log")
    ax.set_xlabel("qubits n")
    ax.set_ylabel("Var[dF/dtheta_0] (exact parameter shift)")
    ax.legend(fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=True)
    ax.grid(True, which="major", color="0.9")
    title = "A4: gradient variance vs width, fixed vs ADAPT-grown, theta ~ N(0, 0.6^2), same draws as S25"
    if partial or not J.get("results", {}).get("slopes"):
        title += "  [PARTIAL FILE]"
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=DPI, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


if __name__ == "__main__":
    which = sys.argv[1:] or ["a2", "a4"]
    if "a2" in which:
        print("wrote", fig_a2())
    if "a4" in which and os.path.exists(os.path.join(RES, "q_var.json")):
        print("wrote", fig_a4())
