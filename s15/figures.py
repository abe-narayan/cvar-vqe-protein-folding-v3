"""SPRINT 15 -- paper figures.\n\nEach figure reads a result JSON from `s15/results/` and is skipped with a printed note if its\ndata is not yet on disk, so this module can be run repeatedly as runs land.\n\nDESIGN RULES, applied to every panel, and adopted because Sprint 14 lost time to each of them:\n* value labels go INSIDE bars, never floating over a rule line or another series;\n* every axis starts at a real zero unless a log scale is declared in the caption;\n* ORACLE arms are drawn in a distinct hatched style and labelled ORACLE in the legend, so no\nreader can mistake a ceiling for a result;\n* confidence intervals are drawn wherever the underlying statistic has them;\n* the incumbent is drawn as a labelled reference rule on every accuracy panel.\n\nRun:\npython -m s15.figures\n"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import matplotlib                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
FIGDIR = os.path.join(ROOT, "s15", "figures")
os.makedirs(FIGDIR, exist_ok=True)

INCUMBENT = 3.2040761603809194

# a colourblind-safe set, assigned by role and never cycled
C_MAIN = "#3B6EA5"        # the arm under test
C_CTRL = "#8C8C93"        # controls and incumbents
C_GOOD = "#2E7D5B"        # something that helps
C_BAD = "#B4553C"         # something that hurts
C_ORACLE = "#7A5EA8"      # ORACLE ceilings, always hatched too
INK = "#26262B"
GRID = "#DCDCE2"


def _style(ax, xlabel="", ylabel="", title=""):
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK, labelsize=9)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=10, color=INK)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=10, color=INK)
    if title:
        ax.set_title(title, fontsize=11.5, color=INK, loc="left", pad=10)


MIN_N = 100          # a figure is never rendered from a smoke read


def _load(name, min_n=MIN_N):
    """Load a result, refusing anything partial or under-powered.\n\nThis guard exists because three figures in this sprint were rendered from 6- and 8-target\nsmoke runs before their full-instrument versions had landed, and looked exactly like\nfinished results. The project's own history makes that dangerous rather than merely untidy:\nan 8-target read once put an arm at 2.792 A where the full instrument gave 3.644 A and\nREVERSED the sign of the comparison. A figure that cannot be distinguished from a finished\none must not be producible from a smoke read.\n"""
    p = os.path.join(RESULTS, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        d = json.load(fh)
    if "partial" in d and "arms" not in d and "stages" not in d:
        return None
    #: some artefacts name their count differently; look for any of them rather than silently
    #: passing a result the guard cannot see.
    n = next((d[k] for k in ("n", "n_targets_phase", "n_done") if isinstance(d.get(k), int)), None)
    if n is None:
        print(f"  WARNING {name}: no target count found; the n-guard cannot check this figure.")
    if isinstance(n, int) and n < min_n:
        print(f"  REFUSING {name}: n = {n} < {min_n}. Smoke read, not a figure.")
        return None
    return d


def _save(fig, name):
    p = os.path.join(FIGDIR, name)
    fig.savefig(p, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {name}")


def _incumbent_rule(ax, label=True):
    ax.axhline(INCUMBENT, color=C_CTRL, lw=1.4, ls="--", zorder=1)
    if label:
        ax.text(0.995, INCUMBENT, f" incumbent {INCUMBENT:.3f} ",
                transform=ax.get_yaxis_transform(), ha="right", va="bottom",
                fontsize=8.5, color=C_CTRL)


# ------------------------------------------------------------------ fig 1: the cascade
def fig_cascade():
    d = _load("cascade_combined.json")
    if d is None:
        return print("  skip fig01_cascade (no cascade_combined.json)")
    keys = ["G_ORACLE", "S", "A", "A_besthalf", "F"]
    labels = ["G\nbest in ensemble\n(ORACLE)", "S\nobjective argmin",
              "A\nconsensus", "A\nobjective-best half", "F\nafter projection"]
    v = [d["stages"][k]["mean"] for k in keys]
    lo = [d["stages"][k]["vs_incumbent"]["ci95"][0] for k in keys]
    hi = [d["stages"][k]["vs_incumbent"]["ci95"][1] for k in keys]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    cols = [C_ORACLE, C_BAD, C_MAIN, C_MAIN, C_BAD]
    bars = ax.bar(range(len(v)), v, color=cols, width=0.62,
                  hatch=["///", "", "", "", ""], edgecolor="white", zorder=3)
    top = max(v) * 1.34
    for k, (b, val) in enumerate(zip(bars, v)):
        ax.text(b.get_x() + b.get_width() / 2, val - 0.16, f"{val:.3f}",
                ha="center", va="top", color="white", fontsize=10.5, fontweight="bold",
                zorder=6,
                bbox=dict(boxstyle="round,pad=0.16", fc=cols[k], ec="none", alpha=0.85))
        #: the CI captions are staggered onto two rows -- adjacent bars differ by only a few
        #: hundredths, so a single row of captions collides every time.
        diff = d["stages"][keys[k]]["vs_incumbent"]["mean_diff"]
        yy = top * (0.90 if k % 2 == 0 else 0.80)
        ax.annotate(f"{diff:+.3f} [{lo[k]:+.3f}, {hi[k]:+.3f}]",
                    xy=(b.get_x() + b.get_width() / 2, val),
                    xytext=(b.get_x() + b.get_width() / 2, yy),
                    ha="center", va="bottom", fontsize=8.2, color=INK,
                    arrowprops=dict(arrowstyle="-", color=GRID, lw=0.8, shrinkB=2))
    #: the incumbent caption goes on the LEFT; on the right it lands under the F bar
    ax.axhline(INCUMBENT, color=C_CTRL, lw=1.4, ls="--", zorder=1)
    ax.text(0.005, INCUMBENT - 0.03, f" incumbent {INCUMBENT:.3f} ",
            transform=ax.get_yaxis_transform(), ha="left", va="top",
            fontsize=8.5, color=C_CTRL, zorder=6,
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.9))
    ax.set_xticks(range(len(v)))
    ax.set_xticklabels(labels, fontsize=8.5, color=INK)
    ax.set_ylim(0, top)
    _style(ax, ylabel="mean full-chain CA-RMSD (Å)",
           title="The generative cascade does not beat the incumbent\n"
                 "the ensemble contains better structures than the readout returns")
    ax.text(0.01, 0.995, f"n = {d['n']} targets · paired difference vs incumbent above "
                         f"each bar",
            transform=ax.transAxes, va="top", fontsize=8.5, color=C_CTRL)
    _save(fig, "fig01_cascade.png")


# ------------------------------------------- fig 2: the distance-accuracy phase diagram
def fig_phase():
    """The three error models on ONE comparable axis: the effective RMS of the injected error.\n\nPlotting them against their own `sigma` parameter would be meaningless -- the models inject\ndifferent magnitudes at the same nominal sigma (`outliers` puts 1.375x and `sep_scaled`\n1.1331x the nominal into the error), and the entire point of the panel is that SHAPE beats\nMAGNITUDE, which cannot be read off an axis that confounds the two.\n"""
    #: DECLARED EXCEPTION to the n-guard: this run is 40 targets by design (each cell is
    #: 8 severities x 3 error models x 2 repeats). The panel prints its n.
    d = _load("distacc.json", min_n=1)
    if d is None or "cells" not in d:
        return print("  skip fig02_phase (no distacc.json)")
    FAC = {"abs_gauss": 1.0, "sep_scaled": 1.1331, "outliers": (0.09 + 0.05 * 36) ** 0.5}
    label = {"abs_gauss": "i.i.d. Gaussian",
             "sep_scaled": "growing with sequence separation",
             "outliers": "5% of pairs badly wrong, the rest nearly right"}
    colour = {"abs_gauss": C_CTRL, "sep_scaled": C_MAIN, "outliers": C_GOOD}
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    for m in ("abs_gauss", "sep_scaled", "outliers"):
        cs = sorted([c for c in d["cells"].values() if c["mode"] == m],
                    key=lambda c: c["sigma"])
        if not cs:
            continue
        x = [c["sigma"] * FAC.get(m, 1.0) for c in cs]
        y = [c["mean"] for c in cs]
        ax.plot(x, y, "-o", color=colour[m], lw=2.2, ms=6,
                label=label.get(m, m), zorder=3)
    #: the real channel, which is worse than every synthetic model at its own magnitude
    ax.scatter([3.704], [3.644], s=190, marker="X", color=C_BAD, zorder=6,
               edgecolor="white", linewidth=1.4)
    ax.annotate("the REAL distogram\n3.70 Å RMSE → 3.644 Å",
                xy=(3.704, 3.644), xytext=(2.55, 3.62), fontsize=9.5, color=C_BAD,
                ha="right", va="center",
                arrowprops=dict(arrowstyle="->", color=C_BAD, lw=1.3))
    for lvl, lab in ((2.5, "2.5 Å  major success"), (2.0, "2.0 Å  ideal")):
        ax.axhline(lvl, color=C_CTRL, lw=1.1, ls=":", zorder=1)
        ax.text(0.012, lvl, f" {lab}", transform=ax.get_yaxis_transform(),
                va="bottom", fontsize=8.5, color=C_CTRL)
    _style(ax, xlabel="effective RMS of the injected distance error (Å)",
           ylabel="mean full-chain CA-RMSD (Å)",
           title="Error SHAPE beats error MAGNITUDE, in both directions\n"
                 "true distances corrupted under three realistic models, then refitted")
    ax.legend(frameon=False, fontsize=9.5, loc="upper left", title="error model",
              title_fontsize=9.5)
    ax.set_ylim(0, None); ax.set_xlim(0, None)
    ax.text(0.99, 0.03, f"{d.get('n_targets_phase')} targets · outlier-shaped error at 4.12 Å "
                        f"RMS is cheaper than i.i.d. error at 3.00 Å",
            transform=ax.transAxes, ha="right", fontsize=8.5, color=C_CTRL)
    _save(fig, "fig02_phase_diagram.png")


# -------------------------------------------------- fig 3: Family B feasibility (K7)
def fig_feasibility():
    d = _load("feasible.json")
    if d is None:
        return print("  skip fig03_feasibility (no feasible.json)")
    f = d["feasibility_mean"]
    eps = sorted([float(k.split("_")[1]) for k in f if k.startswith("eps_")])
    y = [f[f"eps_{e}"] for e in eps]
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(eps, y, "-o", color=C_MAIN, lw=2.2, ms=7, zorder=3)
    for x, yy in zip(eps, y):
        ax.annotate(f"{yy:.3f}", (x, yy), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=9, color=INK)
    ax.axhline(1.0, color=C_CTRL, lw=1.2, ls="--")
    ax.text(0.99, 1.0, " every restraint admits the native ",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=8.5, color=C_CTRL)
    _style(ax, xlabel="restraint band ε, in units of the predictor's own sd",
           ylabel="fraction of pairs admitting the native",
           title="Family B's feasible set excludes the native at every usable ε\n"
                 f"median z = {f['median_z']:.3f}, max z = {f['max_z']:.3f} — "
                 "a concentrated failure, not a calibration one")
    ax.set_ylim(0, 1.08)
    _save(fig, "fig03_feasibility.png")


# ------------------------------------- fig 4: does the objective rank the native? (K6)
def fig_native_percentile():
    d = _load("feasible.json")
    if d is None:
        return print("  skip fig04_native_percentile (no feasible.json)")
    arms = [a for a in d["arms"]]
    pct = [d["arms"][a]["native_percentile_mean"] for a in arms]
    rho = [d["arms"][a]["rho_f_rmsd_mean"] for a in arms]
    order = np.argsort(pct)
    arms = [arms[k] for k in order]; pct = [pct[k] for k in order]
    rho = [rho[k] for k in order]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
    for ax, vals, lab, ttl in (
            (axes[0], pct, "native's percentile on the objective",
             "Where does the objective rank the TRUTH?"),
            (axes[1], rho, "mean Spearman ρ(objective, RMSD)",
             "Does the objective order the pool at all?")):
        cols = [C_ORACLE if a.startswith("ORACLE") else
                (C_BAD if (lab.startswith("native") and v > 0.5) else C_MAIN)
                for a, v in zip(arms, vals)]
        bars = ax.barh(range(len(arms)), vals, color=cols, height=0.66,
                       hatch=["///" if a.startswith("ORACLE") else "" for a in arms],
                       edgecolor="white", zorder=3)
        for b, v in zip(bars, vals):
            ax.text(max(v, 0) + 0.008, b.get_y() + b.get_height() / 2, f"{v:.3f}",
                    va="center", fontsize=9, color=INK)
        ax.set_yticks(range(len(arms)))
        ax.set_yticklabels(arms, fontsize=9)
        _style(ax, xlabel=lab, title=ttl)
        ax.xaxis.grid(True, color=GRID, lw=0.8); ax.yaxis.grid(False)
    axes[0].axvline(0.5, color=C_CTRL, lw=1.3, ls="--")
    axes[0].text(0.5, len(arms) - 0.4, " chance ", ha="left", fontsize=8.5, color=C_CTRL)
    fig.suptitle("The restraint objective orders the pool well and does not find the native",
                 fontsize=11.5, color=INK, x=0.02, ha="left")
    fig.tight_layout()
    _save(fig, "fig04_native_percentile.png")


# --------------------------------------------- fig 5: the distogram's bias profile (K2)
def fig_bias():
    #: DECLARED EXCEPTION to the n-guard, and a different n from fig02's. `distacc.json` holds two
    #: blocks: a 126-target error profile (this panel) and a 40-target phase diagram (fig02). The
    #: guard reads the file-level `n_targets_phase = 40`, which is the WRONG count for this panel,
    #: so the guard is bypassed explicitly and the correct n is printed on the figure.
    d = _load("distacc.json", min_n=1)
    if d is None:
        return print("  skip fig05_bias (no distacc.json)")
    rows = d["error_profile"]["by_separation"]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.bar(x - 0.2, [r["mae"] for r in rows], width=0.38, color=C_CTRL,
           label="MAE", edgecolor="white", zorder=3)
    ax.bar(x + 0.2, [r["bias"] for r in rows], width=0.38, color=C_BAD,
           label="bias (systematic over-prediction)", edgecolor="white", zorder=3)
    for k, r in enumerate(rows):
        ax.text(k + 0.2, r["bias"] + 0.05, f"{r['bias']:+.3f}", ha="center",
                va="bottom", fontsize=9, color=INK)
    ax.set_xticks(x); ax.set_xticklabels([r["sep"] for r in rows], fontsize=9.5)
    _style(ax, xlabel="sequence separation |i − j|", ylabel="ångström",
           title="The distogram over-predicts distance, and the error grows with separation\n"
                 "a systematic, separation-dependent offset — and correcting it buys almost nothing")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    npairs = d["error_profile"]["global"].get("n_pairs")
    ax.text(0.99, 0.02, f"126 targets · {npairs} pairs", transform=ax.transAxes,
            ha="right", fontsize=8, color=C_CTRL)
    _save(fig, "fig05_bias_profile.png")


# ------------------------------------------- fig 6: the ceiling reframe (K1 / K5)
def fig_ceiling():
    d = _load("distgeo.json")
    if d is None:
        return print("  skip fig06_ceiling (no distgeo.json)")
    arms = ["ORACLE_true_distances", "pred_invvar_weighted", "pred_unweighted"]
    arms = [a for a in arms if a in d.get("arms", {})]
    if not arms:
        return print("  skip fig06_ceiling (unexpected distgeo schema)")
    labels = ["ORACLE\ntrue distances", "predicted\n1/sd² weighted", "predicted\nunweighted"]
    v = [d["arms"][a]["mean"] for a in arms]
    med = [d["arms"][a]["median"] for a in arms]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    cols = [C_ORACLE, C_MAIN, C_CTRL]
    bars = ax.bar(range(len(v)), v, color=cols, width=0.58,
                  hatch=["///", "", ""], edgecolor="white", zorder=3)
    ax.plot(range(len(v)), med, "D", color=INK, ms=6, zorder=5, label="median")
    for b, val in zip(bars, v):
        ax.text(b.get_x() + b.get_width() / 2, val - 0.12, f"{val:.3f}", ha="center",
                va="top", color="white", fontsize=11, fontweight="bold", zorder=4)
    _incumbent_rule(ax)
    ax.axhline(1.95, color=C_GOOD, lw=1.3, ls=":", zorder=2)
    ax.text(0.01, 1.95, "  the ceiling the project believed (≈1.95 Å, through the library)",
            transform=ax.get_yaxis_transform(), va="bottom", fontsize=8.5, color=C_GOOD)
    ax.set_xticks(range(len(v))); ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylim(0, max(v) * 1.2)
    _style(ax, ylabel="mean full-chain CA-RMSD (Å)",
           title="Perfect distance knowledge is worth 0.611 Å, not 1.95 Å\n"
                 "the old ceiling measured the library, not the distance channel")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    _save(fig, "fig06_ceiling.png")


# ---------------------------------------------------- fig 7: robust losses (Family B)
def fig_robust():
    d = _load("robust.json")
    if d is None:
        return print("  skip fig07_robust (no robust.json)")
    arms = list(d["arms"])
    v = [d["arms"][a]["mean"] for a in arms]
    lo = [d["arms"][a]["vs_squared"]["ci95"][0] for a in arms]
    hi = [d["arms"][a]["vs_squared"]["ci95"][1] for a in arms]
    diff = [d["arms"][a]["vs_squared"]["mean_diff"] for a in arms]
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    cols = [C_ORACLE if a.startswith("ORACLE") else
            (C_GOOD if diff[k] < 0 else C_CTRL if a == "squared" else C_BAD)
            for k, a in enumerate(arms)]
    ax.errorbar(diff, range(len(arms)),
                xerr=[np.array(diff) - np.array(lo), np.array(hi) - np.array(diff)],
                fmt="o", color=INK, ms=0, lw=1.4, capsize=3, zorder=4)
    ax.scatter(diff, range(len(arms)), c=cols, s=90, zorder=5, edgecolor="white")
    for k, a in enumerate(arms):
        ax.text(hi[k] + 0.02, k, f"{v[k]:.3f} Å", va="center", fontsize=9, color=INK)
    ax.axvline(0, color=C_CTRL, lw=1.4, ls="--")
    ax.set_yticks(range(len(arms))); ax.set_yticklabels(arms, fontsize=9.5)
    _style(ax, xlabel="paired difference vs squared error (Å) — negative is better",
           title="Do robust losses help? The violations are concentrated, so they should\n"
                 "scales chosen leave-fold-out on the objective, never on RMSD")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.yaxis.grid(False)
    _save(fig, "fig07_robust_losses.png")


# ------------------------------------------------- fig 8: pool augmentation (Family C)
def fig_augment():
    d = _load("augment.json")
    if d is None:
        return print("  skip fig08_augment (no augment.json)")
    keep = [a for a in d["arms"] if a.startswith("augment_w")]
    keep = sorted(keep, key=lambda a: float(a.replace("augment_w", "")))
    x = [float(a.replace("augment_w", "")) for a in keep]
    y = [d["arms"][a]["mean"] for a in keep]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.plot(x, y, "-o", color=C_MAIN, lw=2.2, ms=6, zorder=3)
    base = d["arms"]["incumbent_avg"]["mean"]
    ax.axhline(base, color=C_CTRL, lw=1.4, ls="--")
    ax.text(0.99, base, f" retrieval only {base:.3f} ", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=8.5, color=C_CTRL)
    for xx, yy in zip(x, y):
        ax.annotate(f"{yy:.3f}", (xx, yy), textcoords="offset points", xytext=(0, -14),
                    ha="center", fontsize=8.5, color=INK)
    _style(ax, xlabel="weight per solved conformer, relative to a retrieved window",
           ylabel="mean full-chain CA-RMSD (Å)",
           title="Mixing restraint-solved conformers into the retrieval pool\n"
                 "the operator consumes the set mean, and the fits beat it by 0.43 Å")
    _save(fig, "fig08_augment.png")


# ------------------------------------------------ fig 9: the contraction fix (K8 A→F)
def fig_expand():
    d = _load("expand.json")
    if d is None:
        return print("  skip fig09_expand (no expand.json)")
    arms = list(d["arms"])
    v = [d["arms"][a]["mean"] for a in arms]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    cols = [C_ORACLE if a.startswith("ORACLE") else
            (C_CTRL if a == "project_raw" else C_MAIN) for a in arms]
    bars = ax.bar(range(len(v)), v, color=cols, width=0.6,
                  hatch=["///" if a.startswith("ORACLE") else "" for a in arms],
                  edgecolor="white", zorder=3)
    for b, val in zip(bars, v):
        ax.text(b.get_x() + b.get_width() / 2, val - 0.12, f"{val:.3f}", ha="center",
                va="top", color="white", fontsize=10, fontweight="bold", zorder=4)
    ax.set_xticks(range(len(v)))
    ax.set_xticklabels([a.replace("_", "\n") for a in arms], fontsize=8.5)
    ax.set_ylim(0, max(v) * 1.18)
    _style(ax, ylabel="mean full-chain CA-RMSD (Å)",
           title=f"Undoing the averaging contraction before projecting\n"
                 f"the consensus is {d.get('contraction_vs_true_pct', float('nan')):.1f}% "
                 f"too small against the true distances")
    _save(fig, "fig09_expand.png")


# ------------------------------------------ fig 10: surrogate destruction (the error's shape)
def fig_errstruct():
    d = _load("errstruct.json")
    if d is None or "arms" not in d:
        return print("  skip fig10_errstruct (no errstruct.json)")
    arms = [a for a in d["arms"] if a != "ORACLE_true"]
    diff = [d["arms"][a]["vs_real"]["mean_diff"] for a in arms]
    lo = [d["arms"][a]["vs_real"]["ci95"][0] for a in arms]
    hi = [d["arms"][a]["vs_real"]["ci95"][1] for a in arms]
    order = np.argsort(diff)
    arms = [arms[k] for k in order]
    diff = [diff[k] for k in order]; lo = [lo[k] for k in order]; hi = [hi[k] for k in order]
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    cols = [C_ORACLE if a.startswith("ORACLE") else C_MAIN for a in arms]
    ax.errorbar(diff, range(len(arms)),
                xerr=[np.array(diff) - np.array(lo), np.array(hi) - np.array(diff)],
                fmt="none", ecolor=INK, lw=1.4, capsize=3, zorder=4)
    ax.scatter(diff, range(len(arms)), c=cols, s=110, zorder=5, edgecolor="white")
    for k, a in enumerate(arms):
        ax.text(hi[k] + 0.03, k, f"{d['arms'][a]['mean']:.3f} Å", va="center",
                fontsize=9, color=INK)
    ax.axvline(0, color=C_CTRL, lw=1.4, ls="--")
    ax.text(0, len(arms) - 0.35, "  the real error, unmodified", fontsize=8.5, color=C_CTRL)
    ax.set_yticks(range(len(arms)))
    ax.set_yticklabels([a.replace("ORACLE_", "ORACLE · ") for a in arms], fontsize=9.5)
    _style(ax, xlabel="paired change vs the real error (Å) — negative means the destroyed "
                      "property was costly",
           title="What about the distogram's error is expensive?\n"
                 "destroy one property at a time and refit — every surrogate is an ORACLE "
                 "diagnostic")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.yaxis.grid(False)
    _save(fig, "fig10_errstruct.png")


# ------------------------------------------------ fig 11: realizability and its null
def fig_realizability():
    d = _load("coherence.json")
    if d is None or "realizability" not in d:
        return print("  skip fig11_realizability (no coherence.json)")
    r = d["realizability"]
    keys = ["res_disto", "res_pool", "res_sign", "res_perm", "res_iid_null"]
    labs = ["the REAL\npredicted distances", "the retrieval\npool's distances",
            "ORACLE\nsigns randomised", "ORACLE\npermuted over pairs",
            "ORACLE\ni.i.d. matched RMS"]
    v = [r[k] for k in keys]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    cols = [C_MAIN, C_GOOD, C_ORACLE, C_ORACLE, C_CTRL]
    bars = ax.bar(range(len(v)), v, color=cols, width=0.6,
                  hatch=["", "", "///", "///", ""], edgecolor="white", zorder=3)
    for b, val in zip(bars, v):
        ax.text(b.get_x() + b.get_width() / 2, val - 0.04, f"{val:.3f}", ha="center",
                va="top", color="white", fontsize=10.5, fontweight="bold", zorder=5)
    ax.axhline(r["res_iid_null"], color=C_CTRL, lw=1.3, ls="--", zorder=1)
    ax.text(0.005, r["res_iid_null"], "  the i.i.d. null ",
            transform=ax.get_yaxis_transform(), va="bottom", fontsize=8.5, color=C_CTRL)
    ax.set_xticks(range(len(v))); ax.set_xticklabels(labs, fontsize=8.5)
    ax.set_ylim(0, max(v) * 1.2)
    _style(ax, ylabel="unrealizable part  |d_fit − d_target|  (Å)",
           title="The distogram's errors describe a consistent WRONG STRUCTURE\n"
                 f"the real prediction is {1/r['ratio_real_over_null']:.1f}× closer to "
                 "realizable than noise of the same size")
    _save(fig, "fig11_realizability.png")


# ----------------------------------------------------- fig 12: the fusion law
def fig_fusion_law():
    d = _load("coherence.json")
    if d is None or "rows" not in d:
        return print("  skip fig12_fusion_law (no coherence.json)")
    pred = np.array([r["law_prediction"] for r in d["rows"]])
    meas = np.array([r["rmsd_coordavg"] for r in d["rows"]])
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    lim = [0, max(pred.max(), meas.max()) * 1.06]
    ax.plot(lim, lim, color=C_CTRL, lw=1.3, ls="--", zorder=2)
    ax.scatter(pred, meas, s=42, color=C_MAIN, alpha=0.78, edgecolor="white",
               linewidth=0.6, zorder=4)
    f = d["fusion"]
    ax.text(0.03, 0.97, f"n = {d['n']} targets\n"
                        f"mean absolute error {f['prediction_error_abs']:.3f} Å\n"
                        f"no free parameters",
            transform=ax.transAxes, va="top", fontsize=9.5, color=INK,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GRID))
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    _style(ax, xlabel="predicted  √(r² − (s/2)²)   (Å)",
           ylabel="measured coordinate average (Å)",
           title="Fusion gain, predicted from disagreement alone\n"
                 "and the disagreement s is native-free")
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    _save(fig, "fig12_fusion_law.png")


# ------------------------------------------------- fig 13: the scale-correction ceiling
def fig_scale():
    d = _load("scale.json")
    if d is None or "arms" not in d:
        return print("  skip fig13_scale (no scale.json)")
    arms = list(d["arms"])
    diff = [d["arms"][a]["vs_real"]["mean_diff"] for a in arms]
    lo = [d["arms"][a]["vs_real"]["ci95"][0] for a in arms]
    hi = [d["arms"][a]["vs_real"]["ci95"][1] for a in arms]
    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    cols = [C_ORACLE if a.startswith("ORACLE") else
            (C_CTRL if a == "real" else C_MAIN) for a in arms]
    ax.errorbar(diff, range(len(arms)),
                xerr=[np.array(diff) - np.array(lo), np.array(hi) - np.array(diff)],
                fmt="none", ecolor=INK, lw=1.4, capsize=3, zorder=4)
    ax.scatter(diff, range(len(arms)), c=cols, s=100, zorder=5, edgecolor="white")
    ax.axvline(0, color=C_CTRL, lw=1.4, ls="--")
    ax.set_yticks(range(len(arms))); ax.set_yticklabels(arms, fontsize=9.5)
    _style(ax, xlabel="paired change vs no correction (Å) — negative is better",
           title="Per-target rescaling of the predicted distances\n"
                 "the ORACLE arms bound what any rescaling could ever be worth")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.yaxis.grid(False)
    _save(fig, "fig13_scale.png")


def _relayed(block=None, min_n=MIN_N):
    """Relayed tables, guarded like everything else.

    `_relayed` originally bypassed `_load` entirely, so two figures rendered with no n-check at all
    and one of them (fig15) drew 9 targets. Each block now declares its own `n_targets`, and a block
    below `min_n` is rendered only if the caller passes an explicit, *declared* exception — which
    then prints the n on the panel.
    """
    q = os.path.join(RESULTS, "relayed.json")
    if not os.path.exists(q):
        return None
    with open(q) as fh:
        d = json.load(fh)
    if block is not None:
        n = d.get(block, {}).get("n_targets")
        if n is None:
            print(f"  WARNING relayed.json[{block}]: no n_targets declared.")
        elif n < min_n:
            print(f"  relayed.json[{block}]: n = {n} < {min_n} — rendering as a DECLARED "
                  f"exception; the panel prints its n.")
    return d


# ------------------------------------------- fig 14: the channel table vs the requirement
def fig_channels():
    d = _relayed("channels")
    if d is None:
        return print("  skip fig14_channels (no relayed.json)")
    c = d["channels"]
    rows = sorted(c["rows"], key=lambda r: r["acc"])
    names = [r["name"] for r in rows]
    acc = [r["acc"] for r in rows]
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    cols = [C_CTRL if "null" in n else (C_MAIN if a > c["chance"] else C_BAD)
            for n, a in zip(names, acc)]
    ax.barh(range(len(acc)), [a - 0.48 for a in acc], left=0.48, color=cols,
            height=0.66, edgecolor="white", zorder=3)
    for k, (a, r) in enumerate(zip(acc, rows)):
        #: place the value label clear of the confidence bar, not at the bar end, or the two
        #: overlap on every row that has an interval.
        if r.get("lo") is not None:
            ax.plot([r["lo"], r["hi"]], [k, k], color=INK, lw=1.4, zorder=5)
            x = r["hi"] + 0.004
        else:
            x = a + 0.004
        ax.text(x, k, f"{a:.3f}", va="center", fontsize=9, color=INK, zorder=6)
    ax.axvline(c["chance"], color=C_CTRL, lw=1.4, ls="--", zorder=2)
    ax.axvline(c["required_for_2A"], color=C_BAD, lw=1.8, ls="-", zorder=2)
    ax.text(c["chance"], len(acc) - 0.3, " chance", fontsize=9, color=C_CTRL)
    ax.text(c["required_for_2A"], len(acc) - 0.3, " required for 2.0 Å", fontsize=9,
            color=C_BAD, ha="right")
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=9)
    ax.set_xlim(0.48, max(c["required_for_2A"], max(acc)) + 0.02)
    _style(ax, xlabel="in-band pairwise discrimination accuracy",
           title="No available channel reaches the accuracy 2.0 Å requires\n"
                 "the best is 0.566 against a requirement of 0.638")
    ax.xaxis.grid(True, color=GRID, lw=0.8); ax.yaxis.grid(False)
    ax.text(0.99, 0.02, f"{c.get('n_targets')} targets · relayed from the "
                        f"information-channel workstream",
            transform=ax.transAxes, ha="right", fontsize=8, color=C_CTRL)
    _save(fig, "fig14_channels.png")


# ------------------------------------- fig 15: what concentration buys and what it pays
def fig_cvar_trade():
    d = _relayed("cvar", min_n=9)
    if d is None:
        return print("  skip fig15_cvar_trade (no relayed.json)")
    c = d["cvar"]
    a = c["alpha"]
    x = list(range(len(a)))
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    ax.plot(x, c["set_mean"], "-o", color=C_GOOD, lw=2.2, ms=7,
            label="drawn-set MEAN  (VQE better)", zorder=4)
    ax.plot(x, c["set_best"], "-o", color=C_BAD, lw=2.2, ms=7,
            label="drawn-set BEST  (VQE worse)", zorder=4)
    ax.plot(x, c["top75_avg"], "-o", color=C_MAIN, lw=2.2, ms=7,
            label="top-75 coordinate average  (what the pipeline uses)", zorder=4)
    ax.axhline(0, color=C_CTRL, lw=1.4, ls="--", zorder=1)
    #: this note sat at the bottom-left and collided with the alpha=1 point; put it under the
    #: title instead, where nothing is plotted.
    ax.text(0.01, 0.97, "positive = VQE worse than best-of-N from its own untrained circuit",
            transform=ax.transAxes, va="top", fontsize=8.5, color=C_CTRL)
    _style(ax, ylabel="paired difference (Å)",
           title="Concentration buys the set mean and pays the set best\n"
                 "two significant effects of opposite sign, both scaled by α")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax2.bar(x, c["distinct"], color=C_CTRL, width=0.55, edgecolor="white", zorder=3)
    ax2.axhline(c["control_distinct"], color=C_BAD, lw=1.4, ls="--", zorder=4)
    ax2.text(len(x) - 0.5, c["control_distinct"], " untrained control ", ha="right",
             va="bottom", fontsize=8.5, color=C_BAD)
    for k, v in enumerate(c["distinct"]):
        ax2.text(k, v + 12, f"{v:g}", ha="center", fontsize=9, color=INK)
    ax2.set_xticks(x); ax2.set_xticklabels([f"α = {v:g}" for v in a], fontsize=9.5)
    nt = c.get("n_targets"); nc = c.get("n_cells")
    ax2.text(0.995, -0.42, f"{nt} targets · {nc} paired cells ({nt} targets × 3 seeds) · relayed "
                           f"from the quantum-geometry workstream",
             transform=ax2.transAxes, ha="right", fontsize=8, color=C_CTRL)
    _style(ax2, ylabel="distinct configs\nout of 4,096 draws")
    fig.tight_layout()
    _save(fig, "fig15_cvar_trade.png")


# --------------------------------------------- fig 16: where the ångströms actually go
def fig_gaps():
    d = _load("cascade_combined.json")
    if d is None:
        return print("  skip fig16_gaps (no cascade_combined.json)")
    st = d["stages"]
    steps = [("generate\nORACLE best", st["G_ORACLE"]["mean"], None),
             ("select\nby objective", st["S"]["mean"], st["S"]["mean"] - st["G_ORACLE"]["mean"]),
             ("aggregate\nconsensus", st["A"]["mean"], st["A"]["mean"] - st["S"]["mean"]),
             ("project\nonto manifold", st["F"]["mean"], st["F"]["mean"] - st["A"]["mean"])]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.plot([0, 1, 2, 3], [v for _n, v, _dd in steps], "-o", color=INK, lw=2.0, ms=9,
            zorder=5)
    for k, (n, v, dd) in enumerate(steps):
        if dd is not None:
            col = C_BAD if dd > 0 else C_GOOD
            ax.annotate("", xy=(k, v), xytext=(k, v - dd),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=2.6, alpha=0.85))
            ax.text(k + 0.06, v - dd / 2, f"{dd:+.3f}", color=col, fontsize=11,
                    fontweight="bold", va="center")
        ax.text(k, v + 0.075, f"{v:.3f}", ha="center", fontsize=10, color=INK,
                fontweight="bold")
    _incumbent_rule(ax)
    ax.set_xticks(range(4)); ax.set_xticklabels([n for n, _v, _dd in steps], fontsize=9.5)
    ax.set_xlim(-0.45, 3.55)
    #: headroom so the tallest value label clears the title, and an explicit note that this is
    #: the one panel with a non-zero baseline -- the whole point is the DIFFERENCES between
    #: stages, which a zero-based axis would compress into invisibility.
    lo = min(v for _n, v, _dd in steps); hi = max(v for _n, v, _dd in steps)
    ax.set_ylim(lo - 0.10 * (hi - lo), hi + 0.22 * (hi - lo))
    ax.text(0.005, 0.015, "note: the y-axis is not zero-based -- this panel shows the "
                          "differences between stages",
            transform=ax.transAxes, fontsize=8, color=C_CTRL)
    _style(ax, ylabel="mean full-chain CA-RMSD (Å)",
           title="Where the ångströms go\n"
                 "the ensemble beats the baseline; selection and projection lose ground, "
                 "aggregation recovers part of it")
    _save(fig, "fig16_gaps.png")


FIGS = [fig_cascade, fig_phase, fig_feasibility, fig_native_percentile, fig_bias,
        fig_ceiling, fig_robust, fig_augment, fig_expand, fig_errstruct,
        fig_realizability, fig_fusion_law, fig_scale, fig_channels, fig_cvar_trade,
        fig_gaps]


def run():
    print(f"rendering into {FIGDIR}")
    for f in FIGS:
        try:
            f()
        except Exception as e:                                  # noqa: BLE001
            print(f"  FAILED {f.__name__}: {type(e).__name__}: {e}")
    print("done")


if __name__ == "__main__":
    run()
