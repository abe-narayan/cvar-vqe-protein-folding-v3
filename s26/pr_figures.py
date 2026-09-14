"""PR lane, Sprint 26: figures for `vqe_research_overview.pptx`, regenerated from artefacts.

Every figure is saved at 190 dpi on a white background (`fig.savefig(path, dpi=190,
facecolor="white")`) so it sits on a white plate inside the dark deck. No number is typed by
hand: coordinates come from the production cache and the universe files, slopes from
`s25/results/q_plateau.json`, ladder values from `s26/pr_values.py`.

  pr_overlay_1S9Z.png   T030 built chain vs native CA trace, Kabsch-superposed (ORACLE use of the
  pr_overlay_9KAR.png   native, for the figure only), through `s12.instrument.superpose_batch`
  pr_width_sweep.png    S25 width sweep (five cells, fitted slopes) and the depth sweep at n = 7
  pr_ladder.png         the accuracy ladder with the basis named on every bar

Run: python s26/pr_figures.py   (writes to s26/figures/)
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
FIG = os.path.join(ROOT, "s26", "figures")
os.makedirs(FIG, exist_ok=True)

from s12 import instrument as I          # the instrument that produced every RMSD in the record
from s26 import pr_values

ACCENT = "#c8781e"      # the deck's amber, darkened for a white plate
NATIVE = "#2b2b33"
GREY = "#8a8a94"
BLUE = "#3b6fb6"
DPI = 190


def _save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, dpi=DPI, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- overlays
def overlay(pdb, tag):
    """Built chain (production cache `ca`) against the native CA trace (`nat_ca`)."""
    with open(os.path.join(I.ROOT, "bench_results", "cache", I.PROD_KEY, f"{pdb}.json")) as fh:
        rec = json.load(fh)
    ca = np.asarray(rec["ca"], float)
    z = np.load(os.path.join(I.UNIV, f"{pdb}.npz"), allow_pickle=True)
    nat = np.asarray(z["nat_ca"], float)                    # ORACLE: read for the figure only
    n = int(rec["n"])
    assert ca.shape == nat.shape == (n, 3)
    rmsd = I.ca_rmsd(ca, nat)
    assert abs(rmsd - rec["rmsd_arm"]) < 1e-9, (rmsd, rec["rmsd_arm"])
    moved = I.superpose_batch(ca[None], nat)[0]             # Kabsch, reflections forbidden
    # view: principal axes of the native, applied to both traces
    natc = nat - nat.mean(0)
    _, _, vt = np.linalg.svd(natc, full_matrices=False)
    P = natc @ vt.T
    Q = (moved - nat.mean(0)) @ vt.T

    views = [(0, 1), (0, 2)]
    both = np.vstack([P, Q])
    # the figure follows the trace's own aspect so a helix is not drawn inside a square of white
    spans = [(np.ptp(both[:, i]) + 2.0, np.ptp(both[:, j]) + 2.0) for i, j in views]
    ratio = max(sy / sx for sx, sy in spans)
    ratio = min(max(ratio, 0.28), 1.0)
    width = 6.8
    fig, axes = plt.subplots(1, 2, figsize=(width, width * 0.5 * ratio + 1.35))
    for ax, (i, j), lab in zip(axes, views, ["view 1 (principal axes 1, 2)", "view 2 (principal axes 1, 3)"]):
        ax.plot(P[:, i], P[:, j], "-o", color=NATIVE, lw=2.2, ms=4.5, label="native CA (model 1)")
        ax.plot(Q[:, i], Q[:, j], "-o", color=ACCENT, lw=2.2, ms=4.5, label="built chain (production)")
        ax.annotate("N", P[0, [i, j]], textcoords="offset points", xytext=(-9, -9), fontsize=8, color=NATIVE)
        ax.annotate("C", P[-1, [i, j]], textcoords="offset points", xytext=(5, 5), fontsize=8, color=NATIVE)
        lo, hi = both[:, [i, j]].min(0) - 1.0, both[:, [i, j]].max(0) + 1.0
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1])
        ax.set_aspect("equal")
        ax.set_title(lab, fontsize=9)
        ax.set_xlabel("A", fontsize=8); ax.set_ylabel("A", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f"{pdb} ({tag}), n = {n}: built chain vs native, CA-RMSD {rmsd:.3f} A "
                 f"(ORACLE-superposed; basis: built chain)", fontsize=9.5)
    fig.tight_layout(rect=(0, 0.09, 1, 0.97))
    return _save(fig, f"pr_overlay_{pdb}.png"), rmsd


# --------------------------------------------------------------------------- width sweep (S25)
def width_sweep():
    with open(os.path.join(ROOT, "s25", "results", "q_plateau.json")) as fh:
        r = json.load(fh)["results"]
    cells = [("linear_alpha1_T0", "alpha = 1, T = 0 (linear cost)", "#444444", "s"),
             ("cvar_alpha025_T0", "alpha = 0.25, T = 0", BLUE, "^"),
             ("cvar_alpha01_T0", "alpha = 0.10, T = 0", "#7aa6d9", "v"),
             ("deployed_a1_T03", "alpha = 1, T = 0.3 (deployed)", ACCENT, "o"),
             ("deployed_a025_T03", "alpha = 0.25, T = 0.3 (deployed)", "#e0a35c", "D")]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.5), gridspec_kw=dict(width_ratios=[1.55, 1]))
    for key, lab, col, mk in cells:
        rows = r[key]["rows"]
        n = np.array([x["n"] for x in rows], float)
        v = np.array([x["var_g0"] for x in rows], float)
        s = r[key]["log2_slope_per_qubit"]
        ax.plot(n, np.log2(v), marker=mk, color=col, lw=1.2, ms=5, label=f"{lab}: slope {s:+.3f}")
        a, b = np.polyfit(n, np.log2(v), 1)
        ax.plot(n, a * n + b, color=col, lw=0.8, ls="--", alpha=0.7)
    # a 2-design reference: the S13 depth-8 decay base (geo_kernel.json), drawn through the linear cell's first point
    with open(os.path.join(ROOT, "s13", "results", "geo_kernel.json")) as fh:
        base = json.load(fh)["scaling"]["lay8_v_w1"]["decay_base"]
    rows = r["linear_alpha1_T0"]["rows"]
    n0, v0 = rows[0]["n"], rows[0]["var_g0"]
    nn = np.array([rows[0]["n"], 10.0], float)
    ax.plot(nn, np.log2(v0) + np.log2(base) * (nn - n0), color="#b03030", lw=1.0, ls=":",
            label=f"reference: depth-8 decay base {base:.3f} per qubit (S13)")
    ax.set_xlabel("qubits n"); ax.set_ylabel("log2 Var[dF/dtheta_0]  (exact parameter shift)")
    ax.set_title("width sweep, depth 3, P = 3n, theta ~ N(0, 0.6^2)", fontsize=9.5)
    ax.legend(fontsize=6.6, loc="lower left", frameon=True)
    ax.grid(alpha=0.25); ax.tick_params(labelsize=8)
    ds = r["depth_sweep_n7"]
    L = [x["layers"] for x in ds]; v = [x["var_g0"] for x in ds]
    ax2.plot(L, v, "-o", color=ACCENT, lw=1.4, ms=5)
    ax2.set_yscale("log"); ax2.set_xlabel("layers L"); ax2.set_ylabel("Var[dF/dtheta_0]")
    ax2.set_title("depth sweep at n = 7, alpha = 0.25, T = 0.3", fontsize=9.5)
    ax2.axvline(3, color=GREY, ls="--", lw=0.9); ax2.annotate("deployed L = 3", (3, max(v)), fontsize=7.5, color=GREY,
                                                                 textcoords="offset points", xytext=(4, -2))
    ax2.grid(alpha=0.25, which="both"); ax2.tick_params(labelsize=8)
    fig.suptitle("Gradient variance of the deployed CVaR free energy (s25/results/q_plateau.json); "
                 "no shot noise, variance over theta only", fontsize=9)
    fig.tight_layout()
    return _save(fig, "pr_width_sweep.png")


# --------------------------------------------------------------------------- the accuracy ladder
def ladder(V):
    bars = [  # (label, token, tier)
        ("constant alpha-helix, zero information [built chain]", "HELIX", "control"),
        ("sequence-only torsion predictor [built chain]", "TORS", "control"),
        ("random 75-subset of the K=500 pool [point cloud]", "RANDOM75", "control"),
        ("shipped pipeline, built chain (the result)", "ARM_MEAN", "production"),
        ("shipped pipeline, point cloud (intermediate)", "AVG_MEAN", "production"),
        ("best member of the shipped top-75 [single window]", "TOPM_BEST", "oracle"),
        ("perfect distance prior, same pipeline [point cloud]", "PRIOR_PERFECT", "oracle"),
        ("best member of the K=500 pool [single window]", "POOL_BEST", "oracle"),
        ("torsion-space ceiling, k = 4 [torsion rebuild]", "TORSION_CEILING", "oracle"),
        ("best window in the whole library [single window]", "UNIVERSE_BEST", "oracle"),
        ("distance geometry from the true distances", "DISTGEO_TRUE", "oracle"),
    ]
    vals = [float(V[t]["value"]) for _, t, _ in bars]
    labs = [("ORACLE: " if tier == "oracle" else "") + lab for lab, _, tier in bars]
    cols = {"control": GREY, "production": ACCENT, "oracle": BLUE}
    fig, ax = plt.subplots(figsize=(7.4, 4.1))
    y = np.arange(len(bars))[::-1]
    ax.barh(y, vals, color=[cols[t] for _, _, t in bars], height=0.68)
    for yi, v in zip(y, vals):
        ax.text(v + 0.05, yi, f"{v:.3f}", va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=7.6)
    ax.set_xlabel("mean CA-RMSD over the 126 dev targets (A); basis in brackets", fontsize=8.5)
    ax.set_xlim(0, max(vals) + 0.9)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="x", alpha=0.25)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=ACCENT, label="PRODUCTION"), Patch(color=GREY, label="CONTROL (native-free)"),
                        Patch(color=BLUE, label="ORACLE (reads the native; a diagnostic, not a method)")],
               fontsize=7.5, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    ax.set_title("Where the accuracy lives: the ladder on one instrument", fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return _save(fig, "pr_ladder.png")


def main():
    V = pr_values.load_values()
    out = {}
    out["overlay_1S9Z"], r1 = overlay("1S9Z", "T030")
    out["overlay_9KAR"], r2 = overlay("9KAR", "hard target")
    out["width_sweep"] = width_sweep()
    out["ladder"] = ladder(V)
    for k, p in out.items():
        print(f"{k:14s} {os.path.relpath(p, ROOT)}")
    print(f"overlay RMSDs reproduce the record: 1S9Z {r1:.6f}  9KAR {r2:.6f}")
    return out


if __name__ == "__main__":
    main()
