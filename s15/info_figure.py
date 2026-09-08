"""SPRINT 15, INFO -- the paper figure: coverage x uncertainty -> predictive RMSD.

Panel A  the primary phase surface with the 3.0 / 2.5 / 2.0 A contours and the regime the
         project actually occupies marked on it.
Panel B  the sigma = 0 slice: WHERE the gaps fall, at identical coverage.
Panel C  the sigma each error model needs at full coverage to reach each level.
Panel D  alignment of each real channel's error with the RMSD-quiet subspace -- the reason
         panel A over-prices every fragment channel.

    python -m s15.info_figure
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s15 import info_lib as L            # noqa: E402

C_LINE = "#111111"
PAL = ["#2166AC", "#B2182B", "#1B7837", "#D6800B", "#762A83"]


def main():
    ph = json.load(open(os.path.join(L.RESULTS, "info_phase.json")))
    rg = json.load(open(os.path.join(L.RESULTS, "info_regime.json")))
    nl = json.load(open(os.path.join(L.RESULTS, "info_null.json")))
    sig = np.asarray(ph["sigmas"], float)
    cov = np.asarray(ph["covs"], float)
    G = np.asarray(ph["surfaces"]["primary_iid_uniform_pool"]["grid"])

    fig = plt.figure(figsize=(13.0, 9.6))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.26)

    # ---------------- Panel A
    ax = fig.add_subplot(gs[0, 0])
    X, Y = np.meshgrid(cov, sig)
    cf = ax.contourf(X, Y, G, levels=np.linspace(0.3, 5.2, 25), cmap="viridis_r")
    cs = ax.contour(X, Y, G, levels=[2.0, 2.5, 3.0], colors=["#FFFFFF", "#FFD166",
                                                             "#EF476F"], linewidths=2.2)
    ax.clabel(cs, fmt="%.1f A", fontsize=9)
    plt.colorbar(cf, ax=ax, label="CA-RMSD (A)")
    ax.set_xlabel("coverage: fraction of residues with a target-specific torsion channel")
    ax.set_ylabel("per-torsion sigma (deg), i.i.d.")
    ax.set_title("A. Phase surface -- 126 targets, uncovered residues filled from the\n"
                 "top-75 retrieval pool (i.i.d. errors: an UPPER BOUND for fragment "
                 "channels)", fontsize=10, loc="left")
    inc = rg["incumbent_sigma_equivalent_deg"]
    ax.plot([1.0], [inc], "o", ms=9, mfc="none", mec="w", mew=2)
    ax.annotate(f"incumbent 3.204 A\nsigma-equiv {inc:.0f} deg", (1.0, inc),
                textcoords="offset points", xytext=(-108, 6), color="w", fontsize=8.5)
    for dy, nm, key in ((18, "top-75 circular mean", "retrieval top-75 circular mean"),
                        (-32, "ORACLE best pool window (emits 1.77 A)",
                         "ORACLE best pool window")):
        v = rg["channels"][key]["rms_deg"]
        ax.plot([1.0], [v], "^", ms=8, mfc="none", mec="w", mew=1.6)
        ax.annotate(f"{nm}\nRMS {v:.0f} deg", (1.0, v), textcoords="offset points",
                    xytext=(-172, dy), color="w", fontsize=8)

    # ---------------- Panel B
    ax = fig.add_subplot(gs[0, 1])
    for k, (nm, lab) in enumerate([
            ("iid_terminal_gap", "gaps at the chain ENDS"),
            ("iid_contig_gap", "one contiguous interior gap"),
            ("primary_iid_uniform_pool", "gaps scattered uniformly"),
            ("iid_mid_gap", "gap centred mid-chain")]):
        g = np.asarray(ph["surfaces"][nm]["grid"])[0]
        ax.plot(cov, g, "-o", ms=4, color=PAL[k], label=lab)
    for lv, c in ((3.0, "#EF476F"), (2.5, "#D6800B"), (2.0, "#1B7837")):
        ax.axhline(lv, ls=":", lw=1, color=c)
        ax.text(0.005, lv + 0.03, f"{lv:.1f} A", fontsize=8, color=c)
    ax.set_xlabel("coverage")
    ax.set_ylabel("CA-RMSD (A)")
    ax.set_title("B. PERFECT torsions on the covered residues (sigma = 0):\n"
                 "where the gaps fall is worth up to 1.79 A at identical coverage",
                 fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25)

    # ---------------- Panel C
    ax = fig.add_subplot(gs[1, 0])
    names, vals = [], []
    for nm in ("primary_iid_uniform_pool", "ar1_0.5", "ar1_0.8", "resclass", "ssbias"):
        s = ph["surfaces"][nm]
        names.append({"primary_iid_uniform_pool": "i.i.d.", "ar1_0.5": "AR(1) rho=0.5",
                      "ar1_0.8": "AR(1) rho=0.8", "resclass": "residue-class bias",
                      "ssbias": "SS-class bias"}[nm])
        vals.append([s["contours"][str(lv)][-1] for lv in (3.0, 2.5, 2.0)])
    vals = np.asarray(vals)
    w = 0.26
    x = np.arange(len(names))
    for k, (lv, c) in enumerate(((3.0, "#EF476F"), (2.5, "#D6800B"), (2.0, "#1B7837"))):
        ax.bar(x + (k - 1) * w, vals[:, k], w, color=c, label=f"{lv:.1f} A")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8.5, rotation=12)
    ax.set_ylabel("per-torsion sigma required (deg)")
    ax.set_title("C. Error SHAPE at full coverage: coherent error is more expensive\n"
                 "than i.i.d. error of the same magnitude", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, title="target RMSD", loc="upper right", ncol=3)
    ax.grid(alpha=0.25, axis="y")
    ax.set_ylim(0, 30)

    # ---------------- Panel D
    ax = fig.add_subplot(gs[1, 1])
    al = nl["alignment"]
    order = sorted(al, key=lambda k: al[k]["mean"])
    y = np.arange(len(order))
    m = [al[k]["mean"] for k in order]
    lo = [al[k]["mean"] - al[k]["ci95"][0] for k in order]
    hi = [al[k]["ci95"][1] - al[k]["mean"] for k in order]
    ax.barh(y, m, xerr=[lo, hi], color="#2166AC", height=0.62, error_kw=dict(lw=1.2))
    ax.axvline(nl["null_random_direction"]["mean"], color=C_LINE, ls="--", lw=1.4)
    ax.text(nl["null_random_direction"]["mean"], -0.9,
            f"random direction {nl['null_random_direction']['mean']:.3f}",
            fontsize=8, ha="center")
    short = {"retrieval pool-500 circular mean": "pool-500 circ. mean",
             "retrieval top-75 circular mean": "top-75 circ. mean",
             "retrieval top-75 sim-weighted": "top-75 sim-weighted",
             "one random top-75 window": "one random top-75 window",
             "top-75 medoid window": "top-75 medoid window",
             "constant alpha-helix (control)": "constant alpha-helix",
             "incumbent projected torsions": "incumbent projection",
             "ORACLE best pool window": "ORACLE best pool window"}
    ax.set_yticks(y); ax.set_yticklabels([short.get(k, k) for k in order], fontsize=8)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("alignment = RMSD damage per unit angular error, vs a random direction")
    ax.set_title("D. Real torsion errors live in the RMSD-QUIET subspace, so panel A\n"
                 "over-prices them; every CI excludes the random-direction null",
                 fontsize=10, loc="left")
    ax.grid(alpha=0.25, axis="x")

    out = os.path.join(L.FIG, "info_phase_diagram.png")
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    print("wrote", out)
    return out


if __name__ == "__main__":
    main()
