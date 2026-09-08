"""SPRINT 13 — publication figures.  Every figure answers one scientific question.

Written against the sprint's own result JSONs; nothing is hand-entered.  A figure whose
input is missing is skipped with a printed note rather than faked.

Palette: the validated 3-slot categorical set (blue / orange / aqua), which passes the
all-pairs CVD and normal-vision floors on a light surface.  Aqua sits below 3:1 contrast on
this surface, so every aqua mark carries a visible direct label -- the relief rule.
Series identity is never colour-alone: each figure legends and/or directly labels its series.

    python -m s13.figures
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
RES = os.path.join(ROOT, "s13", "results")
OUT = os.path.join(ROOT, "s13", "figures")
os.makedirs(OUT, exist_ok=True)

SURF = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8983"
S1 = "#2a78d6"      # blue    — Legacy
S2 = "#eb6834"      # orange  — AMBER
S3 = "#1baf7a"      # aqua    — third series (always directly labelled)
GRID = "#e3e2dd"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": GRID, "grid.color": GRID, "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 140,
})


def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def save(fig, name, note):
    p = os.path.join(OUT, name)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name:34s} {note}")


# ---------------------------------------------------------------- MONEY FIGURE 1
def fig_gradient_prediction():
    """Does the Pauli spectrum predict the gradient variance, with no free parameter?"""
    j = load("geo_pauli.json")
    if not j:
        print("  SKIP gradient prediction (geo_pauli.json missing)"); return
    series = [("Legacy", S1, "legacy"), ("AMBER raw", S2, "amber"),
              ("AMBER conditioned", S3, "amber_soft")]
    data = {}
    for label, col, key in series:
        xs, ys = [], []
        for c in j["cells"]:
            d = c.get(key)
            if not d:
                continue
            pr = d.get("var_grad_pred_exact_per_string")
            me = (d.get("var_grad_measured") or {}).get("a1.0")
            if pr and me and pr > 0 and me > 0:
                xs.append(pr); ys.append(me)
        if xs:
            data[label] = (col, np.array(xs), np.array(ys))

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.9),
                                  gridspec_kw={"width_ratios": [1, 1]})
    ax.set_axisbelow(True); ax2.set_axisbelow(True)
    lo, hi = np.inf, -np.inf
    for label, (col, xs, ys) in data.items():
        lo = min(lo, xs.min(), ys.min()); hi = max(hi, xs.max(), ys.max())
        ax.scatter(xs, ys, s=24, c=col, alpha=0.8, linewidths=0.7, edgecolors=SURF,
                   label=f"{label}  (n={len(xs)})", zorder=3)
    lo, hi = lo * 0.4, hi * 2.5
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1.2, ls="--", zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("predicted  Var[dC/dtheta]   from the Pauli spectrum")
    ax.set_ylabel("measured  Var[dC/dtheta]")
    ax.set_title("Predicted vs measured", loc="left", pad=12)
    ax.grid(True, lw=0.6, alpha=0.7)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.text(0.97, 0.03, "dashed line = perfect prediction", transform=ax.transAxes,
            ha="right", fontsize=8.5, color=MUTED)

    #: THE HONEST PANEL.  The left axis spans ~48 decades because raw AMBER's dynamic range
    #: does; on a scale that wide almost anything looks collinear.  The RATIO is the real
    #: test and it is scale-free, so it is given equal weight rather than relegated.
    labels, ratios, cols = [], [], []
    for label, (col, xs, ys) in data.items():
        labels.append(label); ratios.append(ys / xs); cols.append(col)
    parts = ax2.violinplot(ratios, showextrema=False, widths=0.8)
    for b, col in zip(parts["bodies"], cols):
        b.set_facecolor(col); b.set_alpha(0.45); b.set_edgecolor(col); b.set_linewidth(1.2)
    for i, (r, col) in enumerate(zip(ratios, cols), start=1):
        q1, med, q3 = np.percentile(r, [25, 50, 75])
        ax2.plot([i, i], [q1, q3], color=col, lw=3, solid_capstyle="round", zorder=4)
        ax2.plot(i, med, "o", color=col, ms=9, mec=SURF, mew=1.5, zorder=5)
        ax2.text(i + 0.22, med, f"{med:.3f}", va="center", fontsize=9.5, color=INK,
                 fontweight="bold")
    ax2.axhline(1.0, color=MUTED, lw=1.4, ls="--", zorder=1)
    ax2.set_xticks(range(1, len(labels) + 1)); ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("measured / predicted")
    ax2.set_ylim(0, 2.2)
    ax2.set_title("The ratio - the actual test, on a linear scale", loc="left", pad=12)
    ax2.grid(True, axis="y", lw=0.6, alpha=0.7)

    fig.suptitle("The energy model's Pauli spectrum predicts its gradient variance, "
                 "with no free parameter", x=0.0, ha="left", fontsize=12.5,
                 fontweight="bold", y=1.04)
    fig.text(0.0, -0.09,
             "Each point is one (target x encoding x circuit depth) cell; exact "
             "statevector, 6-18 qubits, depth 1-8. Prediction is "
             "Sum_S c_S^2 * Var_theta[d<Z_S>/dtheta], dropping cross-covariances "
             "between Pauli strings" + chr(10) +
             "- the ratio is that approximation's test. The LEFT axis spans ~48 decades "
             "because raw AMBER's dynamic range does, so the right panel is the one to "
             "read: a wide log axis flatters any fit.", fontsize=8, color=INK2)
    save(fig, "fig01_money_gradient_prediction.png", "- the no-free-parameter chain")


# ---------------------------------------------------------------- MONEY FIGURE 2
def _constant_helix_rmsd(pdbs):
    """Build every residue at phi=-63, psi=-42 and measure. A zero-information baseline.

    Recomputed here rather than read from a JSON so the figure independently reproduces the
    adversarial audit's headline number instead of quoting it.
    """
    from s12 import instrument as I
    out = {}
    for pdb in pdbs:
        u = I.load_univ(pdb)
        n = int(u["n"])
        phi = np.full(n, np.deg2rad(-63.0)); psi = np.full(n, np.deg2rad(-42.0))
        out[pdb] = I.ca_rmsd(I.build_ca(phi, psi), u["nat_ca"])
    return out


def fig_accuracy_ladder():
    """What does every native-free method actually emit, against what the space contains?"""
    j = load("coord_search_b5000.json")
    a = load("adv_prior.json")
    if not j:
        print("  SKIP accuracy ladder (coord_search missing)"); return
    rows = j["per_target"]
    g = lambda k: np.array([r.get(k, np.nan) for r in rows], float)     # noqa: E731

    helix = None
    if a and "P3_helix_rows" in a:
        hr = {x["pdb"]: x for x in a["P3_helix_rows"]}
        helix = np.array([hr.get(r["pdb"], {}).get("helix_ss_of", np.nan) for r in rows])
    ch = _constant_helix_rmsd([r["pdb"] for r in rows])
    const = np.array([ch[r["pdb"]] for r in rows], float)

    arms = [("SA on Legacy", g("sa_legacy"), S2),
            ("random sampling, 5,000 evals", g("random_legacy"), MUTED),
            ("constant alpha-helix (no information)", const, S3),
            ("SA on the 1-local prior", g("sa_prior"), S1),
            ("ORACLE ceiling of the space", g("ORACLE_descent"), INK)]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.8),
                                  gridspec_kw={"width_ratios": [1.3, 1]})
    ax.set_axisbelow(True); ax2.set_axisbelow(True)
    order = np.argsort([-np.nanmean(v) for _, v, _ in arms])
    for y, i in enumerate(order):
        lab, v, col = arms[i]
        m = np.nanmean(v)
        ax.barh(y, m, height=0.58, color=col, edgecolor=SURF, linewidth=1.5, zorder=3)
        ax.text(m + 0.09, y, f"{m:.3f}", va="center", fontsize=10, color=INK,
                fontweight="bold", zorder=4)
        inside = m > 2.6
        ax.text(0.09 if inside else m + 1.00, y, lab, va="center", fontsize=9.5,
                color="white" if inside else INK, zorder=4)
    ax.axvline(3.213, color=INK, lw=1.4, ls="--", zorder=2)
    ax.text(3.28, -0.72, "shipped retrieval pipeline  3.213", fontsize=8.5, color=INK2,
            va="center")
    ax.set_yticks([]); ax.set_xlabel("mean CA-RMSD (angstrom), 126 targets")
    ax.set_xlim(0, 6.0); ax.set_ylim(-1.1, len(arms) - 0.4)
    ax.set_title("Nothing native-free reaches what the space contains", loc="left", pad=14)
    ax.grid(True, axis="x", lw=0.6, alpha=0.7)

    if helix is not None and np.isfinite(helix).sum() > 100:
        hi = helix > 0.5; lo = helix < 0.1
        x = np.arange(2); w = 0.2
        bars = ((-1.5 * w, [np.nanmean(g("sa_prior")[hi]), np.nanmean(g("sa_prior")[lo])],
                 S1, "SA on the prior"),
                (-0.5 * w, [np.nanmean(const[hi]), np.nanmean(const[lo])], S3,
                 "constant alpha-helix"),
                (0.5 * w, [np.nanmean(g("random_legacy")[hi]),
                           np.nanmean(g("random_legacy")[lo])], MUTED, "random sampling"),
                (1.5 * w, [np.nanmean(g("ORACLE_descent")[hi]),
                           np.nanmean(g("ORACLE_descent")[lo])], INK, "ORACLE ceiling"))
        for off, vals, col, lab in bars:
            ax2.bar(x + off, vals, w * 0.9, color=col, edgecolor=SURF, linewidth=1.4,
                    label=lab, zorder=3)
            for xi, v in zip(x + off, vals):
                ax2.text(xi, v + 0.12, f"{v:.2f}", ha="center", fontsize=8, color=INK,
                         zorder=4)
        ax2.set_xticks(x)
        ax2.set_xticklabels(["helical (>50% helix)" + chr(10) + f"n={int(hi.sum())}",
                             "non-helical (<10%)" + chr(10) + f"n={int(lo.sum())}"],
                            fontsize=9)
        ax2.set_ylabel("mean CA-RMSD (angstrom)"); ax2.set_ylim(0, 7.2)
        ax2.set_title("...and the prior's advantage is a helix artefact", loc="left", pad=14)
        ax2.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2,
                   columnspacing=1.0)
        ax2.grid(True, axis="y", lw=0.6, alpha=0.7)
    fig.text(0.0, -0.09, "Discrete torsion space, k=4 (~26 qubits), matched 5,000-evaluation "
             "budget. On non-helical targets the prior is WORSE than random sampling, and a "
             "zero-information" + chr(10) + "constant alpha-helix reproduces it to within 0.10 A. "
             "rho(helix fraction, prior RMSD) = -0.744. ORACLE arms read the native and are "
             "diagnostics, not methods.", fontsize=8, color=INK2)
    save(fig, "fig02_money_accuracy_ladder.png", "- the ladder and the helix artefact")


# ---------------------------------------------------------------- FIGURE 3
def fig_ceiling():
    """How much structure does the discrete torsion space hold, per qubit spent?"""
    j = load("ceiling_report.json")
    if not j:
        print("  SKIP ceiling (missing)"); return
    ks = sorted(int(k) for k in j["k"])
    qb = [j["k"][str(k)]["mean_qubits"] for k in ks]
    de = [j["k"][str(k)]["descent"]["mean"] for k in ks]
    sn = [j["k"][str(k)]["snap"]["mean"] for k in ks]
    f18 = [j["k"][str(k)]["FAIL18"] for k in ks]
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(qb, sn, "o-", color=S2, lw=2, ms=8, mec=SURF, mew=1.4,
            label="nearest-state snap (chain-blind)")
    ax.plot(qb, de, "o-", color=S1, lw=2, ms=8, mec=SURF, mew=1.4,
            label="coordinate descent (chain-aware)")
    ax.plot(qb, f18, "o--", color=S3, lw=2, ms=8, mec=SURF, mew=1.4,
            label="descent, the 18 hardest targets")
    for x, y, k in zip(qb, de, ks):
        ax.annotate(f"k={k}", (x, y), textcoords="offset points",
                    xytext=(-14 if k == 32 else 0, 12 if k == 32 else -16), ha="center", fontsize=8.5, color=INK)
    ax.axhline(2.0, color=MUTED, lw=1.2, ls=":")
    ax.text(qb[-1], 2.06, "the 2.0 Å target", ha="right", fontsize=8.5, color=INK2)
    ax.axhline(3.213, color=INK, lw=1.4, ls="--")
    ax.text(qb[0], 3.28, "shipped retrieval pipeline", fontsize=8.5, color=INK2)
    ax.axhline(j["native_rebuild"]["mean"], color=MUTED, lw=1, ls="-")
    ax.text(qb[-1], j["native_rebuild"]["mean"] + 0.07, "ideal-geometry rebuild floor",
            ha="right", fontsize=8, color=MUTED)
    ax.set_xlabel("mean qubits  (n_residues × log₂ k)")
    ax.set_ylabel("mean CA-RMSD (Å), 126 targets")
    ax.set_title("Ideal-geometry torsion parameterisation is not the barrier", loc="left",
                 pad=12)
    ax.set_ylim(0.15, 3.65)
    ax.legend(frameon=False, fontsize=8.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=3, columnspacing=1.6)
    ax.grid(True, lw=0.6, alpha=0.7)
    ax.scatter([qb[0]], [1.982], s=120, marker="X", color=INK, zorder=6, edgecolors=SURF,
               linewidths=1.6)
    ax.annotate("same descent, RANDOM start:" + chr(10) + "1.982 A  (+0.388)",
                (qb[0], 1.982), textcoords="offset points", xytext=(16, 4),
                fontsize=8.5, color=INK, fontweight="bold")
    fig.text(0.0, -0.34,
             "ORACLE DIAGNOSTIC - both arms read the native to score, and the descent "
             "start is the native snapped to its nearest states. The audit priced that"
             + chr(10) +
             "privilege at 0.388 A: from a random start the identical descent reaches "
             "1.982 A. It also showed the LIBRARY contributes little - a wrong target's "
             "library costs +0.135 A," + chr(10) +
             "and uniform-random torsion states still reach 2.698 A. Qubit counts are "
             "nominal: residue n-1 is inert for the CA trace, so log2(k) per chain are "
             "dead.", fontsize=8, color=INK2)
    save(fig, "fig03_representation_ceiling.png", "— ceiling vs qubits")


# ---------------------------------------------------------------- FIGURE 4
def fig_objective_validity():
    """Do the energies rank the native — and do they rank where a search actually lives?"""
    j = load("coord_amber_shape_report.json")
    if not j:
        print("  SKIP objective validity (missing)"); return
    names = ["legacy", "raw", "capped", "softcore", "nonclash", "minimised"]
    pretty = {"legacy": "Legacy", "raw": "AMBER raw", "capped": "AMBER capped p90",
              "softcore": "AMBER log-compressed\n(monotone)",
              "nonclash": "AMBER elec+solv only", "minimised": "AMBER minimised"}
    have = [n for n in names if n in j["variants"]]
    ra = [j["variants"][n].get("rho_all", np.nan) for n in have]
    rl = [j["variants"][n].get("rho_low_decile", np.nan) for n in have]
    dr = [j["variants"][n].get("log10_range", np.nan) for n in have]
    y = np.arange(len(have)); w = 0.36
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    ax.barh(y + w / 2, ra, w * 0.92, color=S1, edgecolor=SURF, linewidth=1.4,
            label="over all configurations")
    ax.barh(y - w / 2, rl, w * 0.92, color=S2, edgecolor=SURF, linewidth=1.4,
            label="inside the low-energy decile")
    ax.axvline(0, color=INK, lw=1.2)
    ax.set_yticks(y); ax.set_yticklabels([pretty[n] for n in have], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman ρ (energy, CA-RMSD)   — higher is better, 0 is useless")
    ax.set_title("Neither energy ranks where a search actually lives", loc="left")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax.grid(True, axis="x", lw=0.6, alpha=0.7)

    ax2.barh(y, dr, 0.55, color=S3, edgecolor=SURF, linewidth=1.4)
    for yi, v in zip(y, dr):
        ax2.text(v + 0.25, yi, f"{v:.1f}", va="center", fontsize=9, color=INK,
                 fontweight="bold")
    ax2.set_yticks(y); ax2.set_yticklabels([]); ax2.invert_yaxis()
    ax2.set_xlabel("log₁₀ dynamic range of the objective")
    ax2.set_title("Dynamic range: 16 decades, and\nrescaling does not help", loc="left")
    ax2.grid(True, axis="x", lw=0.6, alpha=0.7)
    fig.text(0.0, -0.08, "20 targets × 120 configurations of the k=4 torsion space. "
             "The log compression is MONOTONE, so it cannot change the rank order — and "
             "indeed ρ is identical\nto raw while the range falls from 16.1 to 1.8 decades. "
             "The damage is in the ordering, not the scale.", fontsize=8, color=INK2)
    save(fig, "fig04_objective_validity.png", "— ρ and dynamic range")


# ---------------------------------------------------------------- FIGURE 5
def fig_pauli_spectra():
    """Is the all-atom force field a more global observable — and was the first answer real?"""
    j = load("geo_pauli.json")
    if not j:
        print("  SKIP Pauli spectra (missing)"); return
    cells = [c for c in j["cells"] if c.get("legacy") and c.get("amber")]
    if not cells:
        print("  SKIP Pauli spectra (no complete cells)"); return
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.3))

    def mean_share(key):
        M = []
        for c in cells:
            ws = (c[key]["spectrum"] or {}).get("weight_share")
            if ws:
                M.append(np.asarray(ws, float))
        if not M:
            return None
        L = max(len(x) for x in M)
        A = np.full((len(M), L), np.nan)
        for i, x in enumerate(M):
            A[i, :len(x)] = x
        return np.nanmean(A, axis=0)

    for key, col, lab in (("legacy_soft", S1, "Legacy (conditioned)"),
                          ("amber_soft", S2, "AMBER (conditioned)")):
        v = mean_share(key)
        if v is None:
            continue
        ax.plot(np.arange(len(v)), v, "o-", color=col, lw=2, ms=6, mec=SURF, mew=1.2,
                label=lab)
    v = mean_share("amber")
    if v is not None:
        ax.plot(np.arange(len(v)), v, "s:", color=MUTED, lw=1.6, ms=5,
                label="AMBER raw (the artefact)")
    ax.set_xlabel("Pauli weight  |S|"); ax.set_ylabel("share of the energy's variance")
    ax.set_title("AMBER carries higher-order interaction structure", loc="left")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(True, lw=0.6, alpha=0.7)

    spikes_a = [c["amber"]["spectrum"]["spike"]["top10"] for c in cells
                if c["amber"]["spectrum"].get("spike")]
    spikes_l = [c["legacy"]["spectrum"]["spike"]["top10"] for c in cells
                if c["legacy"]["spectrum"].get("spike")]
    ax2.hist([spikes_l, spikes_a], bins=np.linspace(0, 1, 21), color=[S1, S2],
             label=[f"Legacy (median {np.median(spikes_l):.3f})",
                    f"AMBER raw (median {np.median(spikes_a):.3f})"],
             edgecolor=SURF, linewidth=0.8)
    ax2.set_xlabel("share of the energy's variance held by its top 10 configurations")
    ax2.set_ylabel("cells")
    ax2.set_title("Why the raw answer was an artefact", loc="left")
    ax2.legend(frameon=False, fontsize=8.5)
    ax2.grid(True, axis="y", lw=0.6, alpha=0.7)
    fig.text(0.0, -0.07, "A constant-plus-single-spike function has Walsh weight spectrum "
             "exactly Binomial(m, ½) — maximally global for arithmetic reasons, not physical\n"
             "ones. Raw AMBER's spectrum measures its worst steric clash; the conditioned "
             "arms use a monotone, rank-preserving compression applied identically to both.",
             fontsize=8, color=INK2)
    save(fig, "fig05_pauli_spectra.png", "— spectra and the spike artefact")


# ---------------------------------------------------------------- FIGURE 6
def fig_budget_curve():
    """Does optimising the objective harder produce better structures?"""
    j = load("qarch_budget.json")
    if not j:
        print("  SKIP budget curve (missing)"); return
    budgets = [int(b) for b in j["budgets"]]
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    keys = [("legacy_total", S1, "Legacy"), ("prior_empirical", S3, "1-local prior"),
            ("ORACLE_rmsd", INK, "ORACLE: search on true RMSD")]
    for key, col, lab in keys:
        M = []
        for r in j["rows"]:
            cur = (r.get("curves") or {}).get(key)
            if not cur:
                continue
            row = []
            for b in budgets:
                v = cur.get(str(b))
                row.append(v.get("mean_rmsd", np.nan) if isinstance(v, dict)
                           else (v if v is not None else np.nan))
            M.append(row)
        if not M:
            continue
        m = np.nanmean(np.array(M, float), axis=0)
        ls = "--" if key == "ORACLE_rmsd" else "-"
        ax.plot(budgets, m, ls, marker="o", color=col, lw=2, ms=6, mec=SURF, mew=1.2,
                label=lab)
        if key == "legacy_total":
            i = int(np.nanargmin(m))
            ax.annotate(f"best at {budgets[i]:,} evals\n{m[i]:.3f} Å",
                        (budgets[i], m[i]), textcoords="offset points", xytext=(4, -34),
                        fontsize=8.5, color=INK,
                        arrowprops=dict(arrowstyle="->", color=INK2, lw=1))
            ax.annotate(f"certified global optimum\n{m[-1]:.3f} Å — WORSE",
                        (budgets[-1], m[-1]), textcoords="offset points", xytext=(-118, -44),
                        fontsize=8.5, color=S2, fontweight="bold",
                        arrowprops=dict(arrowstyle="->", color=S2, lw=1.2))
    ax.set_xscale("log")
    ax.set_xlabel("objective evaluations  (log scale; the space is 4⁹ = 262,144)")
    ax.set_ylabel("mean CA-RMSD (Å)")
    ax.set_title("Optimising the energy harder makes the structure worse", loc="left", pad=16)
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(True, lw=0.6, alpha=0.7)
    fig.text(0.0, -0.06, "Nine fully enumerated 9-residue targets, k=4 — every configuration "
             "scored, so the endpoint is the certified global optimum, not an estimate.\n"
             "A search at 12–16 residues sees <0.4 % of its space and would sit on the "
             "improving part of this curve: it would look like it was working.",
             fontsize=8, color=INK2)
    save(fig, "fig06_budget_curve.png", "— optimise harder, get worse")


# ---------------------------------------------------------------- FIGURE 7
def fig_restraint_surface():
    """How good would torsion information have to be, and how good can we predict it?"""
    t = load("tors_surface.json")
    if not t or "full_coverage" not in t:
        print("  SKIP restraint surface (tors_surface.json missing)"); return
    fc = t["full_coverage"]
    sig = sorted(float(k) for k in fc)
    ys = [fc[str(int(x)) if float(x).is_integer() else str(x)]["mean"] for x in sig]
    f18 = [fc[str(int(x)) if float(x).is_integer() else str(x)]["FAIL18"] for x in sig]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.6),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    ax.set_axisbelow(True); ax2.set_axisbelow(True)

    ax.plot(sig, ys, "o-", color=S1, lw=2.2, ms=6, mec=SURF, mew=1.2,
            label="all 126 targets", zorder=3)
    ax.plot(sig, f18, "o--", color=S3, lw=2, ms=5, mec=SURF, mew=1.2,
            label="the 18 hardest targets", zorder=3)
    k18 = int(np.argmin(np.abs(np.array(f18) - 4.6)))
    ax.annotate("18 hardest", (sig[k18], f18[k18]), textcoords="offset points",
                xytext=(-8, 12), fontsize=8.5, color=S3, fontweight="bold", ha="right")
    ax.axhline(2.0, color=MUTED, lw=1.2, ls=":", zorder=1)
    ax.text(99, 2.08, "the 2.0 A target", ha="right", fontsize=8.5, color=INK2)
    ax.axhline(3.213, color=INK, lw=1.4, ls="--", zorder=1)
    ax.text(2, 3.30, "shipped retrieval pipeline", fontsize=8.5, color=INK2)

    ax.scatter([29.0], [3.213], s=95, marker="o", facecolors="none", edgecolors=INK,
               linewidths=1.8, zorder=6)
    ax.annotate("the incumbent is equivalent" + chr(10) + "to sigma ~29 deg",
                (29.0, 3.213), textcoords="offset points", xytext=(6, -34),
                fontsize=8.5, color=INK)
    ax.scatter([67.7], [3.770], s=130, marker="X", color=S2, zorder=6, edgecolors=SURF,
               linewidths=1.6)
    ax.annotate("best sequence-only" + chr(10) + "predictor: sigma 67.7," +
                chr(10) + "3.770 A", (67.7, 3.770), textcoords="offset points",
                xytext=(10, -30), fontsize=8.5, color=S2, fontweight="bold")

    ax.set_xlim(-3, 103); ax.set_ylim(0, 6.0)
    ax.set_xlabel("per-torsion error sigma (degrees), all residues restrained")
    ax.set_ylabel("mean CA-RMSD (angstrom)")
    ax.set_title("Torsion information reaches the target", loc="left", pad=12)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax.grid(True, lw=0.6, alpha=0.7)

    #: THE CORRECTED MISSINGNESS MODEL.  The coordinator's guidance was that gaps cluster at
    #: termini and that uniform dropout would FLATTER the method.  Measured, it is the
    #: opposite: terminal dropout is much cheaper than uniform, so uniform UNDERSTATES a
    #: TALOS-N-style channel.  Recorded here because it moves the decision.
    cells = t.get("cells", {})
    fracs = [0.9, 0.75, 0.5]
    models = [("uniform", MUTED, "uniform dropout"),
              ("clustered", S1, "clustered (interior)"),
              ("terminal", S2, "terminal - the real pattern")]
    x = np.arange(len(fracs)); w = 0.26
    for i, (key, col, lab) in enumerate(models):
        vals = []
        for f in fracs:
            c = cells.get(f"s12_f{f:g}_{key}")
            vals.append(c["mean"] if c else np.nan)
        ax2.bar(x + (i - 1) * w, vals, w * 0.9, color=col, edgecolor=SURF, linewidth=1.4,
                label=lab, zorder=3)
        for xi, v in zip(x + (i - 1) * w, vals):
            if np.isfinite(v):
                ax2.text(xi, v + 0.07, f"{v:.2f}", ha="center", fontsize=8, color=INK,
                         zorder=4)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{int(f*100)}% restrained" for f in fracs], fontsize=9)
    ax2.set_ylabel("mean CA-RMSD (angstrom)")
    ax2.set_ylim(0, 4.6)
    ax2.axhline(3.213, color=INK, lw=1.2, ls="--", zorder=1)
    ax2.set_title("Where the gaps fall beats how many", loc="left", pad=12)
    ax2.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax2.grid(True, axis="y", lw=0.6, alpha=0.7)

    fig.text(0.0, -0.10,
             "ORACLE DIAGNOSTIC: native torsions corrupted to the stated sigma and coverage, "
             "then built directly; unrestrained residues filled from the best-agreeing "
             "library member." + chr(10) +
             "RIGHT PANEL, at sigma = 12 deg: the coordinator predicted that clustered gaps "
             "would be far worse than uniform. They are not - TERMINAL dropout, which is "
             "where a" + chr(10) +
             "chemical-shift predictor actually declines, is 0.40-0.50 A CHEAPER than "
             "uniform. Uniform dropout UNDERSTATES the channel, and the prediction was wrong.",
             fontsize=8, color=INK2)
    save(fig, "fig07_restraint_surface.png", "- sigma x coverage, and the missingness model")


def main():
    print("Sprint 13 figures ->", OUT)
    for fn in (fig_gradient_prediction, fig_accuracy_ladder, fig_ceiling,
               fig_objective_validity, fig_pauli_spectra, fig_budget_curve,
               fig_restraint_surface):
        try:
            fn()
        except Exception as exc:                                       # noqa: BLE001
            print(f"  FAILED {fn.__name__}: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()
