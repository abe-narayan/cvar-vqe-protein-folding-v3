"""SPRINT 15 / QGEOM -- the QNG verdict, computed from the checkpoints.

Turns `qgeom_qng`'s raw cells into the paired statistics the brief requires: for each pair of
optimisers, the paired mean difference in objective gap with a bootstrap CI, W/L, mean/sd and
a null-calibrated concentration read, split by CONDITIONING BAND so the question
"does QNG help where the metric is ill-conditioned?" is answered as a function of the
condition number rather than on a pooled average.

    python -m s15.qgeom_qngreport
"""
from __future__ import annotations

import numpy as np

from s15 import qgeom_lib as G

TAG = "qngrep"
ARMS = ("sgd", "adam", "qng", "qng_diag", "qng_adam", "ng_shots")
BANDS = ((0, 2, "well-conditioned (cond ~ 1, g = I/4)"),
         (2, 100, "mild (2 - 100)"),
         (100, 1e9, "ill-conditioned (> 100)"))


def _rows(key, cellfield="obj_gap"):
    d = G.ck_load("qng").get(key, {})
    out = []
    for k, v in d.items():
        out.append({"cell": k, "cond": v.get("cond", np.nan),
                    "vals": v[cellfield],
                    "n_params": v.get("n_params"),
                    "iters": v.get("iters")})
    return out


def band_of(c):
    for lo, hi, nm in BANDS:
        if lo <= c < hi:
            return nm
    return BANDS[-1][2]


def report(key, label, cellfield="obj_gap", lower_is_better=True):
    rows = _rows(key, cellfield)
    if not rows:
        print(f"  ({key}: no cells yet)")
        return {}
    print()
    print("=" * 112)
    print(f"{label}   [{key}]   n = {len(rows)} cells")
    print("=" * 112)
    out = {}
    # condition numbers are missing in the equal-cost table; take them from A4
    condmap = {r["cell"]: r["cond"] for r in _rows("A4_vqe")}
    for r in rows:
        if not np.isfinite(r["cond"]):
            r["cond"] = condmap.get(r["cell"], np.nan)
    for lo, hi, nm in BANDS:
        sel = [r for r in rows if np.isfinite(r["cond"]) and lo <= r["cond"] < hi]
        if not sel:
            continue
        print(f"\n  --- {nm} --- n = {len(sel)}")
        means = {a: float(np.mean([r["vals"][a] for r in sel])) for a in ARMS}
        print("      mean objective gap: " +
              "  ".join(f"{a}={means[a]:.5f}" for a in ARMS))
        best = min(ARMS, key=lambda a: means[a])
        print(f"      best arm on the mean: {best}")
        for a, b in (("qng", "sgd"), ("qng", "adam"), ("qng_adam", "adam"),
                     ("qng_diag", "sgd"), ("ng_shots", "sgd"), ("adam", "sgd")):
            pr = G.paired([r["vals"][a] for r in sel], [r["vals"][b] for r in sel])
            out[f"{key}|{nm}|{a}-{b}"] = pr
            print(f"      {a + ' - ' + b:22s} {pr['mean']:+10.5f} "
                  f"[{pr['ci_lo']:+.5f}, {pr['ci_hi']:+.5f}]  W/L {pr['win']}/{pr['loss']}"
                  f"  mean/sd {pr['mean_over_sd']:+6.2f}  {pr['verdict']}")
        # null-calibrated concentration: is the effect carried by a few cells?
        for a, b in (("qng", "adam"),):
            d = np.array([r["vals"][a] - r["vals"][b] for r in sel])
            share = float(np.sort(np.abs(d))[::-1][:max(1, len(d) // 10)].sum()
                          / max(np.abs(d).sum(), 1e-30))
            print(f"      CONCENTRATION({a}-{b}): top-10% of cells carry {share:.3f} of "
                  f"|effect|; uniform expectation {max(1, len(d)//10)/len(d):.3f} -> "
                  f"{'CONCENTRATED' if share > 3 * max(1, len(d)//10)/len(d) else 'DIFFUSE'}")
    G.ck(TAG, key, out)
    return out


def main():
    print("=" * 112)
    print("QNG VERDICT.  Lower objective gap is better.  Negative paired diff = first arm")
    print("better.  W = first arm better.  Bands are the MEASURED per-point condition")
    print("number of the Fubini-Study metric at the initialisation.")
    print("=" * 112)
    report("A4_vqe", "A4. THE REAL VQE, EQUAL ITERATIONS (convention A)")
    report("A4b_vqe_equal_cost",
           "A4b. THE REAL VQE, EQUAL HARDWARE COST (convention B)")
    # structural axis
    print()
    print("=" * 112)
    print("A4c. THE STRUCTURAL AXIS -- mode RMSD, and the budget trap")
    print("=" * 112)
    rows = _rows("A4_vqe", "mode_rmsd")
    raw = G.ck_load("qng").get("A4_vqe", {})
    if rows:
        for a in ARMS:
            v = [r["vals"][a] for r in rows]
            print(f"    {a:10s} mean mode RMSD {np.mean(v):6.3f} A")
        cert = [raw[r["cell"]]["certified_rmsd"] for r in rows]
        rnd = [raw[r["cell"]]["random_draw"] for r in rows]
        print(f"    {'certified optimum of the objective':34s} {np.mean(cert):6.3f} A")
        print(f"    {'uniform random draw':34s} {np.mean(rnd):6.3f} A")
        print()
        print("    If the arms beat the CERTIFIED OPTIMUM of the objective they are")
        print("    optimising, then failing to optimise is what produced the structure --")
        print("    the Sprint 14 budget trap, reproduced here on the optimiser axis.")
        gaps = _rows("A4_vqe", "obj_gap")
        for a in ARMS:
            g = np.array([r["vals"][a] for r in gaps])
            m = np.array([r["vals"][a] for r in rows])
            print(f"    rho(objective gap, mode RMSD) for {a:10s} = "
                  f"{G.spearman(g, m):+.3f}")
    print()
    print("done -> s15/results/qgeom_qngrep.json")


if __name__ == "__main__":
    main()
