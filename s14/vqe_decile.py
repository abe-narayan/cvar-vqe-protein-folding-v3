"""SPRINT 14 / VQE -- put the signal family on the coordinator's axis.

The coordinator's structural objective is reported as **in-decile** rank correlation
(+0.370 at w=0.25, against Legacy's +0.043 and raw AMBER's -0.088).  My signal-tunable
family was calibrated on GLOBAL rho, and the two are very different statistics: a global
rho is dominated by the bulk, an in-decile rho by the region a search actually lives in.
Comparing my curve to the coordinator's curve therefore requires the family to be
re-expressed on the in-decile axis, which is what this does.

It also computes the SELECTION GAP the coordinator asked for -- ``RMSD_returned`` minus
``RMSD of the best configuration the arm actually evaluated`` -- for every arm of the Part C
grid, recovered post hoc from the recorded ``rmsd_returned`` and ``rmsd_best_seen``, so no
rerun is needed.  The selection gap separates two failure modes that look identical in a
flat curve: an arm that never SAW anything better (budget-limited) from an arm that saw
something better and could not tell (discrimination-limited).

    python -m s14.vqe_decile
"""
from __future__ import annotations

import json
import os

import numpy as np

from s14 import vqe_lib as V

SIGNALS = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0)


def in_decile_rho(E: np.ndarray, rmsd: np.ndarray, frac=0.1) -> float:
    """rho(objective, RMSD) inside the lowest `frac` of the objective."""
    k = max(3, int(round(frac * len(E))))
    idx = np.argpartition(E, k - 1)[:k]
    return V.spearman(E[idx], rmsd[idx])


def main():
    E9 = V.load_all()
    print("=" * 84)
    print("1. THE SIGNAL FAMILY ON THE IN-DECILE AXIS (the coordinator's statistic)")
    print("   rho(objective, true RMSD) inside the objective's own lowest decile.")
    print("   Sign convention: POSITIVE means the objective ranks correctly.")
    print("=" * 84)
    out = {}
    for base in ("legacy", "prior"):
        print(f"\n--- base = {base} ---")
        print(f"{'signal':>7s} {'global rho':>11s} {'in-decile rho':>14s} "
              f"{'RMSD @ argmin':>14s} {'best in decile':>15s}")
        for s in SIGNALS:
            g, d, am, bd = [], [], [], []
            for p, e in E9.items():
                E = V.blend_objective(getattr(e, base), e.rmsd, s)
                g.append(V.spearman(E, e.rmsd))
                d.append(in_decile_rho(E, e.rmsd))
                am.append(e.rmsd[int(np.argmin(E))])
                k = int(0.1 * len(E))
                bd.append(e.rmsd[np.argpartition(E, k - 1)[:k]].min())
            out.setdefault(base, {})[str(s)] = {
                "global_rho": float(np.mean(g)), "in_decile_rho": float(np.mean(d)),
                "rmsd_at_argmin": float(np.mean(am)),
                "best_in_decile": float(np.mean(bd)),
                "per_target_in_decile": d}
            print(f"{s:7.2f} {np.mean(g):+11.3f} {np.mean(d):+14.3f} "
                  f"{np.mean(am):14.3f} {np.mean(bd):15.3f}")

    # where does the coordinator's +0.370 land on my family?
    print()
    print("=" * 84)
    print("2. LOCATING THE COORDINATOR'S OBJECTIVE ON THIS FAMILY")
    print("=" * 84)
    for base in ("legacy", "prior"):
        xs = [out[base][str(s)]["in_decile_rho"] for s in SIGNALS]
        for tgt, label in ((0.043, "Legacy (s13 measurement)"),
                           (0.370, "coordinator's structural objective, w=0.25")):
            i = int(np.argmin(np.abs(np.array(xs) - tgt)))
            # linear interpolation in signal
            sig = np.interp(tgt, xs, SIGNALS) if xs[0] <= tgt <= xs[-1] else SIGNALS[i]
            print(f"  base={base:7s} in-decile rho {tgt:+.3f} ({label}) "
                  f"-> signal ~ {sig:.3f}  (nearest grid point {SIGNALS[i]:.2f}, "
                  f"rho {xs[i]:+.3f})")
            out.setdefault("locate", {})[f"{base}|{tgt}"] = float(sig)

    # reference: the raw objectives themselves
    print()
    print("  reference in-decile rho of the raw objectives on these nine targets:")
    for name in ("legacy", "prior"):
        v = [in_decile_rho(getattr(e, name), e.rmsd) for e in E9.values()]
        print(f"    {name:8s} {np.mean(v):+.4f}   per target "
              f"{' '.join(f'{x:+.2f}' for x in v)}")
    v = [in_decile_rho(-getattr(e, "legacy"), e.rmsd) for e in E9.values()]
    print(f"    {'-legacy':8s} {np.mean(v):+.4f}  (sign control)")

    V.write("vqe_decile", out)

    # ---------------------------------------------------------------- selection gap
    path = os.path.join(V.RESULTS, "vqe_signal_legacy.json")
    if os.path.exists(path):
        selection_gap(path)
    print("\nwritten -> s14/results/vqe_decile.json")


ARMS = ("random", "greedy", "anneal", "ga", "vqe_a1.0", "vqe_a0.25", "vqe_a0.1",
        "vqe_a0.05")


def selection_gap(path):
    """RMSD returned minus RMSD of the best configuration the arm actually evaluated."""
    d = json.load(open(path))
    cells = d["cells"]
    print()
    print("=" * 84)
    print("3. SELECTION GAP per arm  (RMSD_returned - RMSD_best_EVALUATED)")
    print("   Separates budget-limited from discrimination-limited.  Averaged over")
    print("   whatever cells have completed so far.")
    print("=" * 84)
    sig_of = lambda k: float(k.split("|")[2])
    buds = sorted({int(k.split("|")[3]) for k in cells})
    sigs = sorted({sig_of(k) for k in cells})
    for bud in buds:
        print(f"\n  budget = {bud} evaluations")
        print(f"  {'signal':>7s} " + "".join(f"{a:>11s}" for a in ARMS))
        for s in sigs:
            ks = [k for k in cells if sig_of(k) == s and int(k.split("|")[3]) == bud]
            if not ks:
                continue
            row = []
            for a in ARMS:
                g = [cells[k][a]["rmsd_returned"] - cells[k][a]["rmsd_best_seen"]
                     for k in ks if a in cells[k]
                     and "rmsd_best_seen" in cells[k][a]]
                row.append(np.mean(g) if g else np.nan)
            print(f"  {s:7.2f} " + "".join(f"{v:11.3f}" for v in row))
        print(f"  (n = {len(ks)} targets)")


if __name__ == "__main__":
    main()
