"""SPRINT 14 / VQE -- the hamil w axis, finely, and a LIKE-FOR-LIKE selection gap.

Two questions from the coordinator.

1. The certified argmin is still improving as w falls toward 0.  Where does it stop?  If it
   keeps improving below w=0.25 the shipped arm may be the retrieval torsion prior almost
   alone, which is a much simpler architecture.

2. Is "+1.876 certified (9 targets) vs +2.040 sampled (126 targets) agree to 0.16 A" real
   corroboration, or two different statistics that happen to land close?

   They are NOT the same statistic.  Mine is
       certified argmin RMSD  -  best RMSD inside the objective's own lowest DECILE
   and the coordinator's is
       returned RMSD          -  best RMSD in the SAMPLE actually drawn.
   Different reference sets (a decile of the complete space vs a 20,000-draw sample) and
   different target populations (nine n=9 targets, the length bin s13 shows is the most
   thoroughly explored, vs 126 targets spanning n=9-16).  So this module computes the
   COORDINATOR'S OWN statistic on the enumerated nine, which makes the comparison
   like-for-like and settles whether the agreement is real.

    python -m s14.vqe_wsweep
"""
from __future__ import annotations

import numpy as np

from s14 import vqe_lib as V
from s14.vqe_hamil import tabulate, combine, in_decile_rho

WS = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.35, 0.50, 0.75, 1.0)


def main():
    rows = {p: tabulate(p) for p in V.ENUM_TARGETS}
    E9 = {p: V.Enum(p) for p in V.ENUM_TARGETS}

    print("=" * 104)
    print("1. THE w AXIS, FINELY -- does the certified optimum keep improving below 0.25?")
    print("   Nine fully enumerated targets, all 262,144 configurations, uniform population.")
    print("=" * 104)
    print(f"  {'w':>5s} {'global rho':>11s} {'in-decile':>10s} "
          f"{'CERTIFIED argmin':>17s} {'sd':>6s} {'decile mean':>12s} {'decile best':>12s} "
          f"{'sel gap':>8s}")
    out = {"w_sweep": {}}
    for w in WS:
        g, dr, am, dm, db = [], [], [], [], []
        for p in V.ENUM_TARGETS:
            e = E9[p]
            E = combine(*rows[p], w)
            g.append(V.spearman(E, e.rmsd))
            r, mn, bs = in_decile_rho(E, e.rmsd)
            dr.append(r); dm.append(mn); db.append(bs)
            am.append(float(e.rmsd[int(np.argmin(E))]))
        out["w_sweep"][str(w)] = {"global_rho": float(np.mean(g)),
                                  "in_decile_rho": float(np.mean(dr)),
                                  "argmin": float(np.mean(am)),
                                  "argmin_sd": float(np.std(am)),
                                  "decile_mean": float(np.mean(dm)),
                                  "decile_best": float(np.mean(db)),
                                  "per_target_argmin": {p: a for p, a in
                                                        zip(V.ENUM_TARGETS, am)}}
        print(f"  {w:5.2f} {np.mean(g):+11.3f} {np.mean(dr):+10.3f} "
              f"{np.mean(am):17.3f} {np.std(am):6.3f} {np.mean(dm):12.3f} "
              f"{np.mean(db):12.3f} {np.mean(am)-np.mean(db):8.3f}")
    best_w = min(WS, key=lambda w: out["w_sweep"][str(w)]["argmin"])
    print(f"\n  BEST certified argmin at w = {best_w} "
          f"({out['w_sweep'][str(best_w)]['argmin']:.3f} A).")
    print(f"  reference: random draw {np.mean([e.rmsd.mean() for e in E9.values()]):.3f}, "
          f"space best {np.mean([e.rmsd.min() for e in E9.values()]):.3f}, "
          f"Legacy certified argmin 3.920, s13 torsion prior 3.636")
    print("\n  per-target certified argmin (the mean is a MIXTURE, not a typical target):")
    print(f"  {'w':>5s} " + "".join(f"{p:>8s}" for p in V.ENUM_TARGETS))
    for w in (0.0, 0.10, 0.25, 1.0):
        a = out["w_sweep"][str(w)]["per_target_argmin"]
        print(f"  {w:5.2f} " + "".join(f"{a[p]:8.3f}" for p in V.ENUM_TARGETS))
    print(f"  {'rand':>5s} " + "".join(f"{E9[p].rmsd.mean():8.3f}" for p in V.ENUM_TARGETS))

    # ---------------------------------------------------------- like-for-like gap
    print()
    print("=" * 104)
    print("2. THE COORDINATOR'S OWN STATISTIC, computed on the enumerated nine")
    print("   selection gap = RMSD returned (argmin of the sample) - best RMSD IN THE SAMPLE")
    print("   uniform draws, read as increasing prefixes so budgets are nested (paired).")
    print("=" * 104)
    budgets = (10, 300, 3000, 20000)
    print(f"  {'w':>5s} " + "".join(f"{f'b={b}':>26s}" for b in budgets))
    print(f"  {'':>5s} " + "".join(f"{'ret':>8s}{'best':>9s}{'gap':>9s}"
                                   for b in budgets))
    out["like_for_like"] = {}
    for w in (0.0, 0.10, 0.25, 1.0):
        cells = {b: {"ret": [], "best": []} for b in budgets}
        for p in V.ENUM_TARGETS:
            e = E9[p]
            E = combine(*rows[p], w)
            for s in range(6):
                rng = np.random.default_rng(1000 * s + 3)
                idx = rng.integers(0, e.N, max(budgets))
                for b in budgets:
                    sub = idx[:b]
                    cells[b]["ret"].append(float(e.rmsd[sub[int(np.argmin(E[sub]))]]))
                    cells[b]["best"].append(float(e.rmsd[sub].min()))
        line = ""
        out["like_for_like"][str(w)] = {}
        for b in budgets:
            r = float(np.mean(cells[b]["ret"]))
            bb = float(np.mean(cells[b]["best"]))
            out["like_for_like"][str(w)][str(b)] = {"returned": r, "best_available": bb,
                                                    "gap": r - bb}
            line += f"{r:8.3f}{bb:9.3f}{r-bb:9.3f}"
        print(f"  {w:5.2f} " + line)
    print("\n  The coordinator's 126-target numbers for comparison (hamil, w=0.25):")
    print("    budget      10      300    3,000   20,000")
    print("    returned  3.742    3.565    3.632    3.572")
    print("    best      2.907    1.987    1.697    1.532")
    print("    gap       0.835    1.577    1.935    2.040")
    V.write("vqe_wsweep", out)
    print("\nwritten -> s14/results/vqe_wsweep.json")


if __name__ == "__main__":
    main()
