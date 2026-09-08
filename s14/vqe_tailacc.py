"""SPRINT 14 / VQE -- TAIL-RESTRICTED pairwise accuracy: the statistic nobody reported.

Three proxies for "objective usefulness" have been used this sprint and none of them is one:

  * **global rank correlation** -- dominated by the bulk.
  * **in-decile rank correlation** -- NON-MONOTONE in quality (`s14/vqe_decile.py`): it
    measures how much orderable RMSD spread is left inside the decile, which FALLS as the
    objective improves.
  * **bulk pairwise accuracy** -- measured over random pairs, so it is dominated by the 99%
    of the space a search never visits.  The distogram wins it decisively (0.654 at a 0.25 A
    gap) and has the WORST certified argmin (3.237).

The quantity that actually governs an argmin or a small-m readout is **pairwise accuracy
INSIDE THE REGION THE SEARCH ENDS UP IN**.  That is what this computes: restrict to the
objective's own lowest `frac` of the space, then measure ordering accuracy among pairs whose
true quality differs by at least `gap`.

It is also a direct test of the mechanism that has now been confirmed three ways --
"the distogram orders the bulk, the torsion prior places the optimum".  If that is right:

    at w=0 (pure prior)      tail accuracy should be HIGH and bulk accuracy near chance
    at w=1 (pure distogram)  bulk accuracy should be HIGH and tail accuracy should COLLAPSE

and the crossing between the two is the mechanism made quantitative on a single statistic.

    python -m s14.vqe_tailacc
"""
from __future__ import annotations

import numpy as np

from s14 import vqe_lib as V
from s14.vqe_hamil import tabulate, combine

FRACS = (1.0, 0.1, 0.01, 0.001, 0.0001)          # 1.0 = the bulk statistic
GAPS = (0.0, 0.25, 0.5, 1.0)
WS = (0.0, 0.25, 0.5, 1.0)
NPAIR = 200_000


def pair_accuracy(E, rmsd, idx, gap, rng, npair=NPAIR):
    """Fraction of quality-separated pairs the objective orders correctly.

    Ties in the OBJECTIVE are counted as half, which is the only unbiased convention: a
    tied objective carries no information and scoring it as a win or a loss would import
    the index order (the recorded tie-breaking trap that once invented a 1.386 A winner).
    """
    if idx.size < 2:
        return np.nan, 0
    a = idx[rng.integers(0, idx.size, npair)]
    b = idx[rng.integers(0, idx.size, npair)]
    ok = np.abs(rmsd[a] - rmsd[b]) > gap
    if ok.sum() < 100:
        return np.nan, int(ok.sum())
    a, b = a[ok], b[ok]
    better = rmsd[a] < rmsd[b]                    # a is truly better
    de = E[a] - E[b]
    correct = np.where(de == 0.0, 0.5, (de < 0) == better)
    return float(np.mean(correct)), int(len(a))


def main():
    rows = {p: tabulate(p) for p in V.ENUM_TARGETS}
    E9 = {p: V.Enum(p) for p in V.ENUM_TARGETS}
    rng = np.random.default_rng(0)

    print("=" * 104)
    print("TAIL-RESTRICTED PAIRWISE ACCURACY, nine fully enumerated targets")
    print("  restrict to the objective's own lowest `frac` of the space, then measure")
    print("  ordering accuracy among pairs separated by at least `gap` in true CA-RMSD.")
    print("  frac = 1.0 is the bulk statistic that has been reported so far.")
    print("=" * 104)
    out = {}
    for w in WS:
        lab = {0.0: "hamil w=0 (pure prior)", 1.0: "hamil w=1 (pure distogram)"}.get(
            w, f"hamil w={w}")
        print(f"\n--- {lab} ---")
        print(f"  {'frac of space':>14s} {'n configs':>10s} "
              + "".join(f"{f'gap>{g}':>10s}" for g in GAPS))
        for frac in FRACS:
            accs = {g: [] for g in GAPS}
            npc = []
            for p in V.ENUM_TARGETS:
                e = E9[p]
                E = combine(*rows[p], w)
                k = max(2, int(round(frac * e.N)))
                idx = np.argpartition(E, k - 1)[:k] if k < e.N else np.arange(e.N)
                npc.append(k)
                for g in GAPS:
                    a, _ = pair_accuracy(E, e.rmsd, idx, g, rng)
                    if np.isfinite(a):
                        accs[g].append(a)
            out.setdefault(str(w), {})[str(frac)] = {
                str(g): (float(np.mean(accs[g])) if accs[g] else None) for g in GAPS}
            out[str(w)][str(frac)]["n_configs"] = int(np.mean(npc))
            print(f"  {frac:14.4g} {int(np.mean(npc)):10d} "
                  + "".join(f"{np.mean(accs[g]):10.3f}" if accs[g] else f"{'-':>10s}"
                            for g in GAPS))

    print()
    print("=" * 104)
    print("THE MECHANISM ON ONE STATISTIC: bulk accuracy vs tail accuracy, by w")
    print("=" * 104)
    print(f"  {'w':>5s} {'BULK (frac=1)':>14s} {'TAIL (frac=0.001)':>18s} "
          f"{'tail - bulk':>12s} {'certified argmin':>17s}")
    for w in WS:
        b = out[str(w)]["1.0"]["0.25"]
        t = out[str(w)]["0.001"]["0.25"]
        am = np.mean([E9[p].rmsd[int(np.argmin(combine(*rows[p], w)))]
                      for p in V.ENUM_TARGETS])
        print(f"  {w:5.2f} {b:14.3f} {t:18.3f} {t-b:+12.3f} {am:17.3f}")
    print("\n  (gap > 0.25 A, the ENER threshold.  ENER's floor: no PHYSICAL objective")
    print("  exceeds 0.511 below a 0.25 A gap.)")
    V.write("vqe_tailacc", out)
    print("\nwritten -> s14/results/vqe_tailacc.json")


if __name__ == "__main__":
    main()
