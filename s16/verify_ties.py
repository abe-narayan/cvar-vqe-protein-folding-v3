"""s16/verify_ties.py -- VERIFY workstream.  Does `np.argsort` in `s16/integrate._readouts`
leak the enumeration order on ties?

The programme has been burned once by `np.argmin` on a tied signal reading the ORACLE sort
order and inventing a 1.386 A winner, so the BRIEF makes this a standing check.

THE SPECIFIC RISK HERE.  `_pick` returns `idx` sorted by the restraint objective `en.prior`,
and the objective is correlated with RMSD.  `_readouts` then calls
`np.argsort(legacy, kind="stable")[:keep]` and the same on `amber`.  A stable argsort breaks
ties by POSITION IN `idx`, i.e. by objective rank, i.e. toward better RMSD.  So any tie in the
Legacy or AMBER column silently converts objective skill into apparent filter skill, and the
random control `rnd` (a permutation) gets none of it.

This module censuses the ties without re-running the VQE, by regenerating the `uniform`
generator's ensemble -- which is `SD.stable_rng(pdb, f"unif{seed}")` and therefore reproduces
the integrate run's ensemble bit for bit -- and computing the same Legacy and AMBER columns.

Deliberately small: three targets, one seed, m = 75.  The question is whether ties exist at
all, and if they do at what rate; that does not need the full grid.

Writes `s16/results/verify_ties.json`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
RESULTS = os.path.join(HERE, "results")

from s14 import vqe_lib as V                 # noqa: E402
from s16 import energy_lib as EL             # noqa: E402
from s16.integrate import BUDGET, _pick, _uniform   # noqa: E402


def main(pdbs=("1CS9", "6F3V", "9UV5"), seed=0, m=75):
    from s13 import qarch_lib as QL
    out = {}
    print("=" * 92)
    print("TIE CENSUS in s16/integrate._readouts   (uniform generator, seed 0, m = 75)")
    print("=" * 92)
    print(f"  {'pdb':<8}{'m':>4}{'prior ties':>12}{'Legacy ties':>13}{'AMBER ties':>12}"
          f"{'rho(idx rank, RMSD)':>22}")
    for pdb in pdbs:
        en = V.Enum(pdb)
        sp = QL.Space(pdb, en.k, seq=en.seq, n=en.n, fold=en.fold)
        seen = _uniform(en, BUDGET, seed)
        idx = _pick(en, seen, m)
        S = en.states(idx)
        amber = np.asarray(QL.amber_energies(sp, S), float)
        comp = EL.legacy_components_of_windows(
            en.seq, en.PHI[np.arange(en.n)[None, :], S], en.PSI[np.arange(en.n)[None, :], S])
        legacy = np.asarray(EL.legacy_total_from(comp), float)
        pr = np.asarray(en.prior[idx], float)
        rm = np.asarray(en.rmsd[idx], float)            # ORACLE, for the leak diagnostic only

        def nties(x):
            return int(len(x) - len(np.unique(x)))

        rho = float(np.corrcoef(np.argsort(np.argsort(np.arange(len(idx)))),
                                np.argsort(np.argsort(rm)))[0, 1])
        out[pdb] = {"m": int(len(idx)), "prior_ties": nties(pr),
                    "legacy_ties": nties(legacy), "amber_ties": nties(amber),
                    "rho_idxrank_rmsd": rho}
        print(f"  {pdb:<8}{len(idx):>4}{nties(pr):>12}{nties(legacy):>13}{nties(amber):>12}"
              f"{rho:>+22.3f}")
    tot = sum(v["legacy_ties"] + v["amber_ties"] for v in out.values())
    print(f"\n  total Legacy + AMBER ties across the three targets: {tot}")
    if tot == 0:
        print("  VERDICT: no ties, so the stable-argsort tie-break cannot leak the objective")
        print("  order into `leg` or `amb`.  The tie trap does NOT apply to this result.")
    else:
        print("  VERDICT: ties exist -- the stable argsort breaks them toward better objective")
        print("  rank, which is correlated with RMSD.  The effect must be re-measured with the")
        print("  outcome averaged over the tied argsort set.")
    print("\n  `rho(idx rank, RMSD)` is the size of the channel a tie would leak into: it is the")
    print("  rank correlation between position in `idx` (objective order) and true RMSD.")
    with open(os.path.join(RESULTS, "verify_ties.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote s16/results/verify_ties.json")


if __name__ == "__main__":
    main()
