"""TORSION-PREDICTOR -- the LIBRARY-PRIOR null (the strongest null available).

`s13/coord_objval.prior_logp(seq, k)` is a genuine, deployable, leakage-safe 1-local torsion
prior: `torsion_lib2.library_for` gives k states per residue CLASS (GENERAL / GLY / PRO /
PRE_PRO) clustered over the identity-held-out database, and `prior_logp` recovers each
state's occupancy by assigning every held-out observation of that class to its nearest state.
A residue's class is known from its sequence alone, so this contains **no learned prediction
whatsoever** -- and it is a much stronger null than shuffled labels or composition, because it
already carries the amino-acid-class Ramachandran structure that most of a weak sequence
predictor's apparent skill consists of.

Converted here onto the same 18x18 grid every other arm is scored on, by placing a von Mises
kernel of concentration matched to the 20-degree cell width at each library state and
weighting by its occupancy.  Without the kernel the null would be 32 delta spikes against the
trained models' smooth densities, which would flatter the trained models on every density
score.

    python -m s13.tors_nulls
"""
from __future__ import annotations
import os, sys, time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s13 import tors_common as T          # noqa: E402
from s13.tors_train import post_path      # noqa: E402
from s12 import instrument as I           # noqa: E402

K_LIB = 32
KAPPA = 1.0 / (20.0 * T.D2R) ** 2         # vM concentration matched to the 20-degree cell


def libprior_posterior(seq, k=K_LIB):
    from s13 import coord_objval as CO
    import torsion_lib2 as tl2
    tab = np.asarray(tl2.library_for(seq, k, exclude_seq=seq), float)     # (n, k, 2)
    w = np.exp(CO.prior_logp(seq, k, exclude_seq=seq))                    # (n, k)
    n = len(seq)
    lp = (KAPPA * np.cos(T.GRID_PHI[None, None, :] - tab[:, :, 0:1])
          + KAPPA * np.cos(T.GRID_PSI[None, None, :] - tab[:, :, 1:2]))   # (n, k, G)
    lp -= lp.max(2, keepdims=True)
    comp = np.exp(lp)
    comp /= comp.sum(2, keepdims=True)
    P = (w[:, :, None] * comp).sum(1)
    return (P / P.sum(1, keepdims=True)).astype(np.float32)


def main():
    tg = I.targets()
    out = {}
    t0 = time.time()
    for k, t in enumerate(tg):
        out[t["pdb"]] = libprior_posterior(t["seq"])
        if k % 25 == 0:
            print(f"  {k+1}/{len(tg)} {t['pdb']} [{time.time()-t0:.0f}s] "
                  f"free={I.free_gb():.2f}", flush=True)
    np.savez_compressed(post_path("n_libprior"), **out)
    print("wrote", post_path("n_libprior"))


if __name__ == "__main__":
    main()
