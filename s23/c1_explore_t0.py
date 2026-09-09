"""s23/c1_explore_t0.py -- DATED ADDENDUM, 2026-09-08, NOT PRE-REGISTERED.

`c1_probweight.py`'s pre-registered T in {0.5, 0.2, 0.1} all show entropy pinned near the
9-qubit maximum (~8.9-9.0 of 9 bits) at alpha=0.15 -- unlike the pipeline docstring's cited
alpha=1 mechanism (T=0.1 -> 0.076 bits, full collapse), the CVaR pressure at alpha=0.15 is
already weak enough that even T=0.1's entropy term keeps the tail close to uniform. That makes
every registered T-cell a WEAK test of the weighting hypothesis: if the trained distribution
never departs far from uniform, a probability-weighted average has little room to differ from
a uniform one regardless of whether training helps.

This file adds T=0.0 (pure CVaR, no entropy regulariser) as an EXPLORATORY, POST-HOC condition
to see the weighting effect where the trained distribution is genuinely non-uniform (s22 A2's own
finding: unregularised CVaR collapses hard). It is appended here, run AFTER the pre-registered
primary/secondary numbers were already on disk, and is reported as exploratory throughout --
never substituted for the registered primary. Same alpha, seeds, layers, iters as the
pre-registration; only T changes.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I         # noqa: E402
from s23 import qc_lib as QL            # noqa: E402
from s23 import c1_probweight as C1     # noqa: E402

T_EXPLORE = 0.0


def run(limit=None):
    tg = sorted(I.targets(), key=lambda r: r["pdb"])
    if limit:
        tg = tg[:limit]
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        u = I.load_univ(t["pdb"])
        nat_ca = np.asarray(u["nat_ca"], float)
        import s22.qcand_lib as QC
        pool = QC.build_pool(t["pdb"])
        enc = QC.Encoding(pool["score_dist"])
        W_real = pool["W"]
        seed_rows = []
        for seed in C1.SEEDS:
            tr = C1.one_arm(enc, C1.ALPHA, T_EXPLORE, seed, C1.ITERS, W_real)
            un = C1.one_arm(enc, C1.ALPHA, T_EXPLORE, seed, 0, W_real)
            seed_rows.append({
                "seed": seed, "m_trained": tr["m"], "m_untrained": un["m"],
                "gate1_trained": tr["gate1_pass"], "gate1_untrained": un["gate1_pass"],
                "entropy_bits_trained": tr["entropy_bits"],
                "entropy_bits_untrained": un["entropy_bits"],
                "rmsd_w_trained": I.ca_rmsd(tr["weighted_ca"], nat_ca),
                "rmsd_c_trained": I.ca_rmsd(tr["classical_ca"], nat_ca),
                "rmsd_w_untrained": I.ca_rmsd(un["weighted_ca"], nat_ca),
                "rmsd_c_untrained": I.ca_rmsd(un["classical_ca"], nat_ca),
            })
        rows.append({"pdb": t["pdb"], "n": t["n"], "fold": t["fold"],
                     "by_T": {str(T_EXPLORE): seed_rows}})
        if (c + 1) % 20 == 0:
            print(f"  {c + 1}/{len(tg)}  {time.time() - t0:.1f}s", flush=True)
    out = {"rows": rows, "complete": len(rows) == 126, "n_expected": 126,
          "config": {"alpha": C1.ALPHA, "T_primary": T_EXPLORE, "T_secondary": [],
                     "seeds": list(C1.SEEDS), "layers": C1.LAYERS, "iters": C1.ITERS,
                     "lr": C1.LR, "n_qubits": C1.N_QUBITS},
          "EXPLORATORY": True, "dated": "2026-09-08", "pre_registered": False}
    QL.save_json("c1_explore_t0_raw.json", out)
    print(f"done: {len(rows)}/126 complete={out['complete']} {time.time() - t0:.1f}s")
    return out


if __name__ == "__main__":
    run()
