#!/usr/bin/env python
"""s29/s29_M_F2_gate.py -- F2's reproduction gate, run under S12's OWN scoring convention.

Stage 0 (`s29_M_F2_supply.py`) scores with the SHIPPED per-pair weight 1/(sd+0.5) and got
PROD 3.0624 / ORACLE 2.4254 against `s12/obj_FINDINGS.md`'s 3.078 / 2.402 -- close, not equal.
The cause is identified, not guessed: `s12/obj_common.py:93-97 score_l1` is called from
`s12/obj_profile.py:156 evaluate` with `w=None`, i.e. S12's profile arms are scored with a
UNIFORM per-pair weight. This script re-runs the two gate arms with that convention and nothing
else changed, so the gate either reproduces or the discrepancy is real.

    python s26/jobrun.py --agent S29M --tag CPU --name m_f2_gate --est-ram 0.8 -- \
        python s29/s29_M_F2_gate.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s24 import d_harness as H          # noqa: E402
from s27 import run_pool as RP          # noqa: E402
from s29 import s29_M_F2_supply as F2      # noqa: E402


def main():
    t0 = time.time()
    out = {"weighted": {}, "uniform": {}}
    rows = []
    for k, pdb in enumerate(P.targets()):
        d = F2.target_data(pdb)
        key = RP.rng_for(pdb, "tiekey").random(d["cand"].k)
        kk = np.searchsorted(d["us"], d["sep"])
        base = d["exp"] - d["q_disto"][kk]
        r = {"pdb": pdb, "fold": d["fold"]}
        for nm, q in (("PROD", d["q_disto"]), ("ORACLE_PROF", d["q_true"])):
            tgt = q[kk] + base
            res = np.abs(d["D"] - tgt[None, :])
            for conv, sc in (("weighted", (res * d["w"][None, :]).mean(1)),
                             ("uniform", res.mean(1))):
                top = np.lexsort((key, RP.zr(sc)))[:RP.M]
                C, _ = H.readout_uniform(d["cand"], top)
                r["%s_%s" % (nm, conv)] = float(I.ca_rmsd(C, d["nat"]))
        rows.append(r)
        if (k + 1) % 30 == 0:
            print("  %d/126 (%.1f min)" % (k + 1, (time.time() - t0) / 60), flush=True)
    for conv in ("weighted", "uniform"):
        for nm in ("PROD", "ORACLE_PROF"):
            out[conv][nm] = float(np.mean([r["%s_%s" % (nm, conv)] for r in rows]))
        out[conv]["gap"] = out[conv]["PROD"] - out[conv]["ORACLE_PROF"]
    out["s12_reference"] = {"PROD_plain_L1_vs_Ed": 3.078, "ORACLE_true_profile": 2.402,
                            "gap": 0.676, "source": "s12/obj_FINDINGS.md:306-331 (section 6)"}
    out["n"] = len(rows)
    out["secs"] = time.time() - t0
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    for conv in ("weighted", "uniform"):
        print("%-9s PROD %.4f  ORACLE %.4f  gap %.4f   | S12: 3.078 / 2.402 / 0.676  -> dPROD %+.4f dORACLE %+.4f"
              % (conv, out[conv]["PROD"], out[conv]["ORACLE_PROF"], out[conv]["gap"],
                 out[conv]["PROD"] - 3.078, out[conv]["ORACLE_PROF"] - 2.402))
    ST.save_atomic(os.path.join(HERE, "results", "s29_M_F2_gate.json"),
                   dict(out, rows=rows), module_file=__file__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
