"""SPRINT 14, OBJ -- what ranking quality would actually be needed? The transfer curve.

Sections 1 and 2 measure how well objectives rank.  This module measures how well an
objective WOULD HAVE TO rank to matter, by sweeping a noisy oracle over the enumerated
space:

    e_sigma(x) = rmsd(x) + sigma * N(0, 1)

For each sigma we report, on the identical evaluation sets used everywhere else:
  * in-band pairwise ordering accuracy below 1.5 A (the metric section 1.2 established),
  * mean CA-RMSD of the low-energy decile and of the top-100,
  * the argmin.

That gives a calibration curve `in-band accuracy -> achievable RMSD`, so any measured
objective (Legacy 0.516-0.543, AMBER 0.479-0.501, the learned objective) can be read
against the accuracy the 2.0 A target demands.

ORACLE DIAGNOSTIC throughout -- e_sigma reads the native.  It is a requirement gauge, not
a method.

    python -m s14.obj_noisy
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                    # noqa: E402
from s13.qarch_lib import spearman                 # noqa: E402
from s14 import obj_enum as E                      # noqa: E402
from s14 import obj_train as T                     # noqa: E402
from s14 import obj_exp as X                       # noqa: E402
from s14 import obj_floor as F                     # noqa: E402

SIGMAS = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 1e6]


def run(targets, seed=0):
    rows = []
    t0 = time.time()
    for sg in SIGMAS:
        per = []
        for p in targets:
            en = T.Enum(p)
            idx = X.eval_idx(en)
            r = en.rmsd[idx].astype(np.float64)
            rng = np.random.default_rng(hash((p, sg)) % (1 << 31))
            e = r + sg * rng.standard_normal(len(r))
            o = np.argsort(e, kind="mergesort")
            dec = o[:len(e) // 10]
            pi = F.pair_discrimination(e, r, rng, npair=300_000, sub=dec)
            per.append(dict(
                inband_lt15=float(np.nanmean([pi[k][0] for k in
                                              ("0.25-0.5", "0.5-1.0", "1.0-1.5")
                                              if k in pi])),
                inband_05=float(pi.get("0.25-0.5", (np.nan,))[0]),
                inband_1=float(pi.get("0.5-1.0", (np.nan,))[0]),
                rho_decile=spearman(e[dec], r[dec]),
                dec_mean=float(r[dec].mean()),
                top100=float(r[o[:100]].mean()),
                argmin=float(r[o[0]]),
                space_mean=float(r.mean()), space_min=float(r.min())))
        agg = {k: float(np.nanmean([q[k] for q in per])) for k in per[0]}
        agg["sigma"] = sg
        rows.append(agg)
        print(f"sigma={sg:8.2f}  inband<1.5={agg['inband_lt15']:.3f} "
              f"(0.25-0.5:{agg['inband_05']:.3f})  rho_dec={agg['rho_decile']:+.3f}  "
              f"decile={agg['dec_mean']:.3f}  top100={agg['top100']:.3f}  "
              f"argmin={agg['argmin']:.3f}   [{time.time()-t0:.0f}s]", flush=True)
    print(f"\nspace mean {rows[0]['space_mean']:.3f} A, space min "
          f"{rows[0]['space_min']:.3f} A over {len(targets)} targets")
    I.write("s14_obj_noisy", {"targets": targets, "rows": rows},
            n_expected=len(SIGMAS))
    return rows


if __name__ == "__main__":
    tg = [a for a in sys.argv[1:] if not a.startswith("-")]
    run(tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])])
