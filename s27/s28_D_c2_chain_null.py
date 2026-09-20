#!/usr/bin/env python
"""s27/s28_D_c2_chain_null.py -- lane D: an independent max-over-31 multiplicity null for lane
C2's built-chain recognition audit (S28-L34 caveat (e)).

Reads `s27/results/s28_C2_chain_rows.jsonl` (per target: `scores` = the 16 backbone scorers on
the projected chains; `scores_ca_on_chain` = the 15 CA scorers re-evaluated on them; `names`
gives the structure order; PROD, circ_s0, circ_best, NATIVE, RAND_SIGNED[0], ...). For every
scorer (31 = 16 + 15, the CA ones suffixed `@chain`):

  pref(X vs PROD) = mean over targets of [score(X) < score(PROD)] with ties at 0.5,
  for X in circ_best, circ_s0, NATIVE, RAND_SIGNED[0], with the tie count and the FAIL18 count
  (k of 18, S28-L34 (d));

and the max-over-31 sign-flip null for pref(circ_best): per draw, each target's indicator is
flipped (1 - ind, ties stay 0.5) with probability 1/2, the 31 scorer means are taken and their
maximum recorded; p_max = share of draws whose maximum reaches the observed best pref. Also the
same null for the number of scorers above 0.5 (how many "recognisers" chance gives).

ORACLE: every structure but PROD and the controls' directions is chosen against the native.
Writes `s27/results/s28_D_c2_chain_null.json`.
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

ROWS = os.path.join(HERE, "results", "s28_C2_chain_rows.jsonl")
OUT = os.path.join(HERE, "results", "s28_D_c2_chain_null.json")
STRUCTS = ("circ_best", "circ_s0", "NATIVE", "RAND_SIGNED[0]")


def main(rows_path=ROWS, out_path=OUT, n_null=4000, seed=20260919):
    rows = [json.loads(l) for l in open(rows_path, encoding="utf-8") if l.strip()]
    pdbs = [r["pdb"] for r in rows]
    n = len(pdbs)
    fail = np.array([p in set(I.FAIL18) for p in pdbs])
    scorers = [(nm, "scores") for nm in rows[0]["scores"]] + \
              [(nm + "@chain", "scores_ca_on_chain") for nm in rows[0]["scores_ca_on_chain"]]

    def ind(r, key, base, X):
        names = r["names"]
        a = float(r[key][base][names.index(X)]); b = float(r[key][base][names.index("PROD")])
        if np.isnan(a) or np.isnan(b):
            return np.nan
        return 1.0 if a < b else (0.5 if a == b else 0.0)

    out = {"n": n, "n_scorers": len(scorers), "pref": {}, "max_null": {}}
    IND = {}
    print(f"{n} targets, {len(scorers)} scorers")
    print(f"{'scorer':24s} {'pref(cb)':>8s} {'ties':>5s} {'FAIL18 k/18':>11s} {'pref(s0)':>8s} {'pref(NAT)':>9s} {'pref(RAND)':>10s}")
    for nm, key in scorers:
        base = nm.replace("@chain", "")
        d = {}
        for X in STRUCTS:
            v = np.array([ind(r, key, base, X) for r in rows])
            d[X] = v
        cb = d["circ_best"]
        ok = np.isfinite(cb)
        IND[nm] = cb
        out["pref"][nm] = {X: float(np.nanmean(d[X])) for X in STRUCTS}
        out["pref"][nm]["ties_circ_best"] = int((cb == 0.5).sum())
        out["pref"][nm]["fail18_k"] = float(np.nansum(cb[fail]))
        out["pref"][nm]["n_finite"] = int(ok.sum())
        print(f"{nm:24s} {out['pref'][nm]['circ_best']:8.3f} {out['pref'][nm]['ties_circ_best']:5d} {out['pref'][nm]['fail18_k']:9.1f}/18 "
              f"{out['pref'][nm]['circ_s0']:8.3f} {out['pref'][nm]['NATIVE']:9.3f} {out['pref'][nm]['RAND_SIGNED[0]']:10.3f}")
    # max-over-scorers sign-flip null
    names = [nm for nm, _ in scorers]
    D = np.column_stack([np.nan_to_num(IND[nm], nan=0.5) for nm in names])
    obs = D.mean(0)
    best = names[int(np.argmax(obs))]
    rng = np.random.default_rng(seed)
    mx = np.empty(n_null); n_above = np.empty(n_null)
    for t in range(n_null):
        eps = rng.random(n) < 0.5
        Dp = np.where(eps[:, None], D, 1 - D)
        m = Dp.mean(0)
        mx[t] = m.max(); n_above[t] = (m > 0.5).sum()
    out["max_null"] = dict(best_scorer=best, best_pref=float(obs.max()), n_null=n_null,
                           null_mean=float(mx.mean()), null_p95=float(np.percentile(mx, 95)),
                           p_max=float((mx >= obs.max()).mean()),
                           n_scorers_above_half_observed=int((obs > 0.5).sum()),
                           n_above_half_null_mean=float(n_above.mean()), n_above_half_null_p95=float(np.percentile(n_above, 95)))
    m = out["max_null"]
    print(f"\nmax-over-{len(names)} sign-flip null (n_null {n_null}): best {best} pref {m['best_pref']:.3f}; null mean {m['null_mean']:.3f} p95 {m['null_p95']:.3f}; p_max {m['p_max']:.3f}")
    print(f"scorers above 0.5: observed {m['n_scorers_above_half_observed']} vs null mean {m['n_above_half_null_mean']:.1f} p95 {m['n_above_half_null_p95']:.0f}")
    ST.save_atomic(out_path, out, complete_keys=["n", "pref", "max_null"], module_file=__file__)
    print("wrote", out_path)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default=ROWS)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--n-null", type=int, default=4000)
    a = ap.parse_args()
    main(a.rows, a.out, a.n_null)
