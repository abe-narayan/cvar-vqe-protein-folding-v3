#!/usr/bin/env python
"""s27/s28_D_c2_poolmember.py -- lane D: the pool-member control for lane C2's recognition audit.

C2 (S28-L35) asks which native-free scorer prefers a 0.29 A ORACLE structure to the production
average, with controls that are random signed combinations or Gaussian perturbations. Neither
control is a REAL PROTEIN TRACE. The natural zero-information-about-nativeness structure with
real local geometry is a pool member: the same scorer was computed on all 500 pool members by
S27 (`s27/cache/<pdb>.npz`, the same functions C2's adapters reproduce to 1e-9). For each
scorer this script reports, per target and averaged over 126:

  pct(X)   the fraction of the 500 pool members that score BETTER (lower) than structure X,
           plus half the ties: X's percentile in its own pool (0.5 = a typical member,
           1.0 = worse than every member);
  pref(random pool member vs PROD) = pct(PROD): how often a random real trace beats production;
  pref(circ_best vs PROD)          from C2's rows (ties 0.5), the number S28-L35 quotes;
  the paired contrast of the two (a scorer that "recognises" the ORACLE structure must prefer
  it to production MORE often than it prefers an arbitrary pool member to production);
  head-to-head: 1 - pct(circ_best) = how often the ORACLE structure beats a random pool member
  (production absent), beside NATIVE, RAND_SIGNED[0] and GAUSS_MATCHED[0].

ORACLE: every structure but PROD is chosen against the native (C2's labelling); this script
adds no deployable arm. Writes `s27/results/s28_D_c2_poolmember.json`.
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

ROWS = os.path.join(HERE, "results", "s28_C2_ca_rows.jsonl")
OUT = os.path.join(HERE, "results", "s28_D_c2_poolmember.json")
STRUCTS = ("PROD", "circ_best", "circ_s0", "NATIVE", "RAND_SIGNED[0]", "GAUSS_MATCHED[0]")
#: 2026-09-19 (resume): the same control on the C2 BUILT-CHAIN rows (`s28_C2_chain_rows.jsonl`),
#: whose scores live in two blocks: `scores` (the 16 backbone scorers on the projected chains) and
#: `scores_ca_on_chain` (the 15 CA scorers re-evaluated on them). `--rows`, `--scores-key`, `--out`
#: select them; the defaults reproduce the S28-L36 artefact unchanged. The pool channels in
#: `s27/cache` are on the members' REAL torsions / real CA traces, the rows' structures are
#: PROJECTED (ideal geometry): pct(X) then compares a projected structure against real ones, which
#: is the same mismatch C2's chain audit carries, stated beside every number.


def main(rows_path=ROWS, scores_key="scores", out_path=OUT):
    rows = [json.loads(l) for l in open(rows_path, encoding="utf-8") if l.strip()]
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    fail = np.array([p in set(I.FAIL18) for p in pdbs])
    scorers = list(rows[0][scores_key].keys())
    caches = {p: np.load(os.path.join(HERE, "cache", f"{p}.npz")) for p in pdbs}

    def score(r, sc, name):
        return float(r[scores_key][sc][r["names"].index(name)])

    out = {"n": len(rows), "scorers": {}, "structures": list(STRUCTS),
           "rows": os.path.relpath(rows_path, ROOT).replace(os.sep, "/"), "scores_key": scores_key}
    for sc in scorers:
        if sc not in caches[pdbs[0]].files:
            out["scorers"][sc] = {"note": "no pool channel of this name in s27/cache (pool-relative adapter)"}
            continue
        if any(np.isnan(score(r, sc, k)) for r in rows for k in STRUCTS):
            out["scorers"][sc] = {"note": "NaN score on some structure; skipped"}
            continue
        pct = {k: [] for k in STRUCTS}
        for r in rows:
            pool = np.asarray(caches[r["pdb"]][sc], float)
            for k in STRUCTS:
                v = score(r, sc, k)
                pct[k].append(float((pool < v).mean() + 0.5 * (pool == v).mean()))
        pct = {k: np.array(v) for k, v in pct.items()}
        pref_cb = np.array([1.0 if score(r, sc, "circ_best") < score(r, sc, "PROD")
                            else (0.5 if score(r, sc, "circ_best") == score(r, sc, "PROD") else 0.0) for r in rows])
        pref_s0 = np.array([1.0 if score(r, sc, "circ_s0") < score(r, sc, "PROD")
                            else (0.5 if score(r, sc, "circ_s0") == score(r, sc, "PROD") else 0.0) for r in rows])
        c_cb = ST.compare(pref_cb, pct["PROD"], folds, names=pdbs, label=f"{sc}: pref(circ_best vs PROD) - pref(pool member vs PROD)")
        c_s0 = ST.compare(pref_s0, pct["PROD"], folds, names=pdbs, label=f"{sc}: pref(circ_s0 vs PROD) - pref(pool member vs PROD)")
        res = {
            "pool_percentile_mean": {k: float(v.mean()) for k, v in pct.items()},
            "pool_percentile_median": {k: float(np.median(v)) for k, v in pct.items()},
            "pool_percentile_fail18": {k: float(v[fail].mean()) for k, v in pct.items()},
            "pool_percentile_other108": {k: float(v[~fail].mean()) for k, v in pct.items()},
            "pref_circ_best_vs_prod": float(pref_cb.mean()),
            "pref_circ_s0_vs_prod": float(pref_s0.mean()),
            "pref_pool_member_vs_prod": float(pct["PROD"].mean()),
            "h2h_circ_best_beats_pool_member": float(1 - pct["circ_best"].mean()),
            "h2h_native_beats_pool_member": float(1 - pct["NATIVE"].mean()),
            "h2h_rand_signed_beats_pool_member": float(1 - pct["RAND_SIGNED[0]"].mean()),
            "contrast_circ_best": {k: v for k, v in c_cb.items() if k != "concentration"},
            "contrast_circ_s0": {k: v for k, v in c_s0.items() if k != "concentration"},
            "fmt_circ_best": ST.fmt(c_cb),
            "fmt_circ_s0": ST.fmt(c_s0),
        }
        out["scorers"][sc] = res
        print(f"{sc:16s} pref(cb vs PROD) {pref_cb.mean():.3f}  pref(pool member vs PROD) {pct['PROD'].mean():.3f}  "
              f"diff {c_cb['effect']:+.4f} x{c_cb['effect_over_mde']:+.2f} fold[{c_cb['ci95_fold'][0]:+.3f},{c_cb['ci95_fold'][1]:+.3f}] {c_cb['folds_same_sign']}/5 | "
              f"h2h cb {1-pct['circ_best'].mean():.3f} nat {1-pct['NATIVE'].mean():.3f} rand {1-pct['RAND_SIGNED[0]'].mean():.3f} | "
              f"pct PROD {pct['PROD'].mean():.3f} (med {np.median(pct['PROD']):.2f}) cb {pct['circ_best'].mean():.3f} nat {pct['NATIVE'].mean():.3f}")
    ST.save_atomic(out_path, out, complete_keys=["n", "scorers", "structures"], module_file=__file__)
    print("wrote", out_path)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default=ROWS)
    ap.add_argument("--scores-key", default="scores", help="scores | scores_ca_on_chain")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    main(a.rows, a.scores_key, a.out)
