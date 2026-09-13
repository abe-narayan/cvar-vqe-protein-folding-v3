#!/usr/bin/env python
"""s26/w_tiebreak.py -- the production endpoint's noise floor from BLOSUM tie-breaking at the
K = 500 pool boundary.  Lane W, Sprint 26.  Pre-registration: s26/PREREG_tiebreak_floor.md.

Every target's pool boundary falls inside a BLOSUM tie class (median 115 windows share the
boundary score, ~57 inside the pool), broken by corpus order.  `draws` re-draws the boundary
members uniformly at random (8 seeded draws), runs the production emission on each and stores
the emissions and the native-free triangle bounds; `endpoint` (GATED) reads natives and reports
s_tie, m_tie and the paired MDE between draws.

    python s26/jobrun.py --agent W --tag CPU --name w_tiebreak_probe --est-ram 0.5 -- python s26/w_tiebreak.py draws --probe 1CEK
    python s26/jobrun.py --agent W --tag CPU --name w_tiebreak_draws --est-ram 0.5 -- python s26/w_tiebreak.py draws
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
import w_selfcopy as W                               # noqa: E402  (imported, never edited here)

N_DRAWS = 8
RES = W.RES
K, TOPM = W.K, W.TOPM


def boundary_class(sim, order, k=K):
    """(non-tied prefix of the pool, tie-class members in corpus order, n_in)."""
    sim = np.asarray(sim, float); order = np.asarray(order, int)
    s = sim[order[k - 1]]
    pool = order[:k]
    prefix = pool[sim[pool] > s]
    tied = np.flatnonzero(sim == s)                       # corpus order
    n_in = int(k - len(prefix))
    assert (sim[pool] >= s).all() and n_in == int((sim[pool] == s).sum())
    return prefix, tied, n_in


def draw_pool(sim, order, seed_parts, k=K):
    """One random tie-break: the non-tied prefix, then a uniform subset of the tie class of the
    production size, in corpus order.  Seeded by s15.seed.stable_rng."""
    from s15 import seed as SD
    prefix, tied, n_in = boundary_class(sim, order, k)
    rng = SD.stable_rng(*seed_parts)
    chosen = np.sort(rng.choice(tied, size=n_in, replace=False))
    pool = np.concatenate([prefix, chosen]).astype(int)
    assert len(pool) == k and len(set(pool.tolist())) == k
    return pool, prefix, tied, n_in


def draws_target(t, n_draws=N_DRAWS, verbose=True):
    pdb, seq, fold = t["pdb"], t["seq"], t["fold"]
    u = W.load_blind(pdb)
    Wu = u["W"]; order = np.asarray(u["order"], int); sim = np.asarray(u["sim"], float)
    dg = I.distogram(pdb)
    pool0 = order[:K]
    e0 = W.emit(Wu[pool0], dg, seq, fold)
    gate = W.production_gate(pdb, e0)
    prefix, tied, n_in = boundary_class(sim, order)
    row = {"pdb": pdb, "n": t["n"], "fold": fold, "gate": gate, "n_tied": int(len(tied)), "n_in": n_in,
           "production": {k: e0[k] for k in ("cloud", "chain", "fit")},
           "production_argmin_windows": Wu[pool0][e0["argmin_set"]], "draws": []}
    top0 = set(pool0[e0["top"]].tolist())
    for k in range(n_draws):
        pool, _, _, _ = draw_pool(sim, order, ("w_tiebreak", pdb, k))
        e = W.emit(Wu[pool], dg, seq, fold)
        same = bool(set(pool0[e0["argmin_set"]].tolist()) == set(pool[e["argmin_set"]].tolist()))
        row["draws"].append({"k": k, "pool_overlap": float(len(set(pool.tolist()) & set(pool0.tolist())) / K),
                             "top75_replaced": int(TOPM - len(top0 & set(pool[e["top"]].tolist()))),
                             "argmin_same": same, **W.triangle(e0, e, Wu[pool0], Wu[pool], same),
                             "emissions": {q: e[q] for q in ("cloud", "chain", "fit")},
                             "argmin_windows": Wu[pool][e["argmin_set"]]})
    assert W.finite(row), pdb
    if verbose:
        print("  %s f%d tied %d (in %d) gate=%s | top75 replaced %s | tri_arm %s"
              % (pdb, fold, len(tied), n_in, gate["top75_equals_sub"],
                 " ".join(str(d["top75_replaced"]) for d in row["draws"]),
                 " ".join("%.3f" % d["tri_arm"] for d in row["draws"])), flush=True)
    del u
    return row


def draws(probe=None, verbose=True):
    tg = I.targets()
    if probe:
        row = draws_target([t for t in tg if t["pdb"] == probe][0], verbose=verbose)
        p = W.save("tiebreak_probe_%s" % probe, {"label": "tie-break probe", "rows": [row]})
        print("  ->", p); return [row]
    name = "tiebreak_draws"
    prev = W.load_result(name)
    rows = prev["rows"] if prev else []
    done = {r["pdb"] for r in rows}
    t_last = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(draws_target(t, verbose=verbose))
        if time.time() - t_last > 300 or len(rows) % 10 == 0:
            W.save(name, {"label": "tie-break draws, native-free", "rows": rows}, rows=rows, n_expected=len(tg))
            t_last = time.time()
    tri = np.array([[d["tri_arm"] for d in r["draws"]] for r in rows])
    summ = {"n": len(rows), "top75_equals_sub_count": int(sum(r["gate"]["top75_equals_sub"] for r in rows)),
            "top75_replaced_mean": float(np.mean([d["top75_replaced"] for r in rows for d in r["draws"]])),
            "tri_arm_median_over_targets_of_mean_over_draws": float(np.median(tri.mean(1))),
            "tri_arm_p90": float(np.percentile(tri.mean(1), 90)),
            "argmin_same_frac": float(np.mean([d["argmin_same"] for r in rows for d in r["draws"]]))}
    p = W.save(name, {"label": "tie-break draws, native-free", "summary": summ, "rows": rows}, rows=rows,
               n_expected=len(tg), complete_keys=("pdb", "draws", "production"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


def endpoint(verbose=True):
    W.require_signoff("tiebreak endpoint")
    Z = W.load_result("tiebreak_draws")
    if not (Z and Z.get("complete")):
        raise SystemExit("run `draws` to completion first")
    rows = []
    for r in Z["rows"]:
        N = np.asarray(I.load_univ(r["pdb"])["nat_ca"], float)
        prod = W.rmsd_to_native(r["production"], N); prod["sel"] = W.sel_rmsd(r["production_argmin_windows"], N)
        dr = []
        for d in r["draws"]:
            e = W.rmsd_to_native(d["emissions"], N); e["sel"] = W.sel_rmsd(d["argmin_windows"], N); dr.append(e)
        rows.append({"pdb": r["pdb"], "fold": r["fold"], "production": prod, "draws": dr,
                     "s_tie": {b: float(np.std([e[b] for e in dr], ddof=1)) for b in prod}})
    pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    out = {"label": "GATED tie-break floor", "rows": rows, "stats": {}}
    for b in ("arm", "cloud", "sel", "fit"):
        M = np.array([[e[b] for e in r["draws"]] for r in rows])            # (126, n_draws)
        prod = np.array([r["production"][b] for r in rows])
        means = M.mean(0)
        pair_mde = []; pair_eff = []
        for a in range(M.shape[1]):
            for c in range(a + 1, M.shape[1]):
                res = ST.compare(M[:, a], M[:, c], folds, names=pdbs, label="draw %d vs draw %d [%s]" % (a, c, b))
                pair_mde.append(res["mde"]); pair_eff.append(res["effect"])
        vs_prod = ST.compare(prod, M.mean(1), folds, names=pdbs, label="production convention vs mean over draws [%s]" % b)
        out["stats"][b] = {"m_tie_sd_of_mean_over_draws": float(means.std(ddof=1)), "mean_over_draws": float(means.mean()),
                           "production_mean": float(prod.mean()),
                           "s_tie_median": float(np.median([r["s_tie"][b] for r in rows])),
                           "s_tie_p90": float(np.percentile([r["s_tie"][b] for r in rows], 90)),
                           "paired_mde_between_draws_median": float(np.median(pair_mde)),
                           "paired_mde_between_draws_max": float(np.max(pair_mde)),
                           "paired_effect_between_draws_maxabs": float(np.max(np.abs(pair_eff))),
                           "production_vs_draws": vs_prod}
        if verbose:
            s = out["stats"][b]
            print("  [%s] m_tie %.4f  mean over draws %.4f  production %.4f  s_tie median %.4f p90 %.4f  paired MDE median %.4f max %.4f"
                  % (b, s["m_tie_sd_of_mean_over_draws"], s["mean_over_draws"], s["production_mean"], s["s_tie_median"], s["s_tie_p90"],
                     s["paired_mde_between_draws_median"], s["paired_mde_between_draws_max"]))
            print(ST.fmt(vs_prod))
    p = W.save("tiebreak_endpoint", out)
    print("  ->", p)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("draws", "endpoint"))
    ap.add_argument("--probe", default=None)
    a = ap.parse_args(argv)
    if a.cmd == "draws":
        draws(probe=a.probe)
    else:
        endpoint()
    return 0


if __name__ == "__main__":
    sys.exit(main())
