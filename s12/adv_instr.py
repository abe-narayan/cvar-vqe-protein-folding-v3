"""ADVERSARIAL AUDIT 0 -- the instrument itself.

Every agent's numbers flow through s12/instrument.py.  A defect there contaminates all of
them at once, so this runs first.  Checks, each against the PRODUCTION artefact
(bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json) rather than against another helper:

  A. shipped_score top-75 == production `sub`?  (float32 risk / float32 mean in the
     instrument vs float64 in core.pipeline.score -> s7.debias.score_risk)
  B. coordinate_average(W[pool][sub]) == production `avg_ca`?
  C. project(C, seq, fold) == production `fit_ca` / `ca`?   (sampled: ~5 s per call)
  D. cached s12/cache/disto_<pdb>.npz really the leave-fold-out model for fold(pdb)?
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I


def score_f64(dg, D):
    """What production does: risk upcast to float64, gather, float64 mean."""
    grid = np.asarray(dg["grid"], float)
    risk = np.asarray(dg["risk"], float)
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


def main(limit=None):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    rows = []
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        rec = I.shipped_record(pdb)
        sub = np.asarray(rec["sub"], int)
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(t["n"])
        W64 = u["W"][p]                      # float64 upcast of the float32 npz, as prod
        D = I.pair_dists(W64, i, j)
        sc32 = I.shipped_score(dg, D)                       # instrument path
        sc32r = I.shipped_score(dg, D.astype(np.float32).astype(float))  # selfcheck path
        sc64 = score_f64(dg, D)                             # production path

        o32 = np.argsort(sc32, kind="stable")[:I.M]
        o32r = np.argsort(sc32r, kind="stable")[:I.M]
        o64 = np.argsort(sc64, kind="stable")[:I.M]

        rr = u["rr"][p]
        # A: membership agreement with production
        row = dict(pdb=pdb, fold=t["fold"], n=t["n"],
                   sub_eq_f64=int(np.array_equal(np.sort(o64), np.sort(sub))),
                   sub_eq_f32=int(np.array_equal(np.sort(o32), np.sort(sub))),
                   sub_eq_f32r=int(np.array_equal(np.sort(o32r), np.sort(sub))),
                   n_diff_f32_vs_f64=int(I.M - len(set(o32.tolist()) & set(o64.tolist()))),
                   n_diff_f32r_vs_f64=int(I.M - len(set(o32r.tolist()) & set(o64.tolist()))),
                   argmin_f32=int(np.argmin(sc32)), argmin_f64=int(np.argmin(sc64)),
                   rr_argmin_f32=float(rr[int(np.argmin(sc32))]),
                   rr_argmin_f64=float(rr[int(np.argmin(sc64))]),
                   rr_top75_best_prod=float(rr[sub].min()),
                   rr_top75_best_f64=float(rr[o64].min()),
                   rr_top75_best_f32=float(rr[o32].min()),
                   max_abs_score_diff=float(np.max(np.abs(sc32 - sc64))),
                   n_exact_ties_f32=int(len(sc32) - len(np.unique(sc32))),
                   n_exact_ties_f64=int(len(sc64) - len(np.unique(sc64))))

        # B: coordinate average vs production avg_ca
        C, b = I.coordinate_average(W64[sub])
        avg = np.asarray(rec["avg_ca"], float)
        row["avg_rmsd_to_prod"] = float(I.ca_rmsd(C, avg))
        row["avg_maxabs_to_prod"] = float(np.abs(C - avg).max())
        # D: distogram cache sanity
        row["risk_shape"] = list(np.asarray(dg["risk"]).shape)
        row["npairs_ok"] = int(np.asarray(dg["risk"]).shape[0] == len(i))
        rows.append(row)
        print(json.dumps({k: row[k] for k in
                          ("pdb", "sub_eq_f64", "sub_eq_f32", "sub_eq_f32r",
                           "n_diff_f32_vs_f64", "avg_rmsd_to_prod", "max_abs_score_diff")}),
              flush=True)

    agg = dict(
        n=len(rows),
        n_sub_eq_f64=int(sum(r["sub_eq_f64"] for r in rows)),
        n_sub_eq_f32=int(sum(r["sub_eq_f32"] for r in rows)),
        n_sub_eq_f32r=int(sum(r["sub_eq_f32r"] for r in rows)),
        mean_n_diff_f32_vs_f64=float(np.mean([r["n_diff_f32_vs_f64"] for r in rows])),
        max_n_diff_f32_vs_f64=int(max(r["n_diff_f32_vs_f64"] for r in rows)),
        n_argmin_moved=int(sum(r["argmin_f32"] != r["argmin_f64"] for r in rows)),
        shipped_f32=float(np.mean([r["rr_argmin_f32"] for r in rows])),
        shipped_f64=float(np.mean([r["rr_argmin_f64"] for r in rows])),
        top75_best_prod=float(np.mean([r["rr_top75_best_prod"] for r in rows])),
        top75_best_f64=float(np.mean([r["rr_top75_best_f64"] for r in rows])),
        top75_best_f32=float(np.mean([r["rr_top75_best_f32"] for r in rows])),
        max_avg_rmsd_to_prod=float(max(r["avg_rmsd_to_prod"] for r in rows)),
        mean_avg_rmsd_to_prod=float(np.mean([r["avg_rmsd_to_prod"] for r in rows])),
        n_npairs_bad=int(sum(1 - r["npairs_ok"] for r in rows)),
        mean_exact_ties_f32=float(np.mean([r["n_exact_ties_f32"] for r in rows])),
        mean_exact_ties_f64=float(np.mean([r["n_exact_ties_f64"] for r in rows])),
    )
    print(json.dumps(agg, indent=1))
    I.write("adv_instrument", dict(agg=agg, rows=rows))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
