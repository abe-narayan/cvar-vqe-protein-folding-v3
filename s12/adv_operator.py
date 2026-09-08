"""ADVERSARIAL AUDIT 1c -- is "generation is wasted" about the OBJECTIVE or about the OPERATOR?

The sprint's central claim: every generation-side improvement is wasted because the
objective cannot rank within the improved set.  The alternative it never controlled for:
the terminal coordinate-average operator returns approximately a fixed quantile of its
input set's RMSD distribution, and therefore cannot exploit a good member EVEN WHEN THAT
MEMBER IS RANKED FIRST.

Discriminating measurement.  Build the top-75 as the shipped top-75 with its j
WORST-scoring members replaced by the j ORACLE-BEST members of the same K=500 pool
(j = 0 is exactly production).  This is a *ranking oracle*: it is what a perfect objective
would have handed the operator, on the SAME candidate set.  Then emit through both
terminal operators and watch what each does with it.

  argmin   : if ranking is the binding constraint, j=1 alone must move 3.454 -> 1.711.
  avg->proj: if the OPERATOR is the binding constraint, j=1 moves ~nothing even though
             the set now contains the pool's best member at rank 1.

Also records, per target and per j, the quantile q of the set's own rr distribution at
which the operator's output lands -- the "fixed quantile" hypothesis, stated numerically.

ORACLE/DIAGNOSTIC throughout: `rr` chooses which members are inserted.  Nothing here is
deployable; it is an attribution experiment, not a method.
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I

JS = [0, 1, 2, 5, 10, 25, 50, 75]


def main(project_js=(0, 1, 5, 75)):
    tg = I.targets()
    rows = []
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = u["W"][p]
        rr = u["rr"][p]
        nat = u["nat_ca"]
        dg = I.distogram(pdb, t["seq"], t["fold"])
        i, j = I.pair_index(t["n"])
        sc = I.shipped_score(dg, I.pair_dists(W, i, j))
        order = np.argsort(sc, kind="stable")
        sub0 = order[:I.M]
        oracle_order = np.argsort(rr, kind="stable")

        r = dict(pdb=pdb, n=t["n"], fold=t["fold"], fail=int(pdb in I.FAIL18),
                 pool_best=float(rr.min()), pool_mean=float(rr.mean()))
        for jj in JS:
            # the jj oracle-best members, then fill from the shipped order, size M
            keep = list(oracle_order[:jj])
            for c in sub0:
                if len(keep) >= I.M:
                    break
                if c not in keep:
                    keep.append(int(c))
            keep = np.asarray(keep[:I.M], int)
            rrk = rr[keep]
            C, b = I.coordinate_average(W[keep])
            avg = I.ca_rmsd(C, nat)
            r[f"j{jj}_set_best"] = float(rrk.min())
            r[f"j{jj}_set_mean"] = float(rrk.mean())
            r[f"j{jj}_avg"] = float(avg)
            r[f"j{jj}_argmin"] = float(rr[keep[0]] if jj > 0 else rr[sub0[0]])
            # where does the operator's output land in the set's own rr distribution?
            r[f"j{jj}_q"] = float((rrk < avg).mean())
            if jj in project_js:
                out = I.project(C, t["seq"], t["fold"])
                r[f"j{jj}_fit"] = float(I.ca_rmsd(out["fit_ca"], nat))
        rows.append(r)
        print(json.dumps({k: r[k] for k in ("pdb", "j0_avg", "j1_avg", "j5_avg", "j75_avg",
                                            "j0_set_best", "j1_set_best", "j0_q", "j75_q")}),
              flush=True)

    names = [r["pdb"] for r in rows]
    fail = np.array([r["fail"] for r in rows], bool)
    folds = np.array([r["fold"] for r in rows])
    agg = {"n": len(rows)}
    for jj in JS:
        a = np.array([r[f"j{jj}_avg"] for r in rows])
        agg[f"j{jj}"] = dict(
            avg=float(a.mean()), avg_fail=float(a[fail].mean()), avg_other=float(a[~fail].mean()),
            set_best=float(np.mean([r[f"j{jj}_set_best"] for r in rows])),
            set_mean=float(np.mean([r[f"j{jj}_set_mean"] for r in rows])),
            argmin=float(np.mean([r[f"j{jj}_argmin"] for r in rows])),
            q_mean=float(np.mean([r[f"j{jj}_q"] for r in rows])),
            q_sd=float(np.std([r[f"j{jj}_q"] for r in rows])),
        )
        if f"j{jj}_fit" in rows[0]:
            f = np.array([r[f"j{jj}_fit"] for r in rows])
            agg[f"j{jj}"]["fit"] = float(f.mean())
            agg[f"j{jj}"]["fit_fail"] = float(f[fail].mean())
            agg[f"j{jj}"]["fit_other"] = float(f[~fail].mean())
    # paired stats for the headline contrasts
    a0 = np.array([r["j0_avg"] for r in rows])
    for jj in JS[1:]:
        agg[f"paired_j{jj}_vs_j0_avg"] = I.paired(np.array([r[f"j{jj}_avg"] for r in rows]), a0,
                                                  folds=folds, names=names)
    if "j0_fit" in rows[0]:
        f0 = np.array([r["j0_fit"] for r in rows])
        for jj in project_js:
            if jj == 0:
                continue
            agg[f"paired_j{jj}_vs_j0_fit"] = I.paired(np.array([r[f"j{jj}_fit"] for r in rows]), f0,
                                                      folds=folds, names=names)
    print(json.dumps({k: v for k, v in agg.items() if not k.startswith("paired")}, indent=1))
    I.write("adv_operator", dict(agg=agg, rows=rows))


if __name__ == "__main__":
    main()
