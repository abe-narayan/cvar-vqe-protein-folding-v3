"""SPRINT 14, ENER-9 -- pricing the VALIDATION role in the currency that ships.

E2c and E8 both say the same thing: AMBER's only genuine skill is rejecting bad geometry.
E1's veto arms tested that through the 1-local prior as a selector, which is itself worse
than a random draw, so they could not isolate it.

This module prices the validation role directly, in the currency the project has shown the
terminal operator actually consumes.  The established operator law is

    d_out = 1.16 * d_set_mean + 0.04 * d_set_best,   R^2 0.893

so what a FILTER has to do is move the surviving set's MEAN, not its best.  Measured here:
remove the worst x% of the pool by each objective and report the surviving set's mean and
best, against removing x% AT RANDOM -- the only null that controls for the fact that
shrinking a set changes its statistics.

    python -m s14.ener_filter
"""
from __future__ import annotations

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E
from s14.ener_matrix import ensemble

FRACS = [0.10, 0.25, 0.50, 0.75, 0.90]
GUARDS = ["legacy", "amber", "leg_steric", "amb_nonbonded", "prior", "legacy_nosteric"]
N_NULL = 200


def curve(pdb):
    d = ensemble(pdb)
    R = d["rmsd"]
    n = len(R)
    rng = np.random.default_rng(E.SEED)
    out = {}
    for g in GUARDS:
        v = np.asarray(d[g], float)
        order = np.argsort(v, kind="mergesort")          # best energy first
        rows = []
        for f in FRACS:
            m = max(int(round((1 - f) * n)), 2)          # keep the best (1-f) by energy
            keep = order[:m]
            rows.append(dict(remove_frac=f, n_kept=m,
                             mean=float(R[keep].mean()), best=float(R[keep].min()),
                             operator=float(1.16 * R[keep].mean() + 0.04 * R[keep].min())))
        out[g] = rows
    nulls = []
    for f in FRACS:
        m = max(int(round((1 - f) * n)), 2)
        mm = [R[rng.choice(n, m, replace=False)] for _ in range(N_NULL)]
        nulls.append(dict(remove_frac=f, n_kept=m,
                          mean=float(np.mean([x.mean() for x in mm])),
                          best=float(np.mean([x.min() for x in mm])),
                          operator=float(np.mean([1.16 * x.mean() + 0.04 * x.min()
                                                  for x in mm]))))
    return dict(pdb=pdb, fold=int(d["_z"].fold), full_mean=float(R.mean()),
                full_best=float(R.min()), per_guard=out, null=nulls)


def main():
    rows = [curve(p) for p in E.ENUM_TARGETS]
    folds = [r["fold"] for r in rows]
    print("=== FILTER CURVE: remove the worst x% by each objective ===")
    print("The null removes x% AT RANDOM, which is the only control for set shrinkage.")
    print(f"full pool mean {np.mean([r['full_mean'] for r in rows]):.3f}  "
          f"best {np.mean([r['full_best'] for r in rows]):.3f}\n")
    for fi, f in enumerate(FRACS):
        nl = np.array([r["null"][fi]["mean"] for r in rows])
        nlb = np.array([r["null"][fi]["best"] for r in rows])
        nlo = np.array([r["null"][fi]["operator"] for r in rows])
        print(f"-- remove {int(f*100)}%  (random null: mean {nl.mean():.3f}  "
              f"best {nlb.mean():.3f}  operator {nlo.mean():.3f}) --")
        for g in GUARDS:
            mv = np.array([r["per_guard"][g][fi]["mean"] for r in rows])
            bv = np.array([r["per_guard"][g][fi]["best"] for r in rows])
            ov = np.array([r["per_guard"][g][fi]["operator"] for r in rows])
            p = I.paired(mv, nl, folds=folds)
            po = I.paired(ov, nlo, folds=folds)
            print(f"   {g:18s} mean {mv.mean():6.3f} ({p['mean_diff']:+.3f} "
                  f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']})"
                  f"   best {bv.mean():6.3f} ({bv.mean()-nlb.mean():+.3f})"
                  f"   operator {po['mean_diff']:+.3f} "
                  f"[{po['ci95'][0]:+.3f},{po['ci95'][1]:+.3f}]")
    E.write("ener_filter", dict(
        what="value of each objective as a pool FILTER, against a random-removal null",
        operator_law="d_out = 1.16*d_set_mean + 0.04*d_set_best (project finding)",
        fracs=FRACS, per_target=rows), n_expected=len(E.ENUM_TARGETS))
    return rows


if __name__ == "__main__":
    main()
