"""SPRINT 14, OBJ -- is the learned objective anything more than a compactness axis?

The flip diagnostic showed that the leave-one-target-out learned objective's per-target
global rank correlation with RMSD tracks the target's own rho(contact count, RMSD) at
Spearman +0.951, and the native's z-scored radius of gyration at +0.909.  That is the
signature of a ONE-DIMENSIONAL globular/extended detector wearing a 4069-parameter coat.

This module prices that directly:

  * `rg_fixed`   -- a one-feature objective, score = -(contacts_8A per residue), i.e.
                    "more compact is better", the sign the training majority supports.
                    INFERENCE-LEGAL.
  * `rg_oracle`  -- the same feature with the sign that helps THIS target.
                    ORACLE DIAGNOSTIC; an upper bound on what any compactness axis can do.
  * residualised rho -- partial out [1, Rg, Rg^2, end-to-end, contacts8, contacts8^2]
                    from BOTH the learned score and the true RMSD, within each target,
                    and re-measure the rank correlation of the residuals.
                    If this is ~0 the model carries no structural information beyond a
                    compactness scalar.

    python -m s14.obj_resid [PDB ...]
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
from s14 import obj_model as M                     # noqa: E402
from s14 import obj_train as T                     # noqa: E402
from s14 import obj_exp as X                       # noqa: E402


def shape_basis(fz, pdb_id, S):
    """The global shape covariates a compactness detector could be using."""
    C, phi, psi = fz.geometry(pdb_id, S)
    tf = fz.tf(pdb_id)
    cen = C - C.mean(1, keepdims=True)
    rg = np.sqrt((cen ** 2).sum(2).mean(1))
    e2e = np.linalg.norm(C[:, 0] - C[:, -1], axis=1)
    d = np.linalg.norm(C[:, tf.ii, :] - C[:, tf.jj, :], axis=2)
    far = tf.sep >= 3
    c8 = (d[:, far] < 8.0).sum(1) / tf.n
    Z = np.column_stack([np.ones_like(rg), rg, rg ** 2, e2e, c8, c8 ** 2])
    return Z, rg, c8


def residualise(v, Z):
    coef, *_ = np.linalg.lstsq(Z, v, rcond=None)
    return v - Z @ coef


def run(targets, seed=0, lam=0.03, n_uniform=20000):
    fz = M.Featurizer(cache_dir=True)
    rows = []
    t0 = time.time()
    print(f"{'pdb':6s} {'rho_g':>7s} {'rho_g|shape':>12s} {'rho_dec':>8s} "
          f"{'rho_dec|shape':>14s} {'rgfix_rho':>10s} {'rgorc_rho':>10s} "
          f"{'d_dec_learn':>12s} {'d_dec_rgfix':>12s} {'d_dec_rgorc':>12s}")
    for p in targets:
        train = [q for q in targets if q != p]
        w = T.fit(train, fz, n_uniform=n_uniform, lam=lam, center=True)
        en = T.Enum(p)
        idx = X.eval_idx(en)
        S = T._decode(idx, en.n, en.k)
        e = fz.score_with(w, p, S)
        r = en.rmsd[idx].astype(np.float64)
        Z, rg, c8 = shape_basis(fz, p, S)

        e_res, r_res = residualise(e, Z), residualise(r, Z)
        o = np.argsort(e)
        dec = o[:len(e) // 10]
        rho_g, rho_d = spearman(e, r), spearman(e[dec], r[dec])
        rho_gs = spearman(e_res, r_res)
        rho_ds = spearman(e_res[dec], r_res[dec])

        # one-feature compactness objectives
        fix = -c8                                       # "more compact is better"
        orc = fix * (1.0 if spearman(fix, r) > 0 else -1.0)   # ORACLE sign
        out = dict(pdb=p, fold=en.fold, rho_g=rho_g, rho_g_shape=rho_gs,
                   rho_dec=rho_d, rho_dec_shape=rho_ds,
                   rgfix_rho=spearman(fix, r), rgorc_rho=spearman(orc, r),
                   space_mean=float(r.mean()))
        for nm, v in (("learn", e), ("rgfix", fix), ("rgorc", orc)):
            oo = np.argsort(v)
            out["d_dec_" + nm] = float(r[oo[:len(v) // 10]].mean() - r.mean())
            out["d_t100_" + nm] = float(r[oo[:100]].mean() - r.mean())
        rows.append(out)
        print(f"{p:6s} {rho_g:+7.3f} {rho_gs:+12.3f} {rho_d:+8.3f} {rho_ds:+14.3f} "
              f"{out['rgfix_rho']:+10.3f} {out['rgorc_rho']:+10.3f} "
              f"{out['d_dec_learn']:+12.3f} {out['d_dec_rgfix']:+12.3f} "
              f"{out['d_dec_rgorc']:+12.3f}  [{time.time()-t0:.0f}s]", flush=True)

    A = {k: np.array([r[k] for r in rows], float) for k in rows[0] if k not in ("pdb",)}
    folds = A["fold"].astype(int)
    print("\n-- MEANS --")
    for k in ("rho_g", "rho_g_shape", "rho_dec", "rho_dec_shape", "rgfix_rho",
              "rgorc_rho", "d_dec_learn", "d_dec_rgfix", "d_dec_rgorc",
              "d_t100_learn", "d_t100_rgfix", "d_t100_rgorc"):
        print(f"  {k:16s} {A[k].mean():+.4f}   (sd {A[k].std():.3f}, "
              f"|.|>0.2 on {(np.abs(A[k])>0.2).sum()}/{len(rows)})")
    print("\n-- paired: learned d_decile vs RANDOM(0) --")
    pr = I.paired(A["d_dec_learn"], np.zeros(len(rows)), folds=folds,
                  names=[r["pdb"] for r in rows])
    print(f"   diff={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] "
          f"W/L={pr['n_better']}/{pr['n_worse']} folds={pr.get('per_fold')}")
    print("-- paired: learned d_decile vs rgfix d_decile --")
    pr = I.paired(A["d_dec_learn"], A["d_dec_rgfix"], folds=folds,
                  names=[r["pdb"] for r in rows])
    print(f"   diff={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] "
          f"W/L={pr['n_better']}/{pr['n_worse']}")
    print("-- paired: ORACLE-sign compactness d_decile vs RANDOM(0) --")
    pr = I.paired(A["d_dec_rgorc"], np.zeros(len(rows)), folds=folds,
                  names=[r["pdb"] for r in rows])
    print(f"   diff={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] "
          f"W/L={pr['n_better']}/{pr['n_worse']}")
    I.write("s14_obj_resid", {"lam": lam, "rows": rows}, n_expected=len(targets))
    return rows


if __name__ == "__main__":
    tg = [a for a in sys.argv[1:] if not a.startswith("-")]
    run(tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])])
