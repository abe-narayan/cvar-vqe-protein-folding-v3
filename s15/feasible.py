"""SPRINT 15, coordinator -- IS THE NATIVE STRUCTURE FEASIBLE UNDER ITS OWN PREDICTED RESTRAINTS?

THE QUESTION THIS ANSWERS, AND WHY IT COMES BEFORE ANY MORE OPTIMISATION.

K1-CORRECTED measured that fitting continuous torsions to the raw predicted distogram gives
3.644 A and LOSES to the incumbent by +0.440 [+0.290, +0.592], while fitting to the true
distance matrix through the identical machinery gives 0.611 A. The whole 3.03 A is restraint
error. But "restraint error" is not yet a mechanism, and there are two very different
mechanisms it could be, with opposite engineering consequences:

  (A) THE RESTRAINTS ARE NOISY BUT UNBIASED ABOUT THE TRUTH. The native still sits near the
      objective's minimum; the optimiser lands elsewhere because the basin is wide and flat.
      Consequence: better search, better aggregation, more restraints all help.

  (B) THE RESTRAINTS ARE INCONSISTENT WITH THE TRUTH. The native scores WORSE on the objective
      than structures far from it, so the objective's minimum is somewhere else entirely.
      Consequence: **search is actively harmful**, better optimisation moves further from the
      answer, and nothing but better restraints can help.

This project has seen (B) before, on the other side of the selection/generation divide. The
standing result "the objective does not rank the native" measured the native as the distance
objective's argmin on **3 of 126 targets**, at the **36.8th percentile** -- and that single
fact explained the recurring "optimise harder, get worse" pathology. The generative frame was
adopted precisely to escape it. Whether it actually does is not something to assume.

THE EXPERIMENT IS ALMOST FREE, because it needs no optimisation at all. Every quantity is an
evaluation of the objective at structures we already have:

    f(native)              the objective AT the native CA trace                  ORACLE-scored
    f over the K=500 pool  the objective at every retrieved window
    the pool's own RMSDs   which member is actually closest to the native        ORACLE

From those three we can ask everything that matters:

    * what PERCENTILE of the pool does the native occupy on the objective?  (50 = the objective
      is blind to nativeness; << 50 = it favours the native; >> 50 = it actively disfavours it)
    * does the objective's argmin over the pool have a lower RMSD than a random pool member?
      This is the Roget et al. (arXiv:2606.21241) test, which found the minimum-cost
      conformation has on average a LARGER RMSD than a random feasible one, for a lattice
      contact potential at exactly our chain lengths. Our restraint objective is not a contact
      potential, so their mechanism should not apply -- but that has to be MEASURED, not
      asserted, and it is the single cleanest way to position our work against theirs.
    * is the native FEASIBLE, in the constrained sense that Family B needs: does it satisfy
      |d_ij - dhat_ij| <= eps * sd_ij for a given eps, and on what fraction of pairs? Family B
      is `min E_AMBER subject to C_restraint <= epsilon`, and if the native is infeasible at
      every usable epsilon then Family B's feasible set excludes the answer and the family is
      falsified before it is built.

ARMS. Every restraint channel this sprint has, so the comparison is like-for-like:

    ls_pred     weighted least squares on the distogram's mean, w = 1/sd^2      NATIVE-FREE
    ls_debias   the same, on separation-debiased means, leave-fold-out          NATIVE-FREE
    ml_pred     negative log-likelihood under the full 17-bin distribution      NATIVE-FREE
    ls_pool     weighted least squares on the retrieval pool's median distances NATIVE-FREE
    ml_pool     negative log-likelihood under the pool's empirical histogram    NATIVE-FREE
    ml_comb     the product of the distogram and pool likelihoods               NATIVE-FREE
    ORACLE_true least squares on the true distances -- must put the native at percentile 0

The ORACLE arm is a machinery check with a known answer: if it does not place the native at
the very bottom of the pool, something is wrong with the evaluation, not with the science.

Run:
    python -m s15.feasible
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distml as M                  # noqa: E402
from s15 import pooldist as P                # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

ARMS = ("ls_pred", "ls_debias", "ml_pred", "ls_pool", "ml_pool", "ml_comb", "ORACLE_true")
EPS = (0.5, 1.0, 2.0, 3.0)


def _ls(D, target, w):
    """(m, npairs) distances -> (m,) weighted squared residual."""
    r = D - np.asarray(target, float)[None, :]
    return (r * r * np.asarray(w, float)[None, :]).sum(1)


def evaluate(t, debias=None):
    """Every arm's objective at the native and at all K=500 pool members. One target."""
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    i, j = I.pair_index(n)
    sep = (j - i).astype(float)

    Dpool, _sim, _i, _j = P.pool_distances(pdb, n)                 # (m, npairs)
    W = np.asarray(u["W"], float)[I.pool_idx(u, 500)]
    rms_pool = I.kabsch_rmsd_batch(W, nat)                          # ORACLE, post hoc
    d_nat = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))

    dg = I.distogram(pdb, seq, fold)
    dhat = np.asarray(dg["expected"], float)
    sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
    centres = np.asarray(dg["centres"], float)
    w = 1.0 / sd ** 2
    d_deb = dhat - debias(sep) if debias is not None else dhat

    d_pool = np.median(Dpool, axis=0)
    w_pool = 1.0 / np.maximum(Dpool.std(axis=0), 0.25) ** 2
    tab_dg = M.LogPTable(np.asarray(dg["prob"], float), centres)
    tab_pl = P.hist_table(Dpool, centres)
    tab_cb = P.SumLogP(tab_dg, tab_pl)

    f_nat, f_pool = {}, {}
    f_nat["ls_pred"] = float(_ls(d_nat[None], dhat, w)[0])
    f_pool["ls_pred"] = _ls(Dpool, dhat, w)
    f_nat["ls_debias"] = float(_ls(d_nat[None], d_deb, w)[0])
    f_pool["ls_debias"] = _ls(Dpool, d_deb, w)
    f_nat["ls_pool"] = float(_ls(d_nat[None], d_pool, w_pool)[0])
    f_pool["ls_pool"] = _ls(Dpool, d_pool, w_pool)
    f_nat["ORACLE_true"] = float(_ls(d_nat[None], d_nat, np.ones(i.size))[0])
    f_pool["ORACLE_true"] = _ls(Dpool, d_nat, np.ones(i.size))
    for name, tab in (("ml_pred", tab_dg), ("ml_pool", tab_pl), ("ml_comb", tab_cb)):
        f_nat[name] = float(-tab.rowsum(d_nat[None])[0])
        f_pool[name] = -tab.rowsum(Dpool)

    row = {"pdb": pdb, "n": n, "fold": fold, "m": int(Dpool.shape[0]),
           "rmsd_pool_best": float(rms_pool.min()),
           "rmsd_pool_mean": float(rms_pool.mean()), "arms": {}}
    for a in ARMS:
        fp = np.asarray(f_pool[a], float)
        k = int(np.argmin(fp))
        #: the native's rank among pool members ON THE OBJECTIVE. 0 = the objective's own best.
        pct = float((fp < f_nat[a]).mean())
        row["arms"][a] = {
            "f_native": f_nat[a],
            "native_percentile": pct,
            "argmin_rmsd": float(rms_pool[k]),
            "argmin_minus_random": float(rms_pool[k] - rms_pool.mean()),
            "rho_f_rmsd": float(_spearman(fp, rms_pool)),
            "beats_random": bool(rms_pool[k] < rms_pool.mean()),
        }
    # ---- Family B feasibility: is the native inside the restraint band at all?
    z = np.abs(d_nat - dhat) / sd
    row["feasibility"] = {f"eps_{e}": float((z <= e).mean()) for e in EPS}
    row["feasibility"]["max_z"] = float(z.max())
    row["feasibility"]["median_z"] = float(np.median(z))
    return row


def _spearman(a, b):
    from s14 import vqe_lib as V
    return V.spearman(np.asarray(a, float), np.asarray(b, float))


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)

    # leave-fold-out separation debias, exactly as s15/distcal.py fits it
    data = C.gather(tg)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    path = os.path.join(RESULTS, "feasible.json")
    for c, t in enumerate(tg):
        rows.append(evaluate(t, debias=deb[int(t["fold"])]))
        if (c + 1) % 20 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": rows, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)}", flush=True)

    out = {"n": len(rows), "rows": rows, "arms": {}}
    for a in ARMS:
        pct = np.asarray([r["arms"][a]["native_percentile"] for r in rows])
        amr = np.asarray([r["arms"][a]["argmin_rmsd"] for r in rows])
        adf = np.asarray([r["arms"][a]["argmin_minus_random"] for r in rows])
        rho = np.asarray([r["arms"][a]["rho_f_rmsd"] for r in rows])
        out["arms"][a] = {
            "native_percentile_mean": float(pct.mean()),
            "native_percentile_median": float(np.median(pct)),
            "frac_native_is_argmin": float((pct == 0.0).mean()),
            "frac_native_below_median": float((pct < 0.5).mean()),
            "argmin_rmsd_mean": float(amr.mean()),
            "argmin_minus_random_mean": float(adf.mean()),
            "argmin_beats_random_frac": float(np.mean([r["arms"][a]["beats_random"]
                                                       for r in rows])),
            "rho_f_rmsd_mean": float(rho.mean()),
            "rho_f_rmsd_frac_positive": float((rho > 0).mean()),
            "vs_random_paired": I.paired(
                amr, np.asarray([r["rmsd_pool_mean"] for r in rows]),
                folds=folds, names=pdbs),
        }
    feas = {k: float(np.mean([r["feasibility"][k] for r in rows]))
            for k in rows[0]["feasibility"]}
    out["feasibility_mean"] = feas

    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_feasible", out, n_expected=len(tg))

    print(f"\nn = {len(rows)} targets, K=500 pool members each. "
          f"'native pct' is the native's rank among pool members ON THE OBJECTIVE:\n"
          f"0.5 means the objective is blind to nativeness, above 0.5 means it DISFAVOURS "
          f"the native.\n")
    print(f"{'arm':<14}{'native pct':>12}{'median':>9}{'is argmin':>11}"
          f"{'argmin RMSD':>13}{'vs random':>11}{'wins':>7}{'rho(f,rmsd)':>13}")
    for a in ARMS:
        s = out["arms"][a]
        print(f"{a:<14}{s['native_percentile_mean']:>12.3f}"
              f"{s['native_percentile_median']:>9.3f}"
              f"{s['frac_native_is_argmin']:>11.3f}"
              f"{s['argmin_rmsd_mean']:>13.3f}"
              f"{s['argmin_minus_random_mean']:>+11.3f}"
              f"{s['argmin_beats_random_frac']:>7.2f}"
              f"{s['rho_f_rmsd_mean']:>13.3f}")
    print(f"\nFamily B feasibility -- fraction of pairs with |d_native - dhat| <= eps*sd:")
    for e in EPS:
        print(f"  eps = {e}:  {feas[f'eps_{e}']:.3f}")
    print(f"  median z = {feas['median_z']:.3f}   max z = {feas['max_z']:.3f}")
    print("\n(a restraint set whose band excludes the native at every usable eps falsifies "
          "Family B\n before it is built; a native percentile above 0.5 means optimising the "
          "objective moves AWAY\n from the answer, which is mechanism (B) in the docstring.)")
    return out


if __name__ == "__main__":
    run()
