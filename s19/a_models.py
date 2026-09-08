"""SPRINT 19, AGENT A, Q4 -- THE PREDICTOR ZOO, judged by Ca-RMSD and by error SHAPE.

Six leave-fold-out predictors already exist on disk, all post-fold-repin (Sep 2, after
`peptide_folds.json`), all 183 input features, all trained on the same out-of-fold peptides
plus fold-safe fragments.  They differ in exactly the axes this sprint cares about:

  deployed   dropout 0, uniform per-pair loss weight              <- the incumbent
  d20        dropout 0.2, uniform                                  regularisation
  d20_s1     dropout 0.2, uniform, seed 1                          a second draw
  d20_ens    mean of d20 and d20_s1 probabilities                  seed ensembling
  sw_none    dropout 0.2, `separation_weights("none")` = uniform   the matched control for
  sw_lin     dropout 0.2, `separation_weights("lin")`              THE ANTI-UTILITY LOSS
  pairnet    triangle multiplicative update + axial attention      JOINT CONSISTENCY
  combo      mean of `deployed` and `pairnet` probabilities        DECORRELATED FAMILIES

`sw_lin` vs `sw_none` is a matched pair differing ONLY in the training separation weight,
and `sw_lin` is weighted `min(sep/3, 6)` -- i.e. up to 6x toward LONG-range pairs, which is
the direction the Sprint-18 compass says converts WORSE than uniform.  This is the
pre-registered anti-utility control for A5, and it already exists.

`pairnet` is the pre-registered "joint consistency enforcement" arm: the triangle update
makes each pair's representation depend on the other two sides of every triangle, which is
the architecture-level answer to the realisability defect measured in `a_topo`.

EVERY ARM IS NATIVE-FREE.  The debias is refitted leave-fold-out FOR EACH FAMILY on the FULL
126-target list (the subset trap, brief section 7).  The only ORACLE quantities are the
diagnostics: MAE, kappa, and the residual statistics, each labelled.

Run:  python -m s19.a_models
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

from s19 import a_lib as L                   # noqa: E402
from s19 import a_fit as F                   # noqa: E402
from s19.a_topo import edm_stats, spearman   # noqa: E402
from s12 import instrument as I              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_models.json")

#: (name, kind, checkpoint suffix, dropout).  `kind` = mlp | pairnet | mean-of.
ZOO = [
    ("deployed", "mlp", "_esm_frag.pt", 0.0),
    ("d20", "mlp", "_esm_frag.pt.d20", 0.2),
    ("d20_s1", "mlp", "_esm_frag_s1.pt.d20", 0.2),
    ("sw_none", "mlp", "_esm_frag_sw-none.pt", 0.2),
    ("sw_lin", "mlp", "_esm_frag_sw-lin.pt", 0.2),
    ("pairnet", "pairnet", None, None),
]
BLENDS = [("d20_ens", ["d20", "d20_s1"]), ("combo", ["deployed", "pairnet"])]
NAMES = [z[0] for z in ZOO] + [b[0] for b in BLENDS]


def _load_mlp(fold, suffix, dropout):
    from core import predict as PR
    d_in = PR.feature_width(True)
    return PR.MLP(d_in, seed=0, dropout=dropout).load(
        os.path.join(ROOT, "distogram_models", f"fold{fold}{suffix}"))


def probs_for(tg):
    """(name -> pdb -> prob (npairs, 17)).  Cached, since the models never change."""
    from core import predict as PR
    from core import pipeline as pl
    path = os.path.join(L.CACHE, "zoo_probs.npz")
    if os.path.exists(path):
        z = np.load(path)
        out = {n: {} for n in NAMES}
        for k in z.files:
            n, pdb = k.split("|")
            out[n][pdb] = z[k]
        return out
    pl.guard_esm([t["seq"] for t in tg])
    out = {n: {} for n in NAMES}
    for name, kind, suffix, drop in ZOO:
        models = {}
        for c, t in enumerate(tg):
            f = int(t["fold"])
            if f not in models:
                if kind == "mlp":
                    models[f] = _load_mlp(f, suffix, drop)
                else:
                    models[f] = PR.pairnet_train_fold(f, verbose=False)
            m = models[f]
            if kind == "mlp":
                X, i, j = PR.features(t["seq"], True)
                out[name][t["pdb"]] = m.predict_proba(X).astype(np.float32)
            else:
                dgp = PR.pairnet_distogram(t["seq"], model=m)
                out[name][t["pdb"]] = dgp.prob.astype(np.float32)
        print(f"  {name}: {len(out[name])} targets", flush=True)
    for bname, members in BLENDS:
        for t in tg:
            out[bname][t["pdb"]] = np.mean([out[m][t["pdb"]] for m in members],
                                           axis=0).astype(np.float32)
    np.savez_compressed(path, **{f"{n}|{p}": out[n][p] for n in NAMES for p in out[n]})
    return out


def build_data(tg, probs):
    """name -> pdb -> dict(dhat, sd, ...), with a per-family leave-fold-out sep debias."""
    from core import predict as PR
    centres = np.asarray(PR.CENTRES, float)
    base = C.gather(tg)                       # ORACLE dtrue lives here
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    out = {}
    for name in NAMES:
        fam = {}
        for t in tg:
            p = t["pdb"]
            pr = np.asarray(probs[name][p], float)
            pr = pr / pr.sum(1, keepdims=True)
            exp = (pr * centres).sum(1)
            var = (pr * (centres[None] - exp[:, None]) ** 2).sum(1)
            sd = np.maximum(np.sqrt(np.maximum(var, 1e-6)), 1e-3)
            fam[p] = dict(base[p])
            fam[p]["dhat"] = exp
            fam[p]["sd"] = sd
            fam[p]["prob"] = pr
        deb = L.debias_fns(fam, pdbs, folds)
        for t in tg:
            p = t["pdb"]
            fam[p]["dhat"] = np.maximum(fam[p]["dhat"] - deb[fam[p]["fold"]](fam[p]["sep"]),
                                        2.0)
        out[name] = fam
    return out


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    probs = probs_for(tg)
    fams = build_data(tg, probs)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d0 = fams["deployed"][pdb]
        i, j, nat = d0["i"], d0["j"], d0["nat"]
        dtrue = d0["dtrue"]
        phi0, psi0, avg = F.start(pdb, seq, fold)
        e = {"pdb": pdb, "n": n, "fold": fold, "avg": float(I.ca_rmsd(avg, nat))}
        for name in NAMES:
            d = fams[name][pdb]
            dhat, sd = d["dhat"], d["sd"]
            w = 1.0 / sd ** 2
            rm, fv = F.fit_rmsd(dhat, sd, i, j, phi0, psi0, nat)
            r = dhat - dtrue
            Et = float((w * r ** 2).sum())
            tri, dfc, neg = edm_stats(None, i, j, dhat, n)
            e[name] = rm
            e[name + "_mae"] = float(np.abs(r).mean())          # ORACLE diagnostic
            e[name + "_rms"] = float(np.sqrt((r ** 2).mean()))  # ORACLE diagnostic
            e[name + "_kappa"] = 1.0 - fv / max(Et, 1e-12)      # ORACLE diagnostic
            e[name + "_defect"] = dfc
            e[name + "_tri"] = tri
            e[name + "_sd"] = float(sd.mean())
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_models.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "models", "report")
    g = lambda k: np.array([r[k] for r in rows])                     # noqa: E731
    dep = g("deployed")
    folds = g("fold")
    print(f"\n=== Q4  THE PREDICTOR ZOO   n = {len(rows)} ===")
    print(f"  coordinate average (start)  {g('avg').mean():.3f}")
    print(f"  reproduction gate: deployed = {dep.mean():.3f}  (objceil a0.0 = 3.610)\n")
    print(f"  {'model':<11}{'RMSD':>8}{'med':>7}{'MAE*':>7}{'kappa*':>8}{'EDMdef':>8}"
          f"{'tri%':>7}{'<sd>':>6}   vs deployed (RMSD)")
    tab = {}
    for name in NAMES:
        v = g(name)
        d = v - dep
        m, lo, hi = L.boot(d, rng)
        sg = sum(1 for f in np.unique(folds)
                 if np.sign(d[folds == f].mean()) == np.sign(m) and m != 0)
        print(f"  {name:<11}{v.mean():>8.3f}{np.median(v):>7.3f}{g(name+'_mae').mean():>7.3f}"
              f"{g(name+'_kappa').mean():>8.3f}{g(name+'_defect').mean():>8.3f}"
              f"{100*g(name+'_tri').mean():>7.2f}{g(name+'_sd').mean():>6.2f}"
              f"   {m:+.3f} [{lo:+.3f},{hi:+.3f}] {int((d<0).sum())}W/{int((d>0).sum())}L"
              f" f{sg}/5")
        tab[name] = {"rmsd": float(v.mean()), "median": float(np.median(v)),
                     "mae": float(g(name + "_mae").mean()),
                     "kappa": float(g(name + "_kappa").mean()),
                     "defect": float(g(name + "_defect").mean()),
                     "vs_deployed": {"diff": m, "ci": [lo, hi],
                                     "W": int((d < 0).sum()), "L": int((d > 0).sum()),
                                     "folds_same_sign": sg}}
    print("  * MAE and kappa are ORACLE DIAGNOSTICS (they read the native distances).")

    print("\n  --- THE PRE-REGISTERED MATCHED PAIR: training separation weight ---")
    s, rec = L.report_pair("sw_lin - sw_none", g("sw_lin"), g("sw_none"), rng, folds)
    print("  " + s)
    print("  (sw_lin weights the loss up to 6x toward LONG-range pairs -- the ANTI-utility")
    print("   direction.  A5 predicts sw_lin is WORSE.)")
    tab["_sw_lin_vs_none"] = rec
    s2, rec2 = L.report_pair("MAE sw_lin - sw_none", g("sw_lin_mae"), g("sw_none_mae"), rng)
    print("  " + s2 + "   <- ORACLE, and NOT an outcome")
    tab["_sw_mae"] = rec2

    print("\n  --- DOES RMSD TRACK ACCURACY OR COHERENCE?  (across-model, n = 8 models) ---")
    mm = np.array([tab[n]["mae"] for n in NAMES])
    kk = np.array([tab[n]["kappa"] for n in NAMES])
    vv = np.array([tab[n]["rmsd"] for n in NAMES])
    dd = np.array([tab[n]["defect"] for n in NAMES])
    print(f"    rho(model MAE,   model RMSD) = {spearman(mm, vv):+.3f}")
    print(f"    rho(model kappa, model RMSD) = {spearman(kk, vv):+.3f}")
    print(f"    rho(model defect,model RMSD) = {spearman(dd, vv):+.3f}")
    print("    (n = 8 models -- descriptive, NOT an inferential result)")
    tab["_across_model"] = {"rho_mae": spearman(mm, vv), "rho_kappa": spearman(kk, vv),
                            "rho_defect": spearman(dd, vv)}

    json.dump(tab, open(os.path.join(L.RESULTS, "a_models_report.json"), "w"), indent=1,
              default=float)
    return tab


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
