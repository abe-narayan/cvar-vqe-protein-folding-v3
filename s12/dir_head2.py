"""D5d. Testing the remedy implied by D5b: if the head hurts because its mistakes are
COHERENT (concentrated in the global expand/contract mode), then a head denied every
target-global and row-context feature should make less coherent mistakes.

`local` head: the same model on the pair-local features only -- separation, position,
E[d], sd, the shape of the pair's own predicted distribution, its ESM contact cell and
neighbourhood, the two residues' propensities and embeddings.  Dropped: rg, rg_rel, gexp,
gsd, gent, exp_m_gexp, shres_*, rowexp_*, rowsd_*, deg_*, sh_*.

Reports LFO accuracy, the coherence statistics of D5b, and the emitted ladder with the
matched random-at-same-accuracy null.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)

from s12 import instrument as I
from s12 import obj_common as OC
from s12 import dir_common as DC
from s12 import dir_head as DH
from s12 import dir_why as DW

DROP = ("rg", "rg_rel", "gexp", "gsd", "gent", "exp_m_gexp")
DROP_PREFIX = ("shres_", "rowexp_", "rowsd_", "deg_", "sh_i", "sh_j")
LADDER = (2.0, 3.0, 4.0)


def main():
    z, keys, names = DH.load_bank()
    keep = [k for k, nm in enumerate(names)
            if nm not in DROP and not any(nm.startswith(p) for p in DROP_PREFIX)]
    print(f"local head: {len(keep)} of {len(names)} features", flush=True)
    pdbs = [k[0] for k in keys]; folds = np.array([int(k[3]) for k in keys])
    Xs = {p: np.asarray(z[f"{p}/X"], np.float32)[:, keep] for p in pdbs}
    ys = {p: np.asarray(z[f"{p}/y"], np.int8) for p in pdbs}
    prob = {}
    for f in range(5):
        tr = [k for k in range(len(pdbs)) if folds[k] != f]
        te = [k for k in range(len(pdbs)) if folds[k] == f]
        Xtr = np.vstack([Xs[pdbs[k]] for k in tr]); ytr = np.concatenate([ys[pdbs[k]] for k in tr])
        Xte = np.vstack([Xs[pdbs[k]] for k in te])
        cuts = np.cumsum([len(ys[pdbs[k]]) for k in te])[:-1]
        pv, _ = DH.fit_predict(Xtr, ytr, Xte)
        for kk, part in zip(te, np.split(pv, cuts)):
            prob[pdbs[kk]] = part.astype(np.float32)
        print(f"  fold {f} done free={I.free_gb():.1f}", flush=True)
        del Xtr, Xte

    tgset = {t["pdb"] for t in I.targets()}
    groups = {"corpus787": pdbs, "tuning126": [p for p in pdbs if p in tgset],
              "fail18": [p for p in pdbs if p in set(I.FAIL18)],
              "other108": [p for p in pdbs if p in tgset and p not in set(I.FAIL18)]}
    res = {"n_features": len(keep), "acc": {}}
    for g, sel in groups.items():
        res["acc"][g] = float(np.mean(np.concatenate(
            [(prob[q] > 0.5).astype(int) == ys[q] for q in sel])))
    print("local-head accuracy:", {k: round(v, 4) for k, v in res["acc"].items()}, flush=True)

    rows, t0 = [], time.time()
    for k, t in enumerate(I.targets()):
        pdb, n = t["pdb"], t["n"]
        while I.free_gb() < 1.5:
            print("  waiting on memory", flush=True); time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue = d["exp"], d["dtrue"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        sl = np.where(np.asarray(prob[pdb], float) > 0.5, 1.0, -1.0)
        acc = float((sl == s_true).mean())
        st = DW.field_stats(np.where(sl != s_true, sl, 0.0), d["i"], d["j"], n)
        cells = {"bayes": OC.emit(d, OC.score_bayes(d["D"], d["risk"], d["grid"]), lam=None)["avg_rmsd"]}
        for D in LADDER:
            cells[f"local@{D}"] = DC.emit_target(d, DC.sign_target(exp, sl, D), lam=None)["avg_rmsd"]
            v = []
            for s in range(3):
                rng = np.random.default_rng(5000 + 17 * k + s)
                v.append(DC.emit_target(d, DC.sign_target(exp, DC.corrupt_sign(s_true, 1 - acc, rng), D),
                                        lam=None)["avg_rmsd"])
            cells[f"rand@{D}"] = float(np.mean(v))
        rows.append(dict(pdb=pdb, fold=t["fold"], fail18=pdb in I.FAIL18, acc=acc,
                         field=st, cells=cells))
        if (k + 1) % 30 == 0 or k == 0:
            print(f"  emit {k+1}/126 [{time.time()-t0:.0f}s]", flush=True)
    res["rows"] = rows
    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    res["summary"] = {}
    for g, rs in grp.items():
        by = np.array([r["cells"]["bayes"] for r in rs])
        tab = {"acc": float(np.mean([r["acc"] for r in rs])),
               "scale_frac": float(np.mean([r["field"]["scale_frac"] for r in rs])),
               "top1": float(np.mean([r["field"]["top1"] for r in rs])),
               "eff_rank": float(np.mean([r["field"]["eff_rank"] for r in rs]))}
        for D in LADDER:
            for a in ("local", "rand"):
                x = np.array([r["cells"][f"{a}@{D}"] for r in rs])
                s = I.paired(x, by)
                tab[f"{a}@{D}"] = dict(d=s["mean_diff"], ci=s["ci95"],
                                       wl=[s["n_better"], s["n_worse"]],
                                       drop10=s["drop_top10_mean_diff"])
        res["summary"][g] = tab
        print(f"\n== {g} == acc={tab['acc']:.4f} scale_frac={tab['scale_frac']:.4f} "
              f"top1={tab['top1']:.4f} eff_rank={tab['eff_rank']:.3f}")
        for D in LADDER:
            print(f"  delta={D}: local {tab[f'local@{D}']['d']:+.3f} "
                  f"{tab[f'local@{D}']['wl']}  rand {tab[f'rand@{D}']['d']:+.3f}")
    I.write("dir_head2", res)


if __name__ == "__main__":
    main()
