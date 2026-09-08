"""TORSION-PREDICTOR -- leave-fold-out training of every predictor arm and null.

One process, torch threads = 2, OMP = 2.  Each arm produces, for all 126 tuning targets, a
(n, 324) posterior density on the 18x18 Ramachandran grid, cached to
`s13/cache/tors_post_<arm>.npz`.  Everything downstream reads those caches, so the training
runs once.

ARMS
    p_point    circular point regression (cos/sin), vM density from its own shrinkage
    p_grid     324-way categorical over the 20-degree Ramachandran grid    <- primary
    p_mvm      mixture of 8 bivariate von Mises, trained by NLL
    p_abego4   S12's 4-state ABEGO predictor, reproduced as the reference point
    d_pep      p_grid trained on PEPTIDE residues only     (S7-2 distribution-shift test)
    d_frag     p_grid trained on FRAGMENT residues only
    e_esm      p_grid + ESM-2 PCA-32 (centre + context mean)
    e_esmpep   e_esm on peptide residues only
    n_shuf     MANDATORY NULL: real features, labels permuted
    n_comp     MANDATORY NULL: composition-only features
    n_marg     MANDATORY NULL: the corpus-marginal Ramachandran density, no sequence at all

    python -m s13.tors_train [arm ...]
"""
from __future__ import annotations
import os, sys, json, time, gc

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s13 import tors_common as T          # noqa: E402
from s12 import instrument as I           # noqa: E402

ARMS = {
    "p_point":  dict(family="point",  feat="context", subset="all",  shuffle=False),
    "p_grid":   dict(family="grid",   feat="context", subset="all",  shuffle=False),
    "p_mvm":    dict(family="mvm",    feat="context", subset="all",  shuffle=False),
    "p_abego4": dict(family="abego4", feat="context", subset="all",  shuffle=False),
    "d_pep":    dict(family="grid",   feat="context", subset="pep",  shuffle=False),
    "d_frag":   dict(family="grid",   feat="context", subset="frag", shuffle=False),
    "e_esm":    dict(family="grid",   feat="esm",     subset="all",  shuffle=False),
    "e_esmpep": dict(family="grid",   feat="esm",     subset="pep",  shuffle=False),
    "n_shuf":   dict(family="grid",   feat="context", subset="all",  shuffle=True),
    "n_comp":   dict(family="grid",   feat="comp",    subset="all",  shuffle=False),
    "n_marg":   dict(family="marg",   feat="none",    subset="all",  shuffle=False),
    "a_noPos":  dict(family="grid",   feat="ctxonly", subset="all",  shuffle=False),
    "a_pepPos": dict(family="grid",   feat="ctxonly", subset="pep",  shuffle=False),
}
EPOCHS = 30
SEED = 0


def post_path(arm):
    return os.path.join(T.CACHE, f"tors_post_{arm}.npz")


def fold_path(arm, f):
    return os.path.join(T.CACHE, f"tors_fold_{arm}_f{f}.npz")


def run_arm(arm, seed=SEED, epochs=EPOCHS):
    cfg = ARMS[arm]
    tg = I.targets()
    out, diag = {}, {}
    for f in sorted({t["fold"] for t in tg}):
        fp = fold_path(arm, f)
        if os.path.exists(fp):                       # per-fold resume: the box is shared
            z = np.load(fp, allow_pickle=True)
            for k in z.files:
                if k != "_diag":
                    out[k] = z[k]
            diag[f] = json.loads(str(z["_diag"]))
            print(f"  [resume] {arm} fold {f}", flush=True)
            continue
        T.free_ok(1.5, tag=f"{arm} fold{f}")
        c = T.corpus(f)
        ok = c["ok"]
        sub = {"all": np.ones(len(ok), bool), "pep": c["org"], "frag": ~c["org"]}[cfg["subset"]]
        m = ok & sub
        T.set_abego_map(c["phi"], c["psi"], ok)

        if cfg["family"] == "marg":
            P0 = T.marginal_posterior(c["phi"], c["psi"], m, 1)[0]
            for t in tg:
                if t["fold"] == f:
                    out[t["pdb"]] = np.tile(P0[None], (t["n"], 1)).astype(np.float32)
            diag[f] = {"n_train": int(m.sum()), "val_loss": None, "secs": 0.0}
            print(f"  {arm} fold {f}: marginal over {int(m.sum())} residues", flush=True)
            np.savez_compressed(fold_path(arm, f), _diag=np.array(json.dumps(diag[f])),
                                **{t["pdb"]: out[t["pdb"]] for t in tg if t["fold"] == f})
            del c
            continue

        E = T.esm_matrix_corpus(f)[c["keep"]] if cfg["feat"] == "esm" else None
        X = T.featurise(c["codes"], c["pid"], cfg["feat"], E)[m]
        Y = T.make_labels(cfg["family"], c["phi"][m], c["psi"][m])
        if cfg["shuffle"]:
            Y = np.random.default_rng(seed + 991).permutation(Y)
        t0 = time.time()
        mdl, st = T.train_one(X, Y, cfg["family"], seed=seed, epochs=epochs,
                              groups=c["pid"][m])
        secs = time.time() - t0
        diag[f] = {**st, "secs": round(secs, 1), "n_feat": int(X.shape[1])}
        print(f"  {arm} fold {f}: n={m.sum()} d={X.shape[1]} val={st['val_loss']:.4f} "
              f"({secs:.0f}s, {st['epochs_run']} ep) free={I.free_gb():.2f}", flush=True)
        del X, Y
        for t in tg:
            if t["fold"] != f:
                continue
            cds = T.codes_of(t["seq"])
            Et = T.esm_target(t["seq"]) if cfg["feat"] == "esm" else None
            Xt = T.featurise(cds, None, cfg["feat"], Et)
            out[t["pdb"]] = T.posterior(mdl, Xt, cfg["family"])
        np.savez_compressed(fold_path(arm, f), _diag=np.array(json.dumps(diag[f])),
                            **{t["pdb"]: out[t["pdb"]] for t in tg if t["fold"] == f})
        del mdl, c
        gc.collect()
    np.savez_compressed(post_path(arm), **out)
    return diag


def main(argv):
    arms = argv or list(ARMS)
    disc = T.assert_fold_discipline()
    alld = {}
    dpath = os.path.join(T.RESULTS, "tors_train_diag.json")
    if os.path.exists(dpath):
        alld = json.load(open(dpath)).get("diag", {})
    for a in arms:
        if os.path.exists(post_path(a)) and a in alld:
            print(f"[skip] {a}", flush=True)
            continue
        print(f"[arm] {a}  free={I.free_gb():.2f}", flush=True)
        alld[a] = run_arm(a)
        json.dump({"what": "leave-fold-out torsion predictor training",
                   "epochs": EPOCHS, "seed": SEED, "arms": ARMS,
                   "fold_discipline": disc, "diag": alld},
                  open(dpath, "w"), indent=1, default=str)
    print("wrote", dpath)


if __name__ == "__main__":
    main(sys.argv[1:])
