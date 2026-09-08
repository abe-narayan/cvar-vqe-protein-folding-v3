"""D3c. Train the per-pair error-SIGN head, leave-fold-out.

Label  y = 1  iff  d_true > E[d]   (the distogram UNDER-predicted this distance).
Training corpus: every peptide in `peptide_db` (787), each scored by the distogram fold
model that never saw its fold (`dir_corpus.py`).  The head that scores a target in fold f
is trained on folds != f, so no target's own errors -- and no error from any sequence in
its identity-clustered fold -- ever trains the head that judges it.

Arms
  full      HistGradientBoosting on all deployable pair features        (the head)
  sep       the same model on separation features only                  NULL (b)
  const     always predict the majority sign of the TRAINING folds      trivial baseline
  shellmaj  per-separation-shell majority sign of the TRAINING folds     NULL (b), cheap form
  shuf      full features, labels shuffled WITHIN each training fold    NULL (a)

Writes s12/results/dir_head.json and s12/cache/dir_head_pred.npz (per-target sign
probabilities for the 126 tuning targets, for dir_emit.py).
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
from s12 import dir_feats as DF

CORPUS = os.path.join(ROOT, "s12", "cache", "dir_corpus.npz")
FEATS = os.path.join(ROOT, "s12", "cache", "dir_featbank.npz")
PRED = os.path.join(ROOT, "s12", "cache", "dir_head_pred.npz")


# --------------------------------------------------------------------------- data
def build_featbank(force=False):
    """Features + labels for every corpus peptide, cached."""
    if os.path.exists(FEATS) and not force:
        return
    z = np.load(CORPUS, allow_pickle=True)
    keys = [str(k).split("\t") for k in z["keys"]]
    centres = np.asarray(__import__("core.predict", fromlist=["x"]).CENTRES, float)
    tg = {t["pdb"]: t for t in I.targets()}
    out, names = {}, None
    t0 = time.time()
    for k, (pdb, seq, n, fold) in enumerate(keys):
        n = int(n); fold = int(fold)
        if pdb in tg:                              # inference view: the SHIPPED distogram
            dg = I.distogram(pdb)
            prob, exp, sd = dg["prob"], np.asarray(dg["expected"], float), np.asarray(dg["sd"], float)
        else:
            prob = np.asarray(z[f"{pdb}/prob"], np.float32)
            exp = np.asarray(z[f"{pdb}/exp"], float); sd = np.asarray(z[f"{pdb}/sd"], float)
        X, names = DF.build(seq, prob, exp, sd, centres)
        dtrue = np.asarray(z[f"{pdb}/dtrue"], float)
        out[f"{pdb}/X"] = X
        out[f"{pdb}/y"] = (dtrue > exp).astype(np.int8)
        out[f"{pdb}/err"] = (exp - dtrue).astype(np.float32)
        out[f"{pdb}/exp"] = exp.astype(np.float32)
        if (k + 1) % 100 == 0:
            print(f"  feats {k+1}/{len(keys)} [{time.time()-t0:.0f}s free={I.free_gb():.1f}]", flush=True)
    out["keys"] = np.array(["\t".join(map(str, r)) for r in keys], dtype=object)
    out["names"] = np.array(names, dtype=object)
    np.savez_compressed(FEATS, **out)
    print("wrote", FEATS, flush=True)


def load_bank():
    z = np.load(FEATS, allow_pickle=True)
    keys = [str(k).split("\t") for k in z["keys"]]
    names = [str(s) for s in z["names"]]
    return z, keys, names


# --------------------------------------------------------------------------- models
def fit_predict(Xtr, ytr, Xte, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06, max_depth=6,
                                       min_samples_leaf=60, l2_regularization=1.0,
                                       early_stopping=True, validation_fraction=0.12,
                                       random_state=seed)
    m.fit(Xtr, ytr)
    return m.predict_proba(Xte)[:, 1], m


def main():
    build_featbank()
    z, keys, names = load_bank()
    sepcols = [names.index(c) for c in DF.SEP_COLS]
    tgset = {t["pdb"] for t in I.targets()}

    pdbs = [k[0] for k in keys]; seqs = [k[1] for k in keys]
    ns = np.array([int(k[2]) for k in keys]); folds = np.array([int(k[3]) for k in keys])
    Xs = {p: np.asarray(z[f"{p}/X"], np.float32) for p in pdbs}
    ys = {p: np.asarray(z[f"{p}/y"], np.int8) for p in pdbs}
    npairs = np.array([len(ys[p]) for p in pdbs])
    print(f"corpus: {len(pdbs)} peptides, {npairs.sum()} pairs, "
          f"base rate P(under-predicted) = {np.mean(np.concatenate([ys[p] for p in pdbs])):.4f}",
          flush=True)

    arms = ("full", "sep", "shuf")
    prob = {a: {} for a in arms}
    const = {}; shellmaj = {}
    rng = np.random.default_rng(0)
    for f in range(5):
        tr = [k for k, p in enumerate(pdbs) if folds[k] != f]
        te = [k for k, p in enumerate(pdbs) if folds[k] == f]
        Xtr = np.vstack([Xs[pdbs[k]] for k in tr]); ytr = np.concatenate([ys[pdbs[k]] for k in tr])
        Xte = np.vstack([Xs[pdbs[k]] for k in te])
        cuts = np.cumsum([len(ys[pdbs[k]]) for k in te])[:-1]
        print(f"  fold {f}: train {Xtr.shape} test {Xte.shape} free={I.free_gb():.1f}", flush=True)

        p_full, m = fit_predict(Xtr, ytr, Xte)
        p_sep, _ = fit_predict(Xtr[:, sepcols], ytr, Xte[:, sepcols])
        p_shuf, _ = fit_predict(Xtr, rng.permutation(ytr), Xte)
        for a, pv in (("full", p_full), ("sep", p_sep), ("shuf", p_shuf)):
            for kk, part in zip(te, np.split(pv, cuts)):
                prob[a][pdbs[kk]] = part.astype(np.float32)
        # constant / per-shell majority from the training folds only
        septr = Xtr[:, names.index("sep")]
        cmaj = float(ytr.mean() > 0.5)
        sh = {}
        for s in np.unique(septr):
            m_ = septr == s
            sh[float(s)] = float(ytr[m_].mean() > 0.5) if m_.sum() > 30 else cmaj
        for kk in te:
            p = pdbs[kk]
            const[p] = np.full(len(ys[p]), cmaj, np.float32)
            sept = Xs[p][:, names.index("sep")]
            shellmaj[p] = np.array([sh.get(float(s), cmaj) for s in sept], np.float32)
        del Xtr, Xte

    prob["const"] = const; prob["shellmaj"] = shellmaj
    # --------------------------------------------------------------- metrics
    from sklearn.metrics import roc_auc_score
    def acc(p_, pdbsel):
        c = np.concatenate([(prob[p_][q] > 0.5).astype(int) == ys[q] for q in pdbsel])
        return float(c.mean())
    def auc(p_, pdbsel):
        yy = np.concatenate([ys[q] for q in pdbsel]); pp = np.concatenate([prob[p_][q] for q in pdbsel])
        return float(roc_auc_score(yy, pp)) if len(np.unique(yy)) > 1 else float("nan")

    groups = {"corpus787": pdbs,
              "tuning126": [p for p in pdbs if p in tgset],
              "fail18": [p for p in pdbs if p in set(I.FAIL18)],
              "other108": [p for p in pdbs if p in tgset and p not in set(I.FAIL18)],
              "nontuning": [p for p in pdbs if p not in tgset]}
    res = {"base_rate": float(np.mean(np.concatenate([ys[p] for p in pdbs]))),
           "n_pairs": int(npairs.sum()), "n_peptides": len(pdbs), "groups": {}}
    print(f"\n{'group':12s}{'n_pairs':>9s}" + "".join(f"{a:>10s}" for a in
          ("full", "sep", "shellmaj", "const", "shuf")) + f"{'AUCfull':>9s}{'AUCsep':>9s}")
    for g, sel in groups.items():
        row = {"n_pairs": int(sum(len(ys[q]) for q in sel)),
               "base_rate": float(np.mean(np.concatenate([ys[q] for q in sel])))}
        for a in ("full", "sep", "shellmaj", "const", "shuf"):
            row[f"acc_{a}"] = acc(a, sel)
        row["auc_full"] = auc("full", sel); row["auc_sep"] = auc("sep", sel)
        row["auc_shuf"] = auc("shuf", sel)
        # per separation shell
        sh = {}
        sept = np.concatenate([Xs[q][:, names.index("sep")] for q in sel])
        yy = np.concatenate([ys[q] for q in sel])
        for a in ("full", "sep"):
            pp = np.concatenate([prob[a][q] for q in sel])
            sh[a] = {int(s): float((((pp > 0.5).astype(int) == yy)[sept == s]).mean())
                     for s in np.unique(sept)}
        row["shell_acc"] = sh
        row["shell_n"] = {int(s): int((sept == s).sum()) for s in np.unique(sept)}
        res["groups"][g] = row
        print(f"{g:12s}{row['n_pairs']:9d}" + "".join(f"{row['acc_'+a]:10.4f}" for a in
              ("full", "sep", "shellmaj", "const", "shuf")) +
              f"{row['auc_full']:9.4f}{row['auc_sep']:9.4f}")

    # per-target accuracy on the 126 (for the emit stage's concentration analysis)
    res["per_target"] = {q: {a: float(((prob[a][q] > 0.5).astype(int) == ys[q]).mean())
                             for a in ("full", "sep", "shellmaj", "const", "shuf")}
                         for q in groups["tuning126"]}
    I.write("dir_head", res)
    np.savez_compressed(PRED, **{f"{q}/{a}": prob[a][q] for q in groups["tuning126"]
                                 for a in ("full", "sep", "shellmaj", "const", "shuf")})
    print("wrote", PRED, flush=True)


if __name__ == "__main__":
    main()
