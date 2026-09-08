"""forensics Part 2d: is the target's SIZE predictable, and does a predicted size work as a key?

`forensics_compact` found that an ORACLE size gate on RETRIEVAL (keep universe windows with
|rg - rg_native| < 0.5, then BLOSUM top-500) is worth -0.314 A overall and -1.108 A on FAIL18,
while the same gate applied INSIDE the shipped pool is worth -0.071 A.  The shipped distogram's
implied rg tracks the native rg with slope 0.238 (r = 0.363) -- almost no information.

Stage 1: train a leave-fold-out regressor for the native rg from the TARGET SEQUENCE only
  (ESM-2 pca128 mean/max/first/last pooling + length + composition), labels = library peptides'
  own rg (out-of-fold; these are library structures, a legitimate training label).
Stage 2: run the same retrieval gate with the PREDICTED rg and a widened window, through the
  full chain.  Control: the shipped distogram's implied rg as the gate; and a constant-rg gate
  (predict the training mean) to separate "a size gate helps at all" from "this size is right".
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import pearsonr, spearmanr
from s12 import instrument as I
from s12 import forensics_lib as L

def seq_features(seq, bank):
    E = bank[seq][1]
    comp = np.array([seq.count(c) for c in I.ALPHABET], float) / len(seq)
    return np.concatenate([E.mean(0), E.max(0), E[0], E[-1], [len(seq), 1.0 / len(seq)], comp])

def build():
    from s12 import esm_bank
    import peptide_db as db
    bank = esm_bank.load(); P = L.parents()
    X, y, f, s = [], [], [], []
    folds = P["folds"]
    for k in range(P["n_pep"]):
        q = P["objs"][k]
        X.append(seq_features(q.seq, bank)); y.append(L.rg(q.ca)); f.append(folds[q.seq]); s.append(q.seq)
    return np.array(X, np.float32), np.array(y), np.array(f), s

def main():
    from sklearn.linear_model import RidgeCV
    from sklearn.ensemble import GradientBoostingRegressor
    from sklearn.preprocessing import StandardScaler
    from s12 import esm_bank
    X, y, f, seqs = build()
    tg = I.targets(); bank = esm_bank.load()
    Xt = np.array([seq_features(t["seq"], bank) for t in tg], np.float32)
    ft = np.array([t["fold"] for t in tg])
    yt = np.array([L.rg(I.load_univ(t["pdb"])["nat_ca"]) for t in tg])       # ORACLE label, evaluation only
    tseqs = {t["seq"] for t in tg}
    pred_r, pred_g, pred_c = np.zeros(len(tg)), np.zeros(len(tg)), np.zeros(len(tg))
    for fold in range(5):
        tr = (f != fold) & ~np.isin(seqs, list(tseqs))     # never train on any tuning target's sequence
        te = ft == fold
        sc = StandardScaler().fit(X[tr]); A = sc.transform(X[tr])
        r = RidgeCV(alphas=np.logspace(-1, 4, 20)).fit(A, y[tr])
        g = GradientBoostingRegressor(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, random_state=0).fit(A, y[tr])
        pred_r[te] = r.predict(sc.transform(Xt[te])); pred_g[te] = g.predict(sc.transform(Xt[te])); pred_c[te] = y[tr].mean()
    disto = np.array([json.load(open(os.path.join(I.RESULTS, "forensics_compact.json")))["rows"][k]["rg_pred"] for k in range(len(tg))])
    out = {"n_train_peptides": int((~np.isin(seqs, list(tseqs))).sum()), "leakage": "training excludes every tuning-target sequence and every in-fold peptide"}
    for name, p in [("ridge", pred_r), ("gbt", pred_g), ("constant", pred_c), ("distogram", disto)]:
        F = np.array([t["pdb"] in I.FAIL18 for t in tg])
        out[name] = dict(pearson=float(pearsonr(p, yt)[0]), spearman=float(spearmanr(p, yt).correlation),
                         slope=float(np.polyfit(yt, p, 1)[0]), mae=float(np.abs(p - yt).mean()),
                         mae_fail18=float(np.abs(p - yt)[F].mean()), mae_other=float(np.abs(p - yt)[~F].mean()),
                         within_0_5=float((np.abs(p - yt) < 0.5).mean()), within_1_0=float((np.abs(p - yt) < 1.0).mean()))
    np.savez(os.path.join(I.CACHE, "forensics_rgpred.npz"), ridge=pred_r, gbt=pred_g, const=pred_c, disto=disto, o_true=yt,
             pdb=np.array([t["pdb"] for t in tg], dtype=object))
    I.write("forensics_rgpred_stage1", out)
    for k, v in out.items():
        if isinstance(v, dict):
            print(f"{k:11s} r={v['pearson']:+.3f} rho={v['spearman']:+.3f} slope={v['slope']:.3f} MAE={v['mae']:.3f} (F18 {v['mae_fail18']:.3f} / O108 {v['mae_other']:.3f}) within0.5={v['within_0_5']:.2f} within1.0={v['within_1_0']:.2f}")

if __name__ == "__main__":
    main()
