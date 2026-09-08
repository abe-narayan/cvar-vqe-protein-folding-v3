"""forensics Part 2c: structural alphabet (Foldseek-in-miniature) retrieval key.

Codebook: k-means (K=8) on (cos phi, sin phi, cos psi, sin psi), fitted on FRAGMENT residues
only (never a target).  Labels: nearest state per residue of every parent.
Predictor: per-residue MLP over ESM-2 pca128 of residues i-2..i+2 + position features,
trained LEAVE-FOLD-OUT (fold f model: peptides of folds != f + fold_fragments(f)), i.e. the
same library the retrieval universe uses.  Also a 3-class CA-SS predictor with labels from
the parents' CA geometry (L.ca_ss), giving an SS key that does not depend on the distogram.
Per target cache: agree_pred (hard alphabet agreement per window), ll_norm (soft),
agree_oracle (ORACLE native alphabet), ss_esm_agree, ss_esm_oracle-free accuracy stats.
"""
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from s12 import instrument as I
from s12 import forensics_lib as L

KST = 8
CB = os.path.join(I.CACHE, "forensics_alphabet_codebook.npz")
MODELS = os.path.join(I.CACHE, "forensics_alphabet_models.pkl")

def tors_feat(phi, psi):
    return np.stack([np.cos(phi), np.sin(phi), np.cos(psi), np.sin(psi)], -1)

def codebook():
    if os.path.exists(CB):
        return np.load(CB)["C"]
    from sklearn.cluster import KMeans
    P = L.parents()
    X = np.concatenate([tors_feat(np.asarray(q.phi, float), np.asarray(q.psi, float)) for q, pep in zip(P["objs"], P["is_pep"]) if not pep])
    X = X[np.isfinite(X).all(1)]
    km = KMeans(KST, n_init=10, random_state=0).fit(X)
    C = km.cluster_centers_
    # order states by (phi, psi) angle for readability
    ang = np.arctan2(C[:, 1], C[:, 0]); C = C[np.argsort(ang)]
    np.savez(CB, C=C, counts=np.bincount(np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), 1), minlength=KST))
    return C

def states(phi, psi, C):
    X = tors_feat(np.asarray(phi, float), np.asarray(psi, float))
    d = ((X[..., None, :] - C) ** 2).sum(-1)
    return np.argmin(d, -1).astype(np.int8)

def residue_features(E, n):
    """E (n,128) -> (n, 5*128+3): neighbours i-2..i+2 zero-padded, rel pos, dist to end, length/20."""
    pad = np.zeros((2, E.shape[1]), E.dtype); Ep = np.vstack([pad, E, pad])
    F = np.concatenate([Ep[k:k + n] for k in range(5)], 1)
    pos = np.arange(n)[:, None]
    extra = np.hstack([pos / max(n - 1, 1), np.minimum(pos, n - 1 - pos) / 8.0, np.full((n, 1), n / 20.0)])
    return np.hstack([F, extra]).astype(np.float32)

def train_fold(fold, C):
    from sklearn.neural_network import MLPClassifier
    from s12 import esm_bank
    from core import data as cdata
    bank = esm_bank.load(); P = parents = L.parents()
    members = [P["objs"][k] for k in range(P["n_pep"]) if P["fold"][k] != fold] + list(cdata.fold_fragments(fold, 5))
    X, ya, ys = [], [], []
    for q in members:
        E = bank[q.seq][1]
        X.append(residue_features(E, q.n)); ya.append(states(q.phi, q.psi, C)); ys.append(L.ca_ss(q.ca))
    X = np.vstack(X); ya = np.concatenate(ya); ys = np.concatenate(ys)
    ma = MLPClassifier((256,), alpha=1e-3, batch_size=256, max_iter=60, early_stopping=True, random_state=0, n_iter_no_change=6).fit(X, ya)
    ms = MLPClassifier((256,), alpha=1e-3, batch_size=256, max_iter=60, early_stopping=True, random_state=0, n_iter_no_change=6).fit(X, ys)
    return ma, ms, dict(n_res=int(len(ya)), n_members=len(members), val_alpha=float(ma.best_validation_score_), val_ss=float(ms.best_validation_score_),
                        prior_alpha=np.bincount(ya, minlength=KST).tolist(), prior_ss=np.bincount(ys, minlength=3).tolist())

def models(C):
    import pickle
    if os.path.exists(MODELS):
        with open(MODELS, "rb") as fh:
            return pickle.load(fh)
    out = {}
    for f in range(5):
        t0 = time.time(); ma, ms, info = train_fold(f, C); out[f] = (ma, ms, info)
        print(f"fold {f}: {info} ({time.time()-t0:.0f}s) free={I.free_gb():.1f}", flush=True)
    with open(MODELS, "wb") as fh:
        pickle.dump(out, fh)
    return out

def per_target(t, C, M):
    from s12 import esm_bank
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    path = os.path.join(I.CACHE, f"forensics_alpha_{pdb}.npz")
    u = I.load_univ(pdb)
    ma, ms, _ = M[fold]
    X = residue_features(esm_bank.load()[seq][1], n)
    pa = ma.predict_proba(X); ps = ms.predict_proba(X)
    pred_a = ma.classes_[pa.argmax(1)]; pred_s = ms.classes_[ps.argmax(1)]
    WA = states(u["PHI"], u["PSI"], C)                      # (nw, n) window alphabet
    cass = L.window_cass(pdb, u)
    logp = np.log(np.maximum(pa, 1e-6))
    colmap = {c: k for k, c in enumerate(ma.classes_)}
    idx = np.vectorize(lambda s: colmap.get(int(s), -1))(WA)
    ll = np.where(idx >= 0, logp[np.arange(n)[None, :], np.maximum(idx, 0)], np.log(1e-6)).sum(1)
    ll_norm = (ll + n * np.log(KST)) / np.log(KST)
    o_na, p = None, None
    import peptide_db as db
    p = db.by_pdb(pdb); o_na = states(p.phi, p.psi, C)      # ORACLE native alphabet (label)
    o_nss = L.ca_ss(u["nat_ca"])
    agree_pred = (WA == pred_a[None, :]).sum(1); agree_oracle = (WA == o_na[None, :]).sum(1)
    ss_esm_agree = (cass == pred_s[None, :]).sum(1)
    np.savez_compressed(path, agree_pred=agree_pred.astype(np.int8), ll_norm=ll_norm.astype(np.float32),
                        agree_oracle=agree_oracle.astype(np.int8), ss_esm_agree=ss_esm_agree.astype(np.int8),
                        pred_a=pred_a.astype(np.int8), pred_s=pred_s.astype(np.int8))
    rr = u["rr"]; thr = float(rr[I.pool_idx(u)].min()) + I.BAND; band = rr <= thr
    prior = np.bincount(WA.ravel(), minlength=KST)
    return dict(pdb=pdb, fail18=pdb in I.FAIL18, n=n,
                acc_alpha=float((pred_a == o_na).mean()), acc_alpha_majority=float((o_na == int(np.argmax(prior))).mean()),
                acc_ss_esm=float((pred_s == o_nss).mean()), acc_ss_disto=float((L.distogram_ss(I.distogram(pdb), n) == o_nss).mean()),
                acc_ss_majority=float((o_nss == np.bincount(cass.ravel(), minlength=3).argmax()).mean()),
                acc_alpha_band_windows=float((WA[band] == o_na[None, :]).mean()) if band.any() else None,
                acc_alpha_univ_windows=float((WA == o_na[None, :]).mean()),
                pred_alpha="".join(str(int(x)) for x in pred_a), o_native_alpha="".join(str(int(x)) for x in o_na),
                pred_ss=L.ss_str(pred_s), o_native_ss=L.ss_str(o_nss),
                mean_ll_band=float(ll_norm[band].mean()) if band.any() else None, mean_ll_univ=float(ll_norm.mean()))

if __name__ == "__main__":
    C = codebook(); print("codebook", np.round(np.degrees(np.arctan2(C[:, 1], C[:, 0])), 0), np.round(np.degrees(np.arctan2(C[:, 3], C[:, 2])), 0), flush=True)
    M = models(C)
    log = os.path.join(I.RESULTS, "forensics_build.log")
    while os.path.exists(log) and "done" not in open(log).read():
        print("waiting for forensics_build to finish", flush=True); time.sleep(30)
    rows = []; t0 = time.time()
    for k, t in enumerate(I.targets()):
        rows.append(per_target(t, C, M))
        print(f"[{k+1:3d}/126] {t['pdb']} acc_alpha={rows[-1]['acc_alpha']:.2f} acc_ss_esm={rows[-1]['acc_ss_esm']:.2f} acc_ss_disto={rows[-1]['acc_ss_disto']:.2f} {time.time()-t0:.0f}s", flush=True)
    F = np.array([r["fail18"] for r in rows])
    agg = {}
    for m in ["acc_alpha", "acc_alpha_majority", "acc_ss_esm", "acc_ss_disto", "acc_ss_majority", "acc_alpha_univ_windows", "mean_ll_univ"]:
        v = np.array([r[m] for r in rows], float); agg[m] = dict(all=float(v.mean()), fail18=float(v[F].mean()), other=float(v[~F].mean()))
    for m in ["acc_alpha_band_windows", "mean_ll_band"]:
        v = np.array([r[m] if r[m] is not None else np.nan for r in rows], float); agg[m] = dict(all=float(np.nanmean(v)), fail18=float(np.nanmean(v[F])), other=float(np.nanmean(v[~F])))
    agg["fold_info"] = {f: M[f][2] for f in M}
    I.write("forensics_alphabet", dict(rows=rows, aggregate=agg, codebook_phi=np.degrees(np.arctan2(C[:, 1], C[:, 0])).tolist(), codebook_psi=np.degrees(np.arctan2(C[:, 3], C[:, 2])).tolist()))
    print(json.dumps(agg, indent=1))
