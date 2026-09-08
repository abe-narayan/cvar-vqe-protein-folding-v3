"""E2 -- DEPLOYABLE piece selection: can a 4-8-residue piece near the native sub-trace be found
in the top-k of a deployable key?  Keys (all computable without the native):
    blosum : BLOSUM62 sum between the piece's sequence and the target sub-sequence
    esm    : -sum_t ||pca128(target, s+t) - pca128(parent, start+t)||^2
    lfo    : sum_t log p(state of piece residue t | target residue s+t), an LFO per-residue
             torsion-state classifier (8 circular k-means states) trained on the fold's library
    tors   : -sum_t ||(cos,sin)(phi,psi)_pred(s+t) - (cos,sin)(phi,psi)_piece(t)||^2 (regression head)
    combo  : sum of z-scores of blosum + esm + lfo + tors
    random : control
    ORACLE : local ideal-CA RMSD (the E1 floor, diagnostic only)
Metric: for each (target, s, L) the min local RMSD among the top-k pieces, and hit@k(tau).
"""
from __future__ import annotations
import os, sys, json, time, math
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_bank as AB
from s12 import assembly_common as AC
from s7 import audit

B62 = audit.B62
NS = AB.NSTATE
WIN = 2          # +-2 one-hot neighbours
LS = [4, 5, 6, 7, 8]
KS = [5, 20, 50]
TAUS = [0.5, 1.0]


def onehot_context(code, chain, pos, clen):
    """(R, (2*WIN+1)*20 + 2) one-hot of residues i-WIN..i+WIN within the same chain, + terminal flags."""
    R = len(code); F = np.zeros((R, (2 * WIN + 1) * 20 + 2), np.float32)
    for d in range(-WIN, WIN + 1):
        j = np.arange(R) + d
        ok = (j >= 0) & (j < R)
        ok[ok] &= chain[j[ok]] == chain[ok]
        F[np.where(ok)[0], (d + WIN) * 20 + code[j[ok]].astype(int)] = 1.0
    F[:, -2] = (pos == 0); F[:, -1] = (pos == clen - 1)
    return F


def lib_features(fold):
    R = AB.residues(fold)
    X = np.concatenate([R["esm"].astype(np.float32),
                        onehot_context(R["code"], R["chain"], R["pos"], R["clen"])], 1)
    ang = np.stack([np.cos(R["phi"]), np.sin(R["phi"]), np.cos(R["psi"]), np.sin(R["psi"])], 1).astype(np.float32)
    return X, R["state"].astype(np.int64), ang, R["chain"], R["valid_phi"] & R["valid_psi"]


def target_features(seq, fold):
    from s12 import esm_bank
    e = esm_bank.load()[seq][1].astype(np.float32)
    n = len(seq); code = np.array([I.ALPHABET.index(c) for c in seq], np.int8)
    return np.concatenate([e, onehot_context(code, np.zeros(n, int), np.arange(n), np.full(n, n))], 1)


class LFO:
    def __init__(self, d_in, seed=0):
        import torch, torch.nn as nn
        torch.manual_seed(seed); torch.set_num_threads(1)   # 2 threads spin-wait under contention: 25x slower
        self.torch = torch
        self.net = nn.Sequential(nn.Linear(d_in, 256), nn.ReLU(), nn.Dropout(0.2), nn.Linear(256, 128), nn.ReLU(),
                                 nn.Dropout(0.2), nn.Linear(128, NS + 4))

    def fit(self, X, y, ang, epochs=15, lr=1e-3):
        torch = self.torch
        X = torch.tensor(X); y = torch.tensor(y); ang = torch.tensor(ang)
        opt = torch.optim.Adam(self.net.parameters(), lr=lr, weight_decay=1e-5)
        n = len(X); g = torch.Generator().manual_seed(0)
        self.net.train()
        for ep in range(epochs):
            perm = torch.randperm(n, generator=g)
            for a in range(0, n, 1024):
                idx = perm[a:a + 1024]
                out = self.net(X[idx])
                loss = torch.nn.functional.cross_entropy(out[:, :NS], y[idx]) + ((out[:, NS:] - ang[idx]) ** 2).sum(1).mean()
                opt.zero_grad(); loss.backward(); opt.step()
        self.net.eval()

    def predict(self, X):
        torch = self.torch
        with torch.no_grad():
            out = self.net(torch.tensor(np.asarray(X, np.float32)))
            logp = torch.log_softmax(out[:, :NS], 1).numpy()
            ang = out[:, NS:].numpy()
        # renormalise the regression head onto the unit circles
        ang = ang.copy()
        for a in (0, 2):
            nrm = np.sqrt((ang[:, a:a + 2] ** 2).sum(1, keepdims=True)) + 1e-9
            ang[:, a:a + 2] /= nrm
        return logp, ang


_LFO = {}


def lfo_model(fold, pep_only=False):
    """Train (or load) the fold's LFO classifier; returns (model, info).
    pep_only=True trains on the out-of-fold PEPTIDE chains only (the target distribution;
    S7-2 showed protein fragments shift the distogram, the same may hold for local states)."""
    key = (fold, pep_only)
    if key in _LFO:
        return _LFO[key]
    import torch
    path = os.path.join(I.CACHE, f"asm_lfo_f{fold}{'_pep' if pep_only else ''}.pt")
    X, y, ang, chain, ok = lib_features(fold)
    if pep_only:
        ok = ok & AB.residues(fold)["is_pep"]
    m = LFO(X.shape[1])
    if os.path.exists(path):
        st = torch.load(path); m.net.load_state_dict(st["state"]); m.net.eval(); info = st["info"]
    else:
        rng = np.random.default_rng(fold)
        chains = np.unique(chain[ok]); val_ch = rng.permutation(chains)[: len(chains) // 10]
        isval = np.isin(chain, val_ch)
        tr = ok & ~isval; va = ok & isval
        t0 = time.time(); m.fit(X[tr], y[tr], ang[tr])
        lp, _ = m.predict(X[va]); pred = lp.argmax(1)
        maj = np.bincount(y[tr], minlength=NS).argmax()
        info = {"val_acc": float((pred == y[va]).mean()), "val_top2": float(np.mean([y[va][i] in np.argsort(-lp[i])[:2] for i in range(len(pred))])),
                "majority_acc": float((y[va] == maj).mean()), "n_train": int(tr.sum()), "n_val": int(va.sum()), "seconds": time.time() - t0}
        # then refit on everything valid for deployment
        m = LFO(X.shape[1]); m.fit(X[ok], y[ok], ang[ok])
        torch.save({"state": m.net.state_dict(), "info": info}, path)
    _LFO[key] = (m, info)
    return _LFO[key]


def native_state_accuracy(fold, tg, pep_only=False):
    """ORACLE labels: LFO accuracy on the tuning targets' native torsions (interior residues)."""
    import peptide_db as db
    m, _ = lfo_model(fold, pep_only); R = AB.residues(fold); C = R["centres"]
    acc, top2, nres, maj = [], [], 0, np.bincount(R["state"], minlength=NS).argmax()
    accm = []
    for t in tg:
        if t["fold"] != fold:
            continue
        q = db.by_pdb(t["pdb"]); st = AB.assign_state(q.phi.astype(float), q.psi.astype(float), C)
        lp, _ = m.predict(target_features(t["seq"], fold))
        pred = lp.argmax(1)
        sl = slice(1, t["n"] - 1)
        acc.append((pred[sl] == st[sl]).mean()); accm.append((st[sl] == maj).mean())
        top2.append(np.mean([st[i] in np.argsort(-lp[i])[:2] for i in range(1, t["n"] - 1)]))
    return {"acc": float(np.mean(acc)), "top2": float(np.mean(top2)), "majority": float(np.mean(accm)), "n_targets": len(acc)}


def target_keys(t):
    """Per-target deployable tables: ESM residue-similarity matrix S (n, R), LFO log-probs (n, NS), tors (n, 4)."""
    fold, seq, n = t["fold"], t["seq"], t["n"]
    R = AB.residues(fold)
    from s12 import esm_bank
    e = esm_bank.load()[seq][1].astype(np.float32)
    E = R["esm"].astype(np.float32)
    S = -(((e ** 2).sum(1)[:, None] + (E ** 2).sum(1)[None, :]) - 2.0 * e @ E.T)     # (n, R) negative sq dist
    m, _ = lfo_model(fold)
    logp, ang = m.predict(target_features(seq, fold))
    mp, _ = lfo_model(fold, pep_only=True)
    logp_p, ang_p = mp.predict(target_features(seq, fold))
    tcode = np.array([I.ALPHABET.index(c) for c in seq], int)
    return {"S": S, "logp": logp, "ang": ang, "logp_pep": logp_p, "ang_pep": ang_p, "tcode": tcode}


def piece_scores(t, K, s, L, rng=None, only=None):
    """dict key -> score over all bank pieces (fold, L) at target position s (higher = better).
    only: optional single deployable key name -> compute just what that key needs (E3 speed path)."""
    fold = t["fold"]; R = AB.residues(fold); P = AB.pieces(fold, L)
    idx = P["start"][:, None] + np.arange(L)[None]           # (Pn, L) residue rows
    codes = R["code"][idx].astype(int)
    out = {}
    z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
    if only is not None:
        if only == "random":
            # deployable "diversity" shortlist: a fixed random draw, seeded per (target, s, L)
            g = np.random.default_rng(abs(hash((t["pdb"], s, L))) % (2 ** 32))
            return {"random": g.random(len(codes))}
        if only in ("blosum", "combo", "combo_pep", "blosum_esm"):
            out["blosum"] = B62[codes, K["tcode"][s:s + L][None, :]].sum(1)
        if only in ("esm", "combo", "combo_pep", "blosum_esm", "combo_noblosum"):
            out["esm"] = K["S"][np.arange(s, s + L)[None, :], idx].sum(1)
        if only not in ("blosum", "esm", "blosum_esm"):
            st = R["state"][idx].astype(int)
            A = np.stack([np.cos(R["phi"][idx]), np.sin(R["phi"][idx]), np.cos(R["psi"][idx]), np.sin(R["psi"][idx])], -1)
            if only in ("lfo", "combo", "combo_noblosum"):
                out["lfo"] = K["logp"][np.arange(s, s + L)[None, :], st].sum(1)
            if only in ("tors", "combo", "combo_noblosum"):
                out["tors"] = -((A - K["ang"][s:s + L][None]) ** 2).sum((1, 2))
            if only in ("lfo_pep", "combo_pep"):
                out["lfo_pep"] = K["logp_pep"][np.arange(s, s + L)[None, :], st].sum(1)
            if only in ("tors_pep", "combo_pep"):
                out["tors_pep"] = -((A - K["ang_pep"][s:s + L][None]) ** 2).sum((1, 2))
        if only == "combo":
            out["combo"] = z(out["blosum"]) + z(out["esm"]) + z(out["lfo"]) + z(out["tors"])
        elif only == "combo_noblosum":
            out["combo_noblosum"] = z(out["esm"]) + z(out["lfo"]) + z(out["tors"])
        elif only == "combo_pep":
            out["combo_pep"] = z(out["blosum"]) + z(out["esm"]) + z(out["lfo_pep"]) + z(out["tors_pep"])
        elif only == "blosum_esm":
            out["blosum_esm"] = z(out["blosum"]) + z(out["esm"])
        return out
    out["blosum"] = B62[codes, K["tcode"][s:s + L][None, :]].sum(1)
    out["esm"] = K["S"][np.arange(s, s + L)[None, :], idx].sum(1)
    st = R["state"][idx].astype(int)
    out["lfo"] = K["logp"][np.arange(s, s + L)[None, :], st].sum(1)
    A = np.stack([np.cos(R["phi"][idx]), np.sin(R["phi"][idx]), np.cos(R["psi"][idx]), np.sin(R["psi"][idx])], -1)
    out["tors"] = -((A - K["ang"][s:s + L][None]) ** 2).sum((1, 2))
    out["lfo_pep"] = K["logp_pep"][np.arange(s, s + L)[None, :], st].sum(1)
    out["tors_pep"] = -((A - K["ang_pep"][s:s + L][None]) ** 2).sum((1, 2))
    out["combo"] = z(out["blosum"]) + z(out["esm"]) + z(out["lfo"]) + z(out["tors"])
    out["combo_noblosum"] = z(out["esm"]) + z(out["lfo"]) + z(out["tors"])
    out["combo_pep"] = z(out["blosum"]) + z(out["esm"]) + z(out["lfo_pep"]) + z(out["tors_pep"])
    out["blosum_esm"] = z(out["blosum"]) + z(out["esm"])
    if rng is not None:
        out["random"] = rng.random(len(codes))
    return out


KEYS = ["blosum", "esm", "lfo", "tors", "lfo_pep", "tors_pep", "combo", "combo_noblosum", "combo_pep", "blosum_esm", "random", "ORACLE"]


def run_target(t):
    pdb, n, fold = t["pdb"], t["n"], t["fold"]
    out_path = os.path.join(I.CACHE, f"asm_e2_{pdb}.json")
    if os.path.exists(out_path):
        return json.load(open(out_path))
    u = I.load_univ(pdb); nat = u["nat_ca"]
    K = target_keys(t); rng = np.random.default_rng(0)
    res = {"pdb": pdb, "n": n, "fold": fold, "rows": []}
    for L in LS:
        for s in range(0, n - L + 1):
            r = AC.local_rmsd_all(fold, s, L, nat[s:s + L])
            sc = piece_scores(t, K, s, L, rng); sc["ORACLE"] = -r
            row = {"L": L, "s": s, "base_rate": {str(tau): float((r <= tau).mean()) for tau in TAUS}, "oracle_top1": float(r.min())}
            for key in KEYS:
                top = np.argsort(-sc[key], kind="stable")[:max(KS)]
                rt = r[top]
                row[key] = {str(k): float(rt[:k].min()) for k in KS}
            res["rows"].append(row)
    json.dump(res, open(out_path, "w"))
    return res


def aggregate(rows_by_target, tg):
    pdbs = [r["pdb"] for r in rows_by_target]
    fail = np.array([p in I.FAIL18 for p in pdbs])
    agg = {"keys": KEYS, "L": LS, "k": KS}
    tab = {}
    for L in LS:
        tab[L] = {}
        # per-target mean over positions of min-local-RMSD in top-k, and hit rates
        for key in KEYS:
            for k in KS:
                per_t = np.array([np.mean([row[key][str(k)] for row in r["rows"] if row["L"] == L]) for r in rows_by_target])
                hits = {str(tau): np.array([np.mean([row[key][str(k)] <= tau for row in r["rows"] if row["L"] == L]) for r in rows_by_target]) for tau in TAUS}
                tab[L][f"{key}@{k}"] = {"min_rmsd": AC.group_means(per_t, pdbs),
                                        **{f"hit{tau}": AC.group_means(hits[str(tau)], pdbs) for tau in TAUS}}
        for tau in TAUS:
            br = np.array([np.mean([row["base_rate"][str(tau)] for row in r["rows"] if row["L"] == L]) for r in rows_by_target])
            tab[L][f"base_rate{tau}"] = AC.group_means(br, pdbs)
            for k in KS:
                exp = np.array([np.mean([1 - (1 - row["base_rate"][str(tau)]) ** k for row in r["rows"] if row["L"] == L]) for r in rows_by_target])
                tab[L][f"random_expected_hit{tau}@{k}"] = AC.group_means(exp, pdbs)
    agg["table"] = tab
    # paired: combo vs blosum vs random at L=6, k=20 (min rmsd per target)
    folds = [r["fold"] for r in rows_by_target]
    def per_t(key, L, k):
        return np.array([np.mean([row[key][str(k)] for row in r["rows"] if row["L"] == L]) for r in rows_by_target])
    agg["paired"] = {}
    for L in (5, 6, 8):
        for k in (20,):
            for a, b in (("combo", "blosum"), ("combo", "random"), ("blosum", "random"), ("lfo", "blosum"), ("esm", "blosum"), ("tors", "blosum"), ("combo_noblosum", "combo"),
                         ("lfo_pep", "lfo"), ("tors_pep", "tors"), ("combo_pep", "combo"), ("blosum_esm", "blosum"), ("esm", "random"), ("lfo_pep", "random")):
                agg["paired"][f"L{L}_k{k}_{a}_vs_{b}"] = I.paired(per_t(a, L, k), per_t(b, L, k), folds=folds, names=pdbs)
    return agg


if __name__ == "__main__":
    tg = I.targets()
    print("free GB", I.free_gb(), flush=True)
    lfo_info = {}
    for f in range(5):
        lfo_info[f] = {}
        for pep in (False, True):
            m, info = lfo_model(f, pep)
            nat = native_state_accuracy(f, tg, pep)
            lfo_info[f]["pep" if pep else "lib"] = {"library_holdout": info, "target_native": nat}
            print(f"LFO fold {f} {'PEP' if pep else 'LIB'}: holdout acc {info['val_acc']:.3f} (top2 {info['val_top2']:.3f}, majority {info['majority_acc']:.3f}) | "
                  f"target native acc {nat['acc']:.3f} (top2 {nat['top2']:.3f}, majority {nat['majority']:.3f}) [{info['seconds']:.0f}s]", flush=True)
    rows = []
    t0 = time.time()
    import multiprocessing as mp
    with mp.Pool(2) as pool:
        for k, r in enumerate(pool.imap_unordered(run_target, tg)):
            rows.append(r)
            if k % 10 == 0:
                print(f"[{k+1}/126] {r['pdb']} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
    rows.sort(key=lambda r: r["pdb"])
    agg = aggregate(rows, tg); agg["lfo"] = lfo_info
    print(I.write("assembly_e2_selection", agg))
    for L in LS:
        print(f"L={L}: base0.5={agg['table'][L]['base_rate0.5']['all']:.3f} base1.0={agg['table'][L]['base_rate1.0']['all']:.3f}")
        for key in KEYS:
            d = agg["table"][L][f"{key}@20"]
            print(f"   {key:15s} @20 min_rmsd all/f18/o108 = {d['min_rmsd']['all']:.3f}/{d['min_rmsd']['fail18']:.3f}/{d['min_rmsd']['other108']:.3f}  hit0.5={d['hit0.5']['all']:.3f} hit1.0={d['hit1.0']['all']:.3f}")
