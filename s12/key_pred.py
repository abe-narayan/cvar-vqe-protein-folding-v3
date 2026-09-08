"""The torsion-bin PREDICTOR: sequence -> per-residue posterior over torsion bins.

Leave-fold-out by construction: model for fold f is trained on `key_corpus.corpus(f)`,
which is the legal library for every fold-f target (out-of-fold peptides + this fold's
fragments).  The predictor NEVER sees the pool, the native, or the target's own fold.

Variants (all trained identically, only the data/labels/features change):
    full      peptide + fragment residues, real labels          <- the deployable one
    pep       peptide residues only                              (S7-2 distribution-shift check)
    frag      fragment residues only
    shuf      real data, labels randomly permuted                <- mandatory null
    comp      composition-only features (position independent)   <- composition null
    esm       full + ESM-2 PCA features for the context centre   (optional)
"""
from __future__ import annotations
import os, sys, json, math, time, argparse
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import torch
torch.set_num_threads(2)
import torch.nn as nn

from s12 import instrument as I
from s12 import key_lib as KL
from s12 import key_corpus as KC

W = 7                       # +-W residue context
PAD = 20                    # pad symbol
NAA = 21
PROPS = np.zeros((NAA, 5), np.float32)
from priors import _PROPS as _PP, _DEFAULT as _PD
for a, v in _PP.items():
    PROPS[I.ALPHABET.index(a)] = v
PROPS[PAD] = _PD
PROPS = (PROPS - PROPS.mean(0)) / (PROPS.std(0) + 1e-6)


# ------------------------------------------------------------------------- featurisation
def context_feats(codes, pid=None):
    """(L, (2W+1)*(NAA+5)) sliding-context features.  Context does not cross parents."""
    codes = np.asarray(codes, np.int64)
    L = len(codes)
    if pid is None:
        pid = np.zeros(L, np.int64)
    pid = np.asarray(pid, np.int64)
    idx = np.arange(L)[:, None] + np.arange(-W, W + 1)[None, :]
    ok = (idx >= 0) & (idx < L)
    idxc = np.clip(idx, 0, L - 1)
    ok &= pid[idxc] == pid[:, None]
    ctx = np.where(ok, codes[idxc], PAD)
    oh = np.zeros((L, 2 * W + 1, NAA), np.float32)
    np.put_along_axis(oh, ctx[:, :, None], 1.0, 2)
    pr = PROPS[ctx]
    return np.concatenate([oh, pr], 2).reshape(L, -1)


def comp_feats(codes, pid=None):
    """Composition-only: the parent's aa composition, repeated at every position."""
    codes = np.asarray(codes, np.int64); L = len(codes)
    pid = np.zeros(L, np.int64) if pid is None else np.asarray(pid, np.int64)
    out = np.zeros((L, 20 + 5), np.float32)
    for p in np.unique(pid):
        m = pid == p
        c = np.bincount(codes[m], minlength=NAA)[:20].astype(np.float32)
        c = c / max(c.sum(), 1.0)
        out[m, :20] = c
        out[m, 20:] = (PROPS[codes[m]].mean(0))
    return out


def esm_matrix_corpus(fold):
    """(L, 32) ESM-2 PCA-32 rows aligned to corpus(fold)'s residue order.
    Parents are looked up in the ESM bank by their exact sequence (99.6 % hit rate)."""
    path = os.path.join(I.CACHE, f"key_esm_f{fold}.npy")
    if os.path.exists(path):
        return np.load(path).astype(np.float32)
    from s12 import esm_bank
    b = esm_bank.load()
    c = KC.corpus(fold)
    pid = c["pid"]; codes = c["codes"]
    E = np.zeros((len(codes), 32), np.float32)
    miss = 0
    for p in np.unique(pid):
        m = pid == p
        s = I.decode(codes[m])
        v = b.get(s)
        if v is None:
            miss += 1
            continue
        e = np.asarray(v[0], np.float32)
        if len(e) == m.sum():
            E[m] = e
        else:
            miss += 1
    print(f"  esm fold {fold}: {miss} parents missing", flush=True)
    np.save(path, E.astype(np.float16))
    return E


def esm_feats(codes, pid=None, E=None):
    """context features + centre-residue ESM PCA-32 + context-mean ESM PCA-32."""
    X = context_feats(codes, pid)
    if E is None:
        E = np.zeros((len(codes), 32), np.float32)
    E = np.asarray(E, np.float32)
    L = len(codes)
    pid = np.zeros(L, np.int64) if pid is None else np.asarray(pid, np.int64)
    idx = np.arange(L)[:, None] + np.arange(-W, W + 1)[None, :]
    ok = (idx >= 0) & (idx < L)
    idxc = np.clip(idx, 0, L - 1)
    ok &= pid[idxc] == pid[:, None]
    Em = (E[idxc] * ok[:, :, None]).sum(1) / np.maximum(ok.sum(1, keepdims=True), 1)
    return np.concatenate([X, E, Em.astype(np.float32)], 1)


FEAT = {"context": context_feats, "comp": comp_feats, "esm": esm_feats}


# ------------------------------------------------------------------------------- model
class MLP(nn.Module):
    def __init__(self, d, nb, h=256):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(0.2),
                                 nn.Linear(h, h // 2), nn.ReLU(), nn.Dropout(0.2),
                                 nn.Linear(h // 2, nb))

    def forward(self, x):
        return self.net(x)


def train_one(X, y, nb, seed=0, epochs=25, bs=512, lr=1e-3, val_frac=0.1, groups=None, verbose=False):
    rng = np.random.default_rng(seed)
    g = np.arange(len(X)) if groups is None else np.asarray(groups)
    ug = np.unique(g)
    vg = set(ug[rng.permutation(len(ug))[: max(1, int(val_frac * len(ug)))]].tolist())
    vm = np.array([x in vg for x in g])
    Xt = torch.tensor(X[~vm]); yt = torch.tensor(y[~vm].astype(np.int64))
    Xv = torch.tensor(X[vm]); yv = torch.tensor(y[vm].astype(np.int64))
    torch.manual_seed(seed)
    m = MLP(X.shape[1], nb)
    opt = torch.optim.Adam(m.parameters(), lr=lr, weight_decay=1e-5)
    lossf = nn.CrossEntropyLoss()
    best, bstate, bad = 1e9, None, 0
    n = len(Xt)
    for ep in range(epochs):
        m.train()
        perm = torch.randperm(n)
        for s in range(0, n, bs):
            b = perm[s:s + bs]
            opt.zero_grad()
            l = lossf(m(Xt[b]), yt[b])
            l.backward(); opt.step()
        m.eval()
        with torch.no_grad():
            lv = float(lossf(m(Xv), yv))
            av = float((m(Xv).argmax(1) == yv).float().mean())
        if verbose:
            print(f"  ep{ep} val_loss={lv:.4f} val_acc={av:.3f}", flush=True)
        if lv < best - 1e-4:
            best, bstate, bad = lv, {k: v.clone() for k, v in m.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= 4:
                break
    if bstate:
        m.load_state_dict(bstate)
    m.eval()
    return m, {"val_loss": best, "val_acc": av}


def predict(m, X):
    with torch.no_grad():
        return torch.softmax(m(torch.tensor(np.asarray(X, np.float32))), 1).numpy()


# ------------------------------------------------------------------------------- driver
VARIANTS = {
    "full":  dict(feat="context", subset="all",  shuffle=False),
    "pep":   dict(feat="context", subset="pep",  shuffle=False),
    "frag":  dict(feat="context", subset="frag", shuffle=False),
    "shuf":  dict(feat="context", subset="all",  shuffle=True),
    "comp":  dict(feat="comp",    subset="all",  shuffle=False),
    "esm":     dict(feat="esm", subset="all", shuffle=False),
    "esmpep":  dict(feat="esm", subset="pep", shuffle=False),
}


def build(alph="abego4", variants=tuple(VARIANTS), seed=0, epochs=25):
    nb = KL.nbins(alph)
    tg = I.targets()
    nat = KL.native_torsions()
    folds = sorted({t["fold"] for t in tg})
    post = {v: {} for v in variants}
    diag = {v: {} for v in variants}
    for f in folds:
        c = KC.corpus(f)
        y_all = KL.bins_of(alph, c["phi"].astype(float), c["psi"].astype(float)).astype(np.int64)
        for v in variants:
            cfg = VARIANTS[v]
            m_sub = {"all": np.ones(len(y_all), bool), "pep": c["org"].astype(bool),
                     "frag": ~c["org"].astype(bool)}[cfg["subset"]]
            if cfg["feat"] == "esm":
                E = esm_matrix_corpus(f)
                X = esm_feats(c["codes"][m_sub], c["pid"][m_sub], E[m_sub])
            else:
                X = FEAT[cfg["feat"]](c["codes"][m_sub], c["pid"][m_sub])
            y = y_all[m_sub].copy()
            if cfg["shuffle"]:
                y = np.random.default_rng(seed + 991).permutation(y)
            t0 = time.time()
            mdl, st = train_one(X, y, nb, seed=seed, epochs=epochs, groups=c["pid"][m_sub])
            diag[v][f] = {**st, "n": int(m_sub.sum()), "secs": round(time.time() - t0, 1)}
            print(f"fold {f} {v}: n={m_sub.sum()} val_acc={st['val_acc']:.3f} "
                  f"({time.time()-t0:.0f}s) free={I.free_gb():.2f}", flush=True)
            for t in tg:
                if t["fold"] != f:
                    continue
                cds = np.array([I.ALPHABET.index(a) if a in I.ALPHABET else PAD for a in t["seq"]], np.int64)
                if cfg["feat"] == "esm":
                    from s12 import esm_bank
                    ev = esm_bank.load().get(t["seq"])
                    Et = np.asarray(ev[0], np.float32) if ev is not None else np.zeros((len(cds), 32), np.float32)
                    Xt = esm_feats(cds, None, Et)
                else:
                    Xt = FEAT[cfg["feat"]](cds)
                post[v][t["pdb"]] = predict(mdl, Xt).astype(np.float32)
            del X, y, mdl
        del c
    return post, diag


def score_predictor(post, alph="abego4"):
    """Per-residue accuracy / confusion vs the NATIVE bins (ORACLE evaluation)."""
    nb = KL.nbins(alph)
    nat = KL.native_torsions()
    out = {}
    for v, d in post.items():
        acc, accf, acco = [], [], []
        cm = np.zeros((nb, nb), int)
        top2 = []
        ent = []
        for pdb, P in d.items():
            tb = KL.bins_of(alph, *nat[pdb]).astype(int)
            pb = np.asarray(P).argmax(1)
            a = (pb == tb).mean()
            acc.append(a); (accf if pdb in I.FAIL18 else acco).append(a)
            for x, y in zip(tb, pb):
                cm[x, y] += 1
            o = np.argsort(-np.asarray(P), 1)[:, :2]
            top2.append(float(np.mean([(tb[k] in o[k]) for k in range(len(tb))])))
            ent.append(float(-(P * np.log(np.clip(P, 1e-9, 1))).sum(1).mean()))
        out[v] = {"acc_all": float(np.mean(acc)), "acc_FAIL18": float(np.mean(accf)),
                  "acc_other108": float(np.mean(acco)), "top2_all": float(np.mean(top2)),
                  "mean_entropy": float(np.mean(ent)), "confusion": cm.tolist(),
                  "per_class_recall": (np.diag(cm) / np.maximum(cm.sum(1), 1)).tolist(),
                  "label_freq": (cm.sum(1) / cm.sum()).tolist()}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--alph", default="abego4")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--out", default="key_pred")
    ap.add_argument("--epochs", type=int, default=25)
    a = ap.parse_args()
    post, diag = build(alph=a.alph, variants=tuple(a.variants.split(",")), epochs=a.epochs)
    sc = score_predictor(post, a.alph)
    print(json.dumps({v: {k: x for k, x in s.items() if k != "confusion"} for v, s in sc.items()}, indent=1))
    I.write(a.out + "_acc", {"alph": a.alph, "score": sc, "diag": diag})
    with open(os.path.join(I.RESULTS, a.out + "_post.json"), "w") as fh:
        json.dump({v: {k: np.asarray(x).tolist() for k, x in d.items()} for v, d in post.items()}, fh)
    print("wrote", a.out)
