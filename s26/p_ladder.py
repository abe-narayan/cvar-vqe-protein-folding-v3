"""s26/p_ladder.py -- PROPOSAL C2: THE PRIOR-INPUT LADDER.  Lane P, Sprint 26.

THE QUESTION.  The only steep lever left in the project is the ACCURACY of the distance prior
(S24 L13: -2.15 A per unit gamma at the origin, along the native's own direction only; S25 L12:
a real operator at cos 0.5 travelling 25% is worth +0.024 A).  Every re-reading of the shipped
posterior is closed (S25 L2/L6/L7/L12).  The open item is whether a distogram trained on BETTER
INPUTS is obtainable.  This module retrains the SHIPPED architecture (core.predict.MLP: width 384,
depth 3, 40 epochs, batch 4096, soft-binned CE smooth 0.6, AdamW 2e-3, OneCycle, seed 0, dropout 0)
on the SHIPPED corpus (out-of-fold peptides + core.data.fold_fragments, same order) with the
per-residue block swapped, plus two architecture rungs and one mixture rung, and scores every
rung through the SAME three endpoints against the shipped distogram.

RUNGS
    shipped   the pinned posterior (s12/cache/disto_<pdb>.npz) through this module's own scoring
              path: the reproduction gate for the evaluator (3.4540 / 3.0483 / 3.2148).
    noesm     shipped architecture, ESM block removed: core.predict.features(use_esm=False), 42-d.
    pca32     shipped architecture, shipped inputs rebuilt from s5/esm32.npz (the float32 output
              of esm_features.embed) and s7/repr_cache/esmcon.npz (the exact contact maps).
              REPRODUCTION GATE against distogram_models/fold<f>_esm_frag.pt (`gate`).
    pca32f    as pca32 but the 32-component PCA refit PER FOLD on training sequences only, from
              s5/esmraw.npz.  Isolates "PCA fitting set" from "PCA dimension" for pca128.
    pca128    128-component PCA refit per fold (s7/repr_cache/pca128_fold<f>.npz, S7-11's
              construction) from s5/esmraw.npz.  Lean trainer.
    raw       the full 1280-d ESM-2 650M representation; the MLP's first Linear layer is the
              learned projection.  Lean trainer (residue table gathered per minibatch).
    esm8m     esm2_t6_8M_UR50D (on disk, 30 MB) in place of the 650M model: its 320-d reps ->
              PCA-32 refit per fold, and ITS contact head.  The downward model-size rung; the 3B
              rung is infeasible on this box (s26/results/b1_feasibility.json).
    wide      shipped inputs, MLP width 768 depth 4 (the capacity rung).
    pairnet   the triangle-update PairNet on ESM input: the pinned
              pairnet_models/fold<f>_c64b4_s0.pt (S19 L11 scored it through a distance-geometry
              fit, -0.080 [-0.242,+0.082]; here it goes through the pipeline).  No training.
    mix       (1-lam) * shipped posterior + lam * the K=500 pool's own per-pair 17-bin histogram,
              lam chosen leave-fold-out on the built-chain endpoint.  lam = 0 is the incumbent,
              asserted bit-exactly on the risk table per target.  How this differs from S10-2's
              per-pair Bayes combination (monotone the wrong way) and S19 L14's pool_sep ladder
              (all arms worse than 3.048): those produced a POINT matrix judged by MAE, or moved
              the LOCATION through a distance-geometry fit; this keeps a DISTRIBUTION inside the
              shipped functional and chooses lam on the endpoint, never on the matrix.

ENDPOINTS, per target, every rung, one code path (`endpoint`)
    sel     CA-RMSD of the argmin of the rung's Bayes-risk score over the shipped K=500 pool,
            tie-averaged (ST.argmin_tied).  Selection-side.  Native read AFTER the argmin.
    cloud   top-75 uniform coordinate average in the medoid frame (point cloud, 3.0483 basis).
    arm     s12.instrument.project(cloud, seq, fold)['ca'] (built chain, 3.2148 basis).  PRIMARY.
    fit     the lam=0 chain from the same call ('fit_ca', 3.2041 basis).  Carried, never decided on.
    gam_eff / cos   progress of the rung's posterior toward the native in the S24 MASS
            construction, probability space: <Q-P0, O-P0>/|O-P0|^2 with O = onehot(native bin);
            and in location space: <E_Q-E_P0, D_nat-E_P0>/|D_nat-E_P0|^2 (the s25/loc.py currency).
            NEVER quoted as -2.1496*gam_eff without the cosine beside it (s25 L12).
    mae     mean |E_Q - D_nat|, DIAGNOSTIC ONLY (S7-3: MAE does not price selection).

FALSIFIER (s26/PREREG_C2.md): a rung is a rung only if it beats the shipped distogram on the built
chain by more than its own MDE with the fold-clustered CI excluding zero and 5/5 folds the same
sign.  Anything else is NOT a rung, whatever its MAE does.

MEMORY.  This module never loads esm_cache.npz.  Every ESM input comes from s5/esmraw.npz (487 MB
float32 resident), s5/esm32.npz (11 MB) and s7/repr_cache/esmcon.npz (4 MB), which cover all 6,790
training sequences and all 126 targets (s26/results/p_probe_esm.json, jobrun peak 0.594 GB).
`install_guard()` patches core.data.esm_raw / esm_embed / esm_contacts and esm_features.raw to
those tables, so no code path can reach the bank; an unknown NON-probe sequence is a hard error.

PHASE GATE.  `train`, `probe`, `gate`, `featurise-esm8m` read no native RMSD.  `eval` and
`report` read natives (after selection) and refuse to run unless "PHASE 0 SIGNED OFF" is in
s26/LEDGER.md.

    python s26/p_ladder.py train --rung pca32 --folds 0            (one fold: the time/RAM probe)
    python s26/p_ladder.py gate  --fold 0                          (pca32 vs the pinned fold model)
    python s26/p_ladder.py featurise-esm8m
    python s26/p_ladder.py train --rung <rung> [--seed 1] [--fold-order 4,3,2,1,0]
    python s26/p_ladder.py eval  --rung <rung> [--seed 1]          (after sign-off)
    python s26/p_ladder.py report

Every job checkpoints (one .pt per fold; rows every 10 targets) and is resumable.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import data as D                           # noqa: E402
from core import predict as PR                       # noqa: E402

MODELS = os.path.join(HERE, "models", "p_ladder")
RES = os.path.join(HERE, "results")
os.makedirs(MODELS, exist_ok=True)
os.makedirs(RES, exist_ok=True)

RUNGS = ("shipped", "noesm", "conly", "pca32", "pca32f", "pca128", "raw", "esm8m", "wide", "pairnet", "mix")
TRAINABLE = ("noesm", "conly", "pca32", "pca32f", "pca128", "raw", "esm8m", "wide")
#: conly (PREREG_B2): the physicochemical pair block plus the 13 contact-head columns, NO
#: per-residue embedding block -- the contact head's standalone value as a prior input.
LEAN = ("pca32f", "pca128", "raw", "esm8m")          #: residue-table gather; never materialise X
ARCH = {"wide": dict(width=768, depth=4)}            #: everything else is the shipped 384 x 3
EPOCHS, BATCH, LR, SMOOTH = 40, 4096, 2e-3, 0.6      #: core.predict.MLP.fit defaults, pinned
TOPM, K = 75, 500
NT = int(os.environ.get("NT", "2"))
MIX_LAMS = (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0)
PROBE_SEQ = "ACDEFGHIKL"                             #: core.predict.feature_width's width probe
LEDGER = os.path.join(HERE, "LEDGER.md")
SIGNOFF = "PHASE 0 SIGNED OFF"


def _torch():
    import torch
    torch.set_num_threads(NT)
    return torch


# ============================================================ ESM providers (no bank access)
_T = {}


def _flat_table(path, key="keys", flat="flat", offs="offs", square=False):
    z = np.load(path, allow_pickle=True)
    k, f, o = z[key], z[flat], z[offs]
    out = {}
    for n, s in enumerate(k):
        s = str(s); blk = f[o[n]:o[n + 1]]
        out[s] = blk.reshape(len(s), len(s)) if square else blk
    return out


def esm32():
    """sequence -> (n, 32) float32: the float32 copy of esm_features.embed (s5/esm32.npz)."""
    if "esm32" not in _T:
        _T["esm32"] = _flat_table(os.path.join(ROOT, "s5", "esm32.npz"))
    return _T["esm32"]


def esmraw():
    """sequence -> (n, 1280) float32 (float16 on disk; s5/esmraw.npz).  487 MB resident."""
    if "esmraw" not in _T:
        z = np.load(os.path.join(ROOT, "s5", "esmraw.npz"), allow_pickle=True)
        k, f, o = z["keys"], z["flat"], z["offs"]
        _T["esmraw"] = {str(s): f[o[n]:o[n + 1]].astype(np.float32) for n, s in enumerate(k)}
    return _T["esmraw"]


def contacts():
    """sequence -> (n, n) float32: the exact esm_features.contacts (s7/repr_cache/esmcon.npz)."""
    if "con" not in _T:
        _T["con"] = _flat_table(os.path.join(ROOT, "s7", "repr_cache", "esmcon.npz"), square=True)
    return _T["con"]


def esm8m_tables():
    """sequence -> (reps (n, 320) float32, contacts (n, n) float32) from esm2_t6_8M."""
    if "esm8m" not in _T:
        p = os.path.join(MODELS, "esm8m_raw.npz")
        if not os.path.exists(p):
            raise FileNotFoundError("%s missing: run `featurise-esm8m` first" % p)
        z = np.load(p, allow_pickle=True)
        k, f, o, c, co = z["keys"], z["flat"], z["offs"], z["cflat"], z["coffs"]
        _T["esm8m"] = {str(s): (f[o[n]:o[n + 1]], c[co[n]:co[n + 1]].reshape(len(str(s)), len(str(s))))
                       for n, s in enumerate(k)}
    return _T["esm8m"]


def _pca_apply(raw, z):
    mu, W, sc = np.asarray(z["mu"], np.float64), np.asarray(z["W"], np.float64), np.asarray(z["scale"], np.float64)
    return (((np.asarray(raw, np.float64) - mu) @ W) / np.maximum(sc, 1e-6)).astype(np.float32)


def _pca_fit(X, n_comp, seed=0):
    """Randomised range finder, the s7/repr_select.fit_pca128 construction."""
    X = np.asarray(X, np.float32)
    mu = X.mean(0)
    g = np.random.default_rng(seed).normal(size=(X.shape[1], n_comp + 32)).astype(np.float32)
    Q, _ = np.linalg.qr((X - mu) @ g)
    B = Q.T @ (X - mu)
    _, S, Vt = np.linalg.svd(B, full_matrices=False)
    return {"mu": mu, "W": Vt[:n_comp].T.astype(np.float32),
            "scale": (S[:n_comp] / math.sqrt(len(X))).astype(np.float32)}


def train_entries(fold, n_folds=5):
    """EXACTLY core.predict.train_fold's corpus and order: out-of-fold peptides, then the fold's
    identity-filtered fragments."""
    folds = D.folds(n_folds)
    return [p for p in D.load() if folds[p.seq] != fold] + list(D.fold_fragments(fold, n_folds))


def fold_pca(kind, fold):
    """Per-fold PCA basis for pca32f (from s5/esmraw) and esm8m (from the 8M reps), cached."""
    p = os.path.join(MODELS, "%s_pca_fold%d.npz" % (kind, fold))
    if os.path.exists(p):
        return np.load(p)
    seqs = sorted({e.seq for e in train_entries(fold)})
    if kind == "pca32f":
        src = esmraw(); X = np.concatenate([src[s] for s in seqs], 0)
    elif kind == "esm8m":
        src = esm8m_tables(); X = np.concatenate([src[s][0] for s in seqs], 0)
    else:
        raise ValueError(kind)
    z = _pca_fit(X, 32)
    np.savez(p + ".tmp.npz", **z); os.replace(p + ".tmp.npz", p)
    return np.load(p)


def residue_block(rung, seq, fold):
    """(emb (n, k) float32, contact map (n, n) float32) for one sequence under one rung."""
    if rung in ("pca32", "wide"):
        return esm32()[seq], contacts()[seq]
    if rung == "conly":
        return np.zeros((len(seq), 0), np.float32), contacts()[seq]
    if rung == "pca32f":
        return _pca_apply(esmraw()[seq], fold_pca("pca32f", fold)), contacts()[seq]
    if rung == "pca128":
        z = np.load(os.path.join(ROOT, "s7", "repr_cache", "pca128_fold%d.npz" % fold))
        return _pca_apply(esmraw()[seq], z), contacts()[seq]
    if rung == "raw":
        return esmraw()[seq], contacts()[seq]
    if rung == "esm8m":
        r, c = esm8m_tables()[seq]
        return _pca_apply(r, fold_pca("esm8m", fold)), c
    raise ValueError(rung)


def install_guard():
    """Route every ESM read in core.data / esm_features to the compact tables.  Never the bank.

    The width probe sequence gets zeros of the right shape (as core.pipeline.guard_esm does);
    any other unknown sequence is a hard error, never a silent zero.
    """
    if _T.get("guarded"):
        return
    e32, cn = esm32(), contacts()

    def _raw(seq):
        if seq in cn and seq in esmraw():
            return esmraw()[seq], cn[seq]
        if seq == PROBE_SEQ:
            return np.zeros((len(seq), 1280), np.float32), np.zeros((len(seq), len(seq)), np.float32)
        raise RuntimeError("sequence %r absent from the compact ESM tables; refusing the bank" % seq)

    def _embed(seq):
        if seq in e32:
            return e32[seq]
        if seq == PROBE_SEQ:
            return np.zeros((len(seq), 32), np.float32)
        raise RuntimeError("sequence %r absent from s5/esm32.npz" % seq)

    def _con(seq):
        if seq in cn:
            return cn[seq]
        if seq == PROBE_SEQ:
            return np.zeros((len(seq), len(seq)), np.float32)
        raise RuntimeError("sequence %r absent from s7/repr_cache/esmcon.npz" % seq)

    D.esm_raw, D.esm_embed, D.esm_contacts = _raw, _embed, _con
    D._pipeline_guarded = True
    try:
        import esm_features as ef
        ef.raw = _raw
    except Exception:
        pass
    _T["guarded"] = True


# ============================================================ features (the shipped construction)
def extra_block(seq, emb, cn, i, j):
    """The ESM half of core.predict.features, split into its contact part (static, 13 columns)
    and its embedding part (4k columns), same column order, products in float64 then float32
    exactly as the shipped code does."""
    n = len(seq)
    cn = np.asarray(cn, np.float32); emb = np.asarray(emb, np.float64)
    deg = cn.sum(1, keepdims=True) / max(1, n - 1)
    sep = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    shells = np.stack([np.where((sep >= a) & (sep < b), cn, 0.0).sum(1) / max(1, n)
                       for a, b in ((3, 5), (5, 9), (9, 100))], axis=1)
    nb = np.stack([cn[np.clip(i + a, 0, n - 1), np.clip(j + b, 0, n - 1)]
                   for a, b in ((-1, -1), (1, 1), (-1, 1), (1, -1))], axis=1)
    static = np.column_stack([cn[i, j], nb, deg[i], deg[j], shells[i], shells[j]]).astype(np.float32)
    embb = np.column_stack([emb[i], emb[j], emb[i] * emb[j], np.abs(emb[i] - emb[j])]).astype(np.float32)
    return static, embb


def features_for(rung, seq, fold):
    """(X (npairs, d_in) float32, i, j) for ONE sequence, materialised."""
    base, i, j = PR.pair_features(seq)
    if rung == "noesm":
        return base, i, j
    emb, cn = residue_block(rung, seq, fold)
    static, embb = extra_block(seq, emb, cn, i, j)
    return np.hstack([base, static, embb]).astype(np.float32), i, j


def d_in_for(rung):
    base = PR.pair_features(PROBE_SEQ)[0].shape[1]                       # 42
    k = {"noesm": 0, "conly": 0, "pca32": 32, "wide": 32, "pca32f": 32, "esm8m": 32, "pca128": 128,
         "raw": 1280}[rung]
    return base if rung == "noesm" else base + 13 + 4 * k


def build_training(rung, fold, verbose=True):
    """Training arrays for a fold.  Materialised X for noesm/pca32/wide (<= 183 columns);
    (S static, R residue table, ii, jj) for the lean rungs.  Y = np.digitize(d, BIN_EDGES)."""
    ents = train_entries(fold)
    t0 = time.time()
    if rung in LEAN:
        Ss, Rs, II, JJ, Ys, off = [], [], [], [], [], 0
        for e in ents:
            base, i, j = PR.pair_features(e.seq)
            if not len(base):
                continue
            emb, cn = residue_block(rung, e.seq, fold)
            static, _ = extra_block(e.seq, emb, cn, i, j)
            Ss.append(np.hstack([base, static]).astype(np.float32))
            Rs.append(np.asarray(emb, np.float32))
            II.append((i + off).astype(np.int64)); JJ.append((j + off).astype(np.int64))
            Ys.append(np.digitize(np.linalg.norm(e.ca[i] - e.ca[j], axis=1), PR.BIN_EDGES))
            off += len(e.seq)
        out = dict(S=np.concatenate(Ss, 0), R=np.concatenate(Rs, 0), ii=np.concatenate(II),
                   jj=np.concatenate(JJ), Y=np.concatenate(Ys))
        npairs = len(out["Y"])
    else:
        Xs, Ys = [], []
        for e in ents:
            X, i, j = features_for(rung, e.seq, fold)
            if not len(X):
                continue
            Xs.append(X)
            Ys.append(np.digitize(np.linalg.norm(e.ca[i] - e.ca[j], axis=1), PR.BIN_EDGES))
        out = dict(X=np.vstack(Xs), Y=np.concatenate(Ys))
        npairs = len(out["Y"])
    if verbose:
        print("  rung %-7s fold %d: %d chains, %d pairs, d_in %d, features built in %.0fs"
              % (rung, fold, len(ents), npairs, d_in_for(rung), time.time() - t0), flush=True)
    return out


# ============================================================ the lean trainer
class LeanTrainer:
    """core.predict.MLP.fit semantics with each minibatch's feature rows assembled on the fly:
    X_p = [S[p], R[ii_p], R[jj_p], R[ii_p]*R[jj_p], |R[ii_p]-R[jj_p]|].  Same net, loss, optimiser,
    schedule, permutation stream and seed; the only difference is that the feature moments are
    accumulated in float64 over chunks instead of numpy's float32 X.mean(0)/X.std(0)."""

    def __init__(self, d_in, width=384, depth=3, seed=0):
        self.m = PR.MLP(d_in, width=width, depth=depth, seed=seed, dropout=0.0)

    @staticmethod
    def _rows(S, R, ii, jj, b):
        a, c = R[ii[b]].astype(np.float64), R[jj[b]].astype(np.float64)
        return np.hstack([S[b], a, c, a * c, np.abs(a - c)]).astype(np.float32)

    def _moments(self, S, R, ii, jj, chunk=None):
        n = len(ii); d = S.shape[1] + 4 * R.shape[1]
        #: chunk rows so the float64 chunk stays near 160 MB whatever the width (raw is 5175-d)
        chunk = chunk or max(256, int(2.0e7 / d))
        s1 = np.zeros(d); s2 = np.zeros(d)
        for a in range(0, n, chunk):
            b = np.arange(a, min(n, a + chunk))
            X = self._rows(S, R, ii, jj, b).astype(np.float64)
            s1 += X.sum(0); s2 += (X * X).sum(0)
        mu = s1 / n
        var = np.maximum(s2 / n - mu * mu, 0.0)
        return mu.astype(np.float32), (np.sqrt(var) + 1e-3).astype(np.float32)

    def fit(self, S, R, ii, jj, Y, epochs=EPOCHS, batch=BATCH, lr=LR, smooth=SMOOTH, verbose=False):
        t = _torch()
        m = self.m
        m.mu, m.sd = self._moments(S, R, ii, jj)
        idx = np.arange(PR.NBINS)[None, :]
        T = np.exp(-((idx - Y[:, None]) ** 2) / (2 * smooth ** 2)).astype(np.float32)
        T /= T.sum(1, keepdims=True)
        T = t.as_tensor(T)
        W = t.ones(len(Y)); W = W / W.mean()
        opt = t.optim.AdamW(m.net.parameters(), lr=lr, weight_decay=1e-4)
        sched = t.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr,
                                                total_steps=epochs * max(1, len(Y) // batch + 1))
        n = len(Y)
        g = t.Generator().manual_seed(0)
        m.net.train()
        for ep in range(epochs):
            perm = t.randperm(n, generator=g)
            tot = 0.0
            for s in range(0, n, batch):
                b = perm[s:s + batch]
                Xb = t.as_tensor((self._rows(S, R, ii, jj, b.numpy()) - m.mu) / m.sd)
                loss = (-(T[b] * t.log_softmax(m.net(Xb), 1)).sum(1) * W[b]).mean()
                opt.zero_grad(); loss.backward(); opt.step(); sched.step()
                tot += float(loss.detach()) * len(b)
            if verbose and ep % 10 == 0:
                print("    epoch %d loss %.4f" % (ep, tot / n), flush=True)
        m.net.eval()
        return m


def model_path(rung, fold, seed=0, tag=""):
    return os.path.join(MODELS, "%s_fold%d_s%d%s.pt" % (rung, fold, seed, ("_" + tag) if tag else ""))


def save_model(m, path, meta):
    t = _torch()
    t.save({"sd": m.net.state_dict(), "mu": m.mu, "s": m.sd, **meta}, path + ".tmp")
    os.replace(path + ".tmp", path)


def load_model(path):
    t = _torch()
    z = t.load(path, weights_only=False)
    cfg = z.get("cfg", {})
    m = PR.MLP(int(z["d_in"]), width=int(cfg.get("width", 384)), depth=int(cfg.get("depth", 3)),
               seed=int(z.get("seed", 0)), dropout=0.0)
    m.net.load_state_dict(z["sd"]); m.mu, m.sd = z["mu"], z["s"]; m.net.eval()
    return m, z


def train_rung(rung, folds=range(5), seed=0, epochs=EPOCHS, tag="", verbose=True):
    """Train (or load) the rung's fold models.  One .pt per fold is the checkpoint."""
    assert rung in TRAINABLE, rung
    install_guard()
    out = {}
    for fold in folds:
        path = model_path(rung, fold, seed, tag)
        if os.path.exists(path):
            out[fold] = load_model(path)[0]
            continue
        t0 = time.time()
        data = build_training(rung, fold, verbose=verbose)
        d_in = d_in_for(rung)
        arch = ARCH.get(rung, dict(width=384, depth=3))
        if rung in LEAN:
            tr = LeanTrainer(d_in, seed=seed, **arch)
            m = tr.fit(data["S"], data["R"], data["ii"], data["jj"], data["Y"], epochs=epochs, verbose=verbose)
        else:
            assert data["X"].shape[1] == d_in, (data["X"].shape, d_in)
            m = PR.MLP(d_in, seed=seed, dropout=0.0, **arch).fit(data["X"], data["Y"], epochs=epochs, verbose=verbose)
        secs = time.time() - t0
        save_model(m, path, {"rung": rung, "fold": fold, "seed": seed, "epochs": epochs, "d_in": d_in,
                             "cfg": arch, "n_pairs": int(len(data["Y"])), "train_secs": round(secs, 1),
                             "provenance": ST.provenance(__file__)})
        if verbose:
            print("  rung %-7s fold %d seed %d: trained in %.0fs -> %s" % (rung, fold, seed, secs, path), flush=True)
        out[fold] = m
        del data
    return out


# ============================================================ posteriors per target
def shipped_dg(pdb):
    """The pinned leave-fold-out posterior, cached by s12.instrument (probabilities, risk, grid)."""
    return I.distogram(pdb)


def posterior(rung, pdb, seq, fold, models=None, seed=0, lam=None, D_pool=None):
    """The rung's per-pair 17-bin posterior for one target.  DEPLOYABLE: no native enters."""
    if rung == "shipped":
        dg = shipped_dg(pdb)
        return np.asarray(dg["prob"], np.float64), np.asarray(dg["i"]), np.asarray(dg["j"])
    if rung == "mix":
        dg = shipped_dg(pdb)
        P0 = np.asarray(dg["prob"], np.float64)
        H = pool_histogram(D_pool)
        return (1.0 - lam) * P0 + lam * H, np.asarray(dg["i"]), np.asarray(dg["j"])
    if rung == "pairnet":
        install_guard()
        d = PR.pairnet_distogram(seq, fold=fold, model=PR.pairnet_train_fold(fold, verbose=False))
        return np.asarray(d.prob, np.float64), d.i, d.j
    m = models[fold] if models is not None else load_model(model_path(rung, fold, seed))[0]
    X, i, j = features_for(rung, seq, fold)
    return m.predict_proba(X).astype(np.float64), i, j


def pool_histogram(D_pool, eps=1e-3):
    """(npairs, 17): the K=500 pool's own per-pair distance histogram over the shipped bins."""
    b = np.digitize(np.asarray(D_pool, float), PR.BIN_EDGES)                     # (K, npairs)
    H = np.stack([(b == c).mean(0) for c in range(PR.NBINS)], 1) + eps
    return H / H.sum(1, keepdims=True)


def risk_table(seq, prob, i, j):
    """A genuine core.predict.Distogram from the posterior: the shipped functional, never a
    re-implementation (the s24/priorladder discipline)."""
    d = PR.Distogram(seq, np.asarray(prob, np.float64), np.asarray(i), np.asarray(j))
    return {"grid": d.grid, "risk": np.asarray(d._risk, np.float32), "w": d.w,
            "expected": d.expected, "sd": d.sd}


# ============================================================ endpoints
def score_target(prob, i, j, seq, W_pool):
    """DEPLOYABLE half: scores and the top-75 point cloud.  Reads no native quantity."""
    Dp = I.pair_dists(W_pool, i, j)
    rt = risk_table(seq, prob, i, j)
    sc = np.asarray(I.shipped_score(rt, Dp), float)
    o = np.argsort(sc, kind="stable")
    top = W_pool[o[:TOPM]]
    Pm = I.pairwise_rmsd(top)
    C = I.superpose_batch(top, top[I.medoid(Pm)]).mean(0)
    return sc, o, C, rt, Dp


def endpoint(prob, i, j, seq, fold, u, P0=None, do_project=True):
    """ORACLE-READ-AFTER-SELECTION half.  Native quantities enter only after `score_target`."""
    idx = I.pool_idx(u); W_pool = np.asarray(u["W"][idx], float)
    sc, o, C, rt, Dp = score_target(prob, i, j, seq, W_pool)
    nat = np.asarray(u["nat_ca"], float); rr = np.asarray(u["rr"][idx], float)
    sel, ntie = ST.argmin_tied(sc, rr)
    rec = {"sel": float(sel), "n_tied": int(ntie), "cloud": float(I.ca_rmsd(C, nat)),
           "top75_best": float(rr[o[:TOPM]].min()), "top75_mean": float(rr[o[:TOPM]].mean())}
    if do_project:
        pr = I.project(C, seq, int(fold))
        rec["arm"] = float(I.ca_rmsd(pr["ca"], nat)); rec["fit"] = float(I.ca_rmsd(pr["fit_ca"], nat))
    Dt = I.pair_dists(nat, i, j)
    E = rt["expected"]
    rec["mae"] = float(np.abs(E - Dt).mean())
    nb = np.abs(np.asarray(PR.CENTRES)[None, :] - Dt[:, None]).argmin(1)
    O = np.eye(PR.NBINS)[nb]
    if P0 is not None:
        E0 = (np.asarray(P0) * np.asarray(PR.CENTRES)[None]).sum(1)
        rec.update(progress(np.asarray(prob), np.asarray(P0), O, E, E0, Dt))
    return rec


def progress(Q, P0, O, EQ, E0, Dt):
    """gam_eff and cos in probability space (S24 MASS construction) and location space (s25 loc)."""
    dq = (Q - P0).ravel(); do = (O - P0).ravel()
    dl = EQ - E0; dt = Dt - E0
    def _g(a, b):
        nb2 = float((b * b).sum()); na = float(np.sqrt((a * a).sum())); nbb = math.sqrt(nb2)
        return (float((a * b).sum() / nb2) if nb2 > 0 else float("nan"),
                float((a * b).sum() / (na * nbb)) if na > 0 and nbb > 0 else float("nan"),
                float(na / nbb) if nbb > 0 else float("nan"))
    gp, cp, rp = _g(dq, do); gl, cl, rl = _g(dl, dt)
    return {"gam_prob": gp, "cos_prob": cp, "amp_prob": rp, "gam_loc": gl, "cos_loc": cl, "amp_loc": rl}


def _signed_off():
    try:
        return SIGNOFF in open(LEDGER, encoding="utf-8").read()
    except OSError:
        return False


def result_path(rung, seed=0, tag=""):
    return os.path.join(RES, "p_ladder_%s_s%d%s.json" % (rung, seed, ("_" + tag) if tag else ""))


def eval_rung(rung, seed=0, tag="", verbose=True):
    """All 126 targets through the three endpoints, paired row-by-row with the shipped posterior.
    Resumable: rows are written every 10 targets."""
    if not _signed_off():
        raise SystemExit("PHASE GATE: '%s' is not in %s; eval refuses to run." % (SIGNOFF, LEDGER))
    install_guard()
    tg = I.targets()
    out = result_path(rung, seed, tag)
    rows = json.load(open(out))["rows"] if os.path.exists(out) else []
    done = {r["pdb"] for r in rows}
    models = train_rung(rung, seed=seed, verbose=verbose) if rung in TRAINABLE else None
    need = ["sel", "cloud", "arm", "fit", "mae", "gam_prob", "cos_prob", "gam_loc", "cos_loc"]
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in done:
            continue
        u = I.load_univ(pdb)
        P0, i0, j0 = posterior("shipped", pdb, t["seq"], t["fold"])
        rec = {"pdb": pdb, "n": t["n"], "fold": t["fold"], "seq": t["seq"]}
        if rung == "mix":
            idx = I.pool_idx(u); Dp = I.pair_dists(np.asarray(u["W"][idx], float), i0, j0)
            for lam in MIX_LAMS:
                Q, i, j = posterior("mix", pdb, t["seq"], t["fold"], lam=lam, D_pool=Dp)
                if lam == 0.0:
                    r0 = risk_table(t["seq"], P0, i0, j0)["risk"]; r1 = risk_table(t["seq"], Q, i, j)["risk"]
                    assert np.array_equal(r0, r1), "mix lam=0 is not the incumbent bit-exactly"
                e = endpoint(Q, i, j, t["seq"], t["fold"], u, P0=P0)
                rec.update({"%s@%g" % (k, lam): v for k, v in e.items()})
            rec.update({k: rec["%s@0" % k] for k in need})          # the lam=0 column doubles as the incumbent
        else:
            Q, i, j = posterior(rung, pdb, t["seq"], t["fold"], models=models, seed=seed)
            assert np.array_equal(i, i0) and np.array_equal(j, j0), "pair index mismatch"
            rec.update(endpoint(Q, i, j, t["seq"], t["fold"], u, P0=P0))
        rows.append(rec)
        if verbose:
            print("  [%3d/%d] %s %s sel %.3f cloud %.3f arm %.3f" % (c + 1, len(tg), rung, pdb,
                  rec["sel"], rec["cloud"], rec["arm"]), flush=True)
        if len(rows) % 10 == 0:
            ST.save_atomic(out, {"rung": rung, "seed": seed, "rows": rows}, rows=rows, n_expected=len(tg),
                           module_file=__file__)
        del u
    ST.save_atomic(out, {"rung": rung, "seed": seed, "rows": rows, "lams": list(MIX_LAMS) if rung == "mix" else None},
                   complete_keys=need, rows=rows, n_expected=len(tg), module_file=__file__)
    return rows


# ============================================================ the reproduction gates (native-free)
def gate_pca32(fold, seed=0, tag="", verbose=True):
    """Does the retrained pca32 fold model reproduce the pinned distogram_models/fold<f>_esm_frag.pt?
    Compared on the fold's own held-out targets: max |risk - risk_pinned|, the fraction of pairs
    whose posterior differs by more than 1e-6, and agreement of the K=500 argmin (an index, not an
    RMSD).  NATIVE-FREE."""
    install_guard()
    m, meta = load_model(model_path("pca32", fold, seed, tag))
    pinned = PR.MLP(PR.feature_width(True), seed=0, dropout=0.0).load(
        os.path.join(ROOT, "distogram_models", "fold%d_esm_frag.pt" % fold))
    rows = []
    for t in I.targets():
        if int(t["fold"]) != int(fold):
            continue
        X, i, j = features_for("pca32", t["seq"], fold)
        Xs, i2, j2 = PR.features(t["seq"], True)
        u = I.load_univ(t["pdb"]); idx = I.pool_idx(u); Wp = np.asarray(u["W"][idx], float)
        Dp = I.pair_dists(Wp, i, j)
        pa, pb = m.predict_proba(X), pinned.predict_proba(Xs)
        ra = risk_table(t["seq"], pa, i, j); rb = risk_table(t["seq"], pb, i, j)
        sa = np.asarray(I.shipped_score(ra, Dp)); sb = np.asarray(I.shipped_score(rb, Dp))
        dg = shipped_dg(t["pdb"])
        rows.append({"pdb": t["pdb"], "feat_maxabs": float(np.abs(X - Xs).max()),
                     "feat_frac_diff": float((X != Xs).mean()),
                     "prob_maxabs": float(np.abs(pa - pb).max()),
                     "risk_maxabs": float(np.abs(ra["risk"] - rb["risk"]).max()),
                     "pinned_vs_cache_risk_maxabs": float(np.abs(rb["risk"] - np.asarray(dg["risk"], np.float32)).max()),
                     "argmin_same": bool(int(np.argmin(sa)) == int(np.argmin(sb))),
                     "top75_overlap": float(len(np.intersect1d(np.argsort(sa, kind="stable")[:TOPM],
                                                               np.argsort(sb, kind="stable")[:TOPM])) / TOPM),
                     "score_spearman": float(np.corrcoef(np.argsort(np.argsort(sa)), np.argsort(np.argsort(sb)))[0, 1])})
        del u
    out = {"fold": fold, "seed": seed, "tag": tag, "train_secs": meta.get("train_secs"), "rows": rows,
           "summary": {k: float(np.mean([r[k] for r in rows])) for k in
                       ("feat_maxabs", "feat_frac_diff", "prob_maxabs", "risk_maxabs", "pinned_vs_cache_risk_maxabs",
                        "argmin_same", "top75_overlap", "score_spearman")},
           "provenance": ST.provenance(__file__)}
    p = os.path.join(RES, "p_ladder_gate_pca32_fold%d_s%d%s.json" % (fold, seed, ("_" + tag) if tag else ""))
    ST.save_atomic(p, out, module_file=__file__)
    if verbose:
        print(json.dumps(out["summary"], indent=1)); print("  ->", p)
    return out


def featurise_esm8m(batch=64, verbose=True):
    """esm2_t6_8M over every training sequence and target, from the on-disk checkpoint (no download).
    Writes s26/models/p_ladder/esm8m_raw.npz.  Never touches esm_cache.npz / esm_small.npz."""
    p = os.path.join(MODELS, "esm8m_raw.npz")
    if os.path.exists(p):
        print("exists:", p); return p
    t = _torch()
    import esm
    seqs = sorted({e.seq for e in D.load()} | {f.seq for f in D.load_fragments()} | {x["seq"] for x in I.targets()})
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval(); bc = alphabet.get_batch_converter(); layer = model.num_layers
    keys, reps, cons, offs, coffs = [], [], [], [0], [0]
    with t.no_grad():
        for s in range(0, len(seqs), batch):
            chunk = seqs[s:s + batch]
            _, _, toks = bc([("p%d" % k, q) for k, q in enumerate(chunk)])
            out = model(toks, repr_layers=[layer], return_contacts=True)
            R = out["representations"][layer].numpy(); C = out["contacts"].numpy()
            for m, q in enumerate(chunk):
                n = len(q); keys.append(q)
                reps.append(R[m, 1:n + 1].astype(np.float32)); cons.append(C[m, :n, :n].astype(np.float32).ravel())
                offs.append(offs[-1] + n); coffs.append(coffs[-1] + n * n)
            if verbose and (s // batch) % 20 == 0:
                print("  esm8m %d/%d" % (s + len(chunk), len(seqs)), flush=True)
    np.savez(p + ".tmp.npz", keys=np.array(keys, dtype=object), flat=np.concatenate(reps, 0),
             offs=np.array(offs, np.int64), cflat=np.concatenate(cons, 0), coffs=np.array(coffs, np.int64))
    os.replace(p + ".tmp.npz", p)
    print("wrote", p, "(%d sequences)" % len(keys)); return p


# ============================================================ report
def _rows(rung, seed=0, tag=""):
    p = result_path(rung, seed, tag)
    if not os.path.exists(p):
        return None
    z = json.load(open(p))
    return z["rows"] if z.get("complete") else None


def report(seed=0):
    if not _signed_off():
        raise SystemExit("PHASE GATE: report reads native RMSDs; refused before sign-off.")
    base = _rows("shipped", seed)
    if base is None:
        raise SystemExit("run `eval --rung shipped` first")
    pdbs = [r["pdb"] for r in base]; folds = ST.pinned_folds(pdbs)
    b = {k: np.array([r[k] for r in base]) for k in ("sel", "cloud", "arm", "fit")}
    print("\nC2 PRIOR-INPUT LADDER, n=%d.  PAIRED against the shipped posterior through the SAME path." % len(base))
    print("  shipped through this module: sel %.4f  cloud %.4f  arm %.4f  fit %.4f   (pinned 3.4540 / 3.0483 / 3.2148 / 3.2041)"
          % tuple(b[k].mean() for k in ("sel", "cloud", "arm", "fit")))
    summary = {}
    for rung in RUNGS:
        if rung == "shipped":
            continue
        rows = _rows(rung, seed)
        if rows is None:
            continue
        m = {r["pdb"]: r for r in rows}
        if any(p not in m for p in pdbs):
            continue
        rr = [m[p] for p in pdbs]
        line = {}
        for k in ("arm", "cloud", "sel"):
            r = ST.compare(np.array([x[k] for x in rr]), b[k], folds, names=pdbs, label="%s vs shipped [%s]" % (rung, k))
            print(ST.fmt(r)); line[k] = r
        g = {q: float(np.nanmean([x[q] for x in rr])) for q in ("gam_prob", "cos_prob", "gam_loc", "cos_loc", "mae")}
        print("    gam_eff/cos  prob-space %+.4f / %+.3f   loc-space %+.4f / %+.3f   MAE %.4f (diag only)"
              % (g["gam_prob"], g["cos_prob"], g["gam_loc"], g["cos_loc"], g["mae"]))
        summary[rung] = {"stats": line, "progress": g}
    p = os.path.join(RES, "p_ladder_report_s%d.json" % seed)
    ST.save_atomic(p, {"summary": summary}, module_file=__file__)
    print("  ->", p)
    return summary


# ============================================================ CLI
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("train", "gate", "eval", "report", "featurise-esm8m", "probe"))
    ap.add_argument("--rung", default="pca32")
    ap.add_argument("--folds", default="0,1,2,3,4", help="fold-processing order, e.g. 4,3,2,1,0")
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--epochs", type=int, default=EPOCHS)
    ap.add_argument("--tag", default="")
    a = ap.parse_args(argv)
    folds = [int(x) for x in a.folds.split(",")]
    if a.cmd == "train":
        train_rung(a.rung, folds=folds, seed=a.seed, epochs=a.epochs, tag=a.tag)
    elif a.cmd == "probe":
        t0 = time.time()
        train_rung(a.rung, folds=folds[:1], seed=a.seed, epochs=a.epochs, tag=a.tag or "probe")
        print("probe wall %.0fs" % (time.time() - t0))
        if a.rung == "pca32":
            gate_pca32(folds[0], a.seed, a.tag or "probe")
    elif a.cmd == "gate":
        gate_pca32(a.fold, a.seed, a.tag)
    elif a.cmd == "eval":
        eval_rung(a.rung, seed=a.seed, tag=a.tag)
    elif a.cmd == "report":
        report(a.seed)
    elif a.cmd == "featurise-esm8m":
        featurise_esm8m()
    return 0


if __name__ == "__main__":
    sys.exit(main())
