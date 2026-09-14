"""s26/p_coh.py -- PREREG_coherence: retrain the shipped MLP with (a) a short-favouring per-pair
weight (S19 A5b, `sw_util`) and (b) a per-chain separation-profile coherence penalty (`coh_l<lam>`).

Imports s26/p_ladder.py UNCHANGED (its training jobs may be in flight) and uses its feature
construction, its guard and its eval path; the models are written under s26/models/p_ladder/ with
the rung names `sw_util` and `coh_l0.3` / `coh_l1` / `coh_l3`, so `p_ladder.py eval --rung <name>`
scores them through the identical endpoint code once `L.TRAINABLE` is extended at eval time
(done here by `register()`; nothing in p_ladder.py is edited).

    python s26/p_coh.py train --rung sw_util|coh_l0.3|coh_l1|coh_l3 --folds 0,1,2,3,4
    python s26/p_coh.py eval  --rung <rung>          (after training; phase gate is in p_ladder)
    python s26/p_coh.py modes --rung <rung>          (ORACLE diagnostic: held-out shell residual)
"""
from __future__ import annotations

import argparse
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
from s26 import p_ladder as L                        # noqa: E402
from core import predict as PR                       # noqa: E402

RUNGS = ("sw_util", "coh_l0.3", "coh_l1", "coh_l3")
LAM = {"coh_l0.3": 0.3, "coh_l1": 1.0, "coh_l3": 3.0}
CEN = np.asarray(PR.CENTRES, np.float32)


def register():
    """Make p_ladder treat these rungs as trainable 183-d (pca32-input) rungs at eval time."""
    L.TRAINABLE = tuple(L.TRAINABLE) + RUNGS
    L.RUNGS = tuple(L.RUNGS) + RUNGS
    _rb = L.residue_block; _din = L.d_in_for

    def residue_block(rung, seq, fold):
        return _rb("pca32", seq, fold) if rung in RUNGS else _rb(rung, seq, fold)

    def d_in_for(rung):
        return _din("pca32") if rung in RUNGS else _din(rung)

    L.residue_block = residue_block; L.d_in_for = d_in_for


def util_weight(sep):
    w = np.maximum(6.0 - np.asarray(sep, np.float32) / 3.0, 1.0)
    return (w / w.mean()).astype(np.float32)


def build(fold):
    """Training arrays with the chain id and separation kept, in p_ladder's corpus order."""
    ents = L.train_entries(fold)
    Xs, Ys, Cs, Ss, Ds = [], [], [], [], []
    for c, e in enumerate(ents):
        X, i, j = L.features_for("pca32", e.seq, fold)
        if not len(X):
            continue
        d = np.linalg.norm(e.ca[i] - e.ca[j], axis=1)
        Xs.append(X); Ys.append(np.digitize(d, PR.BIN_EDGES)); Cs.append(np.full(len(X), c, np.int64))
        Ss.append((j - i).astype(np.int64)); Ds.append(d.astype(np.float32))
    return (np.vstack(Xs), np.concatenate(Ys), np.concatenate(Cs), np.concatenate(Ss), np.concatenate(Ds))


def fit(rung, X, Y, chain, sep, dtrue, epochs=L.EPOCHS, batch=L.BATCH, lr=L.LR, smooth=L.SMOOTH, seed=0, verbose=True):
    """core.predict.MLP.fit with (a) a per-example weight or (b) a per-chain shell-coherence penalty.
    Chain-grouped minibatches: the permutation is over CHAINS, each batch takes whole chains until
    >= `batch` pairs, so the penalty sees every pair of its chains."""
    t = L._torch()
    m = PR.MLP(X.shape[1], seed=seed, dropout=0.0)
    m.mu = X.mean(0); m.sd = X.std(0) + 1e-3
    Xn = t.as_tensor(m._norm(X))
    idx = np.arange(PR.NBINS)[None, :]
    T = np.exp(-((idx - Y[:, None]) ** 2) / (2 * smooth ** 2)).astype(np.float32); T /= T.sum(1, keepdims=True)
    T = t.as_tensor(T)
    W = t.as_tensor(util_weight(sep)) if rung == "sw_util" else t.ones(len(X))
    lam = LAM.get(rung, 0.0)
    shell_id = np.full(len(sep), -1, np.int64)
    for k, (a, b) in enumerate(PR.SHELLS):
        shell_id[(sep >= a) & (sep <= b)] = k
    shell_t = t.as_tensor(shell_id); chain_t = t.as_tensor(chain); dt = t.as_tensor(dtrue); cen = t.as_tensor(CEN)
    n_chain = int(chain.max()) + 1
    starts = np.flatnonzero(np.r_[True, chain[1:] != chain[:-1]]); ends = np.r_[starts[1:], len(chain)]
    opt = t.optim.AdamW(m.net.parameters(), lr=lr, weight_decay=1e-4)
    steps_per_epoch = max(1, len(X) // batch + 1)
    sched = t.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, total_steps=epochs * steps_per_epoch)
    g = t.Generator().manual_seed(0)
    m.net.train()
    for ep in range(epochs):
        order = t.randperm(len(starts), generator=g).numpy()
        tot = 0.0; nsteps = 0; pos = 0
        while pos < len(order):
            sel = []; cnt = 0
            while pos < len(order) and cnt < batch:
                c = order[pos]; sel.append(np.arange(starts[c], ends[c])); cnt += ends[c] - starts[c]; pos += 1
            b = t.as_tensor(np.concatenate(sel))
            logit = m.net(Xn[b])
            loss = (-(T[b] * t.log_softmax(logit, 1)).sum(1) * W[b]).mean()
            if lam > 0:
                Ed = (t.softmax(logit, 1) * cen[None]).sum(1)
                dev = Ed - dt[b]
                key = chain_t[b] * 5 + shell_t[b]
                ssum = t.zeros(n_chain * 5).index_add_(0, key, dev)
                scnt = t.zeros(n_chain * 5).index_add_(0, key, t.ones_like(dev))
                used = scnt > 0
                pen = ((ssum[used] / scnt[used]) ** 2).sum() / max(1, len(set(chain[b.numpy()].tolist())))
                loss = loss + lam * pen
            opt.zero_grad(); loss.backward(); opt.step()
            if nsteps < epochs * steps_per_epoch - 1:
                sched.step()
            tot += float(loss.detach()) * len(b); nsteps += 1
        if verbose and ep % 10 == 0:
            print("    epoch %d loss %.4f" % (ep, tot / len(X)), flush=True)
    m.net.eval()
    return m


def train(rung, folds, seed=0, verbose=True):
    assert rung in RUNGS, rung
    L.install_guard(); register()
    for fold in folds:
        path = L.model_path(rung, fold, seed)
        if os.path.exists(path):
            continue
        t0 = time.time()
        X, Y, chain, sep, dtrue = build(fold)
        if verbose:
            print("  rung %-8s fold %d: %d pairs, %d chains, d_in %d" % (rung, fold, len(X), int(chain.max()) + 1, X.shape[1]), flush=True)
        m = fit(rung, X, Y, chain, sep, dtrue, seed=seed, verbose=verbose)
        L.save_model(m, path, {"rung": rung, "fold": fold, "seed": seed, "epochs": L.EPOCHS, "d_in": int(X.shape[1]),
                               "cfg": dict(width=384, depth=3), "n_pairs": int(len(X)), "train_secs": round(time.time() - t0, 1),
                               "lam": LAM.get(rung, 0.0), "provenance": ST.provenance(__file__)})
        if verbose:
            print("  rung %-8s fold %d: %.0fs -> %s" % (rung, fold, time.time() - t0, path), flush=True)
        del X, Y


def modes(rung, seed=0):
    """ORACLE DIAGNOSTIC: held-out per-target shell-profile residual of the posterior mean, rung vs shipped."""
    L.install_guard(); register()
    rows = []
    for t in I.targets():
        u = I.load_univ(t["pdb"]); nat = np.asarray(u["nat_ca"], float)
        P0, i, j = L.posterior("shipped", t["pdb"], t["seq"], t["fold"])
        Q, i2, j2 = L.posterior(rung, t["pdb"], t["seq"], t["fold"], seed=seed)
        Dt = I.pair_dists(nat, i, j); sep = j - i
        def prof(P):
            e = (P * np.asarray(PR.CENTRES)[None]).sum(1) - Dt
            return np.array([abs(e[(sep >= a) & (sep <= b)].mean()) if ((sep >= a) & (sep <= b)).any() else np.nan for a, b in PR.SHELLS])
        rows.append({"pdb": t["pdb"], "fold": int(t["fold"]), "shipped": float(np.nanmean(prof(P0))), "rung": float(np.nanmean(prof(Q)))})
        del u
    pdbs = [r["pdb"] for r in rows]
    c = ST.compare(np.array([r["rung"] for r in rows]), np.array([r["shipped"] for r in rows]), ST.pinned_folds(pdbs), names=pdbs,
                   label="%s vs shipped: held-out |shell-profile residual| of E[d] (ORACLE diagnostic, A)" % rung)
    print(ST.fmt(c))
    ST.save_atomic(os.path.join(L.RES, "p_coh_modes_%s_s%d.json" % (rung, seed)), {"rows": rows, "stats": c}, module_file=__file__)
    return c


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("train", "eval", "modes")); ap.add_argument("--rung", required=True)
    ap.add_argument("--folds", default="0,1,2,3,4"); ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "train":
        train(a.rung, [int(x) for x in a.folds.split(",")], a.seed)
    elif a.cmd == "eval":
        L.install_guard(); register(); L.eval_rung(a.rung, seed=a.seed)
    else:
        modes(a.rung, a.seed)
