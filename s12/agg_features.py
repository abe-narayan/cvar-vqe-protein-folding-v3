"""Feature tensors for the LEARNED SET DECODER (experiment 2).

For each target and a candidate set of size m (default: the shipped top-75), builds

    X  (m, npairs, F)  the FULL per-pair candidate-vs-objective deviation map
    G  (m, Gf)         candidate-level scalars
    P  (m, m)          candidate-vs-candidate CA-RMSD
    Asup (m, n, 3)     candidates superposed on the set medoid  (native-free)
    nat  (n, 3)        ORACLE label
    rr   (m,)          ORACLE per-candidate CA-RMSD

Everything in X and G is computable at inference time from the deployable view
(window coordinates + the shipped leave-fold-out distogram + BLOSUM sim + org flag).
`rr`/`nat` are labels only.
"""
from __future__ import annotations
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_common as A
from s12 import agg_meta

FEAT = os.path.join(ROOT, "s12", "cache", "agg_feat2")   # v2: signed shell/global deviation + rg features
os.makedirs(FEAT, exist_ok=True)

PAIR_NAMES = ["d", "e", "sd", "d_minus_e", "z_de", "risk_d", "mu", "sig", "d_minus_mu",
              "z_dmu", "e_minus_mu", "med_minus_d", "sep", "entropy", "risk_mu", "rank_p",
              "shell_dev", "glob_dev"]
CAND_NAMES = ["z_score", "z_cons", "rank_frac", "org", "z_sim", "rmsd_medoid", "rmsd_avg",
              "z_rg", "mae_de", "corr_de", "z_risk_mean", "frac_pairs_best",
              "rg_abs", "rg_implied", "rg_minus_implied", "n_res", "meand_minus_meane",
              "rg_ratio"]
# columns of X carrying a SIGNED deviation (the coordinate the coherent-direction
# hypothesis lives in).  `abs_cols` ablation takes |.| of exactly these.
SIGNED_PAIR_COLS = [3, 4, 8, 9, 10, 11, 16, 17]
F = len(PAIR_NAMES)
GF = len(CAND_NAMES)


def _z(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / (s if s > 1e-12 else 1.0)


def build_one(pdb, seq, n, m=75, source="top75"):
    d = A.load(pdb, mem=False)
    mt = agg_meta.load(pdb)
    dg = I.distogram(pdb)
    Wp = d["Wp"]; sc = d["sc"].astype(float); rr = d["rr"].astype(float); sub = d["sub"].astype(int)
    if source == "top75":
        idx = sub
    elif source.startswith("topk"):
        k = int(source[4:])
        idx = np.argsort(sc)[:k]
    else:
        raise ValueError(source)
    W = Wp[idx]
    m = len(idx)
    pi, pj = d["pi"].astype(int), d["pj"].astype(int)
    npairs = len(pi)
    D = I.pair_dists(W, pi, pj)                      # (m, npairs)
    e = d["exp"].astype(float); sdv = np.maximum(d["sd"].astype(float), 1e-3)
    grid = dg["grid"]; risk = dg["risk"].astype(float)
    gi = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    rk = risk[np.arange(npairs)[None, :], gi]        # (m, npairs)
    mu = D.mean(0); sig = np.maximum(D.std(0), 1e-3); med = np.median(D, 0)
    gmu = np.clip(((mu - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    rmu = risk[np.arange(npairs), gmu]
    prob = dg["prob"].astype(float)
    ent = -(prob * np.log(prob + 1e-9)).sum(1)
    sep = (pj - pi) / float(n)
    rank_p = np.argsort(np.argsort(D, 0), 0) / max(m - 1, 1)

    ones = np.ones((m, npairs))
    dev = D - e                                            # SIGNED per-pair deviation
    shell_dev = np.zeros_like(dev)
    seps = (pj - pi)
    for s_ in np.unique(seps):
        msk = seps == s_
        shell_dev[:, msk] = dev[:, msk].mean(1, keepdims=True)
    glob_dev = dev.mean(1, keepdims=True) * ones
    X = np.stack([D, e * ones, sdv * ones, dev, dev / sdv, rk, mu * ones, sig * ones,
                  D - mu, (D - mu) / sig, (e - mu) * ones, med - D, sep * ones, ent * ones,
                  rmu * ones, rank_p, shell_dev, glob_dev], -1).astype(np.float32)

    P = I.pairwise_rmsd(W)
    medi = int(np.argmin(P.mean(1)))
    Asup = I.superpose_batch(W, W[medi])
    Cavg = Asup.mean(0)
    cons = P.mean(1)
    scs = sc[idx]
    rg = np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
    mae = np.abs(D - e).mean(1)
    Dc = D - D.mean(1, keepdims=True); ec = e - e.mean()
    corr = (Dc @ ec) / np.maximum(np.linalg.norm(Dc, axis=1) * np.linalg.norm(ec), 1e-9)
    best_at_pair = (rk == rk.min(0)).mean(1)
    # rg implied by a distance matrix:  rg^2 = (1/n^2) sum_{i<j} d_ij^2  (all pairs, incl |i-j|<2)
    def _rg_from(dv):
        full = np.zeros((n, n))
        full[pi, pj] = dv; full[pj, pi] = dv
        for k_ in range(1, 2):
            for a_ in range(n - k_):
                full[a_, a_ + k_] = full[a_ + k_, a_] = 3.8 * k_
        return np.sqrt((full ** 2).sum() / (2.0 * n ** 2))
    rg_impl = _rg_from(e)
    G = np.stack([_z(scs), _z(cons), np.argsort(np.argsort(scs)) / max(m - 1, 1),
                  mt["org"][idx].astype(float), _z(mt["sim"][idx].astype(float)),
                  I.kabsch_rmsd_batch(W, W[medi]), I.kabsch_rmsd_batch(W, Cavg),
                  _z(rg), mae, corr, _z(rk.mean(1)), best_at_pair,
                  rg, np.full(m, rg_impl), rg - rg_impl, np.full(m, float(n)),
                  D.mean(1) - e.mean(), rg / max(rg_impl, 1e-6)], -1).astype(np.float32)

    return {"X": X, "G": G, "P": P.astype(np.float32), "Asup": Asup.astype(np.float32),
            "nat": d["nat"].astype(np.float32), "rr": rr[idx].astype(np.float32),
            "idx": idx.astype(np.int32), "sc": scs.astype(np.float32),
            "pi": pi.astype(np.int32), "pj": pj.astype(np.int32)}


def build(source="top75"):
    for t in I.targets():
        path = os.path.join(FEAT, f"{t['pdb']}_{source}.npz")
        if os.path.exists(path):
            continue
        o = build_one(t["pdb"], t["seq"], t["n"], source=source)
        np.savez_compressed(path, **o)
        print(t["pdb"], flush=True)


def load(pdb, source="top75"):
    z = np.load(os.path.join(FEAT, f"{pdb}_{source}.npz"))
    return {k: z[k] for k in z.files}


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "top75"
    build(src)
    print("feat ok")
