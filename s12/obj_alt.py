"""Step 3: ALTERNATIVE OBJECTIVE REPRESENTATIONS, all judged on EMITTED CA-RMSD.

Every scorer here consumes only the DEPLOYABLE view (the shipped leave-fold-out distogram
for the target sequence) plus, where marked LFO, statistics estimated on the OTHER folds'
targets.  None of them reads the native of the target being scored.

Scorers (score vector over the K=500 pool, lower = better):

  bayes           shipped L1 Bayes risk                                    (baseline)
  l1              mean |d_cand - E[d]|                                     (baseline')
  bayes_M25/M150  the shipped score with a different top-M                 (operator control)
  sdw             1/sd weighted L1 (the distogram's own uncertainty)
  shell_std       per-shell standardised L1: remove the |i-j| mean from both sides first
  scalefree       L1 after the candidate is optimally uniformly rescaled
  contact         cross entropy of the predicted contact map (P(d<8)) vs candidate contacts
  contact_l1      0.5*(z(l1) + z(contact))
  gram3           subspace alignment of the leading 3 eigenvectors of the Gram matrices
  shellw          L1 with per-shell weights, weights fitted LEAVE-FOLD-OUT to emitted quality
  maha            Mahalanobis L2 under an error covariance estimated LEAVE-FOLD-OUT in the
                  (sep_a, sep_b, share) basis -- the cheapest test of "structure beats mean"
  recal           per-shell affine recalibration of E[d], fitted LEAVE-FOLD-OUT, then L1
"""
from __future__ import annotations
import os, json, math
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from s12 import instrument as I
from s12 import obj_common as OC

MAXSEP = 15


# --------------------------------------------------------------------- helpers
def zsc(x):
    x = np.asarray(x, float); s = x.std()
    return (x - x.mean()) / (s if s > 1e-12 else 1.0)


def shell_mean(vec, sep):
    out = np.zeros_like(np.asarray(vec, float))
    for s in np.unique(sep):
        m = sep == s
        out[m] = np.asarray(vec, float)[m].mean()
    return out


def contact_prob(d, cut=8.0):
    """P(d_ij < cut) from the shipped 17-bin distogram (linear interp inside the bin)."""
    from core.predict import BIN_EDGES
    P = d["prob"]
    edges = np.asarray(BIN_EDGES, float)
    k = int(np.searchsorted(edges, cut))            # bins fully below cut: 0..k-1
    p = P[:, :k].sum(1)
    if k < P.shape[1]:
        lo = edges[k - 1] if k > 0 else 0.0
        hi = edges[k] if k < len(edges) else edges[-1] + 4.0
        p = p + P[:, k] * np.clip((cut - lo) / max(hi - lo, 1e-6), 0, 1)
    return np.clip(p, 1e-4, 1 - 1e-4)


def gram_of(D, n, i, j, fill=3.80):
    """Centred Gram matrix from a pair-distance vector (adjacent filled at `fill`)."""
    M = np.zeros((n, n)); M[i, j] = D; M = M + M.T
    for k in range(n - 1):
        M[k, k + 1] = M[k + 1, k] = fill
    J = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * J @ (M ** 2) @ J


# --------------------------------------------------------------------- scorers
def sc_bayes(d, **kw):
    return OC.score_bayes(d["D"], d["risk"], d["grid"])


def sc_l1(d, **kw):
    return OC.score_l1(d["D"], d["exp"])


def sc_sdw(d, **kw):
    return OC.score_l1(d["D"], d["exp"], w=1.0 / np.maximum(d["sd"], 0.5))


def sc_shell_std(d, **kw):
    sep = d["sep"]
    tp = d["exp"] - shell_mean(d["exp"], sep)
    Dc = d["D"] - shell_mean(d["D"].mean(0), sep)[None, :]
    return np.abs(Dc - tp[None, :]).mean(1)


# ------------------------------------------------- ARM 4: separation-residual SCORING
def _shell_masks(sep):
    return [(s, sep == s) for s in np.unique(sep)]


def _partial_sep(R, sep):
    """Remove each candidate's OWN per-shell mean deviation.  R (b, npairs)."""
    out = R.copy()
    for s, m in _shell_masks(sep):
        out[:, m] -= R[:, m].mean(1, keepdims=True)
    return out


def sc_sep_resid(d, target=None, **kw):
    """L1 on the deviation AFTER the candidate's own per-|i-j| mean deviation is removed:
    the score becomes blind to the sequence-separation marginal that every real fragment
    already satisfies.  Deployable (no labels)."""
    t = d["exp"] if target is None else target
    R = d["D"] - t[None, :]
    return np.abs(_partial_sep(R, d["sep"])).mean(1)


def sc_sep_resid_oracle(d, **kw):
    """DIAGNOSTIC: the same scoring rule fed the TRUE matrix (ceiling of the residual channel)."""
    return sc_sep_resid(d, target=d["dtrue"])


def sc_global_offset(d, **kw):
    """Weaker control: remove only ONE global mean deviation per candidate."""
    R = d["D"] - d["exp"][None, :]
    return np.abs(R - R.mean(1, keepdims=True)).mean(1)


def sc_sep_only(d, **kw):
    """The complement of sc_sep_resid: score ONLY the per-shell mean deviations (i.e. keep the
    separation marginal, throw the pair-specific residual away)."""
    R = d["D"] - d["exp"][None, :]
    keep = R - _partial_sep(R, d["sep"])
    return np.abs(keep).mean(1)


def sc_sep_mix(d, alpha=0.5, **kw):
    return (1 - alpha) * zsc(sc_l1(d)) + alpha * zsc(sc_sep_resid(d))


def sc_scalefree(d, **kw):
    D = d["D"]; t = d["exp"]
    a = (D @ t) / np.maximum((D * D).sum(1), 1e-9)          # per-candidate optimal scale
    return np.abs(D * a[:, None] - t[None, :]).mean(1)


def sc_contact(d, cut=8.0, **kw):
    p = contact_prob(d, cut)
    c = (d["D"] < cut).astype(float)
    return -(c * np.log(p)[None, :] + (1 - c) * np.log(1 - p)[None, :]).mean(1)


def sc_contact_l1(d, **kw):
    return 0.5 * (zsc(sc_l1(d)) + zsc(sc_contact(d)))


def sc_gram3(d, k=3, **kw):
    n = d["n"]; i, j = d["i"], d["j"]
    Gt = gram_of(d["exp"], n, i, j)
    w, V = np.linalg.eigh(Gt); Vt = V[:, np.argsort(w)[::-1][:k]]
    out = np.empty(len(d["D"]))
    for a, Dv in enumerate(d["D"]):
        G = gram_of(Dv, n, i, j)
        w2, V2 = np.linalg.eigh(G); V2 = V2[:, np.argsort(w2)[::-1][:k]]
        out[a] = -np.linalg.norm(Vt.T @ V2, "fro") ** 2     # subspace alignment, higher = better
    return out


def sc_shellw(d, weights=None, **kw):
    sep = np.clip(d["sep"], 0, MAXSEP)
    w = np.asarray(weights, float)[sep]
    return OC.score_l1(d["D"], d["exp"], w=np.maximum(w, 1e-6))


def sc_recal(d, recal=None, **kw):
    sep = np.clip(d["sep"], 0, MAXSEP)
    a = np.asarray(recal, float)[sep, 0]; b = np.asarray(recal, float)[sep, 1]
    return OC.score_l1(d["D"], a + b * d["exp"])


def sc_maha(d, Sig=None, ridge=0.15, **kw):
    """Mahalanobis L2 on the deviation vector under a (sep_a,sep_b,share) covariance."""
    sep = np.clip(d["sep"], 0, MAXSEP); i, j = d["i"], d["j"]
    npair = len(sep)
    share = ((i[:, None] == i[None, :]) | (i[:, None] == j[None, :]) |
             (j[:, None] == i[None, :]) | (j[:, None] == j[None, :])).astype(int)
    np.fill_diagonal(share, 2)
    C = Sig[sep[:, None], sep[None, :], share]
    C = 0.5 * (C + C.T)
    C = C + ridge * np.trace(C) / npair * np.eye(npair)
    R = d["D"] - d["exp"][None, :]
    try:
        Y = np.linalg.solve(C, R.T).T
    except np.linalg.LinAlgError:
        return sc_l1(d)
    return (R * Y).sum(1)


def sc_oracle(d, **kw):
    return OC.score_l1(d["D"], d["dtrue"])


def sc_sep_only_oracle(d, **kw):
    R = d["D"] - d["dtrue"][None, :]
    keep = R - _partial_sep(R, d["sep"])
    return np.abs(keep).mean(1)


def sc_random(d, seed=0, **kw):
    return np.random.default_rng(seed).standard_normal(len(d["D"]))


def sc_shellprofile(d, **kw):
    """Deployable: score ONLY the candidate's per-shell mean distance profile against the
    prediction's per-shell mean profile (one number per |i-j|)."""
    sep = d["sep"]; us = np.unique(sep)
    P = np.stack([d["D"][:, sep == s].mean(1) for s in us], 1)
    q = np.array([d["exp"][sep == s].mean() for s in us])
    return np.abs(P - q[None, :]).mean(1)


def sc_shellprofile_oracle(d, **kw):
    sep = d["sep"]; us = np.unique(sep)
    P = np.stack([d["D"][:, sep == s].mean(1) for s in us], 1)
    q = np.array([d["dtrue"][sep == s].mean() for s in us])
    return np.abs(P - q[None, :]).mean(1)


SCORERS = {"bayes": sc_bayes, "l1": sc_l1, "sdw": sc_sdw, "shell_std": sc_shell_std,
           "scalefree": sc_scalefree, "contact": sc_contact, "contact_l1": sc_contact_l1,
           "gram3": sc_gram3, "shellw": sc_shellw, "recal": sc_recal, "maha": sc_maha,
           "sep_resid": sc_sep_resid, "sep_resid_oracle": sc_sep_resid_oracle,
           "global_offset": sc_global_offset, "sep_only": sc_sep_only, "sep_mix": sc_sep_mix,
           "oracle": sc_oracle, "sep_only_oracle": sc_sep_only_oracle, "random": sc_random,
           "shellprofile": sc_shellprofile, "shellprofile_oracle": sc_shellprofile_oracle}


# --------------------------------------------------------------------- LFO fitting
def fit_lfo():
    """Per-fold: per-shell affine recalibration, per-shell |error| scale (-> shellw), and the
    (sep_a,sep_b,share) error covariance.  Trained on the OTHER folds' 126-instrument targets."""
    tg = OC.targets()
    per = {}
    for t in tg:
        d = OC.load(t["pdb"])
        e = d["exp"] - d["dtrue"]
        sep = np.clip(d["sep"], 0, MAXSEP)
        i, j = d["i"], d["j"]
        share = ((i[:, None] == i[None, :]) | (i[:, None] == j[None, :]) |
                 (j[:, None] == i[None, :]) | (j[:, None] == j[None, :])).astype(int)
        np.fill_diagonal(share, 2)
        per[t["pdb"]] = dict(fold=t["fold"], e=e, sep=sep, share=share,
                             exp=d["exp"], dtrue=d["dtrue"])
    folds = sorted({v["fold"] for v in per.values()})
    out = {}
    for f in folds:
        tr = [v for v in per.values() if v["fold"] != f]
        recal = np.zeros((MAXSEP + 1, 2)); recal[:, 1] = 1.0
        wsh = np.ones(MAXSEP + 1)
        num = np.zeros((MAXSEP + 1, MAXSEP + 1, 3)); den = np.zeros_like(num)
        for s in range(MAXSEP + 1):
            X = np.concatenate([v["exp"][v["sep"] == s] for v in tr]) if any((v["sep"] == s).any() for v in tr) else np.array([])
            Y = np.concatenate([v["dtrue"][v["sep"] == s] for v in tr]) if X.size else np.array([])
            if X.size > 30 and X.std() > 1e-6:
                b = np.polyfit(X, Y, 1); recal[s] = [b[1], b[0]]
                wsh[s] = 1.0 / max(np.abs(Y - (b[1] + b[0] * X)).mean(), 0.05)
            elif X.size:
                wsh[s] = 1.0 / max(np.abs(Y - X).mean(), 0.05)
        for v in tr:
            e = v["e"]; sep = v["sep"]; sh = v["share"]
            P = np.outer(e, e)
            np.add.at(num, (sep[:, None].repeat(len(sep), 1), sep[None, :].repeat(len(sep), 0), sh), P)
            np.add.at(den, (sep[:, None].repeat(len(sep), 1), sep[None, :].repeat(len(sep), 0), sh), 1.0)
        Sig = np.where(den > 5, num / np.maximum(den, 1), 0.0)
        # diagonal fallback: never let a shell get a zero variance
        for s in range(MAXSEP + 1):
            if Sig[s, s, 2] <= 1e-6:
                Sig[s, s, 2] = 1.0
        out[int(f)] = dict(recal=recal, wsh=wsh, Sig=Sig)
    return out


if __name__ == "__main__":
    L = fit_lfo()
    np.savez_compressed(os.path.join(OC.CACHE, "obj_lfo.npz"),
                        **{f"{f}/{k}": v[k] for f, v in L.items() for k in ("recal", "wsh", "Sig")})
    for f, v in L.items():
        print("fold", f, "recal slopes", np.round(v["recal"][2:9, 1], 3), "wsh", np.round(v["wsh"][2:9], 3))
