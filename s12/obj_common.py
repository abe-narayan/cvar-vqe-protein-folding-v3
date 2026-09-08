"""Shared apparatus for the OBJECTIVE / ERROR-STRUCTURE agent (sprint 12).

Compact per-target cache (pool only) + the real downstream path as a pure function of a
TARGET PAIR-DISTANCE OBJECTIVE, so that arbitrary (corrupted / alternative) objectives can
be pushed through the production operator chain:

    score K=500 pool  ->  top-M=75  ->  coordinate average  ->  project  ->  CA-RMSD

Nothing here reads benchmark60.  `rr` / `nat_ca` are used for EVALUATION and, where a model
is trained, only as leave-fold-out labels.
"""
from __future__ import annotations
import os, json, math
import numpy as np

from s12 import instrument as I

CACHE = os.path.join(I.ROOT, "s12", "cache")
PACK = os.path.join(CACHE, "obj_pack.npz")


# --------------------------------------------------------------------------- cache
def build_pack(force=False):
    """Pool-only pack: for each target, W(500,n,3), rr(500), nat(n,3), pair index,
    native pair distances, distogram expected/sd/prob/risk/grid."""
    if os.path.exists(PACK) and not force:
        return
    tg = I.targets()
    out = {}
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        i, j = I.pair_index(t["n"])
        W = u["W"][p]
        dg = I.distogram(pdb, t["seq"], t["fold"])
        assert np.array_equal(dg["i"], i) and np.array_equal(dg["j"], j), pdb
        out[f"{pdb}/W"] = W.astype(np.float32)
        out[f"{pdb}/rr"] = u["rr"][p].astype(np.float32)
        out[f"{pdb}/nat"] = u["nat_ca"].astype(np.float32)
        out[f"{pdb}/org"] = u["org"][p]
        out[f"{pdb}/sim"] = u["sim"][p].astype(np.float32)
        out[f"{pdb}/S"] = u["S"][p].astype(np.int8)
        out[f"{pdb}/PHI"] = u["PHI"][p].astype(np.float32)
        out[f"{pdb}/PSI"] = u["PSI"][p].astype(np.float32)
        out[f"{pdb}/D"] = I.pair_dists(W, i, j).astype(np.float32)
        out[f"{pdb}/dtrue"] = I.pair_dists(u["nat_ca"], i, j).astype(np.float32)
        out[f"{pdb}/exp"] = dg["expected"].astype(np.float32)
        out[f"{pdb}/sd"] = dg["sd"].astype(np.float32)
        out[f"{pdb}/prob"] = dg["prob"].astype(np.float32)
        out[f"{pdb}/risk"] = dg["risk"].astype(np.float32)
    out["grid"] = np.asarray(I.distogram(tg[0]["pdb"])["grid"], np.float32)
    out["centres"] = np.asarray(I.distogram(tg[0]["pdb"])["centres"], np.float32)
    np.savez_compressed(PACK, **out)


_Z = None
_TG = None


def targets():
    global _TG
    if _TG is None:
        _TG = I.targets()
    return _TG


def _z():
    global _Z
    if _Z is None:
        build_pack()
        _Z = np.load(PACK)
    return _Z


def load(pdb):
    """dict for one target: W, rr, nat, D, dtrue, exp, sd, prob, risk, i, j, n, seq, fold."""
    z = _z()
    t = {x["pdb"]: x for x in targets()}[pdb]
    n = t["n"]
    i, j = I.pair_index(n)
    d = {k: np.asarray(z[f"{pdb}/{k}"]) for k in
         ("W", "rr", "nat", "D", "dtrue", "exp", "sd", "prob", "risk", "org", "sim", "S", "PHI", "PSI")}
    d["W"] = d["W"].astype(np.float64); d["nat"] = d["nat"].astype(np.float64)
    d["D"] = d["D"].astype(np.float64); d["dtrue"] = d["dtrue"].astype(np.float64)
    d["exp"] = d["exp"].astype(np.float64)
    d.update(pdb=pdb, n=n, seq=t["seq"], fold=t["fold"], i=i, j=j,
             sep=(j - i).astype(int), grid=np.asarray(z["grid"], np.float64),
             centres=np.asarray(z["centres"], np.float64))
    return d


# --------------------------------------------------------------------------- scorers
def score_l1(D, dtarget, w=None):
    """mean_p w_p |D_cand,p - dtarget_p| ; lower better.  D (b,npairs)."""
    r = np.abs(np.asarray(D) - np.asarray(dtarget)[None, :])
    if w is None:
        return r.mean(1)
    w = np.asarray(w, float); w = w / w.sum()
    return r @ w


def score_bayes(D, risk, grid):
    """The shipped Bayes-risk score."""
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


# --------------------------------------------------------------------------- the path
def emit(d, sc, M=75, lam=0.3, sub=None):
    """Push a candidate SCORE vector (lower better) through the production terminal path.
    Returns dict with fit_rmsd (lam=0 synthesis arm), lam_rmsd, argmin_rmsd, sub_best,
    sub (top-M indices into the pool)."""
    if sub is None:
        sub = np.argsort(sc, kind="stable")[:M]
    W = d["W"][sub]
    C, _ = I.coordinate_average(W)
    if lam is None:                      # cheap proxy arm: no projection
        return {"avg_rmsd": I.ca_rmsd(C, d["nat"]), "fit_rmsd": float("nan"),
                "lam_rmsd": float("nan"), "argmin_rmsd": float(d["rr"][int(np.argmin(sc))]),
                "sub_best": float(d["rr"][sub].min()), "sub_mean": float(d["rr"][sub].mean()),
                "sub": np.asarray(sub, int)}
    out = I.project(C, d["seq"], d["fold"], lam=lam)
    return {"avg_rmsd": I.ca_rmsd(C, d["nat"]),
            "fit_rmsd": I.ca_rmsd(out["fit_ca"], d["nat"]),
            "lam_rmsd": I.ca_rmsd(out["ca"], d["nat"]),
            "argmin_rmsd": float(d["rr"][int(np.argmin(sc))]),
            "sub_best": float(d["rr"][sub].min()),
            "sub_mean": float(d["rr"][sub].mean()),
            "sub": np.asarray(sub, int)}


def emit_from_target(d, dtarget, w=None, M=75, lam=0.3):
    return emit(d, score_l1(d["D"], dtarget, w), M=M, lam=lam)


# --------------------------------------------------------------------------- corruption
def make_corruptor(d, kind, rng, **kw):
    """Freeze the randomness, return f(amp) -> corrupted target vector.  Freezing matters:
    amplitude matching must use the SAME draw that is later pushed through the path."""
    dt = d["dtrue"]; n = d["n"]; i, j = d["i"], d["j"]; sep = d["sep"]
    if kind in ("iid", "hetero", "shell_short", "shell_long", "shrink"):
        z = rng.standard_normal(dt.size)
    elif kind == "sign":
        z = rng.choice([-1.0, 1.0], dt.size)
    elif kind == "coord":
        L = kw.get("L", 3.0); t = np.arange(n)
        Kc = np.exp(-0.5 * ((t[:, None] - t[None, :]) / max(L, 1e-6)) ** 2) + 1e-8 * np.eye(n)
        Lc = np.linalg.cholesky(Kc)
        disp0 = Lc @ rng.standard_normal((n, 3))
        return lambda a: np.maximum(np.linalg.norm(
            (d["nat"] + a * disp0)[i] - (d["nat"] + a * disp0)[j], axis=-1), 1.5)
    elif kind == "lowrank":
        r = int(kw.get("rank", 1)); E = np.zeros((n, n))
        for _ in range(r):
            v = rng.standard_normal(n); E += np.outer(v, v)
        z = (E / np.sqrt(r))[i, j]
    elif kind == "scale":
        return lambda a: np.maximum(dt * (1.0 + a), 1.5)
    elif kind == "scale_dn":
        return lambda a: np.maximum(dt * np.maximum(1.0 - a, 0.05), 1.5)
    elif kind == "shuffle":
        # keep the exact error MARGINAL of the real distogram, destroy its pair assignment
        e0 = d["exp"] - dt
        z = rng.permutation(e0)
        return lambda a: np.maximum(dt + a * z, 1.5)
    elif kind in ("shuffle_shell", "err_sep", "err_resid", "interp", "shrink_shell",
                  "pred_shellmean", "err_sep_shuf", "err_resid_shuf"):
        e0 = d["exp"] - dt
        sm = np.zeros_like(e0); pm = np.zeros_like(dt)
        for s in np.unique(sep):
            m = sep == s
            sm[m] = e0[m].mean(); pm[m] = dt[m].mean()
        if kind == "shuffle_shell":
            z = e0.copy()
            for s in np.unique(sep):
                m = sep == s
                z[m] = rng.permutation(e0[m])
        elif kind == "err_sep":
            z = sm
        elif kind == "err_resid":
            z = e0 - sm
        elif kind == "err_sep_shuf":                       # separation part, shells shuffled
            z = rng.permutation(np.unique(sep)); mp = dict(zip(np.unique(sep), z))
            z = np.array([sm[sep == mp[s]][0] for s in sep])
        elif kind == "err_resid_shuf":                     # pair-specific part, shuffled in shell
            z = e0 - sm
            zz = z.copy()
            for s in np.unique(sep):
                m = sep == s
                zz[m] = rng.permutation(z[m])
            z = zz
        elif kind == "interp":
            z = e0
        elif kind == "shrink_shell":                       # regress to the per-shell mean
            s0 = kw.get("slope", 0.376)
            v = np.maximum(pm + s0 * (dt - pm), 1.5)
            return lambda a: v
        elif kind == "pred_shellmean":
            zz = np.zeros_like(dt)
            for s in np.unique(sep):
                m = sep == s
                zz[m] = d["exp"][m].mean()
            return lambda a: np.maximum(zz, 1.5)
        return lambda a: np.maximum(dt + a * z, 1.5)
    else:
        raise ValueError(kind)
    if kind == "hetero":
        z = z * (dt / dt.mean())
    elif kind == "shell_short":
        z = z * (sep <= kw.get("cut", 4))
    elif kind == "shell_long":
        z = z * (sep > kw.get("cut", 4))
    elif kind == "shrink":
        s = kw.get("slope", 0.376); m = dt.mean()
        base = m + s * (dt - m)
        return lambda a: np.maximum(base + a * z, 1.5)
    return lambda a: np.maximum(dt + a * z, 1.5)


def match_amp(f, dtrue, target_mae, lo=1e-5, hi=60.0, it=45):
    """Bisect the amplitude of a frozen corruptor to a requested MAE."""
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if float(np.abs(f(mid) - dtrue).mean()) < target_mae:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def corrupt(d, kind, amp, rng, **kw):
    """Return a corrupted TARGET pair-distance vector derived from the NATIVE (ORACLE arm).

    kinds:
      iid            d + amp * N(0,1)                       (uncorrelated, symmetric)
      sign           d + amp * random_sign                  (uncorrelated, fixed magnitude)
      coord          distances of X + smooth displacement field (length scale kw['L'])
      scale          d * (1 + amp) (global expansion; amp may be negative)
      shrink         mean + s*(d-mean), s chosen to hit slope kw['slope'], + amp noise
      shell_short    iid noise only on |i-j| <= kw['cut']
      shell_long     iid noise only on |i-j| >  kw['cut']
      hetero         iid noise with sd proportional to d
      lowrank        error = sum_k a_k u_k u_k^T style symmetric low-rank pair error
    """
    dt = d["dtrue"].copy(); n = d["n"]; i, j = d["i"], d["j"]; sep = d["sep"]
    if kind == "iid":
        e = amp * rng.standard_normal(dt.size)
    elif kind == "sign":
        e = amp * rng.choice([-1.0, 1.0], dt.size)
    elif kind == "coord":
        L = kw.get("L", 3.0)
        # smooth Gaussian-process displacement field along the sequence
        t = np.arange(n)
        Kc = np.exp(-0.5 * ((t[:, None] - t[None, :]) / max(L, 1e-6)) ** 2) + 1e-8 * np.eye(n)
        Lc = np.linalg.cholesky(Kc)
        disp = (Lc @ rng.standard_normal((n, 3))) * amp
        X = d["nat"] + disp
        return np.maximum(np.linalg.norm(X[i] - X[j], axis=-1), 1.5)
    elif kind == "scale":
        return np.maximum(dt * (1.0 + amp), 1.5)
    elif kind == "shrink":
        s = kw.get("slope", 0.376); m = dt.mean()
        base = m + s * (dt - m)
        e = amp * rng.standard_normal(dt.size)
        return np.maximum(base + e, 1.5)
    elif kind == "shell_short":
        cut = kw.get("cut", 4); e = np.zeros(dt.size)
        msk = sep <= cut
        e[msk] = amp * rng.standard_normal(int(msk.sum()))
    elif kind == "shell_long":
        cut = kw.get("cut", 4); e = np.zeros(dt.size)
        msk = sep > cut
        e[msk] = amp * rng.standard_normal(int(msk.sum()))
    elif kind == "hetero":
        e = amp * rng.standard_normal(dt.size) * (dt / dt.mean())
    elif kind == "lowrank":
        r = int(kw.get("rank", 1))
        E = np.zeros((n, n))
        for _ in range(r):
            v = rng.standard_normal(n)
            E += np.outer(v, v)
        E = E / np.sqrt(r)
        e = amp * E[i, j]
    else:
        raise ValueError(kind)
    return np.maximum(dt + e, 1.5)


def mae_r(dtarget, dtrue):
    dtarget = np.asarray(dtarget, float); dtrue = np.asarray(dtrue, float)
    mae = float(np.abs(dtarget - dtrue).mean())
    if dtarget.std() < 1e-9 or dtrue.std() < 1e-9:
        r = float("nan")
    else:
        r = float(np.corrcoef(dtarget, dtrue)[0, 1])
    return mae, r


if __name__ == "__main__":
    build_pack()
    print("pack ok", os.path.getsize(PACK) / 2 ** 20, "MB")
