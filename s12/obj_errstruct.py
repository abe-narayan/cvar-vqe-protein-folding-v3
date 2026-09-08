"""Step 2: characterise the SHIPPED distogram's error in error-structure coordinates.

Per target, with e = expected - true over the pairs (|i-j| >= 2):
  scalar        MAE, RMSE, r, calibration slope (OLS dtrue ~ exp), bias
  modes         projection of e onto the CONSTANT mode and the SCALING mode (e ~ a*dtrue),
                per-shell bias, terminal-region bias
  spectrum      eigen-spectrum of the symmetric error matrix E; effective rank
                (participation ratio of |eigenvalues|), top-1/top-3 Frobenius share
  sign          autocorrelation of sign(e) between pairs sharing a residue, and along |i-j|
  realisability classical-MDS 3D embedding of the PREDICTED matrix -> X_pred:
                (a) CA-RMSD(X_pred, native) = the objective's own geometric optimum
                (b) residual ||d(X_pred) - exp|| = the non-Euclidean part of the prediction
  null          the same statistics for iid noise and for a coordinate-displacement field
                matched to the same MAE (so "is the real error like a displacement field?")
"""
from __future__ import annotations
import os, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from s12 import instrument as I
from s12 import obj_common as OC


def to_matrix(d, vec, fill_adj=3.80):
    n = d["n"]; M = np.zeros((n, n))
    M[d["i"], d["j"]] = vec; M = M + M.T
    for k in range(n - 1):
        M[k, k + 1] = M[k + 1, k] = fill_adj
    np.fill_diagonal(M, 0.0)
    return M


def err_matrix(d, e):
    n = d["n"]; E = np.zeros((n, n))
    E[d["i"], d["j"]] = e; E = E + E.T
    return E


def mds3(D):
    """Classical MDS to 3D from a (possibly non-Euclidean) distance matrix."""
    n = len(D); J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(G)
    idx = np.argsort(w)[::-1][:3]
    w3 = np.maximum(w[idx], 0.0)
    return V[:, idx] * np.sqrt(w3)[None, :]


def spectrum_stats(E):
    w = np.linalg.eigvalsh(E)
    a = np.abs(w); s = a.sum()
    if s <= 0:
        return dict(eff_rank=float("nan"), top1=float("nan"), top3=float("nan"))
    p = a / s
    eff = float(np.exp(-(p * np.log(p + 1e-30)).sum()))       # participation (entropy) rank
    f = np.sort(w ** 2)[::-1]; fs = f.sum()
    return dict(eff_rank=eff, top1=float(f[0] / fs), top3=float(f[:3].sum() / fs),
                lam_max=float(w[np.argmax(np.abs(w))]))


def sign_share_stats(d, e):
    """Do errors sharing a residue agree in sign more than chance?"""
    i, j = d["i"], d["j"]; s = np.sign(e)
    npair = len(e); agree = []; agree_null = []
    for a in range(npair):
        m = ((i == i[a]) | (i == j[a]) | (j == i[a]) | (j == j[a]))
        m[a] = False
        if m.any():
            agree.append(float((s[m] == s[a]).mean()))
    return float(np.mean(agree)) if agree else float("nan")


def analyse(d, vec, tag):
    """vec = a predicted pair-distance vector."""
    dt = d["dtrue"]; e = np.asarray(vec, float) - dt
    n = d["n"]; sep = d["sep"]
    mae = float(np.abs(e).mean()); rmse = float(np.sqrt((e ** 2).mean()))
    r = float(np.corrcoef(vec, dt)[0, 1]) if np.std(vec) > 1e-9 else float("nan")
    # calibration: dtrue = a + b * pred
    b = np.polyfit(vec, dt, 1) if np.std(vec) > 1e-9 else [np.nan, np.nan]
    b2 = np.polyfit(dt, vec, 1) if np.std(dt) > 1e-9 else [np.nan, np.nan]
    # calibration AFTER removing everything predictable from |i-j| alone (shell means):
    # this is the slope that says how much of the pair-specific signal is real.
    sm_t = np.zeros_like(dt); sm_p = np.zeros_like(dt)
    for s in np.unique(sep):
        m = sep == s
        sm_t[m] = dt[m].mean(); sm_p[m] = np.asarray(vec)[m].mean()
    rt = dt - sm_t; rp = np.asarray(vec, float) - sm_p
    slope_sep = float((rt @ rp) / (rp @ rp)) if (rp @ rp) > 1e-9 else float("nan")
    r_sep = float(np.corrcoef(rt, rp)[0, 1]) if rp.std() > 1e-9 and rt.std() > 1e-9 else float("nan")
    # modes
    const = float(e.mean())
    ec = e - const
    dc = dt - dt.mean()
    a_scale = float((ec @ dc) / (dc @ dc))
    resid = ec - a_scale * dc
    shell = {}
    for s in range(2, min(sep.max(), 8) + 1):
        m = sep == s
        if m.any():
            shell[int(s)] = float(e[m].mean())
    far = sep > 7
    shell["8+"] = float(e[far].mean()) if far.any() else None
    E = err_matrix(d, e)
    sp = spectrum_stats(E)
    # realisability
    Dm = to_matrix(d, vec)
    X = mds3(Dm)
    dpred = np.linalg.norm(X[d["i"]] - X[d["j"]], axis=-1)
    return dict(tag=tag, mae=mae, rmse=rmse, r=r, slope=float(b[0]), intercept=float(b[1]),
                slope_pred_on_true=float(b2[0]), slope_sep=slope_sep, r_sep=r_sep,
                const_frac=float(const ** 2 * len(e) / (e @ e)),
                scale_coef=a_scale,
                scale_frac=float((a_scale ** 2) * (dc @ dc) / (e @ e)),
                resid_frac=float((resid @ resid) / (e @ e)),
                shell_bias=shell, sign_agree=sign_share_stats(d, e), **sp,
                mds_rmsd=float(I.ca_rmsd(X, d["nat"])),
                mds_selfmae=float(np.abs(dpred - vec).mean()),
                mds_mae_to_true=float(np.abs(dpred - dt).mean()))


def main():
    tg = OC.targets()
    rows = []
    rng = np.random.default_rng(7)
    for t in tg:
        d = OC.load(t["pdb"])
        rec = analyse(d, d["exp"], "distogram")
        rec.update(pdb=t["pdb"], n=t["n"], fold=t["fold"], pool_best=float(d["rr"].min()),
                   pool_mean=float(d["rr"].mean()), fail18=t["pdb"] in I.FAIL18)
        rows.append(rec)
        # matched-MAE nulls
        for kind, kw in (("iid", {}), ("sign", {}), ("coord", {"L": 3.0}), ("coord", {"L": 0.5}),
                         ("lowrank", {"rank": 1}), ("shuffle", {})):
            f = OC.make_corruptor(d, kind, rng, **kw)
            amp = OC.match_amp(f, d["dtrue"], rec["mae"])
            v = f(amp)
            rr = analyse(d, v, f"null_{kind}_{kw.get('L', kw.get('rank',''))}")
            rr.update(pdb=t["pdb"], n=t["n"], fold=t["fold"], fail18=t["pdb"] in I.FAIL18)
            rows.append(rr)
        # oracle-MDS control: what does classical MDS of the TRUE matrix give
    I.write("obj_errstruct", {"rows": rows})
    # console summary
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        by[r["tag"]].append(r)
    for k, v in by.items():
        f = lambda q: float(np.nanmean([x[q] for x in v]))
        print(f"{k:20s} mae {f('mae'):.3f} r {f('r'):.3f} slope_t~p {f('slope'):.3f} "
              f"slope_p~t {f('slope_pred_on_true'):.3f} slope_sep {f('slope_sep'):.3f} r_sep {f('r_sep'):.3f} "
              f"effrank {f('eff_rank'):.2f} top1 {f('top1'):.3f} scale_frac {f('scale_frac'):.3f} "
              f"sign_agree {f('sign_agree'):.3f} mds_rmsd {f('mds_rmsd'):.3f}")


def _match_amp(d, kind, target_mae, rng, **kw):
    lo, hi = 1e-4, 30.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        v = OC.corrupt(d, kind, mid, np.random.default_rng(1234), **kw)
        m = float(np.abs(v - d["dtrue"]).mean())
        if m < target_mae:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    main()
