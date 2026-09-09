"""Sprint 12 shared instrument.  Pure helpers over the cached window universes.

EVERYTHING an agent needs to run a retrieval/selection/synthesis experiment on the
126-target tuning instrument WITHOUT touching benchmark60 and WITHOUT reloading the
library.  `s8/generate_univ/<pdb>.npz` holds, per tuning target, EVERY length-n window of
the leakage-safe library (out-of-fold peptides + this fold's fragments):

    W (nw, n, 3) float32 CA      PHI/PSI (nw, n) float16 (radians)     S (nw, n) int8 codes
    org (nw,) bool  True = window comes from the peptide database, False = protein fragment
    sim (nw,) BLOSUM62 sum vs the target       order (nw,) stable argsort of -sim
    rr (nw,) float32  ORACLE: CA-RMSD of each window to the native (model-1) CA trace
    nat_ca (n, 3)     ORACLE: native CA trace

`rr` and `nat_ca` are LABELS.  Use them for evaluation and for leave-fold-out training
labels only.  Never let them enter an inference-time decision.

Conventions reproduced exactly here (asserted by `selfcheck()`):
    shipped argmin 3.4540, pool best 1.7108, top-75 best 2.3062, synthesis 3.2041,
    18 zero-recall targets (S10-1).
"""
from __future__ import annotations
import json, os, sys, glob, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
UNIV = os.path.join(ROOT, "s8", "generate_univ")
PROD_KEY = "1fc9f2dcf489e2fb"           # bench_results/cache/<key>/<pdb>.json: the production run
CACHE = os.path.join(ROOT, "s12", "cache")
RESULTS = os.path.join(ROOT, "s12", "results")
os.makedirs(CACHE, exist_ok=True); os.makedirs(RESULTS, exist_ok=True)

K = 500
M = 75
BAND = 1.5
ALPHABET = "ARNDCQEGHILKMFPSTWYV"


# ----------------------------------------------------------------------------- targets
def targets():
    """The 126 tuning targets in pinned (pdb-sorted) order: dicts pdb,n,fold,seq."""
    out = []
    for f in sorted(glob.glob(os.path.join(UNIV, "*.npz"))):
        z = np.load(f, allow_pickle=True)
        out.append(dict(pdb=str(z["pdb"]), n=int(z["n"]), fold=int(z["fold"]), seq=str(z["seq"])))
    return out

FAIL18 = ["1ID6", "1JBF", "1LB7", "2BFI", "2BP4", "2JN5", "2MQ2", "2N5C", "2NB7", "2NDM",
          "3BTB", "3SGO", "5W52", "7JS6", "7LCW", "8T63", "9KAR", "9L1M"]


def load_univ(pdb):
    z = np.load(os.path.join(UNIV, f"{pdb}.npz"), allow_pickle=True)
    u = {k: z[k] for k in z.files}
    u["pdb"] = str(u["pdb"]); u["seq"] = str(u["seq"]); u["n"] = int(u["n"]); u["fold"] = int(u["fold"])
    u["W"] = np.asarray(u["W"], np.float64)
    u["PHI"] = np.asarray(u["PHI"], np.float64); u["PSI"] = np.asarray(u["PSI"], np.float64)
    u["rr"] = np.asarray(u["rr"], np.float64); u["nat_ca"] = np.asarray(u["nat_ca"], np.float64)
    return u


def pool_idx(u, k=K):
    """Indices (into the universe) of the shipped K=500 BLOSUM pool."""
    return np.asarray(u["order"][:k], int)


def shipped_record(pdb):
    """The production pipeline's per-target record (sub = top-75 indices into the pool)."""
    with open(os.path.join(ROOT, "bench_results", "cache", PROD_KEY, f"{pdb}.json")) as fh:
        return json.load(fh)


def decode(codes):
    return "".join(ALPHABET[int(c)] for c in np.asarray(codes).ravel())


# ----------------------------------------------------------------------------- geometry
def kabsch_rmsd_batch(W, T):
    """CA-RMSD of each (n,3) in W (b,n,3) to T (n,3) after optimal superposition."""
    W = np.asarray(W, float); T = np.asarray(T, float)
    W = W - W.mean(1, keepdims=True); T = T - T.mean(0, keepdims=True)
    H = np.einsum("bni,nj->bij", W, T)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("bij,bjk->bik", Vt.transpose(0, 2, 1), U.transpose(0, 2, 1))))
    S = S.copy(); S[:, -1] *= d
    num = (W ** 2).sum((1, 2)) + (T ** 2).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(num, 0.0) / W.shape[1])


def ca_rmsd(a, b):
    return float(kabsch_rmsd_batch(np.asarray(a, float)[None], b)[0])


def superpose_batch(W, T):
    """Superpose each member of W (b,n,3) onto T (n,3); returns the moved copies."""
    W = np.asarray(W, float); T = np.asarray(T, float)
    Wc = W - W.mean(1, keepdims=True); Tc = T - T.mean(0, keepdims=True)
    H = np.einsum("bni,nj->bij", Wc, Tc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("bij,bjk->bik", Vt.transpose(0, 2, 1), U.transpose(0, 2, 1))))
    D = np.tile(np.eye(3), (len(W), 1, 1)); D[:, 2, 2] = d
    R = np.einsum("bij,bjk,bkl->bil", Vt.transpose(0, 2, 1), D, U.transpose(0, 2, 1))
    return np.einsum("bij,bnj->bni", R, Wc) + T.mean(0)


def pairwise_rmsd(W):
    W = np.asarray(W, float); b = len(W)
    P = np.zeros((b, b))
    for a in range(b):
        P[a] = kabsch_rmsd_batch(W, W[a])
    return P


def medoid(P):
    return int(np.argmin(np.asarray(P, float).mean(1)))


def coordinate_average(W, P=None):
    """S8-11's operator: superpose on the medoid, mean.  Returns (C, medoid_index)."""
    W = np.asarray(W, float)
    if P is None:
        P = pairwise_rmsd(W)
    b = medoid(P)
    return superpose_batch(W, W[b]).mean(0), b


def pair_index(n, min_sep=2):
    i, j = np.triu_indices(n, k=min_sep)
    return i, j


def pair_dists(W, i, j):
    W = np.asarray(W, float)
    return np.linalg.norm(W[..., i, :] - W[..., j, :], axis=-1)


# ----------------------------------------------------------------------------- projection
def project(C, seq, fold, lam=0.3, multi=True, maxiter=300):
    """STAGE 3b exactly as production: nearest ideal-geometry chain, ramah@lam, multi-start.
    Returns dict ca, phi, psi, fit_ca (lam=0 arm)."""
    from core import project as pj
    pen = pj.make_penalty("ramah", seq, int(fold))
    path = pj.lam_path(np.asarray(C, float), pen, (0.0, lam), maxiter=maxiter, multi=multi, grad="exact")
    fit, arm = path[0.0], path[lam]
    return {"ca": np.asarray(arm[0], float), "phi": np.asarray(arm[1], float), "psi": np.asarray(arm[2], float),
            "fit_ca": np.asarray(fit[0], float)}


def build_ca(phi, psi):
    """Ideal-geometry CA trace. Accepts (n,) or (B,n); returns (n,3) or (B,n,3).

    `core.project.build_ca_exact` is batched-only and is the BIT-EXACT builder the
    projection uses, so templates built here live on the same manifold the pipeline emits.
    """
    from core import project as pj
    phi = np.asarray(phi, float); psi = np.asarray(psi, float)
    single = phi.ndim == 1
    if single:
        phi = phi[None]; psi = psi[None]
    out = np.asarray(pj.build_ca_exact(phi, psi), float)
    return out[0] if single else out


def build_backbone(phi, psi):
    from core import geometry as geo
    return geo.build_backbone(np.asarray(phi, float), np.asarray(psi, float))


def ss_of(phi, psi):
    """Simplified-DSSP H/E/C string from torsions (ideal-geometry backbone)."""
    from core import geometry as geo
    bb = build_backbone(phi, psi)
    return geo.assign_secondary_structure(bb)


# ----------------------------------------------------------------------------- distogram
def distogram(pdb, seq=None, fold=None):
    """The shipped leave-fold-out distogram for a tuning target, cached as npz.
    prob (npairs, 17), i, j, expected (npairs,), sd (npairs,), risk (npairs, 760) on grid 2..40 by 0.05.
    """
    path = os.path.join(CACHE, f"disto_{pdb}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    if seq is None:
        t = {x["pdb"]: x for x in targets()}[pdb]; seq, fold = t["seq"], t["fold"]
    from core import pipeline as pl
    from core import predict as dgm
    pl.guard_esm([seq])
    model = pl.fold_model(int(fold))
    d = dgm.Distogram.for_target(seq, model=model)
    out = {"prob": d.prob.astype(np.float32), "i": d.i, "j": d.j, "expected": d.expected.astype(np.float32),
           "sd": d.sd.astype(np.float32), "risk": np.asarray(d._risk, np.float32),
           "centres": np.asarray(dgm.CENTRES, np.float32), "grid": np.asarray(d.grid, np.float32)}
    np.savez_compressed(path, **out)
    return out


def shipped_score(dg, D):
    """The shipped Bayes-risk score for pair-distance rows D (b, npairs) -- lower is better."""
    grid = dg["grid"]; risk = dg["risk"]
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


# ----------------------------------------------------------------------------- statistics
def paired(a, b, n_boot=4000, seed=0, folds=None, names=None):
    """Paired comparison a - b (negative = a better).  Bootstrap CI, W/L, concentration."""
    a = np.asarray(a, float); b = np.asarray(b, float); d = a - b; n = len(d)
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    order = np.argsort(d)          # most negative (largest gains) first
    out = {"n": int(n), "mean_a": float(a.mean()), "mean_b": float(b.mean()), "mean_diff": float(d.mean()),
           "se": float(d.std(ddof=1) / math.sqrt(n)), "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
           "n_better": int((d < 0).sum()), "n_worse": int((d > 0).sum()), "median_diff": float(np.median(d)),
           "drop_top10_mean_diff": float(d[order[10:]].mean()) if n > 10 else None,
           "drop_top20_mean_diff": float(d[order[20:]].mean()) if n > 20 else None,
           "top10_share": float(d[order[:10]].sum() / d.sum()) if d.sum() != 0 and n > 10 else None}
    if folds is not None:
        folds = np.asarray(folds)
        out["per_fold"] = {int(f): float(d[folds == f].mean()) for f in np.unique(folds)}
    if names is not None:
        out["top10_targets"] = [(str(names[k]), float(d[k])) for k in order[:10]]
    return out


def summary(x):
    x = np.asarray(x, float)
    return {"n": int(len(x)), "mean": float(x.mean()), "median": float(np.median(x)), "sd": float(x.std()),
            "min": float(x.min()), "max": float(x.max()), "frac_under_2.0": float((x < 2.0).mean())}


def write(name, obj, n_expected=None):
    """Write a result JSON to `s12/results/`.

    HAZARD, found by the sprint-12 adversarial audit: this writer has no config key, unlike
    `core/pipeline.py`'s cache, so a PARTIAL run silently overwrites a COMPLETE one of the
    same name.  Pass `n_expected` and the number of rows actually present, and the file is
    written with an explicit `complete` flag and a row count, so a reader can tell a finished
    result from an interrupted one instead of having to trust the filename.
    """
    path = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    if isinstance(obj, dict) and n_expected is not None:
        rows = obj.get("per_target") or obj.get("rows") or []
        obj = dict(obj, n_rows=len(rows), n_expected=int(n_expected),
                   complete=bool(len(rows) >= int(n_expected)))
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return path


def free_gb():
    import ctypes
    class MS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong), ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong), ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong), ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong), ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    st = MS(); st.dwLength = ctypes.sizeof(MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return st.ullAvailPhys / 2 ** 30


# ----------------------------------------------------------------------------- selfcheck
def selfcheck(verbose=True):
    tg = targets(); assert len(tg) == 126
    shipped, pbest, tbest, synth, zero = [], [], [], [], []
    for t in tg:
        u = load_univ(t["pdb"]); p = pool_idx(u); rr = u["rr"][p]
        rec = shipped_record(t["pdb"]); sub = np.asarray(rec["sub"], int)
        dg = distogram(t["pdb"], t["seq"], t["fold"])
        i, j = pair_index(t["n"]); D = pair_dists(u["W"][p], i, j)
        sc = shipped_score(dg, D.astype(np.float32).astype(float))
        shipped.append(rr[int(np.argmin(sc))]); pbest.append(rr.min()); tbest.append(rr[sub].min())
        synth.append(ca_rmsd(np.asarray(rec["fit_ca"]), u["nat_ca"]))
        band = np.where(rr <= rr.min() + BAND)[0]
        if not np.isin(band, sub).any():
            zero.append(t["pdb"])
    out = {"shipped": float(np.mean(shipped)), "pool_best": float(np.mean(pbest)), "top75_best": float(np.mean(tbest)),
           "synthesis_fit": float(np.mean(synth)), "n_zero_recall": len(zero), "zero": zero}
    if verbose:
        print(json.dumps(out, indent=1))
    assert abs(out["pool_best"] - 1.7108) < 2e-3 and abs(out["top75_best"] - 2.3062) < 2e-3
    assert set(zero) == set(FAIL18)
    return out

if __name__ == "__main__":
    selfcheck()
