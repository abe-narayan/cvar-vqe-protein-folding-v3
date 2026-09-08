"""SPRINT 15, INFO -- shared helpers for the information-channel study.

Nothing here is an inference-time procedure by itself.  Any function whose name carries
`ORACLE` reads native data and is a diagnostic only.

Geometry caches are built straight off the cached enumerations (`s13/results/qarch_enum_*`
and `s14/cache/obj_enum_*`), which already carry the `PHI`/`PSI` (n, k) tables the
configurations decode against, so no library rebuild is needed and the traces are
bit-identical to the ones the enumerations were scored on.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402

CACHE = os.path.join(ROOT, "s15", "cache")
RESULTS = os.path.join(ROOT, "s15", "results")
FIG = os.path.join(ROOT, "s15", "figures")
for _d in (CACHE, RESULTS, FIG):
    os.makedirs(_d, exist_ok=True)


# ------------------------------------------------------------------ enumerated set
def enum_path(pdb):
    p = os.path.join(ROOT, "s14", "cache", f"obj_enum_{pdb}.npz")
    if os.path.exists(p):
        return p
    q = os.path.join(ROOT, "s13", "results", f"qarch_enum_{pdb}.npz")
    return q if os.path.exists(q) else None


def enum_targets():
    out = []
    for f in (glob.glob(os.path.join(ROOT, "s13", "results", "qarch_enum_*.npz")) +
              glob.glob(os.path.join(ROOT, "s14", "cache", "obj_enum_*.npz"))):
        out.append(os.path.basename(f).split("_")[-1][:-4])
    return sorted(set(out))


def decode(r, n, k):
    r = np.asarray(r, np.int64)
    out = np.empty((len(r), n), np.int8)
    for i in range(n):
        out[:, i] = (r // (k ** (n - 1 - i))) % k
    return out


# ------------------------------------------------------------------ shape scalars
def rg_of(X):
    """(B, n, 3) -> (B,) radius of gyration."""
    X = np.asarray(X, float)
    c = X - X.mean(-2, keepdims=True)
    return np.sqrt((c ** 2).sum(-1).mean(-1))


def shape_scalars(C, min_sep=3, cut=8.0):
    """Rg, end-to-end, contacts-per-residue for a batch of CA traces (B, n, 3)."""
    C = np.asarray(C, float)
    n = C.shape[1]
    rg = rg_of(C)
    e2e = np.linalg.norm(C[:, 0] - C[:, -1], axis=1)
    i, j = np.triu_indices(n, k=min_sep)
    d = np.linalg.norm(C[:, i, :] - C[:, j, :], axis=2)
    c8 = (d < cut).sum(1) / float(n)
    return rg, e2e, c8


def geom_cache(pdb, n_uniform=60000, band_max=2.5, n_band=30000, seed=0):
    """Per-target geometry on a uniform sample and on the near-native band.

    Returns dict with `u_*` (uniform sample) and `b_*` (band) arrays plus the native
    scalars (ORACLE, evaluation only).
    """
    path = os.path.join(CACHE, f"geom_{pdb}.npz")
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return {k: z[k] for k in z.files}
    z = np.load(enum_path(pdb))
    n, k = int(z["n"]), int(z["k"])
    PHI, PSI = np.asarray(z["PHI"], float), np.asarray(z["PSI"], float)
    rmsd = np.asarray(z["rmsd"], np.float32)
    B = len(rmsd)
    rng = np.random.default_rng(seed)
    ui = np.sort(rng.choice(B, min(n_uniform, B), replace=False))
    sel = np.flatnonzero(rmsd <= band_max)
    if len(sel) > n_band:
        sel = np.sort(rng.choice(sel, n_band, replace=False))
    rows = np.arange(n)

    def build(idx):
        S = decode(idx, n, k)
        phi = PHI[rows[None, :], S]
        psi = PSI[rows[None, :], S]
        out = []
        for a in range(0, len(idx), 20000):
            out.append(I.build_ca(phi[a:a + 20000], psi[a:a + 20000]))
        return np.concatenate(out, 0)

    Cu = build(ui)
    Cb = build(sel)
    ru, eu, cu = shape_scalars(Cu)
    rb, eb, cb = shape_scalars(Cb)
    nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)
    out = {"pdb": pdb, "n": n, "k": k, "fold": int(z["fold"]), "seq": str(z["seq"]),
           "u_idx": ui, "u_rmsd": rmsd[ui], "u_rg": ru, "u_e2e": eu, "u_c8": cu,
           "u_legacy": np.asarray(z["legacy"], np.float32)[ui],
           "u_prior": np.asarray(z["prior"], np.float32)[ui],
           "b_idx": sel, "b_rmsd": rmsd[sel], "b_rg": rb, "b_e2e": eb, "b_c8": cb,
           "b_legacy": np.asarray(z["legacy"], np.float32)[sel],
           "b_prior": np.asarray(z["prior"], np.float32)[sel],
           "rg_native_ORACLE": rg_of(nat[None])[0],
           "e2e_native_ORACLE": float(np.linalg.norm(nat[0] - nat[-1])),
           "band_max": band_max, "B": B, "rmsd_min_ORACLE": float(rmsd.min())}
    i2, j2 = np.triu_indices(n, k=3)
    d2 = np.linalg.norm(nat[i2] - nat[j2], axis=1)
    out["c8_native_ORACLE"] = float((d2 < 8.0).sum() / n)
    np.savez_compressed(path, **out)
    return out


# ------------------------------------------------------------------ statistics
def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a - a.mean(); b = b - b.mean()
    den = math.sqrt(float((a * a).sum()) * float((b * b).sum()))
    return float((a * b).sum() / den) if den > 0 else 0.0


def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = _rank(a); rb = _rank(b)
    return pearson(ra, rb)


def _rank(x):
    x = np.asarray(x, float)
    o = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    r[o] = np.arange(len(x), dtype=float)
    # average ties
    xs = x[o]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def boot_ci(a, b, fn=pearson, n_boot=4000, seed=0):
    a = np.asarray(a, float); b = np.asarray(b, float)
    rng = np.random.default_rng(seed)
    n = len(a)
    v = [fn(a[k], b[k]) for k in (rng.integers(0, n, n) for _ in range(n_boot))]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def perm_p(a, b, fn=pearson, n_perm=20000, seed=0, two_sided=True):
    rng = np.random.default_rng(seed)
    b = np.asarray(b, float)
    obs = fn(a, b)
    o = abs(obs) if two_sided else obs
    c = 0
    for _ in range(n_perm):
        v = fn(a, rng.permutation(b))
        c += (abs(v) >= o) if two_sided else (v >= o)
    return float((c + 1) / (n_perm + 1))


def residualise(x, *covs):
    x = np.asarray(x, float)
    A = np.column_stack([np.ones_like(x)] + [np.asarray(c, float) for c in covs])
    beta, *_ = np.linalg.lstsq(A, x, rcond=None)
    return x - A @ beta


def zscore(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / s if s > 0 else x * 0.0


# ------------------------------------------------------------------ concentration
def concentration_verdict(d, n_sim=4000, seed=0):
    """NULL-CALIBRATED concentration check.

    A raw drop-top-10 threshold is invalid when mean/sd is small.  This simulates a
    UNIFORM-EFFECT null with the same n, mean and sd, and asks whether the observed
    `top10_share` is extreme against it.  Emits one PASS/FAIL verdict, not fields.
    """
    d = np.asarray(d, float)
    n = len(d)
    if n <= 10 or d.sum() == 0:
        return {"verdict": "N/A", "reason": "n<=10 or zero total", "n": int(n)}
    o = np.argsort(d)
    share = float(d[o[:10]].sum() / d.sum())
    mean, sd = float(d.mean()), float(d.std(ddof=1))
    rng = np.random.default_rng(seed)
    sim = rng.normal(mean, sd, size=(n_sim, n))
    ss = np.sort(sim, axis=1)[:, :10].sum(1) / sim.sum(1)
    # concentration = observed share MORE extreme (further from 1) than uniform-effect null
    p = float((ss <= share).mean()) if mean < 0 else float((ss >= share).mean())
    ok = p > 0.05
    return {"verdict": "PASS" if ok else "FAIL", "top10_share": share,
            "null_share_median": float(np.median(ss)), "p_concentrated": p,
            "mean_over_sd": float(mean / sd) if sd > 0 else float("inf"),
            "power_note": "mean/sd small: test has little power" if abs(mean / sd) < 0.3
            else "", "n": int(n)}


def report(a, b, folds=None, names=None, label=""):
    """Full paired report with the null-calibrated concentration verdict attached."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    pr = I.paired(a, b, folds=folds, names=names)
    pr["concentration"] = concentration_verdict(a - b)
    pr["label"] = label
    return pr


def fmt_paired(pr):
    c = pr["concentration"]
    return (f"{pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] "
            f"med {pr['median_diff']:+.3f} W/L {pr['n_better']}/{pr['n_worse']} "
            f"conc {c['verdict']}(p={c.get('p_concentrated', float('nan')):.3f}, "
            f"m/sd={c.get('mean_over_sd', float('nan')):+.2f})")


def jwrite(name, obj):
    p = os.path.join(RESULTS, name if name.endswith(".json") else name + ".json")
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return p
