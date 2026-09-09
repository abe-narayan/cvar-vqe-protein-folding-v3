"""s23/qc_lib.py -- WORKSTREAM C shared plumbing: stats (iid + fold-clustered) and rg_z.

Pre-registered in `s23/PREREG_C.md` before either experiment (`c1_probweight.py`,
`c2_rgcond.py`) reads a single RMSD. Nothing here is fit to an outcome.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402

MDE_K = 2.8016  # BRIEF's fixed multiplier: MDE = 2.8016 * SE


# ============================================================================ atomic writes
def save_json(name: str, obj) -> str:
    path = os.path.join(RESULTS, name)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)
    return path


# ============================================================================ paired stats
def _boot_mean(d: np.ndarray, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    n = len(d)
    idx = rng.integers(0, n, size=(n_boot, n))
    return d[idx].mean(1)


def _boot_mean_fold(d: np.ndarray, folds: np.ndarray, rng: np.random.Generator,
                    n_boot: int) -> np.ndarray:
    """Cluster bootstrap: resample the fold IDs with replacement, pool every target in the
    resampled folds, and take the mean -- clusters (folds), not targets, are the resampling
    unit, per BRIEF's 'paired, fold-clustered CIs beside iid'."""
    uniq = np.unique(folds)
    out = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        vals = np.concatenate([d[folds == f] for f in pick])
        out[b] = vals.mean()
    return out


def paired_stats(a: np.ndarray, b: np.ndarray, folds: np.ndarray, n_boot: int = 4000,
                 seed: int = 0, label: str = "") -> dict:
    """`a - b`, negative = a better. iid + fold-clustered CI, SE, MDE, W/L. BRIEF's fixed unit."""
    a = np.asarray(a, float); b = np.asarray(b, float); folds = np.asarray(folds, int)
    d = a - b
    n = len(d)
    se = float(d.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
    mde = MDE_K * se
    rng = np.random.default_rng(seed)
    bs_iid = _boot_mean(d, rng, n_boot)
    bs_fold = _boot_mean_fold(d, folds, rng, n_boot)
    mean_diff = float(d.mean())
    return {
        "label": label, "n": int(n), "mean_diff": mean_diff, "se": se, "mde": float(mde),
        "abs_over_mde": float(abs(mean_diff) / mde) if mde > 0 else float("nan"),
        "ci95_iid": [float(np.percentile(bs_iid, 2.5)), float(np.percentile(bs_iid, 97.5))],
        "ci95_fold": [float(np.percentile(bs_fold, 2.5)), float(np.percentile(bs_fold, 97.5))],
        "n_better": int((d < 0).sum()), "n_worse": int((d > 0).sum()),
        "n_tied": int((d == 0).sum()),
        "worst_degradation": float(d.max()),
        "median_diff": float(np.median(d)),
    }


def verdict(ps: dict) -> str:
    """CI-excludes-zero AND past its own MDE, on BOTH iid and fold-clustered -- else NOT MEASURED."""
    lo_i, hi_i = ps["ci95_iid"]; lo_f, hi_f = ps["ci95_fold"]
    excl = (lo_i * hi_i > 0) and (lo_f * hi_f > 0)
    past_mde = ps["abs_over_mde"] >= 1.0
    if excl and past_mde:
        return "MEASURED"
    if excl and not past_mde:
        return "DIRECTION ESTABLISHED, MAGNITUDE NOT MEASURED"
    return "NULL / NOT MEASURED"


# ============================================================================ rg_z, full n=126
CA_BOND = 3.8046  # trans virtual CA-CA distance, s21/rgsign.py's own constant


def _rg2_from_pairs(d2_sum: float, n: int) -> float:
    return d2_sum / float(n * n)


def _rg_of(P: np.ndarray) -> float:
    P = np.asarray(P, float)
    c = P - P.mean(0, keepdims=True)
    return float(np.sqrt((c * c).sum() / len(P)))


def compute_rg_family(pdb: str, u=None) -> dict:
    """`s21/rgsign.py`'s exact recipe, on the FULL 126-target dev instrument (that file's own
    cache is n=75, D's n<=13 panel only) -- rg_disto, rg_pool_mean, rg_pool_sd, rg_emit, rg_gap,
    rg_z. Native-free except `rg_emit`, which is flagged partly circular exactly as upstream."""
    if u is None:
        u = I.load_univ(pdb)
    n = int(u["n"])
    p = I.pool_idx(u)
    W = np.asarray(u["W"], float)[p]
    dg = I.distogram(pdb, str(u["seq"]), int(u["fold"]))
    i, j = I.pair_index(n)
    dhat = np.asarray(dg["expected"], float)
    d2 = float((dhat ** 2).sum() + (n - 1) * CA_BOND ** 2)
    rg_disto = float(np.sqrt(_rg2_from_pairs(d2, n)))
    rgs = np.array([_rg_of(w) for w in W], float)
    D = I.pair_dists(W, i, j)
    sc = np.asarray(I.shipped_score(dg, D), float)
    order = np.argsort(sc, kind="stable")[:75]
    avg, _b = I.coordinate_average(W[order])
    return {
        "pdb": pdb, "n": n, "fold": int(u["fold"]),
        "rg_disto": rg_disto,
        "rg_pool_mean": float(rgs.mean()),
        "rg_pool_sd": float(rgs.std(ddof=1)),
        "rg_emit": _rg_of(np.asarray(avg, float)),          # PARTLY CIRCULAR, flagged
        "rg_gap": float(rg_disto - rgs.mean()),
        "rg_z": float((rg_disto - rgs.mean()) / max(rgs.std(ddof=1), 1e-9)),
    }


def build_rg_table(force: bool = False) -> list:
    """rg_z family for all 126 dev targets. Cached; recomputed only with `force=True`."""
    path = os.path.join(RESULTS, "rg_table126.json")
    if os.path.exists(path) and not force:
        return json.load(open(path))["rows"]
    rows = []
    tg = I.targets()
    print(f"rg_table126: {len(tg)} targets", flush=True)
    for c, t in enumerate(sorted(tg, key=lambda r: r["pdb"])):
        u = I.load_univ(t["pdb"])
        rows.append(compute_rg_family(t["pdb"], u))
        if (c + 1) % 20 == 0:
            print(f"  {c + 1}/{len(tg)}", flush=True)
    ok = len(rows) == len(tg)
    save_json("rg_table126.json",
              {"rows": rows, "complete": bool(ok), "n_expected": len(tg)})
    return rows
