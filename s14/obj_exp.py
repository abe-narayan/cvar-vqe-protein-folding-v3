"""SPRINT 14, OBJ -- the experiment runner: arms, nulls, ablations, learning curve, OOD.

Every arm is evaluated the same way, on targets held out of its own fit, and reports
OBJECTIVE quality (rho with RMSD, globally and in-band; in-band pairwise ordering
accuracy below 1.5 A) and STRUCTURAL quality (mean RMSD of the low-energy decile and of
the top-100, minus that target's own space mean) on separate axes.

    python -m s14.obj_exp arms      [--loo|--lfo] [--targets A B C]
    python -m s14.obj_exp curve
    python -m s14.obj_exp ood
    python -m s14.obj_exp flip
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                    # noqa: E402
from s13.qarch_lib import spearman                 # noqa: E402
from s14 import obj_enum as E                      # noqa: E402
from s14 import obj_model as M                     # noqa: E402
from s14 import obj_train as T                     # noqa: E402
from s14 import obj_floor as F                     # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
NEVAL = 100_000                     # evaluation subsample; full space used for headline
LAMS = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]


def eval_idx(en, seed=11):
    """Fixed uniform evaluation subsample -- the same configurations for every arm."""
    if en.B <= NEVAL:
        return np.arange(en.B, dtype=np.int64)
    return np.sort(np.random.default_rng(seed).choice(en.B, NEVAL, replace=False))


def score_stats(e, r, rng, snap=None):
    o = np.argsort(e, kind="mergesort")
    dec = o[:max(50, len(e) // 10)]
    st = dict(rho_global=spearman(e, r), rho_decile=spearman(e[dec], r[dec]),
              rmsd_argmin=float(r[o[0]]), rmsd_decile_mean=float(r[dec].mean()),
              top100_mean_rmsd=float(r[o[:100]].mean()),
              top1000_mean_rmsd=float(r[o[:1000]].mean()),
              rmsd_space_mean=float(r.mean()), rmsd_space_min=float(r.min()))
    st["d_decile"] = st["rmsd_decile_mean"] - st["rmsd_space_mean"]
    st["d_top100"] = st["top100_mean_rmsd"] - st["rmsd_space_mean"]
    st["d_argmin"] = st["rmsd_argmin"] - st["rmsd_space_mean"]
    if snap is not None:
        st["snap_pct"] = float((e < e[snap]).mean())
    pi = F.pair_discrimination(e, r, rng, npair=300_000, sub=dec)
    st["inband_lt15"] = float(np.nanmean(
        [pi[k][0] for k in ("0.25-0.5", "0.5-1.0", "1.0-1.5") if k in pi]))
    st["inband_bins"] = {k: v[0] for k, v in pi.items()}
    pg = F.pair_discrimination(e, r, rng, npair=300_000)
    st["global_lt15"] = float(np.nanmean(
        [pg[k][0] for k in ("0.25-0.5", "0.5-1.0", "1.0-1.5") if k in pg]))
    st["global_bins"] = {k: v[0] for k, v in pg.items()}
    return st


def evaluate(pdb_id, w, fz, sequence_blind=False, mask=None, full=False, seed=0):
    en = T.Enum(pdb_id)
    idx = np.arange(en.B, dtype=np.int64) if full else eval_idx(en)
    S = T._decode(idx, en.n, en.k)
    ww = np.asarray(w)[:M.DIM].copy()
    if mask is not None:
        ww = ww * mask
    e = fz.score_with(ww, pdb_id, S, sequence_blind)
    r = en.rmsd[idx].astype(np.float64)
    snap = int(np.where(idx == en.snap)[0][0]) if (idx == en.snap).any() else None
    st = score_stats(e, r, np.random.default_rng(seed), snap)
    st["pdb"] = pdb_id
    st["fold"] = en.fold
    st["n"] = en.n
    return st


def baseline_stats(pdb_id, seed=0, full=False):
    en = T.Enum(pdb_id)
    idx = np.arange(en.B, dtype=np.int64) if full else eval_idx(en)
    r = en.rmsd[idx].astype(np.float64)
    snap = int(np.where(idx == en.snap)[0][0]) if (idx == en.snap).any() else None
    out = {}
    for nm, col in (("legacy", en.legacy), ("prior", en.prior)):
        st = score_stats(col[idx].astype(np.float64), r, np.random.default_rng(seed), snap)
        st["pdb"], st["fold"], st["n"] = pdb_id, en.fold, en.n
        out[nm] = st
    return out


# ---------------------------------------------------------------- feature masks
def feature_mask(part):
    m = np.zeros(M.DIM)
    if part in ("pair", "all"):
        m[:M.NPAIRF] = 1
    if part in ("loc", "all"):
        m[M.NPAIRF:M.NPAIRF + M.NLOCF] = 1
    if part in ("glob", "all"):
        m[M.NPAIRF + M.NLOCF:M.NPAIRF + M.NLOCF + M.NGLOB] = 1
    if part in ("inter", "all"):
        m[M.NPAIRF + M.NLOCF + M.NGLOB:M.DIM - 1] = 1
    m[M.DIM - 1] = 1
    return m


GLOB0 = M.NPAIRF + M.NLOCF


def glob_only_mask(names):
    m = np.zeros(M.DIM)
    for nm in names:
        m[GLOB0 + M.GLOB_NAMES.index(nm)] = 1
    return m


# ---------------------------------------------------------------- arms
def arm_fits(train, fz, arm, seed=0):
    """Return {lam: w} for one arm's training configuration."""
    kw = dict(seed=seed, center=True)
    if arm == "uniform":
        return T.fit(train, fz, n_uniform=20000, n_band=0, **kw)
    if arm == "uniform+band":
        return T.fit(train, fz, n_uniform=20000, n_band=20000, **kw)
    if arm == "band":
        return T.fit(train, fz, n_uniform=0, n_band=20000, **kw)
    if arm == "band01":
        return T.fit(train, fz, n_uniform=0, n_band=20000,
                     idx_fn=lambda en, rng: T.sample_idx(en, 0, 20000, rng,
                                                         band_frac=0.001), **kw)
    if arm == "seqblind":
        return T.fit(train, fz, n_uniform=20000, sequence_blind=True, **kw)
    if arm == "permuted":
        return T.fit(train, fz, n_uniform=20000, permute_label=True, **kw)
    if arm == "randfeat":
        return T.fit(train, fz, n_uniform=20000, random_feature=True, **kw)
    if arm == "leak":
        return T.fit(train, fz, n_uniform=20000, leak=True, **kw)
    raise ValueError(arm)


ARMS = ["uniform", "uniform+band", "band", "band01", "seqblind", "permuted",
        "randfeat", "leak"]
FIELDS = ["rho_global", "rho_decile", "d_decile", "d_top100", "inband_lt15",
          "global_lt15", "snap_pct"]


def pick_lam(train, fz, arm, seed=0, n_inner=3):
    """Select lambda by an inner leave-3-out INSIDE the training set. No test contact."""
    rng = np.random.default_rng(seed + 991)
    inner_val = list(rng.choice(train, min(n_inner, len(train) - 2), replace=False))
    inner_tr = [p for p in train if p not in inner_val]
    ws = arm_fits(inner_tr, fz, arm, seed=seed)
    best, bl = np.inf, LAMS[len(LAMS) // 2]
    sb = arm == "seqblind"
    for l in LAMS:
        if l not in ws:
            continue
        d = [evaluate(p, ws[l], fz, sequence_blind=sb)["d_decile"] for p in inner_val]
        if np.mean(d) < best:
            best, bl = float(np.mean(d)), l
    return bl, best


def run_arms(targets, arms=None, mode="loo", seed=0, full=False, out="obj_arms",
             fixed_lam=None):
    """`fixed_lam` skips nested selection -- use it for the diagnostic sweep, with a
    value PREDEFINED before the sweep.  The headline run uses nested selection."""
    arms = arms or ARMS
    fz = M.Featurizer(cache_dir=True)
    folds = {p: T.Enum(p).fold for p in targets}
    groups = ([(p, [q for q in targets if q != p]) for p in targets] if mode == "loo"
              else [(p, [q for q in targets if folds[q] != folds[p]]) for p in targets])
    res = {a: [] for a in arms}
    base = {"legacy": [], "prior": []}
    t0 = time.time()
    for p in targets:
        b = baseline_stats(p, full=full)
        for nm in base:
            base[nm].append(b[nm])
    print(f"baselines done [{time.time()-t0:.0f}s]", flush=True)
    fitcache = {}
    for p, train in groups:
        key = tuple(sorted(train))
        for a in arms:
            ck = (a, key)
            if ck not in fitcache:
                lam = (fixed_lam if fixed_lam is not None
                       else pick_lam(train, fz, a, seed=seed)[0])
                fitcache[ck] = (lam, arm_fits(train, fz, a, seed=seed)[lam])
            lam, w = fitcache[ck]
            st = evaluate(p, w, fz, sequence_blind=(a == "seqblind"), full=full)
            st["lam"] = lam
            st["arm"] = a
            res[a].append(st)
        print(f"{p} done [{time.time()-t0:.0f}s] "
              + " ".join(f"{a}:{res[a][-1]['d_decile']:+.2f}" for a in arms), flush=True)
    summarise(res, base, targets, out)
    return res, base


def summarise(res, base, targets, out="obj_arms"):
    folds = np.array([T.Enum(p).fold for p in targets])
    print(f"\n{'arm':14s} " + " ".join(f"{f:>12s}" for f in FIELDS) + "   W/L(dec)")
    rows = {}
    for nm, rr in list(base.items()) + list(res.items()):
        v = {f: np.array([s.get(f, np.nan) for s in rr], float) for f in FIELDS}
        wl = f"{int((v['d_decile'] < 0).sum())}/{int((v['d_decile'] > 0).sum())}"
        print(f"{nm:14s} " + " ".join(f"{np.nanmean(v[f]):12.3f}" for f in FIELDS)
              + f"   {wl}")
        rows[nm] = {f: v[f].tolist() for f in FIELDS}
    # paired comparisons of the structural axis against legacy and against random(0)
    print("\n-- paired d_decile vs LEGACY (negative = learned better) --")
    lg = np.array([s["d_decile"] for s in base["legacy"]])
    for nm, rr in res.items():
        a = np.array([s["d_decile"] for s in rr])
        pr = I.paired(a, lg, folds=folds, names=targets)
        print(f"{nm:14s} diff={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},"
              f"{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']} "
              f"drop10={pr['drop_top10_mean_diff']} folds={pr.get('per_fold')}")
    print("\n-- paired d_decile vs RANDOM (0 by construction) --")
    for nm, rr in list(base.items()) + list(res.items()):
        a = np.array([s["d_decile"] for s in rr])
        pr = I.paired(a, np.zeros_like(a), folds=folds, names=targets)
        print(f"{nm:14s} diff={pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},"
              f"{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']} "
              f"drop10={pr['drop_top10_mean_diff']}")
    I.write(out, {"targets": targets, "folds": folds.tolist(),
                  "arms": {k: v for k, v in res.items()},
                  "baselines": base}, n_expected=len(targets))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "arms"
    tg = [a for a in sys.argv[2:] if not a.startswith("-")]
    tg = tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])]
    lfo = "--lfo" in sys.argv
    if mode == "arms":
        # lam = 0.03 was PREDEFINED from the first exploratory leave-one-out sweep
        # (it maximised both mean rho_global and mean d_decile there) and is fixed for
        # every arm in this diagnostic, so no arm gets a selection advantage.
        run_arms(tg, arms=["uniform", "uniform+band", "band01", "seqblind",
                           "permuted", "leak"],
                 mode="lfo" if lfo else "loo", fixed_lam=0.03,
                 out="obj_arms_lfo" if lfo else "obj_arms_loo")
