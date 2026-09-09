"""Is POOL CONSTRUCTION a lever, when the selector is nearly powerless?

THE HYPOTHESIS
Candidates are retrieved into a pool, then one is picked by ranking with the learned
distogram. On the 126-target instrument the pool BEST is ~1.5 A and the SELECTED
candidate is ~3.45 A. `s7.debias` measured why: the predictor's deviations correlate
+0.283 with truth, and ranking by the pool's own mean profile with the sequence thrown
away entirely (`gain0`) selects 3.734 A against the predictor's 3.454 A. The whole
sequence-specific content of the shipped ranker is worth 0.28 A.

If the ranker is nearly powerless then the selected candidate is close to a DRAW from the
pool's own RMSD distribution. Two consequences follow, and they have never been separated
on this project:

  * Under an ORACLE ranker a bigger pool is monotonically better (K=100/500/2000 ->
    1.971/1.782/1.631 A), because a bigger pool has a better BEST.
  * Under a WEAK ranker a bigger pool need not be better and may be WORSE, because what
    a draw from the pool tracks is the pool's MEAN, and retrieval puts its best matches
    first: extending K adds candidates the ranker cannot discriminate, i.e. mostly new
    ways to pick badly.

If the second is what happens, the lever is pool CONSTRUCTION, not pool ranking: build a
smaller, higher-precision pool with a good MEAN and accept a weak selector. That needs no
fix to the predictor and no new model.

What is measured
1. `kcurve` selected RMSD vs K in {25 .. 2000} under the REAL shipped predictor, with
   pool best, pool mean, the oracle-ranked selection, and the ranker's percentile skill
   (where in the pool's own RMSD distribution the selected candidate lands; 50 = a coin).
2. `prune`  native-free pruning rules -- typicality, retrieval score, predicted-confidence
   agreement, medoid tightness, farthest-point diversity, and combinations -- each
   reported with pool best AND pool mean alongside selected RMSD, because pruning that
   improves the mean while destroying the best caps everything downstream.
3. `retr`   BLOSUM-retrieved pool vs an equal-size RANDOM real-fragment pool over the
   same window universe. This tests the recorded claim that random beats retrieval.

DEPLOYABLE vs DIAGNOSTIC -- the split is structural
`select()` and every `prune_*` function take only (pool distances, retrieval score,
predicted distribution). No native coordinate, distance, contact, torsion or RMSD can
reach them; `test_poolsize.py` asserts it by NaN-poisoning `Dnat` and checking the
selected index is bit-identical. Natives are read afterwards, in `_stats()`, to report
the CA-RMSD of whatever was already chosen.

The oracle curve is DIAGNOSTIC ONLY. It lives in `_oracle_select()`, is the only function
in this file that takes `Dnat`, is never called from a pruning or selection path, and is
labelled as an upper bound wherever it is printed.

PROTOCOL
Everything is swept on the 126-target instrument built by `s7.debias` (cluster
representatives, length 9-16, identity-cluster-disjoint from BOTH the 24-target dev set
and the 60-target benchmark; SE 0.147 A against dev's 0.354 A). The 24-target dev set is
touched ONCE, at the end, for the single best configuration only.

    python -m s7.poolsize cache    # per-target predictor arrays (no training)
    python -m s7.poolsize rand     # equal-size random real-fragment pools
    python -m s7.poolsize kcurve   # selected vs K, real / oracle / best / mean
    python -m s7.poolsize prune    # the native-free pruning table
    python -m s7.poolsize retr     # BLOSUM vs random retrieval key
    python -m s7.poolsize dev      # ONE pass on the 24 dev targets
    python -m s7.poolsize report   # print everything already computed
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:
    pass

import distogram as dgm                                                    # noqa: E402
import peptide_db as db                                                    # noqa: E402
from s7 import audit                                                       # noqa: E402
from s7 import debias                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PCACHE = os.path.join(HERE, "poolsize_cache")     #: predictor arrays, per target
RCACHE = os.path.join(HERE, "poolsize_rand")      #: random real-fragment pools

KMAX = audit.KMAX                                 #: 2000, what the pool cache holds
KS = (25, 50, 100, 250, 500, 1000, 2000)
K_REF = debias.K                                  #: 500, the pool every prior number uses
GRID = np.arange(2.0, 40.0, 0.05)                 #: the shipped Bayes-risk grid
RAND_SEED = 20260903

#: pruning sweep: (base pool size, kept size)
PRUNE_BASE = (500, 2000)
PRUNE_KEEP = (25, 50, 100, 250, 500)
CONF_Q = 0.25                                     #: fraction of most-confident pairs used


# ===================================================================== predictor cache
def _pred_path(pdbid):
    return os.path.join(PCACHE, f"{pdbid}.npz")


def _guard_esm():
    """Refuse to materialise the 1.5 GB `esm_cache.npz`.

    `dgm.train_fold` probes feature width with ONE arbitrary training sequence, which is
    not in the 21 MB hot cache; `esm_features.raw` would fall back to the full cache and
    blow the memory budget. Every fold model is already on disk, so nothing trains and the
    probe's VALUES are never used -- only its shape, which fixes `d_in`. So serve the
    probe zeros of the right shape and let every real (cached) sequence through unchanged.
    A sequence we actually predict on that is missing from the hot cache is a hard error,
    never a silent zero.
    """
    import esm_features as ef
    for f in range(5):
        p = dgm._model_path(f, True, True, 0)
        if not os.path.exists(p):
            raise RuntimeError(f"missing fold model {p}; this stage must not train")
    hot = ef._load_small()
    real = _real_seqs()
    _orig = ef.raw

    def raw(sequence):
        if sequence in hot:
            return _orig(sequence)
        if sequence in real:
            raise RuntimeError(f"sequence {sequence!r} absent from the ESM hot cache")
        return np.zeros((len(sequence), 1280), np.float32), \
            np.zeros((len(sequence), len(sequence)), np.float32)

    ef.raw = raw
    return _orig


_REAL = [None]


def _real_seqs():
    """Sequences whose predictions are actually used; these must have true features."""
    if _REAL[0] is None:
        s = {p.seq for p in debias.tuning_targets()} | {p.seq for p in debias.dev_targets()}
        _REAL[0] = s
    return _REAL[0]


def stage_cache(dev=False):
    """Per-target predictor arrays. Reads fold models off disk; trains nothing."""
    os.makedirs(PCACHE, exist_ok=True)
    _guard_esm()
    tg = debias.dev_targets() if dev else debias.tuning_targets()
    todo = [p for p in tg if (dev or debias._has(p.pdb)) and not os.path.exists(_pred_path(p.pdb))]
    print(f"{len(todo)} predictor caches to build ({'dev' if dev else 'tuning'})", flush=True)
    for k, p in enumerate(todo):
        rec = debias._load(p.pdb, tune=not dev)
        dg = debias.distogram_for(rec)
        np.savez_compressed(
            _pred_path(p.pdb) + ".tmp.npz",
            risk=dg._risk.astype(np.float32), expected=dg.expected.astype(np.float32),
            sd=dg.sd.astype(np.float32), w=dg.w.astype(np.float32),
            i=dg.i.astype(np.int32), j=dg.j.astype(np.int32))
        os.replace(_pred_path(p.pdb) + ".tmp.npz", _pred_path(p.pdb))
        print(f"[{k+1:3d}/{len(todo)}] {p.pdb}", flush=True)


def load_pred(pdbid):
    z = np.load(_pred_path(pdbid))
    return {k: z[k] for k in z.files}


# ===================================================================== random pools
def _rand_path(pdbid):
    return os.path.join(RCACHE, f"{pdbid}.npz")


def stage_rand(dev=False):
    """A RANDOM real-fragment pool over the same window universe, same size as BLOSUM's.

    Identical library, identical windowing, identical KMAX -- the ONLY change is that the
    KMAX candidates are drawn uniformly at random instead of by BLOSUM similarity. That
    isolates the retrieval key from everything else.
    """
    os.makedirs(RCACHE, exist_ok=True)
    folds = db.folds(5)
    tg = debias.dev_targets() if dev else debias.tuning_targets()
    tg = [p for p in tg if dev or debias._has(p.pdb)]
    todo = [p for p in tg if not os.path.exists(_rand_path(p.pdb))]
    if not todo:
        print("random pools already built")
        return
    frag_memo = {}
    _orig = dgm._fold_fragments

    def memo(fold, n_folds, threshold=db.IDENTITY_THRESHOLD):
        key = (fold, n_folds, threshold)
        if key not in frag_memo:
            frag_memo[key] = _orig(fold, n_folds, threshold)
        return frag_memo[key]

    dgm._fold_fragments = memo
    try:
        print(f"{len(todo)} random pools to build", flush=True)
        for k, p in enumerate(sorted(todo, key=lambda x: folds[x.seq])):
            t0 = time.time()
            fold, peps, frags = audit.build_pool_members(p, folds)
            W, S, src = audit.windows_of(peps + frags, p.n)
            sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
            rng = np.random.default_rng(RAND_SEED + abs(hash(p.pdb)) % 100000)
            idx = rng.choice(len(W), size=min(KMAX, len(W)), replace=False)
            i, j = audit.pair_index(p.n)
            np.savez_compressed(
                _rand_path(p.pdb) + ".tmp.npz",
                rr=audit.kabsch_rmsd_batch(W[idx], p.ca).astype(np.float32),
                D=audit.pair_dists(W[idx], i, j).astype(np.float32),
                Dnat=audit.pair_dists(p.ca[None], i, j)[0].astype(np.float32),
                sim=sim[idx].astype(np.float32), n_windows=np.int64(len(W)))
            os.replace(_rand_path(p.pdb) + ".tmp.npz", _rand_path(p.pdb))
            print(f"[{k+1:3d}/{len(todo)}] {p.pdb} windows={len(W):6d} "
                  f"best={np.sort(audit.kabsch_rmsd_batch(W[idx], p.ca))[0]:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    finally:
        dgm._fold_fragments = _orig


def load_pool(pdbid, kind="blosum", tune=True):
    """``(D, rr, sim)`` for one target. `rr` is native RMSD and is for REPORTING only."""
    if kind == "rand":
        z = np.load(_rand_path(pdbid))
        return z["D"].astype(float), z["rr"].astype(float), z["sim"].astype(float)
    rec = debias._load(pdbid, tune=tune)
    return rec["D"].astype(float), rec["rr"].astype(float), rec["sim"].astype(float)


# ===================================================================== DEPLOYABLE path
# Nothing below this line until the DIAGNOSTIC banner may touch a native quantity.
def select(pred, D):
    """The shipped ranker's choice. Index into `D` of the minimum Bayes-risk candidate."""
    return int(np.argmin(debias.score_risk(pred["risk"], GRID, D)))


def _typicality(D):
    """Mean |d - pool mean profile| per candidate. Pure pool composition, no predictor."""
    return np.abs(D - D.mean(0)[None, :]).mean(1)


def _conf_disagreement(pred, D, q=CONF_Q):
    """Weighted |d - predicted| over only the `q` most CONFIDENT pairs (lowest dist_sd)."""
    sd = pred["sd"]
    m = max(3, int(round(q * len(sd))))
    keep = np.argsort(sd)[:m]
    return np.abs(D[:, keep] - pred["expected"][None, keep]).mean(1)


def _medoid(D):
    c = D.mean(0)
    return int(np.argmin(np.abs(D - c[None, :]).mean(1)))


def _fps(D, m, seed=0):
    """Greedy farthest-point subset: a maximally DIVERSE pool of size `m`."""
    B = len(D)
    if m >= B:
        return np.arange(B)
    start = _medoid(D)
    picked = [start]
    dmin = np.abs(D - D[start][None, :]).mean(1)
    for _ in range(m - 1):
        nxt = int(np.argmax(dmin))
        picked.append(nxt)
        dmin = np.minimum(dmin, np.abs(D - D[nxt][None, :]).mean(1))
    return np.array(picked)


def prune(name, D, sim, pred, m, seed=0):
    """Return the indices of a size-<= `m` sub-pool. NATIVE-FREE by construction.

    Every criterion is a function of the candidate distance matrices, the retrieval score,
    and the predicted distribution -- all of which exist before any native is consulted.
    """
    B = len(D)
    if m >= B:
        return np.arange(B)
    if name == "none":
        return np.arange(B)
    if name == "sim":                       # keep the best-retrieved (= K truncation)
        return np.argsort(-sim, kind="stable")[:m]
    if name == "typ":                       # most typical: tight around the consensus
        return np.argsort(_typicality(D))[:m]
    if name == "atyp":                      # least typical: the outliers
        return np.argsort(-_typicality(D))[:m]
    if name == "medoid":                    # nearest the medoid candidate
        k = _medoid(D)
        return np.argsort(np.abs(D - D[k][None, :]).mean(1))[:m]
    if name == "fps":                       # maximally diverse
        return _fps(D, m, seed)
    if name == "conf":                      # best agreement on the confident pairs
        return np.argsort(_conf_disagreement(pred, D))[:m]
    if name == "typ+sim":                   # retrieve 2m, then keep the m most typical
        s = np.argsort(-sim, kind="stable")[:min(B, 2 * m)]
        return s[np.argsort(_typicality(D[s]))[:m]]
    if name == "conf+sim":
        s = np.argsort(-sim, kind="stable")[:min(B, 2 * m)]
        return s[np.argsort(_conf_disagreement(pred, D[s]))[:m]]
    if name == "typ+conf":                  # halve by typicality, then by confidence
        s = np.argsort(_typicality(D))[:min(B, 2 * m)]
        return s[np.argsort(_conf_disagreement(pred, D[s]))[:m]]
    if name == "typ_iter":                  # re-estimate the consensus as it tightens
        idx = np.arange(B)
        while len(idx) > m:
            k = max(m, len(idx) // 2)
            idx = idx[np.argsort(_typicality(D[idx]))[:k]]
        return idx
    if name == "rand":                      # the control
        return np.random.default_rng(seed).choice(B, size=m, replace=False)
    raise ValueError(name)


PRUNERS = ("none", "sim", "typ", "atyp", "medoid", "fps", "conf",
           "typ+sim", "conf+sim", "typ+conf", "typ_iter", "rand")


# ===================================================================== DIAGNOSTIC path
# Everything below reads a NATIVE. None of it is deployable and none of it feeds back
# into pool construction or ranking; it exists to draw the upper bound the real curve is
# compared against, and to report RMSD after a choice has already been made.
def _oracle_select(D, Dnat):
    """Diagnostic upper bound. Rank by agreement with the native distance matrix."""
    return int(np.argmin(np.abs(D - Dnat[None, :]).mean(1)))


def _stats(rr, sel):
    """REPORTING ONLY. `rr` is native CA-RMSD, read after the choice was made."""
    s = float(rr[sel])
    return {"sel": s, "best": float(rr.min()), "mean": float(rr.mean()),
            "median": float(np.median(rr)),
            "q10": float(np.percentile(rr, 10)), "q25": float(np.percentile(rr, 25)),
            "q30": float(np.percentile(rr, 30)),
            "pct": float(100.0 * (rr < s).mean()), "n": int(len(rr))}


# ===================================================================== statistics
def paired(a, b):
    return debias.paired(a, b)


def _agg(rows, key):
    v = np.array([r[key] for r in rows], float)
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))


# ===================================================================== stage: kcurve
def stage_kcurve(dev=False):
    tg = [p for p in (debias.dev_targets() if dev else debias.tuning_targets())
          if (dev or debias._has(p.pdb)) and os.path.exists(_pred_path(p.pdb))]
    print(f"K curve on {len(tg)} targets", flush=True)
    per = []
    for k, p in enumerate(tg):
        D, rr, sim = load_pool(p.pdb, "blosum", tune=not dev)
        Dnat = debias._load(p.pdb, tune=not dev)["Dnat"].astype(float)
        pred = load_pred(p.pdb)
        row = {"pdb": p.pdb, "n": int(p.n)}
        for K in KS:
            if K > len(D):
                continue
            Dk, rk = D[:K], rr[:K]
            real = _stats(rk, select(pred, Dk))
            orac = _stats(rk, _oracle_select(Dk, Dnat))       # DIAGNOSTIC
            row[str(K)] = {"real": real, "oracle_diagnostic": orac}
        per.append(row)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tg)}", flush=True)
    summary = []
    for K in KS:
        rows = [r[str(K)] for r in per if str(K) in r]
        if not rows:
            continue
        s = {"K": K, "n_targets": len(rows)}
        s["sel"], s["sel_se"] = _agg([r["real"] for r in rows], "sel")
        s["pool_best"], _ = _agg([r["real"] for r in rows], "best")
        s["pool_mean"], _ = _agg([r["real"] for r in rows], "mean")
        s["pool_median"], _ = _agg([r["real"] for r in rows], "median")
        s["pct"], s["pct_se"] = _agg([r["real"] for r in rows], "pct")
        s["q30"], _ = _agg([r["real"] for r in rows], "q30")
        s["oracle_diagnostic"], _ = _agg([r["oracle_diagnostic"] for r in rows], "sel")
        # Does the selection track the pool's BEST or the pool's BULK? Across targets,
        # regress the selected RMSD on each. A weak ranker returns a draw from the bulk.
        sel = np.array([r["real"]["sel"] for r in rows])
        s["corr_sel_pool_best"] = float(np.corrcoef(
            sel, [r["real"]["best"] for r in rows])[0, 1])
        s["corr_sel_pool_mean"] = float(np.corrcoef(
            sel, [r["real"]["mean"] for r in rows])[0, 1])
        s["corr_sel_pool_q30"] = float(np.corrcoef(
            sel, [r["real"]["q30"] for r in rows])[0, 1])
        ref = [r for r in per if str(K_REF) in r and str(K) in r]
        s["vs_K500"] = paired([r[str(K)]["real"]["sel"] for r in ref],
                              [r[str(K_REF)]["real"]["sel"] for r in ref])
        summary.append(s)
    out = {"n_targets": len(per), "instrument": "dev24" if dev else "tuning126",
           "summary": summary, "per_target": per}
    audit.write_json("poolsize_kcurve_dev.json" if dev else "poolsize_kcurve.json", out)
    _print_kcurve(out)
    return out


def _print_kcurve(out):
    print(f"\n=== SELECTED RMSD vs POOL SIZE ({out['n_targets']} targets, "
          f"{out['instrument']}) ===")
    print(f"{'K':>5} {'sel':>7} {'+-se':>6} {'d(K=500)':>9} {'ci95':>16} {'W/L':>9} "
          f"{'pool_best':>10} {'pool_mean':>10} {'pct':>6} | {'ORACLE*':>8}")
    for s in out["summary"]:
        v = s["vs_K500"]
        print(f"{s['K']:5d} {s['sel']:7.3f} {s['sel_se']:6.3f} {v['mean_diff']:+9.3f} "
              f"[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] {v['n_better']:3d}/{v['n_worse']:<5d} "
              f"{s['pool_best']:10.3f} {s['pool_mean']:10.3f} {s['pct']:6.1f} | "
              f"{s['oracle_diagnostic']:8.3f}")
    print("* ORACLE ranks by the native distance matrix: a DIAGNOSTIC upper bound, "
          "not deployable.")
    print("pct = percentile of the selected candidate in the pool's own RMSD "
          "distribution (50 = a coin flip).")
    print(f"\n--- does the answer track the pool's BEST or its BULK? ---")
    print(f"{'K':>5} {'sel':>7} {'pool q30':>9} {'sel-q30':>8} | "
          f"{'corr(sel,best)':>14} {'corr(sel,mean)':>14} {'corr(sel,q30)':>14}")
    for s in out["summary"]:
        print(f"{s['K']:5d} {s['sel']:7.3f} {s['q30']:9.3f} {s['sel']-s['q30']:+8.3f} | "
              f"{s['corr_sel_pool_best']:14.3f} {s['corr_sel_pool_mean']:14.3f} "
              f"{s['corr_sel_pool_q30']:14.3f}")


# ===================================================================== stage: prune
def stage_prune(dev=False):
    tg = [p for p in (debias.dev_targets() if dev else debias.tuning_targets())
          if (dev or debias._has(p.pdb)) and os.path.exists(_pred_path(p.pdb))]
    cfgs = [(b, m, c) for b in PRUNE_BASE for m in PRUNE_KEEP if m < b
            for c in PRUNERS if c != "none"]
    cfgs += [(b, b, "none") for b in PRUNE_BASE]
    print(f"{len(cfgs)} pruning configurations on {len(tg)} targets", flush=True)
    per = {f"{b}:{m}:{c}": [] for b, m, c in cfgs}
    for k, p in enumerate(tg):
        D, rr, sim = load_pool(p.pdb, "blosum", tune=not dev)
        pred = load_pred(p.pdb)
        for b, m, c in cfgs:
            Db, rb, sb = D[:b], rr[:b], sim[:b]
            idx = prune(c, Db, sb, pred, m, seed=k)
            assert len(idx) <= m and len(np.unique(idx)) == len(idx)
            st = _stats(rb[idx], select(pred, Db[idx]))
            st["pdb"] = p.pdb
            per[f"{b}:{m}:{c}"].append(st)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tg)}", flush=True)
    ref_key = f"{K_REF}:{K_REF}:none"
    ref = [r["sel"] for r in per[ref_key]]
    summary = []
    for key, rows in per.items():
        b, m, c = key.split(":")
        s = {"key": key, "base": int(b), "keep": int(m), "crit": c,
             "n_targets": len(rows)}
        s["sel"], s["sel_se"] = _agg(rows, "sel")
        s["pool_best"], _ = _agg(rows, "best")
        s["pool_mean"], _ = _agg(rows, "mean")
        s["pct"], _ = _agg(rows, "pct")
        s["vs_ref"] = paired([r["sel"] for r in rows], ref)
        summary.append(s)
    summary.sort(key=lambda s: s["sel"])
    # THE number for this experiment. Across every configuration, how much selected RMSD
    # does a unit improvement in the pool's MEAN buy, and how much does a unit loss of the
    # pool's BEST cost? If pool composition were a lever, the first slope would be near 1.
    ref_s = next(s for s in summary if s["key"] == ref_key)
    dm = np.array([s["pool_mean"] - ref_s["pool_mean"] for s in summary])
    dbst = np.array([s["pool_best"] - ref_s["pool_best"] for s in summary])
    ds = np.array([s["sel"] - ref_s["sel"] for s in summary])
    A = np.column_stack([dm, dbst, np.ones(len(dm))])
    coef, *_ = np.linalg.lstsq(A, ds, rcond=None)
    lev = {"n_configs": len(summary),
           "slope_sel_per_pool_mean": float(np.polyfit(dm, ds, 1)[0]),
           "slope_sel_per_pool_best": float(np.polyfit(dbst, ds, 1)[0]),
           "joint_slope_mean": float(coef[0]), "joint_slope_best": float(coef[1]),
           "corr_dsel_dmean": float(np.corrcoef(dm, ds)[0, 1]),
           "corr_dsel_dbest": float(np.corrcoef(dbst, ds)[0, 1]),
           "best_pool_mean_gain": float(-dm.min()),
           "sel_at_best_mean_gain": float(ds[int(np.argmin(dm))]),
           "n_significantly_better": int(sum(s["vs_ref"]["ci95"][1] < 0 for s in summary))}
    # The slopes above are dominated by `atyp`, which wrecks the mean AND the best at once.
    # Restrict to the configurations a pool-construction lever would actually use: those
    # that IMPROVE the pool mean. This is the honest test of the hypothesis.
    m = dm < -0.05
    if m.sum() > 2:
        lev["improving"] = {
            "n": int(m.sum()), "max_mean_gain": float(-dm[m].min()),
            "mean_dsel": float(ds[m].mean()), "best_dsel": float(ds[m].min()),
            "worst_dsel": float(ds[m].max()),
            "slope_per_mean": float(np.polyfit(dm[m], ds[m], 1)[0]),
            "corr_per_mean": float(np.corrcoef(dm[m], ds[m])[0, 1]),
            "slope_per_best": float(np.polyfit(dbst[m], ds[m], 1)[0]),
            "corr_per_best": float(np.corrcoef(dbst[m], ds[m])[0, 1]),
            "corr_meangain_skillloss": float(np.corrcoef(
                -dm[m], np.array([s["pct"] for s in summary])[m] - ref_s["pct"])[0, 1])}
    out = {"n_targets": len(tg), "ref": ref_key, "lever": lev,
           "instrument": "dev24" if dev else "tuning126",
           "summary": summary, "per_target": per}
    audit.write_json("poolsize_prune_dev.json" if dev else "poolsize_prune.json", out)
    _print_prune(out)
    return out


def _print_prune(out, top=None):
    print(f"\n=== NATIVE-FREE POOL PRUNING ({out['n_targets']} targets, "
          f"{out['instrument']}; reference = {out['ref']}) ===")
    print(f"{'base':>5} {'keep':>5} {'criterion':10} {'sel':>7} {'+-se':>6} {'d(ref)':>8} "
          f"{'ci95':>16} {'W/L':>9} {'pool_best':>10} {'pool_mean':>10} {'pct':>6}")
    rows = out["summary"][:top] if top else out["summary"]
    for s in rows:
        v = s["vs_ref"]
        print(f"{s['base']:5d} {s['keep']:5d} {s['crit']:10} {s['sel']:7.3f} "
              f"{s['sel_se']:6.3f} {v['mean_diff']:+8.3f} "
              f"[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] "
              f"{v['n_better']:3d}/{v['n_worse']:<5d} {s['pool_best']:10.3f} "
              f"{s['pool_mean']:10.3f} {s['pct']:6.1f}")
    lv = out.get("lever")
    if lv:
        print(f"\n--- what does pool composition BUY? ({lv['n_configs']} configurations) ---")
        print(f"  d(selected) per d(pool mean), alone : {lv['slope_sel_per_pool_mean']:+.3f} "
              f"(r = {lv['corr_dsel_dmean']:+.3f})")
        print(f"  d(selected) per d(pool best), alone : {lv['slope_sel_per_pool_best']:+.3f} "
              f"(r = {lv['corr_dsel_dbest']:+.3f})")
        print(f"  jointly: mean {lv['joint_slope_mean']:+.3f}, "
              f"best {lv['joint_slope_best']:+.3f}")
        print(f"  the configuration that improves the pool MEAN most improves it by "
              f"{lv['best_pool_mean_gain']:.3f} A and moves selected RMSD by "
              f"{lv['sel_at_best_mean_gain']:+.3f} A")
        print(f"  configurations significantly better than the reference "
              f"(95% CI excludes 0): {lv['n_significantly_better']}")
        im = lv.get("improving")
        if im:
            print(f"\n  restricted to the {im['n']} configurations that IMPROVE the pool "
                  f"mean (by up to {im['max_mean_gain']:.3f} A):")
            print(f"    mean d(selected) {im['mean_dsel']:+.3f} A   "
                  f"best {im['best_dsel']:+.3f}   worst {im['worst_dsel']:+.3f}")
            print(f"    d(sel)/d(pool mean) {im['slope_per_mean']:+.3f} "
                  f"(r {im['corr_per_mean']:+.3f})   "
                  f"d(sel)/d(pool best) {im['slope_per_best']:+.3f} "
                  f"(r {im['corr_per_best']:+.3f})")
            print(f"    corr(pool-mean gain, percentile-skill LOSS) "
                  f"{im['corr_meangain_skillloss']:+.3f}  -- pruning toward typicality "
                  f"eats exactly the skill the ranker had")


# ===================================================================== stage: retr
def stage_retr(dev=False):
    """BLOSUM-retrieved pool vs a RANDOM real-fragment pool of equal size."""
    tg = [p for p in (debias.dev_targets() if dev else debias.tuning_targets())
          if (dev or debias._has(p.pdb)) and os.path.exists(_pred_path(p.pdb))
          and os.path.exists(_rand_path(p.pdb))]
    print(f"retrieval-key comparison on {len(tg)} targets", flush=True)
    per = []
    for k, p in enumerate(tg):
        pred = load_pred(p.pdb)
        Dnat = debias._load(p.pdb, tune=not dev)["Dnat"].astype(float)
        row = {"pdb": p.pdb}
        for kind in ("blosum", "rand"):
            D, rr, sim = load_pool(p.pdb, kind, tune=not dev)
            for K in KS:
                if K > len(D):
                    continue
                Dk, rk = D[:K], rr[:K]
                st = _stats(rk, select(pred, Dk))
                st["oracle_diagnostic"] = float(rk[_oracle_select(Dk, Dnat)])
                row[f"{kind}:{K}"] = st
        per.append(row)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tg)}", flush=True)
    summary = []
    for K in KS:
        keys = [f"blosum:{K}", f"rand:{K}"]
        rows = [r for r in per if all(x in r for x in keys)]
        if not rows:
            continue
        s = {"K": K, "n_targets": len(rows)}
        for kind in ("blosum", "rand"):
            rs = [r[f"{kind}:{K}"] for r in rows]
            s[f"{kind}_sel"], s[f"{kind}_sel_se"] = _agg(rs, "sel")
            s[f"{kind}_best"], _ = _agg(rs, "best")
            s[f"{kind}_mean"], _ = _agg(rs, "mean")
            s[f"{kind}_pct"], _ = _agg(rs, "pct")
            s[f"{kind}_oracle_diagnostic"], _ = _agg(rs, "oracle_diagnostic")
        s["rand_vs_blosum"] = paired([r[f"rand:{K}"]["sel"] for r in rows],
                                     [r[f"blosum:{K}"]["sel"] for r in rows])
        summary.append(s)
    out = {"n_targets": len(per), "instrument": "dev24" if dev else "tuning126",
           "summary": summary, "per_target": per}
    audit.write_json("poolsize_retr_dev.json" if dev else "poolsize_retr.json", out)
    _print_retr(out)
    return out


def _print_retr(out):
    print(f"\n=== RETRIEVAL KEY: BLOSUM vs RANDOM real fragments "
          f"({out['n_targets']} targets, {out['instrument']}) ===")
    print(f"{'K':>5} | {'blosum sel':>10} {'best':>7} {'mean':>7} {'ORC*':>7} "
          f"| {'rand sel':>9} {'best':>7} {'mean':>7} {'ORC*':>7} "
          f"| {'d(rand-blos)':>12} {'W/L':>9}")
    for s in out["summary"]:
        v = s["rand_vs_blosum"]
        print(f"{s['K']:5d} | {s['blosum_sel']:10.3f} {s['blosum_best']:7.3f} "
              f"{s['blosum_mean']:7.3f} {s['blosum_oracle_diagnostic']:7.3f} "
              f"| {s['rand_sel']:9.3f} {s['rand_best']:7.3f} {s['rand_mean']:7.3f} "
              f"{s['rand_oracle_diagnostic']:7.3f} | {v['mean_diff']:+12.3f} "
              f"{v['n_better']:3d}/{v['n_worse']:<5d}")
    print("* ORC = oracle-ranked, DIAGNOSTIC upper bound only.")


# ===================================================================== stage: dev
def stage_dev():
    """ONE pass on the 24 dev targets, for the single best tuning configuration."""
    tune = audit.read_json("poolsize_prune.json")
    kc = audit.read_json("poolsize_kcurve.json")
    best = tune["summary"][0]
    best_k = min(kc["summary"], key=lambda s: s["sel"])
    print(f"tuning winner (pruning): {best['key']}  sel {best['sel']:.3f}")
    print(f"tuning winner (K only) : K={best_k['K']}  sel {best_k['sel']:.3f}")
    stage_cache(dev=True)
    tg = [p for p in debias.dev_targets() if os.path.exists(_pred_path(p.pdb))]
    rows = []
    for k, p in enumerate(tg):
        D, rr, sim = load_pool(p.pdb, "blosum", tune=False)
        pred = load_pred(p.pdb)
        r = {"pdb": p.pdb}
        Dref = D[:K_REF]
        r["ref"] = _stats(rr[:K_REF], select(pred, Dref))
        Db, rb, sb = D[:best["base"]], rr[:best["base"]], sim[:best["base"]]
        idx = prune(best["crit"], Db, sb, pred, best["keep"], seed=k)
        r["best"] = _stats(rb[idx], select(pred, Db[idx]))
        Dk = D[:best_k["K"]]
        r["bestK"] = _stats(rr[:best_k["K"]], select(pred, Dk))
        rows.append(r)
    out = {"n_targets": len(rows), "winner": best["key"], "winner_K": best_k["K"],
           "per_target": rows}
    for tag in ("ref", "best", "bestK"):
        s = {"sel": _agg([r[tag] for r in rows], "sel")[0],
             "sel_se": _agg([r[tag] for r in rows], "sel")[1],
             "pool_best": _agg([r[tag] for r in rows], "best")[0],
             "pool_mean": _agg([r[tag] for r in rows], "mean")[0],
             "pct": _agg([r[tag] for r in rows], "pct")[0]}
        if tag != "ref":
            s["vs_ref"] = paired([r[tag]["sel"] for r in rows],
                                 [r["ref"]["sel"] for r in rows])
        out[tag] = s
    audit.write_json("poolsize_dev.json", out)
    _print_dev(out)
    return out


def _print_dev(out):
    print(f"\n=== DEV SET ({out['n_targets']} targets) -- ONE pass, no iteration ===")
    print(f"winner from tuning: pruning {out['winner']}, K-only {out['winner_K']}")
    print(f"{'arm':22} {'sel':>7} {'+-se':>6} {'d(ref)':>8} {'ci95':>16} {'W/L':>9} "
          f"{'pool_best':>10} {'pool_mean':>10} {'pct':>6}")
    for tag, name in (("ref", f"K={K_REF} full pool"), ("best", out["winner"]),
                      ("bestK", f"K={out['winner_K']} only")):
        s = out[tag]
        v = s.get("vs_ref")
        d = (f"{v['mean_diff']:+8.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] "
             f"{v['n_better']:3d}/{v['n_worse']:<5d}") if v else f"{'--':>8} {'--':>16} {'--':>9}"
        print(f"{name:22} {s['sel']:7.3f} {s['sel_se']:6.3f} {d} "
              f"{s['pool_best']:10.3f} {s['pool_mean']:10.3f} {s['pct']:6.1f}")


# ===================================================================== stage: report
def stage_report():
    for name, printer in (("poolsize_kcurve.json", _print_kcurve),
                          ("poolsize_retr.json", _print_retr),
                          ("poolsize_prune.json", _print_prune),
                          ("poolsize_dev.json", _print_dev)):
        try:
            printer(audit.read_json(name))
        except FileNotFoundError:
            print(f"({name} not run)")


STAGES = {"cache": stage_cache, "rand": stage_rand, "kcurve": stage_kcurve,
          "prune": stage_prune, "retr": stage_retr, "dev": stage_dev,
          "report": stage_report}


def main(argv):
    if len(argv) < 2 or argv[1] not in STAGES:
        print(__doc__)
        print("stages: " + ", ".join(STAGES))
        return 1
    STAGES[argv[1]]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
