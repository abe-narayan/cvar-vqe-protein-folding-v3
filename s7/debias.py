"""Can the shipped distogram be DE-BIASED post-hoc into a ranker? No retraining.

THE DIAGNOSIS THIS TESTS (from `s7.audit`, stage `prior`)
--------------------------------------------------------
Ranking a 500-candidate pool by the NATIVE distance matrix selects 1.782 A. The shipped
learned distogram selects 3.470 A at a dev MAE of 2.219 A. The audit swept three error
models injected into the oracle matrix and found that the SAME MAE buys wildly different
selection depending on the SHAPE of the error:

    iid Gaussian on pair distances      MAE 2.39 -> 2.257 A
    coordinate-space Gaussian           MAE 2.25 -> 2.677 A
    shrinkage toward the pool mean      MAE 2.31 -> 3.830 A   (a = 1.0)

The shipped prior lands on the SHRINKAGE curve essentially exactly (predicted 3.525 at its
own MAE) and is worse than iid noise of equal MAE on 23/24 targets. Inference: the
predictor is not imprecise, it is BIASED TOWARD THE MEAN -- its per-target distance matrix
has collapsed onto the average peptide's profile.

If that is right and the bias is invertible, then rescaling the predicted DEVIATION from
the pool's own mean profile should walk the predictor back down the shrinkage curve
(a ~ 1.0 -> a ~ 0.3 is 3.83 -> 2.21 A), at the cost of a WORSE MAE. MAE is therefore
reported and never steered by.

ARMS
----
1. `base`      the shipped Bayes-risk score, unchanged.
2. `gain{g}`   expected -> Dbar + g*(expected - Dbar), scored by weighted L1.
               `Dbar = D.mean(0)`, the mean distance matrix of the POOL'S OWN candidates:
               no native, no training-set estimate, leakage-free by construction.
               g = 1.0 is the same transform at unit gain and isolates the change of score
               FORM (weighted L1 instead of Bayes risk) from the change of gain.
               `gainVM` is the variance-matched gain, fitted once on the tuning instrument.
3. `pmi{lam}`  sum_p w_p [ -log P(d_p | seq, sep) + lam log P_marg(d_p | sep) ].
               P_marg is the separation-conditional marginal over the fold's own TRAINING
               data (out-of-fold peptides + `distogram._fold_fragments`), never a held-out
               native. lam = 0 is a pure weighted-NLL arm.
4. `sharp{T}`  prob -> prob^(1/T) renormalised, then the shipped Bayes-risk score.
               Two variants: `sharpW` recomputes the inverse-spread weight `w` from the
               sharpened distribution, `sharpF` holds `w` at its shipped value.

DEPLOYABLE vs DIAGNOSTIC -- the split is enforced structurally
--------------------------------------------------------------
`arm_scores()` is the DEPLOYABLE path. Its arguments are the predicted distribution, the
candidate pool's distance matrices, and a training-set marginal. It never receives a
held-out native's coordinates, distances, contacts or RMSD, and `test_debias.py` asserts
that by NaN-poisoning `Dnat` and checking every arm's scores are bit-identical.

Natives enter only afterwards, in `rank_stats()`, to REPORT the CA-RMSD of whatever was
already selected, and in the tuning-only fits (`calibration slope`, variance-matched gain)
which are computed on the 126-target tuning instrument and never on a dev target.

PROTOCOL
--------
Tuning instrument: 126 cluster-representative peptides, length 9-16, whose identity
clusters are disjoint from BOTH the 24-target dev set and the 60-target benchmark. Each is
scored by the model of a fold it is not in, over a pool built exactly as `s7.audit`
builds one. All sweeping happens here. The SINGLE best arm is then run on the 24 dev
targets ONCE, at the end.

    python -m s7.debias pools    # build the 126-target tuning instrument (~10 min)
    python -m s7.debias marg     # separation-conditional marginals, per fold
    python -m s7.debias calib    # THE discriminator: calibration slope on held-out pairs
    python -m s7.debias tune     # every arm on the tuning instrument
    python -m s7.debias dev      # baseline + the one winning arm, on the 24 dev targets
    python -m s7.debias report   # print everything already computed
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.stats import spearmanr

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:
    pass

import distogram as dgm                                                    # noqa: E402
import peptide_db as db                                                    # noqa: E402
from s7 import audit                                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "debias_cache")
K = audit.K_CLAIM                 #: 500, the pool the whole sprint reports on
BAND = 1.5
MIN_BAND = 5                      #: fewer in-band candidates than this and rho is not defined

#: The gain sweep runs BOTH ways. The diagnosis predicted only g > 1 (amplify a shrunken
#: deviation), but `stage_calib` measured a pooled calibration slope of +0.376 -- the
#: least-squares-optimal correction is to shrink the predicted deviation FURTHER, not to
#: amplify it, because the predicted deviations are already nearly full amplitude
#: (sd 2.58 A against a true 3.44 A) and merely uncorrelated with the truth (r = +0.28).
#: g = 0 is the no-sequence control: it ranks by the pool's own mean profile alone.
GAINS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0)
LAMS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25)
TEMPS = (1.0, 0.7, 0.5, 0.35, 0.25)

_EPS = 1e-9


# --------------------------------------------------------------------- target sets
def tuning_targets(n_max=None):
    """Cluster representatives disjoint from dev AND benchmark clusters. PINNED order.

    One peptide per identity cluster, so the instrument is not N views of one fold, and
    disjoint from the benchmark so nothing tuned here can leak into a reported number.
    """
    cl = db.clusters()
    bad = {cl[p.seq] for p in db.dev_set(24)} | {cl[p.seq] for p in db.benchmark()}
    out, seen = [], set()
    for p in sorted(db.load(), key=lambda x: x.pdb):
        if not (9 <= p.n <= 16):
            continue
        c = cl[p.seq]
        if c in bad or c in seen:
            continue
        seen.add(c)
        out.append(p)
    return out[:n_max] if n_max else out


def dev_targets():
    return db.dev_set(24)


# --------------------------------------------------------------------- pools
def stage_pools():
    """Audit-identical pool records for the tuning targets, with a per-fold memo."""
    os.makedirs(CACHE, exist_ok=True)
    folds = db.folds(5)
    frag_memo = {}
    _orig = dgm._fold_fragments

    def memo(fold, n_folds, threshold=db.IDENTITY_THRESHOLD):
        key = (fold, n_folds, threshold)
        if key not in frag_memo:
            frag_memo[key] = _orig(fold, n_folds, threshold)
        return frag_memo[key]

    dgm._fold_fragments = memo
    try:
        tg = tuning_targets()
        print(f"{len(tg)} tuning targets", flush=True)
        for k, p in enumerate(tg):
            path = os.path.join(CACHE, f"{p.pdb}.npz")
            if os.path.exists(path):
                continue
            t0 = time.time()
            rec = audit.build_target(p, folds)
            _save(rec)
            print(f"[{k+1:3d}/{len(tg)}] {rec['pdb']:6} n={rec['n']:2d} "
                  f"fold={rec['fold']} windows={rec['n_windows']:6d} "
                  f"best@{K}={rec['rr'][:K].min():.3f} ({time.time()-t0:.0f}s)", flush=True)
    finally:
        dgm._fold_fragments = _orig


def _save(rec):
    os.makedirs(CACHE, exist_ok=True)
    p = os.path.join(CACHE, f"{rec['pdb']}.npz")
    tmp = p + ".tmp.npz"
    np.savez_compressed(tmp, **{k: np.asarray(v) for k, v in rec.items()})
    os.replace(tmp, p)


def _load(pdbid, tune=True):
    path = os.path.join(CACHE if tune else audit.CACHE, f"{pdbid}.npz")
    z = np.load(path, allow_pickle=True)
    out = {k: z[k] for k in z.files}
    out["pdb"] = str(out["pdb"])
    out["seq"] = str(out["seq"])
    for k in ("n", "fold"):
        out[k] = int(out[k])
    return out


# --------------------------------------------------------------------- marginals
def marg_path(fold):
    return os.path.join(CACHE, f"marg_fold{fold}.npz")


def stage_marg():
    """P(bin | separation) from each fold's own TRAINING corpus. No held-out native."""
    os.makedirs(CACHE, exist_ok=True)
    folds = db.folds(5)
    for fold in range(5):
        if os.path.exists(marg_path(fold)):
            continue
        t0 = time.time()
        entries = ([p for p in db.load() if folds[p.seq] != fold]
                   + list(dgm._fold_fragments(fold, 5)))
        maxsep = 64
        cnt = np.zeros((maxsep + 1, dgm.NBINS), np.float64)
        npair = 0
        for e in entries:
            n = len(e.seq)
            if n < 3:
                continue
            i, j = np.triu_indices(n, k=2)
            d = np.linalg.norm(e.ca[i] - e.ca[j], axis=1)
            b = np.digitize(d, dgm.BIN_EDGES)
            s = np.minimum(j - i, maxsep)
            np.add.at(cnt, (s, b), 1.0)
            npair += len(d)
        np.savez_compressed(marg_path(fold) + ".tmp.npz", cnt=cnt,
                            n_entries=np.int64(len(entries)), n_pairs=np.int64(npair))
        os.replace(marg_path(fold) + ".tmp.npz", marg_path(fold))
        print(f"fold {fold}: {len(entries)} chains, {npair} pairs "
              f"({time.time()-t0:.0f}s)", flush=True)


_MARG = {}


def marginal(fold, alpha=1.0):
    """``(maxsep+1, NBINS)`` smoothed P(bin | sep). Rows sum to 1 by construction."""
    key = (fold, alpha)
    if key not in _MARG:
        z = np.load(marg_path(fold))
        c = z["cnt"] + alpha
        _MARG[key] = c / c.sum(1, keepdims=True)
    return _MARG[key]


# --------------------------------------------------------------------- scoring core
def bin_of(d):
    return np.digitize(np.asarray(d, float), dgm.BIN_EDGES)


def risk_table(prob, w, grid):
    """The shipped Bayes-risk lookup, rebuilt for an arbitrary (prob, w)."""
    r = (prob[:, None, :] * np.abs(grid[None, :, None]
                                   - dgm.CENTRES[None, None, :])).sum(2)
    return (r * w[:, None]).astype(np.float32)


def score_risk(risk, grid, D):
    """Shipped score evaluated straight off a pool's pair-distance matrix ``(B, npairs)``.

    Bit-identical to `Distogram.score(ca)` when D is that batch's pair distances -- the
    same clip-and-gather on the same 0.05 A grid. `test_debias.py` checks that.
    """
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32),
                0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


def score_wl1(D, target, w):
    """Weighted mean |d - target| over pairs."""
    return (np.abs(np.asarray(D, float) - target[None, :]) * w[None, :]).mean(1)


def sharpen(prob, T):
    if T == 1.0:
        return prob
    q = np.power(np.maximum(prob, 1e-30), 1.0 / T)
    return q / q.sum(1, keepdims=True)


def arm_scores(dg, D, fold, arms=None):
    """DEPLOYABLE. Every arm's score for a pool, from prediction + pool only.

    Parameters are the predicted distribution `dg` (sequence -> distances), the pool's own
    pair-distance matrix `D` (B, npairs), and the fold index used ONLY to pick the
    training-set marginal. No native quantity of the target reaches this function.
    Returns ``{arm_name: (B,) score, lower is better}``.
    """
    D = np.asarray(D, float)
    Dbar = D.mean(0)                       # the pool's own mean profile -- no native
    E, w, grid = dg.expected, dg.w, dg.grid
    sep = (dg.j - dg.i).astype(int)
    out = {}

    want = (lambda a: True) if arms is None else (lambda a: a in arms)

    if want("base"):
        out["base"] = score_risk(dg._risk, grid, D)

    for g in GAINS:
        name = f"gain{g:g}"
        if want(name):
            out[name] = score_wl1(D, Dbar + g * (E - Dbar), w)

    # Two principled single values, both scalars fitted once on the tuning instrument:
    # `gainVM` matches the variance of the predicted deviation to the native's, `gainOLS`
    # is the calibration slope itself (the least-squares-optimal rescaling).
    for name, g in (("gainVM", _fitted("gain_variance_matched")),
                    ("gainOLS", _fitted("slope_pooled"))):
        if g is not None and want(name):
            out[name] = score_wl1(D, Dbar + g * (E - Dbar), w)

    if any(want(f"pmi{l:g}") for l in LAMS):
        b = bin_of(D)                                              # (B, npairs)
        lp = np.log(np.maximum(dg.prob, _EPS))                     # (npairs, NBINS)
        nll = -lp[np.arange(len(sep))[None, :], b]                 # (B, npairs)
        M = marginal(fold)
        lm = np.log(np.maximum(M[np.minimum(sep, M.shape[0] - 1)], _EPS))
        lmarg = lm[np.arange(len(sep))[None, :], b]
        for lam in LAMS:
            name = f"pmi{lam:g}"
            if want(name):
                out[name] = ((nll + lam * lmarg) * w[None, :]).mean(1)

    for T in TEMPS:
        nw, nf = f"sharpW{T:g}", f"sharpF{T:g}"
        if not (want(nw) or want(nf)):
            continue
        q = sharpen(dg.prob, T)
        if want(nf):
            out[nf] = score_risk(risk_table(q, w, grid), grid, D)
        if want(nw):
            e2 = (q * dgm.CENTRES).sum(1)
            sd2 = np.sqrt(np.maximum(
                (q * (dgm.CENTRES[None] - e2[:, None]) ** 2).sum(1), 1e-6))
            w2 = dg._shell / np.maximum((sd2 + 0.5) ** dg._gamma, 1e-6)
            w2 = w2 / max(w2.mean(), 1e-12)
            out[nw] = score_risk(risk_table(q, w2, grid), grid, D)
    return out


def arm_expected(dg, D, fold):
    """The distance ESTIMATE each arm effectively ranks against, for MAE reporting."""
    Dbar = np.asarray(D, float).mean(0)
    E = dg.expected
    out = {"base": E}
    for name, g in (("gainVM", _fitted("gain_variance_matched")),
                    ("gainOLS", _fitted("slope_pooled"))):
        if g is not None:
            out[name] = Dbar + g * (E - Dbar)
    for g in GAINS:
        out[f"gain{g:g}"] = Dbar + g * (E - Dbar)
    for lam in LAMS:
        out[f"pmi{lam:g}"] = E
    for T in TEMPS:
        e2 = (sharpen(dg.prob, T) * dgm.CENTRES).sum(1)
        out[f"sharpW{T:g}"] = e2
        out[f"sharpF{T:g}"] = e2
    return out


_CAL = [None]


def _fitted(key):
    """A scalar fitted on the TUNING instrument by `stage_calib`, or None if unfitted.

    Both scalars this serves are tuning-only quantities: they are estimated from the 126
    instrument targets' natives and applied to dev as fixed constants, exactly as any
    hyperparameter would be.
    """
    if _CAL[0] is None:
        p = os.path.join(HERE, "debias_calib.json")
        _CAL[0] = json.load(open(p)) if os.path.exists(p) else {}
    return _CAL[0].get(key)


# --------------------------------------------------------------------- distograms
_MODELS = {}


def distogram_for(rec):
    """The fold model that never saw this target, applied to its sequence."""
    fold = rec["fold"]
    if fold not in _MODELS:
        _MODELS[fold] = dgm.train_fold(fold, True, 5, fragments=True, verbose=False)
    d = dgm.Distogram.for_target(rec["seq"], model=_MODELS[fold])
    # cache the pieces the sharpening arm needs to rebuild `w`
    fws, fg = dgm._score_weights()
    sep = (d.j - d.i).astype(float)
    shell = np.zeros(len(sep))
    for k, (a, b) in enumerate(dgm.SHELLS):
        shell[(sep >= a) & (sep <= b)] = fws[k]
    d._shell, d._gamma = shell, fg
    return d


# --------------------------------------------------------------------- statistics
def rank_stats(score, rr):
    """REPORTING ONLY. `rr` is the native CA-RMSD of each candidate, read after ranking."""
    out = {"sel": float(rr[int(np.argmin(score))]), "pool": float(rr.min()),
           "rho_global": float(spearmanr(score, rr).statistic)}
    m = rr <= rr.min() + BAND
    out["n_band"] = int(m.sum())
    out["rho_band"] = (float(spearmanr(score[m], rr[m]).statistic)
                       if m.sum() >= MIN_BAND else None)
    return out


def paired(a, b):
    """Paired difference a - b with a 95% t interval, and win/loss counts."""
    d = np.asarray(a, float) - np.asarray(b, float)
    n = len(d)
    se = float(d.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    return {"mean_diff": float(d.mean()), "se": se,
            "ci95": [float(d.mean() - 1.96 * se), float(d.mean() + 1.96 * se)],
            "n_better": int((d < -1e-12).sum()), "n_worse": int((d > 1e-12).sum()),
            "n_tied": int((np.abs(d) <= 1e-12).sum())}


def aggregate(rows, arm, base_arm="base"):
    sel = np.array([r[arm]["sel"] for r in rows])
    bse = np.array([r[base_arm]["sel"] for r in rows])
    rg = np.array([r[arm]["rho_global"] for r in rows])
    rb = np.array([r[arm]["rho_band"] if r[arm]["rho_band"] is not None else np.nan
                   for r in rows])
    nb = np.array([r[arm]["n_band"] for r in rows], float)
    ok = ~np.isnan(rb)
    out = {"arm": arm, "n_targets": len(rows),
           "sel": float(sel.mean()),
           "sel_se": float(sel.std(ddof=1) / np.sqrt(len(sel))),
           "pool": float(np.mean([r[arm]["pool"] for r in rows])),
           "mae": float(np.mean([r["mae"][arm] for r in rows])),
           "rho_global": float(rg.mean()),
           "rho_band": float(np.nanmean(rb)) if ok.any() else None,
           "rho_band_sizeweighted": (float((rb[ok] * nb[ok]).sum() / nb[ok].sum())
                                     if ok.any() else None),
           "n_band_scored": int(ok.sum()), "n_band_dropped": int((~ok).sum()),
           "band_size_median": float(np.median(nb))}
    out["vs_base"] = paired(sel, bse)
    return out


# --------------------------------------------------------------------- calibration
def stage_calib():
    """THE discriminator. Regress true deviation on predicted deviation, held-out pairs.

    y = d_true - Dbar, x = expected - Dbar, both taken pair-by-pair over the tuning
    instrument. A slope > 1 says the predicted deviations are systematically too SMALL --
    direct evidence of shrinkage -- and its value is the gain that would fix them in the
    least-squares sense. A slope ~ 1 says the model is not shrunk and the diagnosis behind
    this whole file is wrong.

    DIAGNOSTIC: it reads the tuning targets' natives. It is fitted on the instrument only
    and never on a dev target, and the fitted scalar is the only thing that crosses over.
    """
    recs = [_load(p.pdb) for p in tuning_targets() if _has(p.pdb)]
    print(f"calibrating on {len(recs)} tuning targets", flush=True)
    X, Y, SEP, per = [], [], [], []
    for rec in recs:
        D = rec["D"][:K].astype(float)
        Dn = rec["Dnat"].astype(float)
        dg = distogram_for(rec)
        assert len(dg.expected) == D.shape[1], (rec["pdb"], len(dg.expected), D.shape)
        Dbar = D.mean(0)
        x, y = dg.expected - Dbar, Dn - Dbar
        s = np.polyfit(x, y, 1) if len(x) > 2 else [np.nan, np.nan]
        per.append({"pdb": rec["pdb"], "n": rec["n"], "npairs": int(len(x)),
                    "slope": float(s[0]), "intercept": float(s[1]),
                    "sd_pred_dev": float(x.std()), "sd_true_dev": float(y.std()),
                    "mae": float(np.abs(dg.expected - Dn).mean()),
                    "corr": float(np.corrcoef(x, y)[0, 1])})
        X.append(x); Y.append(y); SEP.append((dg.j - dg.i).astype(int))
    X = np.concatenate(X); Y = np.concatenate(Y); SEP = np.concatenate(SEP)
    sl, ic = np.polyfit(X, Y, 1)
    gvm = float(Y.std() / max(X.std(), 1e-12))
    out = {"n_targets": len(recs), "n_pairs": int(len(X)),
           "slope_pooled": float(sl), "intercept_pooled": float(ic),
           "corr_pooled": float(np.corrcoef(X, Y)[0, 1]),
           "sd_pred_dev": float(X.std()), "sd_true_dev": float(Y.std()),
           "gain_variance_matched": gvm,
           "slope_per_target_median": float(np.median([p["slope"] for p in per])),
           "slope_per_target_iqr": [float(np.percentile([p["slope"] for p in per], 25)),
                                    float(np.percentile([p["slope"] for p in per], 75))],
           "n_targets_slope_gt1": int(sum(p["slope"] > 1 for p in per)),
           "mae_mean": float(np.mean([p["mae"] for p in per])),
           "by_shell": [], "per_target": per}
    for a, b in dgm.SHELLS:
        m = (SEP >= a) & (SEP <= b)
        if m.sum() > 10:
            s2, i2 = np.polyfit(X[m], Y[m], 1)
            out["by_shell"].append({"shell": [a, b], "n_pairs": int(m.sum()),
                                    "slope": float(s2), "intercept": float(i2),
                                    "sd_pred_dev": float(X[m].std()),
                                    "sd_true_dev": float(Y[m].std()),
                                    "corr": float(np.corrcoef(X[m], Y[m])[0, 1])})
    audit.write_json("debias_calib.json", out)
    _CAL[0] = out
    print(f"\npooled calibration slope {sl:+.3f}  (intercept {ic:+.3f}, "
          f"corr {out['corr_pooled']:+.3f})")
    print(f"variance-matched gain      {gvm:.3f}")
    print(f"per-target slope median    {out['slope_per_target_median']:+.3f} "
          f"IQR {out['slope_per_target_iqr'][0]:+.3f}..{out['slope_per_target_iqr'][1]:+.3f}"
          f"   ({out['n_targets_slope_gt1']}/{len(per)} above 1)")
    return out


def _has(pdbid):
    return os.path.exists(os.path.join(CACHE, f"{pdbid}.npz"))


# --------------------------------------------------------------------- sweeps
def run_arms(recs, tune, arms=None, tag=""):
    rows = []
    for k, rec in enumerate(recs):
        D = rec["D"][:K].astype(float)
        rr = rec["rr"][:K].astype(float)
        Dn = rec["Dnat"].astype(float)
        dg = distogram_for(rec)
        sc = arm_scores(dg, D, rec["fold"], arms=arms)
        ex = arm_expected(dg, D, rec["fold"])
        row = {"pdb": rec["pdb"], "n": rec["n"], "fold": rec["fold"],
               "mae": {a: float(np.abs(ex[a] - Dn).mean()) for a in sc}}
        for a, s in sc.items():
            row[a] = rank_stats(s, rr)
        rows.append(row)
        print(f"[{k+1:3d}/{len(recs)}] {rec['pdb']:6} n={rec['n']:2d} "
              f"pool {rr.min():5.3f} base {row['base']['sel']:5.3f}"
              + (f" {tag}" if tag else ""), flush=True)
    return rows


def stage_tune():
    recs = [_load(p.pdb) for p in tuning_targets() if _has(p.pdb)]
    rows = run_arms(recs, True)
    arms = [a for a in rows[0] if a not in ("pdb", "n", "fold", "mae")]
    summary = [aggregate(rows, a) for a in arms]
    summary.sort(key=lambda s: s["sel"])
    out = {"n_targets": len(rows), "K": K,
           "sel_sd_across_targets": float(np.std([r["base"]["sel"] for r in rows], ddof=1)),
           "summary": summary, "per_target": rows}
    audit.write_json("debias_tune.json", out)
    _print_table(summary, f"TUNING INSTRUMENT ({len(rows)} targets, K={K})")
    return out


def stage_dev():
    """ONE pass on the 24 dev targets: the shipped baseline and the single winning arm."""
    tu = audit.read_json("debias_tune.json")
    best = None
    for s in tu["summary"]:
        if s["arm"] != "base":
            best = s["arm"]
            break
    arms = ["base", best]
    print(f"dev run, arms {arms} (winner chosen on the tuning instrument)", flush=True)
    recs = [_load(p.pdb, tune=False) for p in dev_targets()]
    rows = run_arms(recs, False, arms=arms, tag=f"[{best}]")
    summary = [aggregate(rows, a) for a in arms]
    out = {"n_targets": len(rows), "K": K, "winning_arm": best,
           "summary": summary, "per_target": rows}
    audit.write_json("debias_dev.json", out)
    _print_table(summary, f"DEV SET ({len(rows)} targets, K={K}) -- one pass, no iteration")
    return out


def _print_table(summary, title):
    print(f"\n=== {title} ===")
    print(f"{'arm':10} {'sel':>7} {'+-se':>6} {'d(base)':>8} {'ci95':>16} "
          f"{'W/L':>8} {'MAE':>6} {'rho_g':>7} {'in-band':>8} {'wtd':>7} {'nb':>4}")
    for s in summary:
        v = s["vs_base"]
        rb = s["rho_band"]
        rw = s["rho_band_sizeweighted"]
        print(f"{s['arm']:10} {s['sel']:7.3f} {s['sel_se']:6.3f} {v['mean_diff']:+8.3f} "
              f"[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] "
              f"{v['n_better']:3d}/{v['n_worse']:<4d} {s['mae']:6.3f} "
              f"{s['rho_global']:+7.3f} "
              + (f"{rb:+8.3f} " if rb is not None else f"{'--':>8} ")
              + (f"{rw:+7.3f} " if rw is not None else f"{'--':>7} ")
              + f"{s['n_band_scored']:4d}")
    print(f"pool best (a perfect ranker) : {summary[0]['pool']:.3f} A")


def stage_report():
    try:
        c = audit.read_json("debias_calib.json")
        print("=== CALIBRATION SLOPE (regress true deviation on predicted deviation) ===")
        print(f"  {c['n_targets']} tuning targets, {c['n_pairs']} pairs")
        print(f"  pooled slope           {c['slope_pooled']:+.3f}   "
              f"(intercept {c['intercept_pooled']:+.3f}, corr {c['corr_pooled']:+.3f})")
        print(f"  sd(predicted deviation) {c['sd_pred_dev']:.3f} A vs "
              f"sd(true deviation) {c['sd_true_dev']:.3f} A")
        print(f"  variance-matched gain  {c['gain_variance_matched']:.3f}")
        print(f"  per-target median      {c['slope_per_target_median']:+.3f} "
              f"({c['n_targets_slope_gt1']}/{c['n_targets']} above 1)")
        print(f"\n  {'shell':>10} {'pairs':>8} {'slope':>7} {'sd_pred':>8} "
              f"{'sd_true':>8} {'corr':>7}")
        for s in c["by_shell"]:
            print(f"  {str(tuple(s['shell'])):>10} {s['n_pairs']:8d} {s['slope']:+7.3f} "
                  f"{s['sd_pred_dev']:8.3f} {s['sd_true_dev']:8.3f} {s['corr']:+7.3f}")
    except FileNotFoundError:
        print("(calib not run)")
    for name, title in (("debias_tune.json", "TUNING INSTRUMENT"),
                        ("debias_dev.json", "DEV SET")):
        try:
            d = audit.read_json(name)
        except FileNotFoundError:
            continue
        _print_table(d["summary"], f"{title} ({d['n_targets']} targets, K={d['K']})")


STAGES = {"pools": stage_pools, "marg": stage_marg, "calib": stage_calib,
          "tune": stage_tune, "dev": stage_dev, "report": stage_report}


def main(argv):
    if len(argv) < 2 or argv[1] not in STAGES:
        print(__doc__)
        print("stages: " + ", ".join(STAGES))
        return 1
    STAGES[argv[1]]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
