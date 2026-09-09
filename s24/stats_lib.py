"""s24/stats_lib.py -- THE SPRINT'S SHARED STATISTICS. OWNED BY LANE E. USE IT, DO NOT REIMPLEMENT.

Every lane reports the same quantities in the same way, so that two lanes' numbers can be put
side by side without an operator audit first. Import it:

    from s24 import stats_lib as ST
    r = ST.compare(new, old, folds, names=pdbs, label="residual-gen vs incumbent")
    print(ST.fmt(r))

`compare(a, b, ...)` is PAIRED and returns d = a - b, so **negative is better** for an RMSD
endpoint. State that on both sides of every contrast.

WHY EACH PIECE IS HERE (each one is a mistake this project has already made and paid for):

  MDE = 2.8016 * SE, PER COMPARISON.  `mde-is-per-comparison-not-per-instrument`: the quoted
  0.084 A constant was wrong by a factor of 84 in BOTH directions across 26 real comparisons.
  `sd` belongs to the comparison. Never quote a remembered MDE.

  FOLD-CLUSTERED BOOTSTRAP BESIDE THE IID CI.  Targets inside a fold share a distogram model,
  so the iid CI is anti-conservative. Both are reported; a result that needs the iid CI to
  exclude zero is not a result.

  TYPE-M FLAG AT 0.7-1.3x MDE.  An effect near its own MDE is significant only when the noise
  happened to help, so its magnitude is biased upward. `retrodesign` returns the exaggeration
  factor explicitly.

  MEDIAN BESIDE MEAN, AND WORST-TARGET DEGRADATION.  `median-vs-mean-is-the-free-warning`: a
  near-even W/L with a CI excluding zero is SUGGESTIVE of concentration. The warning is the
  GAP, compared against a uniform-effect null -- a raw drop-top threshold is not a valid test
  and has misfired here before, so `concentration()` supplies the null.

  BEST-OF-K AS THE DISTRIBUTION OF THE MAXIMUM.  Sprint 23 caught a "-0.077 A oracle" that was
  101% accounted for by its own best-of-K null. Any arm that takes the best of K things must
  be compared to `best_of_k_null`, which is the min/max order statistic, NOT the mean.

  TIE HANDLING.  `consensus-is-the-only-in-band-discriminator` (the methodological trap):
  `np.argmin` on a tied signal reads the array's ORDER, which on a sorted pool is the oracle
  order, and invented a 1.386 A winner. `argmin_tied` averages the outcome over the tied set.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

Z975, Z80 = 1.959963985, 0.841621234
#: the brief pins the constant at 2.8016 (= z_.975 + z_.80 to 5 s.f.).  Pinned literally so that
#: every lane's MDE is bit-identical to srcdecomp/biasalign, which hard-code the same literal.
MDE_K = 2.8016
NBOOT = 4000


def _rng(*parts):
    from s15 import seed as SD
    return SD.stable_rng(*parts, salt="s15")


def pinned_folds(pdbs=None):
    """The pinned fold assignment for the 126 dev targets, from the instrument itself.

    `benchmark-and-folds-must-be-pinned`: correcting the identity clustering silently moved 13
    benchmark targets and invalidated every fold model. Read the folds, never recompute them.
    """
    from s12 import instrument as I
    m = {t["pdb"]: int(t["fold"]) for t in I.targets()}
    if pdbs is None:
        return np.array([m[t["pdb"]] for t in I.targets()], int)
    return np.array([m[str(p)] for p in pdbs], int)


# ------------------------------------------------------------------ the paired comparison
def compare(a, b, folds=None, names=None, label="", n_boot=NBOOT, seed_parts=("s24stats",)):
    """Paired comparison of two arms over the SAME targets.  d = a - b; negative = a better.

    Returns every quantity the brief requires, in one dict.  `folds` may be omitted only when
    the comparison genuinely is not over the 126 dev targets; say so if you omit it.
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.shape != b.shape:
        raise ValueError("unpaired arms: %s vs %s" % (a.shape, b.shape))
    d = a - b
    n = len(d)
    if not np.isfinite(d).all():
        raise ValueError("non-finite paired differences: %d of %d" % ((~np.isfinite(d)).sum(), n))
    sd = float(d.std(ddof=1))
    se = sd / math.sqrt(n)
    mde = MDE_K * se
    rng = _rng(*seed_parts, label, "iid")
    bi = np.array([d[rng.integers(0, n, n)].mean() for _ in range(n_boot)])

    out = {
        "label": str(label), "n": int(n),
        "mean_a": float(a.mean()), "mean_b": float(b.mean()),
        "median_a": float(np.median(a)), "median_b": float(np.median(b)),
        "effect": float(d.mean()), "median_effect": float(np.median(d)),
        "sd": sd, "se": se, "mde": float(mde),
        "effect_over_mde": float(d.mean() / mde) if mde > 0 else float("nan"),
        "ci95_iid": [float(np.percentile(bi, 2.5)), float(np.percentile(bi, 97.5))],
        "n_better": int((d < 0).sum()), "n_worse": int((d > 0).sum()),
        "n_tied": int((d == 0).sum()),
        "worst_degradation": float(d.max()),
        "best_improvement": float(d.min()),
        "p90_degradation": float(np.percentile(d, 90)),
    }
    if names is not None:
        k = int(np.argmax(d))
        out["worst_target"] = str(names[k])
        o = np.argsort(d)
        out["top5_gain"] = [(str(names[i]), float(d[i])) for i in o[:5]]
        out["top5_loss"] = [(str(names[i]), float(d[i])) for i in o[-5:][::-1]]

    if folds is not None:
        folds = np.asarray(folds)
        if len(folds) != n:
            raise ValueError("folds length %d != n %d" % (len(folds), n))
        F = np.array(sorted(set(folds.tolist())))
        rf = _rng(*seed_parts, label, "fold")
        bf = np.array([np.concatenate([d[folds == q] for q in rf.choice(F, len(F), replace=True)]).mean()
                       for _ in range(n_boot)])
        out["ci95_fold"] = [float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))]
        out["per_fold"] = {int(q): float(d[folds == q].mean()) for q in F}
        out["n_folds"] = int(len(F))
        out["folds_same_sign"] = int(sum(np.sign(v) == np.sign(d.mean())
                                         for v in out["per_fold"].values()))
    else:
        out["ci95_fold"] = None
        out["per_fold"] = None

    rd = retrodesign(d.mean(), se)
    out.update({"power": rd["power"], "type_m": rd["type_m"], "type_s": rd["type_s"]})
    r = abs(out["effect_over_mde"])
    out["type_m_flag"] = bool(0.7 <= r <= 1.3)
    out["verdict"] = _verdict(out)
    out["concentration"] = concentration(d, seed_parts=seed_parts + (label,))
    return out


def _verdict(o):
    """A result must clear its own MDE AND have the FOLD CI exclude zero.  Nothing else counts.

    SIBLING OF THE UNDERPOWERED BUG, found by the s25 audit lane.  The MDE gate was added after
    this function had labelled a 0.39x-MDE effect "MEASURED" because a five-cluster bootstrap CI
    excluded zero.  The same shape survived one line lower: with `folds=None` the verdict fell
    back to the IID CI and could return BETTER/WORSE -- which this module's own docstring
    forbids ("a result that needs the iid CI to exclude zero is not a result").  A caller who
    simply omitted `folds` got the forbidden verdict silently.  It now refuses, and says what
    the iid CI would have claimed so nothing is hidden.
    """
    e, m = o["effect"], o["mde"]
    if abs(e) <= m:
        return "NOT MEASURED (|effect| %.4f <= its own MDE %.4f, %.2fx)" % (abs(e), m, abs(e) / m if m > 0 else float("nan"))
    if o["ci95_fold"] is None:
        ci = o["ci95_iid"]
        would = "WORSE" if ci[0] > 0 else "BETTER" if ci[1] < 0 else "nothing"
        return ("NOT MEASURED (no fold-clustered CI supplied; pass folds=. "
                "the iid CI alone would have said %s)" % would)
    ci = o["ci95_fold"]
    if ci[0] > 0:
        s = "WORSE"
    elif ci[1] < 0:
        s = "BETTER"
    else:
        return "NOT MEASURED (fold CI includes zero)"
    return s + (" [TYPE-M ZONE: magnitude inflated ~%.2fx]" % o["type_m"] if o["type_m_flag"] else "")


def retrodesign(effect, se, alpha=0.05):
    """Gelman/Carlin retrodesign: power, Type-S and Type-M for a TRUE effect equal to `effect`.

    Type-M is the expected |estimate| among significant estimates, divided by the true effect.
    An effect near its own MDE has Type-M well above 1: it is inflated by the very noise that
    made it significant.
    """
    from math import erf, sqrt, exp, pi
    if se <= 0:
        return {"power": float("nan"), "type_s": float("nan"), "type_m": float("nan")}
    ncp = abs(effect) / se
    crit = Z975

    def Phi(x):
        return 0.5 * (1.0 + erf(x / sqrt(2.0)))

    def phi(x):
        return exp(-0.5 * x * x) / sqrt(2.0 * pi)

    p_hi = 1.0 - Phi(crit - ncp)            # significant in the true direction
    p_lo = Phi(-crit - ncp)                 # significant in the WRONG direction
    power = p_hi + p_lo
    if power <= 0:
        return {"power": 0.0, "type_s": float("nan"), "type_m": float("nan")}
    # E|Z| restricted to |Z|>crit, Z ~ N(ncp,1); truncated-normal tail means
    e_hi = ncp * p_hi + phi(crit - ncp)
    e_lo = -(ncp * p_lo - phi(-crit - ncp))
    exp_abs = (e_hi + e_lo) / power
    return {"power": float(power), "type_s": float(p_lo / power),
            "type_m": float(exp_abs / ncp) if ncp > 0 else float("inf")}


# ------------------------------------------------------------------ concentration
def concentration(d, n_boot=NBOOT, seed_parts=("s24stats",)):
    """Is the effect carried by a few targets?  Compared against a UNIFORM-EFFECT null.

    `median-vs-mean-is-the-free-warning`: a raw drop-top threshold is NOT a valid test -- it
    misfired here once.  The null is: the same n targets, a uniform shift of the observed mean
    plus residuals resampled from the observed centred differences.  We report where the
    observed drop-top-10 mean sits in THAT null's distribution.
    """
    d = np.asarray(d, float); n = len(d)
    if n <= 12:
        return None
    o = np.argsort(d)
    obs = float(d[o[10:]].mean())
    mu = float(d.mean())
    resid = d - mu
    rng = _rng(*seed_parts, "conc")
    null = np.empty(n_boot)
    for t in range(n_boot):
        s = mu + resid[rng.integers(0, n, n)]
        null[t] = s[np.argsort(s)[10:]].mean()
    pct = float((null < obs).mean())
    return {"mean": mu, "median": float(np.median(d)),
            "drop_top10_mean": obs,
            "null_p10": float(np.percentile(null, 10)),
            "null_p50": float(np.percentile(null, 50)),
            "null_p90": float(np.percentile(null, 90)),
            "pctile_in_null": pct,
            "flag": bool(pct < 0.05 or pct > 0.95),
            "top10_share": float(d[o[:10]].sum() / d.sum()) if d.sum() != 0 else None}


# ------------------------------------------------------------------ the best-of-K null
def best_of_k_null(values, k, n_boot=NBOOT, minimise=True, seed_parts=("s24stats",)):
    """The distribution of the ORDER STATISTIC, not the mean.

    `values` is the population an arm draws from (e.g. the per-candidate RMSDs a selector
    chooses among, or the per-seed outcomes a "best seed" arm picks from).  Returns the
    distribution of min (or max) over `k` draws WITH replacement.

    Any arm that reports "the best of K" must be scored against this.  A -0.077 A "oracle"
    last sprint was 101% accounted for by exactly this null.
    """
    v = np.asarray(values, float).ravel()
    v = v[np.isfinite(v)]
    if len(v) == 0:
        raise ValueError("no finite values")
    rng = _rng(*seed_parts, "bok", str(k), str(minimise))
    draws = v[rng.integers(0, len(v), size=(n_boot, int(k)))]
    best = draws.min(1) if minimise else draws.max(1)
    return {"k": int(k), "n_pop": int(len(v)), "minimise": bool(minimise),
            "pop_mean": float(v.mean()), "pop_median": float(np.median(v)),
            "expected_best": float(best.mean()), "median_best": float(np.median(best)),
            "ci95": [float(np.percentile(best, 2.5)), float(np.percentile(best, 97.5))],
            "gain_vs_pop_mean": float(best.mean() - v.mean())}


def best_of_k_within(M, n_boot=400, seed_parts=("s24stats",)):
    """The RIGHT null for an oracle taken over K VARIANTS of the same targets.  s25 audit lane.

    `M` is (n_targets, K): every target scored under every one of K settings.  An arm that
    reports `M.min(1)` is a best-of-K arm and BRIEF SS4 requires it be scored against the
    distribution of the MINIMUM.  `best_of_k_null` above is the wrong tool here -- it resamples
    from a POOLED population, which is right for a candidate pool and wrong for a K-column grid.

    THE TRAP THIS EXISTS TO CLOSE.  s25/temper.py rolled its own null and drew, for each target,
    K indices into THAT TARGET'S OWN K columns, with replacement, then took the minimum.  Such a
    draw can only ever return one of the row's own K values and returns its true minimum with
    probability 1-(1-1/K)^K -- 0.656 at K=8.  It is bounded below by the observed statistic BY
    CONSTRUCTION and its "share accounted" is a near-constant function of K, not a measurement.
    It reported 86% and the true transfer was 4%.  That null is computed here too, labelled
    INVALID, so the trap is visible rather than re-inventable.

    Returns `share_accounted`, `residual`, `k_eff`, and -- the number to actually quote --
    `split_half`, from `split_half_transfer`, which nulls itself and needs no null at all.
    """
    M = np.asarray(M, float)
    if M.ndim != 2:
        raise ValueError("best_of_k_within wants (n_targets, K), got %s" % (M.shape,))
    n, K = M.shape
    rng = _rng(*seed_parts, "bokw")
    dev = M - M.mean(1, keepdims=True)
    observed = float((M.min(1) - M.mean(1)).mean())

    #: VALID: deviations resampled ACROSS targets, so any per-target structure in WHICH setting
    #: wins is destroyed while the marginal spread of the deviations is preserved.
    flat = dev.ravel()
    across = float(np.mean([flat[rng.integers(0, len(flat), (n, K))].min(1).mean()
                            for _ in range(n_boot)]))
    #: INVALID, reproduced only so it can be recognised: min of K with-replacement draws from
    #: the target's own row.
    own = float(np.mean([np.minimum.reduce([dev[np.arange(n), rng.integers(0, K, n)]
                                            for _ in range(K)]).mean() for _ in range(n_boot)]))
    sh = split_half_transfer(M, n_boot=n_boot, seed_parts=seed_parts)
    #: effective number of independent settings: exp of the entropy of how often each wins.
    cnt = np.bincount(M.argmin(1), minlength=K).astype(float)
    p = cnt / cnt.sum()
    k_eff = float(np.exp(-(p[p > 0] * np.log(p[p > 0])).sum()))
    return {"n": n, "k": K, "observed_gain": observed,
            "null_across_targets": across,
            "share_accounted": (across / observed) if observed else float("inf"),
            "residual": observed - across,
            "null_own_row_with_replacement_INVALID": own,
            "share_accounted_INVALID": (own / observed) if observed else float("inf"),
            "p_own_min_drawn": float(1.0 - (1.0 - 1.0 / K) ** K),
            "k_eff": k_eff, "argmin_counts": cnt.tolist(),
            "split_half": sh["transfer"], "split_half_frac": sh["frac_of_oracle"],
            "verdict": ("NOT A SIGNAL (split-half transfers %.0f%% of the oracle)"
                        % (100 * sh["frac_of_oracle"]) if abs(sh["frac_of_oracle"]) < 0.25
                        else "residual survives: %+.4f transfers of %+.4f oracle"
                             % (sh["transfer"], observed))}


def split_half_transfer(M, n_boot=400, seed_parts=("s24stats",)):
    """BRIEF SS4's preferred construction for a best-of-K arm: it nulls itself.

    Choose the setting on a random half of the targets, score it on the other half, average both
    directions, repeat.  If the winning setting is target-specific noise the transfer is ~0; no
    separate null is needed, which is exactly why the brief prefers it to a best-of-K null.
    """
    M = np.asarray(M, float)
    n, K = M.shape
    dev = M - M.mean(1, keepdims=True)
    rng = _rng(*seed_parts, "splithalf")
    obs = float((M.min(1) - M.mean(1)).mean())
    out = np.empty(n_boot)
    for t in range(n_boot):
        perm = rng.permutation(n); h1, h2 = perm[:n // 2], perm[n // 2:]
        out[t] = 0.5 * (dev[h2][:, int(np.argmin(dev[h1].mean(0)))].mean()
                        + dev[h1][:, int(np.argmin(dev[h2].mean(0)))].mean())
    return {"oracle_gain": obs, "transfer": float(out.mean()),
            "ci95": [float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))],
            "frac_of_oracle": float(out.mean() / obs) if obs else float("nan")}


def achieved(nominal, measured, label=""):
    """Does the operator do what its parameter is NAMED for?  s25 audit lane.

    The s25 SD(f) arm was labelled "multiply the posterior sd by f" and delivered 1.43x at f=1.5
    while ALSO shifting the posterior mean by -0.375 A -- a confound the size of the defect under
    study, invisible because only the nominal parameter was ever printed.  Print this beside every
    swept parameter: the arm's name is a claim and it is checkable.
    """
    a, b = np.asarray(nominal, float), np.asarray(measured, float)
    r = b / np.where(a == 0, np.nan, a)
    return {"label": str(label), "nominal": a.tolist(), "achieved": b.tolist(),
            "ratio": r.tolist(), "max_shortfall": float(np.nanmax(np.abs(1.0 - r))),
            "flag": bool(np.nanmax(np.abs(1.0 - r)) > 0.10)}


def best_of_k_accounted(observed_arm, baseline, values, k, **kw):
    """What share of a claimed best-of-K gain is the null alone?  >= ~1.0 means the arm is the null.

    observed_arm and baseline are scalars on the SAME basis (state it).  `values` is the
    population the arm searched.
    """
    nul = best_of_k_null(values, k, **kw)
    obs_gain = float(observed_arm - baseline)
    null_gain = float(nul["expected_best"] - baseline)
    return {"observed_gain": obs_gain, "null_gain": null_gain,
            "share_accounted": (null_gain / obs_gain) if obs_gain != 0 else float("inf"),
            "residual_gain": obs_gain - null_gain, "null": nul}


# ------------------------------------------------------------------ ties
def argmin_tied(score, outcome, atol=1e-12, rtol=1e-9):
    """Mean outcome over the ARGMIN SET, not the first index numpy happens to return.

    On a pool held in oracle-sorted order, `np.argmin` on a tied score reads the oracle order
    and manufactures a winner.  Returns (mean outcome over ties, size of the tie set).
    """
    s = np.asarray(score, float); y = np.asarray(outcome, float)
    m = np.nanmin(s)
    tie = np.isclose(s, m, atol=atol, rtol=rtol)
    return float(y[tie].mean()), int(tie.sum())


# ------------------------------------------------------------------ reporting
def fmt(o):
    """One block, identical across lanes.  Print this, not a hand-rolled line."""
    L = []
    L.append("  %s" % (o["label"] or "comparison"))
    L.append("    a %.4f (med %.4f)   b %.4f (med %.4f)   n=%d"
             % (o["mean_a"], o["median_a"], o["mean_b"], o["median_b"], o["n"]))
    L.append("    effect %+.4f   median %+.4f   SE %.4f   MDE %.4f   effect/MDE %+.2f"
             % (o["effect"], o["median_effect"], o["se"], o["mde"], o["effect_over_mde"]))
    L.append("    iid  CI95 [%+.4f, %+.4f]" % tuple(o["ci95_iid"]))
    if o["ci95_fold"] is not None:
        L.append("    fold CI95 [%+.4f, %+.4f]   folds same sign %d/%d   per-fold %s"
                 % (o["ci95_fold"][0], o["ci95_fold"][1], o["folds_same_sign"], o["n_folds"],
                    " ".join("%d:%+.3f" % (k, v) for k, v in sorted(o["per_fold"].items()))))
    L.append("    %dW/%dL/%dT   worst degradation %+.4f%s   p90 %+.4f   power %.2f  Type-M %.2f"
             % (o["n_better"], o["n_worse"], o["n_tied"], o["worst_degradation"],
                (" (%s)" % o["worst_target"]) if "worst_target" in o else "",
                o["p90_degradation"], o["power"], o["type_m"]))
    c = o.get("concentration")
    if c:
        L.append("    concentration: drop-top10 %+.4f vs uniform-effect null p10/p50/p90 "
                 "%+.4f/%+.4f/%+.4f -> pctile %.3f%s"
                 % (c["drop_top10_mean"], c["null_p10"], c["null_p50"], c["null_p90"],
                    c["pctile_in_null"], "  FLAG" if c["flag"] else ""))
    L.append("    VERDICT: %s" % o["verdict"])
    return "\n".join(L)


# ------------------------------------------------------------------ provenance
def provenance(module_file):
    """Stamp a result with the source that produced it.  `provenance(__file__)`.

    A Python process loads its source ONCE, at launch.  Sprint 24 caught `srcdecomp.json`'s P4
    row carrying numbers from `min(N, 2000)` under a script that says `min(N, 300)` and a label
    that says "300 draw" -- a job launched before an edit kept writing pre-edit results, and
    nothing in the file recorded which source made them.  BRIEF.md SS5 already demands the git
    commit; this makes it one call.  Put the result under a "provenance" key.
    """
    import datetime
    import hashlib
    import subprocess
    out = {"module": os.path.basename(str(module_file)),
           "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
    try:
        with open(module_file, "rb") as fh:
            out["source_sha256"] = hashlib.sha256(fh.read()).hexdigest()[:16]
    except Exception as e:
        out["source_sha256"] = "unavailable: %s" % e
    try:
        out["git_commit"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                           capture_output=True, text=True,
                                           timeout=20).stdout.strip()[:12]
        out["git_dirty"] = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                                               capture_output=True, text=True,
                                               timeout=20).stdout.strip())
    except Exception as e:
        out["git_commit"] = "unavailable: %s" % e
    return out


def save_atomic(path, obj, complete_keys=None, rows=None, n_expected=None, module_file=None):
    """Atomic write (tmp + os.replace) with a completion flag on the FULL KEY SET, not a count.

    `complete_keys` is the set of keys every row must carry; `complete` is True only when the row
    count matches AND every row carries every key.  A row count alone has passed a partial file
    off as a finished one in this project before.
    """
    import json as _json
    o = dict(obj)
    if rows is not None and n_expected is not None:
        ok = len(rows) == int(n_expected)
        if complete_keys:
            ok = ok and all(all(k in r for k in complete_keys) for r in rows)
        o["n_rows"] = len(rows); o["n_expected"] = int(n_expected); o["complete"] = bool(ok)
        o["complete_keys"] = list(complete_keys) if complete_keys else None
    if module_file is not None:
        o["provenance"] = provenance(module_file)
    t = str(path) + ".tmp"
    with open(t, "w") as fh:
        _json.dump(o, fh)
    os.replace(t, str(path))
    return path


def selftest():
    rng = np.random.default_rng(0)
    n = 126
    folds = np.repeat(np.arange(5), 26)[:n]
    b = rng.normal(3.05, 0.9, n)
    a = b + rng.normal(-0.05, 0.20, n)
    r = compare(a, b, folds, label="selftest: a true -0.05 shift")
    print(fmt(r))
    assert abs(r["mde"] - 2.8016 * r["se"]) < 1e-12
    r0 = compare(b, b.copy(), folds, label="selftest: identical arms")
    assert r0["effect"] == 0.0 and r0["verdict"].startswith("NOT MEASURED")
    pop = rng.normal(3.0, 0.8, 500)
    bk = best_of_k_null(pop, 20)
    assert bk["expected_best"] < bk["pop_mean"]
    print("\n  best-of-20 null over a 500-population: pop mean %.4f -> expected BEST %.4f "
          "(a 'gain' of %+.4f that is pure order statistic)"
          % (bk["pop_mean"], bk["expected_best"], bk["gain_vs_pop_mean"]))
    y, k = argmin_tied(np.array([1.0, 1.0, 1.0, 2.0]), np.array([9.0, 1.0, 5.0, 0.0]))
    assert k == 3 and abs(y - 5.0) < 1e-12
    print("  argmin_tied: 3-way tie -> mean outcome %.1f (np.argmin would have said 9.0)" % y)

    #: the UNDERPOWERED fix and its sibling.  Both must refuse.
    c = b + rng.normal(-0.3 * r["mde"], 0.20, n)          # a real but sub-MDE effect
    assert compare(c, b, folds, label="sub-MDE")["verdict"].startswith("NOT MEASURED (|effect|")
    assert compare(a, b, label="no folds")["verdict"].startswith("NOT MEASURED (no fold")
    assert compare(a, b, folds, label="with folds")["verdict"] != \
        compare(a, b, label="no folds")["verdict"]
    print("  MDE gate: a sub-MDE effect refuses; and folds=None now refuses instead of")
    print("            reading the verdict off the iid CI (the sibling bug).")

    #: best-of-K over a K-column grid of PURE NOISE.  The valid null must account for ~all of it
    #: and the split-half transfer must be ~0; the invalid own-row null must not.
    M = rng.normal(3.0, 0.5, (126, 1)) + rng.normal(0.0, 0.3, (126, 8))
    w = best_of_k_within(M, n_boot=200)
    print("\n  best_of_k_within on a PURE-NOISE 126x8 grid:")
    print("    observed oracle %+.4f   valid null %+.4f (%.0f%%)   split-half transfer %+.4f (%.0f%%)"
          % (w["observed_gain"], w["null_across_targets"], 100 * w["share_accounted"],
             w["split_half"], 100 * w["split_half_frac"]))
    print("    INVALID own-row-with-replacement null %+.4f (%.0f%%)  <- a near-constant of K,"
          % (w["null_own_row_with_replacement_INVALID"], 100 * w["share_accounted_INVALID"]))
    print("    P(own min drawn) = 1-(1-1/8)^8 = %.3f.  Never quote this null." % w["p_own_min_drawn"])
    assert abs(w["split_half_frac"]) < 0.25, "pure noise must not transfer"
    assert w["share_accounted"] > 0.8, "the valid null must account for a pure-noise oracle"

    ach = achieved([1.5, 3.0], [1.43, 2.32], label="s25 SD(f): nominal vs achieved sd ratio")
    assert ach["flag"]
    print("  achieved(): s25's SD(f) delivers %.2fx at f=1.5 and %.2fx at f=3.0 -> FLAG"
          % (1.43, 2.32))
    print("\n  stats_lib selftest OK")


if __name__ == "__main__":
    selftest()
