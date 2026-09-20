#!/usr/bin/env python
"""s30/s30_D_meter.py -- lane D: THE S30 COST-RMSD METER (charter section 7).

This EXTENDS `s29/s29_D_cost_audit.py` rather than replacing it: the S29 module is imported and
its six diagnostics are re-run unchanged (they reproduce S28/S29's published anchors exactly --
see S30-L1), and everything this file adds is layered on top.  S29's artefacts are never
touched: the ladder cache is redirected to `s30/results/s30_D_ladder_structs/` and the results
to `s30/results/`.

WHAT S29 ALREADY HAD (unchanged, re-run through M):
  (a) LADDER rho(f, ORACLE RMSD) per target over three ladders, fold-clustered CI
  (b) COSINE cos(-grad f at production, direction to the native), rigid body removed
  (c) PCTILE the native's percentile in its own 500-member pool under f
  (d) PREF   share of targets on which f scores a rung below production, + pool-member control

WHAT S30 ADDS (the charter's section 7 list that S29's version lacked):
  (E1) PER-TARGET DISTRIBUTIONS as standard output -- min / p10 / q1 / median / q3 / p90 / max
       and the five extreme targets for every diagnostic.  A mean is not a distribution: two of
       this project's retracted claims were means over a bimodal per-target distribution.
  (E2) FAIL18 vs 108 as a fold-clustered CONTRAST with SE, MDE and a CI, not two bare means.
       WITH THE HONEST CAVEAT the split cannot avoid: fold 0 contains NO FAIL18 target
       (FAIL18 by fold = {1:6, 2:2, 3:4, 4:6, 0:0}), so the cluster bootstrap of the split
       contrast resamples from 4 non-empty clusters and is correspondingly weak.
  (E3) fold-clustered CIs on EVERY diagnostic, including the new ones.
  (E4) MATCHED RANDOM SIGNED STRUCTURES as a FIRST-CLASS control: R=8 independent draws per
       target (S29 used one, seed 0), so the control carries its own draw spread and a
       single-draw control can be seen for what it is.  Reported: pref(circ_best vs PROD) minus
       pref(RAND_SIGNED vs PROD) averaged over draws, the per-draw values, and the CHARTER
       ladder rho recomputed on every draw.
  (E5) AGGREGATE RANDOM-DIRECTION NULL for the cosine.  S29 printed the per-draw mean |cos|
       (0.140), which is the magnitude of ONE random direction, not the null for the
       126-target MEAN.  Both are printed here and labelled.
  (E6) a VERDICT with explicit, pre-stated gates and the 0.7x / 1.0x MDE rule, and
  (E7) a MULTIPLICITY REGISTER: every comparison this meter emits is counted, so the sprint
       coordinator can price a 1x-MDE positive against the number of looks taken.

USAGE
    python s30/s30_D_meter.py build-cache [--chain] [--slice i/n] [--limit N] [--rebuild]
    python s30/s30_D_meter.py meter --f DIS [--basis ca|chain|chain-s28rows] [--tag NAME]
    python s30/s30_D_meter.py meter --f s30.s30_X_cost:my_cost --basis chain
    python s30/s30_D_meter.py verify          # the six S29 anchors, asserted
    python s30/s30_D_meter.py selftest

EVERY NUMBER THE METER PRINTS IS ORACLE: it reads the native to label the ladder, to point the
gradient and to place the native in the pool.  The meter TUNES NOTHING.  A cost function is
given the poisoned context (ctx.nat_ca / ctx.oracle_rr are NaN); a cost that reads them emits
NaN and the meter refuses it.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s29 import s29_D_cost_audit as M      # noqa: E402  the S29 instrument, verified in S30-L1
from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_C2_recog_audit as C2   # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "s30_D_ladder_structs")
M.CACHE = CACHE                            #: S29's cache stays byte-identical

R_RAND = 8                                 #: matched random-signed draws per target (E4)
SALT = "s30D"

#: the shipped cost's metered values, the incumbent every new cost is read against.
#: VERIFIED by `verify` on 2026-09-20 (S30-L1); the source is S29-L2 / s29/results/.
ANCHORS = {
    "ladder_rho_S28_chain": -0.4023, "ladder_rho_S28_ca": -0.1818,
    "cosine": -0.0339, "cosine_random_ref_abs": 0.140,
    "native_pctile_DIS": 0.3676, "native_pctile_DIS_SURR": 0.3688,
    "pref_circ_best_chain": 0.0714, "pref_circ_best_ca": 0.2063,
    "pool_member_chain": 0.0201, "pool_member_ca": 0.1265,
}


def _safe(k):
    return k.replace("[", "_").replace("]", "")


def rs_key(s):
    return "RAND_SIGNED[%d]" % s


# ============================================================================ (E4) cache: extra draws
def extra_path(pdb):
    return M.cache_path(pdb)


def load_extra(pdb, need_chain=False):
    """The R_RAND-1 EXTRA matched random-signed draws (seeds 1..R-1); draw 0 is in M's rungs."""
    f = M.cache_path(pdb)
    if not os.path.exists(f):
        return None
    with np.load(f, allow_pickle=False) as z:
        keys = set(z.files)
        need = ["x_" + _safe(rs_key(s)) for s in range(1, R_RAND)]
        if any(k not in keys for k in need):
            return None
        if need_chain and any("xchain_ca_" + _safe(rs_key(s)) not in keys for s in range(1, R_RAND)):
            return None
        out = {"S": {rs_key(s): np.array(z["x_" + _safe(rs_key(s))]) for s in range(1, R_RAND)},
               "oracle_rmsd": {rs_key(s): float(z["xrmsd_" + _safe(rs_key(s))]) for s in range(1, R_RAND)}}
        if "xchain_ca_" + _safe(rs_key(1)) in keys:
            out["chain"] = {rs_key(s): {q: np.array(z[f"xchain_{q}_{_safe(rs_key(s))}"]) for q in ("ca", "phi", "psi")}
                            for s in range(1, R_RAND)}
            out["oracle_rmsd_chain"] = {rs_key(s): float(z["xrmsdchain_" + _safe(rs_key(s))]) for s in range(1, R_RAND)}
    return out


def _append_npz(path, d):
    """Merge new arrays into an existing .npz without losing what is there."""
    old = {}
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as z:
            old = {k: np.array(z[k]) for k in z.files}
    old.update(d)
    tmp = path + f".tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, **old)
    from s27 import s28_A_amp as A
    A.replace_retry(tmp, path)


def build_cache(pdbs, with_chain=False, rebuild=False):
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        have = M.load_rungs(pdb, need_chain=with_chain)
        have_x = load_extra(pdb, need_chain=with_chain)
        if have is not None and have_x is not None and not rebuild:
            continue
        t1 = time.time()
        cand, dis, top, dg, frame, sur = M.load_target(pdb)
        b = M.build_rungs(pdb, cand, frame, dg, top, with_chain=with_chain)
        M.save_rungs(pdb, cand, b)
        # ---- (E4) the extra matched random-signed draws, same generator, seeds 1..R-1
        d_prod = float(b["meta"]["d_circ_best_prod"])
        d = {}
        for s in range(1, R_RAND):
            Cs, _ = C2.rand_signed(frame, b["S"]["PROD"], d_prod, pdb, s)
            d["x_" + _safe(rs_key(s))] = np.asarray(Cs, float)
            d["xrmsd_" + _safe(rs_key(s))] = np.array(float(C2.rmsd_between(Cs, cand.nat_ca)))   # ORACLE label
            if with_chain:
                pr = I.project(Cs, cand.seq, cand.fold)
                for qq in ("ca", "phi", "psi"):
                    d[f"xchain_{qq}_{_safe(rs_key(s))}"] = np.asarray(pr[qq], float)
                d["xrmsdchain_" + _safe(rs_key(s))] = np.array(float(C2.rmsd_between(np.asarray(pr["ca"], float), cand.nat_ca)))
        d["x_R_RAND"] = np.array(int(R_RAND))
        _append_npz(M.cache_path(pdb), d)
        print(f"  [cache {q+1}/{len(pdbs)}] {pdb} n={cand.n} +{R_RAND-1} matched random-signed draws"
              + (" (+chain)" if with_chain else "")
              + f"  {time.time()-t1:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("cache:", CACHE)


# ============================================================================ per-target, extended
def meter_target(pdb, cost, basis="ca"):
    """S29's row, plus the R_RAND matched random-signed draws (E4)."""
    if basis == "chain-s28rows":
        row = M.meter_target_from_s28_rows(pdb, cost)
        row["rs_available"] = 1            # C2's stored rows carry draw 0 only
        return row
    row = M.meter_target(pdb, cost, basis=basis)
    ex = load_extra(pdb, need_chain=(basis == "chain"))
    if ex is None:
        row["rs_available"] = 1
        return row
    cand, dis, top, dg, frame, sur = M.load_target(pdb)
    ctx = M.Ctx(pdb, cand, frame, dg, top, dis)
    keys = [rs_key(s) for s in range(1, R_RAND)]
    if basis == "ca":
        W = np.stack([ex["S"][k] for k in keys])
        fv = cost.values(W, ctx, sur=sur)
        rr = [ex["oracle_rmsd"][k] for k in keys]
    else:
        W = np.stack([ex["chain"][k]["ca"] for k in keys])
        PHI = np.stack([ex["chain"][k]["phi"] for k in keys])
        PSI = np.stack([ex["chain"][k]["psi"] for k in keys])
        fv = cost.values(W, ctx, sur=sur, PHI=PHI, PSI=PSI)
        rr = [ex["oracle_rmsd_chain"][k] for k in keys]
    fp = row["f"]["PROD"]
    row["rs_f"] = {0: row["f"]["RAND_SIGNED[0]"]}
    row["rs_rmsd"] = {0: row["rmsd"]["RAND_SIGNED[0]"]}
    for s, k in enumerate(keys, start=1):
        row["rs_f"][s] = float(fv[s - 1]); row["rs_rmsd"][s] = float(rr[s - 1])
    row["rs_pref"] = {s: (1.0 if v < fp else (0.5 if v == fp else 0.0)) for s, v in row["rs_f"].items()}
    #: the CHARTER ladder recomputed on every draw: is its rho a property of the cost or of the draw?
    L = [k for k in M.LADDERS["CHARTER"] if k != "RAND_SIGNED[0]"]
    row["rs_rho_charter"] = {}
    for s in row["rs_f"]:
        fs = [row["rs_f"][s]] + [row["f"][k] for k in L]
        ds = [row["rs_rmsd"][s]] + [row["rmsd"][k] for k in L]
        row["rs_rho_charter"][s] = M.spearman(fs, ds)
    row["rs_available"] = len(row["rs_f"])
    return row


# ============================================================================ (E1) distributions
def quantiles(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return None
    p = np.percentile(x, [0, 10, 25, 50, 75, 90, 100])
    return dict(n=int(x.size), min=float(p[0]), p10=float(p[1]), q1=float(p[2]), median=float(p[3]),
                q3=float(p[4]), p90=float(p[5]), max=float(p[6]), mean=float(x.mean()),
                sd=float(x.std(ddof=1)) if x.size > 1 else float("nan"),
                iqr=float(p[4] - p[2]), frac_pos=float((x > 0).mean()), frac_neg=float((x < 0).mean()))


def extremes(x, pdbs, k=5):
    x = np.asarray(x, float)
    ok = np.isfinite(x)
    idx = np.where(ok)[0]
    if idx.size == 0:
        return None
    o = idx[np.argsort(x[idx])]
    return dict(lowest=[(str(pdbs[i]), float(x[i])) for i in o[:k]],
                highest=[(str(pdbs[i]), float(x[i])) for i in o[-k:][::-1]])


# ============================================================================ (E2) FAIL18 vs 108
def split_contrast(x, fail, folds, label, n_boot=2000):
    """mean(FAIL18) - mean(other 108), with a FOLD-CLUSTERED bootstrap CI.

    UNPAIRED: the two groups are different targets, so `ST.compare` (paired) does not apply.
    The cluster bootstrap resamples FOLDS with replacement and recomputes the difference; a
    resample in which either group is empty is dropped and counted.  Fold 0 holds no FAIL18
    target, which is a property of the pinned clustering, not of this code: with 5 folds and
    one empty cell this CI is wide by construction and is reported, not hidden.
    """
    x = np.asarray(x, float); fail = np.asarray(fail, bool); folds = np.asarray(folds)
    ok = np.isfinite(x)
    a, b = x[ok & fail], x[ok & ~fail]
    out = dict(label=str(label), n_fail=int(a.size), n_other=int(b.size))
    if a.size < 2 or b.size < 2:
        out.update(effect=float("nan"), se=float("nan"), mde=float("nan"), ci95_fold=None,
                   note="too few finite targets in one group")
        return out
    eff = float(a.mean() - b.mean())
    se = float(math.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size))        # Welch
    F = np.array(sorted(set(folds[ok].tolist())))
    rng = ST._rng(SALT, label, "splitfold")
    vals, dropped = [], 0
    for _ in range(n_boot):
        pick = rng.choice(F, len(F), replace=True)
        aa = np.concatenate([x[ok & fail & (folds == q)] for q in pick]) if len(pick) else np.array([])
        bb = np.concatenate([x[ok & ~fail & (folds == q)] for q in pick]) if len(pick) else np.array([])
        if aa.size < 1 or bb.size < 1:
            dropped += 1
            continue
        vals.append(aa.mean() - bb.mean())
    per_fold = {}
    for q in F:
        aq = x[ok & fail & (folds == q)]; bq = x[ok & ~fail & (folds == q)]
        per_fold[int(q)] = dict(n_fail=int(aq.size), n_other=int(bq.size),
                                diff=float(aq.mean() - bq.mean()) if aq.size and bq.size else None)
    out.update(mean_fail=float(a.mean()), mean_other=float(b.mean()), effect=eff, se=se,
               mde=float(ST.MDE_K * se), effect_over_mde=float(eff / (ST.MDE_K * se)) if se > 0 else float("nan"),
               median_fail=float(np.median(a)), median_other=float(np.median(b)),
               ci95_fold=[float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))] if len(vals) > 50 else None,
               n_boot_used=len(vals), n_boot_dropped=int(dropped), per_fold=per_fold,
               folds_with_fail=int(sum(1 for v in per_fold.values() if v["n_fail"] > 0)), n_folds=int(len(F)),
               folds_same_sign=int(sum(1 for v in per_fold.values() if v["diff"] is not None and np.sign(v["diff"]) == np.sign(eff))))
    return out


# ============================================================================ the extended summary
DIAGS = [("rho_S28", "ladder rho, S28 rungs"), ("rho_CHARTER", "ladder rho, CHARTER rungs"),
         ("rho_FULL", "ladder rho, FULL rungs"), ("cos", "gradient cosine"),
         ("native_pctile", "native percentile in pool"), ("pref_circ_best", "pref(circ_best vs PROD)")]


def _diag_vectors(rows):
    """The per-target vector of every diagnostic, by a single name, so E1/E2 cover all of them."""
    v = {}
    for nm in ("S28", "CHARTER", "FULL"):
        v["rho_" + nm] = np.array([r["rho"].get(nm, np.nan) for r in rows], float)
    if all("cos" in r for r in rows):
        v["cos"] = np.array([r["cos"] for r in rows], float)
    if all("native_pctile" in r for r in rows):
        v["native_pctile"] = np.array([r["native_pctile"] for r in rows], float)
    for k in ("circ_best", "NATIVE", "circ_opt", "RAND_SIGNED[0]"):
        if all(k in r.get("pref", {}) for r in rows):
            v["pref_" + k] = np.array([r["pref"][k] for r in rows], float)
    if all(k in r.get("pct", {}) for r in rows for k in ("PROD",)):
        v["pct_PROD"] = np.array([r["pct"]["PROD"] for r in rows], float)
    return v


def summarise(rows, cost, basis):
    base = M.summarise(rows, cost, basis)
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    fail = np.array([r["fail18"] for r in rows], bool)
    V = _diag_vectors(rows)
    n_cmp = 0

    # ---------------------------------------------------------------- (E1) per-target distributions
    base["per_target"] = {k: dict(quantiles=quantiles(v), extremes=extremes(v, pdbs)) for k, v in V.items()}
    base["per_target_note"] = ("a mean is not a distribution: read the IQR and the extremes before quoting the mean. "
                               "The five lowest/highest targets are named so a claim can be re-measured on them.")

    # ---------------------------------------------------------------- (E2) FAIL18 vs 108 contrast
    base["split_fail18"] = {k: split_contrast(v, fail, folds, f"{cost.name}/{basis}: {k} FAIL18 - other108")
                            for k, v in V.items()}
    n_cmp += len(base["split_fail18"])
    base["split_fail18_caveat"] = ("fold 0 contains no FAIL18 target (FAIL18 by pinned fold = 1:6, 2:2, 3:4, 4:6, 0:0), "
                                   "so the cluster bootstrap draws from 4 non-empty clusters; the CI is wide by "
                                   "construction and a split effect inside it is NOT a result.")

    # ---------------------------------------------------------------- (E4) matched random-signed control
    have_rs = all(r.get("rs_available", 1) > 1 for r in rows)
    if have_rs and all("pref" in r and "circ_best" in r["pref"] for r in rows):
        S = sorted(rows[0]["rs_f"].keys())
        pref_rs = np.array([[r["rs_pref"][s] for s in S] for r in rows], float)      # (n, R)
        ind_cb = np.array([r["pref"]["circ_best"] for r in rows], float)
        mean_over_draws = pref_rs.mean(1)
        c = ST.compare(ind_cb, mean_over_draws, folds, names=pdbs, seed_parts=(SALT,),
                       label=f"{cost.name}/{basis}: pref(circ_best vs PROD) - pref(matched RAND_SIGNED vs PROD), {len(S)} draws")
        per_draw = []
        for q, s in enumerate(S):
            cq = ST.compare(ind_cb, pref_rs[:, q], folds, seed_parts=(SALT,), label=f"draw {s}")
            per_draw.append(dict(draw=int(s), pref_rs=float(pref_rs[:, q].mean()), effect=cq["effect"],
                                 se=cq["se"], effect_over_mde=cq["effect_over_mde"],
                                 ci95_fold=cq["ci95_fold"], folds_same_sign=cq["folds_same_sign"]))
        rho_rs = np.array([[r["rs_rho_charter"][s] for s in S] for r in rows], float)
        rmsd_rs = np.array([[r["rs_rmsd"][s] for s in S] for r in rows], float)
        base["rand_signed_control"] = dict(
            n_draws=len(S), pref_circ_best=float(ind_cb.mean()),
            pref_rand_signed_mean=float(pref_rs.mean()),
            pref_rand_signed_per_draw=[float(x) for x in pref_rs.mean(0)],
            pref_rand_signed_draw_sd=float(pref_rs.mean(0).std(ddof=1)),
            contrast={k: v for k, v in c.items() if k != "concentration"},
            concentration=c.get("concentration"), fmt=ST.fmt(c), per_draw=per_draw,
            single_draw_range=[float(min(d["effect"] for d in per_draw)), float(max(d["effect"] for d in per_draw))],
            charter_rho_per_draw=[float(np.nanmean(rho_rs[:, q])) for q in range(len(S))],
            charter_rho_draw_sd=float(np.nanmean(rho_rs, 0).std(ddof=1)),
            rmsd_rand_signed_mean=float(rmsd_rs.mean()),
            note="the matched random-signed structure is a NATIVE-FREE direction at the ORACLE displacement of "
                 "circ_best from production: it is the control that asks whether f prefers the 0.29 A oracle "
                 "structure BECAUSE it is near the native, or merely because it is a displacement of production. "
                 "S29 used one draw; the per-draw spread below prices what that single draw was worth.")
        n_cmp += 1 + len(per_draw)
    else:
        base["rand_signed_control"] = None

    # ---------------------------------------------------------------- (E5) aggregate cosine null
    if all("cos_random_ref" in r for r in rows):
        Rf = np.array([r["cos_random_ref"] for r in rows], float)                    # (n, 16)
        tm = Rf.mean(0)                                                               # the 126-target MEAN per draw
        cs = V["cos"]; ok = np.isfinite(cs)
        obs = float(cs[ok].mean())
        base["cosine_null"] = dict(
            n_ref_draws=int(Rf.shape[1]),
            per_draw_abs_mean=float(np.abs(Rf).mean()),
            target_mean_null=dict(mean=float(tm.mean()), sd=float(tm.std(ddof=1)),
                                  p2_5=float(np.percentile(tm, 2.5)), p50=float(np.percentile(tm, 50)),
                                  p97_5=float(np.percentile(tm, 97.5)),
                                  abs_p95=float(np.percentile(np.abs(tm), 95))),
            observed_target_mean=obs,
            observed_z_vs_null=float((obs - tm.mean()) / tm.std(ddof=1)) if tm.std(ddof=1) > 0 else float("nan"),
            note="TWO different nulls, do not mix them. `per_draw_abs_mean` (0.140 for the shipped cost) is the "
                 "magnitude of ONE random direction on ONE target -- it is NOT the null for the 126-target mean. "
                 "`target_mean_null` is that null, and it is ~10x tighter; a cosine of +0.05 is far below 0.140 "
                 "and still several sigma above chance.")
        n_cmp += 1
    else:
        base["cosine_null"] = None

    # ---------------------------------------------------------------- (E6) verdict, (E7) multiplicity
    base["verdict"] = verdict(base, cost, basis)
    n_cmp += len(base["ladder_rho"]) + (1 if base.get("cosine") else 0) + (1 if base.get("native_pctile") else 0) + len(base["pref"])
    base["multiplicity"] = dict(
        comparisons_emitted=int(n_cmp),
        note="every comparison this meter emits, counted. With 8 lanes, a 1x-MDE positive is EXPECTED somewhere; "
             "price any single positive from this instrument against this count and the sprint's running total.")
    base["anchors"] = ANCHORS
    base["text"] = base["text"] + "\n" + render_ext(base)
    return base


# ============================================================================ (E6) the gates
def verdict(o, cost, basis):
    """Pre-stated gates.  The meter never tunes: these thresholds are fixed here, in code, and a
    cost is read against the SHIPPED cost's metered values (the incumbent), not against a
    threshold chosen after seeing the cost's numbers."""
    g, reasons = {}, []
    key = "ladder_rho_S28_chain" if basis.startswith("chain") else "ladder_rho_S28_ca"
    s = o["ladder_rho"]["S28"]
    ci = s.get("ci95_fold")
    if ci is None:
        g["ladder"] = "INDETERMINATE"; reasons.append("ladder rho CI undefined (<3 finite targets)")
    elif ci[1] < 0:
        g["ladder"] = "BLOCK"
        reasons.append(f"S28-ladder rho {s['mean']:+.4f}, fold CI [{ci[0]:+.3f},{ci[1]:+.3f}] ENTIRELY BELOW ZERO: "
                       "lowering this cost RAISES RMSD on the ladder. The charter forbids giving this a day of compute.")
    elif ci[0] > 0:
        g["ladder"] = "PASS"; reasons.append(f"S28-ladder rho {s['mean']:+.4f}, fold CI excludes zero on the good side "
                                             f"(incumbent {ANCHORS[key]:+.4f})")
    else:
        g["ladder"] = "NEUTRAL"; reasons.append(f"S28-ladder rho {s['mean']:+.4f}, fold CI spans zero "
                                                f"(incumbent {ANCHORS[key]:+.4f})")
    c = o.get("cosine"); nl = o.get("cosine_null")
    if c and c.get("ci95_fold") and nl:
        z = nl["observed_z_vs_null"]
        if c["n_nan"]:
            g["cosine"] = "INDETERMINATE"; reasons.append(f"cosine undefined on {c['n_nan']} targets: the survivors are a selected subsample")
        elif c["mean"] > 0 and abs(z) >= 2:
            g["cosine"] = "PASS"; reasons.append(f"cosine {c['mean']:+.4f}, {z:+.1f} sigma above the target-mean random null (incumbent {ANCHORS['cosine']:+.4f})")
        elif c["mean"] < 0 and abs(z) >= 2:
            g["cosine"] = "FAIL"; reasons.append(f"cosine {c['mean']:+.4f}, {z:+.1f} sigma on the WRONG side of the random null")
        else:
            g["cosine"] = "NEUTRAL"; reasons.append(f"cosine {c['mean']:+.4f} within 2 sigma of the random null")
        sh = c.get("shrink_signature")
        if sh and c["mean"] > 0 and sh["rg_ratio_mean"] < 1.0:
            g["cosine"] = "SUSPECT"
            reasons.append(f"the positive cosine comes with a CONTRACTION (Rg ratio {sh['rg_ratio_mean']:.4f} < 1): "
                           "contract addendum 20 -- a cosine bought by shrinking is not information.")
    else:
        g["cosine"] = "N/A"
    p = o.get("native_pctile")
    if p and p.get("ci95_fold"):
        inc = ANCHORS["native_pctile_DIS"]
        g["pctile"] = "PASS" if p["ci95_fold"][1] < inc else ("FAIL" if p["ci95_fold"][0] > inc else "NEUTRAL")
        reasons.append(f"native percentile {p['mean']:.4f} [{p['ci95_fold'][0]:.3f},{p['ci95_fold'][1]:.3f}] vs incumbent {inc:.4f} "
                       "(0 = the native is the pool's best under f; 0.5 = chance)")
    else:
        g["pctile"] = "N/A"
    rc = o.get("rand_signed_control")
    if rc:
        e = rc["contrast"]; r = e["effect_over_mde"]
        lo, hi = rc["single_draw_range"]
        if abs(r) < 0.7:
            g["rand_signed"] = "BELOW 0.7x MDE -- NOT A RESULT"
        elif abs(r) < 1.0:
            g["rand_signed"] = "0.7-1.0x MDE -- NOT DEMONSTRATED"
        else:
            g["rand_signed"] = "PASS" if e["effect"] > 0 else "FAIL"
        reasons.append(f"pref(circ_best) - pref(matched random-signed, {rc['n_draws']} draws) {e['effect']:+.4f} "
                       f"({r:+.2f}x MDE, folds same sign {e['folds_same_sign']}/{e['n_folds']}); a SINGLE draw would "
                       f"have given anything in [{lo:+.4f}, {hi:+.4f}]")
    else:
        g["rand_signed"] = "N/A (no extra draws cached; run build-cache)"
    g["overall"] = "BLOCK" if "BLOCK" in g.values() else ("REVIEW" if any(v in ("SUSPECT", "FAIL") for v in g.values()) else "PROCEED")
    return dict(gates=g, reasons=reasons,
                rule="below 0.7x MDE is not a result; 0.7-1.0x is not a demonstrated improvement; "
                     "a BLOCK on the ladder means the cost moves in the wrong structural direction and gets no endpoint compute.")


# ============================================================================ rendering
def render_ext(o):
    L = ["  --- S30 EXTENSIONS " + "-" * 60]
    L.append("  (E1) PER-TARGET DISTRIBUTIONS (a mean is not a distribution):")
    for k, d in o["per_target"].items():
        q = d["quantiles"]
        if not q:
            continue
        L.append(f"      {k:18s} n {q['n']:3d}  min {q['min']:+.3f}  p10 {q['p10']:+.3f}  q1 {q['q1']:+.3f}  med {q['median']:+.3f}  "
                 f"q3 {q['q3']:+.3f}  p90 {q['p90']:+.3f}  max {q['max']:+.3f}  | mean {q['mean']:+.4f} sd {q['sd']:.3f} IQR {q['iqr']:.3f}  "
                 f"pos {100*q['frac_pos']:.0f}%")
        e = d["extremes"]
        if e:
            L.append(f"        lowest  " + " ".join(f"{p}:{v:+.2f}" for p, v in e["lowest"])
                     + "   highest " + " ".join(f"{p}:{v:+.2f}" for p, v in e["highest"]))
    L.append("  (E2) FAIL18 (18) vs OTHER (108), fold-clustered contrast -- NOT two bare means:")
    for k, s in o["split_fail18"].items():
        if s.get("ci95_fold") is None:
            L.append(f"      {k:18s} {s.get('effect', float('nan')):+.4f}  CI undefined ({s.get('note','')})")
            continue
        L.append(f"      {k:18s} FAIL18 {s['mean_fail']:+.4f} (med {s['median_fail']:+.3f})  other {s['mean_other']:+.4f} (med {s['median_other']:+.3f})  "
                 f"diff {s['effect']:+.4f}  SE {s['se']:.4f}  {s['effect_over_mde']:+.2f}x MDE  fold CI [{s['ci95_fold'][0]:+.3f}, {s['ci95_fold'][1]:+.3f}]  "
                 f"folds same sign {s['folds_same_sign']}/{s['folds_with_fail']}")
    L.append("      CAVEAT: " + o["split_fail18_caveat"])
    rc = o.get("rand_signed_control")
    L.append("  (E4) MATCHED RANDOM-SIGNED CONTROL (first class, %s draws):" % (rc["n_draws"] if rc else "1 -- NOT CACHED"))
    if rc:
        e = rc["contrast"]
        L.append(f"      pref(circ_best vs PROD) {rc['pref_circ_best']:.4f}   pref(matched random-signed vs PROD) {rc['pref_rand_signed_mean']:.4f} "
                 f"(per draw {['%.3f' % x for x in rc['pref_rand_signed_per_draw']]}, draw sd {rc['pref_rand_signed_draw_sd']:.4f})")
        L.append(f"      contrast {e['effect']:+.4f}  SE {e['se']:.4f}  MDE {e['mde']:.4f}  {e['effect_over_mde']:+.2f}x  "
                 f"fold CI [{e['ci95_fold'][0]:+.3f}, {e['ci95_fold'][1]:+.3f}]  folds same sign {e['folds_same_sign']}/{e['n_folds']}  "
                 f"median {e['median_effect']:+.4f}")
        L.append(f"      ONE DRAW would have given anything in [{rc['single_draw_range'][0]:+.4f}, {rc['single_draw_range'][1]:+.4f}] "
                 f"-- S29 reported the seed-0 draw alone")
        L.append(f"      CHARTER ladder rho per draw: {['%+.3f' % x for x in rc['charter_rho_per_draw']]} (draw sd {rc['charter_rho_draw_sd']:.4f}); "
                 f"mean ORACLE RMSD of the random-signed rung {rc['rmsd_rand_signed_mean']:.3f} A")
        cc = rc.get("concentration")
        if cc:
            L.append(f"      concentration vs the UNIFORM-EFFECT null: mean {cc['mean']:+.4f} median {cc['median']:+.4f} "
                     f"drop-top10 {cc['drop_top10_mean']:+.4f} sits at the {100*cc['pctile_in_null']:.0f}th pct of the null "
                     f"[p10 {cc['null_p10']:+.3f}, p90 {cc['null_p90']:+.3f}]  flag {cc['flag']}")
    else:
        L.append("      NOT AVAILABLE: no extra draws in the cache for this basis. Run `build-cache`"
                 " (and `--chain` for the chain basis). S29's single seed-0 draw is all this run has.")
    nl = o.get("cosine_null")
    if nl:
        t = nl["target_mean_null"]
        L.append("  (E5) COSINE NULLS, the two that must not be mixed:")
        L.append(f"      one random direction, one target:   mean |cos| {nl['per_draw_abs_mean']:.4f}   <- NOT the null for a 126-target mean")
        L.append(f"      the 126-target MEAN under random:   {t['mean']:+.4f} +/- {t['sd']:.4f}  [p2.5 {t['p2_5']:+.4f}, p97.5 {t['p97_5']:+.4f}]  |mean| p95 {t['abs_p95']:.4f}")
        L.append(f"      observed target mean {nl['observed_target_mean']:+.4f}  =  {nl['observed_z_vs_null']:+.2f} sigma of THAT null")
    v = o["verdict"]
    L.append("  (E6) VERDICT: " + "  ".join(f"{k}={val}" for k, val in v["gates"].items()))
    for r in v["reasons"]:
        L.append("        - " + r)
    L.append("        rule: " + v["rule"])
    L.append(f"  (E7) MULTIPLICITY: {o['multiplicity']['comparisons_emitted']} comparisons emitted by this run. "
             + o["multiplicity"]["note"])
    return "\n".join(L)


# ============================================================================ driver
def run_meter(spec, basis="ca", limit=0, fd_h=M.FD_H_DEFAULT, tag=None, out=None, pdbs=None, save=True, quiet=False):
    from s25 import phys_lib as P
    cost = M.Cost(spec, fd_h=fd_h)
    if pdbs is None:
        pdbs = P.targets()[:limit] if limit else P.targets()
    rows, t0 = [], time.time()
    for q, pdb in enumerate(pdbs):
        row = meter_target(pdb, cost, basis=basis)
        if row.get("nan_on_rungs") and cost.kind == "callable":
            raise ValueError(f"{pdb}: the cost gave NaN on {row['nan_on_rungs']} under the POISONED context "
                             "(ctx.nat_ca / ctx.oracle_rr are NaN): it reads a native quantity or is undefined; "
                             "the meter refuses it. Errors: " + "; ".join(cost.errors[:4]))
        rows.append(row)
        if not quiet and ((q + 1) % 20 == 0 or q + 1 == len(pdbs)):
            print(f"  [{q+1}/{len(pdbs)}] {pdb} rho S28 {row['rho']['S28']:+.3f}"
                  + (f" cos {row['cos']:+.3f}" if "cos" in row else "")
                  + f"  ({(time.time()-t0)/60:.1f} min)", flush=True)
    summ = summarise(rows, cost, basis)
    summ["rows"] = rows
    summ["fd_h"] = fd_h
    summ["meter"] = "s30/s30_D_meter.py (extends s29/s29_D_cost_audit.py)"
    if save:
        os.makedirs(RESULTS, exist_ok=True)
        tag = tag or cost.name.replace("/", "_")
        path = out or os.path.join(RESULTS, f"s30_D_meter_{tag}_{basis}.json")
        ST.save_atomic(path, summ, rows=rows, n_expected=len(pdbs), complete_keys=("pdb", "f", "rho", "pref"),
                       module_file=__file__)
        summ["path"] = path
    if not quiet:
        print(summ["text"])
        if save:
            print("wrote", summ["path"])
    return summ


# ============================================================================ verify / selftest
def verify():
    """Re-run the shipped cost through the S29 instrument and ASSERT the six published anchors."""
    import tempfile
    bad = []
    with tempfile.TemporaryDirectory() as td:
        a = M.run_meter("DIS", basis="chain-s28rows", out=os.path.join(td, "a.json"), quiet=True)
        b = M.run_meter("DIS", basis="ca", out=os.path.join(td, "b.json"), quiet=True)
        s = M.run_meter("DIS_SURR", basis="ca", out=os.path.join(td, "c.json"), quiet=True)
    checks = [
        ("ladder rho S28, chain", a["ladder_rho"]["S28"]["mean"], ANCHORS["ladder_rho_S28_chain"], 5e-4),
        ("ladder rho S28, CA", b["ladder_rho"]["S28"]["mean"], ANCHORS["ladder_rho_S28_ca"], 5e-4),
        ("gradient cosine", b["cosine"]["mean"], ANCHORS["cosine"], 5e-4),
        ("random-direction |cos|", b["cosine"]["random_ref_mean_abs"], ANCHORS["cosine_random_ref_abs"], 1e-3),
        ("native pctile, DIS", b["native_pctile"]["mean"], ANCHORS["native_pctile_DIS"], 5e-4),
        ("native pctile, DIS_SURR", s["native_pctile"]["mean"], ANCHORS["native_pctile_DIS_SURR"], 5e-4),
        ("pref circ_best, chain", a["pref"]["circ_best"]["mean"], ANCHORS["pref_circ_best_chain"], 5e-4),
        ("pref circ_best, CA", b["pref"]["circ_best"]["mean"], ANCHORS["pref_circ_best_ca"], 5e-4),
        ("pool-member control, chain", a["pool_member_control"]["pref_pool_member_vs_prod"], ANCHORS["pool_member_chain"], 5e-4),
        ("pool-member control, CA", b["pool_member_control"]["pref_pool_member_vs_prod"], ANCHORS["pool_member_ca"], 5e-4),
    ]
    for name, got, want, tol in checks:
        ok = abs(got - want) <= tol
        print(f"  {'OK ' if ok else 'FAIL'} {name:30s} got {got:+.4f}  want {want:+.4f}  dev {abs(got-want):.2e}")
        if not ok:
            bad.append(name)
    if bad:
        raise AssertionError("anchors that did NOT reproduce: " + ", ".join(bad))
    print("  ALL %d S29 ANCHORS REPRODUCE through s30/s30_D_meter.py -> s29/s29_D_cost_audit.py" % len(checks))


def selftest():
    M.selftest()
    q = quantiles([1.0, 2.0, 3.0, 4.0, np.nan])
    assert q["n"] == 4 and q["median"] == 2.5 and q["frac_pos"] == 1.0, q
    assert quantiles([np.nan, np.nan]) is None
    e = extremes([3.0, 1.0, 2.0], ["a", "b", "c"], k=2)
    assert e["lowest"][0][0] == "b" and e["highest"][0][0] == "a", e
    # the split contrast on a known effect
    rng = np.random.default_rng(0)
    x = rng.normal(size=126); fail = np.zeros(126, bool); fail[:18] = True
    x[fail] += 1.0
    folds = np.repeat(np.arange(5), 26)[:126]
    s = split_contrast(x, fail, folds, "selftest")
    assert 0.5 < s["effect"] < 1.5, s["effect"]
    assert s["ci95_fold"] is not None and s["ci95_fold"][0] > 0, s["ci95_fold"]
    # a group with no members anywhere -> graceful
    s2 = split_contrast(x, np.zeros(126, bool), folds, "selftest empty")
    assert not np.isfinite(s2["effect"]), s2
    # the tie convention: an exact tie is half, never an array-order win (S12's 1.386 A lesson)
    assert M.pct_in_pool(np.array([1.0, 1.0, 2.0]), 1.0) == 1.0 / 3.0
    print("  s30_D_meter selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["build-cache", "meter", "verify", "selftest"])
    ap.add_argument("--f", default="DIS")
    ap.add_argument("--basis", default="ca", choices=["ca", "chain", "chain-s28rows"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--slice", default=None, help="i/n: build only shard i of n (0-based)")
    ap.add_argument("--fd-h", type=float, default=M.FD_H_DEFAULT)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--chain", action="store_true")
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()
    if a.mode == "selftest":
        return selftest()
    if a.mode == "verify":
        return verify()
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    if a.slice:
        i, n = (int(x) for x in a.slice.split("/"))
        pdbs = [p for q, p in enumerate(pdbs) if q % n == i]
    if a.mode == "build-cache":
        os.makedirs(CACHE, exist_ok=True)
        build_cache(pdbs, with_chain=a.chain, rebuild=a.rebuild)
    else:
        run_meter(a.f, basis=a.basis, limit=0, fd_h=a.fd_h, tag=a.tag, out=a.out, pdbs=pdbs)


if __name__ == "__main__":
    main()
