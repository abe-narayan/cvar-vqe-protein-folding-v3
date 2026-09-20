#!/usr/bin/env python
"""s29/s29_D_band.py -- lane D: WITHIN-REALISM-BAND ORDERING (pre-registered in `s29/PREREG_S29_D_band.md`).

EVERY NUMBER HERE IS ORACLE (ordering skill is measured against the CA-RMSD to the native).
Nothing is deployable; nothing is tuned; no native quantity chooses any parameter.

The question lane L's S29-L12 leaves open: the perception-distortion theorem forbids a realism
measure from preferring the distortion-optimal answer ACROSS realism levels, and says nothing
about ordering INSIDE one level. So: condition on a native-free realism statistic R, then ask
whether a scorer orders structures by accuracy inside a band of R.

Per target: the 500 pool members (real traces, the homogeneous POOL-ONLY arm) and, reported
beside it, the same with the meter's nine ladder rungs added. Three pre-registered realism
statistics (R1 CAGEO pool-percentile, R2 the geometric bond/Rg band, R3 the leave-one-out
consensus percentile), five fixed bands, the fixed occupancy rule (>= 20 structures and >= 10
distinct RMSDs), the mandatory DEGENERACY check (|Spearman(scorer, R)| > 0.8 inside the band, or
< 20% of the scorer's full range retained), the within-cell shuffle null, the across-band
comparison, and the max-over-scorers sign-flip null.

    python s29/s29_D_band.py probe --limit 1
    python s29/s29_D_band.py run [--limit N] [--basis ca]
    python s29/s29_D_band.py analyse [--basis ca]

Rows `s29/results/s29_D_band_<basis>_rows.jsonl` (one per target, resumable); summary
`s29/results/s29_D_band_<basis>.json`.
"""
from __future__ import annotations

import argparse
import json
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

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_C2_recog_audit as C2   # noqa: E402
from s29 import s29_D_cost_audit as M      # noqa: E402

RESULTS = os.path.join(HERE, "results")
SALT = "s29D"
BAND_EDGES = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)          # fixed in the prereg, never re-cut
MIN_CELL, MIN_DISTINCT_RMSD = 20, 10                 # the occupancy rule, fixed in the prereg
DEGEN_RHO, DEGEN_RANGE = 0.8, 0.20                   # the degeneracy rule, fixed in the prereg
N_SHUFFLE = 200
REALISMS = ("R1_CAGEO", "R2_GEOM", "R3_CONS")
#: ADDENDUM 1 (L's point 2): the width curve is the headline object.  Widths of the R-percentile
#: axis and the fixed window centres, both registered before any analysis.
WIDTHS = (1.0, 0.5, 0.25, 0.10, 0.05)
CENTRES = (0.1, 0.3, 0.5, 0.7, 0.9)
#: the CA-level scorers the S27 cache carries per pool member (so no recomputation on 500 members)
POOL_SCORERS = ["DIS", "DIS_MEAN", "CONTACT_LL", "DISTPOT", "CONTACT", "ENV", "HP", "RG_LAW",
                "RG_UNIV", "EXVOL", "CAGEO", "SS_MATCH"]


def spearman(x, y):
    return M.spearman(x, y)


def _ranks(x):
    from scipy.stats import rankdata
    return rankdata(np.asarray(x, float))


def _pearson_ranks(rx, ry):
    """Spearman as the Pearson correlation of mid-ranks: identical to `scipy.stats.spearmanr`
    including its tie correction, and 200x cheaper in a permutation loop."""
    a = rx - rx.mean(); b = ry - ry.mean()
    d = np.sqrt((a @ a) * (b @ b))
    return float(a @ b / d) if d > 1e-300 else float("nan")


def _shuffle_null(rx, ry, rng, n_draw=N_SHUFFLE):
    """The registered within-cell shuffle null, vectorised: permute the RMSD ranks inside the
    cell and recompute the same statistic (mid-ranks are permutation invariant, so permuting the
    ranks is permuting the values)."""
    a = rx - rx.mean()
    sa = np.sqrt(a @ a)
    b0 = ry - ry.mean()
    sb = np.sqrt(b0 @ b0)
    if sa < 1e-300 or sb < 1e-300:
        return np.full(n_draw, np.nan)
    Bm = np.array([rng.permutation(b0) for _ in range(n_draw)])
    return (Bm @ a) / (sa * sb)


def pct_of(values, ref):
    """Percentile of each value among the reference population (share below + half the ties)."""
    v = np.asarray(values, float)[:, None]
    r = np.asarray(ref, float)[None, :]
    return ((r < v).mean(1) + 0.5 * (r == v).mean(1))


# ============================================================ per target
def target_row(pdb, with_rungs=True):
    t0 = time.time()
    cand, dis, top, dg, frame, sur = M.load_target(pdb)
    n, k = int(cand.n), int(cand.k)
    W = np.asarray(frame.Wp, float)                                   # the pool, posed
    with np.load(os.path.join(M.S27_CACHE, f"{pdb}.npz")) as z:
        ch = {key: np.asarray(z[key], float) for key in z.files if key != "cost_ms"}
    scorers = [s for s in POOL_SCORERS if s in ch]
    rr = np.asarray(cand.oracle_rr, float)                            # ORACLE: per-member RMSD
    # geometry of every member (native-free)
    bond = np.linalg.norm(np.diff(W, axis=1), axis=2).mean(1)
    rg = np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    # leave-one-out consensus (mean CA-RMSD to the other members), native-free
    P = I.pairwise_rmsd(W)
    cons = (P.sum(1) / max(k - 1, 1))
    # the ladder rungs, scored through the meter's own CA adapters
    rung_names, rung_scores, rung_rmsd, rung_geom, rung_cons = [], {}, [], [], []
    if with_rungs:
        rungs = M.load_rungs(pdb)
        if rungs is not None:
            rung_names = list(M.RUNGS)
            Wr = np.stack([rungs["S"][x] for x in rung_names])
            ctx = M.Ctx(pdb, cand, frame, dg, top, dis)
            sc = C2.ca_scores(ctx.cand_poisoned(), Wr, ctx.universe(), dg)
            i2, j2 = I.pair_index(n)
            Dr = I.pair_dists(Wr, i2, j2)
            sc["DIS"] = np.asarray(I.shipped_score(dg, Dr.astype(np.float32).astype(float)), float)
            rung_scores = {s: np.asarray(sc[s], float) for s in scorers if s in sc}
            rung_rmsd = np.array([rungs["oracle_rmsd"][x] for x in rung_names])          # ORACLE
            rung_geom = (np.linalg.norm(np.diff(Wr, axis=1), axis=2).mean(1),
                         np.sqrt(((Wr - Wr.mean(1, keepdims=True)) ** 2).sum(2).mean(1)))
            rung_cons = np.array([float(I.kabsch_rmsd_batch(W, x).mean()) for x in Wr])
    row = dict(pdb=pdb, n=n, k=k, fold=int(cand.fold), fail18=bool(pdb in I.FAIL18),
               scorers=scorers, rungs=rung_names, basis="ca", oracle=True, cells={})
    #: ADDENDUM 2 (the coordinator, from lane L's S29-L19): if the realism statistic is
    #: compactness-loaded, the band removes the very axis S14 measured as carrying in-band skill
    #: (+0.909 with the native's z-scored Rg), and a null is uninterpretable.  rho(R, Rg) on the
    #: pool is native-free and is computed for every R; the ORACLE reference rho(Rg, RMSD) is
    #: computed beside it (labelled ORACLE) because it is the axis in question.
    row["compactness"] = dict(rho_Rg_rmsd_ORACLE=_pearson_ranks(_ranks(rg), _ranks(rr)),
                              rho_bond_rmsd_ORACLE=_pearson_ranks(_ranks(bond), _ranks(rr)),
                              sd_rg=float(rg.std()), mean_rg=float(rg.mean()))
    # ---- the three realism statistics, on the pool and on the rungs
    R = {}
    R["R1_CAGEO"] = (pct_of(ch["CAGEO"], ch["CAGEO"]),
                     pct_of(rung_scores["CAGEO"], ch["CAGEO"]) if rung_names and "CAGEO" in rung_scores else None)
    pb, pr = pct_of(bond, bond), pct_of(rg, rg)
    med = np.array([0.5, 0.5])
    R["R2_GEOM"] = (np.sqrt((pb - med[0]) ** 2 + (pr - med[1]) ** 2),
                    (np.sqrt((pct_of(rung_geom[0], bond) - 0.5) ** 2 + (pct_of(rung_geom[1], rg) - 0.5) ** 2)
                     if rung_names else None))
    R["R3_CONS"] = (pct_of(cons, cons), pct_of(rung_cons, cons) if rung_names else None)
    #: R2 is a distance-to-median, not a percentile: re-express it as its own percentile so the
    #: five fixed band edges mean the same thing for all three statistics (a monotone re-mapping,
    #: registered as "mapped to its percentile" in the prereg's R2 definition).
    if R["R2_GEOM"][1] is not None:
        R["R2_GEOM"] = (pct_of(R["R2_GEOM"][0], R["R2_GEOM"][0]), pct_of(R["R2_GEOM"][1], R["R2_GEOM"][0]))
    else:
        R["R2_GEOM"] = (pct_of(R["R2_GEOM"][0], R["R2_GEOM"][0]), None)
    for rname in REALISMS:
        rp = None
        if rname == "R1_CAGEO":
            rp = pct_of(ch["CAGEO"], ch["CAGEO"])
        elif rname == "R2_GEOM":
            _pb, _pr = pct_of(bond, bond), pct_of(rg, rg)
            rp = pct_of(np.sqrt((_pb - 0.5) ** 2 + (_pr - 0.5) ** 2), np.sqrt((_pb - 0.5) ** 2 + (_pr - 0.5) ** 2))
        else:
            rp = pct_of(cons, cons)
        row["compactness"][rname] = dict(rho_R_rg=_pearson_ranks(_ranks(rp), _ranks(rg)),
                                         rho_R_bond=_pearson_ranks(_ranks(rp), _ranks(bond)),
                                         rho_R_rmsd_ORACLE=_pearson_ranks(_ranks(rp), _ranks(rr)))
    for arm, use_rungs in (("pool", False), ("pool+rungs", True)):
        if use_rungs and not rung_names:
            continue
        for rname in REALISMS:
            rp, rr_rung = R[rname]
            if use_rungs and rr_rung is None:
                continue
            rv = np.concatenate([rp, rr_rung]) if use_rungs else rp
            y = np.concatenate([rr, rung_rmsd]) if use_rungs else rr                      # ORACLE
            for s in scorers:
                sv = np.concatenate([ch[s], rung_scores[s]]) if (use_rungs and s in rung_scores) else ch[s]
                if use_rungs and s not in rung_scores:
                    continue
                full_rng = float(np.ptp(sv))
                # across-band (the record's number)
                rho_ac = _pearson_ranks(_ranks(sv), _ranks(y))
                cells = []
                for b in range(len(BAND_EDGES) - 1):
                    lo, hi = BAND_EDGES[b], BAND_EDGES[b + 1]
                    m = (rv >= lo) & (rv < hi) if b < len(BAND_EDGES) - 2 else (rv >= lo) & (rv <= hi)
                    nb = int(m.sum())
                    if nb < MIN_CELL or len(np.unique(y[m])) < MIN_DISTINCT_RMSD:
                        cells.append(dict(band=b, n=nb, rho=None, reason="thin"))
                        continue
                    rx, ry, rR = _ranks(sv[m]), _ranks(y[m]), _ranks(rv[m])
                    rho = _pearson_ranks(rx, ry)
                    rho_sR = _pearson_ranks(rx, rR)
                    frac_rng = float(np.ptp(sv[m]) / max(full_rng, 1e-30))
                    degen = bool(abs(rho_sR) > DEGEN_RHO if np.isfinite(rho_sR) else False) or frac_rng < DEGEN_RANGE
                    rng = SD.stable_rng(pdb, "s29D_band_null", s, rname, b, salt=SALT)
                    null = _shuffle_null(rx, ry, rng)
                    fin = np.isfinite(null)
                    cells.append(dict(band=b, n=nb, rho=rho, rho_scorer_vs_R=rho_sR, frac_range=frac_rng,
                                      degenerate=degen,
                                      null_mean=float(null[fin].mean()) if fin.any() else float("nan"),
                                      null_p95=float(np.percentile(null[fin], 95)) if fin.any() else float("nan"),
                                      null_p05=float(np.percentile(null[fin], 5)) if fin.any() else float("nan")))
                #: ADDENDUM 1 (L's point 1): the closed form.  rho_SR is NATIVE-FREE; rho_SY and
                #: rho_RY are ORACLE (they read the RMSD).  On ranks, so the Gaussian assumption
                #: enters only through the copula.
                rS, rY, rR = _ranks(sv), _ranks(y), _ranks(rv)
                r_SY = _pearson_ranks(rS, rY); r_SR = _pearson_ranks(rS, rR); r_RY = _pearson_ranks(rR, rY)
                den = np.sqrt(max((1 - r_SR ** 2) * (1 - r_RY ** 2), 1e-30))
                partial = float((r_SY - r_SR * r_RY) / den) if den > 1e-15 else float("nan")
                #: ADDENDUM 1 (L's point 2): the width curve, sliding windows on the R-percentile axis.
                curve = {}
                for w in WIDTHS:
                    vals, ns = [], []
                    for c0 in (CENTRES if w < 1.0 else (0.5,)):
                        lo, hi = c0 - w / 2.0, c0 + w / 2.0
                        mw = (rv >= lo) & (rv <= hi)
                        nb = int(mw.sum())
                        if nb < MIN_CELL or len(np.unique(y[mw])) < MIN_DISTINCT_RMSD:
                            continue
                        vals.append(_pearson_ranks(_ranks(sv[mw]), _ranks(y[mw]))); ns.append(nb)
                    curve[str(w)] = dict(rho=float(np.mean(vals)) if vals else None,
                                         n_windows=len(vals), mean_n=float(np.mean(ns)) if ns else 0.0)
                row["cells"].setdefault(arm, {}).setdefault(rname, {})[s] = dict(
                    rho_across=rho_ac, bands=cells, rho_SY=r_SY, rho_SR=r_SR, rho_RY=r_RY,
                    partial_pred=partial, width_curve=curve)
    row["secs"] = time.time() - t0
    return row


# ============================================================ driver
def rows_path(basis="ca"):
    return os.path.join(RESULTS, f"s29_D_band_{basis}_rows.jsonl")


def run(pdbs, basis="ca"):
    os.makedirs(RESULTS, exist_ok=True)
    path = rows_path(basis)
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    done.add(json.loads(line)["pdb"])
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        row = target_row(pdb)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(row) + "\n")
        c = row["cells"].get("pool", {}).get("R1_CAGEO", {}).get("DIS", {})
        ok = [b["rho"] for b in c.get("bands", []) if b.get("rho") is not None]
        print(f"  [{q+1}/{len(pdbs)}] {pdb} n={row['n']} DIS|R1 across {c.get('rho_across', float('nan')):+.3f} "
              f"in-band {' '.join('%+.3f' % v for v in ok) if ok else '(all thin)'}  {row['secs']:.1f}s "
              f"(elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", path)


# ============================================================ analysis
def analyse(basis="ca", n_max_null=2000):
    rows = [json.loads(l) for l in open(rows_path(basis), encoding="utf-8") if l.strip()]
    rows.sort(key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    fail = np.array([r["fail18"] for r in rows])
    out = dict(check="within-realism-band ordering (PREREG_S29_D_band)", n=len(rows), basis=basis,
               oracle=True, band_edges=list(BAND_EDGES), min_cell=MIN_CELL, arms={})
    for arm in ("pool", "pool+rungs"):
        if arm not in rows[0]["cells"]:
            continue
        out["arms"][arm] = {}
        for rname in REALISMS:
            if rname not in rows[0]["cells"][arm]:
                continue
            scorers = sorted(rows[0]["cells"][arm][rname])
            res = {}
            pooled_ind = {}
            for s in scorers:
                per_t, per_t_ac, occ, deg, nulls = [], [], [], [], []
                band_means = {b: [] for b in range(len(BAND_EDGES) - 1)}
                for r in rows:
                    c = r["cells"][arm][rname].get(s)
                    if not c:
                        continue
                    good = [b for b in c["bands"] if b.get("rho") is not None and np.isfinite(b["rho"])
                            and not b.get("degenerate")]
                    occ.append(len(good))
                    deg.append(sum(1 for b in c["bands"] if b.get("degenerate")))
                    if not good:
                        continue
                    per_t.append(float(np.mean([b["rho"] for b in good])))          # band-averaged within target
                    nulls.append(float(np.mean([b["null_mean"] for b in good])))
                    per_t_ac.append(float(c["rho_across"]) if np.isfinite(c["rho_across"]) else np.nan)
                    for b in good:
                        band_means[b["band"]].append(b["rho"])
                if len(per_t) < 3:
                    res[s] = dict(n_targets=len(per_t), verdict="UNDERPOWERED (fewer than 3 targets)")
                    continue
                idx = [q for q, r in enumerate(rows) if r["cells"][arm][rname].get(s) and
                       any(b.get("rho") is not None and np.isfinite(b["rho"]) and not b.get("degenerate")
                           for b in r["cells"][arm][rname][s]["bands"])]
                f_sub = folds[idx]
                x = np.array(per_t)
                cmp_zero = ST.compare(x, np.zeros(len(x)), f_sub, label=f"{s} | {rname} | {arm}: within-band rho vs 0",
                                      seed_parts=("s29Dband",))
                ac = np.array(per_t_ac)
                ok = np.isfinite(ac)
                cmp_ac = ST.compare(x[ok], ac[ok], f_sub[ok], label=f"{s} | {rname} | {arm}: within-band rho - across-band rho",
                                    seed_parts=("s29Dband",)) if ok.sum() >= 3 else None
                pooled_ind[s] = (x, np.array(idx))
                pp = np.array([r["cells"][arm][rname][s].get("partial_pred", np.nan) for r in rows
                               if r["cells"][arm][rname].get(s)], float)
                wc = {}
                for w in WIDTHS:
                    #: v and its fold labels are built in ONE pass so they cannot drift apart
                    #: (they did once: a finiteness filter on the values and a None filter on the
                    #: indices gave 125 values against 126 folds and ST.compare refused, correctly).
                    v, fw = [], []
                    for q, r in enumerate(rows):
                        c = r["cells"][arm][rname].get(s)
                        if not c or "width_curve" not in c:
                            continue
                        val = c["width_curve"].get(str(w), {}).get("rho")
                        if val is None or not np.isfinite(val):
                            continue
                        v.append(float(val)); fw.append(folds[q])
                    if len(v) >= 3:
                        cw = ST.compare(np.array(v), np.zeros(len(v)), np.array(fw), label=f"width {w}",
                                        seed_parts=("s29Dband",))
                        wc[str(w)] = dict(rho=float(np.mean(v)), n=len(v), ci95_fold=cw["ci95_fold"], se=cw["se"])
                    else:
                        wc[str(w)] = dict(rho=None, n=len(v))
                res[s] = dict(
                    n_targets=len(x), mean_rho_in=float(x.mean()), median_rho_in=float(np.median(x)),
                    partial_pred_mean=float(np.nanmean(pp)) if np.isfinite(pp).any() else None,
                    partial_pred_median=float(np.nanmedian(pp)) if np.isfinite(pp).any() else None,
                    rho_SY_mean=float(np.nanmean([r["cells"][arm][rname][s].get("rho_SY", np.nan) for r in rows if r["cells"][arm][rname].get(s)])),
                    rho_SR_mean=float(np.nanmean([r["cells"][arm][rname][s].get("rho_SR", np.nan) for r in rows if r["cells"][arm][rname].get(s)])),
                    rho_RY_mean=float(np.nanmean([r["cells"][arm][rname][s].get("rho_RY", np.nan) for r in rows if r["cells"][arm][rname].get(s)])),
                    width_curve=wc,
                    se=cmp_zero["se"], ci95_fold=cmp_zero["ci95_fold"], ci95_iid=cmp_zero["ci95_iid"],
                    folds_same_sign=cmp_zero["folds_same_sign"], mde=cmp_zero["mde"],
                    effect_over_mde=cmp_zero["effect_over_mde"], verdict_vs_zero=cmp_zero["verdict"],
                    mean_rho_across=float(np.nanmean(ac)), shuffle_null_mean=float(np.mean(nulls)),
                    mean_cells_per_target=float(np.mean(occ)), mean_degenerate_per_target=float(np.mean(deg)),
                    fail18_rho=float(x[[fail[q] for q in idx]].mean()) if any(fail[q] for q in idx) else float("nan"),
                    other_rho=float(x[[not fail[q] for q in idx]].mean()) if any(not fail[q] for q in idx) else float("nan"),
                    per_band={str(b): (float(np.mean(v)) if v else None) for b, v in band_means.items()},
                    per_band_n={str(b): len(v) for b, v in band_means.items()},
                    in_minus_across=(dict(effect=cmp_ac["effect"], ci95_fold=cmp_ac["ci95_fold"],
                                          effect_over_mde=cmp_ac["effect_over_mde"], verdict=cmp_ac["verdict"],
                                          fmt=ST.fmt(cmp_ac)) if cmp_ac else None),
                    fmt_vs_zero=ST.fmt(cmp_zero))
            # the max-over-scorers sign-flip null on the pooled within-band rho
            names = [s for s in scorers if s in pooled_ind]
            if names:
                common = sorted(set.intersection(*(set(pooled_ind[s][1].tolist()) for s in names)))
                if len(common) >= 10:
                    Dm = np.column_stack([pooled_ind[s][0][[list(pooled_ind[s][1]).index(q) for q in common]] for s in names])
                    obs = Dm.mean(0)
                    rng = SD.stable_rng("s29D_band_maxnull", rname, arm, salt=SALT)
                    mx = np.empty(n_max_null)
                    for t in range(n_max_null):
                        eps = np.where(rng.random(len(common)) < 0.5, 1.0, -1.0)
                        mx[t] = (Dm * eps[:, None]).mean(0).max()
                    best = names[int(np.argmax(obs))]
                    res["_max_null"] = dict(n_scorers=len(names), n_targets=len(common), best_scorer=best,
                                            best_mean_rho=float(obs.max()), null_mean=float(mx.mean()),
                                            null_p95=float(np.percentile(mx, 95)),
                                            p_max=float((mx >= obs.max()).mean()))
            out["arms"][arm][rname] = res
    out["text"] = render(out)
    ST.save_atomic(os.path.join(RESULTS, f"s29_D_band_{basis}.json"), out, module_file=__file__)
    print(out["text"])
    print("wrote", os.path.join(RESULTS, f"s29_D_band_{basis}.json"))
    return out


def render(o):
    L = [f"WITHIN-REALISM-BAND ORDERING (ORACLE DIAGNOSTIC; PREREG_S29_D_band), basis {o['basis']}, n={o['n']} targets",
         f"  bands {o['band_edges']}, occupancy >= {o['min_cell']} structures and >= 10 distinct RMSDs, degenerate cells excluded",
         "  positive rho = the scorer orders by ACCURACY inside a band of matched realism"]
    for arm, per_r in o["arms"].items():
        for rname, res in per_r.items():
            L.append(f"  -- arm {arm}, realism {rname} --")
            L.append(f"     {'scorer':14s} {'rho_in':>8s} {'fold CI':>18s} {'pred rho_SY.R':>13s} {'rho_SY':>7s} {'rho_SR':>7s} "
                     f"{'rho_RY':>7s} | width curve rho at w=1.0/0.5/0.25/0.10/0.05 | {'degen/t':>7s} {'n':>4s}")
            for s, v in sorted(res.items()):
                if s.startswith("_") or "mean_rho_in" not in v:
                    continue
                ia = v["in_minus_across"]
                wc = v.get("width_curve", {})
                cur = " ".join(("%+.3f" % wc[str(w)]["rho"]) if wc.get(str(w), {}).get("rho") is not None else "  --  "
                               for w in WIDTHS)
                L.append(f"     {s:14s} {v['mean_rho_in']:+8.4f} [{v['ci95_fold'][0]:+.3f},{v['ci95_fold'][1]:+.3f}] "
                         f"{(v.get('partial_pred_median') if v.get('partial_pred_median') is not None else float('nan')):+13.4f} "
                         f"{v.get('rho_SY_mean', float('nan')):+7.3f} {v.get('rho_SR_mean', float('nan')):+7.3f} "
                         f"{v.get('rho_RY_mean', float('nan')):+7.3f} | {cur} | {v['mean_degenerate_per_target']:7.2f} {v['n_targets']:4d}")
            mn = res.get("_max_null")
            if mn:
                L.append(f"     max-over-{mn['n_scorers']} sign-flip null on {mn['n_targets']} targets: best {mn['best_scorer']} "
                         f"{mn['best_mean_rho']:+.4f}; null mean {mn['null_mean']:+.4f} p95 {mn['null_p95']:+.4f}; p_max {mn['p_max']:.3f}")
    L.append("  Every number ORACLE; nothing deployable. A positive rho_in is a recognition result only if it also")
    L.append("  beats rho_across (F2) and survives the max-over-scorers null (F4); its Angstrom value is the")
    L.append("  band-restricted selection ceiling, computed separately, never asserted (F3).")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["probe", "run", "analyse"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--basis", default="ca")
    a = ap.parse_args()
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    if a.mode == "analyse":
        analyse(a.basis)
    elif a.mode == "probe":
        t0 = time.time()
        row = target_row(pdbs[0])
        c = row["cells"]["pool"]["R1_CAGEO"]["DIS"]
        print(json.dumps({k: v for k, v in row.items() if k != "cells"}, indent=1))
        print("DIS | R1 across %.3f, bands: %s" % (c["rho_across"], [(b["band"], b["n"], b.get("rho")) for b in c["bands"]]))
        print("%.1f s for one target" % (time.time() - t0))
    else:
        run(pdbs, a.basis)


if __name__ == "__main__":
    main()
