#!/usr/bin/env python
"""s31/s31_F3_prefix.py -- S31 LANE F, F3: is `bestm128 = 2.9027` an order statistic?

Pre-registered in `s31/PREREG_S31_F.md` s11 (fourth amendment), committed before any aggregate
of this family existed.

`bestm128` is a PER-TARGET MINIMUM OVER K = 128 prefix lengths -- the exact shape contract rule 11
exists for.  The built-chain 2.9027 is a CLOUD-SELECTED oracle, projected: `s29_O_ladder.py:542`
picks `m_best` on the cloud curve and projects that one structure.  Said in every sentence.

Arms (prereg s11.2): F3-a ORACLE gain | F3-b zero-skill random m | F3-c ORACLE global m + LFO |
F3-d split-half transfer | F3-e MATCHED RANDOM-VARIANT-FAMILY control | F3-f native-free
fold-held-out m rules | F3-g effective K.

    python s31/s31_F3_prefix.py cheap      # everything except F3-e
    python s31/s31_F3_prefix.py e          # F3-e (recomputes 128 random subsets per target)
    python s31/s31_F3_prefix.py analyse
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
P128 = os.path.join(ROOT, "s29", "results", "s29_O_p128_rows.jsonl")
EROWS = os.path.join(RESULTS, "s31_F3_randfamily_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_F3_prefix.json")
SEED = 31007
K128 = 128
NDRAW = 4                    # matched-random-family draw seeds (contract rule 10)
NCOMP = [0]


def load_p128():
    R = [json.loads(ln) for ln in open(P128) if ln.strip()]
    R.sort(key=lambda r: r["pdb"])
    assert len(R) == 126, len(R)
    return R


def load_chain():
    """(pdb, item) -> built-chain RMSD, from every s29 O-ladder chain shard."""
    D = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "s29", "results", "s29_O_chain_rows*.jsonl"))):
        for ln in open(f):
            if not ln.strip():
                continue
            r = json.loads(ln)
            D[(r["pdb"], r["item"])] = (float(r["rmsd_chain"]), float(r["rmsd_cloud"]))
    return D


def cmp2(a, b, folds, names, label):
    NCOMP[0] += 1
    return ST.compare(np.asarray(a, float), np.asarray(b, float), folds=folds, names=names,
                      label=label, seed_parts=("s31F3", str(SEED)))


def brief(o):
    return dict(label=o["label"], mean_a=o["mean_a"], mean_b=o["mean_b"], effect=o["effect"],
                median_effect=o["median_effect"], se=o["se"], mde=o["mde"],
                effect_over_mde=o["effect_over_mde"], ci95_fold=o["ci95_fold"],
                folds_same_sign=o["folds_same_sign"], W=o["n_better"], L=o["n_worse"], n=o["n"])


def verdict(o):
    r = abs(o["effect_over_mde"]); ci = o["ci95_fold"]
    ex = bool(ci is not None and (ci[0] > 0 or ci[1] < 0))
    if r >= 1.0 and ex:
        return "MEASURED"
    if r >= 1.0:
        return "NOT MEASURED (>=1.0x MDE, fold CI spans zero)"
    if r >= 0.7:
        return "NOT MEASURED (0.7-1.0x MDE)"
    return "NOT A RESULT (<0.7x MDE)"


# --------------------------------------------------------------------------- F3-e
def run_randfamily():
    """MATCHED RANDOM-VARIANT-FAMILY control, in the operator's own space (contract rule 7).

    For each target: the SAME top-128, the SAME operator (superpose on the subset medoid, uniform
    coordinate mean), the SAME K = 128 variants and the SAME size distribution 1..128 -- but each
    variant is a RANDOM subset of that size instead of the score-ordered PREFIX.  So the variant
    index carries no score information and nothing else changes.  NDRAW independent draws, and the
    MEAN over draws is reported, never the maximum (contract rule 10).
    """
    done = set()
    if os.path.exists(EROWS):
        for ln in open(EROWS):
            if ln.strip():
                done.add(json.loads(ln)["pdb"])
    P = {r["pdb"]: r for r in load_p128()}
    for t in I.targets():
        pdb = t["pdb"]
        if pdb in done:
            continue
        u = I.load_univ(pdb)
        dg = I.distogram(pdb)
        n = int(u["n"])
        pool = np.asarray(u["order"], int)[:500]
        W = np.asarray(u["W"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        ii, jj = I.pair_index(n, 2)
        D = np.linalg.norm(W[:, ii, :] - W[:, jj, :], axis=2)
        order = np.argsort(I.shipped_score(dg, D), kind="stable")
        top = order[:K128]
        Wt = W[top]
        Pt = I.pairwise_rmsd(Wt)
        # reproduction gate: the prefix curve must reproduce s29's stored curve
        cur = np.empty(K128)
        for m in range(1, K128 + 1):
            S = np.arange(m)
            b = I.medoid(Pt[np.ix_(S, S)])
            cur[m - 1] = I.ca_rmsd(I.superpose_batch(Wt[S], Wt[S[b]]).mean(0), nat)
        ref = np.asarray(P[pdb]["curve128"], float)
        dev = float(np.abs(cur - ref).max())
        if dev > 1e-3:
            raise RuntimeError("F3 REPRODUCTION GATE FAILED on %s: max|curve-ref| = %.6f"
                               % (pdb, dev))
        rows_rand = []
        for d in range(NDRAW):
            rng = np.random.default_rng(SEED + 1000 * d + abs(hash(pdb)) % 997)
            v = np.empty(K128)
            for m in range(1, K128 + 1):
                S = rng.choice(K128, m, replace=False) if m < K128 else np.arange(K128)
                b = I.medoid(Pt[np.ix_(S, S)])
                v[m - 1] = I.ca_rmsd(I.superpose_batch(Wt[S], Wt[S[b]]).mean(0), nat)
            rows_rand.append(v.tolist())
        row = dict(pdb=pdb, n=n, fold=int(t["fold"]), curve_repro_maxdev=dev,
                   curve_prefix=cur.tolist(), curves_random=rows_rand,
                   disp128=float((Pt.sum() - np.trace(Pt)) / (K128 * (K128 - 1))),
                   rg_sd128=float(np.sqrt(((Wt - Wt.mean(1, keepdims=True)) ** 2)
                                          .sum(2).mean(1)).std(ddof=1)),
                   n_distinct=int(I.shipped_record(pdb).get("n_distinct", -1)))
        with open(EROWS, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        print("%s n=%2d repro %.2e  prefix-min %.3f  rand-min(mean over %d draws) %.3f"
              % (pdb, n, dev, cur.min(), NDRAW,
                 float(np.mean([min(v) for v in rows_rand]))), flush=True)
    print("F3-e complete")


# --------------------------------------------------------------------------- analysis
def analyse():
    R = load_p128()
    CH = load_chain()
    names = [r["pdb"] for r in R]
    folds = np.array([r["fold"] for r in R])
    M = np.array([r["curve128"] for r in R], float)            # 126 x 128, CA POINT CLOUD
    prod_cloud = np.array([r["prod"] for r in R], float)
    m75 = M[:, 74]
    oracle_cloud = M.min(1)
    argm = M.argmin(1) + 1
    out = dict(seed=SEED, n=len(R),
               prereg="s31/PREREG_S31_F.md s11 @ this commit",
               source_rows=P128,
               BASIS_NOTE="the 126x128 matrix is CA POINT CLOUD. The built-chain 2.9027 is a "
                          "CLOUD-SELECTED oracle, projected (s29_O_ladder.py:542 picks m_best on "
                          "the cloud curve and projects that one structure).")

    # ---- chain arms from s29 (same job, same stored clouds)
    have = all((p, it) in CH for p in names for it in ("prod", "bestm128"))
    if have:
        ch_prod = np.array([CH[(p, "prod")][0] for p in names])
        ch_bestm = np.array([CH[(p, "bestm128")][0] for p in names])
        o = cmp2(ch_bestm, ch_prod, folds, names, "F3a.bestm128 - production (BUILT CHAIN, ORACLE)")
        out["F3a_chain"] = dict(cmp=brief(o), verdict=verdict(o),
                                ORACLE="ORACLE / NOT DEPLOYABLE -- a per-target minimum over "
                                       "K = 128 prefix lengths, selected on the cloud",
                                mean_bestm128=float(ch_bestm.mean()),
                                mean_prod=float(ch_prod.mean()))
    else:
        out["F3a_chain"] = dict(error="s29 chain rows do not cover every (pdb, item)")

    o = cmp2(oracle_cloud, m75, folds, names, "F3a.bestm128 - m=75 (CA POINT CLOUD, ORACLE)")
    out["F3a_cloud"] = dict(cmp=brief(o), verdict=verdict(o),
                            mean_oracle=float(oracle_cloud.mean()), mean_m75=float(m75.mean()),
                            mean_prod_record=float(prod_cloud.mean()),
                            argmin_m=dict(mean=float(argm.mean()), median=float(np.median(argm)),
                                          hist=np.histogram(argm, bins=[1, 2, 5, 10, 20, 40, 75,
                                                                        100, 129])[0].tolist(),
                                          bins="[1,2) [2,5) [5,10) [10,20) [20,40) [40,75) "
                                               "[75,100) [100,128]"))

    # ---- F3-b zero-skill random m (its own draw distribution, contract rule 10)
    rng = np.random.default_rng(SEED)
    draws = np.array([M[np.arange(len(R)), rng.integers(0, K128, len(R))].mean()
                      for _ in range(2000)])
    out["F3b_zero_skill_random_m"] = dict(
        basis="CA POINT CLOUD",
        mean_over_draws=float(draws.mean()), sd_over_draws=float(draws.std(ddof=1)),
        pct=[float(np.percentile(draws, q)) for q in (2.5, 50, 97.5)],
        row_mean_over_m=float(M.mean(1).mean()),
        gain_vs_m75=float(draws.mean() - m75.mean()),
        note="a rule with NO skill picking m uniformly; reported as a distribution, not a maximum")

    # ---- F3-c ORACLE global m, and its leave-fold-out twin
    gm = M.mean(0)
    m_star = int(gm.argmin()) + 1
    lfo = np.empty(len(R))
    m_by_fold = {}
    for f in sorted(set(folds.tolist())):
        tr = folds != f
        mk = int(M[tr].mean(0).argmin())
        m_by_fold[int(f)] = mk + 1
        lfo[folds == f] = M[folds == f, mk]
    o_g = cmp2(M[:, m_star - 1], m75, folds, names,
               "F3c.ORACLE global m - m=75 (CA POINT CLOUD)")
    o_l = cmp2(lfo, m75, folds, names, "F3c.leave-fold-out global m - m=75 (CA POINT CLOUD)")
    out["F3c_global_m"] = dict(m_star=m_star, mean_at_m_star=float(gm[m_star - 1]),
                               mean_at_75=float(gm[74]), m_by_fold=m_by_fold,
                               ORACLE_global=dict(cmp=brief(o_g), verdict=verdict(o_g)),
                               LFO=dict(cmp=brief(o_l), verdict=verdict(o_l)))

    # ---- F3-d split-half transfer over the 128 columns
    out["F3d_split_half_transfer"] = ST.split_half_transfer(M, seed_parts=("s31F3", str(SEED)))
    out["F3d_note"] = ("stats_lib.split_half_transfer chooses ONE GLOBAL column on half the "
                       "targets and scores it on the other half: it prices the m axis as a "
                       "GLOBAL setting, which is the only thing a split over targets can price. "
                       "It is NOT a test of a per-target rule -- F3-f is.")

    # ---- F3-e matched random-variant-family control
    if os.path.exists(EROWS):
        E = {json.loads(ln)["pdb"]: json.loads(ln) for ln in open(EROWS) if ln.strip()}
        if len(E) == 126:
            pref = np.array([E[p]["curve_prefix"] for p in names], float)
            rmin = np.array([[min(v) for v in E[p]["curves_random"]] for p in names], float)
            per_draw = [cmp2(rmin[:, d], m75, folds, names,
                             "F3e.random-family min (draw %d) - m=75 (CLOUD)" % d)
                        for d in range(rmin.shape[1])]
            gains = np.array([o["effect"] for o in per_draw])
            g_prefix = float(oracle_cloud.mean() - m75.mean())
            g_rand = float(gains.mean())
            out["F3e_matched_random_family"] = dict(
                construction="same top-128, same operator, same K=128, same size distribution "
                             "1..128; each variant a RANDOM subset of that size instead of the "
                             "score-ordered prefix",
                n_draws=int(rmin.shape[1]),
                gain_prefix=g_prefix, gain_random_mean=g_rand,
                gain_random_per_draw=gains.tolist(),
                gain_random_sd=float(gains.std(ddof=1)),
                share_of_prefix_gain=float(g_rand / g_prefix) if g_prefix else float("nan"),
                REGISTERED_BAR="'the m axis is an order statistic' FIRES if the random family "
                               "reaches >= 80%% of the prefix gain",
                BAR_FIRES=bool(g_prefix and (g_rand / g_prefix) >= 0.80),
                prefix_curve_repro_maxdev=float(max(E[p]["curve_repro_maxdev"] for p in names)))
        else:
            out["F3e_matched_random_family"] = dict(status="INCOMPLETE %d/126" % len(E))
    else:
        out["F3e_matched_random_family"] = dict(status="not run")

    # ---- F3-f native-free per-target m rules, fold-held-out
    feats = {}
    if os.path.exists(EROWS):
        E = {json.loads(ln)["pdb"]: json.loads(ln) for ln in open(EROWS) if ln.strip()}
        if len(E) == 126:
            feats["DISP128"] = np.array([E[p]["disp128"] for p in names], float)
            feats["rg_sd128"] = np.array([E[p]["rg_sd128"] for p in names], float)
            feats["n_distinct"] = np.array([E[p]["n_distinct"] for p in names], float)
    feats["n"] = np.array([r["n"] for r in R], float)
    TROWS = os.path.join(RESULTS, "s31_F_terminal_rows.jsonl")
    if os.path.exists(TROWS):
        T = {json.loads(ln)["pdb"]: json.loads(ln) for ln in open(TROWS) if ln.strip()}
        if len(T) == 126:
            feats["DISP75"] = np.array([T[p]["DISP"] for p in names], float)
            feats["MOVE_AVG"] = np.array([T[p]["move_AVG"] for p in names], float)

    GRID_A = np.arange(1, K128 + 1)
    rules = {}
    for fn, fv in feats.items():
        z = (fv - fv.mean()) / max(fv.std(ddof=1), 1e-12)
        best = None
        pred = np.empty(len(R))
        for f in sorted(set(folds.tolist())):
            tr = folds != f
            bv, ba, bb = None, 75, 0.0
            for a in GRID_A:
                for b in (-48, -32, -20, -12, -6, 0, 6, 12, 20, 32, 48):
                    mm = np.clip(np.round(a + b * z[tr]).astype(int), 1, K128) - 1
                    v = M[np.where(tr)[0], mm].mean()
                    if bv is None or v < bv:
                        bv, ba, bb = v, int(a), float(b)
            mm = np.clip(np.round(ba + bb * z[folds == f]).astype(int), 1, K128) - 1
            pred[folds == f] = M[np.where(folds == f)[0], mm]
            best = (ba, bb)
        o = cmp2(pred, m75, folds, names, "F3f.%s LFO m-rule - m=75 (CA POINT CLOUD)" % fn)
        rules[fn] = dict(cmp=brief(o), verdict=verdict(o), last_fold_fit=best,
                         DEPLOYABLE="native-free feature, threshold fitted LEAVE-FOLD-OUT on the "
                                    "native cloud RMSD -- the FIT touches the native, so this is "
                                    "an upper bound on a deployable rule, not a clean one")
    out["F3f_native_free_rules"] = rules
    if rules:
        bb = min(rules.items(), key=lambda kv: kv[1]["cmp"]["effect"])
        g_or = float(oracle_cloud.mean() - m75.mean())
        out["F3f_summary"] = dict(
            best_rule=bb[0], best_effect=bb[1]["cmp"]["effect"],
            best_over_mde=bb[1]["cmp"]["effect_over_mde"],
            share_of_oracle=float(bb[1]["cmp"]["effect"] / g_or) if g_or else float("nan"),
            REGISTERED_BAR="'no part of it is deployable' FIRES if the best rule is "
                           "> -0.7 x its own MDE",
            BAR_FIRES=bool(bb[1]["cmp"]["effect"] > -0.7 * bb[1]["cmp"]["mde"]))

    # ---- F3-g effective K
    d = np.diff(M, axis=1)
    loc_min = np.array([int(((M[i, 1:-1] < M[i, :-2]) & (M[i, 1:-1] < M[i, 2:])).sum())
                        for i in range(len(R))], float)
    within = {str(t): float((M <= (oracle_cloud[:, None] + t)).sum(1).mean())
              for t in (0.001, 0.01, 0.05, 0.1)}
    ac1 = np.array([float(np.corrcoef(d[i, :-1], d[i, 1:])[0, 1]) for i in range(len(R))])
    out["F3g_effective_K"] = dict(
        n_local_minima=dict(mean=float(loc_min.mean()), median=float(np.median(loc_min)),
                            max=float(loc_min.max())),
        n_m_within_tol_of_min=within,
        lag1_autocorr_of_curve_increments=float(np.nanmean(ac1)),
        curve_range_mean=float((M.max(1) - M.min(1)).mean()),
        note="a SMOOTH curve has few effectively independent choices, so the best-of-K inflation "
             "is priced at K_eff, not at 128; many m within 0.01 of the min means the argmin is "
             "not sharply identified")

    out["multiplicity"] = dict(registered_families=7, comparisons_emitted=int(NCOMP[0]))
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float)[:9000])
    print("\nwrote", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["e", "analyse"])
    a = ap.parse_args()
    if a.phase == "e":
        run_randfamily()
    else:
        analyse()
