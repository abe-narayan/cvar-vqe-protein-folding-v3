"""FAIL18 forensics, step 11 (coordinator request 2 & 3): fusion sweep + native-free router.

Part A: fuse the shipped distogram score with the best contact-agreement forms
        (`topk`, `corr3`) and with the diversity lever, sweep the mixing weight, and
        report top-75 best / recall on all 126 with the foreign-map null alongside.
Part B: a LEAVE-FOLD-OUT native-free ROUTER for "is this target in the FAIL regime?",
        using only deployable features, with the null (does it beat the base rate?).
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I
from s12.fail_esmrescore import load_con, z
from s12.fail_contactforms import all_forms
from s12.fail_diversity import greedy_div
import protein_geometry as geo

WGRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]


def features(pdb, n, fold, seq, Wc, sc, dg, con, rr):
    """Native-free per-target features for the router."""
    prob = np.asarray(dg["prob"], float)
    ent = float(-(prob * np.log(prob + 1e-9)).sum(1).mean())
    sd = float(np.asarray(dg["sd"], float).mean())
    exp = np.asarray(dg["expected"], float)
    i, j = I.pair_index(n)
    Dfull = np.zeros((n, n)); Dfull[i, j] = exp; Dfull[j, i] = exp
    k1 = np.arange(n - 1); Dfull[k1, k1 + 1] = 3.81; Dfull[k1 + 1, k1] = 3.81
    rg_pred = float(np.sqrt((Dfull ** 2).sum() / (2.0 * n * n)))
    c = Wc - Wc.mean(1, keepdims=True)
    rgs = np.sqrt((c ** 2).sum(2).mean(1))
    top = np.argsort(sc, kind="stable")[:I.M]
    P = I.pairwise_rmsd(Wc[top])
    spread75 = float(P[np.triu_indices(I.M, 1)].mean())
    m3 = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    csharp = float(con[m3].max()) if m3.any() else 0.0
    cmass = float(con[m3].sum() / n) if m3.any() else 0.0
    return dict(n=float(n), dg_entropy=ent, dg_sd=sd,
                rg_pred=rg_pred, rg_pool_med=float(np.median(rgs)),
                rg_pred_minus_pool=rg_pred - float(np.median(rgs)),
                rg_pred_minus_top75=rg_pred - float(np.median(rgs[top])),
                spread75=spread75, spread_ratio=spread75 / (float(np.std(rgs)) + 1e-6),
                esm_sharp=csharp, esm_mass=cmass,
                score_gap=float(np.sort(sc)[74] - np.sort(sc)[0]),
                score_sd=float(sc.std()))


def main():
    con_bank = load_con()
    tg = I.targets(); seqs = [t["seq"] for t in tg]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]; band = rr <= rr.min() + I.BAND
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        sc = I.shipped_score(dg, I.pair_dists(Wc, i, j).astype(np.float32).astype(float))
        bb = geo.build_backbone_batch(u["PHI"][p], u["PSI"][p])
        CB = np.asarray(bb.get("CB", bb["CA"]), float)
        Dca = np.linalg.norm(Wc[:, :, None] - Wc[:, None], axis=-1)
        Dcb = np.linalg.norm(CB[:, :, None] - CB[:, None], axis=-1)
        alt = [s for s in seqs if len(s) == n and s != seq]
        own = all_forms(Dca, Dcb, con_bank[seq], n)
        nul = all_forms(Dca, Dcb, con_bank[alt[k % len(alt)]], n)
        div = greedy_div(Wc, sc, m=I.M)
        rec = dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                   pool_best=float(rr.min()), arms={},
                   feat=features(pdb, n, fold, seq, Wc, sc, dg, con_bank[seq], rr))
        for name, v in (("topk", own["topk"]), ("corr3", own["corr3"]),
                        ("topk_null", nul["topk"]), ("corr3_null", nul["corr3"])):
            for w in WGRID:
                s = z(-sc) + w * z(v)
                sub = np.argsort(-s, kind="stable")[:I.M]
                rec["arms"][f"{name}@{w}"] = dict(top75_best=float(rr[sub].min()),
                                                  top75_mean=float(rr[sub].mean()),
                                                  recall=int(band[sub].any()))
        rec["arms"]["divmax"] = dict(top75_best=float(rr[div].min()),
                                     top75_mean=float(rr[div].mean()),
                                     recall=int(band[div].any()))
        rows.append(rec)
        if (k + 1) % 30 == 0:
            print(f"  {k+1}/126 free={I.free_gb():.1f}", flush=True)
    I.write("fail_fusion", rows)

    f = [r for r in rows if r["fail18"]]; o = [r for r in rows if not r["fail18"]]
    base = np.array([r["arms"]["topk@0.0"]["top75_best"] for r in rows])
    folds = np.array([r["fold"] for r in rows])
    print(f"\n== fusion sweep: top-75 best (ORACLE eval), all 126 ==")
    print(f"{'arm':16s}{'all126':>9s}{'FAIL18':>9s}{'other108':>10s}{'recallF':>9s}"
          f"{'d126':>8s}{'ci95':>19s}{'drop10':>8s}")
    summ = {}
    for a in list(rows[0]["arms"]):
        A = np.array([r["arms"][a]["top75_best"] for r in rows])
        AF = np.array([r["arms"][a]["top75_best"] for r in f])
        AO = np.array([r["arms"][a]["top75_best"] for r in o])
        rF = np.mean([r["arms"][a]["recall"] for r in f])
        st = I.paired(A, base, folds=folds)
        summ[a] = dict(all126=float(A.mean()), fail18=float(AF.mean()),
                       other108=float(AO.mean()), recall_fail=float(rF), paired=st)
        print(f"{a:16s}{A.mean():9.3f}{AF.mean():9.3f}{AO.mean():10.3f}{rF:9.2f}"
              f"{st['mean_diff']:8.3f} [{st['ci95'][0]:6.3f},{st['ci95'][1]:6.3f}]"
              f"{(st['drop_top10_mean_diff'] or 0):8.3f}")

    # ---------------- Part B: native-free router, leave-fold-out ------------------
    print("\n== ROUTER: predict FAIL-regime membership out of fold, native-free features ==")
    fkeys = sorted(rows[0]["feat"])
    X = np.array([[r["feat"][k] for k in fkeys] for r in rows], float)
    y = np.array([r["fail18"] for r in rows], float)
    fold = np.array([r["fold"] for r in rows])
    Xz = (X - X.mean(0)) / (X.std(0) + 1e-9)
    print("  univariate |r| with FAIL18 label:")
    for k, kk in enumerate(fkeys):
        print(f"    {kk:22s} r={np.corrcoef(Xz[:,k], y)[0,1]:+.3f}")
    from sklearn.linear_model import LogisticRegression
    pred = np.zeros(len(y))
    for fo in np.unique(fold):
        tr = fold != fo
        m = LogisticRegression(C=0.5, max_iter=2000).fit(Xz[tr], y[tr])
        pred[~tr] = m.predict_proba(Xz[~tr])[:, 1]
    from scipy.stats import rankdata
    r = rankdata(pred); n1 = int(y.sum()); n0 = len(y) - n1
    ra = float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    print(f"  LFO logistic AUC (FAIL18 vs rest) = {ra:.3f}   base rate = {n1/len(y):.3f}")
    rng = np.random.default_rng(0)
    nulls = []
    for _ in range(200):
        ys = rng.permutation(y); pr = np.zeros(len(y))
        for fo in np.unique(fold):
            tr = fold != fo
            if ys[tr].sum() < 2:
                pr[~tr] = 0.5; continue
            m = LogisticRegression(C=0.5, max_iter=2000).fit(Xz[tr], ys[tr])
            pr[~tr] = m.predict_proba(Xz[~tr])[:, 1]
        rr2 = rankdata(pr)
        nulls.append((rr2[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    nulls = np.array(nulls)
    print(f"  label-permutation null AUC = {nulls.mean():.3f} "
          f"[{np.percentile(nulls,2.5):.3f},{np.percentile(nulls,97.5):.3f}]  "
          f"p = {(nulls >= ra).mean():.3f}")
    # precision at the top-18
    top18 = np.argsort(-pred)[:18]
    print(f"  precision@18 = {y[top18].mean():.3f} (base rate {n1/len(y):.3f})")
    I.write("fail_fusion_summary", dict(sweep=summ, router=dict(
        auc=ra, null_mean=float(nulls.mean()),
        null_ci=[float(np.percentile(nulls, 2.5)), float(np.percentile(nulls, 97.5))],
        p=float((nulls >= ra).mean()), precision_at_18=float(y[top18].mean()),
        features=fkeys, pred={r["pdb"]: float(pp) for r, pp in zip(rows, pred)})))


if __name__ == "__main__":
    main()
