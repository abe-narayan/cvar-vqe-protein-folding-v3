"""SPRINT 19, AGENT A, Q1 + Q3 -- pairwise error topology and metric realisability.

DIAGNOSTIC ONLY.  Every statistic on the error axis is ORACLE (it reads `dtrue`); every
feature on the other axis is native-free.  Nothing here is a predictive result.

Q1   where does the predictor fail, as a function of native-free features, and is the error
     ANTI-CORRELATED with the Sprint-18 utility compass (confident + short-range pays)?
Q3   is the predicted distance field even close to a realisable Euclidean distance matrix,
     and does the degree of violation price per-target failure?

Run:  python -m s19.a_topo
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s19 import a_lib as L                   # noqa: E402
from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_topo.json")


def spearman(a, b):
    from scipy.stats import rankdata
    a = rankdata(np.asarray(a, float))
    b = rankdata(np.asarray(b, float))
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else 0.0


def edm_stats(D, i, j, dhat, n):
    """Metric realisability of the predicted distance field.

    The field is only defined for |i-j| >= 2; the chain bonds (|i-j| = 1) are filled with the
    ideal 3.81 A CA-CA spacing, which is what every downstream builder assumes, so the matrix
    is COMPLETE and its Gram spectrum is well defined.  Reported:

      tri_viol   fraction of triples (a,b,c) with D_ab > D_ac + D_cb  (a metric-space test)
      edm_defect sum |lambda_k| for k > 3 over sum |lambda|, of  -1/2 J D^2 J  (rank-3 defect)
      neg_mass   sum of |negative eigenvalues| over sum |lambda|   (non-Euclidean mass)
    """
    M = np.zeros((n, n))
    M[i, j] = dhat
    M[j, i] = dhat
    k = np.arange(n - 1)
    M[k, k + 1] = 3.81
    M[k + 1, k] = 3.81
    # triangle inequality over all triples
    S = M[:, :, None] + M[None, :, :]           # D_ac + D_cb  -> (a, c, b)
    viol = M[:, None, :] > S + 1e-9             # D_ab vs D_ac + D_cb
    mask = np.ones((n, n, n), bool)
    idx = np.arange(n)
    mask[idx, idx, :] = False
    mask[idx, :, idx] = False
    mask[:, idx, idx] = False
    tri = float(viol[mask].mean())
    J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ (M ** 2) @ J
    ev = np.linalg.eigvalsh(G)
    a = np.abs(ev)
    order = np.argsort(-a)
    tot = a.sum()
    defect = float(a[order[3:]].sum() / max(tot, 1e-12))
    neg = float(a[ev < 0].sum() / max(tot, 1e-12))
    return tri, defect, neg


def main():
    rng = SD.stable_rng("s19", "agentA", "topo")
    tg = I.targets()
    data, pdbs, folds = L.gather_all(tg)
    ref = {r["pdb"]: r for r in
           json.load(open(os.path.join(ROOT, "s18", "results", "objceil.json")))["rows"]}

    # -------------------------------------------------- flat per-pair table
    keys = ["sep", "term", "sd", "ent", "maxp", "nmode", "near1", "dhat", "deg_m",
            "pred_contact", "hyd_m", "chg_prod", "vol_m", "hel_m", "she_m",
            "ss_same", "ss_H", "ss_E", "w"]
    flat = {k: np.concatenate([data[p][k] for p in pdbs]) for k in keys}
    flat["absr"] = np.concatenate([data[p]["absr"] for p in pdbs])
    flat["r"] = np.concatenate([data[p]["r"] for p in pdbs])
    flat["nlen"] = np.concatenate([np.full(len(data[p]["sep"]), data[p]["n"], float)
                                   for p in pdbs])
    npairs = len(flat["absr"])

    out = {"n_targets": len(pdbs), "n_pairs": int(npairs)}

    # -------------------------------------------------- Q1: correlations
    print(f"\n=== Q1  ERROR TOPOLOGY   n = {npairs} pairs over {len(pdbs)} targets ===")
    print("  ORACLE axis: |r| = |dhat - dtrue|.  All features native-free.\n")
    print(f"  {'feature':<16}{'rho(|r|,x)':>12}{'rho(r,x)':>12}")
    cors = {}
    for k in keys + ["nlen"]:
        ra = spearman(flat["absr"], flat[k])
        rs = spearman(flat["r"], flat[k])
        cors[k] = {"abs": ra, "signed": rs}
        print(f"  {k:<16}{ra:>12.3f}{rs:>12.3f}")
    out["spearman"] = cors

    # the two compass axes, stated explicitly
    print(f"\n  COMPASS CHECK   rho(|r|, sep)      = {cors['sep']['abs']:+.3f}"
          "   (A1 predicts > +0.30)")
    print(f"                  rho(|r|, 1/sd^2)   = {cors['w']['abs']:+.3f}"
          "   (A1 predicts < -0.30; D9 recorded -0.476)")

    # -------------------------------------------------- Q1: strata
    def strat(name, groups):
        print(f"\n  {name}")
        print(f"    {'stratum':<18}{'n':>7}{'MAE':>8}{'bias':>8}{'RMS':>8}"
              f"{'mean sd':>9}{'share of SSE':>14}")
        sse_tot = (flat["absr"] ** 2).sum()
        rows = []
        for lab, m in groups:
            if m.sum() == 0:
                continue
            a = flat["absr"][m]
            rr = flat["r"][m]
            row = {"stratum": lab, "n": int(m.sum()), "mae": float(a.mean()),
                   "bias": float(rr.mean()), "rms": float(np.sqrt((a ** 2).mean())),
                   "sd": float(flat["sd"][m].mean()),
                   "share_sse": float((a ** 2).sum() / sse_tot)}
            rows.append(row)
            print(f"    {lab:<18}{row['n']:>7}{row['mae']:>8.3f}{row['bias']:>+8.3f}"
                  f"{row['rms']:>8.3f}{row['sd']:>9.3f}{row['share_sse']:>14.3f}")
        return rows

    sep = flat["sep"]
    out["by_sep"] = strat("BY SEQUENCE SEPARATION", [
        ("2-3", (sep <= 3)), ("4-5", (sep >= 4) & (sep <= 5)),
        ("6-7", (sep >= 6) & (sep <= 7)), ("8-10", (sep >= 8) & (sep <= 10)),
        ("11+", sep >= 11)])
    q = np.quantile(flat["sd"], [0.25, 0.5, 0.75])
    out["by_conf"] = strat("BY PREDICTED CONFIDENCE (sd quartile; Q1 = most confident)", [
        ("Q1 sd<%.2f" % q[0], flat["sd"] <= q[0]),
        ("Q2", (flat["sd"] > q[0]) & (flat["sd"] <= q[1])),
        ("Q3", (flat["sd"] > q[1]) & (flat["sd"] <= q[2])),
        ("Q4 sd>%.2f" % q[2], flat["sd"] > q[2])])
    out["by_modes"] = strat("BY NUMBER OF PREDICTIVE MODES", [
        ("unimodal", flat["nmode"] <= 1), ("2 modes", flat["nmode"] == 2),
        (">=3 modes", flat["nmode"] >= 3)])
    out["by_term"] = strat("BY TERMINAL PROXIMITY min(i, n-1-j)", [
        ("0 (terminal)", flat["term"] == 0), ("1", flat["term"] == 1),
        ("2", flat["term"] == 2), (">=3 (core)", flat["term"] >= 3)])
    out["by_ss"] = strat("BY CONSENSUS SS OF THE PAIR (native-free, top-75 circmean)", [
        ("H-H", flat["ss_H"] == 1), ("E-E", flat["ss_E"] == 1),
        ("same class", flat["ss_same"] == 1), ("mixed", flat["ss_same"] == 0)])
    out["by_contact"] = strat("BY PREDICTED CONTACT STATUS", [
        ("dhat < 8 A", flat["pred_contact"] == 1), ("dhat >= 8 A", flat["pred_contact"] == 0)])

    # -------------------------------------------------- Q1: joint utility x error
    print("\n  THE COMPASS CROSS-TAB: mean |r| in the four compass cells")
    short = sep <= 5
    conf = flat["sd"] <= np.median(flat["sd"])
    cells = {}
    print(f"    {'':<22}{'confident':>12}{'unconfident':>14}")
    for lab, sm in (("short (sep<=5)", short), ("long  (sep>=6)", ~short)):
        vals = []
        for cm in (conf, ~conf):
            m = sm & cm
            vals.append(float(flat["absr"][m].mean()) if m.sum() else float("nan"))
            cells[f"{lab.split()[0]}_{'conf' if cm is conf else 'unconf'}"] = vals[-1]
        print(f"    {lab:<22}{vals[0]:>12.3f}{vals[1]:>14.3f}")
    # SSE share landing in the HIGH-utility cell
    hi = short & conf
    out["compass_cells"] = cells
    out["high_utility_pair_share"] = float(hi.mean())
    out["high_utility_sse_share"] = float((flat["absr"][hi] ** 2).sum()
                                          / (flat["absr"] ** 2).sum())
    print(f"\n    high-utility cell (short AND confident) holds "
          f"{hi.mean():.1%} of pairs and {out['high_utility_sse_share']:.1%} of the squared error")

    # -------------------------------------------------- Q1: multivariate
    X = np.column_stack([flat[k] for k in keys] + [flat["nlen"], np.ones(npairs)])
    Xs = (X - X.mean(0)) / np.maximum(X.std(0), 1e-9)
    Xs[:, -1] = 1.0
    y = flat["absr"]
    coef, *_ = np.linalg.lstsq(Xs, y, rcond=None)
    pred = Xs @ coef
    r2 = 1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    print(f"\n  LINEAR MODEL for |r| from native-free features: R^2 = {r2:.3f}")
    ordr = np.argsort(-np.abs(coef[:-1]))
    names = keys + ["nlen"]
    for k in ordr[:8]:
        print(f"    {names[k]:<16}{coef[k]:+.3f}")
    out["linear_absr"] = {"r2": float(r2),
                          "coef": {names[k]: float(coef[k]) for k in range(len(names))}}
    # how well can a native-free model RANK the error?
    out["nf_rank_skill"] = spearman(pred, y)
    print(f"  native-free ranking skill for |r| (Spearman): {out['nf_rank_skill']:+.3f}")

    # -------------------------------------------------- Q3: realisability
    print("\n=== Q3  METRIC REALISABILITY OF THE PREDICTED FIELD ===")
    rows = []
    for p in pdbs:
        d = data[p]
        n = d["n"]
        tri, defect, neg = edm_stats(None, d["i"], d["j"], d["dhat"], n)
        tri_t, def_t, neg_t = edm_stats(None, d["i"], d["j"], d["dtrue"], n)
        rows.append({"pdb": p, "n": n, "fold": d["fold"],
                     "tri": tri, "defect": defect, "neg": neg,
                     "tri_true": tri_t, "defect_true": def_t, "neg_true": neg_t,
                     "a0": float(ref[p]["a0.0"]), "avg": float(ref[p]["avg"]),
                     "resid_rms": float(ref[p]["resid_rms"]),
                     "mae": float(d["absr"].mean())})
    g = lambda k: np.array([r[k] for r in rows])                       # noqa: E731
    print(f"  {'quantity':<28}{'predicted':>12}{'native (ORACLE)':>18}")
    for k, lab in (("tri", "triangle violations"), ("defect", "rank-3 EDM defect"),
                   ("neg", "non-Euclidean mass")):
        print(f"  {lab:<28}{g(k).mean():>12.4f}{g(k + '_true').mean():>18.4f}")
    print(f"\n  {'rho(x, a0.0 RMSD)':<28}")
    for k in ("defect", "neg", "tri", "mae", "resid_rms"):
        print(f"    {k:<20}{spearman(g(k), g('a0')):+.3f}")
    out["realisability"] = {
        "mean": {k: float(g(k).mean()) for k in
                 ("tri", "defect", "neg", "tri_true", "defect_true", "neg_true")},
        "rho_vs_a0": {k: spearman(g(k), g("a0")) for k in
                      ("defect", "neg", "tri", "mae", "resid_rms")},
        "per_target": rows}

    json.dump(out, open(OUT, "w"), indent=1, default=float)
    open(os.path.join(L.RESULTS, "a_topo.COMPLETE"), "w").write("ok\n")
    print(f"\nwrote {OUT}")
    return out


if __name__ == "__main__":
    main()
