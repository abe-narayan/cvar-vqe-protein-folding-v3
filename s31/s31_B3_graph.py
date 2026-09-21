#!/usr/bin/env python
"""s31/s31_B3_graph.py -- lane B: the CHEAP PRE-CHECK on the coordinator's multi-structure
escape, before anything is built on it.

THE OPENING. G1 is a theorem about SINGLE-structure observables. A quantity that is a function of
the SET of candidate distance maps rather than of any one of them escapes it on exactly that
argument. A free-energy BARRIER between two pool basins is such a quantity but costs 10^3-10^5
trajectories per target. The coordinator's proposal: connectivity observables on the candidate
GRAPH -- geodesic, commute time, spectral quantities on the pairwise Kabsch-RMSD matrix `P` --
are functions of the set, escape G1 the same way, and cost one matrix operation on data on disk.

THE TWO WAYS IT CAN DIE, both checked here BEFORE anything is built:

  C1  COLLAPSE. On a DENSE similarity graph the geodesic and the commute time may be monotone
      functions of the DIRECT distance, in which case the "multi-structure" object is a
      relabelling of P and buys nothing new. Measured as the within-target Spearman rho between
      each derived pair quantity and the direct P, over all 75*74/2 pairs.

  C2  REDUCTION TO CONSENSUS. Node-level graph observables (centrality, degree, mean commute
      time) may be monotone in the MEDOID CRITERION mean_j P_ij -- which is the consensus signal
      S12/S30 already measured at -0.172 A. Measured as the in-band partial Spearman
      rho(node stat, ORACLE rmsd | medoid criterion).

If C1 and C2 both fire, the affordable half of the multi-structure escape is closed and only the
unaffordable half (true barriers) remains open. Either outcome is cheap and decisive.

ORACLE: `rr` labels every correlation. Nothing here tunes a deployable parameter.
Basis: per-candidate CA point cloud, on the SHIPPED top-75 band.

USAGE   python s31/s31_B3_graph.py [--limit N] [--knn 8]
"""
from __future__ import annotations

import argparse
import io
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
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s31_B3_graph.json")
S27_CACHE = os.path.join(ROOT, "s27", "cache")


def rankz(x):
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(r.std(), 1e-12)


def spear(x, y):
    return float((rankz(x) * rankz(y)).mean())


def partial(x, y, z):
    a, b, c = rankz(x), rankz(y), rankz(z)
    rxy, rxz, ryz = float((a * b).mean()), float((a * c).mean()), float((b * c).mean())
    den = np.sqrt(max(1e-12, (1 - rxz ** 2) * (1 - ryz ** 2)))
    return float((rxy - rxz * ryz) / den)


def graph_quantities(P, knn):
    """Geodesic and commute-time matrices on a symmetric kNN graph built from P."""
    from scipy.sparse.csgraph import shortest_path
    m = len(P)
    k = min(knn, m - 1)
    A = np.zeros_like(P)
    order = np.argsort(P, axis=1)
    for i in range(m):
        for j in order[i, 1:k + 1]:
            A[i, j] = A[j, i] = P[i, j]          # edge LENGTH = direct RMSD
    G = shortest_path(A, method="D", directed=False, unweighted=False)
    # commute time from the Laplacian of the GAUSSIAN-weighted version of the same graph
    sig = float(np.median(P[P > 0])) if (P > 0).any() else 1.0
    Wt = np.where(A > 0, np.exp(-(A ** 2) / (2.0 * sig ** 2 + 1e-12)), 0.0)
    d = Wt.sum(1)
    L = np.diag(d) - Wt
    Lp = np.linalg.pinv(L)
    vol = float(d.sum())
    diag = np.diag(Lp)
    CT = vol * (diag[:, None] + diag[None, :] - 2.0 * Lp)
    w, V = np.linalg.eigh(L)
    fied = V[:, 1] if len(w) > 1 else np.zeros(m)
    gap = float(w[1]) if len(w) > 1 else float("nan")
    return G, CT, Wt, np.abs(fied), gap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--knn", type=int, default=8)
    a = ap.parse_args()

    tg = I.targets()
    if a.limit:
        tg = tg[: a.limit]

    node_names = ["medoid_crit", "commute_cent", "geo_cent", "degree", "fiedler_abs"]
    c1 = {"geodesic": [], "commute": []}
    rho = {c: [] for c in node_names}
    prho = {c: [] for c in node_names}
    dis_rho, gaps, pdbs, folds = [], [], [], []
    t0 = time.time()

    for ti, t in enumerate(tg):
        pdb, n, fold = t["pdb"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        idx = I.pool_idx(u)
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        sub = sub[sub < len(idx)]
        if len(sub) < 20:
            continue
        W = u["W"][idx][sub]
        rr = u["rr"][idx][sub]                                   # ORACLE
        with np.load(os.path.join(S27_CACHE, f"{pdb}.npz")) as z:
            DIS = np.asarray(z["DIS"], float)[idx.size and slice(None)][sub]

        P = I.pairwise_rmsd(W)
        G, CT, Wt, fied, gap = graph_quantities(P, a.knn)
        iu = np.triu_indices(len(P), 1)
        fin = np.isfinite(G[iu]) & np.isfinite(CT[iu])
        if fin.sum() < 50:
            continue
        c1["geodesic"].append(spear(G[iu][fin], P[iu][fin]))
        c1["commute"].append(spear(CT[iu][fin], P[iu][fin]))

        node = {"medoid_crit": P.mean(1), "commute_cent": np.nanmean(CT, axis=1),
                "geo_cent": np.where(np.isfinite(G), G, np.nanmax(G[np.isfinite(G)])).mean(1),
                "degree": -Wt.sum(1), "fiedler_abs": fied}
        dis_rho.append(spear(DIS, rr))
        for c in node_names:
            x = node[c]
            if not np.isfinite(x).all() or np.std(x) < 1e-12:
                rho[c].append(np.nan); prho[c].append(np.nan)
                continue
            rho[c].append(spear(x, rr))
            prho[c].append(partial(x, rr, node["medoid_crit"]) if c != "medoid_crit"
                           else partial(x, rr, DIS))
        gaps.append(gap); pdbs.append(pdb); folds.append(fold)
        if (ti + 1) % 20 == 0:
            print("  %3d/%d  %.0fs" % (ti + 1, len(tg), time.time() - t0), flush=True)

    folds = np.asarray(folds, int)
    ncomp = 0

    def verdict(v, lab):
        nonlocal ncomp
        v = np.asarray(v, float); m = np.isfinite(v)
        if m.sum() < 20:
            return dict(n=int(m.sum()), note="too few")
        o = ST.compare(-v[m], np.zeros(int(m.sum())), folds=folds[m],
                       names=list(np.array(pdbs)[m]), label=lab, seed_parts=("s31B", "graph"))
        ncomp += 1
        return dict(n=int(m.sum()), mean=float(v[m].mean()), median=float(np.median(v[m])),
                    mde=o["mde"], x_mde=float(v[m].mean() / o["mde"]) if o["mde"] > 0 else float("nan"),
                    ci95_fold=[-o["ci95_fold"][1], -o["ci95_fold"][0]],
                    folds_same_sign=o["folds_same_sign"])

    out = dict(n=len(pdbs), knn=a.knn,
               basis="SHIPPED top-75 band, per-candidate CA point cloud RMSD (ORACLE labels)",
               oracle="ORACLE / NOT DEPLOYABLE", prereg="s31/PREREG_S31_B.md",
               provenance=ST.provenance(__file__))
    out["C1_collapse"] = {k: dict(mean=float(np.mean(v)), median=float(np.median(v)),
                                  p10=float(np.percentile(v, 10)), min=float(np.min(v)),
                                  fired=bool(np.mean(v) > 0.95))
                          for k, v in c1.items() if v}
    out["spectral_gap"] = dict(mean=float(np.mean(gaps)), median=float(np.median(gaps)))
    out["DIS_inband_rho"] = verdict(dis_rho, "DIS_inband")
    out["node_rho"] = {c: verdict(rho[c], "rho_" + c) for c in node_names}
    out["node_partial_on_medoid"] = {c: verdict(prho[c], "prho_" + c) for c in node_names}
    out["comparisons_emitted"] = ncomp
    with io.open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)

    print("\nn %d   %.0fs   comparisons %d" % (len(pdbs), time.time() - t0, ncomp))
    print("\nC1  COLLAPSE -- within-target Spearman(derived pair quantity, DIRECT RMSD) over "
          "75*74/2 pairs")
    for k, v in out["C1_collapse"].items():
        print("    %-12s mean rho = %.4f  median %.4f  p10 %.4f  min %.4f   COLLAPSED(>0.95)=%s"
              % (k, v["mean"], v["median"], v["p10"], v["min"], v["fired"]))
    print("\nC2  NODE-LEVEL in-band skill, and the same PARTIALLED on the medoid criterion")
    print("    %-16s %9s %7s %6s   %9s %7s %6s" %
          ("node stat", "rho", "xMDE", "folds", "rho|medoid", "xMDE", "folds"))
    d0 = out["DIS_inband_rho"]
    if "mean" in d0:
        print("    %-16s %+9.4f %7.2f %5d/5" % ("DIS (incumbent)", d0["mean"], d0["x_mde"],
                                                d0["folds_same_sign"]))
    for c in node_names:
        r, p = out["node_rho"][c], out["node_partial_on_medoid"][c]
        if "mean" not in r:
            continue
        print("    %-16s %+9.4f %7.2f %5d/5   %+9.4f %7.2f %5d/5" % (
            c, r["mean"], r["x_mde"], r["folds_same_sign"],
            p.get("mean", float("nan")), p.get("x_mde", float("nan")),
            p.get("folds_same_sign", 0)))
    print("\n    (medoid_crit's 'partial' column is partialled on DIS, not on itself.)")
    print("\nwrote", OUT, flush=True)


if __name__ == "__main__":
    main()
