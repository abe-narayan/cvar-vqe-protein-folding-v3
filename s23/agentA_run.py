"""s23/agentA_run.py -- WORKSTREAM A: ENSEMBLE AGGREGATION (Sprint 23).

Implements the three experiments pre-registered in `s23/PREREG_A.md`, exactly as declared there:
  1. CLUSTER-RESTRICTED AVERAGING (priority 1)
  2. ENSEMBLE-WIDTH SWEEP under plain vs cluster-restricted aggregation (priority 2)
  3. WEIGHTED AVERAGING, nested-CV fitted (priority 3)

One pass per target computes the full K=500 pairwise-RMSD matrix ONCE and reuses slices of it for
every arm in every experiment, rather than recomputing it three times. Point-cloud basis throughout;
the incumbent `avg_75` is recomputed here (not read from `s21/poolgap.json`) so this file's own
numbers are self-consistent, and it is asserted to match the pinned constant at report time.

Atomic writes (`os.replace` from a `.tmp`), completion flags require the FULL key set (not a row
count), config-derived filename is a single artefact `s23/results/agentA.json` since all three
experiments share one per-target row (no separate configs to disambiguate in the filename).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402

# ------------------------------------------------------------------------------------- config
M_LADDER = (20, 30, 50, 75, 100, 150, 250, 500)
K_LIST = (2, 3, 4, 5)
M_PRIMARY = 75
M_ROBUST = (150, 250)
K_ROBUST = 3
K_DEFAULT_FOR_WIDTH_SWEEP = 3

SCORE_C_GRID = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 1e6)
RANK_P_GRID = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0)
DIST_C_GRID = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 1e6)

FAIL_THRESH = 5.0  # A, pre-registered absolute failure bar


def _save(obj, name):
    path = os.path.join(RES, name)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)


def _vbond(C):
    C = np.asarray(C, float)
    return float(np.linalg.norm(C[1:] - C[:-1], axis=1).mean())


# ------------------------------------------------------------------------------------- clustering
def _cluster_labels(Pm, k):
    Dc = squareform(np.asarray(Pm, float), checks=False)
    Z = linkage(Dc, method="average")
    labels = fcluster(Z, k, criterion="maxclust")
    return labels


def _cluster_stats(Wm, scm, Pm, nat, labels):
    """Per-cluster: coordinate average, RMSD-to-native, size, mean score."""
    ks = sorted(set(labels.tolist()))
    out = []
    for c in ks:
        members = np.where(labels == c)[0]
        if len(members) == 1:
            Cavg = Wm[members[0]]
        else:
            Cavg, _b = I.coordinate_average(Wm[members], Pm[np.ix_(members, members)])
        Cavg = np.asarray(Cavg, float)
        out.append(dict(members=members, rmsd=I.ca_rmsd(Cavg, nat), size=len(members),
                         mean_score=float(scm[members].mean()), C=Cavg))
    return out


def _pick(clusters, rule, rng):
    n = len(clusters)
    sizes = np.array([c["size"] for c in clusters])
    scores = np.array([c["mean_score"] for c in clusters])
    rmsds = np.array([c["rmsd"] for c in clusters])
    if rule == "size":
        return int(np.argmax(sizes))
    if rule == "score":
        return int(np.argmin(scores))
    if rule == "combined":
        size_rank = np.argsort(np.argsort(-sizes))
        score_rank = np.argsort(np.argsort(scores))
        return int(np.argmin(size_rank + score_rank))
    if rule == "oracle":
        return int(np.argmin(rmsds))
    if rule == "random":
        return int(rng.integers(0, n))
    raise ValueError(rule)


def _random_same_size_partition(Wm, scm, Pm, nat, sizes, rng):
    """Shuffle members into groups of the given sizes (matching a real clustering's size
    distribution), then apply the 'combined' pick logic among the resulting groups."""
    m = Wm.shape[0]
    perm = rng.permutation(m)
    groups = []
    i0 = 0
    for s in sizes:
        groups.append(perm[i0:i0 + s])
        i0 += s
    clusters = []
    for members in groups:
        if len(members) == 0:
            continue
        if len(members) == 1:
            Cavg = Wm[members[0]]
        else:
            Cavg, _b = I.coordinate_average(Wm[members], Pm[np.ix_(members, members)])
        Cavg = np.asarray(Cavg, float)
        clusters.append(dict(members=members, rmsd=I.ca_rmsd(Cavg, nat), size=len(members),
                              mean_score=float(scm[members].mean()), C=Cavg))
    pick = _pick(clusters, "combined", rng)
    return clusters[pick]["rmsd"]


# ------------------------------------------------------------------------------------- weighting
def _weighted_avg(Wm, weights, frame):
    S = I.superpose_batch(Wm, frame)
    w = np.asarray(weights, float)
    w = w / w.sum()
    return (S * w[:, None, None]).sum(0)


def _score_weight(scm, c):
    sd = scm.std()
    tau = c * sd if sd > 0 else 1.0
    z = -(scm - scm.min()) / max(tau, 1e-12)
    z = z - z.max()
    return np.exp(z)


def _rank_weight(m, p):
    rank = np.arange(m)  # 0 = best
    return (m - rank) ** p if p > 0 else np.ones(m)


def _dist_weight(d, c):
    sd = d.std()
    sigma = c * sd if sd > 0 else 1.0
    z = -(d - d.min()) / max(sigma, 1e-12)
    z = z - z.max()
    return np.exp(z)


# ------------------------------------------------------------------------------------- main loop
def run():
    tg = I.targets()
    rows = []
    t_start = time.time()
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        rng = SD.stable_rng("s23agentA", pdb)
        u = I.load_univ(pdb)
        idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        order = np.argsort(sc, kind="stable")
        P_full = I.pairwise_rmsd(W[order])   # reorder ONCE so P_full[:m,:m] is always "top-m"
        Wsorted = W[order]
        scsorted = sc[order]

        row = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"])}

        # ---- incumbent, self-consistency ----
        avg75_C, med75_local = I.coordinate_average(Wsorted[:75], P_full[:75, :75])
        avg75_C = np.asarray(avg75_C, float)
        row["avg_75"] = float(I.ca_rmsd(avg75_C, nat))
        row["vb_avg_75"] = _vbond(avg75_C)
        row["medoid_75"] = float(I.ca_rmsd(Wsorted[med75_local], nat))

        # =========================================================== EXPERIMENT 1
        exp1 = {}
        for m in (M_PRIMARY,) + M_ROBUST:
            Wm = Wsorted[:m]; scm = scsorted[:m]; Pm = P_full[:m, :m]
            ks = K_LIST if m == M_PRIMARY else (K_ROBUST,)
            for k in ks:
                labels = _cluster_labels(Pm, k)
                clusters = _cluster_stats(Wm, scm, Pm, nat, labels)
                cell = {"n_clusters": len(clusters),
                        "sizes": sorted([cl["size"] for cl in clusters], reverse=True)}
                for rule in ("size", "score", "combined", "oracle", "random"):
                    pick = _pick(clusters, rule, rng)
                    cell[rule] = float(clusters[pick]["rmsd"])
                    cell["vb_" + rule] = _vbond(clusters[pick]["C"])
                if k == K_ROBUST and m in (M_PRIMARY,):
                    sizes = [cl["size"] for cl in clusters]
                    cell["random_same_size_partition"] = float(
                        _random_same_size_partition(Wm, scm, Pm, nat, sizes, rng))
                exp1["m%d_k%d" % (m, k)] = cell
        row["exp1"] = exp1

        # =========================================================== EXPERIMENT 2
        exp2 = {}
        for m in M_LADDER:
            mm = min(m, Wsorted.shape[0])
            Wm = Wsorted[:mm]; scm = scsorted[:mm]; Pm = P_full[:mm, :mm]
            Cavg, _b = I.coordinate_average(Wm, Pm)
            Cavg = np.asarray(Cavg, float)
            r_plain = float(I.ca_rmsd(Cavg, nat))
            labels = _cluster_labels(Pm, min(K_DEFAULT_FOR_WIDTH_SWEEP, mm))
            clusters = _cluster_stats(Wm, scm, Pm, nat, labels)
            pick = _pick(clusters, "combined", rng)
            r_clust = float(clusters[pick]["rmsd"])
            exp2["m%d" % m] = {"avg_m": r_plain, "clust_m": r_clust,
                                "vb_avg_m": _vbond(Cavg), "vb_clust_m": _vbond(clusters[pick]["C"])}
        row["exp2"] = exp2

        # =========================================================== EXPERIMENT 3
        m75 = 75
        Wm = Wsorted[:m75]; scm = scsorted[:m75]; Pm = P_full[:m75, :m75]
        frame = Wm[med75_local]
        d_med = Pm[med75_local]  # distance of each member to the medoid
        exp3 = {}
        for c in SCORE_C_GRID:
            w = _score_weight(scm, c)
            C = _weighted_avg(Wm, w, frame)
            exp3["score_c%s" % c] = float(I.ca_rmsd(np.asarray(C, float), nat))
        for p in RANK_P_GRID:
            w = _rank_weight(m75, p)
            C = _weighted_avg(Wm, w, frame)
            exp3["rank_p%s" % p] = float(I.ca_rmsd(np.asarray(C, float), nat))
        for c in DIST_C_GRID:
            w = _dist_weight(d_med, c)
            C = _weighted_avg(Wm, w, frame)
            exp3["dist_c%s" % c] = float(I.ca_rmsd(np.asarray(C, float), nat))
        row["exp3"] = exp3

        rows.append(row)
        if (c + 1) % 10 == 0 or (c + 1) == len(tg):
            el = time.time() - t_start
            print("  %d/%d  (%.1fs, %.2fs/target)" % (c + 1, len(tg), el, el / (c + 1)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)}, "agentA.json")

    need_top = ["avg_75", "medoid_75", "exp1", "exp2", "exp3"]
    ok = len(rows) == len(tg) and all(all(k in r for k in need_top) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "m_ladder": M_LADDER, "k_list": K_LIST, "m_primary": M_PRIMARY, "m_robust": M_ROBUST,
           "score_c_grid": SCORE_C_GRID, "rank_p_grid": RANK_P_GRID, "dist_c_grid": DIST_C_GRID,
           "fail_thresh": FAIL_THRESH}, "agentA.json")
    print("COMPLETE" if ok else "INCOMPLETE", flush=True)
    return rows


if __name__ == "__main__":
    run()
