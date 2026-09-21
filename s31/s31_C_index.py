#!/usr/bin/env python
"""s31/s31_C_index.py -- lane C, question C2: candidate-index information allocation (charter §12).

The deployed 7-bit index maps basis state -> the top-128 in SCORE ORDER.  Under `R1-sel` a
re-labelling cannot raise capacity above 7 bits -- it names one of 128 things either way.  What a
re-labelling CAN change is (a) which candidate sets are *subcubes*, hence which distributions a
product-state ansatz can put its mass on, and (b) what a PARTIAL measurement localises.  Both are
measured here, over the SAME 128 candidates, so only the assignment varies (matched-space control).

Maps
  score         deployed: rank r -> binary(r)
  gray          r -> r ^ (r >> 1).  ASSERTED in code: the j-bit prefix partition is IDENTICAL to
                binary's, so Gray coding can only change Hamming geometry, never a prefix cell.
  bisect        recursive 2-medoid bisection of the 128 under pairwise Ca-RMSD; bit j = branch j
  bisect_score  bisect, with branch 0 the better mean-DIS side (coarse-to-fine)
  spectral      Fiedler-vector order of the RMSD similarity graph, then binary
  perm          a seeded random permutation -- the matched zero-information control

    python s31/s31_C_index.py run [--limit N]
    python s31/s31_C_index.py analyse
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

from s12 import instrument as I                    # noqa: E402
from s24 import stats_lib as ST                    # noqa: E402
from s31 import s31_C_cache as CA                  # noqa: E402
from s31 import s31_C_ladder as LD                 # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s31_C_index_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_C_index.json")
TOP, NBITS = 128, 7
MAPS = ("score", "gray", "bisect", "bisect_score", "spectral", "perm")
SALT = "s31Cindex"


# ------------------------------------------------------------------ index maps
def bits_of(code, nbits=NBITS):
    return np.array([[(c >> (nbits - 1 - b)) & 1 for b in range(nbits)] for c in code], np.int8)


def map_score(P, dis, rng):
    return np.arange(TOP)                       # candidate at score-rank r gets code r


def map_gray(P, dis, rng):
    r = np.arange(TOP)
    return r ^ (r >> 1)


def map_perm(P, dis, rng):
    return rng.permutation(TOP)


def _bisect(P, members, dis, by_score, depth, code, out):
    if depth == 0 or len(members) <= 1:
        for j, m in enumerate(members):
            out[m] = code * (1 << depth) + j
        return
    sub = P[np.ix_(members, members)]
    a = int(np.argmax(sub.sum(1)))              # the two most distant medoid seeds
    b = int(np.argmax(sub[a]))
    lab = (sub[a] > sub[b]).astype(int)         # 0 = nearer seed a
    # force a balanced split: assign by the signed distance difference, median cut
    d = sub[a] - sub[b]
    o = np.argsort(d, kind="stable")
    half = len(members) // 2
    lab = np.ones(len(members), int); lab[o[:half]] = 0
    g0 = [members[i] for i in range(len(members)) if lab[i] == 0]
    g1 = [members[i] for i in range(len(members)) if lab[i] == 1]
    if by_score and np.mean(dis[g0]) > np.mean(dis[g1]):
        g0, g1 = g1, g0                         # branch 0 is the better-scoring side
    _bisect(P, g0, dis, by_score, depth - 1, code * 2 + 0, out)
    _bisect(P, g1, dis, by_score, depth - 1, code * 2 + 1, out)


def map_bisect(P, dis, rng, by_score=False):
    out = np.zeros(TOP, int)
    _bisect(P, list(range(TOP)), dis, by_score, NBITS, 0, out)
    return out


def map_spectral(P, dis, rng):
    S = np.exp(-(P / max(np.median(P), 1e-9)) ** 2)
    D = np.diag(S.sum(1))
    L = D - S
    w, v = np.linalg.eigh(np.linalg.inv(np.sqrt(D)) @ L @ np.linalg.inv(np.sqrt(D)))
    f = v[:, 1]
    order = np.argsort(f)
    code = np.empty(TOP, int)
    code[order] = np.arange(TOP)
    return code


MAPFN = {"score": map_score, "gray": map_gray, "perm": map_perm,
         "bisect": lambda P, d, r: map_bisect(P, d, r, False),
         "bisect_score": lambda P, d, r: map_bisect(P, d, r, True),
         "spectral": map_spectral}


# ------------------------------------------------------------------ product-state ceiling
def product_p(theta):
    q = 1.0 / (1.0 + np.exp(-theta))                       # (NBITS,)
    p = np.ones(1)
    for b in range(len(q)):
        p = np.concatenate([p * (1.0 - q[b]), p * q[b]])   # bit b is the b-th MSB
    return p


def oracle_product_ceiling(Wt_by_code, natp, n, rng, starts=4, iters=400):
    """ORACLE. Best RMSD reachable by a PRODUCT distribution over the 7 index bits, read out as the
    p-weighted average.  `Wt_by_code[c]` is the flattened candidate carrying code c."""
    from scipy.optimize import minimize

    def f(theta):
        p = product_p(theta)
        d = p @ Wt_by_code - natp
        return float(d @ d) / n

    best = (float("inf"), None)
    for s in range(starts):
        th0 = np.zeros(NBITS) if s == 0 else rng.normal(0, 1.5, NBITS)
        r = minimize(f, th0, method="L-BFGS-B", options=dict(maxiter=iters))
        if r.fun < best[0]:
            best = (float(r.fun), r.x.copy())
    p = product_p(best[1])
    return best, p


# ------------------------------------------------------------------ one target
def row_for(pdb):
    from s27 import s28_A_amp as A28
    z = CA.load(pdb)
    W = np.asarray(z["W"], float)
    order = np.asarray(z["order"], int)
    dis = np.asarray(z["dis"], float)
    rr = np.asarray(z["rr"], float)                                   # ORACLE
    nat = np.asarray(z["nat"], float)                                 # ORACLE
    top75 = np.asarray(z["top75"], int)
    n = int(z["n"])
    frame = A28.Frame(W, np.sort(top75))
    Wf = frame.Wf
    natp = I.superpose_batch(nat[None], frame.ref)[0].ravel()         # ORACLE, scoring only
    idx = order[:TOP]
    Wt = Wf[idx]                                                      # (128, 3n) in score-rank order
    P = I.pairwise_rmsd(W[idx])
    d = dis[idx]
    rrt = rr[idx]                                                     # ORACLE
    rng = np.random.default_rng(abs(hash((SALT, pdb))) % (2 ** 32))
    iu, ju = np.triu_indices(TOP, 1)
    struct = P[iu, ju]
    good = (rrt <= np.percentile(rrt, 10)).astype(float)              # ORACLE, for the bit MI only

    row = dict(pdb=pdb, n=n, fold=int(z["fold"]), fail18=bool(pdb in I.FAIL18), maps={})
    for name in MAPS:
        code = np.asarray(MAPFN[name](P, d, rng), int)
        assert sorted(code.tolist()) == list(range(TOP)), "%s: %s is not a bijection" % (pdb, name)
        B = bits_of(code)                                             # (128, 7) per candidate
        ham = (B[iu] != B[ju]).sum(1)
        rho = float(ST.__dict__.get("_spearman", lambda a, b: np.nan)(ham, struct)) \
            if "_spearman" in ST.__dict__ else float(
                np.corrcoef(np.argsort(np.argsort(ham)), np.argsort(np.argsort(struct)))[0, 1])

        # the candidate carrying each code
        by_code = np.empty(TOP, int)
        by_code[code] = np.arange(TOP)
        Wc = Wt[by_code]
        rc = rrt[by_code]; dc = d[by_code]

        # ---- M1 product-state ORACLE ceiling
        (fbest, th), p = oracle_product_ceiling(Wc, natp, n, rng)
        X = (p @ Wc).reshape(n, 3)
        m1 = float(I.ca_rmsd(X, nat))                                 # ORACLE

        # ---- M2 partial measurement: j-bit prefix cells
        cells = {}
        for j in range(1, NBITS):
            size = 1 << (NBITS - j)
            best_o, sc_rule, cons_rule = [], [], []
            means = []
            for cidx in range(1 << j):
                sel = np.arange(cidx * size, (cidx + 1) * size)
                Xc = Wc[sel].mean(0).reshape(n, 3)
                means.append((float(I.ca_rmsd(Xc, nat)),               # ORACLE score of the cell
                              float(dc[sel].mean()),
                              float(P[np.ix_(by_code[sel], np.arange(TOP))].mean())))
            arr = np.array(means)
            best_o = float(arr[:, 0].min())                            # ORACLE-best cell
            sc_rule = float(arr[int(np.argmin(arr[:, 1])), 0])         # native-free: best mean DIS
            cons_rule = float(arr[int(np.argmin(arr[:, 2])), 0])       # native-free: most typical
            cells["j%d" % j] = dict(oracle_best=best_o, by_score=sc_rule, by_consensus=cons_rule,
                                    mean_cell=float(arr[:, 0].mean()))

        # ---- M3 bit value: MI between each bit and "in the ORACLE-best decile"
        mi = []
        for b in range(NBITS):
            x = B[:, b]
            m = 0.0
            for xv in (0, 1):
                for yv in (0, 1):
                    pxy = float(((x == xv) & (good == yv)).mean())
                    px = float((x == xv).mean()); py = float((good == yv).mean())
                    if pxy > 0 and px > 0 and py > 0:
                        m += pxy * math.log2(pxy / (px * py))
            mi.append(float(m))

        row["maps"][name] = dict(rho_hamming_struct=rho, product_ceiling=m1,
                                 cells=cells, mi_bits=mi, mi_total=float(sum(mi)))

    # the analytic Gray assertion: the j-bit prefix PARTITION is identical to binary's
    cs = np.asarray(MAPFN["score"](P, d, rng), int)
    cg = np.asarray(MAPFN["gray"](P, d, rng), int)
    same = all(
        {frozenset(np.where(cs >> (NBITS - j) == v)[0].tolist()) for v in range(1 << j)} ==
        {frozenset(np.where(cg >> (NBITS - j) == v)[0].tolist()) for v in range(1 << j)}
        for j in range(1, NBITS))
    row["gray_prefix_partition_identical"] = bool(same)
    row["anchors"] = dict(uniform128=float(I.ca_rmsd(Wt.mean(0).reshape(n, 3), nat)),
                          argmin128=float(rrt.min()))
    return row


def phase_run(limit=0):
    pdbs = CA.all_pdbs()
    if limit:
        pdbs = pdbs[:limit]
    done = set()
    if os.path.exists(ROWS):
        done = {json.loads(l)["pdb"] for l in open(ROWS)}
    t0 = time.time()
    with open(ROWS, "a") as fh:
        for i, p in enumerate(pdbs):
            if p in done:
                continue
            fh.write(json.dumps(row_for(p)) + "\n"); fh.flush()
            if i % 10 == 0 or i == len(pdbs) - 1:
                print("[%3d/%3d] %s %.1fs" % (i + 1, len(pdbs), p, time.time() - t0), flush=True)
    print("rows done in %.1fs" % (time.time() - t0))


def phase_analyse():
    rows = {json.loads(l)["pdb"]: json.loads(l) for l in open(ROWS)}
    pdbs = [p for p in CA.all_pdbs() if p in rows]
    rows = [rows[p] for p in pdbs]
    folds = np.array([r["fold"] for r in rows])
    n = len(rows)
    out = {"n": n, "basis": "CA point cloud", "maps": {},
           "gray_prefix_partition_identical_on_all": bool(all(r["gray_prefix_partition_identical"]
                                                              for r in rows)),
           "anchors": {"uniform128": float(np.mean([r["anchors"]["uniform128"] for r in rows])),
                       "argmin128_ORACLE": float(np.mean([r["anchors"]["argmin128"] for r in rows])),
                       "production_top75_cloud": 3.0483},
           "provenance": ST.provenance(__file__)}

    def g(mp, key):
        return np.array([r["maps"][mp][key] for r in rows])

    for mp in MAPS:
        d = {"rho_hamming_struct": float(g(mp, "rho_hamming_struct").mean()),
             "product_ceiling_ORACLE": float(g(mp, "product_ceiling").mean()),
             "mi_total_ORACLE": float(g(mp, "mi_total").mean()),
             "mi_bits_ORACLE": [float(np.mean([r["maps"][mp]["mi_bits"][b] for r in rows]))
                                for b in range(NBITS)],
             "cells": {}}
        for j in range(1, NBITS):
            k = "j%d" % j
            d["cells"][k] = {q: float(np.mean([r["maps"][mp]["cells"][k][q] for r in rows]))
                             for q in ("oracle_best", "by_score", "by_consensus", "mean_cell")}
        out["maps"][mp] = d

    # ---- F-C2c pre-check, F-C2a ceiling, F-C2b deployable -----------------------------------
    base_rho = out["maps"]["score"]["rho_hamming_struct"]
    out["F_C2c"] = {"score_rho": base_rho,
                    "bisect_rho": out["maps"]["bisect"]["rho_hamming_struct"],
                    "spectral_rho": out["maps"]["spectral"]["rho_hamming_struct"],
                    "perm_rho": out["maps"]["perm"]["rho_hamming_struct"],
                    "bar": 0.15,
                    "fires": bool(min(out["maps"]["bisect"]["rho_hamming_struct"],
                                      out["maps"]["spectral"]["rho_hamming_struct"])
                                  - base_rho >= 0.15)}

    base_c = g("score", "product_ceiling")
    out["F_C2a"] = {}
    for mp in MAPS:
        if mp == "score":
            continue
        c = ST.compare(g(mp, "product_ceiling"), base_c, folds=folds, names=pdbs,
                       label="prodstate_%s_vs_score" % mp, seed_parts=("s31C", "index"))
        out["F_C2a"][mp] = {k: c[k] for k in ("effect", "se", "effect_over_mde", "ci95_fold",
                                              "folds_same_sign", "n_better", "n_worse",
                                              "mean_a", "mean_b")}
    best = min(out["F_C2a"], key=lambda m: out["F_C2a"][m]["effect"])
    out["F_C2a_verdict"] = {"best_map": best, "effect": out["F_C2a"][best]["effect"],
                            "bar": -0.10,
                            "fires": bool(out["F_C2a"][best]["effect"] <= -0.10 and
                                          out["F_C2a"][best]["ci95_fold"][1] < 0)}

    prod = 3.0483
    cands = []
    for mp in MAPS:
        for j in range(1, NBITS):
            for rule in ("by_score", "by_consensus"):
                cands.append((mp, j, rule, out["maps"][mp]["cells"]["j%d" % j][rule]))
    cands.sort(key=lambda t: t[3])
    out["F_C2b"] = {"best_triple": {"map": cands[0][0], "j": cands[0][1], "rule": cands[0][2],
                                    "mean": cands[0][3], "vs_production": cands[0][3] - prod},
                    "bar": -0.10,
                    "fires": bool(cands[0][3] - prod <= -0.10),
                    "top5": [{"map": a, "j": b, "rule": c, "mean": v} for a, b, c, v in cands[:5]]}
    # the grid is (6 maps x 6 j x 2 rules) = 72 cells: price the per-target maximum (rule 11)
    M = np.array([[r["maps"][mp]["cells"]["j%d" % j][rule] for mp in MAPS
                   for j in range(1, NBITS) for rule in ("by_score", "by_consensus")]
                  for r in rows])
    out["F_C2b"]["grid_shape"] = list(M.shape)
    out["F_C2b"]["split_half_transfer"] = ST.split_half_transfer(M)
    out["F_C2b"]["best_of_k_within"] = ST.best_of_k_within(M)

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    print("\nC2 -- CANDIDATE-INDEX INFORMATION ALLOCATION.  Same 128 candidates, only the "
          "assignment varies.  CA point cloud, n=%d\n" % n)
    print("  anchors: production top-75 uniform 3.0483 | uniform over the 128 %.4f | "
          "ORACLE argmin over the 128 %.4f" %
          (out["anchors"]["uniform128"], out["anchors"]["argmin128_ORACLE"]))
    print("\n  %-14s %10s %14s %10s" % ("map", "rho(H,RMSD)", "prod-state ORACLE", "MI bits"))
    for mp in MAPS:
        d = out["maps"][mp]
        print("  %-14s %10.4f %14.4f %10.4f" %
              (mp, d["rho_hamming_struct"], d["product_ceiling_ORACLE"], d["mi_total_ORACLE"]))
    print("\n  F-C2c pre-check (bar +0.15 over score's rho): %s" %
          ("FIRES" if out["F_C2c"]["fires"] else "REFUTED"))
    print("  F-C2a product-state ceiling, best map %s at %+.4f A vs score  ->  %s" %
          (out["F_C2a_verdict"]["best_map"], out["F_C2a_verdict"]["effect"],
           "FIRES" if out["F_C2a_verdict"]["fires"] else "REFUTED"))
    print("\n  PARTIAL MEASUREMENT (j bits measured -> a cell of 2^(7-j), averaged):")
    print("  %-14s %4s %10s %10s %10s %10s" %
          ("map", "j", "ORACLE", "by_score", "by_cons", "cell mean"))
    for mp in MAPS:
        for j in (1, 3, 5):
            c = out["maps"][mp]["cells"]["j%d" % j]
            print("  %-14s %4d %10.4f %10.4f %10.4f %10.4f" %
                  (mp, j, c["oracle_best"], c["by_score"], c["by_consensus"], c["mean_cell"]))
    b = out["F_C2b"]["best_triple"]
    print("\n  F-C2b best native-free triple: %s j=%d %s -> %.4f (%+.4f vs production)  ->  %s" %
          (b["map"], b["j"], b["rule"], b["mean"], b["vs_production"],
           "FIRES" if out["F_C2b"]["fires"] else "REFUTED"))
    print("  grid %s  split-half transfer %s" %
          (out["F_C2b"]["grid_shape"], ST.fmt(out["F_C2b"]["split_half_transfer"])
           if hasattr(ST, "fmt") else out["F_C2b"]["split_half_transfer"]))
    print("  Gray prefix partition identical to binary on all targets: %s" %
          out["gray_prefix_partition_identical_on_all"])
    print("\n  wrote %s" % OUT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=["run", "analyse"])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.phase == "run":
        phase_run(a.limit)
    else:
        phase_analyse()


if __name__ == "__main__":
    main()
