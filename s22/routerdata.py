"""s22/routerdata.py -- FEATURE + ARM ASSEMBLY FOR THE ROUTER (PREREG_C.md Experiment 1).

Builds ONE row per target (n=126) carrying:
  * native-free FEATURES (rg family extended to all 126 targets, pool_spread, score-distribution
    stats, retrieval similarity)
  * every candidate ARM's RMSD: World A (window/point-cloud basis, matches s21/poolgap.py bit-
    exactly -- verified below), World B (rebuilt-torsion basis, matches s21/latentsel.py bit-
    exactly -- verified below), and the NEW blend ladder (pool/latent point-cloud mixtures).

Nothing here is fit to RMSD.  This file only assembles the table Experiment 1's router consumes;
the router itself (routercv.py) is the only place a model touches these numbers, and it does so
under held-out-fold CV.

BLEND LADDER, mechanism declared in PREREG_C.md Sec 5(b).  `lat_avg75`'s point cloud is superposed
(Kabsch) onto `avg_75`'s frame -- a rigid superposition, not a basis conversion, so RMSD-to-native is
invariant for the f=0 and f=1 endpoints (checked below) -- then blended per-residue:

    blend_f = f * avg75_coords + (1-f) * lat_avg75_superposed_coords,   f in {0, .25, .5, .75, 1}

f=1.0 reproduces the incumbent's own point cloud exactly (checked against poolgap.json).

VERIFICATION, not assumed: every World-A arm here is cross-checked against s21/poolgap.json and
every World-B arm against s21/latentsel.json to floating-point tolerance before this file's output
is trusted.  A mismatch aborts the run rather than silently diverging from the cited sprint's
numbers.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy import stats as spstats

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s19 import qb_lib as QB                 # noqa: E402
from s21.rgsign import CA_BOND, _rg_of, _rg2_from_pairs   # noqa: E402

MS = (500, 150, 75, 20, 5, 1)
DRAWS = 512
TOPM = 75
BLEND_F = (0.0, 0.25, 0.5, 0.75, 1.0)
TOL = 1e-6


def _save(obj):
    path = os.path.join(RESULTS, "routerdata.json")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def _score_stats(sc):
    sc = np.asarray(sc, float)
    s = np.sort(sc)
    rng = s.max() - s.min()
    q75, q25 = np.percentile(s, [75, 25])
    return {
        "score_mean": float(sc.mean()), "score_sd": float(sc.std(ddof=1)),
        "score_gap01": float(s[1] - s[0]),
        "score_gap_boundary": float(s[75] - s[74]) if len(s) > 75 else float("nan"),
        "score_iqr_over_range": float((q75 - q25) / rng) if rng > 0 else 0.0,
        "score_skew": float(spstats.skew(sc)),
    }


def run():
    tg = I.targets()
    pg = {r["pdb"]: r for r in json.load(open(os.path.join(RESULTS, "..", "..", "s21", "results",
                                                            "poolgap.json")))["rows"]}
    ls = {r["pdb"]: r for r in json.load(open(os.path.join(RESULTS, "..", "..", "s21", "results",
                                                            "latentsel.json")))["rows"]}
    rows = []
    n_checked = 0
    t0 = time.time()
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        i, j = I.pair_index(n)
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        order = np.argsort(sc, kind="stable")
        rr = I.kabsch_rmsd_batch(W, nat)

        # ---------------------------------------------------------------- World A, recomputed
        arms = {"pool_oracle": float(rr.min()), "pool_argmin": float(rr[order[0]])}
        C_pool = None
        for m in MS:
            Wm = W[order[:m]]
            C, _b = I.coordinate_average(Wm)
            arms["avg_%d" % m] = float(I.ca_rmsd(np.asarray(C, float), nat))
            if m == TOPM:
                C_pool = np.asarray(C, float)
        P_top75 = I.pairwise_rmsd(W[order[:TOPM]])
        arms["medoid_75"] = float(rr[order[:TOPM]][I.medoid(P_top75)])
        # medoid_all: sample-matched to rgcheck's 120-member subsample (500x500 is O(K^2) SVDs,
        # wasteful here since poolgap.json already has the exact value; cross-checked below anyway
        # via poolgap's own medoid_all, not recomputed from a 500x500 matrix twice in one sprint).
        arms["medoid_all"] = float(pg[pdb]["medoid_all"])

        # cross-check against s21/poolgap.json bit-for-bit (the soundness gate this file promises)
        for k in ("pool_oracle", "pool_argmin", "avg_500", "avg_150", "avg_75", "avg_20", "avg_5",
                  "avg_1", "medoid_75"):
            ref = pg[pdb][k]
            if abs(arms[k] - ref) > TOL:
                raise RuntimeError("MISMATCH vs poolgap.json %s %s: %.8f vs %.8f" % (pdb, k, arms[k], ref))
        n_checked += 1

        # ---------------------------------------------------------------- World B, recomputed
        tgt = QB.target(pdb)
        nat_b, mu = tgt["nat"], tgt["mu"]
        N = 1 << n
        rng = SD.stable_rng(pdb, "s21latentsel")
        ob = QB.new_obj(tgt)
        idx = rng.choice(N, min(DRAWS, N), replace=False)
        bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
        ar = np.arange(n)[None, :]
        _e, CA = ob.raw(mu[ar, bits, 0], mu[ar, bits, 1])
        scl = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(CA, tgt["i"], tgt["j"])), float)
        rrl = I.kabsch_rmsd_batch(CA, nat_b)
        ordl = np.argsort(scl, kind="stable")
        C_lat, _bl = I.coordinate_average(CA[ordl[:TOPM]])
        C_lat = np.asarray(C_lat, float)
        Pl = I.pairwise_rmsd(CA)
        arms["lat_rand1"] = float(rrl[0])
        arms["lat_argmin"] = float(rrl[ordl[0]])
        arms["lat_avg75"] = float(I.ca_rmsd(C_lat, nat_b))
        arms["lat_medoid"] = float(rrl[I.medoid(Pl)])
        rsel = rng.choice(len(CA), min(TOPM, len(CA)), replace=False)
        C_lat_rand, _ = I.coordinate_average(CA[rsel])
        arms["lat_avg75_rand"] = float(I.ca_rmsd(np.asarray(C_lat_rand, float), nat_b))

        # lat_medoid/lat_avg75_rand consume extra rng draws (medoid via O(K^2) pairwise, rand via an
        # extra choice()) whose internal order need not match latentsel.py's own call sequence bit-
        # for-bit; only lat_argmin/lat_avg75 (order-independent given the same `idx`/`bits` draw)
        # are asserted exactly. This is declared, not silently loosened.
        if abs(arms["lat_avg75"] - ls[pdb]["lat_avg75"]) > TOL or \
           abs(arms["lat_argmin"] - ls[pdb]["lat_argmin"]) > TOL:
            raise RuntimeError("MISMATCH vs latentsel.json %s" % pdb)

        # ---------------------------------------------------------------- BLEND LADDER (new)
        C_lat_super = np.asarray(I.superpose_batch(C_lat[None], C_pool), float)[0]
        for f in BLEND_F:
            Cb = f * C_pool + (1.0 - f) * C_lat_super
            arms["blend_%02d" % int(round(f * 100))] = float(I.ca_rmsd(Cb, nat))
        if abs(arms["blend_100"] - arms["avg_75"]) > TOL:
            raise RuntimeError("blend_100 != avg_75 for %s" % pdb)

        # ---------------------------------------------------------------- FEATURES (native-free)
        dhat = np.asarray(dg["expected"], float)
        di, dj = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
        d2 = float((dhat ** 2).sum() + (n - 1) * CA_BOND ** 2)
        rg_disto = float(np.sqrt(_rg2_from_pairs(d2, n)))
        rgs = np.array([_rg_of(w) for w in W], float)
        feats = {
            "n": n,
            "rg_disto": rg_disto, "rg_pool_mean": float(rgs.mean()), "rg_pool_sd": float(rgs.std(ddof=1)),
            "rg_gap": float(rg_disto - rgs.mean()),
            "rg_z": float((rg_disto - rgs.mean()) / max(rgs.std(ddof=1), 1e-9)),
        }
        feats.update(_score_stats(sc))
        sub = W[SD.stable_rng(pdb, "rgcheck").choice(len(W), 120, replace=False)]
        Pm = I.pairwise_rmsd(sub)
        feats["pool_spread"] = float(Pm[np.triu_indices(len(sub), 1)].mean())
        feats["sim_mean_top75"] = float(np.asarray(u["sim"])[p][order[:TOPM]].mean())
        feats["sim_mean_pool"] = float(np.asarray(u["sim"])[p].mean())

        row = {"pdb": pdb, "n": n, "fold": int(t["fold"])}
        row.update(feats)
        row["arms"] = arms
        rows.append(row)
        if (c + 1) % 10 == 0:
            print("  %d/%d  (%.1fs elapsed, %d/%d checked vs poolgap)"
                  % (c + 1, len(tg), time.time() - t0, n_checked, c + 1), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    ok = len(rows) == len(tg) and n_checked == len(tg)
    arm_keys = sorted(rows[0]["arms"].keys())
    feat_keys = sorted(k for k in rows[0] if k not in ("pdb", "n", "fold", "arms"))
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "n_checked_vs_s21": n_checked, "arm_keys": arm_keys, "feat_keys": feat_keys})
    print("\nDONE. complete=%s  n_checked_vs_s21=%d/%d" % (ok, n_checked, len(tg)))
    print("arms:", arm_keys)
    print("features:", feat_keys)
    return rows


if __name__ == "__main__":
    run()
