#!/usr/bin/env python
"""s26/w_amberprior.py -- AMBER as a distribution inside the prior: the energy-weighted pool
histogram as the mixture partner of the shipped posterior.  Lane W (lane P's orphaned H_C3a),
Sprint 26.  Pre-registration: s26/PREREG_amber_prior_partner.md.

Stage 1 `clouds` (native-free): for every (lam, beta) cell the mixed posterior
P = (1 - lam) P_shipped + lam H_beta, H_beta the K = 500 pool's per-pair 17-bin histogram
weighted by exp(-beta * zrank(E_AMBER)) (the cached ff14SB/GBn2 single points), through a genuine
core.predict.Distogram risk table, the shipped score, the top-75 and the medoid-frame average;
the cloud and the top-75 stored, no projection.  The rank-permuted-AMBER control (energies
permuted across the 500 candidates, 4 seeded draws) at every beta > 0 cell.  lam = 0 asserted
bit-exact against the shipped risk table; beta = 0 asserted equal to p_ladder's pool_histogram.

Stage 2 `endpoint` (GATED): cloud RMSD per cell; (lam*, beta*) chosen leave-fold-out on the point
cloud; the chosen cell, the beta = 0 chosen cell, the 4 permuted draws at the chosen cell and
the shipped cloud PROJECTED (7 projections per target, checkpointed); ST.compare on the built
chain (PRIMARY) and the cloud.

    python s26/jobrun.py --agent W --tag CPU --name w_amberprior_probe --est-ram 0.5 -- python s26/w_amberprior.py clouds --probe 1A13
    python s26/jobrun.py --agent W --tag CPU --name w_amberprior_clouds --est-ram 0.5 -- python s26/w_amberprior.py clouds
    python s26/jobrun.py --agent W --tag CPU --name w_amberprior_endpoint --est-ram 0.5 -- python s26/w_amberprior.py endpoint
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402
from scipy.stats import rankdata                     # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import predict as PR                       # noqa: E402
import w_selfcopy as W                               # noqa: E402
import p_ladder as PL                                # noqa: E402

K, TOPM = W.K, W.TOPM
LAMS = (0.0, 0.05, 0.1, 0.2, 0.35, 0.5)
BETAS = (0.0, 0.5, 1.0, 2.0, 4.0)
CELLS = [(l, b) for l in LAMS for b in BETAS]          # grid order: lam outer, beta inner
N_PERM = 4
CACHE_AMB = os.path.join(ROOT, "s24", "cache_amber")


def zrank(x):
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / r.std()


def amber_energies(pdb, u):
    z = np.load(os.path.join(CACHE_AMB, pdb + ".npz"), allow_pickle=True)
    idx = np.asarray(z["universe_idx"], int)
    assert np.array_equal(idx, I.pool_idx(u)), "cache_amber pool identity != I.pool_idx"
    return np.asarray(z["e_amber"], float)


def weighted_pool_histogram(D_pool, w, eps=1e-3):
    """(npairs, 17): the pool's per-pair histogram over the shipped bins, member weights w
    (w == 1 reproduces p_ladder.pool_histogram exactly)."""
    b = np.digitize(np.asarray(D_pool, float), PR.BIN_EDGES)                     # (K, npairs)
    w = np.asarray(w, float)
    H = np.stack([((b == c) * w[:, None]).sum(0) / w.sum() for c in range(PR.NBINS)], 1) + eps
    return H / H.sum(1, keepdims=True)


def cell_cloud(P, seq, i, j, Wp, fold):
    rt = PL.risk_table(seq, P, i, j)
    e = W.emit(Wp, rt, seq, fold, project=False)
    return e["cloud"], e["top"], rt


def clouds_target(t, verbose=True):
    from s15 import seed as SD
    pdb, seq, n, fold = t["pdb"], t["seq"], t["n"], t["fold"]
    u = W.load_blind(pdb)
    idx = I.pool_idx(u); Wp = np.asarray(u["W"][idx], float)
    dg = I.distogram(pdb)
    P0 = np.asarray(dg["prob"], float); i, j = np.asarray(dg["i"]), np.asarray(dg["j"])
    Dp = I.pair_dists(Wp, i, j)
    E = amber_energies(pdb, u); zE = zrank(E)
    H0 = PL.pool_histogram(Dp)
    row = {"pdb": pdb, "n": n, "fold": fold, "cells": {}, "perm": {}, "gate": None}
    C0, top0, rt0 = cell_cloud(P0, seq, i, j, Wp, fold)
    row["gate"] = W.production_gate(pdb, {"top": top0, "cloud": C0})
    assert np.array_equal(np.asarray(rt0["risk"]), np.asarray(dg["risk"], np.float32)), "shipped risk table not reproduced"
    for lam, beta in CELLS:
        w = np.exp(-beta * zE)
        H = weighted_pool_histogram(Dp, w)
        if beta == 0.0:
            assert np.allclose(H, H0, atol=0, rtol=0), "beta = 0 != pool_histogram"
        P = (1.0 - lam) * P0 + lam * H
        C, top, rt = cell_cloud(P, seq, i, j, Wp, fold)
        if lam == 0.0:
            assert np.array_equal(np.asarray(rt["risk"]), np.asarray(rt0["risk"])), "lam = 0 is not the incumbent bit-exactly"
        row["cells"]["%g,%g" % (lam, beta)] = {"cloud": C, "top": top, "top75_overlap": float(len(set(top.tolist()) & set(top0.tolist())) / TOPM),
                                               "tri_cloud": float(I.ca_rmsd(C, C0))}
    for k in range(N_PERM):
        rng = SD.stable_rng("w_amberprior", pdb, k)
        zp = zE[rng.permutation(len(zE))]
        for lam, beta in CELLS:
            if beta == 0.0 or lam == 0.0:
                continue
            H = weighted_pool_histogram(Dp, np.exp(-beta * zp))
            P = (1.0 - lam) * P0 + lam * H
            C, top, _ = cell_cloud(P, seq, i, j, Wp, fold)
            row["perm"]["%d:%g,%g" % (k, lam, beta)] = {"cloud": C, "top": top}
    assert W.finite(row), pdb
    if verbose:
        tc = row["cells"]
        print("  %s f%d gate=%s | tri_cloud lam0.2: b0 %.3f b1 %.3f b4 %.3f | lam0.5: b0 %.3f b4 %.3f | overlap lam0.2,b1 %.2f"
              % (pdb, fold, row["gate"]["top75_equals_sub"], tc["0.2,0"]["tri_cloud"], tc["0.2,1"]["tri_cloud"], tc["0.2,4"]["tri_cloud"],
                 tc["0.5,0"]["tri_cloud"], tc["0.5,4"]["tri_cloud"], tc["0.2,1"]["top75_overlap"]), flush=True)
    del u
    return row


def clouds(probe=None, verbose=True):
    tg = I.targets()
    if probe:
        row = clouds_target([t for t in tg if t["pdb"] == probe][0], verbose=verbose)
        p = W.save("amberprior_probe_%s" % probe, {"label": "amber prior partner probe", "rows": [row]})
        print("  ->", p); return [row]
    name = "amberprior_clouds"
    prev = W.load_result(name)
    rows = prev["rows"] if prev else []
    done = {r["pdb"] for r in rows}
    t_last = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(clouds_target(t, verbose=verbose))
        if time.time() - t_last > 300 or len(rows) % 10 == 0:
            W.save(name, {"label": "amber prior partner, stage 1 native-free", "rows": rows}, rows=rows, n_expected=len(tg))
            t_last = time.time()
    summ = {"n": len(rows), "top75_equals_sub_count": int(sum(r["gate"]["top75_equals_sub"] for r in rows)),
            "tri_cloud_median": {c: float(np.median([r["cells"][c]["tri_cloud"] for r in rows])) for c in rows[0]["cells"]}}
    p = W.save(name, {"label": "amber prior partner, stage 1 native-free", "summary": summ, "rows": rows}, rows=rows,
               n_expected=len(tg), complete_keys=("pdb", "cells", "perm"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


def endpoint(verbose=True):
    W.require_signoff("amberprior endpoint")
    Z = W.load_result("amberprior_clouds")
    if not (Z and Z.get("complete")):
        raise SystemExit("run `clouds` to completion first")
    tg = W.targets_by_pdb()
    rows = Z["rows"]; pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    nat = {p: np.asarray(I.load_univ(p)["nat_ca"], float) for p in pdbs}
    # ---- cloud RMSD per cell (all 30) and per permuted cell
    cell_keys = list(rows[0]["cells"].keys())
    Mc = np.array([[I.ca_rmsd(np.asarray(r["cells"][c]["cloud"], float), nat[r["pdb"]]) for c in cell_keys] for r in rows])   # (126, 30)
    # ---- leave-fold-out choice on the point cloud (all cells; and beta = 0 row only)
    def lfo(keys_mask):
        choice = np.empty(len(rows), dtype=object)
        for f in sorted(set(folds.tolist())):
            train = folds != f
            means = Mc[train].mean(0)
            means = np.where(keys_mask, means, np.inf)
            k = int(np.argmin(means))                      # first in grid order among ties: nearest (0, 0)
            choice[folds == f] = cell_keys[k]
        return choice
    all_mask = np.ones(len(cell_keys), bool)
    b0_mask = np.array([c.endswith(",0") for c in cell_keys])
    choice_all = lfo(all_mask); choice_b0 = lfo(b0_mask)
    # ---- projections: shipped, chosen, chosen beta = 0, 4 permuted draws at the chosen cell (or the chosen cell itself if lam = 0 / beta = 0)
    out_path = os.path.join(W.RES, "w_selfcopy_amberprior_endpoint.json")
    prev = json.load(open(out_path)) if os.path.exists(out_path) else None
    done = {r["pdb"]: r for r in (prev["rows"] if prev else [])}
    res_rows = []
    t_last = time.time()
    for k, r in enumerate(rows):
        pdb = r["pdb"]
        if pdb in done:
            res_rows.append(done[pdb]); continue
        t = tg[pdb]; N = nat[pdb]
        ch, chb0 = choice_all[k], choice_b0[k]
        arms = {"shipped": np.asarray(r["cells"]["0,0"]["cloud"], float), "chosen": np.asarray(r["cells"][ch]["cloud"], float),
                "chosen_b0": np.asarray(r["cells"][chb0]["cloud"], float)}
        lam, beta = (float(x) for x in ch.split(","))
        perm = []
        for q in range(N_PERM):
            key = "%d:%g,%g" % (q, lam, beta)
            perm.append(np.asarray(r["perm"][key]["cloud"], float) if key in r["perm"] else arms["chosen"])
        e = {}
        for name, C in list(arms.items()) + [("perm%d" % q, perm[q]) for q in range(N_PERM)]:
            pr = I.project(C, t["seq"], int(t["fold"]))
            e[name] = {"cloud": float(I.ca_rmsd(C, N)), "arm": float(I.ca_rmsd(pr["ca"], N)), "fit": float(I.ca_rmsd(pr["fit_ca"], N))}
        e["perm_mean"] = {b: float(np.mean([e["perm%d" % q][b] for q in range(N_PERM)])) for b in ("cloud", "arm", "fit")}
        rr = {"pdb": pdb, "fold": int(t["fold"]), "chosen": ch, "chosen_b0": chb0, "rmsd": e,
              "cloud_all_cells": {c: float(Mc[k, ci]) for ci, c in enumerate(cell_keys)}}
        res_rows.append(rr)
        if verbose:
            print("  [%3d/126] %s chosen %s (b0: %s): shipped arm %.3f chosen %.3f b0 %.3f perm %.3f" % (k + 1, pdb, ch, chb0, e["shipped"]["arm"], e["chosen"]["arm"], e["chosen_b0"]["arm"], e["perm_mean"]["arm"]), flush=True)
        if time.time() - t_last > 300 or len(res_rows) % 10 == 0:
            W.save("amberprior_endpoint", {"label": "amber prior partner, GATED", "rows": res_rows}, rows=res_rows, n_expected=len(rows))
            t_last = time.time()
    out = {"label": "amber prior partner, GATED", "rows": res_rows, "stats": {}, "prereg": "s26/PREREG_amber_prior_partner.md",
           "chosen_cells_by_fold": {int(f): str(choice_all[folds == f][0]) for f in sorted(set(folds.tolist()))},
           "chosen_b0_by_fold": {int(f): str(choice_b0[folds == f][0]) for f in sorted(set(folds.tolist()))}}
    for b in ("arm", "cloud"):
        ship = np.array([r["rmsd"]["shipped"][b] for r in res_rows]); cho = np.array([r["rmsd"]["chosen"][b] for r in res_rows])
        b0 = np.array([r["rmsd"]["chosen_b0"][b] for r in res_rows]); pm = np.array([r["rmsd"]["perm_mean"][b] for r in res_rows])
        for name, a, c in (("chosen_minus_shipped", cho, ship), ("chosen_minus_b0mix", cho, b0), ("chosen_minus_permutedAMBER", cho, pm), ("b0mix_minus_shipped", b0, ship)):
            res = ST.compare(a, c, folds, names=pdbs, label="%s [%s], LFO cell on the cloud, n=126" % (name, b))
            out["stats"]["%s_%s" % (name, b)] = res
            if verbose:
                print(ST.fmt(res))
    # the per-target oracle over the 30 cells on the cloud, with its valid null
    out["oracle_over_cells_cloud"] = {"observed_min_mean": float(Mc.min(1).mean()), "shipped_cloud_mean": float(Mc[:, 0].mean()),
                                      "best_of_k_within": ST.best_of_k_within(Mc), "split_half_transfer": ST.split_half_transfer(Mc)}
    if verbose:
        o = out["oracle_over_cells_cloud"]
        print("  per-target oracle over 30 cells (cloud): %.4f vs shipped %.4f; best_of_k_within share_accounted %s k_eff %s; split-half transfer %s"
              % (o["observed_min_mean"], o["shipped_cloud_mean"], o["best_of_k_within"].get("share_accounted"), o["best_of_k_within"].get("k_eff"),
                 o["split_half_transfer"].get("transfer")))
    W.save("amberprior_endpoint", out, rows=res_rows, n_expected=len(rows), complete_keys=("pdb", "rmsd", "chosen"))
    print("  ->", out_path)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("clouds", "endpoint"))
    ap.add_argument("--probe", default=None)
    a = ap.parse_args(argv)
    if a.cmd == "clouds":
        clouds(probe=a.probe)
    else:
        endpoint()
    return 0


if __name__ == "__main__":
    sys.exit(main())
