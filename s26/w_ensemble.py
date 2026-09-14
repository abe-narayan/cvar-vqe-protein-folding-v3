#!/usr/bin/env python
"""s26/w_ensemble.py -- test-time ensembling over retrieval variants at FIXED K.  Lane W (lane
P's orphaned idea), Sprint 26.  Pre-registration: s26/PREREG_window_ensembling.md.

Three retrieval keys (BLOSUM45 / 62 / 80) at K = 500, or three K (250 / 500 / 1000) at BLOSUM62;
each pool cut to its own top-75 by the shipped Bayes-risk score, each top-75 averaged in its own
medoid frame, the three clouds superposed onto the production cloud and averaged, one projection.
No shortlist is widened (S17 L12's displacement mechanism cannot act).  Controls: the single
variant clouds, and `boot3` (three medoid-frame means of 75-member resamples of the production
top-75: a zero-information ensembling control matched in cloud and member count).

`clouds` is native-free (universes loaded blind, `rr` and `nat_ca` poisoned to NaN; every stored
array asserted finite) and stores the emissions; `endpoint` (GATED) reads natives and runs
`ST.compare`.  Everything reuses s26/w_selfcopy.py's emission path (production top-75 == `sub`
asserted on every target).

    python s26/jobrun.py --agent W --tag CPU --name w_ensemble_probe --est-ram 0.5 -- python s26/w_ensemble.py clouds --probe 1A13
    python s26/jobrun.py --agent W --tag CPU --name w_ensemble_clouds --est-ram 0.5 -- python s26/w_ensemble.py clouds
    python s26/jobrun.py --agent W --tag CPU --name w_ensemble_endpoint --est-ram 0.3 -- python s26/w_ensemble.py endpoint
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

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import data as D                           # noqa: E402
import w_selfcopy as W                               # noqa: E402  (imported, never edited here)

K, TOPM = W.K, W.TOPM
KEYS = (45, 62, 80)
KS = (250, 500, 1000)
ARMS = ("shipped", "b45", "b80", "ens3", "ensK", "union3", "boot3")
_M = {}


def blosum(key):
    """The BLOSUM<key> matrix re-ordered to core.data.ALPHABET; 62 asserted equal to core's."""
    if key not in _M:
        from Bio.Align import substitution_matrices as SM
        M = SM.load("BLOSUM%d" % key)
        A = "".join(M.alphabet)
        idx = [A.index(c) for c in D.ALPHABET]
        Mr = np.asarray(M, float)[np.ix_(idx, idx)]
        if key == 62:
            assert np.array_equal(Mr, np.asarray(D.BLOSUM62, float)), "Biopython BLOSUM62 != core.data.BLOSUM62"
        _M[key] = Mr
    return _M[key]


def sims(key, q, S):
    """Summed BLOSUM<key> similarity of every window (rows of S, codes) to the query codes."""
    return blosum(key)[np.asarray(q)[None, :], np.asarray(S)].sum(-1)


def cloud_of(Wp, dg, seq, fold):
    """The medoid-frame mean of the shipped score's top-75 of a pool (no projection)."""
    e = W.emit(Wp, dg, seq, fold, project=False)
    return e["cloud"], e["top"]


def combine(C_ref, others):
    """Mean of the reference cloud and the other clouds after superposing each onto it."""
    moved = [I.superpose_batch(np.asarray(C, float)[None], np.asarray(C_ref, float))[0] for C in others]
    return np.mean([np.asarray(C_ref, float)] + moved, axis=0)


def union_cloud(members, ref_medoid):
    """Uniform mean of all members (with multiplicity) superposed onto one reference member."""
    return I.superpose_batch(np.asarray(members, float), np.asarray(ref_medoid, float)).mean(0)


def boot_cloud(top_members, rng):
    """The medoid-frame mean of a 75-member resample WITH replacement of the production top-75."""
    idx = rng.integers(0, len(top_members), len(top_members))
    Wb = np.asarray(top_members, float)[idx]
    P = I.pairwise_rmsd(Wb)
    return I.superpose_batch(Wb, Wb[I.medoid(P)]).mean(0)


def clouds_target(t, verbose=True):
    from s15 import seed as SD
    pdb, seq, n, fold = t["pdb"], t["seq"], t["n"], t["fold"]
    u = W.load_blind(pdb)
    Wu = u["W"]; S = np.asarray(u["S"]); order = np.asarray(u["order"], int); sim62 = np.asarray(u["sim"], float)
    dg = I.distogram(pdb); q = D.encode(seq)
    # ---- the three keys at K = 500
    pools, C, tops = {}, {}, {}
    for key in KEYS:
        s = sims(key, q, S)
        if key == 62:
            assert np.array_equal(s, sim62), "BLOSUM62 sums do not reproduce the universe's sim"
        pool = D.top_k(s, K)
        if key == 62:
            assert np.array_equal(pool, order[:K]), "top_k does not reproduce the pinned pool"
        pools[key] = pool
        C[key], tops[key] = cloud_of(Wu[pool], dg, seq, fold)
    gate = W.production_gate(pdb, {"top": tops[62], "cloud": C[62]})
    # ---- the three K at BLOSUM62
    CK, topsK = {}, {}
    for k in KS:
        pool = D.top_k(sim62, k)
        CK[k], topsK[k] = cloud_of(Wu[pool], dg, seq, fold)
    assert np.allclose(CK[500], C[62])
    # ---- the arms (clouds), then one projection each
    arms = {"shipped": C[62], "b45": C[45], "b80": C[80],
            "ens3": combine(C[62], [C[45], C[80]]),
            "ensK": combine(CK[500], [CK[250], CK[1000]])}
    members = np.concatenate([Wu[pools[key]][tops[key]] for key in KEYS], 0)          # 225, with multiplicity
    top62 = Wu[pools[62]][tops[62]]
    P62 = I.pairwise_rmsd(top62); ref = top62[I.medoid(P62)]
    arms["union3"] = union_cloud(members, ref)
    rng = SD.stable_rng("w_ensemble", pdb, 0)
    arms["boot3"] = combine(boot_cloud(top62, rng), [boot_cloud(top62, rng), boot_cloud(top62, rng)])
    out = {"pdb": pdb, "n": n, "fold": fold, "gate": gate, "arms": {}, "diag": {}}
    chains = {}
    for name in ARMS:
        pr = I.project(arms[name], seq, int(fold))
        chains[name] = pr["ca"]
        out["arms"][name] = {"cloud": arms[name], "chain": pr["ca"], "fit": pr["fit_ca"]}
    # ---- native-free diagnostics
    g = lambda a, b: len(set(a.tolist()) & set(b.tolist()))
    out["diag"] = {
        "pool_overlap_45_62": g(pools[45], pools[62]) / K, "pool_overlap_80_62": g(pools[80], pools[62]) / K,
        "top75_overlap_45_62": g(pools[45][tops[45]], pools[62][tops[62]]) / TOPM,
        "top75_overlap_80_62": g(pools[80][tops[80]], pools[62][tops[62]]) / TOPM,
        "top75_overlap_K250_K500": g(D.top_k(sim62, 250)[topsK[250]], pools[62][tops[62]]) / TOPM,
        "top75_overlap_K1000_K500": g(D.top_k(sim62, 1000)[topsK[1000]], pools[62][tops[62]]) / TOPM,
        "n_distinct_union": int(len({int(x) for key in KEYS for x in pools[key][tops[key]]})),
        "cloud_rmsd_45_62": float(I.ca_rmsd(C[45], C[62])), "cloud_rmsd_80_62": float(I.ca_rmsd(C[80], C[62])),
        "cloud_rmsd_45_80": float(I.ca_rmsd(C[45], C[80])),
        "cloud_rmsd_K250_K500": float(I.ca_rmsd(CK[250], CK[500])), "cloud_rmsd_K1000_K500": float(I.ca_rmsd(CK[1000], CK[500])),
        "tri_cloud": {name: float(I.ca_rmsd(arms[name], C[62])) for name in ARMS},
        "tri_arm": {name: float(I.ca_rmsd(chains[name], chains["shipped"])) for name in ARMS}}
    assert W.finite(out), pdb
    if verbose:
        d = out["diag"]
        print("  %s f%d gate=%s | pool ovl 45/80 %.2f/%.2f top75 ovl %.2f/%.2f | cloud rmsd 45,80 vs 62 %.3f %.3f | tri_arm %s"
              % (pdb, fold, gate["top75_equals_sub"], d["pool_overlap_45_62"], d["pool_overlap_80_62"],
                 d["top75_overlap_45_62"], d["top75_overlap_80_62"], d["cloud_rmsd_45_62"], d["cloud_rmsd_80_62"],
                 " ".join("%s:%.3f" % (a, d["tri_arm"][a]) for a in ARMS if a != "shipped")), flush=True)
    del u
    return out


def clouds(probe=None, verbose=True):
    tg = I.targets()
    if probe:
        row = clouds_target([t for t in tg if t["pdb"] == probe][0], verbose=verbose)
        p = W.save("ensemble_probe_%s" % probe, {"label": "ensemble probe", "rows": [row]})
        print("  ->", p); return [row]
    name = "ensemble_clouds"
    prev = W.load_result(name)
    rows = prev["rows"] if prev else []
    done = {r["pdb"] for r in rows}
    t_last = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(clouds_target(t, verbose=verbose))
        if time.time() - t_last > 300 or len(rows) % 10 == 0:
            W.save(name, {"label": "fixed-K ensembling, native-free", "rows": rows}, rows=rows, n_expected=len(tg))
            t_last = time.time()
    summ = {"n": len(rows), "top75_equals_sub_count": int(sum(r["gate"]["top75_equals_sub"] for r in rows))}
    for k in ("pool_overlap_45_62", "pool_overlap_80_62", "top75_overlap_45_62", "top75_overlap_80_62",
              "cloud_rmsd_45_62", "cloud_rmsd_80_62", "cloud_rmsd_45_80", "n_distinct_union"):
        summ[k + "_median"] = float(np.median([r["diag"][k] for r in rows]))
    summ["tri_arm_median"] = {a: float(np.median([r["diag"]["tri_arm"][a] for r in rows])) for a in ARMS}
    p = W.save(name, {"label": "fixed-K ensembling, native-free", "summary": summ, "rows": rows}, rows=rows,
               n_expected=len(tg), complete_keys=("pdb", "arms", "diag"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


def endpoint(verbose=True):
    W.require_signoff("ensemble endpoint")
    Z = W.load_result("ensemble_clouds")
    if not (Z and Z.get("complete")):
        raise SystemExit("run `clouds` to completion first")
    rows = []
    helix = None
    for r in Z["rows"]:
        N = np.asarray(I.load_univ(r["pdb"])["nat_ca"], float)
        e = {name: {"arm": float(I.ca_rmsd(np.asarray(a["chain"], float), N)),
                    "cloud": float(I.ca_rmsd(np.asarray(a["cloud"], float), N)),
                    "fit": float(I.ca_rmsd(np.asarray(a["fit"], float), N))} for name, a in r["arms"].items()}
        rows.append({"pdb": r["pdb"], "fold": r["fold"], "rmsd": e})
    pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    out = {"label": "GATED fixed-K ensembling endpoints", "rows": rows, "stats": {},
           "basis": "arm PRIMARY (built chain, 3.2148 basis); cloud carried (3.0483 basis)", "prereg": "s26/PREREG_window_ensembling.md"}
    ship = {b: np.array([r["rmsd"]["shipped"][b] for r in rows]) for b in ("arm", "cloud", "fit")}
    out["shipped_means"] = {b: float(v.mean()) for b, v in ship.items()}
    for name in ARMS:
        if name == "shipped":
            continue
        for b in ("arm", "cloud"):
            a = np.array([r["rmsd"][name][b] for r in rows])
            res = ST.compare(a, ship[b], folds, names=pdbs, label="%s minus shipped [%s]" % (name, b))
            out["stats"]["%s_%s" % (name, b)] = res
            if verbose:
                print(ST.fmt(res))
    for b in ("arm", "cloud"):
        a = np.array([r["rmsd"]["ens3"][b] for r in rows]); c = np.array([r["rmsd"]["boot3"][b] for r in rows])
        res = ST.compare(a, c, folds, names=pdbs, label="ens3 minus boot3 (the matched zero-information ensembling control) [%s]" % b)
        out["stats"]["ens3_minus_boot3_%s" % b] = res
        if verbose:
            print(ST.fmt(res))
    p = W.save("ensemble_endpoint", out)
    print("  ->", p)
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
