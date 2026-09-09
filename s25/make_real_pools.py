"""s25/make_real_pools.py -- GENERATE THE REAL PER-TARGET STRUCTURES FOR ALL SEVEN CONFIGURATIONS.

Replaces the quarantined proof build (native + noise) with genuine predictions, produced by exactly
the pipeline the seven-configuration suite measured: for each target and each configuration, the
channel scores are combined under rank normalisation, a GENUINE CVaR-VQE selects the tail, and the
retained set is averaged into the emitted point cloud.

Everything here is the suite's own code path -- `s25.phys_lib` for the channels and the combination,
`s24.d_harness.arm_vqe` for the selector -- with the SAME alpha, T, layers, iters and seed, so every
structure written corresponds cell-for-cell to a row already in `s25/results/phys_suite.json`. The
per-target RMSD of the emitted cloud is asserted against that row before anything is saved.

OUTPUT: one npz per configuration, holding parallel `pdb` and `ca` arrays, plus `production` (the
shipped pipeline's own emitted cloud) and a spec JSON for `s25/resultslab/build.py --mode frozen`.
The results lab then handles the ideal-geometry projection to the BUILT CHAIN, the PDB export with
mandatory provenance headers, the four release gates, and the site.

NOTHING HERE IS SYNTHETIC. Every configuration is marked `is_synthetic: False` explicitly rather
than by default, because `providers.from_npz` treats absent provenance as UNKNOWN and fails the
frozen build -- which is the behaviour we want and must not be satisfied by accident.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s25 import phys_lib as P             # noqa: E402
from s24 import d_harness as H            # noqa: E402

OUT = os.path.join(HERE, "results", "real_pools")
os.makedirs(OUT, exist_ok=True)
TOL = 5e-9


def run(limit=None):
    pdbs = P.targets()
    if limit:
        pdbs = pdbs[:limit]
    #: the suite's own per-target rows, so every structure can be checked against a measured number
    suite = {}
    for p in pdbs:
        f = os.path.join(HERE, "results", "suite_cells", "%s.json" % p)
        if os.path.exists(f):
            for r in json.load(open(f))["rows"]:
                if r.get("norm") == "rank":
                    suite[(p, int(r["config"]))] = float(r["rmsd_vqe"])

    store = {cid: {} for cid, _, _ in P.CONFIGS}
    prod = {}
    bad = []
    print("targets: %d, configs: %d, alpha=%.2f T=%.1f seed=%d"
          % (len(pdbs), len(P.CONFIGS), P.ALPHA, P.TEMP, P.SEED), flush=True)

    for c_i, pdb in enumerate(pdbs):
        cand, ch = P.channels(pdb)                      # asserts the AMBER cache identity
        for cid, name, subset in P.CONFIGS:
            E = P.combine(ch, subset, norm="rank")
            q = H.arm_vqe(cand, E, alpha=P.ALPHA, T=P.TEMP,
                          layers=P.LAYERS, iters=P.ITERS, seed=P.SEED)
            C, _ = I.coordinate_average(cand.W[np.asarray(q["cands"], int)])
            C = np.asarray(C, float)
            store[cid][pdb] = C
            got = float(I.ca_rmsd(C, cand.nat_ca))
            want = suite.get((pdb, cid))
            if want is not None and abs(got - want) > TOL:
                bad.append((pdb, cid, got, want, abs(got - want)))
        #: the shipped pipeline's own emitted cloud, for the `production` row
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        Cp, _ = I.coordinate_average(W[np.argsort(sc, kind="stable")[:75]])
        prod[pdb] = np.asarray(Cp, float)

        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(pdbs)), flush=True)

    if bad:
        print("\n*** %d cells disagree with the suite beyond %.0e ***" % (len(bad), TOL))
        for r in bad[:10]:
            print("    %s cfg%d  got %.9f  want %.9f  d=%.2e" % r)
        raise SystemExit("refusing to write: structures do not reproduce the measured suite")

    slug = {1: "legacy", 2: "amber", 3: "distogram", 4: "legacy_distogram",
            5: "amber_distogram", 6: "legacy_amber", 7: "legacy_amber_distogram"}
    spec = {"basis_note": "emitted point clouds; the results lab projects to the built chain",
            "source": "s25/make_real_pools.py", "configurations": {}}
    for cid, name, subset in P.CONFIGS:
        d = store[cid]
        ks = sorted(d)
        path = os.path.join(OUT, "%s.npz" % slug[cid])
        np.savez_compressed(path, pdb=np.array(ks, object),
                            ca=np.array([d[k] for k in ks], object))
        spec["configurations"][slug[cid]] = {
            "path": os.path.relpath(path, ROOT).replace("\\", "/"),
            "meta": {"is_synthetic": False, "selector": "CVaR-VQE",
                     "hamiltonians": [{"LEG":"Legacy","AMB":"AMBER","DIS":"Distogram"}[c] for c in subset],
                     "distogram_used": ("DIS" in subset), "config_name": name,
                     "alpha": P.ALPHA, "T": P.TEMP, "qubits": 9, "layers": P.LAYERS,
                     "norm": "rank", "candidates": 500, "retained": "CVaR tail"}}
    ks = sorted(prod)
    ppath = os.path.join(OUT, "production.npz")
    np.savez_compressed(ppath, pdb=np.array(ks, object),
                        ca=np.array([prod[k] for k in ks], object))
    spec["configurations"]["production"] = {
        "path": os.path.relpath(ppath, ROOT).replace("\\", "/"),
        "meta": {"is_synthetic": False, "selector": "score-filter top-75 (== CVaR tail, s24 theorem)",
                 "hamiltonians": ["Distogram"], "distogram_used": True, "config_name": "Production", "candidates": 500,
                 "retained": 75}}
    sp = os.path.join(OUT, "spec.json")
    with open(sp, "w") as fh:
        json.dump(spec, fh, indent=1)

    print("\nAll %d x %d cells reproduce the suite to <%.0e." % (len(pdbs), len(P.CONFIGS), TOL))
    for cid, name, _ in P.CONFIGS:
        v = np.array([I.ca_rmsd(store[cid][p], I.load_univ(p)["nat_ca"]) for p in sorted(store[cid])])
        print("  %-26s point cloud mean %.4f" % (name, v.mean()))
    v = np.array([I.ca_rmsd(prod[p], I.load_univ(p)["nat_ca"]) for p in sorted(prod)])
    print("  %-26s point cloud mean %.4f" % ("Production", v.mean()))
    print("\nspec: %s" % sp)


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else None)
