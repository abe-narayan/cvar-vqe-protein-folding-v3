"""s26/q_mde.py -- the reference MDEs for the A1/A3 pre-registrations, from stored artefacts.

Reads the per-target arms of bench_results/cache/464a0ddb5f283e04/<pdb>.json, the 126-target
run of Config(quantum=True, legacy=True) (the only production cache with the quantum arms),
and the per-target arms of s8/integrate_vqe.json (the s8 selection instrument), and puts the
matched contrasts through s24.stats_lib.compare.  No new structure, no new RMSD: every number
here is a re-read of an RMSD that was written before Sprint 26.  MDE = 2.8016 * SE per
comparison.  Output: s26/results/q_mde_reference.json.
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

from s24 import stats_lib as ST                                           # noqa: E402

CACHE = os.path.join(ROOT, "bench_results", "cache", "464a0ddb5f283e04")
OUT = os.path.join(HERE, "results", "q_mde_reference.json")


def main():
    fs = sorted(f for f in os.listdir(CACHE) if f.endswith(".json"))
    rows = [json.load(open(os.path.join(CACHE, f))) for f in fs]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([int(r["fold"]) for r in rows])
    out = {"source_cache": CACHE, "n": len(rows), "contrasts": {}, "means": {}}

    def col(k):
        return np.array([r[k] for r in rows], float)

    for k in ("rmsd_arm", "rmsd_avg", "rmsd_full", "rmsd_q_synth", "rmsd_u_synth", "rmsd_q_avg",
              "rmsd_vqe_sel", "rmsd_vqe_sel_uniform", "shipped", "pool_best"):
        out["means"][k] = float(col(k).mean())
    pairs = [("rmsd_q_synth", "rmsd_u_synth", "built chain: p-weighted vs uniform over the same 128"),
             ("rmsd_q_synth", "rmsd_arm", "built chain: quantum synthesis vs the production top-75 arm"),
             ("rmsd_u_synth", "rmsd_arm", "built chain: uniform-128 synthesis vs the production arm"),
             ("rmsd_q_avg", "rmsd_avg", "point cloud (NOT a structure): weighted vs top-75 average"),
             ("rmsd_vqe_sel", "rmsd_vqe_sel_uniform", "selection: p-weighted medoid vs unweighted medoid"),
             ("rmsd_vqe_sel", "shipped", "selection: VQE_LFO medoid vs the shipped argmin")]
    for a, b, what in pairs:
        r = ST.compare(col(a), col(b), folds, names=pdbs, label=f"{what} [{a} - {b}]")
        print(ST.fmt(r))
        out["contrasts"][f"{a}-{b}"] = r
    # the s8 instrument's circuit-vs-Boltzmann contrast, re-derived from its per-target arrays
    V = json.load(open(os.path.join(ROOT, "s8", "integrate_vqe.json")))
    PT = {k: np.asarray(v, float) for k, v in V["per_target"].items()}
    F8 = np.asarray(V["folds"], int)
    for a, b in (("vqe_a1.0_T0.3", "boltz_T0.3"), ("vqe_a0.25_T0.3", "boltz_T0.3"),
                 ("vqe_LFO", "boltz_T0.3"), ("vqe_LFO", "argmin")):
        r = ST.compare(PT[a], PT[b], F8, names=V["pdbs"], label=f"s8 instrument: {a} - {b}")
        print(ST.fmt(r))
        out["contrasts"][f"s8:{a}-{b}"] = r
    ST.save_atomic(OUT, out, module_file=__file__)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
