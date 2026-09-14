#!/usr/bin/env python
"""s27/run_vqe_chain.py -- the two secondary arms for the finalists (`s27/PREREG.md` H7 and
the built-chain endpoint):

  --vqe    the GENUINE CVaR-VQE selection (`s24.d_harness.arm_vqe`, S25 settings: 9 qubits over
           500 + 12 padding, 3 layers, 80 Adam iterations on the exact parameter-shift gradient,
           alpha 0.18, T 0.5) of the SAME energy vector, seeds 0 and 1, with the set-equality
           gate and the m-ladder decomposition
  --chain  the BUILT-CHAIN endpoint: the top-75 coordinate average projected through the
           production projection (`s12.instrument.project`, ramah 0.3, multi-start)

Configurations are named exactly as in `pool_rows.jsonl` and rebuilt from `s27/cache/<pdb>.npz`
with the same combine/tie rules as `run_pool.py`.  Rows go to `s27/results/vqe_rows.jsonl` /
`chain_rows.jsonl`, resumable per (config, pdb).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

RESULTS = RP.RESULTS
ALPHA, TEMP, LAYERS, ITERS = 0.18, 0.5, 3, 80


def energy_for(config, ch, pdb, key):
    """Rebuild the energy vector of a named configuration from the cached channels."""
    k = len(ch["DIS"])
    if config in ch:
        return RP.zr(ch[config])
    if config.startswith("REJ[") and config.endswith("]->DIS"):
        c = config[4:-6]
        q = int(round(RP.REJECT_Q * k))
        if c == "random":
            keep = RP.rng_for(pdb, "rejrand").choice(k, size=k - q, replace=False)
        else:
            keep = RP.topm(ch[c], k - q, key)
        E = np.full(k, np.inf); E[keep] = ch["DIS"][keep]
        fin = np.isfinite(E)
        E2 = np.empty(k); E2[fin] = RP.zr(E[fin]); E2[~fin] = E2[fin].max() + 1.0
        return E2
    if config.startswith("DIS+adapt("):
        c = config[len("DIS+adapt("):-1]
        dg = I.distogram(pdb)
        pr = np.clip(np.asarray(dg["prob"], float), 1e-12, 1)
        Hbar = float(-(pr * np.log(pr)).sum(1).mean())
        return RP.combine(ch, ["DIS", c], [1.0, 0.5 * Hbar / np.log(17.0)])
    if config in RP.COMPOSITES:
        return RP.combine(ch, RP.COMPOSITES[config])
    if config.startswith("DIS+") and config[4:] in RP.COMPOSITES:
        parts = RP.COMPOSITES[config[4:]]
        return RP.combine(ch, ["DIS"] + parts, [len(parts)] + [1.0] * len(parts))
    if config.startswith("DIS+") and "*" in config:
        w, c = config[4:].split("*", 1)
        return RP.combine(ch, ["DIS", c], [1.0, float(w)])
    if config.startswith("DIS+"):
        parts = config[4:].split("+")
        return RP.combine(ch, ["DIS"] + parts)
    raise KeyError(config)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vqe", action="store_true")
    ap.add_argument("--chain", action="store_true")
    ap.add_argument("--configs", required=True, help="comma-separated configuration names")
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    configs = [c for c in a.configs.split(",") if c]
    seeds = [int(s) for s in a.seeds.split(",")]
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    for mode in (["vqe"] if a.vqe else []) + (["chain"] if a.chain else []):
        path = os.path.join(RESULTS, f"{mode}_rows.jsonl")
        done = set()
        if os.path.exists(path):
            with open(path, encoding="utf-8") as fh:
                for line in fh:
                    try:
                        r = json.loads(line); done.add((r["config"], r["pdb"], r.get("seed", 0)))
                    except Exception:
                        pass
        t0 = time.time()
        for n, pdb in enumerate(pdbs):
            todo = [(c, s) for c in configs for s in (seeds if mode == "vqe" else [0]) if (c, pdb, s) not in done]
            if not todo:
                continue
            cand, ch, _ = RP.channels_for(pdb)
            key = RP.rng_for(pdb, "tiekey").random(cand.k)
            rows = []
            for c, s in todo:
                E = energy_for(c, ch, pdb, key)
                order = np.lexsort((key, E))
                t1 = time.time()
                if mode == "vqe":
                    q = H.arm_vqe(cand, E, alpha=ALPHA, T=TEMP, layers=LAYERS, iters=ITERS, seed=s)
                    m = int(q["m"])
                    gate = H.gate_set_equality(E, q["cands"], m)
                    r_vqe = RP.oracle_rmsd_of_set(cand, q["cands"])
                    r_topm = RP.oracle_rmsd_of_set(cand, order[:m])
                    r_top75 = RP.oracle_rmsd_of_set(cand, order[:RP.M])
                    rows.append(dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), config=c, seed=s,
                                     rmsd_vqe=r_vqe, rmsd_topm=r_topm, rmsd_top75=r_top75, m=m,
                                     eps=float(r_vqe - r_topm), mladder=float(r_topm - r_top75),
                                     entropy_bits=float(q["entropy_bits"]), ess=float(q["ess"]),
                                     cvar=float(q["cvar"]), gate_pass=bool(gate["pass_"]),
                                     gate_equality=bool(gate["equality"]), secs=float(time.time() - t1)))
                else:
                    top = order[:RP.M]
                    C, _ = H.readout_uniform(cand, top)
                    ca = H.readout_projected(cand, C)
                    rows.append(dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), config=c, seed=0,
                                     rmsd_cloud=float(I.ca_rmsd(C, cand.nat_ca)),
                                     rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)), secs=float(time.time() - t1)))
            with open(path, "a", encoding="utf-8") as fh:
                for r in rows:
                    fh.write(json.dumps(r) + "\n")
            print(f"  [{mode} {n+1}/{len(pdbs)}] {pdb} {len(rows)} rows {sum(r['secs'] for r in rows):.1f}s "
                  f"(elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
        print("done:", path)


if __name__ == "__main__":
    main()
