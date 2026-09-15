#!/usr/bin/env python
"""s27/s28_D_reproduce.py -- lane D (Adversary): reproduce ONE S27 number from its artefact.

Picks one row at random (stated seed) from `s27/results/pool_rows.jsonl` (point-cloud rows,
`s27/run_pool.py :: evaluate_target`), recomputes the configuration's top-75 point-cloud RMSD
through the same instrument (`s27.run_pool.channels_for`, `topm`, `P.rmsd_of_set`) and
compares it with the stored value.  Optionally reproduces one `chain_rows.jsonl` row through
`s12.instrument.project` (the built chain, 3 s) and one `vqe_rows.jsonl` row through
`s24.d_harness.arm_vqe` (the genuine CVaR-VQE).

Usage:  python s27/s28_D_reproduce.py --seed 7 [--kind pool|chain|vqe] [--out path]

ORACLE reads happen only inside `P.rmsd_of_set` / `I.ca_rmsd` (the endpoint).
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
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import run_vqe_chain as RVC       # noqa: E402

RESULTS = os.path.join(HERE, "results")


def _rows(name):
    with open(os.path.join(RESULTS, name), encoding="utf-8") as f:
        return [json.loads(l) for l in f if l.strip()]


def _energy_for(config, ch, k, key, pdb):
    """Rebuild the energy vector of a `pool_rows` configuration from the cached channels,
    exactly as `s27/run_pool.py :: evaluate_target` does."""
    import re
    if config in ch:
        return ch[config]
    m = re.fullmatch(r"DIS\+(\d*\.?\d+)\*(.+)", config)
    if m:
        return RP.combine(ch, ["DIS", m.group(2)], [1.0, float(m.group(1))])
    m = re.fullmatch(r"DIS\+(.+)~perm", config)
    if m:
        c = m.group(1)
        perm_rng = RP.rng_for(pdb, "perm")
        perm = {cc: perm_rng.permutation(k) for cc in ch}
        chp = dict(ch); chp[c] = ch[c][perm[c]]
        return RP.combine(chp, ["DIS", c])
    m = re.fullmatch(r"REJ\[(.+)\]->DIS", config)
    if m:
        c = m.group(1)
        q = int(round(RP.REJECT_Q * k))
        if c == "random":
            keep = RP.rng_for(pdb, "rejrand").choice(k, size=k - q, replace=False)
        else:
            keep = RP.topm(ch[c], k - q, key)
        Ek = np.full(k, np.inf); Ek[keep] = ch["DIS"][keep]
        return Ek
    if config in RP.COMPOSITES:
        return RP.combine(ch, RP.COMPOSITES[config])
    m = re.fullmatch(r"DIS\+(ROSETTA_LIKE|PHYSICS_COMP|CONSIST_COMP|STAT_COMP)", config)
    if m:
        parts = RP.COMPOSITES[m.group(1)]
        return RP.combine(ch, ["DIS"] + parts, [len(parts)] + [1.0] * len(parts))
    m = re.fullmatch(r"DIS\+adapt\((.+)\)", config)
    if m:
        raise NotImplementedError("adaptive rows need the distogram entropy; pick another row")
    m = re.fullmatch(r"DIS\+(.+)", config)
    if m:
        return RP.combine(ch, ["DIS", m.group(1)])
    raise KeyError(config)


def reproduce_pool(row):
    pdb, cfg = row["pdb"], row["config"]
    cand, ch, _ = RP.channels_for(pdb)
    k = cand.k
    key = RP.rng_for(pdb, "tiekey").random(k)
    E = _energy_for(cfg, ch, k, key, pdb)
    top = RP.topm(E, RP.M, key)
    r = P.rmsd_of_set(cand, top)
    dis_top = RP.topm(ch["DIS"], RP.M, key)
    r_dis = P.rmsd_of_set(cand, dis_top)
    return dict(kind="pool", pdb=pdb, config=cfg, stored=float(row["rmsd"]), recomputed=float(r),
                abs_diff=abs(float(r) - float(row["rmsd"])),
                stored_dis75=float(row["rmsd_dis75"]), recomputed_dis75=float(r_dis),
                abs_diff_dis75=abs(float(r_dis) - float(row["rmsd_dis75"])))


def reproduce_chain(row):
    pdb, cfg = row["pdb"], row["config"]
    cand, ch, _ = RP.channels_for(pdb)
    k = cand.k
    key = RP.rng_for(pdb, "tiekey").random(k)
    E = RVC.energy_for(cfg, ch, pdb, key)
    order = np.lexsort((key, E))
    top = order[:RP.M]
    C, _ = I.coordinate_average(cand.W[np.asarray(top, int)])
    r_cloud = float(I.ca_rmsd(C, cand.nat_ca))
    pr = I.project(np.asarray(C, float), cand.seq, cand.fold)
    r_chain = float(I.ca_rmsd(np.asarray(pr["ca"], float), cand.nat_ca))
    return dict(kind="chain", pdb=pdb, config=cfg,
                stored_cloud=float(row["rmsd_cloud"]), recomputed_cloud=r_cloud,
                abs_diff_cloud=abs(r_cloud - float(row["rmsd_cloud"])),
                stored_chain=float(row["rmsd_chain"]), recomputed_chain=r_chain,
                abs_diff_chain=abs(r_chain - float(row["rmsd_chain"])))


def reproduce_vqe(row):
    from s24 import d_harness as H
    pdb, cfg, seed = row["pdb"], row["config"], int(row["seed"])
    cand, ch, _ = RP.channels_for(pdb)
    k = cand.k
    key = RP.rng_for(pdb, "tiekey").random(k)
    E = RVC.energy_for(cfg, ch, pdb, key)
    v = H.arm_vqe(cand, E, alpha=RVC.ALPHA, T=RVC.TEMP, layers=RVC.LAYERS, iters=RVC.ITERS, seed=seed)
    r = P.rmsd_of_set(cand, v["cands"])
    return dict(kind="vqe", pdb=pdb, config=cfg, seed=seed,
                stored=float(row["rmsd_vqe"]), recomputed=float(r),
                abs_diff=abs(float(r) - float(row["rmsd_vqe"])),
                stored_m=int(row["m"]), recomputed_m=int(v["m"]), secs=float(v["secs"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--kind", default="pool", choices=["pool", "chain", "vqe"])
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    fn = {"pool": "pool_rows.jsonl", "chain": "chain_rows.jsonl", "vqe": "vqe_rows.jsonl"}[a.kind]
    rows = _rows(fn)
    rng = np.random.default_rng(a.seed)
    # skip adaptive rows (need the distogram entropy) and staged rows for the vqe kind
    for _ in range(50):
        i = int(rng.integers(0, len(rows)))
        row = rows[i]
        if row["config"].startswith("DIS+adapt("):
            continue
        break
    t0 = time.time()
    fn_ = {"pool": reproduce_pool, "chain": reproduce_chain, "vqe": reproduce_vqe}[a.kind]
    res = fn_(row)
    res.update(seed=a.seed, row_index=i, n_rows=len(rows), source=fn, secs=time.time() - t0)
    print(json.dumps(res, indent=1))
    out = a.out or os.path.join(RESULTS, f"s28_D_reproduce_{a.kind}_seed{a.seed}.json")
    ST.save_atomic(out, res, complete_keys=["stored" if a.kind != "chain" else "stored_chain",
                                             "recomputed" if a.kind != "chain" else "recomputed_chain"],
                   module_file=__file__)
    print("wrote", out)


if __name__ == "__main__":
    main()
