#!/usr/bin/env python
"""s27/run_wave2.py -- PREREG addendum 1: H9 (m-ladder), H10 (in-band re-ranking), H11
(non-additive forms), H12 (near-corpus statistical potentials).  Point-cloud endpoint, the
cached channels of `run_pool.py`, per-target rows in `s27/results/wave2_rows.jsonl`.
"""
from __future__ import annotations

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
from s25 import phys_lib as P              # noqa: E402
from s27 import ham_lib as HL              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

ROWS = os.path.join(RP.RESULTS, "wave2_rows.jsonl")
MS = [3, 5, 10, 25, 50, 75, 100, 150]
LADDER_ARMS = ["DIS", "DISTPOT", "CONS", "DIS+DISTPOT", "DIS+CONS", "DIS+ENV", "DIS+SS_MATCH", "DIS+0.5*SS_MATCH"]
INBAND_X = ["DISTPOT", "CONS", "CONTACT_LL", "ENV", "DSSPHB", "POOLGO", "TORS_CONS", "LEG", "AMB", "SS_MATCH"]
NONADD_X = ["DISTPOT", "ENV", "SS_MATCH"]
RAMA = RP.RAMA


def energy(cfg, ch):
    if cfg in ch:
        return RP.zr(ch[cfg])
    if "*" in cfg:
        w, c = cfg[4:].split("*", 1)
        return RP.combine(ch, ["DIS", c], [1.0, float(w)])
    return RP.combine(ch, cfg.split("+"))


def near_corpus_channels(pdb, cand):
    """DISTPOT and CONTACT refitted on the 2,000 most BLOSUM-similar universe windows."""
    u = I.load_univ(pdb)
    order = np.asarray(u["order"], int)[:2000]
    sub = {k: (np.asarray(u[k])[order] if k in ("W", "S", "PHI", "PSI") else u[k]) for k in u}
    dg = I.distogram(pdb, cand.seq, cand.fold)
    cx = HL.Context(cand, sub, dg, RAMA[cand.fold])
    return {"NC_DISTPOT": HL.h_distpot(cx), "NC_CONTACT": HL.h_contact(cx)}


def main():
    done = set()
    if os.path.exists(ROWS):
        for line in open(ROWS, encoding="utf-8"):
            try:
                done.add(json.loads(line)["pdb"])
            except Exception:
                pass
    t0 = time.time()
    pdbs = P.targets()
    for n, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        cand, ch, _ = RP.channels_for(pdb)
        k = cand.k
        key = RP.rng_for(pdb, "tiekey").random(k)
        rows = []

        def rm(idx):
            return RP.oracle_rmsd_of_set(cand, idx)

        base = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold))
        # H9 m-ladder
        rr = RP.rng_for(pdb, "ladder_rand")
        for m in MS:
            rnd = float(np.mean([rm(rr.choice(k, size=m, replace=False)) for _ in range(8)]))
            for cfg in LADDER_ARMS:
                E = energy(cfg, ch)
                rows.append(dict(base, hyp="H9", config=cfg, m=m, rmsd=rm(RP.topm(E, m, key)), rand_m=rnd))
        # H10 in-band re-ranking: DIS top-150 -> best 75 by X
        top150 = RP.topm(ch["DIS"], 150, key)
        rr2 = RP.rng_for(pdb, "inband_rand")
        rnd150 = float(np.mean([rm(rr2.choice(top150, size=75, replace=False)) for _ in range(8)]))
        rows.append(dict(base, hyp="H10", config="DIS150->random75", m=75, rmsd=rnd150))
        rows.append(dict(base, hyp="H10", config="DIS150->DIS75", m=75, rmsd=rm(RP.topm(ch["DIS"], 75, key))))
        for x in INBAND_X:
            Ex = np.full(k, np.inf); Ex[top150] = ch[x][top150]
            rows.append(dict(base, hyp="H10", config=f"DIS150->{x}75", m=75, rmsd=rm(RP.topm(Ex, 75, key))))
            # and the additive in-band form: within the 150, rank by zrank(DIS)+zrank(X)
            Es = np.full(k, np.inf); Es[top150] = RP.zr(ch["DIS"][top150]) + RP.zr(ch[x][top150])
            rows.append(dict(base, hyp="H10", config=f"DIS150->(DIS+{x})75", m=75, rmsd=rm(RP.topm(Es, 75, key))))
        # H11 non-additive
        from scipy.stats import rankdata
        for x in NONADD_X:
            Emax = np.maximum(RP.zr(ch["DIS"]), RP.zr(ch[x]))
            Eprod = RP.zr(rankdata(ch["DIS"]) * rankdata(ch[x]))
            rows.append(dict(base, hyp="H11", config=f"max(DIS,{x})", m=75, rmsd=rm(RP.topm(Emax, 75, key))))
            rows.append(dict(base, hyp="H11", config=f"rankprod(DIS,{x})", m=75, rmsd=rm(RP.topm(Eprod, 75, key))))
        # H12 near-corpus potentials
        nc = near_corpus_channels(pdb, cand)
        ch2 = dict(ch); ch2.update(nc)
        for x in nc:
            rows.append(dict(base, hyp="H12", config=x, m=75, rmsd=rm(RP.topm(nc[x], 75, key))))
            rows.append(dict(base, hyp="H12", config=f"DIS+{x}", m=75, rmsd=rm(RP.topm(RP.combine(ch2, ["DIS", x]), 75, key))))
            rows.append(dict(base, hyp="H12", config=f"DIS+0.5*{x}", m=75, rmsd=rm(RP.topm(RP.combine(ch2, ["DIS", x], [1, 0.5]), 75, key))))
        rows.append(dict(base, hyp="ref", config="DIS", m=75, rmsd=rm(RP.topm(ch["DIS"], 75, key))))
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  [{n+1}/{len(pdbs)}] {pdb} rows={len(rows)} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", ROWS)


if __name__ == "__main__":
    main()
