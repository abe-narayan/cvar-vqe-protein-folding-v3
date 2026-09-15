#!/usr/bin/env python
"""s27/s28_B_train.py -- F5 of `s27/PREREG_S28_B.md`: gradient variance of F vs J at n = 4..9,
checkpointed per target (the governor killed the in-module version at 95.7% box RAM after three
targets; this one appends `s28_B_train_rows.jsonl` per target and resumes).

Property measurement: no RMSD, no native. Uses `s27.s28_B_hop.measure_hop` unchanged.
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

from s22 import qcand_lib as QC            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

ROWS = os.path.join(B.RESULTS, "s28_B_train_rows.jsonl")
NS = [4, 5, 6, 7, 8, 9]
N_THETA = 120


def target_rows(pdb: str):
    cand, ch, _ = RP.channels_for(pdb)
    E_full = RP.zr(ch["DIS"])
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    order = RP.topm(E_full, cand.k, key)
    D = B.pairwise_rmsd_matrix(cand.W)
    rows = []
    for n in NS:
        dim = 1 << n
        if dim <= cand.k:
            sub = order[:dim]
            g = B.kernel_graph(D[np.ix_(sub, sub)])
            A = g["A"]
            E = B.deployed_E(dim)
        else:
            g = B.kernel_graph(D)
            A = B.pad_graph(g["A"], dim)
            E = QC.Encoding(E_full).E
        for J in B.J_GRID:
            r = B.measure_hop(n, B.LAYERS, B.ALPHA, B.TEMP, J, N_THETA, 1009, 0.6, E, A, "full")
            r.update(pdb=pdb, sigma=g["sigma"], lam_max=g["lam_max"]); rows.append(r)
        r = B.measure_hop(n, B.LAYERS, B.ALPHA, B.TEMP, 1.0, N_THETA, 1009, 0.6, E, A, "hop_only")
        r.update(pdb=pdb, sigma=g["sigma"], lam_max=g["lam_max"]); rows.append(r)
        r = B.measure_hop(n, B.LAYERS, 1.0, 0.0, 0.0, N_THETA, 1009, 0.6, E, A, "linear")
        r.update(pdb=pdb, sigma=g["sigma"], lam_max=g["lam_max"]); rows.append(r)
    return rows


def summarise(rows):
    summ = {}
    cells = sorted({(r["mode"], r["J"]) for r in rows})
    for mode, J in cells:
        per_n = {}
        for n in NS:
            v = [r["var_g0"] for r in rows if r["mode"] == mode and r["J"] == J and r["n"] == n]
            m = [r["mean_sq_per_param"] for r in rows if r["mode"] == mode and r["J"] == J and r["n"] == n]
            per_n[str(n)] = dict(var_g0_median=float(np.median(v)), var_g0_min=float(min(v)),
                                 var_g0_max=float(max(v)), msq_median=float(np.median(m)), n_targets=len(v))
        ns_ = np.array(NS, float)
        v = np.array([per_n[str(n)]["var_g0_median"] for n in NS])
        ok = v > 0
        slope = float(np.polyfit(ns_[ok], np.log2(v[ok]), 1)[0]) if ok.sum() > 2 else float("nan")
        m_ = np.array([per_n[str(n)]["msq_median"] for n in NS])
        slope_m = float(np.polyfit(ns_, np.log2(m_), 1)[0])
        summ[f"{mode}|J{J:g}"] = dict(per_n=per_n, log2_slope_per_qubit=slope, log2_slope_msq_per_qubit=slope_m)
    base = summ["full|J0"]["per_n"]
    for kx, s in summ.items():
        s["ratio_to_J0_per_n"] = {n: s["per_n"][n]["var_g0_median"] / base[n]["var_g0_median"] for n in base}
    return summ


def main():
    pdbs = P.targets()
    pick = pdbs[::11][:12]
    done = B._done_pdbs(ROWS)
    t0 = time.time()
    for pdb in pick:
        if pdb in done:
            continue
        rows = target_rows(pdb)
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  {pdb} done ({(time.time()-t0)/60:.1f} min)", flush=True)
    rows = B.load_rows(ROWS)
    have = sorted({r["pdb"] for r in rows})
    summ = summarise(rows)
    out = dict(kind="property measurement, gradient variance vs J, no RMSD", lane="S28B",
               targets=have, n_theta=N_THETA, alpha=B.ALPHA, T=B.TEMP, layers=B.LAYERS, ns=NS,
               j_grid=list(B.J_GRID), n_rows=len(rows), summary=summ)
    ST.save_atomic(B.TRAIN_OUT, out, module_file=__file__)
    print(f"\n{'cell':16s} " + " ".join(f"{'n='+str(n):>11}" for n in NS) + f" {'slope':>8} {'slope msq':>9}")
    for kx, s in summ.items():
        print(f"{kx:16s} " + " ".join(f"{s['per_n'][str(n)]['var_g0_median']:11.3e}" for n in NS)
              + f" {s['log2_slope_per_qubit']:+8.3f} {s['log2_slope_msq_per_qubit']:+9.3f}")
    print("\nratio Var(J)/Var(J=0) per n (full objective):")
    for kx, s in summ.items():
        if kx.startswith("full"):
            print(f"{kx:16s} " + " ".join(f"{s['ratio_to_J0_per_n'][str(n)]:11.3f}" for n in NS))
    print("wrote", B.TRAIN_OUT, "targets", len(have))


if __name__ == "__main__":
    main()
