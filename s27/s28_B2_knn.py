#!/usr/bin/env python
"""s27/s28_B2_knn.py -- S28 LANE B2: the hopping graph with a SPREAD spectrum.

Pre-registered in `s27/PREREG_S28_B.md` ADDENDUM 1 (read it first). Everything but the graph is
`s27/s28_B_hop.py`: the objective, the exact parameter-shift hopping gradient, the J = 0 loop,
the eigensolver, the three readouts, the controls, the trainability measurement.

THE GRAPH
=========
Symmetric k-nearest-neighbour graph on the pool's pairwise CA-RMSD (i ~ j if j is among i's k
nearest or i among j's; binary), degree-normalised D^-1/2 A D^-1/2, scaled to unit spectral
norm (identity for a connected graph, kept so there is one code path), zero diagonal, zero
padding rows. Its spectrum is that of a sparse graph, which is the one respect in which
S28-L8b/S28-L11's mechanism (the Gaussian graph is a near-rank-one typicality projector) names a
different Hamiltonian.

USAGE
=====
    python s27/s28_B2_knn.py --train            F5-B2: gradient variance vs J at n = 4..9, k in {5, 10}
    python s27/s28_B2_knn.py --run --k 10       point-cloud endpoint (gated: only after the S28B chain verdict)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s22 import qcand_lib as QC            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

RESULTS = B.RESULTS
TRAIN_ROWS = os.path.join(RESULTS, "s28_B2_train_rows.jsonl")
TRAIN_OUT = os.path.join(RESULTS, "s28_B2_train.json")
KS = (5, 10)
J_GRID_B2 = (0.0, 0.3, 1.0, 3.0)
NS = [4, 5, 6, 7, 8, 9]
N_THETA = 120


# ======================================================================== the kNN graph
def knn_graph(D: np.ndarray, k: int) -> Dict:
    """Symmetric binary kNN graph on a distance matrix, degree-normalised, unit spectral norm.

    Returns the same keys as `s28_B_hop.finish_graph` plus the kNN report: the binary
    adjacency `A_bin`, degrees (>= k after symmetrisation), connected components, the top
    eigenvector's participation ratio / dim and its overlap with the uniform state.
    """
    D = np.asarray(D, float)
    n = len(D)
    if not (1 <= k < n):
        raise ValueError(f"k must be in [1, n-1], got {k} for n={n}")
    Dm = D.copy()
    np.fill_diagonal(Dm, np.inf)
    nn = np.argsort(Dm, axis=1, kind="stable")[:, :k]         # k nearest of each row
    A_bin = np.zeros((n, n))
    rows = np.repeat(np.arange(n), k)
    A_bin[rows, nn.ravel()] = 1.0
    A_bin = np.maximum(A_bin, A_bin.T)                         # symmetric: i ~ j if either lists the other
    np.fill_diagonal(A_bin, 0.0)
    deg = A_bin.sum(1)
    dinv = 1.0 / np.sqrt(np.maximum(deg, 1e-300))
    A_norm = (dinv[:, None] * A_bin) * dinv[None, :]
    g = B.finish_graph(A_norm, sigma=float("nan"))            # eigvalsh, unit spectral norm, degree of A_norm
    w, v = np.linalg.eigh(g["A"])
    v1 = v[:, -1] * np.sign(v[:, -1].sum() or 1.0)
    u = np.full(n, 1.0 / np.sqrt(n))
    # connected components by BFS on the binary graph
    seen = np.zeros(n, bool)
    ncomp = 0
    for s in range(n):
        if seen[s]:
            continue
        ncomp += 1
        stack = [s]
        seen[s] = True
        while stack:
            x = stack.pop()
            for y in np.flatnonzero(A_bin[x] > 0):
                if not seen[y]:
                    seen[y] = True
                    stack.append(y)
    g.update(A_bin=A_bin, k=int(k), knn_degree=deg, knn_degree_min=int(deg.min()),
             knn_degree_max=int(deg.max()), n_components=int(ncomp),
             lam2_over_lam1=float(w[-2] / w[-1]), perron_pr_over_dim=float(1.0 / np.sum(v1 ** 4) / n),
             perron_uniform_overlap=float((u @ v1) ** 2), density=float(A_bin.sum() / (n * (n - 1))))
    return g


def spectral_report(D: np.ndarray, k: int) -> Dict:
    """The kNN graph's spectral numbers beside the Gaussian graph's, one target."""
    gk = knn_graph(D, k)
    gg = B.kernel_graph(D)
    w, v = np.linalg.eigh(gg["A"])
    v1 = v[:, -1] * np.sign(v[:, -1].sum() or 1.0)
    n = len(D)
    u = np.full(n, 1.0 / np.sqrt(n))
    return dict(k=k,
                knn=dict(lam2_over_lam1=gk["lam2_over_lam1"], perron_pr_over_dim=gk["perron_pr_over_dim"],
                         perron_uniform_overlap=gk["perron_uniform_overlap"], n_components=gk["n_components"],
                         degree_min=gk["knn_degree_min"], degree_max=gk["knn_degree_max"], density=gk["density"]),
                gauss=dict(lam2_over_lam1=float(w[-2] / w[-1]), perron_pr_over_dim=float(1.0 / np.sum(v1 ** 4) / n),
                           perron_uniform_overlap=float((u @ v1) ** 2), sigma=gg["sigma"]))


# ================================================================ F5-B2: trainability
def target_rows(pdb: str, cand, E_full: np.ndarray, order: np.ndarray, D: np.ndarray) -> List[Dict]:
    rows = []
    for k in KS:
        for n in NS:
            dim = 1 << n
            if dim <= cand.k:
                sub = order[:dim]
                g = knn_graph(D[np.ix_(sub, sub)], k)
                A = g["A"]
                E = B.deployed_E(dim)
            else:
                g = knn_graph(D, k)
                A = B.pad_graph(g["A"], dim)
                E = QC.Encoding(E_full).E
            info = dict(pdb=pdb, k=k, lam2_over_lam1=g["lam2_over_lam1"], perron_pr_over_dim=g["perron_pr_over_dim"],
                        perron_uniform_overlap=g["perron_uniform_overlap"], n_components=g["n_components"],
                        knn_degree_min=g["knn_degree_min"], knn_degree_max=g["knn_degree_max"],
                        deg_E_corr=float(np.corrcoef(g["knn_degree"], E[:len(g["knn_degree"])])[0, 1]))
            for J in J_GRID_B2:
                r = B.measure_hop(n, B.LAYERS, B.ALPHA, B.TEMP, J, N_THETA, 1009, 0.6, E, A, "full")
                r.update(info); rows.append(r)
            r = B.measure_hop(n, B.LAYERS, B.ALPHA, B.TEMP, 1.0, N_THETA, 1009, 0.6, E, A, "hop_only")
            r.update(info); rows.append(r)
    return rows


def summarise(rows: List[Dict]) -> Dict:
    summ = {}
    cells = sorted({(r["k"], r["mode"], r["J"]) for r in rows})
    for k, mode, J in cells:
        per_n = {}
        for n in NS:
            v = [r["var_g0"] for r in rows if r["k"] == k and r["mode"] == mode and r["J"] == J and r["n"] == n]
            m = [r["mean_sq_per_param"] for r in rows if r["k"] == k and r["mode"] == mode and r["J"] == J and r["n"] == n]
            per_n[str(n)] = dict(var_g0_median=float(np.median(v)), var_g0_min=float(min(v)), var_g0_max=float(max(v)),
                                 msq_median=float(np.median(m)), n_targets=len(v))
        ns_ = np.array(NS, float)
        v = np.array([per_n[str(n)]["var_g0_median"] for n in NS])
        ok = v > 0
        slope = float(np.polyfit(ns_[ok], np.log2(v[ok]), 1)[0]) if ok.sum() > 2 else float("nan")
        slope48 = float(np.polyfit(ns_[:5][ok[:5]], np.log2(v[:5][ok[:5]]), 1)[0]) if ok[:5].sum() > 2 else float("nan")
        summ[f"k{k}|{mode}|J{J:g}"] = dict(per_n=per_n, log2_slope_per_qubit=slope, log2_slope_4to8=slope48)
    for k in KS:
        base = summ[f"k{k}|full|J0"]["per_n"]
        for kx, s in summ.items():
            if kx.startswith(f"k{k}|"):
                s["ratio_to_J0_per_n"] = {n: s["per_n"][n]["var_g0_median"] / base[n]["var_g0_median"] for n in base}
    spec = {}
    for k in KS:
        for n in NS:
            rs = [r for r in rows if r["k"] == k and r["n"] == n and r["mode"] == "hop_only"]
            spec[f"k{k}|n{n}"] = {key: float(np.median([r[key] for r in rs])) for key in
                                  ("lam2_over_lam1", "perron_pr_over_dim", "perron_uniform_overlap", "n_components",
                                   "knn_degree_min", "knn_degree_max", "deg_E_corr")}
    return dict(summary=summ, spectral=spec)


def train_main():
    from s25 import phys_lib as P
    from s27 import run_pool as RP
    pdbs = P.targets()
    pick = pdbs[::11][:12]
    done = B._done_pdbs(TRAIN_ROWS)
    t0 = time.time()
    for pdb in pick:
        if pdb in done:
            continue
        cand, ch, _ = RP.channels_for(pdb)
        E_full = RP.zr(ch["DIS"])
        key = RP.rng_for(pdb, "tiekey").random(cand.k)
        order = RP.topm(E_full, cand.k, key)
        D = B.pairwise_rmsd_matrix(cand.W)
        rows = target_rows(pdb, cand, E_full, order, D)
        with open(TRAIN_ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  {pdb} done ({(time.time()-t0)/60:.1f} min)", flush=True)
    rows = B.load_rows(TRAIN_ROWS)
    S = summarise(rows)
    gauss = json.load(open(B.TRAIN_OUT, encoding="utf-8"))["summary"] if os.path.exists(B.TRAIN_OUT) else {}
    out = dict(kind="property measurement, gradient variance vs J on the kNN graph, no RMSD", lane="S28B2",
               targets=sorted({r["pdb"] for r in rows}), n_theta=N_THETA, ks=list(KS), ns=NS, j_grid=list(J_GRID_B2),
               n_rows=len(rows), gaussian_reference=dict(
                   hop_only=gauss.get("hop_only|J1", {}).get("per_n"), full_J0=gauss.get("full|J0", {}).get("per_n"),
                   full_J3=gauss.get("full|J3", {}).get("per_n"),
                   slopes={kx: v.get("log2_slope_per_qubit") for kx, v in gauss.items()}), **S)
    ST.save_atomic(TRAIN_OUT, out, module_file=__file__)
    print(f"\n{'cell':20s} " + " ".join(f"{'n='+str(n):>11}" for n in NS) + f" {'slope':>8} {'4..8':>8}")
    for kx, s in S["summary"].items():
        print(f"{kx:20s} " + " ".join(f"{s['per_n'][str(n)]['var_g0_median']:11.3e}" for n in NS)
              + f" {s['log2_slope_per_qubit']:+8.3f} {s['log2_slope_4to8']:+8.3f}")
    if gauss:
        g = gauss["hop_only|J1"]["per_n"]
        print(f"{'gauss hop_only J1':20s} " + " ".join(f"{g[str(n)]['var_g0_median']:11.3e}" for n in NS)
              + f" {gauss['hop_only|J1']['log2_slope_per_qubit']:+8.3f}")
    print("\nspectral (median over targets):")
    for kx, s in S["spectral"].items():
        print(f"  {kx:8s} lam2/lam1 {s['lam2_over_lam1']:.3f}  PR/dim {s['perron_pr_over_dim']:.3f}  |<u|v1>|^2 {s['perron_uniform_overlap']:.3f}  comps {s['n_components']:.0f}  deg {s['knn_degree_min']:.0f}..{s['knn_degree_max']:.0f}  corr(deg,E) {s['deg_E_corr']:+.3f}")
    print("wrote", TRAIN_OUT)


# ================================================================ endpoint (gated)
def run_main(k: int, limit: int = 0):
    """Point-cloud endpoint with the kNN graph. GATED on the S28B built-chain verdict entry."""
    from s25 import phys_lib as P
    from s27 import run_pool as RP
    rows_path = os.path.join(RESULTS, f"s28_B2_rows_k{k}.jsonl")
    real_kernel = B.kernel_graph

    def knn_kernel(D, sigma=None):                       # the graph switch
        return knn_graph(D, k)

    B.kernel_graph = knn_kernel
    B.J_GRID = J_GRID_B2
    done = B._done_pdbs(rows_path)
    pdbs = P.targets()[:limit] if limit else P.targets()
    t0 = time.time()
    try:
        for i, pdb in enumerate(pdbs):
            if pdb in done:
                continue
            t1 = time.time()
            rows = B.run_target(pdb)
            for r in rows:
                r["graph_kind"] = f"knn{k}"
            with open(rows_path, "a", encoding="utf-8") as fh:
                for r in rows:
                    fh.write(json.dumps(r) + "\n")
            j0 = [r for r in rows if r["source"] == "vqe" and r["J"] == 0.0 and r["seed"] == 0][0]
            print(f"  [{i+1}/{len(pdbs)}] {pdb} rows={len(rows)} vqe0={j0['R1']['rmsd']:.3f} "
                  f"{time.time()-t1:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    finally:
        B.kernel_graph = real_kernel
    print("done:", rows_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.train:
        train_main()
    if a.run:
        run_main(a.k, a.limit)


if __name__ == "__main__":
    main()
