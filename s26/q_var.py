"""s26/q_var.py -- A4: gradient variance of ADAPT-grown circuits beside the fixed ansatz.

PROPERTY MEASUREMENT (s26/PREREG_A4.md).  Reads no RMSD, no native.  Extends
s25/q_plateau.py's width sweep (n = 4..13, layers = 3, theta ~ N(0, 0.6^2), the five
(alpha, T) cells, draws 250/250/250/250/200/120/80) to circuits GROWN by qubit-ADAPT on the
deployed energy shape (standardised ranks 1..2^n) at MATCHED parameter count P = 3n.

For the fixed ansatz the S25 code path is reused unchanged (s25.q_plateau.measure) and, before
anything is extended, the n = 7 row of every cell must reproduce s25/results/q_plateau.json
(same seed 1000 + n, same draws); the reproduction is asserted and recorded.

For the grown circuits: for each (pool, n, cell) run_adapt grows the operator list on
deployed_E(2^n) under that cell (seed 0, 50 Adam steps per growth step, eps 1e-3, P = 3n), then
theta ~ N(0, 0.6^2) is drawn for the GROWN circuit's own P parameters (seed 1000 + n, the S25
law) and Var_theta[dF/dtheta_0] and E|g|^2/P are measured with the exact parameter-shift
gradient.  theta_0 is the first parameter in both ansaetze (the first RY angle).  If ADAPT
stopped early (||g_pool|| < eps) the circuit has fewer than 3n parameters; the count is recorded
and the row is labelled.

Artefact: s26/results/q_var.json, resumable per (pool, n, cell).
"""
from __future__ import annotations

import argparse
import json
import math
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

from s26.q_adapt import PauliCircuit, Pauli, free_energy_grad, run_adapt    # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

CELLS = ((1.0, 0.0), (0.25, 0.0), (0.10, 0.0), (1.0, 0.3), (0.25, 0.3))
CELL_KEYS = {(1.0, 0.0): "linear_alpha1_T0", (0.25, 0.0): "cvar_alpha025_T0",
             (0.10, 0.0): "cvar_alpha01_T0", (1.0, 0.3): "deployed_a1_T03",
             (0.25, 0.3): "deployed_a025_T03"}
NS = (4, 6, 7, 8, 10, 12, 13)
DRAWS = {4: 250, 6: 250, 7: 250, 8: 250, 10: 200, 12: 120, 13: 80}
INIT_SD = 0.6


def _s25():
    from s25 import q_plateau as QP
    return QP


def reproduce_s25_n7():
    """The fixed-ansatz n=7 row of every cell must equal s25/results/q_plateau.json."""
    QP = _s25()
    with open(os.path.join(ROOT, "s25", "results", "q_plateau.json")) as fh:
        ref = json.load(fh)["results"]
    out = {}
    worst = 0.0
    for cell in CELLS:
        alpha, T = cell
        r = QP.measure(7, 3, alpha, T, DRAWS[7], seed=1000 + 7, init_sd=INIT_SD)
        row = [x for x in ref[CELL_KEYS[cell]]["rows"] if x["n"] == 7][0]
        d = abs(r["var_g0"] - row["var_g0"]) / abs(row["var_g0"])
        worst = max(worst, d)
        out[CELL_KEYS[cell]] = dict(here=r["var_g0"], s25=row["var_g0"], rel=d)
        print(f"  reproduce n=7 {CELL_KEYS[cell]:20s} here {r['var_g0']:.6e}  s25 {row['var_g0']:.6e}"
              f"  rel {d:.1e}", flush=True)
    assert worst < 1e-9, worst
    return out, worst


def measure_grown(n, alpha, T, pool, n_theta, seed, adam_steps=50, eps=1e-3, grow_seed=0,
                  optimiser="adam_best"):
    QP = _s25()
    E = QP.deployed_E(1 << n)
    t0 = time.perf_counter()
    r = run_adapt(E, alpha, T, n, pool, max_params=3 * n, adam_steps=adam_steps, seed=grow_seed,
                  eps=eps, record_at=(), optimiser=optimiser)
    ops = [Pauli.from_word(w) for w in r["ops"]]
    circ = PauliCircuit(n, ops)
    P = circ.n_params()
    rng = np.random.default_rng(seed)
    g0, gn = [], []
    for _ in range(n_theta):
        th = rng.normal(0.0, INIT_SD, P)
        _, g, _, _, _ = free_energy_grad(circ, th, E, alpha, T)
        g0.append(float(g[0]))
        gn.append(float(np.dot(g, g)))
    g0 = np.asarray(g0)
    gn = np.asarray(gn)
    return dict(n=n, pool=pool, optimiser=optimiser, P=P, P_matched=(P == 3 * n),
                stopped=r["stopped"], dim=circ.dim, alpha=alpha, T=T,
                var_g0=float(g0.var(ddof=1)), mean_g0=float(g0.mean()),
                mean_sq_norm=float(gn.mean()), mean_sq_per_param=float(gn.mean() / P),
                n_theta=n_theta, sequence=r["sequence"], n_distinct_ops=r["n_distinct_ops"],
                consecutive_repeats=r["consecutive_repeats"], F_grown=r["final"]["F"],
                wall=round(time.perf_counter() - t0, 1))


def slope(rows):
    ns = np.array([r["n"] for r in rows], float)
    v = np.array([r["var_g0"] for r in rows], float)
    ok = v > 0
    return float(np.polyfit(ns[ok], np.log2(v[ok]), 1)[0]) if ok.sum() > 2 else float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pools", default="V,L2")
    ap.add_argument("--optimiser", default="adam_best", choices=("adam_best", "lbfgs", "adam_last"))
    ap.add_argument("--ns", default=",".join(str(x) for x in NS))
    ap.add_argument("--out", default=os.path.join(RESULTS, "q_var.json"))
    ap.add_argument("--skip-fixed", action="store_true", help="do not re-run the fixed sweep")
    a = ap.parse_args(argv)
    from s24 import stats_lib as ST
    QP = _s25()
    R = {}
    if os.path.exists(a.out):
        try:
            with open(a.out) as fh:
                R = json.load(fh).get("results", {})
        except Exception:                                                # noqa: BLE001
            R = {}

    def save():
        ST.save_atomic(a.out, dict(kind="property measurement, gradient variance, no RMSD",
                                   lane="Q", sprint=26, prereg="s26/PREREG_A4.md", results=R),
                       module_file=__file__)

    print("0. REPRODUCE s25/results/q_plateau.json at n = 7 before extending anything")
    if "reproduction" not in R:
        rep, worst = reproduce_s25_n7()
        R["reproduction"] = dict(rows=rep, worst_rel=worst, passed=True)
        save()
    print(f"   worst relative deviation {R['reproduction']['worst_rel']:.2e}")

    ns = [int(x) for x in a.ns.split(",")]
    pools = [x for x in a.pools.split(",") if x]

    print("\n1. FIXED ANSATZ (S25 code path), for the side-by-side table")
    R.setdefault("fixed", {})
    if not a.skip_fixed:
        for cell in CELLS:
            key = CELL_KEYS[cell]
            rows = R["fixed"].get(key, {}).get("rows", [])
            have = {r["n"] for r in rows}
            for n in ns:
                if n in have:
                    continue
                r = QP.measure(n, 3, cell[0], cell[1], DRAWS.get(n, 100), seed=1000 + n, init_sd=INIT_SD)
                rows.append(r)
                R["fixed"][key] = dict(rows=sorted(rows, key=lambda x: x["n"]))
                save()
                print(f"  fixed {key:20s} n={n:2d} P={r['P']:3d} Var[dF/dth0] {r['var_g0']:.6e}", flush=True)
            R["fixed"][key] = dict(rows=sorted(rows, key=lambda x: x["n"]),
                                   log2_slope_per_qubit=slope(rows))
            save()

    print("\n2. ADAPT-GROWN circuits at matched P = 3n, own parameters drawn from the same law")
    R.setdefault("grown", {})
    for pool in pools:
        for cell in CELLS:
            key = f"{pool}:{a.optimiser}:{CELL_KEYS[cell]}"
            rows = R["grown"].get(key, {}).get("rows", [])
            have = {r["n"] for r in rows}
            for n in ns:
                if n in have:
                    continue
                r = measure_grown(n, cell[0], cell[1], pool, DRAWS.get(n, 100), seed=1000 + n,
                                  optimiser=a.optimiser)
                rows.append(r)
                R["grown"][key] = dict(rows=sorted(rows, key=lambda x: x["n"]))
                save()
                print(f"  grown {key:26s} n={n:2d} P={r['P']:3d} ({'matched' if r['P_matched'] else r['stopped']})"
                      f" Var[dF/dth0] {r['var_g0']:.6e}  E|g|^2/P {r['mean_sq_per_param']:.6e}  "
                      f"({r['wall']} s)", flush=True)
            R["grown"][key] = dict(rows=sorted(rows, key=lambda x: x["n"]),
                                   log2_slope_per_qubit=slope(rows))
            save()

    print("\n3. SIDE BY SIDE: fitted log2 Var per qubit")
    tab = {}
    for cell in CELLS:
        key = CELL_KEYS[cell]
        row = {"fixed": R["fixed"].get(key, {}).get("log2_slope_per_qubit")}
        for pool in pools:
            row[f"grown_{pool}_{a.optimiser}"] = R["grown"].get(
                f"{pool}:{a.optimiser}:{key}", {}).get("log2_slope_per_qubit")
        tab[key] = row
        print(f"  {key:20s} " + "  ".join(f"{k} {v if v is None else round(v, 4)}" for k, v in row.items()))
    R["slopes"] = tab
    save()
    print("\nwrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
