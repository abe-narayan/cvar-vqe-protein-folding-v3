"""s26/q_var_boot.py -- bootstrap CIs on the A4 slopes (PREREG_A4 addendum 2, the L47 caveat).

Re-measures every (cell, n) row of s26/results/q_var.json with the SAME theta draws (seed
1000 + n, N(0, 0.6^2), the same counts) but stores every per-draw dF/dtheta_0; asserts the
recomputed variance equals q_var.json's var_g0 to relative 1e-12 (fixed rows) and that the
re-grown circuits reproduce their stored operator sequence exactly (grown rows); then a
percentile bootstrap over the draws within each n (2000 resamples, seed 2026) for every
slope and for grown - fixed per cell.  Property measurement, no native, no score.

    python s26/q_var_boot.py    -> s26/results/q_var_boot.json (resumable per row)
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

from core import quantum as Q                                              # noqa: E402
from s26.q_adapt import PauliCircuit, Pauli, free_energy_grad, run_adapt   # noqa: E402
from s26.q_var import CELLS, CELL_KEYS, DRAWS, INIT_SD                     # noqa: E402

SRC = os.path.join(HERE, "results", "q_var.json")
OUT = os.path.join(HERE, "results", "q_var_boot.json")
NBOOT, BSEED = 2000, 2026


def deployed_E(dim):
    r = np.arange(1, dim + 1, dtype=float)
    return (r - r.mean()) / r.std()


def draws_fixed(n, alpha, T, n_theta, seed):
    """Per-draw dF/dtheta_0 for the fixed ansatz, the S25 draw order exactly."""
    circ = Q.StatevectorCircuit(n, 3)
    E = deployed_E(circ.dim)
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    g0 = np.empty(n_theta)
    for i in range(n_theta):
        th = rng.normal(0.0, INIT_SD, P)
        if T == 0.0:
            g = Q.grad_cvar_paramshift(circ, th, E, alpha)
        else:
            _, g, _, _, _ = Q.free_energy(circ, th, E, alpha, T)
        g0[i] = float(g[0])
    return g0


def draws_grown(n, alpha, T, pool, n_theta, seed, expect_seq):
    E = deployed_E(1 << n)
    r = run_adapt(E, alpha, T, n, pool, max_params=3 * n, adam_steps=50, seed=0, eps=1e-3,
                  record_at=(), optimiser="adam_best")
    assert r["sequence"] == expect_seq, (n, alpha, T, pool, r["sequence"], expect_seq)
    circ = PauliCircuit(n, [Pauli.from_word(w) for w in r["ops"]])
    rng = np.random.default_rng(seed)
    P = circ.n_params()
    g0 = np.empty(n_theta)
    for i in range(n_theta):
        th = rng.normal(0.0, INIT_SD, P)
        _, g, _, _, _ = free_energy_grad(circ, th, E, alpha, T)
        g0[i] = float(g[0])
    return g0, P


def slope_of(vars_by_n):
    ns = np.array(sorted(vars_by_n), float)
    v = np.array([vars_by_n[int(k)] for k in ns])
    ok = v > 0
    if ok.sum() < 3:
        return float("nan")
    return float(np.polyfit(ns[ok], np.log2(v[ok]), 1)[0])


def boot_slopes(rows, rng):
    """rows: {n: g0 array}.  Returns NBOOT slope resamples."""
    out = np.empty(NBOOT)
    for b in range(NBOOT):
        vb = {}
        for n, g in rows.items():
            idx = rng.integers(0, len(g), len(g))
            vb[n] = float(np.var(g[idx], ddof=1))
        out[b] = slope_of(vb)
    return out


def main():
    from s24 import stats_lib as ST
    J = json.load(open(SRC))["results"]
    R = {"fixed": {}, "grown": {}, "repro_worst_rel": 0.0}
    if os.path.exists(OUT):
        try:
            R = json.load(open(OUT))["results"]
        except Exception:                                                # noqa: BLE001
            pass

    def save():
        ST.save_atomic(OUT, dict(kind="property, per-draw gradients and bootstrap slope CIs, no RMSD",
                                 lane="Q", sprint=26, prereg="s26/PREREG_A4.md addendum 2",
                                 nboot=NBOOT, boot_seed=BSEED, results=R), module_file=__file__)

    t0 = time.perf_counter()
    print("1. FIXED rows: per-draw g0, reproduction of q_var.json var_g0 asserted")
    for cell in CELLS:
        key = CELL_KEYS[cell]
        R["fixed"].setdefault(key, {})
        for row in J["fixed"][key]["rows"]:
            n = int(row["n"])
            if str(n) in R["fixed"][key]:
                continue
            g0 = draws_fixed(n, cell[0], cell[1], int(row["n_theta"]), seed=1000 + n)
            v = float(g0.var(ddof=1))
            rel = abs(v - row["var_g0"]) / abs(row["var_g0"])
            assert rel < 1e-12, (key, n, v, row["var_g0"], rel)
            R["repro_worst_rel"] = max(R["repro_worst_rel"], rel)
            R["fixed"][key][str(n)] = dict(n=n, P=row["P"], var_g0=v, g0=g0.tolist())
            save()
            print(f"  fixed {key:20s} n={n:2d} var {v:.6e} (rel {rel:.1e})  {time.perf_counter() - t0:.0f} s", flush=True)
    print("\n2. GROWN rows: re-grown (seed 0), sequence asserted, per-draw g0")
    for gkey, block in J["grown"].items():
        pool, opt, key = gkey.split(":")
        if opt != "adam_best":
            continue
        R["grown"].setdefault(gkey, {})
        for row in block["rows"]:
            n = int(row["n"])
            if str(n) in R["grown"][gkey]:
                continue
            g0, P = draws_grown(n, row["alpha"], row["T"], pool, int(row["n_theta"]), seed=1000 + n,
                                expect_seq=row["sequence"])
            v = float(g0.var(ddof=1))
            rel = abs(v - row["var_g0"]) / max(abs(row["var_g0"]), 1e-300)
            R["grown"][gkey][str(n)] = dict(n=n, P=P, P_matched=(P == 3 * n), var_g0=v,
                                            var_q_var=row["var_g0"], rel=rel, g0=g0.tolist())
            save()
            print(f"  grown {gkey:32s} n={n:2d} P={P:3d} var {v:.6e} (rel to q_var {rel:.1e})  "
                  f"{time.perf_counter() - t0:.0f} s", flush=True)
    print("\n3. BOOTSTRAP over the theta draws, within each n")
    rng = np.random.default_rng(BSEED)
    ci = {}
    for cell in CELLS:
        key = CELL_KEYS[cell]
        fixed_rows = {int(k): np.asarray(v["g0"]) for k, v in R["fixed"][key].items()}
        bf = boot_slopes(fixed_rows, rng)
        ci[key] = {"fixed": dict(point=slope_of({n: float(np.var(g, ddof=1)) for n, g in fixed_rows.items()}),
                                 ci95=[float(np.percentile(bf, 2.5)), float(np.percentile(bf, 97.5))],
                                 n_rows=len(fixed_rows))}
        for pool in ("V", "L2"):
            gkey = f"{pool}:adam_best:{key}"
            if gkey not in R["grown"]:
                continue
            rows = {int(k): np.asarray(v["g0"]) for k, v in R["grown"][gkey].items() if v["P_matched"]}
            if len(rows) < 3:
                ci[key][f"grown_{pool}"] = dict(point=float("nan"), ci95=None, n_rows=len(rows), note="degenerate")
                continue
            bg = boot_slopes(rows, rng)
            d = bg - bf
            ci[key][f"grown_{pool}"] = dict(
                point=slope_of({n: float(np.var(g, ddof=1)) for n, g in rows.items()}),
                ci95=[float(np.percentile(bg, 2.5)), float(np.percentile(bg, 97.5))], n_rows=len(rows),
                diff_vs_fixed=dict(point=float(np.mean(bg) - np.mean(bf)),
                                   ci95=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
                                   includes_zero=bool(np.percentile(d, 2.5) <= 0.0 <= np.percentile(d, 97.5))))
        line = f"  {key:20s} fixed {ci[key]['fixed']['point']:+.3f} [{ci[key]['fixed']['ci95'][0]:+.3f},{ci[key]['fixed']['ci95'][1]:+.3f}]"
        for pool in ("V", "L2"):
            g = ci[key].get(f"grown_{pool}")
            if g and g.get("ci95"):
                line += (f"   grown {pool} {g['point']:+.3f} [{g['ci95'][0]:+.3f},{g['ci95'][1]:+.3f}]"
                         f"  diff [{g['diff_vs_fixed']['ci95'][0]:+.3f},{g['diff_vs_fixed']['ci95'][1]:+.3f}]"
                         f"{' incl 0' if g['diff_vs_fixed']['includes_zero'] else ' EXCL 0'}")
            elif g:
                line += f"   grown {pool} degenerate"
        print(line, flush=True)
    R["ci"] = ci
    save()
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
