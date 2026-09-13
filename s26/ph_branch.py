"""s26/ph_branch.py -- PHYSICS CHOOSES THE PROJECTION BRANCH (lane PH, S26).

Pre-registered in `s26/PREREG_branch_select.md`. Runs only if the tournament ranks
`s26/IDEA_branch_select.md` as a survivor.

    python s26/ph_branch.py solutions [--n N]     # CPU: the five production-rung solutions per target
    python s26/ph_branch.py probe --pdb 1A13      # AMBER: relax the five solutions of one target
    python s26/ph_branch.py relax [--n N]         # AMBER: every target, resumable cells
    python s26/ph_branch.py report                # gated: RMSDs and the statistics

The five solutions are exactly the candidates `core.project.lam_path(multi=True)` chooses among
at the production rung (warm start from the lam=0 solution, plus the four generic starts); the
production choice is reconstructed here and asserted equal to `I.project` (gate G1).
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

from s26 import ph_lib as L                                   # noqa: E402
from s26 import ph_reject as PR                               # noqa: E402

LAM = 0.3
MAXITER = 300
SOL_JSON = os.path.join(L.RESULTS, "ph_branch_solutions.json")
RELAX_JSON = os.path.join(L.RESULTS, "ph_branch_relax.json")
REPORT_JSON = os.path.join(L.RESULTS, "ph_branch_report.json")
SOL_KEYS = ("pdb", "n", "fold", "solutions", "prod_index", "g1_max_abs_dca", "n_distinct")
RELAX_KEYS = ("pdb", "n", "fold", "e0", "e1", "converged", "moved", "wall", "ca_relaxed")


# ------------------------------------------------------------------ pure logic (tested)
def production_index(fvals) -> int:
    """`lam_path(multi=True)`: warm (index 0) unless the best generic start (first strictly
    lowest among indices 1..4) is strictly lower than warm."""
    f = [float(x) for x in fvals]
    best = 1
    for k in range(2, len(f)):
        if f[k] < f[best]:
            best = k
    return best if f[best] < f[0] else 0


def energy_pick(e, converged, fallback: int) -> list:
    """Indices of the lowest converged energy (the whole tie set); `fallback` if none converged."""
    e = np.asarray(e, float); c = np.asarray(converged, bool)
    ok = np.flatnonzero(c & np.isfinite(e))
    if len(ok) == 0:
        return [int(fallback)]
    m = e[ok].min()
    return [int(i) for i in ok if np.isclose(e[i], m, rtol=1e-9, atol=1e-12)]


def distinct_count(cas, tol: float = 1e-3) -> int:
    reps = []
    for ca in cas:
        if all(L.ca_rmsd(ca, r) > tol for r in reps):
            reps.append(ca)
    return len(reps)


# ------------------------------------------------------------------ solutions (CPU)
def solutions_target(t) -> dict:
    from core import project as pj
    from s12 import instrument as I
    pool = PR.load_pool(t["pdb"], oracle=False)
    W, sub, seq, fold = pool["W"], pool["sub"], pool["seq"], pool["fold"]
    C, _ = I.coordinate_average(W[sub])
    n = len(C)
    pen = pj.make_penalty("ramah", seq, int(fold))
    cur = pj.fit_multi(C, pen=None, lam=0.0, maxiter=MAXITER, grad="exact")
    sols = [("warm", pj.fit_prior(C, cur[1], cur[2], pen=pen, lam=LAM, maxiter=MAXITER, grad="exact"))]
    for a, b in pj.STARTS:
        ph = np.full(n, math.radians(a)); ps = np.full(n, math.radians(b))
        sols.append((f"start_{a:g}_{b:g}", pj.fit_prior(C, ph, ps, pen=pen, lam=LAM, maxiter=MAXITER, grad="exact")))
    fv = [s[1][3] for s in sols]
    pi = production_index(fv)
    prod_ca = np.asarray(I.project(C, seq, int(fold))["ca"], float)
    g1 = float(np.abs(np.asarray(sols[pi][1][0], float) - prod_ca).max())
    assert g1 < 1e-9, f"{t['pdb']}: reconstructed production choice differs from I.project by {g1}"
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(fold), "seq": seq, "C": C.tolist(),
            "solutions": [{"name": nm, "phi": np.asarray(s[1], float).tolist(),
                           "psi": np.asarray(s[2], float).tolist(), "ca": np.asarray(s[0], float).tolist(),
                           "fval": float(s[3]), "f_start": float(s[4]), "d_to_C": float(s[5])}
                          for nm, s in sols],
            "prod_index": int(pi), "g1_max_abs_dca": g1,
            "n_distinct": distinct_count([np.asarray(s[1][0], float) for s in sols])}


def solutions(n: int = 0) -> None:
    tg = L.targets()
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_branch_sol")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = solutions_target(t)
        L.write_cell("ph_branch_sol", t["pdb"], r)
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} prod={r['prod_index']} distinct={r['n_distinct']} "
              f"fvals={[round(s['fval'], 4) for s in r['solutions']]} g1={r['g1_max_abs_dca']:.1e} ({clk():.0f}s)",
              flush=True)
    rows = list(L.read_cells("ph_branch_sol").values())
    nd = [r["n_distinct"] for r in rows]
    L.save(SOL_JSON, {"what": "the five production-rung projection solutions per target (native-free)",
                      "lam": LAM, "rows": rows,
                      "summary": {"n": len(rows), "n_distinct_mean": float(np.mean(nd)),
                                  "targets_with_2plus_distinct": int(sum(1 for x in nd if x >= 2)),
                                  "prod_index_counts": np.bincount([r["prod_index"] for r in rows], minlength=5).tolist(),
                                  "g1_max": float(max(r["g1_max_abs_dca"] for r in rows))}},
           rows=rows, complete_keys=SOL_KEYS, n_expected=126, module_file=__file__)


# ------------------------------------------------------------------ relaxation (AMBER)
def relax_target(sol_row) -> dict:
    from s26 import ph_c3 as C3
    seq = sol_row["seq"]
    out = {"pdb": sol_row["pdb"], "n": sol_row["n"], "fold": sol_row["fold"],
           "e0": [], "e1": [], "converged": [], "moved": [], "wall": [], "ca_relaxed": [],
           "clash_out": [], "bond_strain_out": []}
    for s in sol_row["solutions"]:
        r = C3.relax_chain(seq, np.asarray(s["phi"], float), np.asarray(s["psi"], float))
        out["e0"].append(r["e0"]); out["e1"].append(r["e1"]); out["converged"].append(r["converged"])
        out["moved"].append(r["moved"]); out["wall"].append(round(r["wall"], 2))
        out["ca_relaxed"].append(np.asarray(r["ca"], float).tolist())
        out["clash_out"].append(r["panel_out"]["n_clash_2A"]); out["bond_strain_out"].append(r["panel_out"]["bond_strain"])
    return out


def probe(pdb: str) -> dict:
    sols = L.read_cells("ph_branch_sol")
    if pdb not in sols:
        t = {x["pdb"]: x for x in L.targets()}[pdb]
        row = solutions_target(t)
        L.write_cell("ph_branch_sol", pdb, row)
    else:
        row = sols[pdb]
    t0 = time.time()
    r = relax_target(row)
    print(json.dumps({"pdb": pdb, "n_distinct": row["n_distinct"], "prod_index": row["prod_index"],
                      "e0": [round(x, 1) for x in r["e0"]], "e1": [round(x, 1) for x in r["e1"]],
                      "converged": r["converged"], "wall_per_solution": r["wall"],
                      "wall_total": round(time.time() - t0, 1)}, indent=1))
    return r


def relax(n: int = 0) -> None:
    sols = L.read_cells("ph_branch_sol")
    tg = [t for t in L.targets() if t["pdb"] in sols]
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_branch")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = relax_target(sols[t["pdb"]])
        L.write_cell("ph_branch", t["pdb"], r)
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} e1={[round(x, 1) for x in r['e1']]} conv={sum(r['converged'])}/5 "
              f"wall={sum(r['wall']):.0f}s ({clk():.0f}s)", flush=True)
    rows = list(L.read_cells("ph_branch").values())
    L.save(RELAX_JSON, {"what": "production relaxation of the five projection solutions per target (native-free)",
                        "rows": rows}, rows=rows, complete_keys=RELAX_KEYS, n_expected=126, module_file=__file__)


# ------------------------------------------------------------------ report (gated)
def report() -> dict:
    L.require_gate("ph_branch report")
    from s24 import stats_lib as ST
    sols = L.read_cells("ph_branch_sol"); rel = L.read_cells("ph_branch")
    pdbs = sorted(p for p in sols if p in rel)
    folds = L.folds_of(pdbs)
    rows = []
    for p in pdbs:
        u = L.univ_oracle(p); nat = u["nat_ca"]
        s, r = sols[p], rel[p]
        cas = [np.asarray(x["ca"], float) for x in s["solutions"]]
        d_built = np.array([L.ca_rmsd(ca, nat) for ca in cas])
        d_relaxed = np.array([L.ca_rmsd(np.asarray(ca, float), nat) for ca in r["ca_relaxed"]])
        pi = s["prod_index"]
        pick_e1 = energy_pick(r["e1"], r["converged"], pi)
        pick_e0 = energy_pick(r["e0"], [True] * 5, pi)
        rows.append({"pdb": p, "n_distinct": s["n_distinct"], "prod": float(d_built[pi]),
                     "prod_relaxed": float(d_relaxed[pi]),
                     "pick_e1": float(np.mean(d_built[pick_e1])), "pick_e1_relaxed": float(np.mean(d_relaxed[pick_e1])),
                     "pick_e1_ties": len(pick_e1), "pick_e0": float(np.mean(d_built[pick_e0])),
                     "random": float(d_built.mean()), "oracle_min": float(d_built.min()),
                     "d_built": d_built.tolist(), "n_converged": int(sum(r["converged"])),
                     "e1_pick_is_oracle": bool(np.isclose(np.mean(d_built[pick_e1]), d_built.min())),
                     "prod_is_oracle": bool(np.isclose(d_built[pi], d_built.min()))})
    A = {k: np.array([x[k] for x in rows], float) for k in ("prod", "pick_e1", "pick_e0", "random", "oracle_min", "pick_e1_relaxed", "prod_relaxed")}
    moved = np.array([x["n_distinct"] >= 2 for x in rows])
    out = {"n": len(rows), "n_moved": int(moved.sum()),
           "frac_e1_pick_is_oracle": float(np.mean([x["e1_pick_is_oracle"] for x in rows])),
           "frac_prod_is_oracle": float(np.mean([x["prod_is_oracle"] for x in rows]))}
    for lab, a, b in (("pick_e1 minus production (built chain)", A["pick_e1"], A["prod"]),
                      ("pick_e1 minus random pick (built chain)", A["pick_e1"], A["random"]),
                      ("pick_e0 minus production (built chain)", A["pick_e0"], A["prod"]),
                      ("random pick minus production (built chain)", A["random"], A["prod"]),
                      ("pick_e1 minus production (relaxed chain)", A["pick_e1_relaxed"], A["prod_relaxed"])):
        for sub_name, m in (("all", np.ones(len(rows), bool)), ("moved", moved)):
            if m.sum() < 12:
                continue
            c = ST.compare(a[m], b[m], folds[m], names=[p for p, o in zip(pdbs, m) if o],
                           label=f"[{sub_name}, n={int(m.sum())}] {lab}")
            print(ST.fmt(c))
            out[f"{lab} [{sub_name}]"] = c
    M = np.array([x["d_built"] for x in rows], float)
    b = ST.best_of_k_within(M)
    out["oracle_over_5"] = b
    print(f"  ORACLE min over the five: {b['observed_gain']:+.4f}, split-half transfer {b['split_half']:+.4f} "
          f"({100 * b['split_half_frac']:.0f}%), k_eff {b['k_eff']:.2f} -> {b['verdict']}")
    L.save(REPORT_JSON, {"what": "branch selection by the relaxed energy vs the objective", "rows": rows, "report": out},
           rows=rows, complete_keys=("pdb", "prod", "pick_e1", "random"), n_expected=126, module_file=__file__)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("solutions", "probe", "relax", "report"))
    ap.add_argument("--pdb", default="1A13")
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "solutions":
        solutions(a.n)
    elif a.mode == "probe":
        probe(a.pdb)
    elif a.mode == "relax":
        relax(a.n)
    else:
        report()
