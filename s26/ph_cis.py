"""s26/ph_cis.py -- THE CIS-PEPTIDE GAP (lane PH, S26).  Pre-registered in s26/PREREG_cis.md.

PART 1 -- THE CENSUS (allowed before the phase gate, ledger L5: reads native OMEGA angles and
pool CA-CA distances, computes no RMSD, selects nothing).  ORACLE DIAGNOSTIC where it reads a
native; native-free where it reads the window bank.

    For each of the 126 dev targets (`s12.instrument.targets()`, natives under `pdbs_ext/`,
    model 1 = the model the instrument scores against):
      omega_i = dihedral(CA_i, C_i, N_{i+1}, CA_{i+1}),  cis iff |omega| < 30 deg
      cross-check: consecutive CA-CA < 3.3 A (cis ~2.9 A, trans ~3.80 A)
      the residue after the bond (proline or not), the deviation |180 - |omega|| per bond
      the same census over EVERY deposited model (NMR ensembles), reported separately
    For the retrieval pool: consecutive CA-CA < 3.3 A over every window of the universe, the
    shipped K=500 pool and the production top-75.

    WHY THE ANSWER MAY BE ZERO BY CONSTRUCTION, stated before the run.  `core/data.py:406`
    drops any peptide whose consecutive CA-CA step is below 3.5 A or above 4.1 A, and
    `core/data.py:698` drops any fragment window containing such a step.  If those filters
    did their job, model-1 natives and pool windows carry no cis bond at all and the census
    measures the FILTER, not the projection.  The census is run anyway, because "the code
    says so" is not a measurement, because the filter reads CA-CA and not omega (a non-planar
    omega between 30 and 150 deg can pass a CA-CA gate), and because the ensemble census and
    the omega-deviation distribution are what the projection-floor question needs.

PART 2 -- THE PROJECTION FLOOR (after the gate).  `core.data.Peptide.rebuild` is, read from
`core/data.py:410-419`, the backbone (N, CA, C stacked) RMSD between the native and the chain
rebuilt from its own phi/psi through `core.geometry.build_backbone` at ideal TRANS geometry,
after Kabsch superposition.  It is exactly the floor of the production projection's
representation on the native's own torsions, on the N/CA/C basis.  Part 2 also computes the
CA-only floor (the basis the endpoint uses), the production built-chain cost
`rmsd_arm - rmsd_avg` per target, and their relation to the omega deviation.

    python s26/ph_cis.py census          # part 1, native-free w.r.t. RMSD
    python s26/ph_cis.py floor           # part 2, gated
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s26 import ph_lib as L                                   # noqa: E402

CENSUS_JSON = os.path.join(L.RESULTS, "ph_cis_census.json")
FLOOR_JSON = os.path.join(L.RESULTS, "ph_cis_floor.json")

CENSUS_KEYS = ("pdb", "n", "fold", "seq_ok", "n_bonds", "omega_deg", "n_cis_omega30",
               "n_cis_ca33", "cis_bonds", "max_omega_dev", "mean_omega_dev",
               "n_nonplanar30", "ens_n_models", "ens_n_models_with_cis",
               "univ_n_windows", "univ_n_cis_windows", "univ_min_step",
               "pool_n_cis_windows", "pool_min_step", "top75_n_cis_windows",
               "top75_min_step")


# ------------------------------------------------------------------ part 1
def census_target(t) -> dict:
    pdb, n, seq = t["pdb"], int(t["n"]), t["seq"]
    nseq, N, CA, C = L.native_backbone(pdb, model_index=0)
    om = L.omega_deg(CA, C, N)
    step = L.consecutive_ca(CA)
    dev = np.abs(180.0 - np.abs(om))
    cis_om = np.abs(om) < L.CIS_OMEGA_DEG
    cis_ca = step < L.CIS_CA_CA
    bonds = []
    for i in np.flatnonzero(cis_om | cis_ca):
        bonds.append({"bond": int(i), "res_after": nseq[i + 1] if i + 1 < len(nseq) else "?",
                      "is_pro": bool(i + 1 < len(nseq) and nseq[i + 1] == "P"),
                      "omega_deg": float(om[i]), "ca_ca": float(step[i]),
                      "cis_by_omega": bool(cis_om[i]), "cis_by_ca": bool(cis_ca[i])})
    # every deposited model, reported separately (the instrument scores model 1 only)
    ens = L.native_ensemble(pdb)
    ens_cis = 0
    for _s, eN, eCA, eC in ens:
        eo = L.omega_deg(eCA, eC, eN)
        if (np.abs(eo) < L.CIS_OMEGA_DEG).any():
            ens_cis += 1
    # the retrieval pool: CA only, native-free
    u = L.univ_nativefree(pdb, keys=("W", "order"))
    W = u["W"]
    st = L.consecutive_ca(W)                     # (nw, n-1)
    wcis = (st < L.CIS_CA_CA).any(1)
    p = L.pool_idx_from_order(u)
    rec = L.prod_record_nativefree(pdb)
    sub = np.asarray(rec["sub"], int)
    top = p[sub]
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]),
            "seq_ok": bool(nseq == seq and len(CA) == n),
            "n_bonds": int(len(om)), "omega_deg": om.tolist(),
            "n_cis_omega30": int(cis_om.sum()), "n_cis_ca33": int(cis_ca.sum()),
            "cis_bonds": bonds,
            "max_omega_dev": float(dev.max()) if len(dev) else 0.0,
            "mean_omega_dev": float(dev.mean()) if len(dev) else 0.0,
            "n_nonplanar30": int((dev > 30.0).sum()),
            "n_pro": int(seq.count("P")),
            "ens_n_models": int(len(ens)), "ens_n_models_with_cis": int(ens_cis),
            "univ_n_windows": int(len(W)), "univ_n_cis_windows": int(wcis.sum()),
            "univ_min_step": float(st.min()),
            "pool_n_cis_windows": int(wcis[p].sum()), "pool_min_step": float(st[p].min()),
            "top75_n_cis_windows": int(wcis[top].sum()),
            "top75_min_step": float(st[top].min())}


def census() -> dict:
    tg = L.targets()
    rows = []
    clk = L.Clock()
    for k, t in enumerate(tg):
        r = census_target(t)
        rows.append(r)
        print(f"[{k + 1}/{len(tg)}] {r['pdb']} n={r['n']} cis(omega<30)={r['n_cis_omega30']} "
              f"cis(CA<3.3)={r['n_cis_ca33']} maxdev={r['max_omega_dev']:.1f} "
              f"ens {r['ens_n_models_with_cis']}/{r['ens_n_models']} "
              f"pool_cis={r['pool_n_cis_windows']} ({clk():.0f}s)", flush=True)
    summ = summarise_census(rows)
    L.save(CENSUS_JSON, {"what": "cis-peptide census, 126 dev natives (model 1) + ensembles + "
                                 "retrieval pool; ORACLE DIAGNOSTIC on the natives, native-free "
                                 "on the pool; no RMSD computed",
                         "criteria": {"cis_omega_deg": L.CIS_OMEGA_DEG, "cis_ca_ca": L.CIS_CA_CA},
                         "rows": rows, "summary": summ},
           rows=rows, complete_keys=CENSUS_KEYS, n_expected=126, module_file=__file__)
    print(json.dumps(summ, indent=1))
    return summ


def summarise_census(rows) -> dict:
    n_cis_om = sum(1 for r in rows if r["n_cis_omega30"] > 0)
    n_cis_ca = sum(1 for r in rows if r["n_cis_ca33"] > 0)
    agree = sum(1 for r in rows if (r["n_cis_omega30"] > 0) == (r["n_cis_ca33"] > 0))
    bonds = [b for r in rows for b in r["cis_bonds"]]
    n_pro = sum(1 for b in bonds if b["is_pro"])
    devs = np.concatenate([np.abs(180.0 - np.abs(np.asarray(r["omega_deg"]))) for r in rows])
    ens_tot = sum(r["ens_n_models"] for r in rows)
    ens_cis = sum(r["ens_n_models_with_cis"] for r in rows)
    return {
        "n_targets": len(rows),
        "n_seq_mismatch": sum(1 for r in rows if not r["seq_ok"]),
        "targets_with_cis_omega30": n_cis_om,
        "targets_with_cis_ca33": n_cis_ca,
        "targets_criteria_agree": agree,
        "cis_bonds_total": len(bonds), "cis_bonds_proline": n_pro,
        "cis_targets": sorted({r["pdb"] for r in rows if r["n_cis_omega30"] > 0}),
        "omega_dev_deg": {"mean": float(devs.mean()), "median": float(np.median(devs)),
                          "p90": float(np.percentile(devs, 90)),
                          "p99": float(np.percentile(devs, 99)), "max": float(devs.max()),
                          "n_bonds": int(len(devs)),
                          "frac_over_20": float((devs > 20).mean()),
                          "frac_over_30": float((devs > 30).mean())},
        "targets_with_nonplanar30": sum(1 for r in rows if r["n_nonplanar30"] > 0),
        "ensemble": {"n_models_total": ens_tot, "n_models_with_cis": ens_cis,
                     "targets_with_cis_in_any_model":
                         sum(1 for r in rows if r["ens_n_models_with_cis"] > 0),
                     "targets_cis_in_ensemble_but_not_model1":
                         sum(1 for r in rows if r["ens_n_models_with_cis"] > 0
                             and r["n_cis_omega30"] == 0)},
        "pool": {"universe_windows": sum(r["univ_n_windows"] for r in rows),
                 "universe_cis_windows": sum(r["univ_n_cis_windows"] for r in rows),
                 "universe_min_step": min(r["univ_min_step"] for r in rows),
                 "pool500_cis_windows": sum(r["pool_n_cis_windows"] for r in rows),
                 "pool500_min_step": min(r["pool_min_step"] for r in rows),
                 "top75_cis_windows": sum(r["top75_n_cis_windows"] for r in rows),
                 "top75_min_step": min(r["top75_min_step"] for r in rows)},
        "filters_cited": {"native_step_gate": "core/data.py:406-407 step.min()<3.5 or step.max()>4.1",
                          "fragment_step_gate": "core/data.py:697-698 (step>3.5)&(step<4.1)",
                          "rebuild_tol": "core/data.py REBUILD_TOL=1.5 (N/CA/C basis)"},
    }


# ------------------------------------------------------------------ part 2 (gated)
FLOOR_KEYS = ("pdb", "n", "fold", "rebuild_bb", "floor_ca", "rmsd_avg", "rmsd_arm",
              "chain_cost", "max_omega_dev", "n_cis_omega30")


def floor_target(t, census_row) -> dict:
    """ORACLE DIAGNOSTIC. Gated. The projection's representation floor on the native's own
    torsions (ideal TRANS geometry), CA basis and N/CA/C basis, beside the production cost."""
    L.require_gate("ph_cis floor")
    import core
    from core import geometry as geo
    pdb = t["pdb"]
    db = core.backend("data")
    pep = db.by_pdb(pdb)
    nseq, N, CA, C = L.native_backbone(pdb, model_index=0)
    phi, psi = geo.extract_torsions(N, CA, C)
    built = geo.build_backbone(phi, psi)
    floor_ca = L.ca_rmsd(built["CA"], CA)
    rec = L.prod_record_oracle(pdb)
    return {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
            "rebuild_bb": float(pep.rebuild) if pep is not None else float("nan"),
            "floor_ca": float(floor_ca),
            "rmsd_avg": float(rec["rmsd_avg"]), "rmsd_arm": float(rec["rmsd_arm"]),
            "chain_cost": float(rec["rmsd_arm"] - rec["rmsd_avg"]),
            "max_omega_dev": float(census_row["max_omega_dev"]),
            "mean_omega_dev": float(census_row["mean_omega_dev"]),
            "n_cis_omega30": int(census_row["n_cis_omega30"])}


def floor() -> dict:
    L.require_gate("ph_cis floor")
    from s24 import stats_lib as ST
    with open(CENSUS_JSON, encoding="utf-8") as fh:
        cen = {r["pdb"]: r for r in json.load(fh)["rows"]}
    tg = L.targets()
    rows = [floor_target(t, cen[t["pdb"]]) for t in tg]
    pdbs = [r["pdb"] for r in rows]
    folds = L.folds_of(pdbs)
    fc = np.array([r["floor_ca"] for r in rows]); fb = np.array([r["rebuild_bb"] for r in rows])
    cost = np.array([r["chain_cost"] for r in rows])
    dev = np.array([r["max_omega_dev"] for r in rows])
    mdev = np.array([r["mean_omega_dev"] for r in rows])
    cis = np.array([r["n_cis_omega30"] > 0 for r in rows])
    from scipy.stats import spearmanr
    summ = {"floor_ca": L.mean_se(fc), "rebuild_bb": L.mean_se(fb),
            "chain_cost": L.mean_se(cost),
            "n_cis_targets": int(cis.sum()),
            "floor_ca_cis": L.mean_se(fc[cis]) if cis.any() else None,
            "floor_ca_noncis": L.mean_se(fc[~cis]),
            "chain_cost_cis": L.mean_se(cost[cis]) if cis.any() else None,
            "chain_cost_noncis": L.mean_se(cost[~cis]),
            "rho_floor_ca_vs_max_omega_dev": float(spearmanr(fc, dev)[0]),
            "rho_floor_ca_vs_mean_omega_dev": float(spearmanr(fc, mdev)[0]),
            "rho_chain_cost_vs_max_omega_dev": float(spearmanr(cost, dev)[0]),
            "rho_chain_cost_vs_floor_ca": float(spearmanr(cost, fc)[0]),
            "share_of_chain_cost_explained_if_additive":
                float(fc.mean() / cost.mean()) if cost.mean() != 0 else float("nan")}
    # the floor is an ORACLE lower bound; compare it to the production cost paired
    cmp_ = ST.compare(cost, fc, folds, names=pdbs,
                      label="production chain cost (rmsd_arm - rmsd_avg) MINUS the ORACLE "
                            "CA floor of ideal trans geometry on the native's own torsions")
    print(ST.fmt(cmp_))
    summ["compare_cost_vs_floor"] = {k: v for k, v in cmp_.items() if k != "concentration"}
    L.save(FLOOR_JSON, {"what": "ORACLE DIAGNOSTIC: projection representation floor on cis vs "
                                "non-cis dev targets; basis stated per key (floor_ca: CA; "
                                "rebuild_bb: N/CA/C; rmsd_*: production basis)",
                        "rows": rows, "summary": summ},
           rows=rows, complete_keys=FLOOR_KEYS, n_expected=126, module_file=__file__)
    print(json.dumps({k: v for k, v in summ.items() if k != "compare_cost_vs_floor"}, indent=1))
    return summ


# ------------------------------------------------------------------ part 2b (gated, addendum 2)
FLOOR2_JSON = os.path.join(L.RESULTS, "ph_cis_floor2.json")
FLOOR2_KEYS = ("pdb", "n", "fold", "floor2_lam03", "floor2_lam0", "floor_ca", "chain_cost")


def floor2() -> dict:
    """ORACLE DIAGNOSTIC (PREREG_cis addendum 2): the native CA trace itself projected through
    the PRODUCTION projection (`I.project`, ramah at 0.3, multi-start, exact gradient) and its
    lam=0 rung, each scored against the native. The tight representation floor of the manifold."""
    L.require_gate("ph_cis floor2")
    from s12 import instrument as I
    from s24 import stats_lib as ST
    with open(FLOOR_JSON, encoding="utf-8") as fh:
        prev = {r["pdb"]: r for r in json.load(fh)["rows"]}
    tg = L.targets()
    done = L.read_cells("ph_cis_floor2")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        u = L.univ_oracle(t["pdb"])
        nat = u["nat_ca"]
        pr = I.project(nat, t["seq"], int(t["fold"]))
        row = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
               "floor2_lam03": L.ca_rmsd(pr["ca"], nat), "floor2_lam0": L.ca_rmsd(pr["fit_ca"], nat),
               "floor_ca": float(prev[t["pdb"]]["floor_ca"]), "chain_cost": float(prev[t["pdb"]]["chain_cost"]),
               "max_omega_dev": float(prev[t["pdb"]]["max_omega_dev"])}
        L.write_cell("ph_cis_floor2", t["pdb"], row)
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} floor2 lam0.3 {row['floor2_lam03']:.3f} lam0 {row['floor2_lam0']:.3f} "
              f"rebuild {row['floor_ca']:.3f} ({clk():.0f}s)", flush=True)
    rows = sorted(L.read_cells("ph_cis_floor2").values(), key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]; folds = L.folds_of(pdbs)
    f03 = np.array([r["floor2_lam03"] for r in rows]); f0 = np.array([r["floor2_lam0"] for r in rows])
    fc = np.array([r["floor_ca"] for r in rows]); cost = np.array([r["chain_cost"] for r in rows])
    from scipy.stats import spearmanr
    summ = {"floor2_lam03": L.mean_se(f03), "floor2_lam0": L.mean_se(f0), "floor_ca": L.mean_se(fc),
            "n_floor2_above_rebuild": int((f03 > fc).sum()),
            "rho_floor2_vs_max_omega_dev": float(spearmanr(f03, [r["max_omega_dev"] for r in rows])[0]),
            "rho_floor2_vs_chain_cost": float(spearmanr(f03, cost)[0])}
    for lab, a, b in (("floor2 (native projected, lam 0.3) MINUS floor_ca (own-torsion rebuild)", f03, fc),
                      ("floor2 lam 0.3 MINUS floor2 lam 0", f03, f0),
                      ("production chain cost MINUS floor2 (lam 0.3)", cost, f03)):
        c = ST.compare(a, b, folds, names=pdbs, label=lab)
        print(ST.fmt(c)); summ[lab] = c
    L.save(FLOOR2_JSON, {"what": "ORACLE DIAGNOSTIC: the native projected through the production "
                                 "projection; the tight representation floor (CA basis)",
                         "rows": rows, "summary": summ},
           rows=rows, complete_keys=FLOOR2_KEYS, n_expected=126, module_file=__file__)
    return summ


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("census", "floor", "floor2"))
    a = ap.parse_args()
    if a.mode == "census":
        census()
    elif a.mode == "floor":
        floor()
    else:
        floor2()
