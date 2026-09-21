#!/usr/bin/env python
"""s31/s31_B1_achirality.py -- lane B, falsifiers F1/F2/F3 of PREREG_S31_B.md.

THE GATE ON B1's COMPUTE. Lemma B1 claims the shipped potential
`app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")` is reflection-invariant,
U(Rx) = U(x) for R the point reflection. If it is, then by G1 the free energy F = E - TS, the
entropy S, and every temperature derivative of F are functions of the candidate's pairwise
distance map and lie inside the class S30 closed -- and no thermodynamics needs to be computed.

This script does not assume the lemma; it measures it.

  F2  ff14SB's PeriodicTorsionForce phase table: every gamma must be 0 or pi (tol 1e-9 rad).
      A CMAPTorsionForce, if present, is reported -- it is NOT even in general.
  F1  |U(x) - U(Rx)| / (|U(x)| + 1) per force group and in total, on real target structures.
      Reflection here is the EXACT point reflection of every atom about the centroid.
  F3  the Legacy potential term by term: T_even/T_odd variance share under (phi,psi) -> (-phi,-psi).
      F3 is the pre-check for B2 -- if LEG_torsion is even, B2 closes immediately by G1.

USAGE
    python s31/s31_B1_achirality.py            # all three
    python s31/s31_B1_achirality.py --n 8      # how many PDBs for F1

ORACLE CONTENT: none. Nothing here reads a native RMSD. F1/F2/F3 are properties of the potential.
"""
from __future__ import annotations

import argparse
import glob
import io
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s31_B1_achirality.json")

FF_XMLS = ("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")   # verbatim from core/amber.py:927


# ---------------------------------------------------------------- F2: the phase table
def f2_phase_table(seq_pdb):
    """Every PeriodicTorsionForce phase in a system built from the shipped forcefield."""
    import openmm as mm
    from openmm import app, unit

    pdb = app.PDBFile(seq_pdb)
    ff = app.ForceField(*FF_XMLS)
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.addHydrogens(ff)
    sysm = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff,
                           constraints=None, rigidWater=False)
    phases, kinds = [], {}
    for f in sysm.getForces():
        kinds[f.__class__.__name__] = kinds.get(f.__class__.__name__, 0) + 1
        if isinstance(f, mm.PeriodicTorsionForce):
            for i in range(f.getNumTorsions()):
                *_, per, ph, k = f.getTorsionParameters(i)
                phases.append(float(ph.value_in_unit(unit.radian)))
    ph = np.asarray(phases, float)
    # distance to the nearest of {0, pi} modulo 2pi
    w = np.mod(ph, 2 * math.pi)
    d = np.minimum(np.minimum(np.abs(w - 0.0), np.abs(w - 2 * math.pi)), np.abs(w - math.pi))
    return dict(n_torsions=int(len(ph)), forces=kinds,
                max_dist_to_0_or_pi_rad=float(d.max()) if len(d) else 0.0,
                unique_phases_rad=sorted({round(float(x), 12) for x in w}),
                has_cmap=bool(kinds.get("CMAPTorsionForce", 0)),
                F2_fired=bool(len(d) and d.max() > 1e-9))


# ---------------------------------------------------------------- F1: U(x) vs U(Rx)
def f1_reflection(pdb_paths):
    import openmm as mm
    from openmm import app, unit

    ff = app.ForceField(*FF_XMLS)
    rows = []
    for p in pdb_paths:
        try:
            pdb = app.PDBFile(p)
            mod = app.Modeller(pdb.topology, pdb.positions)
            mod.addHydrogens(ff)
            sysm = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff,
                                   constraints=None, rigidWater=False)
            # one force group per force so the comparison is term by term
            names = []
            for g, f in enumerate(sysm.getForces()):
                f.setForceGroup(g)
                names.append(f.__class__.__name__)
            integ = mm.VerletIntegrator(1.0 * unit.femtosecond)
            plat = mm.Platform.getPlatformByName("CPU")
            ctx = mm.Context(sysm, integ, plat, {"Threads": "1"})

            X = np.asarray(mod.positions.value_in_unit(unit.nanometer), float)
            Xm = 2.0 * X.mean(0)[None, :] - X            # EXACT point reflection about centroid
            # MATCHED CONTROL (contract rule 7): a pure ROTATION is also ANALYTICALLY exact
            # for this potential, so its numerical residual IS the floating-point noise floor
            # for a coordinate transform of this magnitude. The reflection residual must be
            # read against it, not against zero.
            rng = np.random.default_rng(1234)
            A_ = rng.normal(size=(3, 3))
            Q_, R_ = np.linalg.qr(A_)
            Q_ = Q_ * np.sign(np.diag(R_))[None, :]
            if np.linalg.det(Q_) < 0:
                Q_[:, 0] *= -1.0                          # proper rotation, det = +1
            Xr = (X - X.mean(0)[None, :]) @ Q_.T + X.mean(0)[None, :]

            def energies(P):
                ctx.setPositions(P * unit.nanometer)
                tot = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                    unit.kilocalorie_per_mole)
                per = {}
                for g, nm in enumerate(names):
                    e = ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy()
                    per[nm] = float(e.value_in_unit(unit.kilocalorie_per_mole))
                return float(tot), per

            e0, p0 = energies(X)
            e1, p1 = energies(Xm)
            e2, p2 = energies(Xr)
            rel = abs(e1 - e0) / (abs(e0) + 1.0)
            rel_rot = abs(e2 - e0) / (abs(e0) + 1.0)
            rows.append(dict(pdb=os.path.basename(p)[:-4], n_atoms=int(sysm.getNumParticles()),
                             E=e0, E_mirror=e1, E_rot=e2, rel=rel, rel_rot=rel_rot,
                             per_group={k: dict(E=p0[k], E_mirror=p1[k], E_rot=p2[k],
                                                rel=abs(p1[k] - p0[k]) / (abs(p0[k]) + 1.0),
                                                rel_rot=abs(p2[k] - p0[k]) / (abs(p0[k]) + 1.0))
                                        for k in p0}))
            del ctx, integ
        except Exception as exc:                                   # noqa: BLE001
            rows.append(dict(pdb=os.path.basename(p)[:-4], error=repr(exc)[:300]))
    ok = [r for r in rows if "rel" in r]
    worst = max((r["rel"] for r in ok), default=float("nan"))
    worst_rot = max((r["rel_rot"] for r in ok), default=float("nan"))
    grp, grp_rot = {}, {}
    for r in ok:
        for k, v in r["per_group"].items():
            grp[k] = max(grp.get(k, 0.0), v["rel"])
            grp_rot[k] = max(grp_rot.get(k, 0.0), v["rel_rot"])
    return dict(n_ok=len(ok), n_fail=len(rows) - len(ok), max_rel=float(worst),
                max_rel_rotation_control=float(worst_rot),
                reflection_over_rotation=float(worst / worst_rot) if worst_rot > 0 else float("inf"),
                max_rel_per_group=grp, max_rel_per_group_rotation=grp_rot, rows=rows,
                F1_fired=bool(len(ok) and worst > 1e-6))


# ---------------------------------------------------------------- F3: the Legacy terms
def f3_legacy_parity(n_targets=20, n_members=200, seed=0):
    """Even/odd split of every Legacy component under (phi,psi) -> (-phi,-psi).

    Measured on REAL pool members (the shipped K=500 BLOSUM pool), not on random torsions,
    because the parity share is a property of the term ON THE ENSEMBLE IT SCORES.
    """
    from s12 import instrument as I
    from s16 import energy_lib as EL

    tg = I.targets()
    rng = np.random.default_rng(seed)
    pick = sorted(rng.choice(len(tg), size=min(n_targets, len(tg)), replace=False).tolist())
    acc = {t: {"even": [], "odd": [], "tot": []} for t in EL.LEG_TERMS}
    per_target = []
    for ti in pick:
        t = tg[ti]
        u = I.load_univ(t["pdb"])
        idx = I.pool_idx(u)[:n_members]
        PHI, PSI = u["PHI"][idx], u["PSI"][idx]
        cp = EL.legacy_components_of_windows(t["seq"], PHI, PSI)
        cm = EL.legacy_components_of_windows(t["seq"], -PHI, -PSI)
        row = {"pdb": t["pdb"], "n": t["n"], "fold": t["fold"]}
        for term in EL.LEG_TERMS:
            a, b = np.asarray(cp[term], float), np.asarray(cm[term], float)
            ev, od = 0.5 * (a + b), 0.5 * (a - b)
            vt = float(np.var(a))
            row[term] = dict(var_total=vt, var_even=float(np.var(ev)), var_odd=float(np.var(od)),
                             odd_share=float(np.var(od) / vt) if vt > 1e-18 else float("nan"),
                             mean=float(a.mean()), mean_mirror=float(b.mean()))
            acc[term]["even"].append(np.var(ev)); acc[term]["odd"].append(np.var(od))
            acc[term]["tot"].append(vt)
        per_target.append(row)
    summ = {}
    for term in EL.LEG_TERMS:
        tot = np.asarray(acc[term]["tot"], float)
        od = np.asarray(acc[term]["odd"], float)
        m = tot > 1e-18
        summ[term] = dict(
            mean_odd_share=float(np.mean(od[m] / tot[m])) if m.any() else float("nan"),
            median_odd_share=float(np.median(od[m] / tot[m])) if m.any() else float("nan"),
            n_targets_nondegenerate=int(m.sum()))
    tors = summ.get("torsion", {}).get("mean_odd_share", float("nan"))
    return dict(n_targets=len(pick), n_members=n_members, seed=seed, per_term=summ,
                per_target=per_target,
                F3_fired=bool(np.isfinite(tors) and tors < 0.10))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="PDBs for F1")
    ap.add_argument("--f3-targets", type=int, default=20)
    ap.add_argument("--f3-members", type=int, default=200)
    a = ap.parse_args()

    pdbs = sorted(glob.glob(os.path.join(ROOT, "pdbs_ext", "*.pdb")))[: a.n]
    if not pdbs:
        pdbs = sorted(glob.glob(os.path.join(ROOT, "pdbs", "*.pdb")))[: a.n]

    out = {"forcefield": list(FF_XMLS), "prereg": "s31/PREREG_S31_B.md"}

    print("F2  ff14SB torsion phase table ...", flush=True)
    out["F2"] = f2_phase_table(pdbs[0])
    print("    n_torsions=%(n_torsions)d  max dist to {0,pi} = %(max_dist_to_0_or_pi_rad).3e  "
          "CMAP=%(has_cmap)s  FIRED=%(F2_fired)s" % out["F2"], flush=True)
    print("    forces:", out["F2"]["forces"], flush=True)

    print("F1  U(x) vs U(Rx) on %d structures ..." % len(pdbs), flush=True)
    out["F1"] = f1_reflection(pdbs)
    print("    n_ok=%(n_ok)d n_fail=%(n_fail)d  max relative |dE| = %(max_rel).3e  "
          "FIRED=%(F1_fired)s" % out["F1"], flush=True)
    print("    ROTATION CONTROL (analytically exact -> the FP noise floor): %.3e   "
          "reflection/rotation = %.2fx"
          % (out["F1"]["max_rel_rotation_control"], out["F1"]["reflection_over_rotation"]),
          flush=True)
    print("      %-28s %10s %10s" % ("force group", "reflect", "rotate"), flush=True)
    for k, v in sorted(out["F1"]["max_rel_per_group"].items(), key=lambda kv: -kv[1]):
        print("      %-28s %.3e %.3e" % (k, v, out["F1"]["max_rel_per_group_rotation"][k]),
              flush=True)

    print("F3  Legacy term parity on real pool members ...", flush=True)
    out["F3"] = f3_legacy_parity(a.f3_targets, a.f3_members)
    print("    %-22s %10s %10s" % ("legacy term", "odd share", "median"), flush=True)
    for k, v in sorted(out["F3"]["per_term"].items(), key=lambda kv: -(kv[1]["mean_odd_share"]
                                                                      if np.isfinite(kv[1]["mean_odd_share"]) else -9)):
        print("    %-22s %10.4f %10.4f" % (k, v["mean_odd_share"], v["median_odd_share"]), flush=True)
    print("    F3 FIRED (torsion essentially achiral) = %s" % out["F3"]["F3_fired"], flush=True)

    with io.open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)
    print("\nwrote", OUT, flush=True)


if __name__ == "__main__":
    main()
