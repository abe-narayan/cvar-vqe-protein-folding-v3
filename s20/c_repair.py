"""s20/c_repair.py -- AMBER-FIRST vs PROJECTION-FIRST.  The repair path, in the deployed basis.

WHY THIS MODULE EXISTS, AND WHOSE QUESTION IT IS.  The coordinator corrected `BRIEF.md` section 3
mid-sprint: the **+0.164 A tax is the PROJECTION's, not AMBER's**, and the deployed AMBER arm is
`core.pipeline.Config.amber_k = 10.0, amber_steps = 0` applied to the **BUILT CHAIN**, where it
costs +0.0207 A -- a quarter of the MDE.  Every "AMBER k = 30" number in Sprints 16 and 19 is a
DIFFERENT restraint constant applied to a DIFFERENT input (the point cloud).  The question that
follows is the coordinator's and it has not been asked in this basis:

    the pipeline has TWO repair operators for the same defect -- the ideal-geometry PROJECTION
    and restrained AMBER -- and it currently applies them in one fixed order.  Is the other
    order better ON BOTH AXES?

THE DEPLOYED PATH, READ OUT OF `core/pipeline.py` RATHER THAN ASSUMED (BRIEF section 10):

    average(Wsub, Ps)          -> C, the raw-window coordinate average          POINT CLOUD
    project(C, seq, fold, cfg) -> phi, psi          (lam = 0.3, multi, 300)     BUILT CHAIN
    relax(seq, phi, psi, cfg)  -> refine_coords(k_restraint = 10.0, steps = 0)  EMISSION

TWO POINT CLOUDS, NOT ONE, AND THEY ARE NOT THE SAME OBJECT.  `I.coordinate_average(W)` averages
the RAW WINDOWS and is what the pipeline projects (3.0483 A).  `s15.phys_repl
.averaged_backbone_from` averages the ideal-geometry REBUILDS and is the only one of the two that
carries N/C/O/CB, so it is the only one AMBER can consume (3.0498 A).  They differ by up to
0.373 A per atom, measured.  Every arm here therefore starts from `averaged_backbone_from`, and
the deployed projection is ALSO re-run on that same cloud so the two paths share an input; the
cached deployed start is carried beside it as a labelled reference, never as the matched control.

VALIDITY IS A VECTOR.  Sprint 19's `caonly_k300` is the standing warning: 0.110 A more accurate,
indistinguishable clash count, **42% more cis peptide bonds**.  No arm is promoted on a scalar.

    python -m s20.c_repair --smoke
    python -m s20.c_repair
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                          # noqa: E402
from s14.avgspace import top75_windows                   # noqa: E402
from s15.phys_repl import averaged_backbone_from, random_rigid   # noqa: E402
from s15 import seed as SD                               # noqa: E402
from s16 import energy_lib as EL                         # noqa: E402
from s18 import phys_lib as PL                           # noqa: E402
from s19 import a_fit                                    # noqa: E402  (the cached deployed start)

ATOMS = ("N", "CA", "C", "O", "CB")
K_DEPLOYED = 10.0        # core.pipeline.Config.amber_k -- THE deployed restraint constant
STEPS = 0                # core.pipeline.Config.amber_steps -- unbounded, the deployed protocol
TOL = 1.0
#: the AMBER-first ladder, at and around the deployed constant.
K_FIRST = (10.0, 30.0)
N_FRAME_NULL = 15


def _amber(seq, rep, coords, k):
    from core import amber as am
    return am.refine_coords(seq, rep, {a: np.asarray(coords[a], float) for a in coords
                                       if a in ATOMS},
                            k_restraint=float(k), steps=STEPS, tolerance=TOL,
                            threads=1, memo=False)


def _rec(bb, seq, nat, ref_ca, extra=None):
    d = {a: np.asarray(bb[a], float) for a in ATOMS if a in bb}
    o = {"rmsd": float(I.ca_rmsd(d["CA"], nat)),
         "ca_disp": float(I.ca_rmsd(d["CA"], np.asarray(ref_ca, float))),
         "valid": EL.panel(d, seq)}
    if extra:
        o.update(extra)
    return o


def run_target(t, do_frame_null=True):
    from core import geometry as geo
    from core import amber as am
    import torsion_lib2 as tl2
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float)
    nat = np.asarray(u["nat_ca"], float)
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)

    #: THE SHARED INPUT -- the all-atom coordinate average, the only point cloud AMBER can eat.
    avg, _C, _dev = averaged_backbone_from(W, np.asarray(PHI, float), np.asarray(PSI, float))
    cloud = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}
    cloud_ca = cloud["CA"]

    arms = {}
    arms["cloud"] = _rec(cloud, seq, nat, cloud_ca)
    arms["cloud"]["basis"] = "POINT CLOUD -- not a buildable structure"

    #: --- PATH P: projection first (the DEPLOYED order), on the SHARED input ---
    pr = I.project(cloud_ca, seq, fold)
    ph, ps = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)
    proj = geo.build_backbone(ph, ps)
    arms["proj"] = _rec(proj, seq, nat, cloud_ca)
    arms["proj"]["basis"] = "BUILT CHAIN"
    for k in (K_DEPLOYED, 30.0):
        r = _amber(seq, rep, proj, k)
        arms[f"proj_amber{int(k)}"] = _rec(r["backbone"], seq, nat, cloud_ca, {
            "energy": float(r["energy"]), "energy_initial": float(r["energy_initial"]),
            "converged": bool(r["converged"]),
            "disp_from_own_input": float(I.ca_rmsd(np.asarray(r["ca"], float),
                                                   np.asarray(proj["CA"], float)))})
        arms[f"proj_amber{int(k)}"]["basis"] = "REPAIRED EMISSION"
    #: the DEPLOYED reference, from the cached raw-window projection.  Labelled, not matched.
    dph, dps, _davg = a_fit.start(pdb, seq, fold)
    dproj = geo.build_backbone(np.asarray(dph, float), np.asarray(dps, float))
    arms["deployed_proj"] = _rec(dproj, seq, nat, cloud_ca)
    r = _amber(seq, rep, dproj, K_DEPLOYED)
    arms["deployed_emission"] = _rec(r["backbone"], seq, nat, cloud_ca, {
        "energy": float(r["energy"]), "converged": bool(r["converged"]),
        "disp_from_own_input": float(I.ca_rmsd(np.asarray(r["ca"], float),
                                               np.asarray(dproj["CA"], float)))})

    #: --- PATH A: AMBER first, on the SAME shared input ---
    for k in K_FIRST:
        r = _amber(seq, rep, cloud, k)
        nm = f"cloud_amber{int(k)}"
        arms[nm] = _rec(r["backbone"], seq, nat, cloud_ca, {
            "energy": float(r["energy"]), "energy_initial": float(r["energy_initial"]),
            "converged": bool(r["converged"]),
            "disp_from_own_input": float(I.ca_rmsd(np.asarray(r["ca"], float), cloud_ca))})
        arms[nm]["basis"] = "REPAIRED EMISSION (AMBER-first)"
        #: --- PATH A+: AMBER first, THEN the deployed projection.  Both operators, other order.
        pr2 = I.project(np.asarray(r["ca"], float), seq, fold)
        p2 = geo.build_backbone(np.asarray(pr2["phi"], float), np.asarray(pr2["psi"], float))
        arms[nm + "_proj"] = _rec(p2, seq, nat, cloud_ca)
        arms[nm + "_proj"]["basis"] = "BUILT CHAIN (AMBER-first then projection)"

    rec = {"pdb": pdb, "n": n, "fold": fold, "arms": arms}

    #: the rotated-frame null -- exactly zero by construction; reported with its MAXIMUM.
    if do_frame_null:
        rng = SD.stable_rng(pdb, "s20C_repair_frame")
        Q, tv = random_rigid(rng)
        cr = {a: (cloud[a] @ Q.T) + tv for a in cloud}
        r = _amber(seq, rep, cr, K_DEPLOYED)
        rec["frame_null"] = {
            "d_rmsd": float(I.ca_rmsd(np.asarray(r["ca"], float), nat)
                            - arms["cloud_amber10"]["rmsd"]),
            "d_energy": float(r["energy"] - arms["cloud_amber10"]["energy"])}
    am.clear_cache()
    return rec


def run(limit=0, out="c_repair.json", verbose=True):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    path = os.path.join(RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"K_DEPLOYED": K_DEPLOYED, "STEPS": STEPS, "TOL": TOL,
           "K_FIRST": list(K_FIRST), "lam": 0.3,
           "input": "s15.phys_repl.averaged_backbone_from (the all-atom coordinate average)"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.55)
        nfn = sum(1 for r in rows if "frame_null" in r)
        rows.append(run_target(t, do_frame_null=(nfn < N_FRAME_NULL)))
        if verbose:
            a = rows[-1]["arms"]
            print(f"  {len(rows)}/{len(tg)} {t['pdb']} cloud={a['cloud']['rmsd']:.3f} "
                  f"proj={a['proj']['rmsd']:.3f} "
                  f"proj+A10={a['proj_amber10']['rmsd']:.3f} "
                  f"A10={a['cloud_amber10']['rmsd']:.3f} "
                  f"A10+proj={a['cloud_amber10_proj']['rmsd']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg))
    _write(out, rows, cfg, len(tg))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)


def _write(out, rows, cfg, n_expected):
    fn = [r["frame_null"]["d_rmsd"] for r in rows if "frame_null" in r]
    fe = [r["frame_null"]["d_energy"] for r in rows if "frame_null" in r]
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected))}
    if fn:
        obj["frame_null"] = {"n": len(fn), "max_abs": float(np.max(np.abs(fn))),
                             "mean_abs": float(np.mean(np.abs(fn))),
                             "max_abs_energy": float(np.max(np.abs(fe)))}
    p = os.path.join(RESULTS, out)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.smoke:
        run(limit=3, out="_SMOKE_c_repair.json")
    else:
        run(limit=a.limit)
