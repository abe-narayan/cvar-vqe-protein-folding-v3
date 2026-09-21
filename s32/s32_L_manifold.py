"""LANE L / L-H1 -- the registered test: the DEPLOYED projector, run on the native itself.

PREREG `s32/PREREG_S32_L.md` at commit 88f2da39, hypothesis L-H1.

The L1 census measured the representability floor through the NATIVE-TORSION REBUILD,
which is only an UPPER bound on the distance from a native to the set of emittable
structures: a closer ideal-geometry chain might exist that simply does not use the native's
own torsions.  The registered falsifier is exactly that possibility.

So this module runs the operator the pipeline actually uses -- `core.project.lam_path` at
lambda = 0, multi-start, the identical call stage 3b makes -- on the NATIVE CA TRACE.  Its
residual is the tightest available bound on the floor, produced by the deployed code.

Three arms, all ORACLE / NOT DEPLOYABLE (they read native coordinates as INPUT):

    PROJ_NAT_SHORT   deployed projector on the 126 canonical natives      (L = 9-16)
    PROJ_NAT_LONG    deployed projector on 40-60 residue natives          (L = 40-60)
    AMBER_ESCAPE     stage 4 applied to the projected native, both lengths

AMBER_ESCAPE is the registered second control.  Stage 4 is the only stage that can leave
the ideal-geometry manifold, so if it recovers the floor, the floor does not bind at the
endpoint.  The pre-registration fixes the threshold in advance: **recovery of more than
50% of the floor at L = 45 falsifies the claim that the floor binds.**

    python -m s32.s32_L_manifold short
    python -m s32.s32_L_manifold long
    python -m s32.s32_L_manifold amber
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)

LONG_LO, LONG_HI = 40, 60
N_LONG = int(os.environ.get("N_LONG", "60"))


def _rmsd(a, b):
    from core import geometry as geo
    return float(geo.rmsd(geo.kabsch_superpose(np.asarray(a, float),
                                               np.asarray(b, float)),
                          np.asarray(b, float)))


def _project_native(ca, seq, fold, maxiter=300):
    """The deployed lambda=0 arm on a native CA trace.  Returns (ca_fit, phi, psi, secs)."""
    from core import project as pj
    pen = pj.make_penalty("ramah", seq, int(fold))
    t0 = time.time()
    path = pj.lam_path(np.asarray(ca, float), pen, (0.0, 0.3), maxiter=maxiter,
                       multi=True, grad="exact")
    dt = time.time() - t0
    fit = path[0.0]
    arm = path[0.3]
    return (np.asarray(fit[0], float), np.asarray(fit[1], float),
            np.asarray(fit[2], float), np.asarray(arm[0], float), dt)


def _rebuild_floor(ca, phi, psi):
    from core import geometry as geo
    built = geo.build_backbone(np.nan_to_num(np.asarray(phi, float)),
                               np.nan_to_num(np.asarray(psi, float)))
    return float(geo.rmsd(geo.kabsch_superpose(built["CA"], np.asarray(ca, float)),
                          np.asarray(ca, float)))


# --------------------------------------------------------------------------- short arm
def stage_short(verbose=True):
    """PROJ_NAT_SHORT -- deployed projector on the 126 canonical natives.

    Reads `s8/generate_univ` for coordinates ONLY.  Writes nothing there, re-splits
    nothing, and produces no number about the 3.2105 endpoint.
    """
    sys.path.insert(0, BASE)
    from s12 import instrument as ins
    rows = []
    t0 = time.time()
    for k, t in enumerate(ins.targets()):
        u = ins.load_univ(t["pdb"])
        nat = np.asarray(u["nat_ca"], float)
        fit, ph, ps, arm, dt = _project_native(nat, t["seq"], t["fold"])
        rows.append(dict(pdb=t["pdb"], n=t["n"], fold=int(t["fold"]),
                         proj=_rmsd(fit, nat), proj_lam03=_rmsd(arm, nat), secs=dt))
        if verbose and (k + 1) % 25 == 0:
            print(f"  {k+1}/126 ({time.time()-t0:.0f}s)", flush=True)
    return _finish("L2_proj_native_short", rows, verbose)


# ---------------------------------------------------------------------------- long arm
def _long_targets(n=N_LONG):
    """40-60 residue chains from the L1 census, lowest PDB code first (a fixed rule)."""
    cen = json.load(open(os.path.join(RESULTS, "L1_census.json")))
    rows = [r for r in cen["rows"] if LONG_LO <= r["n"] <= LONG_HI]
    rows.sort(key=lambda r: r["pdb"])
    return rows[:n]


def stage_long(verbose=True):
    """PROJ_NAT_LONG -- the identical operator, identical schedule, only L differs."""
    from core import geometry as geo
    import glob
    rows = []
    t0 = time.time()
    tg = _long_targets()
    for k, r in enumerate(tg):
        path = os.path.join(BASE, "prots", r["pdb"] + ".pdb")
        if not os.path.exists(path):
            cand = glob.glob(os.path.join(BASE, "prots", r["pdb"].lower() + ".pdb"))
            if not cand:
                continue
            path = cand[0]
        seq, coords, phi, psi = geo.native_coords_from_pdb(path)
        nat = np.asarray(coords["CA"], float)
        # fold 0 is an arbitrary label here: the lambda=0 arm never touches the penalty,
        # and only the lambda=0 arm is the floor.  The lambda=0.3 column is reported too
        # and is labelled as penalty-dependent.
        fit, ph, ps, arm, dt = _project_native(nat, seq, 0)
        rows.append(dict(pdb=r["pdb"], n=int(r["n"]), fold=0,
                         proj=_rmsd(fit, nat), proj_lam03=_rmsd(arm, nat),
                         rebuild=_rebuild_floor(nat, phi, psi), secs=dt))
        if verbose:
            print(f"  [{k+1}/{len(tg)}] {r['pdb']} n={r['n']:2d} "
                  f"proj={rows[-1]['proj']:.3f} rebuild={rows[-1]['rebuild']:.3f} "
                  f"({dt:.1f}s)", flush=True)
    return _finish("L2_proj_native_long", rows, verbose)


# --------------------------------------------------------------------------- the escape
def stage_amber(verbose=True):
    """AMBER_ESCAPE -- can stage 4 recover the floor?  Registered threshold: 50%."""
    from core import geometry as geo
    import glob
    sys.path.insert(0, BASE)
    import amber_refine as ar
    import torsion_lib2 as tl2

    def relax(seq, phi, psi):
        BB = geo.build_backbone_batch(np.asarray(phi, float)[None],
                                      np.asarray(psi, float)[None])
        cd = {k: v[0] for k, v in BB.items()}
        tab = tl2.library_for(seq, 9, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        r = ar.refine_coords(seq, rep, cd, k_restraint=10.0, steps=0, components=True)
        ar.clear_cache()
        return np.asarray(r["ca"], float)

    out = {}
    for which, src in (("short", "L2_proj_native_short"),
                       ("long", "L2_proj_native_long")):
        prev = json.load(open(os.path.join(RESULTS, src + ".json")))
        rows = []
        if which == "short":
            from s12 import instrument as ins
            tgt = {t["pdb"]: t for t in ins.targets()}
        for k, r in enumerate(prev["rows"]):
            try:
                if which == "short":
                    u = ins.load_univ(r["pdb"])
                    nat = np.asarray(u["nat_ca"], float); seq = tgt[r["pdb"]]["seq"]
                    fold = int(r["fold"])
                else:
                    p = os.path.join(BASE, "prots", r["pdb"] + ".pdb")
                    if not os.path.exists(p):
                        c = glob.glob(os.path.join(BASE, "prots", r["pdb"].lower() + ".pdb"))
                        if not c:
                            continue
                        p = c[0]
                    seq, coords, _, _ = geo.native_coords_from_pdb(p)
                    nat = np.asarray(coords["CA"], float); fold = 0
                fit, ph, ps, arm, _ = _project_native(nat, seq, fold)
                ca1 = relax(seq, ph, ps)
                rows.append(dict(pdb=r["pdb"], n=int(r["n"]), fold=fold,
                                 proj=_rmsd(fit, nat), after_amber=_rmsd(ca1, nat)))
            except Exception as e:                                       # noqa: BLE001
                rows.append(dict(pdb=r["pdb"], n=int(r["n"]), err=str(e)[:160]))
            if verbose and (k + 1) % 10 == 0:
                print(f"  {which} {k+1}/{len(prev['rows'])}", flush=True)
        good = [r for r in rows if "after_amber" in r]
        pj = np.array([r["proj"] for r in good]); am = np.array([r["after_amber"] for r in good])
        out[which] = dict(n=len(good), proj_mean=float(pj.mean()),
                          after_amber_mean=float(am.mean()),
                          recovered_mean=float((pj - am).mean()),
                          recovered_frac=float((pj - am).mean() / pj.mean()) if pj.mean() else None,
                          n_err=len(rows) - len(good), rows=rows)
    path = os.path.join(RESULTS, "L2_amber_escape.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: {a: b for a, b in v.items() if a != "rows"}
                          for k, v in out.items()}, indent=1))
        print("wrote", path, flush=True)
    return out


def _finish(name, rows, verbose):
    pr = np.array([r["proj"] for r in rows])
    out = dict(n=len(rows), arm=name, oracle=True, deployable=False,
               proj=dict(mean=float(pr.mean()), median=float(np.median(pr)),
                         se=float(pr.std(ddof=1) / max(1, len(pr)) ** 0.5),
                         p90=float(np.percentile(pr, 90)), max=float(pr.max())),
               secs_mean=float(np.mean([r["secs"] for r in rows])), rows=rows)
    if rows and "rebuild" in rows[0]:
        rb = np.array([r["rebuild"] for r in rows])
        out["rebuild"] = dict(mean=float(rb.mean()), median=float(np.median(rb)))
    path = os.path.join(RESULTS, name + ".json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "short"
    {"short": stage_short, "long": stage_long, "amber": stage_amber}[stage]()
