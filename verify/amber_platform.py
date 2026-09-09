"""Does the OpenCL platform move the AMBER/GBn2 block, and at what cost in fidelity?

Four arms on 1A13's native under the production recipe (k=10, steps=0, tolerance=1.0):

  cpu           the golden path: CPU, Threads=1, DeterministicForces
  ocl_double    OpenCL end to end, Precision=double
  ocl_single    OpenCL end to end, default precision -- the trap, since core/amber.py
                sets no Precision property and OpenMM's OpenCL default is single
  hybrid        Hamiltonian BUILT on CPU (so the hydrogen-frame calibration is the
                golden one) and only the minimisation Context on OpenCL/double, which
                isolates "the minimiser ran elsewhere" from "the system was built
                elsewhere"

Reports, per arm: the interaction energy against the pinned -489.9138948277905, the
number of force evaluations the minimiser needed, the CA-RMSD of the emitted structure
against the CPU arm's, and wall time.  Nothing here edits core/amber.py; `refine_coords`
already takes `platform_name` and `builder_for` already keys its cache on it.
"""
import os
import sys
import time

sys.path.insert(0, r"C:\Users\abena\Protein-Folding-Algorithm")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")

import numpy as np
import openmm
import openmm.unit as unit

import protein_geometry as geo
import peptide_db as db
import torsion_lib2 as tl2
import core.amber as A

REF = -489.9138948277905
PID = "1A13"


class _Counter(openmm.MinimizationReporter):
    """Counts what the minimiser actually asked for.  `report` is called once per
    energy+force evaluation, which is the unit `core/amber.py`'s header counts in."""

    def __init__(self):
        super().__init__()
        self.n = 0

    def report(self, iteration, x, grad, args):
        self.n += 1
        return False


def _patched_minimize(counter):
    orig = openmm.LocalEnergyMinimizer.minimize

    def m(ctx, tol=10.0, maxIter=0, reporter=None):
        return orig(ctx, tol, maxIter, counter)
    return orig, m


def _swap_context(H, props):
    """Rebuild the builder's Context on OpenCL, keeping the CPU-built System."""
    plat = openmm.Platform.getPlatformByName("OpenCL")
    H.platform = plat
    H.platform_properties = props
    H.context = openmm.Context(
        H.system, openmm.VerletIntegrator(0.001 * unit.picoseconds), plat, props)
    return H


def run(tag, platform_name, props=None, hybrid=False):
    p = db.by_pdb(PID)
    tab = tl2.library_for(p.seq, 8, p.seq)
    rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
    rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))

    A.clear_cache()
    t0 = time.time()
    if hybrid:
        H = A.builder_for(p.seq, rep, "CPU", 1)
        _swap_context(H, props or {})
    else:
        H = A.builder_for(p.seq, rep, platform_name, 1)
        if props:
            _swap_context(H, props) if platform_name == "OpenCL" else None
    t_build = time.time() - t0

    counter = _Counter()
    orig, patched = _patched_minimize(counter)
    openmm.LocalEnergyMinimizer.minimize = staticmethod(patched)
    try:
        t0 = time.time()
        out = A._run(H, H._heavy_positions(rb, chi1=None),
                     A.K_MODERATE, 0, 1.0, True)
        wall = time.time() - t0
    finally:
        openmm.LocalEnergyMinimizer.minimize = staticmethod(orig)

    cm = out["components"]
    inter = cm["nonbonded"] + cm["solvation"]
    plat_actual = H.context.getPlatform().getName()
    prec = ""
    try:
        prec = H.context.getPlatform().getPropertyValue(H.context, "Precision")
    except Exception:
        prec = "n/a"
    return {"tag": tag, "platform": plat_actual, "precision": prec,
            "interaction": float(inter), "d_ref": float(inter - REF),
            "evals": counter.n, "wall": wall, "build": t_build,
            "energy": float(out["energy"]),
            "ca": np.asarray(out["ca"], float),
            "restraint_rmsd": float(out["restraint_rmsd"])}


def ca_rmsd(a, b):
    import core
    return float(core.backend("numerics").kabsch_rmsd_batch(a[None], b)[0])

if __name__ == "__main__":
    arms = [("cpu", "CPU", None, False),
            ("ocl_single", "OpenCL", None, False),
            ("ocl_double", "OpenCL", {"Precision": "double"}, False),
            ("hybrid_double", "OpenCL", {"Precision": "double"}, True)]
    res = []
    for tag, plat, props, hyb in arms:
        try:
            r = run(tag, plat, props, hyb)
            res.append(r)
            print(f"{r['tag']:14s} {r['platform']:8s} prec={r['precision']:7s} "
                  f"inter={r['interaction']!r}", flush=True)
            print(f"               d(ref)={r['d_ref']:+.6e}  evals={r['evals']:5d}  "
                  f"build={r['build']:6.2f}s  minimise={r['wall']:7.3f}s", flush=True)
        except Exception as e:                                       # noqa: BLE001
            print(f"{tag:14s} FAILED {type(e).__name__}: {str(e)[:400]}", flush=True)
    if res:
        base = res[0]["ca"]
        print("\nCA-RMSD against the CPU arm's emitted structure:")
        for r in res:
            print(f"  {r['tag']:14s} {ca_rmsd(r['ca'], base):.6e} A   "
                  f"restraint_rmsd {r['restraint_rmsd']:.4f} A")
        print("\nspeed against the CPU arm (minimisation only):")
        for r in res:
            print(f"  {r['tag']:14s} {res[0]['wall']/max(r['wall'],1e-9):6.2f}x")
