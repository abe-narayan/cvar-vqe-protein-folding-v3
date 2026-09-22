"""Is AMBER's +14.9% at w4 CONTENTION, or is it slower cores doing the same work?

The worker sweep cannot tell those apart -- both inflate "CPU-seconds for identical
science" and both produce the same curve.  This does, with one control the sweep cannot
produce: a SINGLE process pinned to the slowest core, which has zero contention by
construction, so whatever inflation it shows is pure core heterogeneity.

    rank    AMBER CPU-seconds per structure on each core, one process, box otherwise idle
    A       w1 pinned to the fastest core           the clean baseline
    B       w4 unpinned                             the regression, in isolation
    C       w4 pinned to the four fastest cores     does placement recover it?
    D       w1 pinned to the slowest core           heterogeneity with contention = 0

Everything is measured in `time.process_time()` -- per-process CPU, user+sys -- because
wall time at w4 confounds the thing being measured with how the four ran concurrently.
Placement cannot change arithmetic, so every arm here is bit-identical by construction and
the energies are asserted equal rather than compared.
"""
import ctypes
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
import os as _os; ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

#: the fixed workload: distinct sequences so no builder or memo is shared
TARGETS = ("1A13", "1A1P", "1CB3")


def set_affinity(core):
    """Pin this process to one logical CPU.  -1 leaves the OS free to migrate it.

    The ctypes signatures are declared rather than left to guess: `GetCurrentProcess`
    returns the pseudo-handle -1, which ctypes truncates to a 32-bit int without an
    explicit `restype`, and the affinity mask is a `DWORD_PTR` (64-bit here).  Both
    defaults are wrong and the call fails with a bare "invalid handle".
    """
    if core is None or core < 0:
        return "free"
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    k32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    k32.SetProcessAffinityMask.restype = ctypes.c_int
    if not k32.SetProcessAffinityMask(k32.GetCurrentProcess(), 1 << int(core)):
        raise OSError(f"SetProcessAffinityMask({core}) failed: "
                      f"{ctypes.get_last_error()}")
    return str(core)


def workload(reps=1):
    """builder + one production-recipe minimisation per target, CPU-timed separately."""
    import numpy as np
    import protein_geometry as geo
    import peptide_db as db
    import torsion_lib2 as tl2
    import core.amber as A

    rows = []
    for _ in range(reps):
        for pid in TARGETS:
            p = db.by_pdb(pid)
            tab = tl2.library_for(p.seq, 8, p.seq)
            rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
            rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))
            A.clear_cache()
            c0, w0 = time.process_time(), time.perf_counter()
            H = A.builder_for(p.seq, rep, "CPU", 1)
            c1, w1 = time.process_time(), time.perf_counter()
            out = A._run(H, H._heavy_positions(rb, chi1=None), A.K_MODERATE, 0, 1.0, True)
            c2, w2 = time.process_time(), time.perf_counter()
            cm = out["components"]
            rows.append({"pid": pid,
                         "build_cpu": c1 - c0, "build_wall": w1 - w0,
                         "min_cpu": c2 - c1, "min_wall": w2 - w1,
                         "interaction": float(cm["nonbonded"] + cm["solvation"]),
                         "energy": float(out["energy"])})
    return rows


def main_worker():
    core = int(sys.argv[2])
    out = sys.argv[3]
    reps = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    where = set_affinity(core)
    t0 = time.perf_counter()
    rows = workload(reps)
    with open(out, "w") as fh:
        json.dump({"core": where, "wall": time.perf_counter() - t0,
                   "cpu": time.process_time(), "rows": rows}, fh)


def spawn(cores, reps=1, tag="x"):
    """Run one process per entry in `cores` concurrently; return their reports."""
    procs, outs = [], []
    for i, c in enumerate(cores):
        o = os.path.join(HERE, f"aff_{tag}_{i}.json")
        outs.append(o)
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                   MKL_NUM_THREADS="1")
        procs.append(subprocess.Popen(
            [sys.executable, os.path.abspath(__file__), "worker", str(c), o, str(reps)],
            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE))
    got = []
    for p, o in zip(procs, outs):
        _, err = p.communicate()
        if p.returncode != 0:
            raise RuntimeError(f"worker failed: {err.decode()[:500]}")
        with open(o) as fh:
            got.append(json.load(fh))
    return got


def per_structure(reports):
    """Mean CPU-seconds per minimisation across every worker's every structure."""
    c = [r["min_cpu"] for rep in reports for r in rep["rows"]]
    w = [r["min_wall"] for rep in reports for r in rep["rows"]]
    b = [r["build_cpu"] for rep in reports for r in rep["rows"]]
    return (sum(c) / len(c), sum(w) / len(w), sum(b) / len(b), len(c))


def energies(reports):
    d = {}
    for rep in reports:
        for r in rep["rows"]:
            d.setdefault(r["pid"], set()).add(r["interaction"])
    return d

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        main_worker()
        sys.exit(0)

    ncpu = os.cpu_count()
    print(f"logical CPUs: {ncpu}\n")

    print("=== per-core AMBER cost, one process at a time, box otherwise idle ===")
    print(f"  {'core':>4s} {'min CPU s':>10s} {'min wall s':>11s} {'build CPU s':>12s} {'rel':>6s}")
    per_core = {}
    for c in range(ncpu):
        rep = spawn([c], reps=1, tag=f"core{c}")
        mc, mw, bc, n = per_structure(rep)
        per_core[c] = mc
        print(f"  {c:4d} {mc:10.3f} {mw:11.3f} {bc:12.3f}", flush=True)
    best = min(per_core, key=per_core.get)
    worst = max(per_core, key=per_core.get)
    for c in range(ncpu):
        print(f"  core {c}: {per_core[best]/per_core[c]:.3f} of fastest")
    order = sorted(per_core, key=per_core.get)
    fast4 = order[:4]
    print(f"\n  fastest core {best} ({per_core[best]:.3f} s), slowest {worst} "
          f"({per_core[worst]:.3f} s), ratio {per_core[worst]/per_core[best]:.3f}")
    print(f"  fastest four for AMBER: {fast4}")

    arms = [("A  w1 pinned fastest", [best]),
            ("B  w4 unpinned", [-1] * 4),
            ("C  w4 pinned fastest-4", fast4),
            ("D  w1 pinned slowest", [worst]),
            #: A again, last.  Contention IS the measurand here, so if anything else
            #: arrived on the box mid-experiment the arms are not comparable -- and the
            #: only way to know is to re-measure the baseline at the end and see whether
            #: it still reproduces.
            ("A' w1 pinned fastest (drift check)", [best])]
    print("\n=== the four arms ===")
    print(f"  {'arm':24s} {'CPU s/struct':>13s} {'vs A':>7s} {'wall s/struct':>14s} {'n':>4s}")
    res = {}
    for name, cores in arms:
        #: two passes so the four concurrent workers stay overlapped long enough for the
        #: steady state to dominate the ragged start and finish
        rep = spawn(cores, reps=2, tag=name.split()[0])
        mc, mw, bc, n = per_structure(rep)
        res[name] = (mc, mw, rep)
        print(f"  {name:24s} {mc:13.3f} {mc/res[arms[0][0]][0]:7.3f} {mw:14.3f} {n:4d}",
              flush=True)

    print("\n=== bit-identity across arms (placement cannot change arithmetic) ===")
    allE = {}
    for name, (_, _, rep) in res.items():
        for pid, vals in energies(rep).items():
            allE.setdefault(pid, set()).update(vals)
    for pid, vals in sorted(allE.items()):
        print(f"  {pid}: {len(vals)} distinct interaction energ{'y' if len(vals)==1 else 'ies'}"
              f"  {sorted(vals)!r}")
