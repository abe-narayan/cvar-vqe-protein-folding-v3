"""What does `Threads>1` actually buy on the GBn2 kernel, and what does it cost?

The feasibility model assumes an AMBER thread speedup of 1.5 at t2 and 2.2 at t4.  Those
came from `core/amber.py`'s header, which measured them before this box's core split was
known -- and on a machine with four fast cores and four cores at 0.65x, a thread pool of 4
is not four cores' worth of anything.  So the assumption is measured here.

One process, pinned to a set of cores wide enough for the thread count, so that thread
scaling is not confounded with the OS placing threads on slow cores.  Reports the WALL
speedup of the minimisation (which is what threads are for) and the CPU-seconds burned
(which is what the feasibility bound counts), plus the deviation from the golden -- the
thing that decides whether this is a usable arm at all.
"""
import ctypes
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = r"C:\Users\abena\Protein-Folding-Algorithm"
sys.path.insert(0, ROOT)

REF = -489.9138948277905
TARGETS = ("1A13", "1A1P", "1CB3")
FAST = [1, 0, 2, 3]


def set_mask(cores):
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    k32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    mask = 0
    for c in cores:
        mask |= 1 << int(c)
    if not k32.SetProcessAffinityMask(k32.GetCurrentProcess(), mask):
        raise OSError(f"affinity {cores} failed: {ctypes.get_last_error()}")


def main_worker():
    threads = int(sys.argv[2])
    ncores = int(sys.argv[3])
    out = sys.argv[4]
    set_mask(FAST[:ncores])
    os.environ["OMP_NUM_THREADS"] = str(threads)

    import numpy as np
    import protein_geometry as geo
    import peptide_db as db
    import torsion_lib2 as tl2
    import core.amber as A

    rows = []
    for pid in TARGETS:
        p = db.by_pdb(pid)
        tab = tl2.library_for(p.seq, 8, p.seq)
        rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
        rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))
        A.clear_cache()
        H = A.builder_for(p.seq, rep, "CPU", threads)
        c0, w0 = time.process_time(), time.perf_counter()
        o = A._run(H, H._heavy_positions(rb, chi1=None), A.K_MODERATE, 0, 1.0, True)
        cm = o["components"]
        rows.append({"pid": pid, "cpu": time.process_time() - c0,
                     "wall": time.perf_counter() - w0,
                     "interaction": float(cm["nonbonded"] + cm["solvation"]),
                     "threads_prop": H.platform_properties.get("Threads")})
    json.dump({"threads": threads, "rows": rows}, open(out, "w"))


def run(threads, ncores):
    o = os.path.join(HERE, f"thr_{threads}_{ncores}.json")
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    p = subprocess.run([sys.executable, os.path.abspath(__file__), "worker",
                        str(threads), str(ncores), o],
                       env=env, capture_output=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.decode()[:600])
    return json.load(open(o))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        main_worker()
        sys.exit(0)

    base = None
    print(f"  {'threads':>7s} {'cores':>5s} {'wall s/struct':>14s} {'CPU s/struct':>13s} "
          f"{'wall speedup':>13s} {'CPU cost':>9s} {'1A13 d(golden)':>16s}")
    for threads, ncores in ((1, 1), (2, 2), (4, 4)):
        r = run(threads, ncores)
        w = sum(x["wall"] for x in r["rows"]) / len(r["rows"])
        c = sum(x["cpu"] for x in r["rows"]) / len(r["rows"])
        d = next(x["interaction"] for x in r["rows"] if x["pid"] == "1A13") - REF
        if base is None:
            base = (w, c)
        print(f"  {threads:7d} {ncores:5d} {w:14.3f} {c:13.3f} {base[0]/w:13.3f} "
              f"{c/base[1]:9.3f} {d:+16.3e}", flush=True)
