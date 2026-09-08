"""s12 QUANTUM-ROLE -- build and cache the VQ-A instance families.

Usage:  python -m s12.vq_build <family> [n_targets]
Families:
  A2_64   k=2, m=64,  n_ret=600   -> 4096 assemblies, index encoding 12 qubits   (all 126)
  A2_256  k=2, m=256, n_ret=1200  -> 65536 assemblies, index encoding 16 qubits  (all 126)
  A3_32   k=3, m=32,  n_ret=600   -> 32768 assemblies, index encoding 15 qubits  (n>=12)
  A2_8    k=2, m=8,   n_ret=200   -> 64 assemblies,    ONE-HOT 16 qubits         (all 126)
  A3_5    k=3, m=5,   n_ret=200   -> 125 assemblies,   ONE-HOT 15 qubits         (n>=12)
"""
from __future__ import annotations
import os, sys, time, json
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import vq_lib as V

FAMILIES = {
    "A2_64":  dict(k=2, m=64,  n_ret=600,  minlen=4),
    "A2_256": dict(k=2, m=256, n_ret=1200, minlen=4),
    "A3_32":  dict(k=3, m=32,  n_ret=600,  minlen=4),
    "A2_8":   dict(k=2, m=8,   n_ret=200,  minlen=4),
    "A3_5":   dict(k=3, m=5,   n_ret=200,  minlen=4),
    "A2_5":   dict(k=2, m=5,   n_ret=200,  minlen=4),   # ONE-HOT 10 qubits, 25 feasible
}


def build(fam, limit=None, verbose=True):
    from s12 import asm_lib as A
    cfg = FAMILIES[fam]
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    # fold-major order so each (fold, L) bank is loaded once; banks dropped between folds
    tg = sorted(tg, key=lambda t: (t["fold"], t["n"]))
    cur_fold = None
    ok, skip = 0, 0
    t0 = time.time()
    for q, t in enumerate(tg):
        if t["fold"] != cur_fold:
            A.drop_bank()
            cur_fold = t["fold"]
        if V.segments(t["n"], cfg["k"], cfg["minlen"]) is None:
            skip += 1
            continue
        p = V.instance_path(t["pdb"], cfg["k"], cfg["m"], cfg["n_ret"], cfg["minlen"], "fit_ca")
        if os.path.exists(p):
            ok += 1
            continue
        if I.free_gb() < 1.5:
            print(f"  free_gb {I.free_gb():.2f} < 1.5 -- waiting", flush=True)
            while I.free_gb() < 1.5:
                time.sleep(20)
        V.cached_instance(t, **cfg)
        ok += 1
        if verbose and q % 10 == 0:
            print(f"  [{fam}] {q+1}/{len(tg)} built={ok} skip={skip} "
                  f"{time.time()-t0:.0f}s free={I.free_gb():.2f}", flush=True)
    print(f"[{fam}] done built={ok} skipped={skip} {time.time()-t0:.0f}s", flush=True)
    return ok, skip


if __name__ == "__main__":
    fam = sys.argv[1]
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    build(fam, lim)
