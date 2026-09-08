"""forensics: build the per-target caches (parent map, window SS x2, ESM key) for tuning126."""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from s12 import instrument as I
from s12 import forensics_lib as L

if __name__ == "__main__":
    tg = I.targets()
    t0 = time.time()
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        L.parent_map(pdb, u); L.window_ss(pdb, u); L.window_cass(pdb, u); L.esm_key_cached(pdb, u)
        print(f"[{k+1:3d}/126] {pdb} nw={len(u['W'])} {time.time()-t0:.0f}s free={I.free_gb():.1f}", flush=True)
    print("done", flush=True)
