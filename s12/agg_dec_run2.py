"""Driver for the v2 (signed-deviation, FAIL18-out-of-training) decoder configurations."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import agg_decoder as D

JOBS = [
    ("v2", "v2", None, 0), ("v2_abs", "v2_abs", None, 0), ("v2_null", "v2_null", None, 0),
    ("v2_oracle", "v2_oracle", None, 0), ("v2_s1", "v2", None, 1),
    ("v2_n8", "v2", 8, 0), ("v2_n16", "v2", 16, 0), ("v2_n32", "v2", 32, 0),
    ("v2_n64", "v2", 64, 0), ("v2_101", "v2_101", None, 0), ("v2_hi", "v2_hi", None, 0),
    ("v2_oracle_n8", "v2_oracle", 8, 0), ("v2_oracle_n32", "v2_oracle", 32, 0),
]
if __name__ == "__main__":
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    for tag, cfg, nt, sd in JOBS:
        if only and tag not in only:
            continue
        if os.path.exists(os.path.join(D.ROOT, "s12", "results", f"agg_dec_{tag}.json")):
            print("skip", tag, flush=True); continue
        try:
            D.main(tag, dict(D.CFGS[cfg]), ntrain=nt, seed=sd, verbose=False)
        except Exception as ex:
            import traceback; traceback.print_exc(); print("FAIL", tag, ex, flush=True)
