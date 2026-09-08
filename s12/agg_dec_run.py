"""Sequential driver for the set-decoder configurations."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import agg_decoder as D

JOBS = [
    ("listwise", "listwise", None, 0),
    ("regress", "regress", None, 0),
    ("null", "null", None, 0),
    ("oracle", "oracle", None, 0),
    ("weight101", "weight101", None, 0),
    ("weight_s1", "weight", None, 1),
    ("weight_n8", "weight", 8, 0),
    ("weight_n16", "weight", 16, 0),
    ("weight_n32", "weight", 32, 0),
    ("weight_n64", "weight", 64, 0),
    ("oracle_n8", "oracle", 8, 0),
    ("oracle_n32", "oracle", 32, 0),
    ("weight_hi", "weight_hi", None, 0),
    ("weight_noattn", "weight_noattn", None, 0),
    ("null_A", "null_A", None, 0),
    ("weight_p_only", "weight_p_only", None, 0),
]

if __name__ == "__main__":
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    for tag, cfgname, nt, sd in JOBS:
        if only and tag not in only:
            continue
        if os.path.exists(os.path.join(D.ROOT, "s12", "results", f"agg_dec_{tag}.json")):
            print("skip", tag, flush=True); continue
        try:
            D.main(tag, dict(D.CFGS[cfgname]), ntrain=nt, seed=sd, verbose=False)
        except Exception as ex:
            print("FAIL", tag, ex, flush=True)
