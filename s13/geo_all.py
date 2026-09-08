"""SPRINT 13 GEO -- run every analysis over whatever energy tables currently exist.

One process, sequential, so the compute cap is respected while `s13.geo_build` keeps the
other slot.  Safe to re-run: each stage rewrites its own results file and `geo_pauli`
resumes from what it already computed.

    python -m s13.geo_all [stage ...]
"""
from __future__ import annotations

import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

STAGES = ("anova", "pauli", "metric", "cvar", "land", "grad", "opt", "report")


def main(stages):
    for st in stages:
        t0 = time.time()
        print(f"\n=== {st} ===", flush=True)
        try:
            mod = __import__(f"s13.geo_{st}", fromlist=["main"])
            mod.main()
        except Exception:
            traceback.print_exc()
        print(f"=== {st} done in {time.time() - t0:.0f}s ===", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(STAGES))
