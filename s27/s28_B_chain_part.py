#!/usr/bin/env python
"""s27/s28_B_chain_part.py -- run `s28_B_hop.chain_main` on one half of the targets so the two
halves run as two governed jobs. Part 0 writes the primary `s28_B_chain_rows.jsonl` (and resumes
from it); part 1 writes `s28_B_chain_rows_p1.jsonl`. `s28_B_analyse.py` reads both.

    python s27/s28_B_chain_part.py --part 0|1 --arms A,B,C
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s25 import phys_lib as P              # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

ALL = P.targets()
HALF = len(ALL) // 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", type=int, required=True, choices=(0, 1))
    ap.add_argument("--arms", required=True)
    a = ap.parse_args()
    mine = ALL[:HALF] if a.part == 0 else ALL[HALF:]
    P.targets = lambda: list(mine)                       # chain_main calls P.targets()
    if a.part == 1:
        B.CHAIN_ROWS = os.path.join(B.RESULTS, "s28_B_chain_rows_p1.jsonl")
    print(f"part {a.part}: {len(mine)} targets {mine[0]}..{mine[-1]} -> {B.CHAIN_ROWS}", flush=True)
    B.chain_main([x for x in a.arms.split(",") if x])


if __name__ == "__main__":
    main()
