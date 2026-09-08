"""SPRINT 13 GEO -- build (and cache) the enumerated energy tables.

This is the only expensive step in the study.  Everything else -- metric, gradients,
Hessian, CVaR sweep, optimiser comparison -- is exact numpy over these tables.

THE GRID, AND WHY IT IS SHAPED THIS WAY
    * k = 4 ladder, L = 4..8        -> 8, 10, 12, 14, 16 qubits, peptide length varying
    * k = 2 ladder, L = 8..14       -> 8, 10, 12, 14 qubits at the SAME qubit counts but
                                       DOUBLE the peptide length
    * k = 8 at L = 4, 5             -> 12, 15 qubits, same peptide length, more states
The three ladders cross at 8/10/12/14 qubits, which is what lets "qubit count", "peptide
length" and "states per residue" be separated instead of confounded.

TARGETS: one per leave-fold-out fold, two of them FAIL18 members, all from the 126-target
s12 tuning instrument.  benchmark60 and dev24 are untouched (verified disjoint).

    python -m s13.geo_build [core|extra|all]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import geo_common as G                       # noqa: E402

#: fold 4 / fold 1 / fold 2(FAIL18) -- the full grid runs on these three
CORE = [("1A13", 4), ("1A1P", 1), ("2BFI", 2)]
#: fold 0 / fold 3(FAIL18) / fold 2 -- replication at two sizes only
EXTRA = [("1DEP", 0), ("1ID6", 3), ("1CEK", 2)]

GRID_CORE = ([(L, 4) for L in (4, 5, 6, 7, 8)]
             + [(L, 2) for L in (8, 10, 12)]
             + [(4, 8), (5, 8)])
GRID_EXTRA = [(6, 4), (7, 4)]


def jobs(which: str):
    out = []
    if which in ("core", "all"):
        for pdb, _f in CORE:
            n = len(G.peptide(pdb)["full_seq"])
            for L, k in GRID_CORE:
                if L <= n:
                    out.append((pdb, L, k))
            if n >= 14:
                out.append((pdb, 14, 2))
    if which in ("extra", "all"):
        for pdb, _f in EXTRA:
            for L, k in GRID_EXTRA:
                out.append((pdb, L, k))
    # cheapest first, so the ladder is usable long before the 16-qubit rows land
    return sorted(set(out), key=lambda r: k_states(r))


def I_free():
    from s12 import instrument as I
    return I.free_gb()


def k_states(r):
    _pdb, L, k = r
    return k ** L


def main(which: str = "core"):
    js = jobs(which)
    print(f"{len(js)} (target, L, k) cells; "
          f"{sum(k_states(r) for r in js):,} register states per model", flush=True)
    log = []
    for pdb, L, k in js:
        N = k ** L
        for model in ("legacy", "amber"):
            p = os.path.join(G.CACHE, G._tag(pdb, L, k, model) + ".npy")
            if os.path.exists(p):
                continue
            # Wait for headroom, but do not deadlock on it: this process holds ~0.27 GB
            # (measured, logged as `rss_gb`) against a 1.2 GB cap, so a box that stays under
            # 1.5 GB free is the sibling agents' load. The MemoryError retry below is the
            # real guard -- `core.amber.memory_guard` refuses to build a Context above 92 %
            # and we wait it out there rather than recording a NaN.
            free = G.gate(tries=6 if N < 4096 else 12)
            t0 = time.time()
            for attempt in range(200):
                try:
                    E = G.energy_table(pdb, L, k, model)
                    break
                except MemoryError as exc:
                    # `core.amber.memory_guard` refuses to build an OpenMM Context when the
                    # box is above 92% -- a sibling agent's load, not this process's
                    # (rss_gb below). Wait it out rather than competing for the last page.
                    print(f"  [wait] {exc} (attempt {attempt}, "
                          f"rss={G.rss_gb():.2f} GB)", flush=True)
                    time.sleep(30)
                    free = G.gate(tries=6)
            log.append({"pdb": pdb, "L": L, "k": k, "model": model, "N": int(N),
                        "wall_s": time.time() - t0, "free_gb_before": free,
                        "rss_gb": G.rss_gb(),
                        "min": float(np.min(E)), "max": float(np.max(E)),
                        "median": float(np.median(E)),
                        "n_nonfinite": int((~np.isfinite(E)).sum())})
            G.write("geo_build_log", log)
    G.write("geo_build_log", log)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "core")
