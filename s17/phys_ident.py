"""s17/phys_ident.py -- THE IDENTICAL-CANDIDATE INSTRUMENT.

One candidate set per target; every physical model scored on the SAME structures; then each
model evaluated in all four roles the sprint asks about -- SELECTOR, GATE, RANK FEATURE,
REPAIR OPERATOR.  Sprint 17 sections 26, 27, 32, 34.  Pre-registration in PREREG_phys.md
(E2, E3, E4, E5), written before this module produced a number.

THE CANDIDATE SET.  The shipped K = 500 BLOSUM62 pool per target, in retrieval order.  This
is a **native-free, ranker-neutral** set: it is fixed by the retrieval stage, not by any of
the three scores being compared, so no ranker is handed a set built by its own criterion.
The K = 500 truncation's cost against the 13,000-27,000-window universe is read from the
coordinator's `s17/results/oracle_map.json` and printed BEFORE any physics number.

WHAT IS COMPUTED ON EVERY CANDIDATE (all native-free except `d`):

    d            CA-RMSD to native                       ORACLE -- labels and ceilings only
    dist         the shipped Bayes-risk distogram score  NATIVE-FREE
    legacy       the genuine eleven-component total at DEFAULT_WEIGHTS (never fitted)
    legacy_*     each of the eleven components separately
    amber        genuine ff14SB/GBn2 single point on the ideal-geometry rebuild
    amber_*      bond / angle / torsion / nonbonded / solvation

Legacy and AMBER are read on the SAME ideal-geometry rebuild of the window's (phi, psi), so
the two potentials see one structure and the comparison is exact.

MONOTONE CONDITIONING.  Raw AMBER single points on unrelaxed structures are catastrophically
heavy-tailed (`s16` §2.3: 35.8% above 1e6 kcal/mol, top-10 carrying 99.9% of the variance).
Every statistic here is RANK-based (Spearman, AUROC, argmin, percentile), which
`s16/energy_lib.condition` verified is invariant to the conditioning map at rho = 1.000000
exactly.  No moment of a raw AMBER energy is taken anywhere in this module.

TIE-BREAKING.  `argmin` on a tied score reads the cache's sort order, which on this project
is sometimes the ORACLE order (a recorded trap that once invented a 1.386 A winner).  Every
selection here goes through `energy_lib.argmin_tied`, which averages the outcome over the
full tied argmin set.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                     # noqa: E402
from s15 import seed as SD                          # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s17 import phys_lib as P                       # noqa: E402

K_POOL = 500                    # the candidate set width
AMBER_M = 500                   # AMBER single points per target (0 = skip)
NEAR = 2.0                      # "near-native" threshold, angstrom
GATE_F = (0.10, 0.25, 0.50)     # gate rejection fractions
N_RAND = 3                      # matched-random draws


def _cols(row):
    """The score columns available on a per-target row."""
    return [k for k in row if k.startswith("s_")]


# ------------------------------------------------------------------ per target
def run_target(t, want_amber=True):
    from core import geometry as geo
    import torsion_lib2 as tl2

    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    u = I.load_univ(pdb)
    order = np.asarray(u["order"], int)
    k = min(K_POOL, len(order))
    sel = order[:k]
    W = np.asarray(u["W"], float)[sel]
    PHI = np.asarray(u["PHI"], float)[sel]
    PSI = np.asarray(u["PSI"], float)[sel]
    nat = np.asarray(u["nat_ca"], float)

    #: ORACLE -- labels and ceilings only, never a predictor.
    d = I.kabsch_rmsd_batch(W, nat)

    #: NATIVE-FREE scores.
    i, j = I.pair_index(n)
    dg = I.distogram(pdb, seq, fold)
    dist = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)

    comp = EL.legacy_components_of_windows(seq, PHI, PSI)
    legacy = np.asarray(EL.legacy_total_from(comp), float)

    row = {"pdb": pdb, "n": n, "fold": fold, "k": int(k),
           "d": d.tolist(), "s_dist": dist.tolist(), "s_legacy": legacy.tolist()}
    for term in EL.LEG_TERMS:
        row["s_leg_" + term] = np.asarray(comp[term], float).tolist()

    if want_amber:
        tab = tl2.library_for(seq, 4, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        box = P.ConstrainedBox(seq, rep)
        m = min(AMBER_M, k)
        E = np.full(k, np.nan)
        Et = {c: np.full(k, np.nan) for c in
              ("bond", "angle", "torsion", "nonbonded", "solvation")}
        t0 = time.time()
        for b in range(m):
            c = geo.build_backbone(PHI[b], PSI[b])
            e = box.energy_point({a: np.asarray(v, float) for a, v in c.items()},
                                 components=True)
            E[b] = e["energy"]
            for cn in Et:
                Et[cn][b] = e[cn]
        row["amber_m"] = int(m)
        row["amber_wall"] = round(time.time() - t0, 2)
        row["s_amber"] = E.tolist()
        for cn, v in Et.items():
            row["s_amb_" + cn] = v.tolist()
        box.close()
    return row


def run(targets=None, out="phys_ident.json", want_amber=True, verbose=True):
    tg = targets if targets is not None else I.targets()
    path = os.path.join(RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path))["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for c, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        EL.mem_hold(1.6, tag="ident")
        rows.append(run_target(t, want_amber=want_amber))
        if verbose:
            print(f"  {len(rows)}/{len(tg)} {t['pdb']} ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"rows": rows, "k_pool": K_POOL, "amber_m": AMBER_M},
                  open(path, "w"))
    json.dump({"rows": rows, "k_pool": K_POOL, "amber_m": AMBER_M}, open(path, "w"))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-amber", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="phys_ident.json")
    a = ap.parse_args()
    tg = I.targets()
    if a.limit:
        tg = tg[:a.limit]
    run(tg, out=a.out, want_amber=not a.no_amber)
