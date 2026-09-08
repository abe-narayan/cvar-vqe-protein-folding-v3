"""s19/agentC_reject.py -- Q3: CAN LEGACY REJECT A CLASS OF IMPOSSIBLE CANDIDATE?

Pre-registration: `s19/PREREG_C.md` section 5, written before this module produced a number.

THE ONLY LEGACY ROLE NOT YET REFUTED, and it is narrow.  Sprint 18 closed Legacy as a ranker and
as an objective term, and closed `leg_torsion` as a gate.  What it did NOT test is whether Legacy
can reject a small number of candidates that are not merely bad but PHYSICALLY IMPOSSIBLE -- a
hard steric violation is a different claim from "this one scores worse".

THE RULE IS ABSOLUTE (BRIEF section 8): if Legacy rejects r candidates, the comparison is against
rejecting r candidates AT RANDOM.  And, because Sprint 19's own section-3 result says the damage a
gate does is carried by the DIVERSITY it destroys, a second control is mandatory here:

    rand@r      MATCHED-RANDOM -- r rejected at random, 30 `stable_rng` draws
    divkeep@r   DIVERSITY-PRESERVING -- reject the r most REDUNDANT candidates, i.e. the r whose
                nearest neighbour in the surviving set is closest.  This is the rejection rule
                that costs the ensemble the least ambiguity, and it uses no score at all.

A detector must beat BOTH or the role closes.

The detectors, all native-free:

    steric@r    the r worst by genuine Legacy `steric` (DEFAULT_WEIGHTS, never fitted)
    minheavy@r  the r worst by minimum heavy-atom distance (`s16.energy_lib.panel`) -- an
                independent, force-field-free statement of the same physical claim
    amber@r     the r worst by genuine ff14SB/GBn2 single point

    python -m s19.agentC_reject
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s16 import energy_lib as EL                                  # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402
from s19 import agentC_lib as CL                                  # noqa: E402
from s19 import agentC_kv as KV                                   # noqa: E402

R_LIST = (2, 5, 10, 19)          # 10 is the PRE-REGISTERED PRIMARY
N_RAND = 30


def _reject_worst(score, r):
    """Drop the `r` HIGHEST-scoring (worst) candidates; deterministic index tie-break."""
    s = np.asarray(score, float)
    s = np.where(np.isfinite(s), s, np.inf)
    order = np.lexsort((np.arange(len(s)), -s))
    return np.sort(order[r:])


def _reject_redundant(P, r):
    """Drop the `r` most REDUNDANT candidates -- greedy, each step removes the member whose
    nearest surviving neighbour is closest.  Native-free, score-free, diversity-preserving."""
    P = np.asarray(P, float).copy()
    K = len(P)
    alive = np.ones(K, bool)
    np.fill_diagonal(P, np.inf)
    for _ in range(r):
        sub = np.where(alive[None, :], P, np.inf)
        nn = sub[alive].min(1)
        idx = np.flatnonzero(alive)
        alive[idx[int(np.argmin(nn))]] = False
    return np.sort(np.flatnonzero(alive))


def _panel_minheavy(seq, PHI, PSI):
    """Minimum heavy-atom distance of every candidate's ideal-geometry rebuild."""
    from core import geometry as geo
    c = geo.build_backbone_batch(np.asarray(PHI, float), np.asarray(PSI, float))
    out = np.empty(len(PHI))
    for b in range(len(PHI)):
        bb = {a: np.asarray(c[a], float)[b] for a in ("N", "CA", "C", "O", "CB") if a in c}
        out[b] = EL.panel(bb, seq)["min_heavy"]
    return out


def run_target(t, amber_sp):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(W)
    d_win = np.asarray(I.kabsch_rmsd_batch(W, nat), float)
    comp, legacy = CL.legacy_scores(seq, PHI, PSI)
    mh = _panel_minheavy(seq, PHI, PSI)
    P = np.asarray(I.pairwise_rmsd(W), float)
    rng = SD.stable_rng(pdb, "s19C_reject")

    #: worse = higher.  `min_heavy` is the other way round, so it is negated.
    det = {"steric": np.asarray(comp["steric"], float),
           "minheavy": -np.asarray(mh, float),
           "legacy": np.asarray(legacy, float)}
    if amber_sp is not None:
        det["amber"] = np.asarray(amber_sp, float)

    idxs = {"none": [np.arange(K)]}
    for r in R_LIST:
        for nm, s in det.items():
            idxs[f"{nm}@{r}"] = [_reject_worst(s, r)]
        idxs[f"divkeep@{r}"] = [_reject_redundant(P, r)]
        idxs[f"rand@{r}"] = [np.sort(rng.permutation(K)[r:]) for _ in range(N_RAND)]

    arms = {}
    for nm, lst in idxs.items():
        sub = []
        for idx in lst:
            Y, _Yb, _C, _dv = CL.common_frame_members(W[idx], PHI[idx], PSI[idx])
            k = CL.kv_of(Y, nat, d_win=d_win[idx])
            sub.append(k)
        if len(sub) == 1:
            arms[nm] = sub[0]
        else:
            agg = {k: float(np.mean([s[k] for s in sub]))
                   for k in sub[0] if isinstance(sub[0][k], float)}
            agg["m"] = int(sub[0]["m"]); agg["n_draws"] = len(sub)
            arms[nm] = agg
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]), "K": int(K),
            "min_heavy_pool_min": float(mh.min()), "min_heavy_pool_mean": float(mh.mean()),
            "n_impossible_2p6": int((mh < 2.6).sum()), "n_impossible_2p0": int((mh < 2.0).sum()),
            "arms": arms}


def run(targets=None, out="agentC_reject.json", verbose=True):
    import json
    tg = targets if targets is not None else I.targets()
    asp = KV._load_amber_sp()
    path = os.path.join(CL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"R_LIST": list(R_LIST), "N_RAND": N_RAND, "PRIMARY_r": 10,
           "amber_sp_source": "s18/results/down.json"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.8)
        rows.append(run_target(t, asp.get(t["pdb"])))
        if verbose and len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
        KV._write(out, rows, cfg, len(tg))
    KV._write(out, rows, cfg, len(tg))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    run(I.targets()[:4] if a.smoke else None,
        out="_SMOKE_agentC_reject.json" if a.smoke else "agentC_reject.json")
