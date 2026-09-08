"""s19/agentC_kv.py -- Q1: WHY DOES SCORE-ORDERING DAMAGE AN AVERAGED SET?

Pre-registration: `s19/PREREG_C.md` section 3, written before this module produced a number.

THE QUESTION SPRINT 18 LEFT OPEN.  Every physics filter, placed in front of the operator that
emits the answer, loses to a random gate of the same size (+0.036 to +0.063 A, five scores, every
CI excluding zero).  The recorded operator law

    d_out = 1.16 * d_set_mean + 0.04 * d_set_best

HOLDS for a matched-random gate (miss +0.006 [-0.003, +0.014]) and FAILS for every physics gate
(legacy +0.094 [+0.063, +0.132], leg_torsion +0.079 [+0.060, +0.098]) -- with the law's PREMISE
holding exactly (Legacy improved the set mean 3.551 -> 3.527) and the output still getting worse.

So something a score gate does to a candidate set is invisible to its mean and its best.  The
hypothesis nobody measured is the SURVIVORS' ERROR COVARIANCE: a score gate may keep candidates
whose errors point the same way, so the average no longer cancels them.

THE INSTRUMENT IS AN IDENTITY, NOT A MODEL.  For an average of `m` members in a COMMON FRAME,
the Krogh-Vedelsby / parallel-axis decomposition is exact:

    readout^2 = E_mem^2 - D^2                 and     readout^2 = <d_reb^2> + FRAME^2 - D^2

so the output error of the averaging operator is FULLY accounted for by three measurable
channels, of which the operator law can see only the first.  See `s19/agentC_lib` for the
derivation, the sign conventions, and the two gates that verify both.

WHAT IS ORACLE HERE.  `readout`, `E_mem`, `FRAME` and every per-member error need the native and
are labelled ORACLE diagnostics.  `D` and every gate in this module are native-free.

    python -m s19.agentC_kv --smoke          # 6 targets, writes _SMOKE_*, never a result
    python -m s19.agentC_kv                  # n = 126
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
from s18 import phys_lib as PL                                    # noqa: E402
from s19 import agentC_lib as CL                                  # noqa: E402

F_PRIMARY = 0.50          # pre-registered: keep m = 37 of K = 75
N_RAND = 30               # matched-random draws (Sprint 18 used 3; this tightens the null)

#: The score gates.  `leg_*` are single genuine Legacy components; `legacy` is the eleven-term
#: total at DEFAULT_WEIGHTS, never fitted; `amber_sp` is genuine ff14SB/GBn2, read from the
#: Sprint-18 artefact (gate GC1 verifies that artefact against a live recomputation).
SCORE_ARMS = ("legacy", "leg_torsion", "leg_contact", "leg_steric", "amber_sp", "helix")
DIV_ARMS = ("legacy_spread", "legacy_clust")


def _load_amber_sp():
    """genuine AMBER single points, per candidate, from the Sprint-18 Phase-8 artefact."""
    o = PL.read_complete(os.path.join(PL.RESULTS, "down.json"), need=126)
    return {r["pdb"]: np.asarray(r["s_amber_sp"], float) for r in o["rows"]}


def run_target(t, amber_sp):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(W)
    m = max(2, int(round(K * (1.0 - F_PRIMARY))))

    d_win = np.asarray(I.kabsch_rmsd_batch(W, nat), float)          # ORACLE labels
    comp, legacy = CL.legacy_scores(seq, PHI, PSI)
    scores = {"legacy": legacy,
              "leg_torsion": np.asarray(comp["torsion"], float),
              "leg_contact": np.asarray(comp["contact"], float),
              "leg_steric": np.asarray(comp["steric"], float),
              "helix": CL.helix_score(W, n)}
    if amber_sp is not None:
        scores["amber_sp"] = np.asarray(amber_sp, float)

    P = np.asarray(I.pairwise_rmsd(W), float)                       # native-free, for the gates
    rng = SD.stable_rng(pdb, "s19C_kv")

    idxs = {"none": [np.arange(K)]}
    for nm in SCORE_ARMS:
        if nm in scores:
            idxs[nm] = [CL.keep_lowest(scores[nm], m)]
    idxs["rand"] = [CL.gate_random(K, m, rng) for _ in range(N_RAND)]
    idxs["legacy_spread"] = [CL.gate_score_then_spread(scores["legacy"], P, m, pre=2.0)]
    idxs["legacy_clust"] = [CL.gate_cluster_best(scores["legacy"], P, m,
                                                 SD.stable_rng(pdb, "s19C_clust"))]

    arms = {}
    for nm, lst in idxs.items():
        sub = []
        for idx in lst:
            Y, Ybar, _C, _dev = CL.common_frame_members(W[idx], PHI[idx], PSI[idx])
            k = CL.kv_of(Y, nat, d_win=d_win[idx])
            #: the Sprint-18 diversity column, on the RAW windows, carried for continuity with
            #: `s18.phys_lib.diversity` -- a different object from the KV `D`, never conflated.
            k["div_s18"] = PL.diversity(W[idx])
            k["idx_hash"] = int(np.sum(idx * np.arange(1, len(idx) + 1)))
            sub.append(k)
        if len(sub) == 1:
            arms[nm] = sub[0]
        else:
            #: a matched-random arm is the MEAN over its draws.  Averaging the SQUARED channels
            #: keeps (KV) exact for the averaged arm, which is what the primary statistic needs;
            #: the linear readout is carried beside it for comparison with the Sprint-18 table.
            agg = {k: float(np.mean([s[k] for s in sub]))
                   for k in sub[0] if isinstance(sub[0][k], float)}
            agg["m"] = int(sub[0]["m"])
            agg["kv_resid"] = float(abs(agg["readout2"] - (agg["E_mem2"] - agg["D2"])))
            agg["n_draws"] = len(sub)
            agg["draws"] = sub
            arms[nm] = agg
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]), "K": int(K), "m": int(m),
            "pool_mean": float(d_win.mean()), "pool_best": float(d_win.min()),
            "arms": arms}


def run(targets=None, out="agentC_kv.json", verbose=True):
    tg = targets if targets is not None else I.targets()
    asp = _load_amber_sp()
    path = os.path.join(CL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        import json
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"F_PRIMARY": F_PRIMARY, "N_RAND": N_RAND, "SCORE_ARMS": list(SCORE_ARMS),
           "DIV_ARMS": list(DIV_ARMS), "amber_sp_source": "s18/results/down.json"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(1.0)
        rows.append(run_target(t, asp.get(t["pdb"])))
        if verbose:
            a = rows[-1]["arms"]
            print(f"  {len(rows)}/{len(tg)} {t['pdb']}  none={a['none']['readout']:.3f} "
                  f"legacy={a['legacy']['readout']:.3f} rand={a['rand']['readout']:.3f} "
                  f"D(leg)={a['legacy']['D']:.3f} D(rand)={a['rand']['D']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg))
    _write(out, rows, cfg, len(tg))
    return rows


def _write(out, rows, cfg, n_expected):
    import json
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected))}
    p = os.path.join(CL.RESULTS, out if out.endswith(".json") else out + ".json")
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)
    return p


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="agentC_kv.json")
    a = ap.parse_args()
    tg = I.targets()
    out = a.out
    if a.smoke:
        tg, out = tg[:6], "_SMOKE_agentC_kv.json"
    elif a.limit:
        tg, out = tg[:a.limit], "_PARTIAL_" + a.out
    run(tg, out=out)
