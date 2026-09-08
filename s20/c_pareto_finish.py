"""s20/c_pareto_finish.py -- FINISH the Sprint-19 Q2 AMBER Pareto run, in Sprint 20's own directory.

Sprint 19 left `s19/results/agentC_pareto.json` at **28 of 126 rows**, `complete: false`, running
the FULL pre-registered ladder UNCAPPED (`STEPS = 0`, the deployed protocol) after the
"pathological minimiser" claim was RETRACTED (s19 CLAIMS Z4).  My brief requires me to finish or
supersede it before doing anything new.

WHAT THIS MODULE DOES AND DOES NOT DO.

  * It runs `s19.agentC_pareto.run_target` UNCHANGED -- the same arms, the same restraints, the
    same `_run_ref`, the same `EL.panel` validity vector, the same gates.  No arm is added,
    removed or re-parameterised here.  This is a CONTINUATION, not a new experiment.
  * It writes to `s20/results/c_pareto.json`, NOT to s19's file.  Sprint 19's artefact is left
    exactly as its author left it (brief section 5: never overwrite a prior sprint's files).
  * It SEEDS itself from s19's 28 rows only after checking that s19's persisted `cfg_hash`
    equals the hash the CURRENT module config produces.  If the configs differ the seed is
    refused and the run starts from zero -- rows made by a different operator are never merged.

Usage:
    python -u -m s20.c_pareto_finish --gate     # GC2/GC3 re-run, into s20/results
    python -u -m s20.c_pareto_finish            # continue to n = 126, resumable per target
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I          # noqa: E402
from s18 import phys_lib as PL           # noqa: E402
from s19 import agentC_pareto as CP      # noqa: E402

OUT = os.path.join(RESULTS, "c_pareto.json")
S19 = os.path.join(ROOT, "s19", "results", "agentC_pareto.json")


def _cfg():
    return {"K_LADDER": list(CP.K_LADDER), "STEPS": CP.STEPS, "TOL": CP.TOL,
            "BLEND": list(CP.BLEND),
            "SETS": {k: list(v) for k, v in CP.SETS.items()},
            "K_INCUMBENT": CP.K_INCUMBENT}


def _write(rows, cfg, n_expected, seeded_from=None):
    fn = [r["frame_null"]["d_rmsd"] for r in rows if "frame_null" in r]
    fe = [r["frame_null"]["d_energy"] for r in rows if "frame_null" in r]
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected)),
           "seeded_from": seeded_from,
           "module": "s20.c_pareto_finish (continuation of s19.agentC_pareto, code unchanged)"}
    if fn:
        obj["frame_null"] = {"n": len(fn), "max_abs": float(np.max(np.abs(fn))),
                             "mean_abs": float(np.mean(np.abs(fn))),
                             "max_abs_energy": float(np.max(np.abs(fe)))}
    tmp = OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, OUT)


def seed_rows(cfg, have=()):
    """s19's rows, but ONLY those whose persisted `cfg_hash` matches this module's config.

    A SECOND WORKSTREAM IS RUNNING THE IDENTICAL EXPERIMENT.  Mid-sprint I found another lane's
    `python -m s19.agentC_pareto` writing `s19/results/agentC_pareto.json` while my continuation
    was writing `s20/results/c_pareto.json` -- the same code, the same config hash, the same
    targets, twice, on a box the brief says to keep under 97%.  I stopped MY copy rather than
    theirs (their file is the pre-registration's own artefact and mine is the duplicate) and
    merge their finished rows here instead.  Both runs are deterministic and share a cfg_hash,
    so a merged row is bit-identical to the one I would have computed; the merge is recorded in
    `seeded_from` with the row count taken from each source, and s19's file is never written to.
    """
    if not os.path.exists(S19):
        return [], None
    try:
        o = json.load(open(S19))
    except Exception:
        return [], None
    want = PL.cfg_hash(cfg)
    if o.get("cfg_hash") != want:
        print(f"REFUSING s19 merge: cfg_hash {o.get('cfg_hash')} != {want}", flush=True)
        return [], None
    have = set(have)
    new = [r for r in o.get("rows", []) if r["pdb"] not in have]
    if new:
        print(f"merged {len(new)} rows from s19/results/agentC_pareto.json "
              f"(cfg_hash {want} matches; it now holds {o.get('n_rows')})", flush=True)
    return new, {"path": "s19/results/agentC_pareto.json", "n": len(o.get("rows", [])),
                 "cfg_hash": want, "note": "another lane ran the identical experiment; "
                                           "deterministic, same cfg_hash, merged not recomputed"}


def main(limit=0):
    cfg = _cfg()
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    rows, seeded = [], None
    if os.path.exists(OUT):
        o = json.load(open(OUT))
        if o.get("cfg_hash") == PL.cfg_hash(cfg):
            rows, seeded = o.get("rows", []), o.get("seeded_from")
            print(f"resuming from {len(rows)} rows in s20/results/c_pareto.json", flush=True)
        else:
            print("existing s20 file has a different cfg_hash -- refusing to resume", flush=True)
            sys.exit(1)
    new, s2 = seed_rows(cfg, have={r["pdb"] for r in rows})
    if new:
        rows.extend(new)
        seeded = s2
        _write(rows, cfg, len(tg), seeded)
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        #: re-check the other lane before spending ten AMBER minimisations on a target it may
        #: have finished while this loop was running.
        new, s2 = seed_rows(cfg, have=done)
        if new:
            rows.extend(new)
            seeded = s2
            done |= {r["pdb"] for r in new}
            _write(rows, cfg, len(tg), seeded)
            if t["pdb"] in done:
                continue
        PL.mem_ok(0.55)
        # the frame null is run on the first N_FRAME_NULL rows; s19's seed already carries 15.
        nfn = sum(1 for r in rows if "frame_null" in r)
        rows.append(CP.run_target(t, do_frame_null=(nfn < CP.N_FRAME_NULL)))
        a = rows[-1]["arms"]
        print(f"  {len(rows)}/{len(tg)} {t['pdb']} in={rows[-1]['input_rmsd']:.3f} "
              f"k1={a['k1']['rmsd']:.3f} k30={a['k30']['rmsd']:.3f} "
              f"k300={a['k300']['rmsd']:.3f} ca30={a['caonly_k30']['rmsd']:.3f} "
              f"pb={a['pullback_30_1000']['rmsd']:.3f} ({time.time()-t0:.0f}s)", flush=True)
        _write(rows, cfg, len(tg), seeded)
    _write(rows, cfg, len(tg), seeded)
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.gate:
        #: GC3 (the iteration bound is inert) is **NOT RUN and is recorded as INAPPLICABLE**.
        #: This run uses `STEPS = 0`, the deployed unbounded protocol, so there is NO bound to
        #: certify: GC3 at cap = 0 would compare uncapped against uncapped and pass at
        #: 0.00e+00 with ZERO firings.  Sprint 19's own rule Z6 -- "a pass with zero firings is
        #: not evidence" -- forbids quoting that as a gate.  It is declared here as N/A rather
        #: than run vacuously.
        g2 = CP.gate_GC2()
        with open(os.path.join(RESULTS, "c_gates_amber.json"), "w") as fh:
            json.dump({"GC2": g2,
                       "GC3": {"status": "INAPPLICABLE",
                               "reason": "STEPS=0 (unbounded, deployed protocol): no bound "
                                         "exists to certify; a cap=0 GC3 would pass with zero "
                                         "firings, which s19 Z6 forbids quoting"},
                       "complete": True, "passed": bool(g2["passed"])}, fh, indent=1)
        sys.exit(0 if g2["passed"] else 1)
    main(a.limit)
