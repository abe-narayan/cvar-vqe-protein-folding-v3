"""S15 AUDIT / Part A+D -- rerun the PRODUCTION pipeline from source for a few targets
and compare to the config-keyed cache the instrument reads.

`s12.instrument` reads `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` for `sub` and
`fit_ca`; two of the five pinned constants (`top75_best`, `synthesis_fit`) and one third
(`n_zero_recall`) are therefore CACHE READS.  This module recomputes them.

    python -m s15.audit_pipe run 1A13,1A1P,1CB3,2MK7     # rerun + compare
    python -m s15.audit_pipe amber                        # AMBER cost + determinism
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OPENMM_CPU_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)


def _w(name, obj):
    p = os.path.join(OUT, name if name.endswith(".json") else name + ".json")
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("wrote", p, flush=True)


def run(pdbs):
    import core
    from core import pipeline as pl
    from s12 import instrument as I
    from s15.audit_rmsd import horn_rmsd
    db = core.backend("data")
    cfg = pl.PROD
    key = cfg.key()
    tg = {p.pdb: p for p in db.load()}
    folds = db.folds(cfg.n_folds)
    rows = []
    for pdb in pdbs:
        t = tg[pdb]
        fold = int(folds[t.seq])
        pl.guard_esm([t.seq])
        t0 = time.perf_counter()
        rec, pool, clk = pl.run_target(t, fold, cfg)
        dt = time.perf_counter() - t0
        cached = I.shipped_record(pdb)
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        sub_new = np.asarray(rec["sub"], int)
        sub_old = np.asarray(cached["sub"], int)
        fit_new = np.asarray(rec["fit_ca"], float)
        fit_old = np.asarray(cached["fit_ca"], float)
        row = {"pdb": pdb, "n": int(t.n), "fold": fold, "cfg_key_now": key,
               "cfg_key_cached": cached.get("cfg_key"), "seconds": round(dt, 1),
               "stage_seconds": {k: round(v, 2) for k, v in clk.t.items()},
               "sub_identical": bool(np.array_equal(sub_new, sub_old)),
               "sub_n_diff": int((sub_new != sub_old).sum()) if sub_new.shape == sub_old.shape else -1,
               "fit_ca_max_abs_diff": float(np.abs(fit_new - fit_old).max()),
               "fit_ca_bit": bool(np.array_equal(fit_new, fit_old)),
               "synthesis_fit_cached": horn_rmsd(fit_old, nat),
               "synthesis_fit_rerun": horn_rmsd(fit_new, nat),
               "top75_best_cached": float(u["rr"][I.pool_idx(u)][sub_old].min()),
               "top75_best_rerun": float(u["rr"][I.pool_idx(u)][sub_new].min())}
        row["synthesis_fit_delta"] = row["synthesis_fit_rerun"] - row["synthesis_fit_cached"]
        a = rec.get("amber") or {}
        row["amber"] = {k: a.get(k) for k in ("e0", "e1", "moved", "strain_after", "err")}
        ca_new = np.asarray(rec["ca"], float)
        row["ca_max_abs_diff"] = float(np.abs(ca_new - np.asarray(cached["ca"], float)).max())
        if cached.get("amber_ca") is not None and a.get("ca") is not None:
            row["amber_ca_max_abs_diff"] = float(np.abs(np.asarray(a["ca"], float)
                                                        - np.asarray(cached["amber_ca"], float)).max())
            row["amber_e1_cached"] = cached.get("amber_e1")
        rows.append(row)
        print(json.dumps(row), flush=True)
        _w("audit_pipe_run", {"rows": rows, "n_expected": len(pdbs),
                              "complete": len(rows) == len(pdbs)})
    _w("audit_pipe_run", {"rows": rows, "n_expected": len(pdbs),
                          "complete": len(rows) == len(pdbs)})


# ------------------------------------------------------------------- AMBER cost
def amber(pdb="1A13", n_rep=40):
    """Distinct AMBER single-point cost, and whether repeats hit a memo.

    The 6 ms figure that was live in a shared brief was wrong because the timing loop
    re-scored the SAME structure and hit a result memo.  Here every call gets a distinct
    perturbation, and the memoised timing is measured beside it so the gap is visible.
    """
    import core
    from s12 import instrument as I
    from s13 import qarch_lib as ql
    amb = core.backend("amber")
    u = I.load_univ(pdb)
    p = I.pool_idx(u)
    PHI = np.asarray(u["PHI"], float)[p][:n_rep]
    PSI = np.asarray(u["PSI"], float)[p][:n_rep]
    seq = u["seq"]
    out = {"pdb": pdb, "n": int(u["n"]), "seq": seq, "n_rep": int(n_rep),
           "amber_module": getattr(amb, "__name__", str(amb))}
    fns = [f for f in dir(amb) if "energ" in f.lower() or "single" in f.lower()]
    out["candidate_entrypoints"] = fns
    # distinct structures
    fn = None
    for name in ("single_point", "energy", "energies", "score"):
        if hasattr(amb, name):
            fn = getattr(amb, name); out["entrypoint"] = name; break
    if fn is None:
        out["error"] = "no obvious single-point entrypoint; see candidate_entrypoints"
        _w("audit_pipe_amber", out); return out
    t0 = time.perf_counter()
    vals = [float(fn(seq, PHI[k], PSI[k])) for k in range(n_rep)]
    out["ms_per_distinct_call"] = 1000.0 * (time.perf_counter() - t0) / n_rep
    t0 = time.perf_counter()
    for _ in range(n_rep):
        fn(seq, PHI[0], PSI[0])
    out["ms_per_repeat_call"] = 1000.0 * (time.perf_counter() - t0) / n_rep
    out["values_head"] = vals[:5]
    # determinism: same input twice
    out["deterministic"] = float(fn(seq, PHI[0], PSI[0])) == vals[0]
    _w("audit_pipe_amber", out)
    return out


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    arg = sys.argv[2] if len(sys.argv) > 2 else None
    if cmd == "run":
        run((arg or "1A13,1A1P,1CB3,2MK7").split(","))
    elif cmd == "amber":
        amber(arg or "1A13")
