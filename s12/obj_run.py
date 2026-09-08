"""Parallel driver: run a named set of OBJECTIVE CONDITIONS through the production path
on all (or a subset of) the 126 tuning targets.  Each worker handles one target and
evaluates every condition on it, so the pool pack is opened once per worker.

Usage:  python -m s12.obj_run <spec.json>
where spec.json = {"out": name, "targets": [...] | null, "conds": [ {..}, .. ]}
Condition dicts are interpreted by `make_target()` below.
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
import numpy as np

from s12 import instrument as I
from s12 import obj_common as OC


KW = ("L", "cut", "slope", "rank")


def make_target(d, c, rng):
    """Return (dtarget_vector_or_None, weights_or_None, score_override_or_None)."""
    k = c["kind"]
    if k == "bayes":
        return None, None, OC.score_bayes(d["D"], d["risk"], d["grid"])
    if k == "expected":
        return d["exp"], None, None
    if k == "oracle":
        return d["dtrue"], None, None
    if k == "random":
        return None, None, rng.standard_normal(len(d["D"]))
    if k == "alt":
        from s12 import obj_alt as A
        z = np.load(os.path.join(OC.CACHE, "obj_lfo.npz"))
        f = int(d["fold"])
        kw = dict(c.get("kw", {}))
        for nm, key in (("weights", "wsh"), ("recal", "recal"), ("Sig", "Sig")):
            kw[nm] = z[f"{f}/{key}"]
        return None, None, A.SCORERS[c["scorer"]](d, **kw)
    kw = {q: c[q] for q in KW if q in c}
    f = OC.make_corruptor(d, k, rng, **kw)
    if "mae" in c:
        amp = OC.match_amp(f, d["dtrue"], float(c["mae"]))
    else:
        amp = float(c.get("amp", 0.0))
    return f(amp), None, None


def run_target(args):
    pdb, conds, seed0, lam = args
    d = OC.load(pdb)
    res = []
    for ci, c in enumerate(conds):
        rng = np.random.default_rng(abs(hash((pdb, ci, c.get("seed", 0), seed0))) % (2 ** 32))
        dt, w, sc = make_target(d, c, rng)
        if sc is None:
            sc = OC.score_l1(d["D"], dt, w)
            mae, r = OC.mae_r(dt, d["dtrue"])
        else:
            mae, r = OC.mae_r(d["exp"], d["dtrue"]) if c["kind"] == "bayes" else (float("nan"), float("nan"))
        e = OC.emit(d, sc, lam=lam)
        e.pop("sub")
        e.update(mae=mae, r=r, cond=ci, pdb=pdb, tag=c.get("tag", c["kind"]))
        res.append(e)
    return res


def main(spec_path):
    with open(spec_path) as fh:
        spec = json.load(fh)
    tg = OC.targets()
    names = spec.get("targets") or [t["pdb"] for t in tg]
    conds = spec["conds"]
    lam = spec.get("lam", 0.3)
    jobs = [(p, conds, spec.get("seed", 0), lam) for p in names]
    nproc = int(spec.get("nproc", 3))
    import multiprocessing as mp
    with mp.Pool(nproc) as pool:
        out = []
        for k, r in enumerate(pool.imap_unordered(run_target, jobs, chunksize=1)):
            out.extend(r)
            if (k + 1) % 10 == 0:
                print(f"  {k+1}/{len(jobs)}", flush=True)
    I.write(spec["out"], {"spec": spec, "rows": out})
    print("wrote", spec["out"])


if __name__ == "__main__":
    main(sys.argv[1])
