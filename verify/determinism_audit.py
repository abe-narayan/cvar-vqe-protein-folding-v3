"""Determinism, resumability, and cache-key completeness.

  D1  same target, two SEPARATE processes, bit-identical deployable record
  D2  a run killed mid-way and resumed produces the identical final output
  D3  cache keys discriminate every type that JSON alone would conflate
  D4  a forged sidecar (the signature of a key missing a parameter) makes `load` RAISE
      `CacheCollision` rather than serve a stale array
  D5  the config key covers every science parameter, and the backend set
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

RUNNER = r'''
import json, os, sys
sys.path.insert(0, %r)
import numpy as np
import core
from core import pipeline as P
db = core.backend("data")
folds = db.folds(P.PROD.n_folds)
t = [x for x in P.manifest("smoke8") if x.pdb == %r][0]
rec, _, _ = P.run_target(t, int(folds[t.seq]), P.PROD)
out = {}
for k, v in rec.items():
    if isinstance(v, np.ndarray):
        out[k] = [repr(float(x)) for x in np.asarray(v, float).ravel()]
    elif isinstance(v, dict):
        for k2, v2 in v.items():
            if isinstance(v2, np.ndarray):
                out[k+"."+k2] = [repr(float(x)) for x in np.asarray(v2, float).ravel()]
            elif isinstance(v2, (int, float)):
                out[k+"."+k2] = repr(float(v2))
    elif isinstance(v, (int, float)):
        out[k] = repr(float(v))
print("@@@" + json.dumps(out, sort_keys=True))
'''


def run_child(pdb, env_extra=None):
    env = dict(os.environ)
    env.update(env_extra or {})
    cp = subprocess.run([sys.executable, "-c", RUNNER % (_ROOT, pdb)],
                        capture_output=True, text=True, env=env, cwd=_ROOT, timeout=1800)
    for line in cp.stdout.splitlines():
        if line.startswith("@@@"):
            return json.loads(line[3:])
    raise RuntimeError((cp.stderr or cp.stdout)[-3000:])


def d1_cross_process(out, pdb="1CS9"):
    a = run_child(pdb)
    b = run_child(pdb)
    keys = sorted(set(a) & set(b))
    diff = [k for k in keys if a[k] != b[k]]
    out["D1_pdb"] = pdb
    out["D1_n_quantities"] = len(keys)
    out["D1_differing"] = diff
    out["D1_bit_identical_across_processes"] = not diff
    # and under a different ambient thread count, which must not change a number
    c = run_child(pdb, {"OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "4",
                        "MKL_NUM_THREADS": "4"})
    diff2 = [k for k in keys if a.get(k) != c.get(k)]
    out["D1_differing_under_4_threads"] = diff2
    out["D1_thread_count_invariant"] = not diff2


def d2_resume(out, manifest="smoke8"):
    """Run the manifest, delete half the cache, rerun, and compare the survivors."""
    import core
    from core.pipeline import PROD, cache_dir_for
    cdir = cache_dir_for(PROD, core.backend_report())
    files = sorted(f for f in os.listdir(cdir) if f.endswith(".json"))
    out["D2_cache_dir"] = cdir
    out["D2_n_cached"] = len(files)
    if len(files) < 2:
        out["D2_skipped"] = "not enough cached targets"
        return
    # The per-target record embeds wall-clock `timings`, `wall` and the `resumed` flag,
    # which are properties of the RUN and not of the result. Comparing raw bytes calls a
    # correct resume a failure, so compare the science and report the volatile keys.
    volatile = {"timings", "wall", "resumed"}

    def science_of(path):
        d = json.load(open(path))
        return {k: v for k, v in d.items() if k not in volatile}

    before = {}
    for f in files:
        before[f] = science_of(os.path.join(cdir, f))
    out["D2_volatile_keys_excluded"] = sorted(volatile)

    # simulate a kill: remove the last two entries, then resume the same manifest
    victims = files[-2:]
    bak = tempfile.mkdtemp()
    for v in victims:
        shutil.move(os.path.join(cdir, v), os.path.join(bak, v))
    cp = subprocess.run([sys.executable, "-m", "core.pipeline", "run",
                         "--manifest", manifest, "--workers", "1"],
                        capture_output=True, text=True, cwd=_ROOT, timeout=3600)
    out["D2_resume_rc"] = cp.returncode
    out["D2_resumed_targets"] = victims
    same, differ = [], []
    for v in victims:
        p = os.path.join(cdir, v)
        if not os.path.exists(p):
            differ.append(v + " (missing after resume)")
            continue
        (same if science_of(p) == before[v] else differ).append(v)
    # untouched entries must not have been rewritten either
    for f in files:
        if f in victims:
            continue
        if science_of(os.path.join(cdir, f)) != before[f]:
            differ.append(f + " (untouched entry changed)")
    out["D2_identical_after_resume"] = same
    out["D2_differing_after_resume"] = differ
    out["D2_resume_is_bit_identical"] = not differ
    shutil.rmtree(bak, ignore_errors=True)


def d3_key_types(out):
    from core import cache as C
    ns, ver = "verify_probe", 1
    ks = {
        "int_1": C.key(ns, ver, x=1),
        "float_1": C.key(ns, ver, x=1.0),
        "bool_True": C.key(ns, ver, x=True),
        "str_1": C.key(ns, ver, x="1"),
        "none": C.key(ns, ver, x=None),
        "list_12": C.key(ns, ver, x=[1, 2]),
        "tuple_12": C.key(ns, ver, x=(1, 2)),
        "float_eps": C.key(ns, ver, x=np.nextafter(1.0, 2.0)),
        "arr_12": C.key(ns, ver, x=np.array([1, 2])),
        "arr_12_f8": C.key(ns, ver, x=np.array([1.0, 2.0])),
        "ver2": C.key(ns, 2, x=1),
        "other_ns": C.key("verify_probe2", ver, x=1),
        "extra_param": C.key(ns, ver, x=1, y=0),
    }
    out["D3_keys"] = ks
    out["D3_n_distinct"] = len(set(ks.values()))
    out["D3_n_total"] = len(ks)
    out["D3_all_distinct"] = len(set(ks.values())) == len(ks)
    dupes = {}
    for name, k in ks.items():
        dupes.setdefault(k, []).append(name)
    out["D3_collisions"] = {k: v for k, v in dupes.items() if len(v) > 1}
    # unrepresentable types must be REFUSED, not str()-ed
    try:
        C.key(ns, ver, x=object())
        out["D3_refuses_unhashable_object"] = False
    except TypeError:
        out["D3_refuses_unhashable_object"] = True


def d4_forged_collision(out):
    """Write an entry, corrupt its recorded params, and demand a raise on read."""
    from core import cache as C
    ns, ver = "verify_collide", 1
    arrays = {"a": np.arange(4.0)}
    k = C.store(ns, ver, arrays, alpha=1, beta="x")
    npz, side = C._paths(ns, k)
    out["D4_stored_key"] = k
    out["D4_sidecar_exists"] = os.path.exists(side)
    C._MEM.clear()                       # force the disk path
    meta = json.load(open(side))
    meta["params"] = json.dumps(C._canon({"alpha": 2, "beta": "x"}),
                                sort_keys=True, separators=(",", ":"))
    json.dump(meta, open(side, "w"))
    C._MEM.clear()
    try:
        C.load(ns, ver, alpha=1, beta="x")
        out["D4_raises_CacheCollision"] = False
        out["D4_result"] = "SERVED A STALE ARRAY"
    except C.CacheCollision as exc:
        out["D4_raises_CacheCollision"] = True
        out["D4_message"] = str(exc)[:200]
    except Exception as exc:
        out["D4_raises_CacheCollision"] = False
        out["D4_result"] = f"{type(exc).__name__}: {exc}"
    C.clear(ns)


def d5_config_key(out):
    from dataclasses import replace, fields
    from core.pipeline import PROD
    import core
    base = PROD.key({"x": 1})
    moved, ignored = [], []
    for f in fields(PROD):
        cur = getattr(PROD, f.name)
        if isinstance(cur, bool):
            new = not cur
        elif isinstance(cur, int) and not isinstance(cur, bool):
            new = cur + 1
        elif isinstance(cur, float):
            new = cur + 1.0
        elif isinstance(cur, str):
            new = cur + "_x"
        else:
            continue
        k = replace(PROD, **{f.name: new}).key({"x": 1})
        (ignored if k == base else moved).append(f.name)
    out["D5_fields_that_move_the_key"] = sorted(moved)
    out["D5_fields_that_do_not"] = sorted(ignored)
    # the backend set must be in the key too
    out["D5_backends_move_the_key"] = PROD.key({"x": 1}) != PROD.key({"x": 2})
    out["D5_science_fields"] = sorted(PROD.science().keys())


def main():
    out = {}
    for fn in (d3_key_types, d4_forged_collision, d5_config_key,
               d1_cross_process, d2_resume):
        try:
            fn(out)
        except Exception as exc:
            out[fn.__name__ + "_ERROR"] = f"{type(exc).__name__}: {str(exc)[:600]}"
        print(f"  {fn.__name__} done", flush=True)
    print(json.dumps(out, indent=2, sort_keys=True, default=str))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "determinism_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=str)
    return out


if __name__ == "__main__":
    main()
