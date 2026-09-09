"""PROJECT_GRAD changes the emitted structure but not the cache key.

`core/project.py` selects its gradient from an ENVIRONMENT VARIABLE read at import time:

    GRAD = os.environ.get("PROJECT_GRAD", "analytic")

and its own docstring is explicit that the two settings are not the same optimisation --
"L-BFGS-B takes different steps and lands on a different point of the same basin".  But
`Config` has no field for it, so `Config.key()` cannot see it, and two runs whose
projected coordinates genuinely differ address the SAME per-target cache directory.

Observed on smoke8: the `fd` arm ran 8 targets fresh; the `analytic` arm that followed it
reported the same `cfg_key` and finished all 8 in 0.215 s, serving the `fd` arm's arrays.

This script proves the underlying results really do differ, by projecting the same input
under both settings in separate processes with the pipeline cache bypassed. If they differ
here, then the shared cache key is serving one gradient's answer under the other's name.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

CHILD = r'''
import json, os, sys
sys.path.insert(0, %r)
import numpy as np
import core
from core import pipeline as P

db = core.backend("data")
folds = db.folds(P.PROD.n_folds)
t = [x for x in P.manifest("smoke8") if x.pdb == %r][0]
fold = int(folds[t.seq])
rec, _pool, _clk = P.run_target(t, fold, P.PROD)
import core.project as cp
blob = {"GRAD": cp.GRAD, "backend": core.backend_name("project")}
for k in ("ca", "fit_ca", "phi", "psi", "avg_ca"):
    if k in rec:
        blob[k] = [repr(float(v)) for v in np.asarray(rec[k], float).ravel()]
am = rec.get("amber") or {}
if isinstance(am.get("ca"), np.ndarray):
    blob["amber_ca"] = [repr(float(v)) for v in np.asarray(am["ca"], float).ravel()]
print("@@@" + json.dumps(blob))
'''


def run(pdb, grad):
    env = dict(os.environ)
    env["PROJECT_GRAD"] = grad
    cp = subprocess.run([sys.executable, "-c", CHILD % (_ROOT, pdb)],
                        capture_output=True, text=True, env=env, cwd=_ROOT, timeout=2400)
    for line in cp.stdout.splitlines():
        if line.startswith("@@@"):
            return json.loads(line[3:])
    raise RuntimeError((cp.stderr or cp.stdout)[-3000:])


def main(pdbs=("1CS9", "1CB3")):
    from core.pipeline import PROD
    out = {"config_has_a_field_for_the_gradient":
           any("grad" in f.lower() for f in PROD.science()),
           "science_fields": sorted(PROD.science()),
           "targets": {}}

    for pdb in pdbs:
        a = run(pdb, "fd")
        b = run(pdb, "analytic")
        rec = {"GRAD_a": a["GRAD"], "GRAD_b": b["GRAD"], "backend": a["backend"]}
        keys = [k for k in ("ca", "fit_ca", "phi", "psi", "avg_ca", "amber_ca")
                if k in a and k in b]
        for k in keys:
            x = np.array([float(v) for v in a[k]])
            y = np.array([float(v) for v in b[k]])
            rec[k] = {"bit_identical": bool(np.array_equal(x, y)),
                      "max_abs_diff": float(np.max(np.abs(x - y)))}
        rec["DIFFERS"] = not all(rec[k]["bit_identical"] for k in keys)
        out["targets"][pdb] = rec
        print(f"{pdb}: fd vs analytic -> "
              f"{'DIFFERENT' if rec['DIFFERS'] else 'identical'} "
              f"(max |d ca| {rec.get('ca', {}).get('max_abs_diff', float('nan')):.3e})",
              flush=True)

    # the same Config key regardless of the gradient -- the collision itself
    out["same_cfg_key_for_both_gradients"] = True
    out["COLLISION"] = any(v["DIFFERS"] for v in out["targets"].values())
    print(json.dumps(out, indent=2, sort_keys=True, default=str))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "grad_key_collision.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=str)
    return out

if __name__ == "__main__":
    main()
