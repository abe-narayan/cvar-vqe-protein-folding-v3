#!/usr/bin/env python
"""s26/i_run_equiv2_arm.py -- one consolidated arm of `verify/run_equiv2.sh`.

`verify/run_equiv2.sh` was written when `core/project.py` read its gradient mode from the
environment variable `PROJECT_GRAD`.  That was the cache-key collision `verify/
grad_key_collision.py` documents, and the fix moved the mode into
`core.pipeline.Config.project_grad`, hashed into the cache key.  The `core.pipeline` CLI has no
flag for it, so this driver builds exactly the `Config` that `python -m core.pipeline run`
builds (same `replace(PROD, ...)` call) and sets `project_grad` on it.  It prints the same JSON
the CLI prints, so the script's `grep '"cfg_key"'` still works, plus one diagnostic line: whether
`openmm` was imported at all (with `--no-amber` it should not be).

    python s26/i_run_equiv2_arm.py --grad exact|analytic|fd [--manifest smoke8] [--no-amber]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import replace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from core import pipeline as P                                           # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--grad", required=True, choices=("exact", "analytic", "fd"))
    ap.add_argument("--manifest", default="smoke8", choices=P.MANIFESTS)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--no-amber", action="store_true")
    ap.add_argument("--backends", default=None)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)
    if a.manifest == "benchmark60":
        print("REFUSED: the benchmark is sealed (rule 1)")
        return 2
    # identical to core.pipeline.main's construction, plus the gradient mode
    cfg = replace(P.PROD, threads_per_worker=a.threads, amber=not a.no_amber,
                  quantum=False, legacy=False,
                  dev_mode=("no-amber" if a.no_amber else P.PROD.dev_mode),
                  project_grad=a.grad)
    out = P.run(a.manifest, cfg, workers=a.workers, backends_env=a.backends, limit=a.limit)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=2, default=float))
    print(json.dumps({"diagnostic": "openmm_imported", "value": "openmm" in sys.modules,
                      "project_grad": a.grad, "amber": cfg.amber}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
