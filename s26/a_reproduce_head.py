"""s26/a_reproduce_head.py -- lane A (Adversary): re-run lane E's s26/e_reproduce.py on HEAD.

The Examiner's script hard-codes its output path (s26/results/e_reproduce.json, lane E's
committed artefact).  This wrapper executes e_reproduce.main() unchanged and redirects its
single save_atomic call to s26/results/a_reproduce_head.json, so lane E's file is not
overwritten (lane contract section 0: append-or-create by the lane that owns it).  The
computation is the Examiner's code, byte for byte; only the destination path differs.

    python s26/jobrun.py --agent A --tag CPU --name a_reproduce_head --est-ram 0.1 -- \
        python s26/a_reproduce_head.py
"""
from __future__ import annotations

import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST  # noqa: E402

OUT = os.path.join(ROOT, "s26", "results", "a_reproduce_head.json")
_orig_save = ST.save_atomic


def _redirect(path, obj, *args, **kwargs):
    print(f"[a_reproduce_head] write redirected from {os.path.relpath(path, ROOT)} "
          f"to {os.path.relpath(OUT, ROOT)}")
    return _orig_save(OUT, obj, *args, **kwargs)


ST.save_atomic = _redirect

_spec = importlib.util.spec_from_file_location(
    "e_reproduce", os.path.join(ROOT, "s26", "e_reproduce.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
raise SystemExit(_mod.main())
