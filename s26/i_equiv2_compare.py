#!/usr/bin/env python
"""s26/i_equiv2_compare.py -- compare the arms `verify/run_equiv2.sh` just ran, per target.

Reads each arm's per-target checkpoints straight from `bench_results/cache/<cfg_key>/`, computing
the keys the way `core.pipeline.run` does (the legacy arm's key carries `{"forced": "legacy"}`,
the consolidated arms' keys carry the parent's `core.backend_report()`), and compares the
emitted arrays with `==`: `ca` (the built chain), `fit_ca`, `phi`, `psi`, `avg_ca`, and the
AMBER trace when the arms ran stage 4.  Bit-identity is the expectation for baseline vs
`exact`; `fd` and `analytic` are expected to differ (see `core/project.py`'s docstring).

    python s26/i_equiv2_compare.py [--manifest smoke8] [--no-amber]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import replace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import core                                                              # noqa: E402
from core import pipeline as P                                           # noqa: E402

KEYS = ("ca", "fit_ca", "phi", "psi", "avg_ca")
SCALARS = ("rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full", "n_windows", "n_top")


def arm_cfg(grad, no_amber):
    return replace(P.PROD, threads_per_worker=1, amber=not no_amber, quantum=False, legacy=False,
                   dev_mode=("no-amber" if no_amber else P.PROD.dev_mode), project_grad=grad)


def load_arm(cfg, backends, manifest):
    key = cfg.key(backends)
    d = os.path.join(P.CACHE_ROOT, key)
    want = {p.pdb for p in P.manifest(manifest)}
    rows = {}
    if os.path.isdir(d):
        for f in sorted(os.listdir(d)):
            if f.endswith(".json") and f[:-5] in want:
                with open(os.path.join(d, f)) as fh:
                    r = json.load(fh)
                if r.get("cfg_key") == key:
                    rows[f[:-5]] = r
    return key, rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="smoke8")
    ap.add_argument("--no-amber", action="store_true")
    a = ap.parse_args(argv)
    os.environ.pop("CORE_BACKENDS", None)
    default_backends = core.backend_report()
    arms = {
        "baseline_legacy": load_arm(arm_cfg("exact", a.no_amber), {"forced": "legacy"}, a.manifest),
        "opt_exact": load_arm(arm_cfg("exact", a.no_amber), default_backends, a.manifest),
        "opt_fd": load_arm(arm_cfg("fd", a.no_amber), default_backends, a.manifest),
        "opt_analytic": load_arm(arm_cfg("analytic", a.no_amber), default_backends, a.manifest),
    }
    out = {"manifest": a.manifest, "no_amber": a.no_amber, "backends_default": default_backends,
           "arms": {k: {"cfg_key": v[0], "n_records": len(v[1])} for k, v in arms.items()},
           "pairs": {}}
    print("arm             cfg_key           n")
    for k, (key, rows) in arms.items():
        print(f"{k:15} {key}  {len(rows)}")
    keys = set(v[0] for v in arms.values())
    out["all_four_cfg_keys_distinct"] = len(keys) == 4
    print(f"all four cfg_keys distinct: {out['all_four_cfg_keys_distinct']}")

    base_rows = arms["opt_exact"][1]
    for k in ("baseline_legacy", "opt_fd", "opt_analytic"):
        rows = arms[k][1]
        shared = sorted(set(rows) & set(base_rows))
        rec = {"n_shared": len(shared), "per_key": {}, "scalars": {}}
        for kk in KEYS + ("amber_ca",):
            diffs, ident = [], 0
            for pdb in shared:
                x, y = base_rows[pdb].get(kk), rows[pdb].get(kk)
                if x is None or y is None:
                    continue
                x, y = np.asarray(x, float), np.asarray(y, float)
                if x.shape != y.shape:
                    diffs.append(float("inf"))
                    continue
                d = float(np.max(np.abs(x - y))) if x.size else 0.0
                diffs.append(d)
                ident += int(np.array_equal(x, y))
            if diffs:
                rec["per_key"][kk] = {"n": len(diffs), "n_bit_identical": ident,
                                      "max_abs_diff": max(diffs)}
        for sk in SCALARS:
            xs = [(base_rows[p].get(sk), rows[p].get(sk)) for p in shared]
            xs = [(x, y) for x, y in xs if x is not None and y is not None]
            if xs:
                rec["scalars"][sk] = {"n": len(xs),
                                      "max_abs_diff": max(abs(float(x) - float(y)) for x, y in xs),
                                      "n_equal": sum(int(x == y) for x, y in xs)}
        out["pairs"][f"opt_exact_vs_{k}"] = rec
        print(f"\nopt_exact vs {k}: {len(shared)} shared targets")
        for kk, v in rec["per_key"].items():
            print(f"  {kk:9} bit-identical {v['n_bit_identical']}/{v['n']}  max|d| {v['max_abs_diff']:.3e}")
        for sk, v in rec["scalars"].items():
            print(f"  {sk:9} equal {v['n_equal']}/{v['n']}  max|d| {v['max_abs_diff']:.3e}")
    path = os.path.join(ROOT, "s26", "results", "run_equiv2_compare.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
