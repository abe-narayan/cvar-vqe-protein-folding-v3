"""Shared apparatus for the ERROR-DIRECTION agent (sprint 12).

The question: can the SIGN of the distogram's per-pair error be predicted, and is a
sign-corrected objective worth what the forensics agent's ORACLE arm says it is?

Everything here is a pure function of a per-pair TARGET DISTANCE VECTOR, pushed through
the production terminal path (score K=500 -> top-75 -> coordinate average -> project).
`rr` / `nat_ca` are used for EVALUATION and as leave-fold-out TRAINING LABELS only; every
arm that consumes them at inference time is named `o_*` / ORACLE.

benchmark60 is never read.  dev24 is never run.
"""
from __future__ import annotations
import os, sys, json, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")

from s12 import instrument as I
from s12 import obj_common as OC          # read-only reuse of the pool pack + emit path

CACHE = os.path.join(ROOT, "s12", "cache")
os.makedirs(CACHE, exist_ok=True)


# --------------------------------------------------------------------------- arms
def risk_from_target(grid, tgt):
    """The L1 risk table of a POINT estimate at `tgt` -- exactly what fail_pairs.py's
    `o_sign` arm substitutes for the shipped 17-bin Bayes risk."""
    return np.abs(np.asarray(grid, float)[None, :] - np.asarray(tgt, float)[:, None])


def score_risk(D, risk, grid):
    g = np.clip(((np.asarray(D, float) - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


def emit_target(d, tgt, M=75, lam=0.3):
    """Score the pool with a point-estimate L1 objective at `tgt` and run the real path."""
    sc = score_risk(d["D"], risk_from_target(d["grid"], tgt), d["grid"])
    return OC.emit(d, sc, M=M, lam=lam)


def emit_risk(d, risk, M=75, lam=0.3):
    return OC.emit(d, score_risk(d["D"], risk, d["grid"]), M=M, lam=lam)


def corrupt_sign(s_true, rate, rng):
    """Flip a `rate` fraction of the oracle signs (exactly, not in expectation)."""
    s = np.asarray(s_true, float).copy()
    k = int(round(rate * s.size))
    if k > 0:
        idx = rng.choice(s.size, k, replace=False)
        s[idx] *= -1.0
    return s


def sign_target(exp, s, delta):
    """The sign-corrected objective: shift each predicted distance `delta` in direction s."""
    return np.maximum(np.asarray(exp, float) + float(delta) * np.asarray(s, float), 1.5)


# --------------------------------------------------------------------------- bookkeeping
def controls():
    """The 18 length-matched control targets fail_contrast.py used for MATCH18."""
    p = os.path.join(ROOT, "s12", "results", "fail_contrast.json")
    return [c["ctrl"] for c in json.load(open(p))["controls"]]


def group_stats(rows, key, base_key="base", fail_key="fail18"):
    f = [r for r in rows if r[fail_key]]
    c = [r for r in rows if not r[fail_key]]
    out = {}
    for nm, g in (("fail18", f), ("match18", c)):
        if not g:
            continue
        a = np.array([r["arms"][key] for r in g], float)
        b = np.array([r["arms"][base_key] for r in g], float)
        st = I.paired(a, b)
        out[nm] = dict(mean=float(a.mean()), delta=st["mean_diff"], ci=st["ci95"],
                       wl=[st["n_better"], st["n_worse"]])
    return out
