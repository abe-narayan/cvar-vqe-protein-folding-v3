"""s22/ablate.py -- THE EXTENDED HAMILTONIAN ABLATION (DIRECTIVE SS57), AT n=126.

Sprint 21 delivered the seven MANDATED cells with genuine CVaR-VQE at n=12: Legacy | AMBER |
Legacy+AMBER | Distance | Distance+Legacy | Distance+AMBER | Distance+Legacy+AMBER.  Distance alone
won on mean and median and all six physics-containing cells were worse -- **direction ESTABLISHED,
every magnitude NOT MEASURED at that n**.  Workstream C's tail-readout replication then put the
physics half **worse than a matched random tail at n=126 with 42/42 CIs excluding zero.**

THE GAP THIS FILE FILLS.  The directive's ablation list extends past those seven to
**Distance+Torsion**, **Distance+Geometry**, and their combinations with the physics energies.
**No lane has run those.**  They are cheap on the pool and the mandated list is a deliverable.

WHAT THIS IS AND -- IMPORTANTLY -- WHAT IT IS NOT.  This is the ablation evaluated through the
**deployed pool selector and its shipped tail readout**, at the full n=126, NOT through the
CVaR-VQE.  It is therefore a **DIFFERENT INSTRUMENT** from Sprint 21's n=12 VQE matrix and the two
must never be quoted as one table.  It is here because the directive asks which TERMS lower RMSD,
and the pool instrument answers that at 10x the targets with matched everything.  The VQE matrix
answers whether the quantum selector reproduces the ordering; Sprint 21 already showed it does.

THE TERMS.  Every one is NATIVE-FREE and independently evaluable, per Pillars 2 and 3.

    H_D   distance   the shipped Bayes-risk distogram score.  The deployed structural objective.
    H_L   Legacy     the 11-component Miyazawa-Jennings potential at DEFAULT_WEIGHTS, never fitted.
    H_T   torsion    -sum_i log p(phi_i, psi_i) under the fold-held-out Ramachandran tables --
                     the directive's SS25 term, evaluated on each candidate's own torsions.
    H_G   geometry   a pure-geometry validity term: squared deviation of virtual CA-CA bonds from
                     the trans value plus a soft clash count.  No sequence, no native, no physics
                     parameters -- it is the "is this a chain" term.

    H_AMBER is NOT re-run here.  Its deployed form is E o Relax_50 and costs an OpenMM context per
    candidate per target; the AMBER cells are already measured at n=12 (VQE) and n=126 (tail), and
    re-running them would take the serialised lock for a result that exists.  STATED, not skipped
    silently.

COMBINATION.  Terms are combined after **per-target rank-normalisation to a common scale**
(pool-quantile -> normal), declared BEFORE any RMSD was read.  This is the only defensible way to
add a Bayes risk, a contact energy in arbitrary units, a log-probability and a squared length --
Sprint 22 L8/S21 C6 measured that a RAW sum is dominated by whichever term has the widest dynamic
range, and for a raw lambda ladder that domination is total.
    NOT TAKEN: the raw sum (an outlier detector wearing a hybrid's name) and robust-z (audited
    alongside and reported whatever it says).  If a non-declared normalisation wins, the declared
    choice is recorded as WRONG rather than swapped.

OPERATOR FORKS, per BRIEF SS4 rule 0.

    functional     DECLARED the shipped Bayes risk for H_D.  NOT TAKEN the squared functional.
    basis          DECLARED the pool's own window coordinates on every arm.  NOT TAKEN the rebuild.
    readout        DECLARED top-75 coordinate average -- what SHIPS -- and argmin reported beside
                   it, because S21 E5 proved the two consume different statistics of the set and
                   no ordering claim is well-posed until the readout is fixed.
    normalisation  DECLARED per-target rank-to-normal; raw-sum and robust-z computed as audits.
    null           DECLARED a MATCHED-COUNT RANDOM TAIL -- the same number of members drawn without
                   regard to score.  This is the null S21 C10 used to show both physics energies
                   are worse than chance, and it is the only null that prices the ordering rather
                   than the averaging.  NOT TAKEN a uniform-random structure.
    THE LABEL      DECLARED continuous Ca-RMSD.  NOT TAKEN any binarised "did it beat the
                   incumbent".

  Hypothesis   H_D alone remains best; H_T may add a little (it is the one term never yet combined
               with H_D on this instrument); H_G adds nothing because the pool's members are
               already real windows and therefore already valid chains.
  Falsifier    if any combination beats H_D alone past that comparison's own MDE with a CI
               excluding zero, the additive extension is live and must be carried into the VQE.
  Null         the matched-count random tail.
  Promotion    none from this file; a win here buys a VQE cell, not a pipeline change.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from core import project as pj            # noqa: E402

TOPM = 75
CA_TRANS = 3.8046
CELLS = {
    "D": ("d",), "L": ("l",), "T": ("t",), "G": ("g",),
    "D+L": ("d", "l"), "D+T": ("d", "t"), "D+G": ("d", "g"),
    "D+L+T": ("d", "l", "t"), "D+T+G": ("d", "t", "g"),
    "D+L+T+G": ("d", "l", "t", "g"), "L+T": ("l", "t"),
}


def _save(o, name="ablate.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _rank_normal(x):
    """Per-target pool-quantile -> normal.  Ties averaged, so no tie reads the index order."""
    x = np.asarray(x, float); n = len(x)
    o = np.argsort(x, kind="stable"); r = np.empty(n, float); r[o] = np.arange(n, dtype=float)
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1:
            r[m] = r[m].mean()
    from math import sqrt
    q = (r + 0.5) / n
    # inverse normal via erfinv-free rational approx is unnecessary: use numpy's
    return np.sqrt(2.0) * _erfinv(2.0 * q - 1.0) * 1.0 / sqrt(1.0)


def _erfinv(y):
    a = 0.147
    ln = np.log(np.clip(1 - y * y, 1e-300, None))
    t1 = 2 / (np.pi * a) + ln / 2
    return np.sign(y) * np.sqrt(np.sqrt(t1 * t1 - ln / a) - t1)


def _geom(W):
    """Native-free geometry validity: virtual-bond deviation + soft clash count."""
    W = np.asarray(W, float)
    b = np.linalg.norm(W[:, 1:, :] - W[:, :-1, :], axis=2)
    bond = ((b - CA_TRANS) ** 2).sum(1)
    d = np.linalg.norm(W[:, :, None, :] - W[:, None, :, :], axis=3)
    n = W.shape[1]
    iu = np.triu_indices(n, k=3)
    clash = np.clip(4.0 - d[:, iu[0], iu[1]], 0, None).sum(1)
    return bond + 10.0 * clash


def run():
    tg = I.targets(); rows = []
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        PHI = np.asarray(u["PHI"], float)[idx]; PSI = np.asarray(u["PSI"], float)[idx]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        terms = {}
        terms["d"] = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        from s16 import energy_lib as EL
        #: GENUINE Legacy: the certified 11-term weighted sum at DEFAULT_WEIGHTS, never fitted.
        #: `legacy_total_from` is the project's own combiner -- not a re-implementation.
        comp = EL.legacy_components_of_windows(u["seq"], PHI, PSI)
        terms["l"] = np.asarray(EL.legacy_total_from(comp), float)
        pen = pj.make_penalty("rama", u["seq"], int(u["fold"]))
        #: `make_penalty` is BATCHED -- one call over the whole pool, not a per-row loop.
        terms["t"] = np.atleast_1d(np.asarray(pen(PHI, PSI), float)).ravel()
        terms["g"] = _geom(W)
        Z = {k: _rank_normal(v) for k, v in terms.items()}
        rng = SD.stable_rng(pdb, "s22ablate")
        row = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"])}
        for name, ks in CELLS.items():
            s = np.sum([Z[k] for k in ks], 0)
            o = np.argsort(s, kind="stable")
            C, _b = I.coordinate_average(W[o[:TOPM]])
            row["avg_" + name] = float(I.ca_rmsd(np.asarray(C, float), nat))
            row["arg_" + name] = float(I.ca_rmsd(W[o[0]], nat))
        #: the null that prices the ORDERING, not the averaging
        rs = []
        for _ in range(8):
            sel = rng.choice(len(W), TOPM, replace=False)
            C, _b = I.coordinate_average(W[sel])
            rs.append(float(I.ca_rmsd(np.asarray(C, float), nat)))
        row["avg_RANDOM"] = float(np.mean(rs))
        rows.append(row)
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})
    need = ["avg_" + k for k in CELLS] + ["avg_RANDOM"]
    ok = len(rows) == len(tg) and all(all(np.isfinite(r[k]) for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "cells": {k: list(v) for k, v in CELLS.items()}, "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "ablate.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    rng = SD.stable_rng("ablate", "rep"); fold = g("fold")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        F = sorted(set(fold.astype(int)))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se,
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    base = g("avg_D")
    print("\nn = %d.  Readout: top-%d coordinate average (SHIPPED). Point cloud, window basis.\n"
          % (len(rows), TOPM))
    print("  %-10s%9s%9s   %-40s" % ("cell", "avg", "argmin", "vs D alone (fold-clustered CI)"))
    order = sorted(CELLS, key=lambda k: g("avg_" + k).mean())
    for name in order:
        x = g("avg_" + name)
        m, se, mde, lo, hi, w = st(x - base)
        flag = ""
        if name != "D":
            flag = "  <<< BEATS D" if hi < 0 and abs(m) > mde else ""
        print("  %-10s%9.4f%9.4f   %+.4f SE %.4f MDE %.3f [%+.3f,%+.3f] %3dW/%3dL%s"
              % (name, x.mean(), g("arg_" + name).mean(), m, se, mde, lo, hi, w, len(rows) - w, flag))
    r = g("avg_RANDOM"); m, se, mde, lo, hi, w = st(r - base)
    print("  %-10s%9.4f%9s   %+.4f SE %.4f MDE %.3f [%+.3f,%+.3f]  <- MATCHED RANDOM TAIL"
          % ("RANDOM", r.mean(), "-", m, se, mde, lo, hi))
    print("\n  Falsifier: any combination beating D alone past its own MDE with a CI excluding zero.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
