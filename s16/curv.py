"""s16/curv.py -- the pre-registered follow-up to the quiet-control anomaly.

`s16/steer.py` found, at n = 8, that the arm which the framing REQUIRES to be inert --
a step along the strictly quiet half of the Jacobian spectrum -- moved RMSD by -0.106 A
with a confidence interval excluding zero, while its ORACLE counterpart did nothing
(-0.008). `s16/PREREG_steer.md` fixed the follow-up before the n = 126 run returned.
This module is that follow-up, unchanged from the pre-registration.

Two explanations, and one test each.

(a) CURVATURE. "Quiet" is a property of the Jacobian AT theta_fit. A finite step turns the
    Jacobian, so a direction quiet at the start need not be quiet along the path. The test
    is to shrink the step until the realised change matches the directional derivative. We
    take the derivative by central difference at h = 1e-4 rather than analytically, so no
    superposition convention can quietly enter: it is the same rmsd_of the ladder calls.

        ratio(s) = [RMSD(th + s d) - RMSD(th)] / [s * dRMSD/ds at 0]

    ratio -> 1 as s -> 0 is the linearisation working. If the derivative is POSITIVE (first
    order says the step should hurt) while the finite step HELPS, the gain is curvature and
    the loud/quiet language survives locally only -- it may not be quoted as a finite-step
    budget.

(b) A NON-STEERING SIDE EFFECT. The step toward theta_pool is a pull toward retrieval
    torsions. If a RANDOM direction of the same norm, projected into the same quiet
    subspace, gains the same amount, then the arm measures the STEP and not the DIRECTION,
    and there is no steering in it at all. Three draws, seeded through s15.seed.

A third control is added for completeness and is declared post hoc here rather than
discovered later: a random matched-norm direction in the FULL space, which prices "any
perturbation of this size" without reference to the spectrum.

ORACLE use: the native is read only to evaluate RMSD and to form the true error `e` for
labelled ORACLE arms. No native quantity enters a direction, a step, or a decision.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s16.steer import (pool_torsions, rmsd_of, subspaces,   # noqa: E402
                       _native_torsions)

H = 1e-4
FINE = (0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0)
N_RAND = 3


def _deriv(phi, psi, nat, vec, n):
    """dRMSD/ds at s = 0 along `vec`, by central difference."""
    th = np.concatenate([phi, psi])
    a = th + H * vec
    b = th - H * vec
    return (rmsd_of(a[:n], a[n:], nat) - rmsd_of(b[:n], b[n:], nat)) / (2 * H)


def run(targets=None, n_start=4):
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    data = C.gather(tg)

    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    rows = []
    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, n, fold = d["i"], d["j"], d["n"], d["fold"]
        sd = d["sd"]; nat = d["nat"]
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)

        starts = D.starts(p, d["seq"], n, fold, n_start, SD.stable_rng(p, "s16steer"))
        best = None
        for phi0, psi0, _tag in starts:
            ph, ps, f = A.fit(dhat, sd, i, j, phi0, psi0)
            if best is None or f < best[0]:
                best = (f, ph, ps)
        _f, phi, psi = best
        th = np.concatenate([phi, psi])
        m = 2 * n

        J, V, sv, _CA = subspaces(phi, psi)
        ppl, psl = pool_torsions(p, n, fold)
        u = A.wrap(np.concatenate([ppl, psl]) - th)
        Vq2 = V[:, m // 2:]
        q = Vq2 @ (Vq2.T @ u)                       # the arm under investigation
        nq = float(np.linalg.norm(q))

        nphi, npsi = _native_torsions(p, nat, d["seq"], fold)      # ORACLE, diagnostics
        e = A.wrap(np.concatenate([nphi, npsi]) - th)

        rng = SD.stable_rng(p, "s16curv")
        Vl4 = V[:, :4]
        #: POST HOC, added after the n = 126 flagship returned and labelled as such.  The
        #: pre-registration's fourth rule fired: rank-4 loud carries 0.875 of the RMSD to
        #: first order, yet stepping the TRUE error's loud component buys -0.028 [-0.106,
        #: +0.051] while the full true error buys -3.604.  The live question is therefore no
        #: longer the quiet control (which went inert at n = 126, -0.003 [-0.011, +0.007])
        #: but whether the loud/quiet split is a valid FINITE-STEP object at all.  The
        #: pre-registration permits a new arm only as a diagnostic of a failure; this is one.
        dirs = {"quiet_half": q, "full": u,
                "ORACLE_quiet_half": Vq2 @ (Vq2.T @ e),
                "ORACLE_loud4": Vl4 @ (Vl4.T @ e),
                "ORACLE_true": e}
        for k in range(N_RAND):
            #: (b) matched-norm RANDOM direction inside the SAME quiet subspace.
            z = rng.standard_normal(m)
            z = Vq2 @ (Vq2.T @ z)
            dirs[f"rand_quiet{k}"] = z * (nq / max(float(np.linalg.norm(z)), 1e-12))
            #: post hoc control: matched-norm random direction in the FULL space.
            y = rng.standard_normal(m)
            dirs[f"rand_full{k}"] = y * (nq / max(float(np.linalg.norm(y)), 1e-12))

        base = rmsd_of(phi, psi, nat)
        out = {}
        for name, vec in dirs.items():
            nv = float(np.linalg.norm(vec))
            if nv < 1e-12:
                out[name] = {"deriv": 0.0, "norm": 0.0,
                             "lad": {str(s): base for s in FINE}}
                continue
            g = _deriv(phi, psi, nat, vec, n)
            lad = {}
            for s in FINE:
                x = th + s * vec
                lad[str(s)] = rmsd_of(x[:n], x[n:], nat)
            out[name] = {"deriv": float(g), "norm": nv, "lad": lad}

        rows.append({"pdb": p, "n": int(n), "fold": int(fold), "base": base,
                     "quiet_norm": nq, "full_norm": float(np.linalg.norm(u)),
                     "arms": out})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)
            json.dump({"n": len(rows), "rows": rows},
                      open(os.path.join(RESULTS, "curv.json"), "w"))

    json.dump({"n": len(rows), "rows": rows},
              open(os.path.join(RESULTS, "curv.json"), "w"))
    report(rows)
    return rows


def _boot(x, y, rng, B=4000):
    """paired bootstrap over TARGETS, the unit of analysis everywhere in this project."""
    d = np.asarray(y) - np.asarray(x)
    k = len(d)
    idx = rng.integers(0, k, size=(B, k))
    m = d[idx].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows):
    rng = SD.stable_rng("curv", "report")
    base = np.array([r["base"] for r in rows])
    print(f"\nn = {len(rows)}   plain fit {base.mean():.3f} A\n")

    print("A. IS THE QUIET GAIN FIRST ORDER?  ratio -> 1 as s -> 0 is the linearisation working")
    print(f"  {'arm':<20}{'dRMSD/ds':>10}" + "".join(f"{s:>8}" for s in FINE))
    for name in ("quiet_half", "full", "ORACLE_quiet_half", "ORACLE_loud4", "ORACLE_true"):
        g = np.array([r["arms"][name]["deriv"] for r in rows])
        cells = []
        for s in FINE:
            v = np.array([r["arms"][name]["lad"][str(s)] for r in rows])
            pred = s * g
            ok = np.abs(pred) > 1e-9
            cells.append(np.mean((v - base)[ok] / pred[ok]) if ok.any() else float("nan"))
        print(f"  {name:<20}{g.mean():>10.4f}" + "".join(f"{c:>8.2f}" for c in cells))

    print("\n  (a POSITIVE derivative with a NEGATIVE finite-step effect is curvature, not"
          " steering.)")

    print("\nB. DOES THE DIRECTION MATTER, OR ONLY THE STEP?  RMSD at each step, vs plain fit")
    print(f"  {'arm':<20}" + "".join(f"{s:>8}" for s in FINE) + "     best step   effect [95% CI]")
    groups = {"quiet_half": ["quiet_half"], "full": ["full"],
              "rand_quiet": [f"rand_quiet{k}" for k in range(N_RAND)],
              "rand_full": [f"rand_full{k}" for k in range(N_RAND)],
              "ORACLE_quiet_half": ["ORACLE_quiet_half"],
              "ORACLE_loud4": ["ORACLE_loud4"], "ORACLE_true": ["ORACLE_true"]}
    for label, members in groups.items():
        curve = []
        for s in FINE:
            v = np.mean([[r["arms"][mn]["lad"][str(s)] for r in rows] for mn in members], axis=0)
            curve.append(v)
        mus = [c.mean() for c in curve]
        k = int(np.argmin(mus))
        mu, lo, hi = _boot(base, curve[k], rng)
        print(f"  {label:<20}" + "".join(f"{m:>8.3f}" for m in mus)
              + f"     s={FINE[k]:<5} {mu:+.3f} [{lo:+.3f},{hi:+.3f}]")

    print("\n  (if rand_quiet matches quiet_half, the arm measures the STEP and there is no"
          " steering in it.)")
    print("  CAVEAT: the step here is chosen IN SAMPLE at the ladder's argmin, so every"
          " effect above is\n  optimistic. It is the same rule for every arm, so the"
          " arm-vs-arm contrast is fair, but no\n  number in panel B may be quoted as an"
          " out-of-sample effect -- steer.py's leave-fold-out\n  ladder is the only place"
          " those live.")


if __name__ == "__main__":
    run()
