"""s17/refine.py -- THE ONE MEASUREMENT THAT DECIDES THE FORWARD PLAN.

Sprint 17 closed selection at five levels and showed the shortlist binds: **at B <= 25 the
2.5 A target is unreachable at rho = 1**.  Every ceiling in that analysis -- 2.104 A in the
top-75, 1.711 A at K = 500, 1.313 A over the full universe -- is a ceiling **of the retrieval
pool**.

REFINEMENT ESCAPES THE POOL.  If the output is allowed to be a structure that is not a
retrieved window, none of those ceilings applies.  And refinement is the sprint's one measured
positive lever: the QUANTUM workstream found **local refinement worth ~0.7 A** (3.636 -> 2.95 A
over 6-residue windows), certified by exhaustive enumeration and reached by classical greedy at
1,024 evaluations.

BUT THERE IS A DIRECT TENSION IN THE RECORD, and this module exists to resolve it.

  * On the ENUMERATED instrument (19 targets, coarse k = 4 lattice, 6-residue windows), local
    refinement toward the distance objective buys 0.7 A.
  * On the REAL instrument, a multi-start torsion-space fit to the SAME objective returns
    **3.647 A** -- *worse* than the coordinate average it would replace (3.048 A).

Both cannot be the whole story.  Either the objective's optimum on the real instrument is worse
than the coordinate average (in which case refinement is harmful and the forward plan must not
rest on it), or local refinement from a GOOD start lands somewhere the multi-start global fit
never reaches (in which case it is the plan's foundation).

    THE QUESTION: starting from the best structure the pipeline actually builds, does moving
    toward the distance objective make it better or worse?

ARMS.  All native-free; the native is read only to score.

    avg          the coordinate average of the shipped top-75            the start, ~3.048 A
    proj         its ideal-geometry projection                           the torsion start point
    refine_full  local fit of ALL torsions from `proj`, single start      the global-ish move
    refine_win{w} local fit of one w-residue window at a time, swept      the QUANTUM-shaped move
    CONTROLS
    rand_obj     the identical refinement against a SHUFFLED distogram    zero-information
    rand_move    a random torsion perturbation of the SAME magnitude      matched-magnitude

PRE-REGISTRATION, written before the run.

  HYPOTHESIS.  Local refinement from the coordinate average improves on it, because the
  enumerated instrument says local moves toward this objective are worth ~0.7 A.

  EXPECTED.  Honestly: I expect this to FAIL on the full-torsion arm, because the multi-start
  fit to the same objective returns 3.647 A and there is no reason a local optimum reached from
  a better start should be better than the objective's own global optimum unless the objective
  is badly multi-modal.  The windowed arms are the live ones -- a windowed move changes few
  coordinates and may stop before leaving the good basin.

  SUCCESS.  Some refinement arm beats `avg` at n = 126, target as the unit, paired
  fold-clustered interval excluding zero, AND beats `rand_move` at matched displacement.

  FALSIFIER.  If every refinement arm is worse than `avg`, then **the distance objective's
  optimum is worse than the structure the pipeline already builds**, refinement toward it is
  harmful, and the forward plan cannot rest on it -- the objective must be replaced first.

  THE TRAP.  `rand_obj` is the control that decides whether any gain is the OBJECTIVE or merely
  the RELAXATION.  A torsion fit is also a smoothing operator; if refining against a shuffled
  distogram gains as much, the effect is not the objective's information.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

WINDOWS = (4, 6, 8)
SWEEPS = 3


def _obj(dhat, sd, i, j, phi, psi):
    """the native-free distance objective actually being minimised."""
    CA = I.build_ca(phi, psi)
    d = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    return float((((d - dhat) / sd) ** 2).sum())


def _windowed(dhat, sd, i, j, phi, psi, w, sweeps=SWEEPS):
    """refine one w-residue window at a time, holding the rest fixed; several sweeps.

    This is the QUANTUM workstream's move shape on the real instrument: few coordinates change
    per step, so the structure cannot leave its basin the way a full fit can.
    """
    n = len(phi)
    ph, ps = phi.copy(), psi.copy()
    for _ in range(sweeps):
        for a in range(0, max(n - w + 1, 1)):
            b = min(a + w, n)
            mask = np.zeros(2 * n, bool)
            mask[a:b] = True
            mask[n + a:n + b] = True
            #: `basis` restricts the fit to the window's coordinates
            B = np.eye(2 * n)[:, mask]
            p2, q2, _f = A.fit(dhat, sd, i, j, ph, ps, basis=B)
            ph, ps = p2, q2
    return ph, ps


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    t0 = time.time()
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        rng = SD.stable_rng(pdb, "s17refine")

        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

        e = {"avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
             "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
             "obj_proj": _obj(dhat, sd, i, j, phi0, psi0)}

        #: FULL refinement -- every torsion, single start from the projection
        pf, qf, _ = A.fit(dhat, sd, i, j, phi0, psi0)
        e["refine_full"] = float(I.ca_rmsd(I.build_ca(pf, qf), nat))
        e["obj_full"] = _obj(dhat, sd, i, j, pf, qf)
        e["disp_full"] = float(np.linalg.norm(A.wrap(np.concatenate([pf, qf])
                                                     - np.concatenate([phi0, psi0]))))

        #: WINDOWED refinement -- the shape the enumerated instrument found worth 0.7 A
        for w in WINDOWS:
            pw, qw = _windowed(dhat, sd, i, j, phi0, psi0, w)
            e[f"refine_win{w}"] = float(I.ca_rmsd(I.build_ca(pw, qw), nat))
            e[f"obj_win{w}"] = _obj(dhat, sd, i, j, pw, qw)
            e[f"disp_win{w}"] = float(np.linalg.norm(
                A.wrap(np.concatenate([pw, qw]) - np.concatenate([phi0, psi0]))))

        #: CONTROL 1 -- the identical refinement against a SHUFFLED distogram.  Decides whether
        #: any gain is the OBJECTIVE's information or merely the smoothing a torsion fit applies.
        perm = rng.permutation(len(dhat))
        pr_, qr_, _ = A.fit(dhat[perm], sd[perm], i, j, phi0, psi0)
        e["rand_obj"] = float(I.ca_rmsd(I.build_ca(pr_, qr_), nat))

        #: CONTROL 2 -- a random torsion move of the SAME magnitude as the full refinement
        z = rng.standard_normal(2 * n)
        z *= e["disp_full"] / max(float(np.linalg.norm(z)), 1e-12)
        th = np.concatenate([phi0, psi0]) + z
        e["rand_move"] = float(I.ca_rmsd(I.build_ca(th[:n], th[n:]), nat))

        rows.append({"pdb": pdb, "n": n, "fold": fold, **e})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False}, open(os.path.join(RESULTS, "refine.json"), "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(tg)},
              open(os.path.join(RESULTS, "refine.json"), "w"))
    report(rows)
    return rows


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "refine.json")))["rows"]
    rng = SD.stable_rng("refine", "report")
    g = lambda k: np.array([r[k] for r in rows])          # noqa: E731
    avg = g("avg")
    print(f"\nn = {len(rows)}.  Start = the coordinate average of the shipped top-75.")
    print("Every arm is NATIVE-FREE; the native is read only to score.\n")
    print(f"  {'arm':<16}{'RMSD':>8}{'median':>9}{'vs avg [95% CI]':>26}{'W/L':>9}"
          f"{'objective':>12}{'disp(rad)':>11}")
    order = ["proj", "refine_full"] + [f"refine_win{w}" for w in WINDOWS] + ["rand_obj", "rand_move"]
    print(f"  {'avg (start)':<16}{avg.mean():>8.3f}{np.median(avg):>9.3f}"
          f"{'--':>26}{'--':>9}{'--':>12}{'--':>11}")
    for k in order:
        if k not in rows[0]:
            continue
        v = g(k)
        mu, lo, hi = _boot(v - avg, rng)
        w = int((v < avg).sum()); l = int((v > avg).sum())
        ok = f"{g('obj_' + k.split('refine_')[-1]).mean():>12.1f}" if ("obj_" + k.split("refine_")[-1]) in rows[0] else f"{'--':>12}"
        dk = "disp_" + k.split("refine_")[-1]
        dd = f"{g(dk).mean():>11.3f}" if dk in rows[0] else f"{'--':>11}"
        tag = "  CONTROL" if k.startswith("rand") else ""
        print(f"  {k:<16}{v.mean():>8.3f}{np.median(v):>9.3f}   "
              f"{mu:+.3f} [{lo:+.3f},{hi:+.3f}]{w:>5}/{l}{ok}{dd}{tag}")

    print(f"\n  objective at the projection start: {g('obj_proj').mean():.1f}")
    print("\nREAD.  If every refinement arm is WORSE than the start, the distance objective's")
    print("optimum is worse than the structure the pipeline already builds -- refinement toward")
    print("it is harmful, and any forward plan resting on refinement must replace the objective")
    print("first.  `rand_obj` decides whether a gain is the objective's INFORMATION or merely the")
    print("SMOOTHING a torsion fit applies; `rand_move` prices the displacement itself.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
