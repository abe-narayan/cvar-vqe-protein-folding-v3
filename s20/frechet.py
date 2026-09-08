"""s20/frechet.py -- AVERAGE ON THE MANIFOLD INSTEAD OF AVERAGING THEN PROJECTING.

THE MEASURED LOSS THIS TARGETS.  The pipeline's readout is:

    superpose the top-75 on their medoid  ->  arithmetic mean in ambient R^3  ->  PROJECT onto
    the ideal-geometry manifold

and the projection costs **+0.166 A [+0.130, +0.202]**, 16W/110L, 5/5 folds (L2), paid on **every
target**.  It is the largest measured loss in the pipeline that is not the predictor, and it is
1.9x the instrument's MDE.

WHY IT IS PAID.  The arithmetic mean of superposed structures is a **shrinkage estimator**: it
contracts toward the centroid.  Measured (L20): the coordinate average's mean virtual Ca-Ca bond
is **2.961 A against a physical 3.804 A -- 22.2% contraction, global minimum bond 0.649 A.**  The
mean is not a peptide; the projection's job is to make it one again, and the +0.166 is what that
costs.

    THE QUESTION: if the average is computed ON the manifold in the first place, is the tax
    avoided rather than paid?

THE OPERATOR.  The manifold-constrained Frechet mean: the ideal-geometry structure minimising the
sum of squared superposition distances to the members,

    theta* = argmin_theta  sum_k  RMSD^2( build_ca(theta), W_k )

By construction its output has **exactly ideal bond geometry** -- contraction is impossible, not
merely discouraged -- so if the +0.166 is the price of un-contracting, this operator never incurs
it.

  THIS IS NOT TORSION-SPACE AVERAGING, and the distinction is the whole point.  The recorded
  result that "coordinate averaging beats torsion averaging by 1.024 A" is about averaging the
  ANGLES, which is a different and much worse operator: circular means of torsions ignore the
  chain's geometry entirely.  The Frechet mean minimises a **coordinate** objective and merely
  *constrains* the answer to the manifold.  Same target as the coordinate average, different
  feasible set.

ARMS.  All native-free; the native is read only to score.  Every arm states its BASIS (L20).

    avg          the coordinate average                       POINT CLOUD, ~3.048, contracted 22.2%
    proj_lam0    avg -> projection at lam=0                   BUILT, the incumbent ~3.204
    proj_lam03   avg -> projection at lam=0.3                 BUILT, the shipped arm ~3.215
    frechet      the manifold-constrained Frechet mean        BUILT, the hypothesis
    frechet_w    the same, weighted by the shipped scores     BUILT
    CONTROLS
    frechet_med  Frechet mean started from the MEDOID's torsions, not from the projection
    torsmean     circular mean of the members' torsions       BUILT, the recorded-bad operator

PRE-REGISTRATION, written before the run.

  Hypothesis:          averaging ON the manifold avoids the projection tax rather than paying it.
  Expected mechanism:  the +0.166 is the cost of un-contracting a shrinkage estimator; an operator
                       whose feasible set is the manifold cannot contract, so there is nothing to
                       undo.
  Primary endpoint:    frechet - proj_lam0, paired over targets, fold-clustered bootstrap CI.
  Falsifier:           frechet - proj_lam0 has a CI spanning zero, OR frechet is worse.  Either way
                       the tax is NOT the contraction and Sprint 18's "irreducible" stands.
  Null:                `torsmean` -- a manifold-respecting operator with no coordinate objective.
                       If frechet only matches it, the gain is "being on the manifold", not the
                       Frechet criterion.
  Matched control:     `frechet_med` -- same operator, different start.  Prices start-dependence,
                       since the projection is itself multi-start and could be start-limited.
  Budget:              n = 126, one L-BFGS per arm per target.
  Promotion criterion: beats proj_lam0 with a CI excluding zero AND >= the 0.084 A MDE AND beats
                       `torsmean`, with ideal bond geometry verified on the output.

  EXPECTED, honestly.  I expect a real but PARTIAL recovery.  The projection is already a
  least-squares fit to the contracted cloud, so it is solving a closely related problem; the
  Frechet mean differs in fitting the MEMBERS rather than their contracted mean.  The gap between
  those two is the part I expect to win, and it is not obviously the whole 0.166.  A result of
  3.10-3.15 would support the mechanism without closing the tax.

  THE TRAP I AM WATCHING FOR.  A Frechet mean that lands closer to the MEDOID than to the mean has
  quietly become medoid selection, which is a different (and already-measured) operator.  Distance
  to the medoid is recorded per target for exactly this reason.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from core import project as pj              # noqa: E402
from s12 import instrument as I             # noqa: E402
from s14 import avgspace as AV              # noqa: E402
from s15 import seed as SD                  # noqa: E402

IDEAL_CA = 3.804


def _superpose_to(X, Y):
    """Kabsch-superpose Y onto X (both (n,3)), returning the rotated/translated Y."""
    Xc = X - X.mean(0)
    Yc = Y - Y.mean(0)
    U, S, Vt = np.linalg.svd(Yc.T @ Xc)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    return Yc @ (U @ D @ Vt) + X.mean(0)


def frechet(W, phi0, psi0, w=None, maxiter=300):
    """Manifold-constrained Frechet mean of the member set `W` (m, n, 3).

    Minimises  sum_k w_k * ||superpose(build_ca(theta), W_k) - build_ca(theta)||^2  over torsions.
    Implemented as alternating superposition (members onto the current chain) and an L-BFGS step
    on the torsions against the current superposed target -- the standard generalised-Procrustes
    scheme, with the crucial difference that the estimate is confined to the manifold throughout.
    """
    from scipy.optimize import minimize

    W = np.asarray(W, float)
    m, n, _ = W.shape
    w = np.ones(m) if w is None else np.asarray(w, float)
    w = w / max(w.sum(), 1e-12)

    def chain(x):
        return np.asarray(pj.build_ca_exact(x[:n][None], x[n:][None]), float)[0]

    x = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    prev = None
    for _sweep in range(6):
        C = chain(x)
        #: superpose every member onto the CURRENT estimate, then form the weighted target
        T = np.stack([_superpose_to(C, W[k]) for k in range(m)])
        target = (T * w[:, None, None]).sum(0)

        def fg(v):
            G = pj.frames(v[:n][None], v[n:][None])[0]
            CA = np.asarray(pj.build_ca_exact(v[:n][None], v[n:][None]), float)[0]
            d = CA - target
            return float((d * d).sum()), pj._torsion_grad(G, CA, 2.0 * d)

        res = minimize(fg, x, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-10})
        x = res.x
        if prev is not None and abs(prev - res.fun) < 1e-9:
            break
        prev = res.fun
    return x[:n], x[n:]


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows, t0 = [], time.time()
    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)
        W = np.asarray(AV.top75_windows(pdb)[0], float)
        P = I.pairwise_rmsd(W)
        b = I.medoid(P)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

        def emit(e, key, ca):
            ca = np.asarray(ca, float)
            e[key] = float(I.ca_rmsd(ca, nat))
            d = np.linalg.norm(ca[1:] - ca[:-1], axis=1)
            e["bond_" + key] = float(d.mean())

        e = {"n_mem": int(len(W))}
        emit(e, "avg", avg)
        emit(e, "proj_lam0", pr["fit_ca"])
        emit(e, "proj_lam03", pr["ca"])

        ph, ps = frechet(W, phi0, psi0)
        Cf = I.build_ca(ph, ps)
        emit(e, "frechet", Cf)
        #: THE TRAP: has it become medoid selection?
        e["frechet_to_medoid"] = float(I.ca_rmsd(Cf, W[b]))
        e["proj_to_medoid"] = float(I.ca_rmsd(np.asarray(pr["fit_ca"], float), W[b]))

        #: CONTROL -- same operator, started from the medoid's own torsions
        prm = I.project(np.asarray(W[b], float), seq, fold)
        ph2, ps2 = frechet(W, np.asarray(prm["phi"], float), np.asarray(prm["psi"], float))
        emit(e, "frechet_med", I.build_ca(ph2, ps2))

        #: weighted by the shipped scores (native-free)
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        s = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        wt = np.exp(-(s - s.min()) / max(s.std(), 1e-9))
        ph3, ps3 = frechet(W, phi0, psi0, w=wt)
        emit(e, "frechet_w", I.build_ca(ph3, ps3))

        #: NULL -- circular mean of the members' torsions (the recorded-bad operator)
        PHI, PSI, _sim = __import__("s14.retprior", fromlist=["x"]).windows(pdb, "top75")
        tp = np.arctan2(np.sin(PHI).mean(0), np.cos(PHI).mean(0))
        ts = np.arctan2(np.sin(PSI).mean(0), np.cos(PSI).mean(0))
        emit(e, "torsmean", I.build_ca(tp, ts))

        rows.append({"pdb": pdb, "n": n, "fold": fold, **e})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": False},
                      open(os.path.join(RESULTS, "frechet.json"), "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(tg)},
              open(os.path.join(RESULTS, "frechet.json"), "w"))
    report(rows)
    return rows


def _boot(d, rng, B=4000):
    d = np.asarray(d, float); k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "frechet.json")))["rows"]
    rng = SD.stable_rng("frechet", "report")
    g = lambda k: np.array([r[k] for r in rows])          # noqa: E731
    base = g("proj_lam0")
    print(f"\nn = {len(rows)}.  Baseline = proj_lam0, the INCUMBENT built structure "
          f"({base.mean():.4f} A).  All arms NATIVE-FREE.\n")
    print(f"  {'arm':<14}{'basis':<12}{'RMSD':>8}{'bond':>8}{'vs incumbent [95% CI]':>28}{'W/L':>9}")
    for k, basis in (("avg", "POINT CLOUD"), ("proj_lam0", "built"), ("proj_lam03", "built"),
                     ("frechet", "built"), ("frechet_w", "built"),
                     ("frechet_med", "built"), ("torsmean", "built")):
        v = g(k); bd = g("bond_" + k)
        mu, lo, hi = _boot(v - base, rng)
        w = int((v < base).sum()); l = int((v > base).sum())
        tag = "" if k not in ("frechet_med", "torsmean") else "  CONTROL"
        print(f"  {k:<14}{basis:<12}{v.mean():>8.4f}{bd.mean():>8.3f}   "
              f"{mu:+.4f} [{lo:+.4f},{hi:+.4f}]{w:>5}/{l}{tag}")
    print(f"\n  ideal Ca-Ca bond {IDEAL_CA:.3f} A -- the Frechet arms must sit exactly there")
    print(f"  THE TRAP: dist(frechet, medoid) {g('frechet_to_medoid').mean():.3f} vs "
          f"dist(proj, medoid) {g('proj_to_medoid').mean():.3f}")
    print("    if the Frechet arm is much closer to the medoid than the projection is, it has")
    print("    become medoid selection -- a different, already-measured operator.")
    mu, lo, hi = _boot(g("frechet") - g("torsmean"), rng)
    print(f"\n  vs the NULL (torsion mean): {mu:+.4f} [{lo:+.4f},{hi:+.4f}]")
    print("\nREAD.  The projection tax is +0.166 [+0.130,+0.202].  If `frechet` recovers a")
    print("significant part of it the tax is the CONTRACTION and is avoidable; if not, Sprint 18's")
    print("'irreducible' stands and the tax is the price of the manifold itself.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
