"""SPRINT 15, coordinator -- DISTANCE GEOMETRY: solve for the structure, do not score it.

WHY THIS IS THE RIGHT FIRST EXPERIMENT.

Fourteen sprints have treated peptide prediction as generate-then-select. Sprint 14 measured
that frame to death: selection is discrimination-limited, the selection gap is ~1.9 A and
survives a CERTIFIED global optimum, every objective is at chance inside its own tail, and
running a variational optimiser is worse than best-of-N from an untrained circuit. Meanwhile
GENERATION is fine -- the retrieval pool holds a 1.711 A member and the k=4 torsion space
holds a 1.982 A answer from a random start.

**So the leverage is not in scoring candidates better. It is in not needing to score them.**

The shipped distogram already predicts every CA-CA pair distance with a per-pair uncertainty
`sd`, and **that uncertainty has never been used**. The pipeline consumes the distogram as a
Bayes-RISK SCORE -- a scalar that ranks candidates. That is a selection object. The natural
alternative is the one NMR structure determination has used for forty years: treat the
predicted distances as restraints and **solve** for the conformation that satisfies them.

    f(phi, psi) = sum_{|i-j|>=2}  w_ij * ( d_ij(phi,psi) - dhat_ij )^2 ,   w_ij = 1/sigma_ij^2

This is a different mathematical object from the Bayes-risk score, it is a GENERATION method
rather than a SELECTION method, it uses information (the per-pair sd) that the project has
never touched, and it is solved rather than sampled -- so the tail-ordering pathology that
defeats every search arm does not arise.

WHAT THIS MEASURES, and why both numbers matter.

  ORACLE arm   fit to the TRUE distance matrix. This is the ceiling of the REPRESENTATION:
               how well can ideal-geometry torsions reproduce a distance matrix at all?
               It is an ORACLE DIAGNOSTIC and can never be a headline.

  PREDICTED    fit to the distogram's expected distances, inverse-variance weighted.
               Fully native-free (the distogram is leave-fold-out). This is a real
               predictive arm.

The gap between them is exactly the distogram's error contribution, isolated from every
other stage. If the ORACLE arm is poor, torsion-space distance geometry is closed and we
learn that cheaply. If the ORACLE arm is good and the PREDICTED arm is poor, the distogram
is the bottleneck and its accuracy becomes the target. Either answer redirects the program.

GRADIENTS. `core.project` already carries an exact analytic torsion gradient:
`frames(phi,psi)` returns the cumulative frames, and `_torsion_grad(G, CA, g)` maps
d(objective)/d(CA) to d(objective)/d(phi,psi) in O(n) using the fact that a torsion rotates
its downstream suffix rigidly. So the fit costs one forward build and one reverse sweep per
iteration, not 2n forward builds.

For our objective, with r_ij = CA_i - CA_j and d_ij = |r_ij|:

    df/dCA_k = sum_{j != k} 2 * w_kj * (d_kj - dhat_kj) * r_kj / d_kj

Run:
    python -m s15.distgeo
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from scipy.optimize import minimize          # noqa: E402

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD             # noqa: E402
from core import project as pj               # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)


# ---------------------------------------------------------------- the objective and gradient
def _pair_terms(CA, i, j, dhat, w):
    """Weighted squared distance residual and its gradient w.r.t. CA. O(n^2) in pairs."""
    r = CA[i] - CA[j]                                    # (P, 3)
    d = np.sqrt((r * r).sum(1))
    d = np.maximum(d, 1e-9)
    resid = d - dhat
    f = float((w * resid * resid).sum())
    # dF/dCA_i = 2 w resid * r/d ;  dF/dCA_j = -that
    coef = (2.0 * w * resid / d)[:, None] * r            # (P, 3)
    g = np.zeros_like(CA)
    np.add.at(g, i, coef)
    np.add.at(g, j, -coef)
    return f, g


def fit_distances(dhat, w, i, j, phi0, psi0, maxiter=400):
    """L-BFGS fit of continuous (phi, psi) to a weighted distance target.

    Returns (phi, psi, f_final, n_iter).  Uses the exact analytic torsion gradient.
    """
    n = len(phi0)

    def fg(x):
        phi = x[:n]; psi = x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        f, gCA = _pair_terms(CA, i, j, dhat, w)
        gt = pj._torsion_grad(G, CA, gCA)
        return f, gt

    x0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])
    res = minimize(fg, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return res.x[:n], res.x[n:], float(res.fun), int(res.nit)


# ------------------------------------------------------------------------------ start points
def starts(pdb, seq, n, fold, n_start, rng):
    """Native-free initialisations, in the order they are tried.

    The retrieval-conditioned circular mean is Sprint 14's best torsion channel
    (phi 33.6, psi 59.2 deg) and is the natural warm start; the constant alpha-helix is the
    project's mandated zero-information control; the rest are prior draws.
    """
    from s14 import retprior as R
    out = []
    PHI, PSI, _ = R.windows(pdb, "top75")
    out.append((R.circ_mean(PHI, axis=0), R.circ_mean(PSI, axis=0), "retrieval_circmean"))
    out.append((np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)), "helix"))
    for s in range(max(0, n_start - 2)):
        k = rng.integers(0, PHI.shape[0])
        out.append((PHI[k].copy(), PSI[k].copy(), f"pool_draw_{s}"))
    return out[:n_start]


# ----------------------------------------------------------------------------------- the arm
def run_target(t, n_start=6, seed=0, maxiter=400):
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)
    i, j = I.pair_index(n)                                # |i-j| >= 2
    rng = SD.stable_rng(pdb, seed)

    dg = I.distogram(pdb, seq, fold)                      # leave-fold-out, native-free
    dhat = np.asarray(dg["expected"], float)
    sd = np.asarray(dg["sd"], float)
    sd = np.maximum(sd, 1e-3)

    # ORACLE target: the true distance matrix (diagnostic ceiling only)
    dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))

    arms = {
        "ORACLE_true_distances": (dtrue, np.ones_like(dtrue)),
        "pred_invvar_weighted": (dhat, 1.0 / (sd * sd)),
        "pred_unweighted": (dhat, np.ones_like(dhat)),
    }

    S = starts(pdb, seq, n, fold, n_start, rng)
    row = {"pdb": pdb, "n": n, "fold": fold, "arms": {}}
    for name, (target, w) in arms.items():
        best = None
        for phi0, psi0, tag in S:
            phi, psi, f, it = fit_distances(target, w, i, j, phi0, psi0, maxiter)
            ca = I.build_ca(phi, psi)
            r = float(I.ca_rmsd(ca, nat))
            if best is None or f < best["f"]:
                best = {"f": f, "rmsd": r, "start": tag, "iters": it}
            # the multi-start selection uses the OBJECTIVE only -- never the RMSD
        row["arms"][name] = best
    # the start's own quality, for the "did the fit do anything?" control
    phi0, psi0, _ = S[0]
    row["start_rmsd_retrieval_circmean"] = float(I.ca_rmsd(I.build_ca(phi0, psi0), nat))
    return row


def run(targets=None, n_start=6, out_name="distgeo"):
    tg = targets if targets is not None else I.targets()
    rows = []
    path = os.path.join(RESULTS, out_name + ".json")
    for c, t in enumerate(tg):
        rows.append(run_target(t, n_start=n_start))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):      # checkpoint: S14 lost a study to this
            with open(path, "w") as fh:
                json.dump({"rows": rows, "n_done": len(rows), "n_expected": len(tg)}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    fail = np.isin(pdbs, I.FAIL18)
    from s14 import ladder as L
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)

    out = {"n": len(rows), "n_start": n_start, "rows": rows,
           "incumbent": float(ref.mean()), "arms": {}}
    for name in rows[0]["arms"]:
        v = np.asarray([r["arms"][name]["rmsd"] for r in rows], float)
        out["arms"][name] = {**I.summary(v), "FAIL18": float(v[fail].mean()),
                             "median": float(np.median(v)),
                             "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    st = np.asarray([r["start_rmsd_retrieval_circmean"] for r in rows], float)
    out["start_retrieval_circmean"] = {**I.summary(st), "FAIL18": float(st[fail].mean())}

    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_distgeo", out, n_expected=len(tg))

    print(f"\nincumbent {ref.mean():.3f}   start (retrieval circ-mean) "
          f"{st.mean():.3f}\n")
    print(f"{'arm':<26}{'mean':>8}{'median':>8}{'sd':>7}{'<2A':>7}{'FAIL18':>9}"
          f"{'vs incumbent':>26}")
    for name, a in out["arms"].items():
        v = a["vs_incumbent"]
        print(f"{name:<26}{a['mean']:>8.3f}{a['median']:>8.3f}{a['sd']:>7.3f}"
              f"{a['frac_under_2.0']:>7.2f}{a['FAIL18']:>9.3f}"
              f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    return out


if __name__ == "__main__":
    run()
