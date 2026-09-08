"""SPRINT 15, ALIGN -- THE CONTROLLED TEST OF THE GEOMETRY: same shrinkage, different metric.

`s15/align_fit.py` tests HARD interventions (reweight restraints, restrict to a subspace).
This module tests the single cleanest version of the hypothesis, because it changes ONLY the
geometry and holds the amount of regularisation fixed.

Add to the restraint objective a quadratic penalty on the fit's displacement from its own
start,  `P = lam * delta' A delta`,  with `delta = theta - theta_start`, and use three
different `A`, all normalised to `trace(A) = 2n` so that on an isotropic random displacement
they impose EXACTLY the same expected penalty:

    tik   A = I                                    isotropic -- the NULL.  Ordinary Tikhonov:
                                                   every torsion direction costs the same.
    loud  A = J'J / s_rms^2                        penalises motion in the RMSD-LOUD directions,
                                                   i.e. confines the fit to the quiet subspace.
    quiet A propto (s_max^2 I - J'J)               penalises motion in the RMSD-QUIET directions,
                                                   i.e. confines the fit to the well-conditioned
                                                   subspace.  This is the SMOOTH version of the
                                                   brief's "project onto the top-r singular
                                                   directions" arm, without the hard cut that
                                                   `align_jac` showed is fragile (the bottom
                                                   subspaces of J(native) and J(emitted) share
                                                   their bulk, mean cos^2 0.827, but their
                                                   minimum principal cosine is only 0.170).

If the alignment geometry is a real lever, `loud` or `quiet` must beat `tik` at matched
shrinkage.  If all three land on top of each other, the only thing that mattered was
shrinkage -- ordinary regularisation -- and the geometry is decoration.  That is a clean
falsification and it is the point of the module.

`J` is evaluated at the stage-1 (`squared`) solution, so every arm is NATIVE-FREE.
`lam` is chosen LEAVE-FOLD-OUT afterwards, exactly as in `s15/align_fit.py`.

    python -m s15.align_reg
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from scipy.optimize import minimize          # noqa: E402
from core import project as pj               # noqa: E402
from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import info_lib as L                # noqa: E402
from s15 import seed as SD                   # noqa: E402
import peptide_db as pdb                     # noqa: E402

PATH = os.path.join(ROOT, "s15", "results", "align_reg.json")
#: lam ladder.  The top of the ladder must drive every family back to its START (the fit is
#: frozen), so a U-shape is expected and the leave-fold-out pick must land in the interior.
LAMS = (0.03, 0.3, 3.0, 30.0)
FAMS = ("tik", "loud", "quiet")


def metrics(J, m):
    """The three `A` matrices, each normalised to `trace(A) = m`."""
    JtJ = J.T @ J
    tr = float(np.trace(JtJ))
    smax2 = float(np.linalg.eigvalsh(JtJ)[-1])
    Aq = smax2 * np.eye(m) - JtJ
    trq = float(np.trace(Aq))
    return {"tik": np.eye(m),
            "loud": JtJ * (m / tr) if tr > 0 else np.eye(m),
            "quiet": Aq * (m / trq) if trq > 0 else np.eye(m)}


def fit_reg(dhat, sd, i, j, phi0, psi0, w, Amat, lam, maxiter=400):
    n = len(phi0)
    inv = 1.0 / np.asarray(sd, float)
    t0 = np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)])

    def fg(x):
        phi = x[:n]
        psi = x[n:]
        G = pj.frames(phi[None], psi[None])[0]
        CA = np.asarray(pj.build_ca_exact(phi[None], psi[None]), float)[0]
        rv = CA[i] - CA[j]
        d = np.maximum(np.sqrt((rv * rv).sum(1)), 1e-9)
        r = (d - dhat) * inv
        f = float((w * r * r).sum())
        coef = ((w * 2.0 * r * inv) / d)[:, None] * rv
        gCA = np.zeros_like(CA)
        np.add.at(gCA, i, coef)
        np.add.at(gCA, j, -coef)
        g = pj._torsion_grad(G, CA, gCA)
        dl = x - t0
        Ad = Amat @ dl
        return f + lam * float(dl @ Ad), g + 2.0 * lam * Ad

    res = minimize(fg, t0.copy(), jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    #: report the DATA term only, so the objective-only multi-start selection is comparable
    x = res.x
    dl = x - t0
    fdata = float(res.fun) - lam * float(dl @ (Amat @ dl))
    return x[:n], x[n:], fdata


def run(limit=None, n_start=6):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)
    folds = sorted({data[p]["fold"] for p in pdbs})
    deb = {}
    for f in folds:
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        p = t["pdb"]
        d = data[p]
        n = d["n"]
        i, j, sd = d["i"], d["j"], d["sd"]
        dhat = np.maximum(d["dhat"] - deb[d["fold"]](d["sep"]), 2.0)
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        w0 = 1.0 / (sd * sd)
        S = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p))
        row = {"pdb": p, "n": n, "fold": d["fold"], "arms": {}}

        def emit(name, phi, psi):
            row["arms"][name] = A.score_axes(phi, psi, nphi, npsi, d["nat"], dhat, sd, i, j,
                                             d["dtrue"])

        bb = None
        for phi0, psi0, _tg in S:
            ph, ps, fv = A.fit(dhat, sd, i, j, phi0, psi0, wpair=w0)
            if bb is None or fv < bb[2]:
                bb = (ph, ps, fv, np.asarray(phi0, float), np.asarray(psi0, float))
        emit("squared", bb[0], bb[1])
        emit("START_of_winner", bb[3], bb[4])          # the lam -> infinity limit
        Jsup, _, _ = A.sup_jacobian(bb[0], bb[1])
        Am = metrics(Jsup, 2 * n)

        for fam in FAMS:
            for lam in LAMS:
                best = None
                for phi0, psi0, _tg in S:
                    ph, ps, fv = fit_reg(dhat, sd, i, j, phi0, psi0, w0, Am[fam], lam)
                    if best is None or fv < best[0]:
                        best = (fv, ph, ps)
                emit(f"{fam}_l{lam}", best[1], best[2])
        rows.append(row)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(PATH, "w") as fh:
                json.dump({"rows": rows, "n_done": len(rows)}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed ({(time.time()-t0)/60:.1f} min)",
                  flush=True)
    return summarise(rows)


def summarise(rows):
    from s14 import ladder as LD
    from s15.align_fit import _lfo
    pdbs = [r["pdb"] for r in rows]
    fol = np.asarray([r["fold"] for r in rows], int)
    fail = np.isin(pdbs, I.FAIL18)
    ref = np.asarray([LD.incumbent_rmsd()[p] for p in pdbs], float)
    keys = ("rmsd", "tors_rms_deg", "dist_mae", "resid_medz", "align", "align_self")
    names = sorted(rows[0]["arms"])
    arms = {a: {k: np.asarray([r["arms"][a][k] for r in rows], float) for k in keys}
            for a in names}
    picks = {}
    for fam in FAMS:
        mem = [f"{fam}_l{l}" for l in LAMS]
        v, infold, pk = _lfo(rows, mem)
        picks[fam] = pk
        arms[fam + "_LFO"] = {"rmsd": v}
        arms[fam + "_INFOLD_ORACLE"] = {"rmsd": infold}
        for k in keys[1:]:
            M = {m: np.asarray([r["arms"][m][k] for r in rows], float) for m in mem}
            arms[fam + "_LFO"][k] = np.asarray([M[pk[int(f)]][ii]
                                                for ii, f in enumerate(fol)], float)
            arms[fam + "_INFOLD_ORACLE"][k] = arms[fam + "_LFO"][k]
    base = arms["squared"]["rmsd"]
    out = {"n": len(rows), "picks": picks, "rows": rows, "arms": {}}
    for a in arms:
        v = arms[a]["rmsd"]
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          **{k: float(np.mean(arms[a][k])) for k in keys[1:] if k in arms[a]},
                          "vs_squared": L.report(v, base, folds=fol, names=pdbs),
                          "vs_incumbent": I.paired(v, ref, folds=fol, names=pdbs)}
    with open(PATH, "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"\nn = {len(rows)}   MATCHED-SHRINKAGE TEST (trace(A) = 2n in every family)\n")
    print(f"{'arm':<24}{'RMSD':>8}{'med':>8}{'tors':>7}{'dMAE':>7}{'align':>7}"
          f"{'   vs squared':<24}")
    for a in sorted(out["arms"], key=lambda z: out["arms"][z]["mean"]):
        s = out["arms"][a]
        v = s["vs_squared"]
        print(f"{a:<24}{s['mean']:>8.3f}{s['median']:>8.3f}{s['tors_rms_deg']:>7.1f}"
              f"{s['dist_mae']:>7.3f}{s['align']:>7.3f}   {v['mean_diff']:+.3f} "
              f"[{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] W/L {v['n_better']}/{v['n_worse']}")
    print("\npicks (leave-fold-out): " + json.dumps(picks))
    return out


if __name__ == "__main__":
    run(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
