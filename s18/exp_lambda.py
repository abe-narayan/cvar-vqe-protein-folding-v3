"""RETIRED -- DO NOT READ AS A RESULT.

This is the PROVISIONAL degree-1 implementation written against the interface published in
`s18/COORD_exp_to_math.md` while MATH had not yet landed a module.  MATH landed
`s18/math_anova.py`, whose object is RESIDUE-additive and differs from what this file built,
and every reported Sprint-18 EXPERIMENT arm consumes MATH's object through `s18/exp_obj.py`.
Nothing here contributed a single number to `s18/exp_FINDINGS.md`.  Its only artefact
(`_RETIRED_exp_lambda_provisional.json`, 3 rows from a 3-target smoke whose leave-fold-out
debias was fitted on those 3 targets and is therefore not comparable to anything) is renamed
so it cannot be mistaken for a result.  Kept only as the record of what was built and when.

s18/exp_lambda.py -- THE 126-TARGET TEST THAT DECIDES SPRINT 18.

The direct extension of `s17/refine.py`.  Sprint 17 asked whether moving the coordinate average
toward the DEPLOYED distance objective helps; the answer was no by +0.561 A [+0.407, +0.712] at
a 71% objective reduction.  This module asks whether the DEGREE-1 TRUNCATION of that objective
changes the answer, and prices the whole family in between.

    E_lambda = E_le1 + lambda * E_ge2 = (1-lambda) E_le1 + lambda E_full,
    lambda in {0, 0.25, 0.5, 0.75, 1}   -- PRE-REGISTERED, not extended.

Every arm starts from EXACTLY the same structure: the coordinate average of the shipped top-75
(Control A) and its ideal-geometry projection.  Nothing here reads the native except to score.

ARMS.  See `s18/PREREG_exp.md` for the pre-registered hypothesis, control, criterion and
falsifier of each.

    avg / proj                 Control A and its torsion start point
    lam{0,25,50,75,100}        the ladder, L-BFGS with analytic gradient, same tolerances
    refine_res                 the residue-additive object (adds the intra-residue weight-2 mass)
    argmin_le1 / argmin_res    the EXACT separable global optima -- no search, no budget
    CONTROLS
    rand_obj_le1 / _full       the identical refinement against a SHUFFLED distogram
    rand_move_le1 / _full      a random torsion move of MATCHED displacement
    rand_tab_le1               a random additive surrogate of MATCHED per-coordinate amplitude
                               -- the zero-information control specific to an ADDITIVE object
    SELECTION (Phase 9, protecting the coordinate average)
    sel_full / sel_le1         argmin of each objective OVER THE 75 POOL WINDOWS
    ALIGNMENT (BRIEF sec 5, four named axes, never blurred)
    per-target Spearman, Pearson, pairwise-ordering accuracy and IN-BAND Spearman of each
    objective against RMSD over the same 75 windows.

Run:  python -m s18.exp_lambda            (full, 126 targets, checkpointed every target)
      python -m s18.exp_lambda report
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import scipy.optimize as SO                  # noqa: E402

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s18 import exp_anova as AN              # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "exp_lambda.json")
LAMBDAS = (0.0, 0.25, 0.5, 0.75, 1.0)          # PRE-REGISTERED.  Do not extend.
MAXITER = 400
BAND = 1.5


def _lbfgs(fun, phi0, psi0, maxiter=MAXITER):
    """One L-BFGS-B descent with an analytic gradient.  Returns (phi, psi, f, nfev, nit)."""
    n = len(phi0)

    def fg(x):
        f, g = fun(x[:n], x[n:])
        return f, np.asarray(g, float)

    r = SO.minimize(fg, np.concatenate([phi0, psi0]), jac=True, method="L-BFGS-B",
                    options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return r.x[:n], r.x[n:], float(r.fun), int(r.nfev), int(r.nit)


def _rank(x):
    o = np.argsort(np.argsort(np.asarray(x, float)))
    return o.astype(float)


def _spear(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(_rank(a), _rank(b))[0, 1])


def _pears(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _poa(a, b):
    """Pairwise ordering accuracy: fraction of pairs the objective orders as RMSD does."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    da = a[:, None] - a[None, :]
    db = b[:, None] - b[None, :]
    m = np.triu(np.ones_like(da, bool), 1) & (db != 0) & (da != 0)
    if not m.any():
        return float("nan")
    return float((np.sign(da[m]) == np.sign(db[m])).mean())


def _obj_batch(ob, PHI, PSI):
    """`E_full` for a stack of candidate torsion sets -- used to score the pool."""
    return ob._batch(np.concatenate([np.asarray(PHI, float), np.asarray(PSI, float)], 1))


def run(targets=None, out=OUT):
    tg = targets if targets is not None else I.targets()
    cfg, chash = AN.config_hash()
    src = "math" if AN.load_math() is not None else "s18.exp_anova(provisional)"
    rows = []
    if os.path.exists(out):
        try:
            prev = json.load(open(out))
            if prev.get("config_hash") == chash and not prev.get("complete"):
                rows = prev["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()

    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    deb = {}
    for f in sorted({int(t["fold"]) for t in tg}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    for c, t in enumerate(tg):
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        if pdb in done:
            continue
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        rng = SD.stable_rng(pdb, "s18explambda")

        W, PH, PS, _u = AV.top75_windows(pdb)
        W = np.asarray(W, float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0 = np.asarray(pr["phi"], float)
        psi0 = np.asarray(pr["psi"], float)
        th0 = np.concatenate([phi0, psi0])

        ob = AN.build(dhat, sd, i, j, PH, PS, pdb, mu="pool")
        obu = AN.build(dhat, sd, i, j, PH, PS, pdb, mu="uniform")
        perm = rng.permutation(len(dhat))
        obs = AN.build(dhat[perm], sd[perm], i, j, PH, PS, pdb + "_shuf", mu="pool")

        e = {"pdb": pdb, "n": n, "fold": fold, "seq": seq,
             "avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
             "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
             "E0": ob.E0, "E0_unif": obu.E0,
             "split_rms": ob.split_rms, "split_r": ob.split_r,
             "split_rms_unif": obu.split_rms, "split_r_unif": obu.split_r}

        #: the objective values at the START, both objects, so reduction is measurable
        e["obj_full_proj"] = ob.E_full(phi0, psi0)[0]
        e["obj_le1_proj"] = ob.E_le1(phi0, psi0)[0]
        e["obj_res_proj"] = ob.E_res(phi0, psi0)[0]
        e["obj_le1u_proj"] = obu.E_le1(phi0, psi0)[0]

        #: ---------------------------------------------------------------- THE LAMBDA LADDER
        for lam in LAMBDAS:
            k = f"lam{int(round(lam * 100)):03d}"
            ob.nev = 0
            pl, ql, fl, nf, nit = _lbfgs(lambda p, q, _l=lam: ob.E_lambda(p, q, _l), phi0, psi0)
            e[k] = float(I.ca_rmsd(I.build_ca(pl, ql), nat))
            e[k + "_objlam"] = fl
            e[k + "_objfull"] = ob.E_full(pl, ql)[0]
            e[k + "_objle1"] = ob.E_le1(pl, ql)[0]
            e[k + "_disp"] = float(np.linalg.norm(A.wrap(np.concatenate([pl, ql]) - th0)))
            e[k + "_nfev"] = nf
            e[k + "_nit"] = nit

        #: ------------------------------------------------ the residue-additive object (B3/B5)
        pr_, qr_, fr_, nf, nit = _lbfgs(ob.E_res, phi0, psi0)
        e["refine_res"] = float(I.ca_rmsd(I.build_ca(pr_, qr_), nat))
        e["refine_res_objres"] = fr_
        e["refine_res_objfull"] = ob.E_full(pr_, qr_)[0]
        e["refine_res_disp"] = float(np.linalg.norm(A.wrap(np.concatenate([pr_, qr_]) - th0)))
        e["refine_res_nfev"] = nf

        #: ------------------------------------- the EXACT separable global optima (B4/B5) --
        #: no search, no budget: these are the continuous analogue of the 19-target CERTIFIED
        #: degree-1 argmin that Phase 0 recomputed, and the arm the sprint's premise is about.
        ap, aq = ob.argmin_le1()
        e["argmin_le1"] = float(I.ca_rmsd(I.build_ca(ap, aq), nat))
        e["argmin_le1_objle1"] = ob.E_le1(ap, aq)[0]
        e["argmin_le1_objfull"] = ob.E_full(ap, aq)[0]
        e["argmin_le1_disp"] = float(np.linalg.norm(A.wrap(np.concatenate([ap, aq]) - th0)))
        bp, bq = ob.argmin_res()
        e["argmin_res"] = float(I.ca_rmsd(I.build_ca(bp, bq), nat))
        e["argmin_res_objres"] = ob.E_res(bp, bq)[0]
        e["argmin_res_objfull"] = ob.E_full(bp, bq)[0]

        #: --------------------------------------------------- mu-SENSITIVITY (F1 / PREREG H6)
        pu, qu, fu, _nf, _ni = _lbfgs(obu.E_le1, phi0, psi0)
        e["lam000_unif"] = float(I.ca_rmsd(I.build_ca(pu, qu), nat))
        au, bu = obu.argmin_le1()
        e["argmin_le1_unif"] = float(I.ca_rmsd(I.build_ca(au, bu), nat))

        #: --------------------------------------------------------------------- CONTROLS ----
        #: C1 shuffled distogram, through the IDENTICAL ANOVA machinery, both objects
        ps_, qs_, _f, _nf, _ni = _lbfgs(obs.E_le1, phi0, psi0)
        e["rand_obj_le1"] = float(I.ca_rmsd(I.build_ca(ps_, qs_), nat))
        e["rand_obj_le1_disp"] = float(np.linalg.norm(A.wrap(np.concatenate([ps_, qs_]) - th0)))
        sp_, sq_ = obs.argmin_le1()
        e["rand_obj_argmin_le1"] = float(I.ca_rmsd(I.build_ca(sp_, sq_), nat))
        pf_, qf_, _f, _nf, _ni = _lbfgs(obs.E_full, phi0, psi0)
        e["rand_obj_full"] = float(I.ca_rmsd(I.build_ca(pf_, qf_), nat))

        #: C2 matched-magnitude random move, one per refinement arm it is matched to
        for k, disp in (("le1", e["lam000_disp"]), ("full", e["lam100_disp"]),
                        ("argmin_le1", e["argmin_le1_disp"])):
            z = rng.standard_normal(2 * n)
            z *= disp / max(float(np.linalg.norm(z)), 1e-12)
            th = th0 + z
            e["rand_move_" + k] = float(I.ca_rmsd(I.build_ca(th[:n], th[n:]), nat))

        #: C3 a RANDOM ADDITIVE surrogate of matched per-coordinate amplitude.  This is the
        #: zero-information control an ADDITIVE object needs: same functional form, same
        #: smoothness, same amplitude, no information about the target at all.  If refining
        #: toward THIS moves RMSD as much as refining toward `E_le1`, degree-1 carries nothing.
        amp = np.sqrt((ob.f1 ** 2).mean(1))
        Rt = rng.standard_normal((2 * n, ob.g1))
        Rt -= Rt.mean(1, keepdims=True)
        Rt *= (amp / np.maximum(np.sqrt((Rt ** 2).mean(1)), 1e-12))[:, None]
        cR = np.array([AN._trig_coef1(Rt[q]) for q in range(2 * n)])
        obr = AN.Anova.__new__(AN.Anova)
        obr.__dict__.update(ob.__dict__)
        obr.c1, obr.f1, obr.E0 = cR, Rt, 0.0
        prn, qrn, _f, _nf, _ni = _lbfgs(obr.E_le1, phi0, psi0)
        e["rand_tab_le1"] = float(I.ca_rmsd(I.build_ca(prn, qrn), nat))
        e["rand_tab_le1_disp"] = float(np.linalg.norm(A.wrap(np.concatenate([prn, qrn]) - th0)))
        anp, anq = obr.argmin_le1()
        e["rand_tab_argmin_le1"] = float(I.ca_rmsd(I.build_ca(anp, anq), nat))

        #: ------------------------------------------- PHASE 9 -- selection over the SAME pool
        Rw = I.kabsch_rmsd_batch(W, nat) if hasattr(I, "kabsch_rmsd_batch") else None
        rw = np.asarray(Rw, float) if Rw is not None else np.array(
            [I.ca_rmsd(W[q], nat) for q in range(len(W))])
        ef = _obj_batch(ob, PH, PS)
        el = np.array([ob.E_le1(PH[q], PS[q])[0] for q in range(len(PH))])
        er = np.array([ob.E_res(PH[q], PS[q])[0] for q in range(len(PH))])
        e["pool_best"] = float(rw.min())
        e["pool_mean"] = float(rw.mean())
        e["sel_full"] = float(rw[int(np.argmin(ef))])
        e["sel_le1"] = float(rw[int(np.argmin(el))])
        e["sel_res"] = float(rw[int(np.argmin(er))])
        e["sel_random"] = float(rw.mean())          # expectation of a uniform pick
        #: ALIGNMENT, four named axes, over the same 75 real structures
        band = rw <= (rw.min() + BAND)
        for nm, v in (("full", ef), ("le1", el), ("res", er)):
            e[f"sp_{nm}"] = _spear(v, rw)
            e[f"pe_{nm}"] = _pears(v, rw)
            e[f"poa_{nm}"] = _poa(v, rw)
            e[f"spband_{nm}"] = _spear(v[band], rw[band]) if band.sum() >= 5 else float("nan")
        e["n_band"] = int(band.sum())

        rows.append(e)
        json.dump({"rows": rows, "complete": False, "config": cfg, "config_hash": chash,
                   "anova_source": src, "lambdas": list(LAMBDAS)}, open(out, "w"))
        if (len(rows)) % 5 == 0:
            print(f"  {len(rows)}/{len(tg)}  ({time.time()-t0:.0f}s)  {pdb}", flush=True)

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "config": cfg,
               "config_hash": chash, "anova_source": src, "lambdas": list(LAMBDAS),
               "n_expected": len(tg), "n_rows": len(rows)}, open(out, "w"))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


# --------------------------------------------------------------------------------- statistics
def boot_fold(diff, folds, rng, B=4000):
    """Paired FOLD-CLUSTERED bootstrap: resample folds, not targets.  TARGET is the unit
    inside a fold, but targets in a fold are not independent, so the cluster is the fold."""
    diff = np.asarray(diff, float)
    folds = np.asarray(folds)
    uf = np.unique(folds)
    idx = [np.where(folds == f)[0] for f in uf]
    m = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, len(uf), len(uf))
        m[b] = diff[np.concatenate([idx[p] for p in pick])].mean()
    return float(diff.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def _load(path=OUT):
    o = json.load(open(path))
    if not o.get("complete"):
        print(f"  ** WARNING: {os.path.basename(path)} is PARTIAL "
              f"({len(o['rows'])} rows) -- not a final result **")
    return o


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        from s18 import exp_report
        exp_report.main()
    else:
        run()
