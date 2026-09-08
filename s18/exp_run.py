"""s18/exp_run.py -- THE 126-TARGET TEST THAT DECIDES SPRINT 18.

The direct extension of `s17/refine.py`, with MATH's degree-1 object in place of the deployed
one.  Sprint 17 refined the coordinate average toward the full distance objective: the objective
fell 71% and RMSD ROSE 0.561 A [+0.407, +0.712], 31W/95L.  This asks whether the degree-1
truncation changes that, and prices the whole family in between:

    E_lam = E_le1 + lam * E_ge2 = (1-lam) E_le1 + lam E_full,  lam in {0,.25,.5,.75,1}
    PRE-REGISTERED in `s18/PREREG_exp.md`.  NOT extended.

Every arm starts from EXACTLY the same structure -- the coordinate average of the shipped top-75
(Control A) and its ideal-geometry projection.  The native is read only to score.

The objective, the reference measure, the grid and the quadrature are MATH's
(`s18/math_anova.py`), consumed through `s18/exp_obj.py`, which adds gradients and nothing else.

    python -m s18.exp_run sens        quadrature sensitivity, stratified subset (run FIRST)
    python -m s18.exp_run main        the 126-target experiment
    python -m s18.exp_run shuf        the shuffled-distogram control (a second fit per target)
    python -m s18.exp_run unif        the mu = uniform sensitivity arm
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                # noqa: E402
from s14 import avgspace as AV                 # noqa: E402
from s15 import align_lib as AL                # noqa: E402
from s15 import seed as SD                     # noqa: E402
from s18 import exp_obj as XO                  # noqa: E402
from s18 import math_anova as MA               # noqa: E402
from s18 import math_lib as ML                 # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
LAMBDAS = (0.0, 0.25, 0.5, 0.75, 1.0)          # PRE-REGISTERED.  Do not extend.
GRID = MA.GRID          # MATH owns this; read dynamically, never pinned here
NSAMP = MA.NSAMP        # ditto
BAND = 1.5
#: MATH's `fit` defaults to `batch = 32768`.  On this box that is the SLOW setting: measured
#: 38.8 s at 4096 against 63.1 s at 32768 for the same target, and 10^8 page faults per process
#: under multi-process load.  `batch` only chooses how many rows are chunked together -- every
#: reduction is within a row, and `step = batch // S` keeps each mesh point's S samples in one
#: chunk -- so the coefficients are bit-identical and only the wall clock changes.
FIT_BATCH = 4096
SENS_S = (256, 512, 2048)


# --------------------------------------------------------------------------------- utilities
def start_structure(t, nat):
    """Control A and its projection.  Identical construction to `s17/refine.py`."""
    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    W, PH, PS, _u = AV.top75_windows(pdb)
    W = np.asarray(W, float)
    P = I.pairwise_rmsd(W)
    avg, _b = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    return (W, np.asarray(PH, float), np.asarray(PS, float), np.asarray(avg, float),
            np.asarray(pr["ca"], float), np.asarray(pr["phi"], float),
            np.asarray(pr["psi"], float))


def _rank(x):
    return np.argsort(np.argsort(np.asarray(x, float))).astype(float)


def _sp(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(_rank(a), _rank(b))[0, 1])


def _pe(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _poa(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    da = a[:, None] - a[None, :]
    db = b[:, None] - b[None, :]
    m = np.triu(np.ones_like(da, bool), 1) & (db != 0) & (da != 0)
    return float((np.sign(da[m]) == np.sign(db[m])).mean()) if m.any() else float("nan")


def _target_obj(t, deb, mu, S, grid=GRID):
    """MATH's cached Target -> a gradient-carrying Obj.  Nothing redefined."""
    d = MA.gather_one(t, deb[int(t["fold"])])
    path = os.path.join(ML.CACHE, f"anova_{d['pdb']}_{mu}_{S}_{grid}.npz")
    if os.path.exists(path):
        try:
            return XO.Obj(MA.Target.load(path, d["pdb"])), d
        except IOError:
            os.remove(path)
    T = MA.Target(d["pdb"], d["seq"], d["n"], d["fold"],
                  d["dhat"], d["sd"], d["i"], d["j"]).fit(mu, S, grid, batch=FIT_BATCH)
    T.save(path)
    return XO.Obj(T), d


# ============================================================ QUADRATURE SENSITIVITY (first)
def sens(out=None, k=8):
    """Does the DOWNSTREAM OUTCOME depend on the quadrature size?  Decides S for the main run.

    The comparison is not "do the fields agree" -- it is "does the arm this sprint reports move
    by more than the 0.084 A MDE".  Reported per S, on a length-stratified subset.
    """
    out = out or os.path.join(RESULTS, "exp_sens.json")
    tg = I.targets()
    deb = MA.debias_map(tg)
    ns = np.array([t["n"] for t in tg])
    order = np.argsort(ns, kind="stable")
    pick = [tg[order[q]] for q in np.linspace(0, len(tg) - 1, k).astype(int)]
    rows = []
    t0 = time.time()
    for t in pick:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        _W, _PH, _PS, avg, _ca, phi0, psi0 = start_structure(t, nat)
        e = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
             "avg": float(I.ca_rmsd(avg, nat))}
        for S in SENS_S:
            ob, _d = _target_obj(t, deb, "pool", S)
            p, q, f, _nf, _ni = XO.lbfgs(ob.E_le1, phi0, psi0)
            e[f"le1_S{S}"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
            e[f"le1obj_S{S}"] = f
            ap, aq = ob.argmin_le1(hold=(phi0, psi0))
            e[f"argmin_S{S}"] = float(I.ca_rmsd(I.build_ca(ap, aq), nat))
            e[f"E0_S{S}"] = float(ob.t.E0)
            e[f"fsd_S{S}"] = float(ob.t.f.std())
            print(f"  {pdb} S={S} le1 {e[f'le1_S{S}']:.3f} argmin {e[f'argmin_S{S}']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        rows.append(e)
        json.dump({"rows": rows, "complete": False, "S": list(SENS_S)}, open(out, "w"))
    json.dump({"rows": rows, "complete": len(rows) == len(pick), "S": list(SENS_S),
               "n_expected": len(pick), "n_rows": len(rows)}, open(out, "w"))
    g = lambda kk: np.array([r[kk] for r in rows], float)         # noqa: E731
    print("\n  QUADRATURE SENSITIVITY -- mean over the stratified subset")
    print(f"    {'S':>6}{'refine_le1':>12}{'argmin_le1':>12}{'vs S=2048 (refine)':>21}"
          f"{'vs S=2048 (argmin)':>20}")
    for S in SENS_S:
        print(f"    {S:>6}{g(f'le1_S{S}').mean():>12.3f}{g(f'argmin_S{S}').mean():>12.3f}"
              f"{np.abs(g(f'le1_S{S}') - g('le1_S2048')).mean():>21.3f}"
              f"{np.abs(g(f'argmin_S{S}') - g('argmin_S2048')).mean():>20.3f}")
    print("\n  DECISION RULE (pre-registered): run at the smallest S whose mean |delta| to")
    print("  S = 2048 on BOTH arms is below the 0.084 A MDE.")
    return rows


# ============================================================================== MAIN EXPERIMENT
def main(S=NSAMP, mu="pool", out=None, targets=None, reverse=True, subset=None):
    """`subset` (an int) runs a LENGTH-STRATIFIED subsample -- used only for the mu-sensitivity
    arm, which needs a second full fit per target and does not fit the compute envelope at
    n = 126.  The leave-fold-out debias is ALWAYS fitted on the full 126, never on the subset:
    fitting it on a subset silently changes `dhat` and makes the arm incomparable, which is a
    trap this workstream already walked into once on a 3-target smoke."""
    out = out or os.path.join(RESULTS, f"exp_main_{mu}_{S}.json")
    full = I.targets()
    deb = MA.debias_map(full)
    tg = targets if targets is not None else full
    if subset:
        ns = np.array([t["n"] for t in tg])
        order = np.argsort(ns, kind="stable")
        tg = [tg[q] for q in np.sort(order[np.linspace(
            0, len(tg) - 1, int(subset)).astype(int)])]
        out = out.replace(".json", f"_sub{len(tg)}.json")
    rows = []
    if os.path.exists(out):
        try:
            p = json.load(open(out))
            if p.get("S") == S and p.get("mu") == mu and not p.get("complete"):
                rows = p["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()

    #: MATH's own `build` is populating the SAME coefficient cache in pinned order while this
    #: runs.  So the target order is chosen DYNAMICALLY at each step: take any target whose
    #: coefficients are already on disk (free), and when none is, take the LAST uncached one --
    #: so the two processes converge from opposite ends instead of recomputing each other's
    #: work.  Order affects no result: every target is independent and the artefact is keyed
    #: by pdb.  It only decides the wall clock.
    def _next(remaining):
        for x in remaining:                       # cheapest first: already cached
            if os.path.exists(os.path.join(
                    ML.CACHE, f"anova_{x['pdb']}_{mu}_{S}_{GRID}.npz")):
                return x
        return remaining[-1] if reverse else remaining[0]

    remaining = [x for x in tg if x["pdb"] not in done]
    while remaining:
        t = _next(remaining)
        remaining = [x for x in remaining if x["pdb"] != t["pdb"]]
        pdb, n = t["pdb"], int(t["n"])
        if pdb in done:
            continue
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        W, PH, PS, avg, ca, phi0, psi0 = start_structure(t, nat)
        th0 = np.concatenate([phi0, psi0])
        rng = SD.stable_rng(pdb, "s18exprun", mu, S)
        ob, d = _target_obj(t, deb, mu, S)

        e = {"pdb": pdb, "n": n, "fold": int(t["fold"]),
             "avg": float(I.ca_rmsd(avg, nat)), "proj": float(I.ca_rmsd(ca, nat)),
             "E0": float(ob.t.E0),
             "undetermined": ob.undetermined().tolist(),
             "obj_full_proj": ob.E_full(phi0, psi0)[0],
             "obj_le1_proj": ob.E_le1(phi0, psi0)[0]}

        #: ------------------------------------------------------------- THE LAMBDA LADDER
        for lam in LAMBDAS:
            k = f"lam{int(round(lam * 100)):03d}"
            p, q, f, nf, nit = XO.lbfgs(lambda a, b, _l=lam: ob.E_lambda(a, b, _l), phi0, psi0)
            e[k] = float(I.ca_rmsd(I.build_ca(p, q), nat))
            e[k + "_objlam"] = f
            e[k + "_objfull"] = ob.E_full(p, q)[0]
            e[k + "_objle1"] = ob.E_le1(p, q)[0]
            e[k + "_disp"] = float(np.linalg.norm(AL.wrap(np.concatenate([p, q]) - th0)))
            e[k + "_nfev"] = nf
            e[k + "_nit"] = nit

        #: -------------------------------- the EXACT separable global optimum, both variants
        ap, aq = ob.argmin_le1()                          # MATH's, terminals at a torus corner
        e["argmin_le1"] = float(I.ca_rmsd(I.build_ca(ap, aq), nat))
        e["argmin_le1_objle1"] = ob.E_le1(ap, aq)[0]
        e["argmin_le1_objfull"] = ob.E_full(ap, aq)[0]
        e["argmin_le1_disp"] = float(np.linalg.norm(AL.wrap(np.concatenate([ap, aq]) - th0)))
        bp, bq = ob.argmin_le1(hold=(phi0, psi0))         # terminals held at the start
        e["argmin_le1_hold"] = float(I.ca_rmsd(I.build_ca(bp, bq), nat))
        e["argmin_le1_hold_objfull"] = ob.E_full(bp, bq)[0]
        #: a local descent from the START on E_le1 vs the GLOBAL optimum -- the basin question
        e["argmin_gap"] = e["argmin_le1_objle1"] - e["lam000_objle1"]

        #: --------------------------------- MATH's harmonic truncation (order-1 in ANGLES)
        Th = ob.t.harmonic_truncate(1)
        obh = XO.Obj(Th)
        p, q, f, _nf, _ni = XO.lbfgs(obh.E_le1, phi0, psi0)
        e["refine_h1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
        hp, hq = obh.argmin_le1(hold=(phi0, psi0))
        e["argmin_h1"] = float(I.ca_rmsd(I.build_ca(hp, hq), nat))

        #: ------------------------------ MATH's SUB-RESIDUE (angle-additive) object.
        #: `E_0 + sum_r a_r(phi_r) + sum_r b_r(psi_r)` -- strictly coarser than the
        #: residue-additive object, and the closest legitimate continuous relative of the
        #: lattice's per-qubit truncation.  MATH's own note is that NO continuous object
        #: equals the strict Walsh weight-<=1 projection, which is the content of F5.
        if getattr(ob.t, "A1", None) is not None:
            oba = XO.ObjAng(ob.t)
            p, q, f, _nf, _ni = XO.lbfgs(oba.E, phi0, psi0)
            e["refine_ang"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
            e["refine_ang_objfull"] = ob.E_full(p, q)[0]
            ap2, aq2 = oba.argmin(hold=(phi0, psi0))
            e["argmin_ang"] = float(I.ca_rmsd(I.build_ca(ap2, aq2), nat))
            e["obj_ang_proj"] = oba.E(phi0, psi0)[0]

        #: -------------------------------------------------------------------- CONTROLS ---
        #: C3 zero-information additive field: matched per-residue amplitude AND matched power
        #: spectrum, so form, smoothness and magnitude are held and only information is removed.
        FR = ob.rand_field(rng)
        p, q, f, _nf, _ni = XO.lbfgs(lambda a, b: ob.E_le1(a, b, F=FR, E0=0.0), phi0, psi0)
        e["rand_field_le1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
        e["rand_field_disp"] = float(np.linalg.norm(AL.wrap(np.concatenate([p, q]) - th0)))
        rp, rq = ob.argmin_le1(F=FR, hold=(phi0, psi0))
        e["rand_field_argmin"] = float(I.ca_rmsd(I.build_ca(rp, rq), nat))

        #: C2 matched-magnitude random torsion moves, one per arm they are matched to
        for k, disp in (("le1", e["lam000_disp"]), ("full", e["lam100_disp"]),
                        ("argmin", e["argmin_le1_disp"])):
            z = rng.standard_normal(2 * n)
            z *= disp / max(float(np.linalg.norm(z)), 1e-12)
            th = th0 + z
            e["rand_move_" + k] = float(I.ca_rmsd(I.build_ca(th[:n], th[n:]), nat))

        #: ------------------------- FREE ARMS FROM MATH'S MOMENT TABLES (`with_dhat`).
        #: `_set_dhat` rebuilds E_0 and every field EXACTLY from the conditional moment
        #: tables with NO new chain builds, so these cost nothing and were pre-registered
        #: (PREREG H7) before any 126-target number was read.
        try:
            dtrue = np.sqrt(((nat[d["i"]] - nat[d["j"]]) ** 2).sum(1))
            #: CONTROL C1 -- the deployed distogram with its PAIRS PERMUTED.  The weights w
            #: are left in place, so only the information is destroyed, not the weighting.
            perm = rng.permutation(len(d["dhat"]))
            obs = XO.Obj(ob.t.with_dhat(np.asarray(d["dhat"], float)[perm]))
            p, q, f, _nf, _ni = XO.lbfgs(obs.E_le1, phi0, psi0)
            e["shufobj_le1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
            sp2, sq2 = obs.argmin_le1(hold=(phi0, psi0))
            e["shufobj_argmin"] = float(I.ca_rmsd(I.build_ca(sp2, sq2), nat))
            p, q, f, _nf, _ni = XO.lbfgs(obs.E_full, phi0, psi0)
            e["shufobj_full"] = float(I.ca_rmsd(I.build_ca(p, q), nat))

            #: ORACLE DIAGNOSTIC -- degree-1 fed the NATIVE's own distances.  Answers whether
            #: degree-1 is limited by the DISTANCES or by the TRUNCATION.  Never a result.
            obo = XO.Obj(ob.t.with_dhat(dtrue))
            p, q, f, _nf, _ni = XO.lbfgs(obo.E_le1, phi0, psi0)
            e["ORC_le1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
            op2, oq2 = obo.argmin_le1(hold=(phi0, psi0))
            e["ORC_argmin"] = float(I.ca_rmsd(I.build_ca(op2, oq2), nat))
            p, q, f, _nf, _ni = XO.lbfgs(obo.E_full, phi0, psi0)
            e["ORC_full"] = float(I.ca_rmsd(I.build_ca(p, q), nat))

            #: ORACLE DIAGNOSTIC -- the coordinator's construction read through degree-1:
            #: residual MAGNITUDES kept, DIRECTION destroyed.
            r_ = np.asarray(d["dhat"], float) - dtrue
            sh = np.maximum(dtrue + rng.permutation(r_)
                            * rng.choice([-1.0, 1.0], size=len(r_)), 2.0)
            obh2 = XO.Obj(ob.t.with_dhat(sh))
            p, q, f, _nf, _ni = XO.lbfgs(obh2.E_le1, phi0, psi0)
            e["ORCSHUF_le1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
            p, q, f, _nf, _ni = XO.lbfgs(obh2.E_full, phi0, psi0)
            e["ORCSHUF_full"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
        except (AttributeError, TypeError) as exc:
            e["with_dhat_unavailable"] = str(exc)

        #: ---------------------------------------- PHASE 9 -- selection over the SAME pool
        rw = np.asarray(I.kabsch_rmsd_batch(W, nat), float)
        ef = np.asarray(ob.t.E_full(PH, PS), float)
        el = np.asarray(ob.t.E_le1(PH, PS), float)
        e["pool_best"] = float(rw.min())
        e["pool_mean"] = float(rw.mean())
        e["sel_full"] = float(rw[int(np.argmin(ef))])
        e["sel_le1"] = float(rw[int(np.argmin(el))])
        band = rw <= (rw.min() + BAND)
        e["n_band"] = int(band.sum())
        for nm, v in (("full", ef), ("le1", el)):
            e[f"sp_{nm}"] = _sp(v, rw)
            e[f"pe_{nm}"] = _pe(v, rw)
            e[f"poa_{nm}"] = _poa(v, rw)
            e[f"spband_{nm}"] = _sp(v[band], rw[band]) if band.sum() >= 5 else float("nan")
        #: ensemble diversity of the pool the arms consume
        e["pool_div"] = float(I.pairwise_rmsd(W)[np.triu_indices(len(W), 1)].mean())

        rows.append(e)
        json.dump({"rows": rows, "complete": False, "S": S, "mu": mu, "grid": GRID,
                   "lambdas": list(LAMBDAS), "math_config": MA.CONFIG}, open(out, "w"))
        if len(rows) % 5 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s) {pdb} "
                  f"avg {e['avg']:.2f} lam0 {e['lam000']:.2f} lam1 {e['lam100']:.2f}", flush=True)

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "S": S, "mu": mu, "grid": GRID,
               "lambdas": list(LAMBDAS), "math_config": MA.CONFIG,
               "n_expected": len(tg), "n_rows": len(rows)}, open(out, "w"))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


# ================================================== the shuffled-distogram control (2nd fit)
def shuf(S=NSAMP, out=None, targets=None):
    """CONTROL C1: the identical degree-1 machinery against a SHUFFLED distogram.

    Decides whether any degree-1 gain is the OBJECTIVE'S INFORMATION or merely the smoothing
    that averaging a penalty over a reference ensemble applies to any restraint set at all.
    """
    out = out or os.path.join(RESULTS, f"exp_shuf_{S}.json")
    tg = targets if targets is not None else I.targets()
    deb = MA.debias_map(tg)
    rows = []
    if os.path.exists(out):
        try:
            p = json.load(open(out))
            if not p.get("complete"):
                rows = p["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        pdb, n = t["pdb"], int(t["n"])
        if pdb in done:
            continue
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        _W, _PH, _PS, avg, _ca, phi0, psi0 = start_structure(t, nat)
        rng = SD.stable_rng(pdb, "s18expshuf", S)
        d = MA.gather_one(t, deb[int(t["fold"])])
        perm = rng.permutation(len(d["dhat"]))
        path = os.path.join(ML.CACHE, f"anovaSHUF_{pdb}_pool_{S}_{GRID}.npz")
        if os.path.exists(path):
            T = MA.Target.load(path, pdb)
        else:
            T = MA.Target(pdb, d["seq"], d["n"], d["fold"], d["dhat"][perm], d["sd"][perm],
                          d["i"], d["j"]).fit("pool", S, GRID)
            T.save(path)
        ob = XO.Obj(T)
        e = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "avg": float(I.ca_rmsd(avg, nat))}
        p, q, f, _nf, _ni = XO.lbfgs(ob.E_le1, phi0, psi0)
        e["shufobj_le1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
        ap, aq = ob.argmin_le1(hold=(phi0, psi0))
        e["shufobj_argmin"] = float(I.ca_rmsd(I.build_ca(ap, aq), nat))
        p, q, f, _nf, _ni = XO.lbfgs(ob.E_full, phi0, psi0)
        e["shufobj_full"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
        rows.append(e)
        json.dump({"rows": rows, "complete": False, "S": S}, open(out, "w"))
        if len(rows) % 5 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "complete": len(rows) == len(tg), "S": S,
               "n_expected": len(tg), "n_rows": len(rows)}, open(out, "w"))
    print(f"DONE shuf {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sens"
    Sarg = int(sys.argv[2]) if len(sys.argv) > 2 else NSAMP
    if mode == "sens":
        sens()
    elif mode == "main":
        main(S=Sarg)
    elif mode == "unif":
        main(S=Sarg, mu="uniform", subset=int(sys.argv[3]) if len(sys.argv) > 3 else 30)
    elif mode == "shuf":
        shuf(S=Sarg)
    else:
        raise SystemExit(mode)
