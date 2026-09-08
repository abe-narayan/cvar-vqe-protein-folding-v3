"""s18/phys_lambda.py -- PHASE 6: `leg_contact` AS A TERM IN THE OBJECTIVE.

Pre-registration: `s18/PREREG_phys.md` sections 1-2, written before this module produced a
number.  The headline hypothesis, and my honest prior that it will FAIL, are both recorded there.

    E(theta; lam_c)  =  E_le1_dist(theta)  +  lam_c * s * E_contact(theta)

`E_le1_dist` is consumed through the SINGLE SHARED degree-1 ANOVA definition (`s18.exp_anova`,
which delegates to `s18/math_*.py` the moment MATH exposes `build` + `__anova_interface__`).
**No second definition is built here** -- a second definition makes the comparison void, and the
module records in every artefact which implementation produced it, by config hash.

`E_contact` is genuine Legacy: `core.energy.contact_term` at `DEFAULT_WEIGHTS["contact"] = 1.0`,
never fitted, verified bit-for-bit against `s16.energy_lib` by gate G0 (max |dE| = 1.8e-15).

THE NORMALISATION, and why it is the part that can silently ruin the experiment.  The two terms
have unrelated units, so `lam_c` means nothing until their scales are matched, and matching them
on anything that has seen the native would make every arm ORACLE.  The pre-registered choice is
the per-target standard deviation of each term over the target's own **ranker-neutral K = 500
BLOSUM62 retrieval pool** -- fixed by retrieval, not by any score being compared:

    s = sd_pool(E_le1_dist) / sd_pool(E_contact)      so lam_c = 1 is one pool-sd against one

THE SIGN IS PRE-REGISTERED AND IS NOT TUNED.  Corrected MJ energies are negative for favourable
contacts and Legacy is minimised, so the physically motivated sign is `lam_c > 0`; that is also
the sign of the +0.080 in-band partial correlation.  It is NOT the sign the GLOBAL correlation
(rho = -0.176) would want, and that contradiction is registered as the chief reason the
hypothesis may fail.  Negative `lam_c` is a CONTROL, never a rescue.

THREE PASSES, because Sprint 17 proved these axes dissociate and the brief forbids merging them:

    align    objective ALIGNMENT.  Score the whole K = 500 pool under every objective; global
             and in-band Spearman, in-band selection against a matched-random selection, argmin.
             No optimisation at all.  Cheap, and it can falsify on its own.
    refine   LOCAL REFINEMENT and FINAL RMSD.  One L-BFGS descent per arm from the IDENTICAL
             start -- the ideal-geometry projection of the shipped coordinate average, which is
             the same start `s17/refine.py` used, so the numbers are comparable to L30.
    gen      CANDIDATE GENERATION, ENSEMBLE DIVERSITY, NEAR-NATIVE COVERAGE.  Multi-start
             descent from the first `M_START` top-75 windows in RETRIEVAL order (native-free and
             identical across arms), then the ensemble's best member, its diversity, and the
             RMSD of its coordinate average -- the structure a pipeline built on this objective
             would actually emit.

ONE THING THE CONTACT TERM COSTS, stated before any result.  `E_le1` is separable, so its global
argmin is CERTIFIED with no search.  `E_contact` is a pair term and is not, so **the combined
objective has no certified optimum** and every combined arm is a local-descent result carrying
search error.  The `full` arm has the same limitation, so the comparison is fair, but the
degree-1 arm's certificate is not inherited by the combination and is not claimed for it.
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
CACHE = os.path.join(HERE, "cache")
os.makedirs(CACHE, exist_ok=True)

from s12 import instrument as I                  # noqa: E402
from s14 import avgspace as AV                   # noqa: E402
from s15 import align_lib as A                   # noqa: E402
from s15 import distcal as C                     # noqa: E402
from s15 import seed as SD                       # noqa: E402
from s18 import exp_anova as AN                  # noqa: E402
from s18 import phys_lib as PL                   # noqa: E402

LAM = PL.LAM_LADDER          # (-1, -0.5, 0, +0.5, +1, +2).  NOT extended.
K_POOL = 500                 # the ranker-neutral pool the normalisation is measured on
M_START = 8                  # multi-start count for the generation pass
MAXITER = 400
BAND = 1.5                   # in-band = d <= pool_best + 1.5 A  (s12.instrument.BAND)
NEAR = 2.0


def key(lam):
    return f"lc{int(round(lam * 100)):+04d}"


#: PINNED at driver start.  MATH is landing modules while this runs, and `s18.exp_anova.build`
#: switches to MATH's implementation the instant `s18/math_*.py` exposes the interface.  If that
#: happened mid-run, the first targets would carry the provisional `E_le1` and the rest MATH's --
#: two definitions inside one artefact, which is exactly the failure the brief forbids.  The
#: driver records the source once and `_pin_check` REFUSES to continue if it changes.
_PIN = {"src": None, "hash": None}


def _pin_check():
    src = "math" if AN.load_math() is not None else "s18.exp_anova(provisional)"
    _c, h = AN.config_hash()
    if _PIN["src"] is None:
        _PIN["src"], _PIN["hash"] = src, h
        return
    if (src, h) != (_PIN["src"], _PIN["hash"]):
        raise RuntimeError(
            f"degree-1 definition CHANGED mid-run: {_PIN['src']}/{_PIN['hash']} -> {src}/{h}. "
            f"Refusing to mix two definitions in one artefact; delete the partial and re-run.")


# ------------------------------------------------------------------ the ANOVA object, cached
def build_anova(dhat, sd, i, j, PH, PS, pdb, mu="pool"):
    """`s18.exp_anova.build`, memoised on disk.

    The cache key carries the shared adapter's config hash AND the concrete class's module and
    name, so if MATH lands a different implementation the key changes and everything is rebuilt.
    Only the object's own state is stored; nothing here re-derives the decomposition.
    """
    _cfg, chash = AN.config_hash()
    p = os.path.join(CACHE, f"anova_{pdb}_{mu}_{chash}.npz")
    probe = None
    if os.path.exists(p):
        try:
            z = np.load(p, allow_pickle=True)
            cls = AN.build.__module__          # placeholder; real class recorded below
            meta = json.loads(str(z["_meta"]))
            import importlib
            m = importlib.import_module(meta["mod"])
            cls = getattr(m, meta["cls"])
            ob = object.__new__(cls)
            st = {}
            for k in z.files:
                if k == "_meta":
                    continue
                v = z[k]
                st[k] = v.item() if v.shape == () else v
            ob.__dict__.update(st)
            probe = ob
        except Exception:
            probe = None
    if probe is not None:
        return probe
    ob = AN.build(dhat, sd, i, j, PH, PS, pdb, mu=mu)
    try:
        st = {k: v for k, v in ob.__dict__.items()}
        meta = {"mod": type(ob).__module__, "cls": type(ob).__name__, "cfg": chash}
        np.savez_compressed(p, _meta=json.dumps(meta), **st)
    except Exception:
        pass
    return ob


# ------------------------------------------------------------------ the combined objective
class Combined:
    """`E_le1 + lam_c * s * E_contact`, value and gradient, for L-BFGS."""

    def __init__(self, ob, contact, lam_c, s):
        self.ob, self.ct, self.lam, self.s = ob, contact, float(lam_c), float(s)

    def __call__(self, phi, psi):
        f0, g0 = self.ob.E_le1(phi, psi)
        if self.lam == 0.0:
            return f0, g0
        f1, g1 = self.ct.fg(phi, psi)
        w = self.lam * self.s
        return f0 + w * f1, g0 + w * g1


def combined(base_fun, term, w):
    """`base(theta) + w * term(theta)`, value and gradient.  `w` already carries lam_c AND the
    native-free pool-sd scale, so the caller cannot accidentally apply one without the other."""
    def fun(phi, psi):
        f0, g0 = base_fun(phi, psi)
        if w == 0.0:
            return f0, g0
        f1, g1 = term.fg(phi, psi)
        return f0 + w * f1, g0 + w * g1
    return fun


def lbfgs(fun, phi0, psi0, maxiter=MAXITER):
    from scipy.optimize import minimize
    n = len(phi0)

    def fg(x):
        f, g = fun(x[:n], x[n:])
        return float(f), np.asarray(g, float)

    r = minimize(fg, np.concatenate([np.asarray(phi0, float), np.asarray(psi0, float)]),
                 jac=True, method="L-BFGS-B",
                 options={"maxiter": maxiter, "maxcor": 20, "ftol": 1e-12, "gtol": 1e-10})
    return r.x[:n], r.x[n:], float(r.fun), int(r.nfev), int(r.nit)


# ------------------------------------------------------------------ one target
def run_target(t, data, deb, do_gen=True):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    d = data[pdb]
    i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
    rng = SD.stable_rng(pdb, "s18physlam")

    #: ---- the identical starting candidates: the shipped top-75 and their coordinate average
    W75, PH75, PS75, u = AV.top75_windows(pdb)
    W75 = np.asarray(W75, float)
    PH75 = np.asarray(PH75, float)
    PS75 = np.asarray(PS75, float)
    P = I.pairwise_rmsd(W75)
    avg, _b = I.coordinate_average(W75, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    phi0 = np.asarray(pr["phi"], float)
    psi0 = np.asarray(pr["psi"], float)
    th0 = np.concatenate([phi0, psi0])

    ob = build_anova(dhat, sd, i, j, PH75, PS75, pdb, mu="pool")
    ct = PL.ContactTerm(seq)
    ct_null = PL.ContactTerm(seq, mj=PL.shuffled_mj(seq, SD.stable_rng(pdb, "s18mjnull")))

    #: ---- the ranker-neutral K = 500 pool: the ONLY thing the normalisation is measured on
    ordr = np.asarray(u["order"], int)[:K_POOL]
    Wp = np.asarray(u["W"], float)[ordr]
    PHp = np.asarray(u["PHI"], float)[ordr]
    PSp = np.asarray(u["PSI"], float)[ordr]
    dp = np.asarray(I.kabsch_rmsd_batch(Wp, nat), float)          # ORACLE labels only

    e_full_p = ob._batch(np.concatenate([PHp, PSp], 1))
    e_le1_p = np.array([ob.E_le1(PHp[b], PSp[b])[0] for b in range(len(PHp))])
    e_con_p = np.asarray(ct.of_torsions(PHp, PSp), float)
    e_conN_p = np.asarray(ct_null.of_torsions(PHp, PSp), float)
    sc = PL.pool_scale(e_le1_p, e_con_p)
    scN = PL.pool_scale(e_le1_p, e_conN_p)
    s = sc["s"]

    rec = {"pdb": pdb, "n": n, "fold": fold, "K": int(len(ordr)),
           "avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
           "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
           "scale": sc, "scale_null": scN,
           "pool_best": float(dp.min()), "pool_mean": float(dp.mean()),
           "n_near_pool": int((dp < NEAR).sum()),
           "E0": float(ob.E0), "split_rms": float(getattr(ob, "split_rms", np.nan))}

    #: ---- sensitivity arm: the same scale measured on the shipped TOP-75 instead of the pool.
    e_le1_75 = np.array([ob.E_le1(PH75[b], PS75[b])[0] for b in range(len(PH75))])
    e_con_75 = np.asarray(ct.of_torsions(PH75, PS75), float)
    rec["scale_top75"] = PL.pool_scale(e_le1_75, e_con_75)

    if sc["degenerate"]:
        rec["degenerate"] = True
        return rec

    # =================================================================== PASS 1: ALIGNMENT
    band = dp <= dp.min() + BAND
    al = {}
    series = {"full": e_full_p, "le1": e_le1_p, "contact": e_con_p,
              "contact_null": e_conN_p}
    for lam in LAM:
        series[key(lam)] = e_le1_p + lam * s * e_con_p
    series["null" + key(PL.LAM_PRIMARY)] = e_le1_p + PL.LAM_PRIMARY * scN["s"] * e_conN_p
    rs = rng.standard_normal(len(dp))                     # matched-random score
    series["rand"] = rs
    for nm, v in series.items():
        al[nm] = {
            "rho_global": PL.spearman(v, dp),
            "rho_band": PL.spearman(v[band], dp[band]) if band.sum() >= 4 else float("nan"),
            "sel_global": PL.argmin_tied(v, dp),
            "sel_band": PL.argmin_tied(v[band], dp[band]) if band.sum() >= 1 else float("nan"),
        }
    al["_n_band"] = int(band.sum())
    rec["align"] = al

    # =================================================================== PASS 2: REFINEMENT
    ref = {}

    def _score(nm, ph, ps, fobj, nf, nit):
        ca = I.build_ca(ph, ps)
        ref[nm] = {"rmsd": float(I.ca_rmsd(ca, nat)),
                   "obj_own": float(fobj),
                   "obj_full": float(ob.E_full(ph, ps)[0]),
                   "obj_le1": float(ob.E_le1(ph, ps)[0]),
                   "obj_contact": float(ct.of_torsions(ph, ps)),
                   "disp": float(np.linalg.norm(A.wrap(np.concatenate([ph, ps]) - th0))),
                   "nfev": int(nf), "nit": int(nit)}
        return ref[nm]

    _score("start", phi0, psi0, ob.E_le1(phi0, psi0)[0], 0, 0)

    for lam in LAM:
        ph, ps, f, nf, nit = lbfgs(Combined(ob, ct, lam, s), phi0, psi0)
        _score(key(lam), ph, ps, f, nf, nit)

    ph, ps, f, nf, nit = lbfgs(lambda p, q: ob.E_full(p, q), phi0, psi0)
    _score("full", ph, ps, f, nf, nit)

    #: CONTROL 1 (zero information) -- identical move, MJ table on permuted residue labels.
    ph, ps, f, nf, nit = lbfgs(Combined(ob, ct_null, PL.LAM_PRIMARY, scN["s"]), phi0, psi0)
    _score("null" + key(PL.LAM_PRIMARY), ph, ps, f, nf, nit)

    #: CONTROL 2 (matched random) -- a random torsion move of the PRIMARY arm's own magnitude.
    z = rng.standard_normal(2 * n)
    z *= ref[key(PL.LAM_PRIMARY)]["disp"] / max(float(np.linalg.norm(z)), 1e-12)
    th = th0 + z
    _score("rand_move", th[:n], th[n:], np.nan, 0, 0)

    rec["refine"] = ref

    # =================================================================== PASS 3: GENERATION
    if do_gen:
        gen = {}
        m = min(M_START, len(PH75))
        starts = list(range(m))                    # retrieval order: native-free, identical
        arms = {"le1": Combined(ob, ct, 0.0, s),
                key(PL.LAM_PRIMARY): Combined(ob, ct, PL.LAM_PRIMARY, s),
                "full": (lambda p, q: ob.E_full(p, q)),
                "null" + key(PL.LAM_PRIMARY): Combined(ob, ct_null, PL.LAM_PRIMARY, scN["s"])}
        for nm, fun in arms.items():
            PHo = np.empty((m, n)); PSo = np.empty((m, n))
            for a, b in enumerate(starts):
                p_, q_, _f, _nf, _ni = lbfgs(fun, PH75[b], PS75[b])
                PHo[a], PSo[a] = p_, q_
            Wo = I.build_ca(PHo, PSo)
            do = np.asarray(I.kabsch_rmsd_batch(Wo, nat), float)
            Pw = I.pairwise_rmsd(Wo)
            ca_avg, _bb = I.coordinate_average(Wo, Pw)
            gen[nm] = {"m": int(m),
                       "set_mean": float(do.mean()), "set_best": float(do.min()),
                       "set_median": float(np.median(do)),
                       "n_near": int((do < NEAR).sum()),
                       "diversity": PL.diversity(Wo),
                       "ens_rmsd": float(I.ca_rmsd(np.asarray(ca_avg, float), nat))}
        #: the unrefined starts, as the do-nothing reference for the same m windows
        Ws = W75[:m]
        ds = np.asarray(I.kabsch_rmsd_batch(Ws, nat), float)
        Ps = I.pairwise_rmsd(Ws)
        ca_s, _ = I.coordinate_average(Ws, Ps)
        gen["starts"] = {"m": int(m), "set_mean": float(ds.mean()), "set_best": float(ds.min()),
                         "set_median": float(np.median(ds)), "n_near": int((ds < NEAR).sum()),
                         "diversity": PL.diversity(Ws),
                         "ens_rmsd": float(I.ca_rmsd(np.asarray(ca_s, float), nat))}
        rec["gen"] = gen
    return rec


# ==========================================================================
# PASS `bases` -- DECLARED EXTENSION, added after two things I did not know when
# `PREREG_phys.md` was written.  Both are recorded here rather than folded in silently.
#
#   1. MATH's correction (`s18/MATH_to_EXP.md` section 2): the object the BRIEF's ANOVA bridge
#      actually derives is the RESIDUE-additive `E_res`, not the angle-additive `E_ang` that the
#      provisional adapter calls `E_le1`; and the 2.411 A lattice result belongs to the strict
#      Walsh object, which has NO continuous image.  So "degree-1" is a family, not a thing, and
#      my primary arm's base must be labelled as the coarser member it is.
#   2. The coordinator's alpha ladder: the deployed objective's FUNCTIONAL FORM is sound (alpha=1
#      reaches 1.152 A).  That makes the sharper version of my own hypothesis -- *does
#      `leg_contact` help THE DEPLOYED OBJECTIVE* -- more decision-relevant than helping a
#      truncation of it.
#
# NOTHING ELSE MOVES.  Same contact term, same native-free pool-sd normalisation (recomputed
# against each base's own pool sd), same pre-registered sign, same PRIMARY lam_c = +1.  Only the
# BASE changes, and the ladder is cut to {0, +1} so this cannot become a fishing expedition.
# ==========================================================================
BASE_LAM = (0.0, PL.LAM_PRIMARY)


def run_target_bases(t, data, deb):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    d = data[pdb]
    i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)

    W75, PH75, PS75, u = AV.top75_windows(pdb)
    W75 = np.asarray(W75, float)
    PH75 = np.asarray(PH75, float); PS75 = np.asarray(PS75, float)
    P = I.pairwise_rmsd(W75)
    avg, _b = I.coordinate_average(W75, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    phi0 = np.asarray(pr["phi"], float); psi0 = np.asarray(pr["psi"], float)
    th0 = np.concatenate([phi0, psi0])

    ob = build_anova(dhat, sd, i, j, PH75, PS75, pdb, mu="pool")
    ct = PL.ContactTerm(seq)
    ct_null = PL.ContactTerm(seq, mj=PL.shuffled_mj(seq, SD.stable_rng(pdb, "s18mjnull")))

    ordr = np.asarray(u["order"], int)[:K_POOL]
    PHp = np.asarray(u["PHI"], float)[ordr]; PSp = np.asarray(u["PSI"], float)[ordr]
    e_full_p = ob._batch(np.concatenate([PHp, PSp], 1))
    e_res_p = np.array([ob.E_res(PHp[b], PSp[b])[0] for b in range(len(PHp))])
    e_con_p = np.asarray(ct.of_torsions(PHp, PSp), float)
    e_conN_p = np.asarray(ct_null.of_torsions(PHp, PSp), float)
    sc = {"full": PL.pool_scale(e_full_p, e_con_p), "res": PL.pool_scale(e_res_p, e_con_p)}
    scN = {"full": PL.pool_scale(e_full_p, e_conN_p), "res": PL.pool_scale(e_res_p, e_conN_p)}

    rec = {"pdb": pdb, "n": n, "fold": fold, "scale": sc, "scale_null": scN,
           "avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
           "proj": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)), "arms": {}}
    if sc["full"]["degenerate"] or sc["res"]["degenerate"]:
        rec["degenerate"] = True
        return rec

    def score(nm, fun):
        ph, ps, f, nf, _ni = lbfgs(fun, phi0, psi0)
        rec["arms"][nm] = {
            "rmsd": float(I.ca_rmsd(I.build_ca(ph, ps), nat)),
            "obj_own": float(f), "obj_full": float(ob.E_full(ph, ps)[0]),
            "obj_res": float(ob.E_res(ph, ps)[0]),
            "obj_contact": float(ct.of_torsions(ph, ps)),
            "disp": float(np.linalg.norm(A.wrap(np.concatenate([ph, ps]) - th0))),
            "nfev": int(nf)}

    bases = {"full": (lambda p, q: ob.E_full(p, q)), "res": (lambda p, q: ob.E_res(p, q))}
    for bn, bf in bases.items():
        for lam in BASE_LAM:
            score(f"{bn}_{key(lam)}", combined(bf, ct, lam * sc[bn]["s"]))
        #: CONTROL -- identical move, MJ table on permuted residue labels.
        score(f"{bn}_null{key(PL.LAM_PRIMARY)}",
              combined(bf, ct_null, PL.LAM_PRIMARY * scN[bn]["s"]))
    return rec


# ------------------------------------------------------------------ driver
def run(targets=None, out="lam.json", do_gen=True, verbose=True, bases=False):
    tg = targets if targets is not None else I.targets()
    _cfg, chash = AN.config_hash()
    cfg = {"LAM": list(LAM), "K_POOL": K_POOL, "M_START": M_START, "MAXITER": MAXITER,
           "BAND": BAND, "anova_cfg": _cfg, "anova_hash": chash,
           "anova_source": "math" if AN.load_math() is not None else "s18.exp_anova(provisional)",
           "normalisation": "s = sd_pool(base)/sd_pool(E_contact) over the K=500 retrieval pool",
           "lam_primary": PL.LAM_PRIMARY, "bases_pass": bool(bases),
           "BASE_LAM": list(BASE_LAM) if bases else None}
    path = os.path.join(PL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            prev = json.load(open(path))
            if prev.get("cfg_hash") == PL.cfg_hash(cfg):
                rows = prev.get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    deb = {}
    for f in sorted({int(t["fold"]) for t in tg}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    t0 = time.time()
    _pin_check()
    for t in tg:
        if t["pdb"] in done:
            continue
        _pin_check()
        PL.mem_ok(1.2)
        rows.append(run_target_bases(t, data, deb) if bases
                    else run_target(t, data, deb, do_gen=do_gen))
        if verbose:
            r = rows[-1]
            if bases:
                a = r.get("arms", {})
                print(f"  {len(rows)}/{len(tg)} {t['pdb']} start={r.get('proj', float('nan')):.3f} "
                      f"full={a.get('full_lc+000', {}).get('rmsd', float('nan')):.3f} "
                      f"full+c={a.get('full_lc+100', {}).get('rmsd', float('nan')):.3f} "
                      f"res={a.get('res_lc+000', {}).get('rmsd', float('nan')):.3f} "
                      f"res+c={a.get('res_lc+100', {}).get('rmsd', float('nan')):.3f} "
                      f"({time.time()-t0:.0f}s)", flush=True)
            else:
                rf = r.get("refine", {})
                print(f"  {len(rows)}/{len(tg)} {t['pdb']} start={r.get('proj', float('nan')):.3f} "
                      f"le1={rf.get('lc+000', {}).get('rmsd', float('nan')):.3f} "
                      f"lc+1={rf.get('lc+100', {}).get('rmsd', float('nan')):.3f} "
                      f"full={rf.get('full', {}).get('rmsd', float('nan')):.3f} "
                      f"({time.time()-t0:.0f}s)", flush=True)
        PL.write(out, {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg)},
                 n_expected=len(tg))
    PL.write(out, {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg)}, n_expected=len(tg))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-gen", action="store_true")
    ap.add_argument("--bases", action="store_true",
                    help="the DECLARED EXTENSION: E_full and E_res as the base instead of E_ang")
    ap.add_argument("--out", default="lam.json")
    a = ap.parse_args()
    tg = I.targets()
    if a.limit:
        tg = tg[:a.limit]
    run(tg, out=a.out, do_gen=not a.no_gen, bases=a.bases)
