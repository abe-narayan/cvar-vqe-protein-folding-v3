"""SPRINT 18 / QUANTUM-ADVERSARIAL -- attacks 3 and 4 on the coordinator's `objceil` headline
(`s18/COORD_FINDING.md`).  Attacks 1 and 2 (the separation-stratified shuffle and the
weight-residual pairing) were taken by the coordinator and are NOT duplicated here.

ATTACK 3 -- IS alpha = 1's 1.152 A MEASURING THE OBJECTIVE, OR THE OPTIMISER'S BASIN?

With perfect distances the native is (near) a global optimum of the functional, so a refinement
that lands at 1.152 A is either
  (i)  OBJECTIVE-LIMITED  -- the functional's optimum genuinely is 1.152 A away, because the
       ideal-geometry torsion parameterisation cannot represent the native exactly, or
  (ii) CONVERGENCE-LIMITED -- the optimum is at ~0 and the fit stopped in a basin reachable
       from that particular start.
These have different consequences: (i) makes 1.152 a real ceiling, (ii) makes it a property of
one start and the alpha ladder cannot be read as a requirement curve.

They are separated by THREE measurements, none of which needs a new objective:
  * THE PARAMETERISATION FLOOR.  Project the NATIVE itself into the same ideal-geometry torsion
    parameterisation and measure its RMSD.  Nothing in this parameterisation can beat that, so
    it is a hard lower bound on alpha = 1 that has nothing to do with the optimiser.
  * MULTI-START.  Refine at alpha = 1 from five different starts -- the coordinate average's
    projection (the coordinator's start), a random pool member, a constant ideal alpha-helix
    (the ledger's zero-information start), a random torsion draw, and the NATIVE'S OWN
    PROJECTED TORSIONS (an ORACLE start, the best case).  Spread across starts = basin
    dependence; agreement = a real ceiling.
  * THE OBJECTIVE AT EACH ENDPOINT, printed beside the RMSD, and the objective evaluated AT the
    native's projected torsions.  If the fit's objective is far above the native's, the fit did
    not converge; if it is at or below it, the optimum genuinely is not the native.

ATTACK 4 -- THE REQUIREMENT CURVE ASSUMES THE IMPROVEMENT IS UNIFORM ACROSS PAIRS.

A note first, because it changes what needs testing.  The coordinator's own alpha ladder,

    d_alpha = (1-alpha) dhat + alpha d_true  =  d_true + (1-alpha) * (dhat - d_true)

is EXACTLY a uniform rescaling of the residual vector by (1-alpha).  It therefore already
preserves the residual's direction and its entire correlation structure -- the coordinator's
stated attack-4 concern ("shrink the residual while preserving its correlation structure") is
what the ladder already does.  The live objection is different and sharper:

    a real predictor does NOT improve uniformly across pairs.  It gets better somewhere.

So this builds a family of error models ALL AT THE SAME RESIDUAL RMS as alpha = 0.5 -- the rung
the coordinator quotes as "halve the residual, reach 2.5 A" -- differing only in WHERE the
improvement lands: short-range vs long-range pairs, largest vs smallest residuals, most vs
least confident pairs (the objective weights by 1/sd^2), and a pure de-biasing that removes the
per-separation systematic component.  If the family's outcomes span a wide range, "halve the
residual RMS -> 2.5 A" is one point in a wide band and must not be quoted as a requirement.

EVERY ARM HERE IS AN ORACLE DIAGNOSTIC.  Each builds its optimisation target from native
distances.  Nothing is deployable and nothing may be quoted as a result.

RUN:  python -m s18.q_ceilattack run  [N]
      python -m s18.q_ceilattack report
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
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as AL              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

SALT = "s18qceil"
OUT = os.path.join(RESULTS, "q_ceilattack.json")


def _rmsd_of(phi, psi, nat):
    return float(I.ca_rmsd(I.build_ca(phi, psi), nat))


def _obj(d, phi, psi, dhat, sd, i, j):
    """The functional's own value: sum_p (d_p(theta) - dhat_p)^2 / sd_p^2, at the SAME weights
    the fit used.  Printed beside every endpoint so convergence is visible, not inferred."""
    CA = np.asarray(I.build_ca(phi, psi), float)
    dd = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
    return float((((dd - dhat) / sd) ** 2).sum())


def _matched(r, keep, target_rms):
    """A residual vector that keeps `r` on the pairs where `keep` is False and shrinks it on
    the pairs where `keep` is True, scaled so the whole vector's RMS equals `target_rms`.

    `keep=True` means "the better predictor fixed this pair".  A single global scale `s` is
    then solved for so the RMS matches exactly, which is what makes the family COMPARABLE:
    every member has the identical residual RMS and differs only in WHERE its error sits.
    """
    base = r.copy()
    base[keep] = 0.0
    cur = float(np.sqrt((base ** 2).mean()))
    if cur < 1e-12:
        return base
    return base * (target_rms / cur)


def cell(t, data, deb):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    d = data[pdb]
    i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
    dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
    sep = np.abs(np.asarray(j) - np.asarray(i))
    r = dhat - dtrue
    rms = float(np.sqrt((r ** 2).mean()))
    rng = SD.stable_rng(pdb, "ceilattack", salt=SALT)

    W = np.asarray(AV.top75_windows(pdb)[0], float)
    P = I.pairwise_rmsd(W)
    avg, _b = I.coordinate_average(W, P)
    pr = I.project(np.asarray(avg, float), seq, fold)
    phi0, psi0 = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)

    e = {"pdb": pdb, "n": n, "fold": fold,
         "avg": float(I.ca_rmsd(np.asarray(avg, float), nat)),
         "proj_avg": _rmsd_of(phi0, psi0, nat),
         "resid_rms": rms}

    # ------------------------------------------------ ATTACK 3a: the parameterisation floor
    prn = I.project(nat, seq, fold)                          # ORACLE DIAGNOSTIC
    phin, psin = np.asarray(prn["phi"], float), np.asarray(prn["psi"], float)
    e["param_floor_native_projected"] = _rmsd_of(phin, psin, nat)
    e["obj_at_native_projection"] = _obj(None, phin, psin, dtrue, sd, i, j)

    # ------------------------------------------------ ATTACK 3b: multi-start at alpha = 1
    ideal_phi = np.full(n, np.deg2rad(-57.0)); ideal_psi = np.full(n, np.deg2rad(-47.0))
    kmember = int(rng.integers(0, W.shape[0]))
    prm = I.project(np.asarray(W[kmember], float), seq, fold)
    starts = {
        "avgproj": (phi0, psi0),                              # the coordinator's start
        "member": (np.asarray(prm["phi"], float), np.asarray(prm["psi"], float)),
        "helix": (ideal_phi, ideal_psi),                      # ZERO-INFORMATION start
        "randtors": (rng.uniform(-np.pi, np.pi, n), rng.uniform(-np.pi, np.pi, n)),
        "nativeproj_ORACLE": (phin, psin),                    # the best case
    }
    e["a1_multistart"] = {}
    for nm, (p0, q0) in starts.items():
        p_, q_, f_ = AL.fit(dtrue, sd, i, j, p0, q0)
        e["a1_multistart"][nm] = {
            "rmsd": _rmsd_of(p_, q_, nat),
            "obj": _obj(None, p_, q_, dtrue, sd, i, j),
            "start_rmsd": _rmsd_of(p0, q0, nat),
            "f": float(f_) if np.isscalar(f_) or np.ndim(f_) == 0 else float(np.min(f_)),
        }

    # ------------------------------------------------ ATTACK 4: matched-RMS error families
    #  every member has residual RMS = 0.5 * rms, the rung the coordinator quotes as 2.5 A
    tgt = 0.5 * rms
    fams = {}
    fams["unif"] = r * 0.5                                     # == the coordinator's alpha=0.5
    med_sep = 5
    fams["fix_short"] = _matched(r, sep <= med_sep, tgt)
    fams["fix_long"] = _matched(r, sep > med_sep, tgt)
    big = np.abs(r) >= np.median(np.abs(r))
    fams["fix_big"] = _matched(r, big, tgt)
    fams["fix_small"] = _matched(r, ~big, tgt)
    conf = np.asarray(sd, float) <= np.median(np.asarray(sd, float))
    fams["fix_confident"] = _matched(r, conf, tgt)
    fams["fix_unconfident"] = _matched(r, ~conf, tgt)
    # de-bias only: remove the per-separation-bin mean residual, then match the RMS
    rb = r.copy()
    for lo, hi in ((2, 3), (4, 5), (6, 7), (8, 10), (11, 999)):
        m = (sep >= lo) & (sep <= hi)
        if m.sum() > 0:
            rb[m] = r[m] - r[m].mean()
    cur = float(np.sqrt((rb ** 2).mean()))
    fams["debias_sep"] = rb * (tgt / cur) if cur > 1e-12 else rb

    e["family"] = {}
    for nm, rr in fams.items():
        dd = np.maximum(dtrue + rr, 2.0)
        p_, q_, _f = AL.fit(dd, sd, i, j, phi0, psi0)
        e["family"][nm] = {"rmsd": _rmsd_of(p_, q_, nat),
                           "resid_rms": float(np.sqrt((rr ** 2).mean())),
                           "resid_rms_after_floor": float(
                               np.sqrt(((dd - dtrue) ** 2).mean()))}
    return e


def run(limit=None):
    tg = I.targets()
    if limit:
        tg = tg[:int(limit)]
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    folds = sorted({int(t["fold"]) for t in tg})
    deb = {}
    for f in folds:
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))
    rows = []
    if os.path.exists(OUT):
        try:
            rows = json.load(open(OUT))["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for c, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(cell(t, data, deb))
        if (len(rows)) % 5 == 0:
            json.dump({"rows": rows, "complete": False}, open(OUT, "w"), default=float)
            print(f"  {len(rows)}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    report(rows)


def _boot(dif, rng, B=4000):
    dif = np.asarray(dif, float); k = len(dif)
    m = dif[rng.integers(0, k, size=(B, k))].mean(1)
    return float(dif.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("ceilattack", "report", salt=SALT)
    P = print
    g = lambda k: np.array([r[k] for r in rows])                      # noqa: E731
    P("=" * 100)
    P(f"ATTACK ON `objceil`   n = {len(rows)}   EVERY ARM IS AN ORACLE DIAGNOSTIC")
    P("=" * 100)
    P(f"  coordinate average (start)                 {g('avg').mean():.3f}")
    P(f"  its ideal-geometry projection              {g('proj_avg').mean():.3f}")
    P(f"  distogram residual RMS                     {g('resid_rms').mean():.3f} A")
    P("")
    P("ATTACK 3a -- THE PARAMETERISATION FLOOR (nothing in this torsion space can beat it)")
    P("-" * 100)
    pf = g("param_floor_native_projected")
    P(f"  the NATIVE projected into the same ideal-geometry torsion parameterisation:")
    P(f"     mean {pf.mean():.3f}   median {np.median(pf):.3f}   "
      f"sd {pf.std(ddof=1):.3f}   max {pf.max():.3f}")
    P("")
    P("ATTACK 3b -- MULTI-START AT alpha = 1 (perfect distances)")
    P("-" * 100)
    starts = list(rows[0]["a1_multistart"].keys())
    P(f"{'start':<20s}{'start RMSD':>12s}{'end RMSD':>11s}{'median':>9s}"
      f"{'objective at end':>19s}{'W/L vs avgproj':>17s}")
    base = np.array([r["a1_multistart"]["avgproj"]["rmsd"] for r in rows])
    for s in starts:
        v = np.array([r["a1_multistart"][s]["rmsd"] for r in rows])
        s0 = np.array([r["a1_multistart"][s]["start_rmsd"] for r in rows])
        ob = np.array([r["a1_multistart"][s]["obj"] for r in rows])
        w = int((v < base - 1e-9).sum()); l = int((v > base + 1e-9).sum())
        P(f"{s:<20s}{s0.mean():>12.3f}{v.mean():>11.3f}{np.median(v):>9.3f}"
          f"{ob.mean():>19.2f}{f'{w}/{l}':>17s}")
    on = g("obj_at_native_projection")
    P(f"{'(objective AT the native projection)':<20s}{'':>12s}{pf.mean():>11.3f}"
      f"{np.median(pf):>9.3f}{on.mean():>19.2f}")
    P("")
    for s in starts:
        if s == "avgproj":
            continue
        v = np.array([r["a1_multistart"][s]["rmsd"] for r in rows])
        m, lo, hi = _boot(v - base, rng)
        P(f"  {s:<20s} - avgproj : {m:+.3f} [{lo:+.3f}, {hi:+.3f}]  median "
          f"{np.median(v-base):+.3f}")
    P("")
    P("  READ: if every start lands at the same RMSD, alpha = 1 is OBJECTIVE-LIMITED and the")
    P("  ladder is a real requirement curve.  If the ORACLE native start lands far lower with")
    P("  a LOWER objective, the fit is CONVERGENCE-LIMITED and 1.152 A is a basin, not a")
    P("  ceiling.  The objective column decides it, not the RMSD column.")
    P("")
    P("ATTACK 4 -- MATCHED-RESIDUAL-RMS ERROR FAMILY (all at 0.5x the deployed residual RMS,")
    P("            i.e. the rung quoted as 'halve the residual -> 2.5 A')")
    P("-" * 100)
    fams = list(rows[0]["family"].keys())
    P(f"{'error model':<20s}{'resid RMS':>11s}{'RMSD':>9s}{'median':>9s}"
      f"{'vs unif':>26s}{'W/L':>10s}")
    bu = np.array([r["family"]["unif"]["rmsd"] for r in rows])
    for f in fams:
        v = np.array([r["family"][f]["rmsd"] for r in rows])
        rr = np.array([r["family"][f]["resid_rms"] for r in rows])
        m, lo, hi = _boot(v - bu, rng)
        w = int((v < bu - 1e-9).sum()); l = int((v > bu + 1e-9).sum())
        P(f"{f:<20s}{rr.mean():>11.3f}{v.mean():>9.3f}{np.median(v):>9.3f}"
          f"{f'{m:+.3f} [{lo:+.3f},{hi:+.3f}]':>26s}{f'{w}/{l}':>10s}")
    allv = np.array([[r["family"][f]["rmsd"] for f in fams] for r in rows]).mean(0)
    P("")
    P(f"  SPREAD across error models at the IDENTICAL residual RMS: "
      f"{allv.min():.3f} to {allv.max():.3f} A  (range {allv.max()-allv.min():.3f} A)")
    P("  READ: the coordinator's requirement curve reads the `unif` row only.  If the range")
    P("  above straddles 2.5 A, 'halve the residual RMS' does NOT determine the outcome and")
    P("  the requirement must be quoted with this band.")
    P("")


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "run"
    if m == "run":
        run(sys.argv[2] if len(sys.argv) > 2 else None)
    else:
        report()
