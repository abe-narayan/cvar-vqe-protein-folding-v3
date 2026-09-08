"""s20/c_land_null.py -- THE MATCHED CONTROL FOR Q2's RMSD RESULT.

DECLARED EXTENSION, written while `s20.c_land` was running and BEFORE its table was read, after
its first three per-target lines showed the same pattern three times: minimising **Legacy** in
torsion space moves Ca-RMSD TOWARD the native and minimising **AMBER** moves it AWAY.  A result
of that shape is worthless without a control matched **in the space the operator actually works
in** (`BRIEF` section 6 rule 1, the programme's most repeated error), and `c_land` does not carry
one.  So it is built here rather than asserted there.

THE OPERATOR.  Each minimisation is a move in continuous torsion space of a measured size
`theta_moved` (per-coordinate RMS torus distance, radians) from a start that is a pool member.
The matched control must therefore be **a torsion-space move of the same torus magnitude from the
same start**, and nothing else may differ.

TWO NULLS, because one of them is the wrong kind and this project has been burned by exactly that
(`BRIEF` section 6 rule 4: *uniform-on-the-torus places mass on impossible backbone conformations,
making it a WORSE measure rather than an uninformative one*):

    rand_iso        isotropic random direction on the active coordinates, scaled to the SAME
                    torus magnitude.  This is the STRICT magnitude-matched null and it is
                    expected to be terrible -- it is reported so the scale is visible, not as
                    the null the conclusion rests on.
    toward_member   move the same torus distance ALONG THE GEODESIC toward another randomly
                    chosen member of the same pool.  Realisable, native-free, zero information
                    about the target, and it lands on real peptide conformations.  **This is the
                    null the conclusion rests on.**

FALSIFIER, registered here before the run.  The Q2 claim "minimising Legacy moves toward the
native / minimising AMBER moves away" is **REFUTED as a property of the POTENTIAL** if the
`toward_member` null of the same magnitude reproduces either sign with a CI that does not exclude
the observed effect -- i.e. if any displacement of that size does the same thing, the potential
is not what did it.

No energy is evaluated here at all: this module needs only `build_ca` and the frozen metric.

    python -m s20.c_land_null
"""
from __future__ import annotations

import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
from s12 import instrument as I
from s18 import phys_lib as PL                     # noqa: E402
from s14.avgspace import top75_windows              # noqa: E402
from s15 import seed as SD                          # noqa: E402

N_DRAWS = 30


# ---------------------------------------------------------------------------------------
# A SELF-CAUGHT LABELLING DEFECT.  `s12.instrument.paired`'s `ci95` is a PLAIN i.i.d.
# target-level bootstrap; its `folds` argument only adds a per-fold mean breakdown and does
# NOT cluster the resample.  Sprint 19's numbers were made with `s18.phys_lib.paired`, which
# returns BOTH an i.i.d. `ci` and a FOLD-CLUSTERED `ci_fold`, and "fold-aware" in this
# programme means the latter (BRIEF section 8).  Every CI below is therefore produced by
# `PL.paired` and the FOLD-CLUSTERED interval is the one quoted, with the i.i.d. one printed
# beside it wherever both matter.  Caught by reading `s12/instrument.py:188` rather than
# trusting the parameter name.
def PP(a, b, folds=None):
    st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
    return {"mean_diff": st["mean"], "ci95": st.get("ci_fold", st["ci"]),
            "ci_iid": st["ci"], "median_diff": st["median"],
            "n_better": st["W"], "n_worse": st["L"],
            "mean_a": st["mean_a"], "mean_b": st["mean_b"], "n": st["n"]}


def torus_d(a, b):
    d = np.mod(np.asarray(a, float) - np.asarray(b, float) + np.pi, 2 * np.pi) - np.pi
    return float(np.sqrt((d ** 2).mean()))


def run(src="c_land.json", out="c_land_null.json"):
    o = json.load(open(os.path.join(RESULTS, src)))
    rows_in = o["rows"]
    tg = {t["pdb"]: t for t in I.targets()}
    N_STARTS = o["config"]["N_STARTS"]
    rows = []
    for r in rows_in:
        pdb, n = r["pdb"], int(r["n"])
        W, PHI, PSI, u = top75_windows(pdb)
        W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
        nat = np.asarray(u["nat_ca"], float)
        #: THE IDENTICAL STARTS -- reconstructed by the same deterministic rule c_land used.
        Pm = I.pairwise_rmsd(W)
        med = int(I.medoid(Pm))
        rng0 = SD.stable_rng(pdb, "s20C_land")
        others = [int(x) for x in rng0.permutation(len(W))[:N_STARTS] if int(x) != med]
        order = [med] + others[:N_STARTS - 1]
        starts = [np.concatenate([PHI[b], PSI[b]]) for b in order]
        active = np.arange(1, 2 * n)          # phi[0] moves no atom (EXACT)

        rng = SD.stable_rng(pdb, "s20C_landnull")
        rec = {"pdb": pdb, "n": n, "fold": int(r["fold"]), "null": {}}
        for pot in ("legacy", "amber"):
            mag = float(r["multi"][pot]["theta_moved"])
            rs = float(r["multi"][pot]["rmsd_start"])
            re_ = float(r["multi"][pot]["rmsd_end"])
            iso, twd = [], []
            for th_s in starts:
                for _ in range(N_DRAWS):
                    #: isotropic, matched torus magnitude
                    v = np.zeros_like(th_s)
                    g = rng.standard_normal(len(active))
                    v[active] = g / np.sqrt((g ** 2).mean() * len(th_s) / len(active))
                    th = th_s + mag * v
                    iso.append(float(I.ca_rmsd(I.build_ca(th[:n], th[n:]), nat)))
                    #: geodesic toward another pool member, same torus distance
                    j = int(rng.integers(len(W)))
                    tj = np.concatenate([PHI[j], PSI[j]])
                    d = np.mod(tj - th_s + np.pi, 2 * np.pi) - np.pi
                    dn = float(np.sqrt((d[active] ** 2).mean()))
                    if dn < 1e-9:
                        continue
                    th2 = th_s.copy()
                    th2[active] = th_s[active] + d[active] * (mag / dn)
                    twd.append(float(I.ca_rmsd(I.build_ca(th2[:n], th2[n:]), nat)))
            rec["null"][pot] = {
                "theta_moved": mag,
                "rmsd_start": rs, "rmsd_end": re_,
                "d_real": re_ - rs,
                "rand_iso": float(np.mean(iso)), "d_rand_iso": float(np.mean(iso)) - rs,
                "toward_member": float(np.mean(twd)),
                "d_toward_member": float(np.mean(twd)) - rs,
                "n_draws": int(len(iso))}
        rows.append(rec)
    obj = {"rows": rows, "config": {"N_DRAWS": N_DRAWS, "src": src},
           "n_rows": len(rows), "n_expected": o.get("n_expected", len(rows)),
           "complete": len(rows) >= int(o.get("n_expected", len(rows)))}
    with open(os.path.join(RESULTS, out), "w") as fh:
        json.dump(obj, fh)
    report(obj)
    return obj


def report(o=None, out="c_land_null.json"):
    if o is None:
        o = json.load(open(os.path.join(RESULTS, out)))
    rows = o["rows"]
    n = len(rows)
    tag = "" if o.get("complete") else f"  *** PARTIAL n={n} ***"
    folds = np.array([r["fold"] for r in rows], int)
    print("=" * 100)
    print(f"Q2b  THE MATCHED CONTROL FOR THE TORSION-SPACE MINIMISATION   n = {n} targets{tag}")
    print("     Every arm is a move of the SAME per-coordinate RMS torus magnitude from the SAME")
    print("     starts.  Negative d = moved TOWARD the native.  No energy is evaluated here.")
    print("=" * 100)
    print(f"\n{'potential':<10s}{'|dtheta| (rad)':>15s}{'start':>9s}{'minimised':>11s}"
          f"{'toward_member':>15s}{'rand_iso':>10s}")
    for pot in ("legacy", "amber"):
        f = lambda k: np.array([r["null"][pot][k] for r in rows], float)   # noqa: E731
        print(f"{pot:<10s}{f('theta_moved').mean():15.4f}{f('rmsd_start').mean():9.4f}"
              f"{f('rmsd_end').mean():11.4f}{f('toward_member').mean():15.4f}"
              f"{f('rand_iso').mean():10.4f}")
    print(f"\n{'comparison':<40s}{'mean':>10s}{'CI95':>24s}{'median':>9s}{'W/L':>10s}")
    for pot in ("legacy", "amber"):
        f = lambda k: np.array([r["null"][pot][k] for r in rows], float)   # noqa: E731
        for lab, a, b in (
                (f"{pot}: minimised - start", f("rmsd_end"), f("rmsd_start")),
                (f"{pot}: minimised - toward_member", f("rmsd_end"), f("toward_member")),
                (f"{pot}: minimised - rand_iso", f("rmsd_end"), f("rand_iso")),
                (f"{pot}: toward_member - start", f("toward_member"), f("rmsd_start"))):
            st = PP(a, b, folds=folds)
            star = "  *" if (st["ci95"][0] > 0 or st["ci95"][1] < 0) else ""
            print(f"{lab:<40s}{st['mean_diff']:+10.4f}  [{st['ci95'][0]:+9.4f},"
                  f"{st['ci95'][1]:+9.4f}]{st['median_diff']:+9.4f}"
                  f"{st['n_better']:>5d}/{st['n_worse']:<4d}{star}")
    fl = np.array([r["null"]["legacy"]["rmsd_end"] for r in rows], float)
    fa = np.array([r["null"]["amber"]["rmsd_end"] for r in rows], float)
    st = PP(fl, fa, folds=folds)
    print(f"\n    Legacy-minimised MINUS AMBER-minimised (identical starts, optimiser, budget; "
          f"only H changes): {st['mean_diff']:+.4f} [{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}] "
          f"med {st['median_diff']:+.4f}  {st['n_better']}W/{st['n_worse']}L")
    print("=" * 100)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--report":
        report()
    else:
        run()
