"""SPRINT 15, coordinator -- SURROGATE DESTRUCTION on the distance channel.

WHAT THE PHASE DIAGRAM JUST EXPOSED, AND WHY IT NEEDS EXPLAINING.

`s15/distacc.py` corrupts TRUE distances with i.i.d. Gaussian noise of a known size and fits the
result. The curve, on 40 targets:

    sigma   0.00  0.25  0.50  0.75  1.00  1.50  2.00  3.00
    RMSD    0.554 0.709 0.939 1.306 1.450 1.787 2.075 2.561

The distogram's measured MAE is 2.386 A, which for a Gaussian corresponds to sigma ~ 2.99. The phase
diagram says i.i.d. noise of that size emits about **2.56 A**. The real distogram emits **3.644 A**.

**The real error costs about 1.08 A MORE than i.i.d. noise of the same magnitude.** Something about
its SHAPE, not its size, is expensive.

That is the exact mirror of what the INFO workstream found for the TORSION channel, where a real
64-degree-RMS channel emits 1.77 A while i.i.d. noise of the same magnitude emits 4.71 A -- there,
real error is 2.9 A CHEAPER than i.i.d. Two channels in the same pipeline, both with structured
errors, and the structure helps one and hurts the other. That contrast is worth a paper section on
its own, but only if the mechanism is measured rather than asserted.

THE METHOD: SURROGATE DESTRUCTION. Take the ACTUAL error vector `e = dhat - dtrue`, destroy exactly
one property of it, keep everything else, refit, and measure the change. Each arm isolates one
structural feature, and the arms are nested so their effects can be read off in sequence.

  ORACLE LABELLING, WHICH MATTERS MORE HERE THAN ANYWHERE ELSE IN THE SPRINT. Every surrogate arm
  needs the TRUE error vector `e = dhat - dtrue` in order to destroy a property of it, so every arm
  below except `real` and `debias_sep` is an **ORACLE DIAGNOSTIC**. They measure what the error's
  structure COSTS. None of them is a method, and none may ever appear in a predictive table. If an
  arm shows a large gain, the finding is "this property is expensive", and the follow-up question --
  a separate experiment -- is whether that property can be estimated or removed native-free.

    real              e as measured                    NATIVE-FREE                the predicted arm
    debias_sep        e minus its separation mean, leave-fold-out      NATIVE-FREE. Removes the
                                                                       systematic expansion
    shuffle_pairs     e permuted across pairs within the target  ORACLE. Destroys ALL separation
                                                                       structure, preserves the
                                                                       marginal distribution exactly
    shuffle_signs     |e| kept exactly, signs randomised         ORACLE. Destroys bias and any
                                                                       directional alignment,
                                                                       preserves every magnitude
    winsor_z3         |e| clipped at 3 sd                        ORACLE. Removes the outlier tail
                                                                       (max z = 7.237) and nothing
                                                                       else
    gauss_matched     i.i.d. Gaussian at the target's error RMS  ORACLE. The i.i.d. surface,
                                                                       per target rather than
                                                                       globally
    ORACLE_true       e = 0                                      ORACLE. The ceiling, 0.611 A

READ THE ARMS AS A DECOMPOSITION. If `shuffle_pairs` is much better than `real`, the separation
structure is what costs. If `winsor_z3` is much better, the outliers are. If `shuffle_signs` is
much better, the error has a coherent direction that is being followed. If `gauss_matched` lands on
`real`, there is no structure cost at all and the phase-diagram gap was a mismatch between MAE and
sigma rather than a property of the errors -- which is the null this experiment must be able to
return, and the most likely single explanation for a gap of this size.

**That null is the honest first hypothesis.** MAE-to-sigma conversion assumes Gaussianity, and the
distogram's errors are plainly not Gaussian. `gauss_matched` matches the RMS per target and so
removes that conversion from the comparison entirely. If the gap survives `gauss_matched`, it is
real; if it vanishes, the 1.08 A was an artefact of my own arithmetic and this module says so.

Run:
    python -m s15.errstruct
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

from s12 import instrument as I              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
os.makedirs(RESULTS, exist_ok=True)

#: only `real` and `debias_sep` are native-free; the rest need the true error vector and are
#: ORACLE DIAGNOSTICS, renamed to carry the label into every table row.
ARMS = ("real", "debias_sep", "ORACLE_shuffle_pairs", "ORACLE_shuffle_signs",
        "ORACLE_winsor_z3", "ORACLE_gauss_matched", "ORACLE_true")
N_REP = 3          # the shuffles are stochastic; average over repeats


def surrogate(e, sep, sd, arm, rng, bias_fn=None):
    """Destroy exactly one property of the measured error vector `e = dhat - dtrue`."""
    if arm == "real":
        return e
    if arm == "ORACLE_true":
        return np.zeros_like(e)
    if arm == "debias_sep":
        return e - bias_fn(sep)
    if arm == "ORACLE_shuffle_pairs":
        return e[rng.permutation(len(e))]
    if arm == "ORACLE_shuffle_signs":
        return np.abs(e) * rng.choice([-1.0, 1.0], size=len(e))
    if arm == "ORACLE_winsor_z3":
        z = e / sd
        return np.clip(z, -3.0, 3.0) * sd
    if arm == "ORACLE_gauss_matched":
        return rng.normal(0.0, float(np.sqrt((e * e).mean())), size=len(e))
    raise ValueError(arm)


def run(targets=None, n_start=4):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    fail = np.isin(pdbs, I.FAIL18)

    data = C.gather(tg)
    #: the separation-bias profile is fitted leave-fold-out, exactly as everywhere else, so the
    #: `debias_sep` arm stays native-free and is directly comparable to the distcal arm.
    bias = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        bias[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    res = {a: [] for a in ARMS}
    stats = []
    path = os.path.join(RESULTS, "errstruct.json")

    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, sep, sd = d["i"], d["j"], d["sep"], d["sd"]
        dtrue, dhat = d["dtrue"], d["dhat"]
        e = dhat - dtrue
        w = 1.0 / sd ** 2
        starts = D.starts(p, d["seq"], d["n"], d["fold"], n_start,
                          SD.stable_rng(p, "errstruct"))
        stats.append({"pdb": p, "rms_e": float(np.sqrt((e * e).mean())),
                      "mae_e": float(np.abs(e).mean()), "bias_e": float(e.mean()),
                      "max_z": float(np.abs(e / sd).max())})
        for a in ARMS:
            reps = 1 if a in ("real", "ORACLE_true", "debias_sep",
                              "ORACLE_winsor_z3") else N_REP
            vals = []
            for r in range(reps):
                rng = SD.stable_rng(p, a, r)
                tgt = np.maximum(dtrue + surrogate(e, sep, sd, a, rng, bias[d["fold"]]), 2.0)
                best = None
                for phi0, psi0, _tag in starts:
                    phi, psi, f, _ = D.fit_distances(tgt, w, i, j, phi0, psi0)
                    if best is None or f < best[0]:
                        best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), d["nat"])))
                vals.append(best[1])
            res[a].append(float(np.mean(vals)))
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res["real"], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "arms": {},
           "error_stats": {"rms_mean": float(np.mean([s["rms_e"] for s in stats])),
                           "mae_mean": float(np.mean([s["mae_e"] for s in stats])),
                           "bias_mean": float(np.mean([s["bias_e"] for s in stats])),
                           "max_z_mean": float(np.mean([s["max_z"] for s in stats]))},
           "per_target": {a: dict(zip(pdbs, map(float, res[a]))) for a in ARMS}}
    for a in ARMS:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "FAIL18": float(v[fail].mean()),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_real": I.paired(v, base, folds=folds, names=pdbs),
                          "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_errstruct", out, n_expected=len(tg))

    es = out["error_stats"]
    print(f"\nn = {len(pdbs)}   incumbent {ref.mean():.3f}")
    print(f"the measured distogram error: RMS {es['rms_mean']:.3f} A, MAE "
          f"{es['mae_mean']:.3f} A, mean bias {es['bias_mean']:+.3f} A, "
          f"mean max |z| {es['max_z_mean']:.2f}\n")
    print(f"{'surrogate':<16}{'mean':>8}{'median':>8}{'<2A':>7}{'<2.5A':>8}"
          f"{'FAIL18':>9}{'vs real (negative = the destroyed property was costly)':>56}")
    for a in ARMS:
        s = out["arms"][a]
        v = s["vs_real"]
        print(f"{a:<16}{s['mean']:>8.3f}{s['median']:>8.3f}{s['frac_under_2.0']:>7.2f}"
              f"{s['frac_under_2_5']:>8.2f}{s['FAIL18']:>9.3f}"
              f"{v['mean_diff']:>+31.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print("\n(`gauss_matched` is the decisive control: if it lands on `real`, the distance"
          "\n channel has NO structure cost and the phase-diagram gap was an MAE-to-sigma"
          "\n conversion artefact of my own making. If `real` is much worse, the structure is"
          "\n real and the arms above say which property carries it.)")
    return out


if __name__ == "__main__":
    run()
