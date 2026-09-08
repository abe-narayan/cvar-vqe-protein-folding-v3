"""SPRINT 15, coordinator -- THE DOSE-RESPONSE CURVE FOR DECORRELATING THE ERROR.

THE CLAIM THIS TESTS, AND WHY IT NEEDS A CURVE RATHER THAN A POINT.

K9 established that about 56% of the distogram's error is geometrically REALIZABLE -- it describes a
coherent wrong structure that the fit faithfully reproduces -- and that destroying its sign structure
while keeping every magnitude is worth about an angstrom. K12 established the requirement from the
other side: reaching 2.5 A needs an i.i.d.-class channel at 2.87 A effective RMS, against the
distogram's actual 3.70 A.

Together those say something specific and actionable: **the way to improve this channel is to
DECORRELATE its error, not to shrink it.** But "decorrelate it" is not yet a number. Nobody can act
on it without knowing how much decorrelation is needed, and whether the response is linear, threshold
-like, or saturating. A predictor-development effort that halves the coherence is a very different
project from one that must eliminate it.

THE EXPERIMENT. Interpolate continuously between the real error and a fully sign-destroyed version of
itself, and measure the emitted structure at each step:

    e(t)  =  (1 - t) * e_real  +  t * e_destroyed        t in [0, 1]

At t = 0 this is the distogram as it is (3.644 A on the full instrument). At t = 1 it is the same
error magnitudes with the coherence removed. The curve in between is the dose-response, and
inverting it gives **the fraction of coherence that would have to be removed to reach 2.5 A**.

THE CONTROL THAT MAKES IT MEAN ANYTHING. Mixing two vectors SHRINKS the result: `e(t)` has smaller
RMS than either endpoint whenever they are not parallel, because the components partially cancel. An
uncontrolled interpolation would therefore confound decorrelation with plain error reduction -- and
plain error reduction is exactly the thing K9 says is NOT the lever, so confounding them would
destroy the experiment.

Every arm is therefore run twice:

    raw        e(t) as mixed                          decorrelation AND shrinkage, confounded
    rescaled   e(t) renormalised to ||e_real||        decorrelation ALONE, at constant magnitude

**The `rescaled` ladder is the experiment; the `raw` ladder is the confound, reported beside it so
the size of the confound is visible rather than assumed away.** If `rescaled` is flat, K9 is wrong and
the gains attributed to coherence were really magnitude gains all along.

TWO DESTRUCTION MODES, because "decorrelate" is ambiguous and the two need not agree:

    signs      keep every |e|, randomise the sign          destroys direction, preserves magnitudes
    permute    permute e across pairs within the target    destroys the pair-to-pair pattern

ORACLE STATUS. Every arm needs the true error vector `e = dhat - dtrue` in order to interpolate away
from it, so **every arm here except `t = 0` is an ORACLE DIAGNOSTIC**. This measures what
decorrelation would be WORTH; it is not a method, and no arm may appear in a predictive table. The
follow-up question -- whether a predictor can be trained to have incoherent errors -- is a different
experiment that this one is meant to justify or forestall.

Run:
    python -m s15.decorr
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

TS = (0.0, 0.25, 0.5, 0.75, 1.0)
MODES = ("signs", "permute")
N_REP = 2


def destroy(e, mode, rng):
    if mode == "signs":
        return np.abs(e) * rng.choice([-1.0, 1.0], size=len(e))
    if mode == "permute":
        return e[rng.permutation(len(e))]
    raise ValueError(mode)


def mix(e_real, e_dead, t, rescale):
    """Interpolate, optionally renormalising to the real error's RMS."""
    m = (1.0 - t) * e_real + t * e_dead
    if rescale:
        r0 = float(np.sqrt((e_real * e_real).mean()))
        r1 = float(np.sqrt((m * m).mean()))
        if r1 > 1e-9:
            m = m * (r0 / r1)
    return m


def run(targets=None, n_start=4):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    data = C.gather(tg)

    arms = {}
    for mode in MODES:
        for sc in (True, False):
            for t in TS:
                arms[f"{'rescaled' if sc else 'raw'}_{mode}_t{t:g}"] = (mode, sc, t)
    res = {a: [] for a in arms}
    coh = []
    path = os.path.join(RESULTS, "decorr.json")

    for c, tt in enumerate(tg):
        p = tt["pdb"]; d = data[p]
        i, j, sd = d["i"], d["j"], d["sd"]
        dtrue, dhat = d["dtrue"], d["dhat"]
        e = dhat - dtrue
        w = 1.0 / sd ** 2
        starts = D.starts(p, d["seq"], d["n"], d["fold"], n_start, SD.stable_rng(p, "decorr"))

        for a, (mode, sc, t) in arms.items():
            reps = 1 if t == 0.0 else N_REP
            vals = []
            for r in range(reps):
                rng = SD.stable_rng(p, mode, t, r)
                em = mix(e, destroy(e, mode, rng), t, sc)
                tgt = np.maximum(dtrue + em, 2.0)
                best = None
                for phi0, psi0, _tag in starts:
                    phi, psi, f, _ = D.fit_distances(tgt, w, i, j, phi0, psi0)
                    if best is None or f < best[0]:
                        best = (f, float(I.ca_rmsd(I.build_ca(phi, psi), d["nat"])))
                vals.append(best[1])
            res[a].append(float(np.mean(vals)))
        coh.append({"pdb": p, "rms": float(np.sqrt((e * e).mean()))})
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": res, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    base = np.asarray(res[f"rescaled_{MODES[0]}_t0"], float)
    out = {"n": len(pdbs), "incumbent": float(ref.mean()), "ts": list(TS),
           "arms": {}, "ladders": {}}
    for a in arms:
        v = np.asarray(res[a], float)
        out["arms"][a] = {**I.summary(v), "median": float(np.median(v)),
                          "frac_under_2_5": float((v < 2.5).mean()),
                          "vs_t0": I.paired(v, base, folds=folds, names=pdbs)}
    for sc in ("rescaled", "raw"):
        for mode in MODES:
            out["ladders"][f"{sc}_{mode}"] = [out["arms"][f"{sc}_{mode}_t{t:g}"]["mean"]
                                              for t in TS]
    #: invert the rescaled ladder: what fraction of the coherence must go to reach each target?
    out["required_t"] = {}
    for mode in MODES:
        ys = out["ladders"][f"rescaled_{mode}"]
        for target in (3.204, 2.5, 2.0):
            if target >= ys[0]:
                v = "already"
            elif target <= ys[-1]:
                v = f"<={TS[-1]:g}"
            else:
                v = f"{float(np.interp(-target, [-y for y in ys], list(TS))):.3f}"
            out["required_t"][f"{mode}_{target}"] = v
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s15_decorr", out, n_expected=len(tg))

    print(f"\nn = {len(pdbs)}   incumbent {ref.mean():.3f}")
    print("EVERY ARM EXCEPT t=0 IS AN ORACLE DIAGNOSTIC -- it measures what decorrelation would "
          "be WORTH.\n")
    print(f"{'ladder':<22}" + "".join(f"{'t=' + format(t, 'g'):>10}" for t in TS))
    for sc in ("rescaled", "raw"):
        for mode in MODES:
            y = out["ladders"][f"{sc}_{mode}"]
            tag = f"{sc}/{mode}"
            print(f"{tag:<22}" + "".join(f"{v:>10.3f}" for v in y))
    print("\n('rescaled' holds the error RMS fixed and is THE experiment; 'raw' lets the mix "
          "shrink\n the error and is the confound, shown so its size is visible.)\n")
    print("fraction of the coherence that must be destroyed to reach each target "
          "(rescaled ladder):")
    for mode in MODES:
        row = "  ".join(f"{t} A -> t = {out['required_t'][f'{mode}_{t}']}"
                        for t in (3.204, 2.5, 2.0))
        print(f"  {mode:<10} {row}")
    return out


if __name__ == "__main__":
    run()
