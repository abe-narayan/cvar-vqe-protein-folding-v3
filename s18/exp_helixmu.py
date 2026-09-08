"""s18/exp_helixmu.py -- the ZERO-INFORMATION REFERENCE-MEASURE control for the degree-1 object.

WHY THIS ARM EXISTS.  The coordinator's `priorfit` (n = 126, complete) found that the retrieval
pool's per-residue torsion prior -- **the best torsion channel measured anywhere in this project**
-- is statistically indistinguishable from a CONSTANT IDEAL ALPHA-HELIX at every rung of its
ladder.  ADVERSARIAL separately found a constant-helix START is not worse than the coordinate
average at alpha = 1.  Three zero-information controls in this sprint have now matched their
informative arms, and the programme's own memory already records that on this instrument *a
zero-information constant alpha-helix beats the random control, so "beats random" proves nothing
here*.

That lands squarely on my lane, because **MATH's `E_le1` is DEFINED against a reference measure
`mu`, and the shipped `mu` is the target's own retrieval-pool per-residue torsion marginal** --
a conditioned torsion channel, sitting inside the objective's definition rather than beside it.
So the statement "degree-1 beats its zero-information controls, therefore it carries real
information" is NOT interpretable until the same object is built against a reference measure with
no positional information at all.

THE CONTROL.  `mu_helix` puts every residue's torsions at the ideal alpha-helix
`(-57 deg, -47 deg)` (`core.project.STARTS[1]`), with a small isotropic jitter so the conditional
expectations are not evaluated at a single point mass.  It is native-free, sequence-blind,
position-blind and target-blind: the SAME measure for every residue of every target.

HOW IT IS BUILT WITHOUT TOUCHING MATH'S ESTIMATOR.  The reference measure is the one modelling
knob the BRIEF explicitly says is a choice ("the choice of mu is a modelling decision and must be
native-free.  Report sensitivity to it").  So this module supplies a new NAME to
`math_anova.mu_samples` by wrapping it -- every other name delegates to the original function,
and MATH's `fit`, tabulation, interpolation and argmin are used exactly as shipped.  Nothing in
the estimator is re-implemented.

READING, pre-declared:
  * helix-mu degree-1 lands NEAR pool-mu degree-1  -> the pool's target conditioning contributes
    nothing to this object, and "degree-1 carries information" collapses to "generic backbone
    plausibility carries information".  A fourth zero-information control matching its arm.
  * helix-mu degree-1 lands MUCH WORSE            -> the conditioning is real, and the arm's
    failure is about the truncation rather than about the measure.

Either way the sprint's verdict is unchanged -- degree-1 fails against Control A by +1.287 A
[+1.124, +1.473] under the shipped mu -- so this control settles an INTERPRETATION, not the
falsification.

    python -m s18.exp_helixmu [n_subset]
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

from core import project as PJ                # noqa: E402
from s12 import instrument as I               # noqa: E402
from s18 import exp_obj as XO                 # noqa: E402
from s18 import exp_run as XR                 # noqa: E402
from s18 import math_anova as MA              # noqa: E402
from s18 import math_lib as ML                # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "exp_helixmu.json")
HELIX = np.deg2rad(np.asarray(PJ.STARTS[1], float))     # (-57, -47) degrees
JITTER = np.deg2rad(10.0)

_ORIG = MA.mu_samples


def _mu_samples(pdb, seq, fold, n, name, S, rng):
    """MATH's `mu_samples` plus one extra name.  Every other name delegates unchanged."""
    if name != "helix":
        return _ORIG(pdb, seq, fold, n, name, S, rng)
    phi = HELIX[0] + JITTER * rng.standard_normal((S, n))
    psi = HELIX[1] + JITTER * rng.standard_normal((S, n))
    return phi, psi


MA.mu_samples = _mu_samples


def run(subset=None, out=OUT, S=None):
    S = S or MA.NSAMP
    full = I.targets()
    deb = MA.debias_map(full)                 # ALWAYS the full 126 -- never a subset
    tg = full
    if subset:
        ns = np.array([t["n"] for t in tg])
        order = np.argsort(ns, kind="stable")
        tg = [tg[q] for q in np.sort(order[np.linspace(
            0, len(tg) - 1, int(subset)).astype(int)])]
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
        pdb = t["pdb"]
        if pdb in done:
            continue
        u = I.load_univ(pdb)
        nat = np.asarray(u["nat_ca"], float)
        _W, _PH, _PS, avg, _ca, phi0, psi0 = XR.start_structure(t, nat)
        ob, _d = XR._target_obj(t, deb, "helix", S)
        e = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
             "avg": float(I.ca_rmsd(avg, nat)), "E0": float(ob.t.E0)}
        p, q, f, _nf, _ni = XO.lbfgs(ob.E_le1, phi0, psi0)
        e["helix_le1"] = float(I.ca_rmsd(I.build_ca(p, q), nat))
        ap, aq = ob.argmin_le1(hold=(phi0, psi0))
        e["helix_argmin"] = float(I.ca_rmsd(I.build_ca(ap, aq), nat))
        e["helix_argmin_objfull"] = ob.E_full(ap, aq)[0]
        rows.append(e)
        json.dump({"rows": rows, "complete": False, "S": S, "mu": "helix",
                   "subset": (len(tg) if subset else None)}, open(out, "w"))
        if len(rows) % 5 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s) {pdb}", flush=True)
    json.dump({"rows": rows, "complete": len(rows) == len(tg), "S": S, "mu": "helix",
               "subset": (len(tg) if subset else None), "n_expected": len(tg),
               "n_rows": len(rows)}, open(out, "w"))
    print(f"DONE helix-mu {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


if __name__ == "__main__":
    run(subset=int(sys.argv[1]) if len(sys.argv) > 1 else 30)
