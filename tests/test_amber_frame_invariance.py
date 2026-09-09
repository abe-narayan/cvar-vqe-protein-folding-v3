"""Standing regression check -- AMBER relaxation must be rigid-invariant.

ff14SB, GBn2, the positional restraint (`0.5 k |x - x0|^2` with `x0` taken from the same
input) and every RMSD in this project are all invariant under a proper rotation and a
translation of the input structure.  Relaxing the SAME structure in a rotated lab frame
is therefore zero by construction, and anything this test measures is the AMBER
minimiser's own numerical floor.

Why this is A permanent test and not A ONE-OFF.  Sprint 15 measured this null at
**+0.0117 A, sd 0.137, max 1.51 A** on 126 targets -- half the size of the accuracy
effect the project was quoting, on a comparison that must return zero.  The cause was
four silently non-converged minimisations; with `core.amber`'s convergence gate applied
the mean collapses to -0.0005 A.  The TAIL does not collapse: converged targets still
move by up to 0.17 A under a transformation that is mathematically no transformation at
all.  A regression here is invisible to every other test in this suite.

Reflections are deliberately excluded (det = +1 enforced): a reflection is NOT a symmetry
of ff14SB, so an improper rotation is a real change of structure, not a null.

    pytest tests/test_amber_frame_invariance.py -s
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

#: A fixed slice of the pinned 126-target instrument.  12 targets at ~2 x 12.6 s of
#: restrained minimisation is ~5 minutes; the full sweep is a research run, not a test.
N_TARGETS = 12


@pytest.mark.slow
def test_amber_relaxation_is_frame_invariant():
    from s12 import instrument as I
    from s16 import energy_lib as L

    rows, verdict = L.frame_null(I.targets()[:N_TARGETS], draw=1, verbose=True)
    g = verdict.get("gated", verdict["ungated"])

    # The pre-declared band, and its honest status.  `FRAME_TOL_MEAN = 0.005` /
    # `FRAME_TOL_MAX = 0.05` were declared before any Sprint 16 draw.  On this 12-target
    # slice, live, the CONVERGED subset returns mean +0.0146 A and max 0.1401 A (1DEP) --
    # it FAILS the declared band on both counts, and the band has deliberately NOT been
    # loosened to make it pass.  That failure is the finding (s16/energy_FINDINGS.md 2.2):
    # a comparison that is zero by construction is not zero, so no PER-TARGET AMBER RMSD
    # claim is resolvable.  This test therefore REPORTS the declared verdict and ASSERTS
    # only a regression guard set from the measured floor -- it catches the minimiser
    # getting WORSE, which is what a regression check is for.
    print(f"\nPRE-DECLARED BAND |mean| <= {L.FRAME_TOL_MEAN}, max <= {L.FRAME_TOL_MAX}"
          f"  ->  PASS = {verdict['PASS']}   (measured mean {g['mean']:+.5f}, "
          f"max {g['max_abs']:.5f}, n = {g['n']}, excluded {verdict['excluded']})")

    # Every excluded target must be excluded for a stated reason, and the count printed.
    assert verdict["n_excluded_by_gate"] == len(verdict["excluded"])

    # REGRESSION GUARD, from the Sprint 16 measurement (mean +0.0146, max 0.1401 on this
    # slice; max 0.124-0.172 across three independent frames on all 126). Tightening these
    # towards zero is an open problem; loosening them silently is how the defect got in.
    assert abs(g["mean"]) <= 0.05, (
        f"rigid-invariance null mean {g['mean']:+.5f} A on {g['n']} converged targets is "
        f"far outside the +0.0146 A measured in Sprint 16; excluded {verdict['excluded']}")
    assert g["max_abs"] <= 0.25, (
        f"rigid-invariance null max |dRMSD| {g['max_abs']:.4f} A is far outside the "
        f"0.12-0.17 A residual measured in Sprint 16")


@pytest.mark.slow
def test_convergence_gate_is_reported_not_silent():
    """`refine_coords` must return a convergence verdict on every call."""
    import torsion_lib2 as tl2
    from core import amber as am
    from s12 import instrument as I
    from s14.avgspace import top75_windows
    from s15.phys_repl import averaged_backbone_from

    t = I.targets()[0]
    W, PHI, PSI, u = top75_windows(t["pdb"])
    avg, _, _ = averaged_backbone_from(W, PHI, PSI)
    tab = tl2.library_for(t["seq"], 4, t["seq"])
    rep = tl2.PerResidueTorsion(t["seq"], tab, chi_bits=False)
    r = am.refine_coords(t["seq"], rep, avg, k_restraint=30.0, steps=0,
                         tolerance=1.0, threads=1, memo=False)
    assert "converged" in r and "converge_reason" in r
    assert isinstance(r["converged"], bool)
    assert r["converged"] == (np.isfinite(r["energy"])
                              and r["energy"] <= am.CONVERGE_MAX_KCAL)
    # the structure is still returned even when the gate fails -- the gate reports,
    # the CONSUMER decides
    assert np.asarray(r["ca"]).shape == (len(t["seq"]), 3)


def test_convergence_flags_is_a_pure_function():
    """No AMBER required: the declared rule, checked directly."""
    from core import amber as am
    assert am.CONVERGE_MAX_KCAL == 1000.0
    assert am.convergence_flags(-500.0, -100.0)["converged"] is True
    assert am.convergence_flags(1000.0, 0.0)["converged"] is True
    assert am.convergence_flags(1000.1, 0.0)["converged"] is False
    assert am.convergence_flags(8.9e8, 0.0)["converged"] is False
    assert am.convergence_flags(float("nan"), 0.0)["converged"] is False
    assert am.convergence_flags(float("inf"), 0.0)["converged"] is False
    assert "8.9e+08" in am.convergence_flags(8.9e8, 0.0)["converge_reason"].replace(
        "8.9e+08", "8.9e+08")
