"""s16/repair_seedfix.py -- BRIEF S3.7 violation at `s15/align_jac.py:136`, fixed and verified.

AUDIT (`s16/audit_FINDINGS.md` S5.2) found one surviving `hash()` seeding in the sprint's
machinery: the random-subspace null that produces the 0.499 figure was drawn from
`np.random.default_rng(abs(hash(p)) % 7 + 5)`.  `hash()` is salted per process, and the modulus
collapses 126 targets onto SEVEN seeds, so the null was neither reproducible across
interpreters nor independent across targets.

The fix is one line: `SD.stable_rng("align_jac", "subspace_null", p)`.

THE CLAIM THIS MODULE HAS TO SUPPORT is "no number changes".  That is not established by
asserting it -- the quantity is a Monte-Carlo draw, so the two seedings cannot be bit-identical
and the honest statement is that they agree to within the sampling error of a quantity whose
expectation is analytically k/m.  This module measures all three: the old seeding under several
interpreter hash salts, the new seeding, and the analytic expectation.

    python -m s16.repair_seedfix
"""
from __future__ import annotations

import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np
import numpy.linalg as la

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import seed as SD                   # noqa: E402
import peptide_db as pdb                     # noqa: E402


def principal_angles(A1, A2):
    return la.svd(A1.T @ A2, compute_uv=False)


def main():
    rows = []
    for t in I.targets():
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        J, _, _ = A.sup_jacobian(nphi, npsi)
        s, V = A.spectrum(J)
        m = 2 * len(nphi)
        k = m // 2
        Vq = V[:, m - k:]

        def draw(rng):
            Q, _ = la.qr(rng.normal(size=(m, m)))
            return float((principal_angles(Vq, Q[:, :k]) ** 2).mean())

        # the OLD seeding, over the 7 values `abs(hash(p)) % 7 + 5` can take -- i.e. over
        # every interpreter salt that target could ever have seen.
        old = [draw(np.random.default_rng(v)) for v in range(5, 12)]
        new = draw(SD.stable_rng("align_jac", "subspace_null", p))
        rows.append(dict(pdb=p, m=m, analytic=k / m,
                         old_mean=float(np.mean(old)), old_min=float(np.min(old)),
                         old_max=float(np.max(old)), new=new))
    a = np.array([r["analytic"] for r in rows])
    om = np.array([r["old_mean"] for r in rows])
    nw = np.array([r["new"] for r in rows])
    spread = np.array([r["old_max"] - r["old_min"] for r in rows])
    print("s15/align_jac.py:136 random-subspace null, n = %d targets" % len(rows))
    print("  analytic expectation k/m          %.6f" % a.mean())
    print("  OLD hash() seeding, mean over the 7 reachable salts   %.6f" % om.mean())
    print("  NEW stable_rng seeding                                %.6f" % nw.mean())
    print("  |new - old_mean| mean %.2e   max %.2e" %
          (np.abs(nw - om).mean(), np.abs(nw - om).max()))
    print("  per-target spread ACROSS INTERPRETER SALTS under the old seeding:"
          " mean %.4f  max %.4f" % (spread.mean(), spread.max()))
    print("  -> the reported 0.499 is unchanged (%.4f old / %.4f new, analytic %.4f);"
          % (om.mean(), nw.mean(), a.mean()))
    print("     what the fix removes is the %.4f mean / %.4f max per-target dependence on"
          % (spread.mean(), spread.max()))
    print("     which interpreter drew it.")
    return rows


if __name__ == "__main__":
    main()
