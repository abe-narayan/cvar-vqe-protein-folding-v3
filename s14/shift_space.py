"""SHIFT agent -- STEP 4.  Is the DISTRIBUTION worth more than the point estimate?

Sprint 13 section 5.1c already established that the distributional FAMILY does not matter
(324-cell categorical, 8-component von Mises mixture and circular point regression land
within 3 degrees of one another), so that is not re-litigated here.  The untested question is
the CONSUMPTION one, and Sprint 13 section 6 showed it inverts the ranking: the arm with the
best point estimate had the WORST search space.

So this module measures the shift-derived density in the currency a torsion-constrained VQE
actually consumes -- k candidate states per residue -- using `s13/tors_support.py`'s exact
protocol so the numbers are directly comparable to that report:

    recall    fraction of residues whose native 20-degree cell is in the residue's top-k
              (DEPLOYABLE -- reads no native angle to construct, only to score)
    descent   ORACLE coordinate descent on true CA-RMSD inside the top-k space
              (ORACLE DIAGNOSTIC -- the ceiling a perfect objective could reach)

    python -m s14.shift_space
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                     # noqa: E402
from s13 import tors_common as T                    # noqa: E402
from s13 import tors_eval as EV                     # noqa: E402
from s14.shift_bmrb import RESULTS                  # noqa: E402

KS = (1, 2, 4, 8, 16)
CACHE = os.path.join(ROOT, "s14", "cache")


def space_ceiling(post, pids, nat, ks=KS, sweeps=3):
    """`s13/tors_support.support_ceiling`, restricted to a target list and reading a
    posterior dict rather than a cached arm file.  Same protocol, same 3 sweeps."""
    tgd = {t["pdb"]: t for t in I.targets()}
    rows = []
    for pid in pids:
        P = np.asarray(post[pid], np.float64)
        n = tgd[pid]["n"]
        nphi, npsi = nat[pid]
        _, _, ca = EV.pool_torsions(pid)
        nb = T.grid_bin(nphi, npsi)
        order = np.argsort(-P, 1)
        r = {"pdb": pid, "n": n, "fold": tgd[pid]["fold"]}
        for k in ks:
            top = order[:, :k]
            r["recall%d" % k] = float(np.mean([(nb[i] in top[i]) for i in range(n)]))
            cph = T.GRID_PHI[top]; cps = T.GRID_PSI[top]
            d = (np.abs(T.wrap(cph - nphi[:, None])) ** 2
                 + np.abs(T.wrap(cps - npsi[:, None])) ** 2)
            sel = np.argmin(d, 1)
            ph = cph[np.arange(n), sel].copy(); ps = cps[np.arange(n), sel].copy()
            ph[0] = nphi[0]; ps[-1] = npsi[-1]
            r["snap%d" % k] = float(I.ca_rmsd(I.build_ca(ph, ps), ca))
            cur = float(I.ca_rmsd(I.build_ca(ph, ps), ca))
            for _ in range(sweeps):
                moved = False
                for i in range(n):
                    if k == 1:
                        break
                    PH = np.tile(ph, (k, 1)); PS = np.tile(ps, (k, 1))
                    PH[:, i] = cph[i]; PS[:, i] = cps[i]
                    if i == 0:
                        PH[:, 0] = nphi[0]
                    if i == n - 1:
                        PS[:, -1] = npsi[-1]
                    v = I.kabsch_rmsd_batch(I.build_ca(PH, PS), ca)
                    j = int(np.argmin(v))
                    if v[j] < cur - 1e-9:
                        cur = float(v[j]); ph, ps = PH[j].copy(), PS[j].copy(); moved = True
                if not moved:
                    break
            r["descent%d" % k] = float(cur)
            r["qubits%d" % k] = int(n * np.log2(k)) if k > 1 else 0
        rows.append(r)
    return rows


def main():
    t0 = time.time()
    nat = EV.native_torsions()
    with open(os.path.join(RESULTS, "shift_avail_summary.json")) as fh:
        S = json.load(fh)
    pids = sorted(S["runnable_T3plus"])
    f18 = set(I.FAIL18)

    arms = ["E-SHIFT", "SHUF-SHIFT", "P-SHIFT", "REF-SHIFT", "MASK-ONLY", "SEQ-ONLY"]
    out = {"targets": pids, "n": len(pids), "ks": list(KS)}
    print("STEP 4 -- the density as a SEARCH SPACE, {} runnable targets".format(len(pids)))
    print("{:11s} {:>3s} {:>8s} {:>9s} {:>9s} {:>7s} {:>7s}".format(
        "arm", "k", "recall", "snap", "descent", "<2A", "qubits"))
    for arm in arms:
        p = os.path.join(CACHE, "shift_post_%s.npz" % arm)
        if not os.path.exists(p):
            print("  missing", p)
            continue
        z = np.load(p)
        rows = space_ceiling({k: z[k] for k in z.files}, pids, nat)
        isf = np.array([r["pdb"] in f18 for r in rows])
        a = {}
        for k in KS:
            de = np.array([r["descent%d" % k] for r in rows])
            rc = np.array([r["recall%d" % k] for r in rows])
            sn = np.array([r["snap%d" % k] for r in rows])
            a[str(k)] = {
                "recall": float(rc.mean()), "recall_FAIL18": float(rc[isf].mean()),
                "snap": float(sn.mean()), "descent": float(de.mean()),
                "descent_median": float(np.median(de)),
                "descent_frac_u2": float((de < 2.0).mean()),
                "descent_frac_u1.5": float((de < 1.5).mean()),
                "descent_FAIL18": float(de[isf].mean()),
                "mean_qubits": float(np.mean([r["qubits%d" % k] for r in rows]))}
            print("{:11s} {:3d} {:8.3f} {:9.3f} {:9.3f} {:7.2f} {:7.0f}".format(
                arm, k, a[str(k)]["recall"], a[str(k)]["snap"], a[str(k)]["descent"],
                a[str(k)]["descent_frac_u2"], a[str(k)]["mean_qubits"]))
        out[arm] = a
        out.setdefault("per_target", {})[arm] = rows

    out["what"] = ("STEP 4. recall is deployable; snap/descent are ORACLE DIAGNOSTICS "
                   "(they read native torsions to select inside the predicted space). "
                   "Protocol identical to s13/tors_support.py, 3 descent sweeps.")
    out["secs"] = round(time.time() - t0, 1)
    with open(os.path.join(RESULTS, "shift_space.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote s14/results/shift_space.json  [{:.0f}s]".format(time.time() - t0))


if __name__ == "__main__":
    main()
