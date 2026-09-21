#!/usr/bin/env python
"""s32/s32_Q_verify.py -- lane Q's numbers, asserted against the artefacts that produced them.

Contract: *numbers in the ledger are asserted by a verifier against the artefact that produced
them*.  This is lane Q's half, written so lane V can call `check()` from `s32/s32_verify.py`
without importing anything heavy.

It also re-derives, from the RAW rows rather than from the summary, the three quantities a
reader is most likely to quote:

  * the Q0 block model reproduces `E` on every target (the claim S31 could not test)
  * Q1-T2's in-hull gain is 1 and the orthogonal gain is 0, with the two DISTINGUISHABLE
  * Q1-T3's cardinality statistic -- whether the s-of-K constraint binds

The self-test that can fail (contract rule 5): `check()` is run twice, once on the real rows and
once on rows with one field CORRUPTED, and the corrupted pass MUST fail.  A verifier that cannot
fail is decoration.

    python s32/s32_Q_verify.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def _rows(name):
    with open(os.path.join(RES, name)) as fh:
        return [json.loads(ln) for ln in fh if ln.strip()]


def check(q0_rows=None, q1_rows=None):
    """Return (ok, findings).  Each finding is (name, ok, detail)."""
    out = []
    q0 = q0_rows if q0_rows is not None else _rows("s32_Q0_rows.jsonl")
    q1 = q1_rows if q1_rows is not None else _rows("s32_Q1_rows.jsonl")

    # ---- Q0
    out.append(("Q0 n == 126", len(q0) == 126, len(q0)))
    out.append(("Q0 sc[o] non-decreasing on every target",
                all(r["monotone"] for r in q0),
                sum(1 for r in q0 if not r["monotone"])))
    bm = max(r["blockmodel_err"] for r in q0)
    out.append(("Q0 block model reproduces E exactly (<1e-12)", bm < 1e-12, bm))
    md = max(r["maxdev_ramp"] for r in q0)
    out.append(("Q0 max|E-ramp| == 0.0407 to 4dp", abs(md - 0.0407) < 5e-5, md))
    iff = all((r["maxdev_ramp"] > 1e-12) == (r["n_tied_excess"] > 0) for r in q0)
    out.append(("Q0 E deviates IFF the target has ties", iff, None))
    out.append(("Q0 selection readout identical tie-vs-ramp on 126/126",
                all(r["selidx_tie"] == r["selidx_ramp"] for r in q0),
                sum(1 for r in q0 if r["selidx_tie"] != r["selidx_ramp"])))
    out.append(("Q0 dim == 128 reachable on every target",
                all(r["len_order"] >= 128 for r in q0),
                min(r["len_order"] for r in q0)))
    out.append(("Q0 (alpha,T) takes exactly two values",
                len({(r["alpha"], r["T"]) for r in q0}) == 2,
                sorted({(r["alpha"], r["T"]) for r in q0})))

    # ---- Q1-T1
    out.append(("Q1 n == 126", len(q1) == 126, len(q1)))
    for key, bar in (("err_identity_rel", 1e-10), ("err_a_affine_in_t_rel", 1e-10),
                     ("err_sufficient_w", 1e-10), ("err_sufficient_x_rmsd", 1e-10)):
        v = max(r[key] for r in q1)
        out.append(("Q1-T1 %s < %g" % (key, bar), v < bar, v))
    out.append(("Q1 KKT certificate < 1e-8", max(r["kkt"] for r in q1) < 1e-8,
                max(r["kkt"] for r in q1)))

    # ---- Q1-T2, and the two gains must be DISTINGUISHABLE (else the check is decoration)
    gi = np.array([r["gain_in_mean"] for r in q1 if np.isfinite(r["gain_in_mean"])])
    go = np.array([r["gain_out_mean"] for r in q1])
    out.append(("Q1-T2 in-hull gain == 1 (|g-1| < 1e-5)",
                bool(np.abs(gi - 1.0).max() < 1e-5), float(np.abs(gi - 1.0).max())))
    out.append(("Q1-T2 orthogonal gain == 0 (< 1e-5)",
                bool(go.max() < 1e-5), float(go.max())))
    out.append(("Q1-T2 the two gains are DISTINGUISHABLE (>0.5 apart)",
                bool(gi.mean() - go.mean() > 0.5), float(gi.mean() - go.mean())))
    sr = max(r["sens_rel_max"] for r in q1)
    out.append(("Q1-T2 dx == P_aff(S) dt (rel < 1e-4)", sr < 1e-4, sr))

    # ---- Q1-T3
    if "k500_support" in q1[0]:
        sup = np.array([r["k500_support"] for r in q1], float)
        out.append(("Q1-T3 K=500 convex support recorded", True,
                    {"mean": float(sup.mean()), "median": float(np.median(sup)),
                     "max": int(sup.max()), "p90": float(np.percentile(sup, 90)),
                     "frac_le_10": float((sup <= 10).mean())}))
    return all(o[1] for o in out), out


def main():
    ok, out = check()
    for name, good, detail in out:
        print("%-58s %s   %s" % (name, "OK " if good else "FAIL", detail))
    # the self-test: corrupt one field and demand the verifier NOTICES
    q0 = _rows("s32_Q0_rows.jsonl")
    q1 = _rows("s32_Q1_rows.jsonl")
    q1[0] = dict(q1[0])
    q1[0]["gain_out_mean"] = 0.9                       # an orthogonal gain of 0.9 is impossible
    ok_bad, _ = check(q0_rows=q0, q1_rows=q1)
    print("\nself-test: corrupted rows must FAIL ->", "PASS" if not ok_bad else "BROKEN")
    print("VERDICT:", "ALL CHECKS PASS" if ok else "FAILURES ABOVE")
    sys.exit(0 if (ok and not ok_bad) else 1)


if __name__ == "__main__":
    main()
