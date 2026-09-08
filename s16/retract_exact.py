"""SPRINT 16 / RETRACT -- TASK A part 2: is the ambiguity identity EXACT through Kabsch?

THE QUESTION.  LIT (`s16/lit_FINDINGS.md` B.5) establishes that

    d_avg^2 = (r1^2 + r2^2)/2 - (s/2)^2                                   (*)

is the two-member Krogh-Vedelsby (1995) ambiguity decomposition, EXACT in any inner-product space,
with no coplanarity and no equidistance assumption.  That is a statement about a FIXED FRAME.

Sprint 15 does not evaluate it in a fixed frame.  Every one of the four quantities is measured after
its OWN optimal superposition:
  r1 = ca_rmsd(Xd, nat)          Kabsch rotation R1
  r2 = ca_rmsd(Xp, nat)          Kabsch rotation R2 != R1
  s  = ca_rmsd(Xd, Xp)           Kabsch rotation R3, in general != R2 R1^-1
  d  = ca_rmsd(coordinate_average(Xd,Xp), nat), and `coordinate_average` itself superposes onto
       the medoid before averaging -- a FIFTH rotation.
Kabsch superposition is a minimisation, hence nonlinear and non-commuting with averaging.  So (*)
has no right to hold, and whether it does is an empirical question about real structures.

WHAT THIS MODULE DOES.  It reproduces the two fitted structures of `s15/coherence.py` exactly (same
starts, same seeds, same optimiser) and evaluates the identity under four frame conventions:

  A1  S15_pipeline    exactly Sprint 15's convention (five independent superpositions)
  A2  common_frame    superpose Xd and Xp onto the native ONCE; then s and d are PLAIN RMSDs in
                      that single frame.  (*) must be exact here if LIT is right.
  A3  s_common_only   pipeline d, common-frame s          } single-substitution arms that attribute
  A4  d_common_only   pipeline s, common-frame d          } the pipeline residual to one operation

Both means are reported for every arm: quadratic sqrt((r1^2+r2^2)/2) -- what the identity requires --
and arithmetic (r1+r2)/2 -- what Sprint 15 used.

NOT AN ORACLE ARM IN THE PREDICTIVE SENSE, but it does read natives: this is an AUDIT of a published
number, permitted by BRIEF 3.1 as auditing / post-hoc interpretation.  Nothing here enters a
predictor.

Run:  python -m s16.retract_exact
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s15 import distcal as C               # noqa: E402
from s15 import distgeo as D               # noqa: E402
from s15 import pooldist as P              # noqa: E402
from s15 import seed as SD                 # noqa: E402

RES = os.path.join(ROOT, "s16", "results")
os.makedirs(RES, exist_ok=True)
PATH = os.path.join(RES, "retract_exact.json")


def _fit(d_target, w, d, starts):
    """VERBATIM from s15/coherence.py._fit, minus the residual return."""
    i, j = d["i"], d["j"]
    best = None
    for phi0, psi0, _tag in starts:
        phi, psi, f, _ = D.fit_distances(np.maximum(d_target, 2.0), w, i, j, phi0, psi0)
        if best is None or f < best[0]:
            best = (f, phi, psi)
    return np.asarray(I.build_ca(best[1], best[2]), float)


def _plain_rmsd(a, b):
    """RMSD with NO superposition -- the inner-product-space quantity the identity is about."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.sqrt(((a - b) ** 2).sum() / len(a)))


def run(targets=None, n_start=4):
    tg = targets if targets is not None else I.targets()
    data = C.gather(tg)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        n = d["n"]; nat = np.asarray(d["nat"], float)
        starts = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p, "coherence"))
        Dm, _s, _i, _j = P.pool_distances(p, n)
        d_pool = np.median(Dm, axis=0)
        wp = 1.0 / np.maximum(Dm.std(axis=0), 0.25) ** 2
        Xd = _fit(d["dhat"], 1.0 / d["sd"] ** 2, d, starts)
        Xp = _fit(d_pool, wp, d, starts)

        # --- the four measured quantities, S15's convention -------------------------------
        r1 = I.ca_rmsd(Xd, nat)
        r2 = I.ca_rmsd(Xp, nat)
        s_k = I.ca_rmsd(Xd, Xp)
        Cm, _b = I.coordinate_average(np.stack([Xd, Xp]))
        d_k = I.ca_rmsd(Cm, nat)

        # --- the same quantities in ONE common frame (both superposed onto the native) -----
        Xd2, Xp2 = I.superpose_batch(np.stack([Xd, Xp]), nat)
        natc = nat - nat.mean(0)
        Xd2 = Xd2 - Xd2.mean(0); Xp2 = Xp2 - Xp2.mean(0)
        r1c = _plain_rmsd(Xd2, natc)            # must equal r1
        r2c = _plain_rmsd(Xp2, natc)            # must equal r2
        s_c = _plain_rmsd(Xd2, Xp2)             # >= s_k
        Cc = 0.5 * (Xd2 + Xp2)
        d_c = _plain_rmsd(Cc, natc)             # >= ca_rmsd(Cc, nat)
        d_ck = I.ca_rmsd(Cc, nat)               # common-frame average, then Kabsch to native

        rows.append({"pdb": p, "n": int(n), "fold": int(d["fold"]),
                     "r1": float(r1), "r2": float(r2), "s_kabsch": float(s_k),
                     "d_kabsch": float(d_k), "r1_common": float(r1c), "r2_common": float(r2c),
                     "s_common": float(s_c), "d_common": float(d_c),
                     "d_common_then_kabsch": float(d_ck)})
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            json.dump({"partial": rows, "n_done": c + 1}, open(PATH, "w"))
            print(f"  {c+1}/{len(tg)}  {time.time()-t0:.0f}s", flush=True)

    g = lambda k: np.asarray([r[k] for r in rows], float)      # noqa: E731
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    r1, r2 = g("r1"), g("r2")
    qm2 = 0.5 * (r1 ** 2 + r2 ** 2)
    am2 = (0.5 * (r1 + r2)) ** 2

    def arm(s, dd, m2):
        pred = np.sqrt(np.maximum(m2 - (s / 2.0) ** 2, 0.0))
        e = dd - pred
        return {"pred_mean": float(pred.mean()),
                "measured_mean": float(dd.mean()), "measured_median": float(np.median(dd)),
                "signed_err_mean": float(e.mean()), "signed_err_median": float(np.median(e)),
                "abs_err_mean": float(np.abs(e).mean()), "abs_err_median": float(np.median(np.abs(e))),
                "abs_err_max": float(np.abs(e).max()), "abs_err_sd": float(np.abs(e).std(ddof=1)),
                "n_under": int((e > 0).sum()), "n_over": int((e < 0).sum()),
                "max_abs_err_target": str(pdbs[int(np.argmax(np.abs(e)))])}

    arms = {}
    for mn, m2 in (("quad", qm2), ("arith", am2)):
        arms[f"A1_S15_pipeline_{mn}"] = arm(g("s_kabsch"), g("d_kabsch"), m2)
        arms[f"A2_common_frame_{mn}"] = arm(g("s_common"), g("d_common"), m2)
        arms[f"A3_s_common_only_{mn}"] = arm(g("s_common"), g("d_kabsch"), m2)
        arms[f"A4_d_common_only_{mn}"] = arm(g("s_kabsch"), g("d_common"), m2)
        arms[f"A5_commonavg_then_kabsch_{mn}"] = arm(g("s_common"), g("d_common_then_kabsch"), m2)

    out = {
        "n": len(rows), "arms": arms,
        "reproduction_check": None,
        "kabsch_buys_on_s": I.summary(g("s_common") - g("s_kabsch")),
        "kabsch_buys_on_d": I.summary(g("d_common") - g("d_common_then_kabsch")),
        "medoid_avg_vs_native_frame_avg": I.summary(g("d_kabsch") - g("d_common_then_kabsch")),
        "r1_common_minus_r1": I.summary(g("r1_common") - r1),
        "r2_common_minus_r2": I.summary(g("r2_common") - r2),
        "A1_vs_A2_abs_err_quad": I.paired(
            np.abs(g("d_kabsch") - np.sqrt(np.maximum(qm2 - (g("s_kabsch") / 2) ** 2, 0.0))),
            np.abs(g("d_common") - np.sqrt(np.maximum(qm2 - (g("s_common") / 2) ** 2, 0.0))),
            folds=folds, names=pdbs),
        "rows": rows,
    }

    # determinism: does the reproduction match the persisted S15 artefact?
    try:
        old = {r["pdb"]: r for r in json.load(
            open(os.path.join(ROOT, "s15", "results", "coherence.json")))["rows"]}
        dl = {k: float(np.abs(np.asarray([old[r["pdb"]][ko] for r in rows], float)
                              - g(kn)).max())
              for k, ko, kn in (("rmsd_disto", "rmsd_disto", "r1"),
                                ("rmsd_pool", "rmsd_pool", "r2"),
                                ("disagreement", "disagreement", "s_kabsch"),
                                ("rmsd_coordavg", "rmsd_coordavg", "d_kabsch"))}
        out["reproduction_check"] = {"max_abs_deviation_vs_s15_coherence": dl}
    except Exception as ex:                                     # pragma: no cover
        out["reproduction_check"] = {"error": repr(ex)}

    json.dump(out, open(PATH, "w"), indent=1)
    I.write("s16_retract_exact", out, n_expected=len(tg))

    print(f"\nn = {len(rows)}   reproduction vs s15/results/coherence.json:")
    for k, v in (out["reproduction_check"].get("max_abs_deviation_vs_s15_coherence") or {}).items():
        print(f"    max |delta| {k:16s} {v:.3e}")
    print("\nIS THE IDENTITY EXACT?  (signed err = measured - predicted; quadratic mean)")
    for k in ("A1_S15_pipeline", "A2_common_frame", "A3_s_common_only",
              "A4_d_common_only", "A5_commonavg_then_kabsch"):
        a = arms[f"{k}_quad"]
        print(f"  {k:26s} |err| mean {a['abs_err_mean']:.3e} median {a['abs_err_median']:.3e} "
              f"max {a['abs_err_max']:.3e}  under/over {a['n_under']}/{a['n_over']}")
    print("\n  same arms with S15's ARITHMETIC mean:")
    for k in ("A1_S15_pipeline", "A2_common_frame"):
        a = arms[f"{k}_arith"]
        print(f"  {k:26s} |err| mean {a['abs_err_mean']:.3e} median {a['abs_err_median']:.3e} "
              f"max {a['abs_err_max']:.3e}  under/over {a['n_under']}/{a['n_over']}")
    for k in ("kabsch_buys_on_s", "kabsch_buys_on_d", "medoid_avg_vs_native_frame_avg"):
        v = out[k]
        print(f"  {k:32s} mean {v['mean']:+.4f} median {v['median']:+.4f} max {v['max']:+.4f}")
    return out


if __name__ == "__main__":
    run()
