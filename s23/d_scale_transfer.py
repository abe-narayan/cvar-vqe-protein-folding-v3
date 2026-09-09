"""s23/d_scale_transfer.py -- WORKSTREAM D: DOES THE PER-TARGET ORACLE SCALE TRANSFER?

BRIEF duty (1): "check whether the effect survives if the scale is fitted on one half of a
target's pool and applied to the other, which is the transfer test that saved a similar claim
last sprint" (s22/mreal.py, m-ladder headroom, 65% transferred).

WHY THIS TEST IS NECESSARY AND DIFFERENT FROM gscale.py's "oracle".  d_scale_closedform.py shows
the per-target oracle s* is a ONE-PARAMETER LEAST-SQUARES fit of the SAME coordinate average that
is then scored -- SSE(1)-SSE(s*) = (Cc-Ct)^2/Cc >= 0 ALWAYS (equality iff s*=1 exactly), so a
strict, near-universal "win" is mechanically guaranteed regardless of whether s* reflects a real,
reusable property of the target or pure in-sample fitting.  The only way to tell the two apart is
to fit s on data DISJOINT from the data being scored -- exactly mreal.py's construction, applied to
scale instead of m.

CONSTRUCTION.  Each target's K=500 pool splits into two disjoint random halves of 250.  Within each
half, apply the SAME score-filter + top-75 selection the incumbent uses, build the coordinate
average.  On half A: find the analytic oracle scale s*_A (closed form, fit using half A's own
average against native).  Apply s*_A (UNCHANGED, not refit) to half B's average and score against
native.  Compare to the FIXED s=1.0 on half B (matched data, no scaling).  Average over 8 random
half-splits and both directions (A->B, B->A), as mreal.py does.  Also report the GLOBAL scale
(mean s*_A over all targets in the training direction) applied out-of-sample, as a second point of
comparison against gscale.py's nested-CV global arm.

OPERATOR FORKS (BRIEF Rule 0).  Directional hypothesis (partial, not full, transfer expected, by
analogy with mreal.py's 65%).
    functional     DECLARED shipped Bayes-risk score, matching gscale.py.        NOT TAKEN squared.
    basis          DECLARED pool window coordinates (point cloud), both sides.   NOT TAKEN rebuild.
    readout        DECLARED top-75 coordinate average within each half, then a scalar about the
                   CENTROID -- matches gscale.py exactly.                        NOT TAKEN medoid frame.
    normalisation  DECLARED closed-form analytic scale (no grid truncation --    NOT TAKEN the
                   d_scale_closedform.py audit found gscale.py's [0.90,1.10]     81-pt grid, which
                   grid truncates 50% of targets).                               truncates.
    null           DECLARED fixed s=1.0 on the SAME half-B data.                 NOT TAKEN half-A's
                                                                                   own in-sample oracle
                                                                                   (begs the question).
    THE LABEL      DECLARED continuous held-out Ca-RMSD.                         NOT TAKEN win/loss
                                                                                   binarisation.

  Hypothesis   the analytic oracle scale s* is PARTLY a real per-target property (some transfer,
               like mreal.py's 65%) but the transferred magnitude is much smaller than the
               in-sample oracle -0.34 A, and is concentrated in the FAIL18 subset (pre-registered
               subgroup, since d_scale_closedform.py already found 57% of the in-sample ceiling
               sits in these 18/126 targets).
  Falsifier    if held-out transfer beats fixed s=1.0 past its own MDE with a fold-clustered CI
               excluding zero, on the FULL panel AND on non-FAIL18 alone, the scale correction is
               real and general, not a FAIL18 artefact.
  Null         fixed s=1.0, same half.
  Budget       126 targets x 8 repeats x 2 halves.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402

TOPM = 75
REPS = 8
FAIL18 = {"1ID6", "1JBF", "1LB7", "2BFI", "2BP4", "2JN5", "2MQ2", "2N5C", "2NB7", "2NDM",
          "3BTB", "3SGO", "5W52", "7JS6", "7LCW", "8T63", "9KAR", "9L1M"}


def _save(o, name="d_scale_transfer.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def kabsch_CcCtTt(C, nat):
    """Fix Kabsch rotation at s=1 (scale-invariant), return Cc, Ct, Tt, n for the quadratic."""
    Cc0 = C - C.mean(0); Tc0 = nat - nat.mean(0)
    H = Cc0.T @ Tc0
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.eye(3); D[2, 2] = d
    R = Vt.T @ D @ U.T
    Crot = Cc0 @ R.T
    n = Crot.shape[0]
    Cc = float((Crot ** 2).sum()); Ct = float((Crot * Tc0).sum()); Tt = float((Tc0 ** 2).sum())
    return Cc, Ct, Tt, n


def rmsd_at(Cc, Ct, Tt, n, s):
    sse = max(Cc * s * s - 2 * Ct * s + Tt, 0.0)
    return float(np.sqrt(sse / n))


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, reps: %d" % (len(tg), REPS), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        rng = SD.stable_rng(pdb, "s23dscaletransfer")
        rep = []
        for _r in range(REPS):
            perm = rng.permutation(len(W))
            halves = (perm[:len(W) // 2], perm[len(W) // 2:])
            quad = []
            for h in halves:
                o = h[np.argsort(sc[h], kind="stable")]
                sel = o[:min(TOPM, len(o))]
                C, _b = I.coordinate_average(W[sel])
                quad.append(kabsch_CcCtTt(np.asarray(C, float), nat))
            (CcA, CtA, TtA, nA), (CcB, CtB, TtB, nB) = quad
            sA = CtA / CcA; sB = CtB / CcB
            rep.append({
                "rmsd1_A": rmsd_at(CcA, CtA, TtA, nA, 1.0),
                "rmsd1_B": rmsd_at(CcB, CtB, TtB, nB, 1.0),
                "oracleA": rmsd_at(CcA, CtA, TtA, nA, sA),
                "oracleB": rmsd_at(CcB, CtB, TtB, nB, sB),
                "transfer_AtoB": rmsd_at(CcB, CtB, TtB, nB, sA),   # s from A, scored on B
                "transfer_BtoA": rmsd_at(CcA, CtA, TtA, nA, sB),   # s from B, scored on A
                "s_star_A": sA, "s_star_B": sB,
            })
        agg = lambda k: float(np.mean([x[k] for x in rep]))  # noqa: E731
        rows.append({
            "pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]), "fail18": pdb in FAIL18,
            "fixed_heldout": 0.5 * (agg("rmsd1_B") + agg("rmsd1_A")),
            "sel_heldout": 0.5 * (agg("transfer_AtoB") + agg("transfer_BtoA")),
            "oracle_insample": 0.5 * (agg("oracleA") + agg("oracleB")),
            "s_star_mean": 0.5 * (agg("s_star_A") + agg("s_star_B")),
        })
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("fixed_heldout", "sel_heldout", "oracle_insample", "s_star_mean")
    ok = len(rows) == len(tg) and all(all(np.isfinite(r[k]) for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "reps": REPS, "topm": TOPM})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "d_scale_transfer.json")))["rows"]
    rng = SD.stable_rng("d_scale_transfer", "rep")

    def st(x, fold):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        b = x[rng.integers(0, k, size=(4000, k))].mean(1)
        F = sorted(set(fold.astype(int)))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se, float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)),
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    def block(rs, label):
        g = lambda k: np.array([r[k] for r in rs], float)  # noqa: E731
        fold = g("fold")
        print("\n=== %s (n=%d) ===" % (label, len(rs)))
        print("  %-42s%9s" % ("arm", "RMSD"))
        for k, lab in (("fixed_heldout", "FIXED s=1.0, held out"),
                       ("sel_heldout", "TRANSFER s* (half A) -> half B"),
                       ("oracle_insample", "IN-SAMPLE oracle s* (the suspect)")):
            print("  %-42s%9.4f" % (lab, g(k).mean()))
        m, se, mde, lo, hi, flo, fhi, w = st(g("sel_heldout") - g("fixed_heldout"), fold)
        verdict = "TRANSFERS" if (fhi < 0 and abs(m) > mde) else "DOES NOT TRANSFER"
        print("  PRIMARY  transfer - fixed, HELD OUT : %+.4f SE %.4f MDE %.3f iid[%+.4f,%+.4f] fold[%+.4f,%+.4f] %dW/%dL"
              % (m, se, mde, lo, hi, flo, fhi, w, len(rs) - w))
        print("           -> %s" % verdict)
        m2, se2, mde2, lo2, hi2, _flo2, _fhi2, _w2 = st(g("oracle_insample") - g("fixed_heldout"), fold)
        print("  in-sample oracle - fixed            : %+.4f [%+.4f,%+.4f]  <- the APPARENT ceiling (this half-pool construction)"
              % (m2, lo2, hi2))
        if m2 < 0:
            frac = 100.0 * min(max(m / m2, 0.0), 1.0)
            print("  FRACTION OF THE APPARENT CEILING THAT TRANSFERS: %.0f%%" % frac)
        return m, mde, fhi

    print("\nn = %d targets, %d split-half repeats, half-pool of 250, top-75 within each half." % (len(rows), REPS))
    block(rows, "FULL PANEL")
    fail = [r for r in rows if r["fail18"]]
    ok = [r for r in rows if not r["fail18"]]
    block(fail, "FAIL18 subset (18 known zero-recall targets)")
    block(ok, "non-FAIL18 (108 targets)")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
