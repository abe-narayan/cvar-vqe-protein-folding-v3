"""s24/c_bias.py -- SIGNED DISTANCE BIAS AT n=126, and the tension the coordinator flagged.

THE TENSION.  The coordinator's `referent.py` measures beta = 0.520: about half the distogram's
error reappears in the emitted output, against a mismatched-distogram placebo floor of 0.352.  My
exploratory numbers have the distogram predicting distances too LONG (signed mean +0.367 A) while
the emitted cloud is too SHORT (-0.783 A).  Opposite signs, and yet half the error transfers.

THE RESOLUTION I AM TESTING.  These are not the same statistic.  `beta` and the cosine are
properties of the full PAIR-ERROR VECTOR; the signed mean is one scalar component of it -- the
scale/compactness mode.  A source can inherit the prior's error DIRECTION while a SEPARATE,
operator-induced bias moves the scale the other way.  The candidate operator is already named in
this project: coordinate averaging CONTRACTS the backbone (25.8% measured), and contraction
shortens every pair distance at once.

THE DECISIVE DECOMPOSITION.  Measure the signed pair-distance bias at three stages of the SAME
pipeline, on the SAME pairs:

    members   the 75 SELECTED windows, individually, before any averaging
    average   the emitted uniform coordinate average of exactly those 75
    delta     average - members  =  what the AVERAGING OPERATOR alone contributes

If `members` is near zero or positive and `average` is strongly negative, the contraction is the
operator, not the prior, and there is no contradiction with beta.  Reported alongside the
projection of each error vector onto the distogram's error direction, which is the component beta
actually prices, and onto the scale mode, which is the component the signed mean prices.

n=126, complete.  ORACLE labelling: every quantity here reads the native and is a DIAGNOSTIC,
never an input to a decision.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from s24 import c_probe as CP            # noqa: E402

TOPM = 75
OUT = os.path.join(RES, "c_bias.json")


def run():
    tg = I.targets()
    rows = []
    for t in tg:
        pdb = t["pdb"]; n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n)
        dnat = I.pair_dists(nat[None], i, j)[0]
        g = np.asarray(dg["expected"], float) - dnat            # distogram's signed error

        idx = I.pool_idx(u); Wp = u["W"][idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, i, j)), float)
        top = np.argsort(sc, kind="stable")[:TOPM]
        A = Wp[top]
        Dm = I.pair_dists(A, i, j)                              # (75, npairs)
        em = Dm.mean(0) - dnat                                  # members' mean signed error
        cA = CP._avg(A)
        ea = I.pair_dists(cA[None], i, j)[0] - dnat             # the emitted average's error

        # the whole universe, as the "what real geometry looks like" reference
        eu = I.pair_dists(u["W"], i, j).mean(0) - dnat

        def proj(e):
            """component of e along the distogram's error direction (what beta prices)"""
            ng = float(g @ g)
            return float((e @ g) / ng) if ng > 0 else float("nan")

        rows.append({
            "pdb": pdb, "n": n, "fold": int(u["fold"]),
            "rmsd": float(I.ca_rmsd(cA, nat)),
            "dgram_signed": float(g.mean()), "dgram_mae": float(np.abs(g).mean()),
            "members_signed": float(em.mean()),
            "average_signed": float(ea.mean()),
            "delta_avg_minus_members": float(ea.mean() - em.mean()),
            "universe_signed": float(eu.mean()),
            "proj_members": proj(em), "proj_average": proj(ea),
            "cos_members_dgram": CP._cos(em, g), "cos_average_dgram": CP._cos(ea, g),
            "rms_members": float(np.sqrt((em ** 2).mean())),
            "rms_average": float(np.sqrt((ea ** 2).mean())),
        })
        if len(rows) % 25 == 0:
            print("  %d/%d" % (len(rows), len(tg)), flush=True)
    tmp = OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "complete": len(rows) == len(tg), "n_expected": len(tg)}, fh)
    os.replace(tmp, OUT)
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    F = sorted(set(fold.tolist()))
    rng = SD.stable_rng("c_bias", "rep")
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731

    def st(x):
        x = np.asarray(x, float); se = x.std(ddof=1) / np.sqrt(len(x))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return x.mean(), se, 2.8016 * se, np.percentile(fs, 2.5), np.percentile(fs, 97.5)

    print("\nn=%d.  SIGNED mean pair-distance error, in Angstroms.  Negative = TOO SHORT." % len(rows))
    for k, lab in (("dgram_signed", "the DISTOGRAM's own prediction"),
                   ("universe_signed", "the whole legal universe (real geometry)"),
                   ("members_signed", "the 75 SELECTED members, before averaging"),
                   ("average_signed", "the EMITTED coordinate average"),
                   ("delta_avg_minus_members", "  -> what AVERAGING alone contributes")):
        m, se, mde, lo, hi = st(g(k))
        print("  %-42s %+8.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]" % (lab, m, se, mde, lo, hi))
    print("  distogram MAE %.4f" % g("dgram_mae").mean())

    print("\n  PROJECTION onto the distogram's error direction (the component beta prices)")
    for k, lab in (("proj_members", "members"), ("proj_average", "emitted average")):
        m, se, mde, lo, hi = st(g(k))
        print("    %-16s %+8.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]" % (lab, m, se, mde, lo, hi))
    print("  COSINE with the distogram's error direction")
    for k, lab in (("cos_members_dgram", "members"), ("cos_average_dgram", "emitted average")):
        m, se, mde, lo, hi = st(g(k))
        print("    %-16s %+8.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]" % (lab, m, se, mde, lo, hi))
    print("\n  RMS pair-distance error   members %.4f   emitted average %.4f"
          % (g("rms_members").mean(), g("rms_average").mean()))
    print("\n  READING: if 'members' is near zero and 'average' is strongly negative, the")
    print("  contraction is the AVERAGING OPERATOR, not the prior, and the opposite signs are")
    print("  not in tension with beta -- beta prices the error DIRECTION, the signed mean prices")
    print("  the SCALE mode, and they are different components of the same vector.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
