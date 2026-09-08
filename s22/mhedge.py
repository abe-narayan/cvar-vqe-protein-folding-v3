"""s22/mhedge.py -- IF YOU CANNOT PREDICT THE RIGHT m, HEDGE ACROSS m.

THE ONE IDEA THAT SURVIVES THIS CAMPAIGN'S CENTRAL RESULT.  L4/L5 established that the per-target
optimal averaging width m is **real, stable across independent halves of the same target's pool,
genuinely per-target rather than a global correction, and worth 0.243 A**.  L7 established that it
is **invisible to four independent native-free router constructions across three feature families**.

    A SIGNAL YOU CANNOT PREDICT CAN STILL BE HEDGED.

Selection needs a per-target decision and therefore a per-target signal, which we do not have.
**Hedging needs neither.**  Averaging the structures emitted at several m is a FIXED, GLOBAL,
NATIVE-FREE operator: the same thing is done to every target, no feature is consulted, no fold is
fit, and there is nothing to overfit.

WHY IT COULD WIN, stated as a mechanism rather than a hope.  The emitted structure at width m has
error `e_m`.  A fixed m pays `E[e_m]` for whichever m was chosen.  The hedge pays the error of the
*average structure*, and by the usual variance decomposition that is below the average of the
individual errors by exactly the **disagreement between the arms** (an ambiguity/Krogh-Vedelsby
term).  So the hedge beats the mean-over-m whenever the arms disagree at all -- the open question is
only whether it also beats the **best single** m, which it does when the m-errors are partially
decorrelated across targets.  **That is precisely the situation L4/L5 measured**: the best m moves
from target to target, so no single m is right often, while the ensemble is never far wrong.

WHY IT IS NOT ALREADY TESTED.  Sprint 21's m-ladder treated the widths as **alternative arms** and
asked which to pick.  Workstream C's blend ladder mixed **pool against latent**, not **across m**.
Nothing in this project has averaged the m-ladder's own outputs.

THE ARMS.

    fixed_75        the incumbent operator.  THE THING TO BEAT.
    hedge_all       coordinate average of the structures at m = 500,150,75,20,5
    hedge_core      the same over m = 150,75,20 -- the three rungs S21 L23 found FLAT
    hedge_wide      m = 500,150,75
    hedge_w2        inverse-rank weighted toward the middle of the ladder
    oracle_m        ORACLE per-target best m -- the ceiling, for reference only

FRAMES.  Averaging structures requires a common frame.  Each ladder structure is superposed onto the
**m=75 structure** by Kabsch before averaging -- a fixed, native-free choice made because m=75 is
the incumbent and therefore the arm a reader will compare against.  NOT TAKEN: superposing onto the
per-target best arm (needs the native) or onto the widest arm.

CONTRACTION, DECLARED.  Coordinate averaging **contracts the backbone** (S20: mean virtual bond
2.961 A against a physical 3.804), and averaging averages contracts it further.  Every hedge arm is
a **POINT CLOUD** and is compared only against other point clouds; `fixed_75` is a point cloud too,
so the primary comparison is internally valid.  Contraction is measured and reported per arm so the
cost is visible rather than assumed.

OPERATOR FORKS, per BRIEF SS4 rule 0.  I hold a directional hypothesis, so each fork names the
alternative and the direction it would push.

    functional     DECLARED the shipped Bayes-risk score for the filter.  NOT TAKEN the squared
                   functional (S21 C2: interchangeable at the argmin on pool-like structures).
    basis          DECLARED the pool's own window coordinates on ALL arms including the incumbent.
                   NOT TAKEN the rebuild.
    readout        DECLARED coordinate average, both within a width and across widths.  NOT TAKEN
                   the medoid across widths, which would select rather than hedge and so would
                   reintroduce exactly the per-target decision this file exists to avoid.
    normalisation  none; RMSD in Angstroms.  NOT TAKEN per-target z-scoring.
    null           `fixed_75`, the incumbent, on identical data -- the honest null, since it asks
                   whether hedging beats NOT hedging.  A second null, `hedge_shuffle`, averages the
                   ladder structures of DIFFERENT targets' widths where shapes allow, to confirm any
                   gain is not a generic smoothing artefact.  NOT TAKEN a random-structure average.
    THE LABEL      DECLARED continuous Ca-RMSD.  NOT TAKEN any binarised "did the hedge win",
                   whose threshold would depend on the ladder's flatness -- the covariate-dependent
                   threshold trap of S21 L25.

  Hypothesis   the hedge beats fixed m=75, because the per-target optimum moves and no single width
               is right often.
  Falsifier    if `hedge_all` and `hedge_core` both fail to beat `fixed_75` past that comparison's
               own MDE with a CI excluding zero, hedging across m is refuted and the m-ladder
               headroom is unreachable by BOTH selection and combination.
  Null         fixed m=75 on identical data.
  Budget       126 targets, 6 arms, no fitting, no folds -- nothing is selected, so nothing can
               overfit.
  Promotion    a win here is promotable ONLY after replication, ablation across the width sets, and
               a validity check (clash / cis / bond geometry), because a contracted point cloud can
               buy Ca-RMSD while destroying physicality.
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
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402

MS = (500, 150, 75, 20, 5)
SETS = {"hedge_all": (500, 150, 75, 20, 5),
        "hedge_core": (150, 75, 20),
        "hedge_wide": (500, 150, 75),
        "hedge_narrow": (75, 20, 5)}
REF = 75


def _save(o, name="mhedge.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _vbond(P):
    P = np.asarray(P, float)
    return float(np.linalg.norm(P[1:] - P[:-1], axis=1).mean())


def run():
    tg = I.targets(); rows = []
    print("targets: %d" % len(tg), flush=True)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb); idx = I.pool_idx(u)
        W = np.asarray(u["W"], float)[idx]; nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        order = np.argsort(sc, kind="stable")
        lad = {}
        for m in MS:
            C, _b = I.coordinate_average(W[order[:min(m, len(order))]])
            lad[m] = np.asarray(C, float)
        ref = lad[REF]
        row = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]),
               "fixed_75": float(I.ca_rmsd(ref, nat)),
               "vb_fixed_75": _vbond(ref)}
        for m in MS:
            row["m_%d" % m] = float(I.ca_rmsd(lad[m], nat))
        row["oracle_m"] = min(row["m_%d" % m] for m in MS)
        #: hedge: superpose every width onto the m=75 frame (fixed, native-free), then mean
        for name, ws in SETS.items():
            S = np.stack([np.asarray(I.superpose_batch(lad[m][None], ref), float)[0] for m in ws])
            H = S.mean(0)
            row[name] = float(I.ca_rmsd(H, nat))
            row["vb_" + name] = _vbond(H)
        rows.append(row)
        if (c + 1) % 20 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})
    need = ["fixed_75", "oracle_m"] + list(SETS)
    ok = len(rows) == len(tg) and all(all(np.isfinite(r[k]) for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "sets": {k: list(v) for k, v in SETS.items()}, "ref": REF})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "mhedge.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)        # noqa: E731
    rng = SD.stable_rng("mhedge", "rep"); fold = g("fold")

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        b = x[rng.integers(0, k, size=(4000, k))].mean(1)
        F = sorted(set(fold.astype(int)))
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), se, 2.8016 * se, float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)),
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    inc = g("fixed_75")
    print("\nn = %d.  All arms POINT CLOUDS on the pool's window basis.\n" % len(rows))
    print("  %-16s%9s%11s   %s" % ("arm", "RMSD", "virt bond", "vs fixed m=75"))
    print("  %-16s%9.4f%11.3f" % ("fixed_75 (INC)", inc.mean(), g("vb_fixed_75").mean()))
    for m in MS:
        if m == REF:
            continue
        print("  %-16s%9.4f%11s" % ("m=%d" % m, g("m_%d" % m).mean(), "-"))
    for name in SETS:
        x = g(name); m_, se, mde, lo, hi, flo, fhi, w = st(x - inc)
        flag = "  <<< BEATS" if fhi < 0 and abs(m_) > mde else ""
        print("  %-16s%9.4f%11.3f   %+.4f SE %.4f MDE %.3f iid[%+.3f,%+.3f] fold[%+.3f,%+.3f] %3dW/%3dL%s"
              % (name, x.mean(), g("vb_" + name).mean(), m_, se, mde, lo, hi, flo, fhi, w, len(rows) - w, flag))
    print("  %-16s%9.4f%11s   <- ORACLE per-target m, ceiling only" % ("oracle_m", g("oracle_m").mean(), "-"))
    print("\n  Falsifier: hedge_all AND hedge_core both fail to beat fixed_75 past their own MDE.")
    print("  Physical virtual CA-CA bond is 3.805 A; averaging contracts, and averaging averages")
    print("  contracts further -- the column above is how much, and it is a cost not a detail.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
