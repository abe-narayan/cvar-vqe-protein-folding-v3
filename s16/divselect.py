"""s16/divselect.py -- the ONE constructive experiment the sprint's mechanism suggests.

WHERE THIS COMES FROM.  `s16/diversity.py` established, on the enumerated instrument, that the
terminal coordinate-average operator obeys the Krogh-Vedelsby decomposition exactly,

    ||avg - native||^2  =  mean_i ||x_i - native||^2  -  mean_i ||x_i - avg||^2
                           <-- member error -->          <-- DIVERSITY -->

and that six generators with mean member error spanning 0.056 A produced readouts spanning
0.301 A, the whole difference being the diversity term.  The incumbent pipeline selects its 75
windows by distogram score alone -- that is, it optimises the FIRST term and ignores the
second, which the identity says contributes with equal weight and opposite sign.

THE HYPOTHESIS.  Selecting for score AND diversity should beat selecting for score alone.

THE OBVIOUS OBJECTION, AND WHY IT IS THE POINT.  Admitting more diverse windows admits worse
ones, so the first term rises as the second does.  The identity does not say diversity is free;
it says the two trade at equal weight.  The experiment measures WHERE the trade turns, and
because both terms are computable per arm, a null here is still informative: we will know
whether it failed because diversity did not rise or because member error rose faster.

NATIVE-FREE, AND EXACTLY SO.  Diversity is a pairwise RMSD among candidate windows -- no native.
The score is the shipped distogram score -- no native.  The mixing weight `lam` and the count
`m` are chosen LEAVE-FOLD-OUT on labels, which is a trained hyperparameter and is declared as
such (the same convention as `s15/distcal.fit_correction`), not as a native-free quantity.

CONTROLS, both of the two AUDIT made mandatory after showing this programme's evidence had
read ~5x stronger than it was:
    rand{m}     matched-count random selection from the same candidate pool (2 draws)
    maxdiv{m}   pure diversity, score ignored entirely -- the ZERO-INFORMATION extreme of the
                mixing dial, which must be bad if the score means anything
and the incumbent's own top-75 is the reference every arm is priced against.

FALSIFIER, WRITTEN BEFORE THE RUN.  If no (m, lam) beats the shipped top-75 selection out of
fold with an interval excluding zero, the diversity mechanism -- though exact as arithmetic --
does not yield a usable selection rule, and the sprint reports it as a mechanism without a
lever.  A leave-fold-out selector that lands on lam = 0 is that null in its clearest form.

ONE TRAP THIS MODULE MUST AVOID.  `s16/csteer.py` produced two arms that looked like large
improvements and were endpoint substitutions, detected because their leave-fold-out step sat on
the boundary of the ladder.  Here the analogue is `lam` or `m` pinned at a ladder end, and the
report prints the chosen values per fold for exactly that reason.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

CAND = 200          # native-free pre-filter: the top-CAND of the K=500 pool by score
MS = (20, 50, 75, 120)
LAMS = (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.8, 1.5)


def _zs(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / (s if s > 1e-12 else 1.0)


def _greedy(score_z, P, m, lam):
    """maximise  -score_z(i) + lam * (mean pairwise RMSD to the chosen set).

    Greedy from the best-scoring window.  `P` is the candidate pairwise-RMSD matrix, so the
    diversity term is entirely native-free.
    """
    n = len(score_z)
    sel = [int(np.argmin(score_z))]
    if lam <= 0.0:
        return np.argsort(score_z, kind="stable")[:m]
    chosen = np.zeros(n, bool); chosen[sel[0]] = True
    dsum = P[sel[0]].copy()
    for _ in range(m - 1):
        gain = -score_z + lam * (dsum / len(sel))
        gain[chosen] = -np.inf
        j = int(np.argmax(gain))
        sel.append(j); chosen[j] = True; dsum += P[j]
    return np.asarray(sel, int)


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    rows = []
    for c, t in enumerate(tg):
        p = t["pdb"]; seq = t["seq"]; fold = int(t["fold"])
        u = I.load_univ(p)
        pool = I.pool_idx(u)
        W = np.asarray(u["W"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        n = W.shape[1]
        i, j = I.pair_index(n)
        dg = I.distogram(p, seq, fold)
        score = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)

        keep = np.argsort(score, kind="stable")[:CAND]
        Wc = W[keep]; sc = _zs(score[keep])
        P = I.pairwise_rmsd(Wc)
        rng = SD.stable_rng(p, "s16divsel")

        def readout(sel):
            if len(sel) < 2:
                return float(I.ca_rmsd(Wc[sel][0], nat)), float("nan"), float("nan")
            Ps = P[np.ix_(sel, sel)]
            a, b = I.coordinate_average(Wc[sel], Ps)
            Wm = I.superpose_batch(Wc[sel], Wc[sel][b])
            div = float(np.sqrt(np.mean(((Wm - a) ** 2).sum(axis=(1, 2)) / n)))
            per = np.array([I.ca_rmsd(w, nat) for w in Wc[sel]])   # ORACLE, diagnostic
            return float(I.ca_rmsd(a, nat)), float(np.sqrt((per ** 2).mean())), div

        arms = {}
        for m in MS:
            for lam in LAMS:
                r, e, d = readout(_greedy(sc, P, m, lam))
                arms[f"g{m}_{lam}"] = {"rmsd": r, "err": e, "div": d}
            rr = [readout(rng.permutation(CAND)[:m]) for _ in range(2)]
            arms[f"rand{m}"] = {"rmsd": float(np.mean([x[0] for x in rr])),
                                "err": float(np.mean([x[1] for x in rr])),
                                "div": float(np.mean([x[2] for x in rr]))}
            #: pure diversity, score ignored -- the zero-information extreme.
            r, e, d = readout(_greedy(np.zeros_like(sc), P, m, 1.0))
            arms[f"maxdiv{m}"] = {"rmsd": r, "err": e, "div": d}

        #: the incumbent's own selection, read through the identical operator so the
        #: comparison is of the SELECTION and not of the readout.
        sub = np.asarray(I.shipped_record(p)["sub"], int)
        pos = {int(k): q for q, k in enumerate(keep)}
        ship = np.array([pos[s] for s in sub if int(s) in pos], int)
        r, e, d = readout(ship)
        arms["shipped75"] = {"rmsd": r, "err": e, "div": d, "n": int(len(ship))}

        rows.append({"pdb": p, "fold": fold, "n": int(n), "arms": arms})
        if (c + 1) % 10 == 0:
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)
            json.dump({"rows": rows}, open(os.path.join(RESULTS, "divselect.json"), "w"))

    json.dump({"rows": rows}, open(os.path.join(RESULTS, "divselect.json"), "w"))
    report(rows)
    return rows


def _boot(d, rng, B=4000):
    d = np.asarray(d, float); k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows):
    rng = SD.stable_rng("divselect", "report")
    folds = np.asarray([r["fold"] for r in rows], int)
    ship = np.array([r["arms"]["shipped75"]["rmsd"] for r in rows])
    print(f"\nn = {len(rows)}   shipped top-75 through this operator: {ship.mean():.3f} A")
    print("  (the incumbent's 3.204 A adds the projection and AMBER stages this module omits,")
    print("   so the reference here is the shipped SELECTION read through the SAME operator.)\n")

    print("A. THE TRADE.  Both Krogh-Vedelsby terms per arm; `err` is ORACLE, `div` is native-free.")
    print(f"  {'arm':<14}{'RMSD':>9}{'member err':>12}{'diversity':>11}{'vs shipped':>12}")
    for m in MS:
        for lam in LAMS:
            k = f"g{m}_{lam}"
            v = np.array([r["arms"][k]["rmsd"] for r in rows])
            print(f"  {k:<14}{v.mean():>9.3f}"
                  f"{np.mean([r['arms'][k]['err'] for r in rows]):>12.3f}"
                  f"{np.mean([r['arms'][k]['div'] for r in rows]):>11.3f}"
                  f"{v.mean() - ship.mean():>+12.3f}")
        for k in (f"rand{m}", f"maxdiv{m}"):
            v = np.array([r["arms"][k]["rmsd"] for r in rows])
            print(f"  {k:<14}{v.mean():>9.3f}"
                  f"{np.mean([r['arms'][k]['err'] for r in rows]):>12.3f}"
                  f"{np.mean([r['arms'][k]['div'] for r in rows]):>11.3f}"
                  f"{v.mean() - ship.mean():>+12.3f}   CONTROL")
        print()

    print("B. THE DEPLOYABLE ARM: (m, lam) chosen LEAVE-FOLD-OUT on labels.")
    keys = [f"g{m}_{lam}" for m in MS for lam in LAMS]
    picked, vals = {}, np.zeros(len(rows))
    for f in sorted(set(folds)):
        tr = folds != f
        picked[int(f)] = min(keys, key=lambda k: np.mean(
            [rows[q]["arms"][k]["rmsd"] for q in range(len(rows)) if tr[q]]))
    for q, r in enumerate(rows):
        vals[q] = r["arms"][picked[int(folds[q])]]["rmsd"]
    a, lo, hi = _boot(vals - ship, rng)
    w = int((vals < ship).sum()); l = int((vals > ship).sum())
    print(f"  chosen per fold: {[picked[f] for f in sorted(picked)]}")
    print(f"  vs shipped top-75:  {a:+.3f} [{lo:+.3f},{hi:+.3f}]   median "
          f"{np.median(vals - ship):+.3f}   {w}W/{l}L")
    print("\n  A fold landing on lam = 0.0 is the null in its clearest form: the selector is")
    print("  declining to use diversity at all.  A fold pinned at the LADDER END (lam = 1.5 or")
    print("  m = 120) is the csteer trap -- the arm is choosing an endpoint, not a trade.")


if __name__ == "__main__":
    run()
