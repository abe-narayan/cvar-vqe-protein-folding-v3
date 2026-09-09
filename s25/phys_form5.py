"""S25 / PHYSICS LANE / FORM 5 -- THE PHYSICS AS A PER-TARGET AUDIT, NOT AS A RANKER.

    python s25/phys_form5.py

Pre-registered in `s25/PREREG_PHYS.md` ADDENDUM A1-A7, written to disk before this ran.
Authorised by the coordinator. The architecture is FROZEN: this cannot change what ships unless
it clears its falsifier decisively, in which case the result goes back to the coordinator.

THE QUESTION
============
The functional lever has failed as a FILTER (s24 D1-C), a PARTITION (s24 D1), a FITTED SCORE
(s24 L16, oracle ceiling 0.0148 A) and an EQUAL-WEIGHT SCORE (s25 suite, no free parameter).
Form 5 is the one remaining class: the physics never chooses a candidate, never splits the pool
and never enters the ranking. It supplies a PER-TARGET SIGN AT INFERENCE, which is the only
leverage `in-band-ordering-is-per-target` leaves open.

    s_t = mean_{i in T_t} [ zrank(E_LEG)_i - zrank(E_AMB)_i ]        T_t = C3's own top-75

MECHANISM, STATED BEFORE THE RUN
================================
s25 L1 measured the distogram predicting distances systematically TOO LONG, growing to -0.589 A
at long separation -- its posterior favours over-expanded structures. The s25 landscape measured
the two physics channels as an OPPOSITE-SIGNED compactness pair on this exact pool: Legacy's
top-75 is 0.758 A MORE COMPACT than the pool, AMBER's is 1.103 A MORE EXPANDED, rho(LEG, Rg) =
+0.60 against rho(AMB, Rg) = -0.27. So the signed Legacy-minus-AMBER rank contrast OF THE
DISTOGRAM'S OWN SELECTED SET is a native-free, scale-free statistic that asks whether the
distogram's answer sits where the physics says the pool's feasible region is -- without either
energy ever ranking anything.

THE FALSIFIER, UNCHANGED
========================
The switched arm must beat C3 by MORE THAN ITS OWN MDE with the FOLD-CLUSTERED CI EXCLUDING
ZERO, under nested LOFO CV. If it does not, Form 5 is closed and the functional lever is closed
in all FIVE of its forms.

NO PER-TARGET FITTING. The threshold AND the direction are fitted on the OTHER FOUR FOLDS ONLY.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s15 import seed as SD                 # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s25 import phys_lib as P              # noqa: E402

CELLS = os.path.join(P.RESULTS, "suite_cells")
QGRID = np.arange(0.05, 0.96, 0.05)        # training-fold quantiles for the threshold


def rg_of(W):
    c = W - W.mean(1, keepdims=True)
    return np.sqrt((c ** 2).sum(2).mean(1))


def nested_lofo(sig, a, b, folds):
    """Switch between arm `a` (=C3) and arm `b` (=C5) by thresholding `sig`.

    THRESHOLD AND DIRECTION ARE FITTED ON THE OTHER FOUR FOLDS ONLY. The grid includes the two
    degenerate ends, so the fit is free to decline to switch at all -- without them a "gain"
    could be the grid being unable to express "do nothing".
    """
    out = np.empty(len(a))
    picked = {}
    for f in sorted(set(folds.tolist())):
        te, tr = folds == f, folds != f
        best, arg = None, None
        # degenerate ends first, so ties resolve toward NOT switching
        for cand in [("always_a", None, None), ("always_b", None, None)] + \
                    [("thr", d, float(np.quantile(sig[tr], q)))
                     for d in (+1, -1) for q in QGRID]:
            kind, d, tau = cand
            if kind == "always_a":
                pick = np.zeros(len(a), bool)
            elif kind == "always_b":
                pick = np.ones(len(a), bool)
            else:
                pick = (sig > tau) if d > 0 else (sig < tau)
            v = float(np.where(pick[tr], b[tr], a[tr]).mean())
            if best is None or v < best - 1e-12:
                best, arg = v, cand
        kind, d, tau = arg
        if kind == "always_a":
            pick = np.zeros(len(a), bool)
        elif kind == "always_b":
            pick = np.ones(len(a), bool)
        else:
            pick = (sig > tau) if d > 0 else (sig < tau)
        out[te] = np.where(pick[te], b[te], a[te])
        picked[int(f)] = dict(kind=kind, direction=(int(d) if d else None),
                              tau=(float(tau) if tau is not None else None),
                              train_mean=float(best),
                              n_switched_heldout=int(pick[te].sum()),
                              n_heldout=int(te.sum()))
    return out, picked


def main() -> int:
    # ------------------------------------------------------------- the arms, from the suite
    rows = []
    for f in sorted(glob.glob(os.path.join(CELLS, "*.json"))):
        rows.extend(json.load(open(f))["rows"])
    pdbs = P.targets()
    def arm(cid):
        d = {r["pdb"]: r["rmsd_vqe"] for r in rows if r["norm"] == "rank" and r["config"] == cid}
        return np.array([d[p] for p in pdbs], float)
    C3, C5 = arm(3), arm(5)
    folds = ST.pinned_folds(pdbs)

    # ------------------------------------------------------------- the signals, NATIVE-FREE
    s, r = np.empty(len(pdbs)), np.empty(len(pdbs))
    for i, p in enumerate(pdbs):
        cand, ch = P.channels(p)
        T = np.argsort(P.combine(ch, ("DIS",), "rank"), kind="stable")[:P.M_PROD]
        s[i] = float((P.zrank(ch["LEG"]) - P.zrank(ch["AMB"]))[T].mean())
        g = rg_of(cand.W)
        r[i] = float(g[T].mean() / g.mean() - 1.0)

    print("S25 FORM 5 -- THE PHYSICS AS A PER-TARGET AUDIT OF THE DISTOGRAM")
    print(f"n = {len(pdbs)}   C3 {C3.mean():.4f} A   C5 {C5.mean():.4f} A   "
          f"(rank normalisation, CVaR-VQE arm, point cloud)")
    print(f"signal s_t : mean {s.mean():+.4f}  sd {s.std(ddof=1):.4f}  "
          f"range [{s.min():+.3f}, {s.max():+.3f}]   NATIVE-FREE")
    print(f"signal r_t : mean {r.mean():+.4f}  sd {r.std(ddof=1):.4f}  "
          f"range [{r.min():+.3f}, {r.max():+.3f}]   NATIVE-FREE (declared secondary)")

    out = dict(n=len(pdbs), mean_C3=float(C3.mean()), mean_C5=float(C5.mean()))

    # ------------------------------------------------------------- M1 mechanism, ORACLE DIAG
    from scipy.stats import rankdata
    adv = C5 - C3
    def sp(x, y):
        return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])
    print(f"\n--- M1 MECHANISM CHECK (ORACLE DIAGNOSTIC -- reads the native, can never "
          f"produce a result)")
    print(f"    rho(s_t, C5-C3) = {sp(s, adv):+.4f}      rho(r_t, C5-C3) = {sp(r, adv):+.4f}")
    print(f"    targets where C5 beats C3: {int((adv<0).sum())}/{len(pdbs)}   "
          f"mean advantage when it does {adv[adv<0].mean():+.4f} A, when it does not "
          f"{adv[adv>0].mean():+.4f} A")
    out["mechanism_oracle_diagnostic"] = dict(rho_s_adv=sp(s, adv), rho_r_adv=sp(r, adv),
                                              n_C5_better=int((adv < 0).sum()))

    # ------------------------------------------------------------- N2 degenerate arms
    print(f"\n--- N2 DEGENERATE ARMS")
    d_always5 = ST.compare(C5, C3, folds, names=pdbs, label="always-C5 vs always-C3")
    print(f"    always-C5 vs always-C3   {d_always5['effect']:+.4f}  "
          f"{d_always5['effect_over_mde']:+.2f}x MDE  {d_always5['verdict']}")
    out["always_C5"] = d_always5

    # ------------------------------------------------------------- PRIMARY + N1
    res = {}
    for name, sig, is_primary in (("PRIMARY  s_t = mean[zrank(LEG) - zrank(AMB)] over C3's top-75",
                                   s, True),
                                  ("SECONDARY r_t = Rg contrast of C3's top-75 vs the pool",
                                   r, False)):
        key = "primary" if is_primary else "secondary"
        sw, pick = nested_lofo(sig, C3, C5, folds)
        cmp_ = ST.compare(sw, C3, folds, names=pdbs,
                          label=f"{'PRIMARY' if is_primary else 'SECONDARY'} switched arm "
                                f"(nested LOFO) vs C3")
        # N1: permuted signal, marginal preserved, correspondence to targets destroyed
        rng = SD.stable_rng("form5", key, 0, salt=P.SALT)
        swp, pickp = nested_lofo(sig[rng.permutation(len(sig))], C3, C5, folds)
        cmpp = ST.compare(swp, C3, folds, names=pdbs,
                          label=f"N1 permuted-signal null ({key}) vs C3")
        print("\n" + "=" * 92)
        print(f"{name}")
        print("=" * 92)
        print(f"  per-fold fit (trained on the OTHER FOUR FOLDS ONLY):")
        for f, v in sorted(pick.items()):
            print(f"    fold {f}: {v['kind']:<9} dir {str(v['direction']):>4} "
                  f"tau {('%+.4f' % v['tau']) if v['tau'] is not None else '   --  '}  "
                  f"switched {v['n_switched_heldout']}/{v['n_heldout']} held-out targets")
        print(f"\n{ST.fmt(cmp_)}")
        print(f"\n  N1 permuted-signal null: {cmpp['effect']:+.4f}  "
              f"{cmpp['effect_over_mde']:+.2f}x MDE  "
              f"fold CI [{cmpp['ci95_fold'][0]:+.4f},{cmpp['ci95_fold'][1]:+.4f}]  "
              f"{cmpp['verdict']}")
        fired = (abs(cmp_["effect"]) > cmp_["mde"] and cmp_["ci95_fold"][1] < 0)
        print(f"\n  >>> FALSIFIER {'CLEARED' if fired else 'NOT CLEARED'} "
              f"(needs |effect| > MDE AND fold CI entirely below zero)")
        res[key] = dict(stats=cmp_, per_fold_fit=pick, null_permuted=cmpp,
                        cleared=bool(fired))
    out.update(res)

    # ------------------------------------------------------------- C1 ORACLE ceiling
    orc = np.minimum(C3, C5)
    kk = float(np.corrcoef(C3, C5)[0, 1])
    k_eff = 2.0 / (1.0 + kk)
    print("\n" + "=" * 92)
    print("C1 -- THE ORACLE CEILING ON ANY C3<->C5 SWITCH (FULL LEAKAGE, NOT A RULE)")
    print("=" * 92)
    print(f"  per-target min(C3, C5) = {orc.mean():.4f} A   "
          f"{orc.mean()-C3.mean():+.4f} vs C3")
    print(f"  corr(C3, C5) = {kk:.4f}  ->  k_eff = {k_eff:.2f} independent arms, not 2")
    print(f"  A best-of-K null is INAPPLICABLE at this correlation -- k_eff near 1 is exactly")
    print(f"  how you know. The ceiling stands as ORACLE and is never promoted.")
    print(f"  Nested CV recovers {res['primary']['stats']['effect']:+.4f} of "
          f"{orc.mean()-C3.mean():+.4f}; the gap is the optimism.")
    out["oracle_ceiling"] = dict(mean=float(orc.mean()),
                                 vs_C3=float(orc.mean() - C3.mean()),
                                 corr_C3_C5=kk, k_eff=k_eff)

    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    pc, sc = res["primary"]["cleared"], res["secondary"]["cleared"]
    if pc:
        print("  PRIMARY FALSIFIER CLEARED. This goes to the coordinator; the freeze is theirs.")
    else:
        print("  PRIMARY FALSIFIER NOT CLEARED. **FORM 5 IS CLOSED.**")
        print("  The functional lever is now closed in ALL FIVE of its forms: filter, "
              "partition,\n  fitted score, equal-weight score, and per-target audit.")
    print(f"  (declared secondary, Bonferroni bar at alpha/2: "
          f"{'cleared' if sc else 'not cleared'})")
    out["verdict"] = ("PRIMARY CLEARED" if pc else "FORM 5 CLOSED")

    ST.save_atomic(os.path.join(P.RESULTS, "phys_form5.json"), out,
                   complete_keys=("pdb", "fold", "s", "r", "C3", "C5"),
                   rows=[dict(pdb=p, fold=int(folds[i]), s=float(s[i]), r=float(r[i]),
                              C3=float(C3[i]), C5=float(C5[i]),
                              switched_primary=float(nested_lofo(s, C3, C5, folds)[0][i]))
                         for i, p in enumerate(pdbs)],
                   n_expected=len(pdbs), module_file=__file__)
    print(f"\nwrote s25/results/phys_form5.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
