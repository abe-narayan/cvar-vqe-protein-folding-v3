"""SPRINT 14, ENER-11 -- TAIL-RESTRICTED discrimination, and why bulk accuracy misleads.

The coordinator observed a paradox in the discrimination table: the pure distogram has the
best bulk pairwise accuracy of any objective measured (0.654 at a sub-0.25 A quality gap,
where every physical energy is at 0.51) and the WORST argmin of its own family (3.237
against 2.878 for the 1-local prior, whose bulk accuracy is chance-level at 0.505).  So
neither `pair_accuracy` over uniform pairs nor `decile_rho` at a fixed 10% predicts what an
argmin will return, and the sprint leaned on both.

This measures the regime an argmin actually samples -- pairs drawn from the objective's own
lowest few per cent -- and then decomposes the answer.

## THE PREDICTION, STATED BEFORE THE MEASUREMENT

Finding E2a established that accuracy is a function of the QUALITY GAP, not of anything
else.  Inside a tail the available gaps are small by construction.  Therefore **every
objective's tail accuracy must fall toward 0.5 as q shrinks, whether or not it has lost
skill**, and a raw tail-accuracy table would rank objectives by how heterogeneous their
tails happen to be.  So the informative statistic is the GAP-MATCHED tail accuracy: pairs
restricted to |dRMSD| in a fixed band, which is comparable across objectives and across q.
If the prediction holds, raw tail accuracy collapses and gap-matched accuracy does not.

## THE DECOMPOSITION THAT RESOLVES THE PARADOX

    sel_rmsd = pool_mean + (tail_mean - pool_mean) + (sel_rmsd - tail_mean)
                          \\____ FILTERING ____/   \\____ ORDERING ____/

Exact and additive.  An objective with strong bulk discrimination and a bad argmin has all
its value in FILTERING and a zero-or-positive ORDERING term.  That is a testable statement
about the mechanism rather than a restatement of the paradox.

Legacy and every `leg_*` term are measured on the FULL 262,144-configuration enumeration
(no sampling error at all); AMBER on the uniform stratum only (finding E0).

    python -m s14.ener_tail
"""
from __future__ import annotations

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E

QS = [0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 1.00]
GAP_BAND = (0.5, 1.5)          # a fixed window, so every cell is like-for-like
FULL_OBJ = ["legacy", "prior", "leg_steric", "leg_torsion", "legacy_nosteric",
            "leg_contact", "leg_solvation", "leg_compactness"]


def one_target(pdb):
    z = E.enum(pdb)
    R = z.rmsd
    rng = np.random.default_rng(E.SEED)
    out = {}
    for nm in FULL_OBJ:
        v = z.obj(nm)
        out[nm] = dict(
            tail=[E.tail_accuracy(v, R, q, 200000, rng, GAP_BAND) for q in QS],
            decomp=[E.selection_decomposition(v, R, q) for q in (0.01, 0.10)])
    # AMBER: uniform stratum only
    m = z.uniform_mask
    idx = z.amber_idx[m]
    Ru = R[idx]
    for nm, v in (("amber", z.amber_total[m]), ("amb_nonbonded", z.amb["nonbonded"][m])):
        out[nm] = dict(
            tail=[E.tail_accuracy(v, Ru, q, 200000, rng, GAP_BAND) for q in QS],
            decomp=[E.selection_decomposition(v, Ru, q) for q in (0.01, 0.10)],
            population="uniform stratum, n=%d" % len(idx))
    # a RANDOM-TAIL null: a tail of the same size chosen at random, same second statistic
    nulls = []
    for q in QS:
        mm = max(int(round(q * len(R))), 2)
        acc = []
        for s in range(8):
            sel = np.random.default_rng(E.SEED + s).choice(len(R), mm, replace=False)
            acc.append(E.tail_accuracy(z.legacy[sel], R[sel], 1.0, 50000,
                                       np.random.default_rng(s), GAP_BAND))
        nulls.append(dict(q=q, acc=float(np.nanmean([a["acc"] for a in acc])),
                          acc_gap_matched=float(np.nanmean(
                              [a["acc_gap_matched"] for a in acc])),
                          tail_mean_gap=float(np.mean([a["tail_mean_gap"] for a in acc]))))
    rec = dict(pdb=pdb, fold=int(z.fold), pool_mean=float(R.mean()),
               per_objective=out, legacy_random_tail_null=nulls)
    E._ENUM.pop(pdb, None)          # keep resident memory to one target at a time
    return rec


def main():
    E.wait_for_memory(1.0, "ener_tail")
    rows = [one_target(p) for p in E.ENUM_TARGETS]
    names = list(rows[0]["per_objective"])

    def agg(nm, field, qi):
        v = [r["per_objective"][nm]["tail"][qi][field] for r in rows]
        v = [x for x in v if np.isfinite(x)]
        return float(np.mean(v)) if v else float("nan")

    print("=== RAW tail accuracy: pairs inside the objective's own lowest q ===")
    hdr = "  ".join(f"{q*100:>6.1f}%" for q in QS)
    print(f"{'objective':18s} {hdr}")
    for nm in names:
        print(f"{nm:18s} " + "  ".join(
            f"{agg(nm,'acc',i):7.3f}" if np.isfinite(agg(nm, 'acc', i)) else "    nan"
            for i in range(len(QS))))
    print(f"{'RANDOM-tail null':18s} " + "  ".join(
        f"{np.mean([r['legacy_random_tail_null'][i]['acc'] for r in rows]):7.3f}"
        for i in range(len(QS))))

    print(f"\nmean |dRMSD| available inside each tail (the CONFOUND):")
    print(f"{'objective':18s} {hdr}")
    for nm in names:
        print(f"{nm:18s} " + "  ".join(f"{agg(nm,'tail_mean_gap',i):7.3f}"
                                       for i in range(len(QS))))

    print(f"\n=== GAP-MATCHED tail accuracy (pairs with |dRMSD| in "
          f"[{GAP_BAND[0]}, {GAP_BAND[1]}) A) ===")
    print(f"{'objective':18s} {hdr}")
    for nm in names:
        print(f"{nm:18s} " + "  ".join(
            f"{agg(nm,'acc_gap_matched',i):7.3f}"
            if np.isfinite(agg(nm, 'acc_gap_matched', i)) else "    nan"
            for i in range(len(QS))))
    print(f"{'RANDOM-tail null':18s} " + "  ".join(
        f"{np.mean([r['legacy_random_tail_null'][i]['acc_gap_matched'] for r in rows]):7.3f}"
        for i in range(len(QS))))

    print("\n=== SELECTION DECOMPOSITION at q = 1% and 10% ===")
    print(f"{'objective':18s} {'q':>5s} {'pool':>7s} {'tail_mean':>10s} {'sel':>7s} "
          f"{'FILTERING':>10s} {'ORDERING':>9s} {'CI(ordering)':>20s} {'W/L':>6s}")
    folds = [r["fold"] for r in rows]
    for nm in names:
        for di, q in enumerate((0.01, 0.10)):
            d = [r["per_objective"][nm]["decomp"][di] for r in rows]
            filt = np.array([x["filtering"] for x in d])
            orde = np.array([x["ordering"] for x in d])
            p = I.paired(np.array([x["sel_rmsd"] for x in d]),
                         np.array([x["tail_mean"] for x in d]), folds=folds)
            print(f"{nm:18s} {q*100:4.0f}% {np.mean([x['pool_mean'] for x in d]):7.3f} "
                  f"{np.mean([x['tail_mean'] for x in d]):10.3f} "
                  f"{np.mean([x['sel_rmsd'] for x in d]):7.3f} "
                  f"{filt.mean():+10.3f} {orde.mean():+9.3f} "
                  f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
                  f"{p['n_better']}/{p['n_worse']}")
    E.write("ener_tail", dict(
        what="tail-restricted discrimination and the filtering/ordering decomposition",
        qs=QS, gap_band=list(GAP_BAND), per_target=rows), n_expected=len(rows))
    return rows


if __name__ == "__main__":
    main()
