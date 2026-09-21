"""LANE L -- lane V's hull-capacity control, re-run at 40-60 residues.

THE CLAIM THIS TESTS.  Lane V (`s32_V_hull_capacity.py`) showed that at 9-16 residues the
ORACLE convex hull of 500 fragments is a statement about the EXPRESSIVE CAPACITY of
fragment space, not about retrieval: 500 windows drawn from a DIFFERENT target's universe
reach 1.1626 A against the retrieved pool's 1.1167 A -- NOT A RESULT at 0.49x MDE.  The
mechanism is dimensional.  The target lives in R^{3n}; at n = 13 that is **39 dimensions**,
and a convex hull of 500 points in 39 dimensions is a very large object, so an ORACLE fit of
500 non-negative weights has enough freedom to reach almost anything.

**THE PREDICTION.**  At n ~ 55, 3n ~ **164**.  500 points cannot fill 164 dimensions the way
they fill 39.  So the donor control must FAIL at long length: DONOR500 should be far worse
than BLOSUM500, and the hull should go back to being a statement about retrieval.

If it does, "the pool contains the answer" is capacity at peptide length and retrieval at
protein length -- and lane V's result is length-scoped rather than general.

ARMS, matched exactly to lane V's (same solver, same rounds, same rho, same frame seeding,
same alternation), differing only in the instrument:

  A  BLOSUM500   this target's shipped K=500 BLOSUM bank
  B  RAND500     500 windows drawn UNIFORMLY from this target's own leakage-safe corpus
                 -- retrieval blind, target aware
  C  DONOR500    500 windows drawn uniformly from a DIFFERENT long40 target's leakage-safe
                 corpus -- retrieval blind AND target blind; pure capacity

Every arm additionally excludes the target's own chain, at long length as at short, so no
arm can win by containing the answer verbatim.

ORACLE / NOT DEPLOYABLE throughout (the native sets the weights).
BASIS: **CA point cloud** (an unprojected convex combination). No chain claim is made.

    python -m s32.s32_L_hull_capacity run
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
ROWS = os.path.join(RESULTS, "L5_hull_capacity_rows.jsonl")

K = 500
M = 75
ROUNDS = 5                    # s29_O_ladder.HULL_ROUNDS, as lane V used
RHO = 1e3                     # s29_O_ladder.RHO_SUM, as lane V used
SEED = 32_0032_78


def convex_nnls(A, y, rho=RHO):
    """min |A w - y| over w >= 0, sum w = 1 -- s29_O_ladder.convex_nnls verbatim."""
    from scipy.optimize import nnls
    A = np.asarray(A, float); y = np.asarray(y, float).ravel()
    Aa = np.vstack([A, rho * np.ones((1, A.shape[1]))])
    ya = np.concatenate([y, [rho]])
    w, _ = nnls(Aa, ya, maxiter=50 * max(A.shape))
    s = w.sum()
    return np.full(A.shape[1], 1.0 / A.shape[1]) if s <= 0 else w / s


def hull(W, nat, n, rounds=ROUNDS, ref=None):
    """ORACLE best convex combination -- lane V's `hull`, unchanged."""
    from s12 import instrument as I
    Wf = np.asarray(W, float).reshape(len(W), 3 * n)
    seed = W[0] if ref is None else np.asarray(ref, float)
    natp = I.superpose_batch(np.asarray(nat, float)[None], seed)[0].ravel()
    w = None
    for _ in range(rounds):
        w = convex_nnls(Wf.T, natp)
        X = (w @ Wf).reshape(n, 3)
        natp = I.superpose_batch(np.asarray(nat, float)[None], X)[0].ravel()
    X = (w @ Wf).reshape(n, 3)
    return float(I.ca_rmsd(X, nat)), int((w > 1e-9).sum())


def uniform_windows(corpus, n, k, rng, exclude_seq, exclude_pdb, kmer=9):
    """`k` length-n windows drawn UNIFORMLY from the corpus, leakage-filtered.

    The long-length universe is ~1.7 M windows and is never materialised (see
    `s32_L_corpus`).  Chains are sampled in proportion to how many valid windows they carry,
    so the draw is uniform over WINDOWS and not over chains -- sampling chains uniformly
    would over-weight short chains and would not be the same control.
    """
    tk = {exclude_seq[i:i + kmer] for i in range(len(exclude_seq) - kmer + 1)}
    ok_chains, counts = [], []
    for a in range(len(corpus["pdb"])):
        if corpus["pdb"][a] in exclude_pdb:
            continue
        seq = corpus["seq"][a]
        if len(seq) <= n + 1:
            continue
        if any(x in seq for x in tk):
            continue
        ok = corpus["ok"][a]; fin = corpus["fin"][a]
        cok = np.concatenate([[0.0], np.cumsum(~ok)])
        cfin = np.concatenate([[0.0], np.cumsum(~fin)])
        s0 = np.arange(1, len(seq) - n)
        good = ((cok[s0 + n - 1] - cok[s0]) == 0) & ((cfin[s0 + n] - cfin[s0]) == 0)
        c = int(good.sum())
        if c:
            ok_chains.append((a, s0[good])); counts.append(c)
    if not ok_chains:
        return None
    counts = np.asarray(counts, float)
    p = counts / counts.sum()
    pick = rng.choice(len(ok_chains), size=k, replace=True, p=p)
    out = []
    for a in pick:
        ci, starts = ok_chains[a]
        s = int(starts[rng.integers(0, len(starts))])
        out.append(np.asarray(corpus["ca"][ci], float)[s:s + n])
    return np.stack(out)


def run(verbose=True):
    from core import geometry as geo
    from s12 import instrument as I
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    corpus = cp.load()
    tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
    by_pdb = {t["pdb"]: t for t in tg}
    done = set()
    if os.path.exists(ROWS):
        for line in open(ROWS):
            line = line.strip()
            if line:
                done.add(json.loads(line)["pdb"])
    rng = np.random.default_rng(SEED)
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        t0 = time.time()
        n = int(t["n"])
        W, _, _ = LD._bank_long(t)
        p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
        if not os.path.exists(p):
            p = glob.glob(os.path.join(BASE, "prots", t["pdb"].lower() + ".pdb"))[0]
        _, coords, _, _ = geo.native_coords_from_pdb(p)
        nat = np.asarray(coords["CA"], float)
        ref = W[:M][int(I.medoid(I.pairwise_rmsd(W[:M])))]

        rA, sA = hull(I.superpose_batch(W, ref), nat, n, ref=ref)
        WB = uniform_windows(corpus, n, K, rng, t["seq"], {t["pdb"].upper()})
        rB, sB = hull(I.superpose_batch(WB, ref), nat, n, ref=ref)
        # DONOR: a different long40 target supplies the leakage filter and the sequence
        # blindness; this target's own chain is STILL excluded so no arm can contain the
        # answer verbatim.
        donors = [q["pdb"] for q in tg if q["pdb"] != t["pdb"]]
        dp = donors[int(rng.integers(0, len(donors)))]
        WC = uniform_windows(corpus, n, K, rng, by_pdb[dp]["seq"],
                             {t["pdb"].upper(), dp.upper()})
        rC, sC = hull(I.superpose_batch(WC, ref), nat, n, ref=ref)

        row = dict(pdb=t["pdb"], n=n, fold=int(t["fold"]), three_n=3 * n, donor=dp,
                   blosum500=rA, blosum500_support=sA,
                   rand500=rB, rand500_support=sB,
                   donor500=rC, donor500_support=sC,
                   best_member_blosum=float(I.kabsch_rmsd_batch(W, nat).min()),
                   secs=time.time() - t0)
        with open(ROWS, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        if verbose:
            print(f"  {t['pdb']:6} n={n} 3n={3*n}  BLOSUM500 {rA:.4f} (supp {sA}) | "
                  f"RAND500 {rB:.4f} (supp {sB}) | DONOR500 {rC:.4f} (supp {sC})  "
                  f"{row['secs']:.1f}s", flush=True)
    return analyse(verbose)


def analyse(verbose=True):
    from s24 import stats_lib as ST
    rows = {}
    for line in open(ROWS):
        line = line.strip()
        if line:
            r = json.loads(line); rows[r["pdb"]] = r
    rows = [rows[k] for k in sorted(rows)]
    folds = np.array([r["fold"] for r in rows])
    g = {q: np.array([r[q] for r in rows], float)
         for q in ("blosum500", "rand500", "donor500", "best_member_blosum")}
    out = dict(n=len(rows), basis="CA point cloud; ORACLE / NOT DEPLOYABLE",
               mean_len=float(np.mean([r["n"] for r in rows])),
               mean_three_n=float(np.mean([r["three_n"] for r in rows])),
               arms={q: dict(mean=float(v.mean()), median=float(np.median(v)),
                             se=float(v.std(ddof=1) / len(v) ** 0.5))
                     for q, v in g.items()},
               support={q: float(np.mean([r[q + "_support"] for r in rows]))
                        for q in ("blosum500", "rand500", "donor500")})
    cmp = {}
    for lab, a, b in [("DONOR500 - BLOSUM500", "donor500", "blosum500"),
                      ("RAND500 - BLOSUM500", "rand500", "blosum500")]:
        c = ST.compare(g[a], g[b], folds=folds, label=f"long hull: {lab}")
        cmp[lab] = {kk: c[kk] for kk in
                    ("n", "mean_a", "mean_b", "effect", "median_effect", "se", "mde",
                     "effect_over_mde", "n_better", "n_worse")}
        if "per_fold" in c:
            cmp[lab]["per_fold"] = c["per_fold"]
    out["comparisons"] = cmp
    d = cmp["DONOR500 - BLOSUM500"]
    out["verdict"] = (
        "DONOR CONTROL FAILS AT LENGTH -- the hull is retrieval, not capacity"
        if d["effect"] > 0 and abs(d["effect_over_mde"]) >= 1.0 else
        "NOT MEASURED at this n" if abs(d["effect_over_mde"]) >= 0.7 else
        "NOT A RESULT -- the donor control SURVIVES at length, as at 9-16 residues")
    out["rows"] = rows
    path = os.path.join(RESULTS, "L5_hull_capacity.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(f"\n== long hull capacity  n={out['n']}  L~{out['mean_len']:.1f}  "
              f"3n~{out['mean_three_n']:.0f}   BASIS: CA CLOUD, ORACLE / NOT DEPLOYABLE")
        for q, v in out["arms"].items():
            print(f"   {q:<20}{v['mean']:8.4f}  med {v['median']:7.4f}  SE {v['se']:.4f}")
        print(f"   support (nonzero weights): {out['support']}")
        for lab, c in cmp.items():
            band = ("RESULT" if abs(c["effect_over_mde"]) >= 1 else
                    "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
            print(f"   {lab:<24}{c['effect']:+8.4f}  MDE {c['mde']:.4f}  "
                  f"{abs(c['effect_over_mde']):5.2f}x  W/L {c['n_better']}/{c['n_worse']}"
                  f"  {band}")
        print(f"   VERDICT: {out['verdict']}")
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    {"run": run, "analyse": analyse}[sys.argv[1] if len(sys.argv) > 1 else "run"]()
