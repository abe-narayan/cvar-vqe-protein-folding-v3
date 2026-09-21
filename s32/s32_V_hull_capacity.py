"""S32 LANE V -- is "the pool contains the answer" a fact about THIS pool, or about any 500
fragments of the right length?

Contract rule 17 and REPORT section 2 rest on one number:

    best sparse convex combination, K=500, s=10      1.1139     "2.10 A of headroom"

and conclude *the pool is not the bottleneck; everything downstream of retrieval destroys
2.10 A that is already present.*  The S29 ladder also shows the full hull at **1.1167**, i.e.
the sparsity constraint costs nothing -- so the number is a property of the **convex hull of
500 fragments in R^{3n}**, with 3n ~ 39.  A hull of 500 points in 39 dimensions is a large
object, and an ORACLE fit of 500 non-negative weights summing to one against a 39-dimensional
target has ~39 effective degrees of freedom to play with.

**THE CONTROL NOBODY RAN.**  If a hull of 500 windows chosen WITHOUT looking at the sequence
reaches the same distance from the native, then "the pool contains the answer" is a statement
about the expressive capacity of fragment space at this length, NOT about retrieval -- and the
2.10 A is not headroom that retrieval earned.  Charter section 27: *do not call a bad candidate
pool a selection problem just because the best candidate exists somewhere inside it.*

Three size-matched arms, same operator, same frame construction, same solver, same rounds:

  A  BLOSUM500   the shipped pool (reproduces the ladder's hull number)
  B  RAND500     500 windows drawn uniformly from the SAME target's own universe -- retrieval
                 blind, target aware
  C  DONOR500    500 windows drawn from a DIFFERENT target's universe at the same length --
                 retrieval blind AND target blind; pure capacity

Everything here is ORACLE / NOT DEPLOYABLE (the native sets the weights).
Basis: CA point cloud (an unprojected combination).  Per contract rule 16 the projection price
for a sparse combination of this class is ~0, but that is NOT assumed -- no chain claim is made.
"""
from __future__ import annotations
import argparse, glob, json, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s24 import stats_lib as ST                                          # noqa: E402
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s32_V_hull_capacity_rows.jsonl")
K = 500
ROUNDS = 5                    # s29_O_ladder.HULL_ROUNDS
RHO = 1e3                     # s29_O_ladder.RHO_SUM
SEED = 32_0032_77


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
    """ORACLE best convex combination of the posed windows W (k,n,3), native pose re-solved
    each round (s29's alternating `cf` problem).  `ref` is the medoid frame, used to seed the
    native's pose EXACTLY as `s29_O_ladder.cloud_row` does (`natp0 = superpose_one(nat,
    frame.ref)`), so arm A reproduces the published hull number rather than a different local
    optimum of the same alternation.  Returns (rmsd, n_support)."""
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


def medoid_ref(u, top75):
    W = u["W"]
    return W[top75][int(I.medoid(I.pairwise_rmsd(W[top75])))]


def posed(u, idx, ref):
    """All of `idx` posed on `ref` -- the frame the deployed average uses."""
    return I.superpose_batch(u["W"][np.asarray(idx, int)], np.asarray(ref, float))


def main(shard=None, n_shards=1):
    tg = I.targets()
    by_n = {}
    for t in tg:
        by_n.setdefault(t["n"], []).append(t["pdb"])
    done = set()
    if os.path.exists(ROWS):
        for line in open(ROWS):
            line = line.strip()
            if line:
                done.add(json.loads(line)["pdb"])
    rng = np.random.default_rng(SEED)
    for k, t in enumerate(tg):
        pdb, n = t["pdb"], t["n"]
        if pdb in done or (shard is not None and k % n_shards != shard):
            continue
        t0 = time.time()
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        pool = I.pool_idx(u, K)
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        top75 = pool[sub]
        nw = len(u["rr"])
        r = dict(pdb=pdb, n=n, fold=t["fold"], n_windows=int(nw))

        ref = medoid_ref(u, top75)
        # A -- the shipped BLOSUM 500
        rA, sA = hull(posed(u, pool, ref), nat, n, ref=ref)
        # B -- 500 uniformly at random from the SAME universe
        idxB = rng.choice(nw, size=min(K, nw), replace=False)
        rB, sB = hull(posed(u, idxB, ref), nat, n, ref=ref)
        # C -- 500 from a DIFFERENT target's universe at the same length
        donors = [q for q in by_n[n] if q != pdb]
        if donors:
            dp = donors[int(rng.integers(0, len(donors)))]
            ud = I.load_univ(dp)
            idxC = rng.choice(len(ud["rr"]), size=min(K, len(ud["rr"])), replace=False)
            Wd = ud["W"][idxC]
            #: pose the DONOR windows on THIS target's medoid frame, same operator
            rC, sC = hull(I.superpose_batch(Wd, ref), nat, n, ref=ref)
            r["donor"] = dp
        else:
            rC, sC = float("nan"), 0
            r["donor"] = None
        r.update(blosum500=rA, blosum500_support=sA, rand500=rB, rand500_support=sB,
                 donor500=rC, donor500_support=sC,
                 best_member_blosum=float(u["rr"][pool].min()),
                 best_member_rand=float(u["rr"][idxB].min()), secs=time.time() - t0)
        with open(ROWS if shard is None else ROWS.replace(".jsonl", "_s%d.jsonl" % shard), "a") as fh:
            fh.write(json.dumps(r) + "\n")
        print("  %-6s n=%2d  BLOSUM500 %.4f (supp %d) | RAND500 %.4f (supp %d) | DONOR500 %.4f "
              "(supp %d)  %.1fs" % (pdb, n, rA, sA, rB, sB, rC, sC, r["secs"]), flush=True)


def analyse():
    rows = {}
    for f in sorted(glob.glob(ROWS.replace(".jsonl", "*.jsonl"))):
        for line in open(f):
            line = line.strip()
            if line:
                q = json.loads(line)
                rows[q["pdb"]] = q
    pdbs = sorted(rows)
    G = lambda k: np.array([rows[p][k] for p in pdbs], float)
    folds = np.array([rows[p]["fold"] for p in pdbs], int)
    out = dict(n=len(pdbs), basis="CA point cloud (an unprojected convex combination)",
               oracle="ORACLE / NOT DEPLOYABLE", rounds=ROUNDS, K=K, seed=SEED,
               provenance=ST.provenance(__file__), per_target=rows)
    print("=" * 100)
    print("HULL CAPACITY: is 1.11 A a property of THIS pool or of 500 fragments?  n=%d  ORACLE"
          % len(pdbs))
    for k, lab in (("blosum500", "A  BLOSUM500 (the shipped pool)"),
                   ("rand500", "B  RAND500  (same universe, retrieval-blind)"),
                   ("donor500", "C  DONOR500 (another target's universe, target-blind)")):
        v = G(k)
        out[k] = dict(mean=float(np.nanmean(v)), median=float(np.nanmedian(v)),
                      support_mean=float(np.nanmean(G(k + "_support"))))
        print("  %-44s %.4f   (median %.4f, support %.1f)"
              % (lab, np.nanmean(v), np.nanmedian(v), np.nanmean(G(k + "_support"))))
    for a, b in (("rand500", "blosum500"), ("donor500", "blosum500")):
        m = np.isfinite(G(a)) & np.isfinite(G(b))
        c = ST.compare(G(a)[m], G(b)[m], folds[m], names=list(np.array(pdbs)[m]),
                       label="%s vs %s (CA POINT CLOUD, ORACLE hull)" % (a, b))
        out["%s_vs_%s" % (a, b)] = c
        gate = ("RESULT" if abs(c["effect_over_mde"]) >= 1.0 and c["folds_same_sign"] >= 4
                and min(c["ci95_fold"]) * max(c["ci95_fold"]) > 0
                else "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
        print("  %-20s vs %-12s effect %+0.4f  MDE %.4f  %.2fx  folds %d/5  %dW/%dL  %s"
              % (a, b, c["effect"], c["mde"], c["effect_over_mde"], c["folds_same_sign"],
                 c["n_better"], c["n_worse"], gate))
    print("  READ: if RAND500/DONOR500 land near BLOSUM500, 'the pool contains the answer' is a "
          "statement about\n        fragment-space capacity at 9-16 residues, not about retrieval, "
          "and the 2.10 A is not headroom\n        retrieval earned.")
    print("=" * 100)
    with open(os.path.join(RESULTS, "s32_V_hull_capacity.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    analyse() if a.analyse else main(a.shard, a.n_shards)
