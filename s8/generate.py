"""Can a better-DISTRIBUTED pool move selected RMSD, when the selector is powerless?

The inverted objective
Sprint 7 finding 6 measured that selected CA-RMSD tracks the pool MEAN and moves AGAINST
the pool BEST: from K=25 to K=2000 the pool best improved 2.350 -> 1.504 A while selected
got WORSE, 3.439 -> 3.520, and the pool mean rose 4.177 -> 4.596.  A near-powerless
selector returns approximately a fixed QUANTILE of the pool's own RMSD distribution.

Every previous generation effort maximised pool BEST.  For this pipeline that is the wrong
objective.  If the pool's mean and spread were driven far enough down, selected RMSD would
follow even under a selector barely better than chance.  This module builds pools that
target the distribution instead of the minimum and measures the whole chain:

    pool best -> pool mean -> pool SD -> in-band Spearman -> SELECTED RMSD

with an oracle-ranked column throughout as the ceiling of each pool.

The warning this respects
Sprint 7 also tried native-free PRUNING of a fixed pool: the best rule improved the pool
mean 4.232 -> 3.712 and moved selected only 3.399 -> 3.420.  Pruning cannot create better
candidates and typicality-based pruning eats the ranker's own skill.  So the arms here are
GENERATION arms -- retrieval over the FULL window universe under a different key, torsion
recombination, and torsion-space contraction -- not subsets of the shipped top-500.  A
pool-mean gain that does not reach selected RMSD is reported as a failure, not a result.

What A candidate is
The shipped generator scores every length-n CA window of (out-of-fold peptides + this
fold's fragment library), sorts by BLOSUM62 similarity to the target sequence, and keeps
the top 2000; K=500 is the pool every sprint-7 number reports on.  The universe is only
~13,000 windows, so this module caches ALL of it per target and can retrieve under any
key.  A window carries its parent's real (phi, psi), so torsion-space arms are exact.

DEPLOYABLE vs DIAGNOSTIC -- the split is structural
Every generator lives in `GENERATORS` and is called as `f(view, pred, rng)`.  `view` is
built by `deployable_view()`, which contains the window coordinates, torsions, encoded
sequences, provenance flags and BLOSUM scores -- and NOTHING derived from the target's
native structure.  `pred` is the shipped fold-model prediction for the target sequence.
`test_generate.py` asserts the split by NaN-poisoning every native array and checking each
generator's output is bit-identical.

Native coordinates are read in exactly two places, both after a pool is fixed:
`_chain()` (to report the CA-RMSD of what was already chosen) and `_oracle_select()` /
the `oracleK` arm, which are labelled DIAGNOSTIC and feed nothing.

PROTOCOL
126-target tuning instrument (`s7.debias.tuning_targets`), cluster-disjoint from both the
24-target dev set and the 60-target benchmark.  SE on the mean is 0.147 A.  The 24-target
dev set is reserved for a tuning winner, and this sprint produced none: see RESULT.

    python -m s8.generate univ      # cache the full window universe + Ramachandran tables
    python -m s8.generate diag      # what the selector actually does to a distribution
    python -m s8.generate transfer  # DIAGNOSTIC: what filter SKILL is worth
    python -m s8.generate gen       # every generator arm, full chain to selected RMSD
    python -m s8.generate soft      # sampled rather than truncated quality filters
    python -m s8.generate law       # one law from pool distribution to the answer
    python -m s8.generate report    # print everything already computed

RESULT -- the hypothesis is right and unusable
1. THE LEVER IS REAL.  `stage_transfer` drives a synthetic filter of controlled rank
   correlation against true CA-RMSD across the whole universe.  Selected RMSD falls
   monotonically with it: 3.531 (rho 0.00) -> 3.410 (0.30) -> 3.334 (0.60) -> 3.003
   (0.91) -> 2.406 (perfect), at +0.397 A of selected per A of pool mean, r = +0.872.
   Significance arrives at rho ~ 0.60.  Driving the distribution down DOES move the
   answer.

2. NO NATIVE-FREE KEY PULLS IT.  BLOSUM's own rank correlation with true RMSD over the
   universe is +0.066.  Ramachandran plausibility reaches +0.307, rg agreement +0.321,
   typicality +0.410, the shipped score itself +0.587 -- and every one of them selects
   WORSE than BLOSUM, by +0.11 to +0.37 A.  The mechanism is one column: at matched rank
   correlation the synthetic filter leaves pool best at 1.65 A while the real keys leave
   2.16-3.40 A.  A real key's error is correlated with the property that makes a candidate
   good, so truncating by it deletes the near-native band along with the tail.  Sampling
   instead of truncating (`stage_soft`, Gumbel-top-k at five temperatures on three keys)
   trades the two off exactly as predicted and never wins.

3. 41 DEPLOYABLE ARMS, ZERO WINNERS.  Quality-filtered retrieval, torsion recombination,
   torsion contraction, cluster-centroid contraction, peptide-context matching, soft
   sampling: best is `contract0.5` at -0.021 A [-0.157, +0.114], inside the SE.  The 14
   arms that improve the pool mean move selected +0.054 A, the wrong way.

4. WHY CONTRACTION IN PARTICULAR CANNOT WORK.  `stage_diag`: the pool is spread 4.142 A
   wide (mean pairwise CA-RMSD) about a medoid that is 3.665 A from the native, while the
   weak selector already returns 3.454 A.  Contracting to zero spread lands ON the medoid,
   which is WORSE than what the selector already achieves.  Shrinking the spread is
   strictly harmful unless the pool is re-centred, and re-centring is recognition.

5. THE LAW.  Over all 55 arms measured here -- deployable, synthetic and oracle alike, on
   identical targets with an identical selector -- mean selected RMSD is

       sel = 0.257 * pool_best + 0.311 * pool_mean + 1.673      R^2 0.79, RMSE 0.108 A

   The 1.67 A constant is what no pool improvement touches.  Even the ORACLE best-500
   subset of this library (pool mean 2.350, sd 0.256, 41% of candidates under 2 A) selects
   only 2.406 A, and its selector percentile degrades to 59.9 -- worse than a coin flip.
   So sub-2.4 A is unreachable by pool construction on this library at any filter skill,
   and the pool-distribution axis is closed with its ceiling priced rather than guessed.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np                                                          # noqa: E402
from scipy.stats import spearmanr                                           # noqa: E402

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:
    pass

import distogram as dgm                                                     # noqa: E402
import peptide_db as db                                                     # noqa: E402
import protein_geometry as geo                                              # noqa: E402
from s7 import audit, debias, poolsize                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
UCACHE = os.path.join(HERE, "generate_univ")      #: full window universe, per target
RAMA = os.path.join(HERE, "generate_rama.npz")    #: per-fold (aa, phi, psi) log-densities

K = 500                       #: the pool size every sprint-7 number reports on
BAND = 1.5                    #: in-band = within pool_best + BAND
MIN_BAND = 5
GRID = poolsize.GRID          #: the shipped Bayes-risk grid
CA_CA = 3.80                  #: virtual CA-CA bond, for the rg reconstruction
RB = 36                       #: Ramachandran bins per axis (10 degrees)
SEED = 20260904
ALPHAS = (0.85, 0.7, 0.5, 0.3)        #: torsion-contraction strengths
KEEPS = (0.1, 0.25, 0.5)              #: quality-filter keep fractions of the universe


# ============================================================ stage: universe + rama
def _windows_all(pool, n, is_pep):
    """Every length-n window of every pool member: CA, phi, psi, encoded seq, origin.

    Same iteration order as `audit.windows_of`, which is what pins the BLOSUM argsort and
    therefore reproduces the shipped pool exactly when truncated to the top 2000.
    """
    cas, phis, psis, seqs, org = [], [], [], [], []
    for q, flag in zip(pool, is_pep):
        m = len(q.seq)
        if m < n:
            continue
        e = audit.encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n]); phis.append(q.phi[s:s + n])
            psis.append(q.psi[s:s + n]); seqs.append(e[s:s + n]); org.append(flag)
    return (np.stack(cas), np.stack(phis), np.stack(psis), np.stack(seqs),
            np.array(org, bool))


def _upath(pdbid):
    return os.path.join(UCACHE, f"{pdbid}.npz")


def stage_univ(verbose=True):
    """Cache the FULL window universe per target, plus per-fold Ramachandran tables.

    Verified against `s7/debias_cache`: the BLOSUM top-2000 slice of this universe must
    reproduce the shipped pool's `D` and `rr` to float32, so every arm below is measured
    on provably the same instrument as sprint 7.
    """
    os.makedirs(UCACHE, exist_ok=True)
    folds = db.folds(5)
    memo = {}
    _orig = dgm._fold_fragments

    def frags(fold, n_folds, threshold=db.IDENTITY_THRESHOLD):
        key = (fold, n_folds, threshold)
        if key not in memo:
            memo[key] = _orig(fold, n_folds, threshold)
        return memo[key]

    dgm._fold_fragments = frags
    rama = np.zeros((5, 20, RB, RB), np.float64)
    need_rama = not os.path.exists(RAMA)
    worst_d = worst_r = 0.0
    try:
        tg = debias.tuning_targets()
        todo = [p for p in tg if not os.path.exists(_upath(p.pdb))]
        if verbose:
            print(f"{len(todo)}/{len(tg)} universes to build "
                  f"(rama {'needed' if need_rama else 'cached'})", flush=True)
        if need_rama:
            for fold in range(5):
                t0 = time.time()
                _accum_rama(rama[fold],
                            [q for q in db.load() if folds[q.seq] != fold]
                            + list(frags(fold, 5)))
                print(f"  rama fold {fold}: {rama[fold].sum():.0f} residues "
                      f"({time.time()-t0:.0f}s)", flush=True)
            np.savez_compressed(RAMA, cnt=rama)
            print(f"wrote {RAMA}", flush=True)
        for k, p in enumerate(sorted(todo, key=lambda x: folds[x.seq])):
            fold = folds[p.seq]
            peps = [q for q in db.load() if folds[q.seq] != fold and q.seq != p.seq]
            fr = list(frags(fold, 5))
            t0 = time.time()
            W, PH, PS, S, org = _windows_all(peps + fr, p.n,
                                             [True] * len(peps) + [False] * len(fr))
            sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
            order = np.argsort(-sim, kind="stable")
            i, j = audit.pair_index(p.n)
            z = np.load(os.path.join(debias.CACHE, f"{p.pdb}.npz"), allow_pickle=True)
            top = order[:audit.KMAX]
            ed = float(np.abs(audit.pair_dists(W[top], i, j) - z["D"]).max())
            er = float(np.abs(audit.kabsch_rmsd_batch(W[top], p.ca) - z["rr"]).max())
            assert ed < 1e-3 and er < 1e-3, f"{p.pdb}: universe mismatch {ed} {er}"
            worst_d, worst_r = max(worst_d, ed), max(worst_r, er)
            tmp = _upath(p.pdb) + ".tmp.npz"
            np.savez_compressed(
                tmp, pdb=p.pdb, n=p.n, fold=fold, seq=p.seq,
                W=W.astype(np.float32), PHI=PH.astype(np.float16),
                PSI=PS.astype(np.float16), S=S.astype(np.int8), org=org,
                sim=sim.astype(np.float32), order=order.astype(np.int32),
                rr=audit.kabsch_rmsd_batch(W, p.ca).astype(np.float32),
                nat_ca=np.asarray(p.ca, np.float32))
            os.replace(tmp, _upath(p.pdb))
            if verbose:
                print(f"[{k+1:3d}/{len(todo)}] {p.pdb:6} n={p.n:2d} fold={fold} "
                      f"windows={len(W):6d} dD={ed:.1e} ({time.time()-t0:.0f}s)",
                      flush=True)
    finally:
        dgm._fold_fragments = _orig
    if verbose:
        print(f"worst |dD|={worst_d:.1e}  worst |drr|={worst_r:.1e}", flush=True)


def _accum_rama(tab, entries):
    """Counts into a (20, RB, RB) per-amino-acid (phi, psi) table. Out-of-fold data only."""
    for e in entries:
        a = audit.encode(e.seq)
        bp = np.clip(((np.asarray(e.phi) + math.pi) / (2 * math.pi) * RB).astype(int),
                     0, RB - 1)
        bs = np.clip(((np.asarray(e.psi) + math.pi) / (2 * math.pi) * RB).astype(int),
                     0, RB - 1)
        np.add.at(tab, (a, bp, bs), 1.0)


_RAMA = [None]


def rama_logp(fold):
    """(20, RB, RB) log P(phi, psi | aa) from fold `fold`'s own TRAINING corpus."""
    if _RAMA[0] is None:
        c = np.load(RAMA)["cnt"] + 0.5
        _RAMA[0] = np.log(c / c.sum((2, 3), keepdims=True))
    return _RAMA[0][fold]


def load_univ(pdbid):
    z = np.load(_upath(pdbid), allow_pickle=True)
    return {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
            "seq": str(z["seq"]), "W": z["W"].astype(np.float64),
            "PHI": z["PHI"].astype(np.float64), "PSI": z["PSI"].astype(np.float64),
            "S": z["S"].astype(np.int64), "org": z["org"], "sim": z["sim"].astype(float),
            "order": z["order"].astype(int), "rr": z["rr"].astype(float),
            "nat_ca": z["nat_ca"].astype(np.float64)}


def cached_targets():
    if not os.path.isdir(UCACHE):
        return []
    return sorted(f[:-4] for f in os.listdir(UCACHE) if f.endswith(".npz"))


# ============================================================ DEPLOYABLE path
# Nothing between here and the DIAGNOSTIC banner may touch a native quantity.
def deployable_view(u):
    """The subset of a cached universe a generator is allowed to see.

    `rr` and `nat_ca` are deliberately absent, so a generator cannot read a native even by
    accident; `test_generate.py` NaN-poisons them and asserts bit-identical output.
    """
    return {k: u[k] for k in ("pdb", "n", "fold", "seq", "W", "PHI", "PSI", "S",
                              "org", "sim", "order")}


def pair_D(W, n):
    i, j = audit.pair_index(n)
    return audit.pair_dists(W, i, j)


def shipped_score(pred, D):
    """The shipped Bayes-risk score, evaluated straight off pair distances."""
    return debias.score_risk(pred["risk"], GRID, np.asarray(D, float))


def _wrap(a):
    return (a + math.pi) % (2 * math.pi) - math.pi


def _circmean(A, axis=0):
    return np.arctan2(np.sin(A).mean(axis), np.cos(A).mean(axis))


def rebuild(PHI, PSI):
    """Ideal-geometry CA traces from torsions. (B, n) radians -> (B, n, 3)."""
    return geo.build_backbone_batch(np.asarray(PHI, float), np.asarray(PSI, float))["CA"]


def rg_of(D, n):
    """Radius of gyration from min_sep=2 pair distances alone.

    rg^2 = (1/(2 n^2)) sum_{i,j} d_ij^2.  The cached pair set omits |i-j| < 2; |i-j| = 0
    contributes nothing and |i-j| = 1 is the virtual CA-CA bond, near-constant at 3.80 A,
    so both missing blocks are supplied analytically and no coordinate is needed.
    """
    D = np.asarray(D, float)
    s = 2.0 * (D ** 2).sum(-1) + 2.0 * (n - 1) * CA_CA ** 2
    return np.sqrt(s / (2.0 * n ** 2))


def predicted_rg(pred, n):
    """The rg implied by the PREDICTED distance matrix. Sequence-derived, native-free."""
    return float(rg_of(pred["expected"][None, :], n)[0])


def rama_score(view):
    """Mean log P(phi_i, psi_i | aa_i) per window, under the fold's own training corpus.

    The target's amino acid at each position is used, not the parent window's, so this
    asks 'is this conformation plausible FOR THIS SEQUENCE' -- the native-free
    plausibility filter the sprint brief calls for.
    """
    L = rama_logp(view["fold"])
    a = audit.encode(view["seq"])
    bp = np.clip(((view["PHI"] + math.pi) / (2 * math.pi) * RB).astype(int), 0, RB - 1)
    bs = np.clip(((view["PSI"] + math.pi) / (2 * math.pi) * RB).astype(int), 0, RB - 1)
    return L[a[None, :], bp, bs].mean(1)


def _blosum_top(view, m):
    return view["order"][:m]


def _keep_by(rank_key, cand, m):
    """The `m` candidates of `cand` with the smallest `rank_key` (already restricted)."""
    return cand[np.argsort(rank_key, kind="stable")[:m]]


# --------------------------------------------------------------- generator arms
# Each returns (W, tag) where W is (B, n, 3) CA traces.  `view` carries no native.
def gen_base(view, pred, rng, m=K):
    """The shipped generator: the top-`m` BLOSUM windows, real coordinates."""
    return view["W"][_blosum_top(view, m)]


def gen_universe(view, pred, rng, m=None):
    """No retrieval at all: the whole window universe. The zero-information control."""
    return view["W"]


def gen_rebuilt(view, pred, rng, m=K):
    """The shipped pool rebuilt from its own torsions -- the control for torsion arms."""
    idx = _blosum_top(view, m)
    return rebuild(view["PHI"][idx], view["PSI"][idx])


def _filtered(view, pred, key, keep, m=K, base=None):
    """Retrieve a wide BLOSUM band, then keep by a native-free quality key, then cut to m.

    This is GENERATION, not pruning: the band is `m / keep` wide, so the pool that comes
    out is still `m` candidates -- the filter changes WHICH `m`, it does not shrink the
    pool.  That is the distinction the sprint-7 pruning negative turns on.
    """
    wide = min(len(view["W"]), int(round(m / keep)))
    cand = _blosum_top(view, wide) if base is None else base
    return _keep_by(key[cand], cand, m)


def gen_rama(view, pred, rng, m=K, keep=0.25):
    """Quality-filtered retrieval by Ramachandran plausibility of the target's sequence."""
    return view["W"][_filtered(view, pred, -rama_score(view), keep, m)]


def gen_rg(view, pred, rng, m=K, keep=0.25):
    """Quality-filtered retrieval by agreement with the PREDICTED radius of gyration.

    rg is an average over all pairs, so the predictor's per-pair error largely cancels;
    it is a far more reliable functional of the same prediction than any single distance,
    and a candidate with the wrong compactness is exactly the upper tail of the RMSD
    distribution.  Native-free: rg_pred comes from the sequence model only.
    """
    n = view["n"]
    rgc = rg_of(pair_D(view["W"], n), n)
    return view["W"][_filtered(view, pred, np.abs(rgc - predicted_rg(pred, n)), keep, m)]


def gen_pep(view, pred, rng, m=K):
    """Context matching: retrieve ONLY from isolated peptides, never protein interiors."""
    cand = view["order"][view["org"][view["order"]]][:m]
    if len(cand) < 5:
        cand = _blosum_top(view, m)
    return view["W"][cand]


def gen_frag(view, pred, rng, m=K):
    """The complement of `pep`: protein-interior fragments only. A control."""
    cand = view["order"][~view["org"][view["order"]]][:m]
    if len(cand) < 5:
        cand = _blosum_top(view, m)
    return view["W"][cand]


def gen_rg_rama(view, pred, rng, m=K, keep=0.25):
    """Both native-free quality keys, applied in sequence over a wide retrieval band."""
    n = view["n"]
    wide = min(len(view["W"]), int(round(m / keep)))
    cand = _blosum_top(view, wide)
    half = max(m, len(cand) // 2)
    cand = _keep_by(-rama_score(view)[cand], cand, half)
    rgc = rg_of(pair_D(view["W"][cand], n), n)
    return view["W"][_keep_by(np.abs(rgc - predicted_rg(pred, n)), cand, m)]


def gen_recomb(view, pred, rng, m=K, parents=200, cuts=(2,)):
    """Torsion-space recombination of retrieved fragments.

    Composes each candidate from two or three retrieved windows rather than inheriting a
    whole window wholesale, so a single poor parent contributes only part of a candidate
    and the bad tail is not imported intact.  Torsions are spliced and the chain rebuilt
    with ideal geometry, which is why `rebuilt` is the control this must beat.
    """
    idx = _blosum_top(view, min(parents, len(view["W"])))
    n = view["n"]
    P, Q = view["PHI"][idx], view["PSI"][idx]
    B = len(idx)
    nc = rng.choice(cuts, size=m)
    ph = np.empty((m, n)); ps = np.empty((m, n))
    for c in set(cuts):
        sel = np.where(nc == c)[0]
        if not len(sel):
            continue
        pts = np.sort(rng.integers(2, n - 1, size=(len(sel), c)), axis=1)
        par = rng.integers(0, B, size=(len(sel), c + 1))
        seg = np.zeros((len(sel), n), int)
        for t in range(c):
            seg += (np.arange(n)[None, :] >= pts[:, t:t + 1])
        pick = np.take_along_axis(par, seg, axis=1)
        ph[sel] = P[pick, np.arange(n)[None, :]]
        ps[sel] = Q[pick, np.arange(n)[None, :]]
    return rebuild(ph, ps)


def gen_contract(view, pred, rng, m=K, alpha=0.5):
    """Torsion-space contraction toward the pool's own circular-mean torsion profile.

    The direct test of the sprint's hypothesis: shrink the pool's SPREAD while leaving its
    centre alone, and see whether selected RMSD follows the distribution down.  alpha = 1
    is `rebuilt` exactly; alpha = 0 collapses the pool to a single consensus structure.
    Uses only the pool's own candidates -- no native, no target-specific fit.
    """
    idx = _blosum_top(view, m)
    P, Q = view["PHI"][idx], view["PSI"][idx]
    mp, mq = _circmean(P), _circmean(Q)
    return rebuild(mp + alpha * _wrap(P - mp), mq + alpha * _wrap(Q - mq))


def _medoids(W, k, rng, cap=200):
    """k-medoids in CA-RMSD space over a capped subsample. Native-free by construction."""
    B = min(cap, len(W))
    sub = np.arange(B)
    M = np.stack([audit.kabsch_rmsd_batch(W[sub], W[b]) for b in sub])
    lab = np.zeros(B, int)
    med = [int(np.argmin(M.mean(1)))]
    for _ in range(k - 1):
        d = M[med].min(0)
        med.append(int(np.argmax(d)))
    for _ in range(12):
        lab = np.argmin(M[med], 0)
        new = []
        for c in range(len(med)):
            mm = np.where(lab == c)[0]
            new.append(med[c] if not len(mm) else int(mm[np.argmin(M[np.ix_(mm, mm)].mean(1))]))
        if new == med:
            break
        med = new
    return np.array(med), lab, M


def gen_centroid(view, pred, rng, m=K, k=3, alpha=0.5):
    """Contract toward the nearest of `k` cluster medoids rather than one global mean.

    The near-native band holds 2-3 populated clusters (`s8/repr_ceiling multi`), so a
    single circular mean averages across modes.  This tightens each mode separately.
    """
    idx = _blosum_top(view, m)
    P, Q = view["PHI"][idx], view["PSI"][idx]
    W = view["W"][idx]
    med, lab, M = _medoids(W, k, rng)
    full = np.argmin(np.stack([audit.kabsch_rmsd_batch(W, W[b]) for b in med]), 0)
    ph, ps = P.copy(), Q.copy()
    for c in range(len(med)):
        mm = np.where(full == c)[0]
        if not len(mm):
            continue
        mp, mq = _circmean(P[mm]), _circmean(Q[mm])
        ph[mm] = mp + alpha * _wrap(P[mm] - mp)
        ps[mm] = mq + alpha * _wrap(Q[mm] - mq)
    return rebuild(ph, ps)


def gen_recomb_rg(view, pred, rng, m=K, over=4):
    """Recombine `over * m` candidates, then keep the `m` whose rg matches the prediction.

    Generation followed by native-free quality filtering -- the two levers together.
    """
    W = gen_recomb(view, pred, rng, m=over * m)
    n = view["n"]
    rgc = rg_of(pair_D(W, n), n)
    return W[np.argsort(np.abs(rgc - predicted_rg(pred, n)), kind="stable")[:m]]


def _univ_key(view, pred, name):
    """A native-free ranking key over the WHOLE universe. Lower is better."""
    n = view["n"]
    if name == "blosum":
        return -view["sim"]
    if name == "rama":
        return -rama_score(view)
    D = pair_D(view["W"], n)
    if name == "rg":
        return np.abs(rg_of(D, n) - predicted_rg(pred, n))
    if name == "score":
        return shipped_score(pred, D)
    raise ValueError(name)


def gen_soft(view, pred, rng, m=K, key="score", T=1.0):
    """SOFT quality filtering: sample the pool from the key instead of truncating it.

    `stage_transfer` measures why every hard native-free filter fails. At matched rank
    correlation against true RMSD, a synthetic filter keeps a 1.65 A pool best while the
    real keys keep 2.2-3.4 A: the real keys' errors are correlated with the very property
    that makes a candidate good, so truncating removes the near-native band along with the
    tail. Sampling with weight exp(-z/T) instead of truncating keeps that band alive at
    reduced density while still shifting the mean, which is the only way a biased key can
    buy distribution without paying for it in the pool's best. T -> 0 recovers the hard
    filter, T -> inf recovers uniform sampling. Gumbel-top-k, so it is exact weighted
    sampling without replacement.
    """
    z = _univ_key(view, pred, key)
    z = (z - z.mean()) / max(z.std(), 1e-9)
    g = -z / T + rng.gumbel(size=len(z))
    return view["W"][np.argsort(-g)[:min(m, len(z))]]


GENERATORS = {
    "base":        (gen_base, {}),
    "universe":    (gen_universe, {}),
    "rebuilt":     (gen_rebuilt, {}),
    "pep":         (gen_pep, {}),
    "frag":        (gen_frag, {}),
    "recomb2":     (gen_recomb, {"cuts": (1,)}),
    "recomb3":     (gen_recomb, {"cuts": (2,)}),
    "recomb_mix":  (gen_recomb, {"cuts": (1, 2, 3)}),
    "recomb_rg":   (gen_recomb_rg, {}),
}
for _k in KEEPS:
    GENERATORS[f"rama{_k:g}"] = (gen_rama, {"keep": _k})
    GENERATORS[f"rg{_k:g}"] = (gen_rg, {"keep": _k})
    GENERATORS[f"rgrama{_k:g}"] = (gen_rg_rama, {"keep": _k})
for _a in ALPHAS:
    GENERATORS[f"contract{_a:g}"] = (gen_contract, {"alpha": _a})
    GENERATORS[f"centroid{_a:g}"] = (gen_centroid, {"alpha": _a})
SOFT_T = (0.15, 0.3, 0.6, 1.2, 2.5)
for _key in ("score", "rg", "blosum"):
    for _T in SOFT_T:
        GENERATORS[f"soft_{_key}{_T:g}"] = (gen_soft, {"key": _key, "T": _T})


# ============================================================ DIAGNOSTIC path
# Everything below reads a native.  None of it is deployable and none of it feeds back
# into pool construction or ranking.
def _oracle_select(D, Dnat):
    """DIAGNOSTIC upper bound: rank by agreement with the native distance matrix."""
    return int(np.argmin(np.abs(D - Dnat[None, :]).mean(1)))


def gen_oracleK(u, m=K):
    """DIAGNOSTIC ONLY. The `m` lowest-RMSD windows in the universe.

    The ceiling on any retrieval-based quality filter: no native-free key can select a
    better-distributed `m`-subset of this universe than this one.
    """
    return u["W"][np.argsort(u["rr"], kind="stable")[:m]]


def _chain(rr, score, Dnat=None, D=None):
    """REPORTING ONLY. The full chain for one target's pool, natives read after the fact."""
    sel = int(np.argmin(score))
    s = float(rr[sel])
    best = float(rr.min())
    band = rr <= best + BAND
    out = {"sel": s, "best": best, "mean": float(rr.mean()), "sd": float(rr.std(ddof=1)),
           "median": float(np.median(rr)),
           "q10": float(np.percentile(rr, 10)), "q30": float(np.percentile(rr, 30)),
           "q50": float(np.percentile(rr, 50)), "q90": float(np.percentile(rr, 90)),
           "pct": float(100.0 * (rr < s).mean()), "n": int(len(rr)),
           "f2": float((rr < 2.0).mean()), "f15": float((rr < 1.5).mean()),
           "n_band": int(band.sum()),
           "rho_global": float(spearmanr(score, rr).statistic),
           "rho_band": (float(spearmanr(score[band], rr[band]).statistic)
                        if band.sum() >= MIN_BAND else None)}
    if Dnat is not None:
        out["oracle_diagnostic"] = float(rr[_oracle_select(D, Dnat)])
    return out


def _agg(rows, key):
    v = np.array([r[key] for r in rows if r.get(key) is not None], float)
    if not len(v):
        return None, None
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))


def write_json(name, obj):
    p = os.path.join(HERE, name)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    os.replace(tmp, p)


def read_json(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


# ============================================================ stage: diag
def stage_diag():
    """What does the selector do to a DISTRIBUTION? The measurement that sets the target.

    For each target the shipped scorer picks a candidate at some percentile of the pool's
    own RMSD distribution.  If the selector were a fixed-quantile draw, selected RMSD
    would be exactly predictable from the pool's quantile function -- and then improving
    the distribution would be guaranteed to pay.  This measures how nearly that holds, and
    prices the two ceilings a quality filter faces: the best `K`-subset of the universe
    (`oracleK`) and the universe's own distribution.
    """
    tgs = cached_targets()
    print(f"distribution diagnostics on {len(tgs)} targets", flush=True)
    per = []
    for k, pdbid in enumerate(tgs):
        u = load_univ(pdbid)
        pred = poolsize.load_pred(pdbid)
        n = u["n"]
        Dnat = pair_D(u["nat_ca"][None], n)[0]
        idx = u["order"][:K]
        W, rr = u["W"][idx], u["rr"][idx]
        D = pair_D(W, n)
        sc = shipped_score(pred, D)
        row = {"pdb": pdbid, "n": n, "n_windows": int(len(u["W"]))}
        row["pool"] = _chain(rr, sc, Dnat, D)
        Du = pair_D(u["W"], n)
        row["univ"] = _chain(u["rr"], shipped_score(pred, Du), Dnat, Du)
        o = np.argsort(u["rr"], kind="stable")[:K]
        row["oracleK_diagnostic"] = _chain(u["rr"][o], shipped_score(pred, Du[o]),
                                           Dnat, Du[o])
        # geometry of the pool: how much of the mean is CENTRING error vs SPREAD?
        med, _, M = _medoids(W, 1, np.random.default_rng(SEED))
        row["spread"] = float(M.mean())
        row["medoid_rmsd"] = float(audit.kabsch_rmsd_batch(
            u["nat_ca"][None], W[int(med[0])])[0])
        row["rg_pred"] = predicted_rg(pred, n)
        row["rg_nat"] = float(rg_of(Dnat[None], n)[0])
        row["rg_pool_mean"] = float(rg_of(D, n).mean())
        per.append(row)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    out = {"n_targets": len(per), "K": K, "per_target": per}
    # Is selected RMSD a fixed-quantile draw?  Compare the mean selected RMSD to the mean
    # of the pool's own q-quantile, sweeping q, and report the q that reproduces it.
    sel = np.array([r["pool"]["sel"] for r in per])
    qs = np.arange(1, 100)
    QQ = np.array([[np.percentile(np.sort(x), q) for q in qs] for x in
                   [load_univ(r["pdb"])["rr"][load_univ(r["pdb"])["order"][:K]]
                    for r in per]])
    err = np.abs(QQ.mean(0) - sel.mean())
    out["quantile_law"] = {
        "mean_sel": float(sel.mean()),
        "matching_quantile": int(qs[int(np.argmin(err))]),
        "mean_pct": float(np.mean([r["pool"]["pct"] for r in per])),
        "sd_pct": float(np.std([r["pool"]["pct"] for r in per], ddof=1)),
        "corr_sel_q": float(np.corrcoef(sel, QQ[:, int(np.argmin(err))])[0, 1]),
        "corr_sel_mean": float(np.corrcoef(sel, [r["pool"]["mean"] for r in per])[0, 1]),
        "corr_sel_best": float(np.corrcoef(sel, [r["pool"]["best"] for r in per])[0, 1]),
        "corr_sel_sd": float(np.corrcoef(sel, [r["pool"]["sd"] for r in per])[0, 1])}
    for tag in ("pool", "univ", "oracleK_diagnostic"):
        rows = [r[tag] for r in per]
        out[tag] = {k2: _agg(rows, k2)[0] for k2 in
                    ("sel", "best", "mean", "sd", "q10", "q30", "median", "pct",
                     "f2", "f15", "rho_band", "oracle_diagnostic")}
        out[tag]["sel_se"] = _agg(rows, "sel")[1]
    out["geometry"] = {
        "pool_spread_mean_pairwise_rmsd": float(np.mean([r["spread"] for r in per])),
        "medoid_rmsd_to_native": float(np.mean([r["medoid_rmsd"] for r in per])),
        "rg_pred_err": float(np.mean([abs(r["rg_pred"] - r["rg_nat"]) for r in per])),
        "rg_pool_err": float(np.mean([abs(r["rg_pool_mean"] - r["rg_nat"]) for r in per])),
        "corr_rgerr_rr": _rg_corr(per)}
    write_json("generate_diag.json", out)
    _print_diag(out)
    return out


def _rg_corr(per):
    """Across targets, does a better rg prediction go with a better pool mean?"""
    a = [abs(r["rg_pred"] - r["rg_nat"]) for r in per]
    b = [r["pool"]["sel"] for r in per]
    return float(np.corrcoef(a, b)[0, 1])


def _print_diag(out):
    print(f"\n=== POOL DISTRIBUTION DIAGNOSTICS ({out['n_targets']} targets, "
          f"K={out['K']}) ===")
    print(f"{'pool':22} {'sel':>7} {'best':>7} {'mean':>7} {'sd':>7} {'q10':>7} "
          f"{'q30':>7} {'pct':>6} {'<2A':>6} {'<1.5A':>6} {'ORC*':>7}")
    for tag, name in (("pool", "BLOSUM top-500"), ("univ", "whole universe"),
                      ("oracleK_diagnostic", "ORACLE best-500*")):
        s = out[tag]
        print(f"{name:22} {s['sel']:7.3f} {s['best']:7.3f} {s['mean']:7.3f} "
              f"{s['sd']:7.3f} {s['q10']:7.3f} {s['q30']:7.3f} {s['pct']:6.1f} "
              f"{s['f2']:6.3f} {s['f15']:6.3f} {s['oracle_diagnostic']:7.3f}")
    q = out["quantile_law"]
    print(f"\n--- is the answer a fixed-quantile draw from the pool? ---")
    print(f"  mean selected {q['mean_sel']:.3f} A  ==  the pool's "
          f"{q['matching_quantile']}th percentile on average")
    print(f"  per-target percentile of the selected candidate: "
          f"{q['mean_pct']:.1f} +- {q['sd_pct']:.1f}")
    print(f"  corr(selected, that quantile) {q['corr_sel_q']:+.3f}   "
          f"corr(selected, pool mean) {q['corr_sel_mean']:+.3f}   "
          f"corr(selected, pool best) {q['corr_sel_best']:+.3f}   "
          f"corr(selected, pool sd) {q['corr_sel_sd']:+.3f}")
    g = out["geometry"]
    print(f"\n--- where does the pool mean come from? ---")
    print(f"  mean pairwise CA-RMSD inside the pool (SPREAD)  "
          f"{g['pool_spread_mean_pairwise_rmsd']:.3f} A")
    print(f"  CA-RMSD of the pool MEDOID to the native (CENTRE) "
          f"{g['medoid_rmsd_to_native']:.3f} A")
    print(f"  |rg_predicted - rg_native| {g['rg_pred_err']:.3f} A   "
          f"|rg_poolmean - rg_native| {g['rg_pool_err']:.3f} A   "
          f"corr(rg error, selected) {g['corr_rgerr_rr']:+.3f}")
    print("* ORACLE columns read the native and are upper bounds, not deployable.")


# ============================================================ stage: gen
def _run_arms(tgs, arms, tag=""):
    per = []
    for k, pdbid in enumerate(tgs):
        u = load_univ(pdbid)
        view = deployable_view(u)
        pred = poolsize.load_pred(pdbid)
        n = u["n"]
        Dnat = pair_D(u["nat_ca"][None], n)[0]
        row = {"pdb": pdbid, "n": n}
        for a in arms:
            fn, kw = GENERATORS[a]
            rng = np.random.default_rng(SEED + (abs(hash(pdbid)) % 100000))
            W = fn(view, pred, rng, **kw)
            D = pair_D(W, n)
            rr = audit.kabsch_rmsd_batch(W, u["nat_ca"])      # REPORTING ONLY
            row[a] = _chain(rr, shipped_score(pred, D), Dnat, D)
        o = np.argsort(u["rr"], kind="stable")[:K]
        Do = pair_D(u["W"][o], n)
        row["oracleK_diagnostic"] = _chain(u["rr"][o], shipped_score(pred, Do), Dnat, Do)
        per.append(row)
        print(f"[{k+1:3d}/{len(tgs)}] {pdbid:6} n={n:2d} "
              f"base {row['base']['sel']:5.3f} {tag}", flush=True)
    return per


def _summarise(per, arms, ref="base"):
    base = [r[ref]["sel"] for r in per]
    summary = []
    for a in arms + ["oracleK_diagnostic"]:
        rows = [r[a] for r in per]
        s = {"arm": a, "n_targets": len(rows)}
        for key in ("sel", "best", "mean", "sd", "q10", "q30", "pct", "f2", "f15",
                    "rho_band", "rho_global", "oracle_diagnostic"):
            v, se = _agg(rows, key)
            s[key] = v
            if key == "sel":
                s["sel_se"] = se
        s["vs_base"] = debias.paired([r["sel"] for r in rows], base)
        summary.append(s)
    summary.sort(key=lambda x: (x["arm"] == "oracleK_diagnostic", x["sel"]))
    return summary


def stage_gen(arms=None, out_name="generate_gen.json"):
    tgs = cached_targets()
    arms = arms or list(GENERATORS)
    print(f"{len(arms)} generator arms on {len(tgs)} targets", flush=True)
    per = _run_arms(tgs, arms)
    out = {"n_targets": len(per), "K": K, "instrument": "tuning126",
           "summary": _summarise(per, arms), "per_target": per}
    write_json(out_name, out)
    _print_gen(out)
    return out


def _print_gen(out):
    print(f"\n=== GENERATOR ARMS ({out['n_targets']} targets, {out['instrument']}, "
          f"K={out['K']}) ===")
    print(f"{'arm':14} {'SELECTED':>9} {'+-se':>6} {'d(base)':>8} {'ci95':>17} "
          f"{'W/L':>9} | {'best':>6} {'mean':>6} {'sd':>6} {'q30':>6} "
          f"{'<2A':>6} {'<1.5A':>6} {'band rho':>9} {'pct':>5} | {'ORC*':>6}")
    for s in out["summary"]:
        v = s["vs_base"]
        rb = f"{s['rho_band']:+9.3f}" if s["rho_band"] is not None else f"{'--':>9}"
        print(f"{s['arm']:14} {s['sel']:9.3f} {s['sel_se']:6.3f} "
              f"{v['mean_diff']:+8.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] "
              f"{v['n_better']:3d}/{v['n_worse']:<5d} | {s['best']:6.3f} "
              f"{s['mean']:6.3f} {s['sd']:6.3f} {s['q30']:6.3f} {s['f2']:6.3f} "
              f"{s['f15']:6.3f} {rb} {s['pct']:5.1f} | "
              f"{s['oracle_diagnostic']:6.3f}")
    print("* ORC = the same pool ranked by the NATIVE distance matrix: a DIAGNOSTIC "
          "ceiling for that pool, not deployable.")
    print("oracleK_diagnostic is the best-500 subset of the universe by true RMSD -- "
          "the ceiling on ANY native-free quality filter over this library.")
    _print_lever(out)


def _print_lever(out):
    """THE question: does moving the pool's distribution move selected RMSD?"""
    S = [s for s in out["summary"] if s["arm"] != "oracleK_diagnostic"]
    ref = next(s for s in S if s["arm"] == "base")
    dm = np.array([s["mean"] - ref["mean"] for s in S])
    dq = np.array([s["q30"] - ref["q30"] for s in S])
    dsd = np.array([s["sd"] - ref["sd"] for s in S])
    dbst = np.array([s["best"] - ref["best"] for s in S])
    ds = np.array([s["sel"] - ref["sel"] for s in S])
    print(f"\n--- does the DISTRIBUTION move the ANSWER? ({len(S)} arms) ---")
    for name, x in (("pool mean", dm), ("pool q30", dq), ("pool sd", dsd),
                    ("pool best", dbst)):
        if x.std() < 1e-9:
            continue
        sl = float(np.polyfit(x, ds, 1)[0])
        r = float(np.corrcoef(x, ds)[0, 1])
        print(f"  d(selected) per d({name:9}) : {sl:+.3f}   (r = {r:+.3f})")
    imp = dm < -0.05
    if imp.sum() > 1:
        print(f"  of the {int(imp.sum())} arms that IMPROVE the pool mean (by up to "
              f"{-dm[imp].min():.3f} A): mean d(selected) {ds[imp].mean():+.3f} A, "
              f"best {ds[imp].min():+.3f}, worst {ds[imp].max():+.3f}")
    sig = [s for s in S if s["vs_base"]["ci95"][1] < 0]
    print(f"  arms significantly better than base (95% CI excludes 0): "
          f"{len(sig)}" + (f"  -> {', '.join(s['arm'] for s in sig)}" if sig else ""))


# ============================================================ stage: transfer
#: Spearman targets for the synthetic quality filter. `a` mixes the true RMSD ranking with
#: noise; the ACHIEVED rank correlation is measured, never assumed.
MIXES = (0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0)


def _noisy_key(rr, a, rng):
    """DIAGNOSTIC. A ranking of `rr` degraded to a controlled skill level.

    Models a native-free quality filter of known rank correlation against true RMSD, so
    the payoff of building one can be priced BEFORE building it. Never deployable.
    """
    z = (np.argsort(np.argsort(rr)) + 0.5) / len(rr)
    return a * z + math.sqrt(max(1.0 - a * a, 0.0)) * rng.standard_normal(len(rr)) * z.std()


def stage_transfer():
    """THE price of pool quality: how much selected RMSD does filter SKILL actually buy?

    A synthetic filter of controlled skill selects K candidates out of the full window
    universe. Sweeping its skill from chance to perfect draws the transfer curve from
    'pool distribution' to 'selected RMSD' end to end, and prices every native-free filter
    anyone could build on this library -- including ones that do not exist yet. The real
    keys (BLOSUM, Ramachandran, rg agreement, the shipped score itself) are measured on
    the same universe so their position on that curve can simply be read off.
    """
    tgs = cached_targets()
    print(f"transfer curve on {len(tgs)} targets", flush=True)
    per = []
    for k, pdbid in enumerate(tgs):
        u = load_univ(pdbid)
        view = deployable_view(u)
        pred = poolsize.load_pred(pdbid)
        n = u["n"]
        Dnat = pair_D(u["nat_ca"][None], n)[0]
        Du = pair_D(u["W"], n)
        sc_u = shipped_score(pred, Du)
        rng = np.random.default_rng(SEED + (abs(hash(pdbid)) % 100000))
        row = {"pdb": pdbid, "n": n, "n_windows": int(len(u["W"]))}
        for a in MIXES:
            key = _noisy_key(u["rr"], a, rng)
            idx = np.argsort(key, kind="stable")[:K]
            row[f"mix{a:g}"] = _chain(u["rr"][idx], shipped_score(pred, Du[idx]),
                                      Dnat, Du[idx])
            row[f"mix{a:g}"]["rho_key"] = float(spearmanr(key, u["rr"]).statistic)
        # where the REAL native-free keys sit on that curve: rank correlation of each
        # key against true RMSD over the WHOLE universe.
        rgc = rg_of(Du, n)
        keys = {"blosum": -u["sim"], "rama": -rama_score(view),
                "rg": np.abs(rgc - predicted_rg(pred, n)), "shipped_score": sc_u,
                "typicality": np.abs(Du - Du.mean(0)[None, :]).mean(1)}
        row["real_keys"] = {}
        for name, key in keys.items():
            idx = np.argsort(key, kind="stable")[:K]
            st = _chain(u["rr"][idx], shipped_score(pred, Du[idx]), Dnat, Du[idx])
            st["rho_key"] = float(spearmanr(key, u["rr"]).statistic)
            row["real_keys"][name] = st
        per.append(row)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    out = {"n_targets": len(per), "K": K, "instrument": "tuning126", "per_target": per}
    base = [r["real_keys"]["blosum"]["sel"] for r in per]
    out["curve"] = []
    for a in MIXES:
        rows = [r[f"mix{a:g}"] for r in per]
        s = {"mix": a}
        for key in ("rho_key", "sel", "best", "mean", "sd", "q30", "pct", "f2", "f15",
                    "rho_band", "oracle_diagnostic"):
            s[key] = _agg(rows, key)[0]
        s["sel_se"] = _agg(rows, "sel")[1]
        s["vs_blosum"] = debias.paired([r["sel"] for r in rows], base)
        out["curve"].append(s)
    out["real_keys"] = []
    for name in ("blosum", "rama", "rg", "shipped_score", "typicality"):
        rows = [r["real_keys"][name] for r in per]
        s = {"key": name}
        for key in ("rho_key", "sel", "best", "mean", "sd", "q30", "pct", "f2", "f15",
                    "rho_band", "oracle_diagnostic"):
            s[key] = _agg(rows, key)[0]
        s["sel_se"] = _agg(rows, "sel")[1]
        s["vs_blosum"] = debias.paired([r["sel"] for r in rows], base)
        out["real_keys"].append(s)
    write_json("generate_transfer.json", out)
    _print_transfer(out)
    return out


def _print_transfer(out):
    print(f"\n=== TRANSFER CURVE: what does FILTER SKILL buy? "
          f"({out['n_targets']} targets, K={out['K']} out of the full universe) ===")
    print("A synthetic filter of controlled rank correlation against true CA-RMSD picks "
          "K windows.\nDIAGNOSTIC: it reads natives. It prices native-free filters that "
          "do not exist yet.")
    print(f"{'rho(key,rmsd)':>13} {'SELECTED':>9} {'+-se':>6} {'d(blosum)':>10} "
          f"{'ci95':>17} {'W/L':>9} | {'best':>6} {'mean':>6} {'sd':>6} {'q30':>6} "
          f"{'<2A':>6} {'<1.5A':>6} {'band rho':>9} {'pct':>5} | {'ORC*':>6}")
    for s in out["curve"]:
        v = s["vs_blosum"]
        rb = f"{s['rho_band']:+9.3f}" if s["rho_band"] is not None else f"{'--':>9}"
        print(f"{s['rho_key']:+13.3f} {s['sel']:9.3f} {s['sel_se']:6.3f} "
              f"{v['mean_diff']:+10.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] "
              f"{v['n_better']:3d}/{v['n_worse']:<5d} | {s['best']:6.3f} {s['mean']:6.3f} "
              f"{s['sd']:6.3f} {s['q30']:6.3f} {s['f2']:6.3f} {s['f15']:6.3f} {rb} "
              f"{s['pct']:5.1f} | {s['oracle_diagnostic']:6.3f}")
    print(f"\n--- the REAL native-free keys on the same universe ---")
    print(f"{'key':14} {'rho(key,rmsd)':>13} {'SELECTED':>9} {'+-se':>6} "
          f"{'d(blosum)':>10} {'ci95':>17} {'W/L':>9} | {'best':>6} {'mean':>6} "
          f"{'sd':>6} {'q30':>6} {'<2A':>6} {'band rho':>9}")
    for s in out["real_keys"]:
        v = s["vs_blosum"]
        rb = f"{s['rho_band']:+9.3f}" if s["rho_band"] is not None else f"{'--':>9}"
        print(f"{s['key']:14} {s['rho_key']:+13.3f} {s['sel']:9.3f} {s['sel_se']:6.3f} "
              f"{v['mean_diff']:+10.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}] "
              f"{v['n_better']:3d}/{v['n_worse']:<5d} | {s['best']:6.3f} {s['mean']:6.3f} "
              f"{s['sd']:6.3f} {s['q30']:6.3f} {s['f2']:6.3f} {rb}")
    c = out["curve"]
    dm = np.array([s["mean"] - c[0]["mean"] for s in c])
    ds = np.array([s["sel"] - c[0]["sel"] for s in c])
    if dm.std() > 1e-9:
        print(f"\n  d(selected) per d(pool mean) along the curve: "
              f"{float(np.polyfit(dm, ds, 1)[0]):+.3f}  "
              f"(r = {float(np.corrcoef(dm, ds)[0, 1]):+.3f})")
    print("* ORC = the same pool ranked by the native distance matrix.")


# ============================================================ stage: dev
def stage_dev():
    """ONE pass on the 24 dev targets for the single best tuning arm. No iteration."""
    raise SystemExit("dev pass is gated on a tuning winner; see stage_gen output")


def stage_soft():
    """The soft-filter sweep, paired against `base` on the same 126 targets."""
    arms = ["base"] + [a for a in GENERATORS if a.startswith("soft_")]
    return stage_gen(arms, "generate_soft.json")


# ============================================================ stage: law
def _all_arms():
    """Every arm measured this sprint, from all three sweeps, on identical targets."""
    out = []
    for fname in ("generate_gen.json", "generate_soft.json"):
        try:
            d = read_json(fname)
        except FileNotFoundError:
            continue
        for s in d["summary"]:
            if any(o["name"] == s["arm"] for o in out):    # base / oracleK are in both
                continue
            out.append({"name": s["arm"], "kind": "oracle" if "oracle" in s["arm"]
                        else ("soft" if s["arm"].startswith("soft_") else "generator"),
                        **{k: s[k] for k in ("sel", "best", "mean", "sd", "q30",
                                             "rho_band", "pct", "f2", "f15")}})
    try:
        t = read_json("generate_transfer.json")
        for s in t["curve"]:
            out.append({"name": f"synth_rho{s['rho_key']:+.2f}", "kind": "synthetic",
                        **{k: s[k] for k in ("sel", "best", "mean", "sd", "q30",
                                             "rho_band", "pct", "f2", "f15")}})
        for s in t["real_keys"]:
            out.append({"name": f"key_{s['key']}", "kind": "real_key",
                        **{k: s[k] for k in ("sel", "best", "mean", "sd", "q30",
                                             "rho_band", "pct", "f2", "f15")}})
    except FileNotFoundError:
        pass
    return out


def _fit(rows, cols):
    y = np.array([r["sel"] for r in rows], float)
    A = np.column_stack([[r[c] for r in rows] for c in cols] + [np.ones(len(rows))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    ss = ((y - y.mean()) ** 2).sum()
    return {"cols": list(cols), "coef": [float(c) for c in coef],
            "r2": float(1.0 - (resid ** 2).sum() / ss) if ss > 0 else None,
            "rmse": float(np.sqrt((resid ** 2).mean())), "n": len(rows)}


def stage_law():
    """Is there ONE law from a pool's distribution to the answer it yields?

    Every arm this sprint measured -- 26 generators, 15 soft filters, 8 synthetic filter
    skills, 5 real native-free keys, the oracle subset -- was run on the SAME 126 targets
    with the SAME selector. If mean selected RMSD is a simple function of a pool's own
    summary statistics across all of them, then that function prices any generator anyone
    might build next WITHOUT building it, and the sprint's question has a numerical answer
    rather than a list of failures.
    """
    rows = _all_arms()
    if not rows:
        raise SystemExit("run `gen`, `soft` and `transfer` first")
    fits = {}
    for cols in (("mean",), ("best",), ("q30",), ("sd",), ("best", "q30"),
                 ("best", "mean"), ("best", "mean", "sd"), ("best", "q30", "rho_band")):
        if any(any(r[c] is None for r in rows) for c in cols):
            continue
        fits["+".join(cols)] = _fit(rows, cols)
    best = min(fits.values(), key=lambda f: f["rmse"])
    dep = [r for r in rows if r["kind"] in ("generator", "soft", "real_key")]
    out = {"n_arms": len(rows), "n_deployable": len(dep), "fits": fits,
           "best_fit": best, "arms": rows,
           "deployable_range": {
               "sel_min": float(min(r["sel"] for r in dep)),
               "sel_max": float(max(r["sel"] for r in dep)),
               "mean_min": float(min(r["mean"] for r in dep)),
               "best_min": float(min(r["best"] for r in dep))},
           "oracle": next((r for r in rows if r["kind"] == "oracle"), None)}
    write_json("generate_law.json", out)
    _print_law(out)
    return out


def _print_law(out):
    print(f"\n=== ONE LAW FROM POOL DISTRIBUTION TO SELECTED RMSD "
          f"({out['n_arms']} arms, all on the same 126 targets) ===")
    print(f"{'predictors':28} {'R2':>7} {'RMSE(A)':>8} {'coefficients':>34}")
    for name, f in sorted(out["fits"].items(), key=lambda kv: kv[1]["rmse"]):
        cs = "  ".join(f"{c}{v:+.3f}" for c, v in zip(f["cols"] + ["const"], f["coef"]))
        print(f"{name:28} {f['r2']:7.3f} {f['rmse']:8.3f} {cs:>34}")
    print(f"\n{'arm':22} {'kind':10} {'SELECTED':>9} {'best':>7} {'mean':>7} "
          f"{'q30':>7} {'sd':>7} {'<2A':>6}")
    for r in sorted(out["arms"], key=lambda x: x["sel"]):
        print(f"{r['name']:22} {r['kind']:10} {r['sel']:9.3f} {r['best']:7.3f} "
              f"{r['mean']:7.3f} {r['q30']:7.3f} {r['sd']:7.3f} {r['f2']:6.3f}")
    d = out["deployable_range"]
    print(f"\n  {out['n_deployable']} DEPLOYABLE arms span selected "
          f"{d['sel_min']:.3f}-{d['sel_max']:.3f} A; the best pool mean any of them "
          f"reaches is {d['mean_min']:.3f} A")
    o = out["oracle"]
    if o:
        print(f"  the ORACLE best-{K} subset reaches pool mean {o['mean']:.3f} A and "
              f"selected {o['sel']:.3f} A -- the ceiling of this library")


def stage_report():
    for name, printer in (("generate_diag.json", _print_diag),
                          ("generate_transfer.json", _print_transfer),
                          ("generate_gen.json", _print_gen),
                          ("generate_soft.json", _print_gen),
                          ("generate_law.json", _print_law)):
        try:
            printer(read_json(name))
        except FileNotFoundError:
            print(f"({name} not run)")


STAGES = {"univ": stage_univ, "diag": stage_diag, "gen": stage_gen,
          "transfer": stage_transfer, "soft": stage_soft, "law": stage_law,
          "dev": stage_dev,
          "report": stage_report}


def main(argv):
    if len(argv) < 2 or argv[1] not in STAGES:
        print(__doc__)
        print("stages: " + ", ".join(STAGES))
        return 1
    STAGES[argv[1]]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
