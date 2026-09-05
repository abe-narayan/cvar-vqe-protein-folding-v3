"""When every candidate is already good, what -- if anything -- discriminates?

THE PHENOMENON
--------------
Sprint 8's generation study fitted a law across 55 pool-construction arms:

    selected = 0.257 * pool_best + 0.311 * pool_mean + 1.673      (R^2 0.79)

and the 1.673 A constant is the selector's incompetence.  No pool or representation work
touches it.  Worse, on the ORACLE-tight pool -- the true best 500 windows of the whole
universe, pool mean 2.350 A, 41% of candidates already under 2 A -- the shipped distogram
score still selects 2.406 A and its percentile DEGRADES from 24.9 to 59.9, worse than a
coin flip.  In the near-native regime the shipped objective is not merely uninformative;
it looks anti-informative.

This module asks whether ANY native-free signal has skill inside the near-native band.

THE INSTRUMENT
--------------
126-target tuning set (`s7.debias.tuning_targets`), cluster-disjoint from the 24-target dev
set and the 60-target benchmark, on the cached window universes `s8/generate_univ/*.npz`.
Two pools per target, both size K=500 over the same universe:

    nrm   the shipped BLOSUM top-500                      (deployable construction)
    tgt   the 500 lowest-RMSD windows in the universe      (ORACLE, DIAGNOSTIC)

and inside each, the near-native BAND: candidates within `BAND` = 1.5 A of pool best.

DIAGNOSTIC vs DEPLOYABLE -- the split is structural
---------------------------------------------------
The tight pool and the band are ORACLE constructions: they read `rr` (each window's
CA-RMSD to the native) to decide WHICH candidates are measured on.  They are diagnostics
and feed nothing.  Every SIGNAL is computed by `signal_block()` from a context built by
`deployable_ctx()`, which contains the candidate coordinates, their torsions, the target
sequence, the fold index and the fold's own out-of-fold statistics -- and nothing derived
from the target's native structure.  `stage_leak` NaN-poisons `rr` and `nat_ca` in every
cached universe and asserts every signal is bit-identical.

WHAT IS MEASURED
----------------
For every signal, on both pools:

  * in-band Spearman rho against true CA-RMSD (the diagnostic), and the same rho with
    radius of gyration partialled out -- compactness has produced two false leads in this
    project already;
  * `sel_band`, the CA-RMSD of the candidate the signal picks WITHIN the band, and
    `sel_pool`, what it picks over the whole pool (the deployable number);
  * mean / median / SD / worst / best, fraction below 2.0 and 1.5 A, paired differences
    against the shipped distogram on identical candidates with 95% CIs and W/L counts.

Sprint 7 established that in-band rho is NOT sufficient -- the arm with the best in-band
rho of 28 had a worse selected RMSD -- so every rho is converted to a selected RMSD before
anything is believed.

GEOMETRY NOTE
-------------
A candidate is a real length-n CA window carrying its parent's real (phi, psi).  The
Legacy energy needs N/C/O/CB, which only exist on an ideal-geometry rebuild from those
torsions; that rebuild sits 0.27 A CA-RMSD from the real window on average.  Legacy terms
are therefore evaluated on the rebuild while the outcome `rr` is the real window's, exactly
as sprint 7's Amber comparison did.  `stage_build` records the rebuild displacement so the
size of that caveat is on the record rather than assumed away.

    python -m s8.inband build    # per-target signal cache (both pools)
    python -m s8.inband skill    # the in-band skill table for every signal
    python -m s8.inband anti     # is the anti-correlation real, and is its sign stable?
    python -m s8.inband ens      # ensemble disagreement across independently trained models
    python -m s8.inband comb     # leave-fold-out learned combiners
    python -m s8.inband two      # the DEPLOYABLE two-stage selector the above implies
    python -m s8.inband dev      # THE single dev pass, one pre-registered arm
    python -m s8.inband leak     # NaN-poison assertion for every deployable signal
    python -m s8.inband report   # print everything already computed

RESULT
------
1. THE PREMISE WAS HALF RIGHT.  In-band rho does not go NEGATIVE in the near-native regime,
   it collapses to ZERO and stays there: the shipped score runs +0.568 (whole pool) ->
   +0.126 (1.5 A band) -> -0.020 (0.5 A band) by band width, and +0.126 (pool mean 4.45) ->
   +0.020 (2.85) -> -0.008 (2.35) -> -0.028 (1.94) by pool tightness, with ~48% of targets
   on either side.  It crosses zero at a pool mean of ~2.8 A.

2. AND IT IS STILL WORSE THAN A COIN FLIP.  Against the pool's own mean -- what a random
   draw returns -- the score selects +0.036 to +0.067 A WORSE on the oracle-tight pools at
   m = 50/100/250/500, losing on 84-94 of 126 targets, four CIs excluding zero.  Zero rank
   information plus a biased error still lands on the wrong side of a tight distribution.

3. ONE THING DISCRIMINATES, AND IT HAS NO MODEL IN IT.  Of 38 native-free signals -- the
   distogram, its confidence, its sd-reweightings, its bagged variant, three independently
   trained distograms' disagreement, all 11 Legacy terms and 4 subsets, H-bonding, burial,
   packing, contact order, strain, Ramachandran likelihood, SS self-consistency -- exactly
   one has positive in-band skill in the tight regime: TYPICALITY, agreement with the
   pool's own mean distance profile.  rho +0.125 (partialling rg STRENGTHENS it to +0.177),
   sign stable in all five folds, and better than a random draw by 0.09-0.22 A at every
   pool tightness with CIs excluding zero.

4. LEGACY, PER TERM, IS NOT A DISCRIMINATOR -- and two thirds of its apparent tight-pool
   anti-correlation is radius of gyration (`lg_all` -0.078 -> -0.025 partialled).  That is
   the third compactness false lead in this project.  As a LEARNED 15-column combination it
   does reach 2.228 A on the tight pool against 2.406, and 4.103 A on the normal pool
   against 3.454: Legacy helps only once the pool is already near-native.

5. THE DEPLOYABLE CONSEQUENCE.  The score is a good coarse filter and a bad fine ranker;
   consensus is the reverse.  Compose them (`stage_two`): filter the shipped pool to the
   score's top 75, return the consensus medoid.  3.282 A against 3.454, -0.172
   [-0.316, -0.027], 74W/47L, chosen leave-fold-out over 20 arms with all five folds
   agreeing.  Dev replicates the magnitude (-0.159) and cannot resolve it (SE 0.354).
   The score's ordering INSIDE its own top 75 is worth 0.10 A; consensus over the identical
   set is worth 0.27 A.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np                                                          # noqa: E402
from scipy.stats import spearmanr, rankdata                                 # noqa: E402

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:                                                           # pragma: no cover
    torch = None

import energy_terms as et                                                   # noqa: E402
import legacy_field as lfd                                                  # noqa: E402
import peptide_db as db                                                     # noqa: E402
import protein_geometry as geo                                              # noqa: E402
from s7 import audit, debias, poolsize                                      # noqa: E402
from s8 import generate as gen                                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "inband_cache")
SSTAB = os.path.join(HERE, "inband_ss.npz")

K = 500                #: pool size every sprint-7/8 number reports on
BAND = 1.5             #: in-band = within pool_best + BAND
MIN_BAND = 5           #: fewer in-band candidates than this and rho is undefined
GRID = gen.GRID
CA_CA = 3.80
SEED = 20260905
POOLS = ("nrm", "tgt")

#: Band widths swept by `stage_anti`, in A above pool best.
BANDS = (0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 1e9)
#: Pool TIGHTNESS sweep: the size of the ORACLE-ranked prefix taken as the pool.  Pool
#: size IS the tightness knob here -- the top-25 of the universe is a far tighter pool than
#: its top-8000 -- and it spans the whole range from "every candidate is near-native" to
#: "the pool is the universe".  0 means the deployable BLOSUM top-500 reference.
TIGHTNESS = (25, 50, 100, 250, 500, 1000, 2000, 8000, 0)


# ============================================================ out-of-fold SS table
#: Torsion-space three-state assignment.  Applied identically to the training corpus and
#: to a candidate, so the "self-consistency" signal compares like with like and needs no
#: backbone rebuild on either side.
def ss_from_torsions(PHI, PSI):
    """(B, n) radians -> (B, n) int in {0: helix, 1: sheet, 2: coil}."""
    phi = np.degrees(np.asarray(PHI, float))
    psi = np.degrees(np.asarray(PSI, float))
    out = np.full(phi.shape, 2, np.int8)
    sheet = (phi >= -180) & (phi <= -45) & ((psi >= 90) | (psi <= -150))
    helix = (phi >= -160) & (phi <= -20) & (psi >= -120) & (psi <= 50)
    out[sheet] = 1
    out[helix] = 0
    return out


def _ss_counts(entries, tab):
    for e in entries:
        a = audit.encode(e.seq)
        s = ss_from_torsions(np.asarray(e.phi)[None, :], np.asarray(e.psi)[None, :])[0]
        np.add.at(tab, (a, s), 1.0)


def build_ss_table(verbose=True):
    """(5, 20, 3) out-of-fold P(SS | amino acid).  Never sees a target's structure."""
    if os.path.exists(SSTAB):
        return
    import distogram as dgm
    folds = db.folds(5)
    tab = np.zeros((5, 20, 3))
    memo = {}
    _orig = dgm._fold_fragments

    def frags(fold, n_folds, threshold=db.IDENTITY_THRESHOLD):
        key = (fold, n_folds, threshold)
        if key not in memo:
            memo[key] = _orig(fold, n_folds, threshold)
        return memo[key]

    dgm._fold_fragments = frags
    try:
        for fold in range(5):
            _ss_counts([q for q in db.load() if folds[q.seq] != fold], tab[fold])
            _ss_counts(list(frags(fold, 5)), tab[fold])
            if verbose:
                print(f"  ss fold {fold}: {tab[fold].sum():.0f} residues", flush=True)
    finally:
        dgm._fold_fragments = _orig
    np.savez_compressed(SSTAB, cnt=tab)


_SS = [None]


def ss_logp(fold):
    if _SS[0] is None:
        build_ss_table(verbose=False)
        c = np.load(SSTAB)["cnt"] + 0.5
        _SS[0] = np.log(c / c.sum(2, keepdims=True))
    return _SS[0][fold]


# ============================================================ DEPLOYABLE context
# Nothing between here and the DIAGNOSTIC banner may touch a native quantity.
def deployable_ctx(u, pred, idx):
    """Everything a signal is allowed to see for one target's candidate subset.

    `u` is a cached universe; only its native-free arrays are read.  `idx` selects the
    candidates -- WHICH candidates were chosen may be an oracle decision (that is the
    diagnostic construction), but nothing downstream of this function can tell.
    """
    n = int(u["n"])
    seq = str(u["seq"])
    fold = int(u["fold"])
    W = np.asarray(u["W"], float)[idx]
    PHI = np.asarray(u["PHI"], float)[idx]
    PSI = np.asarray(u["PSI"], float)[idx]
    D = gen.pair_D(W, n)
    BB = geo.build_backbone_batch(PHI, PSI)
    return {"n": n, "seq": seq, "fold": fold, "W": W, "PHI": PHI, "PSI": PSI,
            "D": D, "BB": BB, "pred": pred, "B": len(idx)}


# --------------------------------------------------------------- signal helpers
def _reweight(pred, gamma):
    """The shipped Bayes-risk table re-weighted to a different sd exponent.

    `Distogram` builds `_risk = base_risk * w` with `w = shell / (sd + 0.5) ** gamma`
    normalised to mean 1, and the shipped shell weights are all 1 (no `score_weights.json`
    in this tree, checked).  Dividing the stored risk by the stored `w` recovers the
    unweighted table exactly, so a gamma arm is an exact re-weighting rather than a
    re-derivation.
    """
    w = np.asarray(pred["w"], float)
    base = np.asarray(pred["risk"], float) / np.maximum(w[:, None], 1e-12)
    sd = np.asarray(pred["sd"], float)
    w2 = 1.0 / np.maximum((sd + 0.5) ** gamma, 1e-6)
    w2 = w2 / max(w2.mean(), 1e-12)
    return base * w2[:, None]


def _risk_score(risk, D):
    g = np.clip(((np.asarray(D, float) - GRID[0]) / 0.05).astype(np.int32),
                0, len(GRID) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g].mean(1)


def _risk_per_pair(risk, D):
    g = np.clip(((np.asarray(D, float) - GRID[0]) / 0.05).astype(np.int32),
                0, len(GRID) - 1)
    return risk[np.arange(risk.shape[0])[None, :], g]


def _contact_order(D, i, j, n, cut=8.0):
    c = (D < cut).astype(float)
    sep = (j - i).astype(float)
    tot = c.sum(1)
    return np.where(tot > 0, (c * sep[None, :]).sum(1) / np.maximum(tot, 1) / n, 0.0)


# --------------------------------------------------------------- the signal block
#: Signals are scores: LOWER is the prediction.  A signal whose natural direction is
#: unknown is still stored as-is and its sign is a result, not an assumption.
def signal_block(ctx):
    """(B, n_sig) native-free signal matrix, and the signal names."""
    n, seq, fold = ctx["n"], ctx["seq"], ctx["fold"]
    D, W, BB = ctx["D"], ctx["W"], ctx["BB"]
    PHI, PSI = ctx["PHI"], ctx["PSI"]
    pred = ctx["pred"]
    i, j = audit.pair_index(n)
    sep = (j - i).astype(float)
    exp = np.asarray(pred["expected"], float)
    sd = np.asarray(pred["sd"], float)
    w = np.asarray(pred["w"], float)
    risk = np.asarray(pred["risk"], float)

    out, names = [], []

    def add(name, v):
        v = np.asarray(v, float)
        assert v.shape == (len(D),), (name, v.shape)
        out.append(v)
        names.append(name)

    # ---------------------------------------------------------- distogram family
    add("score", _risk_score(risk, D))                       # the shipped baseline
    add("score_g0", _risk_score(_reweight(pred, 0.0), D))
    add("score_g2", _risk_score(_reweight(pred, 2.0), D))
    add("l1exp", np.abs(D - exp[None, :]).mean(1))
    add("l1exp_w", (np.abs(D - exp[None, :]) * w[None, :]).mean(1))
    add("zdev", (np.abs(D - exp[None, :]) / np.maximum(sd, 1e-6)[None, :]).mean(1))
    # confidence slices: does the score work only where the model is sure?
    q = np.quantile(sd, [0.25, 0.75])
    lo, hi = sd <= q[0], sd >= q[1]
    pp = _risk_per_pair(risk, D)
    add("score_lowsd", pp[:, lo].mean(1) if lo.sum() else np.zeros(len(D)))
    add("score_highsd", pp[:, hi].mean(1) if hi.sum() else np.zeros(len(D)))
    # sd-only: no candidate dependence would be constant, so use the sd-weighted spread
    add("sdmatch", (np.abs(np.abs(D - exp[None, :]) - sd[None, :])).mean(1))
    # bagged rank over random pair subsets: the shipped robustness variant, never
    # measured in-band
    rng = np.random.default_rng(SEED)
    k = max(4, int(0.6 * pp.shape[1]))
    acc = np.zeros(len(D))
    for _ in range(12):
        cols = rng.choice(pp.shape[1], size=k, replace=False)
        acc += rankdata(pp[:, cols].mean(1))
    add("bagged", acc / 12.0)
    # typicality: agreement with the pool's OWN mean profile (no sequence at all)
    Dbar = D.mean(0)
    add("typic", np.abs(D - Dbar[None, :]).mean(1))
    # rg agreement with the rg implied by the predicted matrix
    rg = gen.rg_of(D, n)
    add("rgpred", np.abs(rg - gen.predicted_rg(pred, n)))
    add("rg", rg)                                            # raw compactness

    # ---------------------------------------------------------- Legacy, 11 terms
    bl = lfd.BatchLegacy(seq, _RepShim(PHI, PSI))
    T = bl.terms_from_coords(BB, phi=PHI, psi=PSI)
    for c, t in enumerate(lfd.TERMS):
        add(f"lg_{t}", T[:, c])
    wv = np.array([lfd.FITTED_WEIGHTS.get(t, 0.0) for t in lfd.TERMS], float)
    add("lg_all", T @ wv)
    for sub, keep in lfd.SUBSETS.items():
        if keep is None:
            continue
        ws = np.array([wv[c] if t in keep else 0.0
                       for c, t in enumerate(lfd.TERMS)], float)
        add(f"lg_{sub}", T @ ws)

    # ---------------------------------------------------------- geometry / physics
    # H-bond COUNTS, from the same batched greedy DSSP match the Legacy energies use.
    # The energies are already two Legacy terms; the counts are not, and a count is the
    # quantity a "does this look like a folded peptide" argument actually appeals to.
    bond = bl._hbonds(BB)[2]
    sepm = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    add("nhb", -bond.sum((1, 2)).astype(float))              # more H-bonds = better
    add("nhb_lr", -(bond & (sepm[None] >= 5)).sum((1, 2)).astype(float))
    add("ncontact", -(D < 8.0).sum(1).astype(float))
    add("cordr", _contact_order(D, i, j, n))
    # backbone strain of the REAL window: consecutive CA-CA away from 3.80 A
    step = np.linalg.norm(W[:, 1:] - W[:, :-1], axis=2)
    add("strain", np.sqrt(((step - CA_CA) ** 2).mean(1)))
    # Ramachandran plausibility under the fold's own out-of-fold table
    L = gen.rama_logp(fold)
    a = audit.encode(seq)
    bp = np.clip(((PHI + math.pi) / (2 * math.pi) * gen.RB).astype(int), 0, gen.RB - 1)
    bs = np.clip(((PSI + math.pi) / (2 * math.pi) * gen.RB).astype(int), 0, gen.RB - 1)
    lp = L[a[None, :], bp, bs]
    add("rama", -lp.mean(1))
    add("rama_min", -lp.min(1))
    # secondary-structure self-consistency: does the candidate's SS match the SS this
    # sequence prefers, under an out-of-fold P(SS | aa) table?
    SL = ss_logp(fold)
    ss = ss_from_torsions(PHI, PSI)
    add("ssfit", -SL[a[None, :], ss].mean(1))
    add("sshelix", (ss == 0).mean(1))
    add("sssheet", (ss == 1).mean(1))

    return np.stack(out, 1), names


class _RepShim:
    """The minimal `representation` surface `BatchLegacy` reads: `_phi`, `_psi`, `n_states`.

    `BatchLegacy.__init__` builds a per-residue x per-state Ramachandran penalty table from
    it, which the coords path does not use (we pass `phi`/`psi` explicitly), so one state
    of dummy torsions is enough and costs nothing.
    """

    def __init__(self, PHI, PSI):
        self._phi = np.asarray(PHI, float)[:1].T.copy()
        self._psi = np.asarray(PSI, float)[:1].T.copy()
        self.n_states = 1
        self.n_residues = self._phi.shape[0]


# ============================================================ DIAGNOSTIC path
# Everything below reads a native.  None of it is deployable.
def pool_index(u, kind, m=K):
    """Candidate indices for one pool. `tgt` is an ORACLE construction (reads `rr`).

    `m` is the pool size; for `tgt` it is also the TIGHTNESS knob, since the top-25 of the
    oracle order is a far tighter pool than its top-2000.
    """
    if kind == "nrm":
        return np.asarray(u["order"], int)[:m]
    return np.argsort(np.asarray(u["rr"], float), kind="stable")[:m]


def _cpath(pdbid):
    return os.path.join(CACHE, f"{pdbid}.npz")


def stage_build(verbose=True):
    """Cache the signal matrix and `rr` for both pools of every target."""
    os.makedirs(CACHE, exist_ok=True)
    tgs = gen.cached_targets()
    todo = [p for p in tgs if not os.path.exists(_cpath(p))]
    if verbose:
        print(f"{len(todo)}/{len(tgs)} targets to build", flush=True)
    names = None
    for k, pdbid in enumerate(todo):
        t0 = time.time()
        u = gen.load_univ(pdbid)
        pred = poolsize.load_pred(pdbid)
        blob = {"pdb": pdbid, "n": u["n"], "fold": u["fold"], "seq": u["seq"]}
        for kind in POOLS:
            idx = pool_index(u, kind)
            ctx = deployable_ctx(u, pred, idx)
            S, names = signal_block(ctx)
            blob[f"S_{kind}"] = S.astype(np.float32)
            blob[f"rr_{kind}"] = np.asarray(u["rr"], float)[idx].astype(np.float32)
            blob[f"rg_{kind}"] = gen.rg_of(ctx["D"], ctx["n"]).astype(np.float32)
            # how far the ideal-geometry rebuild the Legacy terms are scored on sits
            # from the real window whose RMSD is the outcome
            R = ctx["BB"]["CA"]
            blob[f"reb_{kind}"] = np.array(
                [geo.rmsd(geo.kabsch_superpose(R[b], ctx["W"][b]), ctx["W"][b])
                 for b in range(len(R))], np.float32)
        blob["names"] = np.array(names)
        tmp = _cpath(pdbid) + ".tmp.npz"
        np.savez_compressed(tmp, **blob)
        os.replace(tmp, _cpath(pdbid))
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {pdbid} n={u['n']} "
                  f"({time.time()-t0:.1f}s)", flush=True)
    if verbose:
        print(f"built {len(todo)}; cache holds "
              f"{len([f for f in os.listdir(CACHE) if f.endswith('.npz')])}", flush=True)


def load_cache(pdbid):
    z = np.load(_cpath(pdbid), allow_pickle=True)
    out = {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
           "seq": str(z["seq"]), "names": [str(x) for x in z["names"]]}
    for kind in POOLS:
        out[f"S_{kind}"] = z[f"S_{kind}"].astype(float)
        out[f"rr_{kind}"] = z[f"rr_{kind}"].astype(float)
        out[f"rg_{kind}"] = z[f"rg_{kind}"].astype(float)
        out[f"reb_{kind}"] = z[f"reb_{kind}"].astype(float)
    return out


def cached():
    if not os.path.isdir(CACHE):
        return []
    return sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))


# --------------------------------------------------------------- statistics
def sel_of(v, rr, invert=False):
    """Expected CA-RMSD of the candidate a signal picks, under RANDOM tie-breaking.

    THIS FUNCTION IS LOAD-BEARING AND ITS FIRST VERSION WAS WRONG.  `np.argmin` returns the
    FIRST index of a tie, and the candidate order is informative in both pools -- BLOSUM
    rank in `nrm`, and true CA-RMSD in the ORACLE-tight pool `tgt`.  Several signals are
    massively tied (`lg_steric` is exactly 0.0 for most candidates, H-bond counts are small
    integers), so plain `argmin` handed them the pool's first element and `lg_steric`
    "selected" 1.386 A on the tight pool against the shipped score's 2.406 -- a pure
    read-out of the oracle sort order, not a signal.

    Averaging the true RMSD over the whole tied argmin set is the expectation of what a
    real selector delivers when it breaks ties at random, is deterministic, and has lower
    variance than drawing one.  A signal that is constant therefore scores the pool mean,
    which is exactly what "carries no information" should cost.
    """
    v = np.asarray(v, float)
    m = v.max() if invert else v.min()
    return float(np.asarray(rr, float)[v == m].mean())


def _rho(x, y):
    if len(x) < MIN_BAND or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(spearmanr(x, y).statistic)


def partial_rho(x, y, z):
    """Spearman rho of x and y with z partialled out, on ranks."""
    if len(x) < MIN_BAND:
        return None
    R = np.column_stack([rankdata(x), rankdata(y), rankdata(z)]).astype(float)
    R -= R.mean(0)
    s = R.std(0)
    if (s < 1e-9).any():
        return None
    C = np.corrcoef((R / s).T)
    den = math.sqrt(max((1 - C[0, 2] ** 2) * (1 - C[1, 2] ** 2), 1e-12))
    return float((C[0, 1] - C[0, 2] * C[1, 2]) / den)


def per_target_stats(S, names, rr, rg, band_w=BAND):
    """Every signal's in-band and whole-pool behaviour on one target's pool."""
    best = float(rr.min())
    m = rr <= best + band_w
    rows = {}
    for c, nm in enumerate(names):
        v = S[:, c]
        r = {"sel_pool": sel_of(v, rr),
             "rho_global": _rho(v, rr)}
        if m.sum() >= MIN_BAND:
            r["sel_band"] = sel_of(v[m], rr[m])
            r["rho_band"] = _rho(v[m], rr[m])
            r["rho_band_prg"] = partial_rho(v[m], rr[m], rg[m])
        rows[nm] = r
    return rows, {"best": best, "mean": float(rr.mean()),
                  "sd": float(rr.std(ddof=1)), "n_band": int(m.sum()),
                  "f2": float((rr < 2.0).mean()), "f15": float((rr < 1.5).mean())}


def _mean_se(v):
    v = np.asarray([x for x in v if x is not None], float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return None, None, int(len(v))
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), int(len(v))


def write_json(name, obj):
    p = os.path.join(HERE, name)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    os.replace(tmp, p)


def read_json(name):
    with open(os.path.join(HERE, name)) as fh:
        return json.load(fh)


# ============================================================ stage: skill
def _paired(a, b):
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[np.isfinite(d)]
    n = len(d)
    if n < 2:
        return None
    se = float(d.std(ddof=1) / math.sqrt(n))
    return {"mean_diff": float(d.mean()), "se": se,
            "ci95": [float(d.mean() - 1.96 * se), float(d.mean() + 1.96 * se)],
            "n_better": int((d < -1e-9).sum()), "n_worse": int((d > 1e-9).sum()),
            "n_tied": int((np.abs(d) <= 1e-9).sum()), "n": n}


def stage_skill():
    """The in-band skill table for every signal, on both pools.

    Each signal is scored in BOTH directions -- a reliably anti-correlated signal is a
    usable selector if and only if its sign is stable, and finding 5 of the brief is
    exactly that question -- so `sel_band_inv` picks the signal's ARGMAX instead.
    """
    tgs = cached()
    print(f"in-band skill on {len(tgs)} targets", flush=True)
    per, names = [], None
    for pdbid in tgs:
        c = load_cache(pdbid)
        names = c["names"]
        row = {"pdb": pdbid, "n": c["n"], "fold": c["fold"]}
        for kind in POOLS:
            S, rr, rg = c[f"S_{kind}"], c[f"rr_{kind}"], c[f"rg_{kind}"]
            sig, pool = per_target_stats(S, names, rr, rg)
            best = float(rr.min())
            m = rr <= best + BAND
            for nm in names:
                if m.sum() >= MIN_BAND:
                    v = S[:, names.index(nm)]
                    sig[nm]["sel_band_inv"] = sel_of(v[m], rr[m], invert=True)
                    sig[nm]["sel_pool_inv"] = sel_of(v, rr, invert=True)
            row[kind] = {"pool": pool, "sig": sig}
            row[f"reb_{kind}"] = float(c[f"reb_{kind}"].mean())
        per.append(row)
    out = {"n_targets": len(per), "K": K, "BAND": BAND, "names": names,
           "per_target": per, "agg": {}}
    for kind in POOLS:
        out["agg"][kind] = _aggregate(per, kind, names)
        out["agg"][kind]["pool"] = {
            k: float(np.mean([r[kind]["pool"][k] for r in per]))
            for k in ("best", "mean", "sd", "n_band", "f2", "f15")}
        out["agg"][kind]["rebuild_displacement"] = float(
            np.mean([r[f"reb_{kind}"] for r in per]))
    write_json("inband_skill.json", out)
    _print_skill(out)
    return out


def _aggregate(per, kind, names):
    agg = {}
    base = np.array([r[kind]["sig"]["score"].get("sel_band", np.nan) for r in per])
    basep = np.array([r[kind]["sig"]["score"].get("sel_pool", np.nan) for r in per])
    for nm in names:
        rows = [r[kind]["sig"][nm] for r in per]
        d = {}
        for key in ("rho_global", "rho_band", "rho_band_prg", "sel_band", "sel_pool",
                    "sel_band_inv", "sel_pool_inv"):
            mu, se, n = _mean_se([x.get(key) for x in rows])
            d[key], d[key + "_se"], d[key + "_n"] = mu, se, n
        v = np.array([x.get("sel_band", np.nan) for x in rows])
        vp = np.array([x.get("sel_pool", np.nan) for x in rows])
        d["vs_score_band"] = _paired(v, base)
        d["vs_score_pool"] = _paired(vp, basep)
        d["f2_band"] = float(np.nanmean(v < 2.0))
        d["f15_band"] = float(np.nanmean(v < 1.5))
        d["f2_pool"] = float(np.nanmean(vp < 2.0))
        d["f15_pool"] = float(np.nanmean(vp < 1.5))
        d["median_band"] = float(np.nanmedian(v))
        d["sd_band"] = float(np.nanstd(v, ddof=1))
        d["worst_band"] = float(np.nanmax(v))
        d["best_band"] = float(np.nanmin(v))
        # SIGN STABILITY of the in-band rho: a sign that flips per target is noise
        rb = np.array([x.get("rho_band", np.nan) for x in rows], float)
        rb = rb[np.isfinite(rb)]
        d["frac_rho_negative"] = float((rb < 0).mean()) if len(rb) else None
        d["rho_band_median"] = float(np.median(rb)) if len(rb) else None
        # per-fold sign, the stability test the brief demands
        folds = np.array([r["fold"] for r in per])
        fr = []
        for f in range(5):
            sel = [r[kind]["sig"][nm].get("rho_band") for r, ff in zip(per, folds)
                   if ff == f]
            mu, se, n = _mean_se(sel)
            fr.append(mu)
        d["rho_band_by_fold"] = fr
        d["sign_stable_folds"] = (
            None if any(x is None for x in fr)
            else bool(all(x > 0 for x in fr) or all(x < 0 for x in fr)))
        agg[nm] = d
    return agg


def _fmt(x, w=7, p=3):
    return f"{x:{w}.{p}f}" if isinstance(x, (int, float)) and x is not None else " " * (w - 1) + "-"


def _print_skill(out):
    names = out["names"]
    for kind, title in (("nrm", "NORMAL pool (BLOSUM top-500)"),
                        ("tgt", "ORACLE-TIGHT pool (best 500 of the universe)")):
        a = out["agg"][kind]
        p = a["pool"]
        print(f"\n=== {title} — {out['n_targets']} targets ===")
        print(f"  pool best {p['best']:.3f}  mean {p['mean']:.3f}  sd {p['sd']:.3f}  "
              f"<2A {p['f2']:.3f}  <1.5A {p['f15']:.3f}  band n {p['n_band']:.1f}  "
              f"rebuild disp {a['rebuild_displacement']:.3f} A")
        print(f"{'signal':18} {'rho_gl':>7} {'rho_bd':>7} {'+-se':>6} {'prg':>7} "
              f"{'selBAND':>8} {'d vs sc':>8} {'W/L':>8} {'selPOOL':>8} {'inv':>7} "
              f"{'neg%':>5} {'stab':>5}")
        order = sorted(names, key=lambda z: (a[z]["sel_band"]
                                             if a[z]["sel_band"] is not None else 9e9))
        for nm in order:
            d = a[nm]
            wl = d["vs_score_band"]
            dv = f"{wl['mean_diff']:+8.3f}" if wl else "       -"
            wlc = f"{wl['n_better']:3d}/{wl['n_worse']:<4d}" if wl else "   -    "
            neg = ("    -" if d["frac_rho_negative"] is None
                   else f"{100 * d['frac_rho_negative']:5.0f}")
            print(f"{nm:18} {_fmt(d['rho_global'])} {_fmt(d['rho_band'])} "
                  f"{_fmt(d['rho_band_se'],6)} {_fmt(d['rho_band_prg'])} "
                  f"{_fmt(d['sel_band'],8)} {dv} {wlc:>8} "
                  f"{_fmt(d['sel_pool'],8)} {_fmt(d['sel_band_inv'],7)} "
                  f"{neg} {str(d['sign_stable_folds']):>5}")


# ============================================================ stage: anti
def stage_anti():
    """Characterise the anti-correlation: where does in-band rho cross zero?

    Two sweeps on identical candidates: band width at fixed pool, and pool tightness at
    fixed band.  Both are ORACLE constructions.  A signal is only usable inverted if its
    sign is stable, so the per-target sign distribution and a paired t against zero are
    reported alongside every mean.
    """
    tgs = cached()
    print(f"anti-correlation characterisation on {len(tgs)} targets", flush=True)
    names = None
    # --- band-width sweep, both pools, all signals
    band_rows = {kind: {b: {} for b in BANDS} for kind in POOLS}
    per = []
    for pdbid in tgs:
        c = load_cache(pdbid)
        names = c["names"]
        for kind in POOLS:
            S, rr, rg = c[f"S_{kind}"], c[f"rr_{kind}"], c[f"rg_{kind}"]
            best = rr.min()
            for b in BANDS:
                m = rr <= best + b
                if m.sum() < MIN_BAND:
                    continue
                for ci, nm in enumerate(names):
                    band_rows[kind][b].setdefault(nm, []).append(
                        (_rho(S[m, ci], rr[m]), sel_of(S[m, ci], rr[m])))
        per.append(pdbid)
    band = {kind: {} for kind in POOLS}
    for kind in POOLS:
        for b in BANDS:
            band[kind][str(b)] = {}
            for nm in names:
                v = band_rows[kind][b].get(nm, [])
                r = [x[0] for x in v]
                s = [x[1] for x in v]
                mu, se, n = _mean_se(r)
                smu, sse, _ = _mean_se(s)
                rr_ = np.array([x for x in r if x is not None], float)
                band[kind][str(b)][nm] = {
                    "rho": mu, "se": se, "n": n, "sel": smu, "sel_se": sse,
                    "frac_neg": float((rr_ < 0).mean()) if len(rr_) else None,
                    "t": (float(mu / se) if (mu is not None and se and se > 0)
                          else None)}

    # --- pool-TIGHTNESS sweep (ORACLE): the top-`m` of the oracle order IS the pool.
    # `m = 0` is the deployable BLOSUM top-500 reference on the same universe.
    # `vs_random` compares selected RMSD against the pool's own MEAN, which is the honest
    # null for the claim "worse than a coin flip".
    tight = {}
    key_sigs = ["score", "l1exp", "zdev", "typic", "lg_all", "lg_hb", "rg", "bagged",
                "rama", "ssfit", "nhb"]
    for m_ in TIGHTNESS:
        acc = {nm: [] for nm in key_sigs}
        pool_stats = []
        for pdbid in tgs:
            u = gen.load_univ(pdbid)
            pred = poolsize.load_pred(pdbid)
            idx = (pool_index(u, "nrm", K) if m_ == 0
                   else pool_index(u, "tgt", min(m_, len(u["rr"]))))
            ctx = deployable_ctx(u, pred, idx)
            S, nm_all = signal_block(ctx)
            rr = np.asarray(u["rr"], float)[idx]
            mk = rr <= rr.min() + BAND
            ok = mk.sum() >= MIN_BAND
            pool_stats.append((float(rr.mean()), float(rr.std(ddof=1)),
                               float((rr < 2.0).mean()), float(rr.min())))
            for nm in key_sigs:
                ci = nm_all.index(nm)
                acc[nm].append((_rho(S[mk, ci], rr[mk]) if ok else None,
                                sel_of(S[mk, ci], rr[mk]) if ok else None,
                                sel_of(S[:, ci], rr), float(rr.mean())))
        ps = np.array(pool_stats)
        tight[str(m_)] = {"pool_mean": float(ps[:, 0].mean()),
                          "pool_sd": float(ps[:, 1].mean()),
                          "pool_f2": float(ps[:, 2].mean()),
                          "pool_best": float(ps[:, 3].mean()), "sig": {}}
        for nm in key_sigs:
            v = acc[nm]
            mu, se, n = _mean_se([x[0] for x in v])
            smu, sse, _ = _mean_se([x[1] for x in v])
            pmu, pse, _ = _mean_se([x[2] for x in v])
            rv = np.array([x[0] for x in v if x[0] is not None], float)
            tight[str(m_)]["sig"][nm] = {
                "rho": mu, "se": se, "n": n, "sel_band": smu, "sel_band_se": sse,
                "sel_pool": pmu, "sel_pool_se": pse,
                "vs_random": _paired([x[2] for x in v], [x[3] for x in v]),
                "frac_neg": float((rv < 0).mean()) if len(rv) else None}
        d = tight[str(m_)]["sig"]["score"]
        lbl = str(m_) if m_ else "BLOSUM500"
        print(f"  pool {lbl:>9}: mean {tight[str(m_)]['pool_mean']:.3f}"
              f"  score rho {d['rho']:+.3f}  sel {d['sel_pool']:.3f}"
              f"  vs random {d['vs_random']['mean_diff']:+.3f}", flush=True)
    out = {"n_targets": len(tgs), "names": names, "band_sweep": band,
           "tightness_sweep": tight, "BANDS": list(BANDS),
           "TIGHTNESS": list(TIGHTNESS)}
    write_json("inband_anti.json", out)
    _print_anti(out)
    return out


def _print_anti(out):
    print(f"\n=== BAND-WIDTH SWEEP: in-band rho vs band width ({out['n_targets']} "
          f"targets) ===")
    show = ["score", "l1exp", "zdev", "typic", "bagged", "lg_all", "lg_hb", "rg",
            "rama", "ssfit", "nhb"]
    for kind, title in (("nrm", "NORMAL pool"), ("tgt", "ORACLE-TIGHT pool")):
        print(f"\n-- {title}: rho (frac of targets negative) --")
        hdr = "  ".join(f"{('all' if b > 100 else f'{b:g}A'):>13}" for b in out["BANDS"])
        print(f"{'signal':16} {hdr}")
        for nm in show:
            cells = []
            for b in out["BANDS"]:
                d = out["band_sweep"][kind][str(b)].get(nm)
                cells.append("            -" if not d or d["rho"] is None else
                             f"{d['rho']:+6.3f}({100*d['frac_neg']:3.0f}%)")
            print(f"{nm:16} " + "  ".join(cells))
    print("\n=== POOL-TIGHTNESS SWEEP (ORACLE): the top-m of the oracle order IS the pool ===")
    print("  `rnd` is selected minus the pool's own MEAN: positive means the signal picks "
          "WORSE than a random draw.")
    cols = ["score", "typic", "lg_all", "rama"]
    print(f"{'pool':>10} {'poolmean':>9} {'poolbest':>9} {'<2A':>6} | " +
          "  ".join(f"{c:>26}" for c in cols))
    for f in out["TIGHTNESS"]:
        t = out["tightness_sweep"][str(f)]
        cells = []
        for c in cols:
            d = t["sig"][c]
            cells.append(f"{d['rho']:+6.3f} sel{d['sel_pool']:6.3f} "
                         f"rnd{d['vs_random']['mean_diff']:+6.3f}")
        print(f"{(f if f else 'BLOSUM500'):>10} {t['pool_mean']:9.3f} "
              f"{t['pool_best']:9.3f} {t['pool_f2']:6.3f} | " +
              "  ".join(f"{c:>26}" for c in cells))


# ============================================================ stage: comb
def _zcols(S):
    mu, sd = S.mean(0), S.std(0)
    return (S - mu) / np.maximum(sd, 1e-9)


def stage_comb(kind="both"):
    """Leave-fold-out weighted combiners over the signal block.

    Sprint 7's ranking-loss experiment was +0.212 A WORSE than plain regression, so the
    first attempt here is the honest simple one: standardise each signal WITHIN a target's
    band (so no per-target scale leaks), fit a linear map to the candidate's true RMSD on
    the training folds, and select by it on the held-out fold.  Two variants: fitted on
    band members only, and on the whole pool.
    """
    tgs = cached()
    kinds = POOLS if kind == "both" else (kind,)
    out = {"n_targets": len(tgs), "arms": {}}
    caches = [load_cache(p) for p in tgs]
    names = caches[0]["names"]
    # a compact, mechanism-diverse subset plus the everything arm
    subsets = {
        "all": names,
        "dist": [n for n in names if n.startswith(("score", "l1exp", "zdev", "sdmatch",
                                                   "bagged", "typic", "rgpred"))],
        "legacy": [n for n in names if n.startswith("lg_")],
        "physgeo": ["nhb", "nhb_lr", "ncontact", "cordr", "strain", "rama", "rama_min",
                    "ssfit", "sshelix", "sssheet", "rg"],
        "score_plus_legacy": ["score", "lg_all", "lg_hb", "lg_burial", "lg_packing"],
        "score_plus_best": ["score", "l1exp", "typic", "bagged", "rg", "rama", "ssfit"],
    }
    for kd in kinds:
        out["arms"][kd] = {}
        for tag, cols in subsets.items():
            ci = [names.index(c) for c in cols]
            for fit_on in ("band", "pool"):
                sel_b, sel_p, rho_b = [], [], []
                for held in range(5):
                    X, Y = [], []
                    for c in caches:
                        if c["fold"] == held:
                            continue
                        S, rr = c[f"S_{kd}"][:, ci], c[f"rr_{kd}"]
                        m = (rr <= rr.min() + BAND) if fit_on == "band" else \
                            np.ones(len(rr), bool)
                        if m.sum() < MIN_BAND:
                            continue
                        X.append(_zcols(S[m]))
                        Y.append(rr[m] - rr[m].mean())
                    if not X:
                        continue
                    X = np.vstack(X)
                    Y = np.concatenate(Y)
                    A = np.hstack([X, np.ones((len(X), 1))])
                    lam = 1e-3 * len(A)
                    G = A.T @ A + lam * np.eye(A.shape[1])
                    beta = np.linalg.solve(G, A.T @ Y)
                    for c in caches:
                        if c["fold"] != held:
                            continue
                        S, rr = c[f"S_{kd}"][:, ci], c[f"rr_{kd}"]
                        m = rr <= rr.min() + BAND
                        v_pool = _zcols(S) @ beta[:-1]
                        sel_p.append(sel_of(v_pool, rr))
                        if m.sum() >= MIN_BAND:
                            v = _zcols(S[m]) @ beta[:-1]
                            sel_b.append(sel_of(v, rr[m]))
                            rho_b.append(_rho(v, rr[m]))
                key = f"{tag}/{fit_on}"
                mb, sb, nb = _mean_se(sel_b)
                mp, sp, npl = _mean_se(sel_p)
                mr, sr, _ = _mean_se(rho_b)
                out["arms"][kd][key] = {
                    "n_cols": len(ci), "sel_band": mb, "sel_band_se": sb, "n_band": nb,
                    "sel_pool": mp, "sel_pool_se": sp, "n_pool": npl, "rho_band": mr}
        # baselines on the same targets
        base_b, base_p = [], []
        for c in caches:
            S, rr = c[f"S_{kd}"], c[f"rr_{kd}"]
            m = rr <= rr.min() + BAND
            v = S[:, names.index("score")]
            base_p.append(sel_of(v, rr))
            if m.sum() >= MIN_BAND:
                base_b.append(sel_of(v[m], rr[m]))
        out["arms"][kd]["SHIPPED score"] = {
            "n_cols": 1,
            "sel_band": _mean_se(base_b)[0], "sel_band_se": _mean_se(base_b)[1],
            "sel_pool": _mean_se(base_p)[0], "sel_pool_se": _mean_se(base_p)[1],
            "n_band": len(base_b), "n_pool": len(base_p), "rho_band": None}
    write_json("inband_comb.json", out)
    _print_comb(out)
    return out


def _print_comb(out):
    for kd in out["arms"]:
        print(f"\n=== LEAVE-FOLD-OUT COMBINERS — "
              f"{'NORMAL' if kd == 'nrm' else 'ORACLE-TIGHT'} pool ===")
        print(f"{'arm':30} {'cols':>5} {'selBAND':>9} {'+-se':>7} {'selPOOL':>9} "
              f"{'+-se':>7} {'rho_bd':>7}")
        rows = sorted(out["arms"][kd].items(),
                      key=lambda kv: (kv[1]["sel_pool"] if kv[1]["sel_pool"]
                                      is not None else 9e9))
        for k, d in rows:
            print(f"{k:30} {d['n_cols']:5d} {_fmt(d['sel_band'],9)} "
                  f"{_fmt(d['sel_band_se'],7)} {_fmt(d['sel_pool'],9)} "
                  f"{_fmt(d['sel_pool_se'],7)} {_fmt(d['rho_band'])}")


# ============================================================ stage: two
#: Filter sizes for the two-stage selector.  `0` means no filter (the whole pool).
QS = (5, 10, 20, 30, 50, 75, 100, 150, 250, 0)


def _medoid_scores(W):
    """Mean CA-RMSD of each candidate to every other. The SUPERPOSITION consensus.

    S8-1 measured that consuming agreement by superposition beats consuming it over pairs
    by 0.28 A at matched information, so if the tight-pool signal really is consensus, the
    superposition medoid should beat the distance-matrix medoid. That is a prediction from
    an existing finding, made before the number was computed.
    """
    B = len(W)
    if B > 400:                      # O(B^2) Kabsch; the filtered subpools are small
        return None
    M = np.zeros((B, B))
    for a in range(B):
        M[a] = audit.kabsch_rmsd_batch(W, W[a])
    return M.sum(1) / max(B - 1, 1)


def stage_two(qs=None):
    """DEPLOYABLE. Two-stage selection: the shipped score FILTERS, consensus SELECTS.

    `stage_anti` measures that the shipped score's in-band skill crosses zero at a pool
    mean of ~2.8 A and is significantly WORSE than a random draw below it, while
    typicality -- agreement with the pool's own mean distance profile, which uses no
    sequence and no model -- is significantly BETTER than random on every tight pool.  The
    two are therefore complementary: the score is a good coarse filter and a bad fine
    ranker, and consensus is the reverse.

    Nothing here reads a native.  The pool is the shipped BLOSUM top-500, the filter is
    the shipped score, and the consensus is computed inside the filtered subpool.  The
    filter size `q` is a hyperparameter, so the headline arm is chosen LEAVE-FOLD-OUT.
    """
    qs = QS if qs is None else tuple(int(x) for x in qs)
    tgs = cached()
    print(f"two-stage selection on {len(tgs)} targets", flush=True)
    per = []
    for k, pdbid in enumerate(tgs):
        u = gen.load_univ(pdbid)
        pred = poolsize.load_pred(pdbid)
        n = int(u["n"])
        idx = pool_index(u, "nrm", K)
        W = np.asarray(u["W"], float)[idx]
        rr = np.asarray(u["rr"], float)[idx]
        D = gen.pair_D(W, n)
        sc = gen.shipped_score(pred, D)
        order = np.argsort(sc, kind="stable")
        row = {"pdb": pdbid, "fold": int(u["fold"]), "best": float(rr.min()),
               "mean": float(rr.mean()), "score": sel_of(sc, rr), "arms": {}}
        zs = (sc - sc.mean()) / max(sc.std(), 1e-12)
        for q in qs:
            sub = order[:q] if q else np.arange(len(rr))
            Ds, Ws, rs = D[sub], W[sub], rr[sub]
            row["arms"][f"typic{q}"] = sel_of(
                np.abs(Ds - Ds.mean(0)[None, :]).mean(1), rs)
            md = _medoid_scores(Ws)
            row["arms"][f"medoid{q}"] = (None if md is None else sel_of(md, rs))
            row["arms"][f"submean{q}"] = float(rs.mean())
        # blends on the FULL pool: z(score) + w * z(typicality)
        ty = np.abs(D - D.mean(0)[None, :]).mean(1)
        zt = (ty - ty.mean()) / max(ty.std(), 1e-12)
        for w in (0.25, 0.5, 1.0, 2.0):
            row["arms"][f"blend{w}"] = sel_of(zs + w * zt, rr)
        per.append(row)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    arms = sorted({a for r in per for a in r["arms"]})
    base = np.array([r["score"] for r in per])
    out = {"n_targets": len(per), "K": K, "qs": list(qs), "per_target": per,
           "base_score": float(base.mean()),
           "base_se": float(base.std(ddof=1) / math.sqrt(len(base))),
           "pool_best": float(np.mean([r["best"] for r in per])),
           "pool_mean": float(np.mean([r["mean"] for r in per])), "arms": {}}
    for a in arms:
        v = np.array([r["arms"].get(a, np.nan) for r in per], float)
        if not np.isfinite(v).all():
            continue
        out["arms"][a] = {
            "sel": float(v.mean()),
            "se": float(v.std(ddof=1) / math.sqrt(len(v))),
            "median": float(np.median(v)), "sd": float(v.std(ddof=1)),
            "worst": float(v.max()), "best": float(v.min()),
            "f2": float((v < 2.0).mean()), "f15": float((v < 1.5).mean()),
            "vs_score": _paired(v, base)}
    # Leave-fold-out choice of the filter size, so no target contributes to picking the
    # arm it is then scored under.  `v` MUST stay aligned with `per` (and therefore with
    # `base`): the first version of this loop appended in fold order and silently paired
    # each held-out value against the wrong target's baseline, which inflated the CI from
    # [-0.316, -0.027] to [-0.608, +0.265] on an identical mean.
    folds = np.array([r["fold"] for r in per])
    fams = {"typic": lambda a: a.startswith("typic") and a[5:].isdigit(),
            "medoid": lambda a: a.startswith("medoid") and a[6:].isdigit(),
            "consensus": lambda a: ((a.startswith("typic") and a[5:].isdigit())
                                    or (a.startswith("medoid") and a[6:].isdigit()))}
    for fam, pick in fams.items():
        cand = [a for a in out["arms"] if pick(a)]
        if not cand:
            continue
        v = np.full(len(per), np.nan)
        picked = []
        for f in range(5):
            tr = folds != f
            if not tr.any() or not (folds == f).any():
                continue
            best_a = min(cand, key=lambda a: float(np.mean(
                [r["arms"][a] for r, t in zip(per, tr) if t])))
            picked.append(f"fold{f}:{best_a}")
            for k2, r in enumerate(per):
                if r["fold"] == f:
                    v[k2] = r["arms"][best_a]
        assert np.isfinite(v).all(), "leave-fold-out coverage is incomplete"
        out["arms"][f"{fam}_LFO"] = {
            "sel": float(v.mean()), "se": float(v.std(ddof=1) / math.sqrt(len(v))),
            "median": float(np.median(v)), "sd": float(v.std(ddof=1)),
            "worst": float(v.max()), "best": float(v.min()),
            "f2": float((v < 2.0).mean()), "f15": float((v < 1.5).mean()),
            "picked_per_fold": picked, "vs_score": _paired(v, base)}
    write_json("inband_two.json", out)
    _print_two(out)
    return out


def _print_two(out):
    print(f"\n=== TWO-STAGE SELECTION (DEPLOYABLE) — {out['n_targets']} targets, "
          f"K={out['K']} ===")
    print(f"  pool best {out['pool_best']:.3f}  pool mean {out['pool_mean']:.3f}  "
          f"SHIPPED score {out['base_score']:.3f} +- {out['base_se']:.3f}")
    print(f"{'arm':16} {'sel':>7} {'+-se':>6} {'med':>7} {'sd':>6} {'<2A':>6} "
          f"{'<1.5A':>6} {'d vs score':>11} {'95% CI':>18} {'W/L':>9}")
    rows = sorted(out["arms"].items(), key=lambda kv: kv[1]["sel"])
    for a, d in rows:
        if a.startswith("submean"):
            continue
        p = d["vs_score"]
        print(f"{a:16} {d['sel']:7.3f} {d['se']:6.3f} {d['median']:7.3f} "
              f"{d['sd']:6.3f} {d['f2']:6.3f} {d['f15']:6.3f} "
              f"{p['mean_diff']:+11.3f} "
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}]".rjust(19) +
              f" {p['n_better']:3d}/{p['n_worse']:<4d}")
    print("  submeanQ (the score-filtered subpool's own MEAN, i.e. what a random draw "
          "from it gives):")
    print("   " + "  ".join(
        f"q={a[7:]}:{d['sel']:.3f}" for a, d in sorted(
            ((a, d) for a, d in out["arms"].items() if a.startswith("submean")),
            key=lambda kv: int(kv[0][7:]))))


# ============================================================ stage: ens
#: The ensemble members.  `s7/repr_models` holds one seed per arm, so these differ by
#: REPRESENTATION rather than by seed.  That is the favourable direction for the
#: disagreement hypothesis -- representation-diverse models err more independently than
#: re-seeded copies of one -- so a null here is the stronger result.  `raw` and `pca128`
#: are excluded: their encoder loads a 243 MB ESM cache and the box is shared.
ENS_ARMS = ("onehot", "phys", "pca32")


def _ens_preds(seq, fold):
    """One `Distogram` per ensemble arm for this sequence. DEPLOYABLE: sequence only."""
    from s7 import repr_select as rs
    out = []
    for arm in ENS_ARMS:
        if not os.path.exists(rs.model_path(arm, 0, fold)):
            continue
        net, enc, mu, sd = rs.train_one(arm, 0, fold, verbose=False)
        out.append(rs.predict(seq, net, enc, mu, sd))
    return out


def stage_ens():
    """Does ENSEMBLE DISAGREEMENT discriminate inside the near-native band?

    Four signals per candidate, all native-free:

      ens_mean   mean of the members' pool-standardised scores  (an ensemble SELECTOR)
      ens_sd     spread of those scores                          (pure DISAGREEMENT)
      ens_pen    ens_mean + ens_sd                               (penalise disagreement)
      ens_conf   L1 against the ensemble-mean predicted matrix, restricted to the pairs
                 the members AGREE on most (the lowest-spread half)

    Scored on both pools against the shipped distogram on identical candidates.
    """
    tgs = cached()
    print(f"ensemble disagreement ({len(ENS_ARMS)} arms) on {len(tgs)} targets", flush=True)
    per = []
    names = ["ens_mean", "ens_sd", "ens_pen", "ens_conf"]
    for k, pdbid in enumerate(tgs):
        u = gen.load_univ(pdbid)
        pred = poolsize.load_pred(pdbid)
        n, fold = int(u["n"]), int(u["fold"])
        dgs = _ens_preds(str(u["seq"]), fold)
        if len(dgs) < 2:
            continue
        row = {"pdb": pdbid, "fold": fold, "n_members": len(dgs), "pools": {}}
        for kind in POOLS:
            idx = pool_index(u, kind)
            W = np.asarray(u["W"], float)[idx]
            rr = np.asarray(u["rr"], float)[idx]
            D = gen.pair_D(W, n)
            Z = []
            for d in dgs:
                v = d.score(W)
                Z.append((v - v.mean()) / max(v.std(), 1e-12))
            Z = np.stack(Z)                                   # (M, B)
            E = np.stack([d.expected for d in dgs])           # (M, npairs)
            spread = E.std(0)
            lo = spread <= np.median(spread)
            sig = {"ens_mean": Z.mean(0), "ens_sd": Z.std(0),
                   "ens_pen": Z.mean(0) + Z.std(0),
                   "ens_conf": np.abs(D[:, lo] - E.mean(0)[None, lo]).mean(1)}
            sc = gen.shipped_score(pred, D)
            rg = gen.rg_of(D, n)
            m = rr <= rr.min() + BAND
            cell = {"score": {"sel_pool": sel_of(sc, rr)}}
            if m.sum() >= MIN_BAND:
                cell["score"]["sel_band"] = sel_of(sc[m], rr[m])
                cell["score"]["rho_band"] = _rho(sc[m], rr[m])
                cell["score"]["rho_band_prg"] = partial_rho(sc[m], rr[m], rg[m])
            for nm, v in sig.items():
                c = {"sel_pool": sel_of(v, rr), "rho_global": _rho(v, rr)}
                if m.sum() >= MIN_BAND:
                    c["sel_band"] = sel_of(v[m], rr[m])
                    c["sel_band_inv"] = sel_of(v[m], rr[m], invert=True)
                    c["rho_band"] = _rho(v[m], rr[m])
                    c["rho_band_prg"] = partial_rho(v[m], rr[m], rg[m])
                cell[nm] = c
            row["pools"][kind] = cell
        per.append(row)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/{len(tgs)}", flush=True)
    out = {"n_targets": len(per), "arms": list(ENS_ARMS), "names": names, "agg": {},
           "per_target": per}
    for kind in POOLS:
        base = np.array([r["pools"][kind]["score"].get("sel_band", np.nan) for r in per])
        basep = np.array([r["pools"][kind]["score"]["sel_pool"] for r in per])
        out["agg"][kind] = {}
        for nm in ["score"] + names:
            rows = [r["pools"][kind][nm] for r in per]
            d = {}
            for key in ("rho_global", "rho_band", "rho_band_prg", "sel_band",
                        "sel_pool", "sel_band_inv"):
                mu, se, nn = _mean_se([x.get(key) for x in rows])
                d[key], d[key + "_se"] = mu, se
            rb = np.array([x.get("rho_band", np.nan) for x in rows], float)
            rb = rb[np.isfinite(rb)]
            d["frac_rho_negative"] = float((rb < 0).mean()) if len(rb) else None
            d["vs_score_band"] = _paired(
                [x.get("sel_band", np.nan) for x in rows], base)
            d["vs_score_pool"] = _paired(
                [x["sel_pool"] for x in rows], basep)
            out["agg"][kind][nm] = d
    write_json("inband_ens.json", out)
    _print_ens(out)
    return out


def _print_ens(out):
    for kind, title in (("nrm", "NORMAL pool"), ("tgt", "ORACLE-TIGHT pool")):
        print(f"\n=== ENSEMBLE DISAGREEMENT, {title} — {out['n_targets']} targets, "
              f"members {out['arms']} ===")
        print(f"{'signal':12} {'rho_gl':>7} {'rho_bd':>7} {'prg':>7} {'selBAND':>8} "
              f"{'d band':>8} {'selPOOL':>8} {'d pool':>8} {'W/L pool':>9} {'inv':>7} "
              f"{'neg%':>5}")
        for nm, d in out["agg"][kind].items():
            wp = d["vs_score_pool"]
            wb = d["vs_score_band"]
            neg = ("    -" if d["frac_rho_negative"] is None
                   else f"{100 * d['frac_rho_negative']:5.0f}")
            wlc = "    -" if not wp else f"{wp['n_better']:3d}/{wp['n_worse']:<4d}"
            print(f"{nm:12} {_fmt(d['rho_global'])} {_fmt(d['rho_band'])} "
                  f"{_fmt(d['rho_band_prg'])} {_fmt(d['sel_band'],8)} "
                  f"{_fmt(wb['mean_diff'] if wb else None,8)} "
                  f"{_fmt(d['sel_pool'],8)} "
                  f"{_fmt(wp['mean_diff'] if wp else None,8)} "
                  f"{wlc:>9} {_fmt(d['sel_band_inv'],7)} {neg}")


# ============================================================ stage: dev
DCACHE = os.path.join(HERE, "inband_devpool")

#: The arm carried to the dev set, fixed by `stage_two`'s leave-fold-out result BEFORE the
#: dev set was touched: filter the shipped pool to the score's top 75, take the consensus
#: medoid.  All five tuning folds chose it independently.  ONE arm, ONE pass, no iteration.
DEV_ARM = ("medoid", 75)


def stage_devpool(verbose=True):
    """Top-K BLOSUM windows WITH COORDINATES for the 24 dev targets.

    `s7/audit_cache` stores each dev pool's distance matrix and RMSD but not its
    coordinates, and the medoid consensus needs coordinates.  This rebuilds the identical
    pool through `audit`'s own code path and ASSERTS it reproduces the cached `D` and `rr`,
    so the dev pass is provably on the same instrument as every other dev number in the
    project.  Nothing here is tuned; it is a cache.
    """
    import distogram as dgm
    os.makedirs(DCACHE, exist_ok=True)
    folds = db.folds(5)
    tg = [q for q in debias.dev_targets()
          if not os.path.exists(os.path.join(DCACHE, f"{q.pdb}.npz"))]
    if verbose:
        print(f"{len(tg)}/24 dev pools to build", flush=True)
    if not tg:
        return
    memo = {}
    _orig = dgm._fold_fragments

    def frags(fold, n_folds, threshold=db.IDENTITY_THRESHOLD):
        key = (fold, n_folds, threshold)
        if key not in memo:
            memo[key] = _orig(fold, n_folds, threshold)
        return memo[key]

    dgm._fold_fragments = frags
    worst_d = worst_r = 0.0
    try:
        for k, q in enumerate(sorted(tg, key=lambda x: folds[x.seq])):
            t0 = time.time()
            fold, peps, fr = audit.build_pool_members(q, folds)
            W, S, src = audit.windows_of(peps + fr, q.n)
            sim = audit.B62[S, audit.encode(q.seq)[None, :]].sum(1)
            idx = np.argsort(-sim, kind="stable")[:K]
            i, j = audit.pair_index(q.n)
            D = audit.pair_dists(W[idx], i, j)
            rr = audit.kabsch_rmsd_batch(W[idx], q.ca)
            rec = audit.load_rec(q.pdb)
            ed = float(np.abs(D - rec["D"][:K]).max())
            er = float(np.abs(rr - rec["rr"][:K]).max())
            assert ed < 1e-3 and er < 1e-3, f"{q.pdb}: dev pool mismatch {ed} {er}"
            worst_d, worst_r = max(worst_d, ed), max(worst_r, er)
            np.savez_compressed(
                os.path.join(DCACHE, f"{q.pdb}.npz"), pdb=q.pdb, n=q.n, fold=fold,
                seq=q.seq, W=W[idx].astype(np.float32), rr=rr.astype(np.float32))
            if verbose:
                print(f"[{k+1:2d}/{len(tg)}] {q.pdb} n={q.n} fold={fold} "
                      f"dD={ed:.1e} drr={er:.1e} ({time.time()-t0:.0f}s)", flush=True)
    finally:
        dgm._fold_fragments = _orig
    if verbose:
        print(f"worst |dD|={worst_d:.1e}  worst |drr|={worst_r:.1e}", flush=True)


def stage_dev():
    """THE SINGLE DEV PASS. One pre-registered arm, no iteration, no selection here.

    The arm is `DEV_ARM`, fixed by the tuning instrument's leave-fold-out result before
    this function was ever run.  The shipped distogram score on the identical pools is the
    comparator.  A secondary row reports the distance-matrix variant of the same consensus,
    clearly labelled as NOT pre-registered.
    """
    stage_devpool()
    fam, q = DEV_ARM
    rows = []
    for pdbid in sorted(f[:-4] for f in os.listdir(DCACHE) if f.endswith(".npz")):
        z = np.load(os.path.join(DCACHE, f"{pdbid}.npz"), allow_pickle=True)
        W = z["W"].astype(float)
        rr = z["rr"].astype(float)
        n = int(z["n"])
        pred = poolsize.load_pred(pdbid)
        D = gen.pair_D(W, n)
        sc = gen.shipped_score(pred, D)
        sub = np.argsort(sc, kind="stable")[:q]
        md = _medoid_scores(W[sub])
        ty = np.abs(D[sub] - D[sub].mean(0)[None, :]).mean(1)
        rows.append({"pdb": pdbid, "n": n, "fold": int(z["fold"]),
                     "best": float(rr.min()), "mean": float(rr.mean()),
                     "score": sel_of(sc, rr),
                     f"{fam}{q}": sel_of(md, rr[sub]),
                     f"typic{q}": sel_of(ty, rr[sub])})
    base = np.array([r["score"] for r in rows])
    out = {"n_targets": len(rows), "K": K, "arm": f"{fam}{q}",
           "pool_best": float(np.mean([r["best"] for r in rows])),
           "pool_mean": float(np.mean([r["mean"] for r in rows])),
           "per_target": rows, "summary": {}}
    for a in ("score", f"{fam}{q}", f"typic{q}"):
        v = np.array([r[a] for r in rows], float)
        out["summary"][a] = {
            "sel": float(v.mean()), "se": float(v.std(ddof=1) / math.sqrt(len(v))),
            "median": float(np.median(v)), "sd": float(v.std(ddof=1)),
            "worst": float(v.max()), "best": float(v.min()),
            "f2": float((v < 2.0).mean()), "f15": float((v < 1.5).mean()),
            "vs_score": _paired(v, base)}
    write_json("inband_dev.json", out)
    _print_dev(out)
    return out


def _print_dev(out):
    print(f"\n=== THE SINGLE DEV PASS — {out['n_targets']} targets, K={out['K']} ===")
    print(f"  pre-registered arm: {out['arm']}  (chosen leave-fold-out on the 126-target "
          f"tuning set, before this ran)")
    print(f"  pool best {out['pool_best']:.3f}  pool mean {out['pool_mean']:.3f}")
    print(f"{'arm':14} {'sel':>7} {'+-se':>6} {'med':>7} {'<2A':>6} {'<1.5A':>6} "
          f"{'d vs score':>11} {'95% CI':>19} {'W/L':>9}")
    for a, d in out["summary"].items():
        p2 = d["vs_score"]
        ci = f"[{p2['ci95'][0]:+.3f},{p2['ci95'][1]:+.3f}]"
        print(f"{a:14} {d['sel']:7.3f} {d['se']:6.3f} {d['median']:7.3f} "
              f"{d['f2']:6.3f} {d['f15']:6.3f} {p2['mean_diff']:+11.3f} {ci:>19} "
              f"{p2['n_better']:3d}/{p2['n_worse']:<4d}")


# ============================================================ stage: leak
def stage_leak():
    """Every deployable signal must be unchanged when the natives are NaN.

    The band and the tight pool are oracle constructions, so the poison test fixes the
    candidate INDEX SET first and then re-runs `signal_block` on a universe whose `rr` and
    `nat_ca` have been replaced with NaN.  Any signal that had read a native would come
    back NaN or different.
    """
    tgs = cached()[:12]
    worst, bad = 0.0, []
    for pdbid in tgs:
        u = gen.load_univ(pdbid)
        pred = poolsize.load_pred(pdbid)
        for kind in POOLS:
            idx = pool_index(u, kind)
            A, names = signal_block(deployable_ctx(u, pred, idx))
            v = dict(u)
            v["rr"] = np.full_like(u["rr"], np.nan)
            v["nat_ca"] = np.full_like(u["nat_ca"], np.nan)
            B, _ = signal_block(deployable_ctx(v, pred, idx))
            d = np.abs(A - B)
            if not np.isfinite(B).all():
                bad.append((pdbid, kind, "nan"))
            if d.max() > 0:
                bad.append((pdbid, kind, float(d.max())))
            worst = max(worst, float(d.max()))
    out = {"n_targets": len(tgs), "n_signals": len(names), "worst_abs_diff": worst,
           "violations": bad, "clean": not bad}
    write_json("inband_leak.json", out)
    print(f"leakage poison test: {len(tgs)} targets x {len(POOLS)} pools x "
          f"{len(names)} signals; worst |diff| = {worst:.3e}; "
          f"{'CLEAN' if not bad else f'{len(bad)} VIOLATIONS'}")
    assert not bad, bad
    return out


# ============================================================ stage: report
def stage_report():
    for name, fn in (("inband_skill.json", _print_skill),
                     ("inband_anti.json", _print_anti),
                     ("inband_two.json", _print_two),
                     ("inband_comb.json", _print_comb)):
        try:
            fn(read_json(name))
        except FileNotFoundError:
            print(f"({name} not computed)")


STAGES = {"build": stage_build, "skill": stage_skill, "anti": stage_anti,
          "comb": stage_comb, "two": stage_two, "devpool": stage_devpool,
          "dev": stage_dev, "ens": stage_ens, "leak": stage_leak,
          "report": stage_report}


def main(argv):
    if not argv or argv[0] not in STAGES:
        print(__doc__)
        print("stages:", " ".join(STAGES))
        return 1
    STAGES[argv[0]](*argv[1:])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
