"""Do the recognition channels make INDEPENDENT errors?  And if so, does fusing them help?

THE QUESTION
------------
Sprint 8 closed representation (S8-1), generation (S8-5) and retrieval (S8-6).  What is
left is recognition, and every channel has been tested ALONE:

    channel                                   native pct   argmin at native
    learned distance objective (distogram)        36.8         3/126
    all-atom Amber, single point, interaction     54.3         0/70
    inverse folding (s8/invfold.py)               30.1         6/126
    Legacy field (11 decomposed terms)            -- never measured this way --

None of them puts its minimum at the native.  But a combination CAN rank the native first
even when no member does -- provided the members fail on DIFFERENT candidates.  If their
errors are correlated, no combination helps, and that closes recognition.

So the first thing this module measures, before any combined scorer is built, is the ERROR
CORRELATION STRUCTURE: for every pair of channels, on byte-identical candidates, the rank
correlation of their scores AND of their errors (score rank minus true-RMSD rank).

THE INSTRUMENT
--------------
The 126-target tuning set on the cached window universes `s8/generate_univ/*.npz`, pool =
the shipped BLOSUM top-500.  Candidates are the real length-n CA windows; `rr` is each
window's CA-RMSD to the native.  This reproduces the shipped baseline exactly:

    selected 3.4540   pool best 1.7108   pool mean 4.4533

asserted in `test_integrate.py`.  The fast channels (distogram family, Legacy's 11 terms
and 4 subsets, geometry/physics) are read from `s8/inband_cache`, which is built by
`s8/inband.py` and is byte-identical to what this module would recompute -- a test asserts
the reproduction.  The two channels that cache does NOT contain, inverse folding and
Amber, are computed here.

AMBER IS SUBSAMPLED, AND THE SUBSAMPLE IS NATIVE-FREE
-----------------------------------------------------
One single-point ff14SB/GBn2 evaluation is ~5-10 s, so 126 x 500 is out of reach on a
contended box.  Amber is therefore evaluated on the first `M_AMBER` candidates in BLOSUM
order -- a subset chosen with no reference to any native and, per S8-5, essentially
uncorrelated with candidate quality (BLOSUM rho vs true RMSD = +0.066).  Every
cross-channel number involving Amber is computed on exactly that subset for ALL channels,
so the comparison is like-for-like.

DIAGNOSTIC vs DEPLOYABLE
------------------------
Every SCORE in this module is native-free.  Natives enter only as (a) `rr`, the label the
correlations are measured against, (b) the reported CA-RMSD of a selection already made,
and (c) arms explicitly named `oracle_*`.  `stage_leak` NaN-poisons the native in every
cached universe and asserts every deployable channel is bit-identical.

Stages, in dependency order.  All are resumable; the three caches live on disk.

    chan       per-target channel matrix cache            (needs s8/inband_cache)
    dmat       per-target 500x500 pairwise CA-RMSD cache
    amber      the Amber channel, ~5 min/target, resumable per target
    native     the native's own score under every channel (ORACLE diagnostic)

    corr       THE ERROR CORRELATION MATRIX                -- decides everything else
    fuse       25 multi-objective scores, leave-fold-out
    filter     15 filter scores x 11 filter sizes in front of the consensus medoid
    critfuse   fusing the channels with the CONSENSUS CRITERION rather than the score
    k25        the same questions in a 25-candidate pool
    natcons    where the NATIVE sits under the consensus criterion (ORACLE diagnostic)

    cvarcheck  the CVaR gradient against two independent references
    vqe        CVaR-VQE over the discrete hypothesis set + its classical ablations
    role       Legacy and AMBER in the roles the sprint's data suggests for them
    refine     an AMBER relaxation of the structure the system selects
    final      the integrated end-to-end system, one ablation per component

The NaN-poison leakage assertion lives in `s8/test_integrate.py` rather than in a stage
here, because it has to re-derive the deployable channels to compare against the cache.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np                                                          # noqa: E402

os.environ.setdefault("NT", "2")
try:
    import torch
    torch.set_num_threads(2)
except Exception:                                                           # pragma: no cover
    torch = None

import protein_geometry as geo                                              # noqa: E402
from s8 import generate as gen                                              # noqa: E402
from s8 import inband as ib                                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ACACHE = os.path.join(HERE, "integrate_amber")
CCACHE = os.path.join(HERE, "integrate_chan")

K = 500                #: pool size, the number every sprint-7/8 result reports on
BAND = 1.5             #: in-band = within pool_best + BAND
M_AMBER = int(os.environ.get("M_AMBER", "40"))   #: Amber candidates per target
MIN_BAND = 5


# ============================================================ pool reconstruction
def pool(u, k=K):
    """The shipped BLOSUM top-`k` of one cached universe, with its rebuild.

    `W` are the real windows (whose RMSD is the outcome); `BB` is the ideal-geometry
    rebuild from the same torsions, which is what any all-atom or torsion-reading channel
    must be scored on.  The end-torsion placeholders are forced to the builder defaults
    for the SAME reason `s8/invfold.build_pool` does it: `geo.extract_torsions` fills
    phi[0] / psi[n-1] with defaults for a peptide-DB native but a window sliced from a
    longer parent keeps the parent's real values, so leaving them alone would let a
    torsion-reading channel discriminate on provenance rather than conformation.
    """
    idx = np.asarray(u["order"], int)[:k]
    PHI = np.asarray(u["PHI"], float)[idx].copy()
    PSI = np.asarray(u["PSI"], float)[idx].copy()
    PHI[:, 0] = geo.DEFAULT_PHI
    PSI[:, -1] = geo.DEFAULT_PSI
    return {"idx": idx, "W": np.asarray(u["W"], float)[idx], "PHI": PHI, "PSI": PSI,
            "BB": geo.build_backbone_batch(PHI, PSI),
            "rr": np.asarray(u["rr"], float)[idx]}


def targets():
    return gen.cached_targets()


# ============================================================ ranking helpers
def _rk(x):
    from scipy.stats import rankdata
    return rankdata(np.asarray(x, float))


def _z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / max(x.std(), 1e-12)


def zrank(x):
    """Standardised rank -- the common currency every channel is converted to."""
    return _z(_rk(x))


def sel_of(score, rr):
    return float(rr[int(np.argmin(np.asarray(score, float)))])


def _spearman(a, b):
    from scipy.stats import spearmanr
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 4 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return None
    v = spearmanr(a, b).statistic
    return None if not np.isfinite(v) else float(v)


def partial_rg(score, rr, rg):
    """Spearman(score, rmsd) with radius of gyration partialled out.

    MANDATORY in this project: sprint 7 finding 10 recorded that Amber's apparent in-band
    skill was carried by compactness, and finding 1 that a generator-agreement lead went
    +0.594 -> -0.07 under exactly this correction.  Same formula as
    `s7/dist_baseline.partial_rg` and `s8/invfold.partial_rg`.
    """
    a, r, g = _z(_rk(score)), _z(_rk(rr)), _z(_rk(rg))
    p_ar, p_ag, p_rg = float(np.mean(a * r)), float(np.mean(a * g)), float(np.mean(r * g))
    den = np.sqrt(max((1 - p_ag ** 2) * (1 - p_rg ** 2), 1e-12))
    return float((p_ar - p_ag * p_rg) / den)


def rg_of(CA):
    CA = np.asarray(CA, float)
    return np.sqrt(((CA - CA.mean(1, keepdims=True)) ** 2).sum(2).mean(1))


# ============================================================ the channel matrix
#: The channels the sprint has to reconcile.  Every one is a SCORE: lower is the
#: prediction.  `dist` is the shipped selector; `lg_*` come from `legacy_field`'s
#: decomposition; `iv_*` are the inverse-folding self-consistency scores; `amber` is the
#: converged interaction-only all-atom energy.  Names in `FAST` are read from
#: `s8/inband_cache`; `iv_*` are computed here; `amber` has its own slow stage.
FAST = ("score", "l1exp", "zdev", "bagged", "typic", "rg",
        "lg_all", "lg_hb", "lg_burial", "lg_packing",
        "lg_steric", "lg_contact", "lg_hbond_local", "lg_hbond_longrange",
        "lg_coop_helix", "lg_coop_sheet", "lg_solvation", "lg_electrostatic",
        "lg_aromatic", "lg_torsion", "lg_compactness",
        "nhb", "nhb_lr", "ncontact", "cordr", "strain", "rama", "rama_min",
        "ssfit", "sshelix", "sssheet")
IV = ("iv_geom", "iv_nbid", "iv_geom200k", "iv_geom_ctr", "iv_nbid_ctr")

#: The four MECHANISM channels the error-correlation question is actually about.  Every
#: other column is a variant of one of these and would inflate the matrix with
#: near-duplicates.
#: `typic` is a fifth channel the brief did not list and S8-8 found afterwards: it is the
#: only one of 38 native-free signals with positive in-band skill, and it carries no
#: sequence and no model at all.  A channel that shares neither a training set nor an input
#: with the other four is the strongest available test of the decorrelation hypothesis, so
#: leaving it out would have made the matrix weaker than the data allows.
PRIMARY = ("dist", "legacy", "invfold", "typic", "amber")
PRIMARY_COL = {"dist": "score", "legacy": "lg_all", "invfold": "iv_nbid",
               "typic": "typic", "amber": "amber"}


def _cpath(pdbid):
    return os.path.join(CCACHE, f"{pdbid}.npz")


def stage_chan(limit=None, verbose=True):
    """Cache the per-candidate channel matrix for every target: FAST + IV + rr + rg.

    FAST columns are lifted from `s8/inband_cache` rather than recomputed.  That cache is
    built on the SAME pool (`order[:500]` of the same universe) and `test_integrate.py`
    asserts the `rr` vectors agree bit for bit, so lifting them is a reproduction, not an
    approximation -- and it keeps two agents from spending the box on the same distogram
    and Legacy evaluations.
    """
    from s8 import invfold as ivf
    os.makedirs(CCACHE, exist_ok=True)
    tgs = targets()
    if limit:
        tgs = tgs[:int(limit)]
    todo = [p for p in tgs if not os.path.exists(_cpath(p))]
    if verbose:
        print(f"chan: {len(todo)}/{len(tgs)} targets to do", flush=True)
    models = {v: ivf.load_model(v) for v in ("geom", "nbid", "geom200k")}
    for k, pid in enumerate(todo):
        t0 = time.time()
        u = gen.load_univ(pid)
        P = pool(u)
        c = ib.load_cache(pid)
        assert np.allclose(c["rr_nrm"], P["rr"], atol=1e-4), pid
        names = list(FAST)
        cols = [c["S_nrm"][:, c["names"].index(n)] for n in FAST]
        for v in ("geom", "nbid", "geom200k"):
            sc = ivf.SelfConsistency(u["seq"], models[v], nbid=ivf.VARIANTS[v][0])
            cols.append(sc.score(P["BB"], P["PHI"], P["PSI"]))
            names.append("iv_" + v)
            if v in ("geom", "nbid"):
                cols.append(sc.score_centred(P["BB"], P["PHI"], P["PSI"]))
                names.append("iv_%s_ctr" % v)
        S = np.stack(cols, 1).astype(np.float32)
        tmp = _cpath(pid) + ".tmp.npz"
        np.savez_compressed(tmp, pdb=pid, n=u["n"], fold=u["fold"], seq=u["seq"],
                            S=S, names=np.array(names),
                            rr=P["rr"].astype(np.float32),
                            rg=rg_of(P["W"]).astype(np.float32))
        os.replace(tmp, _cpath(pid))
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {pid} ({time.time()-t0:.1f}s)", flush=True)


def load_chan(pdbid, with_amber=False):
    z = np.load(_cpath(pdbid), allow_pickle=True)
    out = {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
           "seq": str(z["seq"]), "names": [str(x) for x in z["names"]],
           "S": z["S"].astype(float), "rr": z["rr"].astype(float),
           "rg": z["rg"].astype(float)}
    if with_amber and os.path.exists(_apath(pdbid)):
        a = load_amber(pdbid)
        assert np.allclose(a["rr"], out["rr"][:a["m"]], atol=1e-4), pdbid
        col = np.full(len(out["rr"]), np.nan)
        col[:a["m"]] = a["e"]
        out["S"] = np.hstack([out["S"], col[:, None]])
        out["names"] = out["names"] + ["amber"]
        out["amber_m"] = a["m"]
        out["amber_native"] = a["e_native"]
    return out


def chan_cached():
    if not os.path.isdir(CCACHE):
        return []
    return sorted(f[:-4] for f in os.listdir(CCACHE) if f.endswith(".npz"))


# ============================================================ the Amber channel
def _amber_energy(seq, rep, coords_single, ar):
    """Interaction-only ff14SB/GBn2 energy of ONE structure, or nan on strain failure.

    Identical protocol to `s7/amber_native.amber_energy`, which is the form the project
    recorded as the right way to score Amber (converged under moderate positional
    restraints, `nonbonded + solvation` only -- the bonded terms are dominated by
    ideal-geometry rebuild artefacts and carry no conformational information).  A
    candidate whose bond+angle energy exceeds 1000 kcal/mol did not relax and is recorded
    as a strain rejection rather than silently scored.
    """
    out = ar.refine_coords(seq, rep, coords_single, k_restraint=ar.K_MODERATE,
                           steps=0, tolerance=1.0, components=True)
    cm = out["components"]
    if cm["bond"] + cm["angle"] > 1000.0:
        return float("nan"), True
    return float(cm["nonbonded"] + cm["solvation"]), False


def _apath(pdbid):
    return os.path.join(ACACHE, f"{pdbid}.npz")


def stage_amber(m=None, limit=None, verbose=True):
    """Amber on the first `m` candidates in BLOSUM order, per target, resumable.

    The subsample is NATIVE-FREE by construction: it is a prefix of the shipped retrieval
    order.  S8-5 measured BLOSUM's own rank correlation with candidate quality at +0.066,
    so this prefix is close to an unbiased sample of the pool -- which is what makes it a
    legitimate common candidate set for the cross-channel correlations.
    """
    import amber_refine as ar
    import torsion_lib2 as tl2
    m = M_AMBER if m is None else int(m)
    os.makedirs(ACACHE, exist_ok=True)
    tgs = targets()
    if limit:
        tgs = tgs[:int(limit)]
    todo = [p for p in tgs if not os.path.exists(_apath(p))]
    if verbose:
        print(f"amber: {len(todo)}/{len(tgs)} targets to do, m={m}", flush=True)
    for k, pid in enumerate(todo):
        t0 = time.time()
        u = gen.load_univ(pid)
        P = pool(u)
        seq = u["seq"]
        tab = tl2.library_for(seq, 8, seq)          # excludes self
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        e = np.full(m, np.nan)
        strain = 0
        for b in range(m):
            try:
                e[b], bad = _amber_energy(seq, rep, {kk: v[b] for kk, v in P["BB"].items()}, ar)
                strain += int(bad)
            except Exception:
                pass
        # ORACLE DIAGNOSTIC: the native's own rebuilt energy, for the native-percentile
        # statistic only.  It is never a feature and never enters a ranking.
        e_nat = np.nan
        try:
            import peptide_db as db
            tp = db.by_pdb(pid)
            nb = geo.build_backbone(np.asarray(tp.phi, float), np.asarray(tp.psi, float))
            e_nat, _ = _amber_energy(seq, rep, {kk: np.asarray(v) for kk, v in nb.items()}, ar)
        except Exception:
            pass
        tmp = _apath(pid) + ".tmp.npz"
        np.savez_compressed(tmp, pdb=pid, m=m, e=e.astype(np.float64),
                            e_native=float(e_nat), strain=int(strain),
                            rr=P["rr"][:m].astype(np.float64))
        os.replace(tmp, _apath(pid))
        if verbose:
            ok = int(np.isfinite(e).sum())
            print(f"[{k+1:3d}/{len(todo)}] {pid} ok {ok}/{m} strain {strain} "
                  f"({time.time()-t0:.0f}s)", flush=True)


def amber_cached():
    if not os.path.isdir(ACACHE):
        return []
    return sorted(f[:-4] for f in os.listdir(ACACHE) if f.endswith(".npz"))


def load_amber(pdbid):
    z = np.load(_apath(pdbid), allow_pickle=True)
    return {"pdb": str(z["pdb"]), "m": int(z["m"]), "e": z["e"].astype(float),
            "e_native": float(z["e_native"]), "strain": int(z["strain"]),
            "rr": z["rr"].astype(float)}


# ============================================================ the native's own scores
# ORACLE DIAGNOSTIC.  The native is appended to the pool as a 501st candidate and scored
# by every channel through the identical code path, so `native_pct` -- the fraction of the
# pool that outscores the true structure -- is measured on a like-for-like basis.  Two
# pool-dependent signals (`typic`, `bagged`) are functions of the pool's own distribution,
# which is why the native has to enter as an extra row rather than be scored alone.
#
# This is the statistic S8-1 warned about: it is degenerate for an ORACLE objective (the
# native scores exactly 0 and is trivially first).  Every channel here is LEARNED or
# PHYSICAL and native-free, so it is meaningful for all of them.
NPATH = os.path.join(HERE, "integrate_native.npz")


def stage_native(limit=None, verbose=True):
    from s7 import poolsize
    import peptide_db as db
    from s8 import invfold as ivf
    tgs = chan_cached()
    if limit:
        tgs = tgs[:int(limit)]
    models = {v: ivf.load_model(v) for v in ("geom", "nbid", "geom200k")}
    rows, keep, names_out = {}, [], None
    for k, pid in enumerate(tgs):
        t0 = time.time()
        u = gen.load_univ(pid)
        pred = poolsize.load_pred(pid)
        tp = db.by_pdb(pid)
        n = int(u["n"])
        idx = np.asarray(u["order"], int)[:K]
        W = np.vstack([np.asarray(u["W"], float)[idx],
                       np.asarray(u["nat_ca"], float)[None]])
        PHI = np.vstack([np.asarray(u["PHI"], float)[idx],
                         np.asarray(tp.phi, float)[None]])
        PSI = np.vstack([np.asarray(u["PSI"], float)[idx],
                         np.asarray(tp.psi, float)[None]])
        PHI[:, 0] = geo.DEFAULT_PHI
        PSI[:, -1] = geo.DEFAULT_PSI
        BB = geo.build_backbone_batch(PHI, PSI)
        ctx = {"n": n, "seq": str(u["seq"]), "fold": int(u["fold"]), "W": W,
               "PHI": PHI, "PSI": PSI, "D": gen.pair_D(W, n), "BB": BB,
               "pred": pred, "B": len(W)}
        S, names = ib.signal_block(ctx)
        vec = [S[-1, names.index(nm)] for nm in FAST]
        nm_out = list(FAST)
        for v in ("geom", "nbid", "geom200k"):
            sc = ivf.SelfConsistency(str(u["seq"]), models[v], nbid=ivf.VARIANTS[v][0])
            vec.append(float(sc.score(BB, PHI, PSI)[-1]))
            nm_out.append("iv_" + v)
            if v in ("geom", "nbid"):
                vec.append(float(sc.score_centred(BB, PHI, PSI)[-1]))
                nm_out.append("iv_%s_ctr" % v)
        rows[pid] = np.array(vec, float)
        keep.append(pid)
        names_out = nm_out
        if verbose and (k + 1) % 20 == 0:
            print(f"[{k+1:3d}/{len(tgs)}] {pid} ({time.time()-t0:.1f}s)", flush=True)
    np.savez_compressed(NPATH, pdbs=np.array(keep),
                        V=np.stack([rows[p] for p in keep]),
                        names=np.array(names_out))
    if verbose:
        print(f"native scores for {len(keep)} targets -> {NPATH}", flush=True)


def load_native():
    if not os.path.exists(NPATH):
        return None
    z = np.load(NPATH, allow_pickle=True)
    return {str(p): v for p, v in zip(z["pdbs"], z["V"])}, [str(x) for x in z["names"]]


# ============================================================ 1. ERROR CORRELATION
# The precondition for everything else in this task.  Two channels can only help each
# other if what they get WRONG differs.  Three statistics, and the difference between them
# is the whole subtlety:
#
#   rho_ab      raw rank correlation of the two scores.  Confounded: two channels that
#               both have skill are correlated for a GOOD reason.
#   err_ab      corr(rank_a - rank_true, rank_b - rank_true), literally what the brief
#               asked for.  Its null is NOT zero -- for two independent zero-skill
#               channels it is +0.5, because both errors contain -rank_true.  Reported
#               with that null stated.
#   part_ab     rho_ab partialled for the true RMSD rank.  This is the statistic that is
#               ZERO under conditional independence, and it is the one to read.
#
# All three are algebraically linked (same three Spearman correlations), so this is one
# measurement presented three ways, not three measurements.  The payoff is the last
# column: for z-scored channels the optimal linear fusion has skill
#     rho_comb = sqrt((ra^2 + rb^2 - 2*ra*rb*rho_ab) / (1 - rho_ab^2))
# so the matrix PREDICTS the fusion gain before any combiner is fitted.
def _pair_stats(sa, sb, rr):
    ra, rb, rt = _z(_rk(sa)), _z(_rk(sb)), _z(_rk(rr))
    r_ab, r_at, r_bt = float(np.mean(ra * rb)), float(np.mean(ra * rt)), float(np.mean(rb * rt))
    den = np.sqrt(max((1 - r_at ** 2) * (1 - r_bt ** 2), 1e-12))
    part = (r_ab - r_at * r_bt) / den
    ea, eb = ra - rt, rb - rt
    err = float(np.mean(ea * eb) / max(np.std(ea) * np.std(eb), 1e-12))
    return {"rho_ab": r_ab, "rho_at": r_at, "rho_bt": r_bt, "part_ab": float(part),
            "err_ab": err}


def rho_comb_pred(ra, rb, rab):
    """Skill of the OPTIMAL linear fusion of two z-scored channels.

    Standard two-predictor multiple-correlation identity.  It is deliberately computed
    from the AGGREGATED correlations, never per target and averaged: a per-target |rho|
    is dominated by sampling noise on a 130-member band, and averaging its magnitude
    would manufacture a fusion gain out of nothing.  (An earlier pass of this function did
    exactly that and predicted +0.41 in-band from two channels whose aggregate skills are
    +0.126 and +0.058; the honest figure is +0.130.)
    """
    if ra is None or rb is None or rab is None:
        return None
    num = ra ** 2 + rb ** 2 - 2 * ra * rb * rab
    return float(min(np.sqrt(max(num, 0.0) / max(1 - rab ** 2, 1e-9)), 1.0))


def _agg(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], float)
    if not len(v):
        return None
    return {"mean": float(v.mean()), "median": float(np.median(v)),
            "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "se": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0,
            "n": int(len(v))}


def stage_corr(out=None, verbose=True):
    """THE ERROR CORRELATION MATRIX.  Cheap, and it decides the rest of the task."""
    tgs = [p for p in chan_cached() if os.path.exists(_apath(p))]
    all_t = chan_cached()
    # report the m actually in the cache, not the module default: the stage takes `m` as
    # an argument and the two drifted apart once (the header said 40 while the data was 32)
    m_cached = min([load_amber(p)["m"] for p in tgs], default=M_AMBER)
    res = {"n_targets_full": len(all_t), "n_targets_amber": len(tgs), "K": K,
           "M_AMBER": int(m_cached), "BAND": BAND}

    for scope, names, tlist in (("full", list(FAST) + list(IV), all_t),
                                ("amber", list(FAST) + list(IV) + ["amber"], tgs)):
        if not tlist:
            continue
        cols = [PRIMARY_COL[p] for p in PRIMARY if PRIMARY_COL[p] in names]
        prim = [p for p in PRIMARY if PRIMARY_COL[p] in names]
        acc = {}          # (region, a, b) -> list
        skill = {}        # (region, a) -> list of (rho, rho_partial_rg)
        for pid in tlist:
            c = load_chan(pid, with_amber=(scope == "amber"))
            idx = c["names"]
            rr, rg = c["rr"], c["rg"]
            lim = c.get("amber_m", len(rr)) if scope == "amber" else len(rr)
            regions = {"pool": np.arange(lim)}
            band = np.where(rr[:lim] <= rr[:lim].min() + BAND)[0]
            if len(band) >= MIN_BAND:
                regions["band"] = band
            for reg, sub in regions.items():
                # Amber fails on a handful of strained candidates.  Drop the failed
                # CANDIDATES, not the whole target -- dropping the target would silently
                # restrict every Amber pair to the minority of targets on which no
                # candidate blew up, which is a selection effect on strain.
                ok = {p: c["S"][sub, idx.index(col)] for p, col in zip(prim, cols)}
                for p, v in ok.items():
                    m = np.isfinite(v)
                    if m.sum() < MIN_BAND:
                        continue
                    skill.setdefault((reg, p), []).append(
                        (_spearman(v[m], rr[sub][m]),
                         partial_rg(v[m], rr[sub][m], rg[sub][m])))
                ps = list(ok)
                for i, a in enumerate(ps):
                    for b in ps[i:]:
                        m = np.isfinite(ok[a]) & np.isfinite(ok[b])
                        if m.sum() < MIN_BAND:
                            continue
                        st = _pair_stats(ok[a][m], ok[b][m], rr[sub][m])
                        acc.setdefault((reg, a, b), []).append(st)
        blk = {"channels": prim, "skill": {}, "pairs": {}}
        for (reg, p), v in skill.items():
            blk["skill"].setdefault(reg, {})[p] = {
                "rho": _agg([x[0] for x in v]),
                "rho_partial_rg": _agg([x[1] for x in v])}
        for (reg, a, b), v in acc.items():
            d = {k: _agg([x[k] for x in v]) for k in v[0]}
            sk = blk["skill"][reg]
            d["rho_comb_pred"] = {"mean": rho_comb_pred(
                sk[a]["rho"]["mean"], sk[b]["rho"]["mean"], d["rho_ab"]["mean"])}
            blk["pairs"].setdefault(reg, {})[f"{a}|{b}"] = d
        res[scope] = blk

    path = out or os.path.join(HERE, "integrate_corr.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        _print_corr(res)
    return res


def _fmt(a, d=3):
    return "  n/a " if a is None else f"{a['mean']:+.{d}f}"


def _print_corr(res):
    for scope in ("full", "amber"):
        if scope not in res:
            continue
        blk = res[scope]
        n = res["n_targets_amber"] if scope == "amber" else res["n_targets_full"]
        print(f"\n===== {scope.upper()} ({n} targets, "
              f"{'first %d BLOSUM candidates' % res['M_AMBER'] if scope=='amber' else 'K=500'}) =====")
        for reg in ("pool", "band"):
            if reg not in blk["pairs"]:
                continue
            print(f"\n--- {reg} ---")
            ch = blk["channels"]
            print("skill (Spearman vs true CA-RMSD; second row = rg partialled out)")
            print("       " + "".join(f"{c:>12s}" for c in ch))
            print("rho    " + "".join(f"{_fmt(blk['skill'][reg][c]['rho']):>12s}" for c in ch))
            print("prg    " + "".join(
                f"{_fmt(blk['skill'][reg][c]['rho_partial_rg']):>12s}" for c in ch))
            for key, lbl in (("rho_ab", "score rank corr"),
                             ("part_ab", "ERROR corr (truth-partialled; 0 = independent)"),
                             ("err_ab", "raw error corr (null = +0.5 at zero skill)")):
                print(f"\n{lbl}")
                print("       " + "".join(f"{c:>12s}" for c in ch))
                for a in ch:
                    row = f"{a:<7s}"
                    for b in ch:
                        k = f"{a}|{b}" if f"{a}|{b}" in blk["pairs"][reg] else f"{b}|{a}"
                        row += f"{_fmt(blk['pairs'][reg].get(k, {}).get(key)):>12s}"
                    print(row)
            print("\npredicted skill of the OPTIMAL 2-channel linear fusion")
            print("       " + "".join(f"{c:>12s}" for c in ch))
            for a in ch:
                row = f"{a:<7s}"
                for b in ch:
                    k = f"{a}|{b}" if f"{a}|{b}" in blk["pairs"][reg] else f"{b}|{a}"
                    row += f"{_fmt(blk['pairs'][reg].get(k, {}).get('rho_comb_pred')):>12s}"
                print(row)


# ============================================================ 2. MULTI-OBJECTIVE SCORE
def _paired(a, b, boot=4000, seed=0):
    """Paired difference a-b with a bootstrap 95% CI and win/loss counts."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a - b
    rng = np.random.default_rng(seed)
    bs = d[rng.integers(0, len(d), (boot, len(d)))].mean(1)
    return {"mean_diff": float(d.mean()), "se": float(d.std(ddof=1) / np.sqrt(len(d))),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "wins": int((d < -1e-9).sum()), "losses": int((d > 1e-9).sum()),
            "ties": int((np.abs(d) <= 1e-9).sum()), "n": int(len(d))}


def _arm_summary(sel, nat_pct, nat_first, rho_g, rho_b, rho_bp, base_sel):
    sel = np.asarray(sel, float)
    d = {"sel": float(sel.mean()), "sel_median": float(np.median(sel)),
         "sel_sd": float(sel.std(ddof=1)), "sel_se": float(sel.std(ddof=1) / np.sqrt(len(sel))),
         "best": float(sel.min()), "worst": float(sel.max()),
         "frac_lt2": float((sel < 2.0).mean()), "frac_lt15": float((sel < 1.5).mean()),
         "n": int(len(sel))}
    if base_sel is not None:
        d["vs_dist"] = _paired(sel, base_sel)
    for k, v in (("native_pct", nat_pct), ("rho_global", rho_g),
                 ("rho_band", rho_b), ("rho_band_prg", rho_bp)):
        a = _agg(v)
        if a:
            d[k] = a
    if nat_first:
        d["native_argmin"] = int(sum(nat_first))
    return d


class Fuser:
    """A multi-objective score over `cols`, as a function of standardised RANKS.

    Ranks, not raw values, for a reason this project has already been bitten by: the
    Legacy energies and the Amber energy have heavy tails whose scale varies by orders of
    magnitude between targets, and a z-score of the raw value would let one target's
    steric blow-up dominate a pooled fit.  Ranking within a target also removes every
    per-target scale, so nothing about the target's difficulty can leak into the weights.
    """

    def __init__(self, cols, w, mode="linear"):
        self.cols, self.w, self.mode = list(cols), np.asarray(w, float), mode

    def __call__(self, M):
        """`M` is (B, n_cols) of raw channel scores -> (B,) fused score, lower better."""
        Z = np.stack([zrank(M[:, i]) for i in range(M.shape[1])], 1)
        if self.mode == "borda":
            return (Z * np.sign(self.w)[None, :]).sum(1)
        return Z @ self.w


def _cols_of(c, cols):
    return np.stack([c["S"][:, c["names"].index(x)] for x in cols], 1)


def _eval_fuser(caches, natv, nat_names, cols, wfn, lim=None):
    """Selected RMSD, native percentile and rhos for one fuser over every target."""
    sel, npct, nfirst, rg_, rb_, rbp_ = [], [], [], [], [], []
    used = []
    for ci, c in enumerate(caches):
        L = len(c["rr"]) if lim is None else min(lim, len(c["rr"]))
        M = _cols_of(c, cols)[:L]
        rr, rg = c["rr"][:L], c["rg"][:L]
        if not np.isfinite(M).all():
            continue
        w = wfn(c)
        if w is None:
            continue
        nv = None
        if natv is not None and c["pdb"] in natv:
            nv = np.array([natv[c["pdb"]][nat_names.index(x)] for x in cols], float)
        if nv is not None and np.isfinite(nv).all():
            Mx = np.vstack([M, nv[None]])
            s_all = w(Mx)
            s, s_nat = s_all[:-1], s_all[-1]
            npct.append(100.0 * float((s < s_nat).mean()))
            nfirst.append(bool(s_nat < s.min()))
        else:
            s = w(M)
        sel.append(sel_of(s, rr))
        used.append(ci)
        rg_.append(_spearman(s, rr))
        band = rr <= rr.min() + BAND
        if band.sum() >= MIN_BAND:
            rb_.append(_spearman(s[band], rr[band]))
            rbp_.append(partial_rg(s[band], rr[band], rg[band]))
    return sel, npct, nfirst, rg_, rb_, rbp_, np.array(used, int)


def _fit_lfo(caches, cols, folds=5, lam=1e-2, fit_band=False):
    """Ridge weights per held-out fold.  Training folds only; no target sees its own fit.

    Target variable is the candidate's true CA-RMSD CENTRED WITHIN ITS TARGET, so the fit
    learns to rank inside a pool rather than to predict which targets are hard -- the
    latter is what a pooled regression would spend its capacity on and it is useless for
    selection.
    """
    W = {}
    for held in range(folds):
        X, Y = [], []
        for c in caches:
            if c["fold"] == held:
                continue
            M = _cols_of(c, cols)
            rr = c["rr"]
            if not np.isfinite(M).all():
                continue
            m = (rr <= rr.min() + BAND) if fit_band else np.ones(len(rr), bool)
            if m.sum() < MIN_BAND:
                continue
            Z = np.stack([zrank(M[m][:, i]) for i in range(M.shape[1])], 1)
            X.append(Z)
            Y.append(rr[m] - rr[m].mean())
        if not X:
            W[held] = None
            continue
        A = np.vstack(X)
        y = np.concatenate(Y)
        G = A.T @ A + lam * len(A) * np.eye(A.shape[1])
        W[held] = np.linalg.solve(G, A.T @ y)
    return W


def stage_fuse(out=None, verbose=True, amber=1):
    natv, nat_names = load_native() or (None, None)
    tgs = chan_cached()
    caches = [load_chan(p) for p in tgs]
    base = np.array([sel_of(_cols_of(c, ["score"])[:, 0], c["rr"]) for c in caches])
    res = {"n_targets": len(caches), "K": K, "BAND": BAND,
           "baseline_sel": float(base.mean()), "arms": {}}

    D, LG, IV1 = "score", "lg_all", "iv_nbid"

    def fixed(cols, w):
        f = Fuser(cols, w)
        return lambda c: f

    arms = {"dist": ([D], [1.0])}
    for w in (0.25, 0.5, 1.0):
        arms[f"dist+legacy_w{w}"] = ([D, LG], [1.0, w])
        arms[f"dist+invfold_w{w}"] = ([D, IV1], [1.0, w])
        arms[f"dist+legacy+invfold_w{w}"] = ([D, LG, IV1], [1.0, w, w])
    arms["borda3"] = ([D, LG, IV1], [1.0, 1.0, 1.0])
    arms["legacy"] = ([LG], [1.0])
    arms["invfold"] = ([IV1], [1.0])

    for name, (cols, w) in arms.items():
        mode = "borda" if name.startswith("borda") else "linear"
        f = Fuser(cols, w, mode)
        r = _eval_fuser(caches, natv, nat_names, cols, lambda c, f=f: f)
        res["arms"][name] = _arm_summary(r[0], r[1], r[2], r[3], r[4], r[5], base[r[6]])

    # ---- learned combiners, leave-fold-out ----
    subsets = {
        "lfo_primary3": [D, LG, IV1],
        "lfo_dist_family": ["score", "l1exp", "zdev", "bagged", "typic"],
        "lfo_legacy11": [n for n in FAST if n.startswith("lg_") and n not in
                         ("lg_all", "lg_hb", "lg_burial", "lg_packing")],
        "lfo_iv": list(IV),
        "lfo_all": [n for n in FAST if n != "rg"] + list(IV),
        "lfo_all_plus_rg": list(FAST) + list(IV),
    }
    for tag, cols in subsets.items():
        for fit_band in (False, True):
            Wf = _fit_lfo(caches, cols, fit_band=fit_band)
            nm = tag + ("_bandfit" if fit_band else "")

            def wfn(c, Wf=Wf, cols=cols):
                w = Wf.get(c["fold"])
                return None if w is None else Fuser(cols, w)
            r = _eval_fuser(caches, natv, nat_names, cols, wfn)
            res["arms"][nm] = _arm_summary(r[0], r[1], r[2], r[3], r[4], r[5],
                                           base[r[6]])
            res["arms"][nm]["weights_fold0"] = (
                None if Wf.get(0) is None else [float(x) for x in Wf[0]])
            res["arms"][nm]["cols"] = cols

    path = out or os.path.join(HERE, "integrate_fuse.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        _print_fuse(res)
    return res


def _print_fuse(res, key="arms"):
    print(f"\nn={res['n_targets']}  baseline (shipped distogram) = "
          f"{res['baseline_sel']:.3f} A\n")
    hdr = ("%-28s%9s%8s%8s%9s%9s%8s%8s%9s%9s%9s" %
           ("arm", "sel", "med", "SD", "d vs D", "CI_lo", "CI_hi", "W/L",
            "nat_pct", "rho_g", "rho_bprg"))
    print(hdr)
    rows = sorted(res[key].items(), key=lambda kv: kv[1]["sel"])
    for nm, d in rows:
        v = d.get("vs_dist", {})
        ci = v.get("ci95", [float("nan")] * 2)
        print("%-28s%9.3f%8.3f%8.3f%9.3f%9.3f%8.3f%8s%9.1f%9.3f%9.3f" % (
            nm, d["sel"], d["sel_median"], d["sel_sd"], v.get("mean_diff", float("nan")),
            ci[0], ci[1], "%d/%d" % (v.get("wins", 0), v.get("losses", 0)),
            d.get("native_pct", {}).get("mean", float("nan")),
            d.get("rho_global", {}).get("mean", float("nan")),
            d.get("rho_band_prg", {}).get("mean", float("nan"))))


# ============================================================ 3. FILTER + CONSENSUS
# S8-8 (`s8/inband.py`, another agent, same instrument) established the composition this
# section builds on: the shipped score is a good coarse FILTER and a bad fine RANKER, and
# consensus is the reverse.  `filter to the score's top 75, return the consensus medoid`
# is 3.282 A, -0.172 [-0.316, -0.027].
#
# The question this module asks is the one that composition leaves open and that only a
# multi-channel study can answer: **is the multi-objective score a better FILTER?**  That
# is a different question from whether it is a better ranker -- a filter only has to keep
# the near-native band inside the top q, which is exactly the property S8-5 showed BLOSUM
# has and every better-correlated key lacks.
DPATH = os.path.join(HERE, "integrate_dmat")


def _dpath(pdbid):
    return os.path.join(DPATH, f"{pdbid}.npz")


def stage_dmat(limit=None, verbose=True):
    """Cache the 500x500 pairwise CA-RMSD matrix of every target's pool.

    Stored float16: the consensus medoid is an argmin over row means of a matrix whose
    entries are 0.3-8 A, so 3 significant figures is far more resolution than the decision
    needs, and it keeps the whole instrument under 40 MB on a contended box.
    """
    from s5.lib import kabsch_rmsd_batch
    os.makedirs(DPATH, exist_ok=True)
    tgs = targets()
    if limit:
        tgs = tgs[:int(limit)]
    todo = [p for p in tgs if not os.path.exists(_dpath(p))]
    if verbose:
        print(f"dmat: {len(todo)}/{len(tgs)} targets to do", flush=True)
    for k, pid in enumerate(todo):
        t0 = time.time()
        X = pool(gen.load_univ(pid))["W"]
        B = len(X)
        D = np.zeros((B, B), np.float32)
        for i in range(B):
            D[i, i + 1:] = kabsch_rmsd_batch(X[i + 1:], X[i])
        D = (D + D.T).astype(np.float16)
        tmp = _dpath(pid) + ".tmp.npz"
        np.savez_compressed(tmp, D=D)
        os.replace(tmp, _dpath(pid))
        if verbose and (k + 1) % 20 == 0:
            print(f"[{k+1:3d}/{len(todo)}] {pid} ({time.time()-t0:.1f}s)", flush=True)


_DMEMO = {}


def load_dmat(pdbid):
    if pdbid not in _DMEMO:
        if len(_DMEMO) > 8:
            _DMEMO.clear()
        _DMEMO[pdbid] = np.load(_dpath(pdbid))["D"].astype(np.float32)
    return _DMEMO[pdbid]


def consensus_medoid(D, idx, w=None):
    """The member of `idx` closest to all the others -- the ensemble's representative.

    `w` optionally weights the ensemble (a Boltzmann weight on the filter score).  With
    `w=None` this is exactly S8-8's `medoid`, so the two studies' numbers are comparable.
    """
    sub = D[np.ix_(idx, idx)].astype(float)
    v = sub.mean(1) if w is None else sub @ (w / max(w.sum(), 1e-12))
    return int(idx[int(np.argmin(v))])


def filter_consensus(s, D, q, rr, beta=0.0):
    """Filter to the top `q` of `s`, return the consensus medoid's true RMSD."""
    q = min(int(q), len(s))
    idx = np.argsort(np.asarray(s, float), kind="stable")[:q]
    w = None
    if beta > 0:
        w = np.exp(-beta * zrank(np.asarray(s, float)[idx]))
    return float(rr[consensus_medoid(D, idx, w)])


def _filters(caches):
    """Every FILTER score this study puts on trial, as {name: per-target score vector}.

    Deliberately small and mechanism-diverse.  The four fixed fusions are the mandated
    channels; the two `lfo_*` arms are the learned combiners fitted leave-fold-out; the
    single-channel arms are the controls that say whether fusing bought anything.
    """
    D, LG, IVN, TY = "score", "lg_all", "iv_nbid", "typic"
    fixed = {
        "score": [(D, 1.0)],
        "typic": [(TY, 1.0)],
        "legacy": [(LG, 1.0)],
        "invfold": [(IVN, 1.0)],
        "D+L.25": [(D, 1.0), (LG, 0.25)],
        "D+I.25": [(D, 1.0), (IVN, 0.25)],
        "D+L+I.25": [(D, 1.0), (LG, 0.25), (IVN, 0.25)],
        "D+L+I.5": [(D, 1.0), (LG, 0.5), (IVN, 0.5)],
        "D+T.25": [(D, 1.0), (TY, 0.25)],
        "D+T.5": [(D, 1.0), (TY, 0.5)],
        "D+T+I.25": [(D, 1.0), (TY, 0.25), (IVN, 0.25)],
        "borda_DLI": [(D, 1.0), (LG, 1.0), (IVN, 1.0)],
    }
    out = {}
    for nm, spec in fixed.items():
        out[nm] = {c["pdb"]: sum(w * zrank(c["S"][:, c["names"].index(n)])
                                 for n, w in spec) for c in caches}
    for tag, cols in (("lfo_DLI", ["score", "lg_all", "iv_nbid"]),
                      ("lfo_DLIT", ["score", "lg_all", "iv_nbid", "typic"]),
                      ("lfo_all", [n for n in FAST if n != "rg"] + list(IV))):
        Wf = _fit_lfo(caches, cols)
        d = {}
        for c in caches:
            w = Wf.get(c["fold"])
            if w is not None:
                d[c["pdb"]] = Fuser(cols, w)(_cols_of(c, cols))
        out[tag] = d
    return out


QGRID = (6, 10, 15, 25, 40, 50, 75, 100, 150, 200, 300)


def stage_filter(out=None, verbose=True, beta=0.0):
    """The filter x q grid, with leave-fold-out arm selection on top of it."""
    tgs = chan_cached()
    caches = [load_chan(p) for p in tgs]
    folds = {c["pdb"]: c["fold"] for c in caches}
    rrs = {c["pdb"]: c["rr"] for c in caches}
    base = {c["pdb"]: sel_of(c["S"][:, c["names"].index("score")], c["rr"]) for c in caches}
    F = _filters(caches)
    cells = {}
    for fname, per in F.items():
        for q in QGRID:
            sel = {}
            for pid, s in per.items():
                sel[pid] = filter_consensus(s, load_dmat(pid), q, rrs[pid], beta)
            cells[f"{fname}@{q}"] = sel
    # the S8-8 headline, reproduced here on identical candidates
    order = sorted(cells, key=lambda k: np.mean(list(cells[k].values())))
    res = {"n_targets": len(caches), "beta": beta,
           "baseline_sel": float(np.mean([base[p] for p in tgs])),
           "cells": {}, "QGRID": list(QGRID)}
    b = np.array([base[p] for p in tgs])
    for k in cells:
        v = np.array([cells[k][p] for p in tgs])
        res["cells"][k] = _arm_summary(v, [], [], [], [], [], b)

    # ---- leave-fold-out arm selection: the honest headline ----
    # Picking the best cell of a 120-cell grid on the same 126 targets it is reported on is
    # how a null becomes a claim.  Instead: for each held-out fold, choose the arm with the
    # lowest mean on the OTHER four folds, and use it only on the held-out one.
    lfo, picks = {}, {}
    for held in sorted(set(folds.values())):
        tr = [p for p in tgs if folds[p] != held]
        te = [p for p in tgs if folds[p] == held]
        best = min(cells, key=lambda k: np.mean([cells[k][p] for p in tr]))
        picks[held] = best
        for p in te:
            lfo[p] = cells[best][p]
    v = np.array([lfo[p] for p in tgs])
    res["consensus_LFO"] = _arm_summary(v, [], [], [], [], [], b)
    res["consensus_LFO"]["picks"] = picks
    res["top10"] = order[:10]
    path = out or os.path.join(HERE, "integrate_filter.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nn={res['n_targets']} baseline {res['baseline_sel']:.3f} beta={beta}\n")
        print("%-22s%9s%9s%9s%8s%8s%8s%8s" %
              ("cell", "sel", "d", "CI_lo", "CI_hi", "W/L", "<2A", "<1.5A"))
        for k in order[:25]:
            d = res["cells"][k]
            v_ = d["vs_dist"]
            print("%-22s%9.3f%9.3f%9.3f%8.3f%8s%8.3f%8.3f" % (
                k, d["sel"], v_["mean_diff"], v_["ci95"][0], v_["ci95"][1],
                "%d/%d" % (v_["wins"], v_["losses"]), d["frac_lt2"], d["frac_lt15"]))
        d = res["consensus_LFO"]
        print("\nLEAVE-FOLD-OUT arm selection: %.3f  d %+.3f %s W/L %d/%d" % (
            d["sel"], d["vs_dist"]["mean_diff"], d["vs_dist"]["ci95"],
            d["vs_dist"]["wins"], d["vs_dist"]["losses"]))
        print("picks per fold:", picks)
    return res


# ============================================================ 4. CVaR-VQE
# WHAT THE QUANTUM STATE REPRESENTS, AND WHY CVaR
# ----------------------------------------------
# S8-2 measured that a target's near-native band holds ~2-3 POPULATED structural
# hypotheses.  S8-8 and this study both measured that the right readout from a candidate
# set is not its argmin but its CONSENSUS -- the member closest to the rest of the
# distribution.  Those two facts fix the formulation:
#
#   * the `n`-qubit state's computational basis indexes the top `2**n` candidates of the
#     multi-objective filter, so a basis state IS a structural hypothesis and the state is
#     a distribution over hypotheses -- not a continuous geometric coordinate.  Sprint 7
#     measured that optimising a misspecified continuous geometric objective HARDER makes
#     structures worse (1KVG: objective 1.564 -> 1.495, CA-RMSD 2.276 -> 3.107), so the
#     continuous formulation is ruled out by this project's own data, not by taste.
#   * the objective is the multi-objective score, and CVaR -- not the mean -- is what is
#     optimised, because the quantity we want is the distribution's low-energy TAIL: alpha
#     is exactly the knob that says how much of the candidate set the consensus is taken
#     over.  alpha -> 0 recovers the shipped argmin selector (3.454 A); alpha -> 1 recovers
#     consensus over the whole filtered set.  The measured optimum is strictly inside, and
#     that is the sense in which CVaR is load-bearing here rather than decorative.
#   * the READOUT is the consensus medoid weighted by the optimised p_theta.  CVaR is
#     about a distribution's tail and consensus is about a distribution's centre; the
#     composition is "concentrate the state on the good tail, then return that tail's
#     centre", and both halves are stated because they are different operations.
#
# The circuit is simulated EXACTLY (2**n amplitudes, n <= 10), so p_theta and every
# gradient below are analytic, not sampled -- which is what makes the gradient audit in
# `stage_cvarcheck` a real verification rather than a noise comparison.
class Circuit:
    """`layers` x (RY on every wire, CNOT ring), exact statevector, real amplitudes.

    Real-amplitude RY+CNOT is the standard hardware-efficient ansatz and is expressive
    enough for a multimodal distribution over basis states at depth >= 2, which is the
    property S8-2's 2-3 populated hypotheses actually require.
    """

    def __init__(self, n, layers=3, ring=True):
        self.n, self.layers, self.ring = int(n), int(layers), bool(ring)
        self.dim = 1 << self.n

    def n_params(self):
        return self.n * self.layers

    def _ry(self, psi, q, t):
        n, c, s = self.n, np.cos(t / 2.0), np.sin(t / 2.0)
        psi = psi.reshape(1 << q, 2, 1 << (n - q - 1))
        a, b = psi[:, 0, :].copy(), psi[:, 1, :].copy()
        psi[:, 0, :] = c * a - s * b
        psi[:, 1, :] = s * a + c * b
        return psi.reshape(-1)

    def _cx(self, psi, ctrl, tgt):
        n = self.n
        idx = np.arange(self.dim)
        cb = (idx >> (n - 1 - ctrl)) & 1
        j = np.where(cb == 1, idx ^ (1 << (n - 1 - tgt)), idx)
        return psi[j]

    def state(self, theta):
        theta = np.asarray(theta, float).reshape(self.layers, self.n)
        psi = np.zeros(self.dim)
        psi[0] = 1.0
        for L in range(self.layers):
            for q in range(self.n):
                psi = self._ry(psi, q, theta[L, q])
            for q in range(self.n - 1):
                psi = self._cx(psi, q, q + 1)
            if self.ring and self.n > 2:
                psi = self._cx(psi, self.n - 1, 0)
        return psi

    def probs(self, theta):
        p = self.state(theta) ** 2
        return p / p.sum()


def cvar_exact(E, p, alpha):
    """CVaR of the LOWER `alpha` tail of the energy distribution, from exact `p`.

    CVaR_alpha = (1/alpha) [ sum_{E(x)<q} p(x) E(x) + (alpha - P(E<q)) q ], with `q` the
    alpha-quantile.  Returns (value, q, dCVaR/dp).  The derivative uses the envelope
    theorem -- q's own dependence on p cancels -- so `dp` is exactly (E - q)/alpha on the
    strict tail and 0 elsewhere, with the boundary state carrying its partial mass.
    """
    E, p = np.asarray(E, float), np.asarray(p, float)
    o = np.argsort(E, kind="stable")
    c = np.cumsum(p[o])
    k = int(np.searchsorted(c, alpha, side="left"))
    k = min(k, len(E) - 1)
    q = E[o[k]]
    below = o[:k]
    mass = c[k - 1] if k > 0 else 0.0
    val = (float((p[below] * E[below]).sum()) + (alpha - mass) * q) / alpha
    dp = np.zeros_like(p)
    dp[below] = (E[below] - q) / alpha
    return float(val), float(q), dp


def grad_cvar_paramshift(circ, theta, E, alpha):
    """EXACT gradient by the parameter-shift rule applied to each basis probability.

    p(x) = <psi|Pi_x|psi> is the expectation of a projector, so every p(x) obeys the exact
    two-term shift rule for an RY generator.  Chaining it with dCVaR/dp gives the exact
    analytic gradient with no sampling and no finite-difference step.  This is the
    reference the deployable estimator is audited against.
    """
    theta = np.asarray(theta, float)
    p = circ.probs(theta)
    _, _, dp = cvar_exact(E, p, alpha)
    g = np.zeros(len(theta))
    for k in range(len(theta)):
        tp, tm = theta.copy(), theta.copy()
        tp[k] += np.pi / 2
        tm[k] -= np.pi / 2
        g[k] = float(dp @ (circ.probs(tp) - circ.probs(tm)) / 2.0)
    return g


def grad_cvar_fd(circ, theta, E, alpha, h=1e-5):
    """Central finite differences on the exact CVaR value. Independent of the shift rule."""
    theta = np.asarray(theta, float)
    g = np.zeros(len(theta))
    for k in range(len(theta)):
        tp, tm = theta.copy(), theta.copy()
        tp[k] += h
        tm[k] -= h
        g[k] = (cvar_exact(E, circ.probs(tp), alpha)[0]
                - cvar_exact(E, circ.probs(tm), alpha)[0]) / (2 * h)
    return g


def grad_cvar_score(circ, theta, E, alpha, shots=8192, rng=None, baseline="global"):
    """The SAMPLED score-function estimator -- what a real device would compute.

    `baseline` is the whole point of this function.  A baseline b(x) leaves the estimator
    unbiased only if it is CONSTANT in x, because the correction term is
    b * E_p[grad log p] = 0.  The recorded defect `cvar-gradient-baseline-defect` is that
    the shipped `qansatz.cvar_gradient` subtracts the TAIL MEAN from the tail entries only
    and leaves the rest at zero -- i.e. b(x) = mean * 1[x in tail], a FUNCTION of x, which
    biases the estimator.  `baseline="tail"` reproduces that defect exactly so the audit
    measures it rather than asserting it.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    theta = np.asarray(theta, float)
    p = circ.probs(theta)
    x = rng.choice(len(p), size=int(shots), p=p)
    en = E[x]
    val, q, _ = cvar_exact(E, p, alpha)
    w = np.where(en < q, (en - q) / alpha, 0.0)
    if baseline == "global":
        w = w - w.mean()                    # constant in x -> unbiased
    elif baseline == "tail":
        m = en < q                          # the recorded DEFECT, reproduced verbatim
        if m.any():
            w[m] -= w[m].mean()
    # d log p(x) / d theta_k by the same shift rule, evaluated once per basis state
    G = np.zeros((len(p), len(theta)))
    for k in range(len(theta)):
        tp, tm = theta.copy(), theta.copy()
        tp[k] += np.pi / 2
        tm[k] -= np.pi / 2
        G[:, k] = (circ.probs(tp) - circ.probs(tm)) / 2.0 / np.maximum(p, 1e-15)
    return (w[:, None] * G[x]).mean(0), val


def _cos(a, b):
    return float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-15))


def stage_cvarcheck(n=7, layers=3, trials=12, out=None, verbose=True):
    """Verify the CVaR gradient against two independent references, and price the defect.

    Three estimators of the same quantity, on the SAME theta and the SAME energies:
    parameter-shift (exact), central finite differences on the exact CVaR (exact,
    independent machinery), and the sampled score-function estimator with each of the two
    baselines.  A correct implementation makes the first two agree to ~1e-6 and the
    globally-baselined sampler agree with them up to sampling noise.
    """
    caches = [load_chan(p) for p in chan_cached()[:trials]]
    circ = Circuit(n, layers)
    rng = np.random.default_rng(0)
    rows = []
    for c in caches:
        s = c["S"][:, c["names"].index("score")]
        o = np.argsort(s, kind="stable")[:circ.dim]
        E = zrank(s[o])
        for alpha in (0.1, 0.25, 0.5):
            th = rng.normal(0, 1.0, circ.n_params())
            gp = grad_cvar_paramshift(circ, th, E, alpha)
            gf = grad_cvar_fd(circ, th, E, alpha)
            gg, _ = grad_cvar_score(circ, th, E, alpha, rng=rng, baseline="global")
            gt, _ = grad_cvar_score(circ, th, E, alpha, rng=rng, baseline="tail")
            rows.append({
                "pdb": c["pdb"], "alpha": alpha,
                "cos_shift_fd": _cos(gp, gf),
                "rel_err_shift_fd": float(np.linalg.norm(gp - gf) /
                                          max(np.linalg.norm(gp), 1e-12)),
                "cos_score_global": _cos(gg, gp),
                "cos_score_tailbaseline": _cos(gt, gp),
                "relnorm_global": float(np.linalg.norm(gg) / max(np.linalg.norm(gp), 1e-12)),
                "relnorm_tail": float(np.linalg.norm(gt) / max(np.linalg.norm(gp), 1e-12))})
    res = {"n_qubits": n, "layers": layers, "n_params": circ.n_params(),
           "n_checks": len(rows), "rows": rows,
           "agg": {k: _agg([r[k] for r in rows]) for k in rows[0] if k != "pdb"}}
    path = out or os.path.join(HERE, "integrate_cvarcheck.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nCVaR gradient audit: n={n} qubits, {layers} layers, "
              f"{circ.n_params()} params, {len(rows)} checks")
        for k, v in res["agg"].items():
            if k == "alpha" or v is None:
                continue
            print(f"  {k:26s} mean {v['mean']:+.6f}  median {v['median']:+.6f}  "
                  f"sd {v['sd']:.6f}")
    return res


def free_energy(circ, theta, E, alpha, T):
    """F = CVaR_alpha(E; p_theta) - T * H(p_theta), and its exact gradient.

    WHY AN ENTROPY TERM, AND WHY IT IS NOT A FUDGE
    ----------------------------------------------
    Minimising CVaR alone is degenerate for this task: for ANY alpha the minimiser
    concentrates p on the lowest-energy basis states, so the readout collapses back to the
    argmin -- the shipped selector, 3.454 A.  That was measured before this term existed
    (every `vqe_a*` arm returned the argmin's structure, state entropy 0.01 bits at
    alpha=1), and it is a property of CVaR, not of the optimiser.

    The quantity the consensus readout actually needs is an ENSEMBLE: concentrated on good
    hypotheses, but still broad enough to have a centre.  That is a free energy -- energy
    minus temperature times entropy -- and it is the same object the project's
    ensemble-selection result rests on.  So `alpha` says which part of the energy
    distribution is scored and `T` says how broad the ensemble is, and both are measured.
    """
    p = circ.probs(theta)
    v, q, dp = cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-15))
    H = float(-(p * lp).sum())
    dH = -(lp + 1.0)
    d = dp - T * dH
    g = np.zeros(len(theta))
    for k in range(len(theta)):
        tp, tm = np.array(theta, float), np.array(theta, float)
        tp[k] += np.pi / 2
        tm[k] -= np.pi / 2
        g[k] = float(d @ (circ.probs(tp) - circ.probs(tm)) / 2.0)
    return float(v - T * H), g, p, v, H


def run_cvar_vqe(E, alpha, T=0.0, n=7, layers=3, iters=80, restarts=1, seed=0, lr=0.15):
    """Optimise the CVaR free energy of the hypothesis distribution; return p_theta.

    Adam on the exact parameter-shift gradient.  The gradient is the one audited in
    `stage_cvarcheck`; the sampled estimator agrees with it at cosine 0.994, so using the
    exact one here removes shot noise from a selection experiment without changing what is
    being optimised.
    """
    circ = Circuit(n, layers)
    rng = np.random.default_rng(seed)
    best, best_f = None, np.inf
    for r in range(restarts):
        th = rng.normal(0.0, 0.6, circ.n_params())
        m = np.zeros_like(th)
        v = np.zeros_like(th)
        for t in range(1, iters + 1):
            f, g, _, _, _ = free_energy(circ, th, E, alpha, T)
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g * g
            th = th - lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
        f, _, p, cv, H = free_energy(circ, th, E, alpha, T)
        if f < best_f:
            best_f, best = f, (p, cv, H)
    return best[0], float(best[1]), float(best[2]), circ


ALPHAS = (0.1, 0.25, 1.0)
TEMPS = (0.1, 0.3, 1.0)


def stage_vqe(n=7, layers=3, iters=60, out=None, verbose=True, filt="score", limit=None):
    """CVaR-VQE selection, with the classical ablations that price it honestly."""
    tgs = chan_cached()
    if limit:
        tgs = tgs[:int(limit)]
    caches = [load_chan(p) for p in tgs]
    per = _filters(caches)["score" if filt == "crit" else filt]
    dim = 1 << int(n)
    names = ([f"vqe_a{a}_T{T}" for a in ALPHAS for T in TEMPS]
             + [f"boltz_T{T}" for T in TEMPS]
             + [f"topfrac_{f}" for f in (0.05, 0.1, 0.25, 0.5)]
             + ["argmin", f"medoid{dim}"])
    arms = {k: [] for k in names}
    base = []
    ent = {k: [] for k in names if k.startswith("vqe")}
    for c in caches:
        pid = c["pdb"]
        s = np.asarray(per[pid], float)
        rr, D = c["rr"], load_dmat(pid)
        o = np.argsort(s, kind="stable")[:dim]
        if filt == "crit":
            # The brief's rule: establish the objective first, then let VQE optimise it.
            # The best objective this study measured is the CONSENSUS CRITERION plus half
            # an inverse-folding term, so that -- not the raw shipped score -- is the
            # energy the quantum state's basis is ordered by and the CVaR tail is taken
            # over.  Sprint 7 measured that optimising a MISspecified objective harder
            # makes structures worse, so which objective goes in is the whole question.
            E = zrank(zrank(load_dmat(pid)[np.ix_(o, o)].astype(float).mean(1))
                      + 0.5 * zrank(c["S"][o, c["names"].index("iv_nbid")]))
        else:
            E = zrank(s[o])
        base.append(sel_of(c["S"][:, c["names"].index("score")], rr))
        arms["argmin"].append(float(rr[o[int(np.argmin(E))]]))
        arms[f"medoid{dim}"].append(float(rr[consensus_medoid(D, o)]))
        for f in (0.05, 0.1, 0.25, 0.5):
            # CLASSICAL ABLATION 1: a hard top-fraction of the same filtered set.
            k = max(1, int(round(f * dim)))
            arms[f"topfrac_{f}"].append(float(rr[consensus_medoid(D, o[:k])]))
        for T in TEMPS:
            # CLASSICAL ABLATION 2: the exact Boltzmann ensemble the free energy would
            # reach with an UNCONSTRAINED distribution.  The VQE's p_theta is confined to
            # what a 3-layer RY/CNOT circuit can express, so this is the reference that
            # says whether the circuit's structure helps, hurts, or is neutral.
            w = np.exp(-E / T)
            arms[f"boltz_T{T}"].append(float(rr[consensus_medoid(D, o, w / w.sum())]))
        for a in ALPHAS:
            for T in TEMPS:
                p, _, H, _ = run_cvar_vqe(E, a, T, n=n, layers=layers, iters=iters)
                arms[f"vqe_a{a}_T{T}"].append(float(rr[consensus_medoid(D, o, p)]))
                ent[f"vqe_a{a}_T{T}"].append(H / np.log(2))
    b = np.array(base)
    res = {"n_targets": len(caches), "n_qubits": int(n), "layers": int(layers),
           "dim": dim, "filter": filt, "iters": int(iters), "alphas": list(ALPHAS),
           "temps": list(TEMPS), "baseline_sel": float(b.mean()), "arms": {},
           "pdbs": [c["pdb"] for c in caches],
           "folds": [c["fold"] for c in caches],
           "base_per_target": [float(x) for x in b],
           "per_target": {k: [float(x) for x in v] for k, v in arms.items()},
           "entropy_bits": {k: float(np.mean(v)) for k, v in ent.items()}}
    for k, v in arms.items():
        res["arms"][k] = _arm_summary(np.array(v), [], [], [], [], [], b)
    # ---- leave-fold-out choice of (alpha, T), a 9-cell grid ----
    fo = np.array([c["fold"] for c in caches])
    cells = {k: np.array(v) for k, v in arms.items() if k.startswith("vqe_")}
    lfo, picks = np.zeros(len(caches)), {}
    for held in sorted(set(fo.tolist())):
        tr, te = fo != held, fo == held
        best = min(cells, key=lambda k: cells[k][tr].mean())
        picks[int(held)] = best
        lfo[te] = cells[best][te]
    res["arms"]["vqe_LFO"] = _arm_summary(lfo, [], [], [], [], [], b)
    res["arms"]["vqe_LFO"]["picks"] = picks
    res["per_target"]["vqe_LFO"] = [float(x) for x in lfo]
    path = out or os.path.join(HERE, "integrate_vqe.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nCVaR-VQE  n={n} qubits ({dim} hypotheses), {layers} layers, "
              f"filter={filt}, n_targets={len(caches)}, baseline {b.mean():.3f}\n")
        print("%-18s%9s%9s%9s%8s%9s%10s" %
              ("arm", "sel", "d", "CI_lo", "CI_hi", "W/L", "H(p) bits"))
        for k in sorted(res["arms"], key=lambda z: res["arms"][z]["sel"]):
            d = res["arms"][k]
            v_ = d["vs_dist"]
            print("%-18s%9.3f%9.3f%9.3f%8.3f%9s%10.2f" % (
                k, d["sel"], v_["mean_diff"], v_["ci95"][0], v_["ci95"][1],
                "%d/%d" % (v_["wins"], v_["losses"]),
                res["entropy_bits"].get(k, float("nan"))))
    return res


# ============================================================ 5. THE ROLES
# The four mandated components have to be materially present, and "materially" has to mean
# measured.  This stage puts Legacy and Amber into the two roles the sprint's own data
# suggests for them -- neither of which is "a global selector", which is where both have
# already been measured and failed -- and prices each one.
#
#   Legacy: S8-8 measured its 15-column learned combiner at 2.228 A on an ORACLE-tight pool
#           against the shipped score's 2.406, and 4.103 A on the normal pool against
#           3.454.  So its honest role is LATE, on an already-narrowed set.  The deployable
#           analogue of a tight pool is the neighbourhood of the consensus medoid, and that
#           is what is tested here.
#   Amber:  failed as a single-point global ranker (native percentile 54.3, argmin 0/70).
#           It is affordable only on a shortlist, so it is tested where it is affordable:
#           as a re-ranker and as an ensemble weight inside the BLOSUM top-32.
def _legacy_lfo(caches, cols=None):
    cols = cols or [n for n in FAST if n.startswith("lg_")]
    Wf = _fit_lfo(caches, cols)
    out = {}
    for c in caches:
        w = Wf.get(c["fold"])
        if w is not None:
            out[c["pdb"]] = Fuser(cols, w)(_cols_of(c, cols))
    return out, cols


def stage_role(out=None, verbose=True, q=75):
    """Legacy as a late refiner; Legacy and Amber as ensemble weights."""
    tgs = chan_cached()
    caches = [load_chan(p) for p in tgs]
    per = _filters(caches)["score"]
    lg_lfo, lg_cols = _legacy_lfo(caches)
    base = np.array([sel_of(c["S"][:, c["names"].index("score")], c["rr"]) for c in caches])
    arms = {}

    def add(name, v):
        arms[name] = np.asarray(v, float)

    med, med_lg, med_lgall, med_ty = [], [], [], []
    ref = {k: [] for k in (5, 10, 20)}
    ref_all = {k: [] for k in (5, 10, 20)}
    for c in caches:
        pid = c["pdb"]
        s = np.asarray(per[pid], float)
        rr, D = c["rr"], load_dmat(pid)
        o = np.argsort(s, kind="stable")[:q]
        m = consensus_medoid(D, o)
        med.append(rr[m])
        for nm, col, dst in (("lg", None, med_lg), ("lgall", "lg_all", med_lgall),
                             ("ty", "typic", med_ty)):
            v = lg_lfo[pid][o] if col is None else c["S"][o, c["names"].index(col)]
            w = np.exp(-zrank(v))
            dst.append(rr[consensus_medoid(D, o, w)])
        # LATE REFINEMENT: the k members closest to the consensus medoid form a genuinely
        # tight set (built with no native); rank inside it with Legacy.
        for k in ref:
            nb = o[np.argsort(D[m][o])[:k]]
            ref[k].append(rr[nb[int(np.argmin(lg_lfo[pid][nb]))]])
            ref_all[k].append(rr[nb[int(np.argmin(c["S"][nb, c["names"].index("lg_all")]))]])
    add(f"consensus@{q}", med)
    add(f"consensus@{q}_wLegacyLFO", med_lg)
    add(f"consensus@{q}_wLegacyAll", med_lgall)
    add(f"consensus@{q}_wTypic", med_ty)
    for k in ref:
        add(f"refine_LegacyLFO_top{k}", ref[k])
        add(f"refine_LegacyAll_top{k}", ref_all[k])

    res = {"n_targets": len(caches), "q": q, "legacy_cols": lg_cols,
           "baseline_sel": float(base.mean()), "arms": {}}
    for k, v in arms.items():
        res["arms"][k] = _arm_summary(v, [], [], [], [], [], base)

    # ---------------------------------------------------------------- AMBER
    ata = [p for p in tgs if os.path.exists(_apath(p))]
    if ata:
        cA = [load_chan(p, with_amber=True) for p in ata]
        m = min(c["amber_m"] for c in cA)
        a_arms, npct, nfirst = {}, [], []
        AW = (0.25, 0.5, 1.0)
        rows = {k: [] for k in
                ["argmin_score", "argmin_amber", "consensus8", "consensus16",
                 "consensus16_wAmber", "pool_mean", "pool_best"]
                + [f"argmin_score+{w}amber" for w in AW]
                + [f"consensus16_filt_score+{w}amber" for w in AW]}
        for c in cA:
            rr = c["rr"][:m]
            D = load_dmat(c["pdb"])[:m, :m]
            s = c["S"][:m, c["names"].index("score")]
            e = c["S"][:m, c["names"].index("amber")]
            ok = np.isfinite(e)
            if ok.sum() < 8:
                continue
            idx = np.arange(m)[ok]
            rows["argmin_score"].append(rr[idx[int(np.argmin(s[idx]))]])
            rows["argmin_amber"].append(rr[idx[int(np.argmin(e[idx]))]])
            rows["pool_mean"].append(float(rr[idx].mean()))
            rows["pool_best"].append(float(rr[idx].min()))
            zs, ze = zrank(s[idx]), zrank(e[idx])
            for w in AW:
                f = zs + w * ze
                rows[f"argmin_score+{w}amber"].append(rr[idx[int(np.argmin(f))]])
                of = idx[np.argsort(f, kind="stable")[:16]]
                rows[f"consensus16_filt_score+{w}amber"].append(
                    float(rr[consensus_medoid(D, of)]))
            for qq, nm in ((8, "consensus8"), (16, "consensus16")):
                o = idx[np.argsort(s[idx], kind="stable")[:qq]]
                rows[nm].append(float(rr[consensus_medoid(D, o)]))
                if nm == "consensus16":
                    w = np.exp(-zrank(e[o]))
                    rows["consensus16_wAmber"].append(
                        float(rr[consensus_medoid(D, o, w)]))
            if np.isfinite(c["amber_native"]):
                npct.append(100.0 * float((e[idx] < c["amber_native"]).mean()))
                nfirst.append(bool(c["amber_native"] < e[idx].min()))
        b2 = np.array(rows["argmin_score"])
        res["amber"] = {"n_targets": len(b2), "m": int(m),
                        "native_pct": _agg(npct), "native_argmin": int(sum(nfirst)),
                        "arms": {k: _arm_summary(np.array(v), [], [], [], [], [], b2)
                                 for k, v in rows.items() if v}}
    path = out or os.path.join(HERE, "integrate_role.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nROLES  n={res['n_targets']}  shipped baseline {res['baseline_sel']:.3f}\n")
        print("%-32s%9s%9s%9s%8s%9s" % ("arm", "sel", "d", "CI_lo", "CI_hi", "W/L"))
        for k in sorted(res["arms"], key=lambda z: res["arms"][z]["sel"]):
            d = res["arms"][k]
            v = d["vs_dist"]
            print("%-32s%9.3f%9.3f%9.3f%8.3f%9s" % (
                k, d["sel"], v["mean_diff"], v["ci95"][0], v["ci95"][1],
                "%d/%d" % (v["wins"], v["losses"])))
        if "amber" in res:
            a = res["amber"]
            print(f"\nAMBER, restricted to the BLOSUM top-{a['m']} "
                  f"({a['n_targets']} targets with Amber cached)")
            print(f"  native percentile in the Amber energy distribution: "
                  f"{a['native_pct']['mean']:.1f} (median "
                  f"{a['native_pct']['median']:.1f}), argmin at native "
                  f"{a['native_argmin']}/{a['native_pct']['n']}")
            print("%-32s%9s%9s%9s%8s%9s" % ("arm", "sel", "d", "CI_lo", "CI_hi", "W/L"))
            for k in sorted(a["arms"], key=lambda z: a["arms"][z]["sel"]):
                d = a["arms"][k]
                v = d["vs_dist"]
                print("%-32s%9.3f%9.3f%9.3f%8.3f%9s" % (
                    k, d["sel"], v["mean_diff"], v["ci95"][0], v["ci95"][1],
                    "%d/%d" % (v["wins"], v["losses"])))
    return res


# ============================================================ 6. THE SYSTEM + ABLATIONS
#: The integrated system, stated as a composition so each ablation removes exactly one
#: thing.  Every stage below is DEPLOYABLE: no native coordinate, distance, torsion or
#: RMSD enters any of them.
#:
#:   1 RETRIEVE   BLOSUM top-500                       (S8-6 closed this axis)
#:   2 OBJECTIVE  multi-objective score over the four+one channels
#:   3 HYPOTHESES the objective's top `dim` candidates = the discrete state space
#:   4 ENSEMBLE   CVaR-VQE free energy -> p_theta      (VQE + CVaR)
#:   5 READOUT    p_theta-weighted consensus medoid
#:   6 RERANK     Legacy / AMBER on the narrowed set   (measured, and reported either way)
SYS_FILTER = [("score", 1.0), ("lg_all", 0.25), ("iv_nbid", 0.25)]


def stage_final(out=None, verbose=True, n=7, layers=3, iters=50, q=75):
    """The integrated system against 3.454 A, with one ablation per component."""
    vpath = os.path.join(HERE, "integrate_vqe.json")
    V = json.load(open(vpath)) if os.path.exists(vpath) else None
    tgs = chan_cached()
    caches = [load_chan(p) for p in tgs]
    F = _filters(caches)
    b = np.array([sel_of(c["S"][:, c["names"].index("score")], c["rr"]) for c in caches])
    dim = 1 << int(n)

    def readout(per, mode, qq=q, weights=None):
        v = []
        for c in caches:
            s = np.asarray(per[c["pdb"]], float)
            rr, D = c["rr"], load_dmat(c["pdb"])
            if mode == "argmin":
                v.append(sel_of(s, rr))
            else:
                o = np.argsort(s, kind="stable")[:qq]
                w = None if weights is None else weights(c, o, s)
                v.append(float(rr[consensus_medoid(D, o, w)]))
        return np.array(v)

    arms = {}
    arms["SHIPPED distogram argmin"] = b
    arms["objective argmin (no consensus)"] = readout(F["D+L+I.25"], "argmin")
    arms["FULL: multi-obj filter + consensus"] = readout(F["D+L+I.25"], "cons")
    arms["- Legacy from the objective"] = readout(F["D+I.25"], "cons")
    arms["- invfold from the objective"] = readout(F["D+L.25"], "cons")
    arms["- both (shipped score filter)"] = readout(F["score"], "cons")
    arms["+ typicality in the objective"] = readout(F["D+T+I.25"], "cons")
    arms["- consensus, + Legacy rerank"] = readout(F["D+L+I.25"], "argmin")
    # the composition the criterion-fusion study selects: filter with the shipped score,
    # then rank inside the filter by (consensus criterion + 0.5 * inverse folding)
    cf = []
    for c in caches:
        s = np.asarray(F["score"][c["pdb"]], float)
        D = load_dmat(c["pdb"])
        o = np.argsort(s, kind="stable")[:q]
        f = zrank(D[np.ix_(o, o)].astype(float).mean(1)) + \
            0.5 * zrank(c["S"][o, c["names"].index("iv_nbid")])
        cf.append(float(c["rr"][o[int(np.argmin(f))]]))
    arms["SYSTEM: crit@75 + 0.5*invfold"] = np.array(cf)
    if V and "per_target" in V and V["pdbs"] == tgs:
        for k in ("vqe_LFO", f"medoid{V['dim']}", "argmin"):
            if k in V["per_target"]:
                arms[f"VQE[{k}] over top-{V['dim']}"] = np.array(V["per_target"][k])
        for k in V["per_target"]:
            if k.startswith(("topfrac_", "boltz_")):
                arms[f"ablate {k} (no VQE)"] = np.array(V["per_target"][k])
    res = {"n_targets": len(caches), "q": q, "baseline_sel": float(b.mean()),
           "system_filter": SYS_FILTER, "arms": {},
           "pdbs": [c["pdb"] for c in caches],
           "folds": [c["fold"] for c in caches],
           "per_target": {k: [float(x) for x in v] for k, v in arms.items()}}
    for k, v in arms.items():
        res["arms"][k] = _arm_summary(v, [], [], [], [], [], b)
    # ---- what each component is WORTH ----
    # Stated as `system_without_it - system_with_it`, so a POSITIVE number means the
    # component earns its place and a negative one means it costs.  Each pair differs in
    # exactly one component and is evaluated on identical candidates.
    res["ablations"] = {}
    for lbl, with_, without in (
            ("consensus (vs argmin of the same objective)",
             "FULL: multi-obj filter + consensus", "objective argmin (no consensus)"),
            ("Legacy in the objective",
             "FULL: multi-obj filter + consensus", "- Legacy from the objective"),
            ("inverse folding in the objective",
             "FULL: multi-obj filter + consensus", "- invfold from the objective"),
            ("the whole multi-objective fusion (vs the shipped score alone)",
             "FULL: multi-obj filter + consensus", "- both (shipped score filter)"),
            ("Legacy on top of the shipped score",
             "- invfold from the objective", "- both (shipped score filter)"),
            ("inverse folding on top of the shipped score",
             "- Legacy from the objective", "- both (shipped score filter)"),
            ("typicality on top of the full objective",
             "+ typicality in the objective", "FULL: multi-obj filter + consensus"),
            ("inverse folding on top of the CONSENSUS CRITERION",
             "SYSTEM: crit@75 + 0.5*invfold", "- both (shipped score filter)")):
        if with_ in arms and without in arms:
            p = _paired(arms[without], arms[with_])       # without - with
            res["ablations"][lbl] = p
    path = out or os.path.join(HERE, "integrate_final.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nINTEGRATED SYSTEM  n={len(caches)}  vs the shipped {b.mean():.3f} A\n")
        print("%-40s%9s%8s%9s%9s%8s%9s%8s%8s" %
              ("arm", "sel", "med", "d", "CI_lo", "CI_hi", "W/L", "<2A", "<1.5A"))
        for k in sorted(res["arms"], key=lambda z: res["arms"][z]["sel"]):
            d = res["arms"][k]
            v = d["vs_dist"]
            print("%-40s%9.3f%8.3f%9.3f%9.3f%8.3f%9s%8.3f%8.3f" % (
                k, d["sel"], d["sel_median"], v["mean_diff"], v["ci95"][0], v["ci95"][1],
                "%d/%d" % (v["wins"], v["losses"]), d["frac_lt2"], d["frac_lt15"]))
        print("\nWHAT EACH COMPONENT IS WORTH  (without it) - (with it); "
              "positive = it earns its place)")
        for k, v in res["ablations"].items():
            print("  %-52s %+.3f A  CI [%+.3f, %+.3f]  helps/hurts %d/%d" % (
                k, v["mean_diff"], v["ci95"][0], v["ci95"][1],
                v["losses"], v["wins"]))     # `losses` = targets where removing it HURT
    return res


def mem_pct():
    """Percent of physical RAM in use, via WMI. The box is shared with three agents."""
    try:
        import subprocess
        out = subprocess.run(
            ["powershell", "-NoProfile", "-c",
             "$o=Get-CimInstance Win32_OperatingSystem;"
             "[int](100*(1-$o.FreePhysicalMemory/$o.TotalVisibleMemorySize))"],
            capture_output=True, text=True, timeout=60).stdout.strip()
        return int(out)
    except Exception:
        return -1


def stage_refine(out=None, verbose=True, q=75, limit=None, cap=90):
    """Does an AMBER relaxation of the SELECTED structure improve it?

    The one role for Amber this study can test end to end on its own output: take what the
    system returns, relax it, and measure the CA-RMSD change.  `s8/relax.py` (another
    agent) is measuring the same operation across a k x steps sweep on pool candidates;
    this is the narrow version on the structures this system actually selects, which is the
    only version that can be quoted as this system's Amber ablation.

    Unrestrained convergence (`k=0, steps=0`) and a moderate restraint are both run, since
    the two answer different questions: whether the force field's own minimum is closer to
    the native, and whether a bounded local relaxation is.
    """
    import amber_refine as ar
    import torsion_lib2 as tl2
    from s5.lib import kabsch_rmsd_batch
    tgs = [p for p in chan_cached() if os.path.exists(_apath(p))]
    if limit:
        tgs = tgs[:int(limit)]
    rows = []
    for pid in tgs:
        # The box is shared with three agents and has hit 97% once already.  Wait rather
        # than abandon: a single OpenMM context is ~200 MB, so the right behaviour is to
        # yield until there is room, and to give up only if there never is.
        waits = 0
        while mem_pct() > int(cap) and waits < 20:
            time.sleep(60)
            waits += 1
        if mem_pct() > int(cap):
            print(f"memory above {cap}% for 20 min, stopping after {len(rows)} targets",
                  flush=True)
            break
        u = gen.load_univ(pid)
        P = pool(u)
        c = load_chan(pid)
        s = c["S"][:, c["names"].index("score")]
        D = load_dmat(pid)
        o = np.argsort(s, kind="stable")[:q]
        pick = consensus_medoid(D, o)
        seq = str(u["seq"])
        tab = tl2.library_for(seq, 8, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        cs = {k: v[pick] for k, v in P["BB"].items()}
        # The relaxation starts from the ideal-geometry REBUILD, which sits ~0.3-0.4 A
        # from the real window whose RMSD the pipeline reports.  Differencing against the
        # window's RMSD would charge the rebuild displacement to Amber, so both baselines
        # are recorded and `d_rmsd` below is measured against the rebuild.
        r = {"pdb": pid, "rmsd_window": float(P["rr"][pick]),
             "rmsd0": float(kabsch_rmsd_batch(P["BB"]["CA"][pick][None],
                                              u["nat_ca"])[0]),
             "argmin_rmsd": sel_of(s, P["rr"])}
        for tag, k_r in (("free", 0.0), ("k10", 10.0)):
            try:
                res = ar.refine_coords(seq, rep, cs, k_restraint=k_r, steps=0,
                                       tolerance=1.0, components=True)
                ca = np.asarray(res["ca"], float)
                r[f"rmsd1_{tag}"] = float(kabsch_rmsd_batch(ca[None], u["nat_ca"])[0])
            except Exception:
                r[f"rmsd1_{tag}"] = None
        rows.append(r)
        if verbose:
            print(f"{pid} {r['rmsd0']:.3f} -> free {r.get('rmsd1_free')} "
                  f"k10 {r.get('rmsd1_k10')}", flush=True)
    res = {"n_targets": len(rows), "q": q, "rows": rows, "agg": {}}
    for tag in ("free", "k10"):
        d = np.array([r[f"rmsd1_{tag}"] - r["rmsd0"] for r in rows
                      if r.get(f"rmsd1_{tag}") is not None])
        if len(d):
            res["agg"][tag] = {"d_rmsd": float(d.mean()),
                               "se": float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else 0.0,
                               "improved": int((d < 0).sum()), "n": int(len(d))}
    path = out or os.path.join(HERE, "integrate_refine.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print("\nAMBER relaxation of the SELECTED structure (consensus@%d), n=%d" %
              (q, len(rows)))
        for k, v in res["agg"].items():
            print("  %-6s d_rmsd %+.4f +- %.4f  improved %d/%d" %
                  (k, v["d_rmsd"], v["se"], v["improved"], v["n"]))
    return res


def stage_critfuse(out=None, verbose=True):
    """Fuse the CONSENSUS CRITERION with the channels -- the one fusion the matrix allows.

    The error-correlation matrix says a fusion only pays when the second member's skill is
    comparable to the first's, and every channel measured here is 2-5x weaker than the
    distogram, which is why every score-level fusion was null.  The consensus criterion
    breaks that: inside the score's own top 75 it ranks at +0.224 against the score's
    +0.065 on identical candidates.  So this is the one place the precondition is met, and
    it is tested with leave-fold-out arm selection because a 21-cell grid on a 0.075 A SE
    will otherwise hand back its best cell as a "result".
    """
    tgs = chan_cached()
    caches = [load_chan(p) for p in tgs]
    F = _filters(caches)
    folds = {c["pdb"]: c["fold"] for c in caches}
    base = np.array([sel_of(c["S"][:, c["names"].index("score")], c["rr"]) for c in caches])
    # One pre-declared grid, evaluated once.  Single-channel additions plus the three
    # two-channel ones; the LFO number below is the reported result whatever it says.
    SPEC = ([()] + [((ch, w),) for ch in ("score", "lg_all", "iv_nbid", "typic")
                    for w in (0.25, 0.5)]
            + [(("score", a), ("iv_nbid", b)) for a, b in
               ((0.25, 0.25), (0.25, 0.5), (0.5, 0.5))])
    cells = {}
    for q in (50, 75, 100):
        for spec in SPEC:
            sel = {}
            for c in caches:
                s = np.asarray(F["score"][c["pdb"]], float)
                D = load_dmat(c["pdb"])
                o = np.argsort(s, kind="stable")[:q]
                f = zrank(D[np.ix_(o, o)].astype(float).mean(1))
                for ch, w in spec:
                    v = s[o] if ch == "score" else c["S"][o, c["names"].index(ch)]
                    f = f + w * zrank(v)
                sel[c["pdb"]] = float(c["rr"][o[int(np.argmin(f))]])
            cells[f"crit@{q}" + "".join(f"+{w}*{ch}" for ch, w in spec)] = sel
    res = {"n_targets": len(caches), "baseline_sel": float(base.mean()), "cells": {}}
    for k, v in cells.items():
        res["cells"][k] = _arm_summary(np.array([v[p] for p in tgs]), [], [], [], [],
                                       [], base)
    lfo, picks = {}, {}
    for held in sorted(set(folds.values())):
        tr = [p for p in tgs if folds[p] != held]
        best = min(cells, key=lambda k: np.mean([cells[k][p] for p in tr]))
        picks[int(held)] = best
        for p in [p for p in tgs if folds[p] == held]:
            lfo[p] = cells[best][p]
    res["critfuse_LFO"] = _arm_summary(np.array([lfo[p] for p in tgs]), [], [], [], [],
                                       [], base)
    res["critfuse_LFO"]["picks"] = picks
    path = out or os.path.join(HERE, "integrate_critfuse.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nCONSENSUS-CRITERION FUSION  n={len(caches)}  baseline {base.mean():.3f}\n")
        print("%-24s%9s%9s%9s%8s%9s%8s" %
              ("cell", "sel", "d", "CI_lo", "CI_hi", "W/L", "<2A"))
        for k in sorted(res["cells"], key=lambda z: res["cells"][z]["sel"]):
            d = res["cells"][k]
            v = d["vs_dist"]
            print("%-24s%9.3f%9.3f%9.3f%8.3f%9s%8.3f" % (
                k, d["sel"], v["mean_diff"], v["ci95"][0], v["ci95"][1],
                "%d/%d" % (v["wins"], v["losses"]), d["frac_lt2"]))
        d = res["critfuse_LFO"]
        v = d["vs_dist"]
        print("\nLEAVE-FOLD-OUT over the %d cells: %.3f  d %+.3f [%+.3f, %+.3f] W/L %d/%d"
              % (len(cells), d["sel"], v["mean_diff"], v["ci95"][0], v["ci95"][1],
                 v["wins"], v["losses"]))
        print("picks per fold:", picks)
    return res


def stage_k25(out=None, verbose=True, K25=25):
    """The multi-objective score in the ONE regime S8-6 left open: a 25-candidate pool.

    S8-6 measured that perfect retrieval is worth 3.5x more at K=25 than at K=500 -- "the
    value of retrieval is inversely proportional to how many ways the selector has to pick
    badly" -- and 3.194 A there is the best number in that study.  The same argument
    applies to recognition, and no channel has been tested at K=25 on this instrument.  The
    pool is the BLOSUM top-25, i.e. a genuine retrieval configuration and not a filtered
    slice of the top-500.
    """
    tgs = chan_cached()
    caches = [load_chan(p) for p in tgs]
    F = _filters(caches)
    m = int(K25)
    rows, base = {}, []
    def put(k, v):
        rows.setdefault(k, []).append(v)
    for c in caches:
        rr = c["rr"][:m]
        D = load_dmat(c["pdb"])[:m, :m]
        base.append(sel_of(c["S"][:m, c["names"].index("score")], rr))
        put("pool_best", float(rr.min()))
        put("pool_mean", float(rr.mean()))
        for nm in ("score", "typic", "legacy", "invfold", "D+L.25", "D+I.25",
                   "D+L+I.25", "D+T.25", "lfo_DLI", "lfo_DLIT"):
            if c["pdb"] not in F[nm]:
                continue
            s = np.asarray(F[nm][c["pdb"]], float)[:m]
            put(f"argmin_{nm}", sel_of(s, rr))
            for q in (8, 12, 25):
                o = np.argsort(s, kind="stable")[:q]
                put(f"consensus{q}_{nm}", float(rr[consensus_medoid(D, o)]))
    b = np.array(base)
    res = {"n_targets": len(caches), "K": m, "baseline_sel": float(b.mean()), "arms": {}}
    for k, v in rows.items():
        res["arms"][k] = _arm_summary(np.array(v), [], [], [], [], [], b)
    path = out or os.path.join(HERE, "integrate_k25.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print(f"\nK={m} BLOSUM pool, n={len(caches)}, baseline {b.mean():.3f}\n")
        print("%-30s%9s%9s%9s%8s%9s" % ("arm", "sel", "d", "CI_lo", "CI_hi", "W/L"))
        for k in sorted(res["arms"], key=lambda z: res["arms"][z]["sel"])[:18]:
            d = res["arms"][k]
            v = d["vs_dist"]
            print("%-30s%9.3f%9.3f%9.3f%8.3f%9s" % (
                k, d["sel"], v["mean_diff"], v["ci95"][0], v["ci95"][1],
                "%d/%d" % (v["wins"], v["losses"])))
    return res


def stage_natcons(out=None, verbose=True, qs=(25, 75, 150, 500)):
    """ORACLE DIAGNOSTIC: where does the NATIVE sit under the consensus criterion?

    Native percentile is the sprint's standard test of whether an objective is even
    correctly specified, and it has been reported for every SCORE.  The consensus medoid is
    not a score, so the analogue has to be constructed: drop the native into the filtered
    set and rank it by the same criterion the medoid minimises -- mean CA-RMSD to the other
    members.  0 = the native IS the consensus; 50 = a coin flip.

    This is the one diagnostic that says whether the consensus result is a real signal
    about nativeness or just an outlier-avoidance effect, and the two make opposite
    predictions: a real signal puts the native near 0, outlier avoidance puts it near
    whatever percentile a typical pool member occupies.
    """
    from s5.lib import kabsch_rmsd_batch
    tgs = chan_cached()
    rows = []
    for pid in tgs:
        u = gen.load_univ(pid)
        P = pool(u)
        c = load_chan(pid)
        s = c["S"][:, c["names"].index("score")]
        D = load_dmat(pid)
        dn = kabsch_rmsd_batch(P["W"], np.asarray(u["nat_ca"], float))
        r = {"pdb": pid}
        for q in qs:
            o = np.argsort(s, kind="stable")[:int(q)]
            mem = D[np.ix_(o, o)].astype(float).mean(1)
            # the native's mean distance to the SAME members (it is not one of them, so
            # its own zero self-distance is excluded from both sides by construction)
            nat = float(dn[o].mean())
            r[f"pct@{q}"] = 100.0 * float((mem < nat).mean())
            r[f"medoid_rmsd@{q}"] = float(P["rr"][o[int(np.argmin(mem))]])
            r[f"best_rmsd@{q}"] = float(P["rr"][o].min())
            # ...and how well the SAME criterion ranks the candidates, against how well
            # the shipped score ranks them inside its own filtered set.  The two numbers
            # together are the finding: the criterion ranks candidates several times
            # better than the score does and still puts the native last.
            r[f"crit_rho@{q}"] = _spearman(mem, P["rr"][o])
            r[f"crit_prg@{q}"] = partial_rg(mem, P["rr"][o], c["rg"][o])
            r[f"score_rho@{q}"] = _spearman(s[o], P["rr"][o])
        rows.append(r)
    res = {"n_targets": len(rows), "qs": list(qs), "rows": rows,
           "agg": {k: _agg([r[k] for r in rows]) for k in rows[0] if k != "pdb"}}
    path = out or os.path.join(HERE, "integrate_natcons.json")
    json.dump(res, open(path, "w"), indent=1)
    if verbose:
        print("\nWhere the NATIVE sits under the consensus criterion "
              f"(n={len(rows)}; 0 = the native is the consensus, 50 = coin flip)\n")
        print("%-6s%11s%9s%12s%10s%11s%10s%11s" %
              ("q", "nativePct", "median", "medoidRMSD", "setBest",
               "critRho", "critPrg", "scoreRho"))
        for q in qs:
            a = res["agg"][f"pct@{q}"]
            print("%-6d%11.1f%9.1f%12.3f%10.3f%11.3f%10.3f%11.3f" % (
                q, a["mean"], a["median"], res["agg"][f"medoid_rmsd@{q}"]["mean"],
                res["agg"][f"best_rmsd@{q}"]["mean"], res["agg"][f"crit_rho@{q}"]["mean"],
                res["agg"][f"crit_prg@{q}"]["mean"], res["agg"][f"score_rho@{q}"]["mean"]))
    return res


def main(argv):
    cmd = argv[1] if len(argv) > 1 else "help"
    rest = argv[2:]
    fn = globals().get("stage_" + cmd)
    if fn is None:
        print(__doc__)
        return 1
    kw = {}
    for a in rest:
        if "=" in a:
            k, v = a.split("=", 1)
            try:
                v = int(v)
            except ValueError:
                try:
                    v = float(v)
                except ValueError:
                    pass
            kw[k] = v
    fn(**kw)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
