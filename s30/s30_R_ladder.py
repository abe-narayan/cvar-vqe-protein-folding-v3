#!/usr/bin/env python
"""s30/s30_R_ladder.py -- LANE R: IS NATIVENESS RECOGNISABLE FROM SINGLE-STRUCTURE GEOMETRY?

Pre-registered in `s30/PREREG_S30_R.md` (commit 7eabffee, 2026-09-20 12:40:35), BEFORE any
aggregate existed.  Charter lead L11.

THE INSTRUMENT.  A ladder of structures in which KIND, LOCAL REALISM and PERTURBATION BUDGET are
all matched, and only nativeness varies.  Every rung is an ideal-geometry backbone built from
(phi, psi) by `core.geometry.build_backbone_batch` -- no projection, no coordinate average, no
contraction -- and every perturbed residue draws its torsions from the FOLD'S LEAKAGE-SAFE
Ramachandran table `s8/generate_rama.npz[fold]`, the same table `ham_lib.h_rama` scores against.

  ladder A   anchor = the native's own (phi, psi);   m in M_GRID residues resampled, R draws each
  ladder B   anchor = a random REAL pool member;     same grid, same R
  fixed      the rebuilt native; the projected PROD / circ_best / sub0 chains; the 500 pool members

The ladder's floor is the torsion rebuild, ~0.3 A from `nat_ca`, NOT 0 A.

EVERY RMSD COLUMN IS ORACLE (it reads `cand.nat_ca`) and is used only to LABEL a rung after all
channel values have been computed.  No channel sees a native: audit A3 re-runs the scoring pass
with `nat_ca` NaN'd and requires bit-identical channel values.

    python s30/s30_R_ladder.py prep  [--shard i --of n] [--limit N]     # project the fixed rungs
    python s30/s30_R_ladder.py run   [--shard i --of n] [--limit N]     # generate + score + row
    python s30/s30_R_ladder.py agg                                      # the verdict on F-R1
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
CHAINCACHE = os.path.join(RESULTS, "s30_R_chains")
LADDER_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_D_ladder_structs")
RAMA_PATH = os.path.join(ROOT, "s8", "generate_rama.npz")

M_GRID = (1, 2, 3, 4, 6, 8, 12)
R_DRAWS = 32
NEAR_NATIVE_A = 1.0                     # clause (ii)'s "<= 1 A" band on ladder A
FIXED_RUNGS = ("PROD", "circ_best", "sub0")
DELTA_BINS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 1e9)
DELTA_LABELS = ("0-0.25", "0.25-0.5", "0.5-1", "1-2", "2-4", ">4")

#: Channels whose SIZE-MATCHED twin is well defined (lane L's relay, S30).  A statistical pair
#: potential u(d) = -kT ln[P_obs/P_ref] with a separation-only or quasi-chemical reference
#: (`ham_lib.h_distpot`'s `Pref`, `h_contact`'s `Nexp`) absorbs chain connectivity and composition
#: but NOT size, so it carries a separable +kT sum ln P_ref(d) term that is a pure function of
#: scale.  The twin removes it BY CONSTRUCTION rather than by partialling: the universe AND the
#: candidates are both rescaled to a common radius of gyration, so the potential is fitted and
#: evaluated in REDUCED UNITS and a uniform contraction is exactly invisible.  RG_LAW and RG_UNIV
#: must come out CONSTANT under the twin -- that is the audit that the transform does what it says.
SI_CHANNELS = ("RG_LAW", "RG_UNIV", "EXVOL", "CAGEO", "CONTACT", "DISTPOT", "ENV", "HP",
               "CONS", "DMAP_CONS", "POOLGO", "SS_MATCH")
#: NOT given a twin: the distogram re-readings (DIS, DIS_MEAN, CONTACT_LL) score against a
#: posterior tabulated in absolute Angstroms, so a rescaled candidate is not the same question;
#: and the torsion-space channels (RAMA, DSSPHB, LEG_*, TORS_CONS) are scale-free already because
#: they read (phi, psi), which a uniform rescale does not touch.

_RAMA = None
_PEP = None


# ------------------------------------------------------------------ small helpers
def rama_cnt():
    global _RAMA
    if _RAMA is None:
        _RAMA = np.load(RAMA_PATH)["cnt"]           # (5, 20, 36, 36)
    return _RAMA


def native_torsions(pdb):
    """The native's OWN deposited (phi, psi) from the peptide bank."""
    global _PEP
    if _PEP is None:
        from core import data as cdata
        recs = np.load(cdata.PEPTIDE_CACHE, allow_pickle=True)["records"]
        _PEP = {r["pdb"]: r for r in recs}
    r = _PEP.get(pdb)
    if r is None:
        return None
    return np.asarray(r["phi"], float), np.asarray(r["psi"], float), r["seq"]


def rg_of(W):
    W = np.asarray(W, float)
    c = W.mean(1, keepdims=True)
    return np.sqrt(((W - c) ** 2).sum(-1).mean(1))


def ca_rmsd_to(W, ref):
    from s12 import instrument as I
    return np.asarray(I.kabsch_rmsd_batch(np.asarray(W, float), np.asarray(ref, float)), float)


def spearman(x, y):
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.size < 4 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    r = spearmanr(x, y).correlation
    return float(r) if np.isfinite(r) else float("nan")


def partial_spearman(x, y, Zcols, vshape_of=None):
    """Spearman(x, y) with the rank of every column of Zcols removed (and, for each column named
    in vshape_of, also the rank of |col - median col|).  S29-L50's repaired partial."""
    from scipy.stats import rankdata
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.size < 6 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    cols = []
    for idx, z in enumerate(Zcols):
        z = np.asarray(z, float)
        cols.append(rankdata(z))
        if vshape_of is not None and idx in vshape_of:
            cols.append(rankdata(np.abs(z - np.median(z))))
    Z = np.column_stack(cols)
    Z = Z - Z.mean(0, keepdims=True)
    rx = rankdata(x) - np.mean(rankdata(x))
    ry = rankdata(y) - np.mean(rankdata(y))
    try:
        cx = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
        cy = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    except np.linalg.LinAlgError:
        return float("nan")
    den = math.sqrt(float(cx @ cx) * float(cy @ cy))
    return float(cx @ cy / den) if den > 1e-12 else float("nan")


# ------------------------------------------------------------------ the generator
class Sampler:
    """Draws (phi, psi) for residue i from the FOLD'S LEAKAGE-SAFE Ramachandran table."""

    def __init__(self, fold, seq, pseudo=1.0):
        from core import data as cdata
        aidx = {a: i for i, a in enumerate(cdata.ALPHABET)}
        self.aa = np.array([aidx[c] for c in seq], int)
        cnt = rama_cnt()[int(fold)].astype(float) + pseudo
        P = cnt.reshape(20, -1)
        self.C = np.cumsum(P / P.sum(1, keepdims=True), 1)
        self.step = 2 * math.pi / 36

    def draw(self, i, rng):
        k = int(np.searchsorted(self.C[self.aa[i]], rng.random()))
        pb, sb = k // 36, k % 36
        return (-math.pi + (pb + rng.random()) * self.step,
                -math.pi + (sb + rng.random()) * self.step)


def make_ladder(phi0, psi0, sampler, n, rng, m_grid, R):
    """Return PHI (S, n), PSI (S, n), m (S,), nchanged (S,) -- the budget-matched ladder."""
    PHI, PSI, MS = [], [], []
    for m in m_grid:
        mm = min(int(m), n)
        for _ in range(R):
            p, s = phi0.copy(), psi0.copy()
            idx = rng.choice(n, size=mm, replace=False)
            for i in idx:
                p[i], s[i] = sampler.draw(int(i), rng)
            PHI.append(p); PSI.append(s); MS.append(mm)
    return np.array(PHI), np.array(PSI), np.array(MS, int)


# ------------------------------------------------------------------ stage 1: project fixed rungs
def chain_path(pdb):
    return os.path.join(CHAINCACHE, f"{pdb}.npz")


def prep_target(pdb, rebuild=False):
    from s12 import instrument as I
    from s24 import d_harness as H
    f = chain_path(pdb)
    if os.path.exists(f) and not rebuild:
        with np.load(f) as z:
            if all(f"phi_{k}" in z.files for k in FIXED_RUNGS):
                return "cached"
    cand = H.Candidates.from_universe(pdb)
    z = np.load(os.path.join(LADDER_STRUCTS, f"{pdb}.npz"))
    out = {}
    for k in FIXED_RUNGS:
        key = k.replace("[", "_").replace("]", "").replace(".", "_")
        pr = I.project(np.asarray(z[key], float), cand.seq, int(cand.fold))
        out[f"ca_{k}"] = np.asarray(pr["ca"], float)
        out[f"phi_{k}"] = np.asarray(pr["phi"], float)
        out[f"psi_{k}"] = np.asarray(pr["psi"], float)
        out[f"rmsdca_{k}"] = np.array(float(z["rmsd_" + key]))               # ORACLE, cached label
    os.makedirs(CHAINCACHE, exist_ok=True)
    tmp = f + f".tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, f)
    return "built"


# ------------------------------------------------------------------ stage 2: score
#: The four channels whose value depends on a REFERENCE SET rather than on the structure alone.
#: They must be referenced to the REAL 500-member pool, never to the pool-plus-ladder: a ladder-A
#: rung is surrounded by its own near-native siblings, so a set-referenced channel computed over
#: the union would score it well BECAUSE of the ladder.  That is circular, and at 47% ladder
#: contamination it is not a small effect.  These four are recomputed here against the pool only,
#: which is also the deployable definition (a selector scores a candidate against ITS pool).
POOL_REF = ("CONS", "DMAP_CONS", "TORS_CONS", "POOLGO")


def consistency_pool(W, PHI, PSI, D, sep, K):
    """CONS / DMAP_CONS / TORS_CONS / POOLGO with the reference taken from the first K rows (the
    real pool) only.  Identical to `ham_lib`'s definitions except for the reference index; a pool
    member excludes itself from its own CONS, exactly as `h_cons` does."""
    from s12 import instrument as I
    out = {}
    P = np.asarray(W[:K], float)
    r = np.empty((len(W), K))
    for a in range(len(W)):
        r[a] = I.kabsch_rmsd_batch(P, W[a])
    cons = r.sum(1).copy()
    cons[:K] -= np.diagonal(r[:K, :K])
    out["CONS"] = cons / np.where(np.arange(len(W)) < K, K - 1, K)
    med = np.median(D[:K], axis=0)
    out["DMAP_CONS"] = np.abs(D - med[None, :]).mean(1)

    def cmean(x):
        return np.arctan2(np.sin(x).mean(0), np.cos(x).mean(0))
    mphi, mpsi = cmean(PHI[:K]), cmean(PSI[:K])
    dphi = np.abs(np.angle(np.exp(1j * (PHI - mphi[None, :]))))
    dpsi = np.abs(np.angle(np.exp(1j * (PSI - mpsi[None, :]))))
    out["TORS_CONS"] = (dphi + dpsi).sum(1)
    m = sep >= 3
    C = (D[:, m] < 8.0).astype(float)
    f = C[:K].mean(0)
    out["POOLGO"] = -(C * f[None, :]).sum(1)
    return out


def score_all(cand_like, universe, dg, rcnt, K=None):
    """Every S27 channel on one candidate set.  Returns {name: (k,) float}."""
    from s27 import ham_lib as HL
    from s12 import instrument as I
    cx = HL.Context(cand_like, universe, dg, rcnt)
    ch = {}
    for name, fn in HL.CHANNELS.items():
        if K is not None and name in POOL_REF:
            continue                      # recomputed below against the pool only
        try:
            ch[name] = np.asarray(fn(cx), float)
        except Exception as e:                                   # pragma: no cover
            ch[name] = np.full(cx.k, np.nan)
            print(f"  channel {name} failed: {e}", file=sys.stderr)
    if K is not None:
        ch.update(consistency_pool(cx.W, cx.PHI, cx.PSI, cx.D, cx.sep, K))
    try:
        lt = HL.legacy_terms(cx)
        ch.update({k: np.asarray(v, float) for k, v in lt.items()})
        ch["LEG"] = np.sum([np.asarray(v, float) for v in lt.values()], axis=0)
    except Exception as e:                                       # pragma: no cover
        print(f"  legacy failed: {e}", file=sys.stderr)
    D = I.pair_dists(cx.W, cx.i, cx.j)
    ch["DIS"] = np.asarray(I.shipped_score(dg, D.astype(np.float32).astype(float)), float)
    return ch


def rescale_to(W, r0):
    """Uniformly rescale every structure so its radius of gyration is exactly r0.  A pure change
    of scale: shape, torsions and every dimensionless ratio are untouched."""
    W = np.asarray(W, float)
    c = W.mean(1, keepdims=True)
    g = rg_of(W)[:, None, None]
    return c + (W - c) * (r0 / np.maximum(g, 1e-9))


def score_si(cand_like, universe, dg, rcnt, K=None):
    """The SIZE-MATCHED twins: the universe AND the candidates rescaled to one common Rg, so the
    pair potentials are fitted and evaluated in reduced units.  Scale-invariant BY CONSTRUCTION --
    no partialling, nothing removed post hoc, nothing for a sceptic to argue about.

    This creates NO new information: it is a re-parameterisation of an existing structure->score
    map, so by the S29 section-8 bound it reaches the endpoint only as a cosine.  It is a
    correctness fix for the in-band question, NOT an Angstrom route."""
    from s27 import ham_lib as HL
    from types import SimpleNamespace
    r0 = float(np.median(rg_of(np.asarray(universe["W"], float))))
    u2 = {k: universe[k] for k in ("W", "S", "PHI", "PSI")}
    u2["W"] = rescale_to(np.asarray(universe["W"], float), r0)
    c2 = SimpleNamespace(seq=cand_like.seq, n=cand_like.n, fold=cand_like.fold,
                         k=cand_like.k, W=rescale_to(cand_like.W, r0),
                         PHI=cand_like.PHI, PSI=cand_like.PSI)
    cx = HL.Context(c2, u2, dg, rcnt)
    out = {}
    for name in SI_CHANNELS:
        if K is not None and name in POOL_REF:
            continue
        try:
            out[name + "_SI"] = np.asarray(HL.CHANNELS[name](cx), float)
        except Exception as e:                                   # pragma: no cover
            print(f"  SI channel {name} failed: {e}", file=sys.stderr)
    if K is not None:
        cp = consistency_pool(cx.W, cx.PHI, cx.PSI, cx.D, cx.sep, K)
        out.update({k + "_SI": v for k, v in cp.items() if k in SI_CHANNELS})
    return out, r0


# ------------------------------------------------------------------ stage 3: per-target row
def cv_r2(X, y, m_onehot, folds=5, alphas=(1.0, 10.0, 100.0, 1000.0), seed=0):
    """Held-out R2 of ridge(y ~ [m_onehot | X]) minus ridge(y ~ m_onehot).  Honest, not in-sample."""
    y = np.asarray(y, float)
    n = y.size
    rng = np.random.default_rng(seed)
    part = rng.permutation(n) % folds

    def _r2(Xb):
        pred = np.zeros(n)
        for f in range(folds):
            tr, te = part != f, part == f
            if tr.sum() < 5 or te.sum() < 1:
                return float("nan")
            A = Xb[tr]; mu = A.mean(0); sd = A.std(0) + 1e-9
            A = (A - mu) / sd
            b = y[tr] - y[tr].mean()
            best, bp = None, None
            for al in alphas:
                W = np.linalg.solve(A.T @ A + al * np.eye(A.shape[1]), A.T @ b)
                r = b - A @ W
                if best is None or float(r @ r) < best:
                    best, bp = float(r @ r), W
            pred[te] = ((Xb[te] - mu) / sd) @ bp + y[tr].mean()
        ss = float(((y - pred) ** 2).sum())
        return 1.0 - ss / max(float(((y - y.mean()) ** 2).sum()), 1e-12)

    base = _r2(m_onehot)
    full = _r2(np.column_stack([m_onehot, X]))
    return full, base, (full - base)


def target_row(pdb, seed=0):
    from core import geometry as geo
    from s12 import instrument as I
    from s24 import d_harness as H

    t0 = time.time()
    cand = H.Candidates.from_universe(pdb)
    n, fold = int(cand.n), int(cand.fold)
    nat = np.asarray(cand.nat_ca, float)
    nt = native_torsions(pdb)
    if nt is None:
        return {"pdb": pdb, "skipped": "no native torsions in peptide_db"}
    phi0, psi0, seq = nt
    assert seq == cand.seq, f"{pdb}: peptide_db seq {seq} != cand {cand.seq}"

    sampler = Sampler(fold, seq)
    rng = np.random.default_rng(abs(hash(("s30R", pdb, seed))) % (2 ** 32))
    grid = tuple(m for m in M_GRID if m <= n)

    # --- ladder A: anchored on the native's own torsions
    PA, SA, MA = make_ladder(phi0, psi0, sampler, n, rng, grid, R_DRAWS)
    # --- ladder B: anchored on a REAL pool member
    bi = int(rng.integers(cand.k))
    phiB = np.asarray(cand.PHI[bi], float); psiB = np.asarray(cand.PSI[bi], float)
    PB, SB, MB = make_ladder(phiB, psiB, sampler, n, rng, grid, R_DRAWS)
    # --- the anchors themselves, rebuilt
    PX = np.vstack([phi0[None], phiB[None]]); SX = np.vstack([psi0[None], psiB[None]])

    # --- fixed projected rungs
    zc = np.load(chain_path(pdb))
    PF = np.array([np.asarray(zc[f"phi_{k}"], float) for k in FIXED_RUNGS])
    SF = np.array([np.asarray(zc[f"psi_{k}"], float) for k in FIXED_RUNGS])

    # --- one candidate set, one reference for every pool-consistency channel
    PHI = np.vstack([np.asarray(cand.PHI, float), PA, PB, PX, PF])
    PSI = np.vstack([np.asarray(cand.PSI, float), SA, SB, SX, SF])
    bb = geo.build_backbone_batch(PHI, PSI)
    W = np.asarray(bb["CA"], float)
    # the pool's own rows keep their DEPOSITED coordinates; everything else is the ideal rebuild
    W[:cand.k] = np.asarray(cand.W, float)

    K = cand.k
    iA = np.arange(K, K + len(PA))
    iB = np.arange(iA[-1] + 1, iA[-1] + 1 + len(PB))
    iNATR, iANCHB = iB[-1] + 1, iB[-1] + 2
    iFIX = {k: iB[-1] + 3 + q for q, k in enumerate(FIXED_RUNGS)}

    # A plain namespace, never a Candidates: `nat_ca` and `oracle_rr` are ABSENT by construction,
    # so no channel can read a native even by accident (audit A3).
    from types import SimpleNamespace
    cl = SimpleNamespace(seq=seq, n=n, fold=fold, k=len(W), W=W, PHI=PHI, PSI=PSI)

    universe = I.load_univ(pdb)
    dg = I.distogram(pdb, seq, fold)
    ch = score_all(cl, universe, dg, rama_cnt()[fold], K=K)
    si, r0_si = score_si(cl, universe, dg, rama_cnt()[fold], K=K)
    ch.update(si)

    # ---------------- ORACLE labels, computed ONLY AFTER every channel value exists
    rmsd_nat = ca_rmsd_to(W, nat)
    rg = rg_of(W)
    reb_floor = float(rmsd_nat[iNATR])
    rmsd_anchorB = ca_rmsd_to(W[iB], W[iANCHB])

    row = dict(pdb=pdb, n=n, fold=fold, k_pool=K, n_ladderA=len(PA), n_ladderB=len(PB),
               m_grid=list(grid), R=R_DRAWS, seed=seed,
               rebuild_floor=reb_floor,
               anchorB_rmsd=float(rmsd_nat[iANCHB]),
               # the ladder's share of the candidate set.  It is reported, not hidden -- but the
               # four SET-REFERENCED channels are referenced to the pool ONLY (see POOL_REF), so
               # nothing that depends on a reference set is contaminated by it.
               contamination=float((len(PA) + len(PB) + 5) / len(W)),
               pool_ref_channels=list(POOL_REF),
               pool_best=float(np.min(rmsd_nat[:K])), pool_med=float(np.median(rmsd_nat[:K])),
               rmsd_fixed={k: float(rmsd_nat[i]) for k, i in iFIX.items()},
               secs=0.0, ch={})

    # --- D3: the anchor control's own confound, measured BEFORE the verdict
    cc = [spearman(rmsd_anchorB[MB == m], rmsd_nat[iB][MB == m]) for m in grid]
    row["anchor_confound_per_m"] = [None if not np.isfinite(v) else float(v) for v in cc]
    row["anchor_confound"] = float(np.nanmean(cc))

    # --- A2: realism flatness of ladder A within m (can only weaken my own positives)
    flat = {}
    for nm in ("RAMA", "EXVOL"):
        if nm in ch:
            flat[nm] = float(np.nanmean([spearman(ch[nm][iA][MA == m], rmsd_nat[iA][MA == m])
                                         for m in grid]))
    flat["RG_DEV"] = float(np.nanmean([spearman(np.abs(rg[iA] - np.median(rg[:K]))[MA == m],
                                                rmsd_nat[iA][MA == m]) for m in grid]))
    row["realism_flatness"] = flat

    # --- SI audit: the size-matched twin must make the two PURE functions of Rg constant.
    # This is the proof the transform does what it claims, not an assertion that it does.
    row["si_r0"] = float(r0_si)
    row["si_audit"] = {}
    for nm in ("RG_LAW", "RG_UNIV"):
        k2 = nm + "_SI"
        if k2 in ch:
            v = np.asarray(ch[k2], float)
            sc = float(np.std(v) / (abs(np.mean(v)) + 1e-12))
            row["si_audit"][nm] = dict(rel_sd=sc,
                                       raw_rel_sd=float(np.std(ch[nm]) / (abs(np.mean(ch[nm])) + 1e-12)))
    row["si_rho_rg"] = {}
    for nm in SI_CHANNELS:
        k2 = nm + "_SI"
        if k2 in ch and np.std(ch[k2]) > 0:
            row["si_rho_rg"][k2] = spearman(ch[k2][:K], rg[:K])
        if nm in ch and np.std(ch[nm]) > 0:
            row["si_rho_rg"][nm] = spearman(ch[nm][:K], rg[:K])

    # --- D1: the locality decomposition (held-out R2, m as the baseline block)
    yA = rmsd_nat[iA]
    onehot = np.column_stack([(MA == m).astype(float) for m in grid])
    Xloc = np.column_stack([np.sin(PA), np.cos(PA), np.sin(SA), np.cos(SA),
                            np.abs(np.angle(np.exp(1j * (PA - phi0[None])))),
                            np.abs(np.angle(np.exp(1j * (SA - psi0[None]))))])
    ii, jj = I.pair_index(n)
    Dmap = I.pair_dists(W[iA], ii, jj)
    Xglob = np.column_stack([Dmap, rg[iA][:, None], (Dmap < 8.0).sum(1)[:, None]])
    fl, bl, dl = cv_r2(Xloc, yA, onehot, seed=1)
    fg, bg, dg2 = cv_r2(Xglob, yA, onehot, seed=1)
    row["locality"] = dict(r2_m_only=bl, r2_local=fl, r2_global=fg,
                           d_local=dl, d_global=dg2,
                           n_feat_local=int(Xloc.shape[1]), n_feat_global=int(Xglob.shape[1]))

    # ---------------- per-channel statistics
    near = iA[(rmsd_nat[iA] <= NEAR_NATIVE_A)]
    row["n_near"] = int(near.size)
    if near.size == 0:                       # fall back to the single nearest rung, flagged
        near = iA[[int(np.argmin(rmsd_nat[iA]))]]
        row["near_fallback"] = True
    row["near_rmsd_mean"] = float(np.mean(rmsd_nat[near]))

    for name, v in ch.items():
        v = np.asarray(v, float)
        if not np.isfinite(v).all():
            fin = np.isfinite(v)
            if fin.sum() < 0.9 * v.size or np.std(v[fin]) == 0:
                continue
            v = np.where(fin, v, np.nanmax(v[fin]) + 1.0)
        if np.std(v) == 0:
            continue
        cell = {}
        # LANE L's MINIMUM: rho(score, Rg) travels beside EVERY recognition number, so a positive
        # is falsifiable as "recognition" and a null cannot be the S29-L19 self-fulfilling null.
        cell["rho_rg_pool"] = spearman(v[:K], rg[:K])
        cell["rho_rg_A"] = float(np.nanmean([spearman(v[iA][MA == m], rg[iA][MA == m])
                                             for m in grid]))
        # (1)(2)(3) within-m ordering, raw and partialled on both compactness ranks
        for tag, idx, ms, lab in (("A", iA, MA, rmsd_nat[iA]),
                                  ("B", iB, MB, rmsd_nat[iB]),
                                  ("ANCHOR", iB, MB, rmsd_anchorB)):
            raw, par = [], []
            for m in grid:
                sel = ms == m
                if sel.sum() < 6:
                    continue
                raw.append(spearman(v[idx][sel], lab[sel]))
                par.append(partial_spearman(v[idx][sel], lab[sel], [rg[idx][sel]], vshape_of={0}))
            cell["rho_" + tag] = float(np.nanmean(raw)) if raw else float("nan")
            cell["rho_" + tag + "_part"] = float(np.nanmean(par)) if par else float("nan")
        # D3's confound-immune form: order RMSD-to-native with rank(RMSD-to-anchor) removed
        pb = []
        for m in grid:
            sel = MB == m
            if sel.sum() >= 6:
                pb.append(partial_spearman(v[iB][sel], rmsd_nat[iB][sel],
                                           [rmsd_anchorB[sel], rg[iB][sel]], vshape_of={1}))
        cell["rho_B_part_anchor"] = float(np.nanmean(pb)) if pb else float("nan")
        # (4) PREF against the projected production chain + the S28-L36 pool-member control
        vp = float(v[iFIX["PROD"]])
        cell["pref_near"] = float(np.mean((v[near] < vp) + 0.5 * (v[near] == vp)))
        cell["pref_pool"] = float(np.mean((v[:K] < vp) + 0.5 * (v[:K] == vp)))
        cell["pref_oracle_circbest"] = float((v[iFIX["circ_best"]] < vp)
                                             + 0.5 * (v[iFIX["circ_best"]] == vp))
        cell["pref_native_rebuilt"] = float((v[iNATR] < vp) + 0.5 * (v[iNATR] == vp))
        # S3: the rebuilt native's percentile inside its own ladder A
        cell["pctile_nat_in_A"] = float(np.mean(v[iA] < v[iNATR]))
        # D2: resolution -- pairwise concordance by |delta RMSD|, within m
        conc = np.zeros(len(DELTA_LABELS)); cnt = np.zeros(len(DELTA_LABELS))
        for m in grid:
            sel = np.where(MA == m)[0]
            if sel.size < 4:
                continue
            vv, ll = v[iA][sel], rmsd_nat[iA][sel]
            a, b = np.triu_indices(sel.size, 1)
            d = np.abs(ll[a] - ll[b])
            ok = ((vv[a] < vv[b]) == (ll[a] < ll[b])).astype(float)
            ok[vv[a] == vv[b]] = 0.5
            bi2 = np.clip(np.digitize(d, DELTA_BINS) - 1, 0, len(DELTA_LABELS) - 1)
            for q in range(len(DELTA_LABELS)):
                s2 = bi2 == q
                if s2.any():
                    conc[q] += ok[s2].sum(); cnt[q] += s2.sum()
        cell["conc"] = [float(conc[q] / cnt[q]) if cnt[q] > 0 else None for q in range(len(DELTA_LABELS))]
        cell["conc_n"] = [int(cnt[q]) for q in range(len(DELTA_LABELS))]
        # --- features for the LEAVE-FOLD-OUT combination (prereg item 6).  Rank-z over the union
        # of ladder A, the pool and the fixed rungs, so every target contributes on one scale.
        from s27.ham_lib import zrank
        sub = np.concatenate([iA, np.arange(K), [iFIX["PROD"]]])
        z = zrank(v[sub])
        zA, zP, zprod = z[:len(iA)], z[len(iA):len(iA) + K], float(z[-1])
        nearmask = rmsd_nat[iA] <= NEAR_NATIVE_A
        farmask = rmsd_nat[iA] >= 3.0
        cell["d_near"] = float(np.mean(zA[nearmask]) - zprod) if nearmask.any() else None
        cell["d_far"] = float(np.mean(zA[farmask]) - zprod) if farmask.any() else None
        cell["d_pool"] = float(np.mean(zP) - zprod)
        row["ch"][name] = cell

    row["secs"] = time.time() - t0
    return row


# ------------------------------------------------------------------ drivers
def _pdbs(limit=0, shard=0, of=1):
    from s12 import instrument as I
    ps = [t["pdb"] for t in I.targets()]
    if limit:
        ps = ps[:limit]
    return [p for q, p in enumerate(ps) if q % of == shard]


def cmd_prep(a):
    ps = _pdbs(a.limit, a.shard, a.of)
    t0 = time.time()
    for q, p in enumerate(ps):
        st = prep_target(p, rebuild=a.rebuild)
        if q % 5 == 0 or st == "built":
            print(f"[prep {a.shard}/{a.of}] {q+1}/{len(ps)} {p} {st} ({time.time()-t0:.0f}s)",
                  flush=True)


def cmd_run(a):
    ps = _pdbs(a.limit, a.shard, a.of)
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, f"s30_R_rows.s{a.shard}of{a.of}.jsonl")
    done = set()
    if os.path.exists(out) and not a.rebuild:
        with open(out, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    done.add(json.loads(line)["pdb"])
    t0 = time.time()
    with open(out, "a", encoding="utf-8") as fh:
        for q, p in enumerate(ps):
            if p in done:
                continue
            try:
                r = target_row(p, seed=a.seed)
            except Exception as e:
                import traceback
                traceback.print_exc()
                r = {"pdb": p, "error": repr(e)}
            fh.write(json.dumps(r) + "\n"); fh.flush()
            print(f"[run {a.shard}/{a.of}] {q+1}/{len(ps)} {p} "
                  f"{r.get('secs', 0):.1f}s total {time.time()-t0:.0f}s", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("prep", "run", "both"):
        s = sub.add_parser(nm)
        s.add_argument("--limit", type=int, default=0)
        s.add_argument("--shard", type=int, default=0)
        s.add_argument("--of", type=int, default=1)
        s.add_argument("--rebuild", action="store_true")
        s.add_argument("--seed", type=int, default=0)
    sub.add_parser("agg")
    a = ap.parse_args(argv)
    if a.cmd == "prep":
        cmd_prep(a)
    elif a.cmd == "run":
        cmd_run(a)
    elif a.cmd == "both":
        cmd_prep(a)
        cmd_run(a)
    else:
        from s30 import s30_R_agg
        s30_R_agg.main([])


if __name__ == "__main__":
    main()
