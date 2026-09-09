"""S24 / LANE D / D1 -- HAMILTONIAN-DISAGREEMENT SAMPLING. Directive section 28.

PRE-REGISTERED IN `s24/PREREG_D.md` (D1) BEFORE THIS FILE EXISTED. Fork list sent to the
coordinator before the run. Read that file for the six axes; the summary is in the table below.

WHY THIS IS NOW THE HIGHEST-VALUE ITEM IN THE LANE
==================================================
Two results closed the two obvious levers in the same week:

  * `d_setequality_proof.py` (D0): the CVaR tail's membership is a SUBSET of the classical
    energy-ordered prefix and equals it under full support. **The SELECTION RULE is fixed by
    the algorithm.** No candidate manifold can change it.
  * s24 L3 (coordinator): a retrieval-free, quality-matched source selected by the SAME
    shipped distogram score has bias cosine **+0.9432** against the incumbent -- ABOVE the
    +0.9330 within-source control -- with a flat mixture curve. **PROVENANCE is worth nothing
    once the score does the selecting.**

Selection rule fixed, provenance worthless. What is left is **the functional**, and Legacy and
AMBER are the only two genuinely different, genuinely evaluable functionals this project has.
Hence D1.

THE QUESTION, AND MY REGISTERED PRIOR THAT IT FAILS
===================================================
Does partitioning candidates by Legacy/AMBER DISAGREEMENT isolate subsets whose emitted average
points in a materially different direction from the incumbent's?

I registered, in advance, that I expect NO. Both energies are full-register in torsion space
(`torsion-space-locality-theorem`; "AMBER is less local" is a category error), both are
dominated by compactness and sterics on real protein windows, and L2/L3 have shown that
anything selected to be *good* ends up parallel. Predicted cosines 0.85-0.95: partitions that
differ in QUALITY but not in DIRECTION. If that is what comes back it is the result and it
closes the functional lever the same way L3 closed the provenance lever.

THE SHARED-REFERENT FLOOR, MEASURED FIRST
=========================================
Every bias vector here is `avg - native` after Kabsch, so all of them are measured against a
COMMON referent and therefore correlate BY CONSTRUCTION. Project memory
(`shared-referent-floor`) records that ignoring this turned a "2/3 sequence-independent" claim
into 1/5. So the floor is measured before anything is interpreted: the cosine that two
arbitrary same-size subsets OF THE SAME POOL show against the incumbent. **A disagreement
partition has to beat that floor, not zero.**

GENUINENESS, per the three pillars
==================================
`E_Legacy` = `s16.energy_lib.legacy_components_of_windows` + `legacy_total_from()` at
`DEFAULT_WEIGHTS`, never fitted. `E_AMBER` = genuine ff14SB/GBn2 SINGLE POINT through
`s20.qb2_lib.AmberSP`, which asserts bit-exactness against
`core.amber.refine_coords(k_restraint=0, steps=-1)` on every target before any number is read.

**NO REFINEMENT ANYWHERE.** s23 L11 measured 17 of 17 restrained/unrestrained repair settings
at or worse than no repair, with the ladder bottoming out at the no-op. AMBER is an ENERGY
MEASUREMENT in this module and nothing else: no minimisation, no relaxation, no `AmberRelax`.

SERIALISATION. The OpenMM lane is Lane D's and is serialised under `s24/results/LOCK_AMBER`
(`d_harness.AmberLock`, `os.open(O_CREAT|O_EXCL)`), announced on open and on release, held for
bounded periods.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from s24 import d_harness as H            # noqa: E402

SALT = "s24laneD_D1"
RESULTS = H.RESULTS
M_INC = 75                                # the incumbent's rung; the referent for every cosine
N_FLOOR_DRAWS = 24                        # random same-size subsets for the shared-referent floor


# ============================================================ bias geometry (qmatch's recipe)
# Reproduced from `s24/qmatch.py` VERBATIM so L2/L3's cosines and mine are the same object.
def _kabsch_R(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _bias(C, nat):
    """The emitted average's error vector IN THE NATIVE FRAME -- the only frame in which two
    sources' errors are commensurable (L2's own justification)."""
    R = _kabsch_R(C, nat)
    return (C - C.mean(0)) @ R.T - (nat - nat.mean(0))


def _cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


def _z_raw(x) -> np.ndarray:
    """Moment z-score. DEGENERATE on raw AMBER -- kept as the declared SECONDARY arm."""
    x = np.asarray(x, float)
    return (x - x.mean()) / max(float(x.std()), 1e-12)


def _z(x) -> np.ndarray:
    """RANK standardisation -- the PRIMARY, per PREREG D1-A.

    Genuine ff14SB/GBn2 single points on UNRELAXED retrieval windows are clash-dominated: a
    single 1e18 clash sets the standard deviation and 98.3-99.7% of the pool lands inside
    |z_raw| < 0.1, so the moment z-score degenerates into "which candidate has the worst
    steric clash". Project memory: `pauli-spectrum-delta-spike-artefact` -- condition
    MONOTONICALLY, and 99th-percentile winsorisation is not enough.

    Ranking is strictly monotone, so it changes no ordering, no argmin and no level set and
    cannot manufacture a preference. It is also the currency `core.pipeline._zrank` already
    deploys. Under it 6.0% of candidates sit inside |z| < 0.1, the well-conditioned value.
    """
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(float(r.std()), 1e-12)


# ============================================ alignment with the DISTOGRAM's own error (D1-B)
# VERBATIM from `s24/referent.py::_stats` -- the coordinator's own L5 implementation, reused
# rather than reimplemented so the two numbers are the same object, not nearly the same one.
def _dhat(dg) -> np.ndarray:
    """The distogram's MAP distance prediction: grid value at the argmin of the risk."""
    return np.asarray(dg["grid"], float)[np.argmin(np.asarray(dg["risk"], float), axis=1)]


def _distogram_alignment(Dc: np.ndarray, Dt: np.ndarray, Dh: np.ndarray) -> Dict:
    """`cos` and `beta` of the emitted distance error against the PRIOR's distance error.

    eP = Dh - Dt is the prior's own error; eC = Dc - Dt is the emitted structure's. beta is
    the share of the prior's error that reappears in the output (1 = total inheritance).
    Neither vector is centred, matching referent.py exactly.
    """
    eP, eC = np.asarray(Dh, float) - Dt, np.asarray(Dc, float) - Dt
    den = float((eP * eP).sum())
    nC, nP = float(np.linalg.norm(eC)), float(np.linalg.norm(eP))
    return dict(cos_distogram=float((eC * eP).sum() / (nC * nP)) if nC > 0 and nP > 0
                else float("nan"),
                beta=float((eC * eP).sum() / den) if den > 0 else float("nan"),
                rms_cloud_to_hat=float(np.sqrt(np.mean((Dc - Dh) ** 2))),
                rms_nat_to_hat=float(np.sqrt(np.mean((Dt - Dh) ** 2))))


# ============================================================ the genuine energies
def energy_amber(cand: H.Candidates, sp=None) -> np.ndarray:
    """GENUINE ff14SB/GBn2 SINGLE POINT per candidate, kcal/mol. No minimisation.

    `AmberSP.verify` is called before any number is read -- a gate that fires zero times is
    not evidence, so it is called with the pool's OWN torsions rather than a synthetic probe.
    """
    if not cand.has_torsions:
        raise ValueError("AMBER needs PHI/PSI")
    from s20 import qb2_lib as QB2
    from s13 import qarch_lib as QA
    if sp is None:
        rep = QA.Space(cand.pdb, 4).rep
        sp = QB2.AmberSP(cand.seq, rep)
    v = sp.verify(cand.PHI[:4], cand.PSI[:4], m=4)
    return sp.batch(cand.PHI, cand.PSI), sp, v


# ============================================================ the pre-declared partition
def partition(zL: np.ndarray, zA: np.ndarray) -> Dict[str, np.ndarray]:
    """The partition rule fixed a priori in `s24/PREREG_D.md` D1. Nothing here is tuned.

    `d = zL - zA`; positive means AMBER likes the candidate relatively more.
    """
    d = zL - zA
    ad = np.abs(d)
    q25, q75 = np.quantile(d, 0.25), np.quantile(d, 0.75)
    q75a = np.quantile(ad, 0.75)
    return {
        "AGREE_GOOD": np.flatnonzero((zL < -0.5) & (zA < -0.5)),
        "LEGACY_PREFERS": np.flatnonzero(d < q25),
        "AMBER_PREFERS": np.flatnonzero(d > q75),
        "STRONG_DISAGREE": np.flatnonzero(ad > q75a),
        "AGREE_BAD": np.flatnonzero((zL > 0.5) & (zA > 0.5)),
    }


# ============================================================ one target
def run_target(pdb: str, k: int = 500, m_inc: int = M_INC,
               n_floor: int = N_FLOOR_DRAWS, sp=None) -> Dict:
    cand = H.Candidates.from_universe(pdb, k=k)
    nat = cand.nat_ca
    sc_dist = H.score_shipped(cand)
    sc_leg = H.score_legacy(cand)
    e_amb, sp, ver = energy_amber(cand, sp=sp)

    # --- the distogram's MAP prediction and the native's pair distances, for D1-B
    dg = I.distogram(pdb, cand.seq, cand.fold)
    pi, pj = I.pair_index(cand.n)                     # min_sep=2, the instrument default
    Dh = _dhat(dg)
    Dt = I.pair_dists(nat, pi, pj)

    # --- THE REFERENT: the incumbent's own top-75 by the shipped score
    inc_idx = np.argsort(sc_dist, kind="stable")[:m_inc]
    C_inc = _avg(cand.W[inc_idx])
    e_inc = _bias(C_inc, nat)
    rmsd_inc = float(I.ca_rmsd(C_inc, nat))
    align_inc = _distogram_alignment(I.pair_dists(C_inc, pi, pj), Dt, Dh)

    zL, zA = _z(sc_leg), _z(e_amb)                    # PRIMARY: rank standardisation (D1-A)
    parts = partition(zL, zA)
    # declared SECONDARY: the degenerate moment z-score, run rather than quietly dropped
    parts_raw = partition(_z_raw(sc_leg), _z_raw(e_amb))
    overlap_raw = {nm: (float(len(np.intersect1d(parts[nm], parts_raw[nm]))
                              / max(1, len(parts[nm]))) if len(parts[nm]) else None)
                   for nm in parts}
    frac_degenerate = float(np.mean(np.abs(_z_raw(e_amb)) < 0.1))
    frac_clash = float(np.mean(np.asarray(e_amb, float) > 1e4))

    from scipy.stats import spearmanr
    rho_LA = float(spearmanr(sc_leg, e_amb).statistic)
    rho_Ld = float(spearmanr(sc_leg, sc_dist).statistic)
    rho_Ad = float(spearmanr(e_amb, sc_dist).statistic)

    out_parts = {}
    for name, idx in parts.items():
        if len(idx) < 3:
            out_parts[name] = dict(size=int(len(idx)), rmsd=None, cos_vs_incumbent=None,
                                   q=None, capture_of_top75=None)
            continue
        C = _avg(cand.W[idx])
        e = _bias(C, nat)
        al = _distogram_alignment(I.pair_dists(C, pi, pj), Dt, Dh)
        # --- D1-C, the QUALITY-MATCHED arm. The raw partition is not score-selected, so
        # comparing it to the incumbent conflates "a different functional" with "no selection
        # at all" -- and the raw partitions are 1.4-2.4 A worse, so that confound is large.
        # This arm applies the SHIPPED SCORE *inside* the partition and takes the same rung,
        # which is how a physics functional would actually be deployed: as a PRE-FILTER ahead
        # of the deployed scorer. Score-selected and size-matched by construction.
        sub = idx[np.argsort(sc_dist[idx], kind="stable")[:min(m_inc, len(idx))]]
        Cs = _avg(cand.W[sub])
        es = _bias(Cs, nat)
        als = _distogram_alignment(I.pair_dists(Cs, pi, pj), Dt, Dh)
        out_parts[name] = dict(
            filt_size=int(len(sub)),
            filt_rmsd=float(I.ca_rmsd(Cs, nat)),
            filt_cos_vs_incumbent=_cos(es, e_inc),
            filt_cos_distogram=als["cos_distogram"], filt_beta=als["beta"],
            filt_capture_of_top75=float(len(np.intersect1d(sub, inc_idx))
                                        / max(1, len(inc_idx))),
            size=int(len(idx)),
            rmsd=float(I.ca_rmsd(C, nat)),
            cos_vs_incumbent=_cos(e, e_inc),
            q=float(np.linalg.norm(e) / max(np.linalg.norm(e_inc), 1e-12)),
            capture_of_top75=float(len(np.intersect1d(idx, inc_idx)) / max(1, len(inc_idx))),
            oracle_best=float(np.min(cand.oracle_rr[idx])),
            overlap_with_raw_z=overlap_raw[name],
            **al,
        )

    # --- THE SHARED-REFERENT FLOOR, measured BEFORE anything is interpreted.
    # Random same-size subsets of the SAME pool, through the identical readout. Their cosine
    # against the incumbent is the number a "different direction" claim must beat.
    rng = SD.stable_rng(pdb, "floor", int(m_inc), salt=SALT)
    floor = {}
    for name, idx in parts.items():
        if len(idx) < 3:
            floor[name] = None
            continue
        cs, ds = [], []
        for _ in range(int(n_floor)):
            pick = rng.choice(cand.k, size=len(idx), replace=False)
            Cr = _avg(cand.W[pick])
            cs.append(_cos(_bias(Cr, nat), e_inc))
            # the floor for D1-B too: a RANDOM subset's alignment with the distogram's error.
            # If a random subset already sits at ~0.640, the partition's 0.640 means nothing.
            ds.append(_distogram_alignment(I.pair_dists(Cr, pi, pj), Dt, Dh)["cos_distogram"])
        floor[name] = dict(mean=float(np.mean(cs)), sd=float(np.std(cs, ddof=1)),
                           lo=float(np.percentile(cs, 2.5)),
                           hi=float(np.percentile(cs, 97.5)),
                           cos_distogram_mean=float(np.mean(ds)),
                           cos_distogram_sd=float(np.std(ds, ddof=1)))

    return dict(
        pdb=pdb, n=cand.n, fold=cand.fold, k=cand.k, basis="point_cloud", label="ACHIEVABLE",
        m_inc=int(m_inc), rmsd_incumbent=rmsd_inc,
        incumbent_cos_distogram=align_inc["cos_distogram"], incumbent_beta=align_inc["beta"],
        rho_legacy_amber=rho_LA, rho_legacy_dist=rho_Ld, rho_amber_dist=rho_Ad,
        amber_mean=float(np.mean(e_amb)), amber_sd=float(np.std(e_amb)),
        legacy_mean=float(np.mean(sc_leg)), legacy_sd=float(np.std(sc_leg)),
        frac_amber_z_degenerate=frac_degenerate, frac_amber_clash_gt_1e4=frac_clash,
        amber_verify_max_rel=float(ver[0]) if isinstance(ver, (tuple, list)) else float(ver),
        partitions=out_parts, shared_referent_floor=floor,
        oracle_pool_best=float(np.min(cand.oracle_rr)),
    ), sp


ROW_KEYS = ("pdb", "n", "fold", "k", "basis", "label", "rmsd_incumbent",
            "rho_legacy_amber", "partitions", "shared_referent_floor")


# ============================================================ the panel
def declared_subset(n_targets: int, seed: int = 0) -> List[str]:
    """The target subset, DRAWN BEFORE THE RUN and stratified across the 5 pinned folds.

    Declared rather than chosen: the draw is a pure function of (n_targets, seed) through
    `s15.seed.stable_rng`, so it cannot be re-rolled after seeing a result.
    """
    tg = I.targets()
    by_fold: Dict[int, List[str]] = {}
    for t in tg:
        by_fold.setdefault(int(t["fold"]), []).append(t["pdb"])
    rng = SD.stable_rng("s24", "D1subset", int(n_targets), int(seed), salt=SALT)
    per = int(math.ceil(n_targets / max(1, len(by_fold))))
    out: List[str] = []
    for f in sorted(by_fold):
        pool = sorted(by_fold[f])
        out += list(rng.choice(pool, size=min(per, len(pool)), replace=False))
    return sorted(out)[:n_targets]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-targets", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--k", type=int, default=500)
    ap.add_argument("--out", type=str, default="d1_hamdisagree")
    ap.add_argument("--declare-only", action="store_true",
                    help="print the pre-declared subset and exit, before any AMBER call")
    a = ap.parse_args()

    subset = declared_subset(a.n_targets, a.seed)
    print(f"DECLARED SUBSET n={len(subset)} seed={a.seed}: {' '.join(subset)}", flush=True)
    if a.declare_only:
        H.write(a.out + "_subset", dict(subset=subset, n=len(subset), seed=a.seed,
                                        rows=[]), required_keys=None)
        return 0

    rows, t0 = [], time.time()
    # AMBER is serialised. Announce on open, announce on release, hold for a bounded period.
    with H.AmberLock(tag=f"D1 n={len(subset)} k={a.k}"):
        for ii, pdb in enumerate(subset):
            try:
                r, _sp = run_target(pdb, k=a.k)
                rows.append(r)
            except Exception as e:
                print(f"  !! {pdb}: {type(e).__name__}: {e}", flush=True)
                rows.append(dict(pdb=pdb, error=f"{type(e).__name__}: {e}"))
            el = time.time() - t0
            print(f"  [{ii+1}/{len(subset)}] {pdb}  {el:.0f}s "
                  f"({el/(ii+1):.1f}s/target)", flush=True)
            H.write(a.out, dict(rows=rows, complete=False, n_expected=len(subset),
                                subset=subset, seed=a.seed, k=a.k), required_keys=None)

    ok = [r for r in rows if "error" not in r]
    H.write(a.out, dict(rows=rows, n_expected=len(subset), subset=subset, seed=a.seed,
                        k=a.k, n_ok=len(ok),
                        complete=bool(len(ok) == len(subset))), required_keys=ROW_KEYS)
    report(rows)
    return 0


def report(rows: List[Dict]) -> None:
    ok = [r for r in rows if "error" not in r]
    if not ok:
        print("no successful rows")
        return
    print(f"\n{'='*78}\nD1 -- HAMILTONIAN DISAGREEMENT, n={len(ok)} targets, POINT-CLOUD basis")
    print(f"{'='*78}")
    print(f"incumbent top-{ok[0]['m_inc']} mean RMSD  {np.mean([r['rmsd_incumbent'] for r in ok]):.4f}")
    print(f"rho(Legacy, AMBER) within target: mean {np.mean([r['rho_legacy_amber'] for r in ok]):+.4f}"
          f"  median {np.median([r['rho_legacy_amber'] for r in ok]):+.4f}")
    print(f"rho(Legacy, distogram)          : mean {np.mean([r['rho_legacy_dist'] for r in ok]):+.4f}")
    print(f"rho(AMBER,  distogram)          : mean {np.mean([r['rho_amber_dist'] for r in ok]):+.4f}")
    print(f"AMBER conditioning: frac E>1e4 kcal {np.mean([r['frac_amber_clash_gt_1e4'] for r in ok]):.3f}"
          f" | frac |z_raw|<0.1 {np.mean([r['frac_amber_z_degenerate'] for r in ok]):.3f}"
          f"  (PREREG D1-A: why rank standardisation is PRIMARY)")
    names = list(ok[0]["partitions"].keys())
    print(f"\n--- PRIMARY 1: bias direction vs the incumbent -------------------------------")
    print(f"{'partition':<18}{'size':>6}{'RMSD':>9}{'cos_vs_inc':>12}{'FLOOR':>9}"
          f"{'excess':>9}{'q':>7}{'capt75':>8}")
    print("-" * 78)
    for nm in names:
        s = [r["partitions"][nm] for r in ok if r["partitions"][nm]["rmsd"] is not None]
        f = [r["shared_referent_floor"][nm] for r in ok
             if r["shared_referent_floor"].get(nm)]
        if not s:
            print(f"{nm:<18}{'--':>6}")
            continue
        c = float(np.mean([x["cos_vs_incumbent"] for x in s]))
        fl = float(np.mean([x["mean"] for x in f])) if f else float("nan")
        print(f"{nm:<18}{np.mean([x['size'] for x in s]):>6.0f}"
              f"{np.mean([x['rmsd'] for x in s]):>9.4f}{c:>12.4f}{fl:>9.4f}"
              f"{c - fl:>+9.4f}{np.mean([x['q'] for x in s]):>7.3f}"
              f"{np.mean([x['capture_of_top75'] for x in s]):>8.3f}")
    print("-" * 78)
    print("L2 bars: within-source control +0.9330 | a source must reach cos <= ~0.65 AND")
    print("         <= ~3.9 A standalone to BREAK EVEN. 'excess' = cos - shared-referent floor;")
    print("         a partition carries new DIRECTION only if it is materially BELOW its floor.")

    # ---- PRIMARY 2 (D1-B): does a different FUNCTIONAL escape the distogram's referent?
    print(f"\n--- PRIMARY 2 (D1-B): alignment with the DISTOGRAM'S OWN error ---------------")
    ci = float(np.mean([r["incumbent_cos_distogram"] for r in ok]))
    bi = float(np.mean([r["incumbent_beta"] for r in ok]))
    print(f"{'INCUMBENT top-75':<18}{'':>6}{'':>9}{ci:>12.4f}{'beta':>9}{bi:>9.4f}")
    print(f"{'partition':<18}{'size':>6}{'':>9}{'cos_distog':>12}{'FLOOR':>9}"
          f"{'excess':>9}{'beta':>9}")
    print("-" * 78)
    for nm in names:
        s = [r["partitions"][nm] for r in ok if r["partitions"][nm]["rmsd"] is not None]
        f = [r["shared_referent_floor"][nm] for r in ok
             if r["shared_referent_floor"].get(nm)]
        if not s:
            continue
        c = float(np.mean([x["cos_distogram"] for x in s]))
        fl = float(np.mean([x["cos_distogram_mean"] for x in f])) if f else float("nan")
        print(f"{nm:<18}{np.mean([x['size'] for x in s]):>6.0f}{'':>9}{c:>12.4f}"
              f"{fl:>9.4f}{c - fl:>+9.4f}{np.mean([x['beta'] for x in s]):>9.4f}")
    print("-" * 78)
    print("READING RULE (fixed in advance): if LEGACY_PREFERS and AMBER_PREFERS both sit at the")
    print("incumbent's alignment, both physics functionals inherit the SAME referent and the")
    print("FUNCTIONAL lever closes alongside the provenance lever. A partition that breaks away")
    print("is the sprint's most valuable number. Each floor is a RANDOM same-size subset.")

    # ---- D1-C: the physics functional used as a PRE-FILTER ahead of the deployed scorer.
    # Quality-matched and score-selected, so it removes the "not selected at all" confound.
    print(f"\n--- D1-C: physics as a PRE-FILTER, then the SHIPPED score, same rung ---------")
    print(f"{'INCUMBENT (no filter)':<22}{'':>6}{np.mean([r['rmsd_incumbent'] for r in ok]):>9.4f}"
          f"{1.0:>12.4f}{ci:>12.4f}")
    print(f"{'pre-filter':<22}{'m':>6}{'RMSD':>9}{'cos_vs_inc':>12}{'cos_distog':>12}"
          f"{'capt75':>8}")
    print("-" * 78)
    for nm in names:
        s = [r["partitions"][nm] for r in ok
             if r["partitions"][nm].get("filt_rmsd") is not None]
        if not s:
            continue
        print(f"{nm:<22}{np.mean([x['filt_size'] for x in s]):>6.0f}"
              f"{np.mean([x['filt_rmsd'] for x in s]):>9.4f}"
              f"{np.mean([x['filt_cos_vs_incumbent'] for x in s]):>12.4f}"
              f"{np.mean([x['filt_cos_distogram'] for x in s]):>12.4f}"
              f"{np.mean([x['filt_capture_of_top75'] for x in s]):>8.3f}")
    print("-" * 78)
    print("This is the arm that matters for deployment: it asks whether a physics functional")
    print("used as a FILTER moves the answer, holding the scorer, the rung and the readout fixed.")


if __name__ == "__main__":
    sys.exit(main())
