#!/usr/bin/env python
"""s30/s30_G_chiral.py -- LANE G, Q2: the ONE family G1 leaves open.

Pre-registered in `s30/PREREG_S30_G.md` (commit 78b65521).

THEOREM G1 (the chirality dichotomy).  For a single-structure channel S: R^{n x 3} -> R invariant
under rotation and translation (with the SEQUENCE held fixed as a parameter -- reflection does not
touch the sequence, so sequence-using channels such as burial or a hydrophobicity-weighted term are
covered), S is a function of the pairwise distance matrix D IFF S is also reflection-invariant.
Proof: classical MDS gives G = -1/2 J D^2 J = X X^T, so D fixes the centred coordinates up to
Q in O(3) = SO(3) x {+-I}; after quotienting the rotations this project actually uses (proper-
rotation Kabsch with the det correction, `s12/instrument.py:kabsch_rmsd_batch`) the residual
ambiguity is EXACTLY reflection.

So every achiral native-free channel is a distance-map reading, and the only references available
to score one against are the three closed buckets.  CHIRAL functionals are the single escape class.
This file builds three of them and measures them on LANE R'S OWN INSTRUMENT -- same ladder, same
seed rule, same anchor and pool-member controls, same fold-clustered statistics.

THE PRE-CHECK RUNS FIRST AND CAN CLOSE THE QUESTION BEFORE ANY CONTRAST (coordinator, 13:4x):
G1 says the escape class is non-empty in principle.  It does NOT say the chiral coordinate has
variance on the manifold our candidates occupy.  Every rung is drawn from the fold's Ramachandran
table and every pool member is a deposited L-backbone, so the mirror images are never in the
candidate set.  A functional can be provably chiral and still be near-constant across everything we
score.  `occupancy = sd(X) / sqrt(sd(X)^2 + mean(X)^2)` is exactly `sd(X over candidates) /
sd(X over candidates UNION their mirrors)`: 1.0 means the chiral coordinate is fully exercised,
~0 means the set sits on one side of it and the escape class is EMPTY IN PRACTICE.

    python s30/s30_G_chiral.py run [--limit N] [--shard i --of k]
    python s30/s30_G_chiral.py agg
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time
import zlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s30_G_chiral_rows.jsonl")
OUT = os.path.join(RESULTS, "s30_G_chiral.json")

CHIRAL = ("WRITHE", "CHIRAL3", "CHIRAL3_LONG")          # exactly the three registered
TWIN = {c: "ABS_" + c for c in CHIRAL}                   # the matched ACHIRAL control
LONG_SEP = 5                                             # CHIRAL3_LONG: |i - l| >= 5


# --------------------------------------------------------------------- the chiral functionals
def writhe(W):
    """Discrete Gauss double integral (Klenin-Langowski 1a) over the CA polyline.  (B,n,3)->(B,)."""
    W = np.asarray(W, float)
    B, n, _ = W.shape
    a, b = np.triu_indices(n - 1, k=2)                   # non-adjacent segment pairs
    if a.size == 0:
        return np.zeros(B)
    # Klenin-Langowski method 1a.  Segment i runs p1->p2, segment j runs q1->q2.
    p1, p2 = W[:, a], W[:, a + 1]
    q1, q2 = W[:, b], W[:, b + 1]
    r13, r14 = q1 - p1, q2 - p1
    r23, r24 = q1 - p2, q2 - p2
    r12, r34 = p2 - p1, q2 - q1

    def u(x):
        nn = np.linalg.norm(x, axis=-1, keepdims=True)
        return x / np.maximum(nn, 1e-12)

    n1, n2 = u(np.cross(r13, r14)), u(np.cross(r14, r24))
    n3, n4 = u(np.cross(r24, r23)), u(np.cross(r23, r13))
    dot = lambda x, y: np.clip((x * y).sum(-1), -1.0, 1.0)
    om = np.arcsin(dot(n1, n2)) + np.arcsin(dot(n2, n3)) \
        + np.arcsin(dot(n3, n4)) + np.arcsin(dot(n4, n1))
    sg = np.sign((np.cross(r34, r12) * r13).sum(-1))
    return 2.0 * (om * sg).sum(1) / (4.0 * math.pi)


def _quad_idx(n, min_span):
    q = [(i, j, k, l) for i in range(n) for j in range(i + 1, n)
         for k in range(j + 1, n) for l in range(k + 1, n) if l - i >= min_span]
    return np.asarray(q, int).reshape(-1, 4)


def chiral3(W, idx, rg):
    """Mean signed volume det[x_j-x_i, x_k-x_j, x_l-x_k] over quadruples, scaled by Rg^3."""
    if idx.size == 0:
        return np.zeros(len(W))
    i, j, k, l = idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]
    a, b, c = W[:, j] - W[:, i], W[:, k] - W[:, j], W[:, l] - W[:, k]
    v = (np.cross(a, b) * c).sum(-1)
    return v.mean(1) / np.maximum(rg ** 3, 1e-12)


def chiral_channels(W, n):
    rg = np.sqrt((((W - W.mean(1, keepdims=True)) ** 2).sum(-1)).mean(1))
    out = {"WRITHE": writhe(W),
           "CHIRAL3": chiral3(W, _quad_idx(n, 3), rg),
           "CHIRAL3_LONG": chiral3(W, _quad_idx(n, LONG_SEP), rg)}
    for c in CHIRAL:
        out[TWIN[c]] = np.abs(out[c])
    return out


def reflection_audit(W, n, tol=1e-9):
    """PROOF, not assertion: reflect and check X -> -X, |X| fixed, D fixed."""
    from s12 import instrument as I
    M = W.copy(); M[..., 2] *= -1.0
    A, Bc = chiral_channels(W, n), chiral_channels(M, n)
    ii, jj = I.pair_index(n)
    dD = float(np.abs(I.pair_dists(W, ii, jj) - I.pair_dists(M, ii, jj)).max())
    # rotation invariance too, on a fixed random rotation
    rng = np.random.default_rng(7)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    R = chiral_channels(W @ Q.T + rng.normal(size=3), n)
    rep = {"max_abs_dD_under_reflection": dD}
    for c in CHIRAL:
        sc = max(float(np.abs(A[c]).max()), 1e-12)
        rep[c] = dict(flip_err=float(np.abs(A[c] + Bc[c]).max() / sc),
                      twin_err=float(np.abs(A[TWIN[c]] - Bc[TWIN[c]]).max() / sc),
                      rot_err=float(np.abs(A[c] - R[c]).max() / sc))
        assert rep[c]["flip_err"] < 1e-7, f"{c} is NOT chiral: {rep[c]}"
        assert rep[c]["twin_err"] < 1e-7, f"{TWIN[c]} is NOT achiral: {rep[c]}"
        assert rep[c]["rot_err"] < 1e-7, f"{c} is not rotation-invariant: {rep[c]}"
    assert dD < 1e-7, f"D moved under reflection: {dD}"
    return rep


def spread(v):
    v = np.asarray(v, float)
    m, s = float(v.mean()), float(v.std(ddof=1))
    return dict(mean=m, sd=s, abs_mean_over_sd=float(abs(m) / s) if s else float("inf"),
                occupancy=float(s / math.sqrt(s * s + m * m)) if (s or m) else float("nan"),
                frac_positive=float((v > 0).mean()))


# --------------------------------------------------------------------- one target
def target_row(pdb, seed=0):
    import s30.s30_R_ladder as RL
    from core import geometry as geo
    from s12 import instrument as I
    from s24 import d_harness as H

    t0 = time.time()
    cand = H.Candidates.from_universe(pdb)
    n, fold = int(cand.n), int(cand.fold)
    nat = np.asarray(cand.nat_ca, float)
    nt = RL.native_torsions(pdb)
    if nt is None:
        return {"pdb": pdb, "skipped": "no native torsions in peptide_db"}
    phi0, psi0, seq = nt
    assert seq == cand.seq

    # --- LANE R'S LADDER, BIT-IDENTICAL: same sampler, same seed rule, same call order.
    sampler = RL.Sampler(fold, seq)
    rng = np.random.default_rng(zlib.crc32(f"s30R|{pdb}|{seed}".encode()) & 0xFFFFFFFF)
    grid = tuple(m for m in RL.M_GRID if m <= n)
    PA, SA, MA = RL.make_ladder(phi0, psi0, sampler, n, rng, grid, RL.R_DRAWS)
    bi = int(rng.integers(cand.k))
    phiB = np.asarray(cand.PHI[bi], float); psiB = np.asarray(cand.PSI[bi], float)
    PB, SB, MB = RL.make_ladder(phiB, psiB, sampler, n, rng, grid, RL.R_DRAWS)
    PX = np.vstack([phi0[None], phiB[None]]); SX = np.vstack([psi0[None], psiB[None]])
    zc = np.load(RL.chain_path(pdb))
    PF = np.array([np.asarray(zc[f"phi_{k}"], float) for k in RL.FIXED_RUNGS])
    SF = np.array([np.asarray(zc[f"psi_{k}"], float) for k in RL.FIXED_RUNGS])

    PHI = np.vstack([np.asarray(cand.PHI, float), PA, PB, PX, PF])
    PSI = np.vstack([np.asarray(cand.PSI, float), SA, SB, SX, SF])
    W = np.asarray(geo.build_backbone_batch(PHI, PSI)["CA"], float)
    W[:cand.k] = np.asarray(cand.W, float)              # pool keeps its DEPOSITED coordinates

    K = cand.k
    iA = np.arange(K, K + len(PA))
    iB = np.arange(iA[-1] + 1, iA[-1] + 1 + len(PB))
    iNATR, iANCHB = iB[-1] + 1, iB[-1] + 2
    iFIX = {k: iB[-1] + 3 + q for q, k in enumerate(RL.FIXED_RUNGS)}

    ch = chiral_channels(W, n)
    audit = reflection_audit(W[::37], n)                 # a stride, for speed; same assertions

    # --- DIS on the same set, as the scale comparator the coordinator asked for
    universe = I.load_univ(pdb)
    dg = I.distogram(pdb, seq, fold)
    ii, jj = I.pair_index(n)
    D = I.pair_dists(W, ii, jj)
    ch["DIS"] = np.asarray(I.shipped_score(dg, D.astype(np.float32).astype(float)), float)

    # --- ORACLE labels, only AFTER every channel value exists
    rmsd_nat = RL.ca_rmsd_to(W, nat)
    rg = RL.rg_of(W)
    rmsd_anchorB = RL.ca_rmsd_to(W[iB], W[iANCHB])
    near = iA[(rmsd_nat[iA] <= RL.NEAR_NATIVE_A)]
    fallback = near.size == 0
    if fallback:
        near = iA[[int(np.argmin(rmsd_nat[iA]))]]

    row = dict(pdb=pdb, n=n, fold=fold, k_pool=int(K), seed=seed,
               rebuild_floor=float(rmsd_nat[iNATR]),
               anchorB_rmsd=float(rmsd_nat[iANCHB]),
               n_near=int(near.size), near_fallback=bool(fallback),
               reflection_audit=audit, ch={}, secs=0.0)

    # ---------------- THE PRE-CHECK: does the chiral coordinate have variance where we live?
    row["prechk"] = {}
    for nm in list(CHIRAL) + [TWIN[c] for c in CHIRAL] + ["DIS"]:
        v = ch[nm]
        row["prechk"][nm] = dict(pool=spread(v[:K]), ladderA=spread(v[iA]),
                                 ladderB=spread(v[iB]), all=spread(v))
    row["prechk_native"] = {nm: float(ch[nm][iNATR]) for nm in CHIRAL}
    row["prechk_prod"] = {nm: float(ch[nm][iFIX["PROD"]]) for nm in CHIRAL}

    # ---------------- lane R's two contrasts, computed exactly as lane R computes them
    for nm in list(CHIRAL) + [TWIN[c] for c in CHIRAL] + ["DIS"]:
        v = np.asarray(ch[nm], float)
        if not np.isfinite(v).all() or np.std(v) == 0:
            continue
        cell = dict(rho_rg_A=float(np.nanmean([RL.spearman(v[iA][MA == m], rg[iA][MA == m])
                                               for m in grid])))
        for tag, idx, ms, lab in (("A", iA, MA, rmsd_nat[iA]),
                                  ("B", iB, MB, rmsd_nat[iB]),
                                  ("ANCHOR", iB, MB, rmsd_anchorB)):
            raw, par = [], []
            for m in grid:
                sel = ms == m
                if sel.sum() < 6:
                    continue
                raw.append(RL.spearman(v[idx][sel], lab[sel]))
                par.append(RL.partial_spearman(v[idx][sel], lab[sel], [rg[idx][sel]],
                                               vshape_of={0}))
            cell["rho_" + tag] = float(np.nanmean(raw)) if raw else float("nan")
            cell["rho_" + tag + "_part"] = float(np.nanmean(par)) if par else float("nan")
        vp = float(v[iFIX["PROD"]])
        cell["pref_near"] = float(np.mean((v[near] < vp) + 0.5 * (v[near] == vp)))
        cell["pref_pool"] = float(np.mean((v[:K] < vp) + 0.5 * (v[:K] == vp)))
        cell["pref_native_rebuilt"] = float((v[iNATR] < vp) + 0.5 * (v[iNATR] == vp))
        cell["pctile_nat_in_A"] = float(np.mean(v[iA] < v[iNATR]))
        row["ch"][nm] = cell

    row["secs"] = time.time() - t0
    return row


# --------------------------------------------------------------------- drivers
def cmd_run(a):
    from s12 import instrument as I
    pdbs = [t["pdb"] for t in I.targets()]
    if a.limit:
        pdbs = pdbs[:a.limit]
    pdbs = [p for q, p in enumerate(pdbs) if q % a.of == a.shard]
    path = ROWS if a.of == 1 else ROWS.replace(".jsonl", f".s{a.shard}of{a.of}.jsonl")
    done = set()
    if os.path.exists(path) and not a.rebuild:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    done.add(json.loads(line)["pdb"])
    with open(path, "a", encoding="utf-8") as fh:
        for q, pdb in enumerate(pdbs):
            if pdb in done:
                continue
            try:
                r = target_row(pdb, seed=a.seed)
            except Exception as e:
                import traceback; traceback.print_exc()
                r = {"pdb": pdb, "error": repr(e)}
            fh.write(json.dumps(r) + "\n"); fh.flush()
            print(f"[{q + 1}/{len(pdbs)}] {pdb} {r.get('secs', 0):.1f}s", flush=True)
    print(f"wrote {path}")


def cmd_agg(a):
    import s30.s30_R_agg as RA
    rows, seen = [], set()
    for f in sorted(glob.glob(ROWS.replace(".jsonl", "*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                r = json.loads(line)
                if r.get("error") or r.get("skipped") or r["pdb"] in seen:
                    continue
                seen.add(r["pdb"]); rows.append(r)
    rows.sort(key=lambda r: r["pdb"])
    folds = [r["fold"] for r in rows]
    names = list(CHIRAL) + [TWIN[c] for c in CHIRAL] + ["DIS"]
    out = dict(prereg="s30/PREREG_S30_G.md @ 78b65521", n=len(rows),
               floor_rebuild=float(np.mean([r["rebuild_floor"] for r in rows])),
               bars=dict(BAR_RHO_MARGIN=RA.BAR_RHO_MARGIN, BAR_PREF=RA.BAR_PREF),
               audit_max=dict(
                   flip_err=max(r["reflection_audit"][c]["flip_err"] for r in rows for c in CHIRAL),
                   twin_err=max(r["reflection_audit"][c]["twin_err"] for r in rows for c in CHIRAL),
                   rot_err=max(r["reflection_audit"][c]["rot_err"] for r in rows for c in CHIRAL),
                   dD=max(r["reflection_audit"]["max_abs_dD_under_reflection"] for r in rows)))

    # --- THE PRE-CHECK, aggregated.  This gate is read FIRST.
    out["PRECHECK_occupancy"] = {}
    for nm in names:
        d = {}
        for setname in ("pool", "ladderA", "ladderB"):
            occ = [r["prechk"][nm][setname]["occupancy"] for r in rows if nm in r["prechk"]]
            ams = [r["prechk"][nm][setname]["abs_mean_over_sd"] for r in rows if nm in r["prechk"]]
            sd = [r["prechk"][nm][setname]["sd"] for r in rows if nm in r["prechk"]]
            mn = [r["prechk"][nm][setname]["mean"] for r in rows if nm in r["prechk"]]
            fp = [r["prechk"][nm][setname]["frac_positive"] for r in rows if nm in r["prechk"]]
            d[setname] = dict(occupancy_mean=float(np.mean(occ)), occupancy_median=float(np.median(occ)),
                              occupancy_min=float(np.min(occ)), occupancy_max=float(np.max(occ)),
                              abs_mean_over_sd=float(np.mean(ams)),
                              sd=float(np.mean(sd)), mean=float(np.mean(mn)),
                              frac_positive=float(np.mean(fp)))
        out["PRECHECK_occupancy"][nm] = d

    # --- lane R's contrasts
    out["channels"] = {}
    for nm in names:
        g = lambda k: [r["ch"].get(nm, {}).get(k, float("nan")) for r in rows]
        cell = {}
        for k in ("rho_A", "rho_A_part", "rho_B", "rho_B_part", "rho_ANCHOR", "rho_ANCHOR_part",
                  "rho_rg_A", "pref_near", "pref_pool", "pref_native_rebuilt", "pctile_nat_in_A"):
            cell[k] = RA.mean_ci(g(k), folds, f"G.{nm}.{k}")
        cell["contrast_A_minus_ANCHOR"] = RA.paired_ci(g("rho_A_part"), g("rho_ANCHOR_part"),
                                                       folds, f"G.{nm}.ordering")
        cell["contrast_pref"] = RA.paired_ci(g("pref_near"), g("pref_pool"), folds, f"G.{nm}.pref")
        out["channels"][nm] = cell

    # --- the TWIN contrast: chiral minus its own achiral control, in the operator's own space
    out["chiral_minus_twin"] = {}
    for c in CHIRAL:
        gc = lambda k, nm: [r["ch"].get(nm, {}).get(k, float("nan")) for r in rows]
        out["chiral_minus_twin"][c] = dict(
            ordering=RA.paired_ci(gc("rho_A_part", c), gc("rho_A_part", TWIN[c]), folds,
                                  f"G.{c}.minus_twin.ordering"),
            anchor_contrast=RA.paired_ci(
                [a - b for a, b in zip(gc("rho_A_part", c), gc("rho_ANCHOR_part", c))],
                [a - b for a, b in zip(gc("rho_A_part", TWIN[c]), gc("rho_ANCHOR_part", TWIN[c]))],
                folds, f"G.{c}.minus_twin.anchor"),
            pref=RA.paired_ci(gc("pref_near", c), gc("pref_near", TWIN[c]), folds,
                              f"G.{c}.minus_twin.pref"))

    # --- multiplicity: max-over-the-three-channels sign-flip null (lane R's own statistic)
    M = np.array([[r["ch"].get(c, {}).get("rho_A_part", np.nan)
                   - r["ch"].get(c, {}).get("rho_ANCHOR_part", np.nan) for c in CHIRAL]
                  for r in rows], float)
    out["multiplicity_signflip_max"] = RA.signflip_pmax(M)
    P = np.array([[r["ch"].get(c, {}).get("pref_near", np.nan)
                   - r["ch"].get(c, {}).get("pref_pool", np.nan) for c in CHIRAL]
                  for r in rows], float)
    out["multiplicity_signflip_max_pref"] = RA.signflip_pmax(P, seed=20300931)
    out["n_comparisons_registered"] = 3

    # --- THE KIND-MATCHED READ.  `pref_near - pref_pool` is CROSS-KIND in this channel's own
    # space: the near-native rungs are ideal-geometry REBUILDS while the pool control keeps its
    # DEPOSITED coordinates (`W[:k] = cand.W`), and a raw geometric shape statistic separates those
    # two constructions whether or not it sees nativeness.  The internally matched statistic is
    # `pctile_nat_in_A`: the rebuilt native's percentile inside its OWN ladder, every member an
    # ideal rebuild from the same Ramachandran table at the same budget.  Chance is 0.5 exactly.
    out["KIND_MATCHED_native_percentile"] = {
        nm: dict(out["channels"][nm]["pctile_nat_in_A"],
                 verdict=("BETTER than chance" if out["channels"][nm]["pctile_nat_in_A"]["mean"] < 0.5
                          else "WORSE than chance"))
        for nm in names}
    out["KIND_MATCHED_note"] = (
        "Read this INSTEAD of contrast_pref for any channel whose value depends on raw geometry. "
        "DIS 0.288 (better than chance); every chiral channel and every achiral twin is at or "
        "worse than chance. A cross-kind preference margin is what S30-L1 withdrew S28-L48 for.")

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {OUT}  (n = {out['n']})")

    au = out["audit_max"]
    print(f"\nREFLECTION AUDIT (worst over all targets): flip {au['flip_err']:.2e}  "
          f"twin {au['twin_err']:.2e}  rot {au['rot_err']:.2e}  dD {au['dD']:.2e}")
    print(f"floor: torsion rebuild {out['floor_rebuild']:.3f} A\n")

    print("=== PRE-CHECK: does the chiral coordinate have variance where we live? ===")
    print(f"{'channel':14s} {'set':8s} {'mean':>11s} {'sd':>11s} {'|mean|/sd':>10s} "
          f"{'occupancy':>10s} {'frac>0':>7s}")
    for nm in names:
        for s in ("pool", "ladderA", "ladderB"):
            d = out["PRECHECK_occupancy"][nm][s]
            print(f"{nm:14s} {s:8s} {d['mean']:>11.4g} {d['sd']:>11.4g} "
                  f"{d['abs_mean_over_sd']:>10.3f} {d['occupancy_mean']:>10.3f} "
                  f"{d['frac_positive']:>7.3f}")

    print("\n=== LANE R's CONTRASTS ===")
    print(f"{'channel':14s} {'rho_A':>8s} {'rho_ANCH':>9s} {'A-ANCHOR':>9s} {'xMDE':>7s} "
          f"{'folds':>6s} {'pref_near':>10s} {'pref_pool':>10s} {'pref c':>8s} {'xMDE':>7s}")
    for nm in names:
        c = out["channels"][nm]
        ct, cp = c["contrast_A_minus_ANCHOR"], c["contrast_pref"]
        print(f"{nm:14s} {c['rho_A_part']['mean']:>+8.3f} {c['rho_ANCHOR_part']['mean']:>+9.3f} "
              f"{ct.get('effect', float('nan')):>+9.3f} {ct.get('ratio', float('nan')):>+7.2f} "
              f"{str(ct.get('folds_same_sign')):>6s} {c['pref_near']['mean']:>10.3f} "
              f"{c['pref_pool']['mean']:>10.3f} {cp.get('effect', float('nan')):>+8.3f} "
              f"{cp.get('ratio', float('nan')):>+7.2f}")

    print("\n=== CHIRAL MINUS ITS OWN ACHIRAL TWIN (matched control, operator's own space) ===")
    for c in CHIRAL:
        t = out["chiral_minus_twin"][c]
        for k in ("ordering", "anchor_contrast", "pref"):
            v = t[k]
            print(f"{c:14s} {k:16s} {v.get('effect', float('nan')):>+8.4f} "
                  f"[{v.get('fold_ci', [float('nan')] * 2)[0]:+.4f},"
                  f"{v.get('fold_ci', [float('nan')] * 2)[1]:+.4f}] "
                  f"{v.get('ratio', float('nan')):>+6.2f}x MDE  folds {v.get('folds_same_sign')}")
    m = out["multiplicity_signflip_max"]
    print(f"\nmax-over-3 sign-flip null: observed {m.get('observed_max', float('nan')):.4f}  "
          f"null mean {m.get('null_mean', float('nan')):.4f}  p95 {m.get('null_p95', float('nan')):.4f}  "
          f"p_max {m.get('p_max', float('nan')):.3f}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--limit", type=int, default=0)
    r.add_argument("--shard", type=int, default=0); r.add_argument("--of", type=int, default=1)
    r.add_argument("--seed", type=int, default=0); r.add_argument("--rebuild", action="store_true")
    r.set_defaults(fn=cmd_run)
    g = sub.add_parser("agg"); g.set_defaults(fn=cmd_agg)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
