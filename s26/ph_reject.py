"""s26/ph_reject.py -- AMBER AS A STERIC REJECT FILTER, NOT A RANKER (lane PH, S26).

Pre-registered in `s26/PREREG_amber_reject.md`, written before any RMSD in this file was read.

THE CLAIM UNDER TEST.  Both physics energies rank worse than noise (s25 PHYS section 1.4-1.5),
but AMBER's difficulty is a steric singularity: 58.6% of every pool is above 1e4 kcal/mol and
ten candidates carry 99.7% of its variance (s25/results/phys_landscape.json).  A BINARY reject at
a FIXED PHYSICAL THRESHOLD uses the one thing the energy is known to measure (its worst clash)
and none of the ordering it is known not to carry.  It differs from every earlier filter:

  * s19 Q3 (`s19/agentC_reject.py`) rejected a FIXED COUNT r in {2,5,10,19} of the top-75 by
    rank of Legacy steric / min-heavy / the AMBER single point, no refill, point cloud.  All
    worse than matched-random; AMBER at r=10 +0.0109 [-0.0049,+0.0285] vs random.
  * s24 D1-C rejected a FIXED FRACTION of the K=500 pool by which energy PREFERS a candidate,
    n=30, five partitions, all harmful in sign (+0.26 to +0.55 A), none past its own MDE.
  Here the count is set by physics, not by rank: it is zero on a target whose top-75 has no
  candidate above the threshold and large where the retrieval is clash-ridden; the retained set
  is REFILLED from the next-ranked survivors so m stays 75 (arm R), and reported without refill
  (arm S); the endpoint is reported on the built chain as well as the point cloud.

ARMS, all native-free, per threshold T in (1e3, 1e4, 1e5, 1e6) kcal/mol, PRIMARY T = 1e4:
    anchor    the shipped top-75 (production `sub`), 3.0483 point cloud / 3.2148 built chain
    R_T       reject e > T, refill from the next-ranked survivors of the K=500 pool until 75
              (identically: the top-75 of the pool's survivors in the shipped score order)
    S_T       reject e > T from the shipped top-75, no refill (the set shrinks)
    RANDS_T   CONTROL matched to S_T: reject the SAME COUNT at random from the top-75, 16 draws
    RANDR_T   CONTROL matched to R_T: reject the same count at random from the top-75, refill
              with the next-ranked candidates in score order with no energy judgment, 16 draws
    PERMS_T / PERMR_T   CONTROL: the energy vector permuted within the pool (marginal kept,
              correspondence to candidates destroyed), same threshold, same two operators,
              16 permutations
Zero-information anchor: the shipped top-75 itself.

ENDPOINTS (gated): point-cloud RMSD (`I.coordinate_average` + `I.ca_rmsd`, basis `rmsd_avg`) for
every arm and every draw; built-chain RMSD (`I.project`, basis `rmsd_arm`) for the anchor, R_T
and S_T at every T, and for the four controls at the PRIMARY T with 4 draws (the projection costs
~4.4 s and 16 draws x 4 controls x 126 targets is not affordable; stated here, not discovered).
A retained set identical to the anchor's reuses the anchor's projection (the projection is a
deterministic function of the cloud).

STATISTICS: `ST.compare` paired per target, fold-clustered CI beside iid, MDE, W/L, the
concentration null; the threshold sweep is an ORDER STATISTIC and is quoted through
`ST.best_of_k_within` (split-half transfer), never as its raw minimum.

    python s26/ph_reject.py census                 # native-free, before the gate
    python s26/ph_reject.py cloud [--n N]          # gated: point cloud, every arm
    python s26/ph_reject.py chain [--n N]          # gated: built chain, resumable per target
    python s26/ph_reject.py report                 # the statistics from the cells
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

from s26 import ph_lib as L                                   # noqa: E402

THRESHOLDS = (1e3, 1e4, 1e5, 1e6)
PRIMARY = 1e4
N_DRAWS = 16
N_DRAWS_CHAIN = 4
CENSUS_JSON = os.path.join(L.RESULTS, "ph_reject_census.json")
CLOUD_JSON = os.path.join(L.RESULTS, "ph_reject_cloud.json")
CHAIN_JSON = os.path.join(L.RESULTS, "ph_reject_chain.json")
REPORT_JSON = os.path.join(L.RESULTS, "ph_reject_report.json")


def tkey(T: float) -> str:
    return "1e%d" % int(round(math.log10(T)))


# ------------------------------------------------------------------ the pool, asserted
def load_pool(pdb: str, oracle: bool = False) -> dict:
    """The shipped K=500 pool with the AMBER cache identity ASSERTED bit-for-bit.

    `sub` (the production top-75) is asserted to be the first 75 of the stable argsort of the
    cached distogram score, as a SET; a boundary tie that broke differently would be reported,
    not silently absorbed.
    """
    z = L.cache_amber(pdb)
    if oracle:
        u = L.univ_oracle(pdb)
    else:
        u = L.univ_nativefree(pdb, keys=("W", "order"))
    p = L.pool_idx_from_order(u, k=z["k"])
    assert np.array_equal(p, z["universe_idx"]), f"{pdb}: AMBER cache pool != I.pool_idx"
    assert z["amber_verify_max_rel"] == 0.0, f"{pdb}: cache not verified bit-exact"
    rec = L.prod_record_oracle(pdb) if oracle else L.prod_record_nativefree(pdb)
    sub = np.asarray(rec["sub"], int)
    order = L.argsort_stable(z["score_dist"])
    top = order[:L.M]
    assert set(top.tolist()) == set(sub.tolist()), f"{pdb}: production sub != score top-75 as a set"
    return {"pdb": pdb, "seq": u["seq"], "n": u["n"], "fold": u["fold"],
            "W": u["W"][p], "e": z["e_amber"], "score": z["score_dist"], "order": order,
            "sub": sub, "rec": rec, "nat_ca": u["nat_ca"] if oracle else None}


# ------------------------------------------------------------------ the operators
def reject_refill(order: np.ndarray, e: np.ndarray, T: float, m: int = L.M) -> np.ndarray:
    """Arm R: the first m survivors (e <= T) in score order. Equals 'reject from the top-m and
    refill from the next-ranked survivors'. May return fewer than m if the pool has fewer
    survivors; the caller records that."""
    keep = order[e[order] <= T]
    return np.asarray(keep[:m], int)


def reject_shrink(sub: np.ndarray, e: np.ndarray, T: float) -> np.ndarray:
    """Arm S: the shipped top-m with e > T removed. Empty set falls back to `sub` (recorded)."""
    keep = sub[e[sub] <= T]
    return np.asarray(keep, int)


def random_shrink(sub: np.ndarray, r: int, rng) -> np.ndarray:
    """Control for S: reject r of the top-m at random."""
    sub = np.asarray(sub, int)
    if r <= 0:
        return sub.copy()
    drop = rng.choice(len(sub), size=min(r, len(sub)), replace=False)
    keep = np.ones(len(sub), bool); keep[drop] = False
    return sub[keep]


def random_refill(order: np.ndarray, sub: np.ndarray, r: int, rng, m: int = L.M) -> np.ndarray:
    """Control for R: reject r of the top-m at random, refill with ranks m..m+r-1 in score order
    with no energy judgment."""
    kept = random_shrink(sub, r, rng)
    need = m - len(kept)
    if need <= 0:
        return kept
    tail = [i for i in order[m:] if i not in set(kept.tolist())][:need]
    return np.concatenate([kept, np.asarray(tail, int)])


def permuted_energy(e: np.ndarray, rng) -> np.ndarray:
    return e[rng.permutation(len(e))]


def retained_sets(pool: dict, T: float, n_draws: int = N_DRAWS) -> dict:
    """Every arm's retained index set(s) for one threshold. Native-free."""
    e, order, sub, pdb = pool["e"], pool["order"], pool["sub"], pool["pdb"]
    S = reject_shrink(sub, e, T)
    r = int(len(sub) - len(S))
    R = reject_refill(order, e, T)
    out = {"n_reject": r, "S": S, "S_empty": bool(len(S) == 0), "R": R,
           "R_short": bool(len(R) < L.M),
           "n_survivors_pool": int((e <= T).sum()),
           "refill_depth": int(np.flatnonzero(np.isin(order, R)).max() + 1) if len(R) else 0,
           "R_overlap_anchor": float(len(set(R.tolist()) & set(sub.tolist())) / len(sub))}
    if len(S) == 0:
        S = sub.copy()
    out["S_eff"] = S
    #: PREREG addendum 1 (2026-09-13 00:50, before any RMSD was read): a pool with NO survivor
    #: at T (8 targets at 1e4, 20 at 1e3 in the census) makes R empty.  A deployable operator
    #: must still emit something, and the only native-free choice is "do nothing" (the anchor).
    #: `R_empty` is recorded and the report also shows the finite-only subset.
    out["R_empty"] = bool(len(R) == 0)
    out["R_eff"] = R if len(R) else sub.copy()
    rs = L.stable_rng(pdb, "rands", tkey(T))
    rr = L.stable_rng(pdb, "randr", tkey(T))
    rp = L.stable_rng(pdb, "perm", tkey(T))
    out["RANDS"] = [random_shrink(sub, r, rs) for _ in range(n_draws)]
    out["RANDR"] = [random_refill(order, sub, r, rr) for _ in range(n_draws)]
    perms, ps, pr = [], [], []
    for _ in range(n_draws):
        ep = permuted_energy(e, rp)
        s_ = reject_shrink(sub, ep, T)
        perms.append(int(len(sub) - len(s_)))
        ps.append(s_ if len(s_) else sub.copy())
        pr.append(reject_refill(order, ep, T))
    out["PERMS"] = ps; out["PERMR"] = pr; out["perm_n_reject"] = perms
    return out


# ------------------------------------------------------------------ part 0: the census
CENSUS_KEYS = ("pdb", "n", "fold", "pool_frac_over", "top75_n_reject", "top75_frac_reject",
               "n_survivors_pool", "refill_depth", "R_overlap_anchor", "R_short",
               "rg_anchor", "rg_R", "rg_S", "minsep3_anchor", "minsep3_R", "minsep3_S",
               "perm_n_reject_mean", "e_top75_median", "e_top75_max", "e_pool_median",
               "min_heavy_allatom_top75", "min_heavy_bbcb_top75", "clash_class_top75",
               "rho_e_minheavy_top75", "n_members_minheavy_lt2_allatom",
               "n_members_minheavy_lt2_bbcb")


def _allatom_min_heavy(seq: str, phi: np.ndarray, psi: np.ndarray, min_res_sep: int = 2):
    """Minimum heavy-atom distance of the ideal-geometry rebuild, all heavy atoms placed by
    the reference sidechain builder (`sidechains.build_full_structure`, no OpenMM), over atom
    pairs at residue separation >= min_res_sep; split by class (bb-bb, bb-sc, sc-sc)."""
    from core import geometry as geo
    import sidechains as SC
    bb = geo.build_backbone(np.asarray(phi, float), np.asarray(psi, float))
    full = SC.build_full_structure(seq, bb, add_oxt=False)
    xyz, res, is_bb = [], [], []
    for i, rd in enumerate(full["residues"]):
        for nm, v in rd.items():
            xyz.append(np.asarray(v, float)); res.append(i); is_bb.append(nm in ("N", "CA", "C", "O"))
    X = np.array(xyz); R = np.array(res); B = np.array(is_bb)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    far = np.abs(R[:, None] - R[None, :]) >= min_res_sep
    D = np.where(far, D, np.inf)
    k = int(np.argmin(D)); a, b = divmod(k, len(X))
    cls = "bb-bb" if (B[a] and B[b]) else ("sc-sc" if (not B[a] and not B[b]) else "bb-sc")
    # backbone+CB only, S19's panel definition (N, CA, C, O, CB)
    keep = B | np.array([nm == "CB" for rd in full["residues"] for nm in rd.keys()])
    Dk = D[np.ix_(keep, keep)]
    return float(D.min()), cls, float(Dk.min())


def census_target(t) -> dict:
    pool = load_pool(t["pdb"], oracle=False)
    e, sub, W = pool["e"], pool["sub"], pool["W"]
    u = L.univ_nativefree(t["pdb"], keys=("PHI", "PSI", "order"))
    p = L.pool_idx_from_order(u)
    PHI, PSI = u["PHI"][p], u["PSI"][p]
    row = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
           "e_top75_median": float(np.median(e[sub])), "e_top75_max": float(e[sub].max()),
           "e_pool_median": float(np.median(e)),
           "rg_anchor": float(L.rg_of(W[sub]).mean()),
           "minsep3_anchor": float(L.min_sep3(W[sub]).mean()),
           "pool_frac_over": {}, "top75_n_reject": {}, "top75_frac_reject": {},
           "n_survivors_pool": {}, "refill_depth": {}, "R_overlap_anchor": {}, "R_short": {},
           "rg_R": {}, "rg_S": {}, "minsep3_R": {}, "minsep3_S": {}, "perm_n_reject_mean": {}}
    for T in THRESHOLDS:
        k = tkey(T)
        rs = retained_sets(pool, T)
        row["pool_frac_over"][k] = float((e > T).mean())
        row["top75_n_reject"][k] = rs["n_reject"]
        row["top75_frac_reject"][k] = rs["n_reject"] / len(sub)
        row["n_survivors_pool"][k] = rs["n_survivors_pool"]
        row["refill_depth"][k] = rs["refill_depth"]
        row["R_overlap_anchor"][k] = rs["R_overlap_anchor"]
        row["R_short"][k] = rs["R_short"]
        row["rg_R"][k] = float(L.rg_of(W[rs["R"]]).mean()) if len(rs["R"]) else float("nan")
        row["rg_S"][k] = float(L.rg_of(W[rs["S_eff"]]).mean())
        row["minsep3_R"][k] = float(L.min_sep3(W[rs["R"]]).mean()) if len(rs["R"]) else float("nan")
        row["minsep3_S"][k] = float(L.min_sep3(W[rs["S_eff"]]).mean())
        row["perm_n_reject_mean"][k] = float(np.mean(rs["perm_n_reject"]))
    # where the singularity lives: min heavy-atom distance of every top-75 rebuild
    mh, mk, cls = [], [], []
    for j in sub:
        a, c, b = _allatom_min_heavy(pool["seq"], PHI[j], PSI[j])
        mh.append(a); cls.append(c); mk.append(b)
    mh = np.array(mh); mk = np.array(mk)
    from scipy.stats import spearmanr
    rho = spearmanr(e[sub], mh)[0] if len(sub) > 2 else float("nan")
    row.update({"min_heavy_allatom_top75": {"mean": float(mh.mean()), "min": float(mh.min()),
                                            "median": float(np.median(mh))},
                "min_heavy_bbcb_top75": {"mean": float(mk.mean()), "min": float(mk.min()),
                                         "median": float(np.median(mk))},
                "clash_class_top75": {c: int(sum(1 for x in cls if x == c))
                                      for c in ("bb-bb", "bb-sc", "sc-sc")},
                "clash_class_over_primary": {c: int(sum(1 for x, ee in zip(cls, e[sub])
                                                        if x == c and ee > PRIMARY))
                                             for c in ("bb-bb", "bb-sc", "sc-sc")},
                "rho_e_minheavy_top75": float(rho),
                "n_members_minheavy_lt2_allatom": int((mh < 2.0).sum()),
                "n_members_minheavy_lt2_bbcb": int((mk < 2.0).sum()),
                "n_members_over_primary": int((e[sub] > PRIMARY).sum()),
                "minheavy_allatom_over_primary": float(mh[e[sub] > PRIMARY].mean())
                if (e[sub] > PRIMARY).any() else float("nan"),
                "minheavy_allatom_under_primary": float(mh[e[sub] <= PRIMARY].mean())
                if (e[sub] <= PRIMARY).any() else float("nan")})
    return row


def census(n: int = 0) -> dict:
    tg = L.targets()
    if n:
        tg = tg[:n]
    rows, clk = [], L.Clock()
    for k, t in enumerate(tg):
        r = census_target(t)
        rows.append(r)
        print(f"[{k + 1}/{len(tg)}] {r['pdb']} n={r['n']} reject@1e4 {r['top75_n_reject']['1e4']}/75 "
              f"depth {r['refill_depth']['1e4']} survivors {r['n_survivors_pool']['1e4']} "
              f"minheavy(all) {r['min_heavy_allatom_top75']['min']:.2f} "
              f"rho(e,minheavy) {r['rho_e_minheavy_top75']:+.2f} ({clk():.0f}s)", flush=True)
    summ = summarise_census(rows)
    L.save(CENSUS_JSON, {"what": "native-free census of the steric reject: what each threshold "
                                 "rejects, before any RMSD is read; plus where the AMBER "
                                 "singularity lives in the top-75 rebuilds",
                         "thresholds": list(THRESHOLDS), "primary": PRIMARY,
                         "rows": rows, "summary": summ},
           rows=rows, complete_keys=CENSUS_KEYS, n_expected=126 if not n else n,
           module_file=__file__)
    print(json.dumps(summ, indent=1))
    return summ


def summarise_census(rows) -> dict:
    out = {"n": len(rows), "per_threshold": {}}
    for T in THRESHOLDS:
        k = tkey(T)
        nr = np.array([r["top75_n_reject"][k] for r in rows], float)
        out["per_threshold"][k] = {
            "pool_frac_over": L.mean_se([r["pool_frac_over"][k] for r in rows]),
            "top75_n_reject": L.mean_se(nr),
            "top75_frac_reject": L.mean_se([r["top75_frac_reject"][k] for r in rows]),
            "targets_with_zero_reject": int((nr == 0).sum()),
            "targets_with_all_rejected": int(sum(1 for r in rows if r["top75_n_reject"][k] >= 75)),
            "n_survivors_pool": L.mean_se([r["n_survivors_pool"][k] for r in rows]),
            "targets_R_short": int(sum(1 for r in rows if r["R_short"][k])),
            "refill_depth": L.mean_se([r["refill_depth"][k] for r in rows]),
            "R_overlap_anchor": L.mean_se([r["R_overlap_anchor"][k] for r in rows]),
            "d_rg_R_minus_anchor": L.mean_se([r["rg_R"][k] - r["rg_anchor"] for r in rows]),
            "d_rg_S_minus_anchor": L.mean_se([r["rg_S"][k] - r["rg_anchor"] for r in rows]),
            "d_minsep3_R_minus_anchor": L.mean_se([r["minsep3_R"][k] - r["minsep3_anchor"]
                                                   for r in rows]),
            "d_minsep3_S_minus_anchor": L.mean_se([r["minsep3_S"][k] - r["minsep3_anchor"]
                                                   for r in rows]),
            "perm_n_reject_mean": L.mean_se([r["perm_n_reject_mean"][k] for r in rows]),
        }
    cls = {c: int(sum(r["clash_class_top75"][c] for r in rows)) for c in ("bb-bb", "bb-sc", "sc-sc")}
    clsP = {c: int(sum(r["clash_class_over_primary"][c] for r in rows))
            for c in ("bb-bb", "bb-sc", "sc-sc")}
    out["singularity"] = {
        "closest_contact_class_all_members": cls,
        "closest_contact_class_members_over_1e4": clsP,
        "min_heavy_allatom_mean": L.mean_se([r["min_heavy_allatom_top75"]["mean"] for r in rows]),
        "min_heavy_allatom_worst_per_target": L.mean_se([r["min_heavy_allatom_top75"]["min"]
                                                         for r in rows]),
        "min_heavy_bbcb_mean": L.mean_se([r["min_heavy_bbcb_top75"]["mean"] for r in rows]),
        "min_heavy_bbcb_worst_per_target": L.mean_se([r["min_heavy_bbcb_top75"]["min"]
                                                      for r in rows]),
        "members_minheavy_lt2_allatom_per_target": L.mean_se(
            [r["n_members_minheavy_lt2_allatom"] for r in rows]),
        "members_minheavy_lt2_bbcb_per_target": L.mean_se(
            [r["n_members_minheavy_lt2_bbcb"] for r in rows]),
        "rho_e_minheavy_within_top75": L.mean_se([r["rho_e_minheavy_top75"] for r in rows]),
        "minheavy_allatom_members_over_1e4": L.mean_se(
            [r["minheavy_allatom_over_primary"] for r in rows]),
        "minheavy_allatom_members_under_1e4": L.mean_se(
            [r["minheavy_allatom_under_primary"] for r in rows]),
    }
    return out


# ------------------------------------------------------------------ gated endpoints
def cloud_rmsd(W: np.ndarray, idx: np.ndarray, nat: np.ndarray) -> float:
    """POINT CLOUD basis (`rmsd_avg`): the uniform coordinate average of the retained windows in
    their own medoid frame, scored ORACLE against the native."""
    from s12 import instrument as I
    C, _ = I.coordinate_average(W[np.asarray(idx, int)])
    return float(I.ca_rmsd(C, nat))


def chain_of(W: np.ndarray, idx: np.ndarray, seq: str, fold: int) -> np.ndarray:
    """BUILT CHAIN basis (`rmsd_arm`): the production projection of the cloud (`I.project`,
    ramah at 0.3, multi-start, exact gradient)."""
    from s12 import instrument as I
    C, _ = I.coordinate_average(W[np.asarray(idx, int)])
    return np.asarray(I.project(C, seq, int(fold))["ca"], float)


CLOUD_KEYS = ("pdb", "n", "fold", "anchor", "arms")


def cloud_target(t) -> dict:
    L.require_gate("ph_reject cloud")
    pool = load_pool(t["pdb"], oracle=True)
    W, nat, sub = pool["W"], pool["nat_ca"], pool["sub"]
    anchor = cloud_rmsd(W, sub, nat)
    rec = pool["rec"]
    assert abs(anchor - float(rec["rmsd_avg"])) < 1e-6, \
        f"{t['pdb']}: anchor point cloud {anchor} != production rmsd_avg {rec['rmsd_avg']}"
    arms = {}
    for T in THRESHOLDS:
        k = tkey(T)
        rs = retained_sets(pool, T)
        a = {"n_reject": rs["n_reject"], "S_empty": rs["S_empty"], "R_short": rs["R_short"],
             "R_empty": rs["R_empty"],
             "R": cloud_rmsd(W, rs["R_eff"], nat),
             "S": cloud_rmsd(W, rs["S_eff"], nat),
             "RANDS": float(np.mean([cloud_rmsd(W, s, nat) for s in rs["RANDS"]])),
             "RANDR": float(np.mean([cloud_rmsd(W, s, nat) for s in rs["RANDR"]])),
             "PERMS": float(np.mean([cloud_rmsd(W, s, nat) for s in rs["PERMS"]])),
             "PERMR": float(np.mean([cloud_rmsd(W, s, nat) for s in rs["PERMR"]])),
             "perm_n_reject_mean": float(np.mean(rs["perm_n_reject"]))}
        arms[k] = a
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "basis": "point_cloud",
            "anchor": anchor, "arms": arms}


def cloud(n: int = 0) -> None:
    L.require_gate("ph_reject cloud")
    tg = L.targets()
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_reject_cloud")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = cloud_target(t)
        L.write_cell("ph_reject_cloud", t["pdb"], r)
        a = r["arms"][tkey(PRIMARY)]
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} anchor {r['anchor']:.3f} R {a['R']:.3f} S {a['S']:.3f} "
              f"randS {a['RANDS']:.3f} randR {a['RANDR']:.3f} permS {a['PERMS']:.3f} "
              f"({clk():.0f}s)", flush=True)
    rows = list(L.read_cells("ph_reject_cloud").values())
    L.save(CLOUD_JSON, {"what": "steric reject, POINT CLOUD basis, every arm, 16 draws",
                        "thresholds": list(THRESHOLDS), "primary": PRIMARY, "n_draws": N_DRAWS,
                        "rows": rows},
           rows=rows, complete_keys=CLOUD_KEYS, n_expected=126, module_file=__file__)


CHAIN_KEYS = ("pdb", "n", "fold", "anchor", "anchor_vs_prod", "arms")


def chain_target(t) -> dict:
    L.require_gate("ph_reject chain")
    from s12 import instrument as I
    pool = load_pool(t["pdb"], oracle=True)
    W, nat, sub, seq, fold = pool["W"], pool["nat_ca"], pool["sub"], pool["seq"], pool["fold"]
    rec = pool["rec"]
    t0 = time.time()
    ca0 = chain_of(W, sub, seq, fold)
    anchor = float(I.ca_rmsd(ca0, nat))
    dev_prod = float(I.ca_rmsd(ca0, np.asarray(rec["ca"], float)))
    cache = {tuple(sorted(sub.tolist())): anchor}

    def rmsd_of(idx):
        key = tuple(sorted(np.asarray(idx, int).tolist()))
        if key not in cache:
            cache[key] = float(I.ca_rmsd(chain_of(W, idx, seq, fold), nat))
        return cache[key]

    arms = {}
    for T in THRESHOLDS:
        k = tkey(T)
        rs = retained_sets(pool, T, n_draws=N_DRAWS_CHAIN)
        a = {"n_reject": rs["n_reject"], "S_empty": rs["S_empty"], "R_short": rs["R_short"],
             "R_empty": rs["R_empty"],
             "R": rmsd_of(rs["R_eff"]),
             "S": rmsd_of(rs["S_eff"])}
        if T == PRIMARY:
            a["RANDS"] = float(np.mean([rmsd_of(s) for s in rs["RANDS"]]))
            a["RANDR"] = float(np.mean([rmsd_of(s) for s in rs["RANDR"]]))
            a["PERMS"] = float(np.mean([rmsd_of(s) for s in rs["PERMS"]]))
            a["PERMR"] = float(np.mean([rmsd_of(s) for s in rs["PERMR"]]))
            a["n_draws"] = N_DRAWS_CHAIN
        arms[k] = a
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "basis": "built_chain",
            "anchor": anchor, "anchor_vs_prod": dev_prod, "prod_rmsd_arm": float(rec["rmsd_arm"]),
            "n_projections": len(cache), "wall": round(time.time() - t0, 1), "arms": arms}


def chain(n: int = 0) -> None:
    L.require_gate("ph_reject chain")
    tg = L.targets()
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_reject_chain")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = chain_target(t)
        L.write_cell("ph_reject_chain", t["pdb"], r)
        a = r["arms"][tkey(PRIMARY)]
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} anchor {r['anchor']:.4f} (prod {r['prod_rmsd_arm']:.4f}, "
              f"dev {r['anchor_vs_prod']:.1e}) R {a['R']:.3f} S {a['S']:.3f} randS {a['RANDS']:.3f} "
              f"randR {a['RANDR']:.3f} proj {r['n_projections']} wall {r['wall']}s ({clk():.0f}s)",
              flush=True)
    rows = list(L.read_cells("ph_reject_chain").values())
    L.save(CHAIN_JSON, {"what": "steric reject, BUILT CHAIN basis (I.project); controls at the "
                                "primary threshold with 4 draws",
                        "thresholds": list(THRESHOLDS), "primary": PRIMARY,
                        "n_draws_chain": N_DRAWS_CHAIN, "rows": rows},
           rows=rows, complete_keys=CHAIN_KEYS, n_expected=126, module_file=__file__)


# ------------------------------------------------------------------ the statistics
def report() -> dict:
    from s24 import stats_lib as ST
    out = {}
    for basis, path, name in (("point_cloud", CLOUD_JSON, "cloud"), ("built_chain", CHAIN_JSON, "chain")):
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        rows = sorted(d["rows"], key=lambda r: r["pdb"])
        pdbs = [r["pdb"] for r in rows]
        folds = L.folds_of(pdbs)
        anchor = np.array([r["anchor"] for r in rows])
        res = {"n": len(rows), "anchor_mean": float(anchor.mean()), "complete": d.get("complete")}
        print(f"\n===== {basis}  n={len(rows)}  anchor {anchor.mean():.4f}")
        for T in THRESHOLDS:
            k = tkey(T)
            arms = [r["arms"][k] for r in rows]
            resT = {"n_reject": L.mean_se([a["n_reject"] for a in arms]),
                    "n_R_empty_fallback": int(sum(1 for a in arms if a.get("R_empty"))),
                    "n_S_empty_fallback": int(sum(1 for a in arms if a.get("S_empty"))),
                    "n_targets_R_moved": int(sum(1 for a in arms
                                                 if a["n_reject"] > 0 and not a.get("R_empty"))),
                    "n_targets_S_moved": int(sum(1 for a in arms
                                                 if a["n_reject"] > 0 and not a.get("S_empty")))}
            for arm in ("R", "S", "RANDS", "RANDR", "PERMS", "PERMR"):
                if not all(arm in a for a in arms):
                    continue
                v = np.array([a[arm] for a in arms], float)
                ok = np.isfinite(v)
                c = ST.compare(v[ok], anchor[ok], folds[ok], names=[p for p, o in zip(pdbs, ok) if o],
                               label=f"{basis} {arm}@{k} minus anchor (negative = arm better)")
                print(ST.fmt(c))
                resT[f"{arm}_vs_anchor"] = {q: c[q] for q in c if q != "concentration"}
                resT[f"{arm}_vs_anchor"]["concentration"] = c["concentration"]
            for arm, ctrl in (("S", "RANDS"), ("R", "RANDR"), ("S", "PERMS"), ("R", "PERMR")):
                if not all(arm in a and ctrl in a for a in arms):
                    continue
                v = np.array([a[arm] for a in arms], float); w = np.array([a[ctrl] for a in arms], float)
                ok = np.isfinite(v) & np.isfinite(w)
                c = ST.compare(v[ok], w[ok], folds[ok], names=[p for p, o in zip(pdbs, ok) if o],
                               label=f"{basis} {arm}@{k} minus {ctrl}@{k} (matched control)")
                print(ST.fmt(c))
                resT[f"{arm}_vs_{ctrl}"] = {q: c[q] for q in c}
            res[k] = resT
        # the threshold sweep is an order statistic
        for arm in ("R", "S"):
            M = np.array([[r["arms"][tkey(T)][arm] for T in THRESHOLDS] for r in rows], float)
            ok = np.isfinite(M).all(1)
            if ok.sum() > 10:
                b = ST.best_of_k_within(M[ok])
                res[f"sweep_{arm}"] = b
                print(f"  sweep {arm}: oracle min over {len(THRESHOLDS)} thresholds "
                      f"{b['observed_gain']:+.4f}, split-half transfer {b['split_half']:+.4f} "
                      f"({100 * b['split_half_frac']:.0f}%), k_eff {b['k_eff']:.2f} -> {b['verdict']}")
        out[basis] = res
    L.save(REPORT_JSON, {"what": "steric reject statistics", "report": out}, module_file=__file__)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("census", "cloud", "chain", "report"))
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "census":
        census(a.n)
    elif a.mode == "cloud":
        cloud(a.n)
    elif a.mode == "chain":
        chain(a.n)
    else:
        report()
