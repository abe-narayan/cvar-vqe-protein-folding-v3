"""s26/p_c4.py -- PREREG_C4: routers for the per-target set size m* and scale s* on a feature set
no previous router used.

NEW FEATURE BLOCKS (native-free; `features` caches them to s26/results/p_c4_features.json):
  pax_*   principal-axis spread of the shipped top-75 cloud: eigenvalue fractions and anisotropy of
          the average cloud's coordinate covariance, and the members' per-axis dispersion about the
          average in that frame (the S22 `pool_spread` is a scalar RMSD; this is its shape)
  sim_*   retrieval-score entropy (softmax over the K=500 BLOSUM sims) and normalised gaps
  dg_ent_* per-pair 17-bin entropy of the shipped posterior (mean/max/sd, multimodal fraction)
  con_*   ESM-2 650M contact-map statistics
(sim_*, dg_ent_*, con_* are built by s26/p_b3.features_one; pax_* here.)  The S22 feature set
(`s22/results/routerdata.json`) is carried as the HARNESS POSITIVE CONTROL: through this code it
must reproduce ~0, as S22 L7 measured.

LABELS (read only inside `run`, gated on "PHASE 0 SIGNED OFF"):
  m*   per target from s12/results/agg_surface.json `shipped|m<m>`, m in MS (15 rungs; point cloud)
  s*   per target from s23/results/errdecomp.json `s_star` (closed form)

ROUTERS (linear, nested leave-fold-out, s26/p_stats):
  A  multi-output ridge of the 15 per-m RMSDs -> argmin of the prediction (S12's construction)
  B  ridge on log m* -> nearest grid rung
  S  ridge on s*
Endpoints: routed m -> the persisted surface value (point cloud, no re-emission) and the built
chain through I.project of the routed-m average (gated); routed s -> RMSD(c*s, t) in errdecomp's
frame (gated).  Nulls: fixed m=75 / s=1; label permutation (rows permuted across targets).

    python s26/p_c4.py features   (native-free, allowed now)
    python s26/p_c4.py selftest   (synthetic labels; no RMSD read)
    python s26/p_c4.py run        (after sign-off)
"""
from __future__ import annotations

import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from s26 import p_stats as PS                        # noqa: E402
from s26 import p_ladder as L                        # noqa: E402
from s26 import p_b3 as B3                           # noqa: E402

RES = os.path.join(HERE, "results")
FEAT = os.path.join(RES, "p_c4_features.json")
OUT = os.path.join(RES, "p_c4.json")
MS = [1, 2, 3, 5, 8, 12, 20, 35, 50, 75, 110, 150, 220, 300, 500]
NEW_BLOCKS = {"pax": "pax_", "sim": "sim_", "dgent": "dg_", "con": "con_"}
OLD_FEATS = ["n", "rg_disto", "rg_pool_mean", "rg_pool_sd", "rg_gap", "rg_z", "score_mean", "score_sd",
             "score_gap01", "score_gap_boundary", "score_iqr_over_range", "score_skew", "pool_spread",
             "sim_mean_top75", "sim_mean_pool"]
N_PERM = 200


def pax_features(u, rec):
    """Principal-axis spread of the shipped top-75 cloud.  Native-free."""
    idx = I.pool_idx(u); sub = np.asarray(rec["sub"], int)
    top = np.asarray(u["W"][idx][sub], float)
    P = I.pairwise_rmsd(top); b = I.medoid(P)
    Wm = I.superpose_batch(top, top[b]); C = Wm.mean(0)
    Cc = C - C.mean(0)
    ev, vec = np.linalg.eigh(Cc.T @ Cc / len(Cc))
    ev = ev[::-1]; vec = vec[:, ::-1]
    fr = ev / max(ev.sum(), 1e-12)
    dev = (Wm - C[None]) @ vec                                # members' deviation in the cloud's frame
    ax_sd = np.sqrt((dev ** 2).mean((0, 1)))                  # per-axis dispersion of members
    return {"pax_f1": float(fr[0]), "pax_f2": float(fr[1]), "pax_f3": float(fr[2]),
            "pax_aniso": float(np.sqrt(ev[0] / max(ev[2], 1e-12))),
            "pax_sd1": float(ax_sd[0]), "pax_sd2": float(ax_sd[1]), "pax_sd3": float(ax_sd[2]),
            "pax_sd_ratio": float(ax_sd[2] / max(ax_sd[0], 1e-12)),
            "pax_spread_over_rg": float(np.sqrt((dev ** 2).sum(-1).mean()) / max(np.sqrt(ev.sum()), 1e-12))}


def build_features(verbose=True):
    if os.path.exists(FEAT):
        return json.load(open(FEAT))
    base = B3.build_features(verbose=verbose)
    old = {r["pdb"]: r for r in json.load(open(os.path.join(ROOT, "s22", "results", "routerdata.json")))["rows"]}
    rows = []
    for c, r in enumerate(base["rows"]):
        u = I.load_univ(r["pdb"]); rec = I.shipped_record(r["pdb"])
        row = {k: v for k, v in r.items() if k in ("pdb", "fold") or k.startswith(("sim_", "dg_", "con_", "n"))}
        row.update(pax_features(u, rec))
        row.update({"old_" + k: float(old[r["pdb"]][k]) for k in OLD_FEATS})
        rows.append(row); del u
        if verbose and (c + 1) % 25 == 0:
            print("  pax %d/126" % (c + 1), flush=True)
    names = [k for k in rows[0] if k not in ("pdb", "fold")]
    out = {"rows": rows, "names": names}
    ST.save_atomic(FEAT, out, complete_keys=names, rows=rows, n_expected=126, module_file=__file__)
    return out


def blocks(names):
    out = {"new_all": [k for k in names if not k.startswith("old_")],
           "old_S22": [k for k in names if k.startswith("old_")]}
    for b, pre in NEW_BLOCKS.items():
        out["new_" + b] = [k for k in names if k.startswith(pre) and not k.startswith("old_")]
    return out


def route_m(X, V, folds, rng, n_perm):
    """Routers A and B for the m ladder.  V (n, 15) per-m point-cloud RMSD (a label matrix).
    Returns per-target routed endpoints and the permutation null of the routed mean."""
    n = len(V); fixed = V[:, MS.index(75)]
    predV, aA = PS.nested_predict(X, V, folds, "reg")
    mA = np.argmin(predV, 1); endA = V[np.arange(n), mA]
    logm = np.log(np.array(MS))[np.argmin(V, 1)]
    predL, aB = PS.nested_predict(X, logm, folds, "reg")
    mB = np.abs(np.log(np.array(MS))[None, :] - predL[:, None]).argmin(1); endB = V[np.arange(n), mB]
    nullA, nullB = [], []
    for _ in range(n_perm):
        perm = rng.permutation(n); Vp = V[perm]
        pv, _ = PS.nested_predict(X, Vp, folds, "reg"); nullA.append(float((Vp[np.arange(n), np.argmin(pv, 1)] - Vp[:, MS.index(75)]).mean()))
        lp, _ = PS.nested_predict(X, np.log(np.array(MS))[np.argmin(Vp, 1)], folds, "reg")
        mb = np.abs(np.log(np.array(MS))[None, :] - lp[:, None]).argmin(1); nullB.append(float((Vp[np.arange(n), mb] - Vp[:, MS.index(75)]).mean()))
    return {"A": {"end": endA, "m": mA, "alphas": aA, "null": np.array(nullA)},
            "B": {"end": endB, "m": mB, "alphas": aB, "null": np.array(nullB)}, "fixed": fixed}


def labels_m():
    z = json.load(open(os.path.join(ROOT, "s12", "results", "agg_surface.json")))
    pdbs = [t["pdb"] for t in I.targets()]
    V = np.array([[z[p]["shipped|m%d" % m] for m in MS] for p in pdbs], float)
    return pdbs, V


def labels_s():
    z = json.load(open(os.path.join(ROOT, "s23", "results", "errdecomp.json")))
    return {r["pdb"]: r for r in z["rows"]}


_PREP = {}
_CHAIN = {}


def prep(pdb):
    """Per-target pool in shipped-score order, the native (ORACLE, read once) and the m=75 cloud."""
    if pdb not in _PREP:
        t = {x["pdb"]: x for x in I.targets()}[pdb]
        u = I.load_univ(pdb); idx = I.pool_idx(u); dg = I.distogram(pdb)
        W = np.asarray(u["W"][idx], float); i, j = I.pair_index(int(u["n"]))
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(W, i, j)), float)
        o = np.argsort(sc, kind="stable")
        _PREP[pdb] = {"W": W[o], "nat": np.asarray(u["nat_ca"], float), "seq": t["seq"], "fold": int(t["fold"]),
                      "C75": None}
        _PREP[pdb]["C75"] = cloud_m(pdb, 75)
        del u
    return _PREP[pdb]


def cloud_m(pdb, m):
    """The top-m (by the shipped score) average in the subset's own medoid frame (agg_surface's operator)."""
    p = _PREP[pdb] if pdb in _PREP else prep(pdb)
    top = p["W"][:m]
    if m == 1:
        return top[0]
    P = I.pairwise_rmsd(top)
    return I.superpose_batch(top, top[I.medoid(P)]).mean(0)


def scale_endpoint(pdb, s):
    """RMSD of the shipped m=75 cloud scaled by s about its centroid (errdecomp's construction).  ORACLE read."""
    p = prep(pdb); C = p["C75"]; nat = p["nat"]
    return float(I.ca_rmsd((C - C.mean(0)) * s, nat - nat.mean(0)))


def chain_endpoint(pdb, m):
    """Built chain of the top-m average through the production projection.  Cached per (pdb, m)."""
    key = (pdb, int(m))
    if key not in _CHAIN:
        p = prep(pdb); pr = I.project(cloud_m(pdb, int(m)), p["seq"], p["fold"])
        _CHAIN[key] = float(I.ca_rmsd(pr["ca"], p["nat"]))
    return _CHAIN[key]


def run(n_perm=N_PERM):
    if not L._signed_off():
        raise SystemExit("PHASE GATE: p_c4 run reads persisted RMSDs; refused before sign-off.")
    fz = build_features(); rows = fz["rows"]; names = fz["names"]
    pdbs = [r["pdb"] for r in rows]; folds = np.array([r["fold"] for r in rows])
    pd2, V = labels_m(); assert pd2 == pdbs
    sst = labels_s(); s_star = np.array([sst[p]["s_star"] for p in pdbs]); rmsd1 = np.array([sst[p]["rmsd_1"] for p in pdbs])
    rng = np.random.default_rng(26)
    out = {"n": len(pdbs), "names": names, "MS": MS, "m_router": {}, "s_router": {}}
    fixed = V[:, MS.index(75)]
    bok = ST.best_of_k_within(V, n_boot=300)
    out["m_oracle"] = {"gain": float((V.min(1) - fixed).mean()), "best_of_k_within": {k: v for k, v in bok.items() if k != "argmin_counts"}}
    # anchors: this module's own m=75 cloud vs the persisted surface and errdecomp (same operator, re-derived)
    for p in pdbs:
        prep(p)
    c75 = np.array([I.ca_rmsd(_PREP[p]["C75"], _PREP[p]["nat"]) for p in pdbs])
    out["anchor"] = {"m75_cloud_vs_surface_maxabs": float(np.abs(c75 - fixed).max()),
                     "m75_cloud_vs_errdecomp_rmsd1_maxabs": float(np.abs(c75 - rmsd1).max())}
    print("anchor: rebuilt m=75 cloud vs agg_surface max abs %.2e ; vs errdecomp rmsd_1 %.2e" % (out["anchor"]["m75_cloud_vs_surface_maxabs"], out["anchor"]["m75_cloud_vs_errdecomp_rmsd1_maxabs"]))
    chain75 = np.array([chain_endpoint(p, 75) for p in pdbs])
    out["anchor"]["m75_chain_mean"] = float(chain75.mean())
    print("anchor: rebuilt m=75 built chain mean %.4f (production rmsd_arm 3.2148)" % chain75.mean())
    for bname, cols in blocks(names).items():
        X = np.array([[r[k] for k in cols] for r in rows], float)
        r = route_m(X, V, folds, rng, n_perm)
        rep = {}
        for key in ("A", "B"):
            cmp = ST.compare(r[key]["end"], fixed, folds, names=pdbs, label="C4 m-router %s [%s] routed - fixed m=75 (POINT CLOUD, persisted surface)" % (key, bname))
            rep[key] = {"stats": cmp, "null_mean": float(r[key]["null"].mean()), "null_p05": float(np.percentile(r[key]["null"], 5)),
                        "p_perm": float((r[key]["null"] <= cmp["effect"]).mean()), "alphas": r[key]["alphas"],
                        "m_hist": {int(m): int((np.array(MS)[r[key]["m"]] == m).sum()) for m in MS}}
            print(ST.fmt(cmp)); print("    perm null mean %+.4f p05 %+.4f  p_perm %.3f" % (rep[key]["null_mean"], rep[key]["null_p05"], rep[key]["p_perm"]))
            if bname in ("new_all", "old_S22"):
                routed_chain = np.array([chain_endpoint(p, MS[int(mi)]) for p, mi in zip(pdbs, r[key]["m"])])
                cmpc = ST.compare(routed_chain, chain75, folds, names=pdbs, label="C4 m-router %s [%s] routed - fixed m=75 (BUILT CHAIN, same projection both sides)" % (key, bname))
                rep[key]["stats_chain"] = cmpc; print(ST.fmt(cmpc))
        out["m_router"][bname] = rep
        ps_, a_s = PS.nested_predict(X, s_star, folds, "reg")
        routed = np.array([scale_endpoint(p, float(ps_[k])) for k, p in enumerate(pdbs)])
        cmp = ST.compare(routed, rmsd1, folds, names=pdbs, label="C4 s-router [%s] routed - s=1 (point cloud)" % bname)
        null = []
        for _ in range(min(n_perm, 100)):
            perm = rng.permutation(len(pdbs)); pp, _ = PS.nested_predict(X, s_star[perm], folds, "reg")
            null.append(float(np.mean([scale_endpoint(p, float(pp[k])) for k, p in enumerate(pdbs)]) - rmsd1.mean()))
        out["s_router"][bname] = {"stats": cmp, "rho_pred_vs_sstar": float(np.corrcoef(ps_, s_star)[0, 1]),
                                  "null_mean": float(np.mean(null)), "null_p05": float(np.percentile(null, 5)), "alphas": a_s}
        print(ST.fmt(cmp)); print("    rho(pred, s*) %+.3f ; perm null mean %+.4f" % (out["s_router"][bname]["rho_pred_vs_sstar"], out["s_router"][bname]["null_mean"]))
        ST.save_atomic(OUT, out, module_file=__file__)
    return out


def selftest():
    """Synthetic labels on the real features: a planted per-m structure is routed, a random one is not."""
    fz = build_features(verbose=False); rows = fz["rows"]; names = fz["names"]
    folds = np.array([r["fold"] for r in rows]); n = len(rows)
    cols = blocks(names)["new_all"]; X = np.array([[r[k] for k in cols] for r in rows], float)
    rng = np.random.default_rng(3)
    Z = (X - X.mean(0)) / (X.std(0) + 1e-9)
    t = Z[:, cols.index("pax_f1")] + 0.5 * Z[:, cols.index("sim_ent")]          # a planted latent
    best = np.clip(np.round(7 + 3 * t), 0, len(MS) - 1).astype(int)              # planted best rung
    V = 3.0 + 0.15 * np.abs(np.arange(len(MS))[None, :] - best[:, None]) + 0.05 * rng.normal(size=(n, len(MS)))
    r = route_m(X, V, folds, rng, n_perm=10)
    gA = float((r["A"]["end"] - r["fixed"]).mean()); gB = float((r["B"]["end"] - r["fixed"]).mean())
    print("  planted: routed A %+.3f  B %+.3f  (null A mean %+.3f)  oracle %+.3f" % (gA, gB, r["A"]["null"].mean(), (V.min(1) - r["fixed"]).mean()))
    assert gA < r["A"]["null"].mean() - 0.05 and gB < r["B"]["null"].mean() - 0.05
    Vr = 3.0 + 0.3 * rng.normal(size=(n, len(MS)))
    r = route_m(X, Vr, folds, rng, n_perm=10)
    gA = float((r["A"]["end"] - r["fixed"]).mean())
    print("  random:  routed A %+.3f  null mean %+.3f p05 %+.3f" % (gA, r["A"]["null"].mean(), np.percentile(r["A"]["null"], 5)))
    assert gA >= np.percentile(r["A"]["null"], 5) - 0.05
    print("  p_c4 selftest OK (%d new features, %d old)" % (len(cols), len(blocks(names)["old_S22"])))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "features"
    {"features": lambda: build_features(), "selftest": selftest, "run": run}[cmd]()
