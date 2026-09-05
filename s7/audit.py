"""Adversarial audit of the claim the sprint's forward plan rests on.

THE CLAIM (s6/FINDINGS.md section 13, commit 256539f):

    On 24 cluster-disjoint held-out peptide targets, ranking a 500-candidate pool of real
    protein fragments by agreement with the NATIVE's CA-CA distance matrix gives in-band
    Spearman +0.641 and selects structures averaging 1.836 A against a pool best of
    1.534 A. The learned sequence-to-distance prior gives in-band +0.013 / 3.398 A.
    Inference: "the distance channel is sufficient; the predictor delivers 2% of it."

Nothing here imports the code that produced that number. Pools, RMSD and the distance
comparison are re-derived from `peptide_db`, `fragment_db` and `distogram._fold_fragments`
only; `s5.lib.kabsch_rmsd_batch` is deliberately NOT imported and the local Kabsch is
cross-checked against `protein_geometry.ca_rmsd` (an independent implementation that forms
the rotation explicitly) plus a mirror-image control.

Stages, each resumable and each writing atomically:

    python -m s7.audit pools     # build the 24 pools, cache RMSDs + distance matrices
    python -m s7.audit verify    # RMSD cross-check, mirror control, tie-break sensitivity
    python -m s7.audit leak      # leakage sweep, counts not assurances
    python -m s7.audit analyse   # oracle reproduction, bands, pool size, splits
    python -m s7.audit noise     # THE key experiment: degrade the oracle, watch selection
    python -m s7.audit report    # print everything already computed
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.stats import spearmanr

os.environ.setdefault("NT", "2")

import distogram as dgm                                                    # noqa: E402
import fragment_db as fdb                                                  # noqa: E402
import peptide_db as db                                                    # noqa: E402
import protein_geometry as geo                                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "audit_cache")
KMAX = 2000                      #: cached candidates; K=500 is the claim's pool
K_CLAIM = 500
BANDS = (1.0, 1.5, 2.0, 3.0)
POOL_SIZES = (100, 500, 2000)
SIGMAS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
N_NOISE_SEEDS = 20

AA = "ARNDCQEGHILKMFPSTWYV"
_IDX = {a: i for i, a in enumerate(AA)}

#: BLOSUM62 in the AA order above -- retyped here rather than imported, so a pool built by
#: this file is not the same object the claim's pool was built from.
_B62_ROWS = """4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
4 2 -3 1 0 -3 -2 -1 -3 -1 3
4 -2 2 0 -3 -2 -1 -2 -1 1
5 -1 -3 -1 0 -1 -3 -2 -2
5 0 -2 -1 -1 -1 -1 1
6 -4 -2 -2 1 3 -1
7 -1 -1 -4 -3 -2
4 1 -3 -2 -2
5 -2 -2 0
11 2 -3
7 -1
4"""


def blosum62():
    M = np.zeros((20, 20))
    for i, row in enumerate(r.split() for r in _B62_ROWS.strip().split("\n")):
        for k, v in enumerate(row):
            M[i, i + k] = M[i + k, i] = float(v)
    return M


B62 = blosum62()


def encode(seq):
    return np.array([_IDX.get(c, 0) for c in seq], dtype=np.int64)


# --------------------------------------------------------------------------- geometry
def kabsch_rmsd_batch(P, ref):
    """CA-RMSD of every ``(n,3)`` structure in ``P`` against ``ref``, written from scratch.

    Optimal superposition without forming the rotation: with H = P_c^T ref_c and singular
    values s, the residual is |P_c|^2 + |ref_c|^2 - 2 sum(s), and flipping the sign of the
    smallest singular value when det(V U^T) < 0 is what forbids a reflection -- the step
    that makes a mirror image score as different rather than identical.
    """
    P = np.asarray(P, float)
    ref = np.asarray(ref, float)
    if P.ndim == 2:
        P = P[None]
    Pc = P - P.mean(1, keepdims=True)
    Rc = ref - ref.mean(0, keepdims=True)
    H = np.einsum("bni,nj->bij", Pc, Rc)
    U, S, Vt = np.linalg.svd(H)
    det = np.linalg.det(np.einsum("bji,bkj->bik", Vt, U))     # det(V U^T), per structure
    S = S.copy()
    S[:, -1] *= np.sign(det)
    resid = (Pc ** 2).sum((1, 2)) + (Rc ** 2).sum() - 2.0 * S.sum(1)
    return np.sqrt(np.maximum(resid, 0.0) / P.shape[1])


def pair_index(n, min_sep=2):
    return np.triu_indices(n, k=min_sep)


def pair_dists(W, i, j):
    """(B, npairs) CA-CA distances for the cached pair set."""
    return np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1)


# --------------------------------------------------------------------------- pools
def windows_of(pool, n):
    """Every length-n CA window of every pool member, with provenance."""
    cas, seqs, src = [], [], []
    for q in pool:
        m = len(q.seq)
        if m < n:
            continue
        e = encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n])
            seqs.append(e[s:s + n])
            src.append(q.pdb)
    return np.stack(cas), np.stack(seqs), np.array(src, dtype=object)


def build_pool_members(target, folds, n_folds=5):
    """The leakage-safe library for one target: out-of-fold peptides + fold fragments."""
    fold = folds[target.seq]
    peps = [q for q in db.load() if folds[q.seq] != fold and q.seq != target.seq]
    frags = list(dgm._fold_fragments(fold, n_folds))
    return fold, peps, frags


def build_target(target, folds, kmax=KMAX, tie_seed=None):
    """Cacheable record for one target. `tie_seed` shuffles before the stable argsort,
    which is how the arbitrariness of the BLOSUM tie set is probed."""
    fold, peps, frags = build_pool_members(target, folds)
    W, S, src = windows_of(peps + frags, target.n)
    sim = B62[S, encode(target.seq)[None, :]].sum(1)
    if tie_seed is None:
        order = np.argsort(-sim, kind="stable")
    else:
        perm = np.random.default_rng(tie_seed).permutation(len(sim))
        order = perm[np.argsort(-sim[perm], kind="stable")]
    idx = order[:kmax]
    i, j = pair_index(target.n)
    rec = {
        "pdb": target.pdb, "n": int(target.n), "fold": int(fold), "seq": target.seq,
        "n_windows": int(len(W)), "n_peptides": len(peps), "n_fragments": len(frags),
        "rr": kabsch_rmsd_batch(W[idx], target.ca).astype(np.float32),
        "D": pair_dists(W[idx], i, j).astype(np.float32),
        "Dnat": pair_dists(target.ca[None], i, j)[0].astype(np.float32),
        "sim": sim[idx].astype(np.float32),
        "src": src[idx],
        "sim_cut": float(sim[order[kmax - 1]]),
        "n_tied_at_cut": int((sim == sim[order[K_CLAIM - 1]]).sum()),
        "n_ge_cut": int((sim >= sim[order[K_CLAIM - 1]]).sum()),
    }
    return rec


def cache_path(pdbid, tag=""):
    return os.path.join(CACHE, f"{pdbid}{tag}.npz")


def save_rec(rec, tag=""):
    os.makedirs(CACHE, exist_ok=True)
    p = cache_path(rec["pdb"], tag)
    tmp = p + ".tmp.npz"
    np.savez_compressed(tmp, **{k: np.asarray(v) for k, v in rec.items()})
    os.replace(tmp, p)


def load_rec(pdbid, tag=""):
    z = np.load(cache_path(pdbid, tag), allow_pickle=True)
    out = {k: z[k] for k in z.files}
    for k in ("pdb", "seq"):
        out[k] = str(out[k])
    for k in ("n", "fold", "n_windows", "n_peptides", "n_fragments",
              "n_tied_at_cut", "n_ge_cut"):
        out[k] = int(out[k])
    return out


def dev_targets():
    return db.dev_set(24)


def stage_pools():
    folds = db.folds(5)
    for p in dev_targets():
        if os.path.exists(cache_path(p.pdb)):
            continue
        t0 = time.time()
        rec = build_target(p, folds)
        save_rec(rec)
        print(f"{rec['pdb']:6} n={rec['n']:2d} fold={rec['fold']} "
              f"windows={rec['n_windows']:6d} best@500={rec['rr'][:500].min():.3f} "
              f"tied_at_cut={rec['n_tied_at_cut']:5d} ({time.time()-t0:.0f}s)", flush=True)


# --------------------------------------------------------------------------- rankers
def score_l1(D, Dpred):
    """Mean absolute distance-matrix disagreement -- the oracle's ranking statistic."""
    return np.abs(D - Dpred[None, :]).mean(1)


def rank_stats(score, rr, bands=BANDS):
    """Global rho, per-band rho, and the RMSD of the argmin. Lower score = better, so a
    POSITIVE rho means the score agrees with RMSD."""
    out = {"sel": float(rr[int(np.argmin(score))]),
           "pool": float(rr.min()),
           "rho_global": float(spearmanr(score, rr).statistic)}
    for b in bands:
        m = rr <= rr.min() + b
        out[f"rho_band{b}"] = (float(spearmanr(score[m], rr[m]).statistic)
                               if m.sum() > 4 else None)
        out[f"n_band{b}"] = int(m.sum())
        out[f"sel_band{b}"] = float(rr[m][int(np.argmin(score[m]))]) if m.sum() else None
    return out


def mean_of(rows, key):
    v = [r[key] for r in rows if r.get(key) is not None]
    return float(np.mean(v)) if v else None


# --------------------------------------------------------------------------- analyse
def stage_analyse():
    recs = [load_rec(p.pdb) for p in dev_targets()]
    out = {"per_target": [], "summary": {}}
    for rec in recs:
        D, Dn, rr = rec["D"], rec["Dnat"], rec["rr"]
        row = {"pdb": rec["pdb"], "n": rec["n"], "fold": rec["fold"],
               "n_windows": rec["n_windows"], "n_tied_at_cut": rec["n_tied_at_cut"]}
        for K in POOL_SIZES:
            d, r = D[:K], rr[:K]
            st = rank_stats(score_l1(d, Dn), r)
            for k, v in st.items():
                row[f"K{K}_{k}"] = v
            # is the distance channel just RMSD in disguise?
            s = score_l1(d, Dn)
            row[f"K{K}_pearson_score_rmsd"] = float(np.corrcoef(s, r)[0, 1])
            # dRMSD, the root-mean-square form of the same statistic
            drms = np.sqrt(((d - Dn[None, :]) ** 2).mean(1))
            row[f"K{K}_rho_drmsd_rmsd"] = float(spearmanr(drms, r).statistic)
            row[f"K{K}_sel_drmsd"] = float(r[int(np.argmin(drms))])
            # controls
            row[f"K{K}_sel_random"] = float(np.mean(r))
            mean_mat = d.mean(0)
            row[f"K{K}_sel_poolmean"] = float(r[int(np.argmin(score_l1(d, mean_mat)))])
            row[f"K{K}_rho_poolmean"] = float(spearmanr(score_l1(d, mean_mat), r).statistic)
        out["per_target"].append(row)

    keys = [k for k in out["per_target"][0] if k.startswith("K")]
    for k in keys:
        out["summary"][k] = mean_of(out["per_target"], k)

    # how well behaved is the in-band statistic itself? Its membership is a function of
    # the pool, its size varies by two orders of magnitude across targets, and targets
    # with fewer than five in-band candidates are silently dropped from the mean.
    rows0 = out["per_target"]
    band_diag = {}
    for b in BANDS:
        n = np.array([r[f"K500_n_band{b}"] for r in rows0], float)
        rho = np.array([r[f"K500_rho_band{b}"] if r[f"K500_rho_band{b}"] is not None
                        else np.nan for r in rows0], float)
        ok = ~np.isnan(rho)
        band_diag[f"band{b}"] = {
            "n_targets_scored": int(ok.sum()),
            "n_targets_dropped": int((~ok).sum()),
            "dropped": [r["pdb"] for r, o in zip(rows0, ok) if not o],
            "band_size_min": float(n.min()), "band_size_median": float(np.median(n)),
            "band_size_max": float(n.max()),
            "rho_unweighted": float(np.nanmean(rho)),
            "rho_size_weighted": float(np.nansum(rho[ok] * n[ok]) / n[ok].sum()),
            "rho_min": float(np.nanmin(rho)), "rho_max": float(np.nanmax(rho)),
            "n_targets_negative": int((rho[ok] < 0).sum()),
            "corr_bandsize_rho": float(spearmanr(n[ok], rho[ok]).statistic)}
    out["band_diagnostics"] = band_diag

    # wrong-native control: rank with another target's native matrix, same length
    ctrl = []
    for a in recs:
        alts = [b for b in recs if b["n"] == a["n"] and b["pdb"] != a["pdb"]]
        if not alts:
            continue
        vals = [rank_stats(score_l1(a["D"][:K_CLAIM], b["Dnat"]), a["rr"][:K_CLAIM])
                for b in alts]
        ctrl.append({"pdb": a["pdb"], "n_alts": len(alts),
                     "sel": float(np.mean([v["sel"] for v in vals])),
                     "rho_global": float(np.mean([v["rho_global"] for v in vals])),
                     "rho_band1.5": float(np.mean([v["rho_band1.5"] for v in vals
                                                   if v["rho_band1.5"] is not None]
                                                  or [np.nan]))})
    out["wrong_native_control"] = {
        "per_target": ctrl,
        "sel": float(np.mean([c["sel"] for c in ctrl])),
        "rho_global": float(np.mean([c["rho_global"] for c in ctrl])),
        "rho_band1.5": float(np.nanmean([c["rho_band1.5"] for c in ctrl])),
    }

    # splits: chain length and pool difficulty
    rows = out["per_target"]
    out["splits"] = {}
    med_n = float(np.median([r["n"] for r in rows]))
    for name, sel in (("short", lambda r: r["n"] <= med_n),
                      ("long", lambda r: r["n"] > med_n)):
        g = [r for r in rows if sel(r)]
        out["splits"][f"len_{name}"] = {
            "n_targets": len(g), "mean_n": float(np.mean([r["n"] for r in g])),
            "pool": mean_of(g, "K500_pool"), "sel": mean_of(g, "K500_sel"),
            "rho_band1.5": mean_of(g, "K500_rho_band1.5"),
            "rho_global": mean_of(g, "K500_rho_global"),
            "gap": mean_of(g, "K500_sel") - mean_of(g, "K500_pool")}
    q = np.quantile([r["K500_pool"] for r in rows], [1 / 3, 2 / 3])
    for name, sel in (("easy", lambda r: r["K500_pool"] <= q[0]),
                      ("mid", lambda r: q[0] < r["K500_pool"] <= q[1]),
                      ("hard", lambda r: r["K500_pool"] > q[1])):
        g = [r for r in rows if sel(r)]
        out["splits"][f"pool_{name}"] = {
            "n_targets": len(g), "pool": mean_of(g, "K500_pool"),
            "sel": mean_of(g, "K500_sel"), "rho_band1.5": mean_of(g, "K500_rho_band1.5"),
            "gap": mean_of(g, "K500_sel") - mean_of(g, "K500_pool")}
    write_json("audit_analyse.json", out)
    return out


# --------------------------------------------------------------------------- noise
def stage_noise():
    """Degrade the oracle and watch selection. Two error models, because the answer
    depends on which one a real predictor makes:

      `iid`    Gaussian noise added independently to each predicted pair distance. This
               is the requested sweep. Its MAE is sigma*sqrt(2/pi).
      `coord`  Gaussian displacement of the native CA coordinates, then the distance
               matrix of the displaced structure. Errors are then structurally consistent
               (they come from a real, if wrong, conformation) and correlated across pairs,
               which is what a structure predictor's errors actually look like.
    """
    recs = [load_rec(p.pdb) for p in dev_targets()]
    folds = db.folds(5)
    natives = {p.pdb: p.ca for p in dev_targets()}
    out = {"iid": [], "coord": [], "per_target_iid": []}
    for sigma in SIGMAS:
        agg = {"sel": [], "rho_g": [], "rho_b": [], "mae": [], "sel_b": []}
        per_t = {r["pdb"]: [] for r in recs}
        for rec in recs:
            D, Dn, rr = rec["D"][:K_CLAIM], rec["Dnat"], rec["rr"][:K_CLAIM]
            rng = np.random.default_rng(abs(hash(rec["pdb"])) % (2 ** 31))
            for s in range(N_NOISE_SEEDS if sigma > 0 else 1):
                Dp = Dn + rng.normal(0, sigma, size=Dn.shape) if sigma > 0 else Dn
                st = rank_stats(score_l1(D, Dp), rr)
                agg["sel"].append(st["sel"])
                agg["rho_g"].append(st["rho_global"])
                if st["rho_band1.5"] is not None:
                    agg["rho_b"].append(st["rho_band1.5"])
                if st["sel_band1.5"] is not None:
                    agg["sel_b"].append(st["sel_band1.5"])
                agg["mae"].append(float(np.abs(Dp - Dn).mean()))
                per_t[rec["pdb"]].append(st["sel"])
        out["iid"].append({
            "sigma": sigma, "mae": float(np.mean(agg["mae"])),
            "sel": float(np.mean(agg["sel"])),
            "sel_sd": float(np.std(agg["sel"])),
            "sel_band1.5": float(np.mean(agg["sel_b"])),
            "rho_global": float(np.mean(agg["rho_g"])),
            "rho_band1.5": float(np.mean(agg["rho_b"]))})
        out["per_target_iid"].append(
            {"sigma": sigma,
             **{k: float(np.mean(v)) for k, v in per_t.items()}})
        print(f"iid   sigma {sigma:4.2f}  MAE {out['iid'][-1]['mae']:5.3f}  "
              f"sel {out['iid'][-1]['sel']:5.3f}  "
              f"rho_g {out['iid'][-1]['rho_global']:+.3f}  "
              f"rho_band1.5 {out['iid'][-1]['rho_band1.5']:+.3f}", flush=True)

    for tau in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        agg = {"sel": [], "rho_g": [], "rho_b": [], "mae": [], "rmsd": []}
        for rec in recs:
            D, Dn, rr = rec["D"][:K_CLAIM], rec["Dnat"], rec["rr"][:K_CLAIM]
            ca = natives[rec["pdb"]]
            i, j = pair_index(rec["n"])
            rng = np.random.default_rng(1000 + abs(hash(rec["pdb"])) % (2 ** 31))
            for s in range(N_NOISE_SEEDS if tau > 0 else 1):
                pert = ca + rng.normal(0, tau, size=ca.shape) if tau > 0 else ca
                Dp = pair_dists(pert[None], i, j)[0]
                st = rank_stats(score_l1(D, Dp), rr)
                agg["sel"].append(st["sel"])
                agg["rho_g"].append(st["rho_global"])
                if st["rho_band1.5"] is not None:
                    agg["rho_b"].append(st["rho_band1.5"])
                agg["mae"].append(float(np.abs(Dp - Dn).mean()))
                agg["rmsd"].append(float(kabsch_rmsd_batch(pert[None], ca)[0]))
        out["coord"].append({
            "tau": tau, "mae": float(np.mean(agg["mae"])),
            "self_rmsd": float(np.mean(agg["rmsd"])),
            "sel": float(np.mean(agg["sel"])),
            "sel_sd": float(np.std(agg["sel"])),
            "rho_global": float(np.mean(agg["rho_g"])),
            "rho_band1.5": float(np.mean(agg["rho_b"]))})
        print(f"coord tau {tau:4.2f}  MAE {out['coord'][-1]['mae']:5.3f}  "
              f"selfRMSD {out['coord'][-1]['self_rmsd']:5.3f}  "
              f"sel {out['coord'][-1]['sel']:5.3f}  "
              f"rho_band1.5 {out['coord'][-1]['rho_band1.5']:+.3f}", flush=True)
    write_json("audit_noise.json", out)
    return out


# --------------------------------------------------------------------------- prior
def stage_prior():
    """Where does a REAL predictor sit on the noise curve, and is iid noise the right
    error model for it?

    The learned prior cannot be evaluated here without importing `esm_features`, whose
    cache is 1.5 GB and which this run is forbidden to load. So two things are done
    instead, neither of which needs the model:

      1. a THIRD error model -- shrinkage of the native matrix toward the pool's own mean
         distance profile, `D = (1-a) Dnat + a Dbar`. That is what an underfit regressor
         produces: not noise around the truth but regression toward the average peptide.
         Same MAE axis, structured and biased error.
      2. the recorded per-target numbers for the shipped distogram (`s5/inband.json`,
         data only -- no s5 code is imported) placed on those curves at their own MAE.
    """
    recs = [load_rec(p.pdb) for p in dev_targets()]
    out = {"shrink": [], "prior_vs_curve": []}
    for a in (0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0):
        sel, mae, rg, rb, per = [], [], [], [], {}
        for rec in recs:
            D, Dn, rr = rec["D"][:K_CLAIM], rec["Dnat"], rec["rr"][:K_CLAIM]
            Dbar = D.mean(0)
            Dp = (1 - a) * Dn + a * Dbar
            st = rank_stats(score_l1(D, Dp), rr)
            sel.append(st["sel"]); rg.append(st["rho_global"])
            if st["rho_band1.5"] is not None:
                rb.append(st["rho_band1.5"])
            mae.append(float(np.abs(Dp - Dn).mean()))
            per[rec["pdb"]] = {"mae": mae[-1], "sel": st["sel"]}
        out["shrink"].append({"a": a, "mae": float(np.mean(mae)),
                              "sel": float(np.mean(sel)),
                              "rho_global": float(np.mean(rg)),
                              "rho_band1.5": float(np.mean(rb)), "per_target": per})
        print(f"shrink a={a:4.2f}  MAE {out['shrink'][-1]['mae']:5.3f}  "
              f"sel {out['shrink'][-1]['sel']:5.3f}  "
              f"rho_g {out['shrink'][-1]['rho_global']:+.3f}  "
              f"rho_band1.5 {out['shrink'][-1]['rho_band1.5']:+.3f}", flush=True)

    inband = os.path.join(os.path.dirname(HERE), "s5", "inband.json")
    if os.path.exists(inband):
        with open(inband) as f:
            rows = json.load(f)
        nz = read_json("audit_noise.json")
        iid_mae = np.array([r["mae"] for r in nz["iid"]])
        per_iid = nz["per_target_iid"]
        for r in rows:
            pid = r["pdb"]
            curve = np.array([p[pid] for p in per_iid])
            shrink_mae = np.array([s["per_target"][pid]["mae"] for s in out["shrink"]])
            shrink_sel = np.array([s["per_target"][pid]["sel"] for s in out["shrink"]])
            o = np.argsort(shrink_mae)
            out["prior_vs_curve"].append({
                "pdb": pid, "prior_mae": r["mae"], "prior_sel": r["sel"],
                "prior_rho_global": r["rho_global"],
                "prior_rho_band1.5": r["rho_band1.5"],
                "iid_sel_at_prior_mae": float(np.interp(r["mae"], iid_mae, curve)),
                "shrink_sel_at_prior_mae": float(np.interp(r["mae"], shrink_mae[o],
                                                           shrink_sel[o])),
                "oracle_sel": float(curve[0])})
        pv = out["prior_vs_curve"]
        out["prior_summary"] = {
            "n": len(pv),
            "prior_mae": float(np.mean([p["prior_mae"] for p in pv])),
            "prior_sel": float(np.mean([p["prior_sel"] for p in pv])),
            "iid_sel_at_prior_mae": float(np.mean([p["iid_sel_at_prior_mae"]
                                                   for p in pv])),
            "shrink_sel_at_prior_mae": float(np.mean([p["shrink_sel_at_prior_mae"]
                                                      for p in pv])),
            "oracle_sel": float(np.mean([p["oracle_sel"] for p in pv])),
            "n_prior_worse_than_iid": int(sum(p["prior_sel"] > p["iid_sel_at_prior_mae"]
                                              for p in pv))}
    write_json("audit_prior.json", out)
    return out


# --------------------------------------------------------------------------- verify
def stage_verify():
    """Independent RMSD verification, mirror control, and pool tie-break sensitivity."""
    out = {}
    folds = db.folds(5)
    targets = dev_targets()
    p = targets[0]
    fold, peps, frags = build_pool_members(p, folds)
    W, S, src = windows_of((peps + frags)[:400], p.n)
    W = W[:300]
    mine = kabsch_rmsd_batch(W, p.ca)
    theirs = np.array([geo.ca_rmsd(w, p.ca) for w in W])
    out["kabsch_vs_protein_geometry"] = {
        "n": int(len(W)), "target": p.pdb,
        "max_abs_diff": float(np.abs(mine - theirs).max()),
        "mean_abs_diff": float(np.abs(mine - theirs).mean()),
        "mine_first5": [float(x) for x in mine[:5]],
        "theirs_first5": [float(x) for x in theirs[:5]]}

    # mirror control: reflecting z must NOT superpose onto the original
    mir = p.ca.copy()
    mir[:, 2] *= -1.0
    out["mirror_control"] = {
        "target": p.pdb,
        "self_rmsd_mine": float(kabsch_rmsd_batch(p.ca[None], p.ca)[0]),
        "mirror_rmsd_mine": float(kabsch_rmsd_batch(mir[None], p.ca)[0]),
        "mirror_rmsd_geo": float(geo.ca_rmsd(mir, p.ca)),
        "mirror_distance_matrix_max_diff": float(
            np.abs(pair_dists(mir[None], *pair_index(p.n))[0]
                   - pair_dists(p.ca[None], *pair_index(p.n))[0]).max())}

    # the distance channel is blind to chirality: what does that cost?
    mirror_rows = []
    for t in targets:
        rec = load_rec(t.pdb)
        D, Dn, rr = rec["D"][:K_CLAIM], rec["Dnat"], rec["rr"][:K_CLAIM]
        ca = t.ca.copy()
        ca[:, 2] *= -1.0
        i, j = pair_index(t.n)
        # the mirrored native as a candidate: perfect distance score, real RMSD
        Dmir = pair_dists(ca[None], i, j)[0]
        mirror_rows.append({
            "pdb": t.pdb,
            "score_of_mirror": float(score_l1(Dmir[None], Dn)[0]),
            "best_pool_score": float(score_l1(D, Dn).min()),
            "rmsd_of_mirror": float(kabsch_rmsd_batch(ca[None], t.ca)[0])})
    out["chirality_blindness"] = {
        "per_target": mirror_rows,
        "mean_mirror_rmsd": float(np.mean([m["rmsd_of_mirror"] for m in mirror_rows])),
        "n_mirror_beats_pool": int(sum(m["score_of_mirror"] < m["best_pool_score"]
                                       for m in mirror_rows))}

    # tie-break sensitivity: rebuild every pool with a shuffled tie order
    tie = []
    for t in targets:
        base = load_rec(t.pdb)
        alt = build_target(t, folds, tie_seed=17)
        for rec, tag in ((base, "stable"), (alt, "shuffled")):
            st = rank_stats(score_l1(rec["D"][:K_CLAIM], rec["Dnat"]),
                            rec["rr"][:K_CLAIM])
            rec["_st"] = st
        overlap = len(set(map(tuple, np.round(base["D"][:K_CLAIM], 4)))
                      & set(map(tuple, np.round(alt["D"][:K_CLAIM], 4))))
        tie.append({"pdb": t.pdb,
                    "n_tied_at_cut": base["n_tied_at_cut"],
                    "n_ge_cut": base["n_ge_cut"],
                    "overlap_500": overlap,
                    "pool_stable": base["_st"]["pool"], "pool_shuffled": alt["_st"]["pool"],
                    "sel_stable": base["_st"]["sel"], "sel_shuffled": alt["_st"]["sel"],
                    "band_stable": base["_st"]["rho_band1.5"],
                    "band_shuffled": alt["_st"]["rho_band1.5"]})
        print(f"tie {t.pdb:6} overlap {overlap:4d}/500  pool {base['_st']['pool']:.3f}"
              f"->{alt['_st']['pool']:.3f}  sel {base['_st']['sel']:.3f}"
              f"->{alt['_st']['sel']:.3f}", flush=True)
    out["tie_break"] = {
        "per_target": tie,
        "mean_overlap": float(np.mean([t["overlap_500"] for t in tie])),
        "pool_stable": float(np.mean([t["pool_stable"] for t in tie])),
        "pool_shuffled": float(np.mean([t["pool_shuffled"] for t in tie])),
        "sel_stable": float(np.mean([t["sel_stable"] for t in tie])),
        "sel_shuffled": float(np.mean([t["sel_shuffled"] for t in tie])),
        "band_stable": float(np.nanmean([t["band_stable"] if t["band_stable"] is not None
                                         else np.nan for t in tie])),
        "band_shuffled": float(np.nanmean([t["band_shuffled"]
                                           if t["band_shuffled"] is not None
                                           else np.nan for t in tie]))}
    write_json("audit_verify.json", out)
    return out


# --------------------------------------------------------------------------- leak
def stage_leak():
    """Counts, not assurances."""
    folds = db.folds(5)
    targets = dev_targets()
    all_frags = list(fdb.load())
    out = {"per_target": [], "fragment_db": {
        "n_fragments": len(all_frags),
        "n_source_pdbs": len(set(f.pdb.rsplit("_", 1)[0] for f in all_frags))}}
    fold_filter = {}
    peptide_pdbs = {p.pdb.upper() for p in db.load()}
    out["fragment_db"]["n_source_pdbs_shared_with_peptide_db"] = len(
        {f.pdb.rsplit("_", 1)[0].upper() for f in all_frags} & peptide_pdbs)
    for t in targets:
        fold = folds[t.seq]
        peps = [q for q in db.load() if folds[q.seq] != fold and q.seq != t.seq]
        frags = list(dgm._fold_fragments(fold, 5))
        rec = load_rec(t.pdb)
        # 1. is the target itself in the library?
        self_in_pep = sum(q.seq == t.seq for q in peps)
        self_in_frag = sum(f.seq == t.seq for f in frags)
        # 2. does any library member share the target's source PDB id?
        tid = t.pdb.upper()
        src_pdbs = {f.pdb.rsplit("_", 1)[0].upper() for f in frags}
        share_pdb_frag = int(tid in src_pdbs)
        share_pdb_pep = sum(q.pdb.upper() == tid for q in peps)
        # 3. how many library members exceed the identity threshold to the target?
        id_p = np.array([db.identity(t.seq, q.seq) for q in peps])
        id_f = np.array([db.identity(t.seq, f.seq) for f in frags])
        id_pep = int((id_p >= db.IDENTITY_THRESHOLD).sum())
        id_frag = int((id_f >= db.IDENTITY_THRESHOLD).sum())
        # 4. does my fold filter reproduce _fold_fragments exactly? (computed per fold,
        #    exhaustively -- no k-mer prefilter, so it also tests the production
        #    prefilter's soundness rather than reimplementing it)
        if fold not in fold_filter:
            held = [p.seq for p in db.load() if folds[p.seq] == fold]
            # keyed by (id, sequence): fragment_db reuses one id for several window
            # lengths at the same start, so a set of bare ids silently merges entries
            # and would make this comparison weaker than it looks.
            mine = {(f.pdb, f.seq) for f in all_frags
                    if not any(db.identity(s, f.seq) >= db.IDENTITY_THRESHOLD
                               for s in held)}
            theirs = {(f.pdb, f.seq) for f in frags}
            fold_filter[fold] = (mine, theirs, len(held))
            print(f"  fold {fold}: exhaustive filter keeps {len(mine)}, "
                  f"_fold_fragments keeps {len(theirs)}, "
                  f"held peptides {len(held)}", flush=True)
        mine, theirs, n_held = fold_filter[fold]
        # 5. what actually made the top-500 pool?
        src500 = rec["src"][:K_CLAIM]
        row = {"pdb": t.pdb, "fold": int(fold), "n": t.n,
               "n_peptides": len(peps), "n_fragments": len(frags),
               "self_seq_in_peptides": int(self_in_pep),
               "self_seq_in_fragments": int(self_in_frag),
               "library_shares_target_pdb_id": int(share_pdb_frag + share_pdb_pep),
               "peptides_over_identity_threshold": int(id_pep),
               "fragments_over_identity_threshold": int(id_frag),
               "max_identity_peptides": float(id_p.max()),
               "max_identity_fragments": float(id_f.max()),
               "pool500_best_rmsd": float(rec["rr"][:K_CLAIM].min()),
               "pool500_n_under_0.5A": int((rec["rr"][:K_CLAIM] < 0.5).sum()),
               "fold_filter_matches": bool(mine == theirs),
               "fold_filter_only_mine": len(mine - theirs),
               "fold_filter_only_theirs": len(theirs - mine),
               "pool500_from_peptides": int(sum(not ("_" in str(s)
                                                     and str(s).rsplit("_", 1)[1].isdigit())
                                                for s in src500)),
               "pool500_unique_sources": int(len(set(map(str, src500))))}
        out["per_target"].append(row)
        print(f"{t.pdb:6} fold {fold} lib {len(peps)}+{len(frags)} "
              f"self {self_in_pep + self_in_frag} sharePDB {row['library_shares_target_pdb_id']} "
              f"idHits {id_pep + id_frag} foldfilter_ok {row['fold_filter_matches']}",
              flush=True)
    out["totals"] = {k: int(sum(r[k] for r in out["per_target"]))
                     for k in ("self_seq_in_peptides", "self_seq_in_fragments",
                               "library_shares_target_pdb_id",
                               "peptides_over_identity_threshold",
                               "fragments_over_identity_threshold",
                               "fold_filter_only_mine", "fold_filter_only_theirs")}
    out["totals"]["fold_filter_matches_all"] = all(r["fold_filter_matches"]
                                                   for r in out["per_target"])
    write_json("audit_leak.json", out)
    return out


# --------------------------------------------------------------------------- io
def write_json(name, obj):
    path = os.path.join(HERE, name)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    os.replace(tmp, path)
    print(f"wrote {path}", flush=True)


def read_json(name):
    with open(os.path.join(HERE, name)) as f:
        return json.load(f)


CLAIM = {"rho_band1.5": 0.641, "rho_global": 0.830, "sel": 1.836, "pool": 1.534,
         "prior_rho_band1.5": 0.013, "prior_sel": 3.398}


def stage_report():
    a = read_json("audit_analyse.json")
    s = a["summary"]
    print("\n=== 1. REPRODUCTION, K=500, oracle native distance matrix ===")
    print(f"{'quantity':22} {'claimed':>9} {'reproduced':>11} {'delta':>8}")
    for key, ck in (("K500_rho_global", "rho_global"), ("K500_rho_band1.5", "rho_band1.5"),
                    ("K500_sel", "sel"), ("K500_pool", "pool")):
        print(f"{ck:22} {CLAIM[ck]:9.3f} {s[key]:11.3f} {s[key]-CLAIM[ck]:+8.3f}")
    print("\n=== 2. Is the oracle just RMSD in disguise? ===")
    print(f"  Spearman(distance-L1 score, CA-RMSD) over the pool : {s['K500_rho_global']:+.3f}")
    print(f"  Pearson  (same)                                    : "
          f"{s['K500_pearson_score_rmsd']:+.3f}")
    print(f"  Spearman(dRMSD, CA-RMSD)                           : "
          f"{s['K500_rho_drmsd_rmsd']:+.3f}")
    print(f"  selected by dRMSD instead of L1                    : {s['K500_sel_drmsd']:.3f}")
    print(f"  a PERFECT ranker (= pool best)                     : {s['K500_pool']:.3f}")
    print(f"  random pick from the pool                          : "
          f"{s['K500_sel_random']:.3f}")
    print(f"  pool-mean distance matrix (no sequence info)       : "
          f"{s['K500_sel_poolmean']:.3f}  rho {s['K500_rho_poolmean']:+.3f}")
    w = a["wrong_native_control"]
    print(f"  another target's native matrix (same length)       : {w['sel']:.3f}  "
          f"rho {w['rho_global']:+.3f}  in-band {w['rho_band1.5']:+.3f}")
    print("\n=== 3. Band and pool-size stability ===")
    hdr = f"{'K':>6} " + " ".join(f"{'band'+str(b):>12}" for b in BANDS)
    print(hdr)
    for K in POOL_SIZES:
        print(f"{K:6d} " + " ".join(
            f"{s[f'K{K}_rho_band{b}']:+12.3f}" for b in BANDS))
    print(f"{'sel':>6} " + " ".join(f"{s[f'K500_sel_band{b}']:12.3f}" for b in BANDS)
          + "   (K=500, best-in-band pick)")
    for K in POOL_SIZES:
        print(f"  K={K:5d}  pool {s[f'K{K}_pool']:.3f}  selected {s[f'K{K}_sel']:.3f}  "
              f"global rho {s[f'K{K}_rho_global']:+.3f}")
    bd = a.get("band_diagnostics", {})
    if bd:
        print("\n=== 3b. Does the in-band statistic behave? (K=500) ===")
        print(f"{'band':>6} {'scored':>7} {'dropped':>8} {'size min/med/max':>20} "
              f"{'rho':>7} {'wtd':>7} {'neg':>4} {'rho~size':>9}")
        for b in BANDS:
            d = bd[f"band{b}"]
            print(f"{b:6.1f} {d['n_targets_scored']:7d} {d['n_targets_dropped']:8d} "
                  f"{d['band_size_min']:6.0f}/{d['band_size_median']:.0f}/"
                  f"{d['band_size_max']:.0f}".ljust(48)
                  + f"{d['rho_unweighted']:+7.3f} {d['rho_size_weighted']:+7.3f} "
                  f"{d['n_targets_negative']:4d} {d['corr_bandsize_rho']:+9.3f}")
    print("\n=== 4. Splits ===")
    for k, v in a["splits"].items():
        print(f"  {k:12} n={v['n_targets']:2d} pool {v['pool']:.3f} sel {v['sel']:.3f} "
              f"gap {v['gap']:+.3f} in-band rho "
              f"{v['rho_band1.5'] if v['rho_band1.5'] is not None else float('nan'):+.3f}")
    try:
        nz = read_json("audit_noise.json")
        print("\n=== 5. THE NOISE SWEEP (K=500, 20 seeds/target) ===")
        print(f"{'sigma':>6} {'MAE':>6} {'selected':>9} {'sd':>6} {'rho_glob':>9} "
              f"{'rho_band1.5':>12}")
        for r in nz["iid"]:
            print(f"{r['sigma']:6.2f} {r['mae']:6.3f} {r['sel']:9.3f} {r['sel_sd']:6.3f} "
                  f"{r['rho_global']:+9.3f} {r['rho_band1.5']:+12.3f}")
        print(f"\n{'tau':>6} {'MAE':>6} {'selfRMSD':>9} {'selected':>9} {'rho_band1.5':>12}"
              "   (coordinate-space noise)")
        for r in nz["coord"]:
            print(f"{r['tau']:6.2f} {r['mae']:6.3f} {r['self_rmsd']:9.3f} {r['sel']:9.3f} "
                  f"{r['rho_band1.5']:+12.3f}")
    except FileNotFoundError:
        print("\n(noise stage not run)")
    try:
        pr = read_json("audit_prior.json")
        print("\n=== 5b. A BIASED error model: shrink the native matrix toward the "
              "pool mean ===")
        print(f"{'a':>5} {'MAE':>6} {'selected':>9} {'rho_glob':>9} {'rho_band1.5':>12}")
        for r in pr["shrink"]:
            print(f"{r['a']:5.2f} {r['mae']:6.3f} {r['sel']:9.3f} "
                  f"{r['rho_global']:+9.3f} {r['rho_band1.5']:+12.3f}")
        ps = pr.get("prior_summary")
        if ps:
            print("\n  the SHIPPED distogram (numbers from s5/inband.json, its own pools)")
            print(f"    its MAE                                : {ps['prior_mae']:.3f} A")
            print(f"    what it actually selects               : {ps['prior_sel']:.3f} A")
            print(f"    iid noise at the SAME MAE would select : "
                  f"{ps['iid_sel_at_prior_mae']:.3f} A")
            print(f"    shrinkage at the SAME MAE would select : "
                  f"{ps['shrink_sel_at_prior_mae']:.3f} A")
            print(f"    the oracle on the same targets         : "
                  f"{ps['oracle_sel']:.3f} A")
            print(f"    targets where the prior is WORSE than iid noise of equal MAE: "
                  f"{ps['n_prior_worse_than_iid']}/{ps['n']}")
    except FileNotFoundError:
        pass
    for name, title in (("audit_verify.json", "6. RMSD VERIFICATION / TIE BREAK"),
                        ("audit_leak.json", "7. LEAKAGE")):
        try:
            d = read_json(name)
        except FileNotFoundError:
            continue
        print(f"\n=== {title} ===")
        if "kabsch_vs_protein_geometry" in d:
            k = d["kabsch_vs_protein_geometry"]
            print(f"  my Kabsch vs protein_geometry.ca_rmsd on {k['n']} real windows: "
                  f"max |diff| {k['max_abs_diff']:.2e}")
            m = d["mirror_control"]
            print(f"  mirror: self {m['self_rmsd_mine']:.3e}  reflected "
                  f"{m['mirror_rmsd_mine']:.3f} (geo {m['mirror_rmsd_geo']:.3f})  "
                  f"distance-matrix change {m['mirror_distance_matrix_max_diff']:.2e}")
            c = d["chirality_blindness"]
            print(f"  reflected natives: mean CA-RMSD {c['mean_mirror_rmsd']:.3f} A, "
                  f"score better than every pool member on "
                  f"{c['n_mirror_beats_pool']}/24 targets")
            t = d["tie_break"]
            print(f"  tie-break: mean overlap {t['mean_overlap']:.0f}/500 candidates; "
                  f"pool {t['pool_stable']:.3f}->{t['pool_shuffled']:.3f}  "
                  f"sel {t['sel_stable']:.3f}->{t['sel_shuffled']:.3f}  "
                  f"in-band {t['band_stable']:+.3f}->{t['band_shuffled']:+.3f}")
        else:
            print("  " + json.dumps(d["totals"]))
            print("  " + json.dumps(d["fragment_db"]))


STAGES = {"pools": stage_pools, "verify": stage_verify, "leak": stage_leak,
          "analyse": stage_analyse, "noise": stage_noise, "prior": stage_prior,
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
