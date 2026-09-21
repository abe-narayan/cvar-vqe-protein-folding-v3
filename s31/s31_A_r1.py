"""S31 lane A -- R1 audit + A1 (monotone-reweighting theorem) on the 126 dev targets.

Registered in s31/PREREG_S31_A.md (A1-v, A1-e, A1-d) and in the coordinator's R1 request.

Basis: CA POINT CLOUD throughout.  Nothing here is on the built-chain basis.
ORACLE quantities are named *_oracle and are NOT DEPLOYABLE.

Data path reproduces the deployed filter without re-running retrieval:
    pool  = s8/generate_univ/<pdb>.npz  -> order[:500], W
    score = s27/cache/<pdb>.npz         -> DIS   (the shipped distogram Bayes-risk score)
    top   = np.argsort(DIS, kind='stable')[:128]     == core/pipeline.py:757-759
    Pt    = pairwise Kabsch CA-RMSD of the top-128, float32 round-tripped (pipeline _q)
This reproduces S30 lane T's PREFIX arm to 0.0024 A of the canonical 3.0483; the residual is
the reference's float32 pool round trip, which is not reproduced here.  Declared, not hidden.
"""
import json
import os
import sys
import glob
import time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from core import quantum as qm                                    # noqa: E402
from s8 import consensus2 as cc                                   # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s31", "results", "s31_A_r1.json")
ROWS = os.path.join(ROOT, "s31", "results", "s31_A_r1_rows.jsonl")
VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}
DIM = 128
NQ, LAYERS, ITERS, SEED = 7, 3, 50, 0


# ------------------------------------------------------------------ geometry helpers
def kabsch_rmsd_batch(A, b):
    """CA-RMSD of each A[k] (n,3) onto b (n,3).  Kabsch, reflections forbidden."""
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    n = b.shape[0]
    bc = b - b.mean(0)
    Ac = A - A.mean(1, keepdims=True)
    C = np.einsum("kni,nj->kij", Ac, bc)
    U, S, Vt = np.linalg.svd(C)
    det = np.sign(np.linalg.det(np.einsum("kij,kjl->kil", U, Vt)))
    Ssum = S[:, 0] + S[:, 1] + det * S[:, 2]
    sq = (Ac ** 2).sum((1, 2)) + (bc ** 2).sum()
    return np.sqrt(np.maximum(sq - 2 * Ssum, 0.0) / n)


def medoid_uniform(P):
    return int(np.argmin(P.mean(1)))


def medoid_weighted(P, w):
    return int(np.argmin(P @ (np.asarray(w, float) / max(float(np.sum(w)), 1e-12))))


# ------------------------------------------------------------------ A1 closed form
def p_star(E, alpha, T):
    """THEOREM A1.  The unconstrained argmin of F(p) = CVaR_alpha(E;p) - T H(p), T > 0.

    By Rockafellar-Uryasev, CVaR_alpha(E;p) = max_mu [ mu - (1/alpha) sum_x p_x (mu-E_x)_+ ],
    a max of functions LINEAR in p, so CVaR is CONVEX in p and F is strictly convex.  At the
    optimum the envelope theorem gives dCVaR/dp_x = -(mu-E_x)_+/alpha, so simplex
    stationarity reads

        p*_x  propto  exp( (mu - E_x)_+ / (alpha T) )

    -- a BOLTZMANN law at temperature alpha*T below the threshold mu and FLAT above it --
    with mu the unique scalar fixing sum_{x: E_x < mu} p*_x = alpha.  The optimum generically
    sits at the KINK where that mass equals alpha exactly, which is why mu is a free
    continuous parameter rather than one of the E values; a scan over the atoms of E cannot
    represent it and lands ~1.5e-3 high in p (measured).  At alpha = 1 the threshold is above
    every E_x and the law collapses to the plain Boltzmann exp(-E_x/T).

    p* is a NON-INCREASING function of E_x alone.  Returns (p, mu).
    """
    E = np.asarray(E, float)
    aT = alpha * T
    if alpha >= 1.0 - 1e-12:
        z = -(E - E.min()) / T
        p = np.exp(z)
        return p / p.sum(), float(np.inf)

    def pm(mu):
        z = np.maximum(mu - E, 0.0) / aT
        z = z - z.max()
        p = np.exp(z)
        return p / p.sum()

    def G(mu):
        p = pm(mu)
        return float(p[E < mu].sum())

    lo, hi = float(E.min()) - 1e-9, float(E.max()) + 1e-9
    if G(hi) < alpha:                      # alpha unreachable: the whole set is the tail
        return pm(hi), hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if G(mid) < alpha:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    return pm(mu), float(mu)


def F_of_p(E, p, alpha, T):
    v, _q, _dp = qm.cvar_exact(E, p, alpha)
    lp = np.log(np.maximum(p, 1e-300))
    return float(v + T * float((p * lp).sum()))


def solve_F_numeric(E, alpha, T, iters=6000, lr=0.05, seed=0):
    """Independent check on p_star: mirror descent on the simplex, no closed form used."""
    rng = np.random.default_rng(seed)
    D = len(E)
    g = rng.normal(0, 0.1, D)
    for _ in range(iters):
        p = np.exp(g - g.max())
        p /= p.sum()
        _v, q, dp = qm.cvar_exact(E, p, alpha)
        grad = dp + T * (np.log(np.maximum(p, 1e-300)) + 1.0)
        g = g - lr * (grad - float(grad @ p))
    p = np.exp(g - g.max())
    return p / p.sum()


def kkt_residual(E, p, alpha, T):
    """The DECISIVE A1 verification: simplex stationarity of F at p, solver-independent.

    grad F_x = dCVaR/dp_x + T(log p_x + 1) must be CONSTANT in x at an interior optimum.
    Returns max_x |grad F_x - mean(grad F)|, scaled by |grad F| so it is a relative figure.
    This replaces `maxabs(closed - numeric)` as the verification of record: that statistic
    measures the REFERENCE SOLVER's convergence, not the closed form's correctness, and
    A1-v was registered against it by mistake.
    """
    _v, _q, dp = qm.cvar_exact(E, p, alpha)
    g = dp + T * (np.log(np.maximum(p, 1e-300)) + 1.0)
    return float(np.abs(g - g.mean()).max() / max(np.abs(g).max(), 1e-12))


def kl_bits(p, q):
    p = np.asarray(p, float) + 1e-300
    q = np.asarray(q, float) + 1e-300
    return float((p * (np.log2(p) - np.log2(q))).sum())


def entropy_bits(p):
    p = np.asarray(p, float)
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


# ------------------------------------------------------------------ per target
def one(f):
    z = np.load(f)
    pdb, n, fold = str(z["pdb"]), int(z["n"]), int(z["fold"])
    pool = z["order"].astype(np.int64)[:500]
    nat = z["nat_ca"].astype(np.float64)
    W500 = z["W"][pool].astype(np.float64)
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
    order = np.argsort(DIS, kind="stable")             # == pipeline filter_pool
    top = order[:DIM]
    sub75 = order[:75]
    Wt = W500[top]

    # ---- Pt: pairwise Kabsch CA-RMSD, float32 round trip (pipeline._q)
    Pt = np.empty((DIM, DIM), float)
    for a in range(DIM):
        Pt[a] = kabsch_rmsd_batch(Wt, Wt[a])
    Pt = Pt.astype(np.float32).astype(float)

    # ---- the deployed quantum state
    alpha, T = VQE_LFO[fold % 5]
    E = qm._zrank(DIS[top]) if hasattr(qm, "_zrank") else None
    if E is None:
        from scipy.stats import rankdata
        r = rankdata(DIS[top])
        E = (r - r.mean()) / max(r.std(), 1e-12)
    t0 = time.time()
    p_th, cvar, Hnat, _circ = qm.run_cvar_vqe(E, alpha, T, n=NQ, layers=LAYERS,
                                              iters=ITERS, seed=SEED)
    t_vqe = time.time() - t0
    p_th = np.asarray(p_th, float)
    p_th = p_th / p_th.sum()

    # ---- A1: the closed form and its agreement
    ps, qq = p_star(E, alpha, T)
    p_num = solve_F_numeric(E, alpha, T)
    a1v = float(np.abs(ps - p_num).max())
    a1_kkt_star = kkt_residual(E, ps, alpha, T)
    a1_kkt_theta = kkt_residual(E, p_th, alpha, T)
    a1_F_closed, a1_F_num = F_of_p(E, ps, alpha, T), F_of_p(E, p_num, alpha, T)
    # monotonicity of p* in E (the theorem's structural claim)
    oE = np.argsort(E, kind="stable")
    mono = bool(np.all(np.diff(ps[oE]) <= 1e-12))

    # ---- R1 part 1: the SELECTOR
    b_w = medoid_weighted(Pt, p_th)
    b_u = medoid_uniform(Pt)
    b_s = medoid_weighted(Pt, ps)
    # R1 point 4: every vertex reachable
    reach = int(sum(1 for j in range(DIM) if medoid_weighted(Pt, np.eye(DIM)[j]) == j))

    # ---- R1 part 2: the SYNTHESIS (pipeline.py:880-895), the arm that emits structure
    Sup_w = cc.superpose_batch(Wt, Wt[b_w])
    C_th = np.tensordot(p_th, Sup_w, axes=(0, 0))
    C_st = np.tensordot(ps, Sup_w, axes=(0, 0))
    C_un = Sup_w.mean(0)
    # is the emitted structure a pool member?  distance to the NEAREST one, in Angstrom
    d_member = float(kabsch_rmsd_batch(Wt, C_th).min())

    rr_top = kabsch_rmsd_batch(Wt, nat)
    r_sel_w = float(rr_top[b_w])
    r_sel_u = float(rr_top[b_u])
    r_synth_th = float(kabsch_rmsd_batch(C_th[None], nat)[0])
    r_synth_st = float(kabsch_rmsd_batch(C_st[None], nat)[0])
    r_synth_un = float(kabsch_rmsd_batch(C_un[None], nat)[0])

    # production: uniform mean of the DIS top-75, medoid frame of the 75
    W75 = W500[sub75]
    P75 = np.empty((75, 75), float)
    for a in range(75):
        P75[a] = kabsch_rmsd_batch(W75, W75[a])
    P75 = P75.astype(np.float32).astype(float)
    b75 = medoid_uniform(P75)
    C75 = cc.superpose_batch(W75, W75[b75]).mean(0)
    r_prod75 = float(kabsch_rmsd_batch(C75[None], nat)[0])

    return {
        "pdb": pdb, "n": n, "fold": fold, "alpha": alpha, "T": T, "t_vqe": t_vqe,
        "H_theta_bits": float(Hnat) / np.log(2.0), "H_star_bits": entropy_bits(ps),
        "kl_theta_star_bits": kl_bits(p_th, ps), "kl_star_theta_bits": kl_bits(ps, p_th),
        "tv_theta_star": float(0.5 * np.abs(p_th - ps).sum()),
        "a1v_maxabs": a1v, "a1_kkt_star": a1_kkt_star, "a1_kkt_theta": a1_kkt_theta, "a1_F_closed": a1_F_closed, "a1_F_num": a1_F_num,
        "a1_mono": mono, "q_star": qq,
        "sel_w": b_w, "sel_u": b_u, "sel_star": b_s, "sel_agree_wu": int(b_w == b_u),
        "sel_agree_ws": int(b_w == b_s), "reach_vertices": reach,
        "d_nearest_member": d_member,
        "rmsd_sel_w": r_sel_w, "rmsd_sel_u": r_sel_u,
        "rmsd_synth_theta": r_synth_th, "rmsd_synth_star": r_synth_st,
        "rmsd_synth_uniform": r_synth_un, "rmsd_prod75": r_prod75,
        "rmsd_top128_best_oracle": float(rr_top.min()),
        "rmsd_top128_mean": float(rr_top.mean()),
    }


def main(limit=None):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    files = UNIV if limit is None else UNIV[:limit]
    rows = []
    with open(ROWS, "w") as fh:
        for i, f in enumerate(files):
            r = one(f)
            rows.append(r)
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{len(files)}  {r['pdb']}  kl={r['kl_theta_star_bits']:.4f}",
                      flush=True)
    A = {k: np.array([r[k] for r in rows], float) for k in rows[0]
         if isinstance(rows[0][k], (int, float, bool))}
    folds = np.array([r["fold"] for r in rows], int)

    def fold_se(d):
        m = np.array([d[folds == g].mean() for g in sorted(set(folds.tolist()))])
        return float(m.std(ddof=1) / np.sqrt(len(m)))

    summ = {
        "n": len(rows),
        "A1_verification_maxabs_closed_vs_numeric": float(A["a1v_maxabs"].max()),
        "A1_KKT_residual_closedform_max": float(A["a1_kkt_star"].max()),
        "A1_KKT_residual_p_theta_max": float(A["a1_kkt_theta"].max()),
        "A1_KKT_residual_p_theta_median": float(np.median(A["a1_kkt_theta"])),
        "A1_F_closed_minus_numeric_max": float((A["a1_F_closed"] - A["a1_F_num"]).max()),
        "A1_monotone_all": bool(A["a1_mono"].all()),
        "A1_KL_theta_star_bits_mean": float(A["kl_theta_star_bits"].mean()),
        "A1_KL_theta_star_bits_median": float(np.median(A["kl_theta_star_bits"])),
        "A1_KL_theta_star_bits_p90": float(np.percentile(A["kl_theta_star_bits"], 90)),
        "A1_TV_theta_star_mean": float(A["tv_theta_star"].mean()),
        "A1_dRMSD_theta_minus_star_mean": float((A["rmsd_synth_theta"]
                                                 - A["rmsd_synth_star"]).mean()),
        "A1_dRMSD_theta_minus_star_absmean": float(np.abs(A["rmsd_synth_theta"]
                                                          - A["rmsd_synth_star"]).mean()),
        "A1_dRMSD_theta_minus_star_SE": fold_se(A["rmsd_synth_theta"]
                                                - A["rmsd_synth_star"]),
        "R1_reach_vertices_min": int(A["reach_vertices"].min()),
        "R1_reach_vertices_mean": float(A["reach_vertices"].mean()),
        "R1_sel_agree_weighted_uniform": float(A["sel_agree_wu"].mean()),
        "R1_sel_agree_weighted_star": float(A["sel_agree_ws"].mean()),
        "R1_d_nearest_member_mean": float(A["d_nearest_member"].mean()),
        "R1_d_nearest_member_min": float(A["d_nearest_member"].min()),
        "rmsd_sel_w_mean": float(A["rmsd_sel_w"].mean()),
        "rmsd_sel_u_mean": float(A["rmsd_sel_u"].mean()),
        "rmsd_synth_theta_mean": float(A["rmsd_synth_theta"].mean()),
        "rmsd_synth_star_mean": float(A["rmsd_synth_star"].mean()),
        "rmsd_synth_uniform_mean": float(A["rmsd_synth_uniform"].mean()),
        "rmsd_prod75_mean": float(A["rmsd_prod75"].mean()),
        "rmsd_top128_best_ORACLE_mean": float(A["rmsd_top128_best_oracle"].mean()),
        "H_theta_bits_mean": float(A["H_theta_bits"].mean()),
        "H_star_bits_mean": float(A["H_star_bits"].mean()),
        "t_vqe_total_s": float(A["t_vqe"].sum()),
    }
    with open(OUT, "w") as fh:
        json.dump(summ, fh, indent=2)
    for k, v in summ.items():
        print(f"{k:52s} {v}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
