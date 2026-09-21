"""S31 lane A -- A3: the EXACT readout identity, the Hamiltonian it forces, and its arms.

Registered in s31/PREREG_S31_A.md sections 2-4.  Basis: CA POINT CLOUD throughout.

THE IDENTITY (exact, no approximation, needs only sum_x w_x = 1 -- NOT w >= 0):

    || sum_x w_x W_x - t ||_F^2   =   <w, a>  -  1/2 * w' B w

    a_x = ||W_x - t||_F^2      ORACLE, per-candidate squared error
    B_xy = ||W_x - W_y||_F^2   NATIVE-FREE, and 1/2 w'Bw = tr Sigma_w, the weighted dispersion

Consequences this script measures:
  * the native-free half of the exact readout objective is PAIRWISE and carries a MINUS sign,
    so at fixed quality the readout should MAXIMISE weighted mutual spread;
  * the forced mean-field Hamiltonian is H[w] = diag(a_hat) - B, whose VMC local energy
    E_loc(x) = a_hat_x - (Bw)_x is the exact gradient of the readout error;
  * on the simplex, d/dw of -gamma/2 w'Bw is CONVEX for gamma > 0 and CONCAVE for gamma < 0,
    so the ATTRACTIVE (consensus) sign cannot produce a distribution at all -- its minimiser
    is a vertex.  That is a theorem, checked here, not a measurement.

All ORACLE arms are labelled ORACLE and are NOT DEPLOYABLE.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s8 import consensus2 as cc                                     # noqa: E402
from s31.s31_A_r1 import kabsch_rmsd_batch, medoid_uniform, medoid_weighted, p_star  # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s31", "results", "s31_A_readout.json")
ROWS = os.path.join(ROOT, "s31", "results", "s31_A_readout_rows.jsonl")
DIM = 128
GAMMA_GRID = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 1e9]
NSHUF = 8


# ------------------------------------------------------------------ solvers
def proj_simplex(v):
    u = np.sort(v)[::-1]
    c = np.cumsum(u) - 1.0
    r = np.arange(1, len(v) + 1)
    rho = np.nonzero(u - c / r > 0)[0][-1]
    return np.maximum(v - c[rho] / (rho + 1.0), 0.0)


def min_quad_simplex(Q, c, iters=4000):
    """argmin_{w in simplex} <c,w> + 1/2 w'Qw, accelerated projected gradient.

    Q must be PSD on the simplex tangent space, which is what makes the problem convex.
    """
    D = len(c)
    L = float(np.abs(np.linalg.eigvalsh(Q)).max()) + 1e-9
    w = np.full(D, 1.0 / D)
    y, tk = w.copy(), 1.0
    for _ in range(iters):
        g = c + Q @ y
        wn = proj_simplex(y - g / L)
        tn = 0.5 * (1 + np.sqrt(1 + 4 * tk * tk))
        y = wn + ((tk - 1) / tn) * (wn - w)
        w, tk = wn, tn
    return w


def min_affine(Xf, tf):
    """argmin_{1'w=1} || w'Xf - tf ||, by least squares in the tangent coordinates.

    w = w0 + Z v with w0 uniform and Z a basis of {1'z = 0}, so the constraint is exact and
    the problem is an unconstrained lstsq.  This is READOUT 2 (s27/s28_A_amp.py:105-117):
    an AFFINE combination, negative weights allowed.  ORACLE -- NOT DEPLOYABLE.
    """
    D = Xf.shape[0]
    w0 = np.full(D, 1.0 / D)
    Z = np.eye(D)[:, 1:] - np.eye(D)[:, :1]          # columns sum to 0, rank D-1
    A = (Z.T @ Xf)                                   # (D-1, 3n)
    r = tf - w0 @ Xf
    v, *_ = np.linalg.lstsq(A.T, r, rcond=None)
    return w0 + Z @ v


# ------------------------------------------------------------------ per target
def one(f, rng):
    z = np.load(f)
    pdb, n, fold = str(z["pdb"]), int(z["n"]), int(z["fold"])
    pool = z["order"].astype(np.int64)[:500]
    nat = z["nat_ca"].astype(np.float64)
    W500 = z["W"][pool].astype(np.float64)
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
    order = np.argsort(DIS, kind="stable")
    top, sub75 = order[:DIM], order[:75]
    Wt = W500[top]

    # pairwise Kabsch RMSD block (the DEPLOYED object) -- used only for the frame + A3-frame
    Pt = np.empty((DIM, DIM), float)
    for a in range(DIM):
        Pt[a] = kabsch_rmsd_batch(Wt, Wt[a])
    Pt = Pt.astype(np.float32).astype(float)
    b_u = medoid_uniform(Pt)

    # ---- the EXOGENOUS frame: everything superposed onto the uniform medoid
    Sup = cc.superpose_batch(Wt, Wt[b_u])                       # (128, n, 3)
    t = cc.superpose_batch(nat[None], Wt[b_u])[0]               # native into the same frame
    Xf = Sup.reshape(DIM, -1)
    tf = t.reshape(-1)
    Ef = Xf - tf[None, :]                                       # ORACLE error matrix
    a = (Ef ** 2).sum(1)                                        # ORACLE
    sq = (Xf ** 2).sum(1)
    B = sq[:, None] + sq[None, :] - 2.0 * (Xf @ Xf.T)           # NATIVE-FREE
    np.fill_diagonal(B, 0.0)
    B = np.maximum(B, 0.0)
    G = Ef @ Ef.T                                               # ORACLE Gram

    # ---- A3-i: the identity, on simplex draws AND affine draws with negative weights
    err_i = []
    for kind in ("simplex", "affine"):
        for _ in range(100):
            if kind == "simplex":
                w = rng.dirichlet(np.ones(DIM) * 0.3)
            else:
                w = rng.normal(size=DIM)
                w = w / w.sum()
            lhs = float(((w @ Xf - tf) ** 2).sum())
            rhs = float(w @ a - 0.5 * w @ B @ w)
            err_i.append(abs(lhs - rhs) / max(abs(lhs), 1e-12))
    a3i = float(max(err_i))

    # ---- A3-frame: does the DEPLOYED pairwise-RMSD block satisfy the identity?
    Bdep = (Pt ** 2) * n                                        # n * RMSD^2 = ||.||_F^2 per pair
    frame_err = []
    for _ in range(50):
        w = rng.dirichlet(np.ones(DIM) * 0.3)
        lhs = float(((w @ Xf - tf) ** 2).sum())
        rhs = float(w @ a - 0.5 * w @ Bdep @ w)
        frame_err.append(abs(lhs - rhs) / max(abs(lhs), 1e-12))
    a3frame = float(np.median(frame_err))

    def rms(w):
        C = (np.asarray(w, float) @ Xf).reshape(n, 3)
        return float(kabsch_rmsd_batch(C[None], nat)[0])

    # ---- baselines
    W75 = W500[sub75]
    P75 = np.empty((75, 75), float)
    for aa in range(75):
        P75[aa] = kabsch_rmsd_batch(W75, W75[aa])
    P75 = P75.astype(np.float32).astype(float)
    b75 = medoid_uniform(P75)
    r_prod75 = float(kabsch_rmsd_batch(
        cc.superpose_batch(W75, W75[b75]).mean(0)[None], nat)[0])
    w_unif = np.full(DIM, 1.0 / DIM)
    r_top128u = rms(w_unif)
    rr_top = kabsch_rmsd_batch(Wt, nat)

    # ---- ORACLE ceilings.  NOT DEPLOYABLE.
    w_qp = min_quad_simplex(2.0 * G, np.zeros(DIM))
    r_qp_oracle = rms(w_qp)
    w_af = min_affine(Xf, tf)
    rank_aff = int(np.linalg.matrix_rank(Xf - Xf.mean(0), tol=1e-8))
    r_affine_oracle = rms(w_af)

    # ---- MEB: argmax_simplex 1/2 w'Bw  ==  minimum-enclosing-ball dual.  Native-free.
    w_meb = min_quad_simplex(-B, np.zeros(DIM))                 # min -1/2 w'Bw = MAX dispersion
    r_meb = rms(w_meb)

    # ---- the gamma family.  s = native-free per-target scale.
    s = float(B.sum() / (DIM * (DIM - 1)))
    zc = (DIS[top] - DIS[top].mean()) / max(DIS[top].std(), 1e-12)
    from scipy.stats import rankdata
    rk = rankdata(DIS[top])
    zr = (rk - rk.mean()) / max(rk.std(), 1e-12)                 # the deployed _zrank
    gam = {}
    for g in GAMMA_GRID:
        if g >= 1e8:
            gam[str(g)] = r_meb
        elif g == 0.0:
            gam[str(g)] = float(rr_top[int(np.argmin(zr))])      # vertex: the argmin selector
        else:
            w = min_quad_simplex((g / s) * (-B), zr)
            gam[str(g)] = rms(w)
    # the ATTRACTIVE sign: theorem says the minimiser is a vertex -> argmin of z
    r_gamma_negative_vertex = float(rr_top[int(np.argmin(zr))])

    # ---- oracle calibration slope for CAL (fitted LEAVE-FOLD-OUT in main)
    ac = a - a.mean()
    cal_num = float(ac @ zr)
    cal_den = float(zr @ zr)

    # ---- the PRICE of the quality estimate: how good must a_hat be?  ORACLE sweep.
    zs = zr / max(zr.std(), 1e-12)
    as_ = (a - a.mean()) / max(a.std(), 1e-12)
    price = {}
    for lam in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        ah = (1 - lam) * zs + lam * as_
        rho = float(np.corrcoef(ah, a)[0, 1])
        # gamma fixed at the DERIVED value 1 after putting a_hat into Angstrom^2 via a.std()
        w = min_quad_simplex((1.0 / a.std()) * (-B), ah)
        price[f"{lam:.1f}"] = [rho, rms(w)]

    # ---- controls: shuffled B (spectrum preserved, correspondence destroyed)
    shuf = []
    for _ in range(NSHUF):
        q = rng.permutation(DIM)
        Bs = B[np.ix_(q, q)]
        ws = min_quad_simplex((1.0 / s) * (-Bs), zr)
        shuf.append(rms(ws))
    # ---- control: shuffled score
    q = rng.permutation(DIM)
    r_shufscore = rms(min_quad_simplex((1.0 / s) * (-B), zr[q]))

    return {"pdb": pdb, "n": n, "fold": fold,
            "a3i_max_rel": a3i, "a3frame_median_rel": a3frame,
            "r_prod75": r_prod75, "r_top128u": r_top128u,
            "r_top128_best_ORACLE": float(rr_top.min()),
            "r_argmin_score": float(rr_top[int(np.argmin(zr))]),
            "r_qp_simplex_ORACLE": r_qp_oracle, "r_affine_ORACLE": r_affine_oracle,
            "r_meb": r_meb, "gam": gam, "s_scale": s,
            "cal_num": cal_num, "cal_den": cal_den,
            "r_gamma_negative_vertex": r_gamma_negative_vertex,
            "shufB": shuf, "r_shufscore": r_shufscore,
            "w_meb_support": int((w_meb > 1e-6).sum()),
            "w_qp_support": int((w_qp > 1e-6).sum()),
            "w_affine_negmass": float(np.abs(np.minimum(w_af, 0)).sum()),
            "rank_affine": rank_aff, "price": price,
            "w_gam1_support": int((min_quad_simplex((1.0 / s) * (-B), zr) > 1e-6).sum())}


def main(limit=None):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    files = UNIV if limit is None else UNIV[:limit]
    rows = []
    with open(ROWS, "w") as fh:
        for i, f in enumerate(files):
            rng = np.random.default_rng(abs(hash(("s31A", os.path.basename(f)))) % 2**32)
            r = one(f, rng)
            rows.append(r)
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{len(files)} {r['pdb']}", flush=True)

    folds = np.array([r["fold"] for r in rows], int)
    gl = sorted(set(folds.tolist()))

    def fold_se(d):
        m = np.array([d[folds == g].mean() for g in gl])
        return float(m.std(ddof=1) / np.sqrt(len(m)))

    def cmp(d, tag):
        mu, se = float(d.mean()), fold_se(d)
        mde = 2.8016 * se
        sign = [int(np.sign(d[folds == g].mean())) for g in gl]
        return {"tag": tag, "mean": mu, "SE": se, "MDE": mde,
                "ratio_to_MDE": float(abs(mu) / mde) if mde > 0 else 0.0,
                "median": float(np.median(d)),
                "W": int((d < 0).sum()), "L": int((d > 0).sum()),
                "folds_same_sign": int(max(sign.count(1), sign.count(-1))),
                "verdict": ("NOT A RESULT" if abs(mu) < 0.7 * mde else
                            "NOT MEASURED" if abs(mu) < mde else
                            ("BETTER" if mu < 0 else "WORSE"))}

    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    gamgrid = {str(x): np.array([r["gam"][str(x)] for r in rows], float) for x in GAMMA_GRID}

    # ---- LEAVE-FOLD-OUT gamma (GAM) and LFO calibration (CAL)
    gam_lfo = np.empty(len(rows))
    gam_pick = {}
    for gf in gl:
        tr = folds != gf
        best = min(GAMMA_GRID, key=lambda x: gamgrid[str(x)][tr].mean())
        gam_pick[int(gf)] = best
        gam_lfo[folds == gf] = gamgrid[str(best)][folds == gf]
    # CAL: gamma implied by regressing ORACLE a on the deployed zrank, LFO, rescaled by s
    num, den, sc = g("cal_num"), g("cal_den"), g("s_scale")
    cal_lfo_gamma, cal_r = {}, np.empty(len(rows))
    for gf in gl:
        tr = folds != gf
        beta = float(num[tr].sum() / max(den[tr].sum(), 1e-12))   # a ~ beta * zrank, Angstrom^2
        # objective <w, beta*z> - 1/2 w'Bw  ==  <w,z> - (s/beta)/2/s * w'Bw  -> gamma = s/beta
        gcal = float(np.median(sc[tr] / max(beta, 1e-12)))
        cal_lfo_gamma[int(gf)] = gcal
        near = min(GAMMA_GRID, key=lambda x: abs(x - gcal))
        cal_r[folds == gf] = gamgrid[str(near)][folds == gf]

    summ = {
        "n": len(rows),
        "A3i_identity_max_rel_err": float(g("a3i_max_rel").max()),
        "A3frame_deployed_block_median_rel_err": float(np.median(g("a3frame_median_rel"))),
        "means": {k: float(g(k).mean()) for k in
                  ("r_prod75", "r_top128u", "r_argmin_score", "r_top128_best_ORACLE",
                   "r_qp_simplex_ORACLE", "r_affine_ORACLE", "r_meb", "r_shufscore")},
        "gamma_curve_mean": {k: float(v.mean()) for k, v in gamgrid.items()},
        "gamma_LFO_pick": gam_pick, "gamma_LFO_mean": float(gam_lfo.mean()),
        "CAL_LFO_gamma": cal_lfo_gamma, "CAL_mean": float(cal_r.mean()),
        "meb_support_mean": float(g("w_meb_support").mean()),
        "qp_support_mean": float(g("w_qp_support").mean()),
        "affine_negmass_mean": float(g("w_affine_negmass").mean()),
        "rank_affine_mean": float(g("rank_affine").mean()),
        "gam1_support_mean": float(g("w_gam1_support").mean()),
        "price_curve": {k: [float(np.mean([r["price"][k][0] for r in rows])),
                            float(np.mean([r["price"][k][1] for r in rows]))]
                        for k in rows[0]["price"]},
        "shufB_mean_of_draws": float(np.array([r["shufB"] for r in rows], float).mean()),
        "shufB_per_draw_mean": np.array([r["shufB"] for r in rows],
                                        float).mean(0).round(4).tolist(),
        "comparisons": [
            cmp(cal_r - g("r_prod75"), "PRIMARY  CAL - PROD75"),
            cmp(gam_lfo - g("r_prod75"), "GAM(LFO) - PROD75"),
            cmp(g("r_meb") - g("r_prod75"), "MEB - PROD75"),
            cmp(g("r_meb") - g("r_top128u"), "MEB - TOP128U"),
            cmp(gam_lfo - g("r_top128u"), "GAM(LFO) - TOP128U"),
            cmp(g("r_qp_simplex_ORACLE") - g("r_prod75"), "ORACLE QP(simplex) - PROD75"),
            cmp(g("r_affine_ORACLE") - g("r_prod75"), "ORACLE affine - PROD75"),
        ],
    }
    with open(OUT, "w") as fh:
        json.dump(summ, fh, indent=2)
    print(json.dumps(summ, indent=2)[:6000])


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
