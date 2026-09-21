"""S31 lane A -- (a) R1's realised capacity, and (b) the ansatz's inductive bias.

Both on the CA POINT-CLOUD basis.  ORACLE quantities are labelled.

(a) REALISED CAPACITY of the SELECTION readout `argmin_i (P p)_i`.
    * the alphabet is the number of DISTINCT structures, not 2^k: duplicate candidates give
      P[i,j] = 0 off the diagonal and np.argmin returns the FIRST, so those vertices are
      unreachable.  `core/pipeline.py:766-773` already dedups when BUILDING P; the readout
      does not.
    * value delivered, in S30's currency: bits = k - log2(rank of the selected candidate in
      the ORACLE RMSD order).  ORACLE -- a diagnostic, not a deployable quantity.

(b) INDUCTIVE BIAS.  The candidate index x IS the DIS rank (`o = top[:dim]`), so qubit b is a
    BIT of the rank.  A product state gives p(x) = prod_b q_b^{x_b}(1-q_b)^{1-x_b}; with equal
    per-qubit odds that is a function of POPCOUNT(x) alone, not of the rank's VALUE.  The
    objective's optimum p* is a monotone function of E_x, and E_x is monotone in the rank's
    VALUE.  So the hypothesis is: the ansatz and the objective live in different bases.

    Measured: R^2 of log p on popcount(x) vs on E_x, for p_theta, for p*, and over a SURVEY of
    4000 random draws from the ansatz's own initialisation law -- which characterises the
    reachable MANIFOLD, not one optimised point.

    Intervention (native-free, uses only E): relabel candidates so rank i goes to the i-th
    bitstring in (popcount, value) order, putting the objective into the ansatz's basis.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from core import quantum as qm                                      # noqa: E402
from s8 import consensus2 as cc                                     # noqa: E402
from s31.s31_A_r1 import (kabsch_rmsd_batch, medoid_uniform, medoid_weighted,  # noqa: E402
                          p_star, kl_bits, entropy_bits)

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s31", "results", "s31_A_cap.json")
ROWS = os.path.join(ROOT, "s31", "results", "s31_A_cap_rows.jsonl")
VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}
DIM, NQ, LAYERS, ITERS = 128, 7, 3, 50
NSEED = 8
POP = np.array([bin(x).count("1") for x in range(DIM)], float)
HW_ORDER = np.array(sorted(range(DIM), key=lambda x: (bin(x).count("1"), x)), int)


def r2(y, X):
    """R^2 of an OLS of y on [1, X...]."""
    X = np.atleast_2d(np.asarray(X, float))
    if X.shape[0] != len(y):
        X = X.T
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ beta
    tot = float(((y - y.mean()) ** 2).sum())
    return float(1.0 - (res ** 2).sum() / max(tot, 1e-300))


def bond_dims(psi, n):
    """Schmidt rank of the real amplitude vector across every contiguous cut."""
    out = []
    for c in range(1, n):
        M = np.asarray(psi, float).reshape(1 << c, 1 << (n - c))
        s = np.linalg.svd(M, compute_uv=False)
        out.append(int((s > 1e-10 * max(s[0], 1e-300)).sum()))
    return out


def one(f, survey_rng):
    z = np.load(f)
    pdb, n, fold = str(z["pdb"]), int(z["n"]), int(z["fold"])
    pool = z["order"].astype(np.int64)[:500]
    nat = z["nat_ca"].astype(np.float64)
    W500 = z["W"][pool].astype(np.float64)
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
    top = np.argsort(DIS, kind="stable")[:DIM]
    Wt = W500[top]

    Pt = np.empty((DIM, DIM), float)
    for a in range(DIM):
        Pt[a] = kabsch_rmsd_batch(Wt, Wt[a])
    Pt = Pt.astype(np.float32).astype(float)
    rr = kabsch_rmsd_batch(Wt, nat)                       # ORACLE
    orank = np.empty(DIM, int)
    orank[np.argsort(rr, kind="stable")] = np.arange(1, DIM + 1)

    # ---- (a) the alphabet
    n_distinct = len({Wt[a].tobytes() for a in range(DIM)})
    reach = int(sum(1 for j in range(DIM) if int(np.argmin(Pt[:, j])) == j))

    from scipy.stats import rankdata
    rk = rankdata(DIS[top])
    E = (rk - rk.mean()) / max(rk.std(), 1e-12)
    alpha, T = VQE_LFO[fold % 5]

    p_th, _c, Hnat, circ = qm.run_cvar_vqe(E, alpha, T, n=NQ, layers=LAYERS,
                                           iters=ITERS, seed=0)
    p_th = np.asarray(p_th, float); p_th /= p_th.sum()
    ps, mu_star = p_star(E, alpha, T)
    # the EXACT shape of log p*: a hinge in E.  R^2 on this basis is 1 by construction for p*,
    # so it is the right yardstick for "how much of the circuit's log-p is the right shape".
    hinge = np.maximum((mu_star if np.isfinite(mu_star) else E.max() + 1.0) - E, 0.0)

    b_w, b_u, b_s = medoid_weighted(Pt, p_th), medoid_uniform(Pt), medoid_weighted(Pt, ps)
    b_e = int(np.argmin(E))

    # ---- (b) basis mismatch, on the optimised point and on p*
    lp_th, lp_st = np.log(np.maximum(p_th, 1e-300)), np.log(np.maximum(ps, 1e-300))
    bias = {
        "r2_logptheta_popcount": r2(lp_th, POP),
        "r2_logptheta_E": r2(lp_th, E),
        "r2_logptheta_hinge": r2(lp_th, hinge),
        "r2_logptheta_hinge_pop": r2(lp_th, np.column_stack([hinge, POP])),
        "r2_logpstar_popcount": r2(lp_st, POP),
        "r2_logpstar_E": r2(lp_st, E),
        "r2_logpstar_hinge": r2(lp_st, hinge),
    }
    # the reachable MANIFOLD, from the ansatz's own initialisation law N(0, 0.6)
    TH = survey_rng.normal(0.0, 0.6, (400, circ.n_params()))
    PR = circ.probs_batch(TH)
    LP = np.log(np.maximum(PR, 1e-300))
    bias["r2_survey_popcount_mean"] = float(np.mean([r2(LP[i], POP) for i in range(400)]))
    bias["r2_survey_E_mean"] = float(np.mean([r2(LP[i], E) for i in range(400)]))
    bias["r2_survey_hinge_mean"] = float(np.mean([r2(LP[i], hinge) for i in range(400)]))
    # the CEILING of the manifold on the right shape: best R^2 over the 400 draws
    bias["r2_survey_hinge_max"] = float(np.max([r2(LP[i], hinge) for i in range(400)]))
    bias["survey_entropy_bits_mean"] = float(np.mean([entropy_bits(PR[i]) for i in range(400)]))
    bd = bond_dims(circ.state(TH[0]), NQ)
    bias["bond_dims_random_theta"] = bd

    # ---- (b) intervention: relabel so rank i -> the i-th (popcount, value) bitstring
    Eperm = np.empty(DIM)
    Eperm[HW_ORDER] = np.sort(E)                       # best candidates at low Hamming weight
    inv = np.empty(DIM, int)
    inv[HW_ORDER] = np.argsort(E, kind="stable")       # bitstring -> original candidate
    p_hw, _c2, _H2, _c3 = qm.run_cvar_vqe(Eperm, alpha, T, n=NQ, layers=LAYERS,
                                          iters=ITERS, seed=0)
    p_hw = np.asarray(p_hw, float); p_hw /= p_hw.sum()
    p_hw_back = np.zeros(DIM)
    p_hw_back[inv] = p_hw                              # back to candidate order
    ps_hw = np.zeros(DIM)
    ps_hw[inv] = p_star(Eperm, alpha, T)[0]
    kl_hw = kl_bits(p_hw_back, ps_hw)

    def synth(w):
        b = medoid_weighted(Pt, w)
        S = cc.superpose_batch(Wt, Wt[b])
        return float(kabsch_rmsd_batch((np.asarray(w, float) / w.sum() @ S.reshape(DIM, -1)
                                        ).reshape(n, 3)[None], nat)[0])

    # ---- A1-d: seed sensitivity at the deployed T and at T = 0
    seeds = {}
    for tag, TT in (("T_dep", T), ("T_zero", 0.0)):
        vals, ents = [], []
        for sd in range(NSEED):
            pp, _a, _b, _cc2 = qm.run_cvar_vqe(E, alpha, TT, n=NQ, layers=LAYERS,
                                               iters=ITERS, seed=sd)
            pp = np.asarray(pp, float); pp /= pp.sum()
            vals.append(synth(pp))
            ents.append(entropy_bits(pp))
        seeds[tag] = {"sd_rmsd": float(np.std(vals, ddof=1)),
                      "mean_rmsd": float(np.mean(vals)),
                      "mean_entropy_bits": float(np.mean(ents))}

    return {"pdb": pdb, "n": n, "fold": fold, "alpha": alpha, "T": T,
            "n_distinct": n_distinct, "reach": reach,
            "orank_sel_w": int(orank[b_w]), "orank_sel_u": int(orank[b_u]),
            "orank_sel_star": int(orank[b_s]), "orank_sel_E": int(orank[b_e]),
            "H_theta_bits": float(Hnat) / np.log(2.0), "H_star_bits": entropy_bits(ps),
            "kl_theta_star": kl_bits(p_th, ps), "kl_hw_star": kl_hw,
            "H_hw_bits": entropy_bits(p_hw_back),
            "r_synth_theta": synth(p_th), "r_synth_star": synth(ps),
            "r_synth_hw": synth(p_hw_back), "r_synth_unif": synth(np.ones(DIM)),
            "seeds": seeds, **bias}


def main(limit=None, nseeds_from=0):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    files = UNIV if limit is None else UNIV[:limit]
    rows = []
    with open(ROWS, "w") as fh:
        for i, f in enumerate(files):
            r = one(f, np.random.default_rng(1234 + i))
            rows.append(r)
            fh.write(json.dumps(r) + "\n"); fh.flush()
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{len(files)} {r['pdb']}", flush=True)

    folds = np.array([r["fold"] for r in rows], int)
    gl = sorted(set(folds.tolist()))
    g = lambda k: np.array([r[k] for r in rows], float)              # noqa: E731

    def fold_se(d):
        m = np.array([d[folds == gg].mean() for gg in gl])
        return float(m.std(ddof=1) / np.sqrt(len(m)))

    def cmp(d, tag):
        mu, se = float(d.mean()), fold_se(d)
        mde = 2.8016 * se
        return {"tag": tag, "mean": mu, "SE": se, "MDE": mde,
                "ratio_to_MDE": float(abs(mu) / mde) if mde > 0 else 0.0,
                "verdict": ("NOT A RESULT" if abs(mu) < 0.7 * mde else
                            "NOT MEASURED" if abs(mu) < mde else
                            ("BETTER" if mu < 0 else "WORSE"))}

    bits = lambda k: float(np.log2(DIM) - np.log2(g(k)).mean())      # noqa: E731
    summ = {
        "n": len(rows),
        "ALPHABET_n_distinct_mean": float(g("n_distinct").mean()),
        "ALPHABET_n_distinct_min": float(g("n_distinct").min()),
        "ALPHABET_reachable_vertices_mean": float(g("reach").mean()),
        "ALPHABET_reachable_vertices_min": float(g("reach").min()),
        "CAPACITY_bits_alphabet_mean": float(np.log2(g("reach")).mean()),
        "ORACLE_bits_delivered_sel_theta": bits("orank_sel_w"),
        "ORACLE_bits_delivered_sel_uniform": bits("orank_sel_u"),
        "ORACLE_bits_delivered_sel_pstar": bits("orank_sel_star"),
        "ORACLE_bits_delivered_argmin_E": bits("orank_sel_E"),
        "ORACLE_mean_rank_sel_theta": float(g("orank_sel_w").mean()),
        "ORACLE_mean_rank_sel_uniform": float(g("orank_sel_u").mean()),
        "BIAS_r2_logptheta_popcount": float(g("r2_logptheta_popcount").mean()),
        "BIAS_r2_logptheta_E": float(g("r2_logptheta_E").mean()),
        "BIAS_r2_logpstar_popcount": float(g("r2_logpstar_popcount").mean()),
        "BIAS_r2_logpstar_E": float(g("r2_logpstar_E").mean()),
        "BIAS_r2_survey_popcount": float(g("r2_survey_popcount_mean").mean()),
        "BIAS_r2_survey_E": float(g("r2_survey_E_mean").mean()),
        "BIAS_r2_logptheta_hinge": float(g("r2_logptheta_hinge").mean()),
        "BIAS_r2_logptheta_hinge_pop": float(g("r2_logptheta_hinge_pop").mean()),
        "BIAS_r2_logpstar_hinge": float(g("r2_logpstar_hinge").mean()),
        "BIAS_r2_survey_hinge_mean": float(g("r2_survey_hinge_mean").mean()),
        "BIAS_r2_survey_hinge_max": float(g("r2_survey_hinge_max").mean()),
        "BIAS_survey_entropy_bits": float(g("survey_entropy_bits_mean").mean()),
        "BIAS_bond_dims": rows[0]["bond_dims_random_theta"],
        "HW_kl_to_star_mean": float(g("kl_hw_star").mean()),
        "HW_kl_baseline_mean": float(g("kl_theta_star").mean()),
        "HW_entropy_bits_mean": float(g("H_hw_bits").mean()),
        "means": {k: float(g(k).mean()) for k in
                  ("r_synth_theta", "r_synth_star", "r_synth_hw", "r_synth_unif",
                   "H_theta_bits", "H_star_bits")},
        "A1d_seed_sd_T_deployed": float(np.mean([r["seeds"]["T_dep"]["sd_rmsd"] for r in rows])),
        "A1d_seed_sd_T_zero": float(np.mean([r["seeds"]["T_zero"]["sd_rmsd"] for r in rows])),
        "A1d_entropy_T_deployed": float(np.mean([r["seeds"]["T_dep"]["mean_entropy_bits"]
                                                 for r in rows])),
        "A1d_entropy_T_zero": float(np.mean([r["seeds"]["T_zero"]["mean_entropy_bits"]
                                             for r in rows])),
        "comparisons": [
            cmp(g("r_synth_hw") - g("r_synth_theta"), "HW-relabel - deployed theta"),
            cmp(g("r_synth_theta") - g("r_synth_star"), "theta - p* (both synth)"),
            cmp(g("r_synth_theta") - g("r_synth_unif"), "theta - uniform128"),
        ],
    }
    with open(OUT, "w") as fh:
        json.dump(summ, fh, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
