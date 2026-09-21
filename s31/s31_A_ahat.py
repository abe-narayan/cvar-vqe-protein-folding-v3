"""S31 lane A -- the a_hat arms through the convex readout (coordinator's request).

The identity ||sum w W - t||^2 = <w,a> - 1/2 w'Bw makes ANY per-candidate quality estimate
a_hat convert into an endpoint with NO tuning: solve the convex program
    w = argmin_{w in simplex}  <w, a_hat>  -  1/2 w'Bw          (gamma = 1, the DERIVED value)
a_hat must therefore be in ANGSTROM^2, which is what the leave-fold-out regression supplies.

ARMS
  1  a_hat = DIS z-rank                     -- the baseline (reproduces lane A's CAL)
  2  a_hat = medoid/consensus criterion      -- lane B's `P.mean(1)` (s31_B3_graph.py:141)
  3  a_hat = DIS + consensus, fitted LFO     -- THE ONE THAT MATTERS
  4  a_hat = const                           -- = MEB, pure dispersion maximisation

Every fit is LEAVE-FOLD-OUT on the pinned 5 folds: fold f's coefficients come from the other
four.  A within-fold fit would be leakage, and the identity makes leakage convert EFFICIENTLY,
which is the worst case.  rho(a_hat, a) is reported BEFORE the Angstroms so the mechanism is
visible first.

SCOPE NOTE, declared: lane B measured `medoid_crit` on the SHIPPED top-75 band.  This runs on
the top-128 quantum hypothesis set, which is the set the readout acts on.  Different set, same
statistic.

Basis: CA POINT CLOUD.  `a` is ORACLE and is used ONLY as a training label on the four training
folds and as the evaluation label; no held-out native enters any a_hat.
"""
import json
import os
import sys
import glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s8 import consensus2 as cc                                        # noqa: E402
from s31.s31_A_r1 import kabsch_rmsd_batch, medoid_uniform             # noqa: E402
from s31.s31_A_readout import min_quad_simplex                         # noqa: E402

UNIV = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
OUT = os.path.join(ROOT, "s31", "results", "s31_A_ahat.json")
ROWS = os.path.join(ROOT, "s31", "results", "s31_A_ahat_rows.jsonl")
DIM = 128
NSHUF = 8


def zs(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / max(v.std(), 1e-12)


def load(f):
    z = np.load(f)
    pdb, n, fold = str(z["pdb"]), int(z["n"]), int(z["fold"])
    pool = z["order"].astype(np.int64)[:500]
    nat = z["nat_ca"].astype(np.float64)
    W500 = z["W"][pool].astype(np.float64)
    DIS = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))["DIS"].astype(np.float64)
    order = np.argsort(DIS, kind="stable")
    top, sub75 = order[:DIM], order[:75]
    Wt = W500[top]

    Pt = np.empty((DIM, DIM), float)
    for a in range(DIM):
        Pt[a] = kabsch_rmsd_batch(Wt, Wt[a])
    Pt = Pt.astype(np.float32).astype(float)
    b_u = medoid_uniform(Pt)

    Sup = cc.superpose_batch(Wt, Wt[b_u])
    t = cc.superpose_batch(nat[None], Wt[b_u])[0]
    Xf = Sup.reshape(DIM, -1)
    tf = t.reshape(-1)
    a = ((Xf - tf[None, :]) ** 2).sum(1)                    # ORACLE label, Angstrom^2
    sq = (Xf ** 2).sum(1)
    B = sq[:, None] + sq[None, :] - 2.0 * (Xf @ Xf.T)
    np.fill_diagonal(B, 0.0)
    B = np.maximum(B, 0.0)
    s_t = float(B.sum() / (DIM * (DIM - 1)))                # NATIVE-FREE per-target scale

    from scipy.stats import rankdata
    rk = rankdata(DIS[top])
    f_dis = zs(rk)                                          # the deployed _zrank
    f_con = zs(Pt.mean(1))                                  # lane B's medoid criterion

    # production baseline on this target
    W75 = W500[sub75]
    P75 = np.empty((75, 75), float)
    for aa in range(75):
        P75[aa] = kabsch_rmsd_batch(W75, W75[aa])
    P75 = P75.astype(np.float32).astype(float)
    b75 = medoid_uniform(P75)
    r_prod = float(kabsch_rmsd_batch(
        cc.superpose_batch(W75, W75[b75]).mean(0)[None], nat)[0])

    return dict(pdb=pdb, n=n, fold=fold, Xf=Xf, nat=nat, a=a, B=B, s_t=s_t,
                F=np.column_stack([f_dis, f_con]), r_prod=r_prod)


def main(limit=None):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    files = UNIV if limit is None else UNIV[:limit]
    D = []
    for i, f in enumerate(files):
        D.append(load(f))
        if (i + 1) % 20 == 0:
            print(f"  loaded {i+1}/{len(files)}", flush=True)

    folds = np.array([d["fold"] for d in D], int)
    gl = sorted(set(folds.tolist()))

    # ---- LEAVE-FOLD-OUT regression of (a - mean a)/s_t on the standardised features.
    #      Dividing by the NATIVE-FREE s_t makes the slope dimensionless and the prediction
    #      target-adaptive; multiplying back gives a_hat in Angstrom^2.
    def fit(cols, tr_mask):
        X, y = [], []
        for d, m in zip(D, tr_mask):
            if not m:
                continue
            X.append(d["F"][:, cols])
            y.append((d["a"] - d["a"].mean()) / d["s_t"])
        X = np.vstack(X)
        y = np.concatenate(y)
        A = np.column_stack([np.ones(len(y)), X])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        return beta

    ARMS = {"1_DIS": [0], "2_CONS": [1], "3_DIS+CONS": [0, 1]}
    ahat = {k: [None] * len(D) for k in ARMS}
    betas = {k: {} for k in ARMS}
    for gf in gl:
        tr = folds != gf
        for k, cols in ARMS.items():
            beta = fit(cols, tr)
            betas[k][int(gf)] = [float(b) for b in beta]
            for i in np.nonzero(folds == gf)[0]:
                d = D[i]
                ahat[k][i] = d["s_t"] * (beta[0] + d["F"][:, cols] @ beta[1:])

    rng = np.random.default_rng(31_000_7)
    rows = []
    with open(ROWS, "w") as fh:
        for i, d in enumerate(D):
            def rms(w):
                C = (np.asarray(w, float) @ d["Xf"]).reshape(d["n"], 3)
                return float(kabsch_rmsd_batch(C[None], d["nat"])[0])

            r = {"pdb": d["pdb"], "fold": d["fold"], "r_prod": d["r_prod"]}
            for k in ARMS:
                ah = ahat[k][i]
                r[f"rho_{k}"] = float(np.corrcoef(ah, d["a"])[0, 1])
                r[f"rmsd_{k}"] = rms(min_quad_simplex(-d["B"], ah))
            r["rho_4_CONST"] = 0.0
            r["rmsd_4_CONST"] = rms(min_quad_simplex(-d["B"], np.zeros(DIM)))
            # shuffled-B control on arm 3, and a shuffled-feature control
            sh = []
            for _ in range(NSHUF):
                q = rng.permutation(DIM)
                sh.append(rms(min_quad_simplex(-d["B"][np.ix_(q, q)], ahat["3_DIS+CONS"][i])))
            r["shufB_3"] = sh
            q = rng.permutation(DIM)
            r["rmsd_3_shuffeat"] = rms(min_quad_simplex(-d["B"], ahat["3_DIS+CONS"][i][q]))
            rows.append(r)
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            if (i + 1) % 20 == 0:
                print(f"  solved {i+1}/{len(D)}", flush=True)

    g = lambda k: np.array([r[k] for r in rows], float)                 # noqa: E731

    def fold_se(x):
        m = np.array([x[folds == f].mean() for f in gl])
        return float(m.std(ddof=1) / np.sqrt(len(m)))

    def cmp(x, tag):
        mu, se = float(x.mean()), fold_se(x)
        mde = 2.8016 * se
        sign = [int(np.sign(x[folds == f].mean())) for f in gl]
        return {"tag": tag, "mean": mu, "SE": se, "MDE": mde,
                "ratio_to_MDE": float(abs(mu) / mde) if mde > 0 else 0.0,
                "W": int((x < 0).sum()), "L": int((x > 0).sum()),
                "folds_same_sign": int(max(sign.count(1), sign.count(-1))),
                "verdict": ("NOT A RESULT" if abs(mu) < 0.7 * mde else
                            "NOT MEASURED" if abs(mu) < mde else
                            ("BETTER" if mu < 0 else "WORSE"))}

    S3 = np.array([r["shufB_3"] for r in rows], float)
    summ = {
        "n": len(rows),
        "PROD75_cloud": float(g("r_prod").mean()),
        "CROSSING_rho_required": 0.211,
        "rho": {k: float(g(f"rho_{k}").mean()) for k in
                ("1_DIS", "2_CONS", "3_DIS+CONS", "4_CONST")},
        "rmsd": {k: float(g(f"rmsd_{k}").mean()) for k in
                 ("1_DIS", "2_CONS", "3_DIS+CONS", "4_CONST")},
        "LFO_betas": betas,
        "shufB_3_per_draw_mean": S3.mean(0).round(4).tolist(),
        "comparisons": [cmp(g(f"rmsd_{k}") - g("r_prod"), f"{k} - PROD75")
                        for k in ("1_DIS", "2_CONS", "3_DIS+CONS", "4_CONST")]
        + [cmp(g("rmsd_3_DIS+CONS") - g("rmsd_1_DIS"), "3_DIS+CONS - 1_DIS"),
           cmp(g("rmsd_3_DIS+CONS") - S3.mean(1), "3 - shuffledB(mean of 8)"),
           cmp(g("rmsd_3_DIS+CONS") - g("rmsd_3_shuffeat"), "3 - shuffled-feature")],
    }
    with open(OUT, "w") as fh:
        json.dump(summ, fh, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
