#!/usr/bin/env python
"""s32/s32_P_sign.py -- S32 LANE P, H-P3 follow-up: the common mode is ONE BIT away, per target.

PREREG s32/PREREG_S32_P.md @ 33dfe0d3 (H-P3, P3-6 extension; declared EXPLORATORY in
s32/MULTIPLICITY.md because the sign-rule family was not named in the prereg).

WHY THIS EXISTS.  `s32_P_inband.py` measured cos(mu_hat, mu) for seven native-free direction
estimates and found the MEAN cos indistinguishable from zero for all of them -- while mean cos^2
ran 0.066 to 0.228 against a same-space random control.  A direction with real |cos| and no
resolvable sign is exactly S14's `in-band-ordering-is-per-target` finding in new coordinates:
**the magnitude is there, the per-target SIGN is not.**  This file asks whether any native-free
rule supplies that one bit, and prices what the bit is worth.

    mu = c - t   (the common mode, in the medoid frame of the score's top-128)
    pc1          = the leading right singular vector of the pool deviations d, sign canonicalised
    the bit      = sign(<pc1, mu>)                       ORACLE, 1 bit per target
    the payoff   = |mu| * (1 - sqrt(1 - cos^2))          the common mode a signed pc1 removes

Every sign rule fixes its polarity LEAVE-FOLD-OUT on the pinned folds, so choosing the
convention cannot read the held-out native.  ORACLE / NOT DEPLOYABLE for the bit itself and for
every cos; the RULES are native-free.

    python s32/s32_P_sign.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I                                        # noqa: E402
from scipy.stats import rankdata                                       # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s32_P_sign.json")
ROWS = os.path.join(RESULTS, "s32_P_sign_rows.jsonl")
DIM, SEED, NRAND, NPC = 128, 32_0_9311, 200, 10
EPS = 0.25          # Angstrom step for the two-sided directional probes


def zs(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / max(v.std(), 1e-12)


def rg(C):
    return float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean()))


def risk_of(dg, C):
    grid = np.asarray(dg["grid"], float); risk = np.asarray(dg["risk"], float)
    ii = np.asarray(dg["i"], int); jj = np.asarray(dg["j"], int)
    D = np.linalg.norm(C[ii] - C[jj], axis=1)
    g = np.clip(((D - grid[0]) / 0.05).astype(np.int64), 0, len(grid) - 1)
    return float(risk[np.arange(risk.shape[0]), g].mean())


def one(t, rng, placebo_nats=()):
    pdb = t["pdb"]; n = int(t["n"])
    u = I.load_univ(pdb)
    dg = I.distogram(pdb, t["seq"], t["fold"])
    z27 = np.load(os.path.join(ROOT, "s27", "cache", f"{pdb}.npz"))
    DIS = z27["DIS"].astype(float)
    LEG = z27["LEG"].astype(float); AMB = z27["AMB"].astype(float)
    pool = np.asarray(u["order"], int)[:500]
    W = np.asarray(u["W"], float)[pool]
    nat = np.asarray(u["nat_ca"], float)
    order = np.argsort(DIS, kind="stable"); top = order[:DIM]
    Wt = W[top]
    Pt = I.pairwise_rmsd(Wt).astype(np.float32).astype(float)
    b = I.medoid(Pt)
    Sup = I.superpose_batch(Wt, Wt[b]); tt = I.superpose_batch(nat[None], Wt[b])[0]
    X = Sup.reshape(DIM, -1); tf = tt.reshape(-1)
    c = X.mean(0); mu = c - tf; d = X - c[None, :]
    V = (d ** 2).sum(1)
    nmu = float(np.linalg.norm(mu))

    Us, Sv, Vt = np.linalg.svd(d, full_matrices=False)
    r = int((Sv > 1e-9 * Sv[0]).sum())
    PC = Vt[:NPC]
    #: canonical polarity: the entry of largest magnitude is made positive.  Deterministic,
    #: native-free, and independent of LAPACK's arbitrary sign.
    for k in range(len(PC)):
        PC[k] = PC[k] * np.sign(PC[k][int(np.argmax(np.abs(PC[k])))])
    cos_pc = np.array([float(PC[k] @ mu / max(nmu, 1e-12)) for k in range(len(PC))])
    cum = float(np.cumsum(cos_pc ** 2)[-1])

    #: matched random control (rule 6): a direction drawn uniformly in span(d), the space the
    #: estimators live in -- measured, not the analytic 1/r.
    Q = Vt[:r].T
    cr = []
    for _ in range(NRAND):
        g = Q @ rng.standard_normal(r); g /= max(np.linalg.norm(g), 1e-12)
        cr.append(float(g @ mu / max(nmu, 1e-12)))
    cr = np.array(cr)

    p1 = PC[0]
    proj = d @ p1                                              # (128,) pool coords along pc1
    bit = float(np.sign(p1 @ mu))                              # ORACLE, the thing to predict
    acos = abs(float(p1 @ mu / max(nmu, 1e-12)))

    # ------------------------------------------------------------------ native-free sign rules
    Cn = c.reshape(n, 3); P1 = p1.reshape(n, 3)
    rgp = float(np.sqrt(np.mean(np.asarray(dg["expected"], float) ** 2) / 2.0))

    def two_sided(f):
        return float(f(Cn + EPS * P1) - f(Cn - EPS * P1))

    z = {}
    z["DIS_corr"] = float(np.corrcoef(proj, zs(rankdata(DIS[top])))[0, 1])
    z["LEG_corr"] = float(np.corrcoef(proj, zs(rankdata(LEG[top])))[0, 1])
    z["AMB_corr"] = float(np.corrcoef(proj, zs(rankdata(AMB[top])))[0, 1])
    z["V_corr"] = float(np.corrcoef(proj, V)[0, 1])
    z["SKEW"] = float(((proj - proj.mean()) ** 3).mean() / max(proj.std() ** 3, 1e-12))
    z["RISK_probe"] = two_sided(lambda C: risk_of(dg, C))
    z["RG_probe"] = two_sided(lambda C: abs(rg(C) - rgp))
    w = np.exp(-zs(DIS[top]) / 2.0); w /= w.sum()
    z["SOFT_disp"] = float(p1 @ (c - (w[:, None] * X).sum(0)))
    #: `top` holds indices into the K=500 pool, and the pool IS BLOSUM-ordered, so `top` itself
    #: is the BLOSUM rank.  A first version used the position WITHIN `top`, which is the DIS
    #: rank by construction -- it reproduced DIS_corr to four decimals and was a null input,
    #: not a second question.  Kept as a comment because the defect is the point.
    z["BLOSUM_corr"] = float(np.corrcoef(proj, zs(rankdata(top)))[0, 1])

    # ---------------------------------------------------------------- PLACEBO CONTROL (rule 5)
    #: `pool-error-is-68-percent-common-mode` records that a MISMATCHED same-length native serves
    #: the rescaling correction as well as the true one.  If cos^2(mu, pc1) is a property of the
    #: pool's geometry against ANY plausible compact chain, the alignment below is not about this
    #: target's native and the finding is void.  This control can fail, which is why it is here.
    plac = []
    for q in placebo_nats:
        tp = I.superpose_batch(q[None], Wt[b])[0].reshape(-1)
        mp = c - tp
        plac.append(float(p1 @ mp / max(np.linalg.norm(mp), 1e-12)))
    plac = np.array(plac) if len(plac) else np.zeros(1)

    return dict(pdb=pdb, n=n, fold=int(t["fold"]), fail18=bool(pdb in I.FAIL18),
                rank_d=r, mu_rmsd=float(nmu / np.sqrt(n)), meanV=float(V.mean()),
                bit=bit, abs_cos_pc1=acos, cos_pc=[float(v) for v in cos_pc],
                cum_cos2_10pc=cum, rand_cos2_mean=float((cr ** 2).mean()),
                rand_cos2_sd=float((cr ** 2).std(ddof=1)),
                placebo_cos2_mean=float((plac ** 2).mean()), n_placebo=int(len(plac)), z=z)


def main():
    ts = I.targets()
    if os.path.exists(ROWS):
        R = [json.loads(l) for l in open(ROWS) if l.strip()]
        if len(R) == len(ts):
            print("reusing %d cached rows" % len(R), flush=True)
            return analyse(R)
    by_n = {}
    for t in ts:
        by_n.setdefault(int(t["n"]), []).append(t["pdb"])
    R = []
    for k, t in enumerate(ts):
        rng = np.random.default_rng(SEED + abs(hash(t["pdb"])) % 100003)
        others = [p for p in by_n[int(t["n"])] if p != t["pdb"]][:4]
        pl = [np.asarray(I.load_univ(p)["nat_ca"], float) for p in others]
        R.append(one(t, rng, pl))
        if (k + 1) % 25 == 0:
            print("  %d/126" % (k + 1), flush=True)
    with open(ROWS, "w") as fh:
        for r in R:
            fh.write(json.dumps(r) + "\n")
    return analyse(R)


def analyse(R):
    folds = np.array([r["fold"] for r in R]); gl = sorted(set(folds.tolist()))
    bit = np.array([r["bit"] for r in R])
    ac = np.array([r["abs_cos_pc1"] for r in R])
    mur = np.array([r["mu_rmsd"] for r in R])
    keys = list(R[0]["z"].keys())

    def fold_se(x):
        m = np.array([np.mean(x[folds == f]) for f in gl])
        return float(m.std(ddof=1) / np.sqrt(len(m))), [float(v) for v in m]

    #: THE ACCOUNTING.  Correcting the centroid by a step `alpha` along a unit direction `p`
    #: leaves an error |mu - alpha*p|, i.e. per-target RMSD sqrt(mur^2 - 2*alpha_r*mur*cos +
    #: alpha_r^2) with alpha_r = alpha/sqrt(n).  A WRONG SIGN MAKES IT WORSE; an earlier version
    #: zeroed the loss on wrong-sign targets, which flattered every rule and made a random sign
    #: look worth +0.227 A.  `alpha` must be native-free, so it is a SINGLE global step fitted
    #: LEAVE-FOLD-OUT; only the ORACLE arm is allowed a per-target alpha.
    AGRID = np.linspace(0.0, 3.0, 121)

    def emit(cs, alpha_r):
        return np.sqrt(np.maximum(mur ** 2 - 2.0 * alpha_r * mur * cs + alpha_r ** 2, 0.0))

    def payoff(signs, per_target_alpha=False):
        cs = signs * bit * ac                          # achieved cos, signed
        if per_target_alpha:                           # ORACLE alpha: the textbook optimum
            e = mur * np.sqrt(np.maximum(1.0 - np.maximum(cs, 0.0) ** 2, 0.0))
            e = np.where(cs > 0, e, mur)
            return float(mur.mean() - e.mean()), float(cs.mean()), None
        al = np.empty(len(R))
        for f in gl:                                   # one global step, fitted out of fold
            tr = folds != f
            j = int(np.argmin([emit(cs, a)[tr].mean() for a in AGRID]))
            al[folds == f] = AGRID[j]
        e = np.sqrt(np.maximum(mur ** 2 - 2.0 * al * mur * cs + al ** 2, 0.0))
        return float(mur.mean() - e.mean()), float(cs.mean()), float(al.mean())

    rules = {}
    for k in keys:
        zz = np.array([r["z"][k] for r in R])
        pred = np.sign(zz); pred[pred == 0] = 1.0
        #: polarity fixed LEAVE-FOLD-OUT: the sign convention for fold f is the one that agrees
        #: with the ORACLE bit on the OTHER four folds.  This is the only place the native
        #: enters a deployable-looking rule and it is held out.
        sig = np.empty(len(R))
        for f in gl:
            tr = folds != f
            sig[folds == f] = np.sign(np.mean(pred[tr] * bit[tr])) or 1.0
        agree = (pred * sig == bit).astype(float)
        se, per = fold_se(agree)
        pay, mcos, al = payoff(pred * sig)
        rules[k] = dict(accuracy=float(agree.mean()), se=se, per_fold=per,
                        mde=2.8016 * se, over_mde=float(abs(agree.mean() - 0.5) / (2.8016 * se)),
                        polarity_by_fold={int(f): float(sig[folds == f][0]) for f in gl},
                        achieved_mean_cos=mcos, payoff_A=pay, alpha_LFO=al,
                        verdict=("NOT A RESULT" if abs(agree.mean() - 0.5) < 0.7 * 2.8016 * se else
                                 "NOT MEASURED" if abs(agree.mean() - 0.5) < 2.8016 * se else
                                 "SIGNAL"))

    rng = np.random.default_rng(SEED)
    rnd = [payoff(rng.choice([-1.0, 1.0], len(R)))[0] for _ in range(NRAND)]
    #: DEFECT FIXED: an earlier version passed `np.ones` here, which leaves cs = +-|cos|
    #: -- a random sign wearing an ORACLE label.  The ORACLE arm must pass `bit`, so that
    #: cs = bit*bit*|cos| = +|cos|.  Caught because ORACLE_mean_cos came back 0.0588 when
    #: mean |cos(pc1, mu)| is 0.3943; two numbers that must agree and did not.
    orc, orc_cos, orc_al = payoff(bit)                                  # ORACLE sign, global step
    orc_pt, _, _ = payoff(bit, per_target_alpha=True)                   # ORACLE sign AND step
    assert abs(orc_cos - float(ac.mean())) < 1e-12, (orc_cos, float(ac.mean()))

    cosmat = np.array([r["cos_pc"] for r in R])
    out = dict(
        prereg="s32/PREREG_S32_P.md @ 33dfe0d3 (H-P3 extension, EXPLORATORY family)",
        basis="CA POINT CLOUD, score-top-128 medoid frame; mu = c - t; NOT an endpoint",
        ORACLE="ORACLE / NOT DEPLOYABLE -- `bit`, every cos and every payoff read the native",
        n=126, dim=DIM, eps=EPS, n_rand=NRAND,
        geometry=dict(rank_d_median=float(np.median([r["rank_d"] for r in R])),
                      rank_d_min=int(min(r["rank_d"] for r in R)),
                      rank_d_max=int(max(r["rank_d"] for r in R)),
                      mu_rmsd_mean=float(mur.mean()),
                      mu_rmsd_note="|mu|/sqrt(n): the RMSD the common mode alone would give"),
        pc_spectrum=dict(
            cos2_by_pc=[float(v) for v in (cosmat ** 2).mean(0)],
            cum_cos2_by_pc=[float(v) for v in np.cumsum((cosmat ** 2).mean(0))],
            matched_random_cos2=float(np.mean([r["rand_cos2_mean"] for r in R])),
            matched_random_cos2_sd=float(np.mean([r["rand_cos2_sd"] for r in R])),
            analytic_null_1_over_rank=float(np.mean([1.0 / r["rank_d"] for r in R])),
            pc1_over_random=float((cosmat[:, 0] ** 2).mean() /
                                  max(np.mean([r["rand_cos2_mean"] for r in R]), 1e-12)),
            abs_cos_pc1_mean=float(ac.mean()), abs_cos_pc1_median=float(np.median(ac))),
        the_bit=dict(
            what="sign(<pc1, mu>) -- ONE BIT PER TARGET; payoffs are A of |mu|-RMSD removed, "
                 "baseline mu_rmsd_mean, ORACLE / NOT DEPLOYABLE",
            baseline_mu_rmsd=float(mur.mean()),
            ORACLE_sign_global_step_payoff_A=orc, ORACLE_mean_cos=orc_cos, ORACLE_step=orc_al,
            ORACLE_sign_and_per_target_step_payoff_A=orc_pt,
            random_sign_payoff_A_mean=float(np.mean(rnd)),
            random_sign_payoff_A_sd=float(np.std(rnd, ddof=1)),
            balance=float((bit > 0).mean())),
        placebo_control=dict(
            what="cos^2(pc1, c - t') for up to 4 MISMATCHED same-length natives t'",
            placebo_cos2=float(np.mean([r["placebo_cos2_mean"] for r in R])),
            true_cos2=float((cosmat[:, 0] ** 2).mean()),
            matched_random_cos2=float(np.mean([r["rand_cos2_mean"] for r in R])),
            n_placebo_mean=float(np.mean([r["n_placebo"] for r in R]))),
        rules=rules,
        best_rule=max(rules, key=lambda k: rules[k]["over_mde"]),
        multiplicity_emitted=len(rules) + 1)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
