#!/usr/bin/env python
"""s29/s29_D_theory_check.py -- lane D: an INDEPENDENT check of lane T's registered predictions.

Every number here is ORACLE where it touches `u` (the direction to the native) or the native's
percentile; nothing is deployable and nothing is tuned. The code is written from lane T's
STATEMENTS in `s29/THEORY.md` sections 2.6 and 3e, not from lane T's code: `s29/s29_T_spectra.py`
is never imported, so a shared implementation error cannot produce agreement.

S29-L7, section 2.6 (three clauses):
  identity  (1) <-grad S(c), u> = -<g, r>, with g the shipped per-pair coefficients
                w_a (2 F_a(d_a(c)) - 1) / P and r = Jc u.  Predicted: relative error < 1e-6 on
                126/126.  My implementation of g is from the risk TABLE's own slope (what the
                shipped objective actually differentiates) and, separately, from the posterior
                CDF (what the theory says it is); both are compared with `Surrogate.value_grad`.
  beta      (2) beta = OLS slope of b = m - tau on a = D(c) - tau, m the posterior's per-pair
                median map, tau the length-matched TYPICAL map.  Predicted: beta > 1 on >= 2/3
                of targets and sign(cos_DIS) = -sign(beta - 1) on >= 70%, exceptions concentrated
                at small |beta - 1|.  Two independent tau definitions are used (the theorem's
                sensitivity to tau is the adversary's question, not lane T's): LIB75, the
                sequence-blind 75-window draw S24 made (rebuilt by lane O's own rule), and UNIV,
                the mean map over the whole leakage-safe universe.
  shrink    (3) re-score with the target map shrunk to tau + s (m - tau), s in {1.0, 0.75, 0.5}
                (and a finer grid), re-read the meter's cosine and native percentile.  Predicted:
                the cosine rises monotonically as s falls, crossing zero near s = 1/median(beta),
                while the native percentile gets WORSE.  This is the measured demonstration that
                contract addendum 20's guard is necessary, so the shrink signature (bond, Rg) is
                reported on the same grid.

S29-L11, section 3e prediction (2):
  variance      for ANY unit-spectral-norm observable at n >= 7, Var[dF/dtheta_0] = r_stable/D^2
                within 2x, and J* = D sqrt(0.0305/r_stable).  Falsified by any observable at
                n >= 7 departing by over 3x.  My observables are NOT lane T's (which were the
                pool's Gaussian graph, its double-centered form and the signed agreement matrix):
                random GOE, random diagonal, rank-k projectors, a Pauli string, a k-regular
                circulant and a rank-one projector, whose stable ranks span 1 to ~D/2.

    python s29/s29_D_theory_check.py identity [--limit N]
    python s29/s29_D_theory_check.py beta     [--limit N] [--tau LIB75|UNIV|both]
    python s29/s29_D_theory_check.py shrink   [--limit N]
    python s29/s29_D_theory_check.py variance [--n 7,8,9] [--draws 120]
    python s29/s29_D_theory_check.py all

Results `s29/results/s29_D_theory_{identity,beta,shrink,variance}.json`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A_amp as A             # noqa: E402
from s27 import s28_A2_local as A2         # noqa: E402
from s29 import s29_D_cost_audit as M      # noqa: E402

RESULTS = os.path.join(HERE, "results")
M_BLIND = 75                      # S24's blind draw size, lane O's rule
SHRINK_GRID = (1.0, 0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3)
REGISTERED_S = (1.0, 0.75, 0.5)   # lane T's registered three


# ============================================================ shared pieces
def target_pieces(pdb):
    """(cand, C0 production cloud, dg, i, j, D0 production pair distances, u ORACLE direction)."""
    cand, dis, top, dg, frame, sur = M.load_target(pdb)
    rungs = M.load_rungs(pdb)
    C0 = rungs["S"]["PROD"] if rungs else None
    if C0 is None:
        from s24 import d_harness as H
        C0, _ = H.readout_uniform(cand, top)
    i, j = I.pair_index(cand.n)
    D0 = I.pair_dists(np.asarray(C0, float)[None], i, j)[0]
    u = A2.oracle_direction(C0, cand.nat_ca)                      # ORACLE
    return cand, np.asarray(C0, float), dg, i, j, D0, u, sur


def risk_slope_at(dg, d):
    """d(risk_a)/dd at distance d_a, from the shipped risk TABLE's own linear interpolation
    (this is exactly what `Surrogate.value_grad` differentiates)."""
    risk = np.asarray(dg["risk"], float); grid = np.asarray(dg["grid"], float)
    g0, h = float(grid[0]), float(grid[1] - grid[0])
    G = risk.shape[1]
    idx = np.clip(np.floor((np.asarray(d, float) - g0) / h).astype(int), 0, G - 2)
    rows = np.arange(risk.shape[0])
    return (risk[rows, idx + 1] - risk[rows, idx]) / h


def cdf_coeff_at(dg, d):
    """The theory's coefficient w_a (2 F_a(d) - 1) with F the posterior CDF on the bin centres.

    `w` is not stored in the distogram cache, but the risk table is w_a * E|d - c_b|, so
    w_a = risk_a(d_far) / E|d_far - c_b| at any grid point; taking the mean over the grid is
    exact up to float32 storage.  Returned as (coeff, w).
    """
    prob = np.asarray(dg["prob"], float)
    cen = np.asarray(dg["centres"], float)
    grid = np.asarray(dg["grid"], float)
    risk = np.asarray(dg["risk"], float)
    base = (prob[:, None, :] * np.abs(grid[None, :, None] - cen[None, None, :])).sum(2)
    w = (risk * base).sum(1) / np.maximum((base ** 2).sum(1), 1e-30)     # least squares, exact up to float32
    F = np.cumsum(prob, axis=1)
    k = np.searchsorted(cen, np.asarray(d, float))                        # centres below d
    rows = np.arange(len(cen := cen)) if False else np.arange(prob.shape[0])
    Fd = np.where(k > 0, F[rows, np.clip(k - 1, 0, prob.shape[1] - 1)], 0.0)
    return w * (2.0 * Fd - 1.0), w


def median_map(dg):
    """The posterior's per-pair median map m: the risk's own argmin on the shipped grid (the
    quantity the objective pulls each distance toward), with ties averaged.  The discrete
    median on the bin centres is returned beside it as a cross-check."""
    risk = np.asarray(dg["risk"], float); grid = np.asarray(dg["grid"], float)
    mn = risk.min(1, keepdims=True)
    tie = np.isclose(risk, mn, rtol=1e-12, atol=1e-12)
    m_grid = (tie * grid[None, :]).sum(1) / tie.sum(1)
    prob = np.asarray(dg["prob"], float); cen = np.asarray(dg["centres"], float)
    F = np.cumsum(prob, 1)
    m_disc = cen[np.argmax(F >= 0.5, axis=1)]
    return m_grid, m_disc


def typical_map(pdb, u_univ, i, j, which="LIB75"):
    """tau: the length-matched TYPICAL distance map, native-free.

    LIB75 -- the mean pair-distance map of S24's sequence-blind 75-window draw, rebuilt by lane
    O's rule (`s29/s29_O_ladder.py :: blind_cloud`, `SD.stable_rng("biasalign", pdb)`).
    UNIV  -- the mean pair-distance map over the WHOLE leakage-safe universe for that target.
    """
    W = np.asarray(u_univ["W"], float)
    if which == "UNIV":
        return I.pair_dists(W, i, j).mean(0), dict(source="mean map over the whole leakage-safe universe", n=len(W))
    if which == "POOL":
        #: a THIRD tau the theorem does not ask for, as a sensitivity: the target's own 500-member
        #: retrieval pool's mean map (not target-blind, so not the theorem's tau; reported to show
        #: which "typical" the sign law would need).
        idx = I.pool_idx(u_univ)
        return I.pair_dists(W[idx], i, j).mean(0), dict(source="the target's own 500-member pool mean map (NOT target-blind)", n=int(len(idx)))
    rng = SD.stable_rng("biasalign", pdb)
    c1 = rng.choice(len(W), M_BLIND, replace=False)
    return I.pair_dists(W[c1], i, j).mean(0), dict(source="s24/biasalign.py C1 draw (75 uniform universe windows)", n=M_BLIND)


# ============================================================ clause 1: the identity
def identity(pdbs):
    rows = []
    for q, pdb in enumerate(pdbs):
        cand, C0, dg, i, j, D0, u, sur = target_pieces(pdb)
        P = len(i)
        # LHS: <-grad S(c), u> with grad S from the shipped surrogate
        val, gradC = sur.value_grad(C0)
        lhs = -float(np.asarray(gradC, float).ravel() @ np.asarray(u, float).ravel())
        # RHS: -<g, r>, g from the risk table's slope / P, r = Jc u
        v = C0[i] - C0[j]
        nrm = np.maximum(np.linalg.norm(v, axis=1), 1e-12)
        r = ((v / nrm[:, None]) * (u[i] - u[j])).sum(1)                   # (P,) = Jc u
        g_tab = risk_slope_at(dg, D0) / P
        rhs_tab = -float(g_tab @ r)
        # the same with the theory's CDF coefficients
        g_cdf, w = cdf_coeff_at(dg, D0)
        rhs_cdf = -float((g_cdf / P) @ r)
        den = max(abs(lhs), 1e-30)
        rows.append(dict(pdb=pdb, n=int(cand.n), P=int(P), fold=int(cand.fold),
                         lhs=lhs, rhs_table=rhs_tab, rhs_cdf=rhs_cdf,
                         rel_table=abs(rhs_tab - lhs) / den, rel_cdf=abs(rhs_cdf - lhs) / den,
                         cos_check=A2.cosine(-A2.remove_rigid(gradC, C0), A2.remove_rigid(u, C0)),
                         kernel_dim=int(P - (3 * cand.n - 6))))
        if (q + 1) % 30 == 0:
            print(f"  [{q+1}/{len(pdbs)}]", flush=True)
    rt = np.array([r["rel_table"] for r in rows]); rc = np.array([r["rel_cdf"] for r in rows])
    out = dict(check="S29-L7 clause 1 (the informativeness identity)", n=len(rows), oracle=True,
               rel_table_max=float(rt.max()), rel_table_median=float(np.median(rt)),
               rel_cdf_max=float(rc.max()), rel_cdf_median=float(np.median(rc)),
               n_table_under_1e6=int((rt < 1e-6).sum()), n_cdf_under_1e6=int((rc < 1e-6).sum()),
               n_cdf_under_1e3=int((rc < 1e-3).sum()),
               kernel_dim_mean=float(np.mean([r["kernel_dim"] for r in rows])),
               kernel_frac_mean=float(np.mean([r["kernel_dim"] / r["P"] for r in rows])),
               verdict=("HOLDS: relative error < 1e-6 on %d/%d targets (table coefficients)"
                        % (int((rt < 1e-6).sum()), len(rows))) if (rt < 1e-6).all() else
                       ("FAILS on %d targets (max rel %.2e)" % (int((rt >= 1e-6).sum()), rt.max())),
               rows=rows)
    print(f"  identity: table coefficients rel err max {rt.max():.3e} median {np.median(rt):.3e} "
          f"({int((rt < 1e-6).sum())}/{len(rows)} under 1e-6)")
    print(f"            CDF coefficients   rel err max {rc.max():.3e} median {np.median(rc):.3e} "
          f"({int((rc < 1e-6).sum())}/{len(rows)} under 1e-6, {int((rc < 1e-3).sum())}/{len(rows)} under 1e-3)")
    print(f"            ker(Jc^T) is {100*out['kernel_frac_mean']:.1f}% of pair space (mean over targets)")
    return out


# ============================================================ clause 2: beta and the sign law
def beta_check(pdbs, taus=("LIB75", "UNIV")):
    cosrows = {r["pdb"]: r for r in (json.loads(l) for l in
               open(os.path.join(ROOT, "s27", "results", "s28_A2_cosine_rows.jsonl"), encoding="utf-8") if l.strip())}
    rows = []
    for q, pdb in enumerate(pdbs):
        cand, C0, dg, i, j, D0, u, sur = target_pieces(pdb)
        uu = I.load_univ(pdb)
        m_grid, m_disc = median_map(dg)
        rec = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), P=int(len(i)),
                   cos_DIS=float(cosrows[pdb]["cos"]["DIS"]) if pdb in cosrows else float("nan"),
                   fail18=bool(pdb in I.FAIL18), beta={}, tau_meta={})
        #: kappa_a = 2 p_a(median) (A2's local curvature) and w_a (the shipped shell weight):
        #: the theorem's E[<g,r>] = -sum w kappa var(a)(1 - beta) is a w*kappa-WEIGHTED statement,
        #: so the weighted slope is the version the identity implies and is tested beside the plain one.
        _, w_ship = cdf_coeff_at(dg, D0)
        prob = np.asarray(dg["prob"], float); cen = np.asarray(dg["centres"], float)
        kb = np.clip(np.searchsorted(cen, m_grid) - 1, 0, prob.shape[1] - 1)
        kappa = 2.0 * prob[np.arange(prob.shape[0]), kb]
        ww = np.maximum(w_ship * kappa, 0.0)
        for which in taus:
            tau, meta = typical_map(pdb, uu, i, j, which)
            a = D0 - tau; b = m_grid - tau
            beta = float((a @ b) / max(a @ a, 1e-30))
            b2 = m_disc - tau
            rec["beta"][which] = beta
            rec["beta"][which + "_weighted"] = float((ww * a @ b) / max(ww * a @ a, 1e-30))
            sep = (j - i)
            for lo, hi, tag in ((2, 3, "sep23"), (4, 6, "sep46"), (7, 99, "sep7p")):
                msk = (sep >= lo) & (sep <= hi)
                rec["beta"][f"{which}_{tag}"] = (float((a[msk] @ b[msk]) / max(a[msk] @ a[msk], 1e-30))
                                                 if msk.sum() >= 3 else float("nan"))
            rec["beta"][which + "_disc"] = float((a @ b2) / max(a @ a, 1e-30))
            rec["beta"][which + "_corr"] = float(np.corrcoef(a, b)[0, 1])
            rec["beta"][which + "_sd_a"] = float(a.std()); rec["beta"][which + "_sd_b"] = float(b.std())
            rec["tau_meta"][which] = meta
        rows.append(rec)
        if (q + 1) % 30 == 0:
            print(f"  [{q+1}/{len(pdbs)}]", flush=True)
    out = dict(check="S29-L7 clause 2 (beta > 1 and the sign law)", n=len(rows), oracle=True, taus=list(taus))
    cos = np.array([r["cos_DIS"] for r in rows])
    folds = ST.pinned_folds([r["pdb"] for r in rows])
    for which in taus:
        bet = np.array([r["beta"][which] for r in rows])
        pred = -np.sign(bet - 1.0)
        agree = (np.sign(cos) == pred)
        k = int(agree.sum()); n = len(rows)
        se = float(np.sqrt(0.25 / n))
        small = np.argsort(np.abs(bet - 1.0))
        out[which] = dict(
            beta_median=float(np.median(bet)), beta_mean=float(bet.mean()),
            frac_above_1=float((bet > 1).mean()), n_above_1=int((bet > 1).sum()),
            sign_agreement=k / n, n_agree=k,
            coin_toss_ci=[0.5 - 1.96 * se, 0.5 + 1.96 * se],
            clears_two_thirds=bool((bet > 1).mean() >= 2 / 3),
            clears_seventy=bool(k / n >= 0.70),
            agreement_in_smallest_third=float(agree[small[:n // 3]].mean()),
            agreement_in_largest_third=float(agree[small[-(n // 3):]].mean()),
            corr_median=float(np.median([r["beta"][which + "_corr"] for r in rows])),
            beta_disc_median=float(np.median([r["beta"][which + "_disc"] for r in rows])),
            fold_ci_beta=ST.compare(bet, np.zeros(len(bet)), folds, label="beta " + which,
                                    seed_parts=("s29Dth",))["ci95_fold"],
            fail18_beta=float(np.median(bet[[r["fail18"] for r in rows]])),
            other_beta=float(np.median(bet[[not r["fail18"] for r in rows]])),
            one_over_beta_median=float(1.0 / np.median(bet)),
            beta_weighted_median=float(np.median([r["beta"][which + "_weighted"] for r in rows])),
            frac_weighted_above_1=float(np.mean([r["beta"][which + "_weighted"] > 1 for r in rows])),
            sign_agreement_weighted=float(np.mean(np.sign(cos) == -np.sign(
                np.array([r["beta"][which + "_weighted"] for r in rows]) - 1.0))),
            beta_by_separation={t: float(np.nanmedian([r["beta"].get(f"{which}_{t}", np.nan) for r in rows]))
                                for t in ("sep23", "sep46", "sep7p")})
        print(f"  tau={which:6s} beta median {np.median(bet):.3f}  >1 on {int((bet>1).sum())}/{n} "
              f"({100*(bet>1).mean():.0f}%)  sign agreement {k}/{n} ({100*k/n:.0f}%, coin toss CI "
              f"[{100*out[which]['coin_toss_ci'][0]:.0f}, {100*out[which]['coin_toss_ci'][1]:.0f}]%)  "
              f"agreement in the smallest |beta-1| third {100*out[which]['agreement_in_smallest_third']:.0f}% "
              f"vs largest {100*out[which]['agreement_in_largest_third']:.0f}%   1/median(beta) = "
              f"{out[which]['one_over_beta_median']:.3f}")
        print(f"           WEIGHTED (w*kappa) beta median {out[which]['beta_weighted_median']:.3f}, >1 on "
              f"{100*out[which]['frac_weighted_above_1']:.0f}%, sign agreement {100*out[which]['sign_agreement_weighted']:.0f}%"
              f"   | beta by separation (2-3/4-6/7+): " + "/".join("%.2f" % out[which]["beta_by_separation"][t]
                                                                   for t in ("sep23", "sep46", "sep7p")))
    out["rows"] = rows
    return out


# ============================================================ clause 3: the shrink grid
_SHRINK_CACHE = {}


def _shrunk_score(W, ctx, s, tau_which="LIB75"):
    """The shipped Bayes-risk score with the target map shrunk to tau + s (m - tau).

    NATIVE-FREE: tau is the sequence-blind universe draw, m the posterior's own median map.
    Shifting each pair's risk curve by delta = (1 - s)(tau - m) moves that pair's minimiser to
    tau + s(m - tau) while keeping the risk's shape, which is exactly lane T's construction.
    """
    W = np.asarray(W, float)
    if W.ndim == 2:
        W = W[None]
    key = (ctx.pdb, tau_which)
    if key not in _SHRINK_CACHE:
        i, j = I.pair_index(ctx.n)
        m_grid, _ = median_map(ctx.dg)
        tau, _ = typical_map(ctx.pdb, ctx.universe(), i, j, tau_which)
        _SHRINK_CACHE[key] = (i, j, m_grid, tau)
    i, j, m_grid, tau = _SHRINK_CACHE[key]
    delta = (1.0 - float(s)) * (tau - m_grid)
    D = I.pair_dists(W, i, j) - delta[None, :]
    risk = np.asarray(ctx.dg["risk"], float); grid = np.asarray(ctx.dg["grid"], float)
    g0, h = float(grid[0]), float(grid[1] - grid[0])
    G = risk.shape[1]
    x = np.clip((D - g0) / h, 0, G - 1.0000001)
    i0 = np.floor(x).astype(int); f = x - i0
    rows = np.arange(risk.shape[0])[None, :]
    r0 = risk[rows, i0]; r1 = risk[rows, i0 + 1]
    return (r0 + f * (r1 - r0)).mean(1)


def cost_s1_00(W, ctx):
    return _shrunk_score(W, ctx, 1.00)


def cost_s0_90(W, ctx):
    return _shrunk_score(W, ctx, 0.90)


def cost_s0_80(W, ctx):
    return _shrunk_score(W, ctx, 0.80)


def cost_s0_75(W, ctx):
    return _shrunk_score(W, ctx, 0.75)


def cost_s0_70(W, ctx):
    return _shrunk_score(W, ctx, 0.70)


def cost_s0_60(W, ctx):
    return _shrunk_score(W, ctx, 0.60)


def cost_s0_50(W, ctx):
    return _shrunk_score(W, ctx, 0.50)


def cost_s0_40(W, ctx):
    return _shrunk_score(W, ctx, 0.40)


def cost_s0_30(W, ctx):
    return _shrunk_score(W, ctx, 0.30)


_S_FUN = {1.0: "cost_s1_00", 0.9: "cost_s0_90", 0.8: "cost_s0_80", 0.75: "cost_s0_75", 0.7: "cost_s0_70",
          0.6: "cost_s0_60", 0.5: "cost_s0_50", 0.4: "cost_s0_40", 0.3: "cost_s0_30"}


def shrink(pdbs, grid=SHRINK_GRID):
    out = dict(check="S29-L7 clause 3 (the shrink grid: meter number 2 is gameable)",
               n=len(pdbs), oracle=True, grid=list(grid), cells={})
    for s in grid:
        t0 = time.time()
        r = M.run_meter(__name__ + ":" + _S_FUN[s], basis="ca", pdbs=list(pdbs), quiet=True, save=False)
        c = r["cosine"]; p = r["native_pctile"]; sh = c.get("shrink_signature")
        out["cells"][str(s)] = dict(
            s=float(s), cos=c["mean"], cos_se=c["se"], cos_fold_ci=c["ci95_fold"], cos_n=c["n"],
            cos_pos=c["n_pos"], rand_ref=c["random_ref_mean_abs"],
            pctile=p["mean"] if p else None, pctile_fold_ci=p["ci95_fold"] if p else None,
            pref_circ_best=r["pref"]["circ_best"]["mean"], pref_native=r["pref"]["NATIVE"]["mean"],
            rho_S28=r["ladder_rho"]["S28"]["mean"], rho_CHARTER=r["ladder_rho"]["CHARTER"]["mean"],
            bond_ratio=sh["bond_ratio_mean"] if sh else None, rg_ratio=sh["rg_ratio_mean"] if sh else None,
            n_contracting=sh["n_contracting"] if sh else None,
            mean_f_prod=r["mean_f"]["PROD"], secs=time.time() - t0)
        d = out["cells"][str(s)]
        print(f"  s={s:.2f}  cos {d['cos']:+.4f} [{d['cos_fold_ci'][0]:+.3f},{d['cos_fold_ci'][1]:+.3f}] "
              f"({d['cos_pos']}/{d['cos_n']} positive)   native pctile {d['pctile']:.4f}   "
              f"pref(circ_best) {d['pref_circ_best']:.3f}  rho_S28 {d['rho_S28']:+.3f}  "
              f"bond x{d['bond_ratio']:.4f} Rg x{d['rg_ratio']:.4f} ({d['n_contracting']}/{d['cos_n']} contract)  "
              f"{d['secs']:.0f}s", flush=True)
    cs = [out["cells"][str(s)]["cos"] for s in grid]
    ps = [out["cells"][str(s)]["pctile"] for s in grid]
    out["monotone_cos"] = bool(all(cs[k + 1] >= cs[k] for k in range(len(cs) - 1)))
    out["monotone_pctile_worse"] = bool(all(ps[k + 1] >= ps[k] for k in range(len(ps) - 1)))
    out["cos_at_registered"] = {str(s): out["cells"][str(s)]["cos"] for s in REGISTERED_S}
    zc = [s for k, s in enumerate(grid) if cs[k] >= 0]
    out["zero_crossing_between"] = ([float(grid[[q for q in range(len(grid)) if cs[q] < 0][-1]]),
                                     float(max(zc))] if zc and any(c < 0 for c in cs) else None)
    return out


# ============================================================ S29-L11 prediction 2: the variance law
def observables(n, rng):
    """Unit-spectral-norm observables whose stable ranks span 1 to ~D/2.  NOT lane T's."""
    D = 2 ** n
    out = {}
    G = rng.normal(size=(D, D)); G = (G + G.T) / 2.0
    out["GOE_dense"] = G / np.linalg.norm(G, 2)
    d = rng.normal(size=D)
    out["diag_gauss"] = np.diag(d / np.abs(d).max())
    for k in (1, 4, 32):
        idx = rng.choice(D, k, replace=False)
        P = np.zeros((D, D)); P[idx, idx] = 1.0
        out[f"projector_k{k}"] = P
    Z = np.array([[1.0, 0.0], [0.0, -1.0]]); X = np.array([[0.0, 1.0], [1.0, 0.0]])
    P = np.array([1.0])
    for q in range(n):
        P = np.kron(P, Z if q % 2 == 0 else X)
    out["pauli_string"] = P
    c = np.zeros(D)
    for sh in (1, 2, 3, 5, 8):                      # a 10-regular circulant: r_stable = D/10
        c[sh] = 1.0; c[-sh] = 1.0
    Acirc = np.zeros((D, D))
    for q in range(D):
        Acirc[q] = np.roll(c, q)
    out["circulant_10reg"] = Acirc / np.linalg.norm(Acirc, 2)
    v = rng.normal(size=D); v /= np.linalg.norm(v)
    out["rank_one"] = np.outer(v, v)
    return out


def r_stable(Aop):
    return float((np.linalg.norm(Aop, "fro") ** 2) / (np.linalg.norm(Aop, 2) ** 2))


def variance_law(ns=(7, 8, 9), draws=120, seed=20260920, var_diag_const=0.0305):
    from core import quantum as Q
    out = dict(check="S29-L11 prediction 2 (Var = r_stable/D^2; J* = D sqrt(0.0305/r_stable))",
               ns=list(ns), draws=int(draws), seed=int(seed), var_diag_const=var_diag_const,
               oracle=False, note="a property measurement of the deployed ansatz; no native, no RMSD",
               cells=[])
    for n in ns:
        D = 2 ** n
        circ = Q.StatevectorCircuit(n, 3)
        P = circ.n_params()
        rng = np.random.default_rng(seed + n)
        obs = observables(n, rng)
        TH = np.random.default_rng(seed * 7 + n).normal(0.0, 0.6, (draws, P))
        #: the exact parameter-shift gradient of <psi|A|psi> in theta_0, my own implementation
        e0 = np.zeros(P); e0[0] = np.pi / 2.0
        Sp = circ.states_batch(TH + e0[None, :])
        Sm = circ.states_batch(TH - e0[None, :])
        for name, Aop in obs.items():
            rs = r_stable(Aop)
            gp = np.einsum("bi,ij,bj->b", Sp, Aop, Sp)
            gm = np.einsum("bi,ij,bj->b", Sm, Aop, Sm)
            g = 0.5 * (gp - gm)
            var = float(np.var(g, ddof=1))
            pred = rs / D ** 2
            cell = dict(n=n, D=D, observable=name, r_stable=rs, var_measured=var, var_predicted=pred,
                        ratio=float(var / pred) if pred > 0 else float("nan"),
                        within_2x=bool(0.5 <= var / max(pred, 1e-300) <= 2.0),
                        within_3x=bool(1 / 3 <= var / max(pred, 1e-300) <= 3.0),
                        jstar=float(D * np.sqrt(var_diag_const / max(rs, 1e-30))))
            out["cells"].append(cell)
            print(f"  n={n} D={D:5d} {name:16s} r_stable {rs:9.3f}  Var measured {var:.4e}  predicted {pred:.4e}  "
                  f"ratio {cell['ratio']:.3f} {'OK' if cell['within_2x'] else ('2-3x' if cell['within_3x'] else 'FAIL')}"
                  f"   J* {cell['jstar']:.1f}", flush=True)
    c7 = [c for c in out["cells"] if c["n"] >= 7]
    out["n_cells_ge7"] = len(c7)
    out["n_within_2x"] = int(sum(c["within_2x"] for c in c7))
    out["n_within_3x"] = int(sum(c["within_3x"] for c in c7))
    out["worst_ratio"] = float(max(c7, key=lambda c: abs(np.log(c["ratio"])))["ratio"])
    out["verdict"] = ("HOLDS: %d/%d cells at n >= 7 within 2x (all within 3x), worst ratio %.2f"
                      % (out["n_within_2x"], len(c7), out["worst_ratio"])) if out["n_within_3x"] == len(c7) else \
                     ("FALSIFIED: %d/%d cells depart by more than 3x" % (len(c7) - out["n_within_3x"], len(c7)))
    print(" ", out["verdict"])
    return out


# ============================================================ driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["identity", "beta", "shrink", "variance", "all"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tau", default="both", help="LIB75 | UNIV | POOL | both (LIB75+UNIV) | all")
    ap.add_argument("--n", default="7,8,9")
    ap.add_argument("--draws", type=int, default=120)
    a = ap.parse_args()
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    os.makedirs(RESULTS, exist_ok=True)
    modes = ["identity", "beta", "shrink", "variance"] if a.mode == "all" else [a.mode]
    for mode in modes:
        t0 = time.time()
        print(f"== {mode} ==", flush=True)
        if mode == "identity":
            res = identity(pdbs)
        elif mode == "beta":
            taus = ("LIB75", "UNIV") if a.tau == "both" else (("LIB75", "UNIV", "POOL") if a.tau == "all" else (a.tau,))
            res = beta_check(pdbs, taus)
        elif mode == "shrink":
            res = shrink(pdbs)
        else:
            res = variance_law(tuple(int(x) for x in a.n.split(",")), a.draws)
        res["secs"] = time.time() - t0
        ST.save_atomic(os.path.join(RESULTS, f"s29_D_theory_{mode}.json"), res, module_file=__file__)
        print("wrote", os.path.join(RESULTS, f"s29_D_theory_{mode}.json"), f"({res['secs']:.0f}s)")


if __name__ == "__main__":
    main()
