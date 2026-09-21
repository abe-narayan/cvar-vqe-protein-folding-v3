#!/usr/bin/env python
"""s31/s31_B2_inpool.py -- lane B, B2's core measurement: does ANY torsion channel rank INSIDE
the candidate pool, and does its escape from G1 (its ODD part) carry the skill?

Registered in s31/PREREG_S31_B.md (committed e2fc6257, before the first number). Measurements
M0-M4 and M6 of that document.

THE 2x2 THAT MOTIVATES IT
    For an ideal-geometry backbone the point reflection R acts on torsions EXACTLY as
    (phi, psi) -> (-phi, -psi). So T_even = 1/2[T(phi,psi) + T(-phi,-psi)] is a reflection-
    invariant single-structure observable and by G1 is a function of the distance map -- the
    class S30 closed. ALL of a torsion channel's escape from G1 lives in T_odd. Crossed with
    separability (S30-R: a sum of per-residue terms cannot see a lever arm) the only cell no
    existing theorem closes is CHIRAL AND NON-SEPARABLE. M6 builds one.

CHANNELS
    LEG_torsion   0.15 * sum_i rama_penalty(aa_i, phi_i, psi_i)   -- SEPARABLE, partly chiral
    RAMA          the S27 fold-conditioned Ramachandran channel   -- SEPARABLE
    DIS           the SHIPPED cost, the incumbent to partial out
    HELIX_CONST   RMS angular deviation from a constant alpha-helix (-57, -47).
                  The ZERO-TARGET-INFORMATION control. Contract rule 9: plausible, never
                  uniform-on-the-torus, which is a WORSE measure rather than an uninformative one.
    XTWIST[s]     mean sin(crossing dihedral) over spatially contacting segment pairs at
                  sequence separation >= s   -- CHIRAL and NON-SEPARABLE  (M6, the open cell)
    XTWABS[s]     the same functional with |chi| instead of sin(chi) -- the ACHIRAL TWIN of
                  XTWIST, i.e. a matched control in the operator's own space (contract rule 7)
    RG            radius of gyration, partialled out as the size/compactness confound

ORACLE CONTENT. `rr` (per-candidate CA-RMSD to native) is ORACLE and labels every correlation
reported here. NOTHING in this file tunes a deployable parameter. Basis: CA POINT CLOUD per
candidate -- NOT the 3.2105 A built-chain endpoint and not the 3.0483 A cloud endpoint.

USAGE
    python s31/s31_B2_inpool.py [--limit N] [--members 500]
"""
from __future__ import annotations

import argparse
import io
import json
import math
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s31_B2_inpool.json")
ROWS = os.path.join(RESULTS, "s31_B2_inpool_rows.jsonl")
S27_CACHE = os.path.join(ROOT, "s27", "cache")

HELIX_PHI, HELIX_PSI = -57.0, -47.0
XT_SEPS = (4, 8)
XT_RCUT = 10.0
N_SHUFFLE = 200


# ============================================================================ channels
def torsion_leg(seq, PHI, PSI):
    """0.15 * sum_i rama_penalty -- byte-identical to s27/cache's LEG_torsion (verified)."""
    from core import energy as et
    w = float(et.DEFAULT_WEIGHTS["torsion"])
    B, n = PHI.shape
    out = np.zeros(B)
    for i in range(n):
        aa = seq[i]
        # rama_penalty is memoised on (aa, phi, psi); the pool's torsions are continuous, so
        # the loop is the honest path. Vectorised reimplementation would risk drifting from it.
        out += np.array([et.rama_penalty(aa, float(p), float(q))
                         for p, q in zip(PHI[:, i], PSI[:, i])])
    return w * out


def helix_const(PHI, PSI):
    """RMS angular deviation (deg) from a CONSTANT alpha-helix. Zero TARGET information: it has
    no target argument at all beyond the candidate. Non-uniform, hence a plausible control."""
    dp = np.degrees(np.angle(np.exp(1j * (PHI - math.radians(HELIX_PHI)))))
    ds = np.degrees(np.angle(np.exp(1j * (PSI - math.radians(HELIX_PSI)))))
    return np.sqrt((dp ** 2 + ds ** 2).mean(1))


def xtwist(W, sep, rcut=XT_RCUT):
    """(XTWIST, XTWABS, n_pairs) over a (B, n, 3) CA stack.

    v_i = C[i+1]-C[i], m_i = midpoint. For segment pairs (i, j) with j-i >= sep and
    |m_i - m_j| < rcut:
        chi_ij = atan2( <v_i x v_j, rhat_ij>, <v_i, v_j> )
    Under the point reflection all difference vectors negate, so v x v is INVARIANT while rhat
    NEGATES: chi -> -chi. Hence sin(chi) is EXACTLY ODD and |chi| EXACTLY EVEN. Non-separable by
    construction: every term couples two segments at sequence separation >= sep.
    """
    W = np.asarray(W, float)
    B, n, _ = W.shape
    if n - 1 < sep + 1:
        return np.full(B, np.nan), np.full(B, np.nan), 0
    V = W[:, 1:, :] - W[:, :-1, :]                       # (B, n-1, 3)
    M = 0.5 * (W[:, 1:, :] + W[:, :-1, :])
    ii, jj = np.triu_indices(n - 1, k=sep)
    if len(ii) == 0:
        return np.full(B, np.nan), np.full(B, np.nan), 0
    vi, vj = V[:, ii, :], V[:, jj, :]
    r = M[:, jj, :] - M[:, ii, :]
    dr = np.linalg.norm(r, axis=2)
    rhat = r / np.maximum(dr, 1e-12)[..., None]
    cr = np.cross(vi, vj)
    chi = np.arctan2((cr * rhat).sum(2), (vi * vj).sum(2))
    m = dr < rcut
    cnt = m.sum(1)
    # A member with NO pair inside rcut must not nan the whole target (it did in the first
    # run, which silently dropped XTWIST from the 500-member band entirely). Fall back to that
    # member's own CLOSEST pair, which is always defined.
    bad = cnt == 0
    if bad.any():
        m = m.copy()
        m[bad, np.argmin(dr[bad], axis=1)] = True
        cnt = m.sum(1)
    s = (np.sin(chi) * m).sum(1) / cnt
    e = (np.abs(chi) * m).sum(1) / cnt
    return s, e, int(np.median(cnt))


def rg_of(W):
    W = np.asarray(W, float)
    return np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1))


# ============================================================================ statistics
def rankz(x):
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(r.std(), 1e-12)


def spear(x, y):
    a, b = rankz(x), rankz(y)
    return float((a * b).mean())


def partial(x, y, z):
    """Partial Spearman rho(x, y | z), within one target."""
    a, b, c = rankz(x), rankz(y), rankz(z)
    rxy, rxz, ryz = float((a * b).mean()), float((a * c).mean()), float((b * c).mean())
    den = math.sqrt(max(1e-12, (1 - rxz ** 2) * (1 - ryz ** 2)))
    return float((rxy - rxz * ryz) / den)


def shuffle_null(x, y, n=N_SHUFFLE, seed=0):
    """Mean |rho| under a WITHIN-TARGET permutation of the ORACLE labels."""
    rng = np.random.default_rng(seed)
    a = rankz(x)
    b = rankz(y)
    v = np.empty(n)
    for k in range(n):
        v[k] = abs(float((a * rng.permutation(b)).mean()))
    return float(v.mean()), float(np.percentile(v, 95))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--members", type=int, default=500)
    a = ap.parse_args()

    tg = I.targets()
    if a.limit:
        tg = tg[: a.limit]

    chan_names = (["LEG_tors", "LEG_tors_even", "LEG_tors_odd", "RAMA", "HELIX_CONST", "RG"]
                  + ["XTWIST_%d" % s for s in XT_SEPS] + ["XTWABS_%d" % s for s in XT_SEPS])
    bands = ("p500", "sub75")
    rho = {(c, b): [] for c in chan_names for b in bands}
    prho = {(c, b): [] for c in chan_names for b in bands}     # partialled on DIS
    prho_rg = {(c, b): [] for c in chan_names for b in bands}  # partialled on RG
    null = {(c, b): [] for c in chan_names for b in bands}
    dis_rho = {b: [] for b in bands}
    m0, m1, m2, m5 = [], [], [], []
    pdbs, folds = [], []
    t0 = time.time()
    fh_rows = io.open(ROWS, "w", encoding="utf-8")

    for ti, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        idx = I.pool_idx(u)[: a.members]
        PHI, PSI = u["PHI"][idx], u["PSI"][idx]
        W = u["W"][idx]
        rr = u["rr"][idx]                                   # ORACLE
        nat = u["nat_ca"]
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        sub = sub[sub < len(idx)]

        with np.load(os.path.join(S27_CACHE, f"{pdb}.npz")) as z:
            DIS = np.asarray(z["DIS"], float)[: len(idx)]
            LT = np.asarray(z["LEG_torsion"], float)[: len(idx)]
            RM = np.asarray(z["RAMA"], float)[: len(idx)]

        # -- M1 mirror: (phi, psi) -> (-phi, -psi) ------------------------------------
        LTm = torsion_leg(seq, -PHI, -PSI)
        LTe, LTo = 0.5 * (LT + LTm), 0.5 * (LT - LTm)

        # -- M0 mirror gap: ORACLE, diagnostic of chiral dynamic range ----------------
        rr_m = I.kabsch_rmsd_batch(2.0 * W.mean(1, keepdims=True) - W, nat)   # ORACLE
        gap = np.abs(rr - rr_m)
        q = np.percentile(rr, [25, 75])
        m0.append(dict(pdb=pdb, fold=fold, n=n,
                       mean_gap=float(gap.mean()), median_gap=float(np.median(gap)),
                       pool_iqr=float(q[1] - q[0]), pool_sd=float(rr.std()),
                       gap_over_sd=float(gap.mean() / max(rr.std(), 1e-9))))

        ch = {"LEG_tors": LT, "LEG_tors_even": LTe, "LEG_tors_odd": LTo, "RAMA": RM,
              "HELIX_CONST": helix_const(PHI, PSI), "RG": rg_of(W)}
        xtc = {}
        for s in XT_SEPS:
            xs, xe, cnt = xtwist(W, s)
            ch["XTWIST_%d" % s] = xs
            ch["XTWABS_%d" % s] = xe
            xtc[s] = cnt

        # -- M1 variance shares -------------------------------------------------------
        vt = float(np.var(LT))
        m1.append(dict(pdb=pdb, fold=fold, n=n, var_total=vt,
                       odd_share=float(np.var(LTo) / vt) if vt > 1e-18 else float("nan"),
                       even_share=float(np.var(LTe) / vt) if vt > 1e-18 else float("nan"),
                       xt_pairs={str(k): v for k, v in xtc.items()}))

        # -- M2 saturation: the within-pool spread against the triage gap -------------
        lp = os.path.join(ROOT, "s30", "results", "s30_D_ladder_structs", f"{pdb}.npz")
        if os.path.exists(lp):
            with np.load(lp, allow_pickle=False) as zz:
                if "chain_phi_PROD" in zz.files and "chain_phi_RAND_SIGNED_0" in zz.files:
                    rungs = {k: torsion_leg(seq, np.asarray(zz["chain_phi_" + k], float)[None],
                                            np.asarray(zz["chain_psi_" + k], float)[None])[0]
                             for k in ("PROD", "circ_best", "RAND_SIGNED_0", "NATIVE")
                             if "chain_phi_" + k in zz.files}
                    sd = float(LT.std())
                    m2.append(dict(pdb=pdb, fold=fold, pool_mean=float(LT.mean()), pool_sd=sd,
                                   pool_p5=float(np.percentile(LT, 5)),
                                   pool_p95=float(np.percentile(LT, 95)),
                                   rungs={k: float(v) for k, v in rungs.items()},
                                   z_rand=float((rungs.get("RAND_SIGNED_0", np.nan) - LT.mean()) / max(sd, 1e-9)),
                                   z_prod=float((rungs.get("PROD", np.nan) - LT.mean()) / max(sd, 1e-9)),
                                   z_best=float((rungs.get("circ_best", np.nan) - LT.mean()) / max(sd, 1e-9))))

        # -- M5 COHERENCE with the pool common mode -----------------------------------
        # DERIVATION being tested. Let the readout be an AFFINE combination of pool members,
        # C = sum_m a_m x_m with sum_m a_m = 1 (argmin is a_m = delta, the uniform average is
        # a_m = 1/75). To first order in the pair-distance map,
        #     e_p = d_C,p - d_nat,p = sum_m a_m (d_m,p - d_nat,p) = mu_p + sum_m a_m eta_m,p
        # because sum_m a_m = 1 passes mu_p through with coefficient EXACTLY ONE, whatever the
        # weights are. So coh = corr(e, mu) is a function of the CONCENTRATION of a, and NOT of
        # the ranker that produced it. PREDICTION: every ranker gives the same coh at the same
        # cardinality, the uniform average gives coh = 1 exactly, and only a readout that leaves
        # the pool's affine hull can move it.
        if len(sub) >= 10:
            pi, pj = I.pair_index(n, 2)
            if len(pi) >= 5:
                Dp = I.pair_dists(W[sub], pi, pj)                      # (75, P)
                dn = I.pair_dists(nat[None], pi, pj)[0]                # ORACLE
                mu = (Dp - dn[None, :]).mean(0)                        # ORACLE common mode
                sdmu = float(mu.std())

                def coh_of(e):
                    if sdmu < 1e-12 or float(np.std(e)) < 1e-12:
                        return float("nan")
                    return float(np.corrcoef(e, mu)[0, 1])

                rk = {"DIS": DIS[sub], "LEG_tors": LT[sub], "RAMA": RM[sub],
                      "HELIX_CONST": ch["HELIX_CONST"][sub], "LEG_tors_odd": LTo[sub]}
                ent = dict(pdb=pdb, fold=fold, sd_mu=sdmu)
                for nm, s_ in rk.items():
                    ent["argmin_" + nm] = coh_of(Dp[int(np.argmin(s_))] - dn)
                rng = np.random.default_rng(abs(hash(pdb)) % (2 ** 31))
                ent["argmin_RANDOM"] = float(np.mean(
                    [coh_of(Dp[int(k)] - dn) for k in rng.integers(0, len(sub), 16)]))
                ent["uniform_mean_pairspace"] = coh_of(Dp.mean(0) - dn)
                Cavg, _ = I.coordinate_average(W[sub])
                ent["coordinate_average"] = coh_of(I.pair_dists(Cavg[None], pi, pj)[0] - dn)
                ent["ORACLE_best"] = coh_of(Dp[int(np.argmin(rr[sub]))] - dn)
                m5.append(ent)

        # -- M3 / M4 ------------------------------------------------------------------
        row = dict(pdb=pdb, fold=fold, n=n)
        for b, sel in (("p500", np.arange(len(idx))), ("sub75", sub)):
            if len(sel) < 10:
                continue
            y, d = rr[sel], DIS[sel]
            dis_rho[b].append(spear(d, y))
            for c in chan_names:
                x = ch[c][sel]
                fin = np.isfinite(x) & np.isfinite(y) & np.isfinite(d)
                # drop only the offending MEMBERS, never the whole target
                if fin.sum() < 10 or np.std(x[fin]) < 1e-12:
                    rho[(c, b)].append(np.nan); prho[(c, b)].append(np.nan)
                    prho_rg[(c, b)].append(np.nan); null[(c, b)].append(np.nan)
                    continue
                xf, yf, df, gf = x[fin], y[fin], d[fin], ch["RG"][sel][fin]
                rho[(c, b)].append(spear(xf, yf))
                prho[(c, b)].append(partial(xf, yf, df))
                prho_rg[(c, b)].append(partial(xf, yf, gf))
                mn, _ = shuffle_null(xf, yf, seed=abs(hash((pdb, c, b))) % (2 ** 31))
                null[(c, b)].append(mn)
                row["%s.%s" % (c, b)] = rho[(c, b)][-1]
                row["p.%s.%s" % (c, b)] = prho[(c, b)][-1]
        row["DIS.p500"] = dis_rho["p500"][-1] if dis_rho["p500"] else None
        row["DIS.sub75"] = dis_rho["sub75"][-1] if dis_rho["sub75"] else None
        fh_rows.write(json.dumps(row, default=float) + "\n")
        pdbs.append(pdb); folds.append(fold)
        if (ti + 1) % 10 == 0:
            print("  %3d/%d  %.0fs" % (ti + 1, len(tg), time.time() - t0), flush=True)
    fh_rows.close()

    folds = np.asarray(folds, int)
    nT = len(pdbs)
    ncomp = 0

    def verdict(v, lab):
        """Arm-vs-zero on the per-target statistic. Higher |v| = more skill; the sign is stated."""
        nonlocal ncomp
        v = np.asarray(v, float)
        m = np.isfinite(v)
        if m.sum() < 20:
            return dict(n=int(m.sum()), note="too few finite targets")
        o = ST.compare(-v[m], np.zeros(int(m.sum())), folds=folds[m], names=list(np.array(pdbs)[m]),
                       label=lab, seed_parts=("s31B", "inpool"))
        ncomp += 1
        return dict(n=int(m.sum()), mean=float(v[m].mean()), median=float(np.median(v[m])),
                    se=o["se"], mde=o["mde"],
                    x_mde=float(v[m].mean() / o["mde"]) if o["mde"] > 0 else float("nan"),
                    ci95_fold=[-o["ci95_fold"][1], -o["ci95_fold"][0]],
                    folds_same_sign=o["folds_same_sign"],
                    per_fold={k: -val for k, val in o["per_fold"].items()})

    out = dict(n_targets=nT, members=a.members, bands=list(bands),
               basis="per-candidate CA point cloud RMSD to native (ORACLE labels); NOT the endpoint",
               oracle="ORACLE / NOT DEPLOYABLE -- rr and the mirror gap both read the native",
               prereg="s31/PREREG_S31_B.md", provenance=ST.provenance(__file__))

    # M0
    g = np.array([r["gap_over_sd"] for r in m0], float)
    out["M0_mirror_gap"] = dict(
        mean_gap_over_pool_sd=float(g.mean()), median=float(np.median(g)),
        p10=float(np.percentile(g, 10)), p90=float(np.percentile(g, 90)),
        mean_gap_A=float(np.mean([r["mean_gap"] for r in m0])),
        mean_pool_sd_A=float(np.mean([r["pool_sd"] for r in m0])), rows=m0,
        note=("DIAGNOSTIC of chiral dynamic range, NOT a ceiling: D(x) determines x up to "
              "reflection and the real pool holds no mirror pairs, so on THIS pool an achiral "
              "observable is not information-limited. G1 forbids a new CHANNEL, not a good "
              "function of the old one."))

    # M1
    os_ = np.array([r["odd_share"] for r in m1], float)
    out["M1_parity"] = dict(leg_torsion_odd_share_mean=float(np.nanmean(os_)),
                            median=float(np.nanmedian(os_)),
                            p10=float(np.nanpercentile(os_, 10)),
                            p90=float(np.nanpercentile(os_, 90)),
                            G1_check_fired=bool(np.nanmean(os_) < 0.10), rows=m1)

    # M2
    if m2:
        zr = np.array([r["z_rand"] for r in m2], float)
        zp = np.array([r["z_prod"] for r in m2], float)
        zb = np.array([r["z_best"] for r in m2], float)
        out["M2_saturation"] = dict(
            z_RAND_SIGNED_mean=float(np.nanmean(zr)), z_RAND_SIGNED_median=float(np.nanmedian(zr)),
            z_PROD_mean=float(np.nanmean(zp)), z_circ_best_mean=float(np.nanmean(zb)),
            frac_rand_outside_pool_p95=float(np.nanmean(zr > 1.645)),
            note="z = (rung LEG_torsion - pool mean) / pool sd. ORACLE: circ_best is native-built.",
            rows=m2)

    # M5 -- coherence with the pool common mode
    if m5:
        keys = [k for k in m5[0] if k not in ("pdb", "fold", "sd_mu")]
        out["M5_coherence"] = {"n": len(m5), "admission_bar_S30": 0.6931,
                               "note": ("coh = within-target corr(readout pair error, pool "
                                        "common-mode pair error mu). ORACLE / NOT DEPLOYABLE: "
                                        "both arguments need the native."),
                               "rows": m5}
        for k in keys:
            v = np.array([r[k] for r in m5], float)
            v = v[np.isfinite(v)]
            out["M5_coherence"][k] = dict(mean=float(v.mean()), median=float(np.median(v)),
                                          sd=float(v.std()), n=int(len(v)),
                                          admitted=bool(v.mean() < 0.6931))

    # M3 / M4
    out["M3_rho"] = {}
    out["M4_partial_DIS"] = {}
    out["M4_partial_RG"] = {}
    out["null_absrho"] = {}
    for b in bands:
        out["M3_rho"]["DIS." + b] = verdict(dis_rho[b], "DIS_rho_" + b)
        for c in chan_names:
            out["M3_rho"]["%s.%s" % (c, b)] = verdict(rho[(c, b)], "rho_%s_%s" % (c, b))
            out["M4_partial_DIS"]["%s.%s" % (c, b)] = verdict(prho[(c, b)], "prho_%s_%s" % (c, b))
            out["M4_partial_RG"]["%s.%s" % (c, b)] = verdict(prho_rg[(c, b)], "prhoRG_%s_%s" % (c, b))
            nv = np.asarray(null[(c, b)], float)
            out["null_absrho"]["%s.%s" % (c, b)] = float(np.nanmean(nv)) if np.isfinite(nv).any() else None

    # G5: the zero-information control against each informative arm, PAIRED
    out["G5_vs_helix_control"] = {}
    for b in bands:
        hc = np.asarray(rho[("HELIX_CONST", b)], float)
        for c in chan_names:
            if c in ("HELIX_CONST", "RG"):
                continue
            v = np.asarray(rho[(c, b)], float)
            m = np.isfinite(v) & np.isfinite(hc)
            if m.sum() < 20:
                continue
            o = ST.compare(-np.abs(v[m]), -np.abs(hc[m]), folds=folds[m],
                           names=list(np.array(pdbs)[m]), label="G5_%s_%s" % (c, b),
                           seed_parts=("s31B", "inpool"))
            ncomp += 1
            out["G5_vs_helix_control"]["%s.%s" % (c, b)] = dict(
                mean_abs_channel=float(np.abs(v[m]).mean()), mean_abs_control=float(np.abs(hc[m]).mean()),
                effect_channel_minus_control=float(np.abs(v[m]).mean() - np.abs(hc[m]).mean()),
                mde=o["mde"], x_mde=float((np.abs(v[m]).mean() - np.abs(hc[m]).mean()) / o["mde"]) if o["mde"] > 0 else float("nan"),
                ci95_fold=[-o["ci95_fold"][1], -o["ci95_fold"][0]],
                folds_same_sign=o["folds_same_sign"])

    out["comparisons_emitted"] = ncomp
    with io.open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)

    # ---------------------------------------------------------------- print
    print("\nn targets %d   %.0fs   comparisons %d" % (nT, time.time() - t0, ncomp))
    print("\nM0  MIRROR GAP (ORACLE, diagnostic not a ceiling)")
    print("    mean |rmsd(x) - rmsd(Rx)| = %.4f A;  mean pool sd = %.4f A;  ratio = %.3f"
          % (out["M0_mirror_gap"]["mean_gap_A"], out["M0_mirror_gap"]["mean_pool_sd_A"],
             out["M0_mirror_gap"]["mean_gap_over_pool_sd"]))
    print("\nM1  PARITY of LEG_torsion over the pool")
    print("    odd (chiral) variance share: mean %.4f  median %.4f  [p10 %.4f p90 %.4f]  "
          "G1-check fired = %s"
          % (out["M1_parity"]["leg_torsion_odd_share_mean"], out["M1_parity"]["median"],
             out["M1_parity"]["p10"], out["M1_parity"]["p90"], out["M1_parity"]["G1_check_fired"]))
    if "M2_saturation" in out:
        s = out["M2_saturation"]
        print("\nM2  SATURATION of LEG_torsion (z against the pool's own spread)")
        print("    RAND_SIGNED z = %+.3f (median %+.3f)   PROD z = %+.3f   circ_best z = %+.3f"
              % (s["z_RAND_SIGNED_mean"], s["z_RAND_SIGNED_median"], s["z_PROD_mean"],
                 s["z_circ_best_mean"]))
        print("    share of targets with RAND_SIGNED above the pool's 95th percentile: %.3f"
              % s["frac_rand_outside_pool_p95"])
    for b in bands:
        print("\nM3/M4  BAND = %s      (rho > 0 means LOWER cost -> LOWER RMSD, i.e. correct order)" % b)
        print("    %-16s %8s %8s %7s %6s  %8s %7s %6s  %8s" %
              ("channel", "rho", "xMDE", "folds", "null", "rho|DIS", "xMDE", "folds", "rho|RG"))
        for c in ["DIS"] + chan_names:
            k = "%s.%s" % (c, b)
            r = out["M3_rho"].get(k) or {}
            p = out["M4_partial_DIS"].get(k) or {}
            g2 = out["M4_partial_RG"].get(k) or {}
            if "mean" not in r:
                continue
            print("    %-16s %+8.4f %8.2f %5d/5 %6.3f  %+8.4f %7.2f %5d/5  %+8.4f" % (
                c, r["mean"], r["x_mde"], r["folds_same_sign"],
                out["null_absrho"].get(k) or float("nan"),
                p.get("mean", float("nan")), p.get("x_mde", float("nan")),
                p.get("folds_same_sign", 0), g2.get("mean", float("nan"))))
    if "M5_coherence" in out:
        print("\nM5  COHERENCE with the pool common mode (ORACLE / NOT DEPLOYABLE)")
        print("    S30 admission bar: ADMIT a corrector iff coh < 0.6931; uncorrected = 0.6931")
        for k, v in sorted(((k, v) for k, v in out["M5_coherence"].items()
                            if isinstance(v, dict) and "mean" in v),
                           key=lambda kv: -kv[1]["mean"]):
            print("    %-28s coh = %.4f  (median %.4f, sd %.4f, n %d)  admitted=%s"
                  % (k, v["mean"], v["median"], v["sd"], v["n"], v["admitted"]))
    print("\nG5  |rho| channel MINUS |rho| constant-alpha-helix control (paired, fold-clustered)")
    for k, v in out["G5_vs_helix_control"].items():
        print("    %-24s %+8.4f   %5.2fx MDE   %d/5 folds   CI [%+.4f, %+.4f]" % (
            k, v["effect_channel_minus_control"], v["x_mde"], v["folds_same_sign"],
            v["ci95_fold"][0], v["ci95_fold"][1]))
    print("\nwrote", OUT, "\nwrote", ROWS, flush=True)


if __name__ == "__main__":
    main()
