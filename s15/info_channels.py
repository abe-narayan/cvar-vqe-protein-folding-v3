"""SPRINT 15, INFO, PART A -- every target-specific channel on one common yardstick.

The brief asks what each channel CAN and CANNOT provide, in numbers, so Phase 2 builds only
what the information supports.  Every channel below is scored on the SAME instrument and the
SAME metric the sprint-14 arithmetic is stated in:

    instrument   the shipped K=500 BLOSUM retrieval pool of each of the 126 tuning targets
    in-band      pool members with `rr <= pool_best + 1.5 A` -- the region where the answer
                 lives and where sprint 14 measured every objective at chance
    metric       IN-BAND PAIRWISE ORDERING ACCURACY (ties averaged), against the 0.500 null
                 and against the 0.638 required for 2.0 A through a top-100 operator
    also         global pair accuracy, and the RMSD emitted by a top-24 / top-75 cut

Channels scored (all NATIVE-FREE):
    disto_bayes     the shipped Bayes-risk distogram score -- the incumbent
    disto_l2        plain squared error against the distogram's `expected` column
    disto_l2_sdw    the same, weighted by 1/sd^2 -- THE `sd` COLUMN HAS NEVER BEEN USED
    disto_l2_lowsd  the same, restricted to the half of the pairs with the smallest sd
    esm_contact     agreement with the ESM-2 contact map (via s12/esm_bank.py, 100 MB)
    pool_typicality mean CA-RMSD to the rest of the pool -- the consensus/medoid criterion
    rg_compact      radius of gyration
    ss_consensus    agreement with the pool-consensus secondary structure
    seq_blind_null  a random score, the calibration null

`sd` calibration is measured separately: is the distogram's own uncertainty column
informative about its own error, and does conditioning on it buy in-band ordering?

    python -m s15.info_channels
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s15 import info_lib as L            # noqa: E402

BAND = 1.5
NPAIR = 200_000


def pair_accuracy(score, rmsd, rng, npair=NPAIR):
    """P(the lower-scoring member of a random pair really is the better one), ties at 0.5."""
    n = len(score)
    if n < 2:
        return float("nan")
    a = rng.integers(0, n, npair); b = rng.integers(0, n, npair)
    k = a != b
    a, b = a[k], b[k]
    ds = score[a] - score[b]
    dr = rmsd[a] - rmsd[b]
    k = dr != 0
    a, b, ds, dr = a[k], b[k], ds[k], dr[k]
    if len(ds) == 0:
        return float("nan")
    hit = np.where(ds == 0, 0.5, (np.sign(ds) == np.sign(dr)).astype(float))
    return float(hit.mean())


def emit_topm(score, rmsd, m):
    m = min(m, len(score))
    thr = np.partition(score, m - 1)[m - 1]
    s = score < thr
    ns = int(s.sum())
    t = score == thr
    if ns >= m:
        return float(rmsd[s].mean())
    w = np.zeros(len(score)); w[s] = 1.0; w[t] = (m - ns) / float(t.sum())
    return float((w * rmsd).sum() / w.sum())


def ss_string(phi, psi):
    return I.ss_of(np.asarray(phi, float), np.asarray(psi, float))


def build_scores(t, esm):
    p, n, fold, seq = t["pdb"], int(t["n"]), int(t["fold"]), t["seq"]
    u = I.load_univ(p)
    pi = I.pool_idx(u)
    W = np.asarray(u["W"], float)[pi]
    rr = np.asarray(u["rr"], float)[pi]
    PH = np.asarray(u["PHI"], float)[pi]
    PS = np.asarray(u["PSI"], float)[pi]
    dg = I.distogram(p, seq, fold)
    i, j = I.pair_index(n)
    D = I.pair_dists(W, i, j)
    exp = np.asarray(dg["expected"], float)
    sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
    S = {}
    S["disto_bayes"] = I.shipped_score(dg, D)
    S["disto_l2"] = ((D - exp[None, :]) ** 2).mean(1)
    S["disto_l2_sdw"] = (((D - exp[None, :]) / sd[None, :]) ** 2).mean(1)
    lo = sd <= np.median(sd)
    S["disto_l2_lowsd"] = ((D[:, lo] - exp[None, lo]) ** 2).mean(1)
    lo25 = sd <= np.quantile(sd, 0.25)
    S["disto_l2_lowsd25"] = ((D[:, lo25] - exp[None, lo25]) ** 2).mean(1)
    # the SHIPPED Bayes-risk score restricted to the confident half of the pairs: the
    # cheapest possible use of the `sd` column inside the existing pipeline
    dg_lo = {"grid": dg["grid"], "risk": np.asarray(dg["risk"])[lo]}
    S["disto_bayes_lowsd"] = I.shipped_score(dg_lo, D[:, lo])
    # ESM-2 contact map agreement: reward a structure whose contacts match the predicted map
    if seq in esm:
        con = esm[seq][2]
        cij = con[i, j]
        S["esm_contact"] = -(cij[None, :] * (D < 8.0)).mean(1)
    P = I.pairwise_rmsd(W)
    S["pool_typicality"] = P.mean(1)
    S["rg_compact"] = L.rg_of(W)
    # pool-consensus secondary structure: per-residue majority over the pool
    ss = np.array([list(ss_string(PH[k], PS[k])) for k in range(len(PH))])
    cons = np.array([max("HEC", key=lambda c: int((ss[:, r] == c).sum()))
                     for r in range(ss.shape[1])])
    S["ss_consensus"] = -(ss == cons[None, :]).mean(1)
    return S, rr, {"n": n, "fold": fold, "pool_best": float(rr.min()),
                   "sd_mean": float(sd.mean()), "npairs": int(len(i))}


def sd_calibration():
    """Is the distogram's `sd` column informative about its own error?  ORACLE-EVALUATED."""
    zs, ae, sds, dev = [], [], [], []
    for t in I.targets():
        p, n, seq, fold = t["pdb"], int(t["n"]), t["seq"], int(t["fold"])
        u = I.load_univ(p)
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(p, seq, fold)
        i, j = I.pair_index(n)
        d = np.linalg.norm(nat[i] - nat[j], axis=1)
        exp = np.asarray(dg["expected"], float)
        sd = np.maximum(np.asarray(dg["sd"], float), 1e-3)
        zs.append((d - exp) / sd); ae.append(np.abs(d - exp)); sds.append(sd)
        dev.append(d - exp)
    z = np.concatenate(zs); a = np.concatenate(ae); s = np.concatenate(sds)
    dv = np.concatenate(dev)
    q = np.quantile(s, np.linspace(0, 1, 6))
    rel = []
    for k in range(5):
        m = (s >= q[k]) & (s <= q[k + 1])
        rel.append({"sd_bin": [float(q[k]), float(q[k + 1])], "n": int(m.sum()),
                    "mean_sd": float(s[m].mean()),
                    "rms_error": float(np.sqrt((dv[m] ** 2).mean())),
                    "mean_abs_error": float(a[m].mean()),
                    "z_sd": float(z[m].std())})
    return {"n_pairs": int(len(z)), "z_mean": float(z.mean()), "z_sd": float(z.std()),
            "corr_sd_abserr": L.pearson(s, a),
            "spearman_sd_abserr": L.spearman(s, a),
            "reliability": rel,
            "note": "z_sd = 1 would mean sd is a calibrated 1-sigma; the reliability table "
                    "shows whether rms error grows with the quoted sd."}


def main():
    from s12 import esm_bank
    esm = esm_bank.load()
    rows = []
    for a, t in enumerate(I.targets()):
        S, rr, meta = build_scores(t, esm)
        band = np.flatnonzero(rr <= rr.min() + BAND)
        rng = np.random.default_rng(1000 + a)
        r = {"pdb": t["pdb"], **meta, "band_n": int(len(band)),
             "band_mean": float(rr[band].mean()), "pool_mean": float(rr.mean()),
             "scores": {}}
        rng2 = np.random.default_rng(5000 + a)
        S["seq_blind_null"] = rng2.normal(0, 1, len(rr))
        for nm, sc in S.items():
            sc = np.asarray(sc, float)
            r["scores"][nm] = {
                "inband_acc": pair_accuracy(sc[band], rr[band], rng) if len(band) > 3
                else float("nan"),
                "global_acc": pair_accuracy(sc, rr, rng),
                "emit24": emit_topm(sc, rr, 24), "emit75": emit_topm(sc, rr, 75)}
        rows.append(r)
        if a % 25 == 0:
            print(f"  {a + 1}/126", flush=True)

    names = sorted({k for r in rows for k in r["scores"]})
    out = {"per_target": rows, "summary": {}, "sd_calibration": sd_calibration()}
    print("\nPART A -- CHANNELS ON THE COMMON IN-BAND YARDSTICK "
          f"(126 targets, band = pool_best + {BAND} A)\n")
    print(f"{'channel':<18}{'n':>5}{'in-band acc':>13}{'CI95':>18}{'global':>9}"
          f"{'emit24':>9}{'emit75':>9}")
    null = np.asarray([r["scores"]["seq_blind_null"]["inband_acc"] for r in rows], float)
    for nm in names:
        v = np.asarray([r["scores"][nm]["inband_acc"] for r in rows
                        if nm in r["scores"]], float)
        g = np.asarray([r["scores"][nm]["global_acc"] for r in rows if nm in r["scores"]],
                       float)
        e24 = np.asarray([r["scores"][nm]["emit24"] for r in rows if nm in r["scores"]],
                         float)
        e75 = np.asarray([r["scores"][nm]["emit75"] for r in rows if nm in r["scores"]],
                         float)
        k = ~np.isnan(v)
        ci = L.boot_ci(v[k], v[k], fn=lambda a, b: float(np.mean(a)))
        pr = I.paired(v[k], np.full(k.sum(), 0.5))
        out["summary"][nm] = {"n": int(k.sum()), "inband_acc": float(np.nanmean(v)),
                              "inband_ci95": [0.5 + pr["ci95"][0], 0.5 + pr["ci95"][1]],
                              "inband_median": float(np.nanmedian(v)),
                              "global_acc": float(np.nanmean(g)),
                              "emit24": float(np.nanmean(e24)),
                              "emit75": float(np.nanmean(e75)),
                              "wl_vs_chance": [int((v[k] > 0.5).sum()),
                                               int((v[k] < 0.5).sum())],
                              "concentration": L.concentration_verdict(v[k] - 0.5)}
        s = out["summary"][nm]
        print(f"{nm:<18}{s['n']:>5}{s['inband_acc']:>13.3f}"
              f"   [{s['inband_ci95'][0]:.3f},{s['inband_ci95'][1]:.3f}]"
              f"{s['global_acc']:>9.3f}{s['emit24']:>9.3f}{s['emit75']:>9.3f}")
    print("\n  required for 2.0 A through a top-100 operator: 0.638   "
          "chance: 0.500   sprint-14 best measured (Legacy): 0.539")

    print("\n-- PAIRED against the shipped `disto_bayes`, with the null-calibrated "
          "concentration verdict --")
    ref_a = np.asarray([r["scores"]["disto_bayes"]["inband_acc"] for r in rows], float)
    ref_e = np.asarray([r["scores"]["disto_bayes"]["emit24"] for r in rows], float)
    keep = ~np.isnan(ref_a)
    out["vs_shipped"] = {}
    for nm in names:
        if nm == "disto_bayes" or not all(nm in r["scores"] for r in rows):
            continue
        v = np.asarray([r["scores"][nm]["inband_acc"] for r in rows], float)
        e = np.asarray([r["scores"][nm]["emit24"] for r in rows], float)
        k = keep & ~np.isnan(v)
        # ACCURACY is higher-is-better while `I.paired` reads lower-is-better, so the
        # accuracy arm is negated before the report and the mean flipped back on print;
        # otherwise the W/L columns come out inverted.
        pa = L.report(-v[k], -ref_a[k], names=[r["pdb"] for r, kk in zip(rows, k) if kk])
        pa = dict(pa, mean_diff=-pa["mean_diff"], median_diff=-pa["median_diff"],
                  ci95=[-pa["ci95"][1], -pa["ci95"][0]])
        pe = L.report(e[k], ref_e[k], names=[r["pdb"] for r, kk in zip(rows, k) if kk])
        out["vs_shipped"][nm] = {"inband_acc": pa, "emit24": pe}
        print(f"  {nm:<20} in-band acc {L.fmt_paired(pa)}")
        print(f"  {'':<20} emitted@24  {L.fmt_paired(pe)}")

    c = out["sd_calibration"]
    print("\n-- THE DISTOGRAM's `sd` COLUMN, USED FOR THE FIRST TIME --")
    print(f"  pairs {c['n_pairs']}, z = (true - expected)/sd:  mean {c['z_mean']:+.3f}, "
          f"sd {c['z_sd']:.3f}   (1.000 = calibrated)")
    print(f"  corr(sd, |error|) = {c['corr_sd_abserr']:+.3f}  "
          f"spearman {c['spearman_sd_abserr']:+.3f}")
    print(f"  {'sd bin':<20}{'n':>9}{'mean sd':>10}{'rms error':>11}{'z sd':>8}")
    for b in c["reliability"]:
        print(f"  [{b['sd_bin'][0]:.2f},{b['sd_bin'][1]:.2f}]{'':<8}{b['n']:>9}"
              f"{b['mean_sd']:>10.2f}{b['rms_error']:>11.2f}{b['z_sd']:>8.2f}")
    L.jwrite("info_channels", out)
    return out


if __name__ == "__main__":
    main()
