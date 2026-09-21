#!/usr/bin/env python
"""s31/s31_E1_direction.py -- S31 lane E, measurement E1: is `mu_hat` a usable DIRECTION for `mu`?

    mu      = pool75_mean - d_nat      ORACLE / NOT DEPLOYABLE
    mu_hat  = pool75_mean - expected   NATIVE-FREE  (trace proved in s31_E_lib.native_free_trace)

E1 is the load-bearing measurement of the whole lane: if `mu_hat` does not track `mu`'s direction,
E2/E3 are dead.  Bars registered in `s31/PREREG_S31_E.md` §2 BEFORE this file was run:

    within-target mean cosine  c >= 0.707  STRONG   (c^2 >= 1/2 of mu's energy removed)
                               c >= 0.577  USABLE   (c^2 >= 1/3)
                               c >= 0.30   WEAK
                               c <  0.30   DEAD, or not separated from the matched control

TWO cosines are reported and they answer different questions:
  * UNCENTERED cos -- the quantity the E2 projection operator actually uses.  It carries a
    SHARED-REFERENT FLOOR (`shared-referent-floor`): both vectors share `-d_nat`, and both have
    a large common mean, so a random vector with the same marginal already scores high.  The
    permutation control measures that floor; only the EXCESS over it is information.
  * CENTERED corr -- floor-free, and the quantity S30's `coh` uses.

    python s31/s31_E1_direction.py
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

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s31 import s31_E_lib as L             # noqa: E402

OUT = os.path.join(HERE, "results", "s31_E1_direction.json")
SEED = 310501


def cosu(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 1e-12 and nb > 1e-12 else np.nan


def corrc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return cosu(a, b)


def main():
    d = L.load()
    X, y, mu, mh = d["X"], d["y"], d["mu"], d["mu_hat"]
    tgt, sep, fold = d["tgt"], d["sep"], d["fold"]
    pdbs = [str(p) for p in d["pdbs"]]
    tfold = np.array([fold[tgt == q][0] for q in range(len(pdbs))], int)
    rng = np.random.default_rng(SEED)

    res = {"seed": SEED, "n_targets": len(pdbs), "n_pairs": int(len(y)),
           "native_free_trace": L.native_free_trace(),
           "identity_mu_hat_eq_mu_minus_y_max_abs_dev": float(np.max(np.abs(mh - (mu - y)))),
           "prereg_bars": {"STRONG": 0.707, "USABLE": 0.577, "WEAK": 0.30},
           "masks": {}}

    for mname, gm in (("all_pairs_minsep2", np.ones(len(y), bool)),
                      ("LONGRANGE_sep_ge7", sep >= 7)):
        per = {k: [] for k in ("cos", "cos_perm", "cos_gauss", "corr", "corr_perm",
                               "corr_gauss", "coh0", "r_sd", "corr_pred_from_identity",
                               "slope_mu_on_muhat", "sd_y", "sd_mu", "sd_mh", "npair")}
        keep = []
        for q in range(len(pdbs)):
            m = gm & (tgt == q)
            if m.sum() < 8:
                continue
            a, b, yy = mh[m], mu[m], y[m]
            # permutation control: identical marginal AND identical norm, pairing destroyed
            ap = a[rng.permutation(len(a))]
            # isotropic control, norm-matched
            g = rng.standard_normal(len(a)); g *= np.linalg.norm(a) / np.linalg.norm(g)
            per["cos"].append(cosu(a, b)); per["cos_perm"].append(cosu(ap, b))
            per["cos_gauss"].append(cosu(g, b))
            per["corr"].append(corrc(a, b)); per["corr_perm"].append(corrc(ap, b))
            per["corr_gauss"].append(corrc(g, b))
            c0 = corrc(yy, b)
            per["coh0"].append(c0)
            sy, sm = float(yy.std()), float(b.std())
            r = sy / sm if sm > 1e-12 else np.nan
            per["r_sd"].append(r); per["sd_y"].append(sy); per["sd_mu"].append(sm)
            per["sd_mh"].append(float(a.std()))
            den = 1.0 + r * r - 2.0 * c0 * r
            per["corr_pred_from_identity"].append((1.0 - c0 * r) / np.sqrt(den)
                                                  if den > 1e-12 else np.nan)
            bc = b - b.mean(); ac = a - a.mean()
            per["slope_mu_on_muhat"].append(float(ac @ bc / max(ac @ ac, 1e-12)))
            per["npair"].append(int(m.sum()))
            keep.append(q)
        keep = np.array(keep, int)
        V = {k: np.array(v, float) for k, v in per.items()}
        fo = tfold[keep]; nm = [pdbs[q] for q in keep]
        blk = {"n_targets_used": int(len(keep))}

        def level(name, arr):
            c = ST.compare(arr, np.zeros(len(arr)), folds=fo, names=nm, label="E1|%s|%s" % (mname, name))
            return {"mean": float(np.nanmean(arr)), "median": float(np.nanmedian(arr)),
                    "sd": float(np.nanstd(arr, ddof=1)), "se": c["se"], "mde": c["mde"],
                    "ci95_fold": c["ci95_fold"], "per_fold": c["per_fold"],
                    "min": float(np.nanmin(arr)), "max": float(np.nanmax(arr)),
                    "frac_ge_0577": float(np.mean(arr >= 0.577)),
                    "frac_ge_0707": float(np.mean(arr >= 0.707))}

        for k in ("cos", "cos_perm", "cos_gauss", "corr", "corr_perm", "corr_gauss",
                  "coh0", "r_sd", "slope_mu_on_muhat", "corr_pred_from_identity"):
            blk[k] = level(k, V[k])
        # paired excess over the matched controls -- the registered separation test
        for a_, b_, lab in (("cos", "cos_perm", "cos_vs_perm"), ("cos", "cos_gauss", "cos_vs_gauss"),
                            ("corr", "corr_perm", "corr_vs_perm"),
                            ("corr", "corr_gauss", "corr_vs_gauss")):
            c = ST.compare(V[a_], V[b_], folds=fo, names=nm, label="E1|%s|%s" % (mname, lab))
            blk[lab] = {"excess": c["effect"], "se": c["se"], "mde": c["mde"],
                        "x_mde": c["effect_over_mde"], "ci95_fold": c["ci95_fold"],
                        "folds_same_sign": c["folds_same_sign"],
                        "W": c["n_better"], "L": c["n_worse"], "verdict": c["verdict"]}
        # identity self-check: does the measured corr match the closed form from (coh0, r)?
        blk["identity_check_corr"] = {
            "max_abs_dev": float(np.nanmax(np.abs(V["corr"] - V["corr_pred_from_identity"]))),
            "mean_abs_dev": float(np.nanmean(np.abs(V["corr"] - V["corr_pred_from_identity"])))}
        blk["per_target"] = {nm[i]: {k: float(V[k][i]) for k in V} for i in range(len(keep))}
        # the registered verdict, applied to the PRIMARY statistic (uncentered cos, the operator's)
        c = blk["cos"]["mean"]
        sep_ok = abs(blk["cos_vs_perm"]["x_mde"]) >= 1.0 and blk["cos_vs_perm"]["excess"] > 0
        blk["PREREG_VERDICT_cos"] = (
            "DEAD (c < 0.30)" if c < 0.30 else
            "DEAD (not separated from the matched permutation control at 1.0x MDE)" if not sep_ok
            else "STRONG" if c >= 0.707 else "USABLE" if c >= 0.577 else "WEAK")
        cc = blk["corr"]["mean"]
        sep_ok2 = abs(blk["corr_vs_perm"]["x_mde"]) >= 1.0 and blk["corr_vs_perm"]["excess"] > 0
        blk["PREREG_VERDICT_corr"] = (
            "DEAD (c < 0.30)" if cc < 0.30 else
            "DEAD (not separated from the matched permutation control at 1.0x MDE)" if not sep_ok2
            else "STRONG" if cc >= 0.707 else "USABLE" if cc >= 0.577 else "WEAK")
        res["masks"][mname] = blk

    tmp = OUT + ".%d.tmp" % os.getpid()
    with open(tmp, "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, OUT)

    for mname, b in res["masks"].items():
        print("=== %s  (%d targets)" % (mname, b["n_targets_used"]))
        for k in ("cos", "cos_perm", "cos_gauss", "corr", "corr_perm", "corr_gauss", "coh0",
                  "r_sd", "slope_mu_on_muhat"):
            q = b[k]
            print("  %-22s mean %+.4f  median %+.4f  fold95 [%+.4f,%+.4f]  >=.577 %.2f  >=.707 %.2f"
                  % (k, q["mean"], q["median"], q["ci95_fold"][0], q["ci95_fold"][1],
                     q["frac_ge_0577"], q["frac_ge_0707"]))
        for k in ("cos_vs_perm", "cos_vs_gauss", "corr_vs_perm", "corr_vs_gauss"):
            q = b[k]
            print("  %-22s excess %+.4f  %.2fxMDE  fold95 [%+.4f,%+.4f]  ss%d  %s"
                  % (k, q["excess"], abs(q["x_mde"]), q["ci95_fold"][0], q["ci95_fold"][1],
                     q["folds_same_sign"], q["verdict"]))
        print("  identity check |corr - closed form| max %.2e mean %.2e"
              % (b["identity_check_corr"]["max_abs_dev"], b["identity_check_corr"]["mean_abs_dev"]))
        print("  PREREG VERDICT  cos: %s   |   corr: %s"
              % (b["PREREG_VERDICT_cos"], b["PREREG_VERDICT_corr"]))
    print("native-free trace:", res["native_free_trace"])
    print("identity max|mu_hat-(mu-y)| = %.2e" % res["identity_mu_hat_eq_mu_minus_y_max_abs_dev"])
    print("wrote", OUT)


if __name__ == "__main__":
    main()
