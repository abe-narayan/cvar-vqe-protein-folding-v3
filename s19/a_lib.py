"""SPRINT 19, AGENT A -- shared gather for the predictor's error topology.

Everything here is a pure function of the pinned instrument plus the cached leave-fold-out
distograms.  The ONLY ORACLE quantity is `dtrue` (and anything derived from it), and it is
named `dtrue` / `r` / `absr` everywhere so no arm can consume it by accident.

`dhat` reproduces `s18/objceil.py` EXACTLY: the leave-fold-out separation debias from
`s15/distcal.fit_correction(..., "sep")` fitted on the FULL target list (the subset trap,
brief section 7), then clipped at 2.0 A.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, "cache")
for _d in (RESULTS, CACHE):
    os.makedirs(_d, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s15 import distcal as C                 # noqa: E402
from core import predict as PR               # noqa: E402

HYDRO = {a: PR._PROPS[a][2] for a in PR._PROPS}
CHARGE = {a: PR._PROPS[a][3] for a in PR._PROPS}
VOL = {a: PR._PROPS[a][4] for a in PR._PROPS}
HELIX = {a: PR._PROPS[a][0] for a in PR._PROPS}
SHEET = {a: PR._PROPS[a][1] for a in PR._PROPS}


def debias_fns(data, pdbs, folds):
    """Leave-fold-out additive separation debias, fitted on the FULL list (subset trap)."""
    out = {}
    for f in sorted(set(np.asarray(folds).tolist())):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        out[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(np.asarray(sp, float)), sp))
    return out


def consensus_ss(pdb, n):
    """Native-free H/E/C per residue: simplified DSSP of the top-75 circular-mean torsions."""
    from s14 import retprior as R
    PHI, PSI, _ = R.windows(pdb, "top75")
    phi, psi = R.circ_mean(PHI, axis=0), R.circ_mean(PSI, axis=0)
    s = I.ss_of(phi, psi)
    s = (s + "C" * n)[:n]
    return s


def gather_all(tg=None, with_ss=True):
    """One dict per target with native-free features and the ORACLE residual."""
    tg = tg if tg is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    data = C.gather(tg)
    deb = debias_fns(data, pdbs, folds)

    out = {}
    for t in tg:
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        d = data[pdb]
        i, j, sep = d["i"], d["j"], d["sep"]
        dg = I.distogram(pdb, seq, fold)
        prob = np.asarray(dg["prob"], float)
        centres = np.asarray(dg["centres"], float)
        sd = d["sd"]
        dhat = np.maximum(d["dhat"] - deb[fold](sep), 2.0)
        dtrue = d["dtrue"]

        # --- distribution shape (native-free)
        ent = -(prob * np.log(np.maximum(prob, 1e-12))).sum(1)
        mx = prob.max(1)
        # modes: local maxima over the bin axis carrying > 2% mass
        pad = np.concatenate([np.full((len(prob), 1), -1.0), prob,
                              np.full((len(prob), 1), -1.0)], axis=1)
        ismode = (pad[:, 1:-1] > pad[:, :-2]) & (pad[:, 1:-1] >= pad[:, 2:]) & (prob > 0.02)
        nmode = ismode.sum(1)
        # mass within 1 A of the reported mean -- the coordinator's moment-collapse statistic
        mean_ = (prob * centres).sum(1)
        near = (prob * (np.abs(centres[None] - mean_[:, None]) <= 1.0)).sum(1)

        # --- sequence / position (native-free)
        term = np.minimum(i, n - 1 - j).astype(float)
        hyd = np.array([HYDRO.get(a, 0.0) for a in seq])
        chg = np.array([CHARGE.get(a, 0.0) for a in seq])
        vol = np.array([VOL.get(a, 0.0) for a in seq])
        hel = np.array([HELIX.get(a, 1.0) for a in seq])
        she = np.array([SHEET.get(a, 1.0) for a in seq])

        # --- predicted packing (native-free): neighbours within 8 A in the PREDICTED matrix
        M = np.full((n, n), 99.0)
        M[i, j] = dhat
        M[j, i] = dhat
        np.fill_diagonal(M, 0.0)
        deg = (M < 8.0).sum(1).astype(float)

        rec = {
            "pdb": pdb, "n": n, "fold": fold, "seq": seq,
            "i": i, "j": j, "sep": sep.astype(float),
            "dhat": dhat, "sd": sd, "prob": prob, "centres": centres,
            "ent": ent, "maxp": mx, "nmode": nmode.astype(float), "near1": near,
            "term": term, "w": 1.0 / sd ** 2,
            "hyd_i": hyd[i], "hyd_j": hyd[j], "hyd_m": 0.5 * (hyd[i] + hyd[j]),
            "chg_prod": chg[i] * chg[j], "vol_m": 0.5 * (vol[i] + vol[j]),
            "hel_m": 0.5 * (hel[i] + hel[j]), "she_m": 0.5 * (she[i] + she[j]),
            "deg_m": 0.5 * (deg[i] + deg[j]),
            "pred_contact": (dhat < 8.0).astype(float),
            "nat": d["nat"],
            # ORACLE from here down
            "dtrue": dtrue, "r": dhat - dtrue, "absr": np.abs(dhat - dtrue),
        }
        if with_ss:
            ss = consensus_ss(pdb, n)
            code = {"H": 0.0, "E": 1.0, "C": 2.0}
            sc = np.array([code.get(c, 2.0) for c in ss])
            rec["ss_i"] = sc[i]
            rec["ss_j"] = sc[j]
            rec["ss_same"] = (sc[i] == sc[j]).astype(float)
            rec["ss_H"] = ((sc[i] == 0) & (sc[j] == 0)).astype(float)
            rec["ss_E"] = ((sc[i] == 1) & (sc[j] == 1)).astype(float)
        out[pdb] = rec
    return out, pdbs, folds


def boot(diff, rng, B=4000):
    d = np.asarray(diff, float)
    k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report_pair(name, a, b, rng, folds=None):
    """a - b, negative = a better."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m, lo, hi = boot(a - b, rng)
    d = a - b
    s = (f"{name:<26}{a.mean():>8.3f}  d={m:+.3f} [{lo:+.3f},{hi:+.3f}]"
         f"  med {np.median(d):+.3f}  {int((d<0).sum())}W/{int((d>0).sum())}L")
    if folds is not None:
        folds = np.asarray(folds)
        sg = sum(1 for f in np.unique(folds) if np.sign(d[folds == f].mean()) == np.sign(m))
        s += f"  folds {sg}/{len(np.unique(folds))}"
    return s, {"mean_a": float(a.mean()), "mean_b": float(b.mean()), "diff": m,
               "ci": [lo, hi], "median": float(np.median(d)),
               "W": int((d < 0).sum()), "L": int((d > 0).sum())}
