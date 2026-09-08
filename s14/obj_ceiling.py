"""SPRINT 14, OBJ -- is the in-band failure a TRANSFER problem or a SIGNAL problem?

A leave-fold-out model that cannot order near-native structures has two possible
diseases, and they have completely different prognoses:

  TRANSFER  the features resolve the in-band region, but the mapping does not carry
            across targets.  More targets / better conditioning / a bigger model help.
  SIGNAL    the features do not resolve the in-band region AT ALL.  Nothing helps, and
            the learning curve will be flat by construction.

Three arms separate them, on the near-native band of each enumerated target:

  same_insample   fit on the target's own band, scored on the SAME rows.  Wildly
                  overfit on purpose: it upper-bounds what this feature space can
                  express about this target at all.
  same_heldout    fit on half the target's band, scored on the OTHER half of the SAME
                  target.  No transfer is required -- only within-target generalisation.
  crosstarget     fit on the OTHER targets' bands, scored on this target's band.
                  This is the real problem.

If `same_heldout` is at chance, the disease is SIGNAL and no training-set size fixes it.

`band_max` selects the rows: RMSD below a threshold, i.e. exactly the region where the
1.594 A answer lives and where section 1.2 showed every physical objective is at chance.
Band membership uses the native (ORACLE DIAGNOSTIC) -- this measures a ceiling, it is not
an inference-time procedure.

    python -m s14.obj_ceiling
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                    # noqa: E402
from s13.qarch_lib import spearman                 # noqa: E402
from s14 import obj_enum as E                      # noqa: E402
from s14 import obj_model as M                     # noqa: E402
from s14 import obj_train as T                     # noqa: E402
from s14 import obj_floor as F                     # noqa: E402

BAND_MAX = 2.5
NBAND = 30000
LAMS = [1e-4, 1e-3, 1e-2, 1e-1]


def band_idx(en, nmax=NBAND, band_max=BAND_MAX, seed=0):
    sel = np.flatnonzero(en.rmsd <= band_max)
    if len(sel) > nmax:
        sel = np.random.default_rng(seed).choice(sel, nmax, replace=False)
    return np.sort(sel)


def pair_acc(e, r, rng, upto=1.5, npair=300_000):
    pi = F.pair_discrimination(e, r, rng, npair=npair)
    return float(np.nanmean([pi[k][0] for k in ("0.25-0.5", "0.5-1.0", "1.0-1.5")
                             if k in pi]))


def acc_of(fz, p, idx):
    """One target's centered normal-equation block."""
    en = T.Enum(p)
    a = M.RidgeAccumulator()
    a.add(fz.csr(p, T._decode(idx, en.n, en.k)),
          en.rmsd[idx].astype(np.float64), center=True)
    return a


def acc_sum(parts):
    """Normal equations are ADDITIVE across per-target centered blocks, so the
    leave-one-target-out fit is `total - that target` -- no refitting at all."""
    out = M.RidgeAccumulator()
    for a in parts:
        out.XtX += a.XtX
        out.Xty += a.Xty
        out.n += a.n
    return out


def acc_minus(total, part):
    out = M.RidgeAccumulator()
    out.XtX = total.XtX - part.XtX
    out.Xty = total.Xty - part.Xty
    out.n = total.n - part.n
    return out


def run(targets, seed=0):
    fz = M.Featurizer(cache_dir=True)
    bands = {p: band_idx(T.Enum(p), seed=seed) for p in targets}
    for p in targets:
        print(f"  {p}: band(RMSD<={BAND_MAX}) = {len(bands[p])} configs "
              f"of {T.Enum(p).B}", flush=True)
    rows = []
    t0 = time.time()
    halves = {}
    # MEMORY: each accumulator is a dense DIM x DIM float64 matrix (~136 MB at DIM=4125).
    # Keeping one per target would be >3 GB on a box with ~5 GB free, so only the running
    # TOTAL is retained and each target's own block is recomputed transiently below.
    TOTAL = M.RidgeAccumulator()
    for p in targets:
        idx = bands[p]
        h = np.random.default_rng(seed).permutation(len(idx))
        halves[p] = (idx[h[:len(idx) // 2]], idx[h[len(idx) // 2:]])
        a = acc_of(fz, p, idx)
        TOTAL.XtX += a.XtX
        TOTAL.Xty += a.Xty
        TOTAL.n += a.n
        del a
        print(f"  accumulated {p} [{time.time()-t0:.0f}s]", flush=True)

    print(f"\n{'pdb':6s} {'nband':>6s} "
          f"{'IN-SAMPLE':>22s} {'SAME-TARGET HELD OUT':>22s} {'CROSS-TARGET':>22s}")
    print(f"{'':6s} {'':6s} " + " ".join(f"{'pairacc':>9s} {'rho':>6s} {'d_t100':>6s}"
                                         for _ in range(3)))
    for p in targets:
        en = T.Enum(p)
        idx = bands[p]
        a, b = halves[p]
        rb = en.rmsd[b].astype(np.float64)
        rall = en.rmsd[idx].astype(np.float64)
        out = dict(pdb=p, fold=en.fold, nband=int(len(idx)),
                   band_mean=float(rall.mean()), band_min=float(rall.min()))

        # -- lambda chosen per arm to be MAXIMALLY GENEROUS to the arm --------
        def best(acc, sc_idx, r_sc, tag):
            bb = None
            S = T._decode(sc_idx, en.n, en.k)
            for lam in LAMS:
                w = acc.solve(lam, standardize=True)
                e = fz.score_with(w, p, S)
                pa = pair_acc(e, r_sc, np.random.default_rng(seed))
                if bb is None or pa > bb[0]:
                    o = np.argsort(e)
                    bb = (pa, spearman(e, r_sc),
                          float(r_sc[o[:100]].mean() - r_sc.mean()), lam)
            out[tag + "_pairacc"], out[tag + "_rho"] = bb[0], bb[1]
            out[tag + "_dt100"], out[tag + "_lam"] = bb[2], bb[3]
            return bb

        Af = acc_of(fz, p, idx)
        s1 = best(Af, idx, rall, "insample")
        Ah = acc_of(fz, p, a)
        s2 = best(Ah, b, rb, "heldout")
        del Ah
        Ax = acc_minus(TOTAL, Af)
        del Af
        s3 = best(Ax, idx, rall, "cross")
        del Ax
        rows.append(out)
        print(f"{p:6s} {len(idx):6d} "
              f"{s1[0]:9.3f} {s1[1]:+6.3f} {s1[2]:+5.2f} "
              f"{s2[0]:9.3f} {s2[1]:+6.3f} {s2[2]:+5.2f} "
              f"{s3[0]:9.3f} {s3[1]:+6.3f} {s3[2]:+5.2f}   [{time.time()-t0:.0f}s]",
              flush=True)

    print("\n-- MEANS over targets --")
    for tag, nm in (("insample", "in-sample (upper bound)"),
                    ("heldout", "same-target held out"),
                    ("cross", "cross-target")):
        pa = np.array([r[tag + "_pairacc"] for r in rows])
        rh = np.array([r[tag + "_rho"] for r in rows])
        dt = np.array([r[tag + "_dt100"] for r in rows])
        pr = I.paired(pa, np.full(len(pa), 0.5),
                      folds=np.array([r["fold"] for r in rows]),
                      names=[r["pdb"] for r in rows])
        print(f"  {nm:24s} pairacc<1.5A = {pa.mean():.3f} "
              f"CI[{0.5+pr['ci95'][0]:.3f},{0.5+pr['ci95'][1]:.3f}] "
              f"rho = {rh.mean():+.3f}  d_top100 = {dt.mean():+.3f} A  "
              f"(W/L vs chance {int((pa>0.5).sum())}/{int((pa<0.5).sum())})")
    print(f"\n  band mean RMSD {np.mean([r['band_mean'] for r in rows]):.3f} A, "
          f"band min {np.mean([r['band_min'] for r in rows]):.3f} A")
    I.write("s14_obj_ceiling", {"band_max": BAND_MAX, "rows": rows},
            n_expected=len(targets))
    return rows


if __name__ == "__main__":
    tg = [a for a in sys.argv[1:] if not a.startswith("-")]
    run(tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])][:12])
