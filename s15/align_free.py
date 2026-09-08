"""SPRINT 15, ALIGN -- IS THERE A NATIVE-FREE VERSION OF THE ALIGNMENT DIAGNOSTIC?

The brief's rule: for every ORACLE diagnostic, say whether a native-free version exists and
measure it if so.  `alignment(e) = ||J e|| / (||e|| s_rms)` has two halves and they have
different answers:

  the SUBSPACE (J and its spectrum)  -- NATIVE-FREE.  `s15/align_jac.py` measures the
      bottom-half subspace of J at the incumbent's emitted torsions against the same subspace
      at the native torsions: mean cos^2 0.827 against a 0.499 random-subspace null, spectra
      correlated +0.980.  The geometry does not need the answer.

  the ERROR DIRECTION (e) -- ORACLE by construction.  You cannot know which way you are wrong.
      This module asks whether any native-free surrogate for `e` carries the alignment signal:

        d_fit_start   theta_fit - theta_start            where the fit moved
        d_fit_pool    theta_fit - pool circular mean     disagreement with the retrieval channel
        d_pool_spread the pool's own principal direction of disagreement
        d_multistart  the spread across the multi-start solutions

      Each is scored two ways: (a) does its own alignment correlate with the TRUE alignment
      across the 126 targets, and (b) does it point the same WAY as the true error (cosine)?

A surrogate that fails (a) cannot even diagnose alignment; one that passes (a) but fails (b)
can rank targets but cannot steer a fit.  Both are reported.

    python -m s15.align_free
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import info_lib as L                # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s14 import retprior as RP               # noqa: E402
import peptide_db as pdb                     # noqa: E402


def main(limit=None, n_start=6):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)
    folds = sorted({data[p]["fold"] for p in pdbs})
    deb = {}
    for f in folds:
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rec = []
    for c, t in enumerate(tg):
        p = t["pdb"]
        d = data[p]
        n = d["n"]
        i, j, sd = d["i"], d["j"], d["sd"]
        dhat = np.maximum(d["dhat"] - deb[d["fold"]](d["sep"]), 2.0)
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        w0 = 1.0 / (sd * sd)
        S = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p))
        sols = []
        for phi0, psi0, _tag in S:
            ph, ps, fv = A.fit(dhat, sd, i, j, phi0, psi0, wpair=w0)
            sols.append((fv, ph, ps, np.concatenate([np.asarray(phi0, float),
                                                     np.asarray(psi0, float)])))
        fv, phi, psi, th0 = min(sols, key=lambda z: z[0])

        Jn, _, _ = A.sup_jacobian(nphi, npsi)                 # ORACLE J
        Js, _, _ = A.sup_jacobian(phi, psi)                   # NATIVE-FREE J
        e_true = np.concatenate([A.wrap(phi - nphi), A.wrap(psi - npsi)])

        P7, S7, _ = RP.windows(p, "top75")
        cm = np.concatenate([RP.circ_mean(P7, None, 0), RP.circ_mean(S7, None, 0)])
        th = np.concatenate([phi, psi])
        Wp = np.concatenate([P7, S7], axis=1)
        Wc = A.wrap(Wp - cm[None, :])
        _, _, Vt = np.linalg.svd(Wc - Wc.mean(0), full_matrices=False)
        spread = Vt[0]
        ms = np.asarray([np.concatenate([a, b]) for _, a, b, _s in sols])
        msd = ms.std(0)
        msd = msd / max(np.linalg.norm(msd), 1e-12)

        sur = {"d_fit_start": A.wrap(th - th0),
               "d_fit_pool": A.wrap(th - cm),
               "d_pool_spread": spread,
               "d_multistart": msd}
        r = {"pdb": p, "n": n, "fold": d["fold"],
             "rmsd": float(I.ca_rmsd(I.build_ca(phi, psi), nat := d["nat"])),
             "align_true_ORACLE": A.alignment(Jn, e_true),
             "align_true_selfJ": A.alignment(Js, e_true)}
        ne = float(np.linalg.norm(e_true))
        for k, v in sur.items():
            nv = float(np.linalg.norm(v))
            r["align_" + k] = A.alignment(Js, v) if nv > 1e-12 else float("nan")
            r["cos_" + k] = float(abs(v @ e_true) / (nv * ne)) if nv > 1e-12 else float("nan")
        rec.append(r)
        if (c + 1) % 20 == 0:
            print(f"  {c+1}/{len(tg)}", flush=True)

    out = {"rows": rec}
    at = np.asarray([r["align_true_ORACLE"] for r in rec], float)
    asj = np.asarray([r["align_true_selfJ"] for r in rec], float)
    print("\nIS THERE A NATIVE-FREE VERSION OF THE ALIGNMENT DIAGNOSTIC?  "
          f"n = {len(rec)} targets\n")
    print(f"the ORACLE quantity: alignment of the predicted fit's TRUE torsion error, "
          f"J at native  = {at.mean():.3f}")
    print(f"the same with J at the EMITTED structure (native-free J, ORACLE error) "
          f"= {asj.mean():.3f}   corr {L.pearson(at, asj):+.3f}")
    out["align_true_ORACLE"] = float(at.mean())
    out["align_true_selfJ"] = float(asj.mean())
    out["corr_J_native_vs_emitted"] = L.pearson(at, asj)
    print(f"\n{'native-free surrogate for the ERROR DIRECTION':<38}{'its align':>11}"
          f"{'corr w/ true align':>20}{'|cos| with true e':>20}{'random |cos|':>14}")
    m = 2 * float(np.mean([r["n"] for r in rec]))
    for k in ("d_fit_start", "d_fit_pool", "d_pool_spread", "d_multistart"):
        a = np.asarray([r["align_" + k] for r in rec], float)
        co = np.asarray([r["cos_" + k] for r in rec], float)
        ok = np.isfinite(a) & np.isfinite(co)
        rr = L.pearson(at[ok], a[ok])
        ci = L.boot_ci(at[ok], a[ok])
        out[k] = {"align_mean": float(a[ok].mean()), "corr_with_true_align": rr,
                  "ci95": ci, "cos_mean": float(co[ok].mean()), "n": int(ok.sum())}
        print(f"{k:<38}{a[ok].mean():>11.3f}{rr:>+14.3f} "
              f"[{ci[0]:+.2f},{ci[1]:+.2f}]{co[ok].mean():>20.3f}"
              f"{np.sqrt(2 / (np.pi * m)):>14.3f}")
    print(f"\n  the random-|cos| column is the expected |cos| between two independent "
          f"directions in\n  2n ~ {m:.0f} dimensions, sqrt(2/(pi*2n)) = "
          f"{np.sqrt(2 / (np.pi * m)):.3f}.")
    L.jwrite("align_free", out)
    return out


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
