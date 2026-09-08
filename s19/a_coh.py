"""SPRINT 19, AGENT A, Q2 -- IS THE ERROR COHERENT?

THE CRUX.  Sprint 18 established that permuting the distogram's residuals across pairs takes
the refinement from 3.610 A to 2.560-2.635 A.  Two incompatible mechanisms explain it:

  A3  COHERENCE.  The residual field is structurally coherent -- `dhat` is nearly the distance
      matrix of some OTHER, wrong structure -- so the fit converges confidently onto that wrong
      structure.  An incoherent field of the same size has no consistent structure to converge
      to, the least-squares compromise averages the contradictions out, and the answer stays
      nearer the truth.

  A4  UNREALISABILITY.  The field is so far from any Euclidean distance matrix that the fit is
      dragged onto an impossible compromise.

They pull in OPPOSITE directions: A3 says the field is too realisable, A4 says it is not
realisable enough.  `a_topo` already showed the predicted field carries 0.286 rank-3 EDM defect
against the native's 0.002 -- but that is not a control, because ANY error inflates the defect.
The control is a MAGNITUDE-MATCHED incoherent field.

ARMS (every one is an ORACLE DIAGNOSTIC: the target is built from `dtrue`).

  real         dhat, debiased, clipped -- reproduces objceil `a0.0` = 3.610
  signflip     dtrue + r * eps, eps = +-1 per pair.  THE WHITENING TEST.  Each pair keeps its
               OWN residual magnitude and its OWN weight (no G5/G6a orphaning) and the spatial
               pattern of |r| is untouched.  ONLY the cross-pair sign coherence is destroyed.
  shuffled     objceil's control, reproduced -- residual magnitudes permuted across pairs.
  coherent     THE POSITIVE CONTROL.  Replace the residual with the residual of a REAL
               alternative structure drawn from the retrieval pool, choosing the pool member
               whose residual RMS is closest to the real one.  Perfectly coherent (it IS a
               realisable matrix), matched magnitude.  If coherence is the harm this arm is as
               bad as `real` or worse, and it lands at the pool member's own RMSD.
  coh_scaled   the same pool member's residual, rescaled to match the real residual RMS exactly.
  iso          isotropic noise at the real residual RMS -- objceil's control, reproduced.

STATISTICS ALSO REPORTED (no fit needed):
  kappa = 1 - E_fit / E_true, the fraction of the predicted deviation that IS realisable by
  some ideal-geometry structure, where E_true = sum w (dhat - dtrue)^2 is the objective the
  NATIVE achieves and E_fit is the fit's terminal objective.  kappa near 1 means the field is
  coherently pointing somewhere; kappa near 0 means it is self-contradictory.

  EDM defect of the real, signflip and shuffled fields -- the A3/A4 crux without a fit.

Run:  python -m s19.a_coh
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s19 import a_lib as L                   # noqa: E402
from s19 import a_fit as F                   # noqa: E402
from s19.a_topo import edm_stats, spearman   # noqa: E402
from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_coh.json")
ARMS = ["real", "signflip", "shuffled", "coherent", "coh_scaled", "iso"]


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    data, pdbs, folds = L.gather_all(tg)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        d = data[pdb]
        i, j, sd, nat, n = d["i"], d["j"], d["sd"], d["nat"], d["n"]
        dhat, dtrue = d["dhat"], d["dtrue"]
        r = dhat - dtrue
        rms = float(np.sqrt((r ** 2).mean()))
        w = 1.0 / sd ** 2
        rng = SD.stable_rng(pdb, "s19A_coh")
        phi0, psi0, avg = F.start(pdb, d["seq"], d["fold"])

        # --- the coherent positive control: a REAL alternative structure from the pool
        W = np.asarray(AV.top75_windows(pdb)[0], float)
        DW = np.linalg.norm(W[:, i, :] - W[:, j, :], axis=-1)          # (75, npairs)
        rw = DW - dtrue[None]
        rms_w = np.sqrt((rw ** 2).mean(1))
        k = int(np.argmin(np.abs(rms_w - rms)))
        d_coh = DW[k]
        d_coh_s = dtrue + rw[k] * (rms / max(rms_w[k], 1e-9))

        fields = {
            "real": dhat,
            "signflip": dtrue + r * rng.choice([-1.0, 1.0], size=len(r)),
            "shuffled": dtrue + rng.permutation(r) * rng.choice([-1.0, 1.0], size=len(r)),
            "coherent": d_coh,
            "coh_scaled": d_coh_s,
            "iso": dtrue + rng.standard_normal(len(r)) * rms,
        }
        e = {"pdb": pdb, "n": n, "fold": d["fold"],
             "avg": float(I.ca_rmsd(avg, nat)),
             "resid_rms": rms,
             "pool_member_rmsd": float(I.ca_rmsd(W[k], nat)),
             "pool_member_rms": float(rms_w[k]),
             "E_true": float((w * r ** 2).sum())}
        for a in ARMS:
            fld = np.maximum(fields[a], 2.0)
            rm, fv = F.fit_rmsd(fld, sd, i, j, phi0, psi0, nat)
            e[a] = rm
            e[a + "_f"] = fv
            # the fraction of the field's deviation that a real structure can realise
            Et = float((w * (fld - dtrue) ** 2).sum())
            e[a + "_kappa"] = 1.0 - fv / max(Et, 1e-12)
            e[a + "_rms"] = float(np.sqrt(((fld - dtrue) ** 2).mean()))
            tri, dfc, neg = edm_stats(None, i, j, fld, n)
            e[a + "_defect"] = dfc
            e[a + "_tri"] = tri
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_coh.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "coh", "report")
    g = lambda k: np.array([r[k] for r in rows])                      # noqa: E731
    real = g("real")
    avg = g("avg")
    print(f"\n=== Q2  COHERENCE   n = {len(rows)}   ORACLE DIAGNOSTICS, all arms ===")
    print(f"  coordinate average (start)  {avg.mean():.3f}")
    print(f"  reproduction gate: real = {real.mean():.3f}  (objceil a0.0 = 3.610)")
    print(f"\n  {'arm':<13}{'RMSD':>8}{'med':>8}{'resid RMS':>11}{'kappa':>8}"
          f"{'EDMdef':>8}{'tri%':>7}   vs real")
    tab = {}
    for a in ARMS:
        v = g(a)
        m, lo, hi = L.boot(v - real, rng)
        d = v - real
        print(f"  {a:<13}{v.mean():>8.3f}{np.median(v):>8.3f}{g(a+'_rms').mean():>11.3f}"
              f"{g(a+'_kappa').mean():>8.3f}{g(a+'_defect').mean():>8.3f}"
              f"{100*g(a+'_tri').mean():>7.2f}   {m:+.3f} [{lo:+.3f},{hi:+.3f}]"
              f"  {int((d<0).sum())}W/{int((d>0).sum())}L")
        tab[a] = {"rmsd": float(v.mean()), "median": float(np.median(v)),
                  "resid_rms": float(g(a + "_rms").mean()),
                  "kappa": float(g(a + "_kappa").mean()),
                  "defect": float(g(a + "_defect").mean()),
                  "tri": float(g(a + "_tri").mean()),
                  "vs_real": {"diff": m, "ci": [lo, hi],
                              "W": int((d < 0).sum()), "L": int((d > 0).sum())}}

    print(f"\n  pool member actually used by `coherent`: mean RMSD "
          f"{g('pool_member_rmsd').mean():.3f}, residual RMS {g('pool_member_rms').mean():.3f}"
          f" (real {g('resid_rms').mean():.3f})")

    # THE WHITENING VERDICT
    sf = g("signflip")
    sh = g("shuffled")
    gap_total = float((real - sh).mean())
    gap_sign = float((real - sf).mean())
    print("\n  --- THE WHITENING VERDICT (pre-registered A3b) ---")
    print(f"  real - shuffled  (the whole Sprint-18 gap)     {gap_total:+.3f} A")
    print(f"  real - signflip  (sign coherence ONLY)         {gap_sign:+.3f} A")
    print(f"  fraction of the gap explained by sign coherence alone: "
          f"{gap_sign / gap_total if gap_total else float('nan'):.1%}")
    m, lo, hi = L.boot(sf - sh, rng)
    print(f"  signflip - shuffled  {m:+.3f} [{lo:+.3f},{hi:+.3f}]  "
          "(zero => sign-flipping alone reproduces the whole control)")

    # kappa and failure
    print("\n  --- kappa: how much of the predicted deviation is REALISABLE ---")
    for a in ARMS:
        print(f"    kappa({a:<11}) = {g(a+'_kappa').mean():.3f}")
    print(f"\n  rho(kappa_real, real RMSD)      {spearman(g('real_kappa'), real):+.3f}")
    print(f"  rho(EDM defect real, real RMSD) {spearman(g('real_defect'), real):+.3f}")
    print(f"  rho(resid_rms, real RMSD)       {spearman(g('resid_rms'), real):+.3f}")

    out = {"n": len(rows), "arms": tab,
           "gap_total": gap_total, "gap_sign": gap_sign,
           "signflip_minus_shuffled": {"diff": m, "ci": [lo, hi]},
           "rho_kappa_rmsd": spearman(g("real_kappa"), real),
           "rho_defect_rmsd": spearman(g("real_defect"), real)}
    json.dump(out, open(os.path.join(L.RESULTS, "a_coh_report.json"), "w"), indent=1,
              default=float)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
