"""SPRINT 19, AGENT A, Q2b -- WHICH coherent mode carries the harm?

`a_coh` establishes that destroying the residual's SIGN coherence, with every pair keeping
its own magnitude and its own weight, recovers the whole Sprint-18 gap.  So the harm is a
correlation structure.  This module names it.

PART 1 -- the correlation structure itself (ORACLE diagnostics, no fit).

  * corr(r_ij, r_ik) for pairs SHARING a residue, against pairs sharing none, against a
    within-target permutation null.
  * the variance of r decomposed into nested, interpretable modes, per target:
        offset      a single constant   r_ij ~ c
        scale       a single stretch    r_ij ~ (s - 1) * dtrue_ij
        seprofile   a per-separation profile  r_ij ~ b(|i-j|)   (5 bins)
        additive    a per-residue effect r_ij ~ a_i + a_j       (n parameters)
        residual    whatever is left

PART 2 -- surgery, through the deployed fit.  Each arm REMOVES one mode from the residual
and leaves everything else -- magnitudes, pair assignment, weights -- untouched.  Every arm
is an ORACLE DIAGNOSTIC: the mode is estimated against the native distances.  The question
each answers is "how much of the 1.0 A gap does THIS mode own?"

  real          dhat                                     3.610
  no_offset     r - mean(r)                              per-target constant removed
  no_scale      r - best multiplicative stretch          per-target size error removed
  no_sep        r - b(|i-j|) fitted per target           per-target separation profile removed
  no_additive   r - (a_i + a_j) fitted per target        per-residue misplacement removed
  only_additive keep ONLY (a_i + a_j)
  signflip      the reference whitening from a_coh
  shuffled      the Sprint-18 control

Every removal SHRINKS the residual, so each arm reports its own residual RMS and a
MAGNITUDE-MATCHED twin (`*_m`) that rescales the surviving residual back to the real
residual's RMS.  Without that twin an arm could win purely by having less error, which is
the trap `s18` G2d/G2e was written about.

Run:  python -m s19.a_struct
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
from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_struct.json")
MODES = ["offset", "scale", "seprofile", "additive"]
ARMS = (["real", "signflip", "shuffled", "only_additive"]
        + [f"no_{m}" for m in MODES] + [f"no_{m}_m" for m in MODES])


def _design(i, j, sep, dtrue, n, mode):
    """Design matrix for one coherent mode, weighted-LS fitted against r."""
    p = len(i)
    if mode == "offset":
        return np.ones((p, 1))
    if mode == "scale":
        return dtrue[:, None]
    if mode == "seprofile":
        b = np.clip(np.searchsorted([2, 4, 6, 8, 11, 999], sep, side="right") - 1, 0, 4)
        M = np.zeros((p, 5))
        M[np.arange(p), b] = 1.0
        return M
    if mode == "additive":
        M = np.zeros((p, n))
        M[np.arange(p), i.astype(int)] += 1.0
        M[np.arange(p), j.astype(int)] += 1.0
        return M
    raise ValueError(mode)


def _wls(M, r, w):
    A = M * w[:, None]
    coef, *_ = np.linalg.lstsq(M.T @ A, M.T @ (w * r), rcond=None)
    return M @ coef


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    data, pdbs, folds = L.gather_all(tg)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        d = data[pdb]
        i, j, sep, sd, nat, n = d["i"], d["j"], d["sep"], d["sd"], d["nat"], d["n"]
        dhat, dtrue = d["dhat"], d["dtrue"]
        r = dhat - dtrue
        w = 1.0 / sd ** 2
        rms = float(np.sqrt((r ** 2).mean()))
        rng = SD.stable_rng(pdb, "s19A_struct")
        phi0, psi0, avg = F.start(pdb, d["seq"], d["fold"])

        e = {"pdb": pdb, "n": n, "fold": d["fold"], "npairs": len(r),
             "avg": float(I.ca_rmsd(avg, nat)), "resid_rms": rms}

        # ---------------- PART 1: variance decomposition (weighted, nested)
        var_tot = float((w * r ** 2).sum() / w.sum())
        fitted = {}
        for m in MODES:
            M = _design(i, j, sep, dtrue, n, m)
            fitted[m] = _wls(M, r, w)
            e[f"var_{m}"] = float((w * fitted[m] ** 2).sum() / w.sum() / max(var_tot, 1e-12))
            e[f"k_{m}"] = int(M.shape[1])
        e["var_tot"] = var_tot

        # ---------------- PART 1b: sharing-a-residue correlation
        # correlation of standardised residuals over pair-pairs that share a residue
        z = (r - r.mean()) / max(r.std(), 1e-9)
        P = len(r)
        share = (i[:, None] == i[None, :]) | (i[:, None] == j[None, :]) | \
                (j[:, None] == i[None, :]) | (j[:, None] == j[None, :])
        tri = np.triu(np.ones((P, P), bool), 1)
        prod = z[:, None] * z[None, :]
        e["corr_share"] = float(prod[share & tri].mean()) if (share & tri).any() else 0.0
        e["corr_nonshare"] = float(prod[(~share) & tri].mean()) if ((~share) & tri).any() else 0.0
        zp = rng.permutation(z)
        pp = zp[:, None] * zp[None, :]
        e["corr_share_null"] = float(pp[share & tri].mean()) if (share & tri).any() else 0.0

        # ---------------- PART 2: surgery through the fit
        fields = {"real": dhat,
                  "signflip": dtrue + r * rng.choice([-1.0, 1.0], size=P),
                  "shuffled": dtrue + rng.permutation(r) * rng.choice([-1.0, 1.0], size=P),
                  "only_additive": dtrue + fitted["additive"]}
        for m in MODES:
            rr = r - fitted[m]
            fields[f"no_{m}"] = dtrue + rr
            s = rms / max(float(np.sqrt((rr ** 2).mean())), 1e-9)
            fields[f"no_{m}_m"] = dtrue + rr * s
        for a in ARMS:
            fld = np.maximum(fields[a], 2.0)
            rm, fv = F.fit_rmsd(fld, sd, i, j, phi0, psi0, nat)
            e[a] = rm
            e[a + "_rms"] = float(np.sqrt(((fld - dtrue) ** 2).mean()))
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_struct.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "struct", "report")
    g = lambda k: np.array([r[k] for r in rows])                     # noqa: E731
    real = g("real")
    print(f"\n=== Q2b  WHICH COHERENT MODE?   n = {len(rows)} ===")
    print(f"  reproduction gate: real = {real.mean():.3f}  (objceil a0.0 = 3.610)")

    print("\n  --- PART 1: weighted variance of the residual explained by each mode ---")
    print(f"  {'mode':<14}{'params':>8}{'share of var':>15}")
    for m in MODES:
        print(f"  {m:<14}{g('k_'+m).mean():>8.1f}{g('var_'+m).mean():>15.3f}")
    print(f"\n  correlation of standardised residuals, pairs SHARING a residue: "
          f"{g('corr_share').mean():+.3f}")
    print(f"                                        pairs sharing NONE:        "
          f"{g('corr_nonshare').mean():+.3f}")
    print(f"                                        within-target permutation null: "
          f"{g('corr_share_null').mean():+.3f}")

    print("\n  --- PART 2: surgery through the deployed fit (ALL ORACLE) ---")
    print(f"  {'arm':<16}{'RMSD':>8}{'med':>8}{'residRMS':>10}   vs real")
    tab = {}
    for a in ARMS:
        v = g(a)
        d = v - real
        m, lo, hi = L.boot(d, rng)
        print(f"  {a:<16}{v.mean():>8.3f}{np.median(v):>8.3f}{g(a+'_rms').mean():>10.3f}"
              f"   {m:+.3f} [{lo:+.3f},{hi:+.3f}]  {int((d<0).sum())}W/{int((d>0).sum())}L")
        tab[a] = {"rmsd": float(v.mean()), "median": float(np.median(v)),
                  "resid_rms": float(g(a + "_rms").mean()),
                  "vs_real": {"diff": m, "ci": [lo, hi],
                              "W": int((d < 0).sum()), "L": int((d > 0).sum())}}
    gap = float((real - g("shuffled")).mean())
    print(f"\n  the Sprint-18 gap (real - shuffled) = {gap:+.3f} A.  Share owned by each mode,")
    print("  MAGNITUDE-MATCHED (so the credit is not just 'less error'):")
    for m in MODES:
        own = float((real - g(f"no_{m}_m")).mean())
        print(f"    {m:<14}{own:+.3f} A   ({own / gap if gap else float('nan'):>6.1%} of the gap)")
        tab[f"_own_{m}"] = own
    tab["_gap"] = gap
    json.dump(tab, open(os.path.join(L.RESULTS, "a_struct_report.json"), "w"), indent=1,
              default=float)
    return tab


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
