"""SPRINT 19, AGENT A -- SELF-AUDIT OF THE HEADLINE (brief section 11).

THE DEFECT I FOUND IN MY OWN ARM.  Every field in `a_coh` is clipped at 2.0 A before the fit,
exactly as `s18/objceil.py` does.  Flipping a residual's sign can push a distance below 2.0,
so the clip SHRINKS the whitened field's magnitude: `signflip` carries residual RMS 3.062
against `real`'s 3.205, a 4.5% reduction.  On the alpha ladder a 25% magnitude cut is worth
0.504 A, so 4.5% extrapolates to about 0.091 A -- roughly 7% of the 1.242 A headline.  Small,
but the headline is the sprint's central result and 7% is not nothing.

THE FIX.  `signflip_exact` rescales the whitened residual AFTER the clip so that its residual
RMS equals `real`'s to machine precision, then re-clips and iterates to a fixed point.  If the
effect survives at EXACTLY matched magnitude, the clip is not carrying it.

Also carried: `shuffled_exact` and `iso_exact`, since those controls are clipped harder still
(3.008 and 3.006, a 6.2% reduction), which means the published `signflip - shuffled` gap of
-0.282 was measured with the magnitude confound pointing AGAINST signflip.

Run:  python -m s19.a_audit
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

OUT = os.path.join(L.RESULTS, "a_audit.json")
ARMS = ["real", "signflip", "signflip_exact", "shuffled", "shuffled_exact",
        "iso", "iso_exact"]


def match_rms(r, dtrue, target_rms, iters=60):
    """Scale `r` so that, AFTER clipping at 2.0, the realised residual RMS equals `target_rms`.

    The clip is a nonlinearity, so this is a fixed-point iteration rather than one division.
    Converges in a handful of steps; 60 is belt and braces.
    """
    s = 1.0
    for _ in range(iters):
        fld = np.maximum(dtrue + r * s, 2.0)
        cur = float(np.sqrt(((fld - dtrue) ** 2).mean()))
        if cur < 1e-9:
            break
        s *= target_rms / cur
        if abs(cur - target_rms) < 1e-10:
            break
    return np.maximum(dtrue + r * s, 2.0)


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    data, pdbs, folds = L.gather_all(tg)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], d["nat"]
        dhat, dtrue = d["dhat"], d["dtrue"]
        r = dhat - dtrue
        P = len(r)
        rng = SD.stable_rng(pdb, "s19A_coh")          # SAME seed stream as a_coh
        phi0, psi0, avg = F.start(pdb, d["seq"], d["fold"])
        real_field = np.maximum(dhat, 2.0)
        target = float(np.sqrt(((real_field - dtrue) ** 2).mean()))

        # reproduce a_coh's draws in the same order so the arms are the identical realisations
        sf_r = r * rng.choice([-1.0, 1.0], size=P)
        sh_r = rng.permutation(r) * rng.choice([-1.0, 1.0], size=P)

        e = {"pdb": pdb, "n": d["n"], "fold": d["fold"], "target_rms": target}
        iso_r = rng.standard_normal(P) * float(np.sqrt((r ** 2).mean()))
        fields = {
            "real": real_field,
            "signflip": np.maximum(dtrue + sf_r, 2.0),
            "signflip_exact": match_rms(sf_r, dtrue, target),
            "shuffled": np.maximum(dtrue + sh_r, 2.0),
            "shuffled_exact": match_rms(sh_r, dtrue, target),
            "iso": np.maximum(dtrue + iso_r, 2.0),
            "iso_exact": match_rms(iso_r, dtrue, target),
        }
        for a in ARMS:
            e[a] = F.fit_rmsd(fields[a], sd, i, j, phi0, psi0, nat)[0]
            e[a + "_rms"] = float(np.sqrt(((fields[a] - dtrue) ** 2).mean()))
        rows.append(e)
        if (c + 1) % 20 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_audit.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "audit", "report")
    g = lambda k: np.array([r[k] for r in rows], float)               # noqa: E731
    real = g("real")
    folds = g("fold")
    print(f"\n=== SELF-AUDIT OF THE HEADLINE   n = {len(rows)} ===")
    print(f"  {'arm':<18}{'RMSD':>8}{'med':>8}{'residRMS':>10}   vs real                     folds")
    tab = {}
    for a in ARMS:
        v = g(a)
        d = v - real
        m, lo, hi = L.boot(d, rng)
        sg = sum(1 for f in np.unique(folds)
                 if m != 0 and np.sign(d[folds == f].mean()) == np.sign(m))
        print(f"  {a:<18}{v.mean():>8.3f}{np.median(v):>8.3f}{g(a+'_rms').mean():>10.4f}"
              f"   {m:+.3f} [{lo:+.3f},{hi:+.3f}] {int((d<0).sum()):>3}W/{int((d>0).sum())}L  {sg}/5")
        tab[a] = {"rmsd": float(v.mean()), "resid_rms": float(g(a + "_rms").mean()),
                  "diff": m, "ci": [lo, hi], "W": int((d < 0).sum()),
                  "L": int((d > 0).sum()), "folds": sg}
    s, rec = L.report_pair("signflip_exact - shuffled_exact", g("signflip_exact"),
                           g("shuffled_exact"), rng, folds)
    print("\n  " + s)
    tab["_sf_vs_sh_exact"] = rec
    json.dump(tab, open(os.path.join(L.RESULTS, "a_audit_report.json"), "w"), indent=1,
              default=float)
    return tab


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
