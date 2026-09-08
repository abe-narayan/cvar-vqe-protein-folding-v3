"""SPRINT 16, AUDIT -- the two controls `s16/audit_align2.py` was missing.

1. AN ORACLE LINE SEARCH IS ITSELF A LEVER.  `audit_align2` line-searches along the surrogate
   with the step size chosen on the TRUE RMSD, so `c = 0` is always available and the arm can
   never lose.  The honest null is the same ORACLE line search along a RANDOM direction of the
   same length.  Without it, "the surrogate is worth -1.093 A" is unreadable.

2. STEERING NEEDS THE LOUD COMPONENT.  Cancelling error that already sits in the quiet half buys
   nothing (Part D: RMSD moves 0.002 A per 20 deg of quiet motion).  What a steerer must point at
   is the error's component in the LOUD half of J.  Measure |cos| there, against the same
   matched nulls.

    python -m s16.audit_align3
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
    os.environ[_v] = "2"

import numpy.linalg as la                     # noqa: E402
from s12 import instrument as I               # noqa: E402
from s15 import align_lib as A                # noqa: E402
from s15 import distcal as C                  # noqa: E402
from s15 import distgeo as D                  # noqa: E402
from s15 import seed as SD                    # noqa: E402
from s15.info_regime import ALPHA             # noqa: E402
from s14 import retprior as RP                # noqa: E402
import peptide_db as pdb                      # noqa: E402

PATH = os.path.join(ROOT, "s16", "results", "audit_align3.json")
CS = np.linspace(-1.5, 1.5, 61)


def main(limit=None):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)
    deb = {}
    for f in sorted({data[p]["fold"] for p in pdbs}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    for c_, t in enumerate(tg):
        p = t["pdb"]
        d = data[p]
        n = d["n"]; m = 2 * n
        i, j, sd = d["i"], d["j"], d["sd"]
        dhat = np.maximum(d["dhat"] - deb[d["fold"]](d["sep"]), 2.0)
        nat = d["nat"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        w0 = 1.0 / (sd * sd)
        rng = SD.stable_rng(p, "s16audit3")
        S = D.starts(p, d["seq"], n, d["fold"], 6, SD.stable_rng(p))
        sols = []
        for phi0, psi0, tag in S:
            ph, ps, fv = A.fit(dhat, sd, i, j, phi0, psi0, wpair=w0)
            sols.append((fv, ph, ps))
        fv, phi, psi = min(sols, key=lambda z: z[0])
        th = np.concatenate([phi, psi]); thn = np.concatenate([nphi, npsi])
        e = A.wrap(th - thn); ne = la.norm(e)
        P7, S7, _ = RP.windows(p, "top75")
        cm = np.concatenate([RP.circ_mean(P7, None, 0), RP.circ_mean(S7, None, 0)])
        v = A.wrap(th - cm)

        def line(u):
            u = u / la.norm(u)
            return min(float(I.ca_rmsd(I.build_ca((th - q * ne * u)[:n],
                                                  (th - q * ne * u)[n:]), nat)) for q in CS)

        r0 = float(I.ca_rmsd(I.build_ca(phi, psi), nat))
        rr = float(np.mean([line(rng.normal(size=m)) for _ in range(10)]))

        Jf, _, _ = A.sup_jacobian(phi, psi)
        s_f, V_f = A.spectrum(Jf)
        k = max(1, int(round(0.5 * m)))
        P = V_f[:, :k]                                    # the LOUD half
        el = P.T @ e                                      # error in loud coordinates
        vl = P.T @ v
        cl = float(abs(el @ vl) / (la.norm(el) * la.norm(vl)))
        # matched sign-flip null in the SAME loud coordinates
        nn = float(np.mean([abs(el @ (vl * rng.choice([-1.0, 1.0], k)))
                            / (la.norm(el) * la.norm(vl)) for _ in range(200)]))
        rows.append({"pdb": p, "n": n, "fold": int(d["fold"]), "rmsd": r0,
                     "ORACLE_line_RANDOM_dir": rr,
                     "cos_loud": cl, "cos_loud_signflip_null": nn,
                     "cos_loud_iso_null": float(np.sqrt(2 / (np.pi * k))),
                     "err_frac_loud": float((el @ el) / (e @ e))})
        if (c_ + 1) % 20 == 0:
            print(f"  {c_+1}/{len(tg)}", flush=True)

    folds = np.asarray([r["fold"] for r in rows], int)
    base = np.asarray([r["rmsd"] for r in rows], float)
    rr = np.asarray([r["ORACLE_line_RANDOM_dir"] for r in rows], float)
    pr = I.paired(rr, base, folds=folds)
    print(f"\nORACLE line search along a RANDOM direction (the missing null):")
    print(f"  {rr.mean():.3f} A   vs the fit {pr['mean_diff']:+.3f} "
          f"[{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] W/L {pr['n_better']}/{pr['n_worse']}")
    d2 = json.load(open(os.path.join(ROOT, "s16", "results", "audit_align2.json")))
    pool = np.asarray([r["ORACLE_line_pool_surrogate"] for r in d2["rows"]], float)
    hel = np.asarray([r["ORACLE_line_helix_ZEROINFO"] for r in d2["rows"]], float)
    p1 = I.paired(pool, rr, folds=folds)
    p2 = I.paired(hel, rr, folds=folds)
    print(f"  pool surrogate  MINUS random direction: {p1['mean_diff']:+.3f} "
          f"[{p1['ci95'][0]:+.3f},{p1['ci95'][1]:+.3f}] W/L {p1['n_better']}/{p1['n_worse']}")
    print(f"  alpha-helix ref MINUS random direction: {p2['mean_diff']:+.3f} "
          f"[{p2['ci95'][0]:+.3f},{p2['ci95'][1]:+.3f}] W/L {p2['n_better']}/{p2['n_worse']}")

    cl = np.asarray([r["cos_loud"] for r in rows], float)
    nl = np.asarray([r["cos_loud_signflip_null"] for r in rows], float)
    il = np.asarray([r["cos_loud_iso_null"] for r in rows], float)
    ef = np.asarray([r["err_frac_loud"] for r in rows], float)
    print(f"\nTHE COMPONENT A STEERER MUST POINT AT -- the error's LOUD half:")
    print(f"  share of ||e||^2 in the loud half: {ef.mean():.3f}")
    print(f"  |cos| of the surrogate with the error, INSIDE the loud half: {cl.mean():.3f} "
          f"(median {np.median(cl):.3f})")
    print(f"  its matched sign-flip null {nl.mean():.3f}; its isotropic null {il.mean():.3f}")
    dd = cl - nl
    print(f"  paired gap over the MATCHED null: {dd.mean():+.4f} W/L {(dd>0).sum()}/{(dd<0).sum()}")
    out = {"rows": rows, "line_random": {"mean": float(rr.mean()), "vs_fit": pr},
           "pool_minus_random": p1, "helix_minus_random": p2,
           "cos_loud": float(cl.mean()), "cos_loud_signflip_null": float(nl.mean()),
           "cos_loud_iso_null": float(il.mean()), "err_frac_loud": float(ef.mean())}
    with open(PATH, "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"\nwritten {PATH}")
    return out


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
