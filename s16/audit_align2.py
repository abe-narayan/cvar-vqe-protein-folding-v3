"""SPRINT 16, AUDIT -- follow-ups to `s16/audit_align.py`.

1. `theta_pool` IS `starts[0]`.  `s15/distgeo.starts` makes the retrieval top-75 circular mean
   the FIRST multi-start initialisation, and `s15/align_free.py` builds its headline surrogate
   `d_fit_pool = theta_fit - theta_pool` from that same circular mean.  So on every target where
   the objective-argmin start happens to be start 0, the "channel-disagreement" surrogate and the
   "where the fit moved" surrogate are the SAME VECTOR.  Measure how often.

2. A steering-relevant statement of the surrogate's skill: if you moved along
   `-c * (theta_fit - theta_pool)` for the BEST c (ORACLE line search), how much RMSD would you
   get?  This is not a proposal, it is the ceiling of the one-dimensional steering the surrogate
   supports, and it is the number the flagship needs.  Reported beside the same line search along
   the TRUE error direction (the ORACLE ceiling of a 1-D move) and along the zero-information
   alpha-helix reference direction.

    python -m s16.audit_align2
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

PATH = os.path.join(ROOT, "s16", "results", "audit_align2.json")
CS = np.concatenate([np.linspace(-1.5, 1.5, 61)])


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
    for c, t in enumerate(tg):
        p = t["pdb"]
        d = data[p]
        n = d["n"]
        i, j, sd = d["i"], d["j"], d["sd"]
        dhat = np.maximum(d["dhat"] - deb[d["fold"]](d["sep"]), 2.0)
        nat = d["nat"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        w0 = 1.0 / (sd * sd)
        S = D.starts(p, d["seq"], n, d["fold"], 6, SD.stable_rng(p))
        sols = []
        for phi0, psi0, tag in S:
            ph, ps, fv = A.fit(dhat, sd, i, j, phi0, psi0, wpair=w0)
            sols.append((fv, ph, ps, tag))
        fv, phi, psi, wtag = min(sols, key=lambda z: z[0])
        th = np.concatenate([phi, psi])
        thn = np.concatenate([nphi, npsi])
        e = A.wrap(th - thn)
        P7, S7, _ = RP.windows(p, "top75")
        cm = np.concatenate([RP.circ_mean(P7, None, 0), RP.circ_mean(S7, None, 0)])
        hel = np.concatenate([np.full(n, ALPHA[0]), np.full(n, ALPHA[1])])
        v_pool = A.wrap(th - cm)
        v_hel = A.wrap(th - hel)

        def line(v):
            nv = la.norm(v)
            if nv < 1e-12:
                return np.nan, np.nan
            u = v / nv
            best, bc = np.inf, 0.0
            for cq in CS:
                x = th - cq * la.norm(e) * u
                rr = float(I.ca_rmsd(I.build_ca(x[:n], x[n:]), nat))
                if rr < best:
                    best, bc = rr, float(cq)
            return best, bc

        r0 = float(I.ca_rmsd(I.build_ca(phi, psi), nat))
        rp, cp = line(v_pool)
        rh, ch = line(v_hel)
        re_, ce = line(e)
        rows.append({"pdb": p, "n": n, "fold": int(d["fold"]), "winner_start": wtag,
                     "rmsd": r0,
                     "ORACLE_line_pool_surrogate": rp, "c_pool": cp,
                     "ORACLE_line_helix_ZEROINFO": rh, "c_helix": ch,
                     "ORACLE_line_true_error": re_, "c_true": ce})
        if (c + 1) % 20 == 0:
            print(f"  {c+1}/{len(tg)}", flush=True)

    folds = np.asarray([r["fold"] for r in rows], int)
    base = np.asarray([r["rmsd"] for r in rows], float)
    tags = [r["winner_start"] for r in rows]
    from collections import Counter
    cnt = Counter(tags)
    print(f"\nWINNING MULTI-START (objective argmin), n = {len(rows)}")
    for k, v in cnt.most_common():
        print(f"  {k:<22}{v:>5}  ({v/len(rows)*100:.1f}%)")
    print(f"\n  On the {cnt.get('retrieval_circmean',0)} targets where the winner IS the "
          f"retrieval circular mean, `theta_fit - theta_pool` and `theta_fit - theta_start`\n"
          f"  are the IDENTICAL vector, so the two rows of s15/align_FINDINGS.md §5 are not "
          f"independent surrogates there.")

    print(f"\nORACLE LINE SEARCH along each direction (all ORACLE: c chosen on true RMSD)\n")
    print(f"{'direction moved along':<44}{'RMSD':>8}{'median':>8}{'  vs the fit':<28}")
    out = {"rows": rows, "winner_start_counts": dict(cnt)}
    for k, lab in (("rmsd", "the fit itself (no move)"),
                   ("ORACLE_line_pool_surrogate",
                    "ORACLE line search on theta_fit-theta_pool"),
                   ("ORACLE_line_helix_ZEROINFO",
                    "ORACLE line search on theta_fit-alpha-helix"),
                   ("ORACLE_line_true_error",
                    "ORACLE line search on the TRUE error (ceiling)")):
        v = np.asarray([r[k] for r in rows], float)
        pr = I.paired(v, base, folds=folds)
        out[k] = {"mean": float(v.mean()), "median": float(np.median(v)), "vs_fit": pr}
        print(f"{lab:<44}{v.mean():>8.3f}{np.median(v):>8.3f}"
              f"  {pr['mean_diff']:+.3f} [{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}]"
              f" W/L {pr['n_better']}/{pr['n_worse']}")
    a = np.asarray([r["ORACLE_line_pool_surrogate"] for r in rows], float)
    b = np.asarray([r["ORACLE_line_helix_ZEROINFO"] for r in rows], float)
    pr = I.paired(a, b, folds=folds)
    print(f"\n  pool surrogate MINUS the zero-information alpha-helix direction: "
          f"{pr['mean_diff']:+.3f} [{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}] "
          f"W/L {pr['n_better']}/{pr['n_worse']}")
    out["pool_minus_helix"] = pr
    cvals = np.asarray([r["c_pool"] for r in rows], float)
    print(f"  the ORACLE step size c along the pool surrogate: mean {cvals.mean():+.3f}, "
          f"median {np.median(cvals):+.3f}, sign consistency "
          f"{max((cvals>0).mean(), (cvals<0).mean())*100:.0f}% "
          f"(a native-free steerer must supply this sign)")
    out["c_pool_mean"] = float(cvals.mean())
    out["c_pool_sign_consistency"] = float(max((cvals > 0).mean(), (cvals < 0).mean()))
    with open(PATH, "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print(f"\nwritten {PATH}")
    return out


if __name__ == "__main__":
    main(limit=int(sys.argv[1]) if len(sys.argv) > 1 else None)
