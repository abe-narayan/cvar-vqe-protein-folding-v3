"""s24/c_probe2.py -- LANE C EXPLORATORY PROBE 2.  NOT A CLAIM.  Sizes the mixture effect.

c_probe.py established the (quality, bias) coordinates of four zero-training torsion samplers.
L2's criterion `q < 1/c_eff` is a MEAN-level statement that the library arm satisfies and still
loses, because the targets that violate it violate it badly.  So the criterion is not the
decision -- the endpoint is.  This runs the matched-size mixture ladder directly:

    75 members always.  a from the incumbent's scored top-75, 75-a from the generator's own
    scored top-75.  Identical readout, identical count.  Only composition changes.

Exploratory sizing for the pre-registered n=126 run.  Everything ORACLE-scored post hoc.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from core import project as pj           # noqa: E402
from s24 import c_probe as CP            # noqa: E402

TOPM = 75
NSAMP = 2000
MIX = [(75, 0), (60, 15), (50, 25), (38, 37), (25, 50), (0, 75)]


def run(pdbs=None, nsamp=NSAMP, arms=("T1_blind", "T2_restype", "T3_pool")):
    tg = I.targets()
    sel = tg if pdbs is None else [t for t in tg if t["pdb"] in set(pdbs)]
    rows = []
    for t in sel:
        pdb = t["pdb"]; n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); ij = I.pair_index(n)
        idx = I.pool_idx(u); Wp = u["W"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, *ij)), float)
        top = np.argsort(scp, kind="stable")[:TOPM]
        A = Wp[top]
        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]),
             "incumbent": float(I.ca_rmsd(CP._avg(A), nat))}
        rng0 = SD.stable_rng("c_probe2", pdb)
        for tag in arms:
            fn = {"T0_helix": lambda g: CP.s_helix(n, nsamp, g),
                  "T1_blind": lambda g: CP.s_blind(u, n, nsamp, g),
                  "T2_restype": lambda g: CP.s_restype(u, n, nsamp, g),
                  "T3_pool": lambda g: CP.s_pool(u, n, nsamp, g, idx[top])}[tag]
            ph, ps = fn(rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            sc = np.asarray(I.shipped_score(dg, I.pair_dists(CA, *ij)), float)
            B = CA[np.argsort(sc, kind="stable")[:TOPM]]
            mx = {}
            for a, b in MIX:
                s = np.concatenate([A[:a], B[:b]], 0) if (a and b) else (A[:a] if a else B[:b])
                mx["m%d_%d" % (a, b)] = float(I.ca_rmsd(CP._avg(s), nat))
            # the SS17 union: pool 500 incumbent windows with 500 generated, ONE score, top-75
            gsub = CA[SD.stable_rng("c_probe2", pdb, tag).permutation(len(CA))[:500]]
            Wu = np.concatenate([Wp, gsub], 0)
            src = np.concatenate([np.zeros(len(Wp), int), np.ones(len(gsub), int)])
            scu = np.asarray(I.shipped_score(dg, I.pair_dists(Wu, *ij)), float)
            ou = np.argsort(scu, kind="stable")[:TOPM]
            r[tag] = {"mix": mx, "union": float(I.ca_rmsd(CP._avg(Wu[ou]), nat)),
                      "union_gen_frac": float(src[ou].mean())}
        rows.append(r)
        print("  %s inc %.3f | %s" % (pdb, r["incumbent"], "  ".join(
            "%s u%.3f f%.2f" % (k[:2], r[k]["union"], r[k]["union_gen_frac"]) for k in arms)),
            flush=True)
    p = os.path.join(RES, "c_probe2.json"); tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "nsamp": nsamp, "EXPLORATORY": True}, fh)
    os.replace(tmp, p)
    rep(rows, arms)
    return rows


def rep(rows, arms=("T1_blind", "T2_restype", "T3_pool")):
    inc = np.array([r["incumbent"] for r in rows], float)
    print("\nEXPLORATORY n=%d targets.  Matched-size mixtures, 75 members always, POINT CLOUD."
          % len(rows))
    print("  incumbent %.4f" % inc.mean())
    for tag in arms:
        M = {k: np.array([r[tag]["mix"][k] for r in rows], float) for k in rows[0][tag]["mix"]}
        keys = list(M)
        e0, e1 = M[keys[0]], M[keys[-1]]
        print("  --- %s" % tag)
        for k in keys:
            frac = float(k.split("_")[1]) / TOPM
            line = (1 - frac) * e0 + frac * e1
            d = M[k] - e0
            se = d.std(ddof=1) / np.sqrt(len(d))
            print("      %-10s %8.4f  vs75/0 %+8.4f  vs line %+8.4f  SE %.4f MDE %.4f %3dW/%3dL"
                  % (k, M[k].mean(), d.mean(), (M[k] - line).mean(), se, 2.8016 * se,
                     int((d < 0).sum()), int((d > 0).sum())))
        un = np.array([r[tag]["union"] for r in rows], float)
        gf = np.array([r[tag]["union_gen_frac"] for r in rows], float)
        d = un - inc; se = d.std(ddof=1) / np.sqrt(len(d))
        print("      UNION      %8.4f  vs inc %+8.4f  SE %.4f MDE %.4f %3dW/%3dL  genfrac %.3f"
              % (un.mean(), d.mean(), se, 2.8016 * se, int((d < 0).sum()), int((d > 0).sum()),
                 gf.mean()))


if __name__ == "__main__":
    tg = I.targets()
    run([t["pdb"] for t in tg[::5]])
