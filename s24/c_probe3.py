"""s24/c_probe3.py -- LANE C EXPLORATORY PROBE 3.  WHERE DOES THE COMMON BIAS COME FROM?

c_probe/c_probe2 found something that needs explaining.  A constant alpha-helix -- a source
carrying ZERO information about the target -- emits, after top-75 distogram scoring, a cloud whose
error is 0.72 cosine-aligned with the incumbent's.  Four samplers spanning "knows nothing" to
"knows the retrieved pool" all land between 0.72 and 0.91.  A bias that survives replacing the
entire source is not a property of the source.

HYPOTHESIS.  The common mode is the DISTOGRAM'S OWN PREDICTION ERROR.  Every arm is selected by
the same Bayes-risk functional over the same predicted distance matrix, so every arm is pulled
toward the same wrong geometry.  s23 L9's 68% common-mode fraction, L2's 0.647 cosine and L2(d)'s
precise null would then all be one quantity.

THE TEST, on signed PAIR-DISTANCE errors, which is the space the score actually acts in:

    g_p = dhat_p - d_p(native)          the distogram's signed error on pair p   (native-free input,
                                        ORACLE-scored here)
    e_p = d_p(emitted) - d_p(native)    the emitted cloud's signed error on pair p

If cos(e, g) is large and roughly EQUAL across sources that share nothing but the score, the score
is the bias.  Controls that make it mean something:
  - the UNSCORED uniform-75 average of the same emitted set (removes the selector, keeps the source)
  - the emitted set's own MEAN pair distance error (removes selection entirely)
  - a within-source ceiling, as in biasalign.

Exploratory.  Not a claim.  ORACLE labelling on every quantity that reads the native.
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


def _cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


def run(pdbs=None, nsamp=NSAMP):
    tg = I.targets()
    sel = tg if pdbs is None else [t for t in tg if t["pdb"] in set(pdbs)]
    arms = ("T0_helix", "T1_blind", "T2_restype", "T3_pool")
    rows = []
    for t in sel:
        pdb = t["pdb"]; n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n)
        dnat = I.pair_dists(nat[None], i, j)[0]
        # the distogram's own signed error, on the SAME pair set the score uses
        dhat = np.asarray(dg["expected"], float)
        g = dhat - dnat

        idx = I.pool_idx(u); Wp = u["W"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, i, j)), float)
        top = np.argsort(scp, kind="stable")[:TOPM]
        cA = CP._avg(Wp[top])
        eA = I.pair_dists(cA[None], i, j)[0] - dnat
        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]),
             "incumbent_rmsd": float(I.ca_rmsd(cA, nat)),
             "cos_inc_dgram": _cos(eA, g),
             "dgram_mae": float(np.abs(g).mean()), "dgram_bias": float(g.mean()),
             "inc_dbias": float(eA.mean())}
        rng0 = SD.stable_rng("c_probe3", pdb)
        for tag in arms:
            fn = {"T0_helix": lambda q: CP.s_helix(n, nsamp, q),
                  "T1_blind": lambda q: CP.s_blind(u, n, nsamp, q),
                  "T2_restype": lambda q: CP.s_restype(u, n, nsamp, q),
                  "T3_pool": lambda q: CP.s_pool(u, n, nsamp, q, idx[top])}[tag]
            ph, ps = fn(rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            D = I.pair_dists(CA, i, j)
            sc = np.asarray(I.shipped_score(dg, D), float)
            o = np.argsort(sc, kind="stable")[:TOPM]
            cS = CP._avg(CA[o])                                   # scored
            ru = rng0.choice(len(CA), TOPM, replace=False)
            cU = CP._avg(CA[ru])                                  # unscored
            eS = I.pair_dists(cS[None], i, j)[0] - dnat
            eU = I.pair_dists(cU[None], i, j)[0] - dnat
            eSet = D.mean(0) - dnat                               # no selection at all
            r[tag] = {"cos_scored_dgram": _cos(eS, g), "cos_unif_dgram": _cos(eU, g),
                      "cos_set_dgram": _cos(eSet, g),
                      "cos_scored_inc": _cos(eS, eA), "cos_unif_inc": _cos(eU, eA),
                      "rmsd_scored": float(I.ca_rmsd(cS, nat)),
                      "rmsd_unif": float(I.ca_rmsd(cU, nat))}
        rows.append(r)
        print("  %s inc %.3f  cos(inc,dgram) %+0.3f | %s" % (
            pdb, r["incumbent_rmsd"], r["cos_inc_dgram"],
            " ".join("%s %+0.3f/%+0.3f" % (k[:2], r[k]["cos_scored_dgram"],
                                           r[k]["cos_set_dgram"]) for k in arms)), flush=True)
    p = os.path.join(RES, "c_probe3.json"); tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "nsamp": nsamp, "EXPLORATORY": True}, fh)
    os.replace(tmp, p)
    rep(rows, arms)
    return rows


def rep(rows, arms=("T0_helix", "T1_blind", "T2_restype", "T3_pool")):
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    print("\nEXPLORATORY n=%d.  Cosines are between SIGNED PAIR-DISTANCE error vectors." % len(rows))
    print("  distogram MAE %.3f A   signed mean error %+0.3f A   incumbent's own %+0.3f A"
          % (g("dgram_mae").mean(), g("dgram_bias").mean(), g("inc_dbias").mean()))
    print("\n  cos( emitted distance error , DISTOGRAM's own distance error )")
    print("  %-14s %10s %10s %10s | %10s %10s" % ("arm", "scored75", "unif75", "whole set",
                                                  "cos vs inc", "RMSD"))
    print("  %-14s %10.4f" % ("INCUMBENT", g("cos_inc_dgram").mean()))
    for k in arms:
        f = lambda q: np.array([r[k][q] for r in rows], float)    # noqa: E731
        print("  %-14s %10.4f %10.4f %10.4f | %10.4f %10.4f"
              % (k, f("cos_scored_dgram").mean(), f("cos_unif_dgram").mean(),
                 f("cos_set_dgram").mean(), f("cos_scored_inc").mean(), f("rmsd_scored").mean()))
    print("\n  If the score is the bias, cos(.,dgram) is HIGH under 'scored75' for every arm")
    print("  regardless of source, and DROPS when the selector is removed ('whole set').")


if __name__ == "__main__":
    tg = I.targets()
    run([t["pdb"] for t in tg[::5]])
