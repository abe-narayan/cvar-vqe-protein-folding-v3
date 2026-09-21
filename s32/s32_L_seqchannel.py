"""LANE L / L3c -- does the sequence->structure retrieval channel strengthen with length?

NOT PRE-REGISTERED.  Declared EXPLORATORY and logged as such: it was prompted by the
ladder's own `filter_skill` column, where BLOSUM's within-pool ordering is worth -2.19 A at
L ~ 56 against -0.15 A at L ~ 13.  An effect found after looking is a hypothesis, not a
result, and it is reported that way.

THE FINDING IT BEARS ON.  Project memory records `structure-and-sequence-are-decoupled`:
*the best-matching window has 12% identity; sequence is the wrong retrieval key.*  That was
measured at 9-16 residues.  A 13-residue BLOSUM match is weak evidence of anything; a
55-residue one is close to a homology statement.  So the decoupling may be length-scoped.

THE STATISTIC.  Within each target's K = 500 BLOSUM pool, the Spearman correlation between
the BLOSUM similarity (the deployed retrieval key, native-free) and the ORACLE CA-RMSD of
each window.  Negative is skill: higher similarity should mean lower RMSD.  This is an
IN-BAND measure -- it is computed inside the already-retrieved pool, so it cannot be
inflated by garbage rejection (memory: `decoy-bank-not-a-pool-proxy`), which is the trap a
whole-universe correlation would fall into.

Also reported: the identity of the single best-matching window (the quantity the memory
entry quotes as 12%), at both lengths, computed the same way.

    python -m s32.s32_L_seqchannel run
"""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
K = 500


def run(verbose=True):
    from scipy.stats import spearmanr
    from core import geometry as geo
    from s12 import instrument as I
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    out = {}
    for kind in ("short", "long"):
        rows = []
        if kind == "short":
            tg = I.targets()
        else:
            tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
            corpus = cp.load()
            by = {p: a for a, p in enumerate(corpus["pdb"])}
        for k, t in enumerate(tg):
            if kind == "short":
                u = I.load_univ(t["pdb"])
                p = I.pool_idx(u, K)
                W = np.asarray(u["W"], float)[p]
                sim = np.asarray(u["sim"], float)[p]
                nat = np.asarray(u["nat_ca"], float)
            else:
                W, sim, src = LD._bank_long(t)
                pp = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
                if not os.path.exists(pp):
                    pp = glob.glob(os.path.join(BASE, "prots",
                                                t["pdb"].lower() + ".pdb"))[0]
                _, coords, _, _ = geo.native_coords_from_pdb(pp)
                nat = np.asarray(coords["CA"], float)
            rr = I.kabsch_rmsd_batch(W, nat)                            # ORACLE labels
            rho = float(spearmanr(sim, rr).statistic)
            # rank of the ORACLE best member inside the BLOSUM order (0 = the top hit)
            rows.append(dict(pdb=t["pdb"], n=int(len(nat)), rho=rho,
                             best_rank=int(np.argmin(rr)),
                             rr_top1=float(rr[0]), rr_best=float(rr.min()),
                             rr_mean=float(rr.mean())))
            if verbose and (k + 1) % 25 == 0:
                print(f"  {kind} {k+1}/{len(tg)}", flush=True)
        rh = np.array([r["rho"] for r in rows])
        br = np.array([r["best_rank"] for r in rows], float)
        out[kind] = dict(
            n=len(rows), mean_len=float(np.mean([r["n"] for r in rows])),
            rho_mean=float(rh.mean()), rho_median=float(np.median(rh)),
            rho_se=float(rh.std(ddof=1) / len(rh) ** 0.5),
            frac_rho_negative=float((rh < 0).mean()),
            best_rank_median=float(np.median(br)),
            rr_top1_mean=float(np.mean([r["rr_top1"] for r in rows])),
            rr_best_mean=float(np.mean([r["rr_best"] for r in rows])),
            rr_mean_mean=float(np.mean([r["rr_mean"] for r in rows])),
            rows=rows)
    s, lg = out["short"], out["long"]
    sed = (s["rho_se"] ** 2 + lg["rho_se"] ** 2) ** 0.5
    eff = lg["rho_mean"] - s["rho_mean"]
    out["length_contrast"] = dict(
        effect=float(eff), se=float(sed), mde=float(2.8016 * sed),
        over_mde=float(abs(eff) / (2.8016 * sed)),
        band=("NOT A RESULT" if abs(eff) < 0.7 * 2.8016 * sed else
              "NOT MEASURED" if abs(eff) < 2.8016 * sed else "RESULT"),
        note="EXPLORATORY, not pre-registered; found by looking at the ladder's "
             "filter_skill column. Unpaired across two different instruments.")
    path = os.path.join(RESULTS, "L3c_seqchannel.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: ({a: b for a, b in v.items() if a != "rows"}
                              if isinstance(v, dict) else v)
                          for k, v in out.items()}, indent=1))
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
