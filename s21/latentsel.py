"""s21/latentsel.py -- CAN 512 RANDOM LATENT DRAWS REPLACE THE RETRIEVAL POOL?

THE QUESTION THIS SPRINT ARRIVED AT, AND WHY IT IS THIS ONE.

L14/L17 enumerated the entire latent at n=126: the objective's EXACT argmin ties a zero-evaluation
pool (+0.085, NOT MEASURED) while the ORACLE over the identical modes is 1.72-1.81 A better, with
ZERO reversals in 126 targets.  L18 then priced the readout route with WORKSTREAM D's matched
control and killed it: the objective's ordering skill is -1.184 A at M=1 and is EXHAUSTED BY M=8;
at M=512 its ordered set is no better than 512 RANDOM draws (+0.107, CI spanning zero).

    SO THE ONE THING THE MEASUREMENTS ACTUALLY ESTABLISH IS THIS:
    512 ARBITRARY DRAWS FROM THE LATENT CONTAIN A STRUCTURE AT 1.879 A --
    0.886 A BETTER THAN THE SHIPPED INCUMBENT -- AND THE OBJECTIVE CANNOT FIND IT.

That is an ORACLE statement.  It becomes an ARCHITECTURE only if some NATIVE-FREE operator can
extract it.  This file asks exactly that, and it asks it with the operator the pipeline ALREADY
SHIPS rather than with a new one invented for the occasion:

    THE POOL AND THE LATENT ARE TWO SOURCES OF CANDIDATE STRUCTURES.
    RUN THE IDENTICAL DEPLOYED SELECTION OPERATOR ON BOTH.  WHICH SOURCE IS BETTER?

The incumbent (`pool_tail_avg`) is *already* score-filter + top-75 coordinate average over the K=500
BLOSUM pool.  The matched latent arm is the same operator over 512 random latent draws.  Nothing new
is being proposed; the SOURCE is being swapped and everything downstream is held fixed.  If the
latent wins, the retrieval stage -- not the selector, not the ansatz, not the optimiser -- is what
four sprints of RMSD has been paying for.

ARMS, all native-free except the two labelled ORACLE.  512 draws, matched to the pool's evaluation
budget and to L18's control size.

    lat_rand1        one random draw                            the zero-information null
    lat_argmin       argmin of the deployed objective            the deployed READOUT
    lat_avg75        top-75 by objective -> coordinate average   THE DEPLOYED OPERATOR, on the latent
    lat_avg75_rand   75 RANDOM draws -> coordinate average       averaging WITHOUT the objective
    lat_medoid       consensus medoid of all 512                 the only in-band discriminator with
                                                                 measured positive skill (S12)
    lat_oracle       best of the 512 by true RMSD                ORACLE ceiling of the DRAW SET
    pool_avg75       THE SHIPPED INCUMBENT, identical operator   the thing to beat
    pool_oracle      best of the K=500 pool by true RMSD         ORACLE, for reference

WHY `lat_avg75_rand` IS PRESENT AND IS NOT DECORATION.  `averaging-space-beats-the-objective` records
that coordinate averaging is worth ~1.0 A while the objective contributes ~0.17 A, and L18 has just
shown the objective's ordering is worthless past its top handful.  If `lat_avg75` beats the pool it
is therefore ESSENTIAL to know whether the objective did any of that work.  Without this arm a win
would be attributed to the selector when it belongs to the averaging operator and the source.

THE CONTRACTION HAZARD, DECLARED BECAUSE IT HAS ALREADY COST THIS PROJECT A PUBLISHED NUMBER.
Coordinate averaging CONTRACTS the backbone (mean virtual bond 2.961 A against a physical 3.804;
the "best built 3.048 A" of Sprint 20 was a contracted point cloud, corrected in five places).  Every
averaged arm here is therefore a POINT CLOUD and is labelled so in the output.  A point-cloud RMSD is
NOT comparable to a built-chain RMSD, and the comparison that matters -- `lat_avg75` vs `pool_avg75`
-- is point-cloud against point-cloud, so it is internally valid; the single-structure arms are built
chains and are marked separately.  Do not read across the two groups.

OPERATOR FORKS, DECLARED PER BRIEF SS7 RULE 0 -- the rule the coordinator's own failure produced.

    functional     DECLARED shipped Bayes risk, the deployed selector.  NOT TAKEN squared-distance.
    basis          DECLARED built chains from torsions for BOTH sources -- the pool arm is rebuilt
                   from its windows' own torsions, exactly as `pool_argmin_rb` in L14.  NOT TAKEN
                   the pool's retrieved coordinates, which would give the pool a different manifold
                   from the latent and is the fork that broke L14's first primary.
    readout        DECLARED top-75 coordinate average for the primary, because it is what SHIPS.
                   NOT TAKEN argmin, reported beside it; NOT TAKEN a tail mean, which the pillar
                   does not have (Sprint 21 SS2).
    normalisation  none; RMSD in Angstroms, no per-target scaling.  NOT TAKEN z-scoring per target,
                   which would hide that the mean is set by which targets the arm FAILS on.
    null           TWO.  `lat_rand1` for "is any of this better than a coin flip", and
                   `lat_avg75_rand` for "did the OBJECTIVE do the work, or the AVERAGING".

PRE-REGISTRATION.
  Primary      `lat_avg75 - pool_avg75`, paired over 126 targets, fold-clustered, SE and the
               per-comparison MDE = 2.8016*SE beside it.
  Hypothesis   the latent source is COMPETITIVE but does not clearly win, because the draws are
               unconditioned on sequence beyond the basin fit, while the pool is retrieved.
  Falsifier    if `lat_avg75` beats `pool_avg75` past that comparison's own MDE with a CI excluding
               zero, THE SOURCE IS THE LEVER and the retrieval stage should be replaced.  If it
               loses past the MDE, the latent's ORACLE richness is unreachable by the shipped
               operator and the route needs a discriminator this project does not have.
  Promotion    NONE from this file.  It is a source comparison, not a pipeline change.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402
from s19 import qb_lib as QB                # noqa: E402

DRAWS = 512                      #: matched to L18's control size and the pool's evaluation budget
TOPM = 75                        #: the shipped tail size
TAG = os.environ.get("LS_TAG", "")


def _save(obj):
    path = os.path.join(RESULTS, "latentsel%s.json" % TAG)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def _avg(W, nat):
    """The shipped operator: superpose on the medoid, mean.  Returns (rmsd, medoid_rmsd)."""
    C, b = I.coordinate_average(np.asarray(W, float))
    return float(I.ca_rmsd(np.asarray(C, float), nat)), float(I.ca_rmsd(W[b], nat))


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, %d draws each" % (len(tg), DRAWS), flush=True)
    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        tgt = QB.target(pdb)
        nat, mu = tgt["nat"], tgt["mu"]
        N = 1 << n
        rng = SD.stable_rng(pdb, "s21latentsel")
        ob = QB.new_obj(tgt)

        #: ---- the LATENT source: DRAWS random bitstrings, decoded at their von Mises modes
        idx = rng.choice(N, min(DRAWS, N), replace=False)
        bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
        ar = np.arange(n)[None, :]
        _e, CA = ob.raw(mu[ar, bits, 0], mu[ar, bits, 1])
        sc = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(CA, tgt["i"], tgt["j"])), float)
        rr = I.kabsch_rmsd_batch(CA, nat)
        order = np.argsort(sc, kind="stable")

        lat_avg75, lat_avg75_med = _avg(CA[order[:TOPM]], nat)
        #: the objective-free averaging control: 75 draws chosen WITHOUT looking at the score
        rsel = rng.choice(len(CA), min(TOPM, len(CA)), replace=False)
        lat_avg75_rand, _ = _avg(CA[rsel], nat)
        P = I.pairwise_rmsd(CA)
        lat_medoid = float(rr[I.medoid(P)])

        #: ---- the POOL source: same operator, same basis (REBUILT from the windows' torsions)
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        Wr = np.asarray(I.build_ca(u["PHI"][p], u["PSI"][p]), float)
        scp = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(Wr, tgt["i"], tgt["j"])), float)
        op = np.argsort(scp, kind="stable")
        pool_avg75, _ = _avg(Wr[op[:TOPM]], nat)
        pool_rr = I.kabsch_rmsd_batch(Wr, nat)
        #: and the SHIPPED incumbent exactly as deployed -- window coordinates, not the rebuild
        W = np.asarray(u["W"], float)[p]
        scw = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(W, tgt["i"], tgt["j"])), float)
        ship_avg75, _ = _avg(W[np.argsort(scw, kind="stable")[:TOPM]], nat)

        rows.append({
            "pdb": pdb, "n": n, "fold": int(t["fold"]), "N_latent": N, "n_draws": int(len(CA)),
            "lat_rand1": float(rr[0]), "lat_argmin": float(rr[order[0]]),
            "lat_avg75": lat_avg75, "lat_avg75_rand": lat_avg75_rand,
            "lat_medoid": lat_medoid, "lat_oracle": float(rr.min()),
            "pool_avg75": pool_avg75, "pool_argmin": float(pool_rr[op[0]]),
            "pool_oracle": float(pool_rr.min()), "ship_avg75": ship_avg75})
        if (c + 1) % 10 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg), "draws": DRAWS})

    ok = len(rows) == len(tg) and all(np.isfinite(r["lat_avg75"]) and np.isfinite(r["pool_avg75"])
                                      for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "draws": DRAWS})
    report(rows)
    return rows


def _stat(d, rng, B=4000):
    d = np.asarray(d, float); k = len(d)
    se = float(d.std(ddof=1) / np.sqrt(k))
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return (d.mean(), se, 2.8016 * se,
            float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)), int((d < 0).sum()))


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "latentsel%s.json" % TAG)))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    rng = SD.stable_rng("latentsel", "report")
    fold = g("fold")
    print("\nn = %d targets, %d latent draws each.\n" % (len(rows), DRAWS))
    print("  %-18s%-22s%8s%9s" % ("arm", "basis", "RMSD", "median"))
    for k, b, in (("lat_rand1", "built chain, NULL"), ("lat_argmin", "built chain"),
                  ("lat_medoid", "built chain"), ("lat_oracle", "built chain, ORACLE"),
                  ("pool_argmin", "built chain"), ("pool_oracle", "built chain, ORACLE"),
                  ("lat_avg75", "POINT CLOUD"), ("lat_avg75_rand", "POINT CLOUD, null"),
                  ("pool_avg75", "POINT CLOUD"), ("ship_avg75", "POINT CLOUD, shipped")):
        print("  %-18s%-22s%8.3f%9.3f" % (k, b, g(k).mean(), np.median(g(k))))
    print("\n  PRIMARY and its two nulls (point cloud vs point cloud, same operator):")
    for lab, a, b in (("lat_avg75 - pool_avg75   THE PRIMARY", "lat_avg75", "pool_avg75"),
                      ("lat_avg75 - lat_avg75_rand  (did the OBJECTIVE work?)",
                       "lat_avg75", "lat_avg75_rand"),
                      ("lat_avg75_rand - pool_avg75  (source alone, no objective)",
                       "lat_avg75_rand", "pool_avg75"),
                      ("lat_avg75 - ship_avg75   (vs the SHIPPED window basis)",
                       "lat_avg75", "ship_avg75")):
        m, se, mde, lo, hi, w = _stat(g(a) - g(b), rng)
        fs = " ".join("%+.2f" % (g(a) - g(b))[fold == k].mean() for k in sorted(set(fold.astype(int))))
        print("    %-54s%+.3f SE %.3f MDE %.3f [%+.3f,%+.3f] %3dW/%3dL  folds %s"
              % (lab, m, se, mde, lo, hi, w, len(rows) - w, fs))
    #: WORKSTREAM D: the deployment footnote is a COMPOSITE of source swap AND basis change.
    #: Price the basis change on THIS operator so a reader can decompose it instead of trusting me.
    #: A's +0.011 A is a PER-MEMBER mean; the shift on a coordinate average of 75 is a different
    #: quantity -- the rebuild can change which member is the MEDOID, and the medoid sets the frame.
    m, se, mde, lo, hi, w = _stat(g("pool_avg75") - g("ship_avg75"), rng)
    print("\n  BASIS PRICE on this operator (makes the deployment footnote decomposable):")
    print("    %-54s%+.3f SE %.3f MDE %.3f [%+.3f,%+.3f] %3dW/%3dL"
          % ("pool_avg75(rebuilt) - ship_avg75(window)", m, se, mde, lo, hi, w, len(rows) - w))

    #: WORKSTREAM D: on n=9 targets, 512 draws WITHOUT replacement from a 512-element space is
    #: DETERMINISTIC -- it is the EXHAUSTIVE latent with zero variance, not a sample, and it hands
    #: those targets the ORACLE-richest possible set while every other target sees 6-50% of its own.
    #: The "random" label is wrong there, and wrong in the direction that FAVOURS the latent arm.
    #: A pre-declared sensitivity, not a caveat.
    nn = g("n"); dd = g("lat_avg75") - g("pool_avg75")
    print("\n  D's SENSITIVITY -- the primary, excluding targets where the draw IS the whole space:")
    for lab, msk in (("all targets", np.ones(len(nn), bool)),
                     ("n >= 10 (the draw is a true sample)", nn >= 10),
                     ("n >= 12 (draw <= 12.5% of the latent)", nn >= 12)):
        if msk.sum() > 2:
            m, se, mde, lo, hi, w = _stat(dd[msk], rng)
            print("    %-40s k=%3d  %+.3f SE %.3f MDE %.3f [%+.3f,%+.3f] %3dW/%3dL"
                  % (lab, int(msk.sum()), m, se, mde, lo, hi, w, int(msk.sum()) - w))

    print("\n  ORACLE ceilings of the two SOURCES (scoring only, not achievable):")
    m, se, mde, lo, hi, w = _stat(g("lat_oracle") - g("pool_oracle"), rng)
    print("    lat_oracle - pool_oracle  %+.3f SE %.3f MDE %.3f [%+.3f,%+.3f] %3dW/%3dL"
          % (m, se, mde, lo, hi, w, len(rows) - w))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
