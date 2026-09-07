"""s21/latentrank.py -- WHERE IN THE OBJECTIVE'S RANKING DOES THE BEST STRUCTURE ACTUALLY SIT?

L14 established that enumerating the ENTIRE latent and taking the exact argmin ties a
zero-evaluation pool, while the ORACLE over the identical 2**n modes is 1.60 A better.  L16 then
found that the one non-null lever in the selector is the READOUT Hamiltonian (-0.697 A) and not the
training one (+0.056 A, null).  Those two results meet on one question, which this file answers:

    THE BEST STRUCTURE IN THE LATENT EXISTS AND THE OBJECTIVE DOES NOT RETURN IT.
    HOW FAR DOWN THE OBJECTIVE'S OWN RANKING IS IT?

That number decides by arithmetic, not opinion, whether a better READOUT is a plausible route to the
1.60 A or a hopeless one:

  * ORACLE-best at the 0.01st percentile -> already near the top of the ordering.  A modest
                                           reranking, or a tail average instead of an argmin,
                                           reaches it.  The readout lever is LIVE.
  * ORACLE-best at the 40th percentile   -> it is in the BULK.  No monotone reranking of this
                                           objective can surface it, because the objective assigns
                                           it a typical score.  The readout lever is DEAD ON THIS
                                           OBJECTIVE and only a different H can help.

THREE QUANTITIES, ANSWERING THREE DIFFERENT QUESTIONS.

    oracle_rank_pct   percentile of the ORACLE-best config in the OBJECTIVE's ordering
                      -> "how deep must a reranker dig to reach the answer?"
    argmin_rmsd_pct   percentile of the OBJECTIVE's argmin in the TRUE RMSD ordering
                      -> "how good is what the objective returns, as a structure?"
    top{M}            best TRUE RMSD among the objective's top-M, M = 1, 8, 64, 512
    rand{M}           best TRUE RMSD among M configs drawn AT RANDOM -- THE MATCHED CONTROL
                      -> the CEILING on any reranker allowed to see the top-M and pick perfectly
                         inside it.  M=1 is the deployed argmin; M=1 -> M=512 is the ENTIRE budget
                         available to a smarter readout ON THIS OBJECTIVE.  ORACLE quantity.

THE MATCHED CONTROL, AND WHY THE top{M} LADDER IS MEANINGLESS WITHOUT IT.  `top{M}` is a MINIMUM
OVER M STRUCTURES, so the curve falls from M=1 to M=512 PARTLY BECAUSE 512 DRAWS BEAT 1 DRAW, with
no ordering skill involved whatever.  "A reranker could reach X at M=512" is then not separable
from "the minimum of 512 ARBITRARY latent configurations is X".  WORKSTREAM D registered this as
hazard D9 -- the same class as the min-of-N confound, applied to the new arm -- BEFORE any of these
numbers existed.  `rand{M}` is the min over a seeded random subset of the SAME size M, and

    THE VERTICAL GAP  rand{M} - top{M}  IS THE OBJECTIVE'S ORDERING SKILL,
    IN EXACTLY THE CURRENCY A RERANKER SPENDS.

That gap, not the `top{M}` level, is the quantity any claim about a reranker route must rest on.

THE MATCHED-SET ORACLE, and why it is here.  `latent_oracle` is a MIN OVER 2**n STRUCTURES, so it
improves with n by pure min-of-N with no discrimination content whatever.  Comparing the gap between
an n=13 target (8,192) and an n=16 target (65,536) therefore compares two different set sizes ON THE
ORACLE ARM ALONE, and any cross-n trend in the gap's MAGNITUDE is confounded by that.  WORKSTREAM D
identified this in the coordinator's own task design before running it.  `oracle_8192` re-takes the
minimum over a FIXED, seeded 8,192-configuration subset, identical in size at every n, which makes
the cross-n comparison clean.  At n <= 13 the subset is the whole latent and the two columns
coincide by construction -- which is itself the check that the column is wired correctly.

OPERATOR FORKS, DECLARED PER BRIEF SS7 RULE 0 -- the rule this lane's own failure produced.

    functional     DECLARED shipped Bayes risk (`I.shipped_score`), the deployed selector.
                   NOT TAKEN the squared-distance functional.  Both computed, both reported; if
                   they disagree the DECLARED one stands and the disagreement is recorded as a
                   finding rather than swapped in.
    basis          DECLARED built chain from torsions (`Obj.raw` -> `build_ca_exact`).  The whole
                   file lives INSIDE the latent, so the window-vs-rebuild fork that broke L14's
                   primary cannot arise in this comparison at all.
    readout        DECLARED percentile/rank.  NOT TAKEN mean rank, which an outlier-heavy tail
                   would dominate.  Median across targets is printed beside every mean.
    normalisation  DECLARED percentile of 2**n, so targets of different latent size are
                   comparable.  NOT TAKEN raw rank, which is a function of n by construction.
    null           An objective with NO ordering skill puts the ORACLE-best at the 50th percentile
                   in expectation.  **50.0 is the number every result below is read against**, and
                   it is stated here before the run.

TIES.  An `argmin` on a tied score silently reads the sort order -- the trap that once manufactured
a 1.386 A winner in this project.  The tied set at the objective's minimum is counted and reported,
and `argmin_rmsd_pct` is taken against the MEAN outcome over the tied set rather than against
whichever member the sort happened to place first.

This file re-derives its own structures rather than reading L14's artefact, so `latent_oracle` here
is an INDEPENDENT recomputation of L14's column and must reproduce it exactly.  That agreement is
checked and printed; a mismatch is a defect in one of the two files.
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

MIN_N = int(os.environ.get("LR_MIN_N", 0))
MAX_N = int(os.environ.get("LR_MAX_N", 13))
TAG = os.environ.get("LR_TAG", "")
CHUNK = 2048
MS = (1, 8, 64, 512)
MATCHED = 8192                   #: the matched-set ORACLE size = the samplers' budget
RAND_REP = 32                    #: repeats of the random-M control, averaged


def _save(obj):
    path = os.path.join(RESULTS, "latentrank%s.json" % TAG)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def run():
    tg = [t for t in I.targets() if MIN_N <= int(t["n"]) <= MAX_N]
    print("targets (%d <= n <= %d): %d" % (MIN_N, MAX_N, len(tg)), flush=True)
    rows = []
    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        tgt = QB.target(pdb)
        nat, mu = tgt["nat"], tgt["mu"]
        N = 1 << n
        ob = QB.new_obj(tgt)
        E = np.empty(N); EB = np.empty(N); RR = np.empty(N)
        ar = np.arange(n)[None, :]
        for s0 in range(0, N, CHUNK):
            k = np.arange(s0, min(s0 + CHUNK, N), dtype=np.int64)
            bits = ((k[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
            e, CA = ob.raw(mu[ar, bits, 0], mu[ar, bits, 1])
            E[k] = e
            EB[k] = I.shipped_score(tgt["dg"], I.pair_dists(CA, tgt["i"], tgt["j"]))
            RR[k] = I.kabsch_rmsd_batch(CA, nat)

        #: matched-set ORACLE: identical 8192 draw size at every n, so cross-n is unconfounded
        sub = (np.arange(N) if N <= MATCHED else
               SD.stable_rng(pdb, "latentrank-matched").choice(N, MATCHED, replace=False))
        row = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "N": N,
               "latent_oracle": float(RR.min()),
               "oracle_8192": float(RR[sub].min()),
               "matched_is_full": bool(N <= MATCHED)}
        #: D9's matched control: min over M RANDOM configs, averaged over RAND_REP draws, so the
        #: min-of-M component of the top-M ladder is priced instead of being read as ordering skill.
        rr_rng = SD.stable_rng(pdb, "latentrank-randM")
        for m in MS:
            row["rand%d" % m] = float(np.mean(
                [RR[rr_rng.choice(N, min(m, N), replace=False)].min() for _ in range(RAND_REP)]))
        for tag, sc in (("bayes", EB), ("sq", E)):
            ordr = np.argsort(sc, kind="stable")
            rank_of = np.empty(N, np.int64); rank_of[ordr] = np.arange(N)
            best = int(np.argmin(RR))
            tied = np.flatnonzero(sc == sc[ordr[0]])          # ties at the objective's minimum
            row["oracle_rank_pct_%s" % tag] = 100.0 * (rank_of[best] + 0.5) / N
            row["argmin_rmsd_pct_%s" % tag] = 100.0 * ((RR < RR[tied].mean()).sum() + 0.5) / N
            row["n_tied_at_min_%s" % tag] = int(len(tied))
            for m in MS:
                row["top%d_%s" % (m, tag)] = float(RR[ordr[:m]].min())
        rows.append(row)
        if (c + 1) % 10 == 0:
            print("  %d/%d" % (c + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    ok = len(rows) == len(tg) and len(tg) > 0 and all(
        np.isfinite(r["oracle_rank_pct_bayes"]) and np.isfinite(r["oracle_8192"]) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg)})
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = []
        for t in ([TAG] if TAG else ["", "_hi"]):
            f = os.path.join(RESULTS, "latentrank%s.json" % t)
            if os.path.exists(f):
                rows += json.load(open(f))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    nn = g("n")
    print("\nn = %d targets, n in [%d,%d], EXHAUSTIVE over 2**n modes." % (len(rows), nn.min(), nn.max()))
    print("NULL: an objective with no ordering skill puts the ORACLE-best at percentile 50.0.\n")
    for tag, lab in (("bayes", "DECLARED Bayes risk"), ("sq", "not-taken squared")):
        o = g("oracle_rank_pct_%s" % tag); a = g("argmin_rmsd_pct_%s" % tag)
        print("  %-20s ORACLE-best sits at objective pct : mean %6.2f  median %6.2f  "
              "min %.4f  max %.2f" % (lab, o.mean(), np.median(o), o.min(), o.max()))
        print("  %-20s objective's argmin sits at RMSD pct: mean %6.2f  median %6.2f"
              % ("", a.mean(), np.median(a)))
    print("\n  THE READOUT CEILING -- best true RMSD inside the objective's top-M, against D9's")
    print("  MATCHED CONTROL (min over M RANDOM configs).  GAP = randM - topM = the ORDERING SKILL.")
    print("    %-6s%10s%10s%10s%10s" % ("M", "top-M", "squared", "rand-M", "GAP"))
    for m in MS:
        tb = g("top%d_bayes" % m).mean(); rd = g("rand%d" % m).mean()
        print("    %-6d%10.3f%10.3f%10.3f%10.3f%s" % (
            m, tb, g("top%d_sq" % m).mean(), rd, rd - tb,
            "   <- the DEPLOYED argmin" if m == 1 else ""))
    print("    %-6s%10.3f%10.3f   <- ORACLE over the WHOLE latent" %
          ("all", g("latent_oracle").mean(), g("latent_oracle").mean()))
    print("    %-6s%10.3f%10s   <- ORACLE over a MATCHED 8192 subset (cross-n clean)" %
          ("8192", g("oracle_8192").mean(), ""))
    print("\n  ties at the objective's minimum: median %d (Bayes), %d (squared)"
          % (np.median(g("n_tied_at_min_bayes")), np.median(g("n_tied_at_min_sq"))))
    f = os.path.join(RESULTS, "latentfull.json")
    if os.path.exists(f):
        d = {r["pdb"]: r["latent_oracle"] for r in json.load(open(f))["rows"]}
        b = [(r["latent_oracle"], d[r["pdb"]]) for r in rows if r["pdb"] in d]
        if b:
            x = np.array(b); dd = float(np.abs(x[:, 0] - x[:, 1]).max())
            print("\n  CROSS-CHECK vs L14's INDEPENDENT enumeration, %d shared targets: "
                  "max |diff| = %.2e  %s" % (len(b), dd, "AGREE" if dd < 1e-9 else "*** MISMATCH ***"))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
