"""s24/e_referaudit.py -- LANE E AUDIT OF L5 PRIMARY 2.  IS THE MISMATCHED-SAME-LENGTH DISTOGRAM
THE RIGHT FLOOR FOR beta, AND DOES beta = 0.520 vs 0.352 SURVIVE A HARDER ONE?

WHAT IS UNDER AUDIT.  `s24/referent.py` measures `beta = <eC, eP> / <eP, eP>` with
`eP = Dhat - Dt` (the prior's own pair-distance error) and `eC = Dc - Dt` (the emitted cloud's),
and floors it with a PLACEBO: the same construction using a mismatched same-length target's
`Dhat`.  L5 reports real 0.5201 against floor 0.3521, excess +0.1680 at 2.4x its own MDE.
The coordinator asked Lane E specifically whether that placebo is the right floor and whether the
excess survives a placebo drawn from a different length or from a shuffled pair index.

THE TWO THINGS THE SHIPPED PLACEBO DOES NOT CONTROL.

  (i)  beta IS NOT NORM-FREE.  beta = cos * |eC| / |eP|.  A placebo prior that predicts Dt WORSE
       has a larger |eP| and a mechanically smaller beta even at identical alignment.  The
       norm-free statistic is the cosine, which referent.py already records; the two must be read
       together, and the share of the beta gap that is norm rather than alignment must be stated.

  (ii) A SINGLE RANDOM PEER IS NOT THE TIGHTEST FLOOR.  A floor should keep everything generic --
       length, the sequence-separation profile that dominates a CA distance, protein-likeness --
       and remove only target-specific information.  One random peer is a noisy draw from that
       class.  Four harder floors are built here, all on the same eC and the same Dt:

         P_peer      the shipped floor, reproduced bit-for-bit from referent.py's own RNG stream
         P_meanlen   the MEAN Dhat over all OTHER same-length targets.  Maximal genericness,
                     zero target-specific information, and less noisy than any single peer --
                     the plausible zero-information control per the standing rule that such a
                     control must be plausible rather than degenerate.
         P_bestpeer  the same-length peer whose Dhat is CLOSEST to the real Dhat.  ADVERSARIAL:
                     it maximises the placebo's alignment, so it is the strictest available floor.
         P_difflen   a DIFFERENT-length target's Dhat, restricted to this target's own (i,j) pair
                     set by index mapping, so the pair geometry is preserved and only the length
                     match is broken.  The coordinator's first requested check.
         P_shuffle   the real Dhat with its pair index permuted.  The coordinator's second
                     requested check.  PREDICTED to be a LOOSER floor, because it destroys the
                     sequence-separation structure that is most of what a CA distogram knows --
                     if so, the shipped choice was the conservative one and that should be said.

OPERATOR FORKS (Lane E, no stake).
    functional     DECLARED the shipped Bayes-risk per-pair optimum `Dhat = grid[argmin(risk)]`,
                   identical to referent.py, for the real arm and every placebo.  NOT TAKEN the
                   posterior mean, which referent.py already reports as its declared secondary.
    basis          DECLARED distance space on the score's own pair set (min_sep=2), no structure
                   operations, identical to the file under audit.  NOT TAKEN a coordinate basis.
    readout        DECLARED beta AND the norm-free cosine for every arm, always together.
                   NOT TAKEN beta alone, which is the fork that hides (i).
    normalisation  DECLARED every placebo evaluated against the SAME eC and the SAME Dt as the
                   real arm, so the arms differ in the PRIOR and in nothing else.
                   NOT TAKEN re-deriving the emitted cloud per placebo.
    null           DECLARED five floors spanning looser-than-shipped to adversarially tighter.
                   NOT TAKEN zero, which is the wrong null when both errors share a native.
    THE LABEL      DECLARED continuous beta and cosine.  NOT TAKEN a binarised "does it transfer".

  H0 (Lane E's null, which L5 must beat): the +0.168 excess is an artefact of the placebo prior
     being a WORSE predictor (larger |eP|), or of a single noisy peer, and shrinks to nothing
     against a generic or adversarial floor.
  Falsifier of H0: the excess survives P_meanlen and P_bestpeer, in the cosine as well as in beta.

  Natives are read for EVALUATION ONLY.
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
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from s24 import stats_lib as ST           # noqa: E402

TOPM = 75
OUT = os.path.join(RES, "e_referaudit.json")
ARMS = ("real", "P_peer", "P_meanlen", "P_bestpeer", "P_difflen", "P_shuffle")
NEED = ("real", "P_peer", "P_meanlen", "P_bestpeer", "P_shuffle")


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _dhat(dg):
    """referent.py's functional, unchanged: the per-pair minimiser of the shipped Bayes risk."""
    return np.asarray(dg["grid"], float)[np.argmin(np.asarray(dg["risk"], float), axis=1)]


def _stats(Dc, Dt, Dh):
    eP = Dh - Dt
    eC = Dc - Dt
    den = float((eP * eP).sum())
    if den <= 0:
        return {"beta": float("nan"), "cos": float("nan"), "rms_hat_to_nat": 0.0,
                "norm_eP": 0.0}
    return {"beta": float((eC * eP).sum() / den),
            "cos": float((eC * eP).sum() / np.sqrt(float((eC * eC).sum()) * den)),
            "rms_hat_to_nat": float(np.sqrt(((Dh - Dt) ** 2).mean())),
            #: PRIMARY 1 through the same floor panel: is the emitted cloud closer to THIS
            #: prediction than the native is?  L5 shows it for the real distogram; the question
            #: an auditor must ask is whether it is also true of a prediction that carries no
            #: target-specific information at all.
            "rms_cloud_to_hat": float(np.sqrt(((Dc - Dh) ** 2).mean())),
            "p1_gap": float(np.sqrt(((Dc - Dh) ** 2).mean()) - np.sqrt(((Dt - Dh) ** 2).mean())),
            "norm_eP": float(np.sqrt(den))}


def run():
    tg = I.targets()
    by_len = {}
    for t in tg:
        by_len.setdefault(int(t["n"]), []).append(t["pdb"])
    lens = sorted(by_len)
    dh_cache = {}

    def dhat_of(pdb):
        if pdb not in dh_cache:
            uu = I.load_univ(pdb)
            dh_cache[pdb] = _dhat(I.distogram(pdb, uu["seq"], uu["fold"]))
        return dh_cache[pdb]

    rows = []
    print("targets: %d, lengths %s" % (len(tg), lens), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        n_res = int(u["n"]); i, j = I.pair_index(n_res)
        dg = I.distogram(pdb, u["seq"], u["fold"])
        idx = I.pool_idx(u); Wpool = Wall[idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        A = Wpool[np.argsort(sc, kind="stable")[:TOPM]]
        C = _avg(A)

        Dt = I.pair_dists(nat[None], i, j)[0]
        Dc = I.pair_dists(C[None], i, j)[0]
        Dh = _dhat(dg)

        row = {"pdb": pdb, "n": n_res, "fold": int(u["fold"]),
               "rmsd": float(I.ca_rmsd(C, nat)),
               "real": _stats(Dc, Dt, Dh)}

        peers = [p for p in by_len[n_res] if p != pdb]
        if peers:
            #: P_peer -- referent.py's own choice, same RNG stream, bit-for-bit
            rng = SD.stable_rng("referent", pdb)
            other = peers[int(rng.integers(0, len(peers)))]
            row["placebo_pdb"] = other
            row["placebo_same_fold"] = bool(int(I.load_univ(other)["fold"]) == int(u["fold"]))
            row["P_peer"] = _stats(Dc, Dt, dhat_of(other))

            #: P_meanlen -- the generic same-length prior: mean Dhat over ALL other peers
            row["P_meanlen"] = _stats(Dc, Dt, np.mean([dhat_of(p) for p in peers], axis=0))

            #: P_bestpeer -- ADVERSARIAL: the peer whose Dhat is closest to the real Dhat
            dd = [float(((dhat_of(p) - Dh) ** 2).sum()) for p in peers]
            bp = peers[int(np.argmin(dd))]
            row["bestpeer_pdb"] = bp
            row["P_bestpeer"] = _stats(Dc, Dt, dhat_of(bp))

        #: P_difflen -- a LONGER target's Dhat restricted to this target's own (i,j) pairs
        longer = [p for L in lens if L > n_res for p in by_len[L]]
        if longer:
            rl = SD.stable_rng("e_referaudit", pdb, "difflen")
            o2 = longer[int(rl.integers(0, len(longer)))]
            n2 = int(I.load_univ(o2)["n"])
            i2, j2 = I.pair_index(n2)
            pos = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(i2, j2))}
            sel = np.array([pos[(int(a), int(b))] for a, b in zip(i, j)], int)
            row["difflen_pdb"] = o2; row["difflen_n"] = n2
            row["P_difflen"] = _stats(Dc, Dt, dhat_of(o2)[sel])

        #: P_shuffle -- the real Dhat, pair index permuted
        rs = SD.stable_rng("e_referaudit", pdb, "shuffle")
        row["P_shuffle"] = _stats(Dc, Dt, Dh[rs.permutation(len(Dh))])

        rows.append(row)
        if (c_i + 1) % 20 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            ST.save_atomic(OUT, {"rows": rows}, complete_keys=NEED, rows=rows,
                           n_expected=len(tg), module_file=__file__)

    ST.save_atomic(OUT, {"rows": rows, "topm": TOPM}, complete_keys=NEED, rows=rows,
                   n_expected=len(tg), module_file=__file__)
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rows = [r for r in rows if all(k in r for k in NEED)]
    fold = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]

    def G(a, k):
        return np.array([r[a][k] if a in r else np.nan for r in rows], float)

    print("\n" + "=" * 96)
    print("LANE E AUDIT OF L5 PRIMARY 2 (beta).  n = %d with a same-length peer." % len(rows))
    print("=" * 96)

    print("\n  (1) EVERY FLOOR, ON THE SAME eC AND THE SAME NATIVE.  beta AND the norm-free cosine.")
    print("    %-38s%10s%10s%10s%10s%12s" % ("arm", "beta mn", "beta md", "cos mn", "cos md",
                                             "|eP| (RMS)"))
    LAB = {"real": "REAL distogram",
           "P_peer": "P_peer      shipped floor (1 random peer)",
           "P_meanlen": "P_meanlen   generic same-length mean",
           "P_bestpeer": "P_bestpeer  ADVERSARIAL closest peer",
           "P_difflen": "P_difflen   different length, pairs mapped",
           "P_shuffle": "P_shuffle   real Dhat, pair index permuted"}
    for a in ARMS:
        b = G(a, "beta"); c = G(a, "cos"); e = G(a, "rms_hat_to_nat")
        f = np.isfinite(b)
        if not f.any():
            continue
        print("    %-38s%10.4f%10.4f%10.4f%10.4f%12.4f"
              % (LAB[a], b[f].mean(), np.median(b[f]), c[f].mean(), np.median(c[f]), e[f].mean()))

    print("\n  (2) REAL MINUS EACH FLOOR, PAIRED, fold-clustered.  beta first, then the COSINE.")
    for stat in ("beta", "cos"):
        print("\n    --- %s ---" % stat)
        for a in ARMS[1:]:
            v = G(a, stat); r0 = G("real", stat)
            f = np.isfinite(v) & np.isfinite(r0)
            if f.sum() < 20:
                continue
            res = ST.compare(r0[f], v[f], fold[f], names=[p for p, k in zip(pdbs, f) if k],
                             label="real - %s (%s)" % (a, stat))
            print("      %-12s %+.4f  med %+.4f  SE %.4f  MDE %.4f  e/MDE %+.2f  "
                  "fold[%+.4f,%+.4f]  %3dW/%3dL  n=%d  %s"
                  % (a, res["effect"], res["median_effect"], res["se"], res["mde"],
                     res["effect_over_mde"], res["ci95_fold"][0], res["ci95_fold"][1],
                     res["n_better"], res["n_worse"], int(f.sum()), res["verdict"]))
    print("\n      NOTE ON SIGN: `compare` returns real - floor, so a POSITIVE effect with the")
    print("      fold CI excluding zero is L5's claim SURVIVING.  `verdict` reads WORSE for a")
    print("      positive effect because the module's convention is 'negative is better' for an")
    print("      RMSD endpoint; here positive IS the result.  Read the CI, not the word.")

    print("\n  (3) HOW MUCH OF THE beta GAP IS NORM RATHER THAN ALIGNMENT?")
    print("      beta = cos * |eC| / |eP|, so a floor whose prior predicts WORSE has a larger")
    print("      |eP| and a mechanically smaller beta at identical alignment.")
    br = G("real", "beta"); cr = G("real", "cos")
    for a in ARMS[1:]:
        b = G(a, "beta"); c = G(a, "cos"); f = np.isfinite(b)
        if f.sum() < 20:
            continue
        gb = (br[f] - b[f]).mean(); gc = (cr[f] - c[f]).mean()
        ratio = (br[f].mean() / b[f].mean()) / (cr[f].mean() / c[f].mean()) if b[f].mean() else np.nan
        print("      %-12s beta gap %+.4f   cos gap %+.4f   norm-only multiplier %.3f"
              % (a, gb, gc, ratio))

    print("\n  (3b) PRIMARY 1 THROUGH THE SAME PANEL.  L5 P1 is RMS|Dc-Dhat| - RMS|Dt-Dhat| < 0:")
    print("       'the emitted cloud agrees with the prediction better than the native does'.")
    print("       Is that true of a prediction carrying NO target-specific information?")
    print("    %-38s%12s%12s%12s" % ("arm", "cloud->hat", "nat->hat", "P1 gap"))
    for a in ARMS:
        gp = G(a, "p1_gap"); ch = G(a, "rms_cloud_to_hat"); nh = G(a, "rms_hat_to_nat")
        f = np.isfinite(gp)
        if f.sum() < 20:
            continue
        res = ST.compare(ch[f], nh[f], fold[f], names=[p for p, k in zip(pdbs, f) if k],
                         label="P1 %s" % a)
        print("    %-38s%12.4f%12.4f%+12.4f   SE %.4f MDE %.4f fold[%+.4f,%+.4f] %3dW/%3dL n=%d"
              % (LAB[a], ch[f].mean(), nh[f].mean(), res["effect"], res["se"], res["mde"],
                 res["ci95_fold"][0], res["ci95_fold"][1], res["n_better"], res["n_worse"],
                 int(f.sum())))
    print("      (a NEGATIVE P1 gap on a zero-information arm means P1 is measuring TYPICALITY,")
    print("       not agreement with THIS distogram.)")

    print("\n  (4) THE SHIPPED PLACEBO'S OWN NUISANCE: is the peer in the SAME leave-fold-out fold?")
    sf = np.array([bool(r.get("placebo_same_fold", False)) for r in rows])
    b = G("P_peer", "beta")
    print("      same-fold peers %d/%d.  beta(P_peer) same-fold %.4f  vs  cross-fold %.4f"
          % (int(sf.sum()), len(rows), b[sf].mean() if sf.any() else float('nan'),
             b[~sf].mean() if (~sf).any() else float('nan')))
    print("      (a cross-fold peer's distogram model may have SEEN this target, which would make")
    print("       the floor MORE informative and the excess CONSERVATIVE, not inflated.)")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
