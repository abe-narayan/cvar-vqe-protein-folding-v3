"""s24/e_bowaudit.py -- LANE E ADVERSARIAL AUDIT OF LEDGER L2(b)/(c)/(d).

WHAT IS UNDER AUDIT.  `s24/biasalign.py` reports that the matched-size mixture curve
"bows 0.226 A below the straight line between its endpoints", and the ledger reads that
bow as INDEPENDENT confirmation of bias non-parallelism, corroborated by a fitted
c_eff = 0.655 against a directly measured cosine of 0.647.  Two lanes are now working to a
spec (q* = 1.274 -> "<= 3.9 A standalone") derived from that fit.

THE THREE ATTACKS, each with a MEASURED null rather than an argument.

  (1) THE UNDECLARED COMPOSITION FORK.  biasalign builds every interior mixture as
          sel = concatenate([A[:a], C1[:b]])
      `A` is `Wpool[argsort(score)[:75]]` and is therefore held IN SCORE ORDER, so `A[:a]`
      is the TOP-a BY SCORE, not an a-subset of A.  The interior arms therefore differ from
      the endpoints in TWO ways -- composition AND the selectivity of the retained score
      half -- while the declared `readout` fork asserts they "differ in composition and in
      nothing else".  The alternative not taken (a uniform a-subset of the same 75) is not
      named anywhere in the module.  NULL: run the identical ladder with a uniform a-subset
      and price the difference in the bow.

  (2) THE OPERATOR NULL.  `_avg` recomputes the medoid and the superposition frame at every
      mixture, and the mixed set's spread changes with composition, so the emitted cloud is
      NOT a linear interpolation of the two parent clouds even when the biases are exactly
      parallel.  NULL: the WITHIN-SOURCE ladder C1[:a] + C2[:b] -- two independent uniform
      draws from one source, matched quality, prefixes of a random draw (so attack (1)
      cannot act), and a directly measured cosine of 0.933.  Any bow there, and any gap
      between the fitted c_eff and 0.933, is pure operator.

  (3) THE SHARED REFERENT.  The "two independent routes" both consume the SAME two bias
      vectors.  S(0) = |e_A|/sqrt(n) and S(1) = |e_C|/sqrt(n) are the fit's own endpoints,
      and any interior S(lambda) is, up to the operator, |(1-lambda)e_A + lambda e_C|/sqrt(n).
      Inverting the quadratic model for c is then algebra on the same vectors the cosine is
      computed from, not a second measurement.  NULL: generate the ladder analytically from
      e_A and e_C with NO operator at all, fit c_eff to it, and report how much of the
      claimed agreement is recovered by construction.

OPERATOR FORKS for this audit (enumerated by Lane E, which has no stake in the outcome).

    functional     DECLARED the shipped Bayes-risk distogram score for A, unchanged from
                   biasalign, so the reproduction arm is bit-comparable.  NOT TAKEN any
                   re-scoring.
    basis          DECLARED point cloud throughout, bias vectors in the native frame,
                   RMSD in each set's own medoid frame -- identical to the file under audit,
                   because an audit that changes the basis measures its own change.
                   NOT TAKEN a built chain, NOT TAKEN a shared frame.
    readout        DECLARED uniform mean of exactly 75 members at every ladder point of every
                   ladder.  NOT TAKEN unmatched sizes.
    normalisation  DECLARED the bow reported BOTH in absolute A and as a fraction of that
                   ladder's own S(0), because ladders with different endpoints cannot be
                   compared in absolute A.  NOT TAKEN absolute A alone, which is the
                   comparison the ledger makes.
    null           DECLARED three: the uniform-subset ladder (composition fork), the
                   within-source C1/C2 ladder (operator), and the analytic no-operator ladder
                   (shared referent).  NOT TAKEN "zero bow", which assumes the conclusion.
    THE LABEL      DECLARED the bow in A and the fitted c_eff, the two quantities the ledger
                   actually cites.  NOT TAKEN a binarised "does the curve bow".

  H0 (Lane E's null, which the ledger must beat): the observed -0.226 A bow is reproduced by
     the composition fork plus the operator, with zero bias independence required.
  Falsifier of H0: the within-source ladder recovers its own directly measured cosine, AND
     removing the score prefix leaves the A/C1 bow essentially unchanged.

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

TOPM = 75
MIX = [(75, 0), (60, 15), (50, 25), (38, 37), (25, 50), (0, 75)]
NRND = 4          # independent uniform-subset permutations for the composition-fork null


def _save(o, name="e_bowaudit.json"):
    p = os.path.join(RES, name); t = p + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, p)


def _kabsch_R(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _bias(C, nat):
    R = _kabsch_R(C, nat)
    return (C - C.mean(0)) @ R.T - (nat - nat.mean(0))


def _cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


def _ladder(src1, src2, nat):
    """RMSD at each matched-size mixture point of two 75-member sources."""
    out = {}
    for a, b in MIX:
        if a == 0:
            sel = src2[:b]
        elif b == 0:
            sel = src1[:a]
        else:
            sel = np.concatenate([src1[:a], src2[:b]], axis=0)
        out["m%d_%d" % (a, b)] = float(I.ca_rmsd(_avg(sel), nat))
    return out


def run():
    tg = I.targets()
    rows = []
    print("targets: %d, m=%d, mixtures %s, %d subset perms" % (len(tg), TOPM, MIX, NRND), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        Wall = np.asarray(u["W"], float); nat = np.asarray(u["nat_ca"], float)
        N = len(Wall); n_res = int(u["n"])
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n_res)
        idx = I.pool_idx(u); Wpool = Wall[idx]
        sc = np.asarray(I.shipped_score(dg, I.pair_dists(Wpool, i, j)), float)
        o = np.argsort(sc, kind="stable")
        A = Wpool[o[:TOPM]]                       # HELD IN SCORE ORDER, exactly as biasalign

        #: reproduce biasalign's RNG stream so the reproduction arm is bit-comparable
        rb = SD.stable_rng("biasalign", pdb)
        c1 = rb.choice(N, TOPM, replace=False)
        c2 = rb.choice(N, TOPM, replace=False)
        C1 = Wall[c1]; C2 = Wall[c2]

        eA = _bias(_avg(A), nat); eC1 = _bias(_avg(C1), nat); eC2 = _bias(_avg(C2), nat)

        ra = SD.stable_rng("e_bowaudit", pdb)
        row = {"pdb": pdb, "n": n_res, "fold": int(u["fold"]), "n_universe": int(N),
               "cos_A_C": _cos(eA, eC1), "cos_C1_C2": _cos(eC1, eC2),
               "norm_A": float(np.linalg.norm(eA)), "norm_C1": float(np.linalg.norm(eC1)),
               "norm_C2": float(np.linalg.norm(eC2))}

        #: L_obs -- exact reproduction of the file under audit
        row["L_obs"] = _ladder(A, C1, nat)

        #: L_rnd -- identical, but the score half is a UNIFORM a-subset of the same 75
        rl = []
        for _ in range(NRND):
            perm = ra.permutation(TOPM)
            rl.append(_ladder(A[perm], C1, nat))
        row["L_rnd"] = {k: float(np.mean([d[k] for d in rl])) for k in rl[0]}
        row["L_rnd_sd"] = {k: float(np.std([d[k] for d in rl], ddof=1)) for k in rl[0]}

        #: L_ctl -- the WITHIN-SOURCE operator null: two independent uniform draws
        row["L_ctl"] = _ladder(C1, C2, nat)

        #: the pure m-ladder: is "top-a by score" better or worse than "all 75"?
        row["m_prefix"] = {str(a): float(I.ca_rmsd(_avg(A[:a]), nat)) for a, _ in MIX if a}
        row["m_subset"] = {str(a): float(np.mean([I.ca_rmsd(_avg(A[ra.permutation(TOPM)[:a]]), nat)
                                                  for _ in range(NRND)])) for a, _ in MIX if a}

        #: the ANALYTIC ladder: pure geometry, no operator, no subsetting
        an = {}
        for a, b in MIX:
            lam = b / float(TOPM)
            an["m%d_%d" % (a, b)] = float(np.linalg.norm((1 - lam) * eA + lam * eC1)
                                          / np.sqrt(n_res))
        row["L_ana"] = an
        rows.append(row)

        if (c_i + 1) % 10 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})

    need = ("L_obs", "L_rnd", "L_ctl", "L_ana", "m_prefix", "cos_A_C", "cos_C1_C2")
    ok = len(rows) == len(tg) and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg), "topm": TOPM,
           "mix": MIX, "nrnd": NRND})
    report(rows)
    return rows


# --------------------------------------------------------------------------- analysis
def fit_c(S, lams):
    """Closed-form LS fit of c in  S(l)^2/S(0)^2 = (1-l)^2 + l^2 q^2 + 2 l (1-l) q c.

    S is the ladder in ascending lambda; interior points only enter the fit.  q = S(1)/S(0).
    Linear in c, so there is no optimiser fork to declare.
    """
    S = np.asarray(S, float); lams = np.asarray(lams, float)
    q = S[-1] / S[0]
    y, x = [], []
    for k in range(1, len(S) - 1):
        l = lams[k]
        y.append((S[k] / S[0]) ** 2 - (1 - l) ** 2 - l ** 2 * q ** 2)
        x.append(2 * l * (1 - l) * q)
    x = np.asarray(x); y = np.asarray(y)
    return float((x * y).sum() / (x * x).sum()), float(q)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "e_bowaudit.json")))["rows"]
    keys = ["m%d_%d" % (a, b) for a, b in MIX]
    lams = np.array([b / float(TOPM) for _, b in MIX])
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("e_bowaudit", "rep")
    F = sorted(set(fold.tolist()))

    def st(x):
        x = np.asarray(x, float); k = len(x)
        se = x.std(ddof=1) / np.sqrt(k)
        fs = [np.concatenate([x[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (x.mean(), np.median(x), se, 2.8016 * se,
                float(np.percentile(fs, 2.5)), float(np.percentile(fs, 97.5)), int((x < 0).sum()))

    def bows(tag):
        M = np.array([[r[tag][k] for k in keys] for r in rows], float)   # (n, 6)
        line = (1 - lams)[None, :] * M[:, :1] + lams[None, :] * M[:, -1:]
        return M, M - line

    print("\n" + "=" * 96)
    print("LANE E AUDIT OF L2(b)/(c).  n = %d." % len(rows))
    print("=" * 96)

    LAD = [("L_obs", "A[:a] + C1[:b]        AS SHIPPED (score prefix)"),
           ("L_rnd", "A[perm][:a] + C1[:b]  COMPOSITION-FORK NULL (uniform subset)"),
           ("L_ctl", "C1[:a] + C2[:b]       WITHIN-SOURCE OPERATOR NULL"),
           ("L_ana", "analytic |(1-l)eA + l eC1|  NO OPERATOR AT ALL")]

    print("\n  (1) THE LADDERS, and the bow each one shows")
    for tag, lab in LAD:
        M, B = bows(tag)
        print("\n    %s" % lab)
        print("      %-10s%9s%12s%12s%10s" % ("mixture", "RMSD", "vs 75/0", "vs line", "bow/S(0)"))
        for c, k in enumerate(keys):
            m, md, se, mde, flo, fhi, w = st(B[:, c])
            print("      %-10s%9.4f%+12.4f%+12.4f%+10.4f   SE %.4f MDE %.4f fold[%+.4f,%+.4f]"
                  % (k, M[:, c].mean(), (M[:, c] - M[:, 0]).mean(), m,
                     (B[:, c] / M[:, 0]).mean(), se, mde, flo, fhi))

    print("\n  (2) THE COMPOSITION FORK, priced.  L_obs bow MINUS L_rnd bow, paired per target.")
    _, Bo = bows("L_obs"); _, Br = bows("L_rnd")
    for c, k in enumerate(keys[1:-1], start=1):
        m, md, se, mde, flo, fhi, w = st(Bo[:, c] - Br[:, c])
        v = "SIGNIFICANT" if (abs(m) > mde and (fhi < 0 or flo > 0)) else "ns"
        print("      %-10s %+.4f  median %+.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f] %3dW/%3dL %s"
              % (k, m, md, se, mde, flo, fhi, w, len(rows) - w, v))
    print("      (negative = the SCORE PREFIX manufactures extra bow that composition does not)")

    print("\n  (3) THE PURE m-LADDER: what does keeping only the top-a by score do on its own?")
    for tag, lab in (("m_prefix", "top-a by score"), ("m_subset", "uniform a-subset")):
        v = np.array([[r[tag][str(a)] for a, _ in MIX if a] for r in rows], float)
        print("      %-18s " % lab + "  ".join("a=%d %.4f" % (a, v[:, c].mean())
                                               for c, (a, _) in enumerate([x for x in MIX if x[0]])))
    vp = np.array([[r["m_prefix"][str(a)] for a, _ in MIX if a] for r in rows], float)
    vs = np.array([[r["m_subset"][str(a)] for a, _ in MIX if a] for r in rows], float)
    print("      paired top-a MINUS uniform-a-subset (negative = the prefix is a QUALITY gain):")
    for c, (a, _) in enumerate([x for x in MIX if x[0]]):
        m, md, se, mde, flo, fhi, w = st(vp[:, c] - vs[:, c])
        print("        a=%-4d %+.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f]  %3dW/%3dL  %s"
              % (a, m, se, mde, flo, fhi, w, len(rows) - w,
                 "REAL" if abs(m) > mde and (fhi < 0 or flo > 0) else "ns"))

    print("\n  (4) FITTED c_eff PER LADDER vs THE DIRECTLY MEASURED COSINE FOR THAT PAIR")
    print("      %-38s%10s%10s%10s%10s%9s" % ("ladder", "c_eff mn", "c_eff md", "direct mn",
                                              "direct md", "rho"))
    for tag, lab, ck in (("L_obs", "A/C1  as shipped", "cos_A_C"),
                         ("L_rnd", "A/C1  uniform subset", "cos_A_C"),
                         ("L_ctl", "C1/C2 within-source null", "cos_C1_C2"),
                         ("L_ana", "A/C1  analytic, no operator", "cos_A_C")):
        ce, qq = [], []
        for r in rows:
            c, q = fit_c([r[tag][k] for k in keys], lams)
            ce.append(c); qq.append(q)
        ce = np.asarray(ce); d = np.array([r[ck] for r in rows], float)
        fin = np.isfinite(ce) & np.isfinite(d)
        print("      %-38s%10.4f%10.4f%10.4f%10.4f%9.3f"
              % (lab, ce[fin].mean(), np.median(ce[fin]), d[fin].mean(), np.median(d[fin]),
                 np.corrcoef(ce[fin], d[fin])[0, 1]))
        if tag == "L_obs":
            np.save(os.path.join(RES, "e_ceff_obs.npy"), ce)
            np.save(os.path.join(RES, "e_q_obs.npy"), np.asarray(qq))

    print("\n      READ: the ANALYTIC row is the shared-referent floor.  It contains ZERO")
    print("      information beyond the two bias vectors the direct cosine is computed from,")
    print("      so whatever agreement it shows is agreement BY CONSTRUCTION, not corroboration.")

    print("\n  (5) THE SPEC.  q* = 1/c_eff, per target, from the SHIPPED ladder.")
    ce = np.load(os.path.join(RES, "e_ceff_obs.npy"))
    q = np.load(os.path.join(RES, "e_q_obs.npy"))
    cd = np.array([r["cos_A_C"] for r in rows], float)
    for nm, v in (("c_eff (fitted)", ce), ("cos_A_C (direct)", cd)):
        f = np.isfinite(v)
        print("      %-20s mean %+.4f  median %+.4f  frac<=0 %.3f  p10 %+.3f  p90 %+.3f"
              % (nm, v[f].mean(), np.median(v[f]), float((v[f] <= 0).mean()),
                 np.percentile(v[f], 10), np.percentile(v[f], 90)))
    for nm, v in (("1/c_eff", ce), ("1/cos_A_C", cd)):
        f = np.isfinite(v) & (v > 0)
        qs = 1.0 / v[f]
        print("      q* = %-12s n_pos %3d  mean %12.1f  median %.4f  p25 %.3f  p75 %.3f"
              % (nm, f.sum(), qs.mean(), np.median(qs), np.percentile(qs, 25),
                 np.percentile(qs, 75)))
    print("      1/mean(c_eff) = %.4f      1/mean(cos_A_C) = %.4f" % (1 / ce[np.isfinite(ce)].mean(),
                                                                      1 / cd.mean()))
    print("      1/median(c_eff) = %.4f    1/median(cos_A_C) = %.4f"
          % (1 / np.median(ce[np.isfinite(ce)]), 1 / np.median(cd)))

    print("\n      DOES q < q* PREDICT AN ACTUAL PER-TARGET WIN AT THE SMALLEST MEASURED lambda?")
    M, _ = bows("L_obs")
    win = (M[:, 1] < M[:, 0])                       # m60_15 beats m75_0
    f = np.isfinite(ce) & (ce > 0)
    pred = np.zeros(len(rows), bool); pred[f] = q[f] < (1.0 / ce[f])
    tp = int((pred & win).sum()); fp = int((pred & ~win).sum())
    fn = int((~pred & win).sum()); tn = int((~pred & ~win).sum())
    print("        predicted-help & won %3d   predicted-help & lost %3d" % (tp, fp))
    print("        predicted-hurt & won %3d   predicted-hurt & lost %3d" % (fn, tn))
    base = win.mean()
    ph = tp / max(tp + fp, 1)
    print("        base win rate %.3f -> conditional on the bar %.3f   (lift %+.3f)"
          % (base, ph, ph - base))
    m, md, se, mde, flo, fhi, w = st(M[:, 1] - M[:, 0])
    print("        m60_15 vs m75_0 overall %+.4f SE %.4f MDE %.4f fold[%+.4f,%+.4f] %dW/%dL"
          % (m, se, mde, flo, fhi, w, len(rows) - w))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
