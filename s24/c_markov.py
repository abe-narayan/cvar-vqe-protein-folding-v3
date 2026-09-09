"""s24/c_markov.py -- DOES CHAIN CORRELATION CHANGE THE VERDICT?  The arm my ladder was missing.

WHY THIS EXISTS.  Every arm in `c_ladder.py` samples residues INDEPENDENTLY.  That is the one
obvious objection to Lane C's recommendation not to train a network, because chain correlation is
exactly what an autoregressive torsion model, an HMM, or s19's entangled basin latent would supply
over a product distribution -- s19/BRIEF SS3 names it directly: "a whole region flipping between two
conformer families moves many pairs the same way at once."  If that is the missing ingredient, a
sequence model would be justified and Lane C's verdict would be wrong.

So this adds the cheapest honest version of it and asks whether it changes anything.

    T3_pool     per-residue 2-component von Mises, basin bits drawn INDEPENDENTLY   (the incumbent
                reference: 3 parameters per residue per component, zero training)
    T4_markov   the SAME von Mises basins, but the basin bits drawn from a first-order MARKOV chain
                whose 2x2 transition matrix is FITTED to the retrieved pool's own basin-label
                sequence.  Identical marginals by construction where the chain is stationary; the
                only thing added is CORRELATION BETWEEN NEIGHBOURING RESIDUES' BASINS.
    T5_shuffle  the FALSIFIER'S CONTROL.  The same Markov chain with its transition matrix replaced
                by the product of its own marginals -- i.e. correlation destroyed, marginals held
                EXACTLY fixed.  Any difference between T4 and T5 is chain correlation and nothing
                else.  Without this control, T4 vs T3 would confound correlation with the small
                marginal drift the chain fit introduces.

This is the "control must match the operator's space" rule applied to my own arm: the comparison
that prices correlation must hold the marginals fixed, and T3 does not do that.

Readout, functional, count, basis and nulls are IDENTICAL to `c_ladder.py`.  Nothing is trained.
ORACLE labelling on every quantity that reads the native.
"""
from __future__ import annotations

import json
import os
import sys
import time

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
from s19 import qb_lib as QB             # noqa: E402
from core import project as pj           # noqa: E402
from s24 import c_probe as CP            # noqa: E402
from s24 import c_ladder as CL           # noqa: E402

TOPM = 75
NSAMP = 2000
NMATCH = 500
MIX = [(75, 0), (60, 15), (50, 25), (38, 37), (25, 50), (0, 75)]
OUT = os.path.join(RES, "c_markov.json")


def basin_labels(PHI, PSI, mu):
    """Assign each pool member's residue to its nearest fitted basin, in the sin/cos embedding."""
    m, n = PHI.shape
    lab = np.zeros((m, n), int)
    for i in range(n):
        X = np.column_stack([np.cos(PHI[:, i]), np.sin(PHI[:, i]),
                             np.cos(PSI[:, i]), np.sin(PSI[:, i])])
        C = np.array([[np.cos(mu[i, c, 0]), np.sin(mu[i, c, 0]),
                       np.cos(mu[i, c, 1]), np.sin(mu[i, c, 1])] for c in (0, 1)])
        lab[:, i] = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(2), 1)
    return lab


def fit_markov(lab, eps=1.0):
    """First-order transition matrices T[i] = P(b_{i+1} | b_i), Laplace-smoothed.  Plus the
    marginals, so a correlation-destroying control can hold them exactly fixed."""
    m, n = lab.shape
    T = np.zeros((n - 1, 2, 2))
    for i in range(n - 1):
        for a in (0, 1):
            for b in (0, 1):
                T[i, a, b] = ((lab[:, i] == a) & (lab[:, i + 1] == b)).sum() + eps
            T[i, a] /= T[i, a].sum()
    marg = np.column_stack([(lab == 0).mean(0), (lab == 1).mean(0)])
    return T, marg


def draw_bits(T, marg, N, rng, correlated=True):
    """Sample basin bits.  `correlated=False` replaces every transition row by the NEXT residue's
    marginal, which destroys correlation while holding the marginals exactly fixed."""
    n = len(marg)
    bits = np.zeros((N, n), int)
    bits[:, 0] = (rng.random(N) > marg[0, 0]).astype(int)
    for i in range(n - 1):
        if correlated:
            p1 = T[i][bits[:, i], 1]
        else:
            p1 = np.full(N, marg[i + 1, 1])
        bits[:, i + 1] = (rng.random(N) < p1).astype(int)
    return bits


def run(nsamp=NSAMP, resume=True):
    tg = I.targets()
    rows = []
    if resume and os.path.exists(OUT):
        try:
            rows = json.load(open(OUT)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        pdb = t["pdb"]
        if pdb in done:
            continue
        n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n)
        idx = I.pool_idx(u); Wp = u["W"][idx]; rrp = u["rr"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, i, j)), float)
        top = np.argsort(scp, kind="stable")[:TOPM]
        A = Wp[top]; cA = CP._avg(A); eA = CP._bias(cA, nat)

        PHI = CP.wrap(u["PHI"][idx][top]); PSI = CP.wrap(u["PSI"][idx][top])
        rng0 = SD.stable_rng("c_markov", pdb)
        mu, kap, w = QB.fit_basins(PHI, PSI, rng0)
        lab = basin_labels(PHI, PSI, mu)
        T, marg = fit_markov(lab)
        # how much correlation is actually there?  mutual information of adjacent basin labels
        mi = []
        for q in range(n - 1):
            J = np.zeros((2, 2))
            for a in (0, 1):
                for b in (0, 1):
                    J[a, b] = ((lab[:, q] == a) & (lab[:, q + 1] == b)).mean()
            pa = J.sum(1); pb = J.sum(0); v = 0.0
            for a in (0, 1):
                for b in (0, 1):
                    if J[a, b] > 0 and pa[a] > 0 and pb[b] > 0:
                        v += J[a, b] * np.log(J[a, b] / (pa[a] * pb[b]))
            mi.append(v)

        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]),
             "incumbent": float(I.ca_rmsd(cA, nat)),
             "pool_best_ORACLE": float(rrp.min()),
             "adjacent_basin_MI_nats": float(np.mean(mi))}

        for tag, corr in (("T4_markov", True), ("T5_shuffle", False)):
            bits = draw_bits(T, marg, nsamp, rng0, correlated=corr)
            ph, ps = QB.draw_from_basins(bits, mu, kap, rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            sc = np.asarray(I.shipped_score(dg, I.pair_dists(CA, i, j)), float)
            o = np.argsort(sc, kind="stable")
            B = CA[o[:TOPM]]; cS = CP._avg(B)
            ru = rng0.choice(len(CA), TOPM, replace=False)
            rr = I.kabsch_rmsd_batch(CA, nat)                    # ORACLE
            mx = {}
            for a, b in MIX:
                s = np.concatenate([A[:a], B[:b]], 0) if (a and b) else (A[:a] if a else B[:b])
                mx["m%d_%d" % (a, b)] = float(I.ca_rmsd(CP._avg(s), nat))
            gi = SD.stable_rng("c_markov", pdb, tag).permutation(len(CA))[:NMATCH]
            Wu = np.concatenate([Wp, CA[gi]], 0)
            src = np.concatenate([np.zeros(len(Wp), int), np.ones(len(gi), int)])
            scu = np.asarray(I.shipped_score(dg, I.pair_dists(Wu, i, j)), float)
            ou = np.argsort(scu, kind="stable")[:TOPM]
            r[tag] = {
                "rmsd": float(I.ca_rmsd(cS, nat)),
                "rmsd_unif": float(I.ca_rmsd(CP._avg(CA[ru]), nat)),
                "cos_vs_incumbent": CP._cos(CP._bias(cS, nat), eA),
                "mix": mx,
                "union": float(I.ca_rmsd(CP._avg(Wu[ou]), nat)),
                "union_gen_frac": float(src[ou].mean()),
                "gen_best_ORACLE_matched": float(rr[gi].min()),
                "frac_better_than_pool_p75": float((sc[gi] < np.sort(scp)[TOPM - 1]).mean()),
                "mode": CL.mode_audit(ph, ps, CA),
                "geom": CL.geom_audit(CA, ph, ps),
            }
        rows.append(r)
        if len(rows) % 10 == 0 or len(rows) == len(tg):
            tmp = OUT + ".tmp"
            with open(tmp, "w") as fh:
                json.dump({"rows": rows, "complete": False, "n_expected": len(tg)}, fh)
            os.replace(tmp, OUT)
            print("  %d/%d  %.0fs" % (len(rows), len(tg), time.time() - t0), flush=True)
    tmp = OUT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "complete": len(rows) == len(tg), "n_expected": len(tg),
                   "nsamp": nsamp}, fh)
    os.replace(tmp, OUT)
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    F = sorted(set(fold.tolist()))
    rng = SD.stable_rng("c_markov", "rep")
    inc = np.array([r["incumbent"] for r in rows], float)

    def st(d):
        d = np.asarray(d, float); se = d.std(ddof=1) / np.sqrt(len(d))
        fs = [np.concatenate([d[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (d.mean(), se, 2.8016 * se, float(np.percentile(fs, 2.5)),
                float(np.percentile(fs, 97.5)), int((d < 0).sum()), float(d.max()))

    mi = np.array([r["adjacent_basin_MI_nats"] for r in rows], float)
    print("\nn=%d.  DOES CHAIN CORRELATION HELP?  incumbent %.4f (point cloud)" % (len(rows), inc.mean()))
    print("  adjacent-residue basin mutual information in the retrieved pool: %.4f nats"
          " (max possible %.4f)" % (mi.mean(), np.log(2)))
    print("\n  %-12s %9s %9s %9s %9s %9s" % ("arm", "RMSD", "cos_inc", "union", "genfrac",
                                             "ORACLEgb"))
    for k in ("T4_markov", "T5_shuffle"):
        f = lambda q: np.array([r[k][q] for r in rows], float)    # noqa: E731
        print("  %-12s %9.4f %9.4f %9.4f %9.4f %9.4f"
              % (k, f("rmsd").mean(), f("cos_vs_incumbent").mean(), f("union").mean(),
                 f("union_gen_frac").mean(), f("gen_best_ORACLE_matched").mean()))
    a = np.array([r["T4_markov"]["rmsd"] for r in rows], float)
    b = np.array([r["T5_shuffle"]["rmsd"] for r in rows], float)
    m, se, mde, lo, hi, wn, wt = st(a - b)
    v = ("BEATS" if hi < 0 and abs(m) > mde else "worse" if lo > 0 and abs(m) > mde else "NULL")
    print("\n  THE PRICE OF CHAIN CORRELATION, marginals held EXACTLY fixed:")
    print("    T4_markov - T5_shuffle  %+.4f  SE %.4f  MDE %.4f  fold[%+.4f,%+.4f] %3dW/%3dL"
          "  worst %+.2f  %s  (effect/MDE %.2fx)"
          % (m, se, mde, lo, hi, wn, len(rows) - wn, wt, v, abs(m) / mde if mde else float("nan")))
    for k in ("T4_markov", "T5_shuffle"):
        m, se, mde, lo, hi, wn, wt = st(np.array([r[k]["union"] for r in rows]) - inc)
        v = ("BEATS" if hi < 0 and abs(m) > mde else "worse" if lo > 0 and abs(m) > mde else "NULL")
        print("    %-12s UNION vs incumbent %+.4f  SE %.4f MDE %.4f fold[%+.4f,%+.4f] %3dW/%3dL %s"
              % (k, m, se, mde, lo, hi, wn, len(rows) - wn, v))
    print("\n  %-12s %9s %9s %9s" % ("arm", "uniq_frac", "dup_frac", "pw_RMSD"))
    for k in ("T4_markov", "T5_shuffle"):
        f = lambda q: np.array([r[k]["mode"][q] for r in rows], float)   # noqa: E731
        print("  %-12s %9.4f %9.4f %9.4f" % (k, f("unique_frac").mean(), f("dup_frac").mean(),
                                             f("pw_rmsd_mean").mean()))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
