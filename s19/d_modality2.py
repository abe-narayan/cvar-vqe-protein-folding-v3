"""AGENT D, Sprint 19 -- P2/P3 hardened.

Three things the first pass left open, each of which could reverse a conclusion:

  1. P2a's matched-cell comparison had no CI and no target-level unit.  Added: a target-level
     bootstrap over the cell-matched contrast, plus a regression adjustment on log(sd) and
     separation, plus the same test under the PROMINENCE definition of a mode.
  2. P2b's "ORACLE nearest mode" is a MIN OVER k CANDIDATES.  s19/BRIEF.md section 9: no
     min-of-N ceiling without its min-of-N null.  Added: a control candidate set with the SAME
     |offsets| from the mean and randomised signs, minimised the same way; and a uniformly
     random mode pick, which is what a native-free scheme without an oracle would get.
  3. P3's decisive contrast: is the REAL residual's spatial coherence any different on
     multimodal pairs than on unimodal ones?  If not, multimodality cannot be the mechanism
     of a coherent error.

Run:  python -m s19.d_modality2
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s19 import d_modality as M            # noqa: E402

OUT = os.path.join(HERE, "results", "D_P2P3_hardened")
os.makedirs(OUT, exist_ok=True)
CENTRES, WIDTHS = M.CENTRES, M.WIDTHS


def boot(d, rng, B=4000):
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def target_rows():
    rows = []
    for t in I.targets():
        pdb = t["pdb"]
        z = np.load(os.path.join(ROOT, "s12", "cache", f"disto_{pdb}.npz"))
        p = np.asarray(z["prob"], float)
        i, j = np.asarray(z["i"], int), np.asarray(z["j"], int)
        exp = np.asarray(z["expected"], float)
        sd = np.asarray(z["sd"], float)
        u = np.load(os.path.join(ROOT, "s8", "generate_univ", f"{pdb}.npz"), allow_pickle=True)
        nat = np.asarray(u["nat_ca"], float)
        dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
        Md, Mp = M.dens_modes(p), M.prom_modes(p)
        rng = SD.stable_rng(pdb, "D2")

        def aim(Mask):
            pos = np.where(Mask, CENTRES[None], np.nan)
            off = pos - exp[:, None]                       # signed offsets of the modes
            dd = np.abs(pos - dtrue[:, None])
            has = Mask.any(1)
            bi = np.nanargmin(np.where(np.isnan(dd), np.inf, dd), axis=1)
            best = np.where(has, CENTRES[bi], exp)
            # uniformly random mode (native-free, no oracle)
            rnd = np.empty(len(exp))
            for r in range(len(exp)):
                idx = np.flatnonzero(Mask[r])
                rnd[r] = CENTRES[rng.choice(idx)] if len(idx) else exp[r]
            # MIN-OF-N NULL: same |offsets|, random signs, minimised the same way
            sgn = rng.choice([-1.0, 1.0], size=off.shape)
            cand = exp[:, None] + sgn * np.abs(off)
            cd = np.abs(cand - dtrue[:, None])
            cd = np.where(np.isnan(cd), np.inf, cd)
            nullbest = np.where(has, cand[np.arange(len(exp)), np.argmin(cd, 1)], exp)
            return best, rnd, nullbest, Mask.sum(1)

        b_d, r_d, n_d, k_d = aim(Md)
        b_p, r_p, n_p, k_p = aim(Mp)
        rows.append(dict(pdb=pdb, i=i, j=j, exp=exp, sd=sd, dtrue=dtrue, p=p,
                         sep=(j - i).astype(float), kd=k_d, kp=k_p,
                         best_d=b_d, rnd_d=r_d, null_d=n_d,
                         best_p=b_p, rnd_p=r_p, null_p=n_p,
                         fold=int(t["fold"])))
    return rows


def main():
    rows = target_rows()
    rng = SD.stable_rng("D", "P2P3h")
    cat = lambda k: np.concatenate([r[k] for r in rows])       # noqa: E731
    tid = np.concatenate([np.full(len(r["exp"]), t) for t, r in enumerate(rows)])
    exp, dtrue, sd, sep = cat("exp"), cat("dtrue"), cat("sd"), cat("sep")
    e_mean = np.abs(exp - dtrue)
    sepb = np.clip(np.searchsorted([2, 4, 6, 8, 11, 999], sep, side="right") - 1, 0, 4)

    print(f"\nn = {len(rows)} targets, {len(exp)} pairs\n")

    # ---------------------------------------------------------------- P2a hardened
    print("P2a -- |mean - true| on multimodal vs unimodal pairs, matched")
    for lab, k in [("density modes", "kd"), ("prominence modes", "kp")]:
        multi = cat(k) >= 2
        for nq, qlab in [(5, "sd quintile"), (10, "sd decile")]:
            qs = np.quantile(sd, np.arange(1, nq) / nq)
            sdb = np.searchsorted(qs, sd)
            per_t = []
            for t in range(len(rows)):
                num, den = 0.0, 0.0
                for a in range(5):
                    for b in range(nq):
                        m = (tid == t) & (sepb == a) & (sdb == b)
                        if (m & multi).sum() < 1 or (m & ~multi).sum() < 1:
                            continue
                        w = m.sum()
                        num += w * (e_mean[m & multi].mean() - e_mean[m & ~multi].mean())
                        den += w
                per_t.append(num / den if den else np.nan)
            m_, lo, hi = boot(np.asarray(per_t), rng)
            ok = np.isfinite(per_t)
            print(f"  {lab:<18}{qlab:<14} multi - uni = {m_:+.3f} [{lo:+.3f},{hi:+.3f}]  "
                  f"n_t={int(ok.sum())}  W/L "
                  f"{int((np.asarray(per_t)[ok] < 0).sum())}/"
                  f"{int((np.asarray(per_t)[ok] > 0).sum())}")
        # regression adjustment: residualise on sep dummies x log sd (linear+quadratic)
        X = np.column_stack([np.ones(len(sd)), np.log(sd), np.log(sd) ** 2]
                            + [(sepb == a).astype(float) for a in range(1, 5)]
                            + [((sepb == a) * np.log(sd)) for a in range(1, 5)])
        beta, *_ = np.linalg.lstsq(X, e_mean, rcond=None)
        res = e_mean - X @ beta
        per_t = np.asarray([res[(tid == t) & multi].mean() - res[(tid == t) & ~multi].mean()
                            if ((tid == t) & multi).any() and ((tid == t) & ~multi).any()
                            else np.nan for t in range(len(rows))])
        m_, lo, hi = boot(per_t, rng)
        print(f"  {lab:<18}{'regression-adj':<14} multi - uni = {m_:+.3f} [{lo:+.3f},{hi:+.3f}]")
        print(f"  {lab:<18}{'UNMATCHED':<14} multi - uni = "
              f"{e_mean[multi].mean() - e_mean[~multi].mean():+.3f}   "
              f"(mean sd {sd[multi].mean():.2f} vs {sd[~multi].mean():.2f}, "
              f"mean sep {sep[multi].mean():.1f} vs {sep[~multi].mean():.1f})")

    # ---------------------------------------------------------------- P2b min-of-N null
    print("\nP2b -- mode aim points on MULTIMODAL pairs, with the min-of-N null")
    for lab, kk, bk, rk, nk in [("density", "kd", "best_d", "rnd_d", "null_d"),
                                ("prominence", "kp", "best_p", "rnd_p", "null_p")]:
        multi = cat(kk) >= 2
        e = {"ORACLE_best": np.abs(cat(bk) - dtrue),
             "min-of-N NULL": np.abs(cat(nk) - dtrue),
             "random mode": np.abs(cat(rk) - dtrue)}
        print(f"  --- {lab} modes, {int(multi.sum())} pairs, mean k = "
              f"{cat(kk)[multi].mean():.2f}")
        for name, v in e.items():
            per_t = np.asarray([(v - e_mean)[(tid == t) & multi].mean()
                                if ((tid == t) & multi).any() else np.nan
                                for t in range(len(rows))])
            m_, lo, hi = boot(per_t, rng)
            print(f"    {name:<18} vs mean: {m_:+.3f} [{lo:+.3f},{hi:+.3f}]")
        per_o = np.asarray([(np.abs(cat(bk) - dtrue) - np.abs(cat(nk) - dtrue))[(tid == t) & multi].mean()
                            if ((tid == t) & multi).any() else np.nan for t in range(len(rows))])
        m_, lo, hi = boot(per_o, rng)
        print(f"    {'ORACLE - min-of-N null':<18}  {m_:+.3f} [{lo:+.3f},{hi:+.3f}]   "
              f"<-- the only line that prices MODE INFORMATION")

    # ---------------------------------------------------------------- P3 conditioned
    print("\nP3 -- residual sign coherence (adjacent minus non-adjacent), by modality class")
    agg = {"all": [], "multi": [], "uni": [], "cross": []}
    for r in rows:
        i, j = r["i"], r["j"]
        s = np.sign(r["dtrue"] - r["exp"])
        multi = r["kd"] >= 2
        share = (i[:, None] == i[None, :]) | (i[:, None] == j[None, :]) | \
                (j[:, None] == i[None, :]) | (j[:, None] == j[None, :])
        np.fill_diagonal(share, False)
        A = (s[:, None] * s[None, :]) > 0
        nz = (s[:, None] != 0) & (s[None, :] != 0)
        both_m = multi[:, None] & multi[None, :]
        both_u = (~multi)[:, None] & (~multi)[None, :]
        cross = ~(both_m | both_u)
        for lab, sel in [("all", np.ones_like(both_m)), ("multi", both_m),
                         ("uni", both_u), ("cross", cross)]:
            a = share & nz & sel
            b = (~share) & nz & sel
            np.fill_diagonal(b, False)
            agg[lab].append(float(A[a].mean()) - float(A[b].mean())
                            if a.sum() > 5 and b.sum() > 5 else np.nan)
    for lab in ("all", "multi", "uni", "cross"):
        m_, lo, hi = boot(np.asarray(agg[lab]), rng)
        n_ok = int(np.isfinite(agg[lab]).sum())
        print(f"  {lab:<8} adj - non-adj sign agreement = {m_:+.4f} [{lo:+.4f},{hi:+.4f}]  n_t={n_ok}")

    json.dump({"n": len(rows)}, open(os.path.join(OUT, "meta.json"), "w"))
    open(os.path.join(OUT, "COMPLETE"), "w").write(f"n={len(rows)}\n")


if __name__ == "__main__":
    main()
