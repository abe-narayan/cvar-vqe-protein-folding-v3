#!/usr/bin/env python
"""s29/s29_T_beta_unsat.py -- LANE T's OWN POST-MORTEM of the withdrawn corollary 2b (S29-L29).

THE QUESTION. Lane D's S29-L26 fired both limbs of my registered clause 2: beta is BELOW 1 (median
0.57 to 0.76 under four definitions) and the sign agreement is inside the coin-toss CI. My
post-mortem (`s29/THEORY.md` 2.3b) named two candidate causes and one discriminating test:

  (i) THE LINEARISATION (A2).  phi' = 2F - 1 saturates at +-1, and where it is saturated the
      per-pair coefficient is w*sign(a - b), not w*kappa*(a - b), so the expectation stops being a
      function of the regression slope at all.  If (i) is the cause, the sign law should WORK on
      the pairs where phi' is in its linear regime and FAIL on the saturated ones.
  (ii) (A4) IS AN IDEALISATION.  E[(a-b)(a-n)] = var(a) - cov(a,b) - cov(a,n) + cov(b,n) and I set
      the last two to zero; if they dominate, the sign law fails on BOTH subsets alike and the
      informative reading is cov(a,n) > cov(b,n) -- the POOL tracking the native better than the
      POSTERIOR's median map does.

THE DESIGN.  This file imports lane D's OWN implementation (`s29_D_theory_check`) for every piece
that both entries share -- the production cloud, the posterior's median map, the tau definitions,
the shipped weights, kappa -- and changes exactly one thing: the PAIR SUBSET.  Using the auditor's
code rather than my own is deliberate; a post-mortem written by the author of the failed claim
should not re-derive a favourable variant of the quantity that failed.  The saturated subset is
carried beside the unsaturated one as the matched control, so the test is a CONTRAST, not a
selected half.

ORACLE: cos_DIS comes from `s27/results/s28_A2_cosine_rows.jsonl` and is a diagnostic, as in
S29-L7 and S29-L26.  Nothing here is deployable and nothing selects a parameter.

    python s29/s29_T_beta_unsat.py [--limit N] [--thresh 0.5]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s24 import stats_lib as ST                      # noqa: E402
from s29 import s29_D_theory_check as D              # noqa: E402

OUT = os.path.join(HERE, "results", "s29_T_beta_unsat.json")
TAUS = ("LIB75", "UNIV", "POOL")
THRESHES = (0.3, 0.5, 0.7)


def rows_for(pdbs, taus=TAUS):
    from s12 import instrument as I
    cosrows = {r["pdb"]: r for r in (json.loads(l) for l in
               open(os.path.join(ROOT, "s27", "results", "s28_A2_cosine_rows.jsonl"),
                    encoding="utf-8") if l.strip())}
    out = []
    for q, pdb in enumerate(pdbs):
        cand, C0, dg, i, j, D0, u, sur = D.target_pieces(pdb)
        uu = I.load_univ(pdb)
        m_grid, _ = D.median_map(dg)
        coeff, w_ship = D.cdf_coeff_at(dg, D0)
        phi = coeff / np.maximum(w_ship, 1e-30)              # = 2 F(d) - 1, the saturating factor
        rec = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), P=int(len(i)),
                   cos_DIS=float(cosrows[pdb]["cos"]["DIS"]) if pdb in cosrows else float("nan"),
                   fail18=bool(pdb in I.FAIL18),
                   sat_share={str(t): float((np.abs(phi) >= t).mean()) for t in THRESHES},
                   mean_abs_phi=float(np.abs(phi).mean()), beta={})
        for which in taus:
            tau, _meta = D.typical_map(pdb, uu, i, j, which)
            a = D0 - tau
            b = m_grid - tau
            for t in THRESHES:
                un = np.abs(phi) < t                          # the LINEAR regime of the risk
                sa = ~un                                       # the matched control
                for tag, msk in (("unsat", un), ("sat", sa)):
                    key = f"{which}_{tag}{t}"
                    if msk.sum() >= 5 and float(a[msk] @ a[msk]) > 1e-30:
                        rec["beta"][key] = float((a[msk] @ b[msk]) / (a[msk] @ a[msk]))
                        rec["beta"][key + "_npairs"] = int(msk.sum())
                    else:
                        rec["beta"][key] = float("nan")
                        rec["beta"][key + "_npairs"] = int(msk.sum())
            rec["beta"][which + "_all"] = float((a @ b) / max(a @ a, 1e-30))
        out.append(rec)
        if (q + 1) % 30 == 0:
            print(f"  [{q+1}/{len(pdbs)}]", flush=True)
    return out


def summarise(rows, taus=TAUS):
    cos = np.array([r["cos_DIS"] for r in rows])
    folds = ST.pinned_folds([r["pdb"] for r in rows])
    n = len(rows)
    se = float(np.sqrt(0.25 / n))
    ci = [0.5 - 1.96 * se, 0.5 + 1.96 * se]
    out = dict(check="lane T's own post-mortem of the withdrawn corollary 2b (S29-L29)",
               n=n, oracle=True, taus=list(taus), threshes=list(THRESHES),
               coin_toss_ci=ci,
               sat_share_median={str(t): float(np.median([r["sat_share"][str(t)] for r in rows]))
                                 for t in THRESHES},
               sat_share_fail18={str(t): float(np.median([r["sat_share"][str(t)] for r in rows
                                                          if r["fail18"]])) for t in THRESHES},
               sat_share_other={str(t): float(np.median([r["sat_share"][str(t)] for r in rows
                                                         if not r["fail18"]])) for t in THRESHES},
               mean_abs_phi_median=float(np.median([r["mean_abs_phi"] for r in rows])))
    cells = {}
    for which in taus:
        for t in THRESHES:
            for tag in ("unsat", "sat"):
                key = f"{which}_{tag}{t}"
                bet = np.array([r["beta"][key] for r in rows], float)
                ok = np.isfinite(bet) & np.isfinite(cos)
                if ok.sum() < 30:
                    continue
                agree = (np.sign(cos[ok]) == -np.sign(bet[ok] - 1.0))
                k = int(agree.sum())
                m = int(ok.sum())
                cells[key] = dict(
                    n_targets=m,
                    npairs_median=float(np.median([r["beta"][key + "_npairs"] for r in rows])),
                    beta_median=float(np.median(bet[ok])),
                    frac_above_1=float((bet[ok] > 1).mean()),
                    sign_agreement=k / m,
                    clears_coin_toss=bool(k / m > ci[1]),
                    clears_seventy=bool(k / m >= 0.70),
                    fold_ci_beta=ST.compare(bet[ok], np.zeros(m), folds[ok],
                                            label="beta " + key,
                                            seed_parts=("s29Tbu",))["ci95_fold"])
    out["cells"] = cells
    # the verdict the entry will quote: the primary cell is the registered threshold 0.5
    prim = {}
    for which in taus:
        u = cells.get(f"{which}_unsat0.5")
        sa = cells.get(f"{which}_sat0.5")
        if u and sa:
            prim[which] = dict(
                unsat_agreement=u["sign_agreement"], sat_agreement=sa["sign_agreement"],
                unsat_beta_median=u["beta_median"], sat_beta_median=sa["beta_median"],
                unsat_clears=u["clears_coin_toss"], sat_clears=sa["clears_coin_toss"],
                contrast=u["sign_agreement"] - sa["sign_agreement"])
    out["primary_thresh_0.5"] = prim
    out["verdict"] = ("NOT DECIDED (too few targets for any cell)" if not prim else
                      "(i) SATURATION is the cause" if all(v["unsat_clears"] for v in prim.values())
                      else ("(ii) (A4): the law fails on both subsets alike"
                            if not any(v["unsat_clears"] for v in prim.values())
                            else "MIXED: clears on some tau definitions and not others"))
    return out


# ============================================================ the FULL second-order expression
def full_terms(pdbs, taus=TAUS):
    """The four terms of E[<g,r>] without (A4)'s truncation, and whether their sum predicts the
    measured cosine's sign.

    E[<g,r>] = -( var(a) - cov(a,b) - cov(a,n) + cov(b,n) ) in the w*kappa metric, with
    n = d(t) - tau the NATIVE's deviation from typical (ORACLE, diagnostic only).  My registered
    corollary 2b kept the first two terms and set the last two to zero; this measures all four, so
    the question "was the bookkeeping wrong or was only (A4) wrong" is answered rather than argued.
    """
    from s12 import instrument as I
    cosrows = {r["pdb"]: r for r in (json.loads(l) for l in
               open(os.path.join(ROOT, "s27", "results", "s28_A2_cosine_rows.jsonl"),
                    encoding="utf-8") if l.strip())}
    rows = []
    for q, pdb in enumerate(pdbs):
        cand, C0, dg, i, j, D0, u, sur = D.target_pieces(pdb)
        if cand.nat_ca is None or not np.isfinite(np.asarray(cand.nat_ca, float)).all():
            continue
        uu = I.load_univ(pdb)
        m_grid, _ = D.median_map(dg)
        coeff, w_ship = D.cdf_coeff_at(dg, D0)
        phi = coeff / np.maximum(w_ship, 1e-30)
        prob = np.asarray(dg["prob"], float); cen = np.asarray(dg["centres"], float)
        kb = np.clip(np.searchsorted(cen, m_grid) - 1, 0, prob.shape[1] - 1)
        kappa = 2.0 * prob[np.arange(prob.shape[0]), kb]
        ww = np.maximum(w_ship * kappa, 0.0)
        Dt = I.pair_dists(np.asarray(cand.nat_ca, float)[None], i, j)[0]      # ORACLE
        rec = dict(pdb=pdb, fail18=bool(pdb in I.FAIL18),
                   cos_DIS=float(cosrows[pdb]["cos"]["DIS"]) if pdb in cosrows else float("nan"),
                   terms={})
        for which in taus:
            tau, _ = D.typical_map(pdb, uu, i, j, which)
            a, b, nn = D0 - tau, m_grid - tau, Dt - tau
            for tag, msk in (("all", np.ones(len(a), bool)), ("unsat", np.abs(phi) < 0.5)):
                if msk.sum() < 5:
                    continue
                W = ww[msk]
                va = float(W @ (a[msk] * a[msk])); cab = float(W @ (a[msk] * b[msk]))
                can = float(W @ (a[msk] * nn[msk])); cbn = float(W @ (b[msk] * nn[msk]))
                kept = va - cab                     # what corollary 2b used
                drop = -can + cbn                   # what (A4) set to zero
                rec["terms"][f"{which}_{tag}"] = dict(
                    var_a=va, cov_ab=cab, cov_an=can, cov_bn=cbn,
                    kept=kept, dropped=drop, full=kept + drop,
                    ratio_drop_over_kept=float(abs(drop) / max(abs(kept), 1e-30)))
        rows.append(rec)
        if (q + 1) % 30 == 0:
            print(f"  [{q+1}/{len(pdbs)}]", flush=True)
    cos = np.array([r["cos_DIS"] for r in rows])
    out = dict(check="the full second-order expression (no (A4) truncation)", n=len(rows),
               oracle=True, cells={})
    nn_ = len(rows); se = float(np.sqrt(0.25 / max(nn_, 1)))
    ci = [0.5 - 1.96 * se, 0.5 + 1.96 * se]
    out["coin_toss_ci"] = ci
    for key in sorted({k for r in rows for k in r["terms"]}):
        kept = np.array([r["terms"][key]["kept"] for r in rows if key in r["terms"]])
        full = np.array([r["terms"][key]["full"] for r in rows if key in r["terms"]])
        drop = np.array([r["terms"][key]["dropped"] for r in rows if key in r["terms"]])
        cs = np.array([r["cos_DIS"] for r in rows if key in r["terms"]])
        ok = np.isfinite(cs)
        # E[<g,r>] = -(expression); cos has the sign of -<g,r>, i.e. the sign of the expression
        out["cells"][key] = dict(
            n=int(ok.sum()),
            agree_kept_only=float((np.sign(cs[ok]) == np.sign(kept[ok])).mean()),
            agree_full=float((np.sign(cs[ok]) == np.sign(full[ok])).mean()),
            median_ratio_drop_over_kept=float(np.median(np.abs(drop) / np.maximum(np.abs(kept), 1e-30))),
            frac_dropped_dominates=float((np.abs(drop) > np.abs(kept)).mean()),
            median_cov_an=float(np.median([r["terms"][key]["cov_an"] for r in rows if key in r["terms"]])),
            median_cov_bn=float(np.median([r["terms"][key]["cov_bn"] for r in rows if key in r["terms"]])),
            frac_cov_an_gt_cov_bn=float(np.mean([r["terms"][key]["cov_an"] > r["terms"][key]["cov_bn"]
                                                 for r in rows if key in r["terms"]])))
    out["rows"] = rows
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args(argv)
    from s12 import instrument as I
    pdbs = [t["pdb"] for t in I.targets()]
    if a.limit:
        pdbs = pdbs[:a.limit]
    if a.full:
        outf = full_terms(pdbs)
        ST.save_atomic(os.path.join(HERE, "results", "s29_T_beta_full.json"), outf)
        for k in sorted(outf["cells"]):
            c = outf["cells"][k]
            print(f"  {k:16s} n={c['n']:3d} agree(kept) {c['agree_kept_only']:.3f} "
                  f"agree(FULL) {c['agree_full']:.3f} |drop|/|kept| {c['median_ratio_drop_over_kept']:.2f} "
                  f"drop dominates {c['frac_dropped_dominates']:.3f} cov_an>cov_bn {c['frac_cov_an_gt_cov_bn']:.3f}")
        print("coin CI", [round(x,3) for x in outf["coin_toss_ci"]])
        return 0
    rows = rows_for(pdbs)
    out = summarise(rows)
    out["rows"] = rows
    ST.save_atomic(OUT, out)
    print(json.dumps({k: v for k, v in out.items() if k not in ("rows", "cells")}, indent=1))
    for k in sorted(out["cells"]):
        c = out["cells"][k]
        print(f"  {k:22s} pairs {c['npairs_median']:6.1f}  beta {c['beta_median']:+.3f}  "
              f"agree {c['sign_agreement']:.3f}  clears {c['clears_coin_toss']}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
