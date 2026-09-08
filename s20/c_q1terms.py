"""s20/c_q1terms.py -- Q1, the component-level answer: WHICH PART of Legacy disagrees with AMBER?

Declared EXTENSION, written after `s20/results/c_q1.json` was read and labelled as such
(`PREREG_C.md` allows extensions; each must say what it was written after).  It was written after
section 1 of `c_q1_report.txt` showed `rho(E_Legacy, E_AMBER)` is **negative** (median -0.095),
and it asks the obvious next question, which needs no new AMBER compute.

THE DECOMPOSITION IS AN IDENTITY, NOT A DISCOVERY (label: EXACT).  With `E_Legacy = sum_t w_t c_t`
and `zA` the within-target standardised AMBER single point,

    corr(E_Legacy, E_AMBER) = sum_t  w_t * cov(c_t, zA) / sd(E_Legacy)

so the Pearson correlation between the two potentials decomposes ADDITIVELY over Legacy's eleven
components.  What is scientific is which components carry it, not that it decomposes.

The Spearman columns beside it are NOT additive and are reported as separate diagnostics.

    python -m s20.c_q1terms
"""
from __future__ import annotations

import os
import sys
import json

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
from s12 import instrument as I          # noqa: E402
from s16 import energy_lib as EL         # noqa: E402
from s18 import phys_lib as PL           # noqa: E402
from s20 import c_q1 as Q1               # noqa: E402


# ---------------------------------------------------------------------------------------
# A SELF-CAUGHT LABELLING DEFECT.  `s12.instrument.paired`'s `ci95` is a PLAIN i.i.d.
# target-level bootstrap; its `folds` argument only adds a per-fold mean breakdown and does
# NOT cluster the resample.  Sprint 19's numbers were made with `s18.phys_lib.paired`, which
# returns BOTH an i.i.d. `ci` and a FOLD-CLUSTERED `ci_fold`, and "fold-aware" in this
# programme means the latter (BRIEF section 8).  Every CI below is therefore produced by
# `PL.paired` and the FOLD-CLUSTERED interval is the one quoted, with the i.i.d. one printed
# beside it wherever both matter.  Caught by reading `s12/instrument.py:188` rather than
# trusting the parameter name.
def PP(a, b, folds=None):
    st = PL.paired(np.asarray(a, float), np.asarray(b, float), folds=folds)
    return {"mean_diff": st["mean"], "ci95": st.get("ci_fold", st["ci"]),
            "ci_iid": st["ci"], "median_diff": st["median"],
            "n_better": st["W"], "n_worse": st["L"],
            "mean_a": st["mean_a"], "mean_b": st["mean_b"], "n": st["n"]}


def spearman(a, b):
    return Q1._spearman(a, b)


def run(out="c_q1_terms.json"):
    down = Q1.load_down()
    q1 = json.load(open(os.path.join(RESULTS, "c_q1.json")))
    if not q1.get("complete"):
        raise SystemExit("c_q1.json is not complete")
    w = EL.legacy_weight_vector()
    terms = list(EL.LEG_TERMS)
    rows = []
    for t in I.targets():
        pdb = t["pdb"]
        dr = down[pdb]
        comp = {k: np.asarray(v, float) for k, v in dr["leg_terms"].items()}
        e_leg = np.asarray(EL.legacy_total_from(comp), float)
        e_amb = np.asarray(dr["s_amber_sp"], float)
        d = np.asarray(dr["d"], float)                     # ORACLE (window RMSD)
        sa = e_amb.std()
        zA = (e_amb - e_amb.mean()) / sa if sa > 1e-12 else np.zeros_like(e_amb)
        sl = e_leg.std()
        contrib = {}
        for ti, tn in enumerate(terms):
            c = comp[tn]
            contrib[tn] = float(w[ti] * np.cov(c, zA, bias=True)[0, 1] / sl) if sl > 1e-12 \
                else float("nan")
        rows.append({
            "pdb": pdb, "fold": int(t["fold"]), "n": int(t["n"]),
            "corr_leg_amb": float(np.corrcoef(e_leg, e_amb)[0, 1]),
            "contrib": contrib,
            "sd_term_w": {tn: float(abs(w[ti]) * comp[tn].std()) for ti, tn in enumerate(terms)},
            "rho_term_amb": {tn: spearman(comp[tn], e_amb) for tn in terms},
            "rho_term_d": {tn: spearman(comp[tn], d) for tn in terms},
        })
    obj = {"rows": rows, "terms": terms,
           "weights": {tn: float(w[i]) for i, tn in enumerate(terms)},
           "config": {"source": "s18/results/down.json", "note": "declared extension"},
           "n_rows": len(rows), "n_expected": 126, "complete": len(rows) >= 126}
    obj["cfg_hash"] = PL.cfg_hash(obj["config"])
    with open(os.path.join(RESULTS, out), "w") as fh:
        json.dump(obj, fh)
    report(obj)
    return obj


def report(o=None, out="c_q1_terms.json"):
    if o is None:
        o = json.load(open(os.path.join(RESULTS, out)))
    rows, terms = o["rows"], o["terms"]
    folds = np.array([r["fold"] for r in rows], int)
    tot = np.array([r["corr_leg_amb"] for r in rows], float)
    print("=" * 104)
    print("Q1b  WHICH PART OF LEGACY DISAGREES WITH AMBER?   n =", len(rows))
    print("     EXACT identity: corr(E_Legacy, E_AMBER) = sum_t w_t cov(c_t, zA)/sd(E_Legacy)")
    print("=" * 104)
    st = PP(tot, np.zeros_like(tot), folds=folds)
    print(f"\nPearson corr(E_Legacy, E_AMBER):  mean {tot.mean():+.4f} "
          f"[{st['ci95'][0]:+.4f}, {st['ci95'][1]:+.4f}]  median {np.median(tot):+.4f}  "
          f"{(tot < 0).sum()}/{len(tot)} targets NEGATIVE")
    C = {t: np.array([r["contrib"][t] for r in rows], float) for t in terms}
    print(f"\n{'component':<18s}{'weight':>9s}{'w*sd(c)':>10s}{'contribution':>14s}"
          f"{'CI95':>24s}{'share':>8s}{'rho(c,E_A)':>12s}{'rho(c,ORACLE d)':>17s}")
    order = sorted(terms, key=lambda t: C[t].mean())
    ssum = sum(C[t].mean() for t in terms)
    for t in order:
        s2 = PP(C[t], np.zeros_like(C[t]), folds=folds)
        ra = np.array([r["rho_term_amb"][t] for r in rows], float)
        rd = np.array([r["rho_term_d"][t] for r in rows], float)
        sw = np.array([r["sd_term_w"][t] for r in rows], float)
        print(f"{t:<18s}{o['weights'][t]:9.3f}{sw.mean():10.4f}{C[t].mean():+14.4f}"
              f"  [{s2['ci95'][0]:+9.4f},{s2['ci95'][1]:+9.4f}]"
              f"{C[t].mean()/ssum if ssum else float('nan'):8.2f}"
              f"{np.nanmean(ra):12.3f}{np.nanmean(rd):17.3f}")
    print(f"{'SUM':<18s}{'':>9s}{'':>10s}{ssum:+14.4f}   (identity check vs the measured "
          f"corr {tot.mean():+.4f}: residual {abs(ssum - tot.mean()):.2e})")
    print("\n    rho(c, E_AMBER) and rho(c, ORACLE d) are Spearman and are NOT additive; they "
          "are separate diagnostics.")
    print("=" * 104)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--report":
        report()
    else:
        run()
