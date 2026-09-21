#!/usr/bin/env python
"""s31/s31_B2_sweep_parity.py -- lane B: IS THE 32-COST METER SWEEP'S PREFERENCE CONTRAST
MEASURING CHIRALITY DETECTION RATHER THAN STRUCTURAL SKILL?

The hypothesis, raised by lane B's F3 pre-check and accepted by the coordinator as theirs to
correct: the two chain costs with the largest preference contrast in the S30 32-cost sweep
(`LEG_steric` +0.2515 at 2.88x MDE, `LEG_torsion` +0.2212 at 2.14x) are EXACTLY the two with the
largest ODD (chirality-carrying) variance share (13.02 and 0.301, s31/results/s31_B1_achirality.json).

MECHANISM, construction-level not correlational: `core/geometry.build_backbone_batch` places CB as
    CB = -0.58273431 * cross(b, d) + 0.56802827 * b - 0.54067466 * d + CA
`cross` is a PSEUDOVECTOR, so every CB-dependent Legacy term inherits chirality, and the
Ramachandran basins `core/energy._RAMA_BASINS` sit at negative phi with alpha_L at depth 0.40
against alpha_R's 1.00.

THE TEST. The published statistic is
    contrast(f) = pref(circ_best vs PROD) - mean_r pref(RAND_SIGNED[r] vs PROD)
with pref(X vs PROD) = share of the 126 targets on which f(X) < f(PROD). Decompose f under the
EXACT point reflection R of the rung -- coordinates negated about their centroid, torsions negated:
    f_even = 1/2 [ f(x) + f(Rx) ]     reflection-INVARIANT  -> a distance-map function by G1
    f_odd  = 1/2 [ f(x) - f(Rx) ]     reflection-ANTI-invariant
and recompute the SAME contrast on each part.

  If contrast(f_even) ~ 0 and contrast(f_odd) ~ contrast(f), the row is chirality detection.
  If contrast(f_even) ~ contrast(f), the row stands and the coincidence is a coincidence.

EVERY NUMBER HERE IS ORACLE: `circ_best` is built with the native. The contrast is a diagnostic of
the METER, not a deployable quantity. Basis: BUILT CHAIN rungs (`chain_ca/phi/psi`), matching the
published `_chain` sweep files.

USAGE
    python s31/s31_B2_sweep_parity.py [--limit N] [--draws 8]
"""
from __future__ import annotations

import argparse
import io
import json
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_C2_recog_audit as C2   # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s31_B2_sweep_parity.json")
LADDER = os.path.join(ROOT, "s30", "results", "s30_D_ladder_structs")

CHAIN_COSTS = list(C2.CHAIN_SCORERS)
CA_ASSERT = ["DIS", "CAGEO", "CONTACT"]     # must come out EXACTLY even -- a self-check on the mirror


def mirror_ca(C):
    """EXACT point reflection of a (m, n, 3) stack about each structure's own centroid."""
    C = np.asarray(C, float)
    return 2.0 * C.mean(axis=-2, keepdims=True) - C


def load_rungs(pdb, draws):
    """(names, CA (m,n,3), PHI (m,n), PSI (m,n)) for PROD, circ_best and the RAND_SIGNED draws."""
    with np.load(os.path.join(LADDER, f"{pdb}.npz"), allow_pickle=False) as z:
        keys = set(z.files)
        names, CA, PHI, PSI = [], [], [], []

        def add(nm, pre=""):
            ca, ph, ps = f"{pre}chain_ca_{nm}", f"{pre}chain_phi_{nm}", f"{pre}chain_psi_{nm}"
            if ca not in keys:
                return False
            names.append(nm); CA.append(np.asarray(z[ca], float))
            PHI.append(np.asarray(z[ph], float)); PSI.append(np.asarray(z[ps], float))
            return True

        if not add("PROD") or not add("circ_best"):
            return None
        got = 0
        for s in range(draws):
            if add("RAND_SIGNED_%d" % s, "" if s == 0 else "x"):
                got += 1
        if got == 0:
            return None
        # extra matched-distance controls that are NOT signed affine combinations of the pool
        add("GAUSS_MATCHED_0")
        add("GAUSS_0.3_0")
        add("sub0")
    return names, np.stack(CA), np.stack(PHI), np.stack(PSI)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--draws", type=int, default=8)
    a = ap.parse_args()

    from s27 import s28_A_amp as A

    tg = I.targets()
    if a.limit:
        tg = tg[: a.limit]

    # per-target indicator arrays, one per cost x part
    parts = ("tot", "even", "odd")
    acc = {(c, p): [] for c in CHAIN_COSTS + CA_ASSERT for p in parts}
    # the contrast's two TERMS kept apart, plus two controls that are not signed affine
    # combinations of the pool. A positive contrast can come from the cost liking circ_best OR
    # from it hating the control; only the split says which.
    terms = ("pref_best", "pref_rand", "pref_gaussm", "pref_gauss03", "pref_sub0")
    tacc = {(c, p, q): [] for c in CHAIN_COSTS + CA_ASSERT for p in parts for q in terms}
    mirror_selfcheck = {c: 0.0 for c in CA_ASSERT}
    pdbs, folds, skipped = [], [], []
    t0 = time.time()

    for ti, t in enumerate(tg):
        pdb = t["pdb"]
        rl = load_rungs(pdb, a.draws)
        if rl is None:
            skipped.append(pdb)
            continue
        names, CA, PHI, PSI = rl
        try:
            cand, dis, top, dg = A.load_pool(pdb)
        except Exception as exc:                                   # noqa: BLE001
            skipped.append(pdb + ":" + repr(exc)[:60])
            continue
        univ = I.load_univ(pdb)
        cp = SimpleNamespace(pdb=pdb, seq=t["seq"], n=int(t["n"]), fold=int(t["fold"]),
                             k=int(cand.k), W=np.asarray(cand.W, float),
                             PHI=np.asarray(cand.PHI, float), PSI=np.asarray(cand.PSI, float))

        sc_p = C2.chain_scores(cp, CA, PHI, PSI, univ, dg)
        sc_m = C2.chain_scores(cp, mirror_ca(CA), -PHI, -PSI, univ, dg)
        ca_p = C2.ca_scores(cp, CA, univ, dg)
        ca_m = C2.ca_scores(cp, mirror_ca(CA), univ, dg)
        for c in CA_ASSERT:
            d = float(np.abs(np.asarray(ca_p[c], float) - np.asarray(ca_m[c], float)).max())
            rng = float(np.ptp(np.asarray(ca_p[c], float))) + 1e-12
            mirror_selfcheck[c] = max(mirror_selfcheck[c], d / rng)

        i_best = names.index("circ_best")
        i_prod = names.index("PROD")
        i_rs = [k for k, nm in enumerate(names) if nm.startswith("RAND_SIGNED")]
        i_ex = {"pref_gaussm": "GAUSS_MATCHED_0", "pref_gauss03": "GAUSS_0.3_0",
                "pref_sub0": "sub0"}

        def pref(x, k, kp):
            """The meter's own rule (s29_D_cost_audit.py:454): TIES COUNT 0.5."""
            return 1.0 if x[k] < x[kp] else (0.5 if x[k] == x[kp] else 0.0)

        for c in CHAIN_COSTS + CA_ASSERT:
            v = np.asarray((sc_p if c in CHAIN_COSTS else ca_p)[c], float)
            w = np.asarray((sc_m if c in CHAIN_COSTS else ca_m)[c], float)
            f = {"tot": v, "even": 0.5 * (v + w), "odd": 0.5 * (v - w)}
            for p in parts:
                x = f[p]
                if not np.isfinite(x).all():
                    acc[(c, p)].append(np.nan)
                    for q in terms:
                        tacc[(c, p, q)].append(np.nan)
                    continue
                pb = pref(x, i_best, i_prod)
                pr = float(np.mean([pref(x, k, i_prod) for k in i_rs]))
                acc[(c, p)].append(pb - pr)
                tacc[(c, p, "pref_best")].append(pb)
                tacc[(c, p, "pref_rand")].append(pr)
                for q, nm in i_ex.items():
                    tacc[(c, p, q)].append(pref(x, names.index(nm), i_prod)
                                           if nm in names else np.nan)
        pdbs.append(pdb); folds.append(int(t["fold"]))
        if (ti + 1) % 20 == 0:
            print("  %3d/%d  %.0fs" % (ti + 1, len(tg), time.time() - t0), flush=True)

    folds = np.asarray(folds, int)
    n = len(pdbs)
    print("\nn targets %d   skipped %d   %.0fs" % (n, len(skipped), time.time() - t0), flush=True)
    print("MIRROR SELF-CHECK (CA costs are functions of the distance map -> must be EXACTLY even):")
    for c, v in mirror_selfcheck.items():
        print("  %-10s max|f(x)-f(Rx)| / range(f) = %.3e" % (c, v), flush=True)

    rows, ncomp = [], 0
    for c in CHAIN_COSTS + CA_ASSERT:
        row = {"cost": c, "n": n}
        for p in parts:
            d = np.asarray(acc[(c, p)], float)
            if not np.isfinite(d).all():
                row[p] = {"effect": float("nan"), "note": "non-finite"}
                continue
            # contrast > 0 is the cost DOING WELL, so pass -d as arm a (compare is lower-is-better)
            o = ST.compare(-d, np.zeros(n), folds=folds, names=pdbs,
                           label="%s_%s" % (c, p), seed_parts=("s31B", "sweepparity"))
            ncomp += 1
            row[p] = {"effect": float(d.mean()), "median": float(np.median(d)),
                      "se": o["se"], "mde": o["mde"],
                      "x_mde": float(d.mean() / o["mde"]) if o["mde"] > 0 else float("nan"),
                      "ci95_fold": [-o["ci95_fold"][1], -o["ci95_fold"][0]],
                      "folds_same_sign": o["folds_same_sign"],
                      "per_fold": {k: -v for k, v in o["per_fold"].items()}}
        t_, e_, o_ = row["tot"]["effect"], row["even"]["effect"], row["odd"]["effect"]
        row["odd_share_of_contrast"] = float(o_ / t_) if abs(t_) > 1e-9 else float("nan")
        row["terms"] = {}
        for p in parts:
            row["terms"][p] = {}
            for q in terms:
                d = np.asarray(tacc[(c, p, q)], float)
                m = np.isfinite(d)
                row["terms"][p][q] = float(d[m].mean()) if m.sum() >= 20 else None
            pb = row["terms"][p].get("pref_best")
            for q, lab in (("pref_gaussm", "contrast_vs_GAUSS_MATCHED"),
                           ("pref_gauss03", "contrast_vs_GAUSS_0.3")):
                pv = row["terms"][p].get(q)
                row["terms"][p][lab] = (pb - pv) if (pb is not None and pv is not None) else None
        rows.append(row)

    rows.sort(key=lambda r: -abs(r["tot"]["effect"]))
    print("\n%-22s %9s %9s %9s %8s %8s %8s   %s" %
          ("cost", "tot", "even", "odd", "xMDEtot", "xMDEev", "xMDEodd", "odd/tot"))
    for r in rows:
        print("%-22s %+9.4f %+9.4f %+9.4f %8.2f %8.2f %8.2f   %s" % (
            r["cost"], r["tot"]["effect"], r["even"]["effect"], r["odd"]["effect"],
            r["tot"]["x_mde"], r["even"]["x_mde"], r["odd"]["x_mde"],
            ("%+.2f" % r["odd_share_of_contrast"]) if np.isfinite(r["odd_share_of_contrast"]) else "n/a"))

    print("\nTHE CONTRAST SPLIT INTO ITS TWO TERMS (total cost), plus two controls that are NOT\n"
          "signed affine combinations of the pool. A cost that RECOGNISES nativeness must have a\n"
          "high pref(circ_best vs PROD); a cost that merely DETECTS a broken control has a low one.")
    print("%-22s %10s %10s %10s %10s %10s" %
          ("cost", "pref_best", "pref_rand", "pref_gaussM", "cts_vs_RS", "cts_vs_GM"))
    for r in rows:
        tm = r["terms"]["tot"]
        def g(k):
            return ("%+10.4f" % tm[k]) if tm.get(k) is not None else "       n/a"
        print("%-22s %s %s %s %s %s" % (r["cost"], g("pref_best"), g("pref_rand"),
                                        g("pref_gaussm"), ("%+10.4f" % r["tot"]["effect"]),
                                        g("contrast_vs_GAUSS_MATCHED")))

    out = dict(n=n, skipped=skipped, draws=a.draws, rows=rows,
               mirror_selfcheck=mirror_selfcheck, comparisons_emitted=ncomp,
               basis="BUILT CHAIN rungs (chain_ca/phi/psi), matching the published _chain sweep",
               oracle="ORACLE / NOT DEPLOYABLE -- circ_best is built with the native",
               prereg="s31/PREREG_S31_B.md", provenance=ST.provenance(__file__))
    with io.open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)
    print("\nwrote", OUT, "  comparisons emitted", ncomp, flush=True)


if __name__ == "__main__":
    main()
