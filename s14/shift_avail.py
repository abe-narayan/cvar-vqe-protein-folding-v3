"""SHIFT agent -- STEP 1.  How much backbone chemical-shift information actually exists
for the 126 development targets?

Consumes s14/results/shift_availability.json (built by s14.shift_bmrb) and reports the
availability tiers, per-nucleus presence, per-residue and per-position coverage, terminal
coverage, and the breakdown by FAIL18 / structural class.

Reads no native coordinates.  Chemical-shift assignments and deposition metadata only.

    python -m s14.shift_avail
"""
from __future__ import annotations

import os
import sys
import json

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s14.shift_bmrb import BACKBONE, RESULTS          # noqa: E402

AVAIL = os.path.join(RESULTS, "shift_availability.json")

# BMRB accessions in the 36000-36999 block are PDBj-BMRB depositions; the BMRB REST API
# does not serve them (verified: absent from /v2/list_entries?database=macromolecules).
PDBJ_LO, PDBJ_HI = 36000, 36999

#: Normalised cost weight of a single-residue torsion error by fractional chain position,
#: ten bins N-terminus -> C-terminus.  MEASURED by the coordinator
#: (`s14/position.py` -> `s14/results/position.json`, ORACLE DIAGNOSTIC, 126 targets,
#: single-residue perturbation at sigma 20 deg, 24 reps/residue).  The middle 40 % of a
#: chain carries 62.8 % of the total single-residue cost and the outer 40 % carries 15.7 %.
#: Mechanism: a torsion at position p hinges two rigid segments of length p and n-p, and
#: under Kabsch superposition the displacement scales with the lever arm product p(n-p).
POS_W = np.array([0.012, 0.060, 0.103, 0.144, 0.170,
                  0.169, 0.145, 0.112, 0.070, 0.015])


def position_weights(n):
    """Per-residue cost weight for a chain of length n, normalised to sum to 1."""
    idx = np.minimum((np.arange(n) + 0.5) / n * 10.0, 9.999).astype(int)
    w = POS_W[idx]
    return w / w.sum()


def weighted_coverage(mask):
    """Position-weighted coverage: the fraction of the chain's total single-residue
    torsion-error COST that the restrained residues account for.  This is the axis the
    (sigma, coverage) surface should really be read on -- 75 % coverage concentrated
    mid-chain is worth more than 90 % concentrated at the termini."""
    m = np.asarray(mask, bool)
    w = position_weights(len(m))
    return float(w[m].sum())


def tiers(rec):
    """Which TALOS-class input tier does this deposition support?

    T0  no deposition retrievable
    T1  shifts exist but are 1-H only (no CA/CB/C'/N anywhere)  -> TALOS-N cannot run
    T2  heteronuclear shifts exist but the completeness gate passes on < 75 % of residues
    T3  the TALOS-N completeness gate passes on >= 75 % of residues
    T4  ... on >= 90 % of residues (the regime Sprint 12/13 priced at <= 12 deg)
    """
    if rec is None:
        return 0
    hetero = any(sum(rec["per_nuc"][a]) > 0 for a in ("N", "CA", "CB", "C"))
    if not hetero:
        return 1
    c = rec["cov_talos"]
    if c >= 0.90:
        return 4
    if c >= 0.75:
        return 3
    return 2


def main():
    from s12 import instrument as I

    with open(AVAIL) as fh:
        A = json.load(fh)
    tg = I.targets()
    fail18 = set(I.FAIL18)

    try:
        with open(os.path.join(ROOT, "s12", "results", "forensics_class.json")) as fh:
            cls = {r["pdb"]: r.get("tclass", []) for r in json.load(fh)["rows"]}
    except Exception:                                              # noqa: BLE001
        cls = {}

    rows = []
    for t in tg:
        pdb = t["pdb"]
        a = A[pdb]
        b = a["best"]
        pdbj = [c for c in a["candidates"] if c.isdigit() and PDBJ_LO <= int(c) <= PDBJ_HI]
        r = {
            "pdb": pdb, "n": t["n"], "fold": t["fold"], "fail18": pdb in fail18,
            "tclass": cls.get(pdb, []),
            "n_cand": a["n_candidates"], "n_pdbj_unresolved": len(pdbj),
            "tier": tiers(b),
        }
        if b is None:
            r.update({"bmrb": None, "cov_any": 0.0, "cov_ge3": 0.0, "cov_talos": 0.0,
                      "cov_full6": 0.0,
                      "nuc": {a2: 0.0 for a2 in BACKBONE},
                      "term_cov": 0.0, "int_cov": 0.0, "src": []})
        else:
            n = t["n"]
            pn = {a2: np.array(b["per_nuc"][a2], float) for a2 in BACKBONE}
            nshift = np.sum([pn[a2] for a2 in BACKBONE], axis=0).astype(int)
            m = np.array(b["talos_mask"], float)
            k = min(3, n // 3)
            term = np.zeros(n, bool)
            term[:k] = True
            term[-k:] = True
            r.update({
                "bmrb": b["bmrb"], "src": b["src"], "entity_len": b["entity_len"],
                "entity_name": b["entity_name"], "off": b["off"],
                "cov_any": b["cov_any"], "cov_ge3": b["cov_ge3"],
                "cov_talos": b["cov_talos"], "cov_full6": b["cov_full6"],
                "wcov_talos": weighted_coverage(b["talos_mask"]),
                "wcov_ge3": weighted_coverage(nshift >= 3),
                "nuc": {a2: float(pn[a2].mean()) for a2 in BACKBONE},
                "term_cov": float(m[term].mean()),
                "int_cov": float(m[~term].mean()) if (~term).any() else float("nan"),
                "talos_mask": b["talos_mask"],
                "n_shifts_per_res": [int(x) for x in nshift],
            })
        rows.append(r)

    out = {"rows": rows}
    T = np.array([r["tier"] for r in rows])
    isf = np.array([r["fail18"] for r in rows])
    fibril = np.array([("fibril" in r["tclass"]) or ("lasso" in r["tclass"]) for r in rows])

    def frac(mask):
        return float(mask.mean()), int(mask.sum())

    print("=" * 78)
    print("STEP 1 -- MEASURED BMRB BACKBONE-SHIFT AVAILABILITY, 126 dev targets")
    print("=" * 78)
    print()
    names = {0: "T0 no deposition retrieved",
             1: "T1 1-H only (no CA/CB/C'/N)  -- TALOS-N cannot run",
             2: "T2 hetero, completeness gate < 75 %",
             3: "T3 hetero, gate 75-90 %",
             4: "T4 hetero, gate >= 90 %  (the priced regime)"}
    for k in range(5):
        m = T == k
        print("  {:52s} {:3d}  ({:.3f})".format(names[k], int(m.sum()), m.mean()))
    print()
    print("  any BMRB candidate at all            : {:3d} ({:.3f})".format(
        *reversed(frac(np.array([r["n_cand"] > 0 for r in rows])))))
    print("  usable deposition retrieved (T>=1)   : {:3d} ({:.3f})".format(
        *reversed(frac(T >= 1))))
    print("  heteronuclear backbone shifts (T>=2) : {:3d} ({:.3f})".format(
        *reversed(frac(T >= 2))))
    print("  TALOS-N-runnable at >=75 % (T>=3)    : {:3d} ({:.3f})".format(
        *reversed(frac(T >= 3))))
    print("  TALOS-N-runnable at >=90 % (T4)      : {:3d} ({:.3f})".format(
        *reversed(frac(T == 4))))
    npdbj = sum(1 for r in rows if r["tier"] == 0 and r["n_pdbj_unresolved"] > 0)
    print("  of the T0, PDBj-BMRB (36xxx) not served by the API: {:d}".format(npdbj))
    print()

    print("  per-nucleus presence, mean over residues (T>=2 targets only):")
    sel = [r for r in rows if r["tier"] >= 2]
    for a2 in BACKBONE:
        v = np.array([r["nuc"][a2] for r in sel])
        print("    {:3s}  mean {:.3f}   targets with any {:3d}/{:d}".format(
            a2, v.mean(), int((v > 0).sum()), len(sel)))
    print()

    print("  coverage distribution among T>=2 (the TALOS-N completeness gate):")
    c = np.array([r["cov_talos"] for r in sel])
    for q in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        print("    q{:.2f}  {:.3f}".format(q, float(np.quantile(c, q))))
    print("    mean {:.3f}   frac==1.00 {:.3f}   frac>=0.90 {:.3f}   frac>=0.75 {:.3f}"
          .format(c.mean(), float((c >= 0.999).mean()), float((c >= 0.90).mean()),
                  float((c >= 0.75).mean())))
    print()

    wc = np.array([r["wcov_talos"] for r in sel])
    print("  POSITION-WEIGHTED gate coverage (coordinator's s14/results/position.json law):")
    print("    raw mean {:.3f}   weighted mean {:.3f}   diff {:+.3f}".format(
        c.mean(), wc.mean(), wc.mean() - c.mean()))
    print("    weighted: frac>=0.90 {:.3f}   frac>=0.75 {:.3f}   min {:.3f}".format(
        float((wc >= 0.90).mean()), float((wc >= 0.75).mean()), float(wc.min())))
    lo = np.argsort(c)[:6]
    print("    the six lowest RAW-coverage targets, raw -> weighted:")
    for i in lo:
        print("      {:6s} raw {:.3f} -> weighted {:.3f}".format(
            sel[i]["pdb"], sel[i]["cov_talos"], sel[i]["wcov_talos"]))
    print()

    print("  TERMINAL vs INTERIOR gate coverage (T>=2), termini = first/last min(3,n//3):")
    tc = np.array([r["term_cov"] for r in sel])
    ic = np.array([r["int_cov"] for r in sel])
    ok = ~np.isnan(ic)
    print("    terminal {:.3f}   interior {:.3f}   diff {:+.3f}".format(
        tc[ok].mean(), ic[ok].mean(), tc[ok].mean() - ic[ok].mean()))
    print()

    print("  by subgroup (fraction reaching each tier):")
    for lbl, m in (("all 126", np.ones(126, bool)), ("FAIL18", isf),
                   ("other108", ~isf), ("fibril/lasso", fibril),
                   ("not fibril/lasso", ~fibril)):
        if m.sum() == 0:
            continue
        print("    {:18s} n={:3d}  T>=2 {:.3f}  T>=3 {:.3f}  T4 {:.3f}  mean cov {:.3f}"
              .format(lbl, int(m.sum()), float((T[m] >= 2).mean()),
                      float((T[m] >= 3).mean()), float((T[m] == 4).mean()),
                      float(np.mean([rows[i]["cov_talos"] for i in np.flatnonzero(m)]))))
    print()

    out["summary"] = {
        "n_targets": 126,
        "tier_counts": {str(k): int((T == k).sum()) for k in range(5)},
        "n_pdbj_unresolved_T0": npdbj,
        "cov_talos_mean_T2plus": float(c.mean()),
        "wcov_talos_mean_T2plus": float(wc.mean()),
        "frac_cov_ge_0.90_of_T2plus": float((c >= 0.90).mean()),
        "frac_wcov_ge_0.90_of_T2plus": float((wc >= 0.90).mean()),
        "position_weights_10bin": POS_W.tolist(),
        "term_cov_T2plus": float(tc[ok].mean()),
        "int_cov_T2plus": float(ic[ok].mean()),
        "nuclei_presence_T2plus": {a2: float(np.mean([r["nuc"][a2] for r in sel]))
                                   for a2 in BACKBONE},
        "by_group": {lbl: {"n": int(m.sum()),
                           "T2plus": float((T[m] >= 2).mean()),
                           "T3plus": float((T[m] >= 3).mean()),
                           "T4": float((T[m] == 4).mean())}
                     for lbl, m in (("all", np.ones(126, bool)), ("FAIL18", isf),
                                    ("other108", ~isf), ("fibril_lasso", fibril))},
    }
    out["runnable_T3plus"] = sorted(r["pdb"] for r in rows if r["tier"] >= 3)
    out["runnable_T4"] = sorted(r["pdb"] for r in rows if r["tier"] == 4)
    out["T2plus"] = sorted(r["pdb"] for r in rows if r["tier"] >= 2)
    out["what"] = ("MEASURED availability of BMRB backbone chemical shifts for the 126 dev "
                   "targets. No native structural information read.")
    with open(os.path.join(RESULTS, "shift_avail_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote s14/results/shift_avail_summary.json")
    print("  T>=3 runnable set (n={}): {}".format(
        len(out["runnable_T3plus"]), " ".join(out["runnable_T3plus"])))


if __name__ == "__main__":
    main()
