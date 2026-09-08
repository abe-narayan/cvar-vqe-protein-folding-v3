"""SHIFT agent -- STEP 5.  The honest verdict, and the error breakdowns STEP 3 asked for.

Consumes the saved 324-cell posteriors (s14/cache/shift_post_*.npz), the availability table
and the pricing surface, and answers:

  * for how many targets could a shift-restrained system actually run,
  * what it emits on those targets,
  * how that compares to the incumbent ON THE SAME SUBSET (never the full-set number),
  * where the error sits: by secondary structure, residue type, peptide class, FAIL18,
  * what the full instrument would look like as a mixed arm, clearly labelled as mixed,
  * whether the kill threshold is met.

    python -m s14.shift_verdict
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                     # noqa: E402
from s13 import tors_common as T                    # noqa: E402
from s13 import tors_eval as EV                     # noqa: E402
from s14.shift_bmrb import RESULTS                  # noqa: E402

CACHE = os.path.join(ROOT, "s14", "cache")
R2D = 180.0 / np.pi
D2R = np.pi / 180.0


def ss_of(phi, psi):
    """Coarse three-state assignment from the native torsions.  ORACLE, used only to
    STRATIFY the error, never to make a prediction."""
    p = np.degrees(T.wrap(phi)); s = np.degrees(T.wrap(psi))
    out = np.full(len(p), "C", dtype="<U1")
    out[(p > -160) & (p < -20) & (s > -120) & (s < 50)] = "H"
    out[(p > -180) & (p < -40) & ((s > 90) | (s < -150))] = "E"
    out[p > 0] = "L"          # left-handed / positive phi
    return out


def main():
    t0 = time.time()
    with open(os.path.join(RESULTS, "shift_avail_summary.json")) as fh:
        S = json.load(fh)
    with open(os.path.join(RESULTS, "shift_model.json")) as fh:
        M = json.load(fh)
    with open(os.path.join(RESULTS, "shift_price.json")) as fh:
        PR = json.load(fh)
    avail = {r["pdb"]: r for r in S["rows"]}

    tg = I.targets()
    tgd = {t["pdb"]: t for t in tg}
    f18 = set(I.FAIL18)
    nat = EV.native_torsions()
    inc = EV.incumbent()
    run = sorted(S["runnable_T3plus"])
    notrun = [t["pdb"] for t in tg if t["pdb"] not in set(run)]

    out = {"n_runnable": len(run), "runnable": run,
           "n_not_runnable": len(notrun),
           "incumbent_runnable": float(np.mean([inc[p] for p in run])),
           "incumbent_not_runnable": float(np.mean([inc[p] for p in notrun])),
           "incumbent_all126": float(np.mean([inc[t["pdb"]] for t in tg]))}

    print("=" * 84)
    print("STEP 5 -- VERDICT")
    print("=" * 84)
    print("runnable {}/126   incumbent-on-runnable {:.4f}   incumbent-elsewhere {:.4f}"
          .format(len(run), out["incumbent_runnable"], out["incumbent_not_runnable"]))
    print()

    # ------------------------------------------------- error breakdown on E-SHIFT
    z = np.load(os.path.join(CACHE, "shift_post_E-SHIFT.npz"))
    zs = np.load(os.path.join(CACHE, "shift_post_SEQ-ONLY.npz"))
    rec = []
    for pid in run:
        n = tgd[pid]["n"]
        nphi, npsi = nat[pid]
        seq = tgd[pid]["seq"]
        ss = ss_of(nphi, npsi)
        gate = np.array(avail[pid]["talos_mask"], bool)
        nsh = np.array(avail[pid]["n_shifts_per_res"])
        for arm, P in (("E-SHIFT", np.asarray(z[pid], np.float64)),
                       ("SEQ-ONLY", np.asarray(zs[pid], np.float64))):
            d = T.decode(P)
            ep = np.degrees(T.wrap(d["phi"] - nphi))
            es = np.degrees(T.wrap(d["psi"] - npsi))
            for i in range(n):
                rec.append({"pdb": pid, "arm": arm, "i": i, "n": n, "aa": seq[i],
                            "ss": ss[i], "gate": bool(gate[i]), "nsh": int(nsh[i]),
                            "fail18": pid in f18,
                            "ephi": float(ep[i]) if i > 0 else np.nan,
                            "epsi": float(es[i]) if i < n - 1 else np.nan,
                            "sigma_hat": float(d["sigma_hat"][i]),
                            "gross": bool(np.hypot(ep[i], es[i]) > 60.0)})

    def rms(v):
        v = np.asarray([x for x in v if np.isfinite(x)], float)
        return float(np.sqrt((v ** 2).mean())) if len(v) else float("nan")

    def block(sel, label, arm="E-SHIFT"):
        r = [x for x in rec if x["arm"] == arm and sel(x)]
        if not r:
            return None
        ph = rms([x["ephi"] for x in r]); ps = rms([x["epsi"] for x in r])
        both = [x for x in r if np.isfinite(x["ephi"]) and np.isfinite(x["epsi"])]
        good = [x for x in r if not x["gross"]]
        gph = rms([x["ephi"] for x in good]); gps = rms([x["epsi"] for x in good])
        return {"label": label, "n": len(r), "rms_phi": ph, "rms_psi": ps,
                "sigma": float(np.sqrt((ph ** 2 + ps ** 2) / 2)),
                "sigma_good_only": float(np.sqrt((gph ** 2 + gps ** 2) / 2)),
                "gross_rate": float(np.mean([x["gross"] for x in both])) if both else float("nan")}

    print("E-SHIFT angular error, stratified (ORACLE stratification, DEMONSTRATED errors):")
    print("{:26s} {:>6s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
        "stratum", "n", "RMSphi", "RMSpsi", "sigma", "sig|good", "gross"))
    strata = []
    for lbl, fn in [
        ("ALL", lambda x: True),
        ("gate PASS", lambda x: x["gate"]),
        ("gate FAIL", lambda x: not x["gate"]),
        ("SS: helix", lambda x: x["ss"] == "H"),
        ("SS: extended", lambda x: x["ss"] == "E"),
        ("SS: coil", lambda x: x["ss"] == "C"),
        ("SS: positive-phi", lambda x: x["ss"] == "L"),
        ("GLY", lambda x: x["aa"] == "G"),
        ("PRO", lambda x: x["aa"] == "P"),
        ("not GLY/PRO", lambda x: x["aa"] not in "GP"),
        ("FAIL18", lambda x: x["fail18"]),
        ("other (runnable)", lambda x: not x["fail18"]),
        ("terminal 2 residues", lambda x: x["i"] < 2 or x["i"] >= x["n"] - 2),
        ("interior", lambda x: 2 <= x["i"] < x["n"] - 2),
        ("6 nuclei present", lambda x: x["nsh"] == 6),
        ("5 nuclei", lambda x: x["nsh"] == 5),
        ("4 nuclei", lambda x: x["nsh"] == 4),
        ("<=3 nuclei", lambda x: x["nsh"] <= 3),
    ]:
        b = block(fn, lbl)
        if b:
            strata.append(b)
            print("{:26s} {:6d} {:8.1f} {:8.1f} {:8.1f} {:8.1f} {:8.3f}".format(
                lbl, b["n"], b["rms_phi"], b["rms_psi"], b["sigma"],
                b["sigma_good_only"], b["gross_rate"]))
    out["strata_E_SHIFT"] = strata
    sq = block(lambda x: True, "ALL", arm="SEQ-ONLY")
    out["strata_SEQ_ONLY_all"] = sq
    print("{:26s} {:6d} {:8.1f} {:8.1f} {:8.1f} {:8.1f} {:8.3f}  <- SEQ-ONLY".format(
        "SEQ-ONLY ALL", sq["n"], sq["rms_phi"], sq["rms_psi"], sq["sigma"],
        sq["sigma_good_only"], sq["gross_rate"]))
    print()

    # ------------------------------------------------- the mixed full-instrument arm
    print("MIXED arm (shift-restrained where it runs, incumbent elsewhere) -- LABELLED MIXED:")
    for arm in ("E-SHIFT", "SHUF-SHIFT", "P-SHIFT", "SEQ-ONLY"):
        pt = M[arm]["per_target"]
        mix = np.array([pt[p] if p in pt else inc[p] for p in [t["pdb"] for t in tg]])
        base = np.array([inc[t["pdb"]] for t in tg])
        pr = I.paired(mix, base, folds=[t["fold"] for t in tg],
                      names=[t["pdb"] for t in tg])
        out["mixed_" + arm] = {"mean": float(mix.mean()), "paired_vs_incumbent_all126": pr}
        print("  {:10s} full-126 mixed mean {:.4f}  vs incumbent {:.4f}  diff {:+.4f} "
              "CI [{:+.4f}, {:+.4f}]  W/L {}/{}".format(
                  arm, mix.mean(), base.mean(), pr["mean_diff"],
                  pr["ci95"][0], pr["ci95"][1], pr["n_better"], pr["n_worse"]))
    print()

    # ------------------------------------------------- ORACLE ceiling arms on the subset
    print("Ladder on the {} runnable targets (incumbent {:.4f}):".format(
        len(run), out["incumbent_runnable"]))
    lad = [
        ("ORACLE sigma=0 (build floor)", PR["A_measured_mask"]["measured_sig0"]["mean"], "ORACLE"),
        ("ORACLE sigma=12, measured mask", PR["A_measured_mask"]["measured_sig12"]["mean"], "ORACLE"),
        ("ORACLE TALOS-N mixture, oracle basins", PR["C_talosn"]["oracle2"]["mean"], "ORACLE"),
        ("ORACLE TALOS-N mixture, drop ambiguous", PR["C_talosn"]["drop"]["mean"], "ORACLE"),
        ("ORACLE TALOS-N mixture, collapse", PR["C_talosn"]["collapse"]["mean"], "ORACLE"),
        ("E-SHIFT (measured shifts, real model)", M["E-SHIFT"]["emitted_gated"], "NMR-RESTRAINED"),
        ("SHUF-SHIFT (shifts permuted, null)", M["SHUF-SHIFT"]["emitted_gated"], "NULL"),
        ("P-SHIFT (predicted shifts)", M["P-SHIFT"]["emitted_gated"], "sequence-only"),
        ("REF-SHIFT (zero secondary shift null)", M["REF-SHIFT"]["emitted_gated"], "NULL"),
        ("MASK-ONLY (assignment pattern only)", M["MASK-ONLY"]["emitted_gated"], "NULL"),
        ("SEQ-ONLY (s13 a_pepPos)", M["SEQ-ONLY"]["emitted_gated"], "sequence-only"),
        ("incumbent on the same 54", out["incumbent_runnable"], "shipped"),
    ]
    for lbl, v, tier in sorted(lad, key=lambda x: x[1]):
        print("  {:42s} {:7.3f}   [{}]".format(lbl, v, tier))
    out["ladder"] = [{"label": a, "value": b, "tier": c} for a, b, c in lad]

    # ------------------------------------------- the asymmetry test: E-SHIFT vs each null
    print()
    print("THE ASYMMETRY TEST -- E-SHIFT paired against every null, {} targets:"
          .format(len(run)))
    eg = np.array([M["E-SHIFT"]["per_target"][p] for p in run])
    folds = [tgd[p]["fold"] for p in run]
    out["asymmetry"] = {}
    for other in ("SHUF-SHIFT", "P-SHIFT", "REF-SHIFT", "MASK-ONLY", "SEQ-ONLY"):
        ov = np.array([M[other]["per_target"][p] for p in run])
        pr = I.paired(eg, ov, folds=folds, names=run)
        out["asymmetry"][other] = pr
        print("  E-SHIFT vs {:11s} {:+.4f}  CI [{:+.4f}, {:+.4f}]  W/L {:2d}/{:2d}  "
              "drop10 {:+.4f}  drop20 {:+.4f}".format(
                  other, pr["mean_diff"], pr["ci95"][0], pr["ci95"][1],
                  pr["n_better"], pr["n_worse"],
                  pr["drop_top10_mean_diff"], pr["drop_top20_mean_diff"]))
        print("      per-fold: {}".format(
            {k: round(v, 3) for k, v in pr["per_fold"].items()}))

    # ------------------------------- the availability arithmetic that settles the question
    print()
    print("AVAILABILITY ARITHMETIC -- what the FULL 126-target instrument can reach if the")
    print("shift channel is perfect on every target where it runs and the incumbent is used")
    print("everywhere else.  a = {}/126 = {:.4f}".format(len(run), len(run) / 126.0))
    a = len(run) / 126.0
    rest = out["incumbent_not_runnable"]
    ceil = {}
    for lbl, v in [("ORACLE sigma=0, the ideal-geometry BUILD FLOOR",
                    PR["A_measured_mask"]["measured_sig0"]["mean"]),
                   ("ORACLE sigma=12 on the measured mask",
                    PR["A_measured_mask"]["measured_sig12"]["mean"]),
                   ("ORACLE TALOS-N mixture, oracle basin choice",
                    PR["C_talosn"]["oracle2"]["mean"]),
                   ("ORACLE TALOS-N mixture, TALOS-N behaviour",
                    PR["C_talosn"]["drop"]["mean"]),
                   ("E-SHIFT, the real measured-shift model",
                    M["E-SHIFT"]["emitted_gated"])]:
        full = a * v + (1 - a) * rest
        ceil[lbl] = {"subset": v, "full126": full}
        print("  {:46s} subset {:6.3f} -> full-126 {:6.3f}".format(lbl, v, full))
    # upper bound if the PDBj-BMRB block were accessible (61/126)
    a2 = 61 / 126.0
    v0 = PR["A_measured_mask"]["measured_sig0"]["mean"]
    print("  with the 7 PDBj-BMRB targets recovered (a = 61/126 = {:.4f}):".format(a2))
    print("     perfect torsions on all 61 -> full-126 {:6.3f}".format(
        a2 * v0 + (1 - a2) * rest))
    print("  availability needed for full-126 < 2.0 A with PERFECT torsions: a > {:.3f} "
          "({:.0f} targets)".format((rest - 2.0) / (rest - v0),
                                    126 * (rest - 2.0) / (rest - v0)))
    out["availability_arithmetic"] = {
        "a": a, "incumbent_elsewhere": rest, "ceilings": ceil,
        "a_needed_for_2.0_with_perfect_torsions": float((rest - 2.0) / (rest - v0)),
        "n_targets_needed": float(126 * (rest - 2.0) / (rest - v0))}

    out["what"] = ("STEP 5 verdict. Comparisons are on the 54-target matched subset. "
                   "The MIXED arm is labelled mixed and is not a like-for-like improvement "
                   "to the sequence-only pipeline: E-SHIFT is NMR-restrained.")
    out["secs"] = round(time.time() - t0, 1)
    with open(os.path.join(RESULTS, "shift_verdict.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print()
    print("wrote s14/results/shift_verdict.json  [{:.0f}s]".format(time.time() - t0))


if __name__ == "__main__":
    main()
