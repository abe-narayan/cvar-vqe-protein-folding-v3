"""s17/phys_ca_report.py -- the Ca-preserving repair frontier.  Analysis only, no physics.

Reads `s17/results/phys_ca_p*.json` and reports the ACTUAL Pareto frontier between all-atom
validity and Ca error, with:

  * the convergence gate's exclusion count and excluded ids at EVERY rung;
  * both mandatory controls at every rung -- do-nothing (zero-information) and a
    matched-magnitude random Ca displacement;
  * the constant ideal alpha-helix beside every validity column, because it scores
    Ramachandran 1.000 and zero clashes BY CONSTRUCTION and beat AMBER 0W/65L;
  * the rotated-lab-frame null with its MAXIMUM, against the unchanged pre-declared band.

DECLARED BEFORE READING ANY RESULT (and repeated from PREREG_phys.md): at the `cafix` rung
the Ca coordinates are frozen, so the Ca-RMSD equals the do-nothing RMSD EXACTLY.  That is a
tautology, not a measurement, and the pre-registered accuracy falsifier cannot fire there.
The rung is judged on validity; the accuracy falsifier is applied to the eps > 0 rungs.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I                     # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s17 import phys_lib as P                       # noqa: E402

GATE = 1000.0
VK = ("rama_favoured", "rama_outlier", "n_clash_2A", "n_clash_2p6A", "min_heavy",
      "bond_strain", "angle_strain", "geom_rms_rel_dev", "cis_frac", "omega_dev")


_SEQ = None


def load(tag):
    """Read a pass, and RECOMPUTE the zero-information constant-conformation references.

    They are recomputed here rather than read from the artefact because the first run of
    `phys_ca.py` built them by passing DEGREES to `core.geometry.build_backbone`, which
    takes RADIANS -- so the stored `helix` row is a different constant conformation
    (recomputed torsions -25.9 / -172.9 deg, Ramachandran-favoured 0.000).  The references
    cost nothing to rebuild and no AMBER minimisation depends on them, so the fix is applied
    at read time and the defect is recorded rather than hidden.
    """
    global _SEQ
    p = os.path.join(RESULTS, f"phys_ca_{tag}.json")
    if not os.path.exists(p):
        return None
    d = json.load(open(p))
    if _SEQ is None:
        _SEQ = {t["pdb"]: t["seq"] for t in I.targets()}
    for r in d["rows"]:
        seq = _SEQ[r["pdb"]]
        n = int(r["n"])
        nat = np.asarray(r["nat_ca"], float)
        ca_in = np.asarray(r["input_ca"], float)
        for tag2, bb in (("helix", P.helix_backbone(n)), ("strand", P.strand_backbone(n))):
            sup = I.superpose_batch(np.asarray(bb["CA"], float)[None], ca_in)[0]
            r[tag2] = {"rmsd": float(I.ca_rmsd(bb["CA"], nat)),
                       "validity": P.validity(bb, seq),
                       "ca_disp_rms": float(np.sqrt(((sup - ca_in) ** 2).sum(1).mean())),
                       "ca_disp_max": float(np.sqrt(((sup - ca_in) ** 2).sum(1)).max())}
    return d


def _col(rows, rung, f, default=np.nan):
    out = []
    for r in rows:
        v = r.get(rung)
        out.append(float(v[f]) if v is not None and f in v else default)
    return np.array(out, float)


def _val(rows, rung, key):
    out = []
    for r in rows:
        v = r.get(rung)
        out.append(float(v["validity"].get(key, np.nan)) if v is not None else np.nan)
    return np.array(out, float)


def _gatemask(rows, rung):
    if rung in ("none", "helix", "strand"):
        return np.ones(len(rows), bool)
    e = _col(rows, rung, "energy")
    return np.isfinite(e) & (e <= GATE)


def frontier(tag, rungs, title):
    d = load(tag)
    if d is None:
        print(f"\n[{tag}] not present yet.")
        return None
    rows = d["rows"]
    folds = np.array([r["fold"] for r in rows])
    names = [r["pdb"] for r in rows]
    none = _col(rows, "none", "rmsd")

    print(f"\n{'='*118}\n{title}   n = {len(rows)}   "
          f"(input = the raw all-atom coordinate average; do-nothing = "
          f"{none.mean():.4f} A)\n{'='*118}")

    print("A. THE FRONTIER.  ACCURACY is CA-RMSD to native [ORACLE].  Ca disp is the "
          "REALISED constraint quantity.")
    print(f"  {'rung':>8}{'gate ex':>8}{'Ca disp rms':>12}{'Ca disp max':>12}"
          f"{'RMSD':>9}{'vs DO-NOTHING (gated)':>32}{'vs MATCHED RANDOM':>28}{'wall s':>8}")
    tab = {}
    for rung in ("none", "helix", "strand") + tuple(rungs):
        if rung not in ("none", "helix", "strand") and rung not in rows[0]:
            continue
        ok = _gatemask(rows, rung)
        rm = _col(rows, rung, "rmsd")
        dsp = _col(rows, rung, "ca_disp_rms")
        dmx = _col(rows, rung, "ca_disp_max")
        wall = _col(rows, rung, "wall")
        s = ""
        if rung != "none":
            a = P.paired(rm[ok], none[ok], folds=folds[ok])
            s = f"{a['mean']:+.4f} [{a['ci'][0]:+.4f},{a['ci'][1]:+.4f}] {a['W']:>3}W/{a['L']:<3}L"
        s2 = ""
        rr = _col(rows, rung, "rand_rmsd_mean")
        if np.isfinite(rr).any():
            b = P.paired(rm[ok], rr[ok], folds=folds[ok])
            s2 = f"{b['mean']:+.4f} [{b['ci'][0]:+.4f},{b['ci'][1]:+.4f}] {b['W']:>3}/{b['L']:<3}"
        print(f"  {rung:>8}{int((~ok).sum()):>8}{np.nanmean(dsp):>12.4f}"
              f"{np.nanmean(dmx):>12.4f}{np.nanmean(rm[ok]):>9.4f}{s:>32}{s2:>28}"
              f"{np.nanmean(wall):>8.1f}")
        tab[rung] = {"gate_excluded": int((~ok).sum()),
                     "excluded": [names[i] for i in np.flatnonzero(~ok)],
                     "ca_disp_rms": float(np.nanmean(dsp)),
                     "ca_disp_max": float(np.nanmean(dmx)),
                     "rmsd_gated": float(np.nanmean(rm[ok])),
                     "rmsd_ungated": float(np.nanmean(rm)),
                     "vs_none": (P.paired(rm[ok], none[ok], folds=folds[ok])
                                 if rung != "none" else None),
                     "vs_rand": (P.paired(rm[ok], rr[ok], folds=folds[ok])
                                 if np.isfinite(rr).any() else None),
                     "wall": float(np.nanmean(wall))}

    print("\nB. VALIDITY at each rung, on the SAME structure whose RMSD is quoted.")
    print("   The `helix` row is a ZERO-INFORMATION reference that scores rama 1.000 and 0 "
          "clashes BY CONSTRUCTION;")
    print("   no row below is evidence about the force field except JOINTLY with its Ca disp "
          "column above.")
    hdr = (f"  {'rung':>8}{'ramaFav':>9}{'ramaOut':>9}{'cl<2.0':>8}{'cl<2.6':>8}"
           f"{'minHeavy':>10}{'bond':>8}{'angle':>8}{'geomDev':>9}{'cis':>7}{'omegaDev':>10}"
           f"{'E kcal':>11}")
    print(hdr)
    print("   EACH RUNG IS FOLLOWED BY ITS OWN GATED INPUT (`in@rung`).  The convergence gate")
    print("   selects a DIFFERENT SUBSET per rung, and comparing a gated arm's output against")
    print("   the ungated input's mean is the Sprint-16 failure mode: at `cafix` the 78")
    print("   surviving targets are the LEAST BROKEN ones, whose input already has 0.31")
    print("   sub-2.0 A clashes rather than the instrument-wide 5.61.")
    for rung in ("none", "helix", "strand") + tuple(rungs):
        if rung not in ("none", "helix", "strand") and rung not in rows[0]:
            continue
        ok = _gatemask(rows, rung)
        v = [np.nanmean(_val(rows, rung, k)[ok]) for k in VK]
        e = _col(rows, rung, "energy")
        es = f"{np.nanmean(e[ok]):>11.1f}" if np.isfinite(e).any() else f"{'--':>11}"
        print(f"  {rung:>8}{v[0]:>9.3f}{v[1]:>9.3f}{v[2]:>8.2f}{v[3]:>8.2f}{v[4]:>10.3f}"
              f"{v[5]:>8.4f}{v[6]:>8.4f}{v[7]:>9.4f}{v[8]:>7.3f}{v[9]:>10.2f}{es}")
        tab.setdefault(rung, {})["validity"] = {k: float(x) for k, x in zip(VK, v)}
        if rung not in ("none", "helix", "strand") and int((~ok).sum()):
            iv = [np.nanmean(_val(rows, "none", k)[ok]) for k in VK]
            print(f"  {'in@'+rung:>8}{iv[0]:>9.3f}{iv[1]:>9.3f}{iv[2]:>8.2f}{iv[3]:>8.2f}"
                  f"{iv[4]:>10.3f}{iv[5]:>8.4f}{iv[6]:>8.4f}{iv[7]:>9.4f}{iv[8]:>7.3f}"
                  f"{iv[9]:>10.2f}{'(gated input)':>13}")
            tab[rung]["input_on_gated_subset"] = {k: float(x) for k, x in zip(VK, iv)}

    print("\nC. THE CONJUNCTION -- validity AND staying near the input, paired against the "
          "do-nothing input.")
    print("   A rung 'wins' an axis only if it improves it with a CI excluding zero.")
    print(f"  {'rung':>8}{'d rama fav':>28}{'d clashes<2.0':>28}{'d bond strain':>28}")
    for rung in tuple(rungs):
        if rung not in rows[0]:
            continue
        ok = _gatemask(rows, rung)
        out = []
        for key, sign in (("rama_favoured", -1), ("n_clash_2A", +1), ("bond_strain", +1)):
            a = _val(rows, rung, key)[ok]; b = _val(rows, "none", key)[ok]
            m = np.isfinite(a) & np.isfinite(b)
            pr = P.paired(a[m], b[m], folds=folds[ok][m])
            out.append(f"{pr['mean']:+.4f}[{pr['ci'][0]:+.4f},{pr['ci'][1]:+.4f}]"
                       f"{pr['W'] if sign>0 else pr['L']:>4}/{pr['L'] if sign>0 else pr['W']}")
        print(f"  {rung:>8}{out[0]:>28}{out[1]:>28}{out[2]:>28}")
    return tab


def frame_null(tag="p4"):
    d = load(tag)
    if d is None:
        print("\n[p4] frame null not present yet.")
        return None
    rows = [r for r in d["rows"] if "frame" in r]
    if not rows:
        return None
    out = {}
    print(f"\n{'='*118}\nD. THE ROTATED-LAB-FRAME NULL -- ZERO BY CONSTRUCTION.  Proper "
          f"rotations only (det = +1).\n"
          f"   Pre-declared PASS band, unchanged from s16/energy_lib: |mean| <= "
          f"{EL.FRAME_TOL_MEAN} A AND max <= {EL.FRAME_TOL_MAX} A on the CONVERGED subset.\n"
          f"{'='*118}")
    print(f"  {'rung':>8}{'n':>5}{'gated n':>9}{'mean':>11}{'sd':>9}{'MAX |d|':>11}"
          f"{'n>1e-6':>9}  verdict")
    for rung in sorted({k for r in rows for k in r["frame"]}):
        a = np.array([r[rung]["rmsd"] for r in rows if rung in r and rung in r["frame"]])
        b = np.array([r["frame"][rung]["rmsd"] for r in rows if rung in r and rung in r["frame"]])
        ea = np.array([r[rung]["energy"] for r in rows if rung in r and rung in r["frame"]])
        eb = np.array([r["frame"][rung]["energy"] for r in rows if rung in r and rung in r["frame"]])
        dif = b - a
        ok = np.isfinite(ea) & (ea <= GATE) & np.isfinite(eb) & (eb <= GATE)
        g = dif[ok]
        if len(g) == 0:
            print(f"  {rung:>8}{len(dif):>5}{0:>9}   every target excluded by the gate")
            out[rung] = {"n": int(len(dif)), "gated_n": 0}
            continue
        verdict = ("PASS" if abs(g.mean()) <= EL.FRAME_TOL_MEAN
                   and np.abs(g).max() <= EL.FRAME_TOL_MAX else "FAIL")
        print(f"  {rung:>8}{len(dif):>5}{len(g):>9}{g.mean():>11.5f}{g.std():>9.5f}"
              f"{np.abs(g).max():>11.5f}{int((np.abs(g) > 1e-6).sum()):>9}  {verdict}")
        out[rung] = {"n": int(len(dif)), "gated_n": int(len(g)), "mean": float(g.mean()),
                     "sd": float(g.std()), "max": float(np.abs(g).max()),
                     "n_nonzero": int((np.abs(g) > 1e-6).sum()), "verdict": verdict}
    return out


def reproduce_s16():
    """The instrument's own reproduction check: the k30 rung of pass 1 against
    `s16/results/repair_A.json`, which was itself bit-identical to Sprint 15."""
    d = load("p1")
    if d is None:
        return None
    try:
        s16 = {r["pdb"]: r for r in
               json.load(open(os.path.join(ROOT, "s16", "results", "repair_A.json")))
               ["per_target"]}
    except Exception as ex:
        print(f"  (s16 repair_A not readable: {ex})")
        return None
    dr, de, n = [], [], 0
    for r in d["rows"]:
        s = s16.get(r["pdb"])
        if s is None or "k30" not in r:
            continue
        dr.append(abs(r["k30"]["rmsd"] - s["f0_k30_full"]["rmsd"]))
        de.append(abs(r["k30"]["energy"] - s["f0_k30_full"]["energy"]))
        n += 1
    if not n:
        return None
    print(f"\nE. INSTRUMENT REPRODUCTION.  The k30 rung is built here from an independently "
          f"constructed\n   OpenMM System, not from `core.amber.refine_coords`.  Against "
          f"`s16/results/repair_A.json`\n   (itself bit-identical to Sprint 15) on n = {n} "
          f"targets:")
    print(f"     max |dRMSD|   = {max(dr):.3e} A")
    print(f"     max |dEnergy| = {max(de):.3e} kcal/mol")
    return {"n": n, "max_drmsd": float(max(dr)), "max_denergy": float(max(de))}


def main():
    out = {}
    out["p1"] = frontier("p1", ("cafix", "k30"),
                         "E1 / PASS 1 -- Ca-FIXED AMBER at FULL n, against the incumbent")
    out["p2"] = frontier("p2", ("cafix", "ca1000", "ca300", "ca100", "ca30", "ca10", "ca3",
                                "free", "k30"),
                         "E1 / PASS 2 -- THE FULL LADDER on the pre-registered 30-target "
                         "subsample (SHAPE, not a selection)")
    out["p3"] = frontier("p3", ("eps010", "eps025"),
                         "E1 / PASS 3 -- flat-bottom eps cross-check (both "
                         "parameterisations, one frontier)")
    out["frame"] = frame_null("p4")
    out["reproduce"] = reproduce_s16()
    json.dump(out, open(os.path.join(RESULTS, "phys_ca_report.json"), "w"),
              default=lambda x: float(x) if isinstance(x, np.floating) else str(x))
    return out


if __name__ == "__main__":
    main()
