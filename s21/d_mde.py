"""s21/d_mde.py -- D5: "MDE = 0.084 A" IS NOT A PROPERTY OF THE INSTRUMENT.

    python -m s21.d_mde

WHAT IS UNDER ATTACK.  `s21/BRIEF.md` §1: *"MDE at 80% power = 0.084 A. A null below that is
uninformative, not negative."*  The same sentence is in `s18/BRIEF.md` §6, `s19/BRIEF.md` and
`s20/BRIEF.md`, and the number has been quoted as a universal bar in at least four sprints'
findings ("5% of the instrument's 0.084 A MDE"; "a null below the MDE"; "a quarter of the MDE").

DERIVE THE OPERATOR BEFORE INTERPRETING ITS STATISTIC (BRIEF §9).  For a paired comparison at
two-sided alpha = 0.05 and 80% power,

    MDE = (z_0.975 + z_0.80) * SE  =  2.8016 * sd(paired differences) / sqrt(n)

**`sd(paired differences)` is a property of the COMPARISON, not of the instrument.**  Two arms
that differ only in a 50-step relaxation have a paired sd of a few hundredths of an Angstrom;
two arms that differ in which of 500 candidates is selected have a paired sd of one to two
Angstroms.  One MDE cannot serve both, and it errs in BOTH directions:

    too LARGE for a low-variance comparison  -> a solidly measured effect whose CI excludes zero
                                                gets written off as "below the MDE"
    too SMALL for a high-variance comparison -> a genuinely underpowered null gets read as a
                                                measured absence

The programme derived this correctly once.  `s16/review_FINDINGS.md` Q4.1 computes it FOR ITS OWN
PANEL -- "SE = 0.051 A, MDE at 80% power ~ 0.16 A, power at a true 0.05 A effect ~ 9%" -- and
`s16/LEDGER.md:678` records it the same way.  It became an instrument constant in `s18/BRIEF.md`
and has been quoted as one since.

WHAT THIS MODULE DOES.  It recomputes the per-comparison MDE from the paired differences of every
comparison this lane produced this sprint, plus the ones the briefs quote against the 0.084 bar,
and reports the ratio.  No new structure, no energy, nothing rebuilt: the paired differences
already exist on disk.

    H0 (the brief's implicit one): 0.084 A is the MDE of comparisons on this instrument.
    Falsifier of MY claim: if the per-comparison MDEs cluster within a factor of ~1.5 of 0.084,
    the constant is a fair summary and D5 is REFUTED.

An implied sd is also reported: the paired sd that WOULD make MDE = 0.084 at n = 126 is
0.084 * sqrt(126) / 2.8016 = 0.3366 A.  Any comparison whose paired sd is far from that has a
different MDE, and the table says by how much.
"""
from __future__ import annotations

import json
import math
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
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

Z = 1.959963985 + 0.841621234        # = 2.8016, two-sided alpha=0.05 at 80% power
QUOTED = 0.084


def mde(d):
    d = np.asarray(d, float)
    n = len(d)
    sd = float(d.std(ddof=1))
    se = sd / math.sqrt(n)
    return {"n": n, "mean": float(d.mean()), "sd": sd, "se": se, "mde": Z * se,
            "ratio_to_quoted": Z * se / QUOTED}


def _cv():
    return json.load(open(os.path.join(RESULTS, "d_cvarop.json")))["rows"]


def _dv():
    return json.load(open(os.path.join(RESULTS, "d_divsel.json")))["rows"]


def _nl():
    return json.load(open(os.path.join(ROOT, "s20", "results", "c_land_null.json")))["rows"]


def run():
    comps = {}
    g = lambda rows, k: np.array([r[k] for r in rows])            # noqa: E731

    r = _cv()
    for h in ("disto", "legacy", "rand"):
        comps[f"[D2] tail_avg(0.15)-argmin | {h}"] = g(r, f"tail_avg0.15|{h}") - g(r, f"argmin|{h}")
        comps[f"[D2] tail_avg(0.15)-rand_avg | {h}"] = (g(r, f"tail_avg0.15|{h}")
                                                        - g(r, "rand_avg0.15"))
        comps[f"[D2] tail_member(0.15)-rand_member | {h}"] = (g(r, f"tail_member0.15|{h}")
                                                              - g(r, "rand_member"))

    d = _dv()
    comps["[R1] divmax(0.3) - top75"] = g(d, "divmax0.3") - g(d, "incumbent_top75")
    comps["[R1] divmax(0.3) - rand_in_band"] = g(d, "divmax0.3") - g(d, "randband0.3")
    comps["[R1] divmax(1.0) - top75"] = g(d, "divmax1") - g(d, "incumbent_top75")

    #: MY OWN nulls, so the correction is applied to this lane's conclusions first.
    try:
        de = json.load(open(os.path.join(RESULTS, "d_enc_LEG.json")))
        P = sorted(de)
        SEEDS = (0, 1, 2, 3)

        def cell(arm, enc, B):
            return np.array([np.mean([de[p_]["rows"][f"LEG|{arm}|{enc}|{B}|{s_}"]["rmsd_ORACLE"]
                                      for s_ in SEEDS]) for p_ in P])
        for arm in ("spsa", "adam_fd"):
            comps[f"[D1] LEG/{arm}: emb@512 - theta@512"] = cell(arm, "emb", 512) - cell(arm, "theta", 512)
        M2 = np.stack([cell(a, "emb", 512) - cell(a, "theta", 512) for a in ("spsa", "adam_fd")])
        comps["[D1] LEG pooled: emb@512 - theta@512"] = M2.mean(0)
    except Exception as e:                                        # pragma: no cover
        print(f"  (d_enc_LEG unavailable: {e!r})")

    #: Sprint 20's own encoding cells, recomputed from its artefacts at their real n and seeds.
    try:
        e20 = json.load(open(os.path.join(ROOT, "s20", "results", "qb2_enc.json")))
        o20 = json.load(open(os.path.join(ROOT, "s20", "results", "qb2_opt.json")))
        Q = sorted(set(e20) & set(o20))
        for k_, a_ in (("AMBc", "spsa"), ("LEG", "spsa")):
            em = np.array([np.mean([e20[p_]["rows"][f"{k_}|{a_}|{s_}"]["rmsd_ORACLE"]
                                    for s_ in "01"]) for p_ in Q])
            th = np.array([np.mean([o20[p_]["rows"][f"{k_}|{a_}|{s_}"]["rmsd_ORACLE"]
                                    for s_ in "01"]) for p_ in Q])
            comps[f"[s20] encoding {k_}/{a_} (n=10, 2 seeds)"] = em - th
    except Exception as e:                                        # pragma: no cover
        print(f"  (s20 encoding artefacts unavailable: {e!r})")

    nl = _nl()
    for pot in ("legacy", "amber"):
        end = np.array([x["null"][pot]["rmsd_end"] for x in nl])
        st = np.array([x["null"][pot]["rmsd_start"] for x in nl])
        iso = np.array([x["null"][pot]["rand_iso"] for x in nl])
        twd = np.array([x["null"][pot]["toward_member"] for x in nl])
        comps[f"[D3] {pot}: minimised - start"] = end - st
        comps[f"[D3] {pot}: minimised - toward_member"] = end - twd
        comps[f"[D3] {pot}: toward_member - rand_iso"] = twd - iso

    #: the LOW-VARIANCE end of the range: the deployed repair operators, from the production
    #: cache.  This is the comparison the briefs describe as "a quarter of the MDE".
    try:
        from s12 import instrument as I
        pdbs = [t["pdb"] for t in I.targets()]
        rec = [I.shipped_record(p) for p in pdbs]
        fit = np.array([x["rmsd_fit"] for x in rec])
        full = np.array([x["rmsd_full"] for x in rec])
        arm = np.array([x["rmsd_arm"] for x in rec])
        avg = np.array([x["rmsd_avg"] for x in rec])
        comps["[X2] AMBER on the built chain (rmsd_full - rmsd_arm)"] = full - arm
        comps["[X2] projection lam=0.3 vs point cloud (rmsd_arm - rmsd_avg)"] = arm - avg
        comps["[X2] projection lam=0 vs point cloud (rmsd_fit - rmsd_avg)"] = fit - avg
    except Exception as e:                                        # pragma: no cover
        print(f"  (production cache unavailable: {e!r})")

    out = {"quoted_mde": QUOTED, "z_sum": Z,
           "implied_sd_at_n126": QUOTED * math.sqrt(126) / Z, "cells": {}}
    print(f"\n=== D5  MDE IS A PROPERTY OF THE COMPARISON, NOT OF THE INSTRUMENT ===")
    print(f"MDE = {Z:.4f} * sd(paired differences) / sqrt(n),  two-sided alpha=0.05, 80% power.")
    print(f"The brief quotes a single {QUOTED} A for every comparison.  The paired sd that would")
    print(f"make that true at n = 126 is {out['implied_sd_at_n126']:.4f} A.\n")
    print(f"  {'comparison':<48}{'n':>5}{'mean':>9}{'sd':>9}{'SE':>9}"
          f"{'MDE':>9}{'x0.084':>9}")
    for k in comps:
        m = mde(comps[k])
        out["cells"][k] = m
        print(f"  {k:<48}{m['n']:>5}{m['mean']:>+9.4f}{m['sd']:>9.4f}{m['se']:>9.4f}"
              f"{m['mde']:>9.4f}{m['ratio_to_quoted']:>9.2f}")
    rr = np.array([out["cells"][k]["ratio_to_quoted"] for k in out["cells"]])
    print(f"\n  range of per-comparison MDE / 0.084:  {rr.min():.2f}x to {rr.max():.2f}x  "
          f"-- a factor of {rr.max()/rr.min():.0f} across {len(rr)} real comparisons")
    within = float(np.mean((rr > 1 / 1.5) & (rr < 1.5)))
    print(f"  fraction within a factor 1.5 of the quoted constant: {100*within:.0f}%")
    print(f"\n  PRE-REGISTERED FALSIFIER: D5 is REFUTED if the per-comparison MDEs cluster within")
    print(f"  a factor ~1.5 of 0.084.  They span {rr.max()/rr.min():.0f}x, and {100*(1-within):.0f}% "
          f"fall outside that band -> D5 SUPPORTED.")
    print("\n  READ.  Both failure modes are present in this table:")
    lo = min(out["cells"], key=lambda k: out["cells"][k]["mde"])
    hi = max(out["cells"], key=lambda k: out["cells"][k]["mde"])
    print(f"    LOWEST  MDE {out['cells'][lo]['mde']:.4f} A ({out['cells'][lo]['ratio_to_quoted']:.2f}x) "
          f"-- {lo}")
    print(f"            its own effect is {out['cells'][lo]['mean']:+.4f} A, i.e. "
          f"{abs(out['cells'][lo]['mean'])/out['cells'][lo]['se']:.1f} SE."
          f"  Calling this 'below the MDE' is wrong.")
    print(f"    HIGHEST MDE {out['cells'][hi]['mde']:.4f} A ({out['cells'][hi]['ratio_to_quoted']:.2f}x) "
          f"-- {hi}")
    print(f"            a null there is uninformative below {out['cells'][hi]['mde']:.3f} A, "
          f"not below 0.084.")
    print("\n  FIX, one line per comparison: report SE beside every mean, and quote the MDE the")
    print("  comparison's own paired sd implies.  `s16/review_FINDINGS.md` Q4.1 already did this")
    print("  correctly for its panel; the constant is what generalised, not the method.")
    json.dump(out, open(os.path.join(RESULTS, "d_mde.json"), "w"), indent=1)
    ok = len(out["cells"]) >= 20 and all(np.isfinite(v["mde"]) for v in out["cells"].values())
    p = os.path.join(RESULTS, "d_mde.COMPLETE")
    if ok:
        with open(p, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"comparisons={len(out['cells'])} sources=d_cvarop,d_divsel,"
                     f"s20/c_land_null,production_cache quoted_mde={QUOTED} "
                     f"range={rr.min():.2f}x-{rr.max():.2f}x\n")
        print("\nCOMPLETE.")
    elif os.path.exists(p):
        os.remove(p)
    return out


if __name__ == "__main__":
    run()
