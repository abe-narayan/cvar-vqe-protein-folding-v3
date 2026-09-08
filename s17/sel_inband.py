"""s17/sel_inband.py -- DOES ANY NATIVE-FREE SIGNAL DISCRIMINATE INSIDE THE BAND?

THE REFRAMING THIS MODULE SERVES.  `s17/sel_obj.py` measured, at K = 500 on n = 126, that the
shipped distance objective's rank correlation with the truth is **0.568 globally and 0.131
inside the near-native band**.  Its skill is almost entirely garbage-versus-plausible and
nearly absent where the argmin is actually decided.  That is why 58 re-weightings of the same
distogram returned -0.009 A [-0.128, +0.110] leave-fold-out: **you cannot fix an in-band
problem by re-weighting a global signal.**

So the question is no longer "does it rank".  It is:

    does ANY available native-free signal have in-band discrimination above chance?

and it is asked of everything the programme has, on one instrument, with a chance null and
zero-information references attached.

TWO BAND DEFINITIONS, BOTH REPORTED, BECAUSE THEY ARE DIFFERENT OBJECTS.

    band_oracle    members within best + 1.5 A.  ORACLE-DEFINED: it needs the native to
                   draw, so it is a DIAGNOSTIC band, never a deployable filter.  It is what
                   the argmin decision actually spans, and it is the definition the
                   programme's recorded 0.600 / 0.638 numbers use.
    band_topB      the shipped objective's own top B (B = 25, 75).  NATIVE-FREE: this is what
                   a deployed shortlist actually is, and it is the set a reranker would see.

A signal's in-band number is quoted with the band it was measured in, always.

THE THREE NULLS, all mandatory (BRIEF section 4).

    chance         a seeded random score.  Its measured spread is the null band, printed --
                   not assumed to be zero, because the per-target rho of a random score at
                   band sizes of 20-80 members has a standard deviation of order 0.2.
    alpha / beta   ZERO-INFORMATION references: CA-RMSD of each candidate to a constant ideal
                   alpha helix and to a constant ideal beta strand.  Sprint 16's most-repeated
                   lesson is that these reproduce more than anyone expects.
    randwin        CA-RMSD of each candidate to a RANDOM member of the same pool -- a
                   zero-information "consensus" that knows nothing about which member is good.

UNITS.  Three scales appear in this programme and they are NOT interchangeable:
    pairwise ordering accuracy (null 0.500)  <->  Kendall tau = 2*acc - 1
    Spearman rho_S
    Gaussian-copula rho = sin(pi*tau/2) from accuracy, or 2*sin(pi*rho_S/6) from Spearman
The recorded "in-band ordering 0.600 across targets" is an ACCURACY; as a copula rho it is
0.309.  Every table here prints Spearman and the accuracy equivalent side by side.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from s17 import sel_lib as L             # noqa: E402
from s17 import sel_obj as O             # noqa: E402

K = 500
BANDS = ("oracle1.5", "top25", "top75")
LEG_TERMS = ("steric", "contact", "hbond_local", "hbond_longrange", "coop_helix",
             "coop_sheet", "solvation", "electrostatic", "aromatic", "torsion",
             "compactness")

#: ideal-geometry reference conformations, in radians
ALPHA = (np.deg2rad(-57.0), np.deg2rad(-47.0))
BETA = (np.deg2rad(-139.0), np.deg2rad(135.0))


def acc_from_spearman(rs):
    """Spearman -> Gaussian-copula rho -> pairwise ordering accuracy (null 0.500)."""
    rho = 2.0 * np.sin(np.pi * np.asarray(rs, float) / 6.0)
    tau = (2.0 / np.pi) * np.arcsin(np.clip(rho, -1, 1))
    return 0.5 * (1.0 + tau)


def signals(p, u, sc, W, leg):
    """Every native-free per-candidate signal, lower = predicted better."""
    n, k = p["n"], p["k"]
    D = np.asarray(p["D"], float)
    s = {}
    s["dist_shipped"] = sc
    c = O.context(p)
    s["dist_rankset"] = O._agg(c["risk"], c["w"]["ship"], "rankset")
    s["dist_blend15"] = O.blend_scores(c, betas=(0.15,), scale="pred")["blendpred_0.15"]
    s["dist_z"] = O._agg(c["z"], c["w"]["one"], "mean")

    # --- retrieval channel
    sim = np.asarray(u["sim"], float)[np.asarray(u["order"], int)][:k]
    s["retrieval_rank"] = np.arange(k, dtype=float)
    s["blosum_sim"] = -sim
    s["not_peptide_db"] = (~np.asarray(u["org"], bool)[np.asarray(u["order"], int)][:k]).astype(float)

    # --- consensus / typicality channel (candidate-set-conditional, no model)
    s["typicality"] = np.abs(D - D.mean(0)[None, :]).mean(1)
    s["typicality_med"] = np.abs(D - np.median(D, axis=0)[None, :]).mean(1)

    # --- geometry / compactness channel
    Rg = np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
    exp = p["expected"]
    s2 = float((exp ** 2).sum() * 2) + 2 * (n - 1) * 3.8 ** 2
    rg_pred = float(np.sqrt(max(s2 / (2 * n * n), 1e-6)))
    s["rg_raw"] = Rg
    s["rg_vs_predicted"] = np.abs(Rg - rg_pred)
    cmap = (D < 8.0)
    s["n_contacts"] = -cmap.sum(1).astype(float)
    sep = p["sep"]
    s["contact_order"] = np.where(cmap.sum(1) > 0,
                                  (cmap * sep[None, :]).sum(1) / np.maximum(cmap.sum(1), 1), 0.0)

    # --- Legacy channel
    for t in LEG_TERMS:
        s[f"leg_{t}"] = np.asarray(leg[t], float)[:k]
    s["leg_total"] = np.asarray(leg["_total"], float)[:k]

    # --- NULLS
    rng = SD.stable_rng(p["pdb"], "s17inband")
    s["NULL_chance"] = rng.standard_normal(k)
    for nm, (ph, ps) in (("NULL_alpha", ALPHA), ("NULL_beta", BETA)):
        ref = I.build_ca(np.full(n, ph), np.full(n, ps))
        s[nm] = I.kabsch_rmsd_batch(W, ref)
    s["NULL_randwin"] = I.kabsch_rmsd_batch(W, W[int(rng.integers(k))])
    return s


def one_target(t):
    p = L.pack(t["pdb"], K, want=("D", "W", "T", "org"))
    u = I.load_univ(t["pdb"])
    sc = L.shipped(p)
    rr = p["rr"]
    W = p["W"]
    from s16 import energy_lib as EL
    comp = EL.legacy_components_of_windows(t["seq"], p["PHI"], p["PSI"])
    leg = {q: np.asarray(comp[q], float) for q in LEG_TERMS}
    leg["_total"] = np.asarray(EL.legacy_total_from(comp), float)

    sg = signals(p, u, sc, W, leg)
    masks = {
        "oracle1.5": rr <= rr.min() + 1.5,
        "top25": np.isin(np.arange(len(rr)), np.argsort(sc, kind="stable")[:25]),
        "top75": np.isin(np.arange(len(rr)), np.argsort(sc, kind="stable")[:75]),
    }
    out = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "k": int(p["k"]),
           "oracle": float(rr.min()), "rho": {}, "sel": {}, "band": {}}
    for b, m in masks.items():
        out["band"][b] = {"size": int(m.sum()), "best": float(rr[m].min()),
                          "mean": float(rr[m].mean())}
        out["rho"][b] = {}
        out["sel"][b] = {}
        if m.sum() < 8:
            continue
        for nm, v in sg.items():
            out["rho"][b][nm] = O.spearman(v[m], rr[m])
            out["sel"][b][nm] = L.sel_of(v[m], rr[m])
    return out


def run(targets=None, verbose=True):
    tg = targets if targets is not None else L.targets()
    rows, t0 = [], time.time()
    for q, t in enumerate(tg):
        rows.append(one_target(t))
        if verbose and (q + 1) % 20 == 0:
            print(f"  {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(os.path.join(L.RESULTS, "sel_inband.json"), "w"))
    json.dump({"rows": rows}, open(os.path.join(L.RESULTS, "sel_inband.json"), "w"))
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(L.RESULTS, "sel_inband.json")))["rows"]
    n = len(rows)
    print(f"\n{'='*108}\nIN-BAND DISCRIMINATION SURVEY   n = {n} targets, K = {K}, "
          f"identical candidate set\n{'='*108}")
    print("Spearman rho of each native-free signal against TRUE RMSD, inside the band.")
    print("acc = the same on the pairwise-ordering-accuracy scale (null 0.500), which is the")
    print("scale the programme's recorded 0.600 / 0.638 numbers live on.  They are NOT the")
    print("same number: accuracy 0.600 is Spearman 0.288 is copula rho 0.309.\n")

    for b in BANDS:
        avail = [r for r in rows if r["rho"].get(b)]
        if not avail:
            continue
        keys = sorted(avail[0]["rho"][b].keys())
        sz = np.array([r["band"][b]["size"] for r in avail], float)
        bb = np.array([r["band"][b]["best"] for r in avail], float)
        bm = np.array([r["band"][b]["mean"] for r in avail], float)
        lab = ("ORACLE-DEFINED (diagnostic only)" if b.startswith("oracle")
               else "NATIVE-FREE (a deployable shortlist)")
        print(f"\n--- BAND = {b}   {lab}   n = {len(avail)} targets")
        print(f"    band size {sz.mean():.1f}   ORACLE best in band {bb.mean():.3f}   "
              f"band mean (= random-in-band) {bm.mean():.3f}")
        chance = np.array([r["rho"][b]["NULL_chance"] for r in avail], float)
        print(f"    chance null: mean rho {chance.mean():+.3f}, per-target sd "
              f"{chance.std():.3f}, so the null band on the MEAN is about "
              f"+-{1.96*chance.std()/np.sqrt(len(avail)):.3f}")
        print(f"    {'signal':<22}{'rho':>8}{'acc':>7}{'95% CI on rho':>20}"
              f"{'selected':>10}{'vs band-random':>16}{'W/L':>9}")
        fold = np.array([r["fold"] for r in avail], int)
        res = []
        for kx in keys:
            v = np.array([r["rho"][b][kx] for r in avail], float)
            s = np.array([r["sel"][b][kx] for r in avail], float)
            st = L.boot_target(v)
            pr = L.report_pair(kx, s, bm, fold)
            res.append((abs(v.mean()), kx, v.mean(), st, s.mean(), pr))
        for _a, kx, mu, st, sm, pr in sorted(res, key=lambda z: -z[0]):
            star = "*" if (st["lo"] > 0 or st["hi"] < 0) else " "
            print(f"   {star}{kx:<22}{mu:+8.3f}{acc_from_spearman(mu):>7.3f}"
                  f"   [{st['lo']:+.3f},{st['hi']:+.3f}]{sm:>10.3f}"
                  f"   {pr['diff']:+.3f} [{pr['lo']:+.3f},{pr['hi']:+.3f}]{pr['W']:>5}/{pr['L']}")
        print("    * = the target-level CI on rho excludes zero.  A signal is only useful if")
        print("      its SELECTED column also beats band-random with a CI excluding zero.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
        report()
