#!/usr/bin/env python
"""s32/s32_D_reemit.py -- S32 lane D: RE-EMIT the five artefacts that had no producing script.

Registered in s32/PREREG_S32_D.md, commit 34973b1b.

## Why this file exists -- AUDIT 11, fourth instance sprint-wide, second in this lane

Five lane-D artefacts were produced by INLINE `python -c` commands.  They carry no `provenance`
block, no module, no git commit, no source hash and no pinned seed, and no `.py` in the repository
writes their filenames -- while `s32/REPORT_S32.md` section 4 and `s32/LEDGER.md` quote numbers from
at least three of them.  This is the same defect this lane already fixed once on D1-T
(`s32_D1_signtransfer.py`), and the same defect lane R then found twice in its own lane.

**Two of the five are worse than "no provenance": their HEADLINE NUMBERS ARE NOT IN THE FILE.**
`s32_D0X_circularity.json` and `s32_D4X_rgsign.json` stored only per-target rows; the variance
retention table and the "positive on 81% of targets" figure existed only in stdout.  This file
recomputes and STORES them.

## What is re-emitted, and the load-bearing number in each

  s32_D1_signrandom.json   the within-band LABEL-PERMUTATION NULL.  This is the control that turns
                           E[rho] = 0 into Var(rho) > 0 -- the evidence for the sprint's main
                           reframing.  REPORT section 4 quotes its table.
  s32_D1_signshare.json    cross-scorer in-band SIGN agreement (AMBER vs Legacy 36.5%, z = -3.0)
                           and the ORACLE sign-corrected price mean|rho|.
  s32_D0X_circularity.json the variance-retention table that refuted the cross-lane chirality
                           synthesis, with the size-matched random-75 control and the EXACT
                           achiral twin.  Table was not in the file.
  s32_D4X_rgsign.json      sign(rg_pred - rg_pool) as the per-target bit.  "Positive on 81% of
                           targets" -- the named mechanism for the failure -- was not in the file.
  s32_D5_signprice.json    the MARGINAL a one-bit classifier must beat, per-fold accuracy, and the
                           CLOUD-basis Angstrom screen.  REPORT section 4 quotes the ORACLE-bit price.

## Reproduction discipline

Seeds are pinned to the ORIGINAL inline values and the RNG draw ORDER is reproduced exactly, so the
numbers should come back bit-for-bit.  Every emitted file carries a `reemit_check` block giving the
max absolute difference against the previous artefact on every numeric leaf they share.  **If any
quoted number moves, this script says so loudly and the report must be corrected rather than quietly
keeping the old value.**

ORACLE: `rr` and `nat_ca` are read as EVALUATION LABELS only.  No arm scored here reads them.

    python s32/s32_D_reemit.py
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
    os.environ.setdefault(_v, "1")

from s12 import instrument as I          # noqa: E402
from s24 import stats_lib as ST          # noqa: E402
from s16 import energy_lib as EL         # noqa: E402
from s32.s32_D1_inband import (          # noqa: E402
    spear, spear_partial, rg, ca_pseudo_torsion, _rank, _pear)
from s32.s32_D4_signpred import features_and_rho, lfo_ridge   # noqa: E402

RES = os.path.join(HERE, "results")

# Seeds pinned to the ORIGINAL inline values so the streams reproduce exactly.
SEED_SIGNRANDOM = 32004        # s32_D1_signrandom.json
SEED_D0X = 320041              # s32_D0X_circularity.json
NPERM = 24                     # permutation draws per target per scorer
SC4 = ("AMBER", "DIS", "LEG_total", "LEG_torsion")   # ORDER IS LOAD-BEARING for the RNG stream
HYDRO = set("AVLIMFWCY")


# --------------------------------------------------------------------------- data, one pass
def collect():
    """One pass over the 126 targets.  Everything the five artefacts need."""
    out = []
    for t in I.targets():
        pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
        u = I.load_univ(pdb)
        p = I.pool_idx(u, k=500)
        sub = np.asarray(I.shipped_record(pdb)["sub"], int)
        W = np.asarray(u["W"][p], float)                      # (500, n, 3)
        z = np.load(os.path.join(ROOT, "s24", "cache_amber", pdb + ".npz"))
        assert np.array_equal(np.asarray(z["universe_idx"], int), p), pdb
        comp = EL.legacy_components_of_windows(
            seq, np.asarray(u["PHI"][p], float), np.asarray(u["PSI"][p], float))
        dg = I.distogram(pdb)
        exp = np.asarray(dg["expected"], float)
        out.append({
            "pdb": pdb, "seq": seq, "n": n, "fold": int(t["fold"]),
            "sub": sub, "W": W, "nat": np.asarray(u["nat_ca"], float),
            "rr500": np.asarray(u["rr"][p], float),           # ORACLE label
            "e_amber": np.asarray(z["e_amber"], float),
            "score_dist": np.asarray(z["score_dist"], float),
            "leg": EL.legacy_total_from(comp),
            "leg_torsion": np.asarray(comp["torsion"], float),
            "rg500": rg(W),
            "exp": exp,
            "rg_pred": float(np.sqrt((np.sum(exp ** 2) * 2 + (n - 1) * 2 * 3.8 ** 2)
                                     / (2.0 * n * n)))})
    return out


def band_scores(d):
    """The four scorers on the shipped top-75 band, in the RNG-load-bearing order."""
    s = d["sub"]
    return {"AMBER": d["e_amber"][s], "DIS": d["score_dist"][s],
            "LEG_total": d["leg"][s], "LEG_torsion": d["leg_torsion"][s]}


# --------------------------------------------------------------------------- reemit check
def _leaves(o, pre=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from _leaves(v, pre + "/" + str(k))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from _leaves(v, pre + "/%d" % i)
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        yield pre, float(o)


def reemit_check(path, new):
    """Compare every shared numeric leaf against the previous artefact."""
    if not os.path.exists(path):
        return {"previous_artefact": "ABSENT", "max_abs_diff": None}
    old = json.load(open(path, encoding="utf-8"))
    a, b = dict(_leaves(old)), dict(_leaves(new))
    shared = sorted(set(a) & set(b))
    diffs = [(k, abs(a[k] - b[k])) for k in shared]
    worst = max(diffs, key=lambda x: x[1]) if diffs else (None, 0.0)
    big = sorted([(k, a[k], b[k]) for k, dv in diffs if dv > 1e-9],
                 key=lambda x: -abs(x[1] - x[2]))[:20]
    return {"previous_artefact": os.path.basename(path),
            "n_shared_numeric_leaves": len(shared),
            "max_abs_diff": worst[1], "max_abs_diff_key": worst[0],
            "n_leaves_differing_above_1e-9": len(big),
            "differing": [{"key": k, "old": o, "new": nn} for k, o, nn in big],
            "REPRODUCES": bool(worst[1] <= 1e-9)}


def write(name, obj):
    path = os.path.join(RES, name)
    obj["provenance"] = ST.provenance(__file__)
    obj["prereg_commit"] = "34973b1b"
    obj["reemit_check"] = reemit_check(path, obj)
    ST.save_atomic(path, obj, module_file=__file__)
    rc = obj["reemit_check"]
    flag = "REPRODUCES" if rc.get("REPRODUCES") else "*** CHANGED ***"
    print("  %-30s %-16s max|diff| %s over %s shared leaves" % (
        name, flag, ("%.3e" % rc["max_abs_diff"]) if rc["max_abs_diff"] is not None else "n/a",
        rc.get("n_shared_numeric_leaves")))
    for d in rc.get("differing", [])[:6]:
        print("       %-52s %.6f -> %.6f" % (d["key"][:52], d["old"], d["new"]))
    return obj


# --------------------------------------------------------------------------- the five
def emit_signrandom(data, pdbs, folds):
    """Within-band label-permutation null.  RNG order: per target, per scorer, NPERM draws."""
    rng = np.random.default_rng(SEED_SIGNRANDOM)
    acc = {nm: {"obs": [], "null": []} for nm in SC4}
    for d in data:
        rr = d["rr500"][d["sub"]]
        S = band_scores(d)
        for nm in SC4:                                   # insertion order of the original dict
            v = S[nm]
            acc[nm]["obs"].append(abs(spear(v, rr)))
            acc[nm]["null"].append(
                np.mean([abs(spear(v, rng.permutation(rr))) for _ in range(NPERM)]))
    out = {}
    for nm in SC4:
        o = np.array(acc[nm]["obs"]); nl = np.array(acc[nm]["null"])
        c = ST.compare(o, nl, folds=folds, names=pdbs, label=nm)
        out[nm] = {"obs_mean": float(o.mean()), "null_mean": float(nl.mean()),
                   "ratio": float(o.mean() / nl.mean()),
                   "effect": c["effect"], "se": c["se"], "mde": c["mde"],
                   "x_mde": abs(c["effect"]) / c["mde"], "folds_same_sign": c["folds_same_sign"],
                   "ci95_fold": c["ci95_fold"]}
    out["_meta"] = {
        "seed": SEED_SIGNRANDOM, "n_perm_draws_per_target_per_scorer": NPERM,
        "scorer_order_is_rng_load_bearing": list(SC4),
        "basis": ("per-target |in-band Spearman| against true CA-RMSD inside the shipped top-75 "
                  "band, n=126; the null holds the SCORE VECTOR and BAND SIZE fixed and permutes "
                  "the ORACLE rr labels inside each target's own band"),
        "what_it_establishes": ("Var(rho) > 0 with E[rho] = 0 -- in-band ordering content exists "
                                "on every scorer including the one with exactly zero mean skill"),
        "oracle": "ORACLE rr is an evaluation label only; no arm scored here reads it"}
    return write("s32_D1_signrandom.json", out)


def emit_signshare(data, pdbs, folds):
    """Cross-scorer in-band sign agreement.  Fully deterministic."""
    nm5 = ["AMBER", "DIS", "LEG_total", "LEG_torsion", "RG"]
    R = {k: [] for k in nm5}
    for d in data:
        rr = d["rr500"][d["sub"]]
        S = band_scores(d)
        S["RG"] = d["rg500"][d["sub"]]
        for k in nm5:
            R[k].append(spear(S[k], rr))
    R = {k: np.array(v, float) for k, v in R.items()}
    se = np.sqrt(0.25 / len(pdbs))
    o = np.sign(R["LEG_total"]) * R["AMBER"]
    o2 = np.sign(R["DIS"]) * R["AMBER"]
    c = ST.compare(o, np.zeros(len(pdbs)), folds=folds, names=pdbs, label="AMBER|LEG sign")
    c2 = ST.compare(o2, np.zeros(len(pdbs)), folds=folds, names=pdbs, label="AMBER|DIS sign")
    out = {"note": ("D1-S cross-scorer in-band sign agreement; |rho| rows are ORACLE / NOT "
                    "DEPLOYABLE price"),
           "corr": {a: {b: float(np.corrcoef(R[a], R[b])[0, 1]) for b in nm5} for a in nm5},
           "sign_agree": {a + "|" + b: float((np.sign(R[a]) == np.sign(R[b])).mean())
                          for i, a in enumerate(nm5) for b in nm5[i + 1:]},
           "sign_agree_z_from_chance": {
               a + "|" + b: float(((np.sign(R[a]) == np.sign(R[b])).mean() - 0.5) / se)
               for i, a in enumerate(nm5) for b in nm5[i + 1:]},
           "oracle_sign_corrected_mean_abs_rho": {a: float(np.abs(R[a]).mean()) for a in nm5},
           "amber_oriented_by_leg": {"effect": c["effect"], "se": c["se"],
                                     "x_mde": abs(c["effect"]) / c["mde"],
                                     "folds": c["folds_same_sign"]},
           "amber_oriented_by_dis": {"effect": c2["effect"], "se": c2["se"],
                                     "x_mde": abs(c2["effect"]) / c2["mde"],
                                     "folds": c2["folds_same_sign"]},
           "_meta": {"deterministic": True, "n_pairwise_tests": 10,
                     "bonferroni_bar_abs_z": 2.81,
                     "limitation": ("sign agreement is dominated by targets where |rho| sits in "
                                    "the noise, and the two rho's share a common referent (rr), "
                                    "so memory `shared-referent-floor` applies to their "
                                    "correlation; the ANTI-agreements are the robust half")}}
    return write("s32_D1_signshare.json", out)


def emit_d0x(data, pdbs, folds):
    """Variance retention under the shipped top-75 filter, with a size-matched random-75 control
    and the EXACT achiral twin.  RNG: one rng.choice(500, 75) per target."""
    rng = np.random.default_rng(SEED_D0X)
    keys = ("CHI_odd_sin", "CHI_even_cos", "RG", "E2E", "DIS", "LEG", "RR_ORACLE")
    rows = []
    for d in data:
        W, sub, rr = d["W"], d["sub"], d["rr500"]
        tau = ca_pseudo_torsion(W)
        ch = {"CHI_odd_sin": np.sin(tau).mean(1),          # CHIRAL (odd under reflection)
              "CHI_even_cos": np.cos(tau).mean(1),         # the EXACT achiral twin (even)
              "RG": d["rg500"],
              "E2E": np.linalg.norm(W[:, -1] - W[:, 0], axis=-1),
              "DIS": d["score_dist"], "LEG": d["leg"], "RR_ORACLE": rr}
        r = {"pdb": d["pdb"], "n": d["n"], "fold": d["fold"]}
        ridx = rng.choice(500, 75, replace=False)
        for k in keys:
            v = ch[k]
            s5 = np.std(v)
            r["ret_" + k] = (np.std(v[sub]) / s5) ** 2 if s5 > 1e-12 else np.nan
            r["retrand_" + k] = (np.std(v[ridx]) / s5) ** 2 if s5 > 1e-12 else np.nan
        r["rho_chi"] = spear(ch["CHI_odd_sin"], rr)
        r["rho_chi_even"] = spear(ch["CHI_even_cos"], rr)
        r["rho_chi_pDIS"] = spear_partial(ch["CHI_odd_sin"], rr, ch["DIS"])
        r["rho_chi_pRG"] = spear_partial(ch["CHI_odd_sin"], rr, ch["RG"])
        r["rho_DIS"] = spear(ch["DIS"], rr)
        r["rho_DIS_pCHI"] = spear_partial(ch["DIS"], rr, ch["CHI_odd_sin"])
        r["corr_chi_dis"] = _pear(_rank(ch["CHI_odd_sin"]), _rank(ch["DIS"]))
        r["corr_chi_rg"] = _pear(_rank(ch["CHI_odd_sin"]), _rank(ch["RG"]))
        rows.append(r)
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    zero = np.zeros(len(rows))
    c = ST.compare(g("ret_CHI_odd_sin"), g("ret_CHI_even_cos"), folds=folds, names=pdbs,
                   label="chiral vs its achiral twin")
    agg = {"retention_score75": {k: float(np.nanmean(g("ret_" + k))) for k in keys},
           "retention_random75_control": {k: float(np.nanmean(g("retrand_" + k))) for k in keys},
           "chiral_minus_achiral_twin_retention": {
               "effect": c["effect"], "se": c["se"], "mde": c["mde"],
               "x_mde": abs(c["effect"]) / c["mde"], "folds_same_sign": c["folds_same_sign"],
               "n_better": c["n_better"], "n_worse": c["n_worse"]},
           "global_rho": {}}
    for k, lab in (("rho_chi", "chiral"), ("rho_chi_even", "achiral_twin"),
                   ("rho_chi_pDIS", "chiral_given_DIS"), ("rho_chi_pRG", "chiral_given_RG"),
                   ("rho_DIS", "DIS"), ("rho_DIS_pCHI", "DIS_given_chiral")):
        cc = ST.compare(g(k), zero, folds=folds, names=pdbs, label=lab)
        agg["global_rho"][lab] = {"mean": cc["effect"], "se": cc["se"],
                                  "x_mde": abs(cc["effect"]) / cc["mde"],
                                  "folds_same_sign": cc["folds_same_sign"]}
    agg["rank_corr_chiral_with"] = {"DIS": float(g("corr_chi_dis").mean()),
                                    "RG": float(g("corr_chi_rg").mean())}
    out = {"note": ("D0-X circularity control: variance retention of a CHIRAL coordinate vs its "
                    "EXACT achiral twin (mean cos tau) and achiral controls, under the shipped "
                    "top-75 filter and a size-matched random-75 control.  ORACLE rr used as a "
                    "label only."),
           "AGGREGATES": agg, "rows": rows,
           "_meta": {"seed": SEED_D0X,
                     "finding": ("the filter shrinks EVERYTHING; chirality is shrunk LESS than "
                                 "its own achiral twin; and the chiral coordinate's retention "
                                 "(0.271) equals lane V's ORACLE-quality retention (0.272) -- "
                                 "ONE FACT, NOT TWO.  The cross-lane synthesis is circular."),
                     "aggregates_were_missing": ("the original artefact stored only `rows`; the "
                                                 "retention table existed only in stdout")}}
    return write("s32_D0X_circularity.json", out)


def emit_d4x(data, pdbs, folds):
    """sign(rg_pred - rg_pool_mean) as the per-target bit.  Fully deterministic."""
    rows = []
    for d in data:
        sub = d["sub"]
        rr = d["rr500"][sub]
        S = band_scores(d)
        G = d["rg500"][sub]
        rows.append({"pdb": d["pdb"], "fold": d["fold"],
                     "rho_LEG": spear(S["LEG_total"], rr), "rho_AMB": spear(S["AMBER"], rr),
                     "rho_RG": spear(G, rr),
                     "rg_dis": d["rg_pred"] - float(G.mean()),
                     "rg_dis_med": d["rg_pred"] - float(np.median(G))})
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    sd = np.sign(g("rg_dis"))
    agg = {"frac_rg_disagree_positive": float((g("rg_dis") > 0).mean()),
           "corr_rgdisagree_rhoRG": float(np.corrcoef(g("rg_dis"), g("rho_RG"))[0, 1]),
           "arms": {}}
    for k, lab in (("rho_RG", "RG"), ("rho_LEG", "LEG_total"), ("rho_AMB", "AMBER")):
        y = g(k)
        p = float((y > 0).mean())
        c = ST.compare(sd * y, y, folds=folds, names=pdbs, label=lab)
        agg["arms"][lab] = {
            "marginal": float(max(p, 1 - p)),
            "accuracy_of_sign_rg_disagree": float((sd == np.sign(y)).mean()),
            "accuracy_of_minus_sign": float((-sd == np.sign(y)).mean()),
            "oriented_rho": float((sd * y).mean()),
            "constant_plus_one_control_rho": float(y.mean()),
            "delta_vs_own_constant": c["effect"], "se": c["se"],
            "x_mde": (abs(c["effect"]) / c["mde"]) if c["mde"] else 0.0,
            "folds_same_sign": c["folds_same_sign"]}
    out = {"note": ("D4-X mechanism-motivated single native-free feature "
                    "sign(rg_pred - rg_pool_mean) as the per-target sign; rho columns are ORACLE "
                    "labels used for scoring only"),
           "AGGREGATES": agg, "rows": rows,
           "_meta": {"deterministic": True,
                     "mechanism_of_failure": ("rg_pred - rg_pool_mean is POSITIVE on 81.0% of "
                                              "targets, so sign(rg_disagree) is nearly a constant "
                                              "+1 and its apparent accuracy is the target's "
                                              "marginal wearing a different label.  A one-bit "
                                              "feature that is 81/19 cannot carry a 59/41 label."),
                     "aggregates_were_missing": ("the original artefact stored only `rows`; the "
                                                 "81% and the accuracies existed only in stdout")}}
    return write("s32_D4X_rgsign.json", out)


def emit_d5_signprice(data, pdbs, folds):
    """The MARGINAL a one-bit classifier must beat, per-fold LFO accuracy, and the CLOUD-basis
    Angstrom screen.  Deterministic: lfo_ridge and argsort carry no RNG."""
    rows = features_and_rho()
    fold = np.array([r["fold"] for r in rows], int)
    names = sorted(rows[0]["feat"])
    X = np.array([[r["feat"][k] for k in names] for r in rows], float)
    acc = {}
    for sc in SC4:
        y = np.array([r["rho_" + sc] for r in rows], float)
        yh = lfo_ridge(X, y, fold)
        p = float((y > 0).mean())
        acc[sc] = {"marginal": float(max(p, 1 - p)), "p_pos": p,
                   "lfo_acc": float((np.sign(yh) == np.sign(y)).mean()),
                   "per_fold": [float((np.sign(yh[fold == f]) == np.sign(y[fold == f])).mean())
                                for f in sorted(set(fold))]}
    # --- CLOUD-basis screen (a SCREEN; never differenced against a built-chain number)
    res = {}
    for d in data:
        sub, W, nat = d["sub"], d["W"][d["sub"]], d["nat"]
        rr = d["rr500"][sub]
        S = {"AMBER": d["e_amber"][sub], "LEG_total": d["leg"][sub],
             "DIS": d["score_dist"][sub]}
        C, _ = I.coordinate_average(W)
        res.setdefault("PROD", []).append(I.ca_rmsd(C, nat))
        for sc, v in S.items():
            rho = spear(v, rr)
            for m in (25, 50):
                for lab, sg in (("const", 1.0), ("ORACLE", np.sign(rho) if rho != 0 else 1.0)):
                    idx = np.argsort(sg * v, kind="mergesort")[:m]
                    Cm, _ = I.coordinate_average(W[idx])
                    res.setdefault("%s_m%d_%s" % (sc, m, lab), []).append(I.ca_rmsd(Cm, nat))
    prod = np.array(res["PROD"])
    screen = {"PROD_mean": float(prod.mean())}
    for k in sorted(res):
        if k == "PROD":
            continue
        v = np.array(res[k])
        c = ST.compare(v, prod, folds=folds, names=pdbs, label=k)
        screen[k] = {"mean": float(v.mean()), "delta_vs_prod": c["effect"], "se": c["se"],
                     "x_mde": (abs(c["effect"]) / c["mde"]) if c["mde"] else 0.0,
                     "folds_same_sign": c["folds_same_sign"], "ci95_fold": c["ci95_fold"]}
    out = {"marginal_and_accuracy": acc, "cloud_screen": screen,
           "note": ("CLOUD basis (CA point cloud), a SCREEN for whether an ORACLE per-target sign "
                    "is worth projecting.  NEVER differenced against a built-chain number.  "
                    "ORACLE arms are ORACLE / NOT DEPLOYABLE."),
           "_meta": {"deterministic": True,
                     "the_baseline_that_binds": ("a one-bit classifier must beat the MARGINAL "
                                                 "max(p, 1-p), not 0.50"),
                     "built_chain_followup": ("s32/s32_D5_signchain.py prices the bit on the "
                                              "BUILT CHAIN at -0.2101 A, 1.66x MDE, 5/5 folds, "
                                              "ORACLE / NOT DEPLOYABLE")}}
    return write("s32_D5_signprice.json", out)


def main():
    print("collecting (one pass over 126 targets) ...", flush=True)
    data = collect()
    pdbs = [d["pdb"] for d in data]
    folds = ST.pinned_folds(pdbs)
    assert np.array_equal(folds, np.array([d["fold"] for d in data], int))
    print("re-emitting five artefacts with provenance and pinned seeds:")
    emit_signrandom(data, pdbs, folds)
    emit_signshare(data, pdbs, folds)
    emit_d0x(data, pdbs, folds)
    emit_d4x(data, pdbs, folds)
    emit_d5_signprice(data, pdbs, folds)
    print("\ndone.  Every file now carries provenance, a pinned seed where one is used, and a "
          "`reemit_check` block against the previous artefact.")


if __name__ == "__main__":
    main()
