"""s21/d_lrank.py -- D10: CAN A NATIVE-FREE READOUT TRAVERSE THE OBJECTIVE'S TOP-M ON THE LATENT?

    python -m s21.d_lrank run | report

THE QUESTION, and it is the one left standing after the exhaustive enumeration.  `latentrank.py`
shows an ORACLE ceiling of ~1.99 A inside the objective's TOP-512 against the deployed argmin's
~3.28 A -- about 1.3 A of headroom in a 512-candidate window.  But **every top-M number there is
an ORACLE**: "picks perfectly inside the top-M" consumes the native.

    THE FINDING IS NOT THE CEILING.  IT IS WHETHER ANY NATIVE-FREE READOUT CAN TRAVERSE IT.

That is measured on the POOL in this project (S12's in-band ordering, 0.600 across targets) and is
NOT measured on the LATENT.  And `s21/d_distobj.py` (D7) showed a pool-measured relationship does
not automatically transfer, so the pool's bounds do not carry over in either direction.

THE DESIGN IS DELIBERATELY THE SAME OPERATOR DECOMPOSITION `s21/d_cvarop.py` RAN ON THE POOL, so
the two are directly comparable rather than merely analogous.  Within the objective's top-M window
of the ENUMERATED latent, four readouts, each against its MATCHED-COUNT control:

    argmin        lowest OBJECTIVE value in the window             the deployed selector
    member        a uniformly random member of the top-M          worst case for a converged law
    medoid        the top-M's medoid                              single structure, set-informed
    avg           the top-M's coordinate average                  POINT CLOUD basis, labelled

    rand_argmin / rand_member / rand_medoid / rand_avg
                  the SAME four readouts on a RANDOM M-window (8 draws).  `rand_argmin` is the
                  objective's argmin computed INSIDE the random window, so the random arm is a
                  full native-free selection battery over M arbitrary latent draws -- which is
                  the set where a sub-2.0 A structure demonstrably lives once M = 512.
    ORACLE_topM                            the ceiling -- native-consuming, labelled every time

PRE-REGISTERED IN `s21/PREREG_D.md` D10, before this module produced a number.

  T1  medoid and avg beat argmin on the latent's top-M, as they do on the pool
      (`tail_medoid|disto` -0.447, `tail_avg|disto` -0.377 against matched random).
  T2  THE ONE THAT MATTERS.  Against the MATCHED-COUNT RANDOM WINDOW the gain will be SMALLER on
      the latent than on the pool, and I am not confident it clears the MDE.  Mechanism, from R1's
      falsification: what makes averaging work is ERROR INCOHERENCE, not set spread, and the
      latent's top-M are all built from the SAME per-residue basin mixtures, so their errors should
      be far more COHERENT than independently retrieved pool windows'.  Coherent error survives
      averaging.
  T3  No native-free readout reaches the ORACLE top-M ceiling at any M <= 512.

  FALSIFIER.  If a native-free readout beats its matched-count random control by more than that
  comparison's OWN MDE with a CI excluding zero at EVERY M, the latent's top-M is traversable
  native-free, T2 is REFUTED, and that is a POSITIVE result -- the most valuable outcome available.
  If every readout ties its matched control, the top-M window carries no native-free structure
  beyond what the objective already used, and the reranker route is closed for THIS objective's
  window despite a ceiling 1.3 A below its argmin.

OPERATOR FORKS, DECLARED (BRIEF section 7 rule 0).  Functional: Bayes, matching `pool_argmin_rb`
and the `latentrank` primary; squared NOT taken, and D7 licenses that as immaterial (median
rho 0.973).  Basis: BUILT CHAIN from `build_ca_exact` on both sides -- no window-vs-rebuild
mismatch anywhere; the `avg` arm is a POINT CLOUD and is labelled on its row.  Readout: all four
reported, none privileged.  Null: matched-count random window, same M, same readout, 8 draws.
Normalisation: none, both arms are RMSD in Angstroms.  Ties: stable sort throughout.

DECLARED LIMITATION, IN ADVANCE.  `n <= 13` only (75 of 126) -- the same length-defined subset D8
was run to remove.  Accepted here because this needs the top-M STRUCTURES rather than a summary and
`n >= 14` costs 8x more.  The result is labelled an `n <= 13` result and is NOT generalised.
"""
from __future__ import annotations

import json
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

from s12 import instrument as I      # noqa: E402
from s15 import seed as SD           # noqa: E402
from s18 import phys_lib as PL       # noqa: E402
from s19 import qb_lib as QB         # noqa: E402

MAX_N = 13                  # 2**13 = 8192; see the declared limitation above
MS = (1, 8, 64, 512)
R_DRAW = 8                  # draws for every stochastic readout and every random-window control
CHUNK = 2048
CFG = {"MAX_N": MAX_N, "MS": list(MS), "R_DRAW": R_DRAW,
       "functional": "shipped Bayes-risk (I.shipped_score) on the enumerated modes",
       "basis": "BUILT CHAIN both sides; the avg arm is a POINT CLOUD"}


def _readouts(sub_CA, sub_true, nat, rng, tag, out, score=None):
    """The four readouts on one M-window of built chains.  `sub_true` is ORACLE, scoring only.

    `score` is the window's NATIVE-FREE objective value per member.  For the top-M window it is
    already sorted, so index 0 is the argmin; for a RANDOM window it must be supplied, or
    "argmin" would silently mean "the first member in index order" -- which is a random member,
    not a selection.  That distinction is the whole point of the random-window arm.
    """
    m = len(sub_CA)
    k0 = 0 if score is None else int(np.argmin(np.asarray(score, float)))
    out[f"{tag}_argmin"] = float(sub_true[k0])
    out[f"{tag}_member"] = float(np.mean([sub_true[rng.integers(0, m)] for _ in range(R_DRAW)]))
    if m == 1:
        out[f"{tag}_medoid"] = float(sub_true[0])
        out[f"{tag}_avg"] = float(sub_true[0])
        return
    P = I.pairwise_rmsd(sub_CA)
    out[f"{tag}_medoid"] = float(sub_true[int(I.medoid(P))])
    a, _b = I.coordinate_average(sub_CA, P)
    out[f"{tag}_avg"] = float(I.ca_rmsd(np.asarray(a, float), nat))
    #: NATIVE-FREE geometric spread of the window -- free, P is already built
    out[f"{tag}_spread"] = float(P[np.triu_indices(m, 1)].mean())


def run(targets=None):
    tg = [t for t in (targets if targets is not None else I.targets()) if int(t["n"]) <= MAX_N]
    print(f"targets with n <= {MAX_N}: {len(tg)}", flush=True)
    rows, t0 = [], time.time()
    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        tgt = QB.target(pdb)
        nat = tgt["nat"]
        mu = tgt["mu"]
        N = 1 << n
        ob = QB.new_obj(tgt)
        #: enumerate every latent mode, keeping only the running best-512 by the BAYES score
        keep = max(MS)
        bestE = np.full(keep, np.inf)
        bestCA = np.zeros((keep, n, 3))
        allR_sum, allR_cnt = 0.0, 0
        oracle_all = np.inf
        #: the FULL score and RMSD vectors, 2**n <= 8192 floats each -- 64 KB, so keeping them is
        #: free and makes every downstream statistic EXACT and post-hoc rather than incremental.
        eb_all = np.empty(N); rr_all = np.empty(N)
        for s0 in range(0, N, CHUNK):
            k = np.arange(s0, min(s0 + CHUNK, N), dtype=np.int64)
            bits = ((k[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
            phi = mu[np.arange(n)[None, :], bits, 0]
            psi = mu[np.arange(n)[None, :], bits, 1]
            _e, CA = ob.raw(phi, psi)
            eb = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(CA, tgt["i"], tgt["j"])), float)
            rr = I.kabsch_rmsd_batch(CA, nat)                   # ORACLE, scoring only
            oracle_all = min(oracle_all, float(rr.min()))
            allR_sum += float(rr.sum()); allR_cnt += len(k)
            eb_all[k[0]:k[-1] + 1] = eb; rr_all[k[0]:k[-1] + 1] = rr
            #: merge this chunk into the running best-`keep`
            E = np.concatenate([bestE, eb])
            C = np.concatenate([bestCA, CA], 0)
            o = np.argsort(E, kind="stable")[:keep]
            bestE, bestCA = E[o], C[o]
        top_true = I.kabsch_rmsd_batch(bestCA, nat)             # ORACLE, scoring only
        rng = SD.stable_rng(pdb, "s21_d_lrank")

        e = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "N_latent": N,
             "latent_mean": allR_sum / max(allR_cnt, 1), "latent_oracle": oracle_all}

        #: ---------------------------------------------------------------- NATIVE-FREE per-target
        #: diagnostics from the OBJECTIVE'S OWN SCORE DISTRIBUTION -- nothing here touches `rr`.
        #: The coordinator's L19 shows the objective's containment is BIMODAL: it contains the
        #: ORACLE-best in its top-512 on 59% of targets and is WORSE THAN RANDOM on the other 41%.
        #: If a native-free quantity separates those two regimes the route is worth ~0.4 A; if not,
        #: the aggregate null is what ships.  These are the candidates, and they are free here.
        eo = np.sort(eb_all)
        med, sd = float(np.median(eb_all)), float(eb_all.std())
        e["nf_score_z_top1"] = float((eo[0] - eb_all.mean()) / max(sd, 1e-12))
        e["nf_score_margin"] = float((eo[1] - eo[0]) / max(med - eo[0], 1e-12))
        e["nf_score_conc512"] = float((eo[min(511, N - 1)] - eo[0]) / max(med - eo[0], 1e-12))
        e["nf_score_iqr_over_range"] = float(
            (np.percentile(eb_all, 75) - np.percentile(eb_all, 25))
            / max(eo[-1] - eo[0], 1e-12))
        e["nf_score_skew"] = float(((eb_all - eb_all.mean()) ** 3).mean() / max(sd ** 3, 1e-12))
        #: ORACLE, EVALUATION ONLY -- the label the native-free columns above are asked to predict
        k_or = int(np.argmin(rr_all))
        e["ORACLE_rank_pct"] = float((eb_all < eb_all[k_or]).sum() / N * 100.0)
        for M in MS:
            e[f"ORACLE_contained{M}"] = float(e["ORACLE_rank_pct"] <= 100.0 * M / N)
        for M in MS:
            _readouts(bestCA[:M], top_true[:M], nat, rng, f"top{M}", e)
            e[f"ORACLE_top{M}"] = float(top_true[:M].min())
            #: MATCHED-COUNT control: the SAME readouts on a RANDOM M-window of the same latent.
            #: Re-enumerating for a random window is unaffordable, so the random window is drawn
            #: from a uniformly sampled set of the same size -- declared, and it is the same
            #: object the `latent_mean` column already prices at M = 1.
            acc = {}
            for _ in range(R_DRAW):
                idx = rng.choice(N, M, replace=False)
                bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
                phi = mu[np.arange(n)[None, :], bits, 0]
                psi = mu[np.arange(n)[None, :], bits, 1]
                _e2, CA2 = ob.raw(phi, psi)
                rr2 = I.kabsch_rmsd_batch(CA2, nat)
                #: the window's own NATIVE-FREE objective values, so `r_argmin` is a genuine
                #: score-based selection INSIDE the random window rather than an arbitrary member
                sc2 = np.asarray(I.shipped_score(tgt["dg"],
                                                 I.pair_dists(CA2, tgt["i"], tgt["j"])), float)
                tmp = {}
                _readouts(CA2, rr2, nat, rng, "r", tmp, score=sc2)
                for kk, vv in tmp.items():
                    acc.setdefault(kk, []).append(vv)
                acc.setdefault("r_ORACLE", []).append(float(rr2.min()))
            for kk, vv in acc.items():
                e[f"rand{M}_{kk[2:] if kk.startswith('r_') else kk}"] = float(np.mean(vv))
        rows.append(e)
        if (c + 1) % 5 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            _save(rows, len(tg))
    _save(rows, len(tg))
    report(rows)
    return rows


def _save(rows, n_expected):
    obj = {"rows": rows, "cfg": CFG, "n_expected": int(n_expected)}
    p = os.path.join(RESULTS, "d_lrank.json")
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, p)                       # atomic, per the defect found in latentfull.py
    need = ["latent_mean", "latent_oracle", "ORACLE_rank_pct", "nf_score_z_top1",
            "nf_score_margin", "nf_score_conc512", "nf_score_iqr_over_range", "nf_score_skew"] + [
        f"{a}{M}_{r}" for M in MS for a in ("top", "rand")
        for r in ("argmin", "member", "medoid", "avg")] + [f"ORACLE_top{M}" for M in MS]
    ok = (len(rows) == n_expected == 75
          and all(np.isfinite(r.get(k, np.nan)) for r in rows for k in need))
    fp = os.path.join(RESULTS, "d_lrank.COMPLETE")
    if ok:
        with open(fp, "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                     f"n=75 (every target with n<={MAX_N}) M={list(MS)} R_DRAW={R_DRAW} "
                     f"readouts=argmin,member,medoid,avg controls=matched-count random window "
                     f"keys={len(need)} functional=Bayes basis=built-chain\n")
    elif os.path.exists(fp):
        os.remove(fp)


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RESULTS, "d_lrank.json")))["rows"]
    g = lambda k: np.array([r[k] for r in rows], float)     # noqa: E731
    folds = g("fold")
    n = len(rows)
    print(f"\n=== D10  NATIVE-FREE READOUTS INSIDE THE OBJECTIVE'S TOP-M ON THE LATENT, n = {n} ===")
    print(f"Every target with n <= {MAX_N} (a length-defined subset, declared -- NOT generalised).")
    print("Basis: BUILT CHAIN both sides; the `avg` row is a POINT CLOUD.  ORACLE rows labelled.\n")
    print(f"  latent mean {g('latent_mean').mean():.3f}   "
          f"latent ORACLE (whole space) {g('latent_oracle').mean():.3f}\n")
    print(f"  {'M':>5}  {'readout':<9}{'top-M':>9}{'randM':>9}   "
          f"{'top - rand (the ONLY comparison that means anything)':<52}")
    out = {"n": n, "cfg": CFG, "cells": {}}
    for M in MS:
        print(f"  ORACLE ceiling inside top-{M}: {g(f'ORACLE_top{M}').mean():.3f}   "
              f"(random {M}-window ORACLE {g(f'rand{M}_ORACLE').mean():.3f})")
        for ro in ("argmin", "member", "medoid", "avg"):
            a, b = g(f"top{M}_{ro}"), g(f"rand{M}_{ro}")
            st = PL.paired(a, b, folds=folds)
            ci = st.get("ci_fold", st["ci"])
            mde = 2.8016 * (a - b).std(ddof=1) / np.sqrt(n)
            flag = ("  BEATS its control" if (st["mean"] < 0 and ci[1] < 0)
                    else ("  WORSE" if (st["mean"] > 0 and ci[0] > 0) else "  NOT MEASURED"))
            out["cells"][f"{M}|{ro}"] = {"top": float(a.mean()), "rand": float(b.mean()),
                                         "diff": st["mean"], "ci95_fold": list(ci),
                                         "ci95_iid": list(st["ci"]), "W": st["W"], "L": st["L"],
                                         "own_mde": float(mde)}
            print(f"  {M:>5}  {ro:<9}{a.mean():>9.3f}{b.mean():>9.3f}   "
                  f"{st['mean']:+7.4f} fold[{ci[0]:+.4f},{ci[1]:+.4f}] "
                  f"W/L {st['W']}/{st['L']} MDE {mde:.4f}{flag}")
        print()
    print("  T1 -- do medoid/avg beat the deployed argmin at the same M?")
    for M in MS[1:]:
        for ro in ("medoid", "avg"):
            st = PL.paired(g(f"top{M}_{ro}"), g("top1_argmin"), folds=folds)
            ci = st.get("ci_fold", st["ci"])
            print(f"    M={M:<4}{ro:<8}vs the deployed argmin: {st['mean']:+7.4f} "
                  f"fold[{ci[0]:+.4f},{ci[1]:+.4f}] W/L {st['W']}/{st['L']}")
    print("\n  T3 -- does any native-free readout reach the ORACLE top-M ceiling?")
    for M in MS[1:]:
        best = min(("member", "medoid", "avg"), key=lambda r: g(f"top{M}_{r}").mean())
        st = PL.paired(g(f"top{M}_{best}"), g(f"ORACLE_top{M}"), folds=folds)
        ci = st.get("ci_fold", st["ci"])
        print(f"    M={M:<4}best native-free ({best}) minus ORACLE_top{M}: {st['mean']:+7.4f} "
              f"fold[{ci[0]:+.4f},{ci[1]:+.4f}] W/L {st['W']}/{st['L']}")
    # ---------------------------------------------------------------- THE CONTAINMENT COLUMN
    #: The coordinator's L19: the objective's containment of the ORACLE-best is BIMODAL -- it
    #: contains it in the top-512 on ~59% of targets and is WORSE THAN RANDOM on the other ~41%.
    #: If a NATIVE-FREE quantity separates those regimes the route is worth ~0.4 A; if not, the
    #: aggregate null is what ships.  `ORACLE_rank_pct` / `ORACLE_contained*` consume the native
    #: and are the LABEL; every `nf_*` and `*_spread` column is native-free and is the PREDICTOR.
    print("\n  THE CONTAINMENT COLUMN -- can anything NATIVE-FREE tell the two regimes apart?")
    lab = g("ORACLE_contained512")
    print(f"    ORACLE-best contained in the objective's top-512 on {int(lab.sum())}/{n} = "
          f"{100*lab.mean():.0f}% of targets "
          f"(a random window gives {100*np.mean(512.0/g('N_latent')):.0f}%)")
    print(f"    ORACLE_rank_pct: median {np.median(g('ORACLE_rank_pct')):.2f}  "
          f"mean {g('ORACLE_rank_pct').mean():.2f}   (50.0 is the no-skill null)")
    preds = [k for k in rows[0] if k.startswith("nf_")] + \
            [f"{w}{M}_spread" for M in MS if M > 1 for w in ("top", "rand")]
    preds = [k for k in preds if k in rows[0]]
    out["containment"] = {}
    print(f"\n    {'native-free predictor':<26}{'rho vs rank_pct':>17}{'|AUC| vs contained':>20}")
    ry = np.empty(n); ry[np.argsort(g("ORACLE_rank_pct"), kind="stable")] = np.arange(n)

    def _auc(v, l):
        pos, neg = v[l > 0.5], v[l < 0.5]
        if not len(pos) or not len(neg):
            return float("nan")
        a = float(np.mean([(x < y) + 0.5 * (x == y) for x in pos for y in neg]))
        return max(a, 1.0 - a)                    # folded: direction is free, so must be priced

    for k in preds:
        v = g(k)
        rk = np.empty(n); rk[np.argsort(v, kind="stable")] = np.arange(n)
        rho = float(np.corrcoef(rk, ry)[0, 1])
        a = _auc(v, lab)
        out["containment"][k] = {"rho_vs_rank_pct": rho, "folded_auc_vs_contained": a}
        print(f"    {k:<26}{rho:>17.3f}{a:>20.3f}")

    #: THE BEST-OF-K NULL (PREREG D11, registered BEFORE this was read).  A maximum over K noisy
    #: AUCs is biased upward, so the bar for "the best predictor" is the 95th percentile of the
    #: PERMUTED maximum -- not 0.5, and not a threshold picked by assertion.  Simulated at this
    #: experiment's own n and split, pure-noise labels give a best-of-10 |AUC| whose MEDIAN is
    #: 0.625: the whole 0.55-0.65 band is inside the noise.
    rngp = np.random.default_rng(20210907)
    V = np.array([g(k) for k in preds])
    mx = np.empty(4000)
    for b in range(len(mx)):
        lp = lab[rngp.permutation(n)]
        mx[b] = max(_auc(row, lp) for row in V)
    obs = max(v["folded_auc_vs_contained"] for v in out["containment"].values()
              if np.isfinite(v["folded_auc_vs_contained"]))
    p95 = float(np.percentile(mx, 95))
    out["bestofK_null"] = {"K": len(preds), "observed_best_folded_auc": float(obs),
                           "perm_median": float(np.median(mx)), "perm_p95": p95,
                           "DEMONSTRATED": bool(obs > p95)}
    print(f"\n    BEST-OF-K NULL (PREREG D11): K = {len(preds)} predictors, label permuted 4000x")
    print(f"      permuted best-of-K |AUC|:  median {np.median(mx):.3f}   95th pct {p95:.3f}")
    print(f"      OBSERVED best |AUC|:       {obs:.3f}")
    print(f"      -> the bimodal route is "
          f"{'DEMONSTRATED' if obs > p95 else 'NOT DEMONSTRATED'} by this instrument.")
    print("    Registered rule (D11): at or below the permuted 95th percentile is NOT MEASURED and")
    print("    NOT DEMONSTRATED -- never 'promising' or 'suggestive'.  0.55-0.65 IS the noise band.")

    print("\nREAD.  The ORACLE rows consume the native and are ceilings, never achievements.")
    print("The `top - rand` column is the whole experiment: it is the objective's top-M window's")
    print("native-free content, in the currency a reranker actually spends.")
    json.dump(out, open(os.path.join(RESULTS, "d_lrank_report.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
