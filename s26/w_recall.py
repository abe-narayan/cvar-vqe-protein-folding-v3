#!/usr/bin/env python
"""s26/w_recall.py -- the partial-recall gradient: sequence proximity of each dev target to its
own fold model's training corpus (native-free covariates) against the built-chain RMSD (gated).
Lane W, Sprint 26.  Pre-registration: s26/PREREG_partial_recall_gradient.md.

    python s26/jobrun.py --agent W --tag CPU --name w_recall_cov --est-ram 0.3 -- python s26/w_recall.py covariates
    python s26/jobrun.py --agent W --tag CPU --name w_recall_endpoint --est-ram 0.3 -- python s26/w_recall.py endpoint
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402
from scipy.stats import spearmanr                    # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import data as D                           # noqa: E402
import w_selfcopy as W                               # noqa: E402
import p_ladder as PL                                # noqa: E402

COVS = ("I_long", "I_short", "L_kmer")


def longest_shared_substring(a, corpus_set_by_len):
    """Longest exact substring of `a` present in any training chain (via per-length substring sets)."""
    n = len(a)
    for L in range(n, 0, -1):
        S = corpus_set_by_len.get(L)
        if S is None:
            continue
        for off in range(0, n - L + 1):
            if a[off:off + L] in S:
                return L
    return 0


def covariates(verbose=True):
    tg = I.targets()
    rows = []
    corp_cache = {}
    for t in tg:
        pdb, seq, fold = t["pdb"], t["seq"], t["fold"]
        if fold not in corp_cache:
            ents = PL.train_entries(fold)
            seqs = [e.seq for e in ents]
            C, lens = D.composition_matrix(seqs)
            sets = {}
            for s in seqs:
                for L in range(1, len(s) + 1):
                    if L > 16:
                        break
                    sets.setdefault(L, set()).update(s[o:o + L] for o in range(0, len(s) - L + 1))
            corp_cache[fold] = (seqs, C, lens, sets)
        seqs, C, lens, sets = corp_cache[fold]
        assert seq not in seqs, "the target's own chain is in its fold model's corpus"
        bound_long = np.asarray(D.max_possible_identity_many(seq, C, lens), float)
        cand = np.flatnonzero(bound_long >= 0.3)
        vals_long = D.identity_many(seq, [seqs[k] for k in cand]) if len(cand) else np.zeros(0)
        n = len(seq)
        ca = D.composition_matrix([seq])[0][0]                       # (20,) counts of the target
        shared = np.minimum(ca[None, :], C).sum(1)                   # the admissible multiset bound
        bound_short = shared / np.maximum(np.minimum(lens, n), 1)    # over the SHORTER length
        cand_s = np.flatnonzero(bound_short >= 0.5)
        vals_short = D.identity_many(seq, [seqs[k] for k in cand_s], norm="shorter") if len(cand_s) else np.zeros(0)
        row = {"pdb": pdb, "n": n, "fold": fold, "n_corpus": len(seqs),
               "I_long": float(vals_long.max()) if len(vals_long) else 0.0,
               "I_short": float(vals_short.max()) if len(vals_short) else 0.0,
               "L_kmer": int(longest_shared_substring(seq, sets)),
               "n_candidates_long": int(len(cand)), "n_candidates_short": int(len(cand_s))}
        rows.append(row)
        if verbose:
            print("  %s n=%2d f%d corpus %d  I_long %.3f  I_short %.3f  L_kmer %d" % (pdb, n, fold, len(seqs), row["I_long"], row["I_short"], row["L_kmer"]), flush=True)
    summ = {c: {"mean": float(np.mean([r[c] for r in rows])), "min": float(np.min([r[c] for r in rows])), "max": float(np.max([r[c] for r in rows])),
                "sd": float(np.std([r[c] for r in rows]))} for c in COVS}
    summ["I_long_max_over_targets"] = float(max(r["I_long"] for r in rows))
    p = W.save("recall_covariates", {"label": "recall covariates, native-free (sequences only)", "summary": summ, "rows": rows},
               rows=rows, n_expected=len(tg), complete_keys=("pdb",) + COVS)
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


def fold_boot_rho(x, y, folds, n_boot=4000, seed_parts=("w_recall",)):
    from s15 import seed as SD
    rng = SD.stable_rng(*seed_parts)
    F = np.array(sorted(set(folds.tolist())))
    out = []
    for _ in range(n_boot):
        pick = rng.choice(F, len(F), replace=True)
        idx = np.concatenate([np.flatnonzero(folds == q) for q in pick])
        out.append(spearmanr(x[idx], y[idx]).correlation)
    out = np.array(out)
    return [float(np.nanpercentile(out, 2.5)), float(np.nanpercentile(out, 97.5))]


def perm_p(x, y, rho_obs, n=500, seed_parts=("w_recall", "perm")):
    from s15 import seed as SD
    rng = SD.stable_rng(*seed_parts)
    cnt = 0
    for _ in range(n):
        r = spearmanr(x, rng.permutation(y)).correlation
        if r <= rho_obs:          # one-sided in the direction of HELP (negative rho)
            cnt += 1
    return (cnt + 1) / (n + 1)


def endpoint(verbose=True):
    W.require_signoff("recall endpoint")
    Z = W.load_result("recall_covariates")
    if not (Z and Z.get("complete")):
        raise SystemExit("run `covariates` first")
    rows = Z["rows"]; pdbs = [r["pdb"] for r in rows]; folds = ST.pinned_folds(pdbs)
    recs = {p: I.shipped_record(p) for p in pdbs}
    y = {"rmsd_arm": np.array([recs[p]["rmsd_arm"] for p in pdbs]), "rmsd_avg": np.array([recs[p]["rmsd_avg"] for p in pdbs]),
         "shipped": np.array([recs[p]["shipped"] for p in pdbs])}
    n_len = np.array([r["n"] for r in rows], float)
    out = {"label": "partial recall gradient, GATED", "prereg": "s26/PREREG_partial_recall_gradient.md", "n": len(rows), "mde_rho": 2.8016 / np.sqrt(len(rows) - 3), "stats": {}}
    for c in COVS:
        x = np.array([r[c] for r in rows], float)
        for yk, yv in y.items():
            rho = float(spearmanr(x, yv).correlation)
            ci = fold_boot_rho(x, yv, folds, seed_parts=("w_recall", c, yk))
            p = perm_p(x, yv, rho, seed_parts=("w_recall", "perm", c, yk))
            # partial out chain length: residuals of rank-regressions on n
            from scipy.stats import rankdata
            rx, ry, rn = rankdata(x), rankdata(yv), rankdata(n_len)
            def resid(a, b):
                A = np.column_stack([np.ones_like(b), b]); coef, *_ = np.linalg.lstsq(A, a, rcond=None); return a - A @ coef
            rho_partial = float(np.corrcoef(resid(rx, rn), resid(ry, rn))[0, 1])
            out["stats"]["%s_vs_%s" % (c, yk)] = {"rho": rho, "ci95_fold": ci, "perm_p_one_sided_help": p, "rho_partial_n": rho_partial,
                                                    "gradient_exists": bool(rho <= -0.25 and ci[1] < 0 and p < 0.017),
                                                    "harmful_gradient": bool(rho >= 0.25 and ci[0] > 0)}
            if verbose:
                print("  %-8s vs %-8s rho %+.3f  fold CI [%+.3f, %+.3f]  perm p(help) %.3f  partial(n) %+.3f  -> %s"
                      % (c, yk, rho, ci[0], ci[1], p, rho_partial, "GRADIENT" if out["stats"]["%s_vs_%s" % (c, yk)]["gradient_exists"] else ("HARMFUL" if out["stats"]["%s_vs_%s" % (c, yk)]["harmful_gradient"] else "none at MDE 0.25")))
    print("  MDE in rho at n=%d: %.3f" % (len(rows), out["mde_rho"]))
    p = W.save("recall_endpoint", out)
    print("  ->", p)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("covariates", "endpoint"))
    a = ap.parse_args(argv)
    covariates() if a.cmd == "covariates" else endpoint()
    return 0


if __name__ == "__main__":
    sys.exit(main())
