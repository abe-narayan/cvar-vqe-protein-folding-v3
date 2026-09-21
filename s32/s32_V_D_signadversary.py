"""S32 LANE V -- adversarial replication of lane D's D1-T "the per-target in-band SIGN is a
real latent and it TRANSFERS across a split half".

Lane D reports (s32/results/s32_D1_signtransfer.json, 16 splits/target, label-permutation null):

    scorer       transfer   null     excess   se      xMDE  folds
    AMBER         +0.1125  -0.0056  +0.1181  0.0164   2.57   5/5
    DIS           +0.1890  -0.0033  +0.1923  0.0226   3.04   5/5
    LEG_total     +0.2263  -0.0003  +0.2266  0.0237   3.41   5/5
    LEG_torsion   +0.1499  -0.0008  +0.1507  0.0196   2.75   5/5

**The label-permutation null cannot test the claim that is being made.**  Permuting rr inside B
destroys *all* structure on B, so it asks "is there any relation at all".  The claim is stronger
and more specific: that the relation's SIGN is a property of THE TARGET.  Three things are
therefore added here, none of which is in lane D's artefact:

 1. **The cross-target sign-shuffle null.**  Apply target t' 's sign to target t's held-out half.
    `E[sign_A(t') * rho_B(t)] = E[sign] * E[rho]`, so this null is exactly the value a single
    GLOBAL sign would deliver.  Any excess over it is the per-target part -- which is the claim.

 2. **The global-sign baseline, leave-fold-out.**  What one sign for the whole instrument buys.
    A scorer whose sign is globally consistent (typicality is) will show large "transfer" that is
    not per-target at all, and the label-permutation null will not notice.

 3. **A NATIVE-FREE, GLOBALLY-SIGNED CONTROL SCORER matched to the operator's own space**
    (contract rule 6): band typicality, and Rg.  If a native-free scorer with a known global sign
    transfers as well as AMBER, then "the sign is a latent needing an oracle" is not what the
    number shows.  Plus a **pure-noise scorer**, which must transfer at ~0 -- the test that makes
    this audit able to fail (contract rule 5).

 4. **A duplicate audit of the band.**  The split-half argument needs A and B to be independent
    samples.  The 500-pool is known to contain duplicates (`n_distinct` 66 of 500 on 1A13); an
    exact duplicate landing in both halves couples them.

Only DIS is recomputed here -- AMBER and Legacy are lane D's to re-run -- so this is an
independent replication of ONE of the four rows plus the controls all four need.

ORACLE / NOT DEPLOYABLE: `rr` is read to form every rho.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s24 import stats_lib as ST                                          # noqa: E402
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)
N_SPLITS = 16                 # lane D's protocol, reproduced exactly
SEED = 32_0032_11


def spear(x, y):
    from scipy.stats import rankdata
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 3 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return np.nan
    a, b = rankdata(x), rankdata(y)
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else np.nan


def per_target(pdb, seq, fold, n, rng):
    u = I.load_univ(pdb)
    pool = I.pool_idx(u, I.K)
    W = u["W"][pool]
    rr = u["rr"][pool]
    sub = np.asarray(I.shipped_record(pdb)["sub"], int)
    Wb, rrb = W[sub], rr[sub]
    m = len(sub)

    # ---- scorers on the band
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    D = I.pair_dists(W, i, j)
    dis = I.shipped_score(dg, D.astype(np.float32).astype(float))[sub]
    P = I.pairwise_rmsd(Wb)
    typ = P.mean(1)                                            # native-free: band typicality
    rg = np.sqrt(((Wb - Wb.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))     # native-free
    noise = rng.standard_normal(m)                             # no per-target latent, by design
    scorers = {"DIS": dis, "TYPICALITY": typ, "RG": rg, "NOISE": noise}

    # ---- duplicate audit: exact coordinate duplicates inside the band
    dup = int(m - (P < 1e-6).any(0).sum()) if m else 0
    n_distinct = int(len({tuple(np.round(x.ravel(), 6)) for x in Wb}))

    out = dict(pdb=pdb, fold=fold, m=m, n_distinct_band=n_distinct,
               dup_frac=float(1.0 - n_distinct / max(m, 1)))

    #: DEDUPED band: one representative per exact-coordinate duplicate class.  The split-half
    #: argument needs A and B to be independent samples, and a duplicate landing in both halves
    #: contributes the same (score, rr) pair to each, guaranteeing sign agreement for free.
    seen, keep = set(), []
    for q in range(m):
        key = tuple(np.round(Wb[q].ravel(), 6))
        if key not in seen:
            seen.add(key); keep.append(q)
    keep = np.asarray(keep, int)
    for nm, v in list(scorers.items()):
        scorers[nm + "_DEDUP"] = v[keep]
    rr_of = {nm: (rrb[keep] if nm.endswith("_DEDUP") else rrb) for nm in scorers}
    for nm, v in scorers.items():
        rrb_use = rr_of[nm]
        out["rho_full_%s" % nm] = spear(v, rrb_use)
        tr, nl, half = [], [], []
        mm = len(v)
        for _ in range(N_SPLITS):
            perm = rng.permutation(mm)
            A, B = perm[: mm // 2], perm[mm // 2:]
            rA = spear(v[A], rrb_use[A]); rB = spear(v[B], rrb_use[B])
            if not (np.isfinite(rA) and np.isfinite(rB)):
                continue
            s = 1.0 if rA >= 0 else -1.0
            tr.append(s * rB)
            half.append(rB)
            nl.append(s * spear(v[B], rng.permutation(rrb_use[B])))   # lane D's matched null
        out["transfer_%s" % nm] = float(np.mean(tr)) if tr else np.nan
        out["null_perm_%s" % nm] = float(np.mean(nl)) if nl else np.nan
        out["rhoB_%s" % nm] = float(np.mean(half)) if half else np.nan     # UNSIGNED half rho
        out["absrhoB_%s" % nm] = float(np.mean(np.abs(half))) if half else np.nan
        out["signA_%s" % nm] = float(np.sign(out["rho_full_%s" % nm] or 0.0))
    return out


def main():
    rng = np.random.default_rng(SEED)
    rows = [per_target(t["pdb"], t["seq"], t["fold"], t["n"], rng) for t in I.targets()]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows], int)
    out = dict(n=len(rows), n_splits=N_SPLITS, seed=SEED,
               basis="in-band Spearman vs ORACLE rr inside the shipped top-75; DIAGNOSTIC, "
                     "never differenced against a chain RMSD",
               oracle="ORACLE / NOT DEPLOYABLE", provenance=ST.provenance(__file__), scorers={})

    dupf = np.array([r["dup_frac"] for r in rows])
    out["band_duplicates"] = dict(mean_frac=float(dupf.mean()), max_frac=float(dupf.max()),
                                  n_targets_with_any=int((dupf > 0).sum()),
                                  mean_n_distinct=float(np.mean([r["n_distinct_band"] for r in rows])))
    print("=" * 106)
    print("BAND DUPLICATES: mean %.1f%% of the 75 band members are exact coordinate duplicates; "
          "%d/126 targets have any; mean distinct %.1f"
          % (100 * dupf.mean(), (dupf > 0).sum(), out["band_duplicates"]["mean_n_distinct"]))
    print("  -> a duplicate landing in both halves couples them; the split-half argument assumes "
          "they are independent samples.")
    print("=" * 106)
    print("%-18s %9s %9s %9s %9s %9s %9s %6s  %s"
          % ("scorer", "transfer", "nullPERM", "nullXTGT", "globalSGN", "mean rho",
             "mean|rhoB|", "folds", "verdict on the PER-TARGET claim"))

    for nm in ("DIS", "TYPICALITY", "RG", "NOISE",
               "DIS_DEDUP", "TYPICALITY_DEDUP", "RG_DEDUP", "NOISE_DEDUP"):
        tr = np.array([r["transfer_%s" % nm] for r in rows], float)
        nl = np.array([r["null_perm_%s" % nm] for r in rows], float)
        rhoB = np.array([r["rhoB_%s" % nm] for r in rows], float)
        absB = np.array([r["absrhoB_%s" % nm] for r in rows], float)
        sA = np.array([r["signA_%s" % nm] for r in rows], float)
        rho_full = np.array([r["rho_full_%s" % nm] for r in rows], float)
        ok = np.isfinite(tr) & np.isfinite(rhoB)

        #: NULL 2 -- the cross-target sign shuffle.  200 shuffles, its own distribution.
        xt = []
        for _ in range(200):
            q = rng.permutation(ok.sum())
            xt.append(float((sA[ok][q] * rhoB[ok]).mean()))
        xt = np.array(xt)

        #: the GLOBAL-sign baseline, leave-fold-out so the sign cannot read the held-out target
        gl = np.empty(ok.sum()); rb = rhoB[ok]; rf = rho_full[ok]; fo = folds[ok]
        for f in np.unique(fo):
            s = np.sign(np.nanmean(rf[fo != f])) or 1.0
            gl[fo == f] = s * rb[fo == f]

        c = ST.compare(-tr[ok], -np.full(ok.sum(), float(xt.mean())), fo, names=list(np.array(pdbs)[ok]),
                       label="%s: per-target sign vs the cross-target (global) sign" % nm)
        verdict = ("PER-TARGET" if abs(c["effect_over_mde"]) >= 1.0 and c["folds_same_sign"] >= 4
                   and min(c["ci95_fold"]) * max(c["ci95_fold"]) > 0
                   else "NOT MEASURED" if abs(c["effect_over_mde"]) >= 0.7 else "NOT A RESULT")
        out["scorers"][nm] = dict(
            transfer=float(np.nanmean(tr)), null_perm=float(np.nanmean(nl)),
            null_xtarget_mean=float(xt.mean()), null_xtarget_sd=float(xt.std(ddof=1)),
            global_sign_lfo=float(gl.mean()), mean_rho_full=float(np.nanmean(rho_full)),
            mean_abs_rhoB=float(np.nanmean(absB)),
            excess_over_xtarget=float(np.nanmean(tr) - xt.mean()),
            excess_over_global_sign=float(np.nanmean(tr) - gl.mean()),
            compare_vs_xtarget=c, verdict_per_target=verdict, n=int(ok.sum()))
        print("%-18s %+9.4f %+9.4f %+9.4f %+9.4f %+9.4f %9.4f %4d/5  %s  (%.2fx MDE)"
              % (nm, np.nanmean(tr), np.nanmean(nl), xt.mean(), gl.mean(),
                 np.nanmean(rho_full), np.nanmean(absB), c["folds_same_sign"], verdict,
                 c["effect_over_mde"]))

    print("-" * 106)
    d = out["scorers"]
    print("READ: `nullPERM` is lane D's null -- it destroys ALL structure on B, so it is ~0 for "
          "every scorer INCLUDING the pure-noise one, and cannot distinguish a per-target sign "
          "from a global one.")
    print("      `nullXTGT` applies another target's sign and is therefore exactly what a single "
          "GLOBAL sign delivers.  The per-target claim lives in `transfer - nullXTGT`.")
    print("      NOISE is the falsifier: it has no per-target latent by construction and must "
          "land at ~0 on every column.  It reports transfer %+.4f."
          % d["NOISE"]["transfer"])
    print("      TYPICALITY is native-free and globally signed: its `transfer` %+.4f against a "
          "global-sign baseline of %+.4f prices how much of a large transfer needs no oracle at all."
          % (d["TYPICALITY"]["transfer"], d["TYPICALITY"]["global_sign_lfo"]))
    print("=" * 106)
    with open(os.path.join(RESULTS, "s32_V_D_signadversary.json"), "w") as fh:
        json.dump(dict(out, per_target=rows), fh, indent=1,
                  default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out


if __name__ == "__main__":
    main()
