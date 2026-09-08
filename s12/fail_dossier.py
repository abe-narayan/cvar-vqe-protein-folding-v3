"""FAIL18 forensics, step 2: cheap per-target dossier for ALL 126 tuning targets.

Everything here is computed from the cached universes + the shipped distogram + the
shipped production record.  ORACLE quantities (rr, nat_ca) are used for EVALUATION and
DIAGNOSIS only and are named with an `o_` prefix or documented as such.

Writes s12/results/fail_dossier.json
"""
from __future__ import annotations
import os, sys, json, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I


def rg(x):
    x = np.asarray(x, float)
    return float(np.sqrt(((x - x.mean(0)) ** 2).sum(1).mean()))


def ss_frac(ss):
    n = max(len(ss), 1)
    return dict(H=ss.count("H") / n, E=ss.count("E") / n, C=ss.count("C") / n)


def shell(sep):
    if sep <= 4:
        return "s2_4"
    if sep <= 8:
        return "s5_8"
    return "s9p"


def dossier(t, pep_by_pdb):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    u = I.load_univ(pdb)
    p = I.pool_idx(u)
    rr = u["rr"]                                   # ORACLE
    nat = u["nat_ca"]                              # ORACLE
    rec = I.shipped_record(pdb)
    sub = np.asarray(rec["sub"], int)              # top-75 indices INTO POOL
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    Dp = I.pair_dists(u["W"][p], i, j)
    sc = I.shipped_score(dg, Dp.astype(np.float32).astype(float))
    rrp = rr[p]

    # --- band definitions (ORACLE) -------------------------------------------------
    pool_best = float(rrp.min())
    band = np.where(rrp <= pool_best + I.BAND)[0]          # indices into pool
    band_set = set(band.tolist())
    n_band_in75 = int(len(band_set & set(sub.tolist())))
    # universe-level band, for the query-failure question
    uni_best = float(rr.min())
    uband = np.where(rr <= uni_best + I.BAND)[0]
    uband_in_pool = int(np.isin(uband, p).sum())

    # --- score geometry -----------------------------------------------------------
    order_sc = np.argsort(sc, kind="stable")
    rank_of = np.empty(len(sc), int); rank_of[order_sc] = np.arange(len(sc))
    band_pct = float(np.mean([rank_of[b] for b in band]) / len(sc)) if len(band) else float("nan")
    band_best_rank = int(min(rank_of[b] for b in band)) if len(band) else -1
    # native's own score
    Dn = I.pair_dists(nat[None], i, j)
    sc_nat = float(I.shipped_score(dg, Dn)[0])            # ORACLE diagnostic
    nat_pct = float((sc < sc_nat).mean())
    # in-pool rank correlation of score vs true rmsd
    from scipy.stats import spearmanr
    rho_all = float(spearmanr(sc, rrp).correlation)
    rho_band = float(spearmanr(sc[band], rrp[band]).correlation) if len(band) > 3 else float("nan")

    # --- pool composition ---------------------------------------------------------
    org = np.asarray(u["org"], bool); sim = np.asarray(u["sim"], float)
    seqs_pool = ["".join(I.ALPHABET[c] for c in row) for row in np.asarray(u["S"])[p]]
    comp = dict(
        nw_universe=int(len(rr)), frac_org_universe=float(org.mean()),
        frac_org_pool=float(org[p].mean()),
        frac_org_band=float(org[p][band].mean()) if len(band) else float("nan"),
        frac_org_top75=float(org[p][sub].mean()),
        n_distinct_seq_pool=int(len(set(seqs_pool))),
        sim_pool_mean=float(sim[p].mean()), sim_pool_min=float(sim[p].min()),
        sim_pool_max=float(sim[p].max()),
        sim_band_mean=float(sim[p][band].mean()) if len(band) else float("nan"),
        sim_rank_band_mean=float(np.mean(band)) if len(band) else float("nan"),
    )

    # --- native geometry (ORACLE labels) ------------------------------------------
    pep = pep_by_pdb.get(pdb)
    nat_ss = I.ss_of(pep.phi, pep.psi) if pep is not None else ""
    nat_rg = rg(nat)
    ee = float(np.linalg.norm(nat[-1] - nat[0]))
    contour = float(np.linalg.norm(np.diff(nat, axis=0), axis=1).sum())
    rg_ideal = 2.2 * n ** 0.38   # loose reference; only used for a relative axis

    # pool geometry
    rgp = np.array([rg(w) for w in u["W"][p]])
    rg_pool_mean = float(rgp.mean())
    rg_band_mean = float(rgp[band].mean()) if len(band) else float("nan")
    rg_sel = float(rgp[int(np.argmin(sc))])

    # --- distogram error structure on the NATIVE (ORACLE diagnostic) --------------
    exp = np.asarray(dg["expected"], float)
    dtrue = Dn[0]
    err = exp - dtrue
    sep = (j - i)
    shells = {}
    for s in ("s2_4", "s5_8", "s9p"):
        m = np.array([shell(x) == s for x in sep])
        if m.sum() == 0:
            continue
        shells[s] = dict(n=int(m.sum()), bias=float(err[m].mean()),
                         mae=float(np.abs(err[m]).mean()),
                         r=float(np.corrcoef(exp[m], dtrue[m])[0, 1]) if m.sum() > 2 else float("nan"))
    dgstat = dict(mae=float(np.abs(err).mean()), bias=float(err.mean()),
                  r=float(np.corrcoef(exp, dtrue)[0, 1]),
                  slope=float(np.polyfit(dtrue, exp, 1)[0]),
                  shells=shells,
                  # what the objective would predict for compactness
                  exp_mean=float(exp.mean()), true_mean=float(dtrue.mean()))

    # --- emitted structure --------------------------------------------------------
    out = dict(
        pdb=pdb, n=n, fold=fold, seq=seq, fail18=pdb in I.FAIL18,
        o_pool_best=pool_best, o_uni_best=uni_best,
        o_top75_best=float(rrp[sub].min()), o_pool_mean=float(rrp.mean()),
        o_band_size=int(len(band)), o_band_in_top75=n_band_in75,
        o_uband_size=int(len(uband)), o_uband_in_pool=uband_in_pool,
        o_band_score_pct=band_pct, o_band_best_score_rank=band_best_rank,
        o_native_score=sc_nat, o_native_score_pct=nat_pct,
        o_rho_score_rmsd=rho_all, o_rho_score_rmsd_band=rho_band,
        o_argmin_rmsd=float(rrp[int(np.argmin(sc))]),
        rmsd_fit=float(I.ca_rmsd(np.asarray(rec["fit_ca"]), nat)),
        rmsd_emitted=float(I.ca_rmsd(np.asarray(rec.get("amber_ca", rec["ca"])), nat)),
        rmsd_ca=float(I.ca_rmsd(np.asarray(rec["ca"]), nat)),
        nat_ss=nat_ss, nat_ss_frac=ss_frac(nat_ss),
        nat_rg=nat_rg, nat_rg_over_ref=nat_rg / rg_ideal, end_to_end=ee,
        ee_over_contour=ee / contour, contour=contour,
        rg_pool_mean=rg_pool_mean, rg_band_mean=rg_band_mean, rg_selected=rg_sel,
        rg_gap_nat_pool=nat_rg - rg_pool_mean,
        composition=comp, distogram=dgstat,
    )
    return out


def main():
    import peptide_db
    peps = {p.pdb: p for p in peptide_db.load()}
    tg = I.targets()
    recs = []
    for k, t in enumerate(tg):
        r = dossier(t, peps)
        recs.append(r)
        if (k + 1) % 20 == 0:
            print(f"  {k+1}/126  free={I.free_gb():.1f} GB", flush=True)
    I.write("fail_dossier", recs)
    print("wrote fail_dossier.json")

    f = [r for r in recs if r["fail18"]]; o = [r for r in recs if not r["fail18"]]
    keys = ["n", "o_pool_best", "o_uni_best", "o_top75_best", "o_pool_mean",
            "o_band_size", "o_uband_in_pool", "o_band_score_pct",
            "o_native_score_pct", "o_rho_score_rmsd", "o_rho_score_rmsd_band",
            "o_argmin_rmsd", "rmsd_fit", "rmsd_emitted", "nat_rg", "nat_rg_over_ref",
            "ee_over_contour", "rg_pool_mean", "rg_band_mean", "rg_selected",
            "rg_gap_nat_pool"]
    print(f"{'key':24s} {'FAIL18':>9s} {'other108':>9s}")
    for kk in keys:
        a = np.nanmean([r[kk] for r in f]); b = np.nanmean([r[kk] for r in o])
        print(f"{kk:24s} {a:9.3f} {b:9.3f}")
    for kk in ("mae", "bias", "r", "slope"):
        a = np.nanmean([r["distogram"][kk] for r in f]); b = np.nanmean([r["distogram"][kk] for r in o])
        print(f"{'dg.'+kk:24s} {a:9.3f} {b:9.3f}")
    for kk in ("frac_org_pool", "frac_org_band", "frac_org_top75", "sim_pool_mean",
               "n_distinct_seq_pool", "sim_rank_band_mean"):
        a = np.nanmean([r["composition"][kk] for r in f]); b = np.nanmean([r["composition"][kk] for r in o])
        print(f"{'cmp.'+kk:24s} {a:9.3f} {b:9.3f}")
    for kk in ("H", "E", "C"):
        a = np.nanmean([r["nat_ss_frac"][kk] for r in f]); b = np.nanmean([r["nat_ss_frac"][kk] for r in o])
        print(f"{'ss.'+kk:24s} {a:9.3f} {b:9.3f}")


if __name__ == "__main__":
    main()
