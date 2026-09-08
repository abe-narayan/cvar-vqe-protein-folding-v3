"""TORSION-PREDICTOR -- convert predictor quality into emitted CA-RMSD.

Sprint 12 established three separate times that PREDICTOR ACCURACY DOES NOT PRICE EMITTED
STRUCTURE, so accuracy is reported here only as a comparability anchor against S12's 0.690
ABEGO-4 number.  The currency is:

  (1) EFFECTIVE SIGMA and EFFECTIVE COVERAGE.  `sigma_eff` = RMS circular error over the
      determined angles, in degrees -- the same parameterisation as S12's oracle noise
      (independent N(0, sigma) on every angle).  Coverage = fraction of residues kept by a
      confidence gate on the predictor's own `sigma_hat`.  Those two numbers index S12's
      measured (sigma, coverage) surface, so a predictor converts into an EXPECTED emitted
      RMSD with no new assumption.
  (2) DIRECT BUILD.  The prediction goes through the same `build_ca` path S12 used, and the
      measured RMSD is compared to the surface read-off.  The two must agree.
  (3) An ERROR-SHAPE CONTROL that explains any disagreement: native torsions corrupted by
      the predictor's OWN empirical error distribution (resampled), at full coverage.  If
      the Gaussian surface over-promises, this arm localises the cause in the error shape
      rather than in the pipeline.
  (4) An extended ORACLE sigma sweep (sigma up to 100 deg), because the S12 surface stops at
      20 deg and a real sequence-only predictor lands far outside it.  ORACLE / DIAGNOSTIC.

    python -m s13.tors_eval
"""
from __future__ import annotations
import os, sys, json, math, time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s13 import tors_common as T                       # noqa: E402
from s13.tors_train import ARMS, post_path             # noqa: E402
from s12 import instrument as I                        # noqa: E402
from s12.coord_restraints import torsion_agreement     # noqa: E402
import peptide_db as pdb                               # noqa: E402

D2R = T.D2R
SEEDS = (0, 1, 2)
#: ABEGO class of each of the 324 grid cells, so posterior mass can be aggregated into the
#: four classes S12's predictor reported.
_ABEGO_CELLS = T.abego4(T.GRID_PHI, T.GRID_PSI)


# ------------------------------------------------------------------ ORACLE reference data
def native_torsions():
    """ORACLE, post-hoc evaluation and oracle-diagnostic arms only."""
    out = {}
    for t in I.targets():
        p = pdb.by_pdb(t["pdb"])
        out[t["pdb"]] = (np.asarray(p.phi, float), np.asarray(p.psi, float))
    return out


def determined_mask(n):
    """phi[0] and psi[n-1] are placeholders and do not move a single CA atom."""
    mp = np.ones(n, bool); mp[0] = False
    ms = np.ones(n, bool); ms[-1] = False
    return mp, ms


def incumbent():
    """Per-target emitted CA-RMSD of the shipped synthesis (S12's projection-based `base`)."""
    d = json.load(open(os.path.join(ROOT, "s12", "results", "coord_restraints2.json")))
    return {r["pdb"]: float(r["base"]) for r in d["per_target"]}


SURFACE = json.load(open(os.path.join(ROOT, "s12", "results",
                                      "coord_restraints2_report.json")))["arms"]
S_SIG = (6.0, 12.0, 20.0)
S_FRAC = (1.0, 0.9, 0.75, 0.5, 0.25)


def surface_readoff(sigma, frac, arm="fill_pool"):
    """Bilinear read of S12's measured (sigma, coverage) surface, clamped at its edges."""
    s = float(np.clip(sigma, S_SIG[0], S_SIG[-1]))
    f = float(np.clip(frac, S_FRAC[-1], S_FRAC[0]))
    si = int(np.searchsorted(S_SIG, s) - 1); si = max(0, min(si, len(S_SIG) - 2))
    fr = sorted(S_FRAC)
    fi = int(np.searchsorted(fr, f) - 1); fi = max(0, min(fi, len(fr) - 2))
    s0, s1 = S_SIG[si], S_SIG[si + 1]; f0, f1 = fr[fi], fr[fi + 1]
    def cell(a, b):
        return SURFACE[f"s{a:g}_f{b:g}"][arm]["mean"]
    ws = 0.0 if s1 == s0 else (s - s0) / (s1 - s0)
    wf = 0.0 if f1 == f0 else (f - f0) / (f1 - f0)
    v = ((1 - ws) * (1 - wf) * cell(s0, f0) + (1 - ws) * wf * cell(s0, f1)
         + ws * (1 - wf) * cell(s1, f0) + ws * wf * cell(s1, f1))
    return float(v), bool(sigma > S_SIG[-1] or sigma < S_SIG[0])


# ------------------------------------------------------------------------- pool gap filling
_POOL = {}


def pool_torsions(pdbid):
    if pdbid not in _POOL:
        u = I.load_univ(pdbid)
        p = I.pool_idx(u, I.K)
        _POOL[pdbid] = (np.asarray(u["PHI"][p], float), np.asarray(u["PSI"][p], float),
                        np.asarray(u["nat_ca"], float))
    return _POOL[pdbid]


def fill_pool(pdbid, phi, psi, mask):
    """S12's `fill_pool`: the pool member agreeing best with the KEPT residues supplies gaps."""
    PHI, PSI, nat = pool_torsions(pdbid)
    if not mask.any():
        mask = mask.copy(); mask[0] = True
    ag = torsion_agreement(PHI, PSI, phi, psi, mask)
    d = int(np.argmin(ag))
    return np.where(mask, phi, PHI[d]), np.where(mask, psi, PSI[d]), nat


# =============================================================================== per-arm eval
def eval_arm(arm, nat, inc, thresholds=(1e9, 90, 80, 70, 60, 50, 40, 30, 25, 20, 15)):
    z = np.load(post_path(arm))
    tg = I.targets()
    rows = []
    for t in tg:
        p = z[t["pdb"]]
        n = t["n"]
        d = T.decode(p)
        nphi, npsi = nat[t["pdb"]]
        mp, ms = determined_mask(n)
        ephi = T.ang_err(d["phi"], nphi); epsi = T.ang_err(d["psi"], npsi)
        err = np.concatenate([ephi[mp], epsi[ms]])
        # per-residue error used for the gate diagnostics (both angles where determined)
        per_res = np.sqrt((np.where(mp, ephi, np.nan) ** 2 + np.where(ms, epsi, np.nan) ** 2)
                          / 2.0)
        per_res = np.where(np.isnan(per_res), np.where(mp, ephi, epsi), per_res)
        ab_hat = T.abego4(d["phi"], d["psi"]); ab_nat = T.abego4(nphi, npsi)
        # S12-comparable accuracy: argmax of the posterior mass aggregated into the four
        # ABEGO classes.  This is what S12's 4-way classifier reported; the point-estimate
        # version above additionally pays for the decoder.
        ab_mass = np.stack([P_[:, _ABEGO_CELLS == c].sum(1) for c in range(4)], 1) \
            if (P_ := np.asarray(p, float)) is not None else None
        ab_arg = np.argmax(ab_mass, 1)
        nl = T.nll_native(p, nphi, npsi)
        # -- build, full coverage
        phi_b = d["phi"].copy(); psi_b = d["psi"].copy()
        phi_b[0] = nphi[0]; psi_b[-1] = npsi[-1]      # placeholders: zero effect on the CA trace
        PHI, PSI, ca_nat = pool_torsions(t["pdb"])
        rm_full = I.ca_rmsd(I.build_ca(phi_b, psi_b), ca_nat)
        # -- confidence-gated builds (deployable: the gate is the model's own sigma_hat)
        gate = {}
        for th in thresholds:
            m = d["sigma_hat"] <= th
            ph, ps, _ = fill_pool(t["pdb"], phi_b, psi_b, m)
            gate[f"{th:g}"] = {"cov": float(m.mean()),
                               "rmsd": I.ca_rmsd(I.build_ca(ph, ps), ca_nat),
                               "sig_kept": float(np.sqrt(np.mean(per_res[m] ** 2))) if m.any() else None}
        # -- ORACLE gate: keep residues whose TRUE error is small.  Separates "the predictor
        #    has no accurate residues" from "the predictor cannot tell which ones are".
        ogate = {}
        for th in (12.0, 20.0, 30.0, 45.0, 60.0):
            m = per_res <= th
            ph, ps, _ = fill_pool(t["pdb"], phi_b, psi_b, m)
            ogate[f"{th:g}"] = {"cov": float(m.mean()),
                                "rmsd": I.ca_rmsd(I.build_ca(ph, ps), ca_nat),
                                "sig_kept": float(np.sqrt(np.mean(per_res[m] ** 2))) if m.any() else None}
        for q in (0.9, 0.75, 0.5, 0.25):
            k = max(1, int(round(q * n)))
            m = np.zeros(n, bool); m[np.argsort(per_res)[:k]] = True
            ph, ps, _ = fill_pool(t["pdb"], phi_b, psi_b, m)
            ogate[f"q{q:g}"] = {"cov": float(m.mean()),
                                "rmsd": I.ca_rmsd(I.build_ca(ph, ps), ca_nat),
                                "sig_kept": float(np.sqrt(np.mean(per_res[m] ** 2)))}
        rows.append({
            "pdb": t["pdb"], "n": n, "fold": t["fold"], "base": inc[t["pdb"]],
            "sigma_eff": float(np.sqrt(np.mean(err ** 2))),
            "mae": float(np.mean(err)), "median_err": float(np.median(err)),
            "abego4_acc": float((ab_hat == ab_nat).mean()),
            "abego4_argmax_acc": float((ab_arg == ab_nat).mean()),
            "nll": float(np.mean(nl)),
            "sigma_hat_mean": float(np.mean(d["sigma_hat"])),
            "mass30_mean": float(np.mean(d["mass30"])),
            "rmsd_build": float(rm_full), "gate": gate, "ogate": ogate,
            "frac_err_under_12": float((per_res <= 12).mean()),
            "frac_err_under_20": float((per_res <= 20).mean()),
            "frac_err_under_30": float((per_res <= 30).mean()),
            "frac_err_under_45": float((per_res <= 45).mean()),
            "err_phi": [float(x) for x in ephi[mp]], "err_psi": [float(x) for x in epsi[ms]],
            "sigma_hat": [float(x) for x in d["sigma_hat"]],
            "per_res_err": [float(x) for x in per_res],
        })
    return rows


def arm_report(arm, rows):
    names = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    isf = np.array([r["pdb"] in set(I.FAIL18) for r in rows])
    base = np.array([r["base"] for r in rows], float)
    build = np.array([r["rmsd_build"] for r in rows], float)
    sig = np.array([r["sigma_eff"] for r in rows], float)
    out = {"arm": arm, "n": len(rows),
           "sigma_eff_mean": float(sig.mean()), "sigma_eff_median": float(np.median(sig)),
           "sigma_eff_FAIL18": float(sig[isf].mean()), "sigma_eff_other108": float(sig[~isf].mean()),
           "mae_mean": float(np.mean([r["mae"] for r in rows])),
           "abego4_acc": float(np.mean([r["abego4_acc"] for r in rows])),
           "abego4_acc_FAIL18": float(np.mean([r["abego4_acc"] for r in rows if r["pdb"] in set(I.FAIL18)])),
           "abego4_argmax_acc": float(np.mean([r["abego4_argmax_acc"] for r in rows])),
           "abego4_argmax_acc_FAIL18": float(np.mean([r["abego4_argmax_acc"] for r in rows
                                                      if r["pdb"] in set(I.FAIL18)])),
           "nll": float(np.mean([r["nll"] for r in rows])),
           "build": {**I.summary(build),
                     "frac_under_1.5": float((build < 1.5).mean()),
                     "FAIL18": float(build[isf].mean()), "other108": float(build[~isf].mean()),
                     "paired_vs_incumbent": I.paired(build, base, folds=folds, names=names)},
           "frac_err_under": {k: float(np.mean([r[f"frac_err_under_{k}"] for r in rows]))
                              for k in (12, 20, 30, 45)},
           "gates": {}, "oracle_gates": {}}
    for key, dst in (("gate", "gates"), ("ogate", "oracle_gates")):
        for th in rows[0][key]:
            g = np.array([r[key][th]["rmsd"] for r in rows], float)
            cov = np.array([r[key][th]["cov"] for r in rows], float)
            sk = np.array([r[key][th]["sig_kept"] if r[key][th]["sig_kept"] is not None else np.nan
                           for r in rows], float)
            skm = float(np.nanmean(sk)) if np.isfinite(sk).any() else float("nan")
            pv, pe = surface_readoff(skm if np.isfinite(skm) else 20.0, float(cov.mean()))
            out[dst][th] = {"coverage": float(cov.mean()), "sigma_kept": skm,
                            "rmsd": float(g.mean()), "median": float(np.median(g)),
                            "FAIL18": float(g[isf].mean()), "other108": float(g[~isf].mean()),
                            "frac_under_2.0": float((g < 2.0).mean()),
                            "frac_under_1.5": float((g < 1.5).mean()),
                            "surface_pred": float(pv), "surface_extrapolated": bool(pe),
                            "paired_vs_incumbent": I.paired(g, base, folds=folds, names=names)}
    sp, se = surface_readoff(sig.mean(), 1.0)
    out["surface_pred_full_cov"] = {"sigma": float(sig.mean()), "coverage": 1.0,
                                    "pred_rmsd": sp, "extrapolated": bool(se)}
    return out


# ============================================== ORACLE: extended sigma sweep + error shape
def oracle_sigma_sweep(nat, sigmas=(0, 6, 12, 20, 30, 40, 50, 60, 70, 80, 90, 100), seeds=SEEDS):
    """ORACLE DIAGNOSTIC.  S12's surface stops at 20 deg; a sequence-only predictor lands
    well outside it, so the axis is extended here with the identical noise model and the
    identical builder, at full coverage."""
    tg = I.targets()
    out = {}
    for sg in sigmas:
        vals = []
        for t in tg:
            nphi, npsi = nat[t["pdb"]]
            _, _, ca = pool_torsions(t["pdb"])
            acc = []
            for s in seeds:
                rng = np.random.default_rng(int(1e6 * sg + s + hash(t["pdb"]) % 9973))
                ph = T.wrap(nphi + rng.normal(0, sg * D2R, t["n"]))
                ps = T.wrap(npsi + rng.normal(0, sg * D2R, t["n"]))
                acc.append(I.ca_rmsd(I.build_ca(ph, ps), ca))
            vals.append(float(np.mean(acc)))
        v = np.array(vals)
        isf = np.array([t["pdb"] in set(I.FAIL18) for t in tg])
        out[f"{sg:g}"] = {"mean": float(v.mean()), "median": float(np.median(v)),
                          "frac_under_2.0": float((v < 2.0).mean()),
                          "FAIL18": float(v[isf].mean()), "other108": float(v[~isf].mean()),
                          "per_target": {t["pdb"]: float(x) for t, x in zip(tg, v)}}
        print(f"  ORACLE sigma={sg:g}: {v.mean():.3f}  <2A {(v<2.0).mean():.2f}", flush=True)
    return out


def oracle_errorshape(rows, nat, seeds=SEEDS, label=""):
    """ORACLE DIAGNOSTIC.  Native torsions corrupted by the PREDICTOR'S OWN empirical error
    distribution (errors pooled over all targets and resampled with a random sign), at full
    coverage.  Isolates error SHAPE from error SCALE: if this reproduces the predictor's
    measured build RMSD but the Gaussian surface at the same sigma does not, the shape is
    the explanation."""
    pool_err = np.concatenate([np.array(r["err_phi"] + r["err_psi"]) for r in rows])
    tg = I.targets()
    vals = []
    for t in tg:
        nphi, npsi = nat[t["pdb"]]
        _, _, ca = pool_torsions(t["pdb"])
        acc = []
        for s in seeds:
            rng = np.random.default_rng(int(7717 * s + hash(t["pdb"]) % 9973))
            e1 = rng.choice(pool_err, t["n"]) * rng.choice([-1.0, 1.0], t["n"])
            e2 = rng.choice(pool_err, t["n"]) * rng.choice([-1.0, 1.0], t["n"])
            acc.append(I.ca_rmsd(I.build_ca(T.wrap(nphi + e1 * D2R), T.wrap(npsi + e2 * D2R)), ca))
        vals.append(float(np.mean(acc)))
    v = np.array(vals)
    isf = np.array([t["pdb"] in set(I.FAIL18) for t in tg])
    return {"label": label, "sigma_of_shape": float(np.sqrt(np.mean(pool_err ** 2))),
            "mean": float(v.mean()), "median": float(np.median(v)),
            "frac_under_2.0": float((v < 2.0).mean()),
            "FAIL18": float(v[isf].mean()), "other108": float(v[~isf].mean())}


# ============================================================================== calibration
def calibration(rows, nbin=8):
    sh = np.concatenate([np.array(r["sigma_hat"]) for r in rows])
    er = np.concatenate([np.array(r["per_res_err"]) for r in rows])
    q = np.quantile(sh, np.linspace(0, 1, nbin + 1))
    out = []
    for a, b in zip(q[:-1], q[1:]):
        m = (sh >= a) & (sh <= b)
        if m.sum() < 5:
            continue
        out.append({"sigma_hat_lo": float(a), "sigma_hat_hi": float(b), "n": int(m.sum()),
                    "sigma_hat_mean": float(sh[m].mean()),
                    "actual_rms_err": float(np.sqrt(np.mean(er[m] ** 2))),
                    "actual_median_err": float(np.median(er[m]))})
    r = float(np.corrcoef(sh, er)[0, 1])
    from scipy.stats import spearmanr
    return {"bins": out, "pearson_r": r, "spearman": float(spearmanr(sh, er).statistic)}


# ==================================================================================== driver
def main(argv):
    # `n_libprior` (s13/tors_nulls.py) and `c_kmer` (s13/tors_support.py) are produced
    # outside the training driver but are scored on exactly the same footing.
    arms = argv or [a for a in list(ARMS) + ["n_libprior", "c_kmer"]
                    if os.path.exists(post_path(a))]
    nat = native_torsions()
    inc = incumbent()
    allrows, reports = {}, {}
    for a in arms:
        if not os.path.exists(post_path(a)):
            print(f"[missing] {a}"); continue
        t0 = time.time()
        rows = eval_arm(a, nat, inc)
        allrows[a] = rows
        rep = arm_report(a, rows)
        rep["calibration"] = calibration(rows)
        reports[a] = rep
        g0 = rep["gates"]["1e+09"] if "1e+09" in rep["gates"] else list(rep["gates"].values())[0]
        print(f"[{a}] sigma_eff={rep['sigma_eff_mean']:.1f} deg  abego4={rep['abego4_acc']:.3f}  "
              f"nll={rep['nll']:.3f}  build={rep['build']['mean']:.3f}  "
              f"surf_pred={rep['surface_pred_full_cov']['pred_rmsd']:.3f}  "
              f"({time.time()-t0:.0f}s)", flush=True)
    json.dump({"what": "predictor -> (sigma, coverage) -> emitted RMSD",
               "incumbent_mean": float(np.mean(list(inc.values()))),
               "reports": reports},
              open(os.path.join(T.RESULTS, "tors_eval.json"), "w"), indent=1, default=str)
    np.savez_compressed(os.path.join(T.CACHE, "tors_rows.npz"),
                        **{a: np.array(json.dumps(r)) for a, r in allrows.items()})

    # oracle diagnostics
    op = os.path.join(T.RESULTS, "tors_oracle.json")
    if not os.path.exists(op):
        print("ORACLE extended sigma sweep", flush=True)
        sw = oracle_sigma_sweep(nat)
        shapes = {}
        for a in ("p_grid", "p_point", "n_marg"):
            if a in allrows:
                shapes[a] = oracle_errorshape(allrows[a], nat, label=a)
                print(f"  ORACLE error-shape[{a}]: sigma={shapes[a]['sigma_of_shape']:.1f} "
                      f"-> {shapes[a]['mean']:.3f}", flush=True)
        json.dump({"what": "ORACLE DIAGNOSTIC: extended sigma sweep and predictor error shape",
                   "sweep": sw, "error_shape": shapes}, open(op, "w"), indent=1, default=str)
    print("wrote", os.path.join(T.RESULTS, "tors_eval.json"))


if __name__ == "__main__":
    main(sys.argv[1:])
