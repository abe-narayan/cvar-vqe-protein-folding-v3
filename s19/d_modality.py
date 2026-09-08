"""AGENT D, Sprint 19 -- Blocks P1/P2/P3: attack the coordinator's multimodality census.

Pre-registered in `s19/PREREG_D.md`.  Reads only the 126 cached leave-fold-out distograms
(`s12/cache/disto_*.npz`) and the native CA traces (ORACLE, scoring only).  Touches nothing
in `results/`.

P1  the census, and whether it survives a density correction and a re-binning
P2  is multimodality HARMFUL at the aim point, matched within (separation x sd) cells
P3  is it COHERENT across pairs sharing a residue, against a stratified permutation null

Run:  python -m s19.d_modality
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
    os.environ[_v] = "1"

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402

OUT = os.path.join(HERE, "results", "D_P1P2P3_modality")
os.makedirs(OUT, exist_ok=True)

BIN_EDGES = np.array([4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0,
                      11.0, 12.5, 14.0, 16.0, 19.0, 23.0])
CENTRES = np.concatenate([[BIN_EDGES[0] - 0.5],
                          0.5 * (BIN_EDGES[1:] + BIN_EDGES[:-1]),
                          [BIN_EDGES[-1] + 2.0]])
#: nominal widths; the two catch-all bins get the width their CENTRES convention implies
WIDTHS = np.concatenate([[1.0], np.diff(BIN_EDGES), [4.0]])
NB = len(CENTRES)


# ------------------------------------------------------------------ mode detectors
def local_max(v):
    """(npairs, NB) bool: strict interior local maxima, endpoints compared one-sided."""
    m = np.zeros(v.shape, bool)
    m[:, 1:-1] = (v[:, 1:-1] > v[:, :-2]) & (v[:, 1:-1] > v[:, 2:])
    m[:, 0] = v[:, 0] > v[:, 1]
    m[:, -1] = v[:, -1] > v[:, -2]
    return m


def mass_modes(p, mass_gate=0.02):
    return local_max(p) & (p >= mass_gate)


def dens_modes(p, mass_gate=0.02):
    return local_max(p / WIDTHS[None]) & (p >= mass_gate)


def prom_modes(p, mass_gate=0.05, frac=0.5):
    """Density modes surviving a prominence gate: the valley between two kept modes must
    fall below `frac` x the smaller peak.  Greedy left-to-right merge, weakest peak drops."""
    dens = p / WIDTHS[None]
    cand = local_max(dens) & (p >= mass_gate)
    out = np.zeros(cand.shape, bool)
    for r in range(p.shape[0]):
        idx = list(np.flatnonzero(cand[r]))
        if not idx:
            out[r] = cand[r]
            continue
        changed = True
        while changed and len(idx) > 1:
            changed = False
            for a in range(len(idx) - 1):
                lo, hi = idx[a], idx[a + 1]
                valley = dens[r, lo:hi + 1].min()
                if valley >= frac * min(dens[r, lo], dens[r, hi]):
                    drop = lo if dens[r, lo] < dens[r, hi] else hi
                    idx.remove(drop)
                    changed = True
                    break
        out[r, idx] = True
    return out


def rebin(p, offset):
    """Merge adjacent bins pairwise.  offset 0: (0,1)(2,3)...  offset 1: (0)(1,2)(3,4)..."""
    groups = []
    b = 0
    if offset == 1:
        groups.append([0])
        b = 1
    while b < NB:
        groups.append(list(range(b, min(b + 2, NB))))
        b += 2
    P = np.stack([p[:, g].sum(1) for g in groups], 1)
    W = np.array([WIDTHS[g].sum() for g in groups])
    C = np.array([float((CENTRES[g] * WIDTHS[g]).sum() / WIDTHS[g].sum()) for g in groups])
    return P, W, C


def local_max_gen(v):
    m = np.zeros(v.shape, bool)
    if v.shape[1] > 2:
        m[:, 1:-1] = (v[:, 1:-1] > v[:, :-2]) & (v[:, 1:-1] > v[:, 2:])
    m[:, 0] = v[:, 0] > v[:, 1]
    m[:, -1] = v[:, -1] > v[:, -2]
    return m


def rebin_modes(p, offset, mass_gate=0.02):
    P, W, C = rebin(p, offset)
    return (local_max_gen(P / W[None]) & (P >= mass_gate)), C


# ------------------------------------------------------------------ per-target work
def bayes_median(risk, grid):
    """The selection path's own aim point: argmin_x sum_b p_b |x - c_b| (positive scaling of
    the cached `risk` does not move the argmin).  EXACT reconstruction of the deployed L1 risk."""
    return grid[np.argmin(risk, axis=1)]


def per_target(pdb, seq, n, fold):
    z = np.load(os.path.join(ROOT, "s12", "cache", f"disto_{pdb}.npz"))
    p = np.asarray(z["prob"], float)
    i, j = np.asarray(z["i"], int), np.asarray(z["j"], int)
    exp = np.asarray(z["expected"], float)
    sd = np.asarray(z["sd"], float)
    risk = np.asarray(z["risk"], float)
    grid = np.asarray(z["grid"], float)
    u = np.load(os.path.join(ROOT, "s8", "generate_univ", f"{pdb}.npz"), allow_pickle=True)
    nat = np.asarray(u["nat_ca"], float)
    dtrue = np.sqrt(((nat[i] - nat[j]) ** 2).sum(1))
    sep = (j - i).astype(float)

    M_mass = mass_modes(p)
    M_dens = dens_modes(p)
    M_prom = prom_modes(p)
    M_r0, C0 = rebin_modes(p, 0)
    M_r1, C1 = rebin_modes(p, 1)

    nmass, ndens, nprom = M_mass.sum(1), M_dens.sum(1), M_prom.sum(1)
    nr0, nr1 = M_r0.sum(1), M_r1.sum(1)

    dmed = bayes_median(risk, grid)

    # aim points from the density modes (the honest detector)
    big = np.where(M_dens.any(1), 0.0, 1.0)          # pairs with no detected mode
    mode_pos = np.where(M_dens, CENTRES[None], np.nan)
    with np.errstate(invalid="ignore"):
        # highest-density mode
        dens = np.where(M_dens, p / WIDTHS[None], -np.inf)
        mode1 = CENTRES[np.argmax(dens, axis=1)]
        # nearest mode to the truth (ORACLE)
        dd = np.abs(mode_pos - dtrue[:, None])
        best = CENTRES[np.nanargmin(np.where(np.isnan(dd), np.inf, dd), axis=1)]
    mode1 = np.where(big > 0, exp, mode1)
    best = np.where(big > 0, exp, best)

    e_mean = np.abs(exp - dtrue)
    e_med = np.abs(dmed - dtrue)
    e_mode1 = np.abs(mode1 - dtrue)
    e_best = np.abs(best - dtrue)
    delta = np.abs(best - exp)
    rng = SD.stable_rng(pdb, "D_flip")
    sgn = rng.choice([-1.0, 1.0], size=len(exp))
    e_flip = np.abs(exp + sgn * delta - dtrue)

    # ---- P3 coherence ---------------------------------------------------------
    share = (i[:, None] == i[None, :]) | (i[:, None] == j[None, :]) | \
            (j[:, None] == i[None, :]) | (j[:, None] == j[None, :])
    np.fill_diagonal(share, False)
    flag = (ndens >= 2)
    sepbin = np.clip(np.searchsorted([2, 4, 6, 8, 11, 999], sep, side="right") - 1, 0, 4)

    def assort(f):
        return float((f[:, None] & f[None, :])[share].mean()) if share.any() else np.nan

    obs = assort(flag)
    rng2 = SD.stable_rng(pdb, "D_perm")
    null = []
    for _ in range(500):
        g = flag.copy()
        for b in np.unique(sepbin):
            m = sepbin == b
            g[m] = rng2.permutation(flag[m])
        null.append(assort(g))
    null = np.asarray(null, float)
    z_assort = float((obs - null.mean()) / max(null.std(), 1e-12)) if np.isfinite(obs) else np.nan

    # signed coherence: does (best mode - mean) agree in sign between adjacent pairs?
    s_delta = np.sign(best - exp)
    s_resid = np.sign(dtrue - exp)

    def sign_agree(s):
        A = (s[:, None] * s[None, :]) > 0
        nz = (s[:, None] != 0) & (s[None, :] != 0)
        a = share & nz
        b = (~share) & nz
        np.fill_diagonal(b, False)
        return (float(A[a].mean()) if a.any() else np.nan,
                float(A[b].mean()) if b.any() else np.nan)

    da_adj, da_non = sign_agree(s_delta)
    ra_adj, ra_non = sign_agree(s_resid)

    return {
        "pdb": pdb, "n": int(n), "fold": int(fold), "npairs": int(len(i)),
        "f_mass2": float((nmass >= 2).mean()), "f_dens2": float((ndens >= 2).mean()),
        "f_prom2": float((nprom >= 2).mean()),
        "f_rebin0_2": float((nr0 >= 2).mean()), "f_rebin1_2": float((nr1 >= 2).mean()),
        "mean_nmass": float(nmass.mean()), "mean_ndens": float(ndens.mean()),
        "mass_within1A": float((p * (np.abs(CENTRES[None] - exp[:, None]) <= 1.0)).sum(1).mean()),
        "entropy": float((-(p * np.log(np.maximum(p, 1e-12))).sum(1)).mean()),
        "assort_obs": obs, "assort_null": float(null.mean()), "assort_z": z_assort,
        "delta_sign_adj": da_adj, "delta_sign_non": da_non,
        "resid_sign_adj": ra_adj, "resid_sign_non": ra_non,
        # target-level error means
        "e_mean": float(e_mean.mean()), "e_med": float(e_med.mean()),
        "e_mode1": float(e_mode1.mean()), "e_best": float(e_best.mean()),
        "e_flip": float(e_flip.mean()), "delta": float(delta.mean()),
        "_pair": {  # kept for the pooled matched-cell analysis
            "sep": sep, "sd": sd, "ndens": ndens, "nmass": nmass,
            "e_mean": e_mean, "e_med": e_med, "e_mode1": e_mode1,
            "e_best": e_best, "e_flip": e_flip, "dtrue": dtrue, "exp": exp,
        },
    }


def boot(d, rng, B=4000):
    d = np.asarray(d, float)
    k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def main():
    tg = I.targets()
    rows = [per_target(t["pdb"], t["seq"], t["n"], t["fold"]) for t in tg]
    rng = SD.stable_rng("D", "P1P2P3")
    g = lambda k: np.asarray([r[k] for r in rows], float)      # noqa: E731

    print(f"\nn = {len(rows)} targets, "
          f"{int(sum(r['npairs'] for r in rows))} pairs.\n")

    print("P1 -- THE CENSUS (target-level mean of the per-target multimodal fraction)")
    print(f"  {'detector':<26}{'mean':>8}{'[95% CI]':>22}{'median':>9}")
    for k, lab in [("f_mass2", "mass local max (>=2%)"),
                   ("f_dens2", "DENSITY local max"),
                   ("f_prom2", "density + prominence"),
                   ("f_rebin0_2", "rebinned x2, offset 0"),
                   ("f_rebin1_2", "rebinned x2, offset 1")]:
        m, lo, hi = boot(g(k), rng)
        print(f"  {lab:<26}{m:>8.3f}   [{lo:.3f},{hi:.3f}]{np.median(g(k)):>9.3f}")
    m, lo, hi = boot(g("mass_within1A"), rng)
    print(f"  {'mass within 1A of mean':<26}{m:>8.3f}   [{lo:.3f},{hi:.3f}]")
    m, lo, hi = boot(g("entropy"), rng)
    print(f"  {'predictive entropy':<26}{m:>8.3f}   [{lo:.3f},{hi:.3f}]  (uniform = "
          f"{np.log(NB):.3f})")

    # ---- pooled, matched-cell P2
    P = {k: np.concatenate([r["_pair"][k] for r in rows]) for k in rows[0]["_pair"]}
    tid = np.concatenate([np.full(r["npairs"], t) for t, r in enumerate(rows)])
    sepb = np.clip(np.searchsorted([2, 4, 6, 8, 11, 999], P["sep"], side="right") - 1, 0, 4)
    qs = np.quantile(P["sd"], [0.2, 0.4, 0.6, 0.8])
    sdb = np.searchsorted(qs, P["sd"])
    multi = P["ndens"] >= 2

    print("\nP2a -- is multimodality HARMFUL at the aim point?  |mean - true|, matched cells")
    print(f"  {'cell':<18}{'n_uni':>8}{'n_multi':>9}{'e_uni':>8}{'e_multi':>9}{'diff':>8}")
    diffs, wts = [], []
    for a in range(5):
        for b in range(5):
            m = (sepb == a) & (sdb == b)
            if m.sum() < 40 or multi[m].sum() < 10 or (~multi[m]).sum() < 10:
                continue
            eu = P["e_mean"][m & ~multi].mean()
            em = P["e_mean"][m & multi].mean()
            diffs.append(em - eu)
            wts.append(m.sum())
    diffs, wts = np.asarray(diffs), np.asarray(wts, float)
    print(f"  cells used: {len(diffs)};  weighted mean (multi - uni) = "
          f"{float((diffs * wts).sum() / wts.sum()):+.3f} A;  "
          f"unweighted {diffs.mean():+.3f};  cells worse for multi: "
          f"{int((diffs > 0).sum())}/{len(diffs)}")
    print(f"  UNMATCHED (the confounded comparison): "
          f"{P['e_mean'][multi].mean() - P['e_mean'][~multi].mean():+.3f} A")

    print("\nP2b -- aim points, restricted to MULTIMODAL pairs (ORACLE for e_best)")
    mm = multi
    for k, lab in [("e_mean", "mean (deployed)"), ("e_med", "Bayes L1 median"),
                   ("e_mode1", "highest-density mode"), ("e_best", "ORACLE nearest mode"),
                   ("e_flip", "matched-magnitude flip (CONTROL)")]:
        # target-level aggregation over multimodal pairs
        v = np.asarray([P[k][mm & (tid == t)].mean() if (mm & (tid == t)).any() else np.nan
                        for t in range(len(rows))])
        v = v[np.isfinite(v)]
        m_, lo, hi = boot(v, rng)
        print(f"  {lab:<34}{m_:>7.3f}  [{lo:.3f},{hi:.3f}]   n_t={len(v)}")
    vb = np.asarray([P["e_best"][mm & (tid == t)].mean() - P["e_mean"][mm & (tid == t)].mean()
                     for t in range(len(rows)) if (mm & (tid == t)).any()])
    vf = np.asarray([P["e_flip"][mm & (tid == t)].mean() - P["e_mean"][mm & (tid == t)].mean()
                     for t in range(len(rows)) if (mm & (tid == t)).any()])
    for lab, v in [("ORACLE best - mean", vb), ("random flip - mean (CONTROL)", vf),
                   ("ORACLE best - random flip", vb - vf)]:
        m_, lo, hi = boot(v, rng)
        print(f"  {lab:<34}{m_:+7.3f}  [{lo:+.3f},{hi:+.3f}]  W/L "
              f"{int((v < 0).sum())}/{int((v > 0).sum())}")

    print("\nP3 -- COHERENCE")
    for k, lab in [("assort_obs", "multimodal-multimodal adjacency"),
                   ("assort_null", "  stratified permutation null"),
                   ("assort_z", "  z of observed vs null")]:
        v = g(k)
        v = v[np.isfinite(v)]
        m_, lo, hi = boot(v, rng)
        print(f"  {lab:<38}{m_:>8.3f}  [{lo:.3f},{hi:.3f}]")
    for a, b, lab in [("delta_sign_adj", "delta_sign_non", "sign(bestmode-mean) agreement"),
                      ("resid_sign_adj", "resid_sign_non", "sign(true-mean) agreement")]:
        va, vb_ = g(a), g(b)
        ok = np.isfinite(va) & np.isfinite(vb_)
        m_, lo, hi = boot((va - vb_)[ok], rng)
        print(f"  {lab:<38} adj {va[ok].mean():.3f}  non-adj {vb_[ok].mean():.3f}   "
              f"diff {m_:+.3f} [{lo:+.3f},{hi:+.3f}]")

    out = {"n": len(rows), "rows": [{k: v for k, v in r.items() if k != "_pair"} for r in rows]}
    with open(os.path.join(OUT, "modality.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    np.savez_compressed(os.path.join(OUT, "pairs.npz"), tid=tid, sepb=sepb, sdb=sdb,
                        **{k: P[k] for k in P})
    with open(os.path.join(OUT, "COMPLETE"), "w") as fh:
        fh.write(f"n={len(rows)}\n")
    return rows


if __name__ == "__main__":
    main()
