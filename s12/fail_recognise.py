"""FAIL18 forensics, step 5: THE RECOGNITION QUESTION.

Inside the K=500 pool, what NATIVE-FREE signal separates the near-native band from the
rest (and in particular from what the shipped filter actually selects)?

Signals tested (all deployable: computed from the target sequence, the shipped
distogram, ESM-2 features of the target, and the candidate's own geometry):
  s_score      shipped distogram Bayes-risk score               (reference)
  s_short      the same score restricted to |i-j| in 2..4
  s_long       the same score restricted to |i-j| >= 5
  s_consensus  centrality: -(mean CA-RMSD to the other pool members)
  s_esmcon     agreement between the candidate's 8 A contact map and ESM-2 contacts
  s_legacy     Legacy 11-term knowledge-based energy (FITTED weights)
  s_rgfit      -|rg(cand) - rg_pred|, rg_pred from the distogram's own expected distances
  s_ssfit      agreement of the candidate's DSSP string with a propensity/ESM SS guess

Scored by AUC for discriminating band (rr <= pool_best + 1.5, ORACLE label used for
EVALUATION only) from non-band inside the pool, per target.  Null controls:
  N1 does the signal merely describe difficulty?  -> correlate per-target AUC with
     target difficulty (pool best) and report AUC after partialling it out
  N2 rg-partialling -> AUC of the signal among candidates in a narrow rg band
  N3 concentration  -> drop the 4 best targets
"""
from __future__ import annotations
import os, sys, json, math
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)
from s12 import instrument as I
from s12 import esm_bank
import protein_geometry as geo
from core import energy as en

HELIX_PROP = {  # Pace & Scholtz helix propensity (kcal/mol, lower = more helical)
    "A": 0.0, "L": 0.21, "R": 0.21, "M": 0.24, "K": 0.26, "Q": 0.39, "E": 0.40,
    "I": 0.41, "W": 0.49, "S": 0.50, "Y": 0.53, "F": 0.54, "H": 0.61, "V": 0.61,
    "N": 0.65, "T": 0.66, "C": 0.68, "D": 0.69, "G": 1.00, "P": 3.16}
SHEET_PROP = {  # Chou-Fasman P(beta)
    "V": 1.70, "I": 1.60, "Y": 1.47, "C": 1.19, "W": 1.37, "F": 1.38, "L": 1.30,
    "T": 1.19, "M": 1.05, "A": 0.83, "R": 0.93, "G": 0.75, "D": 0.54, "K": 0.74,
    "S": 0.75, "H": 0.87, "N": 0.89, "Q": 1.10, "P": 0.55, "E": 0.37}


def auc(score, label):
    """P(score of a positive > score of a negative); higher score = more band-like."""
    score = np.asarray(score, float); label = np.asarray(label, bool)
    if label.all() or not label.any():
        return float("nan")
    from scipy.stats import rankdata
    r = rankdata(score)
    n1 = label.sum(); n0 = len(label) - n1
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def rg_batch(W):
    W = np.asarray(W, float)
    c = W - W.mean(1, keepdims=True)
    return np.sqrt((c ** 2).sum(2).mean(1))


def contact_maps(W, thr=8.0):
    W = np.asarray(W, float)
    D = np.linalg.norm(W[:, :, None, :] - W[:, None, :, :], axis=-1)
    return (D < thr)


def signals(t, bank):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    u = I.load_univ(pdb); p = I.pool_idx(u)
    W = u["W"][p]; PHI = u["PHI"][p]; PSI = u["PSI"][p]
    rr = u["rr"][p]                                       # ORACLE label
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n); sep = j - i
    D = I.pair_dists(W, i, j).astype(np.float32).astype(float)
    grid, risk = dg["grid"], dg["risk"]
    g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    R = risk[np.arange(risk.shape[0])[None, :], g]        # (b, npairs)
    s_score = -R.mean(1)
    m_s = sep <= 4; m_l = sep >= 5
    s_short = -R[:, m_s].mean(1)
    s_long = -R[:, m_l].mean(1) if m_l.any() else np.zeros(len(W))

    # consensus / centrality
    P = I.pairwise_rmsd(W)
    s_cons = -P.mean(1)

    # ESM-2 contact agreement
    con = bank[seq][2]                                     # (n,n) probabilities
    mask = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    cm = contact_maps(W)
    w = con[mask]
    cw = cm[:, mask].astype(float)
    # agreement = correlation between the ESM probability and the binary contact
    s_esm = (cw * (w - w.mean())).sum(1) / max(np.sqrt(((w - w.mean()) ** 2).sum()), 1e-9)

    # Legacy energy
    bb = geo.build_backbone_batch(PHI, PSI)
    coords = {k: np.asarray(v, float) for k, v in bb.items()}
    comp = en.components_batch(seq, coords, PHI, PSI)
    s_leg = -np.asarray(en.totals_batch(comp), float)

    # rg predicted from the distogram itself (deployable)
    exp = np.asarray(dg["expected"], float)
    Dfull = np.zeros((n, n)); Dfull[i, j] = exp; Dfull[j, i] = exp
    k1 = np.arange(n - 1)
    Dfull[k1, k1 + 1] = 3.81; Dfull[k1 + 1, k1] = 3.81
    rg_pred = float(np.sqrt((Dfull ** 2).sum() / (2.0 * n * n)))
    rgs = rg_batch(W)
    s_rg = -np.abs(rgs - rg_pred)

    # SS agreement with a propensity-derived guess (no native, no oracle)
    hp = np.array([HELIX_PROP.get(c, 0.7) for c in seq])
    bp = np.array([SHEET_PROP.get(c, 1.0) for c in seq])
    kern = np.ones(5) / 5.0
    hs = np.convolve(1.0 - hp / 3.16, kern, mode="same")
    bs = np.convolve(bp / 1.7, kern, mode="same")
    want_H = hs > bs
    ss = [I.ss_of(PHI[a], PSI[a]) for a in range(len(W))]
    isH = np.array([[c == "H" for c in s] for s in ss], float)
    s_ss = (isH * want_H[None, :]).sum(1) / n + ((1 - isH) * (~want_H)[None, :]).sum(1) / n

    band = rr <= rr.min() + I.BAND
    S = dict(s_score=s_score, s_short=s_short, s_long=s_long, s_consensus=s_cons,
             s_esmcon=s_esm, s_legacy=s_leg, s_rgfit=s_rg, s_ssfit=s_ss)
    out = dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
               pool_best=float(rr.min()), band_frac=float(band.mean()),
               rg_pred=rg_pred, rg_nat=float(rg_batch(u["nat_ca"][None])[0]),
               auc={k: auc(v, band) for k, v in S.items()})

    # N2 rg-partialling: restrict to candidates within 0.5 A of the band's mean rg
    rg_band = float(rgs[band].mean())
    sel = np.abs(rgs - rg_band) < 0.5
    if band[sel].sum() >= 3 and (~band[sel]).sum() >= 10:
        out["auc_rgmatched"] = {k: auc(v[sel], band[sel]) for k, v in S.items()}
        out["n_rgmatched"] = int(sel.sum())
    # what the shipped filter picks vs the band, on each signal (z units)
    top75 = np.argsort(-s_score, kind="stable")[:I.M]
    out["z_gap_band_minus_top75"] = {
        k: float((v[band].mean() - v[top75].mean()) / (v.std() + 1e-9)) for k, v in S.items()}
    return out


def main():
    bank = esm_bank.load()
    ctrl = json.load(open(os.path.join(ROOT, "s12", "results", "fail_contrast.json")))["controls"]
    cset = {c["ctrl"] for c in ctrl}
    tg = [t for t in I.targets() if t["pdb"] in set(I.FAIL18) | cset]
    rows = []
    for k, t in enumerate(tg):
        rows.append(signals(t, bank))
        print(f"  {k+1}/{len(tg)} {t['pdb']} free={I.free_gb():.1f}", flush=True)
    I.write("fail_recognise", rows)

    f = [r for r in rows if r["fail18"]]; c = [r for r in rows if not r["fail18"]]
    keys = list(rows[0]["auc"].keys())
    print(f"\n== AUC (band vs rest, inside the K=500 pool); 0.5 = no skill ==")
    print(f"{'signal':14s} {'FAIL18':>8s} {'MATCH18':>8s} {'F drop4':>8s} {'F rg-matched':>13s} "
          f"{'F z(band-top75)':>16s}")
    summ = {}
    for k in keys:
        a = np.array([r["auc"][k] for r in f]); b = np.array([r["auc"][k] for r in c])
        d4 = np.sort(a)[:-4].mean()
        rgm = np.array([r["auc_rgmatched"][k] for r in f if "auc_rgmatched" in r])
        zg = np.array([r["z_gap_band_minus_top75"][k] for r in f])
        summ[k] = dict(fail18=float(a.mean()), match18=float(b.mean()),
                       fail18_drop4=float(d4),
                       fail18_rgmatched=float(rgm.mean()) if len(rgm) else None,
                       fail18_n_pos=int((a > 0.5).sum()),
                       z_gap=float(zg.mean()))
        print(f"{k:14s} {a.mean():8.3f} {b.mean():8.3f} {d4:8.3f} "
              f"{(rgm.mean() if len(rgm) else float('nan')):13.3f} {zg.mean():16.3f}")

    # N1: does AUC just describe difficulty?
    print("\n== N1 null: per-target AUC vs target difficulty (pool best), FAIL18+MATCH18 ==")
    pb = np.array([r["pool_best"] for r in rows])
    for k in keys:
        a = np.array([r["auc"][k] for r in rows])
        print(f"  {k:14s} r(AUC, pool_best) = {np.corrcoef(a, pb)[0,1]:+.3f}")

    # rg_pred quality
    e = np.array([r["rg_pred"] - r["rg_nat"] for r in rows])
    ef = np.array([r["rg_pred"] - r["rg_nat"] for r in f]); ec = np.array([r["rg_pred"] - r["rg_nat"] for r in c])
    print(f"\n== distogram-implied rg vs native rg: FAIL18 bias {ef.mean():+.3f} "
          f"MAE {np.abs(ef).mean():.3f} | MATCH18 bias {ec.mean():+.3f} MAE {np.abs(ec).mean():.3f}")
    I.write("fail_recognise_summary", dict(auc=summ, rg_pred_bias_fail=float(ef.mean()),
                                           rg_pred_mae_fail=float(np.abs(ef).mean()),
                                           rg_pred_bias_match=float(ec.mean()),
                                           rg_pred_mae_match=float(np.abs(ec).mean())))


if __name__ == "__main__":
    main()
