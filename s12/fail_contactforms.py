"""FAIL18 forensics, step 10 (coordinator request): which FORM of ESM-2 contact
agreement carries in-band skill -- and does ANY of them carry skill that a foreign
contact map does not?

Every form is scored twice: against the target's OWN ESM-2 contact map, and against a
length-matched FOREIGN target's map (the null).  The quantity that matters is the
DIFFERENCE: that is the part of the signal that is about this sequence rather than about
the shape of the statistic.  Judged on in-band AUC BEFORE any emitted-RMSD arm, as asked.

Forms (candidate contact map vs ESM probabilities p):
  corr3      corr(binary CA 8 A map, p) over |i-j| >= 3          [the E4 form]
  corr5      same over |i-j| >= 5
  cb5        same over |i-j| >= 5 using CB-CB at 8 A
  soft5      corr(sigmoid((8-d)/1.0), p) over |i-j| >= 5
  loglik5    sum c log p + (1-c) log(1-p) over |i-j| >= 5
  rank5      Spearman(-d, p) over |i-j| >= 5
  topk       fraction of the n highest-p pairs that are contacts (|i-j| >= 5)
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I
from s12.fail_esmrescore import load_con
import protein_geometry as geo

FORMS = ["corr3", "corr5", "cb5", "soft5", "loglik5", "rank5", "topk"]


def auc(score, label):
    from scipy.stats import rankdata
    label = np.asarray(label, bool)
    if label.all() or not label.any():
        return float("nan")
    r = rankdata(score); n1 = label.sum(); n0 = len(label) - n1
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _corr(X, w):
    wc = w - w.mean()
    Xc = X - X.mean(1, keepdims=True)
    return (Xc @ wc) / (np.sqrt((Xc ** 2).sum(1)) * np.sqrt((wc ** 2).sum()) + 1e-9)


def all_forms(Dca, Dcb, con, n):
    from scipy.stats import rankdata
    out = {}
    m3 = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    m5 = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 5
    p3, p5 = con[m3].astype(np.float64), con[m5].astype(np.float64)
    c3 = (Dca[:, m3] < 8.0).astype(np.float64)
    c5 = (Dca[:, m5] < 8.0).astype(np.float64)
    out["corr3"] = _corr(c3, p3)
    out["corr5"] = _corr(c5, p5) if m5.any() else np.zeros(len(Dca))
    out["cb5"] = _corr((Dcb[:, m5] < 8.0).astype(np.float64), p5) if m5.any() else np.zeros(len(Dca))
    s5 = 1.0 / (1.0 + np.exp(Dca[:, m5] - 8.0)) if m5.any() else np.zeros((len(Dca), 1))
    out["soft5"] = _corr(s5, p5) if m5.any() else np.zeros(len(Dca))
    if m5.any():
        q = np.clip(p5, 1e-4, 1 - 1e-4)
        out["loglik5"] = c5 @ np.log(q) + (1 - c5) @ np.log(1 - q)
        rp = rankdata(p5)
        out["rank5"] = _corr(np.apply_along_axis(rankdata, 1, -Dca[:, m5]), rp)
        k = min(n, len(p5))
        top = np.argsort(-p5)[:k]
        out["topk"] = c5[:, top].mean(1)
    else:
        for f in ("loglik5", "rank5", "topk"):
            out[f] = np.zeros(len(Dca))
    return out


def main():
    con_bank = load_con()
    tg = I.targets(); seqs = [t["seq"] for t in tg]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]
        band = rr <= rr.min() + I.BAND
        bb = geo.build_backbone_batch(u["PHI"][p], u["PSI"][p])
        CB = np.asarray(bb.get("CB", bb["CA"]), float)
        Dca = np.linalg.norm(Wc[:, :, None] - Wc[:, None], axis=-1)
        Dcb = np.linalg.norm(CB[:, :, None] - CB[:, None], axis=-1)
        alt = [s for s in seqs if len(s) == n and s != seq]
        own = all_forms(Dca, Dcb, con_bank[seq], n)
        nul = all_forms(Dca, Dcb, con_bank[alt[k % len(alt)]], n) if alt else own
        rows.append(dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                         auc_own={f: auc(own[f], band) for f in FORMS},
                         auc_null={f: auc(nul[f], band) for f in FORMS}))
        if (k + 1) % 30 == 0:
            print(f"  {k+1}/126 free={I.free_gb():.1f}", flush=True)
    I.write("fail_contactforms", rows)

    f = [r for r in rows if r["fail18"]]; o = [r for r in rows if not r["fail18"]]
    print(f"\n== in-band AUC by form; 'net' = own - foreign-map null (the sequence-specific part) ==")
    print(f"{'form':9s} | {'F own':>7s}{'F null':>7s}{'F net':>7s} | {'O own':>7s}{'O null':>7s}{'O net':>7s} "
          f"| {'F net ci95':>18s}")
    summ = {}
    rng = np.random.default_rng(0)
    for fm in FORMS:
        aw = np.array([r["auc_own"][fm] for r in f]); an = np.array([r["auc_null"][fm] for r in f])
        ow = np.array([r["auc_own"][fm] for r in o]); on = np.array([r["auc_null"][fm] for r in o])
        d = aw - an
        bs = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(4000)])
        ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
        summ[fm] = dict(fail_own=float(aw.mean()), fail_null=float(an.mean()),
                        fail_net=float(d.mean()), fail_net_ci=ci,
                        other_own=float(ow.mean()), other_null=float(on.mean()),
                        other_net=float((ow - on).mean()))
        print(f"{fm:9s} | {aw.mean():7.3f}{an.mean():7.3f}{d.mean():7.3f} | "
              f"{ow.mean():7.3f}{on.mean():7.3f}{(ow-on).mean():7.3f} | [{ci[0]:7.3f},{ci[1]:7.3f}]")
    I.write("fail_contactforms_summary", summ)


if __name__ == "__main__":
    main()
