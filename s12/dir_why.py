"""D5b. Why is the head WORSE than random signs at the head's own accuracy?

The ladder (dir_emitavg) says the trained head, at 0.687 per-pair accuracy, emits +0.291 A
against the incumbent while ORACLE signs corrupted at exactly the same accuracy emit
-0.142 A.  Same accuracy, 0.43 A apart.  The hypothesis is the one obj_FINDINGS 1 already
established for the distogram itself: **coherent error is far more damaging than i.i.d.
error at the same error rate**, because a coherent field moves the whole structure while
independent flips cancel in the coordinate average.

Measured per target, for the head's wrongness field and for a random field of the same
rate (5 draws):
  err_rate        fraction of pairs whose predicted sign is wrong
  agree_share     P(two WRONG pairs that share a residue) / chance  -- coherence of the
                  mistakes (chance = err_rate, so 1.0 = i.i.d.)
  top1            Frobenius share of the leading eigenvector of the symmetric n x n
                  signed-mistake matrix  (i.i.d. ~ 2/n, a displacement field ~ 0.5)
  eff_rank        entropy participation rank of that matrix's eigenvalues
  scale_frac      share of the squared signed mistake in the single global "all pairs
                  wrong the same way" mode
Also the MAE and r_sep of the corrected objective for both, so that "same accuracy, same
MAE, different structure" can be stated directly.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")

from s12 import instrument as I
from s12 import obj_common as OC
from s12 import dir_common as DC

DELTA = 4.0


def field_stats(w, i, j, n):
    """w (npairs,) signed mistake: 0 where right, +-1 where wrong (sign = applied direction)."""
    m = np.abs(w) > 0
    rate = float(m.mean())
    E = np.zeros((n, n)); E[i, j] = w; E[j, i] = w
    ev = np.linalg.eigvalsh(E)
    a2 = ev ** 2
    tot = a2.sum()
    top1 = float(a2.max() / tot) if tot > 0 else float("nan")
    p = a2 / tot if tot > 0 else np.ones_like(a2) / len(a2)
    eff = float(np.exp(-(p * np.log(p + 1e-12)).sum()))
    ones = np.ones(n); ones = ones / np.linalg.norm(ones)
    scale = float((ones @ E @ ones) ** 2 / tot) if tot > 0 else float("nan")
    # coherence: among pairs sharing a residue, P(both wrong) / rate^2
    same = 0.0; both = 0.0
    for a in range(n):
        idx = np.where((i == a) | (j == a))[0]
        if len(idx) < 2:
            continue
        mm = m[idx].astype(float)
        k = len(idx)
        both += (mm.sum() ** 2 - (mm ** 2).sum()); same += k * (k - 1)
    coh = float((both / same) / max(rate ** 2, 1e-9)) if same > 0 else float("nan")
    return dict(err_rate=rate, top1=top1, eff_rank=eff, scale_frac=scale, coh=coh)


def main():
    Z = np.load(os.path.join(ROOT, "s12", "cache", "dir_head_pred.npz"))
    hd = json.load(open(os.path.join(ROOT, "s12", "results", "dir_head.json")))["per_target"]
    rows = []
    for k, t in enumerate(I.targets()):
        pdb, n = t["pdb"], t["n"]
        d = OC.load(pdb)
        exp, dtrue, sep = d["exp"], d["dtrue"], d["sep"]
        i, j = d["i"], d["j"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        sh = np.where(np.asarray(Z[f"{pdb}/full"], float) > 0.5, 1.0, -1.0)
        wh = np.where(sh != s_true, sh, 0.0)                 # signed mistake of the head
        acc = float(hd[pdb]["full"])
        st_h = field_stats(wh, i, j, n)
        st_h["mae"] = float(np.abs(DC.sign_target(exp, sh, DELTA) - dtrue).mean())
        rs = []
        for s in range(5):
            rng = np.random.default_rng(7000 + 13 * k + s)
            sr = DC.corrupt_sign(s_true, 1.0 - acc, rng)
            wr = np.where(sr != s_true, sr, 0.0)
            q = field_stats(wr, i, j, n)
            q["mae"] = float(np.abs(DC.sign_target(exp, sr, DELTA) - dtrue).mean())
            rs.append(q)
        st_r = {kk: float(np.mean([q[kk] for q in rs])) for kk in rs[0]}
        st_o = field_stats(np.where(s_true != s_true, s_true, 0.0), i, j, n)
        st_o["mae"] = float(np.abs(DC.sign_target(exp, s_true, DELTA) - dtrue).mean())
        rows.append(dict(pdb=pdb, n=n, fold=t["fold"], fail18=pdb in I.FAIL18,
                         mae_pt=float(np.abs(exp - dtrue).mean()),
                         head=st_h, rand=st_r, oracle=st_o))
    I.write("dir_why", rows)

    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    out = {}
    ks = ("err_rate", "coh", "top1", "eff_rank", "scale_frac", "mae")
    print(f"{'group':10s}{'arm':7s}" + "".join(f"{k:>12s}" for k in ks) + f"{'mae_pt':>9s}")
    for g, rs in grp.items():
        out[g] = {}
        for a in ("head", "rand"):
            out[g][a] = {k: float(np.nanmean([r[a][k] for r in rs])) for k in ks}
            print(f"{g:10s}{a:7s}" + "".join(f"{out[g][a][k]:12.4f}" for k in ks) +
                  f"{np.mean([r['mae_pt'] for r in rs]):9.3f}")
    I.write("dir_why_summary", out)


if __name__ == "__main__":
    main()
