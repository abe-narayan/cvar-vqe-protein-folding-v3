"""D3d. What has the sign head actually learned?

The head reaches ~0.68 accuracy against a ~0.53 constant baseline.  Two things could
produce that, and they have very different value:

  (i)  ONE BIT PER TARGET -- "this whole peptide is predicted too compact / too extended".
       The distogram's error is known to contain a large global scaling mode (15.7% of the
       squared error, obj_FINDINGS 1), so a head that only learns the target-level bias
       would score well on per-pair accuracy while carrying no pair-specific information.
  (ii) PAIR-SPECIFIC direction, which is what the objective actually needs.

Decomposition reported here:
  acc_const      global majority sign of the training folds
  acc_tmaj       ORACLE per-target majority sign (the ceiling of (i): one oracle bit/target)
  acc_head       the head
  acc_head_res   the head's accuracy on the pairs where the ORACLE sign DISAGREES with the
                 target's own majority -- pure (ii), chance = the minority fraction
  agree_htmaj    fraction of head predictions equal to the head's own per-target majority
Also: accuracy by |error| decile (is the head right where the error is BIG?), by separation
shell, and per target on FAIL18.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")

from s12 import instrument as I
from s12 import dir_head as DH

def main():
    z, keys, names = DH.load_bank()
    Z = np.load(os.path.join(ROOT, "s12", "cache", "dir_head_pred.npz"))
    tgset = {t["pdb"] for t in I.targets()}
    sepc = names.index("sep")
    rows = {}
    for pdb, seq, n, fold in keys:
        if pdb not in tgset:
            continue
        y = np.asarray(z[f"{pdb}/y"], int)
        err = np.asarray(z[f"{pdb}/err"], float)          # exp - dtrue  (ORACLE)
        sep = np.asarray(z[f"{pdb}/X"], float)[:, sepc]
        ph = np.asarray(Z[f"{pdb}/full"], float) > 0.5
        maj = int(y.mean() > 0.5)                          # ORACLE target majority
        hmaj = int(ph.mean() > 0.5)
        res = y != maj                                     # minority (pair-specific) pairs
        rows[pdb] = dict(n=int(n), fold=int(fold), fail18=pdb in set(I.FAIL18),
                         npairs=int(len(y)),
                         acc_head=float((ph == y).mean()),
                         acc_tmaj=float((y == maj).mean()),
                         acc_const=float((y == 0).mean()),
                         acc_head_res=float((ph[res] == y[res]).mean()) if res.any() else float("nan"),
                         minority_frac=float(res.mean()),
                         agree_htmaj=float((ph == bool(hmaj)).mean()),
                         tmaj_correct=int(hmaj == maj),
                         y=y.tolist(), ph=ph.astype(int).tolist(),
                         abserr=np.abs(err).tolist(), sep=sep.astype(int).tolist())

    def agg(sel, k):
        v = np.array([rows[p][k] for p in sel], float)
        return float(np.nanmean(v))

    groups = {"tuning126": list(rows),
              "fail18": [p for p in rows if rows[p]["fail18"]],
              "other108": [p for p in rows if not rows[p]["fail18"]]}
    out = {"groups": {}}
    print(f"{'group':11s}{'acc_head':>10s}{'acc_tmaj':>10s}{'acc_const':>10s}"
          f"{'acc_head_res':>13s}{'minor_frac':>12s}{'agree_htmaj':>13s}{'tmaj_ok':>9s}")
    for g, sel in groups.items():
        r = {k: agg(sel, k) for k in ("acc_head", "acc_tmaj", "acc_const", "acc_head_res",
                                      "minority_frac", "agree_htmaj", "tmaj_correct")}
        # pooled accuracy by |error| decile and by separation
        Y = np.concatenate([rows[p]["y"] for p in sel])
        P = np.concatenate([rows[p]["ph"] for p in sel])
        E = np.concatenate([rows[p]["abserr"] for p in sel])
        S = np.concatenate([rows[p]["sep"] for p in sel])
        q = np.quantile(E, np.linspace(0, 1, 11))
        dec = {}
        for k in range(10):
            m = (E >= q[k]) & (E <= q[k + 1] if k == 9 else E < q[k + 1])
            dec[f"d{k}"] = dict(n=int(m.sum()), lo=float(q[k]), hi=float(q[k + 1]),
                                acc=float((P[m] == Y[m]).mean()) if m.any() else float("nan"))
        r["acc_by_abserr_decile"] = dec
        r["acc_by_sep"] = {int(s): dict(n=int((S == s).sum()),
                                        acc=float((P[S == s] == Y[S == s]).mean()))
                           for s in np.unique(S)}
        # error-WEIGHTED accuracy: does the head get the pairs that matter right?
        r["err_weighted_acc"] = float((E * (P == Y)).sum() / E.sum())
        out["groups"][g] = r
        print(f"{g:11s}{r['acc_head']:10.4f}{r['acc_tmaj']:10.4f}{r['acc_const']:10.4f}"
              f"{r['acc_head_res']:13.4f}{r['minority_frac']:12.4f}"
              f"{r['agree_htmaj']:13.4f}{r['tmaj_correct']:9.3f}")
    print("\nhead accuracy by |error| decile (tuning126):")
    for k, v in out["groups"]["tuning126"]["acc_by_abserr_decile"].items():
        print(f"  {k} |err| {v['lo']:5.2f}-{v['hi']:5.2f}  n={v['n']:6d}  acc={v['acc']:.4f}")
    print("\nhead accuracy by separation (tuning126):")
    for s, v in out["groups"]["tuning126"]["acc_by_sep"].items():
        print(f"  sep {s:2d}  n={v['n']:6d}  acc={v['acc']:.4f}")
    out["per_target"] = {p: {k: rows[p][k] for k in
                             ("n", "fold", "fail18", "npairs", "acc_head", "acc_tmaj",
                              "acc_const", "acc_head_res", "minority_frac", "tmaj_correct")}
                         for p in rows}
    I.write("dir_headanal", out)


if __name__ == "__main__":
    main()
