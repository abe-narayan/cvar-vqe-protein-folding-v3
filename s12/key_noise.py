"""Pool-level corruption ladder for the ORACLE-native key: where does the gain die?

Independent re-measurement of the literature agent's noise study, but on the *native*
key (not the circular best-window key) and on all 126 targets.  Two error models:
  uniform  -- corrupted positions get a uniformly random bin
  biased   -- corrupted positions collapse to the string's modal bin (how a real
              predictor fails: it over-predicts the majority state)
  marginal -- corrupted positions are drawn from the library's global bin marginal
Also a fully random key (p=1) and a random key drawn from the marginal.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import key_lib as KL

MARGINAL = np.array([0.534, 0.396, 0.055, 0.016])      # library abego4 frequencies


def corrupt(bins, p, model, rng, nb=4):
    b = np.asarray(bins).copy()
    m = rng.random(len(b)) < p
    if not m.any():
        return b
    if model == "uniform":
        b[m] = rng.integers(0, nb, m.sum())
    elif model == "biased":
        b[m] = int(np.bincount(bins, minlength=nb).argmax())
    elif model == "marginal":
        b[m] = rng.choice(nb, m.sum(), p=MARGINAL[:nb] / MARGINAL[:nb].sum())
    return b


def run(alph="abego4", ps=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0),
        models=("uniform", "biased", "marginal"), reps=5, out="key_noise"):
    nat = KL.native_torsions()
    nb = KL.nbins(alph)
    rows = []
    for t in I.targets():
        pdb = t["pdb"]
        u = I.load_univ(pdb); rr = u["rr"]
        wb = KL.bins_of(alph, u["PHI"], u["PSI"])
        tb = KL.bins_of(alph, *nat[pdb]).astype(int)
        row = {"pdb": pdb, "fold": u["fold"], "f18": pdb in I.FAIL18,
               "univ_best": float(rr.min()),
               "blosum": KL.pool_stats(rr, KL.retrieve(u["sim"].astype(float)))}
        for mdl in models:
            for p in ps:
                acc, st = [], []
                for r in range(1 if p == 0.0 else reps):
                    rng = np.random.default_rng(hash((pdb, mdl, p, r)) % (2 ** 31))
                    cb = corrupt(tb, p, mdl, rng, nb)
                    acc.append(float((cb == tb).mean()))
                    st.append(KL.pool_stats(rr, KL.retrieve(KL.hard_key_score(wb, cb))))
                row[f"{mdl}_{p}"] = {"acc": float(np.mean(acc)),
                                     **{k: float(np.mean([s[k] for s in st]))
                                        for k in ("best", "mean", "band")}}
        rows.append(row)
        del u
        print(pdb, "done", flush=True)
    I.write(out, rows)
    return rows


def agg(rows):
    keys = [k for k in rows[0] if isinstance(rows[0][k], dict) and "acc" in rows[0][k]]
    out = {}
    for g, sel in [("all", lambda r: True), ("FAIL18", lambda r: r["f18"]),
                   ("other108", lambda r: not r["f18"])]:
        rs = [r for r in rows if sel(r)]
        out[g] = {"n": len(rs),
                  "blosum": {m: float(np.mean([r["blosum"][m] for r in rs]))
                             for m in ("best", "mean", "band")}}
        for k in keys:
            out[g][k] = {m: float(np.mean([r[k][m] for r in rs]))
                         for m in ("acc", "best", "mean", "band")}
    return out


if __name__ == "__main__":
    rows = run()
    A = agg(rows)
    I.write("key_noise_agg", A)
    keys = [k for k in rows[0] if isinstance(rows[0][k], dict) and "acc" in rows[0][k]]
    for g in ("all", "FAIL18", "other108"):
        print(f"\n=== {g} (n={A[g]['n']})  BLOSUM best={A[g]['blosum']['best']:.3f} "
              f"mean={A[g]['blosum']['mean']:.3f} band={A[g]['blosum']['band']:.1f}")
        for k in keys:
            v = A[g][k]
            print(f"  {k:18s} acc={v['acc']:.3f} best={v['best']:.3f} mean={v['mean']:.3f} band={v['band']:.1f}")
