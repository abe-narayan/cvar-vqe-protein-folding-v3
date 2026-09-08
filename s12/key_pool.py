"""Cheap pass: POOL statistics only (no projection) for many key arms at matched K=500.

Used to size the design before spending projection time.  Emitted RMSD is measured
separately by s12/key_emit.py -- pool statistics are NOT a substitute (S7-6).
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import key_lib as KL


def run(post_path, alph="abego4", mixes=(0.25, 0.5, 0.75), out="key_pool"):
    with open(post_path) as fh:
        raw = json.load(fh)
    post = {v: {k: np.asarray(x, float) for k, x in d.items()} for v, d in raw.items()}
    nat = KL.native_torsions()
    rows = []
    for t in I.targets():
        pdb = t["pdb"]
        u = I.load_univ(pdb); rr = u["rr"]
        wb = KL.bins_of(alph, u["PHI"], u["PSI"])
        sim = u["sim"].astype(float)
        arms = {"blosum": sim}
        nb = KL.bins_of(alph, *nat[pdb])
        arms["oracle_native"] = KL.hard_key_score(wb, nb)
        for v, d in post.items():
            if pdb not in d:
                continue
            s = KL.soft_key_score(wb, d[pdb])
            arms[f"pred:{v}"] = s
            arms[f"predhard:{v}"] = KL.hard_key_score(wb, d[pdb].argmax(1))
            if v == "full":
                for w in mixes:
                    arms[f"mix:{v}:{w}"] = (1 - w) * KL.zscore(sim) + w * KL.zscore(s)
        row = {"pdb": pdb, "n": u["n"], "fold": u["fold"], "f18": pdb in I.FAIL18,
               "univ_best": float(rr.min()), "arms": {}}
        for k, s in arms.items():
            row["arms"][k] = KL.pool_stats(rr, KL.retrieve(s))
        rows.append(row)
        del u
        print(pdb, " ".join(f"{k}={v['best']:.2f}" for k, v in row["arms"].items()), flush=True)
    I.write(out, rows)
    return rows


def agg(rows):
    arms = list(rows[0]["arms"])
    out = {}
    for g, sel in [("all", lambda r: True), ("FAIL18", lambda r: r["f18"]),
                   ("other108", lambda r: not r["f18"])]:
        rs = [r for r in rows if sel(r)]
        out[g] = {"n": len(rs)}
        for a in arms:
            out[g][a] = {m: float(np.mean([r["arms"][a][m] for r in rs]))
                         for m in ("best", "mean", "band")}
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--post", default=os.path.join(I.RESULTS, "key_pred_post.json"))
    ap.add_argument("--alph", default="abego4")
    ap.add_argument("--out", default="key_pool")
    a = ap.parse_args()
    rows = run(a.post, a.alph, out=a.out)
    A = agg(rows)
    arms = list(rows[0]["arms"])
    for m in ("best", "mean", "band"):
        print(f"\n--- pool {m} ---")
        for g in ("all", "FAIL18", "other108"):
            print(f"  {g:9s} " + "  ".join(f"{k}={A[g][k][m]:.3f}" for k in arms))
    I.write(a.out + "_agg", A)
