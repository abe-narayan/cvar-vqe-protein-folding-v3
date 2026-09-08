"""E5 mechanism -- interpolate between the PREDICTED and the ORACLE-native bin key.

For each rate r, a fraction r of the target's residues take the ORACLE bin and the rest
keep the predictor's argmax bin.  Locates exactly where between predicted and oracle the
pool gain (and, for a few rates, the emitted gain) appears.  Pool level is cheap; the
emitted version is run through s12/key_emit.py separately.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import key_lib as KL


def run(post_path, variant="full", alph="abego4",
        rates=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), reps=5, out="key_interp"):
    with open(post_path) as fh:
        P = {k: np.asarray(v, float) for k, v in json.load(fh)[variant].items()}
    nat = KL.native_torsions()
    rows = []
    for t in I.targets():
        pdb = t["pdb"]
        u = I.load_univ(pdb); rr = u["rr"]
        wb = KL.bins_of(alph, u["PHI"], u["PSI"])
        ob = KL.bins_of(alph, *nat[pdb]).astype(int)
        pb = P[pdb].argmax(1).astype(int)
        row = {"pdb": pdb, "fold": u["fold"], "f18": pdb in I.FAIL18,
               "blosum": KL.pool_stats(rr, KL.retrieve(u["sim"].astype(float)))}
        for r in rates:
            st, acc = [], []
            for k in range(1 if r in (0.0, 1.0) else reps):
                rng = np.random.default_rng(abs(hash((pdb, r, k))) % (2 ** 31))
                m = rng.random(len(ob)) < r
                b = np.where(m, ob, pb)
                acc.append(float((b == ob).mean()))
                st.append(KL.pool_stats(rr, KL.retrieve(KL.hard_key_score(wb, b))))
            row[f"r{r}"] = {"acc": float(np.mean(acc)),
                            **{m2: float(np.mean([s[m2] for s in st]))
                               for m2 in ("best", "mean", "band")}}
        rows.append(row)
        del u
    I.write(out, rows)
    keys = [k for k in rows[0] if k.startswith("r")]
    for g, sel in [("all", lambda r: True), ("FAIL18", lambda r: r["f18"]),
                   ("other108", lambda r: not r["f18"])]:
        rs = [r for r in rows if sel(r)]
        print(f"\n=== {g} n={len(rs)} BLOSUM best={np.mean([r['blosum']['best'] for r in rs]):.3f} "
              f"mean={np.mean([r['blosum']['mean'] for r in rs]):.3f} "
              f"band={np.mean([r['blosum']['band'] for r in rs]):.1f}")
        for k in keys:
            print(f"  {k:6s} acc={np.mean([r[k]['acc'] for r in rs]):.3f} "
                  f"best={np.mean([r[k]['best'] for r in rs]):.3f} "
                  f"mean={np.mean([r[k]['mean'] for r in rs]):.3f} "
                  f"band={np.mean([r[k]['band'] for r in rs]):.1f}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--post", default=os.path.join(I.RESULTS, "key_pred_post.json"))
    ap.add_argument("--variant", default="full")
    ap.add_argument("--out", default="key_interp")
    a = ap.parse_args()
    run(a.post, a.variant, out=a.out)
