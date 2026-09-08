"""E1 -- reproduce the literature agent's ORACLE ABEGO retrieval, and separate the two
very different oracles it conflates:

  ORACLE-BESTWIN : key = bin string of the universe's TRUE-BEST window   (what lit did)
  ORACLE-NATIVE  : key = bin string of the NATIVE structure's own torsions

Only ORACLE-NATIVE is the ceiling a sequence predictor could ever reach: a predictor
predicts the native's local conformation, not the identity of the best library window.
Also sweeps the alphabet resolution (abego4 / abego5 / km8 / km16).
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import key_lib as KL

ALPHS = ["ss3", "abego4", "abego5", "km8", "km12", "km16"]


def run(pdbs=None, alphs=ALPHS):
    tg = {t["pdb"]: t for t in I.targets()}
    nat = KL.native_torsions()
    pdbs = pdbs or list(tg)
    rows = []
    for pdb in pdbs:
        t = tg[pdb]
        u = I.load_univ(pdb)
        rr = u["rr"]
        p_bl = I.pool_idx(u)
        row = {"pdb": pdb, "n": t["n"], "fold": t["fold"], "f18": pdb in I.FAIL18,
               "nw": int(len(rr)), "univ_best": float(rr.min()),
               "blosum": KL.pool_stats(rr, p_bl)}
        best_w = int(np.argmin(rr))
        for a in alphs:
            wb = KL.bins_of(a, u["PHI"], u["PSI"])
            # ORACLE-BESTWIN
            k1 = KL.hard_key_score(wb, wb[best_w])
            row[f"bestwin_{a}"] = KL.pool_stats(rr, KL.retrieve(k1))
            # ORACLE-NATIVE
            if pdb in nat:
                nb = KL.bins_of(a, nat[pdb][0], nat[pdb][1])
                k2 = KL.hard_key_score(wb, nb)
                row[f"native_{a}"] = KL.pool_stats(rr, KL.retrieve(k2))
                row[f"agree_nat_bestwin_{a}"] = float((nb == wb[best_w]).mean())
                # hybrid: z(blosum) + z(bin match)
                kz = KL.zscore(k2) + KL.zscore(u["sim"].astype(float))
                row[f"nativehyb_{a}"] = KL.pool_stats(rr, KL.retrieve(kz))
        rows.append(row)
        del u
        print(pdb, row["blosum"]["best"], {a: row.get(f"native_{a}", {}).get("best") for a in alphs}, flush=True)
    return rows


def agg(rows, keys):
    out = {}
    for grp, sel in [("all", lambda r: True), ("FAIL18", lambda r: r["f18"]),
                     ("other108", lambda r: not r["f18"])]:
        rs = [r for r in rows if sel(r)]
        out[grp] = {"n": len(rs)}
        for k in keys:
            vals = [r[k] for r in rs if k in r]
            if not vals:
                continue
            out[grp][k] = {"best": float(np.mean([v["best"] for v in vals])),
                           "mean": float(np.mean([v["mean"] for v in vals])),
                           "band": float(np.mean([v["band"] for v in vals]))}
        out[grp]["univ_best"] = float(np.mean([r["univ_best"] for r in rs]))
        for a in ALPHS:
            k = f"agree_nat_bestwin_{a}"
            v = [r[k] for r in rs if k in r]
            if v:
                out[grp][k] = float(np.mean(v))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    tg = [t["pdb"] for t in I.targets()]
    pdbs = tg[: a.n] if a.n else tg
    rows = run(pdbs)
    keys = ["blosum"] + [f"{p}_{al}" for p in ("bestwin", "native", "nativehyb") for al in ALPHS]
    res = {"rows": rows, "agg": agg(rows, keys)}
    print(json.dumps(res["agg"], indent=1))
    I.write("key_oracle", res)
