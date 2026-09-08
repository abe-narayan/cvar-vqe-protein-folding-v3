"""S15 AUDIT / Part A -- is the instrument deterministic under THREAD-COUNT changes?

The distogram is a torch forward pass.  A same-session repeat is bit-identical, but three
different `torch.set_num_threads` values produce three different bit patterns of `prob`.
This module prices that: it recomputes every one of the 126 distograms at a given thread
count and reports the `shipped` constant it implies, so the drift can be compared against
the pinned 3.4540004952559396.

    python -m s15.audit_determinism 1
    python -m s15.audit_determinism 2
    python -m s15.audit_determinism 8
"""
from __future__ import annotations
import os, sys, json, time, hashlib
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)

PINNED = {"shipped": 3.4540004952559396, "pool_best": 1.7108244199364904,
          "top75_best": 2.3061526409453816, "synthesis_fit": 3.2040761603809194,
          "n_zero_recall": 18}


def run(threads):
    import torch
    torch.set_num_threads(int(threads))
    from core import pipeline as pl
    from core import predict as dgm
    from s12 import instrument as I
    tg = I.targets()
    rows, shipped, flips = [], [], []
    t0 = time.time()
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        z = np.load(os.path.join(I.CACHE, f"disto_{pdb}.npz"))
        cached = {kk: z[kk] for kk in z.files}
        pl.guard_esm([t["seq"]])
        model = pl.fold_model(int(t["fold"]))
        d = dgm.Distogram.for_target(t["seq"], model=model)
        fresh = dict(cached)
        fresh["risk"] = np.asarray(d._risk, np.float32)
        fresh["prob"] = np.asarray(d.prob, np.float32)
        fresh["expected"] = np.asarray(d.expected, np.float32)
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        i, j = I.pair_index(t["n"])
        D = I.pair_dists(u["W"][p], i, j).astype(np.float32).astype(float)
        rr = u["rr"][p]
        a = int(np.argmin(I.shipped_score(cached, D)))
        b = int(np.argmin(I.shipped_score(fresh, D)))
        shipped.append(float(rr[b]))
        if a != b:
            flips.append({"pdb": pdb, "cached_idx": a, "fresh_idx": b,
                          "rmsd_cached": float(rr[a]), "rmsd_fresh": float(rr[b])})
        rows.append({"pdb": pdb, "argmin_same": a == b,
                     "rmsd_cached": float(rr[a]), "rmsd_fresh": float(rr[b]),
                     "prob_max_abs_diff": float(np.abs(fresh["prob"].astype(np.float64)
                                                       - cached["prob"].astype(np.float64)).max()),
                     "risk_max_abs_diff": float(np.abs(fresh["risk"].astype(np.float64)
                                                       - cached["risk"].astype(np.float64)).max()),
                     "prob_sha16": hashlib.sha256(fresh["prob"].tobytes()).hexdigest()[:16]})
        if (k + 1) % 20 == 0:
            print(f"[{k+1}/{len(tg)}] running shipped {np.mean(shipped):.6f}", flush=True)
    out = {"torch_threads": int(threads), "n_targets": len(rows),
           "seconds": round(time.time() - t0, 1),
           "shipped_recomputed": float(np.mean(shipped)),
           "shipped_pinned": PINNED["shipped"],
           "shipped_delta": float(np.mean(shipped)) - PINNED["shipped"],
           "n_argmin_flips": len(flips), "flips": flips,
           "max_prob_abs_diff": max(r["prob_max_abs_diff"] for r in rows),
           "max_risk_abs_diff": max(r["risk_max_abs_diff"] for r in rows),
           "rows": rows}
    p = os.path.join(OUT, f"audit_determinism_t{threads}.json")
    with open(p, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in out if k not in ("rows", "flips")}, indent=1))
    print("flips:", json.dumps(flips))
    print("wrote", p)
    return out


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
