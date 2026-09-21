"""LANE L / L3a -- is the pool's common-mode error fraction length-dependent?

PREREG `s32/PREREG_S32_L.md` @ 88f2da39, L-H3a.  EXPLORATORY, declared as such.

THE LENGTH SUSPICION.  Project memory records that 68% of a retrieved pool's squared error
is a bias shared by all members (`pool-error-is-68-percent-common-mode`, exact identity,
verified to 2.7e-14 on 126 targets).  At 9-16 residues a whole target IS one retrieved
window, so every member is a single draw from one library and a shared bias is easy to
believe.  At 40-60 residues a target is still one window here, but drawn from a corpus
1,000x larger and 4x longer, and the question is whether the shared component survives.

THE STATISTIC, copied from the entry that produced it (contract rule 4).  With the m = 75
retained members written in the MEDOID FRAME THE AVERAGE IS ACTUALLY TAKEN IN as
w_k = t + e_k, the average c = mean_k w_k and d_k = w_k - c:

    mean_k |e_k|^2  =  |ebar|^2  +  mean_k |d_k|^2
                       COMMON        IDIOSYNCRATIC

    f = |ebar|^2 / mean_k |e_k|^2 ,   and the i.i.d. model predicts f = 1/m.

THE ENTRY'S OWN WARNING IS CARRIED WITH IT: *f is NOT a screen and must never be used as
one* -- it falls whenever a source is merely more scattered, and |ebar|^2 = n*RMSD^2
exactly, so the numerator is not independent evidence.  It is reported here descriptively,
to answer one question only: does the FRACTION move with length?

SELF-TEST (contract rule 5).  The decomposition is an exact algebraic identity ONLY if the
native is posed into the same frame the average was taken in.  A frame mismatch -- the
single most likely implementation error here -- breaks it.  `check` asserts the identity to
1e-9 on every target, so the test fails precisely when the bug is present.

    python -m s32.s32_L_commonmode run
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
M = 75
K = 500


def decompose(W, nat, m=M):
    """(f, common, idio, identity_residual) for the top-m of `W` in the medoid frame."""
    sys.path.insert(0, BASE)
    from s12 import instrument as I
    Ws = np.asarray(W, float)[:m]
    P = I.pairwise_rmsd(Ws)
    b = I.medoid(P)
    Wp = I.superpose_batch(Ws, Ws[b])          # the frame the average is taken in
    c = Wp.mean(0)
    t = I.superpose_batch(np.asarray(nat, float)[None], c)[0]   # native INTO that frame
    e = Wp - t
    ebar = c - t
    d = Wp - c
    me = float((e ** 2).sum((1, 2)).mean())
    common = float((ebar ** 2).sum())
    idio = float((d ** 2).sum((1, 2)).mean())
    return (common / me if me else np.nan, common, idio, abs(me - common - idio))


def run(verbose=True):
    sys.path.insert(0, BASE)
    from s12 import instrument as I
    from s32 import s32_L_ladder as LD
    from s32 import s32_L_corpus as cp

    out = {}
    worst = 0.0
    for kind in ("short", "long"):
        rows = []
        if kind == "short":
            tg = I.targets()
        else:
            tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
            cp.load()
        for k, t in enumerate(tg):
            if kind == "short":
                u = I.load_univ(t["pdb"])
                W = np.asarray(u["W"], float)[I.pool_idx(u, K)]
                nat = np.asarray(u["nat_ca"], float)
                n = t["n"]
            else:
                W, _, _ = LD._bank_long(t)
                from core import geometry as geo
                import glob
                p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
                if not os.path.exists(p):
                    p = glob.glob(os.path.join(BASE, "prots", t["pdb"].lower() + ".pdb"))[0]
                _, coords, _, _ = geo.native_coords_from_pdb(p)
                nat = np.asarray(coords["CA"], float); n = len(nat)
            f, com, idi, res = decompose(W, nat)
            worst = max(worst, res / max(1e-12, com + idi))
            rows.append(dict(pdb=t["pdb"], n=int(n), f=float(f),
                             common=com, idio=idi))
            if verbose and (k + 1) % 25 == 0:
                print(f"  {kind} {k+1}/{len(tg)}", flush=True)
        fs = np.array([r["f"] for r in rows])
        out[kind] = dict(n=len(rows), m=M,
                         f_mean=float(fs.mean()), f_median=float(np.median(fs)),
                         f_se=float(fs.std(ddof=1) / len(fs) ** 0.5),
                         f_p10=float(np.percentile(fs, 10)),
                         f_p90=float(np.percentile(fs, 90)),
                         iid_prediction=1.0 / M,
                         observed_over_iid=float(fs.mean() * M),
                         mean_len=float(np.mean([r["n"] for r in rows])), rows=rows)
    # SELF-TEST: the decomposition is exact, so a frame error must break it.
    assert worst < 1e-9, f"bias-variance identity violated at {worst:.2e} -- frame mismatch"
    out["identity_worst_relative_residual"] = float(worst)
    path = os.path.join(RESULTS, "L3a_commonmode.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: ({a: b for a, b in v.items() if a != "rows"}
                              if isinstance(v, dict) else v)
                          for k, v in out.items()}, indent=1))
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
