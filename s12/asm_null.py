"""ASM-7: the RIGID ladder on the capacity-null bank.

Same DP, same bank sizes, same minimum piece length as `asm_rigid`, but the windows are
synthetic: each piece's (phi, psi) are drawn i.i.d. from the pooled Ramachandran of the
same legal library and the piece is the ideal-geometry rebuild of those torsions.  The
bank therefore carries the correct marginal residue conformations and NOTHING ELSE -- no
real fragment, no sequence-structure relation, no residue-to-residue correlation.

If the synthetic bank reproduces the real bank's k>=2 ladder, the ladder measures the
capacity of 6k rigid parameters, not the library.

Usage: python -m s12.asm_null
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402
from s12.asm_rigid import dp_ladder                    # noqa: E402
from s12.asm_chain2 import synth_bank                  # noqa: E402

LMIN = 4


def main():
    tg = sorted(I.targets(), key=lambda t: (t["fold"], t["pdb"]))
    res = {}; cur = None; t0 = time.time()
    import s12.asm_chain2 as C2
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); C2._SYNTH.clear(); C2._POOL.clear(); cur = t["fold"]
        u = I.load_univ(t["pdb"]); nat = u["nat_ca"]; del u
        n = t["n"]
        cost = np.full((n, n), np.inf)
        for L in range(LMIN, n + 1):
            nwr = A.bank_ideal(t["fold"], L, lean=True)["nw"]
            b = synth_bank(t["fold"], L, nwr)
            segs = np.stack([nat[a:a + L] for a in range(n - L + 1)])
            M = A.msd_many(b["Ic"], b["In2"], segs) * L
            for a in range(n - L + 1):
                cost[a, L - 1] = M[:, a].min()
        d = dp_ladder(cost, n, kmax=4, lmin=LMIN)
        res[t["pdb"]] = {"n": n, "fold": t["fold"], "nw": int(nwr),
                         **{f"k{k}": d[k] for k in (1, 2, 3, 4)}}
        if (q + 1) % 20 == 0:
            print(f"{q+1}/{len(tg)} {time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)
    I.write("asm_null_rigid", res)
    f18 = np.array([t["pdb"] in I.FAIL18 for t in tg])
    print("\nRIGID ladder, CAPACITY-NULL bank (i.i.d. Ramachandran pieces, matched size)")
    for k in (1, 2, 3, 4):
        v = np.array([res[t["pdb"]][f"k{k}"] for t in tg], float)
        print(f"  k={k}  {np.nanmean(v):.3f}   FAIL18 {np.nanmean(v[f18]):.3f}  "
              f"other {np.nanmean(v[~f18]):.3f}")
    return res


if __name__ == "__main__":
    main()
