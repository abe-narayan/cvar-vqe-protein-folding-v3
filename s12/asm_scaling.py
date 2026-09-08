"""ASM-6: is assembly better than simply having MORE whole windows?

The k-piece oracle floor is a minimum over ~(bank size)^k combinations, so it must be
priced against the coverage-scaling curve of a SINGLE window bank.  For every target we
measure E[ min RMSD | s candidates ] for

    k=1  : s whole windows drawn uniformly from the length-n bank          (effective = s)
    k=2  : s pieces per side, rigid-placement DP over cuts                 (effective ~ C2*s^2)
    k=3  : s pieces per slot, rigid-placement DP over cut pairs            (effective ~ C3*s^3)

where C2/C3 count the admissible cut positions.  If the three families collapse onto ONE
curve against log10(effective candidate count), assembly is buying nothing that raw
candidate count would not.  Everything here is ORACLE.

Usage: python -m s12.asm_scaling
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

LMIN = 4
S1 = [30, 100, 300, 1000, 3000, 10000]
SK = [10, 30, 100, 300, 1000]
R = 6


def main():
    tg = sorted(I.targets(), key=lambda t: (t["fold"], t["pdb"]))
    res = {}; t0 = time.time(); cur = None
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); cur = t["fold"]
        u = I.load_univ(t["pdb"]); nat = u["nat_ca"]; rr = u["rr"]; nw = len(rr); del u
        n = t["n"]
        rng = np.random.default_rng(1000 + q)
        row = {"n": n, "fold": t["fold"], "nw_full": int(nw), "k1": {}, "k2": {}, "k3": {}}
        for s in S1 + [nw]:
            if s > nw:
                continue
            row["k1"][str(s)] = float(np.mean([rr[rng.choice(nw, s, replace=False)].min()
                                               for _ in range(R)]))
        for s in SK:
            v2, v3 = [], []
            for _ in range(R):
                cost = np.full((n, n), np.inf)
                for L in range(LMIN, n - LMIN + 1):
                    b = A.bank(t["fold"], L)
                    sel = rng.choice(len(b["Wc"]), min(s, len(b["Wc"])), replace=False)
                    segs = np.stack([nat[a:a + L] for a in range(n - L + 1)])
                    M = A.msd_many(b["Wc"][sel], b["n2"][sel], segs) * L
                    for a in range(n - L + 1):
                        cost[a, L - 1] = M[:, a].min()
                d = dp_ladder(cost, n, kmax=3, lmin=LMIN)
                v2.append(d[2]); v3.append(d[3] if d[3] is not None else d[2])
            row["k2"][str(s)] = float(np.mean(v2))
            row["k3"][str(s)] = float(np.mean(v3))
        # effective candidate counts
        row["c2"] = int(max(1, n - 2 * LMIN + 1))
        row["c3"] = int(sum(1 for c1 in range(LMIN, n - 2 * LMIN + 1)
                            for c2 in range(c1 + LMIN, n - LMIN + 1)))
        res[t["pdb"]] = row
        if (q + 1) % 10 == 0:
            print(f"{q+1}/{len(tg)} {time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)
    I.write("asm_scaling", res)

    print("\n--- mean oracle floor vs candidate count (126 targets, rigid placement) ---")
    print(f"{'family':6s} {'s':>7s} {'log10 eff':>10s} {'mean RMSD':>10s}")
    pts = []
    for s in S1:
        v = [res[t["pdb"]]["k1"][str(s)] for t in tg if str(s) in res[t["pdb"]]["k1"]]
        if len(v) == len(tg):
            eff = np.mean([np.log10(s) for t in tg])
            print(f"{'k=1':6s} {s:7d} {eff:10.2f} {np.mean(v):10.3f}")
            pts.append(("k1", eff, np.mean(v)))
    v = [res[t["pdb"]]["k1"][str(res[t["pdb"]]["nw_full"])] for t in tg]
    eff = np.mean([np.log10(res[t["pdb"]]["nw_full"]) for t in tg])
    print(f"{'k=1':6s} {'full':>7s} {eff:10.2f} {np.mean(v):10.3f}")
    pts.append(("k1", eff, np.mean(v)))
    for fam, p in (("k=2", 2), ("k=3", 3)):
        for s in SK:
            v = [res[t["pdb"]][f"k{p}"][str(s)] for t in tg]
            eff = np.mean([np.log10(res[t["pdb"]][f"c{p}"]) + p * np.log10(s) for t in tg])
            print(f"{fam:6s} {s:7d} {eff:10.2f} {np.mean(v):10.3f}")
            pts.append((fam, eff, np.mean(v)))
    return res


if __name__ == "__main__":
    main()
