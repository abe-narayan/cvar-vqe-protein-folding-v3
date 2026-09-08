"""ASM-4b: cross the coverage table with native secondary structure and chain position,
and with the shape-dimension of an m-mer CA trace (2m-5 free shape coordinates once the
3.8 A virtual bonds are fixed).  Reads `s12/results/asm_coverage.json`.
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I                        # noqa: E402
import peptide_db as db                                # noqa: E402

MS = (4, 5, 6, 7, 8)


def main():
    cov = json.load(open(os.path.join(ROOT, "s12", "results", "asm_coverage.json")))
    tg = I.targets()
    ss = {}
    for t in tg:
        p = db.by_pdb(t["pdb"])
        ss[t["pdb"]] = I.ss_of(p.phi, p.psi) if p is not None and len(p.seq) == t["n"] else None
    n_ss = sum(1 for v in ss.values() if v)
    out = {"n_with_ss": n_ss, "by_m": {}, "ss_m6": {}, "pos_m6": {}, "rigid_link": {}}
    print(f"native SS available for {n_ss}/126 targets")

    dof = {m: 2 * m - 5 for m in MS}
    print("\nm   shape-DOF  bank-size  mean-of-worst-hole  mean-hole  frac>0.5")
    for m in MS:
        w, a = [], []
        for t in tg:
            c = cov[t["pdb"]]["cov"].get(str(m))
            if c:
                w.append(max(c)); a.extend(c)
        out["by_m"][m] = {"worst": float(np.mean(w)), "mean": float(np.mean(a)),
                          "frac_gt_half": float(np.mean(np.array(a) > 0.5)), "dof": dof[m]}
        print(f"{m}   {dof[m]:8d}   {'~40-80k':>9s}  {np.mean(w):17.3f}  {np.mean(a):9.3f}"
              f"  {np.mean(np.array(a) > 0.5):8.3f}")

    for m in (6, 8):
        acc, pos = {}, {}
        for t in tg:
            s = ss[t["pdb"]]; c = cov[t["pdb"]]["cov"].get(str(m))
            if not s or not c:
                continue
            for a, v in enumerate(c):
                seg = s[a:a + m]
                lab = "H" if seg.count("H") >= m - 1 else ("E" if seg.count("E") >= m - 1
                                                           else ("mixed" if len(set(seg)) > 1 else "C"))
                acc.setdefault(lab, []).append(v)
                p = "Nterm" if a == 0 else ("Cterm" if a + m == len(s) else "mid")
                pos.setdefault(p, []).append(v)
        print(f"\n--- coverage of native {m}-mers by SS class ---")
        for k, v in sorted(acc.items()):
            print(f"  {k:6s} n={len(v):5d} mean {np.mean(v):.3f}  p90 {np.percentile(v,90):.3f}"
                  f"  max {np.max(v):.3f}")
        print(f"--- by chain position ---")
        for k, v in sorted(pos.items()):
            print(f"  {k:6s} n={len(v):5d} mean {np.mean(v):.3f}  p90 {np.percentile(v,90):.3f}")
        out["ss_m6" if m == 6 else "ss_m8"] = {k: [len(v), float(np.mean(v))] for k, v in acc.items()}
        out["pos_m6" if m == 6 else "pos_m8"] = {k: [len(v), float(np.mean(v))] for k, v in pos.items()}

    # does the local hole predict the whole-window floor?
    rl = json.load(open(os.path.join(ROOT, "s12", "results", "asm_rigid_ladder.json")))
    k1 = rl["arms"]["full_lmin4"]["1"]["per_target"]
    for m in MS:
        x = np.array([max(cov[t["pdb"]]["cov"][str(m)]) for t in tg])
        y = np.array([k1[t["pdb"]] for t in tg])
        r = float(np.corrcoef(x, y)[0, 1])
        out["rigid_link"][m] = r
        print(f"corr(worst {m}-mer hole, whole-window floor) = {r:+.3f}")
    I.write("asm_coverage_ss", out)
    return out


if __name__ == "__main__":
    main()
