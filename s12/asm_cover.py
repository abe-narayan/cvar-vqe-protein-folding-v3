"""ASM-4: CONFORMATIONAL COVERAGE of the legal library -- is the deficiency local
(the library does not contain the native's local motifs) or global (it contains every
local motif but no whole window is close)?

For every native sub-segment of length m in {4,5,6,7,8} we take the minimum CA-RMSD over
EVERY length-m window of the fold's legal library (ORACLE / DIAGNOSTIC), and cross it with
the segment's native secondary structure (`I.ss_of` on the native torsions) and its
position (N-terminal / interior / C-terminal).

Usage: python -m s12.asm_cover
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402

MS = (4, 5, 6, 7, 8)


def native_ss(t, nat):
    import peptide_db as db
    p = db.by_pdb(t["pdb"])
    if p is None or len(p.seq) != t["n"]:
        return None
    ss = I.ss_of(p.phi[None, :], p.psi[None, :])
    ss = ss[0] if isinstance(ss, (list, tuple, np.ndarray)) else ss
    return str(ss)


def main():
    tg = I.targets()
    tg = sorted(tg, key=lambda t: (t["fold"], t["pdb"]))
    nat = {t["pdb"]: I.load_univ(t["pdb"])["nat_ca"] for t in tg}
    ss = {}
    for t in tg:
        try:
            ss[t["pdb"]] = native_ss(t, nat[t["pdb"]])
        except Exception:
            ss[t["pdb"]] = None
    out = {}
    cur = None; t0 = time.time()
    for t in tg:
        if t["fold"] != cur:
            A.drop_bank(); cur = t["fold"]
        n = t["n"]; row = {}
        for m in MS:
            if n < m:
                continue
            b = A.bank(t["fold"], m)
            segs = np.stack([nat[t["pdb"]][a:a + m] for a in range(n - m + 1)])
            R = A.rmsd_many(b["Wc"], b["n2"], segs)
            row[str(m)] = R.min(0).tolist()
        out[t["pdb"]] = {"n": n, "fold": t["fold"], "ss": ss[t["pdb"]], "cov": row}
        print(f"{t['pdb']} n={n} " + " ".join(f"m{m}:{max(row[str(m)]):.2f}" for m in MS if str(m) in row)
              + f"  {time.time()-t0:.0f}s", flush=True)
    I.write("asm_coverage", out)

    # ---------------- summaries
    f18 = set(I.FAIL18)
    print("\n--- worst local hole per target (max over segments), by m ---")
    for m in MS:
        w = [max(out[t["pdb"]]["cov"][str(m)]) for t in tg if str(m) in out[t["pdb"]]["cov"]]
        a = [np.mean(out[t["pdb"]]["cov"][str(m)]) for t in tg if str(m) in out[t["pdb"]]["cov"]]
        wf = [max(out[t["pdb"]]["cov"][str(m)]) for t in tg
              if str(m) in out[t["pdb"]]["cov"] and t["pdb"] in f18]
        print(f"m={m}  mean-of-worst {np.mean(w):.3f}  mean-of-mean {np.mean(a):.3f}  "
              f"max {np.max(w):.3f}  frac targets with a hole >0.5A "
              f"{np.mean(np.array(w) > 0.5):.2f}   FAIL18 mean-of-worst {np.mean(wf):.3f}")

    print("\n--- coverage by native SS of the segment (m=6) ---")
    acc = {}
    for t in tg:
        s = out[t["pdb"]]["ss"]; c = out[t["pdb"]]["cov"].get("6")
        if not s or not c:
            continue
        for a, v in enumerate(c):
            lab = max("HEC", key=lambda ch: s[a:a + 6].count(ch))
            acc.setdefault(lab, []).append(v)
            pos = "Nterm" if a == 0 else ("Cterm" if a + 6 == len(s) else "mid")
            acc.setdefault(pos, []).append(v)
    for k, v in acc.items():
        print(f"  {k:6s} n={len(v):5d}  mean {np.mean(v):.3f}  p90 {np.percentile(v, 90):.3f}  "
              f"frac>0.5 {np.mean(np.array(v) > 0.5):.3f}")
    return out


if __name__ == "__main__":
    main()
