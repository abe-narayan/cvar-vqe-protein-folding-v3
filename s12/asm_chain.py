"""ASM-2: the CHAIN-CONSISTENT oracle assembly ladder.

A chain-consistent assembly is ONE ideal-geometry chain: residue i's (phi_i, psi_i) is
taken from whichever piece covers i, the torsion vector is concatenated, and the CA trace
is rebuilt with `core.project.build_ca_exact` (= `I.build_ca`).  There are ZERO free
placement parameters -- the pieces' relative pose is dictated by the junction torsions.

Exactness.  For a 2-piece assembly cut at c, CA_0..CA_{c-1} of the rebuilt chain are a
rigid copy of piece 1's OWN ideal rebuild and CA_c..CA_{n-1} of piece 2's, so

    n * RMSD(assembly, native)^2  >=  SSD_ideal(p1, T[:c]) + SSD_ideal(p2, T[c:])

which is an admissible lower bound.  We branch and bound on it: the reported floor is the
EXACT minimum over the full library x full library product whenever `capped` is false.

Usage: python -m s12.asm_chain [--targets N]
"""
from __future__ import annotations
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402

LMIN = 4
PAIRCAP = 1_500_000
BUILD_CHUNK = 60_000


def _ideal_ssd(fold, L, T_seg):
    b = A.bank_ideal(fold, L)
    return A.msd_many(b["Ic"], b["In2"], T_seg[None])[:, 0] * L, b


def two_piece(t, nat, verbose=False):
    """Exact (or capped) chain-consistent 2-piece floor, over all cuts."""
    n, fold = t["n"], t["fold"]
    best = {"rmsd": np.inf, "cut": None, "i": None, "j": None, "capped": False, "pairs": 0}
    for c in range(LMIN, n - LMIN + 1):
        s1, b1 = _ideal_ssd(fold, c, nat[:c])
        s2, b2 = _ideal_ssd(fold, n - c, nat[c:])
        o1 = np.argsort(s1); o2 = np.argsort(s2)
        S1 = s1[o1]; S2 = s2[o2]
        if S1[0] + S2[0] >= best["rmsd"] ** 2 * n:
            continue
        # seed the incumbent with a small grid
        inc = best["rmsd"] ** 2 * n
        seed = 48
        for pass_ in range(2):
            if pass_ == 0:
                ii, jj = np.meshgrid(np.arange(min(seed, len(S1))),
                                     np.arange(min(seed, len(S2))), indexing="ij")
                ii = ii.ravel(); jj = jj.ravel()
                keep = S1[ii] + S2[jj] < inc
                ii, jj = ii[keep], jj[keep]
            else:
                lim = inc - S1
                cnt = np.searchsorted(S2, lim)
                tot = int(cnt.sum())
                best["pairs"] += tot
                if tot > PAIRCAP:
                    best["capped"] = True
                    ordr = np.argsort(-cnt)          # keep the widest rows first
                    take = np.cumsum(cnt[ordr]) <= PAIRCAP
                    rows = ordr[take]
                else:
                    rows = np.where(cnt > 0)[0]
                if len(rows) == 0:
                    break
                ii = np.repeat(rows, cnt[rows])
                jj = np.concatenate([np.arange(cnt[r]) for r in rows]) if len(rows) else np.array([], int)
            if len(ii) == 0:
                break
            for a in range(0, len(ii), BUILD_CHUNK):
                p1 = o1[ii[a:a + BUILD_CHUNK]]; p2 = o2[jj[a:a + BUILD_CHUNK]]
                phi = np.concatenate([b1["PHI"][p1], b2["PHI"][p2]], 1)
                psi = np.concatenate([b1["PSI"][p1], b2["PSI"][p2]], 1)
                r = A.rmsd_to(A.build_ca(phi, psi), nat)
                k = int(np.argmin(r))
                if r[k] < best["rmsd"]:
                    best.update(rmsd=float(r[k]), cut=c, i=int(p1[k]), j=int(p2[k]))
            inc = best["rmsd"] ** 2 * n
    return best


def one_piece(t, nat):
    b = A.bank_ideal(t["fold"], t["n"])
    r = A.rmsd_many(b["Ic"], b["In2"], nat[None])[:, 0]
    k = int(np.argmin(r))
    return {"rmsd": float(r[k]), "i": k}


def main(limit=None):
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    nat = {}
    for t in tg:
        u = I.load_univ(t["pdb"]); nat[t["pdb"]] = u["nat_ca"]; del u
    res = {}
    t0 = time.time()
    cur_fold = None
    for k, t in enumerate(tg):
        if t["fold"] != cur_fold:
            A.drop_bank(); cur_fold = t["fold"]
        o = one_piece(t, nat[t["pdb"]])
        two = two_piece(t, nat[t["pdb"]])
        res[t["pdb"]] = {"n": t["n"], "fold": t["fold"], "k1_chain": o["rmsd"],
                         "k2_chain": two["rmsd"], "cut": two["cut"],
                         "capped": bool(two["capped"]), "pairs": int(two["pairs"])}
        print(f"{k+1:3d}/{len(tg)} {t['pdb']} n={t['n']} k1={o['rmsd']:.3f} "
              f"k2={two['rmsd']:.3f} cut={two['cut']} pairs={two['pairs']:>8d}"
              f"{' CAP' if two['capped'] else ''} {time.time()-t0:.0f}s "
              f"free {I.free_gb():.2f}", flush=True)
    I.write("asm_chain_ladder", res)
    a = np.array([res[t["pdb"]]["k1_chain"] for t in tg])
    b = np.array([res[t["pdb"]]["k2_chain"] for t in tg])
    f18 = np.array([t["pdb"] in I.FAIL18 for t in tg])
    print(f"\nk=1 chain {a.mean():.3f}  (FAIL18 {a[f18].mean():.3f} / other {a[~f18].mean():.3f})")
    print(f"k=2 chain {b.mean():.3f}  (FAIL18 {b[f18].mean():.3f} / other {b[~f18].mean():.3f})")
    print(f"capped on {sum(res[p]['capped'] for p in res)} targets")
    return res


if __name__ == "__main__":
    lim = None
    if "--targets" in sys.argv:
        lim = int(sys.argv[sys.argv.index("--targets") + 1])
    main(lim)
