"""ASM-3: matched-budget CHAIN-CONSISTENT assembly floors, with the capacity null.

Every arm uses the SAME search procedure: per cut, take the top-M pieces of each side by
their own ideal-rebuild SSD to the corresponding native segment (ORACLE selection), then
exhaustively evaluate the M x M torsion concatenations with `build_ca_exact`.  Arms differ
only in the CANDIDATE BANK:

    real     -- every length-L window of the fold's legal library
    rama     -- CAPACITY NULL: a synthetic bank of the SAME SIZE whose torsions are drawn
                i.i.d. from the pooled (phi,psi) of the same library.  Same number of
                degrees of freedom, same search, zero real fragment structure.
    blosN    -- real bank restricted to the top-N by BLOSUM62 vs the target's sub-sequence
    randN    -- real bank restricted to N uniformly random windows (retrieval null)

`--k3` adds the 3-piece ladder (top-M3 per slot).

Usage: python -m s12.asm_chain2 [--arms real,rama,blos500,rand500] [--M 300] [--limit N]
"""
from __future__ import annotations
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402

LMIN = 4
_SYNTH = {}
_POOL = {}


def rama_pool(fold):
    """All (phi, psi) residue pairs of the fold's legal library."""
    if fold in _POOL:
        return _POOL[fold]
    b = A.bank(fold, 3)
    P = np.stack([b["PHI"].ravel(), b["PSI"].ravel()], 1)
    _POOL[fold] = P
    return P


def synth_bank(fold, L, nw, seed=7):
    """Capacity null bank: nw pieces of length L, torsions i.i.d. from `rama_pool`."""
    key = (fold, L, nw, seed)
    if key in _SYNTH:
        return _SYNTH[key]
    P = rama_pool(fold)
    rng = np.random.default_rng(seed * 1000003 + fold * 101 + L)
    idx = rng.integers(0, len(P), size=(nw, L))
    PHI = P[idx, 0]; PSI = P[idx, 1]
    parts = [A.build_ca(PHI[a:a + 20000], PSI[a:a + 20000]) for a in range(0, nw, 20000)]
    Ideal = np.concatenate(parts, 0)
    d = {"PHI": PHI, "PSI": PSI, "Ic": Ideal - Ideal.mean(1, keepdims=True)}
    d["In2"] = (d["Ic"] ** 2).sum((1, 2))
    _SYNTH[key] = d
    return d


def _candidates(arm, fold, L, seg, codes_sub, rng):
    """(PHI, PSI, ssd) for the arm's candidate bank on this native segment."""
    if arm == "rama":
        nwr = A.bank_ideal(fold, L, lean=True)["nw"]
        b = synth_bank(fold, L, nwr)
        ssd = A.msd_many(b["Ic"], b["In2"], seg[None])[:, 0] * L
        return b["PHI"], b["PSI"], ssd, np.arange(len(ssd))
    b = A.bank_ideal(fold, L, lean=True)
    ssd = A.msd_many(b["Ic"], b["In2"], seg[None])[:, 0] * L
    sel = np.arange(len(ssd))
    if arm.startswith("blos"):
        N = int(arm[4:])
        sim = A.B62()[b["S"], codes_sub[None, :]].sum(1)
        sel = np.argsort(-sim, kind="stable")[:N]
    elif arm.startswith("rand"):
        N = int(arm[4:])
        sel = rng.choice(len(ssd), size=min(N, len(ssd)), replace=False)
    return b["PHI"], b["PSI"], ssd[sel], sel


def two_piece(t, nat, arm, M=300, rng=None, chunk=60000):
    n, fold = t["n"], t["fold"]
    codes = A.encode(t["seq"])
    best = (np.inf, None)
    for c in range(LMIN, n - LMIN + 1):
        PH1, PS1, s1, sel1 = _candidates(arm, fold, c, nat[:c], codes[:c], rng)
        PH2, PS2, s2, sel2 = _candidates(arm, fold, n - c, nat[c:], codes[c:], rng)
        k1 = sel1[np.argsort(s1)[:M]]; k2 = sel2[np.argsort(s2)[:M]]
        ii, jj = np.meshgrid(np.arange(len(k1)), np.arange(len(k2)), indexing="ij")
        ii = ii.ravel(); jj = jj.ravel()
        for a in range(0, len(ii), chunk):
            p1 = k1[ii[a:a + chunk]]; p2 = k2[jj[a:a + chunk]]
            phi = np.concatenate([PH1[p1], PH2[p2]], 1)
            psi = np.concatenate([PS1[p1], PS2[p2]], 1)
            r = A.rmsd_to(A.build_ca(phi, psi), nat)
            k = int(np.argmin(r))
            if r[k] < best[0]:
                best = (float(r[k]), c)
    return best


def three_piece(t, nat, arm, M=40, rng=None, chunk=60000):
    n, fold = t["n"], t["fold"]
    codes = A.encode(t["seq"])
    best = (np.inf, None)
    for c1 in range(LMIN, n - 2 * LMIN + 1):
        for c2 in range(c1 + LMIN, n - LMIN + 1):
            segs = [(0, c1), (c1, c2), (c2, n)]
            P, Q, K = [], [], []
            for (a, b) in segs:
                PH, PS, s, sel = _candidates(arm, fold, b - a, nat[a:b], codes[a:b], rng)
                k = sel[np.argsort(s)[:M]]
                P.append(PH); Q.append(PS); K.append(k)
            g = np.meshgrid(*[np.arange(len(k)) for k in K], indexing="ij")
            g = [x.ravel() for x in g]
            for a0 in range(0, len(g[0]), chunk):
                phi = np.concatenate([P[m][K[m][g[m][a0:a0 + chunk]]] for m in range(3)], 1)
                psi = np.concatenate([Q[m][K[m][g[m][a0:a0 + chunk]]] for m in range(3)], 1)
                r = A.rmsd_to(A.build_ca(phi, psi), nat)
                k = int(np.argmin(r))
                if r[k] < best[0]:
                    best = (float(r[k]), (c1, c2))
    return best


def one_piece(t, nat, arm, rng=None):
    n, fold = t["n"], t["fold"]
    PH, PS, s, sel = _candidates(arm, fold, n, nat, A.encode(t["seq"]), rng)
    k = int(np.argmin(s))
    return float(np.sqrt(s[k] / n))


def main():
    arms = "real,rama,blos500,rand500"
    M, limit, do_k3 = 300, None, False
    av = sys.argv
    if "--arms" in av: arms = av[av.index("--arms") + 1]
    if "--M" in av: M = int(av[av.index("--M") + 1])
    if "--limit" in av: limit = int(av[av.index("--limit") + 1])
    if "--k3" in av: do_k3 = True
    arms = arms.split(",")
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    tg = sorted(tg, key=lambda t: (t["fold"], t["pdb"]))     # bank locality
    nat = {}
    for t in tg:
        u = I.load_univ(t["pdb"]); nat[t["pdb"]] = u["nat_ca"]; del u
    res = {}
    t0 = time.time(); cur = None
    for q, t in enumerate(tg):
        if t["fold"] != cur:
            A.drop_bank(); _SYNTH.clear(); _POOL.clear(); cur = t["fold"]
        rng = np.random.default_rng(20260905 + q)
        row = {"n": t["n"], "fold": t["fold"]}
        for arm in arms:
            row[f"k1_{arm}"] = one_piece(t, nat[t["pdb"]], arm, rng)
            r2, c = two_piece(t, nat[t["pdb"]], arm, M=M, rng=rng)
            row[f"k2_{arm}"] = r2; row[f"cut_{arm}"] = c
            if do_k3 and t["n"] >= 3 * LMIN:
                row[f"k3_{arm}"] = three_piece(t, nat[t["pdb"]], arm, M=40, rng=rng)[0]
        res[t["pdb"]] = row
        print(f"{q+1:3d}/{len(tg)} {t['pdb']} n={t['n']} " +
              " ".join(f"{a}:{row['k1_'+a]:.2f}/{row['k2_'+a]:.2f}" +
                       (f"/{row.get('k3_'+a, float('nan')):.2f}" if do_k3 else "") for a in arms) +
              f"  {time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)
    tag = "asm_chain2_M%d" % M + ("_k3" if do_k3 else "")
    I.write(tag, res)
    f18 = np.array([t["pdb"] in I.FAIL18 for t in tg])
    print("\narm        k1     k2   " + ("  k3" if do_k3 else ""))
    for arm in arms:
        a1 = np.array([res[t["pdb"]][f"k1_{arm}"] for t in tg])
        a2 = np.array([res[t["pdb"]][f"k2_{arm}"] for t in tg])
        s = f"{arm:9s} {a1.mean():.3f} {a2.mean():.3f}"
        if do_k3:
            a3 = np.array([res[t["pdb"]].get(f"k3_{arm}", np.nan) for t in tg])
            s += f" {np.nanmean(a3):.3f}"
        print(s + f"   | FAIL18 k2 {a2[f18].mean():.3f}  other {a2[~f18].mean():.3f}")
    return res


if __name__ == "__main__":
    main()
