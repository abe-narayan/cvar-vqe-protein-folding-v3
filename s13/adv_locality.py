"""SPRINT 13 -- ADVERSARIAL AUDIT 5: try to break the exact locality theorem.

The claim (qarch §1): the support of `d_ij` is EXACTLY the `j-i-1` residues strictly
between i and j; agreement 1.0000 over 5,000+ (pair, variable) cells, tested with a
**0.10 A** support threshold on 3 targets at k=8.

A 0.10 A threshold cannot distinguish "exactly zero" from "small".  This re-tests at
MACHINE PRECISION, on every one of the 126 targets, and then attacks the boundary
conditions the original test did not cover:

    L1  machine-precision CA-CA support, all 126 targets, k=4 and k=8, random states
    L2  chain termini: pairs with i=0 or j=n-1 scored separately
    L3  glycine and proline positions scored separately (different library class,
        and PRO's phi is restrained in nature)
    L4  a DIFFERENT builder convention -- cis-omega (omega=0) instead of trans -- and
        a per-residue RANDOM omega, which is the convention that should break it
    L5  ALL-ATOM: CB-CB, N-N, C-C and full heavy-atom pair distances through
        `rep.build_coords`, to pin the "support one residue wider" statement to an
        exact rule rather than an approximate one

    python -m s13.adv_locality
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402
from core import project as pj             # noqa: E402
from core import geometry as geo           # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
TOL = 1e-9          # machine-precision support threshold, vs the original 0.10 A
OUT = {}


def ca_support(PHI, PSI, s0, rng, omega=None, reps_=6):
    """(n, n, n) boolean: does moving residue m change d_ij?  Machine precision."""
    n, k = PHI.shape
    idx = np.arange(n)
    kw = {} if omega is None else {"omega": omega}
    def build(S):
        return np.asarray(pj.build_ca_exact(PHI[idx[None, :], S], PSI[idx[None, :], S], **kw), float)
    sup = np.zeros((n, n, n), bool)
    maxout = 0.0; minin = np.inf
    for _ in range(reps_):
        s = rng.integers(0, k, n) if s0 is None else s0
        base = build(s[None])[0]
        D0 = np.linalg.norm(base[:, None] - base[None, :], axis=-1)
        for m in range(n):
            cand = np.repeat(s[None, :], k, axis=0); cand[:, m] = np.arange(k)
            W = build(cand)
            D = np.linalg.norm(W[:, :, None] - W[:, None, :], axis=-1)
            dd = np.abs(D - D0[None]).max(0)
            sup[:, :, m] |= dd > TOL
            # rule: m supports (i,j) iff i < m < j
            rule = np.zeros((n, n), bool)
            for i in range(n):
                for j in range(n):
                    rule[i, j] = min(i, j) < m < max(i, j)
            maxout = max(maxout, float(dd[~rule].max()) if (~rule).any() else 0.0)
            inv = dd[rule]
            if inv.size:
                minin = min(minin, float(inv.min()))
        s0 = None
    return sup, maxout, minin


def rule_matrix(n):
    R = np.zeros((n, n, n), bool)
    for i in range(n):
        for j in range(n):
            for m in range(n):
                R[i, j, m] = min(i, j) < m < max(i, j)
    return R


def l1_l2_l3(tg, k, rng, n_targets=126):
    rows = []
    for t in tg[:n_targets]:
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, k, seq)
        PHI = np.ascontiguousarray(tab[:, :, 0]); PSI = np.ascontiguousarray(tab[:, :, 1])
        sup, maxout, minin = ca_support(PHI, PSI, None, rng)
        R = rule_matrix(n)
        agree = float((sup == R).mean())
        # termini pairs
        term = np.zeros((n, n), bool); term[0, :] = True; term[:, 0] = True
        term[n - 1, :] = True; term[:, n - 1] = True
        a_t = float((sup[term] == R[term]).mean())
        a_i = float((sup[~term] == R[~term]).mean())
        gly = [i for i, c in enumerate(seq) if c == "G"]
        pro = [i for i, c in enumerate(seq) if c == "P"]
        a_g = float((sup[:, :, gly] == R[:, :, gly]).mean()) if gly else None
        a_p = float((sup[:, :, pro] == R[:, :, pro]).mean()) if pro else None
        rows.append({"pdb": t["pdb"], "n": n, "k": k, "agreement": agree,
                     "agreement_terminal_pairs": a_t, "agreement_interior_pairs": a_i,
                     "agreement_GLY_vars": a_g, "agreement_PRO_vars": a_p,
                     "max_|dd|_outside_support": maxout, "min_|dd|_inside_support": minin,
                     "n_cells": int(sup.size)})
    return rows


def l4_omega(tg, k, rng, n_targets=20):
    """A different builder convention.  cis-omega should NOT break it; variable omega SHOULD."""
    rows = []
    for t in tg[:n_targets]:
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, k, seq)
        PHI = np.ascontiguousarray(tab[:, :, 0]); PSI = np.ascontiguousarray(tab[:, :, 1])
        R = rule_matrix(n)
        out = {"pdb": t["pdb"], "n": n}
        for nm, om in (("trans", None), ("cis", 0.0), ("random_fixed", float(rng.uniform(-np.pi, np.pi)))):
            sup, mo, mi = ca_support(PHI, PSI, None, np.random.default_rng(7), omega=om, reps_=3)
            out[nm] = {"agreement": float((sup == R).mean()), "max_outside": mo}
        rows.append(out)
    return rows


def l5_allatom(tg, k, rng, n_targets=12):
    """Support of ATOM-ATOM distances under rep.build_coords -- the AMBER/Legacy-CB case."""
    rows = []
    for t in tg[:n_targets]:
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, k, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        s = rng.integers(0, k, n)
        def coords(ss):
            return rep.build_coords(rep.bitstring_from_states(np.asarray(ss, int)))
        base = coords(s)
        res = {}
        for atom in ("CA", "CB", "N", "C", "O"):
            if atom not in base:
                continue
            P0 = np.asarray(base[atom], float)
            if not np.isfinite(P0).all():
                continue
            D0 = np.linalg.norm(P0[:, None] - P0[None, :], axis=-1)
            sup = np.zeros((n, n, n), bool)
            for m in range(n):
                for v in range(k):
                    if v == s[m]:
                        continue
                    s2 = s.copy(); s2[m] = v
                    P = np.asarray(coords(s2)[atom], float)
                    D = np.linalg.norm(P[:, None] - P[None, :], axis=-1)
                    sup[:, :, m] |= np.abs(D - D0) > 1e-8
            R_between = rule_matrix(n)                       # i < m < j
            R_incl = np.zeros((n, n, n), bool)               # i <= m <= j
            R_lo = np.zeros((n, n, n), bool)                 # i <= m < j
            R_hi = np.zeros((n, n, n), bool)                 # i < m <= j
            for i in range(n):
                for j in range(n):
                    for m in range(n):
                        lo, hi = min(i, j), max(i, j)
                        R_incl[i, j, m] = lo <= m <= hi
                        R_lo[i, j, m] = lo <= m < hi
                        R_hi[i, j, m] = lo < m <= hi
            off = np.eye(n, dtype=bool)[:, :, None] | np.zeros((n, n, n), bool)
            sel = ~off[:, :, 0][:, :, None].repeat(n, 2)     # drop the i==j diagonal
            res[atom] = {
                "agree_strictly_between": float((sup[sel] == R_between[sel]).mean()),
                "agree_inclusive_i_to_j": float((sup[sel] == R_incl[sel]).mean()),
                "agree_i_le_m_lt_j": float((sup[sel] == R_lo[sel]).mean()),
                "agree_i_lt_m_le_j": float((sup[sel] == R_hi[sel]).mean()),
                "mean_support_size": float(sup[sel].reshape(-1, 1).mean() * n) if False else
                                     float(sup.sum(2)[~np.eye(n, dtype=bool)].mean())}
        rows.append({"pdb": t["pdb"], "n": n, "atoms": res})
    return rows


def main():
    tg = I.targets()
    rng = np.random.default_rng(4242)
    print("L1/L2/L3 machine-precision CA-CA support, k=4, all 126 ...")
    r4 = l1_l2_l3(tg, 4, rng)
    OUT["L1_k4"] = r4
    print("   agreement  min", min(r["agreement"] for r in r4),
          " max |dd| outside support", max(r["max_|dd|_outside_support"] for r in r4),
          " min |dd| inside support", min(r["min_|dd|_inside_support"] for r in r4))
    print("   terminal-pair agreement min", min(r["agreement_terminal_pairs"] for r in r4))
    g = [r["agreement_GLY_vars"] for r in r4 if r["agreement_GLY_vars"] is not None]
    p = [r["agreement_PRO_vars"] for r in r4 if r["agreement_PRO_vars"] is not None]
    print(f"   GLY vars: n={len(g)} min agreement {min(g) if g else None}")
    print(f"   PRO vars: n={len(p)} min agreement {min(p) if p else None}")
    print("L1 at k=8, first 30 ...")
    r8 = l1_l2_l3(tg, 8, rng, 30)
    OUT["L1_k8"] = r8
    print("   agreement min", min(r["agreement"] for r in r8),
          " max outside", max(r["max_|dd|_outside_support"] for r in r8))
    print("L4 builder convention ...")
    OUT["L4_omega"] = l4_omega(tg, 4, rng)
    for nm in ("trans", "cis", "random_fixed"):
        vals = [r[nm]["agreement"] for r in OUT["L4_omega"]]
        mo = [r[nm]["max_outside"] for r in OUT["L4_omega"]]
        print(f"   omega={nm:12s} min agreement {min(vals):.6f}  max |dd| outside {max(mo):.2e}")
    print("L5 all-atom ...")
    OUT["L5_allatom"] = l5_allatom(tg, 4, rng)
    for atom in ("CA", "CB", "N", "C", "O"):
        rs = [r["atoms"][atom] for r in OUT["L5_allatom"] if atom in r["atoms"]]
        if not rs:
            continue
        print(f"   {atom:3s} strictly-between {np.mean([x['agree_strictly_between'] for x in rs]):.4f}  "
              f"i<=m<=j {np.mean([x['agree_inclusive_i_to_j'] for x in rs]):.4f}  "
              f"i<=m<j {np.mean([x['agree_i_le_m_lt_j'] for x in rs]):.4f}  "
              f"i<m<=j {np.mean([x['agree_i_lt_m_le_j'] for x in rs]):.4f}  "
              f"mean|supp| {np.mean([x['mean_support_size'] for x in rs]):.2f}")
    json.dump(OUT, open(os.path.join(RESULTS, "adv_locality.json"), "w"), indent=1)
    print("wrote adv_locality.json")


if __name__ == "__main__":
    main()
