"""s26/a_dla_check.py -- lane A (Adversary): independent re-derivation of lane Q's A2 (L27).

Written from the construction in s26/PREREG_A2.md and core/quantum.py, without reading
s26/q_dla.py's closure code into this file: Pauli strings as (x, z) bit pairs, CNOT conjugation
X_c -> X_c X_t and Z_t -> Z_c Z_t, the entangler C = CNOT chain (q -> q+1) plus ring (n-1 -> 0),
generators {C^k Y_q C^-k}, breadth-first Lie closure under "commute with a generator".  Both
conjugation conventions and both gate orders are computed.  The numeric cross-check at n = 4, 5
uses a DIFFERENT rank rule from lane Q's (numpy.linalg.matrix_rank on the Gram matrix of the
vectorised real antisymmetric matrices, default tolerance) so the two implementations share no
code and no tolerance.  Property measurement: no native, no RMSD.  Under one minute, under
200 MB, so run directly (contract section 3).

Writes s26/results/a_dla_check.json.
"""
from __future__ import annotations

import itertools
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s26", "results", "a_dla_check.json")


def dim_so(n):
    return (1 << (n - 1)) * ((1 << n) - 1)


def cnot(xz, c, t, n):
    x, z = xz
    bc, bt = 1 << (n - 1 - c), 1 << (n - 1 - t)
    if x & bc:
        x ^= bt
    if z & bt:
        z ^= bc
    return (x, z)


def entangler(n, ring=True):
    pairs = [(q, q + 1) for q in range(n - 1)]
    if ring:
        pairs.append((n - 1, 0))
    return pairs


def conj(xz, n, power, gate_order="forward"):
    pairs = entangler(n)
    if gate_order == "reverse":
        pairs = pairs[::-1]
    for _ in range(power):
        for c, t in pairs:
            xz = cnot(xz, c, t, n)
    return xz


def generators(n, L, convention="right", gate_order="forward"):
    gens = []
    for q in range(n):
        y = (1 << (n - 1 - q), 1 << (n - 1 - q))            # Y_q: x and z bit set
        ks = range(1, L + 1) if convention == "right" else range(0, L)
        for k in ks:
            gens.append(conj(y, n, k, gate_order if convention == "right" else
                             ("reverse" if gate_order == "forward" else "forward")))
    return gens


def anticommute(a, b):
    return (bin(a[0] & b[1]).count("1") + bin(a[1] & b[0]).count("1")) % 2 == 1


def closure(gens, cap=1 << 20):
    seen = set(gens)
    frontier = list(seen)
    while frontier:
        new = []
        for s in frontier:
            for g in gens:
                if anticommute(s, g):
                    p = (s[0] ^ g[0], s[1] ^ g[1])
                    if p not in seen:
                        seen.add(p)
                        new.append(p)
                        if len(seen) > cap:
                            return len(seen), True, seen
        frontier = new
    return len(seen), False, seen


def odd_y(s):
    return bin(s[0] & s[1]).count("1") % 2 == 1


# ---- independent numeric route: real antisymmetric matrices, Gram-matrix rank
_I = np.eye(2); _X = np.array([[0, 1], [1, 0.0]]); _Z = np.array([[1, 0], [0, -1.0]])
_iY = np.array([[0, 1], [-1, 0.0]])       # -i * Y = [[0,-1],[1,0]] real; sign irrelevant to span


def dense_real(s, n):
    """The real matrix proportional to -i P for an odd-Y string P (real antisymmetric)."""
    M = np.eye(1)
    for q in range(n):
        bit = 1 << (n - 1 - q)
        xb, zb = bool(s[0] & bit), bool(s[1] & bit)
        f = _I if not (xb or zb) else (_X if xb and not zb else (_Z if zb and not xb else _iY))
        M = np.kron(M, f)
    return M


def numeric_dim(gens, n):
    A = [dense_real(g, n) for g in gens]
    vecs = [a.ravel() for a in A]

    def rank(vs):
        V = np.stack(vs)
        return int(np.linalg.matrix_rank(V @ V.T))

    basis = list(vecs)
    r = rank(basis)
    while True:
        cands = list(basis)
        for a in A:
            for b in basis:
                B = b.reshape(a.shape)
                cands.append((a @ B - B @ a).ravel())
        # keep an independent subset by incremental Gram-Schmidt (independent of any rtol on svd)
        Q = []
        for v in cands:
            w = v.copy()
            for q in Q:
                w -= (w @ q) * q
            nw = np.linalg.norm(w)
            if nw > 1e-8 * max(1.0, np.linalg.norm(v)):
                Q.append(w / nw)
        r2 = len(Q)
        if r2 == r:
            return r
        basis, r = Q, r2


def main():
    from s24 import stats_lib as ST
    out = {"fixed": {}, "numeric": {}, "conventions": {}}
    t0 = time.perf_counter()
    for n in (4, 5, 6, 7):
        for L in (1, 2, 3, 4):
            g = generators(n, L)
            d, exceeded, S = closure(g)
            out["fixed"][f"n{n}_L{L}"] = {"dim": d, "dim_so": dim_so(n), "full": d == dim_so(n),
                                          "all_odd_y": all(odd_y(s) for s in S), "n_gens": len(set(g))}
            print(f"n={n} L={L} P={n*L:2d} dim {d:6d} so(2^n) {dim_so(n):6d} full={d == dim_so(n)} odd-Y={out['fixed'][f'n{n}_L{L}']['all_odd_y']}")
    for n in (7,):
        for L in (2, 3):
            dims = {}
            for conv, order in itertools.product(("right", "left"), ("forward", "reverse")):
                dims[f"{conv}_{order}"] = closure(generators(n, L, conv, order))[0]
            out["conventions"][f"n{n}_L{L}"] = dims
            print(f"conventions n={n} L={L}: {dims}")
    for n in (4, 5):
        for L in (1, 2, 3):
            g = generators(n, L)
            sym = closure(g)[0]
            num = numeric_dim(g, n)
            out["numeric"][f"n{n}_L{L}"] = {"symbolic": sym, "numeric_gram": num, "agree": sym == num}
            print(f"numeric n={n} L={L}: symbolic {sym} numeric {num} agree {sym == num}")
    out["wall_s"] = round(time.perf_counter() - t0, 2)
    out["reconciliation"] = {"dim_so_128": dim_so(7), "1025_over_8128": 1025 / dim_so(7)}
    ST.save_atomic(OUT, out, module_file=__file__)
    print("wrote", os.path.relpath(OUT, ROOT), "wall", out["wall_s"], "s")


if __name__ == "__main__":
    main()
