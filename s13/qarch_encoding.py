"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 4 -- the encoding study.

Four candidate encodings of a k-state-per-residue torsion space, compared on qubits,
configuration count, feasibility, ansatz EXPRESSIVITY and -- the point the brief insists
on -- how the encoding interacts with LOCALITY of moves.

    binary        log2(k) qubits/residue, state index in plain binary. Whole register
                  feasible, no penalty. A single-qubit flip changes the state index by a
                  power of two, i.e. it can move the residue between torsion states that
                  are arbitrarily far apart on the Ramachandran map.
    gray          same qubits, same feasibility; index -> Gray code, so a single-qubit flip
                  moves to an ADJACENT state index. Only helps if adjacent index means
                  adjacent torsion, which for k-means centres in arbitrary order it does
                  NOT -- so `gray_sorted` reorders the library by a 1-D angular traversal
                  first, and that ordering is the intervention actually under test.
    onehot        k qubits/residue, so n*k qubits; only k^n of 2^(nk) register states are
                  feasible; needs a penalty (or a constrained ansatz). A state change is a
                  2-qubit move, but it never leaves the residue's own block.
    hierarchical  same qubits as binary; the HIGH bits pick a coarse Ramachandran basin
                  (2-means on the k states in the (cos,sin) embedding, recursively), the
                  LOW bits refine inside it. Same feasibility as binary, but the bit-to-
                  geometry map is graded: high bits carry the big moves.

WHAT IS MEASURED

 (a) qubits, configurations, feasible fraction of the register, penalty requirement.
 (b) MOVE LOCALITY, empirically: flip one qubit of a random configuration and measure the
     median torsion change of the affected residue, the median CA-RMSD move of the whole
     structure, and the median |dE| for Legacy and for the prior.
 (c) EXPRESSIVITY of a PRODUCT ansatz, exactly. A single RY layer with no entanglement
     gives an independent Bernoulli per qubit, so the per-residue distributions it can
     represent are exactly the rank-1 (independent-bits) ones. The minimum KL from a target
     per-residue distribution pi to that family is EXACTLY the mutual information between
     the bits of the state index under the encoding's labelling, so the encoding's cost is
     computable in closed form and differs between binary / gray / hierarchical. One-hot
     with post-selection on the feasible subspace is a FULL-RANK family and its minimum KL
     is exactly 0 -- that is a real, provable advantage of one-hot at one layer.
 (d) the same expressivity question for the REAL ansatz at depth: `StatevectorCircuit`
     (RY/CNOT, exact) is fitted to a deliberately MULTIMODAL target distribution over
     configurations of a real 18-qubit instance, and the achieved KL is reported per layer
     count and per encoding.

Output: `s13/results/qarch_encoding.json`.

    python -m s13.qarch_encoding
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402
from core import quantum as qm          # noqa: E402


# ------------------------------------------------------------------ labellings
def gray(k):
    """state index -> gray code word (as an integer)."""
    s = np.arange(k)
    return s ^ (s >> 1)


def angular_order(PHI_i, PSI_i):
    """Order the k states of one residue along a 1-D traversal of the (phi,psi) torus.

    Nearest-neighbour chain from the state with the most negative phi: a cheap, fully
    deterministic ordering that makes 'adjacent index' mean 'adjacent torsion'.
    """
    k = len(PHI_i)
    X = np.column_stack([np.cos(PHI_i), np.sin(PHI_i), np.cos(PSI_i), np.sin(PSI_i)])
    left = list(range(k))
    cur = int(np.argmin(PHI_i))
    order = [cur]; left.remove(cur)
    while left:
        d = ((X[left] - X[cur]) ** 2).sum(1)
        cur = left[int(np.argmin(d))]
        order.append(cur); left.remove(cur)
    return np.array(order)


def hierarchical_order(PHI_i, PSI_i):
    """Recursive 2-means on the (cos,sin) embedding: high bit = coarse basin."""
    X = np.column_stack([np.cos(PHI_i), np.sin(PHI_i), np.cos(PSI_i), np.sin(PSI_i)])

    def rec(idx):
        if len(idx) <= 1:
            return list(idx)
        Y = X[idx]
        c = Y[[0, int(np.argmax(((Y - Y[0]) ** 2).sum(1)))]]
        for _ in range(30):
            lab = ((Y[:, None] - c[None]) ** 2).sum(-1).argmin(1)
            if not lab.any() or lab.all():          # degenerate split -> halve by order
                lab = (np.arange(len(idx)) >= len(idx) // 2).astype(int)
                break
            nc = np.stack([Y[lab == 0].mean(0), Y[lab == 1].mean(0)])
            if np.allclose(nc, c):
                break
            c = nc
        a = [idx[t] for t in np.flatnonzero(lab == 0)]
        b = [idx[t] for t in np.flatnonzero(lab == 1)]
        while len(a) > len(idx) // 2:
            b.insert(0, a.pop())
        while len(b) > len(idx) // 2:
            a.append(b.pop(0))
        return rec(a) + rec(b)

    return np.array(rec(list(range(len(PHI_i)))))


def perm_for(space, kind):
    """(n, k) array: perm[i, code] = the library state that codeword `code` denotes."""
    n, k = space.n, space.k
    P = np.zeros((n, k), int)
    for i in range(n):
        if kind == "binary":
            P[i] = np.arange(k)
        elif kind == "gray":
            P[i, gray(k)] = np.arange(k)
        elif kind == "gray_sorted":
            o = angular_order(space.PHI[i], space.PSI[i])
            P[i, gray(k)] = o
        elif kind == "hierarchical":
            P[i] = hierarchical_order(space.PHI[i], space.PSI[i])
        elif kind == "onehot":
            P[i] = np.arange(k)
        else:
            raise ValueError(kind)
    return P


# ------------------------------------------------------------------ (b) move locality
def move_locality(space, kinds, n_base=64, seed=0):
    """Flip ONE qubit; measure what moves. Binary/gray/hier share a register layout."""
    rng = np.random.default_rng(seed)
    n, k = space.n, space.k
    b = space.bits_per_res
    base_state = space.uniform(n_base, rng)
    Pemp = Q.empirical_prior(space)
    out = {}
    for kind in kinds:
        P = perm_for(space, kind)
        inv = np.argsort(P, axis=1)                      # state -> codeword
        rows = []
        if kind == "onehot":
            # a state change is a 2-qubit move (clear one, set another); a 1-qubit flip is
            # always infeasible. Measure the 2-qubit move, and count the penalty need.
            pairs = [(a, c) for a in range(k) for c in range(k) if a != c]
            for i in range(n):
                for (a, c) in pairs:
                    m = base_state[:, i] == a
                    if not m.any():
                        continue
                    S2 = base_state[m].copy(); S2[:, i] = c
                    rows.append(_move_stats(space, base_state[m], S2, i, Pemp, 2))
        else:
            for i in range(n):
                for bit in range(b):
                    code = inv[i][base_state[:, i]]
                    code2 = code ^ (1 << bit)
                    S2 = base_state.copy(); S2[:, i] = P[i][code2]
                    rows.append(_move_stats(space, base_state, S2, i, Pemp, 1))
        agg = {kk: float(np.median([r[kk] for r in rows])) for kk in rows[0]}
        agg["mean_index_hamming_for_state_change"] = _mean_hamming(P, b)
        out[kind] = agg
    return out


def _move_stats(space, S1, S2, i, Pemp, nflip):
    d_phi = np.abs(Q.wrap(space.PHI[i][S2[:, i]] - space.PHI[i][S1[:, i]]))
    d_psi = np.abs(Q.wrap(space.PSI[i][S2[:, i]] - space.PSI[i][S1[:, i]]))
    ca1, ca2 = space.ca(S1), space.ca(S2)
    from s12 import instrument as I
    rm = np.array([I.ca_rmsd(a, bb) for a, bb in zip(ca1, ca2)])
    e1 = Q.legacy_energy(space, S1); e2 = Q.legacy_energy(space, S2)
    p1 = Q.prior_energy(Pemp, S1); p2 = Q.prior_energy(Pemp, S2)
    return {"qubit_flips": float(nflip),
            "torsion_move_deg": float(np.median(np.degrees(d_phi + d_psi))),
            "structure_move_rmsd_A": float(np.median(rm)),
            "abs_dE_legacy": float(np.median(np.abs(e2 - e1))),
            "abs_dE_prior": float(np.median(np.abs(p2 - p1)))}


def _mean_hamming(P, b):
    """Mean Hamming distance between the codewords of two DIFFERENT states."""
    n, k = P.shape
    inv = np.argsort(P, axis=1)
    tot, cnt = 0.0, 0
    for i in range(n):
        for a in range(k):
            for c in range(k):
                if a == c:
                    continue
                tot += bin(int(inv[i, a]) ^ int(inv[i, c])).count("1"); cnt += 1
    return tot / cnt


# ------------------------------------------------------------------ (c) exact expressivity
def product_ansatz_gap(space, kinds, targets):
    """min KL(pi || independent bits) per residue = mutual information under the labelling.

    EXACT and closed form.  `targets` maps a name to an (n, k) distribution.
    """
    n, k, b = space.n, space.k, space.bits_per_res
    out = {}
    for tname, PI in targets.items():
        out[tname] = {}
        for kind in kinds:
            if kind == "onehot":
                # post-selected product Bernoulli over k wires is full rank on the feasible
                # subspace: p(s) ~ q_s prod_{t!=s}(1-q_t) reaches any pi exactly.
                out[tname]["onehot"] = {"mean_KL_nats": 0.0, "max_KL_nats": 0.0,
                                        "exact": True,
                                        "note": "post-selected product Bernoulli is "
                                                "full-rank on the feasible subspace"}
                continue
            P = perm_for(space, kind)
            inv = np.argsort(P, axis=1)
            kl = []
            for i in range(n):
                pi = PI[i][P[i]]                     # probability of each CODEWORD
                T = pi.reshape((2,) * b)
                m = 0.0
                for bit in range(b):
                    pass
                # mutual information among the b bits = KL(joint || product of marginals)
                marg = [T.sum(axis=tuple(j for j in range(b) if j != bit))
                        for bit in range(b)]
                prod = np.ones_like(T)
                for bit in range(b):
                    shape = [1] * b; shape[bit] = 2
                    prod = prod * marg[bit].reshape(shape)
                mask = T > 0
                m = float((T[mask] * np.log(T[mask] / prod[mask])).sum())
                kl.append(m)
            out[tname][kind] = {"mean_KL_nats": float(np.mean(kl)),
                                "max_KL_nats": float(np.max(kl)), "exact": True}
    return out


# ------------------------------------------------------------------ (d) real ansatz fit
def register_index(space, S, kind, residues=None):
    """Configuration -> computational-basis index, under a given labelling."""
    P = perm_for(space, kind)
    inv = np.argsort(P, axis=1)
    b = space.bits_per_res
    S = np.atleast_2d(np.asarray(S, int))
    res = range(space.n) if residues is None else residues
    idx = np.zeros(len(S), dtype=np.int64)
    for c, i in enumerate(res):
        idx = (idx << b) | inv[i][S[:, c if residues is not None else i]]
    return idx


def ansatz_fit(nq, target_p, layers_list=(1, 2, 3, 4, 6), iters=250, seed=0, lr=0.1):
    """Fit the genuine RY/CNOT `StatevectorCircuit` to a target distribution.

    Minimises KL(target || p_theta) with Adam on the ANALYTIC parameter-shift gradient of
    the exact probabilities -- no sampling anywhere, so an achieved KL is an expressivity
    statement about the ansatz and not about shot noise.
    """
    res = {}
    tgt = np.asarray(target_p, float)
    sup = tgt > 0
    for L in layers_list:
        circ = qm.StatevectorCircuit(nq, layers=L, ring=True)
        rng = np.random.default_rng(seed)
        th = rng.normal(np.pi / 2, 0.3, circ.n_params())
        opt = qm.Adam(circ.n_params(), lr=lr)
        best = np.inf
        for _ in range(iters):
            p = circ.probs(th)
            kl = float((tgt[sup] * np.log(tgt[sup] / np.maximum(p[sup], 1e-300))).sum())
            best = min(best, kl)
            TH = circ._shift_grid(th, np.pi / 2)
            w = -tgt / np.maximum(p, 1e-300)
            g = np.empty(circ.n_params())
            for a in range(0, len(TH), 32):
                PS = circ.probs_batch(TH[a:a + 32])
                for t in range(0, len(PS), 2):
                    g[(a + t) // 2] = float(w @ (0.5 * (PS[t] - PS[t + 1])))
            th = opt.step(th, g)
        res[L] = {"KL_nats": float(best), "n_params": int(circ.n_params())}
    return res


def multimodal_expressivity(pdb_id="1CS9", k=4, n_res=6, n_modes=8, seed=0,
                            kinds=("binary", "gray_sorted", "hierarchical")):
    """Can the RY/CNOT ansatz express a MULTIMODAL distribution over the relevant states?

    ORACLE DIAGNOSTIC in its target only: the target distribution is uniform over the
    `n_modes` lowest-CA-RMSD configurations of a real `n_res`-residue sub-register, which
    is the hardest honest instance of 'the several structural hypotheses worth keeping'.
    The measured quantity -- achieved KL vs depth vs labelling -- is a property of the
    ansatz and the encoding, not a predictive result.
    """
    sp = Q.Space(pdb_id, k)
    res_idx = list(range(n_res))
    grid = np.array(list(itertools.product(range(k), repeat=n_res)), dtype=np.int8)
    base = np.zeros(sp.n, np.int8)
    S = np.repeat(base[None, :], len(grid), 0)
    S[:, res_idx] = grid
    r = sp.rmsd(S)
    modes = np.argsort(r)[:n_modes]
    nq = n_res * sp.bits_per_res
    out = {"pdb": pdb_id, "k": k, "n_res": n_res, "n_qubits": nq, "n_modes": n_modes,
           "mode_rmsd": [float(x) for x in r[modes]], "fits": {}}
    for kind in kinds:
        tgt = np.zeros(1 << nq)
        tgt[register_index(sp, grid[modes], kind, residues=res_idx)] = 1.0 / n_modes
        out["fits"][kind] = ansatz_fit(nq, tgt, seed=seed)
        print(f"    {kind}: " + "  ".join(
            f"L{L}:{v['KL_nats']:.3f}" for L, v in out["fits"][kind].items()), flush=True)
    out["onehot_note"] = (f"one-hot needs {n_res * k} qubits (dim 2^{n_res*k}); the exact "
                          "statevector fit is not affordable, so only the analytic "
                          "one-layer result is reported for it")
    return out


# ------------------------------------------------------------------ driver
def main():
    t0 = time.time()
    kinds = ["binary", "gray", "gray_sorted", "hierarchical", "onehot"]
    rows = []
    for pdb_id in ["1CS9", "1A13"]:
        for k in (4, 8):
            sp = Q.Space(pdb_id, k)
            n, b = sp.n, sp.bits_per_res
            book = {
                "binary": dict(qubits=n * b, configs=float(k) ** n, feasible_frac=1.0,
                               penalty_needed=False, flips_per_state_change="1..log2(k)"),
                "gray": dict(qubits=n * b, configs=float(k) ** n, feasible_frac=1.0,
                             penalty_needed=False, flips_per_state_change="1..log2(k)"),
                "gray_sorted": dict(qubits=n * b, configs=float(k) ** n, feasible_frac=1.0,
                                    penalty_needed=False, flips_per_state_change="1..log2(k)"),
                "hierarchical": dict(qubits=n * b, configs=float(k) ** n, feasible_frac=1.0,
                                     penalty_needed=False, flips_per_state_change="1..log2(k)"),
                "onehot": dict(qubits=n * k, configs=float(k) ** n,
                               feasible_frac=float((k / 2.0 ** k) ** n),
                               penalty_needed=True, flips_per_state_change="2"),
            }
            Pemp = Q.empirical_prior(sp)
            tgt = {"empirical_prior": Pemp,
                   "ORACLE_prior_q0.8": Q.ORACLE_prior(sp, 0.8),
                   "bimodal_alpha_beta": _bimodal(sp)}
            expr = product_ansatz_gap(sp, kinds, tgt)
            mv = move_locality(sp, kinds)
            rows.append({"pdb": pdb_id, "n": n, "k": k, "book": book,
                         "product_ansatz_min_KL": expr, "move_locality": mv})
            print(f"  {pdb_id} k={k}: binary {n*b} qubits, onehot {n*k} qubits, "
                  f"feasible {book['onehot']['feasible_frac']:.2e}; "
                  f"minKL binary {expr['empirical_prior']['binary']['mean_KL_nats']:.4f} "
                  f"gray_sorted {expr['empirical_prior']['gray_sorted']['mean_KL_nats']:.4f} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    print("  multimodal expressivity of the real RY/CNOT ansatz", flush=True)
    mm = multimodal_expressivity()
    Q.write("qarch_encoding", {"what": "encoding study: qubits, feasibility, move locality, "
                                       "exact product-ansatz expressivity, real-ansatz "
                                       "multimodal fit", "rows": rows,
                               "multimodal_expressivity": mm})
    return rows


def _bimodal(space):
    """A deliberately bimodal per-residue target: the two most separated library states."""
    n, k = space.n, space.k
    P = np.full((n, k), 1e-3)
    for i in range(n):
        X = np.column_stack([np.cos(space.PHI[i]), np.sin(space.PHI[i]),
                             np.cos(space.PSI[i]), np.sin(space.PSI[i])])
        D = ((X[:, None] - X[None]) ** 2).sum(-1)
        a, c = np.unravel_index(np.argmax(D), D.shape)
        P[i, a] = P[i, c] = 0.5
    return P / P.sum(1, keepdims=True)


if __name__ == "__main__":
    main()
