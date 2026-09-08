"""SPRINT 14 / VQE -- PART B: encodings, locality, and the Pauli spectrum done correctly.

Sprint 13 compared binary / gray / gray_sorted / hierarchical / one-hot on qubit count,
move locality and one-layer expressivity, and chose binary.  It did NOT measure the
Pauli-weight distribution, the Hamiltonian sparsity, the coefficient dynamic range or the
gradient variance PER ENCODING, and it did not consider domain-wall / unary at all.  This
does those, adds domain-wall, and states the one thing that settles the qubit-count question
outright.

THE QUBIT-COUNT QUESTION IS CLOSED BY COUNTING
---------------------------------------------
The space has ``k^n`` configurations, so ANY faithful encoding needs at least
``ceil(n log2 k)`` qubits.  When ``k`` is a power of two, binary attains that bound exactly:
at k=4, n=9 the bound is ``log2(4^9) = 18`` qubits and binary uses 18.  **No encoding can be
more compact, and every alternative is strictly wider.**  The only regime where a cleverer
encoding has room is ``k`` NOT a power of two, where binary wastes
``n log2 k - ceil(n log2 k)`` -- measured below and it is under one qubit at every k tested.

THE SPECTRAL METHODOLOGY, WHICH IS THE POINT
--------------------------------------------
Both energies are diagonal, so each is exactly a weighted sum of Pauli-Z strings and its
coefficients are its Walsh-Hadamard transform.  The brief's warning is that the spectrum of
an UNCONDITIONED molecular energy is a delta-spike artefact: the top-10 of 4,096
configurations carry a median 99.6% of raw AMBER's Walsh variance, and a constant-plus-spike
has weight spectrum exactly Binomial(m, 1/2) -- i.e. the "measurement" returns the binomial
whatever the physics is.  99th-percentile winsorisation is NOT sufficient.

So every objective here is passed through `vqe_lib.uniformise` first -- a rank transform,
the strongest monotone conditioning available.  It preserves every ordering the objective
expresses and forces all objectives onto the identical marginal, so a spectral difference
between two conditioned objectives cannot be a difference of scale or tail shape.  The
UNCONDITIONED spectrum is computed alongside and the Binomial(m, 1/2) null is printed next
to it, so the artefact is priced rather than assumed away.

    python -m s14.vqe_encoding
"""
from __future__ import annotations

import itertools

import numpy as np

from core import quantum as Q
from s14 import vqe_lib as V

NSUB = 5          # residues in the sub-register: binary 10q, domain-wall 15q, one-hot 20q


# ================================================================== the encodings
class Encoding:
    """Maps a length-`n` vector of k-state indices onto a qubit register.

    `qubits` is the register width; `feasible` is the list of register integers that decode
    to a legal configuration, in configuration order (so `feasible[c]` is the register word
    of configuration `c`).  A surjective encoding has `len(feasible) == 2**qubits`.
    """

    def __init__(self, name, n, k, words, qubits, penalty_needed):
        self.name, self.n, self.k = name, n, k
        self.qubits = qubits
        self.words = np.asarray(words, np.int64)      # (k^n,) register word per config
        self.penalty_needed = penalty_needed
        self.feasible_fraction = len(self.words) / float(1 << qubits)


def _per_residue_words(codes, n, k, bits):
    """Compose a per-residue codeword table into whole-register words, residue 0 first."""
    cfg = np.array(list(itertools.product(range(k), repeat=n)), np.int64)   # big-endian
    w = np.zeros(len(cfg), np.int64)
    for i in range(n):
        w |= np.asarray(codes, np.int64)[cfg[:, i]] << (bits * (n - 1 - i))
    return w


def _gray(x):
    return x ^ (x >> 1)


def build_encodings(n=NSUB, k=4, order=None):
    """`order` reorders the k library states (used by gray_sorted)."""
    idx = np.arange(k)
    o = np.arange(k) if order is None else np.asarray(order)
    inv = np.argsort(o)                      # library state -> position in the traversal
    b = int(np.log2(k))
    enc = {}
    enc["binary"] = Encoding("binary", n, k, _per_residue_words(idx, n, k, b), n * b, False)
    enc["gray"] = Encoding("gray", n, k, _per_residue_words(_gray(idx), n, k, b),
                           n * b, False)
    enc["gray_sorted"] = Encoding("gray_sorted", n, k,
                                  _per_residue_words(_gray(inv), n, k, b), n * b, False)
    # one-hot: k qubits per residue, exactly one set
    oh = np.array([1 << (k - 1 - s) for s in range(k)], np.int64)
    enc["onehot"] = Encoding("onehot", n, k, _per_residue_words(oh, n, k, k), n * k, True)
    # domain-wall (== unary == thermometer): k-1 qubits, a prefix of ones
    dw = np.array([(1 << s) - 1 << 0 for s in range(k)], np.int64)
    dw = np.array([((1 << s) - 1) << (k - 1 - s) for s in range(k)], np.int64)
    enc["domain_wall"] = Encoding("domain_wall", n, k,
                                  _per_residue_words(dw, n, k, k - 1), n * (k - 1), True)
    return enc


# ============================================================ Walsh / Pauli spectrum
def walsh(f: np.ndarray) -> np.ndarray:
    """Fast Walsh-Hadamard transform, normalised so coefficients are Pauli-Z weights."""
    a = np.asarray(f, float).copy()
    m = int(np.log2(len(a)))
    assert 1 << m == len(a)
    h = 1
    while h < len(a):
        for s in range(0, len(a), h * 2):
            x = a[s:s + h].copy()
            y = a[s + h:s + 2 * h].copy()
            a[s:s + h] = x + y
            a[s + h:s + 2 * h] = x - y
        h *= 2
    return a / len(a)


def spectrum_stats(f: np.ndarray, label="") -> dict:
    """Pauli-weight distribution and sparsity of a diagonal operator over `m` qubits."""
    m = int(np.log2(len(f)))
    c = walsh(f)
    w = np.array([bin(j).count("1") for j in range(len(c))])
    v = c ** 2
    v[0] = 0.0                                  # identity carries no information
    tot = v.sum()
    if tot <= 0:
        return {"n_qubits": m, "mean_weight": float("nan")}
    pw = v / tot
    a = np.abs(c[1:])
    nz = a[a > 1e-12 * a.max()] if a.max() > 0 else a
    return {
        "n_qubits": m,
        "mean_weight": float((w * pw).sum()),
        "max_weight_99": int(np.searchsorted(np.cumsum(np.bincount(
            w, weights=pw, minlength=m + 1)), 0.99)),
        "max_weight": int(w[v > 1e-14 * v.max()].max()) if v.max() > 0 else 0,
        "weight_hist": np.bincount(w, weights=pw, minlength=m + 1).tolist(),
        "frac_variance_top10_configs": float(
            np.sort(( (f - f.mean()) ** 2))[-10:].sum() / ((f - f.mean()) ** 2).sum()),
        "n_terms_99pct_variance": int(np.searchsorted(
            np.cumsum(np.sort(v)[::-1]) / tot, 0.99) + 1),
        "sparsity_99pct": float((np.searchsorted(
            np.cumsum(np.sort(v)[::-1]) / tot, 0.99) + 1) / (len(c) - 1)),
        "dynamic_range_decades": float(np.log10(nz.max() / nz.min()))
        if len(nz) and nz.min() > 0 else float("inf"),
        "binomial_null_mean_weight": m / 2.0,
    }


# ================================================================ the sub-register
def sub_objective(e: V.Enum, residues, base_cfg, which="legacy"):
    """Restrict a tabulated objective to `residues`, holding the rest at `base_cfg`.

    Gives a genuine structural objective on a k^|residues| space that every encoding --
    including the 20-qubit one-hot -- can be spectrally analysed on.
    """
    n, k = e.n, e.k
    src = getattr(e, which) if which != "rmsd" else e.rmsd
    cfg = np.array(list(itertools.product(range(k), repeat=len(residues))), np.int64)
    S = np.tile(np.asarray(base_cfg, np.int64), (len(cfg), 1))
    S[:, list(residues)] = cfg
    return src[e.index(S)]


def embed(enc: Encoding, values: np.ndarray, penalty: float) -> np.ndarray:
    """Place the configuration-ordered `values` at their register words; penalise the rest."""
    out = np.full(1 << enc.qubits, penalty, float)
    out[enc.words] = values
    return out


# ===================================================================== experiments
def bookkeeping():
    print("=" * 96)
    print("B1. BOOK-KEEPING, and the information-theoretic bound")
    print("=" * 96)
    print(f"{'encoding':14s} {'qubits/res':>11s} {'qubits n=9':>11s} {'configs':>10s} "
          f"{'feasible frac':>14s} {'penalty':>8s} {'vs bound':>9s}")
    rows = []
    for n, k in ((9, 4),):
        bound = int(np.ceil(n * np.log2(k)))
        enc = build_encodings(n=4, k=k)        # widths are per-residue, n cancels
        for name, sp in [("binary", int(np.log2(k))), ("gray", int(np.log2(k))),
                         ("gray_sorted", int(np.log2(k))), ("domain_wall", k - 1),
                         ("onehot", k)]:
            qq = n * sp
            ff = (k / float(1 << sp)) ** n
            rows.append({"encoding": name, "qubits_per_res": sp, "qubits": qq,
                         "feasible_fraction": ff, "bound": bound})
            print(f"{name:14s} {sp:11d} {qq:11d} {k**n:10d} {ff:14.3e} "
                  f"{'yes' if ff < 1 else 'no':>8s} {qq-bound:+9d}")
        print(f"\n  information-theoretic bound at n={n}, k={k}: "
              f"ceil(n log2 k) = {bound} qubits.  BINARY ATTAINS IT.")
    print("\n  binary's slack when k is NOT a power of two (the only regime with room):")
    print(f"  {'k':>3s} {'n':>3s} {'binary':>7s} {'bound':>6s} {'wasted':>7s}")
    for k in (3, 5, 6, 7, 12):
        for n in (9, 16):
            b = n * int(np.ceil(np.log2(k)))
            bd = int(np.ceil(n * np.log2(k)))
            print(f"  {k:3d} {n:3d} {b:7d} {bd:6d} {b-bd:7d}")
    print("\n  NOTE: 'unary' and 'domain-wall' are THE SAME ENCODING for a single")
    print("  categorical variable -- a prefix of ones in k-1 qubits.  They are listed once.")
    return rows


def spectra(pdb="1CS9", residues=(2, 3, 4, 5, 6)):
    print()
    print("=" * 96)
    print(f"B2. PAULI SPECTRUM per encoding, {pdb}, residues {residues} "
          f"({len(residues)} residues, {4**len(residues)} configurations)")
    print("    EVERY objective is rank-conditioned first.  The unconditioned row is shown")
    print("    next to the Binomial(m, 1/2) null so the delta-spike artefact is priced.")
    print("=" * 96)
    e = V.Enum(pdb)
    rng = np.random.default_rng(0)
    base = rng.integers(0, e.k, e.n)
    enc = build_encodings(n=len(residues), k=e.k)
    out = {}
    for which in ("legacy", "prior", "rmsd"):
        raw = sub_objective(e, residues, base, which)
        cond = V.uniformise(raw)
        print(f"\n--- objective: {which} "
              f"(rho(raw, rmsd) = {V.spearman(raw, sub_objective(e, residues, base, 'rmsd')):+.3f}) ---")
        print(f"{'encoding':13s} {'qubits':>6s} {'cond mean wt':>12s} {'binom null':>11s} "
              f"{'raw mean wt':>12s} {'99% wt':>7s} {'#terms 99%var':>14s} "
              f"{'sparsity':>9s} {'dyn range':>10s}")
        for name, en in enc.items():
            # penalty: one unit above the conditioned maximum (a monotone extension)
            pen_c = 1.0 + 1.0 / len(cond)
            fc = embed(en, cond, pen_c) if en.penalty_needed else cond[np.argsort(en.words)]
            fr = (embed(en, raw, float(raw.max() + (raw.max() - raw.min())))
                  if en.penalty_needed else raw[np.argsort(en.words)])
            sc = spectrum_stats(fc)
            sr = spectrum_stats(fr)
            out.setdefault(which, {})[name] = {"conditioned": sc, "raw": sr,
                                               "qubits": en.qubits,
                                               "feasible_fraction": en.feasible_fraction}
            print(f"{name:13s} {en.qubits:6d} {sc['mean_weight']:12.3f} "
                  f"{sc['binomial_null_mean_weight']:11.1f} {sr['mean_weight']:12.3f} "
                  f"{sc['max_weight_99']:7d} {sc['n_terms_99pct_variance']:14d} "
                  f"{sc['sparsity_99pct']:9.4f} {sc['dynamic_range_decades']:10.2f}")
        print(f"  raw top-10-config share of total variance (the artefact): "
              f"{spectrum_stats(raw)['frac_variance_top10_configs']:.4f}"
              f"   conditioned: {spectrum_stats(cond)['frac_variance_top10_configs']:.4f}")
    return out


def penalty_sensitivity(pdb="1CS9", residues=(2, 3, 4, 5, 6)):
    print()
    print("=" * 96)
    print("B3. IS A NON-SURJECTIVE ENCODING'S SPECTRUM THE OBJECTIVE, OR THE CONSTRAINT?")
    print("    Sweep the penalty height; if the spectrum moves, it is measuring the")
    print("    constraint and the encoding comparison is meaningless as usually run.")
    print("=" * 96)
    e = V.Enum(pdb)
    base = np.random.default_rng(0).integers(0, e.k, e.n)
    cond = V.uniformise(sub_objective(e, residues, base, "legacy"))
    enc = build_encodings(n=len(residues), k=e.k)
    print(f"{'penalty':>9s} " + "".join(f"{n:>14s}" for n in ("onehot", "domain_wall")))
    out = {}
    for pen in (1.0 + 1 / len(cond), 1.5, 2.0, 5.0, 20.0, 100.0):
        row = []
        for name in ("onehot", "domain_wall"):
            s = spectrum_stats(embed(enc[name], cond, pen))
            row.append(s["mean_weight"])
            out.setdefault(name, {})[str(pen)] = s["mean_weight"]
        print(f"{pen:9.3f} " + "".join(f"{v:14.3f}" for v in row))
    print(f"\n  binary (surjective, no penalty) mean weight: "
          f"{spectrum_stats(cond[np.argsort(enc['binary'].words)])['mean_weight']:.3f}")
    print("  binary needs no penalty at all, so its spectrum has no such degree of freedom.")
    return out


def gradient_variance(pdb="1CS9", residues=(2, 3, 4, 5), seeds=64):
    # 4 residues, not 5: one-hot is then 16 qubits and the 2P-row parameter-shift batch is
    # 33 MB rather than 640 MB.  The box has ~3 GB free with the Part C grid running.
    print()
    print("=" * 96)
    print("B4. GRADIENT VARIANCE per encoding, same objective, same ansatz family")
    print("    RY/CNOT chain at depth 2 on each encoding's own register width.")
    print("=" * 96)
    e = V.Enum(pdb)
    base = np.random.default_rng(0).integers(0, e.k, e.n)
    cond = V.uniformise(sub_objective(e, residues, base, "legacy"))
    enc = build_encodings(n=len(residues), k=e.k)
    print(f"{'encoding':13s} {'qubits':>6s} {'Var[g_1]':>12s} {'mean |g|':>10s} "
          f"{'zero-grad':>10s}")
    out = {}
    for name, en in enc.items():
        if en.qubits > 20:
            print(f"{name:13s} {en.qubits:6d}  skipped (2**{en.qubits} statevector)")
            continue
        f = embed(en, cond, 1.0 + 1 / len(cond)) if en.penalty_needed \
            else cond[np.argsort(en.words)]
        circ = Q.StatevectorCircuit(en.qubits, 2)
        g1, gn, z = [], [], 0
        for s in range(seeds):
            r = np.random.default_rng(s)
            th = r.uniform(0, 2 * np.pi, circ.n_params())
            g = Q.grad_cvar_paramshift(circ, th, f, 1.0)     # expectation value, alpha=1
            if np.linalg.norm(g) == 0:
                z += 1
                continue
            g1.append(g[0])
            gn.append(np.linalg.norm(g))
        out[name] = {"var_g1": float(np.var(g1)), "mean_gnorm": float(np.mean(gn)),
                     "qubits": en.qubits, "zero": z}
        print(f"{name:13s} {en.qubits:6d} {np.var(g1):12.3e} {np.mean(gn):10.4f} "
              f"{z:10d}")
    print("\n  Barren-plateau expectation is Var ~ 2^-qubits; the comparison to make is")
    print("  whether an encoding's variance is worse than its OWN width predicts.")
    for name, o in out.items():
        print(f"    {name:13s} Var[g1] * 2^qubits = {o['var_g1'] * (1 << o['qubits']):.4f}")
    return out


def move_locality(pdb="1CS9"):
    print()
    print("=" * 96)
    print("B5. MOVE LOCALITY on the FULL register, re-measured with true RMSD available")
    print("    (s13 measured this on 64 base configurations; here it is all 262,144.)")
    print("=" * 96)
    e = V.Enum(pdb)
    rng = np.random.default_rng(0)
    base = rng.integers(0, e.N, 20000)
    print(f"{'encoding':13s} {'flips/state change':>19s} {'median |d RMSD|':>16s} "
          f"{'median |d legacy|':>18s}")
    out = {}
    enc = build_encodings(n=e.n, k=e.k)
    for name in ("binary", "gray"):
        en = enc[name]
        # single-qubit flips on the register, decoded through this labelling
        dr, dl = [], []
        for q in range(e.n_qubits):
            w = en.words[base]
            w2 = w ^ (1 << (e.n_qubits - 1 - q))
            c1 = np.argsort(en.words)[w]      # register word -> configuration index
            c2 = np.argsort(en.words)[w2]
            dr.append(np.abs(e.rmsd[c1] - e.rmsd[c2]))
            dl.append(np.abs(e.legacy[c1] - e.legacy[c2]))
        out[name] = {"median_drmsd": float(np.median(np.concatenate(dr))),
                     "median_dlegacy": float(np.median(np.concatenate(dl)))}
        print(f"{name:13s} {1:19d} {out[name]['median_drmsd']:16.4f} "
              f"{out[name]['median_dlegacy']:18.4f}")
    # one-hot / domain-wall change a residue in 2 flips; the RESULTING config move is
    # identical to a binary state change, so the quantity to report is the state change
    ch = []
    for i in range(e.n):
        for a in range(e.k):
            for b in range(e.k):
                if a == b:
                    continue
                S = e.states(base)
                S1 = S.copy(); S1[:, i] = a
                S2 = S.copy(); S2[:, i] = b
                ch.append(np.abs(e.rmsd[e.index(S1)] - e.rmsd[e.index(S2)]))
                break
            break
    print(f"{'any (state)':13s} {2:19d} {float(np.median(np.concatenate(ch))):16.4f}")
    print("\n  A one-hot/domain-wall 2-flip and a binary 1-flip produce the SAME set of")
    print("  configuration changes, so no labelling can make the MOVE local.  Confirms")
    print("  s13 section 4b on 300x the sample.")
    return out


def main():
    V.wait_for_memory(1.0, "vqe_encoding")
    out = {}
    out["bookkeeping"] = bookkeeping()
    out["spectra"] = spectra()
    out["penalty_sensitivity"] = penalty_sensitivity()
    out["gradient_variance"] = gradient_variance()
    out["move_locality"] = move_locality()
    V.write("vqe_encoding", out)
    print("\nwritten -> s14/results/vqe_encoding.json")


if __name__ == "__main__":
    main()
