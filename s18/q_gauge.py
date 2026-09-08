"""SPRINT 18 / QUANTUM-ADVERSARIAL -- A3: the degree-1 result is GAUGE-DEPENDENT.

THE ARGUMENT, WHICH IS EXACT AND NOT STATISTICAL.

The lattice encodes k = 4 torsion states of a residue in 2 qubits.  Which bit pattern names
which torsion state is a LABELLING CHOICE.  Nothing physical depends on it: the objective
`E`, its argmin, its Spearman with RMSD, the RMSD of every configuration, and the residue-
additive ANOVA object `RA` are all INVARIANT under relabelling the k states of a residue.

The strict Walsh weight-<=1 object `W1` is NOT.

    per residue i, write  m_i[v] = E_mu[ E | residue i = v ]      (v = 0..k-1)
    * `RA`  picks   v*_i = argmin_v m_i[v]                        <- LABEL-FREE
    * `W1`  keeps only the two SINGLE-QUBIT Walsh coefficients of the 4-vector m_i and picks
            the two bits INDEPENDENTLY by their signs                <- LABEL-DEPENDENT

(That `W1`'s weight-1 coefficients are functions of `m_i` alone is exact: chi_q depends only
on residue i's state, so c_q = (1/k) sum_v m_i[v] chi_q(v).  Verified numerically below
against the full FWHT.)

So `W1` throws away exactly the intra-residue weight-2 coefficient of `m_i` -- the one that
says whether the 4-vector's minimum is at a bit pattern the two independent bit choices can
reach.  Relabelling the states permutes which of the 4 values sits at which corner, and
therefore MOVES `W1`'s argmin, while moving nothing physical at all.

THE TEST.  The gauge group is (S_k)^n -- 24^n relabellings.  The per-residue choice v*_i under
a uniformly random labelling has an exactly computable 24-point distribution, and the choices
are independent across residues, so the induced distribution over CONFIGURATIONS is sampled
exactly.  If the 2.411 A that Sprint 17 reported is a property of the objective, it must be
invariant.  If it moves under a relabelling, it is a property of the ENCODING, has no
continuous-torsion analogue at all, and F5 fires.

BOTH MANDATORY CONTROLS.  Zero-information reference: a uniformly random configuration.
Matched-random: a configuration built by picking each residue's state uniformly at random
from its k options -- the same construction with the objective's information deleted.

RUN:  python -m s18.q_gauge run
      python -m s18.q_gauge report
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s14 import vqe_lib as V                       # noqa: E402
from s15 import seed as SD                         # noqa: E402
from s16 import qphase_lib as QP                   # noqa: E402
from s17 import q_theory as T                      # noqa: E402
from s18 import q_anova as A                       # noqa: E402

RESULTS = os.path.join(ROOT, "s18", "results")
os.makedirs(RESULTS, exist_ok=True)
SALT = "s18quantum"
NSAMP = 4000          # gauge samples per target


def marginals(ins, uniformised=True):
    """m[i, v] = E_mu[E | residue i = v], mu uniform on the enumerated lattice. EXACT."""
    E = ins.hamil().astype(np.float64)
    if uniformised:
        E = V.uniformise(E)
    n, k = ins.n, ins.k
    pw = k ** np.arange(n - 1, -1, -1)
    idx = np.arange(ins.N, dtype=np.int64)
    m = np.empty((n, k))
    for i in range(n):
        s = (idx // pw[i]) % k
        m[i] = np.bincount(s, weights=E, minlength=k) / np.bincount(s, minlength=k)
        del s
    return m, pw, E


def w1_pick(f4, bits=2):
    """Given the k-vector `f4` INDEXED BY CODE, return the code the strict Walsh weight-<=1
    object minimises: each bit chosen independently by the sign of its own coefficient."""
    k = f4.size
    codes = np.arange(k)
    out = 0
    for b in range(bits):
        chi = 1.0 - 2.0 * ((codes >> b) & 1)          # (-1)^{bit b}
        c = float((f4 * chi).mean())
        # minimise c * chi  ->  chi = -sign(c)  ->  bit = 0 if c < 0 else 1
        out |= (0 if c < 0 else 1) << b
    return int(out)


def gauge_cell(pdb, nsamp=NSAMP, uniformised=True):
    ins = QP.inst(pdb)
    rmsd = ins.rmsd
    m, pw, E = marginals(ins, uniformised)
    n, k = ins.n, ins.k
    bits = ins.bits_per_res
    perms = list(itertools.permutations(range(k)))       # the gauge group of one residue
    G = len(perms)

    # ---- for each residue, the state chosen by W1 under each of the G labellings
    #      labelling sigma sends STATE v -> CODE sigma[v]; f4[code] = m[i, sigma^{-1}(code)]
    chosen = np.empty((n, G), int)
    for i in range(n):
        for g, sg in enumerate(perms):
            inv = np.empty(k, int)
            for v, c in enumerate(sg):
                inv[c] = v
            f4 = m[i][inv]
            code = w1_pick(f4, bits)
            chosen[i, g] = int(inv[code])
    ra_state = np.argmin(m, axis=1)
    agree = float((chosen == ra_state[:, None]).mean())

    # ---- identity labelling reproduces Sprint 17's W1 exactly?
    ident = perms.index(tuple(range(k)))
    w1_state_ident = chosen[:, ident]
    idx_ident = int(w1_state_ident @ pw)

    # ---- sample the gauge group: residues independent
    rng = SD.stable_rng(pdb, "gauge", salt=SALT)
    gs = rng.integers(0, G, size=(nsamp, n))
    st = chosen[np.arange(n)[None, :], gs]                      # (nsamp, n)
    idxs = (st * pw[None, :]).sum(1)
    r_gauge = rmsd[idxs]

    # ---- MATCHED-RANDOM control: same construction, objective information deleted
    st_r = rng.integers(0, k, size=(nsamp, n))
    r_rand = rmsd[(st_r * pw[None, :]).sum(1)]

    ra_idx = int(ra_state @ pw)
    ident_r = float(rmsd[idx_ident])
    # MID-RANK percentile: `<` alone reads 0.000 on a degenerate (all-equal) gauge orbit,
    # which is the tie-breaking trap the ledger records.  Ties are split, so a fully
    # degenerate orbit correctly reports 0.500 (the identity is not special) rather than 0.
    pct_mid = float((r_gauge < ident_r).mean() + 0.5 * (r_gauge == ident_r).mean())
    out = {
        "pdb": pdb, "n": n, "k": k, "fold": ins.fold, "N": int(ins.N),
        "uniformised": bool(uniformised), "nsamp": int(nsamp), "gauge_group_per_residue": G,
        "gauge_group_size_log10": float(n * np.log10(G)),
        "W1_identity_argmin_rmsd_ORACLE": float(rmsd[idx_ident]),
        "RA_argmin_rmsd_ORACLE": float(rmsd[ra_idx]),
        "full_argmin_rmsd_ORACLE": float(rmsd[E == E.min()].mean()),
        "space_best_ORACLE": float(rmsd.min()),
        "space_mean_ORACLE": float(rmsd.mean()),
        "frac_residue_labellings_where_W1_agrees_with_RA": agree,
        "gauge": {"mean": float(r_gauge.mean()), "median": float(np.median(r_gauge)),
                  "sd": float(r_gauge.std()),
                  "q": [float(np.percentile(r_gauge, q)) for q in (2.5, 25, 50, 75, 97.5)],
                  "min": float(r_gauge.min()), "max": float(r_gauge.max()),
                  "pct_of_identity_strict": float((r_gauge < ident_r).mean()),
                  "pct_of_identity_midrank": pct_mid,
                  "degenerate_orbit": bool(r_gauge.std() < 1e-9),
                  "n_distinct_configs": int(np.unique(idxs).size),
                  "samples": [float(x) for x in r_gauge]},
        "matched_random": {"mean": float(r_rand.mean()),
                           "median": float(np.median(r_rand)),
                           "sd": float(r_rand.std())},
    }

    # ---- EXACTNESS CHECK: the weight-1 Walsh coefficients from marginals vs from the FWHT
    if ins.N <= (1 << 20):
        c = T.fwht(E)
        nq = ins.n_qubits
        from_fwht = np.array([c[1 << (nq - 1 - q)] for q in range(nq)])
        from_marg = []
        codes = np.arange(k)
        for i in range(n):
            for b in range(bits):
                # qubit index within residue: bit b of the code, MSB first
                bitpos = bits - 1 - b
                chi = 1.0 - 2.0 * ((codes >> bitpos) & 1)
                from_marg.append(float((m[i] * chi).mean()))
        from_marg = np.array(from_marg)
        out["exactness_weight1_coeff_max_abs_diff"] = float(
            np.abs(from_fwht - from_marg).max())
        del c
    return out


def run(targets=None, uniformised=True):
    targets = list(targets or QP.TARGETS19)
    tag = "gauge" if uniformised else "gauge_raw"
    have = A.load(tag)
    for p in targets:
        if p in have:
            print(f"[skip] {p}", flush=True)
            continue
        t0 = time.time()
        r = gauge_cell(p, uniformised=uniformised)
        A.ck(tag, p, r)
        g = r["gauge"]
        print(f"[{p}] {time.time()-t0:5.1f}s  W1(identity) {r['W1_identity_argmin_rmsd_ORACLE']:.3f}"
              f"  RA {r['RA_argmin_rmsd_ORACLE']:.3f}  full {r['full_argmin_rmsd_ORACLE']:.3f}"
              f" | gauge mean {g['mean']:.3f} med {g['median']:.3f}"
              f" 95%[{g['q'][0]:.3f},{g['q'][4]:.3f}]"
              f"  identity at midpct {g['pct_of_identity_midrank']:.3f}"
              f"  agree {r['frac_residue_labellings_where_W1_agrees_with_RA']:.3f}"
              f"  coeffchk {r.get('exactness_weight1_coeff_max_abs_diff', float('nan')):.2e}",
              flush=True)
    A.ck(tag, "_complete", len(A.load(tag)) >= len(QP.TARGETS19))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    rest = sys.argv[2:]
    if mode == "run":
        run(rest or None, True)
    elif mode == "runraw":
        run(rest or None, False)
