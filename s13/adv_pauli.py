"""SPRINT 13 -- ADVERSARIAL AUDIT 6: the Pauli-spectrum -> gradient-variance chain.

Three attacks on the sprint's scientific centrepiece.

  Q1  THE CROSS-VALIDATION GAP.  `s13/walsh_xval.py` validates (a) the TRANSFORM, by running
      walsh's FWHT over geo's OWN cached tables, and (b) the ENERGY TABLE for **legacy only**,
      on 4 cells.  **No AMBER energy table was ever rebuilt independently**, yet every AMBER
      headline rests on those tables.  This rebuilds AMBER tables from scratch through a
      different code path (my own enumeration, my own state->bits mapping, `core.amber` called
      directly) and compares elementwise.  It also checks the INDEXING: table entry x must be
      the energy of the state that MY decoder reads out of x, not of some permutation.

  Q2  THE DROPPED CROSS-COVARIANCES.  `Var_pred = sum_S c_S^2 Var[d<Z_S>/dth]` drops
      `sum_{S!=T} c_S c_T Cov(d<Z_S>, d<Z_T>)`.  The reported ratios (0.993 / 0.915) say the
      dropped term is small FOR THESE HAMILTONIANS.  Constructed adversary: randomise the
      SIGNS of the Walsh coefficients.  The spectrum `|c_S|^2`, and therefore the prediction,
      is IDENTICAL bit for bit, while the true gradient variance is free to move.  If it
      moves a lot, the prediction is a fact about these Hamiltonians and not a theorem.

  Q3  THE ACHIEVABLE RANGE.  For fixed `sum c_S^2`, the true variance
      `c^T Cov c` is maximised/minimised by the extreme eigenvectors of the covariance of
      `d<Z_S>/dth_i` over theta.  Those eigenvalues bound the ratio and say exactly how
      special the real force fields are.

  Q4  EXACTNESS: Parseval and round-trip reconstruction recomputed here, plus a check that
      `geo_common.clean` (which replaces non-finite AMBER entries with the finite maximum)
      never fired on a cell that carries a headline.

    python -m s13.adv_pauli
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import torsion_lib2 as tl2                          # noqa: E402
from s12 import instrument as I                     # noqa: E402
from s13 import geo_common as G                     # noqa: E402
from s13 import geo_pauli as GP                     # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
OUT = {}


# --------------------------------------------------------------------------- Q1
def q1_amber_table(pdb, L, k, max_states=None):
    """Rebuild an AMBER table through an independent path and compare elementwise."""
    import core.amber as AM
    import time
    path = os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")
    if not os.path.exists(path):
        return None
    theirs = np.load(path)
    t = {x["pdb"]: x for x in I.targets()}[pdb]
    seq = t["seq"][:L]
    tab = tl2.library_for(seq, k, t["seq"])          # geo's convention: FULL seq held out
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    b = int(round(np.log2(k)))
    N = k ** L
    idxs = np.arange(N) if max_states is None or max_states >= N else \
        np.random.default_rng(3).choice(N, max_states, replace=False)
    mine = np.empty(len(idxs))
    for _ in range(240):
        if AM.memory_percent() < AM.MEMORY_LIMIT_PERCENT - 2.0:
            break
        time.sleep(10.0)
    for q, j in enumerate(idxs):
        # INDEPENDENT decode: MSB-first per-residue base-k digits
        st = np.array([(int(j) >> (b * (L - 1 - i))) & (k - 1) for i in range(L)], int)
        mine[q] = float(AM.single_point(seq, rep, st, memo=False)["energy"]) \
            if "memo" in AM.single_point.__code__.co_varnames \
            else float(AM.single_point(seq, rep, st)["energy"])
    th = theirs[idxs]
    fin = np.isfinite(th) & np.isfinite(mine)
    rel = np.abs(mine - th) / np.maximum(np.abs(th), 1e-300)
    return {"pdb": pdb, "L": L, "k": k, "n_compared": int(len(idxs)),
            "n_finite_both": int(fin.sum()),
            "max_abs_diff": float(np.abs(mine[fin] - th[fin]).max()),
            "max_rel_diff": float(rel[fin].max()),
            "median_rel_diff": float(np.median(rel[fin])),
            "spearman_mine_vs_theirs": float(np.corrcoef(
                np.argsort(np.argsort(mine[fin])), np.argsort(np.argsort(th[fin])))[0, 1]),
            "n_nonfinite_in_cached": int((~np.isfinite(theirs)).sum())}


# --------------------------------------------------------------------------- Q2/Q3
def q2_q3_cross_terms(pdb, L, k, layers=3, n_theta=128, n_rand=24, seed=11):
    n = int(round(np.log2(k))) * L
    N = 1 << n
    circ = G.circuit(n, layers=layers)
    P = circ.n_params()
    w = GP.popcount(n)
    rng = np.random.default_rng(seed)
    TH = rng.uniform(-np.pi, np.pi, size=(n_theta, P))

    # d<Z_S>/dtheta_i for every S, every i, every theta -- exactly as geo does it
    JH = np.zeros((n_theta, P, N))
    JJ = np.zeros((n_theta, P, N))
    for t in range(n_theta):
        PR = circ.probs_batch(circ._shift_grid(TH[t], np.pi / 2))
        J = (PR[0::2] - PR[1::2]) / 2.0
        JJ[t] = J
        JH[t] = GP.fwht(J)
    VarS = JH.var(0)                              # (P, N)
    out = {"pdb": pdb, "L": L, "k": k, "n_qubits": n, "layers": layers,
           "n_theta": n_theta, "P": P, "variants": {}}

    # Q3: the achievable range of Var_true / Var_pred over ALL unit-norm coefficient vectors
    # Var_true(c) = mean_i c^T Cov_i c ;  Var_pred(c) = sum_S c_S^2 mean_i Var_i[S]
    # -> ratio range = extreme generalised eigenvalues of (mean_i Cov_i, diag(mean_i Var_i))
    C = np.zeros((N, N))
    for i in range(P):
        X = JH[:, i, :] - JH[:, i, :].mean(0)
        C += X.T @ X / (n_theta - 1)
    C /= P
    dg = np.diag(C).copy()                        # == VarS.mean(0)
    keep = dg > dg.max() * 1e-12
    Cs = C[np.ix_(keep, keep)] / np.sqrt(np.outer(dg[keep], dg[keep]))
    ev = np.linalg.eigvalsh(Cs)
    out["Q3_ratio_range"] = {"min_achievable_ratio": float(ev.min()),
                             "max_achievable_ratio": float(ev.max()),
                             "n_strings_kept": int(keep.sum())}

    for v in G.VARIANTS:
        E = G.variant_table(pdb, L, k, v)
        c = GP.fwht(E) / N
        c2 = c ** 2
        # exactness
        var = float(np.var(E)); tot = float(c2[1:].sum())
        recon = GP.fwht(c)                        # inverse (up to no scaling: fwht(fwht(f))=N f)
        rt = float(np.abs(recon - E).max())
        pred = float((c2[None, :] * VarS).sum(1).mean())
        meas = float(np.mean([np.var([JJ[t, i] @ E for t in range(n_theta)])
                              for i in range(P)]))
        rows = {"parseval_rel_err": abs(var - tot) / max(abs(var), 1e-300),
                "reconstruction_max_abs_err": rt,
                "var_grad_pred_exact_per_string": pred,
                "var_grad_measured": meas,
                "ratio_meas_over_pred": meas / max(pred, 1e-300)}
        # Q2: sign randomisation -- identical |c_S|, identical prediction
        ratios = []
        for r in range(n_rand):
            g = np.random.default_rng(1000 + r)
            sg = g.choice([-1.0, 1.0], size=N); sg[0] = 1.0
            c_r = c * sg
            E_r = GP.fwht(c_r)                     # back to configuration space
            m_r = float(np.mean([np.var([JJ[t, i] @ E_r for t in range(n_theta)])
                                 for i in range(P)]))
            ratios.append(m_r / max(pred, 1e-300))
        rows["Q2_sign_randomised_ratio"] = {
            "n": n_rand, "mean": float(np.mean(ratios)), "sd": float(np.std(ratios)),
            "min": float(np.min(ratios)), "max": float(np.max(ratios)),
            "p05": float(np.percentile(ratios, 5)), "p95": float(np.percentile(ratios, 95))}
        out["variants"][v] = rows
    return out


# --------------------------------------------------------------------------- Q4
def q4_clean_fired():
    rows = []
    for p in sorted(glob.glob(os.path.join(G.CACHE, "geo_E_*.npy"))):
        m = re.match(r"geo_E_(\w+)_(\d+)_(\d+)_(\w+)\.npy", os.path.basename(p))
        if not m:
            continue
        E = np.load(p)
        nf = int((~np.isfinite(E)).sum())
        rows.append({"cell": os.path.basename(p), "n": int(E.size), "n_nonfinite": nf,
                     "frac_nonfinite": nf / E.size,
                     "n_distinct": int(len(np.unique(E))),
                     "complete_power_of_k": bool(E.size == int(m.group(3)) ** int(m.group(2)))})
    return rows


def main():
    print("Q4 table completeness / clean() ...")
    OUT["Q4_tables"] = q4_clean_fired()
    bad = [r for r in OUT["Q4_tables"] if r["n_nonfinite"] > 0]
    inc = [r for r in OUT["Q4_tables"] if not r["complete_power_of_k"]]
    print(f"   {len(OUT['Q4_tables'])} cached tables; {len(bad)} with non-finite entries; "
          f"{len(inc)} with the wrong length")
    for r in bad[:8]:
        print("     ", r["cell"], r["n_nonfinite"], "/", r["n"])

    print("Q2/Q3 cross-covariance adversary ...")
    OUT["Q2Q3"] = []
    for (pdb, L, k) in (("1A13", 4, 4), ("2BFI", 5, 4), ("1A1P", 4, 8)):
        if not os.path.exists(os.path.join(G.CACHE, G._tag(pdb, L, k, "amber") + ".npy")):
            continue
        r = q2_q3_cross_terms(pdb, L, k)
        OUT["Q2Q3"].append(r)
        print(f"  {pdb} L={L} k={k} n={r['n_qubits']}  achievable ratio range "
              f"[{r['Q3_ratio_range']['min_achievable_ratio']:.4f}, "
              f"{r['Q3_ratio_range']['max_achievable_ratio']:.4f}]")
        for v, x in r["variants"].items():
            q = x["Q2_sign_randomised_ratio"]
            print(f"     {v:12s} true ratio {x['ratio_meas_over_pred']:.4f} | "
                  f"sign-randomised {q['mean']:.4f} +- {q['sd']:.4f} "
                  f"[{q['min']:.4f},{q['max']:.4f}] | parseval {x['parseval_rel_err']:.1e}")
    json.dump(OUT, open(os.path.join(RESULTS, "adv_pauli.json"), "w"), indent=1)

    print("Q1 independent AMBER table rebuild ...")
    OUT["Q1_amber"] = []
    for (pdb, L, k, mx) in (("1A13", 4, 4, 256), ("2BFI", 4, 4, 256), ("1A1P", 4, 8, 400)):
        try:
            r = q1_amber_table(pdb, L, k, mx)
        except Exception as e:                                      # noqa: BLE001
            r = {"pdb": pdb, "L": L, "k": k, "error": repr(e)}
        if r:
            OUT["Q1_amber"].append(r)
            print("   ", r)
        json.dump(OUT, open(os.path.join(RESULTS, "adv_pauli.json"), "w"), indent=1)
    print("wrote adv_pauli.json")


if __name__ == "__main__":
    main()
