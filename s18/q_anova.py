"""SPRINT 18 / QUANTUM-ADVERSARIAL -- A2: audit the degree-1 OBJECT itself, on the lattice,
before anyone ports it to continuous torsion space.

THE BRIEF'S F5 IS "the continuous mapping is mathematically invalid".  The most likely quiet
death of the sprint is not a bad number, it is a NAME COLLISION: three different objects are
all being called "degree-1" and they are not the same object.

    (W1)  the STRICT WALSH weight-<=1 projection            -- the object that produced 2.411 A
    (RA)  the RESIDUE-ADDITIVE first-order functional ANOVA -- the object the brief's section 4
          E0 + sum_i ( E[E | theta_i] - E0 )                   defines, and the only one that
                                                               survives into continuous torsions
    (RAW) either of the above applied to the RAW objective  -- as against the RANK-UNIFORMISED
                                                               objective Sprint 17 actually used

CLAIM UNDER AUDIT (brief section 4): "When mu is uniform on the enumerated lattice this IS the
Walsh weight-<=1 projection -- that equivalence must be verified numerically, not asserted."

IT IS FALSE AS STATED, and this module proves it in both directions:
  * the residue-additive ANOVA object equals the projection onto Walsh coefficients whose
    support lies INSIDE A SINGLE RESIDUE -- weight 0, weight 1, AND intra-residue weight 2.
    With k = 4 states in 2 qubits a per-residue field NEEDS its intra-residue weight-2 term;
    Sprint 17 measured 95.7% of all weight-2 mass as intra-residue, so this is not a rounding
    difference, it is most of the second-order mass.
  * therefore W1 != RA, and 2.411 A is a number about W1, not about the object the brief
    tells the sprint to port.

And a second name collision that changes the same number: the truncation Sprint 17 reported
is of the RANK-UNIFORMISED objective.  Uniformisation is monotone, so it cannot change the
FULL objective's argmin -- but a projection of a monotone transform is NOT a monotone
transform of the projection, so the truncated object's argmin CAN move.  Both are built here.

ALSO MEASURED, because it decides the quantum question: how SEPARABLE each object is.  A
strictly weight-<=1 objective has ZERO couplings and its global optimum is a closed-form
per-qubit argmin; a residue-additive objective's optimum is a closed-form per-residue argmin.
Either way there is no search problem at all, which is strictly WORSE for a quantum case than
the full objective already was.

RUN:  python -m s18.q_anova run [pdb ...]
      python -m s18.q_anova report
"""
from __future__ import annotations

import gc
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

RESULTS = os.path.join(ROOT, "s18", "results")
os.makedirs(RESULTS, exist_ok=True)
SALT = "s18quantum"
TARGETS19 = QP.TARGETS19
BANDF = 0.01


# ------------------------------------------------------------------ checkpoint
def ck(tag, key, value):
    path = os.path.join(RESULTS, f"q_{tag}.json")
    d = {}
    if os.path.exists(path):
        try:
            d = json.load(open(path))
        except Exception:
            d = {}
    d[key] = value
    d["_written"] = time.strftime("%Y-%m-%d %H:%M:%S")
    tmp = path + f".tmp{os.getpid()}"
    json.dump(d, open(tmp, "w"), indent=1,
              default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)


def load(tag):
    path = os.path.join(RESULTS, f"q_{tag}.json")
    return json.load(open(path)) if os.path.exists(path) else {}


# ------------------------------------------------------------------ the objects
def residue_additive(E, n, k, pw):
    """(RA)  E0 + sum_i ( E_mu[E | residue i = state] - E0 ),  mu = uniform on the lattice.

    This is the brief's section-4 first-order functional ANOVA, computed EXACTLY (not by
    Monte Carlo) because the register is enumerated.  `E_mu[E | theta_i = v]` is the mean of
    E over the k^(n-1) configurations with residue i in state v.
    """
    N = E.size
    E0 = float(E.mean())
    out = np.full(N, E0)
    idx = np.arange(N, dtype=np.int64)
    for i in range(n):
        s = (idx // pw[i]) % k
        m = np.bincount(s, weights=E, minlength=k) / np.bincount(s, minlength=k)
        out += (m[s] - E0)
        del s
    return out


def residue_masks(n, bits, nq):
    """Boolean masks over Walsh index space: which coefficients live inside one residue."""
    idx = np.arange(1 << nq, dtype=np.int64)
    deg = T._popcount(idx)
    # qubit q (0 = MSB) has bit value 1 << (nq-1-q); residue i owns qubits [i*bits, i*bits+bits)
    inside = np.zeros(1 << nq, bool)
    for i in range(n):
        m = 0
        for b in range(bits):
            m |= 1 << (nq - 1 - (i * bits + b))
        inside |= (idx & ~m) == 0
    inside[0] = True                        # the constant
    return idx, deg, inside


def objects(ins, uniformised=True):
    """Build every candidate 'degree-1' object for one target and return them with the
    reference quantities needed to score them."""
    nq, n, k = ins.n_qubits, ins.n, ins.k
    bits = ins.bits_per_res
    pw = k ** np.arange(n - 1, -1, -1)
    Eraw = ins.hamil().astype(np.float64)
    E = V.uniformise(Eraw) if uniformised else Eraw
    idx, deg, inside = residue_masks(n, bits, nq)
    c = T.fwht(E)

    o = {}
    o["full"] = E
    cc = c.copy(); cc[deg > 1] = 0.0
    o["W1"] = T.ifwht(cc)                                    # strict Walsh weight <= 1
    cc = c.copy(); cc[~inside] = 0.0
    o["RA_walsh"] = T.ifwht(cc)                              # residue-block projection
    o["RA"] = residue_additive(E, n, k, pw)                  # ANOVA, computed directly
    cc = c.copy(); cc[deg > 2] = 0.0
    o["W2"] = T.ifwht(cc)                                    # strict Walsh weight <= 2
    del c, cc
    gc.collect()
    return o, deg, inside, pw


def argmin_row(g, rmsd, band):
    tied = np.flatnonzero(g == g.min())
    order = np.argsort(g, kind="mergesort")[:band]
    br = rmsd[order]
    am = float(rmsd[tied].mean())
    # where the argmin sits INSIDE its own band -- the min-of-N null's own coordinate
    pct_in_band = float((br < am).mean())
    return {
        "argmin_rmsd_ORACLE": am,
        "argmin_ties": int(tied.size),
        "argmin_pct_global_ORACLE": float((rmsd < am).mean()),
        "band_mean_ORACLE": float(br.mean()),
        "band_best_ORACLE": float(br.min()),
        "band_sd_ORACLE": float(br.std()),
        "argmin_pct_within_band_ORACLE": pct_in_band,
        "rho_global_ORACLE": float(V.spearman(g, rmsd)),
        "rho_inband_ORACLE": float(V.spearman(g[order], br)),
        "band_rmsd_q": [float(np.percentile(br, q)) for q in (0, 1, 5, 25, 50, 75, 100)],
    }


def separability(g, deg, inside, nq):
    """Variance fractions of a field: at weight <= 1, inside-one-residue, and inter-residue."""
    c = T.fwht(g)
    p = c ** 2
    tot = float(p.sum() - p[0])
    if tot <= 0:
        return {"var_w1": float("nan"), "var_inside": float("nan"),
                "var_inter": float("nan"), "mean_weight": float("nan")}
    w1 = float(p[deg == 1].sum()) / tot
    ins_ = (float(p[inside].sum()) - float(p[0])) / tot
    mw = float((p[1:] * deg[1:]).sum()) / tot
    return {"var_w1": w1, "var_inside_one_residue": ins_,
            "var_inter_residue": 1.0 - ins_, "mean_weight": mw}


# ------------------------------------------------------------------ per-target cell
def cell(pdb, uniformised=True):
    ins = QP.inst(pdb)
    rmsd = ins.rmsd
    band = max(64, ins.N // 100)
    o, deg, inside, pw = objects(ins, uniformised)
    row = {"pdb": pdb, "n": ins.n, "k": ins.k, "nq": ins.n_qubits, "N": int(ins.N),
           "fold": ins.fold, "band": int(band), "uniformised": bool(uniformised),
           "space_best_ORACLE": float(rmsd.min()),
           "space_mean_ORACLE": float(rmsd.mean())}

    # ---- THE EQUIVALENCE AUDIT, both directions
    ra, raw_, w1 = o["RA"], o["RA_walsh"], o["W1"]
    sc = float(np.std(o["full"])) or 1.0
    row["audit"] = {
        "max_abs_RA_minus_residue_block_projection": float(np.abs(ra - raw_).max()),
        "relative_to_objective_sd": float(np.abs(ra - raw_).max() / sc),
        "max_abs_RA_minus_W1": float(np.abs(ra - w1).max()),
        "RA_minus_W1_rel_sd": float(np.abs(ra - w1).max() / sc),
        "sd_RA": float(ra.std()), "sd_W1": float(w1.std()),
        "var_of_W1_over_var_of_RA": float(w1.var() / max(ra.var(), 1e-300)),
        "spearman_RA_W1": float(V.spearman(ra, w1)),
        "frac_configs_same_argmin_RA_W1": float(
            (np.argmin(ra) == np.argmin(w1))),
    }

    # ---- score every object
    row["obj"] = {}
    for name in ("full", "W1", "RA", "W2", "RA_walsh"):
        r = argmin_row(o[name], rmsd, band)
        r.update(separability(o[name], deg, inside, ins.n_qubits))
        r["rho_with_full"] = float(V.spearman(o[name], o["full"]))
        row["obj"][name] = r

    # ---- the CLOSED-FORM optimum of each additive object, and its cost
    #      W1: independent per QUBIT.  RA: independent per RESIDUE.
    idxall = np.arange(ins.N, dtype=np.int64)
    # per-residue table for RA
    E = o["full"]
    E0 = float(E.mean())
    best_state = np.zeros(ins.n, int)
    for i in range(ins.n):
        s = (idxall // pw[i]) % ins.k
        m = np.bincount(s, weights=E, minlength=ins.k) / np.bincount(s, minlength=ins.k)
        best_state[i] = int(np.argmin(m))
        del s
    ci = int(best_state @ pw)
    row["closed_form"] = {
        "RA_argmin_index": ci,
        "RA_argmin_rmsd_ORACLE": float(rmsd[ci]),
        "RA_argmin_matches_enumerated": bool(abs(o["RA"][ci] - o["RA"].min()) < 1e-9),
        "RA_cost_objective_evaluations": int(ins.n * ins.k),
        "W1_cost_objective_evaluations": int(ins.n_qubits * 2),
    }
    del o, ra, raw_, w1
    gc.collect()
    return row


def run(targets=None, uniformised=True):
    targets = list(targets or TARGETS19)
    tag = "anova" if uniformised else "anova_raw"
    have = load(tag)
    for p in targets:
        if p in have:
            print(f"[skip] {p}", flush=True)
            continue
        t0 = time.time()
        r = cell(p, uniformised)
        ck(tag, p, r)
        a = r["audit"]
        print(f"[{p}] n={r['n']} {time.time()-t0:5.1f}s  "
              f"|RA - blockproj| = {a['max_abs_RA_minus_residue_block_projection']:.3e}  "
              f"|RA - W1| = {a['max_abs_RA_minus_W1']:.4f} ({a['RA_minus_W1_rel_sd']:.3f} sd)  "
              f"argmin W1 {r['obj']['W1']['argmin_rmsd_ORACLE']:.3f} "
              f"RA {r['obj']['RA']['argmin_rmsd_ORACLE']:.3f} "
              f"full {r['obj']['full']['argmin_rmsd_ORACLE']:.3f}", flush=True)
    ck(tag, "_complete", len(load(tag)) >= len(TARGETS19))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    rest = sys.argv[2:]
    if mode == "run":
        run(rest or None, uniformised=True)
    elif mode == "runraw":
        run(rest or None, uniformised=False)
