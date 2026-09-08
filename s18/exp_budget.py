"""s18/exp_budget.py -- CONTROL B: is any of this a SEARCH problem?

Sprint 17 found greedy 1-opt certifies the FULL objective's optimum at 1,024 evaluations on the
enumerated instrument, so budget was not binding there.  The BRIEF asks whether that still holds
for degree-1.  This module answers it on the real 126-target instrument.

THE DISCRETE INSTRUMENT.  Greedy 1-opt needs a move set.  Native-free and target-conditioned:
each residue gets `K = 8` candidate `(phi, psi)` states, the K medoids of that residue's 75
values in the target's own shipped retrieval pool under the circular metric.  The start state is
the pool state nearest the projection of the coordinate average -- the SAME Control-A start as
every other arm.

ACCOUNTING CONVENTION, stated once and applied to every method.

    ONE EVALUATION = one call of the objective at one COMPLETE configuration.

L-BFGS is also charged for its gradient: a gradient at a complete configuration is reported
BOTH as 1 evaluation (analytic, what it costs here) and as 2n+1 (what it would cost by finite
differences), because the two conventions give opposite readings of "matched budget" and quoting
one without the other is how a search comparison is rigged.  Wall-clock per call is reported
separately, because an `E_le1` call is a table lookup and an `E_full` call is a chain build --
matched CALLS are not matched COST, and the report says which is which.

METHODS, all from the identical start at the identical budget:
    greedy_full   greedy 1-opt on E_full   (best improving single-residue change, to a local opt)
    greedy_le1    greedy 1-opt on E_le1
    rls_full      random local search (random single-residue change, accept if improving)
    rls_le1       ditto on E_le1
    anneal_full   Metropolis with a geometric temperature schedule
    anneal_le1    ditto
    exact_le1     the EXACT coordinate-wise discrete argmin of E_le1 -- costs n*K evaluations
                  and is CERTIFIED by separability, not by search.

Run:  python -m s18.exp_budget
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s18 import exp_obj as XO                # noqa: E402
from s18 import exp_run as XR                # noqa: E402

RESULTS = os.path.join(HERE, "results")
OUT = os.path.join(RESULTS, "exp_budget.json")
K = 8
BUDGET = 1024


def _circ_medoids(P, S, k, rng):
    """`k` medoid `(phi, psi)` states of one residue's pool values, circular metric."""
    X = np.column_stack([P, S])
    m = len(X)
    if m <= k:
        return X.copy()
    D = np.abs(X[:, None, :] - X[None, :, :])
    D = np.minimum(D, 2 * np.pi - D)
    D = np.sqrt((D ** 2).sum(2))
    #: k-medoids by farthest-point seeding then one PAM-style reassignment sweep
    idx = [int(rng.integers(m))]
    while len(idx) < k:
        idx.append(int(np.argmax(D[:, idx].min(1))))
    for _ in range(6):
        lab = np.argmin(D[:, idx], 1)
        idx = [int(np.array(np.where(lab == c)[0])[
            np.argmin(D[np.ix_(np.where(lab == c)[0], np.where(lab == c)[0])].sum(1))])
            if (lab == c).any() else idx[c] for c in range(k)]
    return X[idx]


class Discrete:
    """A `n x K` state instrument with a counted objective."""

    def __init__(self, ob, states):
        self.ob = ob
        self.S = states                      # (n, K, 2)
        self.n = states.shape[0]
        self.K = states.shape[1]
        self.calls = 0

    def torsions(self, z):
        ph = self.S[np.arange(self.n), z, 0]
        ps = self.S[np.arange(self.n), z, 1]
        return ph, ps

    def E(self, z, which):
        ph, ps = self.torsions(z)
        self.calls += 1
        return (self.ob.E_full(ph, ps)[0] if which == "full" else self.ob.E_le1(ph, ps)[0])


def greedy(D, z0, which, budget):
    """Greedy 1-opt: the best improving single-residue change, repeated to a local optimum."""
    z = np.array(z0, int)
    f = D.E(z, which)
    trace = [(D.calls, f)]
    conv = False
    while D.calls < budget:
        best, bz = f, None
        for a in range(D.n):
            for v in range(D.K):
                if v == z[a] or D.calls >= budget:
                    continue
                w = z.copy(); w[a] = v
                g = D.E(w, which)
                if g < best:
                    best, bz = g, w
        if bz is None:
            conv = True            # a certified LOCAL optimum of the 1-opt neighbourhood
            break
        z, f = bz, best
        trace.append((D.calls, f))
    return z, f, trace, conv


def rls(D, z0, which, budget, rng):
    z = np.array(z0, int)
    f = D.E(z, which)
    while D.calls < budget:
        w = z.copy()
        w[rng.integers(D.n)] = rng.integers(D.K)
        g = D.E(w, which)
        if g < f:
            z, f = w, g
    return z, f


def anneal(D, z0, which, budget, rng, t0=None):
    z = np.array(z0, int)
    f = D.E(z, which)
    bz, bf = z.copy(), f
    T0 = abs(f) * 0.05 + 1e-9 if t0 is None else t0
    while D.calls < budget:
        frac = D.calls / budget
        T = T0 * (0.01 ** frac)
        w = z.copy()
        w[rng.integers(D.n)] = rng.integers(D.K)
        g = D.E(w, which)
        if g < f or rng.random() < np.exp(-(g - f) / max(T, 1e-12)):
            z, f = w, g
            if f < bf:
                bz, bf = z.copy(), f
    return bz, bf


def exact_le1(D, z0):
    """Global discrete argmin of `E_le1`.  CERTIFIED by additivity, exactly `n*K` calls.

    `E_le1` is a sum of per-residue terms, so changing residue `a` changes only its own term:
    the per-residue argmin of the FULL objective (evaluated with only that residue moved) IS
    the global optimum of the whole configuration.  Evaluating the objective itself -- rather
    than reading the nearest mesh cell of the tabulated field, which is what an earlier version
    of this function did and which returned a WORSE value than greedy on 9 of 10 targets --
    keeps the arm exact and keeps the call count honest.

    Residues the field does not touch at all (MATH's theorem: 0 and n-1) are held at the START
    state; an argmin over an identically-flat term would charge an arbitrary structural
    commitment to the objective.  Those residues are RMSD-irrelevant in any case.
    """
    n, Kk = D.n, D.K
    z = np.array(z0, int)
    f = D.ob.t.f
    for a in range(n):
        if np.ptp(f[a]) < 1e-12:
            D.calls += Kk
            continue
        vals = []
        for v in range(Kk):
            w = z.copy(); w[a] = v
            vals.append(D.E(w, "le1"))
        z[a] = int(np.argmin(vals))
    return z


def run(targets=None, out=OUT, S=None, cached_only=False):
    from s18 import math_anova as MA
    from s18 import math_lib as ML
    S = S or MA.NSAMP
    tg = targets if targets is not None else I.targets()
    rows = []
    if os.path.exists(out):
        try:
            p = json.load(open(out))
            if not p.get("complete"):
                rows = p["rows"]
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    data = C.gather(tg)
    pdbs = [t["pdb"] for t in tg]
    deb = {}
    for f in sorted({int(t["fold"]) for t in tg}):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))
    t0 = time.time()

    for t in tg:
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        if pdb in done:
            continue
        d = data[pdb]
        i, j, sd, nat = d["i"], d["j"], d["sd"], np.asarray(d["nat"], float)
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        rng = SD.stable_rng(pdb, "s18budget")
        W, PH, PS, _u = AV.top75_windows(pdb)
        W = np.asarray(W, float)
        P = I.pairwise_rmsd(W)
        avg, _b = I.coordinate_average(W, P)
        pr = I.project(np.asarray(avg, float), seq, fold)
        phi0 = np.asarray(pr["phi"], float)
        psi0 = np.asarray(pr["psi"], float)

        if cached_only and not os.path.exists(os.path.join(
                ML.CACHE, f"anova_{pdb}_pool_{S}_{MA.GRID}.npz")):
            continue
        ob, _dd = XR._target_obj(t, deb, "pool", S)
        St = np.stack([_circ_medoids(PH[:, a], PS[:, a], K, rng) for a in range(n)])
        #: the start state: the discrete state nearest the SAME Control-A projection
        z0 = np.empty(n, int)
        for a in range(n):
            dd = np.abs(St[a] - np.array([phi0[a], psi0[a]]))
            dd = np.minimum(dd, 2 * np.pi - dd)
            z0[a] = int(np.argmin((dd ** 2).sum(1)))

        e = {"pdb": pdb, "n": n, "fold": fold, "K": K, "budget": BUDGET,
             "avg": float(I.ca_rmsd(np.asarray(avg, float), nat))}
        D0 = Discrete(ob, St)
        e["start_objfull"] = D0.E(z0, "full")
        e["start_objle1"] = D0.E(z0, "le1")
        e["start_rmsd"] = float(I.ca_rmsd(I.build_ca(*D0.torsions(z0)), nat))

        for which in ("full", "le1"):
            D = Discrete(ob, St); tw = time.time()
            z, f, tr, conv = greedy(D, z0, which, BUDGET)
            e[f"greedy_{which}"] = float(I.ca_rmsd(I.build_ca(*D.torsions(z)), nat))
            e[f"greedy_{which}_obj"] = f
            e[f"greedy_{which}_calls"] = D.calls
            e[f"greedy_{which}_converged"] = bool(conv)
            e[f"greedy_{which}_trace"] = [[int(a), float(b)] for a, b in tr]
            e[f"greedy_{which}_sec"] = time.time() - tw
            #: the same greedy at a 4x budget -- does more budget buy anything?
            D4 = Discrete(ob, St)
            z4, f4, _tr, c4 = greedy(D4, z0, which, 4 * BUDGET)
            e[f"greedy4x_{which}"] = float(I.ca_rmsd(I.build_ca(*D4.torsions(z4)), nat))
            e[f"greedy4x_{which}_obj"] = f4
            e[f"greedy4x_{which}_calls"] = D4.calls

            D = Discrete(ob, St)
            z, f = rls(D, z0, which, BUDGET, rng)
            e[f"rls_{which}"] = float(I.ca_rmsd(I.build_ca(*D.torsions(z)), nat))
            e[f"rls_{which}_obj"] = f

            D = Discrete(ob, St)
            z, f = anneal(D, z0, which, BUDGET, rng)
            e[f"anneal_{which}"] = float(I.ca_rmsd(I.build_ca(*D.torsions(z)), nat))
            e[f"anneal_{which}_obj"] = f

            #: best-of-N random restarts at the SAME total budget -- the control the record
            #: says every concentration claim must carry
            D = Discrete(ob, St)
            bf, bz = np.inf, None
            while D.calls < BUDGET:
                w = rng.integers(0, K, n)
                g = D.E(w, which)
                if g < bf:
                    bf, bz = g, w
            e[f"bestofN_{which}"] = float(I.ca_rmsd(I.build_ca(*D.torsions(bz)), nat))
            e[f"bestofN_{which}_obj"] = float(bf)

        D = Discrete(ob, St)
        z = exact_le1(D, z0)
        e["exact_le1"] = float(I.ca_rmsd(I.build_ca(*D.torsions(z)), nat))
        e["exact_le1_obj"] = D.E(z, "le1")
        e["exact_le1_calls"] = D.calls
        e["exact_le1_objfull"] = D.E(z, "full")

        rows.append(e)
        json.dump({"rows": rows, "complete": False, "K": K, "budget": BUDGET, "S": S}, open(out, "w"))
        if len(rows) % 5 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s) {pdb}", flush=True)

    json.dump({"rows": rows, "complete": len(rows) == len(tg), "K": K, "budget": BUDGET,
               "S": S, "n_expected": len(tg), "n_rows": len(rows)}, open(out, "w"))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


if __name__ == "__main__":
    run(cached_only=(len(sys.argv) > 1 and sys.argv[1] == "cached"))
