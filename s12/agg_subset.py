"""EXPERIMENT 4 -- SET SELECTION AS A COMBINATORIAL PROBLEM.

Rather than weighting all 75 candidates, choose a SUBSET S whose coordinate average
minimises a NATIVE-FREE objective.  The value of a candidate depends on which others are
chosen, so the problem is genuinely combinatorial.

------------------------------------------------------------------ EXACT objective F(S)
Candidates A_1..A_m are the shipped top-75, pre-superposed on the SET MEDOID (native-free,
fixed once, so the average is a linear function of the indicator x).

    C(S)      = (1/|S|) sum_{c in S} A_c
    d_p(S)    = || C(S)_i - C(S)_j ||                      (p = (i,j), min_sep 2)
    Risk(S)   = mean_p  risk_p( d_p(S) )                   (the SHIPPED distogram Bayes risk,
                                                            identical to I.shipped_score)
    Cons(S)   = mean_{c != c' in S} P_{c c'}               (internal CA-RMSD spread)
    F(S)      = Risk(S) + mu * Cons(S)

------------------------------------------------------- QUADRATIC surrogate (QUBO / Ising)
For a quantum arm the objective must be quadratic in x.  Replace the risk of the average by
the sd-weighted squared deviation of the AVERAGE OF DISTANCES from the distogram mean
(exact to first order for a tight set), at fixed cardinality k = |S|:

    E(S) = (1/npairs) sum_p w_p [ (1/k) sum_{c in S} d_{p,c} - e_p ]^2 ,   w_p = 1/sd_p^2
         = (1/k^2) x^T Q x  -  (2/k) b^T x  +  const
    Q_{c c'} = (1/npairs) sum_p w_p d_{p,c} d_{p,c'}
    b_c      = (1/npairs) sum_p w_p e_p d_{p,c}
    const    = (1/npairs) sum_p w_p e_p^2

    H(x) = (1/k^2) x^T (Q + mu*P) x - (2/k) b^T x + const + A (sum_c x_c - k)^2

with P the candidate-vs-candidate CA-RMSD matrix (zero diagonal).  Written as an Ising
model in the dumped instances.  16-candidate instances are solved EXACTLY by enumerating
all 2^16 subsets, so a quantum arm can be compared honestly.

Usage:
    python -m s12.agg_subset exact            # greedy / local search / SA on F(S), all 126
    python -m s12.agg_subset instances        # dump 16- and 32-candidate QUBO instances
"""
from __future__ import annotations
import os, sys, json, time, itertools
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_features as AF

MU_DEFAULT = 0.0


class Inst:
    """One target's subset-selection instance."""

    def __init__(self, pdb, m=None):
        o = AF.load(pdb)
        self.pdb = pdb
        self.A = o["Asup"].astype(float)                       # (m, n, 3) medoid-superposed
        self.nat = o["nat"].astype(float)
        self.rr = o["rr"].astype(float)
        self.sc = o["sc"].astype(float)
        self.P = o["P"].astype(float)
        self.pi, self.pj = o["pi"].astype(int), o["pj"].astype(int)
        dg = I.distogram(pdb)
        self.e = dg["expected"].astype(float)
        self.sd = np.maximum(dg["sd"].astype(float), 1e-3)
        self.grid = dg["grid"].astype(float)
        self.risk = dg["risk"].astype(float)
        self.m = len(self.A)
        self.npairs = len(self.pi)
        self.Dc = I.pair_dists(self.A, self.pi, self.pj)       # (m, npairs) candidate dists
        w = 1.0 / self.sd ** 2
        self.w = w
        self.Q = (self.Dc * w) @ self.Dc.T / self.npairs
        self.b = (self.Dc @ (w * self.e)) / self.npairs
        self.const = float((w * self.e ** 2).sum() / self.npairs)

    # -------------------------------------------------- exact objective
    def C(self, S):
        return self.A[list(S)].mean(0)

    def risk_of(self, C):
        d = np.linalg.norm(C[self.pi] - C[self.pj], axis=-1)
        g = np.clip(((d - self.grid[0]) / 0.05).astype(np.int32), 0, len(self.grid) - 1)
        return float(self.risk[np.arange(self.npairs), g].mean())

    def cons(self, S):
        S = list(S)
        if len(S) < 2:
            return 0.0
        sub = self.P[np.ix_(S, S)]
        return float(sub.sum() / (len(S) * (len(S) - 1)))

    def risk_dist(self, S):
        """Risk of the AVERAGE PAIR DISTANCE over S -- cardinality-neutral (the coordinate
        average shrinks distances, which biases `risk_of(C(S))` toward tiny subsets)."""
        d = self.Dc[list(S)].mean(0)
        g = np.clip(((d - self.grid[0]) / 0.05).astype(np.int32), 0, len(self.grid) - 1)
        return float(self.risk[np.arange(self.npairs), g].mean())

    def F(self, S, mu=MU_DEFAULT, mode="dist"):
        r = self.risk_dist(S) if mode == "dist" else self.risk_of(self.C(S))
        return r + mu * self.cons(S)

    def rmsd(self, S):
        return I.ca_rmsd(self.C(S), self.nat)

    # -------------------------------------------------- quadratic surrogate
    def H(self, x, k=None, mu=MU_DEFAULT):
        x = np.asarray(x, float)
        k = k or max(x.sum(), 1)
        M = self.Q + mu * self.P
        return float(x @ M @ x / k ** 2 - 2.0 * (self.b @ x) / k + self.const)

    # -------------------------------------------------- fast incremental exact search
    def _cum(self, S):
        return self.A[list(S)].sum(0)

    def greedy(self, mu=MU_DEFAULT, kmax=None, start=None, mode="dist", kmin=3):
        """Forward greedy on F(S) with 1/|S| normalisation; stops when F stops improving."""
        kmax = kmax or self.m
        rest = set(range(self.m))
        S = list(start) if start else []
        for c in S:
            rest.discard(c)
        best = self.F(S, mu, mode) if S else 1e9
        trace = []
        while len(S) < kmax and rest:
            cands = list(rest)
            vals = []
            for c in cands:
                vals.append(self.F(S + [c], mu, mode))
            j = int(np.argmin(vals))
            if vals[j] >= best - 1e-9 and len(S) >= kmin:
                break
            best = vals[j]
            S.append(cands[j]); rest.discard(cands[j])
            trace.append((len(S), best))
        return S, best, trace

    def greedy_fixed_k(self, k, mu=MU_DEFAULT, mode="dist"):
        rest = set(range(self.m)); S = []
        for _ in range(k):
            cands = list(rest)
            vals = [self.F(S + [c], mu, mode) for c in cands]
            j = int(np.argmin(vals))
            S.append(cands[j]); rest.discard(cands[j])
        return S, self.F(S, mu, mode)

    def local_search(self, S, mu=MU_DEFAULT, iters=40, mode="dist"):
        """Swap / add / drop 1-opt until no improvement."""
        S = list(S); best = self.F(S, mu, mode)
        for _ in range(iters):
            improved = False
            rest = [c for c in range(self.m) if c not in S]
            # drop
            for c in list(S):
                T = [x for x in S if x != c]
                if len(T) >= 3:
                    v = self.F(T, mu, mode)
                    if v < best - 1e-9:
                        S, best, improved = T, v, True
                        break
            if improved:
                continue
            # add
            for c in rest:
                v = self.F(S + [c], mu, mode)
                if v < best - 1e-9:
                    S, best, improved = S + [c], v, True
                    break
            if improved:
                continue
            # swap
            for c in list(S):
                for c2 in rest:
                    T = [x for x in S if x != c] + [c2]
                    v = self.F(T, mu, mode)
                    if v < best - 1e-9:
                        S, best, improved = T, v, True
                        break
                if improved:
                    break
            if not improved:
                break
        return S, best

    def anneal(self, mu=MU_DEFAULT, steps=4000, seed=0, T0=None, S0=None, mode="dist"):
        rng = np.random.default_rng(seed)
        S = set(S0) if S0 is not None else set(rng.choice(self.m, self.m // 2, replace=False).tolist())
        cur = self.F(S, mu, mode); best, bestS = cur, set(S)
        T0 = T0 or max(abs(cur) * 0.05, 1e-3)
        for t in range(steps):
            T = T0 * (1 - t / steps) + 1e-6
            c = int(rng.integers(self.m))
            if c in S:
                if len(S) <= 3:
                    continue
                S.discard(c); v = self.F(S, mu, mode)
                if v < cur or rng.random() < np.exp(-(v - cur) / T):
                    cur = v
                else:
                    S.add(c)
            else:
                S.add(c); v = self.F(S, mu, mode)
                if v < cur or rng.random() < np.exp(-(v - cur) / T):
                    cur = v
                else:
                    S.discard(c)
            if cur < best:
                best, bestS = cur, set(S)
        return sorted(bestS), best

    def oracle_greedy(self, kmax=None):
        """ORACLE ceiling for the uniform-average subset problem: greedy on true RMSD."""
        kmax = kmax or self.m
        rest = set(range(self.m)); S = []; best = 1e9
        while len(S) < kmax and rest:
            cands = list(rest)
            vals = [self.rmsd(S + [c]) for c in cands]
            j = int(np.argmin(vals))
            if vals[j] >= best - 1e-9 and len(S) >= 3:
                break
            best = vals[j]; S.append(cands[j]); rest.discard(cands[j])
        return S, best


# ------------------------------------------------------------------ drivers
def run_exact(out="agg_subset_exact"):
    """Classical solution of the subset problem on all 126 targets.

    mode='dist'  : F(S) = risk( mean_c d_c ) + mu*Cons(S)   -- cardinality-neutral
    mode='struct': F(S) = risk( d(mean_c A_c) ) + mu*Cons(S)
    """
    tg = I.targets()
    rows = {}
    t0 = time.time()
    for kk, t in enumerate(tg):
        ins = Inst(t["pdb"])
        full = list(range(ins.m))
        r = {"avg75": ins.rmsd(full), "argmin": ins.rmsd([int(np.argmin(ins.sc))]),
             "F_dist_avg75": ins.F(full, 0.0, "dist"), "F_struct_avg75": ins.F(full, 0.0, "struct")}
        for mode in ("dist", "struct"):
            for mu in (0.0, 0.01, 0.03):
                tag = f"{mode}_mu{mu}"
                S, f, _ = ins.greedy(mu=mu, mode=mode)
                S2, f2 = ins.local_search(S, mu=mu, iters=12, mode=mode)
                S3, f3 = ins.anneal(mu=mu, steps=3000, seed=0, mode=mode)
                S4, f4 = ins.anneal(mu=mu, steps=3000, seed=1, mode=mode)
                cands = [(f2, S2), (f3, S3), (f4, S4)]
                fb, Sb = min(cands, key=lambda z: z[0])
                r[f"greedy_{tag}"] = ins.rmsd(S); r[f"greedy_{tag}_k"] = len(S)
                r[f"gls_{tag}"] = ins.rmsd(S2)
                r[f"sa_{tag}"] = ins.rmsd(S3)
                r[f"best_{tag}"] = ins.rmsd(Sb); r[f"best_{tag}_k"] = len(Sb); r[f"best_{tag}_F"] = fb
                r[f"gap_{tag}"] = float(f2 - fb)
        for k in (10, 25, 40):
            S, f = ins.greedy_fixed_k(k, mu=0.0, mode="dist")
            r[f"greedyk{k}"] = ins.rmsd(S)
            S, f = ins.greedy_fixed_k(k, mu=0.03, mode="dist")
            r[f"greedyk{k}_mu"] = ins.rmsd(S)
        So, bo = ins.oracle_greedy()
        r["ORC_greedy"] = bo; r["ORC_greedy_k"] = len(So)
        # ORACLE fixed-k uniform-average subsets, for the ceiling of the SELECTION problem
        rows[t["pdb"]] = r
        if kk % 10 == 0:
            print(kk, t["pdb"], f"{time.time()-t0:.0f}s", flush=True)
    I.write(out, rows)
    print("done", time.time() - t0)


def _hypo_set(ins, mode, nh):
    """Choose the nh-candidate hypothesis set for a small QUBO instance."""
    if mode == "score":
        return np.argsort(ins.sc)[:nh]
    # structural: nh medoids of an nh-way agglomerative clustering of the top-75
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    Z = linkage(squareform(ins.P, checks=False), method="average")
    lab = fcluster(Z, t=nh, criterion="maxclust") - 1
    out = []
    for c in range(lab.max() + 1):
        ix = np.where(lab == c)[0]
        out.append(int(ix[np.argmin(ins.P[np.ix_(ix, ix)].mean(1))]))
    return np.array(sorted(out))


def _enum_exact(ins, H, bits, card, ok, mu):
    """Vectorised F(S) = risk_dist(S) + mu*Cons(S) over ALL subsets of the hypothesis set."""
    Dh = ins.Dc[H]                                      # (nh, npairs)
    Dm = (bits @ Dh) / np.maximum(card, 1)[:, None]     # (2^nh, npairs) mean pair distance
    g = np.clip(((Dm - ins.grid[0]) / 0.05).astype(np.int32), 0, len(ins.grid) - 1)
    risk = ins.risk[np.arange(ins.npairs)[None, :], g].mean(1)
    Ph = ins.P[np.ix_(H, H)]
    cons = np.einsum("sc,cd,sd->s", bits, Ph, bits) / np.maximum(card * (card - 1), 1)
    F = risk + mu * cons
    F[~ok] = np.inf
    return F, Dm


def dump_instances(n_targets=12, nh=16, mode="cluster", mu=0.01, out="agg_subset_instances",
                   exhaustive=True):
    tg = I.targets()
    pick = [tg[i] for i in np.linspace(0, len(tg) - 1, n_targets).astype(int)]
    insts = []
    t0 = time.time()
    for t in pick:
        ins = Inst(t["pdb"])
        H = np.asarray(_hypo_set(ins, mode, nh), int)
        Q = ins.Q[np.ix_(H, H)]; Ph = ins.P[np.ix_(H, H)]; b = ins.b[H]
        M = Q + mu * Ph
        rec = {"pdb": t["pdb"], "n": t["n"], "fold": t["fold"], "nh": int(nh), "hypo_mode": mode,
               "mu": mu, "hypo_idx_into_top75": H.tolist(),
               "M": np.round(M, 6).tolist(), "b": np.round(b, 6).tolist(), "const": ins.const,
               "objective": "H(x) = x'Mx/k^2 - 2 b'x/k + const,  k = sum(x) >= 2,  x in {0,1}^nh",
               "avg75_rmsd": float(ins.rmsd(range(ins.m))),
               "uniform_all_hypotheses": {"k": int(nh), "rmsd": float(ins.rmsd(H))},
               "ORACLE_best_single_in_H": float(ins.rr[H].min())}
        # --- heuristic solvers on the hypothesis set (the honest classical baseline)
        def Hval_of(sel):
            x = np.zeros(nh); x[list(sel)] = 1
            k = max(x.sum(), 1)
            return float(x @ M @ x / k ** 2 - 2 * (b @ x) / k + ins.const)
        # greedy on the surrogate
        rest = list(range(nh)); S = []
        bestv = np.inf
        while rest:
            vs = [Hval_of(S + [c]) for c in rest]
            j = int(np.argmin(vs))
            if vs[j] >= bestv - 1e-12 and len(S) >= 2:
                break
            bestv = vs[j]; S.append(rest.pop(j))
        rec["greedy_surrogate"] = {"H": bestv, "S": sorted(S), "k": len(S),
                                   "rmsd": float(ins.rmsd(H[S]))}
        if exhaustive:
            allx = np.arange(1 << nh)
            bits = ((allx[:, None] >> np.arange(nh)) & 1).astype(float)
            card = bits.sum(1); ok = card >= 2
            quad = np.einsum("sc,cd,sd->s", bits, M, bits)
            lin = bits @ b
            Hv = np.full(len(bits), np.inf)
            Hv[ok] = quad[ok] / card[ok] ** 2 - 2 * lin[ok] / card[ok] + ins.const
            j = int(np.argmin(Hv))
            rec["exhaustive_surrogate"] = {"H": float(Hv[j]), "S": np.where(bits[j] > 0)[0].tolist(),
                                           "k": int(card[j]), "rmsd": float(ins.rmsd(H[bits[j] > 0]))}
            Fv, _ = _enum_exact(ins, H, bits, card, ok, mu)
            je = int(np.argmin(Fv))
            rec["exhaustive_exact_F"] = {"F": float(Fv[je]), "S": np.where(bits[je] > 0)[0].tolist(),
                                         "k": int(card[je]), "rmsd": float(ins.rmsd(H[bits[je] > 0]))}
            Cs = np.einsum("sc,cnk->snk", bits[ok], ins.A[H]) / card[ok][:, None, None]
            rm = I.kabsch_rmsd_batch(Cs, ins.nat)
            io = np.where(ok)[0][int(np.argmin(rm))]
            rec["ORACLE_best_subset"] = {"rmsd": float(rm.min()), "k": int(card[io]),
                                         "S": np.where(bits[io] > 0)[0].tolist()}
            # how well does the native-free objective rank subsets?
            rec["rank_corr_surrogate_vs_rmsd"] = float(np.corrcoef(Hv[ok], rm)[0, 1])
            rec["rank_corr_exactF_vs_rmsd"] = float(np.corrcoef(Fv[ok], rm)[0, 1])
            rec["rmsd_percentile_of_surrogate_opt"] = float((rm < ins.rmsd(H[bits[j] > 0])).mean())
            rec["greedy_gap_to_exhaustive"] = float(rec["greedy_surrogate"]["H"] - Hv[j])
        insts.append(rec)
        print(t["pdb"], f"{time.time()-t0:.0f}s", flush=True)
    I.write(out, {"instances": insts,
                  "note": "x in {0,1}^nh selects candidates from hypo_idx_into_top75 (indices into "
                          "the shipped top-75); the emitted structure is the uniform coordinate "
                          "average of the selected medoid-superposed candidates."})
    def mn(f):
        return float(np.mean([f(x) for x in insts]))
    print("avg75 %.4f | uniform-H %.4f | greedy %.4f" % (
        mn(lambda x: x["avg75_rmsd"]), mn(lambda x: x["uniform_all_hypotheses"]["rmsd"]),
        mn(lambda x: x["greedy_surrogate"]["rmsd"])), end="")
    if exhaustive:
        print(" | QUBO-opt %.4f | exactF-opt %.4f | ORACLE %.4f | rho(H,rmsd) %.3f" % (
            mn(lambda x: x["exhaustive_surrogate"]["rmsd"]), mn(lambda x: x["exhaustive_exact_F"]["rmsd"]),
            mn(lambda x: x["ORACLE_best_subset"]["rmsd"]), mn(lambda x: x["rank_corr_surrogate_vs_rmsd"])))
    else:
        print()


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "exact":
        run_exact()
    elif cmd == "instances":
        nh = int(sys.argv[2]) if len(sys.argv) > 2 else 16
        mode = sys.argv[3] if len(sys.argv) > 3 else "cluster"
        nt = int(sys.argv[4]) if len(sys.argv) > 4 else 12
        tag = "agg_subset_instances" if (nh == 16 and mode == "cluster") else f"agg_subset_instances_{nh}_{mode}"
        dump_instances(nh=nh, mode=mode, n_targets=nt, out=tag, exhaustive=(nh <= 20))
