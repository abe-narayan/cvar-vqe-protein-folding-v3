"""ASM-1: the RIGID-PLACEMENT oracle assembly ladder + the retrieval/random nulls.

Rigid placement = every piece is superposed onto its own native segment by its own optimal
rigid transform.  It is a pure upper bound and is NOT one chain: k pieces carry 6k free
rigid parameters.  Because the pieces are placed independently, the total is exactly

    RMSD^2 = (1/n) * sum_pieces SSD(piece, native segment)

so the optimum over cut positions is an exact dynamic program over segment costs.

Candidate sets measured (identical DP, different candidate pool per segment):
    full     -- every window of the legal library of that length (coverage ceiling)
    blos<N>  -- top-N by BLOSUM62 sum vs the target's SUB-SEQUENCE (deployable retrieval)
    rand<N>  -- N uniform-random windows of the same length (the capacity null)

Usage:  python -m s12.asm_rigid
"""
from __future__ import annotations
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I, asm_lib as A          # noqa: E402

LMIN = 3
NRAND = 3            # random-null repeats
ARMS = [("full", None), ("blos500", 500), ("blos100", 100), ("rand500", 500), ("rand100", 100)]


def segment_costs(verbose=True):
    """cost[arm][pdb] = (n, n) array; cost[a, L-1] = SSD of the best length-L piece placed
    on native residues a..a+L-1 (np.inf where the segment does not exist)."""
    tg = I.targets()
    by_fold = {}
    for t in tg:
        by_fold.setdefault(t["fold"], []).append(t)
    nat = {}
    for t in tg:
        u = I.load_univ(t["pdb"]); nat[t["pdb"]] = u["nat_ca"]; del u
    cost = {a: {t["pdb"]: np.full((t["n"], t["n"]), np.inf) for t in tg} for a, _ in ARMS}
    best_idx = {t["pdb"]: {} for t in tg}          # (a,L) -> full-bank argmin index
    rng = np.random.default_rng(20260905)
    t0 = time.time()
    for fold, ts in sorted(by_fold.items()):
        codes = {t["pdb"]: A.encode(t["seq"]) for t in ts}
        for L in range(LMIN, 17):
            if not any(t["n"] >= L for t in ts):
                continue
            b = A.bank(fold, L)
            B62 = A.B62()
            nw = len(b["Wc"])
            for t in ts:
                n = t["n"]
                if n < L:
                    continue
                segs = np.stack([nat[t["pdb"]][a:a + L] for a in range(n - L + 1)])
                M = A.msd_many(b["Wc"], b["n2"], segs) * L          # (nw, nseg) SSD
                for a in range(n - L + 1):
                    cost["full"][t["pdb"]][a, L - 1] = M[:, a].min()
                    best_idx[t["pdb"]][(a, L)] = int(np.argmin(M[:, a]))
                for a in range(n - L + 1):
                    sim = B62[b["S"], codes[t["pdb"]][a:a + L][None, :]].sum(1)
                    ordr = np.argsort(-sim, kind="stable")
                    for arm, N in ARMS:
                        if arm == "full":
                            continue
                        if arm.startswith("blos"):
                            sel = ordr[:N]
                            cost[arm][t["pdb"]][a, L - 1] = M[sel, a].min()
                        else:
                            vals = [M[rng.choice(nw, size=min(N, nw), replace=False), a].min()
                                    for _ in range(NRAND)]
                            cost[arm][t["pdb"]][a, L - 1] = float(np.mean(vals))
                del M
            A.drop_bank(fold, L)
        if verbose:
            print(f"  fold {fold} done {time.time()-t0:.0f}s free {I.free_gb():.2f}", flush=True)
    return cost, best_idx


def dp_ladder(cost, n, kmax=8, lmin=4):
    """best[k] = min total SSD using exactly k contiguous pieces, each of length >= lmin."""
    INF = np.inf
    f = np.full((kmax + 1, n + 1), INF)
    f[0, 0] = 0.0
    for k in range(1, kmax + 1):
        for j in range(1, n + 1):
            best = INF
            for i in range(0, j - lmin + 1):
                L = j - i
                if L < lmin or L > n:
                    continue
                c = cost[i, L - 1]
                if np.isfinite(c) and f[k - 1, i] + c < best:
                    best = f[k - 1, i] + c
            f[k, j] = best
    run = np.inf
    out = {}
    for k in range(1, kmax + 1):          # "at most k pieces" -> always feasible, monotone
        run = min(run, f[k, n])
        out[k] = float(np.sqrt(run / n)) if np.isfinite(run) else None
    return out


def main():
    tg = I.targets()
    cost, best_idx = segment_costs()
    out = {"arms": {}, "per_target": {}}
    for lmin in (4, 3):
        for arm, _ in ARMS:
            lad = {k: [] for k in range(1, 9)}
            for t in tg:
                d = dp_ladder(cost[arm][t["pdb"]], t["n"], kmax=8, lmin=lmin)
                for k in range(1, 9):
                    lad[k].append(d[k])
            key = f"{arm}_lmin{lmin}"
            out["arms"][key] = {}
            for k in range(1, 9):
                v = [x for x in lad[k] if x is not None]
                nfeas = len(v)
                if nfeas == 0:
                    continue
                full = [x if x is not None else np.nan for x in lad[k]]
                out["arms"][key][k] = {
                    "n_feasible": nfeas, "mean_feasible": float(np.mean(v)),
                    "mean_all_maxfill": float(np.nanmean(np.where(np.isnan(full),
                                                                  np.nan, full))),
                    "fail18": float(np.nanmean([full[i] for i, t in enumerate(tg)
                                                if t["pdb"] in I.FAIL18])),
                    "other108": float(np.nanmean([full[i] for i, t in enumerate(tg)
                                                  if t["pdb"] not in I.FAIL18])),
                    "per_target": {t["pdb"]: (None if lad[k][i] is None else float(lad[k][i]))
                                   for i, t in enumerate(tg)}}
    # a monotone "best up to k" ladder on the full arm, lmin=4
    I.write("asm_rigid_ladder", out)
    for lmin in (4, 3):
        print(f"\n=== lmin={lmin}  (rigid placement, ORACLE) ===")
        print("arm        " + "".join(f"  k={k}   " for k in range(1, 7)))
        for arm, _ in ARMS:
            key = f"{arm}_lmin{lmin}"
            row = out["arms"][key]
            s = f"{arm:10s} "
            for k in range(1, 7):
                if k in row:
                    s += f" {row[k]['mean_all_maxfill']:6.3f}({row[k]['n_feasible']:3d})"
                else:
                    s += "     -     "
            print(s)
    return out


if __name__ == "__main__":
    main()
