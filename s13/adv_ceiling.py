"""SPRINT 13 -- ADVERSARIAL AUDIT 2: is 1.594 A the LIBRARY's information, or the SEARCH's?

C0 claims "the representation is not the barrier": coordinate descent on true CA-RMSD over
`torsion_lib2.library_for(seq, 4, seq)` reaches 1.594 A in ~26 qubits.  The alternative
reading is that a descent making n * k * sweeps ORACLE CA-RMSD queries is simply a powerful
search that would find a good structure in almost any sufficiently expressive space.

THE DECISIVE CONTROL: run the IDENTICAL descent, same start rule, same sweep cap, same k,
in matched spaces that carry strictly less information than the real library:

    real            tl2.library_for(seq, k, seq)                      -- the claim
    real_seed7      same, seed=7                                      -- seed stability
    wrongseq        the library of a DIFFERENT target's sequence      -- sequence relevance
    marg_class      k-means k over the WHOLE corpus, class-blind,
                    the SAME k states for every residue               -- corpus marginal
    marg_draw       k (phi,psi) pairs drawn at random from the corpus -- marginal, unclustered
    unif_torus      k (phi,psi) pairs uniform on [-pi,pi]^2           -- no Ramachandran at all

Plus two attacks on the descent itself:
    desc_rand_start random start instead of the ORACLE nearest-state snap
    desc_long       24 sweeps and 4 random restarts

Every arm reads the native to score, so every arm is an ORACLE DIAGNOSTIC.

    python -m s13.adv_ceiling [k] [n_targets]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402
import representations as reps             # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)


def _wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# ------------------------------------------------------------------ corpus marginal
_POOL = None


def corpus_pool(exclude_seq=""):
    """Every (phi, psi) observed in the held-out peptide database, class-blind."""
    global _POOL
    if _POOL is not None and exclude_seq == "":
        return _POOL
    ent = pdb.holdout(exclude_seq) if exclude_seq else pdb.load()
    P = np.concatenate([np.column_stack([np.asarray(p.phi, float),
                                         np.asarray(p.psi, float)]) for p in ent])
    if exclude_seq == "":
        _POOL = P
    return P


def _kmeans_states(P, k, seed):
    X4 = np.column_stack([np.cos(P[:, 0]), np.sin(P[:, 0]), np.cos(P[:, 1]), np.sin(P[:, 1])])
    C = tl2._circ_kmeans(X4, np.ones(len(X4)), k, seed=seed)
    return tl2._angles_from_centres(C)


# ------------------------------------------------------------------ counted descent
class Counter:
    def __init__(self):
        self.n = 0


def descent(PHI, PSI, nat, states, max_sweeps=12, cnt=None):
    """ceiling.descent, verbatim, with an ORACLE-query counter bolted on."""
    n, k = PHI.shape
    s = np.asarray(states, int).copy()
    idx = np.arange(n)
    build = lambda ss: I.build_ca(PHI[idx, ss], PSI[idx, ss])          # noqa: E731
    best = I.ca_rmsd(build(s), nat)
    if cnt: cnt.n += 1
    sweeps_used = 0
    for _ in range(max_sweeps):
        moved = False
        sweeps_used += 1
        for i in range(n):
            cur = s[i]
            cand = np.repeat(s[None, :], k, axis=0)
            cand[:, i] = np.arange(k)
            W = I.build_ca(PHI[idx[None, :], cand], PSI[idx[None, :], cand])
            r = I.kabsch_rmsd_batch(W, nat)
            if cnt: cnt.n += k
            b = int(np.argmin(r))
            if r[b] < best - 1e-9:
                best = float(r[b]); s[i] = b; moved = True
            else:
                s[i] = cur
        if not moved:
            break
    return s, best, sweeps_used


def snap(PHI, PSI, phi0, psi0):
    d = np.abs(_wrap(PHI - phi0[:, None])) + np.abs(_wrap(PSI - psi0[:, None]))
    return np.argmin(d, axis=1)


# ------------------------------------------------------------------ spaces
def spaces_for(t, tg, k, rng):
    """dict name -> (n, k, 2) torsion table."""
    seq, n = t["seq"], t["n"]
    out = {}
    out["real"] = tl2.library_for(seq, k, seq)
    out["real_seed7"] = tl2.library_for(seq, k, seq, 7)
    # a DIFFERENT target of the same length, if one exists; else any other target padded
    others = [o for o in tg if o["pdb"] != t["pdb"] and o["n"] >= n]
    o = others[int(rng.integers(len(others)))]
    out["wrongseq"] = tl2.library_for(o["seq"], k, o["seq"])[:n]
    P = corpus_pool()
    C = _kmeans_states(P, k, seed=1)                     # class-blind corpus marginal
    out["marg_class"] = np.repeat(C[None], n, axis=0)
    draw = P[rng.integers(0, len(P), size=(n, k))]       # marginal, unclustered
    out["marg_draw"] = draw
    out["unif_torus"] = rng.uniform(-np.pi, np.pi, size=(n, k, 2))
    return out


def run_target(t, tg, k, rng, sweeps=12):
    u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
    p = pdb.by_pdb(t["pdb"])
    phi0 = np.asarray(p.phi, float); psi0 = np.asarray(p.psi, float)
    n = t["n"]
    row = {"pdb": t["pdb"], "n": n, "fold": t["fold"]}
    sp = spaces_for(t, tg, k, rng)
    for name, tab in sp.items():
        PHI = np.ascontiguousarray(tab[:, :, 0]); PSI = np.ascontiguousarray(tab[:, :, 1])
        s0 = snap(PHI, PSI, phi0, psi0)
        idx = np.arange(n)
        r_snap = I.ca_rmsd(I.build_ca(PHI[idx, s0], PSI[idx, s0]), nat)
        c = Counter()
        _, r_desc, sw = descent(PHI, PSI, nat, s0, sweeps, c)
        row[name] = {"snap": float(r_snap), "descent": float(r_desc),
                     "oracle_queries": int(c.n), "sweeps": int(sw)}
    # --- attacks on the descent itself, on the REAL space only
    tab = sp["real"]
    PHI = np.ascontiguousarray(tab[:, :, 0]); PSI = np.ascontiguousarray(tab[:, :, 1])
    c = Counter()
    best_r = np.inf
    for _ in range(4):
        s0 = rng.integers(0, k, n)
        _, r, _ = descent(PHI, PSI, nat, s0, sweeps, c)
        best_r = min(best_r, r)
    row["desc_rand_start_x4"] = {"descent": float(best_r), "oracle_queries": int(c.n)}
    # the same random-start budget inside the INFORMATION-FREE space
    tabu = sp["unif_torus"]
    PU = np.ascontiguousarray(tabu[:, :, 0]); SU = np.ascontiguousarray(tabu[:, :, 1])
    c = Counter(); best_u = np.inf
    for _ in range(4):
        s0 = rng.integers(0, k, n)
        _, r, _ = descent(PU, SU, nat, s0, sweeps, c)
        best_u = min(best_u, r)
    row["unif_rand_start_x4"] = {"descent": float(best_u), "oracle_queries": int(c.n)}
    c = Counter()
    s0 = snap(PHI, PSI, phi0, psi0)
    _, r_long, sw = descent(PHI, PSI, nat, s0, 40, c)
    row["desc_long40"] = {"descent": float(r_long), "oracle_queries": int(c.n), "sweeps": int(sw)}
    return row


def report(rows, k):
    names = [n for n in ("real", "real_seed7", "wrongseq", "marg_class", "marg_draw",
                         "unif_torus") if n in rows[0]]
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    folds = [r["fold"] for r in rows]
    out = {"k": k, "n": len(rows), "arms": {}}
    hdr = f"{'space':14s} {'snap':>7s} {'descent':>8s} {'<2A':>6s} {'q/target':>9s} {'sweeps':>7s} {'FAIL18':>7s}"
    print(f"\nk={k}, {len(rows)} targets (ORACLE DIAGNOSTIC throughout)\n")
    print(hdr); print("-" * len(hdr))
    base = np.array([r["real"]["descent"] for r in rows], float)
    for nm in names + ["desc_rand_start_x4", "unif_rand_start_x4", "desc_long40"]:
        de = np.array([r[nm]["descent"] for r in rows], float)
        sn = np.array([r[nm].get("snap", np.nan) for r in rows], float)
        q = np.array([r[nm]["oracle_queries"] for r in rows], float)
        sw = np.array([r[nm].get("sweeps", np.nan) for r in rows], float)
        out["arms"][nm] = {"summary": I.summary(de), "mean_snap": float(np.nanmean(sn)),
                           "mean_oracle_queries": float(q.mean()),
                           "mean_sweeps": float(np.nanmean(sw)),
                           "FAIL18": float(de[isf].mean()) if isf.any() else None}
        if nm != "real":
            out["arms"][nm]["paired_vs_real"] = I.paired(de, base, folds=folds)
        print(f"{nm:14s} {np.nanmean(sn):7.3f} {de.mean():8.3f} {(de<2.0).mean():6.2f} "
              f"{q.mean():9.0f} {np.nanmean(sw):7.1f} {de[isf].mean() if isf.any() else float('nan'):7.3f}")
    print("\npaired vs the REAL library's descent (positive = the control is WORSE):")
    for nm in names[1:] + ["desc_rand_start_x4", "unif_rand_start_x4", "desc_long40"]:
        d = out["arms"][nm]["paired_vs_real"]
        print(f"  {nm:18s} {d['mean_diff']:+.4f} [{d['ci95'][0]:+.4f},{d['ci95'][1]:+.4f}] "
              f"{d['n_better']}W/{d['n_worse']}L")
    return out


def main(k=4, n_targets=126, sweeps=12):
    tg = I.targets()
    sel = tg if n_targets >= len(tg) else [tg[i] for i in range(0, len(tg), max(1, len(tg) // n_targets))][:n_targets]
    path = os.path.join(RESULTS, f"adv_ceiling_k{k}.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for i, t in enumerate(sel):
        if t["pdb"] in done:
            continue
        rng = np.random.default_rng(9130905 + i)
        rows.append(run_target(t, tg, k, rng, sweeps))
        if i % 10 == 0 or i == len(sel) - 1:
            json.dump({"what": "matched-space control on the representation ceiling",
                       "k": k, "per_target": rows}, open(path, "w"), indent=1)
            r = rows[-1]
            print(f"  {i+1}/{len(sel)} {t['pdb']} real={r['real']['descent']:.2f} "
                  f"margclass={r['marg_class']['descent']:.2f} "
                  f"unif={r['unif_torus']['descent']:.2f} [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "matched-space control on the representation ceiling",
               "k": k, "per_target": rows}, open(path, "w"), indent=1)
    rep = report(rows, k)
    json.dump(rep, open(os.path.join(RESULTS, f"adv_ceiling_k{k}_report.json"), "w"), indent=1)


if __name__ == "__main__":
    kk = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 126
    main(kk, nt)
