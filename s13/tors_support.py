"""TORSION-PREDICTOR -- what a SEARCH over the predicted density could reach.

THE ARCHITECTURAL POINT.  A torsion-constrained VQE does not consume a point estimate.  It
consumes a SEARCH SPACE: k candidate (phi, psi) states per residue, one qubit register per
residue, and an energy model that picks among them.  So the question that decides whether a
sequence-only torsion channel can support the sprint's architecture is NOT "is the argmax
right" -- it is **"does the top-k support of the predicted per-residue density contain the
native basin, and what is the best structure inside that space?"**

Two measurements per k in {1, 2, 4, 8, 16, 32}:

    recall    fraction of residues whose native 20-degree cell is in the residue's top-k
    snap      ORACLE: per residue independently take the top-k cell nearest the native
              angles, build, measure CA-RMSD.  The greedy ceiling.
    descent   ORACLE: coordinate descent over the same k-per-residue space, minimising
              CA-RMSD directly.  The true ceiling of the restricted space, matching
              `s13/ceiling.py`'s protocol so the two are comparable.

`snap` and `descent` are ORACLE DIAGNOSTICS -- they read native torsions to select inside
the predicted space.  They are the ceiling a perfect downstream objective could reach, not a
result.  The deployable half is `recall` and the size of the space (n_res * log2 k qubits).

Also here, MODEL-FREE and deployable: `c_kmer`, an exact k-mer context lookup in the same
leave-fold-out corpus (no neural network at all).  If the trained model does not beat it,
the model is not the limiting factor; if neither is better than the corpus marginal, the
sequence->local-conformation channel is empty at this corpus size.

    python -m s13.tors_support
"""
from __future__ import annotations
import os, sys, json, time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s13 import tors_common as T                 # noqa: E402
from s13.tors_train import post_path             # noqa: E402
from s13 import tors_eval as EV                  # noqa: E402
from s12 import instrument as I                  # noqa: E402

KS = (1, 2, 4, 8, 16, 32)


# --------------------------------------------------------------- model-free k-mer predictor
def kmer_posterior(fold, seq, orders=(3, 2, 1, 0), min_count=25, alpha=10.0):
    """Back-off empirical p(cell | sequence context) from the fold's legal corpus.

    order 3 = (i-1, i, i+1) triplet, 2 = (i-1, i), 1 = (i), 0 = marginal.  The highest order
    with at least `min_count` observations wins.  No parameters, no training, no tuning
    against any target -- purely the corpus's own conditional frequencies.
    """
    c = T.corpus(fold)
    ok = c["ok"]
    codes = c["codes"]; pid = c["pid"]
    b = T.grid_bin(c["phi"], c["psi"])
    L = len(codes)
    prev = np.where((np.arange(L) > 0) & (np.roll(pid, 1) == pid), np.roll(codes, 1), T.PAD)
    nxt = np.where((np.arange(L) < L - 1) & (np.roll(pid, -1) == pid), np.roll(codes, -1), T.PAD)
    keys = {3: prev * 21 * 21 + codes * 21 + nxt, 2: prev * 21 + codes, 1: codes,
            0: np.zeros(L, np.int64)}
    tabs = {}
    for o in orders:
        k = keys[o][ok]; v = b[ok]
        nk = int(k.max()) + 1
        H = np.zeros((nk, T.NGRID), np.float32)
        np.add.at(H, (k, v), 1.0)
        tabs[o] = H
    cds = T.codes_of(seq); n = len(cds)
    pv = np.concatenate([[T.PAD], cds[:-1]]); nx = np.concatenate([cds[1:], [T.PAD]])
    qk = {3: pv * 21 * 21 + cds * 21 + nx, 2: pv * 21 + cds, 1: cds, 0: np.zeros(n, np.int64)}
    marg = tabs[0][0] / tabs[0][0].sum()
    out = np.zeros((n, T.NGRID), np.float32)
    used = np.zeros(n, np.int64)
    for i in range(n):
        for o in orders:
            H = tabs[o]; kk = int(qk[o][i])
            if kk < len(H) and H[kk].sum() >= min_count:
                # interpolated back-off: an unsmoothed empirical table has hard zeros and
                # would lose the NLL comparison on the smoothing rather than on the signal
                out[i] = (H[kk] + alpha * marg) / (H[kk].sum() + alpha)
                used[i] = o
                break
        else:
            out[i] = marg; used[i] = 0
    return out, used


def build_kmer_arm():
    tg = I.targets()
    out, orders = {}, {}
    for f in sorted({t["fold"] for t in tg}):
        for t in tg:
            if t["fold"] != f:
                continue
            p, u = kmer_posterior(f, t["seq"])
            out[t["pdb"]] = p.astype(np.float32)
            orders[t["pdb"]] = u.tolist()
        print(f"  c_kmer fold {f} done  free={I.free_gb():.2f}", flush=True)
    np.savez_compressed(post_path("c_kmer"), **out)
    json.dump({"backoff_order_used": orders},
              open(os.path.join(T.RESULTS, "tors_kmer_orders.json"), "w"))
    return out


# ------------------------------------------------------------------------ support ceilings
def support_ceiling(arm, nat, ks=KS, sweeps=3):
    z = np.load(post_path(arm))
    tg = I.targets()
    rows = []
    for t in tg:
        P = np.asarray(z[t["pdb"]], np.float64)
        n = t["n"]
        nphi, npsi = nat[t["pdb"]]
        _, _, ca = EV.pool_torsions(t["pdb"])
        nb = T.grid_bin(nphi, npsi)
        order = np.argsort(-P, 1)
        r = {"pdb": t["pdb"], "n": n, "fold": t["fold"]}
        for k in ks:
            top = order[:, :k]                                   # (n, k)
            r[f"recall{k}"] = float(np.mean([(nb[i] in top[i]) for i in range(n)]))
            cph = T.GRID_PHI[top]; cps = T.GRID_PSI[top]         # (n, k)
            # snap: nearest candidate to the native angles, per residue
            d = (np.abs(T.wrap(cph - nphi[:, None])) ** 2
                 + np.abs(T.wrap(cps - npsi[:, None])) ** 2)
            sel = np.argmin(d, 1)
            ph = cph[np.arange(n), sel].copy(); ps = cps[np.arange(n), sel].copy()
            ph[0] = nphi[0]; ps[-1] = npsi[-1]
            r[f"snap{k}"] = I.ca_rmsd(I.build_ca(ph, ps), ca)
            # descent: coordinate descent on the true CA-RMSD inside the same space
            cur = I.ca_rmsd(I.build_ca(ph, ps), ca)
            for _ in range(sweeps):
                moved = False
                for i in range(n):
                    if k == 1:
                        break
                    PH = np.tile(ph, (k, 1)); PS = np.tile(ps, (k, 1))
                    PH[:, i] = cph[i]; PS[:, i] = cps[i]
                    if i == 0:
                        PH[:, 0] = nphi[0]
                    if i == n - 1:
                        PS[:, -1] = npsi[-1]
                    v = I.kabsch_rmsd_batch(I.build_ca(PH, PS), ca)
                    j = int(np.argmin(v))
                    if v[j] < cur - 1e-9:
                        cur = float(v[j]); ph, ps = PH[j].copy(), PS[j].copy(); moved = True
                if not moved:
                    break
            r[f"descent{k}"] = float(cur)
            r[f"qubits{k}"] = int(n * np.log2(k)) if k > 1 else 0
        rows.append(r)
    return rows


def report(arm, rows):
    isf = np.array([r["pdb"] in set(I.FAIL18) for r in rows])
    out = {"arm": arm, "n": len(rows), "k": {}}
    for k in KS:
        sn = np.array([r[f"snap{k}"] for r in rows]); de = np.array([r[f"descent{k}"] for r in rows])
        rc = np.array([r[f"recall{k}"] for r in rows])
        out["k"][str(k)] = {
            "recall_native_cell": float(rc.mean()), "recall_FAIL18": float(rc[isf].mean()),
            "snap_mean": float(sn.mean()), "descent_mean": float(de.mean()),
            "descent_median": float(np.median(de)),
            "descent_frac_under_2.0": float((de < 2.0).mean()),
            "descent_frac_under_1.5": float((de < 1.5).mean()),
            "descent_FAIL18": float(de[isf].mean()), "descent_other108": float(de[~isf].mean()),
            "mean_qubits": float(np.mean([r[f"qubits{k}"] for r in rows]))}
        print(f"  {arm} k={k:2d}: recall={rc.mean():.3f} snap={sn.mean():.3f} "
              f"descent={de.mean():.3f} <2A={(de<2.0).mean():.2f} q={np.mean([r[f'qubits{k}'] for r in rows]):.0f}",
              flush=True)
    return out


def main(argv):
    nat = EV.native_torsions()
    if not os.path.exists(post_path("c_kmer")):
        print("[build] c_kmer (model-free k-mer back-off)", flush=True)
        build_kmer_arm()
    arms = argv or [a for a in ("p_grid", "p_mvm", "e_esm", "c_kmer", "n_marg", "n_shuf")
                    if os.path.exists(post_path(a))]
    out = {}
    for a in arms:
        T.free_ok(1.5, tag=f"support {a}")
        t0 = time.time()
        rows = support_ceiling(a, nat)
        out[a] = report(a, rows)
        print(f"  [{a}] {time.time()-t0:.0f}s", flush=True)
    json.dump({"what": "ORACLE DIAGNOSTIC: ceiling of a search over the predicted torsion "
                       "density; recall is the deployable half",
               "reports": out}, open(os.path.join(T.RESULTS, "tors_support.json"), "w"),
              indent=1, default=str)
    print("wrote s13/results/tors_support.json")


if __name__ == "__main__":
    main(sys.argv[1:])
