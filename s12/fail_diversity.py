"""FAIL18 forensics, step 8: is the ESM-contact gain sequence-specific, or is ANY
perturbation of the top-75 cut worth the same?  (The control E7 demands.)

E7 found that blending ESM-2 contact agreement into the filter lifts FAIL18 band recall
from 0/18 to 10/18 -- but a length-matched FOREIGN contact map does almost as well.  If a
random perturbation also does it, the finding is not "ESM knows the contacts", it is
"the shipped top-75 is a degenerate near-duplicate set and any diversification recovers
the band".

Arms (all deployable, all take the top-75 of the K=500 pool):
  shipped     top-75 by the shipped score
  esm         z(-score) + w z(ESM contact agreement)          [own ESM map]
  esmnull     z(-score) + w z(agreement vs a FOREIGN ESM map) [length-matched]
  ncontact    z(-score) + w z(number of 8 A contacts)         [pure compactness axis]
  random      z(-score) + w z(N(0,1))                         [pure noise]
  divmax      greedy max-min diversity under the score (facility-location style)
Also reports the redundancy of the shipped top-75 (mean pairwise CA-RMSD) and the
in-band AUC of every signal against its own null.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I
from s12.fail_esmrescore import load_con, agree, z

W = 0.75
SEED = 0


def auc(score, label):
    from scipy.stats import rankdata
    label = np.asarray(label, bool)
    if label.all() or not label.any():
        return float("nan")
    r = rankdata(score); n1 = label.sum(); n0 = len(label) - n1
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def greedy_div(Wc, sc, m=I.M, lam=1.0):
    """Greedy: start from the score argmin, then repeatedly add the candidate maximising
    z(-score) + lam * z(min RMSD to what is already selected)."""
    s = z(-sc)
    sel = [int(np.argmax(s))]
    dmin = I.kabsch_rmsd_batch(Wc, Wc[sel[0]])
    while len(sel) < m:
        g = s + lam * z(dmin)
        g[sel] = -1e9
        k = int(np.argmax(g)); sel.append(k)
        dmin = np.minimum(dmin, I.kabsch_rmsd_batch(Wc, Wc[k]))
    return np.array(sel, int)


def main():
    rng = np.random.default_rng(SEED)
    con_bank = load_con()
    tg = I.targets(); seqs = [t["seq"] for t in tg]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        sc = I.shipped_score(dg, I.pair_dists(Wc, i, j).astype(np.float32).astype(float))
        band = rr <= rr.min() + I.BAND

        a_esm = agree(Wc, con_bank[seq])
        alt = [s for s in seqs if len(s) == n and s != seq]
        a_null = agree(Wc, con_bank[alt[k % len(alt)]])
        D = np.linalg.norm(Wc[:, :, None] - Wc[:, None], axis=-1)
        m3 = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
        ncon = (D[:, m3] < 8.0).sum(1).astype(float)
        noise = rng.standard_normal(len(Wc))

        arms = {}
        arms["shipped"] = np.argsort(sc, kind="stable")[:I.M]
        for name, v in (("esm", a_esm), ("esmnull", a_null), ("ncontact", ncon),
                        ("random", noise)):
            arms[name] = np.argsort(-(z(-sc) + W * z(v)), kind="stable")[:I.M]
        arms["divmax"] = greedy_div(Wc, sc)

        # redundancy of the shipped top-75
        P = I.pairwise_rmsd(Wc[arms["shipped"]])
        red = float(P[np.triu_indices(I.M, 1)].mean())
        Pp = I.pairwise_rmsd(Wc[rng.choice(len(Wc), I.M, replace=False)])
        red_rand = float(Pp[np.triu_indices(I.M, 1)].mean())

        rows.append(dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                         pool_best=float(rr.min()),
                         redundancy_top75=red, redundancy_random75=red_rand,
                         auc_esm=auc(a_esm, band), auc_esmnull=auc(a_null, band),
                         auc_ncontact=auc(ncon, band), auc_score=auc(-sc, band),
                         arms={a: dict(top75_best=float(rr[s].min()),
                                       top75_mean=float(rr[s].mean()),
                                       recall=int(band[s].any())) for a, s in arms.items()}))
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/126 free={I.free_gb():.1f}", flush=True)
    I.write("fail_diversity", rows)

    f = [r for r in rows if r["fail18"]]; o = [r for r in rows if not r["fail18"]]
    print(f"\n== top-75 redundancy (mean pairwise CA-RMSD inside the selected 75) ==")
    for tag, g in (("FAIL18", f), ("other108", o)):
        print(f"  {tag:9s} shipped {np.mean([r['redundancy_top75'] for r in g]):.3f}  "
              f"random-75 {np.mean([r['redundancy_random75'] for r in g]):.3f}")
    print(f"\n== in-band AUC vs its own null ==")
    for k in ("auc_score", "auc_esm", "auc_esmnull", "auc_ncontact"):
        print(f"  {k:12s} FAIL18 {np.nanmean([r[k] for r in f]):.3f}   "
              f"other108 {np.nanmean([r[k] for r in o]):.3f}")
    print(f"\n== top-75 best (ORACLE eval), w={W} ==")
    print(f"{'arm':10s} {'all126':>8s} {'FAIL18':>8s} {'other108':>9s} {'recall F':>9s} "
          f"{'d126 vs shipped':>16s} {'ci95':>20s} {'drop10':>8s}")
    base = np.array([r["arms"]["shipped"]["top75_best"] for r in rows])
    folds = np.array([r["fold"] for r in rows])
    summ = {}
    for a in rows[0]["arms"]:
        A = np.array([r["arms"][a]["top75_best"] for r in rows])
        AF = np.array([r["arms"][a]["top75_best"] for r in f])
        AO = np.array([r["arms"][a]["top75_best"] for r in o])
        rF = np.mean([r["arms"][a]["recall"] for r in f])
        st = I.paired(A, base, folds=folds)
        summ[a] = dict(all126=float(A.mean()), fail18=float(AF.mean()),
                       other108=float(AO.mean()), recall_fail=float(rF), paired=st)
        print(f"{a:10s} {A.mean():8.3f} {AF.mean():8.3f} {AO.mean():9.3f} {rF:9.2f} "
              f"{st['mean_diff']:16.3f} [{st['ci95'][0]:7.3f},{st['ci95'][1]:7.3f}] "
              f"{(st['drop_top10_mean_diff'] or 0):8.3f}")
    I.write("fail_diversity_summary", summ)


if __name__ == "__main__":
    main()
