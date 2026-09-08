"""FAIL18 forensics, step 7: is the ESM-2 contact channel USABLE (deployable arm)?

E4 found that ESM-2 contact-map agreement is the only native-free signal with positive
in-band skill on FAIL18 (AUC 0.593) and that it is decorrelated from target difficulty.
Here we ask whether adding it to the shipped filter changes WHAT GETS SELECTED, measured
without any projection cost, on all 126 tuning targets.

Arm:  rescore = z(-shipped_score) + w * z(esm_contact_agreement); take top-75.
Metrics (ORACLE labels for evaluation only):
    top75_best   best true CA-RMSD inside the selected 75  (what synthesis can reach)
    recall       does the top-75 keep ANY near-native band member
    top75_mean   mean true CA-RMSD of the selected 75      (what the average sees)

Controls: w = 0 reproduces the shipped filter exactly; a shuffled-contact-map null
(the same agreement statistic against another target's ESM contacts) prices how much of
any gain is the statistic's shape rather than its content.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I

CON = os.path.join(ROOT, "s12", "cache", "esm_con_targets.npz")
WGRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]


def load_con():
    z = np.load(CON, allow_pickle=True)
    return {str(s): np.asarray(c, np.float32) for s, c in zip(z["seqs"], z["con"])}


def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / (x.std() + 1e-9)


def agree(W, con, thr=8.0):
    """Contact-map agreement: correlation of the candidate's binary 8 A contact map with
    the ESM-2 contact probabilities over |i-j| >= 3."""
    n = W.shape[1]
    m = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    D = np.linalg.norm(W[:, :, None] - W[:, None], axis=-1)
    cc = (D < thr)[:, m].astype(np.float32)
    w = con[m].astype(np.float64)
    wc = w - w.mean()
    return (cc @ wc) / max(np.sqrt((wc ** 2).sum()), 1e-9)


def main():
    con_bank = load_con()
    tg = I.targets()
    seqs = [t["seq"] for t in tg]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        W = u["W"][p]; rr = u["rr"][p]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        sc = I.shipped_score(dg, I.pair_dists(W, i, j).astype(np.float32).astype(float))
        con = con_bank[seq]
        a = agree(W, con)
        # null: another target's ESM contacts, same length
        alt = [s for s in seqs if len(s) == n and s != seq]
        a0 = agree(W, con_bank[alt[k % len(alt)]]) if alt else np.zeros_like(a)
        band = rr <= rr.min() + I.BAND
        rec = dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                   pool_best=float(rr.min()), arms={}, null={})
        for w in WGRID:
            for tag, aa in (("arms", a), ("null", a0)):
                s = z(-sc) + w * z(aa)
                sub = np.argsort(-s, kind="stable")[:I.M]
                rec[tag][str(w)] = dict(top75_best=float(rr[sub].min()),
                                        top75_mean=float(rr[sub].mean()),
                                        recall=int(band[sub].any()))
        rows.append(rec)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/126 free={I.free_gb():.1f}", flush=True)
    I.write("fail_esmrescore", rows)

    f = [r for r in rows if r["fail18"]]; o = [r for r in rows if not r["fail18"]]
    print(f"\n== top-75 best (oracle label, evaluation only) vs blend weight w ==")
    print(f"{'w':>5s} {'all126':>8s} {'FAIL18':>8s} {'other108':>9s} {'recall F':>9s} "
          f"{'recall O':>9s} {'null all':>9s} {'null F':>8s}")
    base = np.array([r["arms"]["0.0"]["top75_best"] for r in rows])
    out = {}
    for w in WGRID:
        k = str(w)
        A = np.array([r["arms"][k]["top75_best"] for r in rows])
        N = np.array([r["null"][k]["top75_best"] for r in rows])
        AF = np.array([r["arms"][k]["top75_best"] for r in f])
        AO = np.array([r["arms"][k]["top75_best"] for r in o])
        NF = np.array([r["null"][k]["top75_best"] for r in f])
        rF = np.mean([r["arms"][k]["recall"] for r in f])
        rO = np.mean([r["arms"][k]["recall"] for r in o])
        out[k] = dict(all126=float(A.mean()), fail18=float(AF.mean()),
                      other108=float(AO.mean()), recall_fail=float(rF),
                      recall_other=float(rO), null_all=float(N.mean()),
                      null_fail=float(NF.mean()))
        print(f"{w:5.2f} {A.mean():8.3f} {AF.mean():8.3f} {AO.mean():9.3f} {rF:9.2f} "
              f"{rO:9.2f} {N.mean():9.3f} {NF.mean():8.3f}")
    # paired stats at the best w on all 126
    best = min(WGRID[1:], key=lambda w: np.mean([r["arms"][str(w)]["top75_best"] for r in rows]))
    A = np.array([r["arms"][str(best)]["top75_best"] for r in rows])
    N = np.array([r["null"][str(best)]["top75_best"] for r in rows])
    folds = np.array([r["fold"] for r in rows]); names = [r["pdb"] for r in rows]
    print(f"\n== paired, w={best}, top75_best vs shipped (all 126) ==")
    st = I.paired(A, base, folds=folds, names=names); print(json.dumps(st, indent=1))
    print(f"\n== same w, NULL contacts vs shipped (all 126) ==")
    print(json.dumps(I.paired(N, base, folds=folds), indent=1))
    I.write("fail_esmrescore_summary", dict(grid=out, best_w=best, paired=st))


if __name__ == "__main__":
    main()
