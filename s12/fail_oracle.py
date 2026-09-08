"""FAIL18 forensics, step 6: ORACLE counterfactuals -- ranking the missing channels.

Every arm is ORACLE / DIAGNOSTIC.  None is deployable.  Each supplies the system with
EXACTLY ONE piece of native information, then runs the REAL downstream path:
    subset of the K=500 pool  ->  top-75 by the shipped distogram score
                              ->  I.coordinate_average  ->  I.project(lam=0.3)
so the answer is in emitted-structure units (CA-RMSD of the projected chain to model-1).

Arms
  base      shipped: top-75 by score                                   (control)
  o_ss      ORACLE native SS string: keep the 150 pool members with the highest
            agreement to it, then top-75 by the shipped score
  o_rg      ORACLE native radius of gyration: keep the 150 closest in rg
  o_cmap    ORACLE native 8 A contact map: keep the 150 best by contact F1
  o_pairs   ORACLE true distance for the 10 pairs the distogram gets most wrong:
            replace those 10 risk rows with |grid - d_true| and rescore the whole pool
  o_best75  ORACLE which pool members are best: top-75 by true CA-RMSD
  o_best1   ORACLE the single best pool member, projected alone
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)
from s12 import instrument as I

KEEP = 150


def rg_batch(W):
    W = np.asarray(W, float); c = W - W.mean(1, keepdims=True)
    return np.sqrt((c ** 2).sum(2).mean(1))


def emit(W, sub, seq, fold, nat):
    C, _ = I.coordinate_average(W[sub])
    out = I.project(C, seq, fold)
    return float(I.ca_rmsd(out["ca"], nat)), float(I.ca_rmsd(out["fit_ca"], nat))


def run_target(t, pep):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    u = I.load_univ(pdb); p = I.pool_idx(u)
    W = u["W"][p]; PHI = u["PHI"][p]; PSI = u["PSI"][p]
    rr = u["rr"][p]; nat = u["nat_ca"]
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    D = I.pair_dists(W, i, j).astype(np.float32).astype(float)
    sc = I.shipped_score(dg, D)
    order = np.argsort(sc, kind="stable")
    res = {}

    def top75_within(keep_idx):
        keep_idx = np.asarray(keep_idx, int)
        s = keep_idx[np.argsort(sc[keep_idx], kind="stable")[:I.M]]
        return s

    res["base"] = emit(W, order[:I.M], seq, fold, nat)

    # (i) ORACLE native secondary structure
    nat_ss = I.ss_of(pep.phi, pep.psi)
    ssc = np.array([sum(a == b for a, b in zip(I.ss_of(PHI[k], PSI[k]), nat_ss)) / n
                    for k in range(len(W))])
    res["o_ss"] = emit(W, top75_within(np.argsort(-ssc, kind="stable")[:KEEP]), seq, fold, nat)

    # (ii) ORACLE native radius of gyration
    rgn = float(rg_batch(nat[None])[0]); rgs = rg_batch(W)
    res["o_rg"] = emit(W, top75_within(np.argsort(np.abs(rgs - rgn), kind="stable")[:KEEP]),
                       seq, fold, nat)

    # (iii) ORACLE native contact map (8 A, |i-j| >= 3)
    m = np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 3
    Dn = np.linalg.norm(nat[:, None] - nat[None], axis=-1)
    cn = (Dn < 8.0)[m]
    Dc = np.linalg.norm(W[:, :, None] - W[:, None], axis=-1)
    cc = (Dc < 8.0)[:, m]
    tp = (cc & cn).sum(1); f1 = 2 * tp / np.maximum(cc.sum(1) + cn.sum(), 1)
    res["o_cmap"] = emit(W, top75_within(np.argsort(-f1, kind="stable")[:KEEP]), seq, fold, nat)

    # (iv) ORACLE true distance for the 10 most-wrong pairs
    dtrue = I.pair_dists(nat[None], i, j)[0]
    exp = np.asarray(dg["expected"], float)
    worst = np.argsort(-np.abs(exp - dtrue))[:10]
    grid = np.asarray(dg["grid"], float)
    risk = np.array(dg["risk"], float, copy=True)
    for w in worst:
        risk[w] = np.abs(grid - dtrue[w])
    g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
    sc2 = risk[np.arange(risk.shape[0])[None, :], g].mean(1)
    res["o_pairs"] = emit(W, np.argsort(sc2, kind="stable")[:I.M], seq, fold, nat)

    # (v) ORACLE which pool members are best
    ob = np.argsort(rr, kind="stable")
    res["o_best75"] = emit(W, ob[:I.M], seq, fold, nat)
    out1 = I.project(W[ob[0]], seq, fold)
    res["o_best1"] = (float(I.ca_rmsd(out1["ca"], nat)), float(I.ca_rmsd(out1["fit_ca"], nat)))

    return dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                pool_best=float(rr.min()),
                arms={k: dict(ca=v[0], fit_ca=v[1]) for k, v in res.items()})


def main():
    import peptide_db
    peps = {p.pdb: p for p in peptide_db.load()}
    ctrl = json.load(open(os.path.join(ROOT, "s12", "results", "fail_contrast.json")))["controls"]
    cset = {c["ctrl"] for c in ctrl}
    tg = [t for t in I.targets() if t["pdb"] in set(I.FAIL18) | cset]
    rows = []
    for k, t in enumerate(tg):
        rows.append(run_target(t, peps[t["pdb"]]))
        r = rows[-1]
        print(f"  {k+1}/{len(tg)} {t['pdb']} F={int(r['fail18'])} " +
              " ".join(f"{a}={v['ca']:.2f}" for a, v in r["arms"].items()) +
              f"  free={I.free_gb():.1f}", flush=True)
        I.write("fail_oracle", rows)

    f = [r for r in rows if r["fail18"]]; c = [r for r in rows if not r["fail18"]]
    arms = list(rows[0]["arms"].keys())
    print(f"\n== ORACLE counterfactuals, emitted CA-RMSD (lam=0.3 arm) ==")
    print(f"{'arm':10s} {'FAIL18':>8s} {'d vs base':>10s} {'MATCH18':>8s} {'d vs base':>10s} "
          f"{'F ci95':>20s} {'W/L':>7s}")
    bf = np.array([r["arms"]["base"]["ca"] for r in f])
    bc = np.array([r["arms"]["base"]["ca"] for r in c])
    summ = {}
    for a in arms:
        af = np.array([r["arms"][a]["ca"] for r in f])
        ac = np.array([r["arms"][a]["ca"] for r in c])
        st = I.paired(af, bf, names=[r["pdb"] for r in f])
        summ[a] = dict(fail18=float(af.mean()), match18=float(ac.mean()),
                       d_fail=float((af - bf).mean()), d_match=float((ac - bc).mean()),
                       ci=st["ci95"], wl=[st["n_better"], st["n_worse"]],
                       per_target={r["pdb"]: r["arms"][a]["ca"] for r in rows})
        print(f"{a:10s} {af.mean():8.3f} {(af-bf).mean():10.3f} {ac.mean():8.3f} "
              f"{(ac-bc).mean():10.3f} [{st['ci95'][0]:7.3f},{st['ci95'][1]:7.3f}] "
              f"{st['n_better']:3d}/{st['n_worse']:<3d}")
    print("\nper-target FAIL18:")
    print(f"{'pdb':6s}{'poolb':>7s}" + "".join(f"{a:>10s}" for a in arms))
    for r in sorted(f, key=lambda r: r["pdb"]):
        print(f"{r['pdb']:6s}{r['pool_best']:7.2f}" +
              "".join(f"{r['arms'][a]['ca']:10.2f}" for a in arms))
    I.write("fail_oracle_summary", summ)


if __name__ == "__main__":
    main()
