"""FAIL18 forensics, step 9: does the diversity lever survive the EMITTED-STRUCTURE test?

E6 showed diversified selection lifts FAIL18 band recall 0/18 -> 14/18 and top-75 best
4.677 -> 3.306 A, at the cost of the selected set's MEAN quality.  The coordinate average
tracks the set mean, so the only test that counts is the emitted CA-RMSD through the real
synthesis path: I.coordinate_average -> I.project(lam=0.3).

Arms (all DEPLOYABLE -- no native information anywhere):
  shipped75   top-75 by the shipped distogram score                       (control)
  divmax75    greedy max-min diversity under the score, 75 members
  random75    z(-score) + 0.75 z(N(0,1)), 75 members
  divmax25    greedy max-min diversity under the score, 25 members
  score25     top-25 by the shipped score                       (cardinality control)
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
from s12.fail_esmrescore import z
from s12.fail_diversity import greedy_div

SEED = 0


def main():
    rng = np.random.default_rng(SEED)
    ctrl = json.load(open(os.path.join(ROOT, "s12", "results", "fail_contrast.json")))["controls"]
    cset = {c["ctrl"] for c in ctrl}
    tg = [t for t in I.targets() if t["pdb"] in set(I.FAIL18) | cset]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]; nat = u["nat_ca"]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        sc = I.shipped_score(dg, I.pair_dists(Wc, i, j).astype(np.float32).astype(float))
        order = np.argsort(sc, kind="stable")
        sels = dict(
            shipped75=order[:I.M],
            score25=order[:25],
            divmax75=greedy_div(Wc, sc, m=I.M),
            divmax25=greedy_div(Wc, sc, m=25),
            random75=np.argsort(-(z(-sc) + 0.75 * z(rng.standard_normal(len(Wc)))),
                                kind="stable")[:I.M],
        )
        rec = dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                   pool_best=float(rr.min()), arms={})
        for a, s in sels.items():
            C, _ = I.coordinate_average(Wc[s])
            out = I.project(C, seq, fold)
            rec["arms"][a] = dict(ca=float(I.ca_rmsd(out["ca"], nat)),
                                  fit_ca=float(I.ca_rmsd(out["fit_ca"], nat)),
                                  set_best=float(rr[s].min()), set_mean=float(rr[s].mean()))
        rows.append(rec)
        print(f"  {k+1}/{len(tg)} {pdb} F={int(rec['fail18'])} " +
              " ".join(f"{a}={v['ca']:.2f}" for a, v in rec["arms"].items()) +
              f"  free={I.free_gb():.1f}", flush=True)
        I.write("fail_emit", rows)

    f = [r for r in rows if r["fail18"]]; c = [r for r in rows if not r["fail18"]]
    base_f = np.array([r["arms"]["shipped75"]["ca"] for r in f])
    base_c = np.array([r["arms"]["shipped75"]["ca"] for r in c])
    print(f"\n== emitted CA-RMSD (lam=0.3) ==")
    print(f"{'arm':11s} {'FAIL18':>8s} {'d':>7s} {'W/L':>7s} {'MATCH18':>8s} {'d':>7s} {'W/L':>7s} {'sum36':>8s}")
    summ = {}
    for a in rows[0]["arms"]:
        af = np.array([r["arms"][a]["ca"] for r in f])
        ac = np.array([r["arms"][a]["ca"] for r in c])
        sf = I.paired(af, base_f, names=[r["pdb"] for r in f])
        scc = I.paired(ac, base_c, names=[r["pdb"] for r in c])
        summ[a] = dict(fail18=float(af.mean()), match18=float(ac.mean()),
                       d_fail=sf["mean_diff"], ci_fail=sf["ci95"],
                       wl_fail=[sf["n_better"], sf["n_worse"]],
                       d_match=scc["mean_diff"], ci_match=scc["ci95"],
                       wl_match=[scc["n_better"], scc["n_worse"]],
                       mean36=float(np.concatenate([af, ac]).mean()))
        print(f"{a:11s} {af.mean():8.3f} {sf['mean_diff']:7.3f} "
              f"{sf['n_better']:3d}/{sf['n_worse']:<3d} {ac.mean():8.3f} {scc['mean_diff']:7.3f} "
              f"{scc['n_better']:3d}/{scc['n_worse']:<3d} {np.concatenate([af,ac]).mean():8.3f}")
    print("\nper-target FAIL18:")
    arms = list(rows[0]["arms"])
    print(f"{'pdb':6s}{'poolb':>7s}" + "".join(f"{a:>11s}" for a in arms))
    for r in sorted(f, key=lambda r: r["pdb"]):
        print(f"{r['pdb']:6s}{r['pool_best']:7.2f}" +
              "".join(f"{r['arms'][a]['ca']:11.2f}" for a in arms))
    I.write("fail_emit_summary", summ)


if __name__ == "__main__":
    main()
