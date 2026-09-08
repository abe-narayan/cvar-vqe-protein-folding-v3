"""FAIL18 forensics, step 14: anatomy of the 10-pair channel.

E12: giving the objective the TRUE distance for the 10 pairs it gets most wrong is worth
1.49 A of emitted RMSD on FAIL18 -- the largest non-selection channel.  Two questions:

  Q1 (deployable)  does the distogram KNOW which pairs it got wrong?  Correlate its own
                   per-pair sd / entropy with |expected - true| (ORACLE evaluation).
  Q2 (ORACLE)      is it the DIRECTION or the MAGNITUDE that carries the value?
                   o_sign    : every pair's expected distance is nudged by a FIXED
                               delta in the correct direction (sign of the error only)
                   o_pairs10 : the 10 worst pairs get their exact true distance (E12)
                   o_pairsK  : the same for K = 3, 5, 20, and for all pairs
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

DELTA = 2.0


def main():
    ctrl = json.load(open(os.path.join(ROOT, "s12", "results", "fail_contrast.json")))["controls"]
    cset = {c["ctrl"] for c in ctrl}
    tg = [t for t in I.targets() if t["pdb"] in set(I.FAIL18) | cset]
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb); p = I.pool_idx(u)
        Wc = u["W"][p]; rr = u["rr"][p]; nat = u["nat_ca"]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n); sep = j - i
        D = I.pair_dists(Wc, i, j).astype(np.float32).astype(float)
        grid = np.asarray(dg["grid"], float); risk0 = np.asarray(dg["risk"], float)
        exp = np.asarray(dg["expected"], float); sd = np.asarray(dg["sd"], float)
        prob = np.asarray(dg["prob"], float)
        ent = -(prob * np.log(prob + 1e-9)).sum(1)
        dtrue = I.pair_dists(nat[None], i, j)[0]
        err = np.abs(exp - dtrue)

        # Q1: does the model know?
        q1 = dict(r_sd_err=float(np.corrcoef(sd, err)[0, 1]),
                  r_ent_err=float(np.corrcoef(ent, err)[0, 1]),
                  r_sep_err=float(np.corrcoef(sep, err)[0, 1]),
                  r_exp_err=float(np.corrcoef(exp, err)[0, 1]))
        worst = np.argsort(-err)[:10]
        q1["worst10_mean_sep"] = float(sep[worst].mean())
        q1["all_mean_sep"] = float(sep.mean())
        q1["worst10_sd_pct"] = float(np.mean([(sd < sd[w]).mean() for w in worst]))
        q1["worst10_frac_overpredicted"] = float((exp[worst] > dtrue[worst]).mean())
        q1["frac_overpredicted"] = float((exp > dtrue).mean())

        def emit_with(risk):
            g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
            sc = risk[np.arange(risk.shape[0])[None, :], g].mean(1)
            sub = np.argsort(sc, kind="stable")[:I.M]
            C, _ = I.coordinate_average(Wc[sub])
            o = I.project(C, seq, fold)
            return dict(ca=float(I.ca_rmsd(o["ca"], nat)),
                        top75_best=float(rr[sub].min()), top75_mean=float(rr[sub].mean()))

        arms = {"base": emit_with(risk0)}
        for K in (3, 5, 10, 20, len(err)):
            R = risk0.copy()
            for w in np.argsort(-err)[:K]:
                R[w] = np.abs(grid - dtrue[w])
            arms[f"o_pairs{K if K < len(err) else 'ALL'}"] = emit_with(R)
        # sign-only: shift every expected distance DELTA in the right direction
        tgt = exp + DELTA * np.sign(dtrue - exp)
        arms["o_sign"] = emit_with(np.abs(grid[None, :] - tgt[:, None]))
        # magnitude-only control: shift by the right amount but a RANDOM sign
        rng = np.random.default_rng(k)
        mag = np.abs(dtrue - exp)
        tgt2 = exp + mag * rng.choice([-1.0, 1.0], len(exp))
        arms["o_magonly"] = emit_with(np.abs(grid[None, :] - tgt2[:, None]))

        rows.append(dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18,
                         pool_best=float(rr.min()), q1=q1, arms=arms))
        print(f"  {k+1}/{len(tg)} {pdb} F={int(rows[-1]['fail18'])} " +
              " ".join(f"{a}={v['ca']:.2f}" for a, v in arms.items()) +
              f" free={I.free_gb():.1f}", flush=True)
        I.write("fail_pairs", rows)

    f = [r for r in rows if r["fail18"]]; c = [r for r in rows if not r["fail18"]]
    print("\n== Q1: does the distogram know which pairs it got wrong? ==")
    for k in sorted(rows[0]["q1"]):
        print(f"  {k:28s} FAIL18 {np.mean([r['q1'][k] for r in f]):7.3f}   "
              f"MATCH18 {np.mean([r['q1'][k] for r in c]):7.3f}")
    print("\n== Q2: emitted CA-RMSD by oracle-distance arm ==")
    bf = np.array([r["arms"]["base"]["ca"] for r in f])
    bc = np.array([r["arms"]["base"]["ca"] for r in c])
    summ = {}
    print(f"{'arm':12s}{'FAIL18':>9s}{'d':>8s}{'W/L':>8s}{'MATCH18':>9s}{'d':>8s}{'t75b F':>8s}")
    for a in rows[0]["arms"]:
        af = np.array([r["arms"][a]["ca"] for r in f]); ac = np.array([r["arms"][a]["ca"] for r in c])
        st = I.paired(af, bf)
        summ[a] = dict(fail18=float(af.mean()), d_fail=st["mean_diff"], ci=st["ci95"],
                       wl=[st["n_better"], st["n_worse"]], match18=float(ac.mean()),
                       d_match=float((ac - bc).mean()),
                       t75b_fail=float(np.mean([r["arms"][a]["top75_best"] for r in f])))
        print(f"{a:12s}{af.mean():9.3f}{st['mean_diff']:8.3f}{st['n_better']:4d}/{st['n_worse']:<3d}"
              f"{ac.mean():9.3f}{(ac-bc).mean():8.3f}{summ[a]['t75b_fail']:8.3f}")
    I.write("fail_pairs_summary", dict(arms=summ,
                                       q1_fail={k: float(np.mean([r["q1"][k] for r in f])) for k in rows[0]["q1"]},
                                       q1_match={k: float(np.mean([r["q1"][k] for r in c])) for k in rows[0]["q1"]}))


if __name__ == "__main__":
    main()
