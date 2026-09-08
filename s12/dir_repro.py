"""D1. Reproduce the o_sign ORACLE arm -- and de-confound it.

fail_pairs.py's `o_sign` arm does TWO things at once:
  (a) it throws away the shipped 17-bin Bayes-risk objective and replaces it with the L1
      risk of a POINT estimate;
  (b) it moves that point estimate 2.0 A in the correct direction.
Only (b) is a sign channel.  The published -1.217 A on FAIL18 is the sum of the two, and
the published `o_magonly` control (correct magnitude, RANDOM sign) gives -1.226 A on the
same targets -- i.e. an arm with ZERO direction information reproduces the whole effect.
That is a strong hint that (a), not (b), is doing the work on FAIL18.

Arms (all on FAIL18 + the 18 length-matched controls, emitted through the real path):
  base      shipped Bayes risk                                    (the incumbent)
  pt        L1 point estimate at E[d], no shift                    <-- THE MISSING CONTROL
  o_signD   L1 point estimate at E[d] + D * sign(true - E[d])      (D = 0.5,1,2,3)
  r_signD   L1 point estimate at E[d] + D * RANDOM sign            (the null at same D)
  o_magonly L1 point estimate at E[d] + |err| * RANDOM sign        (fail_pairs control)
  o_true    L1 point estimate at the true distance                 (perfect distogram)
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
from s12 import obj_common as OC
from s12 import dir_common as DC

DELTAS = (0.5, 1.0, 2.0, 3.0)


def main():
    cset = set(DC.controls())
    names = [t["pdb"] for t in I.targets() if t["pdb"] in set(I.FAIL18) | cset]
    rows = []
    for k, pdb in enumerate(names):
        if I.free_gb() < 1.5:
            print("  waiting on memory", flush=True)
            import time; time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue, grid = d["exp"], d["dtrue"], d["grid"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        rng = np.random.default_rng(k)
        s_rand = rng.choice([-1.0, 1.0], exp.size)

        arms, lam = {}, {}
        e = OC.emit(d, OC.score_bayes(d["D"], d["risk"], grid))
        arms["base"], lam["base"] = e["fit_rmsd"], e["lam_rmsd"]
        sub_best = {"base": e["sub_best"]}

        def run(name, tgt):
            o = DC.emit_target(d, tgt)
            arms[name] = o["fit_rmsd"]; lam[name] = o["lam_rmsd"]
            sub_best[name] = o["sub_best"]

        run("pt", exp)
        for D in DELTAS:
            run(f"o_sign{D}", DC.sign_target(exp, s_true, D))
            run(f"r_sign{D}", DC.sign_target(exp, s_rand, D))
        mag = np.abs(dtrue - exp)
        run("o_magonly", np.maximum(exp + mag * rng.choice([-1.0, 1.0], exp.size), 1.5))
        run("o_true", dtrue)

        rows.append(dict(pdb=pdb, n=d["n"], fold=d["fold"], fail18=pdb in I.FAIL18,
                         pool_best=float(d["rr"].min()), arms=arms, lam=lam,
                         sub_best=sub_best,
                         frac_over=float((exp > dtrue).mean()),
                         mean_abs_err=float(mag.mean())))
        print(f"  {k+1}/{len(names)} {pdb} F={int(rows[-1]['fail18'])} "
              f"base={arms['base']:.2f} pt={arms['pt']:.2f} "
              f"o2={arms['o_sign2.0']:.2f} r2={arms['r_sign2.0']:.2f} "
              f"true={arms['o_true']:.2f} free={I.free_gb():.1f}", flush=True)
        I.write("dir_repro", rows)

    # ------------------------------------------------------------------ report
    keys = list(rows[0]["arms"])
    summ = {}
    print(f"\n{'arm':12s}{'FAIL18':>9s}{'d(base)':>9s}{'d(pt)':>9s}{'W/L vs pt':>11s}"
          f"{'MATCH18':>9s}{'d(base)':>9s}{'d(pt)':>9s}")
    for a in keys:
        row = {}
        for nm, g in (("fail18", [r for r in rows if r["fail18"]]),
                      ("match18", [r for r in rows if not r["fail18"]])):
            x = np.array([r["arms"][a] for r in g])
            b = np.array([r["arms"]["base"] for r in g])
            p = np.array([r["arms"]["pt"] for r in g])
            sb = I.paired(x, b); sp = I.paired(x, p)
            row[nm] = dict(mean=float(x.mean()), d_base=sb["mean_diff"], ci_base=sb["ci95"],
                           wl_base=[sb["n_better"], sb["n_worse"]],
                           d_pt=sp["mean_diff"], ci_pt=sp["ci95"],
                           wl_pt=[sp["n_better"], sp["n_worse"]],
                           sub_best=float(np.mean([r["sub_best"][a] for r in g])),
                           lam=float(np.mean([r["lam"][a] for r in g])))
        summ[a] = row
        f, m = row["fail18"], row["match18"]
        print(f"{a:12s}{f['mean']:9.3f}{f['d_base']:9.3f}{f['d_pt']:9.3f}"
              f"{f['wl_pt'][0]:5d}/{f['wl_pt'][1]:<5d}{m['mean']:9.3f}{m['d_base']:9.3f}{m['d_pt']:9.3f}")
    I.write("dir_repro_summary", summ)


if __name__ == "__main__":
    main()
