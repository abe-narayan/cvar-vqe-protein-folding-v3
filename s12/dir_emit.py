"""D4/D5. End-to-end measurement of the sign-corrected objective, with all three nulls.

Path: score K=500 -> top-75 -> superpose on medoid -> coordinate average -> project.
Reported against TWO references:
  * `bayes` -- the shipped incumbent (3.204 A synthesis on tuning126)
  * `pt`    -- the same L1 point-estimate objective with NO sign correction.  This is the
              honest reference for a sign channel: it holds the objective FORM fixed.

Arms (DELTA fixed by dir_gate.py on the ORACLE sign, before any model was looked at)
  o_sign     ORACLE sign                                            the ceiling of the channel
  head       trained head, hard sign                                the deployable arm
  head_soft  trained head, shift DELTA*(2p-1)                       confidence-weighted
  sep        separation-only head                     NULL (b)  "you relearned the sep prior"
  shellmaj   per-shell majority sign of training folds NULL (b)  cheapest form of the same
  const      global majority sign of training folds               trivial baseline
  shuf       head trained on shuffled labels          NULL (a)  must be exactly zero
  rand_acc   ORACLE sign corrupted to the head's OWN measured accuracy, 3 seeds
                                                      NULL (c)  the sharpest one
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)

from s12 import instrument as I
from s12 import obj_common as OC
from s12 import dir_common as DC

PRED = os.path.join(ROOT, "s12", "cache", "dir_head_pred.npz")
PROJECT = {"bayes", "pt", "o_sign", "head", "head_soft", "sep", "rand_acc0",
           "o_tmaj", "head_tmaj"}


def r_sep(pred, true, sep):
    """Correlation of prediction vs truth after removing the per-shell means of both."""
    a = pred.copy(); b = true.copy()
    for s in np.unique(sep):
        m = sep == s
        a[m] -= a[m].mean(); b[m] -= b[m].mean()
    if a.std() < 1e-9 or b.std() < 1e-9:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def main(delta=None):
    cfg = json.load(open(os.path.join(ROOT, "s12", "results", "dir_gate_pick.json")))
    D = float(delta if delta is not None else cfg["delta"])
    hd = json.load(open(os.path.join(ROOT, "s12", "results", "dir_head.json")))
    acc = hd["per_target"]
    Z = np.load(PRED)
    tg = I.targets()
    rows, t0 = [], time.time()
    for k, t in enumerate(tg):
        pdb = t["pdb"]
        while I.free_gb() < 1.5:
            print("  waiting on memory", flush=True); time.sleep(20)
        d = OC.load(pdb)
        exp, dtrue, sep = d["exp"], d["dtrue"], d["sep"]
        s_true = np.sign(dtrue - exp); s_true[s_true == 0] = 1.0
        p = {a: np.asarray(Z[f"{pdb}/{a}"], float) for a in
             ("full", "sep", "shellmaj", "const", "shuf")}

        tgts = {"pt": exp,
                "o_sign": DC.sign_target(exp, s_true, D),
                "head": DC.sign_target(exp, np.where(p["full"] > 0.5, 1.0, -1.0), D),
                "head_soft": DC.sign_target(exp, 2.0 * p["full"] - 1.0, D),
                "sep": DC.sign_target(exp, np.where(p["sep"] > 0.5, 1.0, -1.0), D),
                "shellmaj": DC.sign_target(exp, np.where(p["shellmaj"] > 0.5, 1.0, -1.0), D),
                "const": DC.sign_target(exp, np.where(p["const"] > 0.5, 1.0, -1.0), D),
                "shuf": DC.sign_target(exp, np.where(p["shuf"] > 0.5, 1.0, -1.0), D)}
        # one bit per TARGET: which way is this whole peptide wrong?
        tm = 1.0 if s_true.mean() > 0 else -1.0
        tgts["o_tmaj"] = DC.sign_target(exp, np.full(exp.shape, tm), D)
        hm = 1.0 if (p["full"] > 0.5).mean() > 0.5 else -1.0
        tgts["head_tmaj"] = DC.sign_target(exp, np.full(exp.shape, hm), D)
        # the head's pair-specific part only: its sign where it disagrees with its own
        # per-target majority, the majority elsewhere -- identical to `head`; instead we
        # test the complement: head sign with the per-target mean removed is not a sign,
        # so the decomposition is done through o_tmaj / head_tmaj above.
        a_head = float(acc[pdb]["full"])
        for s in range(3):
            rng = np.random.default_rng(7000 + 13 * k + s)
            tgts[f"rand_acc{s}"] = DC.sign_target(
                exp, DC.corrupt_sign(s_true, 1.0 - a_head, rng), D)

        arms, lam, avg, meta = {}, {}, {}, {}
        e = OC.emit(d, OC.score_bayes(d["D"], d["risk"], d["grid"]), lam=0.3)
        arms["bayes"], lam["bayes"], avg["bayes"] = e["fit_rmsd"], e["lam_rmsd"], e["avg_rmsd"]
        meta["bayes"] = dict(sub_best=e["sub_best"], argmin=e["argmin_rmsd"],
                             mae=float(np.abs(exp - dtrue).mean()), r_sep=r_sep(exp, dtrue, sep))
        for nm, v in tgts.items():
            o = DC.emit_target(d, v, lam=0.3 if nm in PROJECT else None)
            arms[nm] = o["fit_rmsd"]; lam[nm] = o["lam_rmsd"]; avg[nm] = o["avg_rmsd"]
            meta[nm] = dict(sub_best=o["sub_best"], argmin=o["argmin_rmsd"],
                            mae=float(np.abs(v - dtrue).mean()), r_sep=r_sep(v, dtrue, sep))
        rows.append(dict(pdb=pdb, n=d["n"], fold=d["fold"], fail18=pdb in I.FAIL18,
                         pool_best=float(d["rr"].min()), acc_head=a_head,
                         arms=arms, lam=lam, avg=avg, meta=meta))
        if (k + 1) % 5 == 0 or k == 0:
            print(f"  {k+1}/126 {pdb} bayes={avg['bayes']:.2f} pt={avg['pt']:.2f} "
                  f"o={avg['o_sign']:.2f} head={avg['head']:.2f} ra={avg['rand_acc0']:.2f} "
                  f"[{time.time()-t0:.0f}s free={I.free_gb():.1f}]", flush=True)
            I.write("dir_emit", rows)
    I.write("dir_emit", rows)
    report(rows, D)


def report(rows, D):
    grp = {"all126": rows, "fail18": [r for r in rows if r["fail18"]],
           "other108": [r for r in rows if not r["fail18"]]}
    keys = list(rows[0]["arms"])
    summ = {"delta": D, "groups": {}}
    for g, rs in grp.items():
        folds = np.array([r["fold"] for r in rs]); nm = [r["pdb"] for r in rs]
        tab = {}
        for a in keys:
            row = {}
            for metric in ("avg", "arms", "lam"):
                x = np.array([r[metric][a] for r in rs], float)
                if np.isnan(x).all():
                    continue
                row[{"avg": "avg", "arms": "fit", "lam": "lam"}[metric]] = float(np.nanmean(x))
            for ref in ("pt", "bayes"):
                for metric, tagm in (("avg", "avg"), ("arms", "fit")):
                    x = np.array([r[metric][a] for r in rs], float)
                    y = np.array([r[metric][ref] for r in rs], float)
                    if np.isnan(x).any() or np.isnan(y).any():
                        continue
                    st = I.paired(x, y, folds=folds, names=nm)
                    row[f"{tagm}_vs_{ref}"] = dict(
                        d=st["mean_diff"], ci=st["ci95"], wl=[st["n_better"], st["n_worse"]],
                        drop10=st["drop_top10_mean_diff"], drop20=st["drop_top20_mean_diff"],
                        per_fold=st.get("per_fold"))
            row["mae"] = float(np.mean([r["meta"][a]["mae"] for r in rs]))
            row["r_sep"] = float(np.nanmean([r["meta"][a]["r_sep"] for r in rs]))
            row["sub_best"] = float(np.mean([r["meta"][a]["sub_best"] for r in rs]))
            row["argmin"] = float(np.mean([r["meta"][a]["argmin"] for r in rs]))
            tab[a] = row
        summ["groups"][g] = tab
    I.write("dir_emit_summary", summ)

    for g in grp:
        t = summ["groups"][g]
        print(f"\n== {g} (n={len(grp[g])}), delta={D} ==")
        print(f"{'arm':11s}{'avg':>8s}{'fit':>8s}{'d_avg vs pt':>13s}{'W/L':>9s}"
              f"{'d_fit vs pt':>13s}{'d_fit vs bayes':>16s}{'MAE':>7s}{'r_sep':>7s}{'t75best':>9s}")
        for a in keys:
            r = t[a]
            dv = r.get("avg_vs_pt", {}); df = r.get("fit_vs_pt", {}); db = r.get("fit_vs_bayes", {})
            print(f"{a:11s}{r['avg']:8.3f}{r.get('fit', float('nan')):8.3f}"
                  f"{dv.get('d', float('nan')):13.3f}"
                  f"{dv.get('wl', ['', ''])[0]:>4}/{dv.get('wl', ['', ''])[1]:<4}"
                  f"{df.get('d', float('nan')):13.3f}{db.get('d', float('nan')):16.3f}"
                  f"{r['mae']:7.3f}{r['r_sep']:7.3f}{r['sub_best']:9.3f}")


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else None)
