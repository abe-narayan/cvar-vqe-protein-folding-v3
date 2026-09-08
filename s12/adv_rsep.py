"""ADVERSARIAL AUDIT 1a -- is r_sep = 0.196 an artefact of the partialling?

`s12/obj_errstruct.py:84-90` computes it WITHIN each target by subtracting per-|i-j| SHELL
MEANS from both truth and prediction, then correlating the residuals, then averaging over
the 126 targets.  For a length-n target, shell s holds n-s pairs, so the longest shells
hold 2 pairs and 1 pair.  A 1-pair shell contributes residuals of EXACTLY (0, 0) to both
vectors; a 2-pair shell contributes (+a, -a) and (+b, -b).  Both are degenerate.  The
question is which way that bites.

Six estimators of the same quantity, all on the identical data, plus the matched i.i.d.
null so the CONTRAST (which is what the claim rests on) can be re-checked too:

  E1 as published (per-shell means, all shells)
  E2 drop shells with < 3 pairs
  E3 drop shells with < 5 pairs
  E4 dof-corrected partial correlation (residual dof = npairs - nshells)
  E5 POOLED across targets, with the shell profile removed by a LEAVE-FOLD-OUT regression
     on |i-j| (one-hot over shells + n + sep/n), fitted on the other four folds -- so
     nothing about this target's own pairs is used to define its own baseline
  E6 pooled, within-target z-scored first (removes per-target scale), LFO shell regression
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s12 import instrument as I


def shell_resid(v, sep, minpairs=1):
    r = np.zeros_like(v); keep = np.zeros(len(v), bool)
    for s in np.unique(sep):
        m = sep == s
        if m.sum() >= minpairs:
            r[m] = v[m] - v[m].mean(); keep |= m
    return r, keep


def main(seed=0):
    tg = I.targets()
    rng = np.random.default_rng(seed)
    per = []
    POOL = {"fold": [], "sep": [], "n": [], "t": [], "p": [], "q": []}
    for t in tg:
        u = I.load_univ(t["pdb"]); n = t["n"]
        i, j = I.pair_index(n)
        dt = I.pair_dists(u["nat_ca"][None], i, j)[0]
        dg = I.distogram(t["pdb"], t["seq"], t["fold"])
        vp = np.asarray(dg["expected"], float)
        sep = (j - i).astype(int)
        # matched i.i.d. null at the SAME MAE (the arm the claim is contrasted against)
        e = vp - dt
        z = rng.normal(size=len(dt))
        vq = dt + z * (np.abs(e).mean() / np.abs(z).mean())

        row = dict(pdb=t["pdb"], n=n, fold=t["fold"], npairs=len(dt),
                   nshells=int(len(np.unique(sep))))
        for tag, v in (("dg", vp), ("iid", vq)):
            for mp, nm in ((1, "E1"), (3, "E2"), (5, "E3")):
                rt, k1 = shell_resid(dt, sep, mp); rp, _ = shell_resid(v, sep, mp)
                rt, rp = rt[k1], rp[k1]
                row[f"{tag}_{nm}"] = (float(np.corrcoef(rt, rp)[0, 1])
                                     if rt.std() > 1e-9 and rp.std() > 1e-9 else np.nan)
            # E4: same as E1 but with the dof the shell means consumed
            rt, _ = shell_resid(dt, sep); rp, _ = shell_resid(v, sep)
            r1 = float(np.corrcoef(rt, rp)[0, 1]) if rt.std() > 1e-9 and rp.std() > 1e-9 else np.nan
            dof = row["npairs"] - row["nshells"]
            row[f"{tag}_E4_dof"] = dof
            row[f"{tag}_E4_t"] = (float(r1 * np.sqrt(max(dof - 1, 1) / max(1 - r1 ** 2, 1e-9)))
                                  if np.isfinite(r1) else np.nan)
        POOL["fold"] += [t["fold"]] * len(dt); POOL["sep"] += sep.tolist()
        POOL["n"] += [n] * len(dt); POOL["t"] += dt.tolist()
        POOL["p"] += vp.tolist(); POOL["q"] += vq.tolist()
        per.append(row)
        print(t["pdb"], round(row["dg_E1"], 3), round(row["dg_E2"], 3), round(row["iid_E1"], 3), flush=True)

    P = {k: np.asarray(v, float) for k, v in POOL.items()}
    sep = P["sep"].astype(int); fold = P["fold"].astype(int)
    smax = int(sep.max())
    X = np.column_stack([(sep == s).astype(float) for s in range(2, smax + 1)]
                        + [P["n"], sep / P["n"], np.log(sep)])

    def lfo_resid(y, zscore=False):
        r = np.empty_like(y)
        yy = y.copy()
        if zscore:
            # within-target standardisation, target boundaries implied by npairs
            k = 0
            for row in per:
                m = slice(k, k + row["npairs"])
                seg = yy[m]
                yy[m] = (seg - seg.mean()) / max(seg.std(), 1e-9)
                k += row["npairs"]
        for f in np.unique(fold):
            tr = fold != f; te = fold == f
            A = np.column_stack([X[tr], np.ones(tr.sum())])
            w = np.linalg.lstsq(A, yy[tr], rcond=None)[0]
            r[te] = yy[te] - (np.column_stack([X[te], np.ones(te.sum())]) @ w)
        return r

    out = dict(n=len(per))
    for tag, key in (("dg", "p"), ("iid", "q")):
        for nm in ("E1", "E2", "E3"):
            v = np.array([r[f"{tag}_{nm}"] for r in per], float)
            out[f"{tag}_{nm}"] = dict(mean=float(np.nanmean(v)), median=float(np.nanmedian(v)),
                                      sd=float(np.nanstd(v)), n_nan=int(np.isnan(v).sum()))
        tv = np.array([r[f"{tag}_E4_t"] for r in per], float)
        out[f"{tag}_E4_mean_t"] = float(np.nanmean(tv))
        out[f"{tag}_E4_frac_t_gt2"] = float(np.nanmean(np.abs(tv) > 2))
        rt = lfo_resid(P["t"]); rp = lfo_resid(P[key])
        out[f"{tag}_E5_pooled_lfo"] = float(np.corrcoef(rt, rp)[0, 1])
        rt2 = lfo_resid(P["t"], True); rp2 = lfo_resid(P[key], True)
        out[f"{tag}_E6_pooled_lfo_z"] = float(np.corrcoef(rt2, rp2)[0, 1])
    out["mean_npairs"] = float(np.mean([r["npairs"] for r in per]))
    out["mean_nshells"] = float(np.mean([r["nshells"] for r in per]))
    out["mean_dof_loss_frac"] = float(np.mean([r["nshells"] / r["npairs"] for r in per]))
    print(json.dumps(out, indent=1))
    I.write("adv_rsep", dict(agg=out, rows=per))


if __name__ == "__main__":
    main()
