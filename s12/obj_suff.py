"""Step 4/5: SUFFICIENCY, dimensionality and null controls for the objective channel.

All arms here are ORACLE/DIAGNOSTIC unless marked deployable.  Everything is measured on the
top-75 coordinate average (validated proxy, +0.16-0.20 A below the projected structure).

  dim        how low-dimensional is the usable information: full matrix vs the per-|i-j|
             profile vs its first shells vs a single scalar (Rg, max distance)
  localgeom  is the CA-CA distance matrix the right REPRESENTATION: CA virtual bond angles
             and virtual dihedrals scored against the native values
  shell      shell ablation: swap one shell of the shipped profile for the truth, and the
             reverse (leave-one-shell-out of the true profile)
  pairfrac   replace a fraction of pairs by their true distance, chosen by |error| /
             at random / by longest separation / by shortest separation
  arm2       simulate the coordinator's arm-2 inference form: a training-fold constant
             profile (LFO) plus a PERFECT within-shell residual
  null       partial pool_best and pool_mean out of the profile-MAE <-> emitted relation
"""
from __future__ import annotations
import os, json, collections
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from s12 import instrument as I
from s12 import obj_common as OC


def zs(x):
    s = np.std(x)
    return (x - np.mean(x)) / (s if s > 1e-12 else 1.0)


def virt(X):
    b = X[..., 1:, :] - X[..., :-1, :]
    n = np.linalg.norm(b, axis=-1)
    c = (b[..., :-1, :] * b[..., 1:, :]).sum(-1) / (n[..., :-1] * n[..., 1:])
    ang = np.arccos(np.clip(-c, -1, 1))
    b0, b1, b2 = b[..., :-2, :], b[..., 1:-1, :], b[..., 2:, :]
    n1 = np.cross(b0, b1); n2 = np.cross(b1, b2)
    m = np.cross(n1, b1 / np.linalg.norm(b1, axis=-1, keepdims=True))
    return ang, np.arctan2((m * n2).sum(-1), (n1 * n2).sum(-1))


def main():
    tg = OC.targets()
    prof = {}
    for t in tg:
        d = OC.load(t["pdb"]); sep = d["sep"]; us = list(np.unique(sep))
        prof[t["pdb"]] = dict(fold=t["fold"], us=us,
                              qt=np.array([d["dtrue"][sep == s].mean() for s in us]),
                              qp=np.array([d["exp"][sep == s].mean() for s in us]))
    const = {}
    for f in sorted({v["fold"] for v in prof.values()}):
        acc = collections.defaultdict(list)
        for p, v in prof.items():
            if v["fold"] == f:
                continue                                   # LFO: labels from other folds only
            for k, s in enumerate(v["us"]):
                acc[s].append(v["qt"][k])
        const[f] = {int(s): float(np.mean(a)) for s, a in acc.items()}

    SS = list(range(2, 10)); FR = [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0]
    R = collections.defaultdict(list)
    addin, leave = collections.defaultdict(list), collections.defaultdict(list)
    pairfrac = collections.defaultdict(lambda: collections.defaultdict(list))
    meta = []
    rng = np.random.default_rng(0)
    for t in tg:
        d = OC.load(t["pdb"]); sep = d["sep"]; us = list(np.unique(sep)); v = prof[t["pdb"]]
        emit = lambda sc: OC.emit(d, sc, lam=None)["avg_rmsd"]
        L1 = lambda tv: emit(OC.score_l1(d["D"], tv))
        expand = lambda q: np.array([q[us.index(s)] for s in sep])
        ep, tp = expand(v["qp"]), expand(v["qt"])
        cp = np.array([const[t["fold"]].get(int(s), v["qp"][us.index(s)]) for s in sep])
        r_d, r_t = d["exp"] - ep, d["dtrue"] - tp
        meta.append(dict(pdb=t["pdb"], fold=t["fold"], pool_best=float(d["rr"].min()),
                         pool_mean=float(d["rr"].mean()),
                         prof_mae=float(np.abs(v["qp"] - v["qt"]).mean()),
                         pair_mae=float(np.abs(d["exp"] - d["dtrue"]).mean())))
        # --- dimensionality
        R["full_oracle"].append(L1(d["dtrue"]))
        R["disto"].append(L1(d["exp"]))
        R["random"].append(emit(rng.standard_normal(len(d["D"]))))
        P = np.stack([d["D"][:, sep == s].mean(1) for s in us], 1)
        for k, nm in ((1, "prof_s2"), (2, "prof_s2_3"), (len(us), "prof_all")):
            R[nm].append(emit(np.abs(P[:, :k] - v["qt"][:k][None, :]).mean(1)))
        Wc = d["W"] - d["W"].mean(1, keepdims=True); rg = np.sqrt((Wc ** 2).sum(-1).mean(-1))
        Nc = d["nat"] - d["nat"].mean(0); rgt = float(np.sqrt((Nc ** 2).sum(-1).mean()))
        R["rg_only"].append(emit(np.abs(rg - rgt)))
        R["dmax_only"].append(emit(np.abs(d["D"].max(1) - d["dtrue"].max())))
        # --- local geometry
        aC, dC = virt(d["W"]); aN, dN = virt(d["nat"][None]); aN, dN = aN[0], dN[0]
        s_ang = np.abs(aC - aN[None]).mean(1)
        s_dih = np.abs(np.arctan2(np.sin(dC - dN[None]), np.cos(dC - dN[None]))).mean(1)
        R["loc_ang"].append(emit(s_ang)); R["loc_dih"].append(emit(s_dih))
        R["loc_both"].append(emit(zs(s_ang) + zs(s_dih)))
        R["dist_plus_loc"].append(emit(zs(OC.score_l1(d["D"], d["dtrue"])) + zs(zs(s_ang) + zs(s_dih))))
        # --- arm-2 simulation
        R["trueprof_distoresid"].append(L1(tp + r_d))
        R["distoprof_trueresid"].append(L1(ep + r_t))
        R["constprof_trueresid"].append(L1(cp + r_t))
        R["constprof_distoresid"].append(L1(cp + r_d))
        R["constprof_noresid"].append(L1(cp))
        # --- shell ablation
        for s in SS:
            if s not in us:
                continue
            q = v["qp"].copy(); q[us.index(s)] = v["qt"][us.index(s)]
            addin[s].append(L1(expand(q) + r_d))
            q2 = v["qt"].copy(); q2[us.index(s)] = v["qp"][us.index(s)]
            leave[s].append(L1(expand(q2) + r_d))
        # --- pair fraction
        e = np.abs(d["exp"] - d["dtrue"]); npair = len(e)
        for mode, o in (("worst", np.argsort(-e)), ("random", rng.permutation(npair)),
                        ("longsep", np.argsort(-sep)), ("shortsep", np.argsort(sep))):
            for fr in FR:
                k = int(round(fr * npair)); tv = d["exp"].copy(); tv[o[:k]] = d["dtrue"][o[:k]]
                pairfrac[mode][fr].append(L1(tv))

    pdbs = [m["pdb"] for m in meta]
    f18 = np.array([p in I.FAIL18 for p in pdbs])
    summ = {k: dict(mean=float(np.mean(v)), fail18=float(np.mean(np.array(v)[f18])),
                    other=float(np.mean(np.array(v)[~f18])),
                    frac2=float(np.mean(np.array(v) < 2.0))) for k, v in R.items()}
    # --- null control
    pb = np.array([m["pool_best"] for m in meta]); pm = np.array([m["pool_mean"] for m in meta])
    pmae = np.array([m["prof_mae"] for m in meta]); dmae = np.array([m["pair_mae"] for m in meta])
    emit_d = np.array(R["disto"]); gain = emit_d - np.array(R["full_oracle"])

    def part(y):
        X = np.column_stack([np.ones(len(y)), pb, pm])
        return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]

    null = dict(
        r_profmae_emitted=float(np.corrcoef(pmae, emit_d)[0, 1]),
        r_pairmae_emitted=float(np.corrcoef(dmae, emit_d)[0, 1]),
        r_profmae_emitted_partialled=float(np.corrcoef(part(pmae), part(emit_d))[0, 1]),
        r_pairmae_emitted_partialled=float(np.corrcoef(part(dmae), part(emit_d))[0, 1]),
        r_profmae_gain=float(np.corrcoef(pmae, gain)[0, 1]),
        r_profmae_gain_partialled=float(np.corrcoef(part(pmae), part(gain))[0, 1]))

    out = dict(summary=summ, meta=meta,
               shell_addin={int(s): float(np.mean(v)) for s, v in addin.items()},
               shell_leave={int(s): float(np.mean(v)) for s, v in leave.items()},
               pairfrac={m: {str(f): float(np.mean(v)) for f, v in dd.items()}
                         for m, dd in pairfrac.items()},
               null=null, per_target=({k: list(map(float, v)) for k, v in R.items()}))
    I.write("obj_suff", out)
    for k in sorted(summ, key=lambda k: summ[k]["mean"]):
        s = summ[k]
        print(f"{k:24s} {s['mean']:6.3f}  F18 {s['fail18']:6.3f}  oth {s['other']:6.3f}  <2A {s['frac2']:.2f}")
    print("\nshell ablation  addin / leave-out:")
    for s in sorted(addin):
        print(f"  s={s}  {np.mean(addin[s]):.3f}  {np.mean(leave[s]):.3f}")
    print("\npair fraction:")
    print("  " + "".join(f"{f:>8.2f}" for f in FR))
    for m in pairfrac:
        print(f"  {m:10s}" + "".join(f"{np.mean(pairfrac[m][f]):8.3f}" for f in FR))
    print("\nnull control:", json.dumps(null, indent=1))


if __name__ == "__main__":
    main()
