"""SPRINT 14, ENER-1 -- the decision-grade Legacy-vs-AMBER matrix.

Arms A-E, H, I are measured on ONE fixed conformational ensemble per target so that the
comparison is genuinely controlled: the ensemble is the cached UNIFORM AMBER stratum
(`amber_kind == 0`, ~1,194 configurations per target, KS 0.026 against the full space).
Nothing in the ensemble was selected by any energy or by the label.

    A  Legacy ranking                     argmin of Legacy over the ensemble
    B  AMBER ranking                      argmin of genuine ff14SB/GBn2 single point
    C  Legacy then AMBER                   Legacy keeps the best 10%, AMBER picks
    D  AMBER then Legacy                   AMBER keeps the best 10%, Legacy picks
    E  a structural objective with neither  the leakage-safe 1-local torsion prior;
                                            plus radius-of-gyration as a second control
    H  both combined                       rank(Legacy) + rank(AMBER), and lambda variants
    I  post-hoc validation only            prior selects; Legacy/AMBER only VETO the worst
                                            decile; measures the clash-rejection role alone

    F  Legacy as the search objective       exact certified optimum over the full 262,144
    G  AMBER as the search objective        coordinate descent on genuine AMBER
    J  refinement only                      `s14/ener_refine.py` (9.1 s per structure)

F and J live in their own modules because they cost real compute.

Metrics per objective: global Spearman rho with RMSD; rho INSIDE THE LOWEST-ENERGY DECILE
(the number that matters -- a search never sees the bulk); ORACLE-snap percentile; the
continuous native's percentile; top-10 recovery of the truly-best 1%; pairwise decision
accuracy; RMSD of the selected structure with a tie-averaged argmin; calibration; and the
full energy-vs-RMSD curve.  Reported per target AND pooled, never pooled alone.

    python -m s14.ener_matrix
"""
from __future__ import annotations

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E

RNG = np.random.default_rng(E.SEED)
KEEP = 0.10                      # cascade keep-fraction


# ------------------------------------------------------------------ the ensemble
def ensemble(pdb):
    """The controlled ensemble: the uniform AMBER stratum. Returns dict of aligned arrays."""
    z = E.enum(pdb)
    m = z.uniform_mask
    idx = z.amber_idx[m]
    d = dict(pdb=pdb, idx=idx, n=len(idx), rmsd=z.rmsd[idx],
             legacy=z.legacy[idx], prior=z.prior[idx],
             amber=z.amber_total[m])
    for t in E.LEG_TERMS:
        d["leg_" + t] = z.leg[t][idx]
    for t in E.AMB_TERMS:
        d["amb_" + t] = z.amb[t][m]
    d["legacy_nosteric"] = z.obj("legacy_nosteric")[idx]
    # radius of gyration -- a pure compactness control with no energy in it
    ca = z.ca(idx)
    d["rg"] = np.sqrt(((ca - ca.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))
    d["_z"] = z
    return d


# ------------------------------------------------------------------ scoring
def score(name, Eobj, R, rng, native_e=None, snap_e=None):
    """All metrics for one objective on one ensemble."""
    out = dict(objective=name,
               rho=E.spearman(Eobj, R),
               rho_decile=E.decile_rho(Eobj, R, 0.10),
               rho_quartile=E.decile_rho(Eobj, R, 0.25),
               sel_rmsd=E.argmin_rmsd(Eobj, R),
               top10_recovery=E.topk_recovery(Eobj, R, 10, 0.01),
               top1pct_mean_rmsd=float(R[np.argsort(Eobj, kind="mergesort")
                                        [:max(len(R) // 100, 1)]].mean()),
               pair_acc=E.pair_accuracy(Eobj, R, 200000, rng),
               pool_mean=float(R.mean()), pool_best=float(R.min()))
    out.update({"calib_" + k: v for k, v in E.calibration(Eobj, R).items()})
    if snap_e is not None:
        out["snap_pctile"] = E.percentile_of(snap_e, Eobj)
    if native_e is not None:
        out["native_pctile"] = E.percentile_of(native_e, Eobj)
    return out


def cascade(first, second, R, keep=KEEP):
    """`first` keeps the best `keep` fraction, `second` picks the argmin inside it."""
    m = max(int(round(keep * len(R))), 2)
    sel = np.argsort(first, kind="mergesort")[:m]
    return E.argmin_rmsd(second[sel], R[sel]), float(R[sel].mean()), int(m)


def veto(selector, guard, R, veto_frac=0.10):
    """Post-hoc VALIDATION only: `guard` vetoes its own worst `veto_frac`; `selector` picks
    from what survives.  This isolates the clash-rejection role from the ranking role."""
    m = max(int(round((1 - veto_frac) * len(R))), 2)
    keep = np.argsort(guard, kind="mergesort")[:m]
    return E.argmin_rmsd(selector[keep], R[keep]), float(R[keep].mean())


# ------------------------------------------------------------------ native references
def native_energies(pdb):
    """ORACLE DIAGNOSTIC: Legacy and AMBER of the CONTINUOUS native torsions, and of the
    in-space ORACLE snap.  Used only to place the native on each objective's axis."""
    from core import energy as et
    from core import geometry as geo
    from core import amber as am
    z = E.enum(pdb)
    sp = z.space()
    phi0, psi0 = sp.ORACLE_native_torsions()
    phi0 = np.asarray(phi0, float)[None]; psi0 = np.asarray(psi0, float)[None]
    c = geo.build_backbone_batch(phi0, psi0)
    comp = et.components_batch(sp.seq, c, phi0, psi0)
    leg_nat = float(sum(float(et.DEFAULT_WEIGHTS.get(t, 0.0)) * np.asarray(comp[t])[0]
                        for t in E.LEG_TERMS))
    snap_leg = float(z.legacy[z.snap_index])
    snap_amb = None
    j = np.flatnonzero(z.amber_idx == z.snap_index)
    if len(j):
        snap_amb = float(z.amber_total[j[0]])
    return dict(leg_native_cont=leg_nat, leg_snap=snap_leg, amb_snap=snap_amb,
                snap_rmsd=float(z.rmsd[z.snap_index]),
                native_cont_comp={t: float(np.asarray(comp[t])[0]) for t in E.LEG_TERMS})


# ------------------------------------------------------------------ main
OBJECTIVES = [
    ("legacy", "legacy", 1),
    ("legacy_nosteric", "legacy_nosteric", 1),
    ("amber", "amber", 1),
    ("prior", "prior", 1),
    ("rg", "rg", 1),
] + [("leg_" + t, "leg_" + t, 1) for t in E.LEG_TERMS] \
  + [("amb_" + t, "amb_" + t, 1) for t in E.AMB_TERMS]


def combos(d):
    """Normalised combinations. Rank normalisation first -- AMBER has no finite variance."""
    rl, ra, rp = E.rank_norm(d["legacy"]), E.rank_norm(d["amber"]), E.rank_norm(d["prior"])
    out = {"H_rank_leg+amb": rl + ra,
           "H_rank_leg+amb+prior": rl + ra + rp,
           "H_rank_prior+leg": rp + rl,
           "H_rank_prior+amb": rp + ra}
    for lam in (0.25, 0.5, 1.0, 2.0, 4.0):
        out[f"H_leg+{lam}amb"] = rl + lam * ra
    # robust-z variants: does a non-monotone-across-objectives scale change anything?
    out["H_robustz_leg+amb"] = E.robust_z(d["legacy"]) + E.robust_z(E.softcore(d["amber"]))
    return out


def main():
    per_target = []
    for pdb in E.ENUM_TARGETS:
        d = ensemble(pdb)
        R = d["rmsd"]
        rng = np.random.default_rng(E.SEED)
        nat = native_energies(pdb)
        rows = []
        for name, key, _ in OBJECTIVES:
            snap_e = None
            if key == "legacy":
                snap_e = nat["leg_snap"]
            elif key == "amber":
                snap_e = nat["amb_snap"]
            rows.append(score(name, d[key], R, rng, snap_e=snap_e))
        for name, v in combos(d).items():
            rows.append(score(name, v, R, rng))
        # cascades
        casc = {}
        c_rmsd, c_pool, m = cascade(d["legacy"], d["amber"], R)
        casc["C_legacy_then_amber"] = dict(sel_rmsd=c_rmsd, sub_pool_mean=c_pool, n_kept=m)
        c_rmsd, c_pool, m = cascade(d["amber"], d["legacy"], R)
        casc["D_amber_then_legacy"] = dict(sel_rmsd=c_rmsd, sub_pool_mean=c_pool, n_kept=m)
        c_rmsd, c_pool, m = cascade(d["prior"], d["legacy"], R)
        casc["prior_then_legacy"] = dict(sel_rmsd=c_rmsd, sub_pool_mean=c_pool, n_kept=m)
        c_rmsd, c_pool, m = cascade(d["prior"], d["amber"], R)
        casc["prior_then_amber"] = dict(sel_rmsd=c_rmsd, sub_pool_mean=c_pool, n_kept=m)
        # NULL CONTROL: a random keep-set of the same size, same second stage
        nsel = max(int(round(KEEP * len(R))), 2)
        rr = np.random.default_rng(E.SEED + 7)
        nulls = {"C_null_random_then_amber":
                 float(np.mean([E.argmin_rmsd(d["amber"][s], R[s])
                                for s in (rr.choice(len(R), nsel, False) for _ in range(64))])),
                 "D_null_random_then_legacy":
                 float(np.mean([E.argmin_rmsd(d["legacy"][s], R[s])
                                for s in (rr.choice(len(R), nsel, False) for _ in range(64))]))}
        # arm I: veto only
        vetoes = {}
        for gname, g in (("legacy", d["legacy"]), ("amber", d["amber"]),
                         ("leg_steric", d["leg_steric"]), ("amb_nonbonded", d["amb_nonbonded"])):
            for vf in (0.10, 0.50):
                s, pm = veto(d["prior"], g, R, vf)
                vetoes[f"I_prior_veto_{gname}_{int(vf*100)}"] = dict(sel_rmsd=s, kept_mean=pm)
        # random selection baseline (tie-averaged over the whole pool == pool mean)
        base = dict(random_draw=float(R.mean()), pool_best=float(R.min()),
                    snap_rmsd=nat["snap_rmsd"])
        # degeneracy diagnostic -- a NaN decile rho means the objective is CONSTANT there
        deg = {}
        for name, key, _ in OBJECTIVES:
            v = np.asarray(d[key], float)
            lo = np.nanmin(v)
            m10 = max(len(v) // 10, 10)
            sel = np.argsort(v, kind="mergesort")[:m10]
            deg[name] = dict(tie_frac=float((v == lo).mean()),
                             n_distinct=int(len(np.unique(v[np.isfinite(v)]))),
                             decile_tie_frac=float((v[sel] == lo).mean()))
        per_target.append(dict(pdb=pdb, n=int(d["n"]), fold=int(d["_z"].fold),
                               rows=rows, cascades=casc, nulls=nulls, vetoes=vetoes,
                               degeneracy=deg, base=base, native=nat))
        print(f"{pdb}: n={d['n']} pool_mean={R.mean():.3f} pool_best={R.min():.3f}", flush=True)
        for r in rows:
            if r["objective"] in ("legacy", "amber", "prior", "leg_steric", "amb_nonbonded"):
                print(f"    {r['objective']:16s} rho {r['rho']:+.3f}  decile {r['rho_decile']:+.3f}"
                      f"  sel {r['sel_rmsd']:.3f}  pair {r['pair_acc']:.3f}"
                      f"  top10rec {r['top10_recovery']:.2f}", flush=True)
    pool = pooled(per_target)
    E.write("ener_matrix", dict(
        what="controlled Legacy-vs-AMBER matrix on the uniform stratum ensemble",
        ensemble="amber_kind==0 (uniform), ~1194 configs/target, no energy or label selection",
        keep_fraction=KEEP, per_target=per_target, pooled=pool),
        n_expected=len(E.ENUM_TARGETS))
    report(pool, per_target)
    return per_target, pool


def pooled(per_target):
    names = [r["objective"] for r in per_target[0]["rows"]]
    out = {}
    for nm in names:
        vals = {k: [] for k in ("rho", "rho_decile", "rho_quartile", "sel_rmsd",
                                "top10_recovery", "pair_acc", "top1pct_mean_rmsd")}
        for t in per_target:
            r = [x for x in t["rows"] if x["objective"] == nm][0]
            for k in vals:
                vals[k].append(r[k])
        rec = {}
        for k, v in vals.items():
            m, ci = E.bootstrap_mean_ci(v)
            rec[k] = m; rec[k + "_ci"] = ci
        # delta vs random draw, paired
        sel = np.array(vals["sel_rmsd"]); base = np.array([t["base"]["random_draw"]
                                                           for t in per_target])
        rec["paired_vs_random"] = I.paired(sel, base, folds=[t["fold"] for t in per_target])
        out[nm] = rec
    for key in ("cascades", "vetoes"):
        for cname in per_target[0][key]:
            sel = np.array([t[key][cname]["sel_rmsd"] for t in per_target])
            base = np.array([t["base"]["random_draw"] for t in per_target])
            out[cname] = dict(sel_rmsd=float(sel.mean()),
                              paired_vs_random=I.paired(sel, base,
                                                        folds=[t["fold"] for t in per_target]))
    for cname in per_target[0]["nulls"]:
        v = np.array([t["nulls"][cname] for t in per_target])
        out[cname] = dict(sel_rmsd=float(v.mean()))
    out["_random_draw"] = float(np.mean([t["base"]["random_draw"] for t in per_target]))
    out["_pool_best"] = float(np.mean([t["base"]["pool_best"] for t in per_target]))
    out["_snap"] = float(np.mean([t["base"]["snap_rmsd"] for t in per_target]))
    return out


def report(pool, per_target):
    print("\n=== POOLED (9 targets, uniform ensemble) ===")
    print(f"random draw {pool['_random_draw']:.3f}  pool best {pool['_pool_best']:.3f}  "
          f"ORACLE snap {pool['_snap']:.3f}")
    print(f"{'objective':22s} {'rho':>7s} {'decile':>8s} {'pair':>6s} {'sel':>7s} "
          f"{'d_rand':>8s} {'CI':>20s} {'W/L':>6s}")
    for nm, r in pool.items():
        if nm.startswith("_") or "rho" not in r:
            continue
        p = r["paired_vs_random"]
        print(f"{nm:22s} {r['rho']:+7.3f} {r['rho_decile']:+8.3f} {r['pair_acc']:6.3f} "
              f"{r['sel_rmsd']:7.3f} {p['mean_diff']:+8.3f} "
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
              f"{p['n_better']}/{p['n_worse']}")
    print("\n--- cascades and vetoes ---")
    for nm, r in pool.items():
        if not isinstance(r, dict):
            continue
        if "paired_vs_random" in r and "rho" not in r:
            p = r["paired_vs_random"]
            print(f"{nm:34s} sel {r['sel_rmsd']:.3f}  d {p['mean_diff']:+.3f} "
                  f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}")
    for nm, r in pool.items():
        if isinstance(r, dict) and set(r) == {"sel_rmsd"}:
            print(f"{nm:34s} sel {r['sel_rmsd']:.3f}   (NULL)")
    print("\n--- degeneracy: fraction of the ensemble tied at the objective's minimum ---")
    for nm in sorted(per_target[0].get("degeneracy", {})):
        v = [t["degeneracy"][nm] for t in per_target]
        print(f"{nm:22s} tie-at-min {np.mean([x['tie_frac'] for x in v]):.4f}  "
              f"distinct values {np.mean([x['n_distinct'] for x in v]):.0f}  "
              f"decile tie-frac {np.mean([x['decile_tie_frac'] for x in v]):.4f}")


if __name__ == "__main__":
    main()
