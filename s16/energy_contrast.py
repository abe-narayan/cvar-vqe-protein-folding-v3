"""SPRINT 16, ENERGY -- THE THREE-WAY ENERGY CONTRAST, on ONE instrument at IDENTICAL
chain lengths.

    Legacy / Miyazawa-Jernigan   (genuine, eleven components, DEFAULT_WEIGHTS, UNFITTED)
    AMBER ff14SB / GBn2          (genuine single points, cached, BINDING-MASKED)
    the distance-restraint objective (the shipped leave-fold-out distogram's Bayes risk)

on the SAME exhaustively enumerated configuration spaces: 9 targets at n = 9, k = 4,
262,144 configurations each, every one carrying an exact ORACLE CA-RMSD.  Because the
space is COMPLETE, an argmin here is a CERTIFIED optimum -- no search, no budget, no
optimiser to blame.  Chain length is held at 9 for all three arms, which is the point:
the pathology under test is a property of CONTACT POTENTIALS, not of short chains.

WHAT IS MEASURED, per target per arm
  argmin_minus_mean   certified-optimum RMSD minus the mean RMSD of the SAME population.
                      The mean is the honest "a randomly chosen feasible member" control
                      (uniform random is not a control on this instrument, but a random
                      draw from the enumerated population IS exactly what the mean is).
  rho                 Spearman(score, ORACLE RMSD) over that population.
  native percentile   ORACLE DIAGNOSTIC: fraction of the population scoring BELOW the
                      ORACLE-snapped native configuration.  0 = the objective's argmin is
                      the native; 0.5 = chance.
  in_band             mean RMSD of the score-lowest 100 members.

THE BINDING AMBER DATA RULE is applied to every AMBER statistic:
`amber_kind == 0 AND amber_idx != snap_index`, per-target n printed.  The Legacy and
distogram arms are ALSO reported restricted to exactly those positions, so the three-way
comparison is on one identical config set and not on three different ones.
`s14/results/obj_floor.json` is not opened.

CONDITIONING.  Every statistic above is RANK-BASED and therefore invariant under any
monotone conditioning map, which is why the 2.6e20 kcal/mol tail cannot reach it.  The
moment-based diagnostics that the tail DOES reach are quantified in `energy_gate.py`.

PRIOR ART.  Roget et al., arXiv:2606.21241, enumerate every peptide up to length 15 over
12,446 PDB structures and report that the minimum-COST conformation has on average a
LARGER RMSD than a randomly chosen feasible one, with the mechanism that the
Miyazawa-Jernigan contact potential was parameterised on proteins above 50 residues.
The Legacy row below is that result, and must be cited to them, not claimed.  What is
new here is (a) the same measurement for an all-atom force field, whose parameters have
no such length provenance, and (b) the counterpart arm that BREAKS the pattern.

    python -m s16.energy_contrast
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import glob
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s16 import energy_lib as L            # noqa: E402

CACHE = os.path.join(L.CACHE, "disto")
os.makedirs(CACHE, exist_ok=True)


def enum_states(n, k):
    return np.array(np.meshgrid(*[np.arange(k)] * n, indexing="ij")).reshape(n, -1).T


def disto_scores(pdb, seq, fold, n, k, chunk=32768):
    """The shipped distance-restraint objective over the FULL enumeration.

    Computed here rather than read from `s14/cache/hamil_full_*.npz` so that the index
    alignment with the enumeration is established by construction rather than assumed.
    """
    p = os.path.join(CACHE, f"disto_{pdb}.npy")
    if os.path.exists(p):
        return np.load(p)
    from s13 import qarch_lib as Q
    sp = Q.Space(pdb, k, seq=seq, n=n, fold=fold)
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    S = enum_states(n, k)
    out = np.empty(len(S), np.float32)
    for a in range(0, len(S), chunk):
        ca = sp.ca(S[a:a + chunk])
        D = np.linalg.norm(ca[:, i, :] - ca[:, j, :], axis=-1)
        out[a:a + chunk] = I.shipped_score(dg, D).astype(np.float32)
    np.save(p, out)
    return out


def stats(score, rmsd, snap_pos, tag, pop):
    amin, ntie = L.argmin_tied(score, rmsd)
    return {"arm": tag, "population": pop, "n_pop": int(len(score)),
            "argmin_rmsd": amin, "n_tied_argmin": ntie,
            "mean_rmsd": float(np.mean(rmsd)),
            "argmin_minus_mean": amin - float(np.mean(rmsd)),
            "best_rmsd": float(np.min(rmsd)),
            "rho": L.spearman(score, rmsd),
            "ORACLE_native_pctile": (float((score < score[snap_pos]).mean())
                                     if snap_pos is not None else None),
            "in_band100_rmsd": float(rmsd[np.argsort(score)[:100]].mean())}


def run_target(path, with_disto=True):
    d = L.load_enum(path)
    pdb, n, k = d["pdb"], d["n"], d["k"]
    rmsd = d["rmsd"]
    leg = d["legacy"]
    snap = d["snap_index"]
    pos, brep = L.binding_mask(np.load(path))
    idx = d["amber_idx"][pos]                 # positions in the FULL enumeration
    amb = d["amber_total"][pos]

    rows = [stats(leg, rmsd, snap, "LEGACY_unfitted", "full_enumeration")]
    dis = None
    if with_disto:
        t0 = time.time()
        dis = disto_scores(pdb, d["seq"], d["fold"], n, k).astype(np.float64)
        rows.append(stats(dis, rmsd, snap, "DISTOGRAM_restraint", "full_enumeration"))
        print(f"    disto {pdb} in {time.time()-t0:.0f}s", flush=True)

    # -- the MATCHED population: exactly the binding-masked AMBER positions
    rows.append(stats(amb, rmsd[idx], None, "AMBER_ff14SB_GBn2", "binding_masked"))
    rows.append(stats(leg[idx], rmsd[idx], None, "LEGACY_unfitted", "binding_masked"))
    if dis is not None:
        rows.append(stats(dis[idx], rmsd[idx], None, "DISTOGRAM_restraint",
                          "binding_masked"))

    # -- the native's AMBER percentile against the UNBIASED masked population.
    #    ORACLE DIAGNOSTIC.  The native's own single point is at snap_index, which the
    #    binding rule removes from the REFERENCE population -- correctly, because it is
    #    force-included and would otherwise be compared against itself.
    kind = d["amber_kind"]; aidx = d["amber_idx"]
    hit = np.where(aidx == snap)[0]
    nat_amb = float(d["amber_total"][hit[0]]) if len(hit) else None
    nat_pct = float((amb < nat_amb).mean()) if nat_amb is not None else None

    for r in rows:
        r.update(pdb=pdb, n=n, fold=d["fold"])
    return {"pdb": pdb, "n": n, "k": k, "fold": d["fold"], "rows": rows,
            "binding": brep,
            "ORACLE_native_amber_energy": nat_amb,
            "ORACLE_native_amber_pctile_masked": nat_pct,
            "ORACLE_native_legacy_pctile_masked": float((leg[idx] < leg[snap]).mean()),
            "ORACLE_native_rmsd": float(rmsd[snap]),
            "amber_tail_share_raw": L.tail_share(amb),
            "amber_tail_share_signedlog": L.tail_share(L.condition(amb)),
            "amber_tail_share_winsor99": L.tail_share(L.condition(amb, "winsor99")),
            "amber_tail_share_rank": L.tail_share(L.condition(amb, "rank")),
            "amber_max": float(np.max(amb)), "amber_frac_gt_1e6": float((amb > 1e6).mean())}


def run():
    out = []
    for p in L.ENUM_FILES:
        t0 = time.time()
        out.append(run_target(p))
        print(f"[contrast] {out[-1]['pdb']} n={out[-1]['n']} "
              f"binding n={out[-1]['binding']['n_binding']} ({time.time()-t0:.0f}s)",
              flush=True)
    return out


def report(per):
    import collections
    agg = collections.defaultdict(list)
    for t in per:
        for r in t["rows"]:
            agg[(r["arm"], r["population"])].append(r)
    tab = {}
    for (arm, pop), rs in sorted(agg.items()):
        pdbs = [r["pdb"] for r in rs]
        fold = np.array([r["fold"] for r in rs])
        am = np.array([r["argmin_rmsd"] for r in rs])
        mn = np.array([r["mean_rmsd"] for r in rs])
        tab[f"{arm}|{pop}"] = {
            "n_targets": len(rs), "pdbs": pdbs,
            "mean_n_pop": float(np.mean([r["n_pop"] for r in rs])),
            "argmin_rmsd": float(am.mean()),
            "population_mean_rmsd": float(mn.mean()),
            "argmin_minus_mean": I.paired(am, mn, folds=fold, names=pdbs),
            "rho_mean": float(np.mean([r["rho"] for r in rs])),
            "rho_frac_positive": float(np.mean([r["rho"] > 0 for r in rs])),
            "in_band100_rmsd": float(np.mean([r["in_band100_rmsd"] for r in rs])),
            "ORACLE_native_pctile": (float(np.mean([r["ORACLE_native_pctile"] for r in rs]))
                                     if rs[0]["ORACLE_native_pctile"] is not None else None)}
    cond = {k: float(np.mean([t[k] for t in per])) for k in
            ("amber_tail_share_raw", "amber_tail_share_signedlog",
             "amber_tail_share_winsor99", "amber_tail_share_rank", "amber_frac_gt_1e6")}
    cond["amber_max_over_targets"] = float(np.max([t["amber_max"] for t in per]))
    return {"table": tab, "conditioning": cond, "per_target": per,
            "ORACLE_native_amber_pctile_masked":
                float(np.mean([t["ORACLE_native_amber_pctile_masked"] for t in per
                               if t["ORACLE_native_amber_pctile_masked"] is not None])),
            "ORACLE_native_legacy_pctile_masked":
                float(np.mean([t["ORACLE_native_legacy_pctile_masked"] for t in per])),
            "binding_rule": {t["pdb"]: t["binding"] for t in per}}


if __name__ == "__main__":
    per = run()
    rep = report(per)
    print("\nwrote", L.write("energy_contrast", rep))
    print(f"\n{'arm|population':44s} {'n':>3s} {'argmin':>7s} {'popmean':>8s} "
          f"{'d':>8s} {'CI':>20s} {'rho':>7s} {'natpct':>7s}")
    for k, v in rep["table"].items():
        d = v["argmin_minus_mean"]
        np_ = v["ORACLE_native_pctile"]
        print(f"{k:44s} {v['n_targets']:3d} {v['argmin_rmsd']:7.3f} "
              f"{v['population_mean_rmsd']:8.3f} {d['mean_diff']:+8.3f} "
              f"[{d['ci95'][0]:+.3f},{d['ci95'][1]:+.3f}] {v['rho_mean']:+7.3f} "
              f"{('%.3f' % np_) if np_ is not None else '   -  ':>7s}")
    print("\nconditioning:", {k: round(v, 6) for k, v in rep["conditioning"].items()})
    print("ORACLE native percentile, masked: AMBER "
          f"{rep['ORACLE_native_amber_pctile_masked']:.3f}  LEGACY "
          f"{rep['ORACLE_native_legacy_pctile_masked']:.3f}")
