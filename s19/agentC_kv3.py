"""s19/agentC_kv3.py -- Q1, SECOND DECLARED EXTENSION: WHERE THE GATE MOVES THE MEAN ERROR.

**Written after `agentC_kv2` REFUTED the causal reading of my own pre-registered primary**, and it
says so.  The pre-registration (`PREREG_C.md` section 3.4) stands as written and as reported.

WHAT KV2 KILLED.  Section 3.4 measured that a Legacy gate collapses the ensemble ambiguity `D`,
and that in the EXACT accounting `readout^2 = E_mem^2 - D^2` the diversity channel carries more
than all of the damage.  KV2 then built a gate with the SAME diversity and NO score
(`rand_lowD`, the lowest-`D` of 200 matched-random subsets, D = 1.723 = Legacy's D to three
decimals) and it costs **+0.003 [-0.003, +0.013]** where Legacy costs **+0.063 [+0.010, +0.124]**.

    Diversity collapse is REAL, it is large, and it is NOT the cause.
    An accounting identity is not a mechanism.  This module says so and then asks the right
    question instead.

THE RIGHT QUESTION.  In one FIXED common frame, write each candidate's error `e_i = Y_i - T` and
the pool's mean error `b_pool = mean over all K`.  A gate keeps a subset S, and the operator emits
`b_S = mean over S`.  Then, exactly,

    readout_S^2 = ||b_S||^2 = ||b_pool||^2  +  2 <b_pool, Delta_S>  +  ||Delta_S||^2,
    Delta_S := b_S - b_pool.

A MATCHED-RANDOM gate has `E[Delta_S] = 0`: it pays only the third term, which is the +0.013 A
cost of halving the set, and its cross term averages to zero **whatever its diversity**.  A gate
that hurts more than that must have a NON-ZERO CROSS TERM -- it must displace the ensemble mean
ALONG the direction in which the pool is already wrong.

    ALIGN_S := <b_pool, Delta_S> / ||b_pool||         (angstrom, signed)

is that displacement, in angstrom, measured along the pool's own error direction.  It is zero in
expectation for any score-free gate, and its whole matched-random distribution is measured here
rather than assumed.  **This is the error-covariance question asked in the coordinates where the
answer can be non-trivial**, and unlike `D` it is not forced by an identity to absorb the damage.

ORACLE.  `T`, `e_i`, `b_pool` and `ALIGN` need the native and are labelled ORACLE diagnostics.
Every gate is native-free.

    python -m s19.agentC_kv3
    python -m s19.agentC_kv3 --report
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402
from s19 import agentC_lib as CL                                  # noqa: E402
from s19 import agentC_kv as KV                                   # noqa: E402

R_NULL = 200


def run_target(t, amber_sp):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(W)
    m = max(2, int(round(K * (1.0 - KV.F_PRIMARY))))

    #: ONE fixed common frame for every arm: the full pool's own frame, so `e_i` is the same
    #: object in every arm and subset means are directly comparable.  (The deployed operator
    #: re-medoids per subset; `readout_deployed` from `agentC_kv.json` is the primary, and the
    #: fixed-frame readout is checked against it in the report.)
    Y, Ybar, _C, _dev = CL.common_frame_members(W, PHI, PSI)
    T = I.superpose_batch(nat[None], Ybar)[0]
    e = Y - T                                     # (K, n, 3)   ORACLE
    b_pool = e.mean(0)
    nb = float(np.sqrt(np.sum(b_pool ** 2)))
    u_dir = b_pool / nb if nb > 1e-12 else np.zeros_like(b_pool)

    d_win = np.asarray(I.kabsch_rmsd_batch(W, nat), float)
    comp, legacy = CL.legacy_scores(seq, PHI, PSI)
    scores = {"legacy": legacy,
              "leg_torsion": np.asarray(comp["torsion"], float),
              "leg_contact": np.asarray(comp["contact"], float),
              "leg_steric": np.asarray(comp["steric"], float),
              "helix": CL.helix_score(W, n)}
    if amber_sp is not None:
        scores["amber_sp"] = np.asarray(amber_sp, float)
    P = np.asarray(I.pairwise_rmsd(W), float)
    rng = SD.stable_rng(pdb, "s19C_kv3")

    idxs = {nm: CL.keep_lowest(s, m) for nm, s in scores.items()}
    idxs["legacy_clust"] = CL.gate_cluster_best(scores["legacy"], P, m,
                                                SD.stable_rng(pdb, "s19C_clust"))
    idxs["legacy_spread"] = CL.gate_score_then_spread(scores["legacy"], P, m, pre=2.0)
    idxs["medoid"] = np.sort(np.argsort(P[int(np.argmin(P.mean(1)))])[:m])
    idxs["fps"] = CL._farthest_point(P, m, int(np.argmin(P.mean(1))))

    def stats(idx):
        bS = e[np.asarray(idx, int)].mean(0)
        D = bS - b_pool
        return {"readout_fixed": float(np.sqrt(np.sum(bS ** 2) / n)),
                "align": float(np.sum(D * u_dir) / np.sqrt(n)),
                "dnorm": float(np.sqrt(np.sum(D ** 2) / n)),
                "d_mean": float(d_win[np.asarray(idx, int)].mean())}

    arms = {nm: stats(ix) for nm, ix in idxs.items()}

    #: the MATCHED-RANDOM null, measured rather than assumed: R draws of the same count.
    al, dn, rf = [], [], []
    draws = []
    for _ in range(R_NULL):
        ix = CL.gate_random(K, m, rng)
        s = stats(ix)
        al.append(s["align"]); dn.append(s["dnorm"]); rf.append(s["readout_fixed"])
        draws.append(ix)
    #: the lowest-D score-free subset of the SAME draws, so kv2's arm is present here too.
    Ds = [float(np.sum((Y[ix] - Y[ix].mean(0)) ** 2) / (m * n)) for ix in draws]
    arms["rand_lowD"] = stats(draws[int(np.argmin(Ds))])
    arms["rand_highD"] = stats(draws[int(np.argmax(Ds))])
    arms["rand"] = {"readout_fixed": float(np.mean(rf)), "align": float(np.mean(al)),
                    "dnorm": float(np.mean(dn)),
                    "d_mean": float(np.mean([d_win[ix].mean() for ix in draws]))}
    arms["none"] = stats(np.arange(K))
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]), "K": int(K), "m": int(m),
            "b_pool_rmsd": float(nb / np.sqrt(n)),
            "null_align_mean": float(np.mean(al)), "null_align_sd": float(np.std(al)),
            "null_dnorm_mean": float(np.mean(dn)), "null_readout_mean": float(np.mean(rf)),
            "arms": arms}


def run(targets=None, out="agentC_kv3.json", verbose=True):
    import json
    tg = targets if targets is not None else I.targets()
    asp = KV._load_amber_sp()
    path = os.path.join(CL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"R_NULL": R_NULL, "F_PRIMARY": KV.F_PRIMARY,
           "note": "second declared extension, written after kv2 refuted the causal reading"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.7)
        rows.append(run_target(t, asp.get(t["pdb"])))
        if verbose and len(rows) % 20 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
        KV._write(out, rows, cfg, len(tg))
    KV._write(out, rows, cfg, len(tg))
    return rows


def report():
    o = PL.read_complete(os.path.join(CL.RESULTS, "agentC_kv3.json"), need=126)
    o1 = PL.read_complete(os.path.join(CL.RESULTS, "agentC_kv.json"), need=126)
    o2 = PL.read_complete(os.path.join(CL.RESULTS, "agentC_kv2.json"), need=126)
    rows = o["rows"]
    pdbs = [r["pdb"] for r in rows]
    folds = np.array([r["fold"] for r in rows])
    r1 = {r["pdb"]: r for r in o1["rows"]}
    r2 = {r["pdb"]: r for r in o2["rows"]}
    z = np.zeros(len(rows))

    def c(arm, key):
        return np.array([r["arms"][arm].get(key, np.nan) for r in rows], float)

    def dep(arm):
        """the DEPLOYED readout for the same arm, from the pre-registered artefacts."""
        out = []
        for p in pdbs:
            if arm in r1[p]["arms"]:
                out.append(r1[p]["arms"][arm]["readout"])
            else:
                out.append(r2[p]["arms"][arm]["readout"])
        return np.asarray(out, float)

    arms = [a for a in rows[0]["arms"] if a not in ("rand", "none")]
    print("=" * 104)
    print("Q1 EXTENSION 2 -- WHERE THE GATE MOVES THE ENSEMBLE MEAN.  n = %d, m = %d of K = %d."
          % (len(rows), rows[0]["m"], rows[0]["K"]))
    print("   ALIGN = displacement of the emitted mean along the POOL's own error direction (A).")
    print("   Zero in expectation for ANY score-free gate.  Null measured over %d draws."
          % o["config"]["R_NULL"])
    print("=" * 104)
    print(f"\n   pool mean-error ||b_pool|| = {np.mean([r['b_pool_rmsd'] for r in rows]):.3f} A"
          f"   (= the ungated readout)   matched-random ALIGN null: "
          f"mean {np.mean([r['null_align_mean'] for r in rows]):+.4f} A, "
          f"sd {np.mean([r['null_align_sd'] for r in rows]):.4f} A")
    print(f"\n   {'arm':<15}{'ALIGN':>9}{'|Delta|':>9}{'rd_fixed':>10}{'rd_deployed':>12}"
          f"     ALIGN vs its matched-random null (paired, fold CI)")
    for a in ["legacy", "leg_torsion", "leg_steric", "leg_contact", "amber_sp", "helix",
              "medoid", "rand_lowD", "rand_highD", "fps", "legacy_spread", "legacy_clust"]:
        if a not in rows[0]["arms"]:
            continue
        p = PL.paired(c(a, "align"), c("rand", "align"), folds=folds, names=pdbs)
        ci = p.get("ci_fold", p["ci"])
        print(f"   {a:<15}{np.nanmean(c(a,'align')):>9.4f}{np.nanmean(c(a,'dnorm')):>9.4f}"
              f"{np.nanmean(c(a,'readout_fixed')):>10.3f}{np.nanmean(dep(a)):>12.3f}"
              f"     {p['mean']:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}] {p['W']}W/{p['L']}L")
    print(f"   {'rand':<15}{np.nanmean(c('rand','align')):>9.4f}"
          f"{np.nanmean(c('rand','dnorm')):>9.4f}{np.nanmean(c('rand','readout_fixed')):>10.3f}"
          f"{np.nanmean(dep('rand')):>12.3f}     -- (the null itself)")
    print(f"   {'none':<15}{np.nanmean(c('none','align')):>9.4f}"
          f"{np.nanmean(c('none','dnorm')):>9.4f}{np.nanmean(c('none','readout_fixed')):>10.3f}"
          f"{np.nanmean(dep('none')):>12.3f}     -- (Delta = 0 exactly)")

    print("\n   DOES ALIGN EXPLAIN THE DAMAGE?  Across arms: mean ALIGN vs mean deployed damage.")
    xs, ys, nm = [], [], []
    for a in arms:
        xs.append(float(np.nanmean(c(a, "align"))))
        ys.append(float(np.nanmean(dep(a) - dep("rand"))))
        nm.append(a)
    xs = np.asarray(xs); ys = np.asarray(ys)
    print(f"      corr over {len(xs)} arms = {np.corrcoef(xs, ys)[0,1]:+.4f}")
    for a, x, y in sorted(zip(nm, xs, ys), key=lambda t: -t[1]):
        print(f"      {a:<15} ALIGN {x:+.4f} A    damage vs rand {y:+.4f} A")

    print("\n   THE EXACT TWO-TERM DECOMPOSITION OF THE DAMAGE (A^2, per residue):")
    print("       readout^2 = ||b_pool||^2 + 2*||b_pool||*ALIGN + |Delta|^2      (EXACT)")
    print("       so, against matched-random of the same count, the excess splits into an")
    print("       ALIGNMENT term and a MAGNITUDE term, and into nothing else.")
    b0 = c("none", "readout_fixed")
    ali_r = 2.0 * b0 * c("rand", "align")
    mag_r = c("rand", "dnorm") ** 2
    print(f"      {'arm':<15}{'ALIGN term':>12}{'MAG term':>11}{'total':>10}"
          f"{'align share':>13}     ALIGN term vs its null (fold CI)")
    for a in ["helix", "legacy", "leg_contact", "leg_torsion", "medoid", "amber_sp",
              "leg_steric", "fps", "legacy_spread", "legacy_clust", "rand_lowD",
              "rand_highD"]:
        if a not in rows[0]["arms"]:
            continue
        ali = 2.0 * b0 * c(a, "align") - ali_r
        mag = c(a, "dnorm") ** 2 - mag_r
        tot = ali + mag
        p_ = PL.paired(ali, z, folds=folds, names=pdbs)
        ci = p_.get("ci_fold", p_["ci"])
        sh = np.mean(ali) / np.mean(tot) if abs(np.mean(tot)) > 1e-9 else float("nan")
        print(f"      {a:<15}{np.mean(ali):>12.4f}{np.mean(mag):>11.4f}{np.mean(tot):>10.4f}"
              f"{sh:>13.2f}     {p_['mean']:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]")

    print("\n   POOLED: the five PHYSICS SCORE gates against the four SCORE-FREE gates.")
    sc = [a for a in ("legacy", "leg_torsion", "leg_steric", "leg_contact", "amber_sp")
          if a in rows[0]["arms"]]
    sf = [a for a in ("rand_lowD", "rand_highD", "fps", "medoid") if a in rows[0]["arms"]]
    A_sc = np.mean([c(a, "align") for a in sc], axis=0) - c("rand", "align")
    A_sf = np.mean([c(a, "align") for a in sf], axis=0) - c("rand", "align")
    for lab, v in (("score gates (mean of 5)", A_sc), ("score-FREE gates (mean of 4)", A_sf)):
        q = PL.paired(v, z, folds=folds, names=pdbs)
        ci = q.get("ci_fold", q["ci"])
        print(f"      ALIGN excess, {lab:<30} {q['mean']:+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]"
              f"  med {q['median']:+.4f}  {q['W']}W/{q['L']}L")
    q = PL.paired(A_sf, A_sc, folds=folds, names=pdbs)
    ci = q.get("ci_fold", q["ci"])
    print(f"      score-free MINUS score                       {q['mean']:+.4f} "
          f"[{ci[0]:+.4f}, {ci[1]:+.4f}]  med {q['median']:+.4f}  {q['W']}W/{q['L']}L")

    print("\n   CONCENTRATION CHECK on both pooled statistics (the near-even W/L beside a CI")
    print("   excluding zero is the recorded early warning; drop-top vs a UNIFORM-EFFECT null):")
    for lab, v in (("ALIGN excess, score gates", A_sc),
                   ("readout damage, legacy vs rand", dep("legacy") - dep("rand"))):
        k = CL.concentration(v)
        print(f"      {lab:<32} mean {k['mean']:+.4f}  med {k['median']:+.4f}  "
              f"sd {k['sd']:.4f}  {k['W']}W/{k['L']}L  drop-top5 {k['drop_top_mean']:+.4f}  "
              f"null pct {k['null_pct']:.3f}  -> {k['verdict']}")

    print("\n   IS THERE A NATIVE-FREE HALF?  |Delta| -- how far the gate moved the emitted mean")
    print("   from the UNGATED mean -- needs no native at all.  ALIGN does.  How much of the")
    print("   damage can be seen without the native?")
    dd = np.concatenate([c(a, "dnorm") ** 2 - c("rand", "dnorm") ** 2 for a in arms])
    aa = np.concatenate([2.0 * c("none", "readout_fixed") * (c(a, "align") - c("rand", "align"))
                         for a in arms])
    yy = np.concatenate([dep(a) - dep("rand") for a in arms])
    print(f"      n = {len(yy)} (arm, target) pairs over {len(arms)} arms")
    print(f"      corr(damage, |Delta|^2 excess)  = {np.corrcoef(yy, dd)[0,1]:+.4f}   NATIVE-FREE")
    print(f"      corr(damage, ALIGN term)        = {np.corrcoef(yy, aa)[0,1]:+.4f}   ORACLE")
    print(f"      corr(damage, both summed)       = {np.corrcoef(yy, dd + aa)[0,1]:+.4f}   "
          f"(exact up to the sqrt)")
    print("      per-arm MEANS (the design-time question -- rank the arms without the native):")
    ranks = sorted(((float(np.mean(c(a, 'dnorm'))), float(np.mean(dep(a) - dep('rand'))), a)
                    for a in arms), key=lambda t: t[0])
    from scipy.stats import rankdata as _rk
    rr = np.corrcoef(_rk([t[0] for t in ranks]), _rk([t[1] for t in ranks]))[0, 1]
    print(f"         Spearman(mean |Delta|, mean damage) over {len(ranks)} arms = {rr:+.4f}")

    print("\n   AND PER TARGET, WITHIN THE LEGACY ARM (the causal claim at target level):")
    for a in ("legacy", "helix", "leg_steric"):
        x = c(a, "align") - c("rand", "align")
        y = dep(a) - dep("rand")
        print(f"      {a:<15} corr(ALIGN excess, damage) = {np.corrcoef(x, y)[0,1]:+.4f}"
              f"   n={len(x)}")
    return o


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        run(I.targets()[:5] if a.smoke else None,
            out="_SMOKE_agentC_kv3.json" if a.smoke else "agentC_kv3.json")
