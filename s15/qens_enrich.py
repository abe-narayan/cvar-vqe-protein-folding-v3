"""SPRINT 15 / QENS -- THE CLASSICAL CONTROL THE QRESTRAINT ENRICHMENT RESULT LACKS.

THE CLAIM UNDER TEST (`s15/qrestraint_FINDINGS.md` section 7).  A trained CVaR-VQE state puts
**5.2x to 6.7x** more probability mass below 2.0 A than its own UNTRAINED initialisation
(18.6x under an ORACLE objective), with a clean internal control -- the weakest objective
(`E_ls_pool`) gives **0.98x**, i.e. nothing.  Nine n=9 targets, 3 seeds, alpha = 0.25, a HARD
budget of 8,192 objective evaluations enforced by `s14.vqe_lib.Counter`, exact enumeration of
the 2^18 register so the probabilities are computed and not sampled.

THE GAP, IN THE FINDING AGENT'S OWN WORDS: *"whether a classical sampler at matched budget
reaches the same enrichment is untested here.  Without it, this is a measurement of what the
VQE does, not a claim that only a VQE does it."*

This module runs that control.  Every arm is scored with the IDENTICAL statistic -- the
probability mass its own distribution places below 2.0 / 2.5 / 3.0 A, on the same targets, the
same objectives and the same enumerated register -- and every arm's ENTROPY and EFFECTIVE
SUPPORT are printed beside it, because an enrichment bought by collapsing the distribution is
not an enrichment.

THE ARMS, with their budget class stated:

    uniform      the space's own base rate                                   0 evaluations
    untrained    the ansatz at theta0 -- the VQE's own control               0 evaluations
    untr_matchdiv  draw from theta0 until the trained state's EFFECTIVE SUPPORT
                 distinct configurations have been seen, then weight them uniformly.
                 The coordinator's literal request: the untrained circuit sampled to the
                 same number of distinct configurations.                     0 evaluations
    anneal       classical simulated annealing at the VQE's EXACT budget of 8,192 objective
                 evaluations; its distribution is the empirical visit frequency over the
                 configurations it evaluated.                        **8,192 -- MATCHED**
    anneal_tail  the same annealer, converged tail only (last 8,192/8 evaluations), which is
                 the analogue of reading the VQE's FINAL distribution rather than its whole
                 optimisation history.                               **8,192 -- MATCHED**
    topM         uniform over the M configurations with the LOWEST objective value, with
                 M = the trained state's effective support.  Zero optimisation, matched
                 diversity, native-free.                              262,144 (class B)
    boltz        p ~ exp(-E/T), T solved so the entropy equals the trained state's.
                                                                      262,144 (class B)
    boltz_p0     p ~ p0 * exp(-E/T) at the same matched entropy -- the untrained circuit's own
                 distribution, reweighted by the same objective.      262,144 (class B)

`topM`, `boltz` and `boltz_p0` read the whole tabulated objective, which only an enumerated
instrument allows; they are DIAGNOSTICS of what the objective makes available, not
budget-matched competitors.  `anneal` and `anneal_tail` ARE budget-matched, and they are the
arms that decide the coordinator's question.

ORACLE.  `Enum.rmsd` is read only to score `pmass_below_t` post hoc.  `E_ORACLE_true` is an
ORACLE objective and is labelled as such in every row.

    python -m s15.qens_enrich
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from s15 import qens_lib as Q
from s15 import qrestraint as QR
from s14 import vqe_lib as V
from s14 import vqe_run as R
from s15 import seed as SD

TAG = "enrich"
THRESH = (2.0, 2.5, 3.0)
BUDGET = 8192
ARMS = ("vqe", "uniform", "untrained", "untr_matchdiv", "anneal", "anneal_tail",
        "topM", "boltz", "boltz_p0")


def _stats(p, rmsd, thresh=THRESH):
    p = np.asarray(p, float)
    p = np.maximum(p, 0.0)
    p = p / p.sum()
    m = p > 0
    h = float(-(p[m] * np.log2(p[m])).sum())
    out = {"entropy_bits": h, "eff_support": float(2.0 ** h),
           "n_support": int(m.sum()), "max_prob": float(p.max()),
           "mean_rmsd_under_p": float(p @ rmsd),            # ORACLE, post-hoc
           "mode_rmsd": float(rmsd[int(np.argmax(p))])}     # ORACLE, post-hoc
    for t in thresh:
        out[f"pmass_below_{t}"] = float(p[rmsd < t].sum())  # ORACLE, post-hoc
    return out


def _empirical(idx, N, rmsd):
    cnt = np.bincount(np.asarray(idx, np.int64), minlength=N).astype(float)
    return _stats(cnt, rmsd)


def _entropy_matched_tilt(base, E, h_target):
    """p ~ base * exp(-E/T) with T solved by bisection so H(p) == h_target bits.

    `E` is standardised first so the temperature grid is scale-free.  Returns the tilt whose
    entropy is closest to the target; if the target exceeds H(base) the base is returned.
    """
    E = np.asarray(E, float)
    sd = E.std()
    z = (E - E.min()) / (sd if sd > 1e-12 else 1.0)
    base = np.asarray(base, float)
    base = base / base.sum()

    def at(T):
        w = -z / max(T, 1e-30)
        q = base * np.exp(w - w.max())
        s = q.sum()
        if not np.isfinite(s) or s <= 0:
            q = np.zeros_like(base)
            q[int(np.argmin(z))] = 1.0
            return q
        return q / s

    def H(q):
        m = q > 0
        return float(-(q[m] * np.log2(q[m])).sum())

    lo, hi = 1e-6, 1e6
    if h_target >= H(at(hi)):
        return at(hi)
    for _ in range(120):
        mid = np.sqrt(lo * hi)
        if H(at(mid)) < h_target:
            lo = mid
        else:
            hi = mid
    return at(np.sqrt(lo * hi))


def _topM(E, M, N):
    M = int(max(1, min(N, round(M))))
    idx = np.argpartition(np.asarray(E, float), M - 1)[:M]
    p = np.zeros(N)
    p[idx] = 1.0 / M
    return p


def _untrained_matchdiv(p0, D, rng, N, cap=4_000_000):
    """Draw from `p0` until `D` DISTINCT configurations have been seen; weight them uniformly.

    Drawn in blocks so a very broad `p0` cannot spin; `cap` bounds the total draws and the
    realised distinct count is returned so an unmet match is visible rather than hidden.
    """
    D = int(max(1, min(N, round(D))))
    seen = set()
    drawn = 0
    blk = int(max(64, D // 4))
    while len(seen) < D and drawn < cap:
        b = rng.choice(N, size=min(blk, cap - drawn), p=p0)
        seen.update(b.tolist())
        drawn += len(b)
    u = np.fromiter(seen, np.int64, len(seen))[:D]
    p = np.zeros(N)
    p[u] = 1.0 / len(u)
    return p, int(len(u)), int(drawn)


def run(targets=QR.ENUM9, objs=QR.OBJECTIVES, seeds=(0, 1, 2), budget=BUDGET, alpha=0.25):
    mode = json.load(open(os.path.join(Q.RESULTS, "qrestraint_mode.json")))["targets"]
    prev = Q.ck_load(TAG)
    out = prev.get("cells", {})
    print("=" * 126)
    print("QENS/ENRICH -- the classical control for the QRESTRAINT sub-2 A enrichment result")
    print(f"  {len(targets)} targets x {len(objs)} objectives x {len(seeds)} seeds, "
          f"budget {budget} objective evaluations, alpha {alpha}")
    print("=" * 126)
    for pdb in targets:
        if all(f"{pdb}|{o}|{s}" in out for o in objs for s in seeds):
            print(f"  {pdb} cached", flush=True)
            continue
        V.wait_for_memory(1.6, tag=f"qens_enrich {pdb}")
        t0 = time.time()
        z = QR.Enum(pdb)
        CA = QR.all_ca(z)
        Es, _x = QR.build_objectives(z, CA)
        del CA, _x
        N = z.N
        rmsd = z.rmsd                                        # ORACLE, post-hoc scoring only
        an = R.make_ansatz("mps2f", z.n_qubits)
        from core import quantum as Qq
        bits = Qq.all_bitstrings(z.n_qubits)
        P0 = {}
        for s in seeds:
            th0 = R.init_theta(an, np.random.default_rng(s), 0.8, "random")
            p0 = np.exp(np.asarray(an.logp(th0, bits), float))
            p0 = np.maximum(p0, 0.0)
            P0[s] = p0 / p0.sum()
        del bits
        for name in objs:
            E = np.asarray(Es[name], float)
            h_tr = float(mode[pdb][name]["entropy_bits"])
            eff = float(mode[pdb][name]["eff_support"])
            for s in seeds:
                key = f"{pdb}|{name}|{s}"
                if key in out:
                    continue
                rng = SD.stable_rng("qens_enrich", pdb, name, s)
                p0 = P0[s]
                arms = {}
                arms["uniform"] = _stats(np.full(N, 1.0 / N), rmsd)
                arms["untrained"] = _stats(p0, rmsd)
                pmd, got, drawn = _untrained_matchdiv(p0, eff, rng, N)
                arms["untr_matchdiv"] = {**_stats(pmd, rmsd),
                                         "distinct_target": eff, "distinct_got": got,
                                         "draws_used": drawn}
                c = V.search_anneal(E, z.n, z.k, budget, rng)
                seen = c.all_seen().astype(np.int64)
                arms["anneal"] = {**_empirical(seen, N, rmsd), "evals": int(c.used)}
                arms["anneal_tail"] = {**_empirical(seen[-max(1, len(seen) // 8):], N, rmsd),
                                       "evals": int(c.used)}
                arms["topM"] = _stats(_topM(E, eff, N), rmsd)
                arms["boltz"] = _stats(_entropy_matched_tilt(np.full(N, 1.0 / N), E, h_tr),
                                       rmsd)
                arms["boltz_p0"] = _stats(_entropy_matched_tilt(p0, E, h_tr), rmsd)
                # THE TRAINED ARM, RE-RUN HERE per seed rather than read as a 3-seed mean.
                # Two reasons: the recorded artefact stores only seed-averages, so a paired
                # test against a per-seed classical arm would be mis-paired; and re-running it
                # is an independent reproduction of the enrichment claim itself.
                rv = R.run(E, z.n_qubits, alpha, budget, shots=512, ansatz="mps2f",
                           seed=s, rmsd=rmsd, bits_per_res=z.bits_per_res,
                           exact_dist=True, baseline="const")
                arms["vqe"] = {"entropy_bits": rv["entropy_bits"],
                               "eff_support": rv["eff_support"],
                               "max_prob": rv["max_prob"],
                               "mean_rmsd_under_p": rv["mean_rmsd_under_p"],
                               "mode_rmsd": rv["mode_rmsd"],
                               "evals": int(rv["evals"]),
                               **{f"pmass_below_{t}": float(rv[f"pmass_below_{t}"])
                                  for t in THRESH}}
                arms["_trained_recorded"] = {k: float(v) for k, v in mode[pdb][name].items()
                                             if isinstance(v, (int, float))}
                out[key] = arms
                Q.ck(TAG, "cells", out)
        del Es, z, an, P0
        print(f"  {pdb} done [{time.time()-t0:.0f}s]", flush=True)
    return out


# ================================================================== reporting
def report(thresh=2.0):
    d = Q.ck_load(TAG).get("cells", {})
    if not d:
        print("no cells")
        return {}
    key = f"pmass_below_{thresh}"
    objs = []
    for k in d:
        o = k.split("|")[1]
        if o not in objs:
            objs.append(o)
    print()
    print("=" * 132)
    print(f"ENRICHMENT AT {thresh} A -- every arm's own probability mass below {thresh} A, "
          f"and its ratio to the UNTRAINED circuit")
    print("  ORACLE post-hoc.  `E_ORACLE_true` is an ORACLE DIAGNOSTIC objective in every row.")
    print("=" * 132)
    res = {}
    for name in objs:
        rows = [v for k, v in d.items() if k.split("|")[1] == name]
        base = float(np.mean([r["untrained"][key] for r in rows]))
        rec = float(np.mean([r["_trained_recorded"][key] for r in rows]))
        print(f"\n  --- {name} ---  n = {len(rows)} cells"
              + ("   [ORACLE OBJECTIVE]" if "ORACLE" in name else ""))
        print(f"  {'arm':>16s} {'P(<'+str(thresh)+')':>10s} {'x untrained':>12s} "
              f"{'entropy':>8s} {'eff support':>12s} {'mean RMSD':>10s} {'budget':>9s}")
        cls = {"vqe": "8192", "uniform": "0", "untrained": "0", "untr_matchdiv": "0",
               "anneal": "8192", "anneal_tail": "8192", "topM": "262144",
               "boltz": "262144", "boltz_p0": "262144"}
        cell = {"recorded_qrestraint": {"pmass": rec, "ratio": rec / max(base, 1e-12)}}
        for arm in ARMS:
            if arm not in rows[0]:
                continue
            v = float(np.mean([r[arm][key] for r in rows]))
            cell[arm] = {"pmass": v, "ratio": v / max(base, 1e-12),
                         "entropy_bits": float(np.mean([r[arm]["entropy_bits"]
                                                        for r in rows])),
                         "eff_support": float(np.mean([r[arm]["eff_support"]
                                                       for r in rows])),
                         "mean_rmsd_under_p": float(np.mean([r[arm]["mean_rmsd_under_p"]
                                                             for r in rows]))}
            print(f"  {arm:>16s} {v:10.4f} {v/max(base,1e-12):12.2f} "
                  f"{cell[arm]['entropy_bits']:8.2f} {cell[arm]['eff_support']:12.0f} "
                  f"{cell[arm]['mean_rmsd_under_p']:10.3f} {cls.get(arm,''):>9s}")
        print(f"  {'[recorded VQE]':>16s} {rec:10.4f} {rec/max(base,1e-12):12.2f} "
              f"   (qrestraint's own 3-seed mean, for cross-check)")
        res[name] = cell
    Q.ck(TAG, f"report_{thresh}", res)
    return res


def paired_table(thresh=2.0):
    """The decisive rows: trained VQE minus each classical arm, paired, POSITIVE = VQE better."""
    d = Q.ck_load(TAG).get("cells", {})
    if not d:
        return {}
    key = f"pmass_below_{thresh}"
    objs = []
    for k in d:
        o = k.split("|")[1]
        if o not in objs:
            objs.append(o)
    print()
    print("=" * 132)
    print(f"PAIRED: trained VQE minus each classical arm on P(< {thresh} A).  "
          f"POSITIVE = the VQE puts MORE mass below {thresh} A.")
    print("  cell unit = (target, seed); target unit averages seeds first.  ORACLE post-hoc.")
    print("=" * 132)
    res = {}
    for name in objs:
        rows = {k: v for k, v in d.items() if k.split("|")[1] == name}
        print(f"\n  --- {name} ---" + ("   [ORACLE OBJECTIVE]" if "ORACLE" in name else ""))
        print(f"  {'arm':>16s} {'unit':>7s} {'n':>4s} {'diff':>9s} {'CI95':>21s} "
              f"{'VQEhi/lo':>9s} {'m/sd':>6s} {'verdict':>11s} {'concentration':>20s}")
        for arm in ARMS:
            if arm == "vqe" or arm not in next(iter(rows.values())):
                continue
            tg = sorted({k.split("|")[0] for k in rows})
            a_c = [v["vqe"][key] for v in rows.values()]
            b_c = [v[arm][key] for v in rows.values()]
            a_t = [float(np.mean([v["vqe"][key] for k, v in rows.items()
                                  if k.split("|")[0] == t])) for t in tg]
            b_t = [float(np.mean([v[arm][key] for k, v in rows.items()
                                  if k.split("|")[0] == t])) for t in tg]
            for unit, x, y in (("cell", a_c, b_c), ("target", a_t, b_t)):
                r = Q.verdict_line("", x, y)
                res[f"{name}|{arm}|{unit}"] = r
                print(f"  {arm:>16s} {unit:>7s} {r['n']:4d} {r['mean']:+9.4f} "
                      f"[{r['ci_lo']:+8.4f},{r['ci_hi']:+8.4f}] "
                      f"{r['loss']:4d}/{r['win']:<4d} {r['mean_over_sd']:+6.2f} "
                      f"{r['verdict']:>11s} {r.get('conc_verdict',''):>20s}")
    Q.ck(TAG, f"paired_{thresh}", res)
    return res


def ratio_table(thresh=2.0):
    """The headline `5.2x` is a RATIO OF MEANS over targets, which one target can dominate.

    This prints the per-target ratio distribution instead: the median and the geometric mean
    of `P_arm / P_untrained` computed WITHIN each target, and how many targets show any
    enrichment at all.  A ratio-of-means far above the median ratio means the aggregate is
    carried by the targets with the largest absolute near-native mass.
    """
    d = Q.ck_load(TAG).get("cells", {})
    if not d:
        return {}
    key = f"pmass_below_{thresh}"
    objs = []
    for k in d:
        o = k.split("|")[1]
        if o not in objs:
            objs.append(o)
    print()
    print("=" * 132)
    print(f"PER-TARGET RATIO DISTRIBUTION at {thresh} A (ORACLE post-hoc).  "
          f"`ratio of means` is the headline convention.")
    print("=" * 132)
    res = {}
    for name in objs:
        rows = {k: v for k, v in d.items() if k.split("|")[1] == name}
        tg = sorted({k.split("|")[0] for k in rows})
        print(f"\n  --- {name} ---" + ("   [ORACLE OBJECTIVE]" if "ORACLE" in name else ""))
        print(f"  {'arm':>16s} {'ratio of means':>15s} {'median ratio':>13s} "
              f"{'geo-mean ratio':>15s} {'targets > 1x':>13s}")
        for arm in ARMS:
            if arm not in next(iter(rows.values())):
                continue
            num, den, rat = [], [], []
            for t in tg:
                a = float(np.mean([v[arm][key] for k, v in rows.items()
                                   if k.split("|")[0] == t]))
                b = float(np.mean([v["untrained"][key] for k, v in rows.items()
                                   if k.split("|")[0] == t]))
                num.append(a)
                den.append(b)
                rat.append(a / b if b > 0 else np.nan)
            rat = np.array(rat, float)
            ok = np.isfinite(rat) & (rat > 0)
            rom = float(np.mean(num) / max(np.mean(den), 1e-12))
            med = float(np.nanmedian(rat))
            geo = float(np.exp(np.mean(np.log(rat[ok])))) if ok.any() else np.nan
            frac = float(np.nanmean(rat > 1.0))
            res[f"{name}|{arm}"] = {"ratio_of_means": rom, "median_ratio": med,
                                    "geomean_ratio": geo, "frac_targets_above_1": frac,
                                    "per_target_ratio": rat.tolist(), "targets": tg}
            print(f"  {arm:>16s} {rom:15.2f} {med:13.2f} {geo:15.2f} "
                  f"{frac:13.2f}")
    Q.ck(TAG, f"ratios_{thresh}", res)
    return res


def main():
    import sys
    if "report" not in sys.argv:
        run()
    for t in (2.0, 2.5):
        report(t)
        ratio_table(t)
        paired_table(t)
    print("\nwritten -> s15/results/qens_enrich.json")


if __name__ == "__main__":
    main()
