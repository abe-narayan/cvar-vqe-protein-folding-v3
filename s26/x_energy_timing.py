#!/usr/bin/env python
"""Graph 1: wall-clock time per ENERGY EVALUATION versus peptide length (number of residues) for
the three energy levels the CVaR pipeline has consumed, on identical structures; one line per
energy level, with a fitted power law t ~ n^k as the time-complexity summary.

Levels (each is the project's own evaluation path, unchanged):
  DIST  the shipped distogram Bayes-risk score  (`s12.instrument.shipped_score`; the
        production selector's Hamiltonian is zrank of this score)
  LEG   the genuine 11-term Legacy potential at DEFAULT_WEIGHTS
        (`s16.energy_lib.legacy_components_of_windows` + `legacy_total_from`)
  AMB   the genuine ff14SB/GBn2 single point through OpenMM (`s20.qb2_lib.AmberSP._e`),
        the same call the S25 seven-configuration suite cached 63,000 of

Protocol: 24 dev targets, three per peptide length 9 to 16 (the alphabetically first three of
each length, so every length has the same weight on the x-axis), the first 10 candidates of each
target's shipped top-75 (`sub` in the production cache),
the SAME candidate (its CA window for DIST; its phi/psi for LEG and AMB, rebuilt on ideal
geometry exactly as the S25 suite did) evaluated ONE STRUCTURE PER CALL, 5 timed repeats per
candidate after one untimed warm-up call, `time.perf_counter`. Batched throughput (all 20
candidates in one call, where the API is batched) is recorded beside it in the data file and
is NOT what the graph shows. Nothing native is read; no RMSD is computed. Per-target set-up
(distogram load, OpenMM system build) is excluded from the per-evaluation time and recorded
separately.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s16 import energy_lib as EL           # noqa: E402
from core import geometry as geo           # noqa: E402
from core import energy as et              # noqa: E402

N_PER_LEN, N_CAND, N_REP = 3, 10, 5
OUT_JSON = os.path.join(HERE, "results", "x_energy_timing.json")
OUT_PNG = os.path.join(HERE, "figures", "graph1_time_per_energy_evaluation.png")


def pick_targets():
    by_len = {}
    for t in I.targets():
        by_len.setdefault(int(t["n"]), []).append(t)
    out = []
    for n in sorted(by_len):
        out += sorted(by_len[n], key=lambda t: t["pdb"])[:N_PER_LEN]
    return out


def main():
    from s13 import qarch_lib as QA
    from s20 import qb2_lib as QB2
    targets = pick_targets()
    rows, setup = [], []
    if os.path.exists(OUT_JSON + ".rows.tmp"):
        with open(OUT_JSON + ".rows.tmp", encoding="utf-8") as fh:
            ck = json.load(fh)
        done = {x["pdb"] for x in ck["rows"] if x["level"] == "AMB"}
        rows = [x for x in ck["rows"] if x["pdb"] in done]
        setup = [x for x in ck["setup"] if x["pdb"] in done]
        print(f"  resuming: {sorted(done)} already complete on all three levels")
    for t in targets:
        pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
        if any(x["pdb"] == pdb for x in rows):
            continue
        u = I.load_univ(pdb)
        rec = I.shipped_record(pdb)
        sub = np.asarray(rec["sub"], int)[:N_CAND]          # production top-75 prefix
        pool = I.pool_idx(u)                                 # K=500 pool, universe indices
        idx = pool[sub]
        W = np.asarray(u["W"], float)[idx]                               # (N_CAND, n, 3)
        PHI = np.asarray(u["PHI"], float)[idx]
        PSI = np.asarray(u["PSI"], float)[idx]

        # ---- DIST: the shipped Bayes-risk score on the window's own CA coordinates
        t0 = time.perf_counter()
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        setup.append(dict(pdb=pdb, level="DIST", setup_s=time.perf_counter() - t0))
        D_all = I.pair_dists(W, i, j).astype(np.float32).astype(float)
        I.shipped_score(dg, D_all[:1])                                  # warm-up
        for c in range(N_CAND):
            D = D_all[c:c + 1]
            for r in range(N_REP):
                t0 = time.perf_counter()
                I.shipped_score(dg, D)
                rows.append(dict(pdb=pdb, n=n, level="DIST", cand=c, rep=r,
                                 seconds=time.perf_counter() - t0))
        t0 = time.perf_counter(); I.shipped_score(dg, D_all); tb = time.perf_counter() - t0
        setup[-1]["batched_s_per_structure"] = tb / N_CAND

        # ---- LEG: the 11-term Legacy potential on the ideal-geometry rebuild of phi/psi
        setup.append(dict(pdb=pdb, level="LEG", setup_s=0.0))
        EL.legacy_total_from(EL.legacy_components_of_windows(seq, PHI[:1], PSI[:1]))  # warm-up
        for c in range(N_CAND):
            ph, ps = PHI[c:c + 1], PSI[c:c + 1]
            for r in range(N_REP):
                t0 = time.perf_counter()
                EL.legacy_total_from(EL.legacy_components_of_windows(seq, ph, ps))
                rows.append(dict(pdb=pdb, n=n, level="LEG", cand=c, rep=r,
                                 seconds=time.perf_counter() - t0))
        t0 = time.perf_counter()
        EL.legacy_total_from(EL.legacy_components_of_windows(seq, PHI, PSI))
        setup[-1]["batched_s_per_structure"] = (time.perf_counter() - t0) / N_CAND

        # ---- AMB: the ff14SB/GBn2 single point (OpenMM) on the same rebuild.
        # core.amber.memory_guard refuses a new context above 92% RAM (the user's own load
        # sits at 89 to 93% tonight): wait for headroom rather than change the guard, and
        # drop the previous target's context so contexts do not accumulate.
        from core import amber as am
        while len(am._BUILDERS):
            _, h_old = am._BUILDERS.popitem(last=False)
            am._drop_builder(h_old)
        tw = time.time()
        while am.memory_percent() > 90.5 and time.time() - tw < 1800:
            time.sleep(5)
        t0 = time.perf_counter()
        sp = QB2.AmberSP(seq, QA.Space(pdb, 4).rep)
        setup.append(dict(pdb=pdb, level="AMB", setup_s=time.perf_counter() - t0))
        bb = geo.build_backbone_batch(PHI, PSI)
        coords = [{k: v[c] for k, v in bb.items()} for c in range(N_CAND)]
        sp._e(coords[0])                                                # warm-up
        for c in range(N_CAND):
            for r in range(N_REP):
                t0 = time.perf_counter()
                sp._e(coords[c])
                rows.append(dict(pdb=pdb, n=n, level="AMB", cand=c, rep=r,
                                 seconds=time.perf_counter() - t0))
        setup[-1]["batched_s_per_structure"] = None                     # no batched API
        with open(OUT_JSON + ".rows.tmp", "w", encoding="utf-8") as fh:
            json.dump(dict(rows=rows, setup=setup), fh)
        print(f"  {pdb} n={n} fold={fold}: "
              + "  ".join(f"{lv} {np.mean([x['seconds'] for x in rows if x['pdb']==pdb and x['level']==lv])*1e3:8.3f} ms"
                          for lv in ("DIST", "LEG", "AMB")), flush=True)

    # ---- summary per level: mean over candidates of the per-candidate median of 5 repeats,
    #      with the SD across candidates and a bootstrap 95% CI of the mean
    rng = np.random.default_rng(0)
    summary = {}
    for lv in ("DIST", "LEG", "AMB"):
        per_cand = []
        for pdb in sorted({x["pdb"] for x in rows}):
            for c in range(N_CAND):
                s = [x["seconds"] for x in rows if x["pdb"] == pdb and x["level"] == lv and x["cand"] == c]
                per_cand.append(float(np.median(s)))
        pc = np.asarray(per_cand)
        boots = np.array([pc[rng.integers(0, len(pc), len(pc))].mean() for _ in range(4000)])
        summary[lv] = dict(mean_s=float(pc.mean()), sd_s=float(pc.std(ddof=1)),
                           ci95_s=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                           median_s=float(np.median(pc)), n_candidates=int(len(pc)),
                           n_timings=int(len(pc) * N_REP),
                           batched_s_per_structure_mean=(None if lv == "AMB" else float(np.mean(
                               [s["batched_s_per_structure"] for s in setup if s["level"] == lv]))),
                           setup_s_mean=float(np.mean([s["setup_s"] for s in setup if s["level"] == lv])))
    # ---- per length: mean over candidates (per-candidate medians) of the targets of that length
    by_len = {}
    for lv in ("DIST", "LEG", "AMB"):
        by_len[lv] = {}
        for n in sorted({x["n"] for x in rows}):
            pc = []
            for pdb in sorted({x["pdb"] for x in rows if x["n"] == n}):
                for c in range(N_CAND):
                    s = [x["seconds"] for x in rows if x["pdb"] == pdb and x["level"] == lv and x["cand"] == c]
                    pc.append(float(np.median(s)))
            pc = np.asarray(pc)
            boots = np.array([pc[rng.integers(0, len(pc), len(pc))].mean() for _ in range(2000)])
            by_len[lv][int(n)] = dict(mean_s=float(pc.mean()), sd_s=float(pc.std(ddof=1)),
                                      ci95_s=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                                      n_candidates=int(len(pc)))
    # time complexity: least-squares slope of log t against log n (t ~ n^k), with the R^2
    complexity = {}
    for lv in ("DIST", "LEG", "AMB"):
        ns = np.array(sorted(by_len[lv]), float)
        ts = np.array([by_len[lv][int(n)]["mean_s"] for n in ns])
        k, b = np.polyfit(np.log(ns), np.log(ts), 1)
        pred = k * np.log(ns) + b
        r2 = 1 - np.sum((np.log(ts) - pred) ** 2) / np.sum((np.log(ts) - np.log(ts).mean()) ** 2)
        complexity[lv] = dict(exponent_k=float(k), r2=float(r2), fit="t = exp(b) * n**k",
                              b=float(b), n_range=[int(ns.min()), int(ns.max())])
    out = dict(protocol=__doc__, targets=[t["pdb"] for t in targets], n_per_length=N_PER_LEN,
               by_length=by_len, complexity=complexity,
               n_candidates_per_target=N_CAND, n_repeats=N_REP, summary=summary,
               setup=setup, rows=rows,
               environment=dict(python=sys.version.split()[0], numpy=np.__version__,
                                threads="OMP/MKL/OPENBLAS=1", platform="OpenMM CPU, 1 thread"))
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    # ---- Graph 1
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = {"DIST": "Distogram Bayes-risk score (production selector's H)",
              "LEG": "Legacy 11-term potential",
              "AMB": "AMBER ff14SB/GBn2 single point (OpenMM, CPU)"}
    colors = {"DIST": "#4c72b0", "LEG": "#55a868", "AMB": "#c44e52"}
    lv = ["DIST", "LEG", "AMB"]
    fig, ax = plt.subplots(figsize=(8, 5.2), dpi=190)
    fig.patch.set_facecolor("white")
    for k in lv:
        ns = np.array(sorted(by_len[k]), float)
        m = np.array([by_len[k][int(n)]["mean_s"] for n in ns])
        lo = m - np.array([by_len[k][int(n)]["ci95_s"][0] for n in ns])
        hi = np.array([by_len[k][int(n)]["ci95_s"][1] for n in ns]) - m
        kx = complexity[k]["exponent_k"]
        ax.errorbar(ns, m, yerr=[lo, hi], marker="o", capsize=4, color=colors[k],
                    label=f"{labels[k]}: t ~ n^{kx:.2f} (R2 {complexity[k]['r2']:.2f}); "
                          f"mean {summary[k]['mean_s']*1e3:.3g} ms")
    ax.set_yscale("log")
    ax.set_xlabel("number of residues in the peptide")
    ax.set_ylabel("time per energy evaluation (seconds, log scale)")
    ax.set_title("Time per single-structure energy evaluation versus peptide length\n"
                 f"{len(targets)} dev targets (3 per length) x {N_CAND} candidates x {N_REP} repeats\n"
                 "points = mean of per-candidate medians; error bars = bootstrap 95% CI", fontsize=9.5)
    ax.set_xticks(sorted({int(n) for n in by_len["DIST"]}))
    ax.legend(fontsize=8, loc="center right")
    ax.grid(which="both", alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, facecolor="white")
    print("\nTIME COMPLEXITY (t ~ n^k over the measured lengths)")
    for k in lv:
        print(f"  {k:5s} k = {complexity[k]['exponent_k']:+.3f}  R2 {complexity[k]['r2']:.3f}")
    print("PER LENGTH (ms per evaluation):")
    for n in sorted(by_len["DIST"]):
        print("  n=%2d " % n + "  ".join(f"{k} {by_len[k][n]['mean_s']*1e3:9.4f}" for k in lv))
    print("\nSUMMARY (seconds per single-structure evaluation, all lengths pooled)")
    for k in lv:
        s = summary[k]
        print(f"  {k:5s} mean {s['mean_s']:.6f}  sd {s['sd_s']:.6f}  ci95 [{s['ci95_s'][0]:.6f}, {s['ci95_s'][1]:.6f}]"
              f"  batched/structure {s['batched_s_per_structure_mean']}  setup {s['setup_s_mean']:.3f}")
    print("data:", OUT_JSON)
    print("graph:", OUT_PNG)


if __name__ == "__main__":
    main()
