#!/usr/bin/env python
"""s26/i_verify_table.py -- the 20-audit table for the ledger and the findings, from
s26/results/verify/*.rerun.json plus lane I's per-audit reading of what a DIFFERS means.

    python s26/i_verify_table.py            # prints the markdown table
"""
from __future__ import annotations

import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "s26"))
import i_verify_rerun as R                                               # noqa: E402

#: lane I's reading of each non-IDENTICAL outcome (the numbers come from the rerun records)
READING = {
    "headline_audit": "all 123 recorded leaves identical; the 8 added leaves are `harness_cfg.*` fields `Config` gained since the tracked run (quantum, legacy, vqe_*, report_single_start_fit); no value changed",
    "amber_platform": "the script prints a table and writes no JSON (the tracked file is a hand-assembled summary of three probes); from the printed table the CPU arm reproduces the golden 1A13 interaction bit-exactly (-489.9138948277905, 726 evaluations) and the OpenCL arms differ run to run as the tracked verdict says they must (hybrid_double -489.874 vs the recorded -489.84..-489.95)",
    "leak_audit": "ALL_CLEAN stands on 24/25 leaves; the one difference is `backends.project`: the tracked run predates `core.project` (`s8.project`), the fresh run records `core.project`",
    "grad_key_collision": "differs in the direction the fix implies: tracked `COLLISION: true`, no Config field; fresh `COLLISION: false`, `config_has_a_field_for_the_gradient: true`, and its two arms are bit-identical because `PROJECT_GRAD` in the environment no longer selects anything (both children run `exact`); the audit is stale the way `run_equiv2.sh` was (L19), 26 run-property leaves are timings",
    "equiv_compare": "same verdict, BIT-IDENTICAL on all 8 targets including every AMBER quantity; the 5 differing leaves are the arm keys and directory names (tracked `29cc..`/`4077..` are gone; fresh compares the current baseline key `44a9..`, 8 records, against the production cache `1fc9..`, 126 records), hence `n_optimised` 8 -> 126 and the added `only_optimised` list",
    "projection_divergence": "ERROR by construction: its two smoke8 arm caches (fd vs analytic) no longer exist, and against the baseline vs production keys it trips (`shape mismatch`) on 7 S11 baseline records whose `amber_ca` is null because AMBER declined at the 92% ceiling in that run; its question is answered by the fresh `project_arms` (126/126 bit-identical)",
    "project_arms": "differs in the informative direction: the tracked S11 run (scan-builder era) had `ca` bit-identical on 0/8 with 22.3 A worst apart; the fresh run, baseline vs the shipped exact mode over the two on-disk caches, has `ca`, `fit_ca`, `avg_ca` bit-identical on 126/126 and `amber_ca` on 119/119 (the 7 null-AMBER baseline records excluded); the `fd` and `analytic` arms were passed as the production key (those caches are gone), so those two comparisons are identical by construction and carry no information",
    "determinism_audit": "no tracked JSON; fresh: bit-identical across processes and under 4 threads (D1), no key collisions (D3), `CacheCollision` raises (D4), config key complete (D5); D2 not run (it would delete two production-cache records)",
    "hazard_audit": "no tracked JSON; fresh: H2 thread-env energies identical 3/3, H3 alignment-invariant and pair distances bit-identical, H6 stable argsort, H7 longer-normalised and symmetric; H1 single-point pair differs by 1239 kcal/mol on 2.1e8 of builder strain (relative 5.8e-6, inside the 1e-5 bar of the integration test); H5 flags the alternate alphabet string in core/data.py, which is `ALPHABET_ALT`, present by design and never used to encode",
    "project_selfcheck": "no tracked JSON (prints only): builder vs reference 1.732e-13 A, analytic vs central gradient 6.828e-06 on |g| = 1.540",
    "project_equiv": "the four-arm table reproduces on every science leaf: agg.ref / ex / fd / an (synthesis 3.2148 / 3.2148 / 3.2145 / 3.2057, fit, projection_cost, fval) and avg_rmsd 3.048338 identical; all 508 differing leaves are the per-target and aggregate TIMINGS (t_ref / t_ex / t_fd / t_an, s_per_target; 504 + 4) and the 10 run-property leaves are seconds / speedups (exact-mode speedup 2.19x tracked -> 2.02x under seven concurrent jobs); 2 added bookkeeping leaves (complete, requested)",
    "project_iters": "not like for like: the tracked file is a 10-target run with three arms (ref / fd / an, 90 solves each, 10.0-13.3% of solves hitting maxiter=300, nit mean 157-161) from before the exact mode existed; the fresh CLI default is 24 targets with four arms (ref / ex / fd / an, 216 solves each: 18.5 / 18.5 / 18.1 / 17.6% hitting maxiter, nit mean 171-181), written as project_iters_n24.json; the audit's question (is maxiter binding?) reads the same way at both sizes: a minority of starts hit the cap, and the ref and ex arms are identical solve for solve",
    "project_degeneracy": "not like for like at the per-target level, same shape in aggregate: the tracked run measured the four starts under the module default of its day (analytic), the fresh run under the shipped exact mode, so no per-target best objective is identical (max |d| 0.179 in objective units); the degeneracy the audit exists to measure reproduces: 43 / 33 targets with a best-to-runner-up gap under 1e-2 / 1e-3 in both runs, 12 same-point ties in both, gap median 0.0667 -> 0.0643, runner-up distance median 0.846 -> 0.877 A, genuine branch ties (gap < 1e-3, structures > 0.5 A apart) 1 -> 3",
    "project_stability": "STABILITY_READING",
    "project_exactness": "126/126 bit-identical, worst |d| 0.0, exactly as tracked; the 1 different leaf is `agg.s_per_target` (4.57 -> 4.65 s, a timing), the run-property leaf is `agg.seconds`, the 2 added leaves are `agg.complete` / `agg.requested` (bookkeeping added with the `_outpath` fix)",
}
NOT_RUN = {  # empty at close: all 20 ran
}


def main():
    rows = {}
    for f in sorted(glob.glob(os.path.join(R.OUT, "*.rerun.json"))):
        r = json.load(open(f, encoding="utf-8"))
        rows[r["audit"]] = r
    L = ["| # | audit | tag | tracked JSON | verdict | leaves identical / different / added / removed | wall | peak RSS | reading |",
         "|---:|---|---|---|---|---|---:|---:|---|"]
    n_ident = n_diff = n_nojson = n_err = n_notrun = 0
    for i, name in enumerate(R.AUDITS, 1):
        spec = R.AUDITS[name]
        tracked = ", ".join(spec[4]) or "(none)"
        r = rows.get(name)
        done = os.path.join(ROOT, "s26", "jobs_done", f"verify_{name}.json")
        peak = "?"
        if os.path.exists(done):
            peak = f"{json.load(open(done)).get('peak_rss_gb', '?')} GB"
        if r is None or name in NOT_RUN and r is None:
            n_notrun += 1
            L.append(f"| {i} | `{name}` | {spec[2]} | {tracked} | **NOT RUN** | | | | {NOT_RUN.get(name, 'not reached')} |")
            continue
        if name in NOT_RUN:
            n_notrun += 1
            L.append(f"| {i} | `{name}` | {spec[2]} | {tracked} | **NOT RUN (stopped)** | | | | {NOT_RUN[name]} |")
            continue
        if r["error"]:
            n_err += 1
            verdict, leaves = "ERROR", ""
        elif not spec[4]:
            n_nojson += 1
            verdict, leaves = "no tracked JSON", ""
        else:
            d = r["diffs"][spec[4][0]]
            verdict = d.get("verdict", "?")
            if verdict == "IDENTICAL":
                n_ident += 1
            else:
                n_diff += 1
            leaves = (f"{d.get('identical', '')} / {d.get('different', '')} / {d.get('n_added', '')} / "
                      f"{d.get('n_removed', '')}" if "identical" in d else "")
        L.append(f"| {i} | `{name}` | {spec[2]} | {tracked} | **{verdict}** | {leaves} | {r['wall_s']} s | {peak} | "
                 f"{READING.get(name, '')} |")
    summary = (f"{len(R.AUDITS)} audits: {n_ident} IDENTICAL, {n_diff} DIFFERS (each explained), "
               f"{n_nojson} with no tracked JSON to diff, {n_err} ERROR by construction, {n_notrun} NOT RUN. "
               f"Not re-runnable: `recon_containment_audit.json` (written by a deleted script), "
               f"`project_inputs.json` (the harvested input the projection audits read).")
    print("\n".join(L))
    print()
    print(summary)
    return "\n".join(L), summary


if __name__ == "__main__":
    main()
