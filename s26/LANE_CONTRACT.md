# S26 LANE CONTRACT

Every S26 agent reads this file first and follows it exactly. It is the shared part of every
lane brief. The coordinator (the session that spawned you) is the only writer of
`s26/lanes.json`; everything else in `s26/` is append-or-create by the lane that owns it.

Repository: `C:\Users\abena\Protein-Folding-Algorithm`, branch `s26` (off `ae86a124`).
Python: `C:\Users\abena\miniforge_3\python.exe` (3.13; numpy 2.5, torch 2.13 cpu, pennylane
0.45, openmm 8.5, fair-esm 2.0, psutil 7.2, python-pptx 1.0.2 as `pptx`).

## 0. What the project is, in one paragraph

A structure predictor for 9 to 16 residue peptides, and the measuring instrument around it.
Retrieval (BLOSUM62, K=500 real protein windows) -> distogram Bayes-risk filter (top 75) ->
coordinate average -> multi-start ideal-geometry projection (ramah, lambda 0.3) -> optional
restrained AMBER ff14SB/GBn2 relaxation. A genuine CVaR-VQE (7 qubits, 3 layers, exact
statevector, parameter-shift gradient) can select a subset over the top 128; its Hamiltonian
is `H = diag(zrank(sorted scores))`, a near-constant rank ladder (S25 L17). Production result:
3.2148 A built chain (`rmsd_arm`) on 126 dev targets; the raw average 3.0483 A is a point
cloud, not a structure, and the two are never compared. Sealed benchmark 60: +0.0103 A paired,
CI [-0.160, +0.180], no validated improvement. The findings are the output.

Read `docs/STATE_BRIEF_2026-09-12.md` in full before anything else. Then the files your lane
brief names. Do not skim; a claim you did not read is a claim you will accidentally re-test
or reverse.

## 1. Absolute rules (Part 11 of the campaign prompt)

1. Never open the sealed benchmark. Never pass `--i-am-spending-the-benchmark`. Never read
   `results/benchmark_manifest.json` contents, benchmark PDBs, or any benchmark RMSD.
2. Never delete, move, or regenerate a pinned artefact: `peptide_folds.json`,
   `peptide_clusters.json`, `catrace_prior.npz`, `results/benchmark_manifest.json`,
   `distogram_models/`, the quarantined `distogram_models_large.STALE-...`, and the `pdbs/`
   file set (never add/remove/rename in `pdbs/`).
3. Never exceed 93% RAM or CPU. Never more than two AMBER/OpenMM jobs at once.
4. Never run an experiment without a pre-registered falsifier written first
   (`s26/PREREG_<name>.md`).
5. Never report an iid CI without the fold-clustered CI beside it. Decide on the clustered one.
6. Never state a number without its artefact path.
7. Never let a native quantity into a non-ORACLE operator. `pytest tests/test_pipeline.py`
   (the NaN-poison test) must pass after every production-path change.
8. Never break a tie by array order (use `s24.stats_lib.argmin_tied`).
9. Never delete a superseded claim. Append the retraction.
10. Never call something a quantum advantage. Never call a small gradient a barren plateau
    without the variance-versus-width analysis. The MPS here is exact; never call it a cause.
11. Never soften a proposal to make it survive. Never move an open item to closed on one
    failed experiment.
12. Never leave yourself idle while your queue is non-empty; take the next pre-registered item.
13. In any deliverable meant for the presenter or the report, never write "genuinely",
    "honestly", "leverage", "robust", "delve", "underscore", or an em dash. Short sentences.
    Real numbers with artefact paths.
14. When in doubt about whether something is allowed, it is not. Log the question in the
    ledger and take the conservative path.

## 2. The standing experimental discipline (earned through 23 retractions)

- Pre-register: `s26/PREREG_<name>.md` with hypothesis, exact falsifier, the comparison it is
  measured against, expected effect size against a computed MDE, memory estimate, agent-hours.
  Written BEFORE the first result exists. Never edited after; append addenda.
- Statistics: `from s24 import stats_lib as ST; r = ST.compare(a, b, folds, names=pdbs,
  label=...)` and print `ST.fmt(r)`. It gives SE, MDE = 2.8016*SE, iid CI, fold-clustered CI,
  W/L/ties, median beside mean, the uniform-effect concentration null, and the verdict.
  Folds: `ST.pinned_folds(pdbs)`. A result must clear its own MDE AND have the fold CI exclude
  zero; 0.7 to 1.3x MDE is the Type-M zone and is not a result; below 0.7x is underpowered.
- Any oracle formed as a minimum over K variants: `ST.best_of_k_within(M)` and quote the
  split-half transfer, never the raw minimum. Report `k_eff`. W/L cannot diagnose it.
- Controls must match the operator's space. A zero-information control must be plausible
  (constant alpha-helix at phi=-63, psi=-42 in RADIANS through `core.geometry.build_backbone`),
  never uniform. Every arm carries a zero-information control and a matched random control.
- Basis discipline: point cloud (`rmsd_avg`, 3.0483) vs built chain (`rmsd_arm`, 3.2148) vs
  repaired emission (`rmsd_full`, 3.2355). State the basis on both sides of every contrast.
  Prefer the built chain; it is the production result.
- ORACLE / ACHIEVABLE / PRODUCTION labelled at every appearance. Any operator that reads a
  native coordinate, distance, torsion, RMSD or derived quantity is ORACLE.
- Provenance: `ST.save_atomic(path, obj, complete_keys=NEED, rows=rows, n_expected=126,
  module_file=__file__)`. `complete: true` is gated on the full key set. Results go in
  `s26/results/`.
- Every positive result is re-run on a different seed and a different fold-assignment order
  (fold labels stay pinned; only the order of processing changes) and must replicate within
  its own CI, or it is not a result.
- Every negative result states its power: was the MDE small enough to see the effect if it
  existed? If not, write "underpowered", not "null".
- The instrument: `s12/instrument.py` (`targets()`, `load_univ(pdb)`, `pool_idx`,
  `shipped_record`, `distogram`, `shipped_score`, `project`, `build_ca`, `selfcheck`). The
  production cache is `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` (126 records, keys
  include `sub`, `avg_ca`, `fit_ca`, `ca`, `amber_ca`, `rmsd_avg`, `rmsd_fit`, `rmsd_arm`,
  `rmsd_full`, `shipped`, `pool_best`). Do not edit either.
- Never edit a module while a job launched from it is running. Every job must checkpoint to
  disk at least every 10 minutes and be resumable (the governor may kill it).

## 3. The resource governor (running; pid in `s26/governor.log`)

`s26/governor.py` samples every 5 s, suspends the newest registered job above 93% RAM or CPU,
kills the newest job after 15 s above 95% (CTRL_BREAK first, then terminate, then requeue),
launches queued jobs when the box has been under 88% for 60 s, and caps AMBER jobs at two.
Baseline load with no campaign work is about 67% RAM (the user's own VS Code and other
sessions); you have roughly 4 GB to work with in total across all lanes.

- Any python process that may exceed 200 MB or run longer than one minute goes through:
  `python s26/jobrun.py --agent <LETTER> --tag <CPU|AMBER|ESM|TEST> --name <unique_name>
  --est-ram <GB> -- python <script> [args]`
  It waits for headroom (and an AMBER slot), registers the job, streams stdout/stderr to
  `s26/logs/<name>.log`, and writes `s26/jobs_done/<name>.json` with exit code, wall time and
  PEAK RSS. Quote that peak RSS as the honest memory estimate in your findings.
- Or queue it: `python s26/enqueue.py --agent <LETTER> --tag <TAG> --name <name> --priority
  <int> --est-ram <GB> -- python <script> ...` and poll `s26/jobs_done/<name>.json`.
- If you cannot estimate memory, run a 1-target probe under jobrun first and read the peak.
- Read `s26/governor_state.json` for the live band before launching anything heavy.
- Tag anything that imports openmm or calls `core.amber` as AMBER. Tag anything that loads
  `esm_cache.npz` or runs ESM-2 as ESM. Never load `esm_cache.npz` (1.5 GB) without a probe;
  `core.pipeline.guard_esm` exists to stop that.
- Never run `pytest tests/` yourself unless your brief assigns it; the Infrastructure lane owns
  the test run.

## 4. Files and formats

- `s26/LEDGER.md`: append-only, one `## L<n> -- TITLE (date, lane)` entry per finding or
  decision. Never edit another lane's entry. Retractions are appended as new entries that name
  the entry they retract.
- `s26/STATUS.md`: every 60 minutes each lane appends two lines under its heading: what it is
  running, what it will do next. Timestamped.
- `s26/agent<Letter>_FINDINGS.md`: your findings, S12 to S25 format (tiers DEMONSTRATED /
  ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN; every number with its artefact path; a
  "what damaged my own expectations" section; a "what I did not do and why" section).
- `s26/PREREG_<name>.md`, `s26/IDEA_<name>.md` (tournament entries: hypothesis, why the
  record does not already close it with the sprint cited, exact falsifier, expected effect vs
  computed MDE, memory estimate, agent-hours), `s26/results/*.json` (provenance-stamped),
  `s26/*.py` scripts named `<letter>_<what>.py`.
- Phase gate: no endpoint experiment (anything that reads an RMSD to a native) may run until
  the coordinator posts "PHASE 0 SIGNED OFF" in `s26/LEDGER.md`. Allowed before that: reading,
  writing PREREGs and IDEAs, writing and unit-testing code on synthetic data, and 1-target
  memory probes under jobrun. Ask the coordinator (finish your turn with the question) if
  unsure.

## 5. Git

- Work on branch `s26`. Commit your own files early and often with descriptive messages.
  End every commit message with:
  `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01GpmP5oMKmcGvwZvxDqQ53t`.
- Only `git add` files you created or were assigned. Never `git add -A`. Never rewrite history.
- Every `.py` on the production path (`core/`, the 24 root modules, `s5/ s7/ s8/ s9/`) must
  remain AST-identical to its parent once docstrings are stripped, unless the change is a
  deliberate, pre-registered, tested behaviour change with a ledger entry. Prefer writing new
  code under `s26/`.
- `*.log` and `*.npz` are gitignored globally; `s26/governor.log` is whitelisted. If you need
  to commit a small npz or log under `s26/`, add an explicit `!s26/<path>` whitelist line and
  say so in the ledger.

## 6. Writing for the presenter

The presenter is sixteen, presenting to a postdoc in variational quantum algorithms, and
must be able to say each edited proposal in under two minutes with every number sourced in
the speaker notes. Plain sentences. No stock words (rule 13). Formal but human.
