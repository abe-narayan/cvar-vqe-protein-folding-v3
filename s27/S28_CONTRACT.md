# S28 CONTRACT (read first; binds every S28 lane)

Sprint 28 of the peptide structure-prediction programme, repository
`C:\Users\abena\Protein-Folding-Algorithm`, branch `s26` (the S26 and S27 work is on it; S28 work
also goes on it, under `s27/`). Python `C:\Users\abena\miniforge_3\python.exe`. The goal: real
movement on accuracy, or a real quantum result, or both. The CVaR-VQE stays the spine: what it
sees, how it is encoded, what it optimises and how its output is consumed are all fair game.

Read in this order before any code: `docs/STATE_BRIEF_2026-09-12.md` (the project in one
sitting), `s26/LANE_CONTRACT.md` sections 1 to 3 (the standing discipline and the governor; every
rule there applies with `s27/` in place of `s26/` for new files), `s27/REPORT.md` sections 1, 4, 6
and 10 (what S27 just closed and the mechanism it found), `s27/PREREG.md` section 1 (what is held
fixed), your own brief in `s27/briefs/`, and the files it names.

## Non-negotiables (instrument integrity)
1. Never open benchmark60: never read `results/benchmark_manifest.json`, `s9/final_report.json`,
   any benchmark PDB, native or RMSD. There is no second benchmark.
2. Never regenerate or move `peptide_folds.json`, `peptide_clusters.json`, `catrace_prior.npz`,
   `distogram_models/`, the `pdbs/` file set, or any pinned artefact.
3. Production code (`core/`, the 24 root modules, `s5/ s7/ s8/ s9/`) is frozen at `a15406c`
   (AST-identical modulo docstrings to `ae86a124` on this branch). New work lives in `s27/` as
   `s27/s28_<lane>_<what>.py`, results in `s27/results/s28_*.json`, tests in
   `tests/test_s28_<lane>.py`.
4. ORACLE arms (anything that reads `nat_ca`, `oracle_rr`, a native torsion or RMSD) are labelled
   ORACLE in the same sentence, every time. No native quantity chooses a parameter, a weight, a
   threshold, a subset or a stopping rule of a deployable arm.
5. The built chain (`s12.instrument.project`, the production projection) is the reporting basis;
   the point cloud is an intermediate. Every RMSD names its basis on both sides of a contrast.
6. Statistics: `s24.stats_lib.compare` (paired per target, SE, MDE = 2.8016 x SE, iid CI beside
   the fold-clustered CI, W/L/ties, concentration null). Decide on the fold-clustered CI. Below
   0.7x MDE is NOT A RESULT; never write "trend". 0.7 to 1.3x is the Type-M zone. Any grid choice
   is priced with `ST.best_of_k_within` and the split-half transfer. Any positive is re-run on a
   second seed and a reversed fold order before it is believed.
7. Pre-register the falsifier before every run (`s27/PREREG_S28_<lane>.md`, appended addenda,
   never edited after a result exists). Try to kill your own positives: matched-random controls in
   the operator's space, permuted channels, order-statistic checks, concentration nulls, and the
   NaN-poison test (replace every native quantity with NaN and show the deployable output is
   bit-identical) for any new operator.
8. One AMBER/OpenMM process at a time. Every python job over 200 MB or one minute goes through
   `python s26/jobrun.py --agent S28<X> --tag CPU|AMBER|ESM|TEST --name <name> --est-ram <GB> --
   python <script>`; probe one target first and quote the peak RSS; checkpoint every 10 minutes;
   the governor (`s26/governor.py`, running) may suspend or kill the newest job above 93% RAM.
9. Never call a small gradient a barren plateau, or its absence; never claim a quantum advantage;
   an MPS is exact here and is never a cause.
10. Re-opening a closed question needs a new angle stated up front in the prereg, citing the
    sprint and ledger line that closed it, not a re-run of the old one.
11. Every number carries its artefact path. No number from memory or prose. Keep the unit suite
    green (`pytest tests/ -q -p no:cacheprovider`, AMBER files as a separate TEST job) and add
    tests for anything new.
12. Ties are never broken by array order (`s27/run_pool.py :: topm` with a stable random key, or
    `ST.argmin_tied`).

## Ledger and status
- `s27/LEDGER.md` is the sprint ledger (append-only; S28 entries are `## S28-L<n> -- TITLE (date,
  lane)`, numbered after the last S28 entry; re-read the tail before appending; a colliding number
  takes a b-suffix). Every entry: question, falsifier, result with `ST.fmt` verbatim where a
  contrast is claimed, verdict, artefact path. Retract openly with a new entry.
- `s27/STATUS.md`: two lines per lane per hour (running / next).
- `s27/s28_<lane>_FINDINGS.md`: your findings in the S12 to S25 format (DEMONSTRATED / ORACLE
  DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN; "what damaged my own expectations"; "what I did not
  do and why").
- Commit your own files early and often on branch `s26` with the trailer lines
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01GpmP5oMKmcGvwZvxDqQ53t`. Never `git add -A`.

## The instrument you reuse (do not re-implement)
`s24.d_harness.Candidates.from_universe(pdb, k=500)` (the pool; `W` (500, n, 3) CA windows,
`PHI/PSI`, `nat_ca` and `oracle_rr` are ORACLE labels), `s25.phys_lib.channels(pdb)` (DIS, LEG,
AMB with cache asserts), `s27/cache/<pdb>.npz` via `s27.run_pool.channels_for(pdb)` (the S27
channels, regenerated in 13 min if absent), `s24.d_harness.arm_vqe` (the genuine CVaR-VQE:
`core.quantum.run_cvar_vqe`, exact statevector, exact parameter-shift gradient),
`s24.d_harness.gate_set_equality`, `s24.d_harness.readout_uniform` (the deployed average),
`s12.instrument.project` (the built chain), `s12.instrument.ca_rmsd`, `s24.stats_lib`.
Production anchors: point cloud 3.048338 (DIS top-75), built chain 3.2126 (rebuild basis) /
3.2148 (production cache), S25 VQE arm 3.0580, random-75 null 3.4209.

## Addendum 1 (2026-09-14 19:10, coordinator; binds every lane from now)
13. The sprint's only currency is a mechanism that changes the FINAL BUILT-CHAIN RMSD. The point
    cloud is an intermediate: any arm at 0.7x MDE or better on the point cloud goes to the built
    chain (`s12.instrument.project`) before it is called anything, and every ledger verdict is
    written on the built chain.
14. No cosmetic variants: no finer grids, extra seeds of a null, or re-parameterisations of the
    convex top-75 average. Every positive's ledger entry answers, in the same entry: (a) does the
    gain survive the averaging/readout bottleneck (S27 REPORT section 6) through the built chain,
    and does it hold on the 108 non-FAIL18 targets or only on FAIL18; (b) what is new relative to
    S26/S27, citing the ledger line it goes beyond; (c) for quantum arms, the three-way split:
    exact objective optimum vs the circuit's reached state vs the emitted structure.
15. No positive is accepted, built on, or reported until lane D has posted STANDS or STANDS WITH
    CAVEAT on it; D attacks within the hour of the entry, starting with the built chain if only
    the point cloud was quoted.
