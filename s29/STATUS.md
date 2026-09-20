# S29 STATUS

## Coordinator
- 2026-09-19 23:35 running: sprint opened; governor v2.5 (band 92-94, hard 95.5); cap 6; six lanes spawning (L, M, T, O, D, X).
- next: read O ladder, T sections 2-3, M map before any build; wave 2 (build) on their word.

## L
- 2026-09-19 23:36 running: brief, contract, charter 4/6/12/18, STATE, S28 report 1/2/8, state-brief 5, S28-L48/L49 and the ham_lib scorer list read; topic 1 (native-free QA) searches starting.
- next: s29/lit/L_1_native_free_qa.md with the equation/assumption/information test per family; ledger S29-L1.
- 2026-09-19 23:54 done: topic 1 posted (S29-L1, `s29/lit/L_1_native_free_qa.md`, 15 papers): NO published QA method is trained or benchmarked below 40-50 residues; AF2's pLDDT has no within-target skill on 588 peptides of 10-40 aa (McDonald 2023) = external confirmation of finding 8; nothing importable; the one actionable item is S8's unfinished free-energy stage (1/24 targets), a probe for D or O.
- next: topic 2 (correlated error in ensembles and the methods that break it) -- directly serves H1 and finding 11.

## M
- 2026-09-19 23:48 running: brief, contract, charter 1/9/13, STATE, state-brief, ARCHITECTURE, s25/QUANTUM.md and the code end to end read (core/data, predict, quantum, project, amber, pipeline; s8/generate; s24/d_harness; s22/qcand_lib (the brief's `s25/qcand_lib.py` does not exist); s12/instrument; s24/stats_lib; s27/run_vqe_chain + run_pool); `s29/DATAPATH.md` posted (13 stages, every function cited by line, counts measured: 183 features, 372,881 MLP params, 6.6k training chains of which 90.5% are the fold-shared fragments, universe 7k-39k windows, DIS anchor 3.2126 reproduced from chain_rows).
- next: `s29/CONVENIENCE_CHOICES.md` (20+ choices with the ledger line that tested each or a cheap decisive test), then the harness audit with commands and a ledger entry.

## T
- 2026-09-19 23:57 running: brief + contract + charter 4/7/11/12 + S28 report + S28-L8b/L11/L18b/L21/L23b/L26b/L35/L36/L41/L43/L48 + QUANTUM.md + S13 locality/Pauli + REPORT V.1-V.10 read; job s29T_spectra launched (spectra + hop-only gradient variance for A, A_c, DD^T at n=4..9, 12 targets).
- next: THEORY.md sections 2 (locally informative class) and 3 (spectral condition) written and committed, ledger entries with their predictions.

## O
- 2026-09-19 23:43 running: brief and record read (S10-5, S24 L2/L3/L13, S23 L6/L9, S22 L4, S28-L1b/L26b/L35/L48); PREREG_S29_O.md committed (1e9bb035); building s29_O_ladder.py + tests.
- next: probe one target under jobrun (peak RSS), point-cloud stage on 126, rung 6 LFO step posted first, chain job 28 structures/target.

## D
- 2026-09-19 23:56 done: cost-RMSD meter built and posted (S29-L2; all four S28 anchors reproduce); prereg checks of O (S29-L3) and X (S29-L4); suite gate 378/3/0 on 17 light files (S29-L5).
- next: attack every positive within the hour; O's rung-6 controls and X's D1/P3 controls when their numbers land; heavy test files in the coordinator's quiet window; hourly S28 reproduction.
- 2026-09-20 00:05 done: S29-L6, the meter's first customer -- lane X's pair log-score is +0.200 ladder rho above the shipped cost on the binding S28 ladder (fold CI [+0.100,+0.343], 5/5) and is the first cost in the record that is uninformative rather than anti-informative there; recognition unchanged (native at the 37.8th percentile, 0.25x MDE from the shipped cost's) and its gradient is undefined on 114/126. Reproduction seed 201 (2FBU / DIS+0.5*LEG_compactness) exact, abs diff 0.0.
- next: lanes M, T, O entries as they post; every lane's proposed cost through the meter on request (`--f module:function`).

## X
- 2026-09-19 23:45 running: brief, contract (S29+S28), charter, STATE, S13 dossier, S15 dossier 1/4/5, S20 L1-L4, S21 L1-L28, STATE_BRIEF 4-5, core/quantum.py, d_harness, instrument, stats_lib, s28_B_hop, D's meter draft read; encoding (ii) chosen (3-mer fragments of the DIS top-8, recombined; q = 3*ceil(n/3) <= 18); prereg being written.
- next: s29/PREREG_S29_X.md (nine questions, new angle vs S13/S15/S20/S21/S28-L21, falsifiers, prior), then s29_X_config.py + tests, then the 12-target probe under jobrun.
