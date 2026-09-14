# REPORT S26: The whole project, from first principles to the presentation

Repository `cvar-vqe-protein-folding-v3`. Written against the working tree on branch `s26` (HEAD `b7f1f280`; the prompt names commit `ae86a124`, which is in the history of this branch; see Appendix D, item D1). Production code is frozen at commit `a15406c` (`ARCHITECTURE.md`, header). Every number in this document carries the path of the artefact it was read from. Every number is also listed in Appendix B. Numbers that do not reproduce, or that exist only in prose, are in Appendix D and are not used in the body.

Companion: `docs/REPORT_S26_SUMMARY.md` (one page).

---

## To the presenter, before anything else

You are going to present a project whose most important results are negative, and whose most-quoted early numbers were later corrected. That is not a weakness to hide. It is the substance of what you built: an instrument that can tell a real effect from noise on a 126-peptide test set, and a record of thirty or so ideas that instrument killed. Read this document once, end to end. Then read Part X twice.

Before you trust any number in it, know the following. I reproduced the three headline numbers from artefacts by reading them, as the rules required. The built-chain mean 3.2148 Å lives in `bench_results/compare_tuning126.json` (key `science.rmsd_arm`, value 3.214765154210998). The point-cloud mean 3.0483 Å lives in the same file (key `science.rmsd_avg`, value 3.048338093879532). The per-target 0.182 Å for T030 (PDB 1S9Z) lives in `results/summary/results.csv` line 913 (production row, `rmsd` 0.18242222359245708). All three are where the state brief says they are. But there are discrepancies you must know about, and they are listed in Appendix D: the prompt's commit hash is not the working tree's HEAD (D1); the results-lab leaderboard prints the production mean as 3.2126 while the production cache prints 3.2148, and the results-lab's own docstring says 3.2126 "does not reproduce anywhere" (D2); the T030 figure is 0.18242 in the results lab and 0.18198 in the production cache, both rounding to 0.18 (D3); the slide deck `vqe_research_overview.pptx` is not on disk or in git, so Part X.3 annotates the eleven-slide structure the campaign brief describes, not a file I opened (D4); the `quantum_stage` docstring in `core/pipeline.py` still asserts the withdrawn "+0.113 Å" claim (D5); the "+0.370 rank correlation" the prompt lists as withdrawn does not exist under that value in any artefact I could find, the closest being s14's +0.355 corrected to +0.161 (D6); the machine's core-equivalent figure is 6.43 in one place and 6.22 in another (D7); and several numbers the S26 examination found to be document-only, plus numbers quoted from memory notes, are listed with their status (D8 to D12).

Everything else in this document is sourced. When you are asked a question you cannot answer, the honest sentence is "I do not know; the artefact that would answer it is X, and I can check." That sentence is in Part X.6, and it is a good one.

---

# PART I. THE PROBLEM

## I.1 What a protein is, in the terms this project uses

A **protein** is a chain of small molecules called **amino acids**, joined end to end. There are twenty kinds of amino acid. A protein's **sequence** is the list of which kind sits at each position, written as a string of capital letters, one letter per amino acid (for example `SIRELEARIRELELRI`, the sequence of target 1S9Z; `results/summary/results.csv` line 913). Each amino acid in the chain is called a **residue**. The chain has a repeating **backbone** of three atoms per residue, named N, CA and C. The CA atom is the **alpha carbon** (written Cα); it is the atom the side group hangs from, and it is the atom this project predicts. The other atoms of each residue, the part that differs between the twenty kinds, are the **side chain**. Two consecutive CA atoms are always about 3.8 Å apart, because the atoms between them are rigid. An **Ångström** (Å) is 10⁻¹⁰ m; a hydrogen atom is about 1 Å across.

A **peptide** is a short protein. This project's targets are 9 to 16 residues long (`core/data.py`, `benchmark()`, `len 9-16`).

In water, a protein chain does not stay straight. It curls up into a particular three-dimensional shape, and the same sequence curls up the same way each time. That shape is the **native structure** (or simply the **native**). **Protein folding** is the process of reaching it; the **protein folding problem** is predicting the native shape from the sequence alone. Native shapes are measured experimentally and deposited in the **Protein Data Bank** (PDB), a public database, each with a four-character code such as `1S9Z`. For peptides, the experiment is usually **NMR** (nuclear magnetic resonance spectroscopy), which yields not one structure but a bundle of similar ones; this project uses the first model of the bundle as the native.

Two named local shapes recur. An **alpha helix** is a spiral, 3.6 residues per turn, in which each residue's backbone forms a hydrogen bond (a weak attraction between an N-H group and an O atom) to the residue four positions earlier. A **beta sheet** is made of stretches of chain lying side by side, hydrogen-bonded across. Everything else is called **coil** or **loop**. Each residue's backbone has two rotatable bonds, and the two rotation angles are called **phi** (φ) and **psi** (ψ), the **torsion angles**. The third backbone angle, **omega** (ω), is almost always 180° (a **trans** peptide bond); a **cis** peptide bond (ω near 0°) is rare. A plot of φ against ψ is a **Ramachandran plot**; most real residues fall in two or three dense regions of it. A **distance map** is the table of all pairwise CA to CA distances; a **contact** is a pair of residues closer than some cutoff.

## I.2 How accuracy is measured

The distance between a predicted structure and the native is the **RMSD**: the root mean square of the distances between corresponding CA atoms, after the two structures have been rigidly superposed to minimise that quantity (the **Kabsch** algorithm; `core/geometry.py`, `kabsch_rmsd_batch`, line 1012). RMSD is in Å. Lower is better. 0 Å means identical. For a 9 to 16 residue peptide, under 2 Å is close to the experimental resolution of the NMR bundle itself; 3 Å is "the right general shape"; above 5 Å is wrong.

Two conventions matter, and this project got them confused once and then fixed it (S25 L4, L8, L11; `s25/LEDGER.md`). A **point cloud** is a set of CA positions that need not obey the 3.8 Å neighbour spacing; an average of several real structures is a point cloud. A **built chain** is a set of CA positions that does obey it. Only a built chain can be turned into a protein. The project's primary reporting basis is the built chain (`results/summary/leaderboard.csv`, `basis` column `built_chain_bb`); the point cloud is reported as secondary and is marked "NON-PHYSICAL" in the build log (`s25/results/frozen_build.log`).

## I.3 Why short peptides are hard

Long proteins have many long-range contacts, and modern predictors (AlphaFold2, ESMFold) learn from the patterns those contacts leave in the record of related sequences. A 12-residue peptide has few contacts, few relatives, and is often flexible in solution. Sequence carries less information about structure at this length. The project measured this directly: the best-matching window in its library has only about 12% sequence identity to the target it matches (`docs/FINDINGS.md` §9, "Structure and sequence are decoupled"); φ carries no sequence signal at peptide length, full sequence context predicting φ to 36.1° against 36.4° for a sequence-blind marginal (`s13/cache/tors_rows.npz`, re-derived in S26, `s26/LEDGER.md` L31 and L32; see Appendix D, D8). So the problem the project set itself is hard for reasons that have nothing to do with computing power.

## I.4 What the project does, in one paragraph

Given a peptide sequence, the pipeline retrieves 500 candidate CA traces from a library of real protein and peptide fragments of the same length, scores each candidate against a learned prediction of the target's distance map, keeps the 75 best-scoring, averages them in a common frame, and then projects that average onto a chain with correct 3.8 Å spacing and realistic torsion angles. An optional step relaxes the chain in a physics force field. A variational quantum circuit (a CVaR-VQE on 7 qubits) can be switched on to select the subset instead of a plain "top 75" cut; it is off in the shipped configuration (`core/pipeline.py`, `Config`, `quantum=False`). Parts III to V describe every stage.

## I.5 The headline, correctly stated

Over the 126-target tuning instrument (Part II), the shipped pipeline emits built chains with mean CA RMSD **3.2148 Å** (`bench_results/compare_tuning126.json`, `science.rmsd_arm`) and, on the secondary point-cloud basis, **3.0483 Å** (same file, `science.rmsd_avg`). The results-lab leaderboard, which rebuilds the chains through the same projection, prints 3.2126 Å for the production configuration, median 2.9661 Å, best 0.18242 Å (T030, 1S9Z), worst 8.234 Å (T057, 2MQ2), 28.57% of targets under 2 Å and 50.79% under 3 Å (`results/summary/leaderboard.csv`, production row; the 3.2126 vs 3.2148 difference is Appendix D, D2). The three targets shown as overlays in the presentation are 1S9Z at 0.18 Å (0.18242, `results/summary/results.csv` T030 production row), 2MJQ at 0.40 Å (0.39949916705333455, T053 production row), and 6WPB at 0.42 Å (0.4167616469202071, T098 production row).

This is not a sub-2 Å method. It is a 3.2 Å method whose good half is under 3 Å. The best a perfect chooser could do from the same 500 candidates is 1.7108 Å (`bench_results/compare_tuning126.json`, `science.pool_best`; this is an ORACLE arm: it uses the native to choose and is not achievable). The best a perfect chooser could do inside the 75 the pipeline keeps is 2.3062 Å (same file, `science.top_m_best`; ORACLE, uses the native, not achievable). The gap between 3.2148 and the ORACLE 2.3062 is the selection problem; the gap between the ORACLE 2.3062 and the ORACLE 1.7108 is the shortlist problem; the gap between the ORACLE 1.7108 and 0 is the library and the distance prior. Most of the sprints in Part VI were spent learning that none of these gaps can be closed by anything the project can compute without the native.

## I.6 Why quantum, and what the quantum component turned out to be

The founding idea was that choosing a subset of candidates that are jointly consistent is a combinatorial problem, and that a variational quantum circuit trained on a **CVaR** objective (Conditional Value at Risk: the mean of the lowest α fraction of the energy distribution; Part V.4) might explore that choice better than a greedy cut. The project built that component, verified it to machine precision, and then measured what it does. What it does is reproduce the classical ranking: the set of candidates the trained circuit puts weight on is always a subset of the top of the score's own ordering (the set-equality theorem, verified on 2,592 cells with 0 violations; `s25/results/q_verify.json`). The circuit's contribution to accuracy, measured on the instrument where it was tuned, is −0.1405 Å against the plain argmin at 0.68× the minimum detectable effect, which under the project's fixed rule is not a result (`s25/results/q_alpha.json`, `vs_no_circuit`, "VQE_LFO - argmin (shipped)": mean 3.3135 vs 3.4540, MDE 0.2051, verdict UNDERPOWERED). The quantum work that survives is a set of exact statements about trainability and about how a force field's structure reaches a circuit's gradients (Part V.8 to V.10, Part VII). That is what to present to a VQA group.

---

# PART II. THE INSTRUMENT

The instrument is the thing that separates this project from a demo. It is the set of targets, the cross-validation folds, the statistics, and the rules that decide what counts as a result. Everything else was measured on it.

## II.1 The target sets

**tuning126.** 126 peptides of 9 to 16 residues, chosen to be cluster-disjoint from each other and from the training corpus (a **cluster** here is a group of sequences above a similarity threshold; `core/data.py`, `clusters()`, cached in `peptide_clusters.json`). They are ordered by PDB code (`s12/instrument.py`, `targets()`). Every result in this document, unless it says otherwise, is a mean over these 126. The production configuration key is `1fc9f2dcf489e2fb` (`s12/instrument.py`, `PROD_KEY`; the per-target cache is `bench_results/cache/1fc9f2dcf489e2fb/<PDB>.json`, e.g. `s26/results/e_trace_1S9Z.json`, `target.record_path`).

**Five folds.** The 126 are split into five folds, pinned on disk (`core/data.py`, `folds()`, `peptide_folds.json`). Each fold has its own distance-map model trained with that fold's targets held out (**leave-fold-out**; `core/predict.py`, five models). This matters for the statistics (II.3).

**benchmark60.** A sealed set of 60 targets in `results/benchmark_manifest.json` (sha256 `a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d`, 8002 bytes; `core/data.py`, `benchmark()`, `BENCHMARK_N=60`). It was spent once, in Sprint 9 to 10, and the command line refuses to run on it without the flag `--i-am-spending-the-benchmark` (`core/bench.py`, line 1068 to 1074). The one result: the tuned pipeline scored 2.9610 Å against the baseline's 2.9507 Å on the 60, a difference of +0.0103 Å with 95% CI [−0.1596, +0.1803], 31 wins to 29 losses, SE 0.0867 (`s9/final_report.json`, cited through `verify/headline_audit.json`, `benchmark60`, source "s9/final_report.json (cached, NOT re-run)"). The tuning-set gain of −0.250 Å [−0.383, −0.117] did not transfer (`docs/FINDINGS.md`, S9-10). This report does not open the benchmark and neither should you.

**dev24.** A 24-target development set disjoint from the benchmark (`core/data.py`, `dev_set(n=24, seed=7)`), used when an experiment had to be judged before spending anything.

**FAIL18.** The 18 targets on which the pipeline is worst, kept as a named list for stratified reporting: 1ID6, 1JBF, 1LB7, 2BFI, 2BP4, 2JN5, 2MQ2, 2N5C, 2NB7, 2NDM, 3BTB, 3SGO, 5W52, 7JS6, 7LCW, 8T63, 9KAR, 9L1M (`s12/instrument.py`, `FAIL18`). Ten of the eighteen are amyloid fibrils (peptides whose native form is a stacked, repeating aggregate rather than a single folded chain) or lasso peptides (a chain threaded through a ring formed by its own first residues), a class enrichment with Fisher p = 1.2e-6 (`s12/SPRINT12_DOSSIER.md`).

## II.2 What a "result" is

The project's fixed rule, applied to every comparison since Sprint 21:

1. Every comparison is **paired**: the same 126 targets, arm A minus arm B per target.
2. The **SE** (standard error) of the paired mean is reported beside every mean.
3. The **MDE** (minimum detectable effect) is 2.8016 × SE, per comparison (`s21/LEDGER.md` L8; the constant is the sum of the 1.96 and 0.84 normal quantiles, 80% power at 5% two-sided). An earlier project-wide constant of 0.084 Å was wrong by a factor between 0.11× and 9.54× across real comparisons, a range of 84× (`s21/LEDGER.md` L8).
4. **Verdicts**: |effect| > MDE with the sign favouring A is BETTER; with the other sign is WORSE; |effect| ≤ MDE is NOT MEASURED; below 0.7× MDE the arm is UNDERPOWERED and "NOT A RESULT" (`s25/results/q_verify.json`, "UNDERPOWERED (<0.7x MDE) rule present: True"; example verdict strings throughout `s25/results/phys_suite.json`).
5. **Two confidence intervals**: an i.i.d. paired bootstrap and a **fold-clustered** one. The fold-clustered interval is required because each fold's distance model saw about 100 of the other 125 natives during training, so the 126 per-target errors are not independent; the i.i.d. interval is anticonservative (`s20/LEDGER.md` L-B; `docs/STATE_BRIEF_2026-09-12.md` §6.4).
6. **Wins/losses** (W/L) are printed but never used as the verdict. A W/L split cannot diagnose an order-statistic artefact (`s24/LEDGER.md` L15).
7. **Concentration check.** The mean after dropping the ten largest per-target gains is compared to the same statistic under a uniform-effect null, with percentiles. A raw drop-top threshold is not a test; the null is (`docs/FINDINGS.md`, S9-10; `s25/results/phys_form5.json`, `concentration` blocks with `null_p10/p50/p90`).
8. **Type-M flag.** When power is below about 0.8 the magnitude of a "significant" effect is inflated; the runner prints the inflation factor and flags it (`s25/results/phys_suite.json`, `type_m`, `type_m_flag`).

## II.3 ORACLE arms and controls

An **ORACLE** arm is one that uses the native to make a choice. It answers "how good could this stage be" and is never a result. This document labels every ORACLE number as ORACLE in the same sentence. The main ones, as a ladder (`docs/STATE_BRIEF_2026-09-12.md` §2, each with its own artefact in Appendix B):

- annealing on the Legacy energy alone, 4.624 Å (not ORACLE; a control);
- random 75 from the 500, 3.4251 Å (`s25/results/phys_suite.json`, `random_null_rank.mean`; a control);
- constant alpha helix for every target, 4.065 Å (`s13/SPRINT13_DOSSIER.md`; the zero-information control);
- sequence-only torsion predictor, 3.770 Å (`s13/SPRINT13_DOSSIER.md`);
- shipped pipeline, 3.2148 Å built chain (`bench_results/compare_tuning126.json`);
- ORACLE perfect prior fed to the same pipeline, 2.2261 Å (ORACLE, uses the native; `s24/LEDGER.md` L13);
- ORACLE pool best, 1.7108 Å (ORACLE, uses the native; `bench_results/compare_tuning126.json`, `science.pool_best`);
- ORACLE universe best over all windows, 1.313 Å (ORACLE, uses the native; `docs/STATE_BRIEF_2026-09-12.md` §2);
- ORACLE torsion ceiling, 1.594 Å (ORACLE, uses the native; `s13/SPRINT13_DOSSIER.md`);
- ORACLE distance geometry from the true distance map, 0.611 Å (ORACLE, uses the native; `s15/LEDGER.md`).

Three lessons about controls that the project learned by getting them wrong (`s16/LEDGER.md` L27; `s20/LEDGER.md`; memory `control-must-match-the-operators-space`):

- A control must live in the operator's own space. "Random 75 from the 500" is the control for the top-75 cut; "random windows from the library" is not.
- A zero-information control must be plausible, not uniform. A uniform draw on the torsion torus is a worse structure, not an uninformative one; the constant alpha helix (4.065 Å) beats it by 0.457 Å (helix 4.065 vs random 4.522, `s13/SPRINT13_DOSSIER.md`).
- A best-of-K over correlated arms is an order statistic, not a signal. Any per-target minimum over K variants is mostly best-of-K (`s24/LEDGER.md` L15, L15-A).

## II.4 The leakage guard

The native is a reporting label. The only function that opens it is `label()` in `core/pipeline.py`. Tests poison the native coordinates with NaN and confirm every emitted quantity is bit-identical (`verify/leak_audit.json`: ALL_CLEAN True, 19 quantities on each of 1CB3, 1CEK, 1CS9, n_bit_identical 19 of 19). A declared, unfixed defect stands: the sequence-identity function that built the clusters normalises by the longer sequence, which is leaky at the member level (`core/data.py`, `identity()`, docstring; `s24/LEDGER.md` L4). Measured price on the tuning set: 13 of 126 targets carry a K=500 pool window at ≥0.6 identity, four at exactly 1.0, and removing them moves the mean by +0.0004 Å (`docs/FINDINGS.md`, corrections table, S10-4; the per-target artefact is in git history at `5fa05cd`, `s26/LEDGER.md` L31). The benchmark equivalent is 13 of 60 at ≥0.6 with 8Y3S and 8ZG3 at 1.0, price +0.0030 Å (`docs/FINDINGS.md` S9-10). The S26 fix (`norm="shorter"`) ships dark, because changing the clustering would move the pinned folds (`s26/LEDGER.md` L15).

## II.5 Reproduction

The whole 126-target production run reproduces bit-for-bit between one worker and eight (`bench_results/compare_tuning126.json`: `science_delta` 0 or ±4.4e-16 on every arm; wall 2537.568 s vs 303.137 s, speed-up 8.371). The S26 examination re-ran it fresh and found max per-target disagreement 0.0 on all four bases (`s26/LEDGER.md` L28; `s26/results/e_reproduce.json`). The projection stage is bit-identical on 126 of 126 (`verify/project_exactness.json`, `agg.n_bit_identical` 126, `worst_abs_diff` 0.0). The unit suite is 370 tests, 357 passed, 13 skipped, 0 failed (`s26/results/test_run.json`, `s26/TEST_RUN.md`; the 13 skips are 11 opt-in slow tests and 2 absent artefacts).


---

# PART III. THE PIPELINE, STAGE BY STAGE

The production pipeline is `core/pipeline.py`, `run_target(cfg, seq, fold)`. Its configuration is the `Config` dataclass with defaults `k=500, m=75, penalty="ramah", lam=0.3, amber_k=10.0, amber_steps=0, min_sep=2, maxiter=300, multi_start=True, project_grad="exact", tie_break="stable", reference_precision=True, quantum=False, legacy=False, vqe_qubits=7, vqe_layers=3, vqe_iters=50, vqe_seed=0, legacy_top=5` (`core/pipeline.py`, `Config`). Two targets are traced below at every stage: **T030 / 1S9Z** (`SIRELEARIRELELRI`, 16 residues, fold 1; the best target, 0.18 Å) and **T122 / 9KAR** (`GGWGTVPDWFFNMNW`, 15 residues, fold 3; 7.4336 Å built chain, a FAIL18 member). Both traces were recorded fresh by the S26 examination and match the production cache at 0.0 on every stage (`s26/results/e_trace_1S9Z.json`, `s26/results/e_trace_9KAR.json`, `reporting.*.diff` 0.0).

## Stage 1. Sequence features (ESM-2)

The sequence is embedded with **ESM-2**, a protein language model (a neural network trained to fill in masked letters in millions of protein sequences; its internal representation carries structural information). The 650M-parameter model's per-residue embedding is reduced to 32 principal components (**pca32**) and its contact-prediction head is read as a contact map. The embeddings are cached (`esm_small.npz`; `s26/results/e_trace_9KAR.json`, `stages.esm.hot_cache`). ESM matters: on selection, ESM pca32 buys −0.288 Å [−0.484, −0.092] over one-hot encoding, p = 0.0046, 74W/45L (`docs/FINDINGS.md`, S7-11; the per-target artefact `s7/repr_tune.json` is lost, `s26/LEDGER.md` L11, so the number is document-and-summary only; Appendix D, D9).

## Stage 2. The distance prior (distogram)

A **distogram** is a predicted probability distribution over the distance between every pair of residues. The project's distogram is a small multilayer perceptron per fold (five models, `distogram_models/fold<k>_esm_frag.pt`; `s26/results/e_trace_9KAR.json`, `stages.distogram.model_file`), trained on 787 peptides plus fragments, with the target's own fold held out. It predicts, for each pair (i, j) with |i − j| ≥ 2, a 17-bin probability vector over distance. The bin edges are 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.5, 14.0, 16.0, 19.0, 23.0 Å, and the bin centres are 4.0, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.5, 9.5, 10.5, 11.75, 13.25, 15.0, 17.5, 21.0, 25.0 Å (`core/predict.py`, `BIN_EDGES`, `CENTRES`, `NBINS=17`). The last gap is 4.0 Å wide, which is why the S25 quantisation finding (L6) matters: the per-pair target the score consumes takes only 17 values (`s25/LEDGER.md` L6).

Each pair's posterior also yields an expected distance and a standard deviation `sd`. The prior is over-confident by a factor of two: the standardised residual (true minus expected, over sd) has standard deviation 1.9962 averaged over targets, the 90% interval covers only 61.6% of true distances, and 24.1% of pairs have multimodal posteriors (`s25/results/calib.json`, 126 rows; means computed over the `z_sd`, `cov90`, `multimodal_frac` columns; `s25/LEDGER.md` L1). Calibrating it does not help RMSD, because the score does not use the width (`s25/LEDGER.md` L2, with L7 retracting L2's mechanism claim but not the null).

For 1S9Z the expected distances range 5.05 to 21.07 Å with sd 0.23 to 2.70 Å (`s26/results/e_trace_1S9Z.json`, `stages.distogram`). For 9KAR they range 5.62 to 23.44 Å with sd 0.23 to 7.58 Å (`s26/results/e_trace_9KAR.json`); the prior is much less sure about 9KAR.

## Stage 3. Retrieval (BLOSUM62, K = 500)

The library is every contiguous window of the target's length cut from a bank of real structures: 787 peptides plus 6,003 protein fragments from 1,001 deposits (`README.md`; `docs/STATE_BRIEF_2026-09-12.md` §4). For a 16-mer the universe has 7,193 windows; for a 15-mer, 9,814 (`s26/results/e_trace_1S9Z.json` and `e_trace_9KAR.json`, `stages.retrieve.n_windows`). Each window's sequence is aligned without gaps to the target and scored with **BLOSUM62** (a standard 20 × 20 table of amino-acid substitution scores). The 500 highest-scoring windows form the **pool** (`core/pipeline.py`, `retrieve`). Ties are broken by a stable argsort, so the order of `pdbs/` on disk feeds the tie order (`Config.tie_break="stable"`; the tie-breaking trap is II.3 and S8-8). For 1S9Z the BLOSUM scores in the pool span −1 to 19 and 144 of the 500 windows come from peptides; for 9KAR, −6 to 32 and 160 from peptides (`stages.retrieve.sim_min/sim_max/n_from_peptides` in the two trace files).

Two facts about retrieval that took sprints to establish. BLOSUM retrieval beats random windows at every K on a matched universe (`docs/FINDINGS.md`, S7-12), so retrieval is not the weak stage. But sequence is a poor key at this length: the pool's best member sits at 1.7108 Å on average (ORACLE, uses the native; `bench_results/compare_tuning126.json`), and the pool's error is 68% common-mode, meaning that most of what is wrong about the 500 candidates is wrong in the same direction, so averaging can remove at most the other 32% (`s23/LEDGER.md` L9: common-mode fraction 0.676, 50.7× the i.i.d. prediction).

## Stage 4. Scoring (L1 Bayes risk against the distogram)

Each candidate window w is a CA trace. Its pairwise distances d_p(w) are compared to the distogram with the score

    s(w) = Σ_p w_p Σ_c p_c |d_p(w) − C_c|,

the expected absolute error under the posterior, summed over pairs, where C_c are the 17 bin centres, p_c the posterior, and w_p a per-pair weight `w_p = shell / (sd_p + 0.5)^g` with `shell ≡ 1` and `g = 1` because `score_weights.json` is absent (`core/predict.py`, `Distogram.__init__`, `self.w = shell/np.maximum((self.sd+0.5)**g, 1e-6)`; `_score_weights()` returns ones). The risk is evaluated on a grid `np.arange(2.0, 40.0, 0.05)` and stored as float32 (`core/predict.py`, `_risk`). The sd-derived weights are informative: inside the S18 refinement objective, replacing them with flat weights costs +0.153 Å [+0.066, +0.247], 50W/76L, native-free (`s18/results/q_ceilattack2.json`; `s18/LEDGER.md` L10). Lower score is better. For 1S9Z the 500 scores span 0.693 to 7.269 with 488 distinct values; for 9KAR, 2.242 to 5.179 with 483 distinct (`stages.score` in the two trace files).

## Stage 5. The cut (top-75) and the optional quantum selection

The 75 lowest-scoring candidates are kept (`Config.m=75`; `core/pipeline.py`, `filter_pool`). This is the whole of selection in the shipped configuration: `selector` in the results table reads "score-filter top-75 (== CVaR tail, s24 theorem)" (`results/summary/results.csv`, every production row). When the quantum stage is on, the cut is widened to 128 so that a 7-qubit register can index it, the scores are rank-standardised into an energy vector E (Stage 6), a CVaR-VQE is trained (Part V), and its output distribution p either selects a single medoid (`consensus_medoid(block, p)`) or weights the average (`average_weighted`). In the fresh traces both `sel` and `sel_uniform` are 449 for both targets, i.e. the quantum readout picks the same member the uniform readout picks (`stages.quantum.sel`, `sel_uniform` in both trace files).

Why 75: the operator that follows consumes the set mean, not the set best (`s12/SPRINT12_DOSSIER.md`: d_out = 1.16 × mean + 0.04 × best, R² 0.893 on narrow sets). The per-target optimal m varies and a per-target ORACLE m is worth −0.239 Å (ORACLE, uses the native; `s22/LEDGER.md`), but no native-free router finds it (four routers null, `s22/LEDGER.md`).

For 9KAR this stage is where the target is lost. The pool's best member is at 1.3997 Å from the native, but the best member inside the top-75 is at 5.3436 Å (`s26/results/e_trace_9KAR.json`, `reporting.record_pool_best`, `record_top_m_best`; both ORACLE quantities, they use the native). The score, trained on a prior with sd up to 7.58 Å for this target, ranks the good candidates out. For 1S9Z the pool best and the top-75 best are the same window at 0.1596 Å (`e_trace_1S9Z.json`, `reporting`).

## Stage 6. Energies and the Hamiltonian (only when a selector is on)

When the quantum or Legacy selectors are on, each of the 128 candidates gets an energy. In the shipped quantum path the energy is the **rank-standardised score**: `E = _zrank(scores)`, a ladder from −1.7186 to +1.7186 in 128 equal steps (`core/pipeline.py`, `_zrank`; `s26/results/e_trace_1S9Z.json`, `stages.hamiltonian.E_first5`, `E_last5`, `E_range` 3.4372). Ties in the score make small deviations from the ideal ladder: 0.394% of the range for both traces, and at worst 1.18% over eight targets in S25 (`s25/results/q_gibbs.json`, `spectrum_target_independence.worst_frac_of_range` 0.0118). This is S25 L17: the Hamiltonian barely changes between targets, and Part V.7 explains why that closes several questions at once. The Hamiltonian is diagonal, `H = diag(E)`, so the problem posed to the circuit is: put probability mass on the low-E basis states.

The two physics energies (Legacy and AMBER) can replace or join the score here; Part IV covers them. They are off in production.

## Stage 7. Medoid frame and coordinate average

The 75 (or 128) candidate traces are superposed onto a common frame. The frame is the **medoid**: the member with the smallest summed RMSD to the others (`core/pipeline.py`, `consensus_medoid`). Each member is Kabsch-superposed onto the medoid and the CA coordinates are averaged. For 1S9Z the medoid is local index 44, pool index 449; for 9KAR, local 34, pool 449 (`stages.average.medoid_local`, `medoid_pool_index`). The average is a **point cloud**. Its neighbour spacings are wrong: 3.75 Å mean and 3.36 Å minimum for 1S9Z, 1.97 Å mean, 1.16 Å minimum for 9KAR (`stages.average.bond_mean`, `bond_min`), against 3.80 Å for real chains. Averaging contracts the backbone 22.3% on average (`docs/STATE_BRIEF_2026-09-12.md` §2; `s16/LEDGER.md`). This object scores 3.0483 Å over 126 targets (`bench_results/compare_tuning126.json`, `science.rmsd_avg`). It is labelled `"raw average (illegal)"` in the code (`core/bench.py`, `ARM_NAMES`, line 682). It is the number the project quoted as its result for many sprints, and S25 L4/L8 decided it is not a result but an intermediate.

## Stage 8. Projection onto a built chain

The point cloud is projected onto the nearest chain with ideal geometry: every CA to CA virtual bond 3.804 Å, and torsion angles penalised toward the allowed Ramachandran regions with penalty `ramah` at weight `lam = 0.3` (`Config.penalty="ramah", lam=0.3`; `core/project.py`). The optimisation runs to 300 iterations, from multiple starts (`Config.maxiter=300, multi_start=True`), with exact gradients (`project_grad="exact"`; `core/project.py` docstring: under the analytic mode 126/126 targets move, median 0.042 Å, worst 1.90 Å, so `exact` is the default). The output for 1S9Z has bond mean 3.8039549 Å, min and max equal to 7 decimals (`e_trace_1S9Z.json`, `stages.projection`). This is the **built chain**, the object that scores **3.2148 Å** (`bench_results/compare_tuning126.json`, `science.rmsd_arm`). The projection costs +0.1664 Å against the point cloud, 3.26× its MDE (`s25/LEDGER.md` L8 and L11, `basis_delta` in `results/summary/leaderboard.csv` is 0.1643 on the results-lab rebuild). For 1S9Z the chain lands at 0.18198 Å (`e_trace_1S9Z.json`, `reporting.rmsd_arm`) and the point cloud was 0.19589; for 9KAR the chain lands at 7.4377 Å and the point cloud was 7.0610 (`e_trace_9KAR.json`, `reporting`). The results-lab re-projection (`s25/resultslab/exportlib.py`, `chain_from_ca`, `PRODUCTION_LAM=0.3`) gives 0.18242 for 1S9Z and 7.4336 for 9KAR (`results/summary/results.csv`); the difference between the two routes is Appendix D, D2 and D3.

## Stage 9. Optional restrained AMBER relaxation, then labelling

If `amber=True`, the built chain gets side chains and hydrogens and is minimised in the AMBER ff14SB force field with the GBn2 implicit solvent model through OpenMM, with the backbone atoms N, CA, C restrained by a harmonic spring of k = 10 kcal/mol/Å² and zero dynamics steps (`Config.amber_k=10.0, amber_steps=0`; `core/amber.py`, `RESTRAINED_BACKBONE=("N","CA","C")` line 808, `refine_coords` line 1497). The energy falls from 5,370 to −1,291 kcal/mol for 1S9Z, moving the CA trace 0.191 Å; for 9KAR from 6.0e12 (a clash) to 1,262 kcal/mol, moving 0.610 Å (`stages.amber.record_e0`, `record_e1`, `record_moved` in the two trace files). Over 126 targets the relaxed chain scores 3.2355 Å (`bench_results/compare_tuning126.json`, `science.rmsd_full`), worse than 3.2148 by +0.021. S26 measured that the relaxation's displacement is worse than a random move of its own size (+0.0111, fold CI [+0.0062, +0.0171]) and points away from the native (cos −0.049): it is a validity step, not an accuracy step (`s26/LEDGER.md` L39).

Finally `label()` opens the native and computes RMSDs on every basis (`core/pipeline.py`, `label`). That is the only place the native is read.

## III.10 The two traces side by side

| stage | 1S9Z (T030) | 9KAR (T122) | source |
|---|---|---|---|
| windows in universe | 7,193 | 9,814 | `s26/results/e_trace_*.json`, `stages.retrieve.n_windows` |
| BLOSUM range in pool | −1 to 19 | −6 to 32 | `stages.retrieve.sim_min/max` |
| prior sd range | 0.23 to 2.70 Å | 0.23 to 7.58 Å | `stages.distogram.sd_min/max` |
| score range | 0.693 to 7.269 | 2.242 to 5.179 | `stages.score.sc_min/max` |
| pool best (ORACLE, uses native) | 0.1596 Å | 1.3997 Å | `reporting.record_pool_best` |
| top-75 best (ORACLE, uses native) | 0.1596 Å | 5.3436 Å | `reporting.record_top_m_best` |
| point cloud | 0.1959 Å | 7.0610 Å | `reporting.rmsd_avg` |
| built chain | 0.1820 Å | 7.4377 Å | `reporting.rmsd_arm` |
| after AMBER | 0.3174 Å | 7.5437 Å | `reporting.rmsd_full` |
| wall time (fresh trace) | 13.13 s | 13.21 s | `wall_s` |

The story of the project is in the fourth and fifth rows of the 9KAR column. A good answer was in the pool. The score threw it away. Nothing downstream can get it back.


---

# PART IV. THE TWO PHYSICS HAMILTONIANS

A **force field** is a function that assigns an energy to a set of atom positions, built from terms for bond stretching, angle bending, torsion rotation, electrostatics, and van der Waals contact. Nature's native structure is, to a first approximation, the lowest-free-energy structure, so a force field ought to rank the native lowest among candidates. This project tested two energies in that role. Neither ranks nativeness on its pools. The measurement that settled it is the seven-configuration suite (`s25/results/phys_suite.json`; `s25/LEDGER.md` L16), described in IV.4.

## IV.1 Legacy: an eleven-term coarse energy

"Legacy" is the project's own CA-level energy (`core/energy.py`). Its eleven terms and default weights are steric 4.0, contact 1.0, hbond_local 1.0, hbond_longrange 3.0, coop_helix 2.0, coop_sheet 2.0, solvation 0.5, electrostatic 1.0, aromatic 0.8, torsion 0.15, compactness 0.4 (`core/energy.py`, `DEFAULT_WEIGHTS`, line 244; `FITTED_WEIGHTS = dict(DEFAULT_WEIGHTS)`, line 1038, so nothing was fitted). What it measures, once the mechanism was isolated, is compactness: against 200 matched-random partitions, Legacy-preferred candidates are 0.45 Å more compact in radius of gyration (the root-mean-square distance of the CA atoms from their centre) (−0.4477 [−0.5110, −0.3911], 124W/2L, 5 of 5 folds; `s20/LEDGER.md`, D lane; `s20/agentD_FINDINGS.md`), and AMBER-preferred candidates are expanded with open sterics. Its certified optimum over torsion space (found by exhaustive or provably optimal search) is 0.139 Å worse than a random structure, where a structural objective's certified optimum is 0.885 Å better than random (`s14/LEDGER.md`). Annealing on Legacy alone gives 4.624 Å (`docs/STATE_BRIEF_2026-09-12.md` §2).

## IV.2 AMBER: ff14SB with GBn2 implicit solvent, through OpenMM

**AMBER ff14SB** is a standard all-atom protein force field; **GBn2** is a generalised-Born implicit solvent model that approximates water as a continuum; **OpenMM** is the simulation library that evaluates them (`core/amber.py`). Because a CA trace has no side chains, the project builds all twenty residue types' side chains before evaluation (`sidechains` module; coverage went from 6 of 35 to 35 of 35 targets when that was done, memory `sidechains-all-20-residues`). A raw single-point AMBER energy on a built structure is dominated by steric clashes, so the deployed AMBER channel is a **converged interaction-only** energy: 50 steps of restrained relaxation, then the non-bonded plus solvation terms (`H_AMBER = E ∘ Relax50`; `s20/LEDGER.md`: the relaxed energy explains 58/97/99% of the variance at the three tested depths). Memory: AMBER runs are gated at 92% RAM (`core/amber.py`, `MEMORY_LIMIT_PERCENT=92.0`, line 249) and a chain that starts above `CONVERGE_MAX_KCAL=1000.0` kcal/mol (line 1282) is reported as unconverged. 63,000 AMBER single points are cached in `s24/cache_amber/` (tracked since S26, `s26/LEDGER.md` L1).

Two more facts about AMBER. On matched pools with a radius-of-gyration control it puts the native at the 51st percentile, i.e. it does not tell the native from its decoys (memory `physics-ranks-real-geometry-not-lattice`, overturned entry; `s20/LEDGER.md`). A polarisable force field (AMOEBA) does not fix this (memory `amoeba-does-not-fix-antiranking`). And its Pauli spectrum, once conditioned, is higher-weight than Legacy's (Part V.10).

## IV.3 Normalisation: rank, not moments

Every channel (distogram score, Legacy, AMBER) is put on the same scale by **rank standardisation**: `zrank(x) = (rank(x) − mean)/sd`, which with no ties has mean (K+1)/2 and sd √((K²−1)/12), so every channel lands on the identical marginal. A combined channel is `E_S = zrank(Σ_c zrank(x_c))`; the outer zrank is required because a |S|-channel sum has sd 1.00, 1.41 or 1.73 and CVaR trades energy against T·H, so the effective temperature would differ silently between configurations, worth up to 0.14 Å (`s25/LEDGER.md` L16, "Normalisation, stated mathematically"). Production normalises by rank; every moment z-score variant is research-only (`s26/LEDGER.md` L8). The mean |z_moment| per channel (AMBER 0.1127, distogram 0.7529, Legacy 0.8013) is quoted in L16 but the S26 examination found no artefact for the triple (Appendix D, D10).

## IV.4 The seven-configuration suite (S25 L16)

Seven configurations, every subset of {Legacy, AMBER, Distogram}, run under an identical CVaR-VQE (exact statevector, 9 qubits, 3 layers, 80 Adam iterations on the exact parameter-shift gradient, α = 0.18, T = 0.5, seed 0) with the same K = 500 pool, the same uniform coordinate average, on the **point-cloud** basis (`s25/results/phys_suite.json`: `alpha` 0.18, `T` 0.5, `layers` 3, `iters` 80, `basis` "point_cloud", `incumbent` 3.0483380938795324). α = 0.18 was pinned on a native-free criterion: a 12-target probe gave median realised tail 75, the production rung, so the VQE and the classical top-75 control consume the same number of candidates (`s25/LEDGER.md` L16). Results (`s25/results/phys_suite.json`, `configs_rank.<i>.mean`; all point-cloud basis):

| config | channels | mean Å | vs incumbent 3.0483 | verdict |
|---|---|---|---|---|
| C3 | Distogram | 3.0580 | +0.0097, 0.43× MDE | NOT MEASURED |
| C5 | AMBER + Distogram | 3.1317 | +0.0833, 0.70× MDE | NOT MEASURED |
| C4 | Legacy + Distogram | 3.2151 | +0.1668, 0.91× MDE | NOT MEASURED |
| C7 | Legacy + AMBER + Distogram | 3.2530 | +0.2047, 1.07× MDE | WORSE (type-M zone) |
| random 75 of 500 | none | 3.4251 | | control (`random_null_rank.mean`) |
| C6 | Legacy + AMBER | 3.6742 | +0.6259, 1.89× MDE | WORSE |
| C1 | Legacy | 3.7553 | +0.7070, 2.14× MDE | WORSE |
| C2 | AMBER | 3.8805 | +0.8322, 2.67× MDE | WORSE |

(`vs_incumbent.<i>` blocks give effect, mde, verdict.) The two headline contrasts are against the random-75 control: Legacy alone is +0.3302 Å worse than random (1.99× MDE, 41W/85L) and AMBER alone is +0.4554 Å worse than random (2.42× MDE, 30W/96L) (`controls_rank.1.vs_random`, `controls_rank.2.vs_random`). A rank-permuted physics channel, which carries no information, reproduces them: +0.3274 for Legacy and +0.4712 for AMBER, 26W/100L (`controls_rank.1.vs_perm`, `controls_rank.2.vs_perm`). Both physics energies are measurably worse than noise. Inside a distogram-led configuration AMBER adds +0.0331 at 0.43× MDE, nothing measurable (`controls_rank.5.vs_perm`). In every configuration the VQE matches the classical size-matched top-m at exactly 0.0000 (`controls_rank.<i>.vs_topm`, effect 0.0), which is the set-equality theorem (Part V.7) showing up in the suite.

An ORACLE that picks the better of C3 and C5 per target reaches 2.9400 Å (ORACLE, uses the native; `s25/results/phys_form5.json`, `oracle_ceiling.mean`), with k_eff 1.02: the two arms are 0.9536 correlated, so that oracle is an order statistic over nearly one arm (`oracle_ceiling.corr_C3_C5`, `k_eff`). Form 5 tried to reach it with a native-free switch under nested leave-one-fold-out: +0.0049 Å against C3 with 97 ties, and the permuted-signal null gives +0.0053 (`s25/results/phys_form5.json`, `primary.stats.effect`, `primary.null_permuted.effect`, `verdict` "FORM 5 CLOSED"; `s25/LEDGER.md` L18).

## IV.5 What "the physics does not rank" means and does not mean

It means: on pools of real fragment geometry for 9 to 16 residue peptides, the ordering by either energy is no better than a random ordering at picking the members closest to the native, and worse than the distogram's ordering. It does not mean the force fields are wrong about physics. Every candidate is a real fragment of a real protein, so all of them are physically plausible; the energies are being asked to discriminate between plausible shapes on a different sequence, which is not what they were built for. And AMBER's relaxation is still the only step that turns a built chain into an all-atom model with side chains (Part III, Stage 9), which is a validity role, not an accuracy role (`s26/LEDGER.md` L39).

---

# PART V. THE QUANTUM COMPONENT IN FULL

## V.1 The encoding

The problem posed to the circuit is: given 2ⁿ candidates with energies E_0 … E_{2ⁿ−1}, produce a probability distribution over them that concentrates on low energy. Each candidate is one computational basis state of an n-qubit register: candidate index i is the bit string of i. With n = 7 there are 128 basis states, so the shipped quantum path widens the top-75 cut to 128 (`core/pipeline.py`, `quantum_stage`, `dim = 2**cfg.vqe_qubits`; `Config.vqe_qubits=7`). The suite used n = 9 to index the full pool of 500 with 512 states (`s25/LEDGER.md` L16). The Hamiltonian is diagonal in the computational basis, `H = diag(E)`, and E is the rank-standardised score ladder (Part III, Stage 6). An earlier one-hot encoding (one qubit per candidate, 80 qubits) was infeasible, 0 of 20,000 sampled bit strings valid, and binary indexing (20 qubits for that instance) replaced it (`docs/FINDINGS.md` §6, "VQE over segment choices: encoding and honest limits").

Nothing in this encoding is "quantum" in the sense of a Hamiltonian with off-diagonal terms. The ground state of `diag(E)` is the basis state at the argmin. That is why the project's measured contribution of the circuit reduces to a statement about which classical readout it emulates (V.7).

## V.2 The ansatz

`StatevectorCircuit` in `core/quantum.py`: n qubits, `layers` = 3 repetitions of a layer of RY rotations (one parameter per qubit) followed by a CNOT chain (qubit i controls i+1) with a ring closure (last controls first), then a final layer of RY (`core/quantum.py`, `StatevectorCircuit`, `_entangler_permutation`, `states_batch`, `probs_batch`). Parameters: `layers × n` = 21 (`s25/results/q_verify.json`, "parameter count n_params = layers*n: 21"). Because the gates are RY and CNOT, the amplitudes are real, and the circuit lives in the real orthogonal group; its dynamical Lie algebra is inside so(2ⁿ), dimension 8128 for n = 7 (Part V.10). The circuit is simulated exactly as a dense statevector; a matrix product state (**MPS**, a tensor-network form whose cost is O(n χ³) for bond dimension χ) reproduces it to 3.331e-16 with χ = 2^layers exactly, χ = 2, 4, 8, 16 at 1 to 4 layers, Schmidt rank (the number of non-zero singular values when the state is split between the two halves of the register) across the middle cut 4 at the 2-layer check (`s25/results/q_verify.json`, "max |p_MPS - p_dense| over 24 configs", "chi by layers"; `core/quantum.py`, `MPSAnsatz`, line 412). No sampling occurs anywhere in the trained path: "sampling calls inside run_cvar_vqe/free_energy: NONE" (`s25/results/q_verify.json`); the random number generator is used for the initial angles only, θ ~ N(0, 0.6²) (`core/quantum.py`, `run_cvar_vqe`).

## V.3 The objective: CVaR, then free energy

For a distribution p over energies E, sorted ascending, the **CVaR at level α** is the mean energy of the lowest-α fraction of the probability mass:

    CVaR_α(p) = (1/α) Σ_{i in tail} p_i E_i,   where the tail is the smallest set of lowest-E states with mass ≥ α, the boundary state counted fractionally.

α = 1 is the ordinary expectation; α → 0 is the minimum. CVaR-VQE (Barkoutsos et al. 2020) uses this instead of the expectation so that the optimiser cares about the best shots rather than the average. The exact gradient uses the envelope theorem, ∂CVaR/∂p_i = (E_i − q)/α for tail states, q the boundary energy (`core/quantum.py`, `cvar_exact`).

The deployed objective is a **free energy**, F(p) = CVaR_α(p) − T · H(p), where H is the Shannon entropy of p in nats and T a temperature (`core/quantum.py`, `free_energy`). The entropy term stops the distribution collapsing to a single state. The deployed (α, T) pairs are per fold: `VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}` (`core/pipeline.py`, `VQE_LFO`), i.e. T = 0.3 on every fold and α = 1.0 on three of five, which is a plain entropy-regularised expectation, not a CVaR at all (`s25/LEDGER.md` L3). The table was chosen leave-fold-out on the S8 instrument; the audit reproduces it exactly, per-fold match True on all five, mean 3.3135 (`verify/vqe_lfo_audit.json`, `LFO_TABLE_REPRODUCES` True, `recomputed_LFO_mean` 3.3135).

## V.4 The optimiser and the gradient

Adam, 50 iterations, learning rate 0.15, β₁ = 0.9, β₂ = 0.999, one restart, seed 0 (`core/quantum.py`, `run_cvar_vqe`; `Config.vqe_iters=50, vqe_seed=0`). The gradient of the objective with respect to the 21 angles is the exact **parameter-shift** gradient (for an RY gate, ∂/∂θ ⟨f⟩ = ½[f(θ+π/2) − f(θ−π/2)]), and is verified against an independent finite-difference gradient: cosine 1.000000000 (minimum over 12 checks), relative error 4.597e-10 (`s25/results/q_verify.json`, "cos(param-shift, exact FD)", "rel |g_ps - g_fd| / |g_fd|"; also `verify/cvar_audit.json`, `B_paramshift_vs_INDEPENDENT_fd_cosine` 1.0). The free-energy gradient is verified the same way: cosine 1.000000000, relative error 4.663e-10 over four (α, T) cells (`q_verify.json`). The dense statevector and the reference implementation agree to 5.551e-17 (`q_verify.json`, "max |p_statevector - p_dense|").

**The gradient defect, and its fix.** An earlier score-function (REINFORCE-style) gradient centred its baseline on the tail only, so the correction term b · E_p[∇ log p] was not zero. Its cosine with the true gradient was −0.023 when first priced (`docs/FINDINGS.md`, corrections table, "The CVaR gradient as shipped"; `docs/FINDINGS.md` §6). Measured against the exact gradient on the deployed instrument, the tail-baseline estimator has cosine 0.566586 and 0.519 of the true norm (`s25/results/q_verify.json`, "cos(exact-expectation, TAIL baseline)", "|g_tail| / |g_exact|"); on the audit instrument 0.6847 (`verify/cvar_audit.json`, `D_tail_baseline_cosine_vs_paramshift`). The number depends on the instrument; the sampled tail-only cosine 0.524 quoted in some documents has no artefact (Appendix D, D10). The constant baseline has cosine 1.000000 and is the default everywhere (`q_verify.json`, "cvar_gradient default baseline: const", "grad_cvar_score default baseline: const"; `verify/cvar_audit.json`, `E_all_defaults_const` True).

## V.5 Readout

The trained distribution p over 128 candidates is consumed in one of two ways: `consensus_medoid(block, p)` picks the p-weighted medoid; `average_weighted` forms the p-weighted coordinate average (`core/pipeline.py`). Both have a uniform ablation (`sel_uniform`). The readout is insensitive to distributional differences of nearly half the mass: the trained state at T = 0.3 differs from its Gibbs target by total variation 0.453, and the endpoint difference between the two is −0.0302 Å, 0.24× MDE, NULL (`s25/results/q_gibbs.json`, `divergence_ladder`, T=0.3 row, `tv` 0.4531; `s25/results/q_alpha.json`, `expressivity_vs_boltzmann`, "circuit - exact Boltzmann @T=0.3", 3.2835 vs 3.3138, MDE 0.1261; `s25/LEDGER.md` L15).

## V.6 What the circuit is worth, measured

On the S8 instrument (consensus-medoid selection from a 128-candidate score-filtered set; never comparable to the 3.0483 point-cloud figure), the per-fold VQE_LFO table gives 3.3135 Å against the plain argmin's 3.4540 Å: −0.1405, SE 0.0732, MDE 0.2051, 0.68× MDE, verdict UNDERPOWERED, NOT A RESULT (`s25/results/q_alpha.json`, `vs_no_circuit`, "VQE_LFO - argmin (shipped)"; `verify/vqe_lfo_audit.json`, `shipped_LFO_vs_argmin` −0.1405). Against an exact Boltzmann distribution at T = 0.3 (no circuit), the VQE is −0.0002: NULL (`q_alpha.json`, "VQE_LFO - Boltzmann T=0.3", 3.3135 vs 3.3138). Against uniform top-64, +0.026, NULL. The best cell on the grid, (α = 0.25, T = 0.3), is −0.1715 against argmin (`verify/vqe_lfo_audit.json`, `best_global_vs_argmin`), also below its MDE.

**The "+0.113 Å CVaR contribution" was withdrawn in S25 (L3, L5).** It was the difference between α = 0.1 and α = 1.0 at T = 0.1 (3.3414 vs 3.4540, entropy 6.36 bits vs 0.076 bits; `q_alpha.json`, `inherited_claim`). Two things were wrong with it. T = 0.1 does not ship; at the deployed T = 0.3 the same contrast is +0.0279 with the sign reversed, and the α = 1 arm does not collapse (5.66 bits) (`q_alpha.json`, `primary_alpha_effect_by_T`, "a=0.1 - a=1.0 @T=0.3", 3.3115 vs 3.2835; `s25/LEDGER.md` L3). And even at T = 0.1 it was a marginal mean where a paired statistic was required: paired, it is −0.1126 at 0.51× MDE with median exactly 0.0000, a NULL (`s25/LEDGER.md` L5). The `quantum_stage` docstring in `core/pipeline.py` still carries the withdrawn claim (Appendix D, D5).

**One curve explains everything.** Over the 18 readout arms (9 VQE cells, 3 Boltzmann, 4 uniform top-fractions, argmin, medoid), mean RMSD is a function of readout entropy alone: correlation −0.7423 with H, +0.2700 with α, −0.0234 with T; a quadratic in H fitted on the 9 no-circuit arms has R² 0.704, optimum H* 7.44 bits (above the 7-bit maximum, so "more entropy is better" all the way to uniform over 128); the VQE arms sit on the curve with mean residual +0.0090 against a no-circuit residual sd of 0.0268; α's marginal share of variance after H and H² is 0.032 (`s25/results/q_alpha.json`, `cell_correlations`, `entropy_fit`, `circuit_off_curve`, `alpha_beyond_entropy`). The circuit is a way of choosing an entropy. The entropy is what matters. And 61.9% of targets have no tail constraint at all at the deployed α (`q_alpha.json`, `share_of_targets_with_no_tail_constraint` 0.619).

## V.7 The set-equality theorem, and why the Hamiltonian barely changes

**Theorem (S24, proved; S25 verified on 2,592 cells).** For a diagonal Hamiltonian and any state, the CVaR tail at level α is a subset of the prefix of the energy ordering of length equal to the realised tail. Verified: 2,592 cells, 1,620 containing an exact zero probability (so the tail can skip states), 0 subset-hood violations, 0 holes that were not exactly zero probability, and on the 972 full-support cells the tail equals the prefix exactly, 972 of 972 (`s25/results/q_verify.json`, "SUBSET-hood violations: 0", "full-support cells with EXACT value equality: 972 / 972"). The stronger S24 claim, that the tail is always a prefix, was an over-claim, violated on 1,424 cells (`q_verify.json`, "prefix-hood violations (the s24 REV1 over-claim): 1424"). Consequence: whatever the circuit does, the candidates it can put in the tail are the top of the classical ranking. That is why the `selector` column in the results table says "score-filter top-75 (== CVaR tail, s24 theorem)" and why every suite configuration matches its size-matched classical top-m at 0.0000 (Part IV.4).

**S25 L17: two trained states.** Because E is the rank-standardised score, the energy vector is the same ladder for every target to within 1.18% of its range (worst deviation 0.0406 on a range of 3.4371 over eight targets checked; 0 exactly identical; `s25/results/q_gibbs.json`, `spectrum_target_independence`). The circuit therefore trains toward one of two states, one per (α, T) cell in `VQE_LFO`, regardless of the target; the target enters only through which candidates the indices point at. This is the reason five sprints of "change the Hamiltonian" experiments returned nulls: the Hamiltonian's shape was never changing (`s25/LEDGER.md` L17).

## V.8 Gibbs identity and how well the circuit trains

For F(p) = ⟨E⟩_p − T·H(p) (the α = 1 case), F(p) − F(p*) = T · KL(p ‖ p*), where p* is the Gibbs distribution at T. So the free-energy gap is a divergence, and the ladder of states can be read as a divergence ladder (`s25/results/q_gibbs.json`, `divergence_ladder`; n = 7, 3 layers, 50 iterations, seed 0, `identity_residual` 2.8e-16):

| T | state | F | KL nats | KL bits | TV | H bits |
|---|---|---|---|---|---|---|
| 0.3 (deployed) | random θ | −1.2754 | 3.9274 | 5.666 | 0.7805 | 4.456 |
| 0.3 | 50 Adam steps (deployed) | −2.1831 | 0.9019 | 1.301 | 0.4531 | 5.6706 |
| 0.3 | uniform over 128 | −1.4556 | 3.3269 | 4.800 | 0.7015 | 7.0 |
| 0.3 | point mass at argmin | −1.7186 | 2.4503 | 3.535 | 0.9137 | 0.0 |
| 0.3 | Gibbs optimum | −2.4537 | 0 | 0 | 0 | 4.9135 |
| 0.1 | 50 Adam steps | −1.7176 | 1.4492 | 2.091 | 0.7607 | 0.0747 |
| 0.1 | point mass at argmin | −1.7186 | 1.4392 | 2.076 | 0.7629 | 0.0 |
| 0.1 | Gibbs optimum | −1.8625 | 0 | 0 | 0 | 3.3326 |

(`q_gibbs.json`, `divergence_ladder`, rows as labelled.) Two readings. At the deployed T = 0.3 the trained state closes 78.3% of the gap from the initial mean to the Gibbs optimum, beats the best of 200 random initialisations by a margin of 0.70, and is still 0.90 nats (1.30 bits) from the Gibbs state, so it trains, and it does not reach the target (`q_gibbs.json`, `training_control`, T = 0.3: `frac_gap_closed` 0.783, `beats_best_of_N` True, `margin` 0.6996, N 200). At T = 0.1 the trained state, F = −1.7176, is worse than the point mass at the argmin, F = −1.7186: the entropy term is too small to matter and the circuit is a slow way to find a minimum it could read off (`q_gibbs.json`, T = 0.1 rows). The pre-registered falsifier (KL below 0.1 nats) fired against the draft claim that the circuit reaches Gibbs (`q_gibbs.json`, `falsifier.fired_against_the_draft` True).

## V.9 Gradient variance: the width and depth sweeps, and what may be called a barren plateau

A **barren plateau** is the regime where the variance of the cost gradient over random parameters decays exponentially in the number of qubits, so that a random initialisation sees a flat surface. The project measured Var[∂F/∂θ₀] over 250 random θ (200 at n = 10, 120 at n = 12, 80 at n = 13) at depth 3 for n = 4 to 13 and fitted log₂ Var per qubit (`s25/results/q_plateau.json`, `*.rows`, `log2_slope_per_qubit`; `n_theta` per row):

| cell | slope (log₂ Var per qubit) |
|---|---|
| linear, α = 1, T = 0 | −0.6492 |
| CVaR α = 0.25, T = 0 | −0.2522 |
| CVaR α = 0.10, T = 0 | −0.0472 |
| deployed α = 1, T = 0.3 | −0.3105 |
| deployed α = 0.25, T = 0.3 | −0.2429 |

A 2-design (a circuit whose random outputs match the Haar measure to second moments) on real states gives a slope of −1 per qubit (Var ~ 1/dim). The measured slopes at depth 3 are shallower than that, and the CVaR objectives shallower still. The depth sweep at n = 7, α = 0.25, T = 0.3 gives Var 2.381e-2 at 1 layer, 1.66e-2 at 2, 9.197e-3 at 3, 7.999e-3 at 4, 7.6e-3 at 6, 7.4e-3 at 8, 8.213e-3 at 12: it saturates by 4 layers at n = 7 (`q_plateau.json`, `depth_sweep_n7`, `var_g0`). The ratio of CVaR-gradient variance to linear-gradient variance grows with n: 0.5159 at n = 7, 3.3585 at n = 13 (first entry of `nonlinearity_variance_ratio` per width), which is the tail's nonlinearity in p amplifying the gradient as the register grows.

**What may and may not be said.** The correct sentence is: "at depth 3 and n = 4 to 13, the deployed objective's gradient variance decays at −0.24 to −0.31 log₂ per qubit, well short of the 2-design rate, and the circuit is simulable because n = 7, not because of its structure." The S26 DLA measurement (V.10) is the reason for the "at depth 3" clause. The forbidden sentences are "the CVaR objective avoids barren plateaus" and "a small gradient here is a barren plateau" (`s26/LEDGER.md` L27, last paragraph; `s26/briefs/PR.md` rule "Never call a small gradient a barren plateau"). The parameter count, 21, is tiny against dim so(2⁷) = 8128, so no 2-design claim was ever possible at the deployed cell (`s25/agentQ_FINDINGS.md`).

## V.10 The Pauli spectrum, the locality theorem, and the DLA

**Locality theorem (S13).** In torsion space, the CA to CA distance d_ij depends on exactly the j − i − 1 residues between i and j; all-atom supports are one residue wider. Verified on real geometry with agreement 1.0000 (`s13/results/qarch_locality_geom.json`, `qarch_locality_amber.json`; memory `torsion-space-locality-theorem`). Both energies are full-register, so "AMBER is less local" is a category error at this length.

**Pauli spectrum predicts gradient variance (S13).** Writing a diagonal energy over 2^m torsion states in the Walsh (Pauli-Z) basis, the mean Pauli weight of Legacy is 2.236 and of conditioned AMBER 3.015, AMBER higher on 79 of 79 cells, with 1-plus-2-body share 0.778 for Legacy and 0.629 for AMBER; neither is 2-local (`s13/results/geo_pauli.json`; `s13/SPRINT13_DOSSIER.md` §8). Raw AMBER's spectrum is a delta spike: the top-10 configurations carry 99.6% of its Walsh variance and its mean weight lands on Binomial(m, ½) exactly, 6.001 for m = 12 (`s13/results/walsh_xval.json`; memory `pauli-spectrum-delta-spike-artefact`). The spectrum times the circuit's gradient kernel predicts the measured gradient variance to 0.4% (`s26/PROPOSAL_B_REPLACEMENT.md` §1, citing `s13/results/walsh_predict.json`). Quantum natural gradient was a null here: the metric is g_ii = 0.25 on the diagonal for this ansatz (`s13/SPRINT13_DOSSIER.md`).

**DLA (S26, measured after the report scope but bearing on every sentence above).** The dynamical Lie algebra of the RY/CNOT chain+ring ansatz is the full so(2ⁿ) from depth 2 at n = 4, 5, 7, 8, 10, 11 (8128 at the deployed n = 7), abelian of dimension n at depth 1, and a proper subalgebra only at n = 6 and n = 9 (`s26/results/q_dla.json`; `s26/LEDGER.md` L27). The pre-registered prediction of a proper subalgebra at (7, 3) was falsified. Consequence: nothing algebraic protects this ansatz from a plateau at scale; the S25 result is a depth-3 statement.

## V.11 Entanglement and the continuous encoding (S20)

On the continuous (torsion-window) encoding the VQE does train: −0.210 Å [−0.339, −0.082] against its own untrained circuit, 5 of 5 folds, and −0.299 [−0.407, −0.195] when used for generation (`s20/LEDGER.md` D4; `s20/agentD_FINDINGS.md`). This scoped the earlier lattice-only finding that running the VQE was worse than not running it (memory `concentration-is-wrong-when-discrimination-binds`, S20 scope fix). Entanglement does no work: an entangler ablation is −0.013 [−0.095, +0.077] and the mutual information across the middle cut is 0.045 bits (`s20/LEDGER.md`). The S26 ADAPT experiments (V.12) make the same point from the other side.

## V.12 ADAPT-VQE on this Hamiltonian (S26, in progress at the time of writing)

Because the presenter's audience works on ADAPT-VQE, this is stated carefully. S26 lane Q pre-registered four property measurements (`s26/PREREG_A1.md` to `PREREG_A4.md`). Two are complete. A2 is the DLA result above. A4 grew circuits operator-by-operator on the deployed rank-ladder Hamiltonian from two pools (Tang's minimal 2n−2 pool "V", and all 1- and 2-local odd-Y strings "L2"), by Adam best-iterate, and measured gradient variance versus n (`s26/results/q_var.json`; `s26/LEDGER.md` L35). At α = 1 the grown circuits hold 1 to 3 distinct operators, the rest being repeats of the same single-qubit Y: they are product circuits, their variance does not decay with n (+0.035 and −0.008 log₂ per qubit), and "a large gradient from a product circuit is not trainability; it is the trivial regime". At α = 0.25, T = 0.3, pool L2 grows 4 to 21 distinct 2-local strings and decays at −0.302 per qubit, the fixed ansatz's −0.243 within the error of a 7-point slope. At α = 1 the state the optimiser heads for is a basis state, a product state that seven parameters represent exactly, so an adaptive ansatz grown on it has nothing to train (`s26/LEDGER.md` L35). A1, the RMSD endpoint of ADAPT against the fixed ansatz, is pending its phase gate, expected null, with 30 per-target files on disk and not scored (`s26/agentQ_FINDINGS.md` §1; `s26/results/a1/`). The Tang pools generate so(2^(n−1)+1), about a quarter of so(2ⁿ) (`q_dla.json`, `pools`; L27).

## V.13 What the quantum component is, in one paragraph for a VQA audience

It is an exact, sampling-free, 7-qubit, 21-parameter RY/CNOT ansatz trained by Adam on a free energy F = CVaR_α − T·H over a diagonal rank-ladder Hamiltonian, with exact parameter-shift gradients verified to 1e-10 against finite differences, an exact MPS twin, and a proved-and-verified theorem that its tail is a subset of the classical prefix. It is deployed off. It contributes nothing measurable to accuracy, and it was measured well enough to say so at 0.68× MDE. What it produced that is worth a VQA group's time is the chain "chain geometry → locality → Pauli spectrum → gradient kernel → measured variance to 0.4%", the depth-3 width sweep with its correct scope, the DLA census, and the ADAPT-grown product-circuit observation.


---

# PART VI. THE TWENTY-FIVE SPRINTS AS A STORY

Sprints 1 to 4 built the first lattice model, the first CVaR-VQE, and the first Amber scoring, and are recorded only in memory notes and in the corrections that later sprints made to them (for example, Sprint 4's physics measurement was confounded, `docs/FINDINGS.md` S5-12). The written record begins at Sprint 5. Pre-registration files begin at Sprint 16 (`s16/PREREG_steer.md` is the first); before that, the falsifier line below says "none registered" and quotes the expectation the sprint stated in its own dossier where one exists.

Each sprint: the question, the pre-registered falsifier, the result, what it closed, what it retracted, what it left open.

## Sprint 5: the wall is located

**Question.** Why does the pipeline stop at about 3.3 Å, and which stage is responsible?
**Falsifier.** None registered.
**Result.** The pool was never the problem and search is not the constraint at any budget; the wall is "no resolution at the top": inside the near-native band nothing available separates the candidates (`docs/FINDINGS.md` S5-1 to S5-4). In-band skill is bounded by the prior's own accuracy (S5-7). Structure and sequence are decoupled at this length: the best-matching window has about 12% identity, so retrieval by sequence cannot find the native's relatives (S5-9). The torsion prior is the sixth signal to hit the same wall (S5-10). Physics does rank real geometry once the Sprint 4 confound is removed (S5-12; overturned again in S20 with a radius-of-gyration control, Part IV.2).
**Closed.** Search budget as a lever. Retrieval by predicted structure (S5-13: better selection, worse pool, no net gain).
**Retracted.** The shipped CVaR gradient was found defective (S5-6; priced in S9-5 at cosine −0.023) and fixed (Part V.4).
**Open.** Everything about the prior.

## Sprint 6: selection is measured out

**Question.** Can a learned ranker, a built structure, or an all-atom combiner pick better inside the pool?
**Falsifier.** None registered.
**Result.** Learned rankers fix the global ranking and not the in-band one; properly powered they are exactly null (S6-2, S6-4). Building a structure instead of picking one: null (S6-5). Tail aggregation of pair violations: null (S6-7). VQE/CVaR assembly gives a narrow pool whose narrowness cannot be exploited (S6-8). Decomposed all-atom physics in a learned combiner: null (S6-12). The distance channel is sufficient; the predictor is not (S6-13) (`docs/FINDINGS.md`, Sprint 6 index).
**Closed.** Learned selection from the existing features.
**Retracted.** A correctness bug in the pairwise ranker arm was found before the arm was reported (S6-3; corrections table).
**Open.** The predictor.

## Sprint 7: the predictor is the constraint

**Question.** Is the distance predictor biased, under-trained, or misspecified?
**Falsifier.** None registered; the coordinator named the calibration slope in advance as the discriminator for the shrinkage hypothesis (S7-3).
**Result.** Thirty times more structural training data makes the predictor monotonically worse (S7-2). It is not shrunk toward the average peptide; its error is correlated across pairs and only 28% aligned with truth, calibration slope +0.376 (S7-3). The native is not the objective's minimum (S7-5: native is the argmin on 3 of 126 targets, 36.8th percentile). Pool size has an interior optimum and selection tracks the pool mean, not its best (S7-6). No scoring channel ranks the truth first; AMBER is worse, not better (S7-10). ESM does beat one-hot on selection, −0.288 Å [−0.484, −0.092], p = 0.0046 (S7-11; the per-target artefact is lost, Appendix D, D9). BLOSUM retrieval beats random fragments at every K (S7-12) (`docs/FINDINGS.md`, Sprint 7).
**Closed.** Post-hoc de-biasing; more training data; chirality on the generative path (S7-4).
**Retracted.** The shrinkage diagnosis (refuted by the experiment it commissioned). The memory "ESM adds nothing" (judged by MAE, the wrong metric). The memory "real fragments beat the lattice" read as "random beats BLOSUM" (refuted on a matched test) (corrections table).
**Open.** Whether any architecture consumes the predictor better.

## Sprint 8: full architectural attack

**Question.** Attack every stage at once: representation, aggregation, pool construction, generation, synthesis, AMBER.
**Falsifier.** None registered.
**Result.** The representation is not the bottleneck; the aggregation metric is (S8-1). Multimodality is modest, 2 to 3 populated basins (S8-2). Pool construction is capped at 2.406 Å at any filter skill (S8-5; an ORACLE-type bound, it uses the native to define skill). Inside the near-native band the score is worse than a coin flip and consensus is the only discriminator (S8-8); this is where the tie-breaking trap was found: `np.argmin` on a tied signal read the ORACLE sort order and invented a 1.386 Å winner (S8-8; memory `consensus-is-the-only-in-band-discriminator`). Synthesis (averaging then projecting) beats selection: 3.204 Å, replicated on dev (S8-11). AMBER is a validity stage, not an accuracy stage (S8-12). The training target is misspecified by 0.5 Å per pair and fixing it is worth 0.044 Å (S8-14; did not replicate on dev, +0.070 [−0.118, +0.259], corrections table) (`docs/FINDINGS.md`, Sprint 8).
**Closed.** Generation as a source of candidates (S8-10: worth nothing as candidates).
**Retracted.** S8-4's class decomposition read as a "28% recoverable opportunity" (refuted by S8-6: injecting the library's best window by cheating is worth −0.016 Å [−0.037, −0.001]; the opportunity is 0.004 Å). S8-10's "transfer law" `selected = 0.915 × ref + 0.346` (withdrawn by S8-13: in the operating band the map is the identity, `1.009 × ref + 0.011`, r = 0.981, and passing the best available structure through selection costs +0.049 [−0.009, +0.103]) (corrections table).
**Open.** Whether the synthesis gain transfers to the sealed benchmark.

## Sprint 9: breaking the 2.0 Å barrier, and the benchmark

**Question.** Can refinement, a torsion prior, evolutionary information, or a loop push below 2 Å? Then: spend the benchmark once.
**Falsifier.** None registered; the benchmark pre-registration recorded the configuration and the tuning mean 3.204 (`verify/headline_audit.json`, `benchmark60.preregistration`).
**Result.** The shared bias is real, one interpretable mode, and is not an error (S9-1). A torsion prior buys physical validity and no accuracy (S9-2). Refinement fails and the failure is all selection (S9-5). Evolutionary information exists, is retrievable, and does not help (S9-7). The in-band ceiling is informational, not architectural (S9-8). **The benchmark pass**: 2.9610 vs 2.9507 baseline, +0.0103 [−0.1596, +0.1803], 31W/29L, SE 0.0867; the tuning gain −0.250 [−0.383, −0.117] and dev gain −0.255 did not transfer; dropping the ten largest tuning gains leaves the system +0.2013 Å worse than baseline; the identity leak on the benchmark prices at +0.0030 (S9-10; `s9/final_report.json` via `verify/headline_audit.json`). Concentration, reported as a caveat for four sprints, was the result (`docs/FINDINGS.md`, "Two structural lessons").
**Closed.** The loop (S9-4: a contraction toward its own fixed point near 3.2 Å, never improving on its input). Chemical-shift and evolutionary side channels.
**Retracted.** S9-4's reading "every 1 Å of query improvement buys 0.49 Å" (reframed). S9-5's positive-phi defect attributed to consensus (relocated to the projection by S9-9) (corrections table).
**Open.** The recognition line: why non-replication happened.

## Sprint 10: adjudication and bounds

**Question.** Was the non-replication contamination, and where is the geometric bound?
**Falsifier.** None registered.
**Result.** The loss is set purity and the failure is a minority class (S10-1). The projection path is closed: the required matrix is 2.7× better than anything achievable (S10-2). Physics does not supply decorrelated error (S10-3). The identity-audit conflict was an alphabet permutation between two encoders (`ARNDCQEGHILKMFPSTWYV` vs `ACDEFGHIKLMNPQRSTVWY`, agreeing on three letters): the true figure is 13 of 126 targets with a pool window at ≥0.6 identity, four at 1.0, and the leak prices at +0.0004 Å (S10-4). The reconciled bound ladder: geometry is not the barrier; headroom inside the top-75 is 1.36 Å (0.90 by selection alone, 1.16 by re-weighting, the last 0.20 only with ORACLE per-candidate poses) and inside the 500 is 2.35 Å (1.49 by selection, 2.11 by re-weighting); the top-75 convex hull emits at 1.840 Å and the K=500 hull at 0.853 Å; all of these are ORACLE diagnostics that use the native and none is deployable (S10-5; `docs/FINDINGS.md`, bottom line, "Where the remaining accuracy actually is"; these are the slide-7 numbers) (`docs/FINDINGS.md`, Sprint 10).
**Closed.** Contamination as the explanation. The projection-matrix route.
**Retracted.** S10-1's identity audit, void in full. S10-2's "10/126" revised to 13/126. `s10/forensics` ceilings 2.087/2.160/1.893 revised to 1.987/2.044/1.802 (a solver defect). `s10/mdgen` figures rescoped from 32 targets to 126 (1.802/1.336/0.953). The premise that oracle weighting cannot be worse than the hull (refuted on 9 of 378 cells) (corrections table).
**Open.** Recognition inside the band.

## Sprint 11: performance engineering and consolidation

**Question.** Move everything into `core/`, make it fast, and prove nothing changed.
**Falsifier.** Bit-exactness against the reference implementation.
**Result.** Eight workers reproduce one worker to ±4.4e-16 on every arm, 8.371× faster (`bench_results/compare_tuning126.json`). The four mandated components (VQE weights, Legacy refinement, AMBER, projection) run and three buy nothing measurable: VQE weights vs uniform −0.031 Å at 49W/49L; Legacy refinement +0.030; AMBER +0.021 cost for validity; the whole system +0.017 [−0.019, +0.053] against the like-for-like pipeline (S11-1). The projection is not a function of its input at Ångström resolution: under a rigid motion of the input cloud the reference disagrees with itself by up to 1.62 Å, so bit-exactness is the only defensible standard (S11-2) (`docs/FINDINGS.md`, Sprint 11).
**Closed.** Speed as a constraint.
**Retracted.** Nothing.
**Open.** Everything scientific.

## Sprint 12: adversarial, aggregation, assembly

**Question.** Can a set-transformer learn in-band discrimination? Is multi-piece assembly a route? What does the terminal operator consume?
**Falsifier.** None registered; the dossier states expectations per experiment (`s12/SPRINT12_DOSSIER.md`).
**Result.** The set-transformer (a neural network that consumes an unordered set of candidates) over the full signed deviation map has a flat learning curve, 3.043 → 3.026 across an 8× change in training data, while a leaked label is loud at n = 8 (2.534): in-band discrimination is signal-limited, not sample-limited (`s12/results/agg_dec_v2*.json`; reproduced to the third decimal in S26, `s26/PROPOSAL_C.md` C1). The operator law: d_out = 1.16 × mean + 0.04 × best, R² 0.893 on narrow sets (0.803/0.298 on wide sets); the terminal operator consumes the set mean, not the set best, so a perfect rank-1 is worth −1.74 Å through argmin and −0.03 Å through the m = 75 average. The blind pipeline (no sequence) gives 3.989 vs 3.213 shipped, but on FAIL18 blind 5.425 beats shipped 6.019: sequence conditioning is harmful on the failures. All 204 clusters of 9 to 16-mers are spent; the containment-fresh world supply is 16 targets, 10 of them amyloid fibrils. A QUBO formulation (quadratic unconstrained binary optimisation, the standard form for annealers) over 1,008 instances: simulated annealing finds the exact optimum 90% of the time, greedy 40.5%, and the exact optimum is +0.169 Å worse than the shipped answer (`s12/SPRINT12_DOSSIER.md`).
**Closed.** Learned in-band aggregation. Multi-piece assembly. A fresh benchmark (none exists). QUBO search as a lever.
**Retracted.** Nothing.
**Open.** The prior.

## Sprint 13: torsion space and quantum architecture

**Question.** Can a torsion-space encoding reach 2 Å, and what do the physics energies look like as Hamiltonians?
**Falsifier.** None registered; predictions stated in the dossier.
**Result.** The ORACLE torsion ceiling is 1.594 Å at about 24 live qubits (ORACLE, uses the native), but 88% of it is generic Ramachandran and 0.388 Å is the oracle start (random start 1.982). A constant alpha helix, 4.065 Å, beats a random torus draw, 4.522, by −0.457, so "beats random" proves nothing in torsion space. Legacy's low-decile correlation with RMSD is +0.043 and its certified optimum is worse than random. Sequence-only torsion prediction gives 3.770 Å, +0.558 against shipped; φ carries no sequence signal (36.1° vs 36.4°) and the whole channel is 10.4° of ψ. The locality theorem, the Pauli spectrum (Legacy 2.236, AMBER 3.015, 79/79 cells), the delta-spike artefact of raw AMBER (99.6% of Walsh variance in the top 10 configurations, mean weight Binomial(m, ½) = 6.001), the Pauli-spectrum prediction of gradient variance, and the QNG null (g_ii = 0.25) all land here (`s13/SPRINT13_DOSSIER.md`; `s13/results/geo_pauli.json`, `walsh_xval.json`, `walsh_predict.json`, `qarch_locality_*.json`).
**Closed.** Torsion space as the route to 2 Å. Chemical shifts (later, S14). QNG.
**Retracted.** Nothing.
**Open.** Publishing the trainability chain.

## Sprint 14: chemical-shift restraints and the VQE redesign

**Question.** Can NMR chemical shifts (per-residue experimental numbers that constrain torsions) supply the missing per-target information?
**Falsifier.** None registered.
**Result.** Shifts exist for only 54 of 126 targets; ORACLE-perfect torsions on all of them still leave the instrument at 2.021 Å, and 55 targets would need to be fixed (ORACLE, uses the native). Collapsing a bimodal posterior to its mode costs −2.253 Å against not collapsing it. A certified structural optimum is 0.885 Å better than random where Legacy's is 0.139 Å worse. The VQE is the worst optimiser tested, 0.00 to 0.12 against greedy's 0.71 on the certification rate (`s14/LEDGER.md`).
**Closed.** Chemical shifts as a route to 2.0 Å (by arithmetic).
**Retracted.** The C7 rank correlation +0.355, measured under a proposal drawn from the prior, corrected to +0.161 under a uniform proposal (`s14/LEDGER.md`, C7-CORRECTION; this is the closest artefact to the "+0.370" the prompt lists, Appendix D, D6). The 1.486 Å figure and the coverage gate, superseded.
**Open.** "Do not collapse the posterior" as a design rule.

## Sprint 15: the paper-driven programme

**Question.** Implement the published-literature pipeline the sprint was handed and see whether it reaches its claimed number.
**Falsifier.** None registered.
**Result.** The full-instrument cascade emits 3.321 Å, a significant loss of +0.117 [+0.050, +0.188] against the incumbent (`s15/FINAL_DOSSIER.md`). Torsion distance geometry from the predicted distances gives 3.644 Å, +0.440 [+0.290, +0.592] (`s15/LEDGER.md`, row 1.3, `distgeo.json`). Error shape beats error magnitude: MAE is unusable as a ranking metric, and the prior points the wrong way (calibration slope +0.376, correlation +0.283: full amplitude, wrong direction). The true ORACLE distance-geometry floor from the native's own distances is 0.611 Å (ORACLE, uses the native) (`s15/LEDGER.md`).
**Closed.** The paper's route. Post-hoc calibration.
**Retracted.** Eleven interim claims of the sprint, all listed in `s15/LEDGER.md`.
**Open.** The prior's direction.


## Sprint 16: native-free steering

**Question.** Can a native-free criterion steer the average (choose a subset, a weighting, or a repair) and beat "do nothing"?
**Falsifier.** `s16/PREREG_steer.md`: a steering arm counts only if it beats the size-matched random repair past its own MDE; a "quiet control" (a random displacement of the same magnitude as the proposed repair) is required, and if the gain matches ‖J·δ‖/√n or vanishes as the step shrinks, it is an operator artefact.
**Result.** Every steering arm chose "do nothing" on every fold. The fusion law for combining channels is the Krogh and Vedelsby ambiguity decomposition: fusion gain goes as the square of the weaker channel's skill, so decorrelated errors that exist (truth-partialled 0.04 to 0.26) buy +0.004 to +0.011 Å (`s16/LEDGER.md`; memory `decorrelated-errors-exist-but-are-unusable`). Repair: the projection costs +0.1554, AMBER at k = 30 +0.1333, and a random displacement matched in magnitude gives 3.1606, cosine to the native direction −0.0521: AMBER relaxation's apparent gain dissolves under the matched-random control (`s16/LEDGER.md` L27). The VQE loses to greedy 10 of 10 on certification (`s16/LEDGER.md`).
**Closed.** Native-free steering of the average. AMBER relaxation as an accuracy step (confirmed by S26 L39).
**Retracted.** The `csteer` |c| arm (a sign convention error, `s16/LEDGER.md`).
**Open.** Per-target set size.

## Sprint 17: from averaging to selection readout

**Question.** Does a shortlist or a ranker bind? Does widening K help? Can any of 58 distogram functionals or the VQE beat the shipped cut?
**Falsifier.** `s17/PREREG_select.md`: "If `sel-dist − sel-random` has a CI containing zero at n = 126, the distance objective does not select"; "If the LFO-selected variant's advantage over shipped has a CI containing zero" it is null. `s17/PREREG_quantum.md`: if no VQE configuration reaches a frontier the classical greedy reaches at equal budget, or exhaustive enumeration of the neighbourhood costs no more than the VQE's budget, the VQE is closed; if the mean-field control matches q_θ, the variational family is closed.
**Result.** Widening K from 75 to the full pool costs +0.141 Å; the shortlist binds, not the ranker. A perfect ranker inside the top-25 returns 2.609 Å (ORACLE, uses the native; `s17/results/inband.json`, reproduced 2.6087 in `s26/PROPOSAL_C.md` C1). 58 functionals of the distogram: best −0.009; the distance argmin inside its own top-25 is −0.014 against random. Greedy certifies 100% at 1,024 evaluations where the VQE needs 8,192. The refine_full arm gives 3.610 Å, +0.561 (`s17/LEDGER.md` L1 to L30). ESM's contact head is a filter, not a discriminator (+0.116 in-band, mostly compactness; `s17/LEDGER.md` L23).
**Closed.** Selection inside the pool, at four levels: 58 functionals, 29 in-band signals, 27 target-level calibration features, and the set-transformer (S12). Widening K.
**Retracted.** Nothing.
**Open.** Why the errors are maladaptive.

## Sprint 18: the degree-1 objective

**Question.** Is a degree-1 (residue-additive) objective in torsion space the missing structural prior, as a 19-target 2.411 Å result suggested? Are the distogram's confidences informative?
**Falsifier.** `s18/PREREG_exp.md` F1 to F5: the branch is closed if the full-instrument arm lands ≥ 3.5 Å, if the residue-additive minus full contrast includes zero, if the gain does not survive a random relabelling of the 2-bit torsion code, and so on. `s18/PREREG_phys.md`: "No β beats the permuted control" closes the physics-gate lane.
**Result.** All five falsifiers fired. Degree-1 on the full instrument is 4.335 Å, +1.287. The residue-additive object versus full is −0.004 [−0.390, +0.318]; the only interval excluding zero measures encoding luck, −0.440 [−0.844, −0.133], i.e. which 2-bit code names which torsion state (`s18/LEDGER.md` L5). Physics gates lose to random by +0.036 to +0.063; `leg_contact` +0.722. Shuffled residuals give 2.609 and isotropic noise 2.573 against the real 3.610 (both ORACLE arms, they use the native's residual field): the real errors are worse than random errors of the same size (`s18/LEDGER.md`). 1/sd-derived weights validated native-free, flat weights +0.153 [+0.066, +0.247] (`s18/results/q_ceilattack2.json`; L10).
**Closed.** Degree-1 objectives. Physics gates on the pool.
**Retracted.** L4's "weights are anti-informative" (amended by L10).
**Open.** The structure of the error: coherent, not random.

## Sprint 19: why the errors are maladaptive

**Question.** Why is the pool's error worse than noise of the same size?
**Falsifier.** `s19/PREREG_A.md`: if no native-free criterion beats the size-matched random repair by more than the MDE the repair lane closes; the whitening test (A3b) is decisive; "ρ ≤ 0.1 → the realisability defect does not price failure". `s19/PREREG_C.md`, `PREREG_D.md`: analogous.
**Result.** The pool has a systematic error direction that the predictor reproduces. Flipping the sign of the exact residual field gives 2.408 Å, −1.202, 128% of the gap; coherent errors at identical 0.688 sign accuracy emit +0.31 Å where i.i.d. errors emit −0.14 (`s19/LEDGER.md`; memory `error-coherence-decides-correctors`). The coherent component is shared across predictor families: 0.898/0.833/0.784/0.731/0.575/0.598 across pairs (`s19/LEDGER.md` L11). The shared-referent floor: two quantities measured against a common reference correlate by construction; measuring that floor first turned "two thirds sequence-independent" into "one fifth" (memory `shared-referent-floor`). A CA-only k = 300 arm is −0.110 but carries cis peptide bonds at +0.4242, an invalid structure (`s19/LEDGER.md`).
**Closed.** Correctors trained on the predictor's own features (they inherit its error structure).
**Retracted.** Nothing.
**Open.** Where the coherent error comes from (corpus, BLOSUM, or the score).

## Sprint 20: geometry of the surface, Legacy versus AMBER

**Question.** Is the pool a sampler problem (P1 to P3)? What do Legacy and AMBER measure? Does the fold-clustered CI change verdicts?
**Falsifier.** `s20/PREREG_A.md` §6 F-A1: "if the peptide-only pool's coherent error aligns with the current pool's at ρ ≥ 0.85, the corpus is not the carrier and this closes", in absolute and share-of-ceiling forms. `PREREG_B` to `PREREG_D`: analogous falsifiers per lane.
**Result.** P1, P2, P3 refuted: the K = 500 pool at 3.230 beats every sampler at 8,192 evaluations. Legacy is compactness (−0.4477 Å of radius of gyration, 124W/2L). AMBER's difficulty is a steric singularity; H_AMBER = E ∘ Relax50 explains 58/97/99% at three depths. ρ(Legacy, AMBER) = −0.0886. The Fréchet mean (the structure minimising the summed squared RMSD to all members, instead of the coordinate average) is +0.0048, null. On the continuous encoding the VQE trains: −0.210 [−0.339, −0.082], 5 of 5 folds; entanglement ablation −0.013 [−0.095, +0.077], MI 0.045 bits; α +0.012/+0.013/+0.025. The i.i.d. CI defect is programme-wide and the fold-clustered CI is adopted (`s20/LEDGER.md` L1 to L15, L-B; `s20/agentD_FINDINGS.md`).
**Closed.** The sampler hypothesis. Legacy in generation.
**Retracted.** Scope of the S16 "VQE is worse than not running it" (lattice-only).
**Open.** The exhaustive latent.

## Sprint 21: CVaR-VQE as a selector, and the MDE

**Question.** Enumerate the whole latent space the VQE searches. Does the exact argmin beat the pool? Does the readout entropy matter? Is the MDE constant right?
**Falsifier.** `s21/PREREG_A.md`: "if a Hamiltonian's tail arms do NOT beat their MATCHED-COUNT RANDOM subsets", the Hamiltonian carries no selection information. `s21/PREREG_B.md` F-E1: the gauge falsifier (an embedding at radius r and the angle chart must agree).
**Result.** For n ≤ 16 the register of 65,536 states is dense in 64 ms, so the budget exceeds the 2ⁿ latent on 75 of 126 targets; the exact argmin ties the pool (+0.025 at n = 75, +0.085 at n = 126) and is +1.7202 worse than the ORACLE over the same set, 0W/119L (ORACLE, uses the native); the latent ORACLE minus incumbent is −1.1766 (ORACLE). Readout entropy correlates −0.697 with RMSD while training entropy correlates +0.056; the AMBER readout correlation is −1.299 in the same units. The χ ladder orders nothing (0.37 Å spread); seed sd 0.200 Å. MDE = 2.8016 × SE per comparison; the 0.084 constant is wrong by 0.11× to 9.54× (`s21/LEDGER.md` L1 to L37, L8 for the MDE).
**Closed.** Search inside the latent (exhaustive). The MDE constant.
**Retracted.** The 0.084 Å MDE constant (every earlier "significant at 0.084" is re-read against its own SE).
**Open.** The readout-entropy lever (explained in S22).

## Sprint 22: the final campaign, routing

**Question.** Is there a per-target set size or scale a native-free router can find? Does the readout-entropy lever survive?
**Falsifier.** `s22/PREREG_A.md`: (a) soundness gate on the exact global argmin; (b) if gauge variation exceeds seed variation on a majority of targets the encoding is closed. `s22/PREREG_B.md` §4: single-number falsifier on the routed arm.
**Result.** The routing ceiling, 0.482, is an order statistic. The per-target ORACLE m is worth −0.239 and is real (it transfers split-half; ORACLE, uses the native), but four native-free routers are null, and a hedge is null. Ablation: D (distogram) best, T worst at +0.899. The readout-entropy lever replicates and is then explained by set equality: it is the classical ranking (`s22/LEDGER.md` L1 to L16). A finite-sample bound closes the rest.
**Closed.** Routers. The entropy lever as a quantum effect.
**Retracted.** Nothing.
**Open.** Per-target scale.

## Sprint 23: drive dev below 3.0

**Question.** Is there a per-target scale s* and can it be reached native-free? Is the pool's error common-mode?
**Falsifier.** `s23/PREREG_A.md`: "No achievable rule beats `avg_75` past its own MDE with a CI excluding zero at ANY" size, the lane closes; `clust_m`'s optimal m must differ from 75 to matter. `s23/PREREG_C.md`: primary is `W_trained − avg_matched_trained`; the routed arm must beat fixed m = 75 past its own MDE with both CIs.
**Result.** The per-target ORACLE scale s* is worth −0.3403, 126W/0L, and transfers at 98.9% (ORACLE, uses the native), and the placebo shows it is unreachable in principle: no native-free quantity tracks it. The pool's error is 68% common-mode (fraction 0.676, 50.7× i.i.d.), so averaging removes at most 32% (`s23/LEDGER.md` L9). Uniform weights are optimal; removing four outliers costs +0.142. Repair: 17 of 17 arms at or worse than nothing. Probability-weighted readout +0.0213 null; T = 0 +0.0386 worse (`s23/LEDGER.md` L8). CA-only arms carry cis bonds 0.51 to 0.64 (`s23/LEDGER.md` L1 to L11).
**Closed.** Per-target scale. Repair. Weighting.
**Retracted.** Nothing.
**Open.** Learned candidate generation.

## Sprint 24: learned candidate generation, and the prior ladder

**Question.** Can a generator make better candidates than retrieval? Where does the common-mode error come from? How steep is the prior?
**Falsifier.** `s24/PREREG_A.md`: "zero source chains in any of the 126 universes match a benchmark or dev-set target under any of the five criteria" (contamination gate). `s24/PREREG_C.md` §4: "Every matched-size mixture at or above the pure incumbent" closes the lane. `s24/PREREG_D.md`: subset-hood falsifier, 0 of 3,888 fired.
**Result.** The generator spec was written and retired: every union is null-to-worse and selection's lift is universal (+0.166, L9-A). The shared bias is the score's (L3: selection makes bias parallel, cos 0.9432). The corpus is small, not repetitive: 7,329 windows at L = 16, 2.2% β, 410 peptides (L8). The prior-attribution ladder: −2.1496 Å per unit γ toward a perfect prior at the origin, −0.4168 at the top, concave; at γ = 1 the same pipeline gives 2.2261 Å (ORACLE, uses the native); a MASS variant at γ = 0.1 gives 2.8334, −0.2150, 2.4× MDE, 113W/13L (ORACLE) (`s24/LEDGER.md` L13, L13-A; `s24/results/priorladder.json`). The fitted AMBER weight ceiling is 0.0148 Å; nested CV picked w = −0.5 (L16). Grid oracles are order statistics; only split-half transfer arms survive (L15, L15-A). A leakage defect reaching the benchmark is declared (L4).
**Closed.** Learned generation. The functional lever in all its forms.
**Retracted.** L5's mechanism claim (withdrawn by L11: a generic prior reproduces 102% of β). L15 corrected by L15-A (two of four oracles survive). L9 corrected by L9-A (+0.166 not +0.20). L8 amended by L8-A. L10 amended by L10-A (AMBER selects along a non-parallel direction; it cannot pay for it).
**Open.** Only the prior.

## Sprint 25: the endgame

**Question.** Fix the reporting basis, verify the quantum component to machine precision, run the physics suite, freeze the architecture, build the results lab.
**Falsifier.** `s25/PREREG_Q.md`: A1 any MPS amplitude differing from dense by more than 1e-12; A2 any parameter-shift gradient disagreeing with finite differences beyond truncation error or any sampling code path; A3 one subset violation over adversarial families with exact zeros; A4 the 0.7× MDE rule absent. `s25/PREREG_PHYS.md` §4: F0 cache integrity at relative difference 0.0; F1 the anchor, C3 under classical top-75 must return 3.0483 to four decimals; F2 subset-hood on 882 cells; F3 the registered prior "I expect every AMBER configuration to fail", with the refutation condition written out.
**Result.** L1 the posterior is over-confident by 2× (z sd 1.9962). L2/L7 calibrating it does not help (L7 retracts L2's mechanism). L3/L5 the "+0.113 Å" is withdrawn. L4/L8/L11 the basis decision: built chain, +0.1664 against the point cloud at 3.26× MDE, and a synthetic leaderboard found on disk was quarantined (a `legacy_amber_distogram` at 2.0921 was fake). L6 the score's target is quantised to 17 values. L9 the "86% accounted" null proven to carry no information. L10 the round-trip gate. L12 Phase I closes: moving the score's target 25% toward truth along a realistic direction does not move the endpoint (51 arms, correlation +0.054; +0.024 at cos 0.5). L13 three sprints of record untracked (fixed S26). L14 tests green. L15 the readout is insensitive to half the mass. L16 the seven-configuration suite (Part IV.4). L17 two trained states (Part V.7). L18 Form 5 closed. None of the four Q falsifiers fired; F3's registered prior held (`s25/LEDGER.md` L1 to L18; `s25/results/*.json`).
**Closed.** The quantum component's accuracy claim. Both physics energies as rankers. The functional lever in all five forms.
**Retracted.** L2's mechanism (by L7), the "+0.113" (L3, L5), the 3.0483 as the headline (L4, L8), the s24 REV1 prefix claim (by `q_verify.json`), and L11's own guard ("necessary and not sufficient").
**Open.** The prior. Publishing the trainability chain. The three directions in Part X.5.

## Sprint 26, for context only

S26 is outside this report's scope and was in progress when it was written. It has so far: a governor for the machine (`s26/governor.py`), a Phase 0 examination that reproduced every headline bit-for-bit (`s26/LEDGER.md` L28, L31 to L33), the DLA and ADAPT variance results (Part V.10, V.12), a cis census (no cis peptide bond anywhere on the instrument, L22), the finding that production relaxation is worse than a random move of its own size (L39), a feasibility stop on ESMFold (B1: no checkpoint, no `openfold`, 8.76 GB against 4.4 GB headroom; `s26/PROPOSAL_B.md`), and the beginning of a learned-prior ladder (`s26/PROPOSAL_C.md`).


---

# PART VII. WHAT HAS BEEN ESTABLISHED

These are the statements the record supports, each with its strongest artefact. They are ordered from the most important to the presenter downward.

**1. The binding constraint is the distance prior, and it is the only steep lever.** Moving the prior toward the truth is worth −2.1496 Å per unit at the origin of the attribution ladder and −0.4168 at the top; the same pipeline with a perfect prior reaches 2.2261 Å (ORACLE, uses the native) (`s24/results/priorladder.json`; `s24/LEDGER.md` L13). Every downstream lever measured is flat to a few hundredths: re-consuming the posterior in eleven ways (memory `score-axis-does-not-transfer`), calibrating it (`s25/LEDGER.md` L2/L7), moving its target 25% toward truth along a realistic direction (+0.024 at cos 0.5; `s25/LEDGER.md` L12).

**2. Selection inside the pool is closed at four levels.** 58 functionals, 29 in-band signals, 27 target-level calibration features (`s17/LEDGER.md`), and a set-transformer with a flat learning curve (`s12/results/agg_dec_v2*.json`). A perfect ranker inside the shipped top-25 returns 2.609 Å (ORACLE, uses the native; `s17/results/inband.json`). The terminal operator consumes the set mean, d_out = 1.16 × mean + 0.04 × best (`s12/SPRINT12_DOSSIER.md`).

**3. The pool's error is 68% common-mode.** Fraction 0.676, 50.7× the i.i.d. prediction (`s23/LEDGER.md` L9). Averaging removes at most the other 32%. The coherent component is shared across predictor families (0.575 to 0.898; `s19/LEDGER.md` L11) and is installed by the score's selection (cos 0.9432; `s24/LEDGER.md` L3).

**4. Neither physics energy ranks nativeness on these pools.** Legacy +0.3302 and AMBER +0.4554 worse than a random 75-subset; rank-permuted channels reproduce both; inside a distogram-led configuration AMBER is +0.0331 at 0.43× MDE (`s25/results/phys_suite.json`). Legacy is compactness (`s20/LEDGER.md`). The fitted AMBER weight ceiling is 0.0148 Å (`s24/LEDGER.md` L16).

**5. The quantum component is exact and contributes nothing measurable to accuracy.** Verified: MPS 3.331e-16, statevector 5.551e-17, parameter-shift cosine 1.000000000, 0 subset violations on 2,592 cells (`s25/results/q_verify.json`). Contribution: −0.1405 at 0.68× MDE, NOT A RESULT (`s25/results/q_alpha.json`); matches classical top-m at 0.0000 in all seven suite configurations (`s25/results/phys_suite.json`, `controls_rank.*.vs_topm`). The "+0.113 Å" is withdrawn (`s25/LEDGER.md` L3, L5).

**6. The Hamiltonian barely changes between targets.** Rank-ladder deviation at most 1.18% of range (`s25/results/q_gibbs.json`); the circuit trains toward two states (`s25/LEDGER.md` L17).

**7. The trainability chain is exact and measured.** Locality theorem (agreement 1.0000; `s13/results/qarch_locality_geom.json`); Pauli weights Legacy 2.236, AMBER 3.015, 79 of 79 cells (`s13/results/geo_pauli.json`); raw AMBER is a delta spike (`s13/results/walsh_xval.json`); spectrum × kernel predicts variance to 0.4% (`s13/results/walsh_predict.json`); depth-3 width slopes −0.047 to −0.649 log₂ per qubit (`s25/results/q_plateau.json`); DLA = so(2ⁿ) from depth 2 at n = 7 (`s26/results/q_dla.json`); ADAPT-grown circuits at α = 1 are product circuits (`s26/results/q_var.json`).

**8. The reporting basis is the built chain.** 3.2148 Å; the point cloud 3.0483 is an intermediate labelled illegal in the code (`bench_results/compare_tuning126.json`; `core/bench.py` line 682; `s25/LEDGER.md` L4, L8, L11).

**9. The tuning gain did not transfer to the benchmark.** +0.0103 [−0.1596, +0.1803] on the 60 (`s9/final_report.json` via `verify/headline_audit.json`). Ten targets carried 61% of the tuning gain (`docs/FINDINGS.md` S9-10).

**10. The instrument's own rules are results.** MDE per comparison at 2.8016 × SE (`s21/LEDGER.md` L8); fold-clustered CI beside i.i.d. (`s20/LEDGER.md` L-B); controls in the operator's space; the zero-information control is a helix, not a uniform draw; grid oracles are order statistics (`s24/LEDGER.md` L15-A); tie-breaking by array order invents winners (`docs/FINDINGS.md` S8-8).

**11. Sequence carries little structural information at this length, and is harmful on the failures.** Blind pipeline 3.989 vs 3.213; on FAIL18 blind 5.425 beats shipped 6.019 (`s12/SPRINT12_DOSSIER.md`). φ 36.1° vs 36.4° (`s13/cache/tors_rows.npz`).

**12. No fresh benchmark exists.** All 204 clusters of 9 to 16-mers are spent; the fresh world supply is 16 targets, 10 amyloid (`s12/SPRINT12_DOSSIER.md`).

**13. Everything reproduces.** Bit-identical between worker counts (`bench_results/compare_tuning126.json`), fresh re-run at 0.0 (`s26/results/e_reproduce.json`), projection 126/126 bit-identical (`verify/project_exactness.json`), leak audit clean (`verify/leak_audit.json`), 370 tests 357 passed 13 skipped (`s26/results/test_run.json`).

---

# PART VIII. CLOSED AND OPEN

## VIII.1 Closed (with the closing artefact)

| question | closed by | artefact |
|---|---|---|
| Search harder | S5-2, S12 QUBO, S21 exhaustive latent | `docs/FINDINGS.md` S5-2; `s12/SPRINT12_DOSSIER.md`; `s21/LEDGER.md` |
| Learn a ranker inside the pool | S6, S12, S17 (four levels) | `s12/results/agg_dec_v2*.json`; `s17/LEDGER.md` |
| Widen or shrink K globally | S7-6, S17 | `s17/LEDGER.md` (+0.141) |
| Per-target m or scale by a native-free router | S22, S23 | `s22/LEDGER.md`; `s23/LEDGER.md` |
| Repair or steer the average native-free | S16, S19, S23 (17 of 17) | `s16/LEDGER.md`; `s23/LEDGER.md` |
| AMBER relaxation as an accuracy step | S16 L27, S26 L39 | `s16/LEDGER.md`; `s26/LEDGER.md` L39 |
| Legacy or AMBER as a ranker | S25 L16 | `s25/results/phys_suite.json` |
| A native-free switch between physics configurations | S25 L18 | `s25/results/phys_form5.json` |
| Fusion of decorrelated channels | S16 (Krogh and Vedelsby) | `s16/LEDGER.md` |
| Post-hoc de-biasing or calibration of the prior | S7-3, S15, S25 L2/L7 | `docs/FINDINGS.md`; `s25/LEDGER.md` |
| Re-consuming the posterior (eleven ways) | memory `score-axis-does-not-transfer`; S25 L12 | `s25/LEDGER.md` L12 |
| Chemical shifts | S14 | `s14/LEDGER.md` |
| Torsion space as the route to 2 Å | S13 | `s13/SPRINT13_DOSSIER.md` |
| Degree-1 objectives | S18 | `s18/LEDGER.md` L5 |
| Learned candidate generation | S24 | `s24/LEDGER.md` L9, L9-A |
| Evolutionary information, loops, retrieval by predicted structure | S9-7, S9-4, S5-13 | `docs/FINDINGS.md` |
| The CVaR-VQE's accuracy contribution | S25 L3/L5; suite `vs_topm` | `s25/results/q_alpha.json`; `phys_suite.json` |
| Entanglement as a lever | S20 | `s20/LEDGER.md` |
| QNG | S13 | `s13/SPRINT13_DOSSIER.md` |
| A fresh benchmark | S12 | `s12/SPRINT12_DOSSIER.md` |
| ESMFold on this machine | S26 B1 | `s26/results/b1_feasibility.json` |

## VIII.2 Open (with what would settle it)

| question | state | what settles it |
|---|---|---|
| A better distance prior from inputs this machine can compute | S26 Proposal C ladder running | `s26/results/p_ladder_report_*.json`, a rung beating the shipped prior on the fold-clustered CI |
| The target-specific third of the coherent error | S19 L11 shows two thirds is generic | a predictor whose error decorrelates from the generic reference |
| Per-target sign at inference for in-band ordering | 0.986 within-target, 0.600 across (memory `in-band-ordering-is-per-target`); compactness proxies 0.24 to 0.37 | a native-free proxy above 0.638 |
| The prior/pool disagreement signal | demonstrated, Å value not measured (memory `prediction-pool-disagreement-is-a-native-free-signal`) | a paired n = 126 run |
| The ADAPT endpoint (A1) | pending phase gate, expected null | `s26/results/a1/*.json` scored under `s26/PREREG_A1.md` |
| The n = 6, n = 9 DLA subalgebra | HYPOTHESIS from two cases (`s26/LEDGER.md` L27) | a proof or a third width |
| Publishing the trainability chain | outline exists (`s26/PROPOSAL_B_REPLACEMENT.md`) | a written paper |
| The leakage defect in `identity()` | declared, fix ships dark | re-clustering with pinned folds regenerated and every fold model retrained |
| The 3.2126 vs 3.2148 route difference | Appendix D, D2 | one projection route named as canonical in the results lab |

---

# PART IX. OPERATING MANUAL

## IX.1 Reproduce the headline

    python -m core.pipeline run --manifest tuning126 --workers 6
    python -m core.pipeline stats --manifest tuning126

or, to time both arms and compare them,

    python -m core.bench --arm baseline  --manifest tuning126 --fresh
    python -m core.bench --arm optimised --manifest tuning126 --workers 6 --fresh
    python -m core.bench --compare bench_results/baseline_tuning126.json bench_results/optimised_tuning126.json

(`README.md`, "Running a full-quality experiment"; manifests are `smoke8`, `smoke24`, `tuning126`, `dev24`, `benchmark60`). The run writes per-target records to `bench_results/cache/1fc9f2dcf489e2fb/` and an aggregate whose `science.rmsd_arm` must read 3.214765154210998 and `science.rmsd_avg` 3.048338093879532 (`bench_results/compare_tuning126.json`). The one-worker and eight-worker aggregates must agree to ±4.4e-16. `verify/run_equiv2.sh` runs the baseline-versus-shipped equivalence and must report bit-identical 8/8 (`s26/LEDGER.md` L19).

## IX.2 Rebuild the results lab

    python -m s25.resultslab.build --mode frozen --spec s25/results/real_pools/spec.json

regenerates `results/summary/leaderboard.csv`, `results.csv`, `target_map.json`, the structures under `results/structures/`, and `results/site/index.html` (`s25/resultslab/build.py`; `s25/results/frozen_build.log`). The build requires provenance on every input and passes the pool gate (WARN on two targets, explained in output since S26 L7) and the difficulty gate (`s25/resultslab/schema.py`).

## IX.3 Run the tests

    pytest tests/

370 tests, 357 pass, 13 skip. Set `VERIFY_SLOW=1` for the 11 opt-in slow tests. Run `tests/test_amber.py` and `tests/test_amber_frame_invariance.py` as separate jobs so that only one AMBER process is live (`s26/TEST_RUN.md`). Peak RSS of the non-AMBER files is 1.692 GB.

## IX.4 Do not

- Do not run anything on `results/benchmark_manifest.json`. The CLI refuses without `--i-am-spending-the-benchmark`, and there is no second benchmark (`core/bench.py` line 1068; `s12/SPRINT12_DOSSIER.md`).
- Do not regenerate `peptide_folds.json` or `peptide_clusters.json`. Correcting the identity clustering once moved 13 benchmark targets and invalidated every fold model (memory `benchmark-and-folds-must-be-pinned`).
- Do not quote a point-cloud number as a structure. The built chain is the result (`s25/LEDGER.md` L8).
- Do not quote a mean without its SE and its own MDE (`s21/LEDGER.md` L8).
- Do not run more than one AMBER job at a time on this machine; the governor bands are 88/90/93/95% (`s26/governor.py`; memory `box-baseline-load-is-67-percent`).
- Do not use `np.argmin` on a tied signal without averaging over the tied set (`docs/FINDINGS.md` S8-8).
- Do not cite a memory index line without reading the memory body (memory `read-the-memory-body-not-the-index-line`).

## IX.5 Where things are

| thing | path |
|---|---|
| frozen architecture | `ARCHITECTURE.md` (a15406c) |
| production code | `core/pipeline.py`, `predict.py`, `data.py`, `project.py`, `quantum.py`, `amber.py`, `energy.py`, `geometry.py`, `bench.py` |
| instrument | `s12/instrument.py` |
| per-target production cache | `bench_results/cache/1fc9f2dcf489e2fb/<PDB>.json` |
| aggregate | `bench_results/compare_tuning126.json` |
| results lab | `s25/resultslab/`, outputs in `results/summary/` and `results/structures/` |
| audits | `verify/*.json` |
| S25 quantum and physics artefacts | `s25/results/q_verify.json`, `q_gibbs.json`, `q_alpha.json`, `q_plateau.json`, `phys_suite.json`, `phys_form5.json`, `calib.json` |
| ledgers | `s14/` to `s26/LEDGER.md`; `docs/FINDINGS.md` for S5 to S11; `s12/`, `s13/` dossiers |
| state brief | `docs/STATE_BRIEF_2026-09-12.md` |


---

# PART X. PRESENTING THIS PROJECT

## X.0 Who you are presenting to

Mafalda Ramôa works in the group of Sophia Economou and Edwin Barnes at Virginia Tech, the group that introduced ADAPT-VQE (Grimsley, Economou, Barnes, Mayhall, Nature Communications 2019) and qubit-ADAPT-VQE (Tang et al., PRX Quantum 2021). Public listings at the time of writing describe her as a graduate student in the Barnes group supervised by Barnes and Economou, with co-supervisors Luís Paulo Santos (University of Minho) and Ernesto Galvão (INL, Portugal); your email calls her a postdoc, and either way she is the person in that group who has spent the most time making ADAPT-VQE cheaper to run. Her papers that matter for your ten minutes:

- **Coupled Exchange Operators (CEO) pool**: "Reducing the Resources Required by ADAPT-VQE Using Coupled Exchange Operators and Improved Subroutines" (Ramôa, Anastasiou, Santos, Mayhall, Barnes, Economou; arXiv:2407.08696, 2024). A new operator pool that cuts CNOT count, CNOT depth and measurement cost by up to 88%, 96% and 99.6% on 12 to 14 qubit molecules (LiH, H₆, BeH₂), and beats fixed ansätze such as UCCSD.
- **Hessian recycling**: "Reducing measurement costs by recycling the Hessian in adaptive variational quantum algorithms" (Quantum Science and Technology 10, 015031, 2025; arXiv:2401.05172). Reuses second-derivative information across ADAPT iterations so the optimiser needs fewer measurements.
- **Gradient troughs**: "Strategies for Overcoming Gradient Troughs in the ADAPT-VQE Algorithm" (Stadelmann, Übelher, Ramôa, Sambasivam, Barnes, Economou; arXiv:2512.25004, December 2025). Gradient troughs are ADAPT's version of a flat region: the operator-pool gradients become very small before the energy minimum is reached. The paper diagnoses them during a run and proposes where to insert operators, using the non-commutative algebra of the ansatz.
- **Co-ADAPT-VQE**: "Co-Designed Adaptive Quantum State Preparation Protocols" (Ramôa, Santos, Mayhall, Barnes, Economou; arXiv:2601.20681, January 2026). Grows the ansatz with the hardware's connectivity in the loop; up to 97% CNOT reduction on linear nearest-neighbour devices, over 70% even with all-to-all connectivity.
- Group background you should know exists: the group's 2023 npj Quantum Information paper (Grimsley et al.) arguing that ADAPT-VQE is insensitive to rough parameter surfaces and to barren plateaus, which is the group's claim about trainability that your Part V.9 to V.12 results speak to directly.

What this means for you. She will not be impressed by "we used a VQE". She will be interested in (a) whether your circuit is exact and verified, (b) what your gradient-variance measurements say and what scope they carry, (c) the DLA and product-circuit observations from S26, because they are exactly the kind of algebraic statement her group makes, and (d) whether you know the difference between a gradient trough, a barren plateau, and a small gradient from a trivial circuit. You do. Say so with the numbers.

(Sources for this section: the INSPIRE author page `inspirehep.net/authors/2613045`; the Barnes group pages `www1.phys.vt.edu/~efbarnes/`; arXiv abstracts 2407.08696, 2401.05172, 2512.25004, 2601.20681. These are the only numbers in this document not read from a repository artefact; they are listed in Appendix B under "external".)

## X.1 The ninety-second statement

"I built a pipeline that predicts the backbone shape of short peptides, nine to sixteen residues, from sequence. It retrieves five hundred real fragments, scores them against a learned distance prior, keeps the best seventy-five, averages them and projects the average onto a chain with ideal geometry. On a hundred and twenty-six held-out peptides it gets a mean Cα RMSD of 3.21 Å, a median of 2.97, best 0.18 Å, with about half the targets under 3 Å. That is not state of the art and I do not claim it is.

What I think is worth your time is two things. First, the instrument: paired statistics on all one hundred and twenty-six targets, a per-comparison minimum detectable effect, fold-clustered confidence intervals, pre-registered falsifiers, a sealed benchmark I spent once, and a corrections ledger. That instrument closed about thirty ideas, including most of my own early claims. Second, the quantum component: a seven-qubit, twenty-one-parameter RY-CNOT ansatz trained by Adam on a CVaR free energy over a diagonal Hamiltonian. I verified it to machine precision, exact parameter-shift gradients against finite differences, an exact MPS twin, and a proved theorem that its CVaR tail is always a subset of the classical ranking's prefix. Then I measured what it contributes to accuracy. Nothing, at 0.68 of my minimum detectable effect. What survives is a chain of exact statements: chain geometry fixes the locality of the energy, that fixes its Pauli spectrum, the spectrum times the circuit's gradient kernel predicts the measured gradient variance to 0.4%, the gradient variance at depth three decays at about a quarter of a bit per qubit, the dynamical Lie algebra is the full so(2ⁿ) from depth two, and an ADAPT-grown circuit on this Hamiltonian is a product circuit at α = 1. I would like to talk about what an adaptive ansatz is for when the target state is a product state, and whether the same measurement chain would say something on a Hamiltonian that is not diagonal."

Every number: 3.21 (`bench_results/compare_tuning126.json` rmsd_arm 3.2148); 2.97 (`results/summary/leaderboard.csv` median 2.9661); 0.18 (T030, `results/summary/results.csv`); half under 3 Å (0.5079, `leaderboard.csv`); 126 (`s12/instrument.py`); 7 qubits, 21 parameters (`s25/results/q_verify.json`); 0.68× MDE (`s25/results/q_alpha.json`); 0.4% (`s13/results/walsh_predict.json` via `s26/PROPOSAL_B_REPLACEMENT.md`); a quarter of a bit per qubit (−0.2429 to −0.3105, `s25/results/q_plateau.json`); so(2ⁿ) from depth 2 (`s26/results/q_dla.json`); product circuit at α = 1 (`s26/results/q_var.json`).

**Sixty-word version.** "A retrieval-and-average pipeline for 9 to 16-residue peptides, 3.21 Å mean over 126 held-out targets, with an exact 7-qubit CVaR-VQE that I verified to machine precision and then showed contributes nothing measurable, because its tail is provably a subset of the classical ranking. The result I would publish is the measured chain from force-field locality to gradient variance, plus the DLA and ADAPT product-circuit findings."

## X.2 The discrepancy with the first email, in one sentence

The first email described a 36-qubit lattice model aimed at chignolin and trpzip with a sub-2 Å target; what exists now is a 7-qubit (9 in the physics suite) exact statevector selector over retrieved real fragments, measured on 126 peptides none of which is chignolin or trpzip, at 3.21 Å mean, and the sub-2 Å target was never reached on any instrument (neither the <3.0 Å dev target nor the <2.5 / <2.0 Å targets set in S15 to S25; `docs/STATE_BRIEF_2026-09-12.md` §2).

Every number that differs:

| the email said | what the record says | artefact |
|---|---|---|
| 36 qubits | 7 qubits, 3 layers, 21 parameters in production; 9 qubits in the S25 suite | `core/pipeline.py` `Config.vqe_qubits=7`; `s25/results/q_verify.json`; `s25/LEDGER.md` L16 |
| chignolin, trpzip | 126 cluster-disjoint peptides, 9 to 16 residues; neither chignolin (1UAO) nor trpzip is on the instrument | `s12/instrument.py`; `results/summary/target_map.json` (126 entries, no 1UAO/1LE0/1LE1) |
| sub-2 Å | 3.2148 Å mean built chain; 28.6% of targets under 2 Å | `bench_results/compare_tuning126.json`; `results/summary/leaderboard.csv` |
| a lattice model | real fragment retrieval, K = 500, BLOSUM62 | `core/pipeline.py` `retrieve` |
| AMBER anti-ranking on chignolin | AMBER measurably worse than noise as a ranker on 126 real pools; converged interaction-only protocol on early decoy banks does not transfer to real pools | `s25/results/phys_suite.json`; memory `decoy-bank-not-a-pool-proxy` |
| VQE as the folding engine | VQE off in production; selector reproduces the classical top-75 | `Config.quantum=False`; `results/summary/results.csv` `selector` column |

Say it once, early, without apology: "The project I emailed you about in the summer is not the project I am showing you. Here is what changed and why."

## X.3 Slide by slide

The deck `vqe_research_overview.pptx` is not on disk or in git (Appendix D, D4; `s26/LEDGER.md` L2). The campaign brief describes its structure: title; six "what I built" slides (why peptides are hard; the pipeline; results with structure overlays; the verified quantum component; the barren-plateau result; where the remaining accuracy lives); three "where I want to go" slides (ADAPT-VQE; AlphaFold gaps; learned model plus physics); goal and ask (`s26/briefs/PR.md` §1). The annotation below follows that structure. If your file differs, map by title.

**Slide 1, title.** Say your name, the one-sentence description from X.1's sixty-word version, and that everything on the slides has an artefact path in a document you can send. Thirty seconds.

**Slide 2, why peptides are hard.** Say: few contacts, few relatives, flexible in solution. Numbers you may use: the best-matching library window has about 12% identity (`docs/FINDINGS.md` S5-9); the distance prior is over-confident by 2× (z sd 1.9962, `s25/results/calib.json`); φ carries no sequence signal at this length (36.1° vs 36.4°, `s13/cache/tors_rows.npz`). Do not say "AlphaFold fails on peptides" unless you have measured it; you have not (S26 B1 could not run ESMFold on this machine, `s26/results/b1_feasibility.json`).

**Slide 3, the pipeline.** Nine boxes from Part III. Say the sentence "the terminal operator consumes the set mean" (`s12/SPRINT12_DOSSIER.md`) because it is the one non-obvious design fact and it explains why the quantum selector cannot matter. If the slide shows the quantum stage, label it "off in production".

**Slide 4, results with structure overlays.** The three overlays: 1S9Z 0.18 Å (0.18242, `results/summary/results.csv` T030), 2MJQ 0.40 Å (0.39949916705333455, T053), 6WPB 0.42 Å (0.4167616469202071, T098). The summary numbers: mean 3.21 (3.2126 on the results-lab rebuild, 3.2148 on the production cache; say "3.21" and know both), best 0.18, 28.6% under 2 Å (0.2857), 50.8% under 3 Å (0.5079) (`results/summary/leaderboard.csv`, production row). Random-75 control 3.43 (3.4251, `s25/results/phys_suite.json` `random_null_rank.mean`; point-cloud basis, so strictly it is the control for 3.0483 not 3.2148; say "the random-subset control on the same basis is 3.43"). Say "these three are the best three; the worst is 8.23 Å on 2MQ2 and the median is 2.97" (`leaderboard.csv`). Showing only the best three without saying so is the kind of selective reporting this project's ledger exists to prevent.

**Slide 5, the verified quantum component.** The verification numbers: statevector vs reference 5.551e-17; parameter-shift vs finite difference cosine 1.000000000, relative error 4.597e-10; free-energy gradient 4.663e-10; MPS vs dense 3.331e-16, χ = 2^layers exactly; no sampling in the trained path; set-equality 0 violations on 2,592 cells, 972 of 972 exact equalities on full support (`s25/results/q_verify.json`). Then the contribution: −0.1405 at 0.68× MDE, NOT A RESULT (`s25/results/q_alpha.json`). Say both. The second sentence is what makes the first one credible.

**Slide 6, the barren-plateau result.** Title it "gradient variance at depth 3", not "no barren plateau". Show the five slopes (−0.6492, −0.2522, −0.0472, −0.3105, −0.2429 log₂ per qubit, n = 4 to 13; `s25/results/q_plateau.json`) and the depth saturation at n = 7 (2.381e-2 at 1 layer to 8.2e-3 at 12; same file). Say: "at depth three, with 3n parameters against dim so(2ⁿ) = 8128 at n = 7, the ansatz is nowhere near a 2-design, and the DLA is the full so(2ⁿ) from depth two, so nothing algebraic protects it at scale; this is a shallow-regime measurement" (`s26/results/q_dla.json`; `s26/LEDGER.md` L27). Then the S26 ADAPT result in one line: "an ADAPT-grown circuit on this Hamiltonian at α = 1 selects one to three distinct operators and is a product circuit; its variance does not decay because there is nothing to train" (`s26/results/q_var.json`; L35). That line is written for this audience.

**Slide 7, where the remaining accuracy lives.** The ORACLE ladder (Part II.3), every rung labelled ORACLE on the slide: pool best 1.71, top-75 best 2.31 (`bench_results/compare_tuning126.json`), perfect prior 2.23 (`s24/results/priorladder.json`), hull floors 1.840 and 0.853 (`docs/FINDINGS.md` S10-5). The headroom numbers 1.36 Å inside the top-75 and 2.35 Å inside the 500 are shipped-minus-ORACLE and are diagnostics (`docs/FINDINGS.md`, bottom line). The sentence: "the prior is the only steep lever, −2.15 Å per unit toward truth at the origin; everything downstream is flat to hundredths" (`s24/LEDGER.md` L13). And: "the pool's error is 68% common-mode, so averaging cannot fix it" (`s23/LEDGER.md` L9).

**Slide 8, ADAPT-VQE.** See X.5, direction A, and X.8. Do not put a predicted RMSD on this slide. Put the A2 and A4 measurements and the A1 status (pending, expected null).

**Slide 9, AlphaFold gaps.** See X.5, direction B. State that ESMFold did not run on your machine (8.76 GB against 4.4 GB headroom; no checkpoint; `openfold` absent; `s26/results/b1_feasibility.json`), so this direction is un-measured, and say what you would measure first.

**Slide 10, learned model plus physics.** See X.5, direction C. The physics half has a measured answer: worse than noise as a ranker (`s25/results/phys_suite.json`), and a validity role only (`s26/LEDGER.md` L39). The learned-prior half is the open question and the S26 ladder is running (`s26/PROPOSAL_C.md`).

**Slide 11, goal and ask.** The ask is an independent study with a defined measurement (X.8). Say which direction you favour and why in two sentences, and end with a question for her.

## X.4 Flashcards

Front / back. Learn all of them.

1. Mean built-chain RMSD, 126 targets / 3.2148 Å (`bench_results/compare_tuning126.json` rmsd_arm); leaderboard rebuild 3.2126 (`results/summary/leaderboard.csv`).
2. Point-cloud mean / 3.0483 Å, an intermediate, "raw average (illegal)" (`core/bench.py` 682).
3. Median, best, worst / 2.9661, 0.18242 (1S9Z), 8.234 (2MQ2) (`leaderboard.csv`).
4. Fraction under 2 Å, under 3 Å / 28.57%, 50.79% (`leaderboard.csv`).
5. Benchmark result / +0.0103 [−0.1596, +0.1803], 31W/29L, n = 60, spent once (`s9/final_report.json` via `verify/headline_audit.json`).
6. Tuning gain that did not transfer / −0.250 [−0.383, −0.117] (`docs/FINDINGS.md` S9-10).
7. MDE rule / 2.8016 × SE per comparison; below 0.7× is NOT A RESULT (`s21/LEDGER.md` L8; `s25/results/q_verify.json`).
8. Why fold-clustered CI / each fold model saw about 100 of the other 125 natives (`s20/LEDGER.md` L-B).
9. K, m, qubits, layers, parameters, iterations / 500, 75 (128 with quantum), 7, 3, 21, 50 (`core/pipeline.py` Config; `q_verify.json`).
10. Objective / F = CVaR_α − T·H; VQE_LFO T = 0.3 on all folds, α = 1 on three of five (`core/pipeline.py`; `s25/LEDGER.md` L3).
11. Quantum contribution / −0.1405 vs argmin, 0.68× MDE, UNDERPOWERED (`s25/results/q_alpha.json`).
12. The withdrawn number / +0.113 Å; measured at T = 0.1 which does not ship; a marginal, not paired; sign reverses at T = 0.3 (`s25/LEDGER.md` L3, L5).
13. Set-equality theorem / tail ⊆ prefix; 0 violations on 2,592 cells; 972/972 equal on full support (`q_verify.json`).
14. Two trained states / E is a rank ladder to 1.18% of range; the circuit sees the same Hamiltonian on every target (`q_gibbs.json`; `s25/LEDGER.md` L17).
15. Gradient verification / cosine 1.000000000, rel err 4.597e-10 (`q_verify.json`).
16. Tail-baseline defect / cosine 0.566586 on the deployed instrument, 0.6847 on the audit instrument; original price −0.023; default is const (`q_verify.json`; `verify/cvar_audit.json`; `docs/FINDINGS.md`).
17. Width slopes / −0.649 linear; −0.252, −0.047 CVaR; −0.311, −0.243 deployed; depth 3, n = 4 to 13 (`q_plateau.json`).
18. DLA / so(2ⁿ) from depth 2 at n = 7 (dim 8128); abelian at depth 1; subalgebra only at n = 6, 9 (`s26/results/q_dla.json`).
19. ADAPT on this Hamiltonian / product circuits at α = 1, no decay, nothing to train; L2 pool at α = 0.25 decays −0.302 like the fixed −0.243 (`s26/results/q_var.json`).
20. Pauli weights / Legacy 2.236, AMBER 3.015, 79/79 cells; raw AMBER a delta spike, 99.6% in top-10 (`s13/results/geo_pauli.json`, `walsh_xval.json`).
21. Spectrum predicts variance / to 0.4% (`s13/results/walsh_predict.json`).
22. Locality theorem / d_ij depends on exactly the j − i − 1 residues between; agreement 1.0000 (`s13/results/qarch_locality_geom.json`).
23. Physics suite / Legacy +0.3302 and AMBER +0.4554 worse than random-75; permuted channels reproduce it (`phys_suite.json`).
24. Prior ladder / −2.1496 Å per unit at origin; perfect prior 2.2261 (ORACLE) (`s24/results/priorladder.json`).
25. Common mode / 68% (0.676), 50.7× i.i.d. (`s23/LEDGER.md` L9).
26. Operator law / d_out = 1.16 × mean + 0.04 × best (`s12/SPRINT12_DOSSIER.md`).
27. Calibration / z sd 1.9962; cov90 0.616; multimodal 24.1% (`calib.json`).
28. Leak price / +0.0004 tuning, +0.0030 benchmark; 13/126 at ≥0.6 identity (`docs/FINDINGS.md` S10-4, S9-10).
29. Tests / 370, 357 pass, 13 skip (`s26/results/test_run.json`).
30. Reproduction / bit-identical 1 vs 8 workers; fresh re-run 0.0 (`compare_tuning126.json`; `s26/results/e_reproduce.json`).
31. The three overlay targets / 1S9Z 0.18, 2MJQ 0.40, 6WPB 0.42 (`results.csv`).
32. Why 9KAR fails / pool best 1.40 Å, top-75 best 5.34 Å (both ORACLE, they use the native): the score discards it (`s26/results/e_trace_9KAR.json`).
33. Zero-information control / constant α-helix 4.065; random torus draw 4.522 is worse (`s13/SPRINT13_DOSSIER.md`).
34. Grid oracles / order statistics; only split-half transfer survives (`s24/LEDGER.md` L15-A).
35. What "not a result" means / |effect| below 0.7× its own MDE; say "underpowered", never "trend".

## X.5 The three directions, as the evidence leaves them

**Direction A: ADAPT-VQE.** What the evidence says. The deployed Hamiltonian is diagonal and nearly target-independent (`s25/LEDGER.md` L17); its optimum at α = 1 is a basis state; an ADAPT-grown circuit on it is a product circuit whose gradients do not decay because nothing needs to be entangled (`s26/results/q_var.json`, L35); the fixed ansatz already spans so(2ⁿ) from depth 2 (`q_dla.json`, L27); the Tang minimal pools span only so(2^(n−1)+1) (same file). The endpoint experiment A1 is pending and expected null (`s26/agentQ_FINDINGS.md` §1). So ADAPT cannot improve RMSD here, for the same reason the fixed ansatz cannot: the set-equality theorem bounds every readout by the classical prefix. What ADAPT could do on this project is scientific, not accuracy: it is the cleanest available demonstration that "large gradient" and "trainable" are different things, and that an adaptive method on a product-state target grows nothing. That is a paper-sized observation and it is in her group's language. The honest version of direction A is therefore: "ADAPT as a diagnostic on classical-shaped Hamiltonians", or, if she wants a non-trivial target, a Hamiltonian that is not diagonal, which this project does not have.

**Direction B: the AlphaFold-hardest targets.** What the evidence says. Nothing measured. ESMFold did not run (`s26/results/b1_feasibility.json`). The project's own hard set, FAIL18, is 10/18 amyloid fibril and lasso peptides (`s12/instrument.py`; `s12/SPRINT12_DOSSIER.md`), and on those the sequence-blind pipeline beats the shipped one (5.425 vs 6.019), which is a real and unexplained signal. No fresh benchmark exists (all 204 clusters spent; 16 fresh targets in the world, 10 amyloid). Direction B as "beat AlphaFold where it fails" is unmeasurable here; direction B as "characterise why the shipped pipeline fails on the fibril/lasso class and whether a fragment library without the sequence channel is the right tool for that class" is measurable on the instrument you have.

**Direction C: learn the ranking, then add physics.** What the evidence says. Learning the ranking is closed at four levels (`s17/LEDGER.md`; `s12/results/agg_dec_v2*.json`); adding physics as a ranker is worse than noise (`phys_suite.json`); adding physics as a relaxer is worse than a random move of its own size (`s26/LEDGER.md` L39). What is open is one level up: learn a better distance prior, because the prior is the only steep lever (`s24/results/priorladder.json`) and the S26 ladder is the first attempt (`s26/PROPOSAL_C.md`, rungs C2 to C5). Direction C survives only in that form: "learn the prior", not "learn the ranking", and "physics for validity", not "physics for accuracy".


## X.6 Questions a VQA postdoc will ask, with answers in your voice

Each answer is written to be spoken. Numbers carry their artefact so you can defend them.

**1. What is your Hamiltonian, exactly?**
"Diagonal. H = diag(E), where E is the rank-standardised score of 128 candidate structures, a ladder from −1.72 to +1.72. There are no off-diagonal terms. That is why the ground state is a basis state and why my theorem says the CVaR tail is always inside the classical prefix (`core/pipeline.py` `_zrank`; `s25/results/q_verify.json`)."

**2. Then why use a quantum circuit at all?**
"The founding idea was that choosing a jointly consistent subset was combinatorial. It is not, on this encoding, because the readout consumes the set mean and the tail is a prefix of the ranking. I did not know that when I started; I proved it in Sprint 24 and verified it in Sprint 25 on 2,592 cells. The honest answer is that the circuit is a way of choosing an entropy, and entropy is what the readout responds to: correlation −0.74 with RMSD across 18 readout arms (`s25/results/q_alpha.json`)."

**3. What does the circuit contribute to accuracy?**
"Minus 0.14 Å against the plain argmin, at 0.68 of my minimum detectable effect. By my own fixed rule that is not a result. Against an exact Boltzmann distribution at the same temperature it is minus 0.0002. It ships off (`q_alpha.json`; `core/pipeline.py` `Config.quantum=False`)."

**4. Is it a real CVaR-VQE?**
"The objective is CVaR at level α minus T times the entropy. The deployed table has T = 0.3 on every fold and α = 1 on three of five, so on those folds it is an entropy-regularised expectation, not a CVaR. I found that in Sprint 25 while verifying an inherited claim, and it is in the ledger (`s25/LEDGER.md` L3)."

**5. Statevector or shots?**
"Exact statevector. No sampling anywhere in the trained path; the random number generator is used for the initial angles only. Gradients are exact parameter-shift, cosine 1.000000000 against independent finite differences, relative error 4.6e-10 (`s25/results/q_verify.json`)."

**6. How do you know the parameter-shift rule applies to your objective?**
"It applies to the expectation values; CVaR and entropy are nonlinear in p, so I compute ∂F/∂p analytically, envelope theorem for the tail, and chain it through ∂p/∂θ from parameter shift. The check against finite differences on the full F is the 4.663e-10 figure (`q_verify.json`, free-energy gradient)."

**7. You mention a gradient defect. What was it?**
"The first score-function estimator centred its baseline on the tail only, so the term b·E_p[∇ log p] was not zero. Its cosine with the true gradient was −0.023 when first priced. Measured against the exact gradient today the tail-baseline estimator has cosine 0.57 on the deployed instrument and 0.68 on the audit instrument. The default is the constant baseline, cosine 1.000000 (`docs/FINDINGS.md` corrections table; `q_verify.json`; `verify/cvar_audit.json`)."

**8. Do you see a barren plateau?**
"I measured gradient variance at depth 3 for n = 4 to 13. The deployed cells decay at −0.24 to −0.31 log₂ per qubit, and the linear objective at −0.65; a 2-design would be −1. I do not call that a barren plateau and I do not call it its absence. It is a depth-3 measurement with 3n parameters against dim so(2ⁿ), and the DLA is the full so(2ⁿ) from depth 2, so nothing protects the ansatz at scale (`s25/results/q_plateau.json`; `s26/results/q_dla.json`)."

**9. How did you compute the DLA?**
"Closure on Pauli strings as a set, cross-checked by dense SVD at n = 4 and 5, 12 of 12 cells agreeing. Every generator is an odd-Y real string so the algebra sits in so(2ⁿ). At n = 7 it reaches 8128, the whole thing, at depth 2. My pre-registered guess was a proper subalgebra around 4095; it was wrong and the measured value stands (`s26/results/q_dla.json`; `s26/PREREG_A2.md`)."

**10. Why does the CVaR objective have a shallower decay than the linear one?**
"The tail is nonlinear in p; the variance ratio of the CVaR gradient to the linear gradient grows with n, 0.52 at n = 7 to 3.36 at n = 13. It is the tail's amplification, not a property of the circuit (`q_plateau.json`, `nonlinearity_variance_ratio`)."

**11. Have you tried ADAPT-VQE?**
"Sprint 26 did, on this Hamiltonian. At α = 1 the grown circuit selects one to three distinct operators and the rest are repeats of the same single-qubit Y; it is a product circuit and its variance does not decay with n. That is not trainability; it is the trivial regime, because the target state is a basis state. At α = 0.25 the 2-local pool grows 4 to 21 strings and decays at −0.30 per qubit, the same as the fixed ansatz. The RMSD endpoint is pending and I expect it to be null, for the theorem reason (`s26/results/q_var.json`; `s26/PREREG_A1.md`)."

**12. Which pool?**
"Two: Tang's minimal 2n−2 pool, which I measured generates so(2^(n−1)+1), about a quarter of so(2ⁿ), and the full set of 1- and 2-local odd-Y strings, which generates so(2ⁿ) (`q_dla.json`, `pools`). I did not implement CEOs; I read your paper and I would like to know whether a fermionic-style pool means anything on a diagonal classical Hamiltonian."

**13. What is the entanglement doing?**
"Nothing measurable. An entangler ablation on the continuous encoding is −0.013 [−0.095, +0.077] and the mutual information across the middle cut is 0.045 bits (`s20/LEDGER.md`). The MPS twin reproduces the dense circuit at χ = 2^layers exactly (`q_verify.json`)."

**14. What would make the quantum part matter?**
"A Hamiltonian that is not diagonal, so that the ground state is not a basis state and the set-equality theorem does not apply. I do not have one. If you know a physically motivated one for a discrete conformer choice, that is the conversation I want."

**15. What is your accuracy, and against what?**
"3.21 Å mean Cα RMSD on 126 cluster-disjoint held-out peptides, median 2.97, 51% under 3 Å. Against a random 75-subset of the same pool at 3.43 on the point-cloud basis, a constant helix at 4.07, and a sequence-only torsion predictor at 3.77. Against a perfect chooser from the same pool, 1.71, which is an oracle and not achievable (`bench_results/compare_tuning126.json`; `results/summary/leaderboard.csv`; `s25/results/phys_suite.json`; `s13/SPRINT13_DOSSIER.md`)."

**16. Why not compare to AlphaFold?**
"ESMFold does not run on my machine: no checkpoint, no openfold, 8.76 GB resident against 4.4 GB headroom. I have not measured it and I will not quote a number I have not measured (`s26/results/b1_feasibility.json`)."

**17. Did the improvement replicate on held-out data?**
"No. The tuning gain of −0.25 Å did not transfer to the sealed 60-target benchmark: +0.010 [−0.16, +0.18], 31 wins to 29. Ten targets carried 61% of the tuning gain. I spent the benchmark once and it cannot be spent again (`s9/final_report.json` via `verify/headline_audit.json`)."

**18. So what is the result of the project?**
"That the distance prior is the only steep lever, −2.15 Å per unit toward truth, and that every downstream lever, including the quantum one, is flat to hundredths; and an instrument that can tell those apart. Plus the trainability chain (`s24/results/priorladder.json`; Part VII)."

**19. Why do the physics force fields not rank the native?**
"On pools of real fragments, Legacy is a compactness measure and AMBER puts the native at the 51st percentile. In the seven-configuration suite both are worse than a random subset by 0.33 and 0.46 Å, and a rank-permuted channel reproduces that. AMBER's role is validity: it turns a CA trace into an all-atom model. Its relaxation moves the trace away from the native, cosine −0.05 (`s25/results/phys_suite.json`; `s20/LEDGER.md`; `s26/LEDGER.md` L39)."

**20. What is the minimum detectable effect?**
"2.8016 times the standard error of the paired difference, per comparison. An earlier fixed constant of 0.084 Å was wrong by up to 84× across comparisons. Anything below 0.7 of its own MDE I call underpowered and do not report as a result (`s21/LEDGER.md` L8)."

**21. Why fold-clustered intervals?**
"Each fold's distance model was trained with about 100 of the other 125 natives in its training set, so per-target errors are correlated through the models. The i.i.d. bootstrap is anticonservative; I print both (`s20/LEDGER.md` L-B)."

**22. How do you know there is no leakage?**
"The native is read by one function, `label()`. The tests poison it with NaN and check every emitted quantity is bit-identical; the audit is clean on 19 of 19 quantities. The one declared defect is the identity normalisation used for clustering; its measured price is +0.0004 Å on tuning and +0.0030 on the benchmark (`verify/leak_audit.json`; `core/data.py`; `docs/FINDINGS.md` S10-4)."

**23. What is the biggest mistake you made?**
"Quoting a +0.113 Å CVaR contribution for two sprints. It was measured at a temperature that does not ship, and it was a difference of marginal means where a paired statistic was required; paired, it was 0.51× MDE with median zero. I withdrew it in the ledger (`s25/LEDGER.md` L3, L5)."

**24. And the second biggest?**
"Reporting the point-cloud average, 3.05 Å, as the structure for many sprints. It is contracted 22% and is not a chain. The built chain is 3.21. The code labels the average 'illegal'; I should have read that sooner (`core/bench.py` line 682; `s25/LEDGER.md` L4, L8)."

**25. Why did you pre-register?**
"Because three times in two sprints a control priced a different quantity than its name, and every bias I had not stated pointed my way. From Sprint 16 every lane wrote its falsifier before running; in Sprint 25 my registered prior on the physics suite was 'every AMBER configuration fails', and it held (`s25/PREREG_PHYS.md` §4)."

**26. What would you do with an independent study?**
"Measure the one open lever, the prior, on the instrument I have; and write up the trainability chain with the DLA and ADAPT results as a paper about what a Pauli spectrum does and does not predict. Part X.8 has the order I would do them in."

**27. Is 7 qubits meaningful?**
"No, and I do not claim it is. n = 7 is why it is simulable. The width sweep to n = 13 is the only scaling statement, and it is a depth-3 statement (`q_plateau.json`; `s26/LEDGER.md` L27)."

**28. What is the readout?**
"The trained distribution either picks a p-weighted medoid or weights the coordinate average. Either way it is insensitive to distributional differences of nearly half the mass: the trained and Gibbs states differ by total variation 0.45 and the endpoints differ by 0.03 Å, a quarter of the MDE (`s25/LEDGER.md` L15; `q_gibbs.json`; `q_alpha.json`)."

**29. Does the circuit actually optimise?**
"At T = 0.3 it closes 78% of the free-energy gap from the initial mean to the Gibbs optimum and beats the best of 200 random initialisations by 0.70; it is still 0.90 nats from Gibbs. At T = 0.1 the trained state is worse than the point mass at the argmin, −1.7176 against −1.7186 (`q_gibbs.json`)."

**30. What is the temperature doing?**
"T·H is the entropy bonus. Free energy minus its optimum equals T times the KL divergence to the Gibbs state, so the training curve is a divergence ladder. At the deployed T the trained state has 5.67 bits against the Gibbs state's 4.91 (`q_gibbs.json`)."

**31. Why a rank ladder and not the raw score?**
"So that every channel lands on the same marginal and the effective temperature does not differ between configurations; the outer rank-standardisation is worth up to 0.14 Å of silent temperature shift otherwise. The cost is that the Hamiltonian is then the same for every target to 1.2%, which is why the circuit trains toward the same state every time (`s25/LEDGER.md` L16, L17)."

**32. If I gave you a real quantum device tomorrow, what would you run?**
"Nothing from this pipeline, because the exact statevector at n = 7 already answers every question the pipeline asks. I would want a non-diagonal Hamiltonian first."

**33. What do you not know?**
"Where the target-specific third of the coherent error comes from; whether any native-free proxy reaches the 0.638 in-band correlation that 2 Å needs; whether the S26 prior ladder moves the endpoint; and why n = 6 and n = 9 have a proper DLA subalgebra (Part VIII.2)."

If a question is outside this list: "I do not know. The artefact that would answer it is in the repository and I can check and send you the path."

## X.7 Must-not-claim, with the true replacement

| do not say | say instead | artefact |
|---|---|---|
| "The CVaR-VQE improves accuracy by 0.11 Å" | "The CVaR-VQE's contribution is −0.14 Å at 0.68× MDE, not a result; the 0.113 was withdrawn" | `s25/results/q_alpha.json`; `s25/LEDGER.md` L3/L5 |
| "3.05 Å mean" | "3.21 Å built chain; 3.05 is the point-cloud intermediate" | `bench_results/compare_tuning126.json`; `s25/LEDGER.md` L8 |
| "sub-2 Å" | "28.6% of targets under 2 Å; mean 3.21" | `results/summary/leaderboard.csv` |
| "no barren plateau" | "depth-3 slopes of −0.24 to −0.31 log₂ per qubit; DLA is full so(2ⁿ), so no protection at scale" | `s25/results/q_plateau.json`; `s26/results/q_dla.json` |
| "the circuit trains, so ADAPT would too" | "on this Hamiltonian ADAPT grows a product circuit at α = 1; large gradient is not trainability" | `s26/results/q_var.json` |
| "quantum advantage" or "quantum speed-up" | nothing; the register is 7 qubits and exact | `q_verify.json` |
| "AMBER refines the structure" | "AMBER relaxation is a validity step; its displacement is worse than a random move of its own size" | `s26/LEDGER.md` L39 |
| "physics helps ranking" | "both energies are measurably worse than noise as rankers" | `s25/results/phys_suite.json` |
| "the improvement is significant" (about the benchmark) | "+0.010 [−0.16, +0.18] on the sealed 60; did not transfer" | `s9/final_report.json` via `verify/headline_audit.json` |
| "MDE 0.084 Å" | "MDE is 2.8016 × SE per comparison" | `s21/LEDGER.md` L8 |
| "the rank correlation is +0.37" (or +0.355) | "+0.161 under a uniform proposal; the +0.355 was under a proposal drawn from the prior" | `s14/LEDGER.md` C7-CORRECTION |
| "entanglement matters" | "entangler ablation −0.013 [−0.095, +0.077]; MI 0.045 bits" | `s20/LEDGER.md` |
| "we beat AlphaFold on hard targets" | "not measured; ESMFold did not run here" | `s26/results/b1_feasibility.json` |
| "consensus finds the native" | "consensus is outlier avoidance; the native sits at the 82.8th percentile of the medoid criterion" | memory `consensus-is-outlier-avoidance` |
| "more search would help" | "the exact argmin over the whole latent ties the pool; search is closed" | `s21/LEDGER.md` |
| "a 36-qubit model" | "7 qubits; 9 in the suite" | `core/pipeline.py`; `s25/LEDGER.md` L16 |
| "the pipeline is state of the art" | "3.21 Å on 126 peptides; a fragment-averaging baseline with an exact selector and a measured ceiling" | Part I.5 |

## X.8 Which direction the evidence favours

The evidence favours a version of direction A that is scientific rather than accuracy-seeking, paired with the one measurable accuracy question, which is direction C in its "learn the prior" form.

Reasoning. Directions that promise accuracy through selection, ranking, physics, or search are closed by measurement (Part VIII.1). The only steep lever is the prior (`s24/results/priorladder.json`), and the S26 ladder is already testing whether a better prior is learnable on this machine (`s26/PROPOSAL_C.md`). That is the accuracy study, and it needs no quantum computer. The quantum study that the evidence supports is the one her group is equipped to judge: the measured chain from force-field locality to gradient variance, the DLA census with its two unexplained widths, and the ADAPT product-circuit observation (`s26/PROPOSAL_B_REPLACEMENT.md`). Its natural next question, "what does an adaptive ansatz do when the target Hamiltonian is not diagonal, and can the Pauli-spectrum-times-kernel prediction be extended to that case", is a question this project cannot answer alone and her group can.

So the ask on slide 11 is: an independent study with two threads. Thread one, the prior ladder to completion on the 126-target instrument, pre-registered, reported with fold-clustered intervals. Thread two, the trainability paper, with her group's input on a non-diagonal test Hamiltonian and on whether CEO-style pools change the product-circuit result. Direction B is not proposed as stated; its measurable remnant (the fibril/lasso class where the sequence-blind pipeline wins) can be a sub-question of thread one.

## X.9 The night before: ten items

1. Read X.1 aloud twice with a timer. Ninety seconds is the ceiling.
2. Re-read Appendix D. Know the five things that do not reproduce or exist only in prose, so that none of them is on a slide.
3. Open `results/summary/leaderboard.csv` and read the production row once: 3.2126, 2.9661, 0.18242, 8.234, 0.2857, 0.5079. Open `bench_results/compare_tuning126.json` and read `science.rmsd_arm` 3.2148. Know which is which.
4. Open `s25/results/q_verify.json` and read the ten verification lines (Part V.2 to V.4). These are the numbers she will test you on.
5. Open `s25/results/q_plateau.json` and say the five slopes from memory; then say the sentence in V.9 about scope.
6. Open `s26/results/q_dla.json` and `q_var.json` and say L27 and L35 in one sentence each.
7. Say the set-equality theorem in one sentence and its two numbers (0 violations; 972 of 972).
8. Say X.2 once, plainly, and decide where in the ten minutes it goes (early).
9. Write the ask on a card in your own words: two threads, one measurement each.
10. Sleep. A tired presenter quotes marginals.


---

# APPENDIX A. GLOSSARY

Protein terms.

- **Amino acid**: one of twenty small molecules that chain together to form proteins.
- **Residue**: one amino acid inside a chain.
- **Sequence**: the ordered list of residue types, written as one capital letter per residue.
- **Peptide**: a short protein; here 9 to 16 residues.
- **Backbone**: the repeating N, CA, C atoms of each residue; **CA** (Cα, alpha carbon) is the atom this project predicts; adjacent CA atoms are about 3.8 Å apart.
- **Side chain**: the residue-specific atoms attached to CA.
- **Native structure**: the shape a sequence folds into in nature; the ground truth here, taken from the PDB.
- **PDB**: the Protein Data Bank; four-character codes such as 1S9Z name deposited structures.
- **NMR**: nuclear magnetic resonance spectroscopy; the experimental method behind most peptide structures; yields a bundle of models, of which the first is used.
- **Alpha helix, beta sheet, coil**: the three named local shapes (spiral; side-by-side strands; everything else).
- **Torsion angles φ, ψ, ω**: the rotation angles about the three backbone bonds per residue; φ and ψ are free, ω is 180° (trans) except in rare cis bonds.
- **Ramachandran plot**: φ against ψ; real residues cluster in a few regions; the projection penalty `ramah` pushes toward them.
- **Distance map / contact**: all pairwise CA distances; a contact is a pair closer than a cutoff.
- **Distogram**: a predicted probability distribution over each pairwise distance, here in 17 bins.
- **Fragment / window**: a contiguous stretch of a real structure, cut to the target's length; the candidates come from these.
- **BLOSUM62**: a standard substitution-score table for aligning sequences; used to rank library windows.
- **ESM-2**: a protein language model whose embeddings are used as features.
- **AlphaFold2, ESMFold**: large structure predictors; neither was run in this project.
- **Force field**: a function from atom positions to energy; **AMBER ff14SB** is one; **GBn2** is its implicit water model; **OpenMM** evaluates them.
- **Legacy**: this project's own eleven-term CA-level energy.
- **Radius of gyration**: a measure of how compact a structure is.
- **Amyloid fibril, lasso peptide**: structural classes over-represented among the project's failures.
- **Chignolin, trpzip, Trp-cage**: three famous small model peptides; the first email named the first two; none is on the instrument.
- **Chemical shift**: an NMR observable per residue that constrains torsions; available for 54 of 126 targets.

Measurement terms.

- **RMSD**: root mean square deviation between corresponding CA atoms after optimal rigid superposition (Kabsch); in Å.
- **Kabsch**: the algorithm that finds the superposition minimising RMSD.
- **Point cloud**: a CA set that need not obey 3.8 Å spacing; the average of structures is one.
- **Built chain**: a CA set with ideal spacing and torsions; the reporting basis.
- **Medoid**: the member of a set with smallest total distance to the others.
- **Pool**: the K = 500 retrieved candidates. **Top-75 / top-m**: the kept subset.
- **ORACLE**: any arm that reads the native to make a choice; a diagnostic, never a result.
- **Control**: an arm that removes the information under test while keeping the operator (random 75-subset; permuted channel; constant helix).
- **Paired**: arm A minus arm B on the same target, then averaged.
- **SE**: standard error of the paired mean.
- **MDE**: minimum detectable effect, 2.8016 × SE; 80% power at 5% two-sided.
- **NOT MEASURED / UNDERPOWERED / NOT A RESULT**: |effect| ≤ MDE; below 0.7× MDE.
- **Fold-clustered CI**: a bootstrap that resamples folds, because per-target errors are correlated through the fold models.
- **W/L**: wins and losses per target; printed, never a verdict.
- **Order statistic**: the best of K correlated draws looks like a signal and is not.
- **Common-mode error**: the part of the candidates' error that is shared across them; 68% here.
- **Concentration**: a mean gain carried by a few targets; tested against a uniform-effect null.
- **Type-M**: magnitude inflation of a low-power significant effect.
- **Leave-fold-out**: train on four folds, evaluate on the fifth.
- **tuning126, dev24, benchmark60, FAIL18**: the target sets (Part II.1).
- **Pre-registration / falsifier**: writing before the run what result would count as failure.

Quantum terms.

- **Qubit, basis state, register**: n qubits index 2ⁿ basis states; here a candidate index.
- **Ansatz**: a parameterised circuit; here RY rotations and CNOT chain plus ring, 3 layers, 21 parameters.
- **Statevector**: the exact 2ⁿ-amplitude description of the state; used here directly.
- **MPS / bond dimension χ**: a tensor-network form of the state; χ = 2^layers here, exactly.
- **Expectation, diagonal Hamiltonian**: H = diag(E); ⟨H⟩ = Σ p_i E_i.
- **CVaR_α**: the mean of the lowest-α probability mass of the energy distribution.
- **Free energy F = CVaR − T·H**: the deployed objective; H is Shannon entropy, T a temperature.
- **Gibbs / Boltzmann distribution**: p* ∝ exp(−E/T); the unconstrained minimiser of ⟨E⟩ − T·H.
- **KL divergence, total variation**: two distances between distributions; F − F* = T·KL.
- **Parameter-shift rule**: the exact gradient formula for rotation gates.
- **Score-function (REINFORCE) gradient, baseline**: a sampled gradient estimator and the constant subtracted to reduce its variance; the tail-only baseline was a defect.
- **Adam**: the optimiser; 50 iterations, learning rate 0.15.
- **Barren plateau**: gradient variance decaying exponentially in n; not claimed here in either direction.
- **2-design**: a circuit whose random outputs match the Haar measure to second moments; gives Var ~ 1/dim.
- **DLA (dynamical Lie algebra)**: the Lie algebra generated by the circuit's gate generators; so(2ⁿ) here from depth 2.
- **Pauli spectrum / Walsh basis / Pauli weight**: the expansion of a diagonal energy in products of Z operators; the weight is how many qubits a term touches.
- **Locality**: which residues a distance term depends on; d_ij depends on exactly the j − i − 1 between.
- **ADAPT-VQE**: an ansatz grown one operator at a time from a pool by largest gradient; **qubit-ADAPT**: the Pauli-string pool version; **CEO pool**: Ramôa's coupled-exchange-operator pool; **gradient trough**: ADAPT's small-gradient stall before convergence.
- **Set-equality theorem**: the CVaR tail support is a subset of the classical energy-order prefix.
- **QNG**: quantum natural gradient; null here.
- **VQE_LFO**: the per-fold (α, T) table chosen leave-fold-out.
- **Readout**: how the trained distribution becomes a structure (weighted medoid or weighted average).

---

# APPENDIX B. EVERY NUMBER, WITH ITS PATH

Grouped by where the number first appears. A number that appears more than once is listed once. "state brief" is `docs/STATE_BRIEF_2026-09-12.md`. "external" marks the four items about the audience that were read from the web, not from the repository.

## B.1 Headline and results lab

| number | meaning | path |
|---|---|---|
| 3.214765154210998 | built-chain mean, 126 targets, production cache | `bench_results/compare_tuning126.json`, `science.rmsd_arm` |
| 3.048338093879532 | point-cloud mean | same, `science.rmsd_avg` |
| 3.2354598538973844 | after AMBER relaxation | same, `science.rmsd_full` |
| 3.20407616038092 | projection at lam = 0, no torsion penalty | same, `science.rmsd_fit` |
| 3.4540004952559396 | shipped argmin selection (S8 instrument) | same, `science.shipped` |
| 1.7108244199364904 | ORACLE pool best (uses the native) | same, `science.pool_best` |
| 2.3061526409453816 | ORACLE top-m best (uses the native) | same, `science.top_m_best` |
| 0, ±4.4e-16 | science delta 1 vs 8 workers | same, `science_delta` |
| 8.371; 2537.568 s; 303.137 s | speed-up; wall times | same |
| 3.2126, 2.9661, 0.18242, 8.234 | leaderboard production mean, median, best, worst | `results/summary/leaderboard.csv`, production row |
| 0.2857, 0.5079 | fraction under 2 Å, under 3 Å | same |
| 0.7150 | corr with pool best | same |
| 0.1643 | basis delta (built minus point cloud, rebuild) | same |
| 3.2187/3.0580 (+0.0061, 0.19×) | distogram configuration, rebuild | same |
| 3.3100/3.1317; 3.3732/3.2151; 3.4221/3.2530; 3.8248/3.6742; 3.8844/3.7553; 4.1015/3.8805 | the other six configurations, rebuild | same |
| 0.18242222359245708; 0.19579; 0.15960; 0.17656 | T030 built chain; point cloud; pool best; distogram config | `results/summary/results.csv` line 913 and T030 distogram row |
| 0.39949916705333455; 0.4167616469202071 | T053 2MJQ; T098 6WPB | `results/summary/results.csv` production rows |
| 8.2342; 7.4336; 7.03; 6.9715 | worst four: 2MQ2, 9KAR, 2BFI, 7JS6 | same |
| 126 | targets in target map | `results/summary/target_map.json` |
| a15406c82245, git_dirty True | provenance of the results lab | `results/summary/leaderboard.json`, `provenance` |
| 2142 structures, 6.6 MB | site build | `s25/results/frozen_build.log` |
| 0.18198112330908295; 7.437696305468429 | T030 and T122 built chain, production cache | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json`, `reporting.rmsd_arm` |
| 3.2126 vs 3.2148 | two routes | Appendix D, D2 |

## B.2 Instrument

| number | meaning | path |
|---|---|---|
| 126; 5 folds; 60; 24 | target sets | `s12/instrument.py`; `core/data.py` |
| 1fc9f2dcf489e2fb | production config key | `s12/instrument.py`, `PROD_KEY` |
| a40581ad…8422d, 8002 bytes | benchmark manifest hash | `core/data.py` |
| 9 to 16 | target lengths | `core/data.py`, `benchmark()` |
| 18; 10/18; p = 1.2e-6 | FAIL18 size; fibril/lasso; Fisher | `s12/instrument.py`; `s12/SPRINT12_DOSSIER.md` |
| 2.9610; 2.9507; +0.0103 [−0.1596, +0.1803]; 31W/29L; SE 0.0867 | benchmark pass | `s9/final_report.json` via `verify/headline_audit.json` |
| −0.250 [−0.383, −0.117]; −0.255; +0.2013; +0.3296 | tuning gain; dev gain; drop-10; drop-20 | `docs/FINDINGS.md` S9-10 |
| 3.204; 3.221 | preregistered tuning and dev means | `verify/headline_audit.json`, `benchmark60.preregistration` |
| 2.8016; 0.084; 0.11× to 9.54×; 84 | MDE constant; old constant; error range | `s21/LEDGER.md` L8 |
| 0.7× | underpowered rule | `s25/results/q_verify.json` |
| ~100 of 125 | natives seen by each fold model | `s20/LEDGER.md` L-B |
| 13/126; 4 at 1.0; +0.0004; 13/60; 8Y3S, 8ZG3; +0.0030 | identity leak counts and prices | `docs/FINDINGS.md` S10-4, S9-10; `s26/LEDGER.md` L31 (git `5fa05cd`) |
| 19 of 19; 3 targets | leak audit | `verify/leak_audit.json` |
| 126/126; 0.0; 575.7 s | projection exactness | `verify/project_exactness.json` |
| 370; 357; 13; 11 + 2 | tests | `s26/results/test_run.json`; `s26/TEST_RUN.md` |
| 1.692 GB; 240 s; 261 s; 256 s | test job peaks and times | `s26/TEST_RUN.md` (quoted in state brief §6.1) |
| 4.624; 3.425; 4.065; 3.770; 2.226; 1.711; 1.313; 1.594; 0.611 | ladder: the first four are controls, the last five are ORACLE (use the native) | state brief §2; `s25/results/phys_suite.json`; `s13/SPRINT13_DOSSIER.md`; `s24/results/priorladder.json`; `s15/LEDGER.md` |
| 4.522; −0.457 | random torus draw; helix minus random | `s13/SPRINT13_DOSSIER.md` |

## B.3 Pipeline

| number | meaning | path |
|---|---|---|
| k=500, m=75, lam=0.3, amber_k=10, amber_steps=0, min_sep=2, maxiter=300, vqe 7/3/50/seed 0, legacy_top=5 | Config defaults | `core/pipeline.py`, `Config` |
| 787; 6,003; 1,001 | peptides; fragments; deposits | `README.md`; state brief §4 |
| 7,193; 9,814 | windows for 16-mer and 15-mer | `s26/results/e_trace_1S9Z.json`, `e_trace_9KAR.json` |
| −1 to 19; 144; −6 to 32; 160 | BLOSUM ranges; peptide-origin windows | same |
| 17; edges and centres as listed; 4.0 Å last gap | distogram bins | `core/predict.py` |
| 2.0 to 40.0 step 0.05 | risk grid | `core/predict.py` |
| 1.9962; 0.6159; 0.2824; 0.2412 | z sd; cov90; cov50; multimodal | `s25/results/calib.json` (means over rows) |
| 5.05 to 21.07; 0.23 to 2.70; 5.62 to 23.44; 0.23 to 7.58 | expected distances and sd for the two traces | trace files, `stages.distogram` |
| 0.693 to 7.269; 488; 2.242 to 5.179; 483 | score ranges and distinct counts | trace files, `stages.score` |
| +0.153 [+0.066, +0.247]; 50W/76L; 3.610; 3.763 | flat weights vs sd weights in refinement | `s18/results/q_ceilattack2.json`; `s18/LEDGER.md` L10 |
| 1.16; 0.04; 0.893; 0.803/0.298 | operator law | `s12/SPRINT12_DOSSIER.md` |
| −0.239 | ORACLE per-target m (uses the native) | `s22/LEDGER.md` |
| 1.3997; 5.3436; 0.1596 | 9KAR pool best, top-75 best; 1S9Z pool best (all ORACLE, use the native) | trace files, `reporting` |
| −1.7186 to +1.7186; 3.4372; 0.394%; 1.18%; 0.0406; 3.4371 | rank ladder; deviations | `e_trace_1S9Z.json` `stages.hamiltonian`; `s25/results/q_gibbs.json` |
| 449; 44; 34 | selected pool index; medoid local indices | trace files |
| 3.75/3.36; 1.97/1.16; 3.80 | point-cloud bond mean/min; native | trace files, `stages.average` |
| 22.3%; 2.961; 3.812; 0.649 | contraction; mean virtual bond; native; worst | state brief §2 |
| 3.804; 3.8039549 | ideal bond; achieved | `core/project.py`; trace files, `stages.projection` |
| 0.042; 1.90; 126/126 | analytic-mode movement | `core/project.py` docstring |
| +0.1664; 3.26× | projection cost | `s25/LEDGER.md` L8, L11 |
| 5,370 → −1,291; 0.191; 6.0e12 → 1,262; 0.610 | AMBER energies and movement, two traces | trace files, `stages.amber` |
| 92.0; 1000.0; N, CA, C | AMBER memory limit; convergence cap; restrained atoms | `core/amber.py` lines 249, 1282, 808 |
| +0.0111 [+0.0062, +0.0171]; +0.0385; cos −0.049; 0.220; 58.7%; 125/126 | S26 relaxation findings | `s26/LEDGER.md` L39, L24 |
| 13.13 s; 13.21 s | trace wall times | trace files, `wall_s` |
| 0.676; 50.7× | common-mode | `s23/LEDGER.md` L9 |

## B.4 Physics

| number | meaning | path |
|---|---|---|
| 4.0, 1.0, 1.0, 3.0, 2.0, 2.0, 0.5, 1.0, 0.8, 0.15, 0.4 | Legacy weights | `core/energy.py` line 244 |
| −0.4477 [−0.5110, −0.3911]; 124W/2L; 200 partitions | Legacy is compactness | `s20/LEDGER.md` |
| 0.139; 0.885 | certified optima vs random | `s14/LEDGER.md` |
| 6/35 → 35/35 | side-chain coverage | memory `sidechains-all-20-residues` |
| 58/97/99% | Relax50 variance explained | `s20/LEDGER.md` |
| 63,000 | AMBER single points cached | `s24/cache_amber/`; `s26/LEDGER.md` L1 |
| 51st percentile | AMBER native rank | memory `physics-ranks-real-geometry-not-lattice`; `s20/LEDGER.md` |
| 1.00/1.41/1.73; 0.14 Å | channel-sum sd; temperature shift | `s25/LEDGER.md` L16 |
| 0.1127, 0.7529, 0.8013 | mean abs z_moment per channel (document-only) | `s25/LEDGER.md` L16; Appendix D, D10 |
| 0.18; 0.5; 3; 80; 9 qubits; 12-target probe; 75 | suite settings | `s25/results/phys_suite.json`; `s25/LEDGER.md` L16 |
| 3.0580; 3.1317; 3.2151; 3.2530; 3.4251; 3.6742; 3.7553; 3.8805 | suite means and random null | `s25/results/phys_suite.json`, `configs_rank`, `random_null_rank` |
| +0.0097 0.43×; +0.0833 0.70×; +0.1668 0.91×; +0.2047 1.07×; +0.6259 1.89×; +0.7070 2.14×; +0.8322 2.67× | vs incumbent | same, `vs_incumbent` |
| +0.3302 1.99× 41W/85L; +0.4554 2.42× 30W/96L | Legacy, AMBER vs random | same, `controls_rank.1/2.vs_random` |
| +0.3274; +0.4712 26W/100L | permuted channels | same, `controls_rank.1/2.vs_perm` |
| +0.0331 0.43× | AMBER inside distogram-led | same, `controls_rank.5.vs_perm` |
| 0.0000 | VQE vs classical top-m, all seven | same, `controls_rank.*.vs_topm` |
| 2.9400; −0.1181; 0.9536; 1.0238 | ORACLE min(C3, C5), uses the native; k_eff | `s25/results/phys_form5.json`, `oracle_ceiling` |
| +0.0049 0.14× 97 ties; +0.0053; +0.0737 0.58× | Form 5 primary, null, always-C5 | same |
| 0.0148; w = −0.5 | fitted AMBER weight ceiling | `s24/LEDGER.md` L16 |
| ρ(AMBER, distogram) −0.0270; ρ(Legacy, distogram) +0.3875 | registered prior inputs | `s25/PREREG_PHYS.md` §4 |

## B.5 Quantum

| number | meaning | path |
|---|---|---|
| 7; 3; 21; 128; 512 | qubits, layers, parameters, dims | `s25/results/q_verify.json`; `s25/LEDGER.md` L16 |
| 80 qubits; 0/20,000; 20 qubits | one-hot vs binary encoding | `docs/FINDINGS.md` §6 |
| 5.551e-17; 1.000000000; 4.597e-10; 4.663e-10 | verification | `s25/results/q_verify.json` |
| 3.331e-16; 6.661e-16; χ 2/4/8/16; Schmidt rank 4; 0.3339 | MPS checks | same |
| 0.566586; 0.519; 1.000000; −0.023; 0.6847; 0.6305; 0.685; 0.524 (D10) | baseline cosines | `q_verify.json`; `docs/FINDINGS.md`; `verify/cvar_audit.json` |
| 2592; 1620; 0; 0; 1424; 972/972 | set equality | `q_verify.json` |
| VQE_LFO table; T = 0.3 all; α = 1 on 3/5 | deployed cells | `core/pipeline.py`; `s25/LEDGER.md` L3 |
| 0.15; 0.9; 0.999; N(0, 0.6²); 1 restart | Adam settings; init | `core/quantum.py` |
| 3.3135; 3.4540; −0.1405; SE 0.0732; MDE 0.2051; 0.68×; 66W/48L/12 ties | VQE_LFO vs argmin | `s25/results/q_alpha.json`; `verify/vqe_lfo_audit.json` |
| −0.0002; +0.0262; −0.0308 | vs Boltzmann T=0.3; uniform top-64; top-128 | `q_alpha.json`, `vs_no_circuit` |
| 3.2825; −0.1715 | best cell (0.25, 0.3) | `verify/vqe_lfo_audit.json` |
| 3.3414; 3.4540; 0.1126; 0.0761 bits; 6.3638 bits | inherited claim | `q_alpha.json`, `inherited_claim` |
| −0.1126 0.51× median 0; +0.0279; +0.0126; −0.0011 | α contrasts by T | `q_alpha.json`, `primary_alpha_effect_by_T`; `s25/LEDGER.md` L5 |
| 3.2835; 3.2825; 3.3115 | T = 0.3 marginals | `q_alpha.json`; `s25/LEDGER.md` L3 |
| −0.7423; +0.2700; −0.0234; 0.704; 7.436; 18; +0.0090; 0.0268; 0.032 | entropy curve | `q_alpha.json`, `cell_correlations`, `entropy_fit`, `circuit_off_curve`, `alpha_beyond_entropy` |
| 0.619 | share with no tail constraint | `q_alpha.json` |
| −0.0302 0.24×; +0.0499; +0.0445 | circuit vs Boltzmann by T | `q_alpha.json`, `expressivity_vs_boltzmann` |
| Gibbs ladder rows (F, KL, TV, H) at T = 0.1, 0.3, 1.0 | Part V.8 table | `s25/results/q_gibbs.json`, `divergence_ladder` |
| 0.8915; 0.783; 0.7833; 0.7288; 0.6996; 0.7224; N 200 | training control | `q_gibbs.json`, `training_control` |
| 0.1 nats; fired True | falsifier | `q_gibbs.json`, `falsifier` |
| −0.6492; −0.2522; −0.0472; −0.3105; −0.2429 | width slopes | `s25/results/q_plateau.json` |
| 250/200/120/80 | θ draws by width | same, `n_theta` |
| 2.381e-2; 1.66e-2; 9.197e-3; 7.999e-3; 7.6e-3; 7.4e-3; 8.213e-3 | depth sweep n=7 | same, `depth_sweep_n7` |
| 0.5159 (n=7); 3.3585 (n=13); 0.2712 (n=4) | nonlinearity ratio | same |
| 0.504 | S13 decay base at depth 8 | `s13/results/geo_kernel.json` via `s26/LEDGER.md` L27 |
| 8128; n; 510/1023/2016; 32766/65535/130816; 36…32896; 12/12 | DLA | `s26/results/q_dla.json`; `s26/LEDGER.md` L27 |
| −0.079; +0.006; +0.035; −0.008; −0.246; −0.302; 1.63; 2.5 to 370; 31/32; 1 to 3; 4 to 21 | ADAPT variance | `s26/results/q_var.json`; `s26/LEDGER.md` L35 |
| 1.0000; 2.236; 3.015; 79/79; 0.778/0.629; 99.6%; 6.001; 0.4%; g_ii = 0.25 | S13 quantum architecture | `s13/results/qarch_locality_geom.json`, `geo_pauli.json`, `walsh_xval.json`, `walsh_predict.json`; `s13/SPRINT13_DOSSIER.md` |
| −0.210 [−0.339, −0.082]; −0.299 [−0.407, −0.195]; −0.013 [−0.095, +0.077]; 0.045 bits; +0.012/+0.013/+0.025 | S20 continuous encoding | `s20/LEDGER.md` D4; `s20/agentD_FINDINGS.md` |
| 65,536; 64 ms; 75/126; +0.025; +0.085; +1.7202 0W/119L; −1.1766; −0.697; +0.056; −1.299; 0.37; 0.200 | S21 latent | `s21/LEDGER.md` |
| −0.031 49W/49L; +0.030; +0.021; +0.017 [−0.019, +0.053] | S11 components | `docs/FINDINGS.md` S11-1 |
| 1.62 Å | reference self-disagreement | `docs/FINDINGS.md` S11-2 |
| 0/12; +0.65 to +1.32 | lattice-era VQE vs best-of-N | memory `concentration-is-wrong-when-discrimination-binds` |
| 0.00 to 0.12 vs 0.71; 10/10; 1,024 vs 8,192 | VQE vs greedy | `s14/LEDGER.md`; `s16/LEDGER.md`; `s17/LEDGER.md` |

## B.6 Sprint story and established results

| number | meaning | path |
|---|---|---|
| 12% | best-window identity | `docs/FINDINGS.md` S5-9 |
| +0.376; +0.283; 28% | calibration slope; correlation; alignment | `docs/FINDINGS.md` S7-3; `s15/LEDGER.md` |
| 3/126; 36.8th | native as argmin | `docs/FINDINGS.md` S7-5 |
| −0.288 [−0.484, −0.092]; p 0.0046; 74W/45L | ESM vs one-hot (D9) | `docs/FINDINGS.md` S7-11 |
| 2.406 | pool construction cap | `docs/FINDINGS.md` S8-5 |
| 1.386 | tie-break winner | `docs/FINDINGS.md` S8-8 |
| −0.016 [−0.037, −0.001]; −0.070; 90/126; 0.004; 0.960 | S8-6 | `docs/FINDINGS.md` |
| 0.915/0.346 r 0.934; 1.009/0.011 r 0.981; 1.019/−0.001; +0.049 [−0.009, +0.103]; 2.146 | S8-10/S8-13 | `docs/FINDINGS.md` |
| 0.044; +0.070 [−0.118, +0.259] | S8-14 | `docs/FINDINGS.md` |
| 0.492/1.634 r 0.961 | S9-4 loop | `docs/FINDINGS.md` |
| 0.2513; 0.1973; 0.0558 | positive-phi | `docs/FINDINGS.md` S9-9 |
| 2.7× | projection-matrix requirement | `docs/FINDINGS.md` S10-2 |
| 1.36; 0.90; 1.16; 0.20; 2.35; 1.49; 2.11; 0.24; 0.035; 0.064; 0.853; 1.840; 1.802; 0.953 | S10-5 headroom and hulls (all ORACLE diagnostics, use the native) | `docs/FINDINGS.md` bottom line and S10-5 |
| 2.087/2.160/1.893 → 1.987/2.044/1.802; +0.092/+0.193/+0.296; 2.7e-13; 1.604/1.124/0.828 → 1.802/1.336/0.953; 369/378; 9; 0.085; 400,000 | S10-5 corrections | `docs/FINDINGS.md` corrections table |
| 3.2005; 2.815; 2.954; +0.139 [+0.102, +0.176] | two pairs that are not disagreements | `docs/FINDINGS.md` |
| 3.043 → 3.026; 2.534; 3.0433/3.0491/3.0456/3.0358/3.0258; 2.5342/2.2004/2.1460 | set-transformer | `s12/results/agg_dec_v2*.json`; `s26/PROPOSAL_C.md` C1 |
| 3.989; 3.213; 5.425; 6.019 | blind pipeline | `s12/SPRINT12_DOSSIER.md` |
| 204; 16; 10 | clusters spent; fresh; amyloid | `s12/SPRINT12_DOSSIER.md` |
| 1,008; 90%; 40.5%; +0.169 | QUBO | `s12/SPRINT12_DOSSIER.md` |
| 1.594; ~24; 88%; 0.388; 1.982; +0.043; 3.770; +0.558; 36.1/36.4; 62.4/72.8; 10.4° | S13 torsion (1.594 and 0.388 are ORACLE, use the native) | `s13/SPRINT13_DOSSIER.md`; `s13/cache/tors_rows.npz` |
| 54/126; 2.021; 55; −2.253; 1.486; +0.355 → +0.161 | S14 (2.021 is ORACLE, uses the native) | `s14/LEDGER.md` |
| 3.321 +0.117 [+0.050, +0.188]; 3.644 +0.440 [+0.290, +0.592]; 0.611; 11 | S15 (0.611 is ORACLE, uses the native) | `s15/FINAL_DOSSIER.md`; `s15/LEDGER.md` |
| 0.04 to 0.26; +0.004 to +0.011; +0.1554; +0.1333; 3.1606; −0.0521 | S16 | `s16/LEDGER.md` |
| +0.141; −8.4/−16.4/−17.9%; 2.609 (2.6087; 3.468; 3.454); −0.009; −0.014; +0.116; 58; 29; 27 | S17 (2.609 is ORACLE, uses the native) | `s17/LEDGER.md`; `s17/results/inband.json` |
| 4.335 +1.287; 2.411 (19 targets); −0.004 [−0.390, +0.318]; +0.191; −0.440 [−0.844, −0.133]; +0.108 [+0.011, +0.251]; 2.609; 2.573; +0.036 to +0.063; +0.722 | S18 (2.609 and 2.573 are ORACLE, use the native) | `s18/LEDGER.md` |
| 2.408 −1.202 128%; 0.688; +0.31; −0.14; 0.898/0.833/0.784/0.731/0.575/0.598; −0.110; +0.4242 | S19 | `s19/LEDGER.md` |
| 3.230; 8,192; −0.0886; +0.0048 | S20 | `s20/LEDGER.md` |
| 0.482; +0.899 | S22 | `s22/LEDGER.md` |
| −0.3403 126W/0L 98.9%; +0.142; 17/17; +0.0213; +0.0386; 0.51 to 0.64 | S23 (−0.3403 is ORACLE, uses the native) | `s23/LEDGER.md` |
| cos 0.9432; 7,329; 2.2%; 410; −2.1496; −0.4168; 0.0225 → 3.0; 2.8334 −0.2150 2.4× 113W/13L; 2.2261; 1.6024; +0.024 (cos 0.5); +0.166 (not +0.20); 102%; 7.2% | S24 (the ladder values 2.8334, 2.2261, 1.6024 are ORACLE, use the native) | `s24/LEDGER.md`; `s24/results/priorladder.json` |
| 17 values; ±1.75 Å; 51 arms +0.054; 83%; 2.0921; 350/337 | S25 misc | `s25/LEDGER.md` L6, L12, L9, L11, L14 |
| 0.986; 0.600; 0.638; 0.24 to 0.37; 0.909 | in-band ordering per target | memory `in-band-ordering-is-per-target` |
| 82.8th | native percentile of medoid criterion | memory `consensus-is-outlier-avoidance` |
| 6.43; 6.22 | core-equivalents (D7) | memory `machine-fits-two-heavy-jobs`; `core/pipeline.py` Config docstring |
| 88/90/93/95% | governor bands | `s26/governor.py` |
| 2,771,653,574 + 5,678,116,398 bytes; 8.76 GB; 4.4 GB; 7.4 GB; 16.75 GB; 64.8%; 93% | ESMFold feasibility | `s26/results/b1_feasibility.json`; `s26/PROPOSAL_B.md` B1 |
| 0.347; 1.47; 0.83; 0.08; 0.166 | cis floor | `s26/LEDGER.md` L38 |
| 1e4 kcal/mol; 40 of 75; 11; 8; 96.8% | steric reject census | `s26/LEDGER.md` L23 |

## B.7 External (audience)

| item | source |
|---|---|
| CEO pool: up to 88% / 96% / 99.6% reductions in CNOT count / depth / measurements on LiH, H₆, BeH₂ (12 to 14 qubits) | arXiv:2407.08696 abstract |
| Hessian recycling: Quantum Sci. Technol. 10, 015031 (2025) | arXiv:2401.05172 |
| Gradient troughs: arXiv:2512.25004 (Dec 2025), authors as listed | arXiv abstract |
| Co-ADAPT-VQE: up to 97% CNOT reduction on LNN devices; over 70% all-to-all | arXiv:2601.20681 abstract |
| Supervisors and affiliation | INSPIRE author 2613045; Barnes group people page |

---

# APPENDIX C. CORRECTIONS AND RETRACTIONS LEDGER

Reproduced from `docs/FINDINGS.md` ("Corrections, withdrawals and refutations"), `s14/`, `s16/`, `s18/`, `s21/`, `s24/`, and `s25/LEDGER.md`. Claim as recorded; what happened to it; where.

| claim as recorded | status | where |
|---|---|---|
| S10-1 identity audit: "0/126 tuning pools contain a ≥0.6 window", the 11/126 figure | WITHDRAWN IN FULL, VOID. Two encoders used different alphabet orders agreeing on three letters. True: 13/126 at ≥0.6, four at 1.0. | `docs/FINDINGS.md` S10-4 |
| S10-2's "10/126" leaked and "clean-116" | REVISED to 13/126 and clean-113; leak +0.0004 Å; no conclusion moves | `docs/FINDINGS.md` S10-4 |
| S8-10 transfer law `selected = 0.915 × ref + 0.346` | WITHDRAWN. In the operating band the map is the identity (1.009 × ref + 0.011, r 0.981); passing the best structure through selection costs +0.049 [−0.009, +0.103] | `docs/FINDINGS.md` S8-13 |
| S8-4 class decomposition as a "28% recoverable opportunity" | REFUTED. Injecting the library's best window is worth −0.016 [−0.037, −0.001]; the opportunity is 0.004 Å | `docs/FINDINGS.md` S8-6 and coordinator correction |
| S7-3 shrinkage diagnosis | REFUTED by its own experiment; calibration slope +0.376, not above 1 | `docs/FINDINGS.md` S7-3 |
| "ESM adds nothing for short peptides" | REFUTED; judged by MAE. On selection ESM pca32 is −0.288 [−0.484, −0.092] | `docs/FINDINGS.md` S7-11 |
| "Random fragments beat BLOSUM retrieval" | REFUTED on a matched test; BLOSUM wins at every K except a tie at 2000 | `docs/FINDINGS.md` S7-12 |
| S9-4 loop fit read as "0.49 Å per Å" | REFRAMED; the loop contracts to its own fixed point near 3.2 Å | `docs/FINDINGS.md` S9-4 |
| S9-5 positive-phi defect attributed to consensus | RELOCATED to the projection (0.2513 on the native trace) | `docs/FINDINGS.md` S9-9 |
| S10-3 hull_floor bounds every consensus operator | QUALIFIED; not those ending in the projection | `docs/FINDINGS.md` S10-5 |
| s10/forensics ceilings 2.087/2.160/1.893 | REVISED to 1.987/2.044/1.802 (solver stops short on ill-conditioned Gram) | `docs/FINDINGS.md` S10-5 |
| s10/mdgen 1.604/1.124/0.828 | RESCOPED (32 targets) to 1.802/1.336/0.953 on 126 | `docs/FINDINGS.md` S10-5 |
| "Oracle weighting cannot be worse than the hull" | REFUTED; loses on 9 of 378 cells by up to 0.085 | `docs/FINDINGS.md` S10-5 |
| S8-1 native percentile under an oracle reference | CORRECTED in place; degenerate statistic | `docs/FINDINGS.md` S8-1 |
| S7-10 interim n = 50 | SUPERSEDED by final n = 70 | `docs/FINDINGS.md` |
| S8-14 tuning-instrument sign | DID NOT REPLICATE on dev, +0.070 [−0.118, +0.259] | `docs/FINDINGS.md` |
| The shipped CVaR gradient | DEFECT; tail-only baseline; cosine −0.023; fixed | `docs/FINDINGS.md` S5-6, S9-5 |
| Sprint 6 pairwise ranker arm | CORRECTNESS BUG, found before reporting | `docs/FINDINGS.md` S6-3 |
| 3.2005 vs 3.204; 2.815 vs 2.954 | NOT disagreements: single vs multi-start; raw vs emitted | `docs/FINDINGS.md` |
| S14 C7 rank correlation +0.355 | CORRECTED to +0.161 under a uniform proposal | `s14/LEDGER.md` C7-CORRECTION |
| S14 1.486 Å figure and the coverage gate | SUPERSEDED | `s14/LEDGER.md`; memory `torsion-restraints-reach-the-target` |
| S15 eleven interim claims | RETRACTED | `s15/LEDGER.md` |
| S16 csteer abs(c) arm | RETRACTED | `s16/LEDGER.md` |
| S16 "VQE is worse than not running it" | RESCOPED to lattice-only by S20 D4 | `s20/LEDGER.md` |
| S18 L4 "weights are anti-informative" | AMENDED by L10; 1/sd² validated native-free | `s18/LEDGER.md` L10 |
| Project-wide MDE constant 0.084 Å | WRONG by 0.11× to 9.54×; MDE = 2.8016 × SE per comparison | `s21/LEDGER.md` L8 |
| Project-wide i.i.d. CI | ANTICONSERVATIVE; fold-clustered CI adopted | `s20/LEDGER.md` L-B |
| S24 L5 mechanism claim (two primaries) | WITHDRAWN by L11; generic prior reproduces 102% of β | `s24/LEDGER.md` L11 |
| S24 L9 selection lift +0.20 | CORRECTED to +0.166 | `s24/LEDGER.md` L9-A |
| S24 L15 "one of four oracles survives" | CORRECTED by L15-A: two of four; discriminator is out-of-sample validation | `s24/LEDGER.md` L15-A |
| S24 L10 "AMBER cannot select off-parallel" | AMENDED by L10-A: it does; it cannot pay for it | `s24/LEDGER.md` L10-A |
| S24 L8 "the corpus does not collapse" | AMENDED by L8-A at endpoint resolution | `s24/LEDGER.md` L8-A |
| S24 REV1 "the tail is a prefix" | OVER-CLAIM; subset-hood holds, prefix-hood violated on 1,424 cells | `s25/results/q_verify.json` |
| S25 L2 calibration mechanism | RETRACTED by L7; the null stands | `s25/LEDGER.md` L7 |
| "+0.113 Å CVaR contribution" | WITHDRAWN: wrong temperature; marginal not paired; 0.51× MDE, median 0; sign reverses at T = 0.3 | `s25/LEDGER.md` L3, L5 |
| 3.0483 Å as the headline | RE-BASED: the point cloud is an intermediate; built chain 3.2148 is the result, +0.1664 at 3.26× MDE | `s25/LEDGER.md` L4, L8, L11 |
| S25 L11's adopted guard | "necessary and not sufficient"; a synthetic leaderboard (legacy_amber_distogram 2.0921) was on disk and quarantined | `s25/LEDGER.md` L11 |
| S25 draft "the circuit attains its own optimum" | CUT; KL 0.902 nats, TV 0.453 | `s25/LEDGER.md` L15 |
| S26 L7 two point-cloud numbers | CORRECTED by L9 (written before the query returned) | `s26/LEDGER.md` L9 |
| S26 L28 test count 369/356 | CORRECTED by L32 to 370/357 | `s26/LEDGER.md` L31, L32 |
| S26 PREREG_A2 prediction dim(DLA) < 8128 | FALSIFIED; 8128 | `s26/LEDGER.md` L27 |

---

# APPENDIX D. UNSOURCED AND NON-REPRODUCING NUMBERS

None of these is used as a result in the body. Where the body mentions one, it points here.

**D1. Commit hash.** The prompt and the state brief name `ae86a124` ("refactor: deep clean and reorganize codebase"). The working tree is branch `s26` at `b7f1f280`; `ae86a124` is in its history; production code is frozen at `a15406c` and every S25 artefact carries `git_commit: a15406c82245, git_dirty: True`. Nothing in this report depends on which of the three is "the" commit, but a reader checking out `ae86a124` will not find the S26 artefacts cited in Part V.10, V.12 and the two trace files.

**D2. 3.2126 vs 3.2148.** `results/summary/leaderboard.csv` prints the production built-chain mean as 3.2126 (and `results.csv` per-target values sum to that). `bench_results/compare_tuning126.json` prints `science.rmsd_arm` 3.214765. The results lab re-projects the point cloud through `core.project.lam_path` (`s25/resultslab/exportlib.py`, `chain_from_ca`) and lands 0.0022 lower on average. The exportlib docstring itself says "The brief's Phase-II note quotes a built-chain figure of 3.2126 A. That number does not reproduce anywhere in this repository; the reproducible figures are 3.2041 (lam=0) and 3.2148 (lam=0.3)", which is contradicted by the leaderboard the same module builds. S24 L9 and S18's `exp_main_pool_512` also quote 3.2126. The body uses 3.2148 as the production figure and names 3.2126 as the rebuild. The difference is 0.07× the smallest MDE in the suite and changes no verdict. Which route is canonical is an open item (Part VIII.2).

**D3. T030 0.18242 vs 0.18198.** Same cause as D2: 0.18242 is the results-lab re-projection (`results/summary/results.csv`), 0.18198 the production cache (`s26/results/e_trace_1S9Z.json`, `reporting.rmsd_arm`; `s26/LEDGER.md` L28). Both round to 0.18.

**D4. The slide deck.** `vqe_research_overview.pptx` is not on disk, not in git history, and not in the user's Desktop, Documents, Downloads, OneDrive or profile root (`s26/LEDGER.md` L2). Part X.3 annotates the eleven-slide structure from `s26/briefs/PR.md`, not a file. Any number on the real deck that is not in Appendix B has no artefact.

**D5. The stale docstring.** `core/pipeline.py`, `quantum_stage`, still says the CVaR stage is worth "+0.113 Å by preventing the collapse". Withdrawn in `s25/LEDGER.md` L3/L5. The code does not act on the docstring; the frozen architecture forbids editing production modules; it is listed so that nobody quotes it.

**D6. "+0.370 rank correlation".** The prompt lists a withdrawn "+0.370". No artefact or ledger carries that value. The closest: `s14/LEDGER.md` C7-CORRECTION, +0.355 under a proposal drawn from the prior, corrected to +0.161 under a uniform proposal; and `docs/FINDINGS.md` S7-1's +0.367 "agreement lead", which is a different object (a difficulty thermometer, closed). The body treats +0.355 → +0.161 as the withdrawn correlation.

**D7. Core-equivalents.** The machine's 8 cores are quoted as 6.43 fast-core-equivalents in memory `machine-fits-two-heavy-jobs` and as 6.22 in the `core/pipeline.py` Config docstring. Neither has a persisted measurement file. Not used.

**D8. φ 36.1° vs 36.4°.** Quoted in `s13/SPRINT13_DOSSIER.md` and the state brief without a JSON leaf; the S26 examination listed it as document-only (C26) and the audit re-derived it from `s13/cache/tors_rows.npz` (`s26/LEDGER.md` L31, L32). Used in the body with that path.

**D9. ESM −0.288 Å.** `docs/FINDINGS.md` S7-11; the per-target artefact `s7/repr_tune.json` is lost (`s26/LEDGER.md` L11). Used in the body as a document-level number with that caveat stated.

**D10. Document-only quantum and physics numbers.** The mean |z_moment| triple 0.1127 / 0.7529 / 0.8013 (`s25/LEDGER.md` L16); the sampled tail-only cosine 0.524 (quoted inside S26 claim C24); and the test count 355/13 in older text (superseded by 370/357). Listed by `s26/LEDGER.md` L28 and L31; none is used as a result.

**D11. The identity-leak prices' original artefact.** The +0.0004 / +0.0030 prices cited `s10/idaudit_*.json`, which do not exist in the tree; the audit found them in git history at `5fa05cd` (`s26/LEDGER.md` L31). Used with that reference.

**D12. Numbers quoted from memory files.** Several statements in Parts IV, VI and VII cite `~/.claude/projects/.../memory/*.md` (for example the 82.8th percentile, the 0.986/0.600 in-band figures, the 51st percentile). They summarise sprint artefacts but are not artefacts. Each is marked "memory" inline.

