# S26 RETRACTIONS (lane A, the Adversary)

One entry per retraction, S26's own or of a prior sprint's claim. Format: the claim verbatim,
where it was made (file, ledger line), the artefact that contradicts it, the discriminating
experiment if one was run, and which statement now stands. Nothing superseded is deleted
anywhere; the retraction is appended. When S26 contradicts a prior sprint, both artefacts are
named and the discriminating experiment is run before anything is called retracted (contract
rule 9; campaign Part 4.3, last bullet).

At the close, the entries marked FOR docs/FINDINGS.md are handed to the coordinator for the
corrections-ledger edit; the Adversary does not edit `docs/FINDINGS.md`.

---

## R1 -- THE S26 GOVERNED TEST COUNT: 369 / 356 CORRECTED TO 370 / 357 (S26 internal)

Claim, verbatim: "Combined: 369 tests, 356 passed, 0 failed, 0 errors, 13 skipped"
(`s26/EXAMINATION.md` section D) and "the governed S26 run (`s26/TEST_RUN.md`) is 369 tests, 356
passed, 13 skipped" (`s26/EXAMINATION.md` claim C35; repeated in `s26/agentE_FINDINGS.md` X1 and
`s26/LEDGER.md` L28 item 2).

Contradicted by: `s26/results/test_run.json :: combined` = 370 tests, 357 passed, 0 failed,
0 errors, 13 skipped, 0 memory-guard skips; junit `s26/results/pytest_core_post.xml` 351/338/13
plus `pytest_amber.xml` 16/16 plus `pytest_amber_frame.xml` 3/3.

Discriminating: no experiment; the 369/356 was a read of the 00:41 render of `s26/TEST_RUN.md`,
before job `pytest_core_post` finished at 00:45:56. Found in the Adversary's Phase 0 audit
(`s26/EXAMINATION_AUDIT.md`, finding A1, ledger L31), fixed by the coordinator's appended addenda
to `s26/EXAMINATION.md` and `s26/agentE_FINDINGS.md` and by ledger L32, re-checked in L34.

Now stands: 370 tests, 357 passed, 0 failed, 0 errors, 13 skipped, 0 memory-guard skips.

---

## R2 -- "NO EXPONENTIAL PLATEAU AT n = 4..13" IS SCOPED TO DEPTH 3 (S25 claim, scoped by S26)

Claim, verbatim: the S25 trainability result that there is "no barren plateau at any width
measured" (`docs/STATE_BRIEF_2026-09-12.md` section 5.4), from the fitted log2 Var per qubit of
-0.649 (linear) to -0.311 (deployed) over n = 4..13 (`s25/results/q_plateau.json`).

Both artefacts named: `s25/results/q_plateau.json` (the slopes, which STAND as measured) and
`s26/results/q_dla.json` (A2, the dynamical Lie algebra), cross-checked by the Adversary in
`s26/results/a_dla_check.json`.

Discriminating experiment (run): A2 (ledger L27). dim(DLA) of the deployed ansatz at n = 7 is
8128 = so(128), the full real algebra, already from depth 2. Larocca 2022 / Ragone 2024 give
Var ~ 1/dim(g) only once a circuit is a 2-design over exp(g); with a maximal algebra the absence
of a decay at P = 3n = 21 parameters is a statement about SHALLOWNESS at fixed depth 3, not about
a favourable (small) algebra.

Now stands: no exponential plateau AT DEPTH 3 (P = 3n against dim so(2^n)); the algebra is
maximal, so the S25 slopes are not evidence of a favourable algebra, and the claim must carry
"at depth 3" wherever it is quoted (report, deck, trainability paper). The number is not
withdrawn; its scope is narrowed. FOR docs/FINDINGS.md at the close (a scope correction, S25 ->
S26).

---

## R3 -- LANE Q's PRE-REGISTERED DLA PREDICTION H2b IS FALSIFIED (S26 internal)

Claim, verbatim: "H2b. At the deployed cell (n = 7, L = 3) dim(DLA) < dim so(128) = 8128. Guess
from the n = 6 pattern: 4095 = dim su(64)." (`s26/PREREG_A2.md`, hypothesis H2b).

Contradicted by: `s26/results/q_dla.json :: results/fixed/n7_L3/dim` = 8128 (so(128)), reached
already at L = 2; independently reproduced at 8128 across four conjugation/gate conventions in
`s26/results/a_dla_check.json` (Adversary check, ledger L45).

Discriminating experiment (run): the exact Pauli-string Lie closure of the deployed ansatz's
generators (A2), with a dense-matrix numeric cross-check at n = 4, 5.

Now stands: the deployed RY/CNOT ansatz is controllable on the real sphere (its DLA is the full
so(2^n)) from depth 2 at every width except n = 6 and n = 9. Recorded as falsified by lane Q
itself in L27 and in the `s26/PREREG_A2.md` addendum; no number was softened.

---

## R4 -- LANE PH's PRE-REGISTERED CIS-ENSEMBLE PREDICTION IS FALSIFIED (S26 internal)

Claim, verbatim: "I predict ... non-zero cis bonds in OTHER ensemble models on some targets
(NMR ensembles are not filtered model by model)" (`s26/PREREG_cis.md` section 2, registered
expectation).

Contradicted by: `s26/results/ph_cis_census.json :: summary/ensemble` = 0 of 1,966 deposited
models across the 126 dev targets carry a cis bond (`n_models_with_cis` 0).

Discriminating experiment (run): the omega census over every deposited model (L22), an ORACLE
diagnostic reading omega angles only.

Now stands: no cis peptide bond exists anywhere on the instrument -- model 1 (0/126), every
ensemble model (0/1,966), or the retrieval pool (0/2,352,893 windows, minimum step 3.5045 A) --
by construction of the database step gate. Recorded as falsified by lane PH itself in L22 and in
the `s26/PREREG_cis.md` addendum.

---

## Notes for the close

- C27's +0.0004 / +0.0030 identity-leak prices were flagged document-only in the Phase 0 audit
  (L31, EXAMINATION C27). They are NOT retracted: the +0.0004 dev half was re-derived exactly on
  the production basis by lane W this sprint (L44), and the S10 artefacts turned out to be in git
  history at `5fa05cd` (not absent, as C27 first read). No retraction entry; a sourcing
  correction only, recorded in the audit.
- No prior-sprint accuracy claim has been contradicted by an S26 endpoint result to date; the
  physics-reject direction (L43), the C3 relaxation (L39) and the 2/60 bound (L44) all confirm
  the record rather than overturn it. If an endpoint result later contradicts a prior claim, both
  artefacts are named here and the discriminating experiment is run before the word "retracted".

---

## R5 -- L68's "L-BFGS SELECTS NO OPERATOR AT alpha = 1" IS RETRACTED BY LANE Q (L75); THE GROWTH IS INERT, NOT ABSENT (S26 internal)

Claim, verbatim: "L-BFGS growth at alpha = 1 stops with no operator selected on 78 of 78 targets"
(`s26/LEDGER.md` L68, property half) and "ADAPT with either pool selects no entangling operator on
any of the 78 alpha = 1 targets under L-BFGS" (L68 reading; `s26/PROPOSAL_A.md` sections 3 and 5:
"selects no entangling operator on any of the 78 alpha = 1 targets", "it declines them on all 78
targets; there is nothing to entangle"; `s26/agentQ_FINDINGS.md` 0.1 and 1.3). The Adversary's
L70 caveat 2 repeated it ("on the 78 alpha = 1 targets L-BFGS grows nothing and the 'P21' arm is
the 7-parameter RY layer").

Contradicted by: the A1 records themselves, `s26/results/a1/<pdb>.json :: adapt.<pool>_lbfgs_
zrank.sequence / .trace / .theta` (read correctly by lane PR, L73): under L-BFGS at alpha = 1,
operators ARE appended on 60 of 78 (pool V) and 68 of 78 (pool L2) targets, all multi-qubit, 235
and 393 strings in total; they lower the free energy by at most 1.2e-4 nats, receive angles of at
most 0.018 rad, and leave the state a product state (KL to the product of marginals at most
4.1e-4). Under Adam-best, appended on 78 of 78 (1,092 strings), worth at most 8.6e-4 nats.

Discriminating: the re-read of persisted per-target records (L75); no new experiment. The
ideal-ladder statements (L27, L35: L-BFGS appends nothing on the exact ladder) stand, because on
the ideal ladder E is exactly affine and every pool gradient is exactly zero; on a real target
tie-averaging makes E not exactly affine and the residual gradient exceeds eps = 1e-3.

Now stands: at alpha = 1 the optimum is a product state; an adaptive ansatz free to entangle
appends operators that are inert (worth under 1e-3 nats, angles under 0.02 rad, the state stays a
product to 4e-4 nats); growing anyway moves the built chain by a third of the MDE. Every endpoint
number of L68, the A2 and A4 results and the verdict REPLACE are unchanged. Corrections applied by
appended addenda to `s26/PROPOSAL_A.md` and `s26/agentQ_FINDINGS.md` (section 9), by lane PR to
slide 8 and slide 9 (L76), and by the Adversary to its own L70 caveat 2 (ledger entry below L77).

---

## R6 -- L53's "THE FIRST NATIVE-FREE QUANTITY ABOVE 0.4" IS RETRACTED BY LANE PH (L123) AFTER THE ADVERSARY'S CONTROL (L121) (S26 internal)

Claim, verbatim: "It is the first native-free quantity in this programme's record with a
correlation above 0.4 to the per-target error of the emitted structure" (`s26/LEDGER.md` L53;
`s26/agentPH_FINDINGS.md` section 3b), with the framings "the physics reports when the answer is
untrustworthy" and "how far the relaxation moves the chain predicts its error".

Contradicted by: `s26/results/a_strain_vs_spread.json` (the control L53's prereg did not carry):
the shipped top-75's own pairwise CA-RMSD spread, native-free and available before the
relaxation runs, has partial Spearman +0.452 with `rmsd_arm` (n, Rg partialled; fold CI [+0.280,
+0.609], 5/5 folds); `moved` correlates with that spread at +0.756 and, given it, adds +0.082
(iid CI [-0.103, +0.259], permutation p 0.39).

Discriminating experiment (run): the pool-spread control, job `a_strain_vs_spread` (L121).

Now stands: the pool's own disagreement is the native-free quantity that predicts the emitted
chain's error (rho +0.45); the relaxation's displacement is its proxy and adds nothing given it.
The phenomenon, the quartile table and "not a lever" stand; the sentence and the two framings
are withdrawn (L123; `s26/pr_notes.md` lines 213-216 and 414-415 replaced by lane PR, L126).

---

## R7 -- L56's "BIT-EXACT ON ALL FOUR BASES" IS CORRECTED BY LANE P (L57) (S26 internal)

Claim, verbatim: the L56 title, "the shipped posterior through the ladder's own path is
bit-exact against the production cache on all four bases".

Contradicted by: the block under L56 and `s26/results/p_ladder_shipped_s0.json`: selection and
point cloud reproduce the cache at 0.00e+00 and 7.25e-14, but the built chain lands on the
leaderboard-rebuild basis 3.2126 (the multi-start projection is sensitive to a 1e-14 input
difference), not the cache's 3.2148.

Now stands: the C2 anchor's built chain is the rebuild basis 3.2126 (named on every C2 line,
L57); selection and point cloud are bit-exact.

---

## R8 -- L99's "-0.30" ESM8M SELECTION SLOPE IS CORRECTED BY LANE P (L101) (S26 internal)

Claim, verbatim: L99's note quoting a selection slope of "-0.30" for the 8M-to-650M step.

Contradicted by: `s26/results/p_ladder_report_esm8m_s0.json :: stats/sel`: esm8m vs shipped on
selection is +0.2246, SE 0.0913, MDE 0.2557, 0.88x, fold CI [+0.122, +0.305], NOT MEASURED.

Now stands: +0.225 A on selection (0.88x MDE, not measured); the built-chain +0.242 (1.17x,
WORSE, Type-M) is unchanged (L101).

---

## R9 -- L106's "HELIX CONTENT AND THE LENGTH" IS CORRECTED BY LANE P (L107) (S26 internal)

Claim: L106 named the pool's helix content and the length as the features carrying the size of
the pipeline's gain over a constant helix.

Contradicted by: the descriptive ridge weights and correlations (`s26/agentP_FINDINGS.md`
section 10): the largest standardised weight is the top-75 members' strand (E) fraction
(-0.612; rho(ss_E, d) = -0.638); rho(n, d) = -0.128.

Now stands: the pool's strand content carries the size of the gain; the sign stays
unpredictable (L106 / L107). Sourcing note: the -0.638 is document-sourced (a findings section),
not a results JSON (L120).

---

## R10 -- L7's TWO POINT-CLOUD NUMBERS ARE CORRECTED BY LANE I (L9) (S26 internal)

Claim, verbatim: the raw coordinate average of 1D6X "is already at 1.5921 A" and of 1KWE
"2.2812" (L7).

Contradicted by: `bench_results/cache/1fc9f2dcf489e2fb/{1D6X,1KWE}.json :: rmsd_avg` = 2.0900
and 2.7313 (the entry was composed before the query returned).

Now stands: 1D6X 2.0900 / 1KWE 2.7313 on the point cloud; the L7 mechanism (the average can
land nearer the native than any single member) stands on the corrected numbers (L9). Lane I's
rule for itself is recorded there: a number enters the ledger only after its tool output is read.

---

## Prior-sprint claims that S26 contradicts, scopes or re-sources (for the docs/FINDINGS.md corrections ledger; the coordinator's L130 block carries the matching rows)

| prior claim | where | S26 disposition | S26 artefact |
|---|---|---|---|
| "no barren plateau at any width measured" (S25) | `docs/STATE_BRIEF_2026-09-12.md` 5.4; `s25/results/q_plateau.json` | SCOPED to depth 3 (P = 3n); the DLA is the full so(2^n) from depth 2, so the slopes are a shallowness statement (R2); the grown-vs-fixed slope CI is stored (L119) | `s26/results/q_dla.json`, `a_dla_check.json`, `q_var_boot.json` |
| S13 Pauli mean weights 2.236 (Legacy) / 3.015 (AMBER) | S13 / `docs/CONDENSED_REPORT.md` | NOT re-derived in S26; typed from the claim ledger by lane PR (L61) and kept off the slides (L122); the ratio leaf of `s13/results/geo_pauli.json` is sourced; a sourcing gap, not a contradiction | `s26/PRESENTATION_CHANGES.md` |
| S7-11's ESM value -0.288 A on selection (artefact `s7/repr_tune.json` lost, L11) | `docs/FINDINGS.md` S7-11 | REPLACED by the measured -0.330 A on selection (fold CI [-0.448, -0.212], 5/5) and -0.208 A on the built chain (Type-M, 5/5) (L62) | `s26/results/p_ladder_report_noesm_s0.json` |
| S10-4's identity-leak prices +0.0004 (dev) / +0.0030 (benchmark), artefacts absent (C27) | `docs/FINDINGS.md` S10-4 | the dev half RE-DERIVED exactly on the `fit` basis (+0.0004, fold CI [-0.0001, +0.0010]) and priced on the built chain (+0.0018); the benchmark half stays historical; the 2/60 bound is 0.028 A (mean-CI envelope) to 0.151 A (worst target), MINOR (L44, L55, L58) | `s26/results/w_selfcopy_endpoint.json`, `w_selfcopy_bound.json` |
| the seven-configuration suite and the Legacy +0.330 / AMBER +0.455 verdicts quoted beside 3.2148 | state brief 2, 5.3; professor brief | BASIS NAMED: point cloud; built-chain twins 3.2187 to 4.1015 (L28, L29) | `s25/results/phys_suite.json :: basis`, `results/summary/leaderboard.json` |
| S16 L27's AMBER-relaxation finding (matched random move as accurate) | `s16/` | CONFIRMED and SHARPENED on the production input: worse than a random move (+0.0111, Type-M) and than a move toward a pool member (+0.0385), replicated (L39, L46, L87); the toward-member move is a partial reversal of the projection's displacement (L95) | `s26/results/ph_c3_stage1.json`, `ph_c3_stage1_rep.json` |
| the routers (five constructions fail, S22 L7 / S23 L7) | `s22/`, `s23/` | EXTENDED: twelve m* and six s* routers on four new feature blocks null-to-harmful (L110, L115) | `s26/results/p_c4.json` |
| L38's own-torsion floor 0.347 A (this sprint) | L38 | SCOPED: an upper bound; the tight floor is 0.083 A (L89) | `s26/results/ph_cis_floor2.json` |

---

## R11 -- L68's "NULL AT THE REGISTERED THRESHOLD" IS SCOPED TO SEED 0 (L139, L140) (S26 internal)

Claim, verbatim: "qubit-ADAPT-VQE in place of the fixed ansatz is null at the registered
threshold on the built chain (-0.014 / -0.022 A, 0.23x / 0.36x MDE, fold CIs span zero)" (L68
title) and "'ADAPT is null' fires on both primaries" (L68; `s26/PROPOSAL_A.md` section 3).

Contradicted by: `s26/results/a1s1_stats.json` (seed 1, reversed fold order, L139): the same
two primaries are -0.045 / -0.052 A at 0.71x / 0.79x MDE with fold CIs excluding zero
([-0.096, -0.004], [-0.109, -0.005], 4/5 folds), outside PREREG_A1's +-0.5x "null" band; the
seed-0 points lie inside the seed-1 CIs; the between-seed change is the fixed comparator's
(3.228 -> 3.261 A), the ADAPT arms being seed-stable (3.214 / 3.216).

Discriminating experiment (run): the pre-registered seed-1 replication (PREREG_A1 addendum 2).

Now stands: A1 is NOT MEASURED on either seed (underpowered at seed 0, Type-M zone at seed 1),
its direction is consistent, and its magnitude is governed by the deployed circuit's seed; the
verdict REPLACE rests on the seed-independent facts (the product-state optimum, the full-so(2^n)
DLA, the inert growth). Corrected in `s26/PROPOSAL_A.md` addendum 3 and by the slide-8
qualifier of L140. The number is not withdrawn; the "null" label is scoped to seed 0.
