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
