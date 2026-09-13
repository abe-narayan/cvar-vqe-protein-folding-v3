# SPRINT 26, LANE W (WILDCARD): FINDINGS

Status: Phase 0 not signed off at the time of writing; every number below is native-free or a
count, and every gated measurement is named as waiting. Branch `s26`. Ledger entries (lane W):
L30. Tiers: DEMONSTRATED (measured, artefact on disk) / ORACLE DIAGNOSTIC (reads native-derived
quantities, never selects) / HYPOTHESIS / REFUTED / OPEN. Basis is stated on every RMSD line.
No stock words; no em dashes.

---

## 0. WHAT THIS LANE OWES, AND WHERE EACH ITEM STANDS

| item | state | artefacts |
|---|---|---|
| the 2/60 benchmark self-copy leak bounded from the dev proxy (mandatory, L25) | PRE-REGISTERED; native-free halves running; gated halves wait | `s26/PREREG_selfcopy_bound.md`, `s26/IDEA_selfcopy_proxy_bound.md`, `s26/w_selfcopy.py`, `s26/w_selfcopy_test.py`, `s26/results/w_selfcopy_census.json`, `w_selfcopy_retrieval_probe_1CEK.json`, `w_selfcopy_envelope_probe_1CEK.json` |
| at least three ideas nobody else proposed | FILED (three own, plus the mandatory one) | `s26/IDEA_tiebreak_noise_floor.md`, `s26/IDEA_conformational_identity_floor.md`, `s26/IDEA_window_provenance.md` |
| the top orphaned tournament survivor | WAITING (no `s26/TOURNAMENT.md` yet) | |
| findings, ledger, status, commits | this file; L30; STATUS 09:15; commits below | |

---

## 1. THE SELF-COPY LEAK: WHAT IS ESTABLISHED BEFORE ANY RMSD IS READ

### 1.1 The mechanism has two channels and the proxy matches both. DEMONSTRATED (counts and source).

A benchmark target carried verbatim by a longer database peptide in another pinned fold meets
that peptide twice: its self-window is in the K = 500 pool (channel A) and its native distances
were training labels for the target's fold model (channel B). The retrieval library for a
target in fold k is the out-of-fold peptides plus fold-k-filtered fragments
(`s12/instrument.py` docstring; S9-10 "0 same-fold members"), and the fold-k model trains on
the out-of-fold peptides (`core/predict.py:362`), so a same-fold verbatim relative is withheld
from BOTH and a cross-fold one enters BOTH. Lane I's audit (`s26/results/i_identity_audit.json`,
`dev126_verbatim_containments`) has exactly four cross-fold cases, all in the "carried"
direction: 1CEK (13, fold 2) in a 25-mer of fold 3; 2FBU (12, fold 4) in a 23-mer of fold 0;
2P5H (9, fold 4) in a 17-mer of fold 2; 6B9K (10, fold 0) in a 17-mer of fold 2. The
benchmark's two (S24 L4) are 11-mers carried by one longer peptide in fold 3. Assumption A1 of
the PREREG rests on that.

### 1.2 Channel A cannot touch the emitted structure on two of the four. DEMONSTRATED, native-free.

`s26/results/w_selfcopy_census.json` (job `w_selfcopy_census`, exit 0, 20 s, peak RSS
0.064 GB; reads `S`, `sim`, `order`, `org` of each universe and the production `sub`; `rr` and
`nat_ca` poisoned to NaN on load):

    target  n  fold  exact self-window (universe idx / BLOSUM rank / pool pos / in shipped top-75)   windows >= 0.6 (universe / pool / top-75)
    1CEK   13   2      7 / 0 / 0 / NO                                                                 8 / 3 / 0
    2FBU   12   4   1602 / 0 / 0 / NO                                                                 4 / 1 / 0
    2P5H    9   4   3316 / 0 / 0 / YES                                                                5 / 1 / 1
    6B9K   10   0   1082 / 0 / 0 / YES                                                                7 / 1 / 1

The self-window is BLOSUM rank 0 with the maximal similarity on all four (as S24 L4 stated),
but the distogram filter keeps it in the averaged 75 on 2P5H and 6B9K only. The retrieval
probe on 1CEK (`s26/results/w_selfcopy_retrieval_probe_1CEK.json`, job
`w_selfcopy_retrieval_probe`, 15 s, 0.104 GB) confirms it end to end: with the self-window
dropped and the pool refilled, the cloud, the built chain and the lam = 0 chain are identical
(triangle bounds 0.000 on all three; the refill window, universe index 4575, does not enter
the top-75); dropping all 8 windows at >= 0.6 (3 of them in the pool) leaves the built chain
within 1.3e-7 A. 21 of 126 targets have a >= 0.6 window somewhere in the universe, S10-4's 21
reproduced from sequences alone.

### 1.3 The production emission is reproduced. DEMONSTRATED, native-free.

On 1CEK the reproduced top-75 set equals the production cache's `sub` and the cloud equals
`avg_ca` to 8.9e-16 (`gate` block of both probe files). The built chain matches the cache's
`ca` at 1.1e-4 A CA-RMSD after superposition and the lam = 0 chain matches `fit_ca` at 1.6e-4
A (internal distance matrices agree to 3.2e-4 A); the elementwise difference of up to 3.7 A is
a lab-frame difference of un-superposed chains, as L19 found between projection modes. So the
"with" arm of every part is the production emission to 1e-4 A, three orders below any bound
this lane will quote. (In-session check on the probe file; the gate function in
`s26/w_selfcopy.py` will report the superposed number once its running jobs finish and the
module may be edited again.)

### 1.4 A fold model that trained on 1CEK's own native moves its built chain by about 0.1 A. DEMONSTRATED, native-free, ORACLE-labelled operator.

`s26/results/w_selfcopy_envelope_probe_1CEK.json` (job `w_selfcopy_envelope_probe`, 25 s,
0.298 GB): the four pinned models that had 1CEK's native among their training labels (folds 0,
1, 3, 4) change the posterior by a mean |dE[d]| of 0.70 to 0.79 A per pair, keep 81 to 87% of
the top-75, and move the built chain by 0.124 / 0.079 / 0.137 / 0.122 A (RMSD between the
leaked and the clean chain; this bounds the change in RMSD-to-native by the triangle
inequality). The argmin window changes on all four (bound 0.61 to 0.69 A on the `sel` basis).
One target; the n = 126 run is in flight.

### 1.5 The pool boundary is set by ties on every target. DEMONSTRATED, native-free.

Same census file: 126/126 targets have a BLOSUM tie at the 500th rank; the boundary tie class
holds a median 115 windows (mean 117, max 200 on 7N2I), of which a median 54.5 (mean 57) are
inside the pool and 55.5 outside. About 11% of every pool is chosen by corpus order. Only 31%
of targets keep their BLOSUM rank-0 window in the shipped top-75. This is the basis of
`s26/IDEA_tiebreak_noise_floor.md`.

---

## 2. GATED, WAITING FOR "PHASE 0 SIGNED OFF" (seconds each once the native-free runs are on disk)

- `python s26/w_selfcopy.py endpoint`: signed per-target deltas on sel / cloud / arm / fit and
  the paired gain for the 4 self-copies and the 122 controls (Part A), the n = 126 envelope
  through `ST.compare` with fold-clustered CIs (Part C), Part B if its models exist.
- `python s26/w_selfcopy.py floor`: the same sequence in a different deposit, 18 targets.
- `python s26/w_selfcopy.py report`: Parts D, the bound with its assumptions and class.

Registered predictions (PREREG section 4), so that they can be judged: channel A |delta| < 0.10
A per target on the built chain (exactly 0 on 1CEK and 2FBU unless the refill enters the
top-75); the envelope helps by -0.2 to -0.8 A on the built chain against an MDE of 0.09 to
0.23 A; the cross-deposit median above 1.5 A. The 1CEK probe's 0.08 to 0.14 A triangle bounds
already say the envelope's per-target effect cannot exceed that on 1CEK; whether that is
typical is what the run measures.

---

## 3. IDEAS FILED (tournament entries; each has hypothesis, closure check, falsifier, MDE, memory, hours)

1. `s26/IDEA_selfcopy_proxy_bound.md` (mandatory direction).
2. `s26/IDEA_tiebreak_noise_floor.md`: the built-chain endpoint's sd over random tie-breaks
   of the K = 500 boundary is the pipeline's own noise floor; S17 measured only the ORACLE
   pool best (sd 0.018). 50 min CPU.
3. `s26/IDEA_conformational_identity_floor.md`: same sequence, different deposit, 18 targets
   / 22 partners; the record has four numbers (S24 L4) and no distribution. 1 min, gated.
4. `s26/IDEA_window_provenance.md`: whole-peptide / terminal / interior / fragment windows in
   the top-75; census first; plausibility 0.15.

Rejected after reading, with the closing entry: ensemble-spread floor (S7 / S12, 1.044 A,
does not predict the gap); iterated Procrustes / GPA frames (S23 L4, FINDINGS 4225); distogram
to MDS (S8-10); Rg prediction-vs-pool disagreement as a router (S21 L27/L28 signal, S23 L? router
null); test-time window ensembling (lane P's `IDEA_window_ensembling.md`).

---

## 4. WHAT DAMAGED MY OWN EXPECTATIONS

1. I expected all four self-windows to be in the averaged set because they are BLOSUM rank 0
   with maximal similarity. Two of four are filtered out by the distogram score; channel A is
   exactly zero there. The census cost 20 s and removed half of the question before any RMSD.
2. The synthetic test caught a definitional error in my `sel` triangle bound (a mean over
   cross pairs, which is not zero when the tie sets are identical); the correct bound is zero
   for identical sets and the maximum cross-pair RMSD otherwise. Written down so the next
   person does not average where a maximum is needed.
3. The reproduced built chain differs from the cache's `ca` elementwise (max 3.7 A on 1CEK)
   while the cloud matches to 9e-16. It is a lab-frame difference: 1.1e-4 A after superposition
   (section 1.3). Recorded so the gate is read as "cloud exact; chain compared after
   superposition", and so that nobody writes "chain exact" from the elementwise field.

---

## 5. WHAT I DID NOT DO, AND WHY

- Did not read any benchmark sequence, name, PDB, native or RMSD, and did not parse
  `results/benchmark_manifest.json`, `s9/final_report.json` or any file under `results/`
  naming a benchmark target (L29). The benchmark facts used are the record's 2/60 and the
  mechanism in S24 L4.
- Did not compute the triangle bound on the benchmark itself, although it needs no native
  and no RMSD, because it needs the two benchmark sequences and their universes, which
  L18b/L29 close for this sprint. Named in L30 as a coordinator's option only.
- Did not launch Part B's retrains: they are fold-model TRAINING jobs (1.25 GB, 456 s each),
  which L16b allowed for lane P under three conditions; the ruling for this lane is requested
  in L30 and PREREG section 8. Default: lane P's `pca32` fold models are the reference and
  Part B trains after sign-off.
- Did not use length-matched non-dev control chains for Part B, because naming one risks
  naming a benchmark peptide; dev-target 16-mers are the controls, and the length mismatch is
  stated.
- Did not edit `s26/w_selfcopy.py` while `w_selfcopy_retrieval` and `w_selfcopy_envelope` run
  from it (contract section 2); the superposed chain gate is a separate one-off check until
  they finish.
- Did not predict the results of the runs in flight.

---

## 6. ARTEFACT INDEX

`s26/PREREG_selfcopy_bound.md`; `s26/IDEA_selfcopy_proxy_bound.md`, `IDEA_tiebreak_noise_floor.md`,
`IDEA_conformational_identity_floor.md`, `IDEA_window_provenance.md`; `s26/w_selfcopy.py`,
`s26/w_selfcopy_test.py`; `s26/results/w_selfcopy_census.json`,
`w_selfcopy_retrieval_probe_1CEK.json`, `w_selfcopy_envelope_probe_1CEK.json` (and, when
complete, `w_selfcopy_retrieval.json`, `w_selfcopy_envelope.json`); `s26/jobs_done/w_selfcopy_*.json`;
`s26/logs/w_selfcopy_*.log`; ledger L30; STATUS 09:15.
