# SPRINT 26, LANE W (WILDCARD): FINDINGS

Status: Phase 0 signed off (L33); the mandatory bound is measured and posted (L44). Branch `s26`.
Ledger entries (lane W): L30, L44. Tiers: DEMONSTRATED (measured, artefact on disk) / ORACLE DIAGNOSTIC (reads native-derived
quantities, never selects) / HYPOTHESIS / REFUTED / OPEN. Basis is stated on every RMSD line.
No stock words; no em dashes.

---

## 0. WHAT THIS LANE OWES, AND WHERE EACH ITEM STANDS

| item | state | artefacts |
|---|---|---|
| the 2/60 benchmark self-copy leak bounded from the dev proxy (mandatory, L25) | MEASURED AND POSTED (L44): pre-registered class MINOR by the own-native envelope (0.028 A built chain / 0.048 selection / 0.023 paired gain), IMMATERIAL by every direct dev-proxy measurement (0.002 A or less); Part B 1 of 10 models, the rest wait for headroom and the tournament | `s26/PREREG_selfcopy_bound.md`, `s26/IDEA_selfcopy_proxy_bound.md`, `s26/w_selfcopy.py`, `s26/w_selfcopy_test.py`, `s26/w_train_chain.py`, `s26/w_endpoint_report.py`, `s26/results/w_selfcopy_{census,retrieval,envelope,posterior,floor,endpoint,bound}.json` |
| at least three ideas nobody else proposed | FILED (three own, plus the mandatory one); the tie-break floor has PREREG, code, tests and a one-target probe; the identity floor is measured (Part E) | `s26/IDEA_tiebreak_noise_floor.md` + `PREREG_tiebreak_floor.md` + `w_tiebreak.py`, `s26/IDEA_conformational_identity_floor.md`, `s26/IDEA_window_provenance.md` |
| the top orphaned tournament survivor | WAITING (no `s26/TOURNAMENT.md` yet) | |
| findings, ledger, status, commits | this file; L30, L44; STATUS 09:15, 09:31, 19:4x; commits `8d849504`, `213a5ebb` and the closing one | |

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

### 1.6 Two discrete operators inside the readout amplify a one-member change into a 1 to 2 A move of the emission. DEMONSTRATED, native-free (seen in Part A's control population).

Dropping ONE of the 75 averaged windows (the matched control of Part A: the BLOSUM rank-0 window
removed and the pool refilled) changes the top-75 by that one member (overlap 74/75) and moves
the point cloud by 0.04 to 0.09 A on most targets, as arithmetic says it should. Two things can
turn that into a large move of the emitted structure (`s26/logs/w_selfcopy_retrieval.log`;
mechanism checked in-session on the universes, native-free):

- The MEDOID FRAME. The average is taken after superposing every member onto the medoid of the
  set; on 5H1H and 6EY3 the medoid changed (universe index 430 -> 5397; 8468 -> 7954) and the
  cloud moved 1.294 and 0.736 A, while the same 74-member set averaged in the ORIGINAL medoid's
  frame sits 0.077 and 0.064 A from the production cloud. The frame choice is a discrete
  argmin with jumps of about 1 A. (S23 L4 found iterated Procrustes does not improve the MEAN;
  the frame's discontinuity as a noise source was not priced.)
- The PROJECTION BRANCH. With the medoid unchanged and the cloud moved by 0.04 to 0.09 A, the
  built chain (ramah 0.3, multi-start) moved 0.85 A on 1NIZ, 0.96 on 1CS9, 1.52 on 1M02, 1.58
  on 1RSW, 1.75 on 2LNG and 2.11 on 2BP4 (the lam = 0 chain moved 0.07 to 0.14 A on most of
  them, 1.36 on 2BP4): a branch flip of the penalised projection (lane PH's
  `IDEA_branch_select.md` names the degeneracy; here it is seen as a sensitivity).

Both are native-free facts about the instrument and both feed `s26/IDEA_tiebreak_noise_floor.md`:
an arbitrary tie-break at the pool boundary can flip either. The full incidence over the 122
controls is reported with Part A's gated results (section 3).

### 1.7 The same sequence in a different deposit sits 2.9 A from the native. ORACLE DIAGNOSTIC (Part E, gated; run after L33).

`s26/results/w_selfcopy_floor.json` (job `w_selfcopy_floor`, exit 0, 5 s; complete 22/22,
provenance-stamped). For the 18 dev targets with a verbatim relative in the database (lane I's
list) and their 22 partners, the CA-RMSD between the target's native (model 1) and the copy of
its sequence in the other deposit (the carrier's segment at the substring position, or the
shorter relative against the target's own segment):

    median 2.908 A   mean 2.811   min 0.317 (5V5B)   max 5.502 (7S3O)
    below 1.0 A: 4 of 22 (18%)    below 1.5 A: 6 of 22 (27%)
    the four cross-fold self-copies: 1CEK 0.595, 2FBU 3.278, 2P5H 2.334, 6B9K 4.126 (S24 L4's numbers, reproduced)
    reference scales: the natives' own intra-ensemble spread 1.044 A (record, `docs/FINDINGS.md:2184`);
    the targets' K = 500 pool mean 4.34 A; the copy is WORSE than the target's pool mean on 5 of 22
    (2LNG, 5Z5W, 6B9K, 7S3O, 8ZG2) and beats the pool's best member on 3 of 22 (1ID6, 1NIZ, 5V5B);
    carrier role (n = 15) median 2.98, carried role (n = 7) median 2.04; Spearman(length ratio,
    RMSD) -0.08 (p 0.72, n = 22): no length dependence visible

The registered prediction (median above 1.5 A, fewer than a third below 1.0 A) holds; F5 did
not fire. A verbatim copy of the sequence is, in the median, no better a guess at the native
than a typical unrelated pool window is at its best. This is the mechanism behind a small leak
bound: a copy that sits 2.3 to 4.1 A from the native carries little the pipeline can exploit,
and on some targets it is a worse window than the average one.

---

## 2. THE BOUND, MEASURED (gated after L33; job `w_endpoint_report`, exit 0, 160 s, peak RSS 0.302 GB; ledger L44)

Artefacts: `s26/results/w_selfcopy_endpoint.json` (Parts A, B, C signed against the natives),
`w_selfcopy_bound.json` (Part D), `w_selfcopy_floor.json` (Part E), `w_selfcopy_posterior.json`
(Part B posteriors), built on the native-free `w_selfcopy_retrieval.json` (complete 126/126, job
`w_selfcopy_retrieval`, exit 0, 1318 s, peak 0.116 GB) and `w_selfcopy_envelope.json` (complete
126/126, job `w_selfcopy_envelope`, exit 0, 2711 s, peak 0.318 GB). Basis on every line: sel
3.4540 (argmin over K = 500) / cloud 3.0483 / arm 3.2126 (built chain; the re-projection rebuild
figure of L9, against the cache's 3.2148) / fit 3.2052 (lam = 0 chain), all reproduced by the
"with" arm on 126/126 (top-75 == `sub`, cloud max abs 1.4e-14). Sign: positive delta = the leaked
emission is WORSE. No benchmark file, sequence, native or RMSD read at any point.

### 2.1 Channel A on the four dev self-copies is 0.000 / 0.000 / -0.068 / +0.007 A on the built chain. DEMONSTRATED, PRODUCTION basis.

    target   self-window in top-75   delta arm   cloud     fit       sel       triangle bound arm   percentile in the 122-control distribution (arm)
    1CEK     no                      +0.0000     +0.0000   +0.0000   +0.0000   0.000                0.61
    2FBU     no                      +0.0000     +0.0000   +0.0000   +0.0000   0.000                0.61
    2P5H     yes                     -0.0683     -0.0114   -0.0193   +0.0000   0.202                0.89
    6B9K     yes                     +0.0071     +0.0069   +0.0046   +0.0000   0.034                0.71

On 2P5H the self-window's presence makes the built chain 0.068 A BETTER; on 6B9K 0.007 A worse.
The argmin never changes (sel 0.0000) on the four and on all 122 controls: the BLOSUM rank-0
window is never the score's argmin anywhere on the instrument, which extends S10-4's "no leaked
window is ever the pool's best" to the deployed argmin. Registered F1 holds. Matched control
population (rank-0 window dropped, pool refilled, on the other 122 targets; the dropped window was
in the top-75 on 47 of them): |delta arm| p50 0.000, p90 0.078, p95 0.147, max 0.511; mean signed
+0.006. The four self-copies sit at the 61st to 89th percentile of that distribution: an exact
self-copy is worth about as much as any other rank-0 window.

### 2.2 C27's dev half re-derived: dropping every >= 0.6 window costs +0.0004 A on the lam = 0 chain. DEMONSTRATED.

S10-4's operator on the production basis (all windows at >= 0.6 shorter-normalised identity dropped
and the pool refilled; 21 targets have one in the universe, 13 in the pool, the same 13 targets as
S10-4's table: 1CEK 1FUV 1N9U 1NIZ 1RSW 2BAO 2FBU 2P5H 6B9K 6MK8 6S0N 7N2I 8FLP), CLEAN minus
PRODUCTION, n = 126 (`ST.compare`, computed from `w_selfcopy_endpoint.json :: A/rows[*]/ge06`):

    basis   effect     fold CI95             iid CI95              MDE      W/L/T        verdict
    fit     +0.0004    [-0.0001, +0.0010]    [-0.0003, +0.0013]    0.0012   4/5/117      NOT MEASURED (0.32x MDE)
    arm     +0.0018    [-0.0003, +0.0039]    [-0.0002, +0.0046]    0.0037   5/4/117      NOT MEASURED (0.48x MDE)
    cloud   +0.0003    [-0.0002, +0.0009]                          0.0011   5/4/117      NOT MEASURED
    sel     +0.0000    exact                                       0        0/0/126      identity

S10-4's +0.0004 [-0.0004, +0.0013] on "synthesis fit" reproduces in sign, magnitude and interval,
and its 0.0000 on the shipped argmin reproduces exactly. Claim C27's dev half (L28, L31) now has an
artefact. The built-chain price is +0.0018, also inside its MDE.

### 2.3 The envelope: a fold model that trained on the target's own native emits a chain 0.70 A nearer to it. ORACLE DIAGNOSTIC, DEMONSTRATED.

The four pinned fold models that had the target's native among their training labels, each used
in place of the clean one with everything else fixed (Part C; `w_selfcopy_endpoint.json :: C`):

    basis    leaked mean minus clean   fold CI95              MDE     x MDE   W/L      folds   per leaked model (0 / 1 / 2 / 3 / 4)
    arm      -0.6980 (med -0.2663)     [-0.8312, -0.5764]     0.246   2.83    112/14   5/5     -0.749 / -0.646 / -0.702 / -0.694 / -0.700
    cloud    -0.7094                   [-0.8394, -0.5680]     0.244   2.91    110/16   5/5
    sel      -1.2081 (med -0.9239)     [-1.4377, -0.9745]     0.315   3.83    112/12   5/5     -1.324 / -1.102 / -1.254 / -1.163 / -1.201
    fit      -0.6932                   [-0.8234, -0.5615]     0.245   2.83    110/16   5/5
    arm-sel  +0.5101 (med +0.3497)     [+0.3698, +0.6838]     0.255   2.00    38/88    5/5

Registered F4 (H_C) holds with the effect inside the registered -0.2 to -0.8 A window: the shipped
MLP memorises at the pipeline endpoint. The leaked models change half the top-75 (mean overlap
0.47) and the argmin on 97% of targets; the spread among the four leaked models is 0.07 A (arm),
so the effect is the native's inclusion, not the corpus fifth; concentration sits at the 52nd
percentile of the uniform-effect null. The argmin arm gains more than the built chain, so a leak
of this strength biases the paired gain of the architecture over the shipped argmin AGAINST the
architecture (+0.51 A per leaked target). Consequences beyond the bound: (i) EXAMINATION E4's
"the fold models are not independent" has a magnitude, 0.70 A on the built chain per target whose
native a model saw; (ii) any arm that mixes fold models across targets is leaked by that much.

### 2.4 Channel B, partial (1 of 10 models): the carrier's presence is worth 0.011 A on 1CEK's built chain. DEMONSTRATED, n = 1.

`w_train_chain` was killed by the host after its first model (L40, L41). That model
(`s26/models/w_selfcopy/pca32_fold2_s0_out_1A11.pt`: fold 2 retrained through lane P's exact
pca32 path with 1A11's 276 pairs removed) against lane P's reference `pca32_fold2_s0.pt` (which
reproduces the pinned emission at 0.000 on 1CEK): posterior mean |dE[d]| 1.02 A per pair, top-75
overlap 0.76, and reference minus carrier-out = +0.0114 arm, +0.0204 cloud, -0.0062 sel: removing
the carrier makes the built chain 0.011 A worse, so its presence helped by that much. Both channels
removed on 1CEK: production is 0.0114 A better on arm and 0.0177 on the paired gain. 1CEK is the one
dev case whose copy is near-native (0.595 A, section 1.7). The nine remaining models (three
carrier-out, six control-out) wait for headroom (L42) and for the tournament; F3 is therefore
measured on one target only and its control clause is not measured.

### 2.5 The bound (Part D). PRE-REGISTERED CLASS: MINOR by the envelope, IMMATERIAL by every direct measurement.

(2 / 60) x the per-target quantity, assumptions A1 to A5 of the PREREG (`w_selfcopy_bound.json`):

    source                                                  arm       sel       paired gain   class
    dev-4 channel A, signed max                              0.0023    0.0000    0.0023         IMMATERIAL (< 0.017)
    dev-4 channel A, native-free triangle max                0.0067    0.0000    0.0067         IMMATERIAL
    both channels removed (1CEK only, n = 1)                 0.0004    0.0002    0.0006         IMMATERIAL
    own-native envelope, fold-CI limit, n = 126              0.0277    0.0479    0.0228         MINOR (0.017 to 0.170)

The Part D rule says IMMATERIAL only if all of B_real(arm), B_real(gain) and B_env(arm) are below
0.017 A; the envelope clause fires, so the pre-registered class is MINOR: 0.028 A on the built
chain, 0.048 A on the selection basis (the benchmark's `shipped` argmin arm), 0.023 A on the paired
gain. The envelope is loose by construction (a model trained on the target's own native, where the
benchmark carrier holds a copy that on the dev proxy sits 2.3 to 4.1 A from the native on 3 of 4
cases and whose one measured effect is 60x smaller). What the bound does to the benchmark verdict:
nothing. The un-leaked paired gain lies in +0.0103 +/- 0.023 under the envelope and within 0.003 of
+0.0103 under the dev-proxy measurement, against a CI half-width of 0.170. The caveat for every
benchmark figure: 2/60 self-copies, bounded at 0.028 A (built chain) by the own-native envelope and
0.002 A by the dev proxy (`s26/results/w_selfcopy_bound.json`).

### 2.6 Side measurements from the same runs.

- ORACLE insertion of the withheld same-fold carriers (8 targets, 11 windows,
  `w_selfcopy_endpoint.json :: A/insertion`): 5 of 11 enter the top-75 (F2's "at least half"
  misses by one); |delta arm| < 0.10 on 7 of 8 targets; the moves when the window enters are
  -0.047 (1NIZ, window 0.75 A from the native), -0.0002 / -0.048 / +0.006 (8IL1) and -0.277
  (8ZG2, from a window 5.06 A from the native: a projection-branch flip, section 1.6).
- The triangle bounds are loose by 3 to 10x on the control population (p95 bound 1.44 A against
  p95 signed 0.147 A): the medoid and branch flips move the chain mostly orthogonally to the
  native. The bound remains rigorous; it is not tight.
- Part E (section 1.7): median cross-deposit RMSD 2.908 A; F5 holds.

### 2.7 Falsifier scorecard (PREREG section 4)

F1 holds (max |delta arm| 0.068, max bound 0.202). F2 half: 5/11 insertions enter (needed 6), 7/8
targets under 0.10 A. F3 partial: |delta arm| 0.011 on the one carrier-out model; the control
clause is not measured. F4 holds: -0.698 A, 2.83x MDE, fold CI excludes zero, 5/5 folds. F5 holds:
median 2.91 A. Verdict of the idea: not IMMATERIAL by its own rule (the envelope), MINOR with the
numbers above; the proxy is not REFUTED (no dev self-copy moves by more than 0.10 A with both
channels removed where measured).

---

## 3. IDEAS FILED (tournament entries; each has hypothesis, closure check, falsifier, MDE, memory, hours)

1. `s26/IDEA_selfcopy_proxy_bound.md` (mandatory direction).
2. `s26/IDEA_tiebreak_noise_floor.md`: the built-chain endpoint's sd over random tie-breaks
   of the K = 500 boundary is the pipeline's own noise floor; S17 measured only the ORACLE
   pool best (sd 0.018). PREREG `s26/PREREG_tiebreak_floor.md`, code `s26/w_tiebreak.py`,
   synthetic tests ALL OK (job `w_tiebreak_test`). One-target probe, native-free
   (`s26/results/w_selfcopy_tiebreak_probe_1CEK.json`, job `w_tiebreak_probe`, exit 0, 35 s,
   peak 0.098 GB): 1CEK's boundary tie class holds 130 windows, 44 inside the pool; over 8
   seeded draws 3 to 7 of the 75 averaged members change, the argmin never does, and the
   built chain moves 0.008 to 0.031 A (triangle bounds). On this one target the floor is of
   order 0.02 A on the built chain; the 126-target run (50 min CPU, est-ram 0.5 GB) is not
   launched until the tournament ranks it (L42's headroom rule).
3. `s26/IDEA_conformational_identity_floor.md`: same sequence, different deposit, 18 targets
   / 22 partners; the record has four numbers (S24 L4) and no distribution. 1 min, gated.
4. `s26/IDEA_window_provenance.md`: whole-peptide / terminal / interior / fragment windows in
   the top-75; census first; plausibility 0.15.

Also filed for the record: the medoid-frame and projection-branch flips (section 1.6) are
the mechanism the tie-break floor would measure; the L44 control population already shows the
same one-member change moving the chain 0.5 A on 5% of targets.

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
4. I expected the own-native envelope to be a formality that the direct measurements would
   make redundant. It is the binding clause of my own Part D rule: 0.028 A against the 0.017 A
   line, so the pre-registered class is MINOR even though every direct measurement is 0.002 A
   or less. The rule was written to be conservative and it is; I report the class it gives
   and the reason the envelope is loose, not a softer class.
5. The triangle bound is rigorous and loose by an order of magnitude on this operator (p95
   bound 1.44 A against p95 signed 0.15 A on the controls): the discrete flips move the chain
   orthogonally to the native. A native-free bound that never reads the native cannot know
   the direction of the move; it was the right tool before the gate and the wrong headline
   after it.
6. Part B's chain died with the host's memory kill after one model. The one model it made is
   the informative one (1CEK, the near-native copy), by luck of the pre-registered order, not
   by design.

---

## 5. WHAT I DID NOT DO, AND WHY

- Did not read any benchmark sequence, name, PDB, native or RMSD, and did not parse
  `results/benchmark_manifest.json`, `s9/final_report.json` or any file under `results/`
  naming a benchmark target (L29). The benchmark facts used are the record's 2/60 and the
  mechanism in S24 L4.
- Did not compute the triangle bound on the benchmark itself, although it needs no native
  and no RMSD, because it needs the two benchmark sequences and their universes, which
  L18b/L29 close for this sprint. Named in L30 as a coordinator's option only.
- Launched Part B's retrains once (`w_train_chain`, after sign-off, est-ram 1.5 GB); the host
  killed the governor and the chain after the first model (L40). Nine models remain unbuilt and
  wait for the coordinator's headroom call and the tournament (L42); Part B is reported at
  n = 1 and its control clause as not measured.
- Did not use length-matched non-dev control chains for Part B, because naming one risks
  naming a benchmark peptide; dev-target 16-mers are the controls, and the length mismatch is
  stated.
- Did not edit `s26/w_selfcopy.py` after its jobs started; the gated half runs through the
  driver `s26/w_endpoint_report.py`, and the superposed chain gate stays an in-session check
  recorded in section 1.3.
- Did not run the 126-target tie-break floor (`w_tiebreak.py draws`) or the provenance census:
  neither is ranked, and L42 allows one sub-1 GB job at a time for ranked work only.
- Did not quote the envelope's 0.028 A as the leak's size: it is an upper bound from an
  operator stronger than the leak (own-native training), stated as such beside the direct
  0.002 A.

---

## 6. ARTEFACT INDEX

`s26/PREREG_selfcopy_bound.md` (+ addendum 1), `s26/PREREG_tiebreak_floor.md` (+ addendum 1);
`s26/IDEA_selfcopy_proxy_bound.md`, `IDEA_tiebreak_noise_floor.md`,
`IDEA_conformational_identity_floor.md`, `IDEA_window_provenance.md`; `s26/w_selfcopy.py`,
`s26/w_selfcopy_test.py`, `s26/w_train_chain.py`, `s26/w_endpoint_report.py`, `s26/w_tiebreak.py`,
`s26/w_tiebreak_test.py`; `s26/results/w_selfcopy_{census,retrieval,envelope,posterior,floor,endpoint,bound}.json`,
`w_selfcopy_retrieval_probe_1CEK.json`, `w_selfcopy_envelope_probe_1CEK.json`,
`w_selfcopy_tiebreak_probe_1CEK.json`; `s26/models/w_selfcopy/pca32_fold2_s0_out_1A11.pt`
(untracked, `*.pt` under `s26/models/`); `s26/jobs_done/w_*.json`; `s26/logs/w_*.log`; ledger
L30, L44; STATUS 09:15, 09:31, 19:4x.
