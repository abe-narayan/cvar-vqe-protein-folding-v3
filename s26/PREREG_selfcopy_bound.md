# PREREG_selfcopy_bound -- THE 2/60 BENCHMARK SELF-COPY LEAK, BOUNDED FROM THE DEV-SET PROXY (lane W, Sprint 26)

Written 2026-09-13 09:10, before any endpoint number exists. Code: `s26/w_selfcopy.py`
(synthetic unit tests `s26/w_selfcopy_test.py`). Results: `s26/results/w_selfcopy_*.json`.
Never edited after the first gated run; addenda are appended. Every command that reads a
native coordinate or an RMSD to a native refuses to run unless "PHASE 0 SIGNED OFF" is in
`s26/LEDGER.md` (enforced in code, as `s26/p_ladder.py` does). The native-free commands
(`census`, `retrieval`, `envelope`, `posterior`) poison `rr` and `nat_ca` of every loaded
universe with NaN before use and assert that every stored output is finite.

Absolute constraint honoured throughout: no benchmark sequence, name, PDB, native, RMSD or
manifest is read. The only benchmark-derived facts used are on the record already: the count
2/60 (S24 L4, L15) and the mechanism (S24 L4: two targets carried verbatim by one longer
database peptide sitting in another pinned fold, present in both the training set and the
retrieval library of both). Nothing in this protocol computes anything on a benchmark target,
and per L29 no script here parses `s9/final_report.json`, the manifest or any file under
`results/` naming a benchmark target.

## 0. What the record knows

- Declared defect (state brief section 3, ARCHITECTURE 2.1, EXAMINATION E1): identity
  normalised by the LONGER sequence; 4/126 dev and 2/60 benchmark targets carry a verbatim
  self-copy in their own fold model's training set and retrieval library.
- S10-4 (`docs/FINDINGS.md:4945`): 13/126 dev targets have a K = 500 window at >= 0.6
  identity, four at exactly 1.0 (1CEK, 2FBU, 2P5H, 6B9K); dropping all >= 0.6 windows and
  re-retrieving to K = 500 costs +0.0004 [-0.0004, +0.0013] on the synthesis `fit` basis
  (3.2005 -> 3.2009), +0.0041 on the 13, exactly 0 on the 113 clean. Artefacts
  `s10/idaudit_*.json` NO LONGER EXIST (L28 claim C27: UNSOURCED).
- S9-10 (`docs/FINDINGS.md:4573`): the same operator on the sealed benchmark, 13/60 targets
  with a >= 0.6 window, two at exactly 1.0: +0.0030 [-0.0002, +0.0061] overall, +0.0137 on
  the 13, 0.000000 on the 47 clean. That read is spent and is not repeated.
- S24 L4 (`s24/LEDGER.md:179`): on the four dev cases the self-window is BLOSUM rank #1 and
  sits 0.595 / 3.278 / 2.334 / 4.126 A from the native against pool bests of 0.342 / 2.243 /
  1.840 / 2.077: "the same peptide sequence in a different deposit adopts a different
  conformation". Option (b), quantifying the 2/60 directly, was refused because it needs
  benchmark RMSD.
- Lane I, L15 (`s26/results/i_identity_audit.json`): the four dev self-copies with carriers
  n = 25 / 23 / 17 / 17 in pinned folds 3 / 0 / 2 / 2, longer identity 0.520 / 0.522 / 0.529
  / 0.588, shorter identity 1.0; 14 further dev targets have a verbatim relative in their OWN
  fold (excluded from both training and retrieval by the fold discipline); the minimal fix
  moves exactly the four (plus 8ZG2 by transitivity); benchmark 2/60 (count only, L18).
- `core/predict.py:318` (overfit note): train NLL 0.944 against validation 2.025 at the
  shipped width and depth, 90 epochs worse than 40. The MLP fits training pairs much better
  than held-out pairs. This is the record's only prior on channel B.
- Native-free census run this turn before writing this file (`s26/w_selfcopy.py census` will
  persist it; the numbers were printed in-session from the universes' `S`, `sim`, `order`
  and the production `sub`, no `rr`): the exact self-window is universe index 7 / 1602 /
  3316 / 1082 for 1CEK / 2FBU / 2P5H / 6B9K, BLOSUM rank 0 with sim = max on all four, pool
  position 0, from a database peptide (`org` True); inside the shipped top-75 on 2P5H and
  6B9K only. Also: every one of the 126 targets has a BLOSUM tie at the K = 500 boundary
  (median 115 windows tie at the boundary score, about 57 inside and 60 outside the pool),
  which is the subject of a separate idea (`s26/IDEA_tiebreak_noise_floor.md`).

## 1. The object being bounded

For a leaked target t and a basis b in {sel (argmin over K = 500), cloud (top-75 average),
arm (built chain, PRIMARY), fit (lam = 0 chain), full (AMBER-relaxed, only if the AMBER
addendum runs)}:

    delta_t(b) = RMSD_b(emission WITH the leak) - RMSD_b(emission with BOTH channels removed)

and for the benchmark's headline paired gain g = full(or arm) - sel: delta_t(g) = delta_t(arm)
- delta_t(sel). The benchmark leak's contribution to a benchmark mean is
Delta_bench(b) = (1/60) * sum over the 2 leaked benchmark targets of delta_t(b), which is
bounded in absolute value by (2/60) * max_t |delta_t(b)| over targets whose leak matches.

The two channels, removed as follows:
- Channel A (retrieval): drop every identity-1.0 window from the universe order and take the
  next 500 (the pool is refilled from the next-best BLOSUM windows, exactly S10-4's operator,
  but restricted to EXACT copies because the benchmark's two are exact copies; the >= 0.6
  variant is run beside it to re-derive S10-4 on the production basis).
- Channel B (training): retrain the target's fold model with the carrier chain removed from
  the training corpus, everything else identical (`s26/p_ladder.py`'s pca32 path: shipped
  corpus and order, shipped MLP, 40 epochs, batch 4096, seed 0), and score the target with
  that model.

## 2. Assumptions of the bound, stated

A1. Mechanism match: the benchmark 2 are, like the dev 4, targets CARRIED verbatim by a longer
    database peptide that sits in another pinned fold, so both channels are live (S24 L4). The
    dev 4 therefore have the same channels in the same direction; the benchmark carrier
    carries both targets, the dev carriers carry one each (if two targets share a carrier the
    per-target effect is the same operator applied twice; no interaction is assumed).
A2. Exchangeability of magnitude: the per-target effect on the benchmark 2 is not larger than
    the largest per-target effect on the dev 4 (realistic bound) or than the fold-CI limit of
    the n = 126 envelope (envelope bound). n = 4 cannot estimate a quantile; "max of four" is
    reported as such, and the envelope is the guard against a fifth case being worse.
A3. Same pipeline: the benchmark was scored (S9-10) by the pre-registered pipeline the dev
    instrument runs (K = 500 out-of-fold peptides + fold fragments, top-75, coordinate average,
    ramah@0.3 projection, then AMBER); `s12.instrument` reproduces that pipeline's dev outputs
    at 0.0 (L16b, L28). The AMBER stage is not re-run here except as an optional addendum; the
    production relaxation moves the CA trace 0.220 A RMS (L24) and its effect on a per-target
    leak delta is taken to be second order, stated as an assumption.
A4. Length: the benchmark 2 are 11-mers; the dev 4 are 13 / 12 / 9 / 10-mers. No length
    dependence is assumed or corrected.
A5. The envelope (Part C) bounds channel B from above in the direction of HELP: a model trained
    on the target's own native cannot carry less native information than one trained on a
    carrier's segment of it. It does not bound HARM (a carrier whose segment is far from the
    native can mislead); harm is measured on the dev 4 directly (Part B) and its magnitude is
    bounded by Part A's + Part B's realistic maximum.

## 3. Design, in five parts

### Part A. Channel A on the dev instrument (native-free first, gated second)

Targets: the 4 self-copies (REAL, production situation) and, as the matched control
population, the other 122 targets with their BLOSUM rank-0 window dropped and the pool refilled
(same operator, an unrelated window in the dropped slot). Nine of the 122 carry a partial
(>= 0.6 by shorter-normalised identity, computed from sequences) rank-0 window (S10-4's 13
minus the 4); the control distribution is reported with and without them.

Per target, native-free (before the gate): the shipped emission (reproduced from
`I.distogram` + `I.shipped_score` over the K = 500 pool; gate: the top-75 set equals the
production cache's `sub` and the cloud equals `avg_ca`, both asserted and the max abs
difference stored) and the emission without the dropped window(s): the cloud, the built
chain (`I.project`, ramah 0.3, multi-start) and the lam = 0 chain, stored as arrays; whether
the dropped window was in the top-75; whether the refill window entered it; the top-75
overlap; the argmin tie sets; and the TRIANGLE BOUND tri_b = RMSD_b(with, without) for b in
{cloud, arm, fit} and, for sel, the RMSD between the two argmin windows (0 if the argmin set
is unchanged). Because Kabsch CA-RMSD after optimal superposition is a metric on shapes,
|RMSD(x, nat) - RMSD(y, nat)| <= RMSD(x, y): tri_b is a rigorous per-target upper bound on
|delta_t(b)| from channel A alone that never reads the native.

Gated: delta_t(b) signed, for the 4 and for the 122 controls, from the stored emissions and
`nat_ca` (read only after sign-off).

ORACLE INSERTION extension (labelled ORACLE, the fold discipline deliberately broken; no
target native read): for the 8 dev targets whose verbatim carrier sits in their OWN fold
(1KMR, 1NIZ, 1RG4, 2MD2 x2, 5Z5W, 7S3O, 8IL1 x3, 8ZG2; the carrier's CA segment taken from
the database record at the substring position), the carrier's self-window is INSERTED into
the pool in place of the 500th member and the same quantities are computed. This is what the
fold discipline withholds from those 8, and the price of a verbatim window when the
target is present: a second proxy population for channel A at n = 8 (11 windows), closer
homologs than the benchmark's (longer identity 0.6 to 0.93 against 0.52 or less), reported
beside the 4, never pooled with them.

### Part B. Channel B on the dev 4 by retraining (needs the coordinator's ruling, section 8)

For fold 2 (1CEK's), fold 4 (2FBU's and 2P5H's) and fold 0 (6B9K's), the pca32 fold model is
retrained through `s26/p_ladder.py`'s exact path with ONE chain removed from
`train_entries(fold)`: the carrier (4 models: 1A11 out of fold 2; 2LMF out of fold 4; 2P5J
out of fold 4; 1U6V out of fold 0), and, as matched controls of the training-noise kind, two
DEV-target chains per fold (public ids, chosen by rule before any run as the two longest dev
chains in that fold's training set that are not verbatim relatives of any dev target: fold 0:
9BAF, 8TXS; fold 2: 8TXS, 8T63; fold 4: 9BAF, 8T63; all 16-mers, so 105 pairs removed against
the carriers' 276 / 231 / 120 / 120: the length mismatch is stated as a limitation, chosen so
that no non-dev database member is ever named). The reference model is lane P's
`s26/models/p_ladder/pca32_fold<f>_s0.pt` (the identical function, seed and corpus with
nothing removed) once its chain has produced it, else this lane trains it under the same
function. The pinned production model is carried beside the reference as a second comparator
(lane P's gate: posterior max abs 1.5e-5, argmin 25/25).

Native-free (before the gate, once the models exist): on each of the 4 targets, the posterior
from the reference, the carrier-out and the two control-out models: max |dprob|, mean |dE[d]|
over the target's pairs, top-75 overlap, argmin agreement, and the triangle bounds tri_b
between the emissions (cloud, arm, fit). Gated: delta_t(b) signed, per model.

### Part C. The n = 126 envelope: a fold model that trained on the target's own native

For every dev target t (fold k) and each of the four pinned models j != k
(`distogram_models/fold<j>_esm_frag.pt`, whose training set contains t's native as labels
because t is a database peptide not in fold j), the pipeline is run with model j's posterior
in place of model k's, everything else identical (same K = 500 pool, top-75, average,
projection). ORACLE-LABELLED: the operator has seen the native. Native-free first: the
emissions, top-75 overlap, argmin agreement, posterior differences, and the triangle bounds
against the clean emission. Gated: `ST.compare(leaked, clean, folds, names)` on arm (PRIMARY),
cloud, sel and fit, for the mean over the four leaked models and for each j separately, with
the per-target spread among the four leaked models as the corpus-nuisance control (the four
differ from the clean model by inclusion of t AND by a different fifth of the peptide corpus;
inclusion is the only systematic part). The constant alpha-helix (`s14/results/ladder.json`
`L0_constant_helix`, 4.0648) is carried as the zero-information scale beside the effects.

### Part D. The bound

For each basis and for the paired gain: B_real(b) = (2/60) x max over the dev 4 of
|delta_t(b)| with both channels removed (the "clean" configuration: carrier-out model x
self-window dropped); B_tri(b) = (2/60) x max over the dev 4 of tri_b (native-free, both
channels); B_env(b) = (2/60) x max(|lower|, |upper|) of the fold-clustered CI of Part C's
effect. Materiality is judged against the benchmark's own 95% CI half-width 0.170 A
(S9-10: [-0.1596, +0.1803]): below 0.017 A (one tenth) is IMMATERIAL; 0.017 to 0.170 is
MINOR and is attached to every benchmark figure as a numbered caveat; above 0.170 is MATERIAL
and the benchmark verdict is restated with the leak's range.

### Part E. Why the leak is small (or not): the same sequence in a different deposit (gated)

For the 18 dev targets with a verbatim relative in the database (lane I's
`dev126_verbatim_containments`): the CA-RMSD between the target's native (model 1) and the
relative's segment at the substring position (carrier role: the carrier's segment; carried
role: the target's segment against the relative's whole chain). ORACLE DIAGNOSTIC. Reference
scales: the intra-ensemble spread of the natives (record: mean 1.044 A to model 1,
`docs/FINDINGS.md:2184`) and the target's K = 500 pool mean RMSD. Registered prediction: the
median cross-deposit RMSD exceeds 1.5 A (S24 L4's four are 0.595 to 4.126), i.e. a verbatim
copy is not a near-native answer at peptide length, which is the mechanism that keeps the
bound small. This is entered as its own idea (`s26/IDEA_conformational_identity_floor.md`)
because it also speaks to what any sequence-only predictor can reach on this instrument.

## 4. Exact falsifiers

F1 (Part A, REAL 4): |delta_t(arm)| < 0.10 A on all four and tri_arm < 0.30 A on all four.
    Falsified by any target above either threshold; then channel A is reported as measured
    and the caveat is upgraded. Registered predictions: delta = 0 exactly on 1CEK and 2FBU
    unless the refill window enters the top-75 (probability of order 75/500 per target); on
    2P5H and 6B9K one member of 75 changes and |delta| is of order 1/75 of that member's
    distance from the cloud.
F2 (Part A, ORACLE insertion 8): the inserted verbatim window enters the top-75 on at least
    half of the 11 insertions (it is BLOSUM rank 0 by construction; whether the distogram
    keeps it is the question) and |delta_t(arm)| < 0.10 A on at least 6 of 8 targets.
F3 (Part B, n = 4): native-free, the carrier-out change (max |dprob|, tri_arm) does not exceed
    both control-out changes on more than 2 of the 4 targets; gated, |delta_t(arm)| < 0.10 A
    on all four for the carrier-out model. Falsified otherwise; a falsification here means a
    single training chain moves the endpoint on its verbatim copy and the envelope of Part C
    becomes the operative bound.
F4 (Part C, n = 126): H_C says the leaked model HELPS: effect (leaked mean of four minus clean)
    negative beyond its own MDE with the fold CI excluding zero and 5/5 folds, on arm. The
    registered expectation from the overfit note is -0.2 to -0.8 A on arm. If the fold CI
    includes zero the shipped MLP does not memorise at the pipeline endpoint, channel B is
    bounded by that MDE, and the finding is stated with its power. Either outcome closes the
    channel-B question; neither is "null" without the MDE beside it.
F5 (Part E, n = 18): median cross-deposit RMSD > 1.5 A. Falsified if below 1.0 A, in which
    case verbatim copies ARE near-native on this instrument, the dev 4's small effect is
    luck, and only the envelope bound (Part C) is quoted for the benchmark.
The whole idea's verdict (Part D): IMMATERIAL if B_real(arm), B_real(g) and B_env(arm) are
all below 0.017 A; otherwise the measured class (MINOR / MATERIAL) with its numbers.

## 5. Comparison arms, basis, statistics

- Basis stated on every line: sel 3.4540 / cloud 3.0483 / arm 3.2148 (PRIMARY) / fit 3.2041.
- Part C: `ST.compare(a, b, folds=ST.pinned_folds(pdbs), names=pdbs)`, `ST.fmt`; iid and
  fold-clustered CIs, MDE = 2.8016 x SE, W/L/ties, median beside mean, concentration null;
  verdict from the fold CI and the MDE gate only. Per-model rows beside the mean-of-four row.
- Parts A and B (n = 4, n = 8): per-target signed deltas and triangle bounds tabulated; no CI
  is claimed at n = 4; each dev self-copy delta is placed as a percentile within the n = 122
  control distribution of the same operator (Part A) or against the two control-out models
  (Part B). Ties in any argmin: the tie set is stored and `ST.argmin_tied` is used for `sel`.
- ORACLE labels: Part C and the insertion arm of Part A are ORACLE by construction (an
  operator that saw the native or broke the fold discipline); Part E is an ORACLE DIAGNOSTIC;
  Parts A (real) and B and the bound itself are PRODUCTION-basis measurements of a defect.
- Replication: Part C has no random element (pinned models, deterministic pipeline). Part B's
  seed-1 replication is a pre-declared follow-up only if F3 is falsified. Part A's control
  population is the replication of its own operator.

## 6. Expected effect against the computed MDE

Part C's MDE on arm is expected between 0.09 A (a small prior change, PREREG_C2 section 4:
paired sd 0.207 on the cloud) and 0.23 A (a distogram swap through a distance-geometry fit,
paired sd 0.910); the expected memorisation effect (-0.2 to -0.8 A) would clear the upper
value. Parts A and B at n = 4 have no MDE in the usual sense; their thresholds (0.10 A per
target) are set by materiality: (2/60) x 0.10 = 0.0033 A on the benchmark mean. Expected:
channel A |delta| <= 0.05 A on the two targets whose self-window is in the top-75 and 0 on
the other two; channel B unknown, expected smaller than the envelope by the S24 L4 mechanism.

## 7. Memory, time, agent-hours (measured where a measurement exists)

- `census`: reads `S`, `sim`, `order`, `org` and the production `sub` per target; < 0.2 GB,
  2 min (the in-session run took about 1 min for 126 targets).
- `retrieval`: 126 targets x 3 emissions (with; without exact-or-rank-0; without all >= 0.6)
  + 11 insertions; each emission is one `I.project` at ~3 s: about 20 min, < 0.6 GB
  (lane P measured `I.project` and the tables at 0.594 GB peak). One-target probe first.
- `envelope`: 126 x (1 clean + 4 leaked) emissions = 630 projections, ~32 min; five pinned
  MLPs plus `esm_small.npz` (21 MB): < 0.8 GB expected; probe on one target first.
- `train`: 1.248 GB peak, 456 s per fold model (lane P, `s26/jobs_done/p_probe_train_pca32_f0.json`);
  10 models (4 carrier-out, 6 control-out) = 76 min one at a time, est-ram 1.5 GB.
- `posterior`, `endpoint`, `floor`, `report`: seconds to minutes; endpoints read stored arrays.
- Agent-hours: 3 code and tests, 2 runs and probes, 2 write-up.

## 8. Operator forks declared, and the ruling requested

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| which windows are "the leak" | identity-1.0 windows only (the benchmark's two are exact copies) | all >= 0.6 (S10-4's set; run beside it as the C27 re-derivation, not as the bound) |
| refill | the next-best BLOSUM windows to K = 500 (S10-4's operator, the production rule) | leaving the pool at 499; re-running retrieval on a rebuilt library |
| channel B comparator | lane P's retrained pca32 reference (identical function); the pinned model carried beside it | the pinned model alone (a 1.5e-5 retraining difference would be confounded with the removal) |
| control chains | dev-target chains, public ids, longest available (16-mers) | length-matched non-dev chains (they cannot be named without risking a benchmark id; ruling requested below) |
| bound statistic | max over the 4; fold-CI limit of the envelope; both reported | a mean over 4 with a CI (n = 4 cannot support one) |
| basis | arm PRIMARY; sel, cloud, fit carried; paired gain arm - sel | the AMBER-relaxed emission (optional addendum under tag AMBER if the coordinator wants it; 8 relaxations) |

Ruling requested from the coordinator (rule 14; the lane does not run Part B's training
until answered): may the 10 channel-B retrains (tag CPU, est-ram 1.5 GB, one governor job at
a time, order fixed here: carrier-out 1A11/f2, 2LMF/f4, 2P5J/f4, 1U6V/f0, then control-out
9BAF/f0, 8TXS/f0, 8TXS/f2, 8T63/f2, 9BAF/f4, 8T63/f4) run before sign-off under L16b's three
conditions, and may they interleave with lane P's chain (two training jobs resident at once,
about 2.5 GB, within the band at the 72% baseline read at 08:54)? Also: may a length-matched
non-dev control chain be used if it is referred to only by length and sequence hash? Default
if unanswered: the dev-chain controls, training after sign-off.

## 9. What I will do if the bound is not immaterial

Report the measured class with its numbers on every basis, append the caveat text for the
benchmark figures to `s26/agentW_FINDINGS.md` and the ledger, run Part B at seed 1, and hand
the result to the Adversary. I will not open the benchmark to check.

## ADDENDUM 1 (2026-09-13 19:45) -- what ran, what deviated; nothing above edited

- Native-free halves complete: `retrieval` (126/126, job `w_selfcopy_retrieval`, exit 0, 1318 s,
  peak 0.116 GB), `envelope` (126/126, `w_selfcopy_envelope`, exit 0, 2711 s, peak 0.318 GB),
  `census` (0.064 GB). Gated halves in one job after L33: `w_endpoint_report` (posterior,
  endpoint, report; exit 0, 160 s, peak 0.302 GB); `floor` separately (exit 0, 5 s).
- Deviation 1: Part B ran 1 of its 10 models. `w_train_chain` (est-ram 1.5 GB) built
  `pca32_fold2_s0_out_1A11.pt` and was then killed with the governor by the host for low memory
  (L40); the nine others are not built. Part B is reported at n = 1 (1CEK, carrier-out only), F3's
  control clause is NOT MEASURED, and `w_selfcopy_bound.json :: signed_bounds_gated/both_removed4`
  holds one target despite its key name.
- Deviation 2: the reference model for every fold is lane P's `s26/models/p_ladder/pca32_fold<f>_s0.pt`
  (trained by the same function after this PREREG was written); it reproduces the pinned emission
  at 0.000 on all four self-copy targets, so the "pinned beside reference" comparator is moot.
- Deviation 3: the reproduction gate's chain comparison is made after superposition (1.1e-4 A on
  1CEK; elementwise it is a lab-frame difference of up to 3.7 A). The cloud gate is exact.
- Results and the pre-registered verdict: ledger L44; `s26/agentW_FINDINGS.md` section 2. Class
  MINOR by the envelope clause (0.028 A built chain), IMMATERIAL by every direct measurement.
- Not done: replication of Part B at seed 1 (nothing positive to replicate at n = 1); the AMBER
  addendum (section 8, optional) was not requested.
