# SPRINT 26, LANE W (WILDCARD): FINDINGS

Status: CLOSED at the sprint's end (2026-09-14 04:2x). Branch `s26`. Ledger entries (lane W): L30, L44,
L52, L58, L64, L84, L85, L90, L91, L105, L108, L109 and the closing entry L137 (what ran, what did
not, and the time each open item would need). Every governed job of this lane exited 0 except the
first synthetic test run (fixed before any real run); largest peak RSS 1.250 GB. Tiers: DEMONSTRATED (measured, artefact on disk) / ORACLE DIAGNOSTIC (reads native-derived
quantities, never selects) / HYPOTHESIS / REFUTED / OPEN. Basis is stated on every RMSD line.
No stock words; no em dashes.

---

## 0. WHAT THIS LANE OWES, AND WHERE EACH ITEM STANDS

| item | state | artefacts |
|---|---|---|
| the 2/60 benchmark self-copy leak bounded from the dev proxy (mandatory, L25) | MEASURED AND POSTED (L44): pre-registered class MINOR by the own-native envelope (mean-CI limit 0.028 A built chain / 0.048 selection / 0.023 paired gain; worst single target 0.151 / 0.194 / 0.123 under assumption A2, L55), IMMATERIAL by every direct dev-proxy measurement (0.002 A or less); Part B 1 of 10 models, the rest wait for headroom and the tournament | `s26/PREREG_selfcopy_bound.md`, `s26/IDEA_selfcopy_proxy_bound.md`, `s26/w_selfcopy.py`, `s26/w_selfcopy_test.py`, `s26/w_train_chain.py`, `s26/w_endpoint_report.py`, `s26/results/w_selfcopy_{census,retrieval,envelope,posterior,floor,endpoint,bound}.json` |
| at least three ideas nobody else proposed | FILED (three own, plus the mandatory one); the tie-break floor has PREREG, code, tests and a one-target probe; the identity floor is measured (Part E) | `s26/IDEA_tiebreak_noise_floor.md` + `PREREG_tiebreak_floor.md` + `w_tiebreak.py`, `s26/IDEA_conformational_identity_floor.md`, `s26/IDEA_window_provenance.md` |
| the tournament's own ideas (L51: items 2 and 4) | item 2 conformational_identity_floor MEASURED AND POSTED (L52); item 4 tiebreak_noise_floor MEASURED AND POSTED (L64): floor 0.004 A on the 126-mean, 0.024 A paired MDE between two conventions, every hundredths-level recorded effect inside it | `s26/PREREG_identity_floor.md`, `s26/w_identity_floor_stats.py`, `s26/results/w_identity_floor.json`; `s26/PREREG_tiebreak_floor.md`, `s26/w_tiebreak.py`, `s26/w_tiebreak_report.py`, `s26/results/w_selfcopy_tiebreak_{draws,endpoint}.json`, `w_tiebreak_report.json` |
| the top orphaned survivor: window_ensembling (lane P's idea, L51 item 6) | MEASURED AND POSTED (L84): refuted with power; ens3 minus shipped -0.0005 A (0.02x MDE), against the zero-information resample +0.0042 (0.14x); a gain of 0.028 A or more excluded; the fixed-K form of the mandatory direction is closed | `s26/PREREG_window_ensembling.md`, `s26/w_ensemble.py`, `s26/w_ensemble_test.py`, `s26/results/w_selfcopy_ensemble_{probe_1A13,clouds,endpoint}.json` |
| the remaining orphans and extensions (L77) | window_provenance POSTED in full (L85 census and ORACLE contrast; L109 readout H_P3 refuted with power); amber_prior_partner POSTED (L105: the LFO choice declines to mix on all five folds; every cell worse); Part B complete on the four carrier-out models (L108; F3 falsified by 2P5H in the harmful direction; direct bound 0.008 A, IMMATERIAL), control-out models 1 of 6 built and the rest queued; partial_recall_gradient MEASURED (L90: no gradient at the MDE, I_long rho -0.205 suggestive) and memorisation_on_the_ladder MEASURED (L91: the discounted ladder under-prices the trained operator by 0.35 A, 1.5x MDE) | `s26/PREREG_{window_provenance,amber_prior_partner,partial_recall_gradient,memorisation_on_the_ladder}.md`, `s26/w_{provenance,amberprior,recall,ladder}.py` |
| ideas found on the way (coordinator's request) | FILED: the partial-recall gradient below the 0.6 threshold; the own-native model placed on the S24 prior ladder | `s26/IDEA_partial_recall_gradient.md`, `s26/IDEA_memorisation_on_the_ladder.md` |
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

### 2.4 Channel B on all four carrier-out models: the carrier's presence is worth +0.011 / -0.002 / -0.246 / -0.027 A on the built chain; F3 falsified by 2P5H in the harmful direction. DEMONSTRATED, n = 4 (ledger L108).

The four carrier-out models (1A11 out of fold 2, built before the host kill of L40; 2LMF out of
fold 4, job `w_train_out_2LMF`, 1049 s, peak 1.244 GB; 2P5J out of fold 4, 1063 s, 1.250 GB; 1U6V
out of fold 0, 893 s, 1.248 GB), each through lane P's exact pca32 path with one chain's pairs
removed, against lane P's reference fold models (which reproduce the pinned emission at 0.000);
gated re-run `w_endpoint_report2` (145 s, 0.309 GB), the L55 additions re-applied by
`s26/w_bound_addendum.py`. Delta = reference minus carrier-out; positive = the carrier helped.

    target   carrier   |dE[d]| per pair; top-75 overlap   arm       cloud     fit       sel       gain      copy's distance from the native (L52)
    1CEK     1A11      1.02 A; 0.76                       +0.0114   +0.0204   +0.0115   -0.0062   +0.0177   0.60 A
    2FBU     2LMF      0.79 A; 0.89                       -0.0021   -0.0002   -0.0021   -0.0040   +0.0019   3.28 A
    2P5H     2P5J      1.38 A; 0.77                       -0.2460   -0.1634   -0.2310   +0.0000   -0.2460   2.33 A
    6B9K     1U6V      1.37 A; 0.91                       -0.0265   -0.0318   -0.0188   +0.0000   -0.0265   4.13 A

F3's first clause (|delta arm| < 0.10 on all four) is falsified by 2P5H: training the fold model
on the carrier's copy of 2P5H's sequence, whose geometry is 2.33 A from 2P5H's native, moves the
built chain 0.25 A AWAY from the native. The training channel's sign follows the cross-deposit
distance of the copy (helps at 0.60 A, hurts at 2.3 to 4.1 A): a memorising model with a copy
that is not the native is a bias toward another deposit's conformation, not a gift. Both
channels removed (self-window dropped and carrier out): production minus clean +0.0114 / -0.0021
/ -0.2460 / -0.0243 on the built chain. F3's control clause (against two control-out chains per
fold) is measured on one control so far (9BAF/f0, 1440 s, 1.247 GB); the rest are queued and go
into an addendum.

### 2.5 The bound (Part D). PRE-REGISTERED CLASS: MINOR by the envelope, IMMATERIAL by every direct measurement.

(2 / 60) x the per-target quantity, assumptions A1 to A5 of the PREREG (`w_selfcopy_bound.json`):

    source                                                  arm       sel       paired gain   class
    dev-4 channel A, signed max                              0.0023    0.0000    0.0023         IMMATERIAL (< 0.017)
    dev-4 channel A, native-free triangle max                0.0067    0.0000    0.0067         IMMATERIAL
    both channels removed (n = 4 since L108; 2P5H)           0.0082    0.0002    0.0082         IMMATERIAL
    own-native envelope, fold-CI limit of the MEAN, n = 126  0.0277    0.0479    0.0228         MINOR (0.017 to 0.170)
    the same envelope at the p95 target                      0.0830    0.1154    0.0816         MINOR
    the same envelope at the WORST single target             0.1512    0.1944    0.1232         MINOR on arm and gain; sel crosses 0.170 (9KAR)
                                                             (2BP4)    (9KAR)    (2NBC)         (added per L55, `signed_bounds_gated/C_envelope_per_target`)

The Part D rule says IMMATERIAL only if all of B_real(arm), B_real(gain) and B_env(arm) are below
0.017 A; the envelope clause fires, so the pre-registered class is MINOR: 0.028 A on the built
chain, 0.048 A on the selection basis (the benchmark's `shipped` argmin arm), 0.023 A on the paired
gain. Those three are fold-CI limits of a MEAN effect (Part D as pre-registered, assumption A2);
the Adversary (L55) asked for the per-target reading of the same envelope beside them, and it is:
(2/60) x |leaked-mean minus clean| at the worst single target 0.151 A on the built chain (2BP4,
4.54 A), 0.194 A on selection (9KAR) and 0.123 A on the paired gain (2NBC); at the p95 target
0.083 / 0.115 / 0.082. Under every reading the class on the built chain and on the paired gain is
MINOR; on the selection basis the worst-target reading crosses 0.170 (MATERIAL by the
pre-registered line) for the one benchmark arm most sensitive to memorisation, the `shipped`
argmin. The expected contribution stays 0.008 A (dev proxy, both channels, n = 4 since L108; it was
0.002 at n = 1) to 0.023 A (mean envelope), and where the direct measurement is not zero it is
HARMFUL to the leaked target's built chain and neutral to its argmin, so an un-leaked benchmark
paired gain would if anything be slightly more favourable to the architecture than +0.0103. The
wording to carry: "expected contribution 0.008 to 0.023, mean-CI limit 0.028, worst single
target 0.151, under A2", with the envelope named as the over-bound it is (own-native training, 60x the one
measured carrier effect). The envelope is loose by construction (a model trained on the target's own native, where the
benchmark carrier holds a copy that on the dev proxy sits 2.3 to 4.1 A from the native on 3 of 4
cases and whose one measured effect is 60x smaller). What the bound does to the benchmark verdict:
nothing. The un-leaked paired gain lies in +0.0103 +/- 0.023 under the envelope and within 0.003 of
+0.0103 under the dev-proxy measurement, against a CI half-width of 0.170. The caveat for every
benchmark figure (L55's wording, the dev price updated by L108): 2/60 self-copies; dev-proxy price 0.008 A; own-native envelope
0.028 A (mean CI) to 0.151 A (worst target) on the built chain, under A2; MINOR under every
reading on the built chain and the paired gain; cannot move the benchmark verdict either way
(`s26/results/w_selfcopy_bound.json`, regenerated by `s26/w_bound_addendum.py` with the gain
envelope row and the per-target readings; the original provenance kept under `provenance_original`).

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

## 2b. THE TIE-BREAK NOISE FLOOR, MEASURED (tournament item 4; ledger L64)

`s26/PREREG_tiebreak_floor.md`; jobs `w_tiebreak_draws` (exit 0, 5166 s under a four-job load, peak
0.111 GB; `s26/results/w_selfcopy_tiebreak_draws.json`, complete 126/126, gate top-75 == `sub`
126/126) and `w_tiebreak_endpoint` (exit 0, 25 s, 0.062 GB; `w_selfcopy_tiebreak_endpoint.json`);
ledger numbers by `s26/w_tiebreak_report.py` -> `w_tiebreak_report.json`. Every target's K = 500
boundary falls inside a BLOSUM62 tie class (median 115 windows, 54.5 inside the pool); the
production convention breaks it by corpus order. Eight seeded uniform re-draws of the boundary
members per target, everything downstream identical.

### 2b.1 A random tie-break replaces 3.6 of the 75 averaged members and the argmin on 8.5% of cells. DEMONSTRATED, native-free.

Median 3 replaced (max 15); one target never changes; the chain moves 0.168 A per draw in the
median (triangle bound; p90 1.08 A), seven times more than its RMSD-to-native changes: the moves
are mostly orthogonal to the native, as in the L44 controls. The registered "about 8 of 75" was
an over-estimate: the score's top-75 is drawn mostly from the non-tied prefix.

### 2b.2 The floor: 0.004 A on the 126-mean, 0.024 A as the paired MDE between two conventions (built chain). DEMONSTRATED.

    basis   m_tie (sd of the 126-mean, 8 draws)   s_tie median / p90   paired MDE between two draws, median / max (28 pairs)   production vs mean of draws
    arm     0.0039                                 0.0227 / 0.1290      0.0236 / 0.0321                                          -0.0026, 0.15x MDE, NOT MEASURED
    cloud   0.0017                                 0.0126 / 0.0408      0.0100 / 0.0124                                          -0.0048, 0.55x MDE, NOT MEASURED
    fit     0.0032                                 0.0148 / 0.0745      0.0157 / 0.0187                                          -0.0047, 0.36x MDE, NOT MEASURED
    sel     0.0136                                 0.0000 / 0.0901      0.0466 / 0.0714                                          -0.0201, 0.75x MDE, NOT MEASURED (fold CI [-0.039, -0.002], 107 ties)

Per target on the built chain: s_tie above 0.1 A on 21/126, above 0.5 A on none; range over the 8
draws median 0.065, p90 0.362, max 0.713 (1D6X, the `pool_gate = WARN` target, on a branch flip).
Registered predictions: m_tie 0.003 to 0.010 (0.0039, holds); paired MDE 0.02 to 0.05 (0.0236,
holds); median s_tie 0.03 to 0.08 (0.0227, just below); ~8 of 75 replaced (3.6, below).

### 2b.3 Every hundredths-level effect on the record is inside the paired MDE between two tie-breaks. DEMONSTRATED (the idea's confirmation clause).

All nine recorded effects listed in `w_tiebreak_report.json :: recorded_effects_vs_floor` (C3's
+0.0207 and +0.0111, the S24 restrained relaxation -0.022, the ORACLE functional weight 0.015,
the all-atom reranking +0.004, C27's +0.0004 / +0.0018, the fd-vs-exact projection 0.012, L24's
+0.013) are below 0.0236 A; three (C27's two prices, the reranking) are also below 2 x m_tie =
0.008 A, the spread of the mean that the convention alone produces. What that means and does not
mean: the floor applies to any contrast whose two arms do NOT share the tie-break (a re-run after a
change of the pinned `pdbs/` order, a different K, a re-retrieval, a comparison across sprints or
instruments; S9's 0.004 to 0.005 A stable-vs-plain-argsort residual is this floor on the mean); it
does NOT apply to a paired contrast whose arms share the pool, which is how C3, the functional
weight and the reranking were measured, so their paired SEs already exclude it. The report's
sentence: a hundredths-level effect is real only as a paired contrast with the tie-break held
fixed; quoted across runs it is inside the pipeline's own convention noise (0.024 A at n = 126).

### 2b.4 The production convention is not distinguishable from a random draw. DEMONSTRATED (secondary).

Built chain -0.0026 (0.15x MDE), point cloud -0.0048 (0.55x), lam = 0 chain -0.0047 (0.36x): all
NOT MEASURED. On the selection basis the corpus-order convention's argmin is 0.020 A better than
a random draw's at 0.75x MDE, fold CI [-0.039, -0.002], 13W/6L/107T: suggestive of the S25
finding that the retrieval order is not neutral (rho(pool index, ORACLE RMSD) +0.054), not a
result, and not a lever (a tie-break cannot be chosen native-free).

---

## 2c. THE ORPHAN: FIXED-K WINDOW ENSEMBLING, MEASURED (tournament item 6; ledger L84)

`s26/PREREG_window_ensembling.md` (written before the probe); jobs `w_ensemble_probe` (1A13, 35 s,
0.097 GB), `w_ensemble_clouds` (exit 0, 4351 s wall under a four-job load with two governor
suspensions, peak 0.113 GB; `s26/results/w_selfcopy_ensemble_clouds.json`, complete 126/126,
gate 126/126) and `w_ensemble_endpoint` (exit 0, 10 s, 0.057 GB;
`w_selfcopy_ensemble_endpoint.json`). Basis: built chain PRIMARY (rebuild basis 3.2126), point
cloud carried; no selection basis exists for an ensemble.

### 2c.1 How it differs from widening K, and why the record predicted null. STATED BEFORE RUNNING.

Widening K (S17 L12) draws ONE shortlist from a wider pool and the extra plausible windows
displace near-native members out of that single top-75 (its ORACLE best 2.104 -> 2.572 A).
Fixed-K ensembling keeps three shortlists of K = 500 (BLOSUM45 / 62 / 80), each cut to its own
top-75 by the shipped score (the BLOSUM62 one is the production set, untouched), and averages
the three CLOUDS in the production frame: no shortlist is widened and displacement cannot act.
What can act is whether the three clouds' errors are parallel; S23 L9 (68% common-mode error
shared by the whole universe) and S24 L3 (score selection makes bias parallel across sources,
cosine 0.943 against 0.933 within-source) predicted that they are.

### 2c.2 The three shortlists overlap at 0.83 to 0.88 and their clouds sit 0.14 to 0.20 A apart. DEMONSTRATED, native-free.

Medians over 126: pool overlap 0.83 (45 vs 62) / 0.89 (80 vs 62); top-75 overlap 0.83 / 0.88;
93 distinct members in the union of three; cloud separations 0.18 / 0.14 / 0.20 A; the ensemble
moves the built chain 0.185 A from the production chain, the zero-information resample of the
production top-75 (`boot3`) 0.249 A.

### 2c.3 The ensemble is -0.0005 A on the built chain (0.02x MDE) and +0.0042 against its zero-information control (0.14x). REFUTED with power.

    contrast (built chain)        effect     MDE      fold CI95              folds   W/L      verdict
    ens3 minus shipped            -0.0005    0.0281   [-0.0194, +0.0162]     3/5     64/62    NOT MEASURED (0.02x)
    ens3 minus boot3              +0.0042    0.0297   [-0.0178, +0.0300]     2/5     61/65    NOT MEASURED (0.14x)
    b45 minus shipped             -0.0039    0.0405   [-0.0348, +0.0212]     2/5     60/66    NOT MEASURED
    b80 minus shipped             +0.0237    0.0337   [+0.0052, +0.0402]     4/5     57/69    NOT MEASURED (0.70x; Type-M zone; the BLOSUM80 shortlist alone is if anything worse)
    ensK minus shipped            -0.0153    0.0386   [-0.0380, +0.0047]     3/5     64/62    NOT MEASURED (0.40x)
    union3 minus shipped          -0.0021    0.0191   [-0.0180, +0.0168]     3/5     67/59    NOT MEASURED
    boot3 minus shipped           -0.0047    0.0273   [-0.0161, +0.0045]     3/5     67/59    NOT MEASURED
    point cloud: ens3 minus shipped +0.0016 (0.09x); ens3 minus boot3 -0.0008 (0.04x); ensK minus shipped -0.0116 (0.52x, fold CI [-0.0179, -0.0045], 5/5 folds: suggestive on the cloud only, not carried by its own chain)

H_E is refuted: the design excludes a gain of 0.028 A or more on the built chain (0.019 on the
cloud), and the registered expectation (0.00 to +0.03) held. Mechanism as registered: averaging
parallel errors returns the same error, and the ensemble's move of the emission is the size and
the value of a resample of the production set. The one cell with a fold CI excluding zero in the
helpful direction is the K-variant ensemble on the point cloud (-0.012, 0.52x MDE), the S17 L12
direction (the consensus readout improves with K on the cloud; the built chain does not follow).
The mandatory test-time-ensembling direction is closed in the fixed-K form with this power
statement; the widening-K form was closed by S17 L12. No deviation from the PREREG.

---

## 2d. WINDOW PROVENANCE: THE CENSUS AND THE ORACLE CLASS CONTRAST (tournament item 9; ledger L85; the readout test H_P3 running)

`s26/PREREG_window_provenance.md`; jobs `w_provenance_census` (exit 0, 5 s, 0.004 GB;
`s26/results/w_selfcopy_provenance_census.json`, complete 126/126, 0 unresolved) and
`w_provenance_oracle` (exit 0, 10 s, 0.055 GB; `w_selfcopy_provenance_oracle.json`). Class rule
tested synthetically (`s26/w_provenance_test.py` ALL OK).

### 2d.1 Three quarters of every pool and of every top-75 is fragment windows; whole peptides are 0.6%. DEMONSTRATED, native-free.

    class       K = 500 pool   shipped top-75   targets with one in the top-75
    whole          0.6%            0.7%           30 / 126
    terminal       8.0%            7.2%          113 / 126
    interior      18.7%           18.3%
    fragment      72.7%           73.8%

The score keeps the pool's class mix (it neither favours nor removes peptide windows).

### 2d.2 In the pool a peptide-derived window is 0.42 A nearer the native than a fragment window; inside the top-75 the gap is 0.09 to 0.13 A. ORACLE DIAGNOSTIC (single-window basis).

    where    contrast (first minus second class, per-target mean rr)   effect    MDE     fold CI95             folds   verdict
    pool     peptide_any minus fragment (n = 126)                       -0.416    0.131   [-0.447, -0.376]      5/5     BETTER (3.2x)
    pool     whole+terminal minus fragment (126)                        -0.421    0.136   [-0.453, -0.388]      5/5     BETTER (3.1x)
    pool     interior minus fragment (126)                              -0.411    0.136   [-0.447, -0.368]      5/5     BETTER (3.0x)
    pool     whole+terminal minus interior (126)                        -0.011    0.085   [-0.055, +0.027]      3/5     NOT MEASURED (0.13x)
    top-75   whole+terminal minus fragment (114)                        -0.131    0.116   [-0.165, -0.110]      5/5     BETTER, TYPE-M ZONE (1.13x)
    top-75   peptide_any minus fragment (126)                           -0.086    0.090   [-0.116, -0.056]      5/5     NOT MEASURED (0.96x)
    top-75   interior minus fragment (125)                              -0.061    0.093   [-0.104, -0.022]      4/5     NOT MEASURED (0.66x)
    top-75   whole+terminal minus interior (113)                        -0.056    0.114   [-0.128, +0.007]      4/5     NOT MEASURED (0.49x)

Any peptide window (whole, terminal or interior alike) is 0.41 to 0.42 A nearer the native than a
fragment window in the pool, 5/5 folds, three times its MDE: S19's "the peptide corpus carries
the sequence-structure channel" and S24 L8's "a fragment's conformation is held by contacts
outside the window", measured for the first time on the single-window basis of the shipped
pool; the window's position in its parent does not matter. The score removes most of it: inside
the top-75 the remaining gap is at the edge of what n = 114 to 126 can see. By the PREREG's rule
the achievable readout test H_P3 ran (`s26/w_provenance_readout.py`, job `w_provenance_readout`,
9321 s under load, 0.110 GB; ledger L109): the fragment-class weight chosen leave-fold-out on the
built chain is 0.5 on three folds and 0 on two (never uniform, never up-weighting), and the routed
readout is +0.0066 A against uniform on the built chain (0.15x MDE 0.044; +0.0001 on the cloud) and
-0.0179 against the permuted-weight control (0.38x; -0.0284, 0.68x, on the cloud, 5/5 folds). H_P3
is refuted; a gain of 0.044 A or more is excluded; the registered expectation held. The 0.13 A
single-window gap acts on 8% of a 75-member mean and, as the L44 / L64 controls showed, a few
members' change moves the chain mostly orthogonally to the native.

---

## 2e. THE PARTIAL-RECALL GRADIENT: NONE AT THE MDE; A SUGGESTIVE rho -0.2 THAT RETRIEVAL CAN EXPLAIN (ledger L90)

`s26/PREREG_partial_recall_gradient.md`; jobs `w_recall_cov` (15 s, 0.292 GB;
`s26/results/w_selfcopy_recall_covariates.json`) and `w_recall_endpoint` (45 s, 0.111 GB;
`w_selfcopy_recall_endpoint.json`). Native-free covariates against each target's own fold
model's training corpus: I_long (max pinned identity; 0.33 to 0.59, all below 0.6), I_short (max
containment; 0.50 to 1.00), L_kmer (longest shared substring, 3 to 13). Against `rmsd_arm`:
I_long rho -0.205 (fold CI [-0.315, -0.099], permutation p 0.010, partial on length -0.176),
I_short -0.157, L_kmer -0.088; the registered MDE in rho is 0.253. No gradient exists by the
pre-registered rule (rho <= -0.25); I_long sits at 0.81 to 0.91x the MDE on all three bases with
the fold CI excluding zero: suggestive, not measured. A clean bill for the 0.6 fold threshold at
this resolution; power: |rho| >= 0.25 excluded. Confound stated in the entry: the training
corpus is the retrieval library, so the covariate is also a retrieval-proximity covariate, and a
negative rho is expected from retrieval alone; the mechanism check that would separate the two
was gated on a gradient at the MDE and did not run.

## 2f. THE OWN-NATIVE MODELS ON THE S24 PRIOR LADDER: THE DISCOUNTED CURRENCY UNDER-PRICES A TRAINED OPERATOR (ledger L91; ORACLE DIAGNOSTIC)

`s26/PREREG_memorisation_on_the_ladder.md`; job `w_ladder` (40 s, 0.274 GB;
`s26/results/w_selfcopy_ladder.json`, 504 cells). The four leaked models per target travel
gam_prob 0.281 at cos_prob 0.584 in the S24 probability currency and gam_loc 0.773 at cos_loc
0.911 in the S25 location currency (E[d] MAE 2.34 -> 0.82 A per pair); the measured gain is
-0.709 A on the cloud (L44). The direction-discounted ladder (-2.1496 x gam x cos, the S25 L12
reading) predicts -0.364: disagreement +0.346 per target, 1.5x MDE, fold CI [+0.206, +0.472],
5/5 folds (the registered falsifier fired in the registered direction); the undiscounted
-2.1496 x gam_prob predicts -0.603, inside the MDE. Across the 504 cells corr(prediction, gain)
= +0.41 (S25 L12's re-readings: +0.05) and corr(MAE change, gain) = +0.72. Reading: for a
TRAINED posterior that moves far, the endpoint responds to the full probability-space move, not
to its projection onto the native's direction; the L12 discount was derived from small
re-readings and does not transfer. Consequence for Proposal C's arithmetic: quote -2.15 x
gam_eff (probability space) with the cosine reported beside it, not multiplied in. Not
deployable by construction; a calibration of the currency.

---

## 2g. THE SECOND ORPHAN: AMBER AS A DISTRIBUTION INSIDE THE PRIOR (tournament item 10; ledger L105). REFUTED IN THE STRONGEST FORM THE DESIGN ALLOWS.

`s26/PREREG_amber_prior_partner.md` (+ addendum 1, the staging); jobs `w_amberprior_probe` (5 s),
`w_amberprior_clouds` (135 s, 0.163 GB; 126/126, lam = 0 bit-exact against the shipped risk
table, beta = 0 equal to `p_ladder.pool_histogram`), `w_amberprior_endpoint` (4786 s under load,
0.131 GB). The leave-fold-out choice over 30 (lam, beta) cells on the point cloud is (0, 0) on
every fold: the deployable arm is the incumbent bit-exactly and every registered contrast is an
exact tie (126T). Every mixture cell is worse than the shipped posterior on the cloud, +0.010 A at
lam 0.05 to +0.152 at lam 0.5, monotone in lam (the typicality direction), none past its MDE, 22
of 29 with the fold CI above zero; beta moves a cell by at most 0.02 A with no consistent sign,
so AMBER's ordering carries nothing through the prior that the unweighted histogram does not
(S24 L16's rank-permuted null seen from the prior side). The per-target oracle over the 30 cells
(-0.29 A) is 224% accounted for by its valid across-target null and transfers 8% in split half.
The last AMBER form on the record is closed on this instrument.

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
   order 0.02 A on the built chain. MEASURED at n = 126 (L64, section 2b): m_tie 0.0039 A,
   paired MDE between two draws 0.0236 A (built chain), 3.6 of 75 members replaced per draw;
   every hundredths-level effect on the record is inside the paired MDE; the production
   convention is not distinguishable from a random draw.
3. `s26/IDEA_conformational_identity_floor.md`: same sequence, different deposit, 18 targets
   / 22 partners; the record had four numbers (S24 L4) and no distribution. MEASURED (L52,
   tournament item 2; section 1.7 and `s26/results/w_identity_floor.json`): median 2.908 A over
   22 pairs (3.055 over 18 targets); the copy is +0.97 A worse than the pool's best window
   (fold CI [+0.72, +1.28], 1.06x MDE, Type-M zone) and -1.44 A better than its mean (fold CI
   [-1.72, -1.19], 1.14x MDE, Type-M zone), n = 18, ORACLE on both sides. F5 holds.
4. `s26/IDEA_window_provenance.md`: whole-peptide / terminal / interior / fragment windows in
   the top-75; census first; plausibility 0.15. Ranked 9 (L51); runs as capacity allows, after
   the orphans.
5. `s26/IDEA_partial_recall_gradient.md` (found on the way to L44): does the fold model's
   memorisation (0.70 A at identity 1.0) have a gradient below the 0.6 threshold? Native-free
   covariates (max longer- and shorter-normalised identity to the training corpus, longest
   shared k-mer) against `rmsd_arm`; MDE 0.25 in rho at n = 126; registered expectation: no
   gradient (a clean bill for the fold threshold with a power statement). Cites S22 L3/L7 and
   S23 L3/L7: every router feature so far was a functional of the pool, the posterior or the
   emission, never of the training corpus relative to the query.
6. `s26/IDEA_memorisation_on_the_ladder.md` (found on the way to L44): the own-native models
   are the first real TRAINED operator far from the origin of the S24 prior ladder; measure
   their gam_eff and cos and test whether the ladder's transfer function (-2.15 A per unit
   gamma, with the S25 L12 direction discount) predicts the measured -0.70 A. Registered
   expectation: it under-prices trained priors by 2 to 3 MDE. ORACLE diagnostic, 5 minutes.

Also filed for the record: the medoid-frame and projection-branch flips (section 1.6) are
the mechanism the tie-break floor would measure; the L44 control population already shows the
same one-member change moving the chain 0.5 A on 5% of targets.

Orphan taken (L51 item 6): lane P's `IDEA_window_ensembling.md`. `s26/PREREG_window_ensembling.md`
states how fixed-K ensembling differs from widening K (S17 L12: one shortlist widened, the near-
native members displaced; here three shortlists of K = 500 each cut to their own top-75, no
shortlist widened, three clouds averaged in the production frame) and why the record predicts
null-to-worse anyway (S23 L9: 68% common-mode error shared by all three clouds; S24 L3: score
selection makes bias parallel across sources, cosine 0.943 against a within-source 0.933; S12:
the m-ladder is flat and m = 75 its argmin). Arms: shipped, b45, b80, ens3 (PRIMARY), ensK,
union3, boot3 (the zero-information ensembling control: three resampled clouds of the production
top-75); falsifier: ens3 minus shipped beyond its MDE with the fold CI excluding zero on 5/5 folds
AND beating boot3; expected 0.00 to +0.03 A against an MDE of about 0.05. Code `s26/w_ensemble.py`
(reuses `w_selfcopy.emit`; BLOSUM45/80 from Biopython re-ordered to `core.data.ALPHABET`, the
BLOSUM62 sums asserted equal to the universe's `sim` and `top_k` equal to the pinned pool),
synthetic tests `s26/w_ensemble_test.py` ALL OK. Probe and run after the tie-break job.

Rejected after reading, with the closing entry: ensemble-spread floor (S7 / S12, 1.044 A,
does not predict the gap); iterated Procrustes / GPA frames (S23 L4, FINDINGS 4225); distogram
to MDS (S8-10); Rg prediction-vs-pool disagreement as a router (S21 L27/L28 signal, S23 L7 router
null); BLOSUM rank as an additive score term inside the pool (S17 audit: BLOSUM ordering has
skill, small and decaying; S8-6 bounds fixing retrieval at 0.016 A; and L44's census says the
rank-0 window is never the argmin, so the term would act on a candidate the score already
rejects).

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
   line (0.151 A at the worst single target, L55), so the pre-registered class is MINOR even
   though every direct measurement is 0.008 A or less. The rule was written to be conservative
   and it is; I report the class it gives and the reason the envelope is loose, not a softer
   class.
5. The triangle bound is rigorous and loose by an order of magnitude on this operator (p95
   bound 1.44 A against p95 signed 0.15 A on the controls): the discrete flips move the chain
   orthogonally to the native. A native-free bound that never reads the native cannot know
   the direction of the move; it was the right tool before the gate and the wrong headline
   after it.
6. Part B's chain died with the host's memory kill after one model. The one model it made was
   the near-native case (1CEK) by luck of the pre-registered order; when the other three were
   built one at a time, the largest effect was on 2P5H, and it was HARMFUL to the leaked target
   (-0.246 A, F3 falsified in the direction nobody had registered). I had written F3 as a bound
   on |delta| as if the leak could only help; a memorising model with a copy that sits 2.3 A
   from the native pulls the answer toward the copy.
7. The artefact and the ledger disagreed on the paired gain (L55 caveat 1): `report()` had no
   envelope row for the gain basis, so the JSON said IMMATERIAL where the ledger, from the gain
   row I had computed in-session, said MINOR. A number computed outside the artefact's own
   function is a number the artefact does not carry; fixed by `s26/w_bound_addendum.py`, which
   also had to be re-run after every Part B re-run so the fix is not lost.
8. I over-estimated the tie-break's reach: I predicted about 8 of the 75 averaged members
   would change per draw and measured 3.6, because the score's top-75 is drawn mostly from
   the non-tied prefix. The floor is smaller than I guessed on membership and exactly where I
   guessed on the mean (0.004 A against a registered 0.003 to 0.010).
9. I nearly wrote the floor as a verdict on C3's relaxation. It is not: C3 was measured as a
   paired contrast with the pool held fixed, and that design excludes the tie-break noise by
   construction. The floor speaks to contrasts across runs and instruments, and the entry says
   so in its own words before anyone reads the table the other way (the Adversary's L71 caveat
   is exactly that the title needs the body's qualifier).
10. The recall covariate I proposed cannot separate memorisation from retrieval, because the
   fold model's training corpus is the retrieval library. I named S22/S23's router failures in
   the idea file and missed the confound in my own design; the PREREG addendum says how to fix it
   (partial out a retrieval-only covariate, run the mechanism check unconditionally).
11. I registered that the ladder would under-price the trained operator, and it did, but the
   mechanism is not the one I imagined: the leaked models recall the native's DISTANCES nearly on
   axis (gam_loc 0.77 at cos 0.91) while their 17-bin probability vectors move only 0.28 at cos
   0.58; the two currencies disagree about the same operator, and the endpoint follows the
   location. I had assumed a single "how far toward the native" number exists.
12. The ensembling orphan taught me nothing new and that was the point: I wrote in the PREREG
   that the three clouds would sit within 0.3 A of each other and that averaging them would do
   what a resample of the production set does, and both held to the second decimal (0.14 to 0.20
   A apart; 0.185 against 0.249 A of chain movement; -0.0005 against +0.0042 A of accuracy). A
   registered null that lands where it was registered is worth the 75 minutes because the
   direction is mandatory and is now closed with a power statement instead of an argument.
13. The provenance census surprised me in the direction I did not expect: I thought whole and
   terminal peptide windows would be the good ones; in the pool ANY peptide-derived window is
   0.42 A nearer the native than a fragment window, and the position in the parent does not
   matter (0.13x MDE). The score then removes three quarters of that gap, and the achievable
   readout cannot convert the rest. The interesting fact is about the corpus, not about the
   readout.
14. The amber prior partner refuted itself before its own controls could act: the
   leave-fold-out choice took lam = 0 on every fold, so the permutation seeds were never
   consumed and the "AMBER versus the unweighted histogram" contrast is an exact tie. I had
   planned the controls for a signal that the choice rule never produced; the cell table is the
   result, and it says the mixture is monotonically worse in lam with beta doing nothing.
15. The box, not the science, set the night's pace: the readout took 9,321 s of wall for
   2,016 projections under a six-job load, and one Part B waiter sat at the door for 23 minutes
   before the coordinator's starvation ruling withdrew it. Every job checkpointed and none was
   lost, which is what the contract's rules were for.

---

## 5. WHAT I DID NOT DO, AND WHY

- Did not read any benchmark sequence, name, PDB, native or RMSD, and did not parse
  `results/benchmark_manifest.json`, `s9/final_report.json` or any file under `results/`
  naming a benchmark target (L29). The benchmark facts used are the record's 2/60 and the
  mechanism in S24 L4.
- Did not compute the triangle bound on the benchmark itself, although it needs no native
  and no RMSD, because it needs the two benchmark sequences and their universes, which
  L18b/L29 close for this sprint. Named in L30 as a coordinator's option only.
- Did not finish Part B's control clause: the four carrier-out models are built and measured
  (L108); of the six control-out models one is built (9BAF out of fold 0) and the other five
  wait behind the coordinator's starvation ruling (L111) and the sprint's close. F3's control
  clause is therefore open, and the entry says so. A second-seed retrain of 2P5J (the one
  large value) was not run for the same reason.
- Did not use length-matched non-dev control chains for Part B, because naming one risks
  naming a benchmark peptide; dev-target 16-mers are the controls, and the length mismatch is
  stated.
- Did not edit `s26/w_selfcopy.py` after its first jobs started; every later step runs through
  a driver (`w_endpoint_report.py`, `w_bound_addendum.py`, `w_tiebreak.py`, `w_ensemble.py`,
  `w_provenance*.py`, `w_amberprior.py`, `w_recall.py`, `w_ladder.py`) that imports it.
- Did not quote the envelope's 0.028 A as the leak's size: it is an upper bound from an
  operator stronger than the leak (own-native training), stated as such beside the direct
  0.008 A and, since L55, beside its worst-target reading 0.151 A under A2.
- Did not project the non-chosen cells of the amber prior partner on the built chain (15 h at
  the loaded box's rate); the declared staging chooses on the cloud and projects the chosen
  cell only, and since the chosen cell is the identity on every fold, no built-chain number
  exists for any mixture and none is claimed.
- Did not run the tie-break floor at 16 draws or on the AMBER-relaxed basis, the ensembling
  control at a second seed, or the provenance readout with 8 permutation draws: each is the
  pre-declared replication or extension for a positive, and none of the three was positive.
- Did not run the recall gradient's ORACLE mechanism check (the posterior's MAE against
  I_short): the PREREG gated it on a gradient at the MDE, which did not exist; the addendum
  says it should run unconditionally next time.
- Did not use `docs/REPORT_S26.md` or `docs/REPORT_S26_SUMMARY.md` (L42: not written by any
  lane, untracked) for any number.
- Did not write the report or the deck: lanes E and PR own them; every number this lane
  contributes is in the ledger entries L30, L44, L52, L58, L64, L84, L85, L90, L91, L105, L108
  and L109 with its artefact path, and in this file.

---

## 6. ARTEFACT INDEX

Pre-registrations (each with its addenda): `s26/PREREG_selfcopy_bound.md` (1 to 3),
`PREREG_tiebreak_floor.md` (1, 2), `PREREG_identity_floor.md`, `PREREG_window_ensembling.md` (1),
`PREREG_window_provenance.md` (1, 2), `PREREG_amber_prior_partner.md` (1, 2),
`PREREG_partial_recall_gradient.md` (1), `PREREG_memorisation_on_the_ladder.md` (1).
Ideas: `s26/IDEA_selfcopy_proxy_bound.md`, `IDEA_tiebreak_noise_floor.md`,
`IDEA_conformational_identity_floor.md`, `IDEA_window_provenance.md`,
`IDEA_partial_recall_gradient.md`, `IDEA_memorisation_on_the_ladder.md`.
Code: `s26/w_selfcopy.py` (+ `w_selfcopy_test.py`), `w_train_chain.py`, `w_endpoint_report.py`,
`w_bound_addendum.py`, `w_identity_floor_stats.py`, `w_tiebreak.py` (+ `w_tiebreak_test.py`),
`w_tiebreak_report.py`, `w_ensemble.py` (+ `w_ensemble_test.py`), `w_provenance.py`
(+ `w_provenance_test.py`), `w_provenance_readout.py`, `w_amberprior.py`, `w_recall.py`,
`w_ladder.py`.
Results (`s26/results/`): `w_selfcopy_{census,retrieval,envelope,posterior,floor,endpoint,bound}.json`,
`w_identity_floor.json`, `w_selfcopy_tiebreak_{draws,endpoint}.json`, `w_tiebreak_report.json`,
`w_selfcopy_ensemble_{clouds,endpoint}.json`, `w_selfcopy_provenance_{census,oracle,readout}.json`,
`w_selfcopy_amberprior_{clouds,endpoint}.json`, `w_amberprior_cells_cloud.json`,
`w_selfcopy_recall_{covariates,endpoint}.json`, `w_selfcopy_ladder.json`, and the one-target
probes `w_selfcopy_{retrieval,envelope,tiebreak,ensemble,amberprior}_probe_*.json`.
Models (untracked, `*.pt` under `s26/models/w_selfcopy/`): `pca32_fold2_s0_out_1A11.pt`,
`pca32_fold4_s0_out_2LMF.pt`, `pca32_fold4_s0_out_2P5J.pt`, `pca32_fold0_s0_out_1U6V.pt`,
`pca32_fold0_s0_out_9BAF.pt`.
Job records and logs: `s26/jobs_done/w_*.json`, `s26/logs/w_*.log`.
Ledger (lane W): L30, L44, L52, L58, L64, L84, L85, L90, L91, L105, L108, L109. STATUS lines
under `## W (Wildcard)`. Closing entry: L137.
