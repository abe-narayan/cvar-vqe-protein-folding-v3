# IDEA_selfcopy_proxy_bound -- BOUND THE 2/60 BENCHMARK SELF-COPY LEAK FROM THE DEV-SET PROXY, WITHOUT OPENING THE BENCHMARK (lane W, Sprint 26; mandatory direction, L25)

## Hypothesis

The leak that reaches the sealed benchmark (two 11-mers carried verbatim by one longer database
peptide in another pinned fold, S24 L4; the count 2/60 re-derived in memory by lane I, L15/L18)
enters the pipeline through exactly two channels: (A) the carrier's self-window sits in the
K = 500 retrieval pool at BLOSUM rank 0, and (B) the carrier's native distances were training
labels for the target's own fold model. Both channels exist, with the same direction (target
carried, carrier in another fold), on four dev targets: 1CEK in 1A11 (folds 2 vs 3), 2FBU in
2LMF (4 vs 0), 2P5H in 2P5J (4 vs 2), 6B9K in 1U6V (0 vs 2). The claim: the benchmark leak's
contribution to every benchmark figure, including the paired gain +0.0103 [-0.160, +0.180], can
be bounded from those four dev targets plus an n = 126 upper envelope, using no benchmark
sequence, PDB, native or RMSD, and the bound is immaterial (below 0.017 A, one tenth of the
benchmark CI half-width).

## Why the record does not already close it

- S10-4 priced channel A on dev at +0.0004 A (13 targets with a >= 0.6 window, synthesis `fit`
  basis 3.2005) and S9-10 priced the same operator on the benchmark at +0.0030 [-0.0002,
  +0.0061] over 13 targets. Neither separated the two exact copies from the eleven partial
  ones, neither was measured on the production built chain, and lane E's examination
  (L28, claim C27) finds the S10 artefacts absent from disk: the +0.0004 / +0.0030 are
  document-only numbers the Adversary may strike.
- Channel B was never priced anywhere: S10-4 dropped windows and kept the distogram fixed.
- S24 L4 declared the 2/60 and chose to leave it unquantified because quantifying it directly
  needs benchmark RMSD; the proxy route was named in the S26 campaign prompt and is unclaimed
  (L25).
- A native-free census this turn (`s26/w_selfcopy.py census`, before any RMSD): the four
  self-windows are BLOSUM rank 0 and pool position 0 on all four, but only 2P5H's and 6B9K's
  are inside the shipped top-75 (`sub` of the production cache); on 1CEK and 2FBU channel A
  cannot touch the emitted structure at all.

## Exact falsifier (full protocol in `s26/PREREG_selfcopy_bound.md`)

The bound is computed two ways and both must be reported: (i) realistic, from the four dev
targets with both channels removed (self-window dropped and pool refilled; fold model retrained
without the carrier), per-target signed effect on `sel`, `cloud`, `arm`, `fit` and on the
paired gain `arm - sel`, plus a native-free triangle-inequality bound (RMSD between the leaked
and un-leaked emissions, which bounds the change in RMSD-to-native without reading the native);
(ii) envelope, at n = 126: the built chain from a fold model that had the target's OWN native
among its training labels (the four other pinned fold models) against the clean model, paired,
fold-clustered. The benchmark bound is (2/60) x max per-target effect of (i) and (2/60) x the
fold-CI limit of (ii). The idea is FALSIFIED as "immaterial" if either bound exceeds 0.017 A
on any production basis; it is REFUTED as a proxy if any dev self-copy target moves by more
than 0.10 A on the built chain when both channels are removed (then the S10-4 price was not
representative and the caveat on every benchmark figure is upgraded from "declared" to
"material, bounded at X").

## Expected effect against the computed MDE

Channel A on the built chain: |effect| < 0.10 A per target (S10-4's +0.0041 mean on 13 leaked
targets; two of four have zero possible effect). Channel B: unknown, the record's only prior is
`core/predict.py`'s overfit note (train NLL 0.944 against validation 2.025), which says the MLP
fits its training pairs far better than held-out pairs, so the n = 126 envelope is expected to
be a measurable NEGATIVE effect (the leaked model helps) of order -0.2 to -0.8 A on the built
chain, against an MDE of 0.09 to 0.23 A for a distogram swap (PREREG_C2 section 4). Even at
-0.8 A the benchmark envelope is (2/60) x 0.8 = 0.027 A, and the realistic bound is expected
an order of magnitude smaller because a carrier's segment is not the native (S24 L4: three of
four self-windows sit 1.0 to 2.0 A WORSE than the pool's own best).

## Memory and agent-hours (probe before every run)

Retrieval and envelope parts: one universe at a time (1.3 MB npz, 15 MB resident), the s12
distogram cache and the five pinned fold models; < 0.8 GB expected, 126 x 7 projections at
~3 s = 45 min CPU. Retraining (channel B): 1.248 GB peak and 456 s per fold model measured by
lane P (`s26/jobs_done/p_probe_train_pca32_f0.json`); 4 carrier-out plus 6 control-out models
= 10 jobs, 1.3 h one at a time; the reference models are lane P's `pca32` rung (identical code
path) once its chain reaches them. Agent-hours: 3 code, 2 runs, 2 write-up.

## Information value

Closes open item 2 of the state brief with a number and its assumptions, re-derives claim C27
on the production basis so it can be sourced, and measures for the first time what a
verbatim training copy is worth to the shipped MLP through the whole pipeline (the envelope),
which is also the quantity behind the fold-model non-independence caveat (EXAMINATION E4).

---
Measured (2026-09-13, L44; Adversary L55 STANDS WITH CAVEAT): every direct dev-proxy measurement
is IMMATERIAL (0.002 A or less); the own-native envelope gives the pre-registered class MINOR,
mean-CI limit 0.028 A on the built chain (0.048 selection, 0.023 paired gain), worst single target
0.151 A (2BP4) under assumption A2; MINOR under every reading on the built chain and the paired
gain. `s26/results/w_selfcopy_bound.json`.
