# PREREG A3 -- DOES A TARGET-DEPENDENT HAMILTONIAN GIVE THE SELECTOR SOMETHING TO DO? (lane Q)

Filed 2026-09-13 09:00, before any A3 arm has been built on more than one target. Not edited
after; addenda are appended. Two halves: a PROPERTY half (distinct trained states, operator
sequences; no native) and an ENDPOINT half (gated on "PHASE 0 SIGNED OFF").

## 1. The question (S25 L17, open)

`E = _zrank(sorted scores)` is the same vector on every target up to tie-averaging, so the
deployed selector solves the same variational problem 126 times and there are two trained
states in the whole deployment. Any Hamiltonian that varies between targets while preserving
the CVaR selection semantics (the tail stays a prefix of the score order) gives the state
target-specific structure to represent. Does it, and does the endpoint notice?

## 2. The variants (all strictly increasing in the score, so the classical selection is unchanged)

    zrank   the deployed E (standardised ranks, ties averaged)
    zraw    (s - mean)/sd of the raw Bayes-risk score over the 128: affine, gaps survive
    asinh   S21 L33's `Nt`: asinh of a robust (median/MAD) z-score, re-standardised
    soft    S13's rank-preserving soft compression of the upper tail (geo_common.soft_compress),
            re-standardised
    zraw_Tmatch  zraw at a per-target temperature T' chosen native-free so the Gibbs entropy
            at zraw equals the Gibbs entropy at zrank, T = 0.3 (4.91 bits): separates
            "target-dependent gaps" from "sharper weights"

Monotonicity is asserted at build time (`energy_variants`); S25 L16's finding that a moment
z-score of raw AMBER is not monotone in float64 does not apply to the distogram score (values
of order 1 A, no 1e28 outliers), and the assertion would fire if it did. The squared-risk
functional (posterior mean instead of median) changes the ORDER and is not a variant here; it
is listed in `s26/IDEA_l17_target_dependent_hamiltonian.md` as a different class.

## 3. Arms

    fixed_{ev}_it50, fixed_{ev}_it750       ev in {zrank, zraw, asinh, soft}
    fixed_{ev}_Tmatch_it50                  ev in {zraw, asinh, soft}
    adaptL2_adam_best_{ev}_P7/14/21         ev in {zrank, zraw}
    gibbs_T, uniform128, randH_fixed, randH_adaptL2   (as in A1)

## 4. Hypotheses and falsifiers

H3a (property). Distinct trained states across the 126 targets, per `(alpha, T)` cell, counted
greedily in pdb order with "distinct" = symmetric KL > 0.01 nats (`distinct_states`):
   under zrank, fixed_zrank_it50: at most 3 per cell (2 cells) -- the tie-averaging outliers;
   under zraw, fixed_zraw_it50: at least 100 of 126.
   Falsified if zrank gives more than 10 per cell or zraw gives fewer than 60.
H3b (property). ADAPT-selected operator sequences (pool L2, adam_best): distinct sequences per
   cell, zrank: at most 2 per cell; zraw: more than 40 of the cell's targets.
   Falsified if zrank exceeds 10 per cell or zraw is below 20.
H3c (property). The Gibbs state at zraw, T = 0.3 is sharper than at zrank: mean entropy below
   3 bits against 4.91 (1A13: 1.84 bits, `s26/results/probe/1A13.json`).
H3d (endpoint, secondary, matched readout). `fixed_zraw_it50 - fixed_zrank_it50` on the built
   chain. PREDICTION: POSITIVE (worse), because S25's entropy curve (`s25/results/q_alpha.json`,
   `entropy_curve`: RMSD 3.454 at 0 bits, 3.404 at 3.3 bits, 3.28 near 5.7 to 6.3 bits) puts a
   sharper state above the deployed one, and S23 L8 found sharper weighting harmful at T = 0.
   Expected +0.03 to +0.10 A on the selection basis, smaller on the built chain.
   `fixed_zraw_Tmatch_it50 - fixed_zrank_it50`: PREDICTION null (within 0.5x MDE): once the
   entropy is matched the target-dependent gaps buy nothing the readout can see.
   "The target-dependent Hamiltonian helps" fires only if a Tmatch contrast is negative past
   its MDE with the fold CI excluding zero, 5/5 folds, and replicates on seed 1 and reverse
   order. "It is null" fires within +-0.5x MDE. "It is harmful" fires if positive past MDE with
   the fold CI excluding zero (a result: the entropy mechanism, not a defect).

## 5. Expected effect against the MDE

Reference MDEs as in PREREG_A1 section 7 (`s26/results/q_mde_reference.json`): built chain
0.0958 A (weighted vs uniform), selection 0.1261 to 0.1652 A. The zraw-vs-zrank paired SD may
exceed these (the states differ by design), so the run's own MDE is what counts.

## 6. Memory and time

Same harness as A1 (`s26/q_adapt.py --build --tag a3 --variants zrank,zraw,asinh,soft
--adapt-variants zrank,zraw --pools L2 --optimisers adam_best`): 8 fixed arms + 3 Tmatch + 6
ADAPT snapshots + 4 controls = 21 projections, about 95 s per target, 126 targets about 3.3 h
on one process; peak RSS under 0.4 GB (`s26/jobs_done/q_probe_1A13.json`). Agent-hours: 4 to 6.

## 7. Rule 0, six forks

1. FUNCTIONAL: paired CA-RMSD of the emitted structure. NOT TAKEN: KL, entropy, F.
2. BASIS: built chain via the production projection; selection arm beside it. NOT TAKEN:
   point cloud; AMBER emission.
3. READOUT: as production (`average_weighted`, `project`, `consensus_medoid`). NOT TAKEN: a
   readout that re-ranks by the variant (the variant is monotone; the classical order is the
   same on every arm, which is the point).
4. NORMALISATION: the four monotone variants above, each re-standardised to mean 0, sd 1 over
   the 128 so the temperature keeps its meaning (S25 L16: the outer re-standardisation is
   load-bearing for T), plus the entropy-matched temperature. NOT TAKEN: a variant that
   changes the order (squared risk); winsorisation (S13: 99th-percentile winsorisation is not
   sufficient conditioning and is not monotone at the cap).
5. NULL: fixed_zrank_it50 (the deployed state), gibbs_T, uniform128, randH. NOT TAKEN:
   initialisation means.
6. LABEL: `s24.stats_lib.compare`, fold-clustered, MDE, W/L, median, concentration. NOT TAKEN:
   marginal means.

## 8. What I will not do

No variant is selected by its endpoint. If `zraw_Tmatch` is null and `zraw` is harmful, the
L17 answer is "the spectrum can be made target-dependent; the readout responds to its entropy
and to nothing else that was varied", and that is what is reported.
