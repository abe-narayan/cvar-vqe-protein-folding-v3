# S27 REPORT: ALTERNATIVE HAMILTONIANS FOR THE CANDIDATE-POOL SELECTOR

Sprint 27, 2026-09-14. Pre-registration `s27/PREREG.md` (addenda 1 and 2 written before the
runs they govern). Code `s27/ham_lib.py` (the channels), `s27/run_pool.py`, `run_wave2.py`,
`run_wave3.py`, `run_vqe_chain.py`, `run_trainability.py`; statistics `analyse_pool.py`,
`analyse_arms.py`; tables `make_tables.py` -> `s27/results/tables.md` (reproduced in full as
Appendix A of this file). Every number below is read from `s27/results/*.json`. Nothing in the
production path, the benchmark or any earlier sprint's file was touched.

## 1. The question and the answer

Can a different Hamiltonian, energy model, or principled combination of Hamiltonians give a
meaningfully better final RMSD than the shipped distogram Bayes-risk score, on the established
126-target instrument, with the production readout and the genuine CVaR-VQE selector intact?

**No.** Across 203 configurations on the point cloud (18 new native-free channels in six
families, the eleven Legacy terms one at a time, the three S25 references, 31 equal-weight
pairs with the distogram and their rank-permuted controls, 62 weighted pairs, 32 staged
rejects, 4 composites alone and with the distogram, 7 adaptive mixtures), 47 wave-2 and 9
wave-3 configurations (the m-ladder, in-band re-ranking, non-additive forms, near-corpus
potentials, secondary-structure variants), 80 configurations inside the genuine CVaR-VQE at two
seeds, and 10 on the built chain:

- no configuration beats the shipped score beyond its own minimum detectable effect on any
  endpoint; the best point-cloud effect is -0.026 A at 0.35x MDE, and it vanishes on the built
  chain (+0.0003 A);
- 45 configurations are measurably WORSE (fold-clustered CI above zero and beyond their MDE);
  every new channel used alone is worse than the shipped score, and most are worse than a
  random 75-subset of the same pool;
- the genuine CVaR-VQE arm tracks the classical top-m of every energy it is given (set-equality
  gate 126/126 on all 80 configurations, both seeds), so the selector cannot rescue a Hamiltonian
  the classical ranking does not already favour, and its trainability is the same for every
  Hamiltonian (all within 0.86x to 1.14x of the distogram's gradient variance).

The mechanism that closes the search is measurable and is the sprint's one real finding
(section 6): on this readout, a channel's correct ranking information is anti-useful. The more
a channel knows about which single candidates are near the native beyond what the distogram
knows (its partial rank correlation with the ORACLE candidate RMSD given the distogram), the
more it HURTS the averaged top-75 (Spearman +0.48 across the 18 channels), because such
channels select similar candidates and the readout averages them. The consistency family has
the most information (partial rho +0.25 for CONS) and does the most damage (+0.29 A when added
to the distogram, worse than its own permuted control by +0.27 A at 1.7x MDE).

## 2. What was held fixed

Identical to the S25 seven-configuration suite (`s25/PREREG_PHYS.md`): the 126 dev targets, the
5 pinned folds, the shipped K = 500 BLOSUM pool, the uniform top-75 coordinate average in the
retained set's medoid frame, `s24.stats_lib.compare` (paired, MDE = 2.8016 x SE, fold-clustered
CI beside iid, W/L, concentration null), exact ties broken by a stable per-target random key
(never by array order). Point cloud primary (the shipped score's top-75 reproduces the S25 anchor
3.048338 A on all 126 targets, `s27/results/pool_summary.json :: dis_mean`); built chain through
the production projection for the finalists (the shipped score reproduces the rebuild-basis
anchor 3.2126, `arms_summary.json :: chain/DIS`); the genuine CVaR-VQE of `s24.d_harness.arm_vqe`
(9 qubits, 3 layers, 80 Adam iterations on the exact parameter-shift gradient, alpha 0.18, T 0.5)
reproduces S25's distogram arm 3.0580 (`arms_summary.json :: vqe/DIS/0/mean_vqe`).

Statistical potentials are fitted per target on that target's own leakage-safe window universe
(the peptides of the other four folds plus the fold's identity-filtered fragments, `s8/
generate.py :: stage_univ`), the same rule the shipped distogram trains under; no native enters
any channel; natives are read only in `oracle_*` functions for the endpoint and the labelled
diagnostics.

## 3. The channels (each is `s27/ham_lib.py :: h_<name>`; cost per evaluation in Appendix T2)

| family | channel | what it scores | literature motivation |
|---|---|---|---|
| geometric | RG_LAW | deviation of Rg from 2.2 n^0.38 | folded-protein Rg scaling, nu 0.34 to 0.38 |
| geometric | RG_UNIV | deviation of Rg from the universe median | data-driven compactness target |
| geometric | EXVOL | CA-CA excluded volume below 4.0 A, sep >= 3 | hard-core repulsion |
| geometric | CAGEO | -log p(virtual angle, virtual torsion) from the universe | Levitt-style CA-trace statistics |
| torsion | RAMA | -sum log p_fold(aa, phi, psi), per-fold table `s8/generate_rama.npz` | Ramachandran statistical potential |
| statistical | CONTACT | quasi-chemical contact potential, CA-CA < 7.5 A, fitted on the universe | Miyazawa-Jernigan / Sippl reference state |
| statistical | DISTPOT | distance-dependent pair potential, 1 A bins 3 to 14 A, 3 separation classes | Sippl / DFIRE form at CA level |
| statistical | ENV | one-body burial statistics -log p(neighbour bin \| aa)/p(bin) | Rosetta centroid `env` |
| statistical | HP | -sum KD(aa) x neighbour count | hydrophobic burial, Kyte-Doolittle |
| physics | DSSPHB | Kabsch-Sander H-bond energy on the ideal rebuild, sum of E < -0.5 | DSSP electrostatic H-bond model |
| physics | ELEC | Debye-Hueckel Coulomb between charged residues at CB, eps 80, kappa 0.1/A | screened electrostatics |
| consistency | CONS | mean CA-RMSD to the other pool members | S12 medoid criterion |
| consistency | DMAP_CONS | L1 distance to the pool's median distance map | distance-space consensus |
| consistency | TORS_CONS | circular distance to the pool's per-residue mean torsions | torsion consensus |
| consistency | POOLGO | Go-like reward for the pool's consensus contacts | Go model with the pool as the reference |
| compatibility | SS_MATCH | CA-geometry SS call vs Chou-Fasman propensity | sequence-structure compatibility |
| re-reading | CONTACT_LL | contact-map log-likelihood under the distogram's contact probabilities | coarse re-reading of the posterior |
| re-reading | DIS_MEAN | L1 error to the posterior mean instead of the median | q = 2 location statistic |
| reference | DIS, LEG, AMB, LEG_<11 terms> | the S25 references and Legacy's eleven terms alone | |

Wave 2 added NC_DISTPOT and NC_CONTACT (the same potentials fitted on the target's 2,000 most
BLOSUM-similar universe windows); wave 3 added SS_MATCH2 (the SS call from Kabsch-Sander
H-bonds on the ideal rebuild against helix/strand propensities fitted on the universe).

Sources consulted for the constructions: the Miyazawa-Jernigan contact potential and its
quasi-chemical reference state ([Statistical contact potentials in coarse-grained modeling](https://link.springer.com/chapter/10.1007/978-1-4419-6889-0_6),
[knowledge-based energy functions](https://arxiv.org/pdf/q-bio/0601026)); the Rosetta centroid
terms env, pair, cbeta, vdw, rg ([Rosetta centroid score terms](https://docs.rosettacommons.org/docs/latest/rosetta_basics/scoring/centroid-score-terms));
the Kabsch-Sander H-bond energy and hydrogen placement ([DSSP](https://en.wikipedia.org/wiki/DSSP_(algorithm)),
[mdtraj.kabsch_sander](https://mdtraj.org/1.9.4/api/generated/mdtraj.kabsch_sander.html)); the Rg
scaling exponent of folded proteins ([polymer scaling laws of unfolded and disordered proteins](https://www.pnas.org/doi/full/10.1073/pnas.1207719109)).

## 4. Results by hypothesis (every table in Appendix A)

**H1, singles (T1, T2).** Every channel used alone is worse than the shipped score; 31 of the 36
singles and composites are worse beyond their MDE with the fold CI above zero. The order of
merit among the alternatives: DIS_MEAN 3.078 (a re-reading, +0.030, 0.5x MDE), DISTPOT 3.199
(+0.151, 0.63x; the best of the new Hamiltonians, better than a random subset by -0.222 at
0.9x MDE and better than Legacy 3.754 and AMBER 3.881 by half an Angstrom), STAT_COMP 3.211,
ROSETTA_LIKE 3.350. Below that, every channel is at or worse than a random 75-subset (3.421):
the consistency family (CONS 3.688, TORS_CONS 3.827, DMAP_CONS 3.856, POOLGO 3.949), the
backbone physics (DSSPHB 3.784, ELEC 3.869), the torsion and geometry statistics (RAMA 3.929,
CAGEO 3.967), and the two genuine physics energies (LEG 3.754, AMB 3.881, reproducing S25).
Ranking skill does not predict single-arm value: Spearman between a channel's whole-pool rho
with the ORACLE RMSD and its top-75 effect is 0.04 across the 18 channels
(`s27/results/mechanism.json`).

**H2, equal-weight pairs (T3).** No pair beats the shipped score. Against the rank-permuted
control (X's marginal kept, its correspondence destroyed), only DIS+DISTPOT (-0.040, 0.35x MDE)
and DIS+ENV (-0.035, 0.28x) sit below their own noise, neither beyond its MDE; 11 pairs are
worse than their own noise with the fold CI above zero (DSSPHB, RAMA, TORS_CONS, CONS, CAGEO,
POOLGO, DMAP_CONS, LEG and its helix, local H-bond and torsion terms, RG_UNIV): adding those
channels to the distogram is worse than adding noise of the same amplitude, S25 L16's result
extended to eleven more channels.

**H3, regulariser weights (T4).** Nested leave-fold-out choice of the weight in {0.25, 0.5, 1.0}:
the only held-out effect with the fold CI below zero is SS_MATCH at w = 0.5 (-0.026 A, 0.35x
MDE, chosen on all five folds). Completing the grid to {0.25, 0.5, 0.75, 1.0} in wave 3 the
held-out effect falls to -0.008 A (0.1x MDE, fold CI [-0.044, +0.044]) and `best_of_k_within`
reads "NOT A SIGNAL (split-half transfers 10% of the oracle)" (T4, T7). The half-weight bump was
a grid artefact.

**H4, staged rejects (T1).** Rejecting the worst 10% by any channel before the distogram's
top-75 is null: the best, REJ[POOLGO]->DIS, is -0.016 A (0.37x MDE) and the random-reject
control is +0.007; rejecting by ELEC, CAGEO, LEG_contact, LEG_electrostatic is worse with the fold
CI above zero. S26 L43/L86 (the AMBER reject) generalises to every channel here.

**H5, composites (T1, T2).** The Rosetta-centroid-like composite (ENV + CONTACT + EXVOL +
RG_LAW) is 3.350 alone (+0.302, 0.99x MDE, worse) and +0.012 with the distogram; the physics
composite 3.552 alone and +0.026 with; the consistency composite 3.904 and +0.282; the
statistical composite 3.211 and +0.011. None is a candidate.

**H6, adaptive mixtures (T1).** Weighting the channel by the distogram's mean per-pair entropy
(fixed map, not fitted): +0.009 to +0.3 A for every channel tried. Null to worse.

**H9, the m-ladder (T5).** At no m from 3 to 150 does any arm beat the distogram at the same m
beyond its MDE. The pre-registered prediction that DIS+CONS would win at small m is FALSIFIED:
it is worse at every m (+0.159 at m = 3, +0.286 at m = 75). The distogram's own ladder has its
minimum at m = 75 (3.048; m = 50 3.068, m = 100 3.056), reproducing S22 L4. The per-arm minimum
over m, priced as an order statistic, transfers 19% (DIS+0.5*SS_MATCH) and 29% (DIS+DISTPOT) of
its oracle in split halves: not signals.

**H10, in-band re-ranking (T6).** Inside the distogram's top-150, no channel picks a better 75
than the distogram itself beyond its MDE: SS_MATCH -0.022 (0.25x), everything else null to
worse, and the consistency channels are the worst in-band (CONS +0.184 at 1.08x MDE): the pool's
mode is not the native even among the near-native band. S17 L23 and S25 §1.5 reproduced.

**H11 to H13 (T7).** Non-additive forms (max-rank, rank-product) are null; the near-corpus
potentials are null (DIS+NC_DISTPOT -0.026 at 0.29x); SS_MATCH2 (the better-founded
secondary-structure variant) is no better than SS_MATCH (-0.010 at w = 0.5, 0.16x MDE) and is
4.417 A alone.

**H7, the genuine CVaR-VQE arm (T9).** 80 configurations x 126 targets x 2 seeds. The
set-equality gate passes on every cell; for the 47 configurations whose top-75 has no exact
ties the residual RMSD_VQE - RMSD_top-m is at most 6.3e-4 A (the theorem reading back); for the
33 configurations with heavy ties (EXVOL and LEG_steric at 99% ties, LEG_coop_sheet 96%,
LEG_compactness 84%) the classical prefix is not unique and the residual measures tie
ambiguity, not a quantum effect. Against the distogram's own VQE arm (3.0580, S25's number),
the best configuration is DIS+0.5*SS_MATCH at -0.037 A (0.47x MDE, fold CI [-0.066, -0.010],
seed 1 -0.029, inside seed 0's CI) and DIS+0.25*LEG_hbond_longrange at -0.027 (0.69x); neither
clears its MDE, and 51 configurations are worse with the fold CI above zero. Realised tails are
70 to 75 for every non-degenerate energy (entropy 8.8 bits of 9), 93 with 9.00 bits for the
99%-tied channels (the uniform state: an all-ties Hamiltonian is empty, S25 §2's raw-moment
AMBER seen again). The VQE does not change the answer for any Hamiltonian: the selector's
verdict is the classical ordering's verdict, up to the m-ladder.

**Built chain (T10).** The production projection removes the last of the point-cloud
suggestions: DIS+0.5*SS_MATCH +0.0003 A (0.00x), DIS+DISTPOT -0.003 (0.02x), REJ[POOLGO]->DIS
-0.017 (0.34x, fold CI [-0.037, +0.005]); DISTPOT alone +0.171, DIS+CONS +0.249, CONS +0.553,
all worse with the fold CI above zero.

**H8, trainability (T11).** Under the deployed rank currency every channel's gradient variance
at n = 9, alpha 0.18, T 0.5 lies within 0.87x to 1.09x of the distogram's (the ladder is the
same object up to ties); under the bounded asinh-MAD standardisation that keeps each channel's
own gaps, 0.86x (AMB) to 1.14x (CONTACT). The falsifier (a 2x difference) does not fire: no
Hamiltonian in this family gives the circuit a measurably easier or harder objective.

## 5. Regimes (per length and FAIL18, `s27/results/strata.json`)

Every alternative that is not harmful overall is helpful on the 18 targets where the distogram
is wrong and harmful or null on the 108 where it is right: DIS+0.5*SS_MATCH -0.234 on FAIL18
and +0.009 elsewhere; DIS+DISTPOT -0.191 and +0.010; DIS+ENV -0.401 and +0.064; DISTPOT alone
-0.364 and +0.237; even CONS, LEG and AMB are -0.13 to -0.23 on FAIL18 while costing +0.77 to
+1.00 elsewhere. The net is null because the switch is not available native-free: seven router
constructions across S22, S23 and S26 (L115) fail to identify FAIL18 from the pool, and S26's
adaptive gate here (H6) is null. By length, the top combinations are neutral below 12 residues
and -0.06 to -0.19 at 13 to 14, an interaction inside the noise.

## 6. Mechanism: ranking information is anti-useful on an averaging readout

`s27/results/redundancy.json` and `mechanism.json`. For each channel X: the partial Spearman
correlation of X with the ORACLE candidate RMSD after removing the distogram's rank (what X
knows that the distogram does not), and the effect of DIS+X on the top-75 average.

| channel | partial rho given DIS | DIS+X effect (A) |
|---|---:|---:|
| CONS | +0.246 | +0.286 |
| TORS_CONS | +0.122 | +0.263 |
| DMAP_CONS | +0.119 | +0.356 |
| POOLGO | +0.113 | +0.299 |
| DISTPOT | +0.090 | -0.018 |
| CONTACT_LL | +0.058 | +0.088 |
| CAGEO | +0.048 | +0.290 |
| ENV | +0.038 | -0.002 |
| SS_MATCH | +0.017 | +0.022 |
| AMB | -0.018 | +0.081 |

Across the 18 channels, Spearman(partial rho, DIS+X effect) = +0.48 and Spearman(rho with DIS,
DIS+X effect) = +0.52: the channels that carry the most correct ranking information, and the
channels most redundant with the distogram, are the ones that hurt most when added. The reading
is the S12 set-mean law and S23 L5 in one measurement: the top-75 average is a variance-reduction
operator; a channel that finds near-native members finds SIMILAR members (CONS, by definition,
picks the pool's mode), the retained set loses the diversity the averaging needs, and its mean
moves toward the pool's mode rather than the native. The one family that does not hurt in
combination (DISTPOT, ENV: partial rho 0.04 to 0.09, rho with DIS 0.27 to 0.43) is the one that
adds a little information without collapsing diversity, and it adds too little to measure. This
is why a better ranker is not a better selector on this readout, and it is a property of the
readout, not of any Hamiltonian.

## 7. Verdicts

| category | configurations |
|---|---|
| PROMISING (beyond MDE, fold CI excluding zero, built chain agreeing) | none |
| NULL, underpowered (inside 0.7x MDE; the effect a real lever would show, 0.05 A, is below the MDE of every pair with a real channel) | DIS+DISTPOT, DIS+ENV, DIS+NC_DISTPOT, REJ[POOLGO]->DIS, REJ[LEG_compactness]->DIS, REJ[RG_LAW]->DIS, DIS+0.25*LEG_hbond_longrange, the SS_MATCH family (its 0.5-weight effect is a grid artefact, T4/T7), max(DIS,ENV), rankprod(DIS,DISTPOT), DIS+STAT_COMP, DIS+ROSETTA_LIKE, the adaptive mixtures |
| WORSE (fold CI above zero, beyond MDE) | 45 configurations: every consistency channel alone and with the distogram; DSSPHB, RAMA, CAGEO, TORS_CONS, LEG_coop_helix, LEG_hbond_local with the distogram; every single except DIS_MEAN, DISTPOT, STAT_COMP, ROSETTA_LIKE, LEG_coop_sheet, ENV, LEG_hbond_longrange, EXVOL, LEG_steric (which are worse in the Type-M zone or below) |
| ORACLE-ONLY | the per-arm minimum over m (T5; 19 to 29% split-half transfer); the partial correlations of section 6 |
| MECHANISM SUPPORTED | ranking information reduces retained-set diversity and hurts the average (section 6, 18 channels) |
| MECHANISM UNSUPPORTED | "a target-dependent Hamiltonian changes what the CVaR-VQE can select" (the VQE tracks the classical prefix on all 80 x 2 x 126 cells); "a knowledge-based potential fitted on the retrieval corpus adds corpus information the distogram lacks" (DISTPOT is the best new single and still null in combination) |

## 8. The ranked table (the ten best and the references; the full 203 rows are T1)

| rank | configuration | mean RMSD (A) | vs DIS | x MDE | fold CI | W/L | cost | useful alone? | useful in combination? | inside the CVaR-VQE (vs DIS-VQE) |
|---:|---|---:|---:|---:|---|---|---|---|---|---:|
| 1 | DIS+0.5*SS_MATCH | 3.0223 | -0.026 | 0.35 | [-0.045, -0.004] | 75/51 | 0.001 ms | no (3.553 alone) | null; grid artefact; +0.0003 on the built chain | -0.037 (0.47x), replicates |
| 2 | DIS+DISTPOT | 3.0299 | -0.018 | 0.17 | [-0.097, +0.073] | 62/64 | 0.12 ms | no (3.199 alone) | null; below its own noise by -0.040 (0.35x) | -0.028 (0.25x) |
| 3 | REJ[POOLGO]->DIS | 3.0325 | -0.016 | 0.37 | [-0.041, +0.010] | 25/22 | 0.000 ms | no (3.949 alone) | null; -0.017 on the built chain (0.34x) | -0.024 (0.51x) |
| 4 | REJ[LEG_compactness]->DIS | 3.0340 | -0.014 | 0.30 | [-0.043, +0.018] | 9/13 | 0.2 ms | no | null (moves 22 targets) | -0.016 (0.33x) |
| 5 | DIS+0.25*LEG_hbond_longrange | 3.0354 | -0.013 | 0.39 | [-0.040, +0.025] | 58/58 | 0.2 ms | no (3.470 alone) | null; -0.014 on the built chain | -0.027 (0.69x) |
| 6 | REJ[RG_LAW]->DIS | 3.0356 | -0.013 | 0.27 | [-0.045, +0.021] | 8/11 | 0.000 ms | no | null | -0.022 (0.42x) |
| 12 | DIS+0.5*ENV | 3.0438 | -0.005 | 0.05 | [-0.101, +0.101] | 62/61 | 0.15 ms | no | null | -0.008 (0.09x) |
| 23 | **DIS (production)** | **3.0483** | 0 | | | | 0.014 ms | yes | | 3.0580 |
| 63 | DIS_MEAN | 3.0784 | +0.030 | 0.51 | [-0.017, +0.074] | | 0.000 ms | re-reading | no | +0.019 |
| 96 | DISTPOT | 3.1990 | +0.151 | 0.63 | [+0.019, +0.313] | 52/74 | 0.12 ms | best new single, worse than DIS | see rank 2 | +0.145 (0.60x) |
| 150 | LEG | 3.7543 | +0.706 | 2.13 | [+0.493, +0.840] | | 3 ms | no (S25) | no | S25 |
| 165 | CONS | 3.6878 | +0.639 | 2.02 | [+0.428, +0.810] | | 1.6 ms | no | no: +0.286 with DIS | +0.276 with DIS |
| 170 | AMB | 3.8810 | +0.833 | 2.68 | [+0.687, +0.983] | | 8 ms | no (S25) | no | S25 |

Computational cost is the per-evaluation figure of Appendix T2 (every new channel is 0.001 to
0.15 ms per candidate on this box; CONS 1.6 ms; the S26 timing of Legacy 3 ms and AMBER 8 ms
stands); no configuration is cost-limited.

## 9. Decisively rejected, with the sprint and table that rejects each

- The consistency family as a selector or as a complement (CONS, DMAP_CONS, TORS_CONS, POOLGO,
  their composite): worse alone (+0.64 to +0.90 A, 2x MDE), worse with the distogram (+0.26 to
  +0.36, worse than their own noise at 1.7x MDE), worse at every m and in-band (T1, T3, T5, T6);
  the mechanism is section 6.
- Backbone physics alone or added (DSSPHB, ELEC, the physics composite): worse than random alone;
  DSSPHB with the distogram +0.20 A at 1.15x MDE, worse than noise (T3).
- Torsion and CA-geometry statistics (RAMA, CAGEO, LEG_torsion): worse alone (2.0 to 2.3x MDE)
  and worse than noise added (T3).
- Compactness targets (RG_LAW, RG_UNIV, LEG_compactness): worse alone and added; the reject form
  is null (T1).
- Excluded volume (EXVOL, LEG_steric): 99% of candidates tie at zero, so the channel is empty as
  a selector (the VQE returns the uniform state on it, T9) and null as a reject (T1).
- Every Legacy term alone: all worse than the distogram, ten of eleven worse than random (T2).
- Re-readings of the posterior (DIS_MEAN, CONTACT_LL): worse alone, worse or null added (S25 L12
  extended).
- Adaptive entropy-gated mixtures (T1) and the propensity gate (T7): null; routers stay closed.
- The 0.5-weight SS_MATCH bump: a grid artefact (T4: held-out -0.008 A, split-half 10%).

## 10. The strongest scientifically justified candidate, and what to follow up

There is no candidate for the next architecture: the shipped distogram Bayes-risk score stands
as the best Hamiltonian on this readout, alone and inside the CVaR-VQE, and the sprint's
mechanism (section 6) says why a better ranker cannot help an averaging readout. The follow-ups
worth an experiment, in order, each with the artefact that motivates it:

1. **DISTPOT as the second channel of a NEW readout, not this one.** It is the only new
   Hamiltonian with real, non-redundant information (partial rho +0.09, better than every physics
   energy by 0.5 A alone, below its own noise in combination) and the only one that does not
   collapse diversity. The experiment is not another weight: it is a readout that consumes
   ranking information without averaging (a cluster-then-average or a diversity-preserving
   selection, S21's diversity arms revisited with DISTPOT as the ranker), pre-registered against
   the set-mean law. Cost: one sprint lane; the channel is 0.12 ms per candidate.
2. **A FAIL18 detector with DISTPOT and ENV as inputs** (T-strata: every alternative helps by
   0.2 to 0.4 A on the 18 targets where the distogram is wrong). Seven routers have failed; the
   new inputs are the pool's own statistical-potential distribution (its spread, its agreement
   with the distogram). Prior: null; worth one nested-CV attempt because the prize is 0.03 A on
   the mean if the switch is 50% accurate.
3. **The set-mean law under a ranker with partial rho 0.25.** CONS is the best in-pool ranker on
   record after the distogram and the worst selector; the S12 medoid readout (-0.172 A) was the
   one operator that used it. A pre-registered medoid-plus-average hybrid (the medoid of the
   distogram's top-75 under CONS, averaged with its nearest 20) tests whether any of that ranking
   information can reach the endpoint. Prior: the S23 L5 outlier-removal cost says no.
4. **SS_MATCH2 on a longer instrument only.** Its effect, if real, is 0.01 to 0.03 A, four times
   below what 126 targets resolve; it needs about 1,000 targets, which this project cannot
   assemble (no fresh benchmark exists). Recorded so nobody re-runs it at n = 126.
5. **Nothing quantum.** No Hamiltonian changed the selector's trainability (T11) or its verdict
   (T9); the CVaR-VQE's answer is the classical ordering's answer for every energy tried.

## 11. What damaged the pre-registered expectations

- H1's prior said every single would be null-to-worse; it was worse, and mostly worse than
  random. The prior was too kind to the knowledge-based family: DISTPOT at 3.199 was the only
  one to land clearly better than random.
- H2's prior named DIS+CONS as the pair whose success would surprise least; it was the pair
  that failed hardest (+0.286, worse than noise at 1.7x MDE), and H9's prediction that it would
  win at small m was falsified at every m. The consensus intuition (S12) does not transfer from
  the medoid readout to the average, and the reason is now measured.
- The SS_MATCH bump at w = 0.5 looked like the sprint's one lead through waves 1 and 2 (the same
  weight chosen on all five folds, the effect present at every m and in-band); completing the
  weight grid dissolved it. The grid-oracle rule earned its keep once more.

## 12. What was not done

- No AMBER compute was spent: every AMBER number is the S25 cache (63,000 single points).
- The built-chain endpoint was run for 10 configurations, not 203; every configuration whose
  point-cloud effect was inside 0.4x MDE and positive was not projected (their point-cloud
  verdicts were already null-to-worse, and the projection has never reversed a sign on this
  record).
- The CVaR-VQE arm used the S25 suite's settings (9 qubits over the whole pool), not the
  production selector's (7 qubits over the top-128 with the VQE_LFO table): the production form
  is a special case of the same theorem (tail subset of prefix) and was verified bit-exactly in
  S25 and S26; re-running it on 80 energies would not change a verdict the gate has already
  fixed.
- Learned potentials (a trained pair network as a Hamiltonian) were excluded by the brief's
  "learned-free" spirit; S26's pairnet rung (L103) is the record for that direction (null).

---

## Appendix A. The tables (`s27/results/tables.md`, generated by `s27/make_tables.py`)
## T1. Every configuration, ranked by mean point-cloud RMSD of the top-75 average (n = 126; DIS = 3.0483, random-75 = 3.4209)

| rank | configuration | kind | mean RMSD (A) | vs DIS (A) | x MDE | fold-clustered 95% CI | W/L | vs random-75 (A) | verdict |
|---:|---|---|---:|---:|---:|---|---|---:|---|
| 1 | `DIS+0.5*SS_MATCH` | weighted | 3.0223 | -0.0261 | -0.35 | [-0.045, -0.004] | 75/51 | -0.3986 | NULL |
| 2 | `DIS+DISTPOT` | pair | 3.0299 | -0.0185 | -0.17 | [-0.097, +0.073] | 62/64 | -0.3910 | NULL |
| 3 | `REJ[POOLGO]->DIS` | staged | 3.0325 | -0.0159 | -0.37 | [-0.041, +0.010] | 25/22 | -0.3884 | NULL |
| 4 | `REJ[LEG_compactness]->DIS` | staged | 3.0340 | -0.0143 | -0.30 | [-0.043, +0.018] | 9/13 | -0.3869 | NULL |
| 5 | `DIS+0.25*LEG_hbond_longrange` | weighted | 3.0354 | -0.0130 | -0.39 | [-0.040, +0.025] | 58/58 | -0.3855 | NULL |
| 6 | `REJ[RG_LAW]->DIS` | staged | 3.0356 | -0.0127 | -0.27 | [-0.045, +0.021] | 8/11 | -0.3853 | NULL |
| 7 | `REJ[RG_UNIV]->DIS` | staged | 3.0366 | -0.0118 | -0.25 | [-0.043, +0.021] | 8/14 | -0.3843 | NULL |
| 8 | `REJ[DSSPHB]->DIS` | staged | 3.0404 | -0.0080 | -0.23 | [-0.034, +0.021] | 33/34 | -0.3805 | NULL |
| 9 | `REJ[DMAP_CONS]->DIS` | staged | 3.0404 | -0.0079 | -0.16 | [-0.035, +0.020] | 9/15 | -0.3804 | NULL |
| 10 | `DIS+0.25*SS_MATCH` | weighted | 3.0422 | -0.0061 | -0.15 | [-0.017, +0.004] | 76/49 | -0.3787 | NULL |
| 11 | `DIS+0.25*LEG_coop_sheet` | weighted | 3.0431 | -0.0053 | -0.21 | [-0.019, +0.006] | 40/38 | -0.3778 | NULL |
| 12 | `DIS+0.5*ENV` | weighted | 3.0438 | -0.0045 | -0.05 | [-0.101, +0.101] | 62/61 | -0.3771 | NULL |
| 13 | `REJ[CONTACT_LL]->DIS` | staged | 3.0444 | -0.0039 | -0.24 | [-0.019, +0.006] | 6/7 | -0.3764 | NULL |
| 14 | `REJ[LEG_solvation]->DIS` | staged | 3.0448 | -0.0036 | -0.16 | [-0.024, +0.019] | 29/32 | -0.3761 | NULL |
| 15 | `DIS+0.5*DISTPOT` | weighted | 3.0451 | -0.0033 | -0.05 | [-0.048, +0.044] | 62/64 | -0.3758 | NULL |
| 16 | `DIS+0.25*HP` | weighted | 3.0456 | -0.0028 | -0.09 | [-0.016, +0.011] | 57/60 | -0.3753 | NULL |
| 17 | `DIS+ENV` | pair | 3.0462 | -0.0021 | -0.02 | [-0.128, +0.113] | 64/62 | -0.3747 | NULL |
| 18 | `DIS+0.5*LEG_coop_sheet` | weighted | 3.0466 | -0.0018 | -0.03 | [-0.056, +0.058] | 60/56 | -0.3743 | NULL |
| 19 | `REJ[LEG]->DIS` | staged | 3.0470 | -0.0014 | -0.03 | [-0.023, +0.023] | 21/29 | -0.3739 | NULL |
| 20 | `REJ[TORS_CONS]->DIS` | staged | 3.0474 | -0.0009 | -0.03 | [-0.015, +0.011] | 33/31 | -0.3735 | NULL |
| 21 | `REJ[SS_MATCH]->DIS` | staged | 3.0477 | -0.0006 | -0.03 | [-0.006, +0.005] | 38/41 | -0.3732 | NULL |
| 22 | `REJ[DISTPOT]->DIS` | staged | 3.0480 | -0.0003 | -0.01 | [-0.020, +0.022] | 18/23 | -0.3729 | NULL |
| 23 | `DIS` | single | 3.0483 | +0.0000 | +nan | [+0.000, +0.000] | 0/0 | -0.3726 | reference |
| 24 | `REJ[DIS_MEAN]->DIS` | staged | 3.0483 | +0.0000 | +nan | [+0.000, +0.000] | 0/0 | -0.3726 | reference |
| 25 | `DIS+0.5*LEG_hbond_longrange` | weighted | 3.0487 | +0.0003 | +0.01 | [-0.048, +0.052] | 55/71 | -0.3722 | NULL |
| 26 | `REJ[LEG_hbond_local]->DIS` | staged | 3.0490 | +0.0006 | +0.02 | [-0.017, +0.021] | 39/40 | -0.3719 | NULL |
| 27 | `DIS+0.25*ENV` | weighted | 3.0522 | +0.0039 | +0.08 | [-0.034, +0.048] | 68/53 | -0.3687 | NULL |
| 28 | `DIS+LEG_aromatic~perm` | pair_perm | 3.0544 | +0.0061 | +0.07 | [-0.035, +0.064] | 32/30 | -0.3665 | NULL |
| 29 | `DIS+EXVOL~perm` | pair_perm | 3.0549 | +0.0066 | +0.39 | [-0.000, +0.013] | 54/49 | -0.3660 | NULL |
| 30 | `DIS+0.25*LEG_aromatic` | weighted | 3.0552 | +0.0069 | +0.22 | [-0.008, +0.020] | 31/31 | -0.3657 | NULL |
| 31 | `REJ[random]->DIS` | staged_ctrl | 3.0557 | +0.0073 | +0.30 | [-0.005, +0.022] | 56/70 | -0.3652 | NULL |
| 32 | `DIS+0.25*DISTPOT` | weighted | 3.0562 | +0.0079 | +0.21 | [-0.013, +0.032] | 57/69 | -0.3647 | NULL |
| 33 | `DIS+adapt(DISTPOT)` | adaptive | 3.0569 | +0.0086 | +0.22 | [-0.015, +0.040] | 57/69 | -0.3640 | NULL |
| 34 | `DIS+0.25*RG_UNIV` | weighted | 3.0570 | +0.0087 | +0.18 | [-0.016, +0.047] | 70/56 | -0.3639 | NULL |
| 35 | `DIS+0.25*DIS_MEAN` | weighted | 3.0571 | +0.0087 | +0.38 | [+0.001, +0.019] | 61/65 | -0.3638 | NULL |
| 36 | `REJ[ENV]->DIS` | staged | 3.0577 | +0.0093 | +0.35 | [-0.013, +0.038] | 20/28 | -0.3632 | NULL |
| 37 | `DIS+0.25*LEG_solvation` | weighted | 3.0577 | +0.0094 | +0.23 | [-0.028, +0.051] | 59/67 | -0.3632 | NULL |
| 38 | `DIS+0.25*RG_LAW` | weighted | 3.0581 | +0.0097 | +0.20 | [-0.031, +0.053] | 62/64 | -0.3628 | NULL |
| 39 | `DIS+0.5*HP` | weighted | 3.0581 | +0.0098 | +0.14 | [-0.033, +0.058] | 56/69 | -0.3628 | NULL |
| 40 | `DIS+LEG_coop_sheet` | pair | 3.0588 | +0.0105 | +0.13 | [-0.048, +0.076] | 59/57 | -0.3621 | NULL |
| 41 | `DIS+STAT_COMP` | pair | 3.0591 | +0.0108 | +0.12 | [-0.053, +0.091] | 59/67 | -0.3618 | NULL |
| 42 | `DIS+ROSETTA_LIKE` | pair | 3.0603 | +0.0120 | +0.14 | [-0.057, +0.103] | 66/60 | -0.3606 | NULL |
| 43 | `DIS+0.5*LEG_solvation` | weighted | 3.0606 | +0.0123 | +0.14 | [-0.074, +0.109] | 50/76 | -0.3603 | NULL |
| 44 | `REJ[LEG_aromatic]->DIS` | staged | 3.0608 | +0.0124 | +0.35 | [-0.018, +0.042] | 37/57 | -0.3601 | NULL |
| 45 | `DIS+0.25*LEG_compactness` | weighted | 3.0611 | +0.0128 | +0.30 | [-0.020, +0.047] | 57/66 | -0.3597 | NULL |
| 46 | `REJ[LEG_electrostatic]->DIS` | staged | 3.0612 | +0.0128 | +0.56 | [-0.001, +0.025] | 38/56 | -0.3597 | NULL |
| 47 | `DIS+EXVOL` | pair | 3.0612 | +0.0129 | +0.65 | [+0.005, +0.021] | 19/26 | -0.3597 | NULL |
| 48 | `DIS+0.25*EXVOL` | weighted | 3.0612 | +0.0129 | +0.65 | [+0.005, +0.021] | 19/26 | -0.3597 | NULL |
| 49 | `DIS+0.5*EXVOL` | weighted | 3.0612 | +0.0129 | +0.65 | [+0.005, +0.021] | 19/26 | -0.3597 | NULL |
| 50 | `DIS+0.25*ELEC` | weighted | 3.0620 | +0.0136 | +0.39 | [-0.010, +0.037] | 50/58 | -0.3589 | NULL |
| 51 | `REJ[CAGEO]->DIS` | staged | 3.0620 | +0.0137 | +0.40 | [+0.000, +0.028] | 33/38 | -0.3589 | NULL |
| 52 | `DIS+0.25*LEG_electrostatic` | weighted | 3.0622 | +0.0138 | +0.39 | [-0.011, +0.037] | 54/54 | -0.3587 | NULL |
| 53 | `DIS+0.5*DIS_MEAN` | weighted | 3.0623 | +0.0140 | +0.46 | [-0.001, +0.028] | 57/69 | -0.3586 | NULL |
| 54 | `DIS+0.5*LEG_aromatic` | weighted | 3.0624 | +0.0141 | +0.18 | [-0.017, +0.050] | 30/32 | -0.3585 | NULL |
| 55 | `DIS+ELEC~perm` | pair_perm | 3.0634 | +0.0150 | +0.16 | [-0.024, +0.044] | 59/49 | -0.3575 | NULL |
| 56 | `REJ[LEG_contact]->DIS` | staged | 3.0637 | +0.0153 | +0.47 | [-0.007, +0.038] | 32/50 | -0.3572 | NULL |
| 57 | `REJ[ELEC]->DIS` | staged | 3.0645 | +0.0161 | +0.66 | [+0.000, +0.027] | 38/57 | -0.3564 | NULL |
| 58 | `DIS+CAGEO~perm` | pair_perm | 3.0646 | +0.0163 | +0.17 | [-0.045, +0.097] | 68/58 | -0.3562 | NULL |
| 59 | `DIS+DIS_MEAN` | pair | 3.0673 | +0.0190 | +0.49 | [-0.012, +0.047] | 53/73 | -0.3536 | NULL |
| 60 | `DIS+CONS~perm` | pair_perm | 3.0678 | +0.0195 | +0.20 | [-0.024, +0.071] | 67/59 | -0.3531 | NULL |
| 61 | `REJ[LEG_hbond_longrange]->DIS` | staged | 3.0685 | +0.0201 | +0.72 | [+0.001, +0.040] | 55/71 | -0.3524 | TYPE-M ZONE |
| 62 | `DIS+LEG_steric~perm` | pair_perm | 3.0693 | +0.0210 | +0.62 | [+0.003, +0.036] | 64/61 | -0.3516 | NULL |
| 63 | `DIS+DISTPOT~perm` | pair_perm | 3.0695 | +0.0212 | +0.23 | [-0.032, +0.089] | 67/59 | -0.3514 | NULL |
| 64 | `DIS+0.25*DSSPHB` | weighted | 3.0697 | +0.0214 | +0.43 | [+0.001, +0.042] | 52/74 | -0.3512 | NULL |
| 65 | `DIS+SS_MATCH` | pair | 3.0706 | +0.0222 | +0.20 | [-0.033, +0.090] | 58/68 | -0.3503 | NULL |
| 66 | `DIS+POOLGO~perm` | pair_perm | 3.0706 | +0.0223 | +0.23 | [-0.025, +0.083] | 62/64 | -0.3503 | NULL |
| 67 | `REJ[LEG_coop_sheet]->DIS` | staged | 3.0706 | +0.0223 | +0.81 | [+0.007, +0.038] | 52/74 | -0.3503 | TYPE-M ZONE |
| 68 | `REJ[HP]->DIS` | staged | 3.0708 | +0.0224 | +1.03 | [+0.012, +0.034] | 29/36 | -0.3501 | WORSE |
| 69 | `DIS+0.5*RG_UNIV` | weighted | 3.0708 | +0.0225 | +0.20 | [-0.030, +0.079] | 65/61 | -0.3501 | NULL |
| 70 | `DIS+adapt(DSSPHB)` | adaptive | 3.0708 | +0.0225 | +0.44 | [-0.000, +0.042] | 54/72 | -0.3501 | NULL |
| 71 | `REJ[CONS]->DIS` | staged | 3.0709 | +0.0226 | +0.44 | [+0.000, +0.045] | 13/20 | -0.3500 | NULL |
| 72 | `DIS+0.5*LEG_compactness` | weighted | 3.0714 | +0.0231 | +0.27 | [-0.043, +0.102] | 58/65 | -0.3495 | NULL |
| 73 | `DIS+CONTACT~perm` | pair_perm | 3.0717 | +0.0233 | +0.24 | [-0.027, +0.079] | 65/61 | -0.3492 | NULL |
| 74 | `DIS+SS_MATCH~perm` | pair_perm | 3.0721 | +0.0237 | +0.24 | [-0.023, +0.078] | 66/60 | -0.3488 | NULL |
| 75 | `DIS+0.25*LEG_hbond_local` | weighted | 3.0727 | +0.0243 | +0.48 | [+0.001, +0.045] | 54/72 | -0.3482 | NULL |
| 76 | `DIS+0.25*LEG_contact` | weighted | 3.0737 | +0.0254 | +0.70 | [+0.001, +0.048] | 62/64 | -0.3472 | NULL |
| 77 | `REJ[LEG_coop_helix]->DIS` | staged | 3.0739 | +0.0256 | +0.69 | [+0.009, +0.044] | 36/53 | -0.3470 | NULL |
| 78 | `REJ[EXVOL]->DIS` | staged | 3.0743 | +0.0260 | +0.98 | [+0.019, +0.033] | 49/77 | -0.3466 | TYPE-M ZONE |
| 79 | `DIS+0.25*CONTACT_LL` | weighted | 3.0744 | +0.0261 | +0.59 | [+0.016, +0.037] | 52/64 | -0.3465 | NULL |
| 80 | `DIS+LEG_hbond_local~perm` | pair_perm | 3.0748 | +0.0265 | +0.27 | [-0.017, +0.071] | 65/61 | -0.3461 | NULL |
| 81 | `DIS+PHYSICS_COMP` | pair | 3.0756 | +0.0272 | +0.47 | [-0.011, +0.066] | 57/69 | -0.3453 | NULL |
| 82 | `DIS+LEG_contact~perm` | pair_perm | 3.0756 | +0.0272 | +0.28 | [-0.007, +0.065] | 61/65 | -0.3453 | NULL |
| 83 | `DIS+0.25*LEG_torsion` | weighted | 3.0757 | +0.0273 | +0.51 | [+0.022, +0.034] | 49/77 | -0.3452 | NULL |
| 84 | `DIS+RAMA~perm` | pair_perm | 3.0757 | +0.0274 | +0.28 | [-0.028, +0.078] | 60/66 | -0.3452 | NULL |
| 85 | `DIS+0.25*LEG` | weighted | 3.0757 | +0.0274 | +0.58 | [-0.005, +0.063] | 54/72 | -0.3452 | NULL |
| 86 | `DIS+LEG_electrostatic~perm` | pair_perm | 3.0760 | +0.0276 | +0.30 | [-0.020, +0.059] | 52/56 | -0.3449 | NULL |
| 87 | `DIS+LEG_solvation~perm` | pair_perm | 3.0760 | +0.0276 | +0.30 | [-0.026, +0.091] | 65/61 | -0.3449 | NULL |
| 88 | `REJ[RAMA]->DIS` | staged | 3.0761 | +0.0278 | +0.89 | [+0.018, +0.039] | 38/49 | -0.3448 | TYPE-M ZONE |
| 89 | `DIS+adapt(CONTACT)` | adaptive | 3.0771 | +0.0287 | +0.57 | [-0.005, +0.079] | 58/59 | -0.3438 | NULL |
| 90 | `DIS+LEG_compactness~perm` | pair_perm | 3.0779 | +0.0296 | +0.29 | [-0.012, +0.076] | 64/62 | -0.3430 | NULL |
| 91 | `DIS+0.25*CAGEO` | weighted | 3.0779 | +0.0296 | +0.54 | [+0.009, +0.050] | 60/66 | -0.3430 | NULL |
| 92 | `DIS+0.25*POOLGO` | weighted | 3.0782 | +0.0298 | +0.61 | [+0.008, +0.048] | 52/66 | -0.3427 | NULL |
| 93 | `DIS_MEAN` | single | 3.0784 | +0.0301 | +0.51 | [-0.017, +0.074] | 56/70 | -0.3425 | NULL |
| 94 | `REJ[AMB]->DIS` | staged | 3.0788 | +0.0304 | +0.87 | [+0.016, +0.042] | 41/63 | -0.3421 | TYPE-M ZONE |
| 95 | `REJ[LEG_steric]->DIS` | staged | 3.0788 | +0.0304 | +0.90 | [+0.016, +0.045] | 55/57 | -0.3421 | TYPE-M ZONE |
| 96 | `DIS+0.25*TORS_CONS` | weighted | 3.0788 | +0.0305 | +0.53 | [+0.012, +0.051] | 61/65 | -0.3420 | NULL |
| 97 | `DIS+RG_LAW~perm` | pair_perm | 3.0791 | +0.0308 | +0.32 | [-0.027, +0.084] | 65/61 | -0.3418 | NULL |
| 98 | `DIS+adapt(LEG)` | adaptive | 3.0793 | +0.0310 | +0.65 | [-0.003, +0.063] | 58/68 | -0.3415 | NULL |
| 99 | `DIS+0.25*CONTACT` | weighted | 3.0795 | +0.0312 | +0.66 | [-0.002, +0.078] | 59/59 | -0.3414 | NULL |
| 100 | `DIS+0.5*LEG_contact` | weighted | 3.0800 | +0.0317 | +0.53 | [-0.006, +0.065] | 56/70 | -0.3409 | NULL |
| 101 | `REJ[CONTACT]->DIS` | staged | 3.0803 | +0.0320 | +0.74 | [+0.015, +0.051] | 30/41 | -0.3406 | TYPE-M ZONE |
| 102 | `DIS+LEG_coop_sheet~perm` | pair_perm | 3.0803 | +0.0320 | +0.68 | [-0.002, +0.059] | 48/67 | -0.3406 | NULL |
| 103 | `DIS+0.25*LEG_steric` | weighted | 3.0806 | +0.0323 | +0.86 | [+0.018, +0.047] | 46/60 | -0.3402 | TYPE-M ZONE |
| 104 | `DIS+DIS_MEAN~perm` | pair_perm | 3.0807 | +0.0324 | +0.36 | [-0.011, +0.076] | 63/63 | -0.3402 | NULL |
| 105 | `DIS+ENV~perm` | pair_perm | 3.0807 | +0.0324 | +0.37 | [-0.005, +0.075] | 67/59 | -0.3402 | NULL |
| 106 | `REJ[LEG_torsion]->DIS` | staged | 3.0809 | +0.0325 | +1.01 | [+0.021, +0.046] | 34/54 | -0.3400 | WORSE |
| 107 | `DIS+0.5*LEG_steric` | weighted | 3.0820 | +0.0337 | +0.78 | [+0.017, +0.048] | 47/59 | -0.3389 | TYPE-M ZONE |
| 108 | `DIS+LEG_steric` | pair | 3.0821 | +0.0338 | +0.79 | [+0.017, +0.048] | 47/59 | -0.3388 | TYPE-M ZONE |
| 109 | `DIS+RG_UNIV~perm` | pair_perm | 3.0824 | +0.0340 | +0.37 | [-0.009, +0.087] | 61/65 | -0.3385 | NULL |
| 110 | `DIS+0.5*RG_LAW` | weighted | 3.0826 | +0.0343 | +0.35 | [-0.022, +0.121] | 59/67 | -0.3383 | NULL |
| 111 | `DIS+LEG~perm` | pair_perm | 3.0828 | +0.0344 | +0.35 | [-0.014, +0.085] | 64/62 | -0.3381 | NULL |
| 112 | `DIS+DSSPHB~perm` | pair_perm | 3.0830 | +0.0347 | +0.35 | [-0.019, +0.095] | 65/61 | -0.3379 | NULL |
| 113 | `DIS+adapt(AMB)` | adaptive | 3.0831 | +0.0348 | +0.83 | [+0.015, +0.052] | 56/70 | -0.3378 | TYPE-M ZONE |
| 114 | `DIS+0.25*AMB` | weighted | 3.0850 | +0.0366 | +0.90 | [+0.020, +0.051] | 58/68 | -0.3359 | TYPE-M ZONE |
| 115 | `DIS+LEG_torsion~perm` | pair_perm | 3.0857 | +0.0374 | +0.39 | [-0.018, +0.097] | 65/61 | -0.3352 | NULL |
| 116 | `DIS+0.5*CONTACT` | weighted | 3.0863 | +0.0379 | +0.46 | [-0.009, +0.122] | 59/64 | -0.3346 | NULL |
| 117 | `DIS+LEG_aromatic` | pair | 3.0875 | +0.0392 | +0.31 | [-0.005, +0.091] | 28/34 | -0.3334 | NULL |
| 118 | `DIS+DMAP_CONS~perm` | pair_perm | 3.0878 | +0.0394 | +0.43 | [+0.011, +0.072] | 64/62 | -0.3331 | NULL |
| 119 | `DIS+0.25*RAMA` | weighted | 3.0881 | +0.0397 | +0.80 | [+0.029, +0.056] | 50/76 | -0.3328 | TYPE-M ZONE |
| 120 | `DIS+0.5*LEG_electrostatic` | weighted | 3.0887 | +0.0404 | +0.66 | [-0.011, +0.087] | 46/62 | -0.3322 | NULL |
| 121 | `DIS+TORS_CONS~perm` | pair_perm | 3.0890 | +0.0407 | +0.43 | [-0.004, +0.081] | 58/68 | -0.3319 | NULL |
| 122 | `DIS+HP~perm` | pair_perm | 3.0896 | +0.0413 | +0.44 | [-0.002, +0.087] | 64/62 | -0.3313 | NULL |
| 123 | `DIS+0.5*ELEC` | weighted | 3.0898 | +0.0415 | +0.66 | [-0.011, +0.090] | 46/62 | -0.3311 | NULL |
| 124 | `DIS+adapt(RAMA)` | adaptive | 3.0910 | +0.0427 | +0.75 | [+0.029, +0.059] | 51/75 | -0.3299 | TYPE-M ZONE |
| 125 | `DIS+0.25*LEG_coop_helix` | weighted | 3.0914 | +0.0431 | +0.77 | [+0.015, +0.076] | 38/61 | -0.3295 | TYPE-M ZONE |
| 126 | `DIS+CONTACT` | pair | 3.0933 | +0.0450 | +0.33 | [-0.054, +0.199] | 60/65 | -0.3276 | NULL |
| 127 | `DIS+0.25*DMAP_CONS` | weighted | 3.0936 | +0.0452 | +0.73 | [+0.013, +0.069] | 54/72 | -0.3273 | TYPE-M ZONE |
| 128 | `DIS+adapt(CONS)` | adaptive | 3.0938 | +0.0454 | +0.66 | [+0.012, +0.070] | 59/67 | -0.3271 | NULL |
| 129 | `DIS+0.5*AMB` | weighted | 3.0940 | +0.0456 | +0.61 | [+0.010, +0.077] | 58/68 | -0.3269 | NULL |
| 130 | `DIS+0.25*CONS` | weighted | 3.0944 | +0.0461 | +0.70 | [+0.017, +0.068] | 58/68 | -0.3265 | TYPE-M ZONE |
| 131 | `DIS+HP` | pair | 3.0948 | +0.0465 | +0.35 | [-0.044, +0.164] | 54/72 | -0.3261 | NULL |
| 132 | `DIS+CONTACT_LL~perm` | pair_perm | 3.0948 | +0.0465 | +0.46 | [-0.007, +0.103] | 63/63 | -0.3260 | NULL |
| 133 | `DIS+AMB~perm` | pair_perm | 3.0996 | +0.0513 | +0.52 | [-0.001, +0.099] | 59/67 | -0.3213 | NULL |
| 134 | `DIS+0.5*CONTACT_LL` | weighted | 3.1000 | +0.0517 | +0.83 | [+0.026, +0.074] | 51/70 | -0.3209 | TYPE-M ZONE |
| 135 | `DIS+0.5*DSSPHB` | weighted | 3.1035 | +0.0551 | +0.55 | [+0.005, +0.103] | 47/79 | -0.3174 | NULL |
| 136 | `DIS+0.5*LEG` | weighted | 3.1051 | +0.0568 | +0.58 | [+0.005, +0.112] | 51/75 | -0.3158 | NULL |
| 137 | `DIS+LEG_electrostatic` | pair | 3.1066 | +0.0583 | +0.52 | [-0.007, +0.119] | 48/60 | -0.3143 | NULL |
| 138 | `DIS+ELEC` | pair | 3.1076 | +0.0593 | +0.53 | [-0.007, +0.122] | 47/61 | -0.3133 | NULL |
| 139 | `DIS+LEG_coop_helix~perm` | pair_perm | 3.1100 | +0.0617 | +0.65 | [+0.009, +0.103] | 56/70 | -0.3109 | NULL |
| 140 | `DIS+0.5*LEG_torsion` | weighted | 3.1113 | +0.0630 | +0.68 | [+0.038, +0.083] | 51/75 | -0.3096 | NULL |
| 141 | `DIS+0.5*POOLGO` | weighted | 3.1174 | +0.0691 | +0.59 | [+0.002, +0.136] | 48/76 | -0.3035 | NULL |
| 142 | `DIS+LEG_solvation` | pair | 3.1213 | +0.0729 | +0.51 | [-0.046, +0.206] | 52/74 | -0.2996 | NULL |
| 143 | `DIS+0.5*LEG_hbond_local` | weighted | 3.1260 | +0.0776 | +0.75 | [+0.024, +0.125] | 51/75 | -0.2949 | TYPE-M ZONE |
| 144 | `DIS+AMB` | pair | 3.1296 | +0.0813 | +0.68 | [+0.038, +0.111] | 57/69 | -0.2913 | NULL |
| 145 | `DIS+0.5*CAGEO` | weighted | 3.1316 | +0.0833 | +0.77 | [+0.018, +0.149] | 55/71 | -0.2893 | TYPE-M ZONE |
| 146 | `DIS+0.5*RAMA` | weighted | 3.1318 | +0.0834 | +0.87 | [+0.041, +0.131] | 40/86 | -0.2891 | TYPE-M ZONE |
| 147 | `DIS+LEG_contact` | pair | 3.1343 | +0.0859 | +0.73 | [-0.009, +0.163] | 54/72 | -0.2866 | TYPE-M ZONE |
| 148 | `DIS+CONTACT_LL` | pair | 3.1367 | +0.0884 | +1.09 | [+0.047, +0.125] | 48/73 | -0.2842 | WORSE |
| 149 | `DIS+0.5*LEG_coop_helix` | weighted | 3.1441 | +0.0958 | +0.83 | [+0.030, +0.167] | 38/67 | -0.2768 | TYPE-M ZONE |
| 150 | `DIS+LEG_compactness` | pair | 3.1443 | +0.0960 | +0.53 | [-0.033, +0.246] | 53/70 | -0.2766 | NULL |
| 151 | `DIS+0.5*TORS_CONS` | weighted | 3.1445 | +0.0961 | +0.81 | [+0.029, +0.160] | 49/77 | -0.2764 | TYPE-M ZONE |
| 152 | `DIS+RG_LAW` | pair | 3.1590 | +0.1106 | +0.56 | [+0.009, +0.248] | 56/70 | -0.2619 | NULL |
| 153 | `DIS+0.5*DMAP_CONS` | weighted | 3.1648 | +0.1165 | +0.89 | [+0.023, +0.175] | 48/78 | -0.2561 | TYPE-M ZONE |
| 154 | `DIS+RG_UNIV` | pair | 3.1818 | +0.1335 | +0.74 | [+0.041, +0.222] | 51/75 | -0.2391 | TYPE-M ZONE |
| 155 | `DIS+LEG_torsion` | pair | 3.1841 | +0.1358 | +0.93 | [+0.021, +0.225] | 50/76 | -0.2368 | TYPE-M ZONE |
| 156 | `DIS+LEG_hbond_longrange~perm` | pair_perm | 3.1903 | +0.1419 | +0.85 | [+0.041, +0.223] | 57/69 | -0.2306 | TYPE-M ZONE |
| 157 | `DIS+0.5*CONS` | weighted | 3.1941 | +0.1458 | +1.06 | [+0.050, +0.219] | 49/77 | -0.2268 | WORSE |
| 158 | `DISTPOT` | single | 3.1990 | +0.1507 | +0.63 | [+0.019, +0.313] | 52/74 | -0.2219 | NULL |
| 159 | `STAT_COMP` | composite | 3.2109 | +0.1625 | +0.66 | [-0.010, +0.385] | 54/72 | -0.2100 | NULL |
| 160 | `DIS+LEG` | pair | 3.2131 | +0.1647 | +0.91 | [+0.059, +0.244] | 47/79 | -0.2078 | TYPE-M ZONE |
| 161 | `DIS+LEG_hbond_longrange` | pair | 3.2255 | +0.1772 | +0.83 | [+0.020, +0.365] | 53/73 | -0.1953 | TYPE-M ZONE |
| 162 | `DIS+DSSPHB` | pair | 3.2483 | +0.2000 | +1.15 | [+0.075, +0.296] | 40/86 | -0.1726 | WORSE |
| 163 | `DIS+RAMA` | pair | 3.2690 | +0.2207 | +1.24 | [+0.080, +0.361] | 37/89 | -0.1519 | WORSE |
| 164 | `DIS+TORS_CONS` | pair | 3.3109 | +0.2626 | +1.22 | [+0.099, +0.409] | 40/86 | -0.1100 | WORSE |
| 165 | `DIS+LEG_hbond_local` | pair | 3.3131 | +0.2648 | +1.31 | [+0.089, +0.398] | 44/82 | -0.1078 | WORSE |
| 166 | `DIS+CONS` | pair | 3.3341 | +0.2858 | +1.29 | [+0.082, +0.476] | 43/83 | -0.0868 | WORSE |
| 167 | `DIS+CAGEO` | pair | 3.3386 | +0.2902 | +1.33 | [+0.118, +0.470] | 42/84 | -0.0823 | WORSE |
| 168 | `DIS+CONSIST_COMP` | pair | 3.3427 | +0.2944 | +1.32 | [+0.100, +0.463] | 44/82 | -0.0782 | WORSE |
| 169 | `DIS+LEG_coop_helix` | pair | 3.3458 | +0.2975 | +1.54 | [+0.198, +0.403] | 28/77 | -0.0751 | WORSE |
| 170 | `DIS+POOLGO` | pair | 3.3477 | +0.2993 | +1.38 | [+0.123, +0.435] | 39/87 | -0.0732 | WORSE |
| 171 | `ROSETTA_LIKE` | composite | 3.3499 | +0.3016 | +0.99 | [+0.023, +0.667] | 50/76 | -0.0710 | TYPE-M ZONE |
| 172 | `LEG_coop_sheet` | single | 3.3838 | +0.3355 | +1.45 | [+0.193, +0.493] | 45/81 | -0.0371 | WORSE |
| 173 | `DIS+DMAP_CONS` | pair | 3.4042 | +0.3559 | +1.45 | [+0.146, +0.520] | 41/85 | -0.0167 | WORSE |
| 174 | `EXVOL` | single | 3.4503 | +0.4019 | +1.61 | [+0.252, +0.535] | 45/81 | +0.0294 | WORSE |
| 175 | `LEG_steric` | single | 3.4601 | +0.4117 | +1.64 | [+0.241, +0.545] | 45/81 | +0.0392 | WORSE |
| 176 | `ENV` | single | 3.4669 | +0.4185 | +1.28 | [+0.270, +0.602] | 46/80 | +0.0460 | WORSE |
| 177 | `LEG_hbond_longrange` | single | 3.4695 | +0.4211 | +1.36 | [+0.183, +0.706] | 43/83 | +0.0486 | WORSE |
| 178 | `RG_UNIV` | single | 3.5458 | +0.4974 | +1.62 | [+0.304, +0.638] | 43/83 | +0.1249 | WORSE |
| 179 | `PHYSICS_COMP` | composite | 3.5523 | +0.5040 | +1.81 | [+0.321, +0.686] | 40/86 | +0.1314 | WORSE |
| 180 | `SS_MATCH` | single | 3.5532 | +0.5049 | +1.88 | [+0.381, +0.624] | 39/87 | +0.1323 | WORSE |
| 181 | `CONTACT` | single | 3.5637 | +0.5153 | +1.77 | [+0.287, +0.776] | 42/84 | +0.1428 | WORSE |
| 182 | `CONTACT_LL` | single | 3.5698 | +0.5214 | +1.93 | [+0.366, +0.627] | 48/78 | +0.1489 | WORSE |
| 183 | `LEG_aromatic` | single | 3.5701 | +0.5217 | +1.87 | [+0.397, +0.689] | 39/87 | +0.1492 | WORSE |
| 184 | `RG_LAW` | single | 3.5863 | +0.5380 | +1.69 | [+0.324, +0.777] | 42/84 | +0.1654 | WORSE |
| 185 | `LEG_solvation` | single | 3.6339 | +0.5856 | +1.89 | [+0.324, +0.896] | 40/86 | +0.2130 | WORSE |
| 186 | `CONS` | single | 3.6878 | +0.6394 | +2.02 | [+0.428, +0.810] | 37/89 | +0.2669 | WORSE |
| 187 | `LEG_torsion` | single | 3.6904 | +0.6421 | +2.02 | [+0.464, +0.822] | 41/85 | +0.2695 | WORSE |
| 188 | `LEG_compactness` | single | 3.7129 | +0.6645 | +1.85 | [+0.372, +0.982] | 38/88 | +0.2920 | WORSE |
| 189 | `LEG` | single | 3.7543 | +0.7059 | +2.13 | [+0.493, +0.840] | 38/88 | +0.3334 | WORSE |
| 190 | `DSSPHB` | single | 3.7839 | +0.7356 | +2.14 | [+0.531, +0.883] | 42/84 | +0.3630 | WORSE |
| 191 | `LEG_coop_helix` | single | 3.8120 | +0.7637 | +2.23 | [+0.569, +0.930] | 39/87 | +0.3911 | WORSE |
| 192 | `TORS_CONS` | single | 3.8269 | +0.7785 | +2.09 | [+0.597, +0.948] | 37/89 | +0.4060 | WORSE |
| 193 | `HP` | single | 3.8337 | +0.7853 | +2.15 | [+0.565, +1.032] | 39/87 | +0.4128 | WORSE |
| 194 | `DMAP_CONS` | single | 3.8560 | +0.8077 | +2.12 | [+0.574, +0.999] | 36/90 | +0.4351 | WORSE |
| 195 | `LEG_electrostatic` | single | 3.8577 | +0.8093 | +2.53 | [+0.657, +0.980] | 27/99 | +0.4368 | WORSE |
| 196 | `ELEC` | single | 3.8693 | +0.8210 | +2.57 | [+0.678, +0.995] | 27/99 | +0.4484 | WORSE |
| 197 | `AMB` | single | 3.8810 | +0.8327 | +2.68 | [+0.687, +0.983] | 32/94 | +0.4601 | WORSE |
| 198 | `CONSIST_COMP` | composite | 3.9043 | +0.8560 | +2.24 | [+0.648, +1.034] | 39/87 | +0.4834 | WORSE |
| 199 | `LEG_hbond_local` | single | 3.9262 | +0.8779 | +2.27 | [+0.663, +1.061] | 39/87 | +0.5053 | WORSE |
| 200 | `RAMA` | single | 3.9290 | +0.8807 | +2.31 | [+0.672, +1.054] | 34/92 | +0.5081 | WORSE |
| 201 | `POOLGO` | single | 3.9493 | +0.9010 | +2.28 | [+0.691, +1.074] | 35/91 | +0.5284 | WORSE |
| 202 | `CAGEO` | single | 3.9667 | +0.9184 | +2.33 | [+0.704, +1.104] | 29/97 | +0.5459 | WORSE |
| 203 | `LEG_contact` | single | 4.4028 | +1.3545 | +3.93 | [+1.258, +1.487] | 23/103 | +0.9819 | WORSE |

## T2. The singles and composites with their ORACLE diagnostics (rho = Spearman with the candidate RMSD, whole pool and in-band < 3 A; rho_DIS = rank correlation with the shipped score; overlap = share of DIS's top-75 in the channel's top-75; cost per evaluation)

| channel | family | mean RMSD | vs DIS | x MDE | vs random | x MDE | rho pool | rho in-band | rho_DIS | overlap | tie frac | ms/eval |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `DIS` | reference | 3.0483 | +0.0000 | +nan | -0.3726 | -1.53 | +0.568 | +0.221 | +1.000 | 1.00 | 0.08 | nan |
| `DIS_MEAN` | distogram-rereading | 3.0784 | +0.0301 | +0.51 | -0.3425 | -1.52 | +0.559 | +0.176 | +0.958 | 0.81 | 0.07 | 0.000 |
| `DISTPOT` | statistical | 3.1990 | +0.1507 | +0.63 | -0.2219 | -0.90 | +0.339 | +0.159 | +0.426 | 0.31 | 0.10 | 0.117 |
| `STAT_COMP` | composite | 3.2109 | +0.1625 | +0.66 | -0.2100 | -0.97 | +0.340 | +0.161 | +nan | 0.32 | 0.08 | nan |
| `ROSETTA_LIKE` | composite | 3.3499 | +0.3016 | +0.99 | -0.0710 | -0.27 | +0.301 | +0.042 | +nan | 0.20 | 0.07 | nan |
| `LEG_coop_sheet` | legacy term | 3.3838 | +0.3355 | +1.45 | -0.0371 | -0.30 | -0.061 | +0.037 | -0.076 | 0.15 | 0.96 | 0.204 |
| `EXVOL` | geometric | 3.4503 | +0.4019 | +1.61 | +0.0294 | +0.54 | +0.018 | +0.017 | +0.056 | 0.16 | 0.99 | 0.000 |
| `LEG_steric` | legacy term | 3.4601 | +0.4117 | +1.64 | +0.0392 | +0.70 | +0.020 | +0.064 | +0.051 | 0.16 | 0.99 | 0.204 |
| `ENV` | statistical | 3.4669 | +0.4185 | +1.28 | +0.0460 | +0.16 | +0.225 | +0.049 | +0.268 | 0.19 | 0.46 | 0.151 |
| `LEG_hbond_longrange` | legacy term | 3.4695 | +0.4211 | +1.36 | +0.0486 | +0.18 | -0.077 | -0.043 | -0.107 | 0.12 | 0.36 | 0.204 |
| `RG_UNIV` | geometric | 3.5458 | +0.4974 | +1.62 | +0.1249 | +0.84 | +0.274 | +0.032 | +0.366 | 0.17 | 0.06 | 0.018 |
| `PHYSICS_COMP` | composite | 3.5523 | +0.5040 | +1.81 | +0.1314 | +0.71 | +0.128 | +0.024 | +nan | 0.16 | 0.07 | nan |
| `SS_MATCH` | compatibility | 3.5532 | +0.5049 | +1.88 | +0.1323 | +0.82 | -0.009 | -0.073 | +0.012 | 0.14 | 0.48 | 0.001 |
| `CONTACT` | statistical | 3.5637 | +0.5153 | +1.77 | +0.1428 | +0.59 | +0.099 | +0.044 | +0.121 | 0.19 | 0.38 | 0.125 |
| `CONTACT_LL` | distogram-rereading | 3.5698 | +0.5214 | +1.93 | +0.1489 | +0.57 | +0.314 | +0.202 | +0.589 | 0.44 | 0.60 | 0.000 |
| `LEG_aromatic` | legacy term | 3.5701 | +0.5217 | +1.87 | +0.1492 | +0.75 | +0.100 | +0.024 | +0.098 | 0.15 | 0.53 | 0.204 |
| `RG_LAW` | geometric | 3.5863 | +0.5380 | +1.69 | +0.1654 | +0.70 | +0.305 | +0.048 | +0.329 | 0.12 | 0.07 | 0.000 |
| `LEG_solvation` | legacy term | 3.6339 | +0.5856 | +1.89 | +0.2130 | +0.78 | +0.229 | +0.024 | +0.248 | 0.14 | 0.07 | 0.204 |
| `CONS` | consistency | 3.6878 | +0.6394 | +2.02 | +0.2669 | +1.45 | +0.425 | +0.103 | +0.494 | 0.29 | 0.08 | 1.625 |
| `LEG_torsion` | legacy term | 3.6904 | +0.6421 | +2.02 | +0.2695 | +1.50 | +0.168 | +0.157 | +0.214 | 0.28 | 0.09 | 0.204 |
| `LEG_compactness` | legacy term | 3.7129 | +0.6645 | +1.85 | +0.2920 | +1.03 | +0.277 | +0.009 | +0.274 | 0.10 | 0.84 | 0.204 |
| `LEG` | reference | 3.7543 | +0.7059 | +2.13 | +0.3334 | +2.01 | +0.307 | +0.096 | +0.367 | 0.22 | 0.08 | nan |
| `DSSPHB` | physics | 3.7839 | +0.7356 | +2.14 | +0.3630 | +2.02 | +0.271 | +0.120 | +0.309 | 0.20 | 0.07 | 0.026 |
| `LEG_coop_helix` | legacy term | 3.8120 | +0.7637 | +2.23 | +0.3911 | +2.30 | +0.260 | +0.059 | +0.296 | 0.15 | 0.99 | 0.204 |
| `TORS_CONS` | consistency | 3.8269 | +0.7785 | +2.09 | +0.4060 | +1.75 | +0.364 | +0.129 | +0.419 | 0.24 | 0.07 | 0.001 |
| `HP` | statistical | 3.8337 | +0.7853 | +2.15 | +0.4128 | +1.43 | +0.024 | -0.035 | +0.023 | 0.13 | 0.38 | 0.003 |
| `DMAP_CONS` | consistency | 3.8560 | +0.8077 | +2.12 | +0.4351 | +1.92 | +0.398 | +0.140 | +0.518 | 0.31 | 0.08 | 0.001 |
| `LEG_electrostatic` | legacy term | 3.8577 | +0.8093 | +2.53 | +0.4368 | +1.61 | -0.065 | -0.038 | -0.043 | 0.12 | 0.20 | 0.204 |
| `ELEC` | physics | 3.8693 | +0.8210 | +2.57 | +0.4484 | +1.65 | -0.066 | -0.038 | -0.042 | 0.12 | 0.20 | 0.000 |
| `AMB` | reference | 3.8810 | +0.8327 | +2.68 | +0.4601 | +2.48 | -0.027 | +0.098 | -0.019 | 0.15 | 0.08 | nan |
| `CONSIST_COMP` | composite | 3.9043 | +0.8560 | +2.24 | +0.4834 | +2.12 | +0.398 | +0.118 | +nan | 0.27 | 0.13 | nan |
| `LEG_hbond_local` | legacy term | 3.9262 | +0.8779 | +2.27 | +0.5053 | +2.23 | +0.291 | +0.109 | +0.340 | 0.22 | 0.08 | 0.204 |
| `RAMA` | torsion | 3.9290 | +0.8807 | +2.31 | +0.5081 | +2.28 | +0.244 | +0.135 | +0.299 | 0.26 | 0.10 | 0.001 |
| `POOLGO` | consistency | 3.9493 | +0.9010 | +2.28 | +0.5284 | +2.24 | +0.359 | +0.091 | +0.403 | 0.17 | 0.77 | 0.000 |
| `CAGEO` | geometric | 3.9667 | +0.9184 | +2.33 | +0.5459 | +2.29 | +0.298 | +0.161 | +0.350 | 0.27 | 0.24 | 0.132 |
| `LEG_contact` | legacy term | 4.4028 | +1.3545 | +3.93 | +0.9819 | +3.71 | -0.176 | +0.019 | -0.153 | 0.11 | 0.18 | 0.204 |

## T3. Equal-weight pairs DIS + X: against DIS and against the rank-permuted control (X's marginal kept, its correspondence destroyed)

| pair | mean | vs DIS | x MDE | fold CI | vs permuted | x MDE | fold CI (perm) | reading |
|---|---:|---:|---:|---|---:|---:|---|---|
| `DIS+DISTPOT` | 3.0299 | -0.0185 | -0.17 | [-0.097, +0.073] | -0.0396 | -0.35 | [-0.128, +0.057] | indistinguishable from noise |
| `DIS+ENV` | 3.0462 | -0.0021 | -0.02 | [-0.128, +0.113] | -0.0345 | -0.28 | [-0.157, +0.072] | indistinguishable from noise |
| `DIS+LEG_coop_sheet` | 3.0588 | +0.0105 | +0.13 | [-0.048, +0.076] | -0.0215 | -0.25 | [-0.080, +0.038] | indistinguishable from noise |
| `DIS+DIS_MEAN` | 3.0673 | +0.0190 | +0.49 | [-0.012, +0.047] | -0.0134 | -0.16 | [-0.044, +0.015] | indistinguishable from noise |
| `DIS+SS_MATCH` | 3.0706 | +0.0222 | +0.20 | [-0.033, +0.090] | -0.0015 | -0.01 | [-0.049, +0.054] | indistinguishable from noise |
| `DIS+HP` | 3.0948 | +0.0465 | +0.35 | [-0.044, +0.164] | +0.0052 | +0.04 | [-0.074, +0.101] | indistinguishable from noise |
| `DIS+EXVOL` | 3.0612 | +0.0129 | +0.65 | [+0.005, +0.021] | +0.0063 | +0.40 | [-0.007, +0.016] | indistinguishable from noise |
| `DIS+LEG_steric` | 3.0821 | +0.0338 | +0.79 | [+0.017, +0.048] | +0.0128 | +0.46 | [-0.004, +0.030] | indistinguishable from noise |
| `DIS+CONTACT` | 3.0933 | +0.0450 | +0.33 | [-0.054, +0.199] | +0.0216 | +0.16 | [-0.078, +0.161] | indistinguishable from noise |
| `DIS+AMB` | 3.1296 | +0.0813 | +0.68 | [+0.038, +0.111] | +0.0300 | +0.41 | [-0.021, +0.082] | indistinguishable from noise |
| `DIS+LEG_electrostatic` | 3.1066 | +0.0583 | +0.52 | [-0.007, +0.119] | +0.0307 | +0.24 | [-0.040, +0.101] | indistinguishable from noise |
| `DIS+LEG_aromatic` | 3.0875 | +0.0392 | +0.31 | [-0.005, +0.091] | +0.0331 | +0.37 | [-0.019, +0.086] | indistinguishable from noise |
| `DIS+LEG_hbond_longrange` | 3.2255 | +0.1772 | +0.83 | [+0.020, +0.365] | +0.0353 | +0.18 | [-0.127, +0.211] | indistinguishable from noise |
| `DIS+CONTACT_LL` | 3.1367 | +0.0884 | +1.09 | [+0.047, +0.125] | +0.0418 | +0.41 | [-0.039, +0.101] | indistinguishable from noise |
| `DIS+ELEC` | 3.1076 | +0.0593 | +0.53 | [-0.007, +0.122] | +0.0442 | +0.37 | [-0.024, +0.102] | indistinguishable from noise |
| `DIS+LEG_solvation` | 3.1213 | +0.0729 | +0.51 | [-0.046, +0.206] | +0.0453 | +0.37 | [-0.024, +0.143] | indistinguishable from noise |
| `DIS+LEG_contact` | 3.1343 | +0.0859 | +0.73 | [-0.009, +0.163] | +0.0587 | +0.43 | [-0.018, +0.129] | indistinguishable from noise |
| `DIS+LEG_compactness` | 3.1443 | +0.0960 | +0.53 | [-0.033, +0.246] | +0.0664 | +0.43 | [-0.031, +0.189] | indistinguishable from noise |
| `DIS+RG_LAW` | 3.1590 | +0.1106 | +0.56 | [+0.009, +0.248] | +0.0799 | +0.48 | [-0.002, +0.198] | indistinguishable from noise |
| `DIS+LEG_torsion` | 3.1841 | +0.1358 | +0.93 | [+0.021, +0.225] | +0.0984 | +1.06 | [-0.026, +0.189] | indistinguishable from noise |
| `DIS+RG_UNIV` | 3.1818 | +0.1335 | +0.74 | [+0.041, +0.222] | +0.0995 | +0.62 | [+0.007, +0.186] | worse than its own noise |
| `DIS+LEG` | 3.2131 | +0.1647 | +0.91 | [+0.059, +0.244] | +0.1303 | +1.07 | [+0.060, +0.179] | worse than its own noise |
| `DIS+DSSPHB` | 3.2483 | +0.2000 | +1.15 | [+0.075, +0.296] | +0.1653 | +1.50 | [+0.095, +0.230] | worse than its own noise |
| `DIS+RAMA` | 3.2690 | +0.2207 | +1.24 | [+0.080, +0.361] | +0.1933 | +1.71 | [+0.059, +0.308] | worse than its own noise |
| `DIS+TORS_CONS` | 3.3109 | +0.2626 | +1.22 | [+0.099, +0.409] | +0.2219 | +1.43 | [+0.089, +0.352] | worse than its own noise |
| `DIS+LEG_coop_helix` | 3.3458 | +0.2975 | +1.54 | [+0.198, +0.403] | +0.2358 | +1.69 | [+0.161, +0.316] | worse than its own noise |
| `DIS+LEG_hbond_local` | 3.3131 | +0.2648 | +1.31 | [+0.089, +0.398] | +0.2383 | +1.72 | [+0.106, +0.353] | worse than its own noise |
| `DIS+CONS` | 3.3341 | +0.2858 | +1.29 | [+0.082, +0.476] | +0.2663 | +1.74 | [+0.112, +0.417] | worse than its own noise |
| `DIS+CAGEO` | 3.3386 | +0.2902 | +1.33 | [+0.118, +0.470] | +0.2739 | +1.68 | [+0.114, +0.446] | worse than its own noise |
| `DIS+POOLGO` | 3.3477 | +0.2993 | +1.38 | [+0.123, +0.435] | +0.2770 | +1.74 | [+0.133, +0.407] | worse than its own noise |
| `DIS+DMAP_CONS` | 3.4042 | +0.3559 | +1.45 | [+0.146, +0.520] | +0.3165 | +1.71 | [+0.102, +0.484] | worse than its own noise |

## T4. Regulariser weights (H3): the weight chosen leave-fold-out, the held-out effect, and the grid priced as an order statistic

| complement | held-out effect | x MDE | fold CI | weights chosen per fold | full-leakage best | k_eff |
|---|---:|---:|---|---|---:|---:|
| `SS_MATCH` | -0.0261 | -0.35 | [-0.045, -0.005] | [0.5, 0.5, 0.5, 0.5, 0.5] | -0.0261 | 2.93 |
| `LEG_hbond_longrange` | -0.0130 | -0.39 | [-0.040, +0.025] | [0.25, 0.25, 0.25, 0.25, 0.25] | -0.0130 | 2.69 |
| `RG_UNIV` | +0.0087 | +0.18 | [-0.016, +0.048] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0087 | 2.79 |
| `DISTPOT` | +0.0096 | +0.12 | [-0.035, +0.075] | [1.0, 1.0, 1.0, 1.0, 0.25] | -0.0185 | 2.65 |
| `RG_LAW` | +0.0097 | +0.20 | [-0.031, +0.058] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0097 | 2.73 |
| `HP` | +0.0106 | +0.25 | [-0.015, +0.049] | [0.25, 0.5, 0.25, 0.25, 0.25] | -0.0028 | 2.83 |
| `EXVOL` | +0.0129 | +0.65 | [+0.005, +0.021] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0129 | 1.00 |
| `ELEC` | +0.0136 | +0.39 | [-0.011, +0.037] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0136 | 2.44 |
| `LEG_electrostatic` | +0.0138 | +0.39 | [-0.011, +0.037] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0138 | 2.38 |
| `DIS_MEAN` | +0.0140 | +0.58 | [+0.003, +0.024] | [0.25, 0.25, 0.25, 0.25, 0.5] | +0.0087 | 2.85 |
| `LEG_aromatic` | +0.0146 | +0.33 | [-0.007, +0.044] | [0.25, 0.5, 0.25, 0.25, 0.25] | +0.0069 | 2.14 |
| `DSSPHB` | +0.0214 | +0.43 | [-0.002, +0.042] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0214 | 2.72 |
| `LEG_coop_sheet` | +0.0221 | +0.46 | [-0.015, +0.063] | [0.5, 0.5, 0.25, 0.25, 0.25] | -0.0053 | 2.38 |
| `LEG_hbond_local` | +0.0243 | +0.48 | [+0.002, +0.045] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0243 | 2.59 |
| `LEG_contact` | +0.0254 | +0.70 | [+0.001, +0.048] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0254 | 2.83 |
| `CONTACT_LL` | +0.0261 | +0.59 | [+0.016, +0.037] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0261 | 2.63 |
| `LEG_torsion` | +0.0273 | +0.51 | [+0.022, +0.034] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0273 | 2.86 |
| `LEG` | +0.0274 | +0.58 | [-0.005, +0.063] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0274 | 2.56 |
| `CAGEO` | +0.0296 | +0.54 | [+0.008, +0.050] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0296 | 2.69 |
| `LEG_compactness` | +0.0297 | +0.53 | [-0.020, +0.095] | [0.25, 0.5, 0.25, 0.25, 0.25] | +0.0128 | 2.65 |
| `POOLGO` | +0.0298 | +0.61 | [+0.009, +0.048] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0298 | 2.59 |
| `TORS_CONS` | +0.0305 | +0.53 | [+0.012, +0.051] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0305 | 2.63 |
| `LEG_steric` | +0.0351 | +0.89 | [+0.018, +0.051] | [0.25, 0.25, 0.25, 0.25, 0.5] | +0.0323 | 1.27 |
| `AMB` | +0.0366 | +0.90 | [+0.020, +0.051] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0366 | 2.71 |
| `LEG_solvation` | +0.0370 | +0.58 | [-0.028, +0.110] | [0.5, 0.5, 0.25, 0.25, 0.25] | +0.0094 | 2.83 |
| `RAMA` | +0.0397 | +0.80 | [+0.029, +0.056] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0397 | 2.54 |
| `LEG_coop_helix` | +0.0431 | +0.77 | [+0.014, +0.076] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0431 | 2.13 |
| `ENV` | +0.0432 | +0.51 | [-0.019, +0.114] | [0.5, 1.0, 0.5, 1.0, 0.25] | -0.0045 | 2.65 |
| `DMAP_CONS` | +0.0452 | +0.73 | [+0.011, +0.069] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0452 | 2.74 |
| `CONS` | +0.0461 | +0.70 | [+0.017, +0.068] | [0.25, 0.25, 0.25, 0.25, 0.25] | +0.0461 | 2.68 |
| `CONTACT` | +0.0687 | +0.85 | [-0.002, +0.197] | [0.25, 1.0, 0.25, 0.25, 0.25] | +0.0312 | 2.62 |

SS_MATCH with the four-weight grid (0.25, 0.5, 0.75, 1.0): held-out -0.0083 A (-0.10x MDE), fold CI [-0.044, +0.044], chosen {'0': 0.5, '1': 0.75, '2': 0.5, '3': 0.5, '4': 0.5}, `best_of_k_within`: NOT A SIGNAL (split-half transfers 10% of the oracle) (k_eff 3.77).

## T5. The m-ladder (H9): mean point-cloud RMSD of the top-m average; in parentheses the paired effect against DIS at the same m (a value beyond its MDE with the fold CI excluding zero is marked *)

| m | random-m | `DIS` | `DISTPOT` | `CONS` | `DIS+DISTPOT` | `DIS+CONS` | `DIS+ENV` | `DIS+SS_MATCH` | `DIS+0.5*SS_MATCH` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 3.875 | 3.274 (+0.000) | 3.604 (+0.330) | 3.678 (+0.404) | 3.330 (+0.056) | 3.433 (+0.159) | 3.270 (-0.004) | 3.236 (-0.039) | 3.292 (+0.018) |
| 5 | 3.714 | 3.231 (+0.000) | 3.494 (+0.263) | 3.655 (+0.424) | 3.275 (+0.043) | 3.397 (+0.166) | 3.213 (-0.019) | 3.164 (-0.067) | 3.197 (-0.035) |
| 10 | 3.588 | 3.146 (+0.000) | 3.341 (+0.195) | 3.636 (+0.491) | 3.207 (+0.061) | 3.364 (+0.218) | 3.122 (-0.023) | 3.116 (-0.030) | 3.096 (-0.050) |
| 25 | 3.490 | 3.075 (+0.000) | 3.244 (+0.169) | 3.620 (+0.545) | 3.121 (+0.046) | 3.358 (+0.283) | 3.068 (-0.006) | 3.063 (-0.012) | 3.049 (-0.026) |
| 50 | 3.452 | 3.068 (+0.000) | 3.206 (+0.138) | 3.654 (+0.587) | 3.067 (-0.001) | 3.336 (+0.268) | 3.049 (-0.018) | 3.054 (-0.014) | 3.028 (-0.039) |
| 75 | 3.426 | 3.048 (+0.000) | 3.199 (+0.151) | 3.688 (+0.639) | 3.030 (-0.018) | 3.334 (+0.286) | 3.046 (-0.002) | 3.071 (+0.022) | 3.022 (-0.026) |
| 100 | 3.419 | 3.056 (+0.000) | 3.188 (+0.132) | 3.711 (+0.655) | 3.031 (-0.025) | 3.338 (+0.283) | 3.045 (-0.011) | 3.070 (+0.014) | 3.017 (-0.039) |
| 150 | 3.418 | 3.072 (+0.000) | 3.188 (+0.117) | 3.708 (+0.636) | 3.036 (-0.035) | 3.344 (+0.272) | 3.069 (-0.003) | 3.082 (+0.011) | 3.053 (-0.019) |

## T6. In-band re-ranking (H10): DIS's top-150, then 75 kept by X

| rule | mean | vs DIS top-75 | x MDE | fold CI | W/L |
|---|---:|---:|---:|---|---|
| `DIS150->SS_MATCH75` | 3.0258 | -0.0225 | -0.25 | [-0.057, +0.012] | 68/58 |
| `DIS150->(DIS+SS_MATCH)75` | 3.0263 | -0.0220 | -0.39 | [-0.035, -0.010] | 80/46 |
| `DIS150->ENV75` | 3.0460 | -0.0023 | -0.02 | [-0.110, +0.114] | 59/67 |
| `DIS150->DIS75` | 3.0483 | +0.0000 | +nan | [+0.000, +0.000] | 0/0 |
| `DIS150->(DIS+DISTPOT)75` | 3.0508 | +0.0025 | +0.05 | [-0.021, +0.035] | 62/64 |
| `DIS150->(DIS+ENV)75` | 3.0515 | +0.0031 | +0.05 | [-0.076, +0.088] | 61/63 |
| `DIS150->DISTPOT75` | 3.0535 | +0.0051 | +0.06 | [-0.040, +0.057] | 57/69 |
| `DIS150->random75` | 3.0691 | +0.0207 | +0.25 | [-0.010, +0.054] | 66/60 |
| `DIS150->(DIS+DSSPHB)75` | 3.0811 | +0.0328 | +0.44 | [-0.002, +0.061] | 51/75 |
| `DIS150->(DIS+AMB)75` | 3.0883 | +0.0400 | +0.76 | [+0.018, +0.062] | 54/72 |
| `DIS150->(DIS+CONTACT_LL)75` | 3.0973 | +0.0490 | +0.82 | [+0.021, +0.074] | 49/72 |
| `DIS150->(DIS+LEG)75` | 3.0974 | +0.0491 | +0.63 | [+0.011, +0.097] | 53/73 |
| `DIS150->(DIS+TORS_CONS)75` | 3.1038 | +0.0554 | +0.61 | [+0.012, +0.100] | 50/76 |
| `DIS150->(DIS+POOLGO)75` | 3.1119 | +0.0636 | +0.68 | [+0.011, +0.105] | 47/78 |
| `DIS150->AMB75` | 3.1140 | +0.0657 | +0.66 | [+0.034, +0.094] | 53/73 |
| `DIS150->CONTACT_LL75` | 3.1270 | +0.0786 | +0.90 | [+0.043, +0.110] | 48/78 |
| `DIS150->(DIS+CONS)75` | 3.1352 | +0.0869 | +0.81 | [+0.027, +0.134] | 49/77 |
| `DIS150->DSSPHB75` | 3.1430 | +0.0946 | +0.71 | [+0.023, +0.166] | 48/78 |
| `DIS150->LEG75` | 3.1645 | +0.1161 | +0.80 | [+0.045, +0.194] | 47/79 |
| `DIS150->TORS_CONS75` | 3.1802 | +0.1319 | +0.83 | [+0.044, +0.226] | 49/77 |
| `DIS150->POOLGO75` | 3.2035 | +0.1552 | +0.96 | [+0.041, +0.249] | 46/80 |
| `DIS150->CONS75` | 3.2320 | +0.1837 | +1.08 | [+0.056, +0.293] | 46/80 |

## T7. Non-additive forms (H11), near-corpus potentials (H12) and the SS_MATCH variants (H13)

| configuration | hypothesis | mean | vs DIS | x MDE | fold CI | folds same sign | W/L |
|---|---|---:|---:|---:|---|---|---|
| `DIS+NC_DISTPOT` | H12 | 3.0224 | -0.0259 | -0.29 | [-0.081, +0.026] | 4/5 | 64/62 |
| `DIS+0.5*NC_DISTPOT` | H12 | 3.0309 | -0.0174 | -0.27 | [-0.061, +0.022] | 4/5 | 60/66 |
| `max(DIS,ENV)` | H11 | 3.0372 | -0.0111 | -0.08 | [-0.131, +0.101] | 3/5 | 63/63 |
| `rankprod(DIS,DISTPOT)` | H11 | 3.0406 | -0.0077 | -0.07 | [-0.094, +0.091] | 3/5 | 66/60 |
| `max(DIS,DISTPOT)` | H11 | 3.0504 | +0.0020 | +0.02 | [-0.061, +0.079] | 3/5 | 64/62 |
| `max(DIS,SS_MATCH)` | H11 | 3.0636 | +0.0152 | +0.13 | [-0.031, +0.057] | 3/5 | 56/70 |
| `DIS+0.5*NC_CONTACT` | H12 | 3.0781 | +0.0297 | +0.42 | [+0.012, +0.046] | 4/5 | 64/59 |
| `rankprod(DIS,SS_MATCH)` | H11 | 3.0841 | +0.0357 | +0.31 | [-0.020, +0.090] | 4/5 | 58/68 |
| `DIS+NC_CONTACT` | H12 | 3.0947 | +0.0464 | +0.36 | [-0.012, +0.117] | 3/5 | 58/67 |
| `rankprod(DIS,ENV)` | H11 | 3.1192 | +0.0709 | +0.49 | [-0.015, +0.179] | 4/5 | 54/72 |
| `NC_DISTPOT` | H12 | 3.2120 | +0.1637 | +0.64 | [+0.055, +0.285] | 4/5 | 55/71 |
| `NC_CONTACT` | H12 | 3.7424 | +0.6941 | +2.31 | [+0.504, +0.949] | 5/5 | 40/86 |
| `DIS+0.75*SS_MATCH2` | H13 | 3.0327 | -0.0156 | -0.16 | [-0.064, +0.042] | 3/5 | 65/61 |
| `DIS+0.75*SS_MATCH` | H13 | 3.0376 | -0.0107 | -0.11 | [-0.057, +0.042] | 3/5 | 66/60 |
| `DIS+0.5*SS_MATCH2` | H13 | 3.0380 | -0.0104 | -0.16 | [-0.042, +0.024] | 3/5 | 63/62 |
| `DIS+0.25*SS_MATCH2` | H13 | 3.0422 | -0.0062 | -0.19 | [-0.030, +0.014] | 4/5 | 67/58 |
| `DIS+gate(SS_MATCH)` | H13 | 3.0455 | -0.0028 | -0.07 | [-0.018, +0.011] | 2/5 | 71/54 |
| `DIS+1.0*SS_MATCH2` | H13 | 3.0593 | +0.0110 | +0.08 | [-0.070, +0.099] | 3/5 | 64/62 |
| `DIS+SS_MATCH+SS_MATCH2` | H13 | 3.0598 | +0.0115 | +0.10 | [-0.048, +0.110] | 1/5 | 63/63 |
| `SS_MATCH2` | H13 | 4.4167 | +1.3684 | +3.75 | [+1.234, +1.481] | 5/5 | 22/104 |

## T8. Complementarity and redundancy: partial Spearman of each channel with the ORACLE candidate RMSD given DIS (positive = adds correct ranking information), and its rank correlation with DIS

| channel | partial rho given DIS | rho with DIS | DIS+X effect on the top-75 average |
|---|---:|---:|---:|
| `CONS` | +0.246 | +0.494 | +0.2858 |
| `TORS_CONS` | +0.122 | +0.419 | +0.2626 |
| `DMAP_CONS` | +0.119 | +0.518 | +0.3559 |
| `POOLGO` | +0.113 | +0.403 | +0.2993 |
| `DISTPOT` | +0.090 | +0.426 | -0.0185 |
| `CONTACT_LL` | +0.058 | +0.589 | +0.0884 |
| `CAGEO` | +0.048 | +0.350 | +0.2902 |
| `CONTACT` | +0.039 | +0.121 | +0.0450 |
| `ENV` | +0.038 | +0.268 | -0.0021 |
| `RG_LAW` | +0.033 | +0.329 | +0.1106 |
| `DSSPHB` | +0.022 | +0.309 | +0.2000 |
| `LEG` | +0.018 | +0.367 | +0.1647 |
| `SS_MATCH` | +0.017 | +0.012 | +0.0222 |
| `RAMA` | +0.016 | +0.299 | +0.2207 |
| `HP` | +0.013 | +0.023 | +0.0465 |
| `DIS` | +0.000 | +1.000 |  |
| `EXVOL` | -0.010 | +0.056 | +0.0129 |
| `ELEC` | -0.014 | -0.049 | +0.0593 |
| `AMB` | -0.018 | -0.019 | +0.0813 |

## T9. The genuine CVaR-VQE arm (9 qubits, alpha 0.18, T 0.5, 80 Adam iterations, exact parameter-shift gradient), point cloud: each configuration's VQE selection against DIS's VQE selection, seed 0, with the seed-1 replication, the realised tail m, the set-equality gate, and eps = RMSD_VQE - RMSD_top-m

| configuration | VQE s0 | vs DIS-VQE | x MDE | fold CI | VQE s1 | effect s1 | s1 inside s0's CI | m | max eps | entropy bits | ESS | gate |
|---|---:|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---|
| `DIS+0.5*SS_MATCH` | 3.0214 | -0.0366 | -0.47 | [-0.066, -0.010] | 3.0289 | -0.0291 | yes | 73.6 | 4.1e-14 | 8.81 | 58.5 | pass |
| `DIS+DISTPOT` | 3.0302 | -0.0279 | -0.25 | [-0.094, +0.056] | 3.0404 | -0.0176 | yes | 73.2 | 8.7e-14 | 8.81 | 57.9 | pass |
| `DIS+0.25*LEG_hbond_longrange` | 3.0314 | -0.0267 | -0.69 | [-0.051, +0.011] | 3.0324 | -0.0256 | yes | 74.6 | 3.9e-14 | 8.84 | 60.7 | pass |
| `DIS+0.25*SS_MATCH` | 3.0321 | -0.0259 | -0.56 | [-0.049, -0.005] | 3.0463 | -0.0117 | yes | 74.4 | 6.5e-14 | 8.83 | 60.5 | pass |
| `REJ[POOLGO]->DIS` | 3.0336 | -0.0244 | -0.51 | [-0.045, -0.004] | 3.0307 | -0.0274 | yes | 70.6 | 1.5e-13 | 8.78 | 53.6 | pass |
| `REJ[RG_LAW]->DIS` | 3.0366 | -0.0215 | -0.42 | [-0.042, +0.002] | 3.0421 | -0.0159 | yes | 70.5 | 1.5e-13 | 8.78 | 53.5 | pass |
| `DIS+0.5*DISTPOT` | 3.0402 | -0.0178 | -0.25 | [-0.054, +0.032] | 3.0455 | -0.0125 | yes | 74.3 | 5.7e-14 | 8.84 | 60.1 | pass |
| `REJ[LEG_compactness]->DIS` | 3.0421 | -0.0159 | -0.33 | [-0.037, +0.004] | 3.0399 | -0.0181 | yes | 70.4 | 1.5e-13 | 8.78 | 53.7 | pass |
| `DIS+ENV` | 3.0485 | -0.0095 | -0.07 | [-0.134, +0.099] | 3.0465 | -0.0115 | yes | 73.9 | 8.4e-14 | 8.83 | 59.4 | pass |
| `REJ[random]->DIS` | 3.0500 | -0.0080 | -0.32 | [-0.021, +0.003] | 3.0521 | -0.0059 | yes | 70.1 | 7.2e-14 | 8.78 | 53.6 | pass |
| `DIS+0.5*ENV` | 3.0502 | -0.0078 | -0.09 | [-0.085, +0.083] | 3.0369 | -0.0211 | yes | 73.5 | 2.2e-14 | 8.82 | 58.6 | pass |
| `DIS` | 3.0580 | +0.0000 | +nan | [+0.000, +0.000] | 3.0562 | -0.0018 | no | 74.1 | 1.1e-13 | 8.84 | 60.2 | pass |
| `DIS+STAT_COMP` | 3.0598 | +0.0017 | +0.02 | [-0.055, +0.079] | 3.0608 | +0.0028 | yes | 73.8 | 6.2e-14 | 8.82 | 59.0 | pass |
| `DIS+DIS_MEAN` | 3.0633 | +0.0053 | +0.13 | [-0.015, +0.026] | 3.0669 | +0.0089 | yes | 74.1 | 4.7e-02 | 8.83 | 59.4 | pass |
| `DIS+ROSETTA_LIKE` | 3.0635 | +0.0055 | +0.07 | [-0.057, +0.099] | 3.0591 | +0.0011 | yes | 75.1 | 2.9e-14 | 8.85 | 61.9 | pass |
| `DIS+EXVOL` | 3.0641 | +0.0061 | +0.35 | [-0.006, +0.015] | 3.0630 | +0.0050 | yes | 73.7 | 1.1e-13 | 8.82 | 59.1 | pass |
| `DIS+LEG_coop_sheet` | 3.0662 | +0.0082 | +0.10 | [-0.052, +0.075] | 3.0598 | +0.0018 | yes | 75.4 | 3.7e-14 | 8.86 | 62.3 | pass |
| `DIS+SS_MATCH` | 3.0731 | +0.0151 | +0.13 | [-0.033, +0.097] | 3.0707 | +0.0127 | yes | 72.9 | 2.4e-14 | 8.81 | 57.9 | pass |
| `DIS_MEAN` | 3.0769 | +0.0188 | +0.31 | [-0.024, +0.055] | 3.0826 | +0.0246 | yes | 74.5 | 3.8e-14 | 8.84 | 60.5 | pass |
| `DIS+LEG_steric` | 3.0821 | +0.0241 | +0.50 | [+0.003, +0.044] | 3.0739 | +0.0159 | yes | 73.3 | 2.9e-14 | 8.82 | 58.9 | pass |
| `DIS+PHYSICS_COMP` | 3.0843 | +0.0263 | +0.41 | [-0.019, +0.067] | 3.0793 | +0.0213 | yes | 73.3 | 3.8e-14 | 8.82 | 58.4 | pass |
| `DIS+LEG_aromatic` | 3.0896 | +0.0316 | +0.24 | [-0.015, +0.094] | 3.0793 | +0.0213 | yes | 73.5 | 2.8e-02 | 8.81 | 58.3 | pass |
| `DIS+HP` | 3.0956 | +0.0376 | +0.28 | [-0.052, +0.160] | 3.0960 | +0.0380 | yes | 73.8 | 5.0e-14 | 8.82 | 59.9 | pass |
| `DIS+CONTACT` | 3.0957 | +0.0377 | +0.27 | [-0.058, +0.185] | 3.0969 | +0.0388 | yes | 73.8 | 2.6e-14 | 8.83 | 59.0 | pass |
| `DIS+ELEC` | 3.1092 | +0.0512 | +0.43 | [-0.016, +0.122] | 3.1017 | +0.0437 | yes | 74.2 | 1.3e-01 | 8.83 | 60.3 | pass |
| `DIS+LEG_electrostatic` | 3.1143 | +0.0563 | +0.47 | [-0.015, +0.126] | 3.1001 | +0.0420 | yes | 74.2 | 1.8e-02 | 8.83 | 60.3 | pass |
| `DIS+LEG_solvation` | 3.1248 | +0.0667 | +0.47 | [-0.042, +0.188] | 3.1234 | +0.0654 | yes | 74.6 | 4.6e-02 | 8.83 | 60.4 | pass |
| `DIS+CONTACT_LL` | 3.1300 | +0.0720 | +0.91 | [+0.044, +0.097] | 3.1348 | +0.0767 | yes | 74.0 | 1.0e-13 | 8.84 | 60.6 | pass |
| `DIS+AMB` | 3.1317 | +0.0737 | +0.58 | [+0.038, +0.111] | 3.1228 | +0.0647 | yes | 74.0 | 6.5e-02 | 8.83 | 59.9 | pass |
| `DIS+LEG_contact` | 3.1368 | +0.0788 | +0.67 | [-0.006, +0.154] | 3.1371 | +0.0790 | yes | 73.5 | 6.0e-14 | 8.82 | 58.5 | pass |
| `DIS+LEG_compactness` | 3.1383 | +0.0802 | +0.44 | [-0.041, +0.245] | 3.1343 | +0.0763 | yes | 73.6 | 2.8e-14 | 8.81 | 59.2 | pass |
| `DIS+RG_LAW` | 3.1559 | +0.0979 | +0.49 | [+0.009, +0.238] | 3.1575 | +0.0995 | yes | 74.3 | 3.6e-01 | 8.83 | 59.5 | pass |
| `DIS+LEG_torsion` | 3.1793 | +0.1213 | +0.81 | [+0.014, +0.209] | 3.1864 | +0.1284 | yes | 73.5 | 2.4e-02 | 8.83 | 59.6 | pass |
| `DIS+RG_UNIV` | 3.1876 | +0.1296 | +0.70 | [+0.030, +0.224] | 3.1802 | +0.1222 | yes | 74.5 | 4.6e-02 | 8.84 | 60.1 | pass |
| `STAT_COMP` | 3.2011 | +0.1431 | +0.58 | [-0.031, +0.333] | 3.2104 | +0.1524 | yes | 73.1 | 6.3e-04 | 8.81 | 58.6 | pass |
| `DISTPOT` | 3.2031 | +0.1451 | +0.60 | [+0.021, +0.305] | 3.2161 | +0.1580 | yes | 73.8 | 4.0e-02 | 8.82 | 59.4 | pass |
| `DIS+LEG` | 3.2151 | +0.1571 | +0.85 | [+0.041, +0.227] | 3.2169 | +0.1589 | yes | 75.5 | 6.9e-02 | 8.84 | 61.9 | pass |
| `DIS+LEG_hbond_longrange` | 3.2441 | +0.1861 | +0.83 | [+0.017, +0.359] | 3.2406 | +0.1825 | yes | 75.2 | 1.5e-14 | 8.84 | 61.3 | pass |
| `DIS+DSSPHB` | 3.2543 | +0.1963 | +1.11 | [+0.072, +0.290] | 3.2544 | +0.1963 | yes | 75.1 | 6.1e-14 | 8.83 | 60.5 | pass |
| `DIS+RAMA` | 3.2608 | +0.2028 | +1.12 | [+0.062, +0.340] | 3.2602 | +0.2022 | yes | 73.7 | 6.0e-02 | 8.82 | 59.3 | pass |
| `DIS+TORS_CONS` | 3.3111 | +0.2530 | +1.16 | [+0.090, +0.392] | 3.3048 | +0.2468 | yes | 73.4 | 3.1e-02 | 8.82 | 58.0 | pass |
| `DIS+LEG_hbond_local` | 3.3127 | +0.2547 | +1.25 | [+0.079, +0.383] | 3.3140 | +0.2560 | yes | 75.1 | 9.4e-14 | 8.84 | 61.3 | pass |
| `EXVOL` | 3.3264 | +0.2684 | +1.15 | [+0.109, +0.395] | 3.3246 | +0.2666 | yes | 92.9 | 1.3e+00 | 9.00 | 92.0 | pass |
| `DIS+CAGEO` | 3.3302 | +0.2722 | +1.26 | [+0.110, +0.434] | 3.3326 | +0.2746 | yes | 74.8 | 1.0e-13 | 8.84 | 61.2 | pass |
| `DIS+CONS` | 3.3339 | +0.2759 | +1.23 | [+0.076, +0.457] | 3.3297 | +0.2716 | yes | 74.1 | 3.8e-02 | 8.83 | 59.5 | pass |
| `DIS+CONSIST_COMP` | 3.3404 | +0.2824 | +1.27 | [+0.095, +0.442] | 3.3351 | +0.2771 | yes | 73.4 | 5.4e-04 | 8.81 | 58.8 | pass |
| `DIS+LEG_coop_helix` | 3.3435 | +0.2854 | +1.47 | [+0.175, +0.384] | 3.3367 | +0.2787 | yes | 73.8 | 1.1e-13 | 8.82 | 59.6 | pass |
| `DIS+POOLGO` | 3.3451 | +0.2871 | +1.31 | [+0.104, +0.428] | 3.3494 | +0.2914 | yes | 74.3 | 9.6e-14 | 8.83 | 60.0 | pass |
| `ROSETTA_LIKE` | 3.3570 | +0.2990 | +0.98 | [+0.026, +0.659] | 3.3561 | +0.2980 | yes | 75.0 | 2.1e-14 | 8.83 | 61.1 | pass |
| `LEG_steric` | 3.3692 | +0.3111 | +1.33 | [+0.136, +0.426] | 3.3709 | +0.3129 | yes | 92.8 | 1.2e+00 | 9.00 | 92.1 | pass |
| `DIS+DMAP_CONS` | 3.3990 | +0.3410 | +1.39 | [+0.116, +0.506] | 3.4045 | +0.3465 | yes | 74.9 | 7.9e-02 | 8.84 | 61.0 | pass |
| `LEG_aromatic` | 3.4676 | +0.4096 | +1.47 | [+0.308, +0.597] | 3.4663 | +0.4083 | yes | 82.9 | 1.4e+00 | 8.90 | 74.6 | pass |
| `ENV` | 3.5133 | +0.4553 | +1.39 | [+0.344, +0.591] | 3.5099 | +0.4519 | yes | 71.0 | 3.6e-01 | 8.78 | 55.9 | pass |
| `RG_UNIV` | 3.5341 | +0.4761 | +1.51 | [+0.256, +0.637] | 3.5527 | +0.4947 | yes | 73.8 | 2.1e-14 | 8.83 | 59.6 | pass |
| `PHYSICS_COMP` | 3.5487 | +0.4907 | +1.75 | [+0.317, +0.666] | 3.5455 | +0.4875 | yes | 73.9 | 1.5e-14 | 8.82 | 59.9 | pass |
| `SS_MATCH` | 3.5617 | +0.5037 | +1.82 | [+0.388, +0.626] | 3.5578 | +0.4998 | yes | 73.8 | 4.7e-01 | 8.82 | 59.5 | pass |
| `CONTACT_LL` | 3.5743 | +0.5163 | +1.89 | [+0.359, +0.638] | 3.5850 | +0.5270 | yes | 73.4 | 4.7e-01 | 8.76 | 59.8 | pass |
| `CONTACT` | 3.5830 | +0.5250 | +1.77 | [+0.306, +0.789] | 3.6106 | +0.5526 | yes | 70.6 | 6.3e-02 | 8.77 | 55.0 | pass |
| `RG_LAW` | 3.5960 | +0.5380 | +1.67 | [+0.333, +0.826] | 3.5908 | +0.5328 | yes | 74.3 | 1.5e-14 | 8.82 | 59.3 | pass |
| `LEG_solvation` | 3.6253 | +0.5673 | +1.82 | [+0.322, +0.879] | 3.6477 | +0.5897 | yes | 73.8 | 1.0e-14 | 8.82 | 59.5 | pass |
| `LEG_compactness` | 3.6752 | +0.6172 | +1.73 | [+0.346, +0.933] | 3.6759 | +0.6178 | yes | 85.0 | 1.3e+00 | 8.93 | 78.5 | pass |
| `LEG_torsion` | 3.6867 | +0.6287 | +1.97 | [+0.448, +0.804] | 3.6975 | +0.6395 | yes | 73.8 | 8.2e-14 | 8.82 | 59.3 | pass |
| `CONS` | 3.6929 | +0.6349 | +1.99 | [+0.438, +0.806] | 3.6876 | +0.6296 | yes | 74.3 | 3.2e-14 | 8.83 | 60.5 | pass |
| `LEG` | 3.7553 | +0.6973 | +2.10 | [+0.487, +0.821] | 3.7544 | +0.6964 | yes | 74.5 | 2.6e-14 | 8.83 | 60.1 | pass |
| `DSSPHB` | 3.7898 | +0.7318 | +2.12 | [+0.536, +0.873] | 3.7917 | +0.7336 | yes | 73.9 | 1.9e-14 | 8.82 | 59.9 | pass |
| `LEG_coop_helix` | 3.7947 | +0.7367 | +2.19 | [+0.526, +0.876] | 3.7955 | +0.7375 | yes | 92.2 | 1.1e+00 | 8.99 | 91.0 | pass |
| `LEG_electrostatic` | 3.8187 | +0.7607 | +2.34 | [+0.593, +0.920] | 3.8304 | +0.7724 | yes | 76.5 | 1.1e+00 | 8.85 | 63.5 | pass |
| `TORS_CONS` | 3.8310 | +0.7729 | +2.07 | [+0.602, +0.941] | 3.8309 | +0.7729 | yes | 73.9 | 6.3e-14 | 8.82 | 59.2 | pass |
| `ELEC` | 3.8314 | +0.7734 | +2.39 | [+0.610, +0.926] | 3.8361 | +0.7781 | yes | 76.4 | 1.1e+00 | 8.85 | 63.3 | pass |
| `HP` | 3.8469 | +0.7889 | +2.14 | [+0.573, +1.026] | 3.8494 | +0.7914 | yes | 73.9 | 6.9e-02 | 8.83 | 59.4 | pass |
| `DMAP_CONS` | 3.8551 | +0.7971 | +2.08 | [+0.565, +0.994] | 3.8583 | +0.8003 | yes | 74.2 | 1.1e-13 | 8.83 | 59.9 | pass |
| `AMB` | 3.8805 | +0.8225 | +2.63 | [+0.684, +0.953] | 3.9147 | +0.8567 | yes | 73.2 | 1.5e-14 | 8.83 | 59.4 | pass |
| `CONSIST_COMP` | 3.9049 | +0.8468 | +2.21 | [+0.641, +1.026] | 3.9053 | +0.8473 | yes | 72.7 | 1.7e-03 | 8.80 | 56.7 | pass |
| `RAMA` | 3.9260 | +0.8680 | +2.27 | [+0.659, +1.033] | 3.9375 | +0.8795 | yes | 74.3 | 5.4e-14 | 8.83 | 59.7 | pass |
| `LEG_hbond_local` | 3.9289 | +0.8709 | +2.24 | [+0.652, +1.053] | 3.9283 | +0.8703 | yes | 74.6 | 5.1e-14 | 8.84 | 61.0 | pass |
| `LEG_hbond_longrange` | 3.9593 | +0.9013 | +2.21 | [+0.692, +1.131] | 3.9421 | +0.8840 | yes | 37.2 | 3.5e-02 | 8.11 | 16.6 | pass |
| `CAGEO` | 3.9700 | +0.9120 | +2.29 | [+0.700, +1.088] | 3.9691 | +0.9111 | yes | 73.7 | 2.8e-04 | 8.83 | 58.9 | pass |
| `POOLGO` | 3.9715 | +0.9135 | +2.31 | [+0.710, +1.074] | 3.9673 | +0.9093 | yes | 68.0 | 6.6e-02 | 8.70 | 54.7 | pass |
| `LEG_contact` | 4.4045 | +1.3465 | +3.86 | [+1.241, +1.499] | 4.4397 | +1.3817 | yes | 75.2 | 3.2e-02 | 8.85 | 61.8 | pass |
| `LEG_coop_sheet` | 4.9899 | +1.9319 | +3.20 | [+1.594, +2.273] | 4.9633 | +1.9052 | yes | 14.6 | 1.3e+00 | 7.17 | 9.6 | pass |

## T10. The built-chain endpoint (the production projection of the same top-75 averages) against DIS, with the same sets on the point cloud beside it

| configuration | built chain | vs DIS | x MDE | fold CI | W/L | point cloud | vs DIS | x MDE |
|---|---:|---:|---:|---|---|---:|---:|---:|
| `REJ[POOLGO]->DIS` | 3.1958 | -0.0168 | -0.34 | [-0.037, +0.005] | 24/23 | 3.0325 | -0.0159 | -0.37 |
| `DIS+0.25*LEG_hbond_longrange` | 3.1987 | -0.0139 | -0.33 | [-0.037, +0.006] | 57/59 | 3.0354 | -0.0130 | -0.39 |
| `DIS+DISTPOT` | 3.2100 | -0.0026 | -0.02 | [-0.067, +0.071] | 60/66 | 3.0299 | -0.0185 | -0.17 |
| `DIS` | 3.2126 | +0.0000 | +nan | [+0.000, +0.000] | 0/0 | 3.0483 | +0.0000 | +nan |
| `DIS+0.5*SS_MATCH` | 3.2129 | +0.0003 | +0.00 | [-0.032, +0.030] | 66/60 | 3.0223 | -0.0261 | -0.35 |
| `DIS+ENV` | 3.2361 | +0.0234 | +0.16 | [-0.102, +0.126] | 55/71 | 3.0462 | -0.0021 | -0.02 |
| `DIS_MEAN` | 3.2576 | +0.0450 | +0.68 | [-0.020, +0.105] | 51/75 | 3.0784 | +0.0301 | +0.51 |
| `DISTPOT` | 3.3840 | +0.1714 | +0.70 | [+0.064, +0.304] | 51/75 | 3.1990 | +0.1507 | +0.63 |
| `DIS+CONS` | 3.4619 | +0.2493 | +1.06 | [+0.030, +0.449] | 46/80 | 3.3341 | +0.2858 | +1.29 |
| `CONS` | 3.7656 | +0.5530 | +1.68 | [+0.346, +0.745] | 39/87 | 3.6878 | +0.6394 | +2.02 |

## T11. Trainability (H8): gradient variance of the deployed objective at n = 9, depth 3, alpha 0.18, T 0.5, 120 theta draws, median over 12 targets, under the deployed rank currency and under a bounded standardisation that keeps the channel's own gaps

| channel | Var[dF/dtheta_0], zrank | Var, asinh-MAD | ratio to DIS (asinh) |
|---|---:|---:|---:|
| `DIS` | 3.0505e-02 | 2.5505e-02 | 1.00 |
| `LEG` | 2.9124e-02 | 2.5471e-02 | 1.00 |
| `AMB` | 3.2686e-02 | 2.1941e-02 | 0.86 |
| `CONS` | 3.3257e-02 | 2.4112e-02 | 0.95 |
| `DSSPHB` | 3.0717e-02 | 2.7628e-02 | 1.08 |
| `CONTACT` | 2.7248e-02 | 2.9066e-02 | 1.14 |
| `DISTPOT` | 2.6676e-02 | 2.4774e-02 | 0.97 |
| `RAMA` | 3.0490e-02 | 2.7493e-02 | 1.08 |
| `CAGEO` | 2.8550e-02 | 2.7159e-02 | 1.06 |
| `ENV` | 2.7258e-02 | 2.3503e-02 | 0.92 |
