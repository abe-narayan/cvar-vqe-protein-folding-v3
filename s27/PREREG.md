# S27 PRE-REGISTRATION -- ALTERNATIVE HAMILTONIANS FOR THE CANDIDATE-POOL SELECTOR

Written 2026-09-14 before any endpoint number was read. Sprint 27 asks one question: can a
different Hamiltonian, energy model, or principled combination give a meaningfully better final
RMSD than the shipped distogram Bayes-risk score, on the established 126-target instrument, with
the production readout and the genuine CVaR-VQE selector left intact?

## 1. What is held fixed (identical to the S25 seven-configuration suite, `s25/PREREG_PHYS.md`)

| | value |
|---|---|
| targets | the 126 dev targets, `sorted(pdb)`, all of them, no cherry-picking |
| folds | the 5 pinned folds (`s24.stats_lib.pinned_folds`) |
| pool | the shipped K = 500 BLOSUM pool per target (`s24.d_harness.Candidates.from_universe`) |
| selector | (a) classical top-75 of the energy; (b) the GENUINE CVaR-VQE of `s24.d_harness.arm_vqe` (exact `StatevectorCircuit`, 9 qubits over 500 + 12 padding, 3 layers, 80 Adam iterations on the exact parameter-shift gradient, alpha = 0.18, T = 0.5, seed 0), whose tail is read exactly |
| readout | the deployed uniform coordinate average in the retained set's own medoid frame |
| primary basis | POINT CLOUD (`rmsd_of_set`, anchored: the shipped score's top-75 = 3.048338 A) |
| secondary basis | BUILT CHAIN through the production projection `s12.instrument.project` (anchored 3.2126 on the rebuild basis, L57) for the finalists |
| statistics | `s24.stats_lib.compare`: paired per target, SE, MDE = 2.8016 x SE, iid CI beside the fold-clustered CI, W/L/ties, median, concentration null; the decision is on the fold-clustered CI |
| ties | never broken by array order: a stable per-target random key breaks exact ties in every top-m; the tie fraction is reported |

Nothing native enters any channel. Natives are read only in `oracle_*` functions for the endpoint
and the labelled diagnostics.

## 2. The channels (`s27/ham_lib.py`)

Eighteen new native-free channels in six families (geometric: RG_LAW, RG_UNIV, EXVOL, CAGEO;
torsion: RAMA; statistical pair / burial fitted per target on the leakage-safe universe:
CONTACT, DISTPOT, ENV, HP; backbone physics on the ideal rebuild: DSSPHB, ELEC; pool consistency:
CONS, DMAP_CONS, TORS_CONS, POOLGO; compatibility: SS_MATCH; distogram re-readings: CONTACT_LL,
DIS_MEAN), the eleven Legacy terms one at a time (LEG_*), and the three S25 references (DIS, LEG,
AMB from their caches). Every channel is a vector over the 500 candidates, lower is better.

## 3. Hypotheses, in the order they are tested, with the falsifier and the registered prior

H1 (singles). Some channel, used alone as the selector, beats the shipped score.
  Falsifier: no single channel's top-75 point-cloud RMSD is below DIS's beyond its own MDE with
  the fold CI excluding zero on 5/5 folds. Registered prior: EVERY single is null-to-WORSE.
  Reason: S25 measured the two genuine physics energies at +0.33 / +0.46 A worse than a random
  subset; S12 to S23 closed within-pool ranking at four levels; the consistency channels select
  the most typical 75, which S23 L5 measured as costing +0.142 A on the average. The one family
  with no record is the universe-fitted statistical potentials (CONTACT, DISTPOT, ENV, CAGEO);
  the prior is null, because the pool is drawn from the same corpus they are fitted on, so they
  measure corpus typicality, which is what CONS measures.

H2 (complements). Some channel, added to DIS at equal weight in the S25 rank-standardised
  currency E = zrank(zrank(DIS) + zrank(X)), beats DIS.
  Falsifier: no equal-weight pair beats DIS beyond its MDE with the fold CI excluding zero on
  5/5 folds. Controls, matched in the operator's space: the same pair with X's marginal
  rank-permuted (correspondence destroyed, amplitude kept), and the random-75 null (16 draws).
  Registered prior: null for every pair (S24 L16 and S25 L18 closed the additive form for LEG and
  AMB; the re-readings are closed by S25 L12); the pair whose failure would surprise me least is
  DIS + DSSPHB, and the one whose success would surprise me least is DIS + CONS, because
  consensus is the only in-band discriminator on record (S12).

H3 (regulariser weight). A channel helps only at a small weight. Weights w in {0.25, 0.5} on
  zrank(X) beside zrank(DIS); the weight is CHOSEN leave-fold-out (nested) and the held-out
  effect is what is reported; the full-leakage best is reported beside `best_of_k_within`.
  Falsifier as H2.

H4 (staged). A channel is useful as a REJECT before DIS rather than as a score: reject the
  worst q = 10% of the pool by X, then take DIS's top-75 of the rest. Controls: reject a random
  10%. Falsifier as H2. Prior: null; S26 L43/L86 measured AMBER's reject as harmful at every
  threshold and S24 D1-C the filter form as harmful.

H5 (triples and a composite). The best two complements together, and a Rosetta-centroid-like
  composite (ENV + CONTACT + EXVOL + RG_LAW) and a physics composite (DSSPHB + ELEC + HP +
  EXVOL), each alone and each with DIS. Same falsifier. Prior: null.

H6 (adaptive mixture). A native-free gate: the channel's weight rises with the distogram's mean
  per-pair entropy (a diffuse posterior leans more on X). The map is FIXED a priori
  (w_t = 0.5 * H_t / max_t H_t), not fitted. Same falsifier. Prior: null (routers fail, S22 L7,
  S23 L7, S26 L115).

H7 (the CVaR-VQE arm). For DIS and the finalists (the best four singles and the best four
  combinations by the primary endpoint, plus any that clears its falsifier), the genuine
  CVaR-VQE selection of the same energy, with the set-equality gate, the m-ladder decomposition
  RMSD_VQE = RMSD_top75 + [RMSD_topm - RMSD_top75] + eps, entropy and ESS, and a second seed.
  Falsifier for "the Hamiltonian gives the selector a better objective": the VQE arm on X beats
  the VQE arm on DIS beyond its MDE with the fold CI excluding zero on 5/5 folds, on BOTH seeds.
  Prior: the VQE arm tracks the classical top-m of the same energy to machine precision (S25 L16,
  1764 cells), so H7 reduces to H1 to H6 up to the m-ladder.

H8 (trainability). Gradient variance at n = 9, depth 3, theta ~ N(0, 0.6^2), exact
  parameter-shift, for E under (i) zrank, which makes every channel the same ladder by
  construction, and (ii) the bounded asinh-MAD standardisation, which keeps target-dependent
  gaps; for DIS, LEG, AMB and the finalists. Falsifier for "a Hamiltonian changes trainability":
  the per-channel variance at alpha = 0.18 differs from DIS's by more than 2x under (ii). Prior:
  no difference under (i), possible differences under (ii) that do not reach the endpoint.

## 4. Diagnostics (ORACLE where they read the native; labelled at every appearance)

Per channel: Spearman rho with the ORACLE candidate RMSD over the whole pool and in-band
(< 3 A); rho with DIS; top-75 overlap with DIS's top-75; energy distribution (decades, tie
fraction, top-10 variance share); the channel-by-channel rank-correlation matrix (redundancy);
per-target and per-fold effects; cost per evaluation. Complementarity: partial rho of X with the
ORACLE RMSD given DIS.

## 5. What counts as a result

A candidate is PROMISING only if it clears its falsifier on the point cloud AND the built-chain
endpoint agrees in sign beyond 0.7x its MDE; NULL if inside 0.7x MDE with a stated power; WORSE
if the fold CI is above zero; UNDERPOWERED if the MDE exceeds the effect the prior predicts a
real lever would show (0.05 A); ORACLE-ONLY if its only positive number reads the native. Type-M
zone (0.7 to 1.3x MDE) is not a result. Every grid choice (weights, thresholds) is priced with
`best_of_k_within`. No native RMSD chooses a weight or a threshold.

## 6. Rule 0 forks, the alternative not taken

| fork | taken | not taken |
|---|---|---|
| where statistical potentials are fitted | the target's own leakage-safe universe (13k windows; peptides of other folds + identity-filtered fragments), pseudocount 1 | the whole peptide DB (leaks the target's native); published MJ tables (not reproducible from memory to the digit) |
| currency for combinations | rank standardisation (S25's, monotone, equal marginals) | raw moment (non-monotone in float64 on AMBER, S25 L16) |
| readout | the deployed uniform top-75 average | medoid, weighted, m-tuned readouts (closed S23 L5/L8) |
| basis | point cloud primary, built chain secondary, both named | mixing bases |
| null | random-75 and rank-permuted controls per target | the pinned constant alone |

## ADDENDUM 1 (2026-09-14, written after H1 to H6's pool results, before any wave-2 number)

Wave 1 found no configuration beyond its MDE, and one mechanism worth testing directly: the
consistency channels carry the LARGEST partial rank correlation with the ORACLE candidate RMSD
given DIS (CONS +0.25, TORS_CONS +0.12, DMAP_CONS +0.12, POOLGO +0.11;
`s27/results/redundancy.json`) and are the MOST harmful in the top-75 average (DIS+CONS +0.29 A,
worse than its own permuted control by +0.27 at 1.7x MDE). Ranking skill and set-average value
are decoupled: a channel that finds near-native candidates also finds SIMILAR candidates, and
the readout averages them. Two pre-registered follow-ups:

H9 (the m-ladder). For DIS, DISTPOT, CONS, DIS+DISTPOT, DIS+CONS, DIS+ENV, DIS+SS_MATCH and
  DIS+0.5*SS_MATCH: the top-m point-cloud average at m in {3, 5, 10, 25, 50, 75, 100, 150},
  paired against DIS at the SAME m, and against random m-subsets (8 draws) at the same m.
  Prediction: at small m (<= 10) DIS+CONS beats DIS (ranking skill shows where averaging cannot
  hide it) but no small-m arm beats DIS at m = 75, so no arm is a production improvement;
  the per-arm minimum over m is priced with `best_of_k_within` and never quoted as a gain.
  Falsifier for "a better ranker exists": some arm beats DIS at the same m beyond its MDE with
  the fold CI excluding zero at two adjacent m values.

H10 (in-band re-ranking). Take DIS's top-150, then keep the 75 best by X (X in DISTPOT, CONS,
  CONTACT_LL, ENV, DSSPHB, POOLGO, TORS_CONS, LEG, AMB); controls: keep 75 at random from the
  150 (8 draws), and keep DIS's own top-75 (the incumbent). This asks whether any channel
  discriminates INSIDE the score's near-native band, the one place the record says the gap
  lives. Falsifier as H2. Prior: null-to-worse (S17 L23: in-band discrimination is
  signal-limited; S25 §1.5).

H11 (non-additive forms). For the two least-harmful complements (DISTPOT, ENV) and SS_MATCH:
  max-rank E = max(zrank(DIS), zrank(X)) (a candidate must be good under both) and rank-product
  E = zrank(rank(DIS) * rank(X)). Falsifier as H2. Prior: null.

H12 (near-corpus statistical potential). DISTPOT and CONTACT refitted on the target's 2,000
  most BLOSUM-similar universe windows instead of all 13k (a sequence-conditioned statistical
  potential). Alone, with DIS at equal weight, and at 0.5. Falsifier as H1/H2. Prior: null.

## ADDENDUM 2 (2026-09-14, written after wave 2, before any wave-3 number)

Wave 2 left one channel with a consistent, sub-MDE negative effect in every form tried
(SS_MATCH at half weight: -0.026 A at m = 75, 0.35x MDE, the same weight chosen on all five
folds; -0.02 to -0.05 A across m = 5..100; -0.022 as an in-band re-ranker). Its partial rank
correlation with the ORACLE RMSD is +0.017 and its rank correlation with DIS is 0.012: an
orthogonal, weak signal. Wave 3 asks whether a better-founded version of the same idea
(sequence-conditioned secondary-structure compatibility) carries more of it:

H13 (SS_MATCH variants). (a) SS_MATCH2: the secondary-structure call from Kabsch-Sander
  H-bonds on the ideal rebuild (helix = i -> i+4 bonded on two consecutive residues; strand =
  any |i-j| >= 3 backbone H-bond) against per-residue helix / strand propensities FITTED on
  the target's leakage-safe universe (the fraction of universe windows in which that residue
  type sits in a helix / strand call), instead of Chou-Fasman; (b) SS_MATCH at w = 0.75 to
  complete the weight grid (0.25, 0.5, 0.75, 1.0), nested choice re-run, `best_of_k_within`
  over the four; (c) a FIXED a-priori gate: w_t = 0.5 * c_t / max_t c_t where c_t is the mean
  absolute Chou-Fasman helix-minus-strand contrast of the sequence (a sequence with a strong
  propensity leans more on the channel). Falsifier as H2. Prior: null at the MDE; if (a) is
  larger than SS_MATCH by more than its MDE the mechanism is the H-bond call, if not it is the
  propensity table.
