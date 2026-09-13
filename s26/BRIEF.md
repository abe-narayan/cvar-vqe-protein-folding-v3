# S26 BRIEF: the sprint's questions and their pre-registered falsifiers (lane E)

Written 2026-09-13 from `s26/EXAMINATION.md`, the lane pre-registrations on disk
(`s26/PREREG_C1..C5.md`, `PREREG_B1..B3.md`, `PREREG_c3_control.md`, `PREREG_amber_reject.md`,
`PREREG_cis.md`, `PREREG_identity_null.md`; lane Q's A1 to A4 files are named in `s26/q_*.py`
docstrings and had not landed at assembly) and the eight `s26/IDEA_*.md` tournament entries.
Every falsifier is one sentence. Standard unless stated: built chain (`rmsd_arm`, 3.2148 A,
`bench_results/cache/1fc9f2dcf489e2fb`), paired per target, `s24.stats_lib.compare`, MDE =
2.8016 x SE per comparison, decision on the fold-clustered CI, a positive needs the effect
beyond its own MDE, the fold CI excluding zero and 5 of 5 folds in sign, then a replication at
a second seed and reversed fold order inside its own CI. Nothing runs before "PHASE 0 SIGNED OFF"
except reading, training without evaluation (L16b), and probes.

## 0. What the examination fixed as the ground

- The production number reproduces exactly (E section H): 3.0483 / 3.2041 / 3.2148 / 3.2355,
  max per-target disagreement 0.0 on all four bases; T030 = 1S9Z 0.181981.
- The Hamiltonian the selector sees is the standardised rank ladder of the top-128 scores to
  0.4% of range on both traced targets (E section B); alpha = 1.0 on folds 0, 3, 4, so 78 of
  126 targets carry no tail constraint (`s25/results/q_alpha.json`).
- Four quoted numbers have no artefact (36.1/36.4 deg, +0.0004/+0.0030, the |z_moment| triple,
  355/13) and the seven-configuration suite is quoted on the point-cloud basis (E section I.2).

## 1. Proposal A: ADAPT-VQE on the deployed selector (lane Q)

Question: does growing the circuit operator by operator (qubit-ADAPT, Tang et al. 2021, pools
V / G / L2) change anything the deployed fixed RY+CNOT ansatz produces, on the deployed
objective F = CVaR_alpha(E; p) - T H(p) with E the rank ladder?

| item | hypothesis | falsifier (one sentence) | comparison arm | basis | MDE reference |
|---|---|---|---|---|---|
| A1 null test | ADAPT-grown circuits at matched parameter count (P = 3n = 21) reach a lower F and a different readout than the fixed ansatz | If the ADAPT readout's built-chain RMSD does not beat the fixed-circuit readout beyond its MDE with the fold CI excluding zero on 5/5 folds, A1 is null and the record's "readout is insensitive" stands (S25 L15, 0.24x MDE). | fixed `core.quantum.run_cvar_vqe` at (alpha, T) from `VQE_LFO`, plus uniform-over-128 and best-of-N from the untrained circuit | built chain through the same readout; the S8 selection basis carried | `s26/results/q_mde_reference.json`: p-weighted vs uniform over the same 128, effect -0.0135, SE 0.0342, MDE 0.0958 (0.14x); VQE_LFO vs argmin on the S8 instrument, MDE 0.2051 |
| A2 dynamical Lie algebra | the DLA of the deployed ansatz at depth 3 is a strict subalgebra of so(2^7), i.e. the fixed circuit cannot express the Gibbs state, which is why 0.902 nats remain | If the exact Pauli-string closure (`s26/q_dla.py`) of {C^k Y_q C^-k} reaches dim so(128) = 8128, expressivity is not the cause of the residual and the 0.902 nats is optimisation, not capacity. Lane Q's L27 reports exactly that: the full so(2^n) from depth 2 at n = 7; the prediction is falsified. | the ADAPT pools' closures; numerical nested commutators at n = 4, 5 as the cross-check | property, no RMSD | none needed (an exact count) |
| A3 target-dependent Hamiltonian | ADAPT run on a Hamiltonian that differs between targets (the raw score, not its rank ladder; or a Legacy/AMBER-mixed diagonal) selects differently from ADAPT on the rank ladder | If the per-target operator sequences and readouts on the target-dependent H are identical or within MDE of those on the rank ladder, the Hamiltonian's target dependence buys nothing and S25 L17 ("the Hamiltonian barely changes between targets") extends to the grown circuit. | rank-ladder ADAPT; fixed ansatz on each H | built chain | same reference file |
| A4 gradient variance, grown vs fixed | the log2 Var per qubit of dF/dtheta_0 for ADAPT-grown circuits at P = 3n differs from the fixed ansatz's slopes (-0.649 linear, -0.252 / -0.047 CVaR, -0.311 deployed; `s25/results/q_plateau.json`) | If the grown circuits' slopes lie inside the S25 slopes' spread at n = 4..13 with the n = 7 row reproducing `q_plateau.json`, ADAPT does not change the variance regime; no slope may be called a barren plateau without this width sweep (Rule 10). | S25 fixed-ansatz sweep, reproduced first | property | none (a slope with its fit CI) |

## 2. Proposal B: a learned folding model as the prior (lane P)

| item | hypothesis | falsifier | comparison arm | basis | MDE |
|---|---|---|---|---|---|
| B1 feasibility | ESMFold v1 can run on this box | Measured and failed on three grounds (`s26/results/b1_feasibility.json`, L13): no checkpoint, `omegaconf`/`openfold` missing under Python 3.13 CPU torch, 8.76 GB resident against 4.4 GB headroom; B1 stops and the endpoint half does not run. | none | none | none |
| B2 feasible-scale replacement | the contact head alone (`conly`, 55-d) or the 8M-parameter model (`esm8m`) beats the shipped pca32 prior, and the size axis 8M to 650M is monotone | If no rung beats the shipped posterior on the built chain beyond its MDE with the fold CI excluding zero on 5/5 folds, and `esm8m` is not worse than `pca32` beyond MDE, the language-model size axis is closed from below and "a bigger model" has no measured support. | shipped posterior through `s26/p_ladder.py`; `noesm` as the zero-ESM control | built chain primary; point cloud and selection carried | 0.09 to 0.20 A on the built chain (PREREG_C2 section 4); 0.28 on selection |
| B3 repo-native stand-in | the pipeline's gain over sequence-only (arm 3.2148 vs torsion predictor 3.7705 vs constant helix 4.0648, all built chains) is concentrated on a native-free characterisable subset | If a nested leave-fold-out ridge classifier of sign(d_i) does not exceed the 95th percentile of a 500-draw label-permutation null and the ridge regression of d_i has held-out R^2 <= 0, the subset is not characterisable at n = 126 and B3 is a power statement. | label permutation; FAIL18 reported as an ORACLE stratum only | built chain on all three sides | MDE on d = 0.369 A |
| replacement rule | applies to B1 only | B is replaced by B2 + B3 per `s26/PROPOSAL_B_REPLACEMENT.md`; the arm-choice oracle over {pipeline, predictor, helix} (2.9632, -0.2515) is 92% an order statistic (`ST.best_of_k_within`, k_eff 2.34, L14) and is not a routing result. | | | |

## 3. Proposal C: a better prior and its consumers (lanes P and PH)

| item | hypothesis | falsifier | comparison arm | basis | MDE |
|---|---|---|---|---|---|
| C1 closure reproduction | the two closures (flat learning curve 3.043 to 3.026 vs leaked label 2.534; perfect in-band ranker inside top-25 = 2.609) reproduce from persisted artefacts | If any quoted number is off by more than 0.001 A from the named artefact, the closure is downgraded to "cited, not reproduced" and ledgered. Result on disk: both reproduce exactly (PREREG_C1). | none (a reproduction) | as in the originals (raw cloud; selection) | none |
| C2 prior-input ladder | a distogram trained on inputs the box can compute (pca32f, pca128, raw, esm8m, wide, pairnet, mix, conly) is a better prior and the gain survives the pipeline | If no rung beats the shipped posterior on the built chain beyond its own MDE with the fold CI excluding zero and 5/5 folds, H_C2 is refuted at this ladder and every rung's power is stated. | shipped posterior through the identical path; `pca32` retrained must reproduce the pinned model (gate passed: argmin 25/25, top-75 overlap 1.000, `s26/results/p_ladder_gate_pca32_fold0_s0_probe.json`); `mix` at lam = 0 bit-exact | built chain primary; cloud and selection carried; gam_eff with cos beside it, never -2.1496 x gam_eff alone | 0.09 to 0.20 A built chain; the record predicts every ESM rung inside 0.05 A |
| C3 AMBER refinement with a matched-random control (lane PH; L5, L16b) | the production relaxation's +0.0207 A (3.2148 to 3.2355) is an accuracy step, not the cost of any move of that size | If AMBER does not beat BOTH the S16-form isotropic random displacement of matched CA magnitude and the toward-a-pool-member move of matched magnitude beyond each comparison's MDE with the fold CI excluding zero on 5/5 folds, "refine with physics" is a validity step only. | `rand` (16 draws), `member` (16 draws), do-nothing; stage 2 repeats on the best C2 rung | built chain, CA-RMSD after superposition | SE 0.0034 for a step that moves every target; the registered prediction is AMBER worse than random by +0.01 to +0.02 |
| C4 native-free routers for s* and m | principal-axis spread, retrieval-score entropy, posterior bin entropy and contact-map statistics predict m* or s* held out | If the routed endpoint does not beat fixed m = 75 (or s = 1) beyond its MDE with the fold CI excluding zero on 5/5 folds, or is harmful (fold CI above zero), C4 stops; a positive s* router is a leak to hunt, since s* is a function of the invisible common mode. | fixed incumbent; 200-draw label permutation; the S22/S23 feature set re-run as a harness control (must give ~0) | point cloud (where m*, s* live) and built chain | m router MDE 0.055 (S23 L7 SE 0.0196); s router MDE 0.157 |
| C5 common-mode prediction | the common-mode error (68% of the pool's squared error, `s23/results/errdecomp.json` f_common mean 0.6758) has a predictable component that a training-fold model can subtract | If GLOBAL or RIDGE minus incumbent is not beyond its MDE with the fold CI excluding zero on 5/5 folds AND does not beat RANDOM-MATCHED by the same standard, C5 is refuted at one agent-day and the S16/S19/S23 closure gains a fourth instrument. | ORACLE ceiling (subtract the true ebar), RANDOM-MATCHED, SCALE-ONLY (s*) | built chain primary; point cloud carried | 0.05 to 0.09 A |

## 4. The tournament (mandatory directions and own ideas)

| entry | falsifier in one sentence | control | basis |
|---|---|---|---|
| AMBER as a steric reject filter (`PREREG_amber_reject.md`, PH) | If neither `R_1e4` (reject e > 1e4 kcal/mol, refill to 75) nor `S_1e4` (no refill) beats the anchor beyond its MDE with the fold CI excluding zero on 5/5 folds AND beats its count-matched random reject, the direction is closed in the two forms not previously measured. | `RANDR`, `RANDS` (16 draws), `PERM` (energy permuted within the pool) | point cloud and built chain |
| the cis-peptide gap (`PREREG_cis.md`, PH) | If no model-1 native carries a cis bond (the registered prediction, because `core/data.py` drops steps outside 3.5 to 4.1 A) and the CA floor on any cis target is below the MDE of the chain-cost comparison, the two-bond-length projection is a design for the world supply only. | omega census on every deposited model; the non-planarity tail | ORACLE diagnostic |
| identity null (`PREREG_identity_null.md`, I) | Done: the shorter-normalised criterion at 0.6 is at the null (real 0.5% vs shuffled 0.3% pass rate, ratio 0.62 against a falsifier of < 0.2), so 71/126 is chance and the folds stay pinned (L15). | composition-preserving shuffles | none |
| physics chooses the projection branch (`IDEA_branch_select.md`, PH) | If picking among the multi-start branches by converged AMBER energy does not lower the built-chain RMSD beyond MDE against the objective's own pick and a random branch pick, closed. | random branch; objective pick | built chain |
| strain as difficulty (`IDEA_strain_difficulty.md`, PH) | If `amber_e0`, `amber_e1`, `amber_moved`, `amber_strain_after` do not predict `rmsd_arm` held out beyond a label-permutation null, it is not a calibration output. | permutation null | built chain (a correlation, no emission) |
| rotamer relief (`IDEA_rotamer_relief.md`, PH) | If greedy chi1 relief does not remove most of the above-1e4 mass, the singularity is the pool's, not the builder's; the relieved single point is then re-measured as a ranker before any use. | raw single point | property |
| better prior inputs (`IDEA_better_prior_inputs.md`, P) | The two open inputs (raw attention maps per pair; the pool histogram as a training-time input) are rungs of C2 and fall under its falsifier. | C2's | built chain |
| coherence-penalised training (`IDEA_coherence_penalised_training.md`, P) | If the short-separation-weighted loss (S19 A5b, unblocked by the compact ESM tables) and the shell-coherence penalty do not beat the uniform-loss retrain beyond MDE, coherence is not trainable away. | uniform-loss retrain; the anti-utility `sw_lin` arm | built chain |
| window ensembling (`IDEA_window_ensembling.md`, P) | If the three-key (BLOSUM45/62/80) or three-K ensemble of top-75 clouds, projected once, does not beat the shipped chain beyond MDE with the fold CI excluding zero, fixed-K ensembling adds nothing that widening K did not already price (S17 L12). | shipped single cloud; a same-key three-seed ensemble | built chain |

## 5. What the Adversary should attack first (from the examination)

1. The point-cloud basis of the seven-configuration suite and of the +0.330 / +0.455 physics
   verdicts in presenter documents (E section I.2).
2. The alpha = 1.0 folds: any Proposal A claim about "the CVaR tail" applies to 48 targets, not
   126, unless the arm re-runs the table.
3. The unconverged 9KAR relaxation and the strain-flagged 2BP4 relaxation inside 3.2355, and
   the unversioned production cache (E section I.4, I.5).
4. The four document-only numbers, if any of them reaches a slide.
