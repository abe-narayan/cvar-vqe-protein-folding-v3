# PREREG S29 / LANE X -- THE CONFIGURATION-SPACE CVaR-VQE (the divergent direction, H2)

Written 2026-09-19 23:50 Pacific, before any number of the probe exists. Appended, never edited,
after the first result (S28 contract rule 7). Brief: `s29/briefs/S29X.md`. Contract: `s29/S29_CONTRACT.md`.

## 0. The question, in one sentence

Is there a state space in which the quantum state is a distribution over STRUCTURAL HYPOTHESES
that the pool average cannot represent, with an interacting Hamiltonian, such that the CVaR
tail ENSEMBLE (not the argmin) emits a structure and removing the quantum stage degrades it?

## 1. What the record closed, and the new angle against each ledger line (contract rule 10)

| closed | where | what was closed | why it may not apply here (the new angle) |
|---|---|---|---|
| Torsion-bin space, k = 4 (~24 live qubits) | S13 dossier §1-4 (`s13/SPRINT13_DOSSIER.md`) | ORACLE descent 1.594 A; no native-free objective (Legacy, AMBER, the 1-local prior) finds it; Legacy's certified optimum is +0.139 worse than random; the prior's argmin is a constant helix | The space here is not Ramachandran bins but RECOMBINATION OF RETRIEVED FRAGMENTS (3-mers of real pool members, internal torsion correlations kept); the objective is the posterior, not a physics energy; the readout is a tail ensemble, not the argmin |
| Distance geometry on the predicted posterior | S15 dossier §1.1-1.3, §5 | From true distances 0.611 A; from the predicted posterior 3.644 A (+0.440 worse than production): the distogram's error is a coherent WRONG structure that a faithful optimiser builds | S15's optimiser was unconstrained in torsion space and could realise the wrong structure. A chimera of real fragments cannot realise an arbitrary distance map; the fragment prior constrains what the posterior can pull the structure toward. Whether that constraint helps or hurts is the measurement |
| "Do not collapse the posterior" | S14/S15 (memory `torsion-restraints-reach-the-target`, ORACLE -2.253 A for top-8 support vs argmax) | Collapsing a per-residue posterior to its argmax costs 2.3 A at the ORACLE | The pair term here consumes the FULL 17-bin posterior of every pair (the log score), not its median or mean; the 1-body term is the full Ramachandran table |
| The continuous basin latent (1 qubit / residue, von Mises draws) | S20 L1/L3 (`s20/LEDGER.md`), S21 L11/L14/L17/L20 | The VQE trains (-0.210 vs untrained) but no sampler beats the zero-evaluation pool; the exhaustive latent argmin TIES the pool argmin (+0.025 / +0.085, n = 75 / 126); the latent as a SOURCE through the shipped operator is +0.390 worse than the pool; its ORACLE is 0.56 worse than the pool's | That latent was a parametric prior FITTED to the pool's per-residue marginals (S20 L3: "the least likely thing to decorrelate"). The chimera space is built from the pool members themselves, keeps each fragment's multi-residue correlations, CONTAINS the members it recombines, and its ORACLE ceiling has never been measured |
| The tail average of the pool IS production | S21 L3 | tail_avg(alpha = 0.15) over the pool = the shipped top-75 average to 1e-6; a CVaR-VQE that selects an alpha-tail of the pool and averages it re-derives the incumbent | The tail here is over CHIMERAS: the pool average is not in the space and the average of the tail is not the pool's tail average |
| Search saturates, discrimination binds; in-pool selection is where the Angstroms are | S21 L14/L17/L18/L19/L21/L23 | Ordering skill of the objective is at M = 1 and exhausted by M ~ 8; bimodal containment; consensus family exhausted; 1.338 A in-pool ORACLE gap uncaptured | Not re-opened: I ACCEPT that no in-pool selector captures the gap. The chimera space is a different set; whether the same objective's exact argmin over it ties the pool argmin is the pre-registered classical control (S21's control, transported) |
| The first non-diagonal Hamiltonian on the pool | S28-L21/L22/L41/L43 (`s27/LEDGER.md`) | H = diag(E) - J A with A the Gaussian similarity graph (near rank one); the deployed readout R1 moves <= 0.013 A at any J; the weighted readouts R2/R3 are +0.24 to +0.33 above production at every J; the J = 3 ground state is not representable (overlap 0.80-0.88) | The off-diagonal term here is a TRANSVERSE FIELD on configuration bits: transitions between hypotheses that differ in one fragment index, the QAOA meaning, on a space where a neighbour is structurally meaningful (one segment swapped). Its spectrum is not a similarity graph's and rank collapse of a kernel does not arise (the mixer is the hypercube adjacency, spectrum {q - 2 w : w = 0..q}) |
| The set-equality theorem | S24 `d_harness` docstring, S25 (5.4 of the state brief), S28-L21 | For a diagonal H the realised CVaR tail is a prefix of the energy order; the deployed uniform-tail readout can only choose m | ACCEPTED AS BINDING and stated up front: my readout R1 (uniform over the tail set) is the deployed operator and is prefix-reducible BY THEOREM. The only readouts in which the quantum state can matter beyond m are the WEIGHTED ones (R2, R3); the quantum stage's whole possible contribution beyond m is R2 - R1 on the same set, and it is measured as such |

Charter findings (rule 13): this hypothesis ATTACKS 9 (Hamiltonian variation within a fixed encoding is
a low-dimensional search: the encoding changes), 10 (the pool contains better structures: recombination
may contain better ones still) and 11 (aggregation: an average over chimeras that share fragments is not a
consensus over independently retrieved members). It ACCEPTS 1-7 as binding. It ACCEPTS 8 (recognition)
as the RISK: the objective is a function of the same pair marginals, so if nothing can recognise a good
chimera, the tail ensemble cannot either; the ORACLE percentile of the best chimera under the objective
is measured to say so. Rule 14 (information test): the one method imported is fragment recombination
(Rosetta-style 3-mer assembly restricted to pool fragments). What information does it contain that the
project does not? Structures outside the pool whose every local 3-mer is a real retrieved fragment
with its internal correlations intact -- the pool's near-native LOCAL pieces recombined. Its value is
exactly its ORACLE ceiling relative to the pool's, which is measured first (section 6, D1).

## 2. Encoding: (ii), fragment recombination, and why

(i) per-residue (phi, psi) bins is S13's space (B = 4) and, at B = 2, S19-S21's basin latent: both have
exhaustive-enumeration closures and neither preserves multi-residue correlations. (ii) has never been
run. Chosen: (ii).

- Members: F = 8 pool members = the DIS top-8 of the shipped K = 500 pool (`s24.d_harness.score_shipped`,
  the production functional), exact ties broken by the stable random key of `s27.run_pool.topm`
  (contract rule 12). The top-8 because S21 L18 located the objective's ordering skill at M = 1..8.
- Segments: S = ceil(n / 3) contiguous segments, `np.array_split(range(n), S)` (sizes 3,3,...,2,2);
  n = 9..16 gives S = 3..6. Each segment chooses ONE of the 8 members; the configuration takes that
  member's (phi_i, psi_i) for every residue i in the segment. 3 qubits per segment, binary index in
  DIS-rank order (index 0 = the DIS argmin). Register q = 3 S in {9, 12, 15, 18}; the space has 8^S
  <= 262,144 configurations, exactly enumerable, so every classical control is EXACT.
- Build: `core.project.build_ca_exact(phi, psi)` (the production ideal-geometry builder). A
  configuration whose S choices are all member k is member k's rebuilt trace, bit-exact (test 1).
  No dead qubits: phi[0], phi[n-1], psi[n-1] are inert for the CA trace (S13 correction 10) but no
  segment consists of inert torsions only, and the 1-body term excludes them.
- Basis: every chimera arm is on the REBUILT manifold. Production is quoted on its own basis
  (window point cloud 3.0483 / built chain 3.2126); S21 L20 priced the rebuild on the average at
  +0.016 [-0.004, +0.037], negligible and stated.

## 3. The Hamiltonian

H = H_diag + H_mix on q qubits.

H_diag(x) = E_pair(x) + E_rama(x), in nats:
- E_pair(x) = sum over pairs (i, j), j - i >= 2, of -log(p_ij[bin(d_ij(x))] + eps), p_ij the shipped
  leave-fold-out distogram's 17-bin posterior (`s12.instrument.distogram`, `core.predict.BIN_EDGES`),
  d_ij the realised CA-CA distance of the built configuration, eps = 1e-4 (bounds one pair at 9.2
  nats: the S13 delta-spike lesson). The log score, unweighted: a proper local scoring rule of the
  posterior. Evaluated EXACTLY on the built configuration (the S13 locality theorem makes d_ij a
  function of the j - i - 1 residues between i and j, hence of up to j - i - 1 segment variables; the
  term is (j - i)-local, not 2-local, and no Ising truncation is made).
- E_rama(x) = -sum over live residues i of log P_fold(phi_i, psi_i | aa_i), the per-fold table
  `s8/generate_rama.npz` with the S27 RAMA channel's binning and pseudocount (`s27.ham_lib.h_rama`,
  PSEUDO = 1, 36 x 36 over [-pi, pi)); weight 1 (nats add). Live residues exclude phi[0], phi[n-1],
  psi[n-1] where the whole residue is inert (residue n-1); residue 0's psi is live so residue 0 is
  kept with its phi included as in the channel (a constant per member choice, harmless).
- The Gibbs distribution exp(-H_diag) at T = 1 is the pseudo-posterior over chimeras (pairs treated
  as independent); it is the distribution the quantum state is asked to represent.

H_mix = -Gamma sum_k X_k (the transverse field on every qubit): the transition operator between
hypotheses differing in one bit of one segment's fragment index. Meaning: the state can hold amplitude
coherently across fragment alternatives; its expectation M(theta) = <psi|sum_k X_k|psi> = sum_k
2 sum_x psi(x) psi(x xor e_k) measures amplitude coherence between Hamming-neighbour chimeras.
Gamma_gap = median over the H_diag top-75 configurations and over the q bits of |H_diag(x) -
H_diag(x xor e_k)| (nats; native-free; the scale at which first-order neighbour amplitude
Gamma / dE is order one). Rungs: Gamma in {0, Gamma_gap}. Nothing about Gamma is chosen on RMSD.

The objective (the spine, S28 lane B's form with the graph replaced by the mixer):
F(theta) = CVaR_alpha(H_diag; p_theta) - T S(p_theta) - Gamma M(theta), alpha = 0.15 (production's
tail fraction), T = 1 nat (the posterior's own temperature; S25 L2 says the posterior is ~2x
over-confident and calibrating it makes RMSD worse, so T is NOT tuned; a T = 2 rung exists only on
the classical Gibbs arm as a diagnostic). Ansatz: `core.quantum.StatevectorCircuit(q, layers = 3)`
(RY on every wire, CNOT chain, ring), exact statevector, exact parameter-shift gradient for all three
terms (M is a sum of Pauli expectations, so the shift rule is exact), Adam lr 0.15, 80 iterations,
theta_0 ~ N(0, 0.6), seeds 0 and 1: the production driver's settings (`core.quantum.run_cvar_vqe`).

## 4. Readouts

The CVaR tail = the alpha-tail of p_theta in H_diag order (`core.quantum.cvar_exact`'s tail set:
the prefix whose cumulative mass reaches alpha; realised size m recorded).
- R1 tail-uniform: coordinate average of the tail set, uniform weights, medoid frame
  (`s12.instrument.coordinate_average`): the DEPLOYED operator; prefix-reducible by theorem.
- R2 tail-weighted: the same set, the same medoid frame, weights p_theta(x) renormalised on the
  set: the PRIMARY quantum readout. R2 - R1 on the same set is the quantum stage's whole contribution
  beyond m.
- R3 full-state: p_theta-weighted average over the top-512 configurations by probability (or the
  0.999-mass set if smaller), medoid frame of that set: the state's mean structure, diagnostic.
Then `s12.instrument.project` (the production projection) -> the built chain, the reporting basis.

## 5. Arms and controls (all pre-registered, all on both bases, all 12 targets)

Native-free arms:
- VQE Gamma_gap seed 0 (the arm), VQE Gamma_gap seed 1 (second seed), VQE Gamma = 0 seed 0 (the
  diagonal twin), VQE product-state restriction (CNOT permutation removed, Gamma_gap, seed 0),
  UNTRAINED circuit (theta_0 of seed 0, no training).
- GS: the EXACT ground state of diag(H_diag) - Gamma_gap sum X (Lanczos on 2^q), p = g^2 (the
  diagonalised equivalent); overlap |<g|psi_theta>|^2 recorded (reachability).
- GIBBS(T = 1): p ∝ exp(-H_diag); GIBBS(matched): T solved so that S(p) = S(p_theta) of the seed-0
  Gamma_gap arm (the classical thermal ensemble at matched entropy); GIBBS(T = 2) diagnostic.
- EXACT classical ladder over H_diag: argmin (m = 1), top-8, top-75 (production rung), top-m at the
  VQE's realised m (size-matched, the S24 convention), uniform averages.
- SA: single-flip Metropolis over the q bits, budget 2^q evaluations (the exhaustive cost, S21 L11),
  geometric cooling 20 -> 1 nat, distinct visited configurations with visit counts; tail = top-m
  distinct visited by H_diag at the VQE's realised m; R1 uniform, R2 weights = visit counts.
- PERM: the posterior's rows permuted among pairs of the same sequence separation (seeded,
  native-free); VQE Gamma_gap seed 0 and the exact argmin / top-75 on the permuted H_diag.
- PRODUCTION: `s27/results/chain_rows.jsonl :: DIS seed 0` per target (rmsd_cloud, rmsd_chain).

ORACLE diagnostics (labelled ORACLE in every sentence; never choose anything):
- the chimera space's ORACLE best (min RMSD over all 8^S), point cloud and chain; the DIS top-8's
  best member; the percentile of the ORACLE-best chimera under H_diag (S21 L18's statistic); the
  native's H_diag percentile among chimeras (the meter's (c) analogue); the ORACLE best inside the
  realised tail; the mutual information between segment choices under p_theta vs GIBBS(matched)
  (mechanism beside outcome, rule 18).

NaN-poison (rule 7): every native-free arm is recomputed with nat_ca and oracle_rr replaced by NaN
and must be bit-identical (test 3).

## 6. Falsifiers, registered before the run

D1 (the space). Chimera ORACLE best vs the DIS top-8's best member and vs the pool's ORACLE best
(1.7108 over K = 500, `s12.instrument.selfcheck`). Prior: recombination is worth 0.3 to 0.8 A at the
ORACLE over the top-8's best member on the point cloud; if it is worth < 0.1 A the space is not richer
than its parents and the lane's remaining arms are formalities.

D2 (the classical control, S21's transported). Exact argmin of H_diag over the chimeras vs
production's argmin basis is not available on the chain per target; the registered contrast is the
exact top-75 uniform (production rung) vs production, chain. Prior: ties or worse (S21 L14/L17/L20).

P1 (the endpoint of the probe). R2(VQE Gamma_gap seed 0) - production, built chain, ST.compare with
the 12 targets' pinned folds. GO to the 126 iff effect <= -0.7 x MDE (better). Registered prior:
WORSE than production on the chain by 0.1 to 0.4 A (the S13 record: torsion-space arms lose on the
chain; S28-L21: weighted readouts sit +0.24 to +0.33 above production on the pool).

P2 (SA ties the VQE). R2(VQE) - R2(SA), chain. Prior: |effect| < 0.7 x MDE.

P3 (the quantum stage beyond m). R2(VQE) - R1(VQE, same set), chain. Prior: within 0.7 x MDE; if
anything R2 is worse (concentration, S23's closure of probability weighting on the pool).

P4 (classical equivalent at matched entropy). R2(VQE) - R2(GIBBS matched), chain. Prior: within
0.7 x MDE. The claim "the quantum stage contributes" requires P3 AND P4 to clear 0.7 x MDE on the chain
with the fold CI excluding zero, on both seeds, and the PERM arm not to reproduce it.

P5 (the mixer). R2(VQE Gamma_gap) - R2(VQE Gamma = 0), chain. Prior: within 0.7 x MDE.

The registered prior for the whole probe, verbatim from the brief: ties the pool on the point cloud,
worse on the built chain, and the SA control ties the VQE; the value is settling whether the
tail-ensemble readout changes S21's closure. If P1 fires (an arm better than production at 0.7 x MDE
on 12 targets), the 126 runs on that arm and its controls, after lane D's meter has seen H_diag.

Multiplicity: the probe's endpoint contrasts against production are counted per arm and readout in
the ledger entry (native-free arms x {R1, R2, R3} on the chain); a GO at 1x MDE among K contrasts is
priced against the max-over-K null (rule 17). ORACLE rows are diagnostics and are not counted.

## 7. The nine questions of charter section 11

1. Basis-state meaning: a chimera -- for each of S segments, which of the 8 DIS-top members supplies
   that segment's backbone torsions; 8^S structures, every one a real-fragment recombination built
   with ideal geometry; the 8 parents are the 8 "diagonal" configurations.
2. Hamiltonian meaning: H_diag is the negative log pseudo-posterior of the chimera under the shipped
   distogram (pair log score) plus the Ramachandran log prior; exp(-H_diag) is the pseudo-posterior
   over chimeras. H_mix is the transition operator between chimeras differing in one fragment-index
   bit; Gamma sets the tunnelling scale in nats.
3. Why it should correlate with structure: the distogram orders the bulk (S13: the exact argmin beats
   the space mean by 1.1-1.4 A; S21 L14) and the objective's top-1 is 1.18 A better than a random
   configuration (S21 L18); the prior places the optimum in a Ramachandran-legal region (S13). What the
   record says it CANNOT do is order the top ~8 (S21 L18/L19), and the S15 mechanism says the posterior
   describes a coherent wrong structure. The fragment constraint is the one thing new: whether it stops
   the posterior from building its wrong structure is D1 + P1.
4. What CVaR optimises, over what: CVaR_0.15 of H_diag under p_theta, the computational-basis
   distribution over chimeras, minus T x entropy (an ensemble, not a point) minus Gamma x coherence.
   Its minimiser concentrates the alpha-tail on low-H_diag chimeras (S20 T1: a face on the tail) while
   the entropy and mixer terms spread mass and amplitude to neighbours.
5. What the ansatz can and cannot represent: q x 3 = 27..54 real parameters over 2^q = 512..262,144
   amplitudes. It cannot represent a generic state concentrated on ~10-100 specific chimeras (S28-L41:
   overlap 0.80-0.88 with the J = 3 ground state at 9 qubits); the GS arm and the overlap measure how
   far it is from the exact optimum of the diagonal-plus-mixer problem. Real amplitudes: the mixer's
   ground state is stoquastic (non-negative amplitudes), so real amplitudes are not a restriction for
   it. The product-state restriction removes all correlation between segments.
6. Whether the optimiser reaches the relevant states: measured by the three-way split (contract
   addendum 14c): the exact GS energy vs the reached state's energy vs the emitted structure; the
   overlap; F at convergence vs F of GS and of GIBBS.
7. What the quantum output could contain that the classical control cannot: a coherent superposition
   over mutually compatible chimeras with amplitude on Hamming neighbours of the tail set, i.e.
   correlations between segment choices that a thermal ensemble at matched entropy does not carry
   (measured as the segment-pair mutual information under p_theta vs GIBBS(matched)). The record's
   expectation: a stoquastic ground state's |amplitude|^2 is classically samplable in principle and
   at these register sizes tabulable exactly (S21 L4), so any effect is a property of the OBJECTIVE,
   never a quantum resource (contract rule 9).
8. Whether a classical control reproduces it: GIBBS(matched entropy), SA, the exact top-m, the
   product state, the PERM posterior -- P2, P3, P4, P5.
9. Whether the built chain moves: P1, on the built chain, ST.fmt, both seeds, before any verdict.

## 8. What is fixed by rule (no RMSD anywhere in these choices)

F = 8 (brief), L = 3 (Rosetta 3-mers; the smallest fragment with an internal pair), members = DIS
top-8 with the stable tie key, eps = 1e-4, alpha = 0.15, T = 1 nat, Gamma_gap by the median-gap rule,
layers 3, iters 80, lr 0.15, seeds {0, 1}, SA budget 2^q with cooling 20 -> 1 nat, R3 cap 512, PERM
within-separation permutation seeded by `s15.seed.stable_rng(pdb, "s29X_perm", salt = "s29X")`.
The 12 targets: `s25.phys_lib.targets()[::11][:12]` (S27's trainability set): 1A13 1I6Y 1M02 2BFI
2LWS 2MP9 2P5H 5Z5W 6MBM 7JGX 8HVS 9KAR (n = 14, 11, 12, 12, 12, 12, 9, 12, 16, 13, 10, 15; folds
4, 2, 0, 2, 4, 0, 4, 3, 1, 0, 4, 3; 2BFI and 9KAR in FAIL18).

## 9. Files

`s29/s29_X_config.py` (the space, the Hamiltonian, the arms; the cost callable
`s29.s29_X_config:cost_nll` in lane D's meter contract), `tests/test_s29_X.py`,
`s29/results/s29_X_probe_<pdb>.json` (per target, resumable), `s29/results/s29_X_probe.json`
(aggregate with every ST.fmt block), `s29/s29_X_FINDINGS.md`, ledger `S29-L<n> (date time, X)`.
Jobs: `python s26/jobrun.py --agent S29X --tag CPU --name s29X_probe --est-ram 2 -- python
s29/s29_X_config.py --run --limit 1` first (peak RSS quoted), then the 12.

## 10. What is not run, and why

- The 126 (rule 16: a probe first; GO only on P1).
- AMBER anywhere (no physics energy ranks nativeness: S25 L16; not part of the hypothesis).
- Encoding (i) (closed twice by enumeration: S13, S21).
- Any tuning of alpha, T, Gamma, F, L, the member rule or the readout on RMSD.

---

## ADDENDUM 1 (2026-09-20 00:10, before any probe number exists)

**A. The register cap, and the measured reason for it.** Section 2 said S = ceil(n / 3). Measured
cost of one exact parameter-shift gradient of the deployed objective (`core.quantum.free_energy`,
layers 3, this box, one thread): q = 12 -> 0.10 s, q = 15 -> 1.61 s, q = 18 -> 24.40 s, i.e. 0.1 /
2.1 / 32.5 minutes for the 80-iteration run. At S = ceil(n/3) the n = 16 target needs q = 18 and
one arm costs half an hour. The rule becomes **S = min(ceil(n / 3), 5)**, so q = 3S <= 15 and the
probe's four longest targets cost ~2.5 min per arm instead of ~33. `np.array_split` then makes the
n = 16 target's segments 4,3,3,3,3 rather than 3,3,3,3,3,1. This is a COMPUTE decision taken before
any RMSD existed; it shrinks the space on one target of twelve and is stated in the ledger entry.

**B. Lane D's four holes (S29-L4) are adopted in full, before the run.**
- (a) D1 acquires its matched order-statistic null: a SCRAMBLED-CHIMERA space of identical
  cardinality, identical parents and identical marginal fragment content, built by permuting each
  member's segment blocks across segment POSITIONS (`Space(..., scramble=True)`, seeded
  `stable_rng(pdb, "scramble", salt="s29X")`). **D1's primary becomes chimera ORACLE best minus
  SCRAMBLED ORACLE best**; the contrast against the eight parents is reported as UNMATCHED and
  labelled so in the same table.
- (b) D1's branch is registered: if recombination is worth less than 0.1 A against the SCRAMBLED
  null, the lane reports D1 and the probe's already-computed arms, and starts nothing further on
  this space without saying why in the same entry.
- (c) P3 acquires the participation ratio of the R2 weights (1/sum w^2) and a PR-MATCHED
  RANDOM-WEIGHT control: 8 Dirichlet draws on the same tail set at the same realised PR, each
  through the identical readout and projection (`P3c`). The emitted Rg and mean virtual bond are
  printed for every arm (contract addendum 20(c)).
- (d) Any GO prints its power and its Type-M factor and says "GO, not a result" in the same
  sentence; the fold-clustered CI at n = 12 is descriptive only; the 126-target run re-registers
  its own falsifier rather than inheriting this one.
- (e) The PERM arm reports its seed and its H_diag spectrum (sd, range, top-75 range) beside the
  real one.
- (f) R3 reports the probability mass its top-512 captures; without it R3 is named a top-512
  readout, not a full-state readout.

**C. Contract addendum 21 (the marginal class is bounded, lane T's theorem 2).** H_diag is inside
that class: it is a functional of the shipped posterior's per-pair marginals plus a
sequence-conditioned 1-body prior. What this lane expects to buy is therefore NOT local
informativeness and NOT a gradient (the meter measures the gradient as undefined on 114/126
targets anyway, S29-L6): it is (i) ORDERING over a different candidate SET -- chimeras, which the
pool does not contain -- and (ii) NON-CONTRACTION, an ensemble over recombinations rather than an
average over 75 pool members. Both are measured here. No gradient claim is made, and no cosine
gain is claimed, so addendum 20's shrink audit is answered by reporting (b) the native percentile
(0.378, S29-L6, essentially the shipped cost's 0.368) and (c) the emitted Rg and bond, both of
which are in the per-arm rows.

---

## ADDENDUM 2 (2026-09-20 00:35, before any probe number exists; the coordinator's upgrade and
## lane T's S29-L15 falsifier ladder)

The coordinator raised this lane from "the divergent direction" to the sprint's principal quantum
candidate on lanes L (S29-L13) and T (S29-L15): a formulation is non-classical in the relevant
sense only if the Hamiltonian's terms do not commute AND the prepared object is not an
eigenvector, and T's variance law leaves exactly one operator class with extensive stable rank --
a sum of LOCAL Pauli terms, which is meaningful only where basis states have local structure,
i.e. here. The formulation T names as the smallest qualifying one is the one already registered
in section 3 above (H_diag + Gamma sum_k X_k; a free-energy target, never an eigenvector; the
CVaR tail's coordinate average as the readout). Three things are ADDED, before any number:

**G. The coupling grid.** Gamma in {0, 0.5, 1, 2} x Gamma_gap, seed 0, trained identically. Full
readouts (R1, R2, R3, both bases) for Gamma in {0, 1} x Gamma_gap (the registered arms); R1 and R2
only for 0.5 and 2, which exist to place the gate rather than to emit an endpoint.

**GATE 1 (T's cheapest falsifier, read BEFORE any endpoint number).** TV(p_theta(Gamma),
p_theta(0)) on the sampled distribution must exceed **0.45** -- S25 L15 measured that this readout
cannot resolve a 45%-of-mass distributional difference. Registered branch: **if no grid point's
mean TV clears 0.45, the cell is empty for this instrument, no build follows, and THAT is the
result** -- reported as such, with the endpoint arms that were computed anyway reported beside it
and not used to argue past the gate. If it clears, gate 2 is the correctly-named classical
counterpart (a classical thermal sampler / SA over the same configuration space at matched
evaluations -- the registered SA arm, which becomes the MAIN comparison, not the eigensolver),
and only then gate 3, the endpoint contrast on the built chain.

**M6, lane T's fixed-profile control, in its strongest form.** The EXACT minimiser of the same
objective at Gamma = 0, computed classically with no circuit and no optimiser: by
Rockafellar-Uryasev plus Sion's minimax theorem, p*(x) proportional to exp((t* - E_x)_+/(alpha T))
with t* the maximiser of t - T log sum_x exp((t - E_x)_+/(alpha T)) (`cvar_optimal_law`; verified
to beat 200 random Dirichlet laws on F and to reproduce the Boltzmann law exactly at alpha = 1,
TV 0.0, which is lane T's own assertion). Contrasts M6 (VQE at Gamma = 0 vs p*) and M6g (VQE at
Gamma_gap vs p*) are registered. Prior: the Gamma = 0 arm ties p* or is worse (it is an
under-trained approximation of it, S29-L15's measured 8.836 vs 8.819 bits).

**M5, T's flat-fraction clause.** Every arm reports the fraction of simplex directions on which
the CVaR term is exactly constant, (M - m - 1)/(M - 1). Stated in advance: at Gamma = 0 this
formulation inherits the deployed flatness exactly -- it is a property of CVaR, not of H -- so
the Gamma = 0 arm is NOT expected to break T's M6 control. The mixer term is not a function of p
at all, which is the only structural reason Gamma > 0 could differ, and gate 1 is the test of
whether that difference reaches the readout.

**One scope note on T's reduction, stated so it is not over-read either way.** T's Q1 reduction
("the quantum stage is a target-independent rank-weight profile whose only endpoint channel is
m") rests on E = zrank(scores) being the same standardised rank ladder on every target to 1.18%
of range (S25 L17). H_diag here is in physical nats with target-specific gaps and is NOT
rank-standardised, so the profile is target-dependent by construction. The entry reports the
realised spread of H_diag (sd, range, top-75 range) so the reader can judge how much of T's
reduction survives here rather than taking my word for it. This is a reason the reduction may not
transfer; it is not a claim that anything is gained.

**Ordering.** The run is cheapest-first by register size (q = 9, then 12, then 15), so gate 1's
quantity arrives on the small targets within minutes, as T's ladder requires.
