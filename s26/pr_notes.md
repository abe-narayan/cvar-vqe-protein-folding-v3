# Speaker notes for `vqe_research_overview.pptx` (PR lane, Sprint 26)

`s26/pr_build_deck.py` reads this file. Each `## Slide N` section has a `### spoken` block (the
text the presenter says; the build counts its words) and an optional `### also` block (lines
appended to the notes, not spoken). Every `{TOKEN}` or `{TOKEN:format}` is replaced at build time
by the value read from the artefact named in `s26/pr_values.py`; the build appends a SOURCES
list to every slide's notes with the token, the value and the artefact path, so no number in the
notes is typed by hand. Bases are named where an RMSD is spoken: built chain (the production
result), point cloud (an intermediate), single window (one pool member), selection (the S25
quantum instrument).

Style: short sentences, first person, no stock words, no dashes as punctuation. Never a quantum
advantage. Never a barren plateau from a small gradient.

## Slide 1 -- Title

### spoken
Hello. I built a pipeline that predicts the backbone shape of short peptides, nine to sixteen
residues, from their sequence. On {N_TARGETS} held-out peptides it reaches a mean CA-RMSD of
{ARM_MEAN:.2f} angstroms on the built chain. That is not state of the art, and I do not claim it
is. Inside the pipeline sits a real {Q_QUBITS}-qubit CVaR-VQE with {Q_PARAMS} parameters. I
verified it to machine precision, and then I measured what it adds. I could not detect an
accuracy contribution. What survived is the measuring instrument, and a chain of exact
statements about how a force field's structure reaches a circuit's gradients. Every number on
these slides has a file path in the notes, and I can send the file.

### also
The project I first described in the summer is not the project on these slides. Say that once,
early, without apology: the register is {Q_QUBITS} qubits, exact; the targets are {N_TARGETS} real
peptides; the sub-2 A target was not reached ({ARM_FRAC2:.1%} of targets are under 2 A).

## Slide 2 -- Why 9 to 16 residue peptides are hard

### spoken
Why are short peptides hard? Three reasons I can measure. First, the sequence says almost
nothing about the backbone at this length. Predicting the phi angle from the full sequence
context gives a mean error of {PHI_SEQ:.1f} degrees over {PHI_N_RES:,} residues. A model that
ignores the sequence and uses only the corpus average gives {PHI_BLIND:.1f} degrees. That is the
whole sequence signal for phi: a third of a degree. Second, the candidates I retrieve all share
the same mistake. {COMMON_MODE:.0%} of the pool's error is common to every member, so averaging
cannot remove it. Third, the only steep lever is the distance prior itself. Interpolating the
prior toward the true distances moves the answer by {PRIOR_SLOPE:.2f} angstroms per unit of
progress, on the point-cloud basis, and {PRIOR_GAMMA_FOR_3A:.1%} of the way would reach three
angstroms. That number is an oracle. It says where the accuracy lives, not how to get it. For
scale, a constant alpha-helix, which knows nothing about the sequence, scores {HELIX:.2f}
angstroms on the built chain. The whole pipeline is {HELIX_MINUS_ARM:.2f} angstroms better than
that. Everything is measured on {N_TARGETS} cluster-disjoint targets in five pinned folds,
paired per target, with a minimum detectable effect of {MDE_FACTOR} standard errors computed per
comparison, and a fold-clustered confidence interval next to every mean.

### also
The library behind retrieval: {CORPUS_PEPTIDES} peptides and {CORPUS_FRAGMENTS:,} protein
fragments (ledger L11). The sealed benchmark of 60 was spent once (slide 4). The identity leak:
{SELFCOPY_DEV} dev targets and {SELFCOPY_BENCH} benchmark targets carry a verbatim self-copy in
their fold model's training set; declared. Its dev price, re-derived this sprint (L44):
{LEAK_DEV_PRICE_FIT:+.4f} A on the lam = 0 chain (S10-4's figure reproduced) and
{LEAK_DEV_PRICE_ARM:+.4f} A on the built chain, neither clearing its MDE; the benchmark half is
bounded from the dev proxy without opening the benchmark (slide 4 notes).

## Slide 3 -- The pipeline, stage by stage

### spoken
Here is the pipeline. The sequence goes through ESM-2 into a distogram: a seventeen-bin
posterior over every pairwise CA distance, trained leave-fold-out. In parallel, BLOSUM62
retrieval pulls five hundred real protein windows of the right length from a library of
{CORPUS_PEPTIDES} peptides and {CORPUS_FRAGMENTS:,} fragments. Each window is scored by its Bayes
risk against the posterior, and the best seventy-five are kept. The quantum stage can select a
subset of the top one hundred and twenty-eight, but it is off in production; the {ARM_MEAN:.4f}
result never passes through it. The kept windows are averaged in coordinate space. That average
scores {AVG_MEAN:.4f} angstroms, but it is not a structure: its virtual CA-CA bonds are
{BOND_AVG_MEAN:.2f} angstroms on average, with a worst case of {BOND_AVG_WORST:.2f}. So the
average is projected onto a chain with ideal geometry. That built chain, at {ARM_MEAN:.4f}
angstroms, is the result I report; the projection costs {GAP_ARM_MINUS_AVG:.3f} angstroms and I
pay it because the chain is a real structure. An optional AMBER relaxation gives {FULL_MEAN:.4f}.
I will show later that it is a validity step, not an accuracy step. One design fact matters for
everything that follows. The output is an average, so it tracks the whole set, not its best
member. A selector that can only remove members from a ranked list cannot move it much.

### also
The 1-worker baseline and the 8-worker optimised run agree on every science quantity to
{SCIENCE_DELTA_MAX:.1e} (bit-identical science across the two arms). A fresh run of the
instrument over all {N_TARGETS} stored records disagrees with the record by
{REPRO_MAX_DISAGREE:.1f} on every basis. T030's universe holds {T030_N_WINDOWS:,} windows.

## Slide 4 -- Results on 126 held-out targets, with two overlays

### spoken
These are the results on the built chain. Mean {ARM_MEAN:.2f} angstroms, median
{ARM_MEDIAN:.2f}, best {ARM_MIN:.2f} on {ARM_BEST_PDB}, worst {ARM_MAX:.2f} on {ARM_WORST_PDB}.
{ARM_FRAC2:.0%} of targets are under two angstroms and {ARM_FRAC3:.0%} under three. The upper
figure is the best target, {T030_PDB}, a helix at {T030_ARM:.3f} angstroms. The lower figure is a
hard one, {HARD_PDB}, at {HARD_ARM:.2f}. For the figures I superposed the prediction onto the
native, so the pictures use the answer; nothing in the pipeline does. On {HARD_PDB} the
retrieval pool holds a window at {HARD_POOL_BEST:.2f} angstroms, but the score filter keeps
nothing under {HARD_TOPM_BEST:.2f}, so the average is built from the wrong windows. That is the
failure mode: the score, not the search. For scale on the same targets, a random
seventy-five-window subset scores {RANDOM75:.2f} on the point-cloud basis, against the
production point cloud of {AVG_MEAN:.2f}. The best single window in the pool would give
{POOL_BEST:.2f}; that is an oracle. On the sealed sixty-target benchmark, spent once, the full
system was {BENCH_DELTA:+.4f} angstroms against the shipped baseline, confidence interval
{BENCH_CI_LO:+.3f} to {BENCH_CI_HI:+.3f}, {BENCH_WL}. There is no validated improvement over the
baseline, and I say so.

### also
The overlays are ORACLE-superposed (Kabsch onto the native through `s12.instrument.superpose_batch`,
for the figure only); the RMSD printed on each figure was recomputed by `s26/pr_figures.py` and
equals the production record. The benchmark numbers are typed from claim C06 of
`s26/EXAMINATION.md`; no S26 lane opened `s9/final_report.json`, and the two benchmark means
({BENCH_FULL} full system, {BENCH_SHIPPED} shipped baseline) are asserted by a passing test.
Caveat that attaches to every benchmark figure (L44, L55, L58): {SELFCOPY_BENCH} benchmark targets
carry a verbatim self-copy in their fold model's training set; dev-proxy price
{LEAK_DEV4_BOUND_ARM:.3f} A; own-native envelope {LEAK_BOUND_ARM_MEANCI:.3f} A (mean CI) to
{LEAK_BOUND_ARM_WORST:.3f} A (worst target) on the built chain under assumption A2; class
{LEAK_CLASS} under every reading, against a benchmark CI half-width of {LEAK_BENCH_CI_HALF:.3f};
it cannot move the benchmark verdict either way.

## Slide 5 -- The quantum component: verified, then measured

### spoken
The quantum component is a real CVaR-VQE, and I checked it before I measured it. It is a
{Q_QUBITS}-qubit, {Q_LAYERS}-layer RY and CNOT circuit with {Q_PARAMS} parameters, simulated as
an exact statevector over {Q_DIM} candidate structures. Against a dense simulator written
independently from the gate list, the probabilities agree to {Q_SV_ERR}. The gradient is the
exact parameter-shift rule; against finite differences its cosine is {Q_PS_COS} with relative
error {Q_PS_REL}. Nothing is sampled in the trained path. The Hamiltonian is diagonal: the
rank-standardised scores of the candidates. Because the candidates are already sorted, it is
the same ladder on every target to within {Q_LADDER_WORST_FRAC:.1%} of its range. I proved that
the CVaR tail is always a subset of a prefix of that ranking: {Q_CELLS:,} adversarial cells,
{Q_VIOLATIONS} violations. So the circuit can delete a candidate from the classical top set; it
can never add one. Then I measured. The optimiser does train: it beats the best of
{Q_BEST_OF_N} untrained draws and closes {Q_GAP_MIN:.0%} to {Q_GAP_MAX:.0%} of the free-energy
gap. It stops {Q_KL_NATS:.3f} nats from its own Gibbs optimum, and the endpoint cannot tell:
{Q_CIRC_VS_GIBBS_X:.2f} of the minimum detectable effect. Against the plain argmin the selector
is {Q_VS_ARGMIN:+.4f} angstroms at {Q_VS_ARGMIN_X:.2f} of its MDE, on the selection basis. By my
own rule that is not a result. I claim no quantum advantage.

### also
The deployed table sets alpha = 1 on folds {Q_ALPHA1_FOLDS}, so {Q_SHARE_NO_TAIL:.1%} of targets
carry no tail constraint; there the objective is an entropy-regularised mean. Across 18 readout
arms, mean RMSD tracks the entropy of the readout weights at rho {Q_RHO_H:.2f}; the nine circuit
arms sit {Q_OFF_CURVE:+.4f} A from a curve fitted on the nine no-circuit arms. Against an exact
Boltzmann weighting at the same temperature the selector is {Q_VS_BOLTZ:+.4f} A. The withdrawn
"+0.113 A CVaR contribution" (S25 L5) is not on any slide. The generation lane's exact MPS twin
agrees with the dense circuit to {Q_MPS_ERR}; it is not the selector.

## Slide 6 -- Gradient variance at depth 3, and the algebra

### spoken
I measured gradient variance against width, at depth three, with exact gradients and no shot
noise, for {SWEEP_N_MIN} to {SWEEP_N_MAX} qubits. The fitted slopes, in log2 variance per qubit,
are {SLOPE_LIN:.3f} for the linear cost, {SLOPE_A025_T0:.3f} and {SLOPE_A01_T0:.3f} for the CVaR
tails, and {SLOPE_A1_T03:.3f} and {SLOPE_A025_T03:.3f} for the deployed free energy. For
reference, the deepest circuit I measured, at depth eight, decays with base
{GEO_BASE_DEPTH8:.3f} per qubit, which is {GEO_SLOPE_DEPTH8:.2f} on this scale. So the deployed
circuit decays more slowly than that, and the CVaR non-linearity flattens the decay further. I
do not call this the absence of a barren plateau, and I do not call a small gradient a barren
plateau either. It is a depth-three measurement with {Q_PARAMS} parameters against a Lie algebra
of dimension {DLA_SO128:,}. This sprint I computed that algebra exactly. At seven qubits the
ansatz generates the full so(128) from depth two; depth one is abelian, dimension
{DLA_N7_L1}. So nothing algebraic protects this ansatz at scale, and what I see is the shallow
regime. Depth at seven qubits saturates from four layers on. The draws per width are
{SWEEP_DRAWS_MAX} down to {SWEEP_DRAWS_MIN}, so each variance carries a relative error near ten
to sixteen percent, and every slope is a seven-point fit.

### also
The symbolic Lie closure agrees with the dense SVD rank on {DLA_NUMERIC_AGREE} cells at n = 4, 5.
The CVaR-to-linear variance ratio at n = 7 is {RATIO_N7_A025:.2f} (alpha 0.25) and rises to
{RATIO_N13_A025:.2f} at n = 13. Depth sweep at n = 7: {DEPTH_L1:.2e} at L = 1, {DEPTH_L3:.2e} at
the deployed L = 3, {DEPTH_L4:.2e} at L = 4, {DEPTH_L12:.2e} at L = 12. The A2 figure on the
right is lane Q's, from `s26/results/q_dla.json`; the left figure is rebuilt from
`s25/results/q_plateau.json` by `s26/pr_figures.py`.

## Slide 7 -- Where the remaining accuracy lives

### spoken
Where does the remaining accuracy live? The ladder on the left is measured on one instrument,
with the basis written on every bar. Two native-free controls sit above the pipeline: a
constant helix at {HELIX:.2f} and a sequence-only torsion predictor at {TORS:.2f}, both built
chains. Below the pipeline are oracles that read the native. The best single window in the pool
is {POOL_BEST:.2f} angstroms, and the best inside the shipped top-75 is {TOPM_BEST:.2f}. Those
windows exist, and no native-free ranker finds them; the record closed that route at every
level. The physics does not help either. As selectors, the Legacy potential and AMBER are worse
than a random subset by {LEGACY_VS_RANDOM:+.3f} and {AMBER_VS_RANDOM:+.3f} angstroms, five of
five folds, on the point-cloud basis. And relaxing the built chain with AMBER costs
{C3_AMBER_VS_NONE:+.4f} angstroms; a random move of the same size costs only
{C3_RANDOM_VS_NONE:+.4f}, so the physics is worse than noise, and the difference sits in the
Type-M zone, so its sign is what I have measured, not its size. The one steep lever is the
distance prior. A perfect prior through the same pipeline reaches {PRIOR_PERFECT:.2f} on the
point-cloud basis, from {PRIOR_GAMMA0:.2f}. And the pool's error is {COMMON_MODE:.0%}
common-mode, so nothing downstream of retrieval can remove most of it. The conclusion is plain.
Search, selection, ranking and physics are closed by measurement. If accuracy moves, it moves
through a better prior.

### also
Seven-configuration suite, point-cloud basis (`s25/results/phys_suite.json`): distogram
{SUITE_DIST:.3f}, AMBER+distogram {SUITE_AMBER_DIST:.3f}, Legacy+distogram
{SUITE_LEGACY_DIST:.3f}, all three {SUITE_LAD:.3f}, random-75 {RANDOM75:.3f}, Legacy+AMBER
{SUITE_LEGACY_AMBER:.3f}, Legacy {SUITE_LEGACY:.3f}, AMBER {SUITE_AMBER:.3f}. The same rows on
the built chain (`results/summary/leaderboard.json`): {LB_DISTOGRAM:.3f}, {LB_AMBER_DISTOGRAM:.3f},
{LB_LEGACY_DISTOGRAM:.3f}, {LB_LEGACY_AMBER_DISTOGRAM:.3f}, {LB_LEGACY_AMBER:.3f}, {LB_LEGACY:.3f},
{LB_AMBER:.3f}; the ranking is unchanged. Fold CIs: Legacy vs random {LEGACY_VS_RANDOM_CI},
AMBER vs random {AMBER_VS_RANDOM_CI}. Torsion-space ceiling {TORSION_CEILING:.3f} (ORACLE, k = 4);
distance geometry from the true distances {DISTGEO_TRUE:.3f} (ORACLE); best window in the whole
library {UNIVERSE_BEST:.3f} (ORACLE). AMBER relaxation (L39, Adversary L46): vs the built chain
{C3_AMBER_VS_NONE:+.4f} [fold CI {C3_AMBER_VS_NONE_CI}], {C3_AMBER_VS_NONE_X:.2f} x MDE; vs a
move of the same size toward a random pool member {C3_AMBER_VS_MEMBER:+.4f} [{C3_AMBER_VS_MEMBER_CI}];
vs a random move of the same size {C3_AMBER_VS_RANDOM:+.4f} [{C3_AMBER_VS_RANDOM_CI}],
{C3_AMBER_VS_RANDOM_FOLDS}/5 folds, {C3_AMBER_VS_RANDOM_X:.2f} x MDE (Type-M zone: sign measured,
magnitude an upper bound). A validity step on 124 of 126; on 2BP4 and 9KAR the relaxation breaks a
virtual CA-CA bond and 9KAR does not converge. Physics as a steric reject filter at 1e4 kcal/mol
(L43, L54, point cloud): {REJECT_R_VS_ANCHOR:+.3f} A worse than the shipped top-75 (median
{REJECT_R_VS_ANCHOR_MEDIAN:+.3f}: near zero on the median target, the harm is tail-carried;
Type-M) and {REJECT_R_VS_RANDR:+.3f} vs rejecting the same count at random. A free calibration
flag, not a lever (L53): how far the relaxation moves the chain predicts its error, Spearman
{STRAIN_RHO_MOVED:+.3f} partial on n and Rg, fold CI {STRAIN_RHO_MOVED_CI}; quartile means
{STRAIN_QUARTILE_MEANS} A.

## Slide 8 -- Direction A: ADAPT-VQE (PENDING)

### spoken
PENDING. This slide waits for `s26/PROPOSAL_A.md`.

## Slide 9 -- Direction B: a learned folding model as the prior (PENDING)

### spoken
PENDING. This slide waits for its build step.

## Slide 10 -- Direction C: learn a better distance prior; physics for validity only

### spoken
Direction C says: learn a better distance prior, and keep physics honest. The prior bounds the
accuracy, so the plan has five parts. C1 reproduced the closures this rests on, to the third
decimal. A learned ranker on real features moves from {C1_S12_REAL_N8:.3f} to
{C1_S12_REAL_FULL:.3f} angstroms as its training set grows; the same model with a leaked label
reaches {C1_S12_LEAK_FULL:.3f}. A perfect ranker inside the shipped top-25 returns
{C1_S17_BAND_BEST:.3f} against {C1_S17_RAND:.3f} for a random member. So learning a ranker is
closed, and the proposal is not that. C2 is a ladder of {C2_N_RUNGS} prior inputs, from no
language model to the full 650-million-parameter embedding. Its falsifier was fixed before the
first run: a rung must beat the shipped posterior on the built chain by more than its own MDE,
with a fold-clustered interval excluding zero on five of five folds. Seven rungs are trained;
the evaluations are running. C3 asked whether AMBER relaxation is an accuracy step. It is not.
The relaxation moves the chain {C3_MAG:.3f} angstroms and costs {C3_AMBER_VS_NONE:+.4f}; a
random move of exactly the same size costs {C3_RANDOM_VS_NONE:+.4f}. Physics is worse than
noise, and its displacement points slightly away from the native. So refine-with-physics is
kept as a validity step only: it brings {C3_N_CONVERGED} of 126 chains below a thousand
kilocalories per mole, though on two targets it breaks a virtual bond. C4, routers, and C5, a
learned common-mode correction, are pre-registered and pending. The verdict on C is pending the
ladder.

### also
Verdict form (campaign prompt): KEEP AS STATED / KEEP WITH EDITS / REPLACE; lane P's
`s26/PROPOSAL_C.md` is DRAFT and its verdict is pending C2 to C5. The C1 learning curve spans
{C1_S12_REAL_N8:.4f} (n = 8) to {C1_S12_REAL_FULL:.4f} (full) on the weighted-average basis;
the leaked label reaches {C1_S12_LEAK_N8:.4f} at n = 8. Perfect ranker in the top-25
{C1_S17_BAND_BEST:.4f}, random member {C1_S17_RAND:.4f}, distogram argmin {C1_S17_DIST_PICK:.4f}
(single window). C2 anchor (L56, L57): selection max abs {C2_ANCHOR_SEL_MAXABS:.1e} and point
cloud {C2_ANCHOR_CLOUD_MAXABS:.1e} against the production cache; the built chain through the
ladder's own path is {C2_ANCHOR_ARM:.4f} A, the leaderboard-rebuild basis, so every C2 contrast
is paired on that basis and the 0.002 A offset to {ARM_MEAN:.4f} is a basis difference, not an
effect. Rungs with all five fold checkpoints on disk at build time: {C2_RUNGS_TRAINED}. C3 (L39,
L46): AMBER vs built chain {C3_AMBER_VS_NONE:+.4f} [fold CI {C3_AMBER_VS_NONE_CI}], vs random
{C3_AMBER_VS_RANDOM:+.4f} [{C3_AMBER_VS_RANDOM_CI}] at {C3_AMBER_VS_RANDOM_X:.2f} x MDE (Type-M:
sign measured, magnitude an upper bound), vs toward-member {C3_AMBER_VS_MEMBER:+.4f}
[{C3_AMBER_VS_MEMBER_CI}]; ORACLE cosine of the displacement with the true residual
{C3_COS:+.3f}, positive on {C3_COS_FRAC_POS:.1%} of targets; before relaxation
{C3_FRAC_E0_OVER_1E4:.1%} of built chains sit above 1e4 kcal/mol; after it {C3_N_CONVERGED} of
126 converge (9KAR ends at {C3_E1_MAX:.0f} kcal/mol) and the virtual CA-CA bond stretches from
{C3_BOND_BEFORE:.3f} to {C3_BOND_AFTER:.3f} A. The free calibration flag (L53): the distance the
relaxation moves a chain predicts its error, Spearman {STRAIN_RHO_MOVED:+.3f} partial on n and Rg;
never a selector. Stage 2 of C3 repeats the seven contrasts on the best C2 rung when lane P
delivers it.

## Slide 11 -- Goal and ask

### spoken
My goal is two things. First, find out whether a better distance prior is obtainable from
inputs this machine can compute, on the instrument I already have, and report it with a
pre-registered falsifier and fold-clustered intervals. Second, publish the trainability chain:
the locality theorem, the Pauli spectrum, the kernel, the width sweep, and this sprint's algebra
and product-circuit results. The ask is an independent study with those two threads, and one
question for you: what does an adaptive ansatz do when the target Hamiltonian is not diagonal,
and does the spectrum-times-kernel prediction extend to that case? Three things I will not do.
I will not claim a quantum advantage. I will not call a small gradient a barren plateau. And I
will not report an effect below its own minimum detectable size.

### also
DRAFT (until the coordinator's verdict entry in `s26/LEDGER.md`): if asked "which would you do
first?", the evidence favours direction C in its learn-the-prior form as the accuracy study,
because the prior is the only lever the record measured as steep ({PRIOR_SLOPE:.2f} A per unit,
`s24/results/priorladder.json`) and every selection, ranking, physics and search lever is closed
by measurement; paired with direction A as a scientific study rather than an accuracy study,
because on this diagonal Hamiltonian an ADAPT-grown circuit is a product circuit
(`s26/results/q_var.json`, L35) and the set-equality theorem bounds any readout by the classical
prefix. Direction B as stated cannot be measured on this machine (`s26/results/b1_feasibility.json`).
