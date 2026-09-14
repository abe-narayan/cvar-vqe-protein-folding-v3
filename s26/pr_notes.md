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
Type-M) and {REJECT_R_VS_RANDR:+.3f} vs rejecting the same count at random. A free calibration flag, not a lever, restated per L121 / L123: the pool's own disagreement (the shipped
top-75's pairwise CA-RMSD spread, native-free, available before any relaxation runs) predicts the
emitted chain's error at Spearman {SPREAD_RHO_PARTIAL:+.3f} partial on n and Rg, fold CI
{SPREAD_RHO_PARTIAL_CI}, {SPREAD_FOLDS}/5 folds; the relaxation's displacement tracks that disagreement at
rho {RHO_MOVED_SPREAD:.2f} and adds nothing given it ({MOVED_GIVEN_SPREAD:+.3f}, iid CI {MOVED_GIVEN_SPREAD_CI},
permutation p {MOVED_GIVEN_SPREAD_P:.2f}; `s26/results/a_strain_vs_spread.json`). L53's quartile means by
the displacement, {STRAIN_QUARTILE_MEANS} A, stand as a presentable form of the phenomenon; L53's
"first native-free quantity above 0.4" and "how far the relaxation moves the chain predicts its
error" are retracted (L123).

## Slide 8 -- Direction A: let the circuit grow (qubit-ADAPT-VQE). Verdict: REPLACE

### spoken
We tested letting the circuit grow itself. qubit-ADAPT-VQE starts from one layer of
single-qubit rotations and adds, one at a time, whichever operator would lower the objective
fastest. We ran it on all {A1_N} development peptides with the deployed objective, seed and
readout, and compared the built chains pairwise. The answer is no change: {A1_L2_EFFECT:.3f} to
{A1_V_EFFECT:.3f} angstroms, a quarter to a third of what the comparison can resolve, with
fold-clustered intervals straddling zero. All {A1_N_ARMS} variants point the same way, but they
are one observation, not twelve, and none clears its bar; the resolution is {A1_L2_MDE:.2f}
angstroms. Why? At the deployed setting the optimum the circuit is asked to reach is a product
state: seven single-qubit rotations represent it exactly, to {A1_KL_RY7_MEAN:.4f} nats. When
ADAPT is offered entangling operators it does append them, on {A1_APPENDED_V_LBFGS} to
{A1_APPENDED_ADAM} of the {A1_N_ALPHA1} targets, but they are inert: together they lower the
objective by less than 0.001 nats, their angles stay below 0.02 radians under L-BFGS, and the
state remains a product state to {A1_KL_TO_PRODUCT_MAX:.4f} nats. The deployed 21-parameter
circuit stops {A1_KL_FIXED_MEAN:.2f} nats short of that same optimum, and reaching it exactly
moves the emitted structure by those 0.014 to 0.022 angstroms, a third of the resolution. Its
Lie algebra is already the full real algebra from depth two: shallowness, not structure. So
Proposal A is replaced. The adaptive circuit is a correct, tested tool; it confirmed the
diagnosis instead of curing it. What we publish is on the next slide.

### also
Verdict REPLACE (`s26/PROPOSAL_A.md`, lane Q; L68; accepted by the coordinator in L69, subject
to the Adversary's check of L68): the mechanism is absent, not weak. A1 basis: built chain,
`rmsd_q_synth` (the production projection of the weighted average over the 128 candidates) on
both sides; the deployed selector's arm is {A1_FIXED_MEAN:.4f} A on that basis against the
production top-75 arm {ARM_MEAN:.4f} A, and the whole quantum synthesis is worth
{QSYNTH_VS_ARM:+.4f} A ({QSYNTH_VS_ARM_X:.2f} x MDE) against the classical arm, so the resolution
is about four times the component's own footprint. Primaries (`s26/results/a1_stats.json`):
ADAPT 2-local pool minus fixed {A1_L2_EFFECT:+.4f} (SE {A1_L2_SE:.4f}, MDE {A1_L2_MDE:.4f},
{A1_L2_X:.2f} x, fold CI {A1_L2_CI}, {A1_L2_W}W/{A1_L2_L}L, {A1_L2_FOLDS}/5 folds); ADAPT Tang pool
minus fixed {A1_V_EFFECT:+.4f} (SE {A1_V_SE:.4f}, MDE {A1_V_MDE:.4f}, {A1_V_X:.2f} x, fold CI
{A1_V_CI}, {A1_V_W}W/{A1_V_L}L). All {A1_N_ARMS} ADAPT arms (two pools, two optimisers, 7 / 14 / 21
parameters) are between {A1_ARMS_EFF_MIN:+.4f} and {A1_ARMS_EFF_MAX:+.4f} A at {A1_ARMS_X_MIN:.2f}
to {A1_ARMS_X_MAX:.2f} x MDE; on the selection readout {A1_SEL_EFF_MIN:+.4f} to {A1_SEL_EFF_MAX:+.4f}
A at up to {A1_SEL_X_MAX:.2f} x MDE with 45 to 59 exact ties (the shape of S25's headline: a
direction the instrument cannot size). The twelve arms are one observation, not twelve: their per-target delta vectors correlate
at {A1_ARM_CORR_MEAN:.3f} on average (min {A1_ARM_CORR_MIN:.3f}; Adversary L70 caveat 1). The 21
parameters are a budget: the Adam primaries realise {A1_P21_ADAM_MIN} (repeats of one rotation
merge), L-BFGS {A1_P21_LBFGS_MIN} to {A1_P21_LBFGS_MAX} (L70 caveat 2). Controls: the exact Gibbs
state in place of the circuit {A1_GIBBS_EFFECT:+.4f} ({A1_GIBBS_X:.2f} x) overall, {A1_GIBBS_ALPHA1:+.4f}
on the {A1_N_ALPHA1} alpha = 1 targets where it is the optimum and {A1_GIBBS_ALPHA025:+.4f} on the
{A1_N_ALPHA025} alpha = 0.25 targets where it is not (L70 caveat 3); the matched-entropy random Hamiltonian
{A1_RANDH_EFFECT:+.4f} ({A1_RANDH_X:.2f} x, fold CI {A1_RANDH_CI}, 5/5 folds). Product-state facts
from the {A1_N_RECORDS} per-target records (`s26/results/a1/*.json`, no native): KL(Gibbs ||
product of marginals) mean {A1_KL_PRODUCT_MEAN:.1e}, max {A1_KL_PRODUCT_MAX:.1e}; on the
{A1_N_ALPHA1} alpha = 1 targets the fixed circuit's KL to Gibbs is {A1_KL_FIXED_MEAN:.4f} (max
{A1_KL_FIXED_MAX:.3f}; S25's 0.902 reproduced) and the 7-parameter RY layer's {A1_KL_RY7_MEAN:.1e}
(max {A1_KL_RY7_MAX:.1e}). Growth at alpha = 1, L75's reconciliation of my L73 flag (the two "no operator selected"
sentences of L68 and `s26/PROPOSAL_A.md` are retracted there and corrected by an addendum): the
growth halts by the eps = 1e-3 gradient criterion on {A1_LBFGS_STOP_EPS} of 156 (pool, target)
cells under L-BFGS and at the 21-parameter cap under Adam, but operators ARE appended, on
{A1_APPENDED_V_LBFGS} (pool V) and {A1_APPENDED_L2_LBFGS} (pool L2) of the {A1_N_ALPHA1} targets
under L-BFGS, every one a multi-qubit string, and on {A1_APPENDED_ADAM} of {A1_N_ALPHA1} under
Adam; they are inert: the free energy moves by at most {A1_LBFGS_DF_ABSMAX:.1e} nats under L-BFGS
and {A1_DF_ABSMAX_ALL:.1e} under Adam (median {A1_LBFGS_DF_MEDIAN:.1e} under L-BFGS), the appended
angles stay at or below {A1_LBFGS_ANGLE_MAX:.3f} rad under L-BFGS, and the state stays a product
state to KL {A1_KL_TO_PRODUCT_MAX:.1e} (`s26/results/a1/*.json :: adapt/*/sequence, trace, theta;
arms/adapt*_P21/kl_to_product`). The ideal-ladder statements of L27 and L35 (nothing appended at
alpha = 1) stand: there E is exactly affine and every pool gradient is exactly zero; on a real
target tie-averaging leaves a residual gradient just above eps. The A4 figure
(slide 8) is lane Q's; the product-circuit reading rests on the adapt sets of
`s26/results/q_dla.json` (L47). Slide 9 carries what is published in A's place.
L119's bootstrap (`s26/results/q_var_boot.json`, every q_var.json row reproduced at relative
{A4_BOOT_REPRO}): the alpha = 0.25 grown-L2 minus fixed slope difference is {A4_BOOT_L2_A025_DIFF:+.3f},
95% CI {A4_BOOT_L2_A025_DIFF_CI} (includes zero, so "decays like the fixed ansatz" now carries a stored
interval); the alpha = 1 grown minus fixed differences are {A4_BOOT_A1_DIFF_LO:+.2f} to {A4_BOOT_A1_DIFF_HI:+.2f}
with every interval excluding zero ("no decay" is measured, not read); the A4 figure on this slide
carries those intervals as error bars. A3 (L125, `s26/results/a3_stats.json`, `a3_property.json`): a
target-dependent Hamiltonian (the raw standardised score at a per-target temperature matched to the
deployed entropy) makes the trained states target-dependent, {A3_DISTINCT_TMATCH} of {A3_DISTINCT_N}
distinct against {A3_DISTINCT_ZRANK} under the deployed rank ladder, and the readout does not notice:
{A3_TMATCH_EFFECT:+.4f} A on the built chain, {A3_TMATCH_X:.2f} x MDE, fold CI {A3_TMATCH_CI}.

## Slide 9 -- Direction B (REPLACE), and what we publish: the trainability paper

### spoken
Proposal B wanted to swap our distance prior for a large pretrained folding model. We measured
whether that can even run here: it cannot. ESMFold needs about {B1_RESIDENT_GB:.1f} gigabytes
resident and a GPU-era dependency chain, and this machine had {B1_HEADROOM_GB:.1f} gigabytes of
headroom and no CUDA. So we measured the part of that idea that does fit. Our prior already
reads a 650-million-parameter protein language model; removing it costs {C2_NOESM_ARM:.2f}
angstroms on the built chain and {C2_NOESM_SEL:.2f} on selection, five folds out of five.
An 8-million-parameter model loses all of that, and two thirds of the selection value sits in
the model's contact head. We also asked whether we could tell, without the answer, which
targets the pipeline beats a sequence-only predictor on, to route between them. We cannot: a
classifier is at chance, and the only predictable thing is how much a constant helix loses,
set by how strand-like the retrieved pool is. Our verdict is replace. The replacement is the
trainability paper, and it is what we publish: chain geometry
fixes which residues a distance can depend on, an exact theorem; that fixes the energy's Pauli
spectrum once the raw force field's clash spike is conditioned away; and the spectrum times
the circuit's own kernel predicts the measured gradient variance with no free parameter, a
median ratio of {PAULI_RATIO_LEGACY:.3f} over {PAULI_N_CELLS} cells. This sprint added the
algebra, the product-state target and the ADAPT null. The one untested input, a larger
language model, needs a larger machine.

### also
Verdict REPLACE, final (`s26/PROPOSAL_B.md`, commit a88ea259; coordinator L117; Adversary L120 STANDS). B1 (L13):
"{B1_VERDICT}"; {B1_DOWNLOAD_GB:.2f} GB of checkpoints not on disk; openfold and omegaconf
absent, Python 3.13 vs <= 3.9; {B1_RESIDENT_GB:.2f} GB resident against {B1_HEADROOM_GB:.1f}
GB, fp16 everywhere {B1_FP16_GB:.1f}. B2, the feasible-scale ladder (L62, L63, L99/L101; every
rung through one path, built chain on the rebuild basis {C2_ANCHOR_ARM:.4f}, L57): noesm
{C2_NOESM_ARM:+.3f} A (SE {C2_NOESM_ARM_SE:.3f}, MDE {C2_NOESM_ARM_MDE:.3f}, {C2_NOESM_ARM_X:.2f} x,
fold CI {C2_NOESM_ARM_CI}, {C2_NOESM_ARM_FOLDS}/5, {C2_NOESM_ARM_W}W/{C2_NOESM_ARM_L}L; selection
{C2_NOESM_SEL:+.3f} at {C2_NOESM_SEL_X:.2f} x, fold CI {C2_NOESM_SEL_CI}); conly
{C2_CONLY_ARM:+.3f} ({C2_CONLY_ARM_X:.2f} x, fold CI {C2_CONLY_ARM_CI}, {C2_CONLY_ARM_FOLDS}/5;
selection {C2_CONLY_SEL:+.3f} at {C2_CONLY_SEL_X:.2f} x); esm8m {C2_ESM8M_ARM:+.3f}
({C2_ESM8M_ARM_X:.2f} x, fold CI {C2_ESM8M_ARM_CI}, {C2_ESM8M_ARM_FOLDS}/5; selection
{C2_ESM8M_SEL:+.3f} at {C2_ESM8M_SEL_X:.2f} x); the retrained shipped recipe (pca32) is the
identity on {C2_PCA32_ARM_T} of 126. Isolations: the contact head alone (conly minus noesm)
{B2_CONLY_MINUS_NOESM_SEL:+.3f} on selection (L63: 1.00x MDE, fold CI [-0.374, -0.037], 4/5) and
{B2_CONLY_MINUS_NOESM_ARM:+.3f} on the built chain; the 8M model minus no ESM
{B2_ESM8M_MINUS_NOESM_ARM:+.3f} on the built chain (L99: 0.15x MDE). B3 (L106, L107,
`s26/results/p_b3.json`): sign classifier balanced accuracy {B3_BACC_TORS:.3f} against the
sequence-only predictor (permutation null 95th percentile {B3_NULL95_TORS:.3f}) and
{B3_BACC_HELIX:.3f} against the constant helix (null {B3_NULL95_HELIX:.3f}); the SIZE of the gain
over the helix is partly predictable, held-out R2 {B3_R2_HELIX:.2f}, squared-error reduction
{B3_MSE_X_HELIX:.2f} x MDE, fold CI {B3_MSE_CI_HELIX}, carried by the pool's strand content
(Pearson {B3_RHO_SSE_HELIX:+.3f} between the top-75 strand fraction and arm minus helix); on the
ORACLE stratum FAIL18 the sequence-only predictor wins by {B3_FAIL18_TORS:+.3f} A and nothing
native-free locates it. The two-minute script above is lane P's, with the numbers read from
the artefacts. The trainability paper (`s26/PROPOSAL_B_REPLACEMENT.md`, lane Q): C1 locality
theorem (`s13/results/qarch_locality_geom.json`); C2 raw AMBER's delta-spike spectrum
(`s13/results/walsh_xval.json`); C4 spectrum x kernel, median measured/predicted
{PAULI_RATIO_LEGACY:.4f} over {PAULI_N_CELLS} cells (`s13/results/geo_pauli.json`); C5 the kernel
flat in Pauli weight at depth >= 3 (`s13/results/geo_kernel.json`); C7 the width sweep, slopes
{SLOPE_LIN:.3f} / {SLOPE_A025_T0:.3f} / {SLOPE_A01_T0:.3f} / {SLOPE_A1_T03:.3f} / {SLOPE_A025_T03:.3f}
(`s25/results/q_plateau.json`); C8 the DLA, dim {DLA_SO128} = so(128) from depth 2 at n = 7
(`s26/results/q_dla.json`, slide 6); C9 the grown circuits (`s26/results/q_var.json`, slide 8);
C10 the optimiser trains and the readout cannot tell, {Q_KL_NATS:.3f} nats, {Q_CIRC_VS_GIBBS_X:.2f} x
MDE (`s25/results/q_gibbs.json`, `q_alpha.json`); C11 the set-equality theorem, {Q_CELLS:,} cells,
{Q_VIOLATIONS} violations (`s25/results/q_verify.json`); C12 the product-state target, KL to the
product of marginals at most {A1_KL_PRODUCT_MAX:.1e} on {A1_N_RECORDS} targets and the A1 null
{A1_V_EFFECT:+.4f} A at {A1_V_X:.2f} x MDE (`s26/results/a1/*.json`, `a1_stats.json`, slide 8).
Venue statement: every number is noiseless exact simulation; a submission would add a noise
model, width beyond n = 13, a 2-design control, seeds and error bars on every variance. Scope of
the checking (L120, L122): the S26 additions are Adversary-checked (A2 L45; A4 L47 with the L119
intervals, the alpha = 0.25 grown minus fixed slope difference {A4_BOOT_L2_A025_DIFF:+.3f}, 95% CI
{A4_BOOT_L2_A025_DIFF_CI}, and the alpha = 1 grown minus fixed differences {A4_BOOT_A1_DIFF_LO:+.2f} to
{A4_BOOT_A1_DIFF_HI:+.2f} with every interval excluding zero; the product-state fact L70); the S13
inputs (C1, C2, C4, C5 above) are cited from their artefacts without an S26 re-check, and the S13
Pauli mean weights 2.236 / 3.015 are not on any slide. The lost S7 ESM number (-0.288 A, L11) is
not on the slide; its re-measurement is C2_NOESM_SEL.

## Slide 10 -- Direction C: learn a better distance prior (KEEP WITH EDITS); physics for validity only

### spoken
Our predictor's accuracy is set by its distance prior, and last sprint we measured that moving
that prior {PRIOR_GAMMA_FOR_3A:.1%} of the way toward the truth would take us under three
angstroms. So this sprint we asked whether any input our machine can compute makes the prior
better. We retrained it nine ways: without the language model, with only its contact head,
with a smaller language model, with more of its embedding, a bigger network, a triangle-update
architecture, and a mix with the retrieval pool's own statistics. First, the
retrained recipe reproduces the pipeline exactly on all {C2_PCA32_ARM_T} targets, so every
difference is real. Second, the language model is worth {C2_NOESM_ARM:.2f} angstroms on the
built chain and {C2_NOESM_SEL:.2f} on selection, five folds out of five, and an
8-million-parameter model carries none of that. Third, nothing beats the shipped prior: every
other variant lands within {C2_NULL_EFF_MAX:.2f} angstroms of it, below what 126 targets can
resolve. Routing the set size on new features failed again. The physics step is worse than a
random move of its own size, so it stays as a validity check only: it
turns {VAL_CLASH_TARGETS_BEFORE} emissions with a sub-2-angstrom heavy-atom overlap into
{VAL_CLASH_TARGETS_AFTER}, at a price of {C3_AMBER_VS_NONE:.3f} angstroms of accuracy,
{VAL_BOND_STRAIN_AFTER:.1%} bond strain, and a broken virtual bond on two targets. Our verdict
is keep with edits: the prior is the lever, its derivative is steep, its inputs on this
machine are flat, and the honest next step is a language model too large for this box.

### also
Verdict KEEP WITH EDITS, final (`s26/PROPOSAL_C.md` with addenda 1 and 2, 03:54; coordinator L117; Adversary L120 STANDS). The
four edits: (1) C2 "learn a better prior from better inputs" becomes "the prior's inputs are
measured out at n = 126": the ESM-2 650M channel is the whole of the channel this instrument
can see ({C2_NOESM_ARM:+.3f} A built chain, {C2_NOESM_SEL:+.3f} A selection, {C2_NOESM_ARM_FOLDS}/5
folds), and no reduction, expansion, capacity, joint-consistency or pool-histogram variant of it
moves the endpoint beyond {C2_NULL_EFF_MAX:.2f} A against MDEs of {C2_NULL_MDE_MIN:.2f} to
{C2_NULL_MDE_MAX:.2f} A; the lever that stays open is a larger language model, which does not
fit this box (B1). (2) C3 "refine with physics" becomes "keep AMBER as a validity step only".
(3) C4 "route the set size and scale" becomes "closed": {C4_N_M_ROUTERS} m* routers on four
feature blocks no previous router used, {C4_M_HARMFUL} pointing the harmful way, the largest
{C4_M_MAX_EFF:+.3f} A at {C4_M_MAX_X:.2f} x MDE, none clearing its MDE; {C4_N_S_ROUTERS} s*
routers, all with the wrong sign (rho {C4_S_RHO_MIN:+.2f} to {C4_S_RHO_MAX:+.2f}), costing
{C4_S_EFF_MIN:+.3f} to {C4_S_EFF_MAX:+.3f} A against s = 1 (`s26/results/p_c4.json`, L110, L115);
C5 (predict and subtract the common mode) is CLOSED (`s26/PROPOSAL_C.md` addendum 2, 03:54; ledger
L132; `s26/results/p_c5.json` complete = {C5_COMPLETE}, 126/126): subtracting a predicted common-mode
correction on held-out folds is null to harmful on the built chain (rebuild basis): GLOBAL in the
coordinate frame {C5_GLOBAL_R2:+.4f} A ({C5_GLOBAL_R2_X:.2f} x MDE, fold CI {C5_GLOBAL_R2_CI},
{C5_GLOBAL_R2_VERDICT}), GLOBAL in distance space {C5_GLOBAL_R1:+.4f} ({C5_GLOBAL_R1_X:.2f} x,
{C5_GLOBAL_R1_VERDICT}), RIDGE in distance space on 45 native-free features {C5_RIDGE_R1:+.4f}
({C5_RIDGE_R1_X:.2f} x, fold CI {C5_RIDGE_R1_CI}, {C5_RIDGE_R1_FOLDS}/5 folds, {C5_RIDGE_R1_VERDICT}),
indistinguishable from a random correction of the same size ({C5_RANDOM_R1:+.4f}, {C5_RANDOM_R1_X:.2f} x);
the ORACLE ceilings are {C5_ORACLE_R1:+.3f} A (distance space) and {C5_ORACLE_R2:+.3f} A (coordinate
frame): the common mode is most of the error and nothing native-free touches it (S16, S19 L14, S24 L7
confirmed on the same operator). The raw rung (1280-d input) is NOT RUN: {RAW_FOLDS_TRAINED} of 5 fold
models trained, no evaluation (L133); its inputs are bracketed by the null pca32f and pca128 rungs. (4) "The prior's derivative is steep
and the prior's inputs are flat": every achievable rung moves at cosine 0.2 to 0.5 to the
native's direction (S25 L12). C1 (L26): closures reproduced, {C1_S12_REAL_N8:.4f} to
{C1_S12_REAL_FULL:.4f} real vs {C1_S12_LEAK_N8:.4f} to {C1_S12_LEAK_FULL:.4f} leaked;
{C1_S17_BAND_BEST:.4f} vs {C1_S17_RAND:.4f}. C2, all nine evaluated rungs, built chain (rebuild
basis {C2_ANCHOR_ARM:.4f}), effect (x MDE, fold CI, folds): noesm {C2_NOESM_ARM:+.3f}
({C2_NOESM_ARM_X:.2f} x, {C2_NOESM_ARM_CI}, {C2_NOESM_ARM_FOLDS}/5, WORSE, Type-M); conly
{C2_CONLY_ARM:+.3f} ({C2_CONLY_ARM_X:.2f} x, {C2_CONLY_ARM_CI}); esm8m {C2_ESM8M_ARM:+.3f}
({C2_ESM8M_ARM_X:.2f} x, {C2_ESM8M_ARM_CI}, {C2_ESM8M_ARM_FOLDS}/5, WORSE, Type-M); pca32 identity
({C2_PCA32_ARM_T} ties); pca32f {C2_PCA32F_ARM:+.3f} ({C2_PCA32F_ARM_X:.2f} x, {C2_PCA32F_ARM_CI});
pca128 {C2_PCA128_ARM:+.3f} ({C2_PCA128_ARM_X:.2f} x, {C2_PCA128_ARM_CI}); wide {C2_WIDE_ARM:+.3f}
({C2_WIDE_ARM_X:.2f} x, {C2_WIDE_ARM_CI}); pairnet {C2_PAIRNET_ARM:+.3f} ({C2_PAIRNET_ARM_X:.2f} x,
{C2_PAIRNET_ARM_CI}); mix identity ({C2_MIX_ARM_T} ties, lam* = 0 on 5/5 folds); raw: NOT RUN, {RAW_FOLDS_TRAINED} of 5
folds trained, not evaluated (L133). No rung is better than the shipped prior on the built chain ({C2_ANY_BETTER}).
C3 (L39, L46, L87, L100; `s26/C3_RESULT.md` addenda 1 to 4): AMBER vs the built chain
{C3_AMBER_VS_NONE:+.4f} (fold CI {C3_AMBER_VS_NONE_CI}, {C3_AMBER_VS_NONE_X:.2f} x); vs a random
move of its own size {C3_AMBER_VS_RANDOM:+.4f} ({C3_AMBER_VS_RANDOM_CI}, {C3_AMBER_VS_RANDOM_X:.2f} x,
Type-M: sign measured, size an upper bound; replication {C3_REP_AMBER_VS_RANDOM:+.4f}); vs a
same-size move toward a random pool member {C3_AMBER_VS_MEMBER:+.4f} ({C3_AMBER_VS_MEMBER_CI};
replication {C3_REP_AMBER_VS_MEMBER:+.4f}); the toward-member control itself improves the
chain, {C3_MEMBER_VS_NONE:+.4f}, replicated at {C3_REP_MEMBER_VS_NONE:+.4f}
{C3_REP_MEMBER_VS_NONE_CI}, a statement about the projection's cost (S16 L27), not about
physics; ORACLE cosine {C3_COS:+.3f}. The validity axis (`s26/results/ph_validity.json`, L100):
heavy-atom pairs below 2.0 A {VAL_CLASH_BEFORE:.3f} to {VAL_CLASH_AFTER:.3f} per target
({VAL_CLASH_TARGETS_BEFORE} targets to {VAL_CLASH_TARGETS_AFTER}); closest pair
{VAL_MINHEAVY_BEFORE:.2f} to {VAL_MINHEAVY_AFTER:.2f} A; bond strain {VAL_BOND_STRAIN_AFTER:.1%},
angle strain {VAL_ANGLE_STRAIN_AFTER:.1%}, omega non-planarity {VAL_OMEGA_DEV_AFTER:.1f} deg;
Ramachandran favoured {VAL_RAMA_BEFORE:.3f} to {VAL_RAMA_AFTER:.3f}; {C3_N_CONVERGED} of 126
converge below 1000 kcal/mol (9KAR ends at {C3_E1_MAX:.0f}); a validity step on 124 of 126, on
2BP4 and 9KAR it breaks a virtual bond. Stage 2 reduces to stage 1 because the best C2 rung is
the shipped prior (L112, addendum 4). The calibration flag, restated per L123: the pool's own disagreement predicts the emitted chain's
error (Spearman {SPREAD_RHO_PARTIAL:+.3f} partial on n and Rg, fold CI {SPREAD_RHO_PARTIAL_CI}); the
relaxation's displacement tracks that disagreement at rho {RHO_MOVED_SPREAD:.2f} and adds nothing
given it (`s26/results/a_strain_vs_spread.json`); never a selector.

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
The line for "which would you do first?", verbatim from the coordinator's ruling (ledger L117,
2026-09-14 02:02; verdicts A REPLACE, B REPLACE, C KEEP WITH EDITS):
"If I could do one thing next, I would publish the trainability work first. Every figure in it
already exists as a measured artefact, it needs no new machine, and it is the one part of this
project whose result is exact and complete. The only open accuracy lever is the distance
prior, and the honest next step there is a larger language model than this laptop can hold, so
that comes second and needs a bigger machine. I would not spend more time on the circuit for
accuracy: we now know why it cannot matter here."
("positive" reads "exact" per the coordinator's L122, accepting the Adversary's L120 caveats.) Why this
order (L117, as amended by L122): the paper's inputs are all in hand; its S26 additions are
Adversary-checked (A2, L45; A4 with L119's bootstrap intervals, L47; the product-state fact, L70),
while its S13 inputs (the locality theorem, `s13/results/qarch_locality_geom.json`; the
Pauli-spectrum prediction, `s13/results/walsh_predict.json`, `geo_pauli.json`) are cited from their
artefacts without an S26 re-check, and the S13 Pauli mean weights 2.236 / 3.015 stay off every slide; the prior lever is real ({PRIOR_SLOPE:.2f} A per unit gamma,
`s24/results/priorladder.json`) but every input this machine can compute is measured flat
(nine rungs, none beats the shipped prior, `s26/results/p_ladder_report_*_s0.json`), so it is a
resourcing decision; the circuit is closed as an accuracy lever by three independent facts (the
product-state optimum, the inert ADAPT growth and the A1 null at {A1_V_X:.2f} x MDE,
`s26/results/a1_stats.json`; the full DLA, `s26/results/q_dla.json`; the set-equality theorem,
`s25/results/q_verify.json`).
