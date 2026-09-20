# PREREG S29 LANE B -- THE COMPATIBILITY HAMILTONIAN: A NON-DEGENERATE OFF-DIAGONAL TERM AND WHAT ITS GROUND STATE SELECTS

Registered 2026-09-20 00:10 Pacific, before any number of this lane. Authority `s29/BRIEF.md`;
operational form `s29/S29_CONTRACT.md` (+ `s27/S28_CONTRACT.md`); brief `s29/briefs/S29B.md`.
Nothing below has been measured by this lane at registration time. What I have already read is
listed in section 0.4 so that "pre-registered" means what it says.

---

## 0. THE RE-OPENING (contract rule 10)

### 0.1 The closure being re-opened, by its ledger line
`s27/LEDGER.md` **S28-L8b** (with its adversary check **S28-L11**, its endpoint **S28-L21** /
**S28-L41**, and **S28-L43**): H = diag(E) - J A with A the Gaussian similarity graph
A_ij = exp(-d_ij^2 / 2 sigma^2) at sigma = median off-diagonal CA-RMSD. Measured: the hop-only
gradient variance decays at **-1.7 to -1.8 per qubit** (S28-L11 caveat (a): quote the range, not
-1.844) and at n = 9 sits **7,300x** below the diagonal terms' 3.051e-2 at J = 1 (820x at J = 3);
97% of that variance is the rank-one part v1 v1^T; the Perron vector is 95 to 99% the uniform
state; lambda_2/lambda_1 = 0.113 to 0.139. The endpoint did not move (S28-L41: every deployed-
readout arm within 0.46x MDE of production; R2/R3 +0.25 to +0.37 A worse). The follow-up on a
spread-spectrum kNN graph (**S28-L25**, **S28-L29**, **S28-L46**, **S28-L47**) halved the decay
rate to -1.0 per qubit and raised the n = 9 variance 30 to 60x, and the endpoint still did not
move; B2 closed.

### 0.2 Why the old closure may not apply
S28's own mechanism sentence is a statement about **one similarity measure**, not about
off-diagonality: "the decay is a property of the near-rank-one SPECTRUM ... not of being
off-diagonal" (S28-L8b), and S28-L11(3) strengthened it (a diagonal control with the same
spectrum decays at the same rate on the padding-free registers). The kNN follow-up spread the
spectrum but kept the same *object* -- a non-negative similarity graph whose top eigenvector is
still 0.98 overlapped with the uniform state (S28-L25), i.e. still a typicality projector with a
flatter tail. Neither matrix ever had **structured** eigenvectors.

### 0.3 The licence (lane T's measurement)
On the same 12 trainability targets and the same sub-pools, lane T measured three matrices side
by side (`s29/s29_T_spectra.py`, job `s26/jobs_done/s29T_spectra.json`, rows
`s29/results/s29_T_spectra_rows.jsonl` (72 rows = 12 targets x n 4..9) and
`s29/results/s29_T_grad_rows.jsonl` (216 rows = 12 x 3 matrices x 6 registers), analysis
`s29/results/s29_T_spectra.json`): at n = 9 the raw Gaussian similarity A has
lambda_2/lambda_1 = **0.138**, its double centering A_c = H A H gives **0.465**, and the signed
agreement matrix G = Delta Delta^T gives **0.634** (the three figures quoted in `s29/briefs/S29B.md`
from lane T's rows). Lane T's ledger entry for this measurement was **not yet posted** when this
prereg was written; it is cited here by artefact path and the entry number will be added to this
lane's first ledger entry when it lands. **The S28 closure was about one degenerate similarity
measure.** This lane asks the question S28 could not: with an off-diagonal term whose eigenvectors
carry structure, does the state change, and does what it emits change?

### 0.4 What I have read before registering (honesty about "pre-registered")
Read: the briefs, contracts and charter; S28-L8b / L11 / L21 / L25 / L29 / L41 / L43 / L46 / L47
in full; `s27/s28_B_hop.py`, `s27/s28_B2_knn.py`, `s29/s29_T_spectra.py` in full; the *schema* and
target list of `s29/results/s29_T_grad_rows.jsonl` (216 rows, 12 targets, keys only). **I have not
read a single `var_g0` value of lane T's gradient rows**, and measurement 1's bright line below is
registered blind to them. I have read the n = 4 spectral row of `1A13` and the n = 4 / n = 5
medians in `s29_T_spectra.json` (they scrolled past while I was reading the script) and the three
n = 9 ratios quoted in my brief (0.138 / 0.465 / 0.634). Measurement 1 therefore registers a
bright line on quantities that are **already on disk**, computed by another lane against *its own*
registered prediction; this is what the brief instructs ("extend lane T's rows, do not redo them")
and it is stated here rather than dressed up. Nothing in measurement 2 or 3 exists on disk.

### 0.5 A correction to the brief's gate constant, registered before it bites
The brief's measurement-2 gate reads "not better than the uniform top-75 average (2.954 point
cloud)". **2.954 is the S10-5 ORACLE ladder figure** (`docs/FINDINGS.md` S10-5, quoted in
`s27/PREREG_S28_A.md` line 9 and `s27/LEDGER.md` S28-L0), from a different pool era. On *this*
instrument the DIS top-75 uniform average on the point cloud is **3.048338** at n = 126
(`s27/results/s28_B_rows.jsonl :: rmsd_dis75`, and `s27/results/s28_B_summary.json ::
anchors.dis75_mean`) and **3.252928** on the 12 trainability targets, which are harder than
average. Gating a 12-target probe against 2.954 would compare an operator against a number from
another target set and another pool -- the project's most repeated error (memory:
`control-must-match-the-operators-space`). **The gate below is PAIRED against each target's own
DIS top-75 uniform average** and the 2.954 figure is reported beside it with this provenance.

---

## 1. THE NINE QUESTIONS (charter section 11; contract rule 15)

1. **What a basis state means.** |i> is candidate i of the target's DIS top-500 retrieval pool --
   one real 9-to-16-residue CA window, posed in the pool's frame. The encoding is the deployed
   identity label (`s22.qcand_lib.Encoding`, candidate i at bit-string i), n = 9 qubits,
   dim = 512 = 500 real candidates + 12 padding states whose E is set to max + 10 sd and whose
   rows and columns in every M are zero. Unchanged from S28 deliberately: the point of this lane
   is to change **one** thing, the off-diagonal term.
2. **What the Hamiltonian means.** H = diag(E) - J M. E is the deployed standardised DIS rank
   ladder (`s27.run_pool.zr(ch["DIS"])` through `Encoding`), the shipped cost. M is a
   **compatibility** operator on the pool with the common mode removed:
   - **A_c = H A H / ||.||_2**, H = I - 11^T/N, A the S28 Gaussian graph at unit spectral norm.
     Double centering is classical MDS: A_c's eigenvectors are the pool's **principal coordinates**,
     i.e. its principal modes of *disagreement*, and the typicality mode is gone by construction.
   - **G = Delta Delta^T / n_res**, at unit spectral norm, with Delta_i the deviation of member i
     from the pool mean in the medoid frame (the frame the deployed readout averages in;
     `s12.instrument.superpose_batch`). G_ij is the inner product of two members' deviations from
     the typical structure: positive when they depart the same way.
   - **A** itself is carried as the S28 reference, not as an arm.
   - **G_res** (the brief's optional third: Delta residualised on the posterior's own predicted
     deviation) is **registered as NOT RUN in this pass**, with the reason in section 6.
3. **Why it should correlate with useful structural information.** The identity
   `<psi|G|psi> = || sum_i psi_i Delta_i ||^2 / n_res` (section 3.1) makes the hopping term
   *exactly* the squared norm of the amplitude-weighted departure from the pool's typical
   structure. The record says 68% of the pool's error is common-mode by exact identity (S23 L9,
   memory `pool-error-is-68-percent-common-mode`), so the pool mean carries a bias that no
   averaging removes; anything that can move the answer at all must move it **along the deviation
   modes**, and G is the operator whose eigenvectors are those modes. A_c is the same statement
   made on the similarity measure instead of the coordinates. This is the charter's "Hamiltonian
   encoding the disagreement between the prior and the pool" (section 12) and it enters the
   common-mode error as the thing that is **removed** rather than averaged in.
4. **What the CVaR objective optimises, over what distribution.**
   F(theta) = CVaR_0.18(E; p_theta) - 0.5 H(p_theta) - J <psi_theta|M|psi_theta>, the deployed
   spine (`core.quantum.free_energy`) plus the hopping term with its exact parameter-shift
   gradient, over p_theta = |psi_theta|^2 on the 512-state register. Unchanged from S28 except M.
5. **What the ansatz can and cannot represent.** `core.quantum.StatevectorCircuit(9, 3)`, RY+CNOT,
   **27 parameters**, real amplitudes only: it provably cannot represent a complex state, and
   S28-L43 measured that it reaches the J = 3 Gaussian ground state at squared overlap **0.849**
   (best of 16 starts, 0.796 to 0.880) -- the cap is expressivity of the state; the 2.8 gap in F
   is optimisation by basin selection.
6. **Whether the optimiser can reach the relevant states.** Measured in S28 and re-measured here:
   the circuit lands in the sign-coherent basin on 43/32 of 126 cells (Gaussian, S28-L41) and
   49/33 (kNN, S28-L46). The *bimodality* is the honest description; no trainability word is
   attached to any slope in this lane (contract rule 9).
7. **Whether the quantum output contains information unavailable to the classical baseline.**
   The exact eigensolver ground state of the same H is the classical counterpart and is an arm,
   not a footnote. If the eigensolver reproduces every effect, the answer is "no" and the entry
   will say so in that word.
8. **Whether a classical control can reproduce the effect.** Controls in section 5.4 / 5.5:
   eigensolver; M rank-permuted; M replaced by a random matrix with the **same spectrum**
   (S28-L43's control, which separates spectrum from structure); a classical top-m by v1^2;
   untrained circuit; product-state restriction; matched budget; both seeds.
9. **Whether it changes built-chain RMSD.** Only measurement 3 answers this, only if the
   measurement-2 gate opens, and only through lane D's meter first (measurement 4).

---

## 2. FINDINGS ENGAGED (contract rule 13)

**Attacks** charter finding **7** (was the failure "non-diagonal Hamiltonians don't help" or "that
particular similarity measure was degenerate"? -- this lane answers it), finding **11** (the 68%
common-mode error: G is the operator that removes it rather than averaging it in), and finding
**9** in its weak form (Hamiltonian variation within a fixed encoding is low-dimensional -- true,
and this lane tests exactly how low).

**Accepts as binding** findings **1** (expressivity is not the bottleneck), **2**, **3** (the
objective is the bottleneck), **4**, **5**, **6** and **10**.

**Finding 8 (recognition), engaged explicitly as required.** This lane does **not** claim to
repair recognition and its Hamiltonian contains no new information channel: E is the shipped cost
and M is built from the pool's own geometry, which the pipeline already has. What it proposes is
that the *aggregation* stage, not the recognition stage, is where the 68% common-mode error is
paid, and that a correlated state over mutually compatible members is the operator class that
could exploit it. My own registered prior (section 3.3) is that it **cannot**, for a reason that
is provable rather than empirical, and measurement 2 is designed to say so in one afternoon rather
than after an endpoint run. Under S29-L1 (lane L: no native-free QA method exists at this length)
this is the honest framing: an aggregation experiment, not a recognition claim.

---

## 3. THE MATHEMATICS, AND MY PRIOR STATED AS A DERIVATION

### 3.1 The hopping term of G is the squared coherent departure
With Delta the (k x 3 n_res) matrix of members' deviations from the pool mean in the medoid frame,
G = Delta Delta^T / n_res, so for any real amplitude vector psi

    <psi|G|psi> = psi^T Delta Delta^T psi / n_res = || Delta^T psi ||^2 / n_res
                = || sum_i psi_i Delta_i ||^2 / n_res.

Maximising it (J > 0) asks for a set of members whose departures from typical **reinforce**;
minimising it (J < 0) asks for a set whose departures **cancel**, which is what uniform averaging
does by construction. Both signs are on the grid, and the J < 0 branch is the record's own
consistency mechanism (S27 section 6: consistency mechanisms have been harmful here; memory
`consensus-is-outlier-avoidance`) measured in the same experiment as its opposite.

### 3.2 The sign-mixing lemma (why I expect the term cannot help through a sign-blind readout)
**Claim.** Both A_c and G annihilate the uniform vector: A_c 1 = H A H 1 = 0 because H1 = 0; and
G 1 = Delta (Delta^T 1) / n_res = 0 because deviations from the mean sum to zero
(sum_i Delta_i = 0). Therefore every eigenvector v of either matrix with lambda != 0 satisfies
1^T v = 0, and a nonzero vector orthogonal to 1 **has entries of both signs**.

**Consequence.** As J -> +-infinity the ground state of diag(E) - J M tends to the extremal
eigenvector of M, which is a **signed contrast**: amplitude on *both* extremes of a principal
disagreement mode with opposite signs. The deployed readouts consume p = |psi|^2 and average the
selected members **uniformly and sign-blindly**, so the support is symmetric about the mode and
its coordinate average returns approximately the pool mean -- the very structure the term was
supposed to depart from. The coherence lives in the **phases**, and the readout is where it is
destroyed. The only thing that can tilt the support to one side is diag(E), and E is the shipped
cost, which is anti-correlated with RMSD along the near-native half of the ladder (S29-L2:
rho -0.402 chain / -0.182 CA).

This lemma is the reason measurement 2 exists and is registered **before** it: it converts "does
the term help?" into a cheap exact eigenproblem, and it also names the repair (a **signed**
readout, R4 below), which is measured as a diagnostic in the same pass.

### 3.3 The stable-rank bound on measurement 1's second clause, derived before measuring
Lane T's section-3 derivation gives, as the 2-design leading term,
Var_theta[dF/dtheta_k] ~= ||M - (tr M / D) I||_F^2 / D^2 = r_stable(M) / D^2 at unit spectral norm,
r_stable = ||M||_F^2 / ||M||_2^2 <= rank(M). At n = 9, D = 512 and the diagonal terms' measured
variance is 3.051e-2 (S28-L8b). "Within 30x of the diagonal terms'" therefore requires

    r_stable(M) >= 512^2 * 3.051e-2 / 30 = **266**.

**G cannot satisfy this by construction**: rank(G) <= 3 n_res and, after the mean and the rigid
body, <= 3 n_res - 6, i.e. **at most 42** for the longest peptide here (n_res 9 to 16). So
r_stable(G) <= 42 and the predicted n = 9 variance is at most 42/262144 = 1.6e-4, at least
**190x** below the diagonal terms'. For A_c the bound is rank <= 511, so the clause is not
impossible a priori, but it needs an effectively half-rank flat spectrum, which lambda_2/lambda_1
= 0.465 with a decaying tail does not describe. **I predict the second clause of the bright line
fails for both matrices**, and that the first clause (slopes shallower than -1.0 per qubit) is
straddled or passed. The measurement then checks the 2-design relation itself (the circuit at
depth 3, n = 9 has 27 parameters and is nowhere near a 2-design, S25 `q_plateau.py` scope note),
which lane T's rows already carry as `ratio_meas_over_pred`.

### 3.4 Should a coherent set along a principal disagreement mode be better or worse than the
average? (the brief's explicit question; S27 section 6)
**Worse, in expectation, and for a stated reason.** A one-sided coherent set emits
mean + c * (a principal deviation mode). The native's own deviation from the pool mean is, on any
given target, at some angle to that mode; the record's measurements of that class of operator are
that the per-target **sign** is the unavailable quantity (memory `in-band-ordering-is-per-target`:
0.986 within a target, 0.600 across; `error-shape-not-mae-decides-ranking`: full amplitude, wrong
direction; S16: every native-free arm chose "do nothing"). A zero-mean signed displacement against
a locally convex RMSD costs more when it is wrong than it gains when it is right, so the expected
effect of a sign-blind operator along a disagreement mode is **worse than the average**, not
neutral. Against that: S10-5's ORACLE ladder puts the affine span of the same 500 windows at
0.064 A and the convex hull at 0.853 A, so the *headroom* along these modes is enormous; the
question is entirely whether anything native-free can point along them. This lane's answer will be
a measurement, and its prior is failure.

---

## 4. MEASUREMENT 1 -- SPECTRUM AND TRAINABILITY (no RMSD, no native, its own ledger entry)

**Question.** Is S28's gradient decay a property of the *degenerate spectrum*, and does a
non-degenerate off-diagonal term restore the hopping term's share of the gradient?

**Inputs.** Lane T's rows, consumed, not recomputed (brief: "extend lane T's rows, do not redo
them"): `s29/results/s29_T_spectra_rows.jsonl` (spectra) and `s29/results/s29_T_grad_rows.jsonl`
(hop-only gradient variance, mode `hop_only`, 120 draws, seed 1009, init sd 0.6, depth 3, exact
parameter shift -- S28's estimator line for line via `s27.s28_B_hop.measure_hop`), for
M in {A, A_c, G} at n = 4..9 on S27's 12 trainability targets. Beside them: S28's diagonal cells
from `s27/results/s28_B_train.json :: summary` (full J = 0 at n = 9: 3.051e-2; hop-only A:
4.157e-6) and the kNN rates from `s27/results/s28_B2_train.json`.

**Reported.** Per matrix and per n (median over the 12 targets): lambda_2/lambda_1,
lambda_3/lambda_1, the stable rank and its traceless form, the top eigenvector's participation
ratio and its squared overlap with the uniform state, the hop-only Var[dF/dtheta_0], the ratio to
the diagonal terms', the fitted log2 slope over n = 4..9 **and** over the padding-free n = 4..8
(S28-L11 caveat (a): quote a range, never a third decimal), and the measured/predicted variance
ratio that tests section 3.3's relation.

**B1, the pre-registered bright line** (the brief's wording, disambiguated). The claim
"**the degeneracy was the mechanism of S28's gradient invisibility**" is
- **SUPPORTED** iff, for **both** A_c and G: (i) the fitted log2 slope is **shallower than -1.0
  per qubit** on the n = 4..9 fit, **and** (ii) the median n = 9 hop-only variance is **within 30x**
  of the diagonal terms' 3.051e-2, i.e. **>= 1.017e-3**;
- **REFUTED** (the brief's falsifier) if either clause fails for both matrices;
- **PARTIAL** otherwise, and the entry says which clause failed for which matrix and quotes the
  numbers without a verdict word.
Registered prior (section 3.3): **clause (ii) fails for both; G's failure is a bound, not a
measurement**. Contract rule 9: no slope in this entry is called a plateau or its absence, and no
quantum-advantage word is used.

**Extension this lane adds** (not a redo): the measured-over-predicted variance ratio per matrix
and per n; the stable-rank-versus-variance scatter across all 216 rows; and the *count* of rows
in which G's rank exceeds 3 n_res - 6 (must be zero; it is a correctness check on the builder).

**Multiplicity.** Zero endpoint comparisons. This entry posts whether or not measurement 2 clears.

---

## 5. MEASUREMENT 2 -- WHAT THE EXACT GROUND STATE SELECTS (ORACLE DIAGNOSTIC, LABELLED)

Every RMSD in this section is **ORACLE** (it reads the native to score an achievable selection)
and **nothing in it chooses a deployable parameter**. The selection itself is native-free
(`select` receives W, E, the tie key and M; the native enters only in `oracle_rmsd`).

### 5.1 The object
For each target, each M in {A_c, G, A} and each J on the grid below: the exact ground state of
H = diag(E) - J M by `numpy.linalg.eigh` (512 x 512), p = v0^2, with the S28 global-sign fix.
J = 0 is the one-hot argmin and is **excluded from every comparison** (degenerate; S28's rule).

### 5.2 The J grid, fixed here and not revisited
**J in {-3, -1, -0.3, -0.1, +0.1, +0.3, +1, +3}**, both signs (section 3.1: the negative branch is
the cancelling/consensus direction, the positive branch the coherent one), magnitudes matched to
S28's grid so the A rows are directly comparable to S28-L21. 8 J x 3 M = **24 cells per target**.

### 5.3 Readouts and descriptors, per cell
- **R3 (primary, production-matched)**: the uniform coordinate average over the **m = 75** most
  probable real candidates, ties by the stable key (`s27.run_pool.topm`, contract rule 12). This
  is the brief's "coordinate average over its top-m support" at production's own m.
- **R2 (secondary)**: the p-weighted coordinate average over all 500 real candidates in the
  p-weighted consensus-medoid frame (`s27.s28_B_hop.readout_weighted`).
- **R4 (diagnostic, the mechanism-matched readout)**: the **signed** amplitude-weighted average
  sum_i psi_i W_i / sum_i psi_i, in the medoid frame. This is the readout for which
  <psi|G|psi> is the meaningful quantity (section 3.2). It is **undefined** when the state is
  sign-antisymmetric, so it is recorded only where the sign coherence
  (sum psi)^2/(sum |psi|)^2 >= 0.01, and the number of undefined cells is reported. Because an
  affine combination does not preserve scale, R4 is reported **with** its emitted mean virtual
  CA-CA bond and radius of gyration beside it (memory `averaging-space-beats-the-objective`:
  averaging contracts the backbone 25.8%; contract addendum 1 rule 20's spirit).
- **Descriptors (native-free)**: participation ratio 1/sum p^2; mean pairwise CA-RMSD of the R3
  support (the cluster it selects); Jaccard and overlap of that support with the DIS top-75;
  sign coherence; the spectral gap of H; the realised one-sidedness of the support along M's top
  eigenvector.

### 5.4 Controls in this measurement (the diagnostic set; the full set is section 6)
- **PROD**: the same target's DIS top-75 uniform average -- the paired comparator.
- **PERM**: M rank-permuted (P M P^T, same spectrum and degree multiset, correspondence with E
  destroyed; `s27.s28_B_hop.permuted_graph` generalised to any symmetric M, seed
  `run_pool.rng_for(pdb, "s29B_perm")`).
- **SPEC**: M replaced by Q Lambda Q^T with Q a Haar orthogonal matrix and Lambda M's own
  eigenvalues -- **the same spectrum, no structure** (S28-L43's control; seed
  `rng_for(pdb, "s29B_spec")`). This is the control that separates "spread spectrum" from "these
  particular eigenvectors", which is exactly the claim under test.

### 5.5 THE GATE (pre-registered; the brief's stop rule, made paired and priced)
Run on the **12 trainability targets** first (contract rule 16). Let d(cell) = mean over the 12
targets of [R3 ORACLE RMSD of the cell] - [the same target's DIS top-75 uniform average], paired,
through `s24.stats_lib.compare` with `pinned_folds`.

**GO to measurement 3 iff all four hold:**
1. some cell with M in {A_c, G} has **d <= -0.7 x MDE** of its own paired comparison;
2. that cell's advantage survives the grid as an order statistic: `ST.best_of_k_within` over the
   24-cell matrix reports a **positive split-half transfer** (memory `grid-oracles-are-order-statistics`);
3. the same cell's **PERM and SPEC** controls do not reproduce it (control effect < half the real
   effect), so the claim is about *these eigenvectors*, not about having a spread spectrum;
4. the cell is not the J -> 0 limit (its participation ratio is > 1.5, i.e. it is not the argmin).

**Otherwise: STOP.** The entry says "the term cannot help at the endpoint and here is the number",
the endpoint run is not spent, and measurement 3 is not run. A GO is **a decision to spend
compute, not a result** (S29-L4(d)): the entry prints power and the Type-M factor beside it and
says so in the same sentence; with 12 targets the fold-clustered CI over 5 clusters of 2 to 3
targets is quoted as **descriptive only**.

If the gate opens, measurement 2 is repeated on all **126** targets before any endpoint run, and
the 126-target version is the one quoted.

**Registered prior for measurement 2: the gate does NOT open** (section 3.2 + 3.4). Specifically:
every R3 cell lands within 0.3 A of PROD with the sign mostly positive (worse); R4 has a *larger*
spread than R3 in both directions and a worse mean; and the ORACLE best-of-sign over R4 is much
better than PROD while R4 itself is not, which is the per-target-sign problem restated.

**Multiplicity.** 24 cells x 3 readouts x (REAL, PERM, SPEC) = up to 216 ORACLE diagnostic
comparisons on 12 targets. None is an endpoint comparison; all are counted in the entry and the
grid is priced with `best_of_k_within` before any cell is quoted.

---

## 6. MEASUREMENT 3 -- THE ENDPOINT (ONLY IF THE GATE OPENS)

Registered now so that a GO cannot be followed by a freshly-invented design.

**Arms.** The genuine CVaR-VQE (`s27.s28_B_hop.run_hop_vqe`, i.e. `core.quantum.run_cvar_vqe`
with the hopping gradient added: alpha 0.18, T 0.5, depth 3, 80 iterations, lr 0.15, Adam), at the
**one or two** (M, J) cells the gate opened, seeds 0 and 1, readouts **R1** (the deployed tail)
and **R2**, on **12 targets first, then 126** only if the 12-target result clears 0.7x MDE.

**Comparators.** Production (DIS top-75 uniform) and the same arm at J = 0, same seed, same
readout (F1's shape in S28).

**Controls, all pre-registered (contract rule 15's ten):** J = 0; M rank-permuted; M replaced by a
random matrix with the same spectrum (SPEC); the exact eigensolver ground state (the classical
spectral counterpart); the untrained circuit (best of 16 draws, the same draw law as S28-L43's
representability fit); a **product-state restriction** (the same objective optimised over
separable states, i.e. the depth-0 / no-CNOT circuit at matched parameter count); **matched
budget** (the J = 0 arm given the same number of objective evaluations); both seeds; and the
order-statistic pricing of any grid. A positive that does not survive all of these is reported as
provisional and handed to lane D.

**Falsifiers.** F-B1: the arm's built-chain mean is not better than production by 0.7x MDE with
the fold-clustered CI excluding zero on **both** seeds -> refuted as registered. F-B2: the PERM or
SPEC control reproduces at least half the effect -> the specific eigenvectors contributed nothing.
F-B3: the eigensolver reproduces the effect -> the quantum stage contributed nothing and the entry
says so in that word. **Registered prior: WORSE or null**, as in S28.

**Reporting basis.** Built chain (`s24.d_harness.readout_projected`), point cloud as the screening
basis only, both re-projected in this lane's own process against the S27 anchor
(`s27/results/s28_B_prodcheck.json`'s 0.0-on-126/126 check is the standard to meet).

**G_res, and why it is not run in this pass.** The brief's optional third matrix residualises each
member's deviation on the posterior's own predicted deviation. The posterior here is a distogram,
not a coordinate predictor, so "the posterior's predicted deviation" requires first materialising a
predicted structure and then defining a deviation in the medoid frame -- a second construction with
its own convenience choices, and by lane T's theorem 2 / contract addendum 1 rule 21 the marginal
class it would live in is bounded. It is deferred until G itself has been measured, and if the
gate does not open it is not built at all. Recorded here so its absence is a registered choice.

---

## 7. MEASUREMENT 4 -- THE METER (contract rule 19)

Before any endpoint claim, every arm's cost goes through lane D's cost-RMSD meter
(`s29/s29_D_cost_audit.py`). The object metered is the **structure-level** cost this lane's
Hamiltonian induces, submitted as a callable `s29.s29_B_compat:cost_compat` under the meter's
`f(W, ctx) -> (m,)` contract: for a stack of CA clouds it returns the shipped DIS cost of each
cloud **minus J times the cloud's coherence with the pool's principal disagreement mode**, which
is the term's contribution to how a structure is scored. The four numbers (ladder Spearman on the
S28 ladder, gradient cosine at production with the random-direction reference, the native's
percentile in the pool, the preference for the ORACLE structure with the pool-member control) are
reported in the endpoint entry. Per contract addendum 1 rule 20, any cosine gain is reported with
(a) the implied shrink, (b) the native percentile, (c) the emitted mean virtual bond and Rg.

---

## 8. BUILD, JOBS AND DISCIPLINE

- Code `s29/s29_B_compat.py`, reusing `s27/s28_B_hop.py`'s machinery (pairwise matrix, graph,
  padding, readouts, gate, ORACLE scoring) rather than re-deriving it. Tests
  `tests/test_s29_B.py`: symmetry; zero diagonal where required; unit spectral norm for every
  matrix; the centering identity A_c 1 = 0 and G 1 = 0 to 1e-10; rank(G) <= 3 n_res - 6;
  G's nonzero spectrum equals that of Delta^T Delta / n_res (the pool's deviation covariance);
  PERM and SPEC have M's spectrum to 1e-10; the ground state at J = 0 is the one-hot argmin; and
  **NaN-poison** of `nat_ca` / `oracle_rr` leaves every native-free output bit-identical.
- Results `s29/results/s29_B_*.json(l)`; findings `s29/s29_B_FINDINGS.md`; ledger entries
  `## S29-L<n> -- ... (date time, B)` numbered from the tail in ONE python process, with the date
  read by `date` in the same command as the append (memory `check-the-clock-before-stamping`).
- Every job over a minute or 200 MB through `python s26/jobrun.py --agent S29B --tag CPU --name
  <name> --est-ram <GB> --`, one target probed first with its peak RSS quoted, **per-target
  checkpoints** in an append-only rows file, and the shards sized so the box fills toward 94% RAM
  without breaching it. Jobs are waited on inside a background `until` loop, not by handing back.
- Statistics `s24.stats_lib.compare`; the fold-clustered CI decides; below 0.7x MDE is not a
  result; 0.7 to 1.3x is the Type-M zone and is labelled; `ST.fmt` verbatim wherever a contrast is
  claimed; SE beside every mean (memory `mde-is-per-comparison-not-per-instrument`); ties never by
  array order; every number carries its artefact path.
- Lane D attacks every positive and will demand the **three-way split** -- Hamiltonian quality
  (what the exact ground state selects), optimisation quality (how far the circuit gets toward it),
  emitted structure (what the readout makes of it) -- on any arm that moves. The three legs are
  measured in measurements 2, 3 and 3 respectively and are reported together, not on request.

---

## 9. WHAT WOULD MAKE ME SAY THE WHOLE DIRECTION IS DEAD

If measurement 1 shows both A_c and G still 100x or more below the diagonal terms at n = 9 **and**
measurement 2's gate does not open at any J of either sign, then: the off-diagonal route on the
**candidate-index register** is closed, not because the similarity was degenerate (S28's reason,
now measured to be one measure's property) but because the readout is sign-blind while every
common-mode-removed operator's eigenvectors are sign-mixed (section 3.2) -- a *structural*
statement about the encoding, which is the charter's finding 9 made precise and which points at
lane X's configuration space rather than at another matrix. That is the closure this lane would
post, and it is worth more than another Hamiltonian.

---

# ADDENDUM 1 (2026-09-20 00:09 Pacific, registered before any number of this lane)

The coordinator has directed me to lane T's **S29-L11** (theory section 3, posted 00:06). I had
not read it when sections 0 to 9 above were written and committed (`1c345f07`); the entry
pre-dates my commit on the clock, so I claim independence of **reasoning**, not of clock, and I
note the convergence rather than the priority: S29-L11(d)'s pole symmetry is section 3.2's
sign-mixing lemma, and S29-L11(b)'s `r_stable <= 3 N_res - 6 <= 42` cap is section 3.3's bound.
Three things change, all registered here before measuring.

## A1. Measurement 1 becomes an independent-implementation check of a derived law
S29-L11 supplies the law `Var[dF/dtheta_0] ~= r_stable(M)/D^2` with three no-parameter checks
(S28's Gaussian to 7%, S28-B2's kNN 46x, lane D's J* = 85.7 as 88) and, from lane T's own job,
the answer to my bright line: slopes **-2.305 (A_c)** and **-1.900 (G)** against A's **-1.830**,
n = 9 variances 3.599e-6 / 5.051e-6 / 4.447e-6. **My B1 is therefore already answered NO**, as
section 3.3 predicted it would be, and a build justified by "the spectrum is no longer degenerate"
is dead. B1 stands exactly as registered and will be reported as REFUTED with the numbers; nothing
about it is rewritten after the fact.

What measurement 1 now is: an **independent implementation** of the same quantities in
`s29/s29_B_compat.py` (my own A_c and G builders, my own hop-only gradient-variance estimator with
a finite-difference cross-check of the parameter-shift rule), compared row for row against lane
T's 216 rows, plus the law check. A derived law with three independent checks deserves a fourth
from code that does not share lane T's.
Adopted **verbatim** as additional registered falsifiers (S29-L11 predictions 1 and 2):
- **T1**: on A_c and G the hop-only slope is -2.3 +- 0.3 and -1.9 +- 0.3, not flat, and both n = 9
  variances are within 1.5x of the Gaussian's -- **falsified if either slope is above -1.3 or
  either n = 9 variance exceeds 3e-5**.
- **T2**: for any unit-spectral-norm observable, `Var[dF/dtheta_0] = r_stable/D^2` within 2x at
  n >= 7, and `J* = D sqrt(0.0305/r_stable)` (for G at n = 9: Var 6.4e-6, J* 71) -- **falsified by
  any observable at n >= 7 departing by over 3x**. My independent rows decide both.
- My own prediction, registered: my A_c and G builders will reproduce lane T's `var_g0` to within
  the draw noise of a shared seed law (they use the same 120 draws at seed 1009), and my A rows
  will reproduce `s27/results/s28_B_train.json :: hop_only|J1` at n = 4..8 exactly, as lane T's
  did. A disagreement is a bug in one of the two implementations and will be chased, not averaged.

## A2. Measurement 2's primary readout is SIGNED, and the pole symmetry is the thing under test
S29-L11(d) is right and section 3.2 says the same: every deployed readout is a function of
p = psi^2, the ground state of a centered M is a signed contrast, so a p-readout cancels it to
first order. Measurement 2 is therefore amended:
- **PRIMARY readout: R4, the signed amplitude readout**, reusing S28 lane A's operator unchanged
  (`s27/s28_A_amp.py :: Frame`, `readout`, `weight_diag`, `struct_diag`): affine weights
  w = psi/sum(psi) over the real candidates, in the DIS-top-75 medoid frame, undefined and
  recorded as such when |sum psi| < 1e-9. Lane A's arms were refuted for **accuracy under the
  shipped objective** (S28-L27b), not as a readout, and this is the first Hamiltonian for which
  the signed readout is the mechanism-matched one. Its emitted `rg` and `bond` are reported beside
  every RMSD (an affine combination does not preserve scale).
- **CO-PRIMARY: R3 (p-top-75)**, kept in order to test **S29-L11 prediction 3** directly: with a
  p-readout, the J -> large ground state of diag(E) - J A_c emits the pool mean to within the
  projection floor on **>= 80%** of targets. Registered falsifier for prediction 3: fewer than 80%
  of targets within the floor at the largest |J| on the grid. (I take "the projection floor" to be
  0.05 A of CA-RMSD between the emitted cloud and the top-75 uniform average, and I state that
  choice here rather than after seeing the numbers.)
- **The one-parameter collapse, tested directly.** S29-L11(d) predicts the family emits
  `production +- eta PC1(pool)`. I will regress each cell's R4 emitted cloud on the two-term model
  `c_prod + eta * PC1(pool)` (PC1 = the pool's first shape mode in the same frame, from the
  deviation covariance whose Gram is G) and report the fitted eta, the R^2 of that fit, and the
  residual norm. **Registered prediction: R^2 >= 0.9 on a majority of cells.** If it is low, the
  family is richer than one parameter and that is the finding.

## A3. The sign is the whole question, and it is measured as such
Under the signed readout, everything rests on the sign of eta, which section 2 of lane T's theory
says the marginals do not supply. Registered, before measuring:
- Report the **ORACLE best-of-sign** RMSD (min over +-eta) beside the **realised** sign the state
  itself produces. The gap between them is the price of the sign.
- Test whether the state's own sign carries ORACLE skill at all: the sign of
  `sum_i psi_i <Delta_i, PC1>` against the sign that reduces RMSD, as a binomial test against 0.5
  over targets, with the tie rule of memory `tie-breaking-leaks-the-pool-order` (average the
  outcome over a tied argmin set; never `np.argmin` on a tied signal).
- **A cell that beats production only under the ORACLE-chosen sign is NOT a GO.** It is the
  per-target sign problem restated (memory `in-band-ordering-is-per-target`), it is reported in
  those words, and no endpoint run is spent on it. This is added to the gate of section 5.5 as
  **condition 5**.

## A4. The 126-target run waits for lane O's ceiling
S29-L11 prediction 4 asks lane O for `min_eta mean ORACLE RMSD(c + eta PC1(pool))` at the best
global eta chosen leave-fold-out; predicted **under 0.15 A** better than production. That number
is the ORACLE ceiling of this entire lane's family under the signed readout. Registered reading
rule, before it arrives:
- **above 0.30 A**: lane T's prediction fails, the family is worth much more than the theory
  thinks, and I say so loudly in the entry;
- **under 0.15 A**: my measurement-2 gate is effectively pre-decided, and I say **that** in the
  entry, in the same sentence as any cell that clears, before spending an endpoint run;
- between: read on the measurement.
**Nothing of this lane launches at 126 until that number is on the board.** The 12-target
measurement 2 runs now, because it is minutes and because it decides whether 126 is worth asking.

## A5. What has NOT changed
The gate of section 5.5 (now with condition 5), the J grid, the control set (PERM, SPEC,
eigensolver, untrained, product-state, matched budget, both seeds), the registered priors
(measurement 1's B1 refuted; measurement 2's gate does not open), the ORACLE labelling, the
multiplicity accounting, and section 9's closure statement.

---

# ADDENDUM 2 (2026-09-20 00:33 Pacific, registered before any number of the build it describes)

The coordinator has directed this lane to a second build, derived by lane T in **S29-L15 (Q1)** and
**S29-L17 (section 4b)**: "TAIL THEN AGGREGATE". Both entries read in full before this addendum.
The compatibility-Hamiltonian measurements 1 and 2 (sections 4 and 5 above) are already running and
are finished and posted as registered -- measurement 2 is the ORACLE diagnostic that decides
whether a *centered* off-diagonal term can help at all, and it is worth its two minutes whatever
happens to the new build. Measurement 3 (that lane's endpoint) is gated on it as registered and I
do not expect to spend it. What follows is a **new measurement 5**, with its own falsifiers.

## B2.1 The object
    F(theta) = CVaR_alpha(E; p_theta) - T H(p_theta) + lam f(R_alpha(p_theta)),
    R_alpha(p)  = (1/alpha) [ sum_{x in S(p)} p_x W_x + (alpha - mass(S)) W_{x_q} ],
S(p) the strict CVaR tail along the E order, x_q the boundary state. R_alpha is the tail's own
coordinate average **as a function of p**, so the objective sees *which candidates populate the
tail and with what weight* -- the freedom S29-L15 proves the deployed objective is indifferent to.
Envelope gradient (S29-L17): dR/dp_y = (W_y - W_{x_q})/alpha on the strict tail, so dF/dp_y gains
lam <grad f(R), W_y - W_{x_q}>/alpha there, and the parameter-shift chain costs the same 2P circuit
evaluations as the deployed objective. The tail's ORDER stays a per-state scalar (E), so the cells
are indexed by the SET and R is linear in p on each; a convex f makes F convex on each cell.

## B2.2 Why this is not S28 lane A, stated with both ledger lines
Lane A (`s27/PREREG_S28_A.md`, endpoint **S28-L27b**, refuted) put lam f(C) on the **amplitude**
readout C = sum psi_i W_i / sum psi_i -- signed affine weights over the whole pool -- and was
refuted for accuracy. Here f acts on the **CVaR tail's own average**, which lane A never tested,
the weights are non-negative and confined to the tail, and the first-instance claim is not accuracy
but **FLATNESS**: the object S29-L15 shows the deployed spine cannot see.

## B2.3 A derivation that corrects the flatness gate BEFORE it is measured
The coordinator's gate is "recompute M5 (the fraction of readout-relevant directions along which
the objective is exactly constant at production, deployed 437/511 = 85.5%); if it does not fall
materially the idea is dead". **I register, before measuring, that the raw 85.5% will NOT fall, and
that this does not kill the idea**, because the two numbers are different quantities:

For y strictly above the VaR, p_y does not appear in R_alpha and (generically) perturbing it moves
neither S(p) nor q, so **dR_alpha/dp_y = 0 by the same envelope argument that gives
dCVaR/dp_y = 0**. The f term is therefore flat on *exactly the same* 437-direction subspace as the
CVaR term. Clause (i) below is predicted UNCHANGED at 85.5% for every lam.

What changes is the **overlap**. Under the deployed pair (CVaR + entropy, uniform-average-over-the-
tail-SET readout) the emitted structure is a piecewise-constant function of p: along every
continuous simplex direction its derivative is exactly zero, and it moves only when the SET
changes, i.e. through the single discrete scalar m (S29-L15's reduction). Under TTA the emitted
structure is R_alpha, which is differentiable on each cell and whose non-zero directions are
**exactly** the objective's 74 non-flat ones. So the honest M5 has three clauses and I register all
three:

  (i)   flat_obj: the fraction of the D-1 simplex directions on which dF/ddelta is exactly 0 at
        production. Deployed 85.5%. **Predicted unchanged** under TTA at every lam.
  (ii)  flat_readout: the same for the emitted structure. Deployed: 100% of continuous directions
        (the readout is piecewise constant). **Predicted 85.5% under TTA** (the same 437).
  (iii) **THE ONE THAT MATTERS -- the overlap**: among the directions along which the EMITTED
        STRUCTURE moves, the fraction along which the OBJECTIVE is exactly constant. Deployed:
        undefined over continuous directions, and operationally one discrete scalar, which the
        objective does fix. TTA: **predicted 0%** -- every direction the readout consumes, the
        objective sees.

**Registered gate G5 (replacing "does 85.5% fall").** The mechanism works iff (iii) is 0% or near
it while (i) stays at 85.5%, AND the emitted structure's realised sensitivity to the tail weights
is non-zero (measured as the norm of dC/dp along the tail directions, which is identically zero for
the deployed readout). If (iii) does not fall, the idea is dead and I say so in one line. If clause
(i) *does* fall, my derivation above is wrong and I say **that**, loudly, in the same entry.

## B2.4 The readout question, which decides whether this is a real arm at all
For the objective to steer what is emitted, the **emitted structure must be R_alpha**, the
p-weighted tail average -- NOT production's uniform average over the tail set. That is a change of
readout, and the record is against sharper-than-uniform weights (`s24.d_harness.readout_uniform`'s
own docstring: every departure measured -- sharper weights, clustering, smaller m, amplitude
weighting -- attacks the variance reduction that makes averaging work; S28-L21 / L41 measured the
p-weighted-over-the-whole-pool readout at +0.25 to +0.37 A WORSE than production). At the deployed
optimum p rises 8.4x across the prefix (S29-L15), so R_alpha is a materially sharper average.
Two arms are therefore registered and BOTH reported:

  - **TTA-w** (the honest one): objective and readout are both R_alpha.
  - **TTA-u** (the diagnostic): objective on R_alpha, emitted structure the deployed uniform
    average over the same tail set. The objective then steers a quantity it does not emit, so a
    gain here is weaker evidence; it is carried to separate "the readout change hurt" from "the
    steering did not help".

**Registered zero-lam anchor, measured before anything else:** TTA-w at lam = 0 against production.
If the readout change alone is already worse than production by more than the effect any lam could
plausibly buy, that is reported in the first sentence of the entry and the arm is priced against
*its own* lam = 0, never against production alone.

## B2.5 f, and the price it pays
f is the shipped distogram risk evaluated on R_alpha -- a **marginal-class** objective, so by lane
T's theorem 2 (S29-L7) and contract addendum 1 rule 21 it cannot be locally informative about the
native's deviation from typical, and by S29-L2 its ladder correlation on the near-native rungs is
-0.182 (CA) / -0.402 (chain). f goes through lane D's meter
(`--f s29.s29_B_tta:cost_tail_avg`) **before** the endpoint run, with rule 20's three companions
(implied shrink, native percentile, emitted bond and Rg). If lane D's band experiment finds a
scorer with positive in-band skill, that scorer is registered as an alternative f and the
substitution is reported as a separate arm, not a replacement.

## B2.6 Falsifiers and controls for measurement 5
- **F5a (mechanism)**: G5 above. Minutes, no endpoint run, reported whatever it says.
- **F5b (endpoint, 12 targets first)**: lam on the pre-registered grid
  **lam in {0, 0.1, 0.3, 1.0, 3.0}**, seeds 0 and 1, against (a) production, (b) its own lam = 0,
  (c) lane D's **fixed-profile control M6** (the target-independent rank-weight profile
  p*(alpha, T) applied to the target's own DIS order, no circuit and no optimiser -- the control
  that matters, S29-L15), (d) the untrained circuit at best-of-16, (e) matched budget. Refuted as
  registered if no arm beats BOTH production and its own lam = 0 by 0.7x MDE with the fold CI
  excluding zero on both seeds. The grid is priced with `ST.best_of_k_within`.
- **F5c (the classical counterpart rule 15 requires, from S29-L17)**: greedy plus local search over
  75-subsets minimising f(mean of the subset). T predicts it beats the DIS top-75 on the OBJECTIVE
  by at least 0.10 and is WORSE than production on the built chain by at least 0.1 A. If the second
  clause fails -- the aggregate-optimal subset is BETTER than production -- that is the sprint's
  first genuine opening and goes to 126 immediately. This control runs BEFORE the VQE endpoint,
  because it is cheaper and it bounds what the VQE could find.

## B2.7 My registered prediction for measurement 5
The **mechanism works and the endpoint does not move, or moves the wrong way**: (iii) falls to 0%,
the emitted structure becomes genuinely sensitive to the tail's internal weights, and the endpoint
is WORSE -- partly because f is marginal-class (theorem 2) and steers toward the contracted typical
average, and partly because R_alpha is a sharper-than-uniform average and the record prices that at
+0.25 to +0.37 A. This is the same prediction the coordinator registered, reached from the same two
theory entries, and I record it as agreement rather than as independent support. **If the endpoint
improves, I hand it to lane D the same hour with the three-way split and do not write it up as a
result first.**
