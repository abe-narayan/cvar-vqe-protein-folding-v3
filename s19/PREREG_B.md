# PREREG_B — Sprint 19, AGENT B (VQE / CVaR / QUANTUM)

**Written before any arm ran. Not to be edited after seeing results.**
If a specification here turns out to be wrong, it stays, is flagged as mis-specified, and the
untested regime is recorded OPEN.

Salt: `s19qb`. Seeding: `s15/seed.py::stable_rng` only. `hash()` appears nowhere.

---

## 1. The question, stated so it can fail

> **What should a genuine CVaR-VQE *generate* that classical methods cannot efficiently generate
> at matched budget?**

Sprint 16 closed VQE-as-search. Sprint 17 closed VQE-as-sampler **on the k=4 lattice**, on the
(member error, diversity) plane, with ε-dominance CIs excluding zero for all ten configurations.
Sprint 18 closed VQE-as-scalar-minimiser of the degree-1 objective, and showed the headline that
motivated it was an artefact of a 2-bit labelling.

**What is not closed, and is this lane's only live claim:** those tests all lived on a *discrete
lattice register* with an *arbitrary state labelling*, and Sprint 18's own headline is that
conclusions carried by such a labelling do not port. This lane therefore re-asks the sampler
question in **continuous torsion space**, where the pipeline actually lives, with a latent whose
only arbitrariness is a Z2 per residue (exhaustively testable), and against the full classical
control set at matched budget.

**This is a hypothesis, not a conclusion. A clean closure is the expected outcome and an
acceptable one.**

---

## 2. The prior that damages the hypothesis, recorded BEFORE the arms — as Sprint 18 did

Two facts are already on the record and are stated first because they predict failure.

**(P-a) The CVaR optimum is indifferent to the population.** For an unrestricted family of
distributions `p` over configurations with `E: Z -> R`, `CVaR_alpha(p)` (mean of the lowest
alpha-fraction) attains its minimum `E_min = min_z E(z)` **iff** `p` places mass `>= alpha` on
`argmin E`. Everything `p` does with its remaining `1 - alpha` mass is **free**. So CVaR does not
design a population: it constrains an alpha-tail and is *indifferent* to the rest. Any diversity a
CVaR-trained sampler exhibits is a property of the ansatz's inductive bias and of partial
convergence — not of the objective. **This is a theorem, and it will be verified numerically
(`qb_theory`) before the sampler arms are read.** It generalises Sprint 17 §9b (which proved the
delta-collapse for the `mps2f` family on the lattice) to *any* family and *any* representation,
continuous included.

**(P-b) The terminal operator consumes the set MEAN.** LEDGER: `d_out = 1.16*d_set_mean +
0.04*d_set_best`, R2 0.89. A sampler that raises the *generation ceiling* (set best) is therefore
worth ~0.04 per Angstrom of ceiling through the deployed operator. **Prediction: a generation win,
if it exists at all, will not convert to a realised win.** The brief's §6 requires reporting both
regardless; this pre-registration commits to reporting a generation success *as* a selection
failure if that is what happens.

If both hold, the honest outcome is: the sampler hypothesis is closed, and the reason is
structural rather than an artefact of budget.

---

## 3. The representation — continuous torsion space, no lattice

**Configuration.** `z = (phi_1..phi_n, psi_1..psi_n) in T^{2n}`, the same object
`s15/align_lib.fit` optimises and the same object `core/project.build_ca_exact` builds from.
**No k=4 lattice, no torsion discretisation, no binary encoding of an angle.**

**Objective (deployed, frozen).**

    E(z) = sum_p ( d_p(z) - dhat_p )^2 / sd_p^2

over CA pairs with `|i-j| >= 2`, with `dhat` the leave-fold-out separation-debiased distogram
mean and `sd` its reported sd — i.e. exactly `s15/align_lib.fit`'s `f` at `kind="squared"`,
`wpair=None`. **The subset trap is respected**: `C.gather` and `C.fit_correction` are run on the
full `I.targets()` (n=126) and only then subset.

**The sampler family (this is the object under test).** Per residue `i`, the K=500 BLOSUM
retrieval pool supplies a *position-specific, native-free* empirical distribution over
`(phi_i, psi_i)`. Fit it with a **2-component von Mises mixture** (circular k-means init +
moment-matched concentrations, `kappa` capped so both components keep full support on the torus).
A latent bit `b_i` selects component `b_i`; the angle is then drawn **continuously** from that
component. The sampler's law is

    p_theta(phi, psi)  =  sum_{b in {0,1}^n}  p_theta(b) * prod_i q_{i, b_i}(phi_i, psi_i)

a genuinely continuous density on `T^{2n}` with **full support**, whose only discrete object is
the latent `b`. `p_theta(b)` is the measured distribution of a **genuine, exactly simulated
RY+CNOT circuit** (`core.quantum.MPSAnsatz`, `mps2f`: 2 entangling layers + a final RY layer),
trained by the **genuine sampled CVaR score-function gradient** (`core.quantum.cvar_gradient`,
`baseline="const"` — the corrected estimator, not the recorded `tail` defect) with
`core.quantum.Adam`.

**Why this and not something else.** What an entangled latent adds over a product latent is
*correlation between which conformer basin neighbouring residues occupy*. That is precisely the
structure the sprint brief §3 names — "a whole region flipping between two conformer families
moves many pairs the same way at once". If a structured quantum sampler is ever worth anything on
this problem, this is the construction where it should show. It is the strongest form of the
hypothesis available, and it is being built so that its failure means something.

**Gauge.** The labelling `b_i <-> 1-b_i` (which basin is called 0) is a `Z2^n` gauge. The
objective, the family of representable continuous laws, and every RMSD are invariant. The
*ansatz's inductive bias* need not be. **Pre-registered gauge test** (`qb_gauge`): rerun the
primary quantum arm under random per-residue relabellings and report the identity labelling's
percentile in the orbit, mid-rank tie handling (the Sprint-18 tie trap), on every target run.
**If the identity labelling sits below the 5th percentile of its own gauge orbit, any positive
result is declared encoding-dependent and is not claimable.**

---

## 4. Budget, and what counts as one evaluation

**One objective evaluation = one read of `E` at one continuous configuration.** Circuit
simulation, latent sampling, `grad log p`, and the von Mises draws are free. This is Sprint 14's
accounting rule, kept unchanged so the numbers are comparable across sprints.

**Every budgeted arm gets `B = 8192` evaluations.** Arms costing zero evaluations (the retrieval
pool, the torsion marginals, the constant helix) are reported with their zero and are *not*
credited with the saving; they are there as nulls, not as competitors on cost.

`shots = 512`, `iters = 16`, `lr = 0.15`, `baseline = "const"` — the Sprint-17/18 convention,
unchanged, so no hyperparameter is chosen here.

---

## 5. The arms

**QUANTUM (genuine, entangled, exactly simulated).**

| arm | what |
|---|---|
| `q_a005`, `q_a025`, `q_a100` | CVaR-VQE at alpha = 0.05, 0.25, 1.00 |
| `q_anneal` | alpha annealed 1.00 -> 0.05 (`core.quantum.alpha_schedule`) |
| `q_untrained` | **MANDATORY.** The *same* circuit at theta_0, B samples, **zero gradient steps** — best-of-N from the untrained ansatz |

**`q_untrained` is the control the LEDGER records a prior sprint failing to run** ("running the
VQE is WORSE than not running it, 0/12 cells, +0.65 to +1.32 A, when compared against best-of-N
from the untrained ansatz rather than against an initialisation mean"). **No arm in this lane is
ever compared against an initialisation mean.**

**CLASSICAL, at matched budget.**

| arm | what | why it is the right control |
|---|---|---|
| `c_prod005`, `c_prod025` | identical pipeline, **CNOTs deleted** (`mps2fn`), same CVaR estimator, same seed | isolates entanglement; a product latent has no inter-residue correlation |
| `c_chain` | a **classical learned proposal**: first-order Markov chain over the latent `b`, trained by the *same* CVaR score-function gradient at the same budget | the decisive control — a classical model with exactly the nearest-neighbour correlation the CNOT chain supplies |
| `c_cem` | cross-entropy method on the latent (product), same budget | the standard classical learned proposal |
| `c_metro` (T ladder) | continuous single-residue Metropolis thermostat, 3 temperatures | a classical thermostat that MOVES support |
| `c_anneal` | simulated annealing in continuous torsion space | |
| `c_lbfgs` | multi-start L-BFGS on `E` (the deployed refinement), function evals counted | search, not sampling — the incumbent optimiser |

**ZERO-INFORMATION / plausible nulls (0 evaluations).** Uniform-on-the-torus is **not** used as a
control (brief §8).

| arm | what |
|---|---|
| `c_marg` | i.i.d. draws from the same per-residue mixtures with basins drawn i.i.d. at the pool marginal rate — matched empirical torsion marginals |
| `c_helix` | constant ideal alpha-helix plus a matched isotropic torsion perturbation |

**INCUMBENT.** `pool75` — the shipped top-75 retrieval windows, 0 evaluations, the generation
mechanism actually deployed.

---

## 6. Outcomes — four, reported separately, never collapsed (brief §6)

For every arm, per target:

1. **generation ceiling** — `min` CA-RMSD over the arm's emitted set (ORACLE, post-hoc).
2. **the (M, D) plane** — common-frame member error `M` and diversity `D` via
   `s17/q_lib.md_plane`, whose exact identity `readout^2 = M^2 - D^2` is asserted to `< 1e-9 A^2`
   before any claim. Plus near-native coverage: `#{members within 1.5 A of the pool-best}` and
   `#{members < 2.5 A}`.
3. **selection ceiling** — `min` CA-RMSD over the `m = 75` members the **frozen, native-free**
   shipped Bayes-risk selector (`I.shipped_score`) keeps.
4. **realised** — CA-RMSD of `project(coordinate_average(those 75))`, the frozen terminal
   operator. Projection is the **repair** operator here (valid ideal geometry); AMBER is not run
   and no AMBER claim is made.

Plus **`aug`**: `pool500 ∪ arm's set`, selected and averaged identically — because a generation
success can only reach the pipeline as augmentation.

---

## 7. THE PRE-REGISTERED NUMERICAL DEFINITION OF "UNAVAILABLE TO CLASSICAL METHODS"

Fixed before looking at anything.

**Primary (P1) — realised.** `realised(best quantum arm) - realised(best classical arm)`, paired
over targets, target as the unit, fold-clustered bootstrap CI.

> **SUCCESS requires ALL of:** mean `<= -0.084 A` (the instrument's MDE); 95% CI excluding zero;
> the same sign in `>= 4/5` folds; **and** the quantum arm beating `q_untrained` by the same
> standard; **and** beating both zero-information nulls; **and** the gauge test not firing.

**Secondary (P2) — the Pareto frontier.** Additive epsilon-dominance `eps(p) = min_c max(M_c -
M_p, D_p - D_c)` (`s17/q_lib.eps_dominance`) of each quantum point against the **whole classical
point set**. `eps > 0` = no classical point dominates it.

> **SUCCESS requires:** mean `eps > 0` with a CI excluding zero, **and** the quantum arm's `eps`
> exceeding the leave-one-out `eps` of every classical arm scored against its own siblings.
> Sprint 17's recorded trap: a classical arm scored against its siblings looked exactly like the
> VQE, so a positive `eps` alone proves nothing.

**Tertiary (P3) — generation.** `generation ceiling(quantum) - generation ceiling(best classical)`
with the same standard, **reported separately from P1 even when they disagree**, and with its
min-of-N null (a matched-count uniform draw from the arm's own emitted band) named.

**Any positive on P2 or P3 with a null on P1 is reported as "generation success, selection
failure", per brief §6, and is NOT a pipeline result.**

---

## 8. Kill rules (brief §12) — declared now

Downgrade to CLOSED immediately if the arm:
loses to `q_untrained` · loses to `c_prod` (the CNOT-free control) · loses to `c_chain` (the
classical learned proposal) · loses to `c_marg` or `c_helix` (zero-information) · is matched by
`c_metro`/`c_anneal` (a trivial thermostat) · has a CI spanning the 0.084 A MDE · dissolves under
target-level analysis · fires the gauge test · needs native information anywhere in an inference
decision.

**A clean falsification is the deliverable if that is what the data says.**

---

## 9. Statistics

Target is the unit, n = 126 (or a pre-declared prefix if compute forces it — the prefix is
declared in the artefact and the claim label carries its n). Paired bootstrap CI (4000
resamples) plus the fold-clustered interval and the 5-fold sign vector. Median and W/L printed
beside every mean. `drop_top10` reported **against a uniform-effect null**, never as a raw
threshold (LEDGER: the raw drop-top threshold misfired once here).

**Labels:** EXACT / ORACLE / ESTABLISHED / SUPPORTED / PLAUSIBLE / OPEN / INCONCLUSIVE / REFUTED.
Every RMSD, every `M`, every coverage count and every generation ceiling is **ORACLE** and used
for post-hoc scoring only. `D`, `n_distinct`, the entropy, the objective values, and every
selection decision are **NATIVE-FREE**. No arm's alpha, temperature, budget, basin fit,
stopping rule or seed reads a native quantity.

**The sealed 60-target benchmark is not read, probed, or derived from. Nothing in this lane
touches `results/benchmark_manifest.json`.**

---

## 10. What I expect to happen (recorded so the outcome cannot be re-narrated)

I expect **P1 null or negative** (the theorem in §2 plus the mean-consuming operator), **P2
negative** (Sprint 17's plane result should port), and **P3 the only place with a chance** —
because a correlated basin latent genuinely can emit conformer-family flips a product sampler
cannot, and generation is the one axis this programme has never measured to saturate.

If P3 is positive and P1 is null, the honest headline is *the sampler raises the pool ceiling and
the deployed operator cannot spend it*, and the branch is closed on the pipeline while the
mechanism is recorded OPEN.
