# SPRINT 21 — LEDGER

---

## L1 — THE COORDINATOR'S OPENING EXPERIMENT IS DAMAGED THREE WAYS, ONE OF THEM A HAZARD HIS OWN BRIEF NAMES (2026-09-07, WORKSTREAM D)

**Benchmark seal verified BY HASH, not read**: `results/benchmark_manifest.json` sha256
`a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d`, 8002 bytes — **byte-identical to
Sprint 20's certificate.**

### 1. It did not run

`s21/results/tailprice.json` was 3 rows, **all skipped**, on
`ValueError("could not convert string to float: 'INWKGIAAMAKKLL'")`. Two wrong signatures:
`legacy_components_of_windows` is `(seq, PHI, PSI, chunk)` and was passed `(pdb, W, seq, fold)`;
`amber_energies` is `(space, S, ...)` where `Space` is the **discrete k-state lattice** and `S` an
integer state vector — **that interface cannot score real retrieved windows and reaching for it would
regress to the lattice the brief forbids.** The continuous route is `energy_from_coords` /
`refine_coords(steps<0)`. And the script **never loaded PHI/PSI at all** — only CA coordinates — so
neither physics arm had an input.

*(WORKSTREAM A has since fixed both properly, via the certified `s18/phys_down.py` coordinate path,
and recorded the diff as auditable while leaving the pre-registration unedited.)*

### 2. A VACUOUS COMPLETION FLAG — the exact hazard BRIEF §6 names, committed by the coordinator

`complete = len(rows) == len(tg)` **counts SKIPPED rows**, so **3/3 total failures wrote
`complete: true`.** It also does not require n = 126 — a 3-target smoke satisfies it.

> **The brief says "a completion flag must require the FULL configuration". The coordinator wrote the
> rule and then broke it in the first file he wrote under it.** The flag must require **n = 126 AND
> zero skips AND all seven Hamiltonians present.**

### 3. THE SCIENTIFIC ONE, and it survives fixing the other two

**(a) The primary endpoint was confounded.** `argmin` returns **one member**; `tail{a}` returns the
**coordinate average of a members**. They differ by the **averaging operator**, which Sprint 19
priced at **~1.0 A** ("where you average beats what you rank with"). **`tail(0.15)` would have beaten
`argmin` by a large margin on every Hamiltonian including a random one — and that margin is the
averaging operator, not the tail.** The pre-registered prediction was near-certain to "confirm" for
the wrong reason. **The matched `random{a}` control was the only line in the table worth reading, so
the PRIMARY should have been tail-minus-random, not tail-minus-argmin.**

**(b) A tail MEAN is not what a CVaR-VQE emits.** Verified in source by the coordinator:
`core/quantum.py` returns `vqe_bitstring` (argmin over the final distribution's samples),
`vqe_modal_bitstring` (the mode), `best_seen_bitstring` (argmin over everything seen). **None is an
average of a tail. CVaR is the TRAINING objective; the READOUT is an argmin.** The docstring named an
operator **the deployed pillar does not have.**

**(c) AND THE SHARP CONSEQUENCE — a near-theorem that may answer the sprint's opening question on
paper.** Q6 is EXACT: CVaR's minimiser is a **face supported on the alpha-tail**. *The alpha-tail of
H contains H's global pool minimum.* Therefore:

> **For a pool-restricted selector whose training objective and readout energy are the same H,
> CVaR-VQE and argmin have the SAME optimal answer. CVaR changes the sampling distribution, not the
> selected point.**

**Scope conditions**: requires same-H readout and pool restriction; **breaks** if the readout is an
**average**, or if the **readout H differs from the training H**. WORKSTREAM D is formalising it.

**What this leaves live**: the **averaged** readout — which is the Sprint-19 averaging lever wearing
a CVaR label and must be labelled as such — and a **readout H different from the training H**, which
is genuinely new and untested.

### Disposition

BRIEF §2 corrected in place. The run is **kept** and **re-primaried on `tail − random{a}`**, with
**single-structure readout arms added** (medoid of the tail, and a random member of the tail) so the
**tail operator and the averaging operator are separated**, and the averaged arms labelled as what
they are. **Coordinator error #10, and the second of the "wrote the rule then broke it" kind.**

---

## L2 — THE ENCODING LEVER IS CONFOUNDED BY STEP COUNT, and one of its twelve cells is an identity (2026-09-07, WORKSTREAM D)

Sprint 20 reported `theta` vs `(cos theta, sin theta)` — physically identical — winning on RMSD in
**11/12 cells at p = 0.006**, and the coordinator called it *"the cheapest RMSD lever anyone has
found this programme"*. **Three defects, in descending order of damage.**

### 1. The step-count confound, present in EVERY cell

From the artefacts' own `iters` field:

    arm        embedded   theta    ratio    dRMSD(emb - theta)
    adam_fd         5.1    10.5    0.49     -0.287 / -0.369 / -0.191
    spsa          159.2   207.6    0.77     -0.247 / -0.649 / -0.131
    lbfgs_fd     (scipy, runs to convergence) -0.048 / -0.125 / +0.073

**The embedded arm takes 23-51% fewer steps, and reaches a WORSE objective in 9/12 cells.** On an
instrument whose most reproduced law is that **optimising the deployed objectives harder makes the
structure worse** (S20 Q3: `c_lbfgs`, the arm optimising hardest, is second-worst on the board),
**that is the signature of an arm that optimised less, not of a better coordinate system.**

The Sprint-20 findings **declared** the halved-iteration confound for the FD arms and then wrote that
it *"does not touch the RMSD column"*. **On this instrument that inference is backwards.**
**Supporting sign**: `lbfgs_fd` — the arm least sensitive to step count, since scipy runs to its own
convergence — has **the three smallest effects including the only positive one.**

### 2. SPSA's deficit is a DIAGNOSTIC eating the budget asymmetrically

`track_quality` spends 3 central-FD probes of `2*d` units each, and **`d` doubles in the embedding**,
so the probe consumes **288 of 512 budget units embedded against 144 in theta** — 28% vs 14%, **in
the two arms of the headline comparison.**

### 3. Two of the twelve cells are the SAME measurement

`AMB|nelder` and `AMBc|nelder` are **bit-identical on all 10 targets x 2 seeds in both encodings**.
`AMBc = sign(E)*log1p(|E|)` is **strictly monotone** and Nelder-Mead is **comparison-only** — as is
`best_of_N`. So the sign test counted 12 independent Bernoulli trials where there are **11 distinct
cells: 10 negative, p = 0.0117 not 0.0063** — and **one of the three "CIs excluding zero" is the
duplicate**, leaving two.

> **An identity counted as an empirical trial** — BRIEF §9's own rule. Not a fabrication; the
> transform is documented in the source. **ACTION: exclude every comparison-only arm from any
> AMB-vs-AMBc conditioning table. Those rows are exactly zero by construction and are not nulls.**

### What survives, at the mandated unit

Redone target-level over the 11 distinct cells, bootstrap over 10 targets:

    mean -0.2167   median -0.2586   CI95 [-0.3239, -0.1095]   W/L 9/1   sign p = 0.0215

> **The MEASUREMENT is real and survives the correct unit. Its INTERPRETATION is not established:
> nothing in Sprint 20 separates "a better coordinate system" from "an arm that took 23-51% fewer
> steps".**

### A mechanism for the confound, and a second hazard in the same claim

In the embedding `||u||` is a **pure gauge direction**: `arctan2` is scale-invariant, so the
objective's true gradient has **exactly zero radial component** and **half the embedded parameters
are unphysical**. `adam_fd` normalises per coordinate by `sqrt(v)`, so a direction whose gradient is
**pure finite-difference noise receives a full-size Adam step**. If `||u||` diffuses upward, the
effective **angular** step shrinks as `1/||u||` — **an implicit, unintended step-size annealing
schedule**, which is a mechanism for exactly the observed effect.

### The deciding control, pre-registered by the lane before running it

**Matched EFFECTIVE STEP COUNT, in both directions** — rule 1, matched in the space the operator
works in: `theta_matched` cut to the embedding's iteration count, `emb_matched` raised to theta's,
4 seeds (Sprint 20 used 2). Plus two instruments the artefacts do not carry: **total angular
displacement** `||wrap(z_best - z0)||` and the **gauge-radius drift** `||u||`.

**Prediction**: the embedding's angular displacement is smaller and the gap closes to inside the MDE
once steps are matched both ways. **Falsifier**: if `theta_matched` does not reproduce the
embedding's RMSD and `emb_matched` does not lose its advantage, **the lever is real and the objection
is withdrawn.**

**Disposition**: BRIEF §5 item 3 downgraded in place; WORKSTREAM B told **not to design an ansatz
around the lever until the control reports**, and redirected to the bond-dimension question and QNG,
both unaffected.

---

## L3 — "WHICH HAMILTONIAN IS BEST" IS NOT WELL-DEFINED WITHOUT A READOUT, and the tail-average selector IS the incumbent (2026-09-07, WORKSTREAM D)

`s21/results/d_cvarop.json`, **n = 126, COMPLETE**.

### 1. The strongest soundness gate available here, and it reframes the matrix

    tail_avg(alpha=0.15) | disto  =  3.0483   vs shipped rmsd_avg:  max |diff| = 0.000000 on 126/126, corr 1.000000
    argmin | disto        =  3.454   vs the instrument's pinned `shipped argmin 3.4540`

**Both of the instrument's own pinned constants return exactly.** And they must: **the shipped
pipeline IS "rank the K=500 pool by the Bayes-risk distogram score, keep the top 75 (= 0.15 x 500),
coordinate-average them".**

> **The mandatory matrix's `disto` row at alpha = 0.15 is not a new arm — it is the production
> pipeline wearing a CVaR label.** "argmin vs tail" was settled in production five sprints ago: the
> pipeline already takes the tail-average, **which is exactly why it sits at 3.048 and not at the
> 3.454 argmin.** A CVaR-VQE that selects an alpha-tail of this pool and averages it is **re-deriving
> the incumbent.** *Do not report 3.048 as a CVaR result.*

### 2. The coordinator's original endpoint is maximised by an energy containing NOTHING

`tail_avg(0.15) - argmin`, paired over 126 targets, iid CI with a fold-clustered CI beside it:

    disto   -0.406 [-0.536,-0.276] iid   [-0.513,-0.317] fold   W/L  95/31
    legacy  -0.735 [-1.113,-0.369] iid   [-1.108,-0.402] fold   W/L  92/34
    rand    -1.035 [-1.126,-0.946] iid   [-1.207,-0.914] fold   W/L 123/3

`rand` is a **per-target random permutation of the pool**, averaged over 16 permutations — **zero
information — and it produces the LARGEST gap of the three: 181% of the physics mean and 1.41x the
weakest physics arm.** The lane's pre-registered falsifier (refuted if rand's gap is under half the
physics gaps) **does not fire**.

> **The endpoint measures the READOUT, not the Hamiltonian, and cannot rank Hamiltonians at all.**
> Decomposed with the tail held fixed: the **averaging** half is -0.28 to -0.31 for **all three
> energies including rand** (energy-independent by construction), and the **set-vs-min** half is
> -0.74 for rand because **argmin-of-noise is a random pool member.** *Neither half is physics.*

**This confirms the coordinator's re-primary to `tail - random` was correct, and quantifies how badly
the original endpoint would have misled.**

### 3. THE FINDING WORTH KEEPING — Legacy's sign FLIPS with the readout

Against **matched-count random controls** at every alpha (the only comparison that means anything;
`rand`'s arms sit at 0.00 +/- 0.02 against the same controls, **which validates that the controls are
matched**):

    readout        disto                    legacy
    tail_member   -0.90 [-1.08,-0.73]      -0.41 [-0.56,-0.26]   both BEAT their control
    tail_medoid   -0.45 [-0.67,-0.24]      +0.30 [+0.20,+0.40]   legacy WORSE, 35W/91L
    tail_avg      -0.38 [-0.55,-0.22]      +0.33 [+0.21,+0.45]   legacy WORSE, 41W/85L

Same pattern at alpha = 0.05 and 0.30, all CIs excluding zero.

> **Legacy's tail is BETTER than a random subset if you draw ONE member from it, and WORSE than a
> random subset if you take its medoid or its average.**
>
> **Mechanism, and it is Sprint 20's L2c**: Legacy prefers **compact, pool-typical** geometry, so its
> tail is a **tightly clustered, low-diversity set** — which raises the mean member **and destroys the
> error diversity the averaging operator lives on.**

**CONSEQUENCE, now binding on the mandatory matrix**: *"which Hamiltonian is best" is not well-defined
without fixing the readout*, and **every row must state its readout.** A matrix run at `tail_avg`
reports Legacy as **harmful**; the same matrix at `tail_member` reports it as **helpful**. **Both are
true.**

### 4. The EXACT limit, now measured as well as derived

    tail_min == argmin in all 9 (H, alpha) cells, max |difference| = 0

Q6 says a CVaR-optimal law is supported on the alpha-tail; the alpha-tail **contains** the pool
argmin; and the deployed readouts are **order-based** (`vqe_bitstring`, `vqe_modal_bitstring`,
`best_seen_bitstring`, verified in source). **So a pool-restricted CVaR-VQE with an order-based
readout and one Hamiltonian cannot select anything argmin cannot.**

**Scope, stated by the lane**: needs pool restriction, the same H for training and readout, and an
order-based readout. **It BREAKS if the readout averages** — and then **the averaging operator, not
CVaR, is doing the work**, which is point 2.

> **The lane explicitly does NOT claim this closes the CVaR arm.** What it closes is the
> **pool-restricted, same-H, order-based** case. **Live**: a readout H *different* from the training
> H; an *averaging* readout (where credit belongs to averaging); and the *unconverged* law — for
> which `tail_member` is the right price, which is why it was added.

**Ceilings, ORACLE and labelled**: ORACLE tail-average at alpha = 0.05/0.15/0.30 is **1.644 / 1.963 /
2.275** (POINT CLOUD basis); `argmin|rand` = 4.453.

### Operational

The box has been at **95-96%** while `core.amber.memory_guard` refuses above **92%**, and **two
lanes' runs have already been killed by it.** OpenMM contexts are now **serialised** in BRIEF §11,
priority **A's matrix > C's continuation > D's AMBc half**, with every lane instructed to run its
Legacy/distogram half first and queue the AMBER half.

---

## L4 — NO BOND DIMENSION BREAKS CLASSICAL SIMULABILITY HERE, BECAUSE THE REGISTER IS TOO SMALL (2026-09-07, WORKSTREAM D)

`s21/results/d_sim.json`, **64 cells, COMPLETE**, all normalisation checks pass.

### The question was the wrong one, and that is the answer

Sprint 20's Q7 — the deployed ansatz is **bond-dimension 4, a <=16-state HMM** — is **true and is NOT
the binding constraint.** From source (`s20/qb2_opt.py:322`) the circuit is
`MPSAnsatz(n, layers=2, final_ry=True, entangler="cnot")` with **n = residues**: **one qubit per
residue**, and the instrument is **9 <= n <= 16**. The whole Hilbert space is **512 to 65,536 basis
states.**

The lane enumerated the circuit's **complete** probability distribution **exactly**, for every
register size the 126 targets use, at layers 1, 2, 3, 4, 6, 8, 12, 20:

    n=16, layers=2  (DEPLOYED)                        65,536 states   654 ms via MPS    sum p = 1.0000000000
    n=16, layers=20 (chi SATURATED at chi_max = 256,
                     maximal entanglement for this
                     topology, ring closure included)  65,536 states    64 ms via dense  sum p = 1.0000000000
    every register size, every depth, total                                              17.7 s

> ### The dense route costs 2^n and is INDEPENDENT of depth. There is NO bond dimension at which classical simulability breaks for this circuit — not at chi = 4, not at chi = 256, not at any depth. The register is too small, full stop.

**And the parameter count sharpens it**: **3n = 27 to 48 parameters over 2^n = 512 to 65,536
outcomes.** The ansatz is a **restricted parameterisation of a small classical categorical
distribution**, which has 2^n - 1 free parameters and can be written down and optimised directly.

### Trainability, and the literature that bites

To reach a register where classical simulation is genuinely hard you need **n >~ 50 qubits — 50-residue
chains, outside this benchmark entirely.** There, **Cerezo et al., "Does provable absence of barren
plateaus imply classical simulability?" (arXiv:2312.09121, *Nature Communications*, Aug 2025)**
applies: the structure that **removes** barren plateaus is the same structure that **admits classical
simulation**, because barren plateaus are a curse of dimensionality and current fixes encode the
problem into a small classically simulable subspace — a soft dequantisation.

**The lane stated their caveats rather than hiding behind the headline** (average-case arguments,
smart initialisations, models outside their assumptions, and a possible quantum data-acquisition
phase) **and noted none of them help here — because here there is no data-acquisition phase to need:
the law is tabulable in 64 ms.**

**Separation on a comparable continuous encoding: none found.**

### WHAT THIS DOES NOT SAY — both halves are binding

> It does **not** say the VQE is useless and it does **not** touch Q4: **the VQE genuinely trains**,
> and **a 48-parameter restricted model can be a perfectly good optimiser.**
>
> It says that **nothing measured on this instrument can ever be evidence of quantum resource**, and
> that **no ansatz redesign changes that while the register is one qubit per residue on 9-16-mers.**

**Consequence for WORKSTREAM B**, already routed: RMSD-facing ansatz comparisons proceed unaffected;
**any quantum-resource or advantage component of the night is CLOSED and must be reported as closed
rather than as untested.**

---

## L5 — D1 LEG HALF: the step confound is CONFIRMED as a diagnostic artefact, and the Legacy effect does not replicate (2026-09-07, WORKSTREAM D)

`s21/results/d_enc_LEG.json`, COMPLETE. n = 10 targets, **4 seeds** (Sprint 20 used 2),
`track_quality` **OFF**, budget grid {256, 512, 1024} in both encodings. **Gate: 32/32 bit-for-bit
reproductions of Sprint 20's own stored rows before any number was read.**

**MECHANISM CONFIRMED.** With the diagnostic off, **SPSA's iteration counts are EXACTLY EQUAL in the
two encodings** (256 vs 256 at B = 512). **So SPSA's entire 159-vs-208 step deficit in Sprint 20 was
the `track_quality` diagnostic.** `adam_fd` keeps a 0.5 ratio — **that one is the genuine 4n-vs-2n
dimension effect and cannot be removed at fixed budget**, which is why a grid was measured.

**GAUGE DRIFT IS REAL AND LARGE.** `||u||` starts at 1.000 and ends at **1.94 (max 3.06)** for
SPSA@512, and **2.67 (max 3.79)** at B = 256 — **in a direction whose true gradient is exactly zero.**
The effective angular step therefore **decays ~2x over a run: an unintended annealing schedule.**

**OUTCOME ON LEGACY — the effect does not replicate.**

    A  matched budget   spsa +0.215 [-0.055,+0.507] 5W/5L    adam -0.125 [-0.301,+0.038] 7W/3L
    B1 matched steps    spsa +0.215 (steps already equal)    adam -0.052 [-0.139,+0.040]
    B2 matched steps    spsa +0.215                          adam -0.071 [-0.225,+0.071]
    pooled, target level  +0.045 [-0.109,+0.211]   W/L 5/5   sign p = 1.00

SPSA's point estimate **flips sign**; adam shrinks from -0.191 to -0.05/-0.07, **inside the MDE.**

**HONEST LIMITATION, VOLUNTEERED BY THE LANE AND NOT HIDDEN**: *Legacy was all-ns in Sprint 20 too,
so this is not yet a knockdown of the headline* — the significant Sprint-20 cells were **AMB and
AMBc**. **And the lane's own premise ("fewer steps helps here") is NOT MEASURED on Legacy**: panel C,
the step-count response within each encoding, is +0.134 for adam/theta and -0.164 for spsa/theta,
both CI-spanning. **AMBc decides it, and the lane has committed to report it either way, including if
it falsifies the objection.**

---

## L6 — THE 2026 CONTINUOUS-SPACE PREPRINT, read as method rather than headline (2026-09-07, WORKSTREAM D)

**arXiv:2609.02113**, *"Logarithmic-scale VQE for off-lattice protein structure prediction in
continuous torsional angle space"* (2 Sep 2026). Three points, from the method rather than the number:

**1. Its own numbers carry OUR basis problem, in one sentence.** *"Chignolin reached 0.623 A Ca-RMSD
in RETAINED SNAPSHOTS and 1.199 A in FINAL MODELS; Trp-cage 2.501 A among snapshots, 3.512 A in final
models."* **A 2x and a 1.4x gap between a best-of-set and what the pipeline emits** — the same object
as this programme's ~1.9 A selection gap, reported as two numbers in one sentence.

**2. It independently reproduces this programme's central negative.** *"Energy-ranking imbalances
persisted across sampled landscapes FOR ALL FUNCTIONS."* That is `nothing-ranks-within-the-pool` /
`the-objective-does-not-rank-the-native`, found by an **external group, on different peptides, with
different energies. A third external corroboration.**

**3. The "O(log2 N) qubits" claim is a TRADE, and the paper says so in its last line**: *"By
converting physical qubit constraints into circuit DEPTH constraints."* Reading N torsions off
`2^q = N` basis-state probabilities **does not compress information — it moves the cost into
measurement shots**, needing **O(N^2/eps^2)** to resolve N probabilities each of size ~1/N. And the
simulator readout **extracts torsions from relative PHASES in statevector simulations**, which is
**not measurable on hardware at all** — which is why the hardware runs use a different CDF decoder.
**Two different algorithms under one name**; the hardware number (best **1.758 A** on ibm_cleveland /
ibm_miami) is the only one that prices the hardware path.

**No control against an untrained circuit or best-of-N is reported.**

---

## L7 — R1 (DIVERSITY-MAXIMISING SET SELECTION) IS DEAD, AND KILLING IT RETRACTS THE MECHANISM THE LANE HAD ASSERTED (2026-09-07, WORKSTREAM D)

`s21/results/d_divsel.json`, **n = 126, COMPLETE**, point-cloud basis, fold-clustered CIs. Anchored
to the production number: the incumbent top-75 is **3.0483**, which `d_cvarop` reproduced bit-exactly.

**The proposal.** Every selector in this project ranks candidates **independently** and takes a
top-*k*; **nothing has ever selected the averaging set AS A SET.** Falsifier written before the run:
dead unless `divmax` beats **both** the incumbent top-75 **and** the random-in-band control by more
than the MDE with a CI excluding zero, **and** the sign control (`divmin`) loses.

    band=top30%   divmax 3.0587   vs top75 +0.010 [-0.041,+0.057]   vs rand_in_band -0.015 [-0.044,+0.008]
    band=top50%   divmax 3.0911   vs top75 +0.043 [-0.020,+0.124]   vs rand_in_band -0.046 [-0.115,+0.038]
    band=all      divmax 3.7795   vs top75 +0.731 [+0.571,+0.894]   vs rand_in_band +0.347 [+0.241,+0.445]

**All three bands dead. The falsifier fires exactly as written.**

### The sign control goes the WRONG WAY, and that is the informative part

    band=all   divmin (MINIMISE diversity)  vs rand_in_band  -0.2140 [-0.3974, -0.0400]  66W/60L

> **Minimising set diversity HELPS. So diversity is not the lever — and the mechanism the lane had
> asserted to the coordinator two messages earlier was WRONG. It retracted it here rather than
> leaving it in the record.**

### The corrected mechanism, which reconciles both results

**It is not the SPREAD of the selected set; it is whether its errors are COHERENT.** Averaging
cancels i.i.d. error and **preserves systematic** error — `error-coherence-decides-correctors`: at
identical 0.688 sign accuracy, **coherent mistakes emit +0.31 A and i.i.d. mistakes -0.14 A.**

Legacy's tail shares a **coherent compactness bias** (L2c: Legacy-preferred candidates are 0.45 A
more compact, 124W/2L), **so averaging preserves the bias** — which is why Legacy costs **+0.33 on
`tail_avg`** while **helping `tail_member` by -0.41**. And **maximising spread does not decorrelate
errors, it selects OUTLIERS**, which is why `divmax` collapses to 3.78 on the unbanded pool; while
**minimising spread is outlier avoidance** — `consensus-is-outlier-avoidance` reproduced on a new
operator.

> **L3's finding stands exactly as measured** — the sign of a Hamiltonian flips with the readout,
> n = 126, CIs excluding zero — **but the reason is error COHERENCE, not set diversity.** BRIEF §4
> corrected in place.

### One thing worth keeping from the wreckage

**The score band itself is load-bearing**: a random 75 drawn from the top 50% **loses to the top-75
by +0.0883 [+0.0113, +0.1631]**, CI excluding zero. **The distogram's ordering is doing real work
inside the band** — consistent with L3's `tail_member` at -0.904 against its control.

> **Fourth independent sighting of the same fact: the distogram orders in-band, and nothing else this
> project has tried does.**

### The lane downgrades its own next proposal

Its R2 (per-target readout gating by tail spread) is **weakened by its own data**, because spread is
not the driver. **Recorded as OPEN-but-unmotivated rather than recommended**, with the note that any
such gate would need an **error-coherence** proxy and **no native-free one exists.**

### Compute, and a correct self-throttle

The lane **killed its own AMBc budget grid after one target** — it was running at **~25 min/target
under four-way contention (52 ms per AMBER single point against a 6 ms nominal)** and would have
taken four hours while starving another lane that was logging *"free 0.42 GB < 1.8; waiting"*. The
partial is preserved as `_PARTIAL_d_enc_AMBc_fullgrid_1target.json`.

**It restarted minimal and decisive** — AMBc, spsa, budget 512, 4 seeds — **and the reasoning is
right**: with `track_quality` OFF the SPSA step counts are **exactly equal** in the two encodings, so
**the matched-budget comparison IS the matched-step comparison and no grid is needed.** AMBc/spsa
**-0.649** is the largest Sprint-20 encoding effect and one of only two distinct cells with a CI
excluding zero — **the cell that decides D1.**

---

## L8 — "MDE = 0.084 A" IS NOT A PROPERTY OF THE INSTRUMENT, AND FOUR SPRINTS HAVE QUOTED IT AS ONE (2026-09-07, WORKSTREAM D)

`s21/results/d_mde.json`, **COMPLETE**, **26 real comparisons** drawn from this sprint's artefacts,
Sprint 20's, and the production cache.

### The operator

    MDE = (z_.975 + z_.80) * SE = 2.8016 * sd(paired differences) / sqrt(n)

**`sd(paired differences)` is a property of THE COMPARISON, not of the instrument.** The paired sd
that would make 0.084 correct at n = 126 is **0.3366 A**. *Almost nothing has that sd.*

**Measured: per-comparison MDE / 0.084 ranges 0.11x to 9.54x — a factor of 84. Only 12% fall within
1.5x of the quoted constant.** The lane's pre-registered falsifier (refuted if they cluster within
~1.5x) **does not fire.**

### It errs in BOTH directions, and both errors are already in the record

**TOO LARGE for low-variance comparisons.** The deployed AMBER repair (`rmsd_full - rmsd_arm`,
n = 126): paired sd **0.0385**, SE **0.0034**, **MDE 0.0096 A — one eleventh of the quoted bar.** Its
effect is **+0.0207 = 6.09 SE, power 1.000.**

> **The Sprint-21 brief calls this "a quarter of the MDE", which invites reading a SIX-SIGMA effect as
> negligible.** It is not negligible — it is one of the most precisely measured quantities in the
> programme. Same for the projection tax (MDE 0.051 and 0.048 against effects +0.166 and +0.156).
> **Coordinator error #11**, and it is in the live brief.

**TOO SMALL for high-variance comparisons, and this one is live.** Sprint 20's encoding cell
AMBc/spsa (n = 10, **2 seeds**): paired sd **0.9045**, SE **0.2860**, **MDE 0.8014 A — 9.54x the
quoted bar** — while its reported effect of **-0.649 [-1.205, -0.165]** was flagged significant.

> **The study's own 80%-power MDE is LARGER than the effect it detected.** |d|/SE = **2.27**, post-hoc
> power **0.62**, Type-M magnitude exaggeration conditional on significance **1.27x**. The sign is not
> in doubt (Type-S ~ 0) but **the magnitude is biased upward — exactly what an underpowered design
> produces when it reports a significant result.** The LEG/spsa cell is worse: **power 0.18, Type-M
> 2.39x.**
>
> **That is a third, independent reason the Sprint-20 encoding panel is fragile**, on top of the
> Nelder-Mead duplicate (L2) and the step-count confound (L2, L5).

**AND IT CORRECTS THE LANE ITSELF FIRST.** Its own Legacy null (pooled emb@512 − theta@512, n = 10,
4 seeds) has **MDE 0.242 A, not 0.084.** So its "NOT MEASURED on LEG" label is right **but its bar was
wrong** — that null is uninformative below **0.242**, not below 0.084 — and it weakened its own
conclusion accordingly before telling anyone else.

### The rule, adopted into BRIEF §1

> **Report SE beside every mean, and quote the MDE that comparison's own paired sd implies. Never
> label an effect against the 0.084 constant.**

**The programme derived this correctly once and then lost it**: `s16/review_FINDINGS.md` Q4.1 computes
*"SE = 0.051 A, MDE at 80% power ~ 0.16 A, power at a true 0.05 A effect ~ 9%"* for its own panel, and
`s16/LEDGER.md:678` records it the same way. **The CONSTANT is what generalised into Sprint 18's
brief; the METHOD is what should have.**

### Scope — what this does and does not touch

**Does not touch** the sprint's large results, which are many SE from zero on their own sd: signflip,
the in-manifold contrast, Legacy's compactness axis, the quantum-vs-null gap, the readout sign flip
(all CIs excluding zero on their own paired sd). **Does touch** every claim labelled "below the MDE",
"a quarter of the MDE", or "matched" against the constant — **each must be re-labelled against its own
comparison's SE**, and some will move in *both* directions.

---

## L9 — LEGACY CARRIES NO IN-BAND RANK INFORMATION BEYOND THE DISTOGRAM (2026-09-07, WORKSTREAM D)

`s21/results/d_partial.json`, **COMPLETE**, computed **at zero cost from coefficients already in
WORKSTREAM A's live artefact**.

**The claim it kills.** A's report prints as a headline distributional fact
`ORACLE_rho_dis_d +0.633 | ORACLE_rho_leg_d +0.335 | ORACLE_rho_amb_d -0.019`, which reads as
*"genuine Legacy carries about half the distogram's in-band rank skill"*. **Three lines above, the
same artefact records `rho_dis_leg` = +0.645 median.** Legacy and the **deployed selector** are
strongly rank-correlated, and **a marginal correlation with the truth, for a score already correlated
with the score in production, is not that score's contribution.**

Per target, from that target's own three coefficients, aggregated fold-clustered, n = 42:

    rho(legacy, true)           MARGINAL   +0.3281  fold [+0.3011, +0.3552]    8/42 negative
    rho(disto,  legacy)                    +0.4784  fold [+0.4226, +0.5240]
    rho(legacy, true | disto)   PARTIAL    -0.0076  fold [-0.0616, +0.0681]   25/42 negative

> ### Legacy retains −2% of its marginal rank skill. It carries NO in-band rank information beyond the distogram.

The lane's pre-registered falsifier (refuted if the partial keeps more than half, > +0.16, with a CI
excluding zero) **does not fire**.

**CONTROL, so the partialling is not manufacturing structure**: AMBER's marginal is already ~0 and
its partial **stays** ~0 (-0.0290 [-0.0593, -0.0036] fold).

**INDEPENDENT CROSS-CHECK, already in A's own tables**: `d+l` never beats `d` at any alpha or any
readout. **This is the mechanism for that.**

**Same error class** as Sprint 20's Q11 (every circuit-side landscape metric collapsed once target
difficulty was partialled out) and the shared-referent floor: *a marginal statistic reported where a
partial is the meaningful quantity.*

**CONSEQUENCE, now binding on §4**: the mandatory matrix's **Legacy row must be reported as a
PARTIAL, not a marginal**, or Legacy will appear to contribute when it does not.

---

## L10 — THE READOUT SIGN-FLIP IS A TWO-LANE RESULT, and the pipeline identity is confirmed in three places (2026-09-07)

Two **independent implementations, different n, different code**, agreeing to the third decimal:

                              audit lane (n=126)         WORKSTREAM A (n=42)
    tail_member | legacy      -0.407 [-0.557,-0.259]     -0.406 [-0.682,-0.122]
    tail_medoid | legacy      +0.302 [+0.204,+0.404]     +0.305 [+0.097,+0.530]
    tail_member | disto       -0.904 [-1.078,-0.729]     -0.874 [-1.166,-0.572]

**Legacy beats a matched random subset on a single-member readout and loses to it on the medoid, in
both lanes. L3 is no longer one lane's result.**

**And the pipeline identity is confirmed in THREE places** — A's `tail0.15|disto`, the audit lane's
`tail_avg0.15|disto`, and production's `rmsd_avg` — at **max |diff| = 0.000000 on all 42 targets A
has finished.** *The alpha = 0.15 disto row IS the incumbent.* Adopted as a standing soundness gate.

### AMBER's selection row, settled

From A's table: **AMBER is worse than its matched-count random control at essentially every readout
and every alpha** — averaged tail **+1.614 / +0.998 / +0.485** at alpha 0.01 / 0.05 / 0.15,
**0/5 folds negative**, medoid up to **+2.008** — and its ORACLE rank skill against true RMSD is
**-0.042 [-0.116, +0.032]**, 23/42 targets negative.

> **AMBER contributes nothing to selection over this pool and actively HURTS at small alpha.**
> Sprint 20's L4c ("composites never beat their own first stage") reproduced on a tail selector.

### A live label defect, caught before the table was quoted

`tailprice.py` annotates its rank-skill block *"(negative = lower energy is better)"*. **The sign is
inverted**: `disto` is +0.633 with 2/35 targets negative, so **positive rho means correct ordering**.
As printed, **the legend says the arm that works is the arm that fails.** Routed to A.

### And the lane applied its own MDE finding to its own design

It started a **properly powered LEG+DIST replication at n = 60**, because L8 shows an n = 10 panel's
own MDE on this comparison is **0.24-0.80 A** — so **neither Sprint 20's panel nor its own first one
could ever have resolved the encoding question.** It also killed a stray full-grid process of its own
that `TaskStop` had orphaned and which was competing for memory with A's run, and said so.

---

## L11 — THE BUDGET EXCEEDS THE ENTIRE LATENT SPACE ON 60% OF TARGETS (2026-09-07, WORKSTREAM D)

From source, `s19/qb_lib.py:198`: `draw_from_basins(bits, ...)` takes bits of shape `(B, n)` with
**n = RESIDUES** and bits in {0,1}. **One qubit per residue, two basins each.** So the **entire latent
space the VQE explores is `2**n`**, with n the peptide length.

    |latent| across the 126 tuning targets:   min 512    MEDIAN 8192    max 65,536

Against the budgets this programme actually spends:

    budget   512:  2**n <= budget on    9/126 targets ( 7%)
    budget  4096:  2**n <= budget on   52/126 targets (41%)
    budget  8192:  2**n <= budget on   75/126 targets (60%)

> ### At Sprint 20's 8192-evaluation budget, the budget EQUALS OR EXCEEDS the entire latent space on 75 of 126 targets — and the MEDIAN target's latent space is exactly 8192.

**So on a majority of targets the deployed VQE is not searching a space it cannot enumerate — it is
RESAMPLING a space it could have enumerated within budget.**

**This changes what Sprint 20's Q3 means.** *"No sampler, quantum or classical, at 8192 evaluations
beats the zero-evaluation retrieval pool"* is **not** a statement that search was too hard on those
targets. **Exhaustive enumeration of the latent was affordable and would not have helped.** What binds
is **the objective's ORDERING over the latent**, not the sampler's reach.

> **`search-saturates-discrimination-binds` now has an EXHAUSTIVENESS argument attached, not only six
> empirical instruments.**

### Two consequences, both adopted

1. **Any budget-scaling arm on this latent must state `|latent| = 2**n` beside the budget.** A sweep
   from 512 to 8192 **crosses the point where the budget swallows the whole space on 60% of targets**,
   and the sweep's shape past that point is a **resampling artefact, not a search curve.**
2. **The cheapest possible control for any sampler on this latent is EXHAUSTIVE ENUMERATION** on the
   **n <= 13** targets (75 of 126), costing `2**n <= 8192` evaluations — **no more than the budget the
   samplers already receive.** It yields the **EXACT argmin of the objective over the latent**, which
   is the true optimum of the search half. **If the exact latent argmin does not beat the pool either,
   the search half is closed BY CONSTRUCTION on those 75 targets rather than by six instruments.**
   *Taken by the coordinator; see L12.*

**This reinforces the reachability verdict in the same artefact**: 9-16 qubits, complete law tabulable
in **64 ms at maximal entanglement**, **27-48 circuit parameters over 512-65,536 outcomes**, and now a
**budget that exceeds the whole space on most targets.**


---

## L12 — THE EXHAUSTIVE LATENT CONTROL IS BUILT AND RUNNING (coordinator)

`s21/latentfull.py`, launched over the **75 targets with n <= 13**. It enumerates **every one of the
`2**n` latent configurations** and takes the **EXACT argmin of the deployed objective** — the true
optimum of the search half, at a cost (`2**n <= 8192`) no greater than the budget the samplers
already received.

**One design fact governs the whole file.** `draw_from_basins` (`s19/qb_lib.py:196`) maps a bitstring
to a **von Mises DRAW**, not to a structure, so `2**n` bitstrings are `2**n` **distributions**, and
there is no such thing as "the" structure of a bitstring. The deterministic exhaustive object is
therefore each bitstring decoded at its **von Mises MODE** — which is also exactly the set an
argmin-readout VQE is choosing among. A stochastic arm (one draw per bitstring) runs beside it so
the mode restriction is **priced rather than assumed**.

**The objective is not reimplemented.** The file calls `Obj.raw` — the deployed objective class
itself, unbudgeted, which is what its own docstring provides it for. This removes the possibility of
the control being certified against a lookalike.

**The trap this design had to avoid, because the sprint has already fallen into it twice.**
`latent_argmin` is a SINGLE STRUCTURE; `pool_tail_avg` is an AVERAGE OF 75. Sprint 21 L3 measured
that averaging operator at **-0.28 to -0.31 Å independent of the energy**. Comparing the exhaustive
argmin to the incumbent average would hand the pool that margin for free and manufacture the
conclusion. **The primary is therefore stated against `pool_argmin` — argmin against argmin, the
same readout on both sides** — and the unmatched row is printed only to display the margin it would
otherwise have absorbed.

**Pre-registered, before the run**: hypothesis `latent_argmin >= pool_argmin` with
`latent_oracle << latent_argmin`; **falsifier** — if the exact argmin BEATS the pool past that
comparison's own MDE (`2.8016 * SE`, per L8 an operator and not the quoted constant), the search
half is NOT closed and a better sampler on this latent is worth building; **null** — the latent's own
MEAN, since an objective that cannot beat the average of the space it is searching has no ordering
skill on it at all.

**Two-target smoke (n=9, n=10, 512 and 1024 configurations, exhaustive) already shows the shape:**

| arm | basis | RMSD |
|---|---|---|
| latent_mean (the whole space) | mode set | 3.962 |
| **latent_argmin (EXACT optimum)** | single | **4.017** |
| latent_oracle | single, ORACLE | 1.757 |
| pool_argmin | single | 3.462 |

`argmin - space mean = +0.055`, and `argmin - ORACLE = +2.260`. On these two targets the **exact
optimum of the objective over the entire latent is no better than a random point in it**, while the
space **contains** structures 2.26 Å better. Two targets prove nothing; the full 75 will.

## L13 — WORKSTREAM D CLOSES D1: THE ENCODING LEVER IS NOT SUPPORTED

D replicated the Sprint-20 encoding result with `track_quality` OFF and 4 seeds. **Where step counts
are exactly equal, the effect is zero**, and on DIST/spsa that null is **properly powered**
(`-0.015 [-0.050,+0.022]`, own MDE **0.052**, below the programme's 0.084 bar). The decisive
Sprint-20 cell AMBc/spsa shrinks from `-0.649 [-1.205,-0.165] "sig"` to `-0.295 [-0.636,+0.022]`.
D's own falsifier did not fire.

**The argument that settles it is a symmetry argument, not a statistical one.** At matched budget,
adam_fd gives **-0.165 on LEG (embedding better) and +0.100 on DIST (embedding WORSE), both CIs
excluding zero**. A change of coordinates that provably does not change the physics **cannot help on
one energy and hurt on another**. A change in **step count** can — and the encodings differ 4.6 vs
9.6 steps at fixed budget, because the embedded parameterisation has twice the dimension.

**The durable output is a BUG, not a verdict.** In the embedded arm the gauge radius `||u||` grows
from 1.000 to ~2.5 (max 4.03) along a direction whose true gradient is **exactly zero** — `arctan2`
is scale-invariant, so every unit of radial motion is **injected finite-difference noise**, and
`adam_fd` makes it worse by normalising per coordinate by `sqrt(v)`, giving a pure-noise direction a
full-size step. The effective angular step then decays as `1/||u||`. **`_EmbField` is running a
different OPTIMISER, not a different set of COORDINATES.** Fix (project the update onto the tangent,
or renormalise `u` each iteration) before any future embedded arm is meaningful.

**Sprint 21 §5 item 3 is hereby moved from "determine whether this is robust and exploit it" to
NOT SUPPORTED.** D correctly declines the stronger claim: the AMBc panel's own MDE is ~0.53 Å and
**cannot exclude** an effect of the originally reported size, so the label there is **NOT MEASURED
with a 55%-smaller point estimate** — not "no effect exists". D also declines to claim its own
mechanism is proven, since the step-count response is CI-spanning and oppositely signed on the two
LEG arms. Both refusals are correct and are adopted as written.

---

## L14 — RESULT: THE EXHAUSTIVE LATENT ARGMIN DOES NOT BEAT A ZERO-EVALUATION POOL

`s21/results/latentfull.json`, **complete: true, 75/75 targets**, every one of the `2**n` latent
configurations evaluated at its von Mises mode. Latent sizes 512 / **median 4096** / 8192.

| arm | basis | mean | median |
|---|---|---|---|
| latent_mean — the ENTIRE space averaged | mode set | 4.511 | 4.402 |
| latent_argmin — exact optimum, squared functional | single | 3.396 | 3.403 |
| **latent_argmin_bayes — the EXACT optimum, DEPLOYED selector** | single | **3.281** | 3.294 |
| latent_draw — one stochastic draw per bitstring | single | 3.479 | 3.210 |
| latent_oracle — best of the SAME modes by true RMSD | single, ORACLE | **1.683** | 1.728 |
| pool_argmin — the shipped score's argmin over K=500 | single, window | 3.180 | 3.357 |
| **pool_argmin_rb — the same, REBUILT onto the latent's manifold** | single, rebuild | **3.257** | 3.330 |
| pool_tail_avg — the shipped top-75 average (incumbent) | avg of 75 | 2.765 | 2.678 |

### The falsifier did not fire: the exhaustive optimum TIES the pool

**Primary, on a matched FUNCTIONAL, a matched READOUT and a matched BASIS —
`latent_argmin_bayes - pool_argmin_rb = +0.025, SE 0.072, MDE 0.203, [-0.114,+0.168], 34W/41L`.**
A **dead-flat tie**. Not concentrated (mean-median `-0.047`, 19th percentile of a uniform-effect
null), folds mixed in sign (+0.05, +0.17, +0.03, -0.04, -0.07), the latent winning on **34 of 75**.

**The exact optimum of the objective over the ENTIRE latent is statistically indistinguishable from
the argmin of a K=500 retrieval pool that costs ZERO evaluations.** The pre-registered falsifier
required the exhaustive argmin to **BEAT** the pool past this comparison's own MDE. It does not
beat it; it ties it.

### THE OPERATOR LADDER — and it is the most important thing my lane produced tonight

The first version of this comparison reported **+0.215** and I was ready to write it up as "the
exhaustive argmin LOSES to the pool." It contained **three unstated operator differences**, found in
sequence — one by me before the run, two by Workstream D after it. Removing them one at a time:

| comparison | unstated operator removed | effect |
|---|---|---|
| `latent_argmin - pool_tail_avg` | *(caught by me pre-run)* argmin vs **average of 75** | +0.630 |
| `latent_argmin - pool_argmin` | argmin vs argmin; **squared vs Bayes** functional remains | **+0.215** |
| `latent_argmin_bayes - pool_argmin` | functional matched; **built chain vs window** basis remains | **+0.101** |
| `latent_argmin_bayes - pool_argmin_rb` | **fully matched** | **+0.025** |

**89% of the original effect was operator, not physics.** And the part that should be uncomfortable:
**all three differences pointed the same way — toward the conclusion I expected.** The averaging
operator flatters the pool; the squared functional handicaps the latent; the window basis flatters
the pool. That alignment is not a coincidence to be noted in passing. A pipeline assembled by
someone who already believes a conclusion will, by default, wire the operators so that they agree
with it, because each individual choice is made by reaching for whatever is nearest to hand in the
arm being defended. **The defence is not vigilance. It is the mechanical rule: state the operator on
BOTH sides of every paired comparison, and prefer the arm you did NOT build when they differ.**

> **SELF-CORRECTION, WITHIN THE HOUR, PROMPTED BY WORKSTREAM D's D7.** The first run of this file
> scored the **latent** with `s19/qb_lib.Obj.raw` (squared-distance) and the **pool** with
> `s12/instrument.shipped_score` (Bayes risk). **Two different functionals inside a single paired
> comparison, unstated** — the exact hazard D had just documented in §2, sitting in my own primary.
> It reported `+0.215 [+0.068,+0.370]`, which I would have written up as "the exhaustive argmin
> LOSES to the pool." **On the matched functional the effect is +0.101 and the CI covers zero.**
>
> **The mixture was worth +0.114 SE 0.069 on the latent — more than double the +0.050 D measured on
> pool windows**, same direction. That difference is itself worth recording: **the squared-vs-Bayes
> null is a property of the POOL's structure distribution and does not transfer unchanged to the
> latent's.** A lane invoking D7's interchangeability to skip a re-run must therefore check it is
> operating on pool-like structures.
>
> **SECOND CORRECTION, SAME HOUR, ALSO WORKSTREAM D.** A *second* unstated difference sat in the
> same primary: the `latent_*` arms are **BUILT CHAINS** (torsions through `build_ca_exact`) while
> `pool_argmin` was `W[order[0]]`, the retrieved **window's own coordinates**. `pool_argmin_rb` now
> puts the pool on the latent's own manifold — all 500 windows rebuilt from `u["PHI"]/u["PSI"]`,
> **scored AND measured on the rebuild**, so it is not half-converted. Cost on the pool argmin:
> `+0.076 SE 0.045`.
>
> **I did NOT apply Workstream A's +0.011 Å per-member rebuild cost as a correction term, and that
> refusal is the point.** A per-member mean does not price an **argmin**: rebuilding changes *which
> window wins*, and that effect has no reason to equal the average displacement. Measured here it is
> **7× A's per-member figure**. Correcting an argmin with a mean-shift constant is the same class of
> error as quoting an MDE as a constant (L8).
>
> **The corrections moved the conclusion from "loses" to "ties," which is the STRONGER form of the
> claim** — the search half is closed without needing the objective's optimum to be actively bad,
> and the concentration caveat the mixed comparison required disappears with it. I record all of it
> because every uncorrected version flattered my own conclusion, which is the direction an error is
> hardest to see.
>
> **DISPOSITION, in two parts that must not borrow strength from each other.** The **magnitude** is
> **NOT MEASURED** — `+0.025` against its own MDE of `0.203`. The **falsifier DID NOT FIRE** — it
> was directional, requiring the exhaustive argmin to *beat* the pool past the MDE, and a point
> estimate of `+0.025` with 34W/41L cannot fire it at any power. Low power limits our ability to
> detect an effect in the direction we did not need; it does not weaken the disposition.

### The decomposition, which is the actual finding — and it is not concentrated

All rows below use **ONE functional on BOTH sides**, so none of them was touched by the correction
above. The `_bayes` column is the deployed selector; the squared column is carried for continuity.

| quantity (Bayes readout) | mean | SE | folds |
|---|---|---|---|
| `latent_argmin - latent_ORACLE` (same 2**n modes) | **+1.598** | 0.140 | +1.08 … +2.29, **5/5 same sign** |
| `latent_ORACLE - incumbent` | **-1.082** | 0.141 | -0.853 … -1.487, **5/5 same sign** |
| `latent_argmin - latent_mean` (the null) | **-1.230** | 0.134 | — |
| `latent_draw - latent_argmin` (mode price, squared) | +0.083 | 0.091 | **null** |

(Squared readout, for continuity: `+1.713` / `-1.082` / `-1.116`. Same conclusion, same signs.)

**Read in order, this is the cleanest statement of the wall the programme has produced.**

1. **The objective is NOT noise.** It beats the mean of the very space it is searching by **1.116 Å**,
   ~2.8× that comparison's MDE. Any reading of the form "the distogram score is uninformative" is
   refuted here. **The search works.**
2. **The latent CONTAINS the answer.** Its ORACLE ceiling is **1.082 Å BETTER than the incumbent**,
   on **62 of 75 targets**, all five folds the same sign. Generation is **not** the wall. The von
   Mises basin decoder — a two-state-per-residue model — already spans structures at **1.683 Å** when
   the incumbent ships 2.765.
3. **And the objective, given EXHAUSTIVE access to every one of them, misses by 1.713 Å.** It picks a
   structure that beats the space's average handily and is still 1.7 Å from the best thing in the
   space it just enumerated in full. `latent_argmin` beats the incumbent on only **15 of 75** targets
   while `latent_oracle` beats it on **62**.
4. **The mode restriction is not doing the damage.** The stochastic arm — an actual von Mises draw per
   bitstring, the decoder as deployed — is `+0.083 [MDE 0.255]` against the modes, a null. So (1)-(3)
   are properties of the objective, not of the deterministic simplification used to enumerate.

### What this closes, and in what words

**`search-saturates-discrimination-binds` is no longer an inference from six instruments. On the 75
targets where the latent is enumerable, it is a THEOREM about this instrument**: there is nothing
left in the space for any sampler — quantum, classical, trained, or exhaustive — to find, because
the optimum has been found by exhaustion and it is not good enough. **No ansatz, no optimiser, no
budget, and no amount of quantum resource can recover the 1.598 Å**, because all of them are
different ways of arriving at an argmin that has already been computed exactly.

**The entire remaining leverage on these targets is the 1.598 Å between the objective's exact
optimum and the ORACLE over the identical set of structures — which is a DISCRIMINATION quantity and
is reachable only by changing the objective (or the readout), never by changing the search.**

### Scope, stated so it is not over-read

- **75 of 126 targets** (`n <= 13`). The other 51 have latents up to 65,536 and are **not** covered;
  the theorem is about the enumerable half.
- `latent_oracle` is **ORACLE** — it consumes the native and is a **ceiling, not an achievable arm**.
  It is exactly the quantity that separates generation from discrimination and it is labelled as such
  in every line above.
- The comparison is **argmin against argmin** by construction. The unmatched row
  (`latent_argmin - pool_tail_avg = +0.630`) is recorded only to show the **+0.415 Å the averaging
  operator would have contributed for free** had the primary been stated against the incumbent
  average — the trap named in L12, and the same one that produced two earlier retractions.
- A tempting monotone trend of `d` with `n` (+0.051 / +0.253 / +0.306 across n-bins) **does not
  survive inspection of levels**: the 11-12 bin is simply a harder bin on every arm including the
  pool (`pool_tail_avg` 3.014 vs 2.686 at n=13). **No size trend is claimed.**

---

## L15 — WORKSTREAM D CLOSES O2 WITH A QUANTITATIVE MECHANISM, AGAINST ITSELF

`d_enc_LEG-DIST_n60`, **n=60, 4 seeds, 3840 cells, complete**. Target-level pooled, the brief's unit:
`A (matched budget) -0.0061 [-0.0448,+0.0329]` **W/L 30/30, sign p = 1.0000**, CI half-width 0.039
against an own MDE of ~0.055 — **below the 0.084 bar, so a POWERED flat null.**

**The mechanism is quantitative, not a story.** For `adam`, the embedded arm at B=512 takes the same
number of steps as `theta` at B=256, so the objective's **own step-count response, measured in
theta, PREDICTS the encoding effect as minus itself**:

| cell | step response | predicted | observed | residual |
|---|---|---|---|---|
| LEG / adam_fd | +0.0912 [+.078,+.109] | -0.0912 | -0.1301 | -0.0389 [-0.122,+0.025] |
| DIST / adam_fd | -0.0865 [-.168,-.018] | +0.0865 | +0.0905 | **+0.0040 [-0.023,+0.027]** |

**Both residuals span zero, and DIST's is +0.004 against a predicted +0.087.** This is the part that
converts the argument: the sign-flip symmetry claim said *the coordinate explanation is impossible*;
this says *the step explanation is SUFFICIENT*, predicting the effect's **size** from an
independently measured quantity, and succeeding on the harder arm. **D1 = NOT SUPPORTED now rests on
a mechanism rather than on a null.**

### D closed its own hypothesis AGAINST ITSELF, and this corrects a memory of mine

D had proposed the mechanism as *"fewer steps helps here."* **It is not.** The step response is
**significantly signed in BOTH directions depending on the objective**:

    LEG   +0.091 [+0.078,+0.109]   more optimisation is WORSE
    DIST  -0.087 [-0.168,-0.018]   more optimisation is BETTER

**"Optimising harder hurts" is a property of the LEGACY objective, NOT of this instrument — and the
DEPLOYED distogram objective goes the other way.** This is the sharpest measurement the project has
of the conditional already carried in `search-saturates-discrimination-binds` ("the budget trap is a
property of BAD OBJECTIVES: searching harder hurts on a bad one, is neutral on a mediocre one, and
HELPS on a good one — always state that condition"). **The conditional was written because this trap
had been sprung before; D was quoting it without the condition, and so was I.** Both signs now have
CIs excluding zero, and they are being added to the memory so the conditional cannot be carried
without them again.

**A hypothesis closed against its proposer is the strongest form of closure available.** D's tally
for the sprint: **two of six** registered hypotheses died on their own pre-registered falsifiers
(R1 set-diversity, D7 distance functionals), and O2's mechanism is now a third correction to itself.

## L16 — WORKSTREAM B: THE LEVER IS THE READOUT HAMILTONIAN, NOT THE TRAINING ONE

Complete. **CVaR is the TRAINING objective and the readout is an argmin — nothing forces them to use
the same Hamiltonian.** Separating them on the **identical evaluated set, at ZERO extra budget**:

    changing the READOUT H     -0.697 A [-1.059,-0.352]   5/5 folds
    changing the TRAINING H    +0.056 A [-0.020,+0.129]   null

Replicated independently in two blocks (**-0.697 and -0.708, agreeing to 0.011 Å**), holding on all
17 arm-level tests and on `best_of_N` **where there is no training at all**. It does **not** beat the
incumbent — which already selects with the distogram — so it **prices a design choice** rather than
promoting an arm: *the Hamiltonian inside the training loop is nearly free to change; the one in the
readout is not.* Not free either on validity: the distogram readout carries more close Cα contacts,
significantly on one of three arms.

**Everything else in B is a null, and the nulls are powered where it matters.** The quantum-resource
half is **CLOSED**: measured χ reproduces the analytic `2^layers` exactly on all seven MPS variants,
a ring circuit reaches **χ = 48**, and **the entire ladder from χ=1 to χ=48 spans 0.37 Å** with
gradient norms *rising* with depth — so it is not a trainability null. Adam beats every alternative
on the CVaR loss at **0W/20L on both Hamiltonians, 5/5 folds**, yet **no optimiser beats `best_of_N`
on RMSD** and **training the circuit at all is worth +0.003 Å [-0.136,+0.144]**. QNG fails because
the **MPS Fisher is rank-deficient by construction** — a property of the parameterisation, not of the
landscape — and is significantly worse than no optimisation at all.

### L14 and L16 are the same statement arriving from opposite directions

L14 enumerates the entire latent and finds the **exact** argmin ties a zero-evaluation pool while the
ORACLE over the identical set is **1.60 Å better**: the search half is closed, and all remaining
leverage is discrimination. L16 then asks *where in the selector discrimination actually lives* and
answers: **in the READOUT Hamiltonian (-0.697 Å) and not in the training one (+0.056 Å, null)**,
with the ansatz (0.37 Å over a 48-fold bond-dimension range), the optimiser (no RMSD win over
`best_of_N`), and training itself (+0.003 Å) all null.

**Every knob that belongs to the SEARCH is null; the one knob that belongs to the READOUT is worth
0.7 Å.** That is the sprint's architecture answer, and neither lane could have produced it alone.

**B's declared limits, adopted:** F-Q3 on AMBER is **NOT MEASURED** (56 of 560 rows at hand-off) and
must not be generalised from the Legacy/distogram result, since the badly-conditioned landscape is
where QNG's case was strongest a priori; a 17× conditional quoted from a 15-target snapshot is
**RETRACTED**; F-A3 demanded 0.084 Å precision an n=20 design cannot deliver (achieved MDE 0.22-0.45),
leaving the **0.084-0.30 Å regime OPEN**; and the Block-Q battery is matched on learning rate, not
step norm (Adam moves 5.6× further than SGD).

---

## L17 — THE ENUMERATION IS COMPLETE AT n=126, AND THE SCOPE CAVEAT IS GONE

Workstream D ran the remaining 51 targets (`latentfull_hi.json`, **51/51, complete, flag verified**).
Merged, fold-clustered, per-comparison MDE:

    DISCRIMINATION GAP argmin - ORACLE (squared)  +1.8149 [+1.6129,+2.0594]   0W/122L   MDE 0.340
    DISCRIMINATION GAP argmin - ORACLE (Bayes)    +1.7202 [+1.4991,+1.9752]   0W/119L   MDE 0.326
    argmin - the space's own MEAN                 -1.3674 [-1.5333,-1.2100] 109W/17L    MDE 0.352
    latent ORACLE ceiling - INCUMBENT             -1.1766 [-1.2935,-1.0645] 104W/22L    MDE 0.325
    matched primary: argmin - pool_argmin_rb      +0.0851 [-0.0176,+0.1691]  53W/73L    MDE 0.197

**NOT ONE TARGET IN 126 has the objective's exhaustive argmin beat the latent's own ORACLE best —
0/122 squared, 0/119 Bayes.** L14's statement was about "the enumerable half"; **it is now about the
instrument.** The subset was defined by peptide LENGTH, which is exactly the selection a reader
should distrust, and it is no longer load-bearing.

**Disposition on the matched primary, the two parts kept separate:** falsifier **DID NOT FIRE**
(directional, 53W/73L, positive point estimate); magnitude **NOT MEASURED** (+0.085 against its own
MDE of 0.197).

### P3 answered: budget-versus-latent-size was NOT the operative variable

The 51 new targets are the **only** ones where the deployed VQE's budget was *smaller* than the
latent — the only ones where it was genuinely **searching** rather than resampling. There the
matched primary is `+0.1738 [-0.0079,+0.3715]`, 19W/32L: **exhaustive search still does not beat the
pool, and if anything does slightly worse than on the within-budget half (+0.0247).** Cross-range
difference `+0.1491 [-0.1505,+0.4549]`, NOT MEASURED. **The closure holds precisely where it was
least obvious.**

### D's P1 died, and its MECHANISM died separately — both recorded

D registered that the gap would be **larger** at n>=14 by min-of-N on the ORACLE arm, and declared
the confound in advance. Direction right and inside the guessed band (+0.252 squared, +0.301 Bayes),
but the cross-range CI spans zero — **NOT MEASURED**, the label D pre-committed to whichever way it
came out, now for two independent reasons.

**And the stated mechanism is contradicted by the data.** `latent_oracle` was predicted to *improve*
with n by min-of-N over a 128× larger set; it goes the **other way**, `corr(n, latent_oracle) =
+0.297`, per-length means rising 1.475 -> 2.025 across n=9..16. The arm that actually moves with
length is the space's **MEAN** (`corr = +0.677`, against +0.185 argmin, +0.144 pool). **Longer
peptides are simply harder, and that swamps min-of-N entirely.** Right direction, wrong reason,
unmeasurable anyway — **D's third hypothesis tonight whose mechanism was wrong even where its
direction survived** (after D1's "fewer steps helps" and R1's "diversity drives averaging").

## L18 — THE ANSWER IS IN THE OBJECTIVE'S TOP 5%, AND A TOP-512 READOUT CEILING IS 1.986 Å

`s21/latentrank.py`, n=75 (n<=13; the n>=14 half is running). **Independent recomputation** — its
`latent_oracle` reproduces L14's to `max |diff| = 0.00e+00` on all 75 shared targets.

L14/L17 say all remaining leverage is **discrimination**; L16 says the one non-null lever in the
selector is the **READOUT** Hamiltonian. Those meet on one question with an arithmetic answer:
**the best structure in the latent exists and the objective does not return it — where in the
objective's own ranking is it?** Null stated before the run: **no ordering skill puts it at 50.0.**

**It is at median percentile 5.62** (mean 20.27, min 0.0061, max 99.91):

    < 0.1 pct on 11/75    < 1 pct on 21/75    < 5 pct on 35/75    < 10 pct on 42/75    < 50 pct on 63/75

**The answer is not in the bulk. It is in the objective's top ~5% on the typical target** — so the
objective's failure is a failure of the **last few hundred places of its ordering**, not a failure to
order at all. The objective's own argmin sits at RMSD percentile **14.28** (median): what it returns
is a good structure, just not the best one.

### THE READOUT CEILING — AND WORKSTREAM D's MATCHED CONTROL, WHICH OVERTURNS MY FIRST READING

**I wrote the first version of this entry an hour ago and it was wrong.** It reported the ladder
below without a control and concluded that "a reranker over the objective's top-512 reaches 1.986 Å,
below the mission target." **Workstream D registered hazard D9 before those numbers existed:**

> *The top-M ceiling is a MINIMUM OVER M STRUCTURES. A curve that falls from M=1 to M=512 therefore
> falls PARTLY BECAUSE 512 DRAWS BEAT 1 DRAW, with no ordering skill involved.*

The matched control is the min over **M RANDOM** configs on the same axes; the **vertical gap is the
objective's ordering skill in exactly the currency a reranker spends.** Implemented (`rand{M}`, 32
seeded draws averaged), n=75, incumbent 2.765:

| M | top-M (objective-ordered) | rand-M (D9 control) | **GAP = ordering skill** | CI | W/L |
|---|---|---|---|---|---|
| **1** | 3.281 | 4.466 | **-1.184** | [-1.448,-0.918] | 63/12 |
| 8 | 3.068 | 3.270 | -0.202 | [-0.490,+0.074] | 34/41 |
| 64 | 2.498 | 2.459 | **+0.039** | [-0.204,+0.281] | 39/36 |
| 512 | 1.986 | **1.879** | **+0.107** | [-0.051,+0.300] | 47/28 |

**THE OBJECTIVE'S ORDERING SKILL IS CONCENTRATED ENTIRELY AT M=1 AND IS EXHAUSTED BY M≈8.** At M=1
it is enormous and unambiguous — the argmin is **1.18 Å better than a random configuration**,
[-1.448,-0.918], 63W/12L. By M=64 it is **zero**. At M=512 the objective-ordered set is **if
anything WORSE than 512 random configs** (+0.107, CI spanning zero, 47W/28L).

**So the entire fall of the ladder from 3.281 to 1.986 is min-of-M, not ordering.** My "1.986 Å in
512 evaluations" is true and meaningless: **512 RANDOM latent configs reach 1.879**, and beat the
incumbent by `-0.886 [-1.175,-0.614]` where the objective-ordered 512 manages `-0.779`. Below 2.0 Å
on 40/75 targets at random against 41/75 ordered.

**THE CONCLUSION INVERTS.** "Dig deeper into the objective's ranking" is **DEAD** — there is nothing
down there the ordering knows about. What the numbers actually support is different and much less
comfortable: **the LATENT is rich enough that 512 arbitrary draws from it contain a sub-2.0 Å
structure, and the objective contributes nothing to finding them beyond its top handful.** Any route
to that richness must be a **different selector**, not a deeper readout of this one.

**Reconciling this with the rank statistic above, since they look contradictory and are not.** The
ORACLE-best sits at median percentile **5.62** — genuinely better than the 50.0 null, so the ordering
*does* carry information about where the best structure is. But the **mean** is 20.27 with a max of
99.91: on a large minority of targets the best structure is in the bulk or worse. A **top-512 filter
on a 4096-config latent keeps 12.5%**, and the ordering is not reliable enough to keep the good
structure inside that window more often than a random 12.5% would. **Median rank skill and top-M
containment are different quantities, and only the second is what a reranker spends.**

**Workstream D's own registered prediction on D9 — "the gap is positive and substantial at every M"
— is REFUTED**: it collapses to zero by M=64 and reverses sign. **The hazard it raised was decisive
and its prediction about the hazard was wrong**, which is the fourth time tonight a D hypothesis has
died on its own falsifier while the instrument it motivated proved essential.

### What survives, stated at its real strength

    M=1, the DEPLOYED argmin, vs a random config      -1.184 [-1.448,-0.918]   63W/12L
    the latent's ORACLE ceiling vs the incumbent      -1.177 [-1.294,-1.064]  104W/22L  (n=126, L17)
    512 RANDOM latent draws + perfect selection       -0.886 [-1.175,-0.614]   52W/23L
    the objective's ordering beyond its top ~8        indistinguishable from random

**The selector has exactly one place where it demonstrably works — the very top of its own ordering
— and the mission's remaining 1.7 Å lives in a part of the space its ordering says nothing about.**

### The original ladder, kept for the record


Incumbent on these targets = **2.765** (pool top-75 coordinate average).

| M | best true RMSD in top-M | vs incumbent | W/L | folds |
|---|---|---|---|---|
| **1** (the DEPLOYED argmin) | 3.281 | +0.516 [+0.337,+0.718] | 19/56 | all + |
| 8 | 3.068 | +0.302 [+0.123,+0.486] | 25/50 | all + |
| **64** | 2.498 | **-0.267 [-0.423,-0.119]** | 49/26 | **5/5 same sign** |
| **512** | **1.986** | **-0.779 [-0.972,-0.598]** | 60/15 | **5/5 same sign** |
| all 2**n | 1.683 | -1.082 [-1.370,-0.817] | 62/13 | — |

**A reranker that sees only the objective's top-512 — 512 evaluations, ≤6% of the median 4096-config
latent — and picks perfectly inside it, reaches 1.986 Å: below the mission's 2.0 Å target, on
41 of 75 targets individually, 0.78 Å below the incumbent.** At **M=64** it already beats the
incumbent. **This is the first quantity in four sprints that puts <2.0 Å inside a 512-evaluation
budget.**


### WHAT THIS IS NOT — the caveat that survives the correction and is now the whole point

**Every top-M and rand-M number is an ORACLE**: "picks perfectly inside the M" means consuming the
native. These are **CEILINGS, not achievable arms.** The correction above does not soften that — it
**sharpens** it, because the ceiling that matters is now `rand-512 = 1.879`, reachable by *sampling*,
with **no objective involved at all**.

**The project has measured the achievable fraction on the POOL and it is the binding constraint**
(`in-band-ordering-is-per-target`: 0.986 within a target, **0.600 across targets**, 2.0 Å needs
0.638; `nothing-ranks-within-the-pool`; `in-band-signal-limited-not-sample-limited`). It has **not**
measured it on the LATENT, and D7 showed a pool-measured relationship failing to transfer to the
latent (L14), so those bounds carry over in neither direction.

**THE HIGHEST-VALUE REMAINING EXPERIMENT, restated after the correction:** measure native-free
selection skill over **512 RANDOM latent draws** — not over the objective's top-512, which the
control has just shown is not a privileged set. The enumeration already exists, so it is cheap.
*(Run as L20. It closes negatively.)*

> **SCOPE CORRECTION (Workstream D), so the two numbers are not read as contradicting.** `rand-512 =
> 1.879` is an **n<=13** figure. On the **full n=126 panel** the latent's 512-draw ORACLE is
> **2.280 against the pool's 1.721** — the latent's draws simply *contain* worse structures than the
> pool's windows (L20). Both statements are true of their own subsets. **"512 arbitrary latent draws
> contain a sub-2.0 Å structure" holds where it was measured and does NOT generalise to the full
> panel**, and it was never a comparative claim against the pool, which L20 settles against it.

---

## L19 — THE M=512 NULL IS NOT "NO SKILL". IT IS TWO OPPOSITE REGIMES CANCELLING

Workstream D noticed that L18's numbers do not reconcile: a selector whose ORACLE-best sits at
median percentile 5.62 should keep the answer inside a 6.2% window far more often than chance, which
predicts a **large** gap — yet the measured gap at M=512 is `+0.107`. D proposed latent density as
the resolution. **Measured at n=126, the resolution is different and sharper.**

    ORACLE-best CONTAINED in the objective's top-512 on 74/126 = 59%   (a random window gives 6%)

| regime | k | top-512 | rand-512 | GAP | latent ORACLE |
|---|---|---|---|---|---|
| **containment SUCCEEDS** | 74 | **1.568** | 1.963 | **+0.396** (ordering WINS) | 1.568 |
| **containment FAILS** | 52 | **3.533** | 2.682 | **-0.851** (ordering LOSES) | 2.305 |
| aggregate | 126 | 2.377 | 2.257 | +0.107 | 1.872 |

**The aggregate null is the average of a +0.396 Å win on 59% of targets and a -0.851 Å loss on 41%.**
The ordering is not weak — it is **bimodal**. When it contains the answer it contains it *perfectly*
(top-512 min **equals** the global latent ORACLE, 1.568 = 1.568). When it fails it fails
**catastrophically, and worse than random**, because the top-512 is then a *concentrated wrong
region* — 3.533 against a random window's 2.682, on targets whose own ORACLE is 2.305.

**This is the same shape as `prior-mae-prices-selected-rmsd`: the mean is set by the targets the
selector FAILS on, not by its typical behaviour.** My L18 reading — "ordering skill is exhausted by
M≈8" — was itself an **aggregation artefact**, and it is corrected here. The skill is not exhausted;
it is present and large on 59% of targets and inverted on the rest.

### THE CONSTRAINT ON WHAT THIS CAN CLAIM, and it is severe

**The split variable is ORACLE.** `contained` is defined by `oracle_rank_pct`, which consumes the
native, so this decomposition is **diagnostic of the mechanism and is NOT an achievable arm.** On the
contained subset `top-512 == latent ORACLE` holds *by construction of the split* — that cell is a
restatement of the selection, not evidence of skill.

> **CORRECTED IMMEDIATELY BY WORKSTREAM D, AND THE CORRECTION IS AGAINST THE CELL I HAD CALLED
> LOAD-BEARING.** I wrote that the failure cell was "a genuine, sign-definite statement about the
> objective and not an artefact of the split." **That is wrong.** Conditioning on `C=0` selects
> targets where the top-512 **provably missed**, and imposes **no corresponding condition on
> rand-512**. One arm is selected on its own failure and the other is not, so **THE SIGN IS FORCED**
> — a perfect selector would look bad on the subset defined by its own failures.
>
> **What is NOT forced is the MAGNITUDE, and that is where the finding actually lives.** On the
> failure cell the top-512 sits **1.228 Å** above the latent ORACLE while an arbitrary window sits
> **0.377 Å** above it — **3.26×** further. Conditioning forces "top-512 > ORACLE"; nothing forces it
> to be three times worse than a random draw. **That is the concentrated-wrong-region claim, and it
> survives.**

### The clean version, which needs no conditioning at all (D's formulation, n=126)

    ORACLE_top512 - latent_oracle   the ORDERED window's excess over the global best   +0.507 [+0.359,+0.673]
    rand512       - latent_oracle   the RANDOM  window's excess                        +0.388 [+0.317,+0.460]
    difference    == the D9 gap restated                                               +0.119 [-0.038,+0.290]

**No split, no selection effect, and it recovers the D9 gap exactly.** The bimodal decomposition is
therefore reported as **MECHANISM — "the aggregate null is the average of two regimes" — and NOT as
the evidence**, with the failure cell's sign marked forced wherever it appears.

### And it relocates the question precisely

**If a NATIVE-FREE signal can identify which targets the objective is containing on, the route is
live and worth ~0.4 Å; if it cannot, the aggregate +0.107 is what ships and the route is dead.**
That is a per-target confidence question, and it is exactly what
`in-band-ordering-is-per-target` already identified as the only remaining leverage
(*"the only leverage supplies the per-target SIGN at inference"*). **The programme has now arrived at
that same requirement from a third independent direction.**

**D's density hypothesis is partly right and worth keeping:** a random 512 lands within **0.388 Å**
of the global latent ORACLE on average (the objective-ordered 512 within 0.507 Å), so the latent
*is* dense at M=512 — dense enough that a random window is a strong baseline, which is why the
failure cell can be worse than random at all.

**D's phrasing, adopted:** *the objective is a GOOD FILTER and a USELESS RANKER* — 1.184 Å of skill
at M=1, nothing about the ordering of the next 500. **L19 refines it: it is a good filter on 59% of
targets and an actively harmful one on 41%, and nothing native-free currently tells them apart.**

---

## L20 — THE SOURCE IS NOT THE LEVER: THE RETRIEVAL POOL BEATS THE LATENT THROUGH THE SAME OPERATOR

`s21/latentsel.py`, **n=126, complete**, 512 latent draws per target. The pool and the latent are two
**sources** of candidate structures; the **identical deployed operator** (Bayes score filter -> top-75
coordinate average) runs on both, on the **same built-chain basis**, and nothing downstream varies.

| arm | basis | RMSD | median |
|---|---|---|---|
| lat_rand1 (zero-information null) | built chain | 4.982 | 4.570 |
| lat_argmin | built chain | 3.777 | 3.638 |
| lat_medoid | built chain | 4.437 | 4.072 |
| pool_argmin | built chain | 3.507 | 3.446 |
| lat_oracle (512 draws) | ORACLE | 2.280 | 2.285 |
| **pool_oracle** (K=500) | ORACLE | **1.721** | 1.767 |
| **lat_avg75** | POINT CLOUD | **3.454** | 3.091 |
| lat_avg75_rand (objective-free) | POINT CLOUD | 3.883 | 3.365 |
| **pool_avg75** | POINT CLOUD | **3.064** | 2.815 |
| ship_avg75 (the incumbent) | POINT CLOUD | 3.048 | 2.837 |

**PRIMARY: `lat_avg75 - pool_avg75 = +0.390, SE 0.086, MDE 0.242, [+0.224,+0.558], 47W/79L,
5/5 folds the same sign.** The pre-registered falsifier was *"if the latent beats the pool past the
MDE, the SOURCE is the lever and retrieval should be replaced."* **The opposite fired, past the MDE,
in every fold.** The latent source is **worse**.

**And it is worse on both counts at once.** `lat_oracle - pool_oracle = +0.558 [+0.413,+0.701]`,
35W/91L — at matched set size the pool is not only easier to average, it is **richer**. The latent's
512-draw ORACLE (2.280) still beats the incumbent (3.048), so L18's richness claim survives; **the
pool's ORACLE is simply better still.**

### Three controls, each doing its job

**D's basis price, requested so the deployment footnote is decomposable rather than trusted:**
`pool_avg75(rebuilt) - ship_avg75(window) = +0.016 [-0.004,+0.037]` — **negligible, CI spanning
zero.** So `lat_avg75 - ship_avg75 = +0.405` is a composite whose basis component is ~4%, and the
source swap is essentially the whole of it. **A's per-member +0.011 Å could not have told us this**
(the rebuild can change which member is the MEDOID, and the medoid sets the superposition frame);
measured on this operator it happens to be small, which is a fact and not an assumption.

**D's sensitivity, on targets where "512 draws" is not a draw at all** — at n=9 the space *is* 512,
so sampling without replacement is deterministic and hands those targets the exhaustive latent:

    all targets   k=126  +0.390     n>=10  k=117  +0.427     n>=12  k=93  +0.474

**Removing them makes the latent look WORSE, monotonically.** D's catch was real and pointed exactly
the direction it predicted — the mislabelled arm was flattering the latent.

### The number that looks like it contradicts D9, and does not

**`lat_avg75 - lat_avg75_rand = -0.429 [-0.561,-0.298], 91W/35L, 5/5 folds`** — the objective's
selection of *which* 75 to average is worth **0.43 Å**. D registered the prediction that these would
be **indistinguishable**, reasoning from D9 that the ordering is worthless past its top handful.
**D's prediction is REFUTED — and D9 is untouched, because the two measurements consume different
statistics of the same window.**

    D9's readout is a MIN over the window   -> dominated by the single BEST member
    this readout is an AVERAGE over the window -> dominated by the window's MEAN

`operator-consumes-set-mean` measured that exact decomposition on this instrument:
**`d_out = 1.16 * d_set_mean + 0.04 * d_set_best`, R² 0.89.** So an ordering can be useless for
*containing the best* while being very valuable for *raising the mean*, and this objective is
precisely that — consistent with `structural-objective-beats-the-energies`: **the distogram orders
the bulk.** The apparent contradiction is a readout mismatch, not a measurement error in either lane,
and it is the fourth time this sprint that two numbers reconciled only once the operator each one
consumes was named.

### What this closes

**"Replace retrieval with the generative latent" is CLOSED, negatively, at n=126 with matched
operator, matched basis, two nulls, a priced basis term and a pre-declared sensitivity.** The
pool wins by 0.390 Å through the shipped operator and by 0.558 Å at the ORACLE.

**It is also the first direct evidence for D's T2 mechanism**, which predicted this sign in advance:
*a candidate set generated from a SHARED parametric prior has coherent errors and averages badly,
while an INDEPENDENTLY retrieved set has incoherent errors and averages well — a property of the
GENERATOR, not of the scorer or the set size.* `lat_avg75_rand - pool_avg75 = +0.819` is the source
gap with the objective removed entirely, and it is **twice** the gap with the objective in play.
This is evidence *consistent with* T2, not a measurement of error coherence; the direct test is D's.

---

## L21 — THE CLOSING ARITHMETIC: ALL REMAINING LEVERAGE IS IN-POOL SELECTION, WORTH 1.33 Å

Two numbers from `latentsel.json`, n=126, both ORACLE and both scoring-only, settle where the
mission's remaining Å actually are.

### What the latent adds to the pool, under the most generous possible reading

    pool ORACLE alone        1.721
    latent ORACLE alone      2.280
    UNION of both sources    1.618        the latent adds -0.103 [-0.146,-0.067]

Take **both** sources, give a **perfect** selector free access to the union, and the generative
latent contributes **0.103 Å** — supplying the better structure on **35 of 126** targets. **That is
the ceiling on the latent's entire contribution, not its achieved value**, and L20 has already shown
the shipped operator extracts *negative* value from it (+0.390 Å worse than the pool alone).

**The generative-source branch is CLOSED**: not by an operator failing to extract, but because at the
ORACLE there is only 0.1 Å there to extract.

### Where the Å actually are

    the shipped incumbent                          3.048
    PERFECT selection inside the EXISTING K=500 pool 1.721

    pool ORACLE - incumbent   -1.327 [-1.546,-1.120]   SE 0.107   MDE 0.300   121 of 126 targets

**A perfect in-pool selector is worth 1.327 Å and would land the instrument at 1.721 Å — past the
mission's <2.0 Å target — using the retrieval pool that already ships, with no new source, no
larger budget, no deeper circuit and no different Hamiltonian.**

### The sprint's arithmetic, in one place

| avenue | measured value | status |
|---|---|---|
| search / optimiser / budget | exhaustive argmin ties the pool (+0.085) | **CLOSED by exhaustion**, L14/L17 |
| ansatz / bond dimension χ=1→48 | spans 0.37 Å, gradients rise with depth | **CLOSED**, L16 |
| training the circuit at all | +0.003 [-0.136,+0.144] | **null**, L16 |
| training Hamiltonian | +0.056 [-0.020,+0.129] | **null**, L16 |
| torsion encoding (θ vs sin/cos) | step-count confound; 0 at matched steps | **NOT SUPPORTED**, L13/L15 |
| generative latent as a SOURCE | -0.103 Å at the ORACLE of the union | **CLOSED**, this entry |
| **readout Hamiltonian** | **-0.697 [-1.059,-0.352]**, 5/5 folds | **the one non-null lever**, L16 |
| **in-pool selection** | **-1.327 [-1.546,-1.120]**, 121/126 | **where the Å are** |

**Every avenue this sprint was asked to test is measured at or below 0.1 Å except the two that are
selection.** The mission's 3.2 → 2.5 → 2.0 Å ladder does not require a better sampler, a better
generator, or more quantum resource. **It requires a better discriminator over candidates the
pipeline already retrieves** — and L19 locates the one place that discrimination demonstrably breaks
(the 41% of targets where the objective's window is a concentrated wrong region) without yet
supplying a native-free way to detect it.

---

## L22 — TWO DEFECTS THAT LOOK IDENTICAL AND ARE NOT: READ HAZARD vs WRITE-COLLISION HAZARD

Both found tonight, both in artefact I/O, and **the standard mitigation for the first does nothing
against the second.** Recorded together because they were nearly conflated.

**(a) THE READ HAZARD — non-atomic write.** `latentfull.py` used `json.dump(..., open(path,"w"))`.
Workstream D read a 45-row partial mid-write and computed a 5-target statistic on it before noticing
the file had been 75 rows a minute earlier. **Fix: tmp + `os.replace`**, which every other lane's
writer already did.

**(b) THE WRITE-COLLISION HAZARD — two writers, one path.** Workstream D restarted `d_lrank` twice to
add columns; `TaskStop` killed the first *wrapper* but not its Python child, and the second was never
killed. **Three processes wrote the same path simultaneously**, each clobbering the others' rows. D
noticed only because the row count went **backwards** between two reads.

> **ATOMICITY PROTECTS READERS FROM A HALF-WRITTEN FILE. IT DOES NOTHING AGAINST A SECOND WRITER.**
> Two *atomic* writers to one path produce a **perfectly well-formed file** containing whichever
> process wrote last — here **silently mixing two different configurations**, one build carrying the
> containment columns and two not. A file that passes every structural check and is a chimera.

**Fix: a lock or a per-run output path — not stronger writes.** Derive the artefact filename from the
configuration, as `d_enc` does (`d_enc_LEG-DIST_n60_spsa-adam_fd_B256-512.json` names its own arms),
or refuse to start when the existing file's `cfg_hash` differs from the running one. **`d_enc` got
this right because its tag encodes the arms; `d_lrank` did not, and that is the whole difference.**

**What actually caught it is worth more than the fix.** D's completion flag requires **every row to
carry every `nf_*` and `ORACLE_contained*` key**, not merely the right row count — so a mixed file
could not have passed it. **The flag was written for an unrelated reason and caught this one.** That
is the argument for making completion flags demand the **full key set**, not just `len(rows) ==
n_expected`: a row-count flag is exactly as vacuous here as `tailprice.py`'s was when it counted
skipped rows (Sprint 21 defect #10).

**Coordinator audit following D's report:** every one of my six artefacts
(`latentfull`, `latentfull_hi`, `latentrank`, `latentrank_hi`, `latentsel`, `poolgap`) verified —
each relaunch was strictly sequential with the previous process confirmed finished or dead, so no
collision occurred. **Verified, not assumed**, which is the only acceptable response to this class.

---

## L23 — THE 1.338 Å SELECTION GAP IS ENTIRELY UNCAPTURED. THE INCUMBENT IS THE BEST ARM THERE IS

`s21/poolgap.py`, **n=126, complete**. L21 priced the in-pool selection ceiling; a ceiling table with
no achieved column invites exactly one misreading — that the 1.338 Å is *available*. **It is not.**

    the shipped incumbent (score-filter -> top-75 coordinate average)   3.048
    PERFECT in-pool selection                                           1.711
    the gap                                                             1.338 A

| arm | RMSD | vs incumbent | % of gap captured |
|---|---|---|---|
| top-500 average (selection-free null) | 3.396 | +0.348 [+0.183,+0.520] | -26% |
| top-150 average | 3.072 | +0.023 [-0.035,+0.081] | **-2%, FLAT** |
| **top-75 average — THE INCUMBENT** | **3.048** | — | **0%** |
| top-20 average | 3.091 | +0.042 [-0.015,+0.102] | **-3%, FLAT** |
| top-5 average | 3.231 | +0.183 [+0.094,+0.277] | -14% |
| top-1 average = argmin readout | 3.454 | +0.406 [+0.273,+0.532] | -30% |
| consensus medoid, whole pool | 3.706 | +0.658 [+0.459,+0.854] | -49% |
| consensus medoid of the top-75 | 3.282 | +0.234 [+0.161,+0.302] | -18% |
| ORACLE ceiling | 1.711 | -1.338 [-1.558,-1.138], 123/126 | 100% |

**NOT ONE OF NINE NATIVE-FREE ARMS BEATS THE INCUMBENT**, and the two nearest (top-150 +0.023,
top-20 +0.042) are **flat** against it. So the result is not "we found a weak selector" — it is
**the incumbent is already at the top of everything this programme knows how to do.** The registered
hypothesis was "under 20%"; the registered falsifier ("more than half the gap with a CI excluding
zero") did not fire and could not have.

> **SCOPE LIMIT, IMPOSED BY WORKSTREAM D AND BINDING ON EVERY RESTATEMENT.** *"Captures 0% of the
> 1.338 Å gap"* invites the reading that **the gap is unreachable. That is not what was measured.**
> What was measured is that **NINE SPECIFIC native-free selectors** capture none of it — and **eight
> of the nine are variations on consensus/typicality**, a family `consensus-is-the-only-in-band-
> discriminator` already prices at **-0.172 Å**. **The honest headline is narrower and stronger:
> THE CONSENSUS FAMILY IS EXHAUSTED — not "in-pool selection is closed."** The 1.338 Å is an
> **unclaimed** ceiling, not a **proven-unreachable** one, and every restatement in the report says
> so.

Note the shape: **both the argmin (-30%) and the whole-pool medoid (-49%) are far worse than the
incumbent.** The argmin discards the averaging operator; the whole-pool medoid discards the score.
**The incumbent's value is the CONJUNCTION**, and each half alone loses more than either is worth.

### The m-ladder: a registered forward prediction, resolved NEGATIVE

`operator-consumes-set-mean` predicts **m\* SHRINKS 500 → 75 → 20 → 3-5 as the objective improves.**
That prediction was made from memory *before* the ladder ran, so it is a genuine forward test and not
best-of-K wearing a monotone hat.

**m\* has NOT moved from S8-11's 75.** Therefore, by the prediction's own logic, **the objective has
not improved** — which is precisely what L14/L17/L19/L20 independently establish.

> **WORKSTREAM D's CAVEAT 2 AND 3, both correct and both applied.**
>
> **(2) A monotone shape over a NESTED ladder is weak evidence.** m = 500/150/75/20/5/1 are nested
> subsets, so adjacent rungs are strongly correlated and RMSD-vs-m is smooth almost surely. *"An
> interior minimum exists"* is close to **vacuous** — it is what a smooth curve on a bounded interval
> does. The informative statement is not THAT an m\* exists but **WHERE**, and specifically
> **interior versus boundary**. Stated that way: **m\* is INTERIOR** (neither m=500 nor m=1), which
> is the non-vacuous content.
>
> **(3) Ties bite at the LADDER level, not only in the score.** A per-target `argmin` over a flat
> ladder picks arbitrarily and whichever rung numpy reaches first wins — `Tie-breaking leaks the
> pool order` one level up. Reported as a **tied set** (rungs within 1 SE of each target's own ladder
> minimum) instead: **median tied-set size 2.0 of 6 rungs**, and m=150/75/20 appear in 58/50/48 of
> 126 tied sets respectively. **Pairwise, m=150 and m=20 are both FLAT against m=75 (CIs spanning
> zero, +0.023 and +0.042).** So **"m\* = 75", "m\* = 20" and "m\* = 150" are the same finding, and
> the ladder cannot distinguish them.** D predicted this before the numbers existed.

### A correction to this file's own docstring

It asserted that **no tie handling is needed** because the score is continuous. **The count says
otherwise: 4,434 exact score ties across the 126 pools**, ~35 per 500-member pool. The claim was
wrong; what is *true* is narrower and had to be checked rather than assumed: **the stable sort breaks
those ties toward better BLOSUM retrieval rank, which is NATIVE-FREE, so no oracle order can leak** —
the failure mode that once manufactured a 1.386 Å winner here. The tie-break is **not neutral**, but
it is **not leaky**, and those are different properties. *Counting the ties is what turned an
assertion into a verified statement, and the assertion was false.*

---

## L24 — WORKSTREAM C: AMBER IS GENUINELY HARDER, THE CONTINUATION IS DEGENERATE IN RAW UNITS, AND PRECONDITIONING DOES NOT WORK

### (a) "Is AMBER genuinely a harder optimisation landscape?" — YES, on scale-free measures

Sprint 20 asked this and could not answer it without a normalisation that did not beg the question.
Measured at the medoid start on the shipped top-75 real rebuilds, n=30, using metrics that are
**invariant to the energy's units**:

| scale-free metric | Legacy | AMBER |
|---|---|---|
| condition number (median) | **6.26** | **7,756** |
| participation ratio | 0.497 | **0.067** |
| anisotropy | 4.18 | **18.19** |
| near-zero eigenvalue fraction | **0** | **0.481** |
| spectral skew | 0.190 | 0.797 |

**AMBER's Hessian is ~1,240× worse conditioned, its curvature is concentrated in 6.7% of directions
against Legacy's 50%, and half its spectrum is near-zero where Legacy has none.** The raw
distribution says the same thing from the other side: AMBER's pool skew is 8.13 and kurtosis 68.9
with **16% of members above 1e6 kcal/mol** and a max of **1.8e11** — the steric singularity, now
measured rather than asserted. **Legacy is a compactness model; AMBER is a landscape with a
half-flat, half-vertical spectrum.**

### (b) The mandated `H(λ) = (1−λ)H_L + λH_A` is DEGENERATE in raw units

    raw-unit crossover  lambda* = ||grad E_L|| / (||grad E_L|| + ||grad E_A||)
      median 1.989e-05    min 1.513e-09    max 5.582e-02    spread 36,888,025x

**A uniform λ grid in raw units is AMBER-dominated above λ ≈ 2e-5.** The entire interesting range of
the continuation lives in the first 0.002% of the interval, and every λ a person would naturally
choose (0.1, 0.25, 0.5, 0.75) is **the same Hamiltonian: AMBER.** Worse, λ\* varies by **3.7e7 across
targets**, so **no single λ grid is comparable between two targets.**

**This is not a criticism of the hypothesis — it is the reason it must be run in normalised units,
and it could not have been known without measuring the components separately first.** It is also a
clean instance of the sprint's rule: the operator (here, the units) decides the comparison.

### (c) AMBER preconditioning by Legacy — NOT SUPPORTED

    RMSD    rmsd_L - rmsd_0   +0.154 [-0.032,+0.431]   MDE 0.348   14W/16L    NOT MEASURED

The point estimate runs **against** preconditioning, and the mechanism it was supposed to exploit is
**absent in every diagnostic**:

| after Legacy preconditioning | plain start | preconditioned |
|---|---|---|
| AMBER energy reached (median) | 1.55e5 | **1.74e8** — 1,000× WORSE |
| condition number | 7,756 | **10,840** — worse |
| near-zero fraction | 0.481 | **0.503** — worse |
| ‖grad AMBER‖ | 2.03e7 | **2.38e10** — 1,000× worse |
| participation ratio | 0.0668 | 0.0718 |
| θ moved | — | **0.063 rad** |

**Legacy's optimum is not a good AMBER starting point.** Preconditioning moves the torsions by only
0.063 rad and lands in a region whose AMBER gradient is **three orders of magnitude larger**. The
landscape is not softened; it is entered at a steeper point. **The hypothesis is not merely null on
RMSD — its mechanism is refuted by its own diagnostics.**

### What this does to the sprint's Hamiltonian-design branch

Combined with **L16** (training-H `+0.056`, a null; readout-H `-0.697`), the picture is consistent:
`H(λ)` is a **training-side** construction operating in the dimension measured as null, and its
natural parameterisation is degenerate in the units anyone would use. **Both of the user's named
Hamiltonian-design hypotheses — continuation and preconditioning — are closed or blocked, and
neither closure came from a null on RMSD alone; each came with a mechanism.**
