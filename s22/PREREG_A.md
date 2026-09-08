# PREREG_A — Workstream A, Sprint 22 (CVaR-VQE selector)

Dated 2026-09-07. Frozen before any Sprint-22 number from this workstream exists. Addenda are
dated and appended, never edited in place (BRIEF §5).

Domain: the genuine CVaR-VQE selector — ansatz, CVaR, state preparation, measurement,
parameter optimisation, quantum landscape, candidate selection, quantum-vs-classical controls.
Three items, in the coordinator's stated priority order: (A1) candidate-state quantum encoding,
(A2) CVaR degeneracy / the tail as a distribution, (A3) multi-stage CVaR-VQE.

**What this workstream does NOT redo** (closed in s21, evidence cited): the ansatz ladder
(χ=1→48, best_of_N wins on AMBER, χ orders nothing — s21 L34/B2a); the optimiser comparison (no
optimiser beats best_of_N on RMSD — s21 L16/B3); the θ vs (sinθ,cosθ) encoding (radius≡scale,
NOT SUPPORTED — s21 L30/L35/L37); the training-Hamiltonian sweep (+0.056, null — s21 L16/B7).

**Exact theorem this workstream builds on** (s21 L3/L10, verified in `core/quantum.py`): for a
pool-restricted selector whose training objective and readout energy are the SAME H, with an
order-based readout (`vqe_bitstring`, `vqe_modal_bitstring`, `best_seen_bitstring`), CVaR-VQE and
argmin have the SAME optimal answer — CVaR changes the sampling distribution, not the selected
point. This is expected to reproduce EXACTLY in the candidate-basis encoding (A1) and is used as
a soundness gate, not a novel claim.

---

## A1 — CANDIDATE-STATE QUANTUM ENCODING

**Construction.** `H|i> = E_i|i>`, diagonal in the candidate basis. `n_qubits = ceil(log2(K))`
where `K` is the retrieval pool size (K=500 → 9 qubits, 512 basis states, 12 padding states with
no candidate). Candidate index → qubit computational-basis label is an ARBITRARY injective map
(the "gauge"); padding states get a fixed penalty energy above every real candidate's energy so
they are never selected. This makes CVaR-VQE genuinely and unambiguously the selector over a
finite pool: no basin latent, no continuous torsion, no ansatz standing in for anything but a
distribution over WHICH POOL MEMBER is chosen.

**Hypothesis.** (a) The EXACT optimum (global argmin of H over the K real candidates) is
label-invariant by construction (trivial, a soundness gate). (b) The ACHIEVED (trained,
finite-iteration) selector's outcome is NOT necessarily label-invariant, because the ansatz's
entangling structure (a CNOT chain / MPS of bond dimension χ) imposes a locality bias over
qubit-index Hamming distance that a candidate relabelling scrambles arbitrarily. This is
Sprint 18's encoding-gauge failure, restated for a new register, and this pre-registration exists
so it is tested rather than assumed away a second time.

**Primary endpoint.** Cα-RMSD of the argmin readout and of a tail-average readout (coordinate
average of the trained distribution's realised CVaR-tail-support, `s12.instrument.
coordinate_average`), against the native, ORACLE-scored post-hoc only.

**Falsifier for (a)** [soundness gate, not a claim]: the exact global-argmin candidate identified
by brute-force sort differs from the exact global-argmin candidate recovered by full enumeration
of the trained circuit's support, on ANY of the 16 targets, in ANY permutation. Reported as a
gate with a fired-count (BRIEF §3 rule 5 / §6 rule 7 — a gate that never fires is not evidence,
so the fired-count is printed even though it is expected to be zero).

**Falsifier for (b)**: across `P=8` random relabellings × `2` seeds each, per target, the achieved
outcome (exact CVaR value at convergence; argmin-readout RMSD; tail-average RMSD) does NOT differ
from the canonical-labelling seed-noise distribution (`4` seeds, identity permutation) by more
than that distribution's own spread (i.e. gauge variation is statistically indistinguishable from
seed variation). If falsified — gauge variation exceeds seed variation on a majority of
targets — the candidate encoding has a REAL gauge dependence in its achieved quality and that is
reported as the finding, not hidden.

**Null.** A permutation that maps candidates to qubit labels uniformly at random; the comparator
population is the canonical labelling's own seed-repeat spread, not a theoretical zero.

**Matched controls (Hard Requirements).** Random selection (random subset of size `⌈αK⌉`,
16 draws, averaged); greedy (= the exact classical score-sort, deterministic, IS the classical
control for the order-based readout here — search is not the bottleneck in a 512-state register,
so "greedy" and "exhaustive argmin" coincide and this is stated rather than hidden); simulated
annealing over candidate identity (matched to the VQE's own move class is not meaningful in a
diagonal register, so SA here is a discrete Metropolis walk over bit-flips, run for the SAME
number of unique-candidate evaluations as the VQE training uses non-cache-hit evaluations — since
that number is tiny (≤512) relative to any budget used elsewhere in this programme, this control
is expected to be near-vacuous and is reported as such, not inflated); best_of_N **from the
UNTRAINED circuit** (N draws from the ansatz at its random initial θ, before any training step —
never an initialisation mean, BRIEF §7).

**RMSD decomposition reported for every architecture** (BRIEF's binding requirement):
`R_pool` (oracle ceiling, best candidate in K) / `R_score_argmin` (classical exact sort — the
"search" stage) / `R_VQE_argmin` (trained circuit's argmin readout) / `R_VQE_tailavg` (trained
tail-support coordinate average) / `R_score_tailavg` (classical exact top-α coordinate average,
the incumbent-analogue ceiling for this readout) / `R_repair` (ideal-geometry projection of the
tail average, `s12.instrument.project`).

**Hamiltonians used.** `H_DIST` (native-free shipped Bayes-risk score) primary, per s21 C8f's
established prescription that distance alone beats every physics hybrid and physics energies
should not enter the selector. `H_Legacy` secondary, run in full, independently evaluable
(Pillar 2), cheap (no OpenMM). **`H_AMBER` is NOT run in A1/A2.** This is a declared scope
limit, not an oversight: AMBER is OpenMM-serialised and expensive, s21 B12/L29 already
established the one non-null AMBER role (training-side, argmin taken with something that ranks),
and A1/A2's question is about the ENCODING's gauge and expressivity, which is Hamiltonian-agnostic
by construction — testing it on the cheap, native-free H is sufficient to answer it. If A1/A2
motivate an AMBER cell, it is queued as an addendum and AMBER access is announced/released per
BRIEF §8 before it is taken.

**Ansatz arms.** `MPSAnsatz(n, layers=2, final_ry=True, entangler="cnot")` (deployed default,
χ=4) primary; `entangler="none"` (product state, zero correlation) as the expressivity control —
an arbitrary tail-membership subset generically requires correlation across qubits, so the
product ansatz is predicted to fail at tail-membership fidelity even where it does not fail at
the single-argmin readout (which any full-support distribution recovers via `best_seen_bitstring`
regardless of training quality, since a generic RY layer places non-zero probability on every
computational basis state as an OR-gate-free product of Bernoulli-type marginals — recorded
explicitly here as a mechanism prediction, not a result).

**Targets.** `s20.qb2_lib.subset(n=16)` — deterministic, fold-balanced, drawn from the 126-target
tuning instrument only. The sealed 60-target benchmark is not touched (hash-verified at the top
and bottom of this workstream's run).

**Budget.** Training is EXACT-gradient (`core.quantum.cvar_gradient_exact` via full enumeration
of the ≤512-state register — no sampling noise, no evaluation budget in the s19/s20 sense, since
scoring all K real candidates once is the entire cost of the Hamiltonian and is far cheaper than
any budget used elsewhere in this programme). This is stated as a property of the register size,
not a design choice that favours the encoding — the exact-argmin control is exactly as cheap.
Adam, 300 iterations, lr 0.12 (project default), seeds as specified per cell above.

**Promotion criterion.** A1 promotes to "the candidate encoding is a usable selector primitive"
only if (i) the exactness gate fires zero times (soundness), (ii) the gauge falsifier does NOT
fire (achieved quality is gauge-robust within seed noise) on the primary `cnot` ansatz, AND
(iii) `R_VQE_tailavg` beats `R_score_tailavg`'s matched-count random control (i.e. the trained
tail carries information the exact sort's own random-subset null does not) by more than that
comparison's own MDE. Failing any of the three is reported as a negative with the mechanism
named, not as "the encoding failed" generically.

### Rule 0 — six forks, named before the run

1. **Functional.** `H_DIST` (native-free Bayes risk, primary) vs `H_Legacy` (secondary, run in
   full). AMBER excluded (declared above). The alternative not taken: AMBER-in-training with a
   ranking readout (s21 B12) — not run here because it reopens an AMBER-role question s21 already
   answered; the alternative pushes toward "physics doesn't help", the direction this workstream
   already expects, which is exactly why it is excluded rather than run to pad the case.
2. **Basis.** Point-cloud RMSD throughout (candidates are literal retrieved windows, not a
   built/projected chain, except for `R_repair` which is explicitly a separate, labelled row).
   The alternative not taken — scoring everything post-repair — would inflate every arm's RMSD
   together and hide the encoding question inside a repair-operator effect (s21 L20's
   basis-price control found that effect ≈ 0.016 Å here, small but non-zero, so it is kept
   separate rather than assumed negligible).
3. **Readout.** Both argmin and tail-average are reported, never one alone (s21 L3's central
   lesson). The alternative not taken — reporting only the tail-average because it is what the
   incumbent uses — would hide the argmin soundness gate, which is the more informative number
   for A1's actual question (gauge invariance of the SELECTED SET, not of one summary statistic).
4. **Normalisation.** `H_DIST` is already in native units (Bayes risk, no rescaling); `H_Legacy`
   is used RAW (s21 established Legacy needs no normalisation for a single-Hamiltonian argmin/CVaR
   read, unlike the λ-continuation's cross-Hamiltonian blend, which does not appear here). The
   alternative not taken — standardising both to z-scores before comparison — is unnecessary
   because neither arm is blended with another Hamiltonian in this experiment; flagged so a
   later blend (A3) does not silently inherit an untested normalisation.
5. **Null.** The matched-count random subset draws from the SAME K=500 pool, not a uniform
   torus or a uniform categorical over all `2**n` basis states (which would spend most of its
   mass on padding states and manufacture an easy-looking win for any real selector, BRIEF §7's
   "zero-information controls must be plausible" rule).
6. **THE LABEL.** The candidate↔qubit-index bijection is itself the object under test in A1(b);
   it is generated by `s15.seed.stable_rng` from `(pdb, "qcand_perm", k)` so every permutation is
   reproducible and named, and the canonical (identity) labelling is never silently treated as
   "no label" — it is one specific label among `9!`-ish choices and its seed-noise spread is the
   comparator, not a zero.

---

## A2 — CVaR DEGENERACY AND THE TAIL AS A DISTRIBUTION

**Foundation, verified in source (s21 L3c, reproduced here as a gate).** `cvar_from_probs`
returns `mass`, the exact per-state probability the trained distribution spends inside the
α-tail (ordered by energy, not by probability) — this IS the CVaR-optimal face, exactly, not an
approximation, whenever the full `2**n`-state distribution is enumerable (true throughout A1/A2
at n=9). No sampling is needed to characterise "the face"; it is read off exactly.

**Hypothesis.** The face is not a point: multiple candidates share tail membership, and their
structural (RMSD, coordinate) spread is a measurable quantity the deployed argmin/best-seen
readout discards entirely. An additional structural term MAY break the degeneracy usefully; this
is tested, not assumed.

**Primary endpoint.** For each of 16 targets × each α ∈ {1.0, 0.5, 0.15, 0.05, annealed
(`core.quantum.alpha_schedule`, a0=0.5→a1=0.05), adaptive (α_t = clip(a_min, entropy_t /
max_entropy · a0, a0), a_min=0.05, a0=0.5)} × 4 seeds × 2 ansätze (`cnot`, `none`) × `H_DIST`
(primary) — report: tail size (states with `mass>0`) vs nominal `⌈αK⌉`; effective sample size
`α² / Σmass²`; RMSD spread (sd, range) of tail-support oracle RMSD; pairwise coordinate RMSD
spread within the tail-support set; energy spread within the tail; full-distribution entropy;
gradient variance (across the 4 seeds, at matched iteration); basin/typicality coverage (fraction
of tail-support candidates whose native-free typicality — mean distance to the other 499 pool
members — sits in the pool's own lower half, i.e. whether the tail also selects for
pool-typicality as a side effect, replicating or refuting s21's Legacy-compactness mechanism in
this new substrate); and FINAL Cα-RMSD of the argmin and tail-average readouts.

**The structural-term test.** A candidate-diversity or candidate-compactness penalty
`H' = H_DIST + β·H_struct` (`H_struct` = native-free pairwise-typicality energy, i.e. Legacy's own
compactness proxy re-used here as a STRUCTURAL term rather than a competing Hamiltonian — declared
so it is not confused with A1's `H_Legacy` functional-fork arm) is swept at `β ∈ {0, small,
moderate, large}` (fitted relative to `H_DIST`'s own pool sd, never in raw units — s21 C6's
degeneracy-in-raw-units lesson applied here pre-emptively) and its effect on tail RMSD spread and
on final readout RMSD is reported.

**Falsifier.** If NEITHER the α sweep NOR the structural term changes tail RMSD spread by more
than the per-comparison MDE on a majority of the 16 targets, the degeneracy exists but is
UNSTRUCTURED noise (i.e. members of the CVaR-tied face are not systematically better or worse
along any measured axis), and "breaking the degeneracy usefully" is closed negatively for this
Hamiltonian. This is registered as a legitimate, valuable negative per BRIEF's own instruction.

**Null.** A random-subset-of-matched-size control at each α (same construction as A1's null).

**Matched controls.** As A1 (random, greedy/exact-sort, best_of_N-untrained). Simulated annealing
is not repeated here (A1's SA result is inherited; re-running it under every α cell would not
answer a new question, since SA has no α parameter of its own).

**Budget.** Exact enumeration throughout (as A1); no sampling budget consumed. `β` grid: 4 values
× 16 targets × 4 seeds × 2 ansätze × 6 α settings = 3,072 trainings, each a ≤300-iteration
27-parameter Adam run on a 512-dim exact distribution — estimated well under an hour of CPU on
this box at ~95% utilisation, single-threaded per run, no OpenMM, no serialisation needed.

**Promotion criterion.** The structural term promotes only if it improves `R_VQE_tailavg` beyond
its own random-subset control's MDE on a MAJORITY of targets AND does not merely reproduce the
α-driven averaging-operator effect already priced in s21 (`operator-consumes-set-mean`) — checked
by holding tail SIZE fixed while varying β, so a shrinking-tail artefact cannot masquerade as a
structural effect.

### Rule 0 — six forks, named before the run

1. **Functional.** `H_DIST` alone for the α/degeneracy sweep (functional held fixed so the sweep
   is a pure readout/regularisation study, not a second Hamiltonian comparison — A1 already
   carries the Legacy functional fork).
2. **Basis.** Point-cloud, as A1, with `R_repair` broken out separately.
3. **Readout.** Argmin AND tail-average both reported at every α; α=1 is the degenerate
   full-mean case (equivalent to no CVaR restriction) and is INCLUDED, not dropped, because it is
   the natural zero-restriction control for the whole sweep.
4. **Normalisation.** `β` is expressed as a multiple of `H_DIST`'s own pool sd (never raw units),
   per s21 C6.
5. **Null.** Matched-count random subset, as A1.
6. **THE LABEL.** The α-schedule labels ("annealed", "adaptive") are two SPECIFIC, pre-declared
   functions of iteration/entropy, not tuned post hoc to whichever curve looks best; both formulas
   are fixed above and are not adjusted after seeing results.

---

## A3 — MULTI-STAGE CVaR-VQE

**Design.** Two stages, one circuit, warm-started parameters, on the candidate-basis encoding.
Stage 1 (purpose: enter a useful structural region) trains with `H_Legacy` — motivated by s21
L3's finding that Legacy's own tail is compact/pool-typical, i.e. a plausible "narrow the field to
plausible structures first" role, DISTINCT from s21 C7'/C7'' where staging was tested on the
CONTINUOUS torsion encoding and found harmful; this is a different substrate (discrete candidate
identity) and is tested rather than assumed to inherit that result. Stage 2 (purpose: discriminate
candidates) continues training the SAME parameters with `H_DIST`. The FINAL selection is genuine
CVaR-VQE on `H_DIST` (Pillar 1 preserved) with both argmin and tail-average readouts reported.

**Falsifier.** Staged (`Legacy→DIST`) does not beat single-stage `DIST`-only training (matched
TOTAL iteration budget, i.e. stage 1 + stage 2 iterations for staged ≤ single-stage's iteration
count, so staging is never given a free compute allowance) by more than the comparison's own MDE.
Given s21's C7'/C7'' precedent (staging harmful on the continuous encoding, mechanism: Legacy's
compactness bias fights the discriminating Hamiltonian's own optimum), the registered EXPECTATION
is that this falsifier FIRES (staging does not help) — stated in advance so a negative result is
not read as a surprise requiring a new post-hoc mechanism.

**Null / control.** Single-stage `DIST`-only training at matched total iterations, from the SAME
initial θ per seed (paired design). `best_of_N` from the untrained circuit, as A1/A2.

**Budget.** 16 targets × 4 seeds × {single-stage, staged} = 128 trainings, exact-gradient, cheap.

**Promotion criterion.** Promotes only if staged beats single-stage by more than its own MDE on a
majority of targets AND the argmin readout still ties the exact `H_DIST` argmin (Pillar 1's
"genuine CVaR-VQE is the selector, not a heuristic" requirement — a staged procedure that drifts
the argmin away from the discriminating Hamiltonian's own optimum would be reporting Legacy's
answer wearing a DIST label, which is exactly the mislabelling this programme's Rule 0 exists to
catch).

### Rule 0 — six forks, named before the run

1. **Functional.** Stage 1 = Legacy, Stage 2 = Dist, final = Dist. Alternative not taken: Dist
   first / Legacy second (would test "discriminate then narrow", a different mechanism claim
   with no motivating mechanism in this project's memory — not run, and named so no reader
   assumes it was tried and hidden).
2. **Basis.** Point-cloud, `R_repair` separate, as A1/A2.
3. **Readout.** Argmin and tail-average both reported; the PROMOTION criterion is pinned to the
   argmin (Pillar 1), not the readout most likely to flatter staging.
4. **Normalisation.** Both stages run their own Hamiltonian in its own native units; nothing is
   blended (no `H(λ)` continuation — s21 C6/C6' already closed that construction as degenerate in
   raw units and this workstream does not reopen it).
5. **Null.** Single-stage matched-iteration DIST-only, paired per seed and per initial θ.
6. **THE LABEL.** "Staged" vs "single-stage" is decided by iteration count alone, not by wall
   clock or by a subjective judgement of when stage 1 has "converged" — the total-iteration match
   is fixed above before any run exists.

---

## ADDENDUM 1 (dated 2026-09-07, appended, A2 body above left unedited)

**Pilot finding, before any A2 number was read for the record.** A smoke run of A1's own
machinery (target `1CS9`, α=0.15, 200 unregularised Adam steps on exact CVaR-over-p) collapsed to
`ESS=1.0` — a single-candidate delta, tail-support size 1 against a nominal tail of 77. This is
not a bug: it is exactly what `core.quantum.free_energy`'s own docstring states — *"for ANY alpha
the minimiser concentrates p on the lowest-energy basis states, so the readout collapses back to
the argmin"*. For continuous, generically-distinct candidate energies there is no numerical tie
at the argmin, so the theoretical CVaR-minimising face (s21 L3c) is a single point unless the
objective is regularised toward staying interior.

**Consequence for A2's design, declared before any A2 result exists.** The β-weighted
typicality "structural term" in the frozen A2 body above is REPLACED by the project's own
established mechanism for an interior optimum: the entropy-regularised free energy
`F_alpha,T = CVaR_alpha(p) − T·H(p)` (`core.quantum.free_energy`, lifted from its exact
`StatevectorCircuit` form to the exact `MPSAnsatz` score-function form used throughout this
workstream — both are exact diagonal-Hamiltonian gradients and the derivation is identical). `T`
replaces `β` as the swept axis; `T=0` reproduces the pure-CVaR collapse above and is reported as
the baseline case, not hidden. The falsifier, null, matched controls and promotion criterion in
the frozen A2 body are UNCHANGED — only the mechanism producing an interior distribution changed,
not what counts as success. `T` is expressed relative to the pool's own energy sd (never raw
units, same rule as `β` was given, s21 C6).

**Second consequence, methodological.** The α sweep is now read out at MULTIPLE iteration
checkpoints (5/20/50/100/200/300), not only at final convergence, because the collapse dynamics
— how fast ESS falls from the nominal tail size to 1 — is itself part of "the tail as a
distribution" and was not going to be visible from an endpoint alone. This is an ADDITION to the
frozen reporting list, not a substitution for any of it.

---

## ADDENDUM 2 (dated 2026-09-07, appended, A1 body above left unedited)

**Tie-breaking hazard, caught in the A1 pilot run before its numbers were read for the record.**
Target `5H1H`'s top-2 distogram scores are EXACTLY tied (gap 0.0 to float precision). The first
implementation of `argmin_readout`/`classical_exact_sort` used a bare `np.argmin`/`np.argsort`,
which breaks a tie by ARRAY POSITION — and array position for the trained circuit is bit-index
order, which the label permutation scrambles, while the classical control's array position is
candidate-canonical order, which it does not. Gate 1 (defined as exact-optimum agreement)
therefore fired on `5H1H` under a relabelling purely because the tie broke differently in the two
index spaces, not because of any genuine training or physics difference — the project's
documented "tie-breaking leaks the pool order" trap (`s21` memory), reproduced in this new
substrate before being caught.

**Fix, applied before any A1 or A2 number was read for the record.** Both readouts now return
the FULL set of candidates tied at the minimum energy (`s22/qcand_lib.tied_indices`) and score
ORACLE RMSD as the mean over that set. The tied set is a set of CANDIDATE indices, so it is
gauge-invariant by construction regardless of which bit-index order any given label induces.
Gate 1 is recomputed as SET equality, not index equality. **A first, tie-naive run of A1 is kept
under `s22/results/a1_gauge_PILOT_tiebug.json`** for the record rather than deleted, alongside
its log, exactly as the project's convention is to preserve rather than erase a caught error.

---

## Compute and reporting

AMBER/OpenMM is not used anywhere in A1–A3; no announcement/release cycle is needed for this
workstream's registered work. `s22/results/` artefacts are named from their config (target-count,
α-grid tag, ansatz, Hamiltonian) with atomic (`tmp` + `os.replace`) writes and a `_COMPLETE`
sidecar demanding the full key set, per BRIEF §3 rules 5–6. Every mean is reported with its SE and
its own per-comparison MDE (`2.8016 × SE`), never the retired 0.084 Å constant (BRIEF §3 rule 3).
ORACLE vs ACHIEVED is labelled at every appearance (BRIEF §3 rule 2). The 60-target benchmark is
hash-verified against `results/benchmark_manifest.json` before this workstream's first run and
after its last, unchanged.
