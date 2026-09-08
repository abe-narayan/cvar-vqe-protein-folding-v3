# SPRINT 13 — RESEARCH DOSSIER
## Torsion-constrained VQE/CVaR for peptide structure, and what a molecular energy model does to variational trainability

**The two questions put to the sprint.** Can a torsion-constrained VQE/CVaR architecture
reach < 2.0 Å mean CA-RMSD on unseen peptide structures? And does replacing a coarse-grained
molecular energy with an all-atom force field change the geometry and trainability of the
variational problem?

**The two answers.**

1. **No, and the reason is now measured precisely rather than inferred.** The discrete torsion
   space contains a 1.594 Å answer at ~24 live qubits — but the adversarial audit showed that
   **88% of what the library buys is generic Ramachandran statistics and only 12% is anything
   sequence-related**, and that 0.388 Å of the ceiling is the privileged oracle start. Neither
   a leave-fold-out sequence-only torsion predictor nor either physical energy can find the
   good configurations. Every native-free arm lands at 4.0–4.6 Å, worse than the 3.213 Å
   retrieval pipeline the architecture was meant to replace. **The objective-validity gate
   caught this before a single VQE was built on a broken objective, which is exactly what it
   exists for.**
2. **Yes, measurably — but not for the reason the sprint hypothesised, and the first version
   of the measurement was an artefact that two agents independently caught and corrected.**
   After proper conditioning, AMBER carries higher-order interaction structure than Legacy on
   79 of 79 cells. The chain from energy model → Pauli-weight spectrum → gradient variance is
   closed with **no free parameter** at a median measured/predicted ratio of 0.996–0.998. But
   cost-locality does **not** explain trainability in the regime this project can reach; *n*
   does. And the naive hypothesis — "AMBER is less local" — is a category error: both
   energies are full-register.

Seven agents ran concurrently. **benchmark60 was never read, dev24 was never spent, no tracked
file was modified.** Every arm that touches a native quantity is labelled ORACLE or DIAGNOSTIC.

| findings file | lines | agent |
|---|---|---|
| `s13/coord_FINDINGS.md` | 179 | coordinator — ceiling, cost model, objective validity, search |
| `s13/qarch_FINDINGS.md` | 1,063 | quantum encoding, locality, Hamiltonian, certified enumeration |
| `s13/geo_FINDINGS.md` | 926 | metric, gradients, CVaR, QNG, landscape |
| `s13/walsh_FINDINGS.md` | 716 | Pauli-weight spectra, per-term decomposition, encodings |
| `s13/lit_FINDINGS.md` | 521 | literature, novelty audit, 50-row method table |
| `s13/tors_FINDINGS.md` | 437 | sequence-only torsion prediction |
| `s13/adv_FINDINGS.md` | 470 | adversarial audit |

**Figures: `s13/figures/` (seven, with `s13/figures/README.md` naming the question each
answers; regenerate with `python -m s13.figures`).** The two that carry the sprint are
`fig01_money_gradient_prediction.png` — the spectrum-to-gradient chain, shown as both a
predicted-vs-measured scatter and, because a 48-decade log axis flatters any fit, the ratio
on a linear scale — and `fig02_money_accuracy_ladder.png` — every native-free arm against
what the space contains, with the helix control beside it.

---

# PART I — THE FOLDING RESULT

## 1. The representation is not the barrier

`s13/ceiling.py` → `s13/results/ceiling_report.json`. **ORACLE DIAGNOSTIC.** The architecture
optimises over `torsion_lib2.library_for(seq, k, exclude_seq=seq)`: k states per residue,
sequence-conditioned, target held out, `n_res·log2(k)` qubits. Native torsions snapped to
nearest states, then coordinate descent on true CA-RMSD.

    ideal-geometry rebuild of the NATIVE continuous torsions: 0.347 Å

| k | mean qubits | snap | **descent** | <2 Å | <1.5 Å | FAIL18 | other-108 |
|---|---|---|---|---|---|---|---|
| 4 | 25.9 | 2.860 | **1.594** | 0.73 | 0.48 | 2.078 | 1.513 |
| 8 | 38.9 | 2.218 | **1.184** | 0.91 | 0.70 | 1.539 | 1.125 |
| 16 | 51.8 | 1.844 | **0.876** | 0.95 | 0.89 | 1.296 | 0.806 |
| 32 | 64.8 | 1.375 | **0.634** | 1.00 | 0.97 | 0.783 | 0.609 |

At k=4, in ~26 qubits, the space contains a 1.594 Å answer against a retrieval pipeline at
3.213 Å. The 1.27 Å `snap`→`descent` gap is real interaction structure: the best per-residue
choice is not the best chain, because torsion error compounds.

**Two corrections from the adversarial audit, and they change how this table should be read.**

*The library is not meaningfully sequence-conditioned.* `library_for` defaults to
`mode="class"`: four k-means codebooks over the whole held-out corpus, one per residue class
(GENERAL / GLY / PRO / PRE_PRO). Measured over 126 targets at k=4, a target has a mean of
**2.38 distinct per-residue state tables** — exactly its number of distinct residue classes on
126/126 targets — and 84% of residues fall in GENERAL and therefore share one identical 4-state
table. **The library's entire sequence content is "which of four classes is this residue in".**
So the correct reading of the ceiling is that *a four-state Ramachandran quantisation* suffices
to express 1.594 Å structures — a stronger and more general statement than the sprint first
made, but a much weaker claim about the library being an asset.

*And the specific states barely matter.* The holdout removes a median of **1** entry from 787,
yet the resulting table differs from the no-holdout table by a median 68° matched-state shift —
that is k-means++ reseeding, not information removal. Re-running the entire ceiling at a
different seed moves it by **+0.034 Å [−0.044, +0.117], 61W/60L**.

*The audit then decomposed the ceiling properly, and this is the correction that matters most
in the whole sprint.* Running the identical descent in deliberately degraded spaces, all 126
targets:

| space the descent searches | ceiling |
|---|---|
| the real library | **1.594** |
| a **wrong target's** library | 1.729 (+0.135 [+0.043, +0.236]) |
| a class-blind codebook — the same 4 states at *every* residue | 1.729 (+0.135) |
| **uniform-random torsion states** | **2.698** |
| the real library from a **random** descent start rather than the oracle snap | **1.982** (+0.388 [+0.272, +0.506]) |

**Of the 1.104 Å the real library buys over an information-free space, 88% is generic
Ramachandran statistics and 12% is anything sequence-related.** And note the third row: an
information-free torsion space reaches 2.698 Å, which already beats the 3.213 Å retrieval
pipeline. The honest headline is therefore *ideal-geometry torsion parameterisation* is a
better representation than retrieval — **not** that this particular library is an asset.

## 2. No native-free objective finds it

`s13/coord_search.py`. Simulated annealing over the k=4 space at a **matched 5,000-evaluation
budget**, with uniform random sampling as the mandatory control (the literature survey found
Boulebnane et al., npj QI 2023, reporting QAOA performance on comparable problems "can be
matched by random sampling up to a small overhead").

| arm | mean | median | <2 Å | FAIL18 | other-108 |
|---|---|---|---|---|---|
| random sampling, 5,000 draws | 4.522 | 4.414 | 0.01 | 5.690 | 4.328 |
| **SA on Legacy** | **4.624** | 4.433 | 0.02 | 5.835 | 4.422 |
| SA on the 1-local torsion prior | 3.969 | 3.948 | 0.24 | 5.670 | 3.686 |
| SA on prior + Legacy | 3.980 | 4.030 | 0.25 | 5.763 | 3.683 |
| *ORACLE snap* | *2.860* | *2.297* | *0.42* | *4.262* | *2.627* |
| **ORACLE descent (ceiling)** | **1.594** | 1.548 | 0.73 | 2.078 | 1.513 |
| *shipped retrieval pipeline, for reference* | *3.213* | | | *6.019* | *2.745* |

| paired vs the random control | Δ | W/L | drop-top-10 |
|---|---|---|---|
| SA on the prior | −0.553 [−0.865, −0.238] | 75/51 | −0.243 |
| **SA on Legacy** | **+0.102 [−0.054, +0.255]** | 58/67 | +0.243 |

**Energy and structure are uncorrelated — and that is the correct statement, not "optimising
is worse".** SA drives the Legacy energy well below the native's (on 1A13, −21.6 against the
native's −16.0) and returns worse structures than the best of 5,000 random draws. But the
audit re-ran this with three stronger optimisers that reach 3.7–4.6 kcal/mol lower energy, and
**the sign reverses** (−0.044 to −0.060 Å); every CI spans zero and every drop-top-10 is
positive. Within a target, across six optimisers, ρ(energy reached, RMSD) = **−0.006**. So the
+0.102 Å is not a robust "optimising hurts" effect; what is robust is that **how well you
optimise this energy tells you nothing about the structure you get.** The stronger claim —
that pushing to the certified global optimum actively degrades the answer — is carried by the
exact-enumeration curve in §4, not by this arm. The audit separately confirmed the annealer
itself works: it beats random sampling's minimum energy on 41 of 42 targets by 3.1 kcal/mol.

**And the prior's apparent advantage is entirely a helix artefact** — the control that had to
be run, and was:

| | n | SA on prior | random | ORACLE ceiling |
|---|---|---|---|---|
| helical targets (>50% helix) | 42 | **2.066** | 4.229 | 1.067 |
| non-helical (<10% helix) | 70 | **5.218** | 4.827 | 1.918 |

ρ(helix fraction, SA-prior RMSD) = **−0.744**; ρ(helix fraction, prior's advantage over
random) = −0.646. The 1-local prior's modal state is the α basin, so "optimise the prior"
means "build an α-helix". **On non-helical targets its advantage vanishes** — the audit's
independent CA-only secondary-structure assigner puts the non-helical gap at +0.281
[−0.002, +0.568], so "worse than random" overstates it and "the advantage vanishes" is the
defensible form. It is a helix generator, not a predictor. The artefact itself survived the
audit: ρ = −0.709 under an assigner that never touches torsions, the two assigners agree 88%
at the 0.5 cut, and the "helical targets are just easier" confound explains under a fifth of
the effect (random's swing across strata is 0.60 Å against the prior's 3.15 Å).

**The audit then nailed the mechanism to a single constant, and in doing so weakened the
control.** Building *every* residue at φ = −63°, ψ = −42° — one ideal α-helix, no library, no
search, no budget, no sequence input of any kind:

| arm | all 126 | helical (n=42) | non-helical (n=70) |
|---|---|---|---|
| **ideal α-helix (a zero-information constant)** | **4.065** | **1.793** | 5.567 |
| SA on the torsion prior | 3.969 | 2.066 | 5.218 |
| random sampling on Legacy, 5,000 evaluations | 4.522 | 4.229 | 4.827 |

    constant alpha-helix - random sampling = -0.457 A [-0.822, -0.090], 66W/60L, drop-10 -0.116
    constant alpha-helix - SA on the prior = +0.096 A [-0.031, +0.220]

**A one-line constant reproduces the sprint's best native-free arm to within 0.10 Å with a CI
including zero, and beats the mandatory random-sampling control with a CI excluding zero.**
Two consequences, and the second is a methodological correction that outlives this sprint:
uniform random sampling in a k=4 torsion space is a **weak control**, so "beats random
sampling" is not evidence of anything in this representation; and **the constant α-helix is the
baseline every future arm here must clear.**

A third audit finding corrects the table's own labelling: the torsion prior is exactly 1-local,
so its global minimiser is the closed-form per-residue argmax. The arm labelled "SA on the
prior" returns that closed form to 1e-6 on **126/126 targets** — it performs no search at all,
and calling it a budget-matched search was wrong.

## 3. Neither energy ranks the native — measured three independent ways

The single most important gate in the sprint, and all three measurements agree.

**(a) Coordinator, 20 targets × 120 uniform configurations** (`s13/coord_ambershape.py`).
`rho_low10` is the rank correlation *inside the low-energy decile* — where an optimiser
actually lives, and the number that decides usability:

| variant | ρ over all configs | **ρ in the low-energy decile** | log₁₀ range | native percentile | descent percentile |
|---|---|---|---|---|---|
| Legacy | +0.198 | **+0.043** | 2.14 | 0.319 | 0.308 |
| AMBER raw single-point | −0.036 | −0.088 | **16.06** | 0.382 | 0.485 |
| AMBER capped at p90 | −0.036 | −0.088 | 10.68 | 0.375 | 0.480 |
| AMBER, monotone log compression | −0.036 | −0.088 | **1.81** | 0.382 | 0.485 |
| AMBER, electrostatics + solvation only | +0.015 | +0.027 | 2.42 | 0.399 | 0.483 |
| AMBER after restrained minimisation | +0.103 | — | 2.83 | — | — |

Three readings. **Both energies are garbage detectors, not rankers** — Legacy's low-decile ρ
is +0.043. **The native sits at the 32nd–40th percentile**, and the 1.6 Å descent structure is
scored worse than 31–48% of random garbage. And **the damage is in the rank order, not the
scale**: the log compression is monotone by construction, so ρ is identical to raw while the
dynamic range falls from 16.1 to 1.8 decades. *You cannot fix this by rescaling or softening.*

**(b) Architecture agent, 126 targets × 12,001 configurations.** Legacy ρ +0.190 (median
+0.232); near-native at the **32.5th percentile**; best-of-12,001 is 2.76 Å above the k=4
ceiling. **Remove the steric term and Legacy inverts to +0.395 Å worse than random** — its
entire skill is clash rejection. The 1-local prior beats it by −1.115 Å, 93W/33L.

**(c) Certified ground truth, 9 targets fully enumerated** (262,144 configurations each,
2.36 M structures scored, cached as `s13/results/qarch_enum_<PDB>.npz`). **Legacy's certified
global optimum is +0.139 Å worse than random sampling**, 2W/7L. AMBER's is +0.045 worse, with
in-band ρ = −0.015. No weighting of the two helps: the λ-sweep is monotone with its supremum
at Legacy alone.

## 4. "Optimise harder, get worse" — as a curve, and it is a warning about VQE

The architecture agent's sharpest result. On the enumerated spaces, RMSD versus evaluation
budget is **non-monotone**: 3.764 at 10 evaluations → 3.667 at 300 → **3.920 at the certified
global optimum**. The turn is governed by the *fraction of space explored*, not the evaluation
count: targets with n ≤ 11 (≥1.4% of the space seen) degrade **+0.267 Å** past their minimum
and end worse than random; targets with n ≥ 12 (≤0.36% seen) have not turned yet.

**A VQE at n = 12–16 with 10⁴ evaluations therefore sits on the improving part of a curve
whose limit is known to be bad. It will look like it is working.** Any optimiser comparison in
this problem must report objective quality and structural quality on separate axes.

## 5. Sequence-only torsion prediction is dead

`s13/tors_FINDINGS.md`. A leave-fold-out predictor over discrete Ramachandran cells with
per-residue densities and confidence, trained on the legal corpus.

Fifteen leave-fold-out arms. The best (a 324-cell Ramachandran categorical over sequence
context and residue properties, trained on peptides only):

| quantity | achieved | needed for 2.0 Å | needed for 1.5 Å |
|---|---|---|---|
| σ over all determined angles | **67.7°** | ≤ 17° | ≤ 12° |
| effective coverage at 12° quality | **0.23** | — | ≥ 0.90 |
| effective coverage at 20° quality | **0.36** | ≥ 0.90 | — |
| emitted, direct build | **3.770** | 2.409 | 1.535 |

Against the 3.213 Å incumbent that is **+0.558 [+0.30, +0.82], 36W/90L**, positive in every
fold, and the deficit *grows* when the best targets are dropped (+0.794 at drop-10, +0.915 at
drop-20). **A perfect oracle confidence gate still loses, at 3.478 (+0.265)** — the failure is
not selection, it is that only 23% of residues are ever accurate. The nulls stay dead on
emitted RMSD, including the deployable library-state prior (−0.246 [−0.424, −0.085]), but on
**angular error the composition null is not beaten** (−2.39° [−5.27, +0.81]).

A single number frames the whole route: on the ORACLE σ-surface, **the incumbent's 3.213 Å
corresponds to σ ≈ 29°.** A sequence-only builder has to beat 29° merely to tie the pipeline
it replaces, and the best arm here is at 67.7°.

The mechanism is specific and interpretable:

| arm | MAE φ | MAE ψ |
|---|---|---|
| sequence-blind corpus marginal | **36.4°** | 72.8° |
| residue-class prior, no learning | 36.9° | 68.9° |
| full sequence context + properties | **36.1°** | 62.4° |

**φ is predicted no better by a model that sees the entire sequence than by one that sees no
sequence at all** (36.1° vs 36.4°). The whole measurable sequence→torsion channel is 10.4° of
ψ. The literature agrees independently: SPOT-1D-LM reports φ/ψ MAE 22.2°/40.6° in
single-sequence mode on *full crystallised domains*, which maps onto Sprint 12's measured
surface at roughly 2.4–2.9 Å — and short flexible peptides are harder, not easier.

The residual decomposes cleanly: essentially none of the 3.770 → 1.535 Å gap is library
discreteness, and essentially all of it is angular error.

**And then the redirect, which is the most useful thing in that report.** Judged not as a
builder but as a *search space* — oracle descent over the predicted density, k=8, 39 qubits:

| space | oracle ceiling |
|---|---|
| trained arms | 1.485–1.656 |
| shuffled-label null | 1.442 |
| **sequence-blind corpus marginal** | **1.364** |
| the project's existing `torsion_lib2` | **1.184** |

**A learned density adds nothing as a search space and is worse at every k**, and the
inversion is sharp: the arm with the best emitted RMSD has the *worst* search space.
Per-residue native-cell recall at k=8 is only 0.34, so no per-residue argmax can work — but a
global objective is not bound by that. **The torsion architecture's bottleneck is the objective
that selects inside the space, not the prior that defines it**, which is the same conclusion
Part I reaches from the energy side.

Also confirmed here in a *third* independent model family: **peptide-only training beats
peptide-plus-fragment by +0.315 Å with nine times less data** (S7-2's distribution shift).

---

# PART II — THE TRAINABILITY RESULT

This half does not depend on the folding working. It needs only that the two energy models
are real, the encoding is real, and the measurements are exact.

## 6. An exact locality theorem, and its scope

`s13/qarch_FINDINGS.md` §1. **PROVEN**, 5,000+ (pair, variable) cells, agreement 1.0000, zero
counterexamples. The support of the CA–CA distance `d_ij` is **exactly the `j−i−1` residues
strictly between i and j**, contiguous; non-supporting variables move it by exactly 0.000 Å.
Support size is exactly `s−1` for separation `s`, with no scatter. Mechanism: `d_ij` is an
internal coordinate of the CA_i…CA_j sub-chain and the torsion-independent CA–CA virtual bond
removes the two end residues.

**Scope, and the audit strengthened it.** The architecture agent reported that all-atom
distances do *not* follow the CA rule — supports one residue wider, only 24% matching a clean
interval rule. Re-tested at machine precision (1e-9 rather than a 0.10 Å threshold), the
all-atom case is **also exact, under four atom-type rules**, each at agreement 1.0000:

    CA:  i < m < j          N:   i <= m < j
    C/O: i < m <= j         CB:  i <= m <= j     (two residues wider, not one)

Max |Δd| outside support **1.4e-13 Å** against a minimum inside support of 4.6e-2 Å. The CA
theorem also holds at chain termini, on glycine (74 targets) and proline (52), and under
*cis*-ω and arbitrary fixed ω. So AMBER populates the s=0 and s=1 shells that CA–CA leaves
exactly empty — not because its support is ragged, but because its atoms sit on either side of
the CA that the virtual-bond argument cancels.

**The naive hypothesis is refuted explicitly.** The union of supports over all pairs is the
whole chain, so *both* energies are full-register n-local. "AMBER is less local than Legacy"
is a category error, and the sprint records it as refuted rather than quietly replacing it.

## 7. The artefact, and the correction that makes the result defensible

**This is the sprint's most important methodological event.** The first Pauli-spectrum
measurement reported "AMBER's mean Pauli weight exceeds Legacy's on 14 of 14 cells". The
Pauli agent then found it was an artefact and the geometry agent independently confirmed it
on its own tables:

Over 141 fully enumerated tables, the **top-10 configurations (of 4,096) carry a median 99.6%
of raw AMBER's total Walsh variance**. A constant-plus-single-spike function has Walsh weight
spectrum exactly `Binomial(m, ½)` with mean weight `m/2` — measured 6.001 against a predicted
6.001, L¹ distance 0.0003. **A delta spike is maximally global for arithmetic reasons.** The
numbers were right; the interpretation was measuring steric clashes.

The geometry agent added two things the correction needed: some *Legacy* cells are spiked too
(top-10 share 0.878 on one cell), so it was never one-sided; and 99th-percentile winsorisation
is **not** sufficient conditioning — a winsorised AMBER table still carried 0.986 of its
variance in ten configurations. Its own earlier section was superseded and the raw sweep
preserved at `geo_pauli_v1_rawonly.json`.

## 8. The corrected result, and the money figure

Rank-preserving soft compression applied identically to both models (monotone, so ρ(E,RMSD),
argmin and native percentile are unchanged bit-for-bit):

| quantity | Legacy | AMBER | AMBER higher on |
|---|---|---|---|
| mean Pauli weight | 2.236 | **3.015** | **79 / 79 cells** |
| share of variance at weight ≤ 2 | 0.641 | 0.392 | — |
| ≥4-body Sobol share | 0.069 | **0.122** | 12 / 13 cells |
| cumulative 1+2-body share | 0.778 | 0.629 | 13 / 13 cells |

**Neither model is 2-local.** And the refinement that matters: the non-locality is in
**interaction order, not sequence range** — pair couplings are nearest-neighbour dominated for
both.

**The money figure — the chain closed with no free parameter.** Predicted gradient variance
from `Σ_S c_S² · Var_θ[∂⟨Z_S⟩/∂θ_i]`, against measured. The audit found that the sprint's two
headline ratios came from **two different predictors**, so only the exact per-Pauli-string form
is quoted here (95 cells in the current tables):

| model | median measured / predicted, **exact per-string** | same, coarse **weight-kernel** form |
|---|---|---|
| Legacy | **1.006** | 1.278 |
| AMBER (conditioned) | **1.001** | — |
| AMBER raw | 0.913 | 0.915 |

**The coarse weight-only summary is not good enough and is withdrawn** — it is reproducibly
1.49 on one cell. And raw AMBER's ratio is **not a converging estimator**: its across-θ
standard deviation *rises* with the number of θ samples (0.15 → 4.18, one cell spanning
0.219–12.227), so the conditioned arm is the one to quote everywhere. What survives is the
exact form on conditioned objectives, at 1.006 and 1.001, with the cross-covariance-dropping
approximation stated plainly and the ratio as its test. **This is, as far as the literature
survey found, the first time a real molecular force field's Pauli structure has been measured
and tied to variational gradient behaviour.**

## 9. What generates the structure, and a result nobody predicted

**Per-term decomposition.** AMBER's `nonbonded` term has covariance share **1.000 on all 26
component cells**; every other term contributes zero. Legacy's `steric` term reaches 0.955 at
m=18, rising with length. **Both models are a steric potential plus rounding error.** Strip
AMBER's nonbonded term and it becomes *more* local than Legacy on 6 of 7 ladder rows.

**The spectrum saturates.** Legacy's mean Pauli weight goes 3.89 → 4.10 → 4.17 → 4.27 across
m = 12…18 and the tail beyond weight 3 is flat at 0.60–0.64. **Effective interaction order is
about 3 residues, independent of peptide length** — the objective does not become more global
as the peptide grows. Nobody predicted this.

**Encoding dependence, quantified.** The residue-order spectrum is encoding-invariant to
1e-16, and the qubit-weight spread follows `b·2^(b−1)/(2^b−1)`: predicted 1.000/1.333/1.714 at
k=2/4/8, measured 1.000/1.31–1.40/1.81–1.94. The same Legacy energy has mean weight
**2.70 / 10.00 / 1.27** under binary / one-hot-mean / one-hot-penalty — a factor of 7.9. **A
one-hot Hamiltonian's Pauli spectrum is a free parameter of the implementer**, which is a
caution worth carrying into any published locality claim. (The sprint's "energy-ordered Gray
coding buys −16 to −21% at k=8" is corrected by the audit to **−10.9% against its own
random-relabelling null**, where binary sits 1.62σ *above* the null and Gray loses to the best
of 200 random relabellings on 6 of 6 cells. The effect is real at Stouffer z = −4.84 but far
smaller than first stated.)

## 10. Trainability: the honest negative

**Does the spectrum predict trainability? No.** The measured ansatz kernel `v(w)` is flat in
Pauli weight (`v(w_max)/v(1)` median 1.00) and decays in *n* as 2^(−0.47n)…2^(−0.86n).
Cerezo et al. 2021 requires local-2-design blocks at depth O(log n); a real-amplitude ring
hardware-efficient ansatz is not in that regime, **so the theorem licenses no prediction
here**. What the spectrum does predict is `Var[∂C] ∝ Var(E)` — so raw AMBER's gradient is
literally the gradient of the sampling probability of about ten clash bitstrings.

Exponential and polynomial decay both fit at R² 0.93–1.00 over n = 6–14 and **are not
discriminable**; the sprint says so. The real evidence is that the decay base falls
monotonically with depth onto the 2-design limit: 0.823 → 0.768 → 0.648 → 0.537 → **0.504** at
depth 8. A sharpened, testable prediction: the kernel is flat only at depth ≥ 3; **at depth 2
it falls 0.862 → 0.078 over n = 6 → 14**, so cost-locality would explain behaviour at depth 2
but does not at the depth this project uses.

## 11. Metric, QNG and CVaR — three clean results, two of them refutations

**The metric contains no Hamiltonian, and it is now a unit test.** `g(θ)` is bit-identical
across energy models at matched θ (**0.000e+00**, all seeds), while the two optimisers'
*endpoint* metrics differ by 118–189%. The geometry section is therefore correctly framed as
**trajectory selection**: the energy model does not reshape the manifold, it determines which
region's geometry the optimiser experiences. But the direction is **not stable** — AMBER
reaches the worse-conditioned endpoint on 3 of 6 runs, Legacy on 3, and the largest condition
number is Legacy's. Reported as a null.

**QNG is refuted, and not for the expected reason.** The metric is full rank at every θ, n and
depth; `g_ii = 0.2500` exactly; off-diagonal correlations 0.008–0.037 and *shrinking* with n;
exactly `I/4` at depth 1. The regulariser never bound. "AMBER's Euclidean parameterisation
matches its manifold badly" is **REFUTED** — the manifold is nearly Euclidean, so there is
nothing for a geometry-aware optimiser to fix.

**CVaR at small α is exactly a steric clash filter.** For α ≤ 0.25, `amber` and `amber_soft`
are *literally the same objective* — identical to every printed digit — while differing by
1.7e16× in gradient variance at α = 1. The entire dynamic-range pathology lives in the tail
that CVaR discards. And per the literature, CVaR is **not** a barren-plateau mitigation: the
bottleneck moves from trainability to *estimability*, so effective sample size `α·shots` is
recorded per arm and evaluation-budget parity is explicitly not shot parity.

**Physics honesty, and it is the discipline the brief demanded.** SPSA optimises the AMBER
objective best of five arms (percentile 0.088) and returns the **worst** structure (+0.333 Å
against random), with the RMSD rebuilt from the same bitstring. Better optimisation of a
misaligned objective produces worse physics, measured directly.

---

# PART III — NOVELTY, LIMITS, AND WHAT TO DO NEXT

## 12. The novelty audit, honestly

**Already published, and must be cited rather than claimed.** Per-residue φ/ψ binning with a
binary index and VQE+CVaR on 10–14-residue fragments: `arXiv:2510.06413` (Oct 2025), 127
qubits, 75 fragments, mean 4.89 Å. QFold (2021) did φ/ψ in 2–5 bits. Torsion plus a real force
field under a quantum optimiser: Marchand 2018 (D-Wave), Mato 2022 (HUBO). A lattice CVaR-VQE
reports **1.22–3.11 Å** on 6–20-residue peptides (Kannan et al., `arXiv:2510.15316`, IBM
Heron). **"Torsion encoding" and "VQE for peptides" are not novel and the sprint does not
claim them.**

**What the survey did not find anywhere: varying the ENERGY MODEL as the independent variable
under matched variational conditions.** The existing comparison literature varies the
*optimiser* with a Miyazawa–Jernigan potential fixed. Doga et al. (JCTC 2024) explicitly name
the coarse-grained-versus-all-atom comparison as an open gap and confirm that all-atom force
fields have only ever been used in classical post-processing in quantum protein structure
prediction — never as the objective inside the loop.

**So this is a trainability paper, not a folding paper, and it should say so.**

## 13. Limitations, stated before anyone else states them

- Register sizes 6–18 qubits, depths 1–8, noiseless exact simulation. Nothing here is an
  asymptotic or hardware claim, and no barren-plateau claim is made.
- The certified enumerations are n=9 at k=4 only. The budget-curve turn at n ≥ 12 is
  **untestable at any budget reached**.
- AMBER was never run on the full 126-target instrument, and the restrained-minimisation
  operating point was not run inside the variational loop (2¹⁴ × 5.9 s ≈ 27 h).
- **No AMBER gradient-magnitude claim survives normalisation.** Its interquartile range is
  9.7e+05 on one target and 21.9 on another; no normaliser is robust, and the sprint therefore
  makes no such claim.
- The geometry agent's CA-RMSDs are on 4–8-residue prefixes and are **not** comparable to the
  3.2 Å incumbent.
- Several sections rest on 1–3 cells or a single seed; those are flagged in place and no
  confidence intervals are quoted where none would be honest.
- The 1.594 Å ceiling is partly a property of a privileged oracle start: from a random start
  the identical descent reaches **1.982 Å**, so 0.388 Å of it is the start. It is an upper
  bound on what the space contains, reached by an oracle, and not a bound on any search.
- **The library's sequence content is worth ~0.13 Å of the ceiling.** A wrong target's library
  and a class-blind single-codebook space both cost only +0.135 Å, and uniform-random torsion
  states still reach 2.698 Å. Read the ceiling as a statement about ideal-geometry torsion
  parameterisation, not about this library.
- The "leakage-safe" holdout in `library_for` removes a median of one entry from 787 and does
  not meaningfully change the state tables. Leakage control at the library level is nominal;
  what protects these numbers is that no arm reads a native quantity at inference.
- `library_for` returns a bit-identical torsion table for several targets at short prefixes, so
  some geometry-only results have effective n = 3 rather than 6.

## 14. Corrections made during the sprint

The sprint corrected itself five times, and each correction is recorded rather than absorbed:

1. **My AMBER cost figure of 6 ms was wrong** — it re-evaluated the same state and hit the
   result memo. True cost is 12.5–28 ms per distinct configuration; `refine_coords` is 9.1 s.
   The error was live in the shared brief for an hour.
2. **Raw AMBER's Pauli spectrum was a delta-spike artefact**, caught by one agent and
   independently confirmed by another on its own tables; both superseded their own sections.
3. **99th-percentile winsorisation is insufficient conditioning** — the corrected instrument
   uses rank-preserving soft compression applied identically to both models.
4. **My single-target objective-validity reading (native at percentile 0.000) was optimistic**;
   on 20 targets the native sits at the 32nd–40th percentile, which two independent agents
   reproduced.
5. **`builder_for(...).energy()` minimises 50 steps** and is not a single point — a trap for
   anyone timing or interpreting it.
6. **The torsion library is a 4-class Ramachandran quantiser, not a sequence-conditioned
   library** (adversarial audit): 2.38 distinct state tables per target, 84% of residues in one
   class, and re-seeding moves the states by a median 68° while moving the ceiling by 0.034 Å.
7. **"SA on the torsion prior" performs no search** — the prior is 1-local, so its minimiser is
   the closed-form per-residue argmax, returned identically on 126/126 targets.
8. **Uniform random sampling is a weak control in this space** — a zero-information constant
   α-helix beats it by 0.457 Å [−0.822, −0.090]. Every future arm must clear the constant, not
   the random draw.
9. **The gradient-prediction chain is weakened in its coarse form.** The exact per-Pauli-string
   prediction reproduces at ratio 1.006 (Legacy) and 1.001 (conditioned AMBER); the coarser
   *weight-kernel* summary is 1.278 and is withdrawn, and raw AMBER's ratio does not converge.
   The dossier quotes the exact form on conditioned objectives only.
10. **A DEAD-QUBIT DEFECT in the encoding, found by the audit and previously unreported.**
   `core/project.py` leaves **three** backbone torsions per chain inert for the CA trace, not
   two: `phi[0]` and `psi[n−1]` are the known pair, and **`phi[n−1]` is a third**. Residue
   *n−1* is therefore entirely inert, so **log₂(k) qubits per chain are dead** and every
   `n_res·log₂ k` figure in this sprint is an overcount: "25.9 qubits" at k=4 is 23.9 live,
   and k=32 loses 5 of 64.8. Worse, `phi[n−1]` *does* place C, CB and O — so **Legacy and
   AMBER depend on a state that CA-RMSD provably cannot see.** Any future Hamiltonian should
   drop that variable from the register, and any qubit count quoted from this sprint should
   carry the −log₂(k) correction.
11. **"Energy-ordered Gray coding buys −16 to −21%"** is against the wrong baseline; against
   its own random-relabelling null it is −10.9%.
12. **"On non-helical targets the prior is worse than random"** overstates it; under an
   independent assigner the gap is +0.281 [−0.002, +0.568], so the advantage *vanishes*
   rather than reversing.
13. **My own pre-registered prediction that a sequence-only route would land at 2.4–2.9 Å is
   REFUTED** — measured 3.770 Å. The arithmetic was right and the premise was wrong: I read
   SPOT-1D-LM's MSA-mode row, but the deployable setting is single-sequence (σ ≈ 41–44°, not
   20–30°), and this model is a further 1.4–1.6× worse than SPOT-1D-Single on a chain length
   no published method has been evaluated on. The agent's own pre-registration (σ 45–70°,
   build 3.5–6.0 Å worse) was the accurate one.
14. **My coverage guidance was wrong, and wrong in the favourable direction.** I told the
   torsion agent to model missingness as *clustered* rather than uniform, on the literature's
   report that TALOS-N's gaps cluster at termini and in flexible regions, and predicted that
   uniform dropout would flatter the method. Measured: interior-block clustering is worse than
   uniform by only 0.02–0.13 Å, but **terminal dropout — the literature's actual pattern — is
   0.40–0.50 Å BETTER than uniform** (σ=12°, 90% coverage: 1.667 vs 2.168). **Uniform dropout
   understates a TALOS-N-style channel.** This makes the chemical-shift route more attractive
   than the sprint previously recorded, not less.

An independent leak was also found and handled: four targets whose exact sequence is a
contiguous substring of a training peptide (6B9K, 1CEK, 2FBU, 2P5H) pass the project's
longer-normalised 0.6 identity threshold. This reproduces the coordinator's Sprint 12 finding
by a completely different route, and the torsion agent dropped the containing parents, so its
numbers are strictly cleaner than the project's own standard.

## 15. The architecture the evidence supports

Not the one the sprint set out to build. The evidence says:

- **Keep the torsion representation, for the right reason.** It is the best structural
  representation the project has measured — but the asset is *ideal-geometry torsion
  parameterisation itself*, not the library: an information-free random-state space already
  reaches 2.698 Å against the retrieval pipeline's 3.213 Å, and 88% of what the real library
  adds is generic Ramachandran statistics. Drop the dead residue-(n−1) qubits from any
  register built on it.
- **Do not put Legacy or AMBER in the objective.** Both are clash detectors with ρ ≈ 0 inside
  the low-energy decile; the certified global optimum of each is worse than random sampling.
  Their honest role is a **validity filter** — reject clashing configurations — not a ranker.
- **The missing component is a many-body objective that does not exist.** The 1-local prior's
  argmin is 1.921 Å and the space's best is 0.969 Å on the enumerated targets. That 0.95 Å is
  the prize, and closing it needs an objective with real 2-and-3-body structure aligned with
  nativeness. Nothing in the project or the literature currently supplies one.
- **VQE/CVaR remains genuine and correctly placed**, but on present evidence it would be
  optimising an objective that does not rank the native, and the budget curve shows it would
  *appear* to work at n ≥ 12 while doing so. The honest role is as the optimiser in the
  trainability study, not as the accuracy mechanism.

## 16. Highest-value next steps

1. **Write up the trainability half.** It is complete, corrected, cross-validated, and novel in
   the one dimension the literature leaves open. It does not need the folding to work.
2. **Build the many-body objective** the enumerated ground truth now makes trainable: 2.36 M
   labelled structures with RMSD, 11 Legacy terms, the prior, and AMBER components are cached
   at `s13/results/qarch_enum_<PDB>.npz`. Learning a 2-and-3-body objective *against RMSD* on
   that data, leave-fold-out, is the first thing in this project with a plausible path to the
   0.95 Å prize — and it is a supervised problem, not a physics one.
3. **Measure BMRB chemical-shift coverage** on these exact targets, from shifts alone, with the
   pre-registered kill threshold of 75%. Sprint 12 priced the channel at 1.486 Å at full
   coverage, and this sprint corrected the coverage model in the route's favour: because
   TALOS-N's gaps fall at **termini**, and terminal dropout costs 0.40–0.50 Å *less* than the
   uniform dropout the earlier estimate assumed, the channel is better than previously
   recorded at any given coverage. The kill threshold should be re-derived against the
   terminal-dropout curve rather than the uniform one before the measurement is run.
4. **Replicate the random-sampling control, and re-baseline it.** The one cell where
   variational arms beat random is a single cell × 5 seeds and needs replication, precisely
   because it contradicts Boulebnane. But re-run it against the **constant α-helix**, not
   against uniform random sampling — the audit showed the random control is beaten by a
   zero-information constant, so clearing it demonstrates nothing.
5. **Do not** run a VQE on Legacy or AMBER for accuracy. The budget curve guarantees it will
   look like it is working.
