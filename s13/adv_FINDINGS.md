# Sprint 13 — ADVERSARIAL AUDIT

My job was to destroy this sprint's conclusions. Everything below is a measurement I ran, not
an objection I raised. Code: `s13/adv_instrument.py`, `s13/adv_ceiling.py`, `s13/adv_prior.py`,
`s13/adv_sa.py`, `s13/adv_locality.py`, `s13/adv_pauli.py`, `s13/adv_pauli2.py`.
Results: `s13/results/adv_*.json`. benchmark60 never read; dev24 never run; no tracked file
modified; every arm that reads a native quantity is labelled ORACLE.

**Headline: the two claims I most expected to break — the exact locality theorem and the
exactness of the Pauli/energy-table machinery — survived every attack I could construct,
including attacks their authors had not run. The representation-ceiling claim survives as a
statement but its *interpretation* is substantially weakened. Two instrument defects were
found, neither fatal, one of which changes a qubit count.**

---

## SEVERITY 1 — the shared instrument

### A. Coordinates: the energies see EXACTLY the coordinates the RMSD is measured on

**Status: SURVIVED (my strongest attempt).** `s13/adv_instrument.py` → `results/adv_instrument.json`.

| check | result |
|---|---|
| `I.build_ca(phi,psi)` vs `core.geometry.build_backbone(phi,psi)["CA"]` | max abs deviation **0.0** (15 targets) |
| `I.build_ca` vs `rep.build_coords(rep.bitstring_from_states(s))["CA"]` | max abs deviation **0.0** |
| `I.build_ca` vs `core.amber.single_point(...)["ca"]` | max CA-RMSD **1.9e-07** (Kabsch noise on identical arrays) |
| `rep.state_indices(rep.bitstring_from_states(s)) == s` | True, 15/15 |
| `PerResidueTorsion._phi/_psi` vs `library_for(...)[:,:,0/1]` | max abs error **0.0** |

There is no coordinate mismatch anywhere in the chain `states → bitstring → build_coords →
core.amber / core.energy` versus `states → PHI/PSI → I.build_ca → CA-RMSD`. This was the
defect that would have contaminated everything at once, and it is not there.

### B. DEFECT — three backbone torsions per chain are inert for the CA trace, not two

**Status: DEFECT CONFIRMED AND EXTENDED. Priced below.**

`s13/tors_FINDINGS.md` §0.2 reports two placeholder angles per chain (`phi[0] = -60°`,
`psi[n-1] = -45°`, constant across all 787 `peptide_db` entries — independently confirmed
here: 1 unique value each). That is a *database* fact. The *builder* fact is worse. Perturbing
each angle by 1 rad and measuring the maximum CA coordinate change over 20 targets:

| angle | max |Δ CA| |
|---|---|
| `phi[0]` | **0.000 Å** |
| `psi[n-1]` | **0.000 Å** |
| **`phi[n-1]`** | **0.000 Å** ← not reported anywhere in the sprint |
| `psi[0]` | 22.08 Å |
| `phi[mid]`, `psi[mid]` | 15.67 / 13.33 Å |

`core.project.build_ca_exact` (lines 426–433) places `C[:, i+1]` from `phi[:, i+1]`, and
`C[n-1]` is never read again, so **`phi[n-1]` is a live database angle that the CA builder
discards.** Consequences, in order of importance:

1. **The qubit counts are overstated.** Residue `n-1` contributes *nothing* to the CA trace,
   so its `log2(k)` qubits are inert for every CA-based objective and for CA-RMSD.
   At k=4 the sprint's "~26 qubits" (mean 25.9) is **23.9 live qubits**; 1A13's "28 qubits"
   is 26. At k=32 it is 5 dead qubits of 64.8. Every "n_res·log2(k) qubits" figure in
   `coord_FINDINGS.md` C0/C3 and `qarch_FINDINGS.md` should carry a `−log2(k)` correction.
2. **A genuine asymmetry between the two energy models and the metric.** `phi[n-1]` places
   `C[n-1]`, hence `CB[n-1]` and `O[n-1]`, so **Legacy and AMBER both depend on state[n-1]
   while CA-RMSD provably cannot.** An optimiser on either energy spends budget on a variable
   that cannot change the reported answer. This is small but it is a real free parameter in
   every search arm.
3. It does **not** bias the ceiling. `ceiling.snap` selects terminal states using the dead
   angles; a snap restricted to live angles gives a mean of 3.0654 vs 3.0654 Å over 40 targets
   (Δ = −1.9e-15) even though it picks a different terminal state on 21/40 targets. No effect.

### C. DEFECT — the "sequence-conditioned" torsion library is a 4-class Ramachandran quantiser

**Status: DEFECT CONFIRMED. This is the finding that most changes how C0 should be read.**

`torsion_lib2.library_for(seq, k, exclude_seq)` defaults to `mode="class"`, which is
`_class_library`: four k-means codebooks over the *whole* held-out corpus, one per class
(GENERAL / GLY / PRO / PRE_PRO). Measured over all 126 targets at k=4:

    mean residues per target                         12.96
    mean DISTINCT per-residue state tables            2.38
    equals the number of distinct residue CLASSES    126 / 126 targets
    mean fraction of residues in class GENERAL       0.842

So 84 % of residues in a target share one identical 4-state table, and the library's total
sequence content is "which of 4 classes is this residue in". Calling it
*sequence-conditioned* overstates it by a lot; `torsion_lib2`'s own docstring describes the
per-residue context form, which is **not** the default and is not what any arm in this sprint
used. Quantified consequence in §2 below.

### D. The leakage holdout is nominal, and the library is k-means-unstable

**Status: WEAKENED, but harmless.** `library_for(seq, k, seq)` removes a **median of 1**
database entry from 787, yet the resulting k=4 table differs from the no-holdout table by a
**median 68.0° / max 180.9°** matched-state shift. That is not information removal — it is
k-means++ reseeding (the seed draw is `rng.choice(m, p=w/w.sum())` and `m` changes with the
pool size). Two consequences:

* "leakage-safe" here is a statement about one removed peptide, not about the states;
* the states themselves are a *locally* optimal clustering that jumps between runs.

**It does not matter**, and that is itself the point: re-running the whole ceiling at
`seed=7` moves it by **+0.034 Å [−0.044, +0.117], 61W/60L** (§2). A library whose states can
move 68° without moving the ceiling is a library whose specific states carry little.

The containment leak `tors_FINDINGS.md` §1 reports (4 targets whose exact sequence is a
substring of a training peptide, passing the 0.6 identity threshold) therefore also cannot
move the ceiling: one peptide of 787 does not move a class-level k-means codebook.

---

## SEVERITY 1 — the representation ceiling of 1.594 Å (claim 2)

`s13/adv_ceiling.py` → `results/adv_ceiling_k4.json`, `adv_ceiling_k4_report.json`.
**126 targets, k=4, identical descent, identical start rule, identical sweep cap.**
My reimplementation reproduces `s13/ceiling.py`'s 1.594 Å exactly, so the comparison is clean.

| space | what it knows | snap | **descent** | <2 Å | Δ vs real |
|---|---|---|---|---|---|
| `real` = `library_for(seq,4,seq)` | class + corpus Ramachandran | 2.860 | **1.594** | 0.73 | — |
| `real_seed7` | same, different k-means seed | 2.845 | 1.629 | 0.67 | +0.034 [−0.044,+0.117] |
| `wrongseq` | **a different target's library** | 2.944 | 1.729 | 0.65 | **+0.135 [+0.043,+0.236]** |
| `marg_class` | corpus Ramachandran, class-blind, **same 4 states at every residue** | 2.878 | 1.729 | 0.63 | **+0.135 [+0.039,+0.234]** |
| `marg_draw` | 4 random observed (φ,ψ) pairs per residue | 3.969 | 2.223 | 0.40 | +0.629 [+0.502,+0.757] |
| `unif_torus` | **nothing — uniform on [−π,π]²** | 4.778 | 2.698 | 0.12 | +1.104 [+0.949,+1.262] |

### The coordinator's alternative reading is largely correct, and here is its price

* **The library's sequence content is worth 0.135 Å of the 1.594.** Running the descent in a
  *different target's* library costs +0.135 Å; running it in a class-blind corpus codebook
  that gives **every residue the same 4 states** costs the same +0.135 Å. 91.5 % of the
  ceiling is not the library.
* **Of the 1.104 Å gap between the real library and an information-free space of identical
  size, 88 % is generic Ramachandran shape and 12 % is anything sequence-related.**
* **An information-free space still reaches 2.698 Å** — better than the shipped retrieval
  pipeline's 3.213 Å. So "the space contains a good answer" is a much weaker statement than
  it reads: with an oracle CA-RMSD objective, 4 *random* torsion pairs per residue contain a
  better answer than the entire production pipeline emits.

### The descent itself: converged, but the start is privileged

* **Convergence — the claim SURVIVES.** Running the identical descent for 40 sweeps instead
  of 12 gives a **bit-identical** result on 126/126 targets (Δ = +0.0000, 0W/0L). Mean sweeps
  actually used: **2.5**. The descent is not stopping early; 1.594 Å is its converged value.
* **The snap start IS privileged — the claim is WEAKENED.** Four random restarts of the same
  descent (5.3× the oracle queries: 702 vs 132) reach **1.982 Å, +0.388 [+0.272, +0.506],
  31W/90L**. So ~0.39 Å of the reported ceiling comes from initialising at the native's
  nearest states, which is an oracle operation. The honest "what can a strong search find in
  this space without a native-derived start" number is **1.982 Å**, not 1.594 Å.
* For symmetry, four random restarts in the **information-free** space reach 2.454 Å
  (+0.860 vs real). So even stripped of both the library and the start, the search+space
  combination lands at 2.45 Å against a 3.213 Å pipeline.

### Verdict on claim 2

**The literal claim SURVIVES: the k=4 space genuinely contains 1.594 Å structures and the
descent that finds them is converged.** The *interpretation* — "the representation is not the
barrier", implying the sequence-conditioned library is the asset — is **WEAKENED**: the
library's sequence conditioning is worth 0.135 Å, its specific states are k-means noise, and
0.388 Å of the headline is an oracle initialisation. The honest sentence is: *a discrete
torsion space of 4 states per residue — almost any such space — contains ~2 Å structures, and
an oracle-initialised coordinate descent finds 1.6 Å ones in it.*

Sprint 12's assembly oracle died because the search, not the space, was doing the work. This
is not that failure — the space really is expressive — but the *credit* is misassigned.

---

## SEVERITY 1 — the torsion prior (claim 4) and a constant that beats it

`s13/adv_prior.py` → `results/adv_prior.json`.

### E. "SA on the 1-local torsion prior" performs no search at all

**Status: NEW DEFECT (in the framing, not the number).** The prior is exactly 1-local, so its
global minimiser is the per-residue `argmax` of the occupancy — a closed form. Measured:

    exact per-residue argmin RMSD, all 126 targets    3.9692 Å
    "SA on the prior, 5,000-evaluation budget"        3.9692 Å
    identical to 1e-6 on                              126 / 126 targets

The arm labelled `sa_prior` in `coord_FINDINGS.md` C3 is a **deterministic constant chain**,
not a budget-matched search. Its "−0.553 Å at matched budget vs random" is a comparison
between a closed-form baseline and a 5,000-sample search, which is a fine thing to report but
is not what the table says it is.

### F. A ZERO-INFORMATION CONSTANT BEATS THE RANDOM-SAMPLING CONTROL

**Status: NEW RESULT, strengthens C3(b) and weakens the control.** Building every residue at
φ = −63°, ψ = −42° — one ideal α-helix, no library, no search, no budget, no sequence:

| arm | all 126 | helical > 0.5 (n=42) | non-helical < 0.1 (n=70) |
|---|---|---|---|
| ideal α-helix (a constant) | **4.065** | **1.793** | 5.567 |
| SA on the torsion prior | 3.969 | 2.066 | 5.218 |
| random sampling on Legacy, 5,000 evals | 4.522 | 4.229 | 4.827 |

    ideal alpha-helix - random sampling  = -0.457 A [-0.822, -0.090]  66W/60L  drop-top-10 -0.116
    ideal alpha-helix - SA on the prior  = +0.096 A [-0.031, +0.220]

**The sprint's best native-free arm is reproduced to within 0.10 Å (CI including zero) by a
constant α-helix**, which also beats the mandatory random-sampling control by 0.46 Å with a CI
excluding zero. Two consequences: (i) C3(b)'s "it is a helix generator, not a predictor" is
correct and now has its mechanism nailed to a single constant; (ii) **uniform random sampling
in a k=4 torsion space is a weak control** — it is beaten by a one-line constant — so
"beats random sampling" should not be treated as evidence of anything in this representation.
The constant α-helix is the baseline any future arm must clear.

### G. The helix split survives an independent assigner

**Status: SURVIVED.** `I.ss_of` is simplified DSSP on a backbone rebuilt from native torsions.
I built an assigner that never touches torsions or `core.geometry`: helix from the **native CA
trace alone** (|CA_i−CA_{i+3}| and |CA_{i+1}−CA_{i+4}| in [4.7, 5.8] Å and |CA_i−CA_{i+4}| in
[5.6, 6.8] Å over a 5-residue window).

| | `I.ss_of` | CA-distance assigner |
|---|---|---|
| mean helix fraction | 0.315 | 0.263 |
| Pearson between the two | **0.839** | |
| binary agreement at the 0.5 cut | **0.881** | |
| ρ(helix fraction, SA-prior RMSD) | **−0.744** (p=1.7e−23) | **−0.709** (p=1.6e−20) |
| ρ(helix fraction, prior's advantage over random) | −0.646 | −0.611 |
| helical > 0.5: prior − random | −2.163 [−2.631, −1.678] 37W/5L | −2.312 [−2.872,−1.773] 29W/2L |
| non-helical < 0.1: prior − random | +0.391 [+0.101, +0.676] 28W/42L | +0.281 [−0.002, +0.568] 33W/43L |

Sprint 12's "62 % agreement" figure does not reproduce as a threat here: at the 0.5 helix cut
the two assigners agree on 88 % of targets and the correlation is 0.84. **The −0.744 is not an
artefact of the assigner.**

**One caveat I do owe the coordinator:** the non-helical stratum's "the prior is WORSE than
random" is +0.391 [+0.101, +0.676] under `I.ss_of` but **+0.281 [−0.002, +0.568]** under the
independent assigner — the CI touches zero. The safe statement is *"on non-helical targets the
prior's advantage vanishes"*, not *"the prior is worse than random"*.

### H. The obvious confound, priced

Helical targets *are* easier for everyone, but nowhere near enough to explain the effect:

    rho(helix fraction, ORACLE ceiling)      = -0.585     (helical 1.067 vs non-helical 1.918)
    rho(helix fraction, random sampling)     = -0.203     (helical 4.229 vs non-helical 4.827)
    rho(helix fraction, SA on Legacy)        = -0.175
    rho(helix fraction, SA on the prior)     = -0.744     (helical 2.066 vs non-helical 5.218)

The prior's helical/non-helical swing is 3.15 Å against random sampling's 0.60 Å. The
confound accounts for under a fifth of the effect. **Claim 4 SURVIVES.**

### I. The placeholder-angle attack on the prior FAILS

I expected the 1,574 fabricated terminal angles (φ[0] = −60°, ψ[n−1] = −45°, both α-basin) to
be manufacturing the prior's α mode. They are not:

| quantity | all observations | interior residues only |
|---|---|---|
| fraction of the observation pool that is terminal | 0.117 | — |
| α-basin fraction of the pool | 0.5950 | 0.6021 |
| α-state occupancy of the k=4 GENERAL codebook | 0.6784 | 0.6882 |
| GENERAL table shift when placeholders are dropped | 16.0° | |

Dropping every terminal observation *raises* the α occupancy by 0.010. The prior's helical
mode is a real property of the corpus, not a data-entry artefact. **Attack failed; the prior's
α mode is genuine.**

---

## SEVERITY 2 — the exact locality theorem (claim 6)

`s13/adv_locality.py` → `results/adv_locality.json`. The original test used a **0.10 Å**
support threshold on 3 targets at k=8. I re-tested at **machine precision (1e-9 Å)**.

**Status: SURVIVED, and it is stronger than claimed.**

| attack | cells | agreement with "support = {m : i < m < j}" |
|---|---|---|
| all **126** targets, k=4, 6 random base configurations each | ~285,000 (pair,var) | **1.0000** |
| 30 targets, k=8 | ~68,000 | **1.0000** |
| pairs with i=0 or j=n−1 (**chain termini**) | — | **1.0000** |
| variables at **glycine** positions (74 targets) | — | **1.0000** |
| variables at **proline** positions (52 targets) | — | **1.0000** |
| builder convention `omega = 0` (**cis**) | 20 targets | **1.0000** |
| builder convention `omega =` arbitrary fixed value | 20 targets | **1.0000** |

    max |Delta d| for a variable OUTSIDE the support   1.4e-13 A   (floating-point zero)
    min |Delta d| for a variable INSIDE  the support   4.6e-02 A

Eleven orders of magnitude separate "in support" from "out of support". The rule is exact, and
the mechanism is rigid-body invariance rather than a builder accident, which is why changing
ω does not touch it: φ_i, ψ_i, φ_j and ψ_j all generate rotations about axes **through CA_i or
CA_j**, and a rotation about an axis through an endpoint preserves that endpoint's distances.

### The all-atom case: four exact rules, not one approximate one

`qarch_FINDINGS.md` says the rule "does not extend to all-atom distances: those have support
one residue wider". That is right in spirit and imprecise. Measured through
`rep.build_coords`, at 1e-8 Å, 12 targets:

| atom pair | exact support rule | size at separation s | agreement |
|---|---|---|---|
| CA–CA | `i < m < j` | **s − 1** | **1.0000** |
| N–N | `i ≤ m < j` | s | **1.0000** |
| C–C, O–O | `i < m ≤ j` | s | **1.0000** |
| **CB–CB** | `i ≤ m ≤ j` | **s + 1** | **1.0000** |

So it is a *family* of exact theorems, one per atom type, and CA–CA is the narrowest member.
CB–CB — which is where **Legacy's** pairwise terms live — is **two** residues wider than
CA–CA, not one. That sharpens rather than damages the claim, and it should be stated this way,
because "the CA rule is the Legacy rule" (qarch §1) is only true of Legacy's CA terms; its CB
terms obey `i ≤ m ≤ j`.

---

## SEVERITY 2 — the Pauli spectrum and the gradient prediction (claim 5)

`s13/adv_pauli.py`, `s13/adv_pauli2.py` → `results/adv_pauli.json`, `adv_pauli2.json`.

### J. The cross-validation gap I found, and closed

`s13/walsh_xval.py` validates the **transform** (walsh's FWHT over geo's own cached tables:
`abs_diff_mean_weight = 0.0` on every cell) and the **energy table for `legacy` only**, on 4
cells. **No AMBER energy table had ever been rebuilt independently**, yet every AMBER headline
rests on those tables. I rebuilt them through my own path — my own enumeration, my own
MSB-first per-residue base-k decode, `core.amber.single_point` called directly:

| cell | configurations compared | max abs diff | max rel diff |
|---|---|---|---|
| 1A13 L=4 k=4 | 256 (all) | **0.0** | **0.0** |
| 2BFI L=4 k=4 | 256 (all) | **0.0** | **0.0** |
| 1A1P L=4 k=8 | 400 sampled | **0.0** | **0.0** |

Bit-identical, including the indexing. **The AMBER tables and their bit layout are correct.**

### K. Completeness and exactness

* 71 cached tables: **0** with non-finite entries, **0** with a length ≠ `k**L`. `geo_common.clean`
  (which would silently replace non-finite AMBER entries with the finite maximum) **never fired**.
* Parseval residual recomputed independently on 12 (cell, variant) pairs: **≤ 6.5e-16**.
* Inverse-FWHT round-trip reconstruction error: **≤ 1e-12** relative.

### L. The dropped cross-covariances: a genuine limitation, quantified two ways

`Var_pred = Σ_S c_S² Var[∂⟨Z_S⟩/∂θ]` drops `Σ_{S≠T} c_S c_T Cov(∂⟨Z_S⟩, ∂⟨Z_T⟩)`. Two
constructed attacks:

**(1) Sign randomisation.** Flip the signs of the Walsh coefficients at random. `|c_S|²`, the
weight spectrum, and therefore the prediction are *identical bit for bit*; only the cross terms
change. Over 24 randomisations per cell the ratio `meas/pred` sits at **1.000 ± 0.03–0.05**.
So for a *typical* coefficient vector the dropped term really is negligible.

**(2) The achievable range.** For fixed `Σ c_S²` the true variance is bounded by the extreme
generalised eigenvalues of the `∂⟨Z_S⟩` covariance against its diagonal:

| cell | n | achievable ratio range |
|---|---|---|
| 1A13 L=4 k=4 | 8 | [0.298, 2.662] |
| 2BFI L=5 k=4 | 10 | [0.090, 7.764] |
| 1A1P L=4 k=8 | 12 | **[0.0007, 28.10]** |

**The approximation is not a theorem — an adversarial Hamiltonian breaks it by four orders of
magnitude at 12 qubits — but real force fields are not adversarial, and the sprint's write-up
should say "empirically small for these Hamiltonians", which is what it in fact measures.**

### M. WEAKENED — the ratio is θ-noisy for raw AMBER, and the *weight-kernel* form is off

`s13/adv_pauli2.py`, 8 independent θ draws per cell, at two values of `n_theta`:

| cell | n_theta | variant | `pred_exact` ratio | `pred_weight_kernel` ratio |
|---|---|---|---|---|
| 1A13 L=4 k=4 (n=8) | 256 | legacy | 0.997 ± 0.030 | 0.984 ± 0.025 |
| | 256 | **amber (raw)** | 0.989 ± **0.153** [0.734, 1.199] | 0.991 ± 0.153 |
| | 256 | amber_soft | 1.004 ± 0.039 | 1.023 ± 0.040 |
| 2BFI L=5 k=4 (n=10) | 192 | legacy | 0.995 ± 0.045 | **1.494 ± 0.072** |
| | 192 | **amber (raw)** | 1.074 ± **0.210** [0.824, 1.454] | 1.259 ± 0.244 |
| | 192 | amber_soft | 1.005 ± 0.091 | **1.441 ± 0.145** |
| 1A1P L=4 k=8 (n=12) | 128 | legacy | 0.990 ± 0.035 | **1.367 ± 0.045** |
| | 128 | **amber (raw)** | **1.882 ± 4.181** [0.219, **12.227**] | 1.882 ± 4.181 |
| | 128 | amber_soft | 1.028 ± 0.104 | **1.415 ± 0.134** |

Note the direction of the `n_theta` dependence, which is the diagnostic: doubling `n_theta`
halves the θ-sd for every conditioned quantity (2BFI amber_soft 0.091, legacy 0.061 → 0.045)
but **raises** it for raw AMBER on the spiked k=8 cell (0.226 → 4.181). A variance estimator
whose spread grows with the sample is a heavy-tailed estimator without a usable mean — which
is the same delta-spike pathology the geo agent's own §4b correction identifies, now visible
in the gradient statistic rather than the spectrum.

Three quantitative caveats on "matches measured at median ratio 0.993/0.915 with no free
parameter":

1. **The two headline numbers come from two different predictors.** In the current
   `geo_pauli.json` (95 cells): `ratio_meas_over_pred_exact` medians are legacy **1.006**,
   amber **0.913**; `ratio_meas_over_pred_weight` medians are legacy **1.278**, amber
   **0.915**. The quoted 0.915 matches the *weight-kernel* predictor; the quoted 0.993 matches
   the *exact-per-string* predictor. They should be quoted as one pair from one predictor.
2. **The version that works is not the one the claim names.** `pred_exact` uses the full
   per-string variance array `Var[∂⟨Z_S⟩]`, not the weight spectrum. The
   "spectrum × ansatz kernel" form named in the claim — the one that carries the Cerezo-style
   argument — has a **median ratio of 1.278 for Legacy over 95 cells**, and reproducibly
   **1.50 ± 0.09** on 2BFI L=5 k=4. The weight-averaging step, not the cross-term drop, is
   where the error lives.
3. **Raw AMBER's ratio is not a converging estimator**: θ-sd 0.15–4.18, a single cell spanning
   0.219 to 12.227 across 8 θ draws, and *increasing* with `n_theta`; over the 95 cells in
   `geo_pauli.json`, `ratio_meas_over_pred_exact` for `amber` spans **0.004 to 9.076**. A
   median of 0.913 with that spread is a statement about the median, not about the prediction.
   Conditioning fixes it: `amber_soft`'s θ-sd is 4–40× smaller and shrinks with `n_theta` as it
   should. **Every AMBER gradient-variance number in the sprint should be quoted from
   `amber_soft`, not from raw AMBER** — one more consequence of the geo agent's own §4b
   correction, which should be cited alongside it.

**Verdict: the exactness claims (Parseval, reconstruction, table completeness and indexing,
independent reproduction of the AMBER tables) all SURVIVED. The prediction chain SURVIVES in
its `pred_exact` form with a reproducible ratio of 0.99–1.09, and is WEAKENED in the
weight-spectrum form the write-up names, where Legacy's median ratio is 1.278.**

---

## SEVERITY 2 — the annealer (claim 3)

`s13/adv_sa.py` → `results/adv_sa.json`, `adv_sa_report.json`. **42 targets (every third of
the 126), k=4, matched 5,000 Legacy evaluations except where stated.** My hypothesis was that
`coord_search.anneal` measures a broken annealer: it runs 32 *independent* chains for
`(5000−32)//32 = 155` single-residue moves each, so each chain sees 155 proposals over ~13
residues while the random control sees 5,000 draws.

### N. The annealer works. My attack on claim 3's premise FAILS.

| arm | budget | mean Legacy E reached | ΔE vs random | E lower than random on |
|---|---|---|---|---|
| random (5,000 uniform draws, argmin) | 5,000 | −15.796 | — | — |
| `coord_search.anneal`, batch 32 (**the sprint's arm**) | 5,000 | −18.892 | **−3.095** | **41 / 42** |
| the same annealer, batch 8 | 5,000 | −19.527 | −3.730 | **42 / 42** |
| multi-restart coordinate descent on Legacy | 5,000 | −19.597 | −3.801 | **41 / 42** |
| `deep`: annealer at a 50,000-evaluation budget | 50,000 | **−20.397** | **−4.601** | **42 / 42** |
| `random20x`: 50,000 uniform draws | 50,000 | −17.705 | −1.909 | 39 / 42 |

**`coord_search.anneal` is a functioning optimiser.** It beats the random control's minimum
energy on 41 of 42 targets by 3.1 kcal/mol, and a 10× budget buys a further 1.5. The batch
structure costs something (batch 8 beats batch 32 by 0.6 kcal/mol) but not enough to matter.
**Claim 3's premise — that the arm measures the objective and not the optimiser — SURVIVES.**

### O. WEAKENED — but the directional conclusion does not survive the stronger optimisers

The RMSD side, paired against the same random control:

| arm | ΔRMSD vs random | 95 % CI | W/L | **drop-top-10** |
|---|---|---|---|---|
| annealer batch 32 (**reproduces C3's `sa_legacy`**) | **+0.106** | [−0.175, +0.390] | 17/24 | **+0.491** |
| annealer batch 8 | −0.060 | [−0.310, +0.195] | 22/19 | **+0.281** |
| coordinate descent on Legacy | −0.044 | [−0.322, +0.252] | 22/19 | **+0.345** |
| `deep` (10× budget, lowest energy of all) | −0.059 | [−0.297, +0.192] | 22/19 | **+0.263** |
| `random20x` (10× budget, random) | +0.135 | [−0.152, +0.456] | 18/23 | **+0.501** |

    rho(energy gained below random, RMSD gained above random), 42 targets = -0.133, p = 0.40
    within-target rho(E reached, RMSD emitted) across the six optimisers  = -0.006 (median -0.059)

Three things follow:

1. **C3(a) reproduces in sign and size** on my subset: the same annealer gives +0.106 Å where
   the coordinator reports +0.102 Å over 126. The measurement is sound.
2. **It does not survive being made stronger.** Three optimisers that reach 3.7–4.6 kcal/mol
   *lower* energy than the sprint's arm return *better* structures (−0.044 to −0.060 Å) — and
   those gains are themselves concentration artefacts, since **every arm's drop-top-10 is
   positive** (+0.26 to +0.50). Every CI in the table spans zero. At this budget the arms are
   not distinguishable.
3. **There is no depth-of-optimisation → structural-damage relationship at this budget.**
   ρ = −0.133 (p = 0.40) across targets, and −0.006 within targets across six optimisers that
   span 4.6 kcal/mol. So the +0.102 Å **cannot be attributed to "optimising Legacy harder"**:
   in a 4¹³ ≈ 6.7 × 10⁷ configuration space, 5,000 evaluations is 7 × 10⁻⁵ % of the space.

**This is consistent with `qarch_FINDINGS.md`'s own mechanism**, which says the turn is
controlled by the fraction of the space explored and has not occurred on n ≥ 12 targets at
≤ 0.36 %. My arms sit far inside the not-yet-turned regime. **So the claim "a better optimiser
on this objective returns a worse structure" is carried by qarch's exact-enumeration curve
(262,144 configurations; certified global optimum +0.139 Å worse than random on 9 fully
enumerated targets) and by geo's SPSA-on-AMBER result — not by the SA arm.** C3(a) should be
restated as *"Legacy's energy and structural quality are uncorrelated in this representation
(within-target ρ = −0.006 across six optimisers spanning 4.6 kcal/mol), and the exact optimum
is worse than random"* rather than *"optimising Legacy is worse than random sampling"*, which
is a +0.10 Å effect with a CI spanning zero that reverses under a stronger optimiser.

The **larger conclusion is untouched and, if anything, cleaner**: no arm, at any budget or
optimiser strength tested, gets below 4.36 Å, against a 1.594 Å ceiling and a 3.213 Å shipped
pipeline. The recognition failure is the finding; the sign of the 0.10 Å is not.

---

## SEVERITY 3 — a wrong baseline in `walsh_FINDINGS.md` §4.2

**Status: WEAKENED — the effect is real, and it is about half the size quoted.**
`s13/results/adv_encoding.json`, recomputed from the walsh agent's own `walsh_encode.json`.

§4.2 reports that at k=8 an energy-ordered Gray relabelling "buys **−16 to −21 %** against
plain binary on 3 of 3 targets". The correct baseline for "does the *ordering key* matter" is
the agent's own 200-relabelling random null, which is printed in the same table but not used
for the headline. All **6** k=8 cells in the JSON (the write-up quotes 3):

| baseline | mean effect | range | significance |
|---|---|---|---|
| `energy_gray` vs **binary** (the quoted comparison) | **−17.2 %** | −6.0 % … −21.6 % | — |
| `energy_gray` vs the **random-relabelling null** | **−10.9 %** | −7.7 % … −14.1 % | mean z = −1.98, Stouffer z = **−4.84** |
| plain **binary** vs the same null | **+1.62 σ** | — | binary is an unusually *bad* relabelling |
| `energy_gray` vs the **best of 200 random relabellings** | worse | — | **0 / 6 cells** |

Three corrections: (i) the gain against a fair null is 10.9 %, not 16–21 %, because plain
binary sits 1.6 σ *above* the null and half the quoted effect is binary being bad rather than
`energy_gray` being good; (ii) the range across all six available cells reaches −6.0 %
(1DEP), where plain Gray is also *worse* than binary, so "on 3 of 3 targets" is a
three-cell subset of a six-cell result; (iii) `energy_gray` is beaten by the best of 200
uniform random relabellings on **6 of 6** cells.

The effect nonetheless **survives** at Stouffer z = −4.84 over 6 cells, and the design rule
("optimise the encoding only when log₂k ≥ 3") is unaffected. It should be quoted as −11 % over
a random relabelling.

---

## Cross-cutting recommendation

`geo_FINDINGS.md` §11 records "the variational arms beat uniform random sampling at matched
budget — SUPPORTED on Legacy, one cell". Finding **F** above shows that uniform random
sampling in this representation is beaten by 0.46 Å [−0.82, −0.09] by a **constant ideal
α-helix**. Any arm whose positive result is "beats uniform random sampling" — the VQE arms
included — should be re-reported against the constant-α-helix baseline (4.065 Å mean,
1.793 Å on the 42 helical targets) before it is called a positive. That is a one-line control
and it is now the relevant one.

A minor note on `walsh_FINDINGS.md` §2: P3 is presented in a section headed "PRE-REGISTERED
PREDICTIONS (written before the AMBER run landed)" and the text then says P3 "is already
confirmed at k=4 by §1.1 **before it was written down**". A prediction confirmed before it is
recorded is a derivation, not a pre-registration. P1/P2/P4 appear genuinely prospective.

---

## Full ledger

| claim | verdict |
|---|---|
| energies see the same coordinates the RMSD is measured on | **SURVIVED** — bit-identical (max dev 0.0) through `build_ca`, `build_backbone`, `rep.build_coords` and `core.amber` |
| the exact CA–CA locality theorem (support = `j−i−1` strictly between) | **SURVIVED at machine precision** — 1.0000 on 126 targets, k=4 and k=8, termini, GLY, PRO, cis-ω and arbitrary ω; 1.4e-13 Å outside vs 4.6e-2 Å inside |
| all-atom distances are "one residue wider" | **SURVIVED, SHARPENED** — four exact rules; CB–CB is `i ≤ m ≤ j`, i.e. **two** wider |
| the enumerated energy tables are complete and correctly indexed | **SURVIVED** — 71 tables, 0 non-finite, 0 wrong length, `clean()` never fired; AMBER tables rebuilt independently, **max abs diff 0.0** on 912 configurations |
| Parseval / reconstruction exactness | **SURVIVED** — residuals ≤ 6.5e-16, round-trip ≤ 1e-12 |
| the gradient prediction's dropped cross-covariances | **SURVIVED empirically** (sign-randomised null 1.000 ± 0.03–0.05) but **is not a theorem** — achievable ratio range [0.0007, 28.1] at 12 qubits |
| "spectrum × ansatz kernel predicts gradient variance, ratio 0.993/0.915" | **WEAKENED** — the two numbers come from two different predictors; the *weight-kernel* form has median ratio **1.278** for Legacy over 95 cells and 1.50 ± 0.09 on one cell; raw AMBER's ratio has θ sd 0.15–0.32 and spans 0.004–9.08 across cells |
| the k=4 torsion space contains 1.594 Å structures | **SURVIVED** — reproduced exactly; the descent converges in 2.5 sweeps and is bit-identical at 40 sweeps |
| "the representation is not the barrier" (the library is the asset) | **WEAKENED** — a wrong target's library costs +0.135 Å, a class-blind corpus codebook with the same 4 states at every residue costs +0.135 Å, a *uniform-random* space reaches 2.698 Å; 0.388 Å of the headline is the ORACLE snap start |
| the annealer in `coord_search` is optimising | **SURVIVED** — beats random's minimum energy on 41/42 targets by 3.1 kcal/mol |
| "optimising Legacy is worse than random sampling" (+0.102 Å) | **WEAKENED** — reproduces (+0.106) but reverses under three stronger optimisers (−0.04 to −0.06), every CI spans zero, every drop-top-10 is positive, ρ(depth, damage) = −0.133 (p=0.40), within-target ρ(E, RMSD) = −0.006 |
| "the prior's advantage is entirely a helix artefact" | **SURVIVED** — ρ = −0.709 under an independent CA-only assigner (vs −0.744); assigners agree 88 % at the 0.5 cut; the easiness confound explains under a fifth of the effect |
| "on non-helical targets the prior is WORSE than random" | **WEAKENED** — +0.281 [−0.002, +0.568] under the independent assigner; say "the advantage vanishes" |
| the prior's α mode is a placeholder-angle artefact (**my** hypothesis) | **REFUTED** — dropping every terminal observation *raises* α occupancy 0.678 → 0.688 |
| `sa_prior` is a budget-matched search | **REFUTED** — it equals the closed-form per-residue argmax on **126/126** targets |
| uniform random sampling is an adequate control | **REFUTED** — a constant ideal α-helix beats it by 0.457 Å [−0.822, −0.090] |
| "n_res·log2(k) qubits" | **DEFECT** — residue n−1 is inert for the CA trace; `log2(k)` qubits per chain are dead (25.9 → 23.9 mean at k=4) |
| the library is "sequence-conditioned" | **DEFECT** — 2.38 distinct state tables per 12.96 residues; 84 % of residues share one GENERAL table; = the number of residue classes on 126/126 targets |
| `energy_gray` buys −16 to −21 % at k=8 | **WEAKENED** — −10.9 % against the correct random-relabelling null (half of it is binary being a bad code), 0/6 cells beat the best of 200 random relabellings; survives at Stouffer z = −4.84 |

---

## What I did NOT establish

* I did not attack `qarch_FINDINGS.md`'s enumeration arms, `geo`'s ANOVA/Sobol decompositions,
  the CVaR arms, or `lit_FINDINGS.md` at all.
* I did not run the ceiling control at k = 8/16/32; the k=4 result is the one the sprint
  quotes, but the credit split between "library" and "search" may differ at larger k.
* The annealer audit is on 42 of the 126 targets (every third, `adv_sa.json`); it reproduces
  C3(a)'s effect size but its CIs are correspondingly wider than the coordinator's.
* I did not audit any AMBER-objective *search* arm, only AMBER's energy tables and spectra.
* I did not re-derive the locality theorem analytically for side-chain χ variables
  (`chi_bits=True`); every arm in this sprint uses `chi_bits=False`.
* No AMBER-objective search arm was audited (`sa_amber_nc` was run on a subset only).
