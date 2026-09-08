# Sprint 13 — PAULI-SPECTRUM findings

**Question.** What is the Pauli-weight spectrum of a real molecular force field, and does it
predict variational trainability?

**The object.** A conformation is a bitstring: residue *r* picks one of *k* discrete torsion
states from the sequence-conditioned, leakage-safe library
`torsion_lib2.library_for(seq, k, exclude_seq=seq)`, encoded in *b* = log₂*k* qubits, so an
*L*-residue peptide is *m = L·b* qubits. Both energy models are **functions of the bitstring**,
therefore **diagonal in the computational basis**, therefore *exactly*

    E(z) = Σ_S c_S · Π_{i∈S} z_i ,      z_i = ±1,   S ⊆ {0..m−1}

— the Walsh–Hadamard expansion, whose coefficients are the Pauli-Z string coefficients.

Code: `s13/walsh_lib.py` (transform, spectra, estimators, encodings), `s13/walsh_verify.py`
(W0), `s13/walsh_legacy.py` (W1), `s13/walsh_amber.py` (W2/W3), `s13/walsh_encode.py` (W4),
`s13/walsh_predict.py` (W5). Results: `s13/results/walsh_*.json`.

**Quantities, all under the uniform product measure on configurations:**

| symbol | definition |
|---|---|
| Var | Σ_{S≠∅} c_S² (= the energy's variance, by Parseval) |
| **V_w** | Σ_{\|S\|=w} c_S² / Var — the **Pauli-weight share** |
| W_mean | Σ_w w·V_w — mean Pauli weight |
| W_eff | 1/Σ_w V_w² — participation number over weights |
| **D_d** | Σ_{S : \|res(S)\|=d} c_S² / Var — the **residue-order share** (functional ANOVA) |
| spread | W_mean / D_mean — qubits per interacting residue |

D_d is a property of the energy on the *k*-ary configuration space; **no index encoding can
change it.** V_w is a property of the energy *composed with the encoding*. This separation is
the spine of everything below and is asserted as an exact identity in `walsh_lib.spectrum`
(`order_weight_bound_holds`: every Walsh subset satisfies d ≤ |S| ≤ d·b).

No native coordinates are read anywhere in this strand. Every number is DEMONSTRATED
(measured) unless explicitly marked otherwise.

---

## ★ HEADLINE CORRECTION — raw AMBER's Pauli spectrum is an artefact of ONE configuration

**Status: PROVEN. Read this before `geo_FINDINGS.md` §4b.** It does not contradict any number
in that section; it reinterprets what those numbers measure, and the reinterpretation reverses
the conclusion.

`s13/walsh_xval.py` → `walsh_xval.json`. Over **141 fully enumerated tables** (both agents'
caches, m = 8…16), the share of the total Walsh variance `Σ_{S≠∅} c_S²` carried by the
**single most extreme configuration** and by the **ten** most extreme:

| model | n cells | top-1 share, median [min, max] | **top-10 share, median [min, max]** | L¹ distance of V_w from the pure-delta binomial |
|---|---|---|---|---|
| Legacy | 74 | 0.026 [0.002, 0.323] | 0.186 [0.018, 0.878] | 1.206 |
| **raw AMBER** | 67 | **0.499** [0.016, **1.0000**] | **0.9961** [0.156, **1.0000**] | **0.232** |

**On a median AMBER cell, 99.6 % of the Pauli-Z variance comes from 10 of 4,096 structures**,
against 18.6 % for Legacy.
On 1A1P L=4 k=8 it is 99.98 % from **one**. A function that is a constant plus a spike at a
single point x₀ has Walsh coefficients c_S = ±(E_max−μ)/2^m for **every** S, so its weight
spectrum is exactly Binomial(m, ½): V_w = C(m,w)/(2^m−1), W_mean = m/2. Measured:

| cell | measured W_mean | delta reference m/2 | L¹(V_w, binomial) | top-1 var share |
|---|---|---|---|---|
| 1A1P L=4 k=8 (m=12) | **6.001** | 6.001 | **0.0003** | 0.9998 |
| 2BFI L=6 k=4 (m=12) | 6.002 | 6.001 | 0.0004 | 0.9997 |
| 1A13 L=6 k=4 (m=12) | 5.947 | 6.001 | 0.027 | 0.994 |
| 1CEK L=6 k=4 (m=12) | 4.017 | 6.001 | 0.935 | 0.062 |
| Legacy 1A13 L=8 k=4 (m=16) | 3.778 | 8.000 | 1.496 | 0.002 |

**Consequence.** "AMBER's mean Pauli weight exceeds Legacy's on 14/14 cells" is measured
correctly and is *not* a statement about the force field's interaction structure. On the
concentrated cells it is a statement that AMBER's variance is a delta spike, and **a delta
spike is maximally global for arithmetic reasons that have nothing to do with chain
connectivity or with van der Waals physics.** The r^−12 wall does not make AMBER non-local by
coupling many residues; it makes AMBER's *variance* a point mass, and a point mass has a flat
Pauli spectrum. The measured spectrum of raw AMBER is therefore **not usable as evidence for
the sprint's central hypothesis**, in either direction.

The distinction is testable and I test it in §3: cap the energy (a labelled MODIFIED
OBSERVABLE, `min(E, p99)`) so the point mass is removed, and ask whether AMBER is *still*
higher-weight than Legacy. If it is, the hypothesis survives on the physics; if it is not, the
14/14 result was the spike all along. **First numbers, exact, m=12, 1A13: raw W_mean 5.726 →
capped-p99 W_mean 4.239, against Legacy 3.797.** So it is partly the spike and partly real —
and geo's `amber_soft` sitting between the two is the same effect seen through a different
compression. Full table in §3.

**Practical consequence, independent of the hypothesis.** `Σ_S c_S²` for raw AMBER is
10²³–10³⁰ kcal²/mol², all of it from structures no optimiser would ever accept. A variational
algorithm measures an observable's Pauli terms to a shot precision set by its norm, so **raw
AMBER as written is not a measurable observable**: at 1e23 variance and 1e-3 relative
precision the shot count is astronomically infeasible. Any AMBER-in-the-loop VQE must first
bound the observable, and that decision is a *modelling* decision that changes the spectrum.

---

## ★ INDEPENDENT CROSS-VALIDATION against the geo/ strand  `s13/walsh_xval.py`

Two agents computed this quantity by separate code paths; there is no literature to check it
against. **Level 1 — my transform on geo's own cached energy tables, 34 cells:**

| comparison | result |
|---|---|
| max abs difference in mean Pauli weight | **0.0 (exact)** |
| max abs difference in any weight share V_w | **0.0 (exact)** |
| max relative difference in Var | 8.1e-16 |
| max Parseval residual (mine) | 4.9e-16 |

**Level 2 — energy tables rebuilt from scratch through my own code path** under geo's library
convention, 4 cells, 1,024–4,096 configurations each: max absolute difference
**1.8e-15 kcal/mol** (relative 1.5e-16). The two implementations are bit-identical.

**One convention difference, which is not a disagreement but must be recorded.** geo builds
the library as `library_for(seq[:L], k, exclude_seq=FULL target sequence)`; my own sweeps use
`exclude_seq=seq[:L]`. Both are leakage-safe; they hold out different pools, so they are
different instruments and their absolute numbers differ (e.g. Legacy 1A13 L=6 k=4:
W_mean 3.511 under geo's convention, 3.797 under mine). Every cross-validation number above is
computed under geo's convention. Nothing in this file compares a number from one convention
with a number from the other.

---

## 0. W0 — verification of the transform  `s13/walsh_verify.py` → `walsh_verify.json`

**Status: PROVEN (verification, not a scientific claim).** Nothing below is worth reading
unless these pass, so they lead.

| check | result |
|---|---|
| FWHT vs the naive O(4^m) Walsh matrix, m=10 random data | max abs **2.2e-16** |
| Parseval, real Legacy tables (4 cells) | rel err **0.0 – 1.7e-16** |
| exact reconstruction (inverse transform − table) | max rel err **1.6e-16 – 4.7e-16** |
| **1-local torsion prior: residue-order share D₁** | **1.000000000000000** on every cell |
| 1-local prior: Walsh mass above weight b = log₂k | **2.4e-31 – 8.6e-32** (i.e. zero) |
| order/weight bound d ≤ \|S\| ≤ d·b, all 2^m subsets, m∈{8,10,12,16}, b∈{1,2,3,4} | holds |

**The prior check is not only a check.** The leakage-safe 1-local torsion prior
`−Σ_i log P(s_i)` is 1-local *by construction in the residue variables*, and the transform
recovers D₁ = 1 to fifteen digits — but its **mean Pauli weight is 1.000 at k=2, 1.69 at k=4
and 2.28 at k=8**, i.e. it is *not* a weight-1 observable. **A "1-local" prior is 1-local in
residues and log₂k-local in qubits.** This is the encoding effect appearing already in the
correctness check, and it is the cleanest one-line statement of the whole encoding-dependence
result (§4).

Independent second anchor: Legacy's own `torsion` term comes back with **D_mean = 1.000000
exactly** on all 6 targets at every m — the transform recovers, without being told, that
exactly one of the eleven Legacy terms is residue-separable.

**Sampling estimator.** `walsh_lib.sampled_spectrum` is an unbiased U-statistic,
V̂_w = (N(N−1))⁻¹ Σ_{a≠b} f_a f_b · K_w(ham(x_a,x_b)) with K_w the Krawtchouk polynomial
[t^w](1+t)^{m−a}(1−t)^a, reduced to a Hamming histogram so all m+1 weights come from one
O(N²) pass. Validated against the exact spectrum at m=12 (1A13, L=6, k=4): exact W_mean
3.797 vs estimated 3.849 at N=4096 and 3.982 at N=16384, L¹ error over the whole weight
distribution 0.023 and 0.086. **The convergence is not monotone in N and the block standard
errors are enormous.** That is not an estimator bug: it is the first appearance of the fact
that dominates §2 — the spectrum of these energies is carried by a small number of very
high-energy clash configurations, so any sampling estimator has heavy-tailed error.
**Exact enumeration is not a convenience here, it is a requirement.** (SUPPORTED.)

---

## 1. W1 — the exact Pauli spectrum of the Legacy potential  `s13/walsh_legacy.py`

**Config.** Six s12 tuning targets, one per leave-fold-out fold, two of them FAIL18 members:
1A13, 1A1P, 2BFI, 1DEP, 1ID6, 1CEK. benchmark60 and dev24 untouched. Targets truncated to
the first *L* residues (the geo/qarch convention); the library for the prefix still holds the
**full** target sequence out. Binary encoding. `core.energy.components_batch`, all 11 terms.
Whole register enumerated: 2^m configurations, 0.155–0.24 ms each.

### 1.1 Legacy is a high-weight observable, and it gets worse with size

Mean over the 6 targets, exact, over three crossing ladders so that qubit count, peptide
length and states-per-residue are not confounded (`walsh_legacy_k4_m18_6t.json`,
`walsh_legacy_k2_k8_m16_6t.json`):

| k | m | L | n | V₁ | V₂ | V₃ | **tail V_{>3}** | **W_mean** | D_mean | spread | prior W_mean |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 |  8 | 4 | 6 | 0.327 | 0.460 | 0.141 | 0.071 | 1.978 | 1.471 | 1.345 | 1.680 |
| 4 | 10 | 5 | 6 | 0.170 | 0.291 | 0.230 | 0.309 | 2.862 | 2.154 | 1.328 | 1.684 |
| 4 | 12 | 6 | 6 | 0.061 | 0.164 | 0.200 | 0.575 | 3.891 | 2.924 | 1.331 | 1.686 |
| 4 | 14 | 7 | 6 | 0.031 | 0.130 | 0.202 | 0.637 | 4.103 | 3.143 | 1.306 | 1.688 |
| 4 | 16 | 8 | 6 | 0.080 | 0.130 | 0.168 | 0.622 | 4.165 | 3.104 | 1.342 | 1.511 |
| 4 | **18** | **9** | 6 | 0.073 | 0.136 | 0.190 | **0.602** | **4.274** | **3.204** | 1.334 | 1.360 |
| 2 |  8 |  8 | 6 | 0.560 | 0.286 | 0.114 | 0.040 | 1.647 | 1.647 | 1.000 | 1.000 |
| 2 | 10 | 10 | 6 | 0.042 | 0.189 | 0.327 | 0.442 | 3.465 | 3.465 | 1.000 | 1.000 |
| 2 | 12 | 12 | 6 | 0.017 | 0.190 | 0.364 | 0.429 | 3.513 | 3.513 | 1.000 | 1.000 |
| 2 | 14 | 14 | 3 | 0.016 | 0.209 | 0.380 | 0.395 | 3.455 | 3.455 | 1.000 | 1.000 |
| 8 | 12 |  4 | 6 | 0.187 | 0.275 | 0.246 | 0.293 | 2.972 | 1.607 | 1.849 | 2.251 |
| 8 | 15 |  5 | 6 | 0.031 | 0.069 | 0.121 | 0.779 | 5.008 | 2.880 | 1.739 | 2.254 |

**At 18 qubits, 7 % of the Legacy potential's variance lies at Pauli weight 1 and 60 % lies
above weight 3.** Classification: **PROVEN** for these cells (exact transform, Parseval to
1e-16, six targets, three ladders).

### 1.1b The spectrum SATURATES — the objective does not become more global with length

This is the most consequential number in the strand and it was not predicted by anyone:

> **W_mean at k=4 goes 1.98 → 2.86 → 3.89 → 4.10 → 4.17 → 4.27 across m = 8…18, and D_mean
> goes 1.47 → 2.15 → 2.92 → 3.14 → 3.10 → 3.20. Both flatten from m = 12 onward. The k=2
> ladder saturates at D_mean ≈ 3.5 across L = 10…14. The tail V_{>3} is flat at 0.60–0.64
> from m=12.**

Legacy is not asymptotically global: **its effective interaction order saturates at ≈ 3
residues, independent of peptide length.** Each added residue adds qubits without adding
interaction order, so `V_{>3}` stops growing and W_mean grows only through the spread factor.
That is exactly what §5.2 predicts from the geometry — a long-range CA–CA distance uses its
whole support but behaves like a low-order function of it — and it is what geo's
sequence-separation decay measures in a different basis.

**Trainability implication, and it is good news that the sprint did not expect.** A cost whose
Pauli weight saturates at a constant is, in Cerezo et al.'s taxonomy, in the *local* class
asymptotically, not the global one. If the barren-plateau mechanism this sprint is testing
were the operative one, Legacy would be *trainable* at scale for that reason. It is not the
operative one on this ansatz (§5.4) — but the observable itself is not the obstruction.
Classification: **DEMONSTRATED** to m=18; **HYPOTHESIS** beyond it.

The decomposition W_mean ≈ spread(k) × D_mean(L) is exact to a few per cent throughout, with
spread(k) fixed by the encoding law of §4.1 (1.000 at k=2, 1.333 at k=4, 1.714 at k=8) and
D_mean set by the physics alone.

### 1.2 Which physical term generates the high-weight structure — Legacy

Mean over 6 targets at m=16 (L=8, k=4). `weighted_var_share` is w_t²·Var(E_t)/Var(E_total)
with the shipped `DEFAULT_WEIGHTS`; `cov_share` is w_t·Cov(E_t, E_total)/Var(E_total) and
sums to 1 across terms.

| term | weight | **cov_share** | W_mean | D_mean | V₁ | tail V_{>3} |
|---|---|---|---|---|---|---|
| **steric** | 4.0 | **0.900** | 4.500 | 3.339 | 0.039 | 0.683 |
| compactness | 0.4 | 0.029 | 2.395 | 1.878 | 0.302 | 0.208 |
| contact | 1.0 | 0.021 | 2.456 | 1.849 | 0.187 | 0.169 |
| electrostatic | 1.0 | 0.020 | 2.450 | 1.901 | 0.291 | 0.179 |
| coop_helix | 2.0 | 0.014 | 3.297 | 2.495 | 0.110 | 0.426 |
| hbond_local | 1.0 | 0.013 | 2.362 | 1.800 | 0.343 | 0.203 |
| solvation | 0.5 | 0.004 | 2.258 | 1.656 | 0.308 | 0.159 |
| torsion | 0.15 | 0.003 | 1.636 | **1.000** | 0.364 | **0.000** |
| hbond_longrange | 3.0 | −0.004 | **5.682** | **4.252** | 0.010 | **0.862** |
| coop_sheet | 2.0 | −0.000 | 5.670 | 4.253 | 0.005 | 0.895 |
| aromatic | 0.8 | 0.000 | — | — | — | — (constant on all 6) |

And at m=18 (L=9, k=4), where the concentration is sharper still:

| term | cov_share | W_mean | D_mean | V₁ | tail V_{>3} |
|---|---|---|---|---|---|
| **steric** | **0.955** | 4.406 | 3.290 | 0.062 | 0.623 |
| contact | 0.017 | 2.481 | 1.874 | 0.270 | 0.202 |
| compactness | 0.016 | 2.247 | 1.799 | 0.367 | 0.186 |
| electrostatic | 0.015 | 2.617 | 1.909 | 0.267 | 0.257 |
| coop_helix | 0.006 | 3.230 | 2.477 | 0.123 | 0.396 |
| hbond_local | 0.003 | 2.535 | 1.946 | 0.314 | 0.250 |
| torsion | 0.002 | 1.347 | **1.000** | 0.653 | **0.000** |
| aromatic | 0.000 | 3.187 | 2.449 | 0.112 | 0.385 |
| coop_sheet | −0.000 | 6.289 | 4.603 | 0.004 | 0.927 |
| hbond_longrange | −0.006 | **5.910** | **4.437** | 0.009 | **0.869** |
| solvation | −0.008 | 2.210 | 1.650 | 0.371 | 0.156 |

**The steric share GROWS with peptide length: 0.814 (m=12) → 0.909 (m=16) → 0.955 (m=18).**

**Three findings, and the first one pre-empts the sprint's central hypothesis.**

1. **The Legacy potential is already a steric potential.** One term — `steric`, the r^−12-like
   clash penalty — carries **96 %** of the total variance at m=18 (90 % at m=16, 81 % at
   m=12), and it is
   the highest-weight term that carries any appreciable variance (W_mean 4.50, 68 % of its
   own mass above weight 3). The sprint's expectation was that AMBER would be steric-dominated
   and Legacy spread. **Legacy is steric-dominated too**, and the fraction *grows* with m.
2. **The two highest-weight terms carry no variance.** `hbond_longrange` and `coop_sheet` have
   D_mean ≈ 4.25 and put 86–90 % of their mass above weight 3 — they are the most non-local
   objects in the model — but they contribute ≈ 0 and ≈ −0.0004 of the total variance on
   these 8-residue prefixes. High Pauli weight and high influence are **different axes**, and
   a term list read for locality alone would rank the model exactly backwards.
3. **`torsion` is exactly residue-separable** (D_mean = 1.000000, zero mass above residue-order
   1) and still has W_mean = 1.64 — the encoding fact of §0 again.

Classification: **PROVEN** for these cells; **SUPPORTED** as a statement about the Legacy
potential in general (six targets, one length, one k).

---

## 2. PRE-REGISTERED PREDICTIONS (written before the AMBER run landed)

The architecture agent's locality result (`s13/qarch_FINDINGS.md` §1, PROVEN) says the support
of the CA–CA distance d_ij is *exactly* the j−i−1 residues strictly between i and j,
contiguous, with zero variables outside. Written in the Walsh basis that becomes three
falsifiable predictions, recorded here **before** the corresponding measurements:

**P1 (support).** For a single pair distance d_ij at separation s, every Walsh coefficient c_S
with S touching a residue outside the s−1 interior residues is **exactly zero**, and D_d = 0
for d > s−1. This is a strictly stronger, exact version of qarch's perturbation test.

**P2 (separation-1 is constant).** Var(d_{i,i+1}) = 0 exactly; a separation-1 pair carries no
torsion information at all.

**P3 (weight from order — the encoding law).** If a residue's dependence is generic within
its own b-qubit block, a subset of residue-order d carries mean Pauli weight

    W | D=d  =  d · b·2^{b−1}/(2^b − 1)          → spread = 4/3 (k=4), 1 (k=2), 12/7 (k=8)

**P3 is already confirmed at k=4 by §1.1 before it was written down**: the measured
W|D = 1.454 / 2.679 / 4.043 at d = 1/2/3 against the predicted 1.333 / 2.667 / 4.000, and the
measured spread is 1.31–1.35 against a predicted 1.333, on every cell and every target. The
d=1 block runs ~9 % high (a residue's own k states are *not* quite generic across its b bits);
d≥2 matches to within 1–2 %. Classification: **STRONGLY SUPPORTED**.

**P4 (the sprint's central hypothesis, stated so it can be refuted).** AMBER's spectral mass
sits at higher Pauli weight than Legacy's. §1.2 already makes this **doubtful**: Legacy is 90 %
steric, and if AMBER is likewise dominated by its own r^−12 wall then the two spectra should be
*similar*, and the hypothesis fails. The coordinator's instruction is to say so plainly and
early if that is what the measurement shows. See §3.

---

---

## 3. W2 — the exact spectrum of genuine AMBER, and which term owns it

`s13/walsh_amber.py` → `walsh_amber.json`. `core.amber.single_point` (ff14SB + GBn2, **no
minimisation**, `threads=1`), whole register enumerated, `components=True`. Measured cost on
this box: **8.66 ms** per distinct configuration without components, **12.51 ms** with —
i.e. the AMBER : Legacy ratio is ~50–70×, not the brief's 20× (the brief's 6 ms figure and
qarch's 12.5–13.2 ms are both reproduced here: 8.66 ms plain, 12.51 ms with the five extra
`getState` calls). Memory polled against a 90.5 % ceiling throughout; no `MemoryError` was
caught and continued past.

### 3.1 One term owns AMBER completely

Mean over the completed L=6, k=4 (m=12) cells. `cov_share` = Cov(E_t, E_total)/Var(E_total),
summing to 1 across terms:

Mean over **all six** targets at L=6, k=4 (m=12), exact:

| AMBER term | **cov_share** | W_mean | D_mean | V₁ | tail V_{>3} |
|---|---|---|---|---|---|
| **nonbonded** (LJ + Coulomb) | **1.000244** | **4.965** | **3.730** | 0.015 | **0.792** |
| bond | 0.000000 | 4.838 | 3.609 | 0.060 | 0.730 |
| angle | −0.000000 | 4.124 | 3.126 | 0.184 | 0.634 |
| torsion | −0.000000 | 1.555 | **1.000** | 0.445 | **0.000** |
| solvation (GBn2) | 0.000000 | 2.864 | 2.230 | 0.142 | 0.345 |

**`nonbonded` carries the whole of AMBER's variance — 1.000244 on every one of the six
targets, identically — and every other term contributes zero to six decimal places.** (The
1.000244 is exactly N/(N−1) = 4096/4095: `cov_share` uses `np.cov`'s ddof=1 against a
population variance. The true share is 1.000000 to the precision measured. The same benign
factor sits on the Legacy `cov_share` column, where it is 1.000015 at m=16.) AMBER's `torsion` term reproduces D_mean = 1.000000 exactly — the third
independent confirmation of the transform, and the exact analogue of Legacy's `torsion` term.
(`bond` and `angle` are *not* constant, because the builder's sidechain placement depends on
the backbone state, but their variance is invisible next to the r^−12 wall.)

### 3.2 Removing the nonbonded term collapses the spectrum below Legacy's

All six targets, L=6 k=4, m=12, exact:

| cell | raw AMBER W_mean | capped p99 | capped p95 | **AMBER − nonbonded** | Legacy total |
|---|---|---|---|---|---|
| 1A13 | 5.726 | 4.239 | 3.527 | **3.021** | 3.797 |
| 1A1P | 5.005 | 4.463 | 3.941 | **2.959** | 3.902 |
| 2BFI | 4.033 | 4.525 | 3.541 | **2.519** | 3.966 |
| 1DEP | 5.010 | 4.452 | 3.897 | **2.792** | 3.605 |
| 1ID6 | 4.016 | 4.522 | 3.892 | **3.827** | 4.102 |
| 1CEK | **6.001** (= m/2, pure delta) | 4.362 | 3.938 | **2.810** | 3.973 |
| **mean** | **4.965** | **4.427** | **3.789** | **2.988** | **3.891** |

**Strip the nonbonded term and genuine all-atom AMBER becomes *more* local than the
knowledge-based Legacy potential** — mean W_mean 2.99 against Legacy's 3.89, and lower on
6 of 6 targets individually (W/L 6–0). All of AMBER's excess Pauli weight, and all of its
variance, is one term.

Note also that **capping is not monotone**: on 1A13 and 1CEK it lowers W_mean (5.73→4.24,
6.00→4.36, removing the spike), while on 2BFI and 1ID6 it *raises* it (4.03→4.53) — those are
the two cells whose top-1 variance share is only 6 %, i.e. the ones that were not
spike-dominated to begin with. The capped mean, 4.43, is the most defensible single number for
"AMBER's Pauli weight net of the dynamic-range artefact", and it still exceeds Legacy's 3.89.

### 3.2b The head-to-head, across the ladder (this is the money table)

Mean over 6 targets per row, exact, **138,752 genuine ff14SB/GBn2 single points** at 32.8 ms
mean. Max Parseval residual over every AMBER cell: **3.5e-16**. AMBER's `torsion` term returns
D_mean = 1.0 to **4.4e-16** on all 26 cells — the transform's correctness anchor, reproduced
26 times on the all-atom model.

| k | m | L | n | **Legacy** W_mean | **raw AMBER** | **AMBER capped p99** | **AMBER capped p95** | **AMBER − nonbonded** |
|---|---|---|---|---|---|---|---|---|
| 4 |  8 | 4 | 6 | 1.978 | 3.655 | **3.076** | 2.537 | **1.950** |
| 4 | 10 | 5 | 6 | 2.862 | 4.840 | **3.896** | 3.395 | **2.416** |
| 4 | 12 | 6 | 6 | 3.891 | 4.965 | **4.427** | 3.789 | **2.988** |
| 4 | 14 | 7 | 6 | 4.103 | 5.481 | **4.672** | 4.057 | **3.184** |
| 2 |  8 | 8 | 3 | 1.647 | 3.318 | **2.708** | 2.625 | **2.070** |
| 2 | 12 | 12 | 3 | 3.513 | 5.249 | **4.165** | 3.554 | **2.790** |
| 8 | 12 | 4 | 3 | 2.972 | 5.997 | **5.201** | 4.547 | **3.579** |

Reading the columns:

* **raw AMBER > Legacy on 7 of 7 rows** — but §★ shows most of that gap is the delta spike.
* **capped AMBER > Legacy on 7 of 7 rows**, by 1.10 / 1.03 / 0.54 / 0.57 / 1.06 / 0.65 / 2.23
  units. **So the sprint's hypothesis survives in its bounded form: net of the dynamic-range
  artefact, genuine AMBER really is the higher-weight observable, W/L 7–0.**
* Along the k=4 ladder the gap **narrows** with m (1.10 → 1.03 → 0.54 → 0.57), because Legacy's
  spectrum is still rising while AMBER's capped spectrum rises more slowly. That is the
  opposite of what an "AMBER is more global, therefore harder at scale" story needs.
  HYPOTHESIS (not measured): the two converge around m ≈ 20.
* **AMBER − nonbonded < Legacy on 6 of 7 rows** (the exception is k=2 m=8, 2.070 vs 1.647).
  Everything AMBER has over Legacy is in `nonbonded`.

### 3.3 The answer to the sprint's central question

Putting §1.2 and §3.1 side by side:

| model | term owning the variance | its cov_share | the term's W_mean | model W_mean after removing it |
|---|---|---|---|---|
| Legacy | `steric` (soft clash penalty) | 0.955 at m=18 (0.909 at m=16, 0.814 at m=12) | 4.41 | — |
| AMBER | `nonbonded` (r^−12 LJ + Coulomb) | **1.000** on all 26 cells | 4.78 | 1.95–3.18, **below Legacy on 5 of 5 ladder rows** |

**The coordinator's prediction is confirmed for AMBER and confirmed *a second time*, for
Legacy, where it was not expected.** Both energy models are, in the Pauli basis, a steric
potential plus rounding error. The physics of chain connectivity *is* what generates
cost-function non-locality in torsion space — but it does so through **one interaction type in
both models**, so the naive form of the sprint's hypothesis ("AMBER is less local than Legacy
because it is more realistic") is **REFUTED**. The models differ in the *boundedness* of that
one term, not in the diversity of their interactions:

* Legacy's `steric` is a **bounded, softened** clash penalty. Its worst configuration is +93
  kcal-equivalent; its top-10 configurations carry a median 18.6 % of the variance.
* AMBER's `nonbonded` is an **unbounded r^−12 wall** evaluated on unrelaxed ideal-geometry
  builds. Its worst configuration is 10¹²–10²⁰ kcal/mol; its top-10 carry a median 99.6 %.

That is the honest statement of the difference, and it is a statement about **dynamic range,
not about locality**. Classification: **PROVEN** on these cells for the per-term ownership;
**STRONGLY SUPPORTED** as the general characterisation (four to six targets, m=12).

---

## 4. W4 — encoding dependence: the same physics, different Pauli structure

`s13/walsh_encode.py` → `walsh_encode.json`. Identical energy tables, relabelled.

### 4.1 The invariance, and the law that governs what can change

**Measured: the residue-order spectrum D_d is identical across every encoding to
1.2e-16 – 2.5e-16** (`D_invariance_max_abs_dev`, every cell). That is the exact statement:

> An index encoding cannot change *which residues interact*. It can only change **how many
> qubits each interaction is spread across**, and only within d ≤ |S| ≤ d·log₂k.

**The law for the spread, derived and then measured (this strand's P3):** if a residue's
dependence is generic across its own b qubits, a residue-order-d subset carries mean weight
W|D=d = d·b·2^{b−1}/(2^b−1). Measured against predicted, exact, over every cell:

| k | b | predicted spread W_mean/D_mean | measured (6 targets × 5 lengths) |
|---|---|---|---|
| 2 | 1 | 1.000 | 1.000 (identically — b=1 leaves nothing to spread) |
| 4 | 2 | **1.333** | **1.306 – 1.345** |
| 8 | 3 | 1.714 | (k=8 cells in `walsh_predict.json`) |

and per order, at k=4: predicted W|D = 1.333 / 2.667 / 4.000 / 5.333 / 6.667 for d = 1…5,
measured **1.55–1.59 / 2.59–2.87 / 3.90–4.13 / 5.33–5.56 / 6.47–6.52**. Agreement is within
1–4 % for d ≥ 3 and runs 16–19 % high at d = 1. Classification: **STRONGLY SUPPORTED**.

### 4.2 How much does the encoding actually buy?

Mean Pauli weight of the *same* Legacy energy under different state→bitstring maps
(200 independent random per-residue relabellings give the null):

| cell | binary | gray | prior_gray | **energy_gray** | random null (mean ± sd) | best…worst uniform perm |
|---|---|---|---|---|---|---|
| 1A1P L=5 k=4 | 2.698 | 2.669 | 2.786 | 2.653 | 2.715 ± 0.057 | 2.669 … 2.786 |
| 2BFI L=5 k=4 | 2.710 | 2.684 | 2.946 | 2.910 | 2.768 ± 0.119 | 2.684 … 2.946 |
| 1DEP L=5 k=4 | 3.461 | 3.253 | 3.220 | 3.363 | 3.353 ± 0.071 | 3.253 … 3.461 |
| 1CEK L=5 k=4 | 2.686 | 2.672 | 2.744 | 2.665 | 2.702 ± 0.084 | 2.672 … 2.744 |
| **1A13 L=4 k=8** | **2.905** | 2.730 | 2.448 | **2.282** | 2.552 ± 0.147 | **2.070 … 3.163** |
| **1A1P L=4 k=8** | **3.135** | 2.907 | 2.831 | **2.649** | 2.869 ± 0.143 | **2.354 … 3.469** |
| **2BFI L=4 k=8** | **3.527** | 3.250 | 3.231 | **2.883** | 3.357 ± 0.215 | **2.754 … 4.100** |

(At k=2, b=1, residue-order and Pauli weight coincide by construction and spread ≡ 1
identically; there is no encoding freedom at all.)

**At k=4 the encoding is worth almost nothing** — the entire binary-vs-Gray-vs-random spread
is under 0.3 units, within ~1.5 sd of the random null, because b = 2 leaves only one bit of
freedom per residue. **At k=8 it is worth a lot**: the best and worst uniform relabelling
differ by **47–53 %** in mean Pauli weight, and simply ordering each residue's states by their
marginal mean energy before laying a Gray code over them (`energy_gray`) buys **−16 to −21 %
against plain binary** on 3 of 3 targets (2.282/2.649/2.883 vs 2.905/3.135/3.527), from a
purely classical, native-free, one-line change. Plain Gray buys −6 to −8 %; ordering by the
prior occupancy is worse than ordering by energy on 3 of 3. **PROVEN** on the cells measured;
**SUPPORTED** as a design rule.

The design rule that follows: *the encoding is worth optimising only when log₂k ≥ 3, and the
right ordering key is the residue's own marginal energy, not the prior occupancy.*

### 4.3 One-hot, and why "the" spectrum of a one-hot Hamiltonian does not exist

Under one-hot (L·k qubits, one per (residue, state)) the physical energy is defined **only on
the feasible subspace**; the 2^{Lk} − k^L infeasible strings must be given values, and that
choice is a modelling decision, not a property of the physics. Two extensions were computed
exactly (`walsh_encode.json`, L=4/5 with k=4 → m=16/20):

| target, L=5 k=4 | **binary** (m=10) | **one-hot, mean-extension** (m=20) | **one-hot, penalty-extension** (m=20) |
|---|---|---|---|
| 1A13 | 2.896 | **10.000** | **1.273** |
| 1A1P | 2.698 | **10.000** | **1.273** |
| 2BFI | 2.710 | **10.000** | **1.273** |

* **mean-extension** (infeasible → mean feasible energy): the feasible set is a measure-zero
  scattering of points, and a sum of point masses has a **binomial** Walsh spectrum, so
  W_mean → Lk/2 = **10.000 exactly**, on all three targets, *regardless of the physics*.
  Exactly the artefact of the headline correction, reached from the opposite direction.
* **penalty-extension** (+λ·Σ_r (popcount_r − 1)²): the penalty is exactly **2-local**, so it
  dominates and W_mean → **1.273** — identical at λ = 1× and 10× the energy range, and
  identical across targets, because the physics is already invisible at λ = 1×.

**The same Legacy energy on the same three targets therefore has mean Pauli weight 2.70, or
10.00, or 1.27 — a factor of 7.9 — depending only on the encoding and the extension chosen.
A one-hot Hamiltonian's Pauli-weight spectrum is a free parameter of the implementer, not a
measurement of the force field.** Any locality claim about a one-hot molecular Hamiltonian
must state λ and the extension. This is a genuinely actionable warning for anyone designing such a Hamiltonian,
and it is the strongest form of the encoding-dependence result.

### 4.4 The project's own worked example

Sprint 12's assembly Hamiltonian was **exactly 2-local** on the same physics, because each
fragment was Kabsch-placed onto a fixed anchor and therefore did not depend on the other
fragments' variables. Torsion space removes the anchor, the chain builds sequentially, and
the same pairwise potential becomes full-register. **The same molecular physics is 2-local or
n-local depending only on how coordinates are produced from the decision variables** — that
is the encoding-dependence statement at its sharpest, and it has a worked example inside this
repository.

---

## 5. W5 — from chain geometry to Pauli structure to gradients

`s13/walsh_predict.py` → `walsh_predict.json`.

### 5.1 The support law, exactly (a strengthening of qarch §1)

qarch proved the support law by perturbation: 64 base configurations, a 0.10 Å threshold,
5,000+ (pair, variable) cells, agreement 1.0000. Here it is tested in the Walsh basis over
**every one of the 2^m subsets** of a fully enumerated register, with no threshold:

| quantity | prediction | measured |
|---|---|---|
| variance of d_ij on Walsh subsets touching a residue **outside** {i+1..j−1} | exactly 0 | **≤ 1.1e-29** (all pairs, all cells) |
| variance of a separation-1 pair distance, relative to its mean² | exactly 0 | **≤ 4.2e-32** |
| highest residue-order carrying mass, pair at separation s | s−1 | **s−1, every pair, no exception** |

**The support law is confirmed exactly.** It is a property of ideal-geometry chain building,
and the Walsh basis makes it an identity rather than a threshold test.

### 5.2 …but the mass is *not* generic on its support — the coordinator's prediction is half right

The prediction was that a pair term's weight is "bounded by, and concentrated near,
(s−1)·log₂k". The bound holds exactly. **The concentration does not**: the mass sits far
*below* the bound.

| separation s | support (residues) | max order with mass | W_mean if generic on support | **measured W_mean** |
|---|---|---|---|---|
| 2 | 1 | 1 | 1.333 | **1.06 – 1.26** |
| 3 | 2 | 2 | 2.667 | **1.79 – 2.43** |
| 4 | 3 | 3 | 4.000 | **1.97 – 2.26** |
| 5 | 4 | 4 | 5.333 | **2.49 – 2.71** |

At separation 5 the pair distance *touches* four residues but behaves like a ~2.5-qubit
observable. L¹ distance of the measured residue-order spectrum from the generic-on-support
prediction grows with s (0.12 at s=3 → 0.89 at s=5): **a long-range CA–CA distance is
overwhelmingly a low-order function of the torsions between its endpoints.** Physically this
is chain rigidity — d_ij is dominated by the additive contribution of each intervening
torsion, with high-order corrections small. This is the answer to the coordinator's "if they
disagree, the gap is interesting": the gap is **not** in how a residue's log₂k qubits share
the dependence (§4.1 shows that part obeys its law to 1–4 %); it is at the residue-order
level, and it means **the |i−j|-locality worry is a bound that the physics does not saturate.**

### 5.3 The additive-pair picture is quantitatively exact

Closing the loop end-to-end on a **SURROGATE** pairwise CA potential (12-6 Lennard-Jones on
CA–CA distances — *not* a force field, and named a surrogate everywhere): predict its whole
spectrum as the variance-weighted mixture of the individual pair spectra, dropping all
cross-covariances between pair terms.

| cell | measured W_mean | additive prediction | L¹ over the whole V_w | Σ(pair variances)/Var(E) |
|---|---|---|---|---|
| 1A13 L=6 k=4 | 4.025 | **4.025** | **0.0000** | 1.000 |
| 2BFI L=6 k=4 | 4.031 | 3.982 | 0.161 | 1.028 |

**Cross-covariance between pair terms carries 0–3 % of the variance.** The chain
geometry → exact support → per-pair spectrum → total spectrum sequence is closed, with the
additive approximation costing essentially nothing. Classification: **PROVEN** on the
surrogate; **SUPPORTED** as an approximation for the real models (both of which are dominated
by a single pairwise term, §3.3).

### 5.4 The trainability prediction — and the honest negative

**Regime audit first** (`s13/results/lit_methods.json`, Cerezo 2021 row, verbatim): the
theorem assumes *"V(θ) is an alternating layered ansatz composed of blocks forming local
2-designs"*, and its local-cost guarantee holds *"only up to O(log n) depth"*. The ansatz in
use (`StatevectorCircuit`, real-amplitude RY rotations, ring entanglement, layers 1–5 at
n = 8–16) is **not** composed of local 2-design blocks — real-amplitude blocks cannot 2-design
even locally — and the ring makes the light cone span the register at depth 1. **Cerezo 2021
therefore licenses no prediction here, and saying otherwise is the error the literature agent
warned would end the conversation.**

**What the spectra do predict, using the geo/ strand's validated factorisation.** geo measured
Var[∂C/∂θ_i] = Σ_S c_S² · v(|S|) with the ansatz kernel v(w) computed exactly, and validated
it at median ratio 0.993 (Legacy) / 0.915 (AMBER) over 56 combinations. Extracting v from
`geo_pauli.json`:

| quantity | measured |
|---|---|
| v(w) dependence on **weight** w | **flat**: decay base 0.786–1.013, v(w_max)/v(1) median 1.00, range 0.12–1.25 |
| v̄ dependence on **n**, layers = 1 / 2 / 3 / 5 | v̄ ∝ 2^(−0.47n) / 2^(−0.63n) / 2^(−0.65n) / 2^(−0.86n) |

**So the prediction is:**

> **Var[∂C] ≈ v̄(n) · Var(E), with v̄ ∝ 2^{−c n}, c ≈ 0.47–0.86 rising with depth, and with
> *no* dependence on the Pauli-weight distribution. The energy model enters only through
> Σ_S c_S² = Var(E). This is exponential concentration in n for *every* observable, including
> a weight-1 one — it is the expressibility mechanism (Holmes 2022), not the cost-globality
> mechanism (Cerezo 2021).**

**The direct answer to my assigned question is therefore NO, and it must be reported as such:
the Pauli-weight spectrum of a real molecular force field does *not* predict variational
trainability on this ansatz at these sizes, because this ansatz's gradient kernel is
weight-flat. n predicts it.** Classification: **DEMONSTRATED** for depth 1–5 at n = 8–16 on
this ansatz; **UNKNOWN** for an ansatz that does satisfy Cerezo's assumptions.

### 5.5 The one thing the spectrum *does* predict, and it is sharp

Since only Σ_S c_S² enters, and §★ measured that **99.6 % (median) of raw AMBER's
Σ_S c_S² comes from 10 configurations out of 4,096**, the prediction is mechanical:

> **Raw AMBER's gradient signal at random θ is, to within a fraction of a percent, the
> gradient of the sampling probability of a handful of specific clash bitstrings. It carries
> no information about the other ~4,086 structures.**

Two falsifiable consequences, offered to the geo/ and optimiser strands rather than measured
here (to avoid duplicating their arms):

1. **CVaR at small α should *destroy* AMBER's gradient, not sharpen it.** CVaR_α keeps the
   lowest-α tail and discards the high-energy configurations — i.e. exactly the 10 structures
   that carry 99.6 % of the Pauli norm. Predicted: ‖∇CVaR_α‖/‖∇⟨E⟩‖ falls by orders of
   magnitude as α → 0 for AMBER and stays O(1) for Legacy. This is the opposite of the
   folklore that CVaR helps here, and it is consistent with Qiu et al. 2026 (`lit` C2).
2. **Capping or softening the potential should change AMBER's trainability far more than
   changing its locality does**, because capping removes ~99 % of Σ_S c_S² while moving
   W_mean by only ~1.5 units (§3.2).

---

## 6. SUMMARY — the eight results, tiered

| # | result | tier | where |
|---|---|---|---|
| 1 | **Raw AMBER's Pauli spectrum is a delta-function artefact**: median 99.6 % of Σ_S c_S² from 10 of 4,096 configurations (67 cells); V_w within L¹ 0.0003 of Binomial(m,½) on the worst cells. Legacy: 18.6 % from 10 (74 cells). | PROVEN | §★ |
| 2 | **Two independent implementations agree exactly** — max abs difference in W_mean = 0.0, in any V_w = 0.0, over 34 cells; energy tables agree to 1.8e-15 over 4 cells. | PROVEN | §★ |
| 3 | **One physical term owns each model.** AMBER: `nonbonded` cov_share **1.000 on all 26 component cells**. Legacy: `steric` 0.955 at m=18, rising with length (0.814→0.909→0.955). Remove AMBER's `nonbonded` and it becomes *more* local than Legacy on **6 of 7** ladder rows. | PROVEN | §1.2, §3 |
| 4 | **The naive hypothesis "AMBER is less local because it is more realistic" is REFUTED**; the bounded form survives. Both models are a steric potential plus rounding error; they differ in the *boundedness* of that one term — dynamic range, not interaction diversity. Net of the artefact (capped p99) AMBER is still the higher-weight observable, **W/L 7–0**, but the gap narrows with m (1.10→1.03→0.54→0.57 along the k=4 ladder). | REFUTED (naive form) / SUPPORTED (bounded form) | §3.2b–3.3 |
| 5 | **The spectrum saturates.** W_mean 3.89→4.10→4.17→4.27 and D_mean 2.92→3.14→3.10→3.20 across m = 12…18. Legacy's effective interaction order is ≈ 3 residues *independent of peptide length*. | DEMONSTRATED to m=18 | §1.1b |
| 6 | **The support law is exact in the Walsh basis** (outside-support variance ≤ 1.1e-29; separation-1 pairs constant to 4.2e-32) **but the mass does not saturate the bound** — a separation-5 pair distance touches 4 residues and behaves like a 2.5-qubit observable. | PROVEN | §5.1–5.2 |
| 7 | **The encoding law**: D_d is encoding-invariant (1e-16), and spread = b·2^{b−1}/(2^b−1) — predicted 1.000/1.333/1.714 at k=2/4/8, measured 1.000/1.31–1.40/1.81–1.94. At k=4 the encoding is worth nothing; at k=8, energy-ordered Gray buys −16 to −21 % on 3/3 targets, and the same energy has W_mean 2.70 / 10.00 / 1.27 under binary / one-hot-mean / one-hot-penalty. | PROVEN (invariance), STRONGLY SUPPORTED (law) | §4 |
| 8 | **The spectrum does NOT predict trainability on this ansatz.** The measured gradient kernel v(w) is flat in weight (v(w_max)/v(1) median 1.00) and decays as 2^{−0.47n}…2^{−0.86n} in n. Cerezo 2021's regime (local 2-design blocks, depth O(log n)) does not hold for a real-amplitude ring HEA, so the theorem licenses no prediction here. What the spectrum *does* predict: Var[∂C] ∝ Var(E), so AMBER's gradient is the gradient of the probability of ~10 clash bitstrings. | DEMONSTRATED (negative) | §5.4–5.5 |

**The one-sentence version.** *The Pauli-weight spectrum of a real molecular force field is
owned entirely by its steric term, saturates at an effective interaction order of ≈3 residues
independent of chain length, changes by a factor of 7.9 under different quantum encodings of
the same physics, and — on the hardware-efficient ansatz actually in use — does not predict
gradient variance at all, because that ansatz's gradient kernel is flat in Pauli weight; what
does predict it is the observable's total variance, 99.6 % of which, for unrelaxed all-atom
AMBER, comes from ten clashing structures.*

---

## 7. WHAT THIS STRAND DID **NOT** ESTABLISH

1. **No claim about barren plateaus.** A BP is an asymptotic statement about Var[∂C] as
   n → ∞ under random initialisation. Everything here is n = 8–18, exact, noiseless. The
   fitted v̄ ∝ 2^{−cn} is a three-point fit per depth taken from another agent's data.
2. **The spectra do not transfer to other ansätze.** §5.4's negative is a statement about the
   ansatz whose kernel was measured. An alternating-layered local-2-design ansatz at depth
   O(log n) without a ring might well show the weight dependence Cerezo predicts; that was not
   built and not tested.
3. **Nine-residue prefixes are not proteins.** Every exact spectrum is on L ≤ 9 (k=4) or
   L ≤ 16 (k=2) prefixes of tuning targets. The trend in m is measured, not extrapolated, and
   nothing here licenses extrapolation to a 40-mer.
4. **The sampling estimator is validated but not trusted at scale.** It is unbiased and it
   reproduces the exact spectrum, but its variance is dominated by the same rare clash
   configurations that dominate the spectrum, so above m ≈ 18 the honest position is that
   these spectra are **not** cheaply estimable for AMBER.
5. **"AMBER − nonbonded" is not "AMBER without sterics".** OpenMM's `NonbondedForce` carries
   the Lennard-Jones repulsion *and* the Coulomb term together; there is no way to remove only
   the r^−12 wall through this API. The arm is named for what it actually removes.
6. **The capped and one-hot arms are MODIFIED OBSERVABLES**, labelled as such everywhere. No
   number in this file calls a capped potential "AMBER".
7. **No optimisation was run, no VQE, no CVaR, no gradients were measured here.** Every
   gradient statement in §5.4–5.5 is a prediction derived from another agent's measured
   kernel, offered for them to test.
8. **The k=8 encoding result rests on three cells** (1A13, 1A1P, 2BFI at L=4); the k=4 null
   result is on 6 targets × 3 cells.
9b. **Two exact AMBER cells at m=16** (1A13, 1A1P at L=8, k=4; 65,536 single points each)
   were still running when this file was written; they append to `walsh_amber.json` as further
   `exact: true` rows and nothing in §3 depends on them. The **sampled** AMBER arm at m=16/18
   was deliberately **not** run: §0 shows the estimator's variance is dominated by the same
   clash configurations that dominate the spectrum, so it would have produced a number without
   a usable error bar.
9. **THE SIX TARGETS ARE NOT SIX INDEPENDENT LIBRARIES, and this affects the geo/ strand
   too.** At L=6, k=4, `library_for` returns a **bit-identical** torsion table for 1A1P, 2BFI,
   1ID6 and 1CEK — none of their first six residues is P or G, so all four fall in the general
   residue class and get the same back-off pool. 1A13 (…G…) and 1DEP (…P…) differ. Verified by
   direct array comparison. Consequences:
   * For **geometry-only** quantities — the pair-support law (§5.1), the pair W_mean table
     (§5.2), the LJ surrogate (§5.3) — the effective n is **3, not 6**, and the four
     duplicated rows in `walsh_predict.json` are identical for that reason and not by
     coincidence.
   * For **energy** quantities (§1, §3, §4) the six targets remain distinct, because the
     sequence enters the energy function directly (Legacy's contact/aromatic/hbond typing;
     AMBER's atom types and charges) even on an identical torsion table. The W/L 6–0 in §3.2
     is therefore a real 6-target result, but it is 6 sequences over 3 libraries.
   * The same is true of any per-target count in `geo_FINDINGS.md` at these lengths.

