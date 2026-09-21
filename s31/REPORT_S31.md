# Sprint 31 — CVaR-VQE Protein Folding: the objective, the readout, and where the information is not

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. Assembled as
results land so nothing is reconstructed from memory at the end.

Branch `s26` · instrument: `tuning126`, 126 targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Charter: `s31/BRIEF.md` (verbatim, 1,660 lines) · Ledger: `s31/LEDGER.md` · State: `s31/STATE.md`
Contract: `s31/S31_CONTRACT.md` (32 rules) · Verifier: `s31/s31_verify.py` · Multiplicity:
`s31/MULTIPLICITY.md`

---

## 0. The answer, up front

[PENDING — written last.]

---

## 1. Endpoint and reproducibility

### 1.1 The endpoint did not move

**Production remains 3.2105 Å** (`s29/results/s29_O_chain_rows.jsonl :: item=prod`), mean built-chain
Cα RMSD over the 126 `tuning126` targets. The CA point cloud is **3.0483 Å** and the *set mean* is
**3.5507 Å** — three different objects, and this report names the basis in the same sentence as every
number.

`benchmark60` is **sealed**; its single pre-registered pass is spent and two guard sites enforce it
(`core/pipeline.py:386-388`, `core/bench.py:1074`, `core/pipeline.py:1664`). Folds and clusters were
not regenerated.

### 1.2 The projection is deterministic and chaotic — there is no seed, and the record said otherwise for three sprints

S28, S30's ledger, S30's report, S30's verifier and this sprint's own charter all describe *"an
unpinned multi-start projection seed."* **There is no RNG anywhere on the projection path.**
`core.project.fit_multi` loops over four fixed starts and takes a strict argmin; `fit_prior` is a
deterministic L-BFGS-B call. Reprojecting the same cloud twice is **bit-identical — max |ΔCA|
exactly 0.0**, unchanged under BLAS thread count.

**The real defect is conditioning, and it is worse than a seed would be:**

```
1A13's four lam=0 objectives   0.5091924865 / 0.5091925170 / 0.5091924488 / 0.5091924188
                               -> spread 1e-7, i.e. NOISE
the same four branches at lam=0.3   0.947 / 0.519 / 0.562 / 2.396
                               -> spread 1e-1
```

A decision taken where the objective **cannot** discriminate fixes an outcome where it **can**.
Measured across the benchmark: the median λ=0 branch margin is **3.1e-7**, and **73 of 126 targets
(58%)** have their branch chosen below a 1e-6 margin. A **1e-14 relative** perturbation of the input
cloud moves the emitted chain by a median of **1.6e-3 Å with a 0.511 Å tail** (2LNG), reproducing the
historical 0.517 Å discrepancy exactly. **Amplification ~1e13. Not one target of 126 is unchanged to
1e-9.**

**Three rules follow and were binding on every lane from the moment they were measured:**

1. Reprojection is reproducible **only from bit-identical clouds**. Both sides of any built-chain
   contrast must be projected **in the same job from the same stored clouds**.
2. **Canonical stays 3.2105 Å** — the run the endpoint was declared from, and the only one whose
   clouds are persisted per target. Not retro-fitted.
3. **No built-chain claim below 0.0107 Å** — the spread across the five circulating values.

**A consequence nobody had drawn:** the measured 0.92 cloud→chain transfer is **an average over a map
that is locally chaotic on a substantial minority of targets.** A lane proposing a cloud-level gain
must expect a non-smooth chain response there.

### 1.3 The reference itself is uncertain — and it does not explain the tail

**92.9% of `tuning126` is NMR-determined** (115 solution NMR, 5 electron crystallography, 4 X-ray, 2
solid-state NMR; RCSB GraphQL over all 126). The manifest reference is *"deposited coordinates, model
1"* — an arbitrary member of an ensemble. The ensembles are genuinely wide: mean pairwise CA-RMSD
between deposited models **1.0823 Å**, with model 1 sitting **0.6965 Å** from its own ensemble medoid
on average (max 4.294); 55/111 resolved targets have spread > 1.0 Å and 15/111 > 2.0 Å.

**It does not explain the tail:**

```
corr(production RMSD, ensemble spread)  +0.1118   95% CI [-0.0762,+0.2921]   SPANS ZERO
worst 18    mean RMSD 6.0291   mean ensemble spread 0.9672
other 93    mean RMSD 2.6842   mean ensemble spread 1.1046
tail-minus-rest  -0.1374   SE 0.2421   MDE 0.6783  ->  0.20x MDE   NULL
```

**The tail's targets have, if anything, narrower deposited ensembles.** The same null holds in all
four other arms. FAIL18 is real failure against a reference no worse determined than any other
target's.

**What survives is a caveat on the absolute number only:** a **uniform ~0.70 Å reference term** — a
perfect predictor aiming at the ensemble medoid still scores ~0.70 Å against model 1. In quadrature
`sqrt(3.21² − 0.70²) = 3.13` against 3.21, i.e. **~0.08 Å today and material only near 1 Å**. Because
it is uniform it **cancels in every arm-to-arm delta**, which is this project's actual currency.
**This is explicitly not a reason to re-score against the medoid** — the reference is part of the
sealed instrument. Reported as an uncertainty line; nothing changed.

---

## 2. Scientific hypothesis

The charter's premise was that **the objective is the barrier** — justified by S30's finding that
the same circuit family reaches **0.2516 Å** under an ORACLE objective and **3.4330 Å** under the
deployed native-free one, against a classical average's 3.2071. *"The circuit can express good
solutions; the objective does not point at them."*

The sprint's opening move was to read the code rather than the record, and that moved the diagnosis
twice — first one stage later, then to a place no one had looked:

1. **The readout, not the Hamiltonian, caps the stage.** (§6)
2. **The objective is not a function of the target at all.** (§5)

Both are stronger than the premise they replace, and the second dissolves the charter's central
question rather than answering it.

---

## 3. Mathematical mechanism

Four exact results, each verified numerically, together accounting for every negative in the sprint.

**(i) The deployed CVaR free energy is a convex program with a closed-form global minimiser.**
For `F(p) = CVaR_α(E;p) − T·H(p)` (`core/quantum.py:993`), the Rockafellar–Uryasev lower-tail form
is affine in `p` inside a max, so `F` is convex in `p` and concave in `s`; Sion's minimax gives

```
min_p F = max_s { s - T*log sum_i exp( (s - E_i)_+ / (alpha*T) ) }     [1-D concave]
p*_i  proportional to  exp( (s* - E_i)_+ / (alpha*T) )                  [hinged Gibbs]
```

**The whole 2^n optimisation is pinned by one scalar.** Strictly stronger than T1 — *the entropy
term was added to break T1's degeneracy and the answer is still one number.* Derived twice
independently: by Sion duality (lane L) and by KKT on the simplex (lane A), agreeing exactly.

**(ii) The readout's objective is an exact identity, and half of it is free.** For **any** weights
with `Σ_x w_x = 1` — non-negativity **not** required, so this covers the selection, convex and
affine readouts *and* the uniform average in one formula:

```
|| sum_x w_x W_x - t ||^2_F  =  <w, a>  -  (1/2) w' B w
    a_x  = ||W_x - t||^2_F      ORACLE       (per-candidate quality)
    B_xy = ||W_x - W_y||^2_F    NATIVE-FREE  ((1/2) w'Bw = the weighted dispersion)
```

Verified to **1.66e-11** over 126 × 80 draws, simplex and affine. **The only unknown is `a`.
The readout question and the ranking question are the same problem.**

**(iii) The same identity from the other end, with a measured constant.** In the deployed common
frame, `||X(w) − nat||²/n = c2 + 2⟨c,Dw⟩/n + ||Dw||²/n` with `c = X̄ − nat`, `D = W − X̄`. **`c2` is
a property of the candidate set and no sum-to-one readout can touch it**; a readout beats the set
mean only through the cross term. Measured:

> **A sum-to-one readout helps iff `|cross| > ||Dw||²`. ORACLE convex runs at 2.45×. Every
> native-free rule runs at ≈ 1** — the score's orientation toward the native is real and is
> cancelled, to within **0.01 Å²**, by the dispersion that concentration costs. **The best
> native-free rule captures 3.99% of the ORACLE cross term. That is the size of the entire
> native-free readout channel.**

**(iv) No ranker can pass the coherence admission test.** With `Σa = 1`, in pair-distance space
`e_p = μ_p + Σ_m a_m η_{m,p}` — so **`Σa = 1` passes the pool's common mode through with
coefficient exactly one, whatever the weights are.** `coh` is therefore a function of the readout's
*concentration*, **not of the ranker**:

```
uniform mean in pair space        coh = 1.0000  sd 0.0000   <- exactly, as derived
coordinate average (shipped)            0.9780
argmin by the SHIPPED DIS               0.8288
argmin by a RANDOM pool member          0.8244   <- the shipped cost and a coin, 0.0044 apart
ORACLE best member                      0.6708   <- the only arm under the 0.6931 bar
```

**The entire 43-channel ranking search of S30 was searching a dimension along which the admission
test does not vary**, and the ORACLE ceiling of *all* in-pool ranking barely clears the bar. **The
test can only be passed by leaving the pool's affine hull** — which is exactly what this project's
one confirmed positive, the AMBER relax at k = 30, does.

---

## 4. CVaR-VQE architecture

The deployed stage (`core/pipeline.py:806-877`) is: take the top `2^k` candidates by score, build a
**diagonal** Hamiltonian from their standardised score ranks, minimise a CVaR free energy over a
9-qubit (deployed: **7**-qubit) RY/CNOT statevector circuit, and read the resulting distribution out
in **three** different ways, only one of which emits the answer:

| readout | where | what `p` does | reachable set | capacity |
|---|---|---|---|---|
| **selection** | `core/pipeline.py:795-803` | picks a cell of an argmin arrangement | the distinct pool members | **6.886 bits mean, 6.555 worst** |
| **convex** | `core/pipeline.py:880-895` **SHIPPED, emits the answer** | enters continuously as weights | the **convex hull** | not bit-capped |
| **affine** | `s27/s28_A_amp.py:105-117`, **harness only** | signed amplitudes | **outside** the hull | not bit-capped |

**Two facts about this that the record did not contain.**

**`quantum: bool = False` in `PROD = Config()`** — the 3.2105 Å endpoint is produced with the stage
**off**. And **all three readout classes were conflated in the record**: the charter's *"the circuit
can express good solutions"* rests on numbers from the **affine harness** readout, while a
deployable quantum stage would use the **convex** one. This report distinguishes them everywhere.

---

## 5. Hamiltonian

### 5.1 The finding: it is a target-independent constant

`E = _zrank(pool["sc"][o])`, and `_zrank` (`core/pipeline.py:788-792`) is `rankdata` followed by
standardisation. **The ranks of any 128 distinct values are 1..128, so the standardised vector is a
constant.** Verified on three independent random score vectors:

```
trial 0  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
trial 1  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
trial 2  first5: [-1.718572 -1.691507 -1.664443 -1.637379 -1.610315]   last3: [1.664443 1.691507 1.718572]
```

Since `p*` is a closed-form function of `(E, α, T)` alone and the circuit is seeded at 0, **both the
closed-form optimum and the circuit's output are one fixed weighting curve per α**, identical across
targets. Measured: **`H(p*) = 4.9136` bits at α = 1 with standard deviation `1.4e-4` across the 126
targets**; 6.6392 bits at α = 0.25 with sd 4.5e-3.

> ### The quantum stage carries zero target-specific information.
> The target enters the answer **only** through the readout's own `P` and `W` — never through the
> objective, the Hamiltonian, the CVaR, or the state. The stage answers *"what fixed weight should
> rank `k` receive?"*, which is a **128-number global hyperparameter, not a per-target
> computation.**

**This is strictly stronger than both of the sprint's earlier capacity theorems** — T1 said the
state specifies one integer, R1 said the selection readout's alphabet is 6.886 bits; this says the
state specifies **nothing** — and it explains every other negative at once: why the objective does
not point at good solutions (*it is not a function of the target*), why α is inert beyond reshaping
a fixed curve, why the circuit's 0.930-bit optimisation gap cancels in the mean, and why the convex
optimum over that family converges to the uniform average.

### 5.2 A non-diagonal Hamiltonian is closed by three independent obstructions

1. **CVaR needs an energy per shot.** Only a diagonal `H` gives every measured bitstring a definite
   eigenvalue; an off-diagonal term cannot be attributed to a single outcome. So
   `H = diag(zrank) − λ·W(block)` **is not a CVaR-VQE and cannot be made into one.**
2. **The forced operator is not even linear in the state.** The identity of §3(ii) forces
   `H[w] = diag(â) − B`, whose VMC local energy is the exact gradient of the readout error — a
   **mean-field operator, quartic in ψ**, with no per-shot eigenvalue.
3. **Dimension counting.** A candidate-index register has Hilbert dimension equal to the candidate
   count, so **any** operator on it is 128×128 and its spectrum is a microsecond `eigh`. **No
   Hamiltonian on a candidate-index register can be classically hard.** The clustering literature
   that does this properly uses **one qubit per data point** — 128, not 7 — and shows no advantage
   at 175 points on hardware.

**And the sign was backwards in the proposal I gave lane A.** The identity's native-free half is
**repulsive**: at fixed candidate quality the readout should *maximise* weighted dispersion. Worse,
the attractive branch **cannot produce a distribution at all** — `B` is conditionally negative
definite, so `w ↦ w'Bw` is concave on the simplex and a concave function on a polytope is minimised
**at a vertex**. Consensus-attraction degenerates to *"pick the single best-scoring candidate"*,
**which is the shipped argmin. A theorem, not a measurement.**

---

## 6. Quantum encoding

**R1, in its final form after two falsifications.** The selection readout is
`v = P·(p/Σp); return argmin_i v_i`, with `P` the pairwise Kabsch CA-RMSD matrix. The state enters
only through a linear map then an argmin, so it matters solely through which cell of the arrangement
`{(P_i − P_j)·p = 0}` it lands in.

- **Falsified once (lane C):** the quantifier *"the entire quantum stage"* is false — the **shipped**
  `average_weighted` emits a continuous convex combination. Measured **1.1144 Å (min 0.073) from the
  nearest pool member**.
- **Falsified twice (lane A):** `argmin(P e_j) = j` fails because the pool contains **duplicate
  structures**, for which `P[i,j] = 0` off-diagonal. Reachable vertices **118.45 / 128 mean, 94
  worst** — *equal to the byte-distinct count on every target*. `filter_pool` already dedups
  (`:766-773`); the readout does not.

> **Final form: the SELECTION readout's alphabet is the number of distinct candidates — 6.886 bits
> mean, 6.555 worst — and it is a scored diagnostic arm (`rmsd_vqe_sel`, `:1173`, registered
> `:1609`), not the answer path.**

**Realised value, ORACLE bits delivered:** `p_θ` medoid 1.805, `p*` medoid 1.876, score argmin
1.693, uniform medoid 1.659. The stage moves the selection on **77.8%** of targets and delivers
**0.112 bits over the plain argmin — 1.6% of the alphabet.**

**Index redesign is closed.** Re-labelling cannot raise capacity above `k` bits, and it does not
help in practice either: **Gray coding is a proven no-op** (identical prefix partition on 126/126,
asserted in code; −0.0014 Å at 0.03× MDE); the best index map raises an ORACLE ceiling by −0.0655 Å
and **lowers** the deployable value; **the best map for a deployable partial readout is the random
permutation**; and no index bit carries more than **0.0658 bits** about candidate quality (deployed
map 0.0343, permutation control 0.005). **The best available map doubles a quantity that is 0.03
bits.** The mechanism is a conflict by construction: differentiating the cells makes the best cell
better and the chosen cell worse, because *a coherent cell is a concentrated set* — §3(iii)
arriving at the index.

---

## 7. Objective geometry

**The circuit does not reach the closed form, and the readout cannot tell.** `KL(p_θ‖p*) = 0.930
bits mean, TV 0.378`; `run_cvar_vqe` is **strictly worse in 12/12** verification cells and in
**36/36** targets on the real instrument. **The gap is expressivity, not optimisation**: 80 → 2000
Adam iterations move it by nothing (0.025764 → 0.025944); only depth and restarts close it, and at
layers = 12 / 400 iters / 8 restarts (~60 s) it is still short — against **microseconds** for the
closed form. *21 parameters cannot cover a 127-dimensional simplex.*

**And yet the two readouts differ by −0.0044 Å at 0.088× MDE** while differing **per-target by
0.102 Å in absolute value**. *The optimisation gap is real and cancels in the mean.*

**The inductive bias, characterised.** The reachable set is an **MPS Born machine of bond dimension
exactly `2^layers = 8`** (measured Schmidt ranks 2, 4, 8, 8, 4, 2). Against the objective's own hinge
shape, `R²` is **0.0071** for a random θ, **0.0538** for the best of 400 draws, and **0.467** for the
optimised `p_θ` — *a ceiling at about half the right shape, not a bias toward anything structural.*

**α is inert on most of the benchmark.** The deployed table is
`VQE_LFO = {0:(1.0,0.3), 1:(0.25,0.3), 2:(0.25,0.3), 3:(1.0,0.3), 4:(1.0,0.3)}` — **T = 0.3 on every
fold**, and **α = 1 on three of the five pinned folds, where `CVaR_1 = ⟨E,p⟩` and CVaR is inactive by
construction.**

**A shipped docstring corrected.** `core/quantum.py:997-1002` asserted the entropy collapse *"is a
property of CVaR, not of the optimiser."* Measured at T = 0 on the real pool: **0.258 bits at α = 1
(78 targets), 3.596 bits at α = 0.25 (48 targets)**. At α = 1 the minimiser is the unique argmin
vertex, so the collapse **is** a property of CVaR; at α < 1 the argmin set is `{p : p_x₀ ≥ α}`, a
face of **positive volume**, so the objective does not determine `p` at all and **the entropy
belongs to the optimiser** — the opposite of the sentence. *"For ANY alpha" is false by theorem.*
Corrected in place with the original quoted.

---

## 8. Free-energy work

**Closed by derivation. No thermodynamics was computed**, which is the charter's §7A hierarchy
working as intended.

**The gating question — what is `S` a function of?** `S(x) = Φ_{β,B}[U_seq](x)`. Its arguments are
the candidate, the sequence (only through atom typing in `U`), and `β` plus the basin map, which are
universal. **There is no fourth argument**, so `S` is an *operator* on (sequence, pool), not a third
source.

**Lemma B1, verified rather than assumed:** the shipped `amber14/protein.ff14SB.xml +
implicit/gbn2.xml` is **reflection-invariant** — ff14SB's **1340 `PeriodicTorsionForce` phases sit at
distance 0.000e+00 from {0, π}**, and there is no CMAP.

**Corollary B1:** `F`, `E`, `S` and every `dF/dT` are rotation-, translation- **and
reflection**-invariant single-structure observables ⇒ by G1 each is a **function of the distance
map**, inside the class S30 closed by theorem and measured empty on 43 channels. The enthalpic half
was already closed by measurement; **the entropic half is closed by the same theorem that closed
contact topology and Rg, because `S` is achiral for exactly the reason `E` is.**

**Corollary B1′, which nobody had connected:** an ANM/GNM Hessian is built from **pairwise
distances**, so its spectrum, its log-determinant — *which is the harmonic configurational entropy
`S_harm = const − ½ ln det′ H`* — and every spectral-graph observable derived from it are
distance-map functions. **This closes the elastic-network, normal-mode and energy-landscape-curvature
family of charter §14 without computing one of them.** Basin populations, ensemble reweighting and
temperature-dependent ranking are monotone transforms of single-structure free energies and close
with it; local configurational entropy closes twice.

**A registered falsifier fired and its matched control inverted the reading.** F1's 1e-6 threshold
fired at **3.673e-06** — but a pure **rotation**, equally analytically exact, gives **5.377e-06,
larger** (ratio 0.68×), and the only force groups that could carry chirality are exact to
**2.0e-16**; the residual lives entirely in two terms invariant under *any* isometry. **A registered
threshold set at an absolute value, on a quantity only a matched control can read, is a mis-set
threshold — not a falsification.** Recorded with the original standing; F1 is never quoted without
the rotation control in the same sentence.

**Left open on price, not theory:** a genuine **multi-structure** barrier `F‡(x_a → x_b)` depends on
the *path*, not on `D(x_a)` and `D(x_b)`, so **G1 does not bind it.** Not closed, not run — 10³–10⁵
trajectories per target against a one-AMBER-process box.

---

## 9. Torsion work

**Closed, and the reason is worth more than the closure.**

**The decomposition theorem:** for an ideal-geometry backbone, the point reflection acts on torsions
**exactly** as `(φ,ψ) → (−φ,−ψ)`, so `T_even` is reflection-invariant and by G1 a distance-map
function. **All of a torsion channel's escape from G1 lives in `T_odd`.** Crossed with separability,
the only cell no theorem closes is **chiral AND non-separable**.

`LEG_torsion` is **30.0% odd** — genuinely chiral, so the registered G1 check does not fire. But the
skill sits in **different halves on the two problems**:

| | EVEN half (G1-closed) | ODD half (G1-escaping) |
|---|---|---|
| coarse triage (meter contrast) | +0.1220, 1.35× MDE | **+0.2698, 2.78× MDE** |
| in-pool selection, 500 band | **+0.2428, 2.51× MDE** | +0.0183, 0.28× MDE |
| in-pool selection, top-75 band | +0.0288, 0.40× | +0.0423, 0.66× |

> **The chiral escape from G1 is real, measurable and not degenerate — and it only works on a
> problem the pipeline does not have.** On the problem it does have, the surviving skill is in the
> **G1-closed** half, and even that dies in band.

- **G2 fired as registered:** in-band `LEG_torsion` is **0.65× MDE**, below the 0.7× bar — where the
  **shipped** cost `DIS` is itself **0.83×, NOT MEASURED in its own band**.
- **G5 fired and is the sharpest number:** a **constant α-helix** at (−57, −47) — plausible, never
  uniform-on-the-torus — **beats every torsion channel on both bands**: ρ +0.3413 against
  `LEG_torsion`'s +0.1676 on the 500 band (paired gap **−0.1722 at 3.17× MDE**), and **−0.0901 at
  1.81× MDE, 5/5 folds** in band. **The Ramachandran channel's in-pool skill is a constant-prior
  effect** — *"this benchmark's peptides are helical"* — and one constant does it better than the
  fold-conditioned channel.
- **M6, the only open cell, is empty:** `XTWIST` (chiral and non-separable by construction) scores
  −0.0104 / +0.0320 and **its own achiral twin beats it** (+0.1124). Scoped honestly: **empty at
  9–16 residues**, which is where S30 left chiral functionals too.

**Mechanism, measured:** `LEG_torsion` of the `RAND_SIGNED` rung sits at **z = +2.363** in the pool's
own torsion distribution, with 70.6% of targets above the pool's 95th percentile, against PROD at
+0.844 and ORACLE `circ_best` at +1.129. **The channel's dynamic range is spent outside the pool.**

---

## 10. Sparse / readout work

[PENDING]

---

## 11. Candidate-index information

[PENDING]

---

## 12. FAIL18 analysis

[PENDING]

---

## 13. Controls and nulls

[PENDING]

---

## 14. The 126-target endpoint

[PENDING]

---

## 15. Statistical analysis

[PENDING]

---

## 16. What was closed

[PENDING]

---

## 17. What improved

[PENDING]

---

## 18. What did not improve

[PENDING]

---

## 19. ORACLE ceilings

[PENDING]

---

## 20. The remaining information bottleneck

[PENDING]

---

## 21. Next-sprint recommendation

[PENDING]

---

## Appendix A — every claim withdrawn this sprint

[PENDING]

## Appendix B — the five "worst 18" strata, as a key

Five different strata are called "the worst 18" somewhere in this project's record. They differ by up
to **0.81 Å on the same-named quantity**. Every such number in this report names its stratum in the
same sentence.

| stratum | defined by | `best1_500` |
|---|---|---|
| `FAIL18` | the filter's own recall (`s12/instrument.py:271-278`) | 2.2842 |
| `defn18` | top-18 by widening gain | **= FAIL18, 18 of 18 — an identity** |
| `worst18_poolmean` | pool mean | 2.2286 |
| `worst18_bestpool` | best pool member | 3.0930 |
| "the genuinely worst 18" | production built-chain RMSD | ORACLE best pool member 2.5298 |
| worst 18 by production RMSD, n = 111 matched | ensemble arm | mean RMSD 6.0291 |
