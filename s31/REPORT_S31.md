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
3. **No built-chain MEAN claim below 0.0107 Å** — the spread across the five circulating
   values. **And the per-target floor is far higher.** Lane P measured it directly (S31-L18)
   by running the *same operator in two implementations* — `average_weighted(uniform,
   top-75)` against `coordinate_average(top-75)`, whose **clouds agree to 3.6e-14**:

   ```
   built-chain MEAN     3.2108 vs 3.2105, paired +0.0003, 0.03x MDE  -- a null
   PER TARGET  |d|      mean 0.0134   median 0.0026   p90 0.0329   MAX 0.2285 A
                        exactly zero on 0 of 126; above 0.0107 on 28 of 126 (22%)
   ```

   > **A 1e-14 input difference — the same method written twice — moves the built chain by
   > up to 0.23 Å on individual targets while cancelling to nothing in the mean.** Any
   > per-target chain statement below **~0.03 Å** is indistinguishable from re-running the
   > same computation in a different implementation. This is also a **matched null** for
   > the chaos result above, arriving from shipped code rather than an injected
   > perturbation, and with a distribution rather than a single worst case.

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
targets — the ties in a real pool leave a residual of **max deviation 4.07e-02** across the 126, but the vector is otherwise fixed. Measured: **`H(p*) = 4.9137` bits at α = 1 with standard deviation `3.1e-04` across the 126
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

### 10.1 The readout class question, answered — with its information cost attached

Built chain, n = 126, **ORACLE / NOT DEPLOYABLE**, from `s29/results/s29_O_chain_rows*.jsonl`
(nothing reprojected):

```
production (uniform top-75 average)                3.2105
bestm128   ORACLE prefix-average over the 128      2.9027    T1's entire reach
best1_128  ORACLE argmin over 128                  2.1435    7 BITS
hull_128   ORACLE CONVEX weights over the SAME 128 1.8538    128 REAL NUMBERS
best1_pool ORACLE argmin over 500                  1.7078    9 BITS
```

**Combining reaches further than naming at a fixed candidate set. Naming wins per bit.** Both
sentences are needed and only the pair is honest — and the second is the one that governs, because
the project's scarce resource is information, not expressiveness.

> **Every lever comparison in this report carries its information cost in the same sentence as its
> Ångströms.** Reading `hull_128 = 1.8538` against `best1_128 = 2.1435` as *"the readout class is
> the larger lever"* compares **7 bits against 128 free reals**, and that is the error S30 §9.4
> already named as this project's characteristic one. It is a **larger ceiling at unpriced cost**,
> not a larger lever.

**And S30-L11's conclusion survives its own mismatched support.** S30 claimed *"a plain argmin
dominates at every bit budget"* and supported it with 2-of-**75** against 1-of-**128**. The support
was mismatched — but **S30 computed the correct reference curve in the same file**
(`s30_Q_sparse.json :: argmin_ref`: 1.9383 at 8 bits, **1.7108 at 9**), and the sparse arm prices
out at ~**9.1 bits**, where naming gives 1.7108 — **0.36 Å better at matched information. Nothing
in the ladder beats the argmin at a matched budget.** *An earlier draft of this report stated that
the set-matched ladder "inverts" S30-L11. That was a withdrawal that was itself wrong; it is
retracted here.*

### 10.2 The matched-search-size correction, which is a finding

`2-of-128` is a per-target minimum over **C(128,2) = 8128** ORACLE supports; `argmin over 128` is a
minimum over **128**. Comparing them directly is a **64× oracle-search mismatch**. Measured at
matched `K` (**ORACLE / NOT DEPLOYABLE**, CA cloud):

```
ORACLE argmin over 128 singletons   (K=128)     2.1458        --
2-of-128 uniform, 128 RANDOM pairs  (K=128)     2.2229   +0.0771  WORSE  1.11x MDE, 5/5, 44W/82L
2-of-128 uniform, EXHAUSTIVE        (K=8128)    1.9070   -0.2388  better (3.81x MDE)
```

**And it decomposes exactly:**

```
family LEVEL   singles 3.5847  pairs 3.3469   -0.2378   <- error cancellation IS real
dispersion loss (pairing compresses the family)        +0.3149
net at matched K                                       +0.0771  pairs WORSE
```

> **Error cancellation is real and worth −0.238 Å in level — and is then more than consumed by
> averaging compressing the family's dispersion.** This is §11's prefix law, **second independent
> instance, opposite operator**: one lane decorrelated a family and its minimum **grew**; another
> contracted one and its minimum **shrank**.

### 10.3 The native-free channel, measured exactly

§3(iii)'s identity gives the closure. A sum-to-one readout beats the set mean **only** through the
cross term, and:

```
                        cross    ||Dw||^2   emitted rmsd (CA cloud)
unif_prefix_m128       +0.0000    0.0000        3.0435
soft_T4.00             -0.0353    0.0257        3.0391
soft_T1.00             -0.2025    0.3097        3.0442
typicality T1.00       +0.1844    0.2708        3.0957   <- WRONG SIGN
ORACLE convex128      -10.3193    4.2054        1.7977   <- 2.45x
```

> **ORACLE convex runs at 2.45×; every native-free rule runs at ≈ 1.** The score's orientation
> toward the native is **real** and is cancelled, to within **0.01 Å²**, by the dispersion that
> concentration costs. Sharpening buys cross monotonically and pays dispersion faster every time.
> **The best native-free rule captures 3.99% of the ORACLE convex cross term — that number is the
> size of the entire native-free readout channel.**

### 10.4 The constrained-affine middle has no selection rule

- **No knee by ceiling.** The ridge path from uniform to unconstrained affine is smooth and
  monotone (`||w||₁` 1.01 → 28.3 gives ORACLE 2.7318 → 0.0000), and at `||w||₁ = 2.09` the ORACLE
  already beats the argmin over the whole 500. **There is no ORACLE answer to "what is the right
  constraint".**
- **No answer by fit either.** ORACLE weights are not identifiable: out-of-fold **R² = 0.021**
  against a 5% bar, and **applied it is +0.0309 Å worse**.

**So a constrained-affine readout cannot be selected by ceiling and cannot be selected by fit.**

**And sparsity is an output, not a design choice:** the exact objective's optimum is **already
sparse — 6.54 of 128 members with no sparsity penalty imposed.** Sparse-readout work that does not
improve `â` is spending effort on the half of the objective that is already exact.

### 10.5 The affine readout's ceiling is vacuous by rank

`rank(aff{W_x}) = 32.9 ≥ 3n−3`, so the affine hull of 128 candidates **spans the residual space**
and the ORACLE affine ceiling is **exactly 0.0000 Å on all 126 targets** — 127 weights against ~33
dimensions, an over-parameterised interpolator (median `ess` **1.13 of 128**, `neg_mass` 3.38).

> **Therefore "0.2516 Å under an ORACLE objective" is a statement about a REGULARISER, not a class
> ceiling** — and the simplex constraint is precisely the regulariser the uniform average enjoys for
> free. That *derives* "expressivity without an aligned objective is harmful" instead of observing
> it: the uniform average is the heavily regularised special case of the same affine readout, and
> under the deployed objective it **beats** the unconstrained optimisation.

---

## 11. Candidate-index information

**Closed. Re-labelling cannot buy information, and in practice it does not buy anything else
either.**

- **Gray coding is a proven no-op** — the j-bit prefix partition is **identical to binary's on
  126/126 targets** (asserted in code); measured ceiling difference **−0.0014 Å at 0.03× MDE**.
- **Structure-aware indexing helps the ORACLE and hurts the deployable**, by construction: it
  raises the ORACLE value of a partial measurement (2.4216 → 2.2712) and **lowers** its deployable
  value (3.2369 → 3.3721). *Differentiating the cells makes the best cell better and the chosen
  cell worse, because a coherent cell is a concentrated set* — §3(iii) arriving at the index.
- **The best index map for a deployable partial readout is the random permutation** — the same
  phenomenon §12 finds at the prefix length.
- **No index bit carries meaningful information about candidate quality.** Mutual information
  between an index bit and "in the ORACLE-best decile": deployed map's best bit **0.0343 bits**,
  best map **0.0658**, permutation control 0.005. **The best available map doubles a quantity that
  is 0.03 bits.**

**The 128 → 512 register widening** survives its circularity gate at **−1.1811 Å cloud / −1.1762
chain** on a filter-independent tail (**ORACLE / NOT DEPLOYABLE**, p = 0.0001), **62% of the
published −1.9004** which is retired as a near-tautology. But **it only pays if a selector exists**,
and three facts held by three different lanes compose into a prediction that it would not:

1. **On the tail, the score is worse than uninformative:** FAIL18's ORACLE-best member sits below
   rank 128 at frequency **1.000** against an uninformative-score null of **0.744**.
2. **Widening degrades the very selector it requires:** `Var ~ 16/D` is the predicted global-cost
   barren plateau; 128 → 512 costs **4×** in gradient variance.
3. **`operator-consumes-set-mean`:** admitting 384 more candidates to a set the selector cannot
   order **raises the set mean**, consumed at coefficient 1.16.

> **PREDICTION, labelled as one and not run: a DEPLOYABLE 128 → 512 arm should come out WORSE than
> production, not merely flat.** Caveat carried because it cuts the other way: at n = 7–9 we are
> **not yet gradient-limited**, so fact 2 bounds the direction without quantifying the Ångströms.

---

## 12. FAIL18 analysis

**`FAIL18` is not a stratum. It is an identity, and this sprint proved it.**

A pre-registered definition-matched stratum `defn18` — the 18 targets with the largest widening gain
from the top-75 — turned out to be **`defn18 ∩ FAIL18 = 18 of 18. The same targets, effect identical
to four decimals.** Choosing the 18 targets whose best candidate hides below rank 75, then reporting
that their best candidate also hides below rank 128, **is one statement and not two.**

**The rank diagnostic confirms it independently and says something sharper.** Under a score with no
within-pool skill the ORACLE best member's rank is uniform on 1..500, so `P(rank > 128) = 0.744`:

```
other 108           0.417     the score genuinely pulls the answer into the window
worst18_poolmean    0.667     near the uninformative null
worst18_bestpool    0.722     near the uninformative null
FAIL18              1.000     ABOVE the null -- p = 0.744^18 = 0.005 by chance
```

> **`FAIL18` is the only stratum *worse than an uninformative score*, which cannot happen by
> sampling.** That is the tail's mechanism: **not "hard targets hide their answer deeper" but "the
> score has no within-pool skill on hard targets, so the best member lands roughly uniformly and a
> 128-window misses it ~74% of the time."**

**And the excess is mostly definitional.** Regressing per-target widening effect on `best1(128)`
gives slope −0.4356, R² 0.5564; the residuals are **−0.577 for `FAIL18` against −0.196 for the clean
filter-independent tail**. **Two thirds of `FAIL18`'s excess is the stratum definition.**

**Consequences enforced throughout this report:** filter-independent tails are primary (`worst18` by
pool mean and by best-in-pool), `FAIL18` appears only as a labelled diagnostic, and **Appendix B is
the key** — five different strata are called "the worst 18" somewhere in this project's record and
they differ by up to **0.81 Å on the same-named quantity**.

**The reference is not the explanation** (§1.3): the tail's targets have, if anything, **narrower**
deposited NMR ensembles, and the tail-minus-rest difference in ensemble spread is **0.20× MDE — a
null**. **FAIL18 is real failure.**

**Where the tail's prize is, ORACLE / NOT DEPLOYABLE:** the along-`mu` correction (§20.1) is **4×
concentrated on the tail** — FAIL18 6.0195 → 3.7392, **−2.2803 Å**, against −0.5651 on the other 108.
The tail is where the endpoint's mean is decided and where the only priced route would pay most.

---

## 13. Controls and nulls

Charter item 25. **Several results in this sprint were killed by their own lane's control**, which
is the point of the table.

| control | what it was matched to | what it ruled out |
|---|---|---|
| **Shuffled-`B`** (random relabelling of the pairwise matrix) | quality-blind dispersion maximisation, which beat the shipped argmin by −0.2621 Å at 2.74× MDE | **Killed the reason for its own lane's only positive.** Costs 0.0250 Å at **0.54× MDE — NOT A RESULT** — so the mechanism is *"spread the weights over many candidates"*, **not** *"spread them along the real geometry"*. No part of the gain may be attributed to the pairwise structure |
| **Matched random subsets at fixed `K`** | the per-target ORACLE prefix length | **Destroyed `bestm128` as a lead.** Random subsets reach **149%** of the prefix family's gain: the prefix axis is *worse than an arbitrary 7-bit index* |
| **Matched search size `K`** | 2-of-128 against argmin-over-128 | **Inverted the sign.** A minimum over C(128,2) = 8128 supports read against one over 128; at matched `K` the 0.076 Å gain is a **0.077 Å penalty** |
| **`A2` = the same operator written twice** | every built-chain arm-to-arm comparison | **Recalibrated the whole instrument.** Clouds agree to 3.6e-14, chain means at 0.03× MDE — but per-target |Δ| is 0.0134 mean / 0.0329 p90 / **0.2285 max**. **The 0.0107 Å floor is a MEAN floor** |
| **Pure rotation** (analytically exact, like reflection) | a reflection-invariance test whose 1e-6 threshold had fired at 3.673e-06 | **Inverted a registered falsifier.** Rotation gives **5.377e-06, larger** — the threshold was mis-set at an absolute value only a matched control could read |
| **Constant α-helix** (plausible, never uniform-on-the-torus) | every torsion channel | **Beat all of them on both bands** (−0.1722 at 3.17× MDE) — the Ramachandran channel's in-pool skill is a **constant-prior effect** |
| **Achiral twin** `\|X\|` | chiral functionals | Whatever the chiral channel does, its reflection-invariant shadow already does |
| **Deposited-ensemble spread** | the tail's difficulty | **Ruled out the reference as the explanation** — 0.20× MDE, and the tail's ensembles are *narrower* |
| **Random-feature control** | the predictability split | `y_perp` R² **+0.9403** against −0.0026; `y_along` −0.0051 against −0.0039 |
| **Magnitude-matched shrinkage** | the orthogonal-correction arm | **The perp arm is +0.6169 Å worse than its own magnitude-matched control** (3.25× MDE, 9W/117L) — magnitude is not the explanation |
| **Norm-matched shrinkage** | the projection operation itself | 0.26× MDE — projecting is indistinguishable from making the corrector smaller |
| **Uniform-weight ablation** (`:1116`, already in the code) | the p-weighted convex readout | The shipped ablation, used rather than rebuilt |
| **Shuffled-feature / shuffled-label** | the `â` combination arms | 0.13× and 0.71× MDE — NOT MEASURED |
| **Random permutation index map** | structure-aware indexing | **It is the best map** for a deployable partial readout |
| **Random-18 stratum null** (20,000 draws) | tail claims | Rules out *"an ordinary 18-subset"* — and **does not** rule out *"a stratum defined by the quantity measured on it"*, which is why `FAIL18` needed §12 |
| **Vertex feasibility assert** | the ORACLE convex solver | Caught 4 of 126 targets where both iterative solvers returned a point worse than a feasible vertex; **119W/4L → 119W/0L**, the four losses *were* the four failures |

**Two controls that changed a lane's own conclusion rather than confirming it** deserve naming:
the **shuffled-`B`** control (which demolished the reason for its lane's only positive) and the
**rotation** control (which inverted a registered falsifier that had fired). Both were run by the
lane that stood to lose from them.

---

## 14. The 126-target endpoint

[PENDING]

---

## 15. Statistical analysis

**MDE = 2.8016 × SE, per comparison.** Below **0.7×** is *not a result*; **0.7–1.0×** is *NOT
MEASURED*. Fold-clustered CIs on the pinned folds; `s24.stats_lib.compare` **refuses a verdict** if
`folds` is not passed.

**Multiplicity was tracked sprint-wide and written to as comparisons were emitted**, not
reconstructed at the end — the build item S30 carried forward. `s31/MULTIPLICITY.md` holds the
register; the running total is in the hundreds with a Bonferroni multiplier table computed up front
(k = 32 → 1.43×, k = 912 → 1.74×). **One lane's first total was wrong by 145 and was corrected in
place with the error stated.**

**The verifier** (`s31/s31_verify.py`) recomputes **129 numbers from artefacts, 0 mismatched, 0
missing, 0 flagged.** It **parses the sprint's own documents for every path they name** (28 found, 28
exist) rather than using a hand-kept list — *because a hand-kept list is exactly what fails*. It
checks both reporting bases, asserts `MDE == 2.8016 × SE` from the SEs, and carries forward S30's
sign-convention trap. **It caught two defects in its own author's work.**

### Three places the rules bound against the person applying them

1. **A fold CI that excludes zero does not rescue a sub-MDE effect.** The branch-carry arm came in at
   **0.44× MDE** with a fold CI excluding zero and 5/5 folds agreeing — *and was reported as NOT A
   RESULT*. The widening arm is the same shape at **0.42× MDE** with CI [+0.0028, +0.0557]. This is
   the documented sibling of the underpowered bug (`s24/stats_lib.py:145`), and **the MDE gate
   binds.**
2. **A null in the mean is not "no effect".** `P1`'s per-target spread is **10.7× the implementation
   noise null**, with 74 of 126 targets moving more than that null's p90, while the mean is −0.0112
   at 0.19×. The honest sentence is *"it reshuffles the answer everywhere and buys nothing."*
3. **The per-target floor is not the mean floor.** 0.0107 Å bounds *mean* built-chain claims; the
   per-target floor is **~0.03 Å with a 0.23 Å tail** (§1.2).

### Registered predictions, scored

Ten pre-registrations were committed before their first number. **Predictions that failed are
recorded as failures, including the coordinator's.**

| prediction | by | outcome |
|---|---|---|
| The orthogonal complement is noise (2:1 against E2/E3) | coordinator | **Right for the wrong reason** — it is 94% predictable and *actively harmful*, not noise |
| Long-range R² ≈ 0 (3:1) | coordinator | **Falsified** on the raw statistic (+0.1959) |
| A substantial fraction of `bestm128` survives (2:1) | coordinator | **Falsified** — the surviving fraction is *negative* |
| `LEG_torsion` is predominantly odd (4:1) | lane B | **Failed** — it is 70% even |
| `F-G2` does not fire | lane G (S30) | Held, and its *mechanism* was refuted by its own pre-check |
| The derived readout beats production ([−0.15, +0.10]) | lane A | **Missed at +0.1334**, on the side that says the direction fails |
| Disagreement count 20–45 | lane P | **Wrong — 66** |
| Arm E lies between D and F | lane P | **Falsified in premise and outcome** |
| `p*` does not improve on the VQE (2:1) | lane P | Held |
| A deployable prefix-`m` rule transfers | lane F | **Falsified**, three ways |

> **A pre-registration that only ever confirms is decoration.** Six of these went against the lane
> that wrote them, and three of those six are the coordinator's.

---

## 16. What was closed

Several of these are theorems. Where a direction closed by **derivation before compute was spent**,
that is noted — it is the charter's §24 ladder working as intended, and it is most of this table.

| direction | closed by | how |
|---|---|---|
| **The deployed CVaR objective as a quantum problem** | **theorem** | It is a **convex program with a closed-form global minimiser** pinned by one scalar. `run_cvar_vqe` is strictly worse in **126/126** targets at the deployed settings; `p*` costs 0.0012 s against 0.0985 s. Charter §11 closed for the deployed objective |
| **Any target-specific role for the quantum state** | **derivation + measurement** | `E = _zrank(...)` is a **target-independent constant** (§5). `H(p*)` has sd **3.1e-04** across 126 targets. The stage answers a global hyperparameter question |
| **A non-diagonal Hamiltonian** | **three independent theorems** | CVaR needs a per-shot eigenvalue; the forced operator is mean-field and quartic in ψ; and a candidate-index register's Hilbert dimension *is* the candidate count, so no operator on it can be classically hard |
| **ADAPT-VQE / qubit-ADAPT** | **theorem** | Its selection rule `\|⟨ψ\|[H,A]\|ψ⟩\|` presumes the cost is `⟨H⟩`, a linear functional; **CVaR is not the expectation of any observable**, so the criterion is *undefined*, not merely unhelpful |
| **The free-energy stage (§7A)** | **derivation, no compute** | `S` has no target argument, and the shipped potential is **reflection-invariant** (1340 torsion phases at distance 0.000e+00 from {0,π}), so `F`, `E`, `S` and every `dF/dT` are distance-map functions by G1 |
| **The elastic-network / normal-mode / landscape-curvature family** | **derivation, no compute** | **Corollary B1′** — an ANM/GNM Hessian is built from pairwise distances, so its spectrum and log-determinant (*which is the harmonic entropy*) are distance-map functions |
| **The backbone-torsion channel (§7B)** | measurement | The chiral escape from G1 is real (odd half 2.78× MDE on coarse triage) **and only works on a problem the pipeline does not have** (0.28× in-pool). A **constant α-helix beats every torsion channel on both bands** |
| **Charter §14, the new-observable question** | **provenance** | 92.9% of the benchmark is NMR-determined; for **117/126** the deposited coordinates *are* a fit to the deposited restraints, so any NMR observable is **ORACLE through a different door**. The one genuine escape (VCD/ROA, chiral, works at 9–16 residues) has **zero measured spectra** for these targets |
| **Candidate-index redesign (§12)** | measurement | Gray coding is a **proven no-op**; the best map raises an ORACLE ceiling and **lowers** the deployable value; the best map for a deployable partial readout is the **random permutation**; no index bit carries more than **0.066 bits** about candidate quality |
| **The per-target prefix length `m`** | **matched control** | A per-target minimum over 128 **random subsets** reaches **149%** of the prefix family's gain — the prefix axis is *worse than an arbitrary 7-bit index*. And the transfer was already **FALSIFIED in S29** (−0.0018 global, +0.0079 leave-fold-out) |
| **The constrained-affine readout** | **by exhaustion of selection rules** | No knee by ceiling (the ridge path is smooth and monotone) and no answer by fit (out-of-fold R² 0.021, applied +0.0309 worse). **It cannot be selected by ceiling and cannot be selected by fit** |
| **Consensus as a quality estimate** | **identity** | The consensus criterion **is** the free half of the objective read at uniform weights — `corr(medoid criterion, B·1/D) = 0.9667` over all 126. Setting `â ∝ +consensus` **cancels** the term it was meant to complement |
| **My own opening hypothesis (E1)** | **exact identity** | `mu_hat = mu − y`, so the proposed common-mode estimator's error **is** the prior error it was meant to help predict. Circular |
| **The orthogonal-correction family (E2/E3)** | **ORACLE ceiling** | With perfect magnitude *and* a perfect direction, correcting only the component orthogonal to the common mode is worth **+0.0747 Å — harmful** (§20) |

---

## 17. What improved

[PENDING]

---

## 18. What did not improve

[PENDING]

---

## 19. ORACLE ceilings

Every number here is **ORACLE / NOT DEPLOYABLE**. They bound what an operator could reach with
perfect information; none is achievable, and several are vacuous for reasons worth stating.

| ceiling | basis | value | note |
|---|---|---|---|
| affine readout over 128 | CA cloud | **0.0000** | **vacuous by rank** — 127 weights against ~33 residual dimensions |
| convex readout over 128 | CA cloud | **1.7977** | support 6.54/128 *without* a sparsity penalty; 9.59× MDE, 119W/0L |
| convex readout over 128 | built chain | 1.8538 | S29's rows, quoted unchanged |
| argmin over 500 (9 bits) | built chain | 1.7078 | |
| argmin over 128 (7 bits) | built chain | 2.1435 | **1.067 Å of headroom below production** |
| per-target prefix length `m` | built chain | 2.9027 | **an order statistic — see §12** |
| best pool member, genuinely worst 18 | built chain | 2.5298 | the tail is selection-limited |
| branch choice among the 8 already computed | built chain | **−0.0938 vs production** | 2.92× MDE, 114W/0L; free to compute |

**Two of these deserve their caveat repeated in any quotation.**

**The affine ceiling is not a bound.** `rank(aff{W_x}) = 32.9 ≥ 3n−3`, so the hull spans the space
and the "ceiling" is an interpolation artefact. The source file's own instruction is *"do not quote
any affine-hull ceiling as a bound."*

**`bestm128 = 2.9027` is an order statistic, and the framing that was licensed is narrower than the
one that was used.** See §12. The *ceiling* interpretation stands — **an optimistically biased upper
bound is still a valid upper bound**, so "2.5 Å is unreachable through this architecture" is safe
and **safer than stated**. The *lead* interpretation — "0.308 Å available to a better selector" —
was priced at **−0.0018 global / +0.0079 leave-fold-out** and marked **FALSIFIED in S29**, and the
caveat was dropped in re-quotation.

---

## 20. The remaining information bottleneck

### 20.1 The requirement was stated backwards, and the correction is the sprint's most consequential result

S30 concluded — and this report's own contract rule 29 encoded — that the next channel must supply
an observable whose error is **incoherent** with the pool's common mode. **Measured with ORACLE
magnitude and a perfect direction, that is worth `+0.0747 Å`: harmful.** The decomposition, splitting
the ORACLE ideal correction `y` against the ORACLE common mode `mu` per target (built chain, all
**ORACLE / NOT DEPLOYABLE**):

```
correcting ONLY the along-mu half    -0.8102 A   3.29x MDE, 5/5 folds, 121W/5L   <- the WHOLE prize
correcting ONLY the orthogonal half  +0.0747 A   0.87x MDE, 57W/69L              <- worth less than nothing
along vs perp, paired                -0.8848 A   3.50x MDE, 5/5 folds, 118W/8L
```

And the reason is an **identity, not a fit**. Refitting the same ridge on each half:

```
y_PERP_mu    out-of-fold R2  +0.9403    (random-feature control -0.0026)
y_ALONG_mu   out-of-fold R2  -0.0051    (random-feature control -0.0039)
```

`mu_hat = mu − y` forces `y_perp = −P_perp(mu_hat; mu)` **exactly to 9.3e-15**, and with
`cos(mu_hat, mu) = 0.05` that makes `y_perp ≈ −mu_hat` at corr 0.971. **The predictable half is not
predicted — it is observed.** The unpredictable half is exactly the common mode that S30-L7 proved
non-identifiable from pool data at any K.

> ### What can be predicted is the component ORTHOGONAL to the common mode, and it is harmful.
> ### What would help is the common mode itself, and it is unpredictable.

**S30's slogan survives verbatim; its mechanism was inverted.** The measurement behind it (+0.0554
coherent against −0.2466 i.i.d. at matched R²) is **untouched and still stands** — it was explained
backwards, not measured wrongly.

**Controls, all registered before the numbers:** magnitude is not the explanation (the perfect
correction shrunk to the perp arm's exact norm gives −0.5422, so the perp arm is **+0.6169 worse
than its own magnitude-matched control**, 3.25× MDE, **9W/117L**); the projection is not the
mechanism (against a norm-matched shrinkage, 0.26× MDE); and projecting against the *true* `mu` buys
nothing over not projecting (0.18× / 0.08×).

**Sizing, ORACLE / NOT DEPLOYABLE:** the along-`mu` correction alone puts the built chain at
**2.4025 Å** — which would clear the charter's *ambitious* 2.50 target — and it is **4× concentrated
on the tail** (FAIL18 6.0195 → 3.7392, −2.2803, against −0.5651 on the other 108). It does **not**
beat a perfect prior: ALONG against FULL is −0.0345 at **0.51× MDE, 3/5 folds — a tie, NOT
MEASURED**, and it is reported as a tie.

### 20.2 The second bottleneck, which is a different object

Lane A's closing statement, from the exact readout identity:

> The readout's objective is exact and half of it is free; the conversion from any quality estimate
> to an endpoint is a tuning-free convex program; that program's ORACLE ceiling is **1.7977 Å**
> against production's **3.0483** (CA cloud). **What is missing is a per-candidate quality estimate
> with positive IN-BAND skill.**

Every native-free candidate now measured has in-band skill that is **zero or the wrong sign**:

```
DIS z-rank (the shipped cost)    rho global +0.1176    IN BAND  -0.0262
CONS (the medoid criterion)      rho global +0.4278    IN BAND  -0.2837
DIS + CONS, leave-fold-out       rho global +0.4281    IN BAND  -0.2827
```

**Consensus reaches double the crossing price on the global axis and the wrong sign in band** — its
`ρ_global = +0.428` is **entirely outlier detection**. *The binding axis is in-band `ρ`, and the
global crossing price of 0.211 is necessary-not-sufficient: it holds along an interpolation path and
must not be quoted as a target for a real feature.*

### 20.3 The two bottlenecks are different objects, and whether they are one requirement is open

`â` is *which candidate is better*. `mu` is *how all of them are wrong together*. **Both are things
the pool cannot tell you about itself**, and both would have to come from outside it. Whether they
are two faces of one requirement or two separate ones is **the sharpest question this sprint
produces, and it is stated as open rather than resolved by assertion.**

### 20.4 What is NOT the bottleneck, so the next sprint does not pay for it again

The circuit's expressivity (`p*` is free and closed-form, and substituting it is worth **−0.0112 Å at
0.19× MDE**); the optimiser (80 → 2000 iterations move the gap by nothing); the CVaR α (inactive by
construction on three of five folds); the ansatz depth; the register width (the 75 → 128 widening is
a **+0.0307 Å COST** at **0.42× MDE, NOT MEASURED** — paid only to fill the `2**7` register (`core/pipeline.py:757`) and pointless when `quantum = False`, which is production; *its fold CI [+0.0028, +0.0557] excludes zero while the effect is 0.42×, the exact shape `_verdict` was hardened against, so the MDE gate binds and the CI does not rescue it*); the index encoding; the readout class (exact, and its
native-free channel is **3.99%** of ORACLE); and the candidate set (the set effect is **+0.0049 Å at
0.14× MDE — NOT A RESULT**, while the weighting rule costs **+0.1102 Å at 2.68× MDE**).

---

## 21. Next-sprint recommendation

**One question, stated as a requirement rather than a direction.**

> ### Find a source of information about the POOL'S COMMON MODE — the shared component of the retrieved candidates' error — that does not come from the pool.

This is the inverse of what S30 and this sprint's own opening contract asked for, and the inversion
is measured (§20.1), not argued. Three things make it a well-posed target rather than a wish:

1. **It is priced.** Perfect knowledge of the along-`mu` component alone is worth **−0.8102 Å on the
   built chain (ORACLE / NOT DEPLOYABLE)**, taking the endpoint to **2.4025 Å** — past the charter's
   *ambitious* target — with **4× concentration on the tail**.
2. **The conversion is free and tuning-free.** The readout identity makes any quality estimate
   convert to an endpoint through a **convex program with no hyperparameters** (§3(ii)).
3. **The failure mode is named.** The common mode is **non-identifiable from pool data at any K**
   (S30-L7), so a channel that reads the pool — however cleverly — cannot supply it. *That is the
   one thing we now know for certain about where it must come from.*

**And the parallel requirement, which may or may not be the same one:** a per-candidate quality
estimate with **positive in-band skill**. Every native-free candidate measured has zero or negative
in-band skill (§20.2). **State which of the two any proposal addresses.**

### What NOT to spend the next sprint on

- **Any further readout, index, ansatz, optimiser, α-schedule or register-width work.** All are
  closed above, most by theorem, and §20.4 lists them so the case does not have to be re-made.
- **A deployable 128 → 512 widening**, on the prediction in §11 — three facts held by three
  different lanes compose to say it should come out **worse** than production, not flat.
- **Chiral single-structure functionals at 9–16 residues.** Empty here, twice. **The theorem is
  length-free and only the emptiness is length-dependent**, so this is worth retesting at 40+
  residues on an instrument that does not exist yet.

### Two engineering items that gate future measurement

- **The per-target built-chain floor is ~0.03 Å, not 0.0107** (§1.2). Any S32 claim below that on
  individual targets is indistinguishable from re-running the same code twice. **The branch-carry
  fix (B4) is available, improves the conditioning 1082×, and is a NULL for accuracy at 0.44× MDE**
  — adopt it deliberately at a sprint boundary as a *conditioning change*, never as an improvement,
  and note it moves canonical 3.2105 → 3.2050.
- **The verifier should assert each number against the value on the basis it names**, not merely
  check that a basis is named. An audit for an *unstated* basis does not catch a *misstated* one,
  and this sprint had exactly one of the latter — in its own headline.

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
