# S31 — Lane L (literature). Running file.

Standing role for the whole sprint. Everything here is tagged **THEOREM** (proved, in the
literature or derived here), **PUBLISHED EMPIRICAL** (someone else's measurement), or
**OURS** (measured on this instrument in this sprint). A transfer argument is stated for
every published result, *including where it fails*.

---

## L1.1 — THE DEPLOYED CVaR FREE ENERGY HAS A CLOSED-FORM GLOBAL MINIMISER

**Status: THEOREM (derived here from two standard results) + OURS (verified numerically
against the shipped code).**

The shipped objective is `core/quantum.py:993` `free_energy`:

```
F(p) = CVaR_alpha(E; p) - T * H(p) ,   p = p_theta = |psi_theta|^2  (exact statevector)
```

### The derivation

Rockafellar & Uryasev (*J. Risk* 2:21, 2000) give the variational form of CVaR. For the
**lower** tail — which is what `cvar_exact` computes, and I checked the code matches this
form including its partial-mass boundary term:

```
CVaR_alpha^low(E;p) = max_s { s - (1/alpha) * sum_i p_i (s - E_i)_+ }
```

Inside the max this is **affine in p**. A pointwise max of affine functions is convex, and
`+T*sum p log p` is convex, so

> **F is a CONVEX function of p on the simplex, strictly convex for T > 0.**

It is concave in `s`, the simplex is compact convex, so Sion's minimax applies and the
order swaps. The inner minimisation over `p` is the standard Gibbs variational problem,
giving

```
min_p F  =  max_s { s - T * log sum_i exp( (s - E_i)_+ / (alpha*T) ) }     [1-D, concave]

p*_i  proportional to  exp( (s* - E_i)_+ / (alpha*T) )
```

**Read what `p*` is.** `(s* - E_i)_+` is zero for every candidate at or above the VaR
level `s*`. So the optimal distribution is **uniform on the whole body of the pool and
exponentially tilted only below the VaR threshold** — a *hinged Gibbs* distribution. The
entire `2**n`-dimensional optimisation is pinned by **one scalar `s*`**, obtained by a
golden-section search on a 1-D concave function in microseconds.

This is **strictly stronger than T1**. T1 says the CVaR tail is a prefix of the induced
order. This says: with the entropy term that was added *specifically to break T1's
degeneracy*, the optimum is still determined by a single number — and now that number has
a closed form.

### OURS — the numerical check (`s31/lit_L/lit_L_dequant_check.py`)

`E` = shuffled z-rank ladder (what `_zrank` produces up to the pool's tie structure).
`F_dual` is the 1-D dual value; `F_mirror` is an **independent** mirror-descent optimum
over the full simplex; `F_vqe` is the shipped `run_cvar_vqe` at its deployed settings
(layers 3, iters 80, restarts 1, exact parameter-shift gradient).

```
  n alpha    T |    F_closed      F_dual     F_mirror        F_vqe |  VQE-CF gap   dualgap |  H_cf   H_vqe |   TV
  7  0.10 0.10 |   -2.187057   -2.187057    -2.187054    -2.161293 |   +0.025764  8.27e-11 |  4.685   4.548 | 0.246
  7  0.10 0.05 |   -1.952814   -1.952814    -1.952802    -1.934213 |   +0.018602  3.18e-09 |  4.685   4.393 | 0.314
  7  0.25 0.10 |   -2.146335   -2.146335    -2.146333    -2.111633 |   +0.034702  2.68e-09 |  4.402   4.253 | 0.315
  7  0.25 0.05 |   -1.929407   -1.929407    -1.929403    -1.904976 |   +0.024431  8.72e-10 |  4.272   4.083 | 0.391
  7  1.00 0.10 |   -1.862495   -1.862495    -1.862495    -1.628915 |   +0.233580  1.11e-15 |  2.310   0.671 | 0.901
  7  1.00 0.05 |   -1.762185   -1.762185    -1.762176    -1.597652 |   +0.164533  6.66e-16 |  1.626   0.588 | 0.962
  9  0.10 0.10 |   -2.328500   -2.328500    -2.328500    -2.267516 |   +0.060984  6.51e-11 |  6.056   5.537 | 0.439
  9  0.10 0.05 |   -2.026793   -2.026793    -2.026791    -2.000685 |   +0.026108  4.53e-10 |  6.002   5.460 | 0.385
  9  0.25 0.10 |   -2.286109   -2.286109    -2.286109    -2.210826 |   +0.075283  2.57e-10 |  5.770   5.062 | 0.499
  9  0.25 0.05 |   -2.000873   -2.000873    -2.000869    -1.959409 |   +0.041463  1.19e-10 |  5.618   4.850 | 0.530
  9  1.00 0.10 |   -2.001363   -2.001363    -2.001363    -1.764573 |   +0.236790  4.44e-16 |  3.693   1.181 | 0.854
  9  1.00 0.05 |   -1.832023   -1.832023    -1.832014    -1.751284 |   +0.080739  2.44e-15 |  3.001   0.930 | 0.799
```

- **Strong duality holds to 1e-9 … 1e-16 in all 12 cells** — the closed form *is* the
  global optimum, not an approximation.
- An independent mirror descent over the full simplex reproduces it to ~1e-5.
- **The circuit is strictly worse in every cell**, and its distribution is far away:
  total-variation distance **0.25 – 0.96** from `p*`.
- The closed form is always **more entropic** than what the circuit finds
  (e.g. n=9, alpha=1, T=0.1: `H* = 3.693` bits vs the circuit's `1.181`).

### OURS — it is an EXPRESSIVITY floor, not an optimiser floor (`s31/lit_L/lit_L_gap_budget.py`)

The obvious objection is "you did not optimise hard enough." Priced:

```
alpha     T layers  iters  rst |      F_vqe   F_closed       gap     TV   H_vqe    sec
 0.10  0.10      3     80    1 |  -2.161293  -2.187057 +0.025764  0.246   4.548    0.1
 0.10  0.10      3    400    1 |  -2.164745  -2.187057 +0.022311  0.247   4.521    0.6
 0.10  0.10      3   2000    1 |  -2.161113  -2.187057 +0.025944  0.246   4.567    2.5
 0.10  0.10      3    400    8 |  -2.178638  -2.187057 +0.008418  0.123   4.646    4.3
 0.10  0.10     12   2000    1 |  -2.182959  -2.187057 +0.004098  0.103   4.652   25.8
 0.10  0.10     12    400    8 |  -2.182904  -2.187057 +0.004153  0.105   4.651   43.6
 1.00  0.10      3     80    1 |  -1.628915  -1.862495 +0.233580  0.901   0.671    0.2
 1.00  0.10      3   2000    1 |  -1.629379  -1.862495 +0.233116  0.901   0.673    4.0
 1.00  0.10      3    400    8 |  -1.775067  -1.862495 +0.087428  0.582   1.249    6.2
 1.00  0.10     12    400    8 |  -1.844926  -1.862495 +0.017569  0.098   2.481   59.5
```

**Iterations buy nothing.** 80 -> 2000 at fixed depth leaves the gap flat (0.0258 ->
0.0259; 0.2336 -> 0.2331). What closes it is **depth** (expressivity) and **restarts**
(landscape multimodality). The optimiser has converged; the ansatz cannot reach `p*`.
21 parameters (`n*layers` = 7*3) for a 127-dimensional simplex. Even at 60 s of compute
the gap and a TV of ~0.10 remain, against **microseconds** for the closed form.

### The correction this forces on our own shipped docstring — L3, turned inward

`core/quantum.py:997-1002` says the entropy collapse "**is a property of CVaR, not of the
optimiser**."

- At **T = 0** that is correct, and I can now give the one-line proof: any `p` with mass
  `>= alpha` on the argmin attains `CVaR = E_min` exactly, so the minimiser set is a
  **positive-volume flat face** of the simplex and the optimiser's path — not the
  objective — picks the point inside it. Anything the readout reads there is an Adam
  artefact. **THEOREM.**
- At the **deployed T = 0.1** it is **wrong**, and it is wrong in a way that matters. The
  objective's own optimum at `alpha = 1, T = 0.1` carries **2.310 bits** (n=7); the
  circuit delivers **0.671**. The collapse is a property of **the ansatz**, not of CVaR
  and not of the optimiser.

This also re-prices the sprint's inherited reading of the `alpha` effect. `pipeline.py:825`
reports the state entropy going **0.076 -> 6.36 bits** from `alpha=1.0` to `alpha=0.1` at
`T=0.1`. On my synthetic ladder the circuit shows 0.671 -> 4.548 (3.88 bits) while the
**objective** asks for 2.310 -> 4.685 (2.375 bits). So roughly **60% of the observed
alpha-effect on state entropy is a real property of the objective and ~40% is the ansatz's
failure amplifying it.** *Caveat carried with the number: this is a shuffled z-rank ladder
at n=7, not our measured pool `E`; the split is indicative, not a measurement on the
instrument.*

### What this closes, and the one thing it does not

It **closes charter §11 for the deployed objective**: the CVaR-VQE as shipped is a *lossy
approximate solver for a convex program with an analytic solution*. There is no quantum
mechanism in it — not interference, not spectrum, not entanglement. Whatever it
contributes, it contributes by **failing** to optimise.

It does **not** close the possibility that the failure is *useful*. The circuit's
reachable set is a 21-parameter manifold inside the simplex; that constraint is an
inductive bias, and this project has a standing result (`concentration-is-wrong-when-
discrimination-binds`) that optimising a bad objective harder makes things worse. So:

> **The decisive experiment is a one-line substitution.** Replace `run_cvar_vqe`'s returned
> `p` with the closed-form `p*` and rerun the endpoint on tuning126. It costs nothing.
> - endpoint **improves** -> the quantum layer is replaced by a softmax and a root-find;
> - endpoint **worsens** -> the circuit's *inability* to optimise is the active ingredient,
>   which is S20's law arriving from the objective side and is a genuine, publishable
>   negative about the architecture.
>
> There is no third outcome, and both are results. **Caveat: the readout is
> `consensus_medoid(block, p)`, so a large TV in `p` need not move the medoid.** The
> substitution must be scored at the endpoint, not on `p`.

**Honest limit of the escape hatch.** If the inductive-bias story wins, the bias is a
**Born machine** `p = |psi|^2` from a shallow 1-D RY+CNOT circuit, which is an **MPS Born
machine** with small bond dimension — a classical tensor-network generative model since
Han et al., *PRX* 8:031012 (2018). So even the escape hatch is classical; it would be a
result about a *useful regulariser*, never about quantum advantage.

---

## L1.2 — ALPHA IS A HOMOTOPY PARAMETER IN THE LITERATURE AND A SHAPE PARAMETER HERE

**Kolotouros & Wallden, "Evolving objective function for improved variational quantum
optimization," *Phys. Rev. Research* 4, 023225 (2022)** (arXiv:2105.11766). Ascending-CVaR.

Schedules, verbatim:
```
linear:   alpha_{t+1} = alpha_t + lambda ,   lambda in [0.025, 0.045]
sigmoid:  alpha_t = 1 / (1 + e^(5 - lambda*t)) ,   lambda in [0.3, 0.4]
```
**PUBLISHED EMPIRICAL** (VQE, hardware-efficient ansatz, 15–20 qubits): hard number
partitioning — Ascending-CVaR 95% success, constant CVaR (alpha=0.1) 58%, plain
expectation 0%. Portfolio optimisation (16–20 assets) — 100% vs 1% for expectation.
Max-Cut (15–19 vertices) — 96 vs 53 successful instances.

Their two propositions are what matter:
- **Prop. 1: all `CVaR_alpha` with `alpha <= kappa` share the same ground state.**
- **Prop. 3: a local minimum of `CVaR_{alpha1}` need not be one of `CVaR_{alpha2}`.**

So in the literature **alpha changes the landscape, not the answer** — it is a homotopy
parameter, and the schedule runs **upward**, ending at the plain expectation value.

### The transfer argument, and where it fails

**It does not transfer as-is, and the direction is the reason.** Prop. 1 holds at `T = 0`.
Our deployed objective has `T = 0.1`, and the entropy term makes the optimum genuinely
alpha-dependent: `H*(alpha=0.1) = 4.685` bits vs `H*(alpha=1.0) = 2.310` bits at n=7. So
here **alpha is a distribution-shape knob, not a landscape knob** — and the literature's
schedule ends at `alpha = 1`, which for us is the **narrowest** ensemble, i.e. the
direction the consensus readout least wants.

**The usable translation.** The closed form gives `H*(alpha, T)` analytically, so one can
schedule `alpha` upward *while solving for `T(alpha)` that holds `H*` fixed*. That keeps
Ascending-CVaR's landscape benefit and removes its shape side-effect. This is cheap and
well-founded; it is the only form in which I would recommend the paper's method here.
**Note it is an optimisation-landscape fix, and L1.1 says the landscape is not what limits
us** — so I rank it below the `p*` substitution.

---

## L1.3 — A NON-DIAGONAL HAMILTONIAN IS NOT A CVaR-VQE AT ALL, AND A CANDIDATE-INDEX
## REGISTER CANNOT CARRY A HARD HAMILTONIAN

**Status: THEOREM + literature.** Two independent obstructions to the lane-A plan
`H = diag(zrank) - lambda * W(block)`.

**(a) CVaR is defined only for diagonal `H`.** The CVaR objective needs an energy per
*shot*. For diagonal `H` every bitstring is an eigenstate and each shot carries a definite
eigenvalue; Barkoutsos et al. (*Quantum* 4, 256, 2020) introduce it in exactly that
setting. Off-diagonal elements require coherent, multi-basis estimation and **cannot be
assigned to a single measurement outcome**. So `diag(zrank) - lambda*W` is not a CVaR-VQE
objective — you would have to either pre-diagonalise `W` (a classical eigendecomposition,
which puts the quantum layer *downstream* of the solve) or redefine "CVaR" over some other
random variable, at which point the literature's results no longer apply and neither does
the name.

**(b) Dimension counting kills it at our register size.** The candidate-index encoding
puts `N = 2**n` candidates in an `n`-qubit register, so **the Hilbert space dimension
equals the number of candidates**. *Any* operator on that register — diagonal or not — is
a `128 x 128` or `512 x 512` matrix. Its full spectrum is a microsecond `eigh`. **No
Hamiltonian design on a candidate-index register can be classically hard.** This is not a
statement about our ansatz or our budget; it is arithmetic.

The only encoding with a Hilbert space that is not a small matrix is the **assignment**
encoding — one qubit per candidate, include/exclude — where `n` qubits carry `2**n`
*subsets*. That is what the clustering literature actually uses: **Hamiltonian
formulations of centroid-based clustering** (*Front. Phys.* 2025, arXiv:2502.06542)
formulates k-means/k-medoids as Ising/QUBO with **one qubit per data point**, total qubits
= `N`, tested to 175 points on D-Wave Advantage. Our 128 candidates would need **128
qubits, not 7**. And the authors state no quantum advantage is demonstrated — they frame
it as "opening a path toward" one, with intercluster variants suffering chain breaks that
"produce random solutions."

**Transfer:** the assignment encoding is the only one where the pairwise structure the
readout consumes becomes a genuine `2**N` optimisation — and at `N = 128` the resulting
QUBO is one the project has already priced classically (S30 §9.3: Frank–Wolfe on the
simplex "strictly upper-bounds any circuit on this objective"). So the pairwise Hamiltonian
is either **too small to be hard** (index encoding) or **too big to run and already solved
classically** (assignment encoding). I do not see a middle.

---

## L1.4 — OUR MEASURED GRADIENT LAW IS THE PREDICTED GLOBAL-COST BARREN PLATEAU, AND IT
## PRICES THE 128 -> 512 WIDENING

**Cerezo, Sone, Volkoff, Cincio & Coles, "Cost function dependent barren plateaus in
shallow parametrized quantum circuits," *Nat. Commun.* 12, 1791 (2021).** **THEOREM:**
cost functions defined by **global** observables have exponentially vanishing gradients
**even at O(1) depth**; **local** cost functions vanish at worst polynomially provided
depth is `O(log n)`.

CVaR is maximally global — it depends on the rank ordering of all `2**n` basis energies,
and (L1.1) it is not even a linear functional of the state, so it is not the expectation
of *any* observable.

Our own law `Var[dF/dtheta] = r_stable/D^2` with `r_stable ~ 16D` gives

```
Var ~ 16/D = 16 * 2^(-n)
```

— exponential in the qubit count, which is **exactly the regime Cerezo et al. predict**,
and their theorem says **depth will not fix it**. That is an analytic reason for the
charter's "do not assume deeper circuits automatically help," at no experimental cost.

**The quantitative brake on lane C's widening.** The law is a factor of **4 per two added
qubits**:

```
n = 7  (128)   Var ~ 16/128  = 0.125      deployed
n = 9  (512)   Var ~ 16/512  = 0.031      4x worse  -- affordable
n = 12 (4096)  Var ~ 16/4096 = 0.0039     32x worse -- marginal
```

So 128 -> 512 is **not** blocked by trainability — but it buys a 4x smaller gradient
signal, and anything past ~12 qubits is in the exponential wall. **Caveat, stated because
it cuts the other way:** my L1.1 budget sweep shows that at n = 7–9 we are *not* yet
gradient-limited — depth and restarts still move the objective a lot. The barren plateau is
a bound on where this can go, not a diagnosis of where it is stuck now.

---

## L2 — THE NEW-OBSERVABLE QUESTION (in progress)

Theorem G1 leaves **chiral** and **many-body** observables open. First pass on the
chiral side:

**Vibrational circular dichroism (VCD) and Raman optical activity (ROA)** are genuinely
chirality-sensitive — they measure differential response to circularly polarised light, so
they are *not* functions of the achiral distance map and therefore sit outside G1. Unlike
most of the structural-biology literature this project has had to discount, they work **in
our length band**: VCD is "sensitive to short-range order" and discriminates beta-sheet
from helices and disorder, and ROA has been applied to model peptides of ~6 residues
(*Chem. Rev.* 2020, 120:3381 and refs.).

**But the likely wall is data availability, not physics** — the same wall as chemical
shifts (`torsion-restraints-reach-the-target`: shifts exist for only 54/126 targets, and
ORACLE-perfect torsions on all of them still leave the instrument at 2.021 A). A
forward-computable observable with no measured target values is a *re-reading of the
structure*, not a new channel. Checking coverage next; **I will not report VCD/ROA as a
live channel until I can state how many of the 126 targets have measured spectra.**

---


---

## L1.5 — ADAPT-VQE / qubit-ADAPT: an efficient way to approximate something we can compute exactly

**qubit-ADAPT-VQE** (Tang et al., arXiv:1911.10205) and **ADAPT-QAOA** (Zhu et al.,
arXiv:2005.10258) grow the ansatz by greedily appending the pool operator with the largest
objective gradient. The charter asks whether a problem-inspired pool exists for a
*selection* problem rather than a chemistry one. Two reasons to rank this last:

1. **ADAPT is an expressivity fix, and L1.1 says our expressivity gap has a free
   alternative.** The measured gap between the shipped circuit and `p*` is an expressivity
   floor; ADAPT would climb toward `p*` with fewer parameters. But `p*` is already
   available in microseconds from a 1-D root-find. **ADAPT-VQE would be an efficient way
   to approximate a quantity with a closed form.**
2. **ADAPT's selection rule is not defined for a CVaR objective.** The criterion is
   `|<psi|[H, A]|psi>|`, which presumes the cost is `<H>` — a *linear* functional of the
   state. CVaR is not the expectation of any observable (L1.1: it is a max over `s` of
   affine functionals, i.e. convex and piecewise-linear in `p`, not linear in `rho`). So
   the pool-selection step would have to be redefined before ADAPT could run here, and the
   published guarantees would not carry over.

Transfer verdict: **not worth a lane this sprint.** It would be worth revisiting only if
the `p*` substitution shows the circuit's constrained reachable set is the active
ingredient — in which case the object of interest is the *constraint*, and ADAPT removes
constraints, i.e. points the wrong way.

---

## L2 — THE NEW-OBSERVABLE QUESTION: THE BINDING CONSTRAINT IS PROVENANCE, NOT GEOMETRY

**Status: OURS (measured on the benchmark) + THEOREM (data-processing inequality).**

### The measurement nobody had made: what determined our reference coordinates?

`s31/lit_L/lit_L_expmethod.py`, RCSB GraphQL over all 126 `tuning126` targets:

```
  115  ( 91.3%)  SOLUTION NMR
    5  (  4.0%)  ELECTRON CRYSTALLOGRAPHY
    4  (  3.2%)  X-RAY DIFFRACTION
    2  (  1.6%)  SOLID-STATE NMR
                 -> NMR-determined 117/126 = 92.9%,  X-ray 4/126 = 3.2%
```

### What that does to the question

The charter asks for an observable that is not sequence, not a pool re-reading, not a
distance-map scalar. G1 answers the *geometry* half. The 92.9% answers the half that
actually binds:

1. **Any observable computed at inference from `(sequence, pool)` adds no information,
   whatever equivalence class it occupies.** Data-processing inequality. A new observable
   of that kind can only be a better *estimator*, never a new channel — and the project has
   already measured that ceiling (`in-band-ordering-is-per-target`: 0.600 across targets
   against the 0.638 needed for 2.0 Å). **G1 is therefore not the binding constraint;
   provenance is.**
2. **Any NMR observable of these targets is the data that DETERMINED the reference.** For
   117/126 entries the deposited coordinates are a fit to deposited NOEs, J-couplings and
   torsion restraints. Feeding those back in is **ORACLE through a different door**, not a
   new channel.
3. **This re-prices the chemical-shift direction retroactively.** `torsion-restraints-
   reach-the-target` reports shifts for 54/126 and ORACLE-perfect torsions still landing at
   2.021 Å. TALOS+/TALOS-N derived dihedral restraints from chemical shifts are **standard
   practice** in NMR peptide structure determination, so on NMR-determined targets
   "ORACLE-perfect torsions" was close to a tautology — the torsions that reproduce the
   deposited coordinates are substantially the ones the restraints imposed. *Stated as
   standard practice, NOT verified per-target here.* The direction is already closed; the
   point is that the **reason** should be recorded correctly.

### The only genuine escape, and why it is empty here

A new channel must be **measured on the target molecule and not used in its structure
determination.** Chirality-sensitive vibrational spectroscopy is the best physical
candidate: **VCD and ROA** are differential responses to circularly polarised light, so
they sit outside G1's achiral class by construction, and unlike most structural-biology
methods they work **in our 9–16 band** (VCD is "sensitive to short-range order" and
separates 3-10/alpha helices, beta-sheet and disorder; ROA has been run on ~6-residue
peptides — Keiderling, *Chem. Rev.* 2020, 120(7):3381-3419).

**But there is no repository of measured VCD/ROA spectra keyed to PDB entries.** I searched
for one and found only method and computational-protocol papers. So for these 126 targets
the measured-spectrum count is, to the best of my search, **zero** — against 54/126 for
chemical shifts, which are themselves disqualified by (2).

> **L2 verdict: the new-observable question is CLOSED for this benchmark, and it is closed
> by data availability rather than by physics.** The physics leaves chiral and many-body
> channels open; the benchmark supplies no measurement to put in them. This is the same
> wall as `no-fresh-benchmark-exists`, arriving from the observable side.

---

## L2.1 — IS THE FAIL18 TAIL REFERENCE NOISE? NO. A CHEAP NULL THAT PROTECTS THE SPRINT

The 92.9% NMR finding raises an obvious and alarming follow-up: the manifest's reference is
**"deposited coordinates, model 1"**, an arbitrary member of an ensemble. If the deposited
models disagree, the endpoint is partly scoring against noise.

**ORACLE DIAGNOSTIC of the BENCHMARK, never an inference-time signal.**
`s31/lit_L/lit_L_ensemble_spread.py` pulled every deposited model for all 126 targets and
computed pairwise Kabsch CA-RMSD; 111/126 resolved.

```
mean pairwise CA-RMSD BETWEEN DEPOSITED MODELS : mean 1.0823  median 0.9929  max 4.2648
model 1 -> ensemble medoid                     : mean 0.6965  median 0.4565  max 4.2943
spread > 1.0 A : 55/111        spread > 2.0 A : 15/111        spread > 3.0 A : 2/111
widest: 3BTB 4.265, 6CEJ 4.200, 6GIJ 2.958, 6EY3 2.674, 2MIG 2.629
```

**The reference really is uncertain.** Model 1 sits 0.70 Å from its own ensemble medoid on
average, and 4.29 Å away for 3BTB.

**And it does not explain the tail.** `s31/lit_L/lit_L_tail_vs_spread.py`, production arm,
n = 111 matched:

```
corr(production RMSD, ensemble spread) = +0.1118   95% CI [-0.0762, +0.2921]   spans zero
corr(production RMSD, model1->medoid)  = +0.1341   95% CI [-0.0536, +0.3127]   spans zero

WORST 18 by production RMSD : mean RMSD 6.0291   mean ensemble spread 0.9672
the other 93                : mean RMSD 2.6842   mean ensemble spread 1.1046
difference tail-minus-rest  : -0.1374   SE 0.2421   MDE 0.6783   -> 0.20x MDE  = NULL
```

The tail's targets have, if anything, *slightly narrower* deposited ensembles. **FAIL18 is
real failure against a reference no worse determined than any other target's.** The same
null holds in all four other arms (0.08x–0.90x MDE, every CI spanning zero).

**Why I am reporting a null this prominently.** Stopping at "the reference is uncertain by
1.08 Å" would have been quotable, alarming and wrong, and it would have redirected the
sprint's tail work. This is `control-at-the-decisive-step`: the cheap null sits exactly
where the wrong answer would first have become quotable.

**What survives, as an interpretive caveat only.** The endpoint carries a reference term of
~0.70 Å that is **uniform across targets** — a perfect predictor aiming at the ensemble
medoid would still score ~0.70 Å against model 1. Combined in quadrature that is worth
about `sqrt(3.21^2 - 0.70^2) = 3.13` against 3.21, i.e. **~0.08 Å at the current endpoint**
— small now, and it would matter far more near 1 Å. Because it is uniform, it **cancels in
every arm-to-arm delta**, which is the project's actual currency.

> **This is NOT a reason to re-score against the medoid.** `benchmark-and-folds-must-be-
> pinned`: the reference is part of the sealed instrument and changing it mid-sprint would
> invalidate every comparison. Report it as an uncertainty on the absolute number; change
> nothing.

**A small integrity note for lane D.** `results/summary/results.csv` gives the production
arm mean as **3.2126** over n = 126, against the charter's **3.2105** — a 0.0021
discrepancy. I flag it rather than reconcile it because the charter already lists **the
unpinned projection seed** as an open defect gating every sub-0.01 Å claim, and this is
consistent with exactly that. Lane D owns it.

---


---

## L2.2 — THE COLLAPSE TO ONE SOURCE IS TRUE, AND AS AN INFORMATION STATEMENT IT IS VACUOUS

The coordinator asked me to attack this rather than agree with it. **The premise survives; the
conclusion does not.**

### The premise is right, with one contamination caveat

`pool = f(sequence; library)` with a universal library is processing, so `I(N; pool) <=
I(N; sequence)`. A codebook indexed by a sequence-derived key carries no target argument —
the coordinator's reading is correct, and I could not break it. Enumerating sources is
closed: no field, energy, graph statistic or free energy is a second source.

**One real caveat**: the library is *not* independent of the native for every target.
`containment-threshold-is-at-the-null` records **4/126 targets carrying a verbatim copy in
their own distogram's training set**. For those the Markov chain `N -> S -> pool` is broken
and the DPI does not apply. That is **leakage, not a channel** — it argues for excluding
those 4, never for counting the library as a source.

### The conclusion does not follow, and the reason is Anfinsen

By the thermodynamic hypothesis the native is a function of the sequence, so

```
I(N ; sequence) = H(N)
```

— **the sequence already contains all of it.** The DPI bound therefore reads *"the pool
contains at most everything,"* which is true and empty. **It places no constraint whatever
on achievable accuracy.**

The cleanest way to see that the strong reading is false:

> **If "all target-specific information is the sequence" implied a ceiling, AlphaFold would
> be impossible.** A different computation on the same single source extracts far more than
> ours does. Information-source counting cannot be what limits us, because we are nowhere
> near the source's content.

**And this project has its own internal counterexample.** `esm-adds-nothing-for-short-
peptides` was OVERTURNED: on SELECTION, ESM buys **0.288 Å over one-hot (p = 0.005,
n = 126)**. That is a *pure computation on the sequence* beating a weaker computation on the
same sequence, with no new source and no new measurement. Under the strong reading that gain
could not exist. It does.

### What the argument actually licenses

Not *"you need a measurement of the molecule, not a computation."* That conclusion is
unlicensed, and adopting it would wrongly close the two escapes the field actually used.
What is licensed is the narrower and still very useful:

> **No re-combination of the objects we currently hold can help. The remaining moves are a
> better ESTIMATOR or a better LIBRARY — not a further source.**

### The distinction the collapse invites us to lose, and must not

The library plays **two different roles**, and the enumeration conflates them:

| role | verdict |
|---|---|
| **as an information source** | it is a constant, carries no target argument — **not a source**. Coordinator is right |
| **as an architectural ceiling** | **binding, and measured**: ORACLE distances still give only ~1.95–2.0 Å through this library (`distance-prior-is-the-ceiling`, `sequence-signal-is-the-ceiling`) |

**"Not a source" does not mean "not a constraint."** The ~2.0 Å figure is a *library*
statement, not an information statement, and it is the real ceiling — it survives perfect
information, which is precisely what an informational argument can never produce.

So the two measured constraints that actually bind are both **estimator/architecture**
facts, neither of them informational:

1. in-band ordering **0.600 across targets** against the **0.638** needed for 2.0 Å;
2. ORACLE distances capping at **~1.95–2.0 Å through this library**.

**Proposed central statement, replacing the vacuous form:**

> All target-specific information is the sequence, so enumerating sources is closed — but
> that is a statement about *sources*, not about *ceilings*. Because Anfinsen makes the
> sequence informationally complete, the bound is vacuous and cannot limit accuracy. What
> limits accuracy is measured and architectural: our estimator reaches 0.600 in-band against
> 0.638 needed, and this library caps even ORACLE distances at ~2.0 Å. The open moves are a
> better estimator or a better library. A measurement of the molecule is *one* escape, not
> the only one.

---

## L3.1 — AUDIT OF `bestm128 = 2.9027 Å`: THE TRANSFER ARM EXISTS, AND IT SAYS ZERO

Requested by the coordinator on the suspicion that a per-target ORACLE minimum over 128
prefix lengths is `grid-oracles-are-order-statistics` and was never transfer-tested.

**It was transfer-tested, in the same entry that produced it.** `s29/LEDGER.md:2777`,
S29-L30, in its own heading:

> *"...AND ITS TRANSFERABLE PART IS ZERO — THE ORACLE GLOBAL PREFIX IS m = 72 (WORTH
> −0.0018 Å, i.e. THE SHIPPED 75) AND THE LEAVE-FOLD-OUT PREFIX IS +0.0079 Å WORSE THAN
> PRODUCTION AT 0.30× MDE"*

The 33-arm ladder (`s29/LEDGER.md:4288-4290`) carries a dedicated **`transferable`** column,
and all three prefix-m rungs read `global m = 72 worth -0.0018; LFO +0.0079`.
`s29/REPORT_S29.md:252` lists *"A transferable prefix length m"* in its hypothesis table with
verdict **FALSIFIED**.

So the per-target ORACLE gain of **−0.3079 Å is ~99.4% order statistics**, and S29 both
measured that and said so. `grid-oracles-are-order-statistics` was **honoured here, not
violated** — this is one of the rungs that has a transfer arm.

### The number does not deflate. Its USE does.

This is the part worth being precise about, because the intuition runs the wrong way:

> **Order-statistics inflation makes an ORACLE number optimistically biased — and an
> optimistically biased UPPER BOUND is still a valid upper bound.**

So the inference S29 actually drew — *"2.5 Å is unreachable through this architecture"* — is
**safe, and in fact safer than stated**: the true ceiling is *worse* than 2.9027, which
strengthens unreachability rather than weakening it. **The "architectural ceiling" framing
was licensed.**

What is **not** licensed is reading 2.9027 as *"a −0.3079 Å lead available to a better
selector."* S29 priced that at **−0.0018 Å global / +0.0079 Å LFO** and marked it
**FALSIFIED**. If this sprint is treating "close the gap to 2.9027" as its most promising
remaining lead, that lead was measured to be worth **0.6% of its face value two sprints ago**
and the caveat was dropped in re-quotation — exactly the L3 pattern.

### The one thing genuinely still open, for lane F

The transfer numbers (−0.0018 / +0.0079) were measured on the **POINT CLOUD** — S29-L30's own
arm is 2.7605 against production 3.0483, both cloud. **2.9027 is the BUILT CHAIN.** The
cloud→chain price on that rung is **+0.1422**. So *"the transferable part is zero"* is being
quoted **across bases**. It is very likely fine — the transfer is ~0 on the cloud and the
chain price is near-uniform — but it is **not literally measured on the chain**, and the
charter's basis discipline says not to transfer across bases silently.

> **Recommendation to lane F: do not re-derive the ladder — S29 already did. Pin the one
> open thing: run the global-m and leave-fold-out-m arms ON THE BUILT CHAIN.** That is a
> two-arm measurement, not a sprint, and it either closes the last gap in a two-sprint-old
> claim or finds the only place it could move.

---

*(file continues; entries appended as the sprint runs)*
