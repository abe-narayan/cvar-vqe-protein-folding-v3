# SPRINT 21 — WORKSTREAM B — CVaR-VQE MECHANICS

Pre-registration: `s21/PREREG_B.md`, falsifier-first, written before any arm of this lane ran and
**unedited**. Code: `s21/qb3_lib.py` (charts, ansatz zoo, Fisher/QNG, exact bond dimension), `s21/qb3_run.py`
(runner), `s21/qb3_report.py` (analysis), `s21/qb3_validity.py` (the validity guard on §4),
drivers `s21/run_b.sh` and `s21/run_b2.sh`.
Artefacts, `s21/results/`, ten of them:

    qb3_config          the full declared configuration                       _COMPLETE
    qb3_machinery_gate  validity gates on this lane's own new arithmetic       _COMPLETE
    qb3_gate            F-E5 harness gate, 20 targets, 35,840 comparisons      _COMPLETE
    qb3_gauge_exact     the PREMISE of F-E1: 8 charts identical to 7.8e-16     _COMPLETE
    qb3_a               BLOCK A, ansatz zoo                    DIST+LEG done, AMBc  __PARTIAL_
    qb3_q               BLOCK Q, optimiser battery             DIST+LEG done, AMBc  __PARTIAL_
    qb3_s               BLOCK S, shots x iterations            DIST+LEG done, AMBc  __PARTIAL_
    qb3_concentration   latent concentration vs outcome                        _COMPLETE
    qb3_validity        the validity vector on the readout claim               _COMPLETE
    qb3_ci_fold         fold-clustered CIs on 24 dispositive contrasts         _COMPLETE

Every `__PARTIAL_` names the exact missing cells against the FULL declared configuration, not
against the subset the call happened to run.

Salt `s21qb`. Seeding `s15/seed.py::stable_rng` only. The sealed 60-target benchmark and
`results/benchmark_manifest.json` were not read, probed or derived from; `dev24` was not run.
Target set is `s20.qb2_lib.subset(20)` — the same 20 tuning targets Sprint 20's lane used, so every
arm here pairs against a Sprint-20 arm on the same target.

**Basis, stated once and true of every number below:** built chain from continuous torsions via
`core.project.build_ca_exact`; **readout is the argmin over everything evaluated** (`Field.best_z`),
which is `core/quantum.py`'s `best_seen` readout and *not* a tail mean — per the BRIEF's own
mid-sprint correction. No arm here uses an averaged readout. No point-cloud number appears.

**4 seeds on every variational arm.** Unit of analysis is the target; seeds are averaged *within* a
target before any test, so seeds buy precision and never inflate n.

### THE PRECISION OF THIS LANE, STATED BEFORE ITS RESULTS — and a pre-registration I mis-specified

**The BRIEF's 0.084 Å MDE belongs to the n=126 instrument. This lane runs n=20.** MDE is a property
of a **comparison**, not of an instrument: `MDE = 2.8016 · SE`, and `SE` depends on the paired sd of
the specific contrast. Measured from the paired standard deviations these arms actually produced:

    Ca-RMSD contrasts, Block A     MDE  0.226 - 0.378 A
    Ca-RMSD contrasts, Block Q     MDE  0.244 - 0.452 A
    pooled-arm Ca-RMSD contrasts   MDE  0.200 - 0.349 A
    CVaR TRAINING LOSS contrasts   MDE  0.103 - 0.312   (observed effects 0.214 - 0.524)

> **The objective column of this lane is adequately powered and the RMSD column is not.** Every
> RMSD null below is an exclusion of effects **larger than ~0.2-0.4 Å**, not of effects larger than
> 0.084 Å.
>
> **My own F-A3 is therefore MIS-SPECIFIED**: it says "if no ansatz beats `mps_L2` by ≥ 0.084 Å with
> a CI excluding zero, ansatz expressivity is closed." An n=20 design cannot deliver that test.
> Per BRIEF §8 the pre-registration **stays unedited**, the mis-specification is declared here, and
> **the regime 0.084-0.30 Å is recorded OPEN** rather than reported as closed. What the data *do*
> support is stated at that resolution and no finer.

The one place pooling recovers usable precision is the joint contrast — average every arm **within**
a target, then test, so arm-level noise is averaged rather than stacked and the unit stays the
target:

| pooled contrast (DIST) | difference | CI | MDE |
|---|---|---|---|
| mean of 8 ansatz − `best_of_N` | +0.036 | [−0.133, +0.195] | 0.236 |
| mean of 7 ansatz − `mps_L2` | −0.012 | [−0.148, +0.130] | 0.200 |
| χ ∈ {1,2} − χ = 16 | +0.022 | [−0.163, +0.205] | 0.273 |
| mean of 5 optimisers − `best_of_N` | +0.091 | [−0.083, +0.283] | 0.268 |
| **mean of 5 optimisers − the UNTRAINED circuit** | **+0.003** | **[−0.136, +0.144]** | 0.208 |
| QNG family − `adam` | +0.058 | [−0.185, +0.281] | 0.349 |

**That pooling barely moves the MDE is itself a measurement: the variance here is target-level, not
arm-level.** Which arm you run matters far less than which target you drew.

*(This lane reached the same conclusion independently of, and agrees with, the project-level finding
that **MDE is per-comparison, not per-instrument** — the 0.084 Å constant is a single contrast's
value on a different instrument and must not be carried as a universal threshold. **SE is reported
beside every mean here via the CI**, and each block's own MDE is stated above.)*

---

## 0. THE HEADLINE, BEFORE ANY MECHANISM

**1. The largest effect this lane measured is in the READOUT, not in anything I was asked to tune.**
CVaR is the *training* objective and the readout is an *argmin*, and nothing forces the two to use
the same Hamiltonian. Separating them experimentally, on the identical evaluated set, at zero extra
budget:

    change the READOUT Hamiltonian   -0.697 A [-1.059, -0.352]  15W/5L  5/5 folds
    change the TRAINING Hamiltonian  +0.056 A [-0.020, +0.129]  10W/10L 4/5 folds

**Which Hamiltonian selects is worth twelve times more than which Hamiltonian trains**, it holds on
all ten arms including `best_of_N` where there is no training at all, and it replicates
independently in Blocks A and Q at different shot counts. **The selector carries everything; the
sampler carries nothing.**

**2. The ansatz question's quantum half is CLOSED by register size, not by bond dimension** — which
my own pre-registration named first, as the thing that would damage my framing most. Its RMSD half
is a clean null: the whole ladder from chi = 1 to a measured chi = 48 spans **0.37 A on 20 targets
with an MDE of 0.2-0.4 A**, and gradient norms *rise* with depth, so it is not a trainability null.

**3. The optimizer question separates into two columns that disagree, and that IS the answer.**
Adam beats every alternative on the CVaR training loss at **0W/20L (DIST) and 0W/14L (LEG), 5/5
folds, every CI excluding zero** — and **no optimiser, Adam included, beats `best_of_N` on
Ca-RMSD**. Training the circuit at all is worth **+0.003 A [-0.136, +0.144]**. QNG does not help,
and the reason is a property of the **parameterisation** (the MPS Fisher is rank-deficient by
construction, condition number 2.5e17 at 32 shots) rather than of the landscape — the distinction
the brief said not to conflate.

**4. Two things I got wrong and corrected in place**: I quoted a 17x conditional from a 15-target
snapshot that did not survive to n=20 (§3e, retracted), and my own F-A3 threshold demanded a
precision an n=20 design does not have (§"precision", declared, regime 0.084-0.30 A left OPEN).

---

## 1. VALIDITY GATES ON THIS LANE'S OWN NEW MACHINERY — run before any number was read

BRIEF §10 says read the code before believing the claim. This lane wrote three new pieces of
arithmetic — a per-sample score matrix, a Fisher preconditioner, and an exact bond-dimension
measurement — so each was gated first. `s21/results/qb3_machinery_gate.json`, `_COMPLETE`.

| gate | `mps_L2` n=9 | `mps_L2` n=12 | `sv_ring_L3` n=9 |
|---|---|---|---|
| `max ‖E_p[∇log p]‖` (must be 0 — an identity) | 3.0e-16 | 6.7e-16 | 5.9e-16 |
| `∇log p` vs central FD, max rel. err | 1.1e-9 | 4.9e-10 | 1.0e-10 |
| classical Fisher rank / P | **26 / 27** | **36 / 36** | **27 / 27** |
| Fisher λ_min | 2.4e-8 | 9.9e-5 | 7.9e-3 |
| `cos(∇, natural ∇)` | +0.447 | +0.515 | +0.423 |
| measured χ (amplitude) vs analytic `2^layers` | 4 vs 4 | 4 vs 4 | 16 vs n/a (no analytic χ) |

Two things worth carrying:

* **The score identity holds to machine precision and the FD check passes**, so the CVaR gradient,
  the Fisher and the natural-gradient direction are all arithmetically what they claim to be. The
  new `SVAnsatz` analytic parameter-shift path (`dp/dθ_j = ½[p(θ+π/2 e_j) − p(θ−π/2 e_j)]`, exact
  for an RY generator because `p(b)` is the expectation of a projector) agrees with finite
  differences to 1e-10.
* **The MPS parameterisation carries a genuine redundancy**: its Fisher is rank-deficient by one at
  n=9 and its smallest eigenvalue is 2.4e-8. The ring statevector circuit is full rank with
  λ_min = 7.9e-3. **Natural gradient on the MPS is therefore dominated by the damping in the null
  direction, and that is a property of the parameterisation, not of the landscape.** It is the
  first thing to hold against any Block-Q result.

`cos(∇, natural ∇) ≈ +0.42–0.52` — the preconditioner genuinely rotates the update by ~60°, so a
null in Block Q is not a null because the preconditioner did nothing.

### The harness gate for Block E (F-E5) — it fired, it was not vacuous

`s21/results/qb3_gate.json`, `_COMPLETE`. `best_of_N` draws from the same von Mises basin mixtures
in all 8 parameterisations and is mapped back through the same retraction, so its emitted torsion
set must be identical. Measured at full n: **max |Δ| = 1.11e-16 over 35,840 configuration
comparisons across all 20 targets and all 8 charts**. Not bit-identity — one ulp — and the source is `arctan2(r sin θ, r cos θ)` recovering
`θ` to within a rounding, not a harness mismatch. **The eight charts are physically identical to
machine precision, measured rather than asserted.**

*(Recorded because I got the flag wrong first: `qb3_gate_COMPLETE` was initially written from a
4-target smoke run, which violates the rule that a completion flag must require the FULL
configuration. The flag was removed, `run_gate` was changed to write `_COMPLETE` only when
`len(rows) == N_T` and a `__PARTIAL_` otherwise, and the gate was re-run at n = 20.)*

**And a stronger gate than F-E5 asked for** (`s21/results/qb3_gauge_exact.json`, `_COMPLETE`). F-E5
only checks that `best_of_N` emits the same set. The premise of F-E1 is stronger: *every* chart must
return the *identical objective* at the *identical physical point*. Measured over 8 charts x 24
arbitrary configurations x 3 targets:

    max |objective difference| across charts at the same physical point   7.77e-16   (504 comparisons)
    max |round-trip angle error|, radii 0.25 / 0.5 / 1 / 2 / 4 / 7.3      0.00e+00

**The `arctan2` round trip is EXACT, not approximate, at every radius tested.** So the eight
parameterisations are physically identical to machine precision and **any difference Block E
measures is optimiser geometry, by construction and not by assumption.**

---

## 2. BLOCK A — THE ANSATZ. The quantum-resource half is CLOSED; the RMSD half is a clean null.

### 2a. F-A1 fired, and the coordinator's audit closed it independently

My pre-registration states F-A1 first *because it damages my own lane's framing most*: every target
here has n ≤ 16 residues at **one qubit per residue**, so the full statevector is ≤ 65,536
amplitudes and **every ansatz at this problem size is classically samplable by brute force,
independent of bond dimension**. The coordinator's audit lane enumerated the complete distribution
for every register size the 126 targets use at layers 1–20 — 17.7 s total, `sum p = 1.0000000000`,
**64 ms dense at χ saturated to 256** — and the dense route costs `2^n` independent of depth.

> **There is no bond dimension at which classical simulability breaks for this circuit.**
> "Escapes bond dimension 4" is *achievable* (see 2b) and *not sufficient*: the binding constraint
> is the register, not the bond. Reaching a hard register needs n ≳ 50, i.e. 50-residue chains,
> outside this benchmark entirely — and there Cerezo et al. (arXiv:2312.09121, *Nat. Commun.*,
> Aug 2025) applies: the structure that removes barren plateaus is the structure that admits
> classical simulation.
>
> **CLOSED, not untested.** No ansatz redesign changes it while the register is one qubit per
> residue on 9–16-mers. This does **not** touch the finding that the VQE genuinely trains, and it
> does **not** say a 39-parameter restricted model is a bad optimiser.

### 2b. Bond dimension, MEASURED rather than quoted — and the source claim is exact

`MPSAnsatz.__init__` sets `self.chi = D` with `D` doubling once per entangling layer, i.e.
`χ = 2^layers`. I did not take that on trust: the amplitude tensor is contracted exactly
(`amp(x) = e₀ᵀ ∏_q A[q][:, x_q, :] e₀`) and its Schmidt rank read at every bipartition.

| ansatz | P | χ analytic | **χ measured (amplitude)** | χ measured (probability) | H₀ bits | H₁ bits | ‖∇‖ |
|---|---|---|---|---|---|---|---|
| `prod` (no CNOTs) | 39 | 1 | **1** | 1 | 7.30 | 5.14 | 0.499 |
| `mps_L1` | 26 | 2 | **2** | 3 | 8.44 | 6.70 | 0.376 |
| `mps_L2` **deployed** | 39 | 4 | **4** | 10 | 9.63 | 7.75 | 0.519 |
| `mps_L3` | 51 | 8 | **8** | 30 | 10.16 | 9.34 | 0.676 |
| `mps_L4` | 64 | 16 | **16** | 63 | 10.62 | 9.61 | 0.704 |
| `mps_L2_nofinal` | 26 | 4 | **4** | 5 | 8.43 | 7.13 | 0.429 |
| `share_L2` (3 params) | 3 | 4 | **4** | 10 | 8.39 | 6.89 | 0.619 |
| `sv_ring_L3` (ring, exact SV) | 35 | — | **48** | 66 | 9.02 | 7.69 | 0.566 |

**The analytic χ and the measured Schmidt rank agree exactly on all seven MPS-family variants** — EXACT, and
a check on the source rather than a finding. `sv_ring_L3`'s ring closure is long-range for an open
MPS and it reaches a mean realised Schmidt rank of **48**, i.e. it *does* escape the chain's bond
dimension by a factor of 12. Per 2a that buys no resource claim; it is here as an expressivity rung.

**Trainability does not degrade over the reachable ladder.** ‖∇‖ *rises* monotonically with depth
(0.376 → 0.704) and latent entropy after training rises with depth (6.70 → 9.61 bits). **There is no
barren plateau to escape at n ≤ 16 and 8 gradient steps.** So expressivity and trainability do not
trade here — the pre-registered F-A2 mechanism does not fire, and the null in 2c is not a
trainability null.

### 2c. F-A3 FIRES. The whole expressivity ladder is worth less than the seed noise.

`s21/results/qb3_a.json`. **n = 20 targets, 4 seeds, B = 512, shots = 64, α = 0.25**, objective
`DIST` (the deployed distogram axis). Paired target-level bootstrap, seeds averaged within target.
`sv_ring_L3` is capped at n ≤ 14 (exact statevector) and its comparisons are on the matched
14-target subset, which is stated in the row rather than pooled silently.

| ansatz | χ measured | **− `mps_L2`** | **− `best_of_N`** | absolute |
|---|---|---|---|---|
| `prod` | 1 | +0.011 [−0.195,+0.220] | +0.057 [−0.186,+0.292] | 4.121 |
| `mps_L1` | 2 | −0.102 [−0.257,+0.053] | −0.056 [−0.211,+0.094] | 4.008 |
| **`mps_L2`** deployed | 4 | — | +0.047 [−0.102,+0.197] | 4.110 |
| `mps_L3` | 8 | +0.002 [−0.206,+0.194] | +0.049 [−0.220,+0.295] | 4.112 |
| `mps_L4` | 16 | −0.068 [−0.243,+0.136] | −0.021 [−0.250,+0.215] | 4.042 |
| `mps_L2_nofinal` | 4 | +0.171 [−0.045,+0.401] | +0.217 [−0.036,+0.480] | 4.281 |
| `share_L2` (3 params) | 4 | −0.029 [−0.227,+0.169] | +0.018 [−0.170,+0.208] | 4.081 |
| **`sv_ring_L3`** | **48** | −0.122 [−0.329,+0.092] | −0.081 [−0.261,+0.093] | 3.736 |
| `untrained_L2` | 4 | −0.012 [−0.197,+0.170] | +0.035 [−0.166,+0.232] | 4.099 |
| `best_of_N` | — | −0.047 [−0.197,+0.102] | — | **4.064** |

> **Not one CI excludes zero. A 48× increase in realised bond dimension moves Cα-RMSD by
> −0.122 Å [−0.329, +0.092], and the entire ladder from χ=1 to χ=48 spans 0.37 Å — less than twice
> the 0.200 Å ansatz-seed sd this instrument is known to carry.** — **F-A3 FIRES.** Ansatz
> expressivity is closed on this encoding, from a third direction, after Sprint 20's entanglement
> null and its `x_deep` result. Label: **NOT MEASURED** for every individual arm; **ESTABLISHED**
> for the joint statement that no rung of the ladder reaches the MDE.

**And the null is not a trainability null**, which is the reading that would have let the branch
survive: ‖∇‖ *rises* with depth and the deepest ansatz has the largest gradient. **F-A2's mechanism
does not fire.** Expressivity and trainability do not trade here; both are available, and neither
buys structure.

**The single most useful row is `best_of_N` at 4.064** — the lowest mean of any arm except the
14-target `sv_ring` subset. Zero optimisation, matched budget. *Search saturates, discrimination
binds*, on a seventh instrument.

---

## 3. BLOCK Q — THE OPTIMIZER. The objective separates decisively; RMSD does not move at all.

`s21/results/qb3_q.json`. **n = 20 targets, 4 seeds, B = 512, shots = 32 (16 gradient steps),
α = 0.25**, ansatz fixed at the deployed `mps_L2`. `shots = 32` rather than Sprint 20's 64 is a
declared deviation, made so that an *optimiser* comparison gets 16 steps rather than 8; it is
matched across every arm in the block.

### 3a. On `DIST` — the two columns disagree, and that is the result

**CVaR training loss at the last iteration** (the genuine training objective; *not* `best_std`,
which is the lowest single evaluated energy and is reported separately):

| arm | CVaR − adam | W/L | folds |
|---|---|---|---|
| `sgd` | **+0.237 [+0.121, +0.328]** | 1W/19L | 4/5 |
| `qng` | **+0.235 [+0.110, +0.343]** | 1W/19L | 4/5 |
| `qng_diag` | **+0.214 [+0.146, +0.284]** | 2W/18L | 5/5 |
| `spsa_ansatz` | **+0.420 [+0.301, +0.560]** | 0W/20L | 5/5 |
| `untrained` | **+0.524 [+0.341, +0.764]** | 0W/20L | 5/5 |

**Every CI excludes zero.** Adam is decisively the best optimiser of the CVaR loss, and the circuit
decisively trains: trained-vs-untrained is +0.524 at **0W/20L, 5/5 folds**.

**Cα-RMSD, the same arms, the same runs, the same targets**, against `best_of_N`:

| arm | RMSD − best_of_N | median | W/L | absolute mean / median |
|---|---|---|---|---|
| `adam` | +0.080 [−0.165, +0.355] | −0.069 | 11W/9L | 4.144 / **3.832** |
| `sgd` | +0.047 [−0.190, +0.305] | −0.013 | 10W/10L | 4.111 / 3.966 |
| `qng` | +0.097 [−0.123, +0.335] | +0.052 | 9W/11L | 4.160 / 4.066 |
| `qng_diag` | +0.180 [−0.052, +0.426] | +0.099 | 8W/12L | 4.243 / 4.143 |
| `spsa_ansatz` | +0.052 [−0.111, +0.224] | +0.005 | 10W/10L | 4.116 / 3.990 |
| `untrained` | +0.088 [−0.100, +0.280] | +0.037 | 8W/12L | 4.152 / 4.003 |
| **`best_of_N`** | — | — | — | **4.064** / 4.103 |

Every CI spans zero. **`best_of_N` has the lowest mean of the seven arms** and every optimiser sits
on the wrong side of it.

**The median-vs-mean warning fires on `adam`, and is reported rather than exploited.** Adam's mean is
+0.080 *worse* than `best_of_N` while its median is −0.069 *better*, at 11W/9L — a near-even W/L
with mean and median on opposite sides. That is the signature of an arm that wins narrowly and loses
catastrophically, i.e. concentration, and per `median-vs-mean-is-the-free-warning` it is an early
warning, **not** a result: the CI spans zero, so the label is **NOT MEASURED**, and I am not
entitled to read the median as a win.

> **This is the cleanest instance of `search-saturates-discrimination-binds` the programme has: one
> table where the objective column is significant at 0W/20L with 5/5 folds and the structure column
> is flat, on identical runs.** — **F-Q1 FIRES**, exactly as pre-registered.

**F-Q2 also fires on this Hamiltonian**: no optimiser beats `best_of_N`, so on `DIST` the limitation
is discrimination, not the optimiser, and **no optimiser claim may be made here**.

### 3b. QNG specifically — it is not a null of inaction, and the reason it fails is measurable

The preconditioner rotates the update by ~60° (`cos(∇, natural ∇) = +0.42–0.52`, §1), so `qng` is
doing something. What it is doing is dominated by damping: **the empirical Fisher's median condition
number is 2.5e17**, because 32 shots cannot span 39 parameters, and §1 shows the MPS
parameterisation is *intrinsically* rank-deficient (Fisher rank 26/27 at n=9) on top of that.

> **Geometry-aware preconditioning on a sampler whose Fisher is rank-deficient by construction is
> damping-limited, and that is a property of the PARAMETERISATION, not of the landscape.** State it
> that way rather than as "QNG does not help": the untested regime is `shots ≥ P` with a full-rank
> parameterisation, and it is recorded **OPEN**.

### 3c. Two more instruments for the decoupling, on a design where difficulty CANNOT operate

Both are computed **within each target, across arms**, so target difficulty — the confound L4
showed swallows every raw correlation on this instrument — is constant by construction rather than
partialled out.

    rho(CVaR TRAINING LOSS, Ca-RMSD)  across the 6 optimiser arms   -0.011 [-0.180, +0.163]   9/20 positive
    rho(best evaluated ENERGY, Ca-RMSD) across the 10 ansatz arms   +0.030 [-0.165, +0.239]  10/20 positive

**Coin flips.** Seventh and eighth independent instruments for
`search-saturates-discrimination-binds`, and the first two this programme has that remove the
difficulty confound by design instead of by regression.

### 3d. Latent concentration is a TARGET-DIFFICULTY proxy — L4 replicated on new hardware

`s21/results/qb3_concentration.json`, `_COMPLETE`. L4's standing requirement is that a raw rank
correlation against Ca-RMSD on this instrument is largely a difficulty measurement. Control is
`pool_mean_ORACLE` (an ORACLE difficulty proxy, labelled).

| ansatz | raw rho(H1, RMSD) | **partial** | raw rho(max_prob, RMSD) | **partial** |
|---|---|---|---|---|
| `prod` | +0.380 | +0.055 | -0.338 | -0.027 |
| `mps_L1` | +0.577 | +0.118 | -0.522 | -0.054 |
| `mps_L2` | +0.496 | +0.106 | -0.475 | -0.159 |
| `mps_L3` | +0.346 | -0.023 | -0.403 | -0.071 |
| `mps_L4` | +0.344 | -0.027 | -0.301 | +0.123 |
| `mps_L2_nofinal` | +0.486 | +0.146 | -0.429 | -0.066 |
| `share_L2` | +0.253 | -0.036 | -0.198 | +0.049 |
| `sv_ring_L3` | +0.297 | +0.002 | -0.354 | -0.121 |
| `untrained_L2` | +0.460 | +0.214 | -0.556 | -0.299 |

Raw +0.25 to +0.58 collapses to -0.04 to +0.21. **And within a target, across the eight
expressivity rungs, rho(entropy, RMSD) = -0.106 [-0.312, +0.123]** — a more concentrated latent
does not build a better structure once the confound cannot operate.
**`concentration-is-wrong-when-discrimination-binds`, reproduced with the mechanism attached and
the difficulty confound removed by design.**

### 3e. A CONDITIONAL I OVER-READ AT n=15 AND MUST RETRACT AT n=20 — recorded in full

On Legacy at **15 targets** this lane measured `untrained - mps_L2 = -0.242`, `mps_L2 - best_of_N =
+0.224 with 5/5 folds`, and a cost-of-training of **+0.204 A on Legacy against +0.012 A on the
distogram — a 17x ratio**, and I reported it to the coordinator in those terms. The mechanism was
attractive and pre-existing: Sprint 20 measured Legacy's pool argmin at **5.487 against a pool mean
of 4.739**, so a sampler that gets better at finding low-Legacy configurations should get worse at
building structures, while the distogram's argmin (3.676) sits *below* the pool mean.

**At the full n = 20 it does not hold up.**

| | pool rho(E, RMSD) | pool argmin | untrained - best_of_N | trained - best_of_N | **cost of training** |
|---|---|---|---|---|---|
| **Legacy** | +0.191 | 5.487 | -0.017 [-0.332, +0.303] | +0.122 [-0.124, +0.373] | **+0.138 [-0.071, +0.372], 9W/11L, 3/5** |
| **distogram** | +0.497 | 3.676 | +0.035 [-0.166, +0.232] | +0.047 [-0.102, +0.197] | **+0.012 [-0.170, +0.197], 9W/11L, 3/5** |

The Legacy cost of training fell from +0.204 with **5/5 folds** to +0.138 with **3/5 folds and a CI
spanning zero** on five additional targets. **NOT MEASURED**, and I retract the "17x" framing.

> **What survives is only the ORDERING of two point estimates in the predicted direction, at n=20,
> with both CIs spanning zero.** That is a PLAUSIBLE mechanism and not a measured effect. The
> correct label for "training the CVaR-VQE costs more on a badly-ranking Hamiltonian than on a
> well-ranking one" on this instrument is **OPEN**, and testing it needs the n=126 instrument, not
> more arms on 20 targets.
>
> **Recorded as a self-caught error** in the same class as Sprint 20's `x_warm` min-of-N artefact:
> an interim reading, quoted with a fold count that the remaining targets removed.

### 3f. Block A on Legacy, final, n = 20

Every arm sits within +-0.17 A of `best_of_N` and **no CI excludes zero**; pooled, the mean of the
eight ansatz rungs is `+0.062 [-0.136, +0.258]` against `best_of_N` (MDE 0.295) and
`+0.079 [-0.134, +0.302]` against the untrained circuit. The distogram result (2c) reproduces on the
second Hamiltonian: **the ansatz ladder does nothing on either.**

---

## 4. THE READOUT HAMILTONIAN — the BRIEF's own open question (b), answered, and it is the largest effect this lane measured

The BRIEF was corrected mid-sprint: **CVaR is the TRAINING objective; the readout is an argmin.**
It then named two live regimes, one of which is *"a readout H different from the training H — a
genuinely new and untested design"*. **That question costs nothing to answer**: every arm already
evaluated a set of 512 configurations, so re-taking the argmin under a different Hamiltonian
consumes no budget, builds no new structure, and for `DIST`/`LEG` makes no OpenMM call
(measured: 0.24 s for both readouts over 512 configurations).

**Identity gate first.** For a `K`-trained run, `xread_K` must reproduce the arm's own readout
exactly. Measured **156/156 exact** to all sixteen digits. The diagonal is an identity, not a
finding, and it is verified rather than assumed.

### Pooled over all ten arms, n = 20, 4 seeds — the two knobs, separated experimentally

| what changes | difference | CI | W/L | folds | MDE |
|---|---|---|---|---|---|
| **the READOUT H** (train LEG: read DIST − read LEG) | **−0.697** | **[−1.059, −0.352]** | 15W/5L | **5/5** | 0.518 |
| the TRAINING H (read DIST: train LEG − train DIST) | +0.056 | [−0.020, +0.129] | 10W/10L | 4/5 | 0.111 |

    absolute:  train LEG / read LEG   4.849
               train LEG / read DIST  4.151
               train DIST / read DIST 4.096

> **The readout Hamiltonian is worth −0.697 Å with a CI excluding zero and 5/5 folds. The training
> Hamiltonian is worth +0.056 Å, with the design excluding anything beyond ±0.13 Å.**
>
> **You may train the genuine CVaR-VQE on Legacy and lose 0.056 Å — provided you take the argmin
> with the distogram. Take it with Legacy and you lose 0.70 Å.**

Per-arm the readout effect is −0.522 to −0.840, **every one of the ten CIs excluding zero**,
including `best_of_N`, whose set involves **no training at all** (−0.679 [−1.175, −0.226]). So it
is a property of the **selector**, not of the sampler — `legacy-corrects-the-selector` and Sprint
20's L13 (*"the carrier is the distogram-as-selector, not the pool"*), now measured **inside** the
CVaR-VQE with the two Hamiltonians separated rather than inferred across pipelines.

**What this does and does not license.** It does **not** beat the incumbent: the deployed pipeline
already selects with the distogram, so the finding *confirms* that design rather than improving on
it. What it adds is a **price**: the Hamiltonian in the training loop is nearly free to change, and
the Hamiltonian in the readout is the expensive choice. **Any future arm that puts a physics
Hamiltonian into a CVaR-VQE should be asked which of the two roles it is being put into**, because
the roles differ by more than an order of magnitude in cost. — **SUPPORTED**, n = 20, 4 seeds,
native-free throughout (the distogram readout uses no native information).

### 4a. THE VALIDITY GUARD — and it does NOT come out cleanly in my favour

BRIEF §3: `caonly_k300` was 0.110 A more accurate with an indistinguishable clash count and
`cis_frac` +0.4242 — **a scalar validity score would have promoted a broken structure. Validity is a
VECTOR.** A -0.7 A readout swap is exactly the kind of claim that precedent exists to catch, so the
two readouts' argmins were compared on geometry as well as accuracy.
`s21/qb3_validity.py` -> `s21/results/qb3_validity.json`, `_COMPLETE`, n = 20, 4 seeds.
**Reproducibility gate: 240/240 re-run rows bit-exact against `qb3_a.json`.**

**What is measurable and what is vacuous.** Every configuration here is built by
`build_ca_exact` from continuous (phi, psi) with **ideal** bond geometry, so bond lengths, bond
angles and omega are ideal *by construction* — reporting them as "passing" would be a gate that
never fires. **And the Ramachandran axis cannot referee this particular comparison at all: it is a
TERM INSIDE the Legacy energy** (`core.energy.rama_penalty`), so Legacy-selected structures would
score better on it by construction. *A validity metric that is a term in one of the two competing
objectives is not a neutral referee*, and it is left out rather than reported as a foregone win.
What is free to vary and neutral: non-adjacent Ca contacts, closest approach, and compactness.

| arm | readout | Ca-RMSD | clash <4.0 A | clash <4.5 A | min non-adj Ca | rg − native rg |
|---|---|---|---|---|---|---|
| `mps_L2` | LEG | 4.923 | 0.200 | 0.588 | 4.489 | **−0.528** |
| `mps_L2` | **DIST** | **4.083** | 0.512 | 0.963 | 4.482 | **−0.475** |
| `untrained_L2` | LEG | 4.785 | 0.200 | 0.525 | 4.541 | −0.508 |
| `untrained_L2` | **DIST** | **4.114** | 0.562 | **1.100** | 4.338 | −0.450 |
| `best_of_N` | LEG | 4.801 | 0.237 | 0.637 | 4.451 | −0.546 |
| `best_of_N` | **DIST** | **4.123** | 0.512 | 1.038 | 4.302 | −0.450 |

> **The distogram readout is NOT free.** It buys 0.67-0.84 A of accuracy and costs a **higher count
> of close non-adjacent Ca contacts on all three arms** — +0.28 to +0.36 below 4.0 A and +0.38 to
> +0.58 below 4.5 A — and on `untrained_L2` that increase is **significant**:
> **+0.575 [+0.013, +1.212]**. The other two arms' CIs span zero. Closest approach moves the same
> way (−0.006 to −0.203, all ns).
>
> **In absolute terms it is small** — about one extra Ca pair under 4.5 A per structure, out of
> 60-100 non-adjacent pairs — but it is a real, measured cost and it is stated beside the benefit
> rather than after it.

**And one axis moves the RIGHT way, with a mechanism.** Sprint 20 established that
**Legacy is a compactness model**: Legacy-preferred candidates are 0.45 A more compact. Here
**Legacy's argmin is over-compact against the native by −0.51 to −0.55 A, and the distogram readout
recovers about a tenth of that** (−0.45 to −0.48). The readout swap is partly *undoing Legacy's
compactness bias* — which is what Sprint 20's characterisation predicts it should do, and is a
mechanism rather than a coincidence.

**Scope, stated.** This is a **Ca-level** validity vector on **LEG-trained** sets. All-atom clashes,
`cis_frac` and omega — the axes on which `caonly_k300` actually failed — are **NOT MEASURED here**,
and the readout claim should not be promoted into the pipeline until they are.

### 4c. FOLD-CLUSTERED CIs ON EVERY DISPOSITIVE CONTRAST — 3 of 24 dispositions change

BRIEF §9: `s12.instrument.paired`'s CI is **i.i.d. over targets, not fold-clustered**, and
`s18.phys_lib.paired`'s `ci_fold` must be used where the disposition could depend on it. Every
dispositive contrast in this lane was re-run through it. `s21/results/qb3_ci_fold.json`,
`_COMPLETE`, 24 contrasts.

| contrast | i.i.d. CI | fold-clustered CI | disposition |
|---|---|---|---|
| READOUT swap, Block A (−0.697) | [−1.059, −0.352] | [−1.119, −0.275] | unchanged |
| READOUT swap, Block Q (−0.708) | [−1.047, −0.366] | [−1.149, −0.266] | unchanged |
| `untrained − adam` [CVaR] DIST / LEG | excludes 0 | excludes 0 | unchanged |
| **TRAINING swap, Block A (+0.056)** | [−0.020, +0.129] | **[+0.002, +0.095]** | **gains significance** |
| **`qng − best_of_N` [RMSD] DIST (+0.097)** | [−0.123, +0.335] | **[+0.006, +0.185]** | **gains significance** |
| **Block S LEG, 32−4 steps [CVaR] (−0.096)** | [−0.158, −0.034] | [−0.171, +0.0004] | **loses significance** |

**Three of twenty-four.** Two considerations follow, and the second is the important one.

* The **readout** result — this lane's headline — is **unchanged and slightly wider** under fold
  clustering. It does not depend on the CI construction.
* **The disposition change that matters strengthens a NEGATIVE result, not a positive one.** Under
  the construction the BRIEF prefers, **`qng` is significantly WORSE than zero optimisation on the
  distogram: +0.097 Å [+0.006, +0.185], above the 0.084 Å MDE.** That is a firmer verdict on
  quantum natural gradient than my i.i.d. interval gave, and it is the version I report.
* The Block-A **training swap** gains significance at +0.056 Å — **still below the MDE**, so the
  conclusion ("the training Hamiltonian is worth an order of magnitude less than the readout")
  stands, but the direction is consistent within folds and I flag it rather than quote only the
  interval that spans zero.

### 4b. A CONFOUND IN MY OWN BLOCK Q, and the thing it makes visible

`adam` moves **5.6x further in parameter space than `sgd`** at matched learning rate
(total ||dtheta|| 4.687 vs 0.830; `qng` 1.191, `qng_diag` 1.009, `spsa_ansatz` 2.221,
`untrained` 0.000). So **"adam is the better optimiser of the CVaR loss" is substantially "adam
takes the biggest steps"**, and my optimiser battery is matched on learning rate rather than on
step norm. That is the *same* confound the encoding question turns on, sitting inside my own block,
and it is declared rather than left for someone else to find. A matched-step-norm battery is
recorded **OPEN**.

**And the confound is what makes the decoupling legible.** Within each target, across the six arms:

    rho(parameter displacement, CVaR TRAINING LOSS)  = -0.411 [-0.497, -0.320]   CI excludes zero
    rho(parameter displacement, Ca-RMSD)             = -0.006 [-0.217, +0.211]

> **The one quantity that demonstrably drives the training objective has exactly no relationship to
> the structure.** Moving further in parameter space reliably lowers the loss and does nothing
> whatever to Ca-RMSD.
>
> This also bears directly on the audit lane's "fewer steps helps" mechanism, from a completely
> independent direction: **on the variational arm, moving LESS lowers the objective less and does
> not improve the structure.** Step count is not acting on RMSD here in either direction.

---

## 5. BLOCK Q ON LEGACY, n = 20 — every DIST conclusion replicates, and more sharply

**CVaR training loss vs `adam`** — the objective column, adequately powered:

| arm | CVaR − adam | W/L | folds |
|---|---|---|---|
| `sgd` | **+0.277 [+0.205, +0.358]** | 0W/20L | 5/5 |
| `qng` | **+0.288 [+0.224, +0.352]** | 0W/20L | 5/5 |
| `qng_diag` | **+0.291 [+0.216, +0.372]** | 0W/20L | 5/5 |
| `spsa_ansatz` | **+0.354 [+0.294, +0.417]** | 0W/20L | 5/5 |
| `untrained` | **+0.341 [+0.248, +0.442]** | 1W/19L | 5/5 |

**Adam beats every alternative on 20 targets out of 20, on both Hamiltonians, with every CI
excluding zero and 5/5 folds. QNG is on the wrong side of Adam on both.**

**Ca-RMSD vs `best_of_N`** — the structure column: `adam` +0.162, `sgd` +0.090, `qng` +0.081,
`qng_diag` +0.141, `spsa_ansatz` +0.089, `untrained` −0.007. **Every CI spans zero and every
trained arm is on the wrong side of doing nothing.** Pooled, the five optimisers are
+0.113 [−0.159, +0.394] against `best_of_N` and **+0.120 [−0.041, +0.283] against the untrained
circuit** — on Legacy, training is on the *harmful* side, though not measurably so.

> **F-Q2 fires on the second Hamiltonian too.** No optimiser beats `best_of_N` on Legacy either, so
> **the limitation is the landscape/discrimination and not the optimiser, on both objectives
> tested** — and no optimiser claim may be made from this lane. **F-Q3's contrast (does QNG help on
> one Hamiltonian and not the other?) returns NO on both**: `qng − adam` is +0.235 on `DIST` and
> +0.288 on `LEG`, both on the wrong side, both CIs excluding zero, on the objective; and both null
> on RMSD. The AMBER arm, which is where bad conditioning was supposed to be QNG's case, is queued
> behind the box's OpenMM serialisation and its status is in §9.

### 5b. The readout result REPLICATES INDEPENDENTLY

Block A (10 arms, shots = 64) and Block Q (7 arms, shots = 32) are separate runs with separate
seeds, separate ansatz sets and separate optimisers, on the same 20 targets:

    BLOCK A pooled   read DIST - read LEG   -0.697 [-1.059, -0.352]  15W/5L  5/5 folds
    BLOCK Q pooled   read DIST - read LEG   -0.708 [-1.047, -0.366]  17W/3L  5/5 folds

**They agree to 0.011 A.** The readout effect is not a property of one block's configuration.

---

## 6. BLOCK S — the shots/iterations trade. "Fewer steps helps" is NOT SUPPORTED on the variational arm.

PREREG addendum §6, falsifier fixed before the data. At fixed `B = 512`,
`shots x iterations = B`, so the shot count is an **8x sweep of optimisation effort with the
encoding, chart, dimension, probe size, ansatz, starts, seeds and budget all held exactly fixed**.
Nothing else in this programme isolates step count that cleanly. `s21/results/qb3_s.json`.

### DIST, n = 20, 4 seeds

| shots | gradient steps | CVaR loss | ‖dtheta‖ | RMSD mean | RMSD median | vs `best_of_N` |
|---|---|---|---|---|---|---|
| 16 | **32** | **−1.048** | 5.726 | 4.145 | 4.225 | +0.082 |
| 32 | 16 | −0.998 | 4.771 | **4.090** | 4.011 | **+0.027** |
| 64 | 8 | −0.958 | 3.637 | 4.169 | 4.159 | +0.105 |
| 128 | **4** | **−0.798** | 2.420 | **4.223** | 4.086 | **+0.159** |
| — | 0 | — | 0 | **4.064** | 4.103 | — |

**The objective falls monotonically with step count and the design has clear power for it**:
32 steps beats 4 steps by **−0.250 [−0.438, −0.122], 18W/2L, 5/5 folds**, and every intermediate
step count is significant in the same direction.

**The structure column does not follow it, and what movement there is runs the WRONG WAY for the
hypothesis.** Against the fewest-steps arm, every higher-step arm has a *negative* RMSD point
estimate (−0.077, −0.132, −0.053), and **the 4-step arm is the worst RMSD arm on the board**
(+0.159 against `best_of_N`, 5/5 folds).

### Legacy, n = 18 — the same picture on the second Hamiltonian

| shots | steps | CVaR loss | RMSD mean | vs `best_of_N` |
|---|---|---|---|---|
| 16 | 32 | **0.378** | **4.590** | −0.063 |
| 32 | 16 | 0.338 | 4.656 | +0.003 |
| 64 | 8 | 0.372 | 4.655 | +0.001 |
| 128 | **4** | **0.483** | **4.666** | +0.013 |
| — | 0 | — | 4.653 | — |

CVaR again separates significantly (−0.104 to −0.145 vs the 4-step arm, CIs excluding zero, up to
17W/1L) and RMSD again does not (all four arms within 0.08 A of `best_of_N`, every CI spanning
zero, every one below MDE). **And again all three higher-step arms have negative RMSD point
estimates against the fewest-steps arm.**

> **F-S1 does NOT fire** — RMSD does not improve as the step count falls; it is the fewest-steps
> arm that is worst.
> **F-S3 fires** — the 4-vs-32-step RMSD difference is −0.077 [−0.285, +0.102], below this design's
> MDE, so the sweep is **NOT MEASURED** on RMSD and F-S2 may not be claimed from it either. The
> honest statement is that **an 8x change in optimisation effort moves Ca-RMSD by 0.13 A against an
> MDE of ~0.25 A, with the point estimates opposite in sign to "fewer steps helps".**

**Taken with §4b** — `rho(parameter displacement, Ca-RMSD) = −0.006 [−0.217, +0.211]` while
`rho(parameter displacement, CVaR loss) = −0.411 [−0.497, −0.320]` — **this lane has two
independent instruments saying step count is not an RMSD lever on the variational arm, in either
direction.** That does not settle the encoding question, which lives on the classical torsion arms
and belongs to the audit lane's control; it removes one mechanism from the list of things that
could explain it *here*. — **NOT SUPPORTED**, n = 20, 4 seeds.




---

## 7. WHAT DID NOT FINISH, AND WHAT IS THEREFORE OPEN

Recorded here rather than left to be inferred from a missing table. Every incomplete artefact
carries a `qb3_<tag>__PARTIAL_` file naming the exact shortfall against the FULL declared
configuration, not against the subset the call happened to run.

| block | objective | status |
|---|---|---|
| A (ansatz) | `DIST`, `LEG` | **COMPLETE**, n=20, 4 seeds |
| A (ansatz) | `AMBc` | **NOT RUN** — lowest priority, queued behind everything else |
| Q (optimizer) | `DIST`, `LEG` | **COMPLETE**, n=20, 4 seeds |
| Q (optimizer) | `AMBc` | **running at hand-off — 56 of 560 rows (2 targets) at the time of writing.** Far too few to report; the driver continues and `qb3_q__PARTIAL_` tracks it |
| S (shots/steps) | `DIST`, `LEG` | **COMPLETE**, n=20, 4 seeds |
| S (shots/steps) | `AMBc` | **NOT RUN** |
| E (gauge sweep) | all | **NOT RUN** — de-prioritised on the coordinator's redirect |

**The consequence for F-Q3, stated plainly.** F-Q3 asks whether QNG helps on one Hamiltonian and
not another. **On the two Hamiltonians that completed the answer is NO on both** — `qng − adam` is
+0.235 (`DIST`) and +0.288 (`LEG`) on the objective, both CIs excluding zero and both on the wrong
side, and null on RMSD. **The badly-conditioned AMBER landscape, which is the one place QNG's case
was strongest a priori, is NOT MEASURED by this lane**, and I am not entitled to generalise the
Legacy/distogram result to it. That regime is **OPEN**.

**Block E is OPEN and its design is on record.** The gauge sweep is built, gated (§1: the eight
charts agree to 7.77e-16 and the `arctan2` round-trip is exact) and runnable with one command
(`python -m s21.qb3_run e:DIST 20`). Its distinctive contribution, which the audit lane's
matched-iteration control cannot supply, is that **the angle-scale family `s ∈ {0.5, 1, 2}` varies
effective step SIZE at exactly constant step COUNT and constant probe cost** — every angle chart is
`2n`-dimensional — so size and count are separable in it. `emb_norm` additionally retracts ‖u‖ every
step, which is a direct test of the radius-diffusion annealing the audit lane measured (1.00 → 1.94,
max 3.06).

---

## 8. COMPUTE, AND WHY NO TIMING NUMBER IN THIS REPORT IS A PROPERTY OF THE CODE

BRIEF §10: *measure timings on a quiet box or not at all.* **The box was not quiet.** Through this
lane's run there were consistently **4-6 other heavy Python processes** on it (`s21.c_norm`,
`s21.tailprice` x2, `s21.d_enc`, `s21.a_matrix`, `s21.c_cont`), and physical memory sat at
**83-94%** against `core.amber.memory_guard`'s 92% ceiling. Two consequences, both recorded:

* **Per-arm wall times in the artefacts are contended and are NOT comparable across blocks.** One
  smoke-test cell (`AMBc / mps_L4`, n=16) took **1000 s** where its sibling took 7 s; re-running it
  was impossible because the guard then refused to open an OpenMM context at all. That is memory
  pressure, not a pathological ansatz, and it is exactly the failure mode Sprint 20 recorded when a
  "pathological minimiser" turned out to be CPU starvation.
* **The AMBER-bound blocks were queued behind everything that does not need OpenMM**, and the
  runner was given a retry that waits out a guard refusal instead of losing the block
  (`with_mem_retry`; refusals are counted per target in `mem_guard_refusals`).
* **And I set that wait wrong the first time, which is worth recording because it is the same class
  of mistake as a vacuous gate.** `mem_wait` was initially `limit = 88%, patience = 3600 s` — 3.5
  points *below* the 92% at which OpenMM actually refuses. On a box hovering at **89%** that spends
  up to an hour of wall clock per target waiting for a refusal that was never going to happen. **A
  guard tuned more conservatively than the thing it guards against is not free; it costs more than
  the failure it prevents.** Corrected to `limit = 90.5%, patience = 600 s`, with `with_mem_retry`
  catching a real refusal, and the AMBER step restarted (row-level resume, nothing recomputed).

**A process-race I caused, found, and closed.** Killing the first driver's Python without killing
its shell left the shell's `for` loop alive, so two drivers wrote the same checkpoint for about six
minutes. **It cost rows, and it could not corrupt them**: every row is keyed by
`stable_rng(pdb, kind, arm, seed)`, so both processes computed identical values and a race can only
lose a row, never change one. **Verified rather than argued**: 20 rows were recomputed from the
persisted seeds and reproduced **bit-exactly, 20/20**. Recorded because the next lane to use
`nohup sh` on this box will hit it.

---

## 9. WHAT THIS LANE CLOSED, AND WHAT IT DID NOT

| claim | label | basis |
|---|---|---|
| Bond dimension is not the binding constraint; the register is | **CLOSED** | coordinator's enumeration; my F-A1 pre-registered the same argument |
| `MPSAnsatz` realised chi == analytic `2^layers` | **EXACT** | Schmidt rank of the contracted amplitude tensor, 6/6 variants |
| A ring circuit escapes the chain's bond dimension (chi 4 -> 48) | **EXACT** | measured, n<=14 subset |
| Expressivity and trainability do not trade at n<=16 | **SUPPORTED** | ‖∇‖ rises 0.376 -> 0.704 with depth |
| No ansatz rung beats `mps_L2` or `best_of_N` on Ca-RMSD | **NOT MEASURED** at <=0.30 A; effects >0.3 A excluded | n=20 MDE 0.20-0.38 |
| Ansatz expressivity "closed" as F-A3 specified it | **MIS-SPECIFIED prereg; regime 0.084-0.30 A is OPEN** | §"precision" above |
| Adam optimises the CVaR loss best; every alternative loses | **ESTABLISHED** (DIST **and** LEG) | 0W/20L to 2W/18L on both, all CIs exclude zero, 5/5 folds |
| ...and it is confounded with step norm | **declared**; matched-step battery **OPEN** | ‖dtheta‖ 4.687 vs 0.830 |
| QNG does not beat Adam even on the objective | **SUPPORTED** (DIST **and** LEG) | +0.235 [+0.110,+0.343] and +0.288 [+0.224,+0.352] |
| QNG is significantly WORSE than zero optimisation on Ca-RMSD | **SUPPORTED** under the fold-clustered CI | +0.097 [+0.006, +0.185] (DIST), above MDE |
| QNG is damping-limited because the MPS Fisher is rank-deficient | **EXACT** (rank 26/27, 36/36) + **SUPPORTED** (cond 2.5e17) | machinery gate |
| No optimiser beats `best_of_N` on Ca-RMSD | **SUPPORTED** (DIST, LEG) | every CI spans zero, `best_of_N` lowest mean |
| Training the circuit is worth +0.003 A [-0.136, +0.144] | **SUPPORTED** | pooled 5 optimisers vs untrained, DIST |
| **The READOUT Hamiltonian is worth ~-0.70 A, 5/5 folds** | **SUPPORTED**, replicated in two blocks | A: -0.697 [-1.059,-0.352]; Q: -0.708 [-1.047,-0.366]; 17/17 arm-level CIs exclude zero |
| **The TRAINING Hamiltonian is worth +0.056 A [-0.020, +0.129]** | **SUPPORTED** | same table, MDE 0.111 |
| ...and the readout swap costs close Ca contacts | **SUPPORTED on 1 of 3 arms** (+0.575 [+0.013,+1.212]); all-atom / cis / omega **NOT MEASURED** | §4a |
| ...and it partly undoes Legacy's compactness bias | **SUPPORTED**, mechanism matches s20's L8 | rg vs native -0.53 -> -0.48 |
| Training costs more on a badly-ranking H than a well-ranking one | **OPEN** — I over-read it at n=15 and **retract** | §3e |
| Latent concentration predicts structure | **REFUTED** after the difficulty partial | §3d |
| rho(objective gained, structure gained) is zero | **ESTABLISHED**, two new within-target instruments | §3c |
| Parameter displacement drives the objective and not the structure | **SUPPORTED** | rho -0.411 [-0.497,-0.320] vs -0.006 [-0.217,+0.211] |
| "Fewer gradient steps helps" on the variational arm | **NOT SUPPORTED** — point estimates run the other way | Block S, 8x step sweep, DIST and LEG |
| The 8 charts are physically identical (F-E1's premise) | **EXACT** | 7.77e-16 objective agreement; arctan2 round-trip 0.00e+00 |
| The lane's own new machinery (scores, Fisher, chi) is correct | **EXACT / verified** | E_p[score]=0 to 1e-15; grad vs FD to 1e-10 |
| F-Q3: does QNG help on the badly-conditioned AMBER landscape? | **NOT MEASURED** — AMBER arm did not finish | §7 |

### The three questions I was given, answered

1. **The encoding.** Redirected by the coordinator mid-run; the deciding matched-step control is
   theirs. My contribution is the **gauge sweep** — `r` in the circle chart and `s` in the angle
   chart are exact gauge parameters, and the angle family varies effective step SIZE at **constant
   step COUNT and constant probe cost**, which the matched-iteration control cannot do. Status
   below.
2. **The ansatz.** The quantum-resource half is **CLOSED by register size**, which my own
   pre-registered F-A1 named as the thing that would damage my framing most. The RMSD half is a
   **clean null on two Hamiltonians** at a resolution of ~0.3 A, with the 0.084-0.30 A regime
   explicitly **OPEN**.
3. **The optimizer.** On the distogram and Legacy the limitation is **not the optimiser**: no
   optimiser beats `best_of_N`, so no optimiser claim may be made (F-Q2). QNG specifically fails
   for a reason that is a property of the **parameterisation** (a rank-deficient Fisher), not of
   the landscape — and that distinction is the one the brief said not to conflate.

**And the RMSD connection, which I was told to hold myself to.** Every arm in this lane that lowered
an objective is reported as a negative result on structure, because that is what it is. The
objective column of Block Q is significant at 0W/20L with 5/5 folds. The structure column is flat.


---

## 10. WHAT I WOULD DO NEXT, IN PRIORITY ORDER

1. **Price the readout choice on the n=126 instrument, with the full validity vector.** §4 is this
   lane's only large effect (-0.70 A, replicated in two blocks, 5/5 folds) and it is measured on 20
   targets with a **Ca-level** validity check that already shows a cost. Before anything is built on
   it, it needs all-atom clashes, `cis_frac` and omega — the axes on which `caonly_k300` actually
   failed — and n=126. **It will not improve the incumbent** (which already selects with the
   distogram); its value is that it prices the two ROLES a Hamiltonian can play, and they differ by
   more than an order of magnitude.

2. **Stop asking whether the ansatz or the optimiser can be improved.** Three blocks, two
   Hamiltonians, 4 seeds, ~4,500 budgeted arms: the objective column separates at 0W/20L with 5/5
   folds and the structure column never moves. The remaining untested regimes are named in §7 and
   §9 and none of them is large: matched-step-norm optimisers, `shots >= P` for a full-rank Fisher,
   and the 0.084-0.30 A band this design cannot resolve.

3. **If the AMBER arm matters, run it alone on a quiet box.** F-Q3 is the one assigned question
   this lane could not answer, and the reason is entirely compute: 4-6 competing processes and a
   92% memory ceiling. It is a ~40-minute job unattended and a multi-hour job contended.

4. **Block E is built and gated; run it if the audit lane's control leaves the encoding open.** Its
   distinctive contribution is separating step SIZE from step COUNT — the angle-scale family
   `s in {0.5, 1, 2}` changes effective step size at exactly constant iteration count and constant
   probe cost, which a matched-iteration control cannot do. One command: `python -m s21.qb3_run
   e:DIST 20`.

5. **Carry the two method lessons rather than the numbers.** (a) **MDE is per-comparison**: this
   lane's own contrasts range 0.10-0.45 A and the 0.084 A constant would have licensed three
   conclusions the data do not support. (b) **A guard tuned more conservatively than the thing it
   guards is not free** — an 88% memory wait against a 92% ceiling cost more than the refusal it
   prevented.
