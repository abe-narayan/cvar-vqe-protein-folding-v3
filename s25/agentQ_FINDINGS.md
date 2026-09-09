# LANE Q — CVaR-VQE / QUANTUM. SPRINT 25.

Pre-registered in `s25/PREREG_Q.md`; the Q-B fork list went to the coordinator and was
approved before the run. Q-A, Q-C and the plateau measurement are **property measurements** —
they read no RMSD and have no outcome that can favour a hypothesis — so Rule 0's fork
enumeration does not apply (the s24 Lane-D precedent for `d_setequality_proof.py`) and
registered falsifiers stand in its place. ORACLE / ACHIEVABLE / PRODUCTION labelled at every
appearance. No benchmark inspection anywhere in this lane. `LOCK_TRAIN` and `LOCK_AMBER` were
not taken.

**BASIS NOTICE, BINDING ON EVERY NUMBER BELOW.** Two instruments appear in this lane and they
never share a column.

| instrument | operator | numbers |
|---|---|---|
| **s8** (`s8/integrate_vqe.json`) | consensus-medoid **selection** of ONE member from a 128-candidate score-filtered set | 3.4540 / 3.3414 / 3.2835 / 3.3135 |
| **126-target dev** | score-filtered uniform **top-75 coordinate average** | 3.0483 (point cloud), **3.2148 (built chain)** |

Per the coordinator's L4, the 3.0483 Å arm is a point cloud whose virtual Cα–Cα bonds are
**22.3% short, worst case 0.649 Å — shorter than a covalent C–C bond** — and `core/bench.py:682`
labels it *"raw average (illegal)"*. **This lane does not quote 3.0483 Å as a structural result
anywhere.** Where a system number is needed, the built chain **3.2148 Å** is used and named as
such; the synthesis arm's projection gap is **+0.1664 Å** (not +0.156, which is `rmsd_fit`, a
different arm).

---

## 0. THE RESULT IN SIX LINES

1. **Everything I was handed about the ALGORITHM verified.** Exact MPS at χ = 2^layers
   (3.3e-16 against a dense simulator I wrote from the gate list), exact parameter shift
   (1.000000000 cosine, 4.6e-10 relative error), the set-equality theorem re-derived on my
   own families with exact zeros (0 violations in 2,592 cells), the `baseline="const"` fix
   in place on both estimators.
2. **The one thing I was handed about the CONTRIBUTION did not.** "+0.113 Å" is arithmetically
   exact and reproduces to four decimals, but as a paired contrast it is **0.51× its own MDE
   with a fold CI spanning zero — NULL by the project's own rule.** It was labelled off
   marginal means.
3. **And it is measured at T = 0.1, which is not deployed.** `core/pipeline.py:118` runs
   T = 0.3 on all five folds and **α = 1.0 on three of them — 78 of 126 targets, 61.9%, carry
   no tail constraint at all.** At T = 0.3 the tail is worth **−0.0011 Å, 0.01× MDE**.
4. **One curve explains both knobs.** RMSD tracks the **entropy of the readout weights**
   (ρ = −0.7423 over the 9 VQE cells) and not α (+0.2700) or T (−0.0234); α adds **3.2%** of
   the variance left after H and H². Fit the curve on the **nine no-circuit arms** and the
   nine VQE arms sit on it at **+0.0090 Å against a residual sd of 0.0268 Å**.
5. **The genuinely positive quantum result: the state reaches its own analytic optimum.** At
   α = 1 the objective is `mean(E) − T·H(p)`, whose unconstrained simplex minimiser is exactly
   `exp(−E/T)/Z`. The 21-parameter RY/CNOT state matches that exact optimum to within noise at
   all three temperatures (0.24–0.42× MDE). **That is also the mechanism for why it buys
   nothing: what it approximates is itself classical and cheap.**
   **>>> THIS LINE IS SUPERSEDED BY §2.5b. I flagged the wording myself, measured it on the
   coordinator's ruling, and the falsifier fired.** The optimiser genuinely trains (beats
   best-of-200 from the untrained circuit at all three T, closes 78–89% of the free-energy
   gap), but it lands **0.902 nats and 45% of its mass** from that optimum at the deployed
   temperature — **and the endpoint cannot tell (0.24× MDE). The readout is insensitive to a
   distributional difference of nearly half the mass**, which is a sharper statement than the
   one it replaces.
6. **Five defects of one class are now on the record** (§5) — a default, a fallback or a
   marginal that reads as a PASS. Four are mine to report in someone else's code, one is the
   coordinator's, and the pattern is the most transferable thing this lane produced.

---

## 1. Q-A — THE FOUR STANDING FACTS, INDEPENDENTLY RE-VERIFIED

`s25/q_verify.py` → `s25/results/q_verify.json`, provenance-stamped.

**What "independent" had to mean here.** s24 Lane D's REV1 asserted a *false* theorem and
reported 0/2916 failures because every probability family it generated carried an epsilon
floor, so the assertion could not fire. A green light from a dead assertion is worse than a
red one. So nothing below is checked against a sibling routine in the same module.

### 1.1 `MPSAnsatz` — χ = 2^layers by construction, exact, no truncation. **CONFIRMED.**

    truncation primitives inside the MPSAnsatz class          NONE
    every occurrence of "svd"/"truncat" in core/quantum.py    PROSE ONLY (lines 25, 465, 1591)
    "cutoff" in core/quantum.py                               CVAR_SORT_CUTOFF, the sort/select
                                                              crossover -- unrelated to bonds
    chi by layers (cnot chain, final_ry)      {1: 2, 2: 4, 3: 8, 4: 16}   == 2**layers
    chi with entangler="none"                 1               (product state, the control)
    ring closure applied by MPSAnsatz         False
    built tensor shape at n=8, layers=2       (8, 4, 2, 4)    = (n, chi, 2, chi)
    chi at n=64, layers=2                     4               -- independent of n

Amplitudes were checked against a **dense simulator written in `q_verify.py` from the gate
list**, by explicit Kronecker products — not against `StatevectorCircuit` (a sibling in the
same module, and a different topology: it applies a ring closure, `MPSAnsatz` does not).

    max |p_MPS - p_dense| over 24 configs (n in 3..9, L in 1..3, final_ry both)   3.331e-16
    max | <psi|psi>_MPS - 1 |                                                     6.661e-16

**The MPS claim is verified in the strongest available sense** — not "no truncation call
exists" but "the state equals the exact state". And the bond dimension is not a bound that
happens to hold, it is *saturated*: the Schmidt spectrum across the middle cut at n=6,
layers=2 is `[0.879, 0.462, 0.106, 0.054, 0, 0, 0, 0]` — **rank exactly 4 = χ**, four
non-zero values and four exact zeros. The entangler is doing something: the same angles
through `entangler="none"` move a probability by up to 0.334.

**A scope point that matters for the specification.** `MPSAnsatz` is **not** used by
`core/pipeline.py`. The deployed selector is `StatevectorCircuit` at n=7, layers=3. `MPSAnsatz`
at layers=2 (χ=4) is the *generation-lane* ansatz — s19/s20/s21. The inherited phrasing
"deployed at layers=2, χ=4" is correct about the MPS lane and would be wrong about the
selector, and `s25/QUANTUM.md` separates the two.

### 1.2 `run_cvar_vqe` — Adam on the exact parameter-shift gradient. **CONFIRMED.**

    deployed register                     n=7 qubits, layers=3, dim=128, P=21 parameters
    max |p_statevector - p_dense|         5.551e-17     (against the same independent simulator)
    cos(parameter shift, exact FD)        1.000000000   min over 12 (alpha in {0.1,0.25,1.0})
    rel |g_ps - g_fd| / |g_fd|            4.597e-10     max over 12
    cos(free-energy grad, FD)             1.000000000   min over 4 (alpha,T) incl. the deployed
    rel err of the free-energy grad       4.663e-10     max over 4  -- BOTH terms, CVaR and -TH
    sampling calls in run_cvar_vqe        NONE          rng used for initial angles only

The free-energy row is the one that matters and it had not been separately recorded: the
deployed objective is `CVaR_α − T·H(p)`, and it is the *gradient of the whole thing*, entropy
term included, that had to be exact. It is.

### 1.3 The recorded gradient defect — **reproduces as a defect, NOT at the recorded constant**

    cvar_gradient default baseline        "const"     (the fix)
    grad_cvar_score default baseline      "const"     (the fix)
    cos(exact-expectation, CONST)         1.000000
    cos(exact-expectation, TAIL)          0.566586    <- the defect, zero sampling noise
    |g_tail| / |g_exact|                  0.519

**The defect is real and my instrument reproduces its character exactly** — a large bias at
zero sampling noise, which is what proves it is bias and not variance. **It does not reproduce
the constant.** `core/quantum.py` records +0.655634 at 0.758× (s9, 10 qubits, 1024 amplitudes);
project memory records +0.524 at 0.534× on a third instrument; I get +0.5666 at 0.519× at
n=7, layers=3, α=0.15 on the deployed energy shape. **The number is instrument-dependent and
must always be cited with its instrument.** It is quoted that way in `s25/QUANTUM.md`. This is
not a discrepancy to resolve — three instruments agreeing on the sign, the order of magnitude
and the mechanism *is* the reproduction; a shared constant would have been a coincidence.

### 1.4 The set-equality theorem — **re-derived independently, 0 violations, assertion LIVE**

Own families, built from scratch, **with exact zeros throughout** (`p = 0.0`, never `+1e-30`):
6 register sizes × 9 energy structures × 8 probability structures × 6 alphas, including an
adversarial family that zeroes out *exactly the lowest-energy quarter* of the states.

    n_cells                                              2592
    cells containing an EXACT zero probability           1620      <- the assertion CAN fire
    SUBSET-hood violations   <- THE THEOREM                 0
    holes that were NOT exactly zero-probability            0      <- the mechanism, directly
    prefix-hood violations (the s24 REV1 over-claim)      1424     = 54.9%
    full-support cells                                    972
    full-support cells with EXACT value equality      972 / 972    = 100%

**And I checked the assertion is alive rather than trusting that it is.** A constructed case
with `p[argmin] = 0` at α = 0.30: the realised tail excludes state 0, and is still a subset of
the prefix. Both reported, both true. `tail_indices(energies, alpha)` takes no probability
vector at all — the other half of the theorem, confirmed from the signature.

My prefix-violation rate is 54.9% against Lane D's 29.9% and the coordinator's 58.8%. **The
rates differ because the family mixes differ and that is expected; the three assertions
agree exactly.** Subset-hood is the theorem; equality is the empirical full-support regime.

---

## 2. Q-B — WHAT THE CVaR TAIL ACTUALLY CONTRIBUTES. **THE HEADLINE.**

`s25/q_alpha.py` → `s25/results/q_alpha.json`, n=126, provenance-stamped. Re-analysis of a
fixed on-disk artefact; no new pipeline run. **s8 instrument throughout.**

**DECLARED WEAKENING:** I read the artefact's stored per-arm marginal means before filing the
forks. No paired SE, MDE, CI or verdict had been computed. Declared in the prereg and here.

### 2.1 The inherited claim, re-verified — and then labelled properly

    vqe_a1.0_T0.1   3.4540   entropy 0.0761 bits    identical to the argmin arm PER TARGET
    vqe_a0.1_T0.1   3.3414   entropy 6.3638 bits
    difference     +0.1126   <- the inherited "+0.113 A". Arithmetically exact.

The collapse is real and stronger than reported: at α=1, T=0.1 the arm is not merely *equal in
mean* to the shipped argmin, it is **identical target by target**. The 0.076-bit entropy is
not an approximation to a collapse, it is a collapse.

**But the contrast's own label:**

    a=0.1 - a=1.0 @T=0.1   -0.1126  SE 0.0792  MDE 0.2220  0.51x  55W/44L  med +0.0000  NULL

> **"+0.113 Å of genuine contribution" was never a measured effect, even in the slice where it
> was measured.** 0.51× its own MDE, fold CI spanning zero, a median of exactly 0.0000, and
> W/L barely off even. It was labelled on marginal means where a paired statistic was
> required. **This is the same class as s24 Lane D's `paired_stats` defect, arriving from the
> other direction** — not a harness flattering its user, but a marginal quoted where a paired
> label belonged. The retraction is recorded as the coordinator's, in `s25/LEDGER.md` L5.

### 2.2 PRIMARY — the α effect does not survive the temperature. **MY FALSIFIER DID NOT FIRE.**

    T = 0.1     a=0.10 - a=1.0   -0.1126  SE 0.0792  0.51x  55W/44L  fold 3/5  NULL
                a=0.25 - a=1.0   -0.1067  SE 0.0644  0.59x  36W/25L  fold 4/5  UNDERPOWERED
    T = 0.3     a=0.10 - a=1.0   +0.0279  SE 0.0518  0.19x  50W/50L  fold 4/5  NULL
                a=0.25 - a=1.0   -0.0011  SE 0.0477  0.01x  44W/47L  fold 3/5  NULL
    T = 1.0     a=0.10 - a=1.0   +0.0126  SE 0.0268  0.17x  28W/19L  fold 4/5  NULL
                a=0.25 - a=1.0   +0.0039  SE 0.0240  0.06x  18W/24L  fold 3/5  NULL

**The effect does not reproduce at the other two temperatures and it changes sign.** The
registered falsifier — "if it reproduces, my prediction is wrong and the inherited attribution
stands" — did not fire.

### 2.3 THE DEPLOYED QUESTION, ANSWERED PLAINLY

`core/pipeline.py:118`, verified from source and reconstructed per target:

    VQE_LFO = {0: (1.0, 0.3), 1: (0.25, 0.3), 2: (0.25, 0.3), 3: (1.0, 0.3), 4: (1.0, 0.3)}
    reconstructed from core/pipeline   3.3135   == the stored vqe_LFO arm, per target: True
    targets on an alpha = 1.0 fold      78 / 126 = 61.9%   -- NO TAIL CONSTRAINT AT ALL

> **At the deployed temperature a genuine tail constraint is worth nothing measurable:
> −0.0011 Å at 0.01× MDE (α=0.25), +0.0279 Å at 0.19× MDE (α=0.10). Both NULL.**

**ORACLE counterfactual, labelled ORACLE, and NOT a proposal.** Forcing α<1 on every fold is
worth −0.0311 Å (α=0.25, 0.27× MDE, UNDERPOWERED) or −0.0021 Å (α=0.10). **There is nothing
there worth changing a leave-fold-out table for.** That is the useful half of this finding:
the α=1 folds are not costing accuracy. **The α=1 problem is not a performance problem. It is
a claims problem, and the fix is to change what is said, not what is run.**

### 2.4 ONE CURVE. THE READOUT'S ENTROPY IS THE VARIABLE; α AND T ARE TWO HANDLES ON IT.

`s8/integrate.py:1272-1320` hands **every** arm the same operator, `consensus_medoid(D, o, w)`
— the member of the same 128-candidate set minimising a `w`-weighted mean distance. The arms
differ **only** in `w`. And `E = zrank(s[o])` where `o` is already the score order, so E is the
standardised ranks 1..128, **the same vector for every target**. So the readout entropy is
computable *exactly* for every arm with no re-running: `argmin` → 0 bits, `topfrac_f` →
log₂(k), `boltz_T` → H(exp(−E/T)), `medoid128` → 7 bits, VQE → the stored per-arm mean.

    arm               H_readout (bits)   mean RMSD    family
    argmin                     0.0000      3.4540     degenerate
    vqe_a1.0_T0.1              0.0761      3.4540     VQE (circuit)
    topfrac_0.05               2.5850      3.4106     uniform top-fraction, no circuit
    boltz_T0.1                 3.3326      3.4041     Boltzmann, exact, no circuit
    topfrac_0.1                3.7004      3.3557     uniform top-fraction, no circuit
    boltz_T0.3                 4.9135      3.3138     Boltzmann, exact, no circuit
    topfrac_0.25               5.0000      3.3680     uniform top-fraction, no circuit
    vqe_a1.0_T0.3              5.6605      3.2835     VQE (circuit)
    vqe_a0.25_T0.1             5.9223      3.3473     VQE (circuit)
    topfrac_0.5                6.0000      3.2873     uniform top-fraction, no circuit
    vqe_a0.25_T0.3             6.2993      3.2825     VQE (circuit)
    vqe_a0.1_T0.1              6.3638      3.3414     VQE (circuit)
    boltz_T1.0                 6.4429      3.3008     Boltzmann, exact, no circuit
    vqe_a0.1_T0.3              6.5234      3.3115     VQE (circuit)
    vqe_a1.0_T1.0              6.8071      3.3453     VQE (circuit)
    vqe_a0.1_T1.0              6.9261      3.3579     VQE (circuit)
    vqe_a0.25_T1.0             6.9627      3.3493     VQE (circuit)
    medoid128                  7.0000      3.3443     degenerate

    corr(H, mean RMSD) over the 9 VQE cells       -0.7423
    corr(alpha, mean RMSD)                        +0.2700
    corr(T, mean RMSD)                            -0.0234
    alpha's marginal share of the variance left after H and H^2      0.032
    quadratic RMSD ~ H + H^2 over all 18 arms     R^2 = 0.7045

**The construction that makes it airtight: fit the curve on the NINE NO-CIRCUIT arms only —
arms the circuit is not in — then score the nine VQE arms against it.**

    mean residual, VQE arms          +0.0090 A   (sd 0.0342, n=9)
    residual sd of the fit itself     0.0268 A

> **The circuit sits on a curve fitted without it, one third of a residual standard deviation
> above.** α and T are two handles on one quantity — how concentrated the readout weights are
> — and once that quantity is conditioned on, neither the knobs nor the circuit is visible.

### 2.5 THE EXPRESSIVITY RESULT — the genuinely positive quantum finding

At α = 1 the objective is `mean_p(E) − T·H(p)`, and its **unconstrained minimiser over the
probability simplex is exactly the Boltzmann distribution `exp(−E/T)/Z`**. So the `boltz_T`
arm is not a loose analogy: **it is the exact optimum the 21-parameter RY/CNOT state is
approximating.** The comparison is therefore a clean expressivity measurement.

    circuit - exact Boltzmann @T=0.1   +0.0499  SE 0.0427  0.42x  24W/37L  UNDERPOWERED
    circuit - exact Boltzmann @T=0.3   -0.0302  SE 0.0450  0.24x  43W/44L  NULL
    circuit - exact Boltzmann @T=1.0   +0.0445  SE 0.0511  0.31x  38W/47L  NULL

**THIS IS WHERE I WAS WRONG, AND I CAUGHT IT BY FLAGGING MY OWN SENTENCE.** My draft said the
state "attains its own analytic optimum". That is a claim about DISTRIBUTIONS; what the table
above measures is that two arms' RMSDs are indistinguishable — a claim about a downstream
ENDPOINT, and an underpowered one. **Claiming a mechanism from an endpoint is the exact error
that produced this project's L2.** I flagged it to the coordinator rather than shipping it; the
ruling was *measure the divergence and earn the sentence, and if it is not small that is the
more interesting result*. **It is not small. §2.5b replaces it.**

### 2.5b The divergence, measured. **THE FALSIFIER FIRED AGAINST MY OWN DRAFT.**

`s25/q_gibbs.py` → `s25/results/q_gibbs.json`. Registered falsifier: KL > 0.1 nat makes the
wording unearned. Exact identity `F(p) − F(p*) = T·KL(p‖p*)`, **asserted at runtime to < 1e-9
on every row** so it cross-checks both computations rather than assuming either.

    T = 0.30  (DEPLOYED on folds 0, 3, 4)
      state                                    F        KL nats   KL bits      TV   H bits
      random theta (untrained)             -1.275443    3.927427   5.66608  0.78052  4.4562
      uniform over 2**n                    -1.455609    3.326874   4.79966  0.70154  7.0000
      point mass at argmin (the collapse)  -1.718572    2.450332   3.53508  0.91374  0.0000
      5 Adam steps                         -1.644194    2.698258   3.89276  0.70994  5.1019
      15 Adam steps                        -1.908150    1.818405   2.62340  0.63160  5.4014
      50 Adam steps (DEPLOYED)             -2.183096    0.901916   1.30119  0.45308  5.6706
      the Gibbs optimum itself             -2.453671    0.000000   0.00000  0.00000  4.9135

    KL(trained || optimum) nats:  T=0.1  1.449   T=0.3  0.902   T=1.0  0.373
    total variation:              T=0.1  0.761   T=0.3  0.453   T=1.0  0.351

**The comparator ladder is the point** — a divergence has no units without one, and the
coordinator asked for exactly this. At the deployed temperature the trained state is **0.902
nats from its own optimum and disagrees with it on 45% of its mass.** That is unearned by any
reading of "attains".

**The mandatory control, which is the genuine positive** (`concentration-is-wrong-when-
discrimination-binds`: best-of-N from the UNTRAINED circuit, never an initialisation mean):

      T      F_trained   F_init_mean   F_init_best    F_gibbs   gap closed   beats best-of-200
    0.10     -1.717571     -0.526843     -0.988733  -1.862495        0.891         True
    0.30     -2.183096     -1.207063     -1.483528  -2.453671        0.783         True
    1.00     -4.936715     -3.587831     -4.214344  -5.309822        0.783         True

> **Three readings, in this order.** (i) **The optimiser genuinely trains** — it beats
> best-of-200 from the untrained circuit at every temperature and closes 78–89% of the
> available free-energy gap. (ii) **It does not reach the optimum** — 0.902 nats and 45% of the
> mass away at the deployed temperature, and *broader* than optimal (5.67 bits against 4.91),
> not collapsed. (iii) **The endpoint cannot tell** — 0.24× MDE.
>
> **The readout is insensitive to a distributional difference of nearly half the mass.** That is
> the mechanism behind §2.4 and it is a sharper statement than the one I drafted: it is *why*
> an exact classical Gibbs state, a trained quantum state 45% away from it, and a flat average
> over the top 64 candidates all land within a few hundredths of an Ångström of each other.

**Two details worth keeping.** At T = 0.1 the trained state's free energy (−1.7176) is
marginally *worse* than the plain point mass at the argmin (−1.7186) and its entropy is 0.075
bits — at that temperature the entropy term is too weak to hold the state open and the circuit
converges to what is effectively the collapse. At T = 1.0 the trained state (KL 0.373) is only
modestly better than the **uniform** distribution (KL 0.458).

### 2.5c A source-level consequence, verified on real targets

`top = argsort(sc)` (`core/pipeline.py:768`), so `sc[top[:128]]` is already ascending and
`_zrank` returns the standardised ranks 1..128. **E is therefore very nearly the same vector
for every target** — but *not exactly*, because `rankdata` averages ties and the score does tie.
Measured on 8 real targets: worst deviation **4.06e-02 against an E range of 3.4371 = 1.18% of
the range**, and **0 of 8 exactly identical**. My first draft said "numerically the same vector
for every target"; that was too strong and the check caught it.

> **The spectrum of H is target-independent to within tie-averaging** — not exactly constant,
> but within about 1% of its own range. **ALL of the per-target information enters through which
> candidate occupies which rank, and essentially none through the spectrum of H.**

**Three consequences, and the third is a mechanism the project has been missing.**

1. The deployed selector solves very nearly the **same variational problem** on all 126
   targets. There are effectively **two trained states in the whole deployment** — one per
   `(alpha, T)` cell of `VQE_LFO` — not 126.
2. The quantum stage is therefore **insensitive to the target** by construction: the object it
   optimises barely knows which protein it is looking at.
3. **It independently explains a result measured five times across five sprints and never
   given a mechanism: deeper ansätze and larger χ order nothing.** A more expressive state can
   only pay if `H` carries target-specific structure for it to capture. It carries almost none
   — it is a fixed ladder of standardised ranks. **Extra expressivity has nothing to be
   expressive about.** That is a mechanism rather than a restatement, and it was available from
   four lines of source the whole time.

### 2.6 Against the no-circuit arms

    VQE_LFO - argmin (shipped)    -0.1405  SE 0.0732  MDE 0.2051  0.68x  66W/48L  5/5 folds
                                                                          UNDERPOWERED
    VQE_LFO - Boltzmann T=0.3     -0.0002  SE 0.0457  0.00x  39W/43L  NULL
    VQE_LFO - uniform top-64      +0.0262  SE 0.0279  0.34x  26W/37L  NULL
    VQE_LFO - uniform top-128     -0.0308  SE 0.0590  0.19x  49W/49L  NULL

The headline contrast is **0.68× its own MDE — below its own MDE, not a result by our rule** —
though its direction is consistent across all 5 folds and 66W/48L. **Against an exact
classical Boltzmann weighting at the same temperature the difference is 0.0002 Å.**

### 2.7 Concentration check, run because the medians are zero

`median-vs-mean-is-the-free-warning` fires on `VQE_LFO − argmin`: median −0.0081 against a mean
of −0.1405, 12 of 126 targets tie exactly, 5 targets carry 45.3% of the effect. **But the
drop-top-5 statistic compared against a uniform-effect null sits at the 64.7th percentile of
that null** (observed −0.0403, null mean −0.0668, 95% band [−0.2060, +0.0755]).

> **Concentration is SUGGESTED by the median-vs-mean gap and NOT established by the drop-top
> test.** Reported as both, per the memory entry, which records that a raw drop-top threshold
> is not a valid test and misfired here once before.

---

## 3. GRADIENT VARIANCE — THE BARREN-PLATEAU MEASUREMENT

`s25/q_plateau.py` → `s25/results/q_plateau.json`. Property measurement, no RMSD.
**Results are transcribed into `s25/QUANTUM.md` §7 with their scope conditions.**

The design's own control costs nothing and is exact: **at α = 1 the CVaR reduces identically to
`<ψ|H|ψ>` with `H = diag(E)` — the linear cost of the textbook barren-plateau setting** — and
α < 1 is the same circuit, the same H, the same θ law, with only the non-linearity switched
on. So "does the CVaR non-linearity change the variance picture" is answerable by a ratio at
matched n, rather than by argument.

**Scope, stated so it cannot be over-read.** The deployed register is n = 7; every larger n is
an extrapolation instrument and is labelled so. P = 3n at layers = 3, which is exponentially
smaller than dim so(2ⁿ) = 2ⁿ⁻¹(2ⁿ − 1) — 21 parameters against 8128 at n = 7 — so **this
circuit is nowhere near a 2-design at any n measured** and an observed decay rate is a property
of *this shallow ansatz*, never a 2-design result. RY and CNOT are real, so the reachable
manifold lies in SO(2ⁿ), not SU(2ⁿ).

### 3.1 No exponential plateau at any width measured

    fitted log2 Var[dF/dth_0] per qubit, n = 4..13, layers = 3
      alpha = 1.0, T = 0   (THE LINEAR COST, <psi|H|psi>)      -0.6492
      alpha = 0.25, T = 0                                      -0.2522
      alpha = 0.10, T = 0                                      -0.0472    essentially FLAT
      alpha = 1.0,  T = 0.3   (deployed form)                  -0.3105
      alpha = 0.25, T = 0.3   (deployed form)                  -0.2429

A 2-design would give about −1 per qubit. The steepest column here is −0.65, and it is the
**linear** one. Draws are 250 (n ≤ 8), 200, 120, 80; relative SE of a variance is
`sqrt(2/(m−1))`, so 9% to 16%.

### 3.2 The non-linearity flattens the decay — the answer to the coordinator's question 2

    ratio Var[grad CVaR_alpha] / Var[grad MEAN], identical circuit, identical theta law
      n        alpha=0.25   alpha=0.10
      4          0.2712       0.0877
      7          0.5159       0.0932     <- deployed width
      8          1.4549       0.5241
     10          2.0336       1.1704
     13          3.3585       2.6765

**Two opposite answers depending on the question.** In *magnitude* at the deployed width the
non-linearity shrinks the gradient (0.52× at α=0.25, 0.093× at α=0.10) — unsurprising, since
the tail weight `(q−E)_+/α` is supported on an α fraction of states. In *scaling* it flattens
the decay, the ratio crosses 1 near n ≈ 8, and by n = 13 the CVaR gradient's variance is
**2.7–3.4× the linear one's**. On this ansatz the CVaR objective is better conditioned at
width than the linear cost it reduces to. A mechanism is offered in `QUANTUM.md` §7.2 **as a
reading, explicitly not as a result** — the project has not tested it.

### 3.3 Depth saturates; it does not attenuate

    n = 7, alpha = 0.25, T = 0.3, 250 draws
      L    P    Var[dF/dth_0]
      1    7    2.381e-02
      2   14    1.658e-02
      3   21    9.197e-03    <- DEPLOYED
      4   28    7.999e-03
      6   42    7.644e-03
      8   56    7.440e-03
     12   84    8.213e-03

Depth costs a factor of ~3 and then flattens: L = 4 through L = 12 are equal within sampling
error. **There is no exponential attenuation with depth in the measured range**, and the
deployed L = 3 already sits close to the plateau. (This row's n = 7 value differs from the
width sweep's by ~2 sampling SEs; different seeds, consistent.)

---

## 4. Q-C — THE HARNESS AUDIT. THE FIX IS IN; IT HAS FOUR SIBLINGS.

`s24/d_harness.py` is 18/18 and reusable as delivered. **I measured the MDE fix rather than
reading it:**

    an effect placed at 0.39x its own MDE  ->  UNDERPOWERED (<0.7x MDE, NOT A RESULT)
    at 1.00x                              ->  TYPE-M ZONE (0.7-1.3x MDE, NOT A RESULT)
    at 2.00x                              ->  MEASURED
    MDE == 2.8016 * SE exactly            ->  True

The fix is in and behaves correctly at all three regimes. **The siblings are in §5.**

---

## 5. FIVE DEFECTS OF ONE CLASS. **THE PATTERN IS THE OUTPUT.**

> **A DEFAULT, A FALLBACK, OR A MARGINAL THAT READS AS A PASS.** In every case a quantity that
> was *not measured* — or not measured *the way the label claims* — is rendered as evidence,
> and in every case the error points in the direction the author wanted. This is now five
> instances across three sprints and two authors. It is the most transferable output of this
> lane.

| # | where | the defect | direction |
|---|---|---|---|
| 1 | s24 `d_setequality_proof` REV1 | every probability family carried an epsilon floor, so the prefix-hood assertion **could not fire**; it reported 0/2916 and the claim it "confirmed" was false | toward the author's claim |
| 2 | s24 `d_harness.paired_stats` | an effect at **0.39× its own MDE** labelled `MEASURED` because a 5-cluster bootstrap CI excluded zero; had already mislabelled three D1-C rows | toward the lane's hypothesis |
| 3 | s24 `d_harness.aggregate` | `gate_equality_rate` is computed as `r.get("gate_equality", True)` — **a row that never measured equality is counted as having passed it.** `gate_pass` beside it is a strict `r["gate_pass"]`. Demonstrated: rows with the key absent return a rate of **1.0** | toward "the gate passed" |
| 4 | s24 `d_harness.write` | with `required_keys=None` **no `_COMPLETE` sidecar is written at all**, so a stale sidecar from an earlier complete run survives beside new content and keeps asserting `complete: true` | toward "the artefact is complete" |
| 5 | the `+0.113 Å` claim | a **marginal mean quoted where a paired label was required**; the paired contrast is 0.51× MDE with a CI spanning zero. Carried into two sprint reports and a published artifact | toward "the component contributes" |

A sixth, minor, and reported for completeness rather than as a member: `gate_set_equality`
widens the energy prefix by an **absolute** 1e-12, which can only make subset-hood *easier* to
pass. It is harmless on rank-standardised energies (gaps ~1/k) and vacuous on raw ones (the
ulp at 1e18 is ~100) — but the leniency points the flattering way, and on this list that is
worth writing down. Likewise `paired_stats` drops non-finite pairs silently while `aggregate`
uses `nanmean` for the means, so a mean and its contrast can rest on different target sets;
`n` is reported, so it is visible rather than hidden.

**Recommended, not applied** (I did not edit a Lane-D module): make `gate_equality` a strict
key access, and have `write` refuse `required_keys=None` or delete the stale sidecar.

---

## 6. WHAT DAMAGED MY OWN EXPECTATIONS

1. **I expected to verify four facts and find them all sound. Three were.** The fourth — the
   contribution — was wrong twice over: wrong temperature *and* an unlabelled marginal. I
   found the temperature problem while reading `core/pipeline.py` for the specification, i.e.
   **because I was writing documentation, not because I was auditing.** The audit I *had*
   planned would not have found it: I was going to re-verify the +0.113 arithmetic, and the
   arithmetic is correct.
2. **My own prediction was right for a reason I only half had.** I predicted α and T were
   substitutes acting on concentration. They are. But I predicted it from the *mechanism* and
   had no expectation that the entropy would explain the **no-circuit arms on the same curve**
   — that came from noticing that `s8/integrate.py` gives every arm the identical operator, a
   source fact, not a hypothesis. The strongest construction in this lane was found by reading
   the code that made the artefact.
3. **I nearly quoted the tail-baseline constant as a fact.** It is instrument-dependent
   (+0.5666 / +0.655634 / +0.524 across three instruments). Had my number matched, I would
   have reported a reproduction; because it did not, I had to work out that the *mechanism* is
   what reproduces. **The near miss is the useful part: a matching constant would have taught
   me nothing and I would have believed more than I should.**
4. **I drafted a sentence I could not support, flagged it myself, and it failed.** "The
   variational state attains its own analytic optimum" was an inference from an *endpoint*
   contrast dressed as a statement about *distributions* — and the endpoint contrast was
   itself underpowered at 0.24× MDE. **It is the same error as this project's L2, which I had
   read.** I flagged it rather than shipping it, measured it on the coordinator's ruling, and
   the falsifier fired: 0.902 nats, 45% of the mass. **The replacement is a better result than
   the claim was**, but I should record that the flag came from re-reading my own draft against
   the standard, not from any doubt I had while writing it.
5. **I labelled the headline UNDERPOWERED and it is the number the report wants.**
   `VQE_LFO − argmin = −0.1405` at 0.68× MDE with 5/5 folds agreeing is *exactly* the shape
   that tempts a promotion, and the rule that stops it is the rule Lane D had to fix after it
   flattered Lane D. I am reporting it as UNDERPOWERED with its fold agreement stated beside
   it, and not as a result.

---

## 7. WHAT REMAINS LIMITING, AND WHAT I AM NOT CLAIMING

* **The selector cannot be the lever, and this is now established twice over.** s24 §1 proved
  its selection is the classical top-m set; this lane shows that even the one channel the
  theorem leaves open — where the α-mass cuts, i.e. the readout's spread — is fully explained
  by a scalar entropy that four no-circuit arms trace out for free.
* **NOT MEASURED, and I will not imply otherwise:** whether the circuit would contribute on an
  energy landscape where the classical Boltzmann optimum is *not* cheap to compute. Everything
  here rests on a 128-dimensional diagonal H whose exact Gibbs state is one `exp` call. That
  is a property of this deployment, not of CVaR-VQE.
* **NOT MEASURED:** whether entanglement helps. The standing project result is s20/s21's
  −0.013 [−0.095, +0.077] with the CNOTs deleted — **NOT MEASURED, not refuted** — and I did
  not re-derive it.
* **The binding constraint on this component is the READOUT'S INSENSITIVITY, not the circuit.**
  §2.5b puts a number on it: a trained state and its own analytic optimum, differing by 45% of
  their probability mass, produce RMSDs 0.24× MDE apart. No improvement to the state can
  matter through a readout with that much slack, and that — not expressivity, not the
  optimiser, not the tail — is where the component's contribution is lost.
* **The pillar is satisfied by the algorithm being genuine, not by the tail being
  load-bearing.** That is the sentence I would ship, and §2 is the evidence for both halves.

---

## 8. ARTEFACTS

All provenance-stamped via `s24.stats_lib.save_atomic`. Nothing in `core/` or `s24/` was
edited; the two harness fixes in §5 are recommendations, not applied changes.

    s25/PREREG_Q.md              the fork lists and falsifiers, filed before the runs
    s25/QUANTUM.md               the postdoc-facing specification
    s25/q_verify.py              -> s25/results/q_verify.json     Q-A + Q-C, property
    s25/q_alpha.py               -> s25/results/q_alpha.json      Q-B, directional, n=126
    s25/q_gibbs.py               -> s25/results/q_gibbs.json      divergence + training control
    s25/q_plateau.py             -> s25/results/q_plateau.json    gradient variance
                                    s25/results/q_plateau.log     the run's full console output
