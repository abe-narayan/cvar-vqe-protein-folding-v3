# SPRINT 14 — coordinator findings

Tiering: **DEMONSTRATED** (legitimate inference-time information, CI + null) /
**ORACLE DIAGNOSTIC** (reads the native; prices a ceiling; never a headline) /
**HYPOTHESIS** / **REFUTED**.

Instrument verified at sprint start: `python -m s12.instrument` reproduces
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`.

---

## C0. The shared emission ladder exists, and it has a leakage guard

`s14/ladder.py`. An emitter is a native-free `(pdb, seq, n, fold, rng) -> (phi, psi)`.
Four agents are working in parallel and their arms only compose if every one of them emits
and is scored identically; this is that harness.

`audit_emitter` is the part worth keeping. It monkey-patches `I.load_univ` so every native
trace becomes noise, then demands bit-identical emitter output. Any emitter that reads
`nat_ca`, native torsions or the oracle `rr` column anywhere in its call graph fails
mechanically rather than being caught by review. **Every arm below passed it.** Given that
this project's single most dangerous failure mode is an oracle diagnostic becoming a
headline, the guard should be a precondition for every future arm.

## C1. Levels 0 and 1 reproduce Sprint 13 exactly — the harness is trustworthy

126 targets, three seeds on stochastic arms. **DEMONSTRATED.**

| arm | mean | median | <2 A | FAIL18 | vs incumbent | vs helix |
|---|---|---|---|---|---|---|
| incumbent synthesis pipeline | 3.204 | 2.966 | 0.28 | 6.026 | reference | -0.861 |
| L0 constant alpha-helix, zero information | 4.065 | 4.253 | 0.29 | 5.887 | +0.861 [+0.574,+1.171] | reference |
| L1a uniform random k=4 library state | 5.065 | 4.819 | 0.00 | 6.355 | +1.861 [+1.584,+2.132] | +1.000 |
| L1b sample the class back-off prior | 4.830 | 4.614 | 0.00 | 6.102 | +1.626 [+1.382,+1.872] | +0.765 |
| L1c class prior argmax (its closed-form minimiser) | 3.969 | 3.948 | 0.24 | 5.670 | +0.765 [+0.514,+1.024] | -0.096 |

The constant helix at 4.065 and the prior argmax at 3.969 match Sprint 13 to the third
decimal, and the -0.096 gap between them is the recorded 0.096 A. The harness is sound.

**A distinction that must not be lost.** Sprint 13's "uniform-random torsion states still
reach 2.698 A" is an ORACLE COORDINATE-DESCENT ceiling from a uniform start. Raw uniform
**emission** is 5.065 A. These differ by 2.4 A and it would be easy to quote the wrong one.

**New observation.** The zero-information constant helix has a *lower* FAIL18 mean (5.887)
than the incumbent (6.026), and so does the class prior argmax (5.670). This is the Sprint 12
result — a blind pipeline beats the real one on the failure class — reappearing in torsion
space through a completely different mechanism. On the 18 hardest targets, having no
information is still better than having this system's information.

## C2. The retrieval pool is a better torsion predictor than the trained sequence model — and it still loses

**DEMONSTRATED.** `s14/retprior.py`.

Sprint 13's torsion prior used the residue-CLASS back-off pool, so every GENERAL residue got
the same distribution. But each target's universe cache stores `PHI`/`PSI` for every
retrieved window, so the K=500 BLOSUM pool is already a **position-specific,
sequence-conditioned** empirical distribution over (phi_i, psi_i). Nobody had ever fed the
retrieval half's torsion information to the torsion half. I built it.

Mean absolute torsion error, ORACLE post-hoc, against the established reference points:

| source | phi | psi |
|---|---|---|
| sequence-blind corpus marginal (S13) | 36.4 | 72.8 |
| best learned full-context predictor (S13) | 36.1 | 62.4 |
| **retrieval top-75 circular mean (new)** | **33.6** | **59.2** |
| retrieval pool-500 circular mean | 34.6 | 67.1 |
| retrieval top-75 state argmax at k=4 | 34.0 | 60.7 |

**The retrieval pool beats the trained leave-fold-out sequence predictor on both angles**,
by 2.5 deg of phi and 3.2 deg of psi, with no training at all. It is the best torsion
channel measured anywhere in this project. The pool concentrates well on phi and poorly on
psi: mean resultant length R_phi 0.834, R_psi 0.645, normalised k=4 state entropy 0.548.

And it does not matter. Emitted accuracy:

| arm | mean | vs incumbent |
|---|---|---|
| L2d top-75 similarity-weighted circular mean | 3.514 | +0.310 [+0.128,+0.492] |
| L2e top-75 state argmax, k=4 | 3.842 | +0.638 [+0.401,+0.884] |
| L2f pool-500 state argmax, k=4 | 4.004 | +0.800 [+0.525,+1.080] |
| L2b top-75 circular mean | 4.072 | +0.868 [+0.551,+1.202] |
| L2c top-20 circular mean | 4.124 | +0.920 [+0.591,+1.287] |
| L2a pool-500 circular mean | 4.175 | +0.971 [+0.643,+1.332] |

Every CI excludes zero. **No native-free torsion emitter yet built beats the retrieval
pipeline that consumes the very same windows.** The pipeline's coordinate average of whole
structures extracts more than any per-residue torsion summary of the same pool.

Mechanism, and it follows from geometry rather than from tuning: torsion errors compound
along the chain, coordinate averaging does not. Averaging in torsion space commits every
residue's error into the integration; averaging in coordinate space lets independent errors
cancel in place. This is a concrete, structural reason the retrieval architecture has been
hard to beat, and it was not previously stated.

## C3. Mean absolute torsion error does not order emitted accuracy — REFUTES the axis the restraint surface is indexed on

**DEMONSTRATED** for the ordering; the surface correction is **ORACLE DIAGNOSTIC**.

Two arms built from the same retrieved windows come out backwards:

| arm | phi MAE | psi MAE | emitted |
|---|---|---|---|
| L2b top-75 circular mean | 33.6 deg | 59.2 deg | 4.072 A |
| L2d top-75 similarity-weighted | 38.5 deg | 62.4 deg | **3.514 A** |

4.9 degrees worse on phi, 3.2 worse on psi, and 0.56 A **better** in structure. So sigma
does not order emitted accuracy even within one family of predictors on one target set.

This matters well beyond a curiosity, because the sprint's most attractive route — the
restraint surface at sigma x coverage, from which "sigma 12 deg reaches 1.486 A" and "the
incumbent is equivalent to sigma ~29 deg" are both read — corrupted native torsions with
**i.i.d. Gaussian noise**. Real predictors do not have i.i.d. errors. A retrieval prior
inherits the systematic bias of its library; a shift-based predictor inherits the systematic
bias of its training set; both produce errors correlated along the chain.

What a backbone builder integrates is the *cumulative* error, so at matched marginal sigma:
coherent errors become systematic curvature and the chain leaves the native monotonically,
while anti-correlated errors cancel almost exactly. This is the Sprint 12 error-coherence
law — at identical accuracy, coherent mistakes emit +0.291 A and i.i.d. mistakes -0.142 A —
reappearing in torsion space through chain integration rather than through a corrector.

`s14/coherence.py` measures the real emitters' error autocorrelation and bias, then rebuilds
the restraint surface with coherence as an explicit second coordinate, at matched sigma
across i.i.d. / constant-bias / AR(+0.7) / AR(-0.7) / alternating error fields. Until that
lands, **every sigma threshold in this project should be read as conditional on an i.i.d.
error model that no real predictor satisfies**, and the SHIFT agent's kill threshold must be
re-derived against the coherence mode its actual method produces.

## C4. Coherence is worth as much as accuracy — but real predictors are already in the good regime. My own hypothesis REFUTED.

`s14/coherence.py`, 126 targets, 12 replicates per cell, every arm rescaled to the SAME
marginal sigma so only the along-chain correlation structure differs.

**ORACLE DIAGNOSTIC.** Emitted RMSD at matched sigma, by error-coherence mode:

| sigma | AR(-0.7) | alternating | i.i.d. | constant bias | AR(+0.7) |
|---|---|---|---|---|---|
| 10 deg | 1.094 | 1.258 | 1.382 | 1.565 | 1.911 |
| 12 deg | 1.272 | 1.485 | 1.617 | 1.772 | 2.203 |
| 15 deg | 1.542 | 1.826 | 1.985 | 2.103 | 2.709 |
| 20 deg | 1.981 | 2.347 | 2.528 | 2.586 | 3.237 |

At sigma 12 deg the spread across coherence modes is **0.93 A at identical marginal
accuracy**. Expressed as the sigma needed to reach 2.0 A: AR(-0.7) 20.2 deg, alternating
16.7 deg, i.i.d. 15.1 deg, bias 14.1 deg, AR(+0.7) 10.6 deg. **Coherence is worth nearly a
factor of two in required angular accuracy** — it is not a second-order correction.

**MY HYPOTHESIS WAS THAT REAL PREDICTORS HAVE COHERENT ERRORS. IT IS REFUTED.** Measured
lag-1 error autocorrelation for all ten native-free emitters built so far lies between
**-0.092 and +0.066** — every one of them is in the near-i.i.d. band, and if anything
mildly anti-correlated, which is the *cheap* direction:

| emitter | RMSD | phi MAE | phi lag-1 | phi bias |
|---|---|---|---|---|
| L2d top-75 similarity-weighted | 3.514 | 38.5 | -0.051 | +4.8 |
| L2e top-75 state argmax | 3.842 | 34.0 | -0.035 | +3.6 |
| L1c class prior argmax | 3.969 | 36.5 | -0.028 | +10.5 |
| L0 constant alpha-helix | 4.065 | 39.2 | -0.029 | +19.6 |
| L2b top-75 circular mean | 4.072 | 33.6 | -0.032 | +6.2 |

Consequences, and the second one matters more than the first:

1. **The Sprint 12 restraint surface is approximately valid.** It used i.i.d. corruption,
   and real predictors are i.i.d.-like, so its sigma thresholds stand. My warning in C3 that
   every sigma threshold was priced on the wrong error model is **withdrawn**. The
   i.i.d. surface is if anything mildly conservative. This correction runs in the same
   favourable direction as the Sprint 13 terminal-dropout correction.
   The incumbent's i.i.d. sigma-equivalent measures **27.1 deg** here, a mild refinement of
   the recorded ~29 deg.

2. **Coherence therefore does NOT explain the C3 anomaly.** L2b and L2d have essentially
   identical lag-1 autocorrelation and similar bias, yet 4.9 deg of MAE separates them in
   one direction and 0.56 A of RMSD in the other. The anomaly stands, and the remaining
   candidate is *where* along the chain the errors fall — an error at residue 2 rotates the
   whole downstream chain, the same error at residue 12 moves one residue. `s14/position.py`
   tests it. It can fail.

3. **A coherent-error predictor would be actively dangerous.** Any future method whose
   errors are positively autocorrelated needs sigma 10.6 deg where an i.i.d. one needs 15.1.
   Error coherence should therefore be reported alongside sigma for every torsion channel,
   and the chemical-shift route must report the coherence mode of its actual method.

## C5. The September 2026 preprint is real, it does not threaten the direction, and it corroborates our hardest negative

**LITERATURE-SUPPORTED.** Full ledger and verification status in `s14/lit_FINDINGS.md`
(1,115 lines, primary text fetched and read in full).

The paper is **arXiv:2609.02113**, Cumbo et al. (Cleveland Clinic), submitted 2 September
2026, method **QTF** — a logarithmic-scale VQE for off-lattice protein structure prediction
in continuous torsional space, with chignolin and Trp-cage.

Five things decide our response:

1. **Their headline is oracle-selected.** The advertised 0.623 A is the argmin over 2.68 M
   snapshots *ranked by RMSD to the native*. Their energy-selected result — the predictive
   one — is a median final model of **2.90 A on chignolin and 5.39 A on Trp-cage**, on two
   targets, with an RMSD that excludes terminal residues and is never defined in Methods.
   Our 3.213 A across 126 targets, full-chain, no oracle, is competitive or better. This is
   precisely the oracle-becomes-headline failure this project tiers against. Do not chase
   0.623.

2. **They independently corroborate our central negative.** They find AMBER/OpenMM total
   potential energy **negatively correlated with RMSD on both proteins**, over 7.3 M
   structures — a different group, different targets, same direction as our -0.088. Their
   own tables show snapshot mining finds a sub-2 A structure in ~99% of replicas while their
   energy functions select one in **6% / 4% / 0.3%** of cases. They print no numeric rho, no
   CI and no p-value in 50 pages. Our certified global optimum being +0.139 A worse than
   random on nine fully enumerated targets is strictly stronger and remains unpublished.

3. **There is no VQE in the technical sense and no classical control.** Their energy is
   never a qubit observable — it is computed classically from the reconstructed Cartesian
   structure — and the optimisers are COBYLA/SLSQP. Their parameter count is
   `P = 2n(ceil(N/n)+3) ~ 2N+6n`: **132 classical parameters to produce 43 torsions** for
   chignolin, so the classical search space is about three times *larger* than the torsion
   vector it encodes. The logarithmic-qubit claim is a reparameterisation, not a
   compression. Under this project's rule 1.4 we could not call that VQE. **Do not adopt
   the phase encoding** — inheriting P > N would be a strict architectural regression.

4. **CVaR is entirely absent** from the paper. In-loop CVaR in continuous torsion space
   survives our novelty audit. CVaR-VQE on lattices is taken.

5. **A defect they did not notice, which we have independently measured the hazard for.**
   Their hardware decoder is `theta_i = 2*pi*C_i - pi` with `C_i` a cumulative probability,
   so the decoded torsion vector is **always monotonically sorted with total variation
   <= 2*pi**. That biases it toward near-constant torsion profiles, i.e. regular secondary
   structure — which predicts their chignolin-works/Trp-cage-fails asymmetry from the
   decoder alone. This is the same trap as our constant alpha-helix beating uniform random
   by 0.457 A, and C1 shows the helix is a strong attractor in exactly this way.

**Verdict: no experimental change, one repositioning.** Nothing invalidates the Sprint 14
direction, and the corroboration de-risks our hardest result. The contribution should be
positioned not as a new quantum encoding — that ground is now taken and is worth nothing —
but as **the first properly controlled measurement of whether quantum conformational
optimisation helps, on a certified landscape**. What survives our novelty audit: the exact
locality theorem, the no-free-parameter spectrum-to-gradient-variance chain, the certified
global optimum on an enumerable space, in-loop CVaR in torsion space, the causal
VQE-versus-classical control at matched budget, and the protected statistical instrument.

One open, unverified and consequential item: their custom energy may have an end-to-end
distance bias whose target could be native-derived, which if true would make their
custom-energy chignolin result leaked. Checkable in their source; not a sprint dependency.

## C6. The cost of a torsion error is a symmetric hump peaked at mid-chain — and it explains a Sprint 13 correction

**ORACLE DIAGNOSTIC.** `s14/position.py`, 126 targets, one residue perturbed at a time at
sigma 20 deg, 24 replicates per residue, everything else held native.

Normalised cost weight by fractional chain position, phi and psi combined, N-terminus first:

| 0.0-0.1 | 0.1-0.2 | 0.2-0.3 | 0.3-0.4 | 0.4-0.5 | 0.5-0.6 | 0.6-0.7 | 0.7-0.8 | 0.8-0.9 | 0.9-1.0 |
|---|---|---|---|---|---|---|---|---|---|
| 0.012 | 0.060 | 0.103 | 0.144 | **0.170** | **0.169** | 0.145 | 0.112 | 0.070 | 0.015 |

The middle 40% of the chain carries **62.8%** of the total single-residue cost; the outer
40% carries **15.7%**. That is a **4.0x ratio**, and the terminal tenth at each end is about
**14x cheaper** than mid-chain.

The mechanism is geometric, not empirical. A torsion at position p hinges two rigid segments
of length p and n-p against one another, and under Kabsch superposition the resulting
displacement scales with the lever-arm product p(n-p), which is maximised at mid-chain and
vanishes at both termini. The measured curve is symmetric to within its own noise, exactly
as that predicts.

**This quantitatively explains the Sprint 13 correction that terminal dropout is 0.40-0.50 A
CHEAPER than uniform dropout.** That correction was empirical and surprising; it is now
mechanical. Terminal torsions carry ~1.4% of the cost each against ~17% mid-chain, so
dropping them is nearly free.

**Actionable consequence.** Raw coverage percentage is the wrong axis for any partial-torsion
channel. Coverage must be **position-weighted**: 75% coverage concentrated mid-chain beats
90% concentrated at the termini. The existing "reaches 2.0 A only at >= 90% coverage"
threshold was derived on raw fraction and needs re-deriving on the weighted axis. Sent to
the SHIFT agent, since chemical-shift predictors decline precisely at the termini — the
cheap region — so the channel is more attractive than a raw-coverage reading suggests.

**MY POSITIONAL HYPOTHESIS WAS ALSO WRONG, AND IN A SPECIFIC WAY WORTH RECORDING.** I
predicted N-terminal dominance, on the reasoning that an early error rotates the whole
downstream chain. Measured front/back ratio is **0.886** — if anything the C-terminal half
costs slightly more. Optimal superposition is what defeats the prediction: Kabsch is free to
rotate the whole structure, so there is no privileged end, and the hinge argument replaces
the lever argument.

**And it does not rescue C3.** Re-scoring the ten emitters with these position weights makes
the ordering *worse*, not better: Spearman(error, emitted RMSD) across emitters is **+0.467
for plain MAE and +0.370 position-weighted**. Emitter errors are close to uniform across
positions, so re-weighting adds noise without signal.

**C3 is therefore softened and partly withdrawn.** Plain MAE does order emitters positively
overall (rho +0.467). The top-75 circular mean versus similarity-weighted pair is a local
inversion, not a global failure of the axis. The honest statement is that MAE is a
positively but imperfectly correlated summary with real inversions, so it should not be used
alone to choose between two candidate torsion channels — and the specific inversion that
prompted C3 remains **unexplained** after testing both coherence and position.

## C7. A native-free STRUCTURAL Hamiltonian is eight times better ordered than any physical energy — and its argmin still loses

**DEMONSTRATED.** `s14/hamil.py`, all 126 targets, 4,000 native-free configurations each
(half from the retrieval prior, half uniform), both terms standardised per target.

The objective is deliberately not a molecular energy. It sums two native-free channels the
project already owns, one local and one non-local:
a retrieval-conditioned per-residue torsion prior (1-local, exactly diagonal) and the
shipped leave-fold-out distogram's Bayes-risk score over CA-CA pair distances (many-body,
since by the Sprint 13 locality theorem `d_ij` depends on the `j-i-1` residues between i
and j). `w` is the share on the distogram term.

| w | rho all | **rho low decile** | argmin RMSD | decile mean | FAIL18 | vs incumbent |
|---|---|---|---|---|---|---|
| 0.00 | 0.306 | 0.351 | 3.809 | 3.956 | 6.537 | +0.605 [+0.374,+0.844] |
| 0.25 | 0.390 | **0.370** | 3.605 | 3.740 | 6.268 | +0.401 [+0.282,+0.522] |
| 0.50 | 0.469 | 0.317 | 3.576 | 3.675 | 6.156 | +0.372 [+0.267,+0.482] |
| 0.75 | 0.523 | 0.271 | **3.575** | 3.651 | 6.106 | +0.371 [+0.252,+0.491] |
| 1.00 | 0.539 | 0.109 | 3.615 | 3.663 | 5.989 | +0.411 [+0.263,+0.562] |

Reference: incumbent 3.204; ORACLE best of the same 4,000 sampled configurations **1.683**.

Two things follow, and the second is the important one.

**The objective problem is solvable.** In-decile rank correlation reaches **+0.370**, against
Legacy's **+0.043** and raw AMBER's **-0.088**. That is roughly an eightfold improvement in
the one statistic that governs where a search actually lives, and it was obtained without
learning anything — only by choosing structural channels over energetic ones. The two terms
are genuinely complementary: global ordering rises monotonically with distogram weight
(0.306 to 0.539) while in-decile ordering peaks at a mixture (w = 0.25) and collapses at
pure distogram (0.109). The distogram orders the whole space; the torsion prior orders the
good region.

**And it does not help.** The best argmin is 3.575, still **+0.371 [+0.252,+0.491] worse
than the incumbent**, with every CI excluding zero. The samples contain a 1.683 A structure
and the best-ordered objective this project has ever built selects a 3.575 A one — a
**selection gap of 1.9 A** that an eightfold gain in rank correlation barely dented.

This generalises the Sprint 13 negative from "the physical energies are bad objectives" to
something much stronger: **even a well-ordered native-free objective does not make torsion
search competitive with the retrieval pipeline.** Whether that is a search failure or a
ceiling is exactly what `s14/budgetcurve.py` is now measuring.

## C7-CORRECTION. The control landed and it substantially cuts my own headline

**The "eightfold better ordered" claim in C7 is WITHDRAWN as stated.** It was measured under
a proposal that drew half its configurations from the very prior the w=0 arm scores.
`s14/hamil_control.py` re-ran the identical analysis under three proposals and a
sequence-blind twin, 126 targets, 2,000 configurations each:

| proposal | w | rho all | **rho low decile** | argmin |
|---|---|---|---|---|
| mixed (used in C7) | 0.00 | 0.304 | **0.355** | 3.905 |
| mixed | 0.50 | 0.469 | 0.346 | 3.615 |
| mixed | 1.00 | 0.542 | 0.171 | 3.585 |
| **uniform** | 0.00 | 0.071 | **0.046** | 4.551 |
| **uniform** | 0.50 | 0.324 | 0.120 | 4.079 |
| **uniform** | 1.00 | 0.442 | **0.161** | 3.913 |
| prior | 0.00 | 0.211 | 0.401 | 3.833 |
| prior | 0.50 | 0.386 | **0.416** | 3.553 |
| mixed, SEQUENCE-BLIND twin | 0.00 | 0.105 | 0.073 | 4.228 |

**Under a uniform proposal the torsion-prior term's in-decile correlation collapses from
0.355 to 0.046** — statistically indistinguishable from Legacy's +0.043. The headline number
was substantially measuring agreement between the proposal and the scorer.

The corrected claims:

- **The honest in-decile figure is 0.161**, from the distogram term under a uniform
  proposal, against Legacy's +0.043 and raw AMBER's -0.088. That is roughly **fourfold**
  better than the best physical energy, not eightfold, and it comes from the distogram
  rather than from the torsion prior.
- **The torsion prior's ordering skill is conditional on being in the prior-typical region.**
  Under uniform sampling it is ~0.046; under prior sampling it is 0.401. Both are real
  measurements of different things: sampled uniformly you are mostly in garbage, where prior
  probability separates garbage from garbage. This conditionality should be stated whenever
  the prior term is used, not averaged away.
- **The sequence-blind twin is the reassuring half.** Replacing the retrieval-conditioned
  prior with the class back-off prior drops in-decile from 0.355 to **0.073** and the decile
  mean from 3.958 to 4.658. So the retrieval conditioning carries genuine target-specific
  information and is not merely generic Ramachandran — unlike the Sprint 13 library, 88% of
  whose value was generic.

**What is unaffected.** The argmin conclusion is unchanged and if anything strengthened:
across all three proposals and both twins the argmin lands between 3.54 and 4.55, every one
of them worse than the incumbent 3.204. And **C8's saturation finding stands**, because that
experiment used one proposal consistently throughout and its claim is about the *shape* of
the budget curve, not the absolute rank correlation. C8's prose describing the objective as
"eight times better ordered" should be read as "roughly fourfold, and distogram-driven".

## C8. Searching a GOOD objective harder buys nothing — search saturates at 300 evaluations while the selection gap grows to 2.0 A

**DEMONSTRATED** (native-free search) with an **ORACLE DIAGNOSTIC** reference column.
`s14/budgetcurve.py`, 126 targets, w = 0.25, one pooled sample of 20,000 configurations per
target read as increasing prefixes so every budget is a strict subset of the next.

| budget | objective found | RMSD emitted | <2 A | ORACLE best available | **selection gap** | space seen |
|---|---|---|---|---|---|---|
| 10 | -1.108 | 3.742 | 0.22 | 2.907 | 0.835 | 3.9e-06 |
| 30 | -1.224 | 3.656 | 0.24 | 2.530 | 1.126 | 1.2e-05 |
| 100 | -1.297 | 3.596 | 0.24 | 2.232 | 1.364 | 3.9e-05 |
| **300** | -1.349 | **3.565** | 0.25 | 1.987 | 1.577 | 1.2e-04 |
| 1,000 | -1.391 | 3.599 | 0.25 | 1.818 | 1.781 | 3.9e-04 |
| 3,000 | -1.409 | 3.632 | 0.24 | 1.697 | 1.935 | 1.2e-03 |
| 10,000 | -1.428 | 3.588 | 0.24 | 1.580 | 2.007 | 3.9e-03 |
| 20,000 | -1.436 | 3.572 | 0.25 | 1.532 | **2.040** | 7.9e-03 |

Read the three columns against each other, because separately each one misleads.

**The search works.** The objective improves monotonically and without exception,
-1.108 to -1.436. Anyone reporting only the objective axis would call this a success.

**The candidate pool improves too.** More sampling genuinely puts better structures in
reach: the best available falls 2.907 to 1.532.

**And the emitted structure does not move.** 3.742 at ten evaluations, a minimum of 3.565 at
three hundred, 3.572 at twenty thousand. **A 2,000-fold increase in search buys 0.17 A**, and
all of it arrives by budget 300. Past the minimum the degradation is +0.007 A — nothing.

**The selection gap is the mechanism, and it grows monotonically: 0.835 to 2.040 A.** The
objective's ability to discriminate saturates at about 3.57 A no matter how good the
candidates become. This is a **ceiling of the objective, not of the search**.

**This refines the Sprint 13 result rather than repeating it.** Sprint 13, on Legacy
(in-decile rho +0.043): 3.764 at 10 evaluations, 3.667 at 300, **3.920** at the certified
global optimum — searching harder actively *hurt*. Sprint 14, on a better-ordered objective:
3.742, 3.565, 3.572 — searching harder no longer hurts, and still does not help. **The
objective improvement bought stability, not accuracy.**

**CORRECTION to how that sentence was first written.** I originally called this "an eightfold
in-decile improvement". Two separate problems, both now fixed. The eightfold figure came
from a confounded proposal (see C7-CORRECTION; the honest figure is roughly fourfold and
distogram-driven). And the VQE agent showed that **in-decile rank correlation is not monotone
in objective quality** — on a signal-tunable family it peaks near +0.446 at moderate signal
and falls to +0.340 where the objective is genuinely much better (global rho +0.933), because
a better objective's lowest decile is a narrower and more homogeneous set with less orderable
spread left in it. So part of any in-decile number is ceiling effect rather than skill, and a
single in-decile figure is ambiguous between objectives differing fourfold in real quality.
Both axes are now required in the shared brief. On both axes this objective is: global rho
+0.442 to +0.542, in-decile +0.161 to +0.171. The safe claim is the narrow one — **at this
objective quality the structural readout has stopped responding to search.**

The Sprint 13 fraction-of-space law reproduces exactly, which is a good internal check: the
20 targets that see >= 1% of their space degrade +0.156 A past their minimum, while the 74
that see < 0.1% degrade +0.016 A.

**Consequence for the quantum question, and it is the sprint's central one.** A variational
optimiser is a better search. Search is already saturated at 300 evaluations out of 20,000
on the best native-free objective this project can build. **There is therefore no room for
VQE, CVaR, a better ansatz, a better encoding or a better initialisation to contribute
through search quality in this representation** — whatever they optimise, the structural
readout stopped responding two orders of magnitude earlier.

### C8-CORRECTION. Half of this reproduces exactly and half of it does not. The saturation claim is NOT general.

The VQE agent recomputed **my** statistic on **their** enumeration — uniform draws read as
nested prefixes, same objective, same weight, 6 seeds x 9 targets, independent
implementation, disjoint target set.

| budget | 10 | 300 | 3,000 | 20,000 |
|---|---|---|---|---|
| their returned | 3.297 | 3.173 | 3.177 | **2.920** |
| their ORACLE best available | 2.660 | 1.590 | 1.206 | 0.997 |
| **their selection gap** | 0.637 | 1.583 | 1.971 | 1.923 |
| **my selection gap (126 targets)** | 0.835 | 1.577 | 1.935 | 2.040 |

**THE SELECTION GAP REPRODUCES AT EVERY BUDGET** — to 0.006 A at budget 300 and 0.036 A at
3,000, and within 0.12 A everywhere. Same statistic, disjoint target sets, independent
implementations. That is stronger corroboration than the "within 0.16 A" I first wrote, and
the discrimination finding should be stated at that strength.

**BUT THE SATURATION HALF DOES NOT REPRODUCE, AND MY CLAIM WAS TOO GENERAL.** My returned
column is flat (3.742 to 3.572, 0.17 A of movement). **Theirs falls 0.38 A monotonically**
(3.297 to 2.920). On their nine targets, searching harder genuinely does improve the emitted
structure.

The likely cause is chain length and hence fraction of space explored: their nine are all
n=9, the length bin Sprint 13 identified as the most thoroughly explored, and at budget
20,000 they have seen **7.6%** of their space where my n=9-16 targets have seen a small
fraction of a much larger one. That is the same fraction-of-space law that governed the
Sprint 13 budget curve, reappearing.

**The honest joint statement, which replaces the one above:**

> The discrimination loss is the same everywhere and reproduces exactly. Whether search
> still buys anything on top of it is **target-length dependent** — it does on short chains
> and does not on the full instrument.

"The emitted structure does not move" must not be written as a general claim. It is true on
the 126-target instrument and false on the enumerated nine, and I would have written it as
general without this check.

**What survives unweakened** is the part that matters for the quantum question, because it
does not rest on the flat curve: a ~1.9 A selection gap that reproduces across both target
sets and **survives the certified optimum on the enumerated nine**, where no optimiser can do
better by construction.

### C8b. VQE is a WORSE optimiser than 1-opt greedy, and the structural tie is not a tie

From the Part C grid, on the OBJECTIVE axis, fraction of cells reaching the certified optimum
at budget 30,000:

| method | fraction reaching the certified optimum |
|---|---|
| greedy 1-opt | 0.71 |
| simulated annealing | 0.71 |
| **VQE** | **0.00 - 0.12** |

At signal 0 (the real Legacy objective) greedy and annealing reach the certified optimum in
**100%** of cells and VQE in **0%**.

So VQE is decisively beaten on the thing it is supposed to be good at, while being tied with
those methods on structural accuracy. **That tie is not VQE holding its own.** It is both
methods hitting the same discrimination floor from different distances — the floor is close
enough that a much worse optimiser reaches it too. Reporting only the structural axis would
have made VQE look competitive; reporting both makes it clear it is not.

That closes the search branch. It leaves exactly one branch open, and it is the branch the
project's own history points at: stop reading the argmin, and aggregate a low-energy SET in
coordinate space (`s14/consensus.py`, running). The Sprint 12 terminal-operator law says the
output tracks the set mean rather than the set best, and C2 says coordinate averaging is
where torsion errors cancel.

## C9. Aggregating a low-energy SET recovers 0.29 A over the argmin — and still loses

**DEMONSTRATED.** `s14/consensus.py`, 126 targets, 4,000 native-free configurations each,
w = 0.25, structural Hamiltonian used to rank and the top m coordinate-averaged.

| m | set mean | set best | consensus | <2 A | vs incumbent |
|---|---|---|---|---|---|
| 1 (argmin) | 3.608 | 3.608 | 3.608 | 0.25 | +0.404 [+0.279,+0.536] |
| 20 | 3.608 | 3.066 | 3.370 | 0.27 | +0.166 [+0.077,+0.253] |
| 100 | 3.648 | 2.775 | 3.317 | 0.27 | +0.112 [+0.010,+0.212] |
| **200** | 3.684 | 2.647 | **3.314** | 0.27 | **+0.110 [+0.004,+0.214]** |
| 500 | 3.755 | 2.468 | 3.326 | 0.27 | +0.122 [+0.007,+0.239] |

Reference: incumbent 3.204, ORACLE best of the same 4,000 sampled 1.647.

The mechanism works exactly as C2 predicted. **Set members do not improve** — the set mean is
flat at ~3.6 and actually *worsens* with m as the ranking reaches further down — yet the
consensus improves from 3.608 to 3.314. All 0.29 A of that is error cancellation in
coordinate space, obtained from members that are individually no better.

**And it still loses**, +0.110 [+0.004,+0.214], a CI that barely excludes zero. The optimal
set size is m* = 200, which by the Sprint 12 law (m* shrinks 500 -> 75 -> 20 -> 3-5 as the
objective improves) is itself a readout that this objective is still mediocre.

## C10. The projection onto ideal geometry costs 0.157 A — the price of physical validity, priced for the first time

**DEMONSTRATED.** `s14/avgspace.py`, 126 targets, the same top-75 windows in every arm, so
the candidate set is held constant and only the aggregation operator moves.

| arm | mean | median | FAIL18 |
|---|---|---|---|
| A coordinate average then project (the incumbent route) | 3.205 | 2.966 | 6.034 |
| B coordinate average, no projection | **3.048** | 2.837 | 5.832 |
| C torsion circular mean, built | 4.072 | 3.552 | 6.743 |
| D torsion circular mean, then projected | 4.072 | 3.552 | 6.743 |

Arm A reproduces the incumbent to 0.001 A, which validates the reconstruction. Three
contrasts follow, and each isolates one thing:

**The averaging space alone is worth +1.024 A [+0.698,+1.365], 41W/85L.** With no projection
on either side, coordinate averaging beats torsion averaging of the identical windows by
over an angstrom. **C2's mechanism is now a controlled measurement rather than an
inference**: torsion errors compound through chain integration, coordinate errors cancel in
place, and that is worth more than every torsion-channel improvement in this sprint combined.

**Projection is exactly a no-op on a torsion-built structure**, -0.000 [-0.000,+0.000]. That
is the internal consistency check it should be: a structure built from torsions already lies
exactly on the ideal-geometry manifold the projection targets.

**Projection COSTS 0.157 A [+0.124,+0.191] on a coordinate average, 18W/108L.** The raw
coordinate average at 3.048 beats the shipped pipeline's 3.205, and loses on only 18 of 126
targets.

That last one needs its caveat stated plainly, because it is not a free win: the raw
coordinate average is **not a physically valid backbone** — its bond lengths and angles are
whatever averaging produced. The pipeline projects for a reason. What the number actually
buys is a price tag: **imposing ideal geometry costs 0.157 A of accuracy**, and this project
had never measured it.

The consequence is structural and bears on the whole sprint. **A torsion representation IS
the ideal-geometry manifold by construction**, so every torsion-space method — including any
VQE built on one — pays that 0.157 A before it starts. It is a floor, not an artefact of the
pipeline, and it should be subtracted from the incumbent when judging what a torsion method
must beat: the fair torsion-space target is 3.205, not 3.048.

### C10-CORRECTION. The mechanism is CONTRACTION, not ideal geometry — and an inference of mine was wrong

The ENER agent built the experiment I asked for and it corrected me twice. Both corrections
matter and both are theirs.

**First, my inference was invalid.** I wrote to them that their finding of inert AMBER
refinement (+0.014 A) meant "AMBER cannot recover the 0.157 A either". It does not license
that. **Every structure in their refinement experiment was already on the ideal-geometry
manifold, where AMBER has nothing to fix.** The case at issue is one where the input is
*off*-manifold, which is a different experiment — and is exactly why it was worth running.
My inference is withdrawn.

**Second, and more important: the 0.157 A is not the price of ideal geometry.** Measured on
the averaged full backbone, coordinate averaging **contracts the structure by 16.0%**:

| bond | averaged | correct |
|---|---|---|
| N-CA | 1.238 A | 1.458 A |
| CA-CA | 3.308 A | 3.804 A |

Averaging points that disagree pulls them toward their centroid, so the raw coordinate
average is a systematically shrunken molecule. On that reading the projection's cost is the
price of **re-expanding a contracted structure**, and **any operation that restores valid bond
lengths must pay it** — the ideal-geometry manifold being incidental.

**CONFIRMED AT n=126.** The estimate was unstable in flight — rho = -0.829 at n=6, **-0.224
at n=10**, and **-0.552 at n=126**. Only the last is quotable, and the ENER agent said so
before reporting it. At full n the mechanism holds: **rho(contraction, per-target projection
cost) = -0.552**, and the contraction is far larger than the pilot suggested — **-25.8%**
(N-CA 1.120 against 1.458; C-N 0.864 against 1.329).

So the corrected statement, now demonstrated rather than hypothesised:

> The 0.157 A is **not the price of ideal geometry. It is the price of restoring ANY valid
> bond length to a structure that averaging shrank by a quarter.** The coordinate average
> scores 3.048 A because it is not a molecule; the moment it is required to be one, by either
> route, it costs the same.

That is a better explanation than mine and it changes the claim's scope. The honest statement
is: *coordinate averaging buys accuracy partly by contracting the molecule, and 0.157 A of
that apparent accuracy is not recoverable once physical bond lengths are restored, by any
means.* The number stands; the mechanism I attached to it does not.

### FINAL at n=126, CI 0.03 A wide

| arm | RMSD | vs A | CI | W/L | geom dev | Ramachandran | clashed |
|---|---|---|---|---|---|---|---|
| A: coordinate average then project | 3.205 | — | — | — | 0 by construction | 0.734 | 0.000 |
| B: raw coordinate average | 3.050 | **-0.155** | [-0.191,-0.122] | 99/27 | 0.2480 | 0.836 | **0.452** |
| **B + AMBER at restraint k=10** | **3.207** | **+0.001** | **[-0.014,+0.016]** | 53/73 | **0.0134** | **0.860** | **0.000** |
| B + AMBER k=2 | 3.248 | +0.043 | [+0.020,+0.066] | 51/75 | 0.0152 | 0.863 | 0.000 |
| B + AMBER k=0, free relaxation | 3.477 | +0.272 | [+0.177,+0.370] | 31/95 | 0.0170 | 0.868 | 0.000 |

Their instrument reproduces both of my arms independently: **A at 3.205 exactly**, and B at
3.050 against my 3.048, with the projection cost coming out **-0.155 [-0.191,-0.122], 99W/27L**
against my +0.157 at 108L. Two implementations, same numbers.

**Ideal geometry was not the binding constraint.** A genuine force field — a strictly *weaker*
constraint than the manifold — recovers exactly the same accuracy as the projection: none,
+0.001 with a CI 0.03 A wide.

**Free relaxation is the worst arm (+0.272), and that is the mechanism's signature.**
Unrestrained, the field re-expands to its own preferred geometry and drifts **1.526 A** away,
discarding the ensemble information that made arm B good in the first place.

**This is outcome 2 of the three I pre-stated: ideal geometry was NOT the binding constraint.**
AMBER recovers physical validity at *exactly* the same RMSD cost as the geometric projection.
So the 0.157 A is the price of restoring valid bond lengths to a contracted average, and it is
**irreducible by better physics** — my hypothesis' consequence is confirmed even though its
causal correlation remains unproven.

**But the geometry column contains a genuine positive that RMSD hides, and it is the sprint's
one deployable recommendation.** At equal RMSD, AMBER's output is the better structure:
geometric deviation **0.0134** against the raw average's 0.2480, **zero clashes against 45.2%**,
and **Ramachandran-allowed 0.860 against the PROJECTION's 0.734.**

> If you want a physically valid structure, **AMBER-relax the coordinate average rather than
> project it**: identical accuracy, strictly better stereochemistry, about 12 seconds per
> target. It buys a better structure, never a better number.

That is the clearest instance this sprint of why the brief demanded a geometry audit beside
every RMSD — the accuracy metric is completely blind to it.

### C10-FRONTIER. The frontier does NOT turn — and my "never a better number" is wrong

**The ENER agent's own hypothesis is REFUTED, and so is my closing statement.** I wrote that
AMBER relaxation buys "a better structure, never a better number". Tighter restraint keeps
buying accuracy. 50 targets, arm A 3.171, arm B raw 3.017:

| arm | RMSD | vs A | CI | W/L | geom dev | Ramachandran | clashed |
|---|---|---|---|---|---|---|---|
| B raw (invalid) | 3.017 | -0.154 | [-0.210,-0.099] | 36/14 | 0.2278 | 0.849 | 0.380 |
| B + AMBER k=10 | 3.169 | -0.001 | [-0.023,+0.019] | 17/33 | 0.0130 | 0.854 | 0.000 |
| **B + AMBER k=30** | 3.146 | **-0.025** | **[-0.046,-0.005]** | 26/24 | 0.0181 | **0.883** | 0.000 |
| B + AMBER k=100 | 3.101 | **-0.070** | [-0.098,-0.043] | 35/15 | 0.0374 | 0.894 | 0.000 |
| B + AMBER k=300 | 3.066 | **-0.105** | [-0.140,-0.073] | 40/10 | 0.0665 | 0.857 | 0.000 |

At k=300 AMBER recovers **68% of the projection cost**, CI excluding zero, 40W/10L, **zero
clashes**, Ramachandran 0.857 against the projection's 0.676.

**But the honest question becomes how much of that is bought with strain — and they built a
MEASURED yardstick rather than asserting one.** Free relaxation at k=0 is ff14SB's *own*
equilibrium geometry, and it sits at 0.0170 deviation from the ideal-builder constants.
Against that internal standard:

| restraint | geom deviation | relative to the force field's own equilibrium |
|---|---|---|
| k=10 | 0.0134 | **tighter** than equilibrium |
| **k=30** | **0.0181** | **essentially exactly at it** |
| k=100 | 0.0374 | 2.2x strained |
| k=300 | 0.0665 | 3.9x strained |

**So the defensible claim is k=30: -0.025 A [-0.046, -0.005] at geometry equal to the force
field's own equilibrium.** Beyond that, RMSD is being bought with strain — precisely the
"an RMSD gain bought with invalid geometry is not a gain" trap the brief demanded be guarded
against, and they caught it in their own favourable result.

### AND THE CONCENTRATION CHECK REVERSES IT. Both of us were wrong about the W/L.

The ENER agent asserted — and I published, on their word, in the user-facing report — that
26W/24L with a CI excluding zero is the *good* pattern, a small consistent shift rather than a
few large wins. **Neither of us had run the drop-top curve. It says the opposite.** Same 50
targets, k=30:

| statistic | value |
|---|---|
| mean difference | -0.0248 [-0.0462, -0.0046], 26W/24L |
| **median** | **-0.0089** |
| drop-top-1 | -0.0208 |
| drop-top-3 | -0.0133 |
| drop-top-5 | -0.0066 |
| **drop-top-10** | **+0.0050 — the gain REVERSES** |
| drop-top-20 | +0.0243 |
| **top-10 share of the total gain** | **1.163** |

Five largest gains: 1TOR -0.220, 1N9U -0.210, 2BP4 -0.183, 2LWU -0.171, 2LWS -0.158.
Five largest losses: 1I8E +0.057, 1RG4 +0.060, 1D6X +0.076, 1U62 +0.089, 1KWE +0.146.

**A top-10 share above 1.0 means the ten best targets carry MORE than the entire net gain —
the other forty are net +0.20 A.** With a median of only -0.009, the near-even win/loss is
**many small losses plus a few large wins**, which is precisely the concentrated pattern that
killed `leg_torsion` in C12. Per fold the gain is same-signed but spans a 25x range
(fold 0 -0.009, fold 1 -0.029, fold 2 -0.002, fold 3 -0.054, fold 4 -0.047).

**The methodological lesson, and it is mine as much as theirs: the win/loss ratio was never
the discriminating statistic. The drop-top curve is.** I repeated in a published report a
reassurance that had not been measured, having spent this entire sprint insisting that every
result carry a concentration check. That the claim came from a colleague is not a defence —
it was published under my name and the check takes seconds.

**Their diagnosis of the failure is sharper than mine, and it changes the fix.** The check was
not skipped. `I.paired` returns `drop_top10_mean_diff` and `top10_share` **in the same dict**
the W/L was read from — the answer was already on screen. *"I looked at one field of a
structure that contained the answer and reported the reassuring one."* So the failure mode is
not forgetting to run the check, it is **reading past it**, and the fix is structural rather
than behavioural: the check now emits the three concentration statistics as a single PASS/FAIL
verdict block with the sample-size caveat inline, so the n=126 result cannot be read without
them.

### The free early warning, which generalises beyond this result

Their second point is the most transferable thing to come out of the episode, and it is
available at zero cost in every paired comparison this project runs:

> **A near-even W/L *together with* a CI excluding zero is not weak evidence of a consistent
> effect — it is positive evidence of a CONCENTRATED one.**

If an effect were uniform, a -0.025 A mean would produce a **lopsided** W/L and a **median near
the mean**. Observing 26W/24L with a median of **-0.009** against a mean of **-0.025** already
says the mean is being set by the tail — *before computing a single drop-top number*.

**The median-versus-mean gap is the free early warning**, and it is now printed on the same
line as the W/L. Ratio here: the mean is 2.8x the median. Any future result where that ratio
is large should be treated as concentrated until the drop-top curve says otherwise.

**Not filed REFUTED yet, for one specific and legitimate reason.** At n=50, drop-top-10
removes **20%** of the sample, which is a harsh test the rule was not written for. At n=126 it
removes 8%. The 126-target k=30 run decides it: if drop-top-10 and drop-top-20 hold their sign
there it is real; if the curve decays as it does here, it is two or three targets and it will
be filed REFUTED alongside `leg_torsion`.

### C10-UN-RETRACTION. My concentration verdict was itself wrong, and the error is in the RULE

**The k=30 result is restored.** The "FAILS the concentration check" verdict above is
withdrawn, and the reason matters more than the result: **a raw drop-top threshold is not a
valid concentration test.**

When an effect's mean is small relative to its per-target spread, discarding the ten most
favourable targets removes a large share of the total **even if every target carries an
identical effect**. So at low signal-to-noise the test *must* fail a uniform effect. I mandated
that rule in the sprint brief, enforced it against `leg_torsion`, and then applied it here
without ever calibrating it. Both of us read a guaranteed failure as evidence of concentration.

**The correct test compares the observed statistic to its own null — a UNIFORM effect of the
same mean and sd.** k=30, n=126, mean -0.0221, per-target sd 0.0789 (**mean/sd = 0.28**),
4,000 simulations:

| statistic | observed | null: a UNIFORM effect | percentile |
|---|---|---|---|
| drop-top-10 | -0.0080 | -0.0094 [-0.0234, +0.0048] | **0.62** |
| drop-top-20 | +0.0028 | +0.0006 [-0.0136, +0.0150] | **0.62** |
| top-10 share | 0.668 | 0.698 [0.389, 1.543] | **0.62** |

**Every statistic sits at the 62nd percentile of its own null. There is no evidence of
concentration in any arm** — including the n=50 case I retracted on, where the observed share
of 1.163 is *less* concentrated than the null mean of **1.410**.

**The symmetric caution, which must carry equal weight.** "Not concentrated" does **not** mean
"demonstrated uniform". At mean/sd = 0.28 the drop-top test has almost no power and cannot
settle it either way. What actually carries the result is the CI and **the per-fold table**:
all five pinned folds same sign, **-0.014 to -0.036, a 2.6x range** — down from 25x at n=50,
which is noise narrowing exactly as it should. FAIL18 moves with the rest (-0.011 against
-0.024).

### The final numbers, n=126, full development instrument

| arm | RMSD | vs A | CI | geom dev | x ff14SB equilibrium | clashed | Ramachandran |
|---|---|---|---|---|---|---|---|
| A: projection (incumbent) | 3.205 | — | — | 0 by construction | — | 0.000 | 0.734 |
| B: raw average (invalid) | 3.050 | -0.155 | [-0.191,-0.122] | 0.2480 | 14.6x | **0.452** | 0.836 |
| **B + AMBER k=30** | **3.183** | **-0.022** | **[-0.036,-0.009]** | 0.0192 | **1.13x** | 0.000 | **0.874** |
| k=100 (n=50) | 3.101 | -0.070 | [-0.098,-0.043] | 0.0374 | 2.20x | 0.000 | 0.894 |
| k=300 (n=50) | 3.066 | -0.105 | [-0.140,-0.073] | 0.0665 | 3.91x | 0.000 | 0.857 |

**So "never a better number" is retracted and now stands retracted at VALID geometry** — 1.13x
the force field's own equilibrium — rather than only at strained geometry.

### Caveats, and the third had to be strengthened twice

1. **`tuning126` only.** It has never been near the benchmark, where **nothing has ever
   transferred** in fourteen sprints.
2. **It is 0.7% of the baseline.** A 0.022 A shift on a 3.205 A number.
3. **THE RESTRAINT CONSTANT WAS SELECTED ON DEV-SET RMSD, and there is no native-free rule in
   hand that picks it.** I first wrote "chosen with the frontier visible", implying a
   principled rule merely applied at the wrong time. **The rule is degenerate and the agent
   retracted it themselves:**

| k | strain | x equilibrium | distance from 1.0x | RMSD |
|---|---|---|---|---|
| 0 | 0.0170 | **1.00x (circular)** | **0.00** | 3.477 |
| 2 | 0.0152 | 0.89x | 0.11 | 3.248 |
| 10 | 0.0134 | 0.79x | 0.21 | 3.207 |
| **30** | 0.0192 | 1.13x | 0.13 | **3.183** |
| 100 | 0.0374 | 2.20x | 1.20 | 3.101 |
| 300 | 0.0665 | 3.91x | 2.91 | 3.066 |

   **k=0 is at 1.00x by definition — it IS the equilibrium** — so "closest to equilibrium"
   trivially selects free relaxation, the **worst** arm at 3.477. The non-circular variant,
   "tightest restraint whose strain does not exceed equilibrium", selects **k=10, which gives
   +0.001 A — no gain at all.** Selecting k=30 requires a threshold of ~1.15x, and 1.15 was
   chosen with the frontier visible.

   **What the yardstick does and does not do.** It legitimately *rules out* k=100 and k=300 as
   strained, and that is what stopped -0.105 A becoming the headline. It does **not** select
   k=30 out of {k=2, k=10, k=30} — all three sit at or near equilibrium, and k=30 was picked on
   development-set RMSD from within the acceptable set. **A confirmation must PRE-REGISTER the
   restraint constant**, and a pre-registered k could land anywhere from **+0.001 (k=10) to
   -0.022 (k=30)**.

### Two further checks, one in my favour and one against a claim of mine

**The baseline choice is conservative, not favourable.** `I.project` returns two arms:
`fit_ca` (lambda=0) at 3.2052 and `ca` (lambda=0.3, Ramachandran-penalised) at 3.2126. Arm A
uses `fit_ca` — the **better** of the two. Against the lambda=0.3 arm the gain would be
**larger**: -0.0295 A [-0.0465, -0.0122], 70W/56L. So -0.022 A is the conservative figure.

**"The first native-free accuracy improvement this project has measured" is FALSE, and it was
mine.** The agent declined to carry the superlative, saying they could not verify a claim
spanning fourteen sprints from their own work and asking that it be attributed or checked. It
does not survive the check: **Sprint 12's score-filter + consensus medoid measured
-0.172 A [-0.316, -0.027]** and is recorded in memory as "the first CI excluding zero". That
is both **earlier and roughly eight times larger** than this result.

The correct statement is narrower: *this is a small native-free improvement on the development
instrument, at valid geometry, with a restraint constant selected on that same instrument.* It
is not a first, and I should have checked my own memory index before writing a superlative into
a published report.

### The methodological lesson, now two-sided and built into the tooling

I read a W/L as reassurance (wrong — the median-versus-mean gap had already flagged a
tail-driven mean), then read an uncalibrated drop-top curve as disqualifying (**also wrong** —
that is what uniform looks like at this signal-to-noise).

> **The concentration check is NECESSARY AND NOT SUFFICIENT. It needs its own null, and
> mean/sd must be printed beside it so a reader can see when it has no power.**

The mandated check caught `leg_torsion` correctly and misfired here, and the difference is
signal-to-noise. Both fixes are now in the checker rather than left to judgement.

**My statement "identical accuracy, never a better number" is retracted.** The correct
statement is: *at the force field's own equilibrium geometry, relaxation is worth
-0.025 A [-0.046, -0.005] over the projection at n=50, pending confirmation at n=126; tighter
restraints buy more RMSD but at measurable strain.*

Worth recording for method: their reconstruction asserted agreement with my arm B at run time
rather than assuming it, **and that assertion caught their own first version's frame bug**,
which had produced a spurious "86% bond contraction". A runtime cross-check against another
agent's number is cheap and it worked.

## C11. The ablation ladder closes, and it bounds the VQE contribution to 0.171 A — already obtained classically

**DEMONSTRATED.** `s14/results/vqe_classical_limit.json` plus C9 and C2.

A variational state prepared to exactly reproduce the leakage-safe prior, with **no
optimisation at all**, is the same object as sampling that prior. Measure what it emits when
its samples are aggregated in coordinate space and you have the honest control the brief's
section 23 demands — the number any quantum claim must be differenced against, rather than
against the single torsion vector at 4.072 A.

| B (measurements) | mean | median | <2 A | FAIL18 | vs incumbent |
|---|---|---|---|---|---|
| 1 | 4.192 | 4.138 | 0.13 | 6.154 | +0.988 [+0.765,+1.229] |
| 5 | 3.662 | 3.339 | 0.19 | 5.822 | +0.457 [+0.269,+0.659] |
| 20 | 3.533 | 3.283 | 0.22 | 5.735 | +0.329 [+0.159,+0.506] |
| 75 | 3.515 | 3.274 | 0.24 | 5.799 | +0.311 [+0.139,+0.493] |
| 200 | **3.485** | 3.299 | 0.24 | **5.771** | +0.281 [+0.113,+0.461] |

Set beside C9 and C2, the sprint's ablation ladder is complete on one instrument, one
metric, one candidate space:

| stage | emitted | increment | what it buys |
|---|---|---|---|
| A torsion information, one committed vector | 4.072 | — | the best torsion channel in the project |
| B + sample it 200 times, coordinate-average | 3.485 | **-0.587** | error cancellation, no objective involved |
| C + rank 4,000 with the structural Hamiltonian, average top 200 | 3.314 | **-0.171** | **the entire value of the objective** |
| incumbent retrieval pipeline | 3.204 | -0.110 | still ahead |

**The whole contribution an objective makes in this architecture is 0.171 A, and coordinate
aggregation makes 0.587 A** — three and a half times as much, from an operator with no
information in it whatsoever.

**This bounds the quantum question sharply and quantitatively.** A VQE's only possible
contribution here is to prepare a *better distribution* than the prior — that is, to do the
job that objective-driven selection does. That job is worth **0.171 A**, and it has already
been done classically by ranking 4,000 samples. C8 shows the same search saturates by
evaluation 300 out of 20,000. So the quantum ceiling in this representation is not merely
small, it is **a fraction of 0.171 A**, and it must be won against a classical control that
already banked it.

Stated the way the brief asks: *the structural prior provides 0.587 A through aggregation
alone, objective-driven selection provides a further 0.171 A, and there is no measured
headroom left for VQE or CVaR to provide anything beyond that.*

**A recurring pattern, now on its third independent appearance.** At B = 200 the classical
limit's FAIL18 mean is **5.771 against the incumbent's 6.026**. Zero-information and
low-information arms keep beating the real system on the 18 hardest targets — the constant
helix (5.887), the class prior argmax (5.670), and now prior sampling with consensus. On the
failure class, this system's information is worse than none, exactly as Sprint 12 found by a
different route.

## C12. The sprint's strongest open lead does not replicate — `leg_torsion` REFUTED at n=126

**REFUTED.** `s14/hamil3.py`, 126 targets, 3,000 native-free configurations each, both a
mixed and a strictly uniform proposal, all four required axes reported.

The ENER agent's component sweep left exactly one survivor as a HYPOTHESIS: the Legacy
`torsion` term selected **2.954 A**, -0.835 [-1.607,-0.105], 6W/3L, clearing the constant-helix
baseline and surviving the helix controls. If it had held it would have been the best
native-free selection result in the project. They filed it correctly as a hypothesis, because
it was measured on the **nine enumerated targets only and two of those nine carried 68% of
the effect**.

At fourteen times the sample it is not merely weaker. It is the **worst arm tested**:

| arm | rho all | rho decile | argmin | vs incumbent |
|---|---|---|---|---|
| prior + distogram | 0.388 | **0.367** | **3.571** | +0.367 [+0.251,+0.487] |
| distogram only | 0.539 | 0.095 | 3.633 | +0.429 [+0.283,+0.580] |
| all three terms | 0.376 | 0.323 | 3.670 | +0.466 [+0.333,+0.615] |
| prior only | 0.304 | 0.356 | 3.836 | +0.632 [+0.388,+0.885] |
| `leg_torsion` + distogram | 0.331 | 0.068 | 3.910 | +0.706 [+0.516,+0.913] |
| `leg_torsion` + prior | 0.269 | 0.140 | 4.306 | +1.102 [+0.800,+1.415] |
| **`leg_torsion` only** | 0.230 | **-0.073** | **4.933** | +1.729 [+1.419,+2.052] |

Under a strictly uniform proposal it is worse still: global rho **-0.073**, argmin 5.195,
+1.991 against the incumbent. Its in-decile correlation is **negative under both proposals**
— it anti-ranks where a search lives.

**And it poisons every combination it enters.** Adding it to prior+distogram moves the argmin
from 3.571 to 3.670. Adding it to the distogram alone moves 3.633 to 3.910. There is no
weight at which it contributes.

The 2.954 A figure was concentration on nine targets, and the concentration analysis called
it in advance. This is the strongest possible argument for the sprint rule that every result
carries a drop-top-10 check and a per-target breakdown: **the lead looked like the best
finding of the sprint and was an artefact of two targets.**

One thing does survive from that arm, and it is worth keeping: `prior + distogram` at
in-decile **+0.367** is the best-ordered native-free objective measured, and it beats the
distogram alone (+0.095) by a wide margin on that axis while being slightly worse on global
rho. The two terms are genuinely complementary in the region a search occupies. That is C7's
finding, confirmed on an independent sample with a third term available and rejected.

## C13. CERTIFIED: the structural objective's global optimum is 0.885 A BETTER than random — the first objective in this project of which that is true

**DEMONSTRATED, and certified rather than sampled.** Measured by the VQE agent over the
**full 262,144-configuration enumeration on all nine enumerated targets**, uniform
population, both terms standardised over the complete space. No sampling anywhere.

Reference points on those nine: space best 0.969 A, random draw 3.781 A, Legacy's certified
argmin 3.920 A, the Sprint 13 torsion prior's certified argmin 3.636 A.

| w (distogram share) | global rho | in-decile | **RMSD at CERTIFIED argmin** | decile mean | decile best |
|---|---|---|---|---|---|
| 0.00 | +0.084 | +0.087 | **2.866** | 3.531 | 1.020 |
| 0.25 | +0.176 | +0.105 | 2.896 | 3.410 | 1.020 |
| 0.50 | +0.272 | +0.113 | 3.034 | 3.318 | 1.121 |
| 0.75 | +0.336 | +0.112 | 3.073 | 3.253 | 1.308 |
| 1.00 | +0.360 | +0.115 | 3.237 | 3.244 | 1.513 |

**This is the sprint's clearest positive result, and it corrects the framing of C7 and C8.**
Legacy's certified optimum is **0.139 A WORSE** than a random draw. The Sprint 13 torsion
prior's is 0.145 A better. The structural Hamiltonian's is **0.885 A better** — roughly six
times the torsion prior's margin, and the only objective measured anywhere in this project
whose optimum is somewhere you would actually want to go.

So **Sprint 13's "the certified optimum is in the wrong place" does NOT hold for this
objective.** That claim was specific to the physical energies, and building a structural
rather than energetic objective repaired it. This survives certification, which is the
strongest form the claim can take — it cannot be attributed to insufficient search, because
the optimum is exact over the complete space.

**What did not get repaired, and it is the whole story.** The certified optimum still sits
**1.876 A worse than the best structure inside the objective's own lowest decile** (2.896
against 1.020). That gap is irreducible by definition: no optimiser can do better than the
certified optimum. And it agrees closely with my sampled selection gap of 2.040 A at budget
20,000, so **two independent routes — sampling on 126 targets and exact enumeration on nine —
measure the same ~1.9 A of discrimination loss.**

The refined statement, which is more accurate and more useful than what C8 said alone:

> The objective now points in a genuinely useful direction. It still cannot tell the good
> structures apart from the mediocre ones once it gets there. Optimisation was never the
> binding constraint; discrimination is, and it survives infinite budget.

**AN INVERSION THAT CHANGES THE ARCHITECTURE CHOICE.** Global rho and certified-argmin
quality move in **opposite directions** along w. Adding distogram weight improves global
ordering monotonically (+0.084 to +0.360) while making the certified optimum monotonically
worse (2.866 to 3.237) and the best decile member worse (1.020 to 1.513) — yet the decile
MEAN improves (3.531 to 3.244).

Mechanism: **the distogram orders the bulk; the retrieval torsion prior places the optimum.**

The consequence is not academic. If the terminal operator is an argmin or a small-m
aggregation, w = 0 to 0.25 is correct and more distogram actively hurts. If it is a large-m
set mean, higher w is correct. The Sprint 12 terminal-operator law says the output consumes
the set MEAN — so the two axes give opposite answers and **w must be chosen jointly with the
terminal operator, never on global rho alone.** My choice of w = 0.25 was defensible for the
argmin readout I was using, but I selected it on in-decile rho, which is the wrong criterion
and non-monotone besides.

**HETEROGENEITY, and it is large enough to matter.** Per-target global rho at w = 0.25:
7N2I +0.469, 2P5H +0.366, 2MK7 +0.343, 9UV5 +0.322, 6EY3 +0.306, 6S0N +0.140, 8IS3 +0.048,
6F3V **-0.150**, 1CS9 **-0.262**. **Two of nine targets anti-rank.** Certified argmin ranges
from 0.675 A (2MK7) to 5.256 A (1CS9, worse than that target's own random draw of 4.064).
The 2.896 mean is a mixture of a few excellent targets and two failures, not a typical
target, and it must be reported with that spread.

**In-decile cross-check.** They measure +0.105 at w=0.25 on the uniform full enumeration,
against my +0.046 under uniform sampling and the withdrawn +0.370 under the prior-drawn
proposal. Their figure sits between the two and much nearer mine-under-uniform, which
independently confirms that the +0.370 was a proposal artefact rather than a property of the
objective.

### C13-CORRECTION. I read a point statement as a set statement. The distogram STAYS.

On the strength of "identical certified argmin at w=0 and w=0.25 on 7 of 9 targets" I
proposed dropping the distogram term entirely for any argmin or small-m readout, as a simpler
architecture. **That inference does not follow and the VQE agent refuted it before it was
locked in.** Comparing the top-m SETS rather than only their argmin:

| m | w=0 best | w=0.25 best | w=0 mean | w=0.25 mean | Jaccard(w0, w0.25) |
|---|---|---|---|---|---|
| 1 | 2.866 | 2.896 | 2.866 | 2.896 | 0.667 |
| 5 | 2.776 | **2.527** | 3.048 | 2.988 | 0.496 |
| 20 | 2.441 | **2.349** | 3.045 | 2.991 | 0.533 |
| 75 | 2.138 | 2.127 | 3.068 | 3.003 | 0.580 |
| 200 | 1.956 | 1.908 | 3.108 | 3.056 | 0.612 |

**Half the set turns over.** Jaccard overlap is 0.50-0.61 for m >= 5 even where the single
argmin coincides — a term can leave the minimum untouched while reshaping everything around
it, and the set is what a real terminal operator consumes. And w=0.25 is **better at every
m >= 5 on both readouts**. The only m favouring w=0 is a literal m=1, by +0.030 A, inside the
1.389 A between-target sd.

Two further arguments I should have raised myself and did not:

**The length extrapolation is the wrong way round.** All nine enumerated targets are n=9 — 36
pairs, maximum separation 8. At n=16 the distogram has 120 pairs and separations to 15, which
by the locality theorem is exactly where its many-body content lives. **The enumerated nine
are the shortest bin in the instrument, and no term whose content scales with sequence
separation should be judged on them alone.** That caution generalises beyond this decision.

**Dropping it raises per-target variance by 46%** (between-target sd 0.951 to 1.389). Against
a project whose defining failure mode is a structural class holding 10 of 18 catastrophic
targets, trading consistency for a 0.030 A mean gain is the wrong direction.

**ADOPTED: keep the distogram at w in [0.25, 0.5].** The corrected wording is theirs: *the
distogram does not move the argmin on 7/9 targets, but it turns over half the top-m set and
improves it at every m >= 5.*

## C15. The discrimination floor is UNIVERSAL — the structural objective does not break it either

**DEMONSTRATED.** `s14/discrim.py`, full 262,144-configuration enumeration, no sampling and
no proposal, using the ENER workstream's own `pair_accuracy` so the numbers compose with
theirs. A state-ordering assertion rebuilds the snap configuration's index and checks it
against the cache, because an off-by-one in digit order would pair every objective value with
the wrong RMSD and produce a beautiful, meaningless result.

Pairwise ranking accuracy against quality gap (first two targets; full run in flight):

| objective | 0.00 | 0.25 | 1.00 | 2.00 | 3.00 | reaches 55% at | argmin |
|---|---|---|---|---|---|---|---|
| Legacy | 0.479 | 0.475 | 0.468 | 0.472 | 0.489 | **never** | 4.552 |
| 1-local prior | 0.514 | 0.515 | 0.523 | 0.540 | 0.564 | 3.00 A | 4.685 |
| hamil w=0 | 0.510 | 0.509 | 0.510 | 0.505 | 0.502 | **never** | 3.020 |
| **hamil w=0.25** | 0.515 | **0.515** | 0.516 | 0.515 | 0.519 | **never** | **2.966** |
| hamil w=1 | 0.526 | 0.528 | 0.532 | 0.543 | 0.562 | 3.00 A | 4.426 |

**RETRACTED ON THE FULL NINE TARGETS. The opposite is true.** The two-target reading above was
the worst possible sample — 1CS9 and 2MK7 are precisely the two outliers the VQE agent had
already flagged (1CS9 anti-ranks at rho -0.262; 2MK7 has the best argmin of the nine). Over
all nine enumerated targets:

| objective | 0.00 | 0.10 | 0.25 | 0.50 | 1.00 | 2.00 | 3.00 | 55% at | argmin |
|---|---|---|---|---|---|---|---|---|---|
| Legacy | 0.525 | 0.528 | 0.533 | 0.542 | 0.566 | 0.637 | 0.721 | 1.00 A | 3.920 |
| 1-local prior | 0.504 | 0.504 | 0.505 | 0.506 | 0.513 | 0.532 | 0.555 | 3.00 A | 3.636 |
| hamil w=0 | 0.528 | 0.530 | 0.532 | 0.536 | 0.546 | 0.564 | 0.581 | 1.50 A | **2.878** |
| hamil w=0.25 | 0.560 | 0.563 | **0.569** | 0.579 | 0.603 | 0.656 | 0.707 | **0.00 A** | 2.896 |
| **hamil w=1 (pure distogram)** | **0.633** | 0.641 | **0.654** | 0.675 | 0.716 | 0.782 | 0.816 | **0.00 A** | 3.237 |

**The distogram BREAKS the floor decisively.** ENER's result was that no *physical* objective
exceeds 0.511 pairwise accuracy below a 0.25 A quality gap. The pure distogram reaches
**0.654** there, and **0.633 at zero gap** — it discriminates between structures of
essentially equal quality. The mixed objective reaches 0.569. So the floor is a property of
the **physical energies**, exactly as ENER scoped it, and a learned native-free structural
objective is not bound by it.

**But breaking it does not produce a better optimum, and that is the finding.** The pure
distogram has the **best bulk discrimination and the WORST argmin of the hamil family**
(3.237 against w=0's 2.878). The 1-local prior is its mirror image: the worst discrimination
(0.505 at a 0.25 A gap, indistinguishable from chance) and among the best argmins.

**So neither in-decile rank correlation nor bulk pairwise accuracy predicts selection
quality.** Both were treated during this sprint as proxies for objective usefulness and
**neither is**. In-decile rho is non-monotone in quality (the VQE agent's correction); bulk
pairwise accuracy is measured over random pairs and is dominated by the 99% of the space a
search never visits. The quantity that governs an argmin is accuracy *within the extreme
tail*, and nothing measured this sprint reports it directly.

This is the third independent confirmation of the same mechanism, now from a completely
different statistic: **the distogram orders the bulk, the torsion prior places the optimum.**
It also retrospectively justifies keeping both terms, and is a fourth axis favouring the
mixture the VQE agent argued for.

**And the dissociation is the cleanest statement of what was achieved.** hamil w=0.25 has by
far the best argmin (2.966 against Legacy's 4.552) and among the *worst* pairwise ranking
accuracy. **Its entire advantage is in WHERE its optimum sits, not in any ability to tell
structures apart.** That reconciles the two headline results — a well-placed certified optimum
and a flat discrimination profile are compatible, and this sprint produced both.

*Sample-size note, so the dossier never carries two numbers for one statistic: the figures in
this table are the first two enumerated targets. Over all nine, Legacy's argmin is 3.920 and
hamil w=0.25's is 2.896. Same quantities, different subsets.*

### C15b. A PRE-REGISTERED prediction, recorded before the data exists

The VQE agent read C15 as a mechanism and derived a falsifiable consequence from it. Recording
it here **before the closing experiment runs**, so it is a test rather than a story fitted
afterwards.

**The mechanism.** A best-in-project argmin combined with chance-level pairwise accuracy is
only jointly possible if the objective's entire skill is concentrated in a vanishingly small
fraction of its range: **it knows where the very bottom is and is blind everywhere else.**
Skill lives in the extreme tail and nowhere else.

**Why that matters for CVaR specifically.** From the Part A derivation, `CVaR_alpha` depends
on the objective only through its lowest alpha-fraction, and at `alpha <= p(x*)` only through
the argmin itself. So an objective whose skill is tail-concentrated should be **the first case
in this sprint where small-alpha CVaR has a structural advantage over the expectation value** —
alpha = 1 averages that skill away against 99.5% noise, and small alpha does not.

**The prediction, on hamil w=0.25 at matched budget:**

1. **CVaR at alpha = 0.05-0.1 returns a BETTER structure than alpha = 1.0**, reversing the
   ordering measured on the signal-0.3 blend family (where alpha = 1 won, 2.266 against
   2.459).
2. The margin is **small, a few tenths of an angstrom at most**, because the tail skill is
   real but the floor beneath it is still 0.515.
3. It will **not** beat greedy 1-opt on structure, because greedy also reaches the tail and
   the discrimination floor binds both.

If (1) comes out the other way, the tail-concentration mechanism is wrong and will be recorded
as refuted. If (1) holds, it is **the only place in this sprint where CVaR's defining feature
does something a mean-based classical method cannot** — a genuinely quantum-motivated positive
rather than a null, though (3) would still mean VQE is not the right tool for the job.

## C14. The last route to 2.0 A is closed by arithmetic, not by accuracy

**DEMONSTRATED.** SHIFT workstream, `s14/shift_FINDINGS.md`, eight sections.

Chemical shifts were the only remaining channel capable in principle of supplying the
sigma ~15 deg that 2.0 A requires. The route is now closed, and the closing argument does not
depend on how good a shift-based predictor could be made.

**Availability, unioning four sources including BMRB's own BLAST over deposited polymer
sequences — which no prior sprint used and which alone supplied 15 of the final 54:**

| tier | n of 126 |
|---|---|
| no deposition retrieved | 28 (7 are PDBj-BMRB 36xxx, unresolved rather than absent) |
| **1H-only, TALOS-N cannot run** | **42 (33%)** |
| heteronuclear, gate < 75% | 2 |
| heteronuclear, gate 75-90% | 3 |
| heteronuclear, gate >= 90% | 51 |

**54 of 126 targets runnable (42.9%), upper bound 61.** All 98 sequence matches were exact
substrings, so there is no homolog contamination. The dominant failure is **not missing
depositions** but **1H-only homonuclear peptide NMR on unlabelled synthetic material**. Only
7 of 126 carry the full six-nucleus TALOS-N input.

**THE ARITHMETIC THAT CLOSES IT.** With availability `a = 54/126` and the incumbent at 3.2871
on the other 72 targets, **ORACLE-perfect torsions on every target where shifts exist leave
the full instrument at 2.021 A.** Reaching 2.0 A needs `a > 0.436`, i.e. **55 targets. We
have 54.** The route fails by one target even with a perfect predictor. The realistic
TALOS-N ceiling on the mixed arm is **2.884 A**.

**MY GUIDANCE WAS APPLIED AND THE GATE IT ADDRESSED IS REFUTED.** I sent them the positional
cost law and asked for position-weighted coverage. They did it, and the answer is that
**coverage is not the binding constraint at all**: where the channel runs, coverage is 0.952
raw and 0.969 position-weighted, with only a -0.040 terminal deficit, and measured gaps cost
**0.014 A at sigma 12** — inside seed noise. So Sprint 12's "coverage is the make-or-break
parameter" and the memory note's "close below 75% coverage" are **refuted as the gate**. The
binding quantity is the fraction of *targets* with heteronuclear shifts (0.444), not the
fraction of *residues* within them (0.952). The pre-registered test passed easily and was
simply the wrong test.

**THE 1.486 A FIGURE IS SUPERSEDED.** Under TALOS-N's own published error mixture (87.5%
Strong, 3.9% Generous, 8.6% Ambiguous; 2.8% and 21% gross error rates) the channel lands at
**2.347 A deployable**, 2.029 A with oracle basins. The 1.486 A number describes a channel
with no bad predictions and no ambiguity, which is not a channel that exists.

**Two surprises worth keeping.** Availability on the failure class is *higher*, not lower
(fibril/lasso 0.500, FAIL18 0.444, other 0.398) — solid-state NMR of fibrils requires 13C/15N
labelling. And missing nuclei are **not** a handicap: five nuclei without C' give 48.9 deg,
CA alone gives 55.7 deg, against 50.9 deg for the full six, all within ~5 deg. Sequence-only
through the same architecture is 77.5 deg.

**What the real channel emits.** On the matched 54 targets the incumbent is **3.1133 A**, not
3.213 — the matched-subset comparison the brief required. E-SHIFT emits 3.832 and loses at
**+0.719 [+0.255,+1.194]**, 18W/36L. But it beats every one of its own nulls decisively:
permuted shifts -1.208 [-1.739,-0.648], zero-secondary-shift -0.859 [-1.338,-0.377],
predicted-shifts -0.743 [-1.297,-0.182]. **The experimental channel is real and 2.7x larger
than the entire sequence-to-torsion channel** (-0.86 A against Sprint 13's -0.32 A). The
whole shortfall is the gross-error rate: 23.5% against TALOS-N's 3.5%. Conditional on not
being grossly wrong, E-SHIFT sits at **sigma 19.1 deg**, already in the 2.0-2.2 A regime.

Their measured error coherence is lag-1 **0.079**, near-i.i.d., so my coherence surface
applies unchanged and the 15.1 deg i.i.d. requirement is the right one.

**Torsions from PREDICTED shifts are closed as pre-registered**: +0.116 A from a
zero-secondary-shift null and **1.0 A worse than using the sequence directly**. A clean
negative that closes a route which looked attractive on paper.

### C14b. The best justification for a qubit this project has produced

Two findings from that workstream are transferable and both point the same way.

**Do not collapse the posterior.** ORACLE search over E-SHIFT's own top-8 support, against
its own argmax: **-2.253 A [-2.642,-1.865], 53W/1L**, at 39 qubits. Independently corroborated
on the simulated channel, where collapsing costs +0.30 A over dropping and carrying both
basins is worth +0.62 A. **A bimodal shift posterior is one qubit with a physical
justification** — the first genuine physical argument for a specific qubit anywhere in this
project. Every previous qubit was a discretisation choice.

**Confidently wrong costs 2-3x what absent costs** (at 5% error: 2.459 against 1.890; at 10%:
3.130 against 2.114). So a torsion prior should **abstain rather than guess**, and a
variational search should consume abstentions as free registers. That is a concrete design
rule, and it is the opposite of what a point-estimate predictor does.

**The honest framing, which they state and I endorse.** Every E-SHIFT number is
**NMR-restrained structure determination reported in its own column** — the deposited
coordinates were solved using these very shifts plus NOEs. It is never an improvement to the
sequence-only ladder and must never be compared to it head-to-head without that label.

**Their unprompted suggestion is the best remaining use of the data**, and it inherits none
of the negatives above: 54 targets now have an experimental observable, and a *forward* shift
predictor scores candidate structures against measured shifts **with no native coordinates at
all**. That is a discrimination channel rather than a generation channel — which is precisely
the bottleneck C13 identified.

## C16. THE SELECTION GAP HAS A CLOSED-FORM ACCOUNT — an objective picks a good neighbourhood and then draws at random from it

**DEMONSTRATED.** VQE workstream, `s14/vqe_tailacc.py`, full enumeration on all nine targets.
This is the sprint's central scientific result and it resolves everything else.

**The statistic.** Restrict to the objective's own lowest `frac` of the 262,144
configurations, then measure ordering accuracy among pairs separated by at least 0.25 A in
true CA-RMSD — the ENER threshold. Objective ties are scored 0.5, because scoring them as a
win or a loss would import index order (the recorded tie-breaking trap).

| frac of space | n configs | w=0 | w=0.25 | w=0.5 | w=1.0 |
|---|---|---|---|---|---|
| 1.0 (bulk) | 262,144 | 0.532 | 0.570 | 0.611 | **0.655** |
| 0.1 | 26,214 | 0.534 | 0.542 | 0.545 | 0.549 |
| 0.01 | 2,621 | 0.526 | 0.528 | 0.521 | 0.509 |
| **0.001 (tail)** | 262 | 0.512 | 0.521 | 0.520 | **0.390** |
| **tail − bulk** | | **-0.020** | **-0.049** | **-0.091** | **-0.266** |
| certified argmin | | 2.866 | 2.896 | 3.034 | 3.237 |

**`tail − bulk` is monotone in w and tracks the certified argmin exactly.** The pure distogram
is the best bulk ranker in the family at 0.655 and **orders structures BACKWARDS inside its
own top 0.1% at 0.390**. That is the sharpest form of "the distogram orders the bulk, the
torsion prior places the optimum" produced anywhere this sprint, and it is the fourth
independent confirmation of that mechanism on a statistic unrelated to the first three.

### The mechanism, stated exactly

If an objective is at chance inside its own tail, its argmin should be statistically
indistinguishable from a **uniform draw from that tail**. Tested directly:

| w | tail m | actual argmin | tail MEAN | tail BEST | percentile of the argmin within its own tail |
|---|---|---|---|---|---|
| 0.00 | 26 | 2.866 | 3.000 | 2.335 | 0.419 |
| 0.00 | 262 | 2.866 | 3.113 | 1.889 | 0.434 |
| 0.00 | 2,621 | 2.866 | 3.299 | 1.365 | 0.409 |
| 0.25 | 262 | 2.896 | 3.066 | 1.889 | 0.480 |
| 1.00 | 262 | 3.237 | 3.103 | 2.515 | **0.558** |
| 1.00 | 2,621 | 3.237 | 3.074 | 2.007 | **0.609** |

**The argmin sits at the 41st to 48th percentile of its own tail — indistinguishable from the
50th.** At w=1 it sits at the 56th to 61st, i.e. *worse* than a random draw from its own tail,
exactly as the 0.390 anti-ranking predicts. The actual argmin equals the tail MEAN to within
0.13-0.43 A while the tail BEST is 1.0-1.9 A better.

> **An objective's argmin quality is set by WHERE its tail sits, not by any ability to order
> within it. It picks a good neighbourhood and then draws essentially at random from it.**

w=0 wins the argmin because its extreme tail has the best *mean*, not because it ranks better
inside. That is why the certified optimum improved 0.885 A over random while discrimination
did not improve at all: **those are two different faculties, and only the first one was ever
repaired.**

### And this closes the selection gap quantitatively

`tail mean − tail best` at w=0, m=2,621 is 3.299 − 1.365 = **1.93 A**.

| measurement | value |
|---|---|
| tail mean minus tail best (this account) | **1.93 A** |
| certified selection gap, 9 enumerated targets | **1.876 A** |
| sampled selection gap, 126 targets, budget 20,000 | **2.040 A** |

**These are the same quantity.** The selection gap *is* the price of drawing at random from a
tail you cannot order. Three measurements on two target sets by two workstreams, now with a
closed-form account rather than an empirical value.

### C15b is REFUTED, and for a better reason than the one I supplied

The pre-registered prediction was that small-alpha CVaR would win on a tail-concentrated
objective. **No objective in this family is tail-concentrated.** w=0's tail accuracy is 0.512
— chance, and slightly *below* its own bulk 0.532. Every objective tested is at chance or
worse inside its own tail; the skill is in the bulk for all of them. The substrate the
prediction required does not exist.

The revised prediction, recorded in its place: **small-alpha CVaR should be neutral to harmful
everywhere**, because it concentrates attention on exactly the region where every objective
carries no ordering information. That is consistent with the blend-family sweep, where
alpha = 1 beat alpha = 0.05 (2.266 against 2.459).

Note the refutation was **not** caused by my two-target artefact. The nine-target version kills
it harder and for a deeper reason: the premise itself was false, not merely mismeasured. The
agent registered the prediction, derived its own refutation, and reported it unprompted.

### Methodological recommendation, now settled

Of the four proxies this sprint used for objective usefulness, **only tail-restricted pairwise
accuracy predicts argmin quality**:

| proxy | predicts the argmin? |
|---|---|
| global rank correlation | **No** — w=1 has the best global ordering and the worst argmin |
| in-decile rank correlation | **No** — non-monotone in objective quality |
| bulk pairwise accuracy | **No** — w=1 wins it decisively and has the worst argmin |
| **tail-restricted pairwise accuracy** | **Yes** — `tail − bulk` tracks the certified argmin exactly |

It is cheap (`python -m s14.vqe_tailacc`, seconds on cached tabulations). Caveat: at
frac = 1e-4 (26 configs) it is too noisy to use.

## C17. THE REQUIREMENT IS NOW A NUMBER — 0.638 in-band accuracy for 2.0 A, against 0.539 available

**DEMONSTRATED.** OBJ workstream, `s14/obj_noisy.py`, `s14/obj_headline.py`,
`s14/obj_FINDINGS.md`. This converts every negative in this sprint from "discrimination binds"
into a quantity with a target and a shortfall.

**The enumerated set was extended from 9 to 19 targets** (all ten n=10 targets,
**1.28 x 10^7 exactly-labelled structures**, all five folds covered). Ordering verified
bit-identical to `itertools.product`; energies match the Sprint 13 cache to 1e-4. The n=11
batch was launched and abandoned — the box ran 4-9 competing jobs all sprint and their
processes were measured at **3-5% of one core**. Fold 1 has only 2 held-out targets and that
weakness is stated in their file rather than hidden.

### The requirement, quantified for the first time

A noisy-oracle sweep prices what in-band ordering accuracy is needed to reach a given emitted
RMSD through a given terminal operator:

| | in-band pairwise accuracy |
|---|---|
| **required for 2.0 A through a top-100 operator** | **0.638** |
| Legacy | 0.539 |
| the learned objective | 0.523 |
| the 1-local prior | 0.501 |
| AMBER | **0.463** (below chance) |

**That is an order-of-magnitude gap in ordering information, not a tuning gap.** Every negative
in this sprint reduces to this line.

**And the terminal operator matters as much as the ranker.** A *perfect* ranker emits 2.474 A
through the decile and **1.219 A through top-100** — the same ranking, a 1.26 A difference,
purely from how the set is consumed. That is the Sprint 12 terminal-operator law reappearing
with a number attached, and it independently supports C11's finding that aggregation is worth
3.4x the objective.

### The learned objective is REFUTED, and the nulls are the reason to trust it

A 4,125-parameter sequence-conditioned distance-binned pair potential, leave-fold-out over 19
targets:

| arm | d_decile [CI] | in-band acc | global acc | W/L |
|---|---|---|---|---|
| learned | **-0.166 [-0.467,+0.119]** | 0.523 | 0.565 | 13/6 |
| sequence-BLIND twin | **-0.337 [-0.546,-0.126]** | 0.514 | **0.616** | 16/3 |
| trained in-band | +0.181 | **0.470** | 0.481 | 7/12 |
| permuted-label null | +0.314 [+0.091,+0.547] | 0.450 | — | 7/12 |
| **leaked-label control** | **-1.405 [-1.604,-1.235]** | **0.916** | 0.996 | **19/0** |

The learned arm's CI crosses zero and **drop-top-3 flips it to +0.019**. Its margins over both
nulls have *positive* drop-top-10, and both nulls are themselves worse than random — the
weak-control trap this project has documented before. The leaked-label control in the same
harness reaches -1.405 A on 19/19 with drop-top-10 -1.093, so the harness is loud when a
signal exists.

**The learning curve is not flat — it DECLINES.** Four points now:

| training targets | d_decile | rho_global | rho_decile |
|---|---|---|---|
| 1 | -0.348 | +0.367 | -0.050 |
| **2** | **-0.380** | **+0.411** | -0.059 |
| 4 | -0.231 | +0.229 | -0.003 |
| 8 | -0.230 | +0.229 | +0.033 |
| 12 | -0.162 | +0.216 | +0.010 |
| 15 | -0.166 | +0.166 | +0.048 |

(All five sizes landed; 28 observations per point.)

**More training data makes the objective monotonically worse** from m=2, on both the structural
axis (-0.380 to -0.162) and the objective axis (rho_global +0.411 to +0.216), while in-band
`rho_decile` never leaves zero at any size. **One or two training targets beat twelve.** That is the Sprint 12 flat-curve diagnostic in its
strongest possible form, and it disposes of the residual worry that the curve might be rising
and would keep rising past 15 targets: **it is falling.**

The mechanism was already measured. One target is enough to locate the single geometric
compactness axis; additional targets only average that axis toward the population mean, which
is exactly where the per-target sign disagreement destroys it — the native radius-of-gyration
z-score correlates **+0.909** with per-target skill. **Signal-limited, not sample-limited**, so
neither the 13 unrun targets nor a larger model changes the answer — which retrospectively
justifies abandoning them.

**THE LOUDNESS CALIBRATION LANDED, AND IT REMOVES THE LAST CAVEAT.** The leaked-label arm was
run at every training size, and the two curves **diverge**:

| m | learned d_decile | learned rho_decile | **leak d_decile** | **leak rho_decile** |
|---|---|---|---|---|
| 1 | -0.348 | -0.050 | **-1.273** | **+0.551** |
| 2 | **-0.380** | -0.059 | -1.397 | +0.627 |
| 4 | -0.231 | -0.003 | -1.407 | +0.683 |
| 8 | -0.230 | +0.033 | -1.419 | +0.710 |
| 12 | -0.162 | +0.010 | **-1.423** | **+0.728** |

**The leak arm is already loud at a SINGLE training target (-1.273 A) and RISES monotonically**
to -1.423, with in-band `rho_decile` climbing **+0.551 to +0.728** and global rho saturating at
+0.994. So **the harness is fully capable of expressing a rising learning curve on the in-band
axis** — it simply does not do so when the label is honest.

Against that, the learned arm declines with its in-band `rho_decile` **pinned at zero across a
12-fold change in training data**. No caveat remains on the diagnostic. This is Sprint 12's
decisive negative reproduced on a completely different instrument, and it is now the
**predicted consequence of C19** rather than an independent observation: in-band ordering is
learnable within a target (0.986) and not across targets (0.600), so adding targets averages
the axis away instead of sharpening it.

### Three refutations, two of them of the workstream's own hypotheses

- **Sequence conditioning is HARMFUL.** The sequence-blind twin *beats* the conditioned model
  by +0.172 A (W/L 9/10), and is the only arm clearing random (-0.337 [-0.546,-0.126]) and
  beating Legacy (-0.361 [-0.647,-0.076]). What it buys is **global garbage rejection**
  (0.616 against Legacy's 0.534), **not in-band ranking** (0.514 against Legacy's 0.539) —
  learning buys discrimination in the wrong place.
- **Training in-band does not create in-band skill.** The arm trained exclusively on in-band
  pairs scores **0.470, below chance**.
- **The flip diagnostic.** The model's per-target skill correlates **+0.909 with the native's
  z-scored radius of gyration** and **+0.951 with rho(contacts, RMSD)**. A 4,125-parameter
  potential learned exactly one geometric axis whose correct sign is a per-target property it
  cannot observe.

### Both of my warnings were handled, and one moved a number

The AMBER oracle-conditioning reproduced independently: **-0.405 A on the full subset against
+0.007 for `amber_kind == 0`**, and **AMBER's top-100 skill was 4x overstated** (3.223 to
3.694 A). Their learned objective is audited AMBER-free by grep across all six modules. The
four inert torsions were confirmed **bit-exactly** (rtol=0, atol=0, seven targets, with
residues 1 and n-2 as live controls); masking is now default and the ablation shows it is
worth ~0.00 A.

**They recommend the VQE workstream does NOT adopt their objective as a Hamiltonian**:
`pauli_support(9)` is already full-register (separation 2 needs 2 qubits, separation 8 needs
14), confirming rather than escaping the locality theorem. The one real handle is
`Featurizer(max_sep=s)`, which caps Pauli weight at exactly `2(s-1)`.

## C18. THE PARADOX IS AN EXACT IDENTITY — selection decomposes into FILTERING plus ORDERING

**DEMONSTRATED.** ENER workstream, `s14/ener_tail.py`, `ener_lib.tail_accuracy` and
`ener_lib.selection_decomposition`, written up as their E13. This is the deepest result of the
sprint and it subsumes several others.

**A methodological caution that must travel with every tail number.** They stated the
prediction before running: since discrimination is a function of the quality gap, and gaps
inside a tail are small *by construction*, every objective's tail accuracy must fall toward
chance as the tail shrinks **whether or not it lost skill** — Legacy's available |dRMSD| falls
from 0.920 A over the whole space to 0.416 A at q = 0.2%. So they carry a gap-matched variant
and a **random-tail null**, and **the null sits at 0.524-0.527, not 0.500**. Every tail
accuracy figure in this sprint — mine and the VQE workstream's included — must be read against
that baseline rather than against 0.5.

**The result survives it.** Legacy scores **0.462 gap-matched inside its own lowest 1%**
against **0.523** for Legacy on a random subset of the same size. Gap-matching does not rescue
it, so it is not the shrinking-gap artefact. **An energy is blindest inside its own tail.**

### The decomposition

    sel  =  pool  +  FILTERING  +  ORDERING          exact and additive

At q = 1%:

| objective | FILTERING | ORDERING | net |
|---|---|---|---|
| Legacy | -0.124 | **+0.263 [+0.046,+0.484]** | **+0.139** |
| `leg_steric` | -0.068 | -0.021 [-0.125,+0.082] | -0.089 |
| AMBER | **-0.141** (best measured; `amb_nonbonded` -0.163) | +0.186 | +0.045 |

**Legacy's net +0.139 is exactly the certified-optimum-versus-random number from Sprint 13.
The identity closes on it.** A result that had stood for a sprint as an empirical curiosity is
now the sum of two measured terms.

Three things follow, and each retires a separate open question:

1. **The ORDERING term's sign *is* "optimise harder, get worse", measured directly.** It is
   the single number to carry for any future objective, and it is positive for both physical
   energies — they actively mis-order inside the region a search occupies.
2. **`leg_steric` has a filtering component and no ordering component at all**
   (-0.021, CI spanning zero). Third independent confirmation that Legacy is a gate.
3. **AMBER has the best filtering term measured and a positive ordering term.** Third
   instrument to reach "validator, not ranker", from a third direction.

**And it dissolves the paradox in C15/C16.** Bulk discrimination and argmin quality are
**different terms of one identity**, not competing measures of one thing. The distogram having
the best of one and the worst of the other is not a contradiction and never was; neither term
should ever be reported as predicting the other. That is the cleanest available statement of
why four separate proxies failed to predict selection quality this sprint.

## C15-ATTRIBUTION-CORRECTION. My account of my own error was itself wrong

I wrote that my two-target discrimination sample (1CS9, 2MK7) was "precisely the two outliers
the VQE agent had already flagged". The ENER agent checked the overlap rather than accepting
the tidier claim, and it is wrong: the two targets carrying 68% of the `leg_torsion` effect
were **7N2I (-2.74) and 2MK7 (-2.36)**, not 1CS9. **Only 2MK7 is common to both failures**,
and 1CS9 is in fact where `leg_torsion` performed *worst*.

The defensible statement is the narrow one: **2MK7 is a recurring outlier, and n=2 samples
drawn from these nine targets have now reversed a conclusion twice.** My original error stands
as recorded; my explanation of it did not, and the correction is theirs.

## C19. THE FEATURE SPACE IS NOT THE PROBLEM. In-band discrimination is learnable to near-perfection per target and does not transfer.

**DEMONSTRATED.** `python -m s14.obj_ceiling` — written by the OBJ workstream, left unrun for
lack of an idle core, and run by the coordinator once the box freed up. **12** fully
enumerated targets (the run's own `n_expected`; the enumeration cache holds 19 and the
ceiling arm covers 12 of them). It is the highest-value experiment of the sprint and it answers the question every
other negative left open.

**The question.** Every objective this project has built fails *in band* — inside the
near-native region where a search actually lives. Is that because the feature space cannot
represent in-band discrimination, or because whatever discriminates in one target does not
discriminate in another? The first would mean a richer model is worth building. The second
would mean no model is.

Three arms, identical features, identical training, differing only in what is held out:

| arm | in-band pairwise accuracy (< 1.5 A) | rho | d_top100 | W/L vs chance |
|---|---|---|---|---|
| in-sample (upper bound) | **0.986 [0.974, 0.997]** | +0.869 | **-0.947 A** | 12/0 |
| **same-target held out** | **0.986 [0.973, 0.997]** | +0.868 | **-0.868 A** | 12/0 |
| **cross-target** | **0.600 [0.549, 0.655]** | +0.111 | **-0.047 A** | 10/2 |

(Band mean RMSD 2.195 A, band minimum 1.031 A.)

**The feature space is emphatically not the barrier.** It reaches **0.986** in-band accuracy —
against the **0.638** that 2.0 A requires — and, decisively, **held-out configurations of the
same target score identically to in-sample ones**. The two gaps, computed from the raw rows by
the workstream that owns the module:

| gap | value | CI |
|---|---|---|
| **overfitting** (same-target held out − in-sample) | **-0.0005** | [-0.0011, +0.0000] |
| **transfer** (cross-target − same-target held out) | **-0.3859** | [-0.4366, -0.3320] |

**The overfitting gap is statistically indistinguishable from zero.** That is what makes the
`d_top100 = -0.868 A` figure unarguable: the model genuinely learns a target's in-band ordering
and generalises it perfectly to unseen configurations *of that target*. The entire loss is the
transfer term, and it is 770 times larger.

**And it transfers essentially not at all.** Cross-target accuracy is **0.600**, below the
0.638 requirement, and it converts to **-0.047 A** — nothing. The gap between 0.986 and 0.600
is the entire finding.

### What this settles

**A richer model is not worth building.** A linear pair potential already saturates the
within-target problem at 0.986. Nonlinearity, more parameters, equivariance, graph or
geometric architectures — none of them can improve on 0.986, and none of them addresses
transfer, which is where the loss is. That closes the "train a bigger model" branch cleanly,
and it closes it with a positive measurement rather than a failure to find something.

**It converges with three other measurements onto one mechanism.** The declining learning curve
(one or two training targets beat twelve); the flip diagnostic (per-target skill correlating
**+0.909** with the native's radius-of-gyration z-score and **+0.951** with rho(contacts,
RMSD)); and now the transfer collapse. All three say the same thing: **the in-band
discriminating axis is real, learnable, and per-target — its correct sign is a property of the
individual target that inference cannot observe.** Averaging across targets cancels it, which
is why more data makes the objective worse.

**And it explains the tail result from the other side.** C16 found every objective at chance
inside its own tail. C19 says that is not because tail ordering is unlearnable — it is learnable
at 0.986 — but because the learned ordering is target-specific and a deployed objective must
carry one fixed ordering for all targets.

This is the sprint's deepest finding. The bottleneck is not representation, not capacity, not
search, and not sample size. **It is that the quantity which orders near-native structures
changes sign from target to target, and nothing available at inference reveals which sign
applies.**

### Two comparability cautions, both raised by the workstream against its own numbers

**The 0.524-0.527 tail null does NOT apply here.** The ENER caution is correct for *their*
shape of statistic — unbinned accuracy on the lowest 0.1%, where the available spread is very
tight. The OBJ statistic **conditions on the ΔRMSD bin**, so pairs are matched on separation
and the null returns to chance. Measured with a purely random objective over 8 targets:

| statistic | measured null |
|---|---|
| binned < 1.5 A on the low-energy decile (their headline arm) | **0.5006** (sd 0.0013) |
| binned < 1.5 A on the near-native band (the C19 statistic) | 0.5050 (sd 0.0072) |
| unbinned on the near-native band | 0.5025 (sd 0.0031) |

So Legacy 0.539, learned 0.523, sequence-blind 0.514 and AMBER 0.463 all stand against a true
0.500. **The two nulls must not be conflated**, and the file now says so.

**C19's cross-target 0.600 must NOT be read against the headline 0.523.** Different
populations: C19 trains and evaluates on the oracle-selected near-native band, while the
headline arm trains on uniform samples and evaluates on the objective's own low-energy decile
of the whole space. **The structural number that travels with 0.600 is `d_top100 = -0.047 A`**,
which is nothing — that, not the accuracy figure, is what should be quoted.

**What C19 leaves genuinely open**, recorded as the one live direction: anything supplying the
**per-target sign** of the in-band axis at inference — a *conditioning* signal, not a bigger
objective. Their `noint` arm was the cheap version of exactly that and was null, and it is
underpowered at 19 targets.

## C20. THE DEFINITIVE NEGATIVE — running the VQE is worse than not running it, 0/12

**DEMONSTRATED.** VQE workstream initialisation experiment. This is the strongest single
result in the sprint and it answers the causality question outright.

**The control is the whole point.** "After VQE" is compared not against the initialisation's
*mean* but against **`init_best`: the best of the SAME NUMBER OF SHOTS drawn from the INITIAL
distribution with no optimisation whatsoever.** That is the arm a VQE must beat to have
contributed anything, because best-of-N sampling is free and is what the shipped pipeline
already does. `DELTA = after − init_best`; only a negative DELTA is a contribution.

Signal 0.10, budget 10,240 evaluations, 4 targets x 3 seeds, certified optimum 2.309 A:

| initialisation | init mean | init best | after VQE | **DELTA** | CI95 | W/L |
|---|---|---|---|---|---|---|
| random (native-free) | 4.131 | 1.797 | 2.606 | **+0.808** | [+0.576,+1.072] | **0/12** |
| uniform/zero (native-free) | 4.119 | 1.765 | 2.757 | **+0.992** | [+0.756,+1.242] | **0/12** |
| empirical prior warm start (native-free) | 4.062 | 1.591 | 2.469 | **+0.879** | [+0.539,+1.237] | 1/11 |
| ORACLE prior q=0.50 | 3.893 | 1.309 | 2.626 | **+1.317** | [+0.888,+1.806] | **0/12** |
| ORACLE prior q=0.80 | 3.424 | 1.444 | 2.714 | **+1.271** | [+0.820,+1.796] | **0/12** |
| ORACLE prior q=0.95 | 2.966 | 1.764 | 2.578 | **+0.813** | [+0.274,+1.455] | 2/10 |

Signal 0.30 reproduces it: DELTA +0.654 to +1.118, 0/12 on four of six rows.

> **Taking the best of 10,240 shots from the UNTRAINED circuit beats running the VQE on the
> same budget, every time, by 0.65 to 1.32 A.**

### The mechanism is C16, and the loop closes

The optimisation concentrates the distribution onto the objective's low-energy tail. C16 showed
**every objective is at chance inside that tail** (0.512 at w=0) and **the argmin is a random
draw from it** (41st-48th percentile). So the optimiser trades a broad sample, whose *best*
member is good, for a narrow one whose *mean* it then samples blindly.

**Concentration is exactly the wrong move when discrimination is the binding constraint.**

That is one finding seen from two sides — C16 from the mechanism, C20 from the outcome — and
they should be presented together. It also retroactively explains C11: aggregation beats the
objective 3.4 to 1 because a broad sample is worth more than a concentrated one here.

### The initialisation confound is answered, in the direction least convenient for VQE

The brief's section 21 warns against an initialisation that merely hands VQE the answer.
Measured: **after VQE the result is 2.42-2.76 A regardless of where it started**, from
initialisation means spanning 2.966 to 4.131 A. An ORACLE warm start at q=0.95 (init mean
2.966, init best 1.764) ends at 2.578 — no better than a random start's 2.606. **The optimiser
washes out its own initialisation completely.**

Note the perverse ordering: **the better the oracle initialisation, the worse the DELTA**
(q=0.50 is worst at +1.317), because a better start has a better `init_best` while the VQE
converges to the same place regardless.

So there is no version of this experiment in which a warm start flatters the VQE. It is neither
handing over the answer nor accelerating genuine search — **it is erased.**

### Two process notes, both recorded

The agent **stopped `s14/vqe_ansatz.py` and ran a reduced form instead**, because the box had
fallen to roughly 0.07 cores per Python process with total CPU pinned at 100% by non-Python
work, and that module writes its JSON only on completion — it would have produced nothing.
Its most valuable result was captured standalone (the depth-1-only QNG correction, C-brief).
Also measured there: **the product ansatz is rank 8 with 16 or 24 parameters — exact L-fold
parameter redundancy — while every entangling pattern is full rank.**

Whether QNG *helps* at depth >= 2 was not measured and is marked HYPOTHESIS, not a finding.

## C21. CVaR correctness closes — and a claim I carried in the brief has its SIGN BACKWARDS

**DEMONSTRATED.** VQE workstream final report; 92 correctness tests pass.

**The value estimator is correct**: 1.7e-15 against an independent 60-digit reference,
including fractional partial buckets, tie permutation-invariance, and the `alpha*2^n < 1`
edge. `grad_cvar_paramshift` agrees with finite differences at cosine **0.999999999999998**,
worst relative error 2.85e-09, across n x depth x alpha x seed. `baseline="const"` is unbiased
in direction and magnitude; `baseline="tail"` is biased and `qansatz.cvar_gradient` reproduces
that bias to **0.00e+00**.

### Three things the prior record got wrong or missed

**1. The known defect's bias has a CLOSED FORM.** It is `-c * grad P(E < q)`, verified to
< 1e-9 — **it adds the gradient of the tail PROBABILITY, not noise.** The long-recorded
"cosine -0.023 with the true gradient" is **one draw from a spread measured at +0.06 to
+0.96**. That figure should never again be quoted as a characterisation of the defect; the
closed form replaces it.

**2. A new sampled-CVaR bias** at non-integer `alpha*N`: +0.134 sd at N=13, alpha=0.1,
decaying as 1/N. Recorded in C-brief.

**3. `dCVaR/dp` is identically zero IFF `p(argmin E) >= alpha`** — 0 counterexamples in 3,000
cases.

### THE SIGN CORRECTION, and it is mine to own

The shared brief and the Sprint 13 record both carry: *"CVaR at alpha <= 0.25 makes `amber`
and `amber_soft` literally the same objective"*, filed as a trap in which small alpha destroys
discrimination between objectives. **The sign is backwards.**

Measured: **small alpha DISCRIMINATES BETTER between objectives. It is CONCENTRATION that
collapses them.** The collapse is real, but its cause was mis-attributed to the alpha level
rather than to the distribution narrowing during optimisation. I propagated that framing into
the Sprint 14 brief that four agents worked from, and it is corrected here and in the brief.

Note this is consistent with — and explains — C20: concentration is the destructive operation,
not small alpha per se. The two records now agree.

### One-hot's spectrum measures its constraint, not its objective

Binary attains the information-theoretic bound exactly (18 qubits at n=9, k=4); everything else
is strictly wider. Binary needs **510 Pauli terms** for 99% of variance where one-hot needs
**454,463**. And one-hot's spectrum sits **exactly on the Binomial(m, 1/2) null and is immovable
under a penalty sweep** — it is measuring the feasibility constraint, not the objective. That
invalidates naive cross-encoding spectrum comparisons and is a methodological trap worth
carrying forward.

### Their own refutations, recorded

Gray-coding Pauli advantage (clean on one target, dead over 54 cells); **a structural crossing
at rho ~ 0.72** (CI excluded zero at n=4 and was gone by n=7 — a near-miss positive killed by
its own replication); "alpha is a learning rate"; their pre-registered tail-concentration
prediction, refuted at its premise by their own measurement; and their argmin-based
recommendation to drop the distogram, corrected by a top-m Jaccard test **before the
coordinator acted on it**.

### C21b. The full crossing grid, and a near-miss positive dying monotonically in n

The Part C grid reached all nine targets. **There is no crossing on either axis.**

- **Objective axis:** VQE loses at *every* signal level, with the CI excluding zero on the
  wrong side at **six of seven levels** — W/L 0/7, 1/7, 0/6, 1/6, 0/7, 0/4.
- **Structure axis:** significantly *worse* at signal 0.20 (+0.070 [+0.002,+0.142]) and 1.00
  (+0.030 [+0.004,+0.062]); every other cell spans zero. **No cell favours VQE.**

And the one marginal "crossing" the workstream had briefly written up at rho ~ 0.72 is now
decisively dead, with its death **monotone in n**:

| n targets | signal-0.50 structure difference | CI95 | W/L |
|---|---|---|---|
| 4 | -0.113 | [-0.292, **-0.000**] | 3/1 |
| 6 | -0.053 | [-0.226, +0.071] | 3/2 |
| 7 | -0.044 | [-0.187, +0.062] | 3/3 |
| 8 | -0.037 | [-0.167, +0.055] | 3/4 |
| **9 (complete)** | **-0.022** | [-0.122, +0.051] | **2/5** |

The effect shrank **fivefold** and the win/loss inverted from 3/1 favourable to 2/5
unfavourable. **A real effect does not decay like that — the shape of the decay is the
diagnostic**, and the whole progression is kept in the findings rather than only the corrected
number. **The concentration
hazard landed precisely on the one cell that would have been a positive headline** — a CI whose
upper bound touched zero at n=4. That table is worth more than the corrected number alone, and
it is the fourth time this sprint that an n=4-or-smaller sample produced a conclusion which
reversed on more data.

### C21c. An exact cross-sprint reproduction, and the budget trap is cured by OBJECTIVE QUALITY

The Part C grid completed: all 9 targets, 92 cells. Two things the completion added.

**An exact reproduction that validates the whole instrument.** At signal 0 the certified global
optimum is **3.920 A** against a random draw of **3.781 A** — **+0.139 A worse**. That is
Sprint 13's recorded figure **to three decimal places**, reproduced by independent code on the
same nine spaces. It is the tightest available check that this sprint measured what Sprint 13
measured, and it retroactively underwrites the ENER identity in C18, whose ORDERING term nets
to that same +0.139.

**The budget trap is a property of BAD OBJECTIVES, not of the problem — and this refines C8.**
At full n the non-monotone budget curve **inverts by signal 0.50**, where budget 30,000 reaches
**2.149 A against a certified optimum of 2.153 A** and more search is **strictly better**.

So the complete statement, which supersedes both the Sprint 13 warning and my own C8 framing:

> Searching harder makes structures worse on a bad objective, neither helps nor hurts on a
> mediocre one, and **helps monotonically on a good one**. The trap was never a property of
> search; it was a property of what was being searched.

That does not rescue VQE — at the signal levels achievable with real native-free information
the curve is in the flat-to-harmful regime, and VQE loses to greedy on the objective axis
everywhere regardless. But it means "optimise harder, get worse" must always be stated with
its objective-quality condition attached, and the Sprint 13 memory should carry that condition.

`s14/vqe_FINDINGS.md` is final at 1,231 lines, 14 result JSONs, 9 refutations and retractions
recorded in place.

## C22. The one remaining direction is NOT closed — but the native-free sign signal is weak

**DEMONSTRATED.** `s14/signpred.py`, all 126 targets. This is the sprint's last experiment and
it tests the only direction C19 left with leverage.

**The logic.** C19 showed in-band ordering is learnable to 0.986 within a target and collapses
to 0.600 across targets, so the missing ingredient is a **conditioning signal**: something that
supplies the per-target SIGN of the ordering axis at inference. The flip diagnostic identified
that axis as compactness — per-target skill correlates **+0.909 with the native's z-scored
radius of gyration**. But that radius is a *native* quantity, and the whole question is whether
any native-free predictor of it exists.

**One does, already in the pipeline and never read this way.** The shipped leave-fold-out
distogram predicts every CA-CA pair distance, and a radius of gyration follows in closed form:

    Rg^2 = (1 / (2 N^2)) * sum_ij d_ij^2

So the distogram's `expected` column *is* a native-free Rg predictor, for free. Two further
native-free proxies: the mean Rg of the retrieved window pool, and the Rg of the incumbent's
own emitted structure.

**The decisive test is length-residualised**, because compactness scales with chain length and
length is trivially known — the question is whether a proxy knows *this target's* compactness
beyond how long it is.

| proxy | raw Pearson | **length-residualised** | CI95 |
|---|---|---|---|
| chain length alone (control) | 0.189 | — | [+0.010, +0.375] |
| retrieval pool mean Rg | 0.257 | **0.244** | [+0.031, +0.444] |
| distogram-predicted Rg | 0.361 | **0.313** | [+0.110, +0.518] |
| **incumbent emitted Rg** | 0.404 | **0.365** | [+0.155, +0.577] |

**All three CIs exclude zero.** Native-free proxies carry **real information** about the native
radius of gyration beyond chain length. The direction is **not closed.**

**But the signal is weak, and the gap to what is needed is large.** The best proxy reaches
0.365 where the oracle version of the same quantity reaches **0.909**. Naively chaining
proxy → native Rg → sign gives a composite around 0.33, which for a binary sign corresponds to
roughly **61% accuracy** — real, but far from the 0.986 the within-target model achieves once
it knows the sign, and this project has established that fusion gain goes as the **square** of
the weaker channel's skill.

**Three honest limitations, stated rather than buried.**

1. **The chain was not measured end to end.** I measured proxy → native Rg on 126 targets, and
   the +0.909 skill ↔ native-Rg link was measured by another workstream on a smaller enumerated
   subset. Composing them assumes linearity and independent errors. **The composite number
   above is an estimate, not a measurement** — the direct experiment is proxy → per-target
   sign → emitted RMSD, and it was not run.
2. **The best proxy is partly circular.** The incumbent's emitted Rg is native-free and
   therefore legitimate as a conditioning signal, but it is a property of the very output the
   conditioning is meant to improve. Its 0.365 should be treated as an upper bound among these
   three rather than a clean independent channel.
3. **The `i,i+2` term in the Rg reconstruction is an ideal-geometry approximation** (6.00 A).
   It is applied identically to every arm, so it cannot manufacture a correlation *across*
   targets, but it does add a small common-mode error.

**The verdict, and it is the right note to end the sprint on.** The one direction with leverage
survives its first test. The signal is real, native-free, free to compute, and weak. Whether
0.24-0.37 of compactness information converts into usable ordering is the next sprint's first
question.

**C22 IS SUBSTANTIALLY WEAKENED BY C23 BELOW. Read them together.**

## C23. The direct chain does NOT reproduce — and neither does its ORACLE ceiling

**DEMONSTRATED (a negative).** `s14/signchain.py`, the 12 targets where both halves exist.

C22 composed two separately-measured correlations into an estimate and I flagged that
composition as an assumption. Both halves were on disk, so I measured the chain **directly**
rather than leaving the estimate standing. The result corrects C22 and it corrects it downward.

The question, asked without composition: **does a native-free compactness proxy predict the
per-target sign of the learned objective's cross-target skill?** Skill measured as `cross_rho`
from the ceiling experiment; everything length-residualised; permutation p-values, which are
the right test at n=12.

| arm | Pearson | Spearman | perm p | CI95 |
|---|---|---|---|---|
| **A. ORACLE ceiling — native Rg** | **+0.416** | +0.343 | **0.177** | [-0.251, +0.838] |
| B. distogram-predicted Rg | +0.257 | +0.252 | 0.419 | [-0.287, +0.668] |
| B. retrieval pool mean Rg | -0.044 | +0.049 | 0.888 | [-0.598, +0.490] |
| B. incumbent emitted Rg | +0.196 | +0.189 | 0.537 | [-0.363, +0.658] |

Sign accuracy against a 0.500 majority-class baseline: distogram 0.667, emitted Rg 0.667, pool
0.500, and **the ORACLE native Rg itself only 0.500**. At n=12, 8/12 is not significant.

**Two findings, and the first is the important one.**

**1. The ORACLE ceiling does not reproduce at +0.909 in my hands — it measures +0.416 with
p = 0.177 and a CI spanning zero.** I do not claim the OBJ workstream's number is wrong. The
most likely explanation is that we measured **different quantities**: their flip diagnostic
correlated the *learned model's per-target skill in the headline arm* against the native Rg
**z-score**, while this correlates the *ceiling experiment's `cross_rho`* against the
length-residualised native Rg. Those are different skill measures on different arms, and a
large discrepancy between them is possible without either being in error.

But the consequence for C22 is the same either way: **C22's premise — that per-target skill
tracks native compactness at +0.909 — is not the quantity I then built a chain on.** The
ceiling for *this* chain is +0.416, not significant at n=12, and C22's "roughly 61% sign
accuracy" estimate should not be quoted.

**2. The chain itself is not demonstrated.** Every native-free proxy has a permutation p above
0.4 and a CI spanning zero. The direction is not refuted — n=12 has almost no power, and C22's
126-target proxy-to-native-Rg link stands on its own with all three CIs excluding zero — but
**it is not supported by the only end-to-end measurement available.**

**The honest closing statement, which replaces C22's:**

> Native-free proxies carry real information about a target's compactness (0.24-0.37 on 126
> targets, all CIs excluding zero). Whether compactness predicts the per-target ordering sign
> is **unresolved**: one workstream measures +0.909 on one skill definition, an independent
> implementation measures +0.416 with p = 0.177 on another, and the end-to-end chain from
> native-free proxy to skill is null at n=12. **The next sprint's first job is to resolve that
> discrepancy on a common definition and a larger enumerated set, before building anything.**

This is why the direct experiment was worth running: it took an attractive closing story and
replaced it with an open question that is honestly stated. A composed estimate would have
carried "61% sign accuracy" into the next sprint as though it were measured.

---

## Open at this point

---

## Open at this point

- Does a learned many-body objective have a structurally accurate low-energy region (OBJ).
- At what objective quality does VQE/CVaR start to beat classical search, and does it ever
  beat it on structure rather than energy (VQE).
- Whether shift data exists at all for these targets, and under which coherence mode (SHIFT).
- Each energy model's legitimate role and the decoy-discrimination threshold (ENER).
- Whether the September 2026 continuous-torsion preprint changes the direction (LIT).
