# SPRINT 24 — WORKSTREAM B (THE RESIDUAL GENERATOR) — FINDINGS

**Bottom line, stated first because it damages my own lane.** The residual channel is real, its
ORACLE ceiling is large, and **essentially none of that ceiling is reachable at any torsion accuracy
this project has ever achieved.** With the residual DIRECTION given by the native, the per-residue
confidence gate given by the native, and the amplitude and gate threshold chosen in-sample on the
same 126 targets they are scored on, the arm buys **−0.032 Å at 0.71× its own MDE** at the accuracy
a real model reaches. Every one of those four concessions is an oracle. **No model was trained, and
on this evidence none should be.**

Basis: point cloud on every side of every contrast. Baseline for every residual arm is `P0R` — the
same 75 retained windows rebuilt from their own torsions, no residual — never the real-window
incumbent, so no arm is charged for the rebuild operator. Pre-registered in `s24/PREREG_B.md`.
Artefacts: `s24/results/resid0.json`, `resid0b.json`, `resid0c.json`, all n=126, all `complete`.

**Nothing here consumed Lane A's corpus.** Lane A's exclusion list did not exist at any point during
this work, no training was attempted, and every number below is priced from ORACLE labels already in
the universe files, used for evaluation and ceiling-pricing only. `LOCK_TRAIN` was never taken — all
three runs are pure numpy.

---

## B-0 — THE SUBSTRATE IS FREE (n=126, complete)

The architecture must rebuild each retrieved window from its torsions before a residual can be
applied, and an ideal-geometry rebuild of a real window is not that window. That cost is the floor
of the lane.

    incumbent, real windows              3.0483   <- reproduces the pinned constant exactly
    P0R, rebuilt from own torsions       3.0524
    rebuild - incumbent   +0.0041  SE 0.0047  MDE 0.0131  57W/69L   NULL

Mean |rebuilt window − its own real window| is **0.3730 Å and it cancels in the average** — the same
variance-reduction mechanism s23 L9 identifies. A non-null here would have closed the lane on
arithmetic; it is reported as the null it is.

---

## B-1 — THE DECOMPOSITION. THE COMMON COMPONENT IS THE WHOLE LANE, AND THE STOCHASTIC PART IS WORSE THAN NEUTRAL

Pre-registered derivation: the readout is a uniform mean of 75 members, so to first order it sees
`mean_k d_k`; an i.i.d. residual should therefore contribute ~0 and only the conditional mean can
move the answer. **Measured, and the derivation is confirmed with a correction that makes it worse
for my architecture.**

    arm (all ORACLE)                RMSD      vs P0R      eff/MDE    W/L      cosR    nrm
    ORAC a=0.10                   2.9819    -0.0705       1.77x    93/33    +0.988   0.968
    ORAC a=1.00                   0.3468    -2.7056       6.64x   125/ 1    +0.073   0.160
    COMM a=0.10 (common only)     2.9983    -0.0541       1.30x    87/39    +0.984   0.974
    COMM a=1.00                   2.2357    -0.8167       1.75x    88/38    +0.562   0.759
    IDIO a=0.10 (member-specific) 3.0834    +0.0310       0.87x    50/76    +0.994   1.010
    IDIO a=1.00                   4.0725    +1.0201       2.12x    41/85    +0.779   1.333

`ORAC(1.00) = 0.3468` reproduces the project's known native continuous-torsion rebuild floor
(0.347), so the readout does not destroy a perfect residual and the channel's headroom is genuine.

**COMM carries 77% of the oracle gain at small amplitude** (−0.054 of −0.071). **IDIO — the
member-specific component alone — is strictly WORSE at every amplitude**, monotonically, out to
+1.02 Å.

> **The correction to my own derivation, and it is the finding that matters.** I predicted the
> i.i.d. component would contribute ~0. It does not: it contributes **positive damage that grows
> with amplitude**, and the bias NORM grows with it (nrm 1.010 → 1.333). The mechanism is that
> `IDIO` is mean-free in TORSION space, and the torsion→coordinate map is nonlinear, so a
> **mean-free torsional perturbation is not a mean-free coordinate perturbation** — chain
> integration converts per-member torsional scatter into coherent coordinate distortion. This is
> the same physics as the project's measured +1.024 Å penalty for averaging in torsion space rather
> than coordinate space, arrived at from the perturbation side.
>
> **Consequence for the mandated architecture: the stochastic part of a stochastic torsional
> residual cannot help through this readout and actively hurts.** A model emitting many samples per
> candidate is emitting, in the only component the readout can see, its conditional mean — plus a
> penalty for the spread around it.

**Normalisation fork settled by measurement, not assumption.** `CALT` (residual of the circular
mean) reproduces `COMM` (circular mean of the residuals) to four decimals at every amplitude. The
sin/cos-vs-raw-angle choice the directive proposes is **not load-bearing here**; I tested it rather
than assuming it, and it does not discriminate.

### The matched-accuracy coherence null — an independent replication of a project law

Both arms carry the full oracle direction corrupted at the same per-member angular accuracy,
differing only in whether the mistakes are shared:

    sigma   COH (shared mistakes)      IID (independent mistakes)     gap
    20 deg      -0.654                       -1.971                  1.32
    30 deg      +0.275   <- crosses          -1.140                  1.41
    45 deg      +1.141                       -0.109   <- crosses     1.25
    60 deg      +1.568                       +0.485                  1.08

At identical accuracy, **coherent mistakes cost 1.1–1.4 Å more than i.i.d. mistakes**, and the
crossing points differ by 20°. This is `error-coherence-decides-correctors` reproduced on a
completely different instrument (torsional residuals through a coordinate average, vs per-pair
distance sign correctors), same direction, roughly **3× the effect size** the memory records (1.41
vs 0.44 Å). The zero-information control — the oracle residual's own values with residue positions
permuted, size-matched in the operator's own space — is +0.760 coherent / +0.426 i.i.d., confirming
the ordering holds when the direction is destroyed.

**This is the axis that closes the lane.** A residual model conditioned on TARGET-LEVEL information
— sequence embedding, distogram, fold — makes target-level mistakes, and those mistakes are shared
by all 75 members **by construction**. So a learned residual lives on the COH curve, which crosses
the do-nothing baseline at **≈27°**.

---

## B-1b — THE JOINT (AMPLITUDE × ACCURACY × COHERENCE) GRID: THE CLOSING MEASUREMENT

Full amplitude is not a fair test — a sensible model shrinks under uncertainty. So the honest
question is joint. Noise coherence `rho` interpolates between the two endpoints, variance-matched
throughout; the amplitude is chosen ORACLE, in-sample, per cell.

    ORACLE-best alpha per cell, vs P0R          rho=0 (optimistic)   rho=1 (what a model does)
    sigma = 30 deg                                -1.105  2.89x        -0.065  0.16x   ns
    sigma = 45 deg                                -0.450  1.30x        -0.040  0.79x   TYPE-M
    sigma = 60 deg                                -0.103  0.60x  ns    -0.001  0.03x   ns
    sigma = 70 deg                                -0.065  0.55x  ns    -0.010  0.23x   ns

**The accuracy this project achieves is sigma = 69.7°** over all determined angles, with phi at MAE
36.1° against a **sequence-BLIND** marginal's 36.4° (`phi-carries-no-sequence-signal`). Reading
those as sigma, phi sits near 45° and psi near 78°.

> **At rho = 1 and sigma = 70 — a target-conditioned model at the accuracy this project actually
> reaches — the ORACLE-best amplitude buys −0.0101 Å at 0.23× its own MDE with a fold CI spanning
> zero. A dead null.** Granting the model phi-quality accuracy on BOTH angles (sigma = 45) lifts it
> to −0.0399 Å at 0.79× MDE, which the BRIEF defines as the Type-M zone and not a result.

### The coordinator's L3 diagnostic, answered directly

The question was whether the conditional-mean residual points AGAINST the incumbent's bias
(correcting) or along it (rescaling). Measured with `qmatch._bias` arithmetic, byte-identical:

    COMM a=0.10   cosR +0.984   nrm 0.974     rho=1 sigma=70 best cell   cosR +0.971  nrm 1.027
    COMM a=0.30   cosR +0.926   nrm 0.941
    COMM a=1.00   cosR +0.562   nrm 0.759

**At every achievable amplitude the answer is "rescaling", not "correcting"** — cosine against the
unmodified cloud stays at 0.92–0.99 with the norm essentially unchanged. The bias only genuinely
rotates (cosR 0.56) at oracle amplitudes that require knowing the native. By your own criterion,
this is the outcome that is **not worth training a model for**.

---

## B-1c — THE LAST DOOR: AN ORACLE PER-RESIDUE CONFIDENCE GATE

A perfect confidence gate is the one construction that has ever partly rescued a torsion channel at
this length, so the lane deserved to be closed on its strongest form. The gate is ORACLE — it keeps
exactly the residues where the realised error is under tau — and tau and alpha are chosen in-sample.

    sigma = 70, rho = 1     ungated (tau=180)  +0.026 at a=0.10, +0.380 at a=0.30
                            ORACLE-best cell: tau=15 (keeps 16.2%), a=0.30
                            -0.0317  SE 0.0160  MDE 0.0448  0.71x  fold[-0.039,-0.024]  78W/48L
    sigma = 45, rho = 1     ORACLE-best cell: tau=15 (keeps 27.3%), a=0.30
                            -0.0747  SE 0.0218  MDE 0.0610  1.22x  fold[-0.099,-0.052]  81W/44L

The gate is worth something real — it converts a +0.026 ungated arm into −0.032 — but **the ceiling
of the entire lane, with the direction, the gate, the amplitude and the threshold ALL oracular, is
−0.03 to −0.07 Å in the Type-M zone.** The sprint needs −0.049 Å from 3.0483 to reach 3.0. My
lane's oracle ceiling straddles the sprint's entire target, and every step down from oracle is
worse.

### The stereochemical trap in shrinkage, which is a genuine mechanism and cuts against the lane

    Ramachandran favoured:   P0R 0.955      ORAC a=0.30 0.811      COMM a=0.30 0.859
    clashes per structure:   P0R 0.11       ORAC a=0.30 0.71       COMM a=0.30 0.64

**A partially-applied torsional residual is less chemically valid than either endpoint.** Shrinking
the amplitude interpolates between two Ramachandran basins and lands the residue *between* them —
the same failure mode as a circular mean of a bimodal distribution (s23 L11a). So the two
requirements pull in opposite directions: **amplitude shrinkage is what makes the residual safe in
RMSD, and it is exactly what makes it chemically invalid.** Any residual model that emits a
posterior mean inherits this. Reported per the standing rule that RMSD bought with absurd geometry
is not bought — here it is the *safe* arm that is absurd, which is the opposite of the usual trap.

---

## A NEW FACT WORTH KEEPING: THE TEMPLATE IS A WORSE PHI PREDICTOR THAN NO INFORMATION AT ALL

    mean |wrap(native phi - retained member phi)|  =  40.6 deg     (n=126, 75 members each)
    sequence-BLIND corpus marginal, s13            =  36.4 deg

**The retrieved window's own phi is a worse estimate of the native phi than a Ramachandran marginal
that sees nothing.** So in `template + residual`, the template is not an asset on the phi channel —
conditioning on it starts the model behind a constant prior. This is a new measurement and it
independently supports `phi-carries-no-sequence-signal` from the retrieval side.

---

## THE STRUCTURAL ARGUMENT, RESTATED SO IT SURVIVES

`d*_k = psi_native − psi_k`, and `psi_k` is **fully observed at inference**. Therefore
`E[d*_k | features] = E[psi_native | features] − psi_k` **exactly**. A torsional residual model IS a
native-torsion predictor with a known offset subtracted; the residual framing adds no information of
its own. What it adds is conditioning s13's predictor lacked — the template's torsions and the
distogram — and B-1b prices what that conditioning would have to be worth: it would have to reach
**sigma ≈ 27° with shared mistakes** to break even, against 69.7° achieved and 36° for phi from the
full sequence versus 36.4° from nothing.

## AUDIT HARNESS (built, exercised, and available to Lane C)

`s24/residlib.py` carries the mode-collapse audit (unique structures at 0.10 Å single linkage,
duplicate fraction, ESS from cluster-occupancy entropy, max mode occupancy, pairwise RMSD, circular
torsion entropy) and the validity audit, plus per-comparison MDE, iid **and** fold-clustered CIs,
W/L, worst-target degradation and the drop-top concentration curve. Exercised on the 75-member sets,
which gives the reference values any generator must match:

    P0R   unique 66.8/75   dup 0.110   ESS 63.7   max-mode 0.044   pairRMSD 2.528   tors-ent 0.261

Recorded so the harness is not rebuilt: for structures from `core.geometry.build_backbone`, omega
deviation, cis fraction, bond-length and bond-angle deviation are **0 BY CONSTRUCTION** and Cα–Cα is
3.80 Å by construction. The validity axes that carry information for any torsion-space generator are
**Ramachandran and clashes**, and they are the two that moved here.

## RECOMMENDATION

**Do not train the residual generator.** The lane is closed by an upper bound, not by a failed
attempt: four separate oracles, stacked in the arm's favour, on the full 126-target instrument, land
in the Type-M zone. The compute is better spent elsewhere, and I would rather hand back a closed
lane with a mechanism than a trained model that lands on +0.002 Å.

If any part of this is to be reopened, the informative direction is **not** a better residual model.
It is that `COMM` — a single common correction shared by all 75 members — is where 77% of the
channel lives, which says the object worth predicting is a **per-target common-mode offset**, not a
per-member residual. That is the same quantity s23 L9 proved is invisible from inside the pool and
s23 L6/L9 priced at −0.3403 Å ORACLE and unreachable in principle. My lane rediscovered that
quantity from the generation side and hit the same wall.

### What I did NOT measure, stated so nobody assumes it

* Heteroscedastic accuracy with a **native-free** confidence estimate. B-1c grants an ORACLE gate;
  a real gate is worse, so this loosens nothing, but the *shape* of a real model's error may differ
  from the isotropic Gaussian used here.
* A residual conditioned on per-member information rich enough to make its mistakes genuinely
  i.i.d. (rho ≈ 0). B-1b prices that optimistic bound at −0.065 Å, 0.55× MDE, ns — so it is priced,
  but no such conditioning has been demonstrated to exist.
* Anything downstream of re-scoring. My fork list fixes the top-75 membership BEFORE the residual,
  deliberately; a re-scoring variant is a different experiment and would land on the L3/§1.3 flat
  lever.

---
---

# B-3 — THE PRIOR-ERROR ATTRIBUTION LADDER (new assignment, n=126, complete)

Artefacts `s24/results/priorladder.json`, `poolcheck.json`, both `complete` on the full key set and
both provenance-stamped via `ST.save_atomic`. `resid0/0b/0c.json` retro-stamped, with an explicit
`provenance_note` recording that they are retro-stamps and that the source was unmodified between
run and stamp — stated rather than passed off as an at-write stamp.

**ORACLE TRANSFER FUNCTION, not a system result.** Every arm interpolates the prior toward the
native's own distances. Same K=500 pool, same shipped functional, same top-75, same uniform
coordinate average. **Only the prior changes.**

## The two hazards, addressed before the run

**(1) Is the shipped loss exactly reconstructible? YES — bit-exactly, and verified in-process on
every run.** `core.predict.Distogram` builds `risk[p,t] = w[p]·Σ_c prob[p,c]·|grid[t] − CENTRES[c]|`
with `w[p] = shell/(sd+0.5)^g`, so the functional is bilinear in `prob` with a weight that is itself
a function of `prob` through `sd`. Rather than re-implement it — the silent approximation I was
warned against — the module **constructs a genuine `Distogram` object from the modified `prob` and
uses its own `_risk`**. Reconstruction of the shipped artefact from the shipped `prob` is verified
at **max abs error 0.0**. `_verify_functional()` raises and stops the run if that ever changes.

**(2) Metric realisability — so distance-space interpolation was NOT used.** `Dhat(γ)` in distance
space can be a matrix no structure attains. The ladder therefore interpolates **probability mass**,
whose γ=1 endpoint is the native's own bin and is realisable by construction. **The coordinator's
literal "substitute Dhat" arm is not expressible in the shipped functional at all, because that
functional consumes a distribution, not a point estimate** — writing it as a point-estimate loss
would be a different functional and would make the ladder uninterpretable. Three proper-distribution
arms are used instead, separating the two things "a better prior" can mean: `MASS` (mix toward the
true bin: recentre **and** sharpen), `TILT` (exponential tilt, minimum-KL: recentre **only**),
`DIRAC` (collapse onto the interpolated mean: recentre **and** fully collapse).

## The ladder

    gamma      MASS      TILT     DIRAC   MASSFIXW    beta      top-75 overlap
    0.0      3.0483    3.0483    3.0752    3.0483    +0.5413    1.000
    0.1      2.8334    2.8575    2.9254    2.9261    +0.4855    0.862
    0.2      2.6827    2.7216    2.7992    2.7762    +0.4479    0.769
    0.3      2.5847    2.5871    2.6619    2.6386    +0.4182    0.693
    0.5      2.4199    2.4149    2.4404    2.4235    +0.3711    0.572
    0.7      2.3284    2.3201    2.2799    2.2944    +0.3365    0.493
    1.0      2.2261    2.2261    2.2261    2.2332    +0.2725    0.365
    EXACT (unbinned, different functional, endpoint reference)  2.2367

Every step is significant: γ=0.1 is −0.2150 at 2.45× MDE, 113W/13L; γ=1.0 is −0.8223 at 3.22× MDE,
114W/12L. Fold CIs exclude zero throughout.

**TWO REPRODUCTION GATES, one passed exactly and one that required investigation.**

* **PASSED, to four decimals.** γ=0 through this module's own scoring path emits **3.0483**,
  matching the pinned incumbent and the module's own independently computed `P0`. And β at γ=0 is
  **+0.5413**, reproducing L5's posterior-mean β of **+0.5413 exactly** — an independent code path
  landing on L5's secondary number to four decimals.
* **INVESTIGATED, not explained away.** γ=1 lands at 2.2261 (EXACT 2.2367), about 0.25 Å above the
  standing "perfect distance knowledge caps the instrument at ~1.95–2.0 Å". Per instruction I
  treated the ladder as the suspect. **B-3b (`poolcheck.py`) tests it: same perfect-knowledge
  functional, same top-75, same readout, only the CANDIDATE SET changes.**

      restricted to the shipped K=500 BLOSUM pool   2.2367   <- the ladder's regime
      selecting from the FULL window universe       1.6024   (mean 18,674 windows/target)
      universe - pool500  -0.6342  SE 0.0572  MDE 0.1601  3.96x  fold[-0.696,-0.574]  116W/10L

  **The standing figure lies inside the bracket these two arms define.** The ladder's endpoint is
  pool-restricted *by construction* — which is the correct choice for a transfer function in the
  prior — and the restriction is worth 0.634 Å. The ladder is not wrong; it answers "what is a
  better prior worth **at fixed candidate pool**", which is the decision-relevant question.

## THE SHAPE, WHICH IS THE ACTUAL DELIVERABLE: STRONGLY SATURATING

    dRMSD/dgamma   g=0->0.1  -0.2150      curvature +0.2173  (>0 = SATURATING)
                   g=0->0.3  -0.4636
                   g=0->1.0  -0.8223

**The first 10% of the way to a perfect prior buys 26% of the total available gain; the first 30%
buys 56%.** The curve is concave, not threshold-shaped. This is the opposite of the falsifier and it
is the answer the sprint needs: **incremental prior improvement pays, and pays most at the start.**
A prior 10% better than the shipped one is worth −0.215 Å — over four times the −0.049 Å the sprint
needs to reach 3.0, and 2.45× its own MDE.

## β TRACKS γ MONOTONICALLY — L5's MECHANISM CONFIRMED AT EVERY POINT, NOT ONE

β falls 0.5413 → 0.2725 as the prior improves, monotone across all eleven rungs. L5 measured the
inheritance coefficient at a single point against a placebo floor of 0.3521; the ladder shows it is
a genuine dose-response in the prior's own error. At γ=1, β = 0.2725 sits **below** L5's
shared-referent placebo floor, which is exactly where it should be when there is no prior error left
to inherit.

## A MECHANISM WORTH KEEPING: POSTERIOR COLLAPSE FLIPS SIGN AT γ ≈ 0.55

`DIRAC` is **worse** than `MASS` while the prior is wrong (+0.027 Å at γ=0, +0.092 at γ=0.1) and
**better** once the prior is nearly right (−0.049 at γ=0.7, −0.036 at γ=0.9). Collapsing a
distribution onto its mean is harmful exactly when that mean is wrong and helpful when it is right —
a clean crossover, and a sharpened version of project memory's "do not collapse the posterior" and
"confidently wrong costs 2–3× absent". The practical reading: **confidence should be earned before
it is expressed**, and the shipped prior has not earned it.

`MASSFIXW` tracks `MASS` within 0.01–0.09 Å, so the self-consistent reweighting of `w` is a real but
second-order channel — reported rather than silently folded in.

## THE SYNTHESIS THE SPRINT NOW HAS

Two levers, both measured on the same instrument, same readout, ORACLE on both sides:

    at fixed candidate pool, a PERFECT prior is worth        -0.822 A   (B-3)
    at a perfect prior, the FULL candidate universe is worth -0.634 A   (B-3b)

They are comparable in size — but they are not comparable in *reachability*, and that is the point.
The sprint has measured the candidate lever as **flat for every achievable source** (L2d union
+0.0022; L3's quality-matched retrieval-free source parallel at cosine 0.9432; my own L7). The prior
lever is not flat, is not threshold-shaped, and pays most in its first increment. **On this
evidence the remaining leverage is in the distogram, not in the candidates** — and unlike the
candidate lever, it does not require a large jump to pay.
