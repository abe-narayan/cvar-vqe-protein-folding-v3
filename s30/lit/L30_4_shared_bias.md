# S30 LANE L, TOPIC 5 -- BREAKING A SHARED BIAS GIVEN A PRIOR ON ITS FORM

Commissioned by the coordinator on spawning lane P, whose question is:

> Can the distance prior's systematic **shape** error be predicted and corrected using information
> demonstrably not already in it?

Four sub-items were set: (1) sources provably decorrelated from a sequence-model distance
prediction, (2) shape-versus-scale decompositions of predicted distance matrices, (3) systematic
bias correction where the bias is unidentifiable from the corrupted data alone, (4) a routing
decision on the AMBER-relax split. Ledger: **S30-L18**.

**Method note, and it is the reason item (2) came out the way it did: I read the project's own
record before searching outside it.** That is what caught the S19 prohibition below, and it is the
third time this sprint that an existence-or-precedent check killed a proposal of mine that had
already passed the information test.

---

## 1. ITEM 2 FIRST, BECAUSE IT IS A RE-IMPORT WARNING

**What I had.** The distogram emits per-pair marginals, and the project's DIS score consumes them
pair-independently. A symmetric hollow matrix on `n` residues has `n(n-1)/2` free entries; a
realisable 3-D structure has `3n - 6`. At our lengths:

```
    n =  9    36 entries   21 DOF   ->  42% of the prediction's freedom is metric inconsistency
    n = 13    78 entries   33 DOF   ->  58%
    n = 16   120 entries   42 DOF   ->  65%
```

Schoenberg's theorem gives the exact test (`D` is an EDM iff `-1/2 J D^2 J` is PSD) and the rank-3
condition adds the embedding constraint. Projecting onto that cone is parameter-free, uses no new
information, and would be an "E2 on the prior" in S30-L7's taxonomy.

**Why it is dead.** `s19/agentA_FINDINGS.md` section 1.2, with `s19/a_topo.py` and `s19/a_coh.py`:

| quantity | predicted field | native |
|---|---|---|
| triangle-inequality violations | 4.09 % | 0.02 % |
| rank-3 EDM defect of `-1/2 J D^2 J` | **0.286** | 0.0018 |
| non-Euclidean mass | 0.183 | 0.0009 |

> *"**Independent pairwise marginals are not jointly realisable** is real, measured, and **NOT the
> mechanism**. Nobody should spend on EDM projection, triangle repair, or joint-consistency
> enforcement as a route to RMSD."*

The A3/A4 crux was decided **against** unrealisability: a magnitude-matched *incoherent* field is
worse on both measures (0.340 defect, 10.85% violations) and lands **1.24 A better**. And
`rho(EDM defect, RMSD) = +0.511` against `rho(residual RMS, RMSD) = +0.891` -- the defect is a
**magnitude proxy** and adds nothing over error size. My degrees-of-freedom arithmetic is correct
and irrelevant.

**REJECTED, with the prohibition quoted.** Sources engaged and not imported: Schoenberg (1935);
Gower's EDM rank results; Dattorro's EDM-cone treatise; the rank-constrained EDM least-squares
model for protein conformation (J Glob Optim 2019); low-rank EDM completion (arXiv:1804.04310);
alternating/Dykstra projection onto the EDM cone. **All describe machinery for a defect that is
measured here to be a symptom rather than a cause.**

## 2. WHAT THE RECORD ALREADY KNOWS ABOUT LANE P's TARGET

From the same S19 file, sections 1, 4, 5 and 6 (ORACLE, n = 126):

- **The harm is COHERENCE.** A perfectly coherent matched-magnitude error is the worst arm on the
  board (3.749 vs `real` 3.610). Destroying only cross-pair **sign** coherence at exactly matched
  magnitude buys **-1.202 A [-1.408, -1.004], 112W/14L, 5/5 folds**.
- **Making the field more realisable makes it more coherently wrong** (EDM defect 0.286 -> 0.148,
  kappa -> 0.914, and the arm is worse). S19's own section 2.
- **The mode is named: a per-target separation profile, five numbers per target**, owning
  **+52.5% of the gap, worth -0.525 A** with the native in hand.
- **Variance explained does not price damage.** Nested variance: offset (1) 0.141, stretch (1)
  0.226, separation profile (5) 0.361, per-residue additive (13) 0.430 -- yet removing the offset
  *hurts* (+0.174), removing the stretch *hurts* (+0.069), and the per-residue additive explains
  the most and owns **none** of the harm.
- **The wall.** The best native-free estimator recovers **0.098-0.155 A of the 0.525 A (19-30%)**
  and only by regressing toward what the pipeline already does. *"There is still no native-free
  estimator of coherence."*
- **The bias is SHARED across independent sources** (section 4); the incoherent component does not
  share at all (0.09-0.19 cross-family).

**That last point is S30-L7's escape E3 measured two sprints before I derived it**, and S24's
provenance cosine (0.9432 against a 0.9330 within-source control) is a third independent arrival.

**Composition with S29-L31.** Five numbers per target are five **incidental parameters**. Neyman &
Scott name exactly two escapes: replication within the instance, or a covariate observed at
inference. So lane P's only viable form is:

> a native-free covariate, **observed at inference**, that predicts the per-target separation
> profile **and is generated by neither the distogram nor the retrieval pool** -- because both are
> measured to share the bias being estimated.

**One load-bearing caveat I could not discharge.** S19 predates the benchmark pinning recorded in
memory `benchmark-and-folds-must-be-pinned` (correcting the identity clustering silently moved 13
targets and invalidated every fold model). **The 0.525 A should be re-checked on the pinned
benchmark before anything is built on it.** I flag this rather than assume it; it is the single
most important unverified number in this file.

## 3. ITEM 1 -- WHAT DECORRELATION IS WORTH, DERIVED

Full derivation and table: `s30/lit/s30_L_orthogonality.py`, ledger S30-L18 section 4.

With `rho_max^2 = (r1^2 - 2 c r1 r2 + r2^2)/(1 - c^2)` and the project's anchors (best single
field 0.1128, 0.358 needed for 3.00 A):

- a **perfectly orthogonal** new channel must itself carry `rho = 0.3398`, against `0.3580` if it
  worked alone: **orthogonality is a 5.1% discount on the requirement** (1.6% for 2.50 A);
- the new channel must be **3.01x better than anything the project owns**;
- the gain is **quadratic in the new channel's own skill**: a perfectly orthogonal channel as good
  as everything we own buys **-0.0206 A**;
- reaching 3.00 A by stacking needs **k = 10.1 mutually orthogonal channels** each as good as our
  best; lane F measured the 21 we have at Gram stable rank **2.057**.

> **Decorrelation is not the lever. Skill is.**

**SCOPE LIMIT, and it matters more than the result.** This is the S29 section 5.1 displacement
law, which governs operators in **pool space**. **Lane P works on the PRIOR**, where the
derivative is **-2.15 A per unit** (`prior-derivative-is-the-only-steep-lever`) and the ORACLE
prize is 0.525 A on one named mode. Section 3 prices *"add another field to the 21"*. It does
**not** price *"fix the distance prior's shape"*. The flat lever must not be quoted at the steep
one.

## 4. ITEM 3 -- THE MEASUREMENT-ERROR LITERATURES, AND THE ONE CONDITION THEY ALL NEED

| family | what it needs | why it fails here | verdict |
|---|---|---|---|
| **errors-in-variables** with two error-laden measurements (reliability ratio from their covariance) | the two measurements' errors are **INDEPENDENT** | our two are the distogram profile and the pool profile, and S19 section 4 **measured that they share the bias**. Not unverified -- measured false | REJECTED |
| **Reiersol, Econometrica 18:375-389 (1950)**: EIV identified without an instrument when the latent regressor is **non-normal** (higher moments do the work) | a linear relation between **two observed** error-laden variables | we have one prediction and no second observable of the same quantity in that relation | REJECTED, recorded because it is the standard answer to "identify without an instrument" and will be proposed |
| **multichannel blind deconvolution** (Xu et al. 1995 and successors) | the channels are **COPRIME** -- no common zeros | a shared bias *is* a common factor; coprimeness is exactly what S19 section 4 and S24's 0.9432 measure to be absent | REJECTED |
| **single-channel blind deconvolution** | sparsity / non-negativity / known support | none holds for a distance-profile bias | REJECTED |
| **instrument calibration** | a reference standard, i.e. an anchor | the anchor is the native | REJECTED (= S30-L7 E1/E3; = S29-L4's control-variate rejection) |
| **Kennedy & O'Hagan (2001); Brynjarsdottir & O'Hagan (2014)** | an informative prior on the discrepancy's **shape** | KEPT from S30-L7 -- and S19 section 5 *supplies* the shape (a per-target separation profile), which is why this route is the live one | KEPT |

> **One condition under five names:** every method that separates a shared bias from the truth
> needs **two views whose errors are independent**, and this instrument's two views are measured to
> share the bias. That is S30-L7's E3, reached from the measurement-error literature instead of
> from the likelihood.

**The one constructive reading.** B&O say a discrepancy term helps only given an informative prior
on its *shape*. S19 **has** the shape -- five numbers, a separation profile. So this project is,
unusually, in the regime B&O describe as workable: the form is known, and what is missing is the
per-target *values*. That is a genuinely better position than "we have a bias and no idea of its
structure", and it is why I think lane P is staffed on the right question.

## 5. ITEM 4 -- MY ROUTING RECOMMENDATION, AND A PRIOR I LOWERED

**Hand the k=30 AMBER-relax split (by pool dispersion, and by FAIL18/108) to whoever owns the
artefact, not to me.** It is zero compute and it tests a claim the averaging-artifacts paper
asserts and never measures.

**But I have lowered my own prediction.** S30-L13 put it at 2:1 that the benefit concentrates on
high-dispersion targets. Lane F then measured that pool Rg dispersion does **not** predict filter
failure (-0.128, CI includes zero). Different outcome variable, so not a refutation -- but the
same family, pointing the other way. **Roughly even now**, and not a bet I would defend.

## 6. WHAT DAMAGED MY EXPECTATIONS IN THIS TOPIC

- I expected the EDM-consistency route to be open and possibly strong. It is explicitly
  prohibited in the record, and the reason is the opposite of what I assumed: unrealisability is
  **anti-correlated** with harm.
- I expected "find a decorrelated source" to be a reasonable brief. Priced, orthogonality is worth
  5% and the whole requirement is skill.
- I expected the measurement-error literature to offer at least one route that does not need an
  anchor. Reiersol's non-normality result genuinely is such a route **in its own setting**, and it
  does not transfer, because we lack the second observed variable it operates on.
- I expected to hand lane P a method. I am handing it a **narrowed target and three prohibitions**,
  which is less satisfying and, on this project's record, worth more.
