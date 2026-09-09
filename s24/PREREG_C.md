# PREREG_C — WORKSTREAM C, THE FROM-SCRATCH GENERATOR

Written before `s24/c_ladder.py` runs at n=126. Exploratory sizing runs (`c_probe.py`,
`c_probe2.py`, `c_probe3.py`, `c_probe4.py`, all n=26 on the fixed index-stride subset `tg[::5]`)
are labelled EXPLORATORY in their own artefacts and are **not** promoted; this file registers the
confirmatory run and the readings that will be made of it.

---

## 0. THE PRIOR RESULT I AM REQUIRED TO CONFRONT, AND WHAT I EXPECT TO BE DIFFERENT

`phi-carries-no-sequence-signal`: a model seeing the entire 15-residue sequence context predicts
φ at **36.1° MAE against 36.4° for a distribution that sees no sequence at all** — a 0.8%
improvement. The whole measurable sequence→torsion channel is **10.4° of ψ**. Direct
sequence→torsion→build emits **4.151 Å**, worse than the 3.213 Å built-chain retrieval incumbent,
and a *perfect* confidence gate recovers only 0.5 Å of a 0.94 Å deficit.

**What is different about my object, stated up front.** The earlier result was a deterministic
point predictor scored by MAE. I am building a conditional *distribution*, sampled many times,
with a structural prior in the conditioning set, and read out through a top-75 coordinate average
rather than by building one chain. Three things could in principle differ:

1. ψ's error distribution is measured to be genuinely **bimodal**, so a point estimate cannot
   express what the model knows; a mixture can. This is the only mechanism by which the same
   information content could yield a better structure.
2. The readout is an *average of 75*, not an argmin, so per-sample angular error partially
   cancels; the earlier 4.151 Å was a single built chain.
3. The distogram enters **generation**, not only scoring (directive §26), so the sampler is not
   restricted to the sequence channel the φ result prices.

**What I expect, honestly, before running.** I expect (1) and (2) to be real and (3) to be the
thing that fails. My prior is that the lane clears the spec's RMSD half and fails to move the
endpoint, because s24 L2(d) and s23 L9 both say the binding quantity is common-mode error and
neither (1) nor (2) touches it. I am pre-committing to that expectation so that a null is not
retrofitted as a prediction and a win is not retrofitted as a surprise.

**The experiment that would show this path cannot clear the spec** is §2 arm T2 plus the mixture
ladder: a sequence-conditioned torsion distribution that clears ≤3.9 Å standalone and still leaves
every matched-size mixture at or above the pure incumbent. That is a falsifier for the lane, not
for the model, and it is registered as the primary outcome to look for.

---

## 1. THE LADDER — ordered by how much information the torsion distribution may see

Every arm emits N built chains per target, is scored by the **same** shipped Bayes-risk distogram
functional, and passes through the **same** readout. Only the torsion distribution changes.

| arm | conditioning | status |
|---|---|---|
| `R_rebuild75` | the incumbent's own top-75, rebuilt from their own torsions | CALIBRATION, not a generator |
| `T0_helix` | nothing (constant α-helix −63/−42 + 15° jitter) | ZERO-INFORMATION CONTROL, plausible |
| `T1_blind` | sequence-blind (φ,ψ) *pairs* from this target's legal universe | ACHIEVABLE |
| `T2_restype` | per-residue-type Ramachandran from the same universe | ACHIEVABLE — the sequence channel at its measured ceiling |
| `T3_pool` | per-residue 2-component von Mises fit to the retrieved top-75 torsions | ACHIEVABLE UPPER REFERENCE |

`T3_pool` is the reference a trained model must beat: s14 measured this exact channel at
**φ 33.6° / ψ 59.2°**, better than the trained leave-fold-out sequence predictor's 36.1 / 62.4,
at zero training cost. A learned `p(φ,ψ | sequence, distogram)` that does not beat T3 has not
earned its parameters.

**No model is trained in this run.** Nothing here needs Lane A's corpus: the per-target universe
in `s8/generate_univ/<pdb>.npz` is the project's already-audited leakage-safe library
(out-of-fold peptides + this fold's fragments), and no arm reads `rr` or `nat_ca` except for
post-hoc ORACLE scoring. Lane A's corpus hash will be required before any *trained* arm, and is
cited as absent here rather than assumed.

## 2. THE SIX OPERATOR FORKS (Rule 0). Each names the alternative not taken.

- **functional** — DECLARED the shipped Bayes-risk distogram score, unchanged, for every arm
  including the incumbent. NOT TAKEN a re-tuned score, a temperature, or a generator-specific
  score, either of which would confound source quality with selector change.
- **basis** — DECLARED point cloud for every emitted average, on both sides of every contrast:
  members are built chains, but a uniform average of 75 chains is not itself an ideal-geometry
  chain, so the emitted object is a point cloud and is commensurable with the incumbent's
  3.0483 Å. The built-chain figure (projection onto ideal geometry, comparable to 3.204 Å) is
  reported **separately and never mixed**. NOT TAKEN quoting the built-chain number against the
  point-cloud incumbent, which is the 0.156 Å operator gap the brief forbids.
- **readout** — DECLARED uniform mean of exactly 75 members, selected as top-75 by the shipped
  score, matched in count across every arm. A second, clearly separated **uniform-75 (unscored)**
  readout is reported alongside, because L2's published 0.647 cosine was measured on an *unscored*
  library draw and is not commensurable with a scored one. NOT TAKEN a single readout, which
  would silently pick one of two non-comparable definitions of the spec.
- **normalisation** — DECLARED cosine of raw bias vectors in the native frame, plus the ratio to
  a within-source control built from an independent draw of the *same* sampler. NOT TAKEN the raw
  cosine alone, which is uninterpretable for structures that are all protein-like.
- **null** — DECLARED (i) the within-source cosine as the "shares all its bias" ceiling; (ii)
  `T0_helix`, a plausible zero-information source, not uniform-on-the-torus; (iii) for every
  best-of-N quantity, a matched-count best-of-500 against the shipped K=500 pool, since the
  generator's larger sample budget confers a maximum-order advantage that is not skill.
  NOT TAKEN zero cosine, uniform torsions, or an unmatched sample count.
- **THE LABEL** — DECLARED continuous Cα-RMSD in Å for every endpoint, and the raw signed cosine
  for every angle. NOT TAKEN any binarised "did it decorrelate", any rank statistic, and no
  thresholded pass/fail on the spec beyond quoting the two numbers.

## 3. PRIMARY OUTCOMES, in the order they will be read

1. **standalone mean Cα-RMSD** of each arm, n=126, point cloud, scored top-75 uniform average.
   Spec: ≤ ~3.9 Å to break even.
2. **bias cosine** vs the incumbent, measured exactly as `s24/biasalign.py` measures it (native
   frame, `_bias`/`_cos` reused verbatim), with the within-source control and the ratio.
   Spec: ≤ ~0.65 to break even.
3. **the endpoint** — matched-size mixture ladder (75 members always: 75/0, 60/15, 50/25, 38/37,
   25/50, 0/75) and the §17 union at matched count. Departure from the straight line between the
   endpoints is the geometry-free independence check.
4. **coverage (coordinator's route (a))** — score quantiles of generated vs retrieved at matched
   count, union take-up share, and ORACLE best of each set. Route (a) requires *both* that the
   score prefers generated chains *and* that they are structurally better.
5. **mode-collapse audit** and **geometric validity**, both reported before any RMSD is promoted.

## 4. FALSIFIERS, written before the run

- **Lane falsifier (the one I expect to fire).** Every matched-size mixture at or above the pure
  incumbent AND the union within its own MDE, while the standalone RMSD clears 3.9 Å. That says
  the spec is necessary but not sufficient and the lane does not pay.
- **Mechanism falsifier for L3.** If `cos(emitted distance error, distogram distance error)` does
  NOT rise when the selector is applied — i.e. scored ≈ unscored across arms — then selection is
  not what makes sets parallel and L3's mechanism is wrong.
- **Route (a) falsifier.** If generated samples both take a large share of the union top-75 and
  have a *worse* ORACLE best than the matched-count pool, route (a) is refuted on its own terms.
- **Model-class falsifier.** If `T3_pool` (zero training, 3 parameters per residue per component)
  is not beaten by any trained model on the standalone number, the lane does not justify a network.

## 5. STOPPING RULE

If the lane falsifier fires at n=126, I report it, do not tune, do not search for a variant that
survives, and say so to the coordinator so the sprint's compute goes to Lane B. A negative result
delivered early is the deliverable.
