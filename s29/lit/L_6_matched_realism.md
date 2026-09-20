# L_6 -- RANKING WITHIN A MATCHED-REALISM BAND: WHAT THE LITERATURE OFFERS LANE D (brief topic 6)

Lane L, 2026-09-20. Commissioned by the coordinator after S29-L12. The idea under test: the
perception-distortion theorem (Blau & Michaeli, CVPR 2018) forbids a realism measure from
preferring the distortion-optimal answer ACROSS realism levels, but says nothing about ordering
INSIDE one level. So: fix a realism band, then ask whether anything orders candidates by accuracy
within it.

Four questions were asked. The answers, stated first.

(a) YES, the perception-distortion literature bands -- and its own banded result is a WARNING.
    The PIRM 2018 challenge is exactly this construction (the transpose: fix distortion, rank by
    realism), and it found its realism index correlates with human judgement at Spearman 0.83
    BETWEEN bands but "is not always well-correlated with these scores on a finer scale (rankings
    within the regions)". The skill largely died inside the band. That is the project's own
    global-vs-in-band split, measured independently in another field, on this exact design.
(b) PARTLY. Statistical physics gives the construction and, more usefully, a theorem about what
    banding destroys: inside an exactly constant-energy shell the energy carries ZERO ordering
    information by construction, so every in-band ordering comes from something orthogonal to the
    band variable. Umbrella sampling/WHAM supplies the machinery for undoing a window selection.
    The QA literature has band-like constructions (CASP's best150 / sel20 / BZQ15 subsets) but
    they were built for other reasons.
(c) NO -- and this is the finding the coordinator asked me to state explicitly. I found NO
    published QA evaluation that conditions on a realism statistic and then measures the
    correlation of a score with GDT/RMSD inside the band. The paper that comes closest is
    Hamelryck et al.'s artefacts-and-biases study, which documents the confounds in detail and
    explicitly does NOT do the matching. Lane D's measurement would be novel, not derivative.
(d) YES, there is an exact analogue of the Ueda-Nakano coefficient, and it is the PARTIAL
    CORRELATION. For jointly Gaussian (S, Y, R), the correlation of scorer S with accuracy Y
    conditional on realism R equals the partial correlation rho_SY.R and is CONSTANT in the
    conditioning value. Band width interpolates: as the band narrows, in-band skill goes to
    rho_SY.R; as it widens, it returns to the global rho_SY. The bound is sharp and it has a
    zero: if a scorer's association with accuracy is entirely MEDIATED by realism
    (rho_SY = rho_SR rho_RY), its in-band skill is exactly zero.

---

## (a) Banded evaluation in the perception-distortion literature

### Blau Y, Mechrez R, Timofte R, Michaeli T, Zelnik-Manor L. "The 2018 PIRM Challenge on Perceptual Image Super-resolution." ECCV 2018 Workshops; arXiv:1809.07517

Construction, verbatim. The motivation for banding is that methods at different distortion
levels cannot be compared at all: perceptual algorithms "cannot participate in any challenge or
benchmark based on these standard measures ... and cannot be compared or ranked using these
common metrics". So "the perception-distortion plane was divided into three regions by setting
thresholds on the RMSE values (regions 1/2/3 were defined by RMSE <= 11.5/12.5/16 respectively).
In each region, the goal was to obtain the best mean perceptual quality." The realism axis is
PI = (1/2)((10 - Ma) + NIQE), both no-reference measures -- i.e. native-free, which is the
property we need. Validation: 35 human raters scored outputs 1-4 for "how realistic the image
looked", explicitly without seeing the ground truth, so "this study does not test distortion in
any way, but rather only perceptual quality".

THE RESULT THAT MATTERS FOR LANE D, verbatim: "while the PI is well correlated with the
human-opinion-scores on a coarse scale (in between regions), it is not always well-correlated
with these scores on a finer scale (rankings within the regions)". Spearman 0.83 across bands;
unreliable within them.

Also verbatim, and it is our contracted average in another field: "the outputs of EDSR, a
state-of-the-art algorithm in terms of distortion, are mostly voted as 'definitely fake'. This is
due to the aggressive averaging causing blurriness as a consequence of optimizing for distortion."

Check against this instrument. The design transfers directly and the warning transfers with it.
The project has the same split already measured on its own axis: the distogram is +0.653 global
and +0.091 in-band (`in-band-is-the-only-ranking-metric`), and S28-L48's whole recognition audit
is a between-band comparison (production average vs a 0.25 A ORACLE structure). PIRM says the
expected outcome of banding is that the metric loses most of its apparent skill -- which is not a
reason to skip the experiment, because the residual is exactly the quantity nobody here has
measured, but it IS the reason to power the experiment for a small effect and to pre-register the
band.
KEPT as the design template and as the prior for the effect size (small).

### The tradeoff's own shape, from Blau & Michaeli and confirmed in PIRM

PIRM observed "the tradeoff appears to be stronger in the low distortion regime (Region 1)",
which is Theorem 3's convexity measured. Our operating point is the low-distortion end, so the
project sits where the tradeoff is steepest -- the band must therefore be narrow in realism to be
meaningful, and narrow bands cost sample size. That tension is the experiment's main design
problem, and (d) below prices it.

---

## (b) Equal-energy shells, umbrella sampling, and matched decoys

### The microcanonical statement, which is the cleanest thing in topic 6

Inside a shell of exactly constant energy E, the Boltzmann weight e^{-beta E} is constant, so the
energy induces NO ordering within the shell; all structure inside a constant-energy shell is
carried by the density of states (the entropy). Generalising beyond energy: banding on a
statistic R removes R's own discriminating power by construction, so an in-band ranking
experiment measures ONLY what is orthogonal to R.

This is why the experiment is worth running and also what it can and cannot conclude. A positive
in-band result is, by construction, evidence of information orthogonal to realism -- which is
precisely the quantity H0 in `s29/STATE.md` says the system lacks. A null is evidence that the
library's scorers carry nothing beyond realism. Either way the measurement is interpretable,
which is the property the coordinator wants.

### Umbrella sampling (Torrie & Valleau 1977) and WHAM (Kumar et al. 1992)

Construction. Add a bias potential w(xi) to restrict sampling to a window of a reaction
coordinate xi, sample within the window, then UNBIAS: the unbiased average is recovered as
<A> = <A e^{+beta w}>_biased / <e^{+beta w}>_biased, and WHAM combines overlapping windows by
optimally weighting them.

Assumption and transfer. The machinery assumes you know the bias exactly (you imposed it), so you
can undo it. Lane D's band is an imposed selection with a known rule, so the same logic applies
and gives a concrete discipline: an in-band result is a statement about the band's population,
and any claim about the FULL pool requires reweighting by the (known) selection. Practically this
means D should either (i) confine every claim to the band, or (ii) report the band-conditional
effect and the selection fraction separately. It also warns against the obvious mistake: the
band's members are not a random sample of the pool, so a within-band mean is not an estimate of a
pool mean.
KEPT as discipline, not as method (there is no free energy to reweight here).

### "Energy-matched" / matched decoy constructions in the QA literature

What exists: CASP11's evaluation subsets, used by VoroMQA and others -- best150 (the 150 models
per target selected by a consensus QA algorithm, i.e. a HIGH-QUALITY band), sel20 (20 models
selected by clustering to be as DIFFERENT from each other as possible, i.e. a diversity-spread
set), and VoroMQA's own BZQ15 (up to 15 models from three strong servers, simulating a realistic
selection task). These are restricted sets, so the field does evaluate on subsets.

The result, and it is a methodological warning rather than a template, verbatim from VoroMQA
(Olechnovic & Venclovas 2017): "for sel20 sets, every method that is based solely on analyzing
geometric features and applying statistical potentials (GOAP, DOOP, dDFIRE and all the VoroMQA
variations) achieved worse results than the best-performing composite methods incorporating
evolutionary information in the form of predicted features such as secondary structure or solvent
accessibility: this was definitely not the case for best150 and BZQ15 sets." And, on the same
sets, a trivial "HHpred-agreement" score (TM-score to a homology model, i.e. agreement with an
INDEPENDENT predictor) "performed very similarly as ProQ2-refine and ProQ2 for sel20, but much
worse than all the other tested QA scores for best150 and BZQ15".

Two lessons for lane D, both important. (1) WHICH METHOD WINS DEPENDS ON HOW THE SUBSET WAS
BUILT -- pure-geometry (realism) scorers and information-carrying scorers swap places between
subsets. So the band's construction is not a neutral choice and must be pre-registered with its
rationale, and a result on one band definition does not transfer to another. (2) The thing that
rescued the diverse subset was agreement with an independent predictor, which is the
quasi-single-model signal from L_1 section 2.2 that this instrument does not have. Neither subset
is a realism-matched band: best150 is quality-matched (by consensus, which is the very thing our
pool violates) and sel20 is diversity-selected.
NOTED; no importable construction.

### ANDIS (Yu et al., Bioinformatics 2019) -- a protein-specific echo of the theorem

Recorded because it is the same tension in our own field: the paper's framing is that "native
recognition emphasizes the differences of overall structure quality between native and decoy
structures, while decoy discrimination generally focuses on the backbone differences among decoy
structures ... The potential's abilities of native recognition and decoy discrimination cannot be
optimized simultaneously with the same parameter sets." Between-band ability and within-band
ability require different parameters, and trade against each other. This is the protein-potential
community reaching Blau & Michaeli's conclusion empirically, without the theorem.
KEPT as corroboration, and as a concrete prediction: a scorer tuned for our between-band task
(which is how every statistical potential was fitted) is the WRONG scorer for the in-band task,
so a null from the existing library would be expected and would NOT close the in-band question.

---

## (c) Has any QA evaluation conditioned on a realism statistic first? No.

### Hamelryck T et al. "Artefacts and biases affecting the evaluation of scoring functions on decoy sets for protein structure prediction." (PMC2677743)

What it establishes. Decoy-set evaluation is confounded several ways: for "139 out of 149 of the
decoy sets considered" the native is TRIVIALLY discriminable via individual energy terms
(improper torsions, vdW clashes) rather than by fold quality; MD-generated decoys violate the
i.i.d. assumption, with correlations "aris[ing] primarily as a consequence of differences
_between_ the five trajectories" rather than within them; and enrichment with near-native
structures inflates apparent performance.

What it does NOT do, and this is the answer to (c). It does not propose or implement matching or
conditioning on a confounder -- no compactness-matched comparison, no analysis within a narrow
quality range. Its recommended mitigations are about effective sample size, not about
conditioning. So the most careful methodological paper in this corner of the field identified the
confounding problem and stopped short of the fix lane D is about to apply.

Combined with the CASP subsets above (which are quality- or diversity-selected, not
realism-matched) and with the CASP14/15 assessments (which report per-target correlations and
top-1 losses over full model sets), my conclusion is:

**Nobody has conditioned on a native-free realism statistic and then measured within-band
accuracy ordering.** Lane D's measurement is novel. The nearest published relatives are PIRM's
distortion-banded ranking (the transpose, in vision) and the CASP subsets (restricted, but not on
realism). I flag the usual caveat on any negative literature claim: this is the result of a
targeted search across the QA, decoy-evaluation and perception-distortion literatures, not a
proof of absence.

---

## (d) How much signal can survive inside a band: the exact analogue of Ueda-Nakano

Let Y be accuracy (CA-RMSD), R the native-free realism statistic used to define the band, and S a
candidate scorer. Two standard results give the bound.

**1. The thin-band limit is the partial correlation.** For jointly Gaussian (S, Y, R), the
conditional covariance of (S, Y) given R = r is Sigma_{SY} - Sigma_{SR} Sigma_{RR}^{-1}
Sigma_{RY}, which does not depend on r. The conditional correlation is therefore the partial
correlation, constant across the band:

    rho_SY.R = ( rho_SY - rho_SR rho_RY ) / sqrt( (1 - rho_SR^2) (1 - rho_RY^2) )

**2. The zero of that expression is the whole point of the experiment.** rho_SY.R = 0 exactly
when rho_SY = rho_SR rho_RY, i.e. when the scorer's entire association with accuracy is MEDIATED
by realism. So the matched-realism band is precisely a test of whether any scorer in the library
carries information about accuracy beyond its realism content. This is the sharp form of the
coordinator's question, and it is directly measurable: D needs only rho_SY, rho_SR and rho_RY,
all in-band-free quantities already computable from the S27 channel cache plus one choice of R.

**3. Band width interpolates, and that is the design knob.** Conditioning on a BAND R in [a, b]
is not conditioning on R = r. For a wide band the within-band correlation retains part of the
mediated signal and approaches the global rho_SY; as the band narrows it converges to rho_SY.R.
So in-band skill as a function of band width is a monotone curve between the global correlation
and the partial correlation, and its limit is the quantity of interest. RECOMMENDATION to lane D:
report skill at several band widths and extrapolate, rather than choosing one width -- this also
defuses the obvious multiplicity objection (one pre-registered curve, not k pre-registered bands)
and it converts the sample-size/width tension into a measured trend instead of a single
underpowered point.

**4. The classical name for the attenuation is range restriction** (Thorndike's case II
correction in psychometrics): selecting on a variable correlated with both S and Y predictably
attenuates the observed S-Y correlation, and the correction requires knowing the unrestricted
variance of the selection variable -- which D has, because the band is imposed. This is the
formal reason a raw within-band correlation must never be compared to a raw global correlation
without stating the selection.

**5. The project has already met this logic and recorded the lesson.** The memory note
`shared-referent-floor` says two quantities measured against a COMMON reference correlate by
construction and that the floor must be measured first -- it "turned '2/3 sequence-independent'
into '1/5'". A realism band is the same correction applied deliberately rather than as a
diagnostic. And `decorrelated-errors-exist-but-are-unusable` already reports TRUTH-PARTIALLED
error correlations of 0.04-0.26 with fusion worth only +0.004-0.011 "because the gain goes as
the SQUARE of the weaker channel's skill" -- so the project has measured partial quantities
before, and the squared-skill law says a small partial correlation buys very little. That is the
honest prior for the Angstrom value of any in-band positive.

**Expected magnitude, stated as a method not a number.** I will not invent a figure. But D can
predict the answer before running the endpoint: plug the measured rho_SY (a scorer's global
correlation with RMSD), rho_SR and rho_RY into the formula above. The project's published pair
for the distogram -- +0.653 global, +0.091 in-band -- is consistent with heavy mediation, and if
the same ratio holds for the realism band, the partial correlations will be small and the
squared-skill law will cap the Angstrom value near zero. A positive would therefore be most
interesting as a MECHANISM result (information orthogonal to realism exists) rather than as an
accuracy result.

---

## Verdict table (topic 6)

| # | source | what it gives lane D | verdict |
|---|---|---|---|
| a1 | Blau et al., PIRM 2018 challenge, arXiv:1809.07517 | the banded-evaluation design, with a native-free realism axis and human validation; and the warning that the realism index correlates 0.83 BETWEEN bands but unreliably WITHIN them | KEPT -- the design template and the prior (expect a small in-band effect) |
| a2 | the same, on tradeoff shape | the tradeoff is steepest at low distortion, where we operate, so the band must be narrow | KEPT as a design constraint |
| b1 | the microcanonical / constant-energy-shell argument (textbook) | banding on R destroys R's own ordering power, so an in-band result measures exactly what is orthogonal to R | KEPT -- this is what makes the experiment interpretable either way |
| b2 | Torrie & Valleau (1977); Kumar et al. WHAM (1992) | the discipline for an imposed selection: confine claims to the band, or reweight by the known selection | KEPT as discipline |
| b3 | CASP11 best150 / sel20 / BZQ15 subsets, via Olechnovic & Venclovas (2017) | restricted-set evaluation exists, and the winning method CHANGES with how the subset was built; the sel20 rescue came from agreement with an independent predictor | NOTED -- a warning that the band's construction must be pre-registered; not a realism-matched band |
| b4 | ANDIS (Yu et al., Bioinformatics 2019) | native recognition and decoy discrimination "cannot be optimized simultaneously with the same parameter sets" | KEPT as corroboration and as a prediction: library scorers were fitted for the between-band task |
| c | Hamelryck et al., PMC2677743 | the confounds are documented (139/149 decoy sets trivially discriminable; MD non-independence; near-native enrichment) and the matching fix is explicitly NOT done | KEPT -- the basis for "nobody has done this" |
| d1 | Gaussian conditional correlation = partial correlation | the thin-band limit rho_SY.R, constant across the band, with an exact zero at full mediation | KEPT -- the analogue of the Ueda-Nakano coefficient the coordinator asked for |
| d2 | Thorndike case II range restriction | the name and the correction for selection-induced attenuation | KEPT as reporting discipline |

KEPT 7, NOTED 1. No importable operator; the deliverable is a design, a bound, a prior and a
warning.

## What topic 6 tells lane D, in five lines

1. The design is published (PIRM) and validated against human judgement, so band-conditional
   evaluation is a legitimate methodology, not an invention of this sprint.
2. The measurement itself is NOVEL: no QA evaluation has conditioned on a realism statistic
   before measuring accuracy ordering. The nearest paper documents the confound and declines the
   fix.
3. The quantity being measured has a closed form: the thin-band limit is the partial correlation
   rho_SY.R, which is exactly zero iff the scorer's association with accuracy is fully mediated
   by realism. Compute rho_SY, rho_SR, rho_RY first -- the answer is predictable before any
   endpoint run, and that prediction is the pre-registration.
4. Report skill versus BAND WIDTH as a curve (global correlation at one end, partial correlation
   at the other) rather than picking one width: it answers the multiplicity objection and turns
   the width/power tension into a measured trend.
5. Expect a small effect, and expect the existing library to underperform, because every scorer
   in it was fitted for the between-band task (ANDIS's point). A null from the current library
   would NOT close the in-band question -- it would say the library is the wrong instrument for
   it, which is a different and weaker claim. Say which claim is being made before running.
