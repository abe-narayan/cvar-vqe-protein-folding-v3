# L_8 -- PER-TARGET CONDITIONING, THE SIGN PROBLEM, AND THE POOL'S OWN DISPERSION (brief topic 8)

Lane L, 2026-09-20. The last literature assignment of S29. Aimed at the sprint's unified finding:
what the system lacks is a PER-TARGET SIGN, equivalently a conditioning signal supplied at
inference. Priced exactly by lane T (S29-L23): every native-free operator is a displacement worth
one cosine, RMSD = RMSD_prod sqrt(1 - rho^2), and a per-target sign enters as |rho|(2q-1) -- so
buying the sign is worth as much as buying the direction.

The four answers, stated first.

(a) The field's per-target conditioning signal is MSA DEPTH, and it does not exist for us. Every
    instance-level signal that changes a structure predictor's behaviour per target is an
    alignment statistic (Neff, depth, coupling strength) or a template match. At 9 to 16 residues
    with no MSA and no template, all of them are identically absent. Recycling, on the coordinator's
    specific question, is OPTIMISATION, not information: it conditions on the model's OWN previous
    output and admits no new external evidence, so it cannot supply a sign not already implied by
    the inputs.
(b) THERE IS A THEOREM, AND IT IS EXACTLY OUR SITUATION. The per-target sign is an INCIDENTAL
    PARAMETER in the sense of Neyman & Scott (1948): one nuisance parameter per target with a
    bounded number of observations per target. Pooled maximum likelihood is INCONSISTENT for it,
    and the standard remedy -- conditional or fixed-effects likelihood, which is precisely the
    within-group pairwise structure of topic 7 -- ELIMINATES the nuisance rather than estimating
    it. That closes the loop: cross-target training cannot estimate the sign, and the design that
    handles the confound correctly is the design that throws the sign away. The literature names
    the only two escapes: replication within the target, or a covariate observed at inference.
(c) YES, there is a mature literature on ensemble dispersion as information, and it carries a
    quantitative warning: in numerical weather prediction the spread-skill relationship is known
    to be WEAK EVEN FOR A PERFECT ENSEMBLE, is larger only where the day-to-day variability of
    spread is large, and is most useful when the spread is EXTREME. That is a per-instance signal
    with a structurally limited correlation -- the same shape as everything else the project has
    measured, and it predicts where to look (the tails, i.e. the regime question) rather than
    promising a uniform gain.
(d) My judgement, in section (d), labelled as judgement.

---

## (a) Where the field gets per-target conditioning

### Recycling -- optimisation, not information (the coordinator's specific question)

AlphaFold2 (Jumper et al., Nature 596:583, 2021) runs 4 recycling iterations, feeding the pair
representation, the single representation and the predicted CB coordinates back into the trunk,
with the final loss applied at each iteration (L_2 section 2.7). Every quantity fed back is the
network's OWN OUTPUT. No external evidence enters at recycling time. So recycling is iterative
refinement within a fixed information set: it can improve the use of information already present
in the inputs, and it cannot create a per-target signal that the inputs do not determine.

The distinction matters for this project because it separates two things the record keeps close
together. The project's ORACLE diagnostics already say the information is not the barrier to
EXPRESSION (a 27-parameter family holds a 0.25 A structure on every target, S28-L26b) but is the
barrier to SELECTION. Recycling addresses the former class of problem. It is therefore not a
route to the per-target sign, and a "recycling-like" iterative scheme built here should be
expected to move nothing, for a stated reason rather than by analogy.
REJECTED as a conditioning source; the classification (optimisation, not information) is the
deliverable.

### MSA depth / Neff -- the field's real per-target conditioner, and it is absent here

This is where the field's per-target behaviour actually comes from. AlphaFold2's accuracy scales
with the number of effective sequences and degrades substantially below an MSA depth of roughly
30; residue-level pLDDT correlates with the number of aligned amino acids per residue; orphan
sequences and de novo designs, which lack MSA signal, are the known failure class, and ColabFold
exposes a single-sequence mode precisely because that regime is different. The field's honest
statement of "how well will this target go" is, at bottom, an alignment statistic.

Does it exist for a 9 to 16-mer? **No.** A 13-residue query has no alignable homolog set in any
meaningful sense: any hits are the parent proteins of the fragment, which is the leakage the fold
clustering exists to exclude. The project's own record reaches the same place from the other
direction -- `structure-and-sequence-are-decoupled` (the best-matching window has 12% identity),
and S13's phi result (sequence context predicts phi at 36.1 deg against a sequence-blind 36.4).
There is no coevolution signal in a peptide, because coevolution is a statement about a family
and a peptide of this length has none.
REJECTED as unavailable, and this is the cleanest single statement of why this instrument is
harder than the protein case: the field's per-target conditioner is structurally absent at this
length, not merely weak.

### Template selection and per-target confidence heads

Template selection is the other per-target input, and it is leakage here by the same argument.
Self-reported confidence heads (pLDDT, PAE, pTM) are per-target but are FUNCTIONS OF THE MODEL'S
OWN OUTPUT, not external evidence -- and L_1 section 4.1 measured what that is worth at peptide
length: pLDDT has no within-target ranking skill on 588 peptides of 10 to 40 aa. A confidence head
trained here would be the same object.
REJECTED.

### The one class that is genuinely per-target and genuinely available

Quasi-single-model QA (L_1 section 2.2): condition on the AGREEMENT between the assessed answer
and an INDEPENDENT predictor's answer. CASP15's EMA winners are this class. The signal is
per-target, is not a function of the assessed model alone, and needs no MSA. The project's
blocker is measured and specific: its only second source is 31% angularly independent and 0.76 A
worse (S24 L2/L3), and mixtures sit on a line (cos 0.943). So the class is available in principle
and empty in practice for want of an independent generator. Recorded here because if S30 ever
obtains one, this is the conditioning family that applies, and its value is the independence
rather than the scoring rule.

---

## (b) The sign problem as its own literature

### The theorem: the per-target sign is an incidental parameter

Neyman J, Scott EL (1948), "Consistent estimates based on partially consistent observations",
Econometrica 16:1-32. The canonical example: X_ij ~ N(mu_j, sigma^2) independent, i = 1, 2 and
j = 1, ..., n. As n grows, the number of NUISANCE parameters mu_j grows with it while the number
of observations per nuisance parameter stays at 2. The maximum likelihood estimator of sigma^2
then converges to HALF its true value -- inconsistent, not merely inefficient. The general
statement: when the dimension of the incidental parameters grows with the sample size, ML for the
structural parameters need not be consistent.

Map it onto this project, which I believe is exact. The scorer decomposes as s(x) = f(x) + g(j)
with j the target: f is the STRUCTURAL parameter (shared across targets, what cross-target
training estimates) and g(j) -- the per-target sign of the in-band axis -- is an INCIDENTAL
parameter, one per target, with a bounded number of usable observations per target. S14's
measurement is the incidental-parameters problem observed empirically: within a target the model
reaches 0.986 with an overfitting gap of -0.0005 (g(j) is identifiable when you have that
target's own labels), across targets it delivers 0.600 (g(j) is not estimable from other targets'
answers), and the learning curve DECLINES with more targets because pooling averages over a
quantity whose sign differs per target. That is textbook.

### The sharper half, and it closes the loop with topic 7

The standard remedy for incidental parameters is to ELIMINATE them -- conditional likelihood,
fixed effects, the invariance principle. Topic 7 identified exactly this structure in the sprint's
own design: a within-group pairwise objective cancels any per-group additive term, s_i - s_j =
f(x_i) - f(x_j), which is conditional logistic regression on matched sets. But eliminating a
nuisance parameter is NOT estimating it. The conditional likelihood buys a consistent estimate of
the SHARED f precisely by discarding all information about g(j).

So the two things the project can do with cross-target data are: estimate f consistently while
discarding the sign (pairwise/conditional), or estimate f inconsistently while contaminating it
with the sign (pooled). Neither returns g(j). **The per-target sign is not estimable from other
targets' answers, as a matter of statistical theory and not of model capacity** -- which is the
formal version of S14's "capacity is saturated by a linear model" and of S12's flat learning
curve.

### When IS a sign estimable if the magnitude is not? The two escapes, named

The incidental-parameters literature gives exactly two routes, and they are the two the sprint
should carry into S30.

1. **Replication within the instance.** The inconsistency in Neyman-Scott comes from a fixed,
   small number of observations per nuisance parameter. If the number of INDEPENDENT observations
   of the same target's latent state grows, g(j) becomes estimable. Applied here: an independent
   second view of the SAME target -- a second generator, a second information source, an
   experimental observable. This is (a)'s quasi-single-model class and it is why independence,
   not scoring, is the valuable property.
2. **A covariate observed at inference that predicts g(j).** This converts an incidental parameter
   into a function of observables, and the structural model becomes estimable again. Applied here:
   exactly "a conditioning signal", and the project already has candidates with measured strength
   -- S14's native-free compactness proxies at residualised r = 0.244 (pool mean Rg), 0.313
   (distogram-predicted Rg) and 0.365 (incumbent emitted Rg) against the ORACLE's 0.909, all three
   CIs excluding zero.

### The adjacent literatures, and why they do NOT apply

The coordinator named phase retrieval and sign-ambiguous factor models. I checked both and neither
transfers, which is worth recording so S30 does not spend time on them.
- **Phase retrieval / the importance of phase** (the classical observation that a signal
  reconstructed from phase alone is far more recognisable than one reconstructed from magnitude
  alone) establishes that sign/phase information is disproportionately VALUABLE. It does not
  establish that it is RECOVERABLE without measurements that constrain it; phase retrieval
  algorithms succeed because they have magnitude measurements plus strong priors (support,
  positivity, oversampling). We have no analogous constraint tying the per-target sign to an
  observable.
- **Sign ambiguity in factor models / PCA** is a gauge freedom: the sign of a component is
  arbitrary and conventionally fixed, so the literature is about CONVENTION, not estimation. Our
  sign is not a gauge -- it is a real latent state with a right answer. Lane O's PC1 result
  (|rho| 0.37 to 0.39 per target, sign correct on 52%, leave-fold-out arm +0.0071 A worse,
  S29-L21) is the measured instance: the mode is real and its orientation is at chance.
- **One-bit compressed sensing** recovers a DIRECTION from many sign measurements, which is the
  converse of our problem (we want one sign, and we have no sign measurements).
RECORDED as three negative transfers.

---

## (c) The pool versus the posterior: dispersion as information

### The pointer being tested

Lane T's post-mortem: the pool's deviation from typical tracks the native better than the
posterior's median map does, consistent with the pool being made of REAL structures and the median
map not being one. This is topic 3's manifold argument in a new place -- the per-pair median map
is an unconstrained object that is generally not a realisable distance matrix, while every pool
member is a real peptide window. If that is right, then the pool's own geometry carries
constraints the posterior's summary does not, and the question is whether its DISPERSION (not its
centre) is readable.

### The mature literature on exactly this: spread-skill in ensemble forecasting

Numerical weather prediction has asked "is an ensemble's spread informative about this
forecast's error" for thirty years, which is our question with the names changed.

Whitaker JS, Loughe AF (1998), "The Relationship between Ensemble Spread and Ensemble Mean Skill",
Monthly Weather Review 126:3292-3302. The findings that transfer:
- "even for a perfect ensemble (one in which all sources of forecast error are sampled correctly)
  there need not be a high correlation between spread and skill". The limitation is STRUCTURAL,
  not a defect of the ensemble.
- the correlation "should be larger where the day-to-day variability of spread is large" -- i.e.
  the signal lives in how much the spread ITSELF varies from case to case, not in its level;
- spread "is likely to be most useful as a predictor of skill when it is EXTREME, that is, when
  it is either very large or very small compared to its climatological mean value".

Check against this instrument. The three statements are directly actionable and they are
falsifiable here. (i) A weak pool-spread-to-error correlation is the EXPECTED result even if the
pool is perfectly behaved, so a modest correlation is not evidence of a broken pool and a null is
not evidence of no information. (ii) The quantity to compute is not the pool's spread but the
CROSS-TARGET VARIABILITY of the pool's spread, length-residualised -- if that variability is
small, the literature predicts no usable signal and the route can be closed cheaply. (iii) The
signal should concentrate in the tails, which maps onto the project's standing regime question
(FAIL18 versus the 108) and gives a pre-registered place to look rather than a fishing expedition.
The project's own memory note `prediction-pool-disagreement-is-a-native-free-signal` already says
the contrast is demonstrated but its Angstrom value was never measured; Whitaker & Loughe say what
shape that measurement should have.
KEPT -- the best-matched external literature for the pool-dispersion route, with three concrete
design statements.

### What dispersion cannot do here, stated so it is not over-read

Spread is a MAGNITUDE statistic: it is symmetric under reflection of the error. By S29-L23's
|rho|(2q-1) form, what the system needs is the SIGN. A spread statistic can plausibly say "this
target is hard" (which is a magnitude claim and could gate an abstention or a per-target step
size); it cannot by itself say "the displacement should go this way rather than that way". So the
dispersion route, even if it works, addresses q's magnitude complement rather than q itself,
unless the dispersion is read ANISOTROPICALLY -- i.e. the shape of the spread (its principal axes
and their orientation relative to the candidate), not its scalar size. That anisotropic reading is
the version worth testing, and it is not what the memory note measured.
RECORDED as the distinction that decides whether the route can reach the sign at all.

---

## (d) My judgement: the single most promising direction for S30, and the first cheap measurement

**LABELLED AS JUDGEMENT, not measurement.** Putting (a), (b) and (c) beside what S29 has closed,
the honest picture is this. Lane T's bound (S29-L23) says every native-free operator over the
present information is worth one cosine and that no field built here exceeds |rho| 0.04 against
the 0.628 that 2.5 A requires; the strongest structured field known, the pool's PC1, has |rho|
0.37 to 0.39 with its per-target sign at chance, and would reach 2.98 A with a perfect sign. Topic
8(b) now adds that the missing sign is an INCIDENTAL PARAMETER, so no amount of cross-target
training, capacity, loss design or architecture can estimate it -- that is a theorem, not a
budget problem, and it converts S14's empirical closure into a structural one. Topic 8(a) says the
field's own per-target conditioner (MSA depth) is structurally absent at this length. So the
search for a better native-free SCORER should stop, and I would say so plainly in the report: it
is closed from three directions at once (measurement, derivation, and now statistical theory).
That leaves exactly two doors, which are the two escapes from the incidental-parameters problem,
and I would put S30 through the second one. **The most promising direction is to supply the
per-target sign of the pool's own principal mode from a native-free covariate -- the one operator
class whose value is already priced (2.98 A at a perfect sign, and |rho|(2q-1) at an imperfect
one), whose target is already identified (PC1, |rho| 0.37 to 0.39), and whose candidate
covariates already exist with measured strength (S14's compactness proxies at residualised r
0.244 to 0.365 against an ORACLE 0.909).** It is the only direction I have found where the
quantity to be bought, its price, and a concrete instrument to buy it with are all on the board
simultaneously -- though see the price table below, which caps what it can be worth. The first
door (a genuinely independent second generator) has the higher ceiling and is the one the field's
quasi-single-model winners use, but S24 closed generation on five instruments and the hardware
forbids the obvious sources, so its expected value per sprint-week is lower.

**The first cheap measurement, and it is cheap because it is ORACLE-free on one side.** Take the
126 targets' pool PC1 and its per-target ORACLE sign (already computed for S29-L21). Regress that
sign on the three native-free compactness proxies S14 already validated (pool mean Rg,
distogram-predicted Rg, incumbent emitted Rg), length-residualised, leave-fold-out, and report
the held-out SIGN ACCURACY q with a fold-clustered CI against the 0.50 null. That is one
regression over 126 rows of quantities that already exist -- minutes, no AMBER, no endpoint run.

**And price it before running it, because the price is sobering.** Lane T's formula,
gain = RMSD(1 - sqrt(1 - rho^2 (2q-1)^2)) at rho = 0.37 and RMSD_prod = 3.2126, gives (arithmetic
on T's formula, not a measurement; verified to reproduce T's own 2.98 A at q = 1):

| q (held-out sign accuracy) | 0.55 | 0.60 | 0.65 | 0.70 | 0.80 | 0.90 | 1.00 |
|---|---|---|---|---|---|---|---|
| gain (A) | 0.002 | 0.009 | 0.020 | 0.035 | 0.080 | 0.144 | 0.228 |
| built chain (A) | 3.210 | 3.204 | 3.193 | 3.177 | 3.132 | 3.069 | 2.985 |

**I had these four times too large in my first draft and the corrected table changes my own
conclusion, so I state the correction rather than quietly fixing it.** The gain is QUADRATIC in
(2q-1) -- the same squared-skill law the project already knows from
`decorrelated-errors-exist-but-are-unusable` -- so a sign predictor at q = 0.70, which would be a
genuinely good classifier of a latent state, buys 0.035 A, inside the noise of this instrument.
Even q = 0.80 buys 0.080 A. And the CEILING of this entire direction is 2.985 A at a perfect
sign, which is lane T's 2.98 A.

**So the revised judgement, which is less attractive than the paragraph above implies and is the
honest one.** The sign regression is still the right FIRST measurement -- it is minutes, it uses
only quantities that exist, and it is decisive in both directions. But it should be run as a
MECHANISM measurement, not as a route to the charter's target: on the arithmetic above, no
achievable sign accuracy on PC1 reaches 3.0 A, let alone 2.5 A. If the charter's 2.5 A is the
goal, then by T's own inversion it requires rho = 0.628 against a best-known structured field of
0.37, and the only lever in the record with that slope is the distance prior itself (-2.15 A per
unit, S24 L13; T's section 7 names it as one of three classes that can break assumption B2). My
judgement, stated plainly for the report: **S30's cheapest decisive experiment is the per-target
sign regression, and S30's only plausible route to a materially better number is a better distance
prior.** Those are two different projects and the report should not let the first stand in for
the second. If q comes back at chance, the native-free scorer programme is closed by theorem
(section (b)) plus measurement, and S30 should be a prior project or an honest ceiling paper.

---

## Verdict table (topic 8)

| # | source | what it adds | verdict |
|---|---|---|---|
| a1 | Jumper et al., Nature 596:583 (2021) -- recycling | recycling conditions on the model's OWN output; no external evidence enters | REJECTED as a conditioning source; the classification (optimisation, not information) is the deliverable |
| a2 | AlphaFold2 MSA-depth dependence (accuracy degrades below depth ~30; pLDDT tracks aligned residues per position; single-sequence mode as a distinct regime) | the field's real per-target conditioner is an ALIGNMENT statistic | REJECTED as unavailable -- structurally absent at 9-16 residues; the cleanest statement of why this instrument is harder than the protein case |
| a3 | template selection; self-reported confidence heads (pLDDT/PAE/pTM) | per-target but either leakage or a function of the model's own output | REJECTED (and L_1 4.1 measured pLDDT's in-band skill at peptide length: none) |
| a4 | quasi-single-model QA (CASP15 EMA winners) | agreement with an INDEPENDENT predictor: per-target, no MSA needed | AVAILABLE IN PRINCIPLE, empty in practice (S24 L2/L3); the family that applies if S30 obtains an independent generator |
| b1 | Neyman J, Scott EL, Econometrica 16:1-32 (1948) | the per-target sign is an INCIDENTAL PARAMETER; pooled ML is inconsistent for it; the remedy eliminates rather than estimates it | KEPT -- converts S14's empirical closure into a structural one, and names the only two escapes |
| b2 | phase retrieval / "the importance of phase" | sign/phase information is disproportionately valuable | RECORDED as a negative transfer -- it establishes value, not recoverability; phase retrieval needs magnitude measurements plus strong priors we do not have |
| b3 | sign ambiguity in factor models / PCA | the component sign is a GAUGE fixed by convention | RECORDED as a negative transfer -- our sign is a real latent state, not a gauge (lane O's PC1: 52%) |
| b4 | one-bit compressed sensing | recovers a direction from many sign measurements | RECORDED as a negative transfer -- the converse of our problem |
| c1 | Whitaker JS, Loughe AF, Mon Wea Rev 126:3292-3302 (1998) | spread-skill is weak EVEN FOR A PERFECT ENSEMBLE; larger where day-to-day variability of spread is large; most useful when spread is EXTREME | KEPT -- three concrete design statements for the pool-dispersion route, and a null-interpretation rule |
| c2 | my own distinction, labelled | spread is a MAGNITUDE statistic and is reflection-symmetric, so it cannot reach the sign unless read ANISOTROPICALLY | RECORDED -- decides whether the dispersion route can address q at all |

KEPT 2, AVAILABLE-IN-PRINCIPLE 1, RECORDED 5, REJECTED 3. No importable operator; the deliverables
are a theorem, a classification, three design statements and a judgement.

## What topic 8 tells the sprint, in four lines

1. **The per-target sign is not estimable from other targets' answers.** Neyman-Scott (1948): it
   is an incidental parameter, pooled ML is inconsistent for it, and the correct handling of the
   confound (conditional/fixed-effects, i.e. within-group pairwise) eliminates it rather than
   estimating it. S14's 0.986-versus-0.600 is that theorem measured.
2. **The only two escapes are replication within the target and a covariate observed at
   inference.** Everything the sprint has called "a conditioning signal" is the second escape.
3. **The field's per-target conditioner is MSA depth and it is structurally absent at 9-16
   residues**, which is the cleanest statement of why the peptide case is harder than the protein
   case. Recycling is optimisation, not information, and cannot substitute.
4. **For the pool-dispersion route:** expect a weak spread-skill correlation even in the best case,
   measure the cross-target VARIABILITY of the spread before anything else, look in the tails, and
   note that a scalar spread cannot reach the sign -- only an anisotropic reading can.
