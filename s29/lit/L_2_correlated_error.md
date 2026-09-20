# L_2 -- CORRELATED ERROR IN ENSEMBLES, AND THE METHODS THAT BREAK IT (charter finding 11; brief topic 2)

Lane L, 2026-09-20. Serves the coordinator's H1 (the typicality axis: use the common-mode error
as signal rather than averaging it in) and charter finding 11 (the pool's error is 68%
common-mode, so any method whose power comes from agreement is structurally limited).

Checked against: 126 peptides of 9-16 residues; a 500-member pool of real CA windows retrieved
by BLOSUM sum, top-75 averaged in the medoid frame; S23 L9's exact error decomposition; CPU
only. Notation follows the project: members w_k = t + e_k with t the native, c the average,
d_k = w_k - c, so e_k = (c - t) + d_k with mean_k d_k = 0.

The answer, stated first. The field's central ensemble law is ALREADY IN THIS PROJECT, derived
independently and verified to 2.7e-14 A: S23 L9's identity mean_k|e_k|^2 = |ebar|^2 +
mean_k|d_k|^2 IS the Krogh-Vedelsby ambiguity decomposition in coordinate space. Everything the
literature builds on top of that law to BREAK correlated error requires one of exactly three
inputs, and this instrument has none of them: (i) control over how the members are TRAINED
(negative correlation learning, boosting, diverse ensembles) -- our members are retrieved real
windows, not fitted estimators; (ii) a control variate whose mean is KNOWN, or samples of the
high-fidelity quantity (control variates, multifidelity Monte Carlo) -- the known mean here is
the native; (iii) a known bias RATIO or error exponent (Richardson extrapolation, the two-source
contrast that H1 proposes) -- S24 L2/L3 measured the two sources' geometry and it does not
supply one. A fourth family, self-correction by recycling, is a training-time property of a
structure predictor we cannot retrain.

One quantitative consequence for the sprint, derived below in section 1.3: under the same law,
enlarging the averaged set from 75 members to an INFINITE pool of the same kind is worth about
0.008 A. The averaging route is closed by arithmetic, not by opinion.

---

## 1. The law itself

### 1.1 Krogh & Vedelsby (1995), via Brown, Wyatt and Tino, JMLR 6:1621-1650 (2005), eq (10)

Verbatim: "at a single arbitrary datapoint, the quadratic error of the ensemble estimator is
guaranteed to be less than or equal to the weighted average quadratic error of the component
estimators,

    (f_ens - t)^2 = sum_i c_i (f_i - t)^2 - sum_i c_i (f_i - f_ens)^2                    (10)

where sum_i c_i = 1, c_i >= 0, and f_ens = sum_i c_i f_i." The second term is the AMBIGUITY.

Check against this instrument -- and this is the finding. With uniform c_i = 1/m and the
project's notation, eq (10) reads |c - t|^2 = mean_k|e_k|^2 - mean_k|d_k|^2, which is S23 L9's
identity rearranged:

    mean_k |e_k|^2   =   |ebar|^2   +   mean_k |d_k|^2
                         160.36         63.82
                         COMMON         IDIOSYNCRATIC = the AMBIGUITY

(`s23/errdecomp.py`, `s23/results/errdecomp.json`, exact to 2.7e-14 A on 126/126). The project
called the two terms common-mode and idiosyncratic; the literature calls them ensemble error and
ambiguity. They are the same two terms. The common-mode fraction f = 0.676 is not a pathology of
this pool -- it is where this pool sits on a law that holds for every ensemble.

Information test: the identity is already ours. What the literature adds is the EXPECTATION form
(section 1.2), which turns the identity into a prediction about m, and the vocabulary that
connects our result to a large body of negative results about fixing it.
KEPT as framing, not as method.

### 1.2 The bias-variance-covariance decomposition -- Ueda & Nakano (1996), as given in Brown et al. (2005) eqs (6) to (9)

    bias  = (1/M) sum_i (E{f_i} - t)
    var   = (1/M) sum_i E{(f_i - E{f_i})^2}
    covar = (1/(M(M-1))) sum_i sum_{j != i} E{(f_i - E{f_i})(f_j - E{f_j})}

    E{(fbar - t)^2} = bias^2 + (1/M) var + (1 - 1/M) covar                               (9)

This is the law that matters for finding 11, because of the coefficient on covar: it is
(1 - 1/M), which does NOT decay. Write rho = covar/var for the average member-member
correlation; the reducible factor multiplying var is 1/M + (1 - 1/M) rho, whose limit is rho.
At the project's rho = 0.676:

| M | 1 | 5 | 25 | 75 | 150 | 500 | infinity |
|---|---|---|---|---|---|---|---|
| 1/M + (1-1/M) rho | 1.000 | 0.741 | 0.689 | 0.680 | 0.678 | 0.677 | 0.676 |

Between M = 75 (production's top-75) and M = 500 (the whole pool) the factor moves by 0.0037.
This is the literature's explanation for three separate project results that were previously
only empirical: S17 (widening K makes the realised answer worse -- the averaging gain from extra
members is ~0.4% of a member's variance, so any displacement of good members from the shortlist
dominates it); S24's pool union worth +0.0022 A; and the memory note "the terminal operator
consumes the set MEAN, not the set BEST".

### 1.3 The derived number the coordinator should have: what more members can ever buy

From S23's own two numbers, under exchangeability of the members (the standard assumption that
produces the i.i.d. 1/m term): E|c_m - t|^2 = |c_inf - t|^2 + (1/m) E|w - c_inf|^2 with
E|w - c_inf|^2 = mean_k|d_k|^2 * m/(m-1) = 63.82 * 75/74 = 64.68, so |c_inf - t|^2 =
160.36 - 64.68/75 = 159.50. Converting on the recorded anchor (point cloud 3.0483 A at m = 75,
`s27/results/chain_rows.jsonl`):

    infinite-pool point cloud = 3.0483 * sqrt(159.50/160.36) = 3.040 A,  a gain of 0.008 A.

CAVEATS, stated because this is arithmetic on aggregate numbers, not a measurement: it uses the
aggregate 160.36/63.82 rather than per-target decompositions; the same conversion applied to the
member scale gives 3.604 A against the recorded typical-member 3.7037 A, so the aggregate-to-RMSD
map is approximate at the ~3% level; and any residual correlation among the d_k makes the gain
SMALLER, not larger, so 0.008 A reads as an upper bound on the averaging route. If the
coordinator wants it exact it is a five-minute per-target recomputation in lane O off
`s23/results/errdecomp.json`; I run no experiments.

Consequence: "average more members" is closed at the 0.008 A level. Any remaining leverage in
aggregation must come from WHICH members or from breaking the common mode, not from how many.

---

## 2. Methods that break correlated error, and what each one needs

### 2.1 Negative correlation learning -- Liu & Yao (1999), analysed in Brown, Wyatt & Tino, JMLR 6:1621 (2005)

Construction. Train the members jointly with a penalty that rewards anti-correlated errors. The
diversity-encouraging error the paper derives is

    e_i^div = (1/M) sum_i (1/2)(f_i - t)^2 - kappa (1/M) sum_i (1/2)(f_i - fbar)^2        (17)

with gradient (1/M)[(f_i - t) - kappa (f_i - fbar)] (eq 18): kappa = 0 trains the members
separately, kappa = 1 trains the ensemble as a single unit. NC's penalty coefficient lambda has
a proven upper bound from the Hessian's positive definiteness,

    lambda_upper = M/(M-1),   gamma_upper = M^2 / (2 (M-1)^2)                            (39)

("When lambda or gamma is varied beyond these upper bounds, the Hessian matrix is guaranteed to
be non-positive definite"; asymptotically lambda_upper -> 1, gamma_upper -> 0.5).

Assumption. The members are ESTIMATORS BEING TRAINED, with a loss and a gradient. Our members
are real CA windows retrieved from other proteins. There is no parameter to push. The only
trainable object in the pipeline is the distogram, and there is exactly one of it (not an
ensemble), and its errors are the common mode.

Information test: a mechanism for producing anti-correlated members, which requires a generator
with knobs. The project's generator is retrieval; S24 closed learned/from-scratch generation on
five instruments. REJECTED -- the method's input (trainable members) does not exist here.
NOTE for lane X: if a lane ever builds a GENERATOR with parameters, eq (17) with its kappa is
the right way to ask it for decorrelated members, and eq (39) bounds the knob.

### 2.2 Diverse ensembles, and the modern reversal -- Abe, Buchanan, Pleiss & Cunningham, "Pathologies of Predictive Diversity in Deep Ensembles" (arXiv:2302.00704)

Construction. Decompose the ensemble risk as

    R_ens = R_avg - E[ (1/M) sum_i l(f_i(x), y) - l(fbar(x), y) ]                         (3)

where the bracket (the Jensen gap) is the predictive diversity, non-negative for strictly convex
losses, zero iff all members predict identically. The paper studies ~600 ensembles and finds,
verbatim from the abstract: "these intuitions do not apply to high-capacity neural network
ensembles ... interventions [that trade component performance for diversity] can improve the
performance of small neural network ensembles ... but they harm the performance of the large
neural network ensembles most often used in practice", and "discouraging predictive diversity is
often benign". Even free diversity carries "an opportunity cost": the best members make nearly
identical predictions, and ensembling the BEST members beats ensembling the most diverse ones.

Check. This is the external form of a result the project already has: the consensus/medoid family
(9 arms) closed 0% of the in-pool gap and "diversity-maximising selection: dead"
(`docs/STATE_BRIEF` 5.2). The literature's mechanism is the same trade-off: buying diversity
costs member quality, and when members are already the best available the exchange is negative.
Information test: none for this instrument. REJECTED as a method; KEPT as external confirmation
that "make the pool more diverse" is a route the field has also abandoned.

### 2.3 Control variates -- Glynn & Szechtman, "Some New Perspectives on the Method of Control Variates", in Monte Carlo and Quasi-Monte Carlo Methods 2000, Springer (2002)

Construction, verbatim. "suppose that there exists a random variable Y, jointly distributed with
X, for which EY is KNOWN. Then, the control variate C = Y - EY is guaranteed to be a 'mean zero'
random variable, so that X(lambda) = X - lambda C is an estimator for alpha." The optimal
coefficient is lambda* = cov(X, C)/var C; in the vector case lambda* = (E C C^T)^{-1} E X C
(their eq 1), which the paper identifies as a Hilbert-space projection: minimising the variance
"is equivalent to finding the closest point W in G to X - alpha", characterised by
<X - W*, W> = 0.

Assumption, and it is the whole thing. E Y must be KNOWN. The variance falls by the factor
(1 - R^2) where R is the correlation between the estimator and the control -- but only because
the control's mean is available to subtract.

Check. Here the estimator is the pool average and the candidate controls are other native-free
functionals of the same pool (Rg, typicality, a second retrieval source). We do not know the
expectation of any of them under the correct (native-centred) distribution; knowing it would BE
knowing the native. This is the exact, formal statement of "a bias shared by every member is
invisible to any within-pool statistic": the projection in section 2.3's Hilbert space is onto
the span of ZERO-MEAN controls, and we cannot centre any control without the native.
Information test: the control's known mean. Not available. REJECTED -- and this is the sharpest
way to say why H1's family is hard: an uncentred control does not reduce variance, it moves the
estimator by an unknown amount.

### 2.4 Multifidelity Monte Carlo -- Peherstorfer, Willcox & Gunzburger, SIAM Review 60:550-591 (2018)

Construction. When the low-fidelity model's mean is not known, the estimator is
sMF = sbar_hi(m0) + sum_i alpha_i (sbar_lo^(i)(m_i) - sbar_lo^(i)(m_{i-1})) with m_i > m0. The
paper's key properties: it "is an UNBIASED estimator of E[f_hi]"; the optimal coefficients are
alpha_i* = rho_i sigma_hi / sigma_i (their eq 3.12); and it beats plain Monte Carlo iff

    sqrt(1 - rho_1^2) + sum_i sqrt(c_lo^(i)/c_hi) (rho_i^2 - rho_{i+1}^2) < 1            (3.16)

"both correlation and costs of the models are critical".

Assumption, and it is fatal here. The estimator is anchored by m0 SAMPLES OF THE HIGH-FIDELITY
MODEL. Multifidelity does not remove a low-fidelity model's bias; it reduces the VARIANCE of an
estimator of the high-fidelity quantity, using the low-fidelity model only as a correlated
auxiliary. Our high-fidelity model is the native structure and we have zero samples of it at
deployment.

Information test: samples of the truth. REJECTED. Recorded because it closes a tempting reading
of H1: "combine a biased cheap source with a biased expensive source and extrapolate" is NOT what
multifidelity does -- it always keeps an unbiased anchor.

### 2.5 Richardson-style extrapolation / the two-source contrast (H1's own family)

Construction (classical; no single primary paper, the construction is textbook). If an estimator
has an error expansion A(h) = A + C h^p + O(h^{p+1}) with the exponent p KNOWN, then two
evaluations at h and h/2 cancel the leading term:
A_R = (2^p A(h/2) - A(h)) / (2^p - 1) = A + O(h^{p+1}).

Assumption. The bias of the two evaluations differs by a KNOWN factor (2^{-p}) along a COMMON
direction. That is the only thing that makes the cancellation exact; with an unknown factor the
extrapolation is a one-parameter guess, and with a bias component that is not shared it amplifies
whatever is orthogonal.

Check against H1, and this is the useful part. H1 proposes u = (conditioned answer) - (blind
answer) as the native-free direction and a leave-fold-out scalar step along it. In Richardson's
terms that is assuming (a) the two sources' biases are parallel and (b) their ratio is stable
across targets and learnable from folds. The record measured (a) directly: S24 L2/L3 gives a bias
cosine of 0.647 between the blind library source and the incumbent (31% angularly independent)
with q = 1.231 (the blind source is 0.76 A WORSE), and score-selected mixtures sit on a straight
line at cos 0.943. Parallel-plus-known-ratio is exactly what 0.647 with a worse-quality second
source does not give you: the non-parallel 31% is the part that gets amplified, and it belongs to
the WORSE source. Note also that this is the coordinator's own stated objection in `s29/STATE.md`
-- the literature does not rescue it; it names the condition that is missing (a known ratio along
a shared direction) and confirms the failure mode (amplification of the orthogonal part).
REJECTED as stated; the probe remains worth its 20 minutes because the extrapolation SIGN was
never run, but the literature's verdict is that a two-source contrast without a known bias ratio
is a one-parameter fit, and must be priced as one (leave-fold-out, with the matched control in
the operator's own space).

### 2.6 Boosting / residual fitting

Construction. Fit a second model to the residual of the first: F_{m+1} = F_m + nu h_m with h_m
trained on the negative gradient of the loss at F_m.

Assumption. The residual is LEARNABLE FROM FEATURES available at inference, and the second
learner's errors are not the first learner's errors. Both fail here, and both are measured.
S24 closed the learned residual / from-scratch generation route by oracle upper bound and corpus
census. And the memory note `error-coherence-decides-correctors` is the mechanism in its sharpest
form: at identical 0.688 sign accuracy, COHERENT mistakes emit +0.31 A while i.i.d. mistakes emit
-0.14 A, because a corrector trained on the predictor's own features inherits its error structure.
That is the boosting assumption failing for a stated, measured reason.
REJECTED -- and this is the strongest in-house result in topic 2: the project has measured the
precise condition (error coherence, not accuracy) under which residual fitting helps or hurts.

### 2.7 Self-correction by recycling -- Jumper et al., Nature 596:583-589 (2021)

Construction. AlphaFold2 feeds its own output (the pair representation, the single
representation, and the predicted CB coordinates) back into the trunk and re-runs it; the final
loss is applied to the outputs of each recycling iteration (4 iterations over 48 Evoformer
blocks). The network is TRAINED under recycling, so it learns to improve on its own previous
answer.

Assumption. You can train the predictor with its own output in the loop. Here the distogram is
trained leave-fold-out on 6,800 windows on a CPU box, and the "previous answer" whose systematic
error we would want it to correct is the pool average -- so recycling would require retraining
the distogram conditioned on pipeline output, i.e. a structure-conditioned predictor, which the
project excludes for leakage and cannot afford.
Information test: a training procedure, not an inference-time signal. REJECTED on both counts.

Recorded beside it, because it is the cleanest external version of an in-house lesson: AF2's
pLDDT regresses on true lDDT-Ca with a least-squares fit lDDT-Ca = 0.997 pLDDT - 1.17, Pearson
r = 0.76, ACROSS a broad quality range -- and the same quantity has NO within-target ranking skill
on peptides (L_1 section 4.1). A signal can be strongly calibrated globally and worthless
in-band. This project's own rule, "in-band is the only ranking metric", is here confirmed on the
field's best-known confidence score.

---

## 3. Verdict table (topic 2)

| # | paper | what it needs | verdict here |
|---|---|---|---|
| 1.1 | Krogh & Vedelsby (1995) ambiguity decomposition, via Brown et al. JMLR 6:1621 (2005) eq (10) | nothing -- it is an identity | KEPT as framing: it IS S23 L9's identity; the project re-derived the field's central ensemble law independently |
| 1.2 | Ueda & Nakano (1996) bias-variance-covariance, eq (9) | nothing -- an expectation identity | KEPT: the (1 - 1/M) coefficient on covar is why finding 11 caps aggregation; gives the M-table and the 0.008 A bound of section 1.3 |
| 2.1 | Negative correlation learning, Liu & Yao (1999); bounds in Brown et al. (2005) eq (39) | trainable members with a gradient | REJECTED -- our members are retrieved windows; noted for any future parameterised generator |
| 2.2 | Abe et al., Pathologies of Predictive Diversity (arXiv:2302.00704) eq (3) | nothing; it is a negative result | KEPT as confirmation -- diversity interventions harm good ensembles; matches "diversity-maximising selection: dead" (S12/S17) |
| 2.3 | Glynn & Szechtman, control variates (2002), eq (1) | a control whose mean is KNOWN | REJECTED -- the known mean is the native; formal statement of why within-pool statistics cannot see the common mode |
| 2.4 | Peherstorfer, Willcox & Gunzburger, SIAM Rev 60:550 (2018), eqs (3.12), (3.16) | samples of the HIGH-FIDELITY quantity | REJECTED -- it reduces variance around an unbiased anchor; we have no anchor |
| 2.5 | Richardson extrapolation / two-source contrast (H1's family) | a KNOWN bias ratio along a SHARED direction | REJECTED as stated -- S24 L2/L3 measured cos 0.647 with the second source 0.76 A worse; the non-parallel 31% is what gets amplified |
| 2.6 | Boosting / residual fitting | a residual learnable from features, with decorrelated errors | REJECTED -- S24 closed it; `error-coherence-decides-correctors` gives the exact condition it violates (+0.31 A coherent vs -0.14 A i.i.d. at the same accuracy) |
| 2.7 | Recycling, Jumper et al., Nature 596:583 (2021) | training the predictor with its own output in the loop | REJECTED -- training-time, leakage and hardware; its pLDDT r = 0.76 globally vs no in-band skill at peptide length is the in-band lesson again |

KEPT (4, all as framing or confirmation): the two decompositions; the predictive-diversity
pathology; and (from 2.7) the global-vs-in-band contrast. REJECTED (5 families): negative
correlation learning, control variates, multifidelity, Richardson-style two-source extrapolation
as stated, and residual/boosting correction.

## 4. What topic 2 tells the sprint

1. Finding 11 is a law, not a defect, and the project already has it exactly. Quoting S23 L9 as
   "the ambiguity decomposition, verified to 2.7e-14" makes it citable outside the project.
2. The averaging route is closed by arithmetic: more members is worth <= ~0.008 A (section 1.3).
   Aggregation leverage must come from WHICH members or from breaking the common mode.
3. Every literature method that genuinely breaks correlated error needs one of: trainable
   members, a control with a known mean, samples of the truth, or a known bias ratio. This
   instrument has none. That is a family-level negative result, and it is the honest form of
   "consensus and averaging inherit common-mode error rather than cancelling it".
4. For H1 specifically: the literature does not supply the missing ingredient, and it names it
   precisely -- a KNOWN ratio between the two sources' biases along a SHARED direction. The
   probe should therefore be priced as a one-parameter leave-fold-out fit with its control in the
   operator's own space, and its failure mode (amplifying the orthogonal component of the worse
   source) is predicted by both the literature and S24 L3.
5. The one place this family could ever act here is a parameterised GENERATOR (eq 17's kappa with
   eq 39's bound). Lane X's configuration space is the only S29 direction with that property.
