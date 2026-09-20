# L_3 -- DECISION THEORY OF STRUCTURE POINT ESTIMATES, AND DISTOGRAM CALIBRATION (brief topic 3)

Lane L, 2026-09-20. Serves charter findings 4, 5 and 8, lane T's contraction work (S29-L7), and
contract rule 20 (the gradient cosine is gameable by shrinking toward typicality).

The answer, stated first. There is a theorem in the image-restoration literature that explains
charter finding 8 as a NECESSITY rather than a defect. Blau & Michaeli (CVPR 2018) prove that
distortion and perceptual quality trade off for ANY distortion measure: the estimator that
minimises expected distortion necessarily produces outputs whose DISTRIBUTION differs from the
distribution of real signals, and the deviation grows as distortion falls. In our terms: the
RMSD-optimal structure is necessarily a structure that does not look like a real peptide, so
every native-free scorer that measures realism -- which is what all QA methods, all statistical
potentials and all force fields measure -- must systematically disprefer it. Recognition and
RMSD-optimality are in tension by construction, and the project has measured both ends of that
curve without naming it.

---

## 1. The classical estimators, and which one the project ships

Classical decision theory (Berger, *Statistical Decision Theory and Bayesian Analysis*, 2nd ed.,
ch. 4; uncontroversial, stated here as textbook material rather than a fetched primary source):
under squared loss the Bayes estimator is the posterior MEAN; under absolute loss it is the
posterior MEDIAN; under 0-1 loss it is the MODE. Applied coordinatewise to a vector, each acts
independently per coordinate.

The project ships the L1 version in DISTANCE space: DIS is the L1 Bayes risk against the 17-bin
per-pair posterior, so the per-pair target is the posterior MEDIAN distance
(`docs/STATE_BRIEF` 5.1; `s27/ham_lib.py :: h_dis_mean` is the L2 sibling, the posterior mean).
The readout is then a coordinate AVERAGE over the top-m pool, i.e. an L2 estimator in coordinate
space. So the pipeline mixes an L1 estimator in distance space with an L2 estimator in coordinate
space. That mismatch is not obviously wrong, but it is a choice nobody justified, and it belongs
on lane M's convenience list.

### 1.1 Why the coordinate average contracts, in one line

The project measured that "averaging CONTRACTS the backbone 25.8%"
(memory: `averaging-space-beats-the-objective`). The mechanism is Jensen's inequality on the
norm, and it needs no calibration story: for any two random points X, Y,

    || E[X] - E[Y] ||  =  || E[X - Y] ||  <=  E[ ||X - Y|| ],

with equality only if X - Y is almost surely a fixed direction. So the distance between two
averaged CA positions is at most the average of the distances, strictly less whenever the members
disagree in direction. Averaging structures therefore shrinks every interatomic distance; the
more diverse the pool, the more it shrinks. This is exact, it is not a bias in the distogram, and
it cannot be fixed by reweighting the average -- only by projecting afterwards (which the built
chain does) or by not averaging in coordinate space.

Consequence for rule 20: a shrink toward typicality and the contraction of an average are the
same kind of object, and both improve a distance-space objective while degrading realism. The
meter's new implied-shrink field is measuring the right thing.

---

## 2. The result that explains finding 8

### 2.1 Blau Y, Michaeli T. "The Perception-Distortion Tradeoff." CVPR 2018 (arXiv:1711.06077)

Construction. Distortion is Delta(X, Xhat) = E[delta(X, Xhat)] for any per-sample measure delta;
perceptual quality is a DIVERGENCE d(p_X, p_Xhat) between the distribution of the estimator's
outputs and the distribution of real signals (with the total-variation choice, d is exactly the
best achievable probability of telling an output from a real signal). The
perception-distortion function is

    P(D) = min_{p_{Xhat|Y}}  d(p_X, p_Xhat)   subject to   E[Delta(X, Xhat)] <= D.       (11)

THEOREM 3 (the tradeoff), verbatim: "If d(p,q) of (4) is convex in its second argument, then the
perception-distortion function P(D) of (11) is 1) monotonically non-increasing; 2) convex." And
the paper states: "Theorem 3 requires no assumptions on the distortion measure Delta. This
implies that a tradeoff between perceptual quality and distortion exists for any distortion
measure". Convexity means "the tradeoff is more severe at the low-distortion and at the
high-perceptual-quality extremes".

THEOREM 1, the sharper form: "If p_{X,Y} defines a non-invertible degradation and the estimator
minimizing the mean distortion is unique, then Delta is not a stably distribution preserving
distortion" -- i.e. you cannot have an estimator that both minimises distortion and reproduces
the distribution of real signals.

THE MECHANISM, verbatim, and it is our mechanism: "The MMSE estimate is an average over all
possible explanations to the measured data, weighted by their likelihoods. However the average
of valid images is not necessarily a valid image, so that the MMSE estimate frequently 'falls
off' the natural image manifold."

Assumptions. A non-invertible degradation (many signals consistent with the observation); a
divergence convex in its second argument; nothing about the distortion measure. All hold here:
the "observation" is the sequence, many peptide structures are consistent with a distogram, and
the divergence choice is ours.

Check against this instrument, and this is the payload.
- Our distortion measure is CA-RMSD, the endpoint (charter section 10). Our "perceptual quality"
  is precisely what every native-free scorer in the S27 library computes: a statistical potential
  or force field is, up to a monotone map, a log-density of real structures, so scoring high
  means "looks drawn from p_X".
- Theorem 3 therefore says: as an estimator approaches the RMSD optimum, the divergence between
  its output distribution and the distribution of real peptides must INCREASE. The RMSD-optimal
  answer is the one realism scorers are most certain is fake.
- That is S28-L48 measured, with the sign explained: 20 of 31 scorers prefer the projected
  production average to a 0.25 A ORACLE structure, and CAGEO "prefers any real trace to a
  contracted one" (S28-L36). The scorers are not broken. They are reading the correct axis; the
  axis is orthogonal-to-adversarial to the one the endpoint rewards.
- The project's two operator families sit at the two ends of the curve, exactly as the theorem
  describes: the coordinate AVERAGE is the low-distortion, off-manifold end (contracted 25.8%,
  point cloud 3.0483 A), and the MEDOID / a pool member is the on-manifold, higher-distortion end
  (a real window, which is why realism scorers prefer it, and why consensus is capped at the
  pool's mode -- memory `consensus-is-outlier-avoidance`).
- The built chain's projection is a MEASURED point on this curve: `s12.instrument.project` maps
  the off-manifold average back onto ideal backbone geometry and costs +0.164 A (3.0483 cloud ->
  3.2126 chain, `s27/results/chain_rows.jsonl`). That is the price of realism in this system,
  already paid, already known. The coordinator's wave-2 probe P1 ("the projection price") is
  therefore a measurement of P(D)'s local slope, and should be framed that way.

Information test. What does it contain that the project does not? A THEOREM saying the
recognition failure cannot be engineered away by finding a better scorer, plus the statement that
the tradeoff is steepest exactly where we operate (low distortion). What it does NOT contain is
any way to get both. KEPT -- as the framing result for finding 8, and as a constraint on every
future objective proposal: an objective that is more "realistic" is, by this theorem, not more
RMSD-optimal, and the meter will see that as an anti-correlation.

### 2.2 The converse, also from the same paper, which protects against over-reading

"we could achieve perfect perceptual quality by randomly drawing natural images that have nothing
to do with the original ground-truth images. In this case the distortion would be quite large."
Realism alone buys nothing. This is the project's pool-member control in one sentence (S28-L36/
L37: CAGEO's apparent recognition was a preference for any real trace), and it is why the
pool-member control must stay on every future recognition claim.

### 2.3 What this does and does not license

It does NOT license "so recognition is impossible". The theorem constrains the DISTRIBUTION of an
estimator's outputs, not the per-target ranking of a fixed candidate set. Ranking two structures
of equal realism by accuracy is not forbidden by it. What it forbids is the hope that a realism
score, by itself, identifies the RMSD-optimal answer. Any future native-free ranker must
therefore either (a) be conditioned on the target in a way that is not a realism measure, or
(b) rank WITHIN a realism band (structures of matched plausibility), which is the in-band
discipline the project already enforces for a different reason.

---

## 3. Calibration of the posterior

### 3.1 Guo C, Pleiss G, Sun Y, Weinberger KQ. "On Calibration of Modern Neural Networks." ICML 2017 (arXiv:1706.04599)

Construction. Temperature scaling divides the logits by a single scalar T fitted on a validation
set by NLL. The paper's operative property: it "does not alter which class receives the highest
predicted probability (argmax) or change classification accuracy -- it only adjusts confidence
levels".

Check, and this is a clean explanation of a closed project result. S25 L2 found the shipped
posterior is about 2x over-confident and that CALIBRATING IT MAKES RMSD WORSE. Temperature
scaling is argmax-invariant; it is also MEDIAN-invariant for a symmetric per-pair posterior,
because rescaling the log-probabilities of a symmetric distribution about its centre leaves the
median where it is. So for an L1 Bayes readout, a pure temperature change can only act through
(i) the asymmetry of the 17-bin posterior and (ii) the discretisation (S25 L7: the per-pair
target takes only 17 values, with gaps up to 4 A). Both are artefacts rather than information,
which is why the measured effect had no reason to point the right way. The project's result and
the literature's property agree; the route is closed for a stated reason, not just empirically.
REJECTED as a lever (no new information; argmax/median-invariant), KEPT as the explanation of
why S25 L2 came out the way it did.

### 3.2 Gneiting T, Raftery AE. "Strictly Proper Scoring Rules, Prediction, and Estimation." JASA 102(477):359-378 (2007)

Construction. A scoring rule S(P, y) is PROPER if E_{y~Q}[S(Q,y)] <= E_{y~Q}[S(P,y)] for all P,
strictly proper if equality only at P = Q: honest reporting of the true predictive distribution
is optimal. The logarithmic score and the CRPS are the canonical examples; the CRPS is strictly
proper for distributions with finite first moment and is the generalisation of absolute error to
predictive distributions. The paradigm the paper advocates is "maximising sharpness subject to
calibration".

Check. The distogram is trained by cross-entropy over 17 bins, which IS the logarithmic score and
IS strictly proper -- so the training objective is already correct in the decision-theoretic
sense, and no change of scoring rule will make the posterior better calibrated in a way that
survives (S25 measured this over 51 arms: corr(progress-toward-truth, endpoint) = +0.054).
Information test: none for training. What is worth recording is the DIAGNOSTIC half: CRPS is the
right one-number summary of a per-pair predictive distribution against a realised distance, and
it is distance-sensitive where the log score is not ("the logarithmic score assigns harsh
penalties regardless of closeness"). If any lane wants to compare two posteriors, CRPS is the
defensible metric and MAE is not -- which the project independently found
(`prior-mae-prices-selected-rmsd`: r = 0.19; `error-shape-not-mae-decides-ranking`).
KEPT as the correct diagnostic; REJECTED as a training change.

### 3.3 What AlphaFold / trRosetta actually did about distance uncertainty

Recorded for completeness because the brief asks. AlphaFold1 (Senior et al., Nature 577:706,
2020) fits a smooth potential to the predicted distance histogram and MINIMISES it by gradient
descent in torsion space, i.e. it uses the posterior as an energy rather than taking a per-pair
point estimate; AlphaFold2 (Jumper et al., Nature 596:583, 2021) drops the distogram-as-potential
route entirely in favour of a structure module that emits coordinates directly, and predicts its
own per-residue accuracy (pLDDT, fitted lDDT-Ca = 0.997 pLDDT - 1.17, Pearson r = 0.76). Neither
applies a post-hoc temperature to the distogram. The lesson for this project is the ROUTE, not a
calibration trick: the field's answer to an over-confident per-pair posterior was to stop taking
a per-pair point estimate, not to rescale it. The project's closest analogue to "use the whole
posterior jointly rather than per-pair" is lane X's configuration-space formulation.

---

## 4. Verdict table (topic 3)

| # | paper | what it adds | verdict |
|---|---|---|---|
| 1 | Classical Bayes estimators (L2 -> mean, L1 -> median, 0-1 -> mode); Jensen contraction ||E X - E Y|| <= E||X - Y|| | the one-line proof that coordinate averaging must contract every distance | KEPT as framing: the 25.8% contraction is Jensen, not a distogram bias |
| 2 | Blau & Michaeli, Perception-Distortion Tradeoff, CVPR 2018, Thm 1 and Thm 3 | a THEOREM that the distortion-optimal estimator's output distribution must differ from the real one, for ANY distortion measure, most severely at low distortion | KEPT -- the framing result for charter finding 8; recognition vs RMSD-optimality is a necessity |
| 3 | the same paper's converse ("perfect perceptual quality by randomly drawing natural images") | realism alone buys no accuracy | KEPT -- the pool-member control's theoretical form |
| 4 | Guo et al., temperature scaling, ICML 2017 | argmax-invariance (and, for symmetric posteriors, median-invariance) | REJECTED as a lever; KEPT as the explanation of S25 L2's "calibration makes it worse" |
| 5 | Gneiting & Raftery, JASA 102:359 (2007) | properness of the log score (the distogram's training loss is already strictly proper); CRPS as the distance-sensitive diagnostic | REJECTED as a training change; KEPT as the correct posterior-comparison metric, in place of MAE |
| 6 | Senior et al. Nature 577:706 (2020); Jumper et al. Nature 596:583 (2021) | the field's response to per-pair over-confidence was to stop taking per-pair point estimates | KEPT as a route note supporting lane X, not as an importable method |

KEPT 5 (as framing, diagnostics or route notes), REJECTED 2 as levers (post-hoc calibration;
changing the training scoring rule). No importable operator.

## 5. What topic 3 tells the sprint

1. Finding 8 now has a theorem behind it. The RMSD-optimal estimator is necessarily the one a
   realism scorer rates worst, and the tension is sharpest exactly at the low-distortion end
   where the project operates. Stop looking for a native-free scorer that prefers the ORACLE
   structures; by Theorem 3 a realism-type scorer cannot, and every scorer in the library is a
   realism-type scorer.
2. The corollary is constructive, and it is the one new direction topic 3 yields: rank WITHIN a
   realism band rather than across it. Two structures of matched plausibility are not separated
   by the theorem; the production average and a 0.25 A ORACLE structure are, and that is the
   comparison S28-L48 ran. A recognition test that first matches realism (the pool-member control
   is a crude form of this) is the only version of finding 8's question that is not answered
   "no" in advance.
3. The 25.8% contraction is Jensen's inequality, not a bias; the +0.164 A projection price is a
   measured point on P(D). Probe P1 should be framed as measuring the local slope of the
   perception-distortion curve, which makes its result interpretable whichever way it comes out.
4. Post-hoc calibration is closed for a mechanism reason (argmax/median invariance), which
   upgrades S25 L2 from an empirical dead end to an explained one.
5. Use CRPS, never MAE, to compare two posteriors. The project reached this empirically twice.
