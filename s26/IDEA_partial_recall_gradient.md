# IDEA_partial_recall_gradient -- DOES THE FOLD MODEL'S MEMORISATION HAVE A GRADIENT BELOW THE 0.6 IDENTITY THRESHOLD? (lane W, Sprint 26; found on the way to L44)

## Hypothesis

L44 Part C measured the top of the memorisation curve: a fold model that had the target's OWN
native among its training labels emits a built chain 0.70 A nearer to it (fold CI [-0.83, -0.58],
5/5 folds) and an argmin 1.21 A nearer. The fold discipline removes only the identity-1.0 end of
that curve (every training chain is below 0.6 longer-normalised identity to the target). The
hypothesis is that the curve has a GRADIENT below the threshold: a target whose sequence sits at
0.4 to 0.58 identity to some training chain, or whose sequence is CONTAINED in a training chain
at high shorter-normalised identity (the near-copies lane I's L15 found to be at the null as a
CLUSTERING criterion, which says nothing about their effect on the MLP), receives a posterior
partially recalled from that chain's native, and its built-chain RMSD is lower (if the recalled
structure is near the target's native) or higher (if it is not; L52 says the copy is 2.9 A away
in the median). Native-free covariates, sequences only: (i) I_long = max longer-normalised
identity of the target to any chain in its fold model's training corpus (the pinned criterion,
< 0.6 by construction); (ii) I_short = max shorter-normalised identity (containment); (iii) the
length of the longest exact k-mer shared with any training chain. ORACLE covariate for the
mechanism check only: the cross-deposit RMSD between the target's native and the best-matching
training chain's segment (L52's quantity, generalised past verbatim copies).

## Why the record does not already close it

- S22 L3 (`score_mean`, `pool_spread`: r 0.36 to 0.43 with true difficulty, wrong quartile), S22
  L7 (four router constructions on distogram confidence, score-distribution shape and candidate
  geometry: all fail, one significantly harmful), S23 L3 (eight native-free features against
  s*, all |rho| <= 0.11), S23 L7 (rg_z, a signal from OUTSIDE the objective with skill on the
  ordering task, routes nothing): every native-free feature tried is a functional of the pool,
  the posterior or the emission. None is a property of the TRAINING CORPUS relative to the
  query. Sequence identity to the training set was never entered as a covariate of difficulty
  (`docs/FINDINGS.md:485` records "max identity to any training peptide 0.422 / 0.571" as a
  corpus statistic only).
- S10-4 lists the 13 targets with a >= 0.6 pool WINDOW and their RMSDs (no correlation drawn);
  L15 shows containment at 0.6 is at the null as a clustering threshold (a real dev sequence
  passes against 0.5% of members, a shuffle against 0.3%): that is about chance CONTAINMENT,
  not about whether a chance-level near-copy moves the MLP.
- The memorisation itself (L44 Part C) is measured at identity 1.0 only; L44 Part B (n = 1)
  measured a carrier's presence at identity 0.52 (longer) / 1.0 (shorter) as worth 0.011 A.
  Nothing between 0.52 and 1.0 in the longer measure, and nothing below 1.0 in the shorter
  measure, has been measured.
- This is not a router: nothing is converted into an operator. It is a leak audit of the fold
  discipline's threshold and, if a gradient exists, a difficulty covariate the report must state
  beside every dev number (the memory `benchmark-and-folds-must-be-pinned` says why the threshold
  cannot be moved after the fact: the answer here is a caveat or a clean bill, not a re-pin).

## Exact falsifier

Spearman rho between each covariate (I_long, I_short, shared k-mer length) and `rmsd_arm`
(production cache) over the 126 dev targets, with a fold-clustered bootstrap CI and a 500-draw
label-permutation null, one-sided in the registered direction (higher identity, lower RMSD);
Bonferroni over the three covariates. A gradient EXISTS iff at least one rho <= -0.25 with the
fold CI excluding zero and the permutation p < 0.017. The registered expectation is that it does
NOT exist (rho within +-0.15): the S19 finding "every predictor emits a typical peptide of that
length" and the 0.024 A margin over pool-typicality (S7-8) say the MLP interpolates rather than
recalls below identity 1.0, and L52 says even a verbatim copy is not near-native. If a gradient
does exist, the ORACLE mechanism check runs: the posterior's MAE and gam_eff/cos against the
native (p_ladder's `progress`) as a function of I_short, to say whether the recalled structure is
the target's native or the neighbour's. A positive re-runs on the 60 targets of the largest fold
pair as a split-half transfer (not a second seed: the covariates are deterministic).

## Expected effect against the computed MDE

For a Spearman correlation at n = 126, SE(rho) ~ 1/sqrt(n - 3) = 0.090, so MDE = 2.8016 x 0.090 =
0.25 in rho; the falsifier's threshold is set at that MDE. Expected rho -0.05 to +0.05: below the
MDE, so the deliverable is the power statement ("a gradient of |rho| >= 0.25 is excluded") and
the clean bill for the 0.6 threshold, or, if the covariate reaches it, a caveat with a number on
every dev figure.

## Memory and agent-hours

Sequences of 787 peptides and 6,003 fragments through `core.data.identity_many` with the
composition prefilter (lane I's audit did the same for 126 x 786 in 15 s); the k-mer scan is
trivial; `rmsd_arm` from the cache. Under 0.2 GB, under 5 minutes; one governed job. 1.5
agent-hours including the write-up. Gated (reads `rmsd_arm`), needs its PREREG first.

## Information value

Either a clean bill for the fold discipline's threshold with a power statement (the report's
"the fold models never saw a near-copy that helped them" sentence, currently unsupported), or a
measured gradient that puts a caveat on every dev number and a magnitude on EXAMINATION E4 that
L44's identity-1.0 envelope only bounds from above. Plausibility of a gradient 0.2; value of
either outcome high because the sentence is needed either way.
