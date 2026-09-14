# PREREG_partial_recall_gradient -- DOES THE FOLD MODEL'S MEMORISATION HAVE A GRADIENT BELOW THE 0.6 IDENTITY THRESHOLD? (lane W, Sprint 26; L77 extension)

Written 2026-09-13 22:20, before any covariate was correlated with any RMSD. Idea file:
`s26/IDEA_partial_recall_gradient.md`. Code: `s26/w_recall.py`. Results:
`s26/results/w_recall_*.json`. The covariates are native-free (sequences only) and are computed
first; the correlation with `rmsd_arm` is gated and reads the production cache's RMSDs (Phase 0
signed off, L33).

## 0. What the record knows

L44 Part C: a fold model that saw the target's OWN native (identity 1.0) emits a built chain
0.70 A nearer (fold CI [-0.83, -0.58]). L44 Part B (n = 1): a carrier at longer identity 0.52 /
containment 1.0 is worth 0.011 A. Nothing between. S22 L3 / L7 and S23 L3 / L7: every native-free
difficulty feature tried is a functional of the pool, the posterior or the emission; none is a
property of the training corpus relative to the query. L15: containment at 0.6 is at the null
as a CLUSTERING criterion (0.5% real vs 0.3% shuffled pass rate), which says nothing about the
MLP's response to a chance-level near-copy. `docs/FINDINGS.md:485`: "max identity to any
training peptide 0.422 / 0.571" as a corpus statistic only.

## 1. Hypotheses and exact falsifier

**H_R.** Over the 126 dev targets, higher sequence proximity of the target to its own fold
model's TRAINING chains predicts a lower built-chain RMSD (a partial recall of a training
native that is near the target's native), or a higher one (recall of a training native that is
far from it); either sign is a gradient. Covariates, per target, sequences only:

    I_long   max longer-normalised NW identity (core.data.identity, the pinned criterion) over
             the training corpus of the target's fold model (out-of-fold peptides + fold fragments;
             < 0.6 for peptides by construction; fragments are filtered against the fold's peptides
             at the same threshold)
    I_short  max shorter-normalised identity (containment) over the same corpus
    L_kmer   the longest exact substring shared with any training chain, in residues

Falsifier. Spearman rho between each covariate and `rmsd_arm` (production cache), one-sided in
the direction of HELP (rho <= -0.25, higher proximity lower RMSD), fold-clustered bootstrap CI
(clusters = the 5 pinned folds, 4,000 resamples), and a 500-draw label-permutation p-value;
Bonferroni over three covariates (alpha 0.017 each). A gradient EXISTS iff at least one
covariate has rho <= -0.25 with the fold CI excluding zero and the permutation p < 0.017. A
HARMFUL gradient (rho >= +0.25 by the same standard) is reported as such. Registered
expectation: no gradient, all |rho| within 0.15; then the deliverable is "a gradient of |rho|
>= 0.25 is excluded" and the clean bill for the 0.6 threshold.

Secondary, reported not decided on: the same for `rmsd_avg` and `shipped` (sel); the partial
correlation given chain length n (longer targets have higher I_short by chance, L15's
by-length table); and, if a gradient exists, the ORACLE mechanism check: the posterior's MAE
against the native (`s12/cache/disto_<pdb>.npz` expected vs native distances) against I_short.

## 2. Expected effect against the computed MDE

SE(rho) at n = 126 is about 1/sqrt(123) = 0.090; MDE = 2.8016 x 0.090 = 0.25 in rho, the
threshold above. Expected |rho| < 0.15.

## 3. Memory, time

`core.data.identity_many` with the composition prefilter over 787 + 6,003 chains per target
(lane I's audit: 126 x 786 in 15 s), the substring scan trivial; the cache's 126 records. Under
0.2 GB, about 5 minutes; one governed job; no projection, no native coordinate. No replication
seed (deterministic covariates); the permutation null is seeded by `s15.seed.stable_rng`.

## 4. Operator forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| corpus | the fold model's actual training corpus (`p_ladder.train_entries(fold)`, the same list `core.predict.train_fold` uses) | the whole database (includes same-fold chains the model never saw) |
| identity | the pinned NW identity, both normalisations, gap -1 | BLOSUM similarity (that is retrieval, a different question) |
| statistic | Spearman, one-sided, Bonferroni over 3 | a fitted router (S22/S23: routers are the wrong instrument at n = 126) |
| outcome | `rmsd_arm` (built chain) PRIMARY; `rmsd_avg`, `shipped` carried | the leaked-model envelope (already measured) |

## ADDENDUM 1 (2026-09-14 00:15) -- measured; one confound the file did not name; nothing above edited

- Jobs `w_recall_cov` (15 s, 0.292 GB) and `w_recall_endpoint` (45 s, 0.111 GB); ledger L90.
  No covariate reaches rho <= -0.25; I_long is -0.205 / -0.210 / -0.230 on arm / avg / shipped
  with the fold CI excluding zero and permutation p 0.010 / 0.014 / 0.004: 0.81 to 0.91x the MDE,
  suggestive, not measured. Power: |rho| >= 0.25 excluded.
- Confound not named above: the training corpus is also the retrieval library, so I_long is
  equally a retrieval-proximity covariate; the registered ORACLE mechanism check was gated on a
  gradient at the MDE and did not run. If this idea is taken up again, the mechanism check
  (posterior MAE vs I_short) should run unconditionally and a retrieval-only covariate (the pool
  best's BLOSUM sum) should be partialled out.
