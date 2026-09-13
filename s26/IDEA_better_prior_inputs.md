# IDEA_better_prior_inputs -- WHICH INPUTS THE MACHINE CAN COMPUTE MIGHT MOVE THE GAMMA LADDER (lane P, Sprint 26)

Open item 1 of the state brief: the project measured what a better prior would BUY (S24 L13:
-2.15 A per unit gamma at the origin, gamma = 0.0225 reaches 3.0 A) and never whether one is
OBTAINABLE. This entry lists every input a distogram could consume on this box, with the record's
verdict, and pre-registers the two that are open.

## Hypothesis

A per-pair input the shipped MLP does not see carries target-specific distance information that
the ESM-2 pca32 block does not, and a distogram trained on it moves the endpoint. The two open
candidates are ranked by how much of the input is PAIRWISE, because the shipped prior's pair
channel is 13 scalars out of 183 columns while its failure mode is a pairwise, coherent,
separation-profile error (S19 L14).

| input | pairwise? | status on the record | verdict |
|---|---|---|---|
| ESM-2 650M raw attention-head maps (660 heads, APC-corrected), top-k by contact-head weight, fed per pair | yes | the contact head is one logistic regression over them (`esm_features.py` docstring); only that scalar and 12 derived columns reach the prior; S17 L23 measured the scalar as a RANKER (+0.116 in-band, mostly compactness), never as a prior INPUT | **OPEN, ranked 1** |
| the K = 500 retrieval pool's per-pair distance histogram as an INPUT at training time (retrieval-augmented prior) | yes | as an OUTPUT-side mixture it is typicality: `gain 0` selects 3.734 vs 3.454 (S7-3), the trained predictor beats typicality by 0.024 (S7-8), every move toward the pool is worse (S19 L14); as a training-time INPUT with the native as label it is untested | **OPEN, ranked 2** (foreshadowed null) |
| ESM-2 650M contact head ALONE (no embedding) | yes | never isolated as a prior input | rung `conly`, PREREG_B2 |
| ESM-2 higher-dimensional / raw / per-fold PCA | no (per residue) | S7-11: pca128 / raw indistinguishable from pca32 on selection | rungs pca32f/pca128/raw, PREREG_C2 |
| a larger language model (ESM-2 3B, ESMFold trunk) | partly | infeasible: 5.68 GB fp16 checkpoint, >= 5.7 GB resident vs 4.4 GB headroom (`s26/results/b1_feasibility.json`) | closed on this box; the size axis is esm8m vs 650M (PREREG_B2) |
| MSA / evolutionary features | yes | S9-7: "exists, is retrievable, and does not help" | closed |
| chemical-shift-derived torsions | no | S14: 54/126 coverage; ORACLE-perfect torsions still 2.021 A | closed |
| de-noised / ensemble-mean labels | -- | S8-14: worth 0.044 A by arithmetic | closed |
| more fragment training data | -- | S7-2: monotonically worse (distribution shift) | closed |
| more PEPTIDE training data | -- | the 9-16-mer world supply is spent (memory `no-fresh-benchmark-exists`; corpus census S24 L8: 410 permitted peptides) | closed |
| a physics-refined label or an energy-weighted pool histogram | yes | never as a prior partner; predicted null by S8-14 and S24 L16 | PREREG_C3 draft, PH lane |

## Why the record does not already close the two open ones

Attention maps: every ESM result in the record consumed either the 32-d per-residue projection or
the contact head's single scalar per pair (S7-11, S17 L23, the shipped `features`). No arm has
given the distogram the head maps themselves; S17 L23's "mostly compactness" verdict is about the
logistic summary, whose weights were fitted on protein contacts at |i-j| >= 6 and cannot represent
what a peptide's heads know at |i-j| = 2-5, where the shipped posterior is 2x over-confident
(S25 L1). Retrieval-augmented input: S7-3/S7-8/S19 L14 all put the pool's statistics on the OUTPUT
side (re-ranking or replacing the prediction); a model trained with the pool's histogram as an
INPUT and the native as the LABEL can learn the pool's systematic error, which is exactly the
68% common mode (S23 L9) that no output-side operator can see. The record's prediction is still a
null (the label is 90% model-1 noise-free but the pool's error is shared with the training
fragments, S19 L11: retrieval-aligned 0.784 of the ceiling), and it is stated as such.

## Exact falsifier

PREREG_C2's, applied to two new rungs: `attn` and `ragp`. Built chain, effect beyond its own
MDE, fold-clustered CI excluding zero, 5/5 folds; gam_eff with cos beside it.

## Expected effect vs computed MDE

MDE on the built chain 0.09-0.20 A (PREREG_C2 section 4). Expected: `attn` -0.05 to +0.05
(the head maps' information beyond the contact scalar is unmeasured; S17 L23's ESM-specific
increment is +0.050 in Spearman, i.e. small); `ragp` 0.00 to +0.05. Both are therefore likely
"not measured" at n = 126; the deliverable is the power statement and gam_eff.

## Memory and agent-hours (probe required before running)

`attn`: one forward pass of the 650M model with `need_head_weights=True` over 6,916 sequences.
The checkpoint is 2.60 GB on disk; a probe (1 sequence under jobrun, tag ESM) must measure the
resident model plus the (33 x 20 x n x n) attention tensor before anything else runs; expected
3.0-3.5 GB, i.e. it needs the box to itself and the coordinator's OK. The reduced cache (top-32
head maps per sequence, float16) is ~170 MB. Training: the lean trainer with a 32-column pair
block, ~1.0 GB. 4 agent-hours. `ragp`: K = 500 BLOSUM retrieval for 6,790 training sequences from
the fold-safe library (~1-2 h CPU, one job, < 1 GB), then a 17-column pair block; 3 agent-hours.
