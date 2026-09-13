# PROPOSAL B -- A LEARNED FOLDING MODEL AS THE PRIOR (lane P, Sprint 26)

Status: DRAFT 2026-09-13 09:35. B1 is final. B2 and B3 are filled in as their pre-registered runs
complete (`s26/PREREG_B2.md`, `s26/PREREG_B3.md`). Verdict form per the campaign prompt.

## What the proposal says

Replace or augment the shipped distance prior with the output of a large pretrained folding
model (ESMFold), on the argument that the prior's accuracy is the only steep lever the project
has measured (S24 L13: -2.15 A per unit gamma at the origin) and that a model trained on the
whole PDB should be a better prior than a 250k-parameter MLP trained on 787 peptides.

## What the repository already knows

- The prior's derivative is steep only along the native's own direction; a real operator at
  cos 0.5 travelling 25% is worth +0.024 A (S25 L12). Every re-reading of the existing posterior
  is closed (S25 L2/L6/L7/L12).
- The one part of a large pretrained model this box already consumes, ESM-2 650M, is worth
  -0.288 A [-0.484, -0.092] over one-hot on selection (S7-11), and its structure-supervised
  contact head is a filter, not a discriminator (S17 L23: +0.116 in-band, mostly compactness, no
  argmin gain). The per-target artefact behind S7-11 is lost (`s7/repr_tune.json`; lane P L11).
- The pool's error is 68% common-mode (S23 L9); the harmful coherent error is shared across
  predictor families and two thirds of it is reproduced by zero-information references (S19 L11):
  a better prior must be better in the one third that is target-specific.

## B1 -- can the folding model run here? (final; `s26/results/b1_feasibility.json`, ledger L13)

No, on three independent grounds. (a) No ESMFold or ESM-2 3B checkpoint is on disk; the download
is 2,771,653,574 + 5,678,116,398 bytes (HTTP HEAD, 2026-09-13; network reachable in under 10 s;
nothing downloaded). (b) `import esm.esmfold.v1.pretrained` raises `ModuleNotFoundError:
omegaconf` and `openfold` is absent; fair-esm's README says openfold needs nvcc and Python 3.9 or
lower; this box is Python 3.13 with CPU-only torch 2.13. (c) fair-esm loads the 3B language model
in fp16 and the 690M trunk in fp32 (esmfold.py:43-46): 8.76 GB resident before the load transient,
against 4.4 GB of campaign headroom (governor: 16.75 GB total, 64.8% used, ceiling 93%); fp16
everywhere is 7.4 GB and still does not fit. B1 stops. The 650M contact head is not a stand-in
folding model: it is a contact map and is already 13 of the shipped prior's 183 input columns.

## B2 -- the feasible-scale ladder (pending; rungs `conly` and `esm8m` of `s26/p_ladder.py`)

To be filled from `s26/results/p_ladder_report_conly_s0.json` and `_esm8m_s0.json`: contact head
alone as a prior input; the 8M language model against the 650M; the size axis read as a monotone
test (PREREG_B2). Built chain primary, point cloud and selection carried, ST.fmt verbatim.

## B3 -- where the pipeline gains over sequence-only, and is that set characterisable? (pending; `s26/results/p_b3.json`)

Per-target arms (built chains, persisted): pipeline 3.2148, sequence-only torsion predictor
3.7705, constant helix 4.0648; tors - arm +0.5557 (SE 0.1318, MDE 0.3694, 37W/89L), helix - arm
+0.8500 (SE 0.1494); FAIL18: tors 5.569 / arm 6.032 / helix 5.887 (ledger L14). The
characterisability test (nested ridge classifier of sign(d) against a 300-draw permutation null;
nested ridge regression of d with the fold-clustered MSE-reduction CI) is filled in from
`s26/results/p_b3.json` when `p_b3_run` completes.

## Verdict

Pending B2/B3. Per the coordinator's ruling (2026-09-13): if B2/B3 yield no native-free
characterisable set with a fold-clustered effect above MDE, the verdict is REPLACE, the
replacement is the trainability paper outline written by lane Q (`s26/PROPOSAL_B_REPLACEMENT.md`),
and the feasible-scale ladder result (conly, esm8m) folds into Proposal C as its second candidate.
