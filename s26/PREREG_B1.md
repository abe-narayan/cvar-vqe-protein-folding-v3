# PREREG_B1 -- CAN A LEARNED FOLDING MODEL RUN ON THIS BOX? (lane P, Sprint 26)

Written 2026-09-13. B1 has two halves: a FEASIBILITY MEASUREMENT (no endpoint; taken during the
reading phase as the lane brief instructs, and recorded in `s26/results/b1_feasibility.json`) and,
conditional on feasibility, an ENDPOINT EXPERIMENT (ESMFold's distogram as the prior). The
endpoint half is pre-registered here and has not run.

## Hypothesis

**H_B1 (feasibility).** ESMFold v1 (`esm.pretrained.esmfold_v1`, fair-esm 2.0.0) can be imported,
loaded and run on the 6,790 training sequences and 126 targets inside the campaign's 93% RAM band
(about 4.4 GB of headroom for all lanes, `s26/governor_state.json`), so that a folding-model prior
can be measured on the pipeline.

**H_B1 (endpoint, conditional).** ESMFold's own distogram head, binned onto the shipped 17
`CENTRES`, used as the prior in place of the shipped MLP posterior, beats the shipped posterior on
the built chain by more than its MDE, fold-clustered, 5/5 folds.

## Exact falsifier

Feasibility fails if ANY of: (a) the import fails on a missing dependency that cannot be installed
without CUDA or a Python downgrade; (b) the checkpoints' fp16 resident size exceeds the headroom;
(c) the download exceeds 100 MB without the coordinator's OK. If feasibility fails, B1 STOPS and
the endpoint half is not run; `s26/results/b1_feasibility.json` is the deliverable and Proposal B
is routed to B2/B3 and to `s26/PROPOSAL_B_REPLACEMENT.md`.

Endpoint falsifier: identical to PREREG_C2's (built chain, > MDE, fold CI excludes zero, 5/5 folds).

## Measured 2026-09-13 (artefact `s26/results/b1_feasibility.json`)

- (a) weights: only `esm2_t33_650M_UR50D.pt` (2.60 GB) and `esm2_t6_8M_UR50D.pt` (30 MB) are in
  `~/.cache/torch/hub/checkpoints`; no ESMFold or 3B checkpoint.
- (b) import: `esm.pretrained.esmfold_v1` exists; `import esm.esmfold.v1.pretrained` raises
  `ModuleNotFoundError: No module named 'omegaconf'` (esmfold.py:10); `openfold` is also absent
  (esmfold.py:11-13, trunk.py:11 `StructureModule`). The fair-esm README's install path for
  openfold "requires nvcc" and "python <= 3.9"; this box is Python 3.13 with CPU-only torch 2.13.
- (c) size: README -- esmfold_v1 is 48 (+36) layers, 690M (+3B) parameters. HTTP HEAD:
  `esmfold_3B_v1.pt` 2,771,653,574 bytes, `esm2_t36_3B_UR50D.pt` 5,678,116,398 bytes (8.45 GB
  total). Parameters alone: fp32 14.8 GB, fp16 7.4 GB; as fair-esm loads it (esmfold.py:46
  `self.esm.half()`, trunk fp32): 8.76 GB resident minimum, before the torch.load transient.
- (d) network: `curl -sI https://dl.fbaipublicfiles.com` 200 OK in under 10 s; nothing downloaded.
- (e) band: 16.75 GB total, 64.8% used at 00:19, ceiling 93% = 15.58 GB, headroom 4.4 GB.
  8.76 GB does not fit; 7.4 GB (fp16 everything) does not fit.

**Feasibility verdict: NO on three independent grounds. B1 stops here.**

## The cheaper stand-in the brief asked about

ESM-2 650M's contact head is already exposed (`esm_features.compute` with `return_contacts=True`;
`core.data.esm_contacts`), already cached for every training sequence and target
(`s7/repr_cache/esmcon.npz`), and already INSIDE the shipped prior (13 of the 183 input columns).
S17 L23 measured it as a ranker: in-band Spearman +0.116 [+0.047, +0.184] at top-75, +0.050
[+0.005, +0.097] over its ESM-free twin, mostly compactness, no argmin gain, and a shortlist
significantly worse than matched random. It is a contact map, not a folding model. Its value AS A
PRIOR INPUT on its own (without the embedding block) is not on the record and is carried as rung
`conly` in `s26/p_ladder.py` under PREREG_B2.

## Comparison arm, basis, MDE, memory, agent-hours

Endpoint half (not run): the shipped posterior through `s26/p_ladder.py`'s path; built chain
primary; MDE as PREREG_C2 section 4 (0.09-0.20 A on the built chain). Memory: the measurement
above (0.6 GB of probes); the endpoint half would need >= 8.76 GB and cannot be estimated from a
probe on this box. Agent-hours spent: 1.0 (measurement and this file).

## Rule 0 forks for the endpoint half (kept so the prereg is complete if the box ever changes)

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| functional | ESMFold's distogram logits (64 bins, 2.3125-21.6875 A) re-binned by mass onto the 17 shipped centres, then the shipped Bayes risk from a genuine `core.predict.Distogram` | a Dirac prior on ESMFold's predicted CA coordinates (collapses the posterior; S24 L13-A crossover) |
| basis | built chain primary; cloud and selection carried | ESMFold's own pLDDT/pTM as a metric |
| readout | shipped top-75 average, m = 75 | ESMFold's coordinates as a candidate (that is the generation lane, closed S24) |
| normalisation | ESMFold run once per sequence at chunk_size None; no recycling sweep | tuning num_recycles per target |
| null | the shipped posterior; the shipped posterior mixed with ESMFold's at lam = 0 | a scrambled-sequence ESMFold run (a different question) |
| label | built-chain CA-RMSD | TM-score, pLDDT |
