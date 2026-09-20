# S29 LANE M -- FINDINGS

Standing format: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN, then "what
damaged my own expectations" and "what I did not do and why".

---

## DEMONSTRATED

**M1. The harness is sound, on nine checks.** `s29/s29_M_harness.py`, job `m_harness_audit3`,
`s29/results/s29_M_harness_audit.json`, ledger S29-L9. `ca_rmsd` agrees with an independent Kabsch
to 3.675e-14 and forbids reflections; the native is the pinned universe `nat_ca`, bit-equal to
`peptide_db.npz` and to a fresh parse of the deposit; no benchmark artefact is read anywhere on the
tuning path under an in-process poison stronger than a rename; `pinned_folds` are the pinned folds
on all 126; the 3.2126 / 3.0483 anchors reproduce from `chain_rows.jsonl :: DIS` and from a fresh
re-projection on 6 targets at 0.000e+00 A; `ST.compare`'s MDE is the literal 2.8016 x SE and its
fold CI is a cluster bootstrap over the 5 folds; the 126 are 126 distinct pdb ids, 126 distinct
sequences, lengths 9 to 16.

**M2. The s12 distogram cache is NOT bit-identical to a fresh recomputation, and it changes
nothing that matters.** Over 126 targets: max |prob| 1.97e-06, max |risk| 1.34e-04 (relative
5.63e-06), max |score| 1.91e-06 against a score sd ~1.25. A fresh recomputation is exactly
deterministic in-process (0.0e+00), so this was introduced when the cache was written. The measured
consequence: the 500-candidate score order differs on **2 of 126**, the top-75 set on **0 of 126**,
and the emitted point-cloud RMSD by **0.00e+00 A on 126 of 126**. Actionable: an arm that reads the
score order *below* the top-75 cut should recompute rather than read the cache.

**M3. Leave-fold-out constrains only 9.5% of the distogram's training data.** `fold_fragments`
removes 0 to 15 of 6,003 fragments per fold, so 90.5% of every fold model's ~6,630 training chains
are shared across all five folds. "Leave-fold-out" is a statement about the 787 peptides.

**M4. The production anchor never passes through the quantum stage.** `s27/run_vqe_chain.py:124-130`
(`--chain`) is the classical tie-safe top-75; `arm_vqe` is called only in `--vqe` mode; production
`Config.quantum = False` (`core/pipeline.py:179`) and the production cache record carries
`quantum: null`, `n_top: 75`. Stated because two S29 lanes will otherwise design against a spine
that does not exist as described.

**M5. 17 of 31 recorded modelling choices have never been tested.** `s29/CONVENIENCE_CHOICES.md`.
The five with the largest reach are the 17 bins and their edges, the soft-bin sigma and the unused
separation weights, the retrieval key's untested class (coverage-preserving rather than
quality-ranked), the top-128 truncation's unmeasured ORACLE ceiling, and the medoid frame. `T = 0.5`
in the S25/S27 harness cell has no recorded criterion anywhere.

## ORACLE DIAGNOSTIC

**M6. Leave-fold-out is empirically real, and the distogram memorises its training peptides by 8x.**
Mean NLL of the true distance bins: **3.28777** for the deployed (own-fold, held-out) model against
**1.21275** for the four models that trained on that fold; delta +2.075, SE 0.152, **4.88x MDE**,
correct sign on all five folds. This is the strong form of the leave-fold-out check (a mis-trained
checkpoint would show the opposite sign) and it is simultaneously a statement about the model: the
training-set likelihood is e^2.075 = 8.0x the held-out likelihood. Any in-sample diagnostic on the
corpus -- calibration, sharpness, MAE, a fitted residual correction -- therefore describes
memorisation and must be recomputed out-of-fold. Reads the native; tunes nothing.

## OPEN

**M7. The top-128 ORACLE ceiling is not in the record** and bounds every production quantum arm
(one line: `min(oracle_rr[order[:128]])` per target, meaned). The record has the top-75 (2.3062) and
the pool (1.7108) ceilings but not this one.

**M8. The medoid frame is the one readout choice with no measurement** -- and S23 L9's common-mode
decomposition, which the sprint's leading hypothesis rests on, is *defined* in that frame.
A generalised-Procrustes frame is a one-minute test.

**M9. Every alternative retrieval key ever tested was another native-free quality score.** S17 L12
names the untested class in its own conclusion (a diversity- or coverage-preserving shortlist) and
nobody built one. This is the only untested route to the 1.313 A universe ceiling.

## What damaged my own expectations

1. **I expected the cached posterior to be bit-identical and it is not.** I wrote check 8 with a
   pass criterion of exact equality, it failed, and the honest repair was not to loosen the
   threshold but to *measure the deployed consequence* -- which turned a binary into a statement
   with a number (order 2/126, set 0/126, RMSD 0/126). The first version of that check would have
   been reported as a harness failure and would have been wrong.
2. **I expected "leave-fold-out" to mean the training set changes substantially between folds.** It
   changes by ~30 peptides out of ~6,630 chains. That reframes every fold-clustered CI in the record
   as, if anything, conservative -- and it is an independent reason (beside S25 L17's rank ladder)
   why per-fold quantities carry so little.
3. **I expected the harness audit to find something wrong.** The charter invites "evaluation-harness
   failure" as a hypothesis and the instrument is 25 sprints old. Nine checks found one float32-scale
   cache discrepancy that moves nothing. The instrument is in better shape than the result it
   measures, which is the opposite of the failure mode I was looking for.
4. **I did not expect the memorisation gap to be 8x.** I added check 9 as a leakage detector and it
   returned a model-capacity finding.

## What I did not do, and why

- **I did not rename `results/benchmark_manifest.json`,** as the brief's poison test suggested.
  Contract rule 2 forbids moving a pinned artefact, and an in-process interception of `open`,
  `os.path.exists` and `np.load` is strictly stronger (it catches a cached read and an
  existence probe, which a rename does not). The pinned file was never touched.
- **I re-projected only 6 targets freshly**, the brief's number. The projection is the expensive
  stage (~3 s per target); the cold-distogram half of the path was checked on all 126.
- **I did not re-derive the set-equality theorem or the VQE harness.** S25 verified them on 2,592
  adversarial cells and S24 on 3,888; re-running them would be a cosmetic rerun (contract rule 10).
  What this lane added instead is the statement that the anchor does not pass through them.
- **I did not run `s12/instrument.py:266 selfcheck`**, which asserts pool_best 1.7108 and
  top75_best 2.3062 and re-derives FAIL18: it reads the production cache and is a pool-integrity
  check, not a harness check. A lane that changes the pool must run it.
- **I ran no experiment beyond the audit's checks**, per the brief's discipline section. Every
  "cheap decisive test" in `s29/CONVENIENCE_CHOICES.md` is written as a proposal with its falsifier
  and its cost, not run.
