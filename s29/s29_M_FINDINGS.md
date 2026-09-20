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

**M10. F1 -- the selection functional is not the lever, and the log score's ordering is worth
nothing beyond the swap rate.** S29-L51, `s29/results/s29_M_F1_rows.jsonl` (126 targets x 6 arms).
Swapping ONLY the functional (L1 Bayes risk -> the log score of the same 17-bin posterior), pool,
top-75, readout and projection all held at production: **LOG - PROD +0.0822 (0.55x MDE, 3/5 folds,
NOT MEASURED)**. LOG exchanges 20.6 of 75 members (overlap 0.726); exchanging the same number at
random costs +0.0621, so **LOG - SWAPCTL = +0.0202 at 0.14x MDE**, fold CI [-0.1020, +0.1622],
57W/69L, power 0.07. The functional's ordering buys nothing at the endpoint. C9 in
`s29/CONVENIENCE_CHOICES.md` is closed in the direction the record predicted.

**M11. F1's permutation control proves the log score's pair information is REAL, which is what makes
M10 a strong null rather than a tie.** Scoring pair p's distance against pair perm(p)'s posterior
costs **+0.7295, 2.29x MDE, fold CI [+0.5658, +0.9637], 5/5 folds, power 1.00**, with the top-75
overlap collapsing to 0.199. A working functional whose ordering is simply not better than the
incumbent's.

**M12. Contraction is not the lever, measured a third independent way -- from the selection stage.**
LOG's emitted cloud is **+0.0593 A LESS contracted** (1.34x MDE, fold CI [+0.0203, +0.1039], 4/5
folds) and its endpoint is +0.082 WORSE; the matched functional control L2RISK contracts **MORE**
(2.9207 vs 2.9614) and is **also** worse (+0.0440). Two functionals move contraction in opposite
directions and both move the endpoint the same way. Independent of S29-L10 (the shipped cost's
descent expands) and S29-L12 (the contraction is Jensen on the average).

**M13. F2 -- the shell-profile class closes on a MEASURED SUPPLY GAP, not on an absent
instantiation.** S29-L53. The fitted native-free profile ratio's displacement cosine is **0.0895**,
below lane T's 0.140 line, and it is **beaten by its own zero-information shrink twin**
(RATIO - RSHRINK = **-0.0319**, fold CI [-0.0628, -0.0057], 4/5 folds); pure typicality with no fit
reaches 0.122. The registered stopping rule fired and **no projections were spent**. Supply side:
LFO `corr(r_hat, r_true) = +0.313` against the **incumbent distogram ratio's +0.366** -- the fit
predicts the true profile ratio *worse* than the shipped prior already does, a fifth instance of
`MAE does not price selected RMSD`, now at the level of the profile's ratio.

**M14. The F2 reproduction gate reproduces a seventeen-sprint-old pair of numbers to four decimals,
and the residual 0.02 A was chased to a line of code.** Under S12's own convention: PROD 3.0784 vs
3.078 (+0.0004), ORACLE 2.4023 vs 2.402 (+0.0003), gap 0.6761 vs 0.676. The offset in the weighted
run is `s12/obj_common.py:93-97 score_l1` being called from `s12/obj_profile.py:156` with `w = None`
-- S12's profile arms are uniformly weighted. The S12 comparison is exact rather than approximate.

## ORACLE DIAGNOSTIC

**M6. Leave-fold-out is empirically real, and the distogram memorises its training peptides by 8x.**
Mean NLL of the true distance bins: **3.28777** for the deployed (own-fold, held-out) model against
**1.21275** for the four models that trained on that fold; delta +2.075, SE 0.152, **4.88x MDE**,
correct sign on all five folds. This is the strong form of the leave-fold-out check (a mis-trained
checkpoint would show the opposite sign) and it is simultaneously a statement about the model: the
training-set likelihood is e^2.075 = 8.0x the held-out likelihood. Any in-sample diagnostic on the
corpus -- calibration, sharpness, MAE, a fitted residual correction -- therefore describes
memorisation and must be recomputed out-of-fold. Reads the native; tunes nothing.

**M15. The deployed score is AT CHANCE for locating the best member of its own top-128, and on the
typical target slightly worse than chance.** S29-L53, ORACLE diagnostic.
`bits_delivered = 7 - log2(rank of the ORACLE-best member of the FIXED production top-128 under the
arm's own score)`; the matched null is a uniform random ranking, exactly
`7 - (1/128) sum_{r=1..128} log2 r = 1.4050`. Measured: **shipped 1.4415, effect +0.0366 at 0.10x
MDE**, fold CI [-0.2251, +0.3110], 2/5 folds. **The median is the honest summary and it is
negative: median paired difference -0.575 bits, below the random baseline on 82 of 126 targets,
median rank 72 of 128 against a chance median of 64.5.** The positive mean is carried by a handful
of targets. Even ORACLE profile knowledge supplies only **3.202** of the 7, so the readout's 7-bit
gap is not a profile-shaped gap either. Two readability notes: `ST.compare`'s W/L and
"better/worse" labels INVERT for bits (higher is better), as lane D flagged in S29-L6; and per the
S29-L44 correction, **bits stay in bits** -- nothing here converts a bit into an Angstrom through
the mean-consuming readout.

**M16. F1's arm emits PRODUCTION's answer, not a new one.** The S24-L2/L3 parallel-bias check,
ORACLE and post hoc: the cosine between LOG's and PROD's error vectors against the native, rigid
body removed, is **0.924 mean / 0.969 median** against a random reference of 0.177 at 3n-6 dof. One
number that explains the endpoint result without any statistics.

## REFUTED (by my own measurement, before it reached a claim)

**M17. My own registered prior for F2 was wrong in the flattering direction, on BOTH inputs.**
Addendum 1 predicted a fitted cosine of 0.16-0.24 ("genuinely close to the line") from 0.24-0.37
(native-free compactness proxies) x 0.66 (an ORACLE cosine **inferred by inverting the bound's
identity on the RMSD ratio**). Measured: the ORACLE displacement cosine is **0.483**, not 0.66 --
the inverted identity ran 37% high -- and the fitted arm reached **0.0895**, not 0.16. The shrink
twin is what caught it. **A quantity derived from an identity rather than measured, and wrong in
the direction its author wanted, is the same failure mode lane T hit with the sign formula.**

**M18. The contraction mechanism F1 was assigned to test was withdrawn before F1 ran, and F1's own
data then confirmed the withdrawal** (M12). The prereg's addendum 1 removed the claim on S29-L10
and S29-L12; had it not, F1's +0.059 A less-contracted cloud with a worse endpoint would have been
read as a mechanism that "worked".

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

**M19. F2 does NOT close assumption B2 in general.** It tested ONE quantity (the per-separation mean
distance profile), one parameterisation (the ratio to the pool's own profile) and one model class
(LFO ridge over the named feature blocks). What it removes is the specific hope that the
lowest-dimensional named ORACLE quantity in the repository had a native-free twin. B2's general case
is live -- and more so since lane T's compactness measurement (S29-L50) refuted the coordinator's
expectation and left the general question open. **F2's specific case is closed; B2 is not.** The two
must be kept apart in the report.

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
5. **I expected F2 to be close to the line and it was not close** (M17). Worse, I had computed the
   quantity I was most confident about -- the ORACLE cosine -- by *inverting an identity* rather
   than measuring it, and it came out 37% high. I now treat any number produced by rearranging an
   identity as an estimate that must be measured before it is used to set an expectation.
6. **Twice I wrote a control that was the identity by algebra.** I correctly refused the
   coordinator's monotone re-ranking control (inert through a top-m readout) and then, in the same
   experiment, wrote a permutation control that permuted per-pair values *after* the gather -- which
   for a SUM aggregator is also the identity, and it reproduced LOG to every printed digit. Spotting
   the failure mode in someone else's proposal did not stop me committing it in my own code. The
   12-target probe caught it; a 126-target run would have "confirmed" a null control.
7. **The record predicted F1's endpoint to a hundredth of an Angstrom.** S7's pure-NLL arm at
   +0.093 through argmin; measured +0.082 through the top-75 average and the built chain. I had
   registered "null to worse, +0.00 to +0.10" and still half-expected the different readout to
   matter. It did not.

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
- **I ran no experiment beyond the audit's checks** during the deliverables stage, per the brief's
  discipline section. Every "cheap decisive test" in `s29/CONVENIENCE_CHOICES.md` is written as a
  proposal with its falsifier and its cost, not run. F1 and F2 were later assignments from the
  coordinator, each pre-registered before any number.
- **I did not run F2's deployable projection arm.** The pre-registered stopping rule (addendum 1,
  A1.3) required the measured cosine to clear 0.140 AND beat its shrink twin; it did neither, so no
  box was spent. Reporting the closure was the registered outcome, not a fallback.
- **I did not chase F2 further** -- a second model class, a different structural quantity, a
  per-target sign. Each is a new experiment needing its own pre-registration, and M19 states the
  scope rather than implying the class is exhausted.
- **I did not re-run F1 on a second seed.** No seed enters any F1 arm: the scores are deterministic
  and the only RNG is the stable per-target tie key. The prereg's replication clause applied to a
  positive, and there was none.

## An operational defect of mine, recorded rather than quietly fixed

While F1's run was queued at the 8-job cap I launched the same named job four times (a backgrounded
attempt, a chained waiter, a retry); **`s26/jobrun.py` does not deduplicate by `--name`**, so all
four eventually got slots and ran the same work concurrently, producing 220 duplicate (target, arm)
rows and burning box time other lanes were queued for. Three were killed and the registered pid
kept. No science was affected -- every arm is deterministic given the stable per-target RNG and the
rows dedupe last-wins to identical values -- and F1's analysis is unchanged. **A job name is not a
lock**: check `s26/jobs/<name>.json` and the live process list before launching, and never chain a
launcher behind a waiter that may itself be retried.
