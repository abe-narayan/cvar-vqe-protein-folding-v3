# S29 CONVENIENCE CHOICES -- the modelling decisions made for convenience, and what tested them
(lane M, 2026-09-19; companion to `s29/DATAPATH.md`, same line-number basis, HEAD `b1618a6a`)

The charter's item 7 asks for every place a modelling choice was made **for convenience rather than
for a stated scientific reason**. This file lists them with file and line, what the alternatives
were, whether the record ever tested the choice (with the ledger line), and where it did not, a
**cheap decisive test**.

Three classification notes, because "convenience" is not the same as "arbitrary":

- **TESTED** -- the record contains a measurement of this choice against a named alternative on the
  126-target instrument (or a declared subset), and the entry is cited. A tested choice is not a
  redesign target unless the formulation changes.
- **PINNED** -- the value is load-bearing for reproducibility (folds, clusters, the pool order) and
  changing it invalidates artefacts rather than testing a hypothesis. These are listed because a new
  architecture must either inherit or re-derive them, not because they should be tuned.
- **UNTESTED** -- no measurement exists. These are the redesign targets. **17 of the 31 entries
  below are UNTESTED**, and the five with the largest reach are C1 (the 17 bins and their edges),
  C2 (the soft-bin sigma and the loss), C6 (BLOSUM as the retrieval key, tested once at the wrong
  altitude), C13 (the top-128 / 2^n truncation) and C24 (the medoid frame).

Sprint-numbering convention here: `S25 L6` = `s25/LEDGER.md :: ## L6`; `S28-L21` = `s27/LEDGER.md ::
## S28-L21`; `S10-2` = `docs/FINDINGS.md :: ## S10-2`.

---

## A. The corpus, the folds, the instrument

### C1. 17 distance bins and their edges -- **UNTESTED, and the largest untested choice on the path**
`core/predict.py:54-59`. `BIN_EDGES` = 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0,
12.5, 14.0, 16.0, 19.0, 23.0; `CENTRES` prepends 4.0 (= first edge - 0.5) and appends 25.0 (= last
edge + 2.0). The docstring's reason is qualitative ("fine where CA-CA distances are structurally
informative ... coarse beyond"); no experiment chose the count, the spacing or the two invented
outer centres.
**Alternatives:** more bins (AlphaFold's distogram uses 64 over 2-22 A); uniform spacing; a
mixture-density or a continuous (Gaussian / von Mises-Fisher) head; a quantile binning of the
corpus's own distance distribution; bins per separation shell.
**What the record tested:** only the *consumption* of the existing 17 bins, never the binning
itself. S25 L6 showed the L1 consumer's effective per-pair target takes exactly 17 values with a
largest adjacent gap of 4.0 A -- "up to 2.0 A of pure location error before any estimation error" --
and S25 L12 measured that removing that quantisation exactly (the DEQUANT family, mean-preserving,
verified by `ST.achieved` at 1.000000) is worth **0.003 A with full leakage**. That is a statement
about the consumer, not the representation: DEQUANT smooths a trained 17-bin posterior; it does not
ask what a 64-bin or continuous head would have learned.
**Cheap decisive test:** retrain one fold's MLP with (a) 34 bins on the same range and (b) a
2-parameter continuous head (mean, log sd) under the same features, epochs and seed, and compare on
the *selection* metric (top-75 point cloud) on that fold's 23-30 targets, not on MAE (`MAE does not
price selected RMSD`, memory; S7-3). Cost: 3 x ~10 min training + a 126-target scoring pass at
seconds per target. The falsifier: if the 34-bin arm's selection RMSD is within 0.7x MDE of the
17-bin arm on the held-out fold, the binning is not the limit and C1 closes.

### C2. The soft-binned loss: Gaussian-in-bin-index smoothing at sigma = 0.6, and no per-pair weight -- **UNTESTED**
`core/predict.py:269` (`smooth: float = 0.6`, the default, never overridden by `train_fold`) and
`core/predict.py:276-278` (`sample_w = None` always, so every pair has weight 1).
`separation_weights` (line 208, five modes) exists and is never called on the deployed path; the
checkpoints `fold*_esm_frag_sw-lin.pt` / `sw-none.pt` on disk show it was tried at some point, but
no ledger entry prices it.
**Alternatives:** hard one-hot targets; sigma chosen by validation; smoothing in *distance* rather
than in bin index (the bins are irregular, so a fixed sigma in index units is a different physical
width in every part of the range -- 0.5 A at the short end, 4.0 A at the long end); the separation
weighting the module already implements.
**What the record tested:** nothing on sigma. The docstring at lines 232-237 states the direction
("hard one-hot targets over 17 narrow bins make the model confident about a boundary it cannot
resolve") without a measurement.
**Cheap decisive test:** two extra one-fold trainings at sigma in {0.0, 1.2} and one with
`sample_w = separation_weights(sep, "lin")`, scored on selection on that fold as in C1. Note that
sigma interacts with the score's `1/(sd+0.5)` weight (C5) -- report the realised `sd` distribution
beside the RMSD.

### C3. `IDENTITY_THRESHOLD = 0.6` and the longer-length normalisation -- **TESTED (and a declared defect)**
`core/data.py:404`, `core/data.py:182 identity(norm="longer")`.
**Record:** ARCHITECTURE section 2.1 and `core/data.py:182-200` declare the leak (an 11-mer verbatim
inside a 21-mer scores 0.52 and lands in another fold); 4/126 dev and 2/60 benchmark targets carry a
verbatim self-copy in their own fold model's training set; priced at +0.0004 A dev (S10-4). S26 lane
I measured the corrected ("shorter") normalisation **in memory** and found 0.6 is AT THE NULL for
peptides as a clustering criterion (a real dev sequence passes against 0.5% of members, a shuffled
one against 0.3%; single linkage collapses 470 clusters to 166) -- memory
`containment-threshold-is-at-the-null-for-peptides`. Verdict on record: use the verbatim-substring
test as a per-pair leak test, not a replacement threshold. **PINNED**: `peptide_clusters.json` /
`peptide_folds.json` must never be re-derived (memory `benchmark-and-folds-must-be-pinned`:
correcting the clustering silently moved 13 benchmark targets).
**Open, cheaply:** report the verbatim-substring flag per target beside every S29 result, so any
positive can be re-read on the 122 clean targets. `core/data.py` already ships the flag dark.

### C4. The fragment bank is shared by all five folds -- **UNTESTED, and it changes what "leave-fold-out" means**
`core/data.py:794 fold_fragments` drops a fragment only if it is >= 0.6 identical to a held-out
*peptide*. Measured 2026-09-19: the bank keeps 5,988-6,003 of 6,003 per fold, i.e. **0 to 15
fragments removed**. The distogram's training set is 90.5% fold-shared fragments (6,003 of ~6,630
chains, `s29/DATAPATH.md` stage 3).
**Consequence:** "leave-fold-out" constrains 9.5% of the training data. The five fold models differ
by ~30 peptides each and are therefore much more similar to one another than the phrase suggests --
which is an independent reason (beside S25 L17's rank ladder) to expect per-fold quantities to carry
little information. No ledger entry states this.
**Cheap decisive test:** compute the pairwise agreement of the five checkpoints' predictions on a
fixed sequence set (KL between fold f's and fold g's posteriors on the 126 targets, seconds). If the
cross-fold KL is at the scale of the within-fold prediction entropy, the leave-fold-out structure is
nearly vacuous and every fold-clustered CI in the record is, if anything, conservative; if it is
small, the honest statement is that the five models are near-copies.

### C5. `w_p = 1 / (sd_p + 0.5)`, the +0.5 and the exponent -- **PARTLY TESTED; the +0.5 is UNTESTED**
`core/predict.py:416`, with `_score_weights()` (388) returning `ones, 1.0` because
`score_weights.json` is deliberately absent (the fitted shell weights were RETIRED: "does not
reproduce under corrected conditions ... loses on the benchmark", lines 378-383).
**Tested:** the exponent. S18 L14 swept `sd^-p` for p in {0, 0.5, 1, 1.5, 2, 3, 4} at n = 126 and
found the response monotone up to the deployed value then flat (p = 3 is -0.026 [-0.108, +0.056],
3/5 folds, inside the MDE); S18 L10 measured deleting the weights at **+0.153 [+0.066, +0.247]** and
permuting them at +0.091. NOTE the basis mismatch: S18's arms are `sd^-2` on the s18 refinement
instrument, while the deployed score is `(sd + 0.5)^-1` -- these are different functional forms and
the record's "the deployed exponent is at the optimum" is about the former.
**Untested:** the additive 0.5 A regulariser. It is a pure convenience (it prevents division by a
small sd) and it sets the scale at which confidence stops mattering; sd on this posterior runs
~0.7-4 A, so +0.5 is a 12-70% shift of the weight.
**Cheap decisive test:** re-score the 126 pools with `w = 1/(sd + c)` for c in {0, 0.25, 1.0, 2.0}
(no retraining -- the risk table rebuild is milliseconds per target) and report the top-75 point
cloud with `ST.compare` and `ST.best_of_k_within`. This is 4 x 126 target-scorings, minutes.

### C6. BLOSUM62 sum as the retrieval key, K = 500 -- **TESTED at the wrong altitude; the KEY is the open half**
`core/data.py:155 blosum_similarity` + `167 top_k` (stable argsort); `s8/generate.py:204-205`;
K = 500 (`s7/audit.py:45 K_CLAIM`, `s12/instrument.py:34`).
**Tested (K):** S17 L1 (K = 500 is a truncation artefact, not a ceiling: 1.711 A at K = 500 against
1.313 A over the full universe), S17 L12 (widening K makes the *answer* worse -- the shortlist
ceiling degrades from 2.104 at K = 75 to 2.572 at full, +0.468 [+0.304, +0.642], because extra
windows displace near-native members), S17 L19's matched null (the ceiling gain from widening is
below the order-statistic floor: observed minus null -0.075 [-0.132, -0.022]).
**Tested (key), once:** `s8/generate.py`'s RESULT 2 measured BLOSUM's own rank correlation with
true CA-RMSD over the universe at **+0.066**, and every alternative native-free key (Ramachandran
plausibility +0.307, rg agreement +0.321, typicality +0.410, the shipped score itself +0.587)
selected **worse** by +0.11 to +0.37 A, because a real key's error is correlated with what makes a
candidate good. Memory `real-fragments-beat-the-lattice` adds: BLOSUM wins at every K on a matched
test.
**Still open:** every alternative tested was a *re-ranking of the same universe by another
native-free quality score*. Nobody has tested a key that is **not a quality score at all** -- e.g.
structural-profile matching (predicted secondary structure or the predicted distogram itself as the
retrieval key), or a coverage/diversity-preserving shortlist, which S17 L12's own conclusion names
as the required class ("a diversity- or coverage-preserving shortlist -- a Problem C question").
**Cheap decisive test:** score every universe window with the shipped risk table (already the
`univ` arm in `s8/generate.py`), then build a shortlist by *facility-location* on the top-2000
(maximise coverage of the pool's own geometric diversity subject to score) and compare the top-75
average. One pass, ORACLE-scored, 126 targets, minutes per target with the cached universes.

### C7. `tuning_targets` reads `results/benchmark_manifest.json` at definition time -- **PINNED, audited clean**
`s7/debias.py:103` -> `core/data.py:644 benchmark()` -> `MANIFEST`. The 126 targets are defined as
cluster-disjoint from the sealed 60, so the *definition* touched the manifest once, in S7-S8. The
run-time path never does: audit check 3 (`s29/s29_M_harness_audit.md`) greps every module on the
`chain_rows` path and runs a poison test with the manifest renamed. No action; recorded so that a
new architecture does not re-derive the target list and silently move it.

### C8. FAIL18 is an ORACLE label with an arbitrary 1.5 A band -- **UNTESTED threshold; ORACLE by construction**
`s12/instrument.py:49-50` (the hard-coded list) and `266-286 selfcheck`, which re-derives it:
`band = rr <= rr.min() + BAND` with `BAND = 1.5` (line 36), and a target is FAIL18 iff **no**
in-band member is in the shipped top-75. `rr` is ORACLE, so **FAIL18 is an ORACLE partition of the
instrument** and any arm that conditions on it (gate, weight, route) is an ORACLE arm and must be
labelled so in the same sentence (contract rule 4). Provenance: the "18 zero-recall targets (S10-1)"
(`s12/instrument.py:19`; `docs/FINDINGS.md` S10-1, "the honest target is the 18 zero-recall
targets, where the distogram has no skill at all").
**Untested:** 1.5 A. At 1.0 or 2.0 A the set would be a different size; nobody has published the
sensitivity. This matters because FAIL18-vs-108 splits carry S19's largest effect
(`sequence-conditioning-hurts-the-failures`: blind 5.425 beats shipped 6.019 on FAIL18).
**Cheap decisive test:** print |FAIL_k| for BAND in {1.0, 1.25, 1.5, 2.0, 3.0} and the Jaccard
overlap with FAIL18 (one pass over the 126 universes, seconds). If the set is stable the label is
robust; if it swings, every FAIL18 result needs the band beside it.

---

## B. The score and the energy

### C9. L1 Bayes risk (rather than log-likelihood, L2, or a proper scoring rule) -- **TESTED**
`core/predict.py:420-422`, consumed by `s12/instrument.py:200 shipped_score`.
**Record:** `core/predict.py:13-17` reports the side-by-side measurement that motivated it (the
log-likelihood form scores +0.02 to -0.42 on in-band Spearman, "dominated by the entropy of the
prediction and blind to HOW FAR a wrong distance is wrong", against +0.16 to +0.82 for Bayes risk).
S25 L6/L12 swept the POWER family (`q` from 0.5 to 2): q = 2 (an L2 risk, minimiser = posterior
mean) has `gam_eff` +0.058 at cos +0.253 -- **nearer the truth than the median** -- and the nested-CV
endpoint is +0.0176 (0.96x MDE), i.e. the identity is chosen with full leakage. So: tested, and the
one direction that moves the target toward truth does not move the endpoint. This is the
`better-matrix-worse-ranking` law on the valid instrument.

### C10. The 0.05 A risk grid over 2-40 A and the truncating clip -- **not a convenience worth attacking**
`core/predict.py:419-422`, `s12/instrument.py:203`, `s7/debias.py:238`. The gather uses
`int((d - 2.0)/0.05)` (truncation toward zero, not rounding) clipped to [0, 759). The quantisation
is 0.05 A against a per-pair target quantised to 17 values (C1), so it cannot bind. Recorded for
completeness: the truncation biases every lookup down by a uniform [0, 0.05) A, identically for every
candidate, so it cannot change a ranking.

### C11. `E = zrank(score)` -- the rank ladder -- **TESTED, and it is the mechanism behind five sprints of nulls**
`s27/ham_lib.py:464` / `s24/d_harness.py:274` / `core/pipeline.py:788`; used at
`s27/run_pool.py:68`, `core/pipeline.py:838`.
**Record:** S25 L17 -- because `top = argsort(sc)`, E is the standardised rank ladder 1..128 up to
tie-averaging (worst deviation 1.18% of range on 8 targets, 0 of 8 exactly equal): **the spectrum of
H is target-independent**, there are effectively two trained states in the whole deployment, and
that is why deeper ansatze and larger chi ordered nothing (S21 L16, S21's closure table). The
*monotone* half is justified and tested: raw moment standardisation is not monotone in float64 on
AMBER (breaks argsort on 40/126 targets, ARCHITECTURE section 4), and zrank changes no ordering and
therefore no tail membership. The *choice of ladder* is the convenience.
**Alternatives that preserve the order but not the ladder:** the score's own z (after a monotone
conditioning that is safe on DIS, which has no 1e28 outliers); a gap-preserving transform (rank ->
the score's own quantile spacing); any spectrum with target-dependent gaps.
**Cheap decisive test (already half-run):** S28-L21's spread-spectrum follow-up "altered gradient
scaling without producing an accuracy result" (charter finding 7). The remaining cheap test is
diagnostic, not an endpoint run: measure `corr(E-gap structure, per-target RMSD)` and the realised
tail size m as a function of the spectrum on 12 targets. If m is the only channel (stage 9 of the
map), any spectrum change is a reparameterisation of m and should be priced as such.

### C12. `pad_energy = max(E) + 10 * sd(E)` on the 12 padding states -- **UNTESTED, and provably inert for the SET**
`s22/qcand_lib.py:129, 135`. With 500 candidates in a 512-state register, 12 basis states have no
candidate; they are given a penalty energy so no order-based readout selects them.
**Why it is nearly inert:** the padding states sit at the top of the energy order, so by the
set-equality theorem they are never in the alpha-tail (checked per row: `n_pad_in_tail`,
`s24/d_harness.py:352`). But they are NOT inert for the objective: they carry probability mass that
the entropy term rewards and the CVaR ignores, so `10 sd` sets how much mass the trained state is
willing to waste. 12/512 = 2.3% of the register.
**Cheap decisive test:** re-run `arm_vqe` on 12 targets at `pad_margin` in {1, 10, 100} and report
`m`, `entropy_bits`, `ess` and the point-cloud RMSD. Seconds per cell. If m moves, the padding is a
silent hyperparameter of the readout; if not, the choice is closed.

### C13. Top-128 (production) / 2^n truncation before the VQE -- **UNTESTED as a truncation; PARTLY tested as a shortlist**
`core/pipeline.py:757-759` (`want = max(m, 2^vqe_qubits)`), `core/pipeline.py:836-838`.
Production takes the score's top 128 of 500 and hands the VQE exactly that. In the S25/S27/S28
harness the whole 500 is encoded in 9 qubits instead, so this choice is production-specific.
**What it costs:** the VQE cannot select outside the top-128 (the theorem restricts it further, to a
prefix), so the top-128's own ORACLE ceiling bounds every quantum arm. The record has the adjacent
numbers (top-75 ceiling 2.3062, pool ceiling 1.7108, `s12/instrument.py:19`) but **not the top-128
ceiling**.
**Cheap decisive test (seconds, ORACLE diagnostic):** `min(oracle_rr[order[:128]])` per target over
the 126 universes, mean it, beside the top-75 and pool values. That one number tells every S29 lane
what the 2^7 register is worth as a shortlist before anyone designs a Hamiltonian on it.

### C14. `min_sep = 2` (pairs with j - i >= 2 only) -- **UNTESTED**
`core/pipeline.py:133`, `s12/instrument.py:128`, `core/predict.py:102`. Adjacent CA-CA distances are
excluded from the posterior, the score and the training labels. The stated reason is implicit: the
i,i+1 distance is ~3.8 A in every real chain and carries no information.
**But it does carry one thing:** it is the only pair whose *predicted* value could detect a
contracted candidate, and the emitted average is 22.3% contracted (`ARCHITECTURE.md` section 0). The
score is blind to the one distance that the readout's largest artefact moves.
**Cheap decisive test:** score the 126 pools with min_sep = 1 using the same 17-bin posterior
extended by an analytic i,i+1 term (delta at 3.81 A), and report the top-75 point cloud. Minutes.
Expected null (all candidates are real windows with real 3.8 A bonds, so the term is constant across
candidates) -- which is exactly why it should be stated as measured rather than assumed.

---

## C. The quantum stage

### C15. 7 qubits (production) / 9 qubits (harness) -- **PARTLY TESTED**
`core/pipeline.py:181` (`vqe_qubits = 7`), `s24/d_harness.py:340` (`n = ceil(log2 k)` = 9 at
k = 500), `MAX_QUBITS = 13` (line 125).
**Record:** S25's width sweep (`s25/QUANTUM.md` section 7.1) ran n = 4..13 for gradient variance, and
S21's simulability table ran n = 16 at chi up to 256. Neither is an *accuracy* test of the register
size; the accuracy consequence is C13 (what the register can see).

### C16. RY/CNOT, 3 layers, ring closure -- **TESTED (as depth), and closed**
`core/quantum.py:818-856`; `layers = 3` at `core/pipeline.py:182` and `s27/run_vqe_chain.py:39`.
**Record:** S21's closure table -- "ansatz / bond dimension chi = 1 -> 48: spans 0.37 A, gradients
rise with depth -- CLOSED, L16"; S25 section 7.1's depth sweep (L = 1..12 at n = 7) shows variance
saturating by L = 4 with the deployed L = 3 already close to the plateau; S25 section 7.4 reports
deleting the CNOTs entirely at **-0.013 A [-0.095, +0.077]**. The mechanism for all of it is C11
(there is no target-specific structure for expressivity to capture).

### C17. alpha = 0.18, T = 0.5 (harness) -- **alpha TESTED native-free; T UNTESTED**
`s27/run_vqe_chain.py:39`, `s25/phys_lib.py:71-72`.
**alpha:** pinned by a 12-target native-free probe (`s25/PREREG_PHYS.md` section 1.3): median
realised tail size 75 at alpha = 0.18, against 64 at 0.15 and 84 at 0.20 -- chosen so the VQE and the
classical top-75 control consume the same number of candidates. That is a *matched-control* reason,
which is a scientific reason, and it is stated.
**T = 0.5:** no criterion is recorded anywhere in `PREREG_PHYS.md` or the S25 ledger. It is the one
number in the harness cell with no stated derivation.
**Cheap decisive test:** the entropy/RMSD curve already exists (S25 L15: RMSD tracks readout entropy
at rho -0.7423 and the circuit arms sit on the no-circuit fit at +0.009 A). Re-running T in
{0.25, 0.5, 1.0} on 12 targets and reporting `entropy_bits`, `m` and the point cloud costs minutes
and tells S29 whether T is anything other than a dial on m.

### C18. `VQE_LFO`, with alpha = 1.0 on three of five folds -- **TESTED, and the claim was retracted**
`core/pipeline.py:113`. **Record:** S25 L3/L5 -- alpha = 1.0 on folds 0, 3, 4 means 78 of 126 targets
(61.9%) carry **no tail constraint at all** (at alpha = 1 the CVaR is the mean); the "+0.113 A CVaR
contribution" was withdrawn (measured at T = 0.1, which does not ship; 0.51x MDE as a paired
contrast, fold CI spanning zero, median exactly 0.0000). Forcing alpha < 1 on every fold is an
ORACLE counterfactual worth -0.0311 A at 0.27x MDE (S25 section 6.5). The table is genuinely
leave-fold-out, so reading it is not tuning.

### C19. 80 Adam iterations (harness) / 50 (production), lr 0.15, restarts = 1, theta ~ N(0, 0.6^2), seeds {0, 1} -- **PARTLY TESTED**
`s27/run_vqe_chain.py:39`, `core/pipeline.py:183-184`, `core/quantum.py:1023-1046`.
**Tested:** that the optimiser trains -- S25 section 6.3, against the mandatory control (best-of-200
from the *untrained* circuit, never an initialisation mean): the trained state wins at all three
temperatures and closes 78-89% of the free-energy gap. And that reaching the exact optimum does not
help: the trained state is 0.902 nats / 45% TV from its own Gibbs optimum and the endpoint
difference is 0.24x MDE (S25 L15).
**Untested:** `restarts = 1`, and lr = 0.15. Both are cheap to sweep and, given S25 L15, both are
expected to be inert *through this readout*; the value of saying so is that a new readout (charter
section 12) would change that expectation.
**Note on seeds:** two seeds are run in `s27/run_vqe_chain.py --seeds 0,1` and the S28-L21 anchors
differ by m = 74 vs 71 and by 0.002 A -- the seed moves the rung, not the answer.

### C20. The CVaR tail read as a UNIFORM average over the survivors -- **TESTED, and the alternatives are worse**
`s24/d_harness.py:288 readout_uniform`; the probabilities are discarded at
`s22/qcand_lib.py:352 tail_candidates`.
**Record:** S23 L8 (the probability-weighted readout fails); S28-L21 at n = 126 -- the p-weighted
readout R2 is **+0.32 A** and the p-top-75 readout R3 **+0.26 A** WORSE than production, at every J,
graph and seed, with the fold CI above zero. S23 L5: even a leakage-oracle grid over rank-power and
distance-to-medoid weights selects the uniform point. The mechanism is S23 L5/L9: uniform is the
minimum-variance combination of exchangeable estimators with i.i.d. errors, and every sharpening
attacks the variance reduction that makes averaging work.
**So this is a tested choice, not a convenience** -- but note precisely what it means: the quantum
state's amplitudes are destroyed *because* using them is measurably harmful through this readout.
Any S29 design that wants the amplitudes to matter must change the readout first, and must beat
uniform at matched m.

### C21. m = 75 -- **TESTED (leave-fold-out), and the per-target optimum is unreachable**
`s27/run_pool.py:38`, `core/pipeline.py:127`, from `s8/consensus2.py:92-106`: the leave-fold-out
sweep over 1,120 arms picks the synthesis family in every fold and the cell `sc|75` in three of five
(sc|100 / sc|150 in the other two, within 0.008 A).
**Record:** S22 L4/L5 -- a per-target optimal m is real and transfers across pool halves (-0.239 A)
but **five independent router constructions fail to predict it native-free**, two significantly
harmful held out (S22 L7, S23 L7). S17 L12's mechanism: the readout consumes the set MEAN, so the
rung trades purity against variance reduction.

---

## D. The readout and the projection

### C22. Uniform weights in the coordinate average -- **TESTED** (same evidence as C20: S23 L5/L8, S28-L21).
Recorded separately because it is a different operator (the average over the retained set) from the
tail's weighting. The load-bearing detail: removing the four most geometrically deviant members
costs **+0.142 A** while removing four at random costs nothing (S23 L5) -- the outliers carry the
cancelling error.

### C23. Coordinate averaging at all (rather than selection, or a mode) -- **TESTED, twice, with a correction**
`s12/instrument.py:119`. S17 L17 measured the ORACLE-greedy subset average and found it never beats
the best member on a 5-target smoke; **S17 L18 corrected it at n = 126: the averaging ceiling DOES
beat the best-member ceiling.** Cite L18, not L17. S23 L9 prices the operator: averaging is worth
-0.6554 A (17.7%) on 126/126 targets, and it removes only the 32% idiosyncratic share of the error.

### C24. The MEDOID frame -- **UNTESTED, and it is the one readout choice with no measurement**
`s12/instrument.py:119-125` (`b = medoid(P)`, then superpose all members on member b and mean);
production `core/pipeline.py:934-936` (identical). Every member is aligned to ONE arbitrary member
before averaging.
**Alternatives:** a generalised Procrustes mean frame (iterate to the consensus frame rather than
anchoring on a member); the score-argmin's frame; the frame that minimises the summed squared
deviation (the Frechet mean on SE(3)).
**Why it might matter:** the average of structures superposed on a *member* inherits that member's
rigid placement, and the S23 L9 decomposition is *defined* in the medoid frame ("in the medoid frame
the average is actually taken in"), so the common-mode/idiosyncratic split itself is frame-dependent.
S23 L4 measured iterated Procrustes and found it "does not help, and the contraction is intrinsic" --
that is the closest thing to a test, and it is about the contraction, not the frame choice.
**Cheap decisive test:** replace `medoid` by 3 iterations of generalised Procrustes (initialise at
the medoid) and recompute the 126-target point cloud. `superpose_batch` is already batched; the
whole pass is under a minute. Falsifier: |effect| < 0.7x MDE closes it permanently.

### C25. `ramah` at lam = 0.3 -- **TESTED, and the choice is declared post hoc**
`core/project.py:759`, `s12/instrument.py:139` (`lam=0.3`), `s9/final.py:89-90`.
**Record:** `s8/project.py:118-128` -- the pre-declared rule picked a different arm, and `ramah@0.3`
was carried forward **on a criterion the original rule did not contain** ("This choice is POST HOC
on a criterion the original rule did not contain, and is labelled as such"), at a cost of
+0.004 A [-0.005, +0.013]. The lambda ladder (0, 0.001, ..., 1.0, `s8/project.py:209`) was chosen
geometrically after a single-target pilot showed every penalty saturating by lam ~ 0.3. Dev-24
non-inferiority held an order of magnitude inside its 0.05 A margin. S22 L12 later found the torsion
prior to be the *worst* term in an extended ablation of the objective -- a different objective, but
worth knowing before anyone strengthens it.

### C26. Four multi-start conformations, and multi-start at every rung -- **TESTED, and load-bearing**
`core/project.py:814 STARTS` (extended, alpha, beta, PPII), `933 lam_path(multi=True)`.
**Record:** S9-2 and `core/project.py:912-922`'s docstring -- the projection is degenerate (two
torsion branches at near-equal objective, one plausible), a warm-started optimiser cannot cross
between them, and on 1A1P, 1CS9, 1I6Y the multi-start finds a strictly lower objective with half the
positive-phi rate. The module also measures that the reference **disagrees with itself** under a
rigid motion by up to 1.6 A on some targets, which is why the shipped gradient mode is the bit-exact
one and not the 13x faster analytic one. The count 4 is inherited by value from
`s8/consensus2.FIT_STARTS` and the self-check asserts it is still those four; **whether 4 is enough
is untested** -- a cheap test is to add the medoid's own torsions as a fifth start (the
`extra` argument already exists, `core/project.py:912`) on 126 targets, one pass.

### C27. `maxiter = 300` and `FD_EPS = 1e-5` -- **TESTED (maxiter), inherited (eps)**
`core/project.py:890, 185`; `python -m core.project iters` exists precisely to ask "is maxiter = 300
binding?" (module docstring). `FD_EPS` is the reference's own one-sided step, kept to reproduce
`s8.project` bit-for-bit -- a reproducibility choice, correctly labelled.

### C28. Constant virtual bond 3.804 A (ideal geometry) -- **DECLARED LIMITATION**
`core/project.py:149-160` gives the bonds and angles; the CA-CA distance follows. Consequence: no
cis-peptide (~2.9 A) can be represented. `ARCHITECTURE.md` section 2.6 and `STATE_BRIEF` section 3
declare it; S17 L28 explains the cis-peptide defect and exonerates the force field.

### C29. AMBER k = 10 kcal/mol/A^2, steps = 0, convergence gate 1000 kcal/mol -- **TESTED**
`core/pipeline.py:129-130`, `core/amber.py:1254-1282`.
**Record:** S8-12 priced the stage (k = 10-100 removes ~1e4 kcal/mol of builder strain for +0.011 to
+0.026 A); S16 showed the apparent -0.023 A relaxation gain is against a projection that is itself
0.155 A worse than doing nothing, and a matched-magnitude random displacement is as accurate; S23
L11 measured 17 of 17 restrained/unrestrained settings at or worse than no repair. The gate
threshold 1000 kcal/mol is the same number the codebase already uses for bond+angle strain, declared
before it was applied (`s16/energy_gate.py`). It is a validity stage and is reported as costing
accuracy.

---

## E. The statistics and the harness

### C30. `MDE_K = 2.8016`, `NBOOT = 4000`, the 5-fold cluster bootstrap -- **DERIVED, pinned**
`s24/stats_lib.py:56` (= z_.975 + z_.80 to 5 s.f., pinned literally so every lane's MDE is
bit-identical), `57`, `124-127`. Derived, not convenient. The known limitation is stated in the code:
with only 5 clusters the fold CI is unstable, which is why `_verdict` (145) requires |effect| > MDE
**and** the fold CI to exclude zero, after the audit found a 0.39x-MDE effect labelled "MEASURED".
Memory `mde-is-per-comparison-not-per-instrument` applies.

### C31. The tie key and the lexsort -- **DERIVED from a recorded failure**
`s27/run_pool.py:54-65` (sha256 of `"s27|<pdb>|<tag>"`, not Python's per-process string hash) and
`s27/run_vqe_chain.py:109` (`np.lexsort((key, E))`). Memory: `np.argmin` on a tied signal once read
the ORACLE sort order and invented a 1.386 A winner. Contract rule 12. Not a convenience.

---

## What I would attack first, if the choice were mine

In order of (untested) x (reach), from the map's own accounting:

1. **C1 + C2, the representation of the posterior.** Everything downstream consumes 17 atoms; S25
   L12 proved the *consumption* is closed and explicitly named the open question as "an
   achievability question about a trained predictor". The binning and the loss are the two pieces of
   that predictor nobody has varied.
2. **C6's open half, a coverage-preserving shortlist.** S17 L12 names the class and closes the
   alternative; nobody built one. This is the only untested route to the 1.313 A universe ceiling.
3. **C13, the top-128 ceiling.** A seconds-long ORACLE number that bounds every quantum arm in this
   sprint and is not in the record.
4. **C24, the medoid frame**, because the S23 L9 common-mode decomposition -- the number the whole
   sprint's hypothesis rests on -- is defined in it.
5. **C4, the fold-shared fragment bank**, because it decides what "leave-fold-out" means for every
   CI in the record.
