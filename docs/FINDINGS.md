# FINDINGS

The consolidated research record of this project: **72 findings across Sprints 5-11**, merged
verbatim from the sprint records that used to live in `s5/FINDINGS.md`, `s6/FINDINGS.md`,
`s6/LEAKAGE.md`, `s7/FINDINGS.md`, `s8/FINDINGS.md`, `s9/FINDINGS.md` and `s10/FINDINGS.md`.

**The sprint bodies below are byte-for-byte the originals.** Nothing was paraphrased, rounded,
re-ordered or summarised away: every number, confidence interval, W/L split, p-value and table
stands exactly as it was written by the agent that measured it. Only this front matter is new,
and it exists to do the one job the split files could not do -- **make the corrections impossible
to read past.** Several findings withdraw or revise earlier ones. In separate files a reader
could easily meet a claim and never meet its retraction. The ledger is therefore the first thing
in the document, and every entry names both the claim and the thing that voids it.

**What is not here.** The per-agent working papers (`s8/consensus2_FINDINGS.md`,
`s8/enstarget_FINDINGS.md`, `s8/relax_findings.md`, and the seven `s9/*_FINDINGS.md`) were the
raw material the coordinator records below were written from; each is cited by the finding that
consumed it. They are in git history, as are the Sprint 1-4 reports
(`results/SPRINT2_REPORT.md`, `results/SPRINT3_REPORT.md`, `results/REDESIGN_REPORT.md`,
`docs/*.md`) and the Sprint <=6 working notes from `work/*.md`.
`git log --diff-filter=D --name-only` finds them.

---

## The bottom line

**There is no validated accuracy improvement over the shipped baseline.**

On the 60-target held-out benchmark, which no sprint had touched before its single
pre-registered pass:

| | mean CA-RMSD |
|---|---|
| shipped distogram argmin (baseline) | **2.9507** |
| the Sprint 8/9 synthesis architecture | **2.9610** |
| paired gain | **+0.0103, 95% CI [-0.1596, +0.1803]**, 31W/29L |

against a tuning-instrument gain of **-0.250 [-0.383, -0.117]** and a dev gain of -0.255. The
mean transferred; the *effect* did not. See **S9-10**, and **S10-4** for the adjudication that
the non-replication is explained by concentration (ten targets carried ~60% of the tuning gain)
and **not** by contamination (the leak prices at +0.0004 A).

The architecture is physically better-behaved than the baseline -- exact 3.804 A virtual bonds,
native-like torsion statistics, strain removed -- and no more accurate.

**Where the remaining accuracy actually is** (S10-5; all oracle diagnostics, none deployable):

- Within the shipped top-75: **1.36 A** of headroom -- 0.90 A by selection alone, 1.16 A by
  re-weighting the pipeline's own stack, the last 0.20 A only with oracle per-candidate poses.
- Within the retrieved K=500: **2.35 A** -- 1.49 A by selection, 2.11 A by re-weighting, the
  last 0.24 A with oracle poses.
- **The geometry is not the barrier.** The affine span holds the native at 0.035 A and
  projecting it lands at 0.064 A; the K=500 convex hull emits at **0.853 A**.

The binding constraint is **recognition** -- ranking within the near-native band -- and S9-8
measured that closure as *informational*, not architectural.

---

## Corrections, withdrawals and refutations

Read this table before citing any single finding. The left column is the claim as originally
recorded; the right is what happened to it. Both live in the bodies below.

| Claim as recorded | Status |
|---|---|
| **S10-1's identity audit**: "0/126 tuning pools contain a >=0.6 window (max 0.583)"; "gapped identity exceeds positional by up to 0.417"; the 11/126 whole-library figure | **WITHDRAWN IN FULL -- VOID (S10-4).** Measured against a systematically permuted sequence: the window banks are integer-encoded by `s7.audit.encode` (`ARNDCQEGHILKMFPSTWYV`) and were decoded with `s9.evo.AA` (`ACDEFGHIKLMNPQRSTVWY`). The two orders agree on three letters. True figure: **13/126 targets carry a K=500 pool window at >=0.6 identity, four at exactly 1.0.** S10-1's *conclusion* (leaked windows are never the library best) happens to survive, for an unrelated reason. |
| **S10-2's** "10/126" leaked targets, and its "clean-116" blocks | **REVISED (S10-4)** to **13/126** and clean-**113**. The leak prices at +0.0004 A, so no conclusion moves -- only the counts. |
| **S8-10's transfer law**: `selected = 0.915 * ref_rmsd + 0.346` (r=0.934), read as "accuracy is purchasable at ~0.9 A per A of predictor" and as replacing S8-5's 1.673 incompetence constant | **WITHDRAWN (S8-13).** In the band the pipeline actually occupies the map is the **identity**: `selected = 1.009 * ref + 0.011` (r=0.981), and `1.019 * ref - 0.001` over converged references. The 0.915/0.346 fit spanned 0-8 A anchored at the ref=0 point. Feeding the best available structure through pool selection **costs +0.049 [-0.009, +0.103]**. What survives is cleaner: **the synthesised structure IS the answer; passing it through a selection step adds nothing.** S8-10's oracle datum (best of 100 samples at 2.146 A) stands as measured. |
| **S8-4's class decomposition**, read by the coordinator as an opportunity table -- "the 29 retrieval-limited targets are 28% of the mean and fixable" | **REFUTED (S8-6 and the coordinator correction that follows it).** Class contributions describe where the pool fails, not fractions of the mean available to recover. Injecting the library's best window by cheating is worth **-0.016 A [-0.037, -0.001]**; injecting *every* sub-2 A window is worth -0.070 A and returns the identical structure on 90/126. That 0.960 A share is a **0.004 A** opportunity. The S8-5 law itself over-predicts the value of *injected* pool improvements. |
| **S7-3's shrinkage diagnosis** (the predictor is biased toward the average peptide) -- the coordinator's own, promoted from the audit's noise sweep and used to commission the experiment | **REFUTED by the experiment it commissioned (S7-3, "Correction of record").** The calibration slope, named in advance as the discriminator, came back at **+0.376**; shrinkage requires it well above 1. No synthetic arm reproduces the shipped predictor. The honest description is error correlated across pairs, only 28% aligned with truth. |
| **`esm-adds-nothing-for-short-peptides`**: one-hot 2.216 / pca32 2.133 / pca128 2.212 / raw 2.109, all p >= 0.59 | **REFUTED -- judged by the wrong metric (S7-11).** Every one of those numbers is a distance **MAE**, and S7-3 established MAE is close to useless here (the best-MAE arm of 28 ranked 16th of 28 on selection). Re-judged on selected CA-RMSD, ESM pca32 buys **-0.288 A [-0.484, -0.092], p=0.0046, 74W/45L**. The representation axis is NOT closed. |
| **`real-fragments-beat-the-lattice`** with **`structure-and-sequence-are-decoupled`**, together suggesting BLOSUM retrieval is a worse pool key than random sampling | **REFUTED (S7-12).** On a matched, paired test over an identical window universe, **BLOSUM wins at every K except 2000, where they tie**, and wins on all four measures at once (selected RMSD, pool best, pool mean, oracle ceiling). The original comparison was presumably not matched on window universe or pool size. |
| **S9-4's** loop fit `selected = 0.492 * query + 1.634` (r +0.961), read by the agent as "every 1 A of query improvement buys 0.49 A" | **REFRAMED (coordinator correction).** The derivative is not the relevant comparison; **returning the query itself is, and that has slope 1.0.** The loop is a **contraction toward its own fixed point at ~3.2 A** and never improves on its input. As predictors improve it becomes more useless, not less. |
| **S9-5's** positive-phi "defect", attributed to consensus | **RELOCATED (S9-9).** It belongs to the **projection**. Handed the native CA trace, the arms' own projection returns 0.2513 positive-phi -- higher than every arm. The earlier 0.1973-against-0.0558 comparison set a projected trace against a real backbone. |
| **S10-3's** claim that `hull_floor` upper-bounds every consensus operator "this project has tried or could try" | **QUALIFIED (S10-5).** It does **not** bound an operator that ends in the projection -- and every operator this project ships does. The projection can move a tight hull point *toward* the native. |
| **`s10/forensics`' ceilings** 2.087 / 2.160 / 1.893 | **REVISED (S10-5)** to **1.987 / 2.044 / 1.802**. Cause: a real solver defect -- projected gradient at the 1/L step stops short on an ill-conditioned Gram, and its error grows with set size (+0.092 at m=75, +0.193 at 225, +0.296 at 500). A frame-convention mismatch was excluded by measurement (agreement to 2.7e-13). |
| **`s10/mdgen`'s** 1.604 / 1.124 / 0.828 | **RESCOPED (S10-5).** Those are **32-target** figures; on the full 126-target instrument they are **1.802 / 1.336 / 0.953**. The 32 targets are simply easier. |
| The coordinator's premise that **oracle weighting cannot be worse than the hull** | **REFUTED (S10-5).** The two sets are not nested -- posing each candidate individually at its best destroys the relative geometry that averaging exploits. The hull wins on the mean and on 369/378 cells and **loses on 9**, by up to 0.085 A, verified at 400,000 Frank-Wolfe iterations. |
| **S8-1's** use of native percentile under an oracle reference | **CORRECTED in place (S8-1, "Correction to the brief"):** the statistic is degenerate under an oracle reference. |
| **S7-10's** interim result at n=50 | **SUPERSEDED** by the same finding's FINAL pass at n=70. |
| **S8-14's** tuning-instrument sign | **DID NOT REPLICATE on dev**: +0.070 [-0.118, +0.259]. |
| The **CVaR gradient** as shipped in `qansatz.cvar_gradient` | **DEFECT, found and fixed (S5-6, priced against an exact reference in S9-5).** The baseline is centred on the tail only, so the correction term `b * E_p[grad log p]` is not zero; cosine with the true gradient measured at -0.023. |
| The **pairwise ranker arm** of Sprint 6 | **CORRECTNESS BUG, found before it was reported (S6-3).** |

### Two pairs of numbers that look like disagreements and are not

The bodies below quote both members of each pair, in different sprints, without always
saying which construction produced them. They are **different objects, both correct**.
Neither is an error, and neither supersedes the other.

- **3.2005 vs 3.204 - single-start vs multi-start synthesis.** `3.2005` is
  `s9.synth.fit_w`, one solve started from the medoid's torsions; `s9/synth_probe.json`
  records it under arm `fit`. **The pipeline emits `3.204`, the multi-start value** - it
  is what `s9/synth.py`'s own `EXPECT` states and what benchmark60's pre-registration
  records as `tuning_mean`. All 126 targets differ between the two. Quote **3.204** for
  what the system does; quote 3.2005 only for the single-start probe.
- **2.815 vs 2.954 - raw structure vs emitted structure.** `2.815` is S10-5's bound
  ladder measuring the *raw* uniform mean over the top-75, before projection. `2.954` is
  what that object becomes after the manifold projection, at a price of
  +0.139 [+0.102, +0.176]. **The raw figure is not a pipeline-emitted quantity at all** -
  a raw hull point is ~20% contracted in contour and is not a peptide. S10-5's ladder
  gives both columns for exactly this reason; read the EMITTED column for anything the
  system could return.

Two structural lessons the record itself draws, worth carrying forward:

1. **"This is the fourth time in two sprints that a quantity which looked like an opportunity
   turned out to be a description."** (S8-6.) The others: the inter-generator agreement signal
   (a difficulty proxy), the AMBER in-band advantage (a compactness proxy and a pool-size
   artifact), and the ESM/one-hot tie (an MAE artifact). **The reliable tell is that none of
   them was measured with a control that isolated the causal claim from the descriptive one.**
2. **Concentration was reported as a caveat throughout four sprints. It was in fact the
   result.** (S9-10.) Ten targets carried ~61% of the tuning gain; after dropping the ten
   largest the system is +0.2013 A *worse* than the baseline (+0.3296 after twenty).

### Sprint 26 additions to this ledger (2026-09-14; `s26/RETRACTIONS.md`, `s26/LEDGER.md`)

The S26 record lives in `s26/`; the entries below are the corrections that touch claims
quoted from this file or from the S12 to S25 ledgers, plus S26's own retractions, so that a
reader of this ledger meets them. Nothing below the front matter is edited.

| Claim as recorded | Status |
|---|---|
| **S25 QUANTUM.md / state brief 5.4**: "no barren plateau at any width measured" (log2 Var per qubit -0.649 to -0.311, n = 4..13) | **SCOPED by S26 (R2, ledger L27, L47, L119).** The deployed ansatz's dynamical Lie algebra is the full so(2^n) from depth 2 (dim 8128 at n = 7; `s26/results/q_dla.json`, cross-checked in `s26/results/a_dla_check.json`), so the measured slopes are a statement about shallowness at depth 3, not about a favourable algebra. The numbers stand; every quotation carries "at depth 3". |
| **S7-11**: ESM-pca32 is worth -0.288 A over one-hot on selection (artefact `s7/repr_tune.json`) | **ARTEFACT LOST; RE-MEASURED (ledger L11, L62).** The file is absent from disk and git. S26's ladder measures the ESM-2 650M block at -0.330 A on selection and -0.208 A on the built chain, 5/5 folds (`s26/results/p_ladder_noesm_s0.json`). The direction and size are confirmed; quote the S26 figure with its artefact. |
| **S10-4**: the identity-leak price +0.0004 A (dev) / +0.0030 A (benchmark); the S10 artefacts `s10/idaudit_*.json` | **REPRODUCED / RE-SOURCED (ledger L31, L44, L50).** The dev half reproduces exactly on the same basis (+0.0004 [-0.0001, +0.0010]; +0.0018 [-0.0003, +0.0039] on the production built chain; `s26/results/w_selfcopy_endpoint.json`); the artefacts are in git history at `5fa05cd`, not absent. The benchmark half stays historical; S26 bounds the 2/60 case from the dev proxy at 0.028 A on the mean (0.151 A at the worst target, 2BP4), class MINOR, without opening the benchmark (L44, L58). |
| **S25 L16 / state brief 5.3**: the seven-configuration suite (3.058 ... 3.881), the random-75 null 3.4251 and Legacy +0.330 / AMBER +0.455, quoted beside the built-chain 3.2148 | **BASIS NAMED (ledger L28 item 3, L29).** Those are point-cloud numbers (`s25/results/phys_suite.json :: basis`); the built-chain means of the same rows are 3.2187 to 4.1015 (`results/summary/leaderboard.json`). The ranking is unchanged; every S26 document names the basis on both sides. |
| **S16 L27 / state brief 5.3**: AMBER relaxation's gain dissolves under a matched random displacement | **CONFIRMED AND SHARPENED (ledger L39, L46, L87, L100).** On the production emission the relaxation is worse than a random move of its own size by +0.0111 A (fold CI [+0.0062, +0.0171], Type-M zone at 1.12x MDE) and worse than a move toward a pool member by +0.0385 (2.27x), replicated on new seeds; it is a validity step on 124 of 126 (34 emissions with a sub-2 A heavy-atom overlap become 1) and breaks a virtual bond on 2BP4 and 9KAR. |
| **S22 L7 / S23 L7**: five router constructions fail to predict per-target m* native-free | **EXTENDED (ledger L110, L115).** Twelve more m* routers and six s* routers on six feature blocks no previous router used: none clears its MDE, all but one point the harmful way. |
| **S26 L28 item 2** (lane E): the governed test run "369 / 356 / 13" | **CORRECTED (R1, L31, L32):** 370 / 357 / 13 at the time; with the opt-in tier 370 / 368 / 2 (L113). |
| **S26 L68** (lane Q): "ADAPT under L-BFGS selects no operator at alpha = 1" | **RETRACTED (R5, L75):** operators are appended on 60 / 68 of the 78 alpha = 1 targets and are inert (at most 1.2e-4 nats); the verdict is unchanged. |
| **S26 L53** (lane PH): the relaxation's displacement is "the first native-free quantity above rho 0.4" for the error | **RETRACTED AS WORDED (L121, L123):** the pool's own spread carries the signal (partial rho +0.452); the displacement tracks it at rho 0.76 and adds +0.08 given it. |
| **S26 L99** (lane P): the esm8m selection contrast "-0.30" | **CORRECTED (L101):** +0.225 A (0.88x MDE, not measured). |
| **S26 L106** (lane P): the B3 gain-size feature is the pool's helix content or the length | **CORRECTED (L107):** it is the pool's strand content. |
| **S26 L56** (lane P): the C2 anchor's built chain equals the production cache | **CORRECTED (L57):** selection and point cloud are bit-exact; the built chain lands on the leaderboard-rebuild basis 3.2126 because the multi-start projection is sensitive to a 1e-14 input difference. |
| **S26 L7** (lane I): two point-cloud numbers for 1D6X and 1KWE | **CORRECTED (L9):** written before the query returned; the correct values are 2.0900 and 2.7313. |
| **S26 L38** (lane PH): the cis floor 0.347 A (the native rebuilt from its own torsions) | **SCOPED (L89):** that figure is the own-torsion upper bound; the native projected through the production projection sits 0.083 A from itself, so the constant omega costs about 0.04 A. |
| **S13 geo_FINDINGS**: the Pauli mean weights 2.236 (Legacy) and 3.015 (conditioned AMBER) | **NOT RE-DERIVED BY S26 (L122, L126):** kept off every slide; the paper outline cites them from `s13/results/geo_pauli.json` with that note; the measured/predicted ratio 0.9969 over 104 cells is what S26 reproduced. |

---

## Index

The headings are the verdicts. Bold entries carry a correction.

### Sprint 5 -- the wall is located
1. The pool was never the problem
2. Search is not the constraint at any budget
3. The 3.3 A wall
4. The wall is not "no signal" -- it is no resolution AT THE TOP
5. What this leaves
6. **A defect in the shipped CVaR gradient**
7. In-band skill is bounded by the prior's own accuracy -- confirmed
8. Winner's curse is real, and small
9. Structure and sequence are decoupled at this length -- retrieval by sequence cannot work
10. The torsion prior is null -- the sixth signal to hit the same wall
11. The causal chain, closed
12. **Physics DOES rank real geometry -- Sprint 4's measurement was confounded**
13. Retrieval by predicted structure: better selection, worse pool, no net gain
14. The synthesis is null too -- and the interim numbers were another mirage
15. Sprint 5 verdict

### Sprint 6 -- selection is measured out
1. The benchmark does not have a floor -- the ceiling is ours
2. A learned ranker fixes the global ranking and not the in-band one
3. **A correctness bug in the pairwise arm, found before it was reported**
4. Properly powered: the learned rankers are exactly null
5. Building a structure instead of picking one: null
6. What the failure actually looks like
7. Tail aggregation of pair violations: null
8. VQE/CVaR assembly: a narrow pool, and why its narrowness cannot be exploited
9. Refinement: both hypotheses answered, and the second one explains the first
10. The full ablation, and a usable confidence signal
11. Anchored selection: null
12. Decomposed all-atom physics in a learned combiner: null
13. The distance channel is sufficient; the predictor is not

### Sprint 7 -- the predictor is the constraint
1. The inter-generator agreement lead is a difficulty thermometer -- closed
2. Thirty times more structural data makes the predictor monotonically worse -- closed
3. **The predictor is not shrunk, it is uncorrelated -- post-hoc de-biasing is dead**
4. Chirality is a red herring on the generative path -- closed
5. The objective is misspecified: the native is not its minimum
6. Pool size has an interior optimum, and selection tracks the pool MEAN not its BEST
7. A ceiling worth stating plainly
8. Training against ranking makes it worse, and the predictor is ~= typicality
9. The triangle operator improves the diagnostics significantly and the outcome not at all
10. No available scoring channel ranks the truth first -- AMBER is worse, not better
11. **ESM DOES beat one-hot -- on selection. The recorded ablation was judged by the wrong metric**
12. **BLOSUM retrieval beats random fragments -- the recorded claim is refuted**

### Sprint 8 -- full architectural attack
- S8-1. The representation is not the bottleneck -- the aggregation metric is
- S8-2. Multimodality is real but modest: 2-3 populated basins, not 30
- S8-3. One actionable lead, and it is priced honestly
- S8-4. What the 3.454 A mean is actually made of -- four classes, and a 28% opportunity
- S8-5. Pool construction is capped at 2.406 A at ANY filter skill -- and the law that says why
- S8-6. Retrieval is fixable; fixing it is worth 0.016 A. The 28% "retrieval opportunity" is 0.4%
- **Coordinator correction: class contributions are not recoverable fractions**
- S8-7. Inverse folding is the least misspecified objective measured, and it is still misspecified
- S8-8. Inside the near-native band: the score is worse than a coin flip, and consensus is the only thing that discriminates
- S8-9. The channels' errors ARE decorrelated, and it buys 0.01 A -- plus what consensus actually is
- S8-10. Generation is worth nothing as candidates -- but it exposes a second law with a far lower constant
- S8-11. Synthesis beats selection: 3.204 A, replicated on dev
- S8-12. AMBER is a validity stage, not an accuracy stage
- **S8-13. CORRECTION to S8-10: the "transfer law" is the identity map where we operate**
- S8-14. The training target IS misspecified -- by 0.5 A per pair -- and fixing it is worth 0.044 A

### Sprint 9 -- breaking the 2.0 A barrier
- S9-1. The shared bias is real, is one interpretable mode, and is not an error
- S9-2. A torsion prior buys physical validity for free, and no accuracy
- S9-3. The relational channel is typicality wearing a tournament
- **S9-4. The loop is a contraction toward its own fixed point, not a multiplier**
- S9-5. Refinement fails, and the failure is located exactly: it is all selection
- S9-6. Synthesis is at its ceiling, and the binding constraint is the distance target
- S9-7. Evolutionary information exists, is retrievable, and does not help
- S9-8. The in-band ceiling is informational, not architectural -- the recognition line closes
- **S9-9. Consensus is structurally the wrong operator for a refinement trajectory**
- **S9-10. THE FINAL BENCHMARK PASS: the gain does not replicate**

### Sprint 10 -- adjudication and bounds
- S10-1. The attrition ledger: the loss is set PURITY, and the failure is a minority class *(its identity audit is VOID -- see S10-4)*
- S10-2. The projection path is closed: the required matrix is 2.7x better than anything achievable
- S10-3. Physics does not supply decorrelated error -- the force field agrees with retrieval, target-specifically
- **S10-4. The identity-audit conflict adjudicated: an alphabet permutation, and the leak is worth +0.0004 A**
- **S10-5. The reconciled bound ladder: the geometry is not the barrier**

### Sprint 11 -- performance engineering and consolidation

- **S11-1. The four mandated components run genuinely, and three of them buy nothing measurable** *(VQE weights vs uniform: -0.031 A at 49W/49L; Legacy refinement +0.030; AMBER is the only significant effect and it is a +0.021 A COST for validity; the whole system is +0.017 [-0.019, +0.053] against the like-for-like pipeline)*
- **S11-2. The projection is not a function of its input at angstrom resolution** *(under a rigid motion of the target cloud, which leaves the CA-RMSD exactly invariant, the REFERENCE disagrees with itself by up to 1.62 A -- the same magnitude as the fast arms' divergence, so bit-exactness is the only defensible standard and is what shipped)*

---

## Reproducing any of this

Every finding below names the module that produced it and the JSON it wrote. Those modules were
consolidated into `core/` (see `README.md`); the sprint packages they were written in are in git
history. `python -m core.bench --arm baseline` runs the reference path these numbers were
measured on.

The instrument discipline that makes any of these numbers mean something, and that must be
preserved:

- **The 60-target benchmark (`results/benchmark_manifest.json`) is held out.** It was touched
  exactly once, by S9-10, on pre-registered constants. All optimisation runs on the 126-target
  tuning instrument; `dev_set(24)` is the cluster-disjoint dev split.
- **Folds and clusters are pinned.** Correcting the identity clustering once silently moved 13
  benchmark targets and invalidated every fold model.
- **Leakage is audited per finding**, not per project, and the audit is reported with the result
  (identity thresholds, fold disjointness, NaN-poisoning controls at 0.000e+00).
- **The known filter defect stands** (S9-10, S10-4): `peptide_db.identity` normalises by the
  LONGER sequence, so a long fragment can carry a target's exact k-mer and still align below
  0.6. It is priced at +0.0030 A on the benchmark and +0.0004 A on tuning, and it applies to
  the member-level filter everywhere in this repository.

---

<!-- ==================================================================== -->
<!-- Sprint 5: verbatim from s5/FINDINGS.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s5/FINDINGS.md`. Not edited.*

# Sprint 5 findings

Every number here is on `peptide_db.dev_set(24)` - 24 targets, cluster-disjoint from both
the 60-target and 52-target benchmarks - unless stated otherwise. Sprint 4's held-out
baseline was **pool best 2.355 Å, selected 3.324 Å**.

## 1. The pool was never the problem

| measurement | mean CA-RMSD |
|---|---|
| nearest training window (oracle retrieval over all windows) | **1.175** |
| best single whole window, rebuilt from its own torsions | 1.168 |
| retrieval pool best, top-2000 by BLOSUM | 1.356 |
| retrieval pool best, 100 **random** real fragments | 1.969 |
| **fragment assembly capacity** (L5/s2/K200, greedy+CD+16 restarts) | **0.887** |
| Sprint 4 lattice pool best | 2.355 |

One hundred *arbitrary* real protein fragments beat the entire learned-distogram + lattice
+ VQE pipeline. Not leakage: the 6,003 fragments come from 1,001 PDB files with zero
overlap with the peptide database, the benchmark, or the dev set.

Sequence scoring adds almost nothing on top of random - BLOSUM top-100 is 1.787 against
random-100 at 1.969 - and the reason is that BLOSUM sums over 5-7 columns are
integer-valued, so the top-K out of ~30,000 windows is mostly an arbitrary draw from a
large tied set (two seeds share 7% of their sources).

## 2. Search is not the constraint at any budget

Assembly capacity as a function of search effort, same library throughout:

| search | mean |
|---|---|
| oracle greedy, best fragment per segment independently | 3.781 |
| greedy + coordinate descent | 1.282 |
| random best of 3,000 | 1.513 |
| + 16 restarts | 0.887 |
| **deployable, no search (top-BLOSUM per segment)** | **4.371** |

All but the last use RMSD-to-native as the objective. Any competent search finds ~1 Å
structures when it knows what it is looking for.

Oracle *greedy* (3.781) being far worse than coordinate descent (1.282) on the same
objective is direct evidence that inter-segment coupling carries most of the information -
i.e. that per-position independent choice, the old lattice's factorisation, is wrong.

## 3. The 3.3 Å wall

Five unrelated selectors, all on real-geometry pools whose best member is at 1.5-1.8 Å:

| selector | selected |
|---|---|
| learned distogram | 3.31 (K=100) / 3.47 (K=500) |
| pool density / clustering / weighted density | 3.28 - 3.91 |
| coarse physics `energy_terms` (Sprint 4) | 3.33 |
| all-atom Amber ff14SB/GBn2 (Sprint 4) | 3.33 |
| prior + Amber rank consensus (Sprint 4) | 3.27 |
| *random* | *4.20 - 4.31* |

Every selector recovers about a quarter of the gap between random and pool-best, then
stops - learned or physical, local or global, score-based or density-based. And the
distogram returns the same ~3.4 Å whether the pool's best member is at 1.79 Å or 1.36 Å:
**it is insensitive to the quality of what it is ranking.**

The reading: in a pool where every candidate is already a real peptide conformation,
plausibility is not a discriminating signal, because plausibility is what got them into
the pool. All five selectors measure some form of "is this a typical peptide". None
measures which conformation *this sequence* adopts.

Corollary observed twice: selection gets **worse** as the pool improves (distogram 3.31 →
3.47 from K=100 → K=500; density 3.59 → 3.85). A weak scorer given more candidates has
more chances to be wrong.

## 4. The wall is not "no signal" - it is no resolution AT THE TOP

Preliminary, 6 fold-0 dev targets (the folds whose torsion model finished first), 500-candidate
retrieval pools. Two Spearman correlations against true CA-RMSD: over the whole pool, and
*in-band* - restricted to candidates within 1.5 Å of the pool's own best, which is the only
set a selector's decision actually ranges over.

| signal | global rho | **in-band rho** |
|---|---|---|
| learned distogram | **+0.653** | **+0.091** |
| torsion prior | +0.379 | +0.172 |

The distogram is *excellent* - it separates good structures from garbage at +0.65 - and
essentially blind within the good set. That is this project's own recorded pattern,
"garbage rejection reads as skill", now measured on real fragment pools rather than the
synthetic decoy bank where it was first seen. My earlier reading of the 3.3 Å wall as "no
signal" was wrong: there is a lot of signal, and none of it is where the decision is made.

The torsion prior is *worse* globally and roughly twice as good in-band. Both are weak.

### The shortlist ceiling - why every reranking experiment here failed

The best structure that survives filtering to the distogram's own top-M, i.e. the hard
ceiling on any reranker applied after that filter:

| top-M | ceiling | lost vs full pool |
|---|---|---|
| 1 | 3.110 | +1.853 |
| 5 | 2.663 | +1.406 |
| 25 | **2.437** | +1.180 |
| 50 | 2.059 | +0.801 |
| 100 | 1.958 | +0.701 |
| 250 | 1.482 | +0.225 |
| 500 (all) | 1.257 | 0 |

**The shipped pipeline reranks `--rescore-top 24`.** So Sprint 4's headline negative - that
all-atom Amber reranking is worth +0.004 Å - was measured on a shortlist whose ceiling is
2.44 Å against a pool best of 1.26 Å. The filter discarded 1.18 Å of the available accuracy
*before the force field was ever called*. Amber was not given a chance to fail properly.

That does not overturn the Sprint 4 result as stated (Amber genuinely did not improve the
answer in that configuration) but it does overturn the interpretation. "Physics cannot rank
peptide structures" is not what was measured. What was measured is "physics cannot rank the
24 structures a globally-good, in-band-blind filter hands it".

## 5. What this leaves

The pool is solved. The search is solved. The remaining error is entirely selection, and
the missing quantity is sequence specificity. The one channel this project has never used
is backbone torsions - the published state of the art for this length class (APPTEST,
1.57 Å best-of-100 on 6-12mers from a third of our data) predicts distances *and* φ/ψ,
where we predict only distances.


## 6. A defect in the shipped CVaR gradient

`qansatz.cvar_gradient` - the update that has driven every VQE run in this project -
centres its REINFORCE baseline over the alpha-tail rows and applies it to the tail rows
only. That is not valid: `E[grad log p] = 0` holds over the whole sampled distribution, not
over a data-dependent subset of it, so subtracting a tail-conditioned mean leaves a bias
that does not shrink with more shots.

Verified by enumerating all 256 bitstrings of an 8-qubit ansatz and central-differencing
the exact CVaR objective. A correctly-baselined gradient agrees to **5.8e-10** relative.
The shipped one, averaged over 400 independent batches at alpha=0.2:

| shots | full baseline | tail baseline (shipped) |
|---|---|---|
| 128 | rel. bias 0.076 | 1.013 |
| 2048 | rel. bias 0.041 | **1.128** |

The full baseline's residual is the known O(1/N) plug-in quantile bias and shrinks with
shots; the tail baseline's is the size of the gradient itself and does not. At 2048 shots
the shipped gradient has **cosine -0.023 with the true CVaR gradient** - essentially
orthogonal to the direction it is supposed to be descending.

Two honest qualifications, both measured:

* **Correct is not automatically better.** On a marginal-constrained `chain` ansatz the
  biased tail baseline *outsearched* the unbiased one (argmin hit rate 0.400 vs 0.325) -
  its bias is essentially the cross-entropy method's move, which is a real optimiser. The
  ordering reverses on the default `chain_ry` ansatz (0.650 vs 0.575), which is why the
  corrected version is the default, but correctness alone did not settle it.
* The practical effect on Sprint 4's results is **not** established. The search was not the
  binding constraint there - the objective was - so a wrong gradient may have cost little.

`s5/quantum.py` defaults to the corrected baseline and keeps the shipped one, asserted
bit-identical to `qansatz`'s, so the finding is about `qansatz` and not about a divergent
rewrite.

### VQE over segment choices: encoding and honest limits

One-hot over segments is **unsampleable**, not merely expensive: on the real 1AL1 library
(S=5, K=16, 80 qubits) **0 of 20,000 sampled bitstrings were feasible**, against a
predicted per-segment rate of 2.4e-4. Binary/Gray over a geometry-sorted candidate order
uses 20 qubits for the same problem (40 at K=200) and needs no feasibility repair. Measured
argmin hit rate: one-hot 0.100, binary/Gray 0.650, random 0.050.

On planted chain instances at matched objective evaluations, VQE beats random 0.625 vs
0.300 (S=4, K=8) and 0.650 vs 0.025 (S=5, K=8). But:

* Most of that gain is the trailing RY layer, not entanglement - a product ansatz with no
  entangler scores 0.600 against `chain_ry`'s 0.650. Entanglement is worth <= 0.05 of hit
  rate here, inside 40-seed noise, and should not be claimed to be doing work.
* On an unstructured needle, random search beats VQE 0.133 vs 0.033, as it must. The wins
  above come from real h+J structure in the instance.
* Nothing here is a folding measurement. The one real-library run uses a stand-in
  compaction+clash objective on a single target.


## 7. In-band skill is bounded by the prior's own accuracy - confirmed

All 24 dev targets, 500-candidate real-fragment pools. Spearman of the distogram score
against true CA-RMSD, restricted to candidates within `b` A of the pool's own best, and its
rank correlation with the prior's own MAE on that target:

| band | mean rho | candidates | **corr(rho, MAE)** |
|---|---|---|---|
| 1.0 A | -0.070 | 72 | **-0.464** |
| 1.5 A | +0.025 | 134 | **-0.604** |
| 2.0 A | +0.126 | 185 | **-0.762** |
| 3.0 A | +0.314 | 296 | **-0.791** |
| global | +0.545 | 500 | - |

Mean prior MAE is **2.219 A**. The prediction was that a signal cannot resolve differences
finer than its own noise, so in-band skill should track accuracy and vanish once the band is
narrower than the error. That is exactly what the data shows, and the relationship
strengthens monotonically as the band widens toward the noise scale.

So the wall has an information-theoretic explanation, not an architectural one:
**you cannot rank structures 1.5 A apart using a prior whose distance error is 2.2 A.**

This does not contradict Sprint 4's retraction of the MAE law. That measured MAE against
*selected RMSD* across variants whose MAE all sat in a narrow 1.8-2.1 A band, where nothing
changes. This is a threshold effect: MAE has to fall below the band width before in-band
ranking becomes possible at all. It gives the first concrete target the sequence model has
had - **MAE below ~1.5 A**, from 2.22 A now.

## 8. Winner's curse is real, and small

The argmin of a noisy score is the candidate whose score is most wrong. Reducing the
distogram's top-M by something that does not need fine resolution, all 24 dev targets:

| M | ceiling | mean | random | medoid | trimmed |
|---|---|---|---|---|---|
| 5 | 3.142 | 3.425 | 3.416 | 3.288 | 3.372 |
| 25 | 2.657 | 3.404 | 3.396 | **3.253** | **3.222** |
| 50 | 2.389 | 3.447 | 3.432 | 3.355 | 3.367 |
| 200 | 1.868 | 3.678 | 3.665 | 3.621 | 3.587 |

argmin (the current rule) is **3.470**; the trimmed medoid of the top-25 is **3.222**, a free
0.25 A from the same score and the same pool.

Two things in that table matter more than the 0.25 A. **The mean of the shortlist is flat at
~3.4 for every M from 5 to 200** - filtering harder does not improve the average candidate at
all, it only moves the ceiling. And the ceiling at M=200 is 1.868 A while every reduction of
that same shortlist returns ~3.6. The information is present and no reduction reaches it.

The distogram's whole contribution is a single step: it lifts the pool from mean 4.3 A to
mean 3.4 A and then adds nothing, which is +0.65 global and +0.09 in-band restated in
angstroms.


## 9. Structure and sequence are decoupled at this length - retrieval by sequence cannot work

Audit of the 1.175 Å neighbour ceiling, all 24 dev targets. For each target: the highest
sequence identity to any training peptide it can retrieve from, and the identity of the
window that actually achieves the best RMSD.

| quantity | mean | max |
|---|---|---|
| max identity to any training peptide | 0.422 | 0.571 |
| **identity of the best-RMSD window** | **0.122** | 0.583 |

Windows above 0.5 identity: **1 of 24**. The project's fold-separation threshold is 0.6, so
nothing here is leakage - but the more interesting fact is that the structurally nearest
fragment is, on average, an **essentially unrelated sequence**. 1AL1's best match at 0.48 Å
shares 0% identity with it; 5XER's at 0.31 Å shares 0%; 2N08's at 0.61 Å shares 0%.

This single number explains four separate observations at once:

* **Sequence retrieval barely beats random** (BLOSUM top-100 1.787 vs random-100 1.969). Of
  course it does - the right fragment does not look like the target sequence, so ranking by
  sequence similarity ranks by the wrong thing.
* **The BLOSUM score is degenerate** (7% source overlap between seeds). Not just a tie-breaking
  artefact; there is little to break ties *with*.
* **ESM-2 appears to contribute nothing** to distance prediction here (preliminary: a no-ESM
  model matched the shipped distogram's MAE).
* **"The sequence signal is the ceiling"** - three ESM-2 routes failed, and this says why.

The correct reading is *not* that peptide structure is unpredictable from sequence. It is
that the mapping is many-to-one: many unrelated sequences adopt the same local conformation,
so structural neighbourhood does not imply sequence neighbourhood, and **sequence similarity
is the wrong retrieval key**. A retrieval architecture has to retrieve by *predicted
structure*, not by sequence - or generate rather than retrieve.

That reframes the assembly result. Its library is selected by BLOSUM, which this says is
close to an arbitrary draw; the 0.887 Å capacity is therefore a property of *real fragment
geometry being dense at this length*, not of the retrieval working. Both the capacity and
the failure of the deployable `seq_top1` baseline (4.37 Å) follow from that.


## 10. The torsion prior is null - the sixth signal to hit the same wall

All 24 dev targets, leave-fold-out mixture-of-von-Mises torsion prior scoring each
candidate's own parent torsions (no reconstruction, nothing inferred):

| K | pool best | distogram | **torsion** | dist+tors | random | rho dist | rho tors |
|---|---|---|---|---|---|---|---|
| 100 | 1.787 | 3.309 | **3.320** | 3.514 | 4.176 | +0.513 | +0.337 |
| 500 | 1.534 | 3.470 | **3.544** | 3.420 | 4.312 | +0.545 | +0.349 |
| 2000 | 1.356 | 3.427 | **3.798** | 3.463 | 4.272 | +0.572 | +0.376 |

The torsion channel lands on the wall exactly where everything else does: 3.320 A against
the distogram's 3.309 at K=100. Combining the two is *worse* than either (3.514). Its global
rank correlation is lower than the distogram's at every K.

This was the most strongly motivated idea in the sprint. It is the one architectural
difference between this project and the published state of the art for this length class
(APPTEST predicts distances *and* phi/psi and reaches 1.57 A best-of-100 on 6-12mers), it is
a genuinely independent channel rather than a reparameterisation of distances, and it was
verified correct at the level of the von Mises normalisation, the mixture integrating to 1
on the torus, and the sampler's circular moments. It simply does not rank real peptide
conformations.

One caveat kept honestly: on the 6 fold-0 targets measured early, the torsion prior's
*in-band* rho was +0.172 against the distogram's +0.091 - roughly double. That comparison
was not recomputed in-band across all 24, and the end-to-end number is what decides.
A signal with slightly better in-band rho and clearly worse global rho selects no better,
which is consistent with everything else here.

### Where the wall count stands

| selector | selected |
|---|---|
| learned distogram | 3.31 |
| **torsion prior (new)** | **3.32** |
| pool density / clustering | 3.28 - 3.91 |
| coarse physics `energy_terms` | 3.33 |
| all-atom Amber ff14SB/GBn2 | 3.33 |
| prior + Amber rank consensus | 3.27 |
| *trimmed medoid of distogram top-25* | *3.22* |
| *random* | *4.20 - 4.31* |

Pool best over the same candidates is 1.36 - 1.79 A.


## 11. The causal chain, closed

Per-target Spearman against selected CA-RMSD, 24 held-out targets:

| predictor | rho vs selected |
|---|---|
| **prior MAE** | **+0.856** |
| distogram global rho | -0.789 |
| windows under 2 A in the library (log) | -0.622 |
| library coverage (nearest training window) | +0.568 |
| pool best @500 | +0.560 |
| chain length | +0.307 |

And decomposing the selection gap (selected minus pool best, mean 1.936 A):

    corr(gap, prior MAE)          = +0.770
    corr(gap, library coverage)   = +0.122

**The gap is driven by the prior's accuracy and essentially not at all by how well the
fragment library covers the target.** That completes the chain every measurement in this
sprint has been pointing at:

    prior MAE  ->  in-band ranking skill  ->  selection gap  ->  selected RMSD
      2.219 A        +0.03 at 1.5 A            1.94 A            3.47 A

Each link is measured: MAE 2.219 A (section 7), corr(in-band rho, MAE) -0.60 to -0.79
(section 7), corr(gap, MAE) +0.770 (here), rho(MAE, selected) +0.856 (here).

### The best case is still not good enough

Splitting the dev set at the median library coverage:

| half | pool best | selected |
|---|---|---|
| well-covered (nn <= 1.07 A) | 0.82 | **2.59** |
| poorly-covered (nn > 1.07 A) | 2.25 | 4.35 |

Even where the library contains a 0.82 A structure, selection returns 2.59 A. There is no
subset of this benchmark, however favourable, where the current selection machinery reaches
2.0 A.

### And the prior cannot be improved from sequence - all four arms, complete

| arm | distance MAE | paired vs `onehot` | 95% CI | p | better on |
|---|---|---|---|---|---|
| `onehot` - no ESM at all | 2.216 | - | - | - | - |
| `pca32` - the shipped input | 2.133 | -0.083 | [-0.399, +0.234] | 0.614 | 15/24 |
| `pca128` | 2.212 | -0.004 | [-0.434, +0.426] | 0.986 | 14/24 |
| `raw` - 1280-d, learned projection | **2.109** | -0.107 | [-0.492, +0.278] | 0.591 | 13/24 |
| shipped distogram | 2.219 | - | - | - | - |

**A 1280-dimensional protein language model, given a learned projection and the same
training budget, does not beat a 20-dimensional one-hot encoding at predicting the distances
of a short isolated peptide.** All four arms sit within 0.11 A, no confidence interval
excludes zero, and none is better on more than 15 of 24 targets.

One correction to my own reading. Seeing `pca32 ≈ pca128` I concluded that more ESM
dimensions cannot help. That inference was invalid: PCA ranks directions by variance, not by
relevance to structure, so adding principal components says nothing about whether a
*learned* projection would find informative low-variance directions. The `raw` arm is the
test that settles it, and its fold-0 result (1.773 A, six targets) briefly looked like a
positive before folds 1-4 came in at 2.92 / 2.29 / 2.30 / 1.00. Single folds of 3-7 targets
are not evidence here, which is the same lesson section 8's power analysis records.

Combined with section 9 (the best-matching window has 12% sequence identity), the reading is
that ESM-2's representation is largely orthogonal to which conformation a short *isolated*
peptide adopts. The requirement from section 7 was MAE below ~1.5 A. The best of four
arms reaches 2.109 A.


## 12. Physics DOES rank real geometry - Sprint 4's measurement was confounded

24 dev targets, top-100 BLOSUM retrieval pools, scored on each window's own continuous
torsions via `geo.build_backbone_batch` - real Ramachandran, nothing snapped to a state
library. AMBER through `amber_refine.refine_coords` (not `refine`, which would project the
fragment back onto the discrete lattice and reintroduce the confound under test), k=10,
steps=0, nonbonded+solvation. 6.31 s/structure, 1172/1200 scored, 28 strain-rejected, 0
exceptions.

| scorer | rho whole pool | **rho in-band (<3 A)** | selected |
|---|---|---|---|
| `legacy_all` | +0.365 ± 0.128 (21/24) | **+0.200 ± 0.126** (p=0.007, n=16) | 4.200 |
| `amber` | +0.351 ± 0.160 (19/24) | **+0.320 ± 0.211** (p=0.012, n=13) | 3.837 |
| `rg` (control) | +0.393 ± 0.143 (21/24) | **-0.001 ± 0.207** (p=0.99, n=16) | 4.467 |

The whole-pool column is a trap: radius of gyration, a one-line control, leads it. The
finding is the second column. Inside the near-native band, where compactness carries
*exactly zero* information, both Legacy and Amber stay significantly positive. **The
lattice-build confound was real** - Sprint 4 measured physics against candidates whose local
geometry the lattice had distorted, and on real geometry physics does rank.

Note the comparison that matters: **Amber's in-band rho (+0.320) is more than three times
the distogram's (+0.091)**. It is the best in-band signal found anywhere in this sprint.

**And it still does not convert.** Selected RMSD: Amber 3.837, Legacy 4.200, against the
shipped sequence-only distogram's 3.309 on the same pools. Gain over a random pick is
+0.475 ± 0.487 for Amber (15/24, p=0.069) and +0.112 for Legacy (p=0.77). Converged
all-atom ff14SB/GBn2 at 6.3 s/structure picks *worse* structures than the prior already in
the pipeline, because its global rho (+0.351) is far below the distogram's (+0.545) and
selection ranges over the whole pool.

Legacy and Amber are wrong the same way: between-scorer rho **+0.763 ± 0.041**, and both
≈+0.56 against plain `rg`. The force field is largely re-deriving the cheap model's
ordering, a third of which is compactness.

Per-term: `compactness` +0.410, `rg` +0.393, `solvation` +0.357, `hbond_local` +0.334,
`coop_helix` +0.322. `contact` **anti-ranks** (-0.266, positive on only 5/24). One lead:
`coop_sheet` is far stronger in-band (+0.483 ± 0.102) than overall (+0.007) - the only term
that behaves like near-native discrimination rather than garbage rejection.

Caveats kept: in-band rows are n=13-16 of 24, since targets with <4 in-band candidates
cannot support a correlation, which biases those rows toward easy pools. And per-candidate
scores were summarised rather than persisted, so the one question left unanswered is the
partial correlation - does Amber add anything once `rg` is regressed out?

### The synthesis this implies

The two signals are complementary in exactly the way the architecture needs:

| signal | global rho | in-band rho |
|---|---|---|
| distogram | **+0.545** | +0.091 |
| Amber | +0.351 | **+0.320** |

The distogram is the better *filter* and Amber is the better *ranker*, and Sprint 4 ran them
in that order with a top-24 shortlist whose ceiling was 2.44 A. Widening the filter is the
one combination this project has never tested at a ceiling that permits a good answer.


## 13. Retrieval by predicted structure: better selection, worse pool, no net gain

Section 9 said sequence similarity is the wrong retrieval key. Testing the alternative --
score all ~10-40k leakage-safe windows by the predicted distance matrix and take the top-K:

| key | K | pool best | trimmed medoid |
|---|---|---|---|
| `seq` (BLOSUM) | 100 | **1.787** | 3.624 |
| `seq` | 500 | **1.534** | 3.758 |
| `struct` (distogram) | 100 | 2.519 | 3.280 |
| `struct` | 500 | 2.203 | **3.261** |
| `hybrid` | 500 | 1.582 | 3.465 |

Oracle over all windows: 1.175.

The two keys trade off exactly as the global/in-band split predicts. The structure key
produces a pool the selector picks *better* from (3.261 vs 3.758 at K=500) and a pool
containing *worse* structures (2.203 vs 1.534). Net, it does not beat the best reduction
already measured (3.222).

The more interesting number is that **the distogram's top-500 out of ~30,000 windows has
pool-best 2.203 A, while BLOSUM's top-500 has 1.534 A** -- and BLOSUM is close to an
arbitrary draw. Ranking every window by predicted structure concentrates the pool on
whatever structural type the prior predicts, and systematically misses the native basin more
often than an essentially random selection of real fragments does. A prior with 2.2 A error
used as a retrieval key over a large library is worse than no key at all.


## 14. The synthesis is null too - and the interim numbers were another mirage

Distogram top-50 filter, then rank within it. 24 held-out targets, paired against the
current rule (distogram argmin over the whole pool):

| arm | selected | d vs argmin | 95% CI | p | wins |
|---|---|---|---|---|---|
| distogram argmin (current) | 3.398 | - | - | - | - |
| trimmed medoid of top-25 | 3.285 | -0.113 | [-0.448, +0.222] | 0.516 | 13/24 |
| Legacy argmin in top-50 | 3.288 | -0.110 | [-0.477, +0.257] | 0.563 | 10/24 |
| AMBER argmin in top-50 | 3.329 | -0.069 | [-0.362, +0.225] | 0.650 | 11/24 |
| **dist+AMBER rank-sum in top-50** | **3.269** | -0.129 | [-0.376, +0.119] | 0.319 | 12/24 |

pool best 1.523 · **ceiling of the top-50 filter 2.357** · gap ceiling-to-best-selector
**+0.912**.

No arm is significant. Every confidence interval contains zero.

**And the interim read was wrong again.** At 20 of 24 targets the same table showed 2.890
for the consensus arm and I reported it; the last four targets moved it to 3.269. That is
the fourth time this sprint a partial result has evaporated - after the swlin_both/gauss_g0
pattern in Sprint 4, the raw-arm fold-0 MAE, and the torsion prior's fold-0 in-band rho. The
rule from section 8's power analysis (~11 cells for 0.3 A, ~25 for 0.2 A) exists precisely
for this and I keep reading numbers before they are resolvable.

**What the run does establish** is the sprint's central claim in its most favourable
setting. After filtering 500 real fragments down to 50, of which the best is at 2.357 A,
the best available selector returns 3.269 A. **Selection loses 0.91 A inside a
50-candidate shortlist of real structures.** Physics, learned prior, consensus and density
all fail there together.

## 15. Sprint 5 verdict

**<2.0 A was not achieved, and no intervention produced a statistically significant
end-to-end improvement.** Best held-out result 3.269 A (dist+Amber rank-sum, p = 0.319)
against Sprint 4's 3.324 A. The honest summary is that the deployable number did not move.

What did move is the diagnosis, which is now mechanical rather than empirical:

1. **The pool is not the constraint** - 100 random real fragments (1.969 A) beat the entire
   Sprint 4 pipeline (2.355 A); assembly capacity is 0.887 A.
2. **The search is not the constraint** - coordinate descent reaches 1.282 A under an oracle
   objective; per-segment greedy gives 3.781 A, so the coupling carries the information.
3. **Selection fails for an information-theoretic reason** - in-band rho +0.09, and
   corr(in-band rho, prior MAE) = -0.60 to -0.79 with MAE 2.219 A against a 1.5 A band. You
   cannot rank structures 1.5 A apart with a 2.2 A signal.
4. **That accuracy is not reachable from sequence** - a 1280-d PLM with a learned projection
   ties 20-d one-hot (2.109 vs 2.216, p = 0.591), because structure and sequence are
   decoupled at this length (best-matching window: 12% identity).
5. **Two Sprint 4 conclusions were wrong** - physics was confounded by lattice geometry
   (Amber in-band rho +0.320 on real geometry, 3x the distogram's), and `--rescore-top 24`
   capped every reranking experiment at a 2.44 A ceiling.
6. **A real defect exists in the shipped CVaR gradient** - tail-centred baseline, cosine
   -0.023 with the true gradient at 2048 shots.

<!-- ==================================================================== -->
<!-- Sprint 6: verbatim from s6/FINDINGS.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s6/FINDINGS.md`. Not edited.*

# Sprint 6 findings

## 1. The benchmark does not have a floor - the ceiling is ours

21 of 24 dev references are multi-model NMR ensembles. Measured intra-ensemble CA-RMSD:
**mean 0.956 A, median 0.633 A** (min 0.07, max 4.12) - tighter than the 1.2-1.9 A the
literature reports for short-peptide ensembles.

| correlation | rho |
|---|---|
| ensemble width vs **selection gap** | **+0.119** |
| ensemble width vs selected RMSD | +0.277 |
| ensemble width vs pool best | +0.651 |

Reference uncertainty does not explain the selection gap. Wide ensembles make the *pool*
worse (flexible peptides are harder to match with rigid fragments) but not *selection*. For
16 of 21 targets the pool holds a structure closer to the reference than the reference's own
models are to each other. Two targets are genuinely floppy (8DHZ 4.12, 9U7M 2.54) and are a
minority.

## 2. A learned ranker fixes the global ranking and not the in-band one

117,000 labelled candidates over 234 training targets, 64 features, leave-fold-out at the
target level. First pass (the `rank` arm of this run was invalid - see §3):

| arm | selected | d vs dist | p | wins | rho | **in-band rho** |
|---|---|---|---|---|---|---|
| `dist` (baseline) | 3.435 | - | - | - | 0.551 | **0.013** |
| `reg` (learned) | **3.332** | -0.103 | 0.561 | 13/24 | 0.532 | **0.018** |
| `top` (learned) | 3.593 | +0.158 | 0.392 | 9/24 | 0.476 | 0.097 |
| `legacy` | 4.297 | +0.862 | 0.075 | 8/24 | 0.389 | **0.161** |
| `tors` | 3.544 | +0.109 | 0.631 | 9/24 | 0.349 | 0.098 |
| `density` | 3.648 | +0.213 | 0.509 | 11/24 | 0.473 | -0.015 |
| `medoid` | 3.758 | +0.323 | 0.345 | 11/24 | - | - |

pool best 1.534 · pool mean 4.393 · recall<2 A 10.6%

`reg` is the best arm at 3.332 A and is not significant (p = 0.56). The informative number
is not the mean, it is that **the learned model's in-band rho is 0.018 - the same as the
distogram it was meant to beat.** Fitting a ranker on the deployment distribution, with hard
pairs, did not create in-band skill.

### Corrected pairwise arm

With the Borda-count scoring (§3) and the 0.40 A pair gap, the pairwise arm becomes the best
learned model:

| arm | selected | d vs dist | p | wins | rho | in-band rho |
|---|---|---|---|---|---|---|
| `rank` (corrected) | **3.291** | -0.145 | 0.539 | 12/24 | **0.568** | 0.059 |
| `reg` | 3.332 | -0.103 | 0.561 | 13/24 | 0.532 | 0.018 |
| `dist` | 3.435 | - | - | - | 0.551 | 0.013 |

It has the best global rho of anything measured (0.568, above the distogram's 0.551) and it
does move in-band skill from 0.013 to 0.059 - a real direction, four times the distogram,
and still an eighth of what Amber manages. Neither learned arm is significant at n=24.

**The reason is the feature set, not the model.** A learner cannot extract information its
inputs do not contain, and the inputs here are the same quantities already measured to be
in-band blind. Note the ordering in the last column: the arms with the most in-band skill
are the physical ones (`legacy` 0.161, `top` 0.097 which leans on the Legacy terms), and the
one signal measured at **+0.320 in-band** - converged all-atom Amber - is **absent from the
feature set entirely**, because it costs 6.31 s per structure and 500 candidates x 234
training targets is 200+ hours.

That is the gap to close, and it is a compute-allocation problem rather than a modelling one.

## 3. A correctness bug in the pairwise arm, found before it was reported

The `rank` arm trained a gradient-boosted classifier on feature differences x_a - x_b and
then scored candidates with `decision_function` on raw features. Valid for a LINEAR model,
where f(x_a - x_b) = w.x_a - w.x_b makes w.x a per-candidate utility; invalid for trees,
where f is not linear. The arm was scoring in a space the model was never fitted in. Its
numbers (3.880 A) are discarded, not reported, and archived in
`s6/evaluate_v1_invalid_rank.json`.

Inference is now a Borda count over 48 sampled opponents. The pair-sampling gap also rose
from 0.15 to 0.40 A, with the floor set by §1's measurement rather than by guesswork:
ordering two candidates 0.1 A apart asks the model to reproduce a distinction the experiment
itself does not resolve.


## 4. Properly powered: the learned rankers are exactly null

Evaluated on **234 held-out training targets**, leave-fold-out (a target in fold *f* scored
by the model fitted without fold *f*). Standard error on a mean is **0.113 A**, against
~0.30 A on the 24-target dev set.

| arm | selected | d vs dist | 95% CI | p | wins | rho | in-band |
|---|---|---|---|---|---|---|---|
| `dist` | 3.354 | - | - | - | - | 0.541 | 0.125 |
| `rank` | **3.347** | **-0.007** | [-0.164, +0.151] | 0.93 | 109/234 | 0.539 | 0.124 |
| `reg` | 3.365 | +0.011 | [-0.122, +0.144] | 0.87 | 104/234 | 0.546 | 0.064 |
| `top` | 3.440 | +0.085 | [-0.080, +0.251] | 0.31 | 102/234 | 0.493 | 0.145 |
| `tors` | 3.781 | +0.427 | [+0.230, +0.623] | 0.0000 | 95/234 | 0.308 | 0.111 |
| `density` | 3.617 | +0.263 | [+0.096, +0.429] | 0.0022 | 99/234 | 0.416 | -0.054 |
| `legacy` | 5.122 | +1.768 | [+1.496, +2.039] | 0.0000 | 43/234 | 0.257 | -0.039 |

**The learned rankers are null to within ±0.15 A.** d = -0.007 and +0.011, confidence
intervals tight around zero, and the study resolves 0.22 A effects. Fitting 117,000 labelled
candidates on the deployment distribution, with hard-pair sampling and a corrected Borda
reduction, does not beat the distance prior it was built to beat.

**This retires the dev-set numbers I reported earlier tonight.** `reg` at 3.332 and `rank` at
3.291 against the distogram's 3.435 looked like 0.10-0.15 A improvements; at n=234 they are
+0.011 and -0.007. Those were noise at n=24, exactly as the standard error predicted, and I
reported them anyway. The high-power instrument existed the whole time -- the pools were
already built -- and building it should have preceded reporting any dev comparison.

It also runs the leakage control from `s6/LEAKAGE.md`. Dev flattered the learned arms
relative to held-out training (-0.10 vs +0.011), which is the direction the indirect pathway
would predict but is also exactly what small-n noise produces. Either way the held-out
training numbers are the trustworthy ones, and they are null.

**The in-band column is unchanged by learning**: 0.125 for the distogram, 0.124 for the
pairwise ranker, 0.064 for regression. Learning redistributed nothing. Two arms are
significantly WORSE than the distogram (`tors` +0.427, `density` +0.263, `legacy` +1.768),
which is a useful check that the instrument has power -- it detects real differences when
they exist.


## 5. Building a structure instead of picking one: null

Sprint 4 tried averaging and it failed for a specific reason -- it averaged CARTESIAN
coordinates, broke chain geometry (max |CA-CA - 3.8| = 1.27 A), and snapping back landed
0.31 A worse than the medoid. Torsion-space averaging cannot fail that way: bond lengths and
angles are fixed by `build_backbone`, so any (phi, psi) vector is a valid backbone by
construction. Done with circular means throughout (atan2 of summed unit vectors -- the
arithmetic mean of 170 and -170 degrees is 0, the opposite side of the circle from both).

24 dev targets, averaging the distogram's top-M:

| M | mean | score-weighted | cluster-restricted | concentration R |
|---|---|---|---|---|
| 3 | 3.728 | **3.403** | 3.728 | 0.837 |
| 5 | 3.865 | 3.431 | 3.828 | 0.863 |
| 10 | 4.088 | 3.486 | 4.106 | 0.833 |
| 25 | 3.881 | 3.619 | 3.896 | 0.825 |
| 50 | 3.795 | 3.780 | 3.748 | 0.832 |

distogram argmin 3.435 · medoid 3.758 · pool best 1.534

Nothing beats the argmin. The best arm (score-weighted over the top 3, 3.403) is the argmin
with extra steps, and unweighted averaging gets monotonically *worse* as M grows -- averaging
in more candidates averages in more wrong ones. Concentration R sits near 0.83, so the
inputs do agree reasonably well; the average is a real conformation, it is just not a better
one. The geometry objection that killed this in Sprint 4 is fixed and the idea still does
not work, which is the cleaner way to close it.

## Where Sprint 6 stands

Eight distinct selectors have now been measured on the same pools: the learned distance
prior, coarse physics, all-atom Amber, pool density, torsion prior, trimmed medoid,
three learned rankers, and torsion-space consensus. **Every one lands between 3.29 and 3.44 A**
except the ones that are significantly worse. The only high-powered comparison (n=234,
SE 0.113) says the learned rankers are null to within +-0.15 A.

Against that: the pool's best member is at 1.53 A, the library ceiling is 0.86 A, the
benchmark's own references are tight (0.96 A mean ensemble spread), and reference
uncertainty does not explain the gap (rho +0.119). The candidates are there and nothing
finds them.


## 6. What the failure actually looks like

Every measurement so far is a scalar. This is what the selector is doing, per target.

**The best candidate is ranked at the 39th percentile.** By the distogram, over 24 dev
targets: mean percentile 39.3%, median 36.4%. In the scorer's top 1% on **1 of 24** targets;
top 10% on 5 of 24; top 50% on 14 of 24. That is close to a random position.

This is not "nearly right with a resolution problem at the top". Global rho of +0.55 comes
from pushing genuinely bad candidates down, and the single best candidate is not preferred
at all.

**The pick is structurally unrelated to the best candidate.** Mean CA-RMSD between the
selected structure and the pool's best member: **3.52 A**. The selector is not choosing a
near-miss of the right answer, it is choosing a different conformation.

**Radius of gyration is not the problem.** Native 6.20, pool best 6.30, selected 6.21 -
the pick has the right overall size. It is the wrong shape at the right scale.

### Secondary structure is not the discriminating axis either

The pick matches the native's dominant SS class on 15 of 24 targets against the pool best's
20 of 24, which looks like a gross, fixable failure - SS prediction from sequence is ~80%
Q3, far more accurate than distance prediction, and this project has never used it.

So I bounded it with an oracle: filter candidates to those whose (helix, sheet, other)
vector is within tau of the **native's**, then take the distogram argmin inside.

| tau | candidates kept | ceiling | selected | gain |
|---|---|---|---|---|
| 0.15 | 27.5 | 2.422 | 3.224 | **+0.211** |
| 0.30 | 79.1 | 1.821 | 3.229 | +0.206 |
| 0.50 | 167.1 | 1.703 | 3.350 | +0.085 |
| 0.80 | 273.3 | 1.552 | 3.436 | -0.001 |

against 3.435 unfiltered. **A perfect secondary-structure predictor is worth 0.21 A** - and
a tight filter also costs ceiling (2.42 at tau=0.15, against a pool best of 1.53), so the
usable gain is smaller still. Building an SS predictor is not worth it, and the failure is
not gross SS mismatch even though the class-match counts suggested it might be.

The discrimination that is missing is fine-grained, and every coarse handle -- size, shape
class, secondary structure, compactness, density -- has now been measured and bounded.


## 7. Tail aggregation of pair violations: null

The literature review's specific suggestion: near-native errors are LOCALISED, so a mean
over ~90 pairs dilutes the few that distinguish a wrong basin; score by the worst pairs
instead. Same pools, same predicted distances, only the aggregation changes.

| aggregation | selected | rho | in-band rho |
|---|---|---|---|
| worst 5% of pairs | 3.597 | 0.544 | -0.029 |
| worst 10% | 3.428 | 0.554 | -0.005 |
| worst 25% | **3.366** | 0.559 | 0.018 |
| worst 50% | 3.376 | 0.548 | 0.022 |
| all pairs (mean) | 3.370 | 0.538 | 0.032 |
| single worst pair | 3.669 | - | - |

Nothing. The best tail setting (3.366) and the plain mean (3.370) are the same number, and
in-band rho stays at 0.02-0.03 throughout. The premise -- that near-native structures are
distinguished by having no badly-violated pair -- is not true of this pool, presumably
because with a 2.2 A prior MAE every candidate has badly-violated pairs and which ones they
are carries no signal.

## Selection approaches measured and rejected, Sprints 4-6

Learned distance prior · coarse physics (Legacy) · all-atom Amber argmin · pool density ·
clustering · torsion prior · trimmed medoid · learned regression ranker · learned pairwise
ranker · learned saturating ranker · torsion-space consensus (plain, weighted, clustered) ·
secondary-structure filtering (oracle-bounded at +0.21 A) · tail aggregation of pair
violations · structure-based retrieval · shell reweighting · realizability projection ·
rank consensus.

Seventeen approaches. Every one lands between 3.27 and 3.44 A on held-out targets, or is
significantly worse. The only high-powered test (n=234, SE 0.113) puts the learned rankers
at -0.007 and +0.011 against the distance prior.


## 8. VQE/CVaR assembly: a narrow pool, and why its narrowness cannot be exploited

The end-to-end ablation put CVaR-VQE where the measurements said it might help -- expanding
the candidate pool by searching the segment-assembly space -- rather than on selection. On
10 dev targets:

| pool | best | selected | **selection gap** | candidates |
|---|---|---|---|---|
| retrieval | **1.459** | 3.286 | **1.827** | 500 |
| VQE assembly | 2.636 | **3.034** | **0.397** | 96 |
| union | 1.322 | 3.066 | 1.744 | 596 |

The VQE pool's selection gap is **4.6x smaller**. It is a worse pool that is far easier to
select from, because the search concentrates candidates near the objective's minimum and the
argmin of a homogeneous pool is representative of it. Arithmetic that looked like the
sprint's best lead: a VQE pool as good as retrieval's (1.459) with the VQE pool's gap (0.397)
would select at **1.86 A**.

**It does not work, and the reason is worth more than the lead was.** Raising the search
budget and the segment library, measured directly:

| target | K | evals | objective min | pool best | selected |
|---|---|---|---|---|---|
| 1GJE | 16 | 11,520 | 1.459 | **0.924** | **0.994** |
| 1GJE | 16 | 34,560 | 1.459 | 0.924 | 0.994 |
| 1GJE | 64 | 34,560 | 1.483 | 1.267 | 1.555 |
| 1KVG | 16 | 11,520 | 1.564 | 2.276 | 5.250 |
| 1KVG | 16 | 34,560 | 1.551 | 3.759 | 4.734 |
| 1KVG | 64 | 34,560 | **1.495** | **3.107** | 4.021 |

On 1KVG, tripling the evaluations and quadrupling the library **lowers the objective from
1.564 to 1.495 and makes the structure worse, 2.276 to 3.107 A.** Optimising the objective
better produces worse answers, monotonically. This is Sprint 4's "better search under a bad
objective is worse" measured directly as a dose-response rather than inferred.

So the VQE pool's narrowness cannot be exploited: making it better requires a better
objective, and with the objective we have, searching harder actively moves away from the
native.

**And where the objective happens to be right, the architecture works.** On 1GJE the
objective's minimum sits at a genuinely good structure, and the same machinery returns
**0.994 A selected** against retrieval's 2.97 A on that target. The pipeline is not broken.
It is correctly converting a good objective into a good answer, and the objective is right
on a minority of targets.


## 9. Refinement: both hypotheses answered, and the second one explains the first

8 dev targets, K=100, baseline pool best 1.523 A / mean 4.059 A. Idealisation (fragment ->
rebuild from its own torsions) costs -0.007 A, i.e. free, and is charged to the baseline.

**H1 -- does refinement improve candidates? No.** Amber's best case is **-0.015 A** at
restraint k=10 (6/8 targets, magnitude at noise), and the restraint sweep is a flat null:
mean per-candidate change is +0.009 to +0.023 A at every k, with 57-78% of candidates
getting worse. There is no restraint strength at which all-atom ff14SB/GBn2 becomes an
accuracy operator. Torsion-space refinement makes pool best strictly worse (+0.073 to
+0.449 A).

**H2 -- does refinement narrow the pool? Yes, unanimously, and it backfires.** Free torsion
refinement improves the pool mean on **8/8** targets and shrinks its spread on **8/8** -
a large, consistent effect. And selection gets worse in **every one of the 8 regimes**:

| | none | tors_t10 | tors_free |
|---|---|---|---|
| selected RMSD | **2.549** | 2.918 | 2.787 |
| rho(score, RMSD) | 0.651 | 0.642 | 0.259 |
| **in-band rho** | **0.351** | 0.283 | **0.095** |
| score sd / baseline | 1.000 | 0.799 | **0.139** |
| random - selected | +1.510 | +0.658 | +0.331 |

The mechanism is exact: the pool narrows because every candidate is pulled toward the same
minimum of the same prior, so the selector's own score collapses with it - spread falls to
14% of baseline, in-band rho from 0.351 to 0.095, and the selector's advantage over a random
pick from +1.51 A to +0.33 A. **Narrowing that comes from optimising the selector's own
objective cannot help the selector.**

The diagnostic split names the damage. Mean change by starting quality:

| regime | started <2 A (n=112) | 2-4 A (n=291) | >4 A (n=397) |
|---|---|---|---|
| amber_k10 | -0.002 | +0.000 | +0.022 |
| tors_free | **+0.360** | -0.207 | **-1.847** |

Torsion refinement drags bad candidates in and pushes good ones out - regression toward a
generic prior-preferred conformation. Rank preservation confirms it: rho(pre, post) falls
from 0.995 to 0.374, and the pool's best candidate survives as best in **1 of 8** targets.

### This is the same finding as section 8, reached from the other direction

Both say: **a pool concentrated near the objective's minimum is easier to select from in a
way that is vacuous, because the minimum is in the wrong place.** The VQE search
concentrates candidates there deliberately and gets a 4.6x smaller selection gap that cannot
be cashed in; refinement concentrates them there physically and destroys the selector's
discriminative power outright. Two different mechanisms, one cause.

Discipline note worth recording: the refinement study verified no-native-access rather than
asserting it - replacing the native with noise and checking every regime returns
bit-identical coordinates (max |delta| = 0.0 across all 9 regimes).


## 10. The full ablation, and a usable confidence signal

24 dev targets, each arm differing from the previous by one component:

| arm | pool best | selected | d vs A | 95% CI | p | wins |
|---|---|---|---|---|---|---|
| A retrieval only | 1.523 | 3.398 | - | - | - | - |
| B + VQE/CVaR assembly | **1.465** | **3.310** | -0.088 | [-0.409, +0.233] | 0.597 | 8/24 |
| C + Legacy in the objective | 1.517 | 3.351 | -0.047 | [-0.137, +0.042] | 0.310 | 2/24 |
| *VQE candidates alone* | *2.862* | *3.312* | - | - | - | *gap 0.450* |

**VQE/CVaR is null on the mean and large on a minority.** It improves 5 of 24 targets by
more than 0.5 A -- 1GJE 2.97 -> 0.99, 6GS5 2.50 -> 1.29, 5ZGD 3.42 -> 2.82, 2MI1 4.54 ->
3.77, 8DHZ 8.28 -> 6.67 -- and does nothing or slightly worse elsewhere. That is what
"materially participates" honestly means here: it changes the answer substantially on 21%
of targets, in the right direction on those, and is neutral on average.

Legacy inside the quantum objective is worth -0.047 A on 2 of 24 targets and makes the VQE
pool slightly *worse* (2.862 -> 2.914). It does not earn its place in the search objective.

### Inter-generator agreement predicts accuracy

The two pools are produced by genuinely different processes -- sequence-similarity retrieval
of real fragments, and a CVaR-VQE search over segment assemblies. When their independently
selected structures agree, both are more likely to be right. Over 12 targets, using only the
CA-RMSD between the two picks (no native):

    corr(agreement, error of the better pick)      = +0.594
    corr(agreement, |retrieval error - VQE error|) = +0.699

Routing on it helps a little (3.069 at a 3 A threshold, against retrieval's 3.277 and an
oracle per-target choice of 2.909) but the useful product is the **confidence estimate**:
this is the first quantity in the project that predicts, at inference and without a native,
whether a given prediction is good. Given that AlphaFold2's pLDDT is reported to have *no*
correlation with which of its peptide samples is closest, a +0.594 confidence signal for
peptide structure prediction is worth keeping even though it does not move the mean.


## 11. Anchored selection: null

If the VQE pool says *where* (small selection gap, concentrated near the objective's
minimum) and retrieval says *what* (better structures, but its argmin sits at the 39th
percentile), the obvious combination is to take the retrieval candidate nearest the VQE
pick. 14 dev targets:

| rule | selected | d vs retrieval | p | wins |
|---|---|---|---|---|
| retrieval argmin | 3.384 | - | - | - |
| VQE argmin | **3.103** | -0.282 | 0.210 | 9/14 |
| retrieval candidate nearest the VQE pick | 3.382 | -0.002 | 0.989 | 8/14 |
| distogram argmin among the 25 nearest to it | 3.336 | -0.048 | 0.768 | 5/14 |

Both anchored rules are null. Transplanting the VQE's *region* onto retrieval's *library*
gains nothing -- the retrieval candidates near the VQE pick are no better than the VQE pick
itself, which is another way of saying the objective's minimum is in the wrong place and
every route to it arrives at the same wrong place.

Worth noting the subset effect: VQE-alone is -0.282 here on 14 targets and -0.086 on all 24.
Same direction, different magnitude, n too small to separate them -- the recurring lesson
of this project, recorded again rather than glossed.


## 12. Decomposed all-atom physics in a learned combiner: null

2,480 converged `refine_coords` calls (124 targets x 20 shortlisted), 98.0% scored, 1.98%
strain-rejected, 0 exceptions, 4.97 s/structure, 3.52 h.

**The headroom frames everything.** Full pool best 1.534; **ceiling of the stage-1 top-20
shortlist 2.874**; stage-1 selected 3.332. The shortlist discards 1.34 A of the pool's
1.80 A of headroom *before Amber sees anything*, and an ORACLE re-ranker of that shortlist
gains only **-0.458 A [-0.636, -0.281]**. That is the ceiling on the whole experiment.

| arm | selected | d vs stage-1 | 95% CI | p | in-band rho |
|---|---|---|---|---|---|
| **stage-1 (learned, cheap features)** | **3.332** | - | - | - | +0.013 |
| Amber total | 3.498 | +0.166 | [-0.045, +0.377] | 0.117 | +0.033 |
| Amber interaction-only | 3.508 | +0.176 | [-0.045, +0.397] | 0.113 | -0.074 |
| Amber nonbonded | 3.648 | +0.316 | [+0.105, +0.526] | **0.005** | +0.065 |
| Amber angle | 3.490 | +0.157 | [+0.011, +0.304] | **0.036** | -0.083 |
| stage-2 learned (regression) | 3.409 | +0.077 | [-0.026, +0.179] | 0.135 | +0.014 |
| stage-2 learned (centred) | 3.463 | +0.131 | [-0.020, +0.282] | 0.087 | +0.046 |

Nothing beat stage-1; two physics argmins are significantly worse. In-band rho stays at ~0
for every arm. **Giving a learned combiner decomposed all-atom physics created no in-band
skill.**

### The feature importances contradict a standing conclusion

Out-of-fold within-shortlist rho is +0.056 / +0.081 -- the model has almost no skill to
distribute -- but it does use the physics (Amber is ~33% of total permutation importance).
Per Amber signal, both objectives independently rank:

| | reg | cen |
|---|---|---|
| **torsion** | **+0.0190** | **+0.0168** |
| solvation | +0.0138 | -0.0002 |
| bond | +0.0112 | +0.0047 |
| restraint_rmsd | +0.0085 | **+0.0130** |
| **nonbonded** | **-0.0002** | +0.0049 |

**Torsion first, nonbonded last.** That inverts this project's standing position, recorded
as `amber-converged-interaction-only-ranks`, that the right expression of Amber is
interaction-only (nonbonded + solvation). Given all five components, a learned combiner
ignores the one the project adopted. And under the primary objective `restraint_rmsd` --
how far the backbone moved to relieve strain, not an energy at all -- outranks four of the
five energy terms.

### And the +0.320 did not replicate

Amber's in-band rho of +0.320, the number that motivated this entire experiment, was
measured on a BLOSUM top-50 pool. Under the same <3.0 A band on a learned-ranker shortlist
it is **-0.079 (interaction-only) / +0.117 (total)**. The signal was conditioning-dependent,
and the conditioning that produced it is not the one a two-stage system uses.

**Read this as "the shortlist is the binding constraint", not as "physics has nothing".**
With n=24, SE ~0.3 A, and an oracle ceiling of 0.458 A, the experiment could only ever have
detected a near-oracle re-ranker.


## 13. The distance channel is sufficient; the predictor is not

DeepH3 reports that 6D inter-residue orientations discriminate near-native 8-18 residue
loops better than distances -- the one literature recommendation left untested. Bounded with
an oracle on both sides, so the CHANNEL is isolated from any predictor's accuracy: rank
candidates by agreement with the native's orientations, versus agreement with the native's
distances.

| channel | global rho | in-band rho | selected |
|---|---|---|---|
| oracle distances | **0.830** | **0.641** | **1.836** |
| oracle 6D orientations | 0.542 | 0.511 | 1.904 |

Orientations are worse on every measure. The DeepH3 result does not transfer to isolated
peptides -- plausibly because its discriminating power came from antibody framework context,
which a free 9-15mer does not have.

**And the more important half of this table**: oracle distances rank in-band at **+0.641**
and select at **1.836 A** against a pool best of 1.534. Perfect distance knowledge very
nearly solves the problem with the machinery already built. The learned distance prior
delivers **+0.013** of that +0.641.

That is the sprint's conclusion in one line. The channel is right, the pool is right, the
search is right, the selection machinery is right. **The sequence to distance predictor
delivers 2% of the in-band skill its own channel is capable of**, and nothing downstream
can manufacture the other 98%.

<!-- ==================================================================== -->
<!-- Sprint 6 leakage audit: verbatim from s6/LEAKAGE.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s6/LEAKAGE.md`. Not edited.*

# Leakage audit for the s6 learned ranker

Three disciplines are enforced directly and are easy to check:

1. **Pools.** A target's candidate pool is drawn from peptides outside its identity fold
   plus `distogram._fold_fragments(fold)`, and never contains the target itself. Identical
   for training and dev targets.
2. **Labels.** RMSD labels come only from *training* peptides' own natives. Dev natives are
   read in `s6/evaluate.py` after every score is computed, for reporting only.
3. **Ranker fitting.** `fit_folds` trains one model per identity fold, excluding all
   training targets in that fold. A dev target in fold *f* is scored by the model fitted
   without any fold-*f* training target.

## A fourth, indirect pathway - stated because it is real

The feature extractors are leave-fold-out *with respect to their own target*. For a dev
target in fold *f*, `Distogram.for_target` loads the fold-*f* model, which excluded fold *f*
- so the dev target's own features are clean.

But a **training** target in fold 2 gets its features from the fold-2 distogram, and that
model was trained on folds 0, 1, 3, 4 - which includes the dev peptides in those folds.
So training-target features are computed by a model that has seen some dev structures, and
the ranker is fitted on those features.

**Why this cannot inflate a specific dev prediction.** The contamination is
target-agnostic: a dev peptide in fold 0 contributes ~1 of ~6,630 training chains to the
fold-2 distogram, and its influence appears identically in the features of *every* training
target in fold 2, regardless of which dev target is later scored. There is no pathway by
which the ranker could learn something specific to the dev target it is asked about, because
nothing in its training data is indexed by that target.

**Why it is not simply dismissed.** It could in principle shift the learned
feature→RMSD mapping slightly, in a direction that flatters dev evaluation as a class. The
fully clean alternative is to recompute training-target features five times, once per
excluded fold, which costs 5x the pool build (~4 hours) and was judged not worth it against
a third-order effect.

**The control that tests it.** `s6/train.py` reports ranker performance on *held-out
training targets* (targets in the excluded fold, which share the same indirect pathway)
alongside dev performance. If the ranker were exploiting dev-specific information, dev
performance would exceed held-out-training performance. If the two match, the pathway is
inert. This is reported, not assumed.

## Not done anywhere

- No benchmark (60-target) structure is read by anything in `s6/`.
- No reference RMSD influences any deployed score; `y` enters only as a training label from
  training natives, and as a reporting quantity on dev.
- Oracle quantities (pool best, recall, in-band membership) are diagnostics and are never
  inputs to a selector.

<!-- ==================================================================== -->
<!-- Sprint 7: verbatim from s7/FINDINGS.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s7/FINDINGS.md`. Not edited.*

# Sprint 7 findings

## 1. The inter-generator agreement lead is a difficulty thermometer - closed

Sprint 6 reported that agreement between two structurally different generators predicts
accuracy at **rho +0.594** with no native access, and I called it "the first quantity in
this project that predicts whether a prediction is good". It does not survive scrutiny.

**It shrinks with n.** On all 24 dev targets rather than 12: **rho +0.367, 95% CI
[-0.087, +0.694], permutation p = 0.075.**

**One of the two numbers was a mathematical artefact.** The companion figure, rho +0.699
against |error_a - error_b|, must be discarded outright: the triangle inequality forces
`|e_a - e_b| <= d <= e_a + e_b`, so the agreement distance `d` is *algebraically* bounded by
that quantity. Verified numerically on 200 random pairs. I reported a constraint of the
metric as a discovery.

**Partial correlation removes the rest.** The confound is `dist_sd` - the distogram's own
predicted-distance uncertainty, which is sequence-derived, native-free, and already
available. It correlates +0.693 with the error and +0.590 with the agreement. Controlling
for it:

| generator pair | raw rho | + dist_sd | + length, spread, dist_sd |
|---|---|---|---|
| blosum / VQE (the original) | +0.37 | **-0.07** | -0.10 [-0.56, +0.41] |
| blosum / torsion | +0.52 | +0.10 | +0.14 |
| VQE / torsion | +0.46 | +0.05 | -0.08 |
| VQE seed0 / VQE seed1 | +0.33 | -0.16 | -0.13 |
| *blosum / random (shared library)* | *+0.77* | *+0.50* | *+0.48* |

**Every cross-mechanism pair collapses to zero.** The two that survive share a candidate
library, so their agreement partly measures whether the same window set returned the same
fragment - a property of the pool, not independent evidence.

**And the mechanism hypothesis fails its own prediction.** If agreement worked because
independent generators err independently, more *dissimilar* generators should agree more
informatively. Over 21 generator pairs, rho(dissimilarity, informativeness) = **-0.016**,
and pool-sharing pairs score *higher* (+0.597) than unrelated-mechanism pairs (+0.527).

Committees make the raw number look better and add nothing: raw +0.367 (K=2) -> +0.524
(K=4) -> +0.539 (K=7) while the partial stays flat at -0.10 / +0.02 / +0.11.

**Selection: nothing moved.** Union of 7 generators, 24 targets. Pool best improved
1.465 -> **1.370** by adding five generators; selected did not move - best arm 3.230 against
the distogram's 3.344, -0.11 A on 7/24. Oracle pick over the union would be 2.764.
These become selection approaches 20 through 25 in the same 3.27-3.44 pile.

**What is left of it.** `dist_sd` is a real native-free confidence signal at rho +0.693
against the error - but it is the distogram's own uncertainty, which Sprint 5 already
measured at rho +0.72 against MAE. The agreement work rediscovered an existing signal
indirectly and attributed it to consensus.

Methodological note worth keeping: rank-linear partialling does not fully remove an additive
multi-covariate confound (a constructed pure confound at raw +0.955 still shows partial
+0.522), so a stratified conditional permutation test was added, and its own limit pinned -
with 3 strata a null is strong evidence and a positive is weak. Both results here are null,
which is the favourable direction for this conclusion.

## 2. Thirty times more structural data makes the predictor monotonically worse - closed

The distogram trains on 6,630 chains (92,187 residues, ~638 K pairs). `prots/` holds
**13,751 parseable structures, 2,318,823 residues** - 25x the residue count - and for a
9-15mer only pairs with |i-j| <= 14 matter, which are just as available inside a 100-residue
protein as in a peptide. Nobody had ever changed the amount of data; every previous attempt
attacked the model or its input representation. So this was the one untested axis.

Pair count swept at full source diversity (all 13,751 sources in every arm, so the arms
differ in sample size and nothing else). 600 K is placed deliberately at the production
distogram's own ~638 K pairs as a like-for-like control:

| training pairs | dev MAE | training loss (final) |
|---|---|---|
| 200,000 | 2.385 | 2.8600 |
| 600,000 | 2.630 | 2.8333 |
| 2,000,000 | **3.358** | - |
| *production, 638 K peptide pairs* | *2.219* | - |

**Monotonically worse, and worse than production at matched pair count.** More data is not
merely useless here, it is actively harmful, and the effect is large: 3x the data costs
0.25 A, 10x costs 0.97 A.

**The mechanism is visible in the two losses moving in opposite directions.** Training loss
*falls* (2.8600 -> 2.8333) while dev MAE *rises* (2.385 -> 2.630). This is not overfitting in
the usual sense - the model has 250 K parameters against 2 M examples - it is distribution
shift. A fragment inside a folded protein is held in place by contacts outside the window, so
its conformational distribution is more compact and more regular than an isolated peptide's;
`fragment_db`'s own docstring records this. The asymptotic target of training is therefore the
*fragment* conditional mean, which is the wrong target. The more faithfully the model reaches
it, the further it lands from peptides.

**A caveat stated up front in the experiment's docstring, now resolved.** A null result would
have been ambiguous between "data does not help" and "this data is the wrong distribution".
A *monotonic degradation* is not ambiguous: insufficient sample size plateaus, it does not
reverse. The distribution is the problem.

**Note also how fast the training loss plateaus** - 2.92 at epoch 3, 2.86 at epoch 11, and
essentially the same floor at every data volume. The local-composition features are close to
exhausted almost immediately: the model reaches the conditional mean and stops, because there
is little target-specific signal left in them to fit. That is the same diagnosis finding 3
reaches from an entirely independent direction.

**Status: the data-limited hypothesis is dead.** This does not merely fail to help, it pins
the ceiling to the input representation and the training distribution rather than to sample
size - a stronger and more useful claim than the one previously recorded, which was that
"the sequence signal is the ceiling" inferred from a representation ablation alone.

## 3. The predictor is not shrunk, it is uncorrelated - post-hoc de-biasing is dead

`s7/debias.py`, tests `s7/test_debias.py`, results `s7/debias_calib.json`,
`s7/debias_tune.json`, `s7/debias_dev.json`.

### The hypothesis

The audit (`s7/audit.py`, stage `prior`) established that the same distance MAE buys
selection anywhere between 2.09 and 3.83 A depending only on the SHAPE of the error, and
that the shipped distogram sits on the *shrinkage* curve - errors shaped like regression
toward the average peptide's distance profile - essentially exactly. That predicted a
cheap fix: rescale the predicted DEVIATION from the pool's own mean profile and walk the
predictor back down the shrinkage curve, no retraining, from a ~ 1.0 (3.83 A) toward
a ~ 0.3 (2.21 A). Four arms were built to test it, all scoring off cached pool distance
matrices so the shipped Bayes-risk score is reproduced bit-for-bit at the identity
setting.

### The discriminator says the mechanism is wrong

Regressing true deviation `d_true - Dbar` on predicted deviation `expected - Dbar` over
8,549 held-out pairs from 126 fold-disjoint tuning targets:

| quantity | value |
|---|---|
| pooled calibration slope | **+0.376** |
| correlation | **+0.283** |
| sd(predicted deviation) | 2.584 A |
| sd(true deviation) | 3.437 A |
| variance-matched gain | 1.330 |
| per-target median slope | +0.186 (6/126 above 1) |

A slope **below** 1 is the opposite of what shrinkage predicts. The predicted deviations
are not small - at 2.58 A against a true 3.44 A they are 75% of full amplitude - they
simply **point the wrong way**. The decomposition is exact: slope = corr x
sd_true/sd_pred = 0.283 x 1.330 = 0.376. It holds in every separation shell (slope
+0.29 to +0.61, correlation +0.24 to +0.34), so it is not an artefact of one band.

The two principled single gains therefore disagree by a factor of 3.5: variance matching
says amplify by 1.33, least squares says shrink by 0.376. When a diagnosis implies a
correction and its two natural estimators point in opposite directions, the correction is
not identified, and that was visible before a single ranking number was computed.

### And the arm table confirms it, in both directions

126 tuning targets, K = 500, pool best 1.711 A, per-target sd of selected RMSD 1.646 A,
so **SE on the mean is 0.147 A** (dev's is 0.354 A). Paired differences against the
shipped baseline, 95% t intervals:

| arm | selected | vs base | 95% CI | W/L | MAE | in-band rho |
|---|---|---|---|---|---|---|
| sharpW0.25 (best) | 3.396 | **-0.058** | [-0.150, +0.033] | 30/18 | 2.434 | +0.127 |
| sharpF0.25 | 3.405 | -0.049 | [-0.099, +0.001] | 20/14 | 2.434 | +0.129 |
| gain 0.75 | 3.431 | -0.023 | [-0.130, +0.084] | 40/44 | 2.202 | +0.072 |
| **base (shipped)** | **3.454** | +0.000 | - | - | 2.339 | +0.126 |
| gainVM (1.330) | 3.509 | +0.055 | [-0.033, +0.143] | 33/44 | 2.660 | +0.124 |
| gainOLS (0.376) | 3.530 | +0.076 | [-0.106, +0.258] | 53/69 | 2.163 | -0.026 |
| pmi lam=0 (pure NLL) | 3.547 | +0.093 | [-0.055, +0.241] | 49/40 | 2.339 | **+0.182** |
| gain 1.5 | 3.582 | +0.128 | [+0.025, +0.231] | 35/56 | 2.868 | +0.125 |
| gain 2 | 3.723 | +0.269 | [+0.141, +0.397] | 36/73 | 3.572 | +0.124 |
| **gain 0 (pool mean, no sequence)** | **3.734** | +0.280 | [+0.036, +0.523] | 52/73 | 2.281 | -0.098 |
| gain 3 | 3.779 | +0.325 | [+0.155, +0.496] | 43/75 | 5.176 | +0.131 |
| gain 4 | 3.899 | +0.445 | [+0.236, +0.654] | 42/79 | 6.882 | +0.141 |

Full table (28 arms) in `s7/debias_tune.json`. The gain arm is **monotonically worse in
both directions away from 1.0**: amplifying costs 0.13 / 0.27 / 0.33 / 0.45 A at g = 1.5 /
2 / 3 / 4, and the least-squares-optimal shrink costs 0.08 A. There is no interior optimum
to find. Sharpening is the only family that moves the mean the right way and its best
interval still contains zero.

**The single most damning line is `gain 0`.** That arm throws the sequence away entirely
and ranks by the pool's own mean distance profile. It selects 3.734 A against the shipped
predictor's 3.454 A. **Everything the learned distogram knows about a specific sequence is
worth 0.28 A of selected RMSD**, against a 1.74 A gap to the pool's own best member.

Two secondary observations worth keeping. The pure-NLL arm has the **best in-band Spearman
of any arm (+0.182 against the baseline's +0.126) and a WORSE selected RMSD (3.547)** -
another instance of in-band rho and the decision metric disagreeing, so the protocol of
steering on selected RMSD alone earned its keep. And MAE is uninformative throughout: the
best-MAE arm (gainOLS, 2.163) ranks 16th of 28 on selection, exactly as predicted.

### The one dev number

The single best arm was carried to the 24 dev targets once, with no iteration:

| arm | selected | MAE | vs base | 95% CI | W/L/tie |
|---|---|---|---|---|---|
| base (shipped) | **3.475** | 2.219 | - | - | - |
| sharpW0.25 | **3.423** | 2.302 | -0.052 | [-0.158, +0.053] | 5/2/17 |

Baseline reproduction: **3.475 A at MAE 2.219**, against the 3.470 / 2.219 recorded in
`s5/inband.json` for the shipped scorer. The 0.005 A residual is tie-break noise in the
pool's BLOSUM sort (the audit sorts stably, `s5/inband.py` does not); everything
downstream is trustworthy. Pool best 1.539 A, oracle 1.782 A.

### Read

**This is a real negative and it closes the post-hoc route.** The prediction was
2.2-2.9 A; the measured result is 3.42 A, a statistically indistinguishable 0.05 A from
the shipped 3.475 A, and the calibration slope says why: the bias is not the invertible
kind. Post-hoc de-biasing can rescale a systematically attenuated signal; it cannot rotate
a signal that is only 28% correlated with the truth into alignment. Amplifying a
deviation that points the wrong way just amplifies the error, which is exactly the
monotone degradation the gain sweep shows.

The audit's placement of the prior on the shrinkage curve stands as a *description* of the
damage - a low-correlation full-amplitude error and a shrunken error destroy ranking about
equally at matched MAE - but it is not the *mechanism*, and correcting the second-moment
symptom does not touch the first-moment cause. Combined with finding 2 (more data makes it
monotonically worse) and the earlier representation ablations, the remaining lever is the
input representation and the training distribution: nothing that consumes the current
predictor's output differently is going to recover 1.7 A.

### Correction of record

The shrinkage diagnosis in this finding was **mine**, promoted from the audit's noise sweep
and used to commission this experiment. The experiment refuted it. Recording that with the
same prominence as the claim: I inferred a *mechanism* (bias toward the average peptide) from
what was only a *description* (where the predictor sits on a matched-MAE curve), and the
calibration slope - which I had named in advance as the discriminator - came back at +0.376
when shrinkage requires it well above 1. Note too that no synthetic arm actually reproduces
the shipped predictor: at matched MAE, iid gives 2.26, correlated 2.68, shrinkage 3.83, and
the real predictor 3.47. The honest description is error that is correlated across pairs, so
it does not average out over the ~30-90 pairs in a score, and only 28% aligned with truth.

## 4. Chirality is a red herring on the generative path - closed

`s7/audit.py` proved an exact bound: **reflection leaves the CA-CA distance matrix
EXACTLY unchanged** while moving CA-RMSD to the native by 3.075 A on average, and on 24/24
dev targets a mirrored native scores strictly better under `distogram.score` than every
real pool member. Reproduced here from scratch, per target, with a stricter instrument:
distance-matrix max |change| **0.000e+00 A**, distogram score max |change| **0.000e+00**,
mirror CA-RMSD **3.075 A**. The bound is real and exactly as stated.

The hypothesis was that the GENERATIVE path (`s6/pipeline.py` arms B and C), whose CVaR-VQE
optimises precisely that mirror-blind quantity, drifts into left-handed conformations at zero
objective cost. **It does not, and it structurally cannot.**

### The handedness distribution runs the wrong way

Three independent native-free measures (`s7/chirality.py`): a chirality pseudoscalar `G`
(mean normalised signed volume of consecutive CA bond triples, exactly negated by
reflection), the right-handed fraction of *helical* CA virtual dihedrals, and the
positive-phi (left-handed Ramachandran) rate split by glycine. 24 dev targets, 500
retrieval + 96 VQE candidates each.

| population | G | helix RH frac | phi+ all | phi+ non-Gly | flagged LH |
|---|---|---|---|---|---|
| **native** | +0.291 | 0.900 | 0.1060 | **0.0771** | 0.208 |
| A retrieval pool | +0.338 | 0.937 | 0.0670 | 0.0535 | 0.058 |
| **B VQE candidates** | +0.444 | **0.990** | 0.0715 | **0.0295** | 0.025 |
| **C VQE candidates** | +0.512 | **0.970** | 0.0602 | **0.0187** | 0.000 |

Reference: 12,623 training-peptide residues give positive-phi 5.66% non-Gly, 42.07% Gly.

The VQE candidates are **more right-handed than the retrieval pool and more right-handed
than the natives themselves** - 0.0295 non-Gly positive-phi against the natives' 0.0771 and
real protein's 0.0566, and 99.0% right-handed helices. 2.5% of arm-B and 0.0% of arm-C
candidates are flagged substantially left-handed, against 20.8% of the natives. There is no
left-handed drift to fix; if anything the generator is too right-handed.

### Why: the cheapest possible fix was already in force

`s5/assembly.assemble_batch` in mode "last" **imports real (phi, psi) values wholesale from
real fragment windows and never manufactures a torsion**. Reaching a left-handed
conformation requires the LIBRARY to offer positive-phi torsions, and the library is made of
real protein: 3.04% non-Gly positive phi, and only 0.87% of all its torsions are a
glycine positive-phi value landing on a non-glycine target position - the one mechanism by
which sequence retrieval could import a torsion an L-amino acid could not adopt. The task's
own guess was right: *the library never offered left-handed segments in the first place.*

### The mirror region is unreachable even when the search is paid to go there

The above could be luck, so it was settled adversarially: the same CVaR-VQE, same library,
same encoding, objective replaced by `+chirality_scalar` so the search is explicitly
**rewarded for building the most left-handed structure the library can express**.

* most left-handed `G` reached: **-0.379** (the mirrored native sits at -0.291, so the
  search does overshoot on this scalar);
* but the non-Gly positive-phi rate only rises **0.030 -> 0.088**, still within sight of
  real protein's 0.057 - torsionally it stays L-amino-acid-legal;
* closest CA-RMSD to the **mirrored** native: **4.929 A**, against 5.031 A to the real
  native. It gets no nearer the mirror than it does the original.

So `G` is not a proxy for "is the mirror image": a chain of right-handed backbone torsions
can still coil left-handed globally, and that is what the adversary finds. The actual
mirror image of a real protein requires positive phi nearly everywhere, and the library
cannot express it. **The distance matrix's blindness to reflection is inert on the
generative path for the same structural reason it is inert on the retrieval path.**

### The oracle bound prices the lead at ~0

DIAGNOSTIC ONLY, reads the native, lives in `oracle_reflection_bound` and is called from no
deployable path. Reflect every candidate, keep whichever of the two is closer to the native
- exactly one bit of native information per candidate, unobtainable natively because the
score of a pair is provably identical.

| arm | pool best | selected |
|---|---|---|
| A | 1.523 -> 1.486 (**-0.037 A**, CI [-0.082, +0.000], 3W/0L/21T) | 3.398 -> 3.321 (-0.077 A, 7W/0L/17T) |
| B | 1.465 -> 1.428 (**-0.037 A**, CI [-0.082, +0.000], 3W/0L/21T) | 3.310 -> 3.139 (-0.171 A, 6W/0L/18T) |
| C | 1.517 -> 1.480 (**-0.037 A**, CI [-0.082, +0.000], 3W/0L/21T) | 3.351 -> 3.287 (-0.063 A, 6W/0L/18T) |

**Pool best moves by 0.037 A on 3 of 24 targets** with a CI touching zero. The larger
selected number is a small-screen artefact: arm B's -0.171 A is **73% carried by two targets**
(1KVG -1.60, 7K1M -1.41); excluding them the mean gain over the other 22 is **-0.050 A**, and
the responsible targets are not the same across arms. Note also that the mirror is *closer*
than the original for only **33-36%** of candidates - below the 50% a chirality-agnostic
generator would give, i.e. the candidates are systematically on the correct side.

And the selected gain is **unattainable in principle, not merely in practice**: it requires
adding the mirror of a right-handed candidate to the pool, which is a left-handed
non-protein. No chirality term produces it; only abandoning L-amino-acid stereochemistry does.

### Part 2 run anyway, and null with the mechanism verified

A calibrated chirality term was added to the arm-B objective and genuinely optimised by the
VQE (nothing stubbed, no post-filtering; selection still uses the unmodified distogram
score). Two parts, both structure-only and leakage-free by construction: a Ramachandran NLL
`-log P(sign(phi_i) | aa_i)` from a per-amino-acid table estimated on training peptides
**with the target's identity fold excluded**, plus the left-handed fraction of helical CA
quadruples. Identical targets, library, seeds and retrieval pool.

| w | pool | selected | VQE-only pool | VQE-only selected |
|---|---|---|---|---|
| 0.5 | +0.001 A (0W/1L/23T) | +0.018 A [-0.020, +0.060] (4W/6L/14T) | -0.019 A (9W/8L/7T) | +0.045 A (7W/8L/9T) |
| 2.0 | +0.001 A (0W/1L/23T) | +0.007 A [-0.156, +0.145] (4W/6L/14T) | -0.059 A (12W/5L/7T) | +0.008 A (8W/8L/8T) |

**The term demonstrably works and buys nothing.** It moves the handedness exactly as
designed - non-Gly positive-phi 0.0295 -> 0.0031 (w=0.5) -> 0.0000 (w=2), helical
right-handed fraction 0.990 -> 1.000 - and every accuracy CI straddles zero. This is the
strongest form of a negative: the intervention was verified to act on its target quantity
and still changed nothing downstream. Worth noting that at w=2 it drives positive phi to
*exactly zero*, below real protein's 5.66%, so pushing harder makes the ensemble less
protein-like, not more.

### Verdict

**Chirality is a red herring for the generative path.** The mirror-blindness of the distance
objective is a genuine and exactly-provable defect of the *representation*, but it prices out
at 0.037 A of pool best on 3/24 targets, and the region it fails to police is one the
fragment-assembly parameterisation cannot enter. The Sprint 3 chirality defect and the
"optimising the objective harder makes structures worse" pathology (1KVG: objective
1.564 -> 1.495 while CA-RMSD went 2.276 -> 3.107) both need a different explanation: the
objective is not being fooled by handedness, it is being fooled by something the distance
matrix *can* see. That is where the remaining error is.

`s7/chirality.py` (`measure` | `reach` | `fix`), `s7/test_chirality.py` (15 tests, all
passing: exact distance-matrix invariance, sign flip of all three measures under reflection,
natives right-handed, oracle bound can only improve).

## 5. The objective is misspecified: the native is not its minimum

`s7/whatfools.py`, results `s7/whatfools_native.json`, `_shells.json`, `_sweep.json`,
`_fool.json`. 126-target tuning instrument.

Four sprints have treated the selection gap as a search or ranking problem. It is neither.
Scoring the NATIVE structure itself against the pool under the shipped objective (an ORACLE
DIAGNOSTIC, isolated from any deployable path):

| quantity | value |
|---|---|
| targets where the native is the objective's argmin | **3 / 126** |
| native's mean percentile in the pool's score distribution | **36.8** |
| native's median percentile | 32.8 |

**On 123 of 126 targets at least one candidate scores better than the truth, and on average
~37% of the pool beats it.** No amount of better searching or better selecting can recover a
structure the objective does not rank first. This is the direct explanation of the
"optimise harder, get worse" pathology (1KVG: objective 1.564 -> 1.495 while CA-RMSD went
2.276 -> 3.107): the search was working correctly and the target was wrong.

**My leading hypothesis for WHY was refuted.** I predicted the score is swamped by
short-separation pairs whose distances are fixed by backbone geometry and therefore cannot
discriminate, with the informative long-range pairs down-weighted by the `1/(sd+0.5)^gamma`
rule. The sweep says the sd weighting is *helping*, not hurting:

| arm | selected |
|---|---|
| gamma 0.5 | **3.429** |
| gamma 1 (shipped) | 3.454 |
| gamma 2 | 3.459 |
| longfrom4 (long-range only) | 3.472 |
| gamma 0 (no sd weighting) | 3.528 |
| gamma -0.5 (up-weight uncertain) | 3.535 |
| longfrom6 | 3.544 |

Removing the sd weighting makes it worse; inverting it makes it worse still; ranking on
long-range pairs alone makes it worse. Best arm is -0.025 A against shipped, inside the
0.147 A SE. **Re-weighting cannot fix a misspecified objective**, which is consistent with
the recorded caution that per-separation weights improved a proxy once without transferring.

There is a real tension worth keeping: per shell, the native's percentile is 33.2 / 35.4 /
42.0 / 40.7 / **16.2** across SHELLS 2-3, 4-5, 6-8, 9-13, 14-40. The long-range shell ranks
the native best by a wide margin, yet ranking the POOL on long-range pairs alone is worse.
The long-range channel is the most *discriminating* and the most *noisily predicted* at once,
so its value is real but unreachable through the current predictor.

## 6. Pool size has an interior optimum, and selection tracks the pool MEAN not its BEST

`s7/poolsize.py`, results `s7/poolsize_kcurve.json`, `_prune.json`. Same 126 targets.

| K | selected | pool best | pool mean |
|---|---|---|---|
| 25 | 3.439 | 2.350 | 4.177 |
| 50 | **3.399** | 2.177 | 4.232 |
| 100 | 3.425 | 1.970 | 4.316 |
| 250 | 3.461 | 1.808 | 4.393 |
| 500 | 3.454 | 1.711 | 4.453 |
| 1000 | 3.483 | 1.568 | 4.522 |
| 2000 | 3.520 | 1.504 | 4.596 |

**The predicted inversion is real.** Under an ORACLE ranker bigger pools help monotonically
(K=100/500/2000 -> 1.971/1.782/1.631). Under the REAL ranker they hurt: pool best improves
0.85 A from K=25 to K=2000 while selected gets 0.08 A WORSE. There is an interior optimum
near K=50.

**But the size of the effect is the finding, not the direction.** K=50 beats K=500 by only
0.055 A, inside the 0.147 A SE. Pool construction is not a lever.

**The decomposition is the valuable part.** Selected RMSD tracks the pool MEAN (4.177 ->
4.596 as K grows, selected 3.439 -> 3.520, moving together) and is *anti*-correlated with the
pool BEST (which improves throughout). Pruning shows the same signature from the other side:
the best pruning rule improves the pool mean from 4.232 to 3.712 and moves selected only
3.399 -> 3.420. **A selector that tracked quality would follow the best; this one follows the
mean.** That is a quantitative statement of near-powerlessness, and it agrees exactly with
finding 3's `gain 0` result (sequence discarded entirely costs only 0.28 A) and with the
architecture arms sitting at the zero-information point.

**Not measured:** the BLOSUM-vs-random retrieval-key comparison (`poolsize_retr.json` is
empty). The agent died in an API outage before running it, so the recorded claim that 100
random fragments beat the pipeline's retrieval remains unverified on this instrument.

## 7. A ceiling worth stating plainly

On the same 126 pools, ranking by the target's OWN distance matrix (perfect distance oracle,
diagnostic only) selects **1.994 A**, in-band rho +0.579, against pool best 1.711 and the
shipped predictor's 3.454.

So a *perfect* distance predictor, on the pools this pipeline actually builds, lands at
1.99 A. The sub-2.0 A goal is therefore not merely a predictor problem: it sits exactly at
the ceiling of the distance-ranking approach, and any margin below 2.0 A has to come from
better pools or a scoring channel that is not pairwise CA distances.

## 8. Training against ranking makes it worse, and the predictor is ~= typicality

`s7/rankloss.py`, results `s7/rankloss.json`. n=234 leave-fold-out targets, SE 0.110 A.

I argued this was the most promising remaining lever, on the reasoning that finding 5 showed
the hand-built distance objective is misspecified and a loss trained directly on candidate
RMSD ordering would bypass that rather than inherit it. **It is significantly worse than plain
regression.**

| arm | selected | d vs base | 95% CI | W/L | in-band | MAE |
|---|---|---|---|---|---|---|
| bins (soft-CE distogram) | 3.508 | -0.054 | [-0.164, +0.056] | 119/94 | +0.059 | 2.284 |
| **base (plain regression)** | **3.562** | - | - | - | -0.010 | 2.222 |
| resid (deviation param.) | 3.620 | +0.058 | [-0.052, +0.167] | 108/106 | -0.017 | 2.314 |
| crps | 3.635 | +0.073 | [-0.038, +0.185] | 101/114 | -0.017 | 2.470 |
| rank_p | 3.656 | +0.094 | [-0.025, +0.213] | 101/108 | -0.005 | 2.369 |
| antishrink | 3.657 | +0.095 | [-0.008, +0.198] | 93/109 | -0.029 | 2.701 |
| rank_a | 3.657 | +0.095 | [-0.015, +0.205] | 99/119 | -0.015 | 2.384 |
| rank_w | 3.710 | +0.148 | **[+0.035, +0.260]** | 89/126 | -0.013 | 2.360 |
| **rank** (pairwise ranking loss) | **3.773** | **+0.212** | **[+0.079, +0.344]** | 94/126 | +0.002 | 2.731 |

Both `rank` and `rank_w` are worse with CIs excluding zero. Nothing beats the baseline: the
best arm (`bins`, -0.054) is null. **The training objective was not the problem.** Note the
anti-shrinkage arm is also null-to-worse (+0.095), which independently confirms finding 3 from
the training side: shrinkage was never the defect.

**The reference row is the more important result.** On the same pools:

| reference | selected |
|---|---|
| pool best | 1.628 |
| perfect distance oracle | 1.946 |
| **typicality** (rank by agreement with the pool's own mean profile) | **3.586** |
| trained predictor (base) | **3.562** |
| random | 4.289 |

**The trained predictor beats typicality by 0.024 A.** Typicality uses no sequence and no
model -- it asks only "how average is this candidate". This is a sharper version of finding 3's
`gain 0` result (sequence discarded costs 0.28 A) and it agrees with the architecture arms
sitting at the zero-information point (finding, `s7/pairnet2.py`): a per-pair MLP on one-hot
selects 3.829 against a no-sequence baseline's 3.734.

Three independent measurements now say the same thing: **the deployed predictor contributes
almost nothing over ranking by typicality.** That is not a calibration problem, a training-loss
problem, or a representation problem -- it is an information problem.

**The oracle here is 1.946 A**, consistent with the 1.994 A measured on the other instrument.
Two independent instruments agree that perfect distance knowledge on these pools lands at
~1.95-2.0 A.

## 9. The triangle operator improves the diagnostics significantly and the outcome not at all

`s7/pairnet2.py`, results `s7/pairnet2_select.json`. 126-target tuning instrument, K=500,
fixed one-hot input and fixed data (6,630 chains / 552,199 pairs), so arms differ only in
architecture. `tri` covers 96 targets pending the last fold.

| arm | params | n | pool | SELECTED | +-SE | rho_gl | rho_band | slope | corr | (MAE) |
|---|---|---|---|---|---|---|---|---|---|---|
| mlp (per-pair control) | 41,985 | 126 | 1.711 | 3.808 | 0.159 | +0.390 | -0.060 | 0.102 | 0.028 | 2.347 |
| pair0 (pair tensor, transitions only) | 31,825 | 126 | 1.711 | 3.795 | 0.157 | +0.389 | -0.057 | 0.117 | 0.032 | 2.352 |
| **tri (+ triangle multiplicative update)** | 62,929 | 96 | 1.719 | **3.747** | 0.174 | +0.463 | **+0.001** | **0.231** | **0.146** | 2.399 |
| SHIPPED distogram (ESM input) | - | 126 | 1.711 | **3.454** | 0.147 | +0.568 | +0.126 | 0.376 | 0.283 | 2.339 |

Paired against the per-pair MLP control:

| arm | d SELECTED | 95% CI | p | d rho_band | 95% CI |
|---|---|---|---|---|---|
| pair0 | -0.013 | [-0.071, +0.047] | 0.672 | +0.003 | [-0.002, +0.009] |
| **tri** | **-0.017** | [-0.268, +0.245] | **0.900** | **+0.064** | **[+0.013, +0.115]** |

**Both halves of this are real and they disagree.** The triangle operator produces a
statistically significant in-band rho gain (+0.064, CI excludes zero) and a large calibration
gain (deviation correlation 0.028 -> 0.146, slope 0.102 -> 0.231). It produces **no selection
gain at all** (-0.017, p 0.900). `pair0` ~= `mlp` on everything, which isolates the effect to
the triangle multiplicative operator rather than to pair-tensor plumbing -- so the mechanism
attribution is clean; it is the payoff that is absent.

**This retires the calibration slope as a proxy for the thing we care about.** I promoted that
metric earlier in the sprint as the discriminator, and it did its job diagnostically -- it is
what refuted the shrinkage story in finding 3. But `tri` now improves slope, deviation
correlation AND in-band rho simultaneously while selected RMSD does not budge. Three
consecutive results say the same thing from different directions: better proxies do not buy
better structures. Judge arms on selected RMSD, full stop; use slope and in-band rho only to
explain a selection result, never to predict one.

**And the architecture arms do not reach the shipped model.** All three one-hot arms sit at
3.75-3.81 against the production ESM distogram's 3.454 on identical pools. Whatever the
triangle operator adds, it does not close the gap to the existing system, let alone to the
1.95-2.0 A distance-oracle ceiling.

Sprint 6's memory recorded pairnet as "corr 0.744 -> 0.825, MAE 2.04 -> 1.74" and treated it
as promising. On the correct instrument it is the eighth closed axis.

## 10. No available scoring channel ranks the truth first - Amber is worse, not better

`s7/amber_native.py` + `s7/dist_baseline.py`, results `s7/amber_native.json`,
`s7/dist_baseline.json`. Paired, matched pools (same targets, same BLOSUM retrieval, same
K=100 source, same top-40 scored subset, same ideal-geometry rebuilds). n=25 and continuing.

Finding 5 showed the distance objective is misspecified. All-atom Amber was the strongest
remaining candidate for a correctly-specified one: a genuinely different channel, an
architecturally mandatory component, and recorded in memory at in-band rho **+0.320** against
the distogram's +0.091 with the note that "Sprint 4's negative was a build confound".

| metric | Amber | distogram | diff | 95% CI | Amber wins |
|---|---|---|---|---|---|
| selected RMSD | 3.648 | **3.068** | +0.580 | [-0.002, +1.162] | 11/25 |
| in-band rho | 0.203 | 0.240 | -0.037 | [-0.219, +0.146] | 12/25 |
| **partial(rg) rho** | 0.238 | **0.561** | **-0.323** | **[-0.451, -0.195]** | **2/25** |
| **native percentile** | **51.0** | 33.8 | **+17.2** | **[+4.2, +30.3]** | 9/25 |
| native is the argmin | **0/25** | 2/25 | | | |

**Amber places the rebuilt native at the 51st percentile of its own energy distribution** --
worse than a coin flip among candidates, and the argmin on zero of 25 targets. It is
*significantly worse* than the distance objective at ranking the truth (+17.2 percentile,
CI excludes zero). So the answer to the decisive question is the unfavourable one: **no
scoring channel available to this project ranks the native first.** That is the most
important negative result of the program, because it means the selection gap is not a
selector problem, a search problem, or a predictor problem -- the objectives themselves do
not have their minima at the right structures.

**The +0.320 record does not survive a matched comparison with a compactness control.** Two
methodological points explain the discrepancy, and both were introduced before aggregating
rather than after:

1. **Pool matching.** The recorded distogram figure (+0.091 / +0.126) came from a K=500 pool;
   Amber's came from a small top-K subset. On a matched top-40 subset the distogram scores
   +0.240 in-band, not +0.091. At n=1 I briefly read Amber's +0.656 as ~5x the distogram --
   the same pool scores +0.739 for the distogram. Comparing scorers across pool constructions
   is invalid.
2. **The compactness control.** Amber's nonbonded + solvation terms reward compactness, and a
   compact native makes energy correlate with RMSD through radius of gyration alone. Partialling
   rg out costs Amber 0.203 -> 0.238 raw-to-partial in the mean but destroys its per-target
   advantage: it wins on 2/25 against the distogram's partial +0.561. This is the same confound
   that took finding 1's agreement lead from +0.594 to -0.07.

Honest limitations. (a) `peptide_db.Peptide` exposes only (ca, phi, psi, rebuild) -- no
experimental all-atom coordinates -- so the *deposited* native cannot be force-field scored at
all. The comparator is the native rebuilt from its own torsions through the same builder as
every candidate, which holds geometry provenance constant but shifts the structure by ~0.8 A
CA-RMSD. A physics energy is more geometry-sensitive than a distance score, so some of Amber's
native penalty may be rebuild distortion rather than misspecification; that component cannot be
separated with the data this repo stores. (b) Scoring is converged interaction-only
(`refine_coords`, K_MODERATE, steps=0, nonbonded+solvation) per the recorded recipe, not a free
minimisation. (c) n=25 of a planned 70; the two significant results have CIs excluding zero but
the selected-RMSD difference is borderline.

### Finding 10, FINAL at n=70 (the n=50 interim below is superseded)

The run continued to 50 paired targets. Every effect holds and selected RMSD becomes
significant:

| metric | Amber | distogram | diff | 95% CI | Amber wins |
|---|---|---|---|---|---|
| **selected RMSD** | 3.953 | **3.343** | **+0.609** | **[+0.110, +1.108]** | 18/50 |
| in-band rho | 0.207 | 0.182 | +0.025 | [-0.129, +0.179] | 22/48 |
| **partial(rg) rho** | 0.160 | **0.516** | **-0.356** | **[-0.462, -0.251]** | **4/50** |
| **native percentile** | **52.8** | 40.5 | **+12.4** | **[+2.8, +21.9]** | 18/49 |
| native is the argmin | **0/50** | 5/50 | | | |

Three of the four differences now have confidence intervals excluding zero, all against Amber.
The distogram is the argmin at the native on 5/50 targets; Amber on **zero**. Raw in-band rho
remains the one metric on which the two are indistinguishable (+0.025, CI spans zero) -- and
it is precisely the metric that does not survive the compactness control, where Amber wins
4/50.

**Settled:** all-atom Amber is not a correctly-specified objective for this problem and does
not rank candidates better than the distance prior once its compactness advantage is removed.
Combined with finding 5, neither of the project's two scoring channels places its minimum at
the native. The rebuild-distortion caveat above still applies to the native-percentile figure
and cannot be separated with the data this repo stores; it does NOT apply to the partial(rg)
or selected-RMSD results, which compare the two scorers on identical candidate structures.

### Finding 10 FINAL - all 70 paired targets

| metric | Amber | distogram | diff | 95% CI | Amber wins |
|---|---|---|---|---|---|
| **selected RMSD** | 4.037 | **3.546** | **+0.491** | **[+0.081, +0.901]** | 26/70 |
| in-band rho | 0.181 | 0.170 | +0.011 | [-0.127, +0.149] | 30/66 |
| **partial(rg) rho** | 0.159 | **0.479** | **-0.320** | **[-0.425, -0.216]** | **12/70** |
| **native percentile** | **54.3** (median 61.5) | 42.3 (median 41.2) | **+12.0** | **[+3.4, +20.6]** | 25/69 |
| native is the argmin | **0/70** | 5/70 | | | |

Complete and stable from n=25 through n=70. Amber's median native percentile is **61.5** --
on a typical target, roughly three-fifths of the candidate pool scores better than the truth
under an all-atom force field. It is the argmin at the native on **zero of 70** targets.

Both scoring channels available to this project are therefore misspecified, and the physics
channel is the worse of the two. The only metric on which they tie is raw in-band rho, which
is exactly the metric that does not survive the compactness control.

## 11. ESM DOES beat one-hot - on selection. The recorded ablation was judged by the wrong metric

`s7/repr_select.py`, results `s7/repr_tune.json`. 126-target tuning instrument, paired,
arms differing ONLY in the feature block (verified: the one-hot block is bit-identical to
`s5/dist_ablate`, the scoring path bit-identical to shipped `Distogram.score`).

This project recorded, in memory as `esm-adds-nothing-for-short-peptides`, that a 1280-d
language model TIES 20-d one-hot: one-hot 2.216 / pca32 2.133 / pca128 2.212 / raw 2.109,
all p >= 0.59. **Every one of those numbers is a distance MAE**, and this sprint established
that MAE is close to useless here (finding 3: the best-MAE arm of 28 ranked 16th of 28 on
selection). Re-judged on selected CA-RMSD:

| arm | selected | d vs onehot | 95% CI | p | W/L | in-band | corr |
|---|---|---|---|---|---|---|---|
| onehot (20-d) | 3.754 | - | - | - | - | -0.026 | -0.072 |
| phys (physicochemical, 25-d) | 3.755 | +0.001 | [-0.045, +0.047] | 0.972 | 9/11 | -0.026 | -0.075 |
| **pca32 (ESM, 52-d)** | **3.466** | **-0.288** | **[-0.484, -0.092]** | **0.0046** | **74/45** | **+0.113** | **+0.259** |

**The representation axis is NOT closed.** ESM buys 0.288 A of selected RMSD over one-hot,
p=0.0046 on 126 paired targets, winning on 74 and losing on 45 -- while being invisible to
MAE, the metric the original ablation used. Deviation correlation goes -0.072 -> +0.259.

**`phys` is the control that makes this specific.** A 25-d physicochemical encoding ties
one-hot to +0.001 A. So this is not "any richer encoding helps"; it is the language model
in particular.

**This also re-explains finding 9.** All three one-hot architecture arms sat at 3.75-3.81
against the shipped ESM distogram's 3.454 on identical pools, and I read that as the
architecture failing to reach the production system. It is not an architecture deficit at
all -- 0.288 A of that 0.30-0.35 A gap is the representation. The triangle operator's null
result stands, but the arms were handicapped by their input, and `tri` on ESM input has not
been measured.

**What this does and does not change.** It does NOT break the ceiling: pca32 at 3.466 is
statistically indistinguishable from the shipped distogram's 3.454 (which already uses ESM),
so this is an explanation of why the production system beats the one-hot arms, not an
improvement on it. The ~1.95-2.0 A distance-oracle ceiling (finding 7) is untouched. What it
changes is the recorded claim: the sequence representation carries real, measurable selection
signal that MAE cannot see, and "ESM adds nothing for short peptides" is false.

pca128 and raw 1280-d arms are still running.

### Finding 11 complete - all five representation arms, n=126 paired

| arm | selected | d vs onehot | 95% CI | p | W/L | corr | MAE |
|---|---|---|---|---|---|---|---|
| onehot (20-d) | 3.754 | - | - | - | - | -0.072 | 2.462 |
| phys (25-d) | 3.755 | +0.001 | [-0.045, +0.047] | 0.972 | 9/11 | -0.075 | 2.463 |
| **pca32 (52-d)** | 3.466 | **-0.288** | [-0.484, -0.092] | **0.0046** | 74/45 | +0.259 | 2.111 |
| **pca128 (148-d)** | **3.417** | **-0.337** | [-0.576, -0.099] | **0.0064** | 78/48 | +0.257 | 2.148 |
| **raw (1300-d)** | 3.453 | **-0.301** | [-0.560, -0.041] | **0.0251** | 73/52 | +0.242 | 2.304 |

**All three ESM arms significantly beat one-hot** (p = 0.005 / 0.006 / 0.025) by 0.29-0.34 A,
and are mutually indistinguishable (3.417-3.466, spread well inside the CIs). The
physicochemical control ties one-hot exactly. So the effect is: *a language model, any
reduction of it*, is worth ~0.3 A of selected RMSD over amino-acid identity.

**The MAE column is the sharpest indictment of the old instrument.** The published ablation
ranked `raw` BEST on MAE (2.109). Here `raw` has the WORST MAE of the three ESM arms (2.304)
and mid-pack selection, while `pca128` has mediocre MAE (2.148) and the best selection. MAE
does not merely fail to separate ESM from one-hot -- **it misorders the ESM variants against
each other.** Any conclusion in this repo resting on a distance-MAE comparison should be
regarded as unsupported until re-run on selected RMSD.

pca128's 3.417 is nominally better than the shipped ESM distogram's 3.454, but the difference
is far inside noise and pca128 is not a new system -- it is the same class of input the
production model already uses.

## 12. BLOSUM retrieval beats random fragments - the recorded claim is refuted

`s7/poolsize.py stage_retr`, results `s7/poolsize_retr.json`. 126 targets, paired, identical
window universe; only the retrieval KEY differs.

The memory `real-fragments-beat-the-lattice` recorded that "100 RANDOM real fragments pool at
1.969 A vs the whole pipeline's 2.355", and `structure-and-sequence-are-decoupled` recorded
that the best-matching window has only 12% sequence identity. Together they suggested BLOSUM
retrieval might be an actively worse pool key than random sampling. Measured directly:

| K | BLOSUM sel | random sel | d(rand-blosum) | W/L | BLOSUM best | rand best | BLOSUM oracle | rand oracle |
|---|---|---|---|---|---|---|---|---|
| 25 | 3.439 | 3.516 | +0.077 | 59/67 | 2.350 | 2.512 | 2.617 | 2.744 |
| 50 | **3.399** | 3.586 | +0.186 | 50/76 | 2.177 | 2.340 | 2.445 | 2.544 |
| 100 | 3.425 | 3.563 | +0.138 | 59/67 | 1.970 | 2.172 | 2.241 | 2.390 |
| 250 | 3.461 | 3.545 | +0.085 | 55/70 | 1.808 | 1.903 | 2.091 | 2.235 |
| 500 | 3.454 | 3.543 | +0.089 | 54/69 | 1.711 | 1.767 | 1.994 | 2.130 |
| 1000 | 3.483 | 3.575 | +0.091 | 48/74 | 1.568 | 1.660 | 1.874 | 1.998 |
| 2000 | 3.520 | 3.516 | -0.003 | 57/61 | 1.504 | 1.544 | 1.838 | 1.864 |

**BLOSUM wins at every K except 2000, where they tie**, and it wins on all four measures at
once: selected RMSD, pool best, pool mean (4.453 vs 4.817 at K=500) and the oracle ceiling.
The advantage shrinks as K grows, which is exactly what it should do -- at K=2000 the retrieval
is taking most of the available window universe anyway, so the key stops mattering.

**Retrieval is doing real work and the recorded claim is refuted.** The original comparison
was presumably not matched on window universe or pool size; on a matched, paired test sequence
similarity is a genuinely useful pool key despite the best-matching window sharing only ~12%
identity with the target.

**The oracle column is the most consequential number in this table.** Even at K=2000 with a
PERFECT scorer, BLOSUM pools cap at **1.838 A**. Larger pools do lower the ceiling
(2.617 -> 1.838 from K=25 to K=2000) but with strongly diminishing returns, while the REAL
selector gets monotonically worse over the same range (3.439 -> 3.520). Pool construction
cannot deliver sub-1.5 A, and cannot be exploited at all without a scorer that does not exist.

<!-- ==================================================================== -->
<!-- Sprint 8: verbatim from s8/FINDINGS.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s8/FINDINGS.md`. Not edited.*

# Sprint 8 findings - full architectural attack

Sprint 7 closed twelve axes and left one structural conclusion: neither of the project's two
scoring channels puts its minimum at the native, and perfect distance knowledge caps at
~1.95-2.0 A on the pools this pipeline builds. Sprint 8 attacks the architecture itself.

Instrument throughout: the 126-target tuning set, cluster-disjoint from BOTH the 24-target dev
set and the 60-target benchmark. Neither of those was touched for model selection.

## S8-1. The representation is not the bottleneck - the aggregation metric is

`s8/repr_ceiling.py`, `s8/repr_ceiling_*.json`, 15/15 tests. Pure oracle arithmetic: rank each
pool by agreement with the NATIVE under each representation and read off the ceiling.

Harness validity first: the `s7/debias_cache` npz files store distances and RMSD but not
coordinates, so the pools were rebuilt through the same code path and asserted against the
cached values (worst error 1.9e-6 over all 126). The incumbent arm reproduces sprint 7's table
exactly (2.617 / 2.241 / 1.994 / 1.838 at K=25/100/500/2000). SE 0.098 A at K=500.

| arm | selected | d vs distances | 95% CI | W/L | frac <2.0 |
|---|---|---|---|---|---|
| coords_rmsd (control = pool best) | **1.711** | -0.283 | [-0.349, -0.220] | 89/0 | 0.579 |
| coords_l1 after alignment | 1.727 | -0.267 | [-0.334, -0.204] | 82/6 | 0.571 |
| **local_allcoord (lossless)** | **1.945** | -0.049 | [-0.123, +0.027] | 55/35 | 0.524 |
| residue frames @0.5 | 1.961 | -0.032 | [-0.106, +0.043] | 50/40 | 0.524 |
| **CA distances (incumbent)** | **1.994** | - | - | - | 0.508 |
| CA distances, L2 | 2.006 | +0.013 | [-0.033, +0.058] | 15/13 | 0.508 |
| **6D distance + orientation** | **2.219** | **+0.225** | **[+0.091, +0.378]** | 44/62 | 0.484 |
| CA contacts @8 A | 2.529 | +0.536 | [+0.377, +0.707] | 18/83 | 0.413 |
| backbone torsions | 2.952 | +0.958 | [+0.735, +1.203] | 14/101 | 0.373 |
| CB contacts @8 A | 3.451 | +1.457 | [+1.134, +1.814] | 15/103 | 0.317 |

Ordering is stable at all four K.

**The decisive control.** `local_allcoord` expresses every CA in every interior local frame and
a test reconstructs the structure from it to <1e-9 A - it is **informationally lossless**.
Consumed by the same additive mean-|delta| statistic it selects 1.945, so **83% of the
incumbent's 0.283 A loss survives making the representation lossless.** Only 0.049 A is
attributable to what pairwise CA distances discard, and that is not significant.

The contrast pair is the whole finding: aggregating agreement **over pairs** gives 1.945-2.006
regardless of representation, lossy or lossless; aggregating deviation **over atoms after
superposition** gives 1.711-1.727 regardless of norm. **The 0.28 A is not missing information,
it is discrimination destroyed by averaging.**

**6D is significantly WORSE than plain CA distances at every K** (+0.225). Orientation does add
over its own CB-CB distance channel (2.159 vs a matched 2.532 control), but that channel is far
worse than CA distances and the sum never catches up. And the gain is not chirality: the
virtual-CB construction contains a cross product, so CB-CB distances are already
reflection-sensitive, making the matched control chirality-aware. Building a 6D predictor would
have been the worst available lead.

**Epistemic note the agent added, and it is correct:** any *injective* representation has an
information-theoretic ceiling of exactly pool best, since one could rank by lookup. "The ceiling
of representation R" is only ever the ceiling of R *plus its natural agreement metric*.

**Correction to the brief.** Native percentile is degenerate under an oracle reference - the
native scores exactly 0 and is trivially the argmin at the 0th percentile for every arm. That
test is meaningful only for a LEARNED scorer. The useful substitutes: the native's
ideal-geometry rebuild sits at the 0.0-0.5th percentile on every arm, so all arms are correctly
specified in that sense; and the genuinely best pool candidate sits at the 2.7th percentile
under distances - the ~14/500 candidates that outrank it are merely ~0.3 A worse.

## S8-2. Multimodality is real but modest: 2-3 populated basins, not 30

Same study. In the K=500 near-native band (pool best +1.5 A, mean 134 candidates, mean pairwise
2.52 A), average-linkage at 1.5 A gives 33 clusters (median 22) and **only 1/126 targets is a
single cluster**. But the mass concentrates: only ~2.7 clusters hold >=5% of the band and the
largest holds 41%.

So a single-hypothesis predictor is genuinely misspecified, and multi-hypothesis output is
warranted - at **2-3 hypotheses, not 30**. This sets the size of any ensemble or
quantum-selection state space.

## S8-3. One actionable lead, and it is priced honestly

Since consumption-by-superposition is what wins, the deployable form is: predict a *structure*
and rank by RMSD to it, rather than predict a *matrix* and rank by matrix agreement. Measured
under degradation (`s8/repr_ceiling_noisy.json`): under isotropic noise the advantage grows with
error (-0.283 at 0 A to -0.838 at 4 A), but under the realistic degradation - a real pool
structure standing in for the prediction - it is **-0.139 [-0.207, -0.078] at a 1.72 A
reference and dead by 3 A (+0.019 [-0.092, +0.126])**.

The shipped selector is at 3.454 and nothing in the pipeline predicts near 3 A, so **the lead is
real, worth ~0.1 A, and does not pay at current predictor quality.** It becomes live only if a
predictor reaches ~2.5 A.

**Consequence for the sprint:** abandon representation work. Do not build a 6D predictor, do not
switch to torsions or contacts. What remains is generation and recognition.

## S8-4. What the 3.454 A mean is actually made of - four classes, and a 28% opportunity

`s8/triage.py` (7 resumable stages), `s8/test_triage.py` (24 tests), `s8/triage_*.json`.
ORACLE DIAGNOSTIC throughout; **no target was removed from any mean** and a test asserts the
headline stays 3.454. The instrument reproduces sprint 7 exactly (selected 3.454, pool best
1.711 at K=500), asserted as a test.

### The benchmark is an NMR benchmark

115/126 of the instrument is solution NMR; 4 X-ray; 111/126 are multi-model depositions; 40/126
are excised from a larger protein. Dev(24) and benchmark(60) have the same composition (77/84
NMR), so this transfers. This is precisely the data AF2 excluded, on top of excluding peptides
under 16 residues.

**Intra-ensemble CA-RMSD to model 1 (n=111): mean 1.044 A, median 0.872, max 4.363.** 48/126
targets disagree with themselves by more than 1.0 A. Sprint 6's 0.956 A was typical, not an
outlier.

**But model 1 is not the wrong native, and this was checked rather than assumed.** Rescoring
every rebuilt pool against every deposited model (rebuild reproduces the cache to 1.2e-07 A):

| pool best measured against | mean |
|---|---|
| model 1 (what we score) | 1.778 |
| a randomly chosen model | 1.777 |
| the closest model (min over ~20) | 1.513 |
| the furthest model | 2.273 |

**Model 1 minus a random model is +0.001 A, CI [-0.030, +0.033].** The 0.265 A apparent gain
from scoring against the closest model is a pure order statistic over ~20 draws, not a
correctable bias. **Ill-posedness is real but it is not a free 0.3 A** - it is an irreducible
~1 A noise term on every per-target number.

### Achievability ceilings, stated for the first time

| K | pool best | has a sub-2 A candidate | has a sub-1.5 A one |
|---|---|---|---|
| 25 | 2.350 | 53/126 (42%) | 40/126 |
| 100 | 1.970 | 62/126 (49%) | 48/126 |
| **500** | **1.711** | **73/126 (58%)** | **53/126 (42%)** |
| 2000 | 1.504 | 87/126 (69%) | 64/126 |
| **whole library (~17k windows)** | **1.313** | **102/126 (81%)** | **77/126 (61%)** |

At the K the sprint reports on, **42% of targets have no sub-2 A candidate at all.** The shipped
pipeline returns sub-2 A on 27/126 (21%).

### The class decomposition

| class | n | % | mean sel | pool best | oracle | contributes | % of the mean |
|---|---|---|---|---|---|---|---|
| rankable_missed | 46 | 36.5% | 3.664 | 1.256 | 1.444 | 1.338 A | **38.7%** |
| **retrieval_limited** | **29** | 23.0% | 4.170 | 2.351 | 2.844 | **0.960 A** | **27.8%** |
| library_limited | 24 | 19.0% | 4.651 | 2.869 | 3.275 | 0.886 A | 25.6% |
| solved | 27 | 21.4% | 1.264 | 0.768 | 0.878 | 0.271 A | 7.8% |

Contributions sum to 3.454 exactly (enforced by a test, so nothing can be silently dropped).

**The 29 retrieval-limited targets are the largest actionable opportunity found in two sprints:
a sub-2 A window already exists in the library, at BLOSUM rank 576-21,292 (median ~2,006), and
K=500 never sees it.** Across all targets the best window sits at **BLOSUM percentile 39.3**
against random's 50.0, and is inside the top 500 on only **18/126**. Sequence similarity is
barely better than chance at finding the structurally closest window.

### What mean is attainable

| scenario | mean |
|---|---|
| actual, shipped pipeline | 3.454 |
| perfect selector on the tractable 73, actual elsewhere | **2.469** |
| perfect distance oracle everywhere (s7 finding 7) | 1.994 |
| perfect selector everywhere (pool best, K=500) | 1.711 |
| …credited only down to each target's ensemble noise | 1.813 |
| perfect retrieval AND perfect selection (whole-library best) | 1.313 |
| …credited only down to ensemble noise | 1.502 |

**Sub-2.0 A is not blocked by benchmark composition.** A perfect selector on current pools gives
1.711 - the pools are good enough. It is blocked by the scoring channel, since the best
conceivable distance-based scorer lands at 1.994 and neither available channel puts its minimum
at the native. The realistic near-term target is **~2.5 A** from selection alone on the 73
tractable targets.

### Correlates of failure, with the confounds removed

Partialled for length AND pool best throughout (the two confounds behind this project's earlier
false leads):

| feature | raw rho | partial rho | 95% CI | perm p |
|---|---|---|---|---|
| dist_mae (the distogram's own distance error) | +0.696 | **+0.739** | [+0.638, +0.817] | 0.0002 |
| ss_sheet | +0.389 | **+0.379** | [+0.190, +0.506] | 0.0002 |
| ss_helix | -0.339 | **-0.334** | [-0.498, -0.162] | 0.0006 |
| ensemble spread | +0.021 | +0.013 | [-0.167, +0.185] | - |
| length | +0.076 | +0.052 | [-0.123, +0.239] | - |

Everything else (glycine, proline, Rg, contact order, positive-phi, charge, hydrophobicity,
sequence entropy, library rank) has |partial| <= 0.16 with CIs spanning zero.

Adding dist_mae as a covariate drops ss_sheet to +0.224 [+0.023, +0.372] and kills ss_helix at
-0.105 [-0.292, +0.067], while dist_mae partialled for sheet, helix, length and pool best is
unmoved at +0.704. **Causal chain: sheet content -> the distogram predicts distances badly ->
the selector picks badly**, with a small residual sheet effect. Beta content also independently
hurts retrieval (ss_sheet vs pool best, partialling length: +0.381 [+0.231, +0.540]) - it costs
twice.

**Two negatives worth keeping.** Ensemble spread does NOT predict the selection gap (+0.013):
we do not fail on the ill-posed targets, they merely have noisy scores. And length does not
either (+0.052), which retires "shorter peptides are easier" on this instrument.

Caveat: ss_sheet is read from native torsions, so it describes what we fail on and is not a
runtime confidence signal. The native-free version (sequence-predicted SS) is untested.

### The 24 impossible targets
All solution NMR, whole-library best 2.06-2.91 A. Coil-rich and helix-poor (coil 0.459 vs 0.320
instrument mean; helix 0.310 vs 0.491), slightly longer (14.1 vs 13.0), and disproportionately
ill-posed (ensemble spread 1.477 vs 1.044). None is catastrophically unreachable - a larger or
non-fragment-based library would move most of them.

## S8-5. Pool construction is capped at 2.406 A at ANY filter skill - and the law that says why

`s8/generate.py`, `s8/test_generate.py` (17 tests), `s8/generate_{diag,transfer,gen,soft,law}.json`.
The full length-n window universe was cached per target (7k-39k windows; the shipped top-2000 is
only ~15% of it) and its BLOSUM top-2000 slice reproduces `s7/debias_cache` to 1.9e-06, so every
arm sits provably on the sprint-7 instrument. Baseline reproduces exactly: 3.454 selected,
1.711 pool best, 4.453 pool mean.

### The hypothesis was right, and it does not help

Sprint 7 finding 6 showed selected RMSD tracks the pool MEAN and moves AGAINST the pool best, so
the natural strategy is to drive the whole distribution down rather than maximise the best
member. A synthetic filter of controlled rank correlation against true CA-RMSD prices that
strategy (DIAGNOSTIC - it reads natives; it prices filters that do not exist):

| rho(key, RMSD) | SELECTED | +-SE | pool best | pool mean | pool SD | frac <2 A | selector pct |
|---|---|---|---|---|---|---|---|
| +0.00 | 3.531 | .151 | 1.783 | 4.810 | 1.314 | .060 | 24.9 |
| +0.30 | 3.410 | .143 | 1.649 | 4.014 | 1.130 | .126 | 35.2 |
| +0.60 | 3.334 | .138 | 1.499 | 3.374 | 0.808 | .204 | 48.5 |
| +0.91 | 3.003 | .123 | 1.412 | 2.850 | 0.494 | .298 | 57.2 |
| **+1.00 (oracle)** | **2.406** | .106 | 1.313 | 2.350 | 0.256 | .412 | **59.9** |

**d(selected)/d(pool mean) = +0.397, r = +0.872**, significant from rho ~ 0.60. The transfer
function is real, monotone and measured.

### But no achievable key can pull it

Measured on the same universe, native-free keys and what they actually select:

| key | rho vs true RMSD | vs BLOSUM | pool best left behind |
|---|---|---|---|
| BLOSUM (incumbent) | **+0.066** | - | 1.71 |
| Ramachandran prior | +0.307 | +0.108 A WORSE | 2.16 |
| rg agreement | +0.321 | WORSE | 2.29 |
| typicality | +0.410 | WORSE | 2.81 |
| the shipped distogram score | +0.587 | +0.367 A WORSE | 3.40 |

**Every key better correlated than BLOSUM selects worse.** The mechanism is the last column: at
matched rho the synthetic filter leaves pool best at 1.65 A while the real keys leave 2.16-3.40.
A real key's error is correlated with the property that makes a candidate good, so truncating
deletes the near-native band along with the tail. Softening rather than truncating (Gumbel-top-k,
3 keys x 5 temperatures) trades the two off exactly as predicted and never wins.

Note BLOSUM's own correlation with candidate quality is **+0.066** - essentially nil. It
survives not by being informative but by being unbiased with respect to nativeness.

### 41 deployable arms, zero winners

Quality-filtered retrieval (Ramachandran / rg / combined x 3 keep-fractions), torsion
recombination (2-, 3-, mixed-crossover, +rg re-filter), torsion contraction (4 alpha),
cluster-centroid contraction (4 alpha), peptide-vs-fragment context matching, whole-universe,
soft sampling. **Arms significantly better than base: 0.** Best is `contract0.5` at **-0.021 A
[-0.157, +0.114]**, 68W/58L, inside the 0.147 A SE. **The 14 arms that improved the pool mean
(by up to 0.420 A) moved selected +0.054 A - the wrong way.**

### Why contraction specifically cannot work

The pool is spread **4.142 A** wide (mean pairwise CA-RMSD) about a medoid **3.665 A** from the
native, while the weak selector already returns 3.454 A. **Contracting to zero spread lands on
the medoid, which is worse than what the selector already achieves.** Shrinking spread is
strictly harmful without re-centring, and re-centring is recognition.

Refinement of sprint 7 finding 6: the returned answer sits at the pool's **25th percentile** on
average (per-target 28.9 +- 24.2), and selected tracks the pool's **q30** (slope +0.267,
r +0.643) better than its mean (slope +0.068, r +0.254).

### The law that closes the axis

Across all 55 arms - deployable, synthetic and oracle, same targets, same selector:

**sel = 0.257 * pool_best + 0.311 * pool_mean + 1.673**,  R^2 0.79, RMSE 0.108 A

**The 1.673 A constant is the selector's incompetence, and no pool improvement touches it.**

The sharpest single number: the ORACLE best-500 subset of this library - pool mean 2.350,
SD 0.256, 41% of candidates already under 2 A - still selects only **2.406 A**, and the
selector's percentile *degrades* from 24.9 to **59.9**. On a pool where most candidates are
good, the shipped score is worse than a coin flip. It is not merely uninformative in the
near-native regime; it is anti-informative.

### Verdict

Driving the pool's distribution down does move selected RMSD, but selection failure defeats it
in two independent ways: no achievable native-free key shifts the distribution without biasedly
destroying the near-native band, and even a perfect filter caps the axis at 2.406 A because the
selector's skill collapses on the tighter pool. **Pool construction cannot deliver sub-2.4 A on
this library at any filter skill.**

Nothing was promoted to the 24-target dev set: the best arm is inside the SE, and spending the
held-out set on noise is how a null becomes a claim.

## S8-6. Retrieval is fixable; fixing it is worth 0.016 A. The 28% "retrieval opportunity" is 0.4%

`s8/retrieve.py` (7 resumable stages), `s8/test_retrieve.py` (21 tests, 0 failures),
`s8/retrieve_screen.json`, `_chain.json`, `_analyse.json`. Same 126-target instrument, same
selector. `blosum` reproduces sprint 7 exactly (**3.454** selected, 1.711 pool best, 1.994
distance-oracle at K=500), and a test asserts the cached BLOSUM ranks equal
`s8/triage_libceil.json`'s target by target.

S8-4 credited 29 `retrieval_limited` targets with **0.960 A, 27.8% of the mean**: a sub-2 A
window exists in the leakage-safe library and BLOSUM buries it at median rank ~2,006. This
study went looking for a better retrieval key. **Two of the three questions came back yes; the
third came back no, and the third is the one that decides the axis.**

### The screen: 22 keys over the whole window universe (median 17k windows, 39k max)

| key | best-window pct mean/median | in top-500 | in top-2000 | rec2@500 | pool best@500 | q30@500 | sub-2 count@500 | pool mean@500 |
|---|---|---|---|---|---|---|---|---|
| **blosum (incumbent)** | 41.9 / **39.3** | 18/126 | 35/126 | 73/126 | **1.711** | 3.632 | 54.5 | 4.453 |
| ident (raw identity count) | 41.0 / 39.1 | 17 | 33 | **79** | **1.663** | 3.676 | 49.2 | 4.504 |
| esm_cos (ESM-2 window cosine) | 41.8 / 39.9 | 16 | 31 | 76 | 1.692 | 3.693 | 49.6 | 4.506 |
| rand (control) | 49.1 / 49.4 | 6 | 20 | 70 | 1.798 | 4.116 | 29.6 | 4.816 |
| **disto (shipped distogram)** | **26.3 / 16.5** | **27** | **57** | 60 | 2.161 | **3.342** | **101.5** | **3.605** |
| dconf (its 30% most confident pairs) | 25.1 / 18.3 | 24 | 55 | 70 | 1.968 | 3.346 | 99.9 | 3.726 |
| dshort (sep <= 5 -- a secondary-structure key) | **24.4** / 17.7 | 23 | 46 | 64 | 2.062 | 3.442 | 96.7 | 3.853 |
| fuse_bs (blosum + dshort, rank sum) | 28.2 / 22.3 | 21 | 44 | 73 | 1.852 | 3.410 | 99.1 | 3.899 |

**Yes, retrieval ranking is fixable.** The shipped distogram used as a RETRIEVAL key moves the
structurally best window from BLOSUM percentile 39.3 to **16.5** (median) and from 35/126 to
**57/126** inside the top 2000. It roughly **doubles the near-native band** (84.8 -> 167.2
candidates within library-best + 1.5 A; sub-2 A count 54.5 -> 101.5), improves the pool's q30
from 3.632 to 3.342 and its mean from 4.453 to 3.605 -- the largest pool-distribution move
measured in this project.

**No, that is not the same as a better pool.** Pool BEST gets worse (1.711 -> 2.161) and the
distance-ORACLE ceiling of the pool -- rank it by the target's own true distances, sprint 7
finding 7's statistic -- gets worse with it (**1.994 -> 2.542**). *Random* retrieval builds a
better pool for a perfect ranker (2.135) than the distogram key does. The distogram key
concentrates the pool on one predicted structural type: many near-native candidates on the
targets it gets right, none at all on the targets it does not. **No key in the screen beats
BLOSUM's oracle ceiling; the best is `ident` at 1.983.**

**ESM is not a retrieval key.** Sprint 7 finding 11 established that ESM beats one-hot by
0.288 A as the INPUT to a learned distance predictor. As a similarity kernel over windows it is
indistinguishable from BLOSUM: percentile 41.8 vs 41.9, pool best 1.692 vs 1.711, selected
3.470 vs 3.454. ESM-as-features and ESM-as-metric are different claims and only the first holds.

### The chain: 24 arms, 126 paired targets, SE 0.147 A on the mean and 0.046 A on a difference

| arm | selected | d vs blosum | 95% CI | W/L | pool best | oracle | rho_g | <2 A | <1.5 A |
|---|---|---|---|---|---|---|---|---|---|
| *oracle_injectall* (DIAGNOSTIC) | *3.384* | *-0.070* | *[-0.151, +0.005]* | *25/11* | *1.313* | *1.532* | *0.461* | *0.270* | *0.183* |
| *oracle_inject1* (DIAGNOSTIC) | *3.438* | *-0.016* | *[-0.037, -0.001]* | *4/0* | *1.313* | *1.522* | *0.567* | *0.214* | *0.167* |
| fuse_be (blosum + ESM) | **3.447** | -0.007 | [-0.100, +0.083] | 34/33 | 1.694 | 1.992 | 0.574 | 0.222 | 0.151 |
| **blosum (incumbent)** | **3.454** | - | - | - | 1.711 | 1.994 | 0.568 | 0.214 | 0.159 |
| ident | 3.457 | +0.003 | [-0.087, +0.095] | 38/34 | 1.663 | 1.983 | 0.566 | 0.214 | 0.143 |
| esm_cos | 3.470 | +0.016 | [-0.088, +0.118] | 45/46 | 1.692 | 1.983 | 0.583 | 0.206 | 0.143 |
| aug_dconf | 3.488 | +0.034 | [-0.070, +0.138] | 54/43 | 1.664 | 1.970 | 0.531 | 0.222 | 0.159 |
| aug_esm_cos | 3.512 | +0.058 | [+0.008, +0.117] | 5/13 | 1.685 | 1.968 | 0.567 | 0.206 | 0.143 |
| fuse_bd2 | 3.532 | +0.078 | [-0.015, +0.177] | 42/45 | 1.890 | 2.171 | 0.177 | 0.198 | 0.135 |
| dconf | 3.543 | +0.089 | [-0.027, +0.202] | 55/57 | 1.968 | 2.249 | 0.239 | 0.214 | 0.135 |
| dshort | 3.547 | +0.093 | [-0.016, +0.203] | 55/60 | 2.062 | 2.358 | 0.270 | 0.198 | 0.127 |
| disto | 3.562 | +0.108 | [+0.005, +0.211] | 54/60 | 2.161 | 2.542 | 0.058 | 0.190 | 0.119 |
| aug_disto | 3.562 | +0.108 | [+0.005, +0.211] | 54/60 | 1.692 | 1.972 | 0.533 | 0.190 | 0.119 |
| dz | 3.565 | +0.111 | [+0.007, +0.213] | 54/61 | 2.148 | 2.493 | 0.051 | 0.190 | 0.119 |
| rand (control) | 3.569 | +0.115 | [+0.007, +0.223] | 57/69 | 1.798 | 2.135 | 0.585 | 0.198 | 0.127 |

blosum median 3.478, SD 1.646, best 0.53, worst 8.13.

**Not one deployable arm beats BLOSUM.** The best, `fuse_be`, is -0.007 A on a 0.046 A SE with
34 wins and 33 losses. BLOSUM does significantly beat random (+0.115 [+0.007, +0.223]),
reconfirming sprint 7 finding 12 on this instrument.

### The decisive control: perfect retrieval, purchased by cheating

A retrieval key can do no better than putting the library's structurally best window into the
pool. So put it there. `oracle_inject1` adds the library argmin to BLOSUM's top-500;
`oracle_injectall` adds every sub-2 A window in the library (up to 250). ORACLE-CONSTRUCTED
POOLS, diagnostic only, never mixed into a system number.

| pool | pool best | oracle ceiling | selected | d vs blosum | 95% CI | W/L |
|---|---|---|---|---|---|---|
| blosum K=500 | 1.711 | 1.994 | 3.454 | - | - | - |
| + the single best window | **1.313** | **1.522** | 3.438 | **-0.016** | [-0.037, -0.001] | **4/0** |
| + every sub-2 A window | **1.313** | **1.532** | 3.384 | **-0.070** | [-0.151, +0.005] | 25/11 |

**Perfect retrieval is worth 0.016 to 0.070 A, and the larger figure is not significant.**
`oracle_injectall` returns the *identical structure* to plain BLOSUM on **90/126** targets: up
to 250 sub-2 A candidates were added to the pool and the selector did not notice them.

**On the 29 retrieval-limited targets themselves** -- the entire premise of the study -- perfect
retrieval takes pool best from 2.351 to 1.566 and the pool's oracle ceiling from 2.844 to 1.860,
and moves selected RMSD by **-0.012 A [-0.035, +0.000], improving 2/29** for one window and
**-0.018 A [-0.045, +0.000], improving 3/29** for all of them. Their contribution to the
aggregate mean falls from 0.960 to 0.955.

**So the 0.960 A credited to retrieval is worth 0.004 A of the aggregate mean when retrieval is
made perfect.** The class is real -- the pool genuinely is missing the answer -- but it is not an
opportunity. `retrieval_limited` targets are selector-limited like every other class.

### Why every screening metric misleads here, including the ones the coordinator asked for

Across the 19 arms with both screen and chain numbers:

| pool statistic | r with selected RMSD | reading |
|---|---|---|
| number of sub-2 A candidates | **+0.582** | MORE near-native candidates, WORSE answer |
| near-native band count | **+0.597** | same |
| q30 | -0.359 | a better q30 goes with a worse answer |
| pool best | +0.714 | the only statistic with the intuitive sign |

Those are across-arm correlations over two key families and should not be over-read. The
unconfounded version is per-target: **across the 13 keys, the within-target rank correlation
between pool best and selected RMSD is +0.046 +- 0.039 (n=126)**, and for pool mean
+0.004 +- 0.039. Which pool you hand this selector has, per target, no bearing on what it
returns.

`aug_disto` is the sharpest single case, and it is exactly the additive-augmentation design that
cannot destroy the band: keep BLOSUM's own top-400, ADD the distogram's top-100. Pool best
**improves** 1.711 -> 1.692, the perfect-selector ceiling **improves** 1.994 -> 1.972, and
selected RMSD gets **0.108 A worse**. `aug_esm_cos` changes the answer on only 18/126 targets
and loses 13 of them.

The mechanism for the replacement keys is visible in `rho_g`: pooling by the distogram collapses
the selector's own global rank correlation from **+0.568 to +0.051**, because retrieving on the
selector's score removes the variance the selector ranks on. But that is not the whole story --
`rand` has the HIGHEST rho_g of any arm (+0.585) and the WORST selected RMSD (3.569). No
property of the pool predicts the answer strongly in either direction.

### The one place retrieval does pay, and its price

Repeat the injection at **K=25**, the interior optimum of sprint 7's K curve:

| arm | selected | d vs blosum | 95% CI | W/L | pool best |
|---|---|---|---|---|---|
| *oracle_injectall* K=25 | ***3.194*** | ***-0.245*** | ***[-0.374, -0.121]*** | *37/17* | *1.313* |
| *oracle_inject1* K=25 | ***3.326*** | ***-0.114*** | ***[-0.196, -0.046]*** | *13/3* | *1.313* |
| blosum K=25 | 3.439 | - | - | - | 2.350 |
| ident K=25 | 3.441 | +0.001 | [-0.115, +0.117] | 48/47 | 2.291 |
| fuse_bd2 K=25 | 3.458 | +0.019 | [-0.124, +0.163] | 61/50 | 2.561 |
| esm_cos K=25 | 3.475 | +0.036 | [-0.113, +0.184] | 52/60 | 2.325 |
| dconf K=25 | 3.533 | +0.094 | [-0.046, +0.237] | 57/64 | 2.648 |
| disto K=25 | 3.562 | +0.123 | [-0.033, +0.277] | 56/65 | 2.832 |

**Perfect retrieval is worth 3.5x more in a 25-candidate pool than in a 500-candidate one**, and
at K=25 both injections are significant. The value of retrieval is inversely proportional to how
many ways the selector has to pick badly. **3.194 A is the best number in this study** and the
only place the axis is alive.

It is also out of reach. At K=25 the best real key on pool best is `ident` at 2.291 against the
oracle's 1.313; the best real key on sub-2 recall@25 is BLOSUM itself (53/126); and no real key
beats BLOSUM on selection at K=25 either (`ident` +0.001).


### Diversity-aware top-K is significantly worse, not better

S8-2 measured that the near-native band concentrates in ~2-3 clusters, which is the standard
argument for filling a top-K to cover conformational space rather than to maximise one
similarity score. Maximal-marginal-relevance over the key's own top-4000, trading key rank
against CA-RMSD to what is already chosen, at K=25 where a 25-slot pool is exactly where
redundancy should bite hardest (`mmr_*` arms, 126 paired targets):

| arm | selected | d vs blosum | 95% CI | W/L | pool best |
|---|---|---|---|---|---|
| blosum K=25 | **3.439** | - | - | - | 2.350 |
| mmr_disto K=25 | 3.562 | +0.123 | [-0.033, +0.277] | 56/65 | 2.558 |
| mmr_blosum K=25 | 3.807 | **+0.367** | **[+0.211, +0.531]** | 39/75 | 2.502 |

**Diversifying BLOSUM's top-25 costs 0.367 A, the largest loss of any arm in this study.** It
makes the pool best worse too (2.350 -> 2.502): spreading the 25 slots over the library trades
away candidates near the mode, and the mode is where the near-native mass is. Multimodality
being real does not make diversity a retrieval objective.

### The full-instrument fold-discipline audit

`s8/retrieve_audit.py`, `s8/retrieve_audit.json`. The per-test check runs on one probe target;
this runs the same check on every library member of every target:

**126 targets, 835,247 target-member pairs. 0 identity violations at the 0.6 threshold, 0
members from the target's own identity fold, 0 members equal to the target. Worst identity
anywhere in the instrument: 0.588.**

### What was NOT built, and why

The brief's third suggestion was a learned retriever -- a model trained on training-fold RMSD
labels to score (target sequence, window) pairs. It was not built, and the reason is the oracle
injection rather than a lack of time. A learned retriever's ceiling is perfect retrieval, and
perfect retrieval is worth **0.016-0.070 A at K=500** and **0.114-0.245 A at K=25**, against a
0.147 A SE on the instrument mean. Building a supervised ranker to chase at most a fraction of
0.245 A, on an axis where the six keys that most improved retrieval ranking all made selection
worse, is not a defensible use of a contended box. If the selector is ever fixed, the K=25
injection numbers say what to re-run first.

One further constraint worth recording for whoever tries: there is **no held-out target set left
to train such a model on**. Every 9-16mer in `peptide_db` whose identity cluster is disjoint
from dev(24), benchmark(60) AND the 126-target instrument numbers **zero** -- the instrument
already takes one representative of every remaining cluster. A learned retriever would have to
be leave-fold-out on the instrument itself, or trained on `fragment_db` entries used as
pseudo-targets.

### Verdict

**Retrieval is fixable as a ranking problem and worthless as a lever.** The distogram key more
than halves the best window's percentile and doubles the near-native band; perfect retrieval is
worth 0.016-0.070 A at the shipped K and 0.114-0.245 A at K=25; and no achievable key is worth
anything at all. The 27.8% of the mean that S8-4 attributed to retrieval is 0.4% of the mean in
practice, and S8-4's class labels should be read as descriptions of where the POOL fails, not as
recoverable fractions of the mean.

This closes the retrieval axis and agrees with S8-5 from the opposite direction: that study
found pool construction capped at 2.406 A because the selector collapses on good pools; this one
finds the selector barely responds to a good pool at all. Note the two are not identical --
S8-5's fitted law `sel = 0.257*pool_best + 0.311*pool_mean + 1.673` predicts 0.102 A from
`oracle_inject1`'s 0.398 A pool-best gain, and the direct intervention measures 0.016 A. A
cross-arm regression over-prices pool best; the injection is the causal test.

**Both remaining classes -- `rankable_missed` (46 targets, 1.338 A) and now `retrieval_limited`
(29, 0.960 A) -- are the same defect, and it is recognition.**

Nothing was promoted to the 24-target dev set or the 60-target benchmark; neither was read.

**Leakage audit.** `window_keys()` takes the target's SEQUENCE, its fold index (used only to
choose the fold-disjoint distogram model), and the library windows' own coordinates and
sequences. Library windows come from fold-disjoint training structures, so their geometry is a
legitimate feature. A test NaN-poisons the target's native CA coordinates and asserts all 13
keys are bit-identical over 39,247 windows; a second asserts the signature cannot even accept a
native; a third does the same for the ESM keys. The library is `s7.audit.build_pool_members` --
out-of-fold peptides plus `distogram._fold_fragments` under the standard identity filter --
with **0 identity violations over all 835,247 target-member pairs on the whole instrument**
(`s8/retrieve_audit.py`; worst identity anywhere 0.588 against the 0.6 threshold), **0 members
from any target's own fold**, **0 members equal to a target**, and **0 dev or benchmark targets
in the cache**. Natives enter only as `rr` labels attached after the keys are computed, in the reported
`oracle` ceiling column, and in the two arms named `oracle_*`.

## Coordinator correction: class contributions are not recoverable fractions

I read S8-4's class decomposition as an opportunity table and commissioned the retrieval study on
that basis, telling the record that the 29 retrieval-limited targets were "28% of the mean and
fixable". **That inference was wrong and S8-6 refutes it directly.**

The decomposition says the 29 retrieval-limited targets contribute 0.960 A (27.8%) of the
aggregate mean. It does NOT say that fixing retrieval recovers 0.960 A. The control that settles
it: injecting the library's best window into the pool **by cheating** is worth **-0.016 A
[-0.037, -0.001]**, and injecting *every* sub-2 A window is worth **-0.070 A** while returning the
**identical structure on 90/126 targets**. On the 29 retrieval-limited targets themselves, pool
best improves 2.351 -> 1.566 and selected moves **-0.018 A**, improving 3 of 29. Their 0.960 A
share is a **0.004 A** opportunity.

**The correct reading, which the retrieval agent supplied:** class contributions describe where
the POOL fails, not fractions of the mean that are available to recover. The selector ignores good
candidates handed to it, so improving what is in the pool does not move the answer. Both
"rankable_missed" and "retrieval_limited" are the same defect - recognition - and the labels
distinguish only how the pool happens to fail alongside it.

This is consistent with, and sharper than, the S8-5 law. That law predicts a pool_best improvement
of 0.785 A on those targets should buy 0.257 * 0.785 ~ 0.20 A; the measured value is 0.018 A. **So
the law itself over-predicts the value of pool improvements when the improvement is injected
rather than earned** - the selector is even less able to exploit an added good candidate than its
average behaviour suggests.

Methodological note for the record: this is the fourth time in two sprints that a quantity which
looked like an opportunity turned out to be a description. The others were the inter-generator
agreement signal (a difficulty proxy), the Amber in-band advantage (a compactness proxy and a
pool-size artifact), and the ESM/one-hot tie (an MAE artifact). **The reliable tell is that none of
them was measured with a control that isolated the causal claim from the descriptive one.** Here
the control was an oracle injection, and it cost one experiment to run.

## S8-7. Inverse folding is the least misspecified objective measured, and it is still misspecified

`s8/invfold.py`, `s8/test_invfold.py` (13 tests), `s8/invfold_{leak,verify,train,select,report,chirality,final_dev,final_dev_report}.json`.
126 tuning targets, K=500, every arm scored on byte-identical candidates in one pass. Pool
reproduction verified against `s7/debias_cache` element-by-element (max |rr diff| 2e-7).

### The decisive experiment

Score `P(target sequence | candidate structure)` with a learned inverse-folding model and ask
where the native falls in its own pool's score distribution.

| objective | argmin | mean pct | median pct |
|---|---|---|---|
| **iv_nbid (inverse folding, + neighbour identity)** | **6/126** | **30.1** | **20.1** |
| iv_geom (geometry only) | 3/126 | 36.4 | 32.2 |
| distogram, matched, measured here | 4/126 | 36.0 | 29.0 |
| distogram, s7 finding 5 | 3/126 | 36.8 | 32.8 |
| all-atom Amber, s7 finding 10 | 0/70 | 54.3 | 61.5 |
| shuffled-sequence control | 0/126 | 45.1 | 45.0 |

`iv_nbid` vs the distogram on identical pools: **-5.90 percentile points [-11.98, +0.17],
p = 0.057, better on 75/126**; below the 10th percentile on 34.1% of targets against 23.0%.

**This is the best-specified objective the project has measured, and on 120 of 126 targets
something still scores better than the truth.** Three independent channels now fail to rank the
native first.

### The control that locates the signal exactly

Against a shuffled permutation of the target's OWN sequence (same model, same candidates, same
composition and length - only position assignment destroyed):

| statistic | true vs shuffled | 95% CI | p |
|---|---|---|---|
| native percentile | **-15.03 pct** | [-21.37, -8.68] | **<1e-4** |
| held-out recovery, native - pool mean | **+4.07 pts** | [+2.26, +5.89] | **<1e-4** |
| **SELECTED CA-RMSD** | **+0.121 A** | [-0.239, +0.482] | **0.51** |

The native is explained by its own sequence significantly better than by a shuffled one. But
scoring the POOL with a shuffled sequence selects just as well as with the true one. **The
sequence signal is real at the native and evaporates among the decoys.**

### Selection (primary), paired, SE 0.151 A

| arm | selected | vs dist | 95% CI | p | W/L | <2.0 A | <1.5 A |
|---|---|---|---|---|---|---|---|
| comb dist+iv_nbid w0.5 | 3.419 | -0.087 | [-0.209, +0.034] | 0.16 | 58/59 | 26.2% | 14.3% |
| **distogram** | **3.507** | - | - | - | - | 21.4% | 13.5% |
| iv_nbid alone | 4.127 | +0.620 | [+0.344, +0.896] | <1e-4 | 45/80 | 19.0% | 11.1% |
| pool best | 1.721 | | | | | | |

Nothing beats the distogram; shortlist arms (ceilings 2.21/2.41 A) are all null. **Final dev-24,
run once:** comb 3.222 vs distogram 3.401, -0.179 [-0.596, +0.238], p=0.38 - and the
native-percentile advantage *inverts* (+6.39, p=0.38). The 126-target lead does not replicate and
is not cited as established.

### Controls

- **Not compactness.** rho(score, rg) is +0.05 for every inverse-folding arm against the
  distogram's +0.240 - genuinely not Amber's failure mode. It does not rescue the result:
  rg-partialled rho is +0.212 vs the distogram's +0.512, worse on 106/126.
- **Data scaling is flat, not harmful.** 7.5x more `prots/` data buys +0.133 A selected (p=0.32)
  and -0.77 pct (p=0.68). Not finding 2's monotone degradation, but the hypothesised transfer
  advantage is **not supported** either.
- **Chirality: sensitive but not correct.** On real pool candidates, reflection changes the CA
  distance matrix by exactly 0.000 A and the distogram score by exactly 0.000, while moving the
  inverse-folding score by 0.321 (0.79 of its own SD), leaving 0.9% rank agreement. But the mirror
  scores worse only **48-52% of the time - chance.** The model has never seen a D-backbone and
  returns noise, not a penalty. Chirality-sensitivity is necessary and not sufficient.
- **Leakage.** Zero PDB-id overlap between `prots/` and the peptide DB / dev-24 / bench-60; all
  787 peptide-DB sequences held out; 966 of 412,530 windows dropped by the `_fold_fragments`
  identity rule.

### The mechanism, and why this closes the axis

**Held-out sequence recovery on the targets' own natives is 11.9-13.0%**, against ~9% for
composition alone and ProteinMPNN's ~52% on real proteins. The agent caught that its
training-split figure (43-54%) was inflated - that split divides residues across ~30 overlapping
windows of ONE chain - and reported the honest held-out number instead.

So inverse folding on isolated 9-16mers carries about **three points of recovery above
composition**. That is the entire channel.

**At 9-16 residues there is almost no local structural context to condition on**: a 12-mer has no
core, no burial and no tertiary contacts, so the model is reduced to reading Ramachandran
propensity. The hypothesis that "structure -> residue propensity is a tighter mapping than
sequence -> structure" is **true in general and nearly vacuous at this length.**

And the informative residue: **the same recovery signal that identifies natives is worth zero
among decoys, because the decoys are real protein fragments and explain the sequence about as
well as the truth does.** That is the clearest statement yet of why recognition fails here - it is
not that our scorers are bad, it is that a near-native real fragment and the native are close to
indistinguishable by any protein-likeness criterion at this length.

## S8-8. Inside the near-native band: the score is worse than a coin flip, and consensus is the only thing that discriminates

`s8/inband.py` (8 stages), `s8/test_inband.py` (26 tests, all passing),
`s8/inband_{skill,anti,two,comb,dev,leak}.json`. 126-target tuning instrument, K=500, two
pools per target over the same cached window universe: the shipped BLOSUM top-500 (`nrm`)
and the ORACLE best-500 of the universe (`tgt`). 38 strictly native-free signals.

Instrument validity is asserted as tests: the shipped score reproduces **3.454 A** on the
normal pool and **2.406 A** on the oracle-tight pool, pool best **1.711 / 1.313**, and the
dev pass reproduces sprint 7's **3.475 A / 1.539 A** exactly. `stage_leak` NaN-poisons `rr`
and `nat_ca` and asserts all 38 signals are bit-identical - worst |diff| **0.000e+00**.

### A methodological correction that changed the answer

The first version of the selection statistic used `np.argmin`, which returns the FIRST index
of a tie. The candidate order is informative in **both** pools - BLOSUM rank in one, true
CA-RMSD in the other - and several signals are massively tied (`lg_steric` is exactly 0.0 for
most candidates; H-bond counts are small integers). `lg_steric` therefore "selected"
**1.386 A** on the tight pool against the shipped score's 2.406 and won **123/1**. It was a
pure read-out of the oracle sort order.

`sel_of()` now averages the true RMSD over the whole tied argmin set - the expectation under
random tie-breaking - so a constant signal correctly scores the pool mean.
`t_sel_of_permutation` and `t_argmin_artefact` pin it. **Every number below postdates that
fix**; the pre-fix table would have reported four spurious sub-1.8 A "winners".

### The anti-correlation is a COLLAPSE TO ZERO, not a negative

The brief's premise was that in-band rho goes negative on tight pools. It does not. Measured
on two independent axes:

| band width (normal pool) | all | 5 A | 3 A | 2 A | **1.5 A** | 1 A | 0.75 A | **0.5 A** |
|---|---|---|---|---|---|---|---|---|
| shipped score, in-band rho | +0.568 | +0.504 | +0.315 | +0.192 | **+0.126** | +0.080 | +0.068 | **-0.020** |
| fraction of targets negative | 10% | 15% | 25% | 29% | 35% | 40% | 43% | **52%** |

| pool (ORACLE top-m) | 25 | 50 | 100 | 250 | **500** | 1000 | 2000 | 8000 | BLOSUM500 |
|---|---|---|---|---|---|---|---|---|---|
| pool mean (A) | 1.672 | 1.799 | 1.941 | 2.157 | **2.350** | 2.569 | 2.845 | 3.896 | 4.453 |
| shipped score rho | +0.008 | -0.028 | -0.028 | -0.012 | **-0.008** | +0.003 | +0.020 | +0.088 | +0.126 |
| **selected minus RANDOM** | +0.027 | **+0.058** | **+0.036** | **+0.067** | **+0.056** | +0.027 | -0.050 | -0.535 | -0.999 |
| 95% CI | [-.015,+.069] | **[+.030,+.085]** | **[+.002,+.070]** | **[+.029,+.106]** | **[+.007,+.105]** | [-.040,+.095] | [-.149,+.050] | [-.762,-.308] | [-1.198,-.801] |

**The score's in-band rho crosses zero at a pool mean of ~2.8 A and at a band width of ~0.6 A,
and then stays within ±0.03 of zero** with ~48% of targets on either side. It is not
anti-informative. It is *exactly* uninformative.

**But it is significantly WORSE than a coin flip.** Against the honest null - the pool's own
mean, i.e. what a random draw returns - the shipped score selects 0.036-0.067 A worse on the
tight pools at m = 50, 100, 250 and 500, losing on **84-94 of 126 targets**, with four
confidence intervals excluding zero. That, and not a negative rho, is what the 59.9 percentile
of S8-5 was measuring: a scorer with zero rank information still concentrates its picks on the
wrong side of a tight distribution.

**So the 1.673 A constant of the S8-5 law is not a defect to be corrected. It is the whole
score, evaluated where every candidate is already good.**

### Legacy, per term, in its new role - the mandated test

Sprint 6 measured Legacy as a whole-pool selector at +1.77 A worse (n=234). Per-term in-band
skill, on identical candidates, 126 targets. `prg` partials out radius of gyration; `stab` is
whether the per-target in-band rho keeps the same sign in all five folds.

| term | NORMAL rho_bd | prg | selBAND | d vs score | stab | TIGHT rho_bd | prg | selBAND | d vs score | stab |
|---|---|---|---|---|---|---|---|---|---|---|
| steric | +0.060 | +0.068 | 2.629 | -0.028 | yes | -0.026 | -0.015 | 2.145 | -0.059 | no |
| contact | +0.071 | +0.067 | 2.734 | +0.077 | no | +0.022 | +0.027 | 2.178 | -0.026 | no |
| hbond_local | +0.053 | +0.069 | 2.642 | -0.015 | no | -0.070 | -0.034 | 2.238 | +0.034 | **yes** |
| hbond_longrange | -0.013 | -0.026 | 2.758 | +0.102 | no | +0.000 | -0.011 | 2.147 | -0.057 | no |
| coop_helix | +0.013 | +0.023 | 2.675 | +0.018 | no | -0.071 | -0.047 | 2.158 | -0.046 | **yes** |
| coop_sheet | +0.001 | -0.013 | 2.620 | -0.037 | no | +0.025 | +0.031 | 2.115 | -0.089 | no |
| solvation | +0.020 | +0.032 | 2.814 | +0.157 | no | +0.012 | +0.016 | 2.175 | -0.029 | no |
| electrostatic | +0.011 | +0.002 | 2.776 | +0.120 | no | +0.008 | +0.003 | 2.186 | -0.018 | no |
| aromatic | +0.021 | +0.051 | 2.692 | +0.036 | no | +0.026 | +0.039 | 2.167 | -0.037 | no |
| **torsion** | **+0.101** | **+0.108** | **2.553** | **-0.103** | **yes** | -0.070 | -0.035 | 2.172 | -0.032 | **yes** |
| compactness | -0.003 | +0.042 | 2.766 | +0.109 | no | -0.054 | -0.015 | 2.196 | -0.008 | no |
| **all (fitted weights)** | +0.058 | +0.118 | 2.655 | -0.001 | yes | **-0.078** | -0.025 | 2.182 | -0.022 | **yes** |
| hb subset | +0.061 | +0.077 | 2.635 | -0.022 | yes | -0.071 | -0.044 | 2.180 | -0.024 | **yes** |
| burial subset | +0.067 | +0.075 | 2.721 | +0.065 | yes | +0.031 | +0.044 | 2.182 | -0.023 | no |
| packing subset | +0.028 | +0.096 | 2.794 | +0.137 | no | -0.018 | +0.035 | 2.219 | +0.015 | no |
| *shipped score* | *+0.126* | *+0.206* | *2.657* | - | no | *-0.008* | *+0.044* | *2.204* | - | no |

**No Legacy term reaches the shipped score's in-band rho on the normal pool**, and no term
beats it on the deployable whole-pool number (`lg_torsion`'s -0.103 A in-band comes with a
whole-pool 4.079 A against 3.454). The best single term is `torsion`, i.e. Ramachandran
plausibility, which is the one term with no coupling in it.

**On tight pools the whole Legacy energy is reliably ANTI-correlated** - `lg_all` at -0.078,
negative on 66-67% of targets and in all five folds, with `hb`, `coop_helix`, `hbond_local`,
`torsion` and `rama` the same. Among near-native candidates, **the more protein-like a
structure looks, the worse it is**: the library's most idealised fragments are canonical
helices and strands, and real natives are locally irregular.

**Two thirds of that anti-correlation is compactness.** Partialling rg takes `lg_all` from
-0.078 to **-0.025**, `lg_hb` from -0.071 to -0.044, `lg_coop_helix` from -0.071 to -0.047.
This is the *third* time in this project that an apparent physics signal has turned out to be
radius of gyration (the Amber correlation, the agreement signal, now this). Inverting Legacy
is therefore not a usable selector: the residue after the confound is 0.02-0.05, and the
inverted selections (2.08-2.13 A) do not reach what consensus achieves without inversion.

**Legacy's verdict in its new role: it is not a discriminator, but as a LEARNED combination it
is not worthless in the near-native regime.** The 15-column Legacy combiner, fitted
leave-fold-out, selects **2.228 A** on the tight pool against the shipped score's 2.406 - and
**4.103 A** on the normal pool against 3.454. Legacy only helps once the pool is already
near-native, which is the regime it was never previously tested in and the one it was
architecturally intended for.

### What DOES discriminate: consensus, and only consensus

Of 38 signals, exactly one has positive in-band skill in the tight regime:

| | normal pool | oracle-tight pool |
|---|---|---|
| **typicality** (agreement with the pool's OWN mean distance profile) | rho -0.083, sel 3.709 | **rho +0.125, prg +0.177, selPOOL 2.151** |
| per-fold in-band rho (tight) | - | **+0.190 / +0.157 / +0.118 / +0.103 / +0.070** |
| vs a random draw (tight, m=500) | - | **-0.199 [-0.251, -0.148]**, 99W/27L |
| shipped score | rho +0.126, sel **3.454** | rho -0.008, sel 2.406 |

Typicality uses **no sequence and no model**. Partialling rg *strengthens* it (+0.125 →
+0.177), so it is not the compactness confound. It beats a random draw by 0.09-0.22 A with CIs
excluding zero at **every** pool tightness from m=25 to m=2000, and its sign is stable in all
five folds. It is the exact complement of the shipped score: the score wins on loose pools and
loses on tight ones; typicality does the reverse, and they cross at a pool mean of ~2.9 A.

### Ensemble disagreement: null on the deployable pool, and a weak consensus effect on the tight one

`s8/inband_ens.json`. Three independently trained distograms (`s7/repr_models`, arms
`onehot` / `phys` / `pca32`, one seed each), which differ by REPRESENTATION rather than by
seed - the favourable direction for the disagreement hypothesis, since
representation-diverse models err more independently than re-seeded copies of one, so a
null here is the stronger result. `ens_sd` is the spread of the members' pool-standardised
scores; `ens_conf` is the L1 against the ensemble-mean predicted matrix restricted to the
pairs the members agree on most.

| signal | NORMAL rho_bd | selPOOL | d vs score | 95% CI | TIGHT rho_bd | selPOOL | d vs score | 95% CI |
|---|---|---|---|---|---|---|---|---|
| SHIPPED score | +0.126 | **3.454** | - | - | -0.008 | 2.406 | - | - |
| ens_mean (ensemble selector) | +0.019 | 3.675 | +0.221 | [-0.002, +0.444] | -0.077 | 2.468 | +0.062 | [+0.003, +0.121] |
| **ens_sd (pure disagreement)** | **-0.109** | 4.250 | **+0.796** | **[+0.535, +1.058]** | +0.012 | 2.313 | **-0.093** | **[-0.169, -0.017]** |
| ens_pen (mean + spread) | -0.023 | 3.604 | +0.150 | [-0.081, +0.381] | -0.074 | 2.465 | +0.059 | [-0.003, +0.121] |
| ens_conf (agreed pairs only) | -0.000 | 3.743 | +0.289 | [+0.079, +0.499] | -0.027 | 2.454 | +0.049 | [-0.010, +0.107] |

**On the deployable pool every ensemble arm is worse than the single shipped model**, and
the pure-disagreement arm is worse by 0.796 A with a CI far from zero. Its in-band rho is
**-0.109, negative on 63% of targets**: inside the near-native band the candidates the models
DISAGREE about are the good ones, which is the exact opposite of the confidence reading.
This is the same collapse sprint 7 finding 1 recorded for inter-generator agreement, now
measured on ensemble members of one architecture.

On the tight pool `ens_sd` does beat the score by 0.093 A [-0.169, -0.017] while carrying an
in-band rho of +0.012 - i.e. it is not ranking, it is mildly avoiding the pool's outliers.
That is a weak consensus effect, and typicality delivers **2.7x more of it** (-0.255) from a
statistic with no model in it at all.

### The deployable consequence - the first arm in three sprints with a CI excluding zero

If the score is a good coarse filter and a bad fine ranker, and consensus is the reverse,
compose them. `stage_two`, fully deployable, no native anywhere: **filter the shipped BLOSUM
top-500 to the score's own top q, then return the consensus medoid of that subpool.**

| arm | selected | ±SE | median | <2 A | <1.5 A | d vs score | 95% CI | W/L |
|---|---|---|---|---|---|---|---|---|
| **medoid75 = consensus_LFO** | **3.282** | 0.162 | 3.095 | 0.278 | 0.198 | **-0.172** | **[-0.316, -0.027]** | **74/47** |
| typic20 | 3.315 | 0.156 | 3.299 | 0.262 | 0.175 | -0.139 | [-0.262, -0.015] | 71/43 |
| typic30 | 3.319 | 0.157 | 3.211 | 0.270 | 0.175 | -0.135 | [-0.263, -0.007] | 68/46 |
| medoid100 | 3.330 | 0.160 | 3.222 | 0.286 | 0.198 | -0.124 | [-0.283, +0.034] | 68/54 |
| **SHIPPED score** | **3.454** | 0.147 | 3.478 | 0.214 | 0.159 | - | - | - |
| typicality on the WHOLE pool (no filter) | 3.709 | 0.160 | 3.640 | 0.246 | 0.135 | +0.255 | [+0.026, +0.483] | 50/76 |

The filter is essential: consensus over the unfiltered pool is 0.255 A **worse**. `medoid75`
was chosen **leave-fold-out over 20 arms** and all five folds picked it independently, so the
headline is the honest number, not the best cell. Robustness: negative in **all five folds**
(-0.037 to -0.395), bootstrap CI [-0.317, -0.030] with P(improvement) = **0.991**, sign test
**p = 0.0177**, and a trimmed mean (dropping the 5 largest movers each way) of **-0.153**.
Fraction below 2.0 A rises 0.214 → **0.278**, below 1.5 A 0.159 → **0.198**, median 3.478 →
**3.095**.

The sharpest line in the table is the subpool mean. Filtered to the score's own top 75, the
subpool's mean is **3.551 A** - so the score's ordering *within its own top 75* is worth
3.551 - 3.454 = **0.10 A**, while consensus over that identical set is worth **0.27 A**.

### The single dev pass

One pre-registered arm, no iteration; a test asserts `DEV_ARM` is what the tuning folds chose.
The 24 dev pools were rebuilt with coordinates through `audit`'s own path and reproduce
`s7/audit_cache` to 1.9e-06.

| arm | selected | ±SE | median | <2 A | d vs score | 95% CI | W/L |
|---|---|---|---|---|---|---|---|
| SHIPPED score | 3.475 | 0.354 | 3.526 | 0.208 | - | - | - |
| **medoid75 (pre-registered)** | **3.316** | 0.373 | 3.183 | 0.250 | **-0.159** | [-0.456, +0.139] | 14/10 |
| typic75 (secondary, NOT pre-registered) | 3.215 | 0.340 | 3.230 | 0.208 | -0.260 | [-0.585, +0.065] | 15/9 |

**This replicates the direction and very nearly the magnitude (-0.159 against tuning's -0.172),
and it does not establish significance, because it cannot.** Dev's SE is 0.354 A against a
0.17 A effect; n=24 has no power to resolve it. Reporting it as confirmation would be wrong.
The tuning result is the evidence; dev is a consistency check that passed.

### Learned combinations buy nothing over the single best signal

Leave-fold-out ridge over the signal block, standardised within each target:

| combiner | NORMAL selPOOL | TIGHT selPOOL |
|---|---|---|
| SHIPPED score | **3.454** | 2.406 |
| score + legacy (5 cols) | 3.462 | 2.247 |
| distogram family (12 cols) | 3.504 | **2.120** |
| everything (38 cols) | 3.511 | 2.127 |
| Legacy only (15 cols) | 4.103 | 2.228 |
| physics/geometry only (11 cols) | 4.148 | 2.318 |

**Nothing beats the shipped score on the deployable pool.** On the tight pool the best combiner
(2.120) improves on typicality alone (2.151) by 0.03 A - the gain is typicality, not learning.
And `all/band` has the best in-band rho of the whole study (+0.235) with a *worse* selected
RMSD than `dist/band` (+0.140), which is the fourth independent instance of sprint 7's rule:
in-band rho does not predict selected RMSD.

### The direct answer

**Does anything discriminate within the near-native band? One thing does, and it is not
physics, not sequence, and not the predictor.** It is agreement with the pool's own consensus.
Every sequence-conditioned and physics-based channel - the shipped distogram, its confidence,
its sd-reweightings, its bagged variant, the disagreement of three independently trained
distograms, all 11 Legacy terms and 4 subsets, H-bond counts, burial, packing, contact
order, backbone strain, Ramachandran likelihood and secondary-structure self-consistency -
is at or below zero in-band once the pool is near-native, and the two that look otherwise
(Legacy, rg) are compactness. Where a signal does beat the score on a tight pool
(ens_sd at -0.093) it does so with an in-band rho of zero, by avoiding outliers: that is
consensus again, weakly, through a much more expensive instrument.

The practical statement: **the shipped objective's whole remaining value is as a coarse filter,
and once it has done that job the best available action is to stop asking it and take the
consensus of what it kept.** That is worth -0.172 A [-0.316, -0.027] on the tuning instrument
and -0.159 A on dev - the first selection gain in this project with an interval excluding zero.
It does not approach the 1.994 A distance-oracle ceiling, and it does not touch S8-5's result
that pool construction caps at 2.406 A. It reduces the 1.673 A constant by about a tenth.

### Two caveats kept on the record

**Geometry.** A candidate is a real CA window carrying its parent's real (phi, psi); the Legacy
terms need N/C/O/CB and are therefore scored on an ideal-geometry rebuild sitting **0.399 A**
CA-RMSD from the window whose RMSD is the outcome. That is the same compromise sprint 7's Amber
comparison made, it is measured rather than assumed, and it is small next to the 1.5 A band -
but it applies to every Legacy and geometry row and to none of the distogram or consensus rows.

**The tight pool is an oracle.** Typicality's 2.151 A is not a deliverable number; it is
measured on a pool nobody can build. The deployable consequence is the two-stage arm, whose
subpool is built by the shipped score alone, and that is where the -0.172 A lives.

## S8-9. The channels' errors ARE decorrelated, and it buys 0.01 A - plus what consensus actually is

`s8/integrate.py` (14 resumable stages), `s8/test_integrate.py` (25 tests, 0 failures),
`s8/integrate_{corr,fuse,filter,critfuse,role,vqe,final,natcons,cvarcheck,k25}.json`,
caches `s8/integrate_{chan,dmat,amber}/`. Same 126-target
tuning instrument, same K=500 BLOSUM pools, rebuilt from `s8/generate_univ`. The shipped
selector reproduces **3.454 A** (pool best 1.7108, pool mean 4.4533) and a test asserts it;
a second test asserts this module's pool is bit-identical to `s8/inband_cache`'s, which is
what makes it legitimate to lift the distogram and Legacy columns from that cache rather
than spend the box recomputing them.

**In one paragraph.** The channels' errors ARE decorrelated (truth-partialled 0.04-0.26) and
that buys +0.004 to +0.011 in rank correlation, because fusion gain goes as the SQUARE of the
weaker channel's skill and every channel is 2-5x weaker than the distogram. Twenty-five
combined objectives and a 165-cell filter grid deliver exactly that: nothing. The multi-
objective score is not a better FILTER either - the plain shipped score wins the grid. The one
quantity with comparable skill is the CONSENSUS CRITERION, which out-ranks the score 3.4:1
inside the score's own top 75; fusing with it gives 3.248 A, -0.206 [-0.364, -0.056], and the
channel's marginal credit is -0.034 [-0.100, +0.026]. Consensus is the whole effect - and the
NATIVE sits at the 82.8th percentile of the criterion consensus minimises, so consensus is
outlier avoidance, not nativeness, which caps every member-selecting operator at the pool's
mode. Of the four mandated components, CVaR earns its place mechanically (+0.113 A by
preventing the state collapse that reproduces the shipped selector exactly), VQE is correct
but unnecessary at 128 enumerable hypotheses, inverse folding is worth +0.034 A with a CI
spanning zero, and Legacy is significantly HARMFUL at -0.106 A [-0.199, -0.027] while being
+0.62 error-correlated with Amber - the two physics channels are largely one channel, so the four mandated channels are really three.

### The precondition: measure the errors before building the combination

Every recognition channel has been tested alone and every one is misspecified. The premise
of a combination is that they fail on DIFFERENT candidates. Five channels on byte-identical
candidates, 126 targets:

| statistic | dist \| legacy | dist \| invfold | dist \| typic | legacy \| invfold | legacy \| typic | invfold \| typic |
|---|---|---|---|---|---|---|
| score rank corr | +0.367 | +0.260 | +0.542 | +0.151 | +0.579 | +0.149 |
| **ERROR corr (truth-partialled)** | **+0.260** | **+0.137** | **+0.445** | **+0.074** | **+0.384** | **+0.038** |
| raw error corr (null +0.5) | +0.470 | +0.366 | +0.570 | +0.413 | +0.540 | +0.353 |

with per-channel skill (Spearman vs true CA-RMSD, whole pool) dist **+0.568**, typicality
+0.401, legacy +0.307, invfold +0.212, and the same after partialling radius of gyration
+0.506 / +0.349 / +0.225 / +0.210. The in-band figures are dist +0.126, invfold +0.090,
legacy +0.058, typicality -0.083 - computed here independently and agreeing with S8-8's
table to three decimals on every channel, which is a stronger check on both harnesses than
either module's own tests.

**A methodological point that changes how the matrix reads.** The brief asked for the
correlation of the channels' *errors*, defined as score rank minus true-RMSD rank. That
statistic has a null of **+0.5, not 0** - both error vectors contain `-rank(true)`, so two
independent zero-skill channels correlate at one half. `t_raw_error_correlation_null_is_half`
pins this at +0.500 ± 0.03 on synthetic data. The statistic that is zero under conditional
independence is the score correlation **partialled for the truth**, and that is the row to
read.

**The answer is yes: the errors are largely decorrelated.** Inverse folding against
typicality is +0.038 - two channels that share neither a training set nor an input, and
whose residual errors are essentially independent. Legacy against inverse folding is +0.074.
The only substantial off-diagonal is distogram-against-typicality (+0.445), and those two
are both functions of the same CA distance matrix, so their agreement is structural.

### And it buys nothing, for a reason the same matrix states exactly

For two z-scored channels, the optimal linear fusion has skill
`sqrt((ra^2 + rb^2 - 2 ra rb rho_ab) / (1 - rho_ab^2))`. Evaluated on the aggregated
correlations above:

| fusion | predicted skill | vs the best member |
|---|---|---|
| dist alone | 0.568 | - |
| dist + legacy | 0.578 | **+0.010** |
| dist + typicality | 0.579 | +0.011 |
| dist + invfold | 0.572 | +0.004 |
| legacy + typicality | 0.412 | +0.011 |
| invfold + typicality | 0.429 | +0.028 |

**Decorrelation is necessary and nowhere near sufficient.** A second channel with half the
first's skill contributes almost nothing however independent its errors are, because the
fusion gain goes as the SQUARE of the weaker channel's skill. The matrix predicts a +0.01
rank-correlation gain before any combiner is fitted, and that is what the combiners then
deliver.

*(Care was needed here: an earlier pass computed this identity per target and averaged the
magnitudes, which predicted +0.41 in-band from channels whose aggregate skills are +0.126
and +0.058. Per-target |rho| on a 130-member band is mostly sampling noise, and averaging
its magnitude manufactures a fusion gain out of nothing. `rho_comb_pred` now takes only
aggregated correlations and a test checks it against a direct least-squares fit.)*

### 25 combined objectives, and the honest ceiling on all of them

Weighted z-rank sums, Borda rank fusion, and leave-fold-out ridge combiners fitted on
training folds only (targets standardised within themselves, so no per-target scale can
leak). Native percentile is the fraction of the pool that outscores the true structure -
the native is appended to the pool as a 501st candidate and scored through the identical
code path, so the two pool-dependent signals (`typic`, `bagged`) see the same distribution.

| arm | selected | d vs shipped | 95% CI | W/L | native pct | argmin at native | rho_g | in-band prg |
|---|---|---|---|---|---|---|---|---|
| lfo_primary3 (LFO ridge, 3 channels) | **3.383** | -0.071 | [-0.214, +0.071] | 61/56 | 39.1 | 2/126 | 0.557 | 0.213 |
| dist + 0.25 invfold | 3.419 | -0.035 | [-0.164, +0.095] | 56/48 | 33.9 | 4/126 | 0.562 | 0.220 |
| Borda over 3 channels | 3.444 | -0.010 | [-0.214, +0.184] | 57/69 | 40.1 | 3/126 | 0.471 | 0.182 |
| dist + 1.0 invfold | 3.448 | -0.006 | [-0.141, +0.127] | 53/64 | **30.3** | 4/126 | 0.479 | 0.176 |
| **SHIPPED distogram** | **3.454** | - | - | - | 36.8 | 3/126 | 0.568 | 0.206 |
| lfo_all (36 columns) | 3.471 | +0.017 | [-0.151, +0.181] | 57/61 | 37.6 | 2/126 | 0.545 | 0.168 |
| dist + 1.0 legacy | 3.495 | +0.041 | [-0.134, +0.217] | 59/66 | 50.6 | 1/126 | 0.484 | 0.181 |
| inverse folding alone | 4.130 | +0.676 | [+0.376, +0.978] | 45/81 | **30.1** | **6/126** | 0.212 | 0.105 |
| Legacy alone | 4.490 | +1.036 | [+0.693, +1.393] | 38/88 | 61.7 | 0/126 | 0.307 | 0.117 |

**Not one arm's CI excludes zero.** The best, a three-channel leave-fold-out ridge, is
-0.071 A on a 0.072 A SE, 61W/56L. And native percentile and selected RMSD point in
*opposite* directions: the arm with the best native percentile (dist + invfold, 30.3, from
inverse folding's own 30.1) is 0.065 A worse in selection than the arm with the worst
(39.1). This is sprint 7's lesson again - a diagnostic is not an outcome.

### The one composition that pays, independently reproduced, and what it actually is

S8-8's `medoid75` reproduces here exactly on this module's own pool and code path:
**3.282 A, -0.172 [-0.320, -0.032], 74W/47L**, asserted as a test. The question this study
could add was the one a multi-channel module is uniquely placed to answer: **is the
multi-objective score a better FILTER?** 15 filter scores x 11 filter sizes = 165 cells,
consensus medoid held fixed:

| filter @ q | selected | d vs shipped | 95% CI | W/L | <2 A | <1.5 A |
|---|---|---|---|---|---|---|
| **shipped score @ 75** | **3.282** | **-0.172** | **[-0.320, -0.032]** | 74/47 | 0.278 | 0.198 |
| dist+0.5 legacy+0.5 invfold @ 25 | 3.294 | -0.160 | [-0.335, +0.008] | 69/52 | 0.302 | 0.198 |
| dist+0.25 legacy @ 25 | 3.297 | -0.157 | [-0.305, -0.012] | 68/49 | 0.294 | 0.198 |
| dist+0.25 invfold @ 75 | 3.299 | -0.155 | [-0.311, -0.010] | 72/51 | 0.286 | 0.183 |
| LFO ridge over 3 channels @ 100 | 3.324 | -0.130 | [-0.294, +0.034] | 68/54 | 0.286 | 0.190 |
| shipped score @ 100 | 3.330 | -0.124 | [-0.286, +0.028] | 68/54 | 0.286 | 0.198 |

**The plain shipped score is the best filter of the 165 cells.** Leave-fold-out selection
over the whole grid returns **3.374 (-0.080)** - worse than the pre-registered `score@75`,
because a 165-cell grid on a 0.147 A SE instrument is mostly noise and LFO correctly
refuses to trust it. Adding channels to the filter does not pay; a better filter is not
what is missing.

### What consensus IS - the diagnostic that caps it

Consensus is not a score, so the sprint's standard native-percentile test has to be
constructed for it: drop the native into the filtered set and rank it by the same criterion
the medoid minimises, mean CA-RMSD to the other members. ORACLE DIAGNOSTIC.

| filter size q | native's percentile | median | what the medoid returns | best in the set |
|---|---|---|---|---|
| 25 | **89.1** | 100.0 | 3.369 | 2.609 |
| **75** | **82.8** | **97.3** | **3.282** | 2.306 |
| 150 | 78.8 | 90.3 | 3.383 | 2.106 |
| 500 | 61.8 | 61.8 | 3.706 | 1.711 |

**The native is farther from the filtered set than 83% of the set's own members, and on
half the targets farther than ALL of them.** Consensus therefore carries no information
about nativeness whatever. It works by outlier avoidance: it returns the candidate
distribution's mode, and the mode happens to be closer to the native than the score's
argmin is, because the score's argmin is an outlier the score likes. That is the whole
mechanism, and it caps every operator that RETURNS A POOL MEMBER: no better choice of member
can move toward a structure sitting on the distribution's periphery. (An operator that
SYNTHESISES a structure can leave the pool and is not bounded this way - see the Verdict.)
It also explains why a better ranker does not compound with it: consensus is not ranking.

### The CVaR gradient, audited against two independent references

The project has a recorded defect (`cvar-gradient-baseline-defect`) on the shipped
`qansatz.cvar_gradient`, so the gradient was verified rather than trusted. The circuit here
is a 7-qubit, 3-layer RY/CNOT-ring ansatz simulated exactly (128 amplitudes), which makes
`p_theta` and every gradient analytic and turns the audit into a real verification rather
than a comparison of two noisy estimates. Three estimators of the same quantity, on the
same `theta` and the same energies, 36 checks over 12 targets x 3 alpha levels:

| estimator | cosine with the exact gradient | ‖g‖ / ‖g_exact‖ |
|---|---|---|
| parameter-shift on each basis probability (exact) | - (reference) | 1 |
| central finite differences on the exact CVaR | **+1.000000** (rel. err **0.0e+00**) | 1.000 |
| sampled score-function, CONSTANT baseline | **+0.994** | 0.997 |
| sampled score-function, **tail-only baseline** (the shipped form) | **+0.524** | **0.534** |

**The defect is diagnosed precisely, and it is a one-line error with a clean statement.** A
baseline leaves a score-function estimator unbiased only if it is CONSTANT in `x`, since
the correction term is `b * E_p[grad log p] = 0`. `qansatz.cvar_gradient` subtracts the tail
mean from the tail entries and leaves the non-tail entries at zero, i.e. it uses
`b(x) = mean * 1[x in tail]` - a FUNCTION of x. That biases the estimator: half the
direction and half the magnitude. `p(x) = <psi|Pi_x|psi>` is the expectation of a projector,
so each basis probability obeys the exact two-term shift rule, which is what makes the
reference exact rather than approximate. `t_score_function_gradient_needs_a_constant_baseline`
pins both numbers.

*(The recorded cosine for the defect was -0.023 on the earlier setup and is +0.524 here.
Different ansatz, different objective; the direction of the finding is the same and the
mechanism is now identified rather than observed.)*

### The criterion out-ranks the score inside the score's own filter

The other half of the same diagnostic, and the one that says where a channel could still
pay. Inside the score's top q, on identical candidates:

| q | consensus criterion, rho | rg-partialled | the shipped score's own rho | native's percentile |
|---|---|---|---|---|
| 25 | +0.175 | +0.163 | +0.048 | 89.1 |
| **75** | **+0.224** | **+0.217** | **+0.065** | 82.8 |
| 150 | +0.228 | +0.220 | +0.131 | 78.8 |
| 500 | +0.425 | +0.382 | +0.568 | 61.8 |

**Inside its own top 75 the shipped score ranks at +0.065 and the consensus criterion ranks
at +0.224 - 3.4x better, and not compactness.** So the mode is a shrinkage estimator:
shrinking a candidate toward the pool's centre buys more variance than it costs in bias
*among library windows*, which is why the criterion ranks them well; but the native is not
a library window and lies outside the region being shrunk to, which is why it ends up at
the 83rd percentile. Both facts are the same fact, and together they are the ceiling on
consensus.

### The one fusion the matrix permits, tested with the discipline it needs

The matrix's rule is that a fusion pays only when the second member's skill is comparable
to the first's, and that is exactly why every score-level fusion was null. The consensus
criterion is the one quantity measured here that meets it (+0.224 against +0.065). Fusing
the channels with the CRITERION rather than with the score, leave-fold-out over a
pre-declared grid:

| grid | leave-fold-out selected | d vs shipped | 95% CI | W/L | per-fold picks |
|---|---|---|---|---|---|
| 21 cells (single-channel additions) | **3.248** | **-0.206** | **[-0.364, -0.056]** | 73/50 | unanimous: `crit@75 + 0.5*iv_nbid` |
| 36 cells (+ the three two-channel additions at each q) | 3.274 | -0.180 | [-0.335, -0.036] | 69/49 | three different cells |
| consensus alone (`crit@75`) | 3.282 | -0.172 | [-0.320, -0.032] | 74/47 | - |

**And the marginal credit to inverse folding is not established.** `crit@75 + 0.5*iv_nbid`
minus `crit@75` on identical candidates is **-0.034 A [-0.100, +0.026], 43W/46L/37 ties**,
negative in all five folds (-0.003 / -0.015 / -0.075 / -0.030 / -0.042), fraction below
2 A 0.278 -> 0.310. The sign is stable and the magnitude is not resolvable: the gain comes
from the size of the wins, not their number. Enlarging the grid to 36 cells moved the
LFO number the wrong way and broke the per-fold agreement, which is what a 0.075 A SE does
to a 36-cell grid and the reason the LFO number rather than the best cell is the result.

Every fold's pick contains `iv_nbid` at weight 0.5. That is the only place in this study
where one of the four mandated components shows a stable-signed contribution.

### K=25, the one regime S8-6 left open - also null

S8-6 measured perfect retrieval as worth 3.5x more in a 25-candidate pool than in a
500-candidate one ("the value of retrieval is inversely proportional to how many ways the
selector has to pick badly"), and no channel had been tested there. Baseline 3.439 on the
BLOSUM top-25; 30 arms:

| arm | selected | d | 95% CI | W/L |
|---|---|---|---|---|
| pool best (ORACLE) | 2.350 | -1.089 | [-1.290, -0.908] | 113/0 |
| consensus over the score's top 8 of 25 | **3.386** | -0.054 | [-0.214, +0.104] | 56/49 |
| consensus over dist+0.25 invfold's top 8 | 3.394 | -0.045 | [-0.195, +0.096] | 56/48 |
| **shipped score argmin** | **3.439** | - | - | - |
| dist+0.25 legacy+0.25 invfold argmin | 3.471 | +0.031 | [-0.059, +0.122] | 33/37 |

The same conclusion in a pool 20x smaller: consensus is worth something and no fused
objective is.

### CVaR-VQE over the discrete hypothesis set

**What the quantum state represents.** Each target's 7-qubit computational basis indexes the
**top 128 candidates of the objective**, so a basis state is a structural hypothesis and the
state is a distribution over hypotheses - not a continuous geometric coordinate. That choice
is forced by this project's own data: sprint 7 measured 1KVG going objective 1.564 -> 1.495
while CA-RMSD went 2.276 -> 3.107, so optimising a misspecified *continuous geometric*
objective harder makes structures worse. S8-2's 2-3 populated basins set the scale, and the
readout is the `p_theta`-weighted consensus medoid - CVaR concentrates the state on the
distribution's good TAIL, consensus returns that tail's CENTRE, and the two are different
operations stated separately.

**Minimising CVaR alone is degenerate, and that is a result, not a bug.** For any alpha the
minimiser concentrates `p` on the lowest-energy basis states, so the readout collapses back
to the argmin. Measured before the fix: every pure-CVaR arm returned the argmin's structure,
state entropy 0.01 bits. The objective that has a non-degenerate optimum is the free energy
`F = CVaR_alpha(E; p_theta) - T * H(p_theta)`, which is the same energy-minus-entropy object
the project's ensemble-selection result rests on. `alpha` says which part of the energy
distribution is scored; `T` says how broad the ensemble is.

126 targets, 7 qubits, 3 layers, exact statevector, exact parameter-shift gradient:

| arm | selected | d vs shipped | 95% CI | W/L | H(p), bits of 7 |
|---|---|---|---|---|---|
| **vqe alpha=0.25, T=0.3** | **3.282** | **-0.172** | **[-0.328, -0.027]** | 61/45 | 6.30 |
| vqe alpha=1.0, T=0.3 | 3.284 | -0.170 | [-0.314, -0.028] | 76/46 | 5.66 |
| *topfrac 0.5* - hard tail, NO VQE | *3.287* | *-0.167* | *[-0.315, -0.027]* | 75/47 | - |
| *boltz T=1.0* - exact Boltzmann, NO VQE | *3.301* | *-0.153* | *[-0.300, -0.011]* | 74/45 | - |
| **vqe, (alpha,T) chosen LEAVE-FOLD-OUT** | **3.314** | **-0.140** | **[-0.288, -0.002]** | 66/48 | - |
| vqe alpha=0.1, T=0.1 | 3.341 | -0.113 | [-0.271, +0.036] | 55/44 | 6.36 |
| *medoid128* - uniform, NO VQE | *3.344* | *-0.110* | *[-0.277, +0.044]* | 67/57 | - |
| *boltz T=0.1* | *3.404* | *-0.050* | *[-0.135, +0.029]* | 37/24 | - |
| **vqe alpha=1.0, T=0.1** | **3.454** | **0.000** | - | 0/0 | **0.08** |
| argmin (the shipped selector) | 3.454 | - | - | - | - |

**The last row is the finding.** With the MEAN objective (`alpha = 1`) and a low temperature
the state collapses onto a single hypothesis - 0.08 bits of 7 - and the system becomes
*exactly* the shipped selector, 3.454 A. The CVaR level is what prevents that collapse: at
the same temperature `alpha = 0.1` gives 3.341, so **CVaR is worth +0.113 A against the mean
objective at matched temperature**, and the state entropy column is the mechanism, measured
rather than asserted. The shipped argmin selector and consensus selection are the two ends of
one free-energy family, and the optimum is strictly interior.

**Handing the VQE a better objective does not change any of that.** The brief's rule is to
establish the objective first and then optimise it, so the whole grid was re-run with the
energy set to the best objective this study measured - the consensus criterion plus half an
inverse-folding term - rather than the shipped score (`s8/integrate_vqe_crit.json`). Best
cell 3.298 against 3.282; the best arm is again the classical hard-tail ablation at 3.287;
the leave-fold-out number is identical at 3.314. The readout *is* consensus, so a
consensus-derived objective is redundant with it, and the family collapses to the same
3.29-3.41 band.

**And the variational optimiser itself is correct but not necessary at this problem size.**
Its best cell (3.282) ties the best classical ablation (`topfrac 0.5`, 3.287) and under
leave-fold-out choice of `(alpha, T)` it is 3.314 - *worse* than that ablation. 128
hypotheses is enumerable, the ansatz constrains `p_theta` to what 21 parameters of RY/CNOT
can express, and the constraint costs about what the optimisation buys. Reporting this as a
quantum advantage would be dressing up a null.

### The integrated system, and what each mandated component is worth

Every stage is deployable; natives enter only as the reported CA-RMSD and in rows named
`oracle`. 126 targets, identical candidates throughout, SE on the instrument mean 0.147 A
and on a paired difference 0.072-0.083 A.

| system | selected | median | d vs shipped | 95% CI | W/L | <2 A | <1.5 A |
|---|---|---|---|---|---|---|---|
| **top 75 by the shipped score, then argmin of (consensus criterion + 0.5 invfold)** - the 21-cell grid's LFO choice, unanimous in all five folds | **3.248** | 3.098 | **-0.206** | **[-0.364, -0.056]** | 73/50 | **0.310** | **0.206** |
| the same, LFO over the 36-cell grid (picks split across folds) | 3.274 | - | -0.180 | [-0.335, -0.036] | 69/49 | - | - |
| top 75 by the shipped score, then consensus medoid (= S8-8's `medoid75`) | 3.282 | 3.095 | -0.172 | [-0.320, -0.032] | 74/47 | 0.278 | 0.198 |
| CVaR-VQE alpha=0.25 T=0.3 over the top 128 | 3.282 | 3.152 | -0.172 | [-0.328, -0.027] | 61/45 | 0.278 | 0.198 |
| the same with (alpha, T) chosen leave-fold-out | 3.314 | 3.152 | -0.140 | [-0.288, -0.002] | 66/48 | 0.286 | 0.206 |
| multi-objective filter + consensus | 3.325 | 3.370 | -0.129 | [-0.292, +0.031] | 69/53 | 0.310 | 0.190 |
| multi-objective argmin, no consensus | 3.411 | 3.267 | -0.043 | [-0.185, +0.094] | 59/59 | 0.246 | 0.183 |
| **SHIPPED distogram argmin** | **3.454** | 3.478 | - | - | - | 0.214 | 0.159 |

Full distributional detail for the three that matter (SD, SE, best, worst - 126 targets):

| system | mean | median | SD | SE | best | worst |
|---|---|---|---|---|---|---|
| crit@75 + 0.5 invfold | 3.248 | 3.098 | 1.802 | 0.161 | 0.222 | 7.911 |
| consensus@75 | 3.282 | 3.095 | 1.816 | 0.162 | 0.222 | 8.761 |
| SHIPPED | 3.454 | 3.478 | 1.646 | 0.147 | 0.302 | 7.917 |

Note the SD *rises* with every consensus arm (1.646 -> 1.80-1.82) even as the mean falls:
consensus moves more targets further, in both directions, and 74/47 of those moves are
improvements. It is a higher-variance selector, not a uniformly better one.

Component credits, each pair differing in exactly one thing, on identical candidates
(positive = the component earns its place):

| component | role it was given | credit | 95% CI | helps/hurts |
|---|---|---|---|---|
| **consensus** | the readout | **+0.086** | [-0.060, +0.237] | 68/54 |
| **inverse folding** | added to the CONSENSUS CRITERION | **+0.034** | [-0.026, +0.100] | 43/46 |
| inverse folding | added to the shipped score | -0.017 | [-0.076, +0.046] | 17/29 |
| typicality | added to the multi-objective score | -0.036 | [-0.089, +0.012] | 17/33 |
| the whole multi-objective fusion | replacing the shipped score as filter | -0.043 | [-0.135, +0.042] | 24/36 |
| **Legacy** | added to the shipped score | **-0.106** | **[-0.199, -0.027]** | 17/32 |

**The verdict on the four mandated components, measured rather than asserted:**

* **VQE** - genuinely present: the state is a distribution over 128 *structural
  hypotheses* per target, optimised by exact parameter-shift gradients, and its best cell
  matches the best classical selector in the family. Genuinely *not* necessary: it ties the
  hard-tail ablation and loses to it under leave-fold-out. At 128 enumerable hypotheses a
  variational optimiser has nothing to buy, and the honest statement is that the
  formulation is right and the problem is too small for it.
* **CVaR** - the only component with an unambiguous mechanical role. At matched temperature
  the mean objective collapses the state to 0.08 bits and reproduces the shipped selector
  exactly (3.454); the CVaR tail level prevents that collapse and is worth **+0.113 A**
  there. Its gradient is verified against two independent references, and the project's
  recorded baseline defect is diagnosed as an x-dependent baseline and priced (cosine
  +0.994 -> +0.524).
* **Legacy** - tested in four roles (term in the objective, whole-pool selector, ensemble
  weight on the consensus, late refiner on the consensus neighbourhood). Every one is
  worse than leaving it out, and as a filter term it is the only *significantly* harmful
  component measured: **-0.106 A [-0.199, -0.027]**. Its errors are also +0.599 correlated
  with Amber's, so it is largely the same channel as the other physics term. **On this
  evidence Legacy does not earn a place in the deployed path, and forcing it in costs
  0.106 A.**
* **AMBER** - measured as the fifth channel of the correlation matrix and as a re-ranker
  and ensemble weight on the BLOSUM top-32, which is where a 7 s/structure evaluation is
  affordable.

### AMBER, where it is affordable

One converged interaction-only ff14SB/GBn2 evaluation is ~7 s, so 126 x 500 is out of reach
on a box shared with three agents. Amber was therefore run on the **first 32 candidates in
BLOSUM order** - a native-free prefix of the shipped retrieval order, and per S8-5 an
essentially unbiased sample of the pool since BLOSUM's own correlation with candidate
quality is +0.066. Every cross-channel number involving Amber is computed on exactly that
subset for all five channels, so the comparison is like for like. `s8/integrate_amber/` is
resumable per target; the numbers below are on the targets finished at write-up time and the
cache continues to fill.

| statistic, BLOSUM top-32 | dist | legacy | invfold | typic | **amber** |
|---|---|---|---|---|---|
| rank skill vs true CA-RMSD | +0.557 | +0.294 | +0.299 | +0.267 | **+0.328** |
| the same, rg partialled out | +0.541 | +0.181 | +0.254 | +0.243 | **+0.193** |
| truth-partialled error corr with Amber | +0.257 | **+0.623** | +0.083 | +0.262 | - |
| score rank corr with Amber | +0.380 | **+0.694** | +0.194 | +0.402 | - |

**Amber and Legacy are largely one channel**: score rank correlation +0.694 and error
correlation +0.623, by far the largest off-diagonal anywhere in the study. Two physics
energies over the same coordinates carry the same information, which is the strongest single
argument against including both - and it means the sprint's "four mandated channels" are
really three. Note also that two thirds of Amber's apparent skill is compactness: +0.328
falls to +0.193 once rg is partialled out, the fourth time in this project that a physics
signal has turned out to be radius of gyration.

Selection on the same 32 candidates (19 targets):

| arm | selected | d vs the score's argmin | 95% CI | W/L |
|---|---|---|---|---|
| pool best (ORACLE) | 2.359 | -0.824 | [-1.193, -0.505] | 16/0 |
| **the score's argmin** | **3.183** | - | - | - |
| score + 0.25 Amber | 3.257 | +0.074 | [-0.100, +0.276] | 3/4 |
| score + 0.5 Amber | 3.266 | +0.084 | [-0.099, +0.287] | 5/6 |
| Amber as an ensemble weight on the consensus | 3.315 | +0.132 | [-0.245, +0.455] | 8/10 |
| **Amber's own argmin** | **4.104** | **+0.921** | **[+0.456, +1.432]** | 4/13 |

Native percentile in the Amber energy distribution **52.7 (median 59.4), argmin at the
native 0/19** - a coin flip, and consistent with sprint 7's 54.3 / 0-of-70 on a different
pool. **Not one of the eleven Amber arms beats the score's own argmin**, and Amber alone is
0.92 A worse.

**And watch the drift, because it is the useful methodological point.** Adding Amber to the
score read -0.091 A at 10 targets, -0.003 at 14, +0.043 at 17 and **+0.084 at 19**. A
10-target subset of this instrument will show you whatever you hoped for. The cache
(`s8/integrate_amber/`, 19 targets x 32 candidates) is resumable and was stopped there to
return an OpenMM process to the other three agents on the box; the conclusion has been
stable in sign since n=14.

Amber's remaining role is refinement of the returned structure. Measured on exactly what
this system selects (`s8/integrate.py refine`, `s8/integrate_refine.json`, 17 targets, the
`consensus@75` pick relaxed under ff14SB/GBn2 to convergence):

| relaxation | d(CA-RMSD) | ±SE | improved |
|---|---|---|---|
| unrestrained (`k=0`) | +0.101 | 0.084 | 8/17 |
| moderate restraint (`k=10`) | **-0.0015** | **0.0063** | 9/17 |

**Nothing.** The restrained relaxation changes the answer by one and a half thousandths of
an angstrom; the unrestrained one is a coin flip that trends harmful. `s8/relax.py` measures
the same operation across a k x steps sweep on pool candidates and is the study that can
settle the general case; on the structures this system actually returns, there is no gain to
credit.

**A confound this measurement had to be rerun to remove, and it is worth recording.** The
relaxation starts from the ideal-geometry REBUILD while the pipeline reports the REAL
WINDOW's RMSD. The first pass differenced against the window and read +0.201 A, charging the
rebuild displacement to Amber. The stage now records both baselines and differences against
the rebuild; `t_refine_baseline_is_the_rebuild_not_the_window` pins the convention. (On the
consensus medoids specifically the gap is only 0.085 A - the medoid is a well-built
structure - so most of that +0.201 was real harm from the unrestrained arm rather than the
confound. Both effects were present and only the rerun separates them.)

### Verdict

**Do decorrelated errors exist to be exploited? Yes, they exist. No, they cannot be
exploited.**

The errors are genuinely and measurably independent - inverse folding against typicality at
+0.038, Legacy against inverse folding at +0.074, on 126 targets with SEs of 0.016-0.020.
The hypothesis was right. What kills it is the second half of the same identity: fusion gain
goes as the square of the weaker channel's skill, every channel is 2-5x weaker than the
distogram, and the predicted gain is +0.004 to +0.011 in rank correlation. Twenty-five
combined objectives then deliver exactly that: nothing, with every CI spanning zero.

The one place the precondition is met is the **consensus criterion**, which out-ranks the
shipped score 3.4:1 inside the score's own top 75, and there a fusion does move the number -
`crit@75 + 0.5*invfold` at 3.248 A, -0.206 [-0.364, -0.056]. But the marginal credit to the
channel is -0.034 [-0.100, +0.026] and enlarging the grid to 36 cells dissolves it, so
the honest reading is that **consensus is the whole effect and the channel is decoration on
it**.

And the consensus effect itself is now bounded from a direction nothing had measured: **the
native sits at the 83rd percentile of the very criterion consensus minimises.** Consensus
is outlier avoidance, not a nativeness signal; it returns the candidate distribution's mode,
and the native is on that distribution's periphery. Every gain reported in this section and
in S8-8 is the mode being a better guess than the score's argmin - which caps *member
selection* at the mode, 3.28 A, while the pool best is 1.71 A.

**One qualification, and it matters.** That percentile bounds operators that RETURN A POOL
MEMBER. An operator that SYNTHESISES a structure is not bounded the same way: it can leave
the pool, and the concurrent projected-consensus work (`s8/consensus2.py`) measures exactly
that at 3.204 A. So the right reading of the 82.8th percentile is not "consensus is finished"
but "consensus cannot be improved by picking a better member, only by leaving the set" - and
the measured value of leaving it is so far ~0.08 A.

**Nothing was promoted.** The dev set and the benchmark were not read. The best arm here is
0.034 A from an arm another agent has already taken to dev, and spending a held-out set on
that difference is how a null becomes a claim.

**Leakage audit.** Every deployable score in this module is a function of the candidate
coordinates, their torsions, the target sequence and the fold's own out-of-fold models.
`test_integrate.py` NaN-poisons `rr` and `nat_ca` in the cached universes and asserts the
inverse-folding columns are bit-identical (worst |diff| 0.000e+00); the distogram, Legacy
and geometry columns are lifted from `s8/inband_cache`, whose own leak stage makes the same
assertion for all 38 signals, and a test asserts this module's pool equals that cache's
target by target. Natives enter as `rr` labels, in the two ORACLE diagnostics
(`natcons`, `pool_best`), and nowhere else. The Amber subsample is a prefix of the BLOSUM
order and reads no native. Combiner weights are fitted on training folds only and every
reported combiner number is leave-fold-out.

**Geometry caveat, carried forward.** Candidates are real length-n CA windows; Legacy,
inverse folding and Amber are evaluated on the ideal-geometry rebuild from the same torsions,
which sits ~0.3-0.4 A from the real window whose RMSD is the outcome. Both endpoints
(`phi[0]`, `psi[n-1]`) are forced to the builder defaults so no torsion-reading channel can
discriminate on provenance rather than conformation.

## S8-10. Generation is worth nothing as candidates - but it exposes a second law with a far lower constant

`s8/diffuse.py`, `s8/test_diffuse.py` (33 tests), `s8/diffuse_{data,eval,chain,mds_mirror}.json`.
The first generative model in this project: a sequence-conditioned autoregressive von Mises
mixture over backbone torsions, plus a DDPM arm and Ramachandran baselines. 126-target
instrument; the `base` arm reproduces the shipped 3.454 A and pool best 1.711 A exactly.

### Structure-prediction accuracy

| arm | best@1 | best@10 | best@100 | medoid | circ-mean |
|---|---|---|---|---|---|
| **ar** (AR von Mises mixture, peptides+fragments 8:1) | 4.504 | 2.860 | **2.146** | 3.846 | 4.058 |
| ar_frag / ar_pep | 4.54 / 4.61 | 2.98 / 3.04 | 2.219 / 2.262 | 3.91 / 3.99 | - |
| **torfit** (torsions refined against the shipped distogram) | **3.669** | - | - | **3.616** | 3.971 |
| mds (distogram -> MDS -> torsion manifold) | 4.159 | - | - | 4.159 | - |
| rama / rama3 (no learning) | 4.97 / 4.88 | 3.48 / 3.42 | 2.614 / 2.634 | 4.36 / 4.46 | - |

**The best deployable single best-guess is 3.51-3.62 A.** Best-of-100 reaches 2.146 A, but
choosing which of the 100 requires the native.

### THE SECOND LAW - the most important result of the sprint

Rank the K=500 retrieval pool by CA-RMSD to a PREDICTED STRUCTURE, instead of by agreement with
a predicted matrix. Paired against the 3.454 baseline:

| reference | ref RMSD | selected | d [95% CI] | W/L |
|---|---|---|---|---|
| **best of 100 AR samples - ORACLE** | 2.146 | **2.361** | **-1.093 [-1.336, -0.856]** | 101/25 |
| best of 100 Rama samples - ORACLE | 2.614 | 2.720 | -0.734 [-0.976, -0.500] | 88/38 |
| cross-channel consensus (`s83x_ar`) | 3.324 | 3.331 | -0.123 [-0.273, +0.024] | 68/51 |
| distogram's favourite AR sample | 3.599 | 3.491 | +0.037 [-0.118, +0.193] | 57/67 |
| torfit | 3.616 | 3.560 | +0.106 [-0.005, +0.223] | 48/58 |
| MDS of the distogram | 4.176 | 4.075 | +0.621 [+0.445, +0.798] | 30/91 |

**selected = 0.915 * ref_rmsd + 0.346**,  r = 0.934, over 3,276 (target, reference) pairs,
including the ref=0 anchor which lands exactly on pool best.

Set this beside S8-5's `selected = 0.257*pool_best + 0.311*pool_mean + 1.673`. **This route
replaces the 1.673 A incompetence constant with 0.346 A**, and makes accuracy purchasable at
~0.9 A of answer per A of predictor:

- break-even against 3.454 needs a predicted structure at **3.40 A** (incumbent: 3.51-3.62)
- ~2.5 A selected needs **~2.35 A**
- sub-2.0 A selected needs **~1.8 A**

**-1.093 A on 101/126 wins is the largest effect measured in this project** - and it is an
oracle. But it converts an intractable discrimination problem into a well-posed regression
problem, which is the first time in eight sprints that the objective has had that character.

### Generation adds nothing as candidates - the decisive negative

The S8-8 consensus arm was reimplemented here independently and reproduced to three decimals:
**3.282, -0.172 [-0.316, -0.032], 74/47**. Applying that same operator to the UNION of the
retrieval pool and each generated ensemble:

`s88` 3.282 -> `+rama3` 3.313 -> `+ar` 3.331 -> `+rama` 3.353 -> `+torfit` 3.358.

**Every generated ensemble makes consensus worse** and pushes the interval back across zero.
As raw candidates: `gen_ar` 3.599 (+0.145), `mix_ar` 3.508 (+0.054). Note this also means the
S8-8 gain is specific to the BLOSUM top-500 as constructed, not a general property of
filter-then-consensus.

### Coverage, novelty, validity

- **Coverage**: the AR ensemble holds a sub-2 A structure on 54/126 (sub-1.5 on 40) against the
  K=500 pool's 73/126 and the whole library's 102/126. **Generation does not reach retrieval's
  coverage.**
- **Novelty**: an AR sample sits a median 1.883 A from the nearest pool window and 1.757 A from
  the nearest window of its own training corpus; only 2.7/100 within 0.5 A. It is generating,
  not memorising.
- **Chirality - a concrete harm, measured.** Native positive-phi rate 5.58%; `ar` 6.49%,
  `ar_pep` 6.14%, `mds` 6.72%, but `rama` 2.50% and `ar_frag` 2.67% (the fragment corpus is more
  regular), and **`torfit` 12.29%**. **Refining torsions against a distance objective doubles
  the left-handed rate**, because a distance matrix is exactly mirror-blind - a test shows a
  peptide's own distances reconstruct at 0.008 A with the right mirror and 3.125 A with the
  wrong one, and the Ramachandran mirror tiebreak on MDS is a coin flip (64/126). This is
  sprint 7's chirality mechanism finally caught doing measurable damage.
- Geometry: CA-CA 3.804 +- 0.000 by construction; clash rate 0.8% raw, 0.02% after refinement.

### Data ablation - sprint 7 finding 2, reproduced in a different model class

| corpus | n_train | val NLL (peptides) | best@100 | oracle S8-3 |
|---|---|---|---|---|
| peptides + fragments 8:1 | 22,554 | **0.845** | **2.146** | **2.361** |
| fragments only | 22,000 | 1.138 | 2.219 | 2.434 |
| peptides only | 554 | 0.879 | 2.262 | 2.497 |

**35x more fragment data is 0.26 nats WORSE than 554 peptides** - the distribution-shift result
of sprint 7 finding 2, now reproduced in a completely different model class and loss function.
Fragments-only has the worst likelihood and the second-best coverage: density-fitting and
basin-covering are different objectives.

### Leakage
5 folds, worst identity to any held-out peptide 0.583-0.591 against the 0.6 threshold, 0
violations, 0 held-out peptides in any corpus. A test NaN-poisons all 787 natives and asserts
sampled torsions are bit-identical; the feature signature cannot accept a structure. Dev(24)
and benchmark(60) never read.

### Verdict
**Generation does not add anything retrieval cannot do, at the current selector.** It produces
genuinely novel conformations, covers less than retrieval, and its candidates degrade the one
validated arm. The gap between 2.361 A (oracle pick among its own 100 samples) and 3.49 A (best
native-free pick) is **1.13 A of pure recognition failure** - the same defect as S8-4/5/6, now
isolated inside a 100-member ensemble the model generated itself rather than a 500-member pool
retrieval built.

Not finished on the contended box: the DDPM arm trained on one fold only (9 targets - not a
comparison and not reported as one), and the ESM-vs-one-hot conditioning ablation never ran.
The model-class choice therefore rests on argument plus the Ramachandran floor, not on a
measured AR-vs-diffusion number.

## S8-11. Synthesis beats selection: 3.204 A, replicated on dev

`s8/consensus2.py` (8 resumable stages), `s8/test_consensus2.py` (32 tests),
`s8/consensus2_FINDINGS.md`. Instrument validity asserted before anything ran: shipped score
3.4540, pool best 1.7108, S8-8's `medoid75` 3.2822.

**The result.** Coordinate-average the score-filtered top 75, then PROJECT that average back
onto the manifold of ideal-geometry chains (exact 3.804 A CA-CA).

| | selected | d vs base | 95% CI | W/L | frac <2 A | median |
|---|---|---|---|---|---|---|
| tuning, 126 targets | **3.204** | **-0.250** | [-0.383, -0.117] | 80/46 | 0.278 | 2.966 |
| **dev-24, one pre-registered pass** | **3.221** | **-0.255** | **[-0.501, -0.008]** | 16/8 | | |

Tuning and dev agree to 0.005 A and the dev interval excludes zero. 1.45x S8-8's gain.

**The winning arm abandons selection.** It returns a structure that is in no pool. And the
geometry accounting is the reason it works: the raw coordinate average scores 3.048 but has a
mean CA-CA bond of **2.961 A** against a native 3.812 - it is not a peptide, and **37% of its
apparent gain is geometry rather than conformation.** Rescaling (the naive fix) is 1.04 A
WORSE; projection is the correct inverse. The projection needs no torsions (4 generic starts
reproduce the torsion-seeded form to 0.004 A), which is what made a clean dev pass possible.

**The filter's real job is collapsing the pool onto one mode, not keeping good candidates.**
A RANDOM 75 has a *better* best member (2.096 vs 2.306) and a much worse consensus (3.400 vs
3.048). Filtering 500 -> 75 costs 0.595 A of ceiling and buys 0.424 A of consensus.

**The bias/variance decomposition that explains the whole architecture.** Averageable
independent error is sqrt(3.748^2 - 3.462^2) = **1.437 A**, against **85.3% of the mean square
being bias common to every candidate.** Averaging removes independent error and cannot touch
shared bias. That single figure predicts the flatness of `avg` in m (3.048-3.146 over m=10-250),
the nullity of trimming and weighting, and why synthesis travels 1.145 A from the mode to gain
only 0.234 A.

**Negatives from the same sweep** (7 filters x 8 sizes x 20 operators = 1120 arms):
- **Torsion-space circular mean is +0.618 [+0.265, +0.972] WORSE** - torsion error compounds
  along the chain and exact geometry does not compensate.
- Typicality as a *filter* is a disaster (3.91); score+typicality fusion loses; S8-3's
  RMSD-to-a-reference form ties the score.
- Weighted, trimmed, iterated and cluster-restricted medoids are all null (best 0.023 A).
- Multi-hypothesis: best-of-3 is 2.856 (medoid) / 2.679 (avg) - 0.43 A of ORACLE headroom - but
  all three native-free rules are worse than not splitting, every CI excluding zero the wrong way.
- Refinement does not stack: +0.062/+0.066 A at 5 deg, +0.20/+0.18 at 15 deg, with energy
  falling 6-10 kcal every time.
- LFO over selection operators only gives 3.294: **S8-8 had already found that family's ceiling.**

**No consensus variant carries nativeness information.** The native sits at percentile 81.1
mean / 96.7 median under the medoid criterion, and six variants span 79.9-80.9 - a 1.2 pp spread
on an SE of 2.5. Inside the score's own top 75 the shipped score itself collapses to 79.0
(median 100.0).

**S8-8's rebuild-displacement caveat, corrected:** 0.399 A of displacement is only **+0.009 A**
of error across the pool and **+0.032** on the selected medoid. Real as a displacement, nearly
free as an RMSD penalty.

**Bug found and fixed:** cells were keyed `f"{filter}{m}"`, so `scty`+250 and `scty2`+50 both
produced `"scty250"` and one silently overwrote the other. Fixed with a separator and pinned by
`t_cell_keys_unique`. Neither was a winner, but the first table reported one under the other's
numbers.

## S8-12. AMBER is a validity stage, not an accuracy stage

`s8/relax.py`, `s8/test_relax.py` (12 tests), `s8/relax_findings.md`. 792 relaxations,
14 targets x 12 candidates, on the same 70-target subset and identical pool recipe as
`s7/amber_native.py` - including the plain (not stable) argsort, because BLOSUM ties are large
and a different tie-break gives a different pool. Pool best/mean reproduce sprint 7 to six
decimals and the native interaction energy reproduces bit-for-bit. Genuine ff14SB + GBn2
(`CustomGBForce` verified present).

| k | steps | CA moved | mean dRMSD | 95% CI | W/L | frac improved |
|---|---|---|---|---|---|---|
| 100 | conv | 0.08 | +0.011 | [+0.006, +0.017] | 0/6 | 0.153 |
| 10 (the s7 setting) | conv | 0.22 | +0.012 | [+0.003, +0.022] | 4/10 | 0.351 |
| 1 | conv | 0.51 | +0.026 | [+0.015, +0.037] | 0/6 | 0.431 |
| 0 | 50 | 0.49 | +0.087 | [+0.061, +0.112] | 0/6 | 0.236 |
| 0 | 200 | 0.65 | +0.020 | [-0.045, +0.084] | 4/10 | 0.446 |
| 0 | 1000 | 1.05 | +0.052 | [-0.023, +0.128] | 1/5 | 0.444 |
| 0 | **conv** | 1.39 | **-0.118** | **[-0.221, -0.014]** | **11/3** | 0.595 |

**Converged free minimisation lowers the distribution, but it is COMPRESSION, not refinement.**
By tercile of starting quality: near-native **+0.056** (worse), middle -0.045, far **-0.364**
(73% improved). Within-target rg spread falls 17% while |rg - rg_native| barely moves, and
rho(d_rg, d_RMSD) = +0.217.

**Pool best does not move: +0.013 [-0.197, +0.224], 7W/7L.** And the ceiling check settles it:
against the pipeline's real K=500 pool (best 1.636 A), **0 of 168 relaxed structures is better,
on 0 of 14 targets.** Sharpest case: 1KZ2 candidate 3 starts at 1.121 A - better than the entire
K=500 pool - and converged free minimisation **destroys it to 1.820 A**.

**One mechanism, two symptoms.** A force field whose minima are placed by generic compaction
rather than sequence-specific structure will both fail to refine AND rank by compactness. Sprint
7's ranking failure and sprint 8's refinement failure are the same defect measured twice.

**What AMBER does earn: stereochemical repair.** Restrained relaxation at k=10-100 removes
~10^4 kcal/mol of builder strain for +0.011-0.012 A of CA-RMSD - essentially free. Every
all-atom structure the pipeline emits should pass through it so nothing downstream inherits
builder strain. This matters more now that the winning arm SYNTHESISES a structure by
projection rather than returning a library member.

**Under-powered and declared as such:** the free-energy stage (strain, energy at own minimum,
quasi-harmonic F from a 4 ps Langevin ensemble, Boltzmann-weighted -kT log<e^-E/kT>, basin
width) is built, tested and resumable but completed **1 of 24 targets** after the box sat at
85-95% RAM. **No free-energy claim is made from one target.** The only recordable note: on 1A13
the TOTAL energy puts the native at the 74th percentile against the INTERACTION-ONLY energy's
28th, which independently supports sprint 7's choice to drop bond/angle/torsion terms.

## S8-13. CORRECTION to S8-10: the "transfer law" is the identity map where we operate

`s8/predictor.py`, `s8/predictor_report.log`, 19/19 tests.

S8-10 recorded `selected = 0.915 * ref_rmsd + 0.346` (r=0.934) and I framed it as the sprint's
most important structural result, on the reading that it "replaces the 1.673 A incompetence
constant with 0.346" and makes accuracy purchasable at ~0.9 A per A of predictor. **That
reading is wrong and is withdrawn.**

Re-measured over 2,520 (target, reference) pairs in the band where the pipeline's predictors
actually live:

  **selected = 1.009 * ref + 0.011**  (r = 0.981), and over converged references
  **selected = 1.019 * ref - 0.001**  (r = 0.984)

**The S8-3 arm is the identity map.** The 0.915/0.346 fit spanned 0-8 A and was anchored at the
ref=0 point (where it necessarily returns pool best); inside 3-4 A the slope is 1. Arithmetic
confirms it: `0.940R + 0.248 < R` only for `R > 4.13 A`, and nothing in the pipeline predicts
that badly. Feeding the best available structure through pool selection **costs +0.049
[-0.009, +0.103]** - it hands back a worse pool member than the structure you already had.

**What survives, and it is cleaner than what it replaces: the synthesised structure IS the
answer. Passing it through a selection step adds nothing.** The oracle datum of S8-10 (best of
100 samples at 2.146 A) stands as measured; what does not stand is the extrapolation that a
3.40 A predictor breaks even or that a 1.8 A predictor yields sub-2.0 A.

### The best deployable predictor, and its ceiling

`X_fit` - project onto ideal geometry the coordinate average of (distogram top 15% of the K=500
pool) union (AR model top 15% of 100 samples):

| | selected | d | 95% CI | W/L | <2.0 A | <1.5 A |
|---|---|---|---|---|---|---|
| **X_fit** | **3.217** | −0.237 vs 3.454 | [−0.378, −0.098] | 79/47 | 36/126 | 25/126 |
| | | −0.065 vs S8-8's 3.282 | [−0.123, −0.006] | 74/52 | | |
| P_fit (pool only) | 3.240 | | | | | |

**The cross-channel term is worth 0.023 A and is not significant** - the gain is consensus2's
construction, not the generative channel.

### The decisive ceiling experiment

Same operator, three input sets:

| input set | selected |
|---|---|
| no filter | 3.653 |
| shipped distogram filter | 3.285 |
| **the 75 pool members truly closest to the native (ORACLE)** | **2.220** |

The score is worth **−0.368 [−0.537, −0.204]**; a perfect filter a further **−1.066 [−1.258,
−0.885]**. **Given the right 75, coordinate-average-and-project lands within 0.51 A of pool
best (1.711).**

**So the predictor is not construction-limited.** The whole construction family sits at
3.20-3.34 (consensus2 `fit` 3.204, `P_fit` 3.240, `X_fit` 3.217), and the remaining 1.07 A is
entirely recognition - the axis already closed nine times.

### The ESM ablation diffuse never ran

Trained end to end, 5 folds. **One-hot WINS the held-out peptide NLL on all five folds**
(0.8189 vs 0.8454; ESM is +0.0266 nats worse). **ESM wins every structural estimate**: medoid
−0.161 [−0.326, −0.013], circular mean −0.578, best-of-1 −0.570. **Both effects vanish under
the distogram filter: dgfit −0.025 [−0.127, +0.076].** Real on the raw distribution, gone by the
deployable estimate - and likelihood and structure disagree about which representation is
better, which is a caution against selecting a representation on NLL.

### A clean new negative: a generative model's own likelihood is not a selector

Over its own samples it is worse than the medoid on all five arms (ar +0.096, ar_noesm +0.233,
ar_pep +0.396, ar_frag +0.212, rama +0.671). It is only weakly informative (r = −0.25 with true
RMSD, against the distogram's +0.58), and a weak per-item score loses to consensus. As a filter
over *retrieved* structures it is actively harmful: `Pll_fit` 4.018, **+0.733 [+0.465, +1.006]**.
It reads local torsion typicality and is blind to the fold.

### Chirality: projection does not inherit the mirror defect

X_fit positive-phi **0.0551** against the native instrument's **0.0558**; clash rate 0.004; Rg
6.58 vs native 6.60; CA-CA step exact to 1e-15. It does **not** reproduce `torfit`'s mirror
drift (0.1229), and the reason is precise: **projection minimises CA-RMSD to a chiral object,
not to a mirror-blind distance matrix.**

Note this appears to CONTRADICT the audit's 12.1% positive-phi on consensus2's `fit`. The
implementations differ in convergence - this projection stops 0.92 A from the average, the
audited one 0.79 A - so the over-helical bias is plausibly a function of how tightly the
projection fits an object that is itself compacted (mean CA-CA 2.961 A). Under test.

### Discipline notes
The agent **declined to spend a dev pass**: its best arm is −0.065 over S8-8 with the pool-only
construction accounting for all but 0.023 A, and n=24 at SE 0.354 cannot resolve that. Spending
a pre-registered pass on an effect the instrument cannot distinguish from its incumbent produces
a number that means nothing.

Method note worth keeping: `diffuse.mem_pct` shells out to `Get-CimInstance`, which costs
**minutes per call** on a box at 93% RAM; replaced with a `GlobalMemoryStatusEx` read (0.000 s)
patched in at runtime without modifying the other agent's file.

Leakage: fold-disjoint corpora, worst identity to a held-out peptide 0.583-0.591, 0 violations;
a NaN-poisoning test shows 19 deployable references bit-identical with only the `O_`-prefixed
oracle arm moving.

## S8-14. The training target IS misspecified - by 0.5 A per pair - and fixing it is worth 0.044 A

`s8/enstarget.py` (8 resumable stages), `s8/test_enstarget.py` (44 checks), `s8/enstarget_*.json`.
Instrument reproduces 3.454 / 1.711 / medoid75 3.282 as tests; leakage 160 checks at 0.000e+00.

Every model in this project is trained to predict properties of **PDB model 1**. Since 115/126
instrument targets are solution NMR and 111/126 are multi-model depositions, that means fitting
one arbitrary member of a ~1 A-wide conformational ensemble as if it were exact.

**The defect is real and was never measured before.** 90.3% of training PAIRS come from a
multi-model deposition (698/787 entries, median 20 models). Mean |d(model 1) - d(ensemble mean)|
= **0.497 A**; mean per-pair SD across models **0.635 A**, rising to **1.028 A** in the
|i-j| >= 14 shell. Noise-to-signal is a near-constant 0.135-0.172 per shell.

### The ceiling was computable before training a single model

Label noise enters the error in quadrature. Noise SD **0.635** against model error SD **2.849**
gives a ratio of **0.223**, so a *perfect* de-noising caps at **0.052 A of distance MAE**, and
the ensemble mean of ~20 models should recover ~**0.040** of it.

**Measured: 0.032 A.** Cross-target OLS with length and pool best held fixed (slope reproduced
on three predictors at 0.855 / 0.872 / 0.970) prices 1 A of MAE at ~0.86 A of selected RMSD.

**So the entire axis was worth 0.044 A of selected CA-RMSD against an instrument paired-SE of
0.065** - below the noise floor before any model was trained.

### Arms, all null

| arm | d vs matched model-1 baseline | 95% CI |
|---|---|---|
| pep:ensmean | -0.045 | [-0.172, +0.082] |
| frag:ensmean | -0.061 | [-0.179, +0.058] |
| ensdist | +0.040 | |
| confw (per-pair confidence weighting) | +0.019 | |

Native percentile moves 36.6 -> 37.7-38.1, i.e. slightly **worse**. In-band rho falls +0.099 ->
+0.087. As the FILTER inside S8-8's consensus architecture: 3.4526 vs 3.4533 - worth **0.0007 A**
against an incumbent of 3.282.

Single pre-registered dev pass (`pep:ensmean`, LFO 4/5 folds, rule fixed in code before running):
**+0.070 [-0.118, +0.259]** - does not replicate the tuning sign.

### Four controls kill the mechanism

1. **The differential** (MAE|m1 - MAE|ens, paired against m1 - the statistic that isolates
   "learned a different target"): **+0.0020 [-0.0082, +0.0122]**. The arm sits no closer to the
   ensemble mean than the baseline does.
2. **Noise INJECTION** - add one more draw of the same measured noise instead of removing it.
   It moves the differential by an identical **+0.0022**, improves MAE, and posts the corpus's
   **best** selected RMSD (-0.063). **Adding noise works as well as removing it**, so whatever
   gain exists is generic target smoothing, not de-noising.
3. **Uncertainty calibration** - the one thing only an ensemble can teach. rho(predicted spread,
   target ensemble SD) goes 0.536 -> 0.547, and the **corrupted-label control also reaches
   0.547**.
4. **Dilution.** At 13.4% ensemble coverage against 90.3%, the selection effect is the *same
   size* and the distance effect *reverses sign*. A real mechanism scales with exposure.

### Verdict

**Misspecified? Yes, by 0.5 A per pair. Does fixing it matter? No - and by arithmetic rather
than luck.** Label noise is a second-order term at a noise-to-error ratio of 0.223.

This is the **fifth independent instance of "a better distance matrix is not a better ranking"**,
and the strongest, because this time the better matrix came from repairing a genuine supervision
defect rather than from smoothing or regularisation.

**Do not re-open this with better features or a bigger corpus: better features lower the model
error, which SHRINKS the available gain.**

### Side result

With sequence-only features the shipped `frag` corpus is **worse** than peptides alone
(3.784 / MAE 2.548 against 3.601 / 2.273). **The fragment contribution is feature-dependent** -
nothing in the project had measured that, and it refines sprint 7 finding 2 and S8-10's data
ablation, which both treated the fragment corpus as uniformly harmful or uniformly useful
depending on the model.

<!-- ==================================================================== -->
<!-- Sprint 9: verbatim from s9/FINDINGS.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s9/FINDINGS.md`. Not edited.*

# Sprint 9 findings - breaking the 2.0 A barrier

Objective: mean held-out CA-RMSD below 2.0 A. Entering state: best validated **3.204 A**
(126-target tuning instrument) / 3.221 A (dev-24) against a 3.454 A baseline, from an
architecture that synthesises rather than selects - filter the K=500 BLOSUM pool to the
distogram's top 75, coordinate-average, project onto the manifold of ideal-geometry chains.

## S9-1. The shared bias is real, is one interpretable mode, and is not an error

`s9/bias.py` (9 stages), `s9/test_bias.py` (40 checks), `s9/bias_FINDINGS.md`,
`s9/bias_*.json`, `s9/bias_fitck/` (per-target records behind every table).

**The hypothesis.** S8-11 measured that **85.3% of the candidates' mean square error is bias
common to every candidate**, leaving only 1.437 A of averageable independent error. Common bias
is systematic and systematic error is learnable, so: regress the displacement field from the
synthesised structure to the native on training folds, and apply the correction at inference.
This is a REGRESSION problem, which is the one class that had not repeatedly failed here.

**The conceptual error, caught by the agent.** *"Common to every candidate" is a WITHIN-target
statement.* It does not imply that different targets miss in the same direction. Only the
cross-target component is learnable, and it had never been measured.

| quantity | value |
|---|---|
| cross-target mean pairwise cosine | **+0.1176 +- 0.0027** |
| sign-flip null | +0.0002 +- 0.0029 |
| share of the field common ACROSS targets | **11.7%** (9.3% held out) |
| of which a single constant vector | 5.5 points |
| value of harvesting all of it (ignoring geometry) | **-0.048 A** [-0.096, -0.000] |
| what would matter: 5 per-target deformation numbers | -0.655 A |
| held-out R^2 for those, from 18 target-level features | **0.003-0.018** |

Real at 40 SE, and almost entirely useless: the shared part is small, and the part that would
matter is not predictable from anything available.

**What the shared part actually is: the architecture looking at itself.** The mean local
displacement is **+0.863 +- 0.066 A along the curvature normal** - 13.1 SE, same sign on 91.3%
of targets. The synthesis is **over-curved: CA-CA-CA pseudo-angle 93.21 degrees against the
natives' 103.99**, independently reproducing `s8/audit8`'s 93.5/104.0.

**The mechanism, from the m-sweep.** Individual candidates are native-like (101.3 degrees,
3.81 A bonds). **Averaging contracts the chain** (mean CA-CA 3.810 -> 2.961 A at m=75) and **the
projection buys the lost contour length back as curvature.** The manifold Frechet mean
over-curves identically (98.33 / 95.51 / 93.61 at m = 5 / 25 / 75), so it is **not** the
average-then-project *order* - it is shrinkage in shape space, the same mechanism S8-12 measured
in AMBER's compression.

**Every deployable arm is null, and the geometry control is why.**

| arm | selected | pre-proj | pre-proj CA-CA | d vs 3.204 | 95% CI | W/L |
|---|---|---|---|---|---|---|
| `ridge_auto@1` (deployable) | 3.2055 | 3.0702 | **2.651** | +0.0015 | [-0.0212, +0.0241] | 61/65 |
| `ridge@0.5` | 3.2010 | 3.0854 | **3.134** | -0.0031 | [-0.0128, +0.0066] | 63/63 |
| `dgforce@1` | 3.2060 | 3.2064 | 3.781 | +0.0019 | [-0.0043, +0.0081] | 57/69 |
| `const@1` (control) | 3.2410 | 3.1556 | **2.697** | +0.0369 | [+0.0188, +0.0551] | 45/81 |
| `ridge@0.5~rs` (geometry restored) | 3.5056 | | 3.804 | **+0.3015** | [+0.2088, +0.3943] | 29/97 |

Pre-projection gains of 0.05-0.12 A sit on chains with **2.65-3.13 A bonds** - the S8-11
contraction trap. Restoring legal geometry costs +0.30 to +0.95 A.

**And correcting the bias LEGALLY is monotonically worse.** Driving the scale-free bond angle
toward the natives' own value costs **+0.019 / +0.080 / +0.258 / +0.334 A**. Making the
structure more native-like in that coordinate makes RMSD worse. **The angle is a symptom, not a
lever.**

**Controls, all landed.** Shuffled target (apply A's correction to B): only **+0.0094**
[-0.0010, +0.0197] worse than the real correction - the correction carries almost no
target-specific information. Null features +0.0202. Dropping the ten largest gains turns
**every** arm positive. Leakage 0.000e+00 on 15 held-out targets plus four in-suite poison
tests. No dev-24 pass spent (best effect 0.003 A against a dev SE of ~0.35); benchmark
untouched.

**Verdict.** The shared bias is real at 13 SE and fully explained: it is the deterministic
signature of consensus-then-project, not an error the pipeline is making about the target.
Nine-tenths of the residual is target-specific, and no signal in this project carries it - not
the distogram (0.002 A as a force), not ESM (+0.0009 *worse* than one-hot), not the pool's
spread, not the sequence.

### By-product: the incumbent's own concentration, never previously reported

S8-11's **-0.2499** has median **-0.1085** with 63.5% of targets improved - but **ten targets
carry 61% of it**, and it falls to **-0.013 after dropping twenty**. This must be carried into
any final report: the headline gain is real and highly concentrated.

## S9-2. A torsion prior buys physical validity for free, and no accuracy

`s8/project.py`, `s8/test_project.py` (38 tests). Five penalties, all estimated from training
folds only: a class-conditional (GLY / PRO / PRE-PRO / GENERAL) smoothed Ramachandran
log-density with bilinear interpolation; the same over 20 residue types; `s5/torsion.py`'s von
Mises mixture reused as-is; a one-sided positive-phi barrier; and `ramah`, a hinged log-density
calibrated to the 5th percentile of real residues. Plus an agreement-weighted data term.

**No accuracy.** Chosen arm **3.202 against the incumbent 3.204**, and d(arm, fit) across
m = {25, 50, 75, 150} is +0.001 / +0.002 / -0.002 / -0.000 - flat zero. Top-10 share of the
delta is 1.034 with mean ex-top-10 of exactly +0.000; uniform null across both achievability
strata and all five folds; tuning bound |d| < 0.005 A. **It does not beat 3.204.**

**Validity, free.** Positive-phi for constrained non-glycine residues **17.5% -> 5.5%**, landing
on the 5.40% measured over real library windows and inside the 2.90-5.66% real band. It plateaus
across m (6.5 / 5.7 / 5.5 / 4.6%). Bond geometry exact to 5.3e-15 A; clashes, radius of gyration
and pseudo-angle unchanged in the fourth decimal.

**Dev-24, one pre-registered pass** (endpoints pinned before running; third pass on this family,
declared as such): `fit` reproduces S8-11 exactly (3.221, -0.255 [-0.501, -0.008], 16/8); the arm
is 3.220, d vs fit **-0.001 [-0.008, +0.006]**; positive-phi 17.66% -> 4.50%. Both pre-registered
endpoints met. Leakage bit-identical (0.000e+00) across four penalty families.

### Four findings beyond the brief

**1. The projection is degenerate.** A CA trace admits **two ideal-geometry torsion solutions at
near-equal distance to the average** - one conformationally plausible, one not. The objective
cannot distinguish them, which is exactly why removing the implausible branch costs nothing.
**Warm-starting from lambda=0 cannot cross between the branches; multi-start at every rung
flipped both priors from costing RMSD to free.** This may also explain the 0.036 A spread between
independent implementations of the same operator (`fit` 3.204, `P_fit` 3.240, `X_fit` 3.217) -
branch selection rather than convergence depth.

**2. Sparse support beats sophistication.** The crude one-sided barrier beat both smoothed
density priors, because its gradient is zero where phi is already negative. Mean-log-density
penalties drag *plausible* residues toward basin modes and pay +0.008 to +0.054 A. Prefer a hinge
with zero gradient in the allowed region over a log-density.

**3. `phi[0]` and `psi[n-1]` are never read by the builder.** The audit's 12.1% positive-phi
figure therefore averaged a coin flip into every target. This resolves the apparent contradiction
with `s8/predictor.py`'s 5.51% - the two were measuring different things. It is also why
`peptide_db` stores literal -60/-45 degrees at those positions.

**4. The convergence hypothesis is confirmed and does not help.** Within-target
rho(distance-to-average, positive-phi) = **-0.507**, negative on 98% of 109 targets, with the rate
running 1.1% -> 17.5% as the fit tightens. But early stopping reaches the real rate only at
maxiter 3, costing **3.204 -> 3.395**, while the prior reaches it at full convergence for -0.002.
**At matched plausibility the prior is 0.19 A better than early stopping** - do not trade
convergence for validity, you can have both. Note the BETWEEN-target correlation is **+0.559**,
the opposite sign: a Simpson's paradox that the within-target design was needed to see.

### Two things carried forward

- **The arm to deploy is `ramah@0.3`, not the one the pre-declared rule picked.** That rule ranked
  by RMSD subject to positive-phi - one bit of a 2-D distribution both arms satisfy. On the full
  distribution the hinge dominates: residues in *no* canonical basin fall to 2.90% (barrier
  14.34%, real window 1.59%) and symmetric KL to 6.78 (barrier 8.80), at a cost of +0.004 A
  [-0.005, +0.013]. **This is post hoc and labelled as such; the dev result stands on the pinned
  arm.**
- **No torsion prior fixes the over-helical conformation, and `ramah` makes it worse.** Alpha-R
  goes 63.9% -> 79.4% against a real window's 72.8%; beta stays at 12.7% against 22.4%; the
  pseudo-angle moves 0.74 of the 10.8 degrees it is short. **A per-residue marginal prior
  structurally cannot de-helicalise a chain, because alpha is its own mode.** That defect is
  inherited from projecting a *contracted* average (S9-1) and must be attacked at the averaging
  step, not downstream of it.

## S9-3. The relational channel is typicality wearing a tournament

`s9/relative.py` (10 resumable stages), `s9/test_relative.py` (45 tests). RelNet: a
set-conditioned Siamese twin with an antisymmetrised pair head -- DeepSets context over the
candidate set, a set-conditioned absolute score s_i, and a pair term f(a,b) depending on the
PAIR rather than either member, with logit = (s_a - s_b) + (f(a,b) - f(b,a)). Antisymmetric by
construction so argument order cannot be read; `s_only` prices the pair head exactly.

Instrument validity asserted as tests: 3.4540 / 1.7108 / medoid75 3.2822 / synth 3.2005 /
avg 3.0483; shipped in-band rho +0.131 against S8-8's +0.126; disto native percentile 36.68.

**Marginal skill over typicality is real and significant.** On identical cached pairs, in-band:
`cty_only` **+0.0838 [+0.0485, +0.1191]** over the LFO-best scalar centrality rule (75W/49L);
+0.044 to +0.068 over the shipped score, +0.067 to +0.078 inside its own top-75, every CI
excluding zero. In-band rho +0.147 to +0.210, **positive in all five folds -- nothing else in
the table is.**

**But after partialling radius of gyration the learned arms (+0.20 to +0.23) TIE the shipped
score's +0.212.** Most of the raw rho gain is compactness -- the fourth such instance in this
project.

**Native percentile: negative.** Best arm `s_only` **35.47 +- 2.48** (median 29.6), tied with
the distogram's 36.68 and behind inverse folding's 30.1. Not a better-specified objective.

**Ranker and filter: worth 0.00 A.** Against the 3.200 synthesis incumbent, best filter
`cty_only:filt75_fit` **3.219, +0.018 [-0.043, +0.080], 50W/76L**; stack-25 3.241-3.252. As a
ranker it is WORSE than the shipped score (3.498 pool-wide, 3.368 inside top-75). Against the
3.454 baseline the best is -0.201, but that is the incumbent's gain, not a new one.

**The null is arithmetic, not bad luck.** Priced against the oracle filter (reproduced: 2.095 A,
-1.106 [-1.289, -0.922]): the oracle needs **+0.672 of in-band recall** to buy 1.106 A; the
model delivers **+0.030** (0.328 -> 0.358). At that exchange rate its filter is worth ~0.049 A
against a 0.047 A paired SE. `cty_med75` hands over a strictly better, 22%-tighter set and
still returns 3.229. **Recall is the currency, and nothing measured has moved it by more
than 0.03.**

**Where the relational skill actually lives.** Permutation importance puts all of it in the
TOURNAMENT statistics -- centrality computed relationally (`tourn75` +0.0258, `tourn` +0.0119,
`d_typ` +0.0087). The displacement-profile SHAPE statistics -- where and how two candidates
disagree, the genuinely novel channel -- have importance <= 0.0004. And `full - s_only` is
**-0.003 [-0.020, +0.014]**: pair features add ~2 pp when the model has only centrality inputs,
and nothing once the absolute channels are present.

Leakage 0.000e+00 over 16 targets x 9 quantities. **No dev-24 pass spent** -- running a +0.018
null against a 0.35 A dev SE would spend a twice-used held-out set on noise.

**Two bugs found and fixed before any number was read**, the second serious: the tournament
feature's strict `<` made a self-pair score -0.5 so antisymmetry held only to the tie fraction;
and the "best pure centrality baseline" was LFO-selected over every name starting with `cty_`,
which swept in the LEARNED arms and made the control identical to the model it was controlling.

## S9-4. The loop is a contraction toward its own fixed point, not a multiplier

`s9/loop.py` (6 resumable stages), `s9/test_loop.py` (38 checks). Round 0 reproduces the
instrument exactly on all 126 targets (3.4540 / 1.7108 / avg 3.0483 / fit 3.2005).

**Structural retrieval genuinely moves recall - the first thing in three sprints that does.**
One round from queries of controlled quality (top three rungs ORACLE):

| query | query A | selected | recall_univ | in-band | d vs incumbent |
|---|---|---|---|---|---|
| *native* | *0.000* | ***1.925*** | *0.160* | *0.814* | *-1.275 [-1.490, -1.060]* |
| *library best* | *1.313* | *2.070* | *0.156* | *0.751* | *-1.130* |
| *pool best* | *1.711* | *2.210* | *0.135* | *0.703* | *-0.991* |
| **synthesis (deployable)** | **3.201** | **3.246** | 0.053 | 0.430 | **+0.046 [+0.003, +0.089]** |
| random member | 4.544 | 4.068 | 0.014 | 0.223 | +0.867 |

A native query moves universe recall **0.038 -> 0.160** and returns **1.925 A** - under the
2.0 A objective and under the 1.994 A perfect-distance-oracle ceiling. Nothing else measured in
this project has moved recall past 0.030.

**But closing the loop does not.** Pool-relative in-band recall FALLS 0.328 -> 0.204; the
trajectory is monotone worse (3.201 -> 3.246 -> 3.284 -> 3.307 -> 3.319 -> 3.323), `replace`
**+0.123 [+0.051, +0.194]**, 47W/79L; `fuse` +0.070, `direct` +0.103. This is NOT S8-6's
"the selector ignored the new candidates" null: the filtered set genuinely changed (round-1
overlap 0.095, answer changed on 126/126). No quality stratum pays.

### Coordinator correction to the agent's framing

The agent fitted **selected = 0.492 x query + 1.634** (r +0.961) and read it as "a multiplier on
predictor quality - every 1 A of query improvement buys 0.49 A", advising a re-run once a
predictor reaches 3.0 A. **The derivative is not the relevant comparison; returning the query
itself is, and that has slope 1.0.** Solving 0.492q + 1.634 = q gives q ~ 3.2.

**So the loop is a CONTRACTION toward its own fixed point at ~3.2 A, and is dominated by simply
returning the query for every query better than that.** At a 3.0 A query it returns 3.11; at 2.5
it returns 2.86; at a *native* query it returns 1.925. **It never improves on its input.** As
predictors improve the loop becomes more useless, not less. The agent's own break-even figure
(3.183 A against the pipeline's 3.2005) says the same thing, and is why the deployable arm is
significantly worse.

This is the **third independent instance of the same mechanism**: AMBER relaxation is compression
not refinement (S8-12), averaging shrinks in shape space (S9-1), and structural retrieval
contracts toward the pipeline's fixed point.

### Controls
- **Self-confirmation**: from a 4.54 A random seed the loop travels 0.96 A (92/126 improved) - a
  real search, not pure self-confirmation. But it lands +0.257 [+0.053, +0.461] short of the
  `replace` fixed point, **1.50 A away in structure space**, agreeing to 0.5 A on only 56/126.
  **Half of where it ends up is seed-determined.** `ctl_med` never leaves its seed - the pool
  medoid is its own fixed point.
- **Distortion inheritance is real and is NOT the mechanism.** rho(query curvature, candidate
  curvature) = **+0.907**, yet the loop improves the pseudo-angle 93.47 -> 98.41 while losing
  accuracy, and the distortion-corrected query retrieves **the most native-like candidates in the
  study** (103.55 against the natives' 103.99) and returns **the worst deployable answer**. The
  same inversion S9-1 found when correcting the angle directly.
- Geometry: every emitted structure CA-CA 3.804 +- 0.000, clash 0.00-0.04%, non-Gly positive-phi
  (excluding phi[0]) 5.7-10.1% against real windows' 5.40%.
- Leakage: 720 deployable quantities NaN-poisoned, worst |diff| 0.000e+00; 835,247 library pairs,
  0 violations. No dev pass spent; benchmark untouched.

## S9-5. Refinement fails, and the failure is located exactly: it is all selection

`s9/refine.py`, `s9/refine_FINDINGS.md`, 30/30 tests. Bar: the S8-11 synthesis at 3.2041 A.
**Best of 28 arms 3.1968, -0.0073 [-0.0273, +0.0126], 71/55 - null.**

With the pool removed as a confound, the decomposition is decisive:

| | CA-RMSD |
|---|---|
| ORACLE over the whole move space | **1.9714** (-1.233, 126/0, **55.6% < 2.0 A**) |
| ORACLE over the structures the search ACTUALLY VISITED | **2.4601** (-0.744, 126/0) |
| what the deployable objective SELECTED | 3.2390 (+0.035) |
| **selection loss** | **+0.7789** |

**The search walks past a sub-2.0 A structure on 39.7% of targets and hands back something worse
than where it started.** There is no search, space, geometry or starting-point problem.

**Two independent confirmations.** Exact enumeration of a 16,384-configuration subspace: the
*perfect* argmin of the objective is **+0.0137 worse than not moving**, and the approximate
searches beat it *because they fail to converge* (annealing reaches the exact optimum on 30% of
targets, VQE on 50%). A third of that subspace beats the incumbent on CA-RMSD while the objective
ranks the incumbent in its own top 6%. Separately, the budget curve shows the objective improving
50% while CA-RMSD is worse at **every** budget, including the first 24 evaluations.

**The new quantity.** Against an isotropic-move null at each arm's own step length, every arm
steers **52-139% better than a random direction**, every CI excluding zero. **The objective is
not blind - it is ~0.21 A of accuracy per A of travel short of breaking even**, with
rho(distance moved, damage) = **+0.973**.

### Four-component ablation, now conclusive
- **VQE** ties uniform random sampling (+0.0002 [-0.052, +0.052]) and significantly LOSES to
  annealing (-0.0637, CI excludes zero) - **at 4^n configurations, closing the "128 was
  enumerable" escape S8-9 left open.**
- **CVaR** null both ways (-0.033, +0.007). S8-9's collapse mechanism does not reproduce: at 2n
  parameters over 4^n states the ansatz cannot collapse to a delta (T=0 gives 8.57 bits with the
  tail against 8.66 without, of ~26).
- **Legacy** learns something real in its honest near-native role (+0.1384 rank skill against
  lg_all's +0.0915) and still does not earn a place (-0.0112); as a sole objective it is the worst
  arm measured (3.6087, +0.4047).
- **AMBER** is a validity stage: -0.0264 [-0.0377, -0.0151] for ~10^5-10^6 kcal/mol of strain.

The two positive credits (distogram +0.073, consensus +0.108) earn them only by LIMITING TRAVEL.

### The CVaR gradient defect, fixed and priced against an exact reference
Constant baseline **cos +1.000000 at 1.000x norm**; the shipped tail-only baseline **+0.655634 at
0.758x** - with **zero sampling noise**, so it is bias, not variance. The reimplementation
reproduces qansatz.cvar_gradient to 5.6e-17.

### A validity defect and a recurring inversion
The S8-11 synthesis' non-Gly positive-phi rate is **0.1973 against the natives' 0.0558** - the
torsion-space face of S9-1's 10.78 degree pseudo-angle deficit. And the inversion recurs: **arms
become more native-like geometrically exactly as they become less accurate** (angle 93.21 ->
96.92, positive-phi -> 0.0340). The Ramachandran channel this suggested was declared and killed
(+0.085 / +0.175 / +0.440 A).

Leakage clean (15 targets NaN-poisoned, worst |diff| 0.000e+00). No dev pass spent; benchmark
untouched. Two process faults are on the record: a rewritten chain script raced two runs, and a
Windows share violation on a checkpoint write killed a 12-arm run; both were killed, affected
outputs deleted and re-run, and _write now retries - no number comes from a raced write.

## S9-6. Synthesis is at its ceiling, and the binding constraint is the distance target

`s9/synth.py` (10 stages), `s9/test_synth.py` (29 checks), `s9/synth_FINDINGS.md`,
`s9/synth_{probe,contour,ext,hyp,curve,learn,conc,filt,leak}.json`.

**23 deployable operators on the 126-target instrument with the filter held at `sc|75`. The best
is 3.200 A against the incumbent's 3.201, CI +-0.011. Not one arm's CI excludes zero on the
winning side.**

| operator | d vs incumbent |
|---|---|
| GPA / Karcher mean | +0.015 |
| precision-weighted GPA | +0.021 |
| geometric median | +0.014 |
| manifold-constrained Frechet | +0.011 (robust variant -0.000) |
| weighted projection | +0.006 |
| distance-space consensus | +0.070 |
| equal-weight cluster mixture | +0.130 |
| **oracle per-residue weights** | **+0.027** |
| learned local-frame decoder | +0.014 to +0.024, negative out-of-fold R^2 |

Weighting loses **even when given the true per-residue errors as weights**. That is the sharpest
statement that the operator is not the problem.

### The contour ladder settles the over-curvature question

The coordinate average is **22.7% short of contour** (35.25 against 45.60 A). That missing length
has exactly two places to go:

- into **curvature** - pseudo-angle 121.6 -> 90.9, costing **+0.15 A**
- into **extension** - end-to-end 15.79 -> 20.13 against a native 15.67, costing **+0.72 A**

**The incumbent's 0.156 A is not slack. It is the cheaper of two forced options.** The
over-curvature I had been treating as a fixable artefact is a forced trade-off already resolved
optimally.

### Direction averaging is not a new operator - a coordinator suggestion, refuted with a proof

I proposed averaging unit CA-CA step vectors to preserve contour length by construction. It is
mathematically vacuous here: every candidate is an ideal chain, so all bond vectors share length
L and `mean(x_{j+1} - x_j) = L * mean(u_j)` - normalisation merely removes the L. Pinned to
1e-6 A by `t_dir_average_identity`. It scores **3.924 (+0.724)**, and it carries the *same*
121.62 degree pseudo-angle as the raw average, which independently proves **the over-curvature is
not manufactured by the projection**.

### Two further closures

- **The consensus is the optimum of every direction it can move along**: argmin(t) = 1.0, 1.0,
  0.0 across three native-free travel directions.
- **Multi-hypothesis has 0.51 A of ORACLE headroom and no rule reaches it at any separation.**
  The shipped score's rho with true RMSD among 2-4 A-separated cluster syntheses is +0.056 to
  +0.158, refuting the "coarse separation" escape left open by S8-8 and S8-9.

### What actually limits synthesis

**`O_dfit_nat` = 0.996 A.** Given the NATIVE distance matrix, the same projection machinery
returns **1.0 A with 80% of targets under 2 A**.

So the manifold, the projection, the optimiser and the chirality handling are **2.2 A from
binding**. **The distance target is the constraint** - and it is the only large oracle gap in
this study that is not a ranking problem in disguise. (Note the corresponding deployable path is
already measured and poor: S8-10's `mds` arm, distogram -> MDS -> torsion manifold, returns
4.159 A. The machinery is fine; the predicted matrix is not.)

### Discipline
Fidelity/legality is monotone on both knobs - there is no free legality, and early stopping costs
+0.039. The filter sweep was run last: a plateau, 0.049 A over m = 25-150, with operator ordering
identical at every m. Concentration: the incumbent gain is -0.253 with median -0.116, ten targets
carrying 60.4%, and -0.017 after dropping twenty; every tied arm shares that profile. Leakage:
168 outputs, worst |diff| 0.000e+00. Every answer arm is an exact peptide (CA-CA 3.803954938,
SD < 1e-15, zero out-of-range angles). **dev-24 deliberately not spent** - n=24 at SE 0.35 cannot
resolve a -0.000 effect. The 60-target benchmark was not touched.

**Five bugs caught by their own tests or CIs and documented**: a wrong `CA_CA` constant in the
validity assertion, non-equivariant local frames at chain ends, `mfit` tracking the wrong
objective, a fold-order/cache-order pairing bug in the leave-fold-out loop (the tell was a
+-0.45 A CI on a +0.024 A difference), and a rule-dispatch KeyError in `stage_hyp`.

## S9-7. Evolutionary information exists, is retrievable, and does not help

`s9/evo.py`, `s9/test_evo.py` (35 tests), `s9/evo_*.json`.

### A premise in my brief was false

I briefed this agent that "for 9-16mers a conventional MSA is hard - the sequence is short and
E-values are poor - so be honest if nothing is retrievable." **That was wrong.** Via UniProtKB
exact-peptide search -> parent protein -> UniRef50 -> alignment:

| quantity | value |
|---|---|
| targets occurring verbatim inside a known protein | **71/126** |
| targets with an aligned MSA | **69/126** |
| Neff >= 2 / >= 5 / >= 10 | 46 / 31 / 19 (max 134) |
| mean column entropy | 1.592 nats |

Evolutionary information exists for a majority of this instrument.

### Structural homology is indistinguishable from leakage at this length

The first exhaustive pass over `prots/` - 13,751 files -> 20,929 chains -> **3,519,987 residues**
(`fragment_db` samples this at stride 5) - scoring every 9-16mer window against all 126 targets:

**67/126 targets have a window at >= 0.6 identity (3,190 windows), all excluded by the leakage
filter.** Mean best identity **0.602**, median **0.600** - the median target's best structural
match sits *exactly on the leakage threshold*.

The shuffled-sequence control kills it: real 0.602 against shuffled **0.593**, d **+0.008
[-0.002, +0.018]**, 53W/46L - and **the shuffled sequence crosses 0.6 more often** (28.7 windows
against 25.3). Only **12/126** real sequences beat all five of their own shuffles.

**Homolog detection and leakage detection are the same operation for a 12-mer.**

### All three roles, 126 paired targets

| role | result |
|---|---|
| **filter** (feeding synthesis) | `pblos_sc\|75` **3.249** against the incumbent **3.200**, +0.048 [-0.017, +0.113]. The profile ties plain BLOSUM. Oracle filter 2.095, reproducing S8-13's 1.11 A headroom |
| **retrieval** | `pssm` achieves **the best perfect-scorer ceiling ever measured in this project - 1.963** against BLOSUM's 1.994 and S8-6's best of 22 keys (1.983) - and is **significantly WORSE on the deliverable**, +0.066 [+0.013, +0.118] |
| **conditioning** (LFO) | prof -0.099 [-0.197, +0.000]; against its own shuffled control only **-0.063 [-0.147, +0.020], 62W/64L** |

**The control that kills the mechanism:** the conditioning gain is **larger on the 80
empty-profile targets (-0.111) than on the 46 deep-profile ones (-0.077)** - the exact opposite
of what an evolutionary mechanism predicts.

### Direct answer

Not "sequence was the only thing tried". A second information source exists and was obtained.
**But the limit is the consumer, not the channel.** All three new signals moved pool statistics -
the profile key genuinely retrieves better structures, setting a new best ceiling - and none
moved the answer. **Sequence keys of any kind span 0.18 A in total against a 1.11 A oracle gap.**

### Audit and defects
Identity audit: 1,890 bit-identical checks under NaN-poisoning of every native, worst key change
**0.0**; **0/126** targets' own PDB entry present in the corpus; **0** library identity violations;
profile rows verified to contain sequences only. Dev(24) and benchmark(60) untouched.

**Five defects caught, each of which would have produced a publishable number**: an alphabet
mismatch scoring the wrong residue at every position; a pseudocount that broke the matched-arm
design (worth 0.61 A); its fix collapsing the log-odds key into all-ties; a cache shared between
two key definitions; and a `hash()`-seeded control that differed every run. Plus one conceptual
correction - **Neff = 1 means one effective *cluster*, not an empty MSA** (57 targets, where the
invariant holds at exactly 0.0).

ESM-in-parent-context is built, tested and one command away, deprioritised on my instruction
since ESM has already been tested three ways here.

## S9-8. The in-band ceiling is informational, not architectural - the recognition line closes

`s9/ceiling.py`, `s9/test_ceiling.py` (34 checks), `s9/ceiling_*.json`.

**Instrument validated before anything else.** The pools behind `s8/generate_univ`,
`s8/inband_cache`, `s8/integrate_chan`, `s8/consensus2_cache` and `s9/relative_cache` are the
SAME pool (worst |diff| 0.000e+00 on all 126 targets), which is what makes concatenating their
columns by index legitimate. The instrument reproduces 3.4540 / 1.7108 / shipped in-band rho
+0.1308 / synthesis 3.2005 / oracle filter 2.0947 exactly. Then the **oracle-feature control**:
one leaked true-RMSD column, pipeline otherwise unchanged -> **in-band rho +0.9559, recall
0.937. PASS.** The harness finds a signal that is definitionally there.

**Built:** 105 per-candidate features across seven channels (distogram plus a new per-shell
decomposition and predicted uncertainty; Legacy's 15; geometry/physics including pseudo-angle and
positive-phi; inverse folding; S9-3's set-relative block plus new per-candidate tournament
aggregates; a new per-candidate sequence block) and 13 context columns, each in z and percentile
view -> a 223-column matrix. Trained on in-band rows only, leave-fold-out.

| quantity | result |
|---|---|
| in-band rho | **+0.1960 +- 0.0292 raw, +0.1938 rg-partialled** (permutation floor -0.0015) |
| the shipped score | +0.1308 raw, **+0.2125 partialled** |
| in-band recall gain | **+0.0157** (0.328 -> 0.344) |
| native percentile | **39.45 +- 3.07** (median 28.8), argmin 6/126 |
| selected CA-RMSD | **3.2284, +0.028 [-0.049, +0.105]** vs the 3.200 incumbent, 56W/70L |

**Raw the kitchen sink wins by +0.065; partialled for compactness it LOSES.** Fifth compactness
instance. The sharp control - retrain with all seven rg columns *removed* rather than partialled
- gives +0.1997/+0.1974: it is not the column, it is the shape of what the pool varies in.

**The prediction discipline is the best of the sprint.** The recall gain was converted to an
expected RMSD gain of **-0.026 A using S9-3's exchange rate and written down BEFORE the RMSD was
measured**; measured +0.028, a 0.054 A miss against a 0.039 A paired SE. **The exchange rate held
out of sample** and can now be used to kill selection arms before they are built.

### Why this closes the question rather than adding to the pile

- **The capacity curve is MONOTONE DECREASING.** ridge@100 +0.196 > ridge@10 +0.196 > pair_lin
  +0.191 > gbt_s +0.166 > gbt_xl +0.143 > gbt_l +0.141 > mlp_l +0.119 > mlp_s +0.103 > mlp_m
  +0.079 - with the permutation floor RISING exactly as skill falls (-0.002 -> +0.045). **The most
  heavily regularised ridge beats a 512-512-512 net and a 1500-tree GBT.** `pair_lin` ties ridge,
  so the regression objective is not the limit either.
- **The whole is worse than its best part.** `rel` alone **+0.2244** against all-blocks +0.1960;
  **six of eight blocks have a negative marginal**, and dropping the new per-shell block is the
  largest single improvement available. Exactly what S8-9's fusion identity predicted.
- **The two missing channels were priced, not excused.** Amber (19 targets x 32 candidates
  cached): in-band rho **+0.0645** against its +0.4122 over the whole prefix - the fusion identity
  turns that into ~+0.004. The per-candidate sequence block, the deployable ESM stand-in, scores
  +0.021 alone.
- Concentration both ways: 10 targets carry **127%** of the aggregate rho gain, and dropping the
  10 the filter helps most leaves it **+0.107 A worse**.

**Verdict: in-band discrimination is informationally impossible with the features this project
can compute, not merely unattempted at scale. The recognition line should stop.** S9-4's
better-QUERY result (universe recall 0.038 -> 0.160, 1.925 A) remains the only measured route
past 2.0 A, and it is a generation result.

Two things to carry: `rel` alone is a better in-band ranker than every channel combined, so any
re-opening starts from set-relative centrality alone; and the exchange rate is now validated
out-of-sample. **No dev-24 or benchmark pass spent**, enforced by a test asserting every pdb id in
every output file is one of the instrument's 126.

## S9-9. Consensus is structurally the wrong operator for a refinement trajectory

`s9/traj.py` (8 resumable stages), `s9/test_traj.py` (26 checks), `s9/traj_FINDINGS.md`. The `sa`
generator is bit-identical to `refine.stage_visited`, asserted at 2.4601 / 3.2390 before any
stage runs.

**The precondition passes - the visited set is NOT the incumbent's neighbourhood.** Mean pairwise
spread 1.2076 A; the best visited structure sits **1.5773 A from the incumbent on 77.8% of
targets**; best-at position 0.501 of the walk. As a SET it beats the pool the incumbent averages
(mean **3.306 against 3.551**, -0.245 [-0.304, -0.186], 71.4% of visited structures better than
the top-75's mean) while its best is far worse (2.460 against 1.711).

**All three variants null**, 51 arms over 7 generators. Best `sa01/t_thin75` **-0.0022 [-0.0133,
+0.0090], 63/63**; trajectory alone +0.0097; trajectory union pool +0.0026; damped 4x +0.0331;
objective-free **+0.1493**; uniform **+0.2450**. **Variant 3 fails backwards**: the looser
generators visit *better* structures (rand 2.348, seed 2.380, iso 2.394 against sa's 2.460) and
their consensus is worse.

### The arithmetic that settles it

The S8-11 error decomposition run on both sets for the first time - and it *predicts* every
pre-projection number to 0.06 A:

| set | member | bias (the averaging ceiling) | independent | bias share |
|---|---|---|---|---|
| filtered top-75 (pool) | 3.6148 | **2.8148** | **2.0888** | 0.620 |
| visited (512) | 3.3292 | **3.1006** | 1.0598 | 0.866 |
| visited, objective top-75 | 3.2723 | 3.1927 | 0.5874 | **0.952** |

**Consensus converts INDEPENDENT error into accuracy. The trajectory has almost none, and its
bias floor is WORSE than the pool's.** Refinement narrows a set around one hypothesis - good
members, correlated errors; consensus needs members wrong in different directions. **Each
operator is already on the only set it can use.** The one intervention this endorses - seeding
the 64 chains at 64 different pool members - works mechanically (independent 1.06 -> 1.63, share
0.866 -> 0.768) while raising the bias floor 3.10 -> 3.22, which is exactly what the arms do.

### Two independent confirmations and two corrections

- **The 0.21 A-per-Angstrom exchange rate belongs to the objective, not the readout.** Across all
  51 arms `damage = +0.2182 x (A moved) - 0.0699`, spearman +0.9048 - against refine's **+0.2105**
  for the *argmin* readout over 28 arms. **Two readouts sharing no code pay the same rate.**
- `refine.Archive` had already computed consensus over the objective's top-64 on all 28 arms and
  never reported it (best 3.1993, worst 3.5063); now read out.
- **The positive-phi "defect" belongs to the PROJECTION, not the consensus.** Handed the *native
  CA trace*, the arms' own projection returns **0.2513** positive-phi (0.1737 multi-start) -
  higher than every arm. The natives are CA-only, so the earlier 0.1973-against-0.0558 comparison
  set a projected trace against a real backbone.
- Multi-start finds a nearer manifold point every time (objective spread 0.03-0.07 A) and buys
  **0.0000 A** of accuracy.

Every emitted structure: bond 3.804 A, SD < 0.001. Leakage clean over 12 targets x (21 arms +
trajectory), 0.000e+00. No dev-24 pass spent; benchmark untouched.

**Direct answer: the set contains information and consensus is structurally the wrong operator
for it. The quantity left to attack is the BIAS FLOOR - a retrieval problem, not a refinement,
search, readout or projection problem.** And retrieval is itself already closed (S8-6, S9-7).

## S9-10. THE FINAL BENCHMARK PASS: the gain does not replicate

`s9/final.py` (7 resumable stages), `s9/test_final.py` (35/35), commit `e58b675`.
**The 60-target benchmark had never been touched by any sprint. This was its single
pre-registered pass.** Every constant was fixed by prior work on the tuning instrument and
pinned literally by `t_preregistered_constants`; nothing was tuned, no variants were tried, and
there was no iteration after seeing a number.

Pre-registered system: K=500 BLOSUM62 pool over out-of-fold peptides plus
`distogram._fold_fragments` -> distogram top-75 -> coordinate average + multi-start projection
onto the ideal-geometry manifold with `ramah@0.3` -> AMBER ff14SB/GBn2 restrained relaxation
(k=10, converged).

### Headline, n=60

| | mean | median | SD | min | max | <2.0 A | <1.5 A |
|---|---|---|---|---|---|---|---|
| **FULL SYSTEM (+AMBER)** | **2.9610** | 2.9466 | 1.4806 | 0.465 | 6.005 | 0.333 | 0.233 |
| stage 3 only (declared ablation) | 2.9357 | 2.8957 | 1.4749 | 0.459 | 5.966 | 0.350 | 0.233 |
| shipped distogram argmin | 2.9507 | 3.0092 | 1.4972 | 0.646 | 6.370 | 0.300 | 0.267 |
| top-75 best (filter ceiling) | 2.0267 | 1.8957 | 1.2299 | | | 0.500 | 0.400 |
| pool best (achievability ceiling) | 1.4666 | 1.2298 | 0.8195 | | | 0.717 | 0.583 |
| perfect-distance oracle (diagnostic) | 1.7994 | 1.4870 | 1.0603 | | | 0.600 | 0.500 |

### The result

**full vs shipped: +0.0103, 95% CI [-0.1596, +0.1803], 31W/29L, paired SE 0.0867.**
Ablation vs shipped -0.0150 [-0.1843, +0.1543]. AMBER's own cost +0.0253 [+0.0171, +0.0336],
14W/46L. `ramah@0.3` vs no prior -0.0011 [-0.0138, +0.0116].

Selection gap **+1.4944** to pool best. The shipped baseline's gap is +1.4841 - **the
architecture closes none of it.**

### Did it replicate?

**The mean did. The effect did not, and that is the finding.**

| | tuning n=126 | dev n=24 | benchmark n=60 |
|---|---|---|---|
| system mean | 3.204 | 3.221 | **2.9610** |
| shipped baseline | 3.454 | - | 2.9507 |
| pool best | 1.711 | - | 1.4666 |
| **paired gain vs shipped** | **-0.250** [-0.383, -0.117] | **-0.255** | **+0.0103** [-0.1596, +0.1803] |

The absolute mean is 0.243 A *better* than tuning, but that is a property of the target set, not
the architecture: the benchmark's baseline is 0.50 A lower and its ceiling 0.24 A lower. The
quantity the architecture claims is the **paired gain**, and it is zero. At a 0.0867 paired SE
the pass had ample power to see -0.250, and **both the tuning and dev gains fall outside the
benchmark's CI.**

**This is not a broken pipeline.** Every mechanical quantity transferred exactly: projection cost
+0.168 A (tuning 0.156), the prior null at -0.001 (tuning -0.002), non-Gly positive-phi 4.93%
(tuning 5.5%), and S9-1's over-curvature reproduced at 94.75 against the natives' 104.37. It is a
genuine non-replication of the one accuracy result the programme was resting on, on data that had
never been looked at.

### The concentration was the tell

31 targets gain 13.71 A and 29 lose 14.33 A, so a percentage share of their 0.62 A difference is
not a real quantity and the code says so. The absolute statement is decisive: the top 5 gains
carry -6.27 A, top 10 -9.44 A, top 20 -12.56 A, and **after dropping the ten largest gains the
system is +0.2013 A WORSE than the baseline** (+0.3296 after twenty). By fold the sign flips:
-0.244, -0.100, +0.071, +0.349, -0.127.

Every agent that measured concentration on the tuning instrument reported the same shape - ten
targets carrying ~61% of the gain, falling to -0.013 after dropping twenty. **It was reported
as a caveat throughout. It was in fact the result.**

### Geometry validity

| | full | stage 3 | library window | native |
|---|---|---|---|---|
| CA-CA mean | 3.8654 | **3.8040** | 3.8073 | 3.8078 |
| CA-CA SD within chain | 0.0221 | **9.0e-16** | 0.0204 | 0.0119 |
| pseudo-angle | 94.90 | 94.75 | 100.00 | **104.37** |
| fraction out of 75-150 deg | 0.0026 | 0.0000 | 0.0000 | 0.0000 |
| clash (non-local < 4 A) | 0.0000 | 0.0002 | 0.0015 | 0.0047 |
| rg | 6.143 | 6.084 | 6.111 | 6.199 |

Non-Gly positive-phi excluding `phi[0]`/`psi[n-1]`: **4.93%** pre-AMBER against the 5.40%
real-window reference; 7.80% post-AMBER. AMBER moved 0.246 A mean and removed a median
**7.45e4** kcal/mol (39/60 above 1e4).

### Leakage audit - and a real defect found

- **397,592** library alignments; worst target-vs-member identity **0.583**; **0 violations**;
  0 self-in-library; 0 same-fold members.
- Benchmark intersected with tuning-126 and dev-24: **0 pdb ids, 0 sequences, 0 identity
  clusters** - all six intersections empty.
- **60 targets x 6 quantities NaN-poisoned, worst |diff| 0.000e+00.**
- **NEW DEFECT, found and priced.** The project's filter is member-level and `identity`
  normalises by the LONGER sequence, so a long fragment can carry a target's exact k-mer and
  still align below 0.6. **13 of 60 targets have a pool window at >=0.6 identity, two (8Y3S,
  8ZG3) at exactly 1.0.** The control - drop those windows, re-retrieve to K=500, identical
  downstream - gives **+0.0030 [-0.0002, +0.0061]** overall, +0.0137 on the 13, and **exactly
  0.000000 on the 47 clean targets.** The result is not carried by it, but the defect is real and
  applies to the filter everywhere in this repository.

### Process notes
The memory start gate was raised from 78% to 86% for these runs (an untouchable Chrome holds
~5.5 GB, putting the idle box at 82-83%); the 90% stop gate was unchanged and never fired. One
bug was fixed before any number was read: `distogram.train_fold`'s feature-width probe was
falling through to the 1.5 GB `esm_cache.npz` and taking the box to 94%, fixed with the same
guard `s7/poolsize.py` uses, with a hard error rather than a silent zero.

### What this means for the programme

**There is no validated accuracy improvement over the shipped baseline.** The honest held-out
state is: shipped distogram argmin **2.9507**, the synthesis architecture **2.9610**,
statistically indistinguishable. The architectural work of sprints 8 and 9 produced a system that
is physically better-behaved (exact bond geometry, native-like torsion statistics, strain removed)
and no more accurate.

The mechanistic findings stand - they were measured on the tuning instrument and most are
oracle decompositions that do not depend on the gain replicating. What does not stand is the
claim that filter-then-synthesise beats the baseline.

<!-- ==================================================================== -->
<!-- Sprint 10: verbatim from s10/FINDINGS.md -->
<!-- ==================================================================== -->

> *Merged verbatim from `s10/FINDINGS.md`. Not edited.*

# Sprint 10 findings

Entering state: **no validated accuracy improvement over the shipped baseline.** The Sprint 8/9
synthesis architecture scored 2.9610 against the baseline's 2.9507 on the untouched 60-target
benchmark - paired gain **+0.0103 [-0.1596, +0.1803]** - against a tuning gain of -0.250 (S9-10).

## S10-1. The attrition ledger: the loss is set PURITY, and the failure is a minority class

`s10/forensics.py`, `s10/test_forensics.py` (57/57), `s10/forensics_*.json`.

### The ledger, 126 tuning targets

| stage | mean | median | p10 | p90 | worst |
|---|---|---|---|---|---|
| library best | 1.313 | 1.247 | 0.432 | 2.349 | 2.912 |
| pool best (K=500) | 1.711 | 1.714 | 0.564 | 2.886 | 3.952 |
| filtered-75 best | 2.306 | 2.262 | 0.687 | 4.225 | 6.921 |
| coordinate average | 3.048 | 2.837 | 1.154 | 5.553 | 8.069 |
| projection | 3.201 | 2.992 | 1.162 | 5.626 | 8.266 |
| AMBER k=100 | 3.212 | | | | |

| transition | loss | median | -top10 | -top20 | per-fold range | % worse |
|---|---|---|---|---|---|---|
| retrieval | +0.397 | +0.311 | +0.331 | +0.278 | +0.33…+0.45 | 86% |
| **filter** | **+0.595** | +0.250 | +0.403 | +0.275 | +0.24…+0.75 | 67% |
| **operator** | **+0.742** | +0.490 | +0.603 | +0.509 | +0.65…+0.83 | 96% |
| projection | +0.152 | +0.099 | +0.133 | +0.109 | +0.13…+0.22 | 86% |
| AMBER k=100 | +0.011 | +0.009 | +0.008 | +0.007 | +0.008…+0.014 | 83% |

**The ledger closes to 0.00e+00 - and the agent correctly notes that is arithmetic, not evidence**:
a telescoping sum of best-available differences closes by construction. Stages 4 and 5 were
verified independently: projection +0.1522 [+0.1184, +0.1859] with multi-start buying a null
+0.0025; AMBER +0.0111 [+0.0085, +0.0136] at k=100, reproducing S8-12 on a different object.

### The pool's best member is deleted on 67.5% of targets - and that is NOT the mechanism

Survival into the filtered top-75 is **32.5% (41/126)**; median filter rank **134 of 500**, p90 399,
with 35 targets ranking it in the bottom half. In-band recall into the top 75 is 0.302 against the
0.150 a skill-free filter keeps by accident.

But operator loss is **+0.729 on the 41 survivors and +0.748 on the 85 deleted - identical** - and
re-inserting the pool's best into the shipped top-75 is worth only **-0.026 [-0.033, -0.018]**
(its whole oracle top five, -0.103). **Absence is not the problem.**

### Dilution, decisively - and its ceiling is the closed in-band problem

`keep@k`, replacing the k worst members with copies of the best: **-0.287 [-0.324, -0.250] at
k=10 on 126W/0L**, -0.615 at k=25, bottoming at -0.911. The random-member control gives only
-0.111/-0.155 and turns positive past k=50. The near-native member is present and weighted 1/75.
**But the probe's ceiling is the best member alone (2.305), so un-diluting it IS the in-band
ranking problem S9-8 measured as informationally impossible.**

### The bound that closes the synthesis line

**Oracle per-CANDIDATE weighting: 2.087 A raw, 2.160 A through the projection** (median 1.714,
58% of targets under 2.0). **A perfect re-weighting of the set the pipeline already builds does
not reach 2.0 A on the mean.** Adding oracle per-candidate alignment gives 1.893/1.940; the 2x2 is
additive (weights 0.961 A, alignment 0.234 A, both 1.155 A). The oracle weight is a soft selection
(n_eff 15.7, 0.167 on the best member against a uniform 0.013, rho(w, RMSD) -0.651). The bias
floor reproduces S9-9's 2.815 exactly.

### What actually governs the loss: set purity, not set ceiling

Twelve counterfactual sets through the unchanged pipeline:

**operator loss = 0.741 x (set mean quality) - 1.673,  R^2 0.912**

**Averaging becomes a GAIN operator below a set mean of 2.26 A; the shipped set is at 3.55 A.**
Hence +1.685 (unfiltered K=500), +0.742 (score top-75), +0.252 (oracle top-75), **-0.215 (oracle
top-3)**. And the deployable set containing the best available member anywhere - unfiltered K=500
- returns the **worst** answer (+0.393). Note the constant -1.673 is numerically the same as
S8-5's selection law, which is worth flagging as suspicious coincidence or shared cause.

### The failure is a minority class, not a uniform leak

**18 targets keep NONE of the near-native band and return 6.016 A against a 2.284 A pool best;
36 keep less than chance.** Per-target rho(score, true RMSD) is **+0.107 on those 18 against
+0.645 on the other 108**, and **rho(filter skill, final answer) = -0.811**. The deleted candidate
is **not an outlier** - centrality percentile 56.6, rg percentile 47.4. Separately,
rho(score, pool centrality) = **+0.494**, a sixth compactness instance.

### A mechanism for the benchmark non-replication

`fit - argmin` across recall bands: **+0.008 / -0.089 / -0.213 / -0.426 / -0.321**. **The
architecture's -0.253 tuning gain is produced entirely by targets where the filter works, and is
zero where it does not.** Both termini are explained: pool -> synthesis +1.490 (98% worse),
pool -> argmin +1.743 (100%).

**Every stage cost survives dropping its top ten and top twenty and holds sign in all five folds.
The architecture's own gain does neither** (-0.253 -> -0.167 / -0.116, top-10 share 0.39, per-fold
-0.120 to -0.424). The ledger is robust; the accuracy claim was not.

### Deployable lever tried and killed
The m-plateau is two opposed slopes - widening is **-0.49 A where recall is 0** and **+1.73 A
where recall >= 0.5**. Oracle per-target m is worth **-0.229 [-0.296, -0.163], 80W/0L**, and no
deployable rule reaches it: the best of 39 leave-fold-out thresholds over 13 native-free
statistics is -0.031 [-0.079, +0.017], sign-flipping by fold and turning positive when its ten
best targets are dropped.

### Identity audit -- WITHDRAWN, see S10-4
This section originally reported "0/126 tuning pools contain a >=0.6 window (max 0.583)" and
"gapped identity exceeds positional by up to 0.417". **Both are void.** They were measured against
a permuted alphabet: the cached window banks are encoded with `s7.audit`'s alphabet and were
decoded with `s9.evo`'s, which agree on only three letters. The correct figure is **13/126**
(S10-4). The ledger numbers are unaffected -- the leak prices at +0.0004 A and exactly 0.000000 on
the 113 clean targets -- but the audit itself must not be cited.

### Defects found and fixed, each caught by its own test
A transposed rotation in the weighting solve - it made the "bound" worse than the equal-weight
mean on 5+ targets and moved the headline from 2.648 to **2.087**; `sc75_plus_top5` evicting a
surviving top-five member on 2 targets; a missing retry in `_atomic` that killed a stage at 99/126
on a Windows share violation; and an unsafe positional screen in the identity audit. Leakage
0.000e+00 over 30 deployable arrays. dev-24 and the benchmark untouched.

### Which stage to attack
**The filter - specifically its 18-target zero-recall class - not the operator and not the
average.** The operator's +0.742 is not a property of averaging; it is 0.741 x set purity and goes
negative on a clean set. Give the same operator an oracle top-75 and the answer is 2.094; a
threshold-drawn band-75 gives 2.056. Nothing downstream has that headroom: perfect weighting caps
at 2.160, perfect projection at 0.152, AMBER at 0.011. The caution is that the required
intervention is set purity, i.e. in-band ranking, which S9-8 closed - **so the honest target is
the 18 zero-recall targets, where the distogram has no skill at all and a DIFFERENT signal, not a
better ranker over the same one, is what is missing.**

## S10-2. The projection path is closed: the required matrix is 2.7x better than anything achievable

`s10/dcurve.py`, `s10/test_dcurve.py` (114 checks), `s10/dcurve_*.json`. Instrument reproduced
before anything was believed: shipped argmin 3.454, pool best 1.711, incumbent synthesis 3.201,
`O_dfit_nat` 0.996. Leakage 140 checks, worst |diff| 0.000e+00. Benchmark untouched, dev unspent.

### The quality-transfer curve

`O_dfit_nat` = 0.996 A was the last unexploited oracle gap in the programme: given the NATIVE
distance matrix as a projection target, the existing machinery returns 1.0 A. This measures how
good a matrix has to be.

**Matrix quality required for sub-2.0 A mean output**, as input MAE and Pearson r against the
native:

| error model | MAE @ 2.0 A | r @ 2.0 A | MAE @ 1.5 A |
|---|---|---|---|
| iid on pair distances | 1.328 | 0.892 | 0.637 |
| correlated coordinate field | 1.309 | 0.925 | 0.680 |
| shrinkage toward the pool | 0.961 | 0.926 | 0.446 |
| **the distogram's own error shape** | **0.872** | **0.940** | 0.443 |

**Where the shipped distogram sits: MAE 2.339, r 0.696 -> 3.779 A output, 11.1% under 2 A.**

**The gap is a 2.7x cut in distance MAE** and r 0.696 -> 0.940. The best sequence-conditioned
distance MAE ever recorded in this project is ~2.11 (S7-11), so a **perfect** member of that
family is still **2.4x short**. The best deployable matrix available - the filtered pool's own
mean - is MAE 1.933 -> 3.270 A, needing 2.2x.

### Shape prices this path too, and it was proven causally

At matched MAE 1.0 the four models give **1.784 / 1.756 / 2.054 / 2.178 A** - the realistic error
shape is the **worst** one. Demonstrated causally on a deployable arm: `recal`, a leave-fold-out
shell recalibration, **lowers** the distogram's MAE to 2.222 and **raises** its output to 3.881,
because it buys MAE by shrinking. That is the shrinkage mechanism for the seventh time in this
programme.

### Three converging lines put the constraint on the matrix

- **Operator**: a chirality-guarded multi-start is 0.60 A better at perfect input (0.991 ->
  **0.394 A, 95.2% under 2 A**) and moves the *requirement* by only 0.04 A; at the shipped
  distogram's input it is worth 0.007 A.
- **Hybrid target** (average-plus-matrix, the second question): the mixing weight's optimum is the
  endpoint - 3.201 (alpha=1) -> 3.203 -> 3.208 -> 3.228 -> 3.270, monotone, with the distogram side
  worse. `pool@1` reproduces the incumbent at +0.000 [0, 0], 0W/0L.
- **Eleven deployable matrices** through the same projection: best is the pool mean at 3.270; the
  distogram 3.779; predictive median 3.788; MDS rank-3 truncation 3.965 - **this last is where
  S8-10's 4.159 came from, since the truncation raises the matrix's MAE 2.339 -> 2.681.**

### No deployable gain

Best arm `pool` = 3.270: -0.184 [-0.317, -0.052] against the shipped baseline, 76W/50L, negative
in all five folds - but **drop-top-10 = -0.038 [-0.143, +0.067]**, and it loses to the 3.201
synthesis by +0.069 [+0.019, +0.119]. Not an improvement. Geometry on every emission: CA-CA
3.8040 A, SD 0.0000, zero clashes, pseudo-angle 97.0 (pool) / 105.6 (disto) against the natives'
103.99.

### Where a sufficient matrix DOES exist - and why it is unreachable

The best convex combination of the K=500 pool's **own** matrices is **MAE 0.536 -> 1.682 A, 62.7%
under 2 A**, using about **9.6 of 500** fragments. Restricted to the filtered top-75 it is MAE
1.210 -> 2.447.

**So the required matrix is inside the library's convex hull** - reachable only by weighting 500
candidates correctly, which is the retrieval/ranking axis closed by S8-6 and S9-7. And a per-pair
Bayes combination of pool distances with the distogram's own risk table is **monotone in the wrong
direction** (1.933 -> 2.002 -> 2.230): **the predictor adds no per-pair information the pool
lacks.**

### Two defects, one of which corrects the coordinator

**The distance objective is exactly mirror-blind** (pinned to 1e-9). Raw lowest-objective
multi-start - which I propagated from S9-2 to several agents - **selects enantiomers** and is
worse than a single start at **every** rung (1.363 against 0.991 A at perfect input). S9-2's
prescription is correct for the **coordinate** projection, which is chirality-sensitive, and wrong
for any **distance-based** objective. Handedness must come from the L-signature.

**Identity audit conflict, unresolved.** This study reports **10/126 targets with a K=500 pool
window at >=0.6 identity, four at exactly 1.0** - three of the four sitting 2.3-4.1 A from their
own native, i.e. a sequence leak carrying little structure, with 1CEK the one real case. Every arm
therefore also reports a clean-116 block. **This contradicts `s10/forensics.py`, which audited the
same instrument exhaustively and reported 0/126 pools at >=0.6 (max 0.583).** The two may differ
in pool construction or identity convention - note the benchmark pass (S9-10) found `identity()`
normalises by the LONGER sequence, which makes member-level filtering leaky for long fragments.
**Neither number should be treated as settled until this is resolved.**

### Verdict

**This retires the last unexploited oracle gap.** The 2.2 A between `O_dfit_nat` and the deployed
system is real, but crossing it requires a matrix 2.2-2.7x more accurate than anything in this
project, on the one error shape that prices worst - and the only construction that reaches it is a
correct weighting of 500 pool candidates, which is the closed ranking axis in a different costume.

## S10-3. Physics does not supply decorrelated error - the force field agrees with retrieval, target-specifically

`s10/mdgen.py` (8 stages), `s10/test_mdgen.py` (52 checks), committed through `3621078`.
Genuine `amber14/protein.ff14SB.xml` + `implicit/gbn2.xml` through OpenMM 8.5.2, `CustomGBForce`
verified present. 32 targets (every 4th of the 126-target instrument, declared before running);
per target 5 native-free starts x 3 temperatures (300/450/600 K); restrained minimisation then
Langevin, 4 fs with 3-amu HMR, 2 ps equilibration and 10 snapshots over 14 ps. **480 trajectories,
4.55 CPU-hours.** Benchmark untouched, no dev pass, `DEV_ARM` asserted None.

### The decomposition (ORACLE)

| set | m | member | **FLOOR** | HULL | share | <2A | modes |
|---|---|---|---|---|---|---|---|
| pool75 | 75 | 3.333 | **2.613** | 1.604 | 0.628 | 0.56 | 5.4 |
| pool225 | 225 | 3.580 | 2.612 | **1.124** | 0.520 | 0.59 | 6.4 |
| pool500 | 500 | 4.564 | 2.864 | **0.828** | 0.376 | 0.72 | 6.3 |
| md_all | 150 | 4.562 | **3.172** | 1.589 | 0.459 | 0.47 | 4.2 |
| md_T300 | 50 | 4.465 | 3.082 | 1.815 | 0.458 | 0.41 | 2.6 |
| md_ext | 30 | 6.368 | 5.716 | 4.538 | 0.768 | 0.03 | 3.9 |
| noise_all (control) | 150 | 4.796 | 2.933 | 1.143 | 0.358 | 0.66 | 3.9 |
| union_bal | 225 | 4.041 | **2.747** | 1.261 | 0.426 | 0.59 | 4.6 |
| union_noise | 225 | 4.186 | 2.662 | 1.027 | 0.372 | 0.72 | 4.3 |

The subset reproduces S9-9's reference shape. Paired against the pool's own floor, `md_all` is
**+0.559 [+0.273, +0.845], 10W/22L, 0 of 5 folds lower**, and - the reverse of the pattern that
failed on the benchmark - **dropping the ten targets MD helps most makes it worse** (+0.559 ->
+0.966). A uniform negative, not a concentrated one.

### The union arithmetic

Floor of `a x pool + (1-a) x MD`: **3.172 / 2.933 / 2.747 / 2.633 / 2.613** at a = 0.00 to 1.00.
**Monotone in the pool's weight; the minimum is the pool alone.** Choosing the mixture weight
per-target WITH THE NATIVE IN HAND buys 2.475 against 2.613 - **0.138 A for an oracle.**

### Why: the force field agrees with retrieval about the fold

**cos(pool mean error, MD mean error) = +0.6923 +- 0.2595, positive on 97% of targets**, against a
shuffled-target control (110 cross-target pairs of equal length) of **+0.0287 +- 0.2702**. The
alignment is **target-specific**. And it is not shared contraction: the outward-field direction
explains only ~15% of either mean field, and removing it *raises* the cosine to +0.7319.

**ff14SB/GBn2's basin around a retrieved fragment is centred on very nearly the structure the
retrieval pool is centred on - physics agrees with retrieval about the fold, including where
retrieval is wrong.**

### Two controls that close it

**Matched cardinality.** A convex hull can only shrink as members are added, so a 225-member union
beating a 75-member pool proves nothing. `pool225` - the same deployable filter taken deeper -
reaches hull **1.124** against the union's **1.261**, better on 22/32 targets; `pool500` reaches
**0.828**. **150 more retrieved fragments buy more hull than 150 MD conformers.**

**The force field priced against matched noise.** An MD ensemble differs from its start in two
ways: it moved, and it moved *the way the potential says*. Only the second is physics. The twin -
same start, same count, displaced by the same measured CA-RMSD, force field replaced by white
noise in torsion space - **wins on every axis**: floor 2.933 vs 3.172, hull 1.143 vs 1.589,
sub-2 A coverage 0.66 vs 0.47, and pool-union-noise (2.662) beats pool-union-MD (2.747).
**Independent per-residue noise decorrelates; the force field's motion is collective.**

**The mechanism, measured.** Participation ratio of each set's own spread (cap 32): pool75 **5.4**,
pool225 **6.4**, md_T300 **2.6**, md_T450 3.3, md_T600 4.2. A 14 ps trajectory moves 2.3-3.3 A but
along two to four collective modes, and every trajectory from the same start explores the same
ones. **Consensus can only average away error lying in the directions a set spreads in.**

### Not a budget artefact
A ten-fold truncation sweep moves the floor 3.270 (1.2 ps) -> 3.172 (12 ps); matched-size time
windows (first half vs second half, same member count) differ by 0.09 A. **The ensemble mean is
set by where the trajectories started, not by how long they ran.** The honest exception is
`md_ext`: 14 ps does not fold an extended 12-mer, so that arm IS budget-limited - and it is also
the least aligned with the pool (cos +0.402) and by far the worst (floor 5.716). **Being
differently wrong is not being usefully wrong.**

### The force field buys validity, and validity is not what the floor is made of

| set | CA-CA | pseudo-angle | non-Gly +phi | clash |
|---|---|---|---|---|
| NATIVE (CA trace) | 3.811 | 101.97 | (0.0558 ref) | 0.0062 |
| md_T300 | 3.871 | 102.62 | 0.156 | 0.0001 |
| md_all | 3.871 | **104.20** | 0.165 | 0.0001 |
| noise_all | 3.804 | 106.02 | 0.202 | 0.0107 |
| noise_syn | 3.804 | 95.13 | 0.246 | 0.0116 |

Pseudo-angle **104.20 against the natives' 103.99**, where the synthesis sits at 93.21; clash
0.0001 against the natives' 0.0062. Same lesson as S9-2's torsion prior, now for the potential
itself.

### The floor is a predictive instrument
**spearman(FLOOR, selected) = +1.0000 over eight arms**, with `pool` reproducing the incumbent to
**-0.0022 A** through independently written code, which is what makes every other row paired. Also
newly reported: the floor superposes onto the *native*, which no arm can do, so alignment alone
costs 0.14 A (pool) to 0.22 A (union) - **the floor is a bound, not a target.**

### Leakage and reproducibility
Six targets NaN-poisoned end to end - starting-point selection, the all-atom build, the restrained
minimisation, **the Langevin trajectory itself**, and all eight emissions: worst |diff|
**0.000e+00**, later extended to 8 arms with the result unchanged. Cells regenerated in a fresh
process match disk to 0.000e+00 (per-trajectory OpenMM contexts and crc32 seeds; a shared
integrator's random stream would have made a resumable study depend on its resume boundary, and
`hash()` on a str is salted per process and is not used).

**The identity-filter defect is priced and runs AGAINST the conclusion.** MD retrieves nothing and
is immune; every pool row inherits it, and the pool rows are the ones that beat MD. Over the 20/32
targets with pool best >= 1.0 A: pool75 floor 3.425 / hull 2.182, pool225 3.270 / **1.535**,
md_all 3.687 / 2.010, union_bal 3.433 / 1.719. The conclusion survives - the union moves from
losing to the pool to tying it, while `pool225`'s hull still beats it.

### Verdict and the two instruments worth keeping

The pool's bias floor is binding not because it is a retrieval limit a different distribution could
relieve, but because **the force field's own equilibrium distribution around a retrieved fragment
is centred where the retrieval pool is centred**, and because MD's spread is two-to-four-dimensional
where retrieval's is six. Every route was measured: more temperature (worse), more time (flat), a
genuinely different start distribution (much worse), better weighting (hull, beaten by more
retrieval), and the potential itself priced against matched white noise (beaten by the noise).

Two instruments outlast the negative:
- **`hull_floor`** - the distance from the native to a candidate set's convex hull, which
  **upper-bounds every consensus operator this project has tried or could try**. `pool500` reaches
  **0.828 A**.
- **the mode count** - a set is useful to consensus only in the directions it actually spreads in.

**The gap between `pool500`'s hull (0.828) and its uniform floor (2.864) is where the remaining
accuracy lives, and it is a weighting problem inside retrieval** - which S9-8 measured as
informationally closed. This converges exactly with S10-2's independent result that the best
convex combination of pool *matrices* reaches MAE 0.536 -> 1.682 A using ~9.6 of 500 fragments.

## S10-4. The identity-audit conflict adjudicated: an alphabet permutation, and the leak is worth +0.0004 A

`s10/idaudit.py` (stages cause / pools / audit / price / overlap / report), `s10/test_idaudit.py`
(15 tests), `s10/idaudit_*.json`, commit `49fc708`.

### CORRECTION TO S10-1 - its identity audit is VOID

S10-1 above records "0/126 tuning pools contain a >=0.6 window (max 0.583)" and "gapped identity
exceeds positional by up to 0.417". **All of those numbers are withdrawn.** They were measured
against a systematically permuted sequence.

**The cause is one line.** The cached window banks `S` in `s8/generate_univ/*.npz` are
integer-encoded by `s7.audit.encode`, alphabet `"ARNDCQEGHILKMFPSTWYV"`. `s10/forensics`'s
identity stage encoded the target with `s9.evo.encode` and decoded the windows with `s9.evo.AA` =
`"ACDEFGHIKLMNPQRSTVWY"`. **The two orders agree on three letters (A, S, T)**, so every identity
it computed - positional and gapped alike - compared the target against a permuted string.
`s9/evo.py` already documented this exact hazard in writing, and it was hit anyway. This is the
second alphabet defect of the sprint; the evolutionary study caught the same class of bug scoring
the wrong residue at every position.

```
1CEK, AISVLLAQAVFLL, against its own exact-sequence pool window:
  decoded in audit order   AISVLLAQAVFLL   positional 1.000  gapped 1.000
  decoded in evo   order   ALSYMMAGAYQMM   positional 0.308  gapped 0.462
```
`forensics_identity.json` records 1CEK at `pos_max` 0.3077 - the garbled value to four decimals.

One run of the same 126 K=500 pools, four ways:

| alphabet | measure | >=0.6 | =1.0 | max | reproduces |
|---|---|---|---|---|---|
| audit | gapped | **13** | 4 | 1.000 | **the correct audit** |
| audit | positional | 10 | 4 | 1.000 | S10-2 (dcurve) exactly |
| evo | gapped | 0 | 0 | 0.583 | S10-1 (forensics) exactly |
| evo | positional | 0 | 0 | 0.500 | forensics' `pos_max` column |

Ruled out by measurement rather than argument: pool construction is identical (21 pools rebuilt
from primitives, window banks and top-500 membership bit-identical to cache); `s5/inband.py`'s
non-stable argsort would move 21/500 members, a real difference but not this one; and
normalisation is irrelevant at window level because windows are target-length.

### True leakage status of the 126-target instrument

**13/126 targets carry a K=500 pool window at >=0.6 identity, four at exactly 1.0.** 21/126 have
one somewhere in the universe (86 windows: 52 from database peptides, 34 from fold fragments).

| pdb | n | fold | identity | # >=0.6 | RMSD of that window | pool best |
|---|---|---|---|---|---|---|
| 1CEK | 13 | 2 | 1.000 | 3 | **0.595** | 0.342 |
| 2FBU | 12 | 4 | 1.000 | 1 | 3.278 | 2.243 |
| 2P5H | 9 | 4 | 1.000 | 1 | 2.334 | 1.840 |
| 6B9K | 10 | 0 | 1.000 | 1 | 4.126 | 2.077 |
| 1RSW | 12 | 0 | 0.667 | 2 | 5.672 | 3.183 |
| 6S0N | 9 | 4 | 0.667 | 3 | 2.760 | 1.562 |
| 7N2I | 9 | 1 | 0.667 | 1 | 2.052 | 0.741 |
| 1NIZ | 14 | 0 | 0.643 | 2 | 2.005 | 1.310 |
| 1FUV | 11 | 2 | 0.636 | 2 | 4.155 | 2.414 |
| 6MK8 | 13 | 1 | 0.615 | 1 | 0.863 | 0.641 |
| 8FLP | 13 | 4 | 0.615 | 3 | 1.522 | 1.208 |
| 1N9U | 10 | 4 | 0.600 | 2 | 4.235 | 1.981 |
| 2BAO | 10 | 2 | 0.600 | 1 | 3.930 | 1.223 |

Three of the four exact-identity windows sit 2.3-4.1 A from their own natives - sequence leak
carrying little structure. **1CEK is the one that carries structure.** And **no flagged window is
ever its pool's best member**: dropping all 23 leaves pool best at 1.7108 exactly.

**The mechanism is confirmed on this instrument.** **573 library members across 58 targets clear
the shipped member-level filter while containing the target at >=0.6**, because
`peptide_db.identity` divides by the LONGER sequence. 1CEK's sequence sits verbatim inside the
25-residue peptide 1A11 at a member identity of 0.520.

### The price

Offending windows dropped before the top-K cut, re-retrieved to K=500 from the next-best BLOSUM
windows, downstream unchanged. Baseline reproduces the instrument and is bit-identical to
`s9/synth_cache`.

| arm | base | clean | diff | 95% CI | drop-top-10 |
|---|---|---|---|---|---|
| pool best | 1.7108 | 1.7108 | +0.0000 | [0, 0] | +0.0000 |
| shipped argmin | 3.4540 | 3.4540 | +0.0000 | [0, 0] | +0.0000 |
| synthesis `fit` | 3.2005 | 3.2009 | **+0.0004** | [-0.0004, +0.0013] | +0.0007 |

**Exactly 0.000000 on the 113 clean targets**; +0.0041 [-0.0039, +0.0122] on the 13 leaked ones.
The tuning gain moves -0.2535 [-0.3864, -0.1206] -> **-0.2531 [-0.3859, -0.1202]**.

### The decisive cross-tabulation

The top-10 gain targets - **2RUO, 9BAF, 6BJF, 2MP9, 1M23, 1TOR, 7QZV, 2NBC, 7N2I, 8HVS** - carry
**60%** of the total (S9's ~61% reproduced), with the gain falling to -0.109 after ten and -0.017
after twenty.

**Exactly one (7N2I) is leaked, and the top-20 overlap is also just 7N2I.** The 13 leaked targets
together carry **2%** of the gain - mean per-target **-0.052 against the clean set's -0.277**, i.e.
**the leaked targets gain LESS than average.**

**Concentration is the explanation for S9-10's benchmark non-replication. Contamination is not,
and is not even a contributing one.**

### What needs revising
- **S10-1's identity paragraph is withdrawn in full** (0/126, max 0.583, the 0.417
  gapped-vs-positional claim, and the 11/126 whole-library figure). Its *conclusion* that leaked
  windows are never the library best happens to be true, for an unrelated reason.
- **S10-2's 10/126 should read 13/126**, and its clean-116 blocks should be clean-113. Since the
  leak prices at +0.0004 A this changes no conclusion, only the count.
- **Nothing else moves.** S9-10's benchmark audit, its +0.0030 price, the -0.250 tuning gain, the
  +0.0103 benchmark non-replication and the concentration finding all stand exactly as recorded.

The tests pin the alphabet the window banks are written in, that `peptide_db.identity` normalises
by the longer sequence, that longer-normalisation hides containment (9/20 = 0.45 < 0.6), that the
real database admits a 1CEK-containing member through the holdout, and that positional identity is
a lower bound.

## S10-5. The reconciled bound ladder: the geometry is not the barrier

`s10/bounds.py` (one module, one convention, all 126 targets), `s10/test_bounds.py` (26/26),
`s10/bounds_*.json`, commit `1b662fc`. Every row is an ORACLE DIAGNOSTIC; none is deployable.
Instrument reproduced exactly first (3.4540 / 1.7108 / 3.0483 / 3.2005); leakage 50 checks at
0.000e+00; benchmark untouched.

**Convention**: native = model-1 CA trace; Kabsch CA-RMSD; the two alignments computed separately
and never mixed - `oa` (each candidate posed on the native) and `cf` (the pipeline's medoid frame,
one shared transform solved jointly with the weights); weights on the simplex by Frank-Wolfe with
a duality-gap certificate; projection = `fit_w` from 5 starts, branch chosen by the L-signature.

### The ladder, 126 targets

| kind | bound | raw | **EMITTED** | price | <2A |
|---|---|---|---|---|---|
| R | shipped argmin | - | 3.454 | | 0.21 |
| R | incumbent synthesis | 3.048 | 3.201 | +0.152 | 0.29 |
| N | uniform mean, top-75 | 2.815 | 2.954 | +0.139 | 0.30 |
| P | best member, top-75 | - | 2.306 | | 0.44 |
| N | best re-weighting, top-75 | 1.987 | 2.044 | +0.057 | 0.56 |
| D | best matrix combination, top-75 | - | 2.453 | | 0.43 |
| N | convex hull, top-75 | 1.802 | 1.840 | +0.039 | 0.62 |
| P | perfect distance ranker, K=500 | - | 1.994 | | 0.51 |
| D | best matrix combination, K=500 | - | 1.572 | | 0.70 |
| P | best member, K=500 | - | 1.711 | | 0.58 |
| N | convex hull, K=225 | 1.336 | 1.319 | -0.017 | 0.82 |
| P | best window, whole library | - | 1.313 | | 0.81 |
| D | native distance matrix | - | 0.394 | | 0.95 |
| **N** | **convex hull, K=500** | **0.953** | **0.853** | **-0.101** | **0.96** |
| N | affine span, K=500 | 0.035 | 0.064 | +0.029 | 1.00 |

### The R-1 disagreement, decomposed with no residual

```
forensics, common frame, projected gradient, 126t   2.087
  certified optimum (Frank-Wolfe), same problem     1.987   -0.100   SOLVER (a defect)
  convention: re-weighting -> true convex hull      1.802   -0.185   CONVENTION
  subset: 126 -> mdgen's every-4th 32               1.601   -0.200   SUBSET (benign)
= mdgen's published 1.604 (400 FW iterations)       1.603   +0.002
```

**The solver term is a real defect**: projected gradient at the 1/L step stops short on an
ill-conditioned Gram, and its error **grows with set size** (+0.092 at m=75, +0.193 at 225,
+0.296 at 500). The one candidate cause that would have been a genuine bug - a frame-convention
mismatch - **is excluded by measurement**: `H_hull_frame` and `H_hull` agree to 2.7e-13. R-2 and
R-4 (the bias floors and hull500) are pure subset differences, both reproduced exactly.

### A coordinator premise, refuted

I asserted that "oracle weighting cannot be worse than the hull, since a hull floor is the minimum
over the simplex". **That assumes a containment which does not exist.** The sets
`{R * sum_i w_i Q_i X_i}` and `{sum_i w_i R_i X_i}` are **not nested** - posing each candidate
individually at its best destroys the relative geometry that averaging exploits. The hull wins on
the mean and on 369/378 cells and **loses on 9**, by up to 0.085 A, verified at 400,000
Frank-Wolfe iterations.

### Physical versus non-physical - the question that motivated the study

A hull point is ~20% contracted in contour and is not a peptide. Projecting it re-expands it, and
**for a tight hull that motion is toward the native, not away**:

| object | raw | emitted | price (95% CI) |
|---|---|---|---|
| uniform mean, top-75 | 2.815 | 2.954 | +0.139 [+0.102, +0.176] |
| convex hull, top-75 | 1.802 | 1.840 | +0.039 [+0.009, +0.068] |
| convex hull, K=225 | 1.336 | 1.319 | -0.017 [-0.042, +0.009] |
| **convex hull, K=500** | **0.953** | **0.853** | **-0.101 [-0.128, -0.074]** |

**The post-projection hull does not decay to 1.5-2.0 A as I feared - it improves to 0.853 A**
(median 0.759, 96% under 2 A, worst 2.17), every emission carrying an exact 3.804 A virtual bond,
clash 0.0038 against the natives' 0.0062, and an L-handed backbone. And **the manifold is not the
constraint anywhere**: the affine span holds the native at 0.035 A and projecting that lands at
0.064 A, so an ideal-geometry chain exists within 0.06 A of every native here.

**The geometry is not the barrier, and the closure S9-8 measured is informational.**

### Ceilings that need revising
- `s10/forensics`: **2.087 -> 1.987**, **2.160 -> 2.044**, **1.893 -> 1.802** (solver defect).
- `s10/mdgen`'s 1.604 / 1.124 / 0.828 are **32-target** figures; on the full instrument they are
  **1.802 / 1.336 / 0.953**. The 32 targets are simply easier.
- **S10-3's claim that `hull_floor` upper-bounds every consensus operator "this project has tried
  or could try" must be qualified: it does NOT bound an operator that ends in the projection, and
  every operator this project ships does.** (The projection can move a tight hull point *toward*
  the native.)
- Two smaller corrections: the perfect distance ranker must be ranked on |i-j| >= 2, and the
  Frank-Wolfe duality gap carries the gradient's factor of two - without it the certificate claims
  twice the convergence it has.

### Remaining headroom, authoritative

From the emitted 3.201 A to the best **emittable** structure the candidate set can produce:

- **within the shipped top-75: 1.36 A total** - 0.90 A by selection alone (2.306), 1.16 A by
  re-weighting the pipeline's own stack (2.044), the last 0.20 A only with oracle per-candidate
  poses (1.840).
- **within the retrieved K=500: 2.35 A total** - 1.49 A by selection alone (1.711), 2.11 A by
  re-weighting (1.094), the last 0.24 A with oracle poses (0.853).

**The filter costs about 1 A of reachable-in-principle accuracy, and every remaining lever is a
weighting or ranking problem inside retrieval.**

## S11-1. The four mandated components run genuinely, and three of them buy nothing measurable

`bench_results/fourcomponent_tuning126_w8.json`, 126/126 targets, cold, 8 workers, zero errors,
`quantum=True legacy=True amber=True`, `vqe_qubits=7 vqe_layers=3 vqe_iters=50 vqe_seed=0`,
`legacy_top=5`, K=500, multi-start projection, `maxiter=300`. Benchmark untouched; this is the
tuning instrument.

All four components genuinely execute and materially change the output. The quantum stage costs
9.54 CPU-s and Legacy 0.76 s; enabling the VQE widens the filter from 75 to 128 candidates
(`want = max(m, 2^n)`), which is why AMBER doubles to 3355 CPU-s and projection nearly triples to
2196 s. **The state is not collapsed**: `vqe_entropy_bits` = **5.9129** of 7, so the CVaR free
energy is producing the broad ensemble it is designed to produce rather than degenerating to the
argmin. Nothing here is decorative.

### The ablation table, all 126 targets, paired

| ablation | mean diff | 95% CI | W/L | zero in CI |
|---|---|---|---|---|
| CVaR-VQE selection vs the argmin it collapses to | **-0.1405** | [-0.2840, +0.0030] | 66/48 | **yes** |
| the state's distribution vs a UNIFORM ensemble over the same hypotheses | -0.0308 | [-0.1463, +0.0848] | **49/49** | **yes** |
| quantum-weighted coordinate average vs the uniform one | -0.0135 | [-0.0805, +0.0536] | 62/64 | **yes** |
| quantum synthesis vs the s9 CLASSICAL synthesis | **+0.0133** | [-0.0223, +0.0488] | 57/69 | **yes** |
| Legacy refinement vs the consensus neighbourhood's own centre | **+0.0295** | [-0.0351, +0.0941] | 46/54 | **yes** |
| AMBER validity stage | **+0.0207** | [+0.0140, +0.0274] | 40/86 | **NO** |
| **the whole four-component system vs like-for-like `s9/final.py`** | **+0.0169** | [-0.0192, +0.0530] | 58/68 | **yes** |

### What this says

**The quantum weighting is worth nothing over uniform.** The decisive row is the second: the
optimised state's distribution against a uniform distribution over the *same* 2^7 hypotheses is
-0.031 A at **49 wins and 49 losses** -- a dead heat, and the tightest possible statement that
the VQE's learned weights carry no information the hypothesis set did not already carry. The
same holds when those weights enter the coordinate average (-0.013, 62/64), and the quantum
synthesis is point-estimate *worse* than the classical one (+0.013, 57/69).

**CVaR selection's -0.140 A is real in sign and not significant.** 3.4540 -> 3.3135 at 66W/48L,
CI [-0.2840, +0.0030], missing zero by 0.003. It is the largest quantum effect measured, it
comes from CVaR's tail definition preventing the collapse to argmin, and on this instrument it
cannot be distinguished from chance. It should not be quoted as an improvement.

**Legacy refinement is negative.** +0.0295 A at 46W/54L. Consistent with the standing record
that Legacy is huge on single cells and worthless aggregated.

**AMBER is the only component with a significant effect, and it is a COST.** +0.0207 A
[+0.0140, +0.0274], 40W/86L -- inside the +0.011 to +0.026 A that S8-12 priced. That is the
documented trade: the validity stage buys a physically legal all-atom structure and pays about
0.02 A for it. It is reported, never hidden.

**The system as a whole does not beat the pipeline it extends.** +0.0169 A [-0.0192, +0.0530],
58W/68L against the like-for-like path. The synthesis gain that IS real -- `rmsd_arm` -0.239
[-0.373, -0.106] against the shipped argmin -- is classical, and it is present in both arms.

**So the four components are genuine, and three of the four are accuracy-neutral.** This is the
honest four-component result and it converges with the programme's bottom line rather than
disturbing it: the levers that move CA-RMSD on this instrument are retrieval, filtering and
consensus synthesis, and adding a quantum selector, a quantum weighting or a Legacy refiner on
top of them changes the answer by less than the measurement can resolve.


## S11-2. The projection is not a function of its input at angstrom resolution -- the reference disagrees with itself by 1.62 A

`core/project.py stability`, `verify/project_stability_partial77.json`, commits `f90fd81`,
`bc5cd22`, `304fd9a`. **77 of 126 targets** -- the run was stopped by the RAM gate at 92% and the
table is written to a `_partial77` filename with `complete: false` and `stopped_by` recorded, so it
cannot be mistaken for the full instrument.

### The control

Stage 3b projects a coordinate average `C` onto the ideal-geometry manifold. A CA-RMSD is
**exactly invariant under a rigid motion of the target cloud**, and the Ramachandran penalty does
not see the cloud at all. So projecting onto `R C + t` -- for a random rotation `R` with
`det R > 0`, never a mirror -- is *mathematically the same problem* as projecting onto `C`, while
every floating-point operation inside the Kabsch is differently rounded.

Running the **reference implementation** `s8.project` on both is therefore the reference's own
reproducibility, measured in the same units as any rewrite's divergence.

| quantity | value |
|---|---|
| worst self-disagreement (CA-RMSD) | **1.6246 A** |
| median | 0.0323 A |
| moved > 1e-3 A | 62 / 77 |
| moved > 0.01 A | 46 / 77 |
| moved > 0.1 A | 18 / 77 |
| worst objective difference | 0.0890 |

### What it settles

The consolidated fast arms move the emitted structure on 126/126 targets, worst **1.626 A** (`fd`)
and **1.898 A** (`analytic`). Those looked like a rewrite corrupting the answer. **They are the
same magnitude as the reference disagreeing with itself under an exactly null change.** The stage
is degenerate -- a CA trace admits two ideal-geometry torsion solutions at near-equal objective
distance -- and L-BFGS-B starts from a fully extended or fully helical chain, far from any
minimum, where early steps are large enough to carry a 1e-13 A difference into the other branch.
The fast arm's objective is *lower* on 66 targets and higher on 60: it lands elsewhere, it does
not converge worse.

**So "does the optimised projection agree with the reference?" has no answer at 1 A resolution,
because the reference is not a function at that resolution.** The only defensible standard is
bit-exactness, and that is what shipped: the `exact` arm reproduces `s8.project` at **0.000e+00 on
all 126 targets** for 2.19x, won purely by removing work no value depends on (`numpy.cross`'s
`moveaxis` argument handling was 60% of the builder at 1,589,366 calls; the carbonyl oxygen and CB
are built and never read). The 16.4x `analytic` arm is retained, labelled, excluded from the
headline, and **in the cache key** -- it was previously an import-time environment global, so two
runs under different gradient modes produced the same `cfg_key` and the second served the first's
arrays, yielding a directory of mixed provenance.

### The methodological point

This control exists because an audit caught a claim that the fast path "changes only how fast one
L-BFGS-B iteration is computed". The gradient and the builder had each been validated **in
isolation** -- correctly -- and the emitted structures had never been compared. *Validating the
parts is not validating the whole*, and for a degenerate objective the gap between them is
angstroms, not ulps.


<!-- ==================================================================== -->
<!-- APPENDIX: the agent-level working papers -->
<!-- ==================================================================== -->

# Appendix -- the agent-level working papers

Ten working papers, merged verbatim, that the Sprint 8 and Sprint 9 coordinator records above
were written from. They are kept because they carry per-arm detail, per-target tables and
method notes that the coordinator records cite but do not reproduce.

**The coordinator record governs.** Where a working paper's framing differs from the record
above, the record above is the adjudicated version -- most consequentially for `S9 / loop`,
whose "0.49 A per A of query" multiplier reading is reframed by S9-4 as a *contraction* that
never improves on its input, and for `S9 / refine`, whose positive-phi defect S9-9 relocates
from the consensus to the projection. Read the corrections ledger at the top of this document
before quoting anything from this appendix.


<!-- -------------------------------------------------------------- -->

## Appendix: S8 / consensus2

> *Working paper, merged verbatim from `s8/consensus2_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S8-9. How far filter-then-consensus goes, and what stops it

`s8/consensus2.py` (7 stages), `s8/test_consensus2.py` (32 tests, 0 failures),
`s8/consensus2_{sweep,fit,multi,refine,dev,leak}.json`.  The 126-target tuning instrument,
K=500 BLOSUM pools, one pass over every arm on byte-identical candidates.

**Instrument validity is asserted as a test before anything else runs**: the shipped score
reproduces **3.4540**, pool best **1.7108**, and S8-8's headline arm `medoid75`
**3.2822** - the three numbers this study stands on, to 4 decimal places.

Sprint 8's only surviving intervention was: filter the shipped BLOSUM top-500 to the
score's top 75, return the consensus **medoid** of that set (3.282 vs 3.454, −0.172
[−0.316, −0.027]).  This study asks how far that goes.  **It goes about twice as far as
S8-8 got, and the thing that carries it is not a better selection rule - it is giving up
on selection.**

## The one-line result

| arm | selected | ±SE | median | <2 A | <1.5 A | d vs score | 95% CI | W/L | d vs S8-8 medoid |
|---|---|---|---|---|---|---|---|---|---|
| SHIPPED score | 3.454 | 0.147 | 3.478 | 0.214 | 0.159 | - | - | - | +0.172 |
| S8-8 `sc\|75:medoid` | 3.282 | 0.162 | 3.095 | 0.278 | 0.198 | −0.172 | [−0.316, −0.027] | 74/47 | - |
| **`sc\|75:fit`** (deployable, physical) | **3.204** | 0.152 | 2.966 | 0.278 | 0.190 | **−0.250** | **[−0.383, −0.117]** | 80/46 | **−0.078 [−0.137, −0.019]** |
| *`sc\|75:avg`* (NOT a peptide) | *3.048* | *0.147* | *2.837* | *0.294* | *0.190* | *−0.406* | *[−0.536, −0.275]* | *95/31* | *−0.234* |
| leave-fold-out over all 1120 arms | 3.104 | 0.150 | 2.880 | 0.294 | 0.190 | −0.350 | [−0.484, −0.216] | 92/34 | −0.178 |
| leave-fold-out, SELECTION operators only | 3.294 | 0.158 | 3.077 | 0.254 | 0.190 | −0.160 | [−0.298, −0.023] | 77/45 | +0.011 |
| pool best (ORACLE) | 1.711 | | | | | | | | |

`sc|75:fit` is the coordinate consensus **projected back onto the manifold of
ideal-geometry chains**.  Its CA-CA bond is exactly 3.804 A (sd 9.4e-16); the raw average's
is **2.961 A** against a native 3.812, so the raw average is not a peptide and about
**37% of its apparent −0.406 A is geometry rather than accuracy**.

## 1. The ceiling of filter-then-consensus

For every filtered set: the ORACLE best inside it, the subpool mean (what a random draw
returns), and what consensus actually returns.  Shipped-score filter, 126 targets:

| m | oracle in the filter | subpool mean | medoid | consensus headroom |
|---|---|---|---|---|
| 10 | 2.834 | 3.478 | 3.362 | 0.528 |
| 25 | 2.609 | 3.502 | 3.369 | 0.761 |
| 50 | 2.403 | 3.539 | 3.340 | 0.936 |
| **75** | **2.306** | **3.551** | **3.282** | **0.976** |
| 100 | 2.210 | 3.567 | 3.330 | 1.119 |
| 150 | 2.106 | 3.606 | 3.383 | 1.278 |
| 250 | 1.946 | 3.759 | 3.488 | 1.542 |
| 500 (no filter) | 1.711 | 4.453 | 3.706 | 1.995 |

**The filter destroys achievable accuracy faster than consensus recovers it.**  Tightening
from 500 to 75 costs **0.595 A of ceiling** (1.711 → 2.306) and buys **0.424 A of
consensus** (3.706 → 3.282).  Consensus never gets back what the filter throws away, and
at the operating point ~1.0 A of headroom is left inside a set the method already holds.
That 0.976 A, not the 1.711 A pool best, is what a perfect *in-filter* ranker would be
worth - and it is the same order as the 1.673 A selector constant S8-5 fitted.

**The sharpest single number in the study.**  A RANDOM 75 of the same pool has oracle best
**2.096** - *better* than the shipped score's top 75 at **2.306** - while its subpool mean
is 4.445 against the score's 3.551.  **The score's filter is worse than random at keeping
the pool's best candidate and much better than random at concentrating the subpool**, and
consensus over it returns 3.048 against a random 75's 3.400.  The filter is not working by
retaining good candidates.  It works by collapsing the subpool onto one structural mode,
which is the only condition under which a consensus means anything.  This is S8-5's
"achievable keys destroy the near-native band" mechanism, and here that destruction is the
feature rather than the bug.

## 2. The joint sweep: 7 filter keys x 8 sizes x 20 operators = 1120 arms

The optimum is joint but very flat.  Best cells for the two headline operators:

| filter key | oracle @75 | medoid@75 | avg@75 | medoid@100 | avg@100 |
|---|---|---|---|---|---|
| `sc` shipped score (incumbent) | 2.306 | **3.282** | **3.048** | 3.330 | 3.056 |
| `ref1` CA-RMSD to the score's top-1 (S8-3's form) | 2.405 | 3.311 | 3.068 | **3.265** | **3.043** |
| `refmed` CA-RMSD to the medoid of the score's top-75 | 2.573 | 3.331 | 3.137 | 3.319 | 3.126 |
| `scty2` score + 0.5 x typicality, rank sum | 2.417 | 3.435 | 3.171 | 3.425 | 3.152 |
| `scty` score + typicality, rank sum | 2.602 | 3.608 | 3.340 | 3.595 | 3.330 |
| `ty` typicality alone | 2.759 | 3.910 | 3.672 | 3.933 | 3.701 |
| `rand` control | 2.096 | 3.706 | 3.400 | 3.709 | 3.389 |

* **Filtering by typicality is a disaster** (3.91).  Typicality is a fine *ranker* inside a
  band and a terrible *filter*, which is the exact mirror of the shipped score.  Fusing the
  two also loses: every `scty*` cell is worse than the score alone.  The composition S8-8
  found is not improved by making the filter smarter.
* **S8-3's "rank by RMSD to a predicted structure" form works, and it is a consensus in
  disguise.**  `ref1` (reference = the score's own argmin) ties or slightly beats the score
  as a filter at m=100.  It carries no new information - the reference is chosen by the
  score - so its gain is purely that RMSD-to-a-reference is a better *neighbourhood*
  definition than score rank.
* **Size**: the medoid family has an interior optimum at m=75 (3.282) that falls away in
  both directions; the average family is flat from m=25 to m=250 (3.043-3.133) and only
  degrades at m=500.  Constructing is far less sensitive to the filter size than selecting.
* The `rand` control at m=500 is the whole pool: consensus over the unfiltered pool is
  3.396-3.706 depending on operator.  **The filter is worth 0.35-0.42 A and is essential**,
  reconfirming S8-8 on a much larger arm set.

## 3. The consensus operator: every robust variant of the medoid is null; only CONSTRUCTING moves

At the incumbent filter (score top-75), 126 paired targets:

| operator | selected | d vs score | 95% CI | W/L | d vs S8-8 medoid |
|---|---|---|---|---|---|
| *`avg`* coordinate average (NOT a chain) | *3.048* | *−0.406* | *[−0.536, −0.275]* | *95/31* | *−0.234* |
| *`avg_iter`* Frechet/Karcher mean (NOT a chain) | *3.052* | *−0.402* | *[−0.530, −0.273]* | *95/31* | *−0.230* |
| *`avg_trim`* average of the central 75% | *3.079* | *−0.375* | *[−0.516, −0.234]* | *90/36* | *−0.203* |
| *`avg_wty`* typicality-weighted average | *3.122* | *−0.332* | *[−0.481, −0.183]* | *86/40* | *−0.160* |
| **`fit`** the projection of `avg` onto ideal geometry | **3.200** | **−0.254** | **[−0.386, −0.121]** | 82/44 | **−0.082** |
| `medoid_sq` mean SQUARED RMSD (true Frechet medoid) | 3.260 | −0.194 | [−0.341, −0.048] | 77/43 | −0.023 |
| **`medoid` (S8-8 incumbent)** | **3.282** | **−0.172** | **[−0.316, −0.027]** | 74/47 | - |
| `trim25` trimmed medoid (drop the outer 25%) | 3.287 | −0.167 | [−0.322, −0.013] | 73/50 | +0.005 |
| `snapavg` nearest real candidate to the average | 3.302 | −0.152 | [−0.298, −0.007] | 73/49 | +0.020 |
| `iter2` iterated medoid (re-filter to neighbours, repeat) | 3.306 | −0.148 | [−0.303, +0.007] | 69/55 | +0.024 |
| `reb_med` the medoid rebuilt from its own torsions | 3.314 | −0.140 | [−0.286, +0.007] | 72/54 | +0.032 |
| `wmed_sc` score-weighted medoid | 3.327 | −0.127 | [−0.216, −0.038] | 52/29 | +0.045 |
| `trim50` | 3.343 | −0.111 | [−0.264, +0.041] | 68/55 | +0.060 |
| `typic` distance-matrix typicality inside the filter | 3.363 | −0.091 | [−0.236, +0.054] | 63/53 | +0.081 |
| `clmed3` medoid of the largest of 3 clusters | 3.365 | −0.089 | [−0.247, +0.069] | 69/54 | +0.083 |
| `clmed15` medoid of the largest 1.5 A cluster | 3.398 | −0.056 | [−0.217, +0.105] | 66/57 | +0.116 |
| `wmed_ty` typicality-weighted medoid | 3.461 | +0.007 | [−0.165, +0.180] | 63/60 | +0.179 |
| SHIPPED score (= `score` operator) | 3.454 | - | - | - | +0.172 |
| **`tors_trim`** trimmed circular-mean torsions | **3.922** | **+0.468** | **[+0.153, +0.784]** | 59/67 | +0.640 |
| **`tors`** circular-mean torsions | **4.072** | **+0.618** | **[+0.265, +0.972]** | 57/69 | +0.790 |
| `avg_scale` the average rescaled to a 3.80 A step | 4.086 | +0.632 | [+0.361, +0.903] | 52/74 | +0.804 |

**Everything the brief proposed as a better SELECTION rule is null.**  Weighted medoids
(both weightings), trimmed medoids (both fractions), the iterated medoid, and the
cluster-restricted medoids are all inside or worse than the incumbent's interval.  The one
selection variant that is even nominally better - `medoid_sq`, the mean *squared* RMSD -
buys 0.023 A.  Leave-fold-out over all 12 selection operators x 8 sizes returns **3.294**,
i.e. **the whole family of ways to pick one real candidate by consensus is at its ceiling
and S8-8 already found it.**

**Torsion-space consensus is much worse, not better, and the brief's reasoning for it was
sound but the premise was wrong.**  The circular mean of the torsions (computed as atan2 of
summed unit vectors - `t_circmean_wraparound` pins the 170/−170 case) gives 4.072 A,
**+0.618 [+0.265, +0.972] worse than the shipped score**.  Torsion error compounds along
the chain: averaging the torsions of two structures that agree at the middle swings their
ends apart, so a torsion-space mean is not a structural mean at all.  Exact bond geometry
does not compensate.

**What does work is abandoning selection entirely** and returning a structure no pool
member equals.  That structure sits **1.131 A from the nearest real candidate**, so it is
genuinely new; and `snapavg` - the real candidate closest to it - gives back 3.302, i.e.
**the entire gain is in the construction, not in it being a better pointer**.

## 4. The geometry caveat is now half the result

| arm at `sc\|75` | selected | mean CA-CA | sd | rg / rg_native | distance to the nearest real candidate |
|---|---|---|---|---|---|
| native | - | 3.812 | - | 1.000 | - |
| `medoid` (a real window) | 3.282 | 3.812 | - | 1.018 | 0 |
| *`avg`* | *3.048* | **2.961** | 0.311 | 0.967 | 1.131 |
| *`avg_iter`* | *3.052* | 2.944 | 0.308 | 0.969 | 1.190 |
| `avg_scale` (rescaled) | 4.086 | 3.800 | 0.456 | 1.308 | 1.860 |
| **`fit`** (projected) | **3.200** | **3.804** | **9.4e-16** | - | - |
| `tors` | 4.072 | 3.804 | 0.000 | 1.182 | 1.682 |
| `reb_med` | 3.314 | 3.804 | 0.000 | 1.023 | 0.321 |

**The coordinate average contracts the chain by 22%.**  This is not a bug in the code; it
is what averaging does to a curve - the mean of points scattered about a path lies inside
it.  A shrunken, smoothed object flatters CA-RMSD, so `avg`'s −0.406 A cannot be reported
as an accuracy gain.

Two ways to fix it, and only one works:

* **Rescaling** (`avg_scale`) restores the 3.80 A step by a global scale and inflates the
  radius of gyration to 1.308x native.  It is 1.04 A **worse**.  The contraction is not
  uniform, so a uniform correction is the wrong inverse.
* **Projection** (`fit`): minimise CA-RMSD to the average over the torsions of an
  ideal-geometry chain, i.e. find the nearest thing to the consensus that a real peptide
  could be.  This gives **3.200 A** with exact bond geometry, so **−0.254 A of the −0.406
  survives and −0.152 A was geometry.**

The projection lands **0.791 A** from the average, which is itself a measurement: the
coordinate consensus is 0.79 A off the manifold of physical backbones.

**The projection needs no torsions.**  Started from four fixed generic conformations
(extended, alpha, beta, PPII) rather than the medoid's own torsions it gives **3.204 vs
3.200** on the full instrument, so it runs on any pool of bare CA traces - which is what
made the dev pass possible, since the dev pools store coordinates and no torsions.

### The rebuild displacement, measured

S8-8 carried a 0.399 A rebuild-displacement caveat.  Measured here on the same pools:

* the ideal-geometry rebuild of a candidate sits **0.3985 A** CA-RMSD from the real window
  (reproducing S8-8's figure exactly), **but**
* its CA-RMSD **to the native** is only **+0.0088 A** worse, averaged over all 500x126
  candidates, and **+0.032 A** on the selected medoid specifically (`reb_med` 3.314 vs
  `medoid` 3.282, CI [−0.286, +0.007] against the score).

**So the 0.399 A displacement is real as a displacement and almost free as an error.** The
rebuild moves a structure sideways within its own basin rather than away from the native.
That retires the caveat for every Legacy and geometry row in S8-8: those numbers are not
inflated by 0.4 A, they are inflated by about 0.01 A.

## 5. Multi-hypothesis: 0.43 A of oracle headroom, and no native-free rule reaches any of it

Average-linkage clustering of the score's top 75, consensus per cluster
(`s8/consensus2_multi.json`):

| k | mean clusters | largest cluster's share | best-of-k (ORACLE) | mean-of-k | worst-of-k | pick largest | pick by score | pick tightest |
|---|---|---|---|---|---|---|---|---|
| 1 | 1.00 | 1.000 | 3.282 | 3.282 | 3.282 | 3.282 | 3.282 | 3.282 |
| 2 | 2.00 | 0.861 | **3.003** | 3.572 | 4.141 | 3.371 | 3.480 | 3.599 |
| 3 | 3.00 | 0.762 | **2.856** | 3.646 | 4.346 | 3.365 | 3.567 | 3.691 |
| 5 | 5.00 | 0.643 | **2.752** | 3.695 | 4.520 | 3.350 | 3.591 | 3.771 |

and with the coordinate average as the per-cluster consensus: best-of-2 **2.791**,
best-of-3 **2.679**, best-of-5 **2.587**, against a single average of 3.048.

**S8-2's multimodality is real and returning k answers is worth 0.28-0.53 A - to whoever
can pick among them.  Nothing native-free can.**  All three rules tested are worse than not
splitting at all: picking the largest cluster is +0.089 [+0.029, +0.149] at k=2 and +0.083
[+0.008, +0.157] at k=3; picking the lowest-mean-score cluster is +0.198 and +0.285;
picking the tightest cluster is +0.317 and +0.409.  Every one of those intervals excludes
zero **in the wrong direction**.

The reason is visible in the third column: the largest cluster already holds 86% of the
mass at k=2 and 76% at k=3, so "pick the largest" almost always returns the single-medoid
answer, and the +0.089 is the price of the cases where it does not.  The alternative
hypotheses are small, and being small is exactly what makes them unpickable.

This is the same shape as every other recognition result in sprints 7 and 8: the
information is present in the pool and no native-free channel reads it.

## 6. Does refinement stack?

No, in either form, and the energy falls every time.  `legacy_refine.refine_angles` is a
bounded pattern search in torsion space; started from the consensus output rather than from
the oracle and native starts sprint 6 used (`work/refine_study.py`):

| start | 0 steps | box 5 deg | d vs start | 95% CI | W/L | dE | box 15 deg | d vs start |
|---|---|---|---|---|---|---|---|---|
| the medoid's own rebuild | 3.314 | 3.376 | +0.062 | [-0.005, +0.129] | 58/68 | -7.6 kcal | 3.512 | +0.198 [+0.077, +0.319] |
| the torsion projection `fit` | 3.204 | 3.270 | +0.066 | [+0.004, +0.128] | 58/68 | -6.2 kcal | 3.379 | +0.175 [+0.058, +0.292] |

**Legacy relaxation is neutral-to-harmful from both starts, and the harm grows with the
box.**  The energy falls by 6 to 10 kcal in every cell, so this is the objective's minimum
being in the wrong place, not an optimiser failure -- the same conclusion sprint 6 reached
from oracle and native starts, now measured from the start the pipeline would actually
ship.  A refinement operator whose minima are 1-2 A from the native basin cannot repair a
3.2 A answer.

**Amber was not duplicated.**  The concurrent relaxation study's committed `relax_sweep.json`
shows restrained minimisation at k=100 and k=10 moving a candidate by 0.011 and 0.0003 A of
CA-RMSD at 4-7 s per candidate; there is nothing there for consensus to stack with at this
scale, and the box is shared.

## 7. The single dev pass

One pre-registered arm, run once, no iteration.  `DEV_ARM = ("sc", 75, "fit")` was pinned in
the module and asserted by a test before `stage_dev` was ever executed; the choice is the
leave-fold-out cell (`sc|75`, picked in 3 of 5 folds, the synthesis family in all 5) with
the raw average replaced by its physically-valid projection, which costs 0.156 A of tuning
RMSD and is not negotiable.  The pools are `s8/inband_devpool`, built and asserted against
`s7/audit_cache` by `s8.inband` -- this module does not rebuild or alter them.  Dev pool
best is 1.539 A.

| arm | selected | +-SE | median | <2 A | <1.5 A | d vs score | 95% CI | W/L |
|---|---|---|---|---|---|---|---|---|
| SHIPPED score | 3.475 | 0.354 | 3.526 | 0.208 | 0.167 | - | - | - |
| `sc\|75:medoid` (S8-8's arm, reproduced) | 3.316 | 0.373 | 3.183 | 0.250 | 0.125 | -0.159 | [-0.456, +0.139] | 14/10 |
| **`sc\|75:fit` (pre-registered)** | **3.221** | 0.349 | 3.048 | 0.250 | 0.125 | **-0.255** | **[-0.501, -0.008]** | **16/8** |

`fit` vs `medoid` on dev: -0.096 [-0.236, +0.044], 19/5.

**Tuning said -0.250, dev says -0.255, and that agreement is worth more than the interval.**
The interval does exclude zero and it should NOT be reported as confirmation: n=24 at
SE 0.349 has no power to resolve a 0.25 A effect, so an interval that just clears zero is
noise aligning with the truth rather than evidence.  The tuning instrument is the evidence.
What dev establishes is that the direction and magnitude replicate on a cluster-disjoint set
that no arm selection ever touched, and that S8-8's medoid arm reproduces its own -0.159
exactly.

## 8. What bounds this architecture, quantified two ways

### (a) The consensus criterion contains no nativeness information - ORACLE DIAGNOSTIC

Drop the NATIVE into the score-filtered top 75 and rank the augmented set by the very
criterion each operator minimises.  (The augmented pairwise matrix needs no new work: the
native's row is exactly `rr`, which the cache already holds.)

| criterion, inside the score's top 75 | native's mean percentile | median | +-SE | native below the 10th pct | native is the argmin |
|---|---|---|---|---|---|
| `medoid` (mean RMSD to the others) | **81.1** | 96.7 | 2.5 | 0.048 | 2/126 |
| `medoid_sq` | 80.5 | 97.3 | 2.6 | 0.048 | 2/126 |
| `typic` (distance-matrix typicality) | 80.2 | 95.3 | 2.6 | 0.063 | 3/126 |
| `d_to_avg` (distance to the SYNTHESISED average) | 80.9 | 97.3 | 2.5 | 0.048 | 1/126 |
| `d_to_fit` (distance to the projected consensus) | 80.5 | 96.0 | 2.5 | 0.048 | 2/126 |
| `d_to_avgtrim` (distance to the trimmed average) | **79.9** | 94.7 | 2.5 | 0.040 | 2/126 |
| the shipped score, inside its own top 75 | 79.0 | 100.0 | 3.1 | 0.079 | 3/126 |

This reproduces the integration agent's measurement (82.8 mean / 97.3 median) on this
instrument at 81.1 / 96.7, and answers the question it raises: **no variant of the criterion
places the native lower.**  The whole spread is 79.9 to 81.6 on an SE of 2.5 - squaring the
distance, using the distance matrix instead of superposition, trimming the outliers, and
measuring against the synthesised structure instead of against the other members are all the
same measurement.  The consensus criterion measures centrality and nothing else.

The same table on the UNFILTERED pool and on a random 75 locates the effect exactly:

| set | `medoid` pct | `typic` pct | shipped `score` pct | subpool spread |
|---|---|---|---|---|
| the whole BLOSUM top-500 | 61.5 | 53.1 | **36.8** | 4.247 A |
| a RANDOM 75 of it | 60.6 | 52.2 | **37.2** | 4.263 A |
| the score's top 75 | **81.1** | 80.2 | **79.0** | 2.486 A |
| `ref1`'s top 100 | 81.6 | 75.8 | 66.6 | 2.480 A |

**Two things are visible here and both matter.**  First, the native is *already* a mild
outlier of the library's own distribution - 61.5th percentile of centrality in the raw pool,
where chance is 50 - so consensus starts out slightly anti-correlated with nativeness before
any filtering happens.  Second, the shipped score carries real nativeness information on the
raw pool (36.8th percentile; the native is below the 10th on 26% of targets) and **loses all
of it inside its own top 75** (79.0, median 100.0).  That is S8-8's "good coarse filter, bad
fine ranker" stated at the native rather than in aggregate, and it is the sharpest form of it
the project has.

**The consequence is a hard bound on member selection.**  If the criterion ranks the native
at the 81st percentile of the filtered set, then no better medoid rule, no better filter size
and no reweighting can find the native - what such an operator converges to is the *mode* of
the filtered set.  At `sc|75` that mode is **3.282 A** while the best MEMBER is **2.306 A**,
and the 0.976 A between them is unreachable by any rule that returns a member.  Every null in
section 3 is that bound being hit from a different direction.

### (b) Synthesis is not bounded that way, and it has its own ceiling: 85% of the error is shared

Synthesis returns a structure no member equals, so it is not capped by the mode.  What caps
it instead is how much of the error is *common* to the filtered set.  Decomposing the mean
square at `sc|75`:

| | mean | RMS |
|---|---|---|
| `medoid` (a member) | 3.282 | 3.748 |
| `avg` (synthesis) | 3.048 | 3.462 |
| `fit` (physical synthesis) | 3.200 | 3.634 |
| best member (ORACLE) | 2.306 | 2.679 |

sqrt(3.748^2 - 3.462^2) = **1.437 A** of independent, averageable error, against a **3.462 A**
component shared by every candidate - **85.3% of the medoid's mean square is bias that
averaging cannot touch.**

That figure predicts everything else in the study.  It predicts that the average is flat in
the filter size (3.146 at m=10, 3.048 at 75, 3.056 at 100, 3.072 at 150 - averaging 75
samples has already exhausted the independent component, so more of them buy nothing).  It
predicts that trimming, weighting and iterating the average change it by 0.03-0.07 A.  And it
predicts that synthesis travels **1.145 A away from the mode** to move only **0.234 A closer
to the native**: the direction it goes is variance reduction, not nativeness.

**So the architecture has two ceilings and they are close together.**  Member selection is
capped at the mode, 3.282 A.  Synthesis is capped by the shared bias at ~3.05 A unphysical
and ~3.20 A physical.  The filtered set's own best member, 2.306 A, is below both and
reachable by neither.

## 9. The direct answer: how far this goes, and what stops it

**How far it goes.**  From 3.454 A to **3.204 A** - a deployable, leakage-clean, physically
valid structure whose CA-CA bond is exact to 1e-15 - at -0.250 [-0.383, -0.117] on the 126
tuning targets and -0.255 [-0.501, -0.008] on the 24 dev targets, replicating to 0.005 A.
That is 1.45x S8-8's gain and it reduces S8-5's 1.673 A selector constant by about 15%.
Fraction below 2.0 A goes 0.214 -> 0.278; median 3.478 -> 2.966.

**What stops it, in order of how much it costs:**

1. **The filtered set is biased, not noisy.**  85% of the consensus error is common to every
   candidate the filter kept, and averaging removes only the other 15%.  This is the binding
   constraint and it is a property of the FILTER, not of the operator.
2. **The criterion carries no nativeness information.**  The native is at the 81st percentile
   of every consensus criterion tested, so member selection cannot beat the mode (3.282) even
   though the set contains a 2.306 A member.  Six variants move that percentile by 1.2 points
   on an SE of 2.5.
3. **The filter destroys more ceiling than consensus recovers.**  500 -> 75 costs 0.595 A of
   achievable accuracy and buys 0.424 A of consensus.  Yet the filter cannot be dropped -
   consensus over the unfiltered pool is 0.35-0.42 A worse - because concentrating the
   subpool onto one mode is what makes a consensus mean anything.  A random 75 has a BETTER
   best member (2.096 vs 2.306) and a much worse consensus (3.400 vs 3.048).
4. **Nothing native-free chooses among hypotheses.**  Multi-hypothesis output is worth
   0.28-0.53 A to an oracle and is negative under every rule tested.
5. **Refinement does not stack.**  Legacy relaxation is +0.06 to +0.20 A while lowering the
   energy by 6-10 kcal.

**What would move it.**  Only two things here are not closed.  The first is reducing the
SHARED bias of the filtered set, which is a generation problem: S8-5 and S8-6 both closed it
for the current library, so it means a different library rather than a different filter.  The
second is a criterion that ranks the native below the 81st percentile of its own filtered
set; every native-free signal tried in sprints 7 and 8, including all six variants here,
fails that test, and one that passed it would be worth far more than the 0.25 A this study
bought.

**What this study closes**: the consensus operator (every robust, weighted, iterated,
cluster-restricted and torsion-space variant), the consensus filter key (7 keys x 8 sizes;
the shipped score and RMSD-to-its-own-argmin tie and nothing beats them), multi-hypothesis
selection, refinement stacking, and the rebuild-displacement caveat - which is 0.399 A of
displacement and 0.009-0.032 A of error.

## Leakage

`stage_leak` NaN-poisons `rr` and `nat_ca` in the cache and recomputes **the filter orders
themselves** as well as every operator - so a filter key that had read a native would
change which candidates are even considered, a stronger test than fixing the index set
first.  **12 targets x 7 filters x 8 sizes x 20 operators = 1120 arms per target, worst
absolute difference 0.000e+00.**  `test_consensus2.py` repeats it independently on 3
targets.  Natives enter only as `rr`/`nat_ca` reporting labels, in the columns named
`oracle`/`best_of_k`, and in nothing else.

## Method notes kept on the record

* **A cell-key collision was found and fixed.**  Cells were keyed `f"{filter}{m}"`, so the
  `scty` filter at m=250 and the `scty2` filter at m=50 both produced `"scty250"` and one
  silently overwrote the other.  Keys now carry a separator and `t_cell_keys_unique` asserts
  56 distinct keys where the un-separated form gives 55.  Neither filter was a winner, but
  the first version of the table reported one of them under the other's numbers.
* **The tie trap.**  Every selection operator is scored through `s8.inband.sel_of`, which
  averages the true RMSD over the whole tied argmin set rather than taking `np.argmin`'s
  first index.  `t_ops_tie_safe` permutes the candidate order and asserts all 12 selection
  operators are invariant to 0.0e+00.
* **Leave-fold-out alignment.**  `t_lfo_alignment` re-derives every leave-fold-out mean and
  paired difference from the per-target rows and the recorded per-fold picks, because the
  misalignment bug that inflated an S8-8 CI by 2x on an identical mean is invisible in the
  mean.
* **Circular means.**  `t_circmean_wraparound` asserts the mean of +170 and −170 degrees is
  180, and separately asserts the naive arithmetic mean would have returned 0.


<!-- -------------------------------------------------------------- -->

## Appendix: S8 / enstarget

> *Working paper, merged verbatim from `s8/enstarget_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S8-14. The training target IS misspecified, by 0.5 A per pair, and correcting it is worth 0.04 A

`s8/enstarget.py` (8 resumable stages), `s8/test_enstarget.py` (44 checks, 0 failures),
`s8/enstarget_{noise,eval,filt,dev,leak}.json`, models in `s8/enstarget_models/`,
ensemble label cache `s8/enstarget_ens.npz`.

126-target tuning instrument, K=500, the same cached window universes
(`s8/generate_univ/*.npz`) every sprint-7/8 number reports on. Instrument validity is
asserted as tests: the shipped selector reproduces **3.454 A**, pool best **1.711**, and
S8-8's `medoid75` architecture reproduces **3.282 A** on this module's own code path.
`stage_leak` NaN-poisons `rr` and `nat_ca` on 160 target x arm combinations and asserts the
deployable score is bit-identical - worst |diff| **0.000e+00**.

## The question, and what makes it distinct

Every learned component in this project is trained to predict properties of **PDB model 1**.
115/126 of the instrument is solution NMR and 111/126 are multi-model, so "the structure"
is usually one arbitrary member of a ~1 A-wide conformational ensemble, and a model fitted
to model 1's exact distance matrix is being asked to reproduce conformational noise as
signal.

**This is not the claim S8-4 already tested.** S8-4 measured that scoring pools against a
RANDOM deposited model instead of model 1 changes pool best by +0.001 A [-0.030, +0.033] -
model 1 is not a biased choice of EVALUATION reference, and nothing here disputes that.
The claim tested here is about TRAINING: fitting model 1 exactly forces the model to
memorise the part of model 1 its own ensemble does not reproduce. Both can be true and
both are.

**Legitimacy.** Using the deposited ensemble of a TRAINING structure is ordinary
supervision - the same label file, read differently. No held-out target's native is
touched in any form at inference. Fold discipline is unchanged (`peptide_db.folds(5)` over
identity clusters plus `distogram._fold_fragments`); the identity audit over all five
folds returns a worst train-vs-heldout identity of **0.588** against the 0.600 threshold.

## 1. The irreducible noise floor - the cheap question, answered first

Per-pair CA-CA distances over every deposited model of every training entry. Distances,
not coordinates: a distance matrix is superposition-free, so rigid-body differences
between deposited models contribute nothing and what survives is exactly the
conformational disagreement the label question is about.

**698/787 training entries are multi-model (median 20 models), and 90.3% of training PAIRS
come from one.**

| quantity | value |
|---|---|
| mean \|d(model 1) − d(ensemble mean)\| | **0.497 A** |
| mean per-pair SD across models | **0.635 A** |
| p90 / max of the model-1 deviation | 1.283 / 13.154 A |
| pairs with ensemble SD > 0.5 A / > 1.0 A | 34.2% / 19.1% |

| separation shell | pairs | noise | ens SD | signal SD | noise/signal |
|---|---|---|---|---|---|
| 2-3 | 22,885 | 0.158 | 0.198 | 1.471 | 0.135 |
| 4-5 | 19,737 | 0.297 | 0.387 | 2.731 | 0.142 |
| 6-8 | 23,703 | 0.402 | 0.522 | 3.514 | 0.148 |
| 9-13 | 24,590 | 0.580 | 0.746 | 4.761 | 0.157 |
| **14+** | 12,961 | **1.028** | **1.291** | 7.510 | **0.172** |

The hypothesis survives its own cheap test: the label carries **half an Angstrom** of
conformational noise per pair on average and **a full Angstrom** in the long-range band,
which is the band the whole failure is in. Noise-to-signal is a near-constant 0.135-0.172
per shell.

## 2. But the ceiling was computable before anything was trained

A model that reproduced its noisy labels exactly would carry their noise in quadrature
with its own error, because the conformational displacement of one deposited model from
its ensemble mean is uncorrelated with what a sequence-only predictor gets wrong. With
`MAE = sigma * sqrt(2/pi)`:

| | value |
|---|---|
| label noise SD | 0.635 A |
| model-1 distance MAE of the matched baseline | 2.273 A → model error SD 2.849 A |
| **noise-to-error ratio** | **0.223** |
| MAE floor with the label PERFECTLY de-noised | 2.222 A |
| **maximum recoverable distance MAE** | **0.052 A** |
| …and the ensemble MEAN of ~20 models removes ~78% of the noise, so | **0.040 A expected** |

**The noise is 22% of the error and 5% of the loss.** That is the whole arithmetic of the
axis: a term that enters in quadrature at a fifth of the error is a second-order term, and
no reading of the label can make it first-order.

## 3. The arms - 126 targets, K=500, matched everything

Five arms sharing feature set, architecture, optimiser, seed, folds and corpus; only the
LABEL differs. `m1` reproduces `distogram.MLP.fit` **bit-for-bit** (asserted as a test), so
every delta is attributable to the target.

**Feature-set note, stated up front.** The shipped distogram uses ESM-2 features;
re-training five arms on them needs the 1.5 GB `esm_cache.npz` resident, and the box ran
at 83-96% throughout with three other agents on it. All arms are trained on the
sequence-only feature set, which is a MATCHED comparison - inter-arm deltas are fully
attributable - but means the absolute numbers are not the shipped predictor's. The shipped
ESM score is carried through every table as the external reference so the handicap is
visible.

`pep` = out-of-fold peptides only, **90% of pairs carry an ensemble**: where the
intervention is testable. `frag` = the SHIPPED corpus, peptides plus
`distogram._fold_fragments`; fragments are windows of `prots/` crystal structures, every
one of them single-model, and they are 82% of the training pairs, so only **13.4%** of the
shipped corpus's pairs carry an ensemble at all.

| arm | selected | ±SE | median | <2 A | <1.5 A | MAE\|m1 | MAE\|ens | native pct | d sel vs matched m1 | 95% CI | W/L |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SHIPPED (ESM, external ref) | 3.454 | 0.147 | 3.478 | 0.214 | 0.159 | 2.339 | 2.237 | 36.8 | - | - | - |
| pep:m1noise *(control)* | 3.537 | 0.158 | 3.545 | 0.246 | 0.159 | 2.260 | 2.030 | 37.5 | **-0.063** | [-0.148, +0.021] | 36/35 |
| **pep:ensmean** | **3.556** | 0.155 | 3.545 | 0.262 | 0.159 | **2.241** | **2.008** | 37.7 | -0.045 | [-0.172, +0.082] | 37/41 |
| **pep:m1** *(matched baseline)* | **3.601** | 0.158 | 3.575 | 0.246 | 0.159 | **2.273** | **2.047** | 36.6 | - | - | - |
| pep:confw | 3.620 | 0.155 | 3.633 | 0.238 | 0.143 | 2.293 | 2.078 | 38.1 | +0.019 | [-0.106, +0.145] | 40/48 |
| pep:ensdist | 3.641 | 0.164 | 3.805 | 0.254 | 0.143 | 2.245 | 2.011 | 38.0 | +0.040 | [-0.069, +0.150] | 38/32 |
| frag:ensmean | 3.723 | 0.156 | 3.527 | 0.190 | 0.103 | 2.569 | 2.458 | 42.7 | -0.061 | [-0.179, +0.058] | 47/48 |
| frag:confw | 3.746 | 0.146 | 3.596 | 0.159 | 0.103 | 2.619 | 2.498 | 43.1 | -0.038 | [-0.169, +0.094] | 63/39 |
| **frag:m1** *(shipped corpus baseline)* | **3.784** | 0.151 | 3.676 | 0.151 | 0.087 | 2.548 | 2.437 | 43.1 | - | - | - |
| frag:ensdist | 3.849 | 0.156 | 3.708 | 0.167 | 0.111 | 2.582 | 2.475 | 42.5 | +0.065 | [-0.038, +0.169] | 40/47 |

**Distance MAE, against both references, paired against the same corpus's m1:**

| arm | d MAE\|model 1 | 95% CI | d MAE\|ensemble | 95% CI | d MAE\|m1, long-range only | 95% CI |
|---|---|---|---|---|---|---|
| pep:ensmean | **-0.032** | [-0.076, +0.012] | **-0.039** | [-0.083, +0.005] | -0.252 | [-0.517, +0.013] |
| pep:ensdist | -0.028 | [-0.070, +0.013] | -0.036 | [-0.074, +0.002] | **-0.184** | **[-0.355, -0.012]** |
| pep:m1noise *(control)* | -0.013 | [-0.041, +0.016] | -0.017 | [-0.048, +0.014] | -0.047 | [-0.197, +0.102] |
| pep:confw | +0.019 | [-0.051, +0.090] | +0.031 | [-0.044, +0.106] | +0.013 | [-0.240, +0.266] |
| frag:ensmean | +0.020 | [-0.024, +0.064] | +0.021 | [-0.026, +0.069] | +0.035 | [-0.170, +0.239] |
| frag:ensdist | +0.034 | [-0.013, +0.081] | +0.038 | [-0.009, +0.084] | -0.026 | [-0.277, +0.226] |
| frag:confw | **+0.071** | **[+0.009, +0.133]** | +0.061 | [-0.006, +0.128] | +0.144 | [-0.122, +0.410] |

**On the peptide corpus the measured MAE gain matches the bound predicted before training:
0.032 A against 0.040 A expected, inside the CI. On the shipped corpus it has the opposite
sign.** Exactly one interval in the whole MAE table excludes zero and it is `frag:confw`
being *worse*.

### The dilution test says the effect is not about ensembles at all

The two corpora differ 7-fold in how much of the loss the intervention touches: 90.3% of
`pep`'s training pairs carry an ensemble against 13.4% of `frag`'s. A real noise-removal
effect must scale with that. It does not:

| | `pep` (90.3% ensemble coverage) | `frag` (13.4%) |
|---|---|---|
| ensmean, d selected vs m1 | -0.045 [-0.172, +0.082] | **-0.061** [-0.179, +0.058] |
| ensmean, d MAE vs m1 | -0.032 | **+0.020** |
| ensmean, the ensemble differential | +0.0020 [-0.0082, +0.0122] | -0.0002 [-0.0159, +0.0154] |

**The selection effect is the same size at one seventh of the exposure, and the distance
effect reverses sign.** That is the signature of sampling noise around zero, not of a
mechanism.

## 4. Two controls say it is not noise removal

The brief predicted that an ensemble-trained model should look **worse** against model 1
and **better** against the ensemble mean, and that if it does not, the setup is wrong. It
looks better against both, by nearly the same amount, and the setup is not wrong - the
prediction holds in-sample, not out-of-sample. Out of sample, fitting noise costs error
against *every* reference including the noisy one, so an arm that fits less noise improves
against both. What separates "learned a different target" from "fitted a bit better" is
the **differential**.

`gap = MAE|model 1 − MAE|ensemble` is how much closer to the ensemble mean a prediction
sits than to model 1. Paired against the m1 arm, this cancels the target's own difficulty,
which dominates either MAE alone:

| arm | d gap vs m1 | 95% CI | W/L |
|---|---|---|---|
| pep:ensmean | **+0.0020** | [-0.0082, +0.0122] | 53/58 |
| **pep:m1noise** *(noise ADDED)* | **+0.0022** | [-0.0079, +0.0123] | 55/55 |
| pep:ensdist | -0.0007 | [-0.0095, +0.0082] | 54/57 |
| pep:confw | -0.0108 | [-0.0269, +0.0054] | 54/57 |

**The ensemble-trained model does not sit closer to the ensemble mean than the
model-1-trained model does.** And the control that ADDS one more draw of the same measured
noise moves the differential by an identical +0.0022, improves MAE too
(-0.013 [-0.041, +0.016]) and posts the corpus's best selected RMSD (-0.063 vs m1).

Both interventions widen the effective target distribution of a soft-binned classifier.
That is what is being bought - generic target smoothing, obtainable by adding noise as
readily as by removing it - not a better-specified objective.

**A third control, on the one thing the ensemble target is uniquely supposed to teach.**
`ensdist` predicts the per-pair distribution across models, so it should learn *which
pairs are conformationally undetermined*. Spearman between each arm's predicted per-pair
spread and the held-out target's own deposited ensemble SD (ORACLE DIAGNOSTIC, feeds
nothing):

| arm | rho(predicted sd, ensemble sd) |
|---|---|
| pep:m1 *(baseline)* | 0.536 |
| pep:ensmean | 0.538 |
| **pep:ensdist** | **0.547** |
| **pep:m1noise** *(noise ADDED)* | **0.547** |
| SHIPPED (ESM) | 0.508 |

The model-1-trained distogram **already** predicts which pairs are undetermined, at
rho 0.536 - conformational uncertainty is largely a function of sequence separation and
local propensity, which it learns from model 1 alone. Training on the ensemble
distribution adds +0.011, and the arm trained on a *deliberately corrupted* label adds
exactly the same +0.011.

## 5. Ranking: nothing moves, and the native percentile moves the wrong way

| arm | rho (whole pool) | rho \| rg | in-band rho | sel in-band | native pct |
|---|---|---|---|---|---|
| SHIPPED | +0.568 | +0.506 | +0.126 | 2.657 | 36.8 |
| pep:m1 | +0.514 | +0.476 | +0.099 | 2.665 | **36.6** |
| pep:m1noise | +0.517 | +0.472 | +0.107 | 2.672 | 37.5 |
| pep:confw | +0.521 | +0.471 | +0.124 | 2.685 | 38.1 |
| pep:ensmean | +0.511 | +0.470 | **+0.087** | 2.678 | 37.7 |
| pep:ensdist | +0.509 | +0.463 | +0.091 | 2.662 | 38.0 |
| frag:m1 | +0.500 | +0.447 | +0.063 | 2.752 | 43.1 |
| frag:ensmean | +0.495 | +0.449 | +0.069 | 2.725 | 42.7 |
| frag:ensdist | +0.481 | +0.436 | +0.086 | 2.737 | 42.5 |
| frag:confw | +0.483 | +0.455 | +0.071 | 2.764 | 43.1 |

The native percentile is the sprint's most diagnostic number - the fraction of the pool
that outscores the true structure, benchmarks distogram 36.0-36.8, inverse folding 30.1,
Amber 54.3. **Ensemble training moves it from 36.6 to 37.7-38.1 on `pep` and from 43.1 to
42.5-43.1 on `frag`: 1.1 the wrong way on one corpus and 0.5 the right way on the other,
neither material.** It does not indicate a better-specified objective; if anything it
indicates a marginally smoother one. Partialling radius of gyration changes no ordering,
so this is not the compactness confound that has produced three false leads here.

**A side result worth recording, since the design isolates it.** With sequence-only
features the SHIPPED training corpus is *worse* than peptides alone: `frag:m1` 3.784 /
MAE 2.548 / native pct 43.1 against `pep:m1` 3.601 / 2.273 / 36.6. Fragments help the
ESM-featured production model (that is why they are in it) and hurt the sequence-only one.
Nothing here recommends changing the shipped corpus - the comparison is at the wrong
feature set for that - but it does say the fragment contribution is feature-dependent,
which nothing in the project had measured.

## 6. As a filter inside the validated architecture

S8-8's `medoid75` - filter the shipped top-500 to the score's top 75, return the consensus
medoid - reproduces here at **3.282 A**, and S8-9 already searched 165 filter x size cells
without beating the plain shipped score. The ensemble arm as the filter:

| filter @ q | selected | ±SE | vs incumbent (score@75) | 95% CI |
|---|---|---|---|---|
| **SHIPPED @ 75** | **3.282** | 0.162 | - | - |
| SHIPPED @ 25 | 3.369 | 0.163 | +0.087 | [-0.002, +0.176] |
| SHIPPED @ 150 | 3.383 | 0.160 | +0.101 | [-0.007, +0.209] |
| frag:ensmean @ 150 | 3.452 | 0.168 | +0.170 | [+0.007, +0.333] |
| **pep:ensmean @ 75** | **3.4526** | 0.165 | +0.170 | [-0.038, +0.379] |
| **pep:m1 @ 75** | **3.4533** | 0.168 | +0.171 | [-0.035, +0.377] |
| pep:confw @ 150 | 3.482 | 0.168 | +0.199 | [-0.001, +0.400] |
| **frag:ensmean @ 75** | **3.516** | 0.161 | +0.234 | [+0.075, +0.393] |
| **frag:m1 @ 75** | **3.535** | 0.167 | +0.253 | [+0.091, +0.414] |
| pep:ensdist @ 75 | 3.545 | 0.168 | +0.263 | [+0.067, +0.458] |

**Inside the architecture the intervention is worth 0.0007 A on `pep` and 0.019 A on
`frag`, against an incumbent it misses by 0.17-0.25 A.** Note separately that the
consensus step recovers 0.10-0.19 A for the weaker sequence-only filters (3.601 → 3.453
for `pep:m1`), which is an independent reproduction of S8-8's central mechanism on a
filter it was never tested with.

## 7. The single pre-registered dev pass - and it does not even replicate the sign

One arm, chosen leave-fold-out on the tuning instrument by a rule fixed before the pass
ran (lowest LFO mean selected CA-RMSD among the four `pep` arms; controls excluded by
name). Four of five folds picked `pep:ensmean`. A test asserts the pre-registration
excludes control arms. The 24 dev pools are `s8/inband_devpool`, read-only, and the pass
reproduces sprint 7's shipped 3.475 and S8-8's `medoid75` 3.316 exactly.

The matched `pep:m1` baseline is carried alongside; it is not a second candidate but the
control without which the pre-registered number is uninterpretable, since the `pep` arms
carry a sequence-only feature handicap that has nothing to do with the label.

| arm | selected | ±SE | median | <2 A | vs matched m1 | 95% CI | W/L |
|---|---|---|---|---|---|---|---|
| SHIPPED | 3.475 | 0.354 | 3.526 | 0.208 | +0.080 | [-0.176, +0.336] | 10/14 |
| **pep:m1** | **3.395** | 0.337 | 3.211 | 0.167 | - | - | - |
| **pep:ensmean** *(pre-registered)* | **3.466** | 0.327 | 3.401 | 0.167 | **+0.070** | [-0.118, +0.259] | 9/7 |
| SHIPPED @ 75 + medoid | 3.316 | 0.373 | 3.183 | 0.250 | +0.050 | [-0.270, +0.371] | 9/13 |
| pep:m1 @ 75 + medoid | 3.266 | 0.384 | 3.158 | 0.333 | - | - | - |
| pep:ensmean @ 75 + medoid | 3.389 | 0.405 | 3.217 | 0.333 | **+0.123** | [-0.005, +0.251] | 3/10 |

**Tuning said -0.045; dev says +0.070, and +0.123 inside the consensus architecture at
3W/10L.** With SE 0.33 against a 0.045 A effect, n=24 has no power to resolve either sign,
so this is not evidence of harm - it is evidence that there is nothing here to replicate.
Reporting the tuning direction as a result without this pass would have been wrong.

## 8. What the axis was worth, in the currency the sprint reports

Cross-target OLS of selected CA-RMSD on `dist_mae` with length and pool best held fixed
(the two confounds this project has been caught by):

**1 A of distance MAE is worth 0.855 A [0.728, 0.982] of selected CA-RMSD** for the
shipped predictor, and the slope reproduces on two independently trained ones -
0.872 [0.752, 0.992] for `pep:m1`, 0.970 [0.853, 1.086] for `frag:m1`.

Caveat that makes this an OVER-estimate: it is a cross-TARGET slope inside one predictor,
so `dist_mae` and the selection gap share a common cause - target difficulty - that length
and pool best only partly absorb. The project's complementary measurement on the other
axis (`prior-mae-prices-selected-rmsd`: across achievable PRIORS, r = 0.19) disagrees
exactly as it should. Read as an upper bound, which is the generous direction here.

So the whole axis, priced before it was fitted:

> maximum recoverable MAE **0.052 A** x price **0.855** = **0.044 A of selected CA-RMSD.**

against an instrument SE of 0.147 A on the mean and 0.065 A on a paired difference.
**The effect was below the instrument's resolution before a single model was trained**,
and the measured -0.045 [-0.172, +0.082] is exactly that number, unmeasurably.

## The direct answer

**Was the training target misspecified? Yes, and by more than anyone had checked - half an
Angstrom of conformational noise per pair, a full Angstrom in the long-range band, on 90%
of the training pairs.**

**Does fixing it matter? No, and the reason is arithmetic rather than luck.** Label noise
enters the error in quadrature at a noise-to-error ratio of 0.223, which caps a perfect
de-noising at 0.052 A of distance MAE - 2.3% - and at 0.044 A of selected CA-RMSD, below
the instrument's resolution. The measurement agrees with the cap (MAE -0.032 against
-0.040 predicted; selection -0.045 [-0.172, +0.082]). Four controls say that even that
much is not noise removal:

1. **The differential.** The statistic that would show the arm learned the ensemble target
   - MAE\|model 1 minus MAE\|ensemble, paired against the m1 arm - is +0.0020
   [-0.0082, +0.0122]. The arm sits no closer to the ensemble mean than the model-1 arm
   does.
2. **The noise-injection arm.** ADDING one more draw of the same measured noise moves that
   differential by an identical +0.0022, improves MAE (-0.013), and posts the corpus's
   best selected RMSD (-0.063 vs m1). Whatever the ensemble target buys is target
   smoothing, obtainable in either direction.
3. **Uncertainty calibration.** The one thing the ensemble target uniquely teaches -
   which pairs are undetermined - improves by rho +0.011, and the corrupted-label control
   improves by exactly +0.011.
4. **The dilution test.** At 13.4% ensemble coverage instead of 90.3%, the selection
   effect is the same size (-0.061 vs -0.045) and the distance effect reverses sign
   (+0.020 vs -0.032). A mechanism scales with exposure; noise does not.

And the single pre-registered dev pass does not replicate the tuning sign: +0.070
[-0.118, +0.259] against the matched baseline, +0.123 at 3W/10L inside the consensus
architecture.

This is the fifth independent instance of the rule sprint 7 established and the strongest
version of it, because here the better distance matrix was obtained by fixing a genuine
defect in the supervision rather than by tuning: **a better distance matrix is not a better
ranking, and the objective's misspecification is not in its labels.**

The axis is closed. Do not re-open it by re-training on ensembles with better features or
a bigger corpus: the ceiling is set by the noise-to-error ratio, and better features lower
the error, which *shrinks* the available gain rather than growing it.


<!-- -------------------------------------------------------------- -->

## Appendix: S8 / relax

> *Working paper, merged verbatim from `s8/relax_findings.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# Sprint 8 - AMBER as a refinement stage, not a ranker

`s8/relax.py`, tests `s8/test_relax.py`, results `s8/relax_sweep.json`,
`s8/relax_best.json`, `s8/relax_fe.json`, report snapshot `s8/relax_report.txt`.

## The question

Sprint 7 found all-atom ff14SB/GBn2 to be a worse **ranker** than the learned distance
prior: selected 4.037 vs 3.546 A, native at the 54th percentile, argmin at the native on
0/70 targets. But that scoring was `refine_coords(..., k_restraint=K_MODERATE, steps=0,
components=True)` with energy = nonbonded + solvation. `steps=0` is OpenMM's "minimise to
convergence", and at `k_restraint=10` every backbone N/CA/C is pinned - so **the fold was
never allowed to move**. Three things were therefore untested: relaxation, free energy
rather than single-point energy, and Amber in a generation/refinement role rather than a
ranking one.

## Protocol and reproduction

Instrument: the first 70 of the 126-target tuning instrument (`s7.debias.tuning_targets`),
which is exactly the set `s7/amber_native.py` ran on. Pools are built by the identical
recipe (`s7.amber_native.pool_for` + `s5.lib.windows_full` + BLOSUM62 top-K, ideal-geometry
rebuild through `protein_geometry.build_backbone_batch`), **including the plain - not
stable - argsort**, because the BLOSUM similarity has large tie sets and a different
tie-break gives a different pool (1A1P pool best 3.270 vs 3.336).

* pool best / pool mean reproduce `s7/amber_native.json` to six decimals on 1A13, 1A1P,
  1CB3;
* the native interaction energy reproduces **bit-for-bit**: −489.9138948277905 on 1A13.

Force field: genuine `amber14/protein.ff14SB.xml` + `implicit/gbn2.xml` through OpenMM.
`python -m s8.relax probe` on 1A13 (237 atoms), reproduced here in full:

```
forcefield  ['amber14/protein.ff14SB.xml', 'implicit/gbn2.xml']
openmm platform  CPU
system forces   ['CustomExternalForce', 'CustomGBForce', 'HarmonicAngleForce',
                 'HarmonicBondForce', 'NonbondedForce', 'PeriodicTorsionForce']
GBn2 present    True

k= 100.0 steps=    0    3.17s  E  497.7 ->  -323.0  moved(CA) 0.042 A
k=  10.0 steps=    0    4.95s  E  497.7 ->  -324.9  moved(CA) 0.113 A
k=   1.0 steps=    0    5.64s  E  497.7 ->  -335.4  moved(CA) 0.296 A
k=   0.0 steps=   50    0.38s  E  497.7 ->  -246.4  moved(CA) 0.409 A
k=   0.0 steps=  200    1.39s  E  497.7 ->  -326.3  moved(CA) 0.353 A
k=   0.0 steps= 1000    7.23s  E  497.7 ->  -339.8  moved(CA) 1.295 A
k=   0.0 steps=    0   23.48s  E  497.7 ->  -347.7  moved(CA) 2.409 A
```

`CustomGBForce` is GBn2; the restraint is the `CustomExternalForce`. The displacement
column is what makes the sweep meaningful: k=100 moves the fold 0.04 A and k=0-converged
moves it 2.4 A, so the schedule spans "cannot move" to "goes wherever the force field
wants".

Leakage: no native coordinate, distance, torsion or RMSD enters relaxation or scoring.
Every native read sits below a `# ---- ORACLE DIAGNOSTIC ----` line. A test scrambles the
native coordinates in the database and requires every candidate to come back bit-identical
(`test_pool_construction_reads_no_native_coordinates`). 12/12 tests pass.

Ran: 792 relaxations over 14 targets x 12 candidates for the schedule sweep (6 targets got
all seven schedules, 14 got the three decisive ones); the free-energy stage is expensive
(~15 min/target for 41 structures) and the box was shared with five other agents, so it
covers fewer targets - stated explicitly below rather than aggregated over.

## A. Relaxation moves candidates, and it is not refinement

`refine_coords`, not `refine` (which would project continuous torsions back onto the
discrete state library and confound everything). k_restraint = 0 is free minimisation, the
only setting in which the fold can actually move.

| k | steps | n | targets | CA moved | mean dRMSD | 95% CI (target-level) | W/L | frac improved |
|---|---|---|---|---|---|---|---|---|
| 100 | conv | 72 | 6 | 0.08 | **+0.011** | [+0.006, +0.017] | 0/6 | 0.153 |
| 10 | conv | 168 | 14 | 0.22 | **+0.012** | [+0.003, +0.022] | 4/10 | 0.351 |
| 1 | conv | 72 | 6 | 0.51 | +0.026 | [+0.015, +0.037] | 0/6 | 0.431 |
| 0 | 50 | 72 | 6 | 0.49 | **+0.087** | [+0.061, +0.112] | 0/6 | 0.236 |
| 0 | 200 | 168 | 14 | 0.65 | +0.020 | [−0.045, +0.084] | 4/10 | 0.446 |
| 0 | 1000 | 72 | 6 | 1.05 | +0.052 | [−0.023, +0.128] | 1/5 | 0.444 |
| 0 | conv | 168 | 14 | 1.39 | **−0.118** | **[−0.221, −0.014]** | **11/3** | 0.595 |

Candidates inside one target share a pool, a sequence and a builder, so the
candidate-level interval understates the SE; the target-level interval is the one quoted.

Three things to read off:

1. **Restrained relaxation is RMSD-neutral to slightly harmful**, confirming the docstring
   claim in `amber_refine` on a much larger and much less favourable sample (that claim was
   measured on 12 targets from *oracle* `floor.floor` starts). At the exact sprint-7 setting
   (k=10) the cost is +0.012 A - real but negligible.
2. **A short free minimisation is significantly harmful** (+0.087 A at 50 steps, pool best
   worse on 6/6 targets). The first steps relieve the worst clashes by locally distorting
   the chain.
3. **Converged free minimisation genuinely lowers the RMSD distribution** - −0.118 A, CI
   excluding zero at the target level, better on 11 of 14 targets, 59.5% of individual
   candidates improved.

### But it is a distribution compressor, not a refiner

Point 3 looks like the result the sprint was hoping for. It is not. Two decompositions kill
it.

**By starting quality** (terciles of rmsd0 within each target, converged free minimisation):

| tercile | mean dRMSD | frac improved |
|---|---|---|
| near-native | **+0.056** | 0.518 |
| middle | −0.045 | 0.536 |
| far | **−0.364** | 0.732 |

The entire gain is on the candidates that were furthest from the native, and the
near-native third gets slightly *worse*. That is the signature of collapse toward a generic
basin, not of movement toward each structure's own native.

**The compactness control confirms it.** Within-target spread of the radius of gyration
falls 0.906 → 0.751 A (−17%) under converged free minimisation, while |rg − rg_native|
falls only 1.127 → 1.011. The structures become more like *each other* faster than they
become more like the native, and rho(drg, dRMSD) = +0.217 - the ones that shrink most are
the ones that improve most.

**And pool best does not move**: +0.013 A [−0.197, +0.224], 7 W / 7 L. A distribution that
shifts down while its minimum stays put is exactly what compression looks like.

### The ceiling check that settles it

A relaxation stage would *add* structures to the pool, never replace it, so the honest
generation number is pool best over the union of original and relaxed. That looks
encouraging: 2.336 → 2.216 A, **−0.120 A [−0.251, +0.011]**, with a relaxed structure
taking the best slot on 7/14 targets.

It is an artefact of the sub-pool. Those 12 candidates are the top-12 BLOSUM hits, whose
pool best is 2.336 A. The pipeline's actual K=500 retrieval pool has a best of **1.636 A**
on the same 14 targets (`s7/debias_cache`). Against that:

| k | steps | relaxed structures better than the K=500 pool best | targets with any |
|---|---|---|---|
| 10 | conv | 1 / 168 | 1 / 14 |
| 0 | 200 | 1 / 168 | 1 / 14 |
| 0 | conv | **0 / 168** | **0 / 14** |

**Zero converged-free-relaxed structures out of 168 beat what retrieval already had, on
0 of 14 targets.** The union gain is measured against a weak sub-pool and does not survive
contact with the pool the pipeline actually builds. This is the same shape of error the
project has recorded twice before (`decoy-bank-not-a-pool-proxy`, `better-matrix-worse-ranking`):
an intervention that improves a weak reference and buys nothing against the real one.

**The single instructive case.** Exactly one of the 168 relaxed candidates started better
than the whole K=500 pool: 1KZ2 candidate 3, at 1.121 A against a K=500 best of 1.281 A.
Pinned relaxation (k=10) improves it to **1.077 A**. A short free minimisation holds it at
1.193 A. Converged free minimisation destroys it: **1.820 A**. The one structure in the
experiment that was genuinely better than retrieval is the one the force field's own
minimum is furthest from - which is what "the minima are not near the native" looks like
on a single structure.

**The distribution really does shift, though.** Under converged free minimisation the
fraction of candidates below 2.0 A rises 0.149 -> 0.185 and below 1.5 A rises 0.089 ->
0.143. So the compression is not cosmetic: it moves genuinely more structures under the
thresholds the programme cares about. It just never produces the *best* one, and the
structures it promotes past a threshold are ones that were already close to it.

The result also survives the influence check: dropping the single most favourable target
(1DEP, -0.557) moves the mean from -0.118 to -0.084, still negative.

### Why the ranking failure and the relaxation failure are the same fact

Free minimisation moves a candidate 1.39 A of CA-RMSD, and where it lands is a more compact
structure that is on average slightly closer to *every* native and no closer to *its own*.
A force field whose minima are placed by generic compaction rather than by sequence-specific
structure will (i) fail to refine, and (ii) rank by compactness - which is precisely what
sprint 7 measured when partialling out radius of gyration cut Amber's in-band skill to
+0.159 against the distogram's +0.479. One mechanism, two symptoms.

## B. Free energy instead of single-point energy - under-powered, reported as such

Four native-free scores that price basin **width** as well as depth, each computed on the
candidate through the same builder and relaxation as the native:

* `strain` = E(built) − E(freely relaxed): the candidate's distance from its own force-field
  minimum, which isolates builder strain from conformational quality;
* `E_free`: the energy *at* the candidate's own relaxed minimum, a depth the builder cannot
  corrupt;
* `F_qh` = ⟨E⟩ − kT·S_qh, with S_qh = ½ Σ log λ of the superposed CA covariance of a 4 ps
  300 K Langevin ensemble (quasi-harmonic, up to an additive constant that is identical for
  every candidate of one target and so cancels in the within-target ranking that is all we
  use);
* `F_boltz` = −kT log ⟨exp(−E/kT)⟩ over the same ensemble, which mixes depth and width
  without assuming harmonicity;

plus `width` (the RMS fluctuation itself) and `S_msf`. The ensemble is seeded from the
*restrained*-relaxed candidate so it samples the candidate's own basin, and uses a
hydrogen-mass-repartitioned copy of the same ff14SB/GBn2 system (constraints change the
dynamics that explore the basin, not the potential).

**This stage completed only 1 of its 24 targets before the box filled up.** One target is
not a result and no free-energy claim is made from it. The machinery is committed,
resumable per target, and the report table is wired to both benchmarks (distance objective
3/126 argmin and 36.8th percentile; single-point Amber 0/70 and 54th percentile) so it can
be finished directly.

The one thing worth recording from the single completed target, because it is a *negative*
about a design choice rather than a positive claim: on 1A13 the **total** energy places the
native at the 74th percentile where the interaction-only energy places it at the 28th, which
independently supports sprint 7's choice to drop bond/angle/torsion from the ranking energy.

## C. The integrated chain

Folded into B (the same run re-scores the pinned- and freely-relaxed CA traces with the same
fold-held-out distogram, so `retrieve → relax → score` is measured without a second OpenMM
pass), and therefore equally under-powered. What A already settles is the first link:
**relaxation does not improve pool best**, so there is no pool-best gain for the rest of the
chain to translate. Per the protocol, a chain result is not claimed on one target.

## Verdict - what role should AMBER hold?

**Not a ranker, and not a refiner of backbone conformation either.**

* As a *ranker*, sprint 7 settled it and nothing here rescues it.
* As a *refiner*, this sprint settles it: at the restraint strengths that preserve the fold,
  relaxation is RMSD-neutral to slightly harmful (+0.011 to +0.026 A); at the setting that
  genuinely moves the fold it lowers the RMSD distribution by compressing it toward a
  generic compact basin, improves the far-from-native candidates only, leaves pool best
  untouched, and produces **zero** structures better than the K=500 pool already contains.

What it is genuinely good for is the thing its own module docstring already claims and this
sprint confirms at scale: **stereochemical repair**. An ideal-geometry rebuild sits tens of
thousands of kcal/mol above its own minimum, and restrained relaxation at k=10-100 removes
that for 0.011-0.012 A of CA-RMSD - i.e. essentially free. That is a legitimate material
role: every all-atom structure the pipeline emits should pass through it so that what
leaves the system is physically valid, and no downstream consumer inherits builder strain.
It is a **validity** stage, not an **accuracy** stage.

The remaining honest gap is experiment B. Entropy as a leading term for marginally stable
9-16mers is still a live idea, it has not been refuted here, and the instrument to test it
is committed and resumable.

## What is left on the table, and exactly how to finish it

The box sat at 85-95% RAM for the second half of the sprint (five other agents plus
Chrome), and every stage here refuses to start above 78% and stops cleanly above 90%, so
three queued runs did not get compute. All three are resumable per target - a kill costs
at most one structure - and none needs any code change:

```
# A, sharpest form: relax the K=500 pool's best three members (oracle-selected diagnostic)
NTOP=3 KPOOL=500 PIDS=... python -m s8.relax best        # ~4 min/target, writes relax_best.json

# A, more targets at the three decisive schedules
PIDS=1RG4,1TOR,2BFI,2FBU,2LER,2LU6,2MAA,2MID,2MK7,2MQ2,2N9A,2NDM NCAND=12 \
    python -m s8.relax wide                              # ~6 min/target, appends relax_sweep.json

# B, the free-energy table (the real gap in this sprint)
NSNAP=100 STRIDE=10 FREE_STEPS=200 KA=40 python -m s8.relax fe   # ~15 min/target

python -m s8.relax report                                # aggregates whatever exists
```

`report` prints every table in this document plus the free-energy table paired against
both benchmarks and against `s7/dist_baseline.json` on identical pools. B needs roughly
4 hours of quiet box to reach 16 targets, which is the point at which its native-percentile
column can be compared to 0/70 and 54.0 with any power.

## Caveats stated up front

* 14 targets x 12 candidates for A's three main settings; 6 targets for the four
  sweep-only settings. Rows of the schedule table cover different target sets, so rg0 and
  pool-best columns are comparable *within* a row, not down the table.
* The candidates relaxed are the top-12 BLOSUM hits, not a random sample of the pool. The
  ceiling check against `s7/debias_cache` is what protects the conclusion from that choice;
  the complementary test - relaxing the K=500 pool's *best* members directly - is
  implemented (`python -m s8.relax best`) and was queued behind the memory gate.
* The Langevin ensemble is 4 ps with 100 snapshots, so S_qh over up to 42 CA degrees of
  freedom is a noisy and downward-biased estimate. It is applied identically to every
  candidate and to the native, which is what the within-target ranking needs, but it is not
  a converged conformational entropy.


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / bias

> *Working paper, merged verbatim from `s9/bias_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-bias. The shared bias is real, 13 SE, fully explained - and it is the architecture's own signature, not an error

`s9/bias.py` (9 resumable stages), `s9/test_bias.py` (40 checks, 0 failures),
`s9/bias_{characterise,fit,angle,frechet,controls,leak,lamcurve,lowdim_diag}.json`.

Instrument: the 126-target tuning set, K=500 BLOSUM pools, read through
`s8.consensus2.load_cache`. `stage_build` asserts the three numbers this study stands on
before anything else runs - shipped score **3.4540**, pool best **1.7108**, and S8-11's
synthesis arm `sc|75:fit` in its deployable torsion-free form **3.2041**. **The 60-target
benchmark was not touched by any stage in this file, and no dev-24 pass was spent.**

## The question, and the distinction the study turns on

S8-11 measured the pool's error decomposition: averageable independent error
`sqrt(3.748^2 - 3.462^2) = 1.437 A` against **85.3% of the mean square being bias common to
every candidate**. Averaging harvests the independent 14.7% and cannot touch the rest.

That is a **within-target** statement. It says the 500 candidates of one target miss in the
same direction. It does **not** say that different targets miss in the same direction, and
only the second kind of bias is learnable from anything but target-specific information.
Nobody had measured the second kind.

**Measured: 11.7% of the displacement field is common across targets.** Nine-tenths of it
is target-specific - which is the recognition problem again, re-expressed as a regression.

## 1. The displacement field, in a frame that can carry a correction

For every target: superpose the synthesised structure onto the native, take
`d = nat - F_superposed`, and express `d` in a per-residue **right-handed (tangent, normal,
binormal) frame built from the synthesis itself**. Local components are invariant under the
global superposition RMSD is free to choose, so a correction expressed in them is
equivariant; the frame is a cross product and therefore a rotation and never a reflection,
so a mirrored structure gets a mirrored frame (this repo has a recorded chirality defect,
S7-4, and a mirror-blind representation would let a correction learn the wrong hand).
`t_frames_right_handed`, `t_local_components_rotation_invariant`, `t_frame_is_chiral`.

Each model is scored twice, because the two scores disagree and only one of them is the
project's metric. `r2` is the fraction of the POOLED displacement sum-of-squares removed -
the statistic the 85.3% figure lives in. `rmsd` is the mean over targets of
`sqrt(mean |d - dhat|^2)`, an **upper bound** on what the deployable arm can score (a fresh
Kabsch can only lower it). The pooled `r2` flatters every global mode because it is
dominated by the targets with the largest displacement, while the metric weights all 126
equally.

| model | r2 (pooled SS) | rmsd | d vs 3.204 | 95% CI | W/L |
|---|---|---|---|---|---|
| rigid - one global translation | 0.0000 | 3.2041 | +0.0000 | - | 0/0 |
| radial - `d = alpha (F - centroid)`, one alpha | 0.0496 | 3.1938 | −0.0103 | [−0.0704, +0.0499] | 66/60 |
| **posbin** - mean local field per normalised position | **0.0664** | 3.1558 | **−0.0483** | [−0.0961, −0.0005] | 77/49 |
| posbin + radial, jointly | 0.0864 | 3.1421 | −0.0620 | [−0.1255, +0.0016] | 69/57 |
| *radial_per_target* - ONE number per target | *0.2719* | *2.8288* | *−0.3752* | *[−0.4873, −0.2632]* | *126/0* |
| *pca_5* - FIVE numbers per target | *0.3928* | *2.5494* | *−0.6546* | *[−0.7797, −0.5296]* | *114/12* |

Italic rows are ORACLE. The gap between `posbin` (−0.048, everything that is shared) and
`pca_5` (−0.655, five target-specific numbers) is the whole result in one line: **the field
is low-dimensional per target and the missing information is five numbers, not a
3n-vector.**

**The decisive cross-target statistic.** Mean pairwise cosine between two targets' resampled
local-frame fields: **+0.1176 ± 0.0027**, against a sign-flip null of **+0.0002 ± 0.0029**
- 40 SE above the null, so real; 69.5% of pairs positive; shared fraction of the field
sum-of-squares **0.1171**.

**Two things this rules out as obstacles.**
- The ideal-geometry manifold is **not** the constraint. Applying the true field and
  re-projecting lands at **0.0952 A**.
- The radial mode is an artefact of imperfect superposition, not compaction. After Kabsch,
  `alpha = (|nat'|/|F'|) c - 1` with overlap `c`; the measured 1.019 and 0.859 reproduce
  the fitted −0.125 exactly, and applying it is worth −0.010 A with a CI spanning zero.

## 2. What the shared bias IS: the synthesis is over-curved, at 13 SE

The mean local-frame profile is almost entirely one component:

| | tangent | **normal** | binormal |
|---|---|---|---|
| mean over 126 targets (A) | −0.052 ± 0.044 | **+0.863 ± 0.066** | +0.000 ± 0.040 |
| frac of targets positive | | **0.913** | |

Moving a CA along its local curvature normal moves it toward the chord of its neighbours,
which **straightens** the local bend. So the profile says the synthesis is systematically
over-curved, and the direct measurement agrees:

| | CA-CA-CA interior bond angle |
|---|---|
| projected synthesis | **93.21 deg** |
| native | **103.99 deg** |
| difference | **+10.78 ± 1.00 deg**, 85.7% of targets |

That reproduces `s8/audit8`'s 93.5 / 104.0 to a tenth of a degree, from a displacement field
rather than from torsions - an independent confirmation on a different measurement.

The signal is not a high-frequency zigzag: the adjacent-residue correlation of the normal
component is −0.041. It is a uniform, smooth over-curvature.

## 3. Why it is there - and why that means it is not an error

`s9/bias.py`, the m-sweep over the raw coordinate average:

| set averaged | CA-CA-CA angle | mean CA-CA step |
|---|---|---|
| individual candidates | 101.27 deg | 3.81 A |
| m = 1 (the medoid - a real library fragment) | 103.04 | 3.810 |
| m = 3 | 110.27 | 3.432 |
| m = 10 | 115.95 | 3.205 |
| m = 25 | 119.62 | 3.083 |
| **m = 75 (the incumbent)** | **121.62** | **2.961** |
| NATIVE | 103.99 | 3.812 |
| **projected synthesis** | **93.21** | **3.804** |

**The candidates are fine.** Individually they have native-like bond angles and native-like
bond lengths. Averaging them contracts the chain monotonically in m - the mean of points
scattered about a curve lies inside it - and in raw coordinates that reads as
*straightening*, all the way to 121.6 deg at m=75. The projection then has to restore an
exact 3.804 A step while staying near that contracted object, and it buys the missing 28%
of contour length back **as curvature**, overshooting to 93.2 deg.

So the one systematic, 13-SE, perfectly interpretable component of the shared bias is
**the deterministic signature of the consensus-then-project architecture**, not an error in
what the pipeline believes about the target. It carries no information about the answer,
which is why every attempt to remove it fails - in two different ways, below.

## 4. Removing it by displacement is illegal geometry

The first arm to complete settled it:

| | selected | mean CA-CA step | d vs 3.204 | 95% CI | W/L |
|---|---|---|---|---|---|
| `ridge@0.5` post-projection | 3.2010 | 3.804 | −0.0031 | [−0.0128, +0.0066] | 63/63 |
| `ridge@0.5` **pre**-projection | **3.0854** | **3.134** | −0.119 | | |
| `ridge@0.5~rs` - step restored, then projected | **3.5056** | 3.804 | **+0.3015** | [+0.2088, +0.3943] | 29/97 |
| `const@1` **pre**-projection | 3.1556 | **2.697** | −0.049 | | |
| `const@1~rs` - step restored, then projected | 4.1508 | 3.804 | **+0.9467** | [+0.7880, +1.1054] | 15/111 |

A per-residue displacement that pushes every CA toward the chord of its neighbours *is* a
contraction. The apparent −0.119 A is the S8-11 trap measured a second time: the raw
coordinate average scores 3.048 with 2.961 A bonds and 37% of its apparent gain is geometry
rather than conformation. Here 100% of it is. Every arm in `bias_fit.json` therefore
reports its pre-projection CA-CA step, and `apply_correction` has a `rescale` mode that
restores the step before projecting so the shape change can be read on its own.

## 5. Removing it legally, in the one scale-free representation, is worse

The bond angle does not care about scale, so it is the representation in which the measured
bias can be corrected without leaving the peptide manifold. `stage_angle` re-projects the
same coordinate average under

    minimise  CA-RMSD(build(phi, psi), A)  +  mu * mean (theta_i - theta*)^2

with `theta*` the mean native bond angle over TRAINING FOLDS ONLY (104.22 / 103.69 / 104.48
/ 103.87 / 103.68 deg - stable). `mu = 0` reproduces the incumbent projection bit-for-bit
and is asserted to. Bond length stays exact at 3.804 A at every mu.

| mu | selected | CA-CA-CA angle | d vs mu=0 | 95% CI | W/L |
|---|---|---|---|---|---|
| **0 (the incumbent)** | **3.2040** | 93.22 | - | - | - |
| 8 | 3.2230 | 97.61 | **+0.0189** | [+0.0016, +0.0363] | 52/74 |
| 32 | 3.2844 | 100.61 | **+0.0803** | [+0.0497, +0.1110] | 39/87 |
| 128 | 3.4618 | 103.17 | **+0.2578** | [+0.1922, +0.3233] | 31/95 |
| 512 | 3.5384 | 103.80 | **+0.3343** | [+0.2530, +0.4157] | 27/99 |

**Monotone, in the wrong direction, with every interval excluding zero.** Driving the bond
angle to the natives' own 104.0 deg takes the synthesis from 3.204 past the 3.454 shipped
baseline to 3.538 - it throws away everything S8-11 gained and more.

So the 13-SE bias is **not an error to be corrected**. It is the projection's optimal
response to a contracted average, and forcing the "right" angle only moves the chain away
from the data it is fitted to. The two facts - that the angle is wrong, and that fixing it
hurts - are consistent exactly because the angle is not carrying information about the
target.

## 6. The model, and what it can and cannot learn

**The model.** A per-residue ridge from 55 deployable features to the local-frame
displacement, fitted leave-fold-out, with `lambda` **and** the step size chosen by a nested
LFO *inside* the training folds. Features: normalised position and distance to a terminus;
the synthesis' own pseudo bond angle and pseudo dihedral (which is where predicted
secondary structure enters); the radial vector in local-frame components; the candidate
set's per-residue spread and its mean pairwise RMSD; the distogram's per-residue predicted
SD; **the distogram force** `g_i = sum_j w_ij (d_ij - E[d_ij]) u_ij`, the direction the
out-of-fold distogram wants residue i to move, in local-frame components; five
physico-chemical sequence channels and their 3-window means; and 16 ESM PCA dimensions.

**What it learns, on held-out targets** (`bias_lamcurve.json`, relative field error; 1.000
is predicting zero):

| | relative held-out field error | fraction of the field explained |
|---|---|---|
| predict nothing | 1.00000 | 0% |
| intercept only - a single constant 3-vector | 0.94485 | **5.5%** |
| ridge, one-hot sequence, lambda=1000 | **0.90721** | **9.3%** |
| ridge + ESM, lambda=1000 | 0.90812 | 9.2% |
| lambda=10000 | 0.92213 | | |

`lambda = 1000` is an interior optimum on every fold and in the curve. **ESM is +0.0009
worse than the sequence channels alone**, which is the third independent instance in this
project (S7-11, S8-13, here) of ESM adding nothing to short peptides. So all 55 features
together add **3.8 percentage points** over a single constant vector.

**The low-dimensional arm, and why it cannot work** (`bias_lowdim_diag.json`). `pca_5`
says five per-target numbers are worth −0.655 A, so the missing information is five numbers
and not a 3n-vector. Learn the 5-mode basis on training folds and regress its coefficients
on 18 target-level features:

| lambda | held-out R^2 per mode | field error | (mean-only control in the same basis) |
|---|---|---|---|
| 10 | −0.145 −0.193 −0.057 −0.140 +0.008 | 0.96353 | |
| 100 | −0.027 −0.032 +0.015 −0.057 +0.017 | 0.94437 | |
| 1000 | +0.003 +0.018 +0.018 +0.017 +0.015 | 0.93828 | **0.94125** |

**The coefficients are unpredictable.** At the best regularisation the five modes are
explained to 0.3-1.8%, and at weaker regularisation the R^2 goes *negative* - worse than
predicting the training mean. Predicting them beats the mean-only control in the same basis
by 0.003 of field error. The per-residue model does better (0.907) than the target-level
one (0.938) precisely because the little that is predictable is **local**, not a global
deformation.

**The arm table**, 126 targets, LFO throughout. `pre` is the corrected coordinates before
re-projection and `CA-CA` is that object's mean bond step (ideal 3.804 A).

| arm | selected | median | SD | <2.0 A | <1.5 A | worst | best | `pre` | `CA-CA` | d vs 3.204 | 95% CI | W/L | d vs 3.454 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SHIPPED score | 3.4540 | | | 0.214 | 0.159 | | | | | +0.250 | | 46/80 | - |
| **synthesis `sc\|75:fit` (S8-11)** | **3.2041** | 2.966 | | 0.278 | 0.190 | | | | 3.804 | - | - | - | −0.2499 [−0.383, −0.117] |
| **`ridge@0.5`** | **3.2010** | 2.991 | 1.719 | 0.278 | 0.190 | 8.215 | 0.198 | *3.0854* | *3.134* | **−0.0031** | [−0.0128, +0.0066] | 63/63 | −0.2530 |
| **`ridge_auto@1`** (LFO step size) | **3.2055** | 3.045 | | 0.278 | | | | *3.0702* | *2.651* | **+0.0015** | [−0.0212, +0.0241] | 61/65 | −0.2484 |
| `dgforce@1` | 3.2060 | 2.971 | | 0.270 | | | | 3.2064 | 3.781 | +0.0019 | [−0.0043, +0.0081] | 57/69 | −0.2480 |
| **`const@1`** (the constant control) | **3.2410** | 2.991 | 1.715 | 0.262 | 0.151 | 8.269 | 0.394 | *3.1556* | *2.697* | **+0.0369** | [+0.0188, +0.0551] | 45/81 | −0.2130 |
| `ridge@0.5~rs` (geometry restored) | 3.5056 | 3.051 | 1.923 | 0.238 | 0.151 | 10.591 | 0.528 | 3.5111 | **3.804** | **+0.3015** | [+0.2088, +0.3943] | 29/97 | +0.0516 |
| `const@1~rs` (geometry restored) | 4.1508 | 3.571 | 1.764 | 0.024 | 0.008 | 10.730 | 1.463 | 4.1965 | **3.804** | **+0.9467** | [+0.7880, +1.1054] | 15/111 | +0.6968 |
| *ORACLE - the true field, re-projected* | *0.0952* | | | | | | | | 3.804 | *−3.109* | | *126/0* | |
| *ORACLE - pool best* | *1.7108* | | | | | | | | | | | | |

Read the `CA-CA` column first. **Every pre-projection gain sits on an object that is not a
peptide**, and the `~rs` rows - the same correction with the mean step restored before
projecting - are the price of making it one.

## 7. Controls

**Control 1 - a constant global correction.** `const@1` is exactly that: the training
folds' mean local-frame field, resampled onto the held-out target's length and applied
identically to every one. It is **3.2410, +0.0369 [+0.0188, +0.0551], 45/81** - significantly
worse than doing nothing. Before re-projection it is 3.1556, a −0.049 A gain that reproduces
`stage_characterise`'s `posbin` row (−0.0483) to a thousandth. So **the constant IS the
whole shared field**, the mandate's first control is decisive, and the answer is that there
is nothing for a learned model to be better than: the ridge and the constant are both null,
the ridge only ever beats the constant by contracting the chain further.

**Control 2 - the shuffled target.** Apply target A's predicted correction to target B
(length-matched by resampling), on the best deployable arm `ridge@0.5`.

| | selected | d vs synthesis | 95% CI | W/L | d vs shipped |
|---|---|---|---|---|---|
| `ridge@0.5` (the arm) | 3.2010 | −0.0031 | [−0.0128, +0.0066] | 63/63 | −0.2530 |
| **shuffled target** | **3.2103** | **+0.0062** | [−0.0004, +0.0129] | 55/71 | −0.2437 |
| shuffled **vs the arm** | | **+0.0094** | [−0.0010, +0.0197] | 53/73 | |
| null features (Gaussian noise through the same pipeline) | 3.2243 | +0.0202 | [+0.0089, +0.0316] | 46/80 | −0.2297 |

**The shuffled correction is 0.0094 A worse than the real one, with a CI that touches zero.**
So the correction carries at most nine thousandths of an Angstrom of target-specific
information - which is the same statement as the arm being null, reached from the other
side. The null-feature control (+0.0202) is the floor for "a small step in a learned
direction" and both real arms sit above it, so the pipeline is not broken; there is simply
almost nothing to learn.

**Control 3 - leave-fold-out throughout.** Section 8.

**Control 4 (not asked for, and the one that mattered) - restore the geometry.** The `~rs`
rows. Without them `ridge@0.5`'s pre-projection 3.0854 reads as a −0.119 A result.

**Concentration.** No arm here has a gain to concentrate: every one has a median paired
difference of at least −0.0007 and a fraction improved of at most 0.500, and **dropping the
ten largest per-target gains turns every arm positive.**

| arm | d | median d | frac improved | d with the top 10 gains dropped |
|---|---|---|---|---|
| `ridge@0.5` | −0.0031 | −0.0007 | 0.500 | **+0.0057** [−0.0031, +0.0145] |
| `ridge_auto@1` | +0.0015 | +0.0028 | 0.484 | +0.0229 |
| `dgforce@1` | +0.0019 | +0.0020 | 0.452 | +0.0086 |
| `const@1` | +0.0369 | +0.0363 | 0.357 | +0.0529 |

For `ridge@0.5` the full ladder is +0.0015 [−0.0077, +0.0107] on n=121 after dropping the
top 5, +0.0057 [−0.0031, +0.0145] on n=116 after 10, and +0.0134 [+0.0053, +0.0216] on
n=106 after 20 - the inter-quartile range of its per-target difference is [−0.036, +0.026],
i.e. symmetric noise about zero.

Measured as a by-product, **the incumbent's own concentration** - the number the mandate
warned about, which had not been reported for the arm this study builds on:

| | value |
|---|---|
| S8-11 synthesis vs the shipped baseline | **−0.2499** [−0.383, −0.117], 80/46 |
| median target | **−0.1085** |
| fraction of targets improved | 0.635 |
| ten largest gains as a share of the total | **0.606** |
| drop the top 5 | −0.1598 on n=121 |
| drop the top 10 | **−0.1070** on n=116 |
| drop the top 20 | −0.0127 on n=106 |

S8-11's gain is real at the median (−0.109) and survives dropping ten targets, but 61% of
it sits in ten and it is essentially gone after twenty. That is worth recording next to
the headline.

## 8. Leakage audit

**Clean.** `stage_leak` takes each held-out target, replaces its native with NaN, and
recomputes the entire deployable chain - features, model prediction, correction,
re-projection - asserting the output is bit-identical: **15 held-out targets across all 5
folds, worst |diff| 0.000e+00.** The model itself is trained on the other folds, where the
native is a legitimate label.

Four further poisoning tests run inside the suite: `t_features_ignore_native` (20 targets,
worst 0.000e+00), `t_target_features_ignore_native` (20, 0.000e+00),
`t_model_lfo_holds_out` (0.000e+00) and `t_lowdim_basis_is_training_only` (0.000e+00 - the
deformation basis is fitted without the held-out target's field).

No arm reads a native at inference. The oracle rows (`radial_per_target`, `pca_k`,
`oracle_reprojected`) are labelled ORACLE everywhere they appear and are ceilings, not arms.

## 9. Averaging ON the manifold instead of beside it - the over-curvature is not the order

Section 3 says the incumbent contracts and then repairs. The obvious objection is that
those two steps are done in the wrong order, so `frechet_fit` does the averaging **on** the
manifold:

    minimise  mean_k  CA-RMSD(build(phi, psi), W_k)^p     over (phi, psi)

Every object it compares is a peptide, so there is no contraction to create and none to
repair. `p = 2` is the Frechet mean, `p = 1` the manifold geometric median (the medoid's
analogue). S8-11's 20 consensus operators covered selection, coordinate averages, trimmed
and weighted averages, and the torsion **circular** mean - a different object, and +0.618 A
worse - but not this one. Run on the first 40 targets (each cell is an optimisation over 75
candidates per target, and the question is about an 11-degree angle rather than a
hundredth of an Angstrom).

| construction, m members | CA-CA-CA angle | selected (40 targets) | d vs the incumbent | 95% CI | W/L |
|---|---|---|---|---|---|
| NATIVE | 103.22 | | | | |
| Euclidean average, m=5, **before** projection | 113 (interp.) | | | | |
| Euclidean average, m=75, **before** projection | 121.62 | | | | |
| **incumbent: Euclidean average m=75, projected** | **93.83** | 2.9856 | - | - | - |
| **manifold Frechet mean, m=5** | **98.33** | 3.0509 | +0.0653 | [−0.1002, +0.2309] | 20/20 |
| **manifold Frechet mean, m=25** | **95.51** | 2.9800 | −0.0055 | [−0.0673, +0.0562] | 21/19 |
| **manifold Frechet mean, m=75** | **93.61** | 3.0012 | +0.0156 | [+0.0046, +0.0266] | 13/27 |

**The manifold mean over-curves exactly as much, and in the same way, as
average-then-project.** So the over-curvature is *not* an artefact of doing the two steps in
the wrong order. It is intrinsic to any RMSD-optimal mean of a diverse set of peptides:
minimising the mean squared RMSD to a spread of structures pulls the estimate toward the
centre, and at a fixed 3.804 A bond length the only way for a chain to be central is to curl
up. It is shrinkage in shape space, and it is the same mechanism S8-12 measured when
converged AMBER minimisation compressed structures (rg spread down 17%, near-native cases
made worse and far ones better).

A practical note falls out: **average-then-project is not an approximation to the manifold
mean that costs accuracy - it is as good, and about thirty times cheaper.**

## 10. How much of the shared bias is learnable - the direct answer

**Statistically: about 9%.** Of the displacement field's pooled sum of squares, 11.7% is
common across targets in sample and a 55-feature ridge recovers **9.3% held out** - 5.5
points of which is a single constant 3-vector and 3.8 points everything else.

**In CA-RMSD, harvesting all of it is worth −0.048 A** [−0.096, −0.000] and that is an upper
bound obtained by ignoring the geometry.

**Deployably it is worth nothing, and the reason is not statistical.** The shared component
is one mode - the synthesis is over-curved by 10.78 ± 1.00 deg of CA-CA-CA bond angle - and
that mode is **the deterministic signature of averaging-then-projecting**, not an error in
what the pipeline believes. Individual candidates have native-like angles and native-like
bonds; averaging m of them contracts the chain, and the projection buys the contour length
back as curvature. Correcting it as a displacement contracts the chain to 2.70-3.13 A per
step, and restoring the step costs +0.30 to +0.95 A. Correcting it in the scale-free bond
angle costs +0.019 to +0.334 A, monotonically, all the way to the natives' own value.

**The part that would matter is target-specific and unpredicted.** Five per-target numbers
would be worth −0.655 A [−0.780, −0.530]. Their held-out R^2 from 18 target-level features
is 0.003-0.018, and negative at weaker regularisation. The per-residue model does better
than the target-level one only because the little that is predictable is local.

### The sharpest statement of the limit

This axis was chosen because it sidesteps discrimination, and it did: it never had to say
which candidate was better. It still landed on the same wall, and it landed there with a
number attached. **The residual of the synthesised structure is genuinely target-specific
information that no signal available to this project carries** - not the distogram
(0.002 A as a force), not ESM (+0.0009 worse than one-hot), not the candidate set's own
spread, not the sequence, and not the synthesis' own geometry. Nine-tenths of the error is
per-target, and the one-tenth that is shared is the method looking at itself.

"85.3% of the error is bias common to every candidate" was the premise of this study. The
correct reading is now: **85.3% common to every candidate, 11.7% common to every target.**
The first number invited a regression; the second is what the regression could actually see.

### No dev pass was spent, deliberately

Every arm is null or worse on the 126-target tuning instrument - the best is `ridge@0.5` at
−0.0031 [−0.0128, +0.0066] with 63 wins and 63 losses. dev-24's paired SE on this arm family
is ~0.35 A. Spending a pre-registered pass on an effect the tuning instrument prices at
three thousandths of an Angstrom would produce a number that means nothing (S8-13's
discipline note). `DEV_ARM` is `None` in the code and `stage_dev` asserts it before running.

The 60-target benchmark was not touched.

### What this does and does not license for the next sprint

- **Do not** re-open post-hoc correction of the synthesised structure with better features.
  Better features lower the model error, and the ceiling here is not model error - it is
  that 90% of the field is per-target and the 10% that is not is an artefact of the method.
- **Do not** add a bond-angle or curvature term to the projection expecting accuracy. It is
  measured, monotone, and negative: `s9/bias_angle.json`.
- **Do** read the `CA-CA` step of anything that reports a gain before re-projection. Two
  independent arms in this file gained 0.05-0.12 A that way and both were contraction.
- The open question this study leaves is the one it could not answer: what carries the five
  per-target deformation numbers. It is not in the pool's spread, the distogram, or the
  sequence.


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / ceiling

> *Working paper, merged verbatim from `s9/ceiling_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-6. The ceiling of in-band discrimination: the information is absent, not unexploited

`s9/ceiling.py` (10 resumable stages), `s9/test_ceiling.py` (34 checks, all passing),
`s9/ceiling_{check,oracle,capacity,global,ablate,native,filter,amber,leak}.json`,
`s9/ceiling_cache/` (126 per-target feature blobs).

**The claim under test.** Nine sprints converged on: *within the near-native band, no
native-free signal can discriminate; recognition is the bottleneck and it is unachievable
with the information available at inference.* Every previous instance used ONE channel, a
small combination, or one architecture. This is the kitchen-sink experiment: the most
capable in-band discriminator the repository's caches permit, every channel at once, oracle
RMSD supervision on training folds, leave-fold-out. It is a CEILING measurement, not a
deployment attempt. **No dev-24 or benchmark pass spent; a test asserts the module cannot
reach either, and that every pdb id in every output file is one of the instrument's 126.**

## The instrument, asserted before a single feature was built

`stage_check`. The pools behind `s8/generate_univ`, `s8/inband_cache`, `s8/integrate_chan`,
`s8/consensus2_cache` and `s9/relative_cache` are **the same pool: worst |diff| 0.000e+00 on
all 126 targets, on both `rr` and the raw coordinates**. That is what makes concatenating
columns from four modules BY INDEX a reproduction rather than an approximation. The
instrument reproduces **3.4540 selected / 1.7108 pool best** exactly, the shipped in-band rho
at **+0.1308** (S8-8's +0.126, S9-3's +0.131), the score-filtered synthesis at **3.2005**
(S8-11's 3.204) and the oracle filter at **2.0947, -1.106 [-1.289, -0.922]** (S9-3's
reproduction to three decimals). `signal_block` recomputed here matches `s8/inband_cache` to
1e-5 (float32 storage), so the lifted columns are checked, not trusted.

**105 per-candidate features + 13 target-level context columns**, each given to the model in
two views (within-**pool** z-score and within-**pool** percentile - never within the band,
which is not knowable at inference), for a **223-column** design matrix:

| block | n | what |
|---|---|---|
| `dg` | 11 | distogram family: the shipped Bayes risk, two sd-reweightings, L1, z-dev, confidence slices, bagged |
| `dgs` | 24 | **NEW** - per-shell decomposition (\|i-j\| in 2-4 / 5-8 / 9-16 / 17+, each as MAE, risk, z, in-band fraction), low-sd and high-sd slices, distance-profile rank agreement |
| `lg` | 15 | Legacy's 11 terms + hb / burial / packing subsets + the fitted total |
| `ph` | 15 | H-bond counts, contacts, contact order, strain, Ramachandran, SS self-consistency, rg, rebuild displacement, **CA-CA-CA pseudo-angle mean/sd and positive-phi rate** |
| `iv` | 5 | inverse-folding self-consistency, three variants + two entropy-centred |
| `rel` | 27 | S9-3's set-relative block + **per-candidate tournament aggregates** (Borda over the pool, and inside the shipped filter) |
| `seq` | 8 | **NEW** - the candidate window's own sequence against the target's: identity, BLOSUM, four property distances, hydropathy-profile correlation |
| `ctx` | 13 | target-level context including an ESM summary (constant within a target, so it can only gate) |

**Two channels are documented gaps, both priced rather than excused.**

* **AMBER.** `s8/integrate_amber/` covers 19 targets x the first 32 BLOSUM candidates - 608
  of 63,000 pool members. `stage_amber`: **11.9 of those 32 are in-band**, and on the 14
  targets with enough of them **Amber's IN-BAND rho is +0.0645** against the shipped score's
  +0.131 (its rho over the whole cached prefix is +0.4122 - the same collapse everything
  else shows). S8-9's fusion identity says gain goes as the SQUARE of the weaker channel's
  skill, so a pool-wide Amber column predicts ~+0.004 of rank correlation. Extending it to
  the band is ~40 h of box time for a channel the arithmetic already prices at nothing.
  Legacy - S8-9 measured it at +0.62 error-correlated with Amber, "largely one channel" - is
  in the matrix in full.
* **ESM.** The pool's window sequences are **65/500 covered** by `esm_cache.npz`;
  per-candidate ESM means ESM-2 650M over ~63,000 sequences. ESM enters twice anyway: the
  distogram is ESM-conditioned, so `dg` and `dgs` are its downstream product, and the
  target-level PCA profile is in `ctx`. The per-candidate sequence channel is `seq`, which
  MEMORY records as tying ESM at distance MAE on peptides this short - and `seq` alone
  scores **+0.021**, so the channel is not what is missing.

## The instrument validation, run FIRST

One deliberately-leaked column - the candidate's true RMSD percentile - added, pipeline
otherwise unchanged:

| arm | in-band rho | rg-partialled | recall | selected |
|---|---|---|---|---|
| **gbt_l + LEAK** | **+0.9559** | +0.9526 | 0.937 | 1.775 |
| **ridge@10 + LEAK** | **+0.8625** | +0.8748 | 0.950 | 2.117 |
| gbt_l, no leak | +0.1410 | +0.1533 | 0.235 | - |
| ridge@10, no leak | +0.1956 | +0.1943 | 0.241 | - |
| SHIPPED score | +0.1308 | +0.2125 | 0.328 | 3.454 |

**VERDICT PASS.** The harness recovers a signal that is definitionally present, so its
negatives mean something. The leaked arm is fenced behind `LEAK_ARM` and never enters
another stage.

## S9-6a. The capacity curve is MONOTONE DECREASING. The ceiling is informational.

Every arm trained on in-band rows only, leave-fold-out, each with a label-permutation twin
(RMSD shuffled within each target):

| arm | in-band rho | permuted floor | clears? | rg-partialled | +folds | recall |
|---|---|---|---|---|---|---|
| **ridge@100** | **+0.1960** | -0.0015 | yes | +0.1938 | 5/5 | 0.248 |
| ridge@10 | +0.1956 | -0.0001 | yes | +0.1943 | 5/5 | 0.241 |
| ridge@1 | +0.1918 | -0.0007 | yes | +0.1854 | 5/5 | 0.240 |
| pair_lin (RankNet) | +0.1913 | -0.0066 | yes | +0.1864 | 5/5 | 0.232 |
| gbt_s (200x15) | +0.1656 | +0.0121 | yes | +0.1803 | 5/5 | 0.249 |
| gbt_xl (1500x127) | +0.1427 | +0.0242 | yes | +0.1641 | 5/5 | 0.239 |
| gbt_l (600x63) | +0.1410 | +0.0336 | yes | +0.1533 | 5/5 | 0.235 |
| mlp_l (512-512-512) | +0.1190 | +0.0164 | yes | +0.1251 | 5/5 | 0.206 |
| pair_mlp (256-256) | +0.1110 | +0.0068 | yes | +0.1182 | 5/5 | 0.200 |
| mlp_s (64-64) | +0.1031 | +0.0065 | yes | +0.1147 | 5/5 | 0.190 |
| mlp_m (256-256) | +0.0791 | +0.0446 | **no** | +0.0896 | 5/5 | 0.217 |
| *SHIPPED score* | *+0.1308* | - | - | *+0.2125* | 4/5 | *0.328* |
| *const (control)* | *+0.0000* | *+0.0000* | - | *+0.0000* | 0/5 | *0.175* |
| *ORACLE rr* | *+1.0000* | - | - | *+1.0000* | 5/5 | *1.000* |

**The most heavily regularised ridge is the best model in the table.** It beats a
512-512-512 net, a 1500-tree 127-leaf GBT, and both pairwise arms. The curve does not
saturate and then plateau - it saturates at the lowest capacity tried and then **falls**,
with the permutation floor rising exactly as it falls (ridge -0.002, gbt_s +0.012, gbt_l
+0.034, mlp_m +0.045): what extra capacity buys is overfitting, not skill.

**This is the key diagnostic the brief asked for and it lands unambiguously: the ceiling is
INFORMATIONAL, not architectural.** `pair_lin` ties ridge to 0.005, so the regression
objective is not the limit either - a pairwise ranking loss on in-band pairs, the exact form
S9-3's RelNet used, reaches the same place.

## S9-6b. And the rg-partialled version reverses the sign of the result

Raw, the kitchen sink beats the shipped score: **+0.1960 against +0.1308**, a real +0.065.
With radius of gyration partialled out it is **+0.1938 against the score's +0.2125 - it
loses.** Every single arm in the table has a rg-partialled rho below the shipped score's.

**This is the fifth time in this project that an apparent recognition signal has been
compactness** (the Amber correlation, the agreement signal, Legacy in-band, S9-3's
relational net, and now every arm here). The brief called it "the fourth-most-likely place
for a fifth" and it was.

The sharp form of the control - retrain with all seven rg columns **removed** rather than
partialling rg out of a model that uses it - gives **+0.1997 / +0.1974**, essentially
unchanged. So it is not the literal `rg` column: **whatever the model adds over the score
is rg-shaped by whichever route it takes.** Removing the column does not remove the
information, because 104 other columns carry it.

## S9-6c. The whole is WORSE than its best part

`stage_ablate`, ridge@100 (the capacity winner), each block alone and all-minus-each:

| block | ALONE rho | prg | recall | WITHOUT rho | prg | marginal |
|---|---|---|---|---|---|---|
| **rel** | **+0.2244** | +0.2129 | 0.265 | +0.1689 | +0.1673 | **+0.0271** |
| ph | +0.1744 | +0.1587 | 0.163 | +0.1990 | +0.1939 | -0.0030 |
| lg | +0.1702 | +0.1430 | 0.175 | +0.1963 | +0.1949 | -0.0003 |
| dgs | +0.1548 | +0.1368 | 0.220 | **+0.2079** | +0.2071 | **-0.0119** |
| dg | +0.1174 | +0.1099 | 0.234 | +0.1993 | +0.1985 | -0.0033 |
| iv | +0.0887 | +0.1038 | 0.224 | +0.1871 | +0.1838 | +0.0089 |
| seq | +0.0209 | +0.0258 | 0.168 | +0.1958 | +0.1935 | +0.0002 |
| ctx | +0.0000 | +0.0000 | 0.175 | +0.1982 | +0.1950 | -0.0022 |
| **all blocks** | **+0.1960** | +0.1938 | 0.248 | | | |

**`rel` alone (+0.2244) beats every channel combined (+0.1960).** Six of the eight blocks
have a NEGATIVE marginal - adding them makes the model worse. Dropping the per-shell
distogram decomposition, the study's own new channel, is the single largest improvement
available (+0.0119).

That is precisely what S8-9's error-correlation structure predicts. It is also the fourth
independent confirmation that **the only thing with in-band skill is set-relative
centrality** - S8-8 found it as typicality, S8-9 as the consensus criterion, S9-3 as the
tournament statistics, and here as the one block that survives ablation.

Two channel-level negatives worth keeping. **The per-shell decomposition does not rescue the
distogram**: `dgs` alone (+0.155) beats `dg` alone (+0.117), so the shell structure *is*
informative relative to the aggregate score - and it still adds nothing to the whole, and
subtracts when present. **Sequence is measured out at +0.021**, which retires "the sequence
channel was never given to a discriminator at candidate resolution".

## S9-6d. What band training costs, priced rather than assumed

A band-trained model has never seen an out-of-band row, so its whole-pool order is
meaningless and its whole-pool recall is not a fair reading. `stage_global` runs the
identical models on the whole pool, and `compose()` gives the deployable form S8-8
established - the score picks a coarse top-Q, the model re-ranks it (at Q = 75 the kept set
is provably unchanged, which is why the comparison is made at Q > 75):

| arm | in-band rho | rg-part | recall | d recall | predicted A | selected |
|---|---|---|---|---|---|---|
| SHIPPED score | +0.1308 | +0.2125 | 0.328 | - | - | 3.454 |
| ridge@100, band-trained | +0.1960 | +0.1938 | 0.248 | -0.080 | +0.131 | 6.106 |
| ridge@100, whole-pool | +0.0989 | +0.1276 | 0.336 | +0.008 | -0.013 | 3.431 |
| ridge@100 band @Q150 | +0.1264 | +0.1840 | 0.337 | +0.009 | -0.015 | 3.559 |
| **ridge@100 pool @Q150** | +0.1098 | +0.1773 | **0.344** | **+0.016** | **-0.026** | 3.414 |
| gbt_s pool @Q150 | +0.1136 | +0.1875 | 0.338 | +0.010 | -0.016 | 3.312 |

Band training buys in-band rho (+0.196 vs +0.099) and costs whole-pool ordering, exactly as
the brief anticipated. **But the best in-band rho and the best recall are in different arms,
and neither converts.** The best recall any deployable form reaches is **0.344, a gain of
+0.016** - below S9-3's +0.030, and **1/42 of the +0.672 the oracle filter needs.**

## The four numbers

**1. In-band Spearman rho, leave-fold-out: +0.1960 +- 0.0292 raw, +0.1938 rg-partialled**
(best arm `ridge@100`; permutation floor -0.0015 / -0.0069). Against benchmarks: shipped
+0.126-0.131, relational net +0.147 to +0.210 with a rg-partialled +0.20 to +0.23 that ties
the score's +0.212. **The kitchen sink lands inside the relational net's range and, after
partialling, below the shipped score.** It does not approach the +0.3 the brief set as the
threshold for reopening.

Per-target, not just the mean:

| | p10 | p25 | median | p75 | p90 | frac > 0 |
|---|---|---|---|---|---|---|
| ridge@100 | -0.225 | -0.013 | +0.233 | +0.444 | +0.622 | 0.73 |
| SHIPPED score | -0.327 | -0.104 | +0.153 | +0.388 | +0.578 | 0.66 |

And it is **highly concentrated**: the ten largest per-target rho gains carry **127% of the
aggregate** (i.e. the other 116 targets are net negative), the same shape S9-1 found in the
incumbent's own 61%.

**2. In-band recall gain: +0.0157**, from 0.328 to 0.344, in the best deployable form.
**Predicted RMSD gain at S9-3's exchange rate (1.106 A per 0.672 recall = 1.646 A per unit):
-0.026 A.** *This prediction was made and written down before the RMSD was measured.*

**3. Native percentile: 39.45 +- 3.07, median 28.8, argmin at native 6/126.** The shipped
distogram reproduces on the identical augmented pools at 36.76. Benchmarks: distogram 36.0,
inverse folding 30.1, relational 35.5, Amber 54.3, consensus 82.8. **The kitchen sink is
worse than the distogram it contains and worse than inverse folding. Nothing is materially
below 30. No better-specified objective.**

**4. Selected CA-RMSD through the validated architecture** (the model as the FILTER feeding
S8-11's synthesis: coordinate-average the kept 75, project onto ideal geometry):

| filter | kept_mean | recall | synth fit | SE | median | <2 A | d vs the 3.200 incumbent | W/L |
|---|---|---|---|---|---|---|---|---|
| *ORACLE* | *2.686* | *1.000* | ***2.095*** | *0.094* | *2.001* | *0.500* | *-1.106 [-1.289, -0.922]* | *123/3* |
| **score (incumbent)** | 3.551 | 0.328 | **3.200** | 0.154 | 2.992 | 0.286 | - | - |
| ridge@100 @Q150 | 3.539 | 0.337 | **3.228** | 0.152 | 3.047 | 0.254 | **+0.028 [-0.049, +0.105]** | 56/70 |
| gbt_s @Q150 | 3.546 | 0.332 | 3.231 | 0.152 | 3.017 | 0.262 | +0.031 [-0.029, +0.091] | 53/73 |
| ridge@100 pool @Q150 | 3.485 | 0.344 | 3.236 | 0.158 | 3.045 | 0.278 | +0.036 [-0.022, +0.093] | 59/67 |
| ridge@100 pool | 3.488 | 0.336 | 3.245 | 0.158 | 3.044 | 0.294 | +0.045 [-0.025, +0.114] | 56/70 |
| ridge@100 (band, no filter) | 4.391 | 0.248 | 3.625 | 0.156 | 3.475 | 0.214 | +0.425 [+0.249, +0.600] | 47/79 |

**Best deployable arm 3.2284, +0.028 [-0.049, +0.105] against the incumbent, 56W/70L; -0.226
against the 3.454 baseline, which is the incumbent's gain and not a new one.** The predicted
-0.026 A and the measured +0.028 A differ by 0.054 A against a paired SE of 0.039 - the
exchange rate is, if anything, generous to the model.

Concentration on the RMSD side too, and it is the S9-1 shape: dropping the ten targets on
which `ridge@100@Q150` helps most leaves it at **+0.107 A worse than the incumbent**, against
its headline +0.028. The arm's whole apparent parity is ten targets.

Note the row that settles the mechanism: `ridge@100 pool @Q150` hands the operator a
**strictly better set** (kept_mean 3.485 vs the score's 3.551, recall 0.344 vs 0.328) and
returns a **worse** answer. S9-3 measured the same inversion with `cty_med75`. At this scale
filter quality does not transfer at all.

## Controls

* **Label permutation**, every arm, shuffled within target: floor -0.007 to +0.045, rising
  monotonically with capacity. Every arm except `mlp_m` clears its own floor by more than 2
  SE. `const` scores exactly 0.000 under the corrected tie convention (S8-8's `sel_of`
  defect, re-checked here as a test: a constant signal must cost exactly the pool mean, and
  it does to 1e-6).
* **Feature ablation**, above: the whole is worse than its best part.
* **Oracle-feature control**, run first: PASS at +0.956.
* **Compactness**, two ways: partialled, and retrained with the columns removed.
* **Leakage**: 16 targets x 105 features, `rr` and `nat_ca` NaN-poisoned, **worst |diff|
  0.000e+00**, and no NaN propagates. The band index and the label DO read `rr` - they are
  supervision, not features - and are named in the output as excluded by construction rather
  than silently.
* **Leave-fold-out** asserted as a test, plus a `const` arm that must return zeros.
* **The harness finds planted signals**: ridge recovers a linear one at rho > 0.95, GBT a
  nonlinear one at > 0.7.

34 checks, 34 passing. One test failed during development and **the feature was right**: it
planted a point at the origin and called it the pool medoid, when a finite Gaussian sample's
medoid is the argmin of mean pairwise distance and was a different point.

## The verdict

**In-band discrimination is informationally impossible with the features this project can
compute, not merely unattempted at scale.**

The strong form of the statement is earned by three things that had never been done together:

1. **Capacity was swept and the curve goes the wrong way.** A ridge with the heaviest
   regularisation tried beats every tree and every net. If the limit were architectural, more
   capacity would buy something; it buys overfitting, measured against its own permutation
   floor. There is nothing left to try in model class.
2. **Every channel was present simultaneously and the combination is worse than its best
   member.** Six of eight blocks have a negative marginal. This is not "we tried the wrong
   pair" - the pairwise fusion arithmetic S8-9 wrote down predicted it, and it is now
   confirmed at full width.
3. **The one apparent gain is compactness for the fifth time**, and it survives removing the
   rg columns, which means it is not a feature to delete but a property of what the pool
   varies in.

The recognition line should stop. What remains measurable is on the other side of the ledger:
the oracle filter is worth 1.106 A and the whole in-band signal available to buy it is 1/42
of the recall it needs. **The 1.106 A is not going to be recovered by ranking.** S9-4's one
positive - that a better *query* moves universe recall 0.038 -> 0.160 and returns 1.925 A -
remains the only measured route past 2.0 A on this instrument, and it is a GENERATION result,
not a recognition one.

Two things to carry forward that are not the headline:

* **`rel` alone is a better in-band ranker than every channel combined** (+0.224 vs +0.196).
  If anyone re-opens this, the honest starting point is set-relative centrality alone, not a
  fusion.
* **The exchange rate held.** A prediction of -0.026 A made from recall before the RMSD was
  measured came within 0.054 A of the +0.028 A observed. S9-3's pricing rule is now validated
  out-of-sample and should be used to kill future selection arms before they are built.


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / evo

> *Working paper, merged verbatim from `s9/evo_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-EVO. Is there another information source for a 9-16mer? Yes - and it is not structural

`s9/evo.py`, `s9/test_evo.py` (35 checks, all passing), `s9/evo_*.json`.
Instrument: the 126-target tuning set, reproducing sprint 7/8 exactly (shipped score
**3.4540**, pool best **1.7108**, S8-8 medoid75 **3.2822**), asserted as a test before
anything else ran. Dev(24) and benchmark(60) were not read.

## The question

Eight sprints consumed one input: the target's own sequence. The project's central
conclusion - that sequence alone is insufficient at this length - was an *inference* from
that, never a comparison against an alternative. This study goes looking for the
alternative, in the two places it could exist: **structural homologs** in `prots/`, and
**evolutionary profiles** from public sequence databases.

The answers run in opposite directions, and both are new.

---

## 1. Structural homolog retrieval is dead at this length, and a shuffled sequence proves it

`s9/evo_scan.json`, `s9/evo_null.json`. The whole of `prots/` was parsed to SEQRES -
**13,751 files → 20,929 chains → 3,519,987 residues** - and every 9-16mer window of it was
scored against every instrument target. `fragment_db` samples this corpus at stride 5;
this is the first exhaustive pass over it.

### The availability numbers

| | count |
|---|---|
| targets whose own PDB entry is in the corpus | **0 / 126** |
| targets with an exact copy elsewhere | 1 / 126 |
| targets with a window at identity >= 0.8 elsewhere | 2 / 126 |
| targets with a window at identity >= 0.6 (positional) | 27 / 126 |
| **targets with a window at >= 0.6 in the project's own gapped identity** | **67 / 126** |
| windows at >= 0.6, gapped, over the instrument | 3,190 |
| mean / median / max best gapped identity | **0.602 / 0.600 / 1.000** |

**The median target's best available structural match sits exactly on the leakage
threshold.** Every one of those 3,190 windows is excluded by `peptide_db.folds(5)` and the
0.6 filter, which is the correct behaviour and is the point: *at this length, homolog
detection and leakage detection are the same operation.* A window similar enough to be a
homolog is similar enough to be the answer.

### The control that says the rest is not homology at all

A search over 3.5 M windows for the best 9-16mer match returns a high identity by EXTREME
VALUE whether or not anything homologous exists - and the length dependence gives it away
(mean best positional identity **0.630 at L=9 falling to 0.461 at L=16**, i.e. the shorter
the peptide the "better" its best match). So the same search was run with each target's own
residues SHUFFLED, which destroys homology and preserves length and composition exactly.

| statistic | real sequence | shuffled | difference | 95% CI | W/L |
|---|---|---|---|---|---|
| positional best identity | 0.532 | 0.523 | **+0.008** | [-0.003, +0.020] | 44 / 54 |
| **gapped best identity** (the leakage unit) | **0.602** | **0.593** | **+0.008** | **[-0.002, +0.018]** | 53 / 46 |
| **windows at >= 0.6** | 25.3 | **28.7** | **-3.4** | [-11.2, +4.4] | - |

**A shuffled version of the target finds a match in 13,751 protein structures that is just
as good as the real sequence does**, in both units, and it crosses the leakage threshold
*more often* - 28.7 windows against 25.3. Only 12/126 real sequences beat all five of their
own shuffles. Whatever the search returns, it is reading the size of the database.

**Verdict: there is no structural homolog to retrieve.** The 67 targets with a >=0.6 window
are the ones the leakage filter exists to remove, and everything below that threshold is
indistinguishable from a composition-matched random string. This is a demonstrated fact
about the available data, not an inference.

---

## 2. An evolutionary profile IS retrievable - for 55% of the instrument

`s9/evo_pepsearch.json`, `s9/evo_profile.json`. This is the surprise, and it goes against
the brief's own expectation that "for 9-16mers a conventional MSA is hard".

The chain is: the target's sequence → UniProtKB exact peptide search → parent proteins →
each parent's UniRef50 cluster → a global alignment of every member to the parent → the
columns that align to the peptide. Sequence databases only; no structure is requested at
any step.

| | count |
|---|---|
| **targets whose sequence occurs verbatim inside a known protein** | **71 / 126** |
| median number of parent proteins when it does | 12 |
| **targets with at least one aligned homologous row** | **69 / 126** |
| targets reaching effective depth Neff >= 2 | **46** |
| targets reaching Neff >= 5 | 31 |
| targets reaching Neff >= 10 | 19 (max 134) |
| median rows / median Neff when found | 12 / 4.58 |
| mean profile column entropy | **1.592 nats** (one-hot 0.000, uniform 2.996) |
| mean identity of an aligned row to the target | 0.735 |

**The premise "a 12-mer has no MSA" is false for a majority of this instrument.** The
information source the project has never had does exist, at real depth, for 46 targets.

Neff is Henikoff-weighted, so a profile whose members are all >= 80% identical to the target
counts as ONE effective sequence - which is the honest measure, because such a profile *is*
the target's own one-hot however many rows it has.

---

## 3. The precondition: it is a different key, and the control says how different

`s9/evo_diverge.json`. Before any accuracy number: does a profile actually rank the window
universe differently from a plain substitution matrix?

The primary key is **profile-weighted BLOSUM**, `sum_i sum_a f[i,a] * B62[a, w_i]`, chosen
because of its degenerate case - with an empty MSA the profile is a one-hot and the key
reduces to `blosum_key` **exactly**. Verified on the real instrument, not merely
synthetically: on the **57 targets with no homologous row the two keys agree to 0.0 and
produce the identical filter order**, and `_assert_matched_arms` fails the stage if that
ever stops being true. So the arms are matched on scoring form, on normalisation, and on
every target where evolution has nothing to say.

| comparison | rho (all 126) | rho (Neff >= 2, n=46) |
|---|---|---|
| profile-BLOSUM vs plain BLOSUM | 0.956 | **0.882** |
| profile-BLOSUM vs its own empty-MSA baseline | - | **0.882** (top-500 overlap 0.630) |
| *PSI-BLAST log-odds vs plain BLOSUM* | *0.638* | *0.583* |
| **one-hot log-odds vs plain BLOSUM (CONTROL)** | **0.675** | - |

The last two rows are why the log-odds form is only a secondary arm. It correlates with
BLOSUM at 0.638, and the SAME form with no homologues at all correlates at 0.675 - so
"the profile ranks differently" scored that way is a statement about arithmetic, not about
evolution. The matched key moves the ranking modestly but really: on the 46 deep targets it
changes **38% of the top 500**.

### Five defects found and fixed, every one of which would have produced a publishable number

1. **The alphabet.** `s7.audit` encodes sequences in `ARNDCQEGHILKMFPSTWYV`; this module
   used `ACDEFGHIKLMNPQRSTVWY`. The cached window banks are in `audit`'s order, so a profile
   built here and indexed with a cached window scored the wrong amino acid at every
   position, silently. One filter run and one divergence table were computed that way.
   Pinned by `t_alphabets_are_reconciled`.
2. **The pseudocount, in the matched key.** Caught by an internal consistency check, not a
   test: profile-BLOSUM differed from plain BLOSUM by up to **0.61 A of pool best** on
   targets where it should have differed by zero. With a pseudocount an empty profile is a
   *smoothed* one-hot, so the arms differed everywhere for a reason unrelated to evolution.
   `PSEUDO_KEY = 0`, plus two tests, one through `load_profile` on a real empty-MSA target
   and its real cached window bank.
3. **The pseudocount, in the OTHER key.** Fixing (2) globally broke the log-odds arm: a
   log-odds is undefined at a zero frequency, 19 of 20 columns become `-inf`, every window
   without a perfect match ties, and `argsort` returns the incoming pool order instead of
   the key. So the two forms carry different constants by construction - `PSEUDO_KEY = 0`
   for profile-BLOSUM, which needs none and therefore reduces exactly, and `PSEUDO_LO = 1`
   for the log-odds, which needs one - plus a floor inside `pssm_of` and
   `t_logodds_needs_a_pseudocount`.
4. **A cache shared between two key definitions.** Two chained runs overlapped; the second
   silently reused the first's cached cells because both produced the same cell *names*,
   and the published table was a mix. The cache path now carries a `FILT_VERSION`, so a key
   change makes old cells unreachable rather than reusable.
5. **A control that was a different control on every run.** The random filter arm was
   seeded from `hash(pdb)`, and Python randomises string hashing per process - two runs
   gave 3.619 and 3.658. Now `zlib.crc32`, with a test. (The two runs agreeing to 0.04 A is
   an accidental but welcome replication of the control.)

And one conceptual correction that matters more than any of them: **the matched-arm
predicate is "the MSA is empty", not "Neff = 1"**. Henikoff weighting counts effective
*clusters*, so a target whose homologues are all >= 80% identical to each other has Neff = 1
while its profile still differs from a one-hot wherever those homologues differ - 2RUO is
Neff 1.0 with 0.75 in four columns. Under the wrong predicate the invariant appeared to
fail on 74 targets; under the right one it holds at exactly 0.0 on all 57.

Every table below is from the code after all five.

---

## 4. The FILTER role - where the 1.07 A is, and the profile does not touch it

`s9/evo_filter.json`. 126 paired targets, identical pools, identical synthesis operator
(coordinate-average the filtered subset, project onto ideal geometry); **only the ORDER
differs.** The instrument reproduces S8-11 (`sc|75` = 3.200 against 3.204, the 0.004 A
being BLOSUM tie-break noise) and reproduces the headroom S8-13 measured:

| arm | selected | d vs sc\|75 | 95% CI | W/L | subset best | <2 A |
|---|---|---|---|---|---|---|
| *O_oracle\|25* (DIAGNOSTIC) | *1.769* | *-1.431* | *[-1.634, -1.228]* | *126/0* | *1.711* | *0.659* |
| *O_oracle\|75* (DIAGNOSTIC) | *2.095* | *-1.106* | *[-1.291, -0.921]* | *123/3* | *1.711* | *0.500* |
| **sc\|75 (incumbent)** | **3.200** | - | - | - | 2.306 | 0.286 |
| blos_sc\|25 | 3.215 | +0.014 | [-0.047, +0.076] | 58/68 | 2.480 | 0.294 |
| **pblos_sc\|25** | 3.216 | +0.016 | [-0.047, +0.078] | 55/71 | 2.481 | 0.294 |
| sc\|25 | 3.232 | +0.032 | [-0.026, +0.090] | 55/71 | 2.609 | 0.286 |
| pssm_sc\|75 | 3.241 | +0.040 | [-0.030, +0.110] | 55/71 | 2.156 | 0.270 |
| blos_sc\|75 | 3.244 | +0.043 | [-0.023, +0.109] | 58/68 | 2.155 | 0.262 |
| **pblos_sc\|75 (profile + distogram)** | **3.249** | **+0.048** | **[-0.017, +0.113]** | 57/69 | 2.169 | 0.270 |
| blos\|75 (plain BLOSUM alone) | 3.482 | +0.282 | [+0.120, +0.443] | 45/81 | 2.104 | 0.254 |
| **pblos\|75 (profile alone)** | 3.499 | +0.298 | [+0.138, +0.459] | 43/83 | 2.105 | 0.238 |
| pssm\|75 | 3.509 | +0.309 | [+0.131, +0.487] | 45/81 | **2.057** | 0.222 |
| **pblos\|25** | 3.519 | +0.318 | [+0.145, +0.492] | 43/83 | 2.364 | 0.214 |
| blos\|25 | 3.533 | +0.332 | [+0.154, +0.510] | 46/80 | 2.350 | 0.198 |
| rand\|75 (control) | 3.658 | +0.458 | [+0.276, +0.639] | 41/85 | 2.072 | 0.230 |
| rand\|25 (control) | 3.716 | +0.515 | [+0.328, +0.702] | 39/87 | 2.379 | 0.190 |

**Not one profile arm beats the distogram filter, and the profile does not beat the plain
substitution matrix it was built to improve on** - it is 0.014 A better at m=25 and 0.017 A
worse at m=75, i.e. a coin flip inside a 0.046 A paired SE.

On the pre-registered **Neff >= 2 subgroup** (46 targets, the only ones where a profile can
differ from a one-hot at all), the same nullity with the same sign:

| arm | selected | d vs sc\|75 | 95% CI | W/L |
|---|---|---|---|---|
| *O_oracle\|75* | *2.178* | *-1.041* | *[-1.322, -0.759]* | *46/0* |
| sc\|75 | 3.218 | - | - | - |
| pblos_sc\|75 | 3.282 | +0.064 | [-0.027, +0.154] | 20/26 |
| blos\|75 | 3.648 | +0.430 | [+0.137, +0.723] | 18/28 |
| pblos\|75 | 3.687 | +0.468 | [+0.186, +0.751] | 16/30 |
| **rand\|75 (control)** | **3.805** | +0.587 | [+0.281, +0.892] | 15/31 |

### Why it cannot pay, stated as a number

**Sequence keys of ANY kind are barely better than random as filters.** Over all 126:
BLOSUM 3.482, profile-BLOSUM 3.499, PSI-BLAST 3.509 - against random's 3.619 and the
distogram's 3.200. On the deep subgroup the whole sequence channel collapses to **0.074 A**
(BLOSUM 3.648 against random 3.722). The dynamic range inside which the profile has to find
a gain is smaller than the paired SE of the instrument, while the gap to the oracle is
1.11 A. **Improving the sequence key cannot matter because the sequence key is not doing
the work.**

**And the S8-6 signature repeats exactly.** `pssm|75` keeps a strictly BETTER subset - best
member 2.057 A against the incumbent's 2.306 - and returns a WORSE answer (3.509 vs 3.200).
More near-native candidates, worse output, for the third time in two sprints.

---

## 5. The RETRIEVAL role - the profile builds the best pool ever measured here, and it does not help

`s9/evo_retr.json`. The profile as a retrieval key over the WHOLE window universe (median
17k, max 39k windows), K=500. S8-6 screened 22 keys and none was a profile; this is the one
classical retrieval advance the project had never tried. Instrument check exact: `blosum`
reproduces **1.711 / 1.994 / 3.454**.

| key | pool best | ORACLE ceiling | selected | synthesis | d(synth) | 95% CI |
|---|---|---|---|---|---|---|
| **blosum (incumbent)** | 1.711 | 1.994 | 3.454 | **3.200** | - | - |
| **pssm** | **1.697** | **1.963** | 3.511 | 3.266 | **+0.066** | [+0.013, +0.118] |
| fuse (rank sum) | 1.720 | 2.015 | 3.436 | 3.229 | +0.028 | [+0.004, +0.053] |
| pblos | 1.732 | 2.034 | 3.435 | 3.237 | +0.036 | [+0.005, +0.068] |

**The PSI-BLAST profile key has the best perfect-scorer ceiling of any retrieval key ever
measured on this instrument - 1.963 A, against BLOSUM's 1.994 and S8-6's best of 22 keys
(`ident`, 1.983).** It also improves pool best (1.697 vs 1.711). And its synthesis output is
**significantly WORSE** (+0.066, CI excluding zero).

That is S8-6's law in its cleanest form yet: the profile is a genuinely better *retrieval*
key by every pool statistic and a worse *system*. All three profile arms are significantly
worse on the deliverable, and on the Neff >= 2 subgroup every one is nominally worse again
(pblos +0.039, pssm +0.050, fuse +0.022, all CIs spanning zero).

---

## 6. The CONDITIONING role - the gain is larger where the profile is EMPTY

`s9/evo_cond.json`. Leave-fold-out over `peptide_db.folds(5)`, peptides only, arms
differing ONLY in the feature block - the design S7-11 used to show ESM beats one-hot by
0.288 A. **423/787 training peptides have a profile.** Judged on selected CA-RMSD through
the same synthesis operator, never on MAE, for the reason S7-11 and S8-14 both record.

The profile block is 63 numbers per pair: the two columns' frequency vectors, their
product, the two column entropies, and the profile's effective depth. `shuf` is the same
block with its columns permuted along the sequence - every marginal preserved, the
correspondence between a column and the residue it describes destroyed.

| arm | synthesis | top-1 | d vs base | 95% CI | W/L | sign p |
|---|---|---|---|---|---|---|
| base | 3.336 | 3.451 | - | - | - | - |
| **prof** | **3.237** | 3.447 | **-0.099** | [-0.197, +0.000] | 73/53 | 0.090 |
| shuf (control) | 3.300 | 3.546 | -0.035 | [-0.130, +0.060] | 68/58 | 0.423 |
| **prof - shuf (THE TEST)** | | | **-0.063** | **[-0.147, +0.020]** | **62/64** | **0.929** |

Against its own control the profile is null: 62 wins, 64 losses.

### The decomposition that kills the mechanism

| subgroup | n | prof - base | 95% CI |
|---|---|---|---|
| **Neff >= 2** (real evolutionary depth) | 46 | **-0.077** | [-0.250, +0.095] |
| **Neff == 1** (the profile IS the target's one-hot) | 80 | **-0.111** | [-0.233, +0.012] |

**The nominal gain is LARGER on the 80 targets where the profile can contain no
evolutionary information at all**, and smaller on the 46 where it has the most. A real
mechanism scales with exposure; this one runs backwards. What the block buys is capacity
and a redundant re-encoding of the sequence, which is what `shuf` also buys and what the
Neff = 1 rows confirm independently. This is S8-14's dilution control run in reverse, with
the same verdict.

---

## 7. Leakage audit - the headline, not the housekeeping

`s9/evo_leak.json`, `s9/test_evo.py`. This was the highest-leakage task in the sprint: it
deliberately went looking for other proteins that resemble the target, so the discipline
had to be mechanical.

### The inference-time availability rule, stated in the module docstring before use

At inference on a held-out target the deployable path may read the target's **sequence**;
**UniProtKB and UniRef50** queried with that string (a sequence record is not a structure);
and the fold-filtered structural library the project already uses. It may **not** read the
target's own deposited structure or any NMR model of it, **any structure at >= 0.6 identity
whatever its accession**, or any entry whose inclusion depends on knowing the answer.

### The audit

| check | result |
|---|---|
| **bit-identical deployable quantities under NaN-poisoning of every native array** | **1,890 / 1,890** |
| worst change in any retrieval key when the native is NaN | **0.0** |
| targets whose own PDB entry appears in the 20,929-chain corpus | **0 / 126** |
| library identity violations (`peptide_db.folds(5)`, threshold 0.6) | **0** |
| targets with a >= 0.6 window anywhere in `prots/` (gapped identity) | 67 / 126 |
| such windows in total - **all excluded by the filter** | 3,190 |
| worst identity surviving the filter | 0.583 (positional), reported per target |
| profile rows are sequences only - checked character by character | yes |

Every filter order, every retrieval key and every synthesised structure was recomputed with
`nat_ca`, `rr` and `rr_reb` replaced by NaN and asserted **bit-identical**; the profile and
contact-map builders are asserted unable to *accept* a native in their signatures; and the
cached MSA JSON is asserted to contain only amino-acid letters and gap characters.

### The point worth keeping

**A homolog that is really the same peptide is leakage even under a different accession**,
and at this length that is nearly the whole population of homologs. 67/126 targets have a
window in `prots/` at >= 0.6 gapped identity and every one is removed - while the shuffled
control shows the *remaining* matches are not homologs at all. The leakage filter and the
homolog detector are the same instrument pointed in opposite directions.

---

## 8. The direct answer

**Was sequence alone genuinely the limit, or merely the only thing tried?**

Both halves have now been measured, and they answer differently.

* **It was not merely the only thing tried - a second source exists.** An evolutionary
  profile is retrievable for 69/126 targets and reaches real depth (Neff >= 2) on 46. The
  brief's own expectation that "a conventional MSA is hard" for a 9-16mer is wrong: a
  majority of these peptides are excised from proteins that have large, well-populated
  UniRef50 clusters, and the peptide's own columns can be read straight out of them.
* **And it buys nothing, in every role the brief mandated.** As a FILTER it is +0.053
  [-0.016, +0.121] against the distogram and ties plain BLOSUM. As a RETRIEVAL key it
  improves the pool's perfect-scorer ceiling to the best figure ever measured here (1.963
  against 1.994, beating all 22 keys S8-6 screened) and makes the deliverable
  significantly WORSE (+0.066 [+0.013, +0.118]). As CONDITIONING it is -0.063 [-0.147,
  +0.020] against its own shuffled control, with the effect *larger* where the profile is
  empty.
* **Structural homology, meanwhile, does not exist at this length.** A shuffled sequence
  finds as good a match in 13,751 structures as the real one (+0.005 [-0.007, +0.017]), and
  the 67 targets that do have a >= 0.6 window are precisely the ones the leakage filter
  exists to remove. **Homolog detection and leakage detection are the same operation for a
  12-mer.**

So the honest statement is neither "sequence is the limit" nor "we never tried anything
else". It is: **the limit is not the information channel, it is the consumer.** Every one
of the three new signals moved the pool statistics - the profile key demonstrably retrieves
better structures than BLOSUM - and none of them moved the answer, because the same
recognition defect that S7-5, S7-10, S8-5, S8-6 and S8-13 all measured absorbs whatever is
fed into it. The 1.11 A that a perfect filter is worth (2.095 vs 3.200 here, reproducing
S8-13) remains untouched by better inputs, and this study is the first to demonstrate that
with an input the project did not previously have.

### What was not done, and why

* **ESM contact maps in the parent protein's context** - run ESM-2 on the parent, slice out
  the peptide's block, and use it where the isolated 12-mer's map is used. It is built and
  tested (`stage_esmctx`, a bounded throwaway batch after an unbounded one took the box to
  96.8% RAM) and was deprioritised on the coordinator's instruction: ESM has now been
  measured three ways here - features +0.288 A, similarity metric at exact BLOSUM parity,
  generative conditioning null under the distogram filter - and a fourth variant is not
  where the value lies. The stage is left runnable; `parents` has already cached all 71
  parent sequences, so it is one command away.
* **A learned retriever on the profile.** S8-6 priced perfect retrieval at 0.016-0.070 A at
  K=500 and there is no held-out target set left to train one on; sections 4 and 5 here now
  add that a profile key which *achieves* a better ceiling than any of S8-6's 22 keys still
  makes the deliverable worse.

### What the next sprint should take from this

1. **Stop looking for inputs.** Three genuinely new signals - an exhaustive structural
   search 100x larger than the shipped library, an evolutionary profile at Neff up to 134,
   and a profile retrieval key with a record ceiling - all moved the pool statistics and
   none moved the answer. That is now four independent information channels (sequence, ESM,
   physics, evolution) meeting the same wall.
2. **The one number that still has room is the filter, and it is not an information
   problem.** `O_oracle|75` reaches 2.095 A on pools this pipeline already builds, from
   candidates the distogram already has in front of it. The gap is 1.11 A and every input
   tried lands within 0.05 A of the incumbent.
3. **If anyone re-opens retrieval, the profile key is the one to re-run at K=25**, where
   S8-6 measured perfect retrieval as worth 3.5x more. `pssm` is the only real key ever to
   beat BLOSUM's oracle ceiling, and that regime is the only place S8-6 left the axis alive.


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / loop

> *Working paper, merged verbatim from `s9/loop_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-x. Closing the loop: the retrieval loop works, and the pipeline sits 0.02 A from its break-even

`s9/loop.py` (6 stages), `s9/test_loop.py` (38 checks, 0 failures),
`s9/loop_{repro,run,qcurve,price,audit,leak,natgeom}.json`, `s9/loop_cache/` (per-target
records behind every table). 126-target tuning instrument, K=500, m=75. dev-24 not spent;
the 60-target benchmark not read.

## The gap, and why it was worth attacking

Every stage of this pipeline runs exactly once: retrieve the BLOSUM top-500, filter to the
distogram's top 75, coordinate-average, project. **The output is never fed back into the
input.** It now produces something it did not have before - a synthesised structure at
**3.2005 A that is in no pool** - and a structure is a far better retrieval key than a
sequence: BLOSUM puts the structurally best window at percentile 39.3 against random's
50.0, while the distogram as a key reaches 16.5 (S8-6).

So: retrieve → filter → synthesise → **re-retrieve windows structurally similar to the
synthesis** → re-filter → re-synthesise. Five arms, five rounds each.

Round 0 IS the shipped pipeline and `stage_repro` asserts it on all 126 targets before
anything iterates: **baseline 3.4540 (expect 3.4540), pool best 1.7108 (1.7108), raw
average 3.0483 (3.0483), projected synthesis 3.2005 (3.2005)** - the BLOSUM top-500 slice
of `s8/generate_univ` is byte-identical to `s8/consensus2_cache`'s pool and the recomputed
distogram score matches its cached vector to 0.0, both asserted as tests.

## L-1. The loop converges, and it converges DOWNHILL

| arm | rd0 | rd1 | rd2 | rd3 | rd4 | rd5 | d vs 3.200 incumbent | W/L |
|---|---|---|---|---|---|---|---|---|
| **replace** | 3.201 | 3.246 | 3.284 | 3.307 | 3.319 | **3.323** | **+0.123 [+0.051,+0.194]** | 47/79 |
| **fuse** (BLOSUM+structure rank sum) | 3.201 | 3.233 | 3.254 | 3.263 | 3.269 | **3.270** | **+0.070 [+0.001,+0.139]** | 53/73 |
| **direct** (no score filter) | 3.201 | 3.230 | 3.262 | 3.277 | 3.290 | **3.304** | **+0.103 [+0.013,+0.193]** | 50/76 |

Monotone in every arm, every round, with the CI excluding zero by round 5 in all three.
Against the 3.454 baseline the arms are still ahead (replace -0.131 [-0.269,+0.007], fuse
-0.184 [-0.320,-0.048]) - but that is the incumbent's -0.253 decaying, not a new gain.

**It converges to a fixed point.** Round-to-round filtered-set overlap runs
0.095 → 0.788 → 0.896 → 0.937 → **0.962** (replace), and selected RMSD flattens to
3.319 → 3.323. Damping helps only by slowing the descent: `fuse` overlaps 0.288 at round 1
and loses 0.070 instead of 0.123.

### The primary diagnostic the brief demanded: the input set really did change

S8-6's warning was that injecting the library's best window is worth -0.016 A and injecting
*every* sub-2 A window returns the identical structure on 90/126 targets - i.e. a loop that
only enlarges the pool measures nothing. That is not what happened here.

**Round 1's filtered top-75 overlaps round 0's by 0.095, and the answer changed on
126/126 targets.** The pool overlap is 0.080. This loop replaced 90% of the synthesis input
set and made the answer worse. It is a real intervention with a real negative.

### What the loop does to the pool, and the fifth instance of the same trap

| arm | pool best | pool mean | subset spread |
|---|---|---|---|
| round 0 (BLOSUM) | **1.711** | 4.453 | 2.486 |
| replace round 5 | 2.565 | **3.515** | 1.551 |
| fuse round 5 | 2.231 | 3.587 | 1.857 |
| direct round 5 | 2.917 | **3.366** | 1.058 |

Pool mean improves by up to 1.1 A and pool best degrades by up to 1.2 A - the exact
signature S8-6 measured for the distogram key (1.711 → 2.161), only stronger, because
CA-RMSD to a query is a sharper concentrator than matrix agreement. Evaluating S8-5's
published law `sel = 0.257*pool_best + 0.311*pool_mean + 1.673` on these pools:

| round | pool best | pool mean | law says | actual | residual |
|---|---|---|---|---|---|
| 0 | 1.711 | 4.453 | 3.498 | **3.201** | -0.297 |
| 5 | 2.565 | 3.515 | 3.425 | **3.323** | -0.102 |

**The law predicts the loop should GAIN 0.073 A and it loses 0.123 A.** The residual - the
synthesis operator's advantage over what the pool alone predicts - is itself eaten, falling
from -0.297 to -0.102. Consensus needs a DIVERSE input set; a structurally-retrieved set is
already collapsed onto one mode, so averaging it has nothing left to average away.

## L-2. Recall is the currency, and the loop spends it

| round | recall_pool | recall_univ | selected |
|---|---|---|---|
| 0 | **0.328** | 0.038 | 3.201 |
| 1 | 0.219 | **0.053** | 3.246 |
| 5 | **0.204** | 0.045 | 3.323 |

Two senses, because the loop changes which candidates exist. Pool-relative recall (S9-3's
definition: overlap of the filtered 75 with the pool's own oracle best 75) **falls by
0.124**. Universe-relative recall - overlap with the whole ~18,674-window library's best 75,
the only one a loop can move by construction - rises **+0.015 at round 1 and settles at
+0.007**.

S9-3 priced the exchange: the oracle filter needs **+0.672 of recall to buy 1.106 A**. At
that rate a **-0.124** recall move predicts **+0.203 A**; the measured cost is **+0.123 A**.
Same sign, same order of magnitude. **The loop paid for the recall it destroyed.**

**Direct answer to the brief's question: no. Closing the loop does not move recall. It
moves pool-relative recall backwards by 0.124 and universe-relative recall forwards by
0.007, against the +0.672 that would matter.**

### There is no regime where the deployable loop pays

Stratified by round-0 query quality (round 1 minus round 0, `replace`):

| stratum | n | query | round 1 | d | W/L |
|---|---|---|---|---|---|
| q1 best | 32 | 1.173 | 1.201 | +0.027 [-0.032,+0.086] | 13/19 |
| q2 | 31 | 2.532 | 2.497 | -0.035 [-0.139,+0.070] | 12/19 |
| q3 | 31 | 3.518 | 3.632 | **+0.115 [+0.022,+0.207]** | 10/21 |
| q4 worst | 32 | 5.568 | 5.643 | +0.075 [-0.003,+0.153] | 12/20 |
| query < 3.0 A | 64 | 1.860 | 1.857 | -0.003 [-0.062,+0.056] | 25/39 |

Not one stratum's interval sits below zero - not even the 32 targets whose query was
already 1.17 A. And the loss is broad, not outliers: dropping the ten most-improved targets
takes `replace`'s round-1 delta from +0.046 to **+0.090**.

## L-3. The mechanism is right; the query is 0.018 A too poor. This is the finding.

`stage_qcurve` runs ONE round of structural retrieval from queries of controlled quality.
Three rungs are ORACLE (they read the native) and exist only to price the axis.

| query | query A | selected | pool best | recall_univ | in-band frac | cand angle | d vs incumbent |
|---|---|---|---|---|---|---|---|
| *the native* | *0.000* | ***1.925*** | *1.313* | ***0.160*** | ***0.814*** | *101.77* | ***-1.275 [-1.490,-1.060]*** |
| *library best window* | *1.313* | *2.070* | *1.313* | *0.156* | *0.751* | *101.79* | *-1.130 [-1.331,-0.929]* |
| *pool best member* | *1.711* | *2.210* | *1.349* | *0.135* | *0.703* | *101.75* | *-0.991 [-1.183,-0.799]* |
| **synthesis (DEPLOYABLE)** | **3.201** | **3.246** | 2.470 | 0.053 | 0.430 | 100.60 | **+0.046 [+0.003,+0.089]** |
| rescaled raw average | 4.091 | 3.411 | 2.702 | 0.050 | 0.387 | **103.55** | +0.211 [+0.097,+0.324] |
| pool medoid | 3.658 | 3.621 | 3.162 | 0.016 | 0.368 | 93.84 | +0.420 [+0.206,+0.635] |
| random pool member | 4.544 | 4.068 | 3.317 | 0.014 | 0.223 | 100.07 | +0.867 [+0.597,+1.137] |

**A native query returns 1.925 A - under the sprint's 2.0 A objective, and better than the
perfect-distance-oracle ceiling of 1.994 A that S7 finding 7 and S8-4 put on this library.**
It moves universe recall 0.038 → **0.160** (4.2x) and the in-band fraction 0.430 → **0.814**.
So the mechanism the brief hypothesised is real and large: **retrieval by structure DOES
move recall, by +0.122, when the query is right.** Nothing else in two sprints has moved
recall by more than 0.030.

Every rung is scored on the same 126 targets, so regressing rung means on rung query
quality holds the target set fixed and is causal:

**selected = 0.492 x query + 1.634,  r = +0.9614 over 7 rungs**

**Break-even against not looping (3.201 A) is a query of 3.183 A. The pipeline's own query
is 3.2005 A - 0.018 A short.** The per-target-rung regression (confounded by the fact that
a hard target has both a worse query and a worse answer, so it over-states the slope) puts
break-even at 3.022 A, and S8-3 independently measured that ranking by CA-RMSD to a
reference "is dead by 3 A" (+0.019 [-0.092, +0.126]). **Three estimates, three methods, all
landing between 3.0 and 3.2 A.**

This is not "the loop is wrong". It is **"the loop is premature, and by a measurable
amount."** Every 1 A of query improvement buys 0.49 A through this loop; the first 0.02 A
buys the right to iterate at all.

## L-4. The self-confirmation control: it is a search, and it is path-dependent

Seed round 0 with a deliberately poor structure and run identical dynamics.

| control | seed | final | vs its own seed | vs `replace`'s fixed point |
|---|---|---|---|---|
| `ctl_med` (K=500 pool medoid) | 3.658 | 3.597 | **-0.061 [-0.169,+0.048]**, 64/62 | **+0.274 [+0.044,+0.504]** |
| `ctl_rand` (random pool member) | 4.544 | 3.580 | **-0.964 [-1.220,-0.707]**, 92/34 | **+0.257 [+0.053,+0.461]** |

**It is not a pure fixed-point artefact.** From a 4.54 A random member the loop travels
0.96 A and 92/126 targets improve - it genuinely searches.

**But it is not seed-independent either.** Both controls land ~0.26 A worse than `replace`,
with intervals excluding zero, and comparing the FIXED POINTS THEMSELVES rather than their
errors (equal error is not the same answer):

- CA-RMSD between `ctl_med`'s fixed point and `replace`'s: **1.871 A mean, 0.809 median**,
  agreeing to 0.5 A on **59/126**; their final input sets overlap **0.466**.
- For `ctl_rand`: **1.500 A mean, 0.852 median**, 0.5 A on **56/126**, sets overlap **0.494**.

So the loop converges to a *basin* that is about half seed-determined. On ~45% of targets
the seed is forgotten; on the rest the answer is still the seed's neighbourhood. `ctl_med`
is the sharper case: its round-1 overlap is already 0.564 and it never leaves - **the
pool medoid is its own fixed point**, which is exactly the self-confirmation the brief
warned about, isolated to the arm designed to expose it.

## L-5. Distortion inheritance is real, strong, and not the mechanism

S9-1: averaging contracts the chain and the projection buys the lost contour length back as
curvature, leaving the synthesis at a **93.47 degree** CA-CA-CA pseudo-angle against this
instrument's natives' **103.99**. Does an over-curved query retrieve over-curved windows?

**Yes, decisively.** Spearman rho between the query's pseudo-angle and the retrieved
candidates' mean pseudo-angle, across targets: **+0.907** (replace), +0.942 (direct),
+0.870 (ctl_med), +0.841 (ctl_rand), +0.725 (fuse). And in level: `ctl_med`, whose query
sits at 93.74 degrees, retrieves candidates at **94.62** - six degrees below `replace`'s
100.77.

**And it is not what is costing the accuracy.** Three independent lines:

1. **The aggregate does not drift.** `replace`'s retrieved candidates sit at 100.75 degrees
   at round 0 and **100.77** at round 5. The inheritance is cross-target - which targets get
   curved candidates - not a systematic contraction of the whole instrument.
2. **The loop IMPROVES the query's own geometry while worsening its accuracy.** The emitted
   structure's pseudo-angle goes **93.47 → 98.41**, closing 4.9 of the 10.5 degrees it was
   short, because the retrieved set is tighter (spread 2.486 → 1.551) so averaging contracts
   less. RMSD got worse over the same rounds.
3. **Correcting the distortion at the query is worse, not better.** `synth_raw` - the raw
   average rescaled to a legal bond length, which is the same conformation without the
   contraction-to-curvature conversion - retrieves candidates at **103.55 degrees, the
   closest to the natives' 103.99 of any rung in the study** - and returns **3.411 A,
   +0.211 [+0.097,+0.324] worse than the projected query**. It is a worse query in RMSD
   (4.091 vs 3.201) and RMSD is what the transfer function reads.

**Query CA-RMSD is the only coordinate that matters. This is the third independent
confirmation of S9-1's "the angle is a symptom, not a lever", now from the retrieval side.**

## L-6. Geometry validity

Every emitted structure, every arm, every round (natives on this instrument: CA-CA
**3.8122 A**, pseudo-angle **103.99**, rg **6.601**):

| arm/round | CA-CA | SD | pseudo-angle | pos-phi non-Gly (constrained) | clash | rg |
|---|---|---|---|---|---|---|
| all arms / 0 | 3.804 | 0.000 | 93.47 | 0.1086 | 0.0004 | 6.476 |
| replace / 5 | 3.804 | 0.000 | 98.41 | 0.0908 | 0.0001 | 6.576 |
| fuse / 5 | 3.804 | 0.000 | 96.86 | 0.1006 | 0.0003 | 6.559 |
| direct / 5 | 3.804 | 0.000 | 98.99 | **0.0572** | 0.0000 | 6.582 |
| ctl_med / 4 | 3.804 | 0.000 | 94.02 | 0.0410 | 0.0001 | 6.166 |
| ctl_rand / 4 | 3.804 | 0.000 | 96.96 | 0.0737 | 0.0002 | 6.523 |

Bonds are the builder's exact 3.80395 with SD 0.000 - every answer is a legal peptide.
Non-local clashes are 0.00-0.04%. **The positive-phi rate excludes `phi[0]`**, which S9-2
finding 3 showed the builder never reads, so it is the honest constrained non-glycine
number; it is 10.9% at round 0 and 5.7-10.1% after, against real library windows' **5.40%**.
These are the UNTREATED projections - S9-2's `ramah` hinge is not applied here, and it takes
that rate to ~5% for -0.002 A, so it stacks with anything below at no cost.

## L-7. Leakage and fold discipline

- **`stage_leak`**: the whole loop re-run with `nat_ca` and `rr` NaN-poisoned. **720
  deployable quantities** (6 targets x 3 arms x 3 rounds x geometry, curvature, spread,
  overlap, positive-phi), **worst |diff| 0.000e+00**. Four further in-suite tests poison the
  retrieval key, the synthesis operator, the medoid seed and the full trajectory
  independently.
- **`stage_audit`**: the library is `s7.audit.build_pool_members` - out-of-fold peptides
  plus `distogram._fold_fragments` under the standard identity filter, via `peptide_db`'s
  5 folds. **835,247 target-member pairs, 0 identity violations, 0 members from a target's
  own fold, 0 members equal to a target, worst identity anywhere 0.588** against the 0.6
  threshold.
- The check only a **structural** retriever needs, added here: the closest window in the
  whole universe to each target's own native is **1.313 A on average and 0.160 A at worst**,
  reproducing S8-4's whole-library ceiling exactly. The library does not contain the answer.
- Natives enter only as reporting labels and in the rows marked ORACLE in L-3.
- **No dev-24 pass was spent** - the best deployable arm is +0.070 A the wrong way against a
  0.16 A instrument SE, and spending a thrice-used held-out set on a negative is how a null
  becomes a claim. **The 60-target benchmark was not read.**

## L-8. Protocol table and concentration

| arm | mean | median | SD | SE | best | worst | <2.0 | <1.5 | d vs 3.454 | W/L | d vs 3.200 | W/L |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 3.454 | 3.478 | | 0.147 | | | | | - | - | - | - |
| incumbent (rd 0) | 3.201 | 2.992 | | 0.154 | | | 0.286 | 0.198 | -0.253 [-0.386,-0.121] | 82/44 | - | - |
| replace / 5 | 3.323 | 3.163 | 1.771 | 0.158 | 0.195 | 7.634 | 0.270 | 0.190 | -0.131 [-0.269,+0.007] | 71/55 | +0.123 [+0.051,+0.194] | 47/79 |
| fuse / 5 | 3.270 | 3.131 | 1.774 | 0.158 | 0.199 | 7.889 | 0.294 | 0.190 | -0.184 [-0.320,-0.048] | 73/53 | +0.070 [+0.001,+0.139] | 53/73 |
| direct / 5 | 3.304 | 3.099 | 1.837 | 0.164 | 0.189 | 9.372 | 0.286 | 0.190 | -0.150 [-0.307,+0.006] | 71/55 | +0.103 [+0.013,+0.193] | 50/76 |
| ctl_med / 4 | 3.597 | 3.455 | 2.030 | 0.181 | 0.195 | 8.828 | 0.278 | 0.190 | +0.143 [-0.110,+0.396] | 58/68 | +0.397 [+0.178,+0.615] | 38/88 |
| ctl_rand / 4 | 3.580 | 3.377 | 1.919 | 0.171 | 0.195 | 8.113 | 0.238 | 0.167 | +0.126 [-0.097,+0.349] | 66/60 | +0.380 [+0.186,+0.574] | 38/88 |

Concentration of the paired delta against the 3.454 baseline, reported the way S9-1 required:

| arm | mean | median | frac improved | top-10 share | drop-20 mean |
|---|---|---|---|---|---|
| incumbent | -0.253 | -0.116 | 0.651 | 0.60 | **-0.017** |
| replace / 5 | -0.131 | -0.029 | 0.563 | **1.18** | +0.125 |
| fuse / 5 | -0.184 | -0.069 | 0.579 | 0.86 | +0.062 |
| direct / 5 | -0.150 | -0.048 | 0.563 | 1.11 | +0.135 |

The incumbent's -0.253 already falls to -0.017 after dropping twenty targets (S9-1's
finding, reproduced). After looping, the residual is *more* concentrated still: `replace`'s
ten largest gains carry **118%** of its delta and the remaining 106 targets are net
**+0.125** worse than the baseline. Whatever survives the loop is ten targets, not a mean.

## Verdict

**Closing the loop does not move recall, and the reason is arithmetic rather than
mechanism.** The mechanism works: a native query moves universe-relative recall from 0.038
to 0.160 and returns **1.925 A**, the first sub-2 A number this pipeline has produced on the
tuning instrument. But the transfer function is `selected = 0.492 * query + 1.634`, so the
loop pays only from a query better than **3.183 A**, and the pipeline's own query is
**3.2005 A**. The deployable loop is 0.018 A short of break-even and therefore costs
+0.046 A in one round and +0.123 A at its fixed point, monotonically, on a set it genuinely
replaced (round-1 overlap 0.095, answer changed on 126/126).

The two bounds the brief set were both met and both proved not to be the binding
constraint. The set DID change, so this is not S8-6's "the selector ignored the new
candidates" null. The distortion IS inherited (rho +0.91) but correcting it makes things
worse, because query RMSD is the only coordinate the transfer function reads.

**What this makes actionable.** The loop is a *multiplier on predictor quality with a hard
threshold*, and the threshold is now measured rather than guessed. It is the only mechanism
in three sprints that has moved in-band recall by more than 0.03 - it moves it by 0.122 -
and it turns every 1 A of future query improvement into 0.49 A of answer. Any intervention
that gets a single structure to 2.5 A unlocks a further ~0.35 A for free through one round
of this loop, and to 1.7 A unlocks 0.99 A. **Re-run `s9.loop qcurve` the moment anything
reaches 3.0 A; do not deploy the loop before then.**


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / refine

> *Working paper, merged verbatim from `s9/refine_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-refine. Refinement from a good hypothesis fails too - but not where anyone was looking

`s9/refine.py` (15 resumable stages), `s9/test_refine.py` (30 checks, 0 failures),
`s9/refine_{build,oracle,train,grad,exact,visited,search,budget,null,ablate,amber,valid,leak}.json`.

Instrument: the 126-target tuning set, K=500 BLOSUM pools, read through `s9.bias.load` on
`s8.consensus2.load_cache`. `stage_build` asserts the three numbers this study stands on
before anything else runs - shipped score **3.4540**, pool best **1.7108**, S8-11's synthesis
`sc|75:fit` **3.2041**. **The 60-target benchmark was not touched by any stage in this file,
and no dev-24 pass was spent** (`DEV_ARM` is `None`; `stage_dev` raises before it will run).

---

## The one-line answer

**Refinement from a good hypothesis fails, and it fails for a reason that is now located
exactly: the search VISITS structures at 2.46 A and the objective HANDS BACK 3.24 A.** There
is no search problem, no space problem, no geometry problem and no starting-point problem.
There is **0.78 A of pure selection loss** inside the neighbourhood of a single good structure
- the same quantity S8-5 and S8-8 measured on pools, reappearing unchanged where the pool has
been removed entirely.

---

## 0. The search space, and why it is the right shape

Per residue, four local conformations: option 0 is the **incumbent's own fitted torsion**;
options 1-3 are the three most populated Ramachandran cells of the filtered top-75 candidates
at that residue, as circular means, interpolated a fraction `eta` of the way from the
incumbent. Two qubits per residue.

| | |
|---|---|
| configurations `4^n` | **2.62e5 … 4.29e9**, median 6.71e7 |
| register | **18 … 32 qubits** |
| the all-zero bitstring | rebuilds the incumbent to **9.5e-07 A** (the float32 cache's rounding) |
| CA-CA step of *every reachable structure* | **3.804 A**, worst deviation 4.5e-05 A |

**It is not enumerable.** S8-9 measured CVaR-VQE as correct-but-*unnecessary* over 128
hypotheses, for the stated reason that 128 is enumerable. That experiment is not repeated.

**Geometry is legal by construction, not by repair.** The most repeated failure mode in this
project is a gain that evaporates when validity is restored (S8-11's coordinate average scores
3.048 A on a 2.961 A bond; S9-bias gained 0.05-0.12 A on 2.65-3.13 A chains and paid +0.30 to
+0.95 A to make them peptides). The space here is *torsions* built through
`build_backbone_batch`, so there is no pre-projection object to be tempted by. Every number
below is already a peptide.

**Nothing in the objective refers to the CA-CA-CA pseudo-angle.** S9-bias measured driving it
toward 104 deg as monotonically harmful (+0.019 / +0.080 / +0.258 / +0.334 A). It is reported
for every arm and targeted by none.

## 1. The ORACLE ceiling: the neighbourhood of one good structure reaches 1.97 A

`stage_oracle` searches the same space with the same budget on the **true CA-RMSD** - a
labelled oracle diagnostic and nothing else.

| | reachable | d vs the 3.2041 incumbent | 95% CI | W/L |
|---|---|---|---|---|
| **incumbent (S8-11 synthesis)** | **3.2041** | - | - | - |
| ORACLE, eta = 0.10 | 2.8869 | **−0.3172** | [−0.3736, −0.2608] | **126/0** |
| ORACLE, eta = 0.25 | 2.5808 | **−0.6233** | [−0.7243, −0.5222] | **126/0** |
| ORACLE, eta = 0.50 | 2.2865 | **−0.9176** | [−1.0522, −0.7829] | **126/0** |
| **ORACLE, eta = 1.00** | **1.9714** | **−1.2326** | [−1.3940, −1.0712] | **126/0** |
| *for scale: the K=500 pool's own best member* | *1.7108* | | | |

At eta = 1 the neighbourhood contains, **on 126 of 126 targets**, a legal peptide better than
the incumbent, mean **1.971 A** - below this sprint's 2.0 A objective - with **55.6% of targets
below 2.0 A and 37.3% below 1.5 A**. These are lower bounds: the oracle is a
3072-evaluation annealer, not an exact solver.

## 2. The arm table - 28 arms, 126 targets, and not one of them beats doing nothing

Every arm returns a structure the deployable objective chose; `moved` is CA-RMSD from the
incumbent. Full table in `refine_search.json`; the spine of it:

| arm | mean | median | SD | d vs 3.2041 | 95% CI | W/L | d vs 3.4540 | <2 A | <1.5 A | moved |
|---|---|---|---|---|---|---|---|---|---|---|
| SHIPPED distogram argmin | 3.4540 | 3.478 | 1.646 | +0.250 | | 46/80 | - | 0.214 | 0.159 | |
| **SYNTHESIS incumbent (S8-11)** | **3.2041** | 2.966 | | - | - | - | **−0.2499** | 0.278 | 0.190 | 0 |
| **`sa@0.1`** - the best arm of 28 | **3.1968** | 2.993 | 1.715 | **−0.0073** | [−0.0273, +0.0126] | 71/55 | −0.2572 | 0.270 | 0.198 | 0.277 |
| `vqe@0.5~tr8.0` (trust region so tight it cannot move) | 3.2030 | 2.962 | 1.726 | −0.0011 | [−0.0039, +0.0017] | 61/61 | −0.2510 | 0.278 | 0.190 | **0.010** |
| `vqe@0.1` | 3.2047 | 2.998 | 1.717 | +0.0006 | [−0.0205, +0.0217] | 65/61 | −0.2493 | 0.278 | 0.190 | 0.276 |
| `sa@0.5~consonly` | 3.2113 | 2.951 | 1.739 | +0.0072 | [−0.0105, +0.0250] | 54/69 | −0.2427 | 0.278 | 0.183 | 0.139 |
| `sa@0.25` | 3.2132 | 2.992 | 1.718 | +0.0091 | [−0.0333, +0.0515] | 61/65 | −0.2408 | 0.278 | 0.190 | 0.557 |
| `vqe@0.25` | 3.2206 | 2.970 | 1.721 | +0.0166 | [−0.0281, +0.0612] | 58/68 | −0.2334 | 0.286 | 0.190 | 0.589 |
| `sa@0.5~nolegacy` | 3.2278 | 3.064 | 1.738 | +0.0238 | [−0.0048, +0.0523] | 57/69 | −0.2262 | 0.270 | 0.183 | 0.540 |
| **`sa@0.5`** (the reference arm) | **3.2390** | 3.005 | 1.713 | **+0.0349** | [−0.0185, +0.0883] | 56/70 | −0.2150 | 0.278 | 0.190 | 0.839 |
| `greedy@0.5` | 3.2709 | 3.073 | 1.721 | +0.0668 | [−0.0034, +0.1371] | 56/70 | −0.1831 | 0.286 | 0.190 | 0.895 |
| **`vqe@0.5`** | **3.3027** | 3.153 | 1.712 | **+0.0986** | **[+0.0317, +0.1656]** | 54/72 | −0.1513 | 0.262 | 0.190 | 0.924 |
| **`rand@0.5`** (uniform random) | **3.3029** | 3.150 | 1.738 | +0.0988 | [+0.0255, +0.1722] | 57/69 | −0.1511 | 0.270 | 0.183 | 0.897 |
| `sa@1.0` | 3.3178 | 3.115 | 1.738 | +0.1138 | [+0.0482, +0.1793] | 47/79 | −0.1362 | 0.254 | 0.175 | 1.065 |
| `vqe@1.0` | 3.3370 | 3.131 | 1.751 | +0.1329 | [+0.0598, +0.2060] | 42/82 | −0.1170 | 0.262 | 0.183 | 1.200 |
| `sa@0.5~nocons` | 3.3474 | 3.241 | 1.706 | +0.1433 | [+0.0405, +0.2461] | 56/70 | −0.1066 | 0.254 | 0.175 | 1.290 |
| `sa@0.5~legacyonly` - **Legacy alone** | **3.6087** | 3.348 | 1.858 | **+0.4047** | **[+0.2351, +0.5742]** | 38/88 | **+0.1547** | 0.222 | 0.135 | 1.910 |
| `sa@0.5~ramaonly` - Ramachandran alone | 3.6435 | 3.405 | 2.033 | **+0.4395** | [+0.2537, +0.6253] | 50/76 | +0.1895 | 0.270 | 0.183 | 1.680 |
| *ORACLE, eta=0.5* | *2.2865* | | | *−0.9176* | | *126/0* | | | | |

**The best of 28 arms is `sa@0.1` at 3.1968 A: −0.0073 [−0.0273, +0.0126], 71 wins to 55.**
Null. Second best is the arm whose trust region is so strong it moves 0.010 A. Nothing
material was found.

### The eta sweep is monotone, and exactly anti-correlated with the oracle

| eta | ORACLE reachable | `vqe` | `sa` |
|---|---|---|---|
| 0.10 | 2.8869 | +0.0006 | **−0.0073** |
| 0.25 | 2.5808 | +0.0166 | +0.0091 |
| 0.50 | 2.2865 | **+0.0986** | +0.0349 |
| 1.00 | **1.9714** | **+0.1329** | **+0.1138** |

**The richer the neighbourhood, the worse the deployable objective does in it.** Every step
that adds reachable accuracy adds deliverable error.

### The damage is distance

Across all 28 arms, `spearman(distance travelled, CA-RMSD damage) = **+0.9726**` (pearson
+0.8820), and

    damage = +0.2105 * (A moved) - 0.0842

Every Angstrom the objective moves the structure costs **0.21 A** of accuracy. The arm table
is, to a very good approximation, a ranking of how far each arm walked.

## 3. Where the failure is: the search walks past the answer

`stage_visited` replays the deployable search with the same seed - the trajectory is
bit-identical - while a passive recorder computes, as a labelled ORACLE column that never
touches the objective, the true CA-RMSD of all 3072 structures each run evaluates.

| | CA-RMSD | d vs incumbent | W/L |
|---|---|---|---|
| incumbent | 3.2041 | - | - |
| ORACLE over the **whole** move space | 2.2865 | −0.9176 | 126/0 |
| **ORACLE over the structures actually VISITED** | **2.4601** | **−0.7440** | **126/0** |
| what the deployable objective **SELECTED** | 3.2390 | +0.0349 | 56/70 |
| **selection loss** | | **+0.7789** | |
| targets where the search walked past a sub-2.0 A structure | | **39.7%** | |

**The search visits structures within 0.17 A of the space's own oracle, on 126 of 126 targets,
and hands back one that is worse than where it started.** This is the study's central number.
It also settles the mandate's question about which component to blame: none of the search
components, at any budget or any level of sophistication, because the search is already
walking through the answer.

## 4. Two independent proofs that the objective, not the search, is what fails

### 4a. An EXACT optimiser is worse than not moving

Restricting to 7 free residues makes the subspace 16384 configurations - enumerable
exhaustively. 40 targets, eta = 0.5:

| | CA-RMSD | d vs incumbent | 95% CI | W/L |
|---|---|---|---|---|
| incumbent | 2.9856 | - | - | - |
| **EXACT argmin of the deployable objective** | **2.9993** | **+0.0137** | [−0.0726, +0.1001] | 19/20 |
| simulated annealing, 3072 evaluations | 2.9767 | −0.0088 | [−0.0735, +0.0558] | 22/17 |
| CVaR-VQE, 3072 evaluations | 2.9713 | −0.0142 | [−0.0789, +0.0504] | 18/20 |
| *ORACLE argmin of true CA-RMSD, same subspace* | *2.4881* | *−0.4974* | [−0.6216, −0.3732] | *40/0* |

| | |
|---|---|
| searches reaching the exact optimum | sa **30%**, vqe **50%** |
| rho(objective, true CA-RMSD) over the whole subspace | **+0.4135** |
| fraction of the subspace beating the incumbent - **in the objective** | **0.061** |
| the same - **in true CA-RMSD** | **0.326** |

Three readings, and none of them is "the search is too weak".

1. **A perfect optimiser of the objective is worse than not moving.** The two approximate
   searches beat it, and they beat it *because* they fail to converge.
2. **VQE is again correct and unnecessary.** At 16384 enumerable configurations it ties an
   annealer and both tie the exact answer. S8-9 said this at 128; it holds at 16384.
3. **The mechanism in one line.** A third of the subspace is better than the incumbent in
   CA-RMSD while the objective ranks the incumbent inside its own top 6% - the objective
   already believes the incumbent is nearly optimal, and the 6% it prefers is not the 33% that
   is actually better. rho = +0.41 is real skill and is nowhere near enough.

### 4b. Optimising harder does not help, and the FIRST step already loses

`stage_budget`, eta = 0.5, simulated annealing, 126 targets:

| evaluations | objective E | CA-RMSD | d vs incumbent | 95% CI | moved |
|---|---|---|---|---|---|
| 0 (incumbent) | −3.650 | **3.2041** | - | - | 0.000 |
| 24 | −3.936 | 3.2377 | +0.0337 | [−0.0011, +0.0684] | 0.351 |
| 96 | −4.165 | 3.2369 | +0.0329 | [−0.0077, +0.0734] | 0.528 |
| 384 | −4.332 | 3.2681 | **+0.0641** | [+0.0217, +0.1064] | 0.644 |
| 1536 | −5.001 | 3.2550 | +0.0510 | [−0.0174, +0.1194] | 0.802 |
| 6144 | −5.491 | 3.2341 | +0.0300 | [−0.0299, +0.0899] | 0.915 |

The objective improves by 50% and the structure is **worse than the incumbent at every
budget, including the first 24 evaluations**. Sprint 7's 1KVG result (objective 1.564 → 1.495
while CA-RMSD 2.276 → 3.107) said that optimising a misspecified objective *harder* makes
structures worse. Measured inside a good hypothesis' own neighbourhood the statement is
sharper: **the first move the objective wants is already wrong**, and extra optimisation is
not what does the damage.

## 5. What the objective DOES know - the isotropic-move null

An arm that moves `r` A in a direction uncorrelated with the native adds `sqrt(R^2+r^2) - R` to
a CA-RMSD of `R`, because the displacement is orthogonal on average. That null is computable
per target from the arm's own step length, and it separates "the direction carries nothing"
from "the direction carries something and cannot pay for the travel".

| arm | moved | damage | isotropic null | steering credit | 95% CI | % of null recovered |
|---|---|---|---|---|---|---|
| `sa@0.1` | 0.277 | **−0.0073** | +0.0187 | **−0.0260** | [−0.0472, −0.0048] | **139%** |
| `sa@0.25` | 0.557 | +0.0091 | +0.0782 | −0.0691 | [−0.1123, −0.0258] | 88% |
| `sa@0.5` | 0.839 | +0.0349 | +0.1682 | −0.1333 | [−0.1972, −0.0694] | 79% |
| `sa@1.0` | 1.065 | +0.1138 | +0.2544 | −0.1407 | [−0.2039, −0.0774] | 55% |
| `vqe@0.5` | 0.924 | +0.0986 | +0.2092 | −0.1105 | [−0.1727, −0.0484] | 53% |
| `rand@0.5` | 0.897 | +0.0988 | +0.2068 | −0.1080 | [−0.1972, −0.0187] | 52% |
| `sa@0.5~legacyonly` | 1.910 | +0.4047 | +0.6692 | −0.2646 | [−0.4192, −0.1100] | 40% |

**Every arm beats a random move of its own length, with every interval excluding zero.** The
objective is not uninformative: it steers, and at the shortest step it steers well enough to
recover *more* than the whole random-walk penalty (139%, which is why `sa@0.1` is the one arm
with a negative point estimate). It is simply never informative enough to beat **not moving**,
and its steering decays with distance exactly as the arm table's ordering says.

That is the honest form of this sprint's negative. Not "physics is useless here" - a
50-to-139% recovery against an isotropic null is real information - but "the information is
about a fifth of an Angstrom short of what a step costs, at every step size we can take."

## 6. The four-component ablation, measured, including where it is zero

Each pair differs in exactly one thing, on identical targets. **Negative credit = the
component makes CA-RMSD worse, i.e. it does not earn its place.**

| component | role it was given | with | without | credit | 95% CI | helps/hurts |
|---|---|---|---|---|---|---|
| **VQE** | the search, vs simulated annealing at matched budget | `vqe@0.5` | `sa@0.5` | **−0.0637** | **[−0.1139, −0.0135]** | 54/70 |
| **VQE** | vs coordinate descent at matched budget | `vqe@0.5` | `greedy@0.5` | −0.0318 | [−0.0761, +0.0125] | 60/61 |
| **VQE** | **vs UNIFORM RANDOM at matched budget** | `vqe@0.5` | `rand@0.5` | **+0.0002** | [−0.0516, +0.0520] | 64/59 |
| **CVaR** | the tail level at matched T = 0.3 | `alpha=0.25` | `alpha=1` | −0.0331 | [−0.0786, +0.0124] | 57/60 |
| **CVaR** | the tail level at T = 0 (S8-9's collapse test) | `alpha=0.25` | `alpha=1` | +0.0066 | [−0.0455, +0.0588] | 64/52 |
| **Legacy** | the 15-column learned combiner as a channel | `sa@0.5` | `~nolegacy` | −0.0112 | [−0.0669, +0.0445] | 53/71 |
| **Legacy** | learned combiner vs the shipped total `lg_all` | `sa@0.5` | `~lgall` | +0.0205 | [−0.0399, +0.0809] | 59/63 |
| distogram | the shipped score as a channel | `sa@0.5` | `~nodist` | **+0.0728** | **[+0.0194, +0.1262]** | 73/45 |
| consensus | the consensus criterion as a channel | `sa@0.5` | `~nocons` | **+0.1084** | **[+0.0291, +0.1876]** | 67/46 |
| **AMBER** | restrained ff14SB/GBn2 relaxation of the output | +amber | no amber | **−0.0264** | **[−0.0377, −0.0151]** | 15/45 |

**VQE.** Genuinely present and genuinely unnecessary - this time with a sharper number than
S8-9's. Against a matched-budget annealer it is **significantly worse** (−0.0637, CI excluding
zero); against **uniform random sampling** it is **+0.0002 [−0.0516, +0.0520]**, i.e.
statistically indistinguishable. At 4^n configurations the space is finally too large to
enumerate, which was the escape S8-9 left open, and the answer did not change: what a better
optimiser buys, a misspecified objective spends.

**CVaR.** Null in both directions (−0.033 and +0.007, both CIs spanning zero), and the entropy
column says why. S8-9's mechanism was that the mean objective at low temperature collapses the
state to 0.08 of 7 bits and reproduces the shipped selector exactly, which the tail level
prevents (worth +0.113 A). **That mechanism does not reproduce here.** At 2n = 18-32 qubits
with 2n parameters the ansatz *cannot* concentrate on a single configuration:

| arm | alpha | T | state entropy at convergence (bits of 2n ≈ 26) |
|---|---|---|---|
| `vqe@0.5` | 0.25 | 0.3 | 23.18 |
| `vqe@0.5~mean` | 1.0 | 0.3 | 21.14 |
| `vqe@0.5~cvarT0` | 0.25 | 0.0 | **8.57** |
| `vqe@0.5~meanT0` | 1.0 | 0.0 | **8.66** |
| `vqe@0.5~tr8.0` | 0.25 | 0.3 | 18.12 |

At T = 0 both collapse to ~8.6 bits and **the tail level neither prevents the collapse nor
changes the result** (8.57 vs 8.66 bits; +0.0066 A). CVaR's measured role was specific to a
128-state register where the argmin is one reachable basis state; there is nothing here for it
to prevent.

**Legacy.** Its honest role was the one regime it has ever helped in - a learned 15-column
combiner on an already-near-native set (S8-8: 2.228 A against the shipped score's 2.406) - and
it was given exactly that role, fitted leave-fold-out *inside the refinement neighbourhood*.
The learning is real: the LFO combiner reaches rank skill **+0.1384** against `lg_all`'s
+0.0915 in this space, and using it instead of `lg_all` is worth +0.0205 A. **It still does not
earn its place** (−0.0112 as a channel), and as a *sole* objective it is the worst arm
measured: **3.6087 A, +0.4047 [+0.2351, +0.5742], 38/88** - it destroys the whole S8-11 gain
and then some, ending +0.155 A *worse than the shipped baseline*. Legacy has now been tested in
six roles across three sprints and has never once earned a place.

**AMBER.** Exactly the validity stage S8-12 measured, reproduced in the synthesis regime.
60 targets, one restrained ff14SB/GBn2 minimisation (k = 10) per structure, one OpenMM context
torn down between structures:

| | before | after | d | 95% CI | CA moved | median energy removed | bond after |
|---|---|---|---|---|---|---|---|
| incumbent | 3.2335 | 3.2593 | **+0.0259** | [+0.0135, +0.0383] | 0.300 A | 7.6e5 kcal/mol | 3.858 A |
| the arm | 3.2428 | 3.2692 | +0.0264 | [+0.0151, +0.0377] | 0.285 A | 1.6e5 kcal/mol | 3.857 A |

Same sign and order as S8-12's +0.011-0.012 A. It removes builder strain and costs about
0.026 A. It is a validity stage, it was used as one, and it should stay one.

**And note what the two "positive" credits actually are.** The distogram (+0.073) and the
consensus criterion (+0.108) earn credit only in the sense of *limiting the damage*: dropping
them makes the arm travel further (0.839 → 1.030 and → 1.290 A) and the extra travel is what
costs the accuracy. The consensus criterion is, on this manifold, approximately a trust region
toward the incumbent - S9-bias measured the manifold Frechet mean of the same candidate set as
being the incumbent to within 0.006 A. **Every component of the objective that earns credit
earns it by keeping the structure closer to where it started.**

## 7. Geometry and validity for every arm - and a defect nobody had reported

Bond length, CA-CA-CA pseudo-angle, non-glycine positive-phi rate, clash rate (fraction of
`|i-j| >= 3` CA pairs under 4.0 A) and Legacy steric energy. All 126 targets.

| set | CA-CA bond | CA-CA-CA angle | **posphi (non-Gly)** | clash | steric |
|---|---|---|---|---|---|
| **NATIVE** | 3.812 | **103.99** | **0.0558** | 0.0033 | - |
| **incumbent (3.204)** | **3.804** | **93.21** | **0.1973** | 0.0004 | 0.1 |
| `sa@0.1` (best arm) | 3.804 | 93.41 | 0.1449 | 0.0005 | 0.1 |
| `vqe@0.5~tr8.0` | 3.804 | 93.21 | 0.1895 | 0.0004 | 0.1 |
| `sa@0.5` | 3.804 | 94.43 | 0.1172 | 0.0038 | 0.9 |
| `vqe@0.5` | 3.804 | 94.64 | 0.0880 | 0.0030 | 0.4 |
| `sa@1.0` | 3.804 | 95.72 | 0.1198 | 0.0033 | 0.5 |
| `vqe@1.0` | 3.804 | **96.92** | **0.0543** | 0.0040 | 0.8 |
| `sa@0.5~legacyonly` | 3.804 | 95.00 | 0.1153 | **0.0134** | 2.2 |
| `sa@1.0+rama` | 3.804 | 96.92 | **0.0340** | 0.0030 | 0.3 |

**Bond length is exactly 3.804 A in every row, by construction.** No arm in this file gained
anything by emitting a non-peptide, and none could have.

**The defect.** The S8-11 synthesis' non-glycine **positive-phi rate is 0.1973 against the
natives' 0.0558** - 3.5x. This has not been reported before. It is the same event S9-bias
measured as a 10.78 deg pseudo-angle deficit, seen in torsion space: the projection buys back
the contour length that averaging contracted, and it pays for some of it in left-handed
backbone. Two independent measurements of one artefact of average-then-project.

**And the direction of every arm's geometry is the opposite of its accuracy.** As arms travel
further they become *more* native-like on every geometric axis - the angle climbs 93.21 →
96.92 toward the natives' 103.99, and positive-phi falls 0.1973 → 0.0543, and in
`sa@1.0+rama` all the way to 0.0340, *below* the natives' own 0.0558 - while CA-RMSD gets
monotonically worse. This is S8-8's in-band finding stated in geometry rather than in energy:
**inside the near-native band, the more protein-like a structure looks, the worse it is.**

### The Ramachandran lead, declared and killed

The positive-phi defect suggested a channel, and it was declared as a lead rather than folded
into a headline: a Ramachandran log-likelihood from the **fold's own out-of-fold corpus** read
with the **target's** amino acids - native-free, and the same quantity S8-8 measured as the
best single Legacy term on the deployable pool. Five arms:

| arm | mean | d vs 3.2041 | 95% CI | posphi |
|---|---|---|---|---|
| `sa@0.5` (no rama) | 3.2390 | +0.0349 | [−0.0185, +0.0883] | 0.1172 |
| `sa@0.5+rama` | 3.2895 | **+0.0854** | [+0.0164, +0.1545] | 0.0925 |
| `vqe@0.5+rama` | 3.2880 | +0.0839 | [+0.0212, +0.1467] | 0.0678 |
| `sa@1.0+rama` | 3.3788 | +0.1747 | [+0.0860, +0.2634] | **0.0340** |
| `sa@0.5~consrama` | 3.2781 | +0.0740 | [+0.0201, +0.1279] | 0.0869 |
| `sa@0.5~ramaonly` | **3.6435** | **+0.4395** | [+0.2537, +0.6253] | 0.0709 |

**It fixes the defect and makes the answer worse, monotonically** - the same shape as
S9-bias's pseudo-angle sweep, reached from a different quantity. Its rank skill inside the
move space is **+0.0381**, the lowest of any channel, which is why. Two independent validity
defects of the synthesis have now been identified, and correcting either one directly costs
accuracy.

## 8. Channel rank skill inside the move space

Spearman against the true CA-RMSD over a random sample of the move space (`stage_train`,
192 samples per target, 126 targets):

| channel | rank skill |
|---|---|
| **consensus criterion** | **+0.3292** |
| shipped distogram score | +0.2364 |
| **Legacy, 15-column learned combiner (LFO)** | **+0.1384** |
| Legacy, shipped fitted total `lg_all` | +0.0915 |
| Ramachandran, fold's own corpus | +0.0381 |

Per-fold held-out rho for the Legacy combiner: +0.1223 / −0.0090 / +0.1604 / +0.1657 / +0.2526
(lambda chosen by a nested leave-fold-out inside the training folds).

**Every channel carries positive rank information here.** The failure is not blindness. It is
that rank skill of +0.33 over a broad sample says nothing about the *extreme tail*, which is
the only part an argmin ever sees - the exact-enumeration row makes the same point with
rho = +0.41 and a perfect optimiser that still loses.

## 9. The CVaR gradient: the recorded defect, fixed, and priced against an exact reference

The envelope theorem gives `d CVaR_alpha/d theta = E_p[f(x) grad log p(x)]` with
`f(x) = -(q*-E(x))_+/alpha`, and a baseline leaves a score-function estimator unbiased **iff it
is CONSTANT in x**, since the correction term is `b * E_p[grad log p] = 0`.
`qansatz.cvar_gradient` subtracts the tail mean from the tail entries and leaves the non-tail
entries at zero, i.e. `b(x) = m * 1[x in tail]` - a function of x.
`s9/refine.cvar_grad(..., baseline="tail")` reproduces the shipped routine to **5.6e-17**
(`t_shipped_cvar_gradient_is_the_tail_form`), so what is priced is the shipped code.

The reference is exact: at 10 qubits all 1024 amplitudes of `p_theta` are closed form, so the
CVaR of the exact distribution and its central-finite-difference gradient are exact. 36 checks
over 12 draws x 3 alpha levels:

| estimator | cos with the exact gradient | ‖g‖ / ‖g_exact‖ |
|---|---|---|
| **exact expectation, CONSTANT baseline (the fix)** | **+1.000000** | **1.000** |
| exact expectation, shipped TAIL-ONLY baseline | **+0.655634** | 0.758 |
| sampled (4096 shots), constant baseline | +0.981077 | 1.005 |
| sampled (4096 shots), shipped tail-only baseline | +0.626416 | 0.703 |

The two **exact-expectation** rows are the point: they carry no sampling noise, so the defect
is **bias, not variance** - the shipped estimator points 49 degrees off the true gradient and
is 24% short in magnitude even with infinite shots. S8-9 measured +0.524 on a different
ansatz; the mechanism is now separated from the noise. Every VQE arm in this file uses the
corrected constant-baseline form.

The entropy term is verified too: the CNOT chain and its ring closure form a bijection on
bitstrings (`t_cnot_network_is_a_bijection`, 1024/1024 distinct images), so the state entropy
is the closed-form sum of Bernoulli entropies - matching the explicit 2^nq distribution to
**8.9e-16** with a gradient matching central finite differences at **cos 1.000000**.

## 10. Classical control at every size

| size | classical | quantum | verdict |
|---|---|---|---|
| **4^n = 2.6e5 … 4.3e9** (18-32 qubits), the deployable space | `sa` +0.0349, `greedy` +0.0668, `rand` +0.0988 | `vqe` +0.0986 | VQE **ties uniform random** (+0.0002) and **loses to annealing** (−0.0637, CI excludes zero) |
| **16384** (14 qubits), exactly enumerable | exact 2.9993, `sa` 2.9767 | `vqe` 2.9713 | all three within 0.028 A; the **exact** answer is the worst of them |
| *128 (7 qubits), S8-9* | *`topfrac` 3.287* | *`vqe` 3.282* | *tied there too* |

At every size tested, a classical method matches or beats the variational one. Reporting
anything else would be dressing up a null.

## 11. Per-target concentration

The mandate asks for the incumbent's concentration profile to be reproduced for the new arm.
The incumbent's own (S9-bias, on the same instrument): d = −0.2499, median −0.1085, fraction
improved 0.635, ten largest gains carrying **60.6%** of the total, −0.1070 after dropping ten,
−0.0127 after dropping twenty.

| arm | d | median d | frac improved | top-10 share | drop 5 | drop 10 | drop 20 |
|---|---|---|---|---|---|---|---|
| *S8-11 synthesis vs shipped (for reference)* | *−0.2499* | *−0.1085* | *0.635* | *0.606* | *−0.1598* | *−0.1070* | *−0.0127* |
| **`sa@0.1`** (best arm) | **−0.0073** | **−0.0041** | 0.563 | 2.664 | **+0.0039** | **+0.0132** | **+0.0269** |
| `vqe@0.5~tr8.0` | −0.0011 | +0.0000 | 0.492 | 1.680 | +0.0007 | +0.0008 | +0.0009 |
| `vqe@0.1` | +0.0006 | −0.0008 | 0.516 | −30.507 | +0.0125 | +0.0207 | +0.0348 |
| `sa@0.5` | +0.0349 | +0.0112 | 0.444 | −1.264 | +0.0657 | +0.0859 | +0.1170 |
| `vqe@0.5` | +0.0986 | +0.0198 | 0.429 | −0.391 | +0.1302 | +0.1490 | +0.1792 |

**No arm has a gain to concentrate.** `sa@0.1`'s −0.0073 turns positive after dropping the five
largest gains, its median is −0.004, and its top-10 "share" of 2.66 is the arithmetic of a near-zero
total, not a concentrated win. Contrast the incumbent's profile, which survives dropping ten
targets. This is exactly the shape S9-bias reported for its own null arms.

## 12. Leakage audit

**Clean.** `stage_leak` NaN-poisons the held-out target's native and recomputes the entire
deployable chain - alphabet, standardisation, objective, search - asserting bit-identical
output: **15 targets, worst |diff| 0.000e+00** on both the returned bitstring and the objective
value. Three further poisoning tests run inside the suite: `t_energy_ignores_the_native`
(0.000e+00), `t_search_ignores_the_native` (0.000e+00), `t_legacy_model_is_leave_fold_out`
(five distinct fold models, 15 columns each).

The native enters only as (a) the reported CA-RMSD, (b) a **training label on TRAINING folds**
for the Legacy combiner, and (c) rows explicitly named ORACLE - `stage_oracle`,
`stage_visited`'s `visited_best`, `stage_exact`'s `r_oracle`, and the `pca`/ceiling rows. Every
one is labelled at the point of use and none feeds a deployable decision.

## 13. Discipline

**No dev-24 pass was spent.** dev-24's paired SE on this arm family is ~0.35 A. The best arm on
the 126-target tuning instrument is −0.0073 A. Spending a pre-registered pass on an effect the
tuning instrument prices at seven thousandths of an Angstrom would produce a number that means
nothing (S8-13's and S9-bias's discipline note). `DEV_ARM` is `None`, `stage_dev` refuses to
run, and `t_no_dev_arm` asserts both. **The 60-target benchmark was not touched.**

One process fault is on the record: editing the chain shell script while a `sh` was executing
it started a second chain against the same JSON (POSIX `sh` re-reads a script by byte offset).
Both were killed and every affected output - `search`, `budget`, `exact` - was deleted and
re-run. Separately, `os.replace` raised a Windows share violation when a watcher read a
checkpoint mid-write, killing a 12-arm run; `_write` now retries and the interrupted
checkpoint was recovered intact. **No number in this file comes from a raced write.**

---

## The direct answer to the mandate's question

> *Does refinement from a good hypothesis work where refinement of a pool did not?*

**No - and the reason is not the one the hypothesis proposed.**

The hypothesis was that physics and search had been failing because they operated *before* a
good hypothesis existed. That framing predicted the bottleneck would move once a good starting
structure was available. It did not move. It is in exactly the same place, and it is now
measured with the pool removed as a confound:

- The neighbourhood of the single synthesised structure contains **1.97 A** structures on
  **126 of 126** targets, with 55.6% below 2.0 A. The material is there.
- The search **walks past 2.46 A structures** on 126 of 126 targets, and past a sub-2.0 A
  structure on 39.7% of them. The search is there.
- The geometry is legal by construction at every step; bond length is 3.804 A in every row of
  every table. The validity is there.
- The objective returns **3.2390**. There is **0.78 A of pure selection loss**, and a perfect
  optimiser of the objective is *worse than not moving at all*.

So the last role physics could hold - *refine, once someone else has found the basin* - is
closed, on the same evidence and by the same mechanism as the roles that were closed before it.
What is left is not a ceiling on refinement, a ceiling on search, or a ceiling on the force
field. **It is the ceiling this project has measured in every sprint from a different angle:
nothing available to it can tell a 2.4 A structure from a 3.2 A structure.** S8-5 measured it as
a 1.673 A incompetence constant on pools; S8-8 measured it as a score that selects worse than
its own pool's mean; S9-bias measured it as 90% of the displacement field being target-specific
and unpredicted; this study measures it as 0.78 A of selection loss inside a good hypothesis'
own neighbourhood, where there is no pool, no retrieval, and no discrimination between
candidates left to blame.

The one genuinely new quantity is **how close the objective gets**: it steers 52-139% of the
way against an isotropic null, and needs roughly 0.21 A of accuracy per Angstrom of travel to
break even. That is a target with a number on it, and it is a much smaller number than "build a
better force field". Whether anything can supply it is the question the next sprint inherits.

### What this licenses, and what it forbids

- **Do not** build another refinement objective out of these channels. Five of them were
  measured here (distogram, consensus, Legacy learned and unlearned, Ramachandran), all with
  positive rank skill, all in every combination, and the exact-enumeration control says a
  *perfect* optimiser of the best of them loses.
- **Do not** re-run VQE on a selection or refinement problem in this pipeline. It is now tied
  or beaten by classical search at 128, at 16384, and at 4^n. The formulation is right; there
  is nothing for it to buy.
- **Do not** correct the synthesis' geometric defects directly. Two are now measured - the
  10.78 deg pseudo-angle (S9-bias) and the 3.5x positive-phi rate (here) - and correcting
  either one is monotonically negative.
- **Do** keep AMBER as a validity stage on whatever the pipeline emits: ~10^5-10^6 kcal/mol of
  builder strain removed for +0.026 A.
- **Do** measure any future candidate signal against the **isotropic-move null**, not against
  zero. Every channel in this study beats zero and none beats not moving.


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / synth

> *Working paper, merged verbatim from `s9/synth_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-SYNTH. The synthesis operator is at its ceiling, and the ceiling is geometric

`s9/synth.py` (8 resumable stages), `s9/test_synth.py` (29 checks, all passing),
`s9/synth_{probe,contour,ext,hyp,curve,learn,conc,filt}.json`.

126-target tuning instrument, read-only from `s8/consensus2_cache`. The 60-target
benchmark was not touched and dev-24 was not spent. `stage_probe` asserts five numbers
before anything is believed - shipped score **3.454**, pool best **1.711**, S8-8's
`medoid75` **3.282**, the raw coordinate average **3.048**, S8-11's `fit` **3.201** - and
`stage_leak` NaN-poisons `rr` and `nat_ca` and asserts every deployable output is
bit-identical (worst |diff| 0.000e+00).

---

## The one-line answer

**Synthesis cannot be made better by a better operator. It is already the optimum of the
family, and what stops it is not the estimator, the metric, the manifold or the
projection - it is that the coordinate average is 10.35 A short of contour length, and
every way of paying that debt has now been measured.**

Twenty-three deployable operators. The best is **3.200 A** against the incumbent's
**3.201**. Nothing wins.

---

## 1. The operator comparison - 23 arms, no winner

Filter held at the incumbent `sc|75` throughout, so every difference is attributable to
the operator. Paired against both the 3.454 baseline and the 3.201 incumbent, 126 targets.

| arm | mean | se | median | <2 A | d vs 3.454 | 95% CI | d vs 3.201 | 95% CI | W/L |
|---|---|---|---|---|---|---|---|---|---|
| shipped score | 3.454 | 0.147 | 3.478 | 0.214 | - | - | +0.253 | [+0.121,+0.386] | 44/82 |
| S8-8 medoid | 3.282 | 0.162 | 3.095 | 0.278 | −0.172 | [−0.316,−0.027] | +0.082 | [+0.023,+0.140] | 52/74 |
| **`fit` (incumbent)** | **3.201** | 0.154 | 2.992 | 0.286 | **−0.253** | [−0.386,−0.121] | 0.000 | - | - |
| `mfit_rob` manifold Frechet, L1 | **3.200** | 0.154 | 2.943 | 0.286 | −0.254 | [−0.387,−0.120] | **−0.000** | [−0.011,+0.010] | 10/13 |
| `mfit_w` manifold Frechet, weighted | 3.203 | 0.154 | 2.939 | 0.278 | −0.251 | [−0.384,−0.117] | +0.003 | [−0.010,+0.016] | 8/10 |
| `cfit` contour-constrained | 3.205 | 0.152 | - | - | −0.249 | [−0.379,−0.119] | +0.004 | - | 59/67 |
| `wfit` precision-weighted projection | 3.207 | 0.155 | 2.964 | 0.286 | −0.247 | [−0.381,−0.113] | +0.006 | [−0.007,+0.020] | 67/59 |
| `mix75` fidelity/distance mixture | 3.208 | 0.154 | 2.997 | 0.286 | −0.246 | [−0.380,−0.112] | +0.008 | [−0.003,+0.018] | 63/63 |
| `mfit` manifold Frechet | 3.211 | 0.153 | 2.964 | 0.278 | −0.243 | [−0.374,−0.112] | +0.011 | [−0.007,+0.028] | 55/71 |
| `fit_cl3` cluster mixture, size-weighted | 3.211 | 0.154 | 2.942 | 0.278 | −0.243 | [−0.374,−0.111] | +0.011 | [−0.002,+0.023] | 61/65 |
| `fit_gmed` geometric median (L1) | 3.215 | 0.155 | 2.903 | 0.278 | −0.239 | [−0.376,−0.102] | +0.014 | [−0.007,+0.035] | 51/75 |
| `fit_gpa` GPA / Karcher mean | 3.215 | 0.153 | 2.971 | 0.278 | −0.239 | [−0.370,−0.108] | +0.015 | [−0.006,+0.035] | 51/75 |
| `fit_wgpa` precision-weighted GPA | 3.221 | 0.153 | 3.009 | 0.270 | −0.233 | [−0.363,−0.102] | +0.021 | [−0.001,+0.042] | 49/77 |
| `wfit_gpa` both weighted | 3.231 | 0.156 | 3.006 | 0.278 | −0.223 | [−0.356,−0.090] | +0.031 | [+0.002,+0.059] | 60/66 |
| `fit_clbig` largest cluster only | 3.247 | 0.161 | 2.929 | 0.286 | −0.207 | [−0.360,−0.054] | +0.047 | [−0.020,+0.113] | 53/73 |
| `dfit` distance-space consensus | 3.270 | 0.153 | 3.145 | 0.270 | −0.184 | [−0.316,−0.051] | +0.070 | [+0.020,+0.119] | 63/63 |
| `dfit_med` median distances | 3.294 | 0.157 | 3.113 | 0.278 | −0.160 | [−0.300,−0.020] | +0.093 | [+0.039,+0.148] | 50/76 |
| `dfit_w` precision-weighted distances | 3.297 | 0.154 | 3.075 | 0.270 | −0.157 | [−0.298,−0.017] | +0.096 | [+0.035,+0.157] | 51/75 |
| `fit_cleq3` cluster mixture, equal-weight | 3.331 | 0.142 | 3.117 | 0.222 | −0.123 | [−0.237,−0.010] | +0.130 | [+0.049,+0.212] | 49/77 |
| `dirint` direction average, integrated | 3.924 | 0.213 | - | - | +0.470 | [+0.232,+0.708] | +0.724 | - | 34/92 |

**Not one arm's CI against the incumbent excludes zero on the winning side.** The two that
tie it (`mfit_rob` at −0.000, `mfit_w` at +0.003) are the manifold-constrained Frechet
mean, which is the correctly-posed version of average-then-project - and being correctly
posed is worth nothing measurable.

### Oracle rows, which say where the remaining error is not

| arm | mean | d vs 3.201 |
|---|---|---|
| `O_dfit_nat` - the same projection given the NATIVE distance matrix | **0.996** | −2.205 |
| `O_bestof3` - oracle best of 3 cluster syntheses | 2.847 | −0.354 |
| `O_fit_bestscale` - oracle per-target rescale of the projection | 2.907 | −0.293 |
| `O_wfit_perfectw` - projection weighted by the TRUE per-residue errors | 3.227 | **+0.027** |

`O_dfit_nat` is the important one. **Handed the right distance matrix, this machinery
lands at 1.0 A** - 80.2% of targets below 2 A. The ideal-geometry manifold, the projection,
the optimiser and the chirality handling are all nowhere near binding. What is missing is
the target, not the machine.

`O_wfit_perfectw` is the decisive negative for the whole weighting idea: **per-residue
weighting loses even when the weights are the true errors.**

---

## 2. Where the 0.156 A projection cost actually goes - the contour ladder

The brief called closing the projection cost "a real target", and the coordinator's
diagnosis was that the synthesis over-curves because averaging contracts the chain and the
projection buys the lost length back as curvature. Both are now measured, and the target
is a mirage for a reason that is pure geometry.

| rung | constraint | RMSD | d vs `fit` | contour | pseudo-angle | end-to-end | rg |
|---|---|---|---|---|---|---|---|
| natives | - | - | - | **45.60** | **103.99** | **15.67** | 6.60 |
| a real candidate | - | 3.282 | +0.082 | 45.52 | 99.34 | 16.18 | 6.53 |
| `avg` | none | **3.048** | −0.152 | **35.25** | 121.62 | 15.79 | 6.21 |
| `dirint` | contour only, directions kept | 3.924 | **+0.724** | 45.52 | 121.62 | **20.13** | 8.12 |
| `cfit` | contour only, nearest to the mean | 3.205 | +0.004 | 45.52 | 90.89 | 15.97 | 6.43 |
| `fit` | contour + legal angles | 3.201 | 0.000 | 45.50 | 93.47 | 16.08 | 6.48 |
| `fitms` | as `fit`, 5 starts | 3.203 | +0.003 | 45.50 | 93.17 | 16.01 | 6.48 |

**The average is short by 10.35 A of contour - 22.7%. That length has to go somewhere, and
there are exactly two places it can go.**

- **Into curvature**, holding the positions: the pseudo-angle falls 121.6 → 90.9 (`cfit`)
  or 93.5 (`fit`), and it costs **+0.15 A**.
- **Into extension**, holding the directions: the end-to-end distance rises 15.79 → 20.13
  against a native 15.67, and it costs **+0.72 A**.

So the incumbent's 0.156 A is not slack to be recovered. **It is the cheaper of the only
two ways to pay a 10 A contour debt, and the alternative was measured at 4.8x the price.**

### Direction averaging is not a new operator

The coordinator's priority arm was to average in a representation that cannot contract -
the unit CA-CA step vectors. It is built, and it is **algebraically identical** to
renormalising the steps of the coordinate average:

> every candidate is an ideal-geometry chain, so every bond vector has the same length L,
> hence `mean_i(x_{i,j+1} − x_{i,j}) = L · mean_i u_{i,j}`, and normalisation kills the L.

`t_dir_average_identity` pins the two constructions to 1e-6 A. So direction averaging
inherits the coordinate average's directions exactly and adds one thing: it **integrates**
them, accumulating direction error along the chain instead of anchoring every residue
globally. That is the same compounding that makes S8-11's torsion-space circular mean
+0.618 A worse, and here it is +0.724 A.

**And `dirint` carries the SAME 121.62-degree pseudo-angle as the raw average.** That
settles the attribution the coordinator asked about: the over-curvature is not manufactured
by the projection. It is what restoring contour length at fixed positions costs, and any
operator that restores the length while staying near the mean's positions must pay it.

---

## 3. The consensus is the optimum of every direction it can be moved along

S8-11 records that synthesis travels 1.145 A from the pool's mode to the consensus and
gains 0.234 A. If that were monotone in travel, the mean would be an arbitrary stopping
point. Three native-free directions, each swept past its endpoint, raw constructions so
the projection cannot confound the answer:

| t | A: medoid → consensus75 | B: consensus500 → consensus75 | C: consensus75 → consensus10 |
|---|---|---|---|
| 0.0 | 3.282 | 3.326 | **3.048** |
| 0.5 | 3.104 | 3.111 | 3.056 |
| 0.75 | 3.060 | 3.057 | 3.087 |
| **1.0** | **3.048** | **3.048** | 3.137 |
| 1.25 | 3.066 | 3.092 | 3.205 |
| 1.5 | 3.111 | 3.186 | 3.288 |
| 2.0 | 3.276 | 3.499 | 3.494 |
| 3.0 | 3.830 | 4.407 | 4.024 |

`argmin(t)` = **1.0, 1.0, 0.0**. The consensus sits exactly at the optimum of the mode
direction and the filter direction, and moving toward the score's own top-10 is worse from
the first step. There is no direction left to travel in.

---

## 4. Multi-hypothesis: 0.51 A of headroom and no rule reaches it, at any separation

Split the filtered set into k average-linkage basins, synthesise and project each, choose
by a native-free rule. Every hypothesis is a shippable answer.

| k | ORACLE best | random hypothesis | ORACLE worst | best deployable rule | d vs not splitting |
|---|---|---|---|---|---|
| 2 | 2.950 | 3.432 | 3.914 | `cons75` 3.222 | +0.021 [−0.027,+0.069] |
| 3 | **2.847** | 3.532 | 4.157 | `cons75` 3.237 | +0.036 [−0.027,+0.099] |
| 4 | 2.794 | 3.561 | 4.279 | `size` 3.247 | +0.046 [−0.027,+0.119] |
| 6 | **2.692** | 3.577 | 4.400 | `size` 3.245 | +0.044 [−0.047,+0.135] |

Six rules - cluster size, the shipped distogram score, consensus with the unfiltered
K=500 pool, consensus with the filtered 75, distance-matrix fidelity, centrality. **The
best is statistically tied with not splitting at every k and none is ever better.**

**The regime argument is refuted.** S8-8 and S8-9 closed ranking *inside* the near-native
band; I re-asked the shipped score here because cluster syntheses are 2-4 A apart, the
coarse separation at which the same score is known to work as a filter. Its Spearman rho
with true RMSD among hypotheses is **+0.056 to +0.158**. It does not discriminate at coarse
separation either. The best signal is consensus with the filtered set at rho +0.21, and it
buys nothing.

---

## 5. The fidelity/legality trade-off is monotone - there is no free legality

Two independent knobs, both of which move the geometry toward the natives.

**(a) The mixture:** `alpha` · weighted-RMSD-to-the-average + `(1−alpha)` · distance residual.

| alpha | 1.0 | 0.9 | 0.75 | 0.5 | 0.25 | 0.1 | 0.0 |
|---|---|---|---|---|---|---|---|
| RMSD | **3.201** | 3.203 | 3.208 | 3.228 | 3.259 | 3.267 | 3.270 |
| d(X, avg) | 0.791 | 0.793 | 0.810 | 0.852 | 0.911 | 0.983 | 1.085 |
| pseudo-angle | 93.5 | 93.4 | 93.7 | 94.6 | 95.6 | 96.3 | **97.0** |

**(b) The convergence budget:**

| maxiter | 300 | 80 | 40 | 20 | 10 | 5 |
|---|---|---|---|---|---|---|
| RMSD | 3.201 | 3.203 | 3.200 | 3.202 | 3.225 | 3.239 |
| pseudo-angle | 93.5 | 93.6 | 94.0 | 94.9 | 96.9 | **98.7** |

Both knobs buy pseudo-angle and both cost RMSD, monotonically, and the cost arrives before
the angle has closed a third of its 10.5-degree gap. **The pure-fidelity corner is the
optimum, and early stopping should not be used to buy validity** (+0.039 at maxiter 5) -
which independently agrees with the torsion-prior agent's finding that a prior gets there
at full convergence for free while early stopping costs 0.19 A.

The distance-space arms make the same point from the other side and are the most
informative failure in the study: `dfit_w` lands at **99.3 degrees**, exactly a real
candidate's value, and is **+0.096 A worse**. Making the structure more native-like in the
angle coordinate makes the answer worse. The angle is a symptom, not a lever.

### The projection degeneracy, measured

Five starts (four generic conformations plus the medoid's own torsions). The objective
values span **0.026 A** and the multi-start answer is **+0.003** against the single warm
start (51W/54L). The degeneracy is real, but on this objective - with no prior term -
branch selection is not where the accuracy is.

---

## 6. The shared bias is real, reproducible, 13 SE from zero, and not a lever

The residual to the native, expressed in each residue's own local frame, over 1633
residues on 126 targets:

    mean  [-0.068, +0.739, +0.016] A     SE  [0.048, 0.056, 0.050]     SD  [1.93, 2.27, 2.04]

**+0.739 A along the curvature normal at 13.2 SE**, independently reproducing the bias
agent's +0.863 ± 0.066 (13.1 SE) from a different construction. It is a real, deterministic
signature of the operator.

And it is 0.325 of one standard deviation on that axis. A leave-fold-out ridge on
deployable per-residue features (per-residue spread, radial position, chain position,
distance to the ends, per-residue distance-matrix spread, chain length, mean spread) has
**negative out-of-fold R²** (−0.12 to −0.14), and applying the learned correction and
re-projecting costs **+0.014 [+0.002,+0.027] to +0.024 [+0.008,+0.040]**.

A bias can be significant in the mean and still be swamped by its own variance. This one is.

---

## 7. Per-target concentration - the mandated breakdown

| arm | mean | gain | median gain | top-10 share | drop 10 | 95% CI | W/L | drop 20 |
|---|---|---|---|---|---|---|---|---|
| `fit` incumbent | 3.201 | −0.253 | −0.116 | **0.604** | −0.109 | [−0.214,−0.004] | 72/44 | **−0.017** |
| `mfit_rob` | 3.200 | −0.254 | −0.116 | 0.605 | −0.109 | [−0.215,−0.002] | 71/45 | −0.015 |
| `mfit_w` | 3.203 | −0.251 | −0.110 | 0.608 | −0.107 | [−0.213,−0.001] | 72/44 | −0.011 |
| `wfit` | 3.207 | −0.247 | −0.115 | 0.617 | −0.103 | [−0.211,+0.005] | 72/44 | −0.005 |
| `medoid` | 3.282 | −0.172 | −0.064 | 0.939 | −0.011 | [−0.125,+0.102] | 64/47 | +0.091 |
| **`avg` (illegal)** | 3.048 | −0.406 | −0.236 | **0.386** | −0.270 | [−0.380,−0.161] | 85/31 | **−0.173** |
| `O_bestof3` | 2.847 | −0.607 | −0.416 | 0.317 | −0.450 | [−0.563,−0.338] | 89/27 | −0.337 |

Every arm within 0.01 A of the incumbent has the incumbent's profile: ten targets carry
60% of the gain and there is essentially nothing left after twenty are dropped.

**The exception is the diagnostic one.** The raw average's gain is broadly distributed -
top-ten share 0.386, still −0.173 after dropping twenty - and the projection is what
concentrates it. That is the same fact as section 2 seen a third way: the projection's
cost falls on the median target while the average's gain does not.

---

## 8. The filter-size sweep, run last - a plateau, and a stable operator ordering

Swept only after the operator comparison was complete, so the operators were never chosen
on it. A plateau confirms; a spike would have meant the cell was doing the work.

| m | 25 | 50 | **75** | 100 | 150 | 250 |
|---|---|---|---|---|---|---|
| `avg` (RAW) | 3.075 | 3.068 | **3.048** | 3.056 | 3.072 | 3.133 |
| `fit` | 3.233 | 3.229 | **3.201** | 3.222 | 3.250 | 3.322 |
| `mfit_w` | 3.236 | 3.226 | **3.203** | 3.222 | 3.252 | 3.303 |
| `wfit` | 3.243 | 3.239 | **3.207** | 3.232 | 3.253 | 3.310 |
| `medoid` | 3.369 | 3.340 | **3.282** | 3.330 | 3.383 | 3.488 |

`fit` spans **0.049 A over m = 25 to 150** with a shallow optimum at 75 - a plateau, as
S8-11 reported for the projection family and unlike the medoid's spike. More usefully:
**the operator ordering is identical at every filter size.** The operators are not tied
only in the incumbent's cell; they are tied everywhere, which is what makes the null a
property of the operator family rather than of one choice of m.

---

## 9. What this means

**How much better can synthesis be made? By this route, essentially not at all.**

The estimator is exhausted. GPA/Karcher, precision-weighted GPA, the L1 geometric median,
the manifold-constrained Frechet mean, distance-space consensus, robust and trimmed
variants, cluster restriction and cluster mixtures, extrapolation in three directions, a
weighted projection, a mixture objective, a convergence sweep and a learned local decoder -
23 deployable arms - and the best is 0.001 A from the incumbent with a CI of ±0.011.

**What limits it, stated plainly, in order of size:**

1. **The target, not the machine.** Given the native distance matrix the same projection
   returns **0.996 A**. Every part of the construction is 2.2 A better than the answer it
   currently produces. Nothing about the operator is binding.
2. **Recognition, still.** The oracle-75 filter reaches 2.220 (S8-13) and oracle
   best-of-6 hypotheses reaches 2.692 - both are ranking problems, and ranking has now
   failed in-band (S8-8, S8-9) *and* at coarse separation (section 4).
3. **A geometric debt that cannot be forgiven.** The average is 22.7% short of contour
   length. Paying it in curvature costs 0.15 A; paying it in extension costs 0.72 A. This
   is not a defect to be engineered away - it is what averaging points about a curve does,
   and the incumbent already takes the cheaper option.

**What I would not spend another sprint on:** better means, better metrics on shape space,
better projections, per-residue weighting, torsion or distance representations, or learned
local corrections. Each is now measured with a confidence interval on 126 paired targets.

**What the numbers say is left:** the distance target. `O_dfit_nat` = 0.996 A is the
largest oracle gap in this study by a factor of six, and it is the only one that is not a
ranking problem in disguise.

---

## Discipline notes

- **dev-24 was not spent.** The best arm is −0.000 [−0.011,+0.010] against its incumbent.
  A pre-registered pass at n=24 and SE ~0.35 cannot resolve that, and spending one would
  produce a number that means nothing. This follows S8-13's precedent.
- **Three tests caught real defects while being written**, each of which would have
  produced a wrong published number: `CA_CA` was carried as 3.80 when the builder emits
  3.803954938 (so the validity assertion was checking the wrong constant); `local_frames`
  clamped its triplet at the chain ends, collapsing the second axis onto a fixed world
  vector and destroying the pose-equivariance that makes the residual regression
  meaningful; and `mfit` recorded the inner objective, which moves every round and is not
  the quantity the iteration descends.
- **Two bugs were caught by their own confidence intervals.** The leave-fold-out loop in
  `stage_learn` paired target i's corrected RMSD with target j's incumbent, which showed up
  as a ±0.45 A interval on a +0.024 A difference; corrected, it is ±0.016. And
  `stage_hyp`'s rule dispatch materialised its sort key before the branch that handles the
  `central` rule, crashing the analysis after all 126 targets had been computed - the
  per-target cache made the re-analysis free.
- **Chirality metric.** `frac_dihedral_positive` here is the CA-trace pseudo-dihedral over
  n−3 well-defined values, not backbone phi, so the coordinator's caveat about `phi[0]` and
  `psi[n-1]` never being read by the builder does not apply to this table. The natives
  score 0.318 and `fit` scores 0.218 - the same over-helical signature, measured on a
  quantity with no undefined entries.
- Every arm that is offered as an answer returns an exact ideal-geometry chain: CA-CA
  3.803954938 with SD < 1e-15, zero pseudo-angles outside 75-150 degrees, clash rate
  0.0000-0.0013 against the natives' 0.0033. The raw constructions (`avg`, `gpa`, `wgpa`,
  `gmed`) are reported and labelled RAW because they are diagnostics, not answers - their
  bonds are 2.94-3.04 A and 25% of their pseudo-angles are out of range.


<!-- -------------------------------------------------------------- -->

## Appendix: S9 / traj

> *Working paper, merged verbatim from `s9/traj_FINDINGS.md`. Not edited. Superseded by the coordinator record above wherever they differ.*

# S9-traj. The trajectory is a better set than the pool, and it is the wrong kind of better

`s9/traj.py` (8 resumable stages), `s9/test_traj.py` (26 checks, 0 failures),
`s9/traj_{build,traj,geom,decomp,decomp_seed,synth,multi,valid,leak}.json`.

Instrument: the 126-target tuning set, K=500 BLOSUM pools, read through `s9.bias.load`.
`stage_build` asserts the four numbers this study stands on before anything else runs -
shipped **3.4540**, pool best **1.7108**, S8-11's synthesis **3.2041**, and the refinement
study's own visited-set oracle **2.4601** with its argmin readout **3.2390**.
**The 60-target benchmark was not touched by any stage in this file, and no dev-24 pass was
spent** (`DEV_ARM` is `None`; `stage_dev` raises before it will run).

---

## The one-line answer

**Yes, a set generated by a misdirecting objective still contains information - and no,
consensus cannot extract it, for a reason that is arithmetic rather than empirical.** The
visited set is a *better* set than the retrieval pool it would replace (mean 3.306 A against
the filtered top-75's 3.551, −0.245 [−0.304, −0.186]), and its good structures are genuine
discoveries 1.58 A away from the starting point, not the starting point itself. But **87-95%
of its error is bias common to every member**, because every member is a perturbation of one
structure. Consensus converts *independent* error into accuracy, and the trajectory has
almost none. **51 arms, seven generators, three variants: not one beats the incumbent, and
the best beats the same-code-path control by 0.008 A.**

---

## 0. The precondition, which the mandate required to be measured first

The trajectory is generated BY the objective, so it may carry nothing the objective has not
already spent. `stage_geom` measures its geometry before anything is built from it. The `sa`
generator is bit-identical to `refine.stage_visited` - same seed, same kernel, same schedule
- and is asserted at 2.4601 / 3.2390 before the stage will run.

| | |
|---|---|
| mean pairwise CA-RMSD within the visited set | **1.2076 A** |
| RMS deviation about its own centroid | 0.9268 A |
| mean / max CA-RMSD from the incumbent | 0.9121 / 2.8169 A |
| **CA-RMSD(best visited structure, incumbent)** | **1.5773 A** |
| fraction of targets where that exceeds 0.5 A | **0.778** |
| position of the best structure in the walk | 0.501 of the way through |
| fraction of targets whose best is in the last half | 0.492 |
| running minimum at 25 / 50 / 75 / 100% of the walk | 2.554 / 2.501 / 2.470 / 2.460 |
| rho(objective, true CA-RMSD) over the visited set | +0.1409 |

**The precondition passes.** The set is not a point, and its good structures are not the
incumbent wearing a disguise: they sit 1.58 A away on 77.8% of targets. Most of the
discovery happens early - 75% of the eventual gain is inside the first quarter - but "early"
here means after 768 evaluations, not "at the start". This is the answer to the mandate's
stated failure mode: the trajectory is *not* dominated by structures near the starting
synthesis.

### And the set is BETTER than the pool it would replace

| set | best | p05 | median | mean |
|---|---|---|---|---|
| visited (3072) | 2.460 | 2.864 | 3.273 | **3.306** |
| K=500 retrieval pool | **1.711** | 2.610 | | 4.453 |
| the filtered top-75 the incumbent averages | 2.306 | | | 3.551 |

- visited mean vs top-75 mean: **−0.2447 [−0.3038, −0.1856]**
- visited best vs pool best: **+0.7493 [+0.5797, +0.9188]**
- fraction of visited structures better than the top-75's mean: **0.714**

The trajectory is a *tighter, better-centred* set with a *truncated upper tail*. That is the
opposite trade from the pool, and §1 says why it is the wrong one.

## 1. Why - the S8-11 error decomposition, run on both sets for the first time

Consensus averages away the part of a set's error that is independent between members and
leaves the part common to all of them untouched. With each member optimally superposed onto
the native (an ORACLE construction, diagnostic only):

    e_i     = X_i − native            MSE = mean_i ||e_i||^2 / n
    bias^2  = ||mean_i e_i||^2 / n    share = bias^2 / MSE

`bias` is the **ceiling of any unweighted coordinate average over the set**, before
projection - what infinitely many members would give you.

| set | m | member RMSD | **bias (the ceiling)** | independent | share |
|---|---|---|---|---|---|
| filtered top-75 (what the incumbent averages) | 75 | 3.6148 | **2.8148** | **2.0888** | 0.620 |
| visited, 512 sampled | 512 | 3.3292 | **3.1006** | 1.0598 | **0.866** |
| visited, the objective's own top-75 | 75 | 3.2723 | **3.1927** | 0.5874 | **0.952** |

**The ceiling predicts the measurement.** Observed pre-projection RMSDs are 3.052 (pool
average), 3.160 (trajectory average) and 3.216 (objective top-75 average) - each within
0.06 A of its own set's bias floor, approached from above, as finite-m sampling requires.
This decomposition costs ten minutes and would have priced the entire study before it ran.

Read across and everything below is settled:

- the pool's members are individually **worse** (3.615 vs 3.329) and its average is
  **better**, because 2.09 A of its error is independent and averageable;
- the trajectory's independent error is **halved** (1.06 A) and the objective's own tail of
  it is **quartered** (0.59 A) - the tighter the objective squeezes, the less is left to
  average;
- and the trajectory's **bias floor is worse than the pool's** (3.10 vs 2.81) and worse than
  the incumbent already achieves, so no weighting, union or reweighting can reach where the
  pool already is.

### The one intervention the decomposition endorses, built and measured

If the problem is that every visited structure is a perturbation of **one** structure, then
start the walk somewhere else. `seed` and `seedsa` run the identical kernel with the 64
chains started at 64 **different pool members** - the configuration in the same move space
nearest each candidate's own torsions, native-free - instead of all at the incumbent.

| set (`seed` generator) | m | member | **bias** | independent | share |
|---|---|---|---|---|---|
| filtered top-75 (pool) | 75 | 3.6148 | 2.8148 | 2.0888 | 0.620 |
| visited, 512 sampled | 512 | 3.6810 | **3.2180** | **1.6255** | **0.768** |
| visited, objective's top-75 | 75 | 3.3201 | 3.1136 | 0.9873 | 0.881 |

**The mechanism is confirmed and the trade is bad.** Seeding from independent starts really
does manufacture averageable error - independent 1.06 → 1.63 A, share 0.866 → 0.768 - and it
raises the bias floor by more than it buys (3.10 → 3.22). The arms follow exactly:
`seed/t_all` **+0.2425**, `seed/u_all` +0.0860, best `seed/u_top64` +0.0033.
**You can buy independent error, and it costs more bias than it is worth.**

## 2. The arm table - 51 arms, 7 generators, 3 variants, no winner

Every arm superposes its set onto the incumbent, coordinate-averages, and projects onto the
manifold of ideal-geometry chains (`s9.synth.fit_w`, incumbent-seeded; `stage_multi` repeats
the headline arms multi-start, §5). `pre` is the pre-projection RMSD and `preB` its CA-CA
bond, **reported together always**, because pre-projection gains of 0.05-0.12 A on 2.65-3.13
A chains have cost +0.30 to +0.95 A to legalise twice already in this sprint. `mv` is how far
the emitted structure moved from the incumbent.

Generators: `sa` = the deployable arm at eta 0.5 (bit-identical to the refinement study),
`damp` = the same walk with the objective's authority cut 4x, `iso` = **objective-free**
(every proposal accepted), `sa01` = eta 0.1, `rand` = uniform sampling of the whole space,
`seed`/`seedsa` = objective-free / full-objective walks started at 64 different pool members.
Subset rules: `t_` trajectory alone, `u_` trajectory ∪ retrieval pool, `x_` pool ∪ **all
five** trajectories.

| gen/arm | mean | median | <2 A | d vs 3.2041 | 95% CI | W/L | d vs 3.4540 | pre | preB | mv |
|---|---|---|---|---|---|---|---|---|---|---|
| *ORACLE over the visited set* | *2.4601* | | | *−0.7440* | | *126/0* | | | | |
| **`sa01/t_thin75`** - best of 51 | **3.2019** | 2.949 | 0.278 | **−0.0022** | [−0.0133, +0.0090] | 63/63 | −0.2521 | 3.198 | 3.790 | 0.14 |
| `sa01/t_all` | 3.2021 | 2.959 | 0.278 | −0.0020 | [−0.0130, +0.0091] | 62/64 | −0.2519 | 3.198 | 3.790 | 0.13 |
| `sa01/t_top64` | 3.2033 | 2.969 | 0.278 | −0.0007 | [−0.0171, +0.0156] | 66/60 | −0.2507 | 3.202 | 3.798 | 0.21 |
| `sa01/u_all` | 3.2039 | 2.935 | 0.278 | −0.0002 | [−0.0073, +0.0070] | 61/65 | −0.2501 | 3.093 | 3.312 | 0.09 |
| **incumbent (S8-11)** | **3.2041** | 2.966 | 0.278 | - | - | - | −0.2499 | 3.048 | 2.961 | 0 |
| `sa/u_all` - **variant 2, the union** | 3.2066 | 2.941 | 0.286 | +0.0026 | [−0.0123, +0.0174] | 60/66 | −0.2474 | 3.080 | 3.221 | 0.24 |
| `seed/u_top64` | 3.2074 | 2.958 | 0.286 | +0.0033 | [−0.0221, +0.0288] | 69/57 | −0.2467 | 3.080 | 3.219 | 0.38 |
| `sa/t_boltz1` | 3.2096 | 2.989 | 0.294 | +0.0055 | [−0.0266, +0.0376] | 65/61 | −0.2444 | 3.168 | 3.645 | 0.43 |
| **`sa/pool`** - the incumbent through THIS code path | **3.2099** | 2.979 | 0.278 | +0.0058 | [−0.0040, +0.0156] | 54/72 | −0.2441 | **3.052** | **2.964** | 0.13 |
| `sa/x_top64` - pool ∪ all five trajectories | 3.2105 | 2.986 | 0.286 | +0.0064 | [−0.0271, +0.0399] | 65/61 | −0.2435 | 3.143 | 3.529 | 0.40 |
| `sa/t_all` - **variant 1, the trajectory alone** | 3.2138 | 2.949 | 0.286 | +0.0097 | [−0.0166, +0.0360] | 58/68 | −0.2402 | 3.160 | 3.588 | 0.39 |
| `sa/t_late` (last quarter of the walk) | 3.2179 | 2.973 | 0.286 | +0.0138 | [−0.0197, +0.0473] | 62/64 | −0.2361 | 3.164 | 3.587 | 0.49 |
| `sa/t_top256` | 3.2245 | 3.047 | 0.286 | +0.0204 | [−0.0245, +0.0654] | 57/69 | −0.2295 | 3.199 | 3.703 | 0.56 |
| `seed/t_top64` | 3.2294 | 3.006 | 0.278 | +0.0253 | [−0.0195, +0.0701] | 59/67 | −0.2246 | 3.178 | 3.598 | 0.63 |
| `sa/t_top64` | 3.2306 | 3.042 | 0.278 | +0.0265 | [−0.0211, +0.0740] | 62/64 | −0.2234 | 3.208 | 3.716 | 0.62 |
| `sa/t_top16` | 3.2341 | 3.009 | 0.278 | +0.0301 | [−0.0173, +0.0774] | 58/68 | −0.2199 | 3.217 | 3.740 | 0.66 |
| `damp/t_all` - **variant 3, damped 4x** | 3.2372 | 3.009 | 0.278 | +0.0331 | [−0.0040, +0.0701] | 51/75 | −0.2168 | 3.164 | 3.488 | 0.51 |
| `sa/x_all` | 3.2522 | 3.066 | 0.278 | +0.0481 | [+0.0066, +0.0897] | 52/74 | −0.2018 | 3.154 | 3.365 | 0.47 |
| `sa/x_nopool` - five trajectories, no pool | 3.2666 | 3.102 | 0.270 | +0.0626 | [+0.0125, +0.1127] | 52/74 | −0.1874 | 3.192 | 3.475 | 0.54 |
| `seedsa/t_all` | 3.2960 | 3.035 | 0.270 | +0.0919 | [+0.0211, +0.1627] | 51/75 | −0.1580 | 3.201 | 3.381 | 0.75 |
| `iso/t_all` - **variant 3, objective-FREE** | 3.3534 | 3.131 | 0.262 | +0.1493 | [+0.0562, +0.2425] | 47/79 | −0.1006 | 3.260 | 3.365 | 0.88 |
| `seed/t_all` - seeded at 64 pool members | 3.4465 | 3.208 | 0.254 | +0.2425 | [+0.1082, +0.3767] | 51/75 | −0.0075 | 3.338 | 3.297 | 1.23 |
| `rand/t_all` - uniform sampling | 3.4491 | 3.179 | 0.262 | +0.2450 | [+0.1153, +0.3748] | 44/82 | −0.0049 | 3.341 | 3.308 | 1.18 |

Full 51-row table in `traj_synth.json`. **The best arm of 51 is `sa01/t_thin75` at 3.2019,
d = −0.0022 [−0.0133, +0.0090], 63 wins to 63 - a coin flip.** Against `sa/pool`, which is
the incumbent's own recipe read through this file's code and lands at 3.2099, the best arm is
worth **−0.008 A**. Neither is a result. 51 arms against a paired SE of 0.015-0.05 A means
the best point estimate is selected noise, and the headline is the null, not the −0.0022.

`<2 A` never leaves 0.254-0.294 and `<1.5 A` never leaves 0.167-0.198, against the
incumbent's 0.278 / 0.190. Nothing moved a target across a threshold either.

### The three variants, answered directly

1. **Consensus over the trajectory alone.** `sa/t_all` **+0.0097 [−0.0166, +0.0360]**. Null.
   Every subset rule of it - top-16/64/256/1024, deduplicated top-75, Boltzmann at two
   temperatures, the last quarter, the last tenth, a random thinning - lands between +0.0055
   and +0.0357, ordered monotonically by how far the rule lets the consensus travel.
2. **Trajectory ∪ retrieval pool.** `sa/u_all` **+0.0026 [−0.0123, +0.0174]**,
   `sa/u_thin75` +0.0030, `sa/u_top64` +0.0060. The union is *better than the trajectory
   alone and no better than the pool alone*, which is what a convex combination of a 2.81 A
   bias floor and a 3.10 A bias floor has to be. Unioning **all five** trajectories with the
   pool (`x_top64` +0.0064) does not change it; dropping the pool (`x_nopool` +0.0626) makes
   it worse. **The union hypothesis - "the two sets have different biases and the union may
   beat either" - is measured and false: the biases differ, and one of them is simply
   worse.**
3. **Objective-free or objective-damped search.** Damping 4x (`damp/t_all` +0.0331) is worse
   than not damping; switching the objective off (`iso/t_all` +0.1493 [+0.0562, +0.2425]) is
   significantly worse; uniform sampling (`rand/t_all` +0.2450) worse still; seeding
   elsewhere (`seed/t_all` +0.2425) the same. **Variant 3 fails in the direction opposite to
   the one that motivated it.** And the ORACLE column runs the exact reverse - the looser
   generators visit *better* structures (`rand` 2.348, `seed` 2.380, `iso` 2.394, `damp`
   2.408 against `sa`'s 2.460) and their consensus is worse. The refinement study's finding
   that every arm steers 52-139% better than a random direction is confirmed from the other
   side: **the objective's steering is what keeps the visited set's bias floor from
   drifting, and removing it costs more than the extra spread buys.**

## 3. The exchange rate is a property of the space, not of the readout

Across all 51 arms - seven generators, three variants, eleven subset rules:

    spearman(distance moved, CA-RMSD damage) = +0.9048   (pearson +0.8969)
    damage = +0.2182 * (A moved) − 0.0699

The refinement study measured, for the **objective-argmin** readout over 28 arms:

    damage = +0.2105 * (A moved) − 0.0842        rho +0.9726

**Two readout operators that share no code - a single-structure argmin, and a coordinate
consensus over up to 3072 structures followed by a manifold projection - pay the same
0.21 A per Angstrom.** The exchange rate is not a defect of argmin selection that a better
readout routes around. It is the price of travel in this neighbourhood given the information
available, and changing what you do with the trajectory does not change it. This is the one
genuinely new quantity this study adds, and it closes the readout question the way the
refinement study closed the search question.

It also explains the arm ordering completely. The best arms are the ones that barely move
(`sa01/u_all` 0.09 A, `sa01/t_all` 0.13 A, `sa01/t_thin75` 0.14 A) - the synthesis-readout
counterpart of `vqe@0.5~tr8.0`, the refinement study's second-best arm, whose trust region
held it within 0.010 A of the incumbent.

## 4. Per-target concentration

Against the incumbent (the mandate's requested reproduction of the incumbent's own profile):

| arm | d | median | frac improved | top-10 share | drop 5 | drop 10 | drop 20 |
|---|---|---|---|---|---|---|---|
| *S8-11 synthesis vs shipped, for reference* | *−0.2499* | *−0.1085* | *0.635* | *0.606* | *−0.1598* | *−0.1070* | *−0.0127* |
| `sa01/t_thin75` (best of 51) | −0.0022 | −0.0002 | 0.500 | 4.864 | +0.0048 | +0.0091 | +0.0167 |
| `sa01/t_all` | −0.0020 | +0.0002 | 0.492 | 5.439 | +0.0050 | +0.0095 | +0.0167 |
| `sa/u_all` | +0.0026 | +0.0006 | 0.476 | −4.933 | +0.0115 | +0.0166 | +0.0254 |
| `sa/t_all` | +0.0097 | +0.0032 | 0.460 | −2.297 | +0.0260 | +0.0347 | +0.0498 |
| `sa/pool` (control) | +0.0058 | +0.0001 | 0.429 | −1.368 | +0.0122 | +0.0149 | +0.0183 |

**No arm has a gain to concentrate.** The best arm's median is −0.0002, exactly half its
targets improve, and it turns positive after dropping the five largest gains. Its "top-10
share" of 4.86 is the arithmetic of a near-zero denominator, not a concentrated win - the
same shape the refinement study and S9-bias reported for their own null arms.

Against the **shipped** baseline every arm reproduces the incumbent's profile and adds
nothing to it:

| arm | d vs shipped | median | frac improved | top-10 share | drop 10 | drop 20 |
|---|---|---|---|---|---|---|
| *the record's S8-11 profile* | *−0.2499* | *−0.1085* | *0.635* | *0.606* | *−0.1070* | *−0.0127* |
| `sa/pool` (control) | −0.2441 | −0.1108 | 0.635 | 0.617 | −0.1016 | −0.0079 |
| `sa01/t_thin75` (best) | −0.2521 | −0.1107 | 0.643 | 0.601 | −0.1093 | −0.0154 |
| `sa/u_all` | −0.2474 | −0.1033 | 0.643 | 0.610 | −0.1048 | −0.0109 |

**The gain these arms report against the shipped baseline is S8-11's gain arriving through a
different pipe.** Ten targets carry ~61% of it and it is gone after twenty, exactly as the
record says. Nothing here adds a concentrated win of its own, or dilutes the incumbent's.

## 5. Geometry and validity, for every emitted structure

Bond length, its SD, CA-CA-CA pseudo-angle, non-glycine positive-phi rate (**residue 0
excluded** - `phi[0]` is never read by the builder, so including it averages a coin flip into
every target) and clash rate (fraction of |i−j| ≥ 3 CA pairs under 4.0 A). All 126 targets;
full 45-row table in `traj_valid.json`.

| set | rmsd | CA-CA | sd | angle | posphi | clash |
|---|---|---|---|---|---|---|
| **NATIVE CA trace, refit through the arms' own projection** | - | 3.812 | 0.015 | **103.99** | **0.2513** | 0.0033 |
| **the same, multi-start** | - | 3.812 | 0.015 | 103.99 | **0.1737** | 0.0033 |
| `sa/pool` (the incumbent) | 3.2099 | 3.804 | 0.000 | 93.08 | 0.2197 | 0.0006 |
| `sa/t_all` | 3.2138 | 3.804 | 0.000 | 92.52 | 0.2056 | 0.0004 |
| `sa/u_all` | 3.2066 | 3.804 | 0.000 | 92.74 | 0.2074 | 0.0005 |
| `sa/t_top64` | 3.2306 | 3.804 | 0.000 | 93.66 | 0.1983 | 0.0024 |
| `sa/x_nopool` | 3.2666 | 3.804 | 0.000 | 91.86 | 0.2033 | 0.0004 |
| `iso/t_all` | 3.3534 | 3.804 | 0.000 | 91.34 | 0.2020 | 0.0003 |
| `rand/t_all` | 3.4491 | 3.804 | 0.000 | 91.48 | 0.1885 | 0.0003 |
| *record references, not computed here* | | 3.812 | | *103.99* | *0.0558* (real windows *0.0540*) | *0.0033* |

**Bond length is 3.804 A with SD < 0.001 in every emitted row, by construction.** Nothing in
this file gained anything by emitting a non-peptide, and nothing could have.

**The pseudo-angle is 91.3-93.9 against the natives' 103.99 in every row**, and it moves the
*wrong* way for the worse arms - arms that travel further are *more* over-curved here, the
reverse of the refinement study's arms, which became more native-like as they got worse. The
incumbent's 10.8 degree deficit is untouched by every operator in this file, as S9-bias and
S9-synth both said it would be.

### A defect the record attributes to the wrong operator

`refine`'s validity table reported the synthesis' non-glycine positive-phi rate of 0.1973
against the natives' 0.0558 as a 3.5x defect of average-then-project. **The control was never
run.** The instrument's natives are CA-only traces, so a CA trace does not determine phi/psi;
the right comparison is the same projection operator applied to the **native trace itself**.

It comes back at **0.2513** - higher than every deployable arm in the table - and multi-start
drops it to **0.1737**. So the positive-phi rate is **the projection operator's, not the
consensus'**: inverting a CA trace onto ideal geometry takes a left-handed branch about a
fifth of the time whatever trace it is handed, and the incumbent's 0.1973 is *below* what the
operator does on the right answer. This is S9-2's degeneracy - two near-equidistant torsion
solutions, one plausible and one not - measured for the first time on the native, and it says
the "defect" is a property of CA-to-torsion inversion rather than evidence about the
synthesis. `t_projection_positive_phi_is_the_operators_not_the_consensus` pins it.

### Multi-start changes the projection and not the answer

`stage_multi` repeats five arms with `fit_multi` (four generic starts plus the incumbent's
torsions) instead of the single incumbent-seeded start:

| arm | single-start | multi-start | d | objective spread across starts |
|---|---|---|---|---|
| `pool` | 3.2099 | 3.2112 | **+0.0013** | 0.0282 |
| `t_all` | 3.2138 | 3.2138 | +0.0000 | 0.0649 |
| `t_top64` | 3.2306 | 3.2305 | −0.0000 | 0.0666 |
| `u_top64` | 3.2100 | 3.2099 | −0.0002 | 0.0387 |
| `u_all` | 3.2066 | 3.2067 | +0.0001 | 0.0312 |

The degeneracy is real and measurable - the objective spans 0.03-0.07 A across starts - and
multi-start finds a genuinely nearer point of the manifold every time (mean projection
objective 0.7786 → 0.7783 on `pool`). **It buys 0.0000 A of accuracy, and on the pool arm it
costs 0.0013.** The branch a warm start from the incumbent's own torsions lands on is as good
as the nearest one. Multi-start matters for the *validity* row (0.2513 → 0.1737 positive-phi
on the native) and not for RMSD.

## 6. Leakage audit

**Clean.** `stage_leak` NaN-poisons the held-out target's native and re-runs the **entire**
deployable chain - the 3072-evaluation trajectory, every subset rule, the superposition, the
average and the projection - asserting bit-identical output for **all 21 arms** on each
target: **12 targets, worst |diff| 0.000e+00** on the recorded bitstrings, the objective
values and every emitted structure.

The native enters only as (a) a reported CA-RMSD and (b) columns carrying the `ORACLE_`
prefix, which `t_oracle_columns_are_named` enforces on the trajectory record itself. Every
subset rule reads bitstrings, the deployable objective value and the evaluation order, all
native-free. Three further poison tests run inside the suite:
`t_energy_ignores_the_native` (0.000e+00), `t_emission_ignores_the_native` (0.000e+00 on
21 arms x 3 targets) and the ORACLE-naming check.

One incident is on the record and is itself evidence: `stage_leak` first *crashed*, inside
LAPACK, on the NaN native - in the ORACLE column, the only branch of the recorder that reads
it. Nothing else in the chain noticed the poison. The guard is now explicit and commented.

## 7. Discipline

**No dev-24 pass was spent.** dev-24's paired SE on this family is ~0.35 A and the best arm
on the 126-target tuning instrument is −0.0022 A. `DEV_ARM` is `None`, `stage_dev` refuses to
run, and `t_no_dev_arm_and_no_benchmark` asserts both and greps the module for every
benchmark identifier. **The 60-target benchmark was not touched.**

Resource discipline: one heavy process at a time through three serialised chain scripts,
`ctypes` `GlobalMemoryStatusEx` gating at 78/90%, per-target checkpointing, and 24 MB of
recorded trajectories held out of git. The gate fired repeatedly at 90-97% while other agents
were running, and no stage was lost. Following the refinement study's recorded fault, **no
chain script was edited while a `sh` was executing it** - the two follow-on chains are
separate files.

---

## The direct answer to the mandate's question

> *Does a set generated by a misdirecting objective still contain extractable information?*

**It contains information. Consensus is not the operator that can extract it, and the reason
is not empirical.**

- The set is real: 1.21 A of internal spread, good structures 1.58 A from the start on 77.8%
  of targets, and as a *set* it beats the retrieval pool the incumbent averages by
  **−0.245 [−0.304, −0.186]**.
- The information is real: the ORACLE over it is 2.4601 A, and the objective still ranks
  within it at rho = +0.141.
- **And 87-95% of its error is common to every member**, because every member is a
  perturbation of one structure. Consensus turns independent error into accuracy at a rate
  set by `sqrt(1 − share)`, and the trajectory's independent error is 0.59-1.06 A against the
  pool's 2.09 A. The bias floor of the visited set - the best any coordinate average over it
  could reach with infinitely many members - is **3.10 A, worse than the pool's 2.81 A and
  worse than the incumbent already is.**

So the two halves do not compose, and the failure is *not* the one the mandate feared. The
trajectory is not merely the objective's own preferred region regurgitated: it is a
demonstrably better and genuinely displaced set. The failure is that **consensus and
refinement need opposite things.** Refinement narrows a set around one hypothesis, which is
what makes its members individually good; consensus needs members whose errors point in
different directions, which is what makes a set's members individually bad. The pool is a bad
set of independently-wrong members (share 0.62); the trajectory is a good set of
identically-wrong members (share 0.87-0.95). Each operator is already applied to the only set
it can use, and the seeded generator proves the trade is not negotiable: manufacturing
independent error by starting elsewhere raises the bias floor faster than it lowers the
share.

### What this licenses, and what it forbids

- **Do not** build another readout of a refinement trajectory. Argmin, medoid, the objective's
  top-64 (already computed by `refine.Archive` and reported here for the first time: 28 arms,
  best 3.1993, worst 3.5063), coordinate consensus under eleven subset rules, Boltzmann
  reweighting, time-tailing, and six unions have now all been measured, and every one pays the
  same +0.218 A per A moved.
- **Do not** raise the temperature, remove the objective, or reseed to widen the visited set.
  Measured: damping 4x costs +0.033, removing the objective costs +0.149, uniform sampling
  +0.245, seeding at pool members +0.243 - even though all four produce visited sets whose
  ORACLE is *better* than the deployable arm's.
- **Do** report the bias share of any set before proposing to average it. One cheap ORACLE
  computation predicted every pre-projection number in this file to within 0.06 A.
- **Do not** cite the synthesis' positive-phi rate as a defect of average-then-project. The
  same projection returns 0.2513 on the native trace itself, and 0.1737 multi-start.
- **The quantity to attack is the bias floor, not the spread.** The pool's is 2.81 A and the
  incumbent reaches 3.20 after projection. Nothing in this file, nothing in `s9/refine.py`'s
  28 arms and nothing in `s9/synth.py`'s 23 operators moves either. A set with a lower bias
  floor is a **retrieval** problem - not a refinement problem, not a search problem, not a
  readout problem, and not a projection problem.
