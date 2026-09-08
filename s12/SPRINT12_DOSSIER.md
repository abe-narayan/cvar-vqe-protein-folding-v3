# SPRINT 12 — RESEARCH DOSSIER

**Question put to the sprint:** can this project be transformed into a substantially better
peptide-structure predictor, ideally reaching < 2.0 Å mean CA-RMSD on unseen structures?

**Answer:** No — and the sprint establishes *why* with numbers rather than assertions, on a
tighter and better-controlled instrument than any previous sprint. It also produces one
architectural law that explains four sprints of null results, two corrections to the
project's own ceiling arithmetic, and one clearly-specified direction that could break the
ceiling but was not tested here because it needs a data modality the pipeline does not yet
read.

Ten agents ran concurrently on the 126-target tuning instrument. **benchmark60 was never
read, dev24 was never spent, and no tracked file was modified.** Every agent's numbers were
produced through one shared instrument (`s12/instrument.py`) that reproduces the project's
five pinned constants exactly, and that instrument was itself audited by an adversarial
agent before its outputs were believed.

| findings file | agent |
|---|---|
| `s12/coord_FINDINGS.md` | coordinator — instrument, integrity audits, corpus/null experiments |
| `s12/fail_FINDINGS.md` | failure-target forensics (748 lines) |
| `s12/forensics_FINDINGS.md` | independent failure forensics + structural retrieval (566 lines) |
| `s12/obj_FINDINGS.md` | objective and error structure |
| `s12/agg_FINDINGS.md` | set aggregation, terminal operator, learned set decoder (808 lines) |
| `s12/asm_FINDINGS.md` / `s12/assembly_FINDINGS.md` | two independent multi-piece assembly studies |
| `s12/key_FINDINGS.md` | retrieval keys / local-conformation prediction (680 lines) |
| `s12/dir_FINDINGS.md` | per-pair error-direction head (580 lines) |
| `s12/lit_FINDINGS.md` | 2026 literature and foundation models |
| `s12/vq_FINDINGS.md` | genuine VQE/CVaR combinatorial role |
| `s12/adv_FINDINGS.md` | adversarial audit of everything above |

---

# I. EXECUTIVE CONCLUSION

**What limits the system is not the objective, not the operator, and not the library. It is
that the information needed to tell a near-native candidate from a plausible wrong one is
not present in any channel this project can compute — and the sprint now has a decisive
test of that, not an inference from a pile of failures.**

The previous sprint (S9-8) concluded in-band discrimination was "informationally
impossible" from a monotone-decreasing capacity curve at fixed sample size. The corrections
brief correctly judged that argument insufficient: the learning curve was still rising, and
the untested class was *learned aggregation over the full candidate-vs-objective deviation
map*. That class has now been built and tested — a set-transformer over the
(75 × n_pairs × 18) **signed** deviation tensor plus a (75×75) pairwise-RMSD attention bias,
trained on the aggregate structure's RMSD (differentiable through the coordinate average),
with the 18 catastrophic targets held out of training entirely.

| | emitted CA-RMSD |
|---|---|
| production synthesis | 3.203 |
| **learned set decoder (28 arms, 2 seeds)** | **3.184, d = −0.019 [−0.058, +0.020], 63W/63L, reverses on drop-top-10** |
| learning curve, n = 8 → 16 → 32 → 64 → 75 training targets | 3.043 → 3.049 → 3.046 → 3.036 → **3.026 (FLAT)** |
| identical harness with a **leaked** RMSD label | **2.534 at n = 8**, 2.200 at n = 32, 2.146 at n = 75 |

The harness finds a signal that is definitionally present, immediately, from eight targets.
Given only deployable features it finds nothing, and adding twelve times more training data
does not change that. **The problem is signal-limited, not sample-limited.** That closes
correction C2's remaining hypothesis properly, by measurement.

Three further results bound what is left:

1. **The distance channel's ceiling is ~1.95-2.0 Å, and only with a co-optimised terminal.**
   Scoring the K=500 pool with the *native* distance matrix emits **2.395 Å at the shipped
   m=75** and **~1.95-1.99 Å at m=5**, the cardinality optimal for a perfect objective. The
   adversarial audit forced this correction and it is the least discouraging result of the
   sprint: the statement "no objective upgrade reaches 2.0 Å" is only true *at m=75*. The
   distance line can reach the target — but only if the objective becomes perfect, and no
   deployable objective improvement was found.
2. **The terminal operator consumes a set MEAN, not a set BEST**, which explains why four
   sprints of ranking work measured null (§ III).
3. **On the 18 catastrophic targets, the entire sequence-conditioned machinery is worse
   than nothing.** A pipeline given a *random* 500 windows and a *random* 75 emits 5.425 Å
   on them; the shipped system emits 6.019 Å. On the other 108, sequence conditioning is
   worth 1.004 Å. There is no deployable signal that tells the two regimes apart — three
   independent routers were built and all three are null.

**No deployable change measured in this sprint improves accuracy.** Nine substantive
interventions were built and all nine are negative with their nulls (§ VIII).

**One direction does reach the target, and it needs a data modality the pipeline does not yet
read.** Supplied with backbone torsions at the accuracy a chemical-shift predictor delivers
(σ = 12°, all residues), simply *building* the chain emits **1.486 Å with 86% of targets under
2.0 Å**, and the failure class goes from 6.019 Å to 1.80 Å (§ XVII). The same restraints used
as a *filter* over the retrieval pool are worth only −0.24 Å — so this is a different
architecture, not an extra channel. Its make-or-break parameter is measured: it needs ~90%
residue coverage, and turns negative below 50%.

The honest outcome is the brief's Outcome D/E — *a demonstration that the current
representation is insufficient, with precise ceilings* — plus a specified, priced path out.

---

# II. WHAT CHANGED MY MIND

### Survived
- **BLOSUM retrieval beats random** (S7-12). Confirmed model-free: pool best 1.711 vs 1.815,
  pool mean 4.453 vs 4.817, both CIs excluding zero.
- **AMBER is a validity stage, not an accuracy stage** (S8-12). Untouched; +0.021 Å cost.
- **The CA–CA distance matrix is the right objective representation.** Contact maps, Gram
  eigenvectors, scale-invariant L1, per-shell standardisation, error-covariance Mahalanobis
  and ORACLE CA-level local geometry all lose to it; oracle local geometry *degrades* the
  distance matrix when fused (2.218 → 2.416).
- **The fibril/lasso enrichment in the failure class** (10/18 vs 6/108, Fisher p = 1.2e-06).
  Attacked three ways by the adversarial agent — priced against the whole space of header
  flags actually searched, leave-one-out over targets and classes, and independent
  re-derivation of the labels from the headers — and did not move.
- **The library is saturated by depth.** Reconfirmed from a second direction: coverage holes
  open at 7–8 residues and scale at −0.432 Å per decade of library size.

### Broke or was corrected
- **"The pool cap is 2.406 Å"** (correction C1) — confirmed broken and the reason quantified:
  argmin → average is worth 0.406 Å raw / 0.250 Å emitted, and it is *already banked*.
- **"The objective's own geometric optimum is 4.62 Å, worse than the 3.20 Å the pipeline
  emits."** The number is **refuted (4.62 → 3.53)** by the adversarial agent using a better
  decoder; the conclusion survives at a smaller margin.
- **"r_sep = 0.196"** is estimator-dependent by a factor of two (→ 0.37–0.39 under pooled,
  cross-validated partialling). The *contrast* it rests on survives every estimator.
- **The two-piece assembly floor of 0.53/0.32 Å** — refuted as stated. The ideal-geometry
  rebuild floor is already 0.347 Å, and a matched synthetic bank with i.i.d. Ramachandran
  torsions and zero fragment content reaches 0.328 Å against the real library's 0.309 Å.
  **The library is worth 0.019–0.028 Å as assembly pieces.**
- **The literature agent's ORACLE ABEGO retrieval key was circular** — its key was the bin
  string of the universe's true-best window, which retrieves that window by construction.
  The honest oracle gives FAIL18 pool-best 1.768, not 1.640.
- **My own hypothesis** that the distogram's sequence-separation component was dead weight
  was **refuted by the experiment I commissioned**: separation-partialled scoring emits
  3.642 Å, worse than shipped *and worse than a random score*. The separation profile
  carries 75% of the objective's damage, not 0%.
- **"Perfect filter at m=25 = 1.644 Å"** is a perfect-*RMSD* filter (an oracle on the
  answer). A perfect *distance objective* at m=25 emits **2.123 Å**. The distogram line's
  true ceiling is **1.95 Å**, not 1.64.

### Strengthened
- **S7-2 (more protein-fragment data makes the predictor monotonically worse)** now
  reproduces in a *second, independent model family*: a local-conformation classifier goes
  from 0.690 to 0.618 accuracy when given 8× more fragment data. And the *mechanism* is now
  measured directly (§ III).
- **"Consensus is outlier avoidance, not a nativeness signal."** Confirmed per-target and
  given a mechanism: on the FAIL18 the near-native band is the pool's *outlier* set
  (centrality AUC 0.252), so every centrality-flavoured method is structurally guaranteed to
  fail there.
- **"The predictor ≈ typicality."** The K=500 pool's own mean distance profile — no model, no
  labels, no sequence — matches the trained distogram's profile and *beats it* at separations
  ≥ 7.

---

# III. THE CURRENT INFORMATION BOTTLENECK

Three quantitative statements, each from an independent agent, that together locate it.

### 1. The sequence→structure channel is weak, and 80% of the library is drawn from a corpus where it barely exists

Model-free, every same-length pair in each corpus, sequence identity against CA-RMSD
(`s12/coord_infochannel.py`):

| corpus | ρ(identity, RMSD), n = 9…16 |
|---|---|
| peptide database (787 chains) | **−0.17 to −0.29** |
| protein fragments (6,003, from 13,751 PDB proteins) | **−0.030 to −0.040** |

A 12-residue stretch cut out of a folded protein has the conformation its tertiary context
imposes; its own sequence did not choose it. Inside a target's own retrieval universe —
which is ~80% fragments — the deployable channel measures **ρ(BLOSUM, true RMSD) = −0.066**,
exactly the mixture those two numbers predict, and **+0.001 on the FAIL18**.

This is also the mechanism behind S7-2: the fragment corpus is the wrong conditional
distribution both for *training* the distance predictor and for *populating* the pool.

### 2. The whole architecture's sequence information is worth 0.776 Å — and it is negative where it is needed

The 2×2 that had never been run (`s12/coord_null.py`): BLOSUM-500 × distogram-75 (shipped)
against random-500 × random-75 (told nothing about the target sequence), 3 seeds on the
random arms, everything else identical.

| arm | retrieval | filter | emitted |
|---|---|---|---|
| shipped | BLOSUM-500 | distogram-75 | **3.213** |
| no filter | BLOSUM-500 | random-75 | 3.597 |
| no retrieval | random-500 | distogram-75 | 3.358 |
| **blind** | random-500 | random-75 | **3.989** |

| contrast | value |
|---|---|
| both channels (blind → shipped) | **−0.776 [−1.020, −0.543]**, 90W/36L |
| distogram alone (blind → no-retrieval) | −0.630 [−0.829, −0.431] |
| BLOSUM alone (blind → no-filter) | −0.391 [−0.497, −0.283] |
| distogram *given* BLOSUM | −0.385 |
| BLOSUM *given* the distogram | −0.146 |

The two channels are strongly overlapping (0.391 + 0.630 = 1.021 jointly delivering 0.776),
consistent with the distogram being largely typicality plus a weak sequence signal. And the
subgroup split is the sharpest result in the sprint:

| subgroup | blind | shipped |
|---|---|---|
| other 108 | 3.749 | **2.745** (sequence worth −1.004) |
| **FAIL18** | **5.425** | **6.019** (sequence worth **+0.594 — actively harmful**) |

**The system's sequence conditioning helps by 1.0 Å where the target is in-distribution and
hurts by 0.6 Å where it is not.** The retrieval-key agent found the same thing from the other
side — on the FAIL18 a content-free capacity null beats the real predictor.

### 3. The terminal operator consumes a set MEAN, so ranking improvements are structurally invisible

The adversarial agent's attribution experiment (`s12/adv_operator.py`, 126 targets, 882
perturbations): replace the *j* worst-scoring members of the shipped top-75 with the pool's
*j* oracle-best members, so at j ≥ 1 the pool's best structure is present **and at rank 1** —
exactly what a perfect objective would hand the operator on the same candidate set.

| j | set best | set mean | argmin emits | **average emits** |
|---|---|---|---|---|
| 0 (production) | 2.306 | 3.551 | 3.454 | **3.048** |
| 1 (perfect rank-1) | 1.711 | 3.533 | **1.711** | **3.023** |
| 25 | 1.711 | 3.219 | 1.711 | 2.580 |
| 75 (perfect filter) | 1.711 | 2.686 | 1.711 | **1.963** |

    d_avg = 0.039·d_set_best + 1.162·d_set_mean − 0.028      R² = 0.893
    set_mean alone R² = 0.891        set_best alone R² = 0.187

A **perfect rank-1 decision is worth −1.743 Å through argmin and −0.029 Å through the
average the system ships** — below the sprint's own 0.03 Å ignore threshold. At m=75 the
operator emits at a near-fixed quantile (0.20–0.25) of whatever distribution it is given.
The aggregation agent reproduced the law independently (slope 1.157 vs 1.162).

**The law's own null, run at the coordinator's request, narrows it correctly.** Re-fitted on
three further perturbation designs (random subsets, 22 different scorers, and the cardinality
axis), the **set-mean dominance survives everywhere** — the coefficient on `d_set_mean` is
0.68–1.16 and set-mean alone explains 55–89% of the variance. But the coefficient of 0.039 on
`d_set_best` is **design-specific**; out of design it is 0.18–0.30. The correct statement is
therefore narrower: *the operator is almost blind to a set-best improvement carried by one or
two members out of 75 that leaves the bulk unmoved* — which is exactly the oracle-insertion
case and exactly the assembly case, and is confirmed there by direct measurement (j=1 is
worth −0.029 Å) rather than by a fitted coefficient. The fixed-quantile statement likewise
holds only at fixed cardinality.

**Does this invalidate four sprints of null ranking results?** The question was put to the
adversarial agent and answered by measurement. The terminal **compresses** ranking effects by
a factor of ~1.6–2.0 (perfect ranking: −1.743 Å at m=1 against −1.085 Å at m=75; shipped vs
random: −0.786 against −0.385) but **it does not flip a sign and does not send a real effect
to zero.** Applying the measured compression: a ranker that measured 0.00 Å at m=75 was worth
at most ~0.06 Å through argmin, and no discarded ranker in the record reached the 0.10 Å at
m=75 that would have been worth 0.20 Å through argmin. **S9-8's closure survives this
objection.** The aggregation agent's independent check agrees: its decoder's oracle control is
loud through every terminal (argmin −1.09, avg-25 −0.56, weighted-75 −0.90, 123–125/126 wins),
so a real gain would not have been invisible.

What does change is **future practice**. Improving the ranker moves the optimal cardinality,
so a fixed-m evaluation *understates* a ranking improvement: going from the shipped objective
to ρ = 0.7 is worth −0.967 Å at fixed m=75 and **−1.346 Å at the co-optimised m — a 39%
understatement.** Any future ranker must be scored on the (ranker × m) surface, not at a
point. And the design rule stands: **the filter's job in this architecture is to raise the
mean quality of the retained set, not to find the best member.**

### The corollary that sets the target

Emitted-average ≈ 1.13 × (set mean) − 0.96. **2.0 Å emitted requires a retained-set mean of
about 2.62 Å.** The pool mean is 4.45, the shipped filter delivers 3.55, a perfect filter
delivers 2.69. The filter already captures ~51% of the available set-mean improvement.

---

# IV. THE 18 FAILURE TARGETS

Two agents attacked them independently with different instruments and converged.

### What they are

| flag (PDB header text only) | FAIL18 | other-108 | Fisher p |
|---|---|---|---|
| fibril / amyloid / steric zipper | 6/18 | 4/108 | 5.5e-04 |
| lasso peptide | 4/18 | 2/108 | 3.8e-03 |
| **fibril OR lasso** | **10/18 (56%)** | **6/108 (5.6%)** | **1.2e-06** |
| membrane / micelle / TFE | 5/18 | 45/108 | n.s. — *depleted* |
| disulfide / LINK constrained | 8/18 | 52/108 | n.s. |

The negative controls matter: membrane peptides are the largest class in the database and
are *depleted* among the failures. The enrichment is specific to two classes whose
conformation is imposed by something the sequence does not encode — a covalent thread
through a macrolactam ring, or a cross-β lattice.

### What is not wrong with them

- **Not a label problem.** All deposited models parsed: ensemble spread 0.942 Å on FAIL18
  vs 0.946 on length-matched controls — indistinguishable. Switching to a best-of-ensemble
  convention recovers 0.178 Å on FAIL18 and 0.217 on controls; paired difference −0.039
  [−0.194, +0.113], **null**. Exactly one target (3BTB, spread 4.27 Å, in fast exchange with
  its partner) is a genuine reference problem.
- **Not a retrieval-coverage problem.** The near-native band's in-pool fraction is 0.079 on
  FAIL18 against 0.073 on the other 108, **p = 0.96**. BLOSUM is at chance on *both* groups.
- **Not a projection problem.** Feeding the pool's single best member through the real
  synthesis path emits 2.285 Å against a 2.284 Å pool best — the geometry stage is lossless.

### What is wrong with them

The native sits at the **78th percentile** of the shipped score over its own pool (32nd on
matched controls). ρ(score, true RMSD) is 0.107 vs 0.734, and *negative* inside the band.
The distogram's MAE on these natives is 4.37 Å vs 1.76, concentrated at separations 5–8
(5.13 vs 2.03).

The mechanism is a **coherent global shape error**: the top-75's radius of gyration tracks
the distogram's *implied* rg to within ~0.5 Å on essentially every target, and both miss the
native **in both directions** — 2BFI native rg 11.9 Å against a top-75 at 5.7; 2MQ2 5.4
against 10.3. The rg calibration slope is 0.238. The score is not noisy; it faithfully
executes a confidently wrong predicted shape.

### Reachability

- **Lasso is not reachable by linear-window retrieval.** 6/6 lasso peptides in tuning126
  emit ≥ 6.10 Å (mean 6.71). A mean of **2.3 windows out of ~10,000** lie within 2.0 Å of the
  native and three of the six have **zero**. This is a coverage limit that library growth
  cannot fix, and it costs ≈ 0.17 Å of the tuning mean. It should be declared out of scope.
- **Fibril is reachable in principle and lost at the objective.** 371 windows per target lie
  within 2.0 Å (3.5% of the universe) and the universe best is 1.20 Å; a two-parameter
  constant-(φ,ψ) extended chain scores 0.27–2.41 Å on six of them, beating the best of all
  21,547 windows on 2BFI. The compactness-preferring objective discards them (rg-fit AUC
  0.002–0.023).
- **Every conventional in-band signal is anti-correlated on the FAIL18**: distogram score AUC
  0.301, consensus/centrality 0.252, Legacy energy 0.380, rg-agreement 0.293. The near-native
  band *is* the pool's outlier set there.

---

# V. LITERATURE FINDINGS THAT MATTER

Full survey and 24-row method table in `s12/lit_FINDINGS.md` / `s12/results/lit_methods.json`.

1. **The isolated 9–16-residue NMR peptide is a documented blind spot of the entire
   AlphaFold3-like family**, not a place where this project is behind the field. AF3 fails
   10/12 PDB lasso peptides (8% success on the RODEO set, *Nat Commun* 2025); the family
   fails 67–73% of monomeric NMR structures. ESM-2's training distribution barely contains
   peptides (2.8% of UniRef50 is < 50 residues). All of these need a GPU and none can run here.
2. **The field has retrieved fragments by predicted local backbone conformation since 2011**
   (Rosetta `torsion_bin_probs`), and it is the architecture of both peptide predictors that
   beat this system — PEP-FOLD's 27-state structural-alphabet profile and APPTEST's
   distance+torsion head. **This project retrieves by BLOSUM62 sum alone.** The sprint tested
   that gap directly and it does not close (§ VII, key agent) because the predictor is not
   accurate enough on exactly the targets that need it.
3. **APPTEST reaches 1.96 Å on 42 short peptides** with deployable selection — the number to
   beat, and a reason to believe sub-2.0 Å is attainable by *some* architecture at this length.
4. **Chemical shifts are deposited in the BMRB for 79 of the 126 tuning targets (63%), and
   for 10 of the 18 failures.** TALOS-N converts them to φ/ψ at ~12° accuracy. This is the
   one channel identified that is (a) genuinely orthogonal to sequence, (b) realistically
   available to a predictor rather than extracted from the answer, and (c) not tested here.
5. **Amyloid detectors (WALTZ-DB 2.0, CORDAX) are trained on in-vitro aggregation assays and
   have never seen a PDB coordinate** — zero structural leakage, CPU-free. They are the only
   clean route to routing the fibril class.
6. **Do not fund**: 3Di/ProstT5/SaProt structural retrieval (3Di encodes *tertiary*
   neighbours and is degenerate for a 13-mer), AF3/Boltz/Chai/ESMFold (GPU, and they fail
   this target class), PepFlow/PPFlow/PepGLAD (receptor-conditioned *design*, wrong task),
   RFdiffusion/Genie2/Chroma (designability prior is wrong-signed for floppy peptides),
   MDGen (tetrapeptides only).

---

# VI. ARCHITECTURE LANDSCAPE

Every architecture in the brief's list A–J, priced against what the sprint measured.

| # | architecture | status | evidence |
|---|---|---|---|
| **A** | retrieval → filter → synthesis (incumbent) | **3.204 Å; still the best deployable thing measured** | — |
| **B** | retrieval → multi-piece assembly | **built twice, both negative deployed.** Raises the achievable space (deployable random-shortlist space best 1.438 vs pool 1.711, −0.255 [−0.340,−0.168]) and the ORACLE top-25 emits 1.265 — but the deployable emitted answer is 3.397/3.570, i.e. +0.19/+0.22 **worse**. The oracle advantage also fails a matched synthetic-bank capacity null (0.309 real vs 0.328 content-free) | `asm_`, `assembly_` |
| **C** | template retrieval → template-conditioned refinement | **null.** 307 ideal templates per target appended to the pool: emitted +0.0042 [−0.048, +0.060]. Templates beat the whole 21,547-window universe on only 5/126 and are 1.7 Å *worse* than the library on the FAIL18 | `coord_templates` |
| **D** | foundation-model proposals → selection | **not runnable and not indicated.** GPU-only; documented to fail this exact target class | `lit_` |
| **E** | diffusion/generative proposals → selection | **not indicated.** Receptor-conditioned design models solve a different problem; designability priors are wrong-signed for flexible peptides | `lit_` |
| **F** | distance/torsion prediction → geometric sampling | **capped.** The objective's own best 3-D decoding is 3.53 Å (adversarially improved from a reported 4.62), *worse* than what retrieval+averaging emits. Sampling against this objective cannot beat restricting to real peptide geometry | `obj_`, `adv_` |
| **G** | hybrid retrieval + generative ensemble | **error correlation kills it.** Assembled candidates' error correlates r = 0.928 with the retrieval answer; all 19 retrieval keys correlate r ≥ 0.87 with BLOSUM's per-target error. There is no decorrelated second opinion to fuse | `assembly_`, `forensics_` |
| **H** | direct sequence → structural ensemble | **the channel is 0.776 Å wide in total** and negative on the failure class (§ III.2) | `coord_null` |
| **I** | multiple generators → learned set aggregation | **the sprint's decisive negative.** Flat learning curve against a loud oracle control (§ I) | `agg_` |
| **J** | generator → VQE/CVaR set selection → Legacy → projection → AMBER | **runs genuinely, contributes ~0 accuracy, and now has a defensible role as a certified optimiser benchmark** (§ IX) | `vq_`, `agg_` |

**The incumbent architecture is not mis-specified.** It is a well-chosen local optimum for
the information it has. The one structural criticism that survives is that its terminal
operator makes it blind to its own improvement path — but the sprint also shows there is no
improvement to be blind to, at present.

---

# VII. EXPERIMENTS RUN

Every row: hypothesis → control → result on the 126-target instrument. All emitted CA-RMSD
unless stated. Negative CI values favour the intervention.

| # | hypothesis | control used | result | verdict |
|---|---|---|---|---|
| 1 | Learned aggregation over the full **signed** n×n deviation map recovers in-band skill (correction C2) | leaked-label oracle harness; feature permutation; FAIL18 held out; 3 terminals | 3.184 vs 3.203, **−0.019 [−0.058,+0.020]**, 63/63, reverses on drop-10; learning curve flat 3.043→3.026 while leaked label reaches 2.534 at n=8 | **NEGATIVE, decisive** |
| 2 | A perfect distance objective reaches 2.0 Å | native matrix through the real path | **2.395 at m=75; 1.948 at the optimal m=5** | **CEILING ESTABLISHED** |
| 3 | The terminal operator can cash a ranking improvement | oracle rank-1 insertion, 882 perturbations | −1.743 Å through argmin, **−0.029 Å through the average**; law R² 0.893 | **LAW ESTABLISHED** |
| 4 | Optimal cardinality m is a function of objective quality | 14 objective levels × 15 cardinalities | m* = 500 → 75–110 (today) → 20 → **3–5** (perfect). Crossover at MAE ≈ 0.85; shipped is at 2.34 | **CONFIRMED; m=75 correct today** |
| 5 | Sequence-predicted local conformation (ABEGO) as a **retrieval key** | shuffled-label, composition, and Ramachandran-marginal capacity nulls | predictor 0.690 (peptide-trained) clears the 0.60 gate on the 108, **fails on FAIL18 (0.517)**; emitted **+0.015 [−0.114,+0.137]**; on FAIL18 the *composition null* (−0.646) beats the real predictor (−0.510) | **NEGATIVE** |
| 6 | Per-pair **error-direction** head corrects the objective | pre-registered accuracy gate; matched-accuracy random corruption; separation-only null | head 0.687 LFO accuracy **clears** the 0.684 gate, and emits **+0.291 [+0.124,+0.455]**; oracle signs corrupted to the *same* accuracy emit −0.142 | **NEGATIVE, mechanism identified** |
| 7 | Multi-piece assembly is a better generator | matched synthetic i.i.d.-Ramachandran bank; random-piece draws | oracle 0.309 real vs **0.328 content-free**; deployable **+0.193/+0.218 worse** | **NEGATIVE** |
| 8 | Ideal-geometry templates fill library holes | matched random-window append | **+0.0042 [−0.048,+0.060]** | **NEGATIVE** |
| 9 | Weighting retrieval toward the peptide corpus (7× the sequence→structure channel) | selectivity-matched random sub-universe | peptide-only **−0.020 [−0.097,+0.052]**, drop-10 +0.063; fragment-only **+0.092**; pep vs selectivity control −0.085 [−0.162,−0.012] | **NEGATIVE (direction real, size ~0)** |
| 10 | Separation-partialled scoring (coordinator's hypothesis) | random-score reference | **3.642, +0.436 [+0.276,+0.627]** — worse than a random score | **REFUTED** |
| 11 | ESM-2 contact agreement is a decorrelated in-band channel | foreign length-matched contact map | net +0.026…+0.043 AUC, CIs spanning zero on FAIL18 | **NEGATIVE** |
| 12 | Alternative objective representations (contact map, Gram, orientation, Mahalanobis, local geometry) | shipped objective, matched MAE | all lose; oracle local geometry *degrades* the distance matrix when fused | **NEGATIVE** |
| 13 | Multi-hypothesis / mode selection | 192 native-free rules × 3 clusterings × k=2–5 | oracle headroom 0.367 Å; **best deployable rule +0.013** | **NEGATIVE** |
| 14 | A deployable router (regime, or cardinality) exists | permutation null, LFO | AUC 0.558 vs null 0.511; m-router −0.0021 [−0.062,+0.055] vs a 0.245 Å oracle router | **NEGATIVE ×3** |
| 15 | Error structure matters more than MAE | 15 corruption families × 6 MAE levels × 3 seeds | at the shipped MAE, emitted ranges 2.24–2.85 Å across error shapes; **the real distogram (3.05) is worse than every synthetic corruption**; i.i.d. noise at the same MAE emits 2.443 (−0.762, 100W/26L) | **CONFIRMED** |
| 16 | The sequence→structure channel differs by corpus | length-matched, all pairs | peptides ρ −0.17…−0.29, fragments **−0.030…−0.040** | **CONFIRMED** |
| 17 | The architecture's total sequence information | random-500 × random-75 blind pipeline | **−0.776 [−1.020,−0.543]** overall; **+0.594 (harmful) on FAIL18** | **CONFIRMED** |
| 18 | A genuine combinatorial problem exists for CVaR-VQE | exhaustive enumeration; shuffled-geometry null | 1,008 certified QUBO instances; SA finds the optimum on 90%, greedy is exact on 40.5%; **closing the optimality gap costs +0.022 Å** | **PROBLEM REAL, ACCURACY ZERO** |
| 19 | Training-set containment leaks explain results | per-target price, FAIL18 split | 4/126 verbatim leaks; ρ(containment, RMSD) = −0.129 (p 0.15); FAIL18 are the *cleanest* targets | **NEGATIVE (integrity defect only)** |
| 20 | A fresh, adequately-powered benchmark can be built | RCSB census + full gate/screen pipeline | date-fresh supply 20–28; containment-fresh yield **16**, of which **10 are amyloid fibrils** | **NOT AVAILABLE** |

---

# VIII. NEGATIVE RESULTS, STATED PLAINLY

Nine substantive interventions were built and measured this sprint. **All nine are negative.**
None is reported as promising; none warrants a dev24 pass; dev24 was not spent.

1. Learned set decoder over the signed deviation map — flat learning curve.
2. Torsion-bin / ABEGO retrieval key — predictor fails the gate exactly where it is needed.
3. Per-pair error-direction head — clears its accuracy gate and emits worse.
4. Multi-piece assembly (two independent implementations) — oracle dies to a capacity null,
   deployable is worse.
5. Ideal-geometry template augmentation — null.
6. Peptide-corpus retrieval weighting — null.
7. ESM-2 contact agreement as an in-band channel — dies to a foreign-map null.
8. Separation-partialled scoring — worse than a random score.
9. Multi-hypothesis mode selection, alternative objective representations, and all three
   routers — null.

**Two failure mechanisms recur and should be treated as laws of this problem:**

- **Learned correctors inherit the error structure of what they correct.** The direction
  head's mistakes are *coherent* (0.089 of squared error in the global expand/contract mode
  vs 0.053 for i.i.d. at the same rate; effective rank 6.62 vs 7.28). Coherent wrongness
  moves the coordinate average; i.i.d. wrongness cancels. At *identical* 0.688 accuracy, a
  random corrector emits −0.142 Å and the trained one +0.291 Å — a 0.43 Å swing with opposite
  signs, decided entirely by error correlation. Feature pruning does **not** decohere it.
- **Predictor accuracy does not price emitted structure quality.** Three independent metrics
  now fail: MAE (S7-3, S7-11), in-band rank correlation (S9-8, S11), and per-pair sign
  accuracy (this sprint — thirteen accuracy points over the separation-only null bought
  exactly 0.000 Å). The only intermediate quantity that has ever tracked emitted RMSD here is
  the oracle arm itself.

---

# IX. THE BEST ARCHITECTURE THE EVIDENCE SUPPORTS

The evidence does not support replacing the pipeline. It supports keeping it, fixing three
integrity defects, giving the mandated components scientifically honest roles, and declaring
a domain of applicability.

```
sequence
  |
  +- STAGE 0  CLASS SCREEN  (new)  -- lasso / steric-zipper detection from sequence
  |                                   (WALTZ-DB / CORDAX-style, aggregation-assay-trained,
  |                                    zero structural leakage).  Does NOT re-route the
  |                                    structure; it ATTACHES A CONFIDENCE LABEL and, for
  |                                    the lasso class, declares the target out of scope.
  |
  +- STAGE 1  RETRIEVE     K=500 BLOSUM62 over out-of-fold peptides + fold fragments.
  |                        UNCHANGED.  (Peptide-weighted retrieval is directionally right
  |                        and worth ~0.02 A; not worth a config change.)
  |
  +- STAGE 2  FILTER       the learned distogram's Bayes-risk score, top m.
  |                        m = 75 TODAY, and m is now a DECLARED FUNCTION of objective
  |                        quality: m* = 75 at MAE 2.34, 20 at MAE 1.2, 3-5 at MAE 0.
  |                        Holding m=75 costs +0.386 A if the objective becomes perfect.
  |
  +- STAGE 3  SYNTHESISE   medoid superposition -> coordinate average -> multi-start
  |                        projection onto ideal geometry, ramah lambda=0.3.  UNCHANGED and
  |                        demonstrably at its native-free optimum (30 arms, none better
  |                        by more than 0.004 A).
  |
  +- STAGE 3b CVaR-VQE     SET-LEVEL SUBSET SELECTION over 16 medoid hypotheses:
  |                        F(S) = risk(mean_c d_c) + mu * Cons(S), exact QUBO
  |                        H(x) = x'Mx/k^2 - 2b'x/k + const.  1,008 certified instances,
  |                        all 2^16 enumerated.  REPORTED AS AN OPTIMISER BENCHMARK, not as
  |                        an accuracy method -- closing the classical optimality gap costs
  |                        +0.022 A because corr(H, true RMSD) = +0.077.
  |
  +- STAGE 4  LEGACY       late refiner over the consensus neighbourhood.  UNCHANGED;
  |                        19,000x cheaper than AMBER; measured contribution ~0.
  |
  +- STAGE 5  AMBER        ff14SB/GBn2 restrained relaxation, k=10, converged.
                           VALIDITY stage; costs +0.021 A and earns its place by removing
                           ~10^4 kcal/mol of builder strain.
```

**Why the four mandated components stay, and in these roles.** Each is genuine, each is
measured, and none is forced into a task a classical method does better:

- **VQE/CVaR** now has a real combinatorial problem — subset selection where a candidate's
  value depends on the rest of the set, because the operator is an average and the risk of a
  mean is not the mean of risks. The Hamiltonian is exact (QUBO optimum = true optimum on
  100% of m=2 instances, brute-force verified for N ≤ 16). **Classical optimisation wins**:
  simulated annealing finds the exhaustive optimum on 90% of 1,402 problems and greedy is
  exactly optimal on 40.5%. That is reported as the result. CVaR's one defensible
  contribution is measured rather than asserted: it prevents the collapse to the argmin that
  α=1 produces (0.076 bits of state entropy), returning a distribution over hypotheses.
- **Legacy** is the only physics cheap enough to sit inside a search loop and is retained
  where its skill was measured — a neighbourhood that is already tight.
- **AMBER** is retained as the validity stage its own ablation supports, with its accuracy
  cost reported rather than hidden.
- **The distogram** remains the filter; the sprint bounds its whole line at 1.95 Å.

**Three integrity fixes the sprint found and the architecture should carry:**

1. Write a `folds()` digest into every model checkpoint and refuse to load on mismatch
   (hazard H1 is currently prevented only by a directory rename).
2. Replace the fold-split screen's longer-normalised `identity` with `containment` **plus a
   verbatim-substring test** — 4/126 targets currently have a verbatim copy of themselves in
   their own distogram's training set. Note that a *containment threshold alone* is unusable
   at this length: a random composition-matched sequence scores 0.56–0.63 against the
   fragment bank.
3. Record, with every reported number, which terminal operator it was scored through.

---

# X. IMPLEMENTATION PLAN

| step | module | work | gate |
|---|---|---|---|
| 1 | `core/predict.py`, `core/data.py` | folds digest in checkpoints; containment + substring screen in `clusters()` | must reproduce the pinned manifests bit-identically or the fold-repin hazard fires |
| 2 | `core/pipeline.py` | make `m` a declared function of a measured objective-quality statistic; default unchanged at 75 | the m-router is null, so this ships as a **constant with a documented derivative**, not as an adaptive rule |
| 3 | `core/bench.py` | every results file records the terminal operator and reports both argmin and average arms | no science change |
| 4 | `s12` → `core` promotion | the QUBO instance builder and its certified optima, as a permanent optimiser benchmark for `core.quantum` | Hamiltonian must reproduce the true objective to numerical precision (verified) |
| 5 | new | class screen (aggregation-assay detector) as a **confidence label**, not a router | detector precision > 0.5 measured out-of-fold before it may change any output |
| 6 | **the only accuracy work worth funding** | BMRB chemical-shift ingestion → TALOS-N-style φ/ψ restraints → an additional retrieval/filter channel | see § XIII.1 |

---

# XI. VALIDATION PLAN

**This is the sprint's hardest constraint and it must be stated before any future claim.**

The three existing instruments spend **all 204** identity clusters of 9–16-residue peptides
in the frozen corpus. Outside it, RCSB holds 700 single-protein-entity 9–16-mers in total;
28 released since 2026-01-01; 146 absent from this corpus entirely. Building the
containment-fresh set out in full (downloaded to `s12/newpdbs/`, never to `pdbs_ext/`) gives:

    downloaded 133 -> pass structural gates 40 -> unique sequences 30
    -> clean of leakage 19 -> one per identity cluster -> 16 targets
    of which 10 are amyloid fibril segments, 4 X-ray designed oligomers, 2 solution NMR

**A fresh, adequately-powered, distribution-matched benchmark does not exist, because the
world does not contain one.** Any future validated claim must name which of these it rests on:

1. **dev24**, one pre-registered pass — the only clean confirmation available, spendable once.
2. **A length-extrapolation benchmark** from the 266 unused clusters of 8-mers and 17–26-mers
   (382 peptides — adequately powered, but a different target distribution; it tests
   generalisation across length rather than within it).
3. **Nested cross-validation on tuning126**, with the number of architecture decisions made
   after seeing tuning results declared explicitly.

Leakage discipline for anything new: verbatim-substring test plus identity < 0.6 (null
0.41–0.45), never a bare containment threshold (null 0.56–0.63). Every learned component
leave-fold-out on the pinned folds. Every oracle arm named and labelled.

**Selection pressure incurred this sprint:** twenty hypotheses were tested on tuning126 and
all twenty are negative, so no architecture decision was made *on* the instrument. That is
why the instrument is still usable.

---

# XII. EXPECTED PERFORMANCE

| scenario | expected tuning126 mean | basis |
|---|---|---|
| ship the incumbent unchanged | **3.20–3.24 Å** | measured |
| incumbent + every integrity fix in § IX | **3.20–3.24 Å** | the fixes change reproducibility, not accuracy |
| declare the lasso class out of scope (6 targets) | **~3.04 Å** on the remaining 120 | arithmetic on measured per-target values; **not** an accuracy improvement, a scope change |
| a perfect distance objective, m held at 75 | **2.395 Å** | measured |
| **a perfect distance objective at its co-optimised m = 5** | **~1.95–1.99 Å** | measured (raw 1.832 + the projection's +0.16 Å) |
| a perfect *filter* (oracle on the answer) at m = 75 | 1.963 Å | measured |
| a perfect *filter* at its co-optimised m = 5 | **~1.48–1.64 Å** | measured |
| an objective at ρ = 0.7 with co-optimised m | ~1.70 Å raw | measured |
| **torsion restraints at σ = 12°, ~100% residue coverage** | **1.49 Å (86% under 2 Å; FAIL18 1.80)** | measured, ORACLE (§ XVII) |
| torsion restraints at σ = 12°, 90% coverage | 2.14 Å | measured, ORACLE |
| torsion restraints at σ = 12°, 75% coverage | 2.79 Å | measured, ORACLE |
| torsion restraints at σ = 12°, 50% coverage | 3.29 Å — *worse than the incumbent* | measured, ORACLE |
| torsion restraints at σ = 6°, ~100% coverage | 0.86 Å (99% under 2 Å) | measured, ORACLE |

**The realistic expectation for the next sprint is no accuracy improvement.** Anyone
promising < 2.0 Å from the distance channel is promising something the ceiling forbids.

---

# XIII. HIGHEST-VALUE NEXT EXPERIMENTS, RANKED

Ranked by (expected gain × orthogonality × failure-target relevance) ÷ (cost × leakage risk).

1. **BMRB chemical shifts → φ/ψ restraints. The ceiling diagnostic is DONE (§ XVII) and it
   clears the target: 1.486 Å at σ = 12° with full coverage, 86% of targets under 2.0 Å, and
   the failure class at 1.80 Å instead of 6.019.** What remains is one cheap, decisive
   measurement before any structure code is written: **what fraction of residues will a
   shift-based predictor actually predict on flexible 9–16-mers?** The channel needs ~90%
   coverage to reach 2.0 Å and is *worse than the incumbent* below 50%. That fraction can be
   established from BMRB deposits alone. If it is below ~75%, close the direction. Leakage
   discipline: shifts are an independent observable but the deposited coordinates were solved
   using them, so this is NMR-restrained prediction and belongs in its own column. Note also
   that the value is in BUILDING from the torsions, not in filtering the pool with them — a
   successful version replaces the architecture rather than extending it.
2. **The class screen as a declared domain of applicability.** Aggregation-assay-trained
   detectors carry zero structural leakage. Measure precision/recall on the 126 before any
   routing; one false positive costs +4.7 Å. Report in-scope accuracy as a *separate declared
   number*, as `monomer_manifest.json` already does for the benchmark.
3. **~~Re-price the record's eleven negative ranking results~~ — DONE this sprint, and
   S9-8's closure survives.** The terminal compresses ranking effects by 1.6–2.0× without
   flipping a sign; no discarded ranker reaches the threshold where the compression would
   have mattered. *Replaced by:* score every future ranker on the (ranker × m) surface rather
   than at m=75, since a fixed-m evaluation understates a ranking improvement by ~39%.
4. **The set decoder on the K=500 pool rather than the pre-filtered top-75.** Features are
   cached (`s12/cache/agg_feat/*_topk200.npz`); it is the only variant that could rescue
   members the filter discarded. Expect a negative.
5. **In-band selection on the *assembled* space rather than on retrieval pools.** Every
   in-band result on record was measured on pools of whole library windows; the assembled
   space has different geometry (top-40 pairwise diversity 1.90 Å) and the native sits at the
   23.6th score percentile there — the best of any space measured this sprint. The prize is
   the measured 1.265 Å oracle-25 ceiling.
6. **Deliberate pool diversification on flagged targets.** On the FAIL18, BLOSUM is *worse
   than noise* and a content-free capacity null is worth −0.353 Å. Needs no new channel and
   is nearly free — but requires the § XIII.2 detector to know where to apply it.
7. **A curated irregular-conformation fragment set.** The library's holes are at 7–8 residues
   and concentrated in coil (0.754 Å for pure-coil 8-mers vs 0.199 for helical); the worst
   8-mer hole correlates r = +0.815 with the target's whole-window floor. If data is ever
   collected, select for irregular local conformations, not for more peptides.
8. **An in-distribution subset of tuning126**, defined by experiment class from headers and
   never by performance, so in-scope accuracy can be tracked separately.
9. **Re-tune m whenever the objective moves.** Zero value today, up to 0.386 Å later; the
   derivative is now measured and belongs in any objective work's evaluation.
10. **The length-extrapolation benchmark** (§ XI.2) — the only adequately-powered fresh
    instrument available, worth building before it is needed rather than under pressure.

---

# XIV. THINGS NOT WORTH DOING

- **Any further in-band ranker over retrieval pools.** The untested class is now tested; the
  learning curve is flat against a loud oracle control.
- **Finer structural alphabets, ESM features for local conformation, key fusion, more
  fragment training data.** Resolution degrades the key monotonically; ESM adds 0.002
  accuracy; every key correlates r ≥ 0.87 with BLOSUM's errors; 8× more fragment data costs
  0.072 accuracy.
- **Post-hoc correction of the distance objective in any form.** Two independent deaths now
  (calibration, S7-3; error coherence, this sprint) with different mechanisms.
- **Growing the fragment library.** −0.432 Å per decade; 800× for 0.8 Å.
- **Ideal-geometry template banks.**
- **Generic multi-piece assembly as a generation lever**, unless a ranker with real in-band
  skill exists first.
- **Foundation models requiring a GPU**, and any peptide *design* model.
- **Manufacturing a quantum advantage.** Classical annealing solves the assembly and subset
  problems to optimality; the honest role is the benchmark, and the accuracy is zero either way.
- **Spending benchmark60.** It is gone; nothing in this sprint would justify it.

---

# XV. SCIENTIFIC NOVELTY

1. **The terminal-operator transfer law.** `d_output = 1.16·d_set_mean + 0.04·d_set_best`,
   R² 0.89, reproduced independently by two agents. A consensus operator over a candidate set
   is a near-linear map of the set's *mean* quality and is almost blind to its *best* member.
   This is a general statement about consensus-based structure prediction: it predicts that
   ranking improvements will measure null in any pipeline ending in an average, and it gives
   the design rule that optimal cardinality shrinks as the objective improves
   (m* = 500 → 75 → 20 → 3–5 across the quality range).
2. **Error coherence, not error magnitude, decides whether a learned corrector helps.** At
   identical 0.688 per-pair accuracy, a corrector whose mistakes are coherent with the
   predictor's own error structure emits +0.291 Å and one whose mistakes are i.i.d. emits
   −0.142 Å. A corrector trained on the predictor's own features inherits its error structure
   by construction, and feature pruning does not break that. A concrete, transferable warning
   about learned residual correction in structured-prediction pipelines.
3. **A corpus-level measurement of the sequence→structure channel at peptide length**, and
   the finding that protein-derived fragments carry ~7× less of it than isolated peptides —
   unifying the training-side result (S7-2) and the retrieval-side result under one mechanism.
4. **A worked demonstration that predictor-quality metrics do not price structure quality**,
   with three independent metrics failing (MAE, rank correlation, per-pair sign accuracy) and
   an iso-MAE surface showing a 0.6 Å spread in emitted accuracy at fixed MAE.
5. **A negative result about short-peptide structure prediction that the field's own record
   supports**: the isolated 9–16-residue NMR peptide is a documented blind spot of the
   AlphaFold3-like family, and this sprint measures *why* from the information side rather
   than the model side.
6. **The methodological record itself** — twenty pre-registered hypotheses, twenty negatives,
   every one with its null, and an adversarial audit that refuted two of the coordinator's own
   claims and one of its own agents'. That discipline is unusual and is why these numbers can
   be trusted.

---

# XVI. RESEARCH LEDGER

**PROVEN** — the terminal-operator transfer law; the 1.95 Å ceiling of the distance channel;
the 0.776 Å total sequence-information content of the architecture; the peptide/fragment
corpus asymmetry; the fibril/lasso enrichment; lasso unreachability by linear-window
retrieval; the instrument's fidelity to production.

**STRONGLY SUPPORTED** — in-band discrimination is signal-limited, not sample-limited;
learned correctors inherit the error structure they correct; predictor accuracy does not
price emitted structure; the terminal operator is at its native-free optimum at m = 75 today.

**PROVEN (ORACLE, prices a future modality)** — torsion restraints at σ = 12° with full
coverage emit 1.486 Å, 86% under 2.0 Å, FAIL18 1.80; the channel works by BUILDING, not
filtering; it needs ~90% residue coverage and is negative below 50%.

**OPEN** — whether a shift-based predictor achieves ~90% residue coverage on flexible
9–16-mers; in-band selection on the assembled space; the set decoder on the unfiltered K=500
pool; an aggregation-assay class detector's precision.

**WEAK** — that peptide-weighted retrieval is worth anything (direction real, size ~0.02 Å);
that deliberate diversification helps the failure class (helps recall, hurts the 108).

**REFUTED** — the 2.406 Å pool cap; the 0.53/0.32 Å two-piece assembly floor; the circular
ABEGO oracle; separation-partialled scoring; the ESM contact-agreement lead; the
residual-target distogram retrain; the SS/ABEGO retrieval key; the error-direction head;
ideal-template augmentation; every router.

**INVALIDATED BY AUDIT** — "the objective's own optimum is 4.62 Å" (→ 3.53); "r_sep = 0.196"
as an absolute (estimator-dependent → 0.37–0.39; the contrast survives); "perfect filter at
m = 25 = 1.644 Å" (that is a perfect-*RMSD* filter; a perfect *distance objective* gives
2.123); **"no objective upgrade of any kind reaches 2.0 Å"** (true only at m=75; with a
co-optimised terminal a perfect distance matrix reaches ~1.99 Å — the sprint's most
consequential correction); "the operator is almost blind to the set best" as a general law
(coefficient 0.039 is design-specific; 0.18–0.30 out of design — narrowed, not withdrawn);
"the FAIL18 return 6.0 Å" as an independent fact (near-tautological given the selection rule);
every containment-based leakage screen in this repository used at a 0.6 threshold.

**UNKNOWN** — whether any architecture can reach < 2.0 Å on this target distribution without
experimental restraints.

---

# XVII. THE ONE DIRECTION THAT REACHES THE TARGET, PRICED BEFORE IT IS BUILT

`s12/coord_restraints.py`, `s12/coord_restraints2.py` → `s12/results/coord_restraints*.json`.
**ORACLE / DIAGNOSTIC throughout.** These arms read native torsions to synthesise restraints;
they price a *future modality* and are not a method.

§ XIII ranks BMRB chemical shifts → TALOS-N-style φ/ψ restraints as the top next experiment.
Building it is days of work — NMR-STAR parsing, residue mapping, a torsion predictor — so it
was priced first, without touching the BMRB, by corrupting native torsions to the accuracy a
shift-based predictor actually delivers and running the result through the pipeline.

### The headline: complete torsions at TALOS-N accuracy reach 1.49 Å

| arm at σ = 12°, all residues restrained | mean | FAIL18 | other-108 | < 2 Å |
|---|---|---|---|---|
| shipped synthesis (baseline) | 3.213 | 6.019 | 2.745 | 0.286 |
| **build the chain directly from the restrained torsions** | **1.486** | **1.80** | 1.43 | **0.86** |
| use the restraints as a pool FILTER, then the normal synthesis | 2.972 | 4.78 | 2.67 | 0.31 |
| fuse the restraints with the distogram score as a filter | 2.890 | 5.24 | 2.50 | 0.30 |

Paired: **−1.726 Å [−2.016, −1.464]** against the incumbent. This is the only intervention
measured anywhere in the sprint that goes below 2.0 Å, and it is the only one that fixes the
failure class rather than working around it — 6.019 → 1.80 on the FAIL18.

**Two architectural statements follow immediately.**

1. **The value is in BUILDING, not in FILTERING.** As a retrieval/filter channel the same
   restraints are worth only −0.24 to −0.32 Å. With real local-conformation information the
   retrieval pool is not needed; without it, the pool is all there is. A shift-restrained
   system is therefore a *different architecture*, not an extra channel bolted onto this one.
2. **It confirms the sprint's diagnosis from the opposite direction.** The whole pipeline is
   an elaborate attempt to infer local backbone conformation from sequence. Supplied directly,
   that information is worth 1.7 Å. Predicted from sequence, it is worth nothing — the
   retrieval-key agent's leave-fold-out torsion-bin predictor reached 0.690 accuracy overall
   and 0.517 on the FAIL18, and emitted +0.015 Å.

### The make-or-break parameter is COVERAGE, not accuracy

The first pass filled unrestrained residues with zeros, which is not "no information" but a
specific absurd conformation, and torsion error compounds along a chain. `coord_restraints2`
fills them honestly instead — from the pool member whose torsions best match the *restrained*
residues (`fill_pool`), or from the pool's own Ramachandran mode (`fill_rama`).

| σ | coverage | fill_pool | fill_rama | fill_zero | paired vs base | < 2 Å | FAIL18 |
|---|---|---|---|---|---|---|---|
| 6° | 1.00 | **0.861** | 0.861 | 0.861 | −2.352 [−2.656, −2.080] | 0.99 | 0.97 |
| 6° | 0.90 | 1.672 | 1.983 | 2.562 | −1.540 [−1.801, −1.303] | 0.68 | 2.38 |
| 6° | 0.75 | 2.488 | 2.734 | 3.596 | −0.724 [−1.012, −0.443] | 0.45 | 3.28 |
| **12°** | **1.00** | **1.486** | 1.486 | 1.486 | **−1.726 [−2.016, −1.464]** | **0.86** | **1.80** |
| **12°** | **0.90** | **2.141** | 2.401 | 2.827 | **−1.072 [−1.342, −0.821]** | 0.56 | 2.63 |
| 12° | 0.75 | 2.792 | 2.948 | 3.684 | −0.421 [−0.674, −0.170] | 0.37 | 4.04 |
| 12° | 0.50 | 3.290 | 3.496 | 4.404 | **+0.078 [−0.147, +0.298]** | 0.25 | 4.81 |
| 20° | 1.00 | 2.408 | 2.408 | 2.408 | −0.804 [−1.105, −0.521] | 0.36 | 2.58 |
| 20° | 0.75 | 3.258 | 3.381 | 4.055 | +0.045 [−0.237, +0.309] | 0.15 | 3.95 |

Pool-based gap-filling beats the Ramachandran mode by 0.16–0.31 Å and the zero fallback by
0.9 Å, so the retrieval library does have a real job in a hybrid system — supplying the gaps.
But it cannot rescue sparse coverage: **the channel is worth nothing below ~50% coverage and
turns negative there**, and it only reaches 2.0 Å at roughly **≥ 90% coverage and ≤ 12°**.

**What this means for the build decision.** The experiment is worth doing, and the parameter
that decides it is not the torsion accuracy TALOS-N quotes but **the fraction of residues it
will predict at all on flexible 9–16-mers**. That is the first thing to measure, it can be
measured from BMRB deposits alone without building any structure pipeline, and if the answer
is below ~75% the direction should be closed before anything is built. The leakage caveat
travels with all of it: shifts are an independent observable, but the deposited coordinates
were determined using those same shifts plus NOEs, so any such arm is **NMR-restrained
prediction reported in its own column**, never an improvement to the sequence-only system.

---

# XVIII. THE MANDATED QUANTUM COMPONENT — a verified problem, and a classical win

`s12/vq_FINDINGS.md`. The brief requires a genuine VQE and a genuine CVaR and forbids forcing
them onto a task classical optimisation does better. Both halves were honoured.

**A real 2-local Hamiltonian was found and verified.** The target is cut into k contiguous
segments, each choosing one of m library pieces, each Kabsch-placed onto the shipped `fit_ca`
(a deployable anchor — no native). Because a placed piece depends only on its own segment,
every residue pair falls in exactly one block and the objective decomposes exactly:

    npairs · F(x) = Σ_s h_s(x_s) + Σ_{s<t} J_st(x_s, x_t)

Verified against a from-scratch rebuild at **max |Δ| = 0.0 in float64**, one-hot QUBO ↔ tables
at 1e-15, Ising ↔ QUBO at 1e-13, fully connected. This is a genuine combinatorial problem,
unlike the shipped placement whose energy is a per-candidate score and therefore has no
interaction structure at all. Exact functional ANOVA puts the interaction at **2.4% of
Var(F)** — real, and small.

**Classical optimisation wins decisively, at matched budget:**

| instance | classical | CVaR-VQE, same budget | paired |
|---|---|---|---|
| 12 qubits, 4,096 configs | greedy+1-opt **97.6% optimal, 144 ms** | 61.1% optimal, 9.2 s | +0.0140 [+0.0093, +0.0192], **0W/63L** |
| 16 qubits, 65,536 configs | annealing **100% optimal, 0.96 s** | 32.5% optimal, 51 s | +0.0198 [+0.0141, +0.0263], **0W/85L** |
| one-hot, 0.098% feasible | annealing 75% optimal | **0% optimal** | +2.267, **0W/40L** |

The VQE also loses to uniform random sampling on 4 of 5 budgets, because its state entropy
stays at 9.1 of 12 bits at *every* budget — it never concentrates. This reproduces
`core/quantum.py`'s own recorded finding on two new problems.

**Does CVaR contribute anything a classical optimiser does not?** One thing, and it is not
worth having. Collapse-prevention is confirmed with a convergence control: at T = 0 the α = 1
arm reaches 0.000 bits of entropy against 6.1 bits at α = 0.05, stable under 12× the
iterations and matching the analytic ceiling. But the resulting breadth is **indifference, not
preference** — the p-weighted readout loses to a matched-entropy Boltzmann distribution on
**20 of 20** α×T arms, and the VQE's eight most-probable states never beat the eight
lowest-energy states on band coverage.

**And solving the problem perfectly makes the answer worse.** End to end through the real
projection: anchor 3.204, VQE 3.350, greedy+1-opt 3.369, **the exact certified optimum
3.373 (+0.169 [+0.103, +0.236])**. Every CI excludes zero and all five folds agree. The
optimiser is not what is wrong — the objective is anti-correlated with truth in this family.

**Recommendation, adopted into § IX:** keep the component, re-site it from the top-128
selector onto this verified assembly Hamiltonian, and publish the classical comparison
permanently beside it. Keep CVaR for collapse-prevention, which is real and analytically
explained, and drop the claim that its distribution is a better multi-hypothesis output —
replace that readout with `softmax(−E/T)`.

**One incidental positive, and it is generation-side, not quantum.** Against a shuffled-risk
null, the shipped objective genuinely ranks *inside an anchored family*: correlation +0.449
real vs +0.205 shuffled, optimum at the 24.8th RMSD percentile vs 53.3rd, paired −0.678 Å
[−0.835, −0.531], 104W/22L, surviving rg-partialling (+0.333 real vs −0.017 shuffled). The
record's "the objective does not rank" was measured on *unanchored* candidate sets. On the
FAIL18 it collapses back to chance (+0.119, 45th percentile).

---

# XIX. THE ADVERSARIAL AUDIT'S FINAL LEDGER

`s12/adv_FINDINGS.md` (1,011 lines). The instrument was **cleared** before any output was
believed: `shipped_score` reproduces production's top-75 on **126/126** targets under three
precisions, and `coordinate_average` matches production `avg_ca` to 2.1e-07 Å.

**Refuted:** the 4.62 Å objective optimum (→ 3.532 Å by direct L-BFGS; the conclusion survives
at a quarter the margin, and the attack *strengthens* it — the native scores worst of all
arms on 117/126); assembly as evidence that the objective cannot rank (a capacity artefact —
2.3×10⁵ i.i.d. Ramachandran chains beat the real bank by +0.216 Å); and **"no objective
upgrade reaches 2.0 Å"**, which holds only at frozen m = 75.

**Weakened:** r_sep (0.196 → 0.367–0.391 under pooled leave-fold-out partialling; the contrast
survives every estimator); the FAIL18's fibril/lasso enrichment under an *absolute* rather
than relative band definition (p 0.097/0.244 — the relative-band version survives Bonferroni
over all 66 post-hoc hypotheses at 8.1e-05 and 4,000 max-stat permutations at 0/4000); and my
own set-best coefficient (design-specific).

**Survived:** the fibril/lasso enrichment as defined, with zero label disagreements on an
independent header re-read; the no-fresh-benchmark finding, with all 67 length rejections and
all 16 CA-step rejections individually verified (best defensible relaxation yields ~18, not
40); the aggregation agent's operator optimum, extended along cardinality — **leave-fold-out
picks m = 75 on all five folds, d = 0.0000**; and S9-8's closure.

**The largest single effect the sprint found, and it is a routing effect:** on the FAIL18,
**discarding the distogram entirely** (m = 500, i.e. averaging the whole pool) is worth
**−0.510 Å [−0.857, −0.177], 15W/3L**, against **+0.631 Å [+0.319, +0.983]** on a length- and
fold-matched control — a **1.14 Å interaction** with both CIs excluding zero, confirmed
through the real projection. But a perfect binary router caps at −0.073 Å on the full
instrument, because the failures are only 14% of targets, and per-target oracle cardinality
caps at −0.313 Å. So it is a large effect on a small subgroup with no deployable trigger.

**Defects found in my own instrument and code, all recorded and two fixed:** the λ = 0.3
projection arm carries a ~0.02 Å reproducibility floor (documented; the λ = 0 arm reproduces
to 7e-4 Å and is the one to quote below that threshold); `I.write` had no config key so a
partial run could silently overwrite a complete one (**fixed** — it now records `n_rows`,
`n_expected` and a `complete` flag); and two assembly scripts use an unstable argsort on tied
scores, which is the exact pattern that produced the record's phantom 1.386 Å winner.
