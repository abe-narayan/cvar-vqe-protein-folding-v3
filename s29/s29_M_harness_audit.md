# S29 EVALUATION-HARNESS AUDIT (lane M, 2026-09-20 00:01)

The charter's section 9 lists "evaluation-harness failure (yes, check this too)" among the failure
modes S29 must be able to distinguish. This is that check. Nine checks, every one with its command
and its number; anything that fails is reported as a failure.

**Script:** `s29/s29_M_harness.py` (nine checks, one process, 61 s, peak RSS 0.549 GB).
**Command:** `python s26/jobrun.py --agent S29M --tag CPU --name m_harness_audit3 --est-ram 1.5 --
python s29/s29_M_harness.py`
**Log:** `s26/logs/m_harness_audit3.log` (exit 0, wall 65.2 s).
**Artefact:** `s29/results/s29_M_harness_audit.json`, provenance-stamped by `ST.save_atomic`
(source sha256 `4f10460fce236435`, git `d4b52a6c1ffc`, started 2026-09-20T07:01:07Z).

## VERDICT

**All nine checks pass, and on that evidence the harness is sound**, with one declared,
fully-quantified non-bit-exactness which is stated here rather than in a footnote and which is
proven not to move a number:

> **The s12 distogram cache is NOT bit-identical to a fresh recomputation.** Over all 126 targets,
> max |prob_fresh - prob_cached| = 1.97e-06, max |risk diff| = 1.34e-04 (relative 5.63e-06), max
> |score diff| = 1.91e-06 against a per-target score sd of ~1.25. A fresh recomputation is exactly
> deterministic **in-process** (0.0e+00 between two calls), so the difference was introduced when
> the cache was written -- a different BLAS/thread reduction order or library version, not
> nondeterminism now. **Consequence, measured: the full 500-candidate score order differs on 2 of
> 126 targets (deep in the tail), the top-75 SET differs on 0 of 126, and the emitted point-cloud
> RMSD differs by 0.00e+00 A on 126 of 126.** No result in the record is affected; a future arm that
> depends on the *order* below the cut (not on the set) must recompute rather than read the cache.

---

## Check 1 -- `ca_rmsd` is the Kabsch CA RMSD. **PASS**

An independent five-line implementation written from the definition (centre, SVD of `a.T @ b`,
reflection fixed by `sign(det)`, `sqrt(mean squared deviation)`), compared against
`s12/instrument.py:91 ca_rmsd` on 10 targets x 6 structure pairs each (pool members 0, 1, 7, 100,
499 against the native, plus the top-75 coordinate average against the native), with a rigid-motion
invariance check (a random rotation + translation of one argument).

    worst |I.ca_rmsd - independent Kabsch|, 10 targets x 6 pairs + invariance:  3.675e-14

Sample (1A13 / 1A1P / 1CEK): `ca_rmsd` 3.6077311 / 3.3158968 / 0.5945119, independent
3.6077311 / 3.3158968 / 0.5945119. **Reflections are correctly forbidden**: a z-mirrored copy of the
same window scores differently (1A13 4.6676 vs 3.6077; 1CEK 3.0516 vs 0.5945), and the independent
implementation agrees on the mirrored value too -- so the determinant correction is present in both
and the metric is chirality-sensitive, which is what makes lowest-objective multi-start selection
safe in the projection (`core/project.py` module docstring).

## Check 2 -- the native used is the pinned one. **PASS** (with a correction to the brief)

    universe nat_ca vs peptide_db.npz ca, all 126 targets:        0.000e+00 A
    universe nat_ca vs a FRESH parse of the deposited PDB, 10/10: 0.000e+00 A
    targets absent from peptide_db or with a disagreeing sequence: none

Every target's native is `s8/generate_univ/<pdb>.npz :: nat_ca` (`s12/instrument.py:53 load_univ`,
read as float64 from float32 storage), which is bit-equal to `core.data.load()`'s `ca` and to a fresh
`core/geometry.py:823 native_coords_from_pdb` parse of the deposit (model 1). The exact zeros are
expected, not suspicious: the parser already stores float32-representable values, so the universe's
float32 round trip is exact.

**Correction to the brief.** `catrace_prior.npz` is **not** a native store. It holds a single key,
`nlp`, a (24, 36) float64 pseudo-angle log-density table (`core/geometry.py:884 CATRACE_CACHE`,
`NTHETA = 24`, `NTAU = 36`) -- a CA-trace geometry prior. The pinned native is the universe's
`nat_ca`. Any S29 lane looking for "the pinned native" should read the universe.

## Check 3 -- no benchmark file is read anywhere in the tuning path. **PASS**

Two parts.

**(a) Static grep** over the whole import closure of `s27/run_vqe_chain.py` -- `s27/run_pool.py`,
`s27/ham_lib.py`, `s25/phys_lib.py`, `s24/d_harness.py`, `s24/stats_lib.py`, `s22/qcand_lib.py`,
`s16/energy_lib.py`, `s15/seed.py`, `s12/instrument.py`, `core/{quantum,project,predict,geometry,
cache,data}.py` -- for `benchmark_manifest|final_report\.json|monomer_manifest|\.benchmark\(\)`:

    hits: core/data.py only -- lines 588, 589 (the two path CONSTANTS) and 646 (a docstring).

No call site on the path. The two functions that would read the manifest, `core/data.py:644
benchmark` and `682 monomer_benchmark`, are reached only from `core/pipeline.py:399`
(`targets_of("benchmark60")`, which the CLI refuses without `--i-am-spending-the-benchmark`) and
from `s7/debias.py:103 tuning_targets`, which used the manifest **once, at definition time**, to
make the 126 cluster-disjoint from the sealed 60. The run-time path never calls it: the 126 come
from `s12/instrument.py:41 targets()`, a glob over the pinned universes.

**(b) Poison test, and it is stronger than the brief's rename.** `BenchmarkPoison`
(`s29/s29_M_harness.py:52`) replaces `builtins.open`, `os.path.exists` and `np.load` in-process so
that **any** touch of `benchmark_manifest.json`, `final_report.json` or `monomer_manifest.json`
raises. With the poison live, the deployed path was run end to end on 6 targets (pool construction
from the universe, the shipped score, the tie-safe top-75, the uniform average, the production
projection, the RMSD) and, separately (check 8), on all 126 targets through a cold distogram:

    forbidden reads with the poison live, 6 fresh rebuilds + 126 cold distograms:  0

Why not the brief's rename: a rename is defeated by a value already in memory or by a cached read,
and copying the tree would copy 1.1 GB of pinned caches; more importantly, **contract rule 2 forbids
moving a pinned artefact**, and `results/benchmark_manifest.json` is one. The interception is
strictly stronger and touches nothing on disk.

## Check 4 -- `pinned_folds` are the pinned folds, and the leave-fold-out models exclude. **PASS**

    ST.pinned_folds()            == the universes' `fold` field, all 126:      True
    the universes' `fold`        == peptide_folds.json[seq], all 126:          True
    ST.pinned_folds([10 pdbs])   == the same 10 entries, by pdb:               True
    10 targets (the first 5 and last 5) served by fold_model(own fold),
      checkpoint present on disk, and their own sequence ABSENT from that
      fold's training entries:                                                 True
    verbatim fragment self-copies among those 10:                              none

`s24/stats_lib.py:65 pinned_folds` reads the instrument (`I.targets()`), never `data.folds()`, which
is the right thing (memory `benchmark-and-folds-must-be-pinned`). The structural half of this check
says the right file is served; check 9 is the empirical half, and it is the one that could have
caught a mis-trained checkpoint.

## Check 5 -- the production anchor 3.2126 reproduces. **PASS**

From the stored rows (`s27/results/chain_rows.jsonl`, config `DIS`, n = 126):

    mean rmsd_chain  3.2126   (anchor 3.2126)
    mean rmsd_cloud  3.0483   (anchor 3.0483)

From a **fresh re-projection** on the first 6 targets (1A13, 1A1P, 1CB3, 1CEK, 1CS9, 1D0W), rebuilt
from `s27/cache/<pdb>.npz` through `RP.zr(DIS)` -> `np.lexsort((tiekey, E))[:75]` ->
`H.readout_uniform` -> `H.readout_projected` -> `I.ca_rmsd`, with the poison of check 3 live:

    worst | fresh - stored row | over 6 targets, cloud AND chain:  0.000e+00 A   (18 s)

The anchor is reproducible bit-for-bit, and the path that reproduces it reads no benchmark artefact.
**Note for every lane:** this path contains no VQE (`s27/run_vqe_chain.py:124-130`), which is a
property of the anchor, not of this audit -- see `s29/DATAPATH.md`, "Read this first".

## Check 6 -- `ST.compare`'s MDE and its fold CI. **PASS**

Source read: `s24/stats_lib.py:56` carries the literal `MDE_K = 2.8016` (pinned so every lane's MDE
is bit-identical), `s24/stats_lib.py:119-127` builds the fold CI by resampling **FOLDS** with
replacement (`rf.choice(F, len(F), replace=True)`, then concatenating those folds' per-target
differences) -- a genuine cluster bootstrap over the 5 clusters, 4,000 draws, beside a separate iid
bootstrap. Reproduced on a synthetic paired sample (n = 126, 5 folds):

    ST.compare SE   0.026133   recomputed sd(d)/sqrt(n)  0.026133   (identical)
    ST.compare MDE  0.073215   2.8016 x SE               0.073215   (identical)
    n_folds 5, ci95_fold present beside ci95_iid

`_verdict` (line 145) requires **both** |effect| > its own MDE **and** the fold CI to exclude zero,
and refuses to return BETTER/WORSE at all when `folds` is omitted -- both guards were added after the
audit that found a 0.39x-MDE effect labelled "MEASURED". The known limitation is in the code and
stands: with only 5 clusters the fold CI is unstable.

## Check 7 -- the 126 are the 126. **PASS**

    targets 126 | distinct pdb ids 126 | distinct sequences 126 | lengths 9..16
    universe files in s8/generate_univ 126 | order is pdb-sorted: True
    length histogram  9:9  10:11  11:13  12:19  13:23  14:14  15:16  16:21

## Check 8 -- the cached posterior is not stale in any way that moves a number. **PASS, with the declared non-bit-exactness**

A cold recomputation of the distogram through the MLP for **all 126 targets** (`pl.fold_model` ->
`dgm.Distogram.for_target`), with the poison live, compared against the `s12/cache/disto_<pdb>.npz`
the deployed score actually reads:

    forbidden reads:                                   0
    in-process determinism (two calls, one process):   0.0e+00        <- exact
    max |prob_fresh - prob_cached|, 126 targets:       1.97e-06
    max |risk_fresh - risk_cached|:                    1.34e-04   (relative 5.63e-06)
    max |score_fresh - score_cached| over 500 cands:   1.91e-06   (score sd ~1.25)
    500-candidate score ORDER differs on:              2 / 126
    top-75 SET differs on:                             0 / 126
    max |point-cloud RMSD difference|:                 0.00e+00 A on 126 / 126

The pass criterion is the deployed consequence (no set, no RMSD moves), which is the only way a
stale cache could change a result. The bit-equality sub-question answers NO and is reported as such
in the verdict above. **Actionable consequence:** an S29 arm whose operator depends on the score's
order *below the top-75 cut* -- a ranking correlation over the full 500, a reranker's tail, an
order-statistic null -- should recompute the posterior rather than read the cache, or accept a
1-in-63 chance that its order differs from another lane's.

## Check 9 -- EMPIRICAL leave-fold-out, and it is the strongest evidence here. **PASS** (ORACLE DIAGNOSTIC)

A structural check (check 4) proves only that the right file is served. This one asks whether the
checkpoint *behaves* like a model that never saw its fold. For every fold model f and every target t
it computes the mean negative log-likelihood the model assigns to the target's **true** distance bin
(ORACLE, a diagnostic; it tunes nothing and chooses nothing), then compares, per target, the
**deployed** model (the one for that target's own fold) against the mean of the four models that
**did** train on that fold:

    mean NLL, DEPLOYED (own-fold, held out):                 3.28777
    mean NLL, the four models that SAW this fold:            1.21275
    delta (deployed - seen):                                +2.07502   SE 0.15166   4.88x MDE
    per-fold deployed NLL:  0: 3.0636  1: 3.1248  2: 3.6206  3: 3.4722  4: 3.1808

The sign is correct on the aggregate and on every fold: **the model that ships is decisively worse
on its own targets than the models that trained on them.** A checkpoint trained with the wrong
exclusion would show the opposite sign. Leave-fold-out is real, on all 126 targets, empirically.

**And a finding that falls out of it, for the other lanes.** The gap is 2.075 nats -- the models that
saw a peptide assign its true bins about **e^2.075 = 8.0x** more probability than the held-out model
does. The distogram MLP **memorises its training peptides heavily** (372,881 parameters against 787
peptides plus 6,003 fragments, no dropout by design, `core/predict.py:335`). Two consequences:
(i) any in-sample diagnostic computed on the corpus -- calibration, sharpness, MAE, a fitted
correction -- describes memorisation, not the model the pipeline deploys, and must be computed
out-of-fold; (ii) the "sequence signal" the record prices (memory `sequence-signal-is-the-ceiling`)
is the 3.288-nat held-out model, not the 1.213-nat in-sample one, and the 8x gap is an upper bound
on what more capacity of this kind could buy if the generalisation gap could be closed.

---

## What this audit did NOT check, stated so nobody over-reads it

- **The VQE arm's own harness.** `arm_vqe`, `exact_face`, `gate_set_equality` and the set-equality
  theorem were verified in S25 (`s25/q_verify.py`, 2,592 adversarial cells; `s25/QUANTUM.md`
  section 5) and re-run by S24 and the coordinator; this audit did not re-derive them. It did
  establish (check 5, `s29/DATAPATH.md`) that **the 3.2126 anchor does not pass through them**.
- **The AMBER stage** (`core/amber.py`), not on the `chain_rows` basis.
- **`s12/instrument.py:266 selfcheck`** was not run here: it asserts `pool_best 1.7108` and
  `top75_best 2.3062` and re-derives FAIL18, and it reads the production cache
  (`bench_results/cache/1fc9f2dcf489e2fb`). A lane that changes the pool must run it.
- **The projection's own degeneracy.** `core/project.py`'s docstring records that the reference
  disagrees with itself by up to 1.6 A under a rigid motion of its input on some targets. That is a
  property of the stage, it is measured and declared, and this audit reproduced the *emitted* values
  bit-for-bit (check 5) rather than re-opening it.
- **Only 6 targets were freshly re-projected** (the brief's number). The projection is the expensive
  stage; 126 would be ~20 min. The cold-distogram half of the path was checked on all 126.

## Corrections to the brief, for the coordinator

1. `s25/qcand_lib.py` does not exist. The `Encoding` class is **`s22/qcand_lib.py:120`** (and there
   is an older one at `s14/vqe_encoding.py:50`). `s24/d_harness.py:113` imports the s22 one.
2. `catrace_prior.npz` is a pseudo-angle prior (`nlp`, 24 x 36), not a native store. See check 2.
3. The brief assigns the harness audit to lane M (`s29/briefs/S29M.md` deliverable 3) while
   `s29/STATE.md` ("Standing measurements before any build") assigns it to lane D. This file is M's;
   if D also ran one, the two should be reconciled rather than duplicated.
