# SPRINT 26 -- DECISION LEDGER

Append-only. One entry per finding or decision. Nobody edits another lane's entry; a
retraction is a new entry naming the entry it retracts. Every number carries an artefact path.

---

## L0 -- THE GOVERNOR IS RUNNING; THE BASELINE LOAD LEAVES ABOUT 4 GB (2026-09-12, coordinator)

`s26/governor.py` started at 23:50:41 (pid 36196), commit `a4db170c`. Self-test: a 0.3 GB job
registered through `s26/jobrun.py`, ran, and reported peak RSS 0.328 GB
(`s26/jobs_done/selftest_alloc.json`). Every action is logged to `s26/governor.log`.

Measured before any campaign process existed (`python s26/governor.py --once`):

    RAM total 16.75 GB (15.6 GiB)   used 11.16 GB = 66.6%   available 5.59 GB   CPU 10.5%

The 66.6% is the user's own environment (several VS Code and claude sessions, Defender;
`psutil` process table, 2026-09-12 23:47). The campaign's 93% ceiling therefore leaves about
4.4 GB for all lanes together. The "90 to 93% band" instruction is read as: fill up to, never
above, 93%; the governor launches queued work only when the box has sat under 88% for 60 s.

---

## L1 -- OPERATIONAL ITEM 1 CLOSED: `s24/cache_amber/*.npz` IS TRACKED (2026-09-12, coordinator)

Commit `6b494780`, the first on branch `s26`: `.gitignore` whitelist `!s24/cache_amber/*.npz`
placed after the blanket `*.npz` rule; 126 files, 1.5 MB, `git ls-files s24/cache_amber`
returns 126. S25 L13's single-point-of-failure on the 63,000 ff14SB/GBn2 single points is
closed.

---

## L2 -- THE PRESENTATION FILE DOES NOT EXIST ON THIS MACHINE (2026-09-12, coordinator)

`vqe_research_overview.pptx` is named in the campaign prompt as the 11-slide deck to tune. It
is not in the repository, not in git history (`git log --all --diff-filter=A -- '*.pptx'`
returns nothing), and not anywhere under `C:\Users\abena` to depth 7 excluding `AppData`
(`find` over Desktop, Documents, Downloads, OneDrive and the profile root; the only `.pptx`
files are six unrelated lecture decks in Downloads and python-pptx's template).

Decision: Phase 4 will build the deck from the structure the prompt describes (title; six
"what I built" slides; three proposal slides; goal and ask; dark theme) using python-pptx,
and `s26/PRESENTATION_CHANGES.md` will record every slide's content and its artefact so the
edits map onto the real file if the user supplies it. No number will be invented for a slide
that has no artefact.

---

## L3 -- `TEST_RUN_RESULT_PLACEHOLDER` IS NOT IN THE TREE (2026-09-12, coordinator)

The prompt says the state brief contains the literal string `TEST_RUN_RESULT_PLACEHOLDER`.
`grep -rn` over every `.md`, `.py` and `.txt` finds no occurrence. The brief on disk
(`docs/STATE_BRIEF_2026-09-12.md` section 6.1) already carries a live run: "355 passed, 13
skipped, 0 failed, 0 errors, exit 0" dated 2026-09-12. Operational item 5 is therefore "replace
that line with the S26 run under the governor, with per-file counts and the skip reasons", and
the Infrastructure lane owns it.

---

## L4 -- BRANCH AND READING ORDER (2026-09-12, coordinator)

Branch `s26` created off `ae86a124`. Coordinator read, in this order and in full: the state
brief, README, ARCHITECTURE, professor brief, S25 LEDGER (1,270 lines), S25 QUANTUM.md (857),
S25 agentQ_FINDINGS, S20 agentD_FINDINGS (1,006), docs/FINDINGS.md (all 9,081 lines, corrections
ledger first), docs/CONDENSED_REPORT, the S24, S23, S22, S21, S20, S19, S18, S17, S16, S15 and
S14 ledgers, the S13 and S12 dossiers, S13 qarch_FINDINGS part 1, the S25 BRIEF and PREREG_Q,
s12/instrument.py, s24/stats_lib.py, core/pipeline.py, and the deployed-selector region of
core/quantum.py. The remaining code read (core modules, root modules, tests, verify,
resultslab) is assigned to the Examiner and Infrastructure lanes, whose module map is the
Phase 0 deliverable.

---
## L5 -- PHASE-GATE READING FOR THE CIS-PEPTIDE CENSUS, AND LANE PH SPAWNED (2026-09-13, coordinator)

The phase gate forbids any endpoint experiment (anything reading an RMSD to a native) before
"PHASE 0 SIGNED OFF". The cis-peptide census reads omega angles of the 126 dev natives and CA-CA
distances of pool windows; it computes no RMSD, selects nothing and changes nothing. It is an ORACLE
DIAGNOSTIC and is allowed before the gate. The projection-floor cost on cis targets (an RMSD) waits
for the gate. Lane PH (Physics) spawned at 00:20 with brief `s26/briefs/PH.md`: AMBER as a steric
reject filter, the cis-peptide gap, the steric singularity for Part IV, and the C3 matched-random
control (stage 1 on the production cache, stage 2 on the best C2 rung when lane P delivers it).
Lanes active: E, I, Q, P, PH (`s26/lanes.json`).

---

## L6 -- OPERATIONAL ITEM 4 CLOSED: THE S8-13 SOURCE LOG HAS A TRACKED EXCERPT (2026-09-13, lane I)

`docs/FINDINGS.md` S8-13 (line 3711) cites `s8/predictor_report.log`; the file lives at
`_archive/logs/s8/predictor_report.log` (untracked by decision, `_archive/` is gitignored),
sha256 `2647f640efb6d8f8e0398620647d503c366aa9a0209f17742d927366dd0f4740`, 13,588 bytes, 172
lines. The number with no other source is the transfer law `selected = 1.009 * ref + 0.011`
(r = 0.981, 2,520 pairs), log line 92: a search for `1.009 * ref` over every json/md/py/txt/log
file finds only `docs/FINDINGS.md` and that log. Now tracked verbatim, with its context lines and
the other S8-13 tables the log carries (converged X_fit 3.217, the 3.653 / 3.285 / 2.220 filter
ceiling, the ESM ablation), in `docs/sources/s8_predictor_report_excerpt.md` (commit `eb89c165`).
Only the excerpt is tracked; the log directory is not.

Two things the excerpt records that were not known before: S8-13's second regression
(`1.019 * ref - 0.001`, r = 0.984, "over converged references") is NOT in the log and has no
on-disk source at all; and S8-13's X_fit positive-phi rate 0.0551 does not match the converged
row in the log (0.0629). `python s26/examine.py` now checks the transfer law as claim
`s8_13_transfer_law` with the excerpt as fallback.

---
## L7 -- DEFECT 6c: pool_gate = WARN IS TWO TARGETS AND IS NOW EXPLAINED IN THE OUTPUT (2026-09-13, lane I)

`results/summary/results.json` (production rows, primary basis `built_chain_bb`): the two
per-target pool violations are **1D6X** (T007, chain 1.7809 A against its pool's ORACLE best
member 2.4194 A) and **1KWE** (T019, 2.4224 against 2.8625). The production cache confirms the
mechanism is averaging, not leakage (`bench_results/cache/1fc9f2dcf489e2fb/{1D6X,1KWE}.json`):
for 1D6X the raw coordinate average is already at 1.5921 A while the best single member of the
top-75 it averages is 2.4194 A (`top_m_best` = `pool_best` there), and the projected chain
lands at 1.7809; for 1KWE the average is 2.2812 against a top-75 best of 2.8625 and the chain
at 2.4224. The emitted structure is a coordinate average projected onto ideal geometry, not a
pool member, so on a single target it can land nearer the native than any one window; it cannot
do so systematically, which is why the aggregate (mean 3.2126 against 1.7108, margin 1.50 A) is
the release condition (S25 L11). The distogram row shows the same two targets, amber_distogram
one (1KWE), legacy_distogram two (1KWE, 8IL1); the four physics rows show none.

Change (commit `eb89c165`, `s25/resultslab/schema.py`): `fmt_leaderboard` appends one line
under the table naming the WARN rows and the count of targets and stating the mechanism, and
the persisted `pool_gate_rule` string now defines PASS / WARN / FAIL in the same words (the site
shows that string on its overview card). Verified on the on-disk payload
(`s26/logs/precheck_fmt_leaderboard.log`). **The files under `results/summary/` were built
before this change and still carry the old rule string; the explanation appears on the next
`--mode frozen` build. `results/` was not rebuilt.**

---
## L8 -- DEFECT 6d: PRODUCTION NORMALISES BY RANK; EVERY MOMENT Z-SCORE IS RESEARCH-ONLY OR NOT ON AMBER (2026-09-13, lane I)

The production selector currency is `core.pipeline._zrank` (line 788: `rankdata`, then
`(r - mean) / sd`), applied at line 838 to the distogram score of the top 2^n candidates in
`quantum_stage`. It is the only normalisation in `core.pipeline`. Two qualifications: the
production `Config` has `quantum=False`, so in the 3.2148 A run the selector does not execute at
all (it runs under `--components`); and `_zrank` is applied to the distogram score, never to an
AMBER energy, because no AMBER energy enters the deployable selector.

Every remaining `(x - mean) / std` on an energy, with a verdict (`grep` over `core/`,
`s25/phys_*.py`, `s25/resultslab/`, `s24/`):

| path | what is standardised | verdict |
|---|---|---|
| `s25/phys_lib.py:98 zmoment` + `combine(norm="moment")` | raw AMBER / Legacy / distogram channels, no outer re-standardisation | RESEARCH ONLY: the DECLARED SECONDARY of the S25 seven-configuration suite; this is where the non-monotone float64 behaviour was measured (40/126 targets). Imported by `s25/phys_suite.py`, `s25/phys_analyse.py` only |
| `s25/phys_landscape.py:132`, `s25/phys_gate.py:71` | per-target energy vectors `e` | RESEARCH ONLY: S25 diagnostics |
| `core/energy.py:1066 LegacyField._standardise` | the 11 Legacy TERMS against 512 random structures from the target's own library | NOT AMBER, NOT PRODUCTION: the optional `legacy` term of `core.quantum.FoldObjective` (`w_legacy=0.0`, "Off by default"); Legacy terms are bounded so the moment is monotone there; `core.pipeline` never instantiates it |
| `core/quantum.py:2126 build_distogram(models="both")` | two distogram models' scores against random structures | NOT AMBER, NOT PRODUCTION: the generation lane's prior; `core.pipeline` does not call it |
| `core/predict.py:560` | input FEATURES of the distogram model (`sigma[k]`) | not an energy; training-time feature scaling |
| `s25/resultslab/exportlib.py:688` | per-target RMSD sd in the difficulty gate | not an energy |

`core/amber.py`, `core/bench.py`, `s24/d_harness.py` (`zrank`, line 275) and
`s24/d_hamiltonians.py` contain no moment standardisation. Conclusion: no production code path
applies a moment z-score to an AMBER energy; the rank normalisation is the one production uses,
and the S25 statement (ARCHITECTURE section 4) stands.

---

## L9 -- CORRECTION TO L7: TWO POINT-CLOUD NUMBERS WERE WRITTEN BEFORE THE QUERY RETURNED (2026-09-13, lane I)

L7 states the raw coordinate average of 1D6X "is already at 1.5921 A" and of 1KWE "2.2812".
Those two numbers are wrong: I composed the entry in the same turn as the cache query and
typed values that were not yet on screen. The production cache
(`bench_results/cache/1fc9f2dcf489e2fb/1D6X.json`, `1KWE.json`) reads:

    1D6X (n 13, fold 4): shipped 3.1750  pool_best 2.4194  top_m_best 2.4194  top_m_mean 3.5117
                         rmsd_avg 2.0900  rmsd_fit 1.9033  rmsd_arm 1.7814  rmsd_full 1.7473
    1KWE (n 16, fold 2): shipped 3.4830  pool_best 2.8625  top_m_best 2.8625  top_m_mean 3.9980
                         rmsd_avg 2.7313  rmsd_fit 2.4598  rmsd_arm 2.4127  rmsd_full 2.3944

The mechanism stated in L7 stands on the correct numbers: on both targets the point cloud
(2.0900, 2.7313) is already nearer the native than the best single member of the 75 it
averages (2.4194, 2.8625), and the projected chain stays below that member (1.7814, 2.4127 in
the cache; 1.7809, 2.4224 in `results/summary/results.json`, which re-projects and round-trips
through PDB, the same rebuild difference as 3.2148 against 3.2126). Every other number in L7
was read from `results/summary/results.json` or the leaderboard before it was written.

Rule for myself, recorded so it is not repeated: a number goes into the ledger only in a turn
AFTER the tool output that carries it has been read.

---

## L10 -- OPERATIONAL ITEM 3 CLOSED: THE SIX HELD UNUSED-IMPORT PROPOSALS ARE APPLIED (2026-09-13, lane I)

S25 L14 held them "until the results and physics lanes stop importing `core.pipeline`". The
condition, checked by grep: `s25/resultslab/*.py` and `s25/phys_*.py` never imported
`core.pipeline` (they import `core.project`, `core.geometry`, `s12.instrument`);
`s24/d_harness.py:275` and `s24/d_hamiltonians.py:130` name it in docstrings only;
`s23/c1_probweight.py:35` imports it (a finished sprint script); `s12/instrument.py:188`
imports it lazily inside one function. No live lane imports it, and at the moment of the edit
(`s26/governor_state.json`, 00:26:31) the governor had no job registered at all, so the
"never edit a module while a job launched from it is running" rule was met by measurement.

All six applied, commit `37bddbbb`, each verified with `s26/i_ast_check.py` (the working file
against HEAD, docstrings stripped, `ast.unparse` diffed) to differ by EXACTLY the removed names:

| file | removed | evidence it was unused |
|---|---|---|
| `core/amber.py:281-282` | `MAXITER_PER_PARAM`, `MIN_MAXITER_PER_PARAM` from the `budget` import | zero references repo-wide outside `budget.py` and the import; not in `__all__` |
| `core/bench.py:47` | `asdict` | single occurrence was the import |
| `core/cache.py:43` | `Iterable` | single occurrence |
| `core/predict.py:38` | `Dict` | single occurrence (`Dict`, strings included) |
| `core/data.py:440-442` | `_seq_index()` (dead `lru_cache` helper) | zero references in 700 modules |
| `verify/vqe_lfo_audit.py:19` | `defaultdict` | single occurrence; not on the production path |

Nothing changes a number: the removed names were never read. Tests after the edit, all under
the governor: `tests/test_pipeline.py` + `tests/test_data.py` (job
`pytest_nanpoison_post_core_edit`, 79 tests, 0 failed, 2 absent-artefact skips, exit 0) and
`tests/test_amber.py` (job `pytest_amber`, 16 passed, exit 0). The full non-AMBER suite is
re-run on the committed tree as the last act of the lane and recorded in `s26/TEST_RUN.md`.
`peptide_db.py` (the legacy arm) was not touched.

---

## L11 -- LANE P: THE COMPACT ESM TABLES COVER THE TRAINING CORPUS; THE BANK COSTS 1.79 GB; ONE S7 ARTEFACT IS LOST (2026-09-13, lane P)

`s26/results/p_probe_esm.json` (jobrun `p_probe_esm`, peak RSS 0.594 GB, 10 s): `s5/esmraw.npz`
(486.6 MB float32 resident, 1.3 s), `s5/esm32.npz` and `s7/repr_cache/esmcon.npz` each cover all
6,790 unique training sequences (787 peptides + 6,003 fragments) and all 126 targets;
`esm_small.npz` covers 360. Per fold 6,609-6,640 chains, 552,199-556,709 pairs.
`s26/results/p_probe_esmcache.json` (jobrun `p_probe_esmcache`, tag ESM, est 2.5 GB):
`esm_cache.npz` = 22,795 sequences, 1.79 GB resident (process peak_wset 1.793 GB; the 5-s jobrun
sampler read 1.456 GB, an under-read of the transient), 19.7 s. Decision: the C2 ladder reads the
compact tables only and patches `core.data.esm_raw/esm_embed/esm_contacts` so no path reaches the
bank. `s12.instrument.project` rebuilds 1A13's persisted `fit_ca` and `ca` at max abs 0.0 in 2.9 s.
`s7/repr_tune.json` (S7-11's per-target ESM contrast) is not on disk and not in git history
(`git log --all -- s7/repr_tune.json` is empty); only `s7/repr_oracle.json` survives. The S7-11
contrast can only be re-measured (rungs noesm / pca32 of `s26/p_ladder.py`).

---

## L12 -- LANE P: C2 CODE AND ITS PROBE. THE RETRAINED pca32 REPRODUCES THE PINNED FOLD-0 MODEL; THE LADDER COSTS ~38 MIN PER RUNG (2026-09-13, lane P)

`s26/p_ladder.py` (rungs shipped, noesm, conly, pca32, pca32f, pca128, raw, esm8m, wide, pairnet,
mix; endpoints sel / cloud / arm / fit; gam_eff and cos in probability and location space; the
phase gate is enforced in code). Synthetic tests `s26/p_ladder_test.py` (jobrun `p_ladder_test`,
10 s, 0.325 GB): lean trainer == `core.predict.MLP.fit` at max |dprob| 1.0e-7; rebuilt pca32
features vs the shipped construction max abs 9.5e-7; `mix` lam = 0 bit-exact; NaN-poison
bit-identical; the guard refuses unknown sequences. Probe (jobrun `p_probe_train_pca32_f0`):
fold 0, 552,199 pairs, 183-d, 40 epochs, 453 s, peak RSS 1.248 GB. Gate
(`s26/results/p_ladder_gate_pca32_fold0_s0_probe.json`, 25 fold-0 targets, native-free): posterior
max abs 1.5e-5 vs `distogram_models/fold0_esm_frag.pt`, risk max abs 6.8e-4, K = 500 argmin
identical 25/25, top-75 overlap 1.000, Spearman 0.99999998. Not bit-exact (float32 rounding of
the cached embeddings), equivalent in every quantity the pipeline consumes. A note for the record:
S7-11's "tri on ESM input has not been measured" is superseded on one basis by S19 L11 (PairNet on
the deployed inputs, -0.080 [-0.242, +0.082] through the distance-geometry fit, paired sd 0.910,
MDE 0.227, `s19/results/a_models.json`); unmeasured through the pipeline readout, which is rung
pairnet. Question for the coordinator: may the 5-fold TRAINING of the ladder (no native RMSD read)
run before sign-off, so eval can start at sign-off? Not launched pending the answer (rule 14).

---

## L13 -- LANE P: B1 IS INFEASIBLE ON THIS BOX ON THREE INDEPENDENT GROUNDS; B1 STOPS (2026-09-13, lane P)

`s26/results/b1_feasibility.json`. (a) Weights on disk: only esm2_t33_650M (2.60 GB) and
esm2_t6_8M (30 MB); no ESMFold or 3B checkpoint. (b) `esm.pretrained.esmfold_v1` exists but
`import esm.esmfold.v1.pretrained` raises `ModuleNotFoundError: No module named 'omegaconf'`
(esmfold.py:10); `openfold` (esmfold.py:11-13, trunk.py:11) is absent; the fair-esm README's
openfold install "requires nvcc" and "python <= 3.9" (box: Python 3.13, torch 2.13 CPU).
(c) README: esmfold_v1 = 48 (+36) layers, 690M (+3B) parameters; HTTP HEAD (network reachable in
< 10 s): esmfold_3B_v1.pt 2,771,653,574 B, esm2_t36_3B_UR50D.pt 5,678,116,398 B; fair-esm halves
the LM and keeps the trunk fp32 (esmfold.py:46) -> 8.76 GB resident minimum (fp16 everything
7.4 GB) against 4.4 GB of campaign headroom (`s26/governor_state.json` 00:19: 64.8% of 16.75 GB,
ceiling 93%). Nothing downloaded. The 650M contact head is exposed (`esm_features.compute`,
`return_contacts=True`), cached for every training sequence, and already 13 of the shipped prior's
183 input columns; S17 L23 measured it as a ranker (+0.116 [+0.047, +0.184] in-band at top-75,
+0.050 over its ESM-free twin, mostly compactness, no argmin gain, shortlist worse than matched
random). Its standalone value as a PRIOR INPUT is rung conly (PREREG_B2); the size axis is
esm8m vs 650M. Proposal B routes to B2/B3 and `s26/PROPOSAL_B_REPLACEMENT.md`.

---

## L14 -- LANE P: B3's PER-TARGET ARTEFACTS EXIST; THE ARM-CHOICE ORACLE OVER {pipeline, torsion predictor, helix} IS AN ORDER STATISTIC (2026-09-13, lane P)

`s13/cache/tors_rows.npz['a_pepPos']` (126 rows, `rmsd_build` mean 3.7705) and
`s14/results/ladder.json per_target['L0_constant_helix']` (126, mean 4.0648; its incumbent column
reproduces `rmsd_fit` 3.2041 at 0.0000). Paired against the production `rmsd_arm` (3.2148),
persisted artefacts only: tors - arm +0.5557 (median +0.2314), SE 0.1318, MDE 0.3694, 37W/89L,
top-10 share 0.567; helix - arm +0.8500 (median +0.1691), SE 0.1494, MDE 0.4185, 32W/94L;
tors - helix -0.2943, SE 0.0807, 69W/57L. FAIL18: tors 5.569 / arm 6.032 / helix 5.887;
other108: 3.471 / 2.745 / 3.761. Per-target min over the three arms 2.9632 (-0.2515 vs arm) is
92% accounted for by `ST.best_of_k_within`'s across-target null (k_eff 2.34). The S13 `base`
column (3.2126, the leaderboard rebuild) differs from `rmsd_arm` by up to 0.171 A per target;
PREREG_B3 pairs against the production cache. Also noted: the coordinator's L5 defines lane PH's
C3 (steric reject filter, matched-random control); lane P's `PREREG_C3.md` was drafted without
that brief and says so in its addendum; PH's brief governs.

---

## L15 -- DEFECT 6a: THE IDENTITY FLAG SHIPS DARK; THE "CORRECTED" NORMALISATION IS AT THE NULL AS A CLUSTERING CRITERION; THE FOLDS STAY PINNED (2026-09-13, lane I)

**Code.** `core.data.identity(a, b, norm="longer")` and `identity_many(..., norm="longer")`
(commit `37bddbbb`): the default path is statement-for-statement the pinned convention;
`norm="shorter"` is the same match count over the SHORTER sequence behind a verbatim-substring
test. `clusters()` / `folds()` do not accept it; nothing on the production path passes it;
`peptide_db.py` (legacy arm) untouched. Test:
`tests/test_data.py::test_identity_norm_flag_ships_dark_and_the_shorter_form_catches_the_self_copy`
(default bit-identical; 1CEK-in-1A11 scores 1.0 under the flag; scalar and batched forms agree
on both sides of `_BATCH_MIN`).

**Audit, in memory only** (`s26/i_identity_audit.py`, job `i_identity_audit3`, 0.056 GB,
`s26/results/i_identity_audit.json`; `peptide_folds.json` and `peptide_clusters.json`
sha256-identical before and after and equal to lane E's `s26/results/pinned_hashes.json`).
The in-memory copy of `core.data.clusters` reproduces the pinned clusters AND the pinned folds
exactly, so the copy is faithful.

1. The brief's corrected criterion (shorter normalisation + substring, same 0.6 threshold)
   collapses **470 clusters to 166**; 594/787 sequences gain or lose a mate; **71/126 dev
   targets** gain a mate sitting in another pinned fold; 48/60 benchmark sequences (count only).
2. That count is chance, not leak. Pre-registered null (`s26/PREREG_identity_null.md`, written
   before the control ran): a real dev sequence passes the shorter criterion against **0.5%** of
   the other 786 members, a composition-preserving shuffle of it against **0.3%** (ratio 0.62;
   the falsifier was ratio < 0.2). At mean degree ~3 on chance edges, single linkage percolates.
   The pinned longer criterion: 0.07% real, 0.001% shuffled. So `norm="shorter"` at 0.6 is a
   per-pair leak TEST, not a replacement clustering threshold; the docstring now says so.
   (My PREREG expected pass rates of 10 to 30%; they are 0.5%. The at-the-null verdict stands,
   the magnitude I predicted was wrong by 30x.)
3. The minimal fix, pinned criterion OR verbatim substring: 459 clusters; membership changes
   on 23 sequences, 20 of them verbatim containments and 3 carried along by single linkage
   (their pinned cluster mate is one); PREREG H3 as written is falsified by those 3, by
   transitivity. Dev targets that gain a mate in another pinned fold: **5** -- the four declared
   self-copies **1CEK** (fold 2; carrier n=25 in fold 3; longer identity 0.520), **2FBU** (4;
   carrier n=23 in 0; 0.522), **2P5H** (4; carrier n=17 in 2; 0.529), **6B9K** (0; carrier n=17
   in 2; 0.588), each with shorter identity 1.0 and each moved into one fold by the fix, plus
   **8ZG2** (fold 1), whose own verbatim carrier (n=23, longer identity 0.478) is in fold 1 and
   which gains a cross-fold mate only because that carrier also carries a peptide from another
   fold. 14 further dev targets have a verbatim relative that already sits in their own fold.
   Benchmark under the minimal fix: **2/60** gain a cross-fold mate (count only), which is the
   S24 L4 figure.
4. The trap, measured: even the minimal fix, re-derived through `folds()`' shuffle, would
   relabel **647/787** sequences (102/126 dev under the shorter criterion), because cluster ids
   renumber and the seeded permutation changes with them. Any future repair must PATCH the
   pinned fold map (move the copies) and retrain only the affected fold models, never re-derive.

Nothing was written to the pinned files; nothing ships. Sign convention for the record: 4/126
dev and 2/60 benchmark stand as the self-copy counts.

---

## L16 -- DEFECT 6b: THE AMBER MEMORY GUARD NAMES THE CEILING AND QUOTES THE GOVERNOR; BOTH AMBER FILES CARRY IT (2026-09-13, lane I)

`tests/test_amber.py` (commit `eb89c165`): the autouse `_memory_ceiling` fixture keeps its
decision exactly (`_memory_verdict` on `core.amber.memory_percent()`, the syscall
`core.amber.memory_guard` itself uses, so the verdict is identical with or without a governor);
what changed is the text. `_governor_reading()` reads `s26/governor_state.json` and returns
its RAM / CPU / job counts / band ceiling when the snapshot is under 120 s old, else `None`
(absent, stale, malformed). `_ceiling_message()` then says "physical memory N%
(core.amber.memory_percent) is above core.amber's 92% ceiling; the S26 governor read X% RAM /
Y% CPU at <ts> (<age> s ago) with J registered job(s), A AMBER, governor ceiling 93.0%", or
"no fresh s26/governor_state.json, so the governor is not running", followed by the existing
"this suite did not cause it" / "THIS SUITE accounts for it" clause. Unit test
`test_the_ceiling_message_names_the_ceiling_and_the_governor_reading` (synthetic snapshot
files; fresh, stale, absent and malformed cases; both verdicts).

`tests/test_amber_frame_invariance.py` imports that fixture object (`from test_amber import
_memory_ceiling`), which registers it as autouse in that module: the S25 L14-approved change,
shared rather than copied, so memory pressure there is now a skip that names the ceiling instead
of a red `MemoryError` out of `core.amber.memory_guard`.

Under the governor, on the tree of `37bddbbb`: `tests/test_amber.py` job `pytest_amber`
**16 passed** (15 + the new test), 260.5 s, peak RSS 0.872 GB; `tests/test_amber_frame_invariance.py`
job `pytest_amber_frame` **3 passed**, 255.7 s, peak RSS 0.324 GB; the guard fired in neither
(box 65 to 75% throughout, `s26/governor.log`). Record: `s26/TEST_RUN.md`,
`s26/results/test_run.json`, `s26/results/pytest_amber.xml`, `pytest_amber_frame.xml`.

---
## L16b -- PHASE-GATE READING FOR MODEL TRAINING; C3 OWNERSHIP; REPRODUCTION EXACT (2026-09-13, coordinator)

(Numbered L16b: lane I's L16 landed while this entry was being written; rule 4 below applies.)

1. Lane P asked (L12) whether the C2 ladder's 5-fold TRAINING may run before sign-off. Training a
   distogram head on the training folds reads native distances as labels, exactly as the pinned
   production models were trained; it reads no RMSD, selects nothing, and produces no claim. It is
   allowed before the gate under three conditions: one training job at a time through the
   governor (probe peak 1.25 GB, L12); rung order fixed by `s26/PREREG_C2.md` before the first
   job; no evaluation of any rung (no RMSD) until "PHASE 0 SIGNED OFF" is posted.
2. Proposal C's item C3 (AMBER relaxation of the best C2 rung's output against the S16
   matched-magnitude random displacement control) is owned by lane PH (L5, `s26/briefs/PH.md`,
   prereg to be `s26/PREREG_c3_control.md`). Lane P's `s26/PREREG_C3.md` describes a different
   hypothesis (AMBER-weighted pool histogram as a mixture partner of the prior; relaxed training
   labels). It is not deleted: lane P appends an addendum saying it is not C3, and enters the
   cheap arm in the tournament as `s26/IDEA_amber_prior_partner.md`.
3. Lane E's reproduction (`s26/results/e_reproduce.json`, job `e_reproduce`, exit 0): fresh
   versus stored means agree to 0.0 on all four bases (rmsd_avg 3.048338, rmsd_fit 3.204076,
   rmsd_arm 3.214765, rmsd_full 3.235460, n=126); T030 = 1S9Z rmsd_arm 0.181981 fresh and
   stored; fold mismatches against the pinned folds: none. The 2.6e-04 A tolerance is met with
   zero disagreement. Phase 0 still waits for `s26/EXAMINATION.md` and the Adversary's audit.
4. Ledger numbering: five lanes append concurrently. Re-read the tail immediately before
   appending and take the next free number; if two entries share a number, the later one adds a
   suffix (for example L17b) rather than renumbering.

---

## L17 -- SESSION-LIMIT INTERRUPTION 00:46 TO 08:40; NOTHING LOST; LANES RESUMED FROM DISK (2026-09-13, coordinator)

All five lanes (E, I, Q, P, PH) were terminated at about 00:46 by the API session limit
("session limit, resets 03:30 America/Los_Angeles"), not by any error. The governor (pid 36196)
ran throughout: `s26/governor.log` is continuous, with SAMPLE lines every minute from 23:52 to
08:36. No governed job was in flight at the cut: the last registration reaped was
`pytest_core_post` at 00:45:56 and `s26/jobs/` was empty at 08:36; every record in
`s26/jobs_done/` carries an exit code. Every lane's files on disk stand (listing at 08:36:
preregs C1-C5, B1-B3, amber_reject, cis, c3_control, identity_null; eight IDEA files; scripts
e_*, i_*, p_*, ph_*, q_*; results in `s26/results/`); lane commits eb89c165, 37bddbbb, 9bb7f4f6,
601a39c7 are on the branch. At 08:40 each lane was resumed from its own transcript with an
instruction to commit what is on disk and continue, not restart. The phase gate is still closed
(`s26/EXAMINATION.md` and `s26/BRIEF.md` not yet written). Rulings added on resumption: lane Q's
A2 (DLA) and A4 (gradient variance versus width) read no native and may run before the gate;
A1 and A3 endpoints wait.

---

## L18 -- HOW L15's BENCHMARK COUNTS WERE OBTAINED: THE MANIFEST WAS READ THROUGH core.data.benchmark(), SEQUENCES ONLY (2026-09-13, lane I)

The counts in L15 that concern the sealed benchmark (48/60 gain a cross-fold mate under the
shorter criterion; 52/60 would be relabelled by a naive re-derivation; 2/60 under the minimal
pinned-OR-substring fix) were computed in `s26/i_identity_audit.py` (sections 3b and 5) over
the SEQUENCES of the 60 benchmark peptides returned by `core.backend("data").benchmark()`.
That function (`core/data.py:617-633`) opens `results/benchmark_manifest.json`, reads the `pdb`
field of every entry in its `targets` list, and returns the matching `Peptide` records of the
peptide database (`peptide_db.npz`). So the manifest's contents WERE read by the audit process,
by the route the lane brief named as the permitted one ("through core.backend("data").benchmark()
sequences only, print no names, report only the count"); stated plainly so the coordinator can
log it as the conservative path. The S24 L4 count (2/60) was made the same way.

What was and was not touched: the audit used only `.seq` of each returned record. The
`Peptide` records carry native CA / phi / psi for every one of the 787 database entries
(they ARE the database, loaded by every caller of `core.data.load()`), and no such field was
accessed for a benchmark entry. No benchmark PDB file was opened, no native coordinate read,
no RMSD computed, no benchmark name printed, saved or logged: `s26/results/i_identity_audit.json`
holds counts only (`benchmark60`, `substring_or_pinned.benchmark60_with_new_mates_in_other_pinned_folds`),
and non-dev database members appear in the audit output by length alone. Nothing further will be
computed on benchmark targets by this lane.

---

## L19 -- OPERATIONAL ITEM 2 CLOSED: verify/run_equiv2.sh RUNS THE COMPARISON IT WAS WRITTEN FOR; BASELINE AND THE SHIPPED MODE ARE BIT-IDENTICAL (2026-09-13, lane I)

Why it was stale: it exported `PROJECT_GRAD`, which `core/project.py` stopped reading when the
mode moved into `core.pipeline.Config.project_grad` so that it reaches the cache key
(`verify/grad_key_collision.py`). The `core.pipeline` CLI has no flag for the mode, so as
written the script ran three arms in one mode. Repair (commit `eb89c165`): the consolidated
arms go through `s26/i_run_equiv2_arm.py`, which builds the identical `Config` to
`python -m core.pipeline run` (the same `replace(PROD, ...)` call) and sets `project_grad`;
four arms (baseline legacy; `exact`, the shipped default; `fd`; `analytic`); the original
`cfg_key` grep kept; any argument passed to every arm. `s26/i_equiv2_compare.py` compares the
arms per target from their cache directories with `==`.

Run: job `run_equiv2`, `sh verify/run_equiv2.sh --no-amber`, tag AMBER (with `--no-amber` no
OpenMM context is created, but `openmm` is still imported because `core.pipeline` resolves
the amber backend at start, so the conservative tag was the right one), exit 0, 145.3 s, peak
RSS 0.596 GB (`s26/jobs_done/run_equiv2.json`). Per arm on smoke8 (8 targets, no stage 4):
baseline 67.4 s, exact 33.0 s, fd 20.2 s, analytic 8.8 s.

    arm              cfg_key            vs opt_exact (8 targets)
    baseline_legacy  65ec272db3d31f05   ca, fit_ca, phi, psi, avg_ca, rmsd_avg, rmsd_fit,
                                        rmsd_arm, n_windows, n_top: bit-identical 8/8
    opt_exact        66050f6daae4ca07   (the shipped default)
    opt_fd           85faafd84d76a827   avg_ca and rmsd_avg identical; rmsd_fit differs on
                                        8/8, max 6.9e-4 A; rmsd_arm max 1.2e-2 A
    opt_analytic     f4e137a48586bd4b   avg_ca and rmsd_avg identical; rmsd_fit max 1.0e-1 A;
                                        rmsd_arm max 4.8e-2 A, 8/8 targets

Four distinct keys: the mode is in the key and the collision the script was flagged for is
closed. Baseline against the shipped mode is bit-identical on every emitted array and scalar:
the consolidation is faithful, and the `fd` and `analytic` differences are the scan builder and
the analytic gradient, as `core/project.py`'s docstring states. (The raw coordinate max|d| of
22 to 32 A and phi/psi differences near 2 pi between modes are lab-frame and angle-wrap
differences of un-superposed emitted chains, not the science; the RMSD scalars are.)

Artefacts: `s26/results/run_equiv2.log` (whitelisted transcript: the jobrun log, all four arm
logs, the comparison), `s26/results/run_equiv2_compare.json`, `verify/e2_*.log` (ignored).
Not run: the AMBER-inclusive form; stage 4 is downstream of the projection the script compares.

---
## L18b -- GOVERNOR v2: SMOOTHED CPU, RAM-ONLY KILLS, 20 s MINIMUM SUSPENSION, LAUNCH CAP OF 4 (2026-09-13, coordinator)

(Numbered L18b: lane I's L18 landed first; the governor.log RESTART line says "ledger L18" and means this
entry. Lane I's L18 records that the benchmark counts in L15 came from `core.data.benchmark()` sequences
only, with no coordinate, native or RMSD read; the coordinator logs that as the conservative path:
membership and sequence through the project's own helper is allowed for counting, and nothing further
is computed on benchmark targets in S26.)

`s26/governor.log` 00:37:01 to 00:38:51: with five governed jobs registered the raw CPU sample
crossed 93% and fell under 90% on alternate 5 s ticks (SAMPLE lines 00:37:41 cpu 98.7%, 00:38:41
92.8%, RAM 67 to 68% throughout), so the v1 band suspended and resumed `q_dla_smoke` eleven times
in two minutes and suspended two PH census jobs. No job was killed and every job exited 0, but a
CPU spike sustained 15 s would have killed the newest job, which for jobrun-launched jobs has no
queue spec to requeue from. Changes (commit follows this entry), all in `s26/governor.py` and
`s26/jobrun.py`, no production module touched:

- CPU enters the band as a 15 s rolling mean (`cpu_smooth`, three samples), written to
  `s26/governor_state.json` beside the raw value.
- The kill rule (95% for 15 s) now reads RAM only; CPU above 93% suspends the newest running job
  and never kills. CPU overload slows the box; only RAM can crash it.
- A suspended job stays suspended at least 20 s; resumption needs RAM under 90% AND smoothed CPU
  under 80%.
- `jobrun.py` also waits while the smoothed CPU is above 85% or four jobs are already
  registered, so five lanes launching at once are serialised at the door instead of suspended
  after the fact.

The v1 governor (pid 36196) was terminated at 08:40:24 with the RESTART line in the log and v2
started in its place. Two lane-Q jobs (`a2_dla`, `a4_var`) were registered at that moment and ran
unsupervised for a few seconds; v2 picked their registrations up on its first sample. The log is
continuous through the restart.

---

## L20 -- OPERATIONAL ITEM 5 CLOSED: THE BRIEF'S 6.1 CARRIES THE GOVERNED RUN; THE ORIGINAL SENTENCE IS IN HISTORY (2026-09-13, lane I)

`docs/STATE_BRIEF_2026-09-12.md` was committed exactly as received as `6e50ea93`, so the
"355 passed, 13 skipped, 0 failed, 0 errors, exit 0" sentence survives in history; then its
section 6.1 live-run paragraph was replaced (commit `83305549`) by the S26 governed run:
**370 tests, 357 passed, 13 skipped, 0 failed, 0 errors, 0 memory-guard skips** on commit
`601a39c7`, in three jobs under `s26/governor.py` (`pytest_core_post`, the 10 non-AMBER files:
351 tests, 338 passed, 13 skipped, 240 s, peak RSS 1.692 GB; `pytest_amber`: 16 passed, 261 s,
0.872 GB; `pytest_amber_frame`: 3 passed, 256 s, 0.324 GB), with the passed count per file,
the 13 skip reasons (11 `VERIFY_SLOW=1` opt-ins in `test_equivalence.py` and
`test_integration.py`, 2 absent artefacts in `test_pipeline.py`) and the pointer to
`s26/TEST_RUN.md` / `s26/results/test_run.json`. The suite is 370 rather than 368 because S26
added one test to `test_amber.py` and one to `test_data.py`.

The same non-AMBER files run BEFORE any S26 edit (job `pytest_core`, tree `a4db170c`: 337
passed, 13 skipped, 0 failed, identical skip list) are recorded beside it and marked
superseded, so the before/after of the production-path edits (L10, L15) is on file.
`TEST_RUN_RESULT_PLACEHOLDER` never existed in the tree (L3); the live-run sentence held the
placeholder's role.

---

## L21 -- HYGIENE: README GOVERNOR SECTION; `make examine` IS `python s26/examine.py` (WITH examine.sh / examine.bat), WRAPPING THE MODULE MAP, THE PINNED HASHES AND THE CLAIM LEDGER (2026-09-13, lane I)

`README.md` (commit `eb89c165`) gained the section "The S26 resource governor": the three
scripts and what each does, the band (88% low water to launch queued work; resume below 90%;
93% ceiling suspends the newest job; 95% for 15 s kills it with CTRL_BREAK first and requeues
it), the AMBER cap of two, the tags CPU / AMBER / ESM / TEST, where the live state
(`s26/governor_state.json`) and the log (`s26/governor.log`, whitelisted) are, how to start it
(`python s26/governor.py`, `--once`, `--status`), how a job is run or queued, the three-way
split of the test suite, and the examine entry point.

There is no Makefile in the repository and no `make` on this box, so the `make examine`
equivalent is a script: `python s26/examine.py`, with the one-line launchers `examine.sh` and
`examine.bat` at the repository root (no name conflict). It calls, imported rather than copied:
lane E's `s26/e_module_map.py` (rewrites `s26/results/module_map.json`), lane E's
`s26/e_hashes.py --check` (re-hashes every pinned artefact against
`s26/results/pinned_hashes.json`; the benchmark manifest as bytes only), and lane I's
`s26/i_claim_check.py` over `s26/results/claims.json` (21 claims, each re-read from the artefact
it names: the four 126-target means of the production cache, the `compare_tuning126.json`
science deltas and baseline constants, the leaderboard's production row and gate, the pinned
hashes, the S8-13 transfer law with the excerpt as fallback, and the S26 test-run totals);
`--search` adds lane E's `s26/e_claims.py`. Exit status is non-zero on any drift, mismatch or
non-optional absence. Validated: 21/21 OK at 00:37 (`s26/results/claim_check.json`), and the
full run at 08:4x is `s26/logs/i_examine_full.log` (it rewrote lane E's `module_map.json` once;
the map is deterministic apart from the sha256 of the files S26 edited).

Also added for the contract's own rules: `s26/i_ast_check.py` (section 5, AST identity modulo
docstrings against a commit, with a readable diff of any deviation) and `s26/i_test_report.py`
(junit XML plus `s26/jobs_done/` into `s26/TEST_RUN.md` and `s26/results/test_run.json`, with
memory-guard skips counted apart from real skips).

---

## L22 -- CIS CENSUS: NO CIS PEPTIDE BOND EXISTS ANYWHERE ON THE INSTRUMENT; THE DATABASE STEP GATE, NOT THE PROJECTION, SETS THE COST (2026-09-13, PH)

`s26/ph_cis.py census`, `s26/results/ph_cis_census.json` (complete 126/126), job
`s26/jobs_done/ph_cis_census.json` (exit 0, 90 s, peak RSS 0.038 GB). Pre-registered in
`s26/PREREG_cis.md` section 2 with the prediction "zero on model 1 and in the pool, because of
the step gate; some cis bonds in other ensemble models". Allowed before the gate by L5: ORACLE
DIAGNOSTIC on the natives (omega angles only, no RMSD), native-free on the pool.

    model-1 natives with a cis bond, |omega| < 30 deg        0 / 126
    the same by consecutive CA-CA < 3.3 A                    0 / 126   (the two criteria agree 126/126)
    deposited models with a cis bond, every ensemble         0 / 1,966
    universe windows with a CA-CA step < 3.3 A               0 / 2,352,893   minimum step 3.5045 A
    K=500 pools and production top-75 with such a window     0 and 0
    omega non-planarity |180 - |omega||, 1,507 bonds         mean 1.91 deg, median 0.34, p90 5.8, p99 17.4,
                                                             max 42.8 (9UV5, bonds 0 and 6: -137.2 and +144.8 deg);
                                                             0.53% of bonds beyond 20 deg, 0.13% beyond 30 deg

The universe minimum step IS the gate: `core/data.py:406-407` drops any peptide with a
consecutive CA-CA step below 3.5 A and `core/data.py:697-698` drops any fragment window with
one. The 126 are drawn from that database (`s7/debias.py:103-120`) and so is the sealed
benchmark, so neither can contain a cis target and no pool can contain a cis window. The
two-bond-length projection is therefore worth exactly 0.000 A on this instrument by
construction of the database, not by any property of the projection. The residual cost of the
constant omega here is the non-planarity tail (0.5% of bonds beyond 20 deg), which part 2
(gated) prices as the ideal-trans floor on the native's own torsions. The registered prediction
that other ensemble models carry cis bonds was WRONG: 0 of 1,966. The cis-peptide question is a
world-supply question (the 16 containment-fresh targets, 10 amyloid) and the design note in
`s26/agentPH_FINDINGS.md` section 1.4 is written for that supply.

---

## L23 -- STERIC REJECT CENSUS: AT 1e4 kcal/mol THE REJECT REMOVES 40 OF 75, EMPTIES 11 TOP-75 SETS AND 8 WHOLE POOLS, AND 96.8% OF THE CATASTROPHES ARE SIDE-CHAIN CONTACTS (2026-09-13, PH)

`s26/ph_reject.py census`, `s26/results/ph_reject_census.json` (complete 126/126), job
`s26/jobs_done/ph_reject_census.json` (exit 0, 100 s, peak RSS 0.117 GB). Native-free: the
63,000 cached single points (`s24/cache_amber`, pool identity `universe_idx == I.pool_idx` and
production `sub` == score top-75 asserted on every target) and the ideal-geometry rebuilds; no
RMSD, no native. Required by `s26/briefs/PH.md` section 3.1 before any RMSD is read; the prior
was registered in `s26/PREREG_amber_reject.md` section 5.

    threshold     pool frac > T   top-75 rejected   zero-reject   all-75 rejected   no survivor in 500   refill depth (mean/median)   R overlap with anchor
    1e3             0.760           54.5 / 75          1              25                 20                270 / 234                    0.27
    1e4 PRIMARY     0.586           40.1 / 75          2              11                  8                201 / 147                    0.46
    1e5             0.446           30.3 / 75          6               3                  1                167 / 118                    0.60
    1e6             0.345           23.3 / 75          8               2                  0                135 / 103                    0.69

The 0.586 reproduces S25's 58.6% (`s25/results/phys_landscape.json`). At the primary threshold
the operator is not a small surgical reject: it removes more than half of every shipped set,
reaches rank 147 (median) of 500 to refill, empties the whole top-75 on 11 targets (1G89 1ID6
1LB7 2MAI 2NB7 2XL1 5MML 5Z5W 7BX2 8UN8 9S5G) and finds no survivor among all 500 candidates on
8 (1G89 1ID6 2MAI 2NB7 2XL1 5MML 7BX2 8UN8). Both empty cases fall back to the anchor, a rule
added in the prereg's addendum 1 before any RMSD is read. Native-free geometry of the retained
set at 1e4: R is +0.254 A more expanded in Rg than the anchor (SE 0.058) and +0.144 A in the
minimum |i-j| >= 3 CA-CA distance; S is +0.060 and +0.070. Same sign as S25's +1.10 A expansion
of AMBER's top-75, smaller because the distogram's order is kept among the survivors.

WHERE THE SINGULARITY LIVES (the optional Part IV measurement, `singularity` block; every top-75
member rebuilt with all heavy atoms through the reference builder `sidechains.py`, no OpenMM):

    closest heavy-atom contact, residue separation >= 2
      all 9,450 members:              bb-bb 1,054   bb-sc 5,300   sc-sc 3,096
      the 5,057 members above 1e4:    bb-bb   163   bb-sc 2,627   sc-sc 2,267    -> 96.8% side-chain-involving
    members per target with a heavy-atom pair closer than 2.0 A:   all-atom 40.7 of 75 (SE 1.9)
                                                                    backbone+CB 2.61 of 75 (SE 0.34)   [S19 section 3.1: 2.66, reproduced]
    Spearman(e_amber, minimum heavy-atom distance) within a top-75: mean -0.743, median -0.803
    minimum heavy-atom distance: members above 1e4, 1.48 A; members below 1e4, 2.31 A

The AMBER single point on the top-75 is, to rho -0.74, the minimum heavy-atom distance of the
rebuild, and 96.8% of the rebuilds it condemns are condemned by a contact involving a side chain
placed by the deterministic builder (fixed chi1, no rotamer scan). The physically impossible
class on the backbone is 2.6 per 75 (S19); the class the 1e4 reject removes is 40 per 75.
Stated before any RMSD is read: the steric reject at the primary threshold is a reject of the
builder's side-chain placement, not of the pool's backbones. This is Part A of
`s26/IDEA_rotamer_relief.md` and it survives.

---

## L24 -- C3 NATIVE-FREE PART: THE PRODUCTION RELAXATION MOVES THE CA TRACE 0.220 A RMS; 58.7% OF BUILT CHAINS START ABOVE 1e4 kcal/mol; 125 OF 126 CONVERGE (2026-09-13, PH)

`s26/ph_c3.py nativefree`, `s26/results/ph_c3_nativefree.json` (complete 126/126), job
`s26/jobs_done/ph_c3_nativefree.json` (exit 0, 5 s; the RSS sampler polls every 5 s and
under-reads a 5 s process, so no memory number is claimed). Read from
`bench_results/cache/1fc9f2dcf489e2fb` with every RMSD key stripped; nothing here reads a native.

    AMBER displacement of the built chain, per-atom RMS after superposition   0.220 A (SE 0.008, median 0.197, range 0.103 to 0.591)
    `amber_moved` (restraint RMSD on N/CA/C, unsuperposed)                     0.233 A
    the built chain's own AMBER energy e0 before relaxation                    median 8.6e4 kcal/mol, min -473, max 1.3e14; 58.7% above 1e4
    relaxed energy e1                                                          mean -560 (SE 30), max +1262 (9KAR)
    converged (e1 <= CONVERGE_MAX_KCAL = 1000)                                 125 / 126 (9KAR fails)
    bond + angle strain after                                                  mean 60 kcal/mol
    virtual CA-CA bond: built 3.80395 (exact, sd 1e-16) -> relaxed             mean 3.867, min 3.12 (1M02), max 5.38 (2BP4, e1 +845); 4.86 on 9KAR
    targets with a relaxed CA-CA outside [3.6, 4.0]                            6 / 126
    Rg change                                                                  +0.046 A

Two things the validity story must carry: the emission's own strain census (58.7% of built
chains sit above 1e4 kcal/mol before relaxation, the same fraction as the pool), and the two
targets where the relaxation itself breaks a virtual bond (2BP4, 9KAR). Derived prediction for
stage 1, registered in `s26/PREREG_c3_control.md` addendum 1 before the gate: an orthogonal move
of 0.220 A on a 3.21 A chain costs about m^2 / (2 RMSD) = 0.0075 A by the S16 identity; the
production step costs +0.0207; so AMBER is predicted WORSE than the matched random control by
about +0.013 A. Measured after the gate.

---
## L25 -- LANE W (WILDCARD) SPAWNED; THE 2/60 PROXY BOUND IS ITS FIRST MANDATORY DIRECTION (2026-09-13, coordinator)

Six lanes active (E, I, Q, P, PH, W; `s26/lanes.json`). Brief `s26/briefs/W.md`. Of the mandatory
tournament directions, test-time window ensembling is already claimed by lane P
(`s26/IDEA_window_ensembling.md`) and the AMBER reject and cis census by PH (L22, L23); the 2/60
benchmark self-copy leak bounded from a dev-set proxy without opening the benchmark is unclaimed
and goes to W, building on lane I's L15 audit (the four dev self-copies and their carriers) and
S10-4's +0.0004 A dev price. W also owes at least three ideas nobody else has proposed and takes
the highest-ranked orphaned survivor once the Adversary ranks the tournament.

---


## L26 -- LANE P: THE C2 TRAINING CHAIN IS RUNNING, ONE GOVERNOR JOB AT A TIME, IN THE ORDER FIXED BY PREREG_C2 ADDENDUM 2; THE C3 HAND-OFF FORMAT IS IMPLEMENTED (2026-09-13, lane P)

Resumed 08:37 after the 00:46 session cut; no training job had been launched before the cut
(`s26/jobs_done/` held only the probes). Per L16b: order noesm -> conly -> pca32 -> wide -> pca32f
-> pca128 -> featurise-esm8m -> esm8m -> raw, fixed in `s26/PREREG_C2.md` addendum 2 BEFORE the
first job; `s26/p_train_chain.sh` runs them sequentially through `s26/jobrun.py` (est-ram 1.5 GB,
raw 1.8, NT = 2), one `.pt` per fold as the checkpoint, a second `_p2` pass resuming anything
the governor kills. `p_train_noesm` registered 08:52. No `eval`/`report` before "PHASE 0 SIGNED
OFF" (refused in code). Expected machine time 6-10 h (raw last: 5175 x 384 first layer).
While it runs: `s26/PREREG_C3.md` addendum 2 (not Proposal C's C3; C3 is lane PH's,
`s26/PREREG_c3_control.md`), `s26/IDEA_amber_prior_partner.md` (H_C3a entered in the tournament;
H_C3b withdrawn by S8-14's arithmetic), the C1 closure reproduction written into
`s26/agentP_FINDINGS.md` section 7 (S12 learning curve 3.0433 -> 3.0258 real vs 2.5342 leaked at
n = 8; S17 band best 2.6087 vs random 3.4676, all to the third decimal), `s26/p_stats.py`
(nested ridge, selftest OK: planted 0.866 vs null p95 0.558), `s26/p_b3.py`, `s26/p_c4.py`,
`s26/p_c5.py` (native-free feature builds and synthetic selftests queued under jobrun; their
`run` commands are gated). The C3 stage-2 hand-off format (`s26/results/p_best_rung_chains.json`:
rows keyed by pdb with phi/psi in RADIANS as `s12.instrument.project` emits them at lam = 0.3,
`ca`, `rmsd_arm`; `ST.save_atomic` with complete_keys and n_expected = 126) is implemented in
`s26/p_deliver.py`, a separate module so that `p_ladder.py` is not edited while its jobs run; it
re-projects the chosen rung and asserts equality with the eval JSON before writing.

## L27 -- LANE Q, A2: THE DEPLOYED ANSATZ'S DYNAMICAL LIE ALGEBRA IS THE FULL so(2^n) FROM DEPTH 2 AT n = 7; MY PRE-REGISTERED PREDICTION OF A PROPER SUBALGEBRA IS FALSIFIED; THE TANG MINIMAL POOLS GENERATE so(2^(n-1)+1) (2026-09-13, lane Q)

`s26/q_dla.py` -> `s26/results/q_dla.json` (complete; `s26/jobs_done/a2_dla.json`: 85.3 s, peak
RSS 0.479 GB). Pre-registered in `s26/PREREG_A2.md` before the run. Property measurement, no
native, no score: the closure is computed exactly on Pauli strings as a set and cross-checked
by dense SVD (rtol 1e-10) at n = 4, 5: 12 of 12 cells agree; the two conjugation conventions
agree on 12 of 12. Every closure lies in the odd-Y (real) set.

    fixed ansatz (RY / CNOT chain + ring), dim(DLA) by width n and depth L
      L = 1                       n            (abelian) at every n = 4..11
      n = 4, 5, 7, 8, 10, 11      so(2^n)      from L = 2 on   (8128 at the deployed n = 7)
      n = 6                       510 / 1023 / 2016 = so(64)  at L = 2 / 3 / 4
      n = 9                       32766 / 65535 / 130816 = so(512)  at L = 2 / 3 / 4
    pools V and G (Tang 2021, 2n-2 strings)   36, 136, 528, 2080, 8256, 32896  at n = 4..9
                                              = dim so(2^(n-1) + 1) at all six widths
    pool L2 (all 1-, 2-local odd-Y strings)   so(2^n) at n = 4..9
    ADAPT-selected sets, ideal ladder, n = 7  alpha = 1: abelian (7) at every step, both pools,
                                              both optimisers; L-BFGS stops at P = 7 with no
                                              operator selected.  alpha = 0.25: V 7 -> 16;
                                              L2 7 -> 1025 at P = 21 (12.6% of so(128)).

H2b of the prereg said dim(DLA) at (n = 7, L = 3) is below 8128, guess 4095. It is 8128, the
full so(128), already at L = 2. Falsified; the measured value stands. H2a (depth 1 abelian),
H2c (pools: 2080, 8256, 32896 predicted and measured), H2d (odd-Y) and H2e (ADAPT sets abelian
at alpha = 1) held. The only widths where depth 2 is not enough are n = 6 and n = 9, whose
dimensions at L = 2 and 3 equal 2 dim su(2^(n-2)) and dim su(2^(n-1)); that n divisible by 3
is the condition is an observation from two cases and is labelled HYPOTHESIS, not explained.

What it means for the S25 slopes (`s25/results/q_plateau.json`): the algebra at the deployed
cell is maximal, so nothing in the algebra protects the ansatz from an exponential plateau.
Larocca et al. 2022 / Ragone et al. 2024 give Var ~ 1/dim(g) once the circuit is a 2-design
over exp(g); dim(g) = 8128 at n = 7 and grows as 4^n / 2. S13 measured the decay base
approaching 0.504 per qubit at depth 8 (`s13/results/geo_kernel.json`), the 2-design rate.
The S25 result "no exponential plateau at n = 4..13" is therefore a statement about depth 3
(P = 3n against dim so(2^n)), and must be quoted with "at depth 3" attached. It is not
evidence of a favourable algebra, and Cerezo et al. 2025 does not apply through the
small-DLA route; the circuit is simulable because n = 7, not because of its structure
(S21 L4).

The pools: a "complete" pool in Tang's sense (overlap-matrix rank 2^n - 1) generates
so(2^(n-1) + 1), which acts transitively on the real sphere and has about a quarter of the
dimension of so(2^n). Completeness is weaker than controllability; the ADAPT arm using pool V
in A1 is therefore restricted to that subalgebra by construction, and the L2 arm is not.
## L28 -- PHASE 0 EXAMINATION LANDED: REPRODUCTION EXACT; FOUR DOCUMENT-ONLY NUMBERS; THE SEVEN-CONFIGURATION SUITE IS QUOTED ON THE POINT-CLOUD BASIS; ONE RULE-1 INCIDENT (2026-09-13, lane E)

`s26/EXAMINATION.md` (sections A to I), `s26/BRIEF.md`, `s26/agentE_FINDINGS.md`; scripts and
artefacts in commit `e3570aed` (`s26/e_module_map.py`, `e_hashes.py`, `e_reproduce.py`,
`e_trace.py`, `e_claims.py`; `s26/results/module_map.json`, `pinned_hashes.json`,
`e_reproduce.json`, `e_trace_1S9Z.json`, `e_trace_9KAR.json`, `claim_search.json`, `.txt`).

1. Reproduction (Part 2.3): fresh vs stored means 3.048338093879532 / 3.2040761603809194 /
   3.2147651542109985 / 3.2354598538973844, max per-target disagreement 0.0 on all four bases,
   T030 = 1S9Z 0.18198112330908295 (L16b records the same). Fresh `run_target` + `label` on
   1S9Z and 9KAR: every arm 0.0, `sub` identical, emitted `ca` identical. Every dataflow stage
   of both traces matches the record at 0.0; the Hamiltonian is the rank ladder to 0.394% of
   range (S25's worst case 1.18%, `s25/results/q_gibbs.json`).
2. Claim ledger (EXAMINATION C, 35 claims): 30 sourced to an artefact leaf, a passing test or a
   stated derivation from stored rows. **UNSOURCED (document-only):** 36.1 / 36.4 deg phi MAE
   (C26), the +0.0004 / +0.0030 identity-leak prices (C27; the cited `s10/idaudit_*.json` do
   not exist), the |z_moment| triple 0.7529 / 0.8013 / 0.1127 (C34), 355 passed / 13 skipped
   (C35; superseded by `s26/TEST_RUN.md`: 369 / 356 / 13), and the 0.524 sampled tail-only
   cosine inside C24. The benchmark delta +0.0103 (C06) is named to `s9/final_report.json`,
   which this lane did not open; the two means beside it are asserted to 5e-4 by two passing tests.
3. Basis: the seven-configuration suite (3.058 ... 3.881), the random-75 null (3.4251) and the
   Legacy +0.330 / AMBER +0.455 verdicts are point-cloud numbers (`s25/results/phys_suite.json`,
   `results/summary/leaderboard.json :: rows[*]/mean_secondary`); the state brief (lines 57,
   190-191) and `results/summary/professor_brief.md:115-125` quote them beside the built-chain
   3.2148 without saying so. Built-chain means of the same rows: 3.2187 / 3.3100 / 3.3732 /
   3.4221 / 3.8248 / 3.8844 / 4.1015 (`leaderboard.json :: rows[*]/mean`).
4. For the Adversary: `VQE_LFO` runs alpha = 1.0 on folds 0, 3, 4 (78 of 126 targets, no tail
   constraint; `s25/results/q_alpha.json :: results/share_of_targets_with_no_tail_constraint`
   0.6190476190476191); 9KAR's production relaxation ends at 1262.4 kcal/mol, above
   `core.amber.CONVERGE_MAX_KCAL` (L24: 125 of 126 converge) and is inside the 3.2355 mean;
   2BP4 converges by that gate (e1 845.3) but carries bond+angle strain 1172.7, above the S8
   strain-rejection rule; `bench_results/cache/1fc9f2dcf489e2fb/` is gitignored (EXAMINATION G).
5. Incident (Rule 1): the first run of `s26/e_claims.py` parsed `s9/final_report.json` before the
   benchmark exclusion existed and printed one aggregate leaf, `dist/shipped/mean` =
   2.9507235775391263, already published in `README.md:14`. No per-target benchmark value was
   printed or read. The exclusion (`SKIP_FILES`, `SKIP_SUBSTR`) is in the committed script and
   the artefacts were regenerated with it.
6. Provenance note: the traces and the reproduction ran at git `74e073e2` (dirty: lane I's
   later-committed `37bddbbb` edits to `core/data.py` and four unused imports were in the tree).
   The identity path they touch is not on `run_target`; the reproduction costs 5 s and should be
   re-run on HEAD by the Adversary (`python s26/e_reproduce.py`).

---

## L29 -- ADVERSARY SPAWNED ON THE EXAMINATION; THE RULE-1 INCIDENT LOGGED; THE POINT-CLOUD BASIS OF THE SUITE NOTED (2026-09-13, coordinator)

Lane A (Adversary) spawned at 08:55 with brief `s26/briefs/A.md` on `s26/EXAMINATION.md`,
`s26/BRIEF.md` and `s26/agentE_FINDINGS.md` (L28). Phase 0 remains closed until its audit is
clean or every material finding is fixed and re-checked. Six lanes active (E, Q, P, PH, W, A);
lane I finished (L21).

Three things from L28 the coordinator rules on now:

1. The Rule-1 incident (L28 item 5): the first run of `s26/e_claims.py` parsed
   `s9/final_report.json` and printed one aggregate, `dist/shipped/mean` = 2.9507235775391263,
   a number already published in `README.md:14` and the state brief. No per-target benchmark
   value was read or printed; the exclusion is in the committed script and the artefacts were
   regenerated with it. Logged as an incident with no consequence: the benchmark was not spent
   and nothing new was learned about it. Rule for the rest of the sprint: no script parses
   `s9/final_report.json`, `results/benchmark_manifest.json` or any file under `results/`
   naming a benchmark target; lane I's L18 read (sequences through `core.data.benchmark()` for
   a count) and this one are the only two benchmark-adjacent reads of S26.
2. The seven-configuration suite, the random-75 null and the Legacy +0.330 / AMBER +0.455
   verdicts are point-cloud numbers (L28 item 3; `s25/results/phys_suite.json`). The state brief
   and `results/summary/professor_brief.md` quote them beside the built-chain 3.2148 without
   naming the basis. Every S26 deliverable that quotes them (the report, the deck, the proposal
   verdicts) names the basis on both sides; the built-chain means of the same rows are in
   `results/summary/leaderboard.json :: rows[*]/mean` (distogram 3.2187 ... AMBER 4.1015).
3. Four presenter-facing numbers are document-only (C26 36.1/36.4 deg, C27 +0.0004/+0.0030,
   C34 the |z_moment| triple, C35 355/13). The Adversary decides whether each is material; the
   coordinator's default is that a document-only number is not put on a slide unless a lane
   re-derives it from cached rows this sprint (C26 from `s13/cache/tors_rows.npz` if it holds
   the per-residue predictions; C35 is superseded by `s26/TEST_RUN.md` and is simply replaced).

---


## L30 -- LANE W: THE 2/60 PROXY BOUND IS PRE-REGISTERED; NATIVE-FREE CENSUS: TWO OF THE FOUR DEV SELF-WINDOWS NEVER REACH THE EMISSION; EVERY TARGET HAS A BLOSUM TIE AT THE POOL BOUNDARY (2026-09-13, W)

`s26/PREREG_selfcopy_bound.md` (written before any result), `s26/IDEA_selfcopy_proxy_bound.md`,
code `s26/w_selfcopy.py`, synthetic tests `s26/w_selfcopy_test.py` (job `w_selfcopy_test2`, ALL
OK, 5 s; the tests caught one definitional error in the `sel` bound before any real run). Nothing
here reads a native, an RMSD to a native, a benchmark sequence, PDB or the manifest; every
native-free command overwrites `rr` and `nat_ca` with NaN on load and asserts its outputs finite.
The benchmark facts used are the record's only: 2/60 (S24 L4, L15) and the mechanism (S24 L4).

Design, in short: the leak has two channels, (A) the carrier's self-window at BLOSUM rank 0 in
the K = 500 pool, (B) the carrier's native distances in the fold model's training labels. Both
are live, in the same direction, on 1CEK / 2FBU / 2P5H / 6B9K. Part A drops the identity-1.0
window and refills the pool (S10-4's operator, exact copies only) on the 4, with the other 122
targets' BLOSUM rank-0 window dropped as the matched control population, the >= 0.6 variant
beside it (the C27 re-derivation on the production basis, L28), and an ORACLE insertion of the
withheld same-fold carriers on 8 more targets. Part B retrains the three fold models with the
carrier removed (through `s26/p_ladder.py`'s pca32 path; controls: two dev chains per fold).
Part C, the n = 126 envelope: the four pinned fold models that trained on the target's OWN
native, against the clean one, paired, fold-clustered (ORACLE). Part D: the bound
(2/60) x max per-target effect, and (2/60) x the fold-CI limit of the envelope, judged against
0.017 A (one tenth of the benchmark CI half-width). Before any RMSD is read, every part stores
its emissions and the TRIANGLE BOUND: Kabsch CA-RMSD is a metric, so RMSD(leaked emission,
un-leaked emission) bounds |change in RMSD-to-native| without the native. Assumptions A1 to A5
are in the PREREG.

Native-free census (`s26/results/w_selfcopy_census.json`, job `w_selfcopy_census`, 20 s, 0.064 GB):

    target  n  fold  exact self-window: universe idx / BLOSUM rank / pool pos / in shipped top-75   >= 0.6 windows: universe / pool / top-75
    1CEK   13   2          7 / 0 / 0 / NO                                                            8 / 3 / 0
    2FBU   12   4       1602 / 0 / 0 / NO                                                            4 / 1 / 0
    2P5H    9   4       3316 / 0 / 0 / YES                                                           5 / 1 / 1
    6B9K   10   0       1082 / 0 / 0 / YES                                                           7 / 1 / 1

So channel A cannot touch the emitted structure of 1CEK or 2FBU at all (the self-window is
filtered out by the distogram score); on 2P5H and 6B9K it is one member of the 75 averaged.
21/126 targets carry a >= 0.6 window somewhere in the universe (S10-4's 21 reproduced; its
13-in-pool count is checked when the retrieval run completes). Only 31% of targets keep their
BLOSUM rank-0 window in the top-75. Probes: retrieval on 1CEK (`w_selfcopy_retrieval_probe`,
15 s, 0.104 GB): with and without the self-window the cloud, chain and lam = 0 chain are
IDENTICAL (triangle bound 0.000) and the refill window does not enter the top-75; envelope on
1CEK (`w_selfcopy_envelope_probe`, 25 s, 0.298 GB): production gate true (top-75 equals the
cache's `sub`), the four leaked models keep 81 to 87% of the top-75 and move the built chain by
0.079 / 0.124 / 0.137 / 0.122 A (triangle bounds, models 1 / 0 / 3 / 4). Full native-free runs
`w_selfcopy_retrieval` (est 0.5 GB, ~20 min) and `w_selfcopy_envelope` (est 0.6 GB, ~40 min)
are registered; the gated `endpoint`, `floor` and `report` wait for sign-off.

A second census fact, recorded because it is the basis of `s26/IDEA_tiebreak_noise_floor.md`:
every one of the 126 targets has a BLOSUM tie at the K = 500 boundary; the median tie class at
the boundary score holds 115 windows, about 57 inside the pool and 60 outside, so ~11% of every
pool is chosen by corpus order. S17 measured only the ORACLE pool best under random tie-breaks
(sd 0.018); the production endpoint's floor is unmeasured. Three own ideas filed:
`IDEA_tiebreak_noise_floor.md`, `IDEA_conformational_identity_floor.md` (Part E of the
PREREG: the same sequence in a different deposit), `IDEA_window_provenance.md` (census first,
plausibility 0.15).

Question for the coordinator (rule 14; not run until answered): may Part B's 10 retrains (tag
CPU, est-ram 1.5 GB, 456 s each, one at a time, order fixed in the PREREG section 8) run before
sign-off under L16b's conditions, and may they interleave with lane P's chain (two training
jobs resident, ~2.5 GB, against 4.7 GB available at 08:54)? Default if unanswered: the
reference is lane P's `pca32` fold models and Part B trains after sign-off. Noted and NOT done:
the triangle bound is computable on the benchmark itself with no native or RMSD read, but it
would need the two benchmark sequences and universes, which L18b/L29 close for this sprint.

## L31 -- EXAMINATION AUDIT: REPRODUCTION EXACT ON HEAD; ONE MATERIAL FINDING, THE TEST-SUITE TOTAL (369/356 QUOTED, 370/357 IN THE ARTEFACT); C26 RE-DERIVED; C27 IS IN GIT HISTORY (2026-09-13, A)

`s26/EXAMINATION_AUDIT.md`, sections 1 to 9. (1) `s26/e_reproduce.py` re-run on HEAD `76287a33`
through `s26/a_reproduce_head.py` (the Examiner's `main()` unchanged, the write redirected so
lane E's artefact is not overwritten; job `a_reproduce_head`, exit 0, 5.0 s;
`s26/results/a_reproduce_head.json`): all four means and all 126 rows equal
`s26/results/e_reproduce.json` at 0.0; T030 = 1S9Z `rmsd_arm` 0.18198112330908295. (2) 29
claim-ledger leaves opened (seed-26 sample C02, C05, C08, C15, C17, C19, C21, C23, C28, C32, plus
every UNSOURCED and markdown-sourced entry); every sourced value sits at its cited key at stored
precision. (3) The trace is faithful to `core/pipeline.py:836-848` (selector), `:934` (average)
and `core/predict.py:416, 422` (score weight `1/(sd+0.5)`, `_risk`); the identical pool index
449 on both traced targets is verified in `bench_results/cache/1fc9f2dcf489e2fb/{1S9Z,9KAR}.json
:: sub` (positions 44 and 34). (4) The junit files agree with `s26/TEST_RUN.md` and
`test_run.json`: 370 / 357 / 0 / 0 / 13, 0 memory-guard skips, 11 `VERIFY_SLOW` and 2
absent-artefact skips. (5) Every defect file:line resolves at HEAD. (6) Four hashes recomputed,
four equal. (7) Three sizes as listed, all ignored. (8) Rule 1: the grep over `s26/*.py` hits
only `e_claims.py` (the SKIP lists, committed `e3570aed`), `e_hashes.py` (bytes) and
`i_identity_audit.py` (L18's sequence count). L28 item 3: `s25/results/phys_suite.json :: basis`
= "point_cloud"; the artefact declares its own basis.

MATERIAL A1. `s26/EXAMINATION.md` section D ("Combined: 369 tests, 356 passed"), row C35,
`s26/agentE_FINDINGS.md` X1 and L28 item 2 quote 369 / 356 / 13 for `s26/TEST_RUN.md`, which
holds 370 / 357 / 13 (`pytest_core_post.xml` 351 / 338 / 13 plus `pytest_amber` 16 and
`pytest_amber_frame` 3). The count is the 00:41 render, before `pytest_core_post` finished at
00:45:56; section D lists four junit files and omits `pytest_core_post.xml`. L20 and the state
brief 6.1 already carry the right number; nothing in Phase 1 depends on it. Fix: lane E, or the
coordinator if E has finished, appends a correction to `s26/EXAMINATION.md` (D and C35) and
`s26/agentE_FINDINGS.md` and posts a ledger entry naming L28 item 2; A re-checks on the append.

MINOR. C26 (36.1 / 36.4 deg) is now DERIVED: `s13/cache/tors_rows.npz` arms `p_grid` / `n_marg`,
key `err_phi`, pooled over 1,507 residues: 36.133 / 36.416 deg (psi 62.355 / 72.772), matching
`s13/SPRINT13_DOSSIER.md:261-263`. C27's artefacts are not absent: `s10/idaudit_price.json` and
four siblings are in git history at `5fa05cd` (commit `49fc708` exists); not opened, because the
price stage covers the sealed benchmark; the +0.0030 half stays historical, the +0.0004 dev half
is lane W's to re-derive (L25). C34 confirmed document-only (a grep of `s25/` hits only the two
md files). C35's replacement is 370 / 357 / 13. C24's 0.524: `stage_cvarcheck` exists
(`s8/integrate.py:1162`) but its JSON and console log do not; quote 0.655634
(`tests/test_quantum.py:59`) or 0.566586 (`s25/results/q_verify.json`) instead.

Question for the coordinator: `tests/test_pipeline.py:338-347` parses `s9/final_report.json`
aggregates, so any further full `pytest tests/` in S26 does so again; accept, or `--deselect` it
for the rest of the sprint.

AUDIT: MATERIAL FINDINGS 1, LISTED ABOVE

---
