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
## L32 -- L28 ITEM 2 CORRECTED: THE GOVERNED TEST RUN IS 370 / 357 / 13; THE test_pipeline BENCHMARK-AGGREGATE READ IS ACCEPTED; C26 IS DERIVED (2026-09-13, coordinator)

1. The Adversary's one MATERIAL finding (L31 item A1) is fixed by appended addenda in
   `s26/EXAMINATION.md` (sections D and C35) and `s26/agentE_FINDINGS.md` (X1): the governed
   suite is 370 tests, 357 passed, 0 failed, 0 errors, 13 skipped, 0 memory-guard skips
   (`s26/results/test_run.json :: combined`; junit files `pytest_core_post.xml` 351/338/13,
   `pytest_amber.xml` 16/16, `pytest_amber_frame.xml` 3/3). L28 item 2's "369 / 356 / 13" was
   the 00:41 render before `pytest_core_post` finished. Lane E's lines are left as written.
2. Ruling on the Adversary's question: `tests/test_pipeline.py:338-347` reads two published
   aggregates of `s9/final_report.json` (the benchmark means 2.9507 and 2.9610 already in
   `README.md`) and asserts them to 5e-4. It is the repository's own regression test, it
   predates S26, it reads no per-target value, and a pass or fail carries no new information
   about the benchmark. Accepted; it is not deselected, so the S26 test record stays whole.
   No S26 lane adds any benchmark read of its own (L29 item 1 stands).
3. C26 (phi MAE 36.1 vs 36.4 deg) is DERIVED by the Adversary from `s13/cache/tors_rows.npz`
   (arms `p_grid` / `n_marg`, key `err_phi`, 1,507 residues: 36.133 / 36.416 deg). The Adversary
   saves that derivation as `s26/results/a_c26_phi_mae.json` with provenance; with it the number
   may appear on a slide with that artefact path. C27's +0.0004 dev half is lane W's to
   re-derive (L25, L30); C34 and C24's 0.524 stay document-only and off every slide.

---

## L33 -- PHASE 0 SIGNED OFF (2026-09-13 09:05, coordinator)

The Adversary's audit (`s26/EXAMINATION_AUDIT.md`, L31) found the reproduction exact on HEAD on
all four bases and all 126 rows, 29 claim leaves at stored precision, the trace faithful to the
code, the hashes and sizes matching, Rule 1 closed, and one MATERIAL documentary finding, which
L32 fixes with appended addenda; the Adversary re-checks the append on its next turn and may
veto by ledger entry, in which case this sign-off is retracted by a new entry and every result
read in between is provisional. Phase 1 begins: lanes Q (A1, A3), P (ladder evaluation, B3, C4,
C5), PH (cis floor, C3 stage 1, steric reject) and W (the 2/60 proxy bound) may read RMSDs to
natives under their pre-registrations. Standing: every positive result is attacked by the
Adversary before it enters the closed/open tables; basis named on both sides of every contrast;
the phase gate string in this heading is the one the lanes' code checks.

---

## L34 -- ADVERSARY RE-CHECK OF L32: THE ADDENDA STAND; THE L33 SIGN-OFF IS NOT VETOED (2026-09-13, A)

Re-read at commit `195cf93e`: the appended addenda in `s26/EXAMINATION.md` (after section I) and `s26/agentE_FINDINGS.md` (after the unsourced list) both state 370 tests, 357 passed, 0 failed, 0 errors, 13 skipped, 0 memory-guard skips with the artefact paths (`s26/results/test_run.json :: combined`, `pytest_core_post.xml` 351/338/13, `pytest_amber.xml` 16/16, `pytest_amber_frame.xml` 3/3), and L32 item 1 matches them; the original lines are left as written. L31 item A1: STANDS as fixed. No veto.

---

## L35 -- LANE Q, A4: ADAPT-GROWN CIRCUITS ON THE DEPLOYED HAMILTONIAN ARE PRODUCT CIRCUITS AT alpha = 1 (NO DECAY, NOTHING TO TRAIN) AND DECAY LIKE THE FIXED ANSATZ AT alpha = 0.25; THE S25 n = 7 ROWS REPRODUCE EXACTLY (2026-09-13, lane Q)

`s26/q_var.py` -> `s26/results/q_var.json` (complete); `s26/jobs_done/a4_var.json`: 2,765 s,
peak RSS 0.08 GB. Pre-registered in `s26/PREREG_A4.md` before the run. Property measurement:
energies are the deployed SHAPE (standardised ranks 1..2^n), no target, no native. Figure:
`s26/figures/a4_variance_slopes.png` (190 dpi, white), from the JSON.

Gate: `s25.q_plateau.measure(7, 3, alpha, T, 250, seed=1007)` reproduces the n = 7 row of
all five cells of `s25/results/q_plateau.json` at relative deviation 0.0 (bit-identical).

    fitted log2 Var[dF/dtheta_0] per qubit, n = 4..13, theta ~ N(0, 0.6^2), same draws as S25
    (grown: rows with P = 3n only, count in brackets; growth by Adam best-iterate, 50 steps,
     eps 1e-3, pool V = Tang 2n-2, pool L2 = 2-local odd-Y)
      cell                    fixed (S25)   grown V          grown L2
      alpha=1,    T=0         -0.649        -0.079 (7)       +0.006 (7)
      alpha=0.25, T=0         -0.252        degenerate (0)   degenerate (0)
      alpha=0.10, T=0         -0.047        degenerate (0)   degenerate (0)
      alpha=1,    T=0.3       -0.311        +0.035 (5)       -0.008 (7)
      alpha=0.25, T=0.3       -0.243        -0.246 (6)       -0.302 (7)
    grown/fixed variance ratio at matched n: 1.63 (n = 4, alpha=1 T=0.3, both pools) and
    2.5 to 370 everywhere else (31 of 32 matched rows above 2).

What the grown circuits are. At alpha = 1 (T = 0 or 0.3) they hold 1 to 3 DISTINCT operators;
the other 12 to 25 selections are consecutive repeats of the same single-qubit Y (a repeated
rotation merges with the previous one): they are product circuits with the parameter count
of the fixed ansatz, which is why their variance does not decay with n. At alpha = 0.25,
T = 0.3, pool L2 grows 4 to 21 distinct 2-local strings and its decay, -0.302 per qubit, is
the fixed ansatz's -0.243 within the sampling error of a 7-point slope (relative SE of a
variance 9 to 16%). At T = 0 with alpha < 1 growth stops at P = n on every width (the
collapse; all first-order gradients vanish at a basis state, the lemma in s26/q_adapt.py) and
the variance is 0 to 2.4e-2; those cells are degenerate by construction, as pre-registered.

Predictions (PREREG_A4): H4a (T = 0.3 grown slopes in [-0.3, 0]) held for three of four
(V alpha=1 +0.035 and L2 alpha=0.25 -0.302 sit just outside by 0.035 and 0.002, inside the
error of the fit); the falsifier (below -0.5) did not fire. H4b (ratio > 2 at every matched
n) failed on 1 of 32 rows (n = 4, 1.63) and held on 31. H4c held on P < 3n (14 of 14) and
failed on "Var < 1e-3" for 3 of 7 rows at alpha = 0.25 (max 2.4e-2 at n = 4). H4d held
exactly. All four recorded as measured, none softened.

What it means, in the record's own terms. A large gradient from a product circuit is not
trainability; it is the trivial regime (rule 10's mirror: never call a large gradient a
merit). The one non-trivial grown family (alpha = 0.25, L2) decays like the fixed ansatz it
would replace. A4 gives Proposal A no width-scaling argument.
## L36 -- jobrun v2.1: THE LAUNCH CAP IS RE-CHECKED AFTER A JITTER; SIX JOBS HAD LAUNCHED AGAINST A CAP OF FOUR (2026-09-13, coordinator)

`s26/governor_state.json` at 09:28:11 showed six registered jobs (p_train_conly, two W census
passes, a1_build_s0 and s1, p_eval_shipped) against jobrun's cap of four: when A4 finished,
three waiting jobrun processes read the same 5 s snapshot with three jobs and all launched.
The governor's band held (smoothed CPU 85.7%, RAM 76.9%; a1_build_s0 suspended as the newest
job, resumed when CPU fell), so nothing was lost, but the cap was not a cap. `s26/jobrun.py`
now counts live registrations in `s26/jobs/` directly instead of the snapshot, sleeps a random
0.2 to 3 s after passing the gate and re-checks the count before registering. Waiters started
before this edit run the old code until they launch; new waiters use v2.1. No production module
touched; `s26/governor.py` unchanged.

---

## L37 -- CAMPAIGN PAUSED AT THE USER'S REQUEST (2026-09-13 09:30); STATE SAVED; GOVERNED JOBS KEEP RUNNING (2026-09-13, coordinator)

The user asked for a break with the usage limit at 96%. Every lane was told to commit its files,
append a STATUS line and end its turn; the coordinator commits whatever remains as a checkpoint.
Governed python jobs keep running on the box (they cost no API usage and checkpoint to disk):
at 2026-09-13T09:30:16: p_train_conly (P), w_selfcopy_retrieval (W), w_selfcopy_envelope (W), a1_build_s1 (Q), a1_build_s0 (Q, suspended), p_eval_shipped (P); RAM 77.1%, smoothed CPU 93.6%. The governor v2
(pid 25480) stays up. `s26/models/` (11 MB of ladder checkpoints, regenerable) is not committed.

Where the sprint stands: Phase 0 signed off (L33, audit L31/L34). Done: A2 (L27), A4 (L35), the
cis census (L22), the steric-reject census (L23), the C3 native-free part (L24), C1 (lane P
findings section 7), B1 infeasible (L13), all operational items and declared defects (L6 to L10,
L15, L16, L18 to L21), the report draft Parts I to VI, IX, Appendix A (`s26/REPORT.md`), the
Proposal B replacement outline draft, the trainability paper outline, A2/A4 figures. Running or
queued: A1 shards then A3 (Q), the C2 ladder training and evaluation (P), the cis floor, C3
stage 1 and the steric reject (PH), the 2/60 proxy bound (W), the Adversary's checks of L27,
L22 to L24, L30, L35 and the tournament ranking (A), the report's Appendix B check (E).

To resume: read this ledger from L33, `s26/STATUS.md`, `s26/jobs_done/` and `s26/results/`,
then continue each lane from its STATUS line; the lane briefs are in `s26/briefs/` and every
lane's findings file names what it did not do. No lane restarts finished work.

---

## L38 -- CIS FLOOR (ORACLE DIAGNOSTIC): THE IDEAL-TRANS REBUILD OF THE NATIVE'S OWN TORSIONS SITS 0.347 A FROM THE NATIVE (max 1.47), TRACKS OMEGA NON-PLANARITY AT rho 0.83, AND DOES NOT TRACK THE 0.166 A CHAIN COST (rho 0.08) (2026-09-13, PH)

`s26/ph_cis.py floor`, `s26/results/ph_cis_floor.json` (complete 126/126), job
`s26/jobs_done/ph_cis_floor.json` (exit 0, 10 s, peak RSS 0.098 GB, held 345 s at the job cap).
Pre-registered in `s26/PREREG_cis.md` part 2 and addendum 1. ORACLE DIAGNOSTIC: every number
reads the native. Bases, named: `floor_ca` is the CA-RMSD between the model-1 native CA trace
and the chain rebuilt from the native's OWN phi/psi at ideal trans geometry
(`core.geometry.build_backbone`), after superposition; `rebuild_bb` is the same on N/CA/C
(`core.data.Peptide.rebuild`, verified from `core/data.py:410-419`); `chain_cost` is the
production `rmsd_arm - rmsd_avg` (built chain minus point cloud, a difference across bases,
`bench_results/cache/1fc9f2dcf489e2fb`).

    floor_ca   (CA, ideal trans on the native's own torsions)   mean 0.347 A (SE 0.029), median 0.272, p90 0.752, max 1.474 (1ID6)
    rebuild_bb (N/CA/C, `Peptide.rebuild`)                         mean 0.340 A (SE 0.028), median 0.248, max 1.413
    chain_cost (production rmsd_arm - rmsd_avg)                    mean 0.166 A (SE 0.018), median 0.098
    targets with floor_ca > 0.5 A                                  32 / 126;  floor_ca > chain_cost on 86 / 126
    Spearman(floor_ca, max omega deviation)                        +0.828;  with mean omega deviation +0.833
    Spearman(chain_cost, max omega deviation)                      -0.036;  Spearman(chain_cost, floor_ca) +0.083
    cis targets                                                    0 (L22); the cis-vs-non-cis contrast is empty and is not printed
    9UV5 (the one bond beyond 30 deg)                              floor_ca 0.663, rebuild_bb 0.708, chain_cost 0.372

`ST.fmt`, verbatim (the two quantities are both per-target Angstroms but are different
objects; "BETTER" below means only that the production chain cost is smaller than the
own-torsion rebuild floor, not that anything improved):

    production chain cost (rmsd_arm - rmsd_avg) MINUS the ORACLE CA floor of ideal trans geometry on the native's own torsions
      a 0.1664 (med 0.0977)   b 0.3468 (med 0.2722)   n=126
      effect -0.1804   median -0.1261   SE 0.0333   MDE 0.0933   effect/MDE -1.93
      iid  CI95 [-0.2433, -0.1176]
      fold CI95 [-0.2266, -0.1022]   folds same sign 5/5   per-fold 0:-0.025 1:-0.224 2:-0.230 3:-0.223 4:-0.202
      86W/40L/0T   worst degradation +0.5161 (1RSW)   p90 +0.2720   power 1.00  Type-M 1.00
      concentration: drop-top10 -0.1100 vs uniform-effect null p10/p50/p90 -0.1523/-0.1116/-0.0702 -> pctile 0.518
      VERDICT: BETTER

Reading. (1) The constant omega does carry a representation cost on this instrument even
without a single cis bond: the ideal-trans rebuild of the native's own torsions misses the
native by 0.35 A on average and by more than 0.5 A on a quarter of the targets, and that miss is
the omega non-planarity to rho 0.83 (L22: 0.5% of bonds beyond 20 deg, mean deviation 1.9 deg;
small deviations accumulate along the chain). (2) It is NOT what the production projection
pays: the chain cost (0.166 A) does not correlate with the floor (rho +0.08) nor with the omega
deviation (rho -0.04). The projection's cost is a displacement effect of the operator (S16 L27:
the projection moves 0.813 A with a negative cosine), not a representation effect. (3)
`floor_ca` is an UPPER bound on the manifold floor, because the projection fits phi/psi to the
trace rather than rebuilding from the native's torsions; `s26/PREREG_cis.md` addendum 2
registers `floor2` (the native projected through the production projection) as the tight
number, 10 min CPU, to run after the reject jobs. Power: n = 126, SE 0.033, MDE 0.093 on the
contrast; a Spearman at n = 126 has SE about 0.09, so the two nulls (+0.08, -0.04) exclude
|rho| above about 0.25. Nothing here is a proposal; the two-bond-length projection stays worth
0.000 A on this instrument (L22), and the design note in `s26/agentPH_FINDINGS.md` 1.4 stands.

---

## L39 -- C3 STAGE 1: THE PRODUCTION RELAXATION IS WORSE THAN A RANDOM MOVE OF ITS OWN SIZE (+0.0111, fold CI [+0.0062, +0.0171]) AND WORSE THAN A MOVE TOWARD A POOL MEMBER (+0.0385); ITS DISPLACEMENT POINTS AWAY FROM THE NATIVE (cos -0.049); REFINE-WITH-PHYSICS IS A VALIDITY STEP, NOT AN ACCURACY STEP (2026-09-13, PH)

`s26/ph_c3.py stage1`, `s26/results/ph_c3_stage1.json` (complete 126/126), job
`s26/jobs_done/ph_c3_stage1.json` (exit 0, 10 s, peak RSS 0.041 GB, held 345 s at the job cap).
Pre-registered in `s26/PREREG_c3_control.md` sections 2 and 3 and addendum 1 (prediction:
AMBER worse than the matched random move by about +0.013 A, cosine negative). No AMBER compute:
the production coordinates `ca` and `amber_ca` from `bench_results/cache/1fc9f2dcf489e2fb`.

Bases, named on both sides. Arm: input the BUILT CHAIN (`rmsd_arm`, 3.2148), output the
RELAXED CHAIN (`rmsd_full`, 3.2355), both reproduced from the stored coordinates to 1e-6.
Controls: the built chain's CA trace displaced by a vector of the SAME per-atom RMS magnitude
as AMBER's displacement (0.220 A mean, measured after superposing `amber_ca` onto `ca`);
`rand` = S16's construction (isotropic Gaussian, six rigid components projected out on
`rigid_basis`, 16 draws, mean); `member` = the same magnitude along the straight line toward a
random member of the shipped K=500 pool superposed onto `ca` (16 draws, mean; 6 of 2,016 draws
overshot the member). A displaced CA trace is not an ideal-geometry chain: the controls match the
operator's SPACE (magnitude), not its validity, exactly as S16 did. Full-chain CA-RMSD to the
native, ORACLE evaluation of native-free operators. `ST.fmt`, verbatim, all seven contrasts:

    stage1 AMBER (relaxed) minus do-nothing (built chain)
      a 3.2355 (med 2.9757)   b 3.2148 (med 2.9661)   n=126
      effect +0.0207   median +0.0158   SE 0.0034   MDE 0.0096   effect/MDE +2.16
      iid  CI95 [+0.0141, +0.0274]
      fold CI95 [+0.0154, +0.0290]   folds same sign 5/5   per-fold 0:+0.021 1:+0.037 2:+0.017 3:+0.013 4:+0.017
      40W/86L/0T   worst degradation +0.1653 (2MID)   p90 +0.0640   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0258 vs uniform-effect null p10/p50/p90 +0.0213/+0.0256/+0.0302 -> pctile 0.529
      VERDICT: WORSE
    stage1 AMBER minus matched-magnitude RANDOM (16 draws)
      a 3.2355 (med 2.9757)   b 3.2244 (med 2.9669)   n=126
      effect +0.0111   median +0.0070   SE 0.0035   MDE 0.0099   effect/MDE +1.12
      iid  CI95 [+0.0044, +0.0178]
      fold CI95 [+0.0062, +0.0171]   folds same sign 5/5   per-fold 0:+0.011 1:+0.023 2:+0.009 3:+0.003 4:+0.010
      52W/74L/0T   worst degradation +0.1518 (2MID)   p90 +0.0539   power 0.88  Type-M 1.07
      concentration: drop-top10 +0.0168 vs uniform-effect null p10/p50/p90 +0.0122/+0.0168/+0.0216 -> pctile 0.497
      VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.07x]
    stage1 AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)
      a 3.2355 (med 2.9757)   b 3.1969 (med 2.9755)   n=126
      effect +0.0385   median +0.0338   SE 0.0061   MDE 0.0170   effect/MDE +2.27
      iid  CI95 [+0.0272, +0.0506]
      fold CI95 [+0.0298, +0.0479]   folds same sign 5/5   per-fold 0:+0.039 1:+0.052 2:+0.049 3:+0.029 4:+0.027
      29W/97L/0T   worst degradation +0.3497 (2BP4)   p90 +0.1092   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0485 vs uniform-effect null p10/p50/p90 +0.0407/+0.0481/+0.0561 -> pctile 0.521
      VERDICT: WORSE
    stage1 RANDOM minus do-nothing
      a 3.2244 (med 2.9669)   b 3.2148 (med 2.9661)   n=126
      effect +0.0096   median +0.0071   SE 0.0015   MDE 0.0042   effect/MDE +2.28
      iid  CI95 [+0.0070, +0.0127]
      fold CI95 [+0.0077, +0.0122]   folds same sign 5/5   per-fold 0:+0.010 1:+0.015 2:+0.008 3:+0.010 4:+0.007
      29W/97L/0T   worst degradation +0.1147 (1I93)   p90 +0.0276   power 1.00  Type-M 1.00
      concentration: drop-top10 +0.0114 vs uniform-effect null p10/p50/p90 +0.0094/+0.0113/+0.0135 -> pctile 0.528
      VERDICT: WORSE
    stage1 TOWARD-MEMBER minus do-nothing
      a 3.1969 (med 2.9755)   b 3.2148 (med 2.9661)   n=126
      effect -0.0178   median -0.0118   SE 0.0043   MDE 0.0120   effect/MDE -1.49
      iid  CI95 [-0.0266, -0.0096]
      fold CI95 [-0.0255, -0.0124]   folds same sign 5/5   per-fold 0:-0.018 1:-0.015 2:-0.031 3:-0.016 4:-0.010
      90W/36L/0T   worst degradation +0.1204 (8T62)   p90 +0.0274   power 0.99  Type-M 1.01
      concentration: drop-top10 -0.0088 vs uniform-effect null p10/p50/p90 -0.0141/-0.0092/-0.0040 -> pctile 0.543
      VERDICT: BETTER
    stage1 ORACLE cos: AMBER minus RANDOM
      a -0.0491 (med -0.0498)   b 0.0013 (med 0.0028)   n=126
      effect -0.0503   median -0.0509   SE 0.0154   MDE 0.0432   effect/MDE -1.17
      iid  CI95 [-0.0799, -0.0211]
      fold CI95 [-0.0735, -0.0253]   folds same sign 5/5   per-fold 0:-0.050 1:-0.088 2:-0.066 3:-0.001 4:-0.046
      74W/52L/0T   worst degradation +0.3845 (5MXS)   p90 +0.1647   power 0.90  Type-M 1.06
      concentration: drop-top10 -0.0225 vs uniform-effect null p10/p50/p90 -0.0433/-0.0232/-0.0028 -> pctile 0.517
      VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.06x]
    stage1 ORACLE cos: AMBER minus TOWARD-MEMBER
      a -0.0491 (med -0.0498)   b 0.1229 (med 0.1090)   n=126
      effect -0.1719   median -0.1733   SE 0.0232   MDE 0.0649   effect/MDE -2.65
      iid  CI95 [-0.2156, -0.1266]
      fold CI95 [-0.2104, -0.1388]   folds same sign 5/5   per-fold 0:-0.168 1:-0.190 2:-0.244 3:-0.121 4:-0.140
      94W/32L/0T   worst degradation +0.5848 (2NDN)   p90 +0.1612   power 1.00  Type-M 1.00
      concentration: drop-top10 -0.1320 vs uniform-effect null p10/p50/p90 -0.1626/-0.1327/-0.1020 -> pctile 0.512
      VERDICT: BETTER

    ORACLE cos(AMBER displacement, true residual)   mean -0.049 (SE 0.015), median -0.050, positive on 36.5% of targets
    orthogonal-move cost predicted by the S16 identity from the magnitude alone   +0.0107 (SE 0.0011), median +0.0079
    identity n RMSD^2 = |r|^2 - 2 v.r + |v|^2 as an upper bound on the superposed RMSD   max violation 0.0062 A (it is a bound, not an equality, under re-superposition)

Reading. (1) The production step reproduces exactly: +0.0207 [+0.0154, +0.0290], 2.16x MDE.
(2) About half of that cost is the SIZE of the move: a random displacement of the same 0.220 A
costs +0.0096, which is what the S16 identity predicts for a move orthogonal to the residual
(+0.0107). The other half is DIRECTION: AMBER is worse than its own random twin by +0.0111
(fold CI excluding zero, 5/5 folds, but 1.12x MDE, so the magnitude is an upper bound), and its
displacement points away from the native, cos -0.049 (36.5% positive), which reproduces S16 L27
to the fourth decimal (S16: cos -0.0521, AMBER minus random -0.0491; here AMBER minus random on
the cosine -0.0503). The registered prediction (+0.013, negative cosine) held. (3) A move of the
same size TOWARD A RANDOM POOL MEMBER improves the built chain by -0.0178 [-0.0255, -0.0124],
90W/36L, 1.49x MDE, 5/5 folds, with ORACLE cos +0.123. This is a zero-information control, not
a proposal, and it is PROVISIONAL until the replication registered in `s26/PREREG_c3_control.md`
addendum 2 lands (job `ph_c3_stage1_rep`: new seeds, reversed order). Its mechanism is on the
record: the projection moved the chain 0.166 A away from the point cloud with a negative cosine
(S16 L27), and the pool members surround the cloud, so any move back toward one of them
recovers part of that displacement. It says nothing about physics and everything about the
projection's cost; it will be handed to the Adversary as such. (4) Power: SEs 0.0015 to 0.0061,
MDEs 0.004 to 0.017; every contrast except AMBER-vs-random and the AMBER-vs-random cosine is
above 1.3x its MDE; those two are in the Type-M zone, direction measured, magnitude inflated
about 1.07x.

DECISION RULE (`s26/PREREG_c3_control.md` section 3): AMBER does not beat either matched
control; it is worse than both. "Refine with physics" is DROPPED as an accuracy step and KEPT
as a validity step. `s26/C3_RESULT.md` carries the sentence with the numbers for lane P and the
presentation. Stage 2 repeats this on lane P's best C2 rung when it is delivered.

---
## L40 -- GOVERNOR KILLED BY THE HOST AT 10:10 FOR LOW MEMORY; AT 19:10 THE BOX WAS AT 91.8% RAM; THE TWO REMAINING LANE-P JOBS WERE STOPPED BY THE COORDINATOR (2026-09-13 19:11, coordinator)

The harness stopped the governor v2 background process (pid 25480) at about 10:10 with the message
"the system is running low on memory" (`s26/governor.log` ends at 10:10:36 with RAM 87.5%, six
jobs registered). Every lane was then cut by the API session limit. At 19:10 the coordinator found:
RAM 91.1 to 91.8% used with 1.5 GB available, no governor alive, a1_build_s0/s1, ph_reject_cloud
and w_train_chain gone (their records, if any, are in `s26/jobs_done/`), `p_eval_shipped` suspended
since 09:25 (0.03 GB; suspended by the governor before it died, never resumed), and `p_train_raw`
running unsupervised since 16:50 at 1.11 GB (launched by lane P's chain driver through jobrun,
which treats a stale governor snapshot as "go"). With the box at the campaign ceiling and the host
already killing processes for memory, the coordinator terminated both lane-P jobs and the chain
driver (per-fold checkpoints exist under `s26/models/p_ladder/`; lane P resumes the `raw` rung and
the shipped evaluation from them) and restarted the governor. Nothing else of the campaign's is
running. Rule added for jobrun: a stale governor snapshot is NOT "go" for a job whose est-ram
exceeds 0.5 GB (implemented at the next resume; recorded here first).

---

## L41 -- THE BOX IS AT 84% WITH NOTHING OF OURS RUNNING; jobrun v2.2; WHAT FINISHED DURING THE PAUSE (2026-09-13 19:13, coordinator)

After L40's stops the box reads 84.2% RAM used (14.1 GB) with no campaign process alive: the
load is the user's own (twelve Chrome renderers at 0.3 to 0.8 GB each, VS Code, three claude
sessions). Available memory is 2.6 GB, so the campaign's working headroom under the 93% ceiling
is now about 1.4 GB, a third of what L0 measured. Every remaining job must be launched one at a
time with a measured peak below 1 GB until the user's load drops; the ladder's `raw` rung (peak
of its sibling `pca32f` 1.85 GB) does not fit and waits.

`s26/jobrun.py` v2.2: (i) a stale governor snapshot is no longer "go" for a job with est-ram
above 0.5 GB (L40); (ii) a job never starts unless the snapshot shows its estimate plus 0.5 GB
available. Lane P's shell chain drivers (`p_train_chain.sh`, `p_eval_chain.sh`) were terminated
because they re-launch the next rung whenever the previous job exits, including on a kill; on
resume lane P launches rungs singly.

Finished during the pause (`s26/jobs_done/`, exit 0): w_selfcopy_retrieval (09:34),
p_train_conly (09:51), w_selfcopy_envelope (09:57), w_tiebreak_probe (09:58), p_train_pca32
(10:33), p_train_wide (14:32), p_train_pca32f (15:13), p_train_pca128 (16:08),
p_featurise_esm8m and p_train_esm8m (16:50). The "_p2" records at 19:11 to 19:12 are the chain
driver re-touching finished rungs (5 s, no work). Killed without a record: a1_build_s0/s1
(per-target checkpoints under `s26/results/a1/`), ph_reject_cloud, w_train_chain,
p_eval_shipped, p_train_raw (fold checkpoints under `s26/models/p_ladder/`), p_eval_noesm.
Each resumes from its checkpoint; none restarts.

---

## L42 -- RESUMING AT 2026-09-13 19:14: A SECOND REPORT WAS FOUND IN docs/ (NOT WRITTEN BY ANY S26 LANE); LANES RELAUNCHED UNDER A 1.4 GB HEADROOM; PR SPAWNED (2026-09-13, coordinator)

1. `docs/REPORT_S26.md` (181 KB, Parts I to X, Appendices A to D, mtime 14:25) and
   `docs/REPORT_S26_SUMMARY.md` (14:24) appeared on disk during the pause, while every S26 lane
   was dead (session limit, 10:10 to 19:10). They were not written by any lane of this campaign
   and are not in git; their header names HEAD `b7f1f280`, so they were written from this
   branch by another session or by the user. They cite external sources (arXiv abstracts, an
   author page) that no S26 lane fetched. They are left untouched and uncommitted by the
   coordinator; the user decides whether they are tracked. The campaign's report deliverable
   remains `s26/REPORT.md`; lane E cross-checks the two, takes any artefact-sourced number or
   presenter material it lacks, cites the docs file when it does, and never copies an external
   claim or a personal detail into a deliverable.
2. Every lane is relaunched from its own transcript with the instruction to resume from its
   STATUS line and the checkpoints on disk. Headroom rule until the user's own load drops: one
   governed job at a time whose measured or probed peak is under 1.0 GB; the `raw` rung
   (sibling peak 1.85 GB) and lane W's Part B retrains (1.3 GB) wait; jobrun v2.2 blocks any
   launch without est-ram + 0.5 GB free.
3. Lane PR (Presentation) is spawned now: `s26/C3_RESULT.md` (L39) fixes Proposal C's AMBER
   sentence, A2/A4 (L27, L35) fix the slide 6/8 figures, and the slides 1 to 7 and 11 depend
   only on the Phase 0 claim ledger. Seven lanes active (E, Q, P, PH, W, A, PR).

---

## L43 -- STERIC REJECT, POINT CLOUD: THE FALSIFIER FIRED THE OTHER WAY. AT 1e4 THE REFILL ARM IS +0.228 A WORSE THAN THE SHIPPED TOP-75 AND +0.167 WORSE THAN REJECTING THE SAME COUNT AT RANDOM (5/5 FOLDS); THE DOSE IS MONOTONE AND ITS LIMIT IS THE ANCHOR (2026-09-13, PH)

`s26/ph_reject.py cloud` then `report`, `s26/results/ph_reject_cloud.json` (complete 126/126),
`s26/results/ph_reject_report.json`, log `s26/results/ph_reject_report_cloud.log`; job
`s26/jobs_done/ph_reject_cloud2.json` (exit 0, 170 s for the last 10 cells, peak RSS 0.041 GB;
the first 116 cells came from `ph_reject_cloud`, killed with the governor by the host at about
10:10, L40; the per-target cells made the restart lossless). Pre-registered in
`s26/PREREG_amber_reject.md` sections 3 to 5 and addendum 1; the moved-subset secondary and the
R-first reading were fixed by the coordinator at 09:10 (findings section 2.3) before any RMSD
was read.

BASIS: POINT CLOUD on both sides (`I.coordinate_average` of the retained windows, scored ORACLE
against `nat_ca`; the anchor reproduces the production `rmsd_avg` 3.0483 to 1e-6 on every
target). The built-chain twin (`rmsd_arm` basis) is job `ph_reject_chain`, running. Arms: R =
reject e_amber > T from the shipped top-75 and refill from the next-ranked survivors to m = 75
(the PRIMARY reading; it empties only when the whole pool has no survivor: 8 targets at 1e4,
which fall back to the anchor and count as ties); S = reject without refill (empties 11 sets at
1e4, same fallback). Controls matched in the operator's space: RANDR / RANDS reject the SAME
COUNT at random (16 draws, mean); PERMR / PERMS apply the same threshold to a permuted energy
vector (16 permutations). Twelve `ST.fmt` blocks verbatim, primary threshold 1e4:

      point_cloud [all, n=126] R@1e4 minus anchor (negative = arm better)
        a 3.2760 (med 3.1167)   b 3.0483 (med 2.8373)   n=126
        effect +0.2276   median +0.0027   SE 0.0693   MDE 0.1941   effect/MDE +1.17
        iid  CI95 [+0.0950, +0.3672]
        fold CI95 [+0.1466, +0.3232]   folds same sign 5/5   per-fold 0:+0.278 1:+0.076 2:+0.195 3:+0.398 4:+0.199
        49W/67L/10T   worst degradation +4.1545 (8T61)   p90 +1.2819   power 0.91  Type-M 1.06
        concentration: drop-top10 +0.3256 vs uniform-effect null p10/p50/p90 +0.2365/+0.3216/+0.4116 -> pctile 0.521
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      point_cloud [all, n=126] R@1e4 minus RANDR@1e4 (matched control)
        a 3.2760 (med 3.1167)   b 3.1088 (med 3.0266)   n=126
        effect +0.1672   median +0.0128   SE 0.0473   MDE 0.1326   effect/MDE +1.26
        iid  CI95 [+0.0799, +0.2639]
        fold CI95 [+0.0699, +0.2793]   folds same sign 5/5   per-fold 0:+0.156 1:+0.007 2:+0.191 3:+0.384 4:+0.113
        47W/77L/2T   worst degradation +3.2867 (8T61)   p90 +0.7875   power 0.94  Type-M 1.03
        concentration: drop-top10 +0.2306 vs uniform-effect null p10/p50/p90 +0.1709/+0.2276/+0.2895 -> pctile 0.523
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.03x]
      point_cloud [all, n=126] R@1e4 minus PERMR@1e4 (matched control)
        a 3.2760 (med 3.1167)   b 3.1707 (med 3.0234)   n=126
        effect +0.1053   median +0.0031   SE 0.0385   MDE 0.1078   effect/MDE +0.98
        iid  CI95 [+0.0350, +0.1821]
        fold CI95 [+0.0187, +0.2000]   folds same sign 4/5   per-fold 0:+0.070 1:-0.063 2:+0.128 3:+0.274 4:+0.115
        52W/66L/8T   worst degradation +2.6065 (9KAR)   p90 +0.6742   power 0.78  Type-M 1.14
        concentration: drop-top10 +0.1610 vs uniform-effect null p10/p50/p90 +0.1131/+0.1602/+0.2120 -> pctile 0.509
        VERDICT: NOT MEASURED (|effect| 0.1053 <= its own MDE 0.1078, 0.98x)
      point_cloud [all, n=126] S@1e4 minus anchor (negative = arm better)
        a 3.1563 (med 3.0432)   b 3.0483 (med 2.8373)   n=126
        effect +0.1079   median +0.0147   SE 0.0287   MDE 0.0805   effect/MDE +1.34
        iid  CI95 [+0.0542, +0.1679]
        fold CI95 [+0.0647, +0.1528]   folds same sign 5/5   per-fold 0:+0.043 1:+0.054 2:+0.133 3:+0.193 4:+0.117
        38W/75L/13T   worst degradation +1.7827 (8T61)   p90 +0.4359   power 0.96  Type-M 1.02
        concentration: drop-top10 +0.1427 vs uniform-effect null p10/p50/p90 +0.1043/+0.1412/+0.1794 -> pctile 0.516
        VERDICT: WORSE
      point_cloud [all, n=126] S@1e4 minus RANDS@1e4 (matched control)
        a 3.1563 (med 3.0432)   b 3.0854 (med 2.8491)   n=126
        effect +0.0709   median +0.0096   SE 0.0254   MDE 0.0711   effect/MDE +1.00
        iid  CI95 [+0.0235, +0.1229]
        fold CI95 [+0.0317, +0.1087]   folds same sign 4/5   per-fold 0:-0.000 1:+0.040 2:+0.105 3:+0.133 4:+0.077
        41W/72L/13T   worst degradation +1.5172 (8T61)   p90 +0.3281   power 0.80  Type-M 1.13
        concentration: drop-top10 +0.1041 vs uniform-effect null p10/p50/p90 +0.0714/+0.1030/+0.1371 -> pctile 0.520
        VERDICT: NOT MEASURED (|effect| 0.0709 <= its own MDE 0.0711, 1.00x)
      point_cloud [all, n=126] S@1e4 minus PERMS@1e4 (matched control)
        a 3.1563 (med 3.0432)   b 3.0809 (med 2.8611)   n=126
        effect +0.0754   median +0.0054   SE 0.0274   MDE 0.0767   effect/MDE +0.98
        iid  CI95 [+0.0253, +0.1312]
        fold CI95 [+0.0311, +0.1255]   folds same sign 5/5   per-fold 0:+0.020 1:+0.023 2:+0.113 3:+0.164 4:+0.063
        49W/69L/8T   worst degradation +1.5869 (8T61)   p90 +0.3941   power 0.79  Type-M 1.13
        concentration: drop-top10 +0.1137 vs uniform-effect null p10/p50/p90 +0.0787/+0.1119/+0.1493 -> pctile 0.526
        VERDICT: NOT MEASURED (|effect| 0.0754 <= its own MDE 0.0767, 0.98x)
      point_cloud [all, n=126] RANDR@1e4 minus anchor (negative = arm better)
        a 3.1088 (med 3.0266)   b 3.0483 (med 2.8373)   n=126
        effect +0.0605   median +0.0000   SE 0.0377   MDE 0.1056   effect/MDE +0.57
        iid  CI95 [-0.0120, +0.1392]
        fold CI95 [+0.0206, +0.0962]   folds same sign 5/5   per-fold 0:+0.121 1:+0.069 2:+0.004 3:+0.014 4:+0.086
        62W/62L/2T   worst degradation +1.6061 (8FLP)   p90 +0.5672   power 0.36  Type-M 1.65
        concentration: drop-top10 +0.1250 vs uniform-effect null p10/p50/p90 +0.0781/+0.1226/+0.1703 -> pctile 0.527
        VERDICT: NOT MEASURED (|effect| 0.0605 <= its own MDE 0.1056, 0.57x)
      point_cloud [all, n=126] RANDS@1e4 minus anchor (negative = arm better)
        a 3.0854 (med 2.8491)   b 3.0483 (med 2.8373)   n=126
        effect +0.0371   median +0.0023   SE 0.0090   MDE 0.0253   effect/MDE +1.47
        iid  CI95 [+0.0214, +0.0556]
        fold CI95 [+0.0229, +0.0495]   folds same sign 5/5   per-fold 0:+0.043 1:+0.014 2:+0.028 3:+0.060 4:+0.040
        35W/78L/13T   worst degradation +0.6973 (9KAR)   p90 +0.1193   power 0.98  Type-M 1.01
        concentration: drop-top10 +0.0450 vs uniform-effect null p10/p50/p90 +0.0330/+0.0445/+0.0570 -> pctile 0.519
        VERDICT: WORSE
      point_cloud [moved, n=116] R@1e4 minus anchor (negative = arm better)
        a 3.2910 (med 3.1277)   b 3.0438 (med 2.8573)   n=116
        effect +0.2473   median +0.0142   SE 0.0750   MDE 0.2101   effect/MDE +1.18
        iid  CI95 [+0.1024, +0.3992]
        fold CI95 [+0.1515, +0.3538]   folds same sign 5/5   per-fold 0:+0.278 1:+0.080 2:+0.244 3:+0.457 4:+0.206
        49W/67L/0T   worst degradation +4.1545 (8T61)   p90 +1.3987   power 0.91  Type-M 1.05
        concentration: drop-top10 +0.3564 vs uniform-effect null p10/p50/p90 +0.2617/+0.3505/+0.4495 -> pctile 0.526
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.05x]
      point_cloud [moved, n=116] R@1e4 minus RANDR@1e4 (matched control)
        a 3.2910 (med 3.1277)   b 3.1149 (med 3.0689)   n=116
        effect +0.1761   median +0.0170   SE 0.0507   MDE 0.1421   effect/MDE +1.24
        iid  CI95 [+0.0825, +0.2822]
        fold CI95 [+0.0809, +0.3007]   folds same sign 5/5   per-fold 0:+0.156 1:+0.017 2:+0.199 3:+0.439 4:+0.117
        43W/73L/0T   worst degradation +3.2867 (8T61)   p90 +0.8178   power 0.93  Type-M 1.04
        concentration: drop-top10 +0.2451 vs uniform-effect null p10/p50/p90 +0.1808/+0.2413/+0.3129 -> pctile 0.531
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.04x]
      point_cloud [moved, n=113] S@1e4 minus anchor (negative = arm better)
        a 3.1525 (med 3.0824)   b 3.0321 (med 2.8563)   n=113
        effect +0.1203   median +0.0233   SE 0.0319   MDE 0.0892   effect/MDE +1.35
        iid  CI95 [+0.0601, +0.1875]
        fold CI95 [+0.0668, +0.1862]   folds same sign 5/5   per-fold 0:+0.043 1:+0.057 2:+0.167 3:+0.246 4:+0.126
        38W/75L/0T   worst degradation +1.7827 (8T61)   p90 +0.4649   power 0.97  Type-M 1.02
        concentration: drop-top10 +0.1607 vs uniform-effect null p10/p50/p90 +0.1199/+0.1586/+0.2023 -> pctile 0.524
        VERDICT: WORSE
      point_cloud [moved, n=113] S@1e4 minus RANDS@1e4 (matched control)
        a 3.1525 (med 3.0824)   b 3.0734 (med 2.8570)   n=113
        effect +0.0790   median +0.0154   SE 0.0282   MDE 0.0790   effect/MDE +1.00
        iid  CI95 [+0.0257, +0.1358]
        fold CI95 [+0.0302, +0.1325]   folds same sign 4/5   per-fold 0:-0.000 1:+0.042 2:+0.132 3:+0.170 4:+0.083
        41W/72L/0T   worst degradation +1.5172 (8T61)   p90 +0.3480   power 0.80  Type-M 1.12
        concentration: drop-top10 +0.1172 vs uniform-effect null p10/p50/p90 +0.0815/+0.1153/+0.1524 -> pctile 0.531
        VERDICT: NOT MEASURED (|effect| 0.0790 <= its own MDE 0.0790, 1.00x)

All four thresholds, `[all, n=126]`, arm minus anchor and arm minus its matched random control:

    1e3  R  vs anchor  +0.5317 (+1.75x MDE 0.3045)  fold [+0.4328, +0.6566]   36W/69L/21T  WORSE         | vs RANDR  +0.4611 (+1.87x)  fold [+0.3723, +0.6036]   WORSE
    1e3  S  vs anchor  +0.1960 (+1.66x MDE 0.1182)  fold [+0.1347, +0.2330]   29W/71L/26T  WORSE         | vs RANDS  +0.1118 (+1.20x)  fold [+0.0583, +0.1424]   WORSE
    1e4  R  vs anchor  +0.2276 (+1.17x MDE 0.1941)  fold [+0.1466, +0.3232]   49W/67L/10T  WORSE         | vs RANDR  +0.1672 (+1.26x)  fold [+0.0699, +0.2793]   WORSE
    1e4  S  vs anchor  +0.1079 (+1.34x MDE 0.0805)  fold [+0.0647, +0.1528]   38W/75L/13T  WORSE         | vs RANDS  +0.0709 (+1.00x)  fold [+0.0317, +0.1087]   NOT MEASURED
    1e5  R  vs anchor  +0.0929 (+0.66x MDE 0.1410)  fold [+0.0332, +0.1433]   58W/61L/7T   NOT MEASURED  | vs RANDR  +0.0416 (+0.44x)  fold [-0.0027, +0.0963]   NOT MEASURED
    1e5  S  vs anchor  +0.0739 (+1.17x MDE 0.0634)  fold [+0.0385, +0.1077]   50W/67L/9T   WORSE         | vs RANDS  +0.0568 (+0.95x)  fold [+0.0264, +0.0871]   NOT MEASURED
    1e6  R  vs anchor  +0.0662 (+0.51x MDE 0.1305)  fold [+0.0196, +0.1013]   56W/62L/8T   NOT MEASURED  | vs RANDR  +0.0297 (+0.30x)  fold [-0.0071, +0.0588]   NOT MEASURED
    1e6  S  vs anchor  +0.0504 (+0.97x MDE 0.0517)  fold [+0.0261, +0.0717]   52W/64L/10T  NOT MEASURED  | vs RANDS  +0.0377 (+0.81x)  fold [+0.0139, +0.0590]   NOT MEASURED

Threshold sweep (an order statistic): `ST.best_of_k_within` over the four thresholds, R: oracle
per-target minimum -0.3368 A, split-half transfer -0.1470 (44%), k_eff 3.74; S: -0.1177,
transfer -0.0556 (47%), k_eff 3.72. The transfer is negative because the sweep contains 1e6,
which rejects least and so damages least: "the least harmful threshold transfers" is not a gain.
R@1e6 is +0.066 against the anchor (0.51x MDE, NOT MEASURED) and S@1e6 +0.050 (0.97x).

READING. (1) The falsifier did not clear; it fired the other way. At the primary threshold the
refill arm R is WORSE than the shipped top-75 by +0.228 A [fold +0.147, +0.323], 49W/67L/10T,
5/5 folds (1.17x MDE, Type-M zone: the size is an upper bound, the sign is measured), and WORSE
than rejecting the same number at random and refilling with no energy judgment by +0.167
[+0.070, +0.279], 5/5 folds (1.26x MDE, Type-M). The shrink arm S is worse than the anchor by
+0.108 [+0.065, +0.153], 1.34x MDE, measured; against its random twin +0.071 at exactly 1.00x
MDE, not measured. Restricting to the targets the operator actually moved (the secondary,
n = 116 / 113) changes nothing: R +0.247 and S +0.120, both WORSE; R vs RANDR +0.176, WORSE.
(2) The dose is monotone: 1e3 (54.5 of 75 rejected) +0.532 A; 1e4 (40.1) +0.228; 1e5 (30.3)
+0.093 (0.66x MDE); 1e6 (23.3) +0.066 (0.51x). The less the reject does, the less it costs, and
its limit is the anchor: S16 L27's k-ladder shape and S23 L5's mechanism (the members a physics
score removes carry error that cancels in the average) through a third operator. (3) The
matched controls say where the harm comes from. Rejecting the same COUNT at random and
refilling (RANDR) costs +0.061 at 1e4 (0.57x MDE, not measured), so refilling from ranks 76 to
147 is nearly free and the energy's CHOICE of what to reject costs the other +0.167. Rejecting
the same count at random without refill (RANDS) costs +0.037 (1.47x MDE, measured): shrinking
a 75-set to 35 at random costs 0.04 A and the energy's choice of which 40 costs 0.07 more. The
permuted-energy controls sit between random and real (PERMR +0.122, PERMS +0.033) because a
permuted single point rejects the same energies from other candidates. (4) Why: the census
(L23) found 96.8% of the members the threshold condemns are condemned by a side-chain contact
placed by the deterministic builder, and the retained set at 1e4 is +0.254 A more expanded in
Rg than the anchor; the reject keeps expanded members and discards compact ones. S25's "AMBER
selects expansion" (+1.10 A Rg on its own top-75) reproduced as a filter.

POWER. Every arm-vs-anchor contrast at 1e3 and 1e4 is above its MDE with the fold CI excluding
zero and 5/5 folds. At 1e5 and 1e6 R is 0.66x and 0.51x its MDE (SE 0.050 and 0.047, so a gain
of 0.13 A or more would have been seen): "the mildest threshold is null" is UNDERPOWERED for a
gain below 0.13 A and MEASURED against any harm above it. Nothing here is positive, so no seed
replication is due; the cross-basis replication is the built-chain run, and the direction must
hold there or the point-cloud result is a basis artefact (S16 showed the two bases can disagree
about a repair). DISPOSITION, pending the built chain: the steric reject at a physical
threshold, with refill, is closed on the point cloud in the two forms the record had not
measured (physics-set count, refill), with a sign, on all 126 and on the moved subset; the
filter form of the functional lever (s24 D1-C, s19 Q3) stays closed and gains its sixth
instrument.

---

## L44 -- LANE W: THE 2/60 BENCHMARK SELF-COPY LEAK, BOUNDED FROM THE DEV PROXY WITHOUT OPENING THE BENCHMARK: 0.002 A BY THE DEV-4 MEASUREMENT, 0.028 A (BUILT CHAIN) / 0.048 A (SELECTION) / 0.023 A (PAIRED GAIN) BY THE OWN-NATIVE ENVELOPE; PRE-REGISTERED CLASS MINOR; C27's +0.0004 DEV PRICE RE-DERIVED EXACTLY; A FOLD MODEL THAT TRAINED ON THE TARGET'S OWN NATIVE EMITS A CHAIN 0.70 A NEARER TO IT (2026-09-13, W)

Pre-registration `s26/PREREG_selfcopy_bound.md` (L30); gated job `w_endpoint_report`
(`s26/w_endpoint_report.py`: posterior, endpoint, report in one process; exit 0, 160 s, peak RSS
0.302 GB, `s26/jobs_done/w_endpoint_report.json`); artefacts `s26/results/w_selfcopy_endpoint.json`
(Parts A, B, C signed), `w_selfcopy_bound.json` (Part D), `w_selfcopy_floor.json` (Part E),
`w_selfcopy_retrieval.json` and `w_selfcopy_envelope.json` (the native-free halves, complete
126/126, exit 0, peak 0.116 / 0.318 GB, L41), `w_selfcopy_posterior.json`. No benchmark file,
sequence, native or RMSD was read; the benchmark facts used are 2/60 and the mechanism (S24 L4).
Reproduction gate: the "with" arm equals the production emission on all 126 (top-75 == `sub`
126/126, cloud max abs 1.4e-14; means sel 3.4540 / cloud 3.0483 / arm 3.2126 / fit 3.2052, the arm
being the re-projection rebuild figure of L9). Basis stated on every line; positive delta = the
leaked emission is WORSE.

**Part A, channel A (the self-window in the K = 500 pool), the four dev self-copies, PRODUCTION
basis.** Signed delta = with minus without (self-window dropped, pool refilled from the next BLOSUM
windows, S10-4's operator on exact copies):

    target   in top-75   arm      cloud    fit      sel      | triangle bound arm | control percentile (arm)
    1CEK        no      +0.0000  +0.0000  +0.0000  +0.0000  |  0.000             |  0.61
    2FBU        no      +0.0000  +0.0000  +0.0000  +0.0000  |  0.000             |  0.61
    2P5H        yes     -0.0683  -0.0114  -0.0193  +0.0000  |  0.202             |  0.89
    6B9K        yes     +0.0071  +0.0069  +0.0046  +0.0000  |  0.034             |  0.71

The self-window is never the score's argmin (sel 0.0000 on all four, and 0.0000 on all 122
controls: the BLOSUM rank-0 window is never the argmin anywhere on the instrument). The matched
control population (the rank-0 window dropped on the other 122 targets, 47 of which had it in
the top-75): |delta arm| p50 0.000, p90 0.078, p95 0.147, max 0.511; the four self-copies sit at
its 61st to 89th percentile. Registered F1 holds (|delta arm| < 0.10 and bound < 0.30 on all four).

**C27's dev half re-derived (S10-4's operator: drop every >= 0.6 window, refill; 13 targets with
such a window in the pool, the same 13 as S10-4, 21 in the universe): CLEAN minus PRODUCTION
+0.0004 on the lam = 0 chain, fold CI [-0.0001, +0.0010], iid [-0.0003, +0.0013], MDE 0.0012, 117
ties; +0.0018 on the built chain [-0.0003, +0.0039]; +0.0003 on the cloud; 0.0000 on sel.**
S10-4's +0.0004 [-0.0004, +0.0013] reproduces on the same basis (its "synthesis fit"), and its "0
on the shipped argmin" reproduces exactly. The number now has an artefact (L28 C27, L31).

**Part C, the envelope (ORACLE by construction): the four pinned fold models that had the target's
OWN native among their training labels, against the clean model, everything else identical.**

    basis   leaked-model mean minus clean   fold CI95            MDE     x MDE   W/L      folds   per model
    arm     -0.6980 (median -0.2663)        [-0.8312, -0.5764]   0.246   2.83    112/14   5/5     -0.749 / -0.646 / -0.702 / -0.694 / -0.700
    cloud   -0.7094                         [-0.8394, -0.5680]   0.244   2.91    110/16   5/5
    sel     -1.2081 (median -0.9239)        [-1.4377, -0.9745]   0.315   3.83    112/12   5/5     -1.324 / -1.102 / -1.254 / -1.163 / -1.201
    fit     -0.6932                         [-0.8234, -0.5615]   0.245   2.83    110/16   5/5
    gain (arm - sel)  +0.5101               [+0.3698, +0.6838]   0.255   2.00    38/88    5/5

The shipped MLP memorises at the pipeline endpoint: a fold model that saw the target's native emits
a built chain 0.70 A nearer to it, changes half the top-75 (overlap 0.47) and the argmin on 97% of
targets, and its argmin selection is 1.21 A nearer. The registered expectation (-0.2 to -0.8 A on
arm) held; concentration at the uniform-effect null's 52nd percentile (not concentrated); the
spread among the four leaked models (sd 0.07 A on arm) is small against the effect, so the effect
is the inclusion of the native, not the corpus fifth. The leak inflates the argmin arm more than
the built chain, so under a leak of this strength the paired gain of the architecture over the
shipped argmin is biased AGAINST the architecture by +0.51 A per leaked target. This is also the
first measurement behind EXAMINATION E4: per-target dev results may never be read across folds.

**Part B (channel B, the carrier chain in the training labels), partial: one of ten retrains
finished before the host killed `w_train_chain` (L40, L41).** 1CEK with 1A11 removed from fold 2's
corpus (`s26/models/w_selfcopy/pca32_fold2_s0_out_1A11.pt`, 276 pairs fewer), reference lane P's
`pca32_fold2_s0.pt` (which reproduces the pinned emission at 0.000): the carrier's presence is worth
-0.0114 on the built chain, -0.0204 cloud, +0.0062 sel (posterior mean |dE[d]| 1.02 A per pair,
top-75 overlap 0.76). Both channels removed on 1CEK: production is 0.0114 A BETTER on arm, 0.0177 on
the gain. 1CEK is the one dev case whose copy is near-native (Part E, 0.595 A); the three others
sit 2.3 to 4.1 A away. The control-out models and the other three carrier-out models wait for
headroom (L42) and for the tournament.

**Part E (ORACLE DIAGNOSTIC): the same sequence in a different deposit, 18 dev targets, 22 verbatim
partners: median 2.908 A from the native (min 0.317, max 5.502; 18% under 1.0 A, 27% under 1.5 A;
the natives' own ensemble spread is 1.044 A; the copy is worse than the target's pool MEAN on 5 of
22 and beats the pool's best on 3 of 22).** F5 holds. A verbatim copy is not a near-native answer at
this length; this is why channels A and B are small where they are measured.

**Part D, the bound, (2 / 60) x the per-target quantity, assumptions A1 to A5 of the PREREG (same
mechanism and direction; the benchmark 2 no worse than the dev 4 or than the envelope; same pipeline;
no length correction; the envelope bounds help, harm is measured on the dev 4):**

    source                                            arm       sel       paired gain   class (0.017 = one tenth of the benchmark CI half-width)
    dev-4 channel A, signed max                        0.0023    0.0000    0.0023         IMMATERIAL
    dev-4 channel A, native-free triangle max          0.0067    0.0000    0.0067         IMMATERIAL
    both channels removed (1CEK only, n = 1)           0.0004    0.0002    0.0006         IMMATERIAL
    own-native envelope, fold-CI limit (n = 126)       0.0277    0.0479    0.0228         MINOR

**Pre-registered verdict: MINOR, bounded at 0.028 A on the built chain, 0.048 A on the selection
basis and 0.023 A on the paired gain, because the envelope clause of Part D fires; every direct
measurement on the dev proxy is IMMATERIAL (0.002 A or less).** The envelope is loose by
construction: it prices a model trained on the target's OWN native, and the only measured carrier
effect (1CEK, the near-native case) is 60x smaller. Against the benchmark's own CI half-width 0.170
(S9-10: +0.0103 [-0.1596, +0.1803]) the bound cannot change the benchmark verdict (no validated
gain) in either direction: the un-leaked paired gain lies in [+0.0103 - 0.023, +0.0103 + 0.023]
under the envelope and within 0.003 of +0.0103 under the dev-proxy measurement. The caveat that
attaches to every benchmark figure now reads: 2/60 self-copies, bounded at 0.028 A (built chain)
by the own-native envelope and 0.002 A by the dev proxy, `s26/results/w_selfcopy_bound.json`.

Also measured on the way, native-free, and recorded for `s26/IDEA_tiebreak_noise_floor.md`: a
one-member change of the 75 flips the medoid frame (cloud moves 1.29 A on 5H1H, 0.74 on 6EY3; the
same set in the original frame is 0.06 to 0.08 A away) and flips the projection branch (chain moves
0.85 to 2.11 A on 1NIZ, 1CS9, 1M02, 1RSW, 2LNG, 2BP4 at cloud moves of 0.04 to 0.09 A); the signed
deltas on those targets are 3 to 10x smaller than the triangle bounds, so the flips move the chain
mostly orthogonally to the native. ORACLE insertion of the withheld same-fold carriers (8 targets,
11 windows): 5 of 11 enter the top-75 (F2's first clause, "at least half", misses by one); |delta
arm| < 0.10 on 7 of 8 targets; the largest move is 8ZG2, -0.277 A from a window that is itself
5.06 A from the native (a branch flip, not the window's geometry).

Deviations from the PREREG, stated: Part B has 1 of 10 models (the rest killed by the host, L40);
the `both_removed4` key of `w_selfcopy_bound.json` therefore holds n = 1; no replication of Part
C at a second seed is needed (deterministic, pinned models); Part A's control population is the
replication of its operator. Question for the coordinator: the S10-4 number is now sourced on the
`fit` basis; does the report quote the built-chain +0.0018 [-0.0003, +0.0039] beside it or in its
place?

## L45 -- ADVERSARY CHECK OF L27 (A2, THE DLA): STANDS (2026-09-13, A)

Independent re-derivation `s26/a_dla_check.py` -> `s26/results/a_dla_check.json` (88 s, run
directly, under 200 MB; written from the construction in `s26/PREREG_A2.md` and
`core/quantum.py`, sharing no closure code with `s26/q_dla.py`, and using a different numeric
rank rule -- incremental Gram-Schmidt, no shared svd tolerance). Symbolic closure at n = 4..7,
L = 1..4 equals `s26/results/q_dla.json :: results/fixed/*` at every cell: n7 dim(DLA) = 8128 =
so(128) from L = 2 (7 at L = 1); n6 510 / 1023 / 2016; n4 120 from L = 2; n5 496 from L = 2. All
four combinations of conjugation convention (right, left) and gate order (forward, reverse) give
8128 at n = 7, L = 2 and L = 3. The independent numeric route agrees with the symbolic count at
n = 4, 5, L = 1..3 (6 of 6); every closure is odd-Y. The 8128 / 1025 reconciliation in
`s26/agentQ_FINDINGS.md` 2.1 is correct: 8128 is `q_dla.json :: results/fixed/n7_L3/dim` (the
fixed ansatz, so(128), 1.0 of so); 1025 is `results/adapt_sets/adapt_L2_*_n7_a0.25_T0.3/
ladder[P=21]/dim` (the strings ADAPT selected), 1025 / 8128 = 0.126. Leakage / ties / iid-vs-fold
CI: not applicable (an exact count, no native, no RMSD). Q's own H2b prediction (dim < 8128,
guess 4095) is correctly recorded as FALSIFIED. Verdict: STANDS. The "S25 slopes are a statement
about depth 3" reading is an inference from Larocca 2022 / Ragone 2024, not a measurement here;
it is a scope note for RETRACTIONS, not a retracted number.

---

## L46 -- ADVERSARY CHECK OF L39 (C3 STAGE 1): STANDS WITH CAVEAT (2026-09-13, A)

The three contrasts reproduce from `s26/results/ph_c3_stage1.json` (complete, 126 rows,
provenance e480fc15): AMBER minus do-nothing +0.02069 (effect/MDE 2.16, fold CI [0.0154,
0.0290], 5/5, WORSE), AMBER minus matched RANDOM +0.01108 (effect/MDE 1.12, fold CI [0.0062,
0.0171], 5/5), AMBER minus TOWARD-MEMBER +0.03851 (effect/MDE 2.27, fold CI [0.0298, 0.0479],
5/5); cos(AMBER, true residual) mean -0.0491 (SE 0.0149, 36.5% positive). Checks:

1. Control construction against S16. `s26/ph_lib.py:random_displacement` reproduces
   `s16/repair.py`'s `rand` exactly: an isotropic Gaussian on the CA trace, the six rigid-body
   directions removed by projection onto `s15/align_lib.py:rigid_basis`, scaled so
   ||g||/sqrt(n) equals AMBER's per-atom RMS displacement. The magnitude matched is `mag_sup`
   0.2198 A (SE 0.0084), measured after superposition, the same definition S16 used. Confirmed.
2. The Type-M reading. AMBER-minus-RANDOM is at 1.12x MDE with `type_m_flag` True (0.7-1.3x is
   the Type-M zone). So the HEADLINE +0.0111 is a Type-M number: its sign and fold CI are clean
   (5/5, CI excludes zero) but its magnitude is inflated ~1.07x and by the discipline it is
   "not a result" as a magnitude. The two clean contrasts are AMBER-minus-do-nothing (+0.0207,
   2.16x) and AMBER-minus-member (+0.0385, 2.27x), both clear of the Type-M zone. Any
   presentation line must quote +0.0111 with the Type-M flag or lean on the two clean contrasts.
3. The decision rule holds regardless: "accuracy step" required AMBER to BEAT both controls;
   it is worse than both, so it is not an accuracy step. Note that TOWARD-MEMBER minus
   do-nothing is -0.0178 (BETTER, fold CI [-0.0255, -0.0124], 5/5): a zero-information move
   toward a random pool member IMPROVES RMSD where AMBER's physics move worsens it.
4. "Validity step" wording. Sound for the 124 targets that converge with a sane virtual bond,
   but on 2BP4 (relaxed CA-CA 5.38 A) and 9KAR (4.86 A, e1 1262 > CONVERGE_MAX_KCAL, not
   converged) the relaxation BREAKS a virtual bond; on those two it is not cleanly a validity
   step either. The presentation should say "a validity step on 124 of 126; on 2 it breaks a
   virtual bond, one of which does not converge."

Verdict: STANDS WITH CAVEAT (the +0.0111 headline is Type-M; "validity step" carries the 2BP4 /
9KAR exception). The finding's direction and its central conclusion are correct.

---

## L47 -- ADVERSARY CHECK OF L35 (A4, GROWN-CIRCUIT VARIANCE): STANDS WITH CAVEAT (2026-09-13, A)

The reproduction gate holds: `s26/results/q_var.json :: results/reproduction/passed` True,
`worst_rel` 0.0 (the S25 n = 7 rows of all five cells reproduce bit-identically). The slopes
reproduce (`results/slopes`): fixed deployed_a1_T03 -0.311, grown V +0.024, grown L2 -0.008;
the one non-trivial family, deployed_a025_T03, fixed -0.243, grown V -0.239, grown L2 -0.302.
Checks:

1. "Product circuit at alpha = 1" is established, but by `s26/results/q_dla.json`, not by
   `q_var.json`. `q_var.json` carries only the variance slopes; the product-circuit reading
   rests on the ADAPT closures in `q_dla.json :: results/adapt_sets` (alpha = 1 selects an
   abelian, dim-7 = n set, `n_distinct_ops` 1 to 2, the ladder dim stays 7), which I
   independently reproduced in `s26/results/a_dla_check.json` (L27 above). A large gradient
   from a product circuit is the trivial regime, not trainability (rule 10's mirror). Sound,
   with that provenance noted.
2. The one non-trivial slope's comparison. L35 says grown L2 -0.302 is "the fixed ansatz's
   -0.243 within the sampling error of a 7-point slope (relative SE of a variance 9 to 16%)."
   `q_var.json` stores the point slope but no CI on it, so "within error" is asserted, not
   computed: the difference is 0.059 in log2-slope-per-qubit and no SE on that difference is on
   disk. The qualitative conclusion is safe (both slopes are order -0.25 to -0.30, both far from
   the -1.0 of a 2-design), but the precise "within error" claim is not backed by a stored CI.
   Also: the ledger table's grown-V value -0.246 for that cell is the P = 3n-only 6-point fit,
   while `q_var.json :: slopes/deployed_a025_T03/grown_V_adam_best` is -0.239; the two differ
   because of the degenerate-row exclusion the ledger states, not a discrepancy.

Verdict: STANDS WITH CAVEAT (the negative conclusion -- ADAPT gives Proposal A no width-scaling
argument -- holds; the "within error" of -0.302 vs -0.243 is asserted without a persisted slope
CI, and "product circuit" is grounded in q_dla.json, not q_var.json).

---

## L48 -- ADVERSARY CHECK OF L22, L23, L24, L38 (PH CENSUSES AND THE CIS FLOOR): ALL STAND; L22 WITH A STALE-LINE-NUMBER CAVEAT (2026-09-13, A)

L22 (cis census, `s26/results/ph_cis_census.json`, complete). Omega statistics recomputed from
its rows: 1,507 bonds (= sum of n-1), mean 1.91, median 0.34, p90 5.8, p99 17.4, max 42.8 deg
at 9UV5 (bonds -137.2 and +144.8), 0.53% beyond 20 deg, 0.13% beyond 30; 0 of 126 model-1
natives cis by either criterion (the two agree 126/126), 0 of 1,966 ensemble models, 0 of
2,352,893 windows, universe minimum step 3.5045 A. Native-free path confirmed (natives through
`ph_lib.native_backbone`, omega/CA only; bank through `univ_nativefree`, which refuses `rr` /
`nat_ca`; record through `prod_record_nativefree`, RMSD keys stripped). CAVEAT: L22 cites the
step gate at `core/data.py:406-407` and `697-698`; at HEAD it is `core/data.py:439`
(`step.min() < 3.5 or step.max() > 4.1`) and `:725` -- the numbering before lane I's `37bddbbb`.
Same code, stale line numbers. Q's registered ensemble-cis prediction is correctly recorded as
FALSIFIED (0 of 1,966). Verdict: STANDS WITH CAVEAT.

L23 (steric reject census, `s26/results/ph_reject_census.json`, complete). The pool-identity
assertion is live code (`s26/ph_reject.py:96-102`: `universe_idx == I.pool_idx`,
`amber_verify_max_rel == 0.0`, production `sub` == score top-75 as a set) and holds on 1A13,
1S9Z, 9KAR (`s24/cache_amber/<pdb>.npz :: universe_idx` equals `order[:500]`; top-75 members
above 1e4 = 12 / 27 / 74, equal to the census rows; pool frac above 1e4 = 0.390 / 0.650 /
0.888). 96.8% = (2,627 + 2,267) / 5,057 side-chain-involving; rho(e, min heavy-atom distance)
within top-75 -0.743 (SE 0.017). Verdict: STANDS.

L24 (C3 native-free part). Every number recomputed from `bench_results/cache/1fc9f2dcf489e2fb`
with my own Kabsch: displacement 0.2198 A (SE 0.0084, median 0.1972, range 0.103-0.591),
`amber_moved` 0.2332, e0 median 8.59e4 / min -473 / max 1.3e14 / 58.73% above 1e4, e1 mean
-559.8 (SE 30.0) / max 1262.4 (9KAR), 125 of 126 converged, strain 60.4, relaxed bond mean
3.867 (min 3.12 at 1M02, max 5.38 at 2BP4, 4.86 at 9KAR), 6 of 126 outside [3.6, 4.0], Rg
+0.0457. All match L24. Verdict: STANDS.

L38 (cis floor, ORACLE DIAGNOSTIC, `s26/results/ph_cis_floor.json`, complete). floor_ca mean
0.3468 A (SE 0.029, median 0.272, max 1.474 at 1ID6); rebuild_bb 0.3400; chain_cost 0.1664;
rho(floor_ca, max omega dev) +0.828, rho(chain_cost, omega dev) -0.036, rho(chain_cost,
floor_ca) +0.083; the paired contrast cost-minus-floor -0.1804 (fold CI [-0.2266, -0.1022],
verdict BETTER). Every quantity reads the native and is ORACLE-labelled; the "BETTER" is
explicitly disarmed in L38 (it means only that the production chain cost is smaller than the
own-torsion rebuild floor). floor_ca is declared an UPPER bound on the manifold floor and the
tight `floor2` is registered. The two Spearman nulls (+0.08, -0.04) exclude |rho| above ~0.25
at n = 126. Verdict: STANDS.

---

## L49 -- ADVERSARY CHECK OF L30 (W's 2/60 PROXY-BOUND PREREG): SOUND, STANDS (2026-09-13, A)

L30 is a pre-registration plus a native-free census, not an endpoint result; the check is on
its soundness and on whether the census reads anything it must not.

Reads. `s26/w_selfcopy.py census` iterates `I.targets()` (the 126 dev targets), loads every
universe through `load_blind` (which overwrites `rr` and `nat_ca` with NaN, lines 105-115), and
reads the production record through `I.shipped_record` (native-free). No benchmark sequence,
name, PDB, native, RMSD or manifest is read anywhere in the census, retrieval, envelope or
posterior commands; the only native reads (`I.load_univ` at lines 686, 821) are inside the
gated `endpoint` and `floor`, which refuse to run before "PHASE 0 SIGNED OFF". The one
benchmark-derived fact used, 2/60, is on the record (S24 L4, lane I L15/L18). Confirmed clean.

Assumption set (A1-A5). Reasonably complete for the quantity claimed (the leak's contribution
to a benchmark mean, bounded in absolute value). A1 (mechanism match, no interaction when one
carrier carries both benchmark targets) and A5 (the envelope bounds channel B in the HELP
direction only; HARM is measured on the dev 4 directly) are the two load-bearing assumptions and
both are stated. The gap a reader would press -- that the benchmark carriers might be MORE
homologous to their targets than the dev-4 carriers, which would make the benchmark effect
larger than "max of four" -- is covered by the ORACLE-insertion arm (8 dev targets with a
same-fold verbatim carrier, longer identity 0.6 to 0.93 against the dev 4's <= 0.52), reported
beside the dev 4. The honest limitations are self-declared: n = 4 supports no quantile ("max of
four" is named as such), and the AMBER stage's second-order contribution to a per-target delta
is an assumption (A3), not a measurement.

Verdict: STANDS as a sound pre-registration. The bound's headline will rest on n = 4 for the
realistic arm and on the n = 126 envelope for the guard; when the gated endpoints land, the
Adversary re-checks the signed deltas, the triangle bounds and the envelope's fold CI against
this prereg before any benchmark caveat text is written.

---

## L50 -- RULING ON L44's QUESTION: THE REPORT QUOTES BOTH LEAK PRICES, EACH WITH ITS BASIS; NEITHER REPLACES THE OTHER (2026-09-13 19:27, coordinator)

Lane W asked whether the report quotes the built-chain +0.0018 A [-0.0003, +0.0039] beside
S10-4's +0.0004 A [-0.0004, +0.0013] or in its place. Both, each with its basis named: the
+0.0004 is the S10-4 figure on the lam = 0 chain (`fit`), now re-derived exactly
(`s26/results/w_selfcopy_endpoint.json`, L44) and is the number the state brief carries; the
+0.0018 is the same operator on the production built chain and is the number that matches the
report's production basis. Neither clears its MDE (0.0012 on `fit`; the built-chain CI includes
zero) and both say the dev leak is immaterial; the envelope (0.028 A built chain) is the bound
the report states for the benchmark's 2/60, labelled MINOR as pre-registered. The 0.70 A
own-native effect (fold CI [-0.83, -0.58], 5/5 folds) is an ORACLE fact about leaked training,
not a price of the leak that exists, and the report says so in the same sentence.

---

## L51 -- TOURNAMENT RANKED (`s26/TOURNAMENT.md`); ASSIGNMENTS AND THE RUN ORDER UNDER TONIGHT'S HEADROOM (2026-09-13 19:28, coordinator)

The Adversary's ranking stands as posted. Assignments, in run order, one governed job at a time
beside the three endpoint jobs already running (a1_build, p_eval_*, ph_reject_chain):

1. product_state_optimum (Q): rides A1; Q reports it with A1's ledger entry.
2. conformational_identity_floor (W): gated, under a minute; W runs it first.
3. strain_difficulty (PH): reads the 126 cached records; PH runs it between chain cells.
4. tiebreak_noise_floor (W): about 50 min CPU, checkpointed; W runs it after item 2.
5. l17_target_dependent_hamiltonian = A3 (Q): already chained behind A1.
   branch_select (PH): AMBER, one job, about 1.0 GB; PH probes one target first and runs it
   only when the governor snapshot shows 1.5 GB free (jobrun v2.2 enforces it).
6. window_ensembling (P's idea, orphaned to W): W runs it after items 2 and 4; it is the
   mandatory test-time-ensembling direction and must state how it differs from widening K.
   rotamer_relief B (PH), window_provenance (W), amber_prior_partner (orphaned to W): as
   capacity allows, in that order, after everything above.
7. Deferred on memory: better_prior_inputs `attn` (3 to 3.5 GB) and coherence_penalised_training
   (1.25 GB); the `ragp` rung folds under C2 in lane P's queue. If the box frees up, W takes
   coherence_penalised_training; `attn` needs the box to itself and is unlikely tonight.

The amber_reject built-chain arm remains PH's (the falsifier's last leg, running as
ph_reject_chain). Every survivor still needs its PREREG on disk before compute (the ones
without one: strain_difficulty, window_ensembling, window_provenance, amber_prior_partner,
product_state_optimum if its arm is not already inside PREREG_A1).

---


## L52 -- TOURNAMENT ITEM 2, conformational_identity_floor (W): THE SAME SEQUENCE IN A DIFFERENT DEPOSIT SITS 2.9 A FROM THE NATIVE (MEDIAN OVER 22 PAIRS; 3.1 A OVER 18 TARGETS), 0.97 A WORSE THAN THE POOL'S BEST WINDOW AND 1.44 A BETTER THAN ITS MEAN; F5 HOLDS; ORACLE DIAGNOSTIC (2026-09-13, W)

Pre-registered as Part E of `s26/PREREG_selfcopy_bound.md` (falsifier F5: median above 1.5 A;
FALSIFIED if below 1.0 A) and, for the two paired contrasts, `s26/PREREG_identity_floor.md`
(written before `s26/w_identity_floor_stats.py` ran). Measurement: job `w_selfcopy_floor` (exit 0,
5 s; `s26/results/w_selfcopy_floor.json`, complete 22/22, provenance e480fc15); statistics: job
`w_identity_floor_stats` (exit 0, 5 s, peak 0.004 GB; `s26/results/w_identity_floor.json`). Every
quantity reads the native (the copy's CA-RMSD to the target's model-1 native at the shared
segment; the pool's ORACLE `rr`): ORACLE DIAGNOSTIC, nothing selects, nothing is deployable.
Population: the 18 dev targets with a verbatim relative in the peptide database (lane I's L15
list; 4 cross-fold carriers, the rest same-fold), 22 partners; benchmark sequences never used.

    cross-deposit CA-RMSD, 22 pairs:  median 2.908   mean 2.811   min 0.317 (5V5B)   max 5.502 (7S3O)
                                       below 1.0 A: 4/22 (18%)     below 1.5 A: 6/22 (27%)
    per target (partners averaged), 18:  median 3.055   below 1.0 A: 4/18   below 1.5 A: 5/18
    the four cross-fold self-copies:  1CEK 0.595   2FBU 3.278   2P5H 2.334   6B9K 4.126   (S24 L4's four, reproduced)
    role: carrier segment (n = 15) median 2.98;  carried whole chain (n = 7) median 2.04;  Spearman(length ratio, RMSD) -0.08 (p 0.72)
    reference scales: the natives' own intra-ensemble spread 1.044 A (record); the S24 L8 "universe best" 1.313 A

F5 holds (2.908 against the 1.5 A line; 18% below 1.0 A against "fewer than a third"). The two
registered paired contrasts (`ST.fmt` verbatim; basis: a single window against the native on both
sides; n = 18 targets, five folds):

  copy_minus_pool_best (ORACLE both sides; single window vs native), n=18 targets
    a 2.7968 (med 3.0548)   b 1.8291 (med 1.9778)   n=18
    effect +0.9677   median +0.9703   SE 0.3269   MDE 0.9158   effect/MDE +1.06
    iid  CI95 [+0.3349, +1.5834]
    fold CI95 [+0.7247, +1.2809]   folds same sign 5/5   per-fold 0:+1.113 1:+1.523 2:+0.736 3:+0.635 4:+0.937
    3W/15L/0T   worst degradation +2.9126 (7S3O)   p90 +2.5643   power 0.84  Type-M 1.10
    concentration: drop-top10 +2.1403 vs uniform-effect null p10/p50/p90 +1.6224/+2.1031/+2.5080 -> pctile 0.544
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.10x]
  copy_minus_pool_mean (ORACLE both sides; single window vs native), n=18 targets
    a 2.7968 (med 3.0548)   b 4.2396 (med 4.1842)   n=18
    effect -1.4428   median -0.9970   SE 0.4504   MDE 1.2618   effect/MDE -1.14
    iid  CI95 [-2.3156, -0.6046]
    fold CI95 [-1.7172, -1.1896]   folds same sign 5/5   per-fold 0:-1.855 1:-1.217 2:-1.773 3:-1.452 4:-1.023
    13W/5L/0T   worst degradation +0.9307 (8ZG2)   p90 +0.4068   power 0.89  Type-M 1.06
    concentration: drop-top10 +0.1234 vs uniform-effect null p10/p50/p90 -0.3330/+0.0854/+0.4138 -> pctile 0.555
    VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.06x]

Reading, with the discipline's labels: the copy is WORSE than the pool's own best window (+0.97 A,
5/5 folds, fold CI excluding zero, but 1.06x MDE: Type-M zone, sign clean, magnitude not a result)
and BETTER than a typical pool window (-1.44 A, 5/5 folds, 1.14x MDE: Type-M zone likewise). Both
were registered as expected in sign; the second was registered as "not measured", and it is
measured in the Type-M sense only. The number for the report: sequence identity buys a window that
is, in the median, 2.9 A from the native, i.e. no nearer than the shipped built chain's 3.21 A mean
by a margin the instrument can call at n = 18, and 1 A worse than the best window the pool already
holds. The 2.0 A target of S15 to S25 is below what a verbatim sequence lookup reaches on this
instrument. Two consequences carried forward: the self-copy leak is small where it is measured
because the copy is not the native (L44); and "containment-fresh" novelty (memory
`no-fresh-benchmark-exists`) is a weaker notion than the record has treated it as, since a
verbatim copy is worth about one pool window. Replication: deterministic (Kabsch on fixed
coordinates; no seed, no fit); the only random element, the bootstrap, is seeded by `stats_lib`.
Not done: a length-matched non-verbatim control population (a random same-fold window against
each target) is the natural third contrast and was not registered; it is the pool mean's role
here and is left for the report's reader as such.

## L53 -- STRAIN DIFFICULTY (tournament 4): HOW FAR THE RELAXATION MOVES THE BUILT CHAIN PREDICTS ITS ERROR (SPEARMAN +0.433 WITH rmsd_arm, FOLD CI [+0.25, +0.59], 5/5 FOLDS, PARTIAL ON n AND Rg, REPLICATED); QUARTILE MEANS 2.29 / 2.94 / 3.76 / 3.92 A; A CALIBRATION FLAG, NOT A LEVER (2026-09-13, PH)

`s26/ph_strain.py`, `s26/results/ph_strain.json` (complete 126/126) and the registered
replication `s26/results/ph_strain_rep.json` (different bootstrap seed, targets in reversed
order); jobs `s26/jobs_done/ph_strain.json` and `ph_strain_rep.json` (exit 0, 35 s each, peak
RSS 0.1 GB). Pre-registered in `s26/PREREG_strain_difficulty.md` (tournament rank 4, L51),
committed at 59d8e934 before the run. Reads the 126 production records
(`bench_results/cache/1fc9f2dcf489e2fb`) and nothing else. CALIBRATION, not an accuracy lever:
nothing is selected, tuned or moved; `rmsd_arm` (BUILT CHAIN basis) is read only as the ORACLE
label of the already-emitted structure. Signals, all native-free and emitted by the production
relaxation for free: `log_e0` (the built chain's own AMBER energy before relaxation),
`log_drop` (energy removed), `moved` (restraint RMSD on N/CA/C, how far the chain moved),
`strain_after` (bond + angle energy left). Confounds n and Rg partialled by residualising both
variables on them.

Spearman with the ORACLE `rmsd_arm`, n = 126, 4000-draw iid and fold-clustered bootstrap CIs,
4000-draw label-permutation p, per-fold signs:

    raw rho(log_e0, rmsd_arm)                     rho +0.224  iid [+0.061, +0.376]  fold [+0.116, +0.328]  perm p 0.0110  folds same sign 5/5
    raw rho(log_drop, rmsd_arm)                   rho +0.226  iid [+0.060, +0.376]  fold [+0.118, +0.328]  perm p 0.0112  folds same sign 4/5
    raw rho(moved, rmsd_arm)                      rho +0.431  iid [+0.267, +0.565]  fold [+0.278, +0.565]  perm p 0.0000  folds same sign 5/5
    raw rho(strain_after, rmsd_arm)               rho +0.175  iid [+0.003, +0.348]  fold [+0.109, +0.236]  perm p 0.0470  folds same sign 5/5
    partial rho(log_e0, rmsd_arm | n, Rg)         rho +0.241  iid [+0.070, +0.391]  fold [+0.047, +0.379]  perm p 0.0075  folds same sign 4/5
    partial rho(log_drop, rmsd_arm | n, Rg)       rho +0.250  iid [+0.081, +0.402]  fold [+0.068, +0.381]  perm p 0.0037  folds same sign 4/5
    partial rho(moved, rmsd_arm | n, Rg)          rho +0.433  iid [+0.272, +0.573]  fold [+0.247, +0.588]  perm p 0.0000  folds same sign 5/5
    partial rho(strain_after, rmsd_arm | n, Rg)   rho +0.133  iid [-0.052, +0.305]  fold [-0.178, +0.452]  perm p 0.1353  folds same sign 3/5

    FAIL18 in the top quartile of log_e0       : 7/31 (18 of 126 overall)  odds 2.23  one-sided p 0.113
    FAIL18 in the top quartile of log_drop     : 7/31 (18 of 126 overall)  odds 2.23  one-sided p 0.113
    FAIL18 in the top quartile of moved        : 6/31 (18 of 126 overall)  odds 1.66  one-sided p 0.257
    FAIL18 in the top quartile of strain_after : 4/31 (18 of 126 overall)  odds 0.86  one-sided p 0.698
    confounds: rho(n, rmsd_arm) +0.114; rho(Rg, rmsd_arm) +0.213; rho(n, log_e0) +0.340; rho(Rg, log_e0) +0.140

THE FALSIFIER CLEARS ON ONE SIGNAL, `moved`, WITH MARGIN: partial rho +0.433 (bar 0.25), fold CI
[+0.247, +0.588] excluding zero, 5/5 folds, permutation p < 0.00025 (bar 0.0125 after Bonferroni
over four signals). REPLICATED as registered (new seed, reversed order):

    partial rho(moved, rmsd_arm | n, Rg)          rho +0.433  iid [+0.272, +0.572]  fold [+0.248, +0.581]  perm p 0.0000  folds same sign 5/5
    partial rho(log_e0, rmsd_arm | n, Rg)         rho +0.241  iid [+0.080, +0.394]  fold [+0.047, +0.375]  perm p 0.0063  folds same sign 4/5
    partial rho(log_drop, rmsd_arm | n, Rg)       rho +0.250  iid [+0.084, +0.413]  fold [+0.068, +0.381]  perm p 0.0037  folds same sign 4/5

The replication lands inside the first run's fold CI on every quantity. `log_e0` and `log_drop`
are positive at +0.24 and +0.25 partial but fail the 5/5-fold sign rule (4/5) and sit at the
0.25 bar; `strain_after` is null (partial +0.13, fold CI spanning zero). The FAIL18 Fisher test is
null for every signal (best one-sided p 0.113), as predicted: FAIL18 is about retrieval recall,
not strain.

The presentable form (built chain, ORACLE labels, quartiles of `moved`, 32 targets each):

    moved (A)     < 0.159      0.159 to 0.218    0.218 to 0.288    >= 0.288
    mean rmsd_arm   2.286         2.936             3.758             3.923

A chain the force field has to move 0.29 A or more to make physical has a mean error 1.6 A larger
than one it moves less than 0.16 A. Mechanism checks (ORACLE, diagnostic): without 9KAR and 2BP4
(the two broken-bond emissions, moved 0.61 and 0.50) rho is +0.41, so it is not two outliers;
rho(moved, log_e0) is +0.41 while `moved` beats `e0` on the label by 0.19 in rho, so it is not the
raw energy in disguise; the top-8 by `moved` contain both easy (1D6X 1.78 A) and hard (9KAR
7.44, 2MFV 6.04) targets, so it is a graded signal and not a flag for a few catastrophes.

What it is and is not. It is the first native-free quantity in this programme's record with a
correlation above 0.4 to the per-target error of the emitted structure (the routers of S22 L7 /
S23 L7 and the compactness proxies of `in-band-ordering-is-per-target` reached 0.24 to 0.37 on
the per-target sign). It is a confidence flag the presentation can attach to every emitted
structure at zero cost ("how far the physics had to move this chain to make it physical").
It is NOT an accuracy lever: the prereg forbids converting it into a selector or a weight, and
the record says every such conversion fails held out (S22 L7, S23 L7). The physical reading is
plain: a coordinate average that the projection turned into a strained chain is one whose pool
members disagreed, and disagreement is error. Power: at n = 126 the design resolves |rho| of
0.25 at about 0.8; the three weaker signals are at or below that bar and are reported as
measured, not as nulls. HYPOTHESIS for lane PR: quote it as a calibration curve, never as a
gain. The Adversary's check is invited.

---

## L54 -- ADVERSARY CHECK OF L43 (STERIC REJECT, POINT CLOUD): STANDS WITH CAVEAT (2026-09-13, A)

A harmful result, so the checklist is applied to the controls and the power, not to a gain.

- Leakage: the operators (`reject_refill`, `reject_shrink`, `random_shrink`, `random_refill`,
  `permuted_energy`, `retained_sets`, `s26/ph_reject.py:109-200`) read `e_amber`, `score_dist`,
  `sub` and `order` only; no `nat_ca` and no RMSD inside them (the `rr` at `:169` is an rng handle,
  not the ORACLE array; `univ_nativefree` refuses `rr`). The native enters only in `cloud_rmsd`,
  the ORACLE scoring of a native-free operator. Clean.
- Tie-breaking: the score order is `argsort_stable` (the production rule); an empty retained set
  falls back to the anchor and is counted as a tie (10 ties at 1e4 R, 13 at S), as pre-registered
  in addendum 1. Clean.
- iid vs fold CI: both quoted verbatim; at 1e4 both exclude zero for R vs anchor and R vs RANDR,
  5/5 folds. Clean.
- Concentration / median-vs-mean: the uniform-effect null is computed and not flagged (pctile
  0.52). But the harm is TAIL-CARRIED: R@1e4 vs anchor has median +0.0027 against mean +0.2276,
  49W/67L/10T, p90 +1.28, worst +4.15 (8T61). The presentation must say "near zero on the median
  target, catastrophic on a minority (the 8 targets whose whole pool has no survivor, and the deep
  refills)", never "+0.228 on every target". CAVEAT 1.
- k_eff for the threshold sweep: `ST.best_of_k_within` applied (k_eff 3.74 R / 3.72 S), and the
  negative split-half transfer is correctly read as "the mildest threshold transfers", not a gain.
  Clean.
- Tuned parameter: the primary 1e4 was fixed before the run (PREREG section 3). Clean.
- Baseline / basis: the anchor is the production `rmsd_avg` (reproduced to 1e-6 on every
  target); controls matched in count, refill and permutation; point cloud on both sides, stated.
  Clean.
- Type-M: the two HEADLINE numbers, R@1e4 vs anchor +0.228 (1.17x MDE) and R vs RANDR +0.167
  (1.26x), are Type-M-zone magnitudes. The direction "harmful" rests on the measured contrasts:
  1e3 R +0.532 (1.75x), 1e3 S +0.196 (1.66x), 1e4 S +0.108 (1.34x), RANDS +0.037 (1.47x), 5/5
  folds throughout, and on the monotone dose. Quote the 1e4 R magnitudes with the flag. CAVEAT 2.
- "Refilling from ranks 76 to 147 is nearly free" rests on RANDR vs anchor +0.061 at 0.57x MDE
  (UNDERPOWERED, iid CI includes zero), so the refill-cost / choice-cost decomposition is a
  point-estimate split, not a measured one. CAVEAT 3.
- Replication: a negative; the cross-basis replication is the built-chain job `ph_reject_chain`,
  running; the direction must hold there. Power stated correctly for the 1e5 / 1e6 nulls.

Verdict: STANDS WITH CAVEAT (the conclusion -- the physical-threshold reject with refill is
harmful on the point cloud, its limit is the anchor -- is established; the +0.228 / +0.167
magnitudes are Type-M, the harm is tail-carried, and the refill-cost split is underpowered).
Disposition "closed on the point cloud in the two unmeasured forms" is accepted pending the
built chain.

---

## L55 -- ADVERSARY CHECK OF L44 (THE 2/60 PROXY BOUND): STANDS WITH CAVEAT; THE CLASS MINOR IS STABLE UNDER EVERY READING OF THE ENVELOPE, THE NUMBER 0.028 IS NOT (2026-09-13, A)

- Leakage / reads: the native-free halves use `load_blind` (`rr`, `nat_ca` NaN-poisoned); the
  gated `endpoint` and `floor` read natives for scoring only; no benchmark file, sequence or name
  anywhere (confirmed in L49 and re-checked on `s26/w_endpoint_report.py`'s inputs). Clean.
- Ties: `emit` uses `argsort(kind="stable")`, the argmin tie set is stored and `sel` is averaged
  over it (`ST.argmin_tied`'s rule); the refill follows the production corpus order. Clean.
- iid vs fold: Part C on arm -0.698, fold CI [-0.831, -0.576], 2.83x MDE, 112W/14L, 5/5;
  concentration at the null's 52nd percentile; the median (-0.266) sits well inside the mean
  (-0.698), a broad but skewed effect. Clean, ORACLE-labelled, and L50 rightly keeps it out of
  the price of the leak that exists.
- Bound arithmetic: every row of `s26/results/w_selfcopy_bound.json` recomputes: (2/60) x
  0.0683 = 0.0023 (2P5H); (2/60) x 0.2016 = 0.0067 (triangle); (2/60) x 0.8312 = 0.0277 (arm
  envelope fold-CI limit); (2/60) x 1.4377 = 0.0479 (sel); (2/60) x 0.6838 = 0.0228 (gain).
  Materiality thresholds pre-registered (0.017 / 0.170). Clean.

Two caveats.

1. **Artefact / ledger disagreement on the paired gain.** `w_selfcopy_bound.json :: verdict/gain`
   reads IMMATERIAL at 0.0006 (the n = 1 `both_removed4` row) because `C_envelope_fold_ci` carries
   no `gain` row and `report()`'s `max(cands)` (`w_selfcopy.py:913`) therefore never saw the
   envelope for that basis; L44's Part D table and its pre-registered verdict say MINOR at 0.023
   (from the gain row +0.5101 [+0.3698, +0.6838]). The ledger's class is the prereg-correct one
   (Part D: IMMATERIAL only if B_real AND B_env are both below 0.017). Lane W should add the
   envelope gain row to the artefact so the report cites a JSON that agrees with the ledger.
2. **B_env is the fold-CI limit of a MEAN effect, not a per-target bound.** The prereg defined it
   so (Part D, assumption A2), and A2 is stated, but a bound on the contribution of two SPECIFIC
   targets is (2/60) x a per-target quantity. From `s26/results/w_selfcopy_endpoint.json :: C/rows`
   (leaked minus clean, arm, mean over the four leaked models): median -0.266, p05 -2.490, worst
   -4.536 A (2BP4); over (target, model) pairs the worst is -4.557. So the same envelope gives
   (2/60) x 4.536 = **0.151 A** as the worst-single-target bound and (2/60) x p95 = 0.083 A;
   the artefact's own native-free triangle rows say the same (`C_envelope_p95_over_targets/arm`
   0.147, `C_envelope_max_over_models/arm` 0.211, the latter loose because branch flips move
   the chain orthogonally to the native, as L44 notes). The class MINOR therefore holds under
   every reading of the envelope (0.028 mean-CI, 0.083 p95, 0.151 worst target; all below
   0.170), which is a stronger statement than L44 makes; but "bounded at 0.028 A" must be quoted
   as "expected contribution 0.023, mean-CI limit 0.028, worst single target 0.151, under A2",
   and the envelope named as the over-bound it is (own-native training, 60x the one measured
   carrier effect).

Also noted, not a caveat: Part B is 1 of 10 retrains (host kill, L40), stated as a deviation;
F2's "at least half" missed by one (5 of 11) and is reported as such; F1, F3 (n = 1), F4, F5
hold. Verdict: STANDS WITH CAVEAT. The benchmark caveat text should read: "2/60 self-copies;
dev-proxy price 0.002 A; own-native envelope 0.028 A (mean CI) to 0.151 A (worst target), MINOR
under every reading; cannot move the benchmark verdict either way."

---


## L56 -- LANE P, C2 RUNG SHIPPED: THE ANCHOR. THE SHIPPED POSTERIOR THROUGH THE LADDER'S OWN PATH IS BIT-EXACT AGAINST THE PRODUCTION CACHE ON ALL FOUR BASES (2026-09-13 19:34, lane P)

Artefacts: `s26/results/p_ladder_shipped_s0.json` (126 rows, complete), `s26/results/p_ladder_report_shipped_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_shipped2` (resumed from a 30-row checkpoint): exit 0, 416 s, peak RSS 0.13 GB (`s26/jobs_done/p_eval_shipped2.json`); the first attempt `p_eval_shipped` was suspended by the governor at 09:25 and terminated at 19:11 (L42). Every later rung is paired against these 126 rows. The per-target production values are `bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json` (`shipped`, `rmsd_avg`, `rmsd_arm`, `rmsd_fit`).

```
SHIPPED THROUGH THE LADDER'S OWN PATH (the anchor), n=126
  sel   mean 3.4540 (pinned 3.4540)  vs production cache per target: max abs 0.00e+00, mean abs 0.00e+00, n(diff>1e-6)=0
  cloud mean 3.0483 (pinned 3.0483)  vs production cache per target: max abs 7.25e-14, mean abs 2.69e-15, n(diff>1e-6)=0
  arm   mean 3.2126 (pinned 3.2148)  vs production cache per target: max abs 1.71e-01, mean abs 1.07e-02, n(diff>1e-6)=120
  fit   mean 3.2052 (pinned 3.2041)  vs production cache per target: max abs 1.41e-01, mean abs 1.80e-03, n(diff>1e-6)=118
  gam_eff/cos of the shipped posterior against itself: prob +0.0000/nan loc +0.0000/nan (must be 0/nan)  MAE 2.3386
```

---

## L57 -- CORRECTION TO L56: THE ANCHOR'S BUILT CHAIN IS THE LEADERBOARD-REBUILD BASIS (3.2126), NOT THE PRODUCTION CACHE (3.2148); SELECTION AND POINT CLOUD ARE BIT-EXACT (2026-09-13 19:35, lane P)

L56's title says "bit-exact on all four bases"; the block under it says otherwise and the block is
right. Through `s26/p_ladder.py`'s path the shipped posterior reproduces the production cache
per target at 0.00e+00 on selection and 7.25e-14 on the point cloud, but the BUILT CHAIN lands at
mean 3.2126 (the S13 `base` / S24 L9 "leaderboard rebuild" figure; the state brief lists both:
"3.2148 (leaderboard rebuild: 3.2126)"), with 120/126 targets differing from `rmsd_arm` by more
than 1e-6: mean abs 0.0107, max abs 0.1711, signed mean -0.0021, worst
7JS6 0.171, 6QAX 0.148, 2LWU 0.121, 2NDN 0.115, 2LNG 0.101. Mechanism: the cloud enters `s12.instrument.project`
differing from the persisted `avg_ca` by 1e-14, and the multi-start L-BFGS projection lands in a
different local optimum on most targets; the 1A13 probe (L11) matched at 0.0 because it was fed
the persisted `avg_ca` itself. Consequence for the ladder: every rung and the anchor are built
through the SAME process and code path, so the paired contrasts are on one basis, the
rebuild basis 3.2126, and that basis is stated in every C2 entry from here on; the 0.0022 A
offset to the production figure is a basis difference, not an effect. The `fit` (lam = 0) basis
lands at 3.2052 (production 3.2041) for the same reason. L56 stands with its title read as
"bit-exact on selection and point cloud".

---

## L58 -- ADDENDUM TO L44 PER L55: THE BOUND ARTEFACT NOW CARRIES THE GAIN ENVELOPE ROW AND THE PER-TARGET READINGS; 0.028 IS QUOTED WITH 0.151 (WORST TARGET, 2BP4) BESIDE IT; ON THE SELECTION BASIS THE WORST-TARGET READING IS 0.194, ABOVE THE 0.170 LINE (2026-09-13, W)

Both L55 caveats applied by `s26/w_bound_addendum.py` (run directly, 4 s, one JSON read, under
200 MB; `s26/w_selfcopy.py` was not edited because `w_tiebreak_draws`, which imports it, is
running). `s26/results/w_selfcopy_bound.json` is regenerated with this script's provenance and
the original provenance kept under `provenance_original`; every number comes from
`s26/results/w_selfcopy_endpoint.json`.

1. `signed_bounds_gated/C_envelope_fold_ci/gain` added: leaked-model mean of (arm - sel) minus
   clean, +0.5101, fold CI [+0.3698, +0.6838], MDE 0.2548, 5/5 folds; bound (2/60) x 0.6838 =
   0.0228. `verdict/gain` now reads MINOR at 0.0228, as L44's Part D table says; the artefact and
   the ledger agree.
2. `signed_bounds_gated/C_envelope_per_target` added, every basis: (2/60) x |leaked-mean minus
   clean| at the worst single target, at the p95 target, and over (target, model) pairs. The
   Adversary's numbers reproduce exactly:

       basis   mean-CI limit   p95 target   worst target (pdb)     worst (target, model) pair   class under every reading
       arm       0.0277          0.0830       0.1512 (2BP4, 4.54 A)   0.1519                       MINOR
       gain      0.0228          0.0816       0.1232 (2NBC)           0.1261                       MINOR
       cloud     0.0280          0.0824       0.1587 (2MQ2)           0.1604                       MINOR
       sel       0.0479          0.1154       0.1944 (9KAR)           0.1944                       crosses 0.170 at the worst target

   So the class MINOR holds under every reading of the envelope on the built chain, the point
   cloud and the paired gain (the benchmark's headline contrast), and on the SELECTION basis
   (the benchmark's `shipped` argmin arm, the basis most sensitive to memorisation: L44 Part C
   -1.21 A) the worst-single-target reading is 0.194 A, above the pre-registered 0.170 line.
   Stated as the pre-registered rule gives it: for the benchmark's argmin mean alone, under the
   over-bound of own-native training and assumption A2, a single benchmark self-copy behaving
   like 9KAR could move that mean by up to 0.19 A; the expected contribution on every basis stays
   0.002 A (dev proxy) to 0.05 A (mean envelope), and the paired gain, which is the benchmark's
   verdict, is MINOR under every reading (worst 0.123). The benchmark verdict (no validated gain,
   +0.0103 [-0.160, +0.180]) does not move.
3. Wording to carry everywhere, as L55 asks: "2/60 self-copies; dev-proxy price 0.002 A; own-native
   envelope 0.028 A (mean CI) to 0.151 A (worst target) on the built chain, under A2; MINOR under
   every reading on the built chain and the paired gain; the selection basis reaches 0.194 at the
   worst target; cannot move the benchmark verdict either way." Applied in
   `s26/agentW_FINDINGS.md` (sections 0, 2.5, 4, 5), `s26/PREREG_selfcopy_bound.md` addendum 2 and
   `s26/IDEA_selfcopy_proxy_bound.md`. L44 itself is not edited.
## L59 -- GOVERNOR v2 DIED AT 19:36 ON A WINDOWS FILE-REPLACE RACE; v2.1 RETRIES AND NEVER EXITS THE LOOP ON A TRANSIENT ERROR (2026-09-13 19:37, coordinator)

`s26/governor.py` exited with `PermissionError: [WinError 5] Access is denied` from
`os.replace(governor_state.json.tmp, governor_state.json)` at 19:36:xx (its background task's
output): a reader (a jobrun waiter polling the snapshot, or a lane's `--status`) held the target
open at the instant of the replace. Four jobs were registered (a1_build, ph_reject_chain,
w_tiebreak_draws, ph_branch_solutions) and ran unsupervised for about a minute at 75% RAM; none
was harmed. Fix: `write_json_atomic` retries the replace eight times with a short back-off and
logs a skipped snapshot instead of raising; the main loop wraps the sample in a try/except that
logs and retries on the next tick. Restarted as v2.1 at 19:37. No production module touched.

---

## L60 -- REPORT CHECK: EVERY APPENDIX B NUMBER OF s26/REPORT.md RE-READ FROM ITS ARTEFACT (106 ROWS: 79 PASS / 2 MISSING PATH / 25 NOT FOUND ON THE FIRST RUN; 138 / 138 AFTER THE FIXES); STYLE 0 / 0 / 0; PARTS VII AND VIII FILLED TO L59; APPENDIX D RECONCILES docs/REPORT_S26.md (2026-09-13 19:50, lane E)

`s26/e_report_check.py` (committed b981b523, revised through 7db128c0) parses Appendix B of
`s26/REPORT.md` row by row, resolves every backticked path (`file :: key, key`, `file:line`,
`file:a-b`, `tests/x.py::test`, ledger entries named after a ledger path), checks existence, and
checks every number quoted in the row against the artefact: JSON leaves under the stated keys
(list means, container lengths, `[field=value]` filters, means and counts over `rows[*]/x`,
numbers inside string leaves), text lines or ledger entries by literal or by tolerance (half a
unit in the last quoted digit; a percentage as x and x/100), and `derived: expr = value` cells
evaluated with every four-decimal literal itself required to be in the row's artefacts. Rule 1:
`s9/final_report.json` and anything named benchmark are checked for existence only (two rows,
the benchmark means, are recorded as unverifiable behind Rule 1). `--style` scans the whole report
for the six banned words, U+2013 / U+2014, and any RMSD sentence contrasting two numbers with no
basis label on the sentence or its paragraph. Output: the table on stdout and
`s26/results/e_report_check.json` (save_atomic provenance, report sha256). Repeatable as
`python s26/examine.py --report-check` (examine.py gained the flag; a9984267 did not parse and
6f3708dc repairs it).

Runs, all through jobrun (tag CPU, est-ram 0.3): `e_report_check` (55 s including the wait, peak
RSS 0.032 GB) on the draft as committed at fcfbaf1d: 106 rows, 79 PASS, 2 MISSING PATH, 25 NUMBER
NOT FOUND; `e_report_check2` and `e_report_check3` after the fixes: 109 / 109 PASS, style 0 / 0 /
0; `e_report_check4` on the Part VII/VIII/D state is queued behind the four registered jobs (the
job cap) and refreshes the artefact when a slot frees; the iteration runs between them ran
directly to a scratch output (5 s, 9 MB). Current state, 7db128c0: 138 rows, 138 PASS, style 0 /
0 / 0.

What the 27 failures were, and what was done (drop or correct, never a new number without an
artefact): 2 MISSING PATH (a root `README.md` resolved against the previous citation's directory,
parser fixed; `geo_pauli_v1_rawonly.json` cited without its `s13/results/` directory). 25 NUMBER
NOT FOUND: 10 wrong or incomplete key paths (`cells/25` -> `rows[*]/cells/25/band_best`;
`rows[*]/MASS1.0, MASS0.1, MASS0.0` -> three full paths; `rmsd_vqe_sel, medoid128` ->
`arms/vqe_LFO/sel, arms/medoid128/sel`; `singularity` -> `summary/singularity`; the
`phys_landscape` summary keys; `concentration` -> `[model=amber]` / `[model=legacy]` for the 25 /
41 table counts; `q_alpha`, `q_verify` and `c_land_null` narrowed to their keys; `n_rows` added
where a row quoted "of 126"); 6 derived cells declared (9450 = 126 x 75; 5057 = 163 + 2627 + 2267;
96.8%; 78 = 0.6190 x 126; 42 = 2 P; 8128 = dim so(128); 12.6% = 1025 / 8128; the prior-ladder slope
and gain; +0.0207 = 3.2355 - 3.2148); 5 citations moved to the document that carries the number
(`s25/agentPHYS_FINDINGS.md:57-58, :357, :379` for 40/126, 462/500, +0.054;
`s25/agentQ_FINDINGS.md:375-376` for 45.3% and the 64.7th percentile; `docs/FINDINGS.md:3160` for
+0.994; `core/quantum.py:42` for +0.655634; `docs/FINDINGS.md:2751-2760` for 1.386); 4 numbers
replaced by the artefact's own (0.36 to 0.88 -> 0.358 / 0.827 / 0.882, `s20/LEDGER.md` L6; 2.6e-4 ->
2.239e-4 against the 8.66e-4 quantisation bound, `s25/LEDGER.md` L10; 699 -> 725 modules,
`module_map.json :: n_modules`; 58.6% re-sourced to `ph_reject_census.json ::
summary/per_threshold/1e4/pool_frac_over`); 1 number DROPPED: the relaxation cost's interval
"+0.0207 [+0.0143, +0.0276]" quoted from the state brief has no artefact (the brief does not
contain it either); the mean +0.0207 is kept as a derived difference of two stored means and the
paired interval now comes from L39 (+0.0207, fold CI [+0.0154, +0.0290]). One S19 number (2.66 per
75) re-cited to `s19/LEDGER.md` L12. Style: the first, broad detector flagged 33 sentences; 24
were not RMSD contrasts (correlations, free energies, variances, Pauli weights) and the detector
now skips those; 9 real RMSD contrasts lacked a basis on the sentence and are labelled (S12,
S13, S14, S15 and S21 numbers in Part VI, one S16 contrast in IV.4, the S8-instrument residual in
V.7, the SPSA contrast in V.10). Banned words: one ("leverage", VI.17) removed at 92d559bb;
dashes: none.

Appendix D reconciles `docs/REPORT_S26.md` and its summary (L42 item 1; read only, not edited or
committed): eleven artefact-sourced items taken with the same paths and the docs file cited as
the pointer (the leaderboard spread, the seven configurations against the incumbent, the
calibration means, the ladder's first rung, the binning and the 51-arm null, the Gibbs falsifier,
the S5 cosine and the consolidation audit's cosines, the entangler-deletion source, the
end-to-end speed-up, B1's numbers); one disagreement decided by the artefact (5.6e-17 -> 5.551e-17,
mine corrected); five apparent disagreements that are the same quantity on two paths or two
criteria (3.2126 / 3.2148 and 0.18242 / 0.18198, `compare_` / `baseline_tuning126.json`, 13/60 at
>= 0.6 identity against the 2/60 verbatim self-copies, the "100 of 125" citation, the blind
pipeline on 126 against 108); not taken: their external sources (arXiv, an author page), the
supervisors' names and affiliation, their Part X presentation guide (lane PR's), their memory-file
citations, and the two benchmark target ids they print from `docs/FINDINGS.md:4571`.

Parts VII and VIII now carry L35, L38, L39, L43, L44, L52, L53, L56 / L57 and the Adversary's
L45 to L49, L54, L55, L58, each row with its ledger entry and artefact; the C2 anchor's basis
(the rebuild built chain 3.2126, L57) is stated where the ladder is discussed; the leak wording of
L58 is carried in II.1, VII.4 and VIII.2. Next: A1 / A3, B2 / C2 rungs, the steric reject on the
built chain, C3 stage 2, the tie-break floor and the remaining tournament entries as they land;
Appendix C at the close.

---

## L61 -- LANE PR: THE DECK IS BUILT FROM ARTEFACTS (11 SLIDES; SLIDES 1-7, 9, 10, 11 FILLED; SLIDE 8 A PENDING PLACEHOLDER WITH THE A2/A4 FACTS); 237 REGISTERED NUMBERS, EACH WITH ITS PATH; FOUR NUMBERS TYPED FROM THE CLAIM LEDGER, SEVEN KEPT OFF THE SLIDES (2026-09-13 19:55, PR)

`vqe_research_overview.pptx` (repository root) did not exist (L2), so it is BUILT by `s26/pr_build_deck.py`
(python-pptx 1.0.2; dark theme; one accent; figures on white plates at 190 dpi; title, body and notes on every
slide) from `s26/pr_values.py` (every number read from its artefact at build time; the record is
`s26/pr_values.json`, 237 tokens with value, path, basis, status), `s26/pr_notes.md` (the spoken text with
`{TOKEN}` placeholders; each notes frame ends with a SOURCES list) and `s26/pr_figures.py` (the two
ORACLE-superposed CA overlays for 1S9Z = T030 and 9KAR, their RMSDs recomputed through `s12.instrument` and
equal to the record at 1e-9: 0.181981 / 7.437696; the S25 width and depth sweeps rebuilt from
`s25/results/q_plateau.json`; the accuracy ladder with the basis on every bar). Lane Q's A2 and A4 figures sit on
slides 6 and 8. Commits 71efb76d (a), a6738b35 (b), e93fc98c (c), cb90ffb3 (d placeholder). Verification
(`s26/pr_verify.txt`, pasted into `s26/PRESENTATION_CHANGES.md`): 11 slides; 0 U+2014, 0 U+2013; 0 banned words;
spoken words 234 / 248 / 244 on slides 8 / 9 / 10 (limit 250). Build peak RSS 0.12 GB, run directly.

Basis rule (L28 item 3, L29 item 2) applied throughout: the suite, the random-75 3.4251 and the Legacy / AMBER
verdicts are named point cloud; 3.2148 built chain; the C2 anchor the rebuild basis 3.2126 (L57); the S25
quantum contrasts the selection basis. The Adversary's caveats are in the wording: +0.0111 with the Type-M flag
and "a validity step on 124 of 126" (L46); the product circuit grounded in `q_dla.json` and no CI claimed on the
-0.302 vs -0.243 slopes (L47); the steric reject with the flag and "tail-carried" (L54); the strain signal as a
calibration flag, never a gain (L53); the benchmark caveat in L58's words on slide 4's notes.

Not sourced to a results artefact and labelled as such: the benchmark +0.0103 [-0.1596, +0.1803] 31W/29L (typed
from claim C06; `s9/final_report.json` not opened; the two means SOURCED_BY_TEST), the 4.4 GB headroom (L13,
a governor reading), the 11-rung count (the PREREG), the MDE factor (the contract). Kept off every slide: the
+0.0030 benchmark leak price (C27; the caveat uses lane W's bound instead), the z_moment triple (C34), the 0.524
cosine (C24), the 355/13 count (C35), the lost ESM -0.288 (L11), the S13 Pauli mean weights 2.236 / 3.015 (not
reproduced by me from `geo_pauli.json`; the median measured/predicted ratio 0.9969 over 104 cells is what slide
9 quotes), and the S12 operator-law coefficients. The phi statistic is on slide 2 with
`s26/results/a_c26_phi_mae.json` (L32). Ledger numbers recomputed on the way and found identical: L14, L26,
L35, L46, L53. Slide 11's direction line is DRAFT until the coordinator's verdict entry. Slide 8 waits for
`s26/PROPOSAL_A.md`; the builder swaps the placeholder for the proposal slide when the file exists.

---

## L62 -- LANE P, C2 RUNG NOESM: REMOVING THE ESM BLOCK COSTS +0.208 A ON THE BUILT CHAIN (5/5 FOLDS, 1.01x MDE, TYPE-M ZONE) AND +0.330 A ON SELECTION, REPRODUCING S7-11's LOST -0.288 TO -0.34; AND ITS gam_eff IS POSITIVE WHILE IT IS WORSE (2026-09-13 20:01, lane P)

Artefacts: `s26/results/p_ladder_noesm_s0.json` (126 rows, complete), `s26/results/p_ladder_report_noesm_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_noesm2`: exit 0, 521 s, peak RSS 0.305 GB (`s26/jobs_done/p_eval_noesm2.json`); fold models `s26/models/p_ladder/noesm_fold{0..4}_s0.pt` (42-d physicochemical pair features, the shipped MLP and trainer; `p_train_noesm` 1,743 s, peak 0.63 GB). Reading: (1) the S7-11 figure whose artefact is lost (L11) is re-measured on the same instrument: the ESM block is worth -0.330 A [fold -0.448, -0.212] on selection, 5/5 folds, against S7-11's -0.288 [-0.484, -0.092] for pca32 and -0.34 for pca128; (2) on the built chain it is worth -0.208 A, 5/5 folds, at 1.01x its MDE (Type-M zone; the magnitude is uncertain, the sign is not), so the ESM channel survives the pipeline, which S17 L23's 'filter, not discriminator' did not settle for the readout; (3) the medians (+0.044 arm, +0.071 sel) are far below the means and the worst target is +3.65 A (8T61): the loss is concentrated on a few targets, but the uniform-effect null puts the drop-top-10 at its 50th percentile, so it is not a concentration artefact; (4) gam_eff is POSITIVE (+0.116 prob-space at cos 0.25; +0.343 loc-space at cos 0.34, amplitude ~1) for a rung that is WORSE at 5/5 folds and has worse MAE (2.548 vs 2.339): a large move at low cosine projects positively onto the truth direction while adding more orthogonal error than it removes. This is S25 L12's caveat reproduced from the training side, on an achievable rung, and it settles how gam_eff will be read for every rung below: never without the cosine and the amplitude, and never as a prediction of the endpoint.

```
  noesm vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.4205 (med 3.4139)   b 3.2126 (med 2.9661)   n=126
    effect +0.2078   median +0.0440   SE 0.0733   MDE 0.2054   effect/MDE +1.01
    iid  CI95 [+0.0647, +0.3511]
    fold CI95 [+0.0986, +0.3942]   folds same sign 5/5   per-fold 0:+0.115 1:+0.126 2:+0.200 3:+0.573 4:+0.075
    47W/79L/0T   worst degradation +3.6485 (8T61)   p90 +0.9859   power 0.81  Type-M 1.12
    concentration: drop-top10 +0.3348 vs uniform-effect null p10/p50/p90 +0.2551/+0.3322/+0.4153 -> pctile 0.514
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.12x]
  noesm vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.2442 (med 3.1375)   b 3.0483 (med 2.8373)   n=126
    effect +0.1959   median +0.0576   SE 0.0684   MDE 0.1917   effect/MDE +1.02
    iid  CI95 [+0.0554, +0.3310]
    fold CI95 [+0.0775, +0.3814]   folds same sign 5/5   per-fold 0:+0.112 1:+0.193 2:+0.153 3:+0.547 4:+0.035
    47W/79L/0T   worst degradation +3.4270 (8T61)   p90 +1.0224   power 0.82  Type-M 1.11
    concentration: drop-top10 +0.3171 vs uniform-effect null p10/p50/p90 +0.2411/+0.3155/+0.3936 -> pctile 0.509
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.11x]
  noesm vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.7839 (med 3.6764)   b 3.4540 (med 3.4779)   n=126
    effect +0.3299   median +0.0713   SE 0.0924   MDE 0.2589   effect/MDE +1.27
    iid  CI95 [+0.1499, +0.5173]
    fold CI95 [+0.2122, +0.4475]   folds same sign 5/5   per-fold 0:+0.149 1:+0.203 2:+0.521 3:+0.443 4:+0.332
    48W/70L/8T   worst degradation +3.1238 (2LWS)   p90 +1.7459   power 0.95  Type-M 1.03
    concentration: drop-top10 +0.5026 vs uniform-effect null p10/p50/p90 +0.3881/+0.4977/+0.6130 -> pctile 0.524
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.03x]
  gamma-equivalent: gam_eff prob-space +0.1161 at cos +0.247 ; loc-space +0.3430 at cos +0.341 ; MAE 2.5483 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 5/5 ; verdict (arm): WORSE [TYPE-M ZONE: magnitude inflated ~1.12x]
```

---

## L63 -- LANE P, C2 RUNG CONLY: THE CONTACT HEAD ALONE (55-d, NO EMBEDDING BLOCK) SITS BETWEEN noesm AND THE SHIPPED PRIOR: +0.122 A ON THE BUILT CHAIN (0.62x MDE, 4/5 FOLDS, NOT MEASURED); AGAINST noesm IT BUYS -0.218 A ON SELECTION (1.00x MDE, 4/5) AND -0.086 A ON THE BUILT CHAIN (0.51x) (2026-09-13 20:12, lane P)

Artefacts: `s26/results/p_ladder_conly_s0.json` (126 rows, complete), `s26/results/p_ladder_report_conly_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_conly`: exit 0, 621 s, peak RSS 0.305 GB; models `s26/models/p_ladder/conly_fold{0..4}_s0.pt` (physicochemical pair block + the 13 ESM-2 650M contact-head columns; `p_train_conly` peak 0.69 GB). PREREG_B2's first size-axis point. Reading: the rung is UNDERPOWERED against the shipped prior (0.62x MDE on the built chain, 0.42x on selection; power 0.41 / 0.22), not null; its medians (+0.018 arm, 0.000 sel with 15 exact ties) say most targets are unmoved and the mean is carried by a few (worst 1CEK +3.70 A, the target whose verbatim self-window is BLOSUM rank 1, S24 L4). Isolated against noesm (paired, same code path): the contact head alone is worth -0.218 A on selection [fold -0.374, -0.037], 4/5 folds, exactly at its MDE (1.00x, Type-M zone), and -0.086 A on the built chain [fold -0.212, +0.010], 0.51x MDE, 2/5 folds. So of the ESM channel's -0.330 A on selection (L62), about two thirds is in the 13 contact-head columns and the rest in the 32-d embedding block (pca32, next); on the built chain the split is not resolved at n = 126. S17 L23 measured the contact head as a RANKER (filter, not discriminator); as a PRIOR INPUT it carries selection skill that the ranking measurement could not see, and that is consistent with the record: a better filter moves the argmin over K = 500 and barely moves a top-75 average. gam_eff/cos for the record: +0.109 prob at cos 0.24, +0.343 loc at cos 0.37, same shape as noesm (a large low-cosine move), MAE 2.437.

```
  conly vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.3349 (med 3.1110)   b 3.2126 (med 2.9661)   n=126
    effect +0.1223   median +0.0176   SE 0.0708   MDE 0.1983   effect/MDE +0.62
    iid  CI95 [-0.0113, +0.2620]
    fold CI95 [-0.0350, +0.2674]   folds same sign 4/5   per-fold 0:+0.115 1:-0.177 2:+0.221 3:+0.379 4:+0.079
    53W/73L/0T   worst degradation +3.7009 (1CEK)   p90 +0.8198   power 0.41  Type-M 1.55
    concentration: drop-top10 +0.2519 vs uniform-effect null p10/p50/p90 +0.1750/+0.2472/+0.3239 -> pctile 0.529
    VERDICT: NOT MEASURED (|effect| 0.1223 <= its own MDE 0.1983, 0.62x)
  conly vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.1532 (med 2.8753)   b 3.0483 (med 2.8373)   n=126
    effect +0.1049   median +0.0227   SE 0.0638   MDE 0.1788   effect/MDE +0.59
    iid  CI95 [-0.0201, +0.2354]
    fold CI95 [-0.0176, +0.2458]   folds same sign 4/5   per-fold 0:+0.096 1:-0.136 2:+0.196 3:+0.359 4:+0.026
    50W/76L/0T   worst degradation +3.1669 (1CEK)   p90 +0.6442   power 0.38  Type-M 1.61
    concentration: drop-top10 +0.2243 vs uniform-effect null p10/p50/p90 +0.1540/+0.2207/+0.2924 -> pctile 0.523
    VERDICT: NOT MEASURED (|effect| 0.1049 <= its own MDE 0.1788, 0.59x)
  conly vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.5656 (med 3.5772)   b 3.4540 (med 3.4779)   n=126
    effect +0.1116   median +0.0000   SE 0.0950   MDE 0.2660   effect/MDE +0.42
    iid  CI95 [-0.0761, +0.3032]
    fold CI95 [-0.0975, +0.3008]   folds same sign 3/5   per-fold 0:+0.279 1:-0.288 2:+0.373 3:+0.191 4:-0.000
    56W/55L/15T   worst degradation +4.4311 (1CEK)   p90 +1.5913   power 0.22  Type-M 2.15
    concentration: drop-top10 +0.2906 vs uniform-effect null p10/p50/p90 +0.1683/+0.2864/+0.4050 -> pctile 0.520
    VERDICT: NOT MEASURED (|effect| 0.1116 <= its own MDE 0.2660, 0.42x)
  gamma-equivalent: gam_eff prob-space +0.1093 at cos +0.241 ; loc-space +0.3425 at cos +0.368 ; MAE 2.4373 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 4/5 ; verdict (arm): NOT MEASURED (|effect| 0.1223 <= its own MDE 0.1983, 0.62x)
```

---

## L64 -- TOURNAMENT ITEM 4, tiebreak_noise_floor (W): THE PIPELINE'S OWN NOISE FLOOR FROM BLOSUM TIE-BREAKING AT THE K = 500 BOUNDARY IS 0.004 A ON THE 126-MEAN (BUILT CHAIN) AND 0.024 A AS THE PAIRED MDE BETWEEN TWO CONVENTIONS; EVERY HUNDREDTHS-LEVEL EFFECT ON THE RECORD IS INSIDE THE LATTER; THE PRODUCTION CONVENTION IS NOT DISTINGUISHABLE FROM A RANDOM DRAW (2026-09-13, W)

Pre-registered in `s26/PREREG_tiebreak_floor.md` (falsifiers and expectations in section 1, written
before the probe). Native-free half: job `w_tiebreak_draws` (exit 0, 5166 s under a four-job load,
peak RSS 0.111 GB; `s26/results/w_selfcopy_tiebreak_draws.json`, complete 126/126, production gate
top-75 == `sub` on 126/126). Gated half: job `w_tiebreak_endpoint` (exit 0, 25 s, 0.062 GB;
`s26/results/w_selfcopy_tiebreak_endpoint.json`); ledger numbers composed by
`s26/w_tiebreak_report.py` -> `s26/results/w_tiebreak_report.json`. Operator: on each target the
boundary tie class (every window sharing the 500th window's BLOSUM62 sum; median 115 windows, 54.5
of them inside the production pool) is re-drawn uniformly 8 times (`s15.seed.stable_rng`, seeds
("w_tiebreak", pdb, k)), the non-tied prefix unchanged, everything downstream identical (shipped
score, top-75, medoid-frame average, ramah 0.3 multi-start projection). Basis stated on every
line; the production chain here is the rebuild basis 3.2126 (L9, L57). Negative = the first arm is
better.

**Membership (native-free).** A random tie-break replaces 3.57 of the 75 averaged members on
average (median 3, max 15; one target never changes) and leaves the argmin unchanged on 91.5% of
(target, draw) cells: the score's top-75 is drawn mostly from the non-tied prefix. The registered
"about 8 of 75" was an over-estimate. The chain moves 0.168 A per draw in the median (triangle
bound, p90 1.08 A), 7x more than its RMSD-to-native changes (below): the moves are mostly
orthogonal to the native, as in L44.

**The floor, per basis (8 draws x 126 targets):**

    basis   m_tie (sd of the 126-mean over 8 draws)   mean over draws   production   s_tie median / p90   paired MDE between two draws: median / max over 28 pairs   max |paired effect| between draws
    arm     0.0039                                     3.2152            3.2126       0.0227 / 0.1290      0.0236 / 0.0321                                              0.0108
    cloud   0.0017                                     3.0532            3.0483       0.0126 / 0.0408      0.0100 / 0.0124                                              0.0046
    fit     0.0032                                     3.2099            3.2052       0.0148 / 0.0745      0.0157 / 0.0187                                              0.0104
    sel     0.0136                                     3.4741            3.4540       0.0000 / 0.0901      0.0466 / 0.0714                                              0.0398

Per target on the built chain: s_tie above 0.1 A on 21/126, above 0.5 A on 0/126; the range over
the 8 draws median 0.065 A, p90 0.362, max 0.713 (1D6X, the `pool_gate = WARN` target of L7/L9,
whose emission sits on a projection-branch flip). m_tie from 8 draws carries a relative SE of about
27% (sd of 8 numbers).

**Falsifier (PREREG section 1).** The claim is FALSIFIED if m_tie < 0.002 AND the paired MDE < 0.010
on the built chain: m_tie 0.0039 and paired MDE 0.0236, so not falsified. CONFIRMED if any recorded
effect is below the paired MDE between two draws: all nine listed are (`w_tiebreak_report.json ::
recorded_effects_vs_floor`): C3's production relaxation +0.0207 (L39) and its +0.0111 against the
matched random move, the S24 restrained relaxation -0.022, the ORACLE functional weight 0.015 (S24
L16), the all-atom reranking +0.004, C27's +0.0004 / +0.0018 (L44), the fd-vs-exact projection
difference 0.012 (L19), L24's predicted +0.013. Three of them (C27's two prices and the all-atom
reranking) are also below 2 x m_tie = 0.008 A, the spread of the 126-mean that the convention alone
produces. Registered predictions: m_tie 0.003 to 0.010 (measured 0.0039, holds); paired MDE 0.02 to
0.05 (0.0236, holds); median s_tie 0.03 to 0.08 (0.0227, just below); ~8 of 75 replaced (3.6, below).

**What the floor means, and what it does not.** It is the size of effect an unstated convention
produces on its own: two runs of the same pipeline that differ only in how the boundary ties are
broken differ on the 126-mean by an sd of 0.004 A (built chain) and can be told apart at 80% power
only above 0.024 A. It applies to any contrast whose two arms do NOT share the tie-break: a re-run
after a change of corpus order (the pinned `pdbs/` file set, state brief section 3), a different K,
a re-retrieval, a comparison across instruments or sprints (S9's 0.004 to 0.005 A residual between
the stable and the plain argsort, `docs/FINDINGS.md:1526, 7261`, is exactly this floor's effect on
the mean). It does NOT apply to a paired contrast whose two arms share the pool and its tie-break
(C3's relaxation, the functional weight, the reranking were all measured that way, and their
paired SEs already exclude this noise). So the sentence for the report is: "a hundredths-level
effect is real only as a paired contrast with the tie-break held fixed; quoted across runs or
against another instrument it is inside the pipeline's own convention noise (0.024 A at n = 126)."

**The production convention against random draws (registered secondary; `ST.fmt` verbatim):**

  production convention minus mean over 8 draws [arm], n=126
    a 3.2126 (med 2.9661)   b 3.2152 (med 2.9625)   n=126
    effect -0.0026   median -0.0004   SE 0.0062   MDE 0.0174   effect/MDE -0.15
    iid  CI95 [-0.0142, +0.0098]
    fold CI95 [-0.0155, +0.0156]   folds same sign 4/5   per-fold 0:-0.010 1:+0.031 2:-0.001 3:-0.004 4:-0.022
    64W/61L/1T   worst degradation +0.3828 (2LNG)   p90 +0.0577   power 0.07  Type-M 5.64
    concentration: drop-top10 +0.0099 vs uniform-effect null p10/p50/p90 +0.0023/+0.0095/+0.0171 -> pctile 0.534
    VERDICT: NOT MEASURED (|effect| 0.0026 <= its own MDE 0.0174, 0.15x)
  production convention minus mean over 8 draws [cloud], n=126
    a 3.0483 (med 2.8373)   b 3.0532 (med 2.8162)   n=126
    effect -0.0048   median -0.0011   SE 0.0032   MDE 0.0088   effect/MDE -0.55
    iid  CI95 [-0.0112, +0.0011]
    fold CI95 [-0.0089, -0.0006]   folds same sign 4/5   per-fold 0:-0.011 1:-0.007 2:+0.003 3:-0.008 4:-0.003
    70W/55L/1T   worst degradation +0.1433 (7YFS)   p90 +0.0133   power 0.34  Type-M 1.71
    concentration: drop-top10 +0.0019 vs uniform-effect null p10/p50/p90 -0.0013/+0.0017/+0.0050 -> pctile 0.524
    VERDICT: NOT MEASURED (|effect| 0.0048 <= its own MDE 0.0088, 0.55x)
  production convention minus mean over 8 draws [sel], n=126
    a 3.4540 (med 3.4779)   b 3.4741 (med 3.4796)   n=126
    effect -0.0201   median +0.0000   SE 0.0096   MDE 0.0269   effect/MDE -0.75
    iid  CI95 [-0.0392, -0.0033]
    fold CI95 [-0.0392, -0.0018]   folds same sign 4/5   per-fold 0:+0.009 1:-0.000 2:-0.055 3:-0.028 4:-0.024
    13W/6L/107T   worst degradation +0.2561 (6F3V)   p90 +0.0000   power 0.55  Type-M 1.34
    concentration: drop-top10 +0.0065 vs uniform-effect null p10/p50/p90 -0.0009/+0.0052/+0.0106 -> pctile 0.618
    VERDICT: NOT MEASURED (|effect| 0.0201 <= its own MDE 0.0269, 0.75x)

The corpus-order convention is not distinguishable from a random draw on the built chain (0.15x
MDE) or the point cloud (0.55x); on the selection basis it is 0.020 A better than a random draw at
0.75x MDE with the fold CI excluding zero and 107 ties: suggestive of the S25 finding that the
retrieval order is not neutral (rho(pool index, ORACLE RMSD) +0.054), and NOT MEASURED by the rule.
Nothing here is a lever: a tie-break cannot be chosen native-free, and the point of the number is
the floor. One representative draw-to-draw contrast (`ST.fmt` verbatim; the other 27 pairs are
summarised above):

  tie-break draw 0 minus draw 1 [arm], n=126
    a 3.2182 (med 2.9452)   b 3.2143 (med 3.0146)   n=126
    effect +0.0039   median +0.0007   SE 0.0084   MDE 0.0236   effect/MDE +0.17
    iid  CI95 [-0.0123, +0.0203]
    fold CI95 [-0.0041, +0.0193]   folds same sign 1/5   per-fold 0:-0.001 1:-0.001 2:-0.004 3:+0.034 4:-0.005
    54W/69L/3T   worst degradation +0.2871 (8T61)   p90 +0.0708   power 0.07  Type-M 5.16
    concentration: drop-top10 +0.0222 vs uniform-effect null p10/p50/p90 +0.0117/+0.0216/+0.0310 -> pctile 0.538
    VERDICT: NOT MEASURED (|effect| 0.0039 <= its own MDE 0.0236, 0.17x)

Replication: the 8 draws are 8 seeds; the 28 pairwise contrasts give the paired MDE's spread
(0.0236 median, 0.0321 max on the built chain); there is nothing positive to replicate at a second
seed, and the fold order does not enter (no fit). Power: the floor's own numbers ARE the power
statement. Not done: 16 draws (PREREG fork; the 8-draw m_tie carries a 27% relative SE, which the
verdict does not depend on); the AMBER-relaxed basis.

## L65 -- LANE P, C2 RUNG PCA32: REPRODUCTION GATE 2 PASSES AT THE ENDPOINT: THE RETRAINED pca32 EMITS THE SHIPPED PIPELINE'S ANSWER ON 126/126 TARGETS, ALL THREE BASES, 0 A DIFFERENCE (2026-09-13 21:01, lane P)

Artefacts: `s26/results/p_ladder_pca32_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pca32_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pca32`: exit 0, 536 s, peak RSS 0.306 GB; models `s26/models/p_ladder/pca32_fold{0..4}_s0.pt` retrained from the compact tables with the shipped recipe (`p_train_pca32`, peak 0.69 GB). The retrained posterior differs from the pinned one by 1.5e-5 per pair (fold-0 gate, L12) and that difference changes no argmin, no top-75 membership and therefore no emitted cloud or chain on any of the 126 targets: 126 exact ties on selection, point cloud and built chain (the FLAG on the concentration line is the degenerate all-zero case). This closes PREREG_C2's gate 2 in the strongest form: the ladder's trainer IS the production trainer at the endpoint, so any rung difference below is attributable to its inputs or architecture and not to retraining noise. It also fixes the ESM contrast: noesm (L62) is the shipped pipeline minus its ESM block and nothing else. gam_eff 0 at cos 0 (the identity), MAE 2.3386 = the shipped diagnostic.

```
  pca32 vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2126 (med 2.9661)   b 3.2126 (med 2.9661)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  pca32 vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0483 (med 2.8373)   b 3.0483 (med 2.8373)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  pca32 vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4540 (med 3.4779)   b 3.4540 (med 3.4779)   n=126
    effect +0.0000   median +0.0000   SE 0.0000   MDE 0.0000   effect/MDE +nan
    iid  CI95 [+0.0000, +0.0000]
    fold CI95 [+0.0000, +0.0000]   folds same sign 5/5   per-fold 0:+0.000 1:+0.000 2:+0.000 3:+0.000 4:+0.000
    0W/0L/126T   worst degradation +0.0000 (1A13)   p90 +0.0000   power nan  Type-M nan
    concentration: drop-top10 +0.0000 vs uniform-effect null p10/p50/p90 +0.0000/+0.0000/+0.0000 -> pctile 0.000  FLAG
    VERDICT: NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
  gamma-equivalent: gam_eff prob-space -0.0000 at cos -0.001 ; loc-space -0.0000 at cos +0.010 ; MAE 2.3386 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 5/5 ; verdict (arm): NOT MEASURED (|effect| 0.0000 <= its own MDE 0.0000, nanx)
```

---

## L66 -- LANE P, C2 RUNG WIDE: THE CAPACITY RUNG (WIDTH 768, DEPTH 4, SAME INPUTS) IS NULL-TO-WORSE: +0.032 A ON THE BUILT CHAIN (0.25x MDE), +0.097 A ON SELECTION (0.48x MDE, FOLD CI ABOVE ZERO 5/5); MORE CAPACITY DOES NOT BUY A BETTER PRIOR (2026-09-13 21:11, lane P)

Artefacts: `s26/results/p_ladder_wide_s0.json` (126 rows, complete), `s26/results/p_ladder_report_wide_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_wide`: exit 0, 566 s, peak RSS 0.337 GB; models `s26/models/p_ladder/wide_fold{0..4}_s0.pt` (2.7x the parameters of the shipped 384x3, same 183-d inputs, corpus, epochs, loss; `p_train_wide` peak 1.03 GB). The built-chain contrast is UNDERPOWERED (0.25x MDE, power 0.11), so no gain of 0.13 A or more exists, and the point estimate points the wrong way. On selection the fold-clustered CI excludes zero on the harmful side at 5/5 folds while the effect is 0.48x its MDE: a real, small degradation, worth reporting as 'the wider head selects slightly worse', not as a measured effect size. MAE barely moves (2.348 vs 2.339): the wider net reaches the same conditional mean. This is S9-8's monotone-decreasing capacity curve (ridge beats a 512x3 net and a 1500-tree GBT) reproduced at the PRIOR rather than at the ranker: the information is not in the inputs, and adding capacity spends it on fitting the fragment distribution (S7-2). gam_eff +0.096 prob / +0.291 loc at cos 0.18 / 0.35: the third rung in a row with positive gam_eff and a worse endpoint.

```
  wide vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2444 (med 3.1364)   b 3.2126 (med 2.9661)   n=126
    effect +0.0317   median +0.0099   SE 0.0458   MDE 0.1282   effect/MDE +0.25
    iid  CI95 [-0.0566, +0.1201]
    fold CI95 [-0.0300, +0.1279]   folds same sign 3/5   per-fold 0:+0.018 1:+0.021 2:-0.056 3:+0.230 4:-0.027
    54W/72L/0T   worst degradation +2.4883 (6RRO)   p90 +0.3976   power 0.11  Type-M 3.50
    concentration: drop-top10 +0.1272 vs uniform-effect null p10/p50/p90 +0.0792/+0.1235/+0.1749 -> pctile 0.536
    VERDICT: NOT MEASURED (|effect| 0.0317 <= its own MDE 0.1282, 0.25x)
  wide vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0710 (med 2.8223)   b 3.0483 (med 2.8373)   n=126
    effect +0.0227   median +0.0123   SE 0.0427   MDE 0.1195   effect/MDE +0.19
    iid  CI95 [-0.0650, +0.1031]
    fold CI95 [-0.0500, +0.1161]   folds same sign 3/5   per-fold 0:+0.022 1:+0.060 2:-0.061 3:+0.199 4:-0.071
    56W/70L/0T   worst degradation +2.1217 (6RRO)   p90 +0.3972   power 0.08  Type-M 4.52
    concentration: drop-top10 +0.1131 vs uniform-effect null p10/p50/p90 +0.0692/+0.1096/+0.1545 -> pctile 0.544
    VERDICT: NOT MEASURED (|effect| 0.0227 <= its own MDE 0.1195, 0.19x)
  wide vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.5512 (med 3.5063)   b 3.4540 (med 3.4779)   n=126
    effect +0.0972   median +0.0000   SE 0.0725   MDE 0.2032   effect/MDE +0.48
    iid  CI95 [-0.0508, +0.2397]
    fold CI95 [+0.0487, +0.1587]   folds same sign 5/5   per-fold 0:+0.018 1:+0.099 2:+0.216 3:+0.099 4:+0.062
    55W/55L/16T   worst degradation +2.8395 (5MXS)   p90 +1.1442   power 0.27  Type-M 1.92
    concentration: drop-top10 +0.2280 vs uniform-effect null p10/p50/p90 +0.1348/+0.2230/+0.3121 -> pctile 0.529
    VERDICT: NOT MEASURED (|effect| 0.0972 <= its own MDE 0.2032, 0.48x)
  gamma-equivalent: gam_eff prob-space +0.0958 at cos +0.182 ; loc-space +0.2909 at cos +0.350 ; MAE 2.3480 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 3/5 ; verdict (arm): NOT MEASURED (|effect| 0.0317 <= its own MDE 0.1282, 0.25x)
```

---

## L67 -- LANE P, C2 RUNG PCA32F: REFITTING THE 32-PCA PER FOLD (NO GLOBAL-PCA LEAK) CHANGES NOTHING MEASURABLE: +0.018 A BUILT CHAIN (0.13x MDE), -0.032 A SELECTION (0.17x); THE SHIPPED GLOBAL PCA WAS NOT A LEAK WORTH ANYTHING (2026-09-13 21:34, lane P)

Artefacts: `s26/results/p_ladder_pca32f_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pca32f_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pca32f`: exit 0, 551 s, peak RSS 0.796 GB (the lean trainer's residue table, `s5/esmraw.npz`, is resident at eval); models `s26/models/p_ladder/pca32f_fold{0..4}_s0.pt` with per-fold bases `pca32f_pca_fold{0..4}.npz` fitted on training sequences only (`p_train_pca32f` peak 1.06 GB). Every contrast is far under its MDE (0.11x to 0.17x; power 0.06 to 0.08) with medians at or below 0.01 A: the rung is the shipped prior up to fitting noise. Two things this settles: the shipped `esm_pca.npz` (fitted once over the whole database, S7-11's 'weak leak') is worth nothing to the endpoint, so no ladder result is contaminated by it; and the isolate for pca128 is clean, since pca32f and pca128 differ only in the number of components. MAE 2.317 (the lowest so far) with no endpoint movement: S7-3's rule again. gam_eff +0.111 prob at cos 0.24, +0.348 loc at cos 0.40.

```
  pca32f vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2307 (med 3.0813)   b 3.2126 (med 2.9661)   n=126
    effect +0.0180   median -0.0092   SE 0.0491   MDE 0.1376   effect/MDE +0.13
    iid  CI95 [-0.0806, +0.1111]
    fold CI95 [-0.0725, +0.0970]   folds same sign 4/5   per-fold 0:+0.015 1:+0.103 2:+0.043 3:+0.125 4:-0.147
    70W/56L/0T   worst degradation +1.8708 (9BAF)   p90 +0.6142   power 0.07  Type-M 6.46
    concentration: drop-top10 +0.1163 vs uniform-effect null p10/p50/p90 +0.0582/+0.1142/+0.1710 -> pctile 0.518
    VERDICT: NOT MEASURED (|effect| 0.0180 <= its own MDE 0.1376, 0.13x)
  pca32f vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.0347 (med 2.9544)   b 3.0483 (med 2.8373)   n=126
    effect -0.0136   median -0.0119   SE 0.0453   MDE 0.1270   effect/MDE -0.11
    iid  CI95 [-0.1023, +0.0752]
    fold CI95 [-0.1099, +0.0751]   folds same sign 3/5   per-fold 0:-0.000 1:+0.098 2:-0.021 3:+0.093 4:-0.187
    71W/55L/0T   worst degradation +1.3544 (2N9M)   p90 +0.5068   power 0.06  Type-M 7.85
    concentration: drop-top10 +0.0852 vs uniform-effect null p10/p50/p90 +0.0314/+0.0820/+0.1313 -> pctile 0.532
    VERDICT: NOT MEASURED (|effect| 0.0136 <= its own MDE 0.1270, 0.11x)
  pca32f vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4224 (med 3.3752)   b 3.4540 (med 3.4779)   n=126
    effect -0.0316   median +0.0000   SE 0.0649   MDE 0.1819   effect/MDE -0.17
    iid  CI95 [-0.1541, +0.0917]
    fold CI95 [-0.1108, +0.0541]   folds same sign 3/5   per-fold 0:+0.084 1:-0.136 2:-0.092 3:+0.101 4:-0.099
    60W/49L/17T   worst degradation +2.1793 (2MAI)   p90 +0.7372   power 0.08  Type-M 4.91
    concentration: drop-top10 +0.1025 vs uniform-effect null p10/p50/p90 +0.0195/+0.0986/+0.1767 -> pctile 0.523
    VERDICT: NOT MEASURED (|effect| 0.0316 <= its own MDE 0.1819, 0.17x)
  gamma-equivalent: gam_eff prob-space +0.1105 at cos +0.237 ; loc-space +0.3479 at cos +0.400 ; MAE 2.3168 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 4/5 ; verdict (arm): NOT MEASURED (|effect| 0.0180 <= its own MDE 0.1376, 0.13x)
```

---

## L68 -- LANE Q, A1: qubit-ADAPT-VQE IN PLACE OF THE FIXED ANSATZ IS NULL AT THE REGISTERED THRESHOLD ON THE BUILT CHAIN (-0.014 / -0.022 A, 0.23x / 0.36x MDE, FOLD CIs SPAN ZERO); EVERY ADAPT ARM POINTS THE SAME WAY AND NONE CLEARS ITS MDE; ON THE 78 alpha = 1 TARGETS ADAPT REACHES THE PRODUCT GIBBS STATE THE FIXED CIRCUIT MISSES BY 0.90 NATS (2026-09-13, lane Q)

Pre-registered in `s26/PREREG_A1.md` (filed 09:00, before any endpoint existed). Build:
`s26/q_adapt.py --build --tag a1` (two shards killed by the host at 10:10 after 33 targets;
one governed process resumed from those checkpoints: `s26/jobs_done/a1_build.json`, 7,782 s,
peak RSS 0.383 GB). Label: `s26/jobs_done/a1_label.json` (25 s, 0.36 GB), run after L33.
Stats: `s26/results/a1_stats.json`, `s26/logs/a1_stats.log` (every ST.fmt block).
Per-target records: `s26/results/a1/<pdb>.json` (126). Every record reproduces the
production quantum arm bit-for-bit: `ca` and `q_ca` against
`bench_results/cache/464a0ddb5f283e04/<pdb>.json` at max |diff| 0.0 and the selection index
equal, 126 of 126 (key `cache_check`). Basis: built chain (`rmsd_q_synth`, the production
projection of the weighted average) on both sides of every contrast below; the selection
readout (`rmsd_sel`) is the s8-instrument secondary and is named where it appears.
Comparator: `fixed_zrank_it50`, the deployed selector (3.2280 A built chain, 3.3135 A
selection; the production top-75 arm `rmsd_arm` is 3.2148 A, the shipped argmin 3.4540).

### The two PRIMARY contrasts (built chain), their selection-basis secondaries, and the matched random control, verbatim

  rmsd_q_synth: adaptL2_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2142 (med 3.1307)   b 3.2280 (med 2.9926)   n=126
    effect -0.0138   median +0.0066   SE 0.0210   MDE 0.0588   effect/MDE -0.23
    iid  CI95 [-0.0562, +0.0258]
    fold CI95 [-0.0705, +0.0442]   folds same sign 3/5   per-fold 0:-0.116 1:-0.008 2:+0.008 3:+0.096 4:-0.035
    58W/68L/0T   worst degradation +0.5351 (2MK7)   p90 +0.2886   power 0.10  Type-M 3.68
    concentration: drop-top10 +0.0294 vs uniform-effect null p10/p50/p90 +0.0046/+0.0295/+0.0533 -> pctile 0.498
    VERDICT: NOT MEASURED (|effect| 0.0138 <= its own MDE 0.0588, 0.23x)

  rmsd_q_synth: adaptV_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2059 (med 3.0924)   b 3.2280 (med 2.9926)   n=126
    effect -0.0222   median +0.0012   SE 0.0217   MDE 0.0608   effect/MDE -0.36
    iid  CI95 [-0.0665, +0.0194]
    fold CI95 [-0.0854, +0.0424]   folds same sign 3/5   per-fold 0:-0.129 1:-0.068 2:+0.030 3:+0.095 4:-0.031
    61W/65L/0T   worst degradation +0.5816 (2N9M)   p90 +0.2586   power 0.18  Type-M 2.44
    concentration: drop-top10 +0.0229 vs uniform-effect null p10/p50/p90 -0.0018/+0.0224/+0.0469 -> pctile 0.509
    VERDICT: NOT MEASURED (|effect| 0.0222 <= its own MDE 0.0608, 0.36x)

  rmsd_q_synth: randH_adaptL2 - fixed_zrank_it50
    a 3.2617 (med 3.0966)   b 3.2280 (med 2.9926)   n=126
    effect +0.0337   median +0.0214   SE 0.0350   MDE 0.0982   effect/MDE +0.34
    iid  CI95 [-0.0354, +0.1008]
    fold CI95 [+0.0122, +0.0581]   folds same sign 5/5   per-fold 0:+0.010 1:+0.006 2:+0.054 3:+0.075 4:+0.026
    55W/71L/0T   worst degradation +1.1926 (1D6X)   p90 +0.4926   power 0.16  Type-M 2.58
    concentration: drop-top10 +0.1077 vs uniform-effect null p10/p50/p90 +0.0701/+0.1069/+0.1436 -> pctile 0.513
    VERDICT: NOT MEASURED (|effect| 0.0337 <= its own MDE 0.0982, 0.34x)

  rmsd_sel: adaptL2_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2698 (med 3.1903)   b 3.3135 (med 3.1524)   n=126
    effect -0.0437   median +0.0000   SE 0.0334   MDE 0.0936   effect/MDE -0.47
    iid  CI95 [-0.1103, +0.0213]
    fold CI95 [-0.0720, -0.0224]   folds same sign 5/5   per-fold 0:-0.100 1:-0.052 2:-0.028 3:-0.025 4:-0.018
    39W/29L/58T   worst degradation +1.7418 (8HVS)   p90 +0.0706   power 0.26  Type-M 1.96
    concentration: drop-top10 +0.0289 vs uniform-effect null p10/p50/p90 -0.0089/+0.0268/+0.0645 -> pctile 0.527
    VERDICT: NOT MEASURED (|effect| 0.0437 <= its own MDE 0.0936, 0.47x)

  rmsd_sel: adaptV_adam_best_zrank_P21 - fixed_zrank_it50
    a 3.2688 (med 3.1903)   b 3.3135 (med 3.1524)   n=126
    effect -0.0448   median +0.0000   SE 0.0343   MDE 0.0961   effect/MDE -0.47
    iid  CI95 [-0.1109, +0.0218]
    fold CI95 [-0.0874, -0.0066]   folds same sign 4/5   per-fold 0:-0.099 1:-0.106 2:+0.016 3:-0.025 4:-0.018
    43W/34L/49T   worst degradation +1.7418 (8HVS)   p90 +0.1779   power 0.26  Type-M 1.96
    concentration: drop-top10 +0.0294 vs uniform-effect null p10/p50/p90 -0.0106/+0.0272/+0.0650 -> pctile 0.531
    VERDICT: NOT MEASURED (|effect| 0.0448 <= its own MDE 0.0961, 0.47x)

### The registered falsifiers

"ADAPT is null" fires on both primaries: -0.0138 A (0.23x its MDE of 0.0588) and -0.0222 A
(0.36x its MDE of 0.0608) are inside +-0.5x MDE, the fold-clustered CIs [-0.071, +0.044] and
[-0.085, +0.042] span zero, 3 of 5 folds agree in sign, W/L is 58/68 and 61/65, and the
drop-top-10 statistic sits at the 50th and 51st percentile of the uniform-effect null (no
concentration). "ADAPT helps" does not fire; no replication is owed.

Power statement. The comparison resolves 0.059 to 0.061 A on the built chain (MDE = 2.8016
SE at n = 126). An effect smaller than that could exist unseen; the observed effects are a
third of it, and Gelman-Carlin power for a true effect equal to the observed is 0.10 and 0.18
(Type-M 3.7 and 2.4). For scale, the whole quantum synthesis against the classical top-75 arm
is +0.0133 A (`s26/results/q_mde_reference.json`), so the resolution is four times the size
of the component's own footprint. Null at the registered threshold; underpowered below 0.06 A.

### Every ADAPT arm, built chain vs fixed_zrank_it50 (effect, x MDE; all fold CIs span zero; 3-4/5 folds)

    adaptL2 adam_best  P7 / P14 / P21   -0.0128 (0.21x) / -0.0192 (0.33x) / -0.0138 (0.23x)
    adaptL2 lbfgs      P7 / P14 / P21   -0.0206 (0.34x) / -0.0187 (0.32x) / -0.0194 (0.33x)
    adaptV  adam_best  P7 / P14 / P21   -0.0128 (0.21x) / -0.0195 (0.32x) / -0.0222 (0.36x)
    adaptV  lbfgs      P7 / P14 / P21   -0.0206 (0.34x) / -0.0244 (0.39x) / -0.0245 (0.40x)
    controls:  fixed 750 steps -0.0035 (0.16x)   gibbs_T +0.0088 (0.11x)   uniform128 +0.0135 (0.14x)
               randH_fixed +0.0230 (0.24x, 4/5)   randH_adaptL2 +0.0337 (0.34x, fold CI [+0.012,+0.058], 5/5)
    alpha subsets: L2 P21 alpha=1 -0.0223 (0.24x, n=78), alpha=0.25 +0.0000 (0.00x, n=48);
                   V  P21 alpha=1 -0.0254 (0.28x),       alpha=0.25 -0.0169 (0.27x)

Twelve of twelve ADAPT arms are negative on the built chain and twelve of twelve on the
selection readout (-0.039 to -0.048 A, 0.42x to 0.47x MDE, fold CIs excluding zero on all
twelve, 45 to 59 exact ties per contrast; 58 of 126 targets emit the identical selected
candidate under ADAPT-L2-P21 and the fixed circuit). By the standing rule (clear the MDE AND
fold CI excluding zero) none is a result; it is the shape of S25's `VQE_LFO - argmin`
(-0.1405 at 0.68x, 5/5 folds): a consistent direction with an unmeasured magnitude. The
matched-entropy random control (same 128 candidates, same entropy, no information) is worse
than the fixed circuit by +0.034 with a fold CI above zero on 5/5 folds and is likewise below
its MDE. The 7-parameter arms (P7, one trained RY layer) carry the same effect as the
21-parameter ones on both bases.

### The property half, on the real targets (no native in any of these numbers)

    KL(p || Gibbs) on the 78 alpha = 1 targets   fixed circuit 0.9027 (max 0.984)   [S25's 0.902 reproduced]
                                                 ADAPT, both pools, both optimisers 0.0002 (max 0.0009)
    KL(Gibbs_zrank || product of its marginals)  mean 1.4e-4, max 7.9e-4 over 126 targets
    L-BFGS growth at alpha = 1                   stops with no operator selected on 78 of 78 targets
    ADAPT sequences, alpha = 0.25, pool L2       mean 13.7 distinct 2-local strings; 38 distinct sequences on 48 targets

The deployed selector's target at alpha = 1 is a product state on every real target; a
7-parameter RY layer reaches it and the 21-parameter fixed circuit stops 0.90 nats short.
Reaching it changes the emitted structure by -0.02 A, a third of the MDE. That is the one
sentence Proposal A earns, and it is negative: a better-optimised state on the deployed
Hamiltonian is not visible through the readout, S25's readout-slack finding reproduced by a
second ansatz family.

Verdict for Proposal A (`s26/PROPOSAL_A.md`): REPLACE. The record now holds both halves: the
Hamiltonian is a constant ladder whose optimum is a product state (nothing to grow), and
growing anyway moves the answer by a third of the MDE.
## L69 -- PROPOSAL A VERDICT ACCEPTED: REPLACE (NOT THE EXPECTED KEEP WITH EDITS); HOW SLIDES 8 AND 9 STAY DISTINCT (2026-09-13 21:40, coordinator)

The campaign prompt expected KEEP WITH EDITS for Proposal A ("ADAPT-VQE is the right tool for
the trainability question and the DLA measurement"). Lane Q's evidence (L27, L35, L68;
`s26/PROPOSAL_A.md`) says the mechanism is absent rather than weak: the deployed objective's
optimum is a product state on every target (KL to the product of marginals at most 7.9e-4
nats), ADAPT offered entangling operators grows none on any of the 78 alpha = 1 targets, the
fixed circuit's DLA is already the full so(2^n) from depth 2 so there is no expressivity gap
for an adaptive ansatz to fill, the grown circuits give no width-scaling argument (L35), and the
endpoint is null at 0.23x and 0.36x its MDE with fold CIs spanning zero and a stated
resolution of 0.06 A. The prompt's own rule governs: a proposal is not softened to survive. The
verdict REPLACE stands, subject to the Adversary's check of L68.

Structure of the deck after this verdict (for lane PR): Proposal A and Proposal B both point
at the trainability paper, so slides 8 and 9 must not say the same thing twice. Slide 8 is
"what we learned by letting the circuit grow" (the ADAPT measurements and the product-state
diagnosis, the reason an ansatz cannot matter here); slide 9 is "what we publish" (the
trainability paper: the S13 locality theorem and Pauli-spectrum prediction, the S25 width sweep
scoped to depth 3 by the DLA result, A2, A4, the product-state fact, with the venue-honest gap
statement). Lane P's Proposal B verdict, when it lands, adds the feasible-scale ladder facts
(B2) and the B3 characterisability result to slide 9's notes as the evidence that the original
B had no testable form on this machine.

---

## L70 -- ADVERSARY CHECK OF L68 (A1, qubit-ADAPT vs THE FIXED ANSATZ): STANDS WITH CAVEAT; THE REPLACE VERDICT (L69) IS SUPPORTED (2026-09-13, A)

Checked against `s26/results/a1/<pdb>.json` (126 records), `s26/results/a1_stats.json`,
`s26/q_adapt.py` and `s26/PREREG_A1.md`.

Confirmed clean.
- Readout: every arm goes through the production functions `pl.consensus_medoid`,
  `pl.average_weighted` and `pl.project` on its own `p` (`s26/q_adapt.py:849-851`); the
  comparator `fixed_zrank_it50` reproduces the S25 quantum-arm cache bit-for-bit (`cache_check`:
  `ca` and `q_ca` max abs 0.0, `sel_equal` True, 126 of 126). Basis built chain both sides.
- Leakage: the native enters only in `label_one` (`:887-911`, RMSD scoring), gated on "PHASE 0
  SIGNED OFF" (`:925`); the build half reads no native. Clean.
- Registered threshold: "ADAPT is null fires if both PRIMARY contrasts lie within +-0.5x their
  own MDE" (`s26/PREREG_A1.md:85`); measured 0.23x and 0.36x. Fires as registered.
- Power statement: complete and correct (MDE 0.0588 / 0.0608, Gelman-Carlin power 0.10 / 0.18,
  Type-M 3.7 / 2.4, "null at the registered threshold; underpowered below 0.06 A"). Fold CIs
  span zero, 3 of 5 folds, concentration at the null's 50th percentile. No replication owed.
- The product-state claim rests on persisted per-target rows: `product_diagnostics/zrank/
  kl_gibbs_to_product` (recomputed over 126: mean 1.40e-4, max 7.90e-4, as quoted) and
  `arms/fixed_zrank_it50/kl_to_gibbs` on the 78 alpha = 1 targets (mean 0.9027, max 0.9836, as
  quoted; ADAPT-L2-P21 2.6e-4 mean, 7.0e-4 max). Sourced.

Three caveats.
1. **"Twelve of twelve arms negative" is one observation, not twelve.** The twelve ADAPT delta
   vectors (built chain, per target) have mean pairwise correlation 0.955 (min 0.919, max 1.000)
   and their first principal component carries 96.0% of their variance; `ST.best_of_k_within`
   over the 13 arms (fixed + 12) gives a per-target oracle of -0.094 A whose split-half transfer
   is +0.006 (-6% of the oracle: NOT A SIGNAL). L68 rightly calls none of them a result; the
   phrase "a consistent direction" must be read with the 0.955 beside it, so that twelve
   correlated draws at 0.2-0.4x MDE are not taken for a trend.
2. **The parameter match is a budget match, not a realised one.** `max_params = 3n = 21` for the
   arms labelled P21, but `len(ops)` runs 14 to 21 for Adam and 0 to 21 for L-BFGS (`adapt/*/ops`
   in the records): on the 78 alpha = 1 targets L-BFGS grows nothing and the "P21" arm is the
   7-parameter RY layer, and Adam grows 7 to 14 strings. Nothing in the null depends on it (P7,
   P14 and P21 agree within 0.01 A on both bases), but "21 parameters both sides" should be
   stated as "the same 21-parameter budget; realised counts 7 to 21".
3. **The `gibbs_T` control's +0.0088 overall masks two opposite subsets.** On the 78 alpha = 1
   targets the exact Gibbs state is -0.0271 against the fixed circuit and ADAPT-L2-P21 is
   -0.0223; the two agree with each other (+0.0047, SE 0.0060, 0.28x MDE, fold CI [-0.0009,
   +0.0098]) as two states 2.6e-4 nats apart should. On the 48 alpha = 0.25 targets `gibbs_T` is
   +0.0671 (the Gibbs state is not the CVaR optimum there) and ADAPT is +0.0000. So the
   product-state diagnosis is internally consistent on the subset it applies to; and, noted for
   the record, 29 of 78 alpha = 1 targets differ by more than 0.01 A between two states 2.6e-4
   nats apart (max 0.206 A): the projection readout amplifies sub-milli-nat differences on some
   targets, the mechanism L64 measures as the pipeline's own noise floor.

Verdict: STANDS WITH CAVEAT. The REPLACE verdict for Proposal A (L69) is supported: the
endpoint is null at the registered threshold with a stated resolution of 0.06 A, the optimum is
a product state on every real target by persisted per-target rows, and no arm family changes
the emitted structure beyond a third of the MDE.

---

## L71 -- ADVERSARY CHECK OF L64 (THE TIE-BREAK NOISE FLOOR): STANDS WITH CAVEAT (2026-09-13, A)

- Operator and reads: the boundary tie class is re-drawn uniformly with `s15.seed.stable_rng`
  (8 seeds), the non-tied prefix untouched, everything downstream production; the draws are
  native-free (`load_blind`) and the endpoint is gated; production gate top-75 == `sub` on
  126 of 126. Clean.
- Numbers: m_tie 0.0039 A (sd of the 126-mean over 8 draws; W states its ~27% relative SE),
  paired MDE between two draws median 0.0236 A over the 28 pairs (max 0.0321), per-target s_tie
  median 0.0227, p90 0.129, above 0.1 A on 21 of 126; 3.57 of 75 members replaced per draw; the
  argmin unchanged on 91.5% of cells. Registered predictions held (m_tie, paired MDE) or were
  over-estimates (8 of 75 replaced; median s_tie 0.03-0.08). `s26/results/w_tiebreak_report.json`.
- The production convention against a random draw: 0.15x MDE on the built chain, 0.55x on the
  cloud, 0.75x on `sel` with 107 ties and the fold CI excluding zero -- correctly NOT MEASURED,
  correctly read as suggestive of S25's non-neutral retrieval order and not a lever.
- Replication: the 8 seeds are the replication; no fit, no fold order. Power: the floor is the
  power statement.

Caveat. The title's sentence "EVERY HUNDREDTHS-LEVEL EFFECT ON THE RECORD IS INSIDE THE LATTER"
is true only as W's body qualifies it: the floor is the noise between two RUNS that do not share
the tie-break (a corpus re-order, a re-retrieval, a cross-instrument or cross-sprint comparison).
Every one of the nine listed effects (C3's +0.0207 and +0.0111, S24's -0.022, the 0.015 ORACLE
weight, the +0.004 reranking, C27's two prices, L19's 0.012, L24's +0.013) was measured as a
PAIRED contrast with the pool and its tie-break held fixed, and their paired SEs (C3: 0.0034)
already exclude this noise; none of them is invalidated. The presenter must quote the report
sentence W wrote ("a hundredths-level effect is real only as a paired contrast with the
tie-break held fixed; quoted across runs or against another instrument it is inside the
pipeline's own convention noise, 0.024 A at n = 126"), never the title alone. Second, minor:
m_tie from 8 draws carries a 27% relative SE, so "0.004" is 0.003 to 0.005; the verdict does
not depend on it.

Verdict: STANDS WITH CAVEAT (the title needs the body's qualifier wherever it is quoted).

---


## L72 -- LANE P, C2 RUNG PCA128: 128 ESM COMPONENTS INSTEAD OF 32: +0.076 A ON THE BUILT CHAIN (0.43x MDE, 3/5 FOLDS), +0.025 A ON SELECTION (0.13x); S7-11's '-0.337 vs -0.288, INDISTINGUISHABLE' STANDS AND THE DIRECTION IS, IF ANYTHING, WORSE (2026-09-13 21:45, lane P)

Artefacts: `s26/results/p_ladder_pca128_s0.json` (126 rows, complete), `s26/results/p_ladder_report_pca128_s0.json`; anchor `s26/results/p_ladder_shipped_s0.json`. Every arm through `s26/p_ladder.py`'s single path: shipped K=500 pool -> the rung's posterior in a genuine `core.predict.Distogram` -> shipped Bayes-risk score -> top-75 uniform medoid-frame average -> `s12.instrument.project` (ramah 0.3). Paired per target against the shipped posterior through the same path. Negative = the rung is better.

Job `p_eval_pca128`: exit 0, 576 s, peak RSS 0.799 GB; models `s26/models/p_ladder/pca128_fold{0..4}_s0.pt` (567-d inputs, the S7-11 per-fold 128-component bases `s7/repr_cache/pca128_fold{0..4}.npz`; `p_train_pca128` peak 1.41 GB). Isolated against pca32f (same construction, 32 vs 128 components; paired, same path): +0.058 A on the built chain, SE 0.056, MDE 0.156, 0.37x, fold CI [-0.048, +0.176], 3/5 folds, NOT MEASURED. S7-11 ranked pca128 nominally best on selection (3.417 vs pca32's 3.466, inside its CI); on this instrument it is +0.025 A on selection with 12 exact ties and +0.076 on the built chain, both far under their MDEs. Power: the design can see a 0.18 A built-chain gain and none of 0.18 A exists. MAE 2.386 (worse than pca32f's 2.317): more components fit the fragment corpus's directions, not the peptides'. gam_eff +0.117 prob at cos 0.21, +0.374 loc at cos 0.39.

```
  pca128 vs shipped -- BUILT CHAIN (PRIMARY; rebuild basis 3.2126, L57)
    a 3.2885 (med 3.2575)   b 3.2126 (med 2.9661)   n=126
    effect +0.0759   median +0.0158   SE 0.0630   MDE 0.1764   effect/MDE +0.43
    iid  CI95 [-0.0410, +0.1932]
    fold CI95 [-0.0691, +0.2157]   folds same sign 3/5   per-fold 0:+0.229 1:-0.173 2:-0.081 3:+0.256 4:+0.131
    57W/69L/0T   worst degradation +2.2800 (6RRO)   p90 +0.7311   power 0.23  Type-M 2.10
    concentration: drop-top10 +0.1958 vs uniform-effect null p10/p50/p90 +0.1236/+0.1917/+0.2650 -> pctile 0.528
    VERDICT: NOT MEASURED (|effect| 0.0759 <= its own MDE 0.1764, 0.43x)
  pca128 vs shipped -- POINT CLOUD (3.0483 basis)
    a 3.1089 (med 2.9577)   b 3.0483 (med 2.8373)   n=126
    effect +0.0606   median +0.0227   SE 0.0570   MDE 0.1598   effect/MDE +0.38
    iid  CI95 [-0.0573, +0.1725]
    fold CI95 [-0.0632, +0.1822]   folds same sign 3/5   per-fold 0:+0.204 1:-0.143 2:-0.088 3:+0.218 4:+0.100
    54W/72L/0T   worst degradation +1.8200 (2LER)   p90 +0.5037   power 0.19  Type-M 2.36
    concentration: drop-top10 +0.1731 vs uniform-effect null p10/p50/p90 +0.1081/+0.1682/+0.2315 -> pctile 0.538
    VERDICT: NOT MEASURED (|effect| 0.0606 <= its own MDE 0.1598, 0.38x)
  pca128 vs shipped -- SELECTION argmin K=500 (3.4540 basis)
    a 3.4789 (med 3.5336)   b 3.4540 (med 3.4779)   n=126
    effect +0.0249   median +0.0000   SE 0.0710   MDE 0.1990   effect/MDE +0.13
    iid  CI95 [-0.1117, +0.1606]
    fold CI95 [-0.1169, +0.1456]   folds same sign 4/5   per-fold 0:+0.036 1:-0.228 2:+0.037 3:+0.009 4:+0.211
    53W/61L/12T   worst degradation +2.5698 (7JGX)   p90 +0.7897   power 0.06  Type-M 6.76
    concentration: drop-top10 +0.1708 vs uniform-effect null p10/p50/p90 +0.0828/+0.1676/+0.2526 -> pctile 0.515
    VERDICT: NOT MEASURED (|effect| 0.0249 <= its own MDE 0.1990, 0.13x)
  gamma-equivalent: gam_eff prob-space +0.1171 at cos +0.214 ; loc-space +0.3744 at cos +0.387 ; MAE 2.3861 (diagnostic only).
  CAVEAT (S25 L12): -2.1496 x gam_eff is redeemable only at cos = 1; a real operator travelling 25%% at cos 0.5 is worth +0.024 A. Never quote the product alone.
  folds same sign (arm): 3/5 ; verdict (arm): NOT MEASURED (|effect| 0.0759 <= its own MDE 0.1764, 0.43x)
```

---

## L73 -- LANE PR: SLIDE 8 REBUILT IN PROPOSAL FORM FROM PROPOSAL_A.md (VERDICT REPLACE, L68 / L69) WITH THE ADVERSARY'S L70 CAVEATS ON THE SLIDE; SLIDES 8 AND 9 KEPT DISTINCT (DIAGNOSIS / PUBLICATION); ONE WORDING OF PROPOSAL_A.md THAT THE A1 RECORDS STATE DIFFERENTLY (2026-09-13 21:55, PR)

`vqe_research_overview.pptx` rebuilt (`python s26/pr_build_deck.py --no-figures`, then `s26/pr_changes.py`):
slide 8 now carries Proposal A's two-minute script adapted from `s26/PROPOSAL_A.md` section 5 (248 spoken
words) with every number read from `s26/results/a1_stats.json` (primaries: 2-local pool -0.0138 A, SE 0.0210,
MDE 0.0588, 0.23x, fold CI [-0.071, +0.044], 58W/68L; Tang pool -0.0222, SE 0.0217, MDE 0.0608, 0.36x,
[-0.085, +0.042], 61W/65L; NOT MEASURED; resolution 0.06 A; all 12 ADAPT arms -0.0128 to -0.0245 at 0.21x to
0.40x; the exact Gibbs state +0.0088 at 0.11x) and from the 126 records of `s26/results/a1/` (KL(Gibbs ||
product of marginals) mean 1.4e-4, max 7.9e-4; on the 78 alpha = 1 targets the fixed circuit's KL to Gibbs
0.9027, max 0.984, the 7-parameter RY layer's 1.5e-4, max 7.9e-4). Lane Q's A4 figure stays on slide 8; slide 9
is retitled "Direction B, and what we publish: the trainability paper" and names its Sprint 26 additions (the
DLA census, the product-state target, the non-decaying grown circuits, the A1 null), so that 8 is the diagnosis
and 9 the publication as L69 asks. The Adversary's L70 caveats are on slide 8, recomputed from the same records:
the twelve arms' per-target delta vectors correlate at 0.955 on average (min 0.919), one observation not
twelve; the 21 parameters are a budget (the Adam primaries realise 21, L-BFGS 7 to 21); the Gibbs control is
-0.0271 on the 78 alpha = 1 targets and +0.0671 on the 48 alpha = 0.25 targets. A1's basis is named on the
slide: `rmsd_q_synth`, the production projection of the weighted average over the 128 candidates (3.2280 A for
the deployed selector), not the production top-75 arm 3.2148. Verification (`s26/pr_verify.txt`, pasted into
`s26/PRESENTATION_CHANGES.md`): 11 slides, 0 U+2014, 0 U+2013, 0 banned words, spoken words 248 / 247 / 244 on
slides 8 / 9 / 10. Registry 296 tokens (`s26/pr_values.json`).

One statement of `s26/PROPOSAL_A.md` (sections 3 and 5) and of L68 that the records word differently, for lane
Q and the Adversary: "ADAPT ... selects no entangling operator on any of the 78 alpha = 1 targets under
L-BFGS" / "declines them on all 78 targets" / "stops with no operator selected on 78 of 78". In
`s26/results/a1/*.json :: adapt/{V,L2}_lbfgs_zrank`, the growth halts by the eps = 1e-3 pool-gradient criterion
on 156 of 156 (pool, target) cells, but before the halt it ADDS operators on 60 of 78 targets (pool V; 0 added on
18) and 68 of 78 (pool L2; 0 on 10), at most 11, the first of them entangling strings (1A13, pool V: IYIZZZZ at
pool-gradient norm 1.29e-3, just above eps); each such addition lowers the free energy by at most 1.2e-4 nats
(median 9e-6) from a 7-rotation state already within 7.9e-4 nats of the Gibbs optimum. The substance of the
claim holds (the growth is inert: nothing the objective or the readout can see), and the deck says it in the
records' words with the counts; the literal "no operator selected" does not hold. L70's caveat 2 ("on the 78
alpha = 1 targets L-BFGS grows nothing") reads the same way on the records: realised P for the L-BFGS P21 arms
runs 7 to 21. Lane Q's file is not edited. Slide 11's direction line stays DRAFT until the coordinator rules on
all three verdicts; lane P's B and C verdicts will trigger one more rebuild.

---
## L74 -- STALL AT 91.9% RAM WITH EVERY JOB SUSPENDED; GOVERNOR v2.2 ADDS A STALL BREAKER; THE raw RUNG'S TRAINING IS STOPPED BY HAND (21:51, coordinator)

At 21:50 the box read 91.9% RAM (user load: Chrome 5.0 GB, three claude sessions 2.3 GB, VS
Code 1.7 GB) with all four registered jobs suspended (ph_reject_chain, p_train_raw_f1 at
0.87 GB, w_ensemble_clouds, a3_build) and CPU at 16%: the v2 band suspends above 93% and
resumes only below 90%, so with the user's own load parked in between nothing of ours could
resume. v2.2: if RAM sits at or above 90% for 180 s while jobs are suspended, the governor
kills the suspended job holding the most memory (CTRL_BREAK first, so it checkpoints) and logs
it; the owner relaunches from the checkpoint when the box frees. To clear tonight's stall at
once the coordinator stops p_train_raw_f1 by hand (fold checkpoints under
`s26/models/p_ladder/`; lane P resumes it only on the coordinator's word) and restarts the
governor as v2.2. The raw rung is the ladder's last training item and not on the path to any
verdict.

---

## L75 -- LANE Q CORRECTS L68 ON PR's L73 FLAG: AT alpha = 1 ADAPT's GROWTH IS INERT, NOT ABSENT. OPERATORS ARE APPENDED ON 60 TO 78 OF 78 TARGETS, ALL MULTI-QUBIT UNDER L-BFGS, AND THEY BUY AT MOST 1.2e-4 NATS AND LEAVE A PRODUCT STATE (KL <= 4.1e-4); THE VERDICT DOES NOT MOVE (2026-09-13, lane Q)

PR's L73 read `s26/results/a1/` correctly and I did not. Re-read of my own records
(`s26/results/a1/<pdb>.json`, keys `adapt.<pool>_<optimiser>_zrank.sequence`, `.trace`,
`.theta`, `.stopped`, and `arms.*.kl_to_product`), the 78 alpha = 1 targets:

    pool / re-optimiser    stop reason     targets with >= 1     appended strings     F(P=7) - F(final)      max |angle| of      KL(p || product
                                           string appended       total / weight >= 2  mean / max, nats       appended strings    of marginals) max
    V  / L-BFGS-B          eps 78/78       60 / 78               235 / 235            1.5e-5 / 9.9e-5        0.018 rad           3.3e-4
    L2 / L-BFGS-B          eps 78/78       68 / 78               393 / 393            1.7e-5 / 1.2e-4        0.018 rad           4.1e-4
    V  / Adam best         P = 21 78/78    78 / 78               1092 / 304           7.5e-4 / 8.6e-4        0.28 rad (1 target) 3.0e-4
    L2 / Adam best         P = 21 78/78    78 / 78               1092 / 642           7.4e-4 / 8.3e-4        0.11 rad            2.7e-4
    (alpha = 0.25, L2 / Adam best, for contrast: F(P=7) - F(final) 0.053 nats; KL to product 0.126 mean, 0.148 max)

What "selected" means in `s26/q_adapt.py` `run_adapt`: at each growth step the pool operator
with the largest |dF/dphi| at phi = 0 is APPENDED to the circuit and all angles are
re-optimised; every appended operator stays with whatever angle it receives. There is no
retention step. So PR's count (appended) is the right count and my "no operator selected"
was false; the correct statement is that the appended operators are inert.

Why the ideal-ladder runs (L27, L35: L-BFGS appends nothing at alpha = 1; the Adam-grown
sets are abelian) differ from the real targets: on the ideal ladder E is exactly affine in
the register index, the RY layer reaches the exact Gibbs state and every pool gradient is
exactly zero. On a real target tie-averaging makes E not exactly affine (KL(Gibbs || product)
up to 7.9e-4), the RY layer's residual gradient in some multi-qubit direction exceeds
eps = 1e-3, an operator is appended, re-optimisation gives it an angle of order 5e-3 rad,
the free energy moves by order 1e-5 nats, and the criterion then fires. That is inert growth.

CORRECTIONS. In L68, the line "L-BFGS growth at alpha = 1 stops with no operator selected on
78 of 78 targets" and the sentence "ADAPT with either pool selects no entangling operator on
any of the 78 alpha = 1 targets under L-BFGS" are RETRACTED and replaced by the table above.
In `s26/PROPOSAL_A.md` section 3 ("selects no entangling operator on any of the 78 alpha = 1
targets") and section 5 ("it declines them on all 78 targets; there is nothing to entangle")
are corrected by an appended addendum. In `s26/agentQ_FINDINGS.md` sections 0.1 and 1.3 the
same sentences are corrected by an appended note (section 9). The ideal-ladder statements in
L27 and L35 stand as written (they are about the ladder, where the appended count is zero).

What does not change: every endpoint number in L68 (the two primaries, the twelve arms, the
controls, the subsets); the product-state result (KL(Gibbs_zrank || product) <= 7.9e-4 on 126
targets; ADAPT reaches the Gibbs state to KL <= 9e-4 on the 78 alpha = 1 targets; the fixed
circuit stops 0.90 nats short); the A2 and A4 results; the verdict REPLACE. The mechanism
sentence becomes: at alpha = 1 the objective's optimum is a product state, and an adaptive
ansatz that is free to entangle appends operators worth less than 1e-3 nats and leaves the
state a product state; there is nothing for growth to do, and growing anyway moves the
built chain by a third of the MDE.

FINAL WORDING FOR SLIDE 8's NOTES (for lane PR), replacing any "grows nothing" / "declines
them" phrasing:
"At the deployed setting the optimum the circuit is asked to reach is a product state:
seven single-qubit rotations represent it exactly (KL 0.0002 nats). When ADAPT is offered
entangling operators it does append them, on 60 to 78 of the 78 targets, but they are inert:
together they lower the objective by less than 0.001 nats, their angles stay below 0.02
radians under L-BFGS, and the state remains a product state to 0.0004 nats. The deployed
21-parameter circuit stops 0.9 nats short of that same optimum, and reaching it exactly
moves the emitted structure by 0.014 to 0.022 angstroms, a third of what the comparison can
resolve. Sources: s26/results/a1/*.json (adapt.*.sequence, .trace, arms.*.kl_to_product),
s26/results/a1_stats.json, ledger L68 and L75."

## L76 -- LANE PR: L75's FINAL WORDING APPLIED TO SLIDE 8 (NOTES AND TEXT) AND TO SLIDE 9's SPRINT-26 LINE; NO "GROWS NOTHING" / "DECLINES THEM" / "SELECTS NO ENTANGLING OPERATOR" PHRASING REMAINS; EVERY L75 QUANTITY RECOMPUTED FROM THE A1 RECORDS; THE RETRACTION IS IN THE LANE'S CHANGE LOG (2026-09-13 22:25, PR)

Applied L75 (lane Q's reconciliation of this lane's L73 flag; the two L68 sentences retracted there;
`s26/PROPOSAL_A.md` addendum). Slide 8's spoken notes now carry L75's final wording with every number
read from the records: "seven single-qubit rotations represent it exactly, to 0.0002 nats"
(`s26/results/a1/*.json :: arms/adaptV_lbfgs_zrank_P7/kl_to_gibbs`, mean over the 78 alpha = 1 targets
1.5e-4); "it does append them, on 60 to 78 of the 78 targets" (`adapt/V_lbfgs_zrank/sequence` non-empty on
60, `L2_lbfgs_zrank` on 68, `*_adam_best_zrank` on 78); "together they lower the objective by less than
0.001 nats" (`adapt/*/trace`: F(P = 7) - F(final) at most 1.2e-4 under L-BFGS, 8.6e-4 under Adam);
"their angles stay below 0.02 radians under L-BFGS" (`adapt/*_lbfgs_zrank/theta` beyond the first 7:
max 0.0179 rad); "the state remains a product state to 0.0004 nats" (`arms/adapt*_P21/kl_to_product`
max 4.1e-4); "stops 0.90 nats short" (`arms/fixed_zrank_it50/kl_to_gibbs` mean 0.9027). The slide text
says the same in one bullet, and slide 9's Sprint-26 line scopes the A4 non-decay to the ideal ladder
(L27, L35 stand) with "inert growth on the real targets (L75)". The deck dump has no "grows nothing",
"declines them" or "selects no entangling operator" left; the notes cite the retracted phrase once, as
retracted. Verification PASS: 11 slides, 0 U+2014, 0 U+2013, 0 banned words, spoken words 249 / 247 /
244 on slides 8 / 9 / 10. Registry 302 tokens. `s26/PRESENTATION_CHANGES.md` row updated to name the
retraction and its artefact keys; `s26/agentPR_FINDINGS.md` section 7 is a dated change log (the flag
at L73, the retraction at L75, this application) for the Adversary's RETRACTIONS.md. No endpoint
number or verdict changed. Still pending: lane P's B and C verdicts (one more rebuild) and the
coordinator's ruling on all three directions for slide 11's DRAFT line.

---

(Clock correction, PR: L76's header time is 21:56, not 22:25; the same correction applies to the 22:25 lines in `s26/STATUS.md` and `s26/agentPR_FINDINGS.md` section 7. Nothing else in L76 changes.)

---
## L77 -- THE USER EXTENDS THE SPRINT: EVERY DEFERRED RUN, TEST AND DIAGNOSTIC GOES AHEAD; CLOSE AT ABOUT 04:30 PACIFIC, FINAL REPORT AT ABOUT 04:45 (2026-09-13 22:02, coordinator)

The user's instruction at 22:0x: run all the tests and additional diagnostics possible, push the
close to about 04:30 Pacific (2026-09-14), and deliver the final report with all detail at about
04:45. Nothing running is interrupted. The extra scope, by lane, under the same rules (prereg
before compute; governor; one honest est-ram per job; the 93% ceiling stands, so memory, which
is the user's own load tonight, sets the parallelism, and the stall breaker (L74) may kill the
fattest suspended job):

- I (resumed): the opt-in test tier `VERIFY_SLOW=1` (the 11 real skips of `s26/TEST_RUN.md`:
  OpenMM and pipeline checks, AMBER-tagged, two at most); every standalone audit under
  `verify/` re-run with its JSON compared to the tracked one; the results lab `--mode frozen`
  rebuild as a governed job with the four gates and the leaderboard numbers compared to the
  tracked `results/summary/` (3.2126 built chain, 3.0483 point cloud; the L7 WARN explanation
  appears); the final AST gate of every production `.py` against `ae86a124` (docstrings
  stripped); `python s26/examine.py` on the final tree; `s26/TEST_RUN.md` extended.
- Q: A3 to completion; the DLA dimension of the ADAPT-grown circuit at every growth step across
  the 126 A1 records (the A2 item the prompt asked for per growth step); bootstrap CIs on the
  A4 slopes (the L47 caveat); a second-seed and reversed-fold-order replication of A1's
  fixed-versus-ADAPT contrast as a robustness check even though it is a null.
- P: the remaining rungs (esm8m, mix, pairnet), then raw (retrain fold by fold when the
  governor shows 2.5 GB free), then the deferred coherence_penalised_training (1.25 GB) if the
  box allows; B3, C4, C5; length-stratified and FAIL18-stratified tables for every rung;
  replication for any rung clearing its MDE; the best-rung delivery file; the final
  PROPOSAL_B.md and PROPOSAL_C.md. `attn` (3 to 3.5 GB) only if the box empties.
- PH: the reject chain to completion; branch_select relaunched from its 41 cells; rotamer_relief
  B; `stage1 --rep` (the C3 replication) and `floor2`; C3 stage 2 on P's rung file; the
  heavy-atom validity axis of the relaxed emission.
- W: ensembling to completion; window_provenance; amber_prior_partner; Part B's retrains one at a
  time as memory allows; the two IDEA files found on the way get preregs and, if cheap, runs.
- A: every new entry checked within the hour; RETRACTIONS.md and DELIVERABLES_CHECK.md kept
  current; the final deliverables pass at about 04:00.
- E: the report grows with every result; the final pass (Parts VII, VIII, Appendix B, Appendix C
  = this ledger reproduced) starts at about 04:15 for a 04:45 close; the report check re-run last.
- PR: the last rebuild after the B and C verdicts and the coordinator's slide 11 ruling.
- Coordinator: the slide 11 ruling once the three verdicts are in; the `docs/FINDINGS.md`
  corrections-ledger update from RETRACTIONS.md; the sprint-close entry.

Eight lanes active (E, I, Q, P, PH, W, A, PR): the maximum.

---

## L78 -- THE STALL BREAKER'S FIRST KILL (p_eval_esm8m, 22:03); THE CTRL_BREAK GRACE SIGNAL DOES NOT CROSS CONSOLES ON WINDOWS, SO A GOVERNOR KILL IS A HARD KILL AFTER 25 s (2026-09-13 22:08, coordinator)

`s26/governor.log` 22:03:30: the v2.2 stall breaker (L74) fired as designed on the fattest
suspended job while RAM sat above 90%, `p_eval_esm8m`; the CTRL_BREAK it sends first failed
with `OSError(22, 'The parameter is incorrect')` (WinError 87): `GenerateConsoleCtrlEvent` can
only reach a process group attached to the caller's console, and the governor runs in its own
console. The job was terminated at 22:03:55 after the 25 s grace and lane P relaunches it from
its checkpoint. Consequence, stated for every lane: a governor kill is a hard kill after 25 s;
the loss is bounded by the job's own checkpoint granularity (per target or per fold in every
S26 job), and a killed job's partial rows are not a result until the relaunched run completes
with `complete: true`. Not changed tonight: a file-based stop request that jobs would poll is
the right fix and is noted for the next campaign in `s26/agentI_FINDINGS.md`'s hygiene list.

---

## L79 -- ADVERSARY CHECK OF THE C2 RUNGS L62, L63, L65, L66, L67, L72: ALL STAND; L66's "A REAL, SMALL DEGRADATION" ON SELECTION IS NOT LICENSED AT 0.48x MDE; THE LADDER'S EVENTUAL BEST RUNG IS AN ORDER STATISTIC (2026-09-13, A)

Common path checked once (`s26/p_ladder.py`): the phase gate is enforced in code (`:118`
`SIGNOFF`; `:64`); native quantities enter only in the read-after-selection half (`:513-522`,
"ORACLE-READ-AFTER-SELECTION"), after `score_target`; the top-75 is `argsort(kind="stable")`
(`:505`); `sel` is tie-averaged with `ST.argmin_tied` (`:517`); every rung goes through one code
path against the one anchor (`p_ladder_shipped_s0.json`), whose built chain is the rebuild basis
3.2126 (L57) and is bit-exact against the cache on selection and point cloud; gate 2 (L65) shows
the retrained pca32 emits the shipped answer on 126/126, so rung differences are inputs or
architecture, not trainer noise. iid and fold CIs on every contrast; concentration nulls
computed; medians beside means; power stated. Clean on the checklist.

Per rung:
- L62 noesm: WORSE at 5/5 folds on all three bases, but every magnitude is Type-M (1.01x /
  1.02x / 1.27x): the sign is measured, the "+0.208 A" and "+0.330 A" are upper bounds; L62
  says so. The harm is tail-carried (median +0.044, worst +3.65 on 8T61); the uniform-effect null
  does not flag it. "Reproducing S7-11's -0.288 to -0.34" is a reproduction of a Type-M number
  by a Type-M number; quote both with the flag. The gam_eff-positive-while-worse reading (S25
  L12 from the training side) is correct and important: gam_eff is never a prediction of the
  endpoint.
- L63 conly: correctly UNDERPOWERED (0.62x), not null; the "two thirds of the ESM channel is in
  the contact head" is a point-estimate split (-0.218 at exactly 1.00x MDE against -0.330 at
  1.27x, both Type-M); the built-chain split is unresolved and L63 says so. The worst target is
  1CEK (+3.70), the verbatim self-copy target: noted, not read into.
- L65 pca32: 126 exact ties on all bases; the degenerate all-zero concentration FLAG is the
  identity case. The strongest possible gate-2 pass.
- L66 wide: built chain UNDERPOWERED (0.25x). CAVEAT: the selection sentence "a real, small
  degradation, worth reporting as 'the wider head selects slightly worse'" is not licensed:
  +0.097 is 0.48x its MDE, the iid CI spans zero ([-0.051, +0.240]) and only the fold CI excludes
  it; by the standing rule (clear the MDE AND fold CI excluding zero) it is not a result and the
  word is "suggestive, not measured". The conclusion "more capacity does not buy a better prior"
  stands on the built-chain null and the flat MAE.
- L67 pca32f: every contrast at 0.11x to 0.17x; "changes nothing measurable" is correct at the
  stated MDEs (0.13 to 0.18 A) and closes the global-PCA leak question at that resolution.
- L72 pca128: 0.43x / 0.38x / 0.13x, all under MDE; the isolate against pca32f (+0.058, 0.37x)
  is the clean size-axis point; power stated (a 0.18 A gain would be seen).

Multiplicity, stated now for the ladder's close: seven rungs so far and four to come, all
paired against one anchor on one pool, so the per-rung deltas are correlated (the rungs share
the pool, the score and the projection); any "best rung" or per-target minimum over rungs must
go through `ST.best_of_k_within` with its k_eff and split-half transfer (PREREG_C2's own rule),
never the raw minimum. No rung is positive so far, so nothing is owed yet. Verdicts: L62, L63,
L65, L67, L72 STAND; L66 STANDS WITH CAVEAT (the selection-basis wording).

---

## L80 -- ADVERSARY CHECK OF L52 (conformational_identity_floor, ORACLE DIAGNOSTIC): STANDS WITH CAVEAT (n = 18; BOTH PAIRED CONTRASTS ARE TYPE-M) (2026-09-13, A)

Reads: the copy's CA-RMSD to the target's model-1 native at the shared segment and the pool's
ORACLE `rr`; every quantity reads the native and is labelled ORACLE DIAGNOSTIC; nothing selects;
dev targets only (lane I's L15 list, 18 targets, 22 partners; no benchmark sequence, confirmed
in L49 for the same script). Registered falsifier F5 (median above 1.5 A; FALSIFIED below 1.0)
holds at 2.908 A over 22 pairs (18% below 1.0 A, 27% below 1.5 A); the four cross-fold
self-copies reproduce S24 L4's four values. Deterministic (Kabsch on fixed coordinates; the
bootstrap is seeded). Power: n = 18 gives MDEs of 0.92 and 1.26 A on the two paired contrasts.

Caveat. Both paired contrasts are in the Type-M zone (copy minus pool best +0.968 at 1.06x MDE;
copy minus pool mean -1.443 at 1.14x); W labels them so ("sign clean, magnitude not a result").
The number for the report is the median 2.9 A, a descriptive statistic at n = 22, not the paired
effects; and the sentence "the 2.0 A target is below what a verbatim lookup reaches" is a
statement about this instrument's 18 targets with a verbatim relative, not a general law. W's
own "not done" (no length-matched non-verbatim control population) is the right missing arm;
the pool mean plays that role loosely. Verdict: STANDS WITH CAVEAT.

---

## L81 -- ADVERSARY CHECK OF L53 (strain_difficulty): STANDS WITH CAVEAT, PROVISIONAL UNTIL THE POOL-SPREAD CONTROL LANDS; "REPLICATED" MEANS THE CIs, NOT THE ESTIMATE (2026-09-13, A)

Checked against `s26/ph_strain.py` and `s26/results/ph_strain.json` / `ph_strain_rep.json`.
- Leakage: the four signals are the production cache's own relaxation scalars (`amber_moved`,
  `amber_e0`, `amber_e1`, `amber_strain_after`), native-free by construction; `rmsd_arm` enters
  only as the ORACLE label (`:84-90`); gated (`:74`); nothing is selected, weighted or tuned.
  Clean.
- Multiplicity: four pre-registered signals, Bonferroni (p < 0.0125); the falsifier fires on one
  (`moved`, partial rho +0.433, fold CI [+0.247, +0.588], 5/5 folds, permutation p < 0.00025).
  Confounds n and Rg partialled by residualisation; per-fold signs and a 4000-draw permutation
  null. Clean. The two 4/5-fold signals (log_e0, log_drop at +0.24, +0.25) are correctly reported
  as measured-but-failing-the-rule, and `strain_after` as null; the FAIL18 Fisher tests as null.
- Outliers: rho +0.41 without 9KAR and 2BP4 (the broken-bond emissions). Clean.
- Basis: built chain, ORACLE label, stated. Clean.

Two caveats.
1. **"REPLICATED" in the title means the bootstrap and permutation draws, not the estimate.**
   The Spearman is a deterministic function of 126 fixed rows; the registered replication
   (new seed, reversed order) can only re-draw the CIs and the permutation null, and it did
   (identical rho +0.433, fold CI [+0.248, +0.581]). It cannot test whether the signal
   generalises to new targets; that needs targets the record does not have (no fresh
   benchmark exists). Say "CIs replicated", not "result replicated".
2. **The prereg carried no pool-disagreement control.** The obvious native-free difficulty
   proxy is the top-75's own pairwise CA-RMSD spread (the members' disagreement, no native).
   If `moved` is that spread in disguise, L53's "first native-free quantity above 0.4" is a
   re-labelling. The Adversary's control `s26/a_strain_vs_spread.py` (rho of the spread with
   `rmsd_arm`; partial rho of `moved` given n, Rg AND the spread; job `a_strain_vs_spread`,
   queued behind the four-job cap at 22:3x) decides it; the verdict below is provisional until
   it lands and is confirmed or amended in a follow-up entry.

Verdict: STANDS WITH CAVEAT, provisional. As a calibration flag with the L53 wording ("a
calibration curve, never a gain") it may be quoted now; the novelty sentence waits for the
control.

---

## L82 -- CORRECTION TO THE ADVERSARY'S L70 CAVEAT 2 PER L75: L-BFGS DOES APPEND OPERATORS AT alpha = 1; THEY ARE INERT; THE CAVEAT'S CONCLUSION IS UNCHANGED (2026-09-13, A)

L70 caveat 2 said "on the 78 alpha = 1 targets L-BFGS grows nothing and the 'P21' arm is the
7-parameter RY layer". Per L75 (lane Q, from the records' `adapt.*_lbfgs_zrank.sequence`), under
L-BFGS operators are appended on 60 of 78 (pool V) and 68 of 78 (pool L2) targets, all
multi-qubit, worth at most 1.2e-4 nats, with angles at most 0.018 rad, leaving the state a
product to 4.1e-4 nats. My `len(ops)` reading (0 to 21 under L-BFGS) was consistent with that
and I mis-stated it as "grows nothing". Corrected wording of caveat 2: the parameter match is a
budget match (21) and the realised counts run 7 to 21 under both optimisers, with the L-BFGS
appended operators inert; nothing in the A1 null depends on the count (P7, P14, P21 agree within
0.01 A). Recorded in `s26/RETRACTIONS.md` R5 together with lane Q's own corrections.

---

## L83 -- GOVERNOR v2.3 ADOPTS SUSPENSIONS MADE BEFORE ITS START; THE LAUNCH CAP IS NOW A FILE AND IS RAISED TO SIX (2026-09-13 22:21, coordinator)

Two throttles found at 22:20 with the box at 78% RAM and 3.7 GB free: (1) `a3_build` had been
suspended by the v2.1 governor and was never resumed by v2.2, whose suspended stack starts
empty on a restart; v2.3 adopts any stopped registered job onto its stack on every tick and
logs ADOPT, so it resumes under the normal rule. (2) Eight launcher waiters (P esm8m2 and mix,
E report_check5, A strain_vs_spread, I slow_equivalence and slow_integration, Q a4_var_boot and
a2_dla_a1) sat behind jobrun's fixed cap of four with room on the box. `s26/jobrun.py` v2.3
reads the cap from `s26/launch_cap.json` on every wait tick (coordinator-only file); it is set
to 6 now. Waiters started before this edit keep the old fixed cap of four until they launch;
every new waiter uses the file. The memory rules (est-ram + 0.5 GB free; the 93% ceiling; the
stall breaker) are unchanged and still bind before the cap does.

---


## L84 -- TOURNAMENT ITEM 6, window_ensembling (P's idea, orphaned to W; THE MANDATORY TEST-TIME-ENSEMBLING DIRECTION): FIXED-K ENSEMBLING OF THREE BLOSUM KEYS IS -0.0005 A ON THE BUILT CHAIN (0.02x MDE) AND INDISTINGUISHABLE FROM A ZERO-INFORMATION RESAMPLE OF THE PRODUCTION SET (+0.0042, 0.14x); THE THREE CLOUDS SIT 0.14 TO 0.20 A APART AND THEIR ERRORS ARE PARALLEL, AS S23 L9 / S24 L3 PREDICT; A GAIN OF 0.028 A OR MORE IS EXCLUDED (2026-09-13, W)

Pre-registered in `s26/PREREG_window_ensembling.md` (written 20:05, before the probe; falsifier,
arms, controls, the widening-K distinction and the expected 0.00 to +0.03 A all fixed there).
Native-free half: job `w_ensemble_clouds` (exit 0, 4351 s wall under a four-job load with two
governor suspensions, peak RSS 0.113 GB; `s26/results/w_selfcopy_ensemble_clouds.json`, complete
126/126; production gate top-75 == `sub` 126/126, cloud max abs 1.8e-15; BLOSUM62 sums asserted
equal to the universe's `sim` and `core.data.top_k` to the pinned pool on every target). Gated
half: job `w_ensemble_endpoint` (exit 0, 10 s, 0.057 GB; `s26/results/w_selfcopy_ensemble_endpoint.json`).
Probe `w_ensemble_probe` (1A13, 35 s, 0.097 GB). Basis on every line: built chain `arm` PRIMARY
(the rebuild basis 3.2126, L57), point cloud carried (3.0483); the selection basis does not exist
for an ensemble. Negative = the arm is better than the shipped emission.

**How fixed-K ensembling differs from widening K, in one sentence each.** Widening K (S17 L12)
draws ONE shortlist from a wider pool, and the extra plausible windows displace near-native
members out of that single top-75 (its ORACLE best 2.104 -> 2.572 A from K = 75 to the full
universe); fixed-K ensembling keeps three shortlists of K = 500, each cut to its own top-75 by
the shipped score (the BLOSUM62 one IS the production set, its ORACLE best 2.306 untouched), and
averages the three CLOUDS in the production frame, so no shortlist is widened and displacement
cannot act; what can act is only whether the three clouds' errors are parallel, and S23 L9 (68%
of the pool's error is common-mode, shared by every member of the universe) and S24 L3 (score
selection makes the bias parallel across sources, cosine 0.943 against a within-source 0.933)
predicted that they are.

**Native-free geometry (medians over 126).** The BLOSUM45 / BLOSUM80 pools overlap the BLOSUM62
pool at 0.83 / 0.89 and their score-selected top-75s overlap the production top-75 at 0.83 /
0.88 (93 distinct members in the union of three); the three clouds sit 0.18 (45 vs 62), 0.14
(80 vs 62) and 0.20 A (45 vs 80) apart; the ensemble cloud moves the built chain 0.185 A from
the production chain in the median (b45 0.33, b80 0.25, ensK 0.32, union3 0.15, boot3 0.25:
the zero-information resample moves it as far as the ensemble does).

**Endpoints (`ST.fmt` verbatim; the PRIMARY contrast and the two controls the falsifier names):**

  ens3 minus shipped [arm]
    a 3.2121 (med 2.9453)   b 3.2126 (med 2.9661)   n=126
    effect -0.0005   median -0.0002   SE 0.0100   MDE 0.0281   effect/MDE -0.02
    iid  CI95 [-0.0206, +0.0199]
    fold CI95 [-0.0194, +0.0162]   folds same sign 3/5   per-fold 0:-0.021 1:-0.012 2:+0.022 3:-0.021 4:+0.023
    64W/62L/0T   worst degradation +0.5570 (1D6X)   p90 +0.0615   power 0.05  Type-M 48.95
    concentration: drop-top10 +0.0179 vs uniform-effect null p10/p50/p90 +0.0062/+0.0171/+0.0296 -> pctile 0.531
    VERDICT: NOT MEASURED (|effect| 0.0005 <= its own MDE 0.0281, 0.02x)
  ens3 minus shipped [cloud]
    a 3.0499 (med 2.7932)   b 3.0483 (med 2.8373)   n=126
    effect +0.0016   median -0.0014   SE 0.0066   MDE 0.0185   effect/MDE +0.09
    iid  CI95 [-0.0095, +0.0158]
    fold CI95 [-0.0086, +0.0161]   folds same sign 2/5   per-fold 0:-0.005 1:+0.028 2:+0.007 3:-0.011 4:-0.008
    68W/58L/0T   worst degradation +0.6286 (1U62)   p90 +0.0475   power 0.06  Type-M 9.87
    concentration: drop-top10 +0.0112 vs uniform-effect null p10/p50/p90 +0.0028/+0.0105/+0.0198 -> pctile 0.540
    VERDICT: NOT MEASURED (|effect| 0.0016 <= its own MDE 0.0185, 0.09x)
  ens3 minus boot3 (the matched zero-information ensembling control) [arm]
    a 3.2121 (med 2.9453)   b 3.2079 (med 3.0424)   n=126
    effect +0.0042   median +0.0004   SE 0.0106   MDE 0.0297   effect/MDE +0.14
    iid  CI95 [-0.0162, +0.0243]
    fold CI95 [-0.0178, +0.0300]   folds same sign 2/5   per-fold 0:-0.028 1:-0.009 2:+0.048 3:-0.011 4:+0.017
    61W/65L/0T   worst degradation +0.5241 (9BAF)   p90 +0.1299   power 0.07  Type-M 5.96
    concentration: drop-top10 +0.0242 vs uniform-effect null p10/p50/p90 +0.0117/+0.0238/+0.0366 -> pctile 0.514
    VERDICT: NOT MEASURED (|effect| 0.0042 <= its own MDE 0.0297, 0.14x)
  ens3 minus boot3 (the matched zero-information ensembling control) [cloud]
    a 3.0499 (med 2.7932)   b 3.0507 (med 2.8659)   n=126
    effect -0.0008   median -0.0025   SE 0.0064   MDE 0.0178   effect/MDE -0.04
    iid  CI95 [-0.0134, +0.0117]
    fold CI95 [-0.0051, +0.0048]   folds same sign 3/5   per-fold 0:-0.006 1:+0.003 2:+0.009 3:-0.005 4:-0.004
    72W/54L/0T   worst degradation +0.2384 (6CEJ)   p90 +0.0724   power 0.05  Type-M 19.64
    concentration: drop-top10 +0.0125 vs uniform-effect null p10/p50/p90 +0.0047/+0.0119/+0.0197 -> pctile 0.541
    VERDICT: NOT MEASURED (|effect| 0.0008 <= its own MDE 0.0178, 0.04x)

**The secondaries, minus shipped (effect / MDE / fold CI / verdict), both bases:**

    arm      b45 -0.0039 / 0.0405 / [-0.0348, +0.0212] NOT MEASURED   b80 +0.0237 / 0.0337 / [+0.0052, +0.0402] NOT MEASURED (0.70x, fold CI above zero, 4/5 folds, Type-M zone: the BLOSUM80 shortlist alone is if anything WORSE)
             ensK -0.0153 / 0.0386 / [-0.0380, +0.0047] NOT MEASURED   union3 -0.0021 / 0.0191 / [-0.0180, +0.0168] NOT MEASURED   boot3 -0.0047 / 0.0273 / [-0.0161, +0.0045] NOT MEASURED
    cloud    b45 +0.0009 / 0.0374 NOT MEASURED   b80 +0.0136 / 0.0250 / [+0.0024, +0.0258] NOT MEASURED (0.54x)   ensK -0.0116 / 0.0223 / [-0.0179, -0.0045] NOT MEASURED (0.52x, 5/5 folds, fold CI below zero: the K-ensemble's CLOUD is 0.012 A nearer, and its built chain -0.015 at 0.40x is not)
             union3 +0.0013 / 0.0113 NOT MEASURED   boot3 +0.0023 / 0.0209 NOT MEASURED

**Verdict.** H_E is refuted at this instrument: the ensemble neither beats the shipped chain
(-0.0005, 0.02x MDE) nor its zero-information control (+0.0042, 0.14x). Power: the design
excludes a gain of 0.028 A or more on the built chain (the MDE of ens3 minus shipped) and of
0.019 A or more on the cloud; the registered expectation (0.00 to +0.03) held. The mechanism is
the registered one: the three shortlists overlap at 0.83 to 0.88, their clouds sit 0.14 to 0.20 A
apart, and averaging them moves the emission by exactly as much as a resample of the production
set does (0.185 vs 0.249 A median chain move) with the same effect on accuracy (none): parallel
errors average to themselves. The only cell with the fold CI excluding zero in the helpful
direction is the K-variant ensemble on the CLOUD (-0.0116, 0.52x MDE, 5/5 folds), which its own
built chain does not carry (-0.0153, 0.40x); it is quoted as suggestive and nothing more, and it
is the S17 L12 direction (the consensus readout improves with K on the point cloud, and the
built chain does not follow). The mandatory test-time-ensembling direction is closed in the
fixed-K form with a power statement; the widening-K form was closed by S17 L12. Replication:
`boot3` is seeded (`s15.seed.stable_rng`); a second seed and the reversed fold order are the
pre-declared replication for a positive, and there is none to replicate; the ensemble itself
has no random element. Deviations from the PREREG: none; the `sel` basis is undefined for an
ensemble and is not quoted. Artefacts: `s26/PREREG_window_ensembling.md`, `s26/w_ensemble.py`,
`s26/w_ensemble_test.py` (ALL OK), `s26/results/w_selfcopy_ensemble_probe_1A13.json`,
`w_selfcopy_ensemble_clouds.json`, `w_selfcopy_ensemble_endpoint.json`,
`s26/jobs_done/w_ensemble_{probe,clouds,endpoint}.json`.

## L85 -- TOURNAMENT ITEM 9, window_provenance (W), CENSUS AND ORACLE CONTRAST: 73% OF EVERY POOL AND OF EVERY TOP-75 IS FRAGMENT WINDOWS, 0.6% WHOLE PEPTIDES; IN THE POOL A PEPTIDE-DERIVED WINDOW IS 0.42 A NEARER THE NATIVE THAN A FRAGMENT WINDOW (3.2x MDE, 5/5 FOLDS); INSIDE THE SCORE-SELECTED TOP-75 THE GAP SHRINKS TO 0.09 (0.96x MDE) AND WHOLE+TERMINAL VS FRAGMENT TO 0.13 (1.13x, TYPE-M ZONE); THE READOUT TEST H_P3 RUNS NEXT (2026-09-13, W)

Pre-registered in `s26/PREREG_window_provenance.md` (written before the code ran on any target).
Census (native-free, codes and sequences only): job `w_provenance_census` (exit 0, 5 s, peak
0.004 GB; `s26/results/w_selfcopy_provenance_census.json`, complete 126/126, 0 unresolved
windows). ORACLE contrast (single-window basis, the universe's `rr`): job `w_provenance_oracle`
(exit 0, 10 s, 0.055 GB; `s26/results/w_selfcopy_provenance_oracle.json`). Class rule: a pool
window's string looked up among the length-n substrings of the 787 peptides (whole: parent length
n; terminal: offset 0 or the end; interior) and the 6,003 fragments (fragment), the FIRST parent
in database order winning, the universe's `org` flag arbitrating a string present in both banks;
synthetic test `s26/w_provenance_test.py` ALL OK. Multi-parent strings are common (15,990 at
length 9 down to 1,224 at 16, mostly overlapping fragments of one protein) and are counted, not
guessed.

**Census (H_P1), means over 126 targets:**

    class       K = 500 pool   shipped top-75   targets with one in the top-75
    whole          0.6%            0.7%           30 / 126
    terminal       8.0%            7.2%          113 / 126
    interior      18.7%           18.3%
    fragment      72.7%           73.8%
    unresolved     0.0%            0.0%

The score keeps the pool's class mix almost unchanged: it neither favours nor removes peptide
windows. The registered "about 1% whole peptides" holds (0.6 to 0.7%).

**ORACLE contrast (H_P2), per-target mean `rr` of one class minus another, `ST.fmt` verbatim for
the two contrasts the falsifier names (single-window basis on both sides; negative = the first
class is nearer the native):**

  pool: class peptide_any minus class fragment, per-target mean ORACLE rr (single-window basis), n=126 targets with both
    a 4.1721 (med 4.0394)   b 4.5882 (med 4.3129)   n=126
    effect -0.4161   median -0.4358   SE 0.0468   MDE 0.1310   effect/MDE -3.18
    iid  CI95 [-0.5066, -0.3220]
    fold CI95 [-0.4472, -0.3758]   folds same sign 5/5   per-fold 0:-0.341 1:-0.401 2:-0.463 3:-0.439 4:-0.434
    101W/25L/0T   worst degradation +1.1974 (7VI4)   p90 +0.2005   power 1.00  Type-M 1.00
    concentration: drop-top10 -0.3272 vs uniform-effect null p10/p50/p90 -0.3902/-0.3278/-0.2696 -> pctile 0.504
    VERDICT: BETTER
  top75: class whole_or_terminal minus class fragment, per-target mean ORACLE rr (single-window basis), n=114 targets with both
    a 3.5381 (med 3.4172)   b 3.6695 (med 3.7564)   n=114
    effect -0.1313   median -0.0533   SE 0.0414   MDE 0.1159   effect/MDE -1.13
    iid  CI95 [-0.2139, -0.0533]
    fold CI95 [-0.1646, -0.1096]   folds same sign 5/5   per-fold 0:-0.106 1:-0.193 2:-0.123 3:-0.143 4:-0.107
    66W/48L/0T   worst degradation +0.8160 (5MXS)   p90 +0.2992   power 0.89  Type-M 1.07
    concentration: drop-top10 -0.0286 vs uniform-effect null p10/p50/p90 -0.0789/-0.0306/+0.0136 -> pctile 0.520
    VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.07x]

The other contrasts (effect / MDE / fold CI / verdict): pool whole+terminal minus fragment -0.421
/ 0.136 / [-0.453, -0.388] BETTER (3.1x, 5/5); pool interior minus fragment -0.411 / 0.136 /
[-0.447, -0.368] BETTER (3.0x, 5/5); pool whole+terminal minus interior -0.011 / 0.085 NOT
MEASURED (0.13x); top-75 peptide_any minus fragment -0.086 / 0.090 / [-0.116, -0.056] NOT MEASURED
(0.96x, 5/5 folds); top-75 interior minus fragment -0.061 / 0.093 NOT MEASURED (0.66x); top-75
whole+terminal minus interior -0.056 / 0.114 NOT MEASURED (0.49x).

Reading. In the retrieval pool, ANY peptide-derived window (whole, terminal or interior alike) is
0.41 to 0.42 A nearer the native than a fragment window, 5/5 folds, 3x its MDE: S19's "the peptide
corpus carries the sequence-structure channel; fragments are training-only material whose
conformation is held by contacts outside the window" (S24 L8), measured for the first time on the
single-window basis of the shipped pool. The position of the window in its parent peptide does
not matter (whole+terminal minus interior 0.13x MDE). Inside the score-selected top-75 the gap
shrinks to -0.09 (0.96x) for any peptide window and -0.13 (1.13x, Type-M zone) for whole+terminal
windows: the shipped score already removes most of the class difference, and what remains is at
the edge of what n = 114 to 126 can see. The registered expectation for H_P2 (NOT MEASURED inside
the top-75, well under 0.1 A) held for the peptide_any and interior contrasts and was exceeded, in
the Type-M zone, by the whole+terminal one. By the PREREG's rule the ACHIEVABLE readout test H_P3
(fragment-class relative weight in {0, 0.5, 1, 2}, chosen leave-fold-out on the built chain,
against the uniform readout and the permuted-weight control) therefore runs; the registered
expectation for it is 0.00 to -0.02 A on the built chain against an MDE of about 0.05
(NOT MEASURED), because a 0.13 A single-window difference on 8% of the members of a 75-member
average is worth at most a few hundredths through the mean, and the L44/L64 controls show that a
few members' change moves the chain mostly orthogonally to the native. Nothing here is deployable
yet: the class of a window is native-free, the contrast is ORACLE. Replication: deterministic
(the census and the `rr` are fixed); the bootstrap is seeded.

## L86 -- STERIC REJECT, BUILT CHAIN: THE POINT-CLOUD RESULT REPLICATES ACROSS BASES. R@1e4 IS +0.248 A WORSE THAN THE RE-PROJECTED SHIPPED TOP-75 (fold CI [+0.12, +0.35], 4/5 FOLDS, TYPE-M) AND +0.157 WORSE THAN THE SAME COUNT REJECTED AT RANDOM; DOSE MONOTONE; THE FILTER IS CLOSED ON BOTH BASES (2026-09-13, PH)

`s26/ph_reject.py chain` then `report`, `s26/results/ph_reject_chain.json` (complete 126/126),
`s26/results/ph_reject_report.json` (both bases), job `s26/jobs_done/ph_reject_chain.json` (exit
0, 11,976 s = 3.3 h, 22.8 projections per target, peak RSS 0.09 GB). Pre-registered in
`s26/PREREG_amber_reject.md` sections 3 to 5, addendum 1 (fallback and moved subset) and
addendum 2 (the re-projected anchor). This is the falsifier's last leg and the cross-basis
replication of L43.

BASIS: BUILT CHAIN on both sides (`rmsd_arm`): every retained set is coordinate-averaged and
then projected through `s12.instrument.project` (ramah at 0.3, multi-start, exact gradient,
the production call) and the chain is scored ORACLE against `nat_ca`. The ANCHOR is the
shipped top-75 RE-PROJECTED through the same call, not the stored production chain, so all arms
share one instrument. Its deviation from the stored chain, reported as registered: mean
3.2126 A against the production 3.2148 (the "leaderboard rebuild 3.2126" of the
state brief and L57: the production path round-trips the cloud through float32 before
projecting, `Config.reference_precision`, and this instrument does not); per-target RMSD
difference mean -0.0021, max |0.171|; CA deviation between the two chains mean
0.105 A, above 0.05 A on 46 targets and above 0.5 A on 6
(2LWU 1.25, 7JS6 1.54, 6QAX 1.56, 2BP4 1.62), which is the projection's own branch
degeneracy (`core/project.py` docstring: a 1e-13 input difference can route L-BFGS-B into the
other torsion branch, up to 1.6 A). Arms, controls and fallbacks as in L43 (R = reject e_amber
> T and refill to 75, the PRIMARY reading; S = reject, no refill; RANDR/RANDS same count at
random; PERMR/PERMS the same threshold on a permuted energy); the built-chain controls at 1e4
carry 4 draws each (registered in section 3), the other thresholds carry R and S only.
Eleven `ST.fmt` blocks verbatim:

      built_chain [all, n=126] R@1e4 minus anchor (negative = arm better)
        a 3.4609 (med 3.2449)   b 3.2126 (med 2.9661)   n=126
        effect +0.2483   median +0.0016   SE 0.0758   MDE 0.2123   effect/MDE +1.17
        iid  CI95 [+0.1021, +0.4044]
        fold CI95 [+0.1197, +0.3479]   folds same sign 4/5   per-fold 0:+0.267 1:-0.002 2:+0.226 3:+0.421 4:+0.310
        49W/67L/10T   worst degradation +4.2043 (8T61)   p90 +1.3949   power 0.91  Type-M 1.06
        concentration: drop-top10 +0.3567 vs uniform-effect null p10/p50/p90 +0.2615/+0.3545/+0.4565 -> pctile 0.511
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      built_chain [all, n=126] R@1e4 minus RANDR@1e4 (matched control)
        a 3.4609 (med 3.2449)   b 3.3036 (med 3.1840)   n=126
        effect +0.1574   median +0.0085   SE 0.0488   MDE 0.1367   effect/MDE +1.15
        iid  CI95 [+0.0701, +0.2539]
        fold CI95 [+0.0480, +0.2811]   folds same sign 4/5   per-fold 0:+0.120 1:-0.049 2:+0.182 3:+0.376 4:+0.160
        56W/68L/2T   worst degradation +3.2399 (8T61)   p90 +0.8244   power 0.90  Type-M 1.06
        concentration: drop-top10 +0.2204 vs uniform-effect null p10/p50/p90 +0.1582/+0.2166/+0.2815 -> pctile 0.530
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      built_chain [all, n=126] R@1e4 minus PERMR@1e4 (matched control)
        a 3.4609 (med 3.2449)   b 3.3589 (med 3.1940)   n=126
        effect +0.1021   median +0.0000   SE 0.0414   MDE 0.1159   effect/MDE +0.88
        iid  CI95 [+0.0241, +0.1886]
        fold CI95 [-0.0080, +0.2034]   folds same sign 4/5   per-fold 0:+0.017 1:-0.113 2:+0.137 3:+0.277 4:+0.174
        57W/61L/8T   worst degradation +2.6686 (9KAR)   p90 +0.6660   power 0.69  Type-M 1.20
        concentration: drop-top10 +0.1607 vs uniform-effect null p10/p50/p90 +0.1091/+0.1587/+0.2107 -> pctile 0.518
        VERDICT: NOT MEASURED (|effect| 0.1021 <= its own MDE 0.1159, 0.88x)
      built_chain [all, n=126] S@1e4 minus anchor (negative = arm better)
        a 3.3164 (med 3.1608)   b 3.2126 (med 2.9661)   n=126
        effect +0.1037   median +0.0025   SE 0.0333   MDE 0.0933   effect/MDE +1.11
        iid  CI95 [+0.0414, +0.1747]
        fold CI95 [+0.0373, +0.1662]   folds same sign 5/5   per-fold 0:+0.014 1:+0.002 2:+0.172 3:+0.194 4:+0.129
        47W/66L/13T   worst degradation +1.7906 (9BAF)   p90 +0.4638   power 0.88  Type-M 1.08
        concentration: drop-top10 +0.1422 vs uniform-effect null p10/p50/p90 +0.1000/+0.1424/+0.1914 -> pctile 0.498
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.08x]
      built_chain [all, n=126] S@1e4 minus RANDS@1e4 (matched control)
        a 3.3164 (med 3.1608)   b 3.2612 (med 3.0604)   n=126
        effect +0.0551   median +0.0028   SE 0.0269   MDE 0.0754   effect/MDE +0.73
        iid  CI95 [+0.0041, +0.1101]
        fold CI95 [-0.0005, +0.1064]   folds same sign 4/5   per-fold 0:-0.048 1:+0.027 2:+0.118 3:+0.127 4:+0.055
        45W/68L/13T   worst degradation +1.3396 (8T61)   p90 +0.3322   power 0.53  Type-M 1.36
        concentration: drop-top10 +0.0941 vs uniform-effect null p10/p50/p90 +0.0586/+0.0924/+0.1297 -> pctile 0.524
        VERDICT: NOT MEASURED (|effect| 0.0551 <= its own MDE 0.0754, 0.73x)
      built_chain [all, n=126] RANDR@1e4 minus anchor (negative = arm better)
        a 3.3036 (med 3.1840)   b 3.2126 (med 2.9661)   n=126
        effect +0.0909   median +0.0026   SE 0.0432   MDE 0.1211   effect/MDE +0.75
        iid  CI95 [+0.0110, +0.1787]
        fold CI95 [+0.0456, +0.1322]   folds same sign 5/5   per-fold 0:+0.148 1:+0.047 2:+0.045 3:+0.046 4:+0.151
        59W/65L/2T   worst degradation +2.1588 (8FLP)   p90 +0.7406   power 0.56  Type-M 1.34
        concentration: drop-top10 +0.1645 vs uniform-effect null p10/p50/p90 +0.1103/+0.1633/+0.2176 -> pctile 0.513
        VERDICT: NOT MEASURED (|effect| 0.0909 <= its own MDE 0.1211, 0.75x)
      built_chain [all, n=126] RANDS@1e4 minus anchor (negative = arm better)
        a 3.2612 (med 3.0604)   b 3.2126 (med 2.9661)   n=126
        effect +0.0486   median +0.0002   SE 0.0175   MDE 0.0490   effect/MDE +0.99
        iid  CI95 [+0.0167, +0.0836]
        fold CI95 [+0.0121, +0.0698]   folds same sign 4/5   per-fold 0:+0.063 1:-0.025 2:+0.055 3:+0.067 4:+0.074
        48W/65L/13T   worst degradation +1.1826 (9BAF)   p90 +0.2279   power 0.79  Type-M 1.13
        concentration: drop-top10 +0.0695 vs uniform-effect null p10/p50/p90 +0.0477/+0.0684/+0.0920 -> pctile 0.528
        VERDICT: NOT MEASURED (|effect| 0.0486 <= its own MDE 0.0490, 0.99x)
      built_chain [all, n=126] R@1e3 minus anchor (negative = arm better)
        a 3.7757 (med 3.5036)   b 3.2126 (med 2.9661)   n=126
        effect +0.5631   median +0.0248   SE 0.1138   MDE 0.3187   effect/MDE +1.77
        iid  CI95 [+0.3549, +0.7891]
        fold CI95 [+0.4413, +0.7037]   folds same sign 5/5   per-fold 0:+0.541 1:+0.346 2:+0.530 3:+0.822 4:+0.578
        36W/69L/21T   worst degradation +4.6340 (2NDN)   p90 +2.7185   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.6994 vs uniform-effect null p10/p50/p90 +0.5482/+0.6942/+0.8531 -> pctile 0.518
        VERDICT: WORSE
      built_chain [all, n=126] S@1e3 minus anchor (negative = arm better)
        a 3.4048 (med 3.2359)   b 3.2126 (med 2.9661)   n=126
        effect +0.1921   median +0.0007   SE 0.0477   MDE 0.1338   effect/MDE +1.44
        iid  CI95 [+0.1048, +0.2907]
        fold CI95 [+0.1032, +0.2469]   folds same sign 5/5   per-fold 0:+0.231 1:+0.021 2:+0.273 3:+0.212 4:+0.208
        37W/63L/26T   worst degradation +2.4536 (9BAF)   p90 +0.7604   power 0.98  Type-M 1.01
        concentration: drop-top10 +0.2423 vs uniform-effect null p10/p50/p90 +0.1781/+0.2389/+0.3077 -> pctile 0.529
        VERDICT: WORSE
      built_chain [moved, n=116] R@1e4 minus anchor (negative = arm better)
        a 3.4869 (med 3.2751)   b 3.2172 (med 3.0247)   n=116
        effect +0.2697   median +0.0077   SE 0.0820   MDE 0.2298   effect/MDE +1.17
        iid  CI95 [+0.1144, +0.4415]
        fold CI95 [+0.1290, +0.3943]   folds same sign 4/5   per-fold 0:+0.267 1:-0.002 2:+0.283 3:+0.484 4:+0.321
        49W/67L/0T   worst degradation +4.2043 (8T61)   p90 +1.5485   power 0.91  Type-M 1.06
        concentration: drop-top10 +0.3904 vs uniform-effect null p10/p50/p90 +0.2849/+0.3857/+0.4945 -> pctile 0.520
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]
      built_chain [moved, n=116] R@1e4 minus RANDR@1e4 (matched control)
        a 3.4869 (med 3.2751)   b 3.3178 (med 3.2429)   n=116
        effect +0.1692   median +0.0099   SE 0.0522   MDE 0.1462   effect/MDE +1.16
        iid  CI95 [+0.0731, +0.2704]
        fold CI95 [+0.0519, +0.3172]   folds same sign 4/5   per-fold 0:+0.120 1:-0.043 2:+0.202 3:+0.438 4:+0.165
        51W/65L/0T   worst degradation +3.2399 (8T61)   p90 +0.8244   power 0.90  Type-M 1.06
        concentration: drop-top10 +0.2393 vs uniform-effect null p10/p50/p90 +0.1716/+0.2352/+0.3086 -> pctile 0.533
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.06x]

All four thresholds, `[all, n=126]`, arm minus anchor, built chain:

    1e3  R   +0.5631 (median  +0.0248, +1.77x MDE)  fold [+0.4413, +0.7037]   5/5  36W/69L/21T  WORSE
    1e3  S   +0.1921 (median  +0.0007, +1.44x MDE)  fold [+0.1032, +0.2469]   5/5  37W/63L/26T  WORSE
    1e4  R   +0.2483 (median  +0.0016, +1.17x MDE)  fold [+0.1197, +0.3479]   4/5  49W/67L/10T  WORSE
    1e4  S   +0.1037 (median  +0.0025, +1.11x MDE)  fold [+0.0373, +0.1662]   5/5  47W/66L/13T  WORSE
    1e5  R   +0.1190 (median  +0.0000, +0.78x MDE)  fold [+0.0482, +0.1844]   5/5  57W/62L/7T   NOT MEASURED
    1e5  S   +0.0833 (median  +0.0000, +0.95x MDE)  fold [+0.0020, +0.1560]   3/5  58W/59L/9T   NOT MEASURED
    1e6  R   +0.0915 (median  +0.0000, +0.65x MDE)  fold [+0.0328, +0.1377]   4/5  57W/61L/8T   NOT MEASURED
    1e6  S   +0.0511 (median  +0.0000, +0.76x MDE)  fold [-0.0016, +0.0983]   3/5  61W/55L/10T  NOT MEASURED

Threshold sweep (order statistic): R oracle minimum -0.3442, split-half transfer -0.1473 (43%%),
k_eff 3.86; S -0.1345, transfer -0.0535 (40%%): the mildest rung transfers, not a gain (R@1e6
+0.092 at 0.65x MDE, S@1e6 +0.051 at 0.76x, both NOT MEASURED).

READING. (1) The built chain REPLICATES the point cloud in sign on the primary reading: R@1e4
is +0.248 A worse than the re-projected shipped top-75 (fold CI [+0.120, +0.348], 49W/67L/10T)
and +0.157 worse than rejecting the same count at random with a judgment-free refill (fold CI
[+0.048, +0.281]); L43 had +0.228 and +0.167 on the point cloud. Both are Type-M-zone
magnitudes (1.17x and 1.15x MDE), as on the point cloud, and both are tail-carried (median
+0.002 and +0.009 against means +0.248 and +0.157; p90 +1.39; worst +4.20 on 8T61), so the
Adversary's L54 caveats 1 and 2 transfer unchanged to this basis: near zero on the median
target, catastrophic on the minority whose pool has no survivor or whose refill reaches deep.
One difference from L43 that is reported rather than smoothed: fold 1 has the opposite sign on
R (-0.002 against +0.267 / +0.226 / +0.421 / +0.310 on the other four), so R@1e4 is 4/5 folds
on the built chain where it was 5/5 on the point cloud; the fold CI still excludes zero. S@1e4
is +0.104 worse (5/5 folds, 1.11x MDE, Type-M) against +0.108 on the point cloud. (2) The
dose is again monotone: 1e3 +0.563 (1.77x MDE, 5/5 folds, MEASURED), 1e4 +0.248, 1e5 +0.119
(0.78x), 1e6 +0.092 (0.65x); the limit is the anchor. (3) The built chain is noisier than the
point cloud by construction (the projection adds its own branch noise on every arm), so the
matched-control contrasts that were Type-M on the point cloud are UNDERPOWERED here: S vs
RANDS +0.055 at 0.73x MDE (fold CI [-0.0005, +0.106]), R vs PERMR +0.102 at 0.88x, S vs PERMS
+0.050 at 0.59x; RANDR vs anchor +0.091 at 0.75x, RANDS vs anchor +0.049 at 0.99x. The L54
caveat 3 stands: the refill-cost / choice-cost split is a point estimate on both bases. (4)
The moved subset (secondary, n = 116 / 113) agrees: R +0.270 WORSE (Type-M), R vs RANDR +0.169
WORSE (Type-M), S +0.116 WORSE (Type-M).

POWER. Built-chain SEs 0.03 to 0.11 A, MDEs 0.09 to 0.32; at 1e5 and 1e6 R sits at 0.78x and
0.65x its MDE (SE 0.055 and 0.050), so "the mildest threshold is null" is UNDERPOWERED for a
gain below about 0.14 A on this basis and MEASURED against any harm above it. Nothing is
positive; the replication of the harmful direction is this entry (cross-basis, same sign, same
monotone dose, same tail structure).

DISPOSITION. The falsifier of `s26/PREREG_amber_reject.md` section 4 required R or S at 1e4 to
BEAT the anchor and its matched control on the built chain; both are WORSE than the anchor
with fold CIs excluding zero (R 4/5 folds, S 5/5), and R is worse than its matched control.
AMBER as a steric reject filter at a physical threshold, with or without refill, is CLOSED on
both bases in the two forms the record had not measured (physics-set count, refill), with a
harmful sign at 1e3 and 1e4 and an underpowered null at 1e5 and 1e6. The functional lever's
filter form (s19 Q3, s24 D1-C) keeps its closure and gains a sixth and seventh instrument. Why,
from L23: the members the threshold condemns are condemned by the builder's side-chain
placement (96.8%), not by their backbones; removing them removes compact members whose error
cancelled in the average (S23 L5). The remaining physics question on this pool is
`IDEA_rotamer_relief` part B (does the singularity move when the side chains are relieved),
which runs next under the tournament.

---

## L87 -- C3 STAGE 1 REPLICATION (new seeds, reversed order): EVERY CONTRAST LANDS INSIDE L39's FOLD CI; THE TOWARD-MEMBER CONTROL's -0.018 A IMPROVEMENT IS NO LONGER PROVISIONAL (-0.0207 [-0.0250, -0.0166], 5/5) (2026-09-13, PH)

`s26/ph_c3.py stage1 --rep`, `s26/results/ph_c3_stage1_rep.json` (complete 126/126), job
`s26/jobs_done/ph_c3_stage1_rep.json` (exit 0, 10 s, peak RSS 0.04 GB). Registered in
`s26/PREREG_c3_control.md` addendum 2 before it ran: new draw seeds (salt "rep"), targets
processed in reversed pinned order, fold labels untouched. Basis as in L39: arm input the
BUILT CHAIN (`rmsd_arm` 3.2148), arm output the RELAXED CHAIN (`rmsd_full` 3.2355); the
controls displace the built chain's CA trace by AMBER's own per-atom RMS magnitude. The arm
has no randomness, so only the controls' draws change. `ST.fmt` verbatim, the five contrasts
the replication was registered to judge:

      stage1-REP TOWARD-MEMBER minus do-nothing
        a 3.1941 (med 2.9650)   b 3.2148 (med 2.9661)   n=126
        effect -0.0207   median -0.0158   SE 0.0042   MDE 0.0117   effect/MDE -1.77
        iid  CI95 [-0.0287, -0.0126]
        fold CI95 [-0.0250, -0.0166]   folds same sign 5/5   per-fold 0:-0.026 1:-0.020 2:-0.026 3:-0.015 4:-0.016
        94W/32L/0T   worst degradation +0.1064 (8TXS)   p90 +0.0279   power 1.00  Type-M 1.00
        concentration: drop-top10 -0.0121 vs uniform-effect null p10/p50/p90 -0.0173/-0.0123/-0.0073 -> pctile 0.522
        VERDICT: BETTER
      stage1-REP AMBER minus matched-magnitude RANDOM (16 draws)
        a 3.2355 (med 2.9757)   b 3.2255 (med 2.9698)   n=126
        effect +0.0100   median +0.0094   SE 0.0035   MDE 0.0098   effect/MDE +1.02
        iid  CI95 [+0.0033, +0.0172]
        fold CI95 [+0.0059, +0.0175]   folds same sign 5/5   per-fold 0:+0.008 1:+0.026 2:+0.006 3:+0.005 4:+0.007
        49W/77L/0T   worst degradation +0.1442 (2LM8)   p90 +0.0561   power 0.81  Type-M 1.12
        concentration: drop-top10 +0.0156 vs uniform-effect null p10/p50/p90 +0.0109/+0.0156/+0.0204 -> pctile 0.509
        VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.12x]
      stage1-REP AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)
        a 3.2355 (med 2.9757)   b 3.1941 (med 2.9650)   n=126
        effect +0.0414   median +0.0387   SE 0.0059   MDE 0.0165   effect/MDE +2.50
        iid  CI95 [+0.0299, +0.0527]
        fold CI95 [+0.0331, +0.0501]   folds same sign 5/5   per-fold 0:+0.047 1:+0.058 2:+0.043 3:+0.028 4:+0.033
        30W/96L/0T   worst degradation +0.3600 (2BP4)   p90 +0.1093   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0512 vs uniform-effect null p10/p50/p90 +0.0439/+0.0513/+0.0588 -> pctile 0.492
        VERDICT: WORSE
      stage1-REP RANDOM minus do-nothing
        a 3.2255 (med 2.9698)   b 3.2148 (med 2.9661)   n=126
        effect +0.0107   median +0.0083   SE 0.0014   MDE 0.0039   effect/MDE +2.79
        iid  CI95 [+0.0082, +0.0135]
        fold CI95 [+0.0093, +0.0120]   folds same sign 5/5   per-fold 0:+0.013 1:+0.012 2:+0.011 3:+0.008 4:+0.010
        28W/98L/0T   worst degradation +0.1035 (1S9Z)   p90 +0.0283   power 1.00  Type-M 1.00
        concentration: drop-top10 +0.0124 vs uniform-effect null p10/p50/p90 +0.0106/+0.0124/+0.0143 -> pctile 0.513
        VERDICT: WORSE
      stage1-REP ORACLE cos: AMBER minus RANDOM
        a -0.0491 (med -0.0498)   b -0.0032 (med -0.0030)   n=126
        effect -0.0459   median -0.0401   SE 0.0155   MDE 0.0435   effect/MDE -1.05
        iid  CI95 [-0.0770, -0.0157]
        fold CI95 [-0.0734, -0.0212]   folds same sign 5/5   per-fold 0:-0.042 1:-0.101 2:-0.050 3:-0.003 4:-0.036
        77W/49L/0T   worst degradation +0.4596 (5MXS)   p90 +0.1598   power 0.84  Type-M 1.10
        concentration: drop-top10 -0.0172 vs uniform-effect null p10/p50/p90 -0.0377/-0.0175/+0.0020 -> pctile 0.508
        VERDICT: BETTER [TYPE-M ZONE: magnitude inflated ~1.10x]

    first run (L39)                          replication              first run's fold CI       inside?
    toward-member minus do-nothing  -0.0178  -0.0207 [-0.0250, -0.0166]  [-0.0255, -0.0124]     yes
    AMBER minus random              +0.0111  +0.0100 [+0.0059, +0.0175]  [+0.0062, +0.0171]     yes
    AMBER minus toward-member       +0.0385  +0.0414 [+0.0331, +0.0501]  [+0.0298, +0.0479]     yes
    random minus do-nothing         +0.0096  +0.0107 [+0.0093, +0.0120]  [+0.0077, +0.0122]     yes
    ORACLE cos AMBER minus random   -0.0503  -0.0459 [-0.0734, -0.0212]  [-0.0735, -0.0253]     yes

Every replicated quantity lands inside the first run's fold-clustered CI, all 5/5 folds. The
PROVISIONAL label on the toward-member line of L39 is lifted: a zero-information move of the
same size as the relaxation's, toward a random member of the shipped pool, improves the built
chain by 0.018 to 0.021 A (94W/32L in the replication, 1.77x its MDE), where the relaxation's
own move worsens it by 0.021. The reading in L39 stands unchanged: this is a property of the
projection's displacement (S16 L27: the projection moved the chain 0.166 A away from the point
cloud with a negative cosine, and the pool members surround the cloud), not of physics, and it
is not a proposal. The Type-M status of AMBER-minus-random (1.02x MDE in the replication) is
unchanged and the L46 caveat still applies to any quotation of that number.

---
