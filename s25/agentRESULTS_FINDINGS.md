# SPRINT 25 — RESULTS / RENDERER LANE. DESIGN NOTES AND FINDINGS.

Status: **INFRASTRUCTURE COMPLETE. FOUR RELEASE GATES LIVE AND TESTED.
NO REAL RESULT DATA GENERATED — `results/` carries no RMSD of any kind.**
Awaiting the physics lane's seven-configuration output format before the frozen build.

---

## 0. WHAT EXISTS, AND HOW TO RUN IT

| file | what it is |
|---|---|
| `s25/resultslab/exportlib.py` | export layer: T-map, PDB writer/reader, chain builder, round-trip verifier, **all four gates** |
| `s25/resultslab/providers.py` | architecture-agnostic registry: configuration LABEL → callable returning a Cα trace |
| `s25/resultslab/schema.py` | results data model, leaderboard, the three build-time refusals |
| `s25/resultslab/site.py` | static-site generator |
| `s25/resultslab/site/{index.html,style.css,app.js}` | browse flow and 3D viewer |
| `s25/resultslab/build.py` | the one command |
| `s25/resultslab/test_export.py` | 20 python tests |
| `s25/resultslab/site/test_site.js` | browser-side geometry, basis-key and header tests |

```
python -m s25.resultslab.exportlib                          # build/verify the T-number map
python -m s25.resultslab.test_export                        # 20 tests, ~90 s
python -m s25.resultslab.test_export --full                 # + the 126-target round trip
python -m s25.resultslab.build --mode selftest --limit 20   # SANDBOXED end-to-end proof
node    s25/resultslab/site/test_site.js                    # browser-side tests
python -m s25.resultslab.build --mode frozen --spec <spec>  # the real build (needs a spec)
```

`--mode selftest` writes to `s25/resultslab/_selftest/` and to the system temp directory.
**It cannot write to `results/`.** `--mode frozen` is the only mode that can, and it requires
a spec file that does not exist yet.

---

## 1. THE FOUR RELEASE GATES — ALL LIVE, ALL TESTED, ALL DEMONSTRATED FIRING

### (a) Pool-oracle gate — `exportlib.pool_oracle_gate`
Mean RMSD must be ≥ the mean ORACLE best member of each target's own K=500 candidate pool,
1.7108 Å (recomputed here from `u["rr"][pool_idx].min()`, cached to
`results/summary/pool_best.json`). Aggregate violation is HARD and stops the build.

Per-target violations are reported as a WARNING only, and the distinction is not squeamishness:
the emitted structure is a coordinate **average** over retained candidates and is therefore not
itself a pool member, so it *can* legitimately beat every single window on one target. It
cannot do that systematically across 126. The incumbent shows exactly this signature — `WARN`,
with a handful of per-target violations and a healthy aggregate margin.

**Confirmed weakest of the four, as flagged.** At σ = 2.4 the synthetic arm clears it with only
a WARN.

### (b) Difficulty gate — `exportlib.difficulty_gate` — THE ONE THAT DISCRIMINATES
`corr(per-target RMSD, that target's ORACLE pool_best) ≥ 0.30` **and** `sd(per-target RMSD) ≥ 0.70`.

Calibration re-measured independently on this instrument, n=126:

| arm | corr | sd | mean |
|---|---|---|---|
| incumbent, point cloud (`avg_ca`) | **+0.6988** | 1.6466 | 3.0483 |
| incumbent, fit arm (`fit_ca`) | +0.7126 | 1.7264 | 3.2041 |
| **incumbent, PRODUCTION chain (`ca`)** | **+0.7157** | **1.7361** | **3.2148** |
| synthetic σ=0.9 | −0.1201 | 0.1951 | 1.4259 |
| synthetic σ=1.6 | +0.0765 | 0.3313 | 2.5246 |
| synthetic σ=2.4 | +0.0193 | 0.4656 | 3.8043 |

(The coordinator's 1.6400 is the same statistic at `ddof=0`; I report `ddof=1`.)

The separation is roughly 5× the floor in both statistics, so fixed floors are defensible.
**Its limitation is stated in the code beside it**: it is calibrated on ONE genuine arm, and a
real method that were uniformly accurate across easy and hard targets would trip it. A failure
here is a claim to investigate, not a verdict to hide — but nothing ships until a human has looked.

### (c) Provenance-in-the-header gate — `exportlib.verify_header` — THE PRIMARY GATE
Every exported PDB carries `REMARK 999 PROVENANCE / MODULE_SHA / GIT_COMMIT / BASIS`, plus
target, configuration, RMSD, referent, Cα-bond statistics, and — on anything not genuine —
a `WARNING *** NOT A RESULT — TEST FIXTURE, SYNTHETIC ***` line. `verify_header` reads the
file **back off disk** (not the dict that was meant to be written) and the frozen build accepts
only `genuine`. Absent provenance is `unknown` and fails. `providers.from_json`/`from_npz` no
longer default to `False`.

The test also strips the `PROVENANCE` line out of a written file and asserts the gate notices.

### (d) Fixture quarantine
`FIXTURE_DIR` is the **system temp directory**, not `results/structures/_fixtures/` — a
directory beside the real exports is one `cp -r` from shipping. Fixture labels are
`_zero`, `_noise`, `_incumbent`, and `providers.register` **refuses outright** to bind a
synthetic provider to a real configuration name, so the proof build's failure mode is no
longer expressible. Sandboxed builds also keep their *natives* out of `results/`.

### The gates demonstrated at FULL SCALE, n=126, sandboxed selftest
```
label        provenance pool   diff        corr       sd  expected
_incumbent   genuine    WARN   PASS     +0.7157   1.7361  PASS OK
_zero        synthetic  FAIL   FAIL     +0.3692   0.2278  FAIL OK
_noise       synthetic  WARN   FAIL     -0.1504   0.5215  FAIL OK
gates discriminate as designed
```
and the leaderboard the same run produced:
```
configuration     n      mean   median  pt-cloud    effect     MDE  pool  diff    corr
_zero        126/126   0.9904   0.9777    1.4354   -2.2244  0.4193  FAIL  FAIL  +0.369
_noise       126/126   2.3955   2.3126    3.7622   -0.8193  0.4842  WARN  FAIL  -0.150
_incumbent   126/126   3.2148   2.9662    3.0483        --      --  WARN  PASS  +0.716
```
**The end-to-end path reproduces the frozen headline exactly**: the incumbent replay comes out
of export -> schema -> leaderboard -> site at **3.2148 A built chain / 3.0483 A point cloud**.
And it is the clearest statement of why gate (a) alone is not enough: `_noise` is a bag of
isotropic noise, it posts a **-0.82 A "improvement" that clears the pool gate with a WARN and
reads BETTER on the paired statistics**, and only the difficulty gate calls it (corr -0.150).
At n=126 `_noise`'s sd rises to 0.5215 -- still below the 0.70 floor, but closer than at n=20,
which is the honest reason the gate needs BOTH statistics rather than sd alone.
`build.py::_assert_gates_discriminate` makes this an assertion, not a printout: **a gate that
never fires is not evidence**, so the selftest fails if the gates stop separating the replay
from the noise. Selftest gates run SOFT on purpose — a hard gate would abort before the table
that demonstrates them exists. The frozen build runs them hard.

---

## 2. THE ROUND TRIP — FULL INSTRUMENT, PRODUCTION BASIS, n=126

`python -m s25.resultslab.test_export --full`:

```
all_ok                              True
worst |file - instrument|  chain    2.641e-04 A   (tolerance 2.0e-03)
worst |file - instrument|  cloud    2.443e-04 A
reprojection vs the cached arm      0.000e+00 A   <- bit-exact stage 3b
mean RMSD  BUILT CHAIN (production) 3.2148 A
mean RMSD  point cloud (non-phys.)  3.0483 A
basis gap (chain - cloud)           +0.1664 A
mean virtual CA-CA bond   native    3.8122 A
                          chain     3.8039 A
                          cloud     2.9614 A   <- 22.3% contracted
```

Every frozen number reproduced independently, through this lane's own code path. The
quantisation bound is `sqrt(3) × 5e-4 = 8.66e-04 Å` (Kabsch is 1-Lipschitz in the point set),
so the observed error is PDB's `%8.3f` field and nothing else.

**The projection is bit-exact against the shipped arm** — `chain_from_ca(lam=0.3)` reproduces
the cached `ca` field to 0.000e+00 on all 126. The chain builder in this package *is*
production stage 3b, not a re-implementation of it.

Browser side, over both bases: worst |Δ| between the browser's independent Horn-quaternion
recomputation and Python's = **3.0e-07 Å** (the header's `%.6f`), and a mirrored native does
not superpose (3.4114 Å) — proper rotations only.

Independent measurement of the basis ruling's evidence, n=126: point-cloud mean virtual bond
**2.9614 Å** against native **3.8122 Å**, worst single bond **0.6492 Å**, **80/126** targets
under 3.4 Å. Every figure matches L4/L8 exactly.

---

## 3. BASIS — RULED, AND MADE STRUCTURALLY UNMIXABLE

**Production = BUILT CHAIN at lam=0.3 (SYNTHESIS stage 3), 3.2148 Å.** Point cloud 3.0483 Å is
retained everywhere as an explicitly labelled non-physical intermediate. Reported always paired.

`export_pair()` emits both bases from one call, so **there is no code path that yields one
number without the other**. Beyond that:

* separate directories — `results/structures/chain/` vs `results/structures/`;
* `REMARK 999 BASIS` and `BASIS_ROLE` in every file (`PRODUCTION RESULT` /
  `NON-PHYSICAL INTERMEDIATE -- not a molecule`);
* `basis` is a **required** schema field — `schema._require_basis` raises on a present record
  that carries an RMSD without one;
* the site's structure keys are `T001__legacy__built_chain_bb`, so the viewer physically
  cannot load one basis and print the other's number;
* the viewer prints the basis in the header badge, in the large RMSD label, in the meta row,
  in every table heading, and on the basis toggle, which moves the structure and the number
  together.

**Two corrections to numbers I was given.** The brief's 3.2126 Å does not reproduce anywhere in
this repository and nothing here quotes it. The +0.156 gap is the **lam=0 fit arm**; the
SYNTHESIS arm's gap is **+0.1664**.

---

## 4. EXISTING CODE FOUND AND REUSED

**Reused**
* `core/geometry.py:847 write_pdb` — `exportlib._atom_lines` is asserted **byte-identical** to
  its ATOM body; a local writer exists only because `write_pdb` takes a single `remark` string
  and the lab needs a parseable `REMARK 999` block.
* `core/geometry.py:148 build_backbone`; `core/project.py` `lam_path`/`make_penalty`/
  `build_ca_exact` — production stage 3b, called directly (rather than via
  `instrument.project`) because the arm's **torsions** are needed to build a backbone.
* `s12/instrument.py` — `targets`, `load_univ`, `ca_rmsd`, `pool_idx`, `shipped_record`.
  **`ca_rmsd` is the only RMSD in this package.**
* `s24/stats_lib.py` — `compare`, `provenance`, `save_atomic`. No statistic is hand-rolled.

**Deliberately not reused:** `core/geometry`'s reader drops any residue missing N, CA or C by
design, which is every residue of a Cα-only trace; `exportlib.read_pdb` is a minimal column
reader, asserted against `parse_pdb_ensemble` where both apply.

**There was no visualisation layer in the repo** — no HTML, JS or CSS anywhere; only
`s13/figures.py` (matplotlib). The 3D renderer is new; nothing was duplicated.

---

## 5. THE DATA MODEL

`results/summary/results.{json,csv}`, one record per (target, configuration), full key set
enforced row by row by `ST.save_atomic(complete_keys=RECORD_KEYS)`. Each row carries **both
bases**: `rmsd`/`basis`/`prediction_path` and `rmsd_secondary`/`basis_secondary`/
`prediction_path_secondary`, plus `basis_delta`, the Cα-bond statistics for both bases and the
native, `pool_best`, `provenance`, `git_commit`, `module_hash`, and the metadata the viewer shows.

**Absent means absent.** Every grid cell produces a record; a cell with no prediction gets
`status: "absent"`, `rmsd: null` and an `absent_reason`. A provider that raises produces an
absent record carrying the exception text — not a crash, not a silent skip. `export_prediction`
rejects a wrong-length or non-finite trace rather than writing something plausible.

`leaderboard.{json,csv}` adds, from `stats_lib.compare` paired against the baseline: effect,
median effect, SE, **MDE = 2.8016 × SE**, effect/MDE, iid CI, **fold-clustered CI beside it**,
per-fold effects, folds-same-sign, W/L/T, worst degradation, Type-M, verdict — plus both gate
statuses and `corr_with_pool_best`. A row whose fold CI includes zero renders as `NOT MEASURED`.

**Three build-time refusals** in `schema.build`: undeclared provenance; a synthetic record
outside a selftest; and the gates.

---

## 6. THE RENDERER

Static, self-contained, no server; 3Dmol.js **2.5.5 from cdnjs, pinned exactly**, with a
graceful message if it cannot be reached (no number on the page depends on it). Structures are
inlined into `structures.js` because a `file://` page cannot `fetch()` a sibling file.

* Prediction as a thick blue Cα tube; native as a thin grey tube, visibly subordinate; legend.
* **Camera consistency**: each pane frames an invisible copy of the *native*, so every
  configuration of a target gets an identical camera — framing on the drawn shapes would
  rescale each pane by how far its own prediction wanders, which is exactly the comparison
  being made by eye.
* **RMSD is the largest thing on a target page**, with its basis, the paired other-basis value,
  and target / configuration / provenance / length / fold / selector / Hamiltonians / distogram /
  CVaR α / qubits / layers / candidates / retained beneath it.
* Browse flow: overview (headline, leaderboard with both bases and both gates, architecture,
  "How to check this", target index) → method (distribution, per-target table) → target
  (3D + overlay + basis toggle + provenance) → **`#/compare/<T>`, the same-target side-by-side**
  of every configuration, each with its own RMSD and the best cell outlined.
* Framing: the page is built so a reader can **check** the work — the "How to check this" card
  states that every number is recomputed from exported bytes, every structure carries its own
  provenance, every basis is stated, absent means absent, and the statistics come from
  `stats_lib`. Not a victory lap.
* Accessibility: `prefers-reduced-motion`, `:focus-visible`, skip link, `role="img"` +
  `aria-label` on every 3D pane and the histogram, wide tables scroll in their own container,
  light/dark via `prefers-color-scheme`.
* Verified in Chrome across overview / method / target / compare, both bases; no console errors.

---

## 7. DEFECTS FOUND AND FIXED IN MY OWN CODE

1. **`export_prediction` was not atomic.** It wrote twice — once to obtain the quantised
   coordinates, once to stamp the RMSD computed from them — so an interruption left
   `RMSD_CA PENDING` on disk (`T006__amber.pdb` did). Now `quantise()` obtains those
   coordinates **by formatting and reparsing in memory**, so the file is written once, complete.
   `np.round(x, 3)` would not do: numpy rounds half-to-even on the binary value while `%8.3f`
   rounds the decimal representation, and they disagree on ties. Asserted bit-equal to the file.
2. **A hash-seeded RNG is not deterministic.** The synthetic provider seeded from
   `hash((pdb, tag, seed))`; Python **salts string hashing per process**, so two runs of the
   "same" build emitted different structures. Now `zlib.crc32`, with a test that launches
   subprocesses under two `PYTHONHASHSEED` values and asserts bit-identical output.
   Worth generalising: **any seed derived from `hash()` of a string is a reproducibility hole**,
   and this repo should assume it is present wherever that pattern appears.
3. **Per-record timestamps made determinism unverifiable.** One timestamp per build now;
   content hashes identically across consecutive builds.
4. **`os.replace` loses a race on Windows.** A build writing ~2,000 small files hits the
   on-access virus scanner holding a just-created file; `os.replace` then fails with WinError 5.
   Bounded retry, which is correct because the write succeeded and only the rename raced.
5. **The fixture labels were `distogram` and `legacy`** — so this module's own test output
   contained `T001  distogram  14  0.000000  0.000000`. Renamed to `_zero`/`_noise` and
   quarantined; see gate (d).

---

## 8. THE QUARANTINED PROOF BUILD — WHAT I TAKE FROM IT

The defect was narrow and exactly as diagnosed: **the label did not reach the artefact.**
`is_synthetic` was correct in `results.json`, correct in the leaderboard, printed on the page
in three places — and absent from all 1,134 PDB headers. Separated from its JSON, every one of
those files read as a genuine result, and `T001__distogram.pdb` read as a perfect one.

Three structural changes, in increasing order of how little they depend on anyone remembering:

1. the label now travels **inside the file** (gate c), so separation cannot strip it;
2. a synthetic provider **cannot be bound** to a real configuration name, so the artefact
   cannot acquire a misleading name in the first place;
3. the **difficulty gate** does not depend on labels at all — it reads the numbers and asks
   whether they behave like a method.

---

## 9. NOTES FOR THE COORDINATOR

* **Waiting on the physics lane's output format** before the frozen build, as instructed. What
  the lab needs per configuration: `{"<PDB id>": [[x,y,z], ...]}` — the **emitted point cloud**,
  residue order matching the sequence — plus a `meta` block that must include
  `is_synthetic: false`. The lab does the stage-3b projection itself at `lam=0.3`, so every
  configuration is chained by the same operator. Anything the lane emits differently, tell me
  and I will add a loader; nothing else in the package changes.
* **I have invented no configurations.** `CONFIGURATIONS` is a list of eight permitted LABELS
  with no meaning attached; a label with no registered provider produces nothing at all.
* **`results/` currently contains no RMSD**: `benchmark_manifest.json` and
  `monomer_manifest.json` (pre-sprint, untouched), `summary/target_map.json` (pure mapping),
  `summary/pool_best.json` (ORACLE pool minima, gate input, no prediction), and
  `summary/professor_brief.md` (another lane's, untouched).
* **`s24/stats_lib.py` is being edited by another lane while I import it** — it grew from 384 to
  ~490 lines mid-session and `_verdict` now emits a more detailed NOT MEASURED string. My use
  (`compare`, `provenance`, `save_atomic`) is unaffected and I have not touched the file.
* I have taken **no lock**, touched **no benchmark60 target**, and generated **no result data**.

---

## 10. OPEN DEFECT FOUND WHILE A FROZEN BUILD WAS IN FLIGHT — NOT YET FIXED

At 20:37 another lane launched `--mode frozen` against this lab (log
`s25/results/frozen_build.log`, empty because Python buffers a redirected stdout). It is
running now. Its files verify: `MODULE_SHA 90e160cc06020ac1` matches `exportlib.py` on disk
exactly, so no source drift; `PROVENANCE genuine`; `BASIS built_chain_bb`; `CHAIN_LAM 0.300`.
Read-only spot check of four of its own outputs, recomputed from the bytes against the
exported native:

```
T001  header 2.738867  from files 2.738867  |d| 9.3e-08
T005  header 3.993981  from files 3.993981  |d| 3.0e-08
T050  header 6.592727  from files 6.592727  |d| 2.1e-07
T100  header 1.492276  from files 1.492276  |d| 2.6e-07
```

**I have stopped editing every module in `s25/resultslab` for the duration** (BRIEF SS4: never
edit a module while a job launched from it is still running; an artefact was written from
vanished source here once already).

### THE DEFECT: the gates protect the summary, not the files

`schema.collect()` writes structures to `results/structures/` as it goes; the pool-oracle and
difficulty gates run afterwards, in `leaderboard()`. So a configuration that FAILS a gate has
already left 252 provenance-stamped PDBs in `results/`, and the build then raises — leaving
files on disk with no `results.json` beside them. **That is precisely the quarantined failure
mode**: a directory of plausible-looking structures separated from the JSON that would have
disqualified them. Being `PROVENANCE genuine` does not help; the gate's whole point is that a
file can be genuinely produced and still be wrong.

**Proposed fix, needs the coordinator's go-ahead because the module is in use:** collect into a
STAGING directory (`results/.staging/`), run every gate, and only then promote to
`results/structures/` with an atomic directory move. A build that fails a gate then leaves
nothing behind but a log. Second-best if staging is unwanted: gate each configuration
immediately after collecting it and before the next one starts, which still leaves the failing
configuration's own files behind but stops the run one configuration earlier.

Until that is fixed, **a failed frozen build must be treated as leaving debris**: check
`results/structures/` against `results/summary/results.json` and delete any structure the JSON
does not account for.
