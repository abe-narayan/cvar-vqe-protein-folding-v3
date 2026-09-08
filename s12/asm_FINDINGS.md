# ASSEMBLY AND RETRIEVAL — s12 findings

**Question.** Is whole-window retrieval the wrong granularity? Can multi-piece fragment
assembly beat it, and is the placement problem tractable without the native?

**Headline.** All three reconnaissance numbers replicate (0.45 / 0.29 Å two-piece rigid
floor, +0.338 Å placement cost) — and the direction is **refuted anyway**. A synthetic bank
of the same size, built from i.i.d. Ramachandran draws with no fragment content at all,
reaches the same floor (rigid 0.328 vs 0.309; chain-consistent 0.675 vs 0.647, CI crosses
zero). The k ≥ 2 oracle ladder measures the fitting capacity of the placement parameters,
not the library. The deployable version is **+0.218 Å worse** than the shipped synthesis
(CI [+0.106, +0.336], 46W/80L), for the record's own reason: assembly improves the candidate
set (oracle best 1.711 → **1.175 Å**) and the shipped objective then chooses worse from it
(top-75 best 2.306 → 2.871 Å; the best assembly sits at the **74th** score percentile on
FAIL18). One positive worth keeping: **assembly is the only lever in the record that moves
the deployable pool's achievable floor**, and it becomes valuable the moment — and not
before — a ranker with in-band skill exists.

---

## 0. Instrument validation (run first, before any claim)

| check | result |
|---|---|
| sub-length bank machinery at L = n vs `s8/generate_univ/<pdb>.npz` | coordinates **bit-identical** on 12 spot targets (`max|ΔW| < 1e-4`), window counts identical |
| my analytic-SVD Kabsch vs `I.kabsch_rmsd_batch` | max |Δ| **6.9e-7 Å**, 7× faster |
| my `k=1` full-universe oracle over the rebuilt banks | **1.313 Å** — reproduces the record's "whole-universe best window" exactly |
| bank recipe | `peptide_db.load()` filtered by `folds(5)[seq] != fold(target)` + `core.predict._fold_fragments(fold, 5)`; identical iteration order to `s8/generate.py:_windows_all` |

Banks are keyed by (fold, L) only. That is leakage-safe: the target's own sequence lives in
its own fold and is removed by the fold filter, which is why the L = n bank reproduces the
shipped universe exactly.

Code: `s12/asm_lib.py`. Banks cached in `s12/cache/asm_banks/` (5 folds × L = 3…16).

---

## 1. The rigid-placement oracle ladder (ORACLE, upper bound, NOT one chain)

Each piece is placed by its **own** optimal rigid transform onto its own native segment.
Because the pieces are independent, `n·RMSD² = Σ_pieces SSD(piece, native segment)`, so the
optimum over cut positions is an exact DP over segment costs — no approximation anywhere.

Candidate set per segment (the DP is identical across arms; only the bank differs):
`full` = every library window of that length; `blos500/100` = top-N by BLOSUM62 vs the
target's **sub**-sequence (deployable retrieval); `rand500/100` = N uniform-random windows.

**Minimum piece length 4, "at most k pieces", mean CA-RMSD over the 126 tuning targets:**

| arm | k=1 | k=2 | k=3 | k=4 |
|---|---|---|---|---|
| full     | **1.313** | **0.309** | 0.102 | 0.087 |
| blos500  | 1.711 | 0.541 | 0.255 | 0.233 |
| blos100  | 1.970 | 0.674 | 0.365 | 0.341 |
| **rand500** | 1.795 | **0.575** | 0.254 | 0.228 |
| rand100  | 2.155 | 0.752 | 0.380 | 0.343 |

FAIL18 / other-108 split:

| arm | k=1 F18 / o108 | k=2 F18 / o108 | k=3 F18 / o108 |
|---|---|---|---|
| full    | 1.640 / 1.259 | 0.453 / 0.285 | 0.119 / 0.099 |
| blos500 | 2.284 / 1.615 | 0.771 / 0.502 | 0.293 / 0.248 |
| rand500 | 2.281 / 1.714 | 0.806 / 0.537 | 0.292 / 0.248 |

With minimum piece length 3 the ladder collapses further: full k=4 = 0.018 Å, k=5 = 0.011 Å.

**Does the reconnaissance claim replicate?** Approximately yes, on the rigid convention.
Claimed ≈ 0.53 Å on the FAIL18 and ≈ 0.32 Å on the rest; I measure **0.453 / 0.285** with
the full library and **0.771 / 0.502** with a 500-window retrieval per segment. The claim
is in the right band and is real *as arithmetic*. What it is not is evidence of anything —
see next.

---

## 2. The decisive control: random pieces reach the same floor

Hypothesis under test: the two-piece gain comes from shorter pieces having better matches
(retrievable information). Null: the gain is the capacity of extra rigid degrees of freedom.

| | k=1 | k=2 | k=1→k=2 gain |
|---|---|---|---|
| BLOSUM top-500 per segment | 1.711 | 0.541 | **−1.170** |
| **500 uniformly random windows** | 1.795 | 0.575 | **−1.220** |
| difference (retrieval − random) | −0.084 | **−0.034** | +0.050 |

Random pieces reproduce **96 %** of the two-piece "gain". Under chain-consistency
(§3b) the same contrast is 1.030 vs 1.044, paired d = −0.014 Å, CI [−0.053, +0.023],
64W/62L — an exact null. Sequence-based retrieval of *pieces* is worth 0.034 Å at k=2 — below the sprint's 0.03 Å ignore threshold once you note
that random pieces give a *larger* gain than retrieved ones (they start worse at k=1 and
end in the same place).

### 2b. The capacity null, run properly (`asm_null.py`)

The §2 contrast varies *which* real windows are offered. The stronger control replaces the
library entirely: a synthetic bank of the **same size** for every (fold, length), whose
pieces are the ideal-geometry rebuilds of (phi, psi) drawn **i.i.d. from the pooled
Ramachandran of the same legal library**. It has the right residue-conformation marginal
and nothing else — no real fragment, no sequence-structure relation, no residue-to-residue
correlation. Same DP, same minimum piece length.

| k | real library | capacity null | what the library is worth |
|---|---|---|---|
| 1 | **1.313** | 1.617 | **−0.304** |
| 2 | 0.309 | 0.328 | −0.019 |
| 3 | 0.102 | 0.111 | −0.009 |
| 4 | 0.087 | 0.094 | −0.007 |

The real fragment library is worth 0.304 Å as a source of whole windows and **0.019 Å** as
a source of assembly pieces. Structureless random torsion pieces reproduce 94 % of the
two-piece oracle floor and 99 % of the three-piece one. That is the refutation: the k ≥ 2
oracle ladder measures the fitting capacity of the placement parameters, not the library.
(The same holds under chain-consistency, §3b: 0.647 real vs 0.675 null.)

Parameter counting says the same thing: k rigid pieces carry 6k parameters, a target has
3n ≈ 39 coordinates and 6 are global, so k = 4 pieces (24 parameters) is already
over-parameterised — and indeed k = 4 reads 0.087 Å from the full library and 0.228 Å from
*random* 500-window banks. **The rigid-placement ladder is a capacity curve, not a
structural result, and must not be quoted as an achievable bound** (same failure mode as
correction C4, the convex-hull ladder).

---

## 3. Candidate-count scaling: cardinality is the wrong currency (`asm_scaling.py`)

E[min CA-RMSD | s candidates], averaged over the 126 targets, oracle, rigid placement.
`k=1` subsamples the length-n bank; `k=2`/`k=3` subsample every piece bank to s and run the
same exact DP. "log10 eff" is log10 of the number of distinct assemblies searched.

| family | s | log10 eff | mean floor |
|---|---|---|---|
| k=1 | 30 | 1.48 | 2.493 |
| k=1 | 100 | 2.00 | 2.120 |
| k=1 | 300 | 2.48 | 1.882 |
| k=1 | 1000 | 3.00 | 1.650 |
| k=1 | 3000 | 3.48 | 1.493 |
| k=1 | full (≈16k) | 4.21 | **1.313** |
| k=2 | 10 | 2.74 | 1.043 |
| k=2 | 100 | 4.74 | 0.684 |
| k=2 | 1000 | 6.74 | 0.486 |
| k=3 | 100 | — | 0.346 |
| k=3 | 1000 | — | 0.212 |

Two things follow.

1. **Growing the window library is hopeless.** The k=1 coverage curve is −0.432 Å per
   decade of windows. From today's 1.313 Å, reaching 0.5 Å needs ≈ 10^6.1 windows — ~800×
   the library. And this curve is *optimistic*: it subsamples an existing library, so its
   samples are less redundant than the windows a 3.2× bigger corpus actually adds (which
   moved universe best by only 0.049 Å, S-record). Depth is closed, and this prices it.
2. **Assembly is not just more candidates.** At a *matched* 10^4.2 assemblies searched,
   k=2 interpolates to ≈ 0.77 Å against k=1's 1.313 Å. So the k≥2 families are not on the
   k=1 curve. That is exactly what you expect from a product set with 12 free rigid
   parameters rather than 10^4.2 independent samples — cardinality is the wrong axis;
   degrees of freedom is the right one. Which is why the null in §2 and §4 is the test
   that matters, not the count.

---

## 3b. Chain-consistency: the placement cost, and the null again (`asm_chain2.py`)

A **chain-consistent** assembly is ONE ideal-geometry chain: residue i takes (phi_i, psi_i)
from whichever piece covers it, the torsion vector is concatenated, and the CA trace is
rebuilt by `build_ca_exact`. There are ZERO free placement parameters. Search: per cut,
the top-M = 150 pieces per side by their own ideal-rebuild SSD (an admissible lower bound
on any assembly containing that piece), then all M² concatenations rebuilt and measured
exactly. Sanity: exact branch-and-bound (`asm_chain.py`) gives 1A13 = 0.782 vs M=300's
0.78 and 1A1P = 0.978 vs M=150's 1.04, so M=150 loses ≲ 0.06 Å and loses it in every arm.

| arm | k=1 | k=2 | k=2 FAIL18 | k=2 other-108 |
|---|---|---|---|---|
| **real** (full library) | 1.313 | **0.647** | 0.875 | 0.609 |
| **rama** (capacity null, matched size) | 1.617 | **0.675** | 1.050 | 0.613 |
| blos500 (deployable-sized retrieval) | 1.721 | 1.030 | 1.376 | 0.972 |
| rand500 (retrieval null) | 1.803 | 1.044 | 1.389 | 0.986 |

Three readings.

1. **The placement cost is +0.338 Å** (rigid k=2 0.309 → chain-consistent k=2 0.647). The
   reconnaissance's "placement floor ≈ 0.36 Å" replicates. So the two-piece floor does NOT
   collapse under chain-consistency — it roughly doubles but stays far under the 1.313 Å
   whole-window floor.
2. **The capacity null survives chain-consistency.** A synthetic bank of the same size
   whose pieces are i.i.d. Ramachandran draws — no real fragment, no sequence-structure
   relation, no residue-residue correlation — reaches 0.675 Å against the real library's
   0.647 Å. The whole library is worth **0.028 Å** at k=2 (it is worth 0.304 Å at k=1).
   Whatever the two-piece floor measures, it is almost entirely the expressive power of
   "two torsion blocks joined at a free junction", not the contents of the library.
3. **Retrieval is worth nothing at the piece level**: blos500 1.030 vs rand500 1.044.
   BLOSUM retrieval of *pieces* buys 0.014 Å where it buys 0.082 Å on whole windows.
   Sub-sequences are too short to carry retrieval signal — consistent with the record's
   "structure and sequence are decoupled".

**k = 3 chain-consistent** (40-target subset, M = 40 per slot, real vs capacity null;
29 targets have n ≥ 12 so three ≥ 4-residue pieces fit): k=2 0.675 / k=3 0.748 (real),
0.738 / 0.775 (null); best-of-(k2,k3) 0.600 real vs 0.634 null, paired d = −0.035 Å,
CI [−0.111, +0.040], 19W/10L. k=3 does **not** beat k=2 here — at M = 40 the 3-piece
search is far weaker than the 2-piece one, so this is a search-budget statement, not a
statement about the 3-piece space. What does carry through is the null: real and synthetic
banks are again indistinguishable.

### 3c. A second join convention agrees (`asm_overlap.py`)

Overlapping pieces with a shared-residue join (programme item d), done properly: piece 1
covers [0, c+3), piece 2 covers [c, n), they share 3 residues, and piece 2 is placed by
superposing its first 3 residues on piece 1's last 3. Both pieces keep their **real**
deposited geometry; o = 3 makes the join fully determined and it needs no native.

| join convention | k=2 oracle | FAIL18 | other-108 |
|---|---|---|---|
| free rigid placement (not a chain) | 0.309 | 0.453 | 0.285 |
| torsion concatenation (ideal chain) | 0.647 | 0.875 | 0.609 |
| **3-residue overlap superposition** | **0.649** | 0.873 | 0.612 |

Two independent, physically realisable join rules land on the same number, so the +0.34 Å
placement cost is a property of the problem, not of one convention. One caveat in favour
of torsion concatenation: the overlap join's emitted virtual CA-CA bonds run 3.36–4.54 Å
(ideal 3.80; the database gate is 3.5–4.1), i.e. it sometimes emits a broken chain, whereas
torsion concatenation is ideal-geometry by construction.

Incidental: the ideal-geometry rebuild costs **nothing** at k=1 — the best rebuilt window
is 1.3134 Å against the best real-coordinate window's 1.3134 Å (mean paired difference
−0.0000 Å, max |Δ| 0.235 Å; different windows win). Torsion-space working is free here.

## 4. Coverage, not size (`asm_cover.py`, `asm_cover_ss.py`) — the crux

For every native sub-segment of length m, the minimum CA-RMSD over EVERY length-m window
of the fold's legal library (ORACLE). Shape-DOF = 2m−5: the free coordinates of an m-mer CA
trace once the 3.8 Å virtual bonds are fixed.

| m | shape-DOF | bank size | mean *worst* hole per target | mean hole | frac of segments > 0.5 Å |
|---|---|---|---|---|---|
| 4 | 3 | 72–79k | 0.070 | 0.035 | 0.000 |
| 5 | 5 | 66–72k | 0.224 | 0.129 | 0.000 |
| 6 | 7 | 59–66k | 0.393 | 0.253 | 0.085 |
| 7 | 9 | 52–59k | 0.539 | 0.383 | 0.312 |
| 8 | 11 | 46–52k | 0.696 | 0.520 | 0.481 |

Coverage by native secondary structure (`I.ss_of` on the native torsions — DIAGNOSTIC):

| class | 6-mers (n, mean) | 8-mers (n, mean) |
|---|---|---|
| pure helix | 352, **0.119** | 236, **0.199** |
| pure coil  | 342, 0.342 | 219, **0.754** |
| mixed      | 309, 0.307 | 296, 0.601 |
| N-term / mid / C-term (m=8) | — | 0.581 / 0.496 / 0.552 |

**Answer to the crux.** It is *not* true that the library covers every local motif and only
assembly is missing. The library covers 4- and 5-mers essentially exactly (0.04–0.13 Å) —
but that is a dimension-counting fact (3–5 shape DOF sampled 70,000 times), not a property
of this corpus. Real holes open at m = 7–8 (9–11 DOF), they are **concentrated in coil and
irregular segments** (0.754 Å for pure-coil 8-mers vs 0.199 Å for helical ones), they are
**not** at the termini, and the worst 8-mer hole per target correlates **r = +0.815** with
that target's whole-window oracle floor. Helices are solved; irregular local structure is
not, and the missing whole-window accuracy is largely the missing *local* accuracy at
7–8 residues. FAIL18 targets have systematically worse holes (mean-of-worst 0.873 Å at
m = 8 vs 0.696 overall).

So: "assembly, not coverage" is the wrong dichotomy. At the length where holes appear
(7–8), covering the motif *is* the assembly problem, and it is a coverage problem in an
11-dimensional space that no library of 10^4–10^5 fragments from this corpus will fill.

## 5. Deployable assembly — no native in any decision (`asm_deploy.py`)

Generation: for every cut c, the top-200 pieces per side by BLOSUM62 against the target's
**sub**-sequence, all 200² torsion concatenations rebuilt with `build_ca_exact`
(80k–360k assemblies per target, mean 2.3 × 10^5). Scoring: the **shipped** leave-fold-out
distogram Bayes-risk score, i.e. exactly the production objective. Terminal operators:
argmin, and `I.coordinate_average` of the 75 best-scoring → `I.project`.
Instrument check: `pool_argmin` reproduces the record's 3.454 exactly.

| arm | mean | FAIL18 | other-108 |
|---|---|---|---|
| pool_argmin (shipped baseline) | **3.454** | 6.008 | 3.028 |
| pool1_argmin (same pool, ideal-rebuilt) | 3.507 | 6.039 | 3.085 |
| **asm_argmin** | **3.570** | 6.124 | 3.145 |
| mix_argmin (pool ∪ assemblies) | 3.559 | 6.124 | 3.132 |
| **asm_avg75 → project** | **3.422** | 6.017 | 2.990 |
| mix_avg75 → project | 3.421 | 6.017 | 2.989 |
| shipped synthesis (`fit_ca`) | 3.204 | — | — |
| *oracle* pool_best | 1.711 | 2.284 | 1.615 |
| *oracle* **asm_best** | **1.175** | 1.591 | 1.106 |
| *oracle* top-75-by-score best, pool | 2.306 | — | — |
| *oracle* top-75-by-score best, assemblies | **2.871** | 5.503 | 2.432 |

Paired against the shipped argmin (3.454), and against the synthesis (3.204):

| arm | vs 3.454 | 95 % CI | W/L | drop-10 | vs 3.204 | 95 % CI | W/L |
|---|---|---|---|---|---|---|---|
| asm_argmin | **+0.116** | [−0.002, +0.235] | 58/68 | +0.224 | **+0.366** | [+0.226, +0.510] | 38/88 |
| mix_argmin | +0.105 | [−0.008, +0.221] | 58/68 | +0.212 | +0.355 | — | — |
| asm_avg75 | −0.032 | [−0.126, +0.062] | 66/60 | **+0.066** | **+0.218** | [+0.106, +0.336] | 46/80 |
| mix_avg75 | −0.033 | [−0.127, +0.061] | 66/60 | +0.065 | +0.217 | [+0.104, +0.335] | 46/80 |

**No deployable version beats the incumbent.** The −0.03 Å of `asm_avg75` against the
distogram argmin is inside its CI and reverses to +0.066 Å when the top ten targets are
dropped; measured against the real production terminal (synthesis, 3.204) the assembly
path is **+0.218 Å worse, CI excluding zero, 46W/80L**.

The mechanism is the record's own: **the candidate set improves and the objective cannot
use it.**

- deployable assembly improves the oracle best of the candidate set from 1.711 → **1.175 Å**
  (−0.536 Å, and −0.69 Å on FAIL18) — assembly *is* a real generation improvement;
- but the top-75 the shipped filter keeps gets **worse**, 2.306 → 2.871 Å (5.503 on FAIL18);
- the best assembly sits at the **30.5th** score percentile overall and the **74.2nd** on
  FAIL18 — on the hard targets the objective ranks the best structure worse than chance
  among 10^5 candidates.

Enlarging the hypothesis space under an objective with calibration slope 0.376 and in-band
ρ ≈ 0.13 costs more than the extra reach earns. This is the same law as the record's
"the objective does not rank the native / optimise harder, get worse", now measured at
10^5 hypotheses instead of 500.

## 6. Answers to the three questions I was set

**Does the two-piece oracle claim replicate?** *Yes, arithmetically.* Rigid placement over
the full library gives 0.453 Å on FAIL18 and 0.285 Å on the other 108 (claimed 0.53 / 0.32),
and the placement cost — the gap between free rigid placement and a single ideal chain — is
+0.338 Å (claimed ≈ 0.36 Å). All three reconnaissance numbers land.

**Does it survive the random-piece null?** *No.* A synthetic bank of the same size, built
from i.i.d. Ramachandran draws with no fragment, no sequence relation and no residue-residue
correlation, reaches 0.328 Å rigid (real 0.309) and 0.675 Å chain-consistent (real 0.647,
paired d = −0.028, CI [−0.066, +0.009], drop-10 +0.015). The library is worth 0.304 Å as a
source of whole windows and 0.019–0.028 Å as a source of assembly pieces. **The k ≥ 2 oracle
ladder is a capacity curve.**

**Does it survive chain-consistency?** *Partly.* It roughly doubles (0.309 → 0.647) but
stays far below the 1.313 Å whole-window floor, and two independent join conventions agree.
So chain-consistency is not what kills assembly — the null is.

**Does any deployable version beat the incumbent?** *No.* +0.218 Å against the synthesis
baseline, CI [+0.106, +0.336], 46W/80L; the argmin arm is +0.366 Å worse. Assembly makes the
candidate set genuinely better (oracle best 1.711 → 1.175 Å) and the shipped objective then
makes a worse choice from it (top-75 best 2.306 → 2.871 Å; best assembly at the 74th score
percentile on FAIL18).

## 7. What the coordinator should do with this — ranked

1. **Retire the assembly direction as a generation lever, and retire the reconnaissance
   two-piece floor as a bound.** Add it to the corrections list beside C4 (hull capacity):
   *any oracle over k independently-placed or freely-jointed pieces must be priced against a
   matched-size i.i.d.-Ramachandran bank before it is quoted.* One line, `asm_null.py`.
2. **Assembly is worth keeping as a POOL-BUILDER if and only if the ranking problem is
   solved first.** It is the only intervention in this sprint that moved the achievable
   floor of a deployable candidate set by half an ångström (1.711 → 1.175, and 2.284 → 1.591
   on FAIL18). Nothing else on the record moves pool best at all. The moment a ranker with
   real in-band skill exists, this is where the headroom is — but not one hour before.
3. **The library's real deficiency is 7–8-residue irregular local structure, not whole-window
   coverage and not size.** Holes are ≈ 0.20 Å for helical 8-mers and 0.754 Å for coil
   8-mers, and the worst 8-mer hole predicts the whole-window floor at r = +0.815. If new
   data is ever collected (the record says a fresh benchmark needs new data anyway), it
   should be selected for *irregular* local conformations, not for more peptides.
4. **Stop pricing generation by candidate count.** The whole-window coverage curve is
   −0.432 Å per decade and is itself optimistic; 800× the library buys 0.8 Å. Depth is
   closed, and now it is priced rather than asserted.
5. **Do not run this on dev24.** Nothing here earns a confirmation pass: the one arm that
   is not a null is negative.

## Negatives, in full

| claim tested | result |
|---|---|
| two-piece assembly beats whole-window retrieval (oracle) | true but null-explained |
| the gain comes from better sequence matches on shorter pieces | **refuted** (0.014–0.034 Å) |
| the gain is uniform vs concentrated on FAIL18 | roughly proportional; FAIL18 keeps its 1.4–1.5× penalty at every k |
| chain-consistency destroys the two-piece floor | **refuted** (+0.338 Å, not collapse) |
| overlapping joins beat torsion concatenation | **refuted** (0.649 vs 0.647) |
| k = 3 chain-consistent beats k = 2 | not shown (search-budget limited) |
| deployable assembly beats the shipped argmin | **refuted** (+0.116) |
| deployable assembly beats the synthesis | **refuted** (+0.218, CI excludes 0) |
| the library covers every local motif (so only assembly is missing) | **refuted** at m ≥ 7 |

## Leakage audit

- Every bank is built from `peptide_db.load()` filtered by `folds(5)[seq] != fold(target)`
  plus `core.predict._fold_fragments(fold, 5)` — the same library the shipped universe uses;
  verified bit-identical at L = n on 12 targets. Banks are keyed by (fold, L) only, which is
  safe because a target's own sequence is in its own fold.
- `nat_ca` enters only ORACLE/DIAGNOSTIC arms (§1, §2, §3, §4) and the evaluation of §5.
  The §5 generation and scoring read only the target sequence, BLOSUM62, the fold's library
  and the shipped leave-fold-out distogram. Native SS in §4 is labelled DIAGNOSTIC.
- Three independent reproductions of shipped numbers: whole-universe best 1.313, shipped
  distogram argmin 3.454, shipped synthesis 3.2041.
- No benchmark60 file, no `s9/final_cache`, no dev24 target, no `esm_cache.npz` was opened.
  No file outside `s12/` was modified. No git write command was run.

## Files

| file | what |
|---|---|
| `s12/asm_lib.py` | sub-length banks, fast many-vs-many Kabsch, batched CA builder |
| `s12/asm_rigid.py` | §1–2 rigid oracle ladder + retrieval/random arms |
| `s12/asm_chain.py` | exact branch-and-bound 2-piece chain floor (validation) |
| `s12/asm_chain2.py` | §3b matched-budget chain-consistent ladder + capacity null + k=3 |
| `s12/asm_null.py` | §2b rigid ladder on the capacity-null bank |
| `s12/asm_overlap.py` | §3c overlapping-piece join |
| `s12/asm_scaling.py` | §3 candidate-count scaling curves |
| `s12/asm_cover.py`, `s12/asm_cover_ss.py` | §4 conformational coverage + SS crossing |
| `s12/asm_deploy.py` | §5 deployable assembly |
| `s12/asm_report.py` | pooled tables and paired statistics |
| `s12/results/asm_*.json`, `s12/results/asm_*.log` | results |
