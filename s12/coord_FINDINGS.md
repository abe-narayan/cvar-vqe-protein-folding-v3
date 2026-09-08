# Sprint 12 — COORDINATOR findings

The coordinator's own work: instrument construction, integrity audits, and the
architecture-level null control. Agent findings are in `s12/{fail,asm,obj,agg,lit,forensics}_FINDINGS.md`.

---

## C0. The instrument (`s12/instrument.py`)

One module over the 126 cached window universes (`s8/generate_univ/*.npz`), so every agent
scores against the same objects and the benchmark is never touched. `python -m s12.instrument`
reproduces, exactly:

| quantity | pinned | reproduced |
|---|---|---|
| shipped distogram argmin | 3.4540 | 3.4540004952559396 |
| K=500 pool best (oracle) | 1.7108 | 1.7108244199364904 |
| top-75 best (oracle) | 2.3062 | 2.3061526409453816 |
| synthesis, λ=0 (`fit_ca`) | 3.2041 | 3.2040761603809194 |
| zero-recall targets | 18 | 18, identical set |

The 126 leave-fold-out distograms are cached to `s12/cache/disto_*.npz`, and a compact
per-residue ESM-2 bank (PCA-32, a fresh PCA-128, and contact maps for all 22,795 banked
sequences) to `s12/cache/esm_bank.npz` — 100 MB against the 1.5 GB original, so four agents
can hold ESM features at once on a 16 GB box. Coverage: 787/787 peptides, 6,003/6,003 small
fragments, 22,007/22,007 large fragments.

---

## C1. The TRAINING-SET containment leak is real but small — and 4 targets carry a verbatim copy

`s12/coord_contain.py` → `s12/results/coord_containment.json`.

Hazard H3 says `peptide_db.identity` normalises by the LONGER sequence, so a target sitting
verbatim inside a longer library member scores ~0.52, lands in a different identity cluster,
and can therefore land in a different FOLD — which puts it in the training set of the
distogram that scores the target. That pathway had never been measured. It is real:

| | count |
|---|---|
| tuning targets with a ≥0.6-**containment** member in their own distogram's training set | 74/126 |
| tuning targets with a **VERBATIM** copy in their own distogram's training set | **4/126** |

The four: `1CEK ⊂ 1A11`, `2FBU ⊂ 2LMF`, `2P5H ⊂ 2P5J`, `6B9K ⊂ 1U6V` — every one at
containment 1.00 and identity 0.52–0.59, i.e. every one passing the production 0.6 filter.

**It does not explain performance.** ρ(max training containment, emitted RMSD) = −0.129
(p = 0.15); the five targets above 0.9 containment average 3.115 Å against 3.407 for the
rest, on n = 5. And the direction is *against* the failure story: the FAIL18 have **lower**
training containment than the rest (0.301 vs 0.417), so the catastrophic targets are the
*clean* ones, not the contaminated ones.

**Status:** a genuine integrity defect worth fixing (write a containment-based screen and a
folds digest into every checkpoint), not an explanation of any result.

## C2. The 0.6 containment threshold is AT THE NULL for peptides — any audit that used it is uninterpretable

`s12/results/coord_containment_null.json`. Twelve composition-matched RANDOM sequences per
length, scored against the 6,003-fragment library and the 787-peptide database:

| n | max containment vs fragments (mean / max) | max identity vs fragments | max containment vs peptides |
|---|---|---|---|
| 9 | **0.62** / 0.78 | 0.43 | 0.55 |
| 11 | **0.58** / 0.67 | 0.44 | 0.54 |
| 13 | **0.56** / 0.60 | 0.45 | 0.52 |
| 16 | **0.63** / 0.80 | 0.41 | 0.45 |

A sequence with no relationship to anything scores 0.56–0.63 containment against the
fragment bank. **`containment ≥ 0.6` against a bank of this size therefore rejects a random
sequence about half the time and is evidence of nothing.** Identity is usable (null 0.41–0.45
against a 0.6 threshold); verbatim substring is unambiguous (p ≈ 0 under the null at n ≥ 9).

This is the sixth instance of the project's recurring failure mode — a quantity treated as
informative without its null — and it applies to my own first pass at C3, which I withdrew.

## C3. There is NO adequately-powered fresh benchmark. The supply is 16, and 10 of the 16 are amyloid fibrils

`s12/coord_benchN.py` → `s12/results/coord_benchN_manifest.json`.

Starting point: the three instruments between them spend **all 204** identity clusters of
9–16-residue peptides in the frozen corpus. Zero cluster-disjoint supply remains inside it.
So the question is what exists outside it. Measured against RCSB on 2026-09-05:

| supply | count |
|---|---|
| 9–16mer single-protein-entity entries, all time | 700 (472 solution NMR) |
| … released after 2026-01-01 (**date-fresh**) | 28 |
| … released after 2026-06-01 | 20 |
| … single-entity, no nucleic acid, absent from `pdbs/`+`pdbs_ext/` (**containment-fresh**) | 146 |

The date-fresh supply cannot reach the pre-registered minimum of ~40 *before* any quality
gate. The containment-fresh supply was built out in full — downloaded to `s12/newpdbs/`
(never to `pdbs_ext/`, which would re-permute `_scan` and re-draw the whole manifest):

| stage | n |
|---|---|
| downloaded | 133 |
| pass the structural gates (length 9–16, finite torsions, CA step 3.5–4.1, rebuild ≤ 1.5) | 40 |
| unique sequences | 30 |
| clean of leakage (no verbatim copy, identity < 0.6, containment < 0.85 — all clear of the C2 null) | 19 |
| one per identity cluster | **16** |

**16 targets, against a pre-registered minimum of 40.** And the composition is worse than the
count: **10 of the 16 are amyloid fibril / nanofibril segments**, 4 are X-ray designed
oligomers or bicyclics, and **2 are ordinary solution-NMR peptides**. The 9–16mer supply this
corpus never fetched is dominated by exactly the class `results/monomer_manifest.json` already
excludes as not folded in isolation.

**Consequence for the sprint, stated plainly: no architecture developed here can be validated
on a fresh, adequately-powered, distribution-matched benchmark, because the world does not
contain one.** The remaining options are (i) `dev24`, one pre-registered pass, (ii) a
length-extrapolation benchmark drawn from the 266 unused clusters of 8-mers and 17–26-mers
(powered — 382 peptides — but a different target distribution), or (iii) nested
cross-validation on tuning126 with the decision count declared. Any claim of a validated
improvement must name which of these it rests on.

---

## C4. Ideal-geometry template augmentation — NULL (`s12/coord_templates.py`)

The literature agent measured that on 2BFI a chain built from a single constant (φ,ψ) scores
0.27 Å while the best of all 21,547 windows in that target's universe scores 1.35 Å. If
canonical conformations are systematically absent from an empirical window bank, appending a
few hundred ideal templates to the pool is free and needs no detector. Tested: 307 templates
per target (10 textbook conformations, a 20°-step Ramachandran sweep including the
left-handed basin, and ideal β-hairpins over every interior turn position and four turn
types), appended to the K=500 pool, everything else untouched.

| ORACLE coverage | mean |
|---|---|
| best window in the whole legal universe | 1.313 |
| **best ideal template** | **2.496** |
| templates beat the whole universe on | **5/126** |
| templates beat the K=500 pool on | 18/126 |
| FAIL18: template best vs universe best | 3.336 vs 1.640 |

| emitted (λ=0.3) | mean | paired vs base | W/L | drop-10 |
|---|---|---|---|---|
| base (shipped pool) | 3.2126 | — | — | — |
| **+ 307 ideal templates** | **3.2169** | **+0.0042 [−0.048, +0.060]** | 53/51 | +0.050 |
| + 307 random universe windows (null) | 3.2365 | +0.0239 [−0.008, +0.057] | 65/61 | +0.052 |

**Null, and the mechanism is informative.** The filter is not rejecting the templates — 104
of 126 targets keep at least one in the top-75, 8.4 on average, and on 8 targets a template
is the score's outright argmin. They simply are not better: in aggregate the ideal
conformations are 1.2 Å *worse* than the empirical library, and on the FAIL18 they are 1.7 Å
worse. The literature agent's steric-zipper result is real but is a property of ~6 amyloid
targets, not a general coverage gap, and capturing it requires a class detector with the
false-positive risk that entails (+4.7 Å on one mis-routed target). **Generic canonical
conformations are not a missing channel.** This also independently corroborates the assembly
agent's finding that the library's holes are in *irregular* local structure, not in the
regular conformations an ideal template can express.

## C5. The sequence→structure channel, measured model-free (`s12/coord_infochannel.py`)

Every same-length pair in each corpus, sequence identity against CA-RMSD. No model, no
targets' natives, no training.

| corpus | n=9 | n=11 | n=13 | n=16 | mean pair RMSD |
|---|---|---|---|---|---|
| **peptide database** (787 chains) | −0.285 | −0.187 | −0.231 | −0.242 | 4.1–5.5 Å |
| **protein fragments** (6,003) | −0.036 | −0.030 | −0.031 | −0.030 | 3.6–5.5 Å |

(ρ(identity, RMSD); negative = more similar sequence gives more similar structure. Every
value p < 1e-5 — the fragment channel is significant and *tiny*, not absent.)

**The peptide corpus carries ~7× the sequence→structure correlation of the protein-fragment
corpus.** The mechanism is not mysterious and `fragment_db`'s own docstring states it: a
12-residue stretch cut out of a folded protein has the conformation its tertiary context
imposes, and its own sequence did not choose it. The record already knew the *training*-side
version of this (S7-2: adding protein-fragment pairs made the distogram monotonically worse).
The retrieval-side version had not been asked, and fragments are ~80% of every pool.

Inside a target's own universe, the deployable channel measures:

| quantity | value |
|---|---|
| ρ(BLOSUM similarity, true CA-RMSD) over the whole universe | **−0.066** (median −0.072) |
| … wrong sign on | 21/126 targets |
| … on the FAIL18 | **+0.0013** — nothing at all |
| … on the other 108 | −0.077 |
| BLOSUM top-500 vs a random 500: pool best | 1.711 vs 1.815 (−0.104 [−0.181, −0.033]) |
| BLOSUM top-500 vs a random 500: pool mean | 4.453 vs 4.817 (−0.364 [−0.427, −0.299]) |

So BLOSUM retrieval is real (S7-12 stands: it beats random) and is operating at a rank
correlation of about 0.07 — exactly the mixture the two corpus numbers predict for a pool
that is 80% fragments. **On the catastrophic targets the retrieval key carries literally zero
information about nativeness.**

## C6. A hypothesis of mine, refuted by the objective agent — recorded because it was wrong

From C5 and the objective agent's r_sep = 0.196 I inferred that the distogram's
sequence-separation component was dead weight — every pool window is a real fragment and
already satisfies the separation prior — and commissioned a deployable arm that scored
candidates on the separation-partialled deviation only. **It emits 3.642 Å, +0.436
[+0.276, +0.627] worse than shipped and worse than a random score (3.601).** The premise is
false: pool windows span helix to extended, so their per-shell distance profiles differ
enormously and the profile *is* what discriminates them. The agent's decomposition makes it
quantitative — the separation profile carries 75% of the objective's damage (0.607 Å of
0.813) and the within-shell residual only 25% (0.206 Å) — which is the opposite of what I
predicted. A low within-shell correlation does not mean the separation channel is
uninformative; it means the *within-shell* channel is nearly empty while the *profile*
channel is both informative and badly estimated.

## C7. Retrieval weighted toward the peptide corpus — NULL (`s12/coord_corpus.py`)

C5 says the peptide corpus carries ~7× the sequence→structure correlation of the protein
fragments that make up ~80% of every pool, and the forensics agent found the shipped top-75
*under*-selects peptide windows (0.177) exactly where the near-native band is *enriched* in
them (0.389). The historical `gen_pep` arm measured a peptide-only pool as null — but on the
ARGMIN terminal, which responds to the pool best, while the same file records the peptide
pool's better MEAN (4.321 vs 4.453), and the shipped operator responds to the mean. So the
arm had been scored against the terminal least able to use what it improved. Re-measured on
the synthesis arm, with a selectivity-matched control:

| arm | emitted | argmin | pool best | pool mean | top-75 best | peptide fraction of pool |
|---|---|---|---|---|---|---|
| base (shipped) | 3.2126 | 3.4540 | 1.7108 | 4.4533 | 2.3062 | 0.273 |
| **peptide-only** | **3.1922** | 3.4451 | 1.6714 | 4.3208 | 2.1730 | 1.000 |
| fragment-only | 3.3042 | 3.5617 | 1.8571 | 4.6055 | 2.3796 | 0.000 |
| selectivity control | 3.2769 | 3.5165 | 1.7153 | 4.6201 | 2.2694 | 0.215 |
| provenance boost ×10 | 3.1835 | 3.4515 | 1.6336 | 4.3382 | 2.2144 | 0.758 |

| paired vs base | value | W/L | drop-10 |
|---|---|---|---|
| peptide-only | **−0.0204 [−0.0970, +0.0519]** | 65/61 | +0.063 |
| fragment-only | **+0.0916 [+0.0453, +0.1398]** | 54/72 | +0.127 |
| boost ×10 | −0.0292 [−0.0734, +0.0158] | 62/64 | +0.024 |
| **peptide-only vs the selectivity control** | **−0.0847 [−0.1623, −0.0117]** | 70/56 | — |

**The direction is real and the size is nil.** Provenance beats a selectivity-matched random
sub-universe by −0.085 Å with a CI excluding zero, and dropping to fragments alone costs
+0.092 Å — so the corpus asymmetry measured in C5 does propagate into the pipeline. But
against the *actual* baseline, which already contains 27% peptide windows, peptide-only is
−0.020 Å with a CI spanning zero and a drop-top-10 that reverses the sign. Not worth a
config change. Recorded because the mechanism is sound and the null is informative: the
corpus asymmetry is real, and it is already largely captured.

## C8. The architecture's total sequence information is 0.776 Å — and it is NEGATIVE on the failures (`s12/coord_null.py`)

The 2×2 that had never been run. The pipeline has exactly two sequence-conditioned stages:
BLOSUM retrieval (which 500 windows enter the pool) and the distogram filter (which 75
survive). Everything downstream is sequence-blind geometry. Random arms averaged over 3
seeds; identical downstream path in every cell.

| arm | retrieval | filter | emitted | median | <2 Å |
|---|---|---|---|---|---|
| shipped argmin (reference) | BLOSUM-500 | distogram | 3.4540 | 3.478 | 0.214 |
| **shipped** | BLOSUM-500 | distogram-75 | **3.2126** | 2.966 | 0.286 |
| no filter | BLOSUM-500 | random-75 | 3.5974 | 3.515 | 0.222 |
| no retrieval | random-500 | distogram-75 | 3.3583 | 3.245 | 0.262 |
| **blind** | random-500 | random-75 | **3.9887** | 3.936 | 0.159 |

| contrast | mean diff | W/L | drop-10 |
|---|---|---|---|
| **both channels** (blind → shipped) | **−0.7760 [−1.0200, −0.5429]** | 90/36 | −0.522 |
| distogram alone (blind → no-retrieval) | −0.6304 [−0.8293, −0.4314] | 87/39 | −0.425 |
| BLOSUM alone (blind → no-filter) | −0.3912 [−0.4973, −0.2832] | 92/34 | −0.282 |
| distogram *given* BLOSUM | −0.3848 [−0.5632, −0.2186] | 84/42 | −0.191 |
| BLOSUM *given* the distogram | −0.1457 [−0.2262, −0.0722] | 81/45 | −0.047 |

**Everything the architecture knows about the target sequence is worth 0.776 Å**, and the two
channels overlap heavily (0.391 + 0.630 = 1.021 delivering 0.776 jointly) — consistent with
the distogram being largely typicality plus a weak sequence term. Both contrasts survive
dropping their ten best targets, so this is a broad effect, not a concentrated one.

The subgroup split is the sharpest single result of the sprint:

| subgroup | blind | shipped | sequence conditioning is worth |
|---|---|---|---|
| other 108 | 3.749 | **2.745** | **−1.004 Å** |
| **FAIL18** | **5.425** | **6.019** | **+0.594 Å — actively harmful** |

**A pipeline that is told nothing about the target sequence beats the real one on the 18
catastrophic targets.** The retrieval-key agent reached the same conclusion from the other
direction (on FAIL18 a content-free capacity null outperforms the trained predictor, and
BLOSUM is worse than noise). Three independently built routers are null, so the two regimes
cannot be told apart at inference time. This is the quantitative core of the sprint's
conclusion: the system's sequence conditioning is a genuine 1 Å channel where the target is
in-distribution and an active liability where it is not.
