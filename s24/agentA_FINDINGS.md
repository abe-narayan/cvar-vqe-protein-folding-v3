# SPRINT 24 — WORKSTREAM A — DATA, CORPUS AND LEAKAGE

**CORPUS_HASH = `29e3b67e8ca0c03d`** — cite this on every downstream artefact.

Pre-registration: `s24/PREREG_A.md`, sent to and approved by the coordinator before either run.
Code: `s24/a_corpus.py` (Run 1), `s24/a_ha2.py` (H-A2), `s24/a_charac.py` (Run 2).
Artefacts: `s24/results/a_corpus_permitted_29e3b67e8ca0c03d.json` (the permitted corpus + exclusion
list), `a_corpus_audit.json` (full audit, every hit record, per-target impact),
`a_ha2_fragment_leak.json`, `a_ha2_selfwindow.json`, `a_charac_cheap_29e3b67e8ca0c03d.json`,
`a_charac_heavy_29e3b67e8ca0c03d.json`.

**Benchmark discipline observed throughout.** Benchmark target IDENTIFIERS and SEQUENCES were read,
for exclusion only. No benchmark result file, no benchmark structure, and no benchmark RMSD was
opened at any point. Where sizing a finding would have required that, I stopped and handed the
decision up (§H-A2.5); the coordinator ruled the benchmark stays sealed.

---

## 0. HEADLINE, INCLUDING THE PARTS THAT DAMAGE MY OWN EXPECTATIONS

1. **The exclusion list exists and is materialised.** 443 of 6,790 source chains excluded (388
   deposits), 6,347 permitted, hash `29e3b67e8ca0c03d`, `complete=true` on the full key set.
2. **My own H-A2 was FALSIFIED.** I predicted the leak would enter through the fragment half. It
   does not; the fragment filter is sound in practice. It enters through the **peptide** half, by a
   normalisation defect the codebase already documents.
3. **A real leak reaches the sealed benchmark: 2 of 60 targets.** 8Y3S and 8ZG3 each sit verbatim
   inside 8XTO, which is in a different fold and therefore in both of their distogram training sets
   and both of their retrieval libraries. Sized on the dev instrument where sizing is permitted: the
   self-copy is BLOSUM rank #1 with certainty, but is **worse than the pool's own best on 3 of 4
   dev cases**. Real, structural, small.
4. **My H-A3 was also FALSIFIED, and in the opposite direction to the one I expected.** I predicted
   the corpus would prove to be a small number of folds repeated. At τ = 1.0 Å it is not: the
   near-duplicate leader curve is still climbing near-linearly at the largest sample I can afford,
   and at L=16 over the *whole* permitted corpus, 4,485 of 7,329 windows are mutually >1.0 Å apart.
   **The corpus is not repetitive — it is SMALL.** That is a different diagnosis and it points
   somewhere else. *(The one place I was right: at τ = 2.0 Å — the resolution the endpoint actually
   lives on — the 9-residue band does collapse and does saturate, to ~250 distinct shapes out of
   38,316 windows. It is a short-window phenomenon that weakens with length.)*
5. **The answer to Lane C's decisive question is NO, by my own pre-registered criterion**, and the
   binding number is not redundancy — it is that the permitted **peptide** bank is 410 chains and
   8,354 residues, and the 16-residue band contains 7,329 windows in total, below my pre-registered
   10^4 floor before any independence collapse is applied. §D2.6 states this precisely, with the
   one caveat that would make it a qualified NO rather than a flat one.

---

## D1. THE CORPUS, THE EXCLUSION LIST, AND THE HASH

### D1.1 What the corpus actually is — read from the builder, not assumed

`s8/generate.py:stage_univ` writes one npz per tuning target. For target `p` in fold `f`:

```python
peps  = [q for q in peptide_db.load() if folds[q.seq] != f and q.seq != p.seq]
frags = distogram._fold_fragments(f, 5, threshold=0.6)
W     = every length-n window of (peps + frags), IN THAT ORDER
```

`W`, `PHI`, `PSI`, `S`, `org`, `sim`, `order`, `rr`, `nat_ca` are all derived from that. There is no
third source. The universe is therefore a deterministic function of exactly two banks:

| bank | module | cache | chains | deposits | lengths | origin |
|---|---|---|---|---|---|---|
| peptides | `peptide_db` | `peptide_db.npz` | 787 | 787 | 8–26 | `pdbs/` (61) + `pdbs_ext/` (1,463), RCSB single protein chains 8–26 aa, queried 2026-08-31, gated on torsion-rebuild fidelity ≤1.5 Å and exact-sequence dedup |
| fragments | `fragment_db` | `fragment_db.npz` | 6,003 | 1,001 | 9–20 | contiguous windows of `prots/` (13,751 files scanned), stride 5, ≤6 per protein, rebuild tol 1.0 Å |

**Deposit overlap between the two banks: ZERO.** 6,790 source chains total.

**Provenance is proven, not inferred.** `s24/a_corpus.py:universe_provenance` reconstructs the exact
`(source chain → window index range)` map for every target from these two banks and asserts it
against the npz files. **Window counts match exactly on 126/126** — 2,352,893 windows, mean 18,674
per target, range 7,016 (1RG4) to 39,410 (7N2I). Any lane can now map any window index
back to its source chain.

### D1.2 The instruments, and why containment was not one of them

Pre-registered functional fork. Leak = union of {exact sequence equality, **verbatim substring in
either direction**, Needleman–Wunsch **identity ≥ 0.6 normalised by the longer sequence**,
source-deposit-id equality}.

**Containment ≥ 0.6 was NOT used, and this is deliberate — say so to anyone who later proposes to
"fix" the leak with it.** Project memory (`containment-threshold-is-at-the-null`) measures
composition-matched random 9–16mers at **0.56–0.63 containment** against this exact bank. A
containment-0.6 rule is a threshold at chance and is evidence of nothing.

Both instruments I did use were re-nulled here, on this corpus, per the pre-registered null fork.
12 composition-matched random sequences per length 9–16, 96 total, seeded from
`s15.seed.stable_rng`, scored against the identical bank by the identical criteria:

| n | max NW identity (mean) | max NW identity (max) | verbatim-substring hits | exact hits |
|---|---|---|---|---|
| 9 | 0.491 | 0.556 | 0/12 | 0/12 |
| 10 | 0.485 | 0.545 | 0/12 | 0/12 |
| 11 | 0.492 | 0.571 | 0/12 | 0/12 |
| 12 | 0.463 | 0.533 | 0/12 | 0/12 |
| 13 | 0.461 | 0.500 | 0/12 | 0/12 |
| 14 | 0.475 | 0.571 | 0/12 | 0/12 |
| 15 | 0.444 | **0.600** | 0/12 | 0/12 |
| 16 | 0.431 | 0.500 | 0/12 | 0/12 |

**identity ≥ 0.6 sits at roughly the 1% tail of the null** (1 of 96 random sequences touched it, and
only at n=15). **Verbatim substring is at p ≈ 0** — 0 of 96, at every length. Both instruments are
clear of chance. Neither number is quoted anywhere in this document without this table behind it.

### D1.3 The exclusion list

Pre-registered basis fork: the verdict is taken at **source-chain level and propagated to every
window derived from it**; for a fragment the source chain is its parent deposit, so one contaminated
fragment removes all six fragments cut from that deposit. That is the strict reading and it is what
shipped. Both the strict and the loose count are in the artefact so the operator is visible.

```
corpus source chains                        6,790
chains HIT by at least one criterion           389
chains EXCLUDED after deposit propagation      443     <- the shipped exclusion list
deposits excluded                              388
chains PERMITTED                             6,347     <- the shipped training corpus
```

Per-criterion, published for every criterion including the ones that fire least:

| held-out set | exact | verbatim substring | identity ≥0.6 | same deposit |
|---|---|---|---|---|
| tuning126 | 126 | 148 | 202 | 126 |
| dev24 | 24 | 32 | 47 | 24 |
| bench60 | 60 | 76 | 129 | 60 |

Every target is itself a `peptide_db` entry, so "exact" necessarily equals the set size — that is a
sanity check on the audit, not a finding.

**Impact on the shipped universes:** 107,761 of 2,352,893 windows (4.58%) come from an excluded
source. Per target: mean 4.35%, median 4.43%, **worst 5.19%**. 82.0% of the excluded windows are
peptide-derived. So a generator trained on the permitted corpus loses roughly one window in 23, and
no target loses more than one in 19.

### D1.4 A hazard found while building the key: fragment ids are not unique

`fragment_db` ids are `<PDB>_<start>`, but one start yields one fragment per length in
`LENGTHS = 9..20`, so **803 of the 6,003 ids collide**. Any downstream code that keys fragments by
`f.pdb` silently merges 9-mers with 20-mers. The corpus primary key used here is `kind:id:len` and is
asserted injective. **Lanes B and C must key on that, not on `f.pdb`.**

---

## H-A2. THE LEAKAGE MECHANISM — MY HYPOTHESIS WAS WRONG AND THE TRUE ONE REACHES THE BENCHMARK

### H-A2.1 What I predicted, and why it is false

I predicted that `distogram._fold_fragments` filters fragments only against the *held-out fold's*
peptides while the 60 benchmark targets sit in *other* folds, so the fragment filter never sees them.
**Read from source, that is false.** `core/bench.py:207` and `core/pipeline.py:1315` both set
`fold = db.folds(cfg.n_folds)[target.seq]`, so a target's **own** fold is the held-out fold at
inference, and the fragment filter *is* run against the target's own sequence. H-A2 as written is
falsified. I record it because a falsified hypothesis of my own is the cheapest thing in this report
and because the true mechanism was only found by chasing it.

### H-A2.2 The fragment half is clean — measured, not argued

Across all 210 targets (126 + 24 + 60), against the fragments that actually survive into each
target's own fold model:

| test | tuning126 | dev24 | bench60 |
|---|---|---|---|
| fragment contains the target verbatim, or vice versa | 0/126 | 0/24 | 0/60 |
| fragment ≥0.6 identical that the **unsound 3-mer prefilter** would have skipped | 0/126 | 0/24 | 0/60 |
| fragment's parent deposit id equals the target's deposit id | 0/126 | 0/24 | 0/60 |

The filter removes only 0–15 of 6,003 fragments per fold (fold sizes 6001/5988/6001/6003/6002) — and
the raw bank does contain 12 related fragments, so it is catching real things. Spot-checked directly:
`1EY9_91` (`MVNEALVRQGLA`, identity 0.923 and a verbatim substring of dev target 2FXZ `KMVNEALVRQGLA`)
**is absent** from the shipped fold-1 fragment set. **The fragment filter works.**

**Known defect with measured zero impact.** `distogram._fold_fragments` still carries the 3-mer
prefilter (`if fk and k and not (fk & k): continue`) that `peptide_db.clusters` removed from the
clustering as *unsound* — its premise, "a pair above 0.6 identity must share a 3-mer", is false, and
it once let a benchmark target train on a 0.70-identity homolog. It survives in the fragment filter.
Measured impact today: **exactly zero** missed homologs across all 210 targets. Coordinator ruled:
logged as a known defect, **not patched this sprint**.

### H-A2.3 The peptide half is where the leak is, and the codebase already says so

`peptide_db.identity` normalises by the **longer** sequence. So an 11-mer target sitting **verbatim**
inside a 21-mer library chain scores 11/21 = 0.52, falls below the 0.60 threshold, lands in a
**different identity cluster**, therefore a **different fold**, therefore is in that target's fold
model's training set *and* in `library_members`, which excludes only exact-sequence matches.

This is not a subtle discovery. `core/data.py:188` documents the convention as leaky at the member
level in its own docstring. And `s9/final.py:210` states the opposite as a guarantee —

> "anything at or above `db.IDENTITY_THRESHOLD` to the target is in the target's own cluster, hence
> its own fold, hence excluded here"

— which is exactly what the longer-normalisation breaks. **The consequence was documented and never
counted.** The count is this workstream's contribution.

| set | targets carrying a verbatim self-copy in their own fold model's training set |
|---|---|
| tuning126 | **4 / 126** — 1CEK in 1A11, 2FBU in 2LMF, 2P5H in 2P5J, 6B9K in 1U6V (reproduces the standing project fact exactly) |
| dev24 | **0 / 24** |
| bench60 | **2 / 60** — see below |

### H-A2.4 The two benchmark targets, from identifiers and sequences only

```
8Y3S  11 aa  GVAFRAPSIHG              fold 1
8ZG3  11 aa  LGGGSVRFGPG              fold 4
8XTO  21 aa  LGGGSVRFGPGVAFRAPSIHG    fold 3
8ZUG  23 aa  YRQSSATSSFGGLGGGSVRFGPG  fold 1
```

`8XTO[0:11] == 8ZG3` and `8XTO[10:21] == 8Y3S`. **8XTO is literally the chain spanning both
benchmark targets, and it sits in fold 3 — i.e. in the training set and the retrieval library of
BOTH of them.** 8ZUG additionally contains 8ZG3 verbatim (8ZUG is fold 1, so it is correctly
excluded from 8Y3S's model, but not from 8ZG3's).

### H-A2.5 The size, measured only where measuring is permitted

On the 4 affected **dev** targets I reconstructed the exact window index of the self-copy in the
shipped universe and read its BLOSUM rank and its ORACLE `rr`:

| target | n | BLOSUM rank of the self-window | its `rr` (Å) | pool best K=500 | universe best |
|---|---|---|---|---|---|
| 1CEK | 13 | **0 (rank #1)** | 0.595 | 0.342 | 0.310 |
| 2FBU | 12 | **0** | 3.278 | 2.243 | 1.154 |
| 2P5H | 9 | **0** | 2.334 | 1.840 | 1.386 |
| 6B9K | 10 | **0** | 4.126 | 2.077 | 1.515 |

**The exposure is certain** — a verbatim self-sequence window is BLOSUM rank #1 and is therefore in
the K=500 pool with probability 1. **It is not a near-native answer**: on 3 of the 4 it is *worse*
than the pool's own best, by 1.0–2.0 Å. The same sequence in a different deposit adopts a different
conformation. This is consistent with the standing project note that this does not explain
performance.

**Where I stopped.** Whether 8Y3S/8ZG3's self-windows are near-native, and therefore whether the
2/60 moved any reported number, requires reading benchmark structures' RMSD. **I did not do it and I
did not estimate it.** The coordinator ruled: the benchmark stays sealed, the 2/60 is carried as a
declared caveat, nobody quantifies it. Re-pinning folds with a substring test is logged as a
standing decision for the user, not taken, because `peptide_db.folds` warns that renumbering
silently reassigns every sequence and invalidates every trained model on disk.

**Read on the incumbent:** 4/126 affected, self-window worse than pool best on 3 of 4. The 3.0483 Å
incumbent is **declared, not retracted**.

---

## D2. CORPUS CHARACTERISATION

All numbers on the **permitted** corpus at `CORPUS_HASH 29e3b67e8ca0c03d`.

### D2.1 Size, and the number that actually binds

```
permitted source chains          6,347
permitted independent deposits   1,400
  peptides                         410 chains / 410 deposits /  8,354 residues
  fragments                      5,937 chains / 990 deposits / 80,738 residues (exactly 6 per deposit)
```

**Nearly half the peptide bank is gone.** 410 of 787 peptides survive exclusion — the 210 held-out
targets are themselves peptides, and their ≥0.6-identity neighbourhoods take the rest.

**This is the finding Lane C most needs and it is easy to miss.** The generator's *training* corpus
is strictly smaller than the retrieval library it must beat. Retrieval draws from ~630 out-of-fold
peptides per target; a legitimately-trained generator gets **410**. And project memory records that
the peptide corpus carries **7× the sequence–structure channel** of the protein fragments that make
up 80% of every pool. The generator is therefore handicapped precisely in the half that carries the
signal it needs.

Windows per length, permitted corpus:

| L | peptide windows | fragment windows | total | peptide share |
|---|---|---|---|---|
| 9 | 5,074 | 33,242 | 38,316 | 0.132 |
| 10 | 4,667 | 27,305 | 31,972 | 0.146 |
| 11 | 4,260 | 22,051 | 26,311 | 0.162 |
| 12 | 3,854 | 17,456 | 21,310 | 0.181 |
| 13 | 3,449 | 13,507 | 16,956 | 0.203 |
| 14 | 3,047 | 10,174 | 13,221 | 0.230 |
| 15 | 2,649 | 7,378 | 10,027 | 0.264 |
| 16 | 2,258 | 5,071 | **7,329** | 0.308 |

**The 16-residue band has 7,329 windows in total, before any independence collapse.** That is already
below the 10^4 floor I pre-registered as a NO condition.

### D2.2 Vocabulary, torsions and multimodality — the label's second half PASSES

- 89,092 permitted residues; full 20-letter vocabulary, no unknown characters; vocabulary entropy
  4.19 bits of a possible 4.32.
- Global Ramachandran basin composition: **αR 0.536, β/PPII 0.394, αL 0.055, other 0.016.**
- **Residue types whose φ/ψ is >80% concentrated in one basin: 0 of 20.** Every residue type is
  genuinely bimodal between αR and β. The pre-registered ">80% in one basin per residue type" NO
  condition **does not fire**.
- The φ channel is near-constant across sequence, consistent with the standing project result that
  φ carries no sequence signal at peptide length; the variance lives in ψ.

### D2.3 Secondary structure — and a real structural bias

Simplified DSSP over the permitted chains (exact per chain, then sliced to windows):

```
H 0.375     E 0.022     C 0.604
```

**β-sheet is essentially absent (2.2%).** That is mechanical, not accidental: DSSP assigns E only
with a *paired* strand, and both banks are isolated short chains — a fragment's partner strand was
left outside the window by construction. `fragment_db`'s own docstring states the caveat: a
fragment's conformation is partly held by contacts outside the window, so it is "a source of
geometry, not a source of what does this peptide do on its own", and fragments were declared
training-only for exactly this reason.

**Consequence for both generator designs: a model trained on this corpus will have seen almost no
extended-sheet conformations in their native pairing context, and will be biased toward helix and
coil.** Any evaluation that reports a mean over targets should also report whether the sheet-bearing
targets are the losses.

Distinct H/E/C window strings, permitted corpus:

| L | distinct SS strings | of 3^L possible | ESS exp(H) | most common string's share |
|---|---|---|---|---|
| 9 | 355 | 19,683 | 21.7 | 0.279 |
| 10 | 467 | 59,049 | 29.1 | 0.249 |
| 11 | 591 | 177,147 | 38.9 | 0.221 |
| 12 | 711 | 531,441 | 52.0 | 0.195 |
| 13 | 819 | 1.6e6 | 68.7 | 0.172 |
| 14 | 891 | 4.8e6 | 88.8 | 0.153 |
| 15 | 913 | 1.4e7 | 112.0 | 0.137 |
| 16 | 895 | 4.3e7 | 137.5 | 0.122 |

**At the level of secondary-structure pattern the corpus has an effective size of 20–140.** Note the
distinct-string count *plateaus* at ~900 from L=14 onward while the number of possible strings grows
by 3× per residue — the corpus stops producing new SS patterns and starts producing new *geometry
within* the same patterns.

### D2.4 The near-duplicate structure: THE DECISIVE MEASUREMENT

Greedy leader clustering under Kabsch CA-RMSD, length-stratified, on nested random subsamples;
`complete=true`. The Rg prefilter used to make it affordable is **admissible** (RMSD ≥ |Rg_a − Rg_b|
follows from `sum(S) ≤ sqrt(|a|²|b|²)`) and is self-tested against the unpruned computation. The
L=16 arm reproduced an earlier independent pilot to the digit, so the walk is deterministic.

**Leader counts at the largest sample available per band** (N = 10,000, except L=16 where the whole
permitted corpus is 7,329):

| L | n windows | τ=0.5 Å | τ=1.0 Å | τ=1.5 Å | τ=2.0 Å |
|---|---|---|---|---|---|
| 9 | 38,316 | 4,837 | 2,549 | 903 | **213** |
| 10 | 31,972 | 5,234 | 3,316 | 1,610 | 468 |
| 11 | 26,311 | 5,671 | 4,070 | 2,379 | 887 |
| 12 | 21,310 | 5,908 | 4,488 | 3,018 | 1,464 |
| 13 | 16,956 | 6,169 | 4,872 | 3,597 | 2,026 |
| 14 | 13,221 | 6,493 | 5,173 | 4,098 | 2,663 |
| 15 | 10,027 | 6,736 | 5,454 | 4,398 | 3,138 |
| 16 | **7,329** | 5,375 | **4,485** | 3,717 | 2,815 |

**Against the two pre-registered plausible zero-information controls, at MATCHED N = 5,000:**

| L | corpus τ=1.0 | Ramachandran-marginal control τ=1.0 | constant α-helix control (any τ, any L) |
|---|---|---|---|
| 9 | 1,648 | 3,312 | 1 |
| 10 | 2,038 | 4,330 | 1 |
| 11 | 2,430 | 4,837 | 1 |
| 12 | 2,658 | 4,960 | 1 |
| 13 | 2,818 | 4,993 | 1 |
| 14 | 2,975 | 4,999 | 1 |
| 15 | 3,105 | 5,000 | 1 |
| 16 | 3,271 | 5,000 | 1 |

The corpus sits roughly **2× more concentrated** than a plausible uncorrelated draw and vastly above
the degenerate floor. It carries real, non-trivial structure. **It does not carry enough of it.**

**H-A3 IS FALSIFIED, AND THE DIRECTION MATTERS.** I predicted the τ=1.0 Å independent count would be
one to two orders of magnitude below the raw window count. It is not. Doubling the sample from 5,000
to 10,000 multiplies the leader count by 1.55 (L=9) to 1.76 (L=15), and at L=16 over the *entire*
permitted corpus, 4,485 of 7,329 windows are mutually >1.0 Å apart. **There is no saturation. The
corpus is not repetitive — it is small.** At L=9–14 my leader count is bounded by my sample size, so
what I report there is a **lower bound on independence, not the value.**

**The one place my prediction was right, and it is the resolution that matters most.** At τ = 2.0 Å —
which is the scale the endpoint actually lives on, given a dev pool best of 2.355 Å and an incumbent
of 3.048 Å — the 9-residue band **does** collapse and **does** saturate: 186 leaders at N=5,000 and
213 at N=10,000, a ratio of 1.15 for a doubled sample. **At the resolution the project's metric cares
about, a 9-mer corpus of 38,316 windows contains on the order of 250 distinct shapes.** That collapse
weakens with length (L=16 at τ=2.0 is 2,815 of 7,329 and still climbing), so it is a short-window
phenomenon.

**Secondary, pre-registered, in torsion space** (leader clustering on `[cosφ,sinφ,cosψ,sinψ]` at
0.35 rad RMS): 4,345 / 4,604 / 4,829 / 5,074 / 5,257 / 5,431 / 5,670 leaders at N=10,000 for L=9–15,
and 3,354 at N=5,000 for L=16. **It agrees with the CA-space primary in shape** — near-linear growth,
no saturation — at uniformly higher counts, i.e. torsion space is the less forgiving metric. **No
disagreement to report**, which is the outcome I committed to publishing either way.

### D2.5 Representations already in the repo

| representation | where | shape / convention | note for a generator |
|---|---|---|---|
| CA trace | universe `W` `(nw, n, 3)` float32 | deposited coordinates, **not centred, not scaled** | every RMSD in the project is Kabsch, rotation+translation only — **scale is never fitted**, so a systematic scale error goes straight into the endpoint (s23 L1) |
| torsions | `PHI`, `PSI` `(nw, n)` **float16**, radians | parent's real φ/ψ; φ of residue 0 is undefined and fragment extraction skips it | float16 is ~0.001 rad — adequate for a target but **not** for a residual whose scale is smaller than that; Lane B should recompute φ/ψ in float64 from source rather than train on the npz arrays |
| sequence codes | `S` `(nw, n)` int8 over `ARNDCQEGHILKMFPSTWYV` | `s12.instrument.decode` | |
| ideal-geometry builder | `core/project.build_ca_exact(phi, psi)`, batched | `N-CA 1.458, CA-C 1.525, C-N 1.329; ∠N-CA-C 111.0°, ∠CA-C-N 116.2°, ∠C-N-CA 121.7°; ω = π (trans, fixed)` | **bit-identical** to `core/geometry.build_backbone_batch(...)["CA"]` and that identity is load-bearing — the fast scan builder differs at 1e-13 Å and moves the emitted structure on 126/126 targets, median 0.031 Å, worst 1.63 Å. Use `build_ca_exact`. ω is fixed trans, so **cis-peptide conformations are unreachable by construction**. |
| distogram | `s12.instrument.distogram(pdb)` → `prob (npairs,17)`, `expected`, `sd`, `risk (npairs,760)` on grid 2–40 Å step 0.05 | leave-fold-out MLP ensemble, `core/pipeline.fold_model(fold)` | a **trained predictor with its own training set** — see D3 |
| ESM | `esm_features.py`, `esm_cache.npz` (22,795 seqs), `esm_small.npz` (360), `esm_pca.npz` | ESM-2 `esm2_t33_650M_UR50D`; final-layer per-residue embedding reduced to **32 PCs**; plus **attention contact logits** | two separate exposure paths — see D3 |
| pair index | `s12.instrument.pair_index(n, min_sep=2)` | upper triangle, |i−j| ≥ 2 | |
| readout | `s12.instrument.coordinate_average` | superpose on the **medoid**, then mean, m=75 | contracts the virtual CA-CA bond by 22% while the radius of gyration is only 6% short (s23 L1) — a generator that emits *single* structures is not on the same manifold as the incumbent's output and must not be compared to it without stating the basis |

**Basis discipline reminder for both lanes:** point-cloud and built-chain RMSD differ by 0.156 Å of
pure operator choice. State the basis on both sides of every contrast.

### D2.6 THE DECISIVE ANSWER FOR LANE C

**NO — the permitted corpus does not contain enough genuinely independent structural information to
support a from-scratch generative model for 9–16mers, by the criterion I stated before I saw the
data. And the reason is not the one I expected.**

The pre-registered label had two halves and a set of NO conditions. Scoring them honestly:

| pre-registered NO condition | fires? | evidence |
|---|---|---|
| fewer than ~10^4 independent windows in ANY length band | **YES** | L=16 has **7,329 windows in total**, before any independence collapse; L=15 has 10,027 raw of which only 2,649 are peptide-derived |
| >80% of windows in one Ramachandran basin class per residue type | **no** | 0 of 20 residue types exceeds 0.80; the most concentrated is Ala at 0.713 |
| independent count does not clearly exceed the plausible zero-information control | **no** | corpus is ~2× more concentrated than the Ramachandran-marginal control and ~3,000× above the α-helix floor — it clearly carries structure |

One of three conditions fires, and it is the one about size. **The corpus is torsionally multimodal
and structurally non-redundant. It is simply too small in the length bands that matter.**

**The number that actually binds is worse than the one my label named.** The permitted corpus is
1,400 independent deposits, and the half that is *isolated-peptide* geometry — the thing a 9–16mer
generator is being asked to model — is **410 chains and 8,354 residues**. The other 91% of residues
are protein fragments that `fragment_db`'s own docstring declares training-only, because their
conformation is held by contacts outside the window. And retrieval, which the generator must beat,
draws from ~630 out-of-fold peptides per target while the generator gets 410. **The generator is
handicapped, relative to the incumbent, precisely in the half that carries the sequence–structure
signal** — and project memory prices that half at 7× the fragments' channel.

**What I recommend, holding no stake now that the answer is NO:**

1. **Lane B's residual design is not bound by this finding in the same way, and this is the important
   corollary.** A residual model learns a correction over ~2.35M (candidate, target) pairs; its
   sample count is the *pool* size, not the corpus size. The 7,329-window ceiling constrains
   *p(structure)*, not *p(correction | structure)*.
2. If a from-scratch model is built anyway, build the **smallest** one that shows genuine
   multimodality (as BRIEF §3 already says) and quote its parameter count against **7,329**, not
   against the 18,674-windows-per-target figure the corpus advertises.
3. Either way, **stratify by native secondary structure**. At E = 2.2% the corpus cannot teach a
   model what a paired strand looks like.

**Five caveats against my own conclusion, stated so this is not one-sided.**

1. The 10^4 threshold is mine. It was set in advance and rests on a naive ten-samples-per-parameter
   rule; a 10^3-parameter model would clear it. I said so in the pre-registration, before the data.
2. Windows of different lengths share residues, so the eight bands are **not** eight independent
   pieces of evidence.
3. A leader count is a covering number at one resolution. It is not information content, and a
   generative model does not have to memorise windows — it can generalise across them. "10× the
   parameter count" is a rule of thumb, not a theorem.
4. At L=9–14 the count is bounded by my sample size, so those are lower bounds; the true independent
   counts are larger and my NO is *conservative in the wrong direction* there. The NO rests on L=15
   and L=16, where the raw corpus is smaller than the threshold and no sampling argument can help.
5. The finding is about *this* corpus. Nothing here says a from-scratch generator is impossible —
   only that it is not supported by 410 permitted peptides. Project memory
   (`no-fresh-benchmark-exists`) records that the containment-fresh world supply of new 9–16mer
   targets is 16, ten of them amyloid fibrils, so enlarging the peptide bank is not a cheap fix
   either.

---

## D3. THE LEAKAGE-RISK REGISTER — WRITTEN BEFORE EITHER GENERATOR IS BUILT

### D3.0 The shared exposure both designs inherit: the distogram's own training set

**What it was.** `core/pipeline.fold_model(fold)` → `distogram.train_fold(fold, use_esm=True,
n_folds=5, fragments=True)`, which trains on

```
[p for p in peptide_db.load() if folds[p.seq] != fold]   +   _fold_fragments(fold, 5, 0.6)
```

i.e. **~630 out-of-fold peptides and ~6,000 fragments**, with the target's own identity cluster held
out.

**Does it overlap the dev targets? Yes, massively — and this is the shipped condition, not a new
defect.** The 126 dev targets are spread across the 5 folds (25/23/25/23/30), so:

```
other dev126 targets inside a dev target's own fold model's training set   mean 100.5  (min 96, max 103)
other bench60 targets inside a bench target's own fold model's training set mean  45.8
dev126 targets inside a benchmark target's fold model's training set        mean 101.2
```

The target itself and its cluster are correctly excluded. But **each fold model saw ~100 of the other
125 dev natives**, and there are only **five** distinct models covering all 126 targets.

**Three consequences, all of which bind this sprint:**

1. **A generator conditioned on the distogram inherits ~100 dev natives of exposure per target.** It
   is not conditioned on the target's own native — but the distogram is a *learned summary of the
   dev distribution*, and a generator that learns to trust it is learning the dev distribution
   through a legitimate-looking channel.
2. **Per-target results are not independent across the 126.** Five models, ~100 shared natives each.
   iid CIs over 126 targets will be anticonservative. **Fold-clustered CIs are not optional here** —
   they are the only honest interval, and BRIEF §4 already requires them.
3. **The 4/126 (and 2/60) verbatim self-copies of §H-A2 are exposure to the target's OWN native**,
   which is categorically different from the above and is the only part that is a defect.

### D3.1 ESM-2 — the exposure nobody in this project can audit

Two distinct channels, and they carry different risk:

| channel | what it is | risk |
|---|---|---|
| final-layer embedding → 32 PCs | ESM-2 650M, trained on UniRef50 — **sequence only** | The model has almost certainly seen sequences homologous to these peptides. Sequence-only pretraining is not structural leakage, but it is not auditable and cannot be excluded from. **Declare it; do not claim to have excluded it.** |
| **attention contact logits** | ESM-2's contact head is a **logistic regression trained on real PDB contacts** — `esm_features.py`'s own docstring says it "is the only part of ESM-2 that was ever supervised on structure" | **This channel is structure-supervised on a training set that certainly includes deposits related to these targets, and that training set is not available to audit.** It is the single largest un-auditable structural exposure in the pipeline. |
| the PCA itself | `esm_pca.npz` — 1280→32 projection "fitted once over the whole database" | The projection directions were chosen using held-out targets' **sequences** (not structures). Low severity, but it is a global fit across the split and should be declared, not defended. |

`esm_features.py` says "Nothing here reads a structure, so ESM features are available for a target
whose native is being held out." That is true **of inference** and false of **training**. Both
generator designs should state which of these three channels they use.

The counter-consideration, stated so this is not one-sided: memory
(`esm-adds-nothing-for-short-peptides`, OVERTURNED) records that on SELECTION, ESM buys 0.288 Å over
one-hot (p=0.005, n=126). Dropping ESM to close this exposure has a measured price.

### D3.2 Design (a) — the RESIDUAL generator: `retrieved candidate + learned (Δφ, Δψ)`

| # | path by which native information could enter | the check that catches it |
|---|---|---|
| a1 | **Training on the target's own native residual.** The residual's natural label is `native − retrieved`; if the training pair is ever formed for a held-out target, the model is trained on the answer. | Assert at dataset-construction time that no training pair's target chain appears in the exclusion list. Cite `CORPUS_HASH 29e3b67e8ca0c03d`. Hard assert, not a filter. |
| a2 | **`rr` and `nat_ca` in the universe npz.** Every universe carries ORACLE arrays. A conditioning tensor assembled by `{k: z[k] for k in z.files}` picks them up silently. | NaN-poison `rr` and `nat_ca` before building any conditioning tensor and assert the model output is bit-identical. `s8/test_generate.py` already does exactly this for the deployable/diagnostic split — reuse that test, do not rewrite it. |
| a3 | **Selecting the retrieved candidate to condition on using `rr`.** "Condition on the best retrieved candidate" is native selection wearing an innocent name. | The conditioning candidate must be chosen by a native-free rule (BLOSUM rank, shipped score) and that rule must be stated in the module docstring. |
| a4 | **A loss that consumes native RMSD, native contacts, native compactness or native scale.** BRIEF §4 forbids it. Note the specific trap: the incumbent's output is 22% contracted in bond length, so a loss with *any* scale term is fitting a known native-derived quantity. | Enumerate every loss term and state, per term, what it reads. Reject any term whose value changes when the native is permuted. |
| a5 | **Early stopping / hyperparameters chosen on dev RMSD.** Nested CV is required by BRIEF §4. | Nested CV over the pinned folds. Any hyperparameter fitted on the folds it is scored on is reported separately as optimism. |
| a6 | **Distogram conditioning** — inherits D3.0 in full (~100 dev natives per target, 5 models). | Report the distogram-conditioned and distogram-free arms. If the model only works with the distogram, say so. |
| a7 | **The float16 trap, which is a correctness risk not a leakage one but will be mistaken for one.** `PHI`/`PSI` are float16 (~0.001 rad). A residual smaller than that is quantisation noise, and a model that appears to learn a tiny residual may be learning the quantisation grid. | Recompute φ/ψ in float64 from the source chains. Report the residual's magnitude against 0.001 rad. |
| a8 | **The 4 dev / 2 bench self-copies of §H-A2** are the highest-BLOSUM-rank candidates for those targets, so a residual model conditioned on "the top-ranked retrieved candidate" conditions on a verbatim self-copy on exactly those targets. | Report those 4 dev targets separately. Do not let them into a headline mean without a footnote. |

**Mode-collapse audit precedes any RMSD claim** (BRIEF §4): unique-structure count, torsion entropy,
pairwise RMSD, effective sample size, duplicate fraction. **Geometric validity with every pool:**
Ramachandran, ω, bond length, bond angle, clashes, chain breaks.

### D3.3 Design (b) — the FROM-SCRATCH generator: `p(φ,ψ | sequence, distogram)`

Everything in D3.2 applies except a3 and a7, plus:

| # | path | the check |
|---|---|---|
| b1 | **The distogram is the dominant conditioning signal and it carries D3.0's exposure in full.** A from-scratch generator conditioned on it is not "from scratch" with respect to the dev distribution — it is a decoder for a predictor that saw ~100 dev natives. | State it in the module docstring. Run a **distogram-free** arm (sequence only). The gap between them is the size of the inherited exposure, and it should be reported as such rather than as an ablation win. |
| b2 | **Training on windows of the corpus while evaluating on peptide targets.** The corpus is 91% fragment residues, and fragments are context-held geometry that `fragment_db` itself declares training-only. | Report peptide-only and fragment-included training arms. A win that only exists with fragments is a win on the wrong distribution. |
| b3 | **The sheet deficit (§D2.3, E = 2.2%).** A generator that never emits extended-paired conformations will look good on a helix/coil-dominated panel and fail silently on sheet targets. | Stratify every result by the target's native SS class. |
| b4 | **Sampling temperature / number of samples tuned on dev RMSD.** Best-of-K nulls are the distribution of the MAXIMUM (BRIEF §4); a "−0.077 Å oracle" last sprint was 101% accounted for by its own best-of-K null. | Every best-of-K claim carries its size-matched best-of-K null. |
| b5 | **The `s12` targets list is derived from the universe directory.** `s12.instrument.targets()` globs `s8/generate_univ/*.npz`. Adding or removing a file silently changes the instrument. | Assert `len(targets()) == 126` and that the pdb set matches `debias.tuning_targets()` at the top of every run. |
| b6 | **ω is fixed trans in the builder.** A generator that learns φ/ψ from deposited structures containing cis-prolines learns torsions its own builder cannot realise. | Report the fraction of training residues whose deposited ω is non-trans, and exclude or flag them. |

### D3.4 What would falsify "the corpus is clean"

Stated now so it can be checked later rather than argued: re-run `s24/a_corpus.py`; the hash must
come back `29e3b67e8ca0c03d`. If `peptide_db.npz`, `fragment_db.npz`, `peptide_clusters.json` or
`peptide_folds.json` changes, the hash changes and **every downstream result that cited the old hash
is void**, because `peptide_db.folds` warns that renumbering clusters silently reassigns every
sequence.

---

## E. REPRODUCTION, DISCIPLINE AND WHAT LANE E SHOULD RE-DERIVE

```
python -m s24.a_corpus            # census, exclusion audit, null, artefact + hash   (~45 s)
python -m s24.a_ha2               # the fragment/peptide leak measurement            (~25 s)
python -m s24.a_charac cheap      # vocabulary, torsions, SS, window counts          (~2 s)
python -m s24.a_charac heavy      # leader curves + both controls; TAKES LOCK_TRAIN  (~33 min)
```

**Discipline observed.** Atomic writes (`tmp` + `os.replace`) on every artefact. Completion flags
gated on the **full key set** — all 126 target ids, all four criteria per held-out set, all eight
null lengths, all eight length bands, all four τ, both controls, and
`len(exclude) + len(permitted) == len(corpus)` — never a row count. All four artefacts report
`complete=true`. Seeds from `s15.seed.stable_rng`. Config-derived filenames carrying the corpus hash.
Results under `s24/results/`. `LOCK_TRAIN` taken with `os.open(..., O_CREAT|O_EXCL)` for the heavy
run only, released in a `finally` block; verified absent afterwards. `LOCK_AMBER` never touched.

**Statistics discipline.** Everything in D1 and D2 is a **census**, not an estimate: the whole corpus
and the whole target sets are enumerated. Censuses are reported without confidence intervals and are
labelled as such. **No comparison in this workstream has an SE**, so MDE = 2.8016 × SE is not quoted
anywhere here — quoting it would be the error `mde-is-per-comparison-not-per-instrument` warns
about. The one place a CI *will* be needed is §D3.0: any downstream comparison over the 126 dev
targets must use fold-clustered CIs, because five distogram models each saw ~100 of the other 125
dev natives and the per-target results are therefore not independent.

**Stake, restated.** I declared before running that I had a stake in D2's answer, since "the corpus is
adequate" would keep Lane C alive. The answer came back NO, against my own interest as declared. Lane
E has been tasked to re-enumerate Run 2's forks; **the specific things worth re-deriving are (i) the
10^4 threshold in my label, which is the only arbitrary number in the NO, and (ii) whether the leader
count is the right readout at all, since it is a covering number and not an information content.** I
have no stake in D1 — a longer exclusion list costs me nothing and costs the sprint compute — and
that enumeration stands.

**Falsifiers I named in advance and how they landed.**

| pre-registered hypothesis | outcome |
|---|---|
| H-A1 the shipped universes contain held-out-derived material | **CONFIRMED**, at 4.58% of windows; and 4/126 + 2/60 targets carry a verbatim self-copy |
| H-A2 the contamination enters through the FRAGMENT half | **FALSIFIED.** Fragments 0/210. It enters through the peptide half. |
| H-A3 the corpus is a small number of folds repeated | **FALSIFIED at τ=1.0 Å** (no saturation at any length); confirmed only at τ=2.0 Å for L=9 |

Three hypotheses, two falsified. The exclusion list, which is the deliverable two lanes were blocked
on, does not depend on any of them.
