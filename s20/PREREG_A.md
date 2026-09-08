# SPRINT 20 — PRE-REGISTRATION, WORKSTREAM A
## Can any candidate generator produce structures whose COHERENT error decorrelates from the retrieval pool's?

Written **before any arm ran**. Not edited after seeing data (brief §7). Deviations are recorded in
`s20/agentA_FINDINGS.md` as deviations, never by editing this file.

---

## 0. INSTRUMENT AND REPRODUCTION GATE (fixed before the question)

126 cluster-disjoint tuning targets, pinned folds, `s12/instrument.py`. Full-chain Cα-RMSD, frozen
implementation, all residues, proper rotations. Sealed 60-target benchmark untouched.

**Verified before writing this file** (so the reference numbers below are mine, not quoted):

| object | mean Cα-RMSD |
|---|---|
| raw coordinate average of the shipped top-75 (`s19/cache/start_*.npz:avg`) | **3.0483** — "best built" |
| projection of that average, λ=0.3 (`I.project`) | 3.2126 |
| shipped production structure (`shipped_record(pdb)["fit_ca"]`) | 3.2041 — incumbent |

My reconstruction of the shipped top-75 — `pool_idx(u)` → `pair_dists` → `shipped_score` →
`argsort[:75]` — is **set-identical to `shipped_record(pdb)["sub"]`** on the target tested.
A full 126/126 set-identity check is a gate on the run and is reported.

> **Consequence recorded in advance**: the programme's best built object, 3.048 Å, is the
> **unprojected** coordinate average. The projection *costs* +0.164 Å. The **primary emitted object
> for every source in this pre-registration is therefore the raw coordinate average of that source's
> top-75**, referenced to 3.048. The projected variant is reported as the secondary
> "through-the-real-pipeline" number, referenced to 3.213.

---

## 1. HYPOTHESIS

**H-A.** The harmful coherent error component identified in Sprint 19 (shared 0.87 with the
retrieval pool, 0.81 cross-architecture, 0.64–0.67 with sequence-blind references) is carried by the
**corpus the candidate windows come from**. A candidate source drawn from a *different corpus* —
specifically the length-matched peptide database rather than the protein fragments that are ~70–85%
of every pool — will emit a structure whose coherent error is **less aligned** with the incumbent
pool's than a candidate source drawn from the *same* corpus at matched operation.

**H-A2 (the payoff).** If a source's coherent error decorrelates, coordinate-averaging that source's
emitted structure with the incumbent's cancels part of the coherent component, and integrated
Cα-RMSD falls.

---

## 2. EXPECTED MECHANISM

Every stage of the pipeline is fed by one window universe. The distogram is trained on fragments
from that universe; the candidate pool is retrieved from it; the selector scores against the
distogram. Sprint 19 measured the resulting error as coherent and shared, and named the corpus as
the untested carrier. If the peptide sub-corpus carries a materially different structural prior
(project memory: 7× the sequence–structure channel of the protein fragments), the two sub-corpora
should be *differently wrong* — their emitted coherent errors should point in different directions —
and the average of two differently-wrong structures is closer to the native than either.

If instead the coherent error is a property of **the target's own length and the generic peptide
manifold**, any two disjoint subsets of the universe are wrong in the same direction and the corpus
identity buys nothing over a random partition.

---

## 3. THE OPERATOR, STATED EXACTLY (so controls can be matched in its space — brief §6 rule 1)

For a candidate index set `Ω ⊆ universe`:

    rank Ω by the universe's own stable BLOSUM order `u["order"]` (filtered, so ties break identically)
    take the first 500                       -> the SOURCE POOL
    score all 500 with the shipped Bayes-risk distogram score `I.shipped_score`
    take the 75 lowest                       -> the SOURCE TOP-75
    X_S = I.coordinate_average(top-75)       -> THE EMITTED STRUCTURE (raw; primary)
    X_S^proj = I.project(X_S, seq, fold)     -> the projected variant (secondary)

**Every arm uses this identical operator.** Only `Ω` moves. The two non-retrieval arms (`tors`,
`helix`) substitute a generated 500-member candidate set for the retrieval step and are otherwise
identical from `shipped_score` onward.

### Error fields

    e_S = d(X_S) - d_true      over pairs (i,j), min_sep = 2      [ORACLE, evaluation only]

`X_S` is a real structure, so `e_S` **is** its coherent error — the Sprint-19 split `r = r_coh +
r_inc` is degenerate here with `r_inc = 0`. This is an identity, not a finding, and is stated so no
reader mistakes it for one.

    rho(S, T) = mean-centred Pearson correlation of e_S and e_T across that target's pairs,
                averaged over the 126 targets   (identical semantics to `s19/a_source._corr`)

`e_pool` := `e` of the incumbent arm (`pool`).

---

## 4. ARMS

**Retrieval partitions (disjoint window subsets of one universe, identical operator):**

| arm | Ω |
|---|---|
| `pool` | `order[:500]` — the incumbent, exact reproduction |
| `pep` | windows with `org == True` (peptide database) |
| `prot` | windows with `org == False` (protein fragments) |
| `halfA`, `halfB` | a **random 50/50 partition** of the universe, `stable_rng(pdb, "s20A_half")` |

`pep`/`prot` and `halfA`/`halfB` are **both** disjoint partitions of the same universe consumed by
the same operator. They differ *only* in whether the partition is by corpus or at random. This is
the matched control, in the operator's own space.

**Non-partition sources:**

| arm | what it is | role |
|---|---|---|
| `rand500` | 500 windows chosen uniformly at random from the universe | changes the *selection criterion*, not the corpus |
| `tors` | 500 chains from a classical continuous-torsion generator: a 3-state (H/E/C) first-order Markov chain with per-state Gaussian (φ,ψ), parameters from standard backbone geometry, **not fitted to this project's corpus**; sequence-blind, length is the only input | "not made of the pool" — brief candidate 3 |
| `helix` | the constant ideal α-helix (φ=−63°, ψ=−42°) | zero-information reference (brief §6 rule 4) |

**Fusion arms (the payoff test):**

| arm | construction |
|---|---|
| `fuse(A,B)` | superpose `X_B` on `X_A`, take the midpoint; primary and projected variants |
| `merge(A,B)` | coordinate average of the **union of the two top-75 sets** (150 members) |

Fusion is run for `(pool,pep)`, `(pool,prot)`, `(pool,tors)`, `(pool,rand500)` and — as the matched
control — **`(halfA,halfB)`**, which fuses two same-corpus halves at the same member count.

---

## 5. PRIMARY ENDPOINT

**Integrated Cα-RMSD of the emitted structure through the real pipeline**, n = 126, paired
target-level bootstrap CI, fold-stratified, median and W/L beside every mean.

The **decision metric for the hypothesis** is `rho(e_S, e_pool)`, reported for every source, with
member RMSD, geometric diversity, error diversity and near-native coverage beside it. A source that
decorrelates but does not lower integrated RMSD is reported as a decorrelated source that did not
convert, not as a win.

**Ceilings reported separately for every source, always four numbers:**

| ceiling | definition |
|---|---|
| generation | `min rr` over the source's 500-member pool (ORACLE) |
| selection | `min rr` over the source's top-75 (ORACLE); i.e. the best any native-free selector inside the band could reach |
| repair | the emitted structure's raw average vs its projected variant (the +0.164 Å repair tax, measured per source) |
| realised | Cα-RMSD of the emitted structure |

---

## 6. FALSIFIER (pre-registered; brief §6 rule 5)

**F-A1 — the corpus is not the carrier.**
The brief states the bar in absolute terms: *if the peptide-only pool's coherent error aligns with
the current pool's at ρ ≥ 0.85 — the cross-architecture level — the corpus is not the carrier and
this closes.* I adopt that, and because my ρ is measured between **structures'** error fields and
Sprint 19's 0.784/0.898 were measured between **distogram families'** coherent components, the two
are not on the same scale. I therefore pre-register **both** forms and state now which governs:

* **F-A1a (absolute)**: `rho(e_pep, e_pool) ≥ 0.85`.
* **F-A1b (share of ceiling)**: `rho(e_pep, e_prot) ≥ 0.85 × rho(e_halfA, e_halfB)`, i.e. the corpus
  partition decorrelates no more than a random partition of the same universe under the same
  operator.

**If the measured same-corpus ceiling `rho(e_halfA, e_halfB)` is itself below 0.85, F-A1a is
unusable as stated and F-A1b governs.** This is written before the data because I expect the
structure-level ceiling to sit lower than the distogram-family ceiling, and I refuse to be able to
choose afterwards.

**F-A1 firing closes candidate source 1 and source 2** (both are corpus swaps) and demotes the whole
"decorrelate by changing the corpus" branch.

**F-A2 — decorrelation exists but does not convert.**
If some source reaches `rho ≤ 0.6 × ceiling` **and** its fusion arm fails to beat `pool` by more than
the 0.084 Å MDE with a CI excluding zero, then decorrelation is real and worthless, and the lever is
reported as measured-and-closed rather than open.

**F-A3 — the matched control eats the effect.**
If `fuse(pool,pep) − pool` and `fuse(halfA,halfB) − halfA` are statistically indistinguishable, any
fusion gain is "averaging two structures", not "averaging two *decorrelated* structures", and must be
reported as the former.

---

## 7. NULL AND MATCHED CONTROLS

| control | what it nulls |
|---|---|
| `halfA`/`halfB` random partition | the corpus partition — matched in count, operator, selector and universe |
| `rand500` | BLOSUM retrieval as the selection criterion, at matched count |
| `helix` | zero-information reference structure (a *reference*, not a gate — the role-dependence in brief §6 rule 4 is noted) |
| `fuse(halfA,halfB)` | the fusion operator itself, at matched member count and matched corpus |
| set-identity of my `pool` top-75 with `shipped_record["sub"]`, 126/126 | my reimplementation of the operator |

---

## 8. BUDGET

Phase 1 (no projection): all arms, all 126 targets — seconds per target, well under one core-hour.
Phase 2 (projection, 5.5 s/call): only the arms Phase 1 makes decisive, ≤ 10 arms × 126 targets
≈ 2 h wall on one process. One heavy process, `OMP_NUM_THREADS=1`, load checked before launch
(measured at prereg time: CPU 61%, 3.89 GB free — so **one** process, not two).

---

## 9. PROMOTION CRITERION

A source or fusion arm is **promoted** only if, at n = 126:

1. integrated Cα-RMSD beats the `pool` arm's 3.048 Å by **more than 0.084 Å (the MDE)**, and
2. the paired bootstrap CI excludes zero, and
3. the sign holds in **≥ 4/5 folds**, and
4. it beats its own matched control (`fuse(halfA,halfB)` for fusion arms; `rand500` for source arms)
   past the MDE, and
5. it needs no native information anywhere in generation or selection.

Anything short of all five is reported at its own label (`NOT MEASURED`, `OPEN`, `REFUTED`), never as
an improvement. Per the brief: an argmin over arms on this tuning instrument is a hyperparameter
chosen here, and is labelled as such.

---

## 10. WHAT WOULD MAKE ME ABANDON THE LANE

If F-A1 fires **and** F-A2 fires, the "decorrelated generation" lever is closed for corpus swaps and
for the classical torsion generator, and I move to the secondary mission (ensemble-first
architecture, directive §30) without further spending here. I state in advance that this is the
outcome I consider most likely, given that the incoherent component is what averaging suppresses and
Sprint 19 measured the coherent component as shared even with references that know nothing about the
sequence.
