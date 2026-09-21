# PREREG — S32 LANE L (LENGTH GENERALISATION)

**Committed before the first comparison. Lane L owns the charter's closing question,
*"maybe test on longer proteins"*, and S31's recommendation that the chiral-functional
emptiness "is worth retesting at 40+ residues on an instrument that does not exist yet".**

Contract rule 11 is absolute and is restated here as this lane's first constraint:
**nothing in this lane opens, alters, re-splits or contaminates `tuning126` or
`benchmark60`.** Every object lane L builds is separately named, carries its own fold
structure, and **no S32 headline about the 3.2105 Å endpoint may be computed on it.**

---

## 0. What was already measured before this document existed, and why that is legitimate

L1 is a **census plus one geometric property**. It contains no arm, no comparison, no
selection and no hypothesis test, so no pre-registration is possible or meaningful for it:
it is the inventory that decides whether the rest of the lane can exist at all. Those
numbers are listed in §1 as INPUTS. Everything from §2 onward is pre-registered.

Files already written: `s32/s32_L_feasibility.py`, `s32/s32_L_floor.py`,
`s32/results/L1_census.json`, `L1_peptide_floor.json`, `L1_windows.json`,
`L1b_floor_vs_length.json`.

---

## 1. INPUTS — the L1 census (no hypothesis, no comparison)

### 1a. Candidate supply

Length-*n* windows the retrieval stage can serve, by source:

| n | deployed library (peptide_db + fragment_db) | `prots/` (13,751 chains) |
|---|---|---|
| 13 | 17,853 | 2,153,224 |
| 26 | **0** | 1,974,465 |
| 35 | **0** | 1,850,715 |
| 45 | **0** | 1,713,271 |
| 55 | **0** | 1,576,535 |

The deployed library is capped at 25 residues (`peptide_db`) and 20 (`fragment_db`,
`LENGTHS = range(9, 21)`). **At n ≥ 26 it serves exactly zero candidates**, so no longer
instrument can reuse the deployed retrieval bank; one must be built from `prots/`, where
supply is ~96× the deployed supply at n = 13 and is not a constraint at any length ≤ 55.

### 1b. Target supply

`prots/` chains passing the same quality gates the 9–16mer library applies (contiguous
CA 3.5–4.1 Å, finite interior torsions, no X):

| band | chains |
|---|---|
| 26–35 | 2 |
| 36–45 | 25 |
| 46–55 | 144 |
| 56–65 | 200 |
| 66–80 | 228 |

599 in 26–80 total; 191 rejected for chain breaks. Redundancy clustering will reduce this
and the surviving count is a §2 deliverable, not an input.

### 1c. The representability floor — the number that reframes the lane

**DEFINITION, carried with the value (contract rule 4):** the *representability floor* is
the CA-RMSD between a native CA trace and the ideal-geometry chain
`build_backbone(φ_nat, ψ_nat)` built from **that same native's own torsions**. Fixed bond
lengths, fixed bond angles and fixed ω are the only structures this architecture can emit,
at every stage: retrieval windows are re-expressible, the projection's output is
`build_ca(φ, ψ)` by construction, and the AMBER stage is restrained to it. The floor is
therefore an **upper bound on the distance from the native to the set of emittable
structures**, and it is an ORACLE quantity used as a diagnostic only.

It is **not** S29's "architectural ceiling" (2.9027 Å), which is the best ORACLE prefix-*m*
average over the top-128 — a readout quantity on `tuning126`. The two are different
objects and are never differenced.

Measured, paired within **the same 371 proteins** at every length (so corpus, molecule,
resolution and admission policy are all held fixed and only *L* moves):

| L | mean | median | SE |
|---|---|---|---|
| 9 | 0.430 | 0.433 | 0.0036 |
| 13 | 0.700 | 0.699 | 0.0070 |
| 16 | 0.926 | 0.926 | 0.0099 |
| 20 | 1.260 | 1.257 | 0.0148 |
| 26 | 1.797 | 1.776 | 0.0223 |
| 35 | 2.619 | 2.560 | 0.0338 |
| 45 | 3.474 | 3.383 | 0.0471 |
| 55 | 4.225 | 4.088 | 0.0583 |
| 70 | 5.212 | 5.047 | 0.0724 |

`floor ≈ 0.0299 · L^1.239`, R² = 0.9958.

On the canonical 126 targets themselves the same quantity is **0.347 Å mean / 0.272
median** — 11% of the 3.2105 Å endpoint, which is why it has never mattered.

**SELF-TEST (contract rule 5).** `s32_L_floor.selftest` runs the same estimator on chains
*built by* the ideal-geometry builder from random torsions. These lie exactly on the
manifold, so a correct estimator must return ~0 at every length; an estimator that grew
with L for any other reason — the failure mode that would manufacture this result — fails.
Measured: 4e-15 (L=9) to 9e-14 (L=70), worst 1.8e-13, and it still moves for a single
0.5 Å displacement. The test can fail and does not.

---

## 2. HYPOTHESES

### L-H1 — THE FLOOR HYPOTHESIS (primary, and it reframes the lane)

**HYPOTHESIS.** The dominant obstacle to running this architecture at 40–60 residues is not
candidate supply, not target supply and not selection skill. It is that the *set of
structures the architecture can emit* recedes from the native as L^1.24, so at 40–60
residues the representation alone consumes more error than the entire current endpoint.

**MECHANISM.** Real backbones vary their N–Cα–C angle (τ) by several degrees and their ω
from 180° by a few degrees. Each is a small local error, but it acts as a lever arm on
everything downstream of it in the chain, so the CA displacement accumulates
super-linearly in L. At 9–16 residues there is no lever arm to speak of; at 45 there is.

**PREDICTION.** The **deployed projector** (`core.project.lam_path` at λ=0, multi-start,
the exact operator stage 3b uses), applied to the **native CA trace itself**, will return a
residual that rises with L on the same power-law scale, and at L ∈ [40, 60] will exceed
0.5 × the 3.2105 Å endpoint.

**FALSIFIER.** If the deployed projector on the native returns < 1.0 Å at L = 45 — i.e. if
the native-torsion rebuild is a loose bound and a much closer ideal-geometry chain exists —
then L-H1 is false, the representation is not the obstacle, and the lane proceeds to L-H2
with no reframing.

**CONTROL (rule 6 — matched to the operator's own space).** The identical projector, the
identical multi-start schedule and the identical penalty, run on the canonical 126 natives.
Both arms are the *same operator on a native*; only L differs. Reported as a pair.

**SECOND CONTROL — the AMBER escape.** Stage 4 is a restrained ff14SB/GBn2 minimisation and
is the one stage that can leave the ideal-geometry manifold. If it recovers the floor, the
floor does not bind at the endpoint. TEST: AMBER-relax the projected native under the
deployed schedule (k = 10, steps = 0) and re-measure. **Registered in advance: a recovery
of more than 50% of the floor at L = 45 falsifies the claim that the floor binds at the
endpoint**, and the lane will say so.

**DEPLOYMENT CONDITION.** None — this is a diagnostic, ORACLE / NOT DEPLOYABLE. Its output
is a *scope statement* on every conclusion this project has drawn at 9–16 residues.

### L-H2 — THE LADDER-SHAPE HYPOTHESIS (the scientific core, L2)

**HYPOTHESIS.** The *shape* of the loss ladder — pool contains the answer, selection
destroys most of it, projection costs in proportion to how non-physical the readout's
output is — is a property of the architecture and not of peptide length, and therefore
survives at 40–60 residues once the ladder is expressed **relative to the representability
floor of its own instrument**.

**MECHANISM.** Each rung is a different operator (retrieval truncation, prefix truncation,
argmin-vs-average readout). None of them references length. If they are length-free, the
ladder's *proportions* transfer even though every absolute number degrades.

**PREDICTION, registered in advance as the falsifiable claim.**
On the new instrument, with every rung measured on the **built chain**:
- **P1.** Pool headroom survives: `best single member (K=500)` is at least 0.75 Å below
  `PRODUCTION-equivalent`, i.e. generation is still not the bottleneck.
- **P2.** Selection is still the largest single rung: `selection / readout` is the biggest
  of the three deployable rungs (retrieval filter, prefix, selection).
- **P3.** The projection cost still tracks how far the readout's output is from a valid
  chain: a real member projects for < 0.05 Å and a dense 75-member average projects for
  more.

**FALSIFIER.** P1 fails (pool headroom < 0.75 Å) ⇒ generation *is* the bottleneck at
length and the 9–16mer conclusion does not transfer. P2 fails ⇒ the "selection binds"
conclusion is a peptide artefact. P3 fails ⇒ the projection-cost law (contract rule 16) is
length-scoped. **Each is recorded as a separate verdict; a mixed outcome is reported as
mixed, not rounded to the hypothesis.**

**CONTROL.** Every rung on the new instrument is computed by the **same code path** used to
compute it on `tuning126` wherever the code is length-agnostic, and any rung that needs new
code ships its own zero-information control (constant α-helix, per memory
`zero-information-control-must-be-plausible`).

**STATISTICAL RULE.** `s24/stats_lib.compare`, LOWER IS BETTER, `folds` supplied from the
new instrument's own fold structure (§3). **MDE = 2.8016 × SE, per comparison, stated
beside every mean.** < 0.7× MDE = NOT A RESULT; 0.7–1.0× = NOT MEASURED.

**DEPLOYMENT CONDITION.** None. Every ORACLE rung is labelled ORACLE / NOT DEPLOYABLE.
Lane L proposes no deployable change to the canonical pipeline.

### L-H3 — LENGTH-SUSPICIOUS PRIOR FINDINGS (L3, exploratory, declared as such)

Registered as **exploratory**, one comparison each, logged to `MULTIPLICITY.md`:
- **L3a.** Is the common-mode fraction of pool error length-dependent? (At 9–16 a whole
  target is one retrieved window; at 45 it is not.) Measured as the same statistic the
  `pool-error-is-68-percent-common-mode` entry defines, on the new pool.
- **L3b.** Does a constant α-helix still beat the torsion channel at 40–60 residues?

**No L3 result may be quoted as overturning a `tuning126` finding**, because it is measured
on a different instrument. It can only say "this finding is or is not length-portable".

### L-H4 — DEPLOYABILITY AT LENGTH (L4)

Descriptive, not a hypothesis: projection wall-clock and multi-start branch degeneracy at
L = 45 vs L = 13, candidate-bank cost, and the MDE the new instrument can actually resolve.

---

## 3. THE NEW INSTRUMENT — defined here, BEFORE any result is computed on it

**Name: `long40`.** It is never merged with, compared against, or substituted for
`tuning126`.

**Target admission**, in this fixed order, all criteria native-free except the quality
gates which read only the deposited backbone:
1. chain from `prots/`, first model, first chain, **40 ≤ n ≤ 60**;
2. contiguous CA (all steps in 3.5–4.1 Å), finite interior torsions, no `X`;
3. single-chain deposition (`nchains == 1`) — a monomer, so the deposited conformation is
   not held by an interface the prediction cannot see;
4. **leakage:** the target's own PDB contributes nothing to `peptide_db` or `fragment_db`,
   and no member of either bank, nor of the distogram's training set, has sequence identity
   ≥ 0.4 to the target;
5. **redundancy:** single-linkage cluster all survivors at identity ≥ 0.4; keep exactly one
   representative per cluster, chosen by the **lowest PDB code** (a rule fixed here, before
   any RMSD is known, so it cannot be chosen on outcome).

**Fold structure.** 5 folds assigned by **sequence-identity clustering of the admitted
targets only**, balanced by cluster size, seeded at 0, **written to
`s32/results/long40_folds.json` and frozen before the first arm runs.** The canonical
`peptide_folds.json` is not read and not written.

**Candidate bank.** All length-n windows of `prots/` chains, excluding the target's own
chain and any chain with identity ≥ 0.4 to the target, scored by BLOSUM62 sum with a
**stable** argsort (the pinned convention), top K = 500.

**Power, stated before any arm runs (contract rule: state the MDE it can resolve).** The
instrument admits an expected 40–90 targets. Taking the canonical instrument's paired
per-target sd as the scale, a paired comparison at n = 60 resolves
**MDE ≈ 2.8016 × sd/√60**. The realised per-target sd is unknown until the first ladder
column exists, so **the first thing reported after the ladder's PRODUCTION column is the
instrument's realised SE and MDE**, and any rung smaller than that is reported as NOT
MEASURED. Lane L commits in advance to reporting the MDE even when it is embarrassing: if
`long40` resolves only ~0.5 Å, the L2 predictions P1–P3 are gated at that resolution and
the sub-threshold ones are reported as NOT MEASURED, not as negatives.

---

## 4. WHAT LANE L WILL NOT CLAIM

- It will not claim any change to the 3.2105 Å endpoint.
- It will not claim a `tuning126` finding is refuted; only that it is or is not portable.
- It will not present an ORACLE rung as a deployment candidate.
- **It will report "the instrument cannot support a verdict" if that is the outcome**, with
  the specific quantity that blocks it, rather than manufacturing a weaker instrument that
  produces unresolvable numbers.
