# SPRINT 24 — PREREG A — DATA, CORPUS AND LEAKAGE

Written **before** any Workstream-A run. Sent to the coordinator before execution, per BRIEF §4
Rule 0. Enumerator: Workstream A. **STAKE DECLARED:** A has a stake in D2's decisive question
(a "corpus is adequate" answer keeps Lane C alive). Per Rule 0 that is a weakening; Lane E should
re-enumerate the D2 forks. A has **no** stake in D1 (the exclusion list): a longer exclusion list
costs A nothing and costs the sprint compute, so the incentive there is neutral-to-adverse and the
enumeration is safe.

---

## RUN 1 — THE CORPUS CENSUS AND EXCLUSION AUDIT (`s24/a_corpus.py`)

**Not a directional hypothesis.** It is an enumeration: what source structures exist, which of
them are related to a held-out target, and by which criterion. Forks are still named because the
CRITERION is an operator choice and the project's most repeated error is an unstated one.

### The six axes

| axis | TAKEN | ALTERNATIVE NOT TAKEN |
|---|---|---|
| **functional** | leak criterion = union of {exact sequence equality, verbatim substring in EITHER direction, Needleman–Wunsch identity ≥ 0.6 normalised by the longer sequence, source-PDB-ID equality, source-PDB-ID equality after chain stripping} | `containment ≥ 0.6`. **REJECTED WITH EVIDENCE**, not by taste: `containment-threshold-is-at-the-null` measures composition-matched random 9–16mers at 0.56–0.63 containment against this exact bank, i.e. the threshold is AT the chance level. Also not taken: structural homology by TM-score/CA-RMSD to the benchmark natives — that would require opening benchmark structures, which the mandate forbids. |
| **basis** | audit at the level of the **SOURCE CHAIN** (a `peptide_db` entry, or the parent deposit of a `fragment_db` fragment), then propagate the verdict to **every window derived from it** | window-level audit only. NOT TAKEN because a clean 9-mer window of a contaminated parent still carries that parent's conformation, and a window-only audit would pass it. |
| **readout** | the exclusion list is the set of **source identifiers removed**; reported as (a) source count, (b) window count, (c) per-target fraction of that target's universe removed, (d) worst-target fraction | "number of contaminated targets" alone. NOT TAKEN: a count of targets hides how much corpus each one removes, which is what Lane B/C actually need. |
| **normalisation** | identity normalised by the **LONGER** sequence — `peptide_db.identity`, the production rule that `folds()` and `clusters()` already use | normalising by the shorter sequence (= containment). NOT TAKEN: it is the same quantity memory records as sitting at the null, and it would also make the audit incommensurable with the fold assignment that everything else in the repo depends on. |
| **null** | every threshold I quote is accompanied by its **measured chance level**: 12 composition-matched random sequences per length 9–16, seeded `s15.seed.stable_rng("s24A","null",n)`, scored against the identical bank by the identical criteria | quoting a threshold without its null. This is the project's sixth-recorded instance of that error and it is not repeated here. |
| **THE LABEL** | **"leak" = a training-corpus source chain that IS a held-out target, CONTAINS one verbatim, is ≥0.6-identical to one, or is a window of the SAME DEPOSITED ENTRY as one.** Held-out = the 126 dev targets ∪ the 60 sealed benchmark ∪ the 24 dev_set. | "leak = the generator's samples resemble the native." NOT TAKEN: that is an outcome label, it can only be evaluated after training, it is circular (a good generator would be labelled a leak), and it cannot gate a corpus that must exist before training starts. |

### Pre-registered outcome-readings (written before the run)

* **H-A1: the shipped window universes already contain material derived from held-out targets.**
  FALSIFIER: zero source chains in any of the 126 universes match a benchmark or dev-set target
  under any of the five criteria. If the falsifier fires I will say so and the exclusion list will
  be empty, which is a real and reportable outcome.
* **H-A2 (mechanism-specific): the contamination, if it exists, enters through the FRAGMENT half,
  not the peptide half** — the peptide half is fold-filtered by identity cluster and so is already
  protected, whereas `_fold_fragments` filters fragments only against *the held-out fold's
  peptides*, and the 60 benchmark targets sit in **other** folds. FALSIFIER: contamination is
  found in the peptide half at a comparable or greater rate.
* I commit to publishing the count under **every** criterion, including criteria that return zero,
  so the reader can see which operator produced the headline number.

---

## RUN 2 — CORPUS CHARACTERISATION AND THE DECISIVE INDEPENDENCE QUESTION (`s24/a_charac.py`)

**This one is directional.** Lane C exists only if the answer is "yes".

### The six axes

| axis | TAKEN | ALTERNATIVE NOT TAKEN |
|---|---|---|
| **functional** | independence measured as **near-duplicate collapse under CA-RMSD after Kabsch** on the permitted windows: a window is a duplicate of an earlier one if CA-RMSD < τ. Report the whole curve τ ∈ {0.5, 1.0, 1.5, 2.0} Å rather than one τ | independence by torsion-space distance (circular RMSD on φ/ψ). NOT TAKEN as the primary because the generator's endpoint is Cα-RMSD and `control-must-match-the-operators-space` says the control must live in the operator's space; the torsion-space version is reported as a SECONDARY and any disagreement is reported. |
| **basis** | window level, **length-stratified** (9,10,…,16 separately) — windows of different length are never compared, mirroring `basis discipline` | pooling all lengths and comparing across them, which would invent similarity from length alone. |
| **readout** | **effective independent count** = number of greedy-clustering leaders at each τ, plus the entropy-based effective sample size `exp(H)` of the cluster-size distribution | raw window count (18,674/target). NOT TAKEN: it is the number the corpus advertises and it is the number most likely to be wrong. |
| **normalisation** | per length, per **permitted** corpus (after Run 1's exclusion), and reported both as an absolute count and as a **ratio to the parameter count** of the smallest plausible generator | reporting only the absolute count with no model-size referent. |
| **null** | two zero-information controls in the operator's space: (i) windows built from a **constant α-helix** at the same lengths, (ii) windows built from torsions drawn from the **per-residue-type marginal** Ramachandran table (`s8/generate_rama.npz`), which is plausible rather than degenerate. NOT uniform-on-the-torus. | uniform random torsions. NOT TAKEN — `zero-information-control-must-be-plausible` records that uniform-on-the-torus is a WORSE measure, not an uninformative one. |
| **THE LABEL** | **"the corpus supports a from-scratch generator" = the permitted corpus contains ≥ 10× as many τ=1.0 Å-independent windows per length band as the parameter count of the smallest generator Lane C would plausibly train (taken here as 10^5 parameters ⇒ 10^6 independent windows needed for the naive rule, 10^4 for a 10^3-parameter model), AND the φ/ψ distribution shows genuine per-residue-type multimodality rather than a single dominant basin.** Both halves are stated now so neither can be quietly dropped. | "supports a generator" = "there is a lot of data". NOT TAKEN: it has no referent and cannot fail. |

### Pre-registered outcome-readings

* **H-A3: the permitted corpus is effectively a small number of folds repeated.** Direction stated
  in advance: I EXPECT the τ=1.0 Å independent count to be one to two orders of magnitude below
  the raw window count, and I expect the α-helix basin to dominate. If it is not, I report the
  falsification.
* **I state now what would make me answer NO to the decisive question**, so that answering NO later
  is not a post-hoc rescue: fewer than ~10^4 independent windows in ANY length band, or >80% of
  windows in one Ramachandran basin class per residue type, or an independent count that does not
  exceed the plausible zero-information control's by a clear margin.
* Every reported difference carries SE and effect/MDE with MDE = 2.8016 × SE. Counts are censuses,
  not estimates, and are reported without a CI, labelled as such.

---

## RUN 3 — THE LEAKAGE-RISK REGISTER (`s24/agentA_FINDINGS.md` §D3)

No hypothesis; it is a design document written **before** either generator is built, per the
mandate. Its one measurable component — *what the distogram's own training set was, and whether it
overlaps the dev targets* — is executed as part of Run 1 and inherits Run 1's forks.

---

## DISCIPLINE COMMITMENTS

* Atomic writes (`tmp` + `os.replace`) everywhere. Config-derived filenames.
* Completion flags gated on the **full key set** (all 126 target ids present, all criteria present),
  never a row count.
* Seeds from `s15.seed.stable_rng`.
* Results under `s24/results/`.
* No heavy job without `s24/results/LOCK_TRAIN` taken via `os.open(..., O_CREAT|O_EXCL)`. Run 1 is
  not heavy (a few thousand alignments, seconds). Run 2's all-pairs RMSD **is** heavy and takes the
  lock.
* **Benchmark discipline:** benchmark TARGET IDENTIFIERS and SEQUENCES are read for the sole purpose
  of exclusion. No benchmark result file, no benchmark RMSD, no benchmark structure is opened. If an
  exclusion cannot be made without reading results, A stops and says so rather than guessing.
