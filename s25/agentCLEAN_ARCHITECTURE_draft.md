<!--
DRAFT — S25 CLEANUP LANE. NOT INSTALLED.

This becomes ARCHITECTURE.md after the freeze. Sections marked
<!-- PLACEHOLDER: FREEZE --> depend on decisions the research lanes have not made yet;
each one names the specific fact it is waiting for. Everything NOT so marked is measured
from the code as it stands on 2026-09-08 and should survive the freeze unchanged.
-->

# ARCHITECTURE

How a sequence becomes a structure, what every component actually computes, and where the
system is known to be limited.

Companion documents: `README.md` (install and run), `FINDINGS.md` (what the project
learned, negative results included), and each sprint's `LEDGER.md` (the evidence for a
particular decision).

---

## 1. DATA FLOW

```
    sequence (8–20 residues)
        │
        │  core.data.retrieve
        ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │ STAGE 1  RETRIEVE                                                     │
  │   every length-n window of the leakage-safe library                   │
  │     = out-of-fold peptides ∪ this fold's protein fragments            │
  │   scored by BLOSUM62 similarity sum, top K = 500 kept                 │
  └───────────────────────────────────────────────────────────────────────┘
        │  pool: (500, n, 3) Cα, plus φ/ψ and residue codes
        ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │ STAGE 2  FILTER            core.predict                               │
  │   a sequence-conditioned distance posterior over every (i, j) pair    │
  │   scores each candidate by L1 Bayes risk; top M = 75 survive          │
  └───────────────────────────────────────────────────────────────────────┘
        │
        ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │ STAGE 2b EVALUATE          core.energy (H_Legacy) · core.amber (H_AMBER)│
  │   two physical energies, rank-standardised into a common scale        │
  └───────────────────────────────────────────────────────────────────────┘
        │
        ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │ STAGE 2c SELECT            core.quantum                               │
  │   CVaR-VQE over a Hamiltonian-derived objective; the CVaR tail of the │
  │   measured distribution is the selected ensemble                      │
  └───────────────────────────────────────────────────────────────────────┘
        │
        ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │ STAGE 3  SYNTHESISE        core.geometry · core.project               │
  │   3a  coordinate-average the ensemble about its medoid                │
  │   3b  multi-start projection onto the ideal-geometry manifold under a │
  │       fold-disciplined Ramachandran prior                             │
  └───────────────────────────────────────────────────────────────────────┘
        │
        ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │ STAGE 4  RELAX             core.amber                                 │
  │   ff14SB / GBn2 restrained relaxation on the OpenMM CPU platform      │
  └───────────────────────────────────────────────────────────────────────┘
        │
        ▼
    structure  +  a per-target record  +  (labels only, after the fact) RMSD
```

> <!-- PLACEHOLDER: FREEZE -->
> **Component order.** The brief fixes what must be present, not the sequence. If an
> experiment justifies moving the selector before the filter, or the relaxation inside the
> ensemble loop, this diagram and §5 change together. *Needs: the frozen order.*

The evaluation path is deliberately separate from all of this. `s12/instrument.py` reads
the cached window universes, scores candidates, and computes every RMSD in the project.
**It never participates in a decision.** Native coordinates enter through `label()` only,
after the structure is final; `tests/test_pipeline.py` NaN-poisons every native quantity
and asserts that every emitted coordinate is bit-identical.

---

## 2. THE DISTANCE POSTERIOR (`core.predict`)

### What it predicts

For a sequence of length *n* and every pair (i, j) with |i − j| ≥ 2, a discrete
distribution over Cα–Cα distance bins. Inputs are a sequence representation (ESM-2
embeddings reduced to 32 principal components, plus one-hot identity and separation
features) and retrieved-fragment statistics. Five models, one per pinned fold; a target is
always scored by the model that did not see its cluster.

### How a candidate is scored — the form matters more than the model

Not by likelihood. By the **L1 Bayes risk** of the predicted distribution, evaluated on a
fine distance grid:

```
    r_ij(d)  =  Σ_b  p_ij(b) · | d − c_b |            c_b = bin centre
    grid     =  2.00, 2.05, …, 39.95 Å
    score(X) =  mean over pairs of  w_ij · r_ij( ‖x_i − x_j‖ )      lower is better
```

precomputed as a table so that scoring a batch of candidates is a gather rather than a
convolution.

Two consequences are load-bearing and both are measured:

* **The minimiser of an L1 risk is the posterior MEDIAN.** This is a *location*-based
  ranker. Posterior WIDTH is very nearly irrelevant to it.
* Consequently, calibrating the posterior's width does not help the endpoint. S25-L1
  measured the posterior as almost perfectly centred (z_mean −0.052) and **twice too
  narrow** (z_sd 1.996 against a nominal 1.0; 90% coverage actually 0.6159). S25-L2 then
  fixed that — a calibration-fitted widening drove held-out z_sd to 0.9834, essentially
  perfect — and RMSD got **worse by +0.0538 Å**. Variance calibration, confidence
  calibration, sharpening, heavy-tail treatment and entropy approaches are closed by
  measurement, not by argument.

A third measured property that constrains what can be done about it: the signed error
grows monotonically with separation, −0.048 Å at |i−j| = 2 to −0.589 Å at 9–15, and 24.1%
of pairs have a multimodal posterior. But S24 Workstream C measured that this offset is
**not a distogram defect** — real protein windows carry the same offset against these
natives — so it is a corpus/native scale mismatch that the candidates share, and
correcting only the prior may break a cancellation that currently helps.

---

## 3. THE TWO HAMILTONIANS

### 3.1 `H_Legacy` — the 11-term knowledge-based potential (`core.energy`)

Eleven decomposed terms at `DEFAULT_WEIGHTS`, **never fitted and never proxied**:

| term | weight | term | weight |
|---|---:|---|---:|
| `steric` | 4.0 | `electrostatic` | 1.0 |
| `contact` | 1.0 | `aromatic` | 0.8 |
| `hbond_local` | 1.0 | `torsion` | 0.15 |
| `hbond_longrange` | 3.0 | `compactness` | 0.4 |
| `coop_helix` | 2.0 | `solvation` | 0.5 |
| `coop_sheet` | 2.0 | | |

`BatchLegacy` builds the per-sequence tables once and scores a whole batch. `SUBSETS`
names term families (`hb`, `burial`, `packing`) used to test independence from the
distance prior — `hb` is pure backbone geometry with no amino-acid identity in it at all,
which makes it the cleanest such test.

### 3.2 `H_AMBER` — ff14SB / GBn2 through real OpenMM (`core.amber`)

`amber14/protein.ff14SB.xml` + `implicit/gbn2.xml`, CPU platform, `DeterministicForces`
on. There is no fallback and no proxy: AMBER is the validity stage. Side chains are built
for all 20 residues (`sidechains.py`), which is what took benchmark coverage from 6/35 to
35/35 targets. Bit-exactness against the reference implementation is asserted by
`tests/test_amber.py`, and frame invariance — that relaxation does not depend on the input
coordinate frame — by `tests/test_amber_frame_invariance.py`.

`core.amber` refuses to open a context above 92% physical memory (`memory_guard`). That is
a real operating constraint on a 15.6 GB box, not defensive decoration.

### 3.3 Normalisation between them — stated mathematically, because it has to be

Raw AMBER on unrelaxed retrieved windows is **clash-dominated**: 57.5% of a pool is above
1e4 kcal/mol and the worst single value measured is 7.1e18. A z-score over that
distribution puts 99.4% of the pool inside |z| < 0.1 — the normalisation destroys the
ranking it was supposed to preserve. So the pipeline uses **rank standardisation**: each
Hamiltonian's raw values are replaced by their rank within the pool, mapped to a common
scale, before any combination. Every comparison between the two energies in this
repository is on ranks.

This is also why CVaR's translation-equivariance and positive-homogeneity matter
(`tests/test_cvar.py`): they are what licenses the rescaling at all.

---

## 4. THE VQE AND ITS CVaR (`core.quantum`)

This is a real variational quantum eigensolver, not enumeration wearing its name.

### 4.1 The circuit

`layers × (RY on every wire → CNOT chain → optional ring closure)` applied to |0…0⟩, in
three simulators that are each **exact** for their regime:

| simulator | regime | what it is for |
|---|---|---|
| PennyLane `lightning.qubit` | n ≤ 30 | the production device; `tests/test_quantum.py` asserts probabilities against it with `==` |
| `StatevectorCircuit` | n ≲ 20 | a batched dense statevector, for gradients |
| `MPSAnsatz` | no qubit ceiling | an **exact** matrix-product state — the chain topology caps the bond dimension at 2^layers, so nothing is truncated |

Probabilities come out of the simulated state, bitstrings are drawn from that
distribution, and the objective is evaluated on what was drawn. `all_bitstrings`,
`cvar_gradient_exact` and `grad_cvar_fd` enumerate the register; they are **verification
instruments**, they are named so, and `test_the_search_does_not_enumerate` asserts the
search touches a vanishing fraction of 2^n.

### 4.2 The CVaR objective

The conditional value at risk of the objective distribution over measured bitstrings —
the mean of its lower α tail, not the mean:

```
    CVaR_α(θ)  =  max_t  [ t − (1/α) · E_{x ~ p_θ} (t − E(x))_+ ]
```

At the optimal *t* (the α-quantile *q*) the envelope theorem kills the *t* dependence and
leaves the score-function form

```
    ∇ CVaR  =  −(1/α) · E_p[ (q − E(x))_+ ∇ log p_θ(x) ]
            =   E_p[ f(x) ∇ log p_θ(x) ],     f(x) = −(q − E(x))_+ / α
```

α is annealed geometrically over a run (`alpha_schedule`, 0.5 → 0.05). Early on, a narrow
tail concentrates the distribution before the landscape has been sampled — the premature
concentration failure. Late on, a broad tail is just the mean, and the mean is not what we
want to minimise.

### 4.3 The gradient-baseline defect, and why the test suite names it

A baseline *b* may be subtracted from *f* **if and only if it is constant in x**, because
the correction term is `b · E_p[∇ log p] = 0`. The original `qansatz.cvar_gradient`
subtracted the **tail mean from the tail entries only** and left the rest at zero — that is
`b(x) = m · 1[x ∈ tail]`, a *function of x*, and the identity it relies on does not hold on
a data-dependent subset. **It is a bias, not extra variance**, and it measured a cosine of
−0.023 against the true gradient.

`test_cvar_gradient_baseline_is_constant_not_tail_only` does not read the source to check
which baseline is in use. It measures both against an exact classical reference on an
enumerable register and asserts the shipped default lands on +1.000000 rather than
+0.655634. Reintroducing the defect fails with the number that names the bug.

### 4.4 An honest measurement about the VQE

At these problem sizes the VQE ties uniform random sampling and loses to annealing (S9
`stage_search`). On the lattice encoding, running it was measured as **worse than not
running it** — 0/12 cells against best-of-N from the untrained circuit. On the continuous
encoding it genuinely trains (−0.210 Å, 5/5 folds the same sign). The durable lesson is the
control discipline: **always control against best-of-N from the untrained circuit, never
against an initialisation mean.**

The component is a pillar and stays genuine. The honest number is reported rather than
engineered away.

---

## 5. AGGREGATION AND OUTPUT

### 5.1 Coordinate averaging (Stage 3a)

Superpose every ensemble member on the **medoid** (the member minimising mean pairwise
RMSD), then take the coordinate mean.

Two measured facts govern this stage:

* **Where you average beats what you rank with.** Coordinate averaging beats torsion
  averaging by 1.024 Å, and the choice of ranking objective contributes only 0.171 Å.
* **Averaging contracts the backbone by 25.8%**, which is what Stage 3b then has to undo.
  `tests/test_instrument.py` controls for this: averaging rigid copies of one structure
  must return that structure exactly, so a test cannot pass merely by collapsing.

The medoid is `np.argmin` over row means, so a tie resolves to the order the pool arrived
in. That is deterministic and therefore right, but it means anything ranking on a signal
with large tie sets inherits the candidate order — this project once had an `argmin` over
a tied native-free signal silently read the *oracle* sort order and invent a 1.386 Å
winner. Average the outcome over the tied argmin set.

### 5.2 The manifold projection (Stage 3b, `core.project`)

The averaged cloud is not a protein: bond lengths are wrong and the torsions are not
physical. Stage 3b finds the nearest ideal-geometry chain under a fold-disciplined
Ramachandran penalty, over a λ path `(0.0 → λ)` with multi-start and exact gradients.
The λ = 0 arm is the pure geometric fit and is reported beside the penalised one, because
the difference between them is the price of the prior.

### 5.3 Output and provenance

Every target writes one atomic JSON checkpoint keyed on a SHA-1 over every parameter that
can change a number — including **the live backend set**, so a consolidated result can
never be served out of a reference-arm cache. Anything promoted carries a provenance stamp
with the FULL KEY SET; `complete` is never inferred from a row count or a filename.

> <!-- PLACEHOLDER: FREEZE -->
> **The results layout.** `results/summary/{results.json, results.csv, leaderboard.csv,
> final_report.md, professor_brief.md}`, `results/structures/T001__method.pdb`,
> `results/site/`. *Needs: the results-lab lane's schema and the renderer's entry point.*

---

## 6. KNOWN LIMITATIONS

Stated plainly, because they are the most useful part of this document.

1. **Discrimination binds; search does not.** A ~1.9 Å selection gap survives the
   *certified* optimum on two disjoint target sets. Exhaustive enumeration exceeds the 2^n
   latent on 75 of 126 targets and exact argmin merely ties the pool, while an oracle over
   the same set is 1.60 Å better. The problem is not that the optimiser is too weak.
   Searching harder helps only on a good objective, is neutral on a mediocre one, and
   actively hurts on a bad one — always state which.

2. **Nothing ranks reliably within the pool.** The dev pool's best is 2.355 Å against
   3.324 Å returned. All-atom AMBER reranking is worth +0.004 Å. Thirty-eight native-free
   signals were tested; only typicality has positive in-band skill, and the best validated
   in-band operator (score-filter + consensus medoid) is −0.172 Å [−0.316, −0.027].

3. **The distance prior is the ceiling.** Oracle distances give 0.36 Å pool / 0.98 Å
   selected through the same library. Perfect distance knowledge caps the instrument at
   ~1.95–2.0 Å, and the trained predictor beats pool-typicality by 0.024 Å.

4. **Neither energy ranks the native lowest.** Not `H_Legacy`, not `H_AMBER`, and not
   AMOEBA — polarizable physics does not fix it either. On matched pools with an Rg
   control, AMBER puts the native at the **51st percentile**. The native is the distance
   objective's argmin on 3 of 126 targets (36.8th percentile). This is not a single-term
   defect; it is distributed.

5. **The posterior is over-confident by 2×, and fixing it makes RMSD worse** (§2). The
   shipped Bayes risk is a location-based ranker and width is nearly irrelevant to it.

6. **Sequence conditioning is worth 0.776 Å overall and is HARMFUL on the failures.** On
   the 18 hardest targets the blind pipeline (5.425 Å) *beats* the shipped one (6.019 Å).
   The peptide corpus carries ~7× the sequence–structure channel of the protein fragments
   that are 80% of every pool.

7. **φ carries no sequence signal at peptide length.** Full sequence context predicts φ at
   36.1° against 36.4° for a sequence-blind marginal. The whole channel is 10.4° of ψ.

8. **In-band ordering is learnable but per-target.** 0.986 within a target with zero
   overfitting, 0.600 across targets; reaching 2.0 Å needs 0.638. Capacity is saturated by
   a *linear* model. The only remaining leverage is supplying the per-target sign at
   inference, and native-free compactness proxies reach 0.24–0.37 against the oracle's
   0.909 — open, but weak.

9. **No fresh benchmark exists.** All 204 clusters of 9–16mers are spent. The
   containment-fresh world supply is 16 targets, 10 of them amyloid fibrils. And a 0.6
   containment threshold is *at the null* for peptides — random sequences score 0.56–0.63
   against the fragment bank — so containment must be checked by identity or verbatim
   substring, not by that threshold.

10. **A known filter defect stands.** `identity` normalises by the longer sequence, so a
    long fragment can carry a target's exact k-mer and still align below 0.6. Priced at
    +0.0030 Å on the benchmark. It is documented rather than fixed because fixing it would
    re-pin the folds, and re-pinning the folds once already moved 13 benchmark targets and
    invalidated every trained model.

11. **The instrument's numerical floor is ~1.3e-7 Å.** RMSD is reconstructed from singular
    values as `|P|² + |T|² − 2Σs`; for near-identical structures that is a cancellation of
    two ~1e4 Å² quantities. Four orders below anything reported, but a self-RMSD is never
    exactly 0.0 and code should not assume it is.

12. **`s12/instrument.py` duplicates four functions from `core.geometry`.** They agree
    bit-for-bit today and `tests/test_instrument.py` asserts it at `==`. One of them,
    `pairwise_rmsd`, is behaviourally different from its `core` counterpart (`core`
    symmetrises, the instrument does not). The published findings were measured with the
    instrument's version.

---

## 7. THE MODULE MAP

| module | responsibility |
|---|---|
| `core/__init__.py` | the backend switch: consolidated module if it satisfies `CONTRACT`, else the root module it replaces. Records the choice in `ACTIVE`, and that set enters every cache key. |
| `core/data.py` | peptide database, pinned folds, identity clusters, manifests, retrieval |
| `core/geometry.py` | backbone build, Kabsch, RMSD, secondary-structure assignment |
| `core/predict.py` | the distance posterior and its Bayes-risk table |
| `core/energy.py` | the 11 Legacy terms and their weights |
| `core/amber.py` | ff14SB/GBn2 through OpenMM; the memory guard; the builder cache |
| `core/quantum.py` | CVaR, the ansatz families, the three simulators, SPSA, the drivers |
| `core/project.py` | the manifold projection and the fold-disciplined Ramachandran prior |
| `core/cache.py` | on-disk caches and subset extraction from the large banks |
| `core/pipeline.py` | the stages, parallel and resumable |
| `core/bench.py` | the harness that measures them, per stage, against the `s9/final.py` reference |
| `s12/instrument.py` | the evaluation instrument — every RMSD in the project |

The 24 modules at the repository root are the `CORE_BACKENDS=legacy` reference arm: the
pre-consolidation implementations that every equivalence claim is proven against. The
import graph and the four reachability closures are in `s25/agentCLEAN_reachability.json`,
regenerable with `s25/agentCLEAN_reach.py`.
