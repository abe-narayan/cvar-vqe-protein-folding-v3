# SPRINT 31 LEDGER

Entries are `## S31-L<n> -- TITLE (date time, lane)` with the verdict in the heading.
Run `date` in the same command as the append. Corrections are annotated in place with the
original wording left standing.

Ledger numbers are a read-modify-write with no lock and **collided three times in S30**. Take the
next free number and re-check it immediately before writing; if you collide, renumber the
later-stamped entry and annotate in place.

---

## S31-L0 -- THE CHARTER (SAVED AND **VERIFIED**), THE BENCHMARK SEAL, AND THE ARCHITECTURAL READING THAT DETERMINED THE LANE MAP: **A 512-DIMENSIONAL STATE IS BEING USED TO SPECIFY ONE INTEGER** (2026-09-20 23:44, coordinator)

### The charter

Saved verbatim as `s31/BRIEF.md` — **1,660 lines, 47,348 bytes, existence and size verified by
reading the file back, not by asserting it.** S30-L0 made this same claim and the file was never
written; it stood false for 77 minutes and was the *fifth* instance in this project of prose naming
a path that does not exist — recorded, as it happens, by the very entry that committed it.
`s31/s31_verify.py` (lane D) asserts every path this ledger claims to have written.

One hard constraint (CVaR-VQE is the spine and the main scientific object), one endpoint
(**mean built-chain Cα RMSD on 126 targets, currently 3.2105 Å**), primary target < 3.00,
ambition < 2.50.

### Benchmark seal — verified, not assumed

`benchmark60` is sealed. `core/pipeline.py:386-388`: *"ITS SINGLE PRE-REGISTERED PASS IS SPENT — the
harness refuses to run it without `--i-am-spending-the-benchmark`, and it is never a tuning
instrument."* Two independent guard sites: `core/bench.py:1074` and `core/pipeline.py:1664`. The
working instrument this sprint is **`tuning126`**. Folds and clusters are pinned and will not be
regenerated.

### THE ARITHMETIC (carried from S30, unchanged — the mean is a tail statistic)

```
production, n = 126        mean 3.2105   median 2.9661
the worst 18 targets       mean 6.2758
the other 108              mean 2.6997

cap the worst 10 at 3.00 A  ->  2.9074  (-0.3031)   BEATS the primary target
cap the worst 18 at 3.00 A  ->  2.7426  (-0.4680)
cap the worst 30 at 3.00 A  ->  2.5778  (-0.6328)
```

And the tail is **selection-limited, not pool-limited**: on the genuinely worst 18 the ORACLE best
pool member is **2.5298 Å** with 11 of 18 under 3.00. The material is already in the candidate sets.

### THE ARCHITECTURAL READING, made from the code before any lane was briefed

`core/pipeline.py` ~838:

```python
E = _zrank(np.asarray(pool["sc"], float)[o])     # DIAGONAL: rank-standardised candidate scores
p, cvar, H, _circ = qm.run_cvar_vqe(E, alpha, T, n=cfg.vqe_qubits, ...)
block = np.asarray(Pt, float)[:dim, :dim]        # PAIRWISE distances between candidates
local = consensus_medoid(block, p)               # readout: p-weighted 1-median over `block`
```

**Two structural facts follow, and between them they set the sprint:**

1. **Theorem T1 (S30): for a diagonal `H`, the CVaR tail is always a prefix of the order induced by
   the objective.** The optimised state's entire structural output is therefore **one integer** —
   the prefix length `m`.

   > **A 9-qubit, 512-dimensional state is being used to specify a single number.** That is the
   > information-allocation defect at the centre of the architecture, and it is not a defect of
   > expressivity, optimisation or circuit depth.

2. **The Hamiltonian and the readout are decoupled.** `H` encodes *how good is each candidate
   alone* (diagonal, scores). The readout consumes `block` — *how structurally similar are
   candidates to each other* (off-diagonal). **The optimiser never sees the geometry it will be
   read out through.** This matters because the project's one established in-band discriminator is
   consensus/typicality (S12: score-filter + consensus medoid, −0.172 Å, the project's first CI
   excluding zero), and consensus is **irreducibly pairwise** — it cannot be written as a diagonal
   function of candidate index without precomputing it, which collapses it back to a diagonal.

And the fact that makes the whole thing worth attacking rather than abandoning: the same circuit
family reaches **0.2516 Å** built chain under an ORACLE objective and **3.4330 Å** under the
deployed native-free one, against a classical average's 3.2071. **The circuit can express excellent
states. The objective does not point at them.**

### THE LANE MAP — five launched, three slots held

| lane | remit | first question |
|---|---|---|
| **A** | quantum/CVaR theory | Characterise the objectives whose CVaR argmin is *not* a prefix. Which hypothesis of T1 fails off-diagonal — and does a generalised prefix theorem still apply? Then: is `H = diag(zrank) − λ·W(block)` physically meaningful, and is it classically reducible (charter §11)? |
| **B** | physical model | §7A free energy and §7B torsion. **Theoretical pre-check first:** `F = E − TS`; the enthalpic half is already closed (AMBER puts the native at the 51st percentile and is *worse* than chance as a ranker), so the novelty must live in `S` — and `S` must have a target argument that (sequence, pool) does not already carry |
| **C** | encoding / readout | §7C native-free sparse support; §12 candidate-index allocation; and the **128→512 widening's undischarged circularity gate**, which is cheap and either revives or closes a direction |
| **D** | verification / integrity | The four defects: the withdrawn positive in shipped code at `core/pipeline.py:821`; **the unpinned projection seed** (gates every sub-0.01 Å claim); the launcher/governor contradiction; and the verifier + sprint-wide multiplicity register |
| **E** | the incoherence hypothesis | See below — the sprint's sharpest new idea |

### THE NEW HYPOTHESIS, AND WHY I THINK IT IS THE SPRINT'S BEST SHOT (lane E)

S30's deepest result is that **what can be predicted is coherent with the pool's common mode and
therefore harmful**: at identical out-of-fold R² = 0.2355 a fitted corrector emits **+0.0554 Å
worse** while a synthetic i.i.d. one emits **−0.2466 Å better**. Every corrector the project owns
*raises* the residual's coherence with the common mode (0.6931 → 0.786 / 0.783 / 0.917).

S30 tested *"fit a corrector and apply it."* **It never tested "apply only the component orthogonal
to the common mode."**

And the enabling observation, which I believe S30 missed: since `mu = pool75_mean − d_nat` and the
distogram's `expected` is an estimate of `d_nat`, the quantity

```
mu_hat = pool75_mean - expected
```

**is computable at inference with no native at all** — and both terms are already features in
`s30_P_lr.py`. This is the project's standing "prediction/pool disagreement is a native-free signal"
arriving where it can actually be spent.

**My registered prior: roughly 2:1 against a deployable gain.** If the fitted corrector's residual
has coh 0.9172 it is *almost entirely* in the `mu` direction, so the orthogonal complement may be
nearly pure noise — giving the i.i.d. arm's *magnitude* without its *information*. **That is still a
clean result**, because it would show the i.i.d. benefit comes from genuine orthogonal signal rather
than from incoherence per se, which sharpens what a fourth observable must supply.

### WHAT I AM CARRYING FORWARD AS BINDING (until a lane breaks it)

- **Basis discipline.** Built chain **3.2105**; CA cloud **3.0483**; *set mean* 3.5507. Three
  different objects. The measured cloud→chain transfer for a correction is **0.92** — the 1.16
  coefficient in `operator-consumes-set-mean` maps *set mean* to output and **must not** be used to
  transfer a cloud delta (its own inputs give 4.21 Å against production's actual 3.21).
- **MDE = 2.8016 × SE, per comparison.** Below 0.7× is not a result; 0.7–1.0× is NOT MEASURED.
- **G1 (S30):** every achiral rotation/translation-invariant single-structure observable is a
  distance-map reading. Renaming one is not a new channel.
- **ORACLE labels travel inside the sentence carrying the number**, not a paragraph later.
