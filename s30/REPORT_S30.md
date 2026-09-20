# Sprint 30 — CVaR-VQE Protein Folding: the search for the first real accuracy breakthrough

**Status: FINAL.** Every section is written from artefacts, not from memory. Two independent
adversaries audited it under contract rule 28 (*the main team may not approve its own positive*):
lane V on §1/§3/§4/Appendix A (**18 defects**, two of which inverted a headline) and lane Z on
§0/§2/§5/§6/§9–§13 (**23 defects**, four severe, one of which mislabelled the sprint's deliverable).
**All are fixed in place with the original error stated.** `s30/s30_verify.py` recomputes **36/36**
headline numbers from the artefacts with 0 mismatches.

Branch `s26` · benchmark: 126 dev targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Ledger: `s30/LEDGER.md` · Running state: `s30/STATE.md` · Contract: `s30/S30_CONTRACT.md` (28 rules)
Charter: `s30/BRIEF.md`

---

## 0. The answer, up front

The charter asked for the first real accuracy breakthrough, with a primary target of < 3.0 Å and an
ambitious target of < 2.5 Å against production's **3.2105 Å**. It also said, explicitly, that *"a
well-evidenced ceiling argument that redirects the next three sprints is worth more than a fragile
2.98 Å."*

### 0.1 The endpoint did not move

| # | question | answer |
|---|---|---|
| 1 | Strongest new architecture? | **None was deployed.** Every candidate was closed — by measurement, by theorem, or by price — before it reached the pipeline |
| 2 | Final mean built-chain RMSD? | **3.2105 Å**, unchanged |
| 3 | Did it beat 3.21? | **No** |
| 4 | Did it beat 3.0? | **No** |
| 5 | Did it beat 2.5? | **No** |
| 6 | Exact paired effect vs production? | **Zero — nothing was deployed.** Every *applied* prior correction is **NOT MEASURED or worse**: the three fitted arms emit **+0.0554, +0.1461, +0.1513 Å** at **0.68×, 0.91×, 0.95× MDE** (CA point cloud, against 3.0483). The AMBER relax at k = 30 is **−0.0221 Å**, the one confirmed effect — **0.69% of baseline**, *not deployable* (its restraint constant has no native-free selection rule, §11.2), and **not closed** (§3, §7.4) |
| 7 | Fold-clustered CI? | Per comparison; every one is reported beside its effect. No sprint-wide constant |
| 8 | MDE? | Per comparison, MDE = 2.8016 × SE. **Below 0.7× is not a result; 0.7–1.0× is NOT MEASURED** |
| 9 | What happened on FAIL18? | **Unchanged.** And one of the sprint's own headlines died here: FAIL18 is defined by the filter's own recall, so claims about the filter's effect on it are partly forced arithmetic (§Appendix A) |
| 10 | What happened on the other 108? | **Unchanged** |
| 11 | What information source produced the gain? | **There was no gain.** The sprint's contribution is the reason why |

### 0.2 What the sprint established instead

**The target is reachable. Nothing we own can reach it.** Both halves are measured, and they are
worth more than a fragile 2.98 would have been.

> **Five ORACLE signs on long-range pairs are worth −0.3259 Å on the built chain — 3.2126 → 2.8867,
> 2.05× MDE, 5/5 folds. That clears the charter's primary target.** *(Lane P's arm is measured
> against the S27 projection's 3.2126 rather than the S29 `prod` row's 3.2105; the two differ by
> 0.0021 Å because the projection seed is not pinned — §1.1. The **delta** is what transfers.)* The prize is **five bits per
> target**, concentrated ~4× on the tail and 68% in `|i−j| ≥ 7`.

And the reason those five bits are unobtainable is **mechanistic, not a failed search**:

1. **The source enumeration collapses to two** — sequence and library — with physics an *operator*
   on either, not a third source (theorem G1, repaired by lane P).
2. **Both are measured.** Sequence: out-of-fold R² **0.83%** against a pre-registered 1.96% bar.
   Pool: genuinely informative at **+0.0758**, and **strictly harmful when applied.**
3. **The mechanism is a theorem about *which part* is identifiable.** The predictable part of the
   prior's error *is* the pool's common mode — and the common mode is exactly the component proved
   non-identifiable from pool data at any K.

> ### What can be predicted is coherent and therefore harmful; what would help is incoherent and therefore unpredictable.

A corrector fitted on the distogram and pool emits **+0.0554 Å worse (0.68× MDE — NOT MEASURED)**;
a synthetic i.i.d. one at **identical out-of-fold R² = 0.2355** emits **−0.2466 Å better (1.77×
MDE)**. *A 0.30 Å swing at matched accuracy* — and the sharper statement is that **all three fitted
arms sit at zero-to-worse (0.68×, 0.91×, 0.95× MDE) while both matched-R² i.i.d. arms clear their
MDE in the other direction.** The i.i.d. arms are **ORACLE-constructed and not deployable** — the
price of a channel nobody has, never an achievement. All four numbers are **CA point cloud**
against 3.0483; the measured cloud→chain transfer is 0.92.

### 0.3 The deliverable

Not a number — **a development-time admission test** that replaces an endpoint run with one
correlation, and that **no corrector this project has built would have passed**:

```
coh = corr(a corrector's RESIDUAL, the pool's common-mode pair error)
      -- both arguments require the native (s30_P_lr.py:58,78,202). THIS GATE IS ORACLE.

    uncorrected                               0.6931
    the three fitted correctors                0.7857 / 0.7827 / 0.9172   <- ALL RAISE IT
    imposed-structure ORACLE arms only         0.5865 / 0.5365

    ADMIT iff coh < 0.6931.
```

> **It is ORACLE, and I first wrote "native-free" here. It is not.** `y = expected − d_nat` and
> `mu = pool75_mean − d_nat`; the `d_nat` is printed in the definition in §12.2. **It therefore
> cannot screen a corrector at inference.** What it *can* do is decide, on the 126 labelled dev
> targets, whether a candidate is worth an endpoint experiment at all — one correlation instead of
> a pipeline run — which is the expensive step this project keeps paying.

**What passing it would be worth**, and the label travels with the number: the two matched-R²
arms that clear the gate are **ORACLE-constructed and NOT DEPLOYABLE** — a price for a channel
nobody has, never an achievement — and neither reaches the 2× MDE this sprint used as its survival
standard:

```
R2 0.16 -> -0.126 A   1.23x MDE   [TYPE-M ZONE flagged in the artefact]
R2 0.24 -> -0.247 A   1.77x MDE
                      both CA point cloud; the measured cloud->chain transfer is 0.92
```

**This is not "decorrelated from the distogram"** — orthogonality was priced at a 5.1% discount on
the requirement. Incoherence is a different condition, and it separates a −0.25 Å arm from a
+0.06 Å one at identical accuracy.

### 0.4 How much of this report is negative, and why that is the point

Thirteen lanes, 29 ledger entries, ten pre-registrations. **Four of the ten registered falsifiers
fired against the lane that wrote them.** Directions closed this sprint: the field combination, the
sparse weighted readout, subset objectives through an averaging readout, the second-moment escape
(twice, independently), generative spaces, torsion encodings, common-mode correction, filter width,
achiral single-structure channels, and recognition from single-structure geometry.

**Twenty-six claims were withdrawn in Appendix A's tables — ten of them mine — plus lane G's
self-killed 2.17×-MDE positive, and one withdrawal that was itself later withdrawn.**
The report adversary then found **eighteen** further defects in this document — two of which inverted a
headline — and they are fixed in place with the original error stated. That process is the reason
the surviving results are worth anything, and §Appendix A records all of it.

**The one thing the sprint got wrong about itself:** it claimed, in its own opening ledger entry,
to have saved the charter to disk, and did not — while that same entry was recording four prior
instances of exactly that failure as a reason to check rather than assume. The file was written 77
minutes later. The fix is now mechanical: the verifier asserts that every path the ledger claims to
have written exists.


---

## 1. The arithmetic that governed the sprint

Computed at the open, before any experiment, from S29's stored rows:

```
production, n = 126      mean 3.2105   median 2.9661
the worst 18 targets     mean 6.2758
the other 108            mean 2.6997

cap the worst 10 at 3.00 Å  ->  mean 2.9074  (−0.3031)   beats the primary target
cap the worst 18 at 3.00 Å  ->  mean 2.7426  (−0.4680)
cap the worst 30 at 3.00 Å  ->  mean 2.5778  (−0.6328)

versus improving EVERY ONE of the 126 by 0.20 Å  ->  mean 3.0105  (−0.2000)
```

**Fixing ten targets is worth more than improving all 126 by 0.20 Å.** The endpoint is a mean over
a distribution whose median already clears 3.0, so a mechanism that works on typical targets and
leaves the tail alone is close to worthless here. This is why the sprint staffed the tail first.

The caveat carried with it throughout: capping is an ORACLE operation. It bounds the prize; it does
not deliver it.

---


### 1.1 The endpoint's own reproducibility, which nobody had pinned

Four values for "production built chain" were in circulation this sprint. They are now reconciled,
and the reconciliation is itself a finding:

```
s29/results/s29_O_chain_rows.jsonl :: item=prod      chain 3.2105   cloud 3.0483   <- canonical
s27/results/chain_rows.jsonl       :: config=DIS     chain 3.2126   cloud 3.0483
lane D's independent re-projection                   chain 3.2071
lane P's independent re-projection                   chain 3.2126
```

> **The CA point cloud is identical to four decimals in both records — 3.0483. Only the built chain
> differs.** The projection from cloud to chain is **multi-start and its seed is not pinned**, so
> the same cloud projected twice gives means differing by ~0.002 Å, and independent re-projections
> range to 3.2071. The known cause is documented: a multi-start branch flip on a single target
> (S28-L18/L27b/L43), which lane D reproduced with **median per-target |Δ| of 0.0012 Å and one
> 0.513 Å outlier**.

**This report uses 3.2105** (the S29 `prod` row, which is what the charter quotes). The spread
matters for exactly one reason: **the full range across instruments is 0.0107 Å, and the sprint's
one confirmed deployable effect — the AMBER relax at k=30 — is 0.0221 Å.** That effect is only
**2×** the endpoint's own reconstruction spread. It does not threaten anything at the 3.00 or 2.50
scale, but it means **any future claim below ~0.01 Å on the built chain is inside the noise of the
instrument that measures it**, and the projection seed should be pinned before one is made.

## 2. What was tested

### 2.1 The lanes

Thirteen lanes ran, never more than eight concurrently, against the charter's floor of four
whenever meaningful work was available. Ledger entries `S30-L0` … `S30-L28`, 29 in all.

| lane | remit | entries | outcome |
|---|---|---|---|
| **D** | adversary, permanent; owns and extends the cost/RMSD meter | 7 | Closed the field-combination question; withdrew two published numbers, one of them S29's |
| **L** | literature, permanent | 5 | Three theorems (reference-state compactness, common-mode non-identifiability, set-selection); found the one paper that contradicts our E2 |
| **T** | theory: the bit accounting, and when the CVaR tail stops being a prefix | 3 | **T1** (the tail is always a prefix); the codebook reframe; the value-of-a-bit law. Retracted its own headline |
| **F** | the failure tail — highest leverage by the opening arithmetic | 3 | Tail is **selection-limited, not pool-limited**; shape not scale; closed filter width by ceiling |
| **X** | divergent, permanent: should the quantum stage select at all, or generate? | 2 | The set-mean decomposition — generation closed **jointly with the readout** |
| **R** | is nativeness recognisable from one structure at all (L11) | 2 | Ordering survives, **preference fails on all 43 channels**; the local-feature null is a theorem |
| **Q** | L5 and L6 together: subset objectives and sparse readouts | 2 | Both closed; named the null that made lane T retract |
| **P** | the prior — the one measurement the whole bound reduced to | 2 | The registered null, **and** the coherence result that explains it |
| **G** | the two questions nobody had measured | 1 | **G1**, the chirality dichotomy; confirmed a claim two lanes had talked themselves out of |
| **V** | report adversary (rule 28: the main team may not approve its own positive) | — | Auditing sections 1, 3, 4 and Appendix A against the artefacts |
| **W** | what the CVaR-VQE actually contributed (charter items 12–19) | — | Investigating |
| **Y** | synthesis over the sprint's own record (items 23, 25, 26) | — | §7 and §8; also found seven ledger-vs-ledger contradictions |
| **Z** | second report adversary, on the sections lane V never saw | — | Auditing |
| *(coordinator)* | integration, the report, and a fair share of the errors | 1 | S30-L0 |

### 2.2 Pre-registration

**Ten pre-registration files**, every one committed before the number it predicted existed:
`PREREG_S30_{D_gram, F1, F2, F3, G, P, Q_sparse, R, T, X}.md`. Three lanes' registered falsifiers
then **failed**, and all three reported the failure as the result:

- **Lane D** — its own falsifier for the combination question was not met (S30-L21).
- **Lane T** — retracted its own headline on a null lane Q told it to run (**S30-L15 §3c**).
- **Lane G** — **both** registered priors were directionally wrong, and it said so first (S30-L26).
- **Lane F** — F2a and F2c, both refuted by its own measurements (§Appendix A.5).

*Lane P's P1 is **not** in this list: it **held** (registered bar 1.96%, best arm 0.83%). The
falsified prediction in lane P's work was **mine** — I put 3:1 on long-range R² also coming back
≈ 0 and it is +0.1959.*

Lane P registered its bars at **1.96% / 12.82% / 39.44%** R² before the regression existed and came
in at **0.83%**. Lane G registered a power note stating a half-split *could not* reach its own MDE
unless more than 100% of the effect sat in one half — then measured that more than 100% does.

> **A pre-registration that only ever confirms is decoration.** Four of ten fired against the lane
> that wrote them, which is the evidence that they were real.

### 2.3 The permanent roles, honestly assessed

The charter required three permanent roles. Two were honoured as written; one was not.

- **Adversary — honoured in substance, NOT in the form the charter asked.** The charter says
  *"rotate the adversary role so no agent is only ever attacking its own prior conclusions."*
  **Lane D held the named role for essentially the whole sprint and it was not rotated.** What
  happened instead is that adversarial action distributed itself: lane D attacked lane F's headline
  (S30-L23), lane Q named the null that made lane T retract, lane F closed a mechanism I had
  synthesised one minute earlier, lane G attacked lane L's prior and its own, and lane V was stood
  up at the end specifically to audit *my* report. The charter's purpose was served; its mechanism
  was not followed, and a sprint that relied on distribution happening spontaneously got lucky.
- **Divergent — honoured.** Lane X held it permanently and used it to overturn my premise rather
  than to explore variations.
- **Literature — honoured.** Lane L read for the whole sprint, not only at the start, and returned
  mid-sprint as §6 of the charter asks.

### 2.4 Controls

This sprint's controls are unusually load-bearing — several results died to their own control, and
one control **refuted its lane's stated mechanism while strengthening its conclusion**. The full
register of every control run and what each one ruled out is **§7.3**, compiled from the ledger.


---

## 3. What was closed, and by what

The sprint's dominant output. Each row is a direction that is now shut, with the instrument that
shut it. Several are theorems rather than measurements.

| direction | closed by | how |
|---|---|---|
| combining the 21 native-free displacement fields | D + T | ORACLE global ρ = **0.169**, LFO **0.012** (worse than the best single field) against **0.358** needed for 3.00 Å; Gram stable rank 2.057 — the fields span ~2 directions, not 11 |
| sparse weighted readout (S29's last open ladder class) | Q | **by price** — a plain argmin dominates at every bit budget; 2-of-75 with free weights is 2.1683 Å at 11.4+ bits against argmin-over-128 at **2.1435 Å for 7.0 bits** |
| second-moment / quadric escape | Q **and** T, independently | +0.2059 Å (Q, 1.81× MDE) and +0.086 (T, 2.95×); the structured `disp2` form scores 3.3585 against production's 3.0483 |
| subset objective through an averaging readout | T | **T1**: the tail is *always* a prefix — of the order induced by ∇V at the optimum. S29 §4.3 and §4.5 are incompatible and §4.3 wins |
| generative structural spaces | X | closed **jointly with the readout**: at ρ≈0 the endpoint is near-unit-slope in the set **mean**, and width buys ceiling while costing mean |
| torsion / configuration encodings | T | arithmetically infeasible — **48 bits** (2 bits/torsion x 2 torsions/residue x 12 residues) against **7** deployed |
| common-mode correction from pool data | L | **non-identifiable**: the likelihood depends on (t, μ) only through t+μ at any K; **m_eff = 1.4** of 75 members |
| E1 — a prior on the shared bias's form | L | three independent ways |
| E2 — constraint repair — **NOT CLOSED; listed here only because its *field-scale* form is** (see §7.4) | L, **revised by G** | field-scale: removing 98% of clashes costs **+0.08 Å**; works here only because 2/n is 15.4% at n=13. **G then CONFIRMED it is real and divergence-graded** (−0.0406, 3.56× MDE, 5/5 folds length-matched; 97.2% of the gain in the high-dispersion half; divergence primary, contraction its shadow at ρ 0.947). **It stays closed on SIZE — 0.69% of baseline, and its restraint constant has no native-free selection rule** (§A.7) |
| filter width as a free lunch | F | no knee — benefit/harm cross smoothly at ~55%/55%; the **ORACLE global argmin over k IS the shipped 75** |
| recognising nativeness from single-structure geometry | R | ordering survives on 2 of 43 channels; **preference fails on all 43**. And the mechanism is a **theorem on this instrument**: local features give ΔR² **−0.089** against global **+0.600**, so **a sum of per-residue terms cannot see a lever arm** |
| **every achiral single-structure channel** | G | **by theorem (G1).** Classical MDS: a rotation/translation-invariant channel is a function of the distance map **iff** it is reflection-invariant — so contact topology, Rg profile, inertia/shape tensors, separation distributions, contact order, burial/SASA are all *distance-map re-readings*, not new channels. Sequence enters as a fixed parameter and does not affect the argument |
| the one family that theorem leaves open — **chiral functionals** | G | **built and empty, and not degenerate.** The pre-check refuted the easy explanation: chiral-axis occupancy is **0.685–0.965 against the shipped cost's 0.345** — *more* exercised than DIS. Best channel's anchor contrast **+0.0405 against a max-over-3 null mean of 0.0408** (p = 0.430); `WRITHE − \|WRITHE\|` is 0.29× MDE with the CI including zero. **Length-scoped:** the theorem is length-free, the emptiness is not |
| correcting the prior from native-free features | P | **the registered null** — out-of-fold R² **0.83%** against a pre-registered 1.96% bar, and applying the fitted displacement emits **3.0519 against production's 3.0483**, i.e. marginally *worse* |
| correcting the prior at all, from {sequence, library} | P | **by theorem about which part is identifiable.** The predictable component **is** the pool's common mode, which is non-identifiable from pool data at any K. At identical out-of-fold R² = 0.2355 a fitted corrector emits **+0.0554 worse** and a synthetic i.i.d. one **−0.2466 better** — *a 0.30 Å swing at matched accuracy* |

---

## 4. What was established positively

### 4.1 The tail is selection-limited, not pool-limited

The claim holds on **two** tail definitions, and the stratum label matters, so both are given.

**On `FAIL18`** — the *filter-defined* stratum (`s12/instrument.py:271-278`: the targets whose
score-top-75 retained **zero** members within 1.5 Å of the pool optimum) — the ORACLE best pool
member is **2.2842 Å** (built chain 2.2845), with **13 of the 18** holding a member under 3.00 and
the worst bottoming out at **3.5436**.

**On the genuinely worst 18** (by production built-chain RMSD, an order statistic, overlap with
FAIL18 **13 of 18**) the same figures are **2.5298 Å**, **11 of 18** under 3.00, worst **3.9523**.

Either way the direction survives and it is the point: **2.53 < 3.00, so the material to fix the
tail is already inside the candidate sets the pipeline is handed.** The larger claim I first wrote
— that *nothing* here is conditioned on what it measures — was too strong: the 2.2842 figure is
computed on a stratum defined by the filter's own recall, and only the worst-18 figure is free of
that.

**Retrieval is NOT MEASURED at any stratum — which is not the same as exonerated.** BLOSUM's 500
against a random 500 of the same universe fails to clear its MDE everywhere, most pointedly on the
tail (0.01× its own MDE). But the measurement is **ORACLE best-of-pool on the CA cloud** with an
MDE of 0.29 Å on the tail, power 0.05 and type-M 71, whereas the record it appears to revise
(blind 5.425 against shipped 6.019 on FAIL18) is about the **emitted built-chain endpoint**. *A
different operator, on a different basis, at a power that could not have detected the effect.*
**The prior record stands; this adds an absence of evidence, not evidence of absence** — and it is
this sprint's one instance of the control-space error (§6.1).

**The tail ordering result, with the three strata separated rather than merged:** the shipped
score's Spearman with ORACLE in-pool RMSD is **+0.6446 on the easy 108** (fold CI [+0.5905,
+0.7063], 5/5 folds) against

```
FAIL18 (filter-defined)   +0.1066   fold CI includes zero, 3/5 folds
worst18 by pool mean      +0.3798   no fold CI stored
worst18 by best-in-pool   +0.3438   no fold CI stored
random-18 null, lower bound         +0.4067
```

All three sit below the random-18 null's lower bound, so the degradation is real — but the
filter-defined number is **3.2–3.6× lower than the other two**, and only it has a CI spanning zero.
The honest statement is therefore **"the score orders its own pool markedly worse on hard
targets,"** not that it cannot order it at all.

It is the distogram's own error that drives this (ρ = −0.799, −0.819 length-residualised) and the
chain is fully mediated. On shape versus scale, lane F attached a quotation condition to its own
row and it travels here verbatim: **shape error is 83.1% of it, and `|scale error|` is 3.23× on
FAIL18 against shape's 1.93× — both facts belong in any quotation of this row.** Scale is refuted
as the *mechanism* (partial ρ = −0.067, p = 0.46, against shape's −0.641 at p = 6.5e−16) while
being the larger *absolute* deviation. Lane F also records `F3c.fires = False` with the caveat
*"statistic guessed by lane F; lane L was unreachable"* — an unregistered self-chosen statistic,
which is why this is stated as a refuted mechanism and not as a headline.

### 4.2 The pool is a codebook, not a channel

The charter asked where the 1.44 usable bits went and where the other 5.56 were spent. The question
is not well-posed: **the 500 deposited backbones carry the structure and the index only names it**,
so bits are not conserved across an index. Seven index bits move the ORACLE ladder 4.108 → 1.898 Å
**on the CA point cloud**,
which through the displacement bound is ρ = 0.887.

**The multiplier that follows is a free parameter and this sprint corrected it downward, so it is
reported with its formula rather than as a headline.** The information is
`I = −(d/2)·log₂(1 − ρ²)`, and everything depends on `d`:

```
d = 3n - 6 = 32.88  (the full coordinate dimension)   ->  36.6 bits, a 5.2x ratio
d = 6               (the effective dimension)         ->   6.69 bits, a 0.95x ratio
```

**At the corrected `d` the multiplier inverts to slightly under one.** S30-L14 §6 says so in lane
T's own words — *"S30-L15 corrects the `d` this should be evaluated at, and the corrected numbers
are smaller"* — and Appendix A records the companion 12 → 3.78 bit withdrawal made on exactly this
ground. **I originally published the `d = 32.88` figure with neither the formula nor the `d`.**

> What survives, and it does not depend on `d` at all: **bits are not conserved across an index,
> because the 500 deposited backbones carry the structure and the index only names it.** The
> *qualitative* claim is robust; the multiplier is not, and the system still cannot supply the
> seven bits either way.

The **value-of-a-bit law** follows and makes allocations comparable: `D(R) = a + c·2^(−R/γ)` fits at
**R² = 0.9983** (a = 1.3312 Å, γ = 3.1636), so `−dD/dR = 0.219·(D − 1.331)` Å per bit. Candidate
indexing beats subset cardinality by **3×** — which *explains* the 3.5× S29 measured — and
**torsion/configuration encodings are arithmetically infeasible** at this width (**2 bits per
torsion x 2 torsions per residue x 12 residues = 48 bits**, against 7 deployed; the ledger's
per-torsion figure is 51.8 bits). *As first worded — "4 basins per residue" — the arithmetic gives
24, not 48; the factor of two is the two backbone torsions.* The charter called the encoding the least-examined
component and the likely hidden bottleneck. **It was examined and it is not.**

### 4.3 The field library spends its rank on the wrong direction

Lane L derived that a fixed-reference statistical potential contains a **separable term that is a
pure function of scale** (`+kT·Σ ln P_ref(d_ij)`; a uniform 10% contraction of a 13-mer with zero
shape change moves it ~14 kT, while a size-matched reference moves by exactly 0.0000). It predicted,
falsifiably and with its concession pre-stated, that **one of the Gram's two effective directions
must therefore be the radial one**. Measured independently:

```
radial share of the Gram trace                        0.5798  [+0.545, +0.603]
cos(dominant principal direction, radial)             0.947 mean / 0.984 median
radial share of lambda_1 (ratio of trace shares, NOT a projection)   91.2%
stable rank with radial removed                       1.705 -> 2.642
   (ALL of the above are PAIR-DISTANCE space; coordinate space is 3.4-3.6 -- A.1)

cos(direction to the native, radial), all targets    -0.0675
cos(direction to the native, radial), FAIL18         -0.2524
```

> **The library spends the majority of its two available directions on a component that is
> orthogonal to the answer** — and lane F had independently shown that shape, not scale, is 83.1%
> of the tail's error. **Two lanes, one chain.**

**The anti-alignment half is withdrawn.** The FAIL18 figure of −0.2524 is **0.67× its own MDE**
(n = 18, SE 0.1337) — *below this sprint's own "not a result" line* — has no CI anywhere in the
artefact, and one of four non-empty folds (−0.6491) carries it. I bolded it in the first draft. The
orthogonal-in-general half stands, because at 0.74× MDE it is claiming a **null**, which is what a
sub-MDE number can legitimately support.

It does **not** follow that deflating scale helps: the residual 42% has no large eigenvalue, the
combination's ceiling already prices the whole span, and lane L withdrew its own deployable proposal
on exactly this point — the size-matched field *is* the deflated field, and deflation reallocates
rank without creating any.

### 4.4 Ordering exists; preference does not

On a ladder where **kind, local realism and perturbation budget are all matched** — every rung an
ideal-geometry backbone built from torsions drawn from the fold's leakage-safe Ramachandran table,
no projection, no averaging, no contraction — the verdict splits:

- **Ordering survives.** Two of 43 channels clear the bar, DIS at +0.347 with an anchor contrast of
  +0.134 (max-over-channels sign-flip null p = 0.000).
- **Preference fails on all 43.** The best `pref_near` in the library is 0.640, under the bar, and
  that channel prefers a *random pool member* to production on 0.790.
- **The largest preference effect in the library points the wrong way:** DIS prefers production to a
  0.55 Å structure on **94.2%** of targets. The project's most-cited negative, with its confound
  removed, comes out **sharper**.
- The leave-fold-out combination prefers the near-native rung to production on **93.0%** of held-out
  targets — and prefers an **arbitrary pool member on 100%** and a **3 Å rung on 100%**. Margin
  **−0.070 [−0.110, −0.028]** — *1.04× MDE at 4/5 folds, so it only just clears its own bar.*
  *It learned "is this production?", not "is this near-native".*

**And the mechanism makes the null a theorem on this instrument.** Held-out R² at matched capacity:
local features give **ΔR² −0.089**, global features **+0.600** — and the local block is
**ORACLE-advantaged**, handed per-residue deviations from the native anchor, and still adds nothing
beyond knowing the perturbation budget. **A sum of per-residue terms cannot see a lever arm.**

Resolution, which bounds what any support rule could ever select on: DIS concordance inside the near
band is **0.520** at |Δ| = 0–0.25 Å; concordance reaches 0.822 only above 4 Å, **and that figure
comes from the all-pairs table rather than the near band**. **Coarse triage and nothing
else.**

### 4.5 Generation is closed jointly with the readout

For a coordinate-average terminal, `set_mean² ≈ B² + S²` where **B is the endpoint itself** and S is
the set's spread. The terminal's entire value is spread extraction (`avg_gain = 0.4143·S − 0.1559`,
r = 0.885), and within target `corr(S, B) = +0.085` **after dropping the `T0_helix` arm** — the
artefact's primary value over all arms is **−0.4744**, which lane X itself flagged as driven by that
one arm and "must not be quoted". Both belong here: on the arms that matter, **concentration is
orthogonal to bias**. So a
set mean improved purely by concentration is worth **zero by algebra**, and the half that pays *is*
the endpoint, which is non-identifiable from pool data at any K.

The extreme case settles it: the most concentrated source ever built here (spread 0.568 against the
pool's 1.577) is **the worst endpoint in the record, 3.789 on the built chain**. **For an averaging terminal, spread is
the raw material, not a defect.**

---

## 5. What the CVaR-VQE contributed

*(Charter items 12–19. Audited by lane W against the code and the artefacts, not against the
sprint's own prose; full working in `s30/QUANTUM_W.md`.)*

### 5.1 The honest summary

Sprint 30 ran **no CVaR-VQE**. It *did* run a 9-qubit depth-3 statevector circuit on all 126
targets with exact parameter-shift Jacobians and 300 Adam iterations — but only to regenerate the
meter's ORACLE rungs (`s30/results/s30_D_ladder_structs/`, 126 files, asserted to 1e-6 against
S28's stored values). **"No VQE" and "no quantum compute" are different claims and only the first
is true**; a report sentence saying the sprint spent no quantum compute would be false as written.

The production anchor does not pass through the quantum stage at all: `core/pipeline.py:241`
(`PROD = Config()`) with `core/pipeline.py:179` (`quantum: bool = False`), and the governing
comment at `:173-178` says it outright — *"AMBER participates, Legacy is only imported for scoring,
and VQE/CVaR do not participate at all."*

What the sprint contributed to the quantum spine is **three closures and a theorem**: **T1**, which
shows the CVaR tail is always a prefix of the order induced by `∇V` at the optimum, so the deployed
stage's structural output is the single integer `m`; **T1b plus two independent quadric
measurements**, which price the only two Hamiltonian classes that could have changed that order and
find both reachable but non-transferable — **closed before either was built**, the first time this
project's cheap pre-check fired *ahead* of the spend rather than after; and **the bit accounting**,
which shows the encoding is not where the loss is.

The reason none of it was carried into a circuit is **S30-L19**: on a kind-matched, budget-matched
ladder, ordering survives on 2 of 43 native-free channels and **preference fails on all 43**, with a
mechanism (local ΔR² −0.089 against global +0.600) that makes the null a theorem on this instrument.

> **An objective cannot be pointed at a target nothing can see.**

### 5.2 The single row that answers items 12 and 16 at once

From this sprint's own cache, means over n = 126:

```
rung                                                    built chain      CA cloud
circ_best   same circuit, ORACLE objective, best of 5       0.2516        0.2884
circ_s0     ORACLE objective, 1 start (regenerated S30)     0.3175        0.3854
PROD        the deployed uniform average                    3.2071        3.0483
circ_opt    THE SAME CIRCUIT, DEPLOYED native-free score    3.4330        3.3850
```

*(`PROD` here is the meter's own re-projection; see §1.1 on the ±0.002 Å projection spread.)*

> **The same circuit reaches 0.2516 Å when its objective is the native and 3.4330 Å when its
> objective is the shipped score — worse than the classical average it was meant to improve.** The
> circuit is not the problem. The objective is, and it has been all along.

### 5.3 Item by item

| # | item | answer |
|---|---|---|
| 12 | What did CVaR-VQE contribute? | **Nothing as an endpoint contribution.** Its residual role, from T1: **α supplies τ; it is not doing selection.** Lane L's finite-shot CVaR bias hypothesis was tested and **refuted, in the opposite direction to its guess** (S30-L8) |
| 13 | What Hamiltonian was used? | **NONE.** Charter §14 withdrew the obligation. Two candidate classes (T1b halfspace; second-moment quadric) were priced and **closed before either was built** — S30-L12 and S30-L15 independently, two samplers, same sign |
| 14 | Why different from the diagonal rank-ladder? | **No new Hamiltonian to be different** — and S30 established the rank-ladder's diagonality **is not the defect**. `s30/THEORY.md:450-453` removes Hamiltonian structure, ansatz, optimiser, pool and search from the list. Referent: `H = diag(zrank(top-128 scores))`, `core/pipeline.py:838`, the same ladder on every target to **1.18% of range** |
| 15 | What did the quantum state represent? | **Partly applicable.** T1's trichotomy: under the deployed formulation the state determines **exactly one integer, `m`**. The codebook result says the 7-bit register realises 36.6 bits at `d = 32.88` (0.95× at the corrected `d = 6`; §4.2) and **delivers under one** |
| 16 | Could a classical control reproduce it? | **The honest form is "is there any positive for a classical control to reproduce?" and the answer is NO.** Every decisive classical control runs *against* the quantum stage: one sort reproduces the tail exactly; Frank–Wolfe is the correct counterpart and is polynomial; argmin dominates the sparse readout at every budget; `PROD` beats `circ_opt` |
| 17 | Trainability | **NOT APPLICABLE, not measured this sprint.** S28/S29 positions cited in `QUANTUM_W.md` rather than restated |
| 18 | Gradient variance | **NOT APPLICABLE, not measured** — but S30 measured the governing input of `Var[∂F/∂θ] = r_stable/D²`, namely `r_stable`, **on 126 real pools in 76 seconds** over data already on disk. Carry lane Q's scope fix: **`r_stable = 1.859` is pair-distance space; coordinate space is 3.404** and does not fire the threshold |
| 19 | DLA / ansatz structure | **NOT APPLICABLE, not measured.** S29's scope correction stands: the algebra is **not** maximal at n = 9, L = 3 |

### 5.4 A defect in shipped code, found while answering item 12

`core/pipeline.py:821` — the `quantum_stage` docstring — still asserts *"the CVaR tail is worth
+0.113 Å by preventing the collapse, and that is the component's measured role."* **S25-L5 withdrew
that number** (`s25/LEDGER.md:233-244`, with the replacement figure at `:280`: *"never a measured effect … by
this project's own fixed rule that is a NULL"*) and replaced it with **−0.1405 Å at 0.68× MDE**.

> Anyone answering item 12 from the source file rather than the ledger gets a **withdrawn positive**
> presented as the component's measured role. This is the **sixth** instance of prose asserting a
> state that does not hold, and **the first located in shipped code rather than in a report.** It
> is not fixed here — `core/` was read-only for this audit — and it is listed in §11 as open.


## 6. Statistical discipline and multiplicity

### 6.1 The rules actually enforced

| rule | where it bit |
|---|---|
| **MDE = 2.8016 x SE, per comparison** - never a sprint-wide constant | the S29 constant 0.084 A is wrong by up to 84x in both directions |
| **Below 0.7x MDE is not a result; 0.7-1.0x is NOT MEASURED** | killed S29's built-chain preference row at 0.56x (§6.3) and lane P's headline at 0.75x |
| **Fold-clustered CI, on the pinned folds** | `s24.stats_lib._verdict` refuses a verdict when `folds=None` rather than falling back to the IID CI |
| **Pre-register the falsifier before the number exists** | every lane did; three lanes' falsifiers then failed and were reported as failures |
| **A control must match the operator's own space** | the project's most repeated error, and **there WAS an instance this sprint**: the "retrieval is exonerated" claim in §4.1 compared an ORACLE best-of-pool CA-cloud measurement against a record about the emitted built-chain endpoint. Caught by the report adversary, not by the lane |
| **The native never tunes a deployable parameter** | every ORACLE arm is labelled ORACLE and none is deployable |

### 6.2 Multiplicity, counted rather than asserted

The meter counts every comparison it emits - **30 on the built chain, 33 on the CA cloud** - and
states the consequence in the artefact itself: *"With 8 lanes, a 1x-MDE positive is EXPECTED
somewhere; price any single positive from this instrument against this count and the sprint's
running total."*

Sprint-wide the count is in the **hundreds**: 43 channels (lane R) + 21 fields (lane D) + 8 feature
arms x 5 strata (lane P) + the width sweep, the quadric/halfspace ladder, the bit ladder, and the
63 above. This is why **the surviving positives are the ones at 2-10x MDE with 5/5 folds**, not the
ones at 1.2x. Two lanes ran max-over-channels nulls rather than per-channel ones (lane R's
sign-flip null across 43, p = 0.000; lane Q's across-target null on the search ladder) and that is
the standard the next sprint should inherit for any swept family.

**The honest limitation:** the sprint-wide total is a *lower bound*. Lanes counted their own
comparisons; nobody maintained a single register, and by the time that was obvious the lanes had
closed. **A sprint-wide multiplicity register, written to as comparisons are emitted, is a build
item for S31** - not a discipline problem to exhort about.

### 6.3 Three statistic/null mismatches, which is the sprint's real methodological finding

Three times a statistic was compared against a null belonging to a *different* statistic. They are
the same defect wearing three costumes:

1. **`7 - log2 r` is right-skewed**, so its median sits 0.416 bits below its own mean. Reading the
   median against the analytic mean says "worse than chance" **by construction**. (Mine, reported
   to the user twice.)
2. **`abs(signed mean)` scored against a null for `mean |cos|`** - `s29/s29_D_fields.py:243,263`
   labels one statistic and computes the other.
3. **A 0/0.5/1 preference indicator has modal value 0**, so its median is uninformative by
   construction and a median paired difference of 0.0000 is *not* the median-vs-mean warning firing.

And a fourth, latent, caught by the verifier rather than by a reader: **`s24.stats_lib.compare` is
lower-is-better** (`d = a - b`, negative = a better) because its native statistic is RMSD. Fed a
**preference rate**, which is higher-is-better, its `verdict` string reads `WORSE` for an effect of
**+0.1716 that is the good direction**. Lane D's gate handles it correctly
(`s30_D_meter.py:451`, `PASS if effect > 0`) and no S30 document quotes the inverted field -
`s30/s30_verify.py` now asserts the trap's shape so it stays known.

### 6.4 The control that deserves more credit than it got

Lane G's **chirality occupancy pre-check** is the only control this sprint that **refuted its own
lane's registered mechanism and strengthened the conclusion by doing so.** Lane G predicted the
chiral coordinate would be near-degenerate at peptide length; occupancy came back at **0.685–0.965
against the shipped cost's 0.231–0.345** — *more* exercised, not less. Had it come back degenerate,
the null would have been uninterpretable.

> **A pre-check that can only remove your own excuse is worth running before the contrast, not
> after.** It costs one pass and it converts a quiet null into a real negative.

---

> **The generalisation, and it is the one to carry forward:** a statistics library encodes a
> *direction* as well as a test. Every quantity handed to it - preference rates, accuracies,
> correlations, R2, win rates, concordances - must be checked against that direction, because the
> failure is silent and the field that is easiest to quote is the one that is wrong.


## 7. Every hypothesis entertained and killed

Charter item 26. Ninety-one hypotheses were entertained across ten lanes. They died three
different deaths, and the distinction is the most reusable thing in the table: a hypothesis
**killed by measurement** can be reopened by a better instrument; one **killed by theorem** cannot
be reopened at all within the model class named; one **killed by price** is not wrong — it is
correct and unaffordable, and reopening it requires a change to the budget, not to the biology.

The sprint's own distribution over ninety-one: **58 killed by measurement, 19 by theorem, 8 by
price**, with **6 entertained and still open** (§7.4). Roughly half of the 58 were killed by the
lane that raised them, usually before anyone asked, and **every one of the ten lanes killed at
least one of its own**.

### 7.1 Killed by measurement

Ordered by lane. "predicted" is what the lane or the coordinator said would happen, taken from the
pre-registration where one exists.

| # | hypothesis | lane | what was predicted | what killed it | effect / ×MDE |
|---|---|---|---|---|---|
| 1 | The failure tail is **pool-limited** — the candidate sets do not contain good structures | F | F1a fires only if ORACLE best-of-pool on FAIL18 is **below 3.00 Å** (`PREREG_S30_F1.md`) | Direct ORACLE measurement of the 18 pools | **2.2842 Å** point cloud (2.2845 chain), 13/18 under 3.00, worst pool 3.54 (S30-L2 §2) — one-sided threshold, no MDE by design |
| 2 | **BLOSUM retrieval** is what harms hard targets (project memory `sequence-conditioning-hurts-the-failures`) | F | Retrieval degrades the tail | BLOSUM-500 vs a random 500 of the same universe, 2,000 draws | FAIL18 **−0.0034, 0.01× MDE**; all-126 −0.0722, 0.87×; other-108 −0.0837, 0.99× — **NOT MEASURED everywhere** (S30-L2 §3) |
| 3 | A **wider filter is a free lunch** — keep the body's benefit, shed the tail's harm (the coordinator's hypothesis, registered as F2a) | F | Some k>75 retains **≥80%** of the 108's set-mean benefit while shedding **≥50%** of FAIL18's harm | The two curves measured across k ∈ {75…500} | No k fires. Benefit/harm cross **smoothly at ~55%/55% near k≈275**; no knee (S30-L16 §2) |
| 4 | **Widening the filter rescues hard targets** | F | FAIL18 improves monotonically to k=400 (−0.6193, random-18 null p=0) | Its own filter-independent-tail control | Does **not replicate**: ≈0 on worst-18-by-ORACLE-best, and **reversed at +0.7313** on worst-18-by-pool-mean. One mechanism (ρ = **+0.81** with the filter's own set-mean benefit) explains all four strata (S30-L16 §4) |
| 5 | FAIL18 is simply where production was **unluckiest** (regression to the mean) | F | A linear RTM control explains the k=400 gain | Regression of the per-target effect on production RMSD over the 108, extrapolated | Predicted **+0.069**, observed **−0.619**, residual −0.688 — RTM does *not* explain it (S30-L16 §4) |
| 6 | **Lane L's scale-separable term** is the mechanism behind the filter's collapse (F3c) | F | Partial Spearman(ρ_pool, scale error \| shape) **≤ −0.30** with fold CI excluding zero | Shape/scale decomposition of the distogram's error | **−0.067, p = 0.46** against shape's **−0.641, p = 6.5e−16**. Scale error *is* elevated 3.23× on the tail and is still not the mechanism; shape is 83.1% of the tail's error (S30-L17 §3). **Provisional** — lane F guessed the statistic from a second-hand relay and said so |
| 7 | The distogram is **confidently** wrong on hard pools, not merely wrong (F3d — the coordinator's framing, adopted by lane F) | F | error×confidence beats error alone by **≥ 0.10 of Spearman** | Both predictors of ρ_pool measured | error×confidence **−0.729**, *worse* than error alone **−0.799**; gain **−0.069** against a +0.10 bar (S30-L17 §4) |
| 8 | **Pool Rg dispersion** predicts where the filter fails | F | A compactness-disagreement route | The F4 pool-Rg-sd column | **−0.128, fold CI [−0.316, +0.036] — NOT MEASURED** (S30-L17 §2) |
| 9 | "The median target is **worse than chance**" under the deployed score | D (correcting the coordinator) | S29 §11's supporting paragraph, reported to the user twice | A simulated uniform-ranking null for the **exact** statistic, 20,000 draws, plus KS and χ² | All five statistics at chance, **p = 0.16 to 0.75**; KS D=0.0685 p=0.571; χ² 8.35 p=0.303. "82 of 126 below 1.405" is what a random ranking does to itself (null expects **78.8, z = +0.60**) (S30-L4) |
| 10 | "**No field's mean \|cos\| clears the random reference 0.140**" — S29's printed verdict, the basis of assumption B2's argument | D | The 21-field class is empty | Reading `s29/s29_D_fields.py:243,263` and recomputing both nulls | The named quantity clears 0.1398 on **21 of 21** (0.2505–0.3249); against the correct signed null (+0.0014 ± 0.0144) **11 of 21** have fold CIs excluding zero, best **+7.7 σ**. `any_beats_reference` was `[]`; the correct count is 11 (S30-L5) |
| 11 | The 21 fields **combine** past ρ = 0.358 (the coordinator's proposal after S30-L5) | D | √(Σρᵢ²) ≈ 0.33–0.36 if the fields are near-orthogonal | The Gram, then the combination itself | ORACLE **global** ρ = **0.1693** (0.69 bits vs 3.22 needed, **0.046 Å** on the chain); leave-fold-out **0.0948**, and the 21-parameter fitted arm **0.0124** — all **below 0.7× MDE**. Best honest arm vs best single field: **−0.0267, −0.36× MDE** (S30-L21) |
| 12 | The ORACLE **per-target** combination's ρ = 0.9491 is a ceiling of the field class | D | 0.949 → 1.011 Å | Its own matched random 21-dimensional subspace control | Control reaches **0.8095** by dimension counting. Excess **+0.1396, 6.83× MDE, 5/5 folds** is real; the raw 0.949 is not a ceiling. *"Without that control this would have gone into the record as a 0.97 Å result"* (S30-L21) |
| 13 | **Deflating the radial/scale direction** would create rank in the field library | D (against lane L's constructive half) | A size-matched reference builds a *different* field | Lane L's own algebra, checked by lane D's spectrum | The score decomposes additively, so **the size-matched field IS the deflated field**; residual 42% has no large eigenvalue and the whole span is already priced at ρ = 0.1693 (S30-L22; conceded in S30-L18 §8) |
| 14 | S29's **built-chain preference contrast** (+0.0635, "f prefers the 0.29 Å ORACLE structure to production more often than a matched random displacement") | D | 0.92× MDE — a borderline positive | Eight matched random-signed draws instead of one | **+0.0357, 0.56× MDE**, fold CI [−0.011, +0.082], folds 3/5. **S29's seed-0 draw was the maximum of its own eight** (range [+0.0159, +0.0635]). Relative draw noise 59% on the chain (S30-L24) |
| 15 | Lane F's **F1c FAIL18 row, +1.7674 Å** — "the harmful thing the pipeline does on hard targets is located in one stage" | F, killed by D | 3.87× MDE, 0W/18L, random-18 null p = 0 | The adversary read the definition: `s12/instrument.py:271-278` defines FAIL18 as the targets whose top-75 retained **zero** in-band members | Reproduced at −1.7622 (opposite sign convention). **≥ 0.869 Å (49%) is forced arithmetic**; `keep = 1` shows **−0.0363**; outside the 18 the filter is **+0.0245 Å, better** than random; all-126 effect is **100% the 18** (S30-L23) |
| 16 | The instrument (`s29_D_cost_audit.py`) could gate new costs on the reporting basis | D | `--basis chain` works | Cold run | **0 of 126 ladder-cache files carried chain projections**; `chain-s28rows` serves only the 31 names S28 stored. The two bases differ in the **sign** of the CHARTER ladder ρ (+0.2603 CA, −0.0921 chain) (S30-L3 §3) |
| 17 | "Run `selftest` and confirm it reproduces the baselines" (the coordinator's instruction to lane D) | D | A reproduction check | Reading the code | `selftest` is a **synthetic 8-residue check** running in 0.4 s that would pass on a meter whose every dev-set number had drifted. The instruction was unsatisfiable; a real `verify` now exists (S30-L3 §2a; 22/22 matched, 0 mismatched, 0 missing) |
| 18 | **F-R1** — nativeness is recognisable from single-structure geometry, on a kind-matched ladder | R | Fires only if some channel clears **both** (i) anchor contrast ≥ +0.10 with fold CI excluding zero **and** (ii) pref(≤1 Å rung vs production) ≥ 0.65 with margin ≥ +0.10 over the pool-member control | 43 channels on an ideal-geometry, budget-matched, leakage-safe ladder | **Clause (i) fires on 2/43** (DIS +0.347, contrast +0.134, p_max 0.000). **Clause (ii) fires on nothing** — best `pref_near` in the library is RAMA at **0.640**, under the bar, and RAMA prefers a random pool member on **0.790**. F-R1 **does not fire** (S30-L19) |
| 19 | The **leave-fold-out combination over 40 channels recognises nativeness** — it prefers a 0.55 Å structure to production on 93.0% of held-out targets | R | A headline | Its own pool-member and 3 Å-rung controls | Prefers an **arbitrary pool member on 100%** and a **3 Å rung on 100%**. Margin **−0.070 [−0.110, −0.028], 1.04× MDE, 4/5 folds, n = 114**. *It learned "is this production?"* (S30-L19) |
| 20 | **LEG_torsion** — S29 §12.0's last outside-class-M hope | R | In-band skill +0.181 outside the bounded scorer class | The matched ladder | Anchor contrast **+0.024, fold CI [+0.010, +0.038]** — a **quarter** of the registered +0.10 margin; pref 0.503 against a 0.698 pool-member control; ranks the rebuilt native at the **exact median (0.503)** of the native's own perturbations (S30-L19) |
| 21 | **Per-residue / local channels** can see nativeness | R | RAMA, LEG_torsion, DSSPHB, CAGEO, HP | Held-out R² at matched capacity (84 local vs 80 global features), with the local block **ORACLE-advantaged** | local **ΔR² −0.089** [−0.122, −0.052]; global **+0.600** [+0.585, +0.617]. It only costs capacity (S30-L19 D1) — this is the measurement that makes the class's null a theorem (§7.2 #8) |
| 22 | Lane R's own **0.500/0.500/0.500** combination cell | R | Three answers to three questions | Three identical values across three different questions | A **null-input artefact**: the channel filter required `d_near` on every row and 3 targets have no rung under 1 Å, so the channel set was empty (S30-L19, disclosure) |
| 23 | Removing the **forced scale term** decorrelates a channel from Rg (lane L's construction, tested by lane R) | R | ρ(score, Rg) → 0 | Size-invariant twins built and audited (RG_LAW relative sd 0.5009 → **3.2e−16**) | Several twins are **more** Rg-correlated: CONS −0.374 → **−0.807**, DMAP_CONS −0.050 → **−0.729**, POOLGO +0.611 → **−0.678**. At peptide length shape and size are genuinely correlated; scale-invariance by construction is a **correctness fix, not a decorrelation** (S30-L19) |
| 24 | The **second-moment (quadric) escape** beats the halfspace class it was derived to escape into | Q, then T independently | T1b: VC dim 34 → ~595, bit cap 50.4 → 152.1 | Enumeration on the real pools with the readout held fixed | Q: **−0.0105, 0.16× MDE (NOT MEASURED)** at ORACLE m; **+0.2059, 1.81× MDE, 37W/89L (WORSE)** at the shipped m=75 (S30-L12). T at matched budget K=5,000: **+0.0860, 2.95× MDE, 5/5 folds, 21W/105L** (S30-L15 M4). Two samplers, two lanes, same sign |
| 25 | The **dispersion rule** the second-moment construction points at (S29's own post-mortem: "operators that read the pool's own dispersion") | Q | The one structured form the derivation specifies | `disp2(x) = ‖W_x − mean(W)‖²` scored directly | **3.3585 Å at m=75 against production's 3.0483** — 0.31 Å *worse* than what ships. Every native-free rule direction tested (rg, dist_to_medoid, PC1–PC3) is worse at both set sizes (S30-L12) |
| 26 | The **halfspace/quadric ceiling** (1.735/1.725 Å) is a property of the class | Q | A reachable ceiling | `best_of_k_within` against an across-target null + split-half transfer | LINEAR: **225% accounted** by the null, split-half **+0.0323 (−3%)**; QUADRIC 203%, **−0.0204 (2%)**. The null *exceeds* the observed gain in both (S30-L12) |
| 27 | The cut classes would degenerate to m = 1, making them a relabelling of the argmin | Q | Its own worry | Best prefix size per direction | **Median m = 402**; 81% of directions peak at m ≥ 50, only 4% at m ≤ 2. The class genuinely selects large sets — *"it simply cannot select them from outside the target"* (S30-L12) |
| 28 | A **native-free support rule** for a sparse readout (F3) | Q | It would be the contested question | Best of {score-prefix, farthest-point diversity} vs the mean of 8 random supports | Confirmed **and worthless**: +0.364 Å over random at s=2, and **every** native-free support arm is worse than production (3.198 vs 3.048) even with ORACLE continuous weights handed to it free (S30-L11 F3) |
| 29 | **"Measure the ceiling first, always"** — the admission gate the coordinator handed lane X | X | A gate on generative spaces | Scored against the endpoint's actual direction on 10 cells, two independent families | **1 of 10** sign-correct against **5/6** for the frozen source law; mean \|residual\| **0.333 Å vs 0.0153 — 22×**. On the S7-6 ladder the ceiling improves 0.846 Å while the endpoint gets **0.081 Å worse** (S30-L10 §1) |
| 30 | **H-X2** — the FAIL18 tail is conformational ambiguity in the reference (multi-model NMR depositions) | X | Registered dead if mean floor < 1.0 Å **and** it explains < 25% of the tail excess; lane X predicted ~0.6 Å and <5% | RMSD(ensemble medoid, model 1) over 111 multi-model targets | Floor mean **0.6136 Å**; the tail's floor is **−0.1776 lower**, i.e. **−5.0%** of the excess — the sign is reversed. **PRICED AND DEAD**; the hard targets have *tighter* deposited ensembles (S30-L10 §5) |
| 31 | **H-X3** — a generator whose *typical* member is better (the coordinator's construction) | X | Δ(set mean) ≈ −0.2 Å | Four arms plus the algebra `set_mean² ≈ B² + S²` | Every `d_set_mean` is **positive**; the best is a tie — T3_pool **+0.0409, 0.57× MDE, NOT MEASURED**. The counting falsifier **tied its bar (24.80% vs 25%)** and decided nothing; the verdict rests on the split-half transfer of `d_B`: **+0.1597 [+0.0676, +0.2687]** — wrong sign, CI excluding zero (S30-L20) |
| 32 | "Nobody has built a **typical-good generator**" (the coordinator's premise) | X | An unexplored design | `T0_helix` — constant α-helix plus 15° jitter | Somebody had. It is the most concentrated source ever built here (S **0.5676** vs the pool's 1.5770) and **the worst endpoint in the record, 3.7892**. Its averaging gain collapses to 0.0771 against the pool's 0.5023 (S30-L20 §1(3)) |
| 33 | Lane X's own within-target **corr(S,B) = −0.4744** | X | A strong negative worth quoting | Dropping one arm | **Driven by T0_helix alone**; without it, **+0.0851**. Lane X flagged it as an artefact that must not be quoted (S30-L20 §1(2)) |
| 34 | "**Candidate generation is closed on five instruments**" (project memory, used to rule proposals out) | X | A closure | Reading the four artefacts | **Four of five measured only the endpoint and carry no oracle arm.** An untrained per-residue-type Ramachandran sampler beats the whole K=500 pool's best on **45 of 126** targets and the endpoint still worsens. The closure should read *"no achievable source produces candidates the shipped score can convert into Ångströms"* (S30-L10 §4) |
| 35 | "**The pool's error is 68% common-mode, 50.7× the i.i.d. prediction**" is a cap on generation | X | The premise of lane X's own remit | Reading the memory's **body** | `‖ē‖² = n·RMSD²` **exactly** — the numerator *is* the output RMSD, not independent evidence about it; and the same file says *"f is NOT a screen and must never be used as one"* (S30-L10 §6) |
| 36 | **Finite-shot CVaR estimator bias** gives the VQE an incentive to shape its distribution for the estimator | L | A candidate mechanism for "the objective improves and the structure does not move" | Four checks at the deployed cell (n=9, α=0.18, 2,048 shots, tail k=369) | **REFUTED, in the opposite direction**: the empirical lower-tail CVaR is *pessimistically* biased; bias +0.0000…+0.0025, no trend with support, argmin unmoved (penalty +0.0000), and concentration is **penalised by less than one sd of the shot noise** (S30-L8 §4). Closure lapses if any lane proposes a low-α or low-shot arm |
| 37 | **EDM projection / triangle repair** on the predicted distance matrix | L | 58% of the prediction's freedom is metric inconsistency — parameter-free, no new information needed | Grepping the codebase before proposing it | S19 had measured it and left a **named prohibition**. Refuted in the *unexpected* direction: a magnitude-matched **incoherent** field is worse on both realisability measures (defect 0.340, 10.85% violations) and lands **1.24 Å better** (S30-L18 §1) |
| 38 | **E3** in its obvious form — a candidate source whose bias is not the same μ | L | Change the generator, change the shared bias | S24's quality-matched provenance cosine | **0.9432 against a 0.9330 within-source control** — two candidates from completely different sources are *more* aligned than two draws from the same source. μ is a property of **the prior candidates are scored against**, not of where they come from (S30-L7) |
| 39 | Lane L's **"the AMBER-relax benefit concentrates on divergent pools", 2:1** | L, then G | Lane L lowered it to roughly even on lane F's Rg-dispersion null; lane G registered **2:1 against** | Lane G measured it | **CONFIRMED on every form.** ρ(d, DISP_rmsd) **−0.316**, fold-preserving permutation null **p₂ = 0.000**; high-dispersion half −0.0429, low half −0.0012 — **97.2% of the gain in one half**; **length-matched split −0.0406 [−0.0481, −0.0319], 3.56× MDE, 5/5 folds**. **Both lanes' downgrades are withdrawn** (S30-L26 Q1; report §A.7) |
| 40 | Lane G's own registered **2:1 against** concentration, reasoned from `corr(S,B)=+0.085` and the 68% common-mode result | G | It does not concentrate | Its own measurement | Wrong, and the reasoning was wrong too: **the relax does not repair the bias, it repairs the geometry the averaging destroyed**, and the destruction scales with spread — ρ(DISP, contraction) = **+0.947**, and partialling DISP collapses contraction's effect from −0.276 to **+0.078** (S30-L26) |
| 41 | The confirmed divergence effect is a **tail intervention** (the coordinator's inference on commissioning Q1) | G | "If it concentrates on divergent pools it is a tail intervention, and the tail is where the prize is" | Divergence and the production tail measured as separate objects | **The premise holds and the conclusion does not.** FAIL18 dispersion vs the 108: **+0.4583 [−0.3845, +1.0472], 0.47× MDE — NOT MEASURED**; relax gain on FAIL18 **−0.0113** against the 108's **−0.0239** — *the wrong direction*; targeting the divergent half buys **−0.0215 against −0.0221**, i.e. nothing (S30-L26) |
| 42 | The divergence effect is **ten big winners relabelled** | G | The real risk — 10 targets carry 66.8% of the whole-sample gain | A drop-top ladder | Survives dropping the ten biggest winners (**−0.0310, 5/5 folds**) and dies only at drop-top-20, where the arm has **no gain left to localise** (+0.0028). A gradient, not a relabelling (S30-L26) |
| 43 | The divergence effect is **chain length** | G | ρ(n, DISP) = +0.267 made it necessary to check | Partialling and a length-matched split | Partialling n **strengthens** it (−0.316 → −0.338, p₂ 0.000); ρ(n, d) = +0.036; the length-matched split is the strongest form of the result (S30-L26) |
| 44 | **F-G2** — a chiral functional carries nativeness signal (the one family theorem G1 leaves open) | G | Fires at anchor contrast **≥ +0.10** with fold CI excluding zero **and** exceeding its own achiral twin by ≥ 0.7× the paired MDE | WRITHE, CHIRAL3, CHIRAL3_LONG on lane R's bit-identical ladder | Best: WRITHE **+0.0405 [−0.0255, +0.1015], 0.41× MDE, 3/5 folds**, against a max-over-3 sign-flip null whose **mean is +0.0408, p_max 0.430**. CHIRAL3 **−0.0363**, CHIRAL3_LONG **−0.0394** (wrong sign). WRITHE minus its achiral twin: **+0.0170, 0.29× MDE**, CI includes zero (S30-L26 Q2) |
| 45 | Lane G's own registered **mechanism** for why F-G2 would fail — the chiral coordinate is near-constant at n = 9–16 | G | 4:1 that F-G2 does not fire, *because the class is degenerate* | The occupancy pre-check the coordinator required **before** any contrast | **REFUTED.** Occupancy WRITHE 0.685, CHIRAL3 **0.953**, CHIRAL3_LONG **0.965**, against **DIS 0.345 (pool) / 0.231 (ladder)**. The chiral axis is *more* fully exercised than the shipped cost is. **The negative is stronger for it** (S30-L26; `s30/results/s30_G_chiral.json :: PRECHECK_occupancy`) |
| 46 | **WRITHE's preference contrast +0.1641** — the only channel all sprint to clear *both* of lane R's registered preference clauses | G | 2.17× MDE, 5/5 folds, max-null p = 0.000, `pref_near` 0.771 | Lane G's own kind-matched statistic, applied before it was quoted anywhere | **Cross-kind**: `pref_pool`'s control keeps **deposited** coordinates while the near rungs are **ideal rebuilds**. Kind-matched read: the rebuilt native sits at the **0.6061** percentile of its own ladder (\|WRITHE\| **0.6732**), both **worse than the 0.5 chance line**, against DIS's 0.2876. ****Second** instance of the cross-kind confound this sprint (S30-L26; the widening result was circular, which §A.7 records as a different defect)** |
| 47 | **P1** — some native-free feature set prices the ORACLE displacement out of fold | P | Bars registered at **1.96% / 12.82% / 39.44%** R² before the regression existed | Everything thrown at `e` at once, 8 arms × 5 strata, sign-equivariant, leave-fold-out | Best arm **R² 0.0083**, excess over the matched-dimension control **+0.0095 at 0.75× MDE — NOT MEASURED**. Nothing clears the bar on any stratum. **P1–P5 all hold** (S30-L25 §1) |
| 48 | A **global 5-number** separation profile is deployable | P | The ORACLE profile is worth −0.5740 Å | Applying the leave-fold-out global profile | **+0.0036 Å, 0.11× MDE**. The profile's per-target dispersion is **~9× its mean**, so a global constant is arithmetically a no-op (S30-L25 §3) |
| 49 | The **pool** can estimate the prior's separation profile | P | The only genuinely non-tautological arm in the design | `NF_POOLPROF5` | **+0.0395 Å, wrong sign** — regressing the prior toward the pool *hurts*. `pool-error-is-68-percent-common-mode` arriving at the prior, a third independent site (S30-L25 §3) |
| 50 | The **implied-endpoint conversion** from ρ is a fair summary of a correction's value | P (correcting every lane, itself first) | ρ → 3.0338 Å | Applying the fitted displacement out of fold | Emits **3.0519 Å against production's 3.0483 — slightly worse**. The conversion flatters a no-op by 0.018 Å. The sprint's fourth checklist entry (S30-L25 §1) |
| 51 | **Per-residue** prior corrections beat the 5-number profile | P | n parameters > 5 | `ORACLE_PERRES` | **−0.5714 vs −0.5740** — extra parameters add **nothing** (S30-L25 §3) |
| 52 | A single **offset** or a single **stretch** corrects the prior | P | The 1-parameter nestings | Both measured ORACLE | `ORACLE_OFFSET1` **−0.1443, 0.87×**; `ORACLE_STRETCH1` **−0.1501, 0.85×** — both **NOT MEASURED** (S30-L25 §3) |
| 53 | The **long-range** R² is also ≈ 0 (the coordinator's registered prediction, at 3:1) | P | ≈ 0 | Per-pair out-of-fold regression restricted to \|i−j\| ≥ 7 | **FALSIFIED on the raw statistic: +0.1959.** The coordinator's objection to S30-L25's aggregate statistic was right and his prior about it was wrong (S30-L27 §2) |
| 54 | …and therefore **long-range R² is a route** | P | +0.1959 numerically clears the 12.82% primary bar | Nesting, plus the category check | **+0.1406 is calibration** (already closed: `error-shape-not-mae-decides-ranking`; S25-L2 measured calibration making RMSD worse); **−0.0184** for the whole uncollapsed posterior; only **+0.0758 is genuinely new, from the pool**. And the bar itself does not apply — it was registered for the **aggregate displacement** R², not a long-range-restricted **pair-space** R² (S30-L27 §2, flagged by lane P as a category difference, not a clearance) |
| 55 | A **fitted prior corrector** helps once it is accurate enough | P | R² 0.2355 out of fold against a −0.04 control is a good predictor | Applying it through production's own downstream | **+0.0554 Å WORSE.** At *identical* R², an ORACLE-constructed i.i.d. corrector emits **−0.2466 BETTER** — a **0.30 Å swing at matched accuracy** (S30-L27 §3) |
| 56 | Any **native-free sign channel** supplies the five signs per target | P | The five-bit question | Seven channels against the **plausible** always-positive baseline (the sign of the leave-fold-out global profile), per bin | Baseline is **0.556 / 0.571 / 0.579 / 0.643 / 0.627** — **0.627 is free at long range**. Best channel POOLDIS75 **0.667, excess +0.040, 0.25× MDE**; **every channel emits worse than production** (+0.0401 to +0.1492) (S30-L27 §1) |
| 57 | The **"two qubits" synthesis** — a selecting readout over a wider register, aimed at the tail (the coordinator's, 13:07) | coordinator | Six lanes converging on one intervention | Lane F's own measurements, one minute later | Filter width is closed by ceiling (the **ORACLE global argmin over k IS the shipped 75**), and widening-rescues-the-tail does not replicate on either filter-independent tail. **Corrected at 13:08, before it was quoted** (STATE, "CORRECTION TO THE SYNTHESIS BELOW"). *The register half is not closed — see §7.4 #3* |
| 58 | "**LEG_torsion is at chance**" | coordinator | Lane R's result, in the coordinator's words | Lane R recomputing from its own artefact | Anchor contrast **+0.024, fold CI [+0.010, +0.038] — excludes zero.** It *does* order nativeness slightly above its control; it fails by being a quarter of the margin, not by being noise (S30-L19; report §A.1) |

### 7.2 Killed by theorem

These cannot be reopened by a better instrument. Each is stated with the model class it is relative
to, because a theorem without its quantifier is an overclaim — a distinction lane L insisted on and
the coordinator asked for explicitly.

| # | hypothesis | lane | the theorem that killed it | quantifier / scope |
|---|---|---|---|---|
| 1 | **The CVaR tail stops being a prefix** under some condition (the charter's own question, and S29's) | T | **T1**: at every KKT point of `min_{Λ(p)} V`, `∇V(λ*)_x < μ ⇒ λ*_x = p_x` and `> μ ⇒ 0`. The tail is **always** a prefix — of the order induced by `∇V` at the optimum | Any `V` differentiable on the box-simplex. The question is mis-posed; the load-bearing one is *"the tail stops being reproducible by one classical sort iff `λ → ∇V(λ)` has more than one fixed point"* (S30-L9) |
| 2 | S29 §4's "**the first formulation whose classical counterpart genuinely goes away**" | T | §4.3 (keep the order exogenous, for well-posedness) **provably guarantees** what §4.5 claims to escape: with an exogenous order the emitted set is still `argsort(E)[:m]` and the endpoint channel is still the integer `m` | S29 proved a **free** subset optimum is non-prefix; it did not show the lifted CVaR readout can reach it. §4.5's 3-state witness needs `p₃ = 0` exactly (S30-L9 §2) |
| 3 | **Submodularity guarantees** (Nemhauser–Wolsey–Fisher `1−1/e`, Das & Kempe `1−e^{−γ}`) apply to set selection here | L | Both require **monotone**. S29-L25 states in its own words that `V` is *"neither additive nor monotone"* | **Monotonicity, not submodularity, is what we fail first** — every submodularity guarantee in the plan was void before it was invoked (S30-L8) |
| 4 | **DPPs / facility location / diversity-aware selection** | L | `V(S) = f(mean_S W)` factors through the **centroid**, so `V` is constant on centroid-equivalence classes | There is no diversity structure to exploit. Rejected by algebra, not by trial (S30-L8) |
| 5 | The **common mode is estimable from the pool** with enough members | L | Under `w_k = t + μ + d_k`, the likelihood depends on `(t, μ)` **only through `t + μ`**, so `(t, μ)` and `(t−c, μ+c)` induce identical data distributions for every `c` and every K | Relative to model class **M**, and **from the pool only** — *not* "invisible to any method". Exactly three escapes exist (E1/E2/E3), and the record prices all three. `m_eff = 224.18/160.36 = 1.398` (S30-L7) |
| 6 | **PCA / factor models / ICA** recover the shared bias | L | Centring removes exactly the component in question before the decomposition runs; factor models recover the direction only up to the sign-and-scale gauge | Explains why S29-L47's ORACLE global η for the PC1 family is **`+0.0000` exactly** — an exact zero is something being identically zero for a reason (S30-L7) |
| 7 | **E1** — a prior on the shared bias's *form* | L | Three independent closures: (a) the optimal rescale `s* = ⟨c,t⟩/‖c‖²` is a function of the invisible component and of nothing else; (b) a corrector trained on the predictor's own features inherits its error structure; (c) S29-L47 measured four global scalars whose ORACLE optima are 0.0–0.6% of their per-target gains, **two exactly zero** | Brynjarsdottir & O'Hagan (2014) require an informative prior on the discrepancy's *shape*; (a), (b), (c) are the three shapes available here (S30-L13) |
| 8 | A **per-residue sum** can order nativeness | R | A structure can be locally perfect everywhere and globally wrong. Matched local statistics make local channels blind **by construction**, and the measurement (§7.1 #21) confirms it | *On this instrument* — a kind-matched, budget-matched torsion ladder with floor 0.347 Å. **"A sum of per-residue terms cannot see a lever arm"** (S30-L19 D1) |
| 9 | A **non-chiral single-structure channel** can carry information the distance map does not | G | **G1** (classical MDS): `G = −½JD²J = XXᵀ`, so `D` fixes the centred coordinates up to `O(3)`; after proper-rotation Kabsch the residual ambiguity is **exactly reflection**. Hence `S` is a function of `D` **iff** `S` is reflection-invariant | **Single-structure** channels only, sequence held fixed. Does **not** bound set-referenced channels (bucket b, closed separately) or a fourth reference this project does not have. The theorem does not depend on chain length; only the emptiness does (S30-L26) |
| 10 | The charter's candidate list — contact topology, Rg profile, inertia/asphericity, end-to-end distributions, the separation profile, contact order, excluded volume, packing density, burial/SASA — are **new channels** | G | Corollary G1a: each is reflection-invariant, hence a function of `D` | *"They are not new channels; they are coordinate systems on the distance map."* The list falls on **one** ground, not one at a time (S30-L26) |
| 11 | There is a **fourth reference** to score a shape descriptor against | G | Corollary G1b: this project has exactly three — the predicted distogram (class M, bounded by theorem 2), the pool (typicality; anchor contrasts CONS −0.056, DMAP_CONS −0.053, POOLGO −0.050), and universal physics (**target-independent by construction**, so it cannot supply a per-target sign) | Explains lane R's split verdict rather than merely surviving it: **all 43 of its channels are either per-residue sums or `D`-functionals against one of those three references. There was no fourth kind in the library to test** (S30-L26) |
| 12 | S28-L48 (**"20 of 31 scorers prefer production to a 0.25 Å ORACLE structure"**) establishes that nativeness is unrecognisable | R | Blau & Michaeli's perception–distortion theorem predicts the **cross-kind** half from the construction alone: every rung differs in kind (a projected coordinate average vs a circuit amplitude readout vs a least-squares fit vs perturbations of PROD) | **Posted before lane R's own numbers existed**, on the coordinator's instruction. It survives as *"no native-free scorer prefers a near-native structure of a DIFFERENT CONSTRUCTION to the production average"* — which is the real selection setting (S30-L1) |
| 13 | A **dispersion term** can see the pool's common-mode error | T | `tr(Σ_λ)` is invariant to a common shift of the pool, so by the S23-L9 identity it operates on the **idiosyncratic 32% only** | Stated in advance of the measurement that then closed the class anyway (S30-L9 §4) |
| 14 | **Errors-in-variables** identifies the reliability ratio from two error-laden measurements | L | EIV needs the two measurements' errors to be **independent**; S19 §4 measured that our two (the distogram profile and the pool profile) **share** the bias | Not merely unverified — measured false (S30-L18 §6) |
| 15 | **Reiersøl (1950)** — EIV identified without an instrument when the latent is non-normal | L | Needs a linear relation between two *observed* error-laden variables; we have one prediction and **no second observable of the same quantity** | Recorded because it is the standard answer to "identify without an instrument" and someone will propose it (S30-L18 §6) |
| 16 | **Multichannel blind deconvolution** | L | Requires the channels to be **coprime** — no common zeros. A shared bias **is** a common factor | Coprimeness is exactly what S19 §4 and S24's 0.9432 measure to be absent (S30-L18 §6) |
| 17 | **Instrument calibration** against a reference standard | L | The anchor is the native — this is E1/E3 again, and S29-L4's control-variate rejection | The pattern across all five measurement-error families is **one condition under five names**: two views whose errors are independent (S30-L18 §6) |
| 18 | The **reference state's compactness loading** is a measured empirical pattern (S29-L50) | L | `u(d) = −kT ln[P_obs/P_ref]` contains the **separable** term `+kT Σ ln P_ref(d_ij)`, a pure function of scale. For DOPE's ball reference `P_ref(d;a) = (1/a)g(d/a)` exactly, so it is scale-invariant **iff** the reference tracks the candidate's size | Forced, not fitted: a uniform 10% contraction of a 13-mer moves a fixed-reference score **−0.256 kT/pair (~14 kT)** with zero shape change; a size-matched reference moves **exactly 0.0000** (S30-L6). Confirmed independently at 58.0% of the Gram trace (S30-L22) |
| 19 | The readout's **7 bits are a channel** through which 1.44 arrive — the charter's premise | T | Bits are not conserved across an index: **the 500 deposited backbones ARE the information and the index only names one** | 7 index bits move the ORACLE ladder 4.1080 → 1.8978 Å, which through the displacement bound is ρ = 0.8869 = **36.63 displacement bits, a 5.23× ratio**. *"The readout's 7 bits are worth five times their face value, and the system cannot supply one"* (S30-L14 §4) |

### 7.2b Killed by price — correct, and arithmetically unaffordable

Nothing in this block is wrong. Each is a construction that works and costs more than the budget.

| # | hypothesis | lane | the price | verdict |
|---|---|---|---|---|
| 1 | **The sparse weighted readout** — S29's last ladder class not closed by ceiling (2 members with ORACLE weights emit 1.4315 Å) | Q | Support `log₂C(\|S\|,s)` + weights `log₂C(L+s−1,s−1)`, counted in the same currency the argmin spends | At **every** budget B=3…9 the plain argmin over the top-2^B beats the best fully-priced sparse arm by **+0.1622 to +0.7267 Å**. The sparse family needs ~15 bits to reach what the argmin gives at 9 (S30-L11 F4) |
| 2 | The famous **1.4315 Å** sparse arm is a route | Q | 2 of 500 = **16.93 bits** *plus* an unbounded continuous weight channel | The same pool's argmin reaches **1.7078 Å for 8.97 bits and nothing else**. On the built chain, 2-of-top-75 with free weights is **2.1683 Å at 11.4+ bits** against the top-128 argmin's **2.1435 Å at 7.0 bits** — better while spending 4.4 fewer bits (S30-L11) |
| 3 | **Torsion / configuration encodings** for the quantum register | T | `2n·log₂(k)` bits: at mean n = 12.96, **25.9 bits at k=2, 41.1 at k=3, 51.8 at k=4, 77.8 at k=8** | Against **7** deployed qubits (`core/pipeline.py:181`) or 9 in the harness — **0.27 bits per torsion where 1 bit names a single Ramachandran basin.** 4–8× the register this project has ever run, before expressivity is discussed (S30-L14 §5) |
| 4 | **Subset cardinality** is a competitive bit allocation | T | The value-of-a-bit law `−dD/dR = 0.2191·(D − 1.3312)` makes allocations comparable only through floor and tail index | Candidate identity **0.132 Å/bit** at R=7 against subset cardinality's **0.044 — dominated 3.0×**, which *explains* the 3.5× S29 measured empirically (S30-L14 §5) |
| 5 | **Decorrelation** — find a source whose errors are orthogonal to the distogram's (the coordinator's brief to lane P) | L | `ρ_max² = (r₁² − 2c·r₁r₂ + r₂²)/(1 − c²)` in closed form | A **perfectly orthogonal** new channel must itself carry ρ = 0.3398 where alone it would need 0.3580 — **orthogonality is a 5.1% discount on the requirement** (1.6% for 2.50 Å), and the channel must be **3.01× better than anything we own**. A perfectly orthogonal channel merely as good as our best buys **0.02 Å**. Reaching 3.00 Å by stacking needs **10.1 mutually orthogonal** channels against a measured Gram stable rank of 2.057. **Decorrelation is not the lever; skill is** (S30-L18 §4) |
| 6 | The **tail-then-aggregate lift** with an endogenous order opens a combinatorial search | T | **T1b**: reachable tails are halfspace cuts of VC dim `d+1`; by Sauer–Shelah the class carries `log₂ Σ_{i≤d+1} C(500,i)` bits | Against **300.6 bits** of free 75-subset choice: **24.3 bits** at `d_eff = 2` (stable rank), 50.4 at k90 = 6, 126.7 at k99 = 21, 175.5 at the nominal 33. A **250-to-275-bit collapse**. The registered rank-collapse rule (stable rank < 2.0) fired on a **76-second** pass over data already on disk, rather than a lane-week (S30-L9 §3) |
| 7 | **Ensemble / consensus / bagging** methods transfer at their published gains | L | Every such method prices its gain in K, and this pool's effective independent size is **m_eff = 1.398** | *"The honest conversion factor from any such paper to this instrument is (1.4 / K_theirs)."* The 75-member pool carries the statistical content of about **1.4 independent members** (S30-L7) |
| 8 | **Solving the set-selection problem better** is worth something | L | By Maurey's empirical method any hull point is within `R/√m`; with the project's own dispersion `√63.82 = 7.99` that is **0.22 Å at m=75**, 0.86 at m=5, 1.36 at m=2. Frank–Wolfe on the simplex runs in **seconds** and strictly upper-bounds any circuit on this objective | *"Exhaustive search already solved it (S29-L25) and the answer was worse than production. What is scarce is the information needed to SPECIFY a good set, and no solver supplies information"* (S30-L8) |

### 7.3 Every control run and what it ruled out

Charter item 25. One row per control. **Matched to** names the arm the control was matched to;
**ruled out** is what it eliminated, and where a control ruled out *less* than it was credited
with, the row says so in bold. Several rows are controls that killed the result they were run to
support; one refuted its own lane's stated mechanism and made the conclusion stronger for it.

#### Order-statistic and best-of-K controls

| control | matched to | draws | what it ruled out |
|---|---|---|---|
| **Best-of-random-75 from the same pool** (`PREREG_S30_F1`, seed 30001; the registered 2,000-draw sampler was replaced by the exact closed form `P(min rank ≥ i) = C(n−i,k)/C(n,k)`, cross-checked against the sampler to max \|Δ\| **0.0258 Å**) | the 500→75 filter, in its own operator space | exact + 2,000-draw cross-check | That the raw ORACLE-best drop from 500 to 75 is evidence of anything. **It did NOT rule out the stratum being defined by the outcome** — see S30-L23 below |
| **Best-of-random-500 from the universe** (seed 30001) | BLOSUM retrieval | 2,000/target | Retrieval as the tail's culprit: FAIL18 **0.01× MDE**, all strata NOT MEASURED (S30-L2 §3). This is the row that revised project memory |
| **`best_of_k_within` across-target null + split-half transfer** (lane Q) | the 128-direction halfspace/quadric search | 128 directions × 126 targets | The cut classes' apparent ceiling: **225% / 203% accounted**, transfer **+0.0323 (−3%) / −0.0204 (2%)**, k_eff 78.8. *"To reopen it, someone must exhibit a rule that produces a direction, not a larger search over directions"* (S30-L12) |
| **Common direction bank** — columns made comparable by a deterministic native-free sign convention, so column k is the same *rule* on every target (lane T's M6, run at lane Q's request) | lane T's own K = 5,000 search | 5,000 | Made `best_of_k_within` well posed and **retracted lane T's headline**: 196%/188% accounted, transfer 19%/20%, and the transferable rule lands **+0.4724 Å WORSE than production, 1.88× MDE, 5/5 folds, 39W/87L** (S30-L15 §3c) |
| **`FREE` = best-of-5,000 random 75-subsets at matched budget** (lane T's registered null for M4) | the halfspace class | 5,000 | HALFSPACE − FREE = −0.8568, 4.93× MDE — lane T read this as *"the class is real, not best-of-K"*. **This control was credited with more than it ruled out, and lane T said so:** it asks whether a structured class beats an unstructured one at equal budget, **not whether the winning direction is the same direction twice** (S30-L15 §3c) |
| **Random-support arm scored as the mean of 8 draws, never the per-target min** (lane Q, contract rule 8) | the sparse readout's support rule | 8 | The per-target minimum as an arm: min-of-8 gain −1.1056 Å is **108% accounted** by the across-target null, split-half transfer **+0.0007 (0%)** (S30-L11) |
| **Greedy-gap control** — exhaustive s = 2 over the top-64, `C(64,2) = 2016` pairs enumerated per target | the forward-greedy support search | exhaustive | That greedy supports are lower bounds: **2.1180 exhaustive against 2.1427 greedy, +0.0247 Å** (S30-L11) |
| **Split-half transfer of `d_B`, 400 splits** | lane X's "wins on both set mean and endpoint" counting arm | 400 | The counting arm, which had itself **tied its bar (24.80% vs 25%)**: transfer **+0.1597 [+0.0676, +0.2687]**, wrong sign, CI excluding zero (S30-L20) |

#### Matched-dimension and matched-space controls

| control | matched to | what it ruled out |
|---|---|---|
| **Matched random 21-dimensional subspace** (lane D) | the ORACLE **per-target** field combination | That ρ = **0.9491** is a property of the fields. The control reaches **0.8095** by dimension counting; the excess is **+0.1396, SE 0.0073, 6.83× MDE, 5/5 folds**. *"Without that control this would have gone into the record as a 0.97 Å result"* (S30-L21) |
| **Rank-matched random subspace at every rank r** (lane D) | the rank curve | That the excess is an artefact of the near-degenerate tail: at rank 4 the real set reaches **0.6627 against 0.3299, +0.333, 6.31× MDE** (`s30_D_gram.json :: rank_curve`) |
| **C-DIM — matched-dimension random orthonormal K-frame** in the rigid-free space, features recomputed and model refit (lane P) | the 14-direction native-free basis | Two things at once. On the regression: control **−0.0012**, excess **+0.0095 at 0.75× MDE — NOT MEASURED**. On the ORACLE capture: control **0.4444** against the basis's **0.8273**, excess **+0.3829, 10.17× MDE, 5/5 folds, CI [+0.369, +0.392]** — the sprint's cleanest positive. Lane P checked rather than assumed that the frame is inside the pool's span (at n = 9–16 the 75 deviations span all 21–42 dimensions), so this is **not** an instance of `control-must-match-the-operators-space` (S30-L25) |
| **Matched controls in the operator's own space for the source law** (lane X): two-term / mean-only / best-only | the ceiling gate | The gate. \|resid\| **0.0153 / 0.1198 / 0.3330 Å**; sign correct **5/6 / 5/6 / 1/6** (S30-L10 §2) |
| **The achiral twin `\|X\|`** — same functional, same scale, same ladder, differing in **nothing but chirality** (lane G; `\|X\|` is reflection-invariant, hence a `D`-functional by G1) | each chiral channel | That chiral content contributes anything: WRITHE − \|WRITHE\| = **+0.0170, 0.29× MDE**, CI includes zero; CHIRAL3 and CHIRAL3_LONG are **worse** than their twins at ordering (−0.0760, −0.0808, 5/5 folds). *"Whatever WRITHE can do, its own reflection-invariant shadow already does"* (S30-L26) |
| **Matched local/global feature counts** — 84 local vs 80 global, with the local block **ORACLE-advantaged** (handed per-residue circular deviation from the native anchor) (lane R) | the D1 locality decomposition | Capacity as the explanation for the local null. With that advantage and at matched capacity the local block still gives **ΔR² −0.089** (S30-L19 D1) |
| **The random m-subset**, 8 seeds per cell, common random numbers across k, seed 30003 (lane F) | the averaging operator | Seed noise in the width sweep; the seed sd is printed beside every cell (`PREREG_S30_F2`) |

#### Plausible zero-information baselines (not coin flips)

| control | matched to | what it ruled out |
|---|---|---|
| **"Always predict positive"** — the sign of the leave-fold-out global profile, per separation bin (lane P) | every native-free sign channel | The whole sign channel. The baseline is **0.556 / 0.571 / 0.579 / 0.643 / 0.627** and **0.627 is free at long range**; no channel clears 1× MDE over it in any bin, and every one emits worse than production. The matched-accuracy i.i.d. ladder crosses zero at **~0.63–0.65 realised accuracy — exactly where the free baseline already sits** (S30-L27 §1) |
| **`T0_helix` — a constant α-helix plus 15° jitter** (lane X) | the "typical-good generator" hypothesis | That concentration pays. It is the most concentrated source ever built (S **0.5676** vs the pool's 1.5770) and **the worst endpoint in the record, 3.7892**. It also *created* an artefact: the raw within-target corr(S,B) = −0.4744 is driven by this arm alone, and lane X withdrew it (S30-L20) |
| **C-ZERO — a plausible zero-information field** (constant contraction/expansion to the pool's mean Rg), explicitly *not* an isotropic random field (lane P, `PREREG_S30_P`) | the scale mode | Crediting a trivially available scale move as signal |
| **Three structured fields that sit at the signed null** — MEDOID −0.0097, MSET_50 −0.0187, EXPAND −0.0208 (lane D) | the +7.7 σ field result | The `zero-information-control-must-be-plausible` objection, *from inside the survey*: a signed mean of +0.11 is **not** a generic property of any structured displacement field (S30-L5) |
| **The two nulls for a cosine, named and separated** (meter extension E5) | any new cost's cosine | Reading a 126-target mean against the magnitude of **one** random direction on **one** target. Per-draw \|cos\| **0.1398**; target-mean null sd **0.0144**, 95% [−0.0169, +0.0307] — ~10× tighter. The shipped cosine of **−0.0339 is z = −2.45** against the second and unremarkable against the first. *"A lane reading its new cost's +0.09 cosine against 0.140 would discard a field that is 6 σ from the signed null"* (S30-L5, S30-L3 E5) |

#### Nulls for the exact statistic

| control | matched to | draws | what it ruled out |
|---|---|---|---|
| **Simulated uniform-ranking null for each of five statistics** plus KS and χ² on the full rank distribution (lane D) | S29's "worse than chance" paragraph | 20,000 | All three "worse than chance" claims. Null mean **1.4050**, null **median 0.9888**, null win rate **62.5%**; observed p = 0.16–0.75, KS p = 0.571, χ² p = 0.303 (S30-L4) |
| **Random-18 null** (`PREREG_S30_F1`, seed 30002) | every FAIL18 claim in lanes F and Q | 20,000 | That FAIL18 behaves like an ordinary random 18-subset. **It does NOT rule out FAIL18 being defined by the quantity measured on it** — the distinction S30-L23 had to introduce, and the reason `p = 0` was not a defence |
| **Max-over-channels per-target sign-flip null across all 43 channels** (lane R) | DIS's +0.347 ordering result | 500 | Multiplicity as the explanation: null mean 0.062, p95 0.105, **p_max 0.000** on both `rho_A_part` and `pref_near`. The standard the report's §6.2 says the next sprint should inherit |
| **Max-over-3 sign-flip null** (lane G, lane R's statistic) | the three chiral channels | as lane R's | The chiral positive: observed **0.0405** against a null **mean of 0.0408**, p95 0.0852, **p_max 0.430** — on its own null to three decimal places (S30-L26) |
| **Fold-preserving permutation null** — dispersion labels permuted **within fold** (lane G) | ρ(relax gain, pool divergence) | 10,000 | Fold structure in either variable manufacturing the correlation: **p₂ = 0.000** on all three dispersion variables (S30-L26) |
| **Uniform-effect null** for concentration (lane G; also lane D on the meter) | the high-minus-low dispersion gap | — | Reading concentration off a raw drop-top threshold. Lane G's observed gap sits at **percentile 0.001** of [−0.0264, +0.0289]; lane D's CA preference contrast sits at the **51st** percentile, no flag (S30-L26, S30-L24) |
| **C-PERM — target-label permutation within length strata** (lane P) | the feature→target link | — | Spurious association: `CPERM R² = −0.0388` (S30-L25) |
| **Row-permutation and random-feature controls beside every long-range block** (lane P) | the per-pair regression | — | Row-permutation **−0.0254**, random-feature **−0.0381**. **Lane P flagged that the random-feature control does NOT catch the shared-referent floor**, because random features contain no `expected` — which is why the +0.1959 had to be nested rather than compared to the control (S30-L27 §2) |
| **Label-shuffled fit** (lane R) | the leave-fold-out combination | — | Leakage: shuffled pref 0.439 against a 0.325 control, margin +0.114 (S30-L19) |
| **R = 8 matched random-signed draws** (meter extension E4) | S29's single seed-0 control | 8 | **Withdrew one S29 number and confirmed another.** Built chain **+0.0635 → +0.0357**, 0.92× → **0.56× MDE**, fold CI now spanning zero; CA **+0.1746 → +0.1716**, still 1.84× MDE, 5/5 folds. Relative draw noise **59% on the chain vs 27% on CA** (S30-L24) |

#### Stratum, circularity and definition controls

| control | matched to | what it ruled out |
|---|---|---|
| **Two filter-independent tails** — worst-18 by pool mean, worst-18 by ORACLE best-in-pool (lane F, carried by lanes F, P and G) | FAIL18 | Three different things, and the pattern is the point. The **filter effect survives** (+0.6708 p = 0.0148; +0.6320 p = 0.0300, against FAIL18's +1.7674 — inflated ~2.7× by its own definition). The **widening effect does not** (≈0, and reversed at **+0.7313**). The **score's skill collapse does** (+0.3798, +0.3438 against the 108's +0.6446) — *"the first tail statement in this lane that survives on all three definitions."* Lane P's separation-profile prize also survives (**−1.405 / −1.164** against the 108's −0.359) (S30-L2 §6, S30-L16 §4, S30-L17 §1, S30-L25 §4) |
| **Stratification by `keep`** — the number of in-band members the filter actually retained (lane D, 200 draws/target, seed 3030) | lane F's F1c FAIL18 row | **The headline itself.** `keep = 0` (which *is* FAIL18) −1.7622; `keep = 1` **−0.0363**; `keep ≥ 9` **+0.0380**; Spearman(keep, effect) over the 108 = **−0.046**. There is **no gradient** — the effect is a step at the selection boundary, and one retained in-band member removes 98% of it. By definition `top_best − pool_best ≥ 1.5` on these 18, so the predicate alone forces **\|effect\| ≥ 0.869 Å = 49% of the headline** (S30-L23). **This is the control that named the failure mode: *a matched control in the right space does not rescue a stratum defined by the outcome.*** |
| **Difficulty control** — `G_filt` regressed on production RMSD over the 108, 95% prediction band (lane F, F1b clause 3) | the filter's tail failure | That the filter's failure is a restatement of difficulty: all 18/18 sit above the band, mean z = **+4.99**. **This control shares the defect S30-L23 found** — `G_filt` on FAIL18 is the quantity the predicate forces (see §7.4 #1) |
| **Linear regression-to-the-mean control** (lane F) | the widening effect on FAIL18 | That widening's tail gain is production having been unluckiest there: predicted **+0.069**, observed **−0.619** (S30-L16 §4) |
| **44/126 targets with a filter worse than 90% of random 75-subsets, 26 of them outside FAIL18** (lane F) | the FAIL18-only reading | That the pathology is confined to the defined stratum: those 26 have a production mean of **3.1360** against **2.4095** for the other 82 non-FAIL18 targets (S30-L2 §6) |
| **The kind-matched percentile of the rebuilt native inside its own ladder** (lane G; chance exactly 0.5 by construction) | WRITHE's +0.1641 preference positive | **Lane G's own positive.** WRITHE **0.6061 [0.5609, 0.6522]**, \|WRITHE\| **0.6732**, both worse than chance at 5/5 folds, against DIS's **0.2876**. Withdrawn by lane G **before it was quoted anywhere** (S30-L26) |
| **`pctile_nat_in_A`** — the share of the native's own Rama-resampled perturbations that score better than the rebuilt native (lane R) | all 43 channels | That anything puts the native first. DIS 0.291, LEG 0.283, LEG_torsion 0.503, **RAMA 0.639** — under the shipped cost, 29% of the native's own perturbations score better; under the Ramachandran term, 64% do (S30-L19) |
| **Pool-member control** — `pref` against a random *real* pool member (lane R, from S28-L36's veto) | every preference number | The whole of clause (ii). The LFO combination prefers a pool member on **100%**; RAMA's 0.640 becomes a **−0.150** contrast; LEG_torsion's 0.503 against 0.698 becomes **−0.195** (S30-L19) |
| **Anchor control** — ladder B anchored on a random real pool member, scored as the contrast `rho_A − rho_ANCHOR` (lane R; the control in the operator's own space, reused by lane G) | "orders distance from an anchor" vs "orders nativeness" | Magnitude-ordering masquerading as nativeness. **Its own reliability was measured and is low**: draw-to-draw correlation **0.147**, mean \|diff\| 0.235, from an accidental duplicate shard — so the *per-target* contrast carries that noise even though its mean is stable (`s30_R_stability.json`). Any future use should average several anchors. See §7.4 #5 |
| **D3 — the anchor control's own confound, measured before the verdict was read** (lane R) | the anchor control | That clause (i)'s +0.10 margin was unreachable for reasons unrelated to the channels: median corr(distance-to-anchor, distance-to-native) = **+0.259** (mean +0.214), and the prereg's caveat triggers only above 0.5, so **it does not trigger** (S30-L19 D3) |
| **A2 — realism flatness audit** (lane R), explicitly an audit that can only *weaken* its own positives | the ladder's construction | That the ladder's residual realism gradient explains DIS's +0.347: RAMA +0.116, EXVOL +0.161, \|Rg − median pool Rg\| +0.190 — DIS is about **twice the largest gradient**, and the gradient makes every null stronger and every positive smaller than it looks (S30-L19) |
| **The declared FAIL18 caveat** (lane G, `PREREG_S30_G`) | its own F-G1c stratum split | Pre-emptively: reported as a **descriptive** split, not a clean causal stratum, because FAIL18 is defined by the filter's recall (S30-L23) and **fold 0 contains no FAIL18 target** (1:6, 2:2, 3:4, 4:6, **0:0**), so every FAIL18/108 CI in the sprint draws from four clusters and is wide by construction (S30-L3 E2) |

#### Construction controls — properties asserted in code rather than argued in prose

| control | matched to | what it ruled out |
|---|---|---|
| **A3 — the channel object has no `nat_ca` and no `oracle_rr` attributes at all**; `ham_lib.Context` consumes the universe through exactly four keys (W, S, PHI, PSI); locked by `tests/test_s30_R.py::test_t4` (lane R) | every one of 43 channels | Native leakage, **by construction rather than by a poison run** — which is the stronger form. Every RMSD column is computed only *after* every channel value exists (S30-L19) |
| **The reflection/rotation audit asserted in code** (lane G): reflect every structure and assert `X → −X`, `\|X\|` unchanged, `D` unchanged | the three chiral channels | **It caught a real bug.** The first Klenin–Langowski sign term was **not rotation-invariant (rot_err 2.03)** and the assertion stopped the run **before any number existed**. Final audit, worst over 126 targets: `X` flips at 0.00e+00, `\|X\|` and `D` at 0.00e+00, rotation invariance 4.19e−13 (S30-L26) |
| **The sign-equivariant model with no intercept** (lane P) | the out-of-fold regression | That the fit could launder an ORACLE sign. A PC's sign is arbitrary, so every feature is either signed (flips with `U_k`) or unsigned-with-interaction; the fit is equivariant by construction. *"Without this the arm is lane D's ρ = 0.949 that was 0.809 dimension"* (S30-L25 §1) |
| **The internal consistency floor** — `w = e_j` must reproduce field *j* exactly, so a true maximum cannot fall below the best single field (lane D) | lane D's own ORACLE global arm | Its own objective mismatch: a least-squares arm returned ρ **below** the best single field, which is impossible for a true maximum. Replaced with direct mean-cosine maximisation, 22 starts (S30-L21) |
| **Nested leave-one-fold-out ridge selection inside the training folds** (lane D) | lane D's own LFO arm | Its own strawman: an untuned ridge on a stable-rank-2 Gram returned a **negative cosine (−0.218)**. The tuned version selects the largest grid value on 4/5 folds — the fitted covariance contributes nothing (S30-L21) |
| **Leakage-free field selection inside the training folds** (lane D) | lane D's own EQ11 arm | Selection leakage, **priced**: ρ **0.1139 → 0.0948**, i.e. the leakage was worth **0.019 in ρ** (S30-L21) |
| **Identifiability check registered in advance** — `r(pool mean, pool best) = −0.9898` across the K-ladder (lane X) | fitting the source law's coefficients | Fitting `(a, b)` on the family being tested. **No coefficients were fitted**; the Sprint-20 wide-set values were frozen (S30-L10 §2) |
| **The disclosed non-out-of-sample cell** — K25→K2000, computed by hand *before* the prereg (lane X, contract rule 14) | the residual bar | Its own motivating cell contaminating the test: excluded from the bar and carried separately as `disclosed_cell_NOT_out_of_sample` (S30-L10 §2) |
| **Exogenous frame** — all 500 candidates superposed onto the pool medoid **once** (lane T, M4) | the direction search | A tail-dependent frame acting as an unstated operator (S30-L15 §1) |
| **Tie handling** — no `np.argmin` on a tied signal anywhere; ties score 0.5 and preferences average over the tied set (lane R, A5) | every preference statistic | The `tie-breaking leaks the pool order` failure mode (S30-L19) |
| **Bit-identical ladder reuse** — `s30_G_chiral.py` imports lane R's `Sampler`, `make_ladder`, seed rule `crc32("s30R\|<pdb>\|<seed>")` in the same call order (lane G) | lane G's F-G2 numbers | That lane G and lane R measured different instruments: **the 0.347 Å torsion-rebuild floor reproduces exactly** (S30-L26) |
| **Reproduction gate** — the score order recomputed from the posterior, aborting unless the recomputed top-75 equals production's `sub` **set-wise on 126/126** (lane F, F2 and F3) | every k-sweep and ρ_pool row | That lane F was scoring a different order than production. Declared residual caveat: deep order *below* the cut carries a 2-in-126 chance of differing from a fresh recomputation, and every k > 75 row inherits it (S30-L16, S30-L17) |
| **Cold re-run of all ten S29 anchors** (lane D, job 1) | the meter itself | Instrument drift: all ten reproduce to **< 5e−4**; now asserted as `s30_D_meter.py verify` (S30-L3 §1) |
| **`s30/s30_verify.py`** — every quoted anchor asserted against its artefact, plus an assertion on the shape of the `ST.compare` sign-convention trap | the report's own numbers | Prose drifting from artefacts: **22 matched, 0 mismatched, 0 keys/files not found** (re-run 2026-09-20 by lane Y) |
| **A1 — independent chain reprojection** (lane R) | production's own record | That the 17-sprint-old chain numbers are stale: **3.2071 Å against the record's 3.2126**, median per-target \|Δ\| **0.0012 Å**, with one 0.513 Å outlier whose cause is the documented multi-start branch flip (S28-L18/L27b/L43) (S30-L19) |

#### A control applied before any compute was spent

| control | matched to | what it ruled out |
|---|---|---|
| **The occupancy pre-check** — `sd(X over candidates) / sd(X over candidates ∪ their mirrors)`, required by the coordinator *before* any contrast (lane G) | the chiral escape class | **This is the sprint's most instructive control.** Lane G predicted the chiral coordinate would be near-constant at peptide length and registered that as its *mechanism* for expecting F-G2 to fail. Measured: WRITHE **0.685**, CHIRAL3 **0.953**, CHIRAL3_LONG **0.965**, against **DIS 0.345 (pool), 0.231 (ladder)**, with both signs occupied (frac>0 = 0.62–0.78). **The control refuted the lane's own stated mechanism and made its conclusion stronger**: the class is genuinely exercised, so the null that follows is a real negative rather than a degenerate coordinate. *"The easy explanation is gone"* (S30-L26) |
| **Lane L's E3 falsifier, applied to lane X's remit before any endpoint run** — require the generated source's provenance cosine against the incumbent pool to fall **below** the 0.9330 within-source control | generative structural spaces | The measured value is **0.9432** (`s24/results/qmatch.json`, S24-L3, replicated by S24 lane E). **The falsifier had already failed, so no endpoint compute was spent** (S30-L7 "FOR LANE X", applied in S30-L20 §2) |
| **The rank-collapse pre-check** — stable rank of the centred pair-distance matrix, the check the charter's L9 demanded *"before spending compute"* (lane T, registered rule: < 2.0 ⇒ closed at the encoding level) | the tail-then-aggregate lift and lane Q's direction | **1.859** (median 1.865, min 1.247, max 2.721), PC1 = 55.4%, k90 = 5.61. The rule fired. *"That verdict cost one 76-second pass over data that was already on disk rather than a lane-week"* — and it is the same failure mode that killed S29's first non-diagonal Hamiltonian, caught this time **before** the spend (S30-L9 §3) |
| **The feature-space qualifier on that same number** (lane Q, adopted in full by lane T) | the rank-collapse rule itself | **Over-application of the control.** 1.859 is the **pair-distance** space; in **coordinate** space it is **3.404 (T) / 3.619 (Q)** with **k90 = 11.2, not 5.6** — *above* lane T's own 2.0 threshold. Sparse/weighted readouts and second-moment constructions act on **coordinates**. The verdict is unchanged but now rests on the direct measurements (M4, M6), not on the threshold (S30-L15 §3d) |

#### The gate and the register

| control | matched to | what it ruled out |
|---|---|---|
| **The meter's PASS/BLOCK gate at the 0.7× / 1.0× MDE rule** (E6) | every cost function proposed this sprint | *"A BLOCK on the ladder means the cost moves in the wrong structural direction and gets no endpoint compute."* **The shipped cost BLOCKs on both bases** — the chain on the ladder alone (ρ −0.4023), the cloud on the ladder *and* the cosine (S30-L3, S30-L24) |
| **The multiplicity register** (E7) | the meter's own output | Uncounted comparisons: **30 on the built chain and 33 on the CA cloud** (`s30_D_meter_DIS_{chain,ca}.json :: multiplicity.comparisons_emitted`, verified by lane Y), with the consequence stated in the artefact itself. *Flagged: S30-L24's closing line says "30 comparisons from the chain run, 24 from the CA run" — the artefact says 33; the report's §6.2 already quotes the artefact.* The honest limitation the report records stands: the sprint-wide total is a **lower bound**, because lanes counted their own and nobody maintained one register |

### 7.4 Entertained and NOT killed — what is honestly still open

Six hypotheses were entertained and are neither confirmed nor closed. Recording them here is the
difference between a ceiling argument and a tidy one.

| # | hypothesis | status | why it is not closed |
|---|---|---|---|
| 1 | **The filter is differentially harmful on hard targets** (lane F's F1b, all three clauses) | **NOT ESTABLISHED, and not withdrawn either** | S30-L23 withdrew F1c's FAIL18 row because the stratum is the outcome. **F1b's clauses 1 and 2 rest on the same quantity** — `G_filt = top75_best − pool_best` on FAIL18, which the selection predicate forces to be ≥ 1.5 (measured 2.393). No entry extends the withdrawal to F1b, and no entry re-derives F1b on a filter-independent tail. S30-L23's own words apply verbatim: *"the honest test defines the stratum on information the measurement does not reuse."* The surviving, non-circular statements are the +0.6708/+0.6320 filter effects and S30-L17's skill collapse |
| 2 | **Lane P's built-chain endpoint** for the ORACLE-sign / leave-fold-out-magnitude arm | **RESOLVED — and this row is lane Y's own false positive, kept because the cause is instructive** | Lane Y read `s30_P_chain_rows.jsonl` **while the job was still writing it** (the log prints progress every 10 targets, so it showed 110/126) and concluded the run had died. The completed artefact holds **126 rows and 126 unique PDBs**, `s30/results/s30_P_chain.json` **does exist**, and every arm is at n = 126: PROD 3.2126, ORACLE_SEPPROF5 2.6791, ORACLEsign_LFOmag **2.8867**. My registered −0.35 to −0.45 prediction is **resolved at −0.3259** — a 0.024 Å shortfall at 0.15× MDE, which is not a falsification at any resolution this instrument has. **The lesson is real and new: checking a file is the right instinct, but a file being written by a live job is not evidence about that job.** |
| 3 | **Widening the quantum register 128 → 512 helps the tail** (lane Q's incidental finding; the live half of the "two qubits" synthesis) | **OPEN, with a known circularity the sprint never discharged** | −1.9004 Å on FAIL18 against −0.1907 on the 108, difference of stratum means −1.7097, SE 0.2207, **2.77× MDE** (S30-L11). The coordinator flagged at 13:08 that FAIL18 is defined by production, production is the filtered set's average, and a bad filter is exactly what pushes good candidates down the ranking — so *"the good candidate is outside the window"* may be **produced by** the thing it is offered as evidence about. He called the filter-independent-tail check *"cheap"* and *"the gate on that whole direction."* **No ledger entry runs it.** The finding is ORACLE in any case |
| 4 | **E2 — constraint repair** | **OPEN, doubly capped, and NOT a tail intervention** | The report's §3 lists E2 in the "What was closed" table. It is not closed: S30-L7 has it *"OPEN at −0.022 Å"*, S30-L13 *"open but doubly capped"*, and S30-L26 **confirms** the divergence half at −0.0406, 3.56× MDE, 5/5 folds. What is established is narrower and sharper: it is real, mechanistically explained, divergence-graded, **0.69% of baseline**, it does **not** reach the production tail, and its restraint constant still has **no native-free selection rule** |
| 5 | **The A−ANCHOR contrast** as the sprint's discriminator of record | **VALID IN THE MEAN, NOISY PER TARGET** | Lane R measured draw-to-draw correlation **0.147** on `rho_ANCHOR_part` (mean \|diff\| 0.235) against 0.912 for `rho_A_part` and 0.982 for `pref_near`, and said future uses should average several anchors (`s30_R_stability.json`). Lane G's F-G2 verdict and lane R's clause (i) both turn on this contrast. Both verdicts are safe — lane G's best value sits on its own null mean to three decimals, and lane R's DIS clears by 2× the largest realism gradient — but the per-target column should not be reused without averaging anchors |
| 6 | **A native-free covariate, observed at inference, predicting the separation profile, generated by neither the distogram nor the pool** (lane L's design constraint, S30-L18 §3) | **UNTESTED** | Lane P tested two candidates (a global constant and the pool) and both are zero or negative, and says so explicitly: *"Lane L's design constraint survives this entry untested, and it is now a five-bit question rather than an open-ended one."* S30-L27 then converted it from an adjective into a measurable admission rule: **a prior corrector is worth building iff it LOWERS `coh` below the uncorrected error's own +0.6931**, where `coh` is the within-target correlation of the corrector's residual with the pool's common-mode pair error. Every channel this project owns **raises** it (0.693 → 0.78–0.92); the ORACLE arms that lower it to 0.54–0.59 are worth **−0.126 Å at R² = 0.16 and −0.247 Å at R² = 0.24** |

---

## 8. The twelve S29 leads: pursued, rejected, and why

Charter item 23. `BRIEF.md` §8 is explicit that the twelve are *"a leads register, not a task
list"*, that *"a sprint ignoring ten of them and finding the mechanism is a success, and one
executing all twelve and finding nothing is not"*, and that *"inventing a thirteenth lead that
beats all of these is an entirely acceptable outcome."* Six were pursued to a verdict, four were
pursued inside another lane's question, and **two were not pursued at all**. The reasons for the
two are given as they actually are, not retrofitted.

| lead | status | lane(s) | what was found |
|---|---|---|---|
| **L1 — Calibrate the posterior** | **Pursued indirectly, and closed from three sides** | L, P | The charter asked specifically to *"test the common-mode hypothesis directly — and note that if common-mode error does dominate, that is a finding about the pool, not just about calibration."* It was tested directly and it does dominate: S30-L7's non-identifiability theorem says no pool-only estimator of the shared bias exists at any K, with `m_eff = 1.398` of 75. Lane P then measured the calibration arms: a leave-fold-out global profile emits **+0.0036 Å (0.11× MDE)**, a pool-estimated one **+0.0395 Å, wrong sign**, and a fitted per-pair corrector at **R² 0.2355 out of fold** emits **+0.0554 Å worse** (S30-L25, S30-L27). Lane L had already placed calibration as row 4 of the five-instance realism law (S25-L2: calibration improves, RMSD worse). **No temperature/per-bin calibration arm was run as such**, because all three of its possible forms were priced or measured first |
| **L2 — Remove the contraction blind spot** | **Pursued, and answered in both directions** | L, D, R, G | The charter's exact question was *"can the objective distinguish a physically valid structure from an artificially contracted one without destroying legitimate structural variation?"* **Yes, by construction, and it buys nothing.** Lane L derived that a fixed-reference pair potential contains a separable pure-scale term worth **−0.256 kT/pair (~14 kT over a 13-mer)** for a 10% contraction with zero shape change, where a size-matched reference moves **exactly 0.0000** (S30-L6). Lane R built the size-invariant twins and audited them (RG_LAW relative sd 0.5009 → **3.2e−16**): on the channel lane L named it helps as predicted (DISTPOT_SI beats DISTPOT at +0.193 vs +0.136 ordering, +0.112 vs +0.056 anchor contrast) — **and several twins become *more* Rg-correlated** (CONS −0.374 → −0.807), because at peptide length shape and size are genuinely correlated. Lane D then measured the consequence that matters: the radial direction is **58.0% of the Gram trace** and **94.8% of λ₁**, while the direction to the native is **orthogonal to it (−0.0675) and anti-aligned on FAIL18 (−0.2524)** (S30-L22). Lane L withdrew the deployable half of its own proposal: **the size-matched field IS the deflated field, and deflation reallocates rank without creating any** (S30-L18 §8). Lane G added the mechanism from another side: contraction is a **shadow** of pool divergence (ρ 0.947; partialling divergence collapses contraction's effect from −0.276 to +0.078) |
| **L3 — Residual prior** | **Pursued as the sprint's decisive measurement** | P, with L narrowing it and T reformulating it | Lane T restated assumption B2 as one estimable number, `ρ_max = √(R²(e ~ S))`; lane L priced the charter's *"demonstrably different information"* requirement out of the picture (**orthogonality is a 5.1% discount; skill is the lever**, S30-L18 §4). Lane P measured it: **R² 0.0083 out of fold, excess over a matched-dimension control +0.0095 at 0.75× MDE — NOT MEASURED**, against a registered 1.96% bar, on 8 arms × 5 strata (S30-L25). The charter's explicit warning — *"do not train a second model on the same feature representation"* — was honoured and then surpassed: S30-L27 measured **why** it fails, at identical accuracy. **What can be predicted is coherent and therefore harmful; what would help is incoherent and therefore unpredictable.** Every fitted corrector raises its residual's coherence with the pool's common mode (0.693 → 0.78–0.92); only an ORACLE-imposed i.i.d. one lowers it (→ 0.54–0.59), and the swing at matched R² is **0.30 Å** |
| **L4 — Backbone-torsion channel** | **Pursued, closed negatively with a mechanism** | R, T | The charter said *"do not assume it is useful because an audit left it open. Prove it or kill it."* Killed. On a kind-matched, budget-matched ladder LEG_torsion's anchor contrast is **+0.024 [+0.010, +0.038]** — a quarter of the registered +0.10 margin — its preference is **0.503 against a 0.698 pool-member control**, and it ranks the rebuilt native at the **exact median (0.503)** of the native's own perturbations (S30-L19). **The charter's menu of alternative formulations (torsion distributions, transitions, pairwise couplings, smoothness, secondary-structure sectors, non-additive interactions) was not worked through, and the honest reason is that D1 makes the whole class a theorem rather than a search:** the RMSD signal is absent from local features (**ΔR² −0.089**, with the local block ORACLE-advantaged) and abundant in global ones (**+0.600**), so a sum of per-residue terms cannot see a lever arm. Separately, lane T closed torsion space as an *encoding*: **48 bits for 4 basins/residue at n = 12 against 7 deployed** |
| **L5 — Free-energy / subset-selection CVaR-VQE** | **Pursued and closed, three independent ways** | T, Q, L | S29 flagged this as the highest-priority quantum reformulation. **T1**: the tail is *always* a prefix — of the order induced by `∇V` at the optimum — so the charter's question ("under what conditions does it stop") is mis-posed, and S29 §4.3 and §4.5 are incompatible with §4.3 winning. **T1b** prices the endogenous-order repair at **24.3–175.5 bits against 300.6** of free set choice. **The one surviving escape** (a second-moment term, reached independently by lane T's algebra and S29's own post-mortem) was then measured by two lanes with two samplers: **+0.2059 Å, 1.81× MDE at m=75** (Q) and **+0.0860, 2.95× MDE, 5/5 folds** (T), with the structured `disp2` form emitting **3.3585 against production's 3.0483**. Lane L supplied the classical counterpart the contract requires: **Frank–Wolfe on the simplex, seconds, strictly upper-bounds any circuit on this objective.** The charter's own instruction — *"if the classical problem is better, say so clearly"* — is answered: it is, and the reason is that **no solver supplies information** |
| **L6 — Sparse weighted readout** | **Pursued and closed by price** | Q | The charter noted *"S29 estimated the relevant selection task may need substantially more than 7 bits"* and permitted larger registers *"where justified"*. Lane Q counted both channels in the readout's own currency. Three of four registered falsifiers confirmed (support carries 76.7% of the gain at s=2; weights cost only **2.32 bits**; a native-free support rule does beat random by +0.364 Å) and the fourth **refuted**: at every budget B = 3…9 the plain argmin over the top-2^B beats the best fully-priced sparse arm by **+0.1622 to +0.7267 Å**. S29's famous **1.4315 Å** arm costs ≥ 16.93 bits *plus* an unbounded continuous channel where the same pool's argmin reaches **1.7078 for 8.97 bits** (S30-L11). **Closed not by ceiling but by price** — and the confirmed falsifier F3 is the informative one: the skill exists, is real, and every native-free support arm still sits **above production** |
| **L7 — FAIL18** | **Pursued hardest of the twelve — six lanes touched it** | F, D, X, Q, P, G | **Two candidate causes excluded.** It is not pool-limited (ORACLE best of the FAIL18 pools is **2.2842 Å**, 13/18 under 3.00, S30-L2) and it is not conformational ambiguity in the reference (**the hard 18 have *tighter* deposited NMR ensembles**, floor difference −0.1776 = **−5.0%** of the excess, S30-L10 §5). **One mechanism located and measured on all three tail definitions:** the score cannot order its own pool there — ρ_pool **+0.1066 with the fold CI including zero** against **+0.6446** on the 108, random-18 null **p = 0**, replicating at +0.3798/+0.3438 on the filter-independent tails; and it is the distogram's own error (**−0.799**, −0.819 length-residualised), **83.1% of it shape, with scale refuted as the mechanism** (S30-L17). **The charter's warning was honoured to the letter:** *"a native-free regime detector is allowed but must be genuinely different from the closed router family. Do not build increasingly elaborate routers from the same posterior-derived features."* **No router was built.** Lane F enumerated the seven closed feature families in its prereg, measured all of them, and reported that its three best native-free correlates (+0.615, −0.528, −0.428) are **all inside them**. One economic fact did change: the ORACLE gate prize is **−0.1193 Å at 1.29× MDE**, about **2.5×** the −0.0474 Å against which lane C's detector was dismissed, so *"even a perfect detector is too small to measure"* is specific to lane C's operator and must not be quoted generally (S30-L16 §5). Lane G closed the last open route into the tail: E2's divergence effect is real and **does not reach it** |
| **L8 — Information-budget study** | **Pursued and answered by inverting the question** | T | The charter's premise was that 1.44 of 7 bits arrive and asked where the other 5.56 went. **They did not go anywhere — the pool is a codebook, not a channel.** 7 index bits move the ORACLE ladder 4.1080 → 1.8978 Å, which through the displacement bound is **36.63 displacement bits, 5.23×** (S30-L14 §4). The **value-of-a-bit law** `D(R) = a + c·2^(−R/γ)` fits at **R² = 0.9983** (a = 1.3312 Å, γ = 3.1636), giving `−dD/dR = 0.2191·(D − 1.3312)` and making allocations comparable through floor and tail index alone, with cardinality irrelevant. **The charter's own hypothesis is refuted:** *"the encoding may be the hidden bottleneck, and it is the least-examined component."* It was examined; **candidate indexing wins by 3×** (0.132 vs 0.044 Å/bit, which *explains* the 3.5× S29 measured), its ORACLE floor is 1.3312 Å against a 2.50 Å target, and **it is not the bottleneck.** What is missing instead has a dimension and a price: **six per-target coefficients, 0.855 bits for 3.00 Å and 3.782 for 2.50 Å inside a subspace the pool hands you for nothing** (S30-L15) |
| **L9 — Non-diagonal Hamiltonians, only if actually new** | **NOT PURSUED as a build — closed upstream by the check the lead itself demanded** | T, Q | The lead's own text says *"check for rank collapse before spending compute — that check is cheap and would have saved the previous attempt."* It was run first, as the second sentence of lane T's theory work: stable rank of the centred pair-distance matrix **1.859** (median 1.865, max 2.721), PC1 55.4%, k90 5.61 — **lane T's pre-registered rule (< 2.0 ⇒ closed at the encoding level) fired, on a 76-second pass over data already on disk.** The one structurally motivated interaction the record actually pointed at — a second-moment/dispersion term, which is a genuine off-diagonal coupling and was reached independently by two derivations — was then **measured** rather than assumed, by two lanes, and is **worse** (§7.1 #24–25). The lead's real question — *"does the off-diagonal structure encode information that materially changes the optimization problem?"* — is answered **no**, and answered before a Hamiltonian existed. The scope qualifier belongs with it: **stable rank 1.86 is the pair-distance space; in coordinate space it is 3.404/3.619 with k90 = 11.2**, above the threshold, so the verdict rests on the direct measurements and not on the rule that fired |
| **L10 — ADAPT-VQE** | **NOT PURSUED. Its precondition never arrived** | — (lane W audited the absence) | The lead is explicitly conditional: *"determine whether adaptive ansätze become meaningful **under a new formulation** ... whether a genuinely correlated structural Hamiltonian creates a meaningful adaptive growth process."* No such Hamiltonian survived. L5 closed by theorem (T1, T1b), its one surviving escape closed by measurement (M4, S30-L12), and L9's rank pre-check fired — so an adaptive pool had nothing correlated to grow into, and S26's product-circuit finding at α = 1 still governs the unchanged diagonal problem. **No CVaR-VQE was trained, no ansatz designed, and no trainability, gradient-variance or DLA measurement was taken** (`s30/QUANTUM_W.md` §0, with the absence verified seven ways). **The precise statement matters, because the obvious one is false:** a parameterised **9-qubit statevector circuit did run on all 126 targets**, with exact parameter-shift Jacobians and 300 Adam iterations per target — as a **regeneration check for the meter's ladder rungs**, asserted to 1e−6 against S28, on an ORACLE objective rather than a CVaR one. *"A report sentence saying 'no quantum compute was spent this sprint' would be false as written."* And that regeneration produced the sharpest available statement about where the barrier is **not**: the same circuit reaches **0.2516 Å** on the ORACLE objective and **3.4330 Å** on the deployed native-free one, against production's 3.2071 — **the ansatz can express the answer and the objective cannot point at it** |
| **L11 — Off-pool recognition** | **Pursued as the sprint's foundational question, and answered NO** | R, then G | The charter judged this *"may sit underneath several of the other leads"* and that *"a negative answer would be one of the most important results the project could produce."* It is negative, on an instrument built to remove the confound that invalidated the project's previous attempt (every rung an ideal-geometry backbone from leakage-safe Ramachandran torsions; kind, local realism and perturbation budget all matched; floor **0.347 Å**, not 0). **The verdict splits and the split is more informative than a flat null: ordering survives (DIS +0.347, +0.134 above its anchor control, p_max 0.000 over 43 channels), preference fails on all 43**, and the leave-fold-out combination that prefers a 0.55 Å structure to production on 93.0% prefers a **random pool member on 100%** and a **3 Å rung on 100%** — margin **−0.070 [−0.110, −0.028]**. Resolution bounds what any support rule could select on: DIS concordance in the near band is **0.520** at \|Δ\| = 0–0.25 Å, reaching 0.822 only above 4 Å — **coarse triage and nothing else** (S30-L19). Lane G then closed the class rather than the library: by **G1** every reflection-invariant single-structure channel **is** a distance-map reading, the only three references available are the three closed buckets, and the one family the theorem leaves open was built, is **genuinely exercised** (occupancy 0.69–0.97 against DIS's 0.23–0.35) and is **empty** (S30-L26). Lane R's own scope caveat travels with this: D1 says the information **is** in the global shape (held-out R² 0.911, ORACLE); what is missing is a **native-free globally-reaching channel**, which is *"a supply problem, not an impossibility"* |
| **L12 — Combinations and anything not listed** | **Pursued; the sprint invented two thirteenth leads and killed one of them itself** | all | The lead asks for *"the formulation that dissolves several at once rather than the one that patches each in turn"*, and names two candidate identifications: **L3 ≈ L11** and **L5 ≈ L8**. Both were confirmed and both dissolved downward rather than upward. L5/L8: lane T's **T1 + the value-of-a-bit law** answer both in one currency — the register question is settled and the selector question is not. L3/L11: lane R's D1 (**local −0.089 vs global +0.600**) and lane P's separation-bin split (**68% of the prize in \|i−j\| ≥ 7**) are *"two instruments, one statement"*, and lane G's G1 explains why lane R's library had no fourth kind to test. **The genuinely new lead the sprint invented and pursued** was the coordinator's: *what ρ does the best **combination** of the 21 fields reach?* — arithmetic that landed close enough to the 0.358 threshold to require measurement. Answer: **ORACLE global 0.1693, leave-fold-out 0.0948, fitted-21 0.0124**, and the reason is **rank, not count** — the fields span ~2 directions, so `√(Σρᵢ²) = 0.33 was the right formula applied to the wrong rank` (S30-L21). Lane T had registered `ρ_comb = ρ₀/√s₁` with the field count cancelling, and it predicted the measurement to within **0.008**. **The second new lead — the "two qubits" synthesis — the coordinator killed himself one minute after posting it** (§7.1 #57), and its surviving half is §7.4 #3 |

**Summary against the charter's own standard.** Six leads pursued to a verdict (L5, L6, L7, L8,
L11, and L3), four pursued inside another lane's question (L1, L2, L4, L12), two not pursued (L9,
L10) — one because the cheap check the lead itself demanded closed it, one because its precondition
never existed. The charter said a sprint that ignores ten and finds the mechanism is a success, and
one that executes all twelve and finds nothing is not. This sprint executed roughly ten and found
**no mechanism that moves the endpoint** — but it converted five of the twelve from open directions
into theorems, which is the outcome the charter's §18 ranks above *"a fragile 2.98 Å."*

---

## 9. Literature relied on and rejected

The charter asked for literature read *"for most of the sprint, not only at the start"* and to
**"read the equations, not the abstracts."** Lane L held the role permanently. The distinguishing
feature of this sprint's literature use is that **three papers closed project directions that no
experiment here had the power to close**, and one **contradicted a result we had measured** — which
turned out to be the most useful thing any of them did.

### 9.1 Relied on, and what each one settled

| source | what it settled here |
|---|---|
| **Blau & Michaeli — the perception–distortion tradeoff** | Predicts the *cross-kind* half of S28-L48's "20 of 31 scorers prefer production to a 0.25 Å ORACLE structure". Lane R used it to withdraw the project's most-cited negative **before its own numbers existed** (S30-L1). If perception–distortion predicts a result from the construction alone, the result is not evidence about nativeness |
| **Kennedy & O'Hagan (2001); Brynjarsdottir & O'Hagan, *Inverse Problems* 30:114007 (2014)** | Calibration parameter and discrepancy are **not jointly identifiable**, and a discrepancy term helps *only* given a strongly informative prior on its shape — with a *wrong* prior worse than none. This closes escape **E1** and is this project's own "confidently wrong costs 2–3× absent" arriving from a **second, independent literature** |
| **Maehara, *Oper. Res. Lett.* 43:526 (2015)** | The CVaR of a stochastic submodular set function **is not submodular**, and **no polynomial-time multiplicative approximation exists unless P = NP**. The sharpest available statement about the sprint's one hard constraint |
| **Wilder (AAAI 2018); Ohsaka & Yoshida (2017)** | The escape from Maehara is to stop asking for a single set and relax to a **portfolio — a distribution over sets** — restoring a `1 − 1/e` guarantee via continuous DR-submodular maximisation |
| **Nemhauser–Wolsey–Fisher; Das & Kempe** | Both guarantees (`1 − 1/e`, `1 − e^(−γ)`) require **monotone**. S29-L25 states in its own words that `V` is *"neither additive nor monotone"*. **Monotonicity, not submodularity, is what we fail first** — so every submodularity guarantee in the plan was void |
| **Goldberg (1984); Maurey's empirical method** | Fixed-`m` is the hard formulation (densest-`k`-subgraph); free-`m` is poly-time by max-flow. Maurey bounds any hull point within `R/√m`, which with the project's own dispersion `√63.82 = 7.99` gives **0.22 Å at m = 75** — i.e. the grid is not the constraint |
| **ANDIS (Yu et al. 2019)** | *"Native recognition and decoy discrimination cannot be optimized simultaneously with the same parameter sets."* An explicit statement, from the field, of the trade this project keeps rediscovering |
| **DOPE (Shen & Sali), incl. the authors' own DOPE-24 ablation** | The reference state's entire support for a 13-mer is [0, 15.05 Å] against a 15 Å table cutoff — the top 13% has **zero reference density**. The authors say plainly that DOPE *"is less accurate for smaller proteins"*; lane L turned that into the mechanism for the field's 40–50 residue wall |
| **Klenin & Langowski** | The writhe/Gauss double-integral formulation lane G used to build the only channel class that theorem G1 leaves open |
| **Abe et al., arXiv:2302.00704** | Already in the S29 index as the negative result on diversity interventions; used to reject diversity-aware selection without re-running it |

### 9.2 The one that contradicted us, which was worth more than the ones that agreed

*Improving consensus structure by eliminating averaging artifacts* (PMC2662860), 2090
non-homologous single-domain proteins under 200 residues:

```
baseline averaging (COMBO)   63.0 %  of atoms in clashes < 3.6 A
MCORE (their repair)          1.09 %                              <- a 58x reduction
PULCHRA                       3.64 %
RMSD, MCORE refined           3.36 A  against  3.28 A original    <- +0.08 A, WORSE
```

**Our E2 result (−0.022 Å) has the opposite sign from the published one at scale.** Lane L did not
explain the discrepancy away — it derived why both can be true: **the bonded fraction of all pairs
goes as ~2/n**, so a constraint repair is **~15% of the geometry at n = 13 and ~1% at n = 200**.
The published result is the large-`n` limit of ours.

> This is the sprint's model for how to use literature: a published number that **disagrees** with
> a measured one is a constraint on the mechanism, not a reason to doubt either. It also re-caps E2
> immediately after it escapes the identification invariance — perception–distortion catches what
> the non-identifiability theorem lets through.

### 9.3 Read and rejected, with reasons

| rejected | why |
|---|---|
| **DPPs / facility location / diversity-aware selection** | `V` is constant on centroid-equivalence classes, so **there is no diversity structure to exploit**. Rejected by algebra, not by trial |
| **Submodular maximisation guarantees** | Void — they require monotone, and `V` is not (above) |
| **QUBO / Ising formulations** | Available and exact if `f` were quadratic in the centroid, but **moot**: the continuous relaxation is cheap, and Frank–Wolfe on the simplex runs in seconds and **strictly upper-bounds any circuit on this objective** |
| **Portfolio/CVaR transfer from Wilder and Ohsaka & Yoshida** | **Mismatch stated rather than hidden:** their CVaR is over *exogenous* randomness; ours is over a distribution **the optimiser controls**, so the theorems do not transfer. What transfers is only the design lesson — relax the set, do not search it |
| **Fine-grained sub-region Ramachandran tables** | No sequence channel to exploit at peptide length; `phi` is predicted at 36.1° by full sequence context against 36.4° sequence-blind |

### 9.4 The reframe the literature produced

Lane L's reading, set beside the project's own numbers — 2 members with ORACLE **weights** reach
1.4315 Å against 75 members with ORACLE **membership** at 2.3055 Å, while choosing 2 of 500 costs
~17.9 bits against 7 for the top-128 argmin — gives:

> **Solving the set problem better is worth nothing. Exhaustive search already solved it
> (S29-L25) and the answer was worse than production. What is scarce is the information needed to
> SPECIFY a good set, and no solver supplies information.**

Which is the same sentence the sprint's measurements arrived at independently, from the other side.


## 10. The cost/RMSD meter's baselines for every cost tested

Built first, as §7 of the charter demanded, and extended mid-sprint with a `verify` verb and an
8-draw control. `python s30/s30_D_meter.py meter --f <name-or-module:function> --basis chain`.

### 10.1 The charter's four anchors all recompute

| anchor | charter's value | measured | basis |
|---|---|---|---|
| ladder rho (S28 ladder) | ~ -0.40 | **-0.4023** [-0.477, -0.322], 5/5 | built chain |
| cosine | ~ -0.03 | **-0.0339**, z = -2.45 vs its own null | CA cloud |
| native percentile | ~ 36.9 | **0.3676** [0.306, 0.405] | either |
| ORACLE preference | ~ 0.07 | **0.0714** | built chain |

Verified by `s30/s30_verify.py`, **36/36 matched, 0 mismatched, 0 missing**. The extension did not
move the instrument.

### 10.2 The shipped cost, on both bases

```
                          built chain      CA point cloud
ladder rho, S28            -0.4023          -0.1818       lowering the cost RAISES RMSD
ladder rho, CHARTER        -0.1964          +0.2603
ladder rho, FULL           -0.3187          +0.1179
cosine                      n/a             -0.0339       (-2.45 sigma, WRONG side)
native percentile           0.3676           0.3676       (0 = native is the pool's best; 0.5 = chance)
pref(ORACLE best vs PROD)   0.0714           0.2063
pref(pool member vs PROD)   0.0201           0.1265
production RMSD             3.2071           3.0483
GATE                        BLOCK            BLOCK
```

**Both bases BLOCK, and they block for different reasons** - the chain on the ladder alone, the
cloud on the ladder *and* the cosine. The gate's rule: *"a BLOCK on the ladder means the cost moves
in the wrong structural direction and gets no endpoint compute."*

### 10.3 The two costs whose percentiles were compared, and the mislabel I made

`native_pctile_DIS = 0.3676` and `native_pctile_DIS_SURR = 0.3688`. I quoted the **DIS_SURR** value
in lane D's brief beside five `DIS` numbers. They differ by 0.0012 so nothing downstream moved, but
the anchors block now carries both, named.

### 10.4 The cosine null, corrected from S29

Two different nulls, which S29 mixed:

- `per_draw_abs_mean = 0.1398` - the magnitude of **one** random direction on **one** target.
- `target_mean_null` - the null for the **126-target mean**, ~10x tighter: sd 0.0144, 95% interval
  [-0.0169, +0.0307].

The shipped cosine of **-0.0339** is z = **-2.45** against the second. Against the first it looks
unremarkable. **A cosine of +0.05 is far below 0.1398 and still several sigma above chance** - the
two nulls answer different questions and S29's paragraph used the wrong one.

### 10.5 What the meter could not see, which is a finding about the meter

Per contract rule 24, *a cost the meter cannot see is a finding about the meter.* Two surfaced:

1. **The meter was blind on the reporting basis.** It ran on the CA point cloud while the endpoint
   is the **built chain**. Both bases now run and, as §6.3 shows, they do not agree - the
   preference contrast is a *result* on one and *not a result* on the other.
2. **Its random-signed control was a single draw.** S29's seed-0 draw turned out to be the
   **maximum of its own eight**. The relative draw noise is **59% on the chain against 27% on CA**,
   because on the chain the cost prefers production to nearly everything, so the control is a rare
   event. The meter now takes `R = 8` draws and reports the per-draw distribution and the
   single-draw range.

> **A single-draw control is least trustworthy exactly where the effects are smallest - which is
> where this project's remaining effects live.**

### 10.6 Every cost in the library, on the endpoint basis

Charter item 27 asks for the meter's baselines *"for every serious cost function tested."* All **32**
costs the meter can address were run on the **built chain** — the basis the endpoint is reported on
— emitting **912 comparisons**.

**Read this as a descriptive baseline table, not a search for a winner.** The largest |z| is **8.06**
against an expected max-of-32 under the null of ~2.5–2.9, so the top rows are *not* order-statistic
artefacts — but what they measure is **coarse triage**, not what the pipeline needs (§10.6.2).

```
cost                    kind   ladder  lfold  contrast   xMDE   pctile  gate
LEG_steric             chain  +0.1553    3/5   +0.2515  +2.88      n/a  PROCEED
CAGEO                     ca  +0.1075    4/5   +0.2341  +2.49   0.6874  REVIEW
LEG_torsion            chain  +0.0502    3/5   +0.2212  +2.14      n/a  PROCEED
RAMA                   chain  +0.0760    2/5   +0.1498  +1.48      n/a  PROCEED
LEG_coop_helix         chain  -0.3828    5/5   +0.0645  +1.18      n/a  BLOCK
CONTACT                   ca  +0.1081    5/5   +0.1305  +1.09   0.4138  PROCEED
EXVOL                     ca  -0.1736    5/5   +0.0456  +0.93   0.5680  BLOCK
DISTPOT                   ca  -0.1849    5/5   +0.1002  +0.89   0.4321  BLOCK
ENV                       ca  -0.0520    3/5   +0.0918  +0.87   0.3932  PROCEED
CONTACT_LL                ca  -0.3197    5/5   +0.0655  +0.68   0.4025  BLOCK
DIS_SURR                  ca  -0.3890    5/5   +0.0367  +0.58   0.3688  BLOCK
DIS  (SHIPPED)            ca  -0.4023    5/5   +0.0357  +0.56   0.3676  BLOCK
LEG                    chain  -0.1560    4/5   +0.0546  +0.55      n/a  BLOCK
LEG_coop_sheet         chain  +0.4829    4/5   +0.0074  +0.47      n/a  PROCEED
DIS_MEAN                  ca  -0.4141    5/5   +0.0218  +0.36   0.3836  BLOCK
LEG_compactness        chain  +0.0604    3/5   +0.0441  +0.36      n/a  PROCEED
HP                        ca  -0.1265    4/5   +0.0273  +0.23   0.4644  BLOCK
TORS_CONS_POOL         chain  -0.1735    5/5   +0.0159  +0.17      n/a  BLOCK
LEG_hbond_longrange    chain  +0.4850    5/5   +0.0055  +0.15      n/a  PROCEED
RG_UNIV                   ca  -0.1430    5/5   +0.0169  +0.15   0.5151  BLOCK
RG_LAW                    ca  -0.0386    3/5   +0.0169  +0.14   0.5169  REVIEW
LEG_solvation          chain  -0.1330    4/5   -0.0179  -0.15      n/a  BLOCK
LEG_aromatic           chain  -0.0917    4/5   -0.0139  -0.17      n/a  PROCEED
LEG_contact            chain  -0.0157    2/5   -0.0253  -0.21      n/a  PROCEED
LEG_electrostatic      chain  -0.0544    3/5   -0.0327  -0.30      n/a  PROCEED
LEG_hbond_local        chain  -0.2310    4/5   -0.0263  -0.33      n/a  BLOCK
ELEC                   chain  -0.0468    3/5   -0.0417  -0.39      n/a  PROCEED
POOLGO_POOL               ca  -0.2552    5/5   -0.0714  -0.68   0.5144  BLOCK
DMAP_CONS_POOL            ca  -0.3784    5/5   -0.0933  -1.03   0.5865  BLOCK
SS_MATCH                  ca  -0.2248    5/5   -0.0967  -1.12   0.4440  BLOCK
CONS_POOL                 ca  -0.5633    5/5   -0.1290  -1.52   0.6175  BLOCK
DSSPHB                 chain  -0.2527    5/5   -0.1835  -1.93      n/a  BLOCK

ladder   = Spearman(cost, RMSD) over the S28 ORACLE ladder; POSITIVE is the good direction
contrast = pref(ORACLE best vs PROD) - pref(matched random-signed vs PROD), 8 draws
pctile   = the native's percentile in its own 500-pool (0 = best, 0.5 = chance)
```

### 10.6.1 A defect in this table, which the table's own signature caught

**The first version of this sweep printed `pctile = 0.0000` for sixteen costs**, making it look as
though the native were each one's argmin. It is not: **`native_pctile` is only computed for CA-kind
costs**, and my collector's `r.get('pctile') or 0` turned a **missing key into a fabricated zero**.

> Sixteen identical values across sixteen different costs is **exactly the signature that caught
> lane R's null-input artefact** (three identical 0.500s across three different questions,
> §Appendix A). I wrote the same bug three hours after recording that lesson in this report.
> **`or 0` on a possibly-absent numeric key manufactures data.** The column now reads `n/a`.

### 10.6.2 What the table says, and what it does not

**Three descriptive facts, none of which is a route:**

1. **The shipped cost ranks 12th of 32 on the preference contrast** (+0.0357, **0.56× MDE — below
   the not-a-result line**) and has **the second-worst ladder ρ in its own library** (−0.4023; only
   `CONS_POOL` is worse). *Lowering the shipped cost raises RMSD on the ORACLE ladder, and eleven
   other costs do better on a criterion it was never selected for.*
2. **Coarse triage is widespread.** Four costs clear 1× MDE on the contrast with a *positive*
   ladder ρ — `LEG_steric` (+2.88), `CAGEO` (+2.49), `LEG_torsion` (+2.14), `RAMA` (+1.48) — and
   `CONTACT` manages +1.09 with 5/5 ladder folds. Several are physics terms carrying no distogram
   information at all.
3. **And not one of them ranks the native inside its own pool.** Every computed percentile is
   0.37 or worse, i.e. at or below chance-adjacent, and **`CAGEO` — the second-strongest contrast
   in the library — puts the native at the 69th percentile, *worse* than chance.** The two
   abilities are **decoupled**.

**Why this is not a licence to swap the cost.** The contrast asks *"does this cost prefer a
near-native structure to a randomly displaced one?"* The pipeline needs *"does this cost order
candidates **within the top-75**?"* — which lane R measured directly on a kind-matched ladder and
found **preference fails on all 43 channels**. The fold agreement on the top rows is also weak
(`LEG_steric` 3/5 ladder folds, `LEG_torsion` 3/5, `RAMA` 2/5), and the binding stage is the filter,
not the score's coarse behaviour (§4.1).

> **The sweep reproduces lane R's split verdict across 32 costs on the endpoint basis: coarse
> ordering is common and cheap; in-band selection is absent everywhere.** That the shipped cost is
> near-worst on the coarse criterion while being the one selected for in-band use is consistent —
> and it is the clearest single illustration that *these are two different abilities and this
> project needs the one nothing has.*


## 11. What remains open

### 11.1 The one scientific question

**Does an observable exist whose error is incoherent with the pool's common mode?** Everything else
below is housekeeping by comparison. It has an admission test and a price (§12).

**The named gap in the closure, stated by the lane that produced it:** neither G1 nor lane P has
shown the *mutual information* between the native-free feature set and the prior's error is zero.
What is shown is that **the extractable part is the wrong part, by a theorem about which part is
identifiable.** That is a stronger claim than an exhausted search and a weaker one than an
information-theoretic zero, and it should be quoted as exactly that.

### 11.2 Open because it was never measured at the right scale

**Chiral functionals at 40+ residues.** G1 proves every achiral single-structure channel is a
distance-map reading, leaving chiral functionals as the only escape family; lane G built them and
found them empty. But **the theorem does not depend on chain length — only the emptiness does.** A
13-mer barely crosses itself, so writhe is dominated by local helical handedness. The family is
worth retesting where a chain can actually knot. *This project's instrument cannot ask the question:
the benchmark is 9–16mers and all 204 clusters are spent.*

**E2's restraint constant.** The AMBER relax at k = 30 has a whole-sample gain of **−0.0221 Å** (0.69% of baseline). What lane G confirmed at **−0.0406, 3.56× MDE, 5/5 folds** is the **dispersion contrast** — high-minus-low divergence within chain-length
tertiles — i.e. that the effect *concentrates*, not its size. **Both numbers describe the arm and they are different quantities; an earlier draft used the contrast as the gain.** It
**still has no native-free rule that selects its restraint constant.** It was picked on dev-set RMSD. Until that rule exists it is not deployable, independent
of its size.

### 11.3 Open defects, all located, none fixed here

| defect | location | status |
|---|---|---|
| A **withdrawn positive asserted in shipped code** — "the CVaR tail is worth +0.113 Å … the component's measured role" | `core/pipeline.py:821` | S25-L5 withdrew it and replaced it with −0.1405 Å at 0.68× MDE. `core/` was read-only for this audit. **Sixth instance of prose asserting a state that does not hold; first in shipped code** |
| **The cloud→chain projection seed is not pinned** | the multi-start projection | Two records give 3.2105 and 3.2126 from an *identical* cloud (§1.1). Harmless at the 3.00 Å scale; **fatal for any future sub-0.01 Å claim** |
| **The launch gate contradicts the governor** | `s26/jobrun.py:39` (`CPU_START = 85.0`, refused at `:137`) against `s26/governor.py:61` (`CPU_CEILING = 101.0`, CPU suspension deliberately **disabled**) | The governor is content up to 101% CPU; `jobrun` refuses to *launch* above 85%. During exactly the 94–95% the charter asks for, launches are blocked while the governor is happy, so lanes launch detached and bypass the governor. **Open since S29; this is the user's call, not mine** — the one-line fix is to raise `CPU_START` to match the governor's disabled-CPU policy and leave RAM as the binding constraint, which it is |
| **No sprint-wide multiplicity register** | — | Lanes counted their own comparisons; the sprint-wide total is a *lower bound* in the hundreds. A register written to as comparisons are emitted is a build item, not an exhortation |

### 11.4 Open because the sprint chose not to spend on it

The charter's twelve leads are a **leads register, not a task list**, and it says so. §8 records which
were pursued, which were not, and — where nobody touched one — the honest reason, without
retrofitted justifications.


## 12. The next highest-value scientific question

> ### Find an observable of *this molecule* whose error is **incoherent with the pool's common
> mode** — not merely decorrelated from the distogram.

This is not a restatement of "get a better prior". It is a specific, measurable, and cheaply
falsifiable condition that the sprint derived and then priced.

### 12.1 Why this and not something else

Three results compose into it, each independently established:

1. **The source enumeration collapses to two** — sequence and library — with physics as an
   *operator* on either rather than a third source (G1, repaired by lane P: a universal function
   evaluated on a target-specific argument yields target-specific output, so physics has no target
   argument of its own).
2. **Both sources are measured.** Sequence-derived features: out-of-fold R² **0.83%** against a
   registered 1.96% bar; the entire uncollapsed posterior is worth **−0.018**. The pool:
   genuinely informative at **+0.0758** out of fold at long range, and **strictly harmful when
   applied**, raising the residual's coherence with the pool's common mode from 0.6931 to **0.9172**.
3. **The mechanism says why, and it is a theorem, not a search result.** The predictable part of the
   prior's error *is* the common mode, and the common mode is exactly the component S30-L7 proves
   non-identifiable from within the pool.

> **What can be predicted is coherent and therefore harmful; what would help is incoherent and
> therefore unpredictable.**

### 12.2 The admission test — use this, not decorrelation

Define `coh` as the **within-target correlation of a candidate corrector's residual with the pool's
common-mode pair error** `μ_{t,p}` (the mean over the retained 75 of `d_{m,p} − d_nat,{t,p}`; the
S30-L7 quantity).

```
uncorrected                                     coh = 0.6931
every fitted corrector this project owns        coh = 0.783 / 0.786 / 0.917   <- ALL RAISE IT
imposed-structure ORACLE arms only              coh = 0.587 / 0.537           <- only these lower it

ADMIT a prior corrector iff it LOWERS coh below 0.6931.
```

**This is not "decorrelated from the distogram."** Lane L priced orthogonality at a **5.1% discount
on the requirement** — decorrelation is not the lever. Incoherence is a different condition and it
is the one that separates a −0.25 Å corrector from a +0.06 Å one **at identical accuracy**.

### 12.3 The price, which is the encouraging part

The required quality is **far below the displacement bound**, because this corrects the *prior* and
acts through the filter rather than being a displacement added to the answer:

```
R2 = 0.16  ->  -0.126 A        ORACLE-CONSTRUCTED, NOT DEPLOYABLE:
R2 = 0.24  ->  -0.247 A        a PRICE for a channel nobody has, never an achievement
```

And the target is sized: **five ORACLE signs on long-range pairs are worth −0.3259 Å on the built
chain (3.2126 → 2.8867, 2.05× MDE, 5/5 folds)** — which would clear the charter's primary target.
**The prize is five bits per target** — one sign per separation bin, `SEP_EDGES` =
[(2,2),(3,3),(4,4),(5,6),(7,99)]. And the deployable sign is **not merely weak, it is free**:
matched like-for-like it scores **0.600 against an always-positive baseline of 0.556 on all pairs,
and 0.643 against 0.627 at long range** — *above* the baseline in both spaces, but by a margin that
buys nothing. ≈ 0.8 is needed. Lane P's own wording is the right one: **essentially all of it is
free.** *(An earlier draft compared the all-pairs accuracy against the long-range-only baseline and
reported a shortfall that does not exist in either space.)*

### 12.4 The secondary recommendation, and three measurements agree on it

**If the next sprint has one place to spend bits, it is the readout.**

```
readout index bit      0.376 A/bit   (lane T's value-of-a-bit law at D = 3.05)
prior-sign bit         0.0713 A/bit  (lane P, cloud)  ->  5.3x
prior-sign bit         0.0652 A/bit  (lane P, chain)  ->  5.8x
candidate indexing over subset cardinality            ->  3x  (lane T, independent)
```

**Flagged unregistered, and the currencies are not exchangeable** — one names a candidate, one
shifts a prior, and lane T's law was fitted on an ORACLE ladder. But three independent measurements
agree in direction, which is more than any single one of them is worth.

### 12.5 What NOT to spend on

Closed by theorem, by price, **or by measurement** — and §7 keeps the three apart deliberately,
because a measurement-closed direction is the reopenable kind. With the closure named: the field combination (rank, not count);
sparse weighted readouts (argmin dominates at every bit budget); subset objectives through an
averaging readout (T1); the second-moment/quadric escape (twice, independently); generative spaces
(closed *jointly* with the readout); torsion encodings (48 bits against 7); common-mode correction
from pool data (non-identifiable at any K); achiral single-structure channels (G1); and
recognition from single-structure geometry (ordering survives, preference fails on all 43).


## 13. Architecture diagram of what was actually built

### 13.1 The production path — unchanged by this sprint

```
  sequence
     |
     v
  [1] BLOSUM62 retrieval, K = 500            <- NOT MEASURED as harmful at any stratum
     |                                          (and not exonerated either -- see 4.1)
     v
  [2] leave-fold-out ESM-2 650M distogram, 17 bins
     |                                          <- the ceiling lives here: its error is
     |                                             83.1% shape, and 68% of the recoverable
     v                                             prize is in |i-j| >= 7
  [3] L1 Bayes-risk score  ->  top-75 filter
     |                                          <- rho +0.6446 on the easy 108,
     |                                             +0.1066 to +0.3798 on the tail
     v
  [4] coordinate average of the 75            <- extracts spread (r = 0.885); cancels
     |                                           i.i.d. error, NOT coherent error
     v
     CA point cloud   3.0483 A                <- exact, reproduces to 4 decimals
     |
     v
  [5] multi-start ideal-geometry projection   <- SEED NOT PINNED: +-0.002 A (section 1.1)
     |
     v
     BUILT CHAIN      3.2105 A   <- THE ENDPOINT (charter section 2)


  [Q] quantum stage ....................... OFF THE PATH
      core/pipeline.py:179  quantum = False
      core/pipeline.py:241  PROD = Config()
      "VQE/CVaR do not participate at all"  (:173-178)
```

### 13.2 What S30 actually built

No new production architecture. **Nothing in §13.1 changed**, and that is the result, not an
omission — every candidate replacement was closed by measurement, by theorem or by price before it
reached the pipeline.

What was built is **instrumentation**:

```
  s30/s30_D_meter.py     the cost/RMSD meter, extended from s29's
                         + a `verify` verb        (the old `selftest` was an 8-residue synthetic
                                                   check that would pass on a drifted meter)
                         + BOTH reporting bases   (it had been blind on the built chain)
                         + an R = 8 draw control  (the single draw was the maximum of its own 8)
                         + a multiplicity counter (30 chain / 33 cloud comparisons per run)

  s30/s30_verify.py      36 headline numbers recomputed from artefacts, 0 mismatches
                         + the endpoint's own reproducibility
                         + an assertion that every path the ledger claims to have written EXISTS
                           (the mechanical fix for the S30-L0 failure)

  s30/THEORY.md          T1, T1b, the codebook reframe, the value-of-a-bit law
  s30/BRIEF.md           the charter, verbatim -- written 77 minutes after being claimed
  s30/AUDIT_V.md         the report adversary's findings (rule 28)
  s30/QUANTUM_W.md       items 12-19, audited against the code
```

### 13.3 The one structural change the sprint argues for

Not a new stage — a **gate on an existing one**:

```
  a candidate prior corrector
     |
     v
  compute coh = corr(residual, pool common-mode pair error)     [S30-L7 quantity]
     |
     +-- coh >= 0.6931  ->  REJECT.  It will emit WORSE at any accuracy.
     |                      (every corrector this project owns lands here)
     |
     +-- coh <  0.6931  ->  admit and measure the endpoint.
                            ORACLE arms reaching this: R2 0.16 -> -0.126 A (1.23x MDE, TYPE-M)
                                                       R2 0.24 -> -0.247 A (1.77x MDE)
                            NOT DEPLOYABLE -- a price, not an achievement. CA cloud basis.
```

**This gate is the sprint's deliverable — and it is ORACLE** (`s30_P_lr.py:58,78,202`: both
arguments of the correlation need `d_nat`). It is therefore a **development-time** test, run on the
126 labelled dev targets, **not** an inference-time screen. Within that scope it is cheap, it
discriminates a −0.25 Å arm from a +0.06 Å one at identical out-of-fold accuracy, and **none of the
three correctors built here passes it** — all three raise `coh` to 0.78–0.92.

*The two arms that do pass are ORACLE-constructed and not deployable; their gains are 1.23× and
1.77× MDE on the CA cloud, below the 2× this sprint treated as the survival standard, and the
artefact flags the smaller one as TYPE-M.*


---

## Appendix A — claims withdrawn during this sprint

The sprint's most distinctive feature was lanes destroying their own results, usually before anyone
asked. Recorded in full because the charter requires it and because it is the reason the surviving
results are worth anything.

### A.1 The coordinator's (mine)

| claim | how it died |
|---|---|
| **The "two qubits" synthesis** — a selecting readout over a wider register aimed at the tail | Written at 13:07, corrected at 13:08. Lane F closed filter width by ceiling (the ORACLE global argmin over k **is** the shipped 75) and showed the widening-rescues-the-tail result does not replicate on either filter-independent tail (≈0, and *reversed* at +0.73 on one) |
| **"The median target is worse than chance"** — reported to the user twice | `7 − log₂r` is right-skewed. Null mean 1.4050, null **median 0.9888**, null win rate 62.5%. "Below random on 82 of 126" is what a random ranking does to itself (null expects 78.8, z = +0.60). I read a median against a mean's baseline |
| **"LEG_torsion is at chance"** | My wording, not lane R's. Its anchor contrast is **+0.024 with the fold CI excluding zero** — it *does* order slightly above its control. It fails by being a quarter of the required margin, not by being noise |
| **The ceiling gate I handed lane X** — "measure the ceiling first, always" | An **anti-predictor**: 1 of 10 cells correct on predicting the endpoint's direction, against 5/6 for the set-mean law; mean residual 0.333 Å vs 0.0153 — **22×** |
| **"Nobody has built a typical-good generator"** | Somebody had, and **diversity is the correct design**. For an averaging terminal spread is the raw material; the most concentrated source ever built is the worst endpoint in the record |
| **"Run `selftest` and confirm it reproduces the baselines"** | `selftest` is a synthetic 8-residue check that runs in 0.4 s and **would pass on a meter whose every number had drifted**. The instruction was unsatisfiable. A real `verify` now exists |
| **The 0.3688 anchor in lane D's brief** | That is **DIS_SURR's**; the shipped cost is 0.3676. I listed it beside five `DIS` numbers |
| **"Find a source with decorrelated errors"** | Priced by lane L: a *perfectly orthogonal* channel is worth a **5.1% discount** on the requirement and must be **3.01× better** than anything owned. **Decorrelation is not the lever; skill is** |
| **Quoting stable rank 1.86 without its feature space** | 1.86 is **pair-distance** space; **coordinate space is 3.4–3.6** with k90 = 11.2 not 5.6 — and sparse readouts and second-moment constructions act on coordinates |
| **Citing S28-L48 as established** | Every rung in that ladder differs in **kind** as well as nativeness, so it measured a cross-kind preference that perception–distortion already predicts (lane R, S30-L1, posted before its own numbers existed) |

### A.2 Lane T's — four, three toward its own hypothesis

| claim | how it died |
|---|---|
| **P2b: retrieval worth +0.10–0.30 Å** | Measured **0.07**; its copula reasoning overestimated by 2–4× |
| **P4c: searched classes beat prefix by < 0.6 Å** | Measured **0.982** — missed by 64%, in the flattering direction |
| **"12 bits ≈ six coefficients — two independent routes"** | An **unregistered aside**; a coincidence, not a derivation. Withdrawn and corrected *downward* to 3.78 bits, and flagged louder precisely because it had no bar attached |
| **"HALFSPACE − PREFIX = −0.9816 Å, BETTER"** — its own headline | Lane Q named the null; lane T ran it with a common direction bank: **196% accounted**, split-half transfer 19%, and the **transferable rule lands +0.472 Å WORSE than production** (1.88× MDE, 5/5 folds). Its own diagnosis: the matched random-subset null asked whether a structured class beats an unstructured one, not whether **the winning direction is the same direction twice** |

### A.3 Lane L's

| claim | how it died |
|---|---|
| **EDM projection / triangle repair as a route** | Grepped the codebase after deriving it: S19 had already measured it and left a **named prohibition**. Refuted in the *unexpected* direction — an incoherent magnitude-matched field is worse on realisability and lands **1.24 Å better** |
| **The deployable half of its own reference-state proposal** | Its own algebra refutes it: the score decomposes additively, so the size-matched field **is** the deflated field, and deflation reallocates rank without creating any |
| **"The AMBER-relax benefit concentrates on divergent pools", 2:1** | Lowered to **roughly even** by lane L itself, on lane F's dispersion null (−0.128, CI includes zero) — adjacent evidence pointing the other way. **THIS WITHDRAWAL IS ITSELF WITHDRAWN (S30-L26): the claim is CONFIRMED at −0.0406, 3.56× MDE, 5/5 folds on a length-matched split, with 97.2% of the gain in the divergent half.** Lane L moved on a null about a *different outcome* and lane G's matching prior rested on a *different mechanism*; two revisions agreeing was not evidence. See §A.7 |
| **"Two independent arrivals" on confidently-wrong** | Downgraded to one literature result plus one project result on a *different object*, after lane F measured error×confidence (−0.729) as **worse** than error alone (−0.799) |

### A.4 Lane X's

| claim | how it died |
|---|---|
| **Its own S30-L10 admission condition** | Amended to be **stricter**: Δ(set mean) is two channels, not one, so `ADMIT iff Δ(bias B) + 0.298·Δ(set best) < 0` |
| **Its counting falsifier** | Volunteered as having **essentially tied its bar — 24.80% against 25%** — and therefore deciding nothing; the verdict rests on the transfer arm alone |
| **A −0.4744 correlation it could have quoted** | Flagged as **driven by one arm**; drop it and the value is +0.0851. "The strong negative is an artefact and must not be quoted" |

### A.5 Lane F's, lane D's, lane R's, lane Q's

| claim | how it died |
|---|---|
| **Lane F's F1c FAIL18 row** (+1.7674 Å) | The adversary: FAIL18 is *defined* as the targets whose top-75 retained **zero** in-band members, so **49% is forced arithmetic** and the rest is truncation-inflated. Lane F *had* checked circularity against production; the predicate is the **filter's own recall**. **"A matched control in the right space does not rescue a stratum defined by the outcome"** |
| **Lane F's own F2a and F2c** | Both refuted by its own measurements — no knee in filter width, and its widening result flagged by itself as circular |
| **Three defects in lane D's own combination arms** | An untuned ridge returning a negative cosine; a least-squares objective violating its own consistency floor (caught because ρ fell below the best single field, impossible for a true maximum); selection leakage worth 0.019 in ρ |
| **Lane R's 0.500/0.500/0.500 cell** | A **null-input artefact** — three identical values across three different questions is what caught it |
| **Lane Q's framing of my point 4** | "The framing survives; the operator does not" — identity is not only the better question but the **cheaper** one, so a sparse weighted readout is the expensive way to ask it |

### A.7 The one withdrawal that was itself wrong, and lane G's self-kill

**A withdrawal is a claim too, and it can be wrong in the same ways.** Appendix A recorded lane L
lowering "the AMBER-relax benefit concentrates on divergent pools" from 2:1 to roughly even. Lane G
then measured it and **confirmed the concentration**: the high-minus-low dispersion contrast within
chain-length tertiles is **−0.0406 Å, 3.56× MDE, 5/5 folds**, with **97.2%** of the gain in the
high-dispersion half and a permutation null at p = 0.000. *(The arm's own whole-sample gain is
−0.0221 Å — the contrast measures where the gain sits, not how big it is.)* Lane G had
*also* registered 2:1 against, and lost.

> **Two lanes agreeing did not make a prior.** Lane L moved on lane F's null, which was about a
> different *outcome* (filter failure, not relax gain). Lane G moved on the common-mode argument,
> which was about a different *mechanism* (bias repair, not geometry repair). Each rested on an
> object that was not the one under test, and their agreement read as convergence.

What survives the confirmation is narrower than the claim: **the effect is too small to matter,
which is not the same as failing to reach the tail, and my first draft argued it on the wrong
stratum.**

I originally quoted two FAIL18 rows. **FAIL18 is the wrong tail for a dispersion-graded question**
— it is the *filter-defined* stratum and is barely enriched in high dispersion (11 of 18 against
9.0 expected) — and `s30/results/s30_G_disp2.json` holds the filter-independent rows **right beside
the ones I quoted**, where both legs reverse:

```
                                  FAIL18 (quoted)        worst18 by pool mean (not quoted)
tail dispersion vs the rest       +0.4583  0.47x MDE     +1.6116  1.44x MDE  CI [+0.583,+2.076]
relax gain, tail vs rest          -0.0113 vs -0.0239     -0.0386  vs  -0.0193
```

On the filter-independent tail the relax gain is **twice as large on the tail**, not half. The
ledger printed both; the report printed one. Worse, the FAIL18 relax-gain difference I called *"the
wrong direction"* is **0.22× MDE** — I used a number below my own not-a-result line as evidence of
absence.

**The third leg is sound and is sufficient on its own:** targeting the divergent half buys
**−0.0215 against −0.0221**, i.e. nothing, because the other half contributed nothing to begin
with. And the whole arm is **0.69% of the 3.2126 Å baseline**. So **E2 stays out of the mechanism
column on SIZE, not on reach**, and its restraint constant still has no native-free selection
rule.

**And lane G killed its own positive before anyone quoted it.** WRITHE's preference contrast is
**+0.1641, 2.17× MDE, 5/5 folds, max-null p = 0.000**, clearing *both* of lane R's preference
clauses that none of 43 channels cleared. It is **cross-kind** — the control keeps deposited
coordinates while the near rungs are ideal rebuilds. The kind-matched statistic settles it: the
rebuilt native's percentile inside its own ladder is **0.6061 for WRITHE and 0.6732 for |WRITHE|,
both worse than the 0.5 chance line**, against DIS's 0.2876. **Second instance of the cross-kind confound this sprint**, after S28-L48 (withdrawn by lane R).
*My first draft said "third" and counted the widening result, which Appendix A correctly records as
circular rather than cross-kind — a different defect.*

### A.6 The three checklist entries this sprint earned

1. **A matched control in the right space does not rescue a stratum defined by the outcome.**
2. **When a statistic is a nonlinear transform, its mean, median and win-rate have three different
   nulls** — never read one against another's baseline.
3. **"The winning direction is the same direction twice"** is a failure mode distinct from the
   order-statistic one, and needs a common direction bank to detect.
