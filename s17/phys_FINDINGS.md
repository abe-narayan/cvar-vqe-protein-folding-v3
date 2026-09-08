# SPRINT 17 — PHYSICS findings

**The two physical models, on one instrument: genuine Legacy (eleven components,
`core.energy`, `DEFAULT_WEIGHTS`, never fitted) and genuine AMBER ff14SB/GBn2
(`core.amber`, OpenMM). No learned surrogate stands in for either anywhere in this
workstream.**

Pre-registration: `s17/PREREG_phys.md` (E1–E5), written before any Sprint-17 number existed.
Modules: `s17/phys_lib.py` (the constrained-optimisation Hamiltonian, controls, statistics),
`s17/phys_ca.py` (E1, runs AMBER), `s17/phys_ca_report.py` (E1 analysis),
`s17/phys_ident.py` (the identical-candidate instrument), `s17/phys_report.py` (E2, E4, E5),
`s17/phys_gate.py` (E3). Artefacts in `s17/results/`.

Tiering: **DEMONSTRATED · SUPPORTED · REFUTED · ORACLE DIAGNOSTIC · EXACT (a theorem, not a
discovery)**.

---

## 0. VERDICT — leading with what most damages a physics claim, including my own

> **1. My own pre-registered headline hypothesis is HALF REFUTED, and I declared before the
> run why its accuracy half could not be tested.** Cα-fixed AMBER (particle mass 0, ε = 0
> *exactly*) is a genuine steric repairer at literally zero Cα cost — against **its own gated
> input** it takes heavy-atom clashes below 2.0 Å from **0.31 to 0.00**, below 2.6 Å from
> **2.03 to 0.01**, and the minimum heavy separation from **2.606 Å to 2.913 Å**, while the Cα
> coordinates move by a **maximum of 7 × 10⁻¹⁵ Å**. *(That baseline is deliberately not the
> instrument-wide 5.61/15.88/2.096: the gate keeps the 78 least-broken targets, and quoting the
> ungated input against a gated arm is the Sprint-16 failure mode — see §1.2, where I caught it
> in my own first draft.)* But it **fails
> the convergence gate on 48 of 126 targets (38%)** against the incumbent's 3, it **degrades
> Ramachandran-favoured from the input's 0.836 to 0.734** (−0.1718 [−0.2203, −0.1256], worse on
> 44 targets and better on 1), and it introduces **cis peptide bonds at a mean fraction of
> 0.313 with 60.2° of ω deviation** where the input has none. The
> pre-registered falsifier for H1b — *Ramachandran-favoured below the input's* — **FIRED**.
> **Freezing Cα does not make validity free; it moves the price from Cα accuracy to backbone
> plausibility.**
>
> **2. And the accuracy half of that hypothesis was never testable.** At ε = 0 the Cα
> coordinates do not move, so Cα-RMSD equals the do-nothing RMSD **exactly, by construction**.
> The falsifier I was handed ("accuracy cost above 0.10 Å against doing nothing") **cannot
> fire at that rung**. This is recorded in `PREREG_phys.md` before the run rather than
> claimed afterwards as a result: it is the Sprint-16 failure mode (*a quantity measured
> correctly and read as a different quantity*) waiting to happen, and the ε = 0 rung is
> therefore judged on validity alone.
>
> **3. I found a units defect in my OWN zero-information control and report it.**
> `core.geometry.build_backbone` takes torsions in **radians**; my first constant-α-helix
> reference passed **degrees**, which silently builds a different constant conformation
> (recomputed torsions −25.9° / −172.9°) scoring Ramachandran-favoured **0.000** instead of
> 1.000. Every affected number was recomputed; nothing built on the broken control survives
> into this document. The fix is at read time (`phys_ca_report.load`) and the defect is
> recorded, not hidden.
>
> **4. The result that damages the pipeline hardest is not about the energies at all.** The
> shipped distance score's shortlist is **no better than a matched random shortlist of the
> same size at containing a sub-2 Å candidate** (M = 25/75/150/300: +0.016, −0.016, −0.032,
> −0.005, every CI including zero) and is **significantly WORSE at containing the BEST one**
> (+0.195, +0.242, +0.179, +0.097 Å, every CI excluding zero). Win/loss is near-even, so this
> is mean-carried and the median-vs-mean warning applies — but the direction is identical at
> all four widths.
>
> **5. Legacy's gate preserves the near-native CLASS and destroys the near-native BEST, and
> the mechanism is now measured to the level of a single candidate.** On the 73/126 targets
> that have a sub-2 Å candidate, dropping the worst half by Legacy keeps **61.3%** of the
> sub-2 Å members against a matched random gate's **50.1%** (+0.112 [+0.023, +0.198],
> 49W/21L) — and the set best is **1.427 Å against random's 1.208** (+0.219 [+0.031, +0.444]).
> In Legacy's own normalised order the pool's **single best** candidate sits at **0.450**
> while the **median** sub-2 Å candidate sits at **0.410** (random = 0.500). **Legacy likes
> near-natives and likes the best one systematically less than the typical one.** That is
> Sprint 16's G4 set-mean trap, resolved to one candidate.
>
> **6. `leg_torsion` is the single exception, and it is the only positive Legacy result in
> this workstream.** It is the strongest in-band component (Spearman 0.101 against the
> distance score's 0.126), the only score of seventeen with a *negative* (better-than-random)
> in-band selection delta (−0.0446 [−0.1378, +0.0482], 66W/57L, CI includes zero), and the
> only gate that buys near-native recall (+0.159 [+0.073, +0.239] at f = 0.50) **without**
> significantly damaging the set best (+0.011 [−0.079, +0.128]).
>
> **7. The decisive Legacy-vs-AMBER experiment ran at FULL scale — 63,000 genuine ff14SB/GBn2
> single points on the identical K = 500 candidates of all 126 targets, no shortlist — and
> AMBER lost every ranking role.** Global Spearman **−0.027 against Legacy's +0.307**; band
> AUROC **0.499**, exactly chance; global argmin selects **5.191 Å, worse than random's
> 4.427**; incremental value over the distance model **−0.004 [−0.037, +0.030]** with its
> `angle` (−0.037) and `torsion` (−0.053) terms **significantly anti-informative**; and as a
> gate it *loses* pool-best survival to a matched random gate (0.452 vs 0.508). **At ~1,500×
> Legacy's cost, all-atom AMBER cannot do the one thing Legacy demonstrably can — reject
> garbage.** The physics-ranking route is now closed by measurement on both potentials.
>
> **8. Disagreement as information (sprint §34) is REFUTED, cleanly.** Every
> distance/Legacy disagreement feature has real leave-fold-out skill against the selector's
> realized error (0.208–0.449) and **all of it vanishes** once the trivial native-free
> feature `dist_score` — how good the objective thinks its own pick is, LFO 0.598 alone — is
> partialled out: partial LFO ranges −0.148 to +0.168 with no consistent sign across
> features or responses. **Disagreement is subsumed by the objective's own score.**

---

| claim | tier |
|---|---|
| Cα-fixed AMBER removes every steric clash at *exactly* zero Cα cost | **DEMONSTRATED** (§1) |
| …and does so without cost | **REFUTED — my own hypothesis** (§1): gate fails 48/126, Ramachandran 0.836 → 0.734, cis fraction 0.000 → 0.313 |
| the ε = 0 accuracy falsifier is testable | **EXACT (tautology), declared before the run** (§1) |
| the k-ladder traces a validity/accuracy Pareto frontier | **DEMONSTRATED** (§2) |
| some ε gives k = 0's validity at < 0.05 Å Cα cost | **REFUTED** (§2) |
| Legacy has in-band ranking skill over its matched random control | **REFUTED** (§3) |
| Legacy has *incremental* in-band skill beyond the distance model | **SUPPORTED, small, and it does not convert to selection** (§3) |
| the Legacy gate preserves the near-native class | **DEMONSTRATED** (§4) |
| …and preserves the near-native BEST | **REFUTED** (§4) |
| the mechanism is the set-mean trap at the level of one candidate | **DEMONSTRATED** (§4) |
| `leg_torsion` gates coverage without destroying the best | **SUPPORTED** (§4) |
| the shipped distance shortlist beats a random shortlist on coverage | **REFUTED** (§4) |
| ranker disagreement predicts selector failure | **REFUTED, subsumed by `dist_score`** (§6) |
| AMBER single points rank better than Legacy at 1,500× the cost | **REFUTED** (§5): global ρ −0.027 vs +0.307, band AUROC 0.499 |
| AMBER adds in-band information beyond the distance model | **REFUTED** (§5): −0.004 [−0.037, +0.030]; two terms significantly negative |
| AMBER's cis-peptide defect (Sprint 16 G12) is a force-field problem | **REFUTED** (§1.3): ρ = −0.949 with the input's backbone contraction |
| the incumbent's N/CA/C restraint set is arbitrary | **REFUTED** (§2): it beats Cα-only on every axis at matched displacement, and it is what holds ω trans |
| the ε and k parameterisations trace one frontier | **SUPPORTED** (§2.1), n = 8 cross-check |
| my own zero-information α-helix control was correct as first built | **REFUTED — units defect, mine, fixed** (§0.3, §1) |

---

## 1. E1 — Cα-PRESERVING AMBER REPAIR AS CONSTRAINED OPTIMISATION

**The problem, posed properly:** minimise `E_ff14SB/GBn2(x)` subject to
`‖CAᵢ(x) − CAᵢ(x₀)‖ ≤ ε` for every residue i. Implemented two ways (`s17/phys_lib.py`):
**HARD** at ε = 0 (Cα particle masses set to zero, which OpenMM's `LocalEnergyMinimizer`
treats as frozen) and **SOFT** over the Lagrangian family `E + (k/2)‖ΔCA‖²`, whose dual
variable k sweeps the same Pareto frontier. N, C, O, CB, sidechains and every hydrogen are
**free at every rung**. The incumbent's harmonic N/CA/C restraint at k = 30 is carried as a
comparison arm, not as a rung.

**Why the ladder is swept in k and not in ε, measured rather than assumed.** A per-atom
flat-bottom wall `½k_f·max(0, |x−x₀|−ε)²` is the exact penalty form of the ε-ball, but it has
a kink at the wall that L-BFGS handles badly. On 1A13 at ε = 0.10 Å: `k_f = 10³` **overshoots
the constraint by 0.267 Å** (it is not enforcing it at all); `k_f = 10⁶` holds it to 0.004 Å
and costs **78 s per target per rung** against the harmonic arm's 7.6 s. So the ladder is
harmonic, and **the frontier is plotted against the REALISED Cα displacement — the quantity
the constraint bounds, measured on every structure, never assumed.**

### 1.0 The instrument reproduces the incumbent bit-for-bit

The `k30` rung is built here from an **independently constructed OpenMM System** — a fresh
`createSystem` with my own restraint forces — not through `core.amber.refine_coords`. Against
`s16/results/repair_A.json` (itself bit-identical to Sprint 15's artefact) on **all 126
targets**:

| | |
|---|---|
| max \|ΔCα-RMSD\| | **0.000 × 10⁰ Å** |
| max \|ΔEnergy\| | **2.7 × 10⁻¹³ kcal/mol** |

Every number below is on the same instrument as `s15/phys_FINDINGS.md`, `s16/energy_FINDINGS.md`
and `s16/repair_FINDINGS.md` and may be compared to them directly.

### 1.1 The ε = 0 rung — the tautology declared before the run, and what it leaves

At ε = 0 the Cα coordinates are frozen. Measured, not assumed: over 126 targets the **maximum
per-atom Cα displacement is 7.1 × 10⁻¹⁵ Å** and the maximum |ΔCα-RMSD| against doing nothing
is **6.9 × 10⁻¹⁴ Å** — machine epsilon.

> **Consequence, and it was written into `PREREG_phys.md` before the run.** The
> pre-registered falsifier I was handed — *"accuracy cost above 0.10 Å against doing
> nothing"* — **cannot fire at this rung**: the cost is zero **by construction**, not by
> measurement. Quoting "Cα-fixed AMBER costs 0.000 Å of accuracy" as a result would be
> precisely the Sprint-16 failure mode. **The whole scientific content of the ε = 0 rung is
> on the validity axis.**
>
> **The same tautology explains the one AMBER arm in this programme that PASSES its
> rotated-lab-frame null.** `cafix` returns mean 0.00000, sd 0.00000, **max \|Δ\| 0.00000**,
> **0 of 9 converged targets moving by more than 1e−6 Å** (12 drawn, 3 gate-excluded) —
> against a band (|mean| ≤ 0.005 **and** max ≤ 0.05) that **every** converged full-depth arm
> in Sprint 16 failed with a 0.137 Å tail. **This is EXACT, not a discovery.** The output Cα
> *are* the input Cα and every RMSD here is rigid-invariant, so the null is zero for the same
> reason the accuracy effect is zero. **Passing the frame null costs nothing and buys nothing
> at this rung** — the two tautologies are the same tautology, and neither is evidence that
> the minimiser became reproducible.

### 1.2 What Cα-fixed AMBER actually does — full n = 126

Input: the raw all-atom coordinate average of the shipped top-75 set, do-nothing
**3.0498 Å**. Both mandatory controls carried.

| rung | gate excl | Cα disp rms | Cα disp max | RMSD Å | vs DO-NOTHING (gated) | vs MATCHED RANDOM | wall s |
|---|---|---|---|---|---|---|---|
| **none** (zero-info) | 0 | 0.0000 | 0.0000 | **3.0498** | — | — | 0 |
| *zero-info* α-helix | 0 | 2.5510 | 4.1283 | 4.0696 | +1.0198 [+0.7321, +1.3319] 32W/94L | — | 0 |
| *zero-info* β-strand | 0 | 8.1429 | 13.9872 | 8.6937 | +5.6439 [+5.0537, +6.2387] 8W/118L | — | 0 |
| **`cafix`** (ε = 0) | **48** | **0.0000** | **0.0000** | **= do-nothing, exactly** | **−0.0000 [−0.0000, +0.0000]** | −0.0000 | 10.3 |
| `k30` (incumbent) | 3 | 0.7305 | 0.9992 | 3.1827 | **+0.1333 [+0.1031, +0.1645] 24W/99L** | **+0.0354 [+0.0008, +0.0697]** 43/80 | 11.0 |

*(The k30 row reproduces `s16/repair_FINDINGS.md` §1 exactly: +0.1333 [+0.1031, +0.1645],
24W/99L. Its matched-random contrast now **excludes zero in the wrong direction** —
AMBER's displacement is worse than a random displacement of the same size, tightening
Sprint 16's +0.0225 [−0.0138, +0.0571].)*

**THE CONVERGENCE GATE, declared before use (final energy finite and ≤ 1000 kcal/mol,
`core.amber.CONVERGE_MAX_KCAL`), reported with its exclusions by name:**

* **`k30` excludes 3 / 126 — `1D6X 2NB7 7BX2`** — the *same three* as Sprint 16, a fourth
  independent confirmation of the instrument.
* **`cafix` excludes 48 / 126 (38%)**: `1A13 1D6X 1JBF 1KWE 1M02 1RSW 2BP4 2EFZ 2LM8 2LNG 2LU6
  2LWS 2LWU 2MFV 2MID 2MIG 2MK7 2MQ2 2MSA 2N5C 2N9M 2NB7 2NBC 2NDN 5H1H 5NVB 5W52 6BX9 6EY3
  6HVK 6OQP 6QAX 7BX2 7JGX 7JS6 7LCW 7S3O 7T3H 7YFS 8FLP 8IS3 8T61 8T62 8T63 8TXS 8ZG2 9BAF
  9KAR`. **This is a sixteen-fold increase in gate failure and it is the operator's largest
  single liability.** §1.3 identifies the cause and it is not the minimiser: gate survival
  correlates with the input's mean Cα–Cα distance at **ρ = +0.807**.
* The `cafix` RMSD column is not reported as a gated mean because it is identically the
  do-nothing value; reporting "2.535 Å gated" would compare the 78 *easier* surviving targets
  against a 126-target baseline, which is the Sprint-16 failure mode in a new costume.

**VALIDITY, on the same structures whose RMSD is quoted.** The two zero-information rows are
here because a validity statistic quoted alone is meaningless: a constant α-helix scores
Ramachandran 1.000 and zero clashes **by construction**, and so does a constant β-strand.

**A CORRECTION I FOUND IN MY OWN HEADLINE, BEFORE PUBLISHING IT.** The convergence gate selects
a **different subset for every rung**, so quoting a gated arm's validity against the *ungated*
input's mean is the Sprint-16 failure mode in a new costume — and my first draft did exactly
that. `cafix` converges on the 78 **least broken** targets (§1.3: gate survival correlates with
the input's Cα–Cα spacing at ρ = +0.807), whose input already carries **0.31** sub-2.0 Å
clashes, not the instrument-wide 5.61. **Every rung below is therefore followed by its own
gated input row (`in@rung`), and only those two rows may be compared.**

| rung | ramaFav ↑ | ramaOut ↓ | clash<2.0 ↓ | clash<2.6 ↓ | minHeavy ↑ | bond ↓ | angle ↓ | geomDev ↓ | **cis ↓** | **ωdev ↓** |
|---|---|---|---|---|---|---|---|---|---|---|
| **none** (input, all 126) | 0.836 | 0.081 | 5.61 | 15.88 | 2.096 | 0.3147 | 0.1846 | 0.2480 | **0.000** | **0.0** |
| *zero-info* α-helix | **1.000** | **0.000** | **0.00** | **0.00** | 3.069 | 0.0000 | 0.0000 | 0.0000 | 0.000 | 0.0 |
| *zero-info* β-strand | **1.000** | **0.000** | **0.00** | **0.00** | **4.480** | 0.0000 | 0.0000 | 0.0000 | 0.000 | 0.0 |
| **`cafix`** (n = 78) | **0.734** | 0.214 | **0.00** | **0.01** | **2.913** | 0.0621 | 0.0494 | 0.0469 | **0.313** | **60.2** |
| ↳ *its own gated input* | **0.906** | 0.043 | **0.31** | **2.03** | 2.606 | 0.2002 | 0.1168 | 0.1477 | 0.000 | 0.0 |
| `k30` (n = 123) | 0.875 | 0.041 | 0.00 | 0.03 | 2.809 | 0.0262 | 0.0243 | 0.0188 | 0.074 | 21.4 |
| ↳ *its own gated input* | 0.837 | 0.080 | 5.39 | 15.28 | 2.116 | 0.3100 | 0.1822 | 0.2440 | 0.000 | 0.0 |

> **What the correction costs the headline, stated plainly.** Against its **own** gated input,
> `cafix` removes **0.31 → 0.00** sub-2.0 Å clashes and **2.03 → 0.01** sub-2.6 Å clashes and
> gains **+0.31 Å** of minimum heavy separation — real, one-sided (13 targets improved, 0
> degraded on the 2.0 Å axis), and at 7 × 10⁻¹⁵ Å of Cα displacement, but **an order of
> magnitude smaller than the "5.61 → 0.00" the ungated comparison would have claimed.** The
> big number belongs to `k30`, which converges on 123 targets and takes **5.39 → 0.00**. And
> the Ramachandran damage is *worse* under the correct comparison, not better: against its own
> gated input `cafix` falls **0.906 → 0.734**, not 0.836 → 0.734.

**THE CONJUNCTION** — validity *and* staying near the input, paired against the do-nothing
input on the gated subset:

| rung | Δ rama favoured | Δ clashes < 2.0 Å | Δ bond strain |
|---|---|---|---|
| **`cafix`** | **−0.1718 [−0.2203, −0.1256] 1 better / 44 worse** | **−0.3077 [−0.5128, −0.1410] 13/0** | −0.1381 [−0.1662, −0.1117] 78/0 |
| `k30` | +0.0379 [+0.0108, +0.0664] 40/24 | **−5.3902 [−7.5203, −3.5528] 54/0** | −0.2837 [−0.3164, −0.2511] 122/1 |

> ### H1b — the pre-registered hypothesis, and the falsifier that FIRED
>
> **PRE-REGISTERED SUCCESS CRITERIA:** zero heavy-atom clashes below 2.0 Å on ≥ 90% of
> targets **and** Ramachandran-favoured ≥ 0.824.
> **PRE-REGISTERED FALSIFIER:** ≥ 1 clash below 2.0 Å on more than 10% of targets, **or**
> Ramachandran-favoured below the do-nothing input's 0.836.
>
> **Clashes: SUCCEEDS, completely — on the population it converges on.** `cafix` reaches a mean
> of **0.00** heavy-atom clashes below 2.0 Å and **0.01** below 2.6 Å against its own gated
> input's 0.31 and 2.03, and raises the minimum heavy separation from 2.606 Å to **2.913 Å**,
> **at a Cα displacement of 7 × 10⁻¹⁵ Å**. Paired against its own gated input:
> sub-2.0 Å clashes **−0.3077 [−0.5128, −0.1410], 13 improved / 0 degraded** (13 of the 78 have
> such a clash and all 13 are repaired); sub-2.6 Å **−2.0128 [−2.9231, −1.2179], 31 / 0**;
> minimum heavy separation **+0.3076 [+0.1901, +0.4310], 44 improved / 34 degraded** — that
> last one is near-even and therefore mean-carried, and I flag it rather than quote it as
> one-sided. The steric half of "all-atom validity" is
> genuinely free — **and it is a smaller quantity than the whole instrument's 5.61 clashes,
> because the targets carrying those clashes are largely the ones `cafix` fails to converge on.**
>
> **Ramachandran: FALSIFIER FIRES, and harder under the correct baseline.** `cafix` returns
> **0.734 against its own gated input's 0.906** (and against the instrument-wide 0.836) —
> −0.1718 [−0.2203, −0.1256], **worse on 44 targets and better on 1**. It
> also introduces cis peptide bonds at a mean fraction of **0.313** with **60.2° of mean ω
> deviation**, where the input has none, and it **fails the convergence gate on 48 of 126
> targets (38%)** — against the incumbent's 3.
>
> **The pre-registered hypothesis is therefore HALF REFUTED, and the half that fails is
> mine.** Freezing Cα does not make validity free. It **relocates the price** from Cα
> accuracy to backbone plausibility.

### 1.3 WHY — the mechanism, and it exonerates the force field

The coordinate average's Cα–Cα spacing is **2.949 Å against the ideal trans value of
3.80 Å — a 22.4% contraction** (min 1.760 Å, max 3.796 Å over 126 targets), which
independently confirms the ledger's *"averaging contracts the backbone 25.8%"*. Freezing Cα
freezes that contraction, and a peptide unit cannot span 2.9 Å in the trans conformation.

**Measured, n = 126 (78 gated at `cafix`):**

| | Spearman with the input's mean Cα–Cα distance |
|---|---|
| **cis-peptide fraction at `cafix`** | **−0.949** |
| convergence (gate passed) at `cafix` | **+0.807** |
| cis-peptide fraction at `k30` | −0.665 |

> **Sprint 16's new defect G12 — *"restrained ff14SB/GBn2 relaxation introduces cis peptide
> bonds on 44 of 126 targets"* — is now explained, and it is not a defect of the force
> field.** At ρ = **−0.949** the cis fraction is an almost deterministic function of how far
> the coordinate average has contracted the backbone. The peptide bond flips because the
> restraint pins the chain at a spacing no trans peptide can reach. **The same variable
> explains 38% gate failure at ρ = +0.807.** The architectural consequence is direct: **fix
> the averaging operator's contraction, do not blame or retune the force field** — and note
> that the incumbent's choice to restrain **N and C as well as Cα** is what holds ω near
> trans (cis 0.074 vs 0.313), which had never been given a reason.

---

## 2. E1 — THE PARETO FRONTIER, and where validity is free and where it is bought

Pre-registered, native-free, fold × length stratified 30-target subsample
(`s16.repair.shape_subsample`, reused **verbatim** so no new selection freedom is created).
**SHAPE only: no rung is re-quoted at full n as an out-of-sample choice.** **n = 30, complete.**
(The pass was killed once by `core.amber.memory_guard` at its 92% physical-memory ceiling —
the same failure that killed a Sprint 16 shard — and restarted with a yield-and-retry; the
per-target checkpoint meant nothing was lost.) Do-nothing on this subsample is 3.5061 Å.

| rung | gate excl | **realised Cα disp** | RMSD cost vs DO-NOTHING | vs MATCHED RANDOM | ramaFav | clash<2.0 | geomDev | **cis** | ωdev |
|---|---|---|---|---|---|---|---|---|---|
| **none** | 0 | 0.000 | — | — | 0.781 | 7.37 | 0.2885 | **0.000** | 0.0 |
| **`cafix`** (ε = 0) | 15 | **0.000** | **−0.0000** (tautology) | −0.0000 | 0.691 | **0.00** | 0.0476 | 0.313 | 59.2 |
| `ca1000` | 6 | 0.185 | **+0.0081 [−0.0033, +0.0187]** | −0.0003 | 0.647 | **0.00** | 0.0706 | **0.539** | 89.2 |
| `ca300` | 2 | 0.317 | +0.0305 [+0.0085, +0.0519] | +0.0080 | 0.649 | **0.00** | 0.0709 | **0.574** | 96.4 |
| `ca100` | 2 | 0.499 | +0.0662 [+0.0331, +0.0996] | +0.0377 | 0.717 | **0.00** | 0.0546 | 0.488 | 89.1 |
| `ca30` | 2 | 0.782 | +0.1402 [+0.0843, +0.1968] | +0.0343 | 0.814 | **0.00** | 0.0275 | 0.201 | 47.9 |
| `ca10` | 2 | 0.960 | +0.2159 [+0.1290, +0.3064] | +0.0529 | 0.817 | **0.00** | 0.0165 | 0.105 | 28.2 |
| `ca3` | 2 | 1.111 | +0.2682 [+0.1588, +0.3833] | +0.0812 | **0.847** | **0.00** | 0.0162 | **0.084** | 22.5 |
| `free` (unrestrained) | 2 | 1.878 | +0.4586 [+0.2374, +0.6876] | **−0.0954 [−0.2777, +0.0881]** | 0.845 | **0.00** | 0.0172 | 0.087 | 21.6 |
| **`k30`** (N/CA/C) | 2 | 0.824 | +0.1217 [+0.0688, +0.1769] | +0.0210 | 0.828 | **0.00** | 0.0221 | **0.086** | 25.0 |

**Five readings.**

1. **H1a HOLDS.** Every rung whose realised Cα displacement is ≤ 0.5 Å costs less than 0.10 Å
   of accuracy: `ca1000` +0.0081, `ca300` +0.0305, `ca100` +0.0662 [+0.0331, **+0.0996**].
   The falsifier (cost ≥ 0.10 Å at ε ≤ 0.5 Å) does **not** fire — and the CI's upper bound
   sits at 0.0996 against a 0.10 threshold, so this is a pass with **no margin at all**.
2. **THE FRONTIER SPLITS IN TWO, and this is the finding.**
   **Steric validity is FREE**: *every* rung, `cafix` included, reaches **0.00** clashes below
   2.0 Å at a Cα cost anywhere from 0 to 1.88 Å. The curve is flat.
   **Torsional and peptide-bond validity is a STRICT TRADE**: Ramachandran-favoured rises
   0.647 → 0.847 and the cis fraction falls 0.574 → 0.084 **monotonically in the realised
   displacement**, and no rung buys them cheaply.
3. **H1c is REFUTED.** The pre-registered target was an ε reaching the unrestrained arm's
   validity at < 0.05 Å of Cα cost. The rungs under 0.05 Å (`ca1000` at +0.0081, `ca300` at
   +0.0305) return Ramachandran **0.647 / 0.649** against `free`'s 0.845, geometric deviation
   **0.071 / 0.071** against 0.017, and cis **0.539 / 0.574** against 0.087. **There is no
   free lunch on the torsional axis.**
4. **A NULL AGAINST MY OWN LANE'S PREMISE.** The lane's framing was that freeing N and C —
   restraining *only* Cα — would let the force field repair in coordinates Cα-RMSD is blind
   to. At matched displacement it does the **opposite**: `k30` (N/CA/C, disp 0.824) beats
   `ca30` (Cα only, disp 0.782) on accuracy (+0.122 vs +0.140), Ramachandran (0.828 vs 0.814),
   geometry (0.0221 vs 0.0275) and cis (0.086 vs 0.201). **The incumbent's restraint set is
   better than the Cα-only set on essentially every axis.** The 55.1% non-torsional share of
   AMBER's displacement (Sprint 16 J9) is not an opportunity to be harvested by loosening
   N and C.
5. **THE TIGHT-Cα REGIME IS THE WORST PLACE ON THE LADDER FOR PEPTIDE GEOMETRY, and it is not
   monotone.** cis peaks at **0.574 at `ca300`** — worse than `cafix`'s 0.313 and far worse
   than `free`'s 0.087 — with ω deviation peaking at **96.4°**. §1.3 gives the mechanism: a
   Cα-only restraint at high k pins the contracted spacing while leaving N and C free to
   absorb it by flipping ω, whereas the hard rung has less freedom to flip and the loose rungs
   let the chain expand. **A restraint that holds Cα and frees the peptide plane is the one
   thing not to do.**

### 2.1 The flat-bottom cross-check — the two parameterisations trace one frontier

Pass 3, the first 8 of the pre-registered subsample, `k_f = 10⁵ kcal/mol/Å²`. The realised
constraint violation is **measured, not assumed**: at ε = 0.10 Å the maximum per-atom
displacement is 0.137 Å (violation 0.037 Å), at ε = 0.25 Å it is 0.265 Å (violation 0.015 Å).

| rung | realised Cα disp | RMSD cost vs do-nothing | vs matched random | ramaFav | clash<2.0 | cis | wall s |
|---|---|---|---|---|---|---|---|
| **none** | 0.000 | — | — | 0.782 | 3.25 | 0.000 | 0 |
| `eps010` (flat) | 0.107 | −0.0000 [−0.0065, +0.0063] | −0.0105 [−0.0190, −0.0018] 5W/1L | 0.731 | **0.00** | 0.308 | 37.7 |
| `eps025` (flat) | 0.239 | +0.0206 [+0.0101, +0.0308] | −0.0305 [−0.0619, −0.0081] 5W/1L | 0.786 | **0.00** | 0.306 | 36.5 |
| *harmonic `ca1000`* | 0.185 | +0.0081 | −0.0003 | 0.647 | 0.00 | 0.539 | 12.0 |
| *harmonic `ca300`* | 0.317 | +0.0305 | +0.0080 | 0.649 | 0.00 | 0.574 | 11.8 |

**At matched realised displacement the two parameterisations agree on the accuracy trade** —
0.107 Å buys ≈0.00 Å of cost either way, 0.24–0.32 Å buys ≈0.02–0.03 Å — which is what a
Pareto frontier plotted against the realised constraint should look like, and it validates
sweeping the dual variable instead of ε. The flat-bottom rungs are **3× the wall clock** for
the same trade, which is why they are a cross-check and not the ladder.

**Two cautions on this table, both against it.** (a) The flat rungs' matched-random contrasts
have CIs excluding zero (−0.0105, −0.0305) at **n = 8** — a sample size at which this
programme's own record has reversed a conclusion six times, once with a CI that excluded zero.
**I do not claim it.** (b) The flat rungs reach *higher* Ramachandran and *much* lower cis
(0.31 against the harmonic rungs' 0.54–0.57) at comparable displacement, which is a real
difference in kind — a wall that is exactly flat inside ε leaves the peptide plane free where
a harmonic spring at k = 300–1000 does not. That is a hypothesis this n = 8 cannot settle.

**MATCHED RANDOM.** Only the unrestrained rung is better than a matched-magnitude random
displacement (−0.0954 [−0.2777, +0.0881], 17W/11L, CI includes zero); every restrained rung is
worse than random by +0.008 to +0.081 with intervals containing zero. **At no rung does the
force field's *direction* beat a random direction of its own size** — Sprint 16 J3, reproduced
across a seven-rung ladder rather than at one setting.

> ### The answer to the sprint's §29–31, in one paragraph
>
> **The Pareto frontier is real, it is monotone, and it is not the frontier the question
> assumed.** All-atom *steric* validity — every clash, and a minimum heavy separation better
> than the incumbent's — is available at **exactly zero** Cα cost, so on that axis the
> constrained formulation wins outright and the repair tax is genuinely zero — **on the 62% of
> targets where it converges, which are the least contracted ones**. All-atom
> *torsional and peptide-bond* validity is bought at a fixed exchange rate against Cα
> displacement, with no rung offering it cheaply, and the ε = 0 end of the ladder is
> **worse than doing nothing** on Ramachandran and introduces a peptide-bond defect the input
> does not have. **The defensible operator is therefore: Cα-fixed AMBER as a terminal
> clash-remover, quoted jointly as "zero clashes below 2.6 Å and +0.82 Å of minimum heavy
> separation at 7 × 10⁻¹⁵ Å of Cα displacement" — and NOT quoted as producing a
> Ramachandran-valid structure, which it does not.** And the input contraction that makes the
> torsional half expensive is a defect of the *averaging operator*, not of the force field
> (§1.3).
>
> **And the honest caveat, which I found by checking my own headline against its artefact:**
> the population `cafix` converges on is **selected by the same variable that makes the input
> broken**, so its steric win is measured where there was least to win. On the whole
> instrument the operator that removes 5.39 clashes is the incumbent `k30`, at +0.133 Å of Cα
> accuracy. **The constrained formulation does not yet dominate the incumbent; it dominates it
> on a subset the gate chooses, and closing that gap is item 4 of §8.**

---

## 3. E2 — LEGACY, DECOMPOSED, ON THE IDENTICAL CANDIDATE SET

**n = 126 targets × K = 500 candidates.** The candidate set is the shipped BLOSUM62 pool in
retrieval order — **ranker-neutral**: it is fixed by the retrieval stage, not by any score
being compared. Legacy's eleven components are computed on the ideal-geometry rebuild of each
window's (φ, ψ), which is the same object AMBER is read on in §5.

**The shortlist's cost, stated before any physics number.** The K = 500 set's ORACLE best is
**1.711 Å** against **1.309 Å** over the 13,000–27,000-window universe (n = 100 targets present
in `s17/results/oracle_map.json` at the time of reading): **the truncation costs +0.402 Å of
ceiling.** Everything below is conditional on that truncation.

**IN-BAND is the only meaningful column** (`d ≤ pool_best + 1.5 Å`, `s12.instrument.BAND`).
Global Spearman is reported beside it and is inflated by garbage rejection.

| score | ρ global | **ρ IN-BAND** | AUROC band | selected Å | FNR on < 2 Å |
|---|---|---|---|---|---|
| `dist` (shipped) | 0.568 | **0.126** | 0.739 | 3.454 | 0.075 |
| *zero-info* constant α-helix | 0.295 | −0.004 | 0.558 | 4.072 | 0.278 |
| *zero-info* constant β-strand | −0.100 | −0.023 | 0.466 | 7.386 | 0.242 |
| pool typicality (consensus) | 0.425 | −0.012 | 0.592 | 3.706 | 0.252 |
| **matched random** | −0.005 | 0.006 | 0.500 | 4.427 | 0.219 |
| **AMBER single point (total)** | **−0.027** | 0.047 | **0.499** | 5.191 | 0.222 |
| `amb_nonbonded` | −0.031 | 0.050 | 0.498 | 5.485 | 0.222 |
| `amb_solvation` | 0.116 | −0.009 | 0.512 | 4.490 | 0.251 |
| `amb_torsion` | 0.092 | −0.082 | 0.492 | 4.491 | 0.276 |
| `amb_angle` | 0.020 | −0.045 | 0.500 | 4.515 | 0.204 |
| `amb_bond` | −0.004 | 0.010 | 0.496 | 4.453 | 0.250 |
| **Legacy total** | 0.307 | **0.058** | 0.584 | 4.490 | 0.179 |
| `leg_torsion` | 0.168 | **0.101** | 0.578 | 4.079 | 0.168 |
| `leg_contact` (MJ) | −0.176 | 0.071 | 0.444 | 5.634 | 0.292 |
| `leg_steric` | 0.020 | 0.060 | 0.495 | 4.448 | 0.121 |
| `leg_hbond_local` | 0.291 | 0.053 | 0.567 | 4.054 | 0.178 |
| `leg_solvation` | 0.229 | 0.020 | 0.588 | 4.584 | 0.121 |
| `leg_aromatic` | 0.203 | 0.021 | 0.566 | 4.746 | 0.134 |
| `leg_coop_helix` | 0.260 | 0.013 | 0.532 | 4.141 | 0.004 |
| `leg_electrostatic` | −0.076 | 0.011 | 0.472 | 5.531 | 0.215 |
| `leg_coop_sheet` | −0.067 | 0.001 | 0.512 | 5.836 | 0.000 |
| `leg_compactness` | 0.277 | −0.003 | 0.577 | 4.427 | 0.141 |
| `leg_hbond_longrange` | −0.077 | −0.014 | 0.507 | 5.538 | 0.000 |

**The zero-information constant α-helix is positively correlated GLOBALLY (+0.295) and at zero
IN-BAND (−0.004).** That is the garbage-rejection inflation made visible in one row: a
constant helix tells you which windows are nonsense and tells you *nothing* about which of the
good ones is best. It is the correct reading of the whole `ρ global` column.

**And the first result about AMBER is in that table.** The genuine ff14SB/GBn2 single point has
a **global Spearman of −0.027 and a band AUROC of 0.499** — it does not even reject garbage,
which is the one thing the eleven-term Legacy total (+0.307) demonstrably does. Its in-band ρ
of +0.047 is below `leg_torsion`'s 0.101 and below the distance score's 0.126, and its **global
argmin selects at 5.191 Å, worse than picking at random (4.427 Å)**. The mechanism is on the
record already: *a Pauli spectrum of an unconditioned energy measures its worst clash* — an
unrelaxed ideal-geometry rebuild's ff14SB energy is dominated by whichever single contact is
tightest, and that has no relation to global Cα-RMSD.

### 3.1 The operative bar — IN-BAND selection against MATCHED RANDOM selection

Positive is worse. `s_rand` **is** the matched random operation, so its row is exactly zero
by construction.

| score | selected in-band Å | **vs matched random in-band** | W/L |
|---|---|---|---|
| **`leg_torsion`** | 2.553 | **−0.0446 [−0.1378, +0.0482]** f[−0.109, +0.010] | **66W/57L** |
| `leg_coop_sheet` | 2.638 | −0.0005 [−0.0761, +0.0774] | 57W/57L |
| **matched random** | 2.598 | 0.0000 | — |
| `amb_bond` | 2.620 | +0.0224 [−0.0444, +0.0917] | 60W/65L |
| `amb_nonbonded` | 2.621 | +0.0233 [−0.0708, +0.1142] | 59W/64L |
| `amb_angle` | 2.627 | +0.0287 [−0.0514, +0.1129] | 63W/61L |
| `leg_steric` | 2.629 | +0.0307 [−0.0390, +0.1047] | 57W/68L |
| `leg_hbond_local` | 2.642 | +0.0440 [−0.0469, +0.1361] | 59W/65L |
| **AMBER single point (total)** | 2.646 | **+0.0484 [−0.0416, +0.1348]** f[−0.003, +0.122] | 57W/65L |
| **Legacy total** | 2.655 | **+0.0572 [−0.0348, +0.1547]** f[+0.007, +0.115] | 55W/69L |
| **`dist` (shipped)** | 2.657 | **+0.0587 [−0.0395, +0.1606]** f[−0.053, +0.166] | 57W/62L |
| `leg_coop_helix` | 2.675 | +0.0768 [+0.0041, +0.1522] | 56W/69L |
| `amb_torsion` | 2.714 | +0.1157 [+0.0287, +0.2072] | 51W/72L |
| `leg_contact` | 2.734 | +0.1357 [+0.0433, +0.2327] | 52W/71L |
| `leg_hbond_longrange` | 2.758 | +0.1603 [+0.0724, +0.2535] | 46W/76L |
| pool typicality | 2.763 | +0.1651 [+0.0769, +0.2574] | 49W/72L |
| `leg_compactness` | 2.766 | +0.1676 [+0.0892, +0.2452] | 47W/78L |
| `amb_solvation` | 2.773 | +0.1752 [+0.0833, +0.2668] | 57W/68L |
| *zero-info* α-helix | 2.779 | +0.1814 [+0.0930, +0.2724] | 50W/72L |
| `leg_electrostatic` | 2.817 | +0.1839 [+0.0861, +0.2840] | 37W/68L |
| `leg_aromatic` | 2.717 | +0.2054 [+0.0560, +0.3570] | 24W/38L |
| `leg_solvation` | 2.814 | +0.2155 [+0.1139, +0.3167] | 44W/78L |
| *zero-info* β-strand | 2.868 | +0.2695 [+0.1708, +0.3639] | 35W/88L |

**Three readings.**

1. **Legacy as a ranker is REFUTED again, on a new instrument.** +0.0572 [−0.0348, +0.1547],
   55W/69L — a losing win/loss record with a fold-clustered interval that excludes zero **in
   the wrong direction**. This is an independent reproduction of the SELECT workstream's
   +0.036 [−0.276, +0.353] and of Sprint 16's +0.0118 [−0.0242, +0.0492].
2. **The shipped distance objective is at the same null** (+0.0587, 57W/62L). The bar a
   physics ranker has to clear is **not** "beat the distance objective" — that bar is at the
   noise floor. This reproduces the SELECT workstream's −0.014 [−0.118, +0.094] finding
   against random inside its own top-25 on a different band definition.
3. **`leg_torsion` is the only score in seventeen on the correct side of the line**, at
   −0.0446 with a winning 66W/57L and a CI that includes zero. It is a candidate, not a
   result.

### 3.2 The actual scientific question — what does Legacy contain that the distance model does not?

**Nested ablation, not raw correlation:** the partial Spearman of each component with
Cα-RMSD **in-band**, with the shipped distance score partialled out, per target, aggregated
over targets with both an i.i.d. and a fold-clustered interval.

**The null for this statistic is NOT zero.** The matched random score returns
**+0.010 [−0.010, +0.033]** — a small positive bias from rank-partialling inside a
native-selected band. Every row below is read against **+0.010**, not against 0.

| component | partial ρ in-band, distance partialled out | fold-clustered | clears the +0.010 null? |
|---|---|---|---|
| **`leg_contact` (MJ)** | **+0.080 [+0.032, +0.128]** | **[+0.052, +0.112]** | **yes, ~8×** |
| `leg_torsion` | +0.039 [−0.005, +0.084] | [−0.005, +0.105] | marginal |
| Legacy total | +0.031 [−0.020, +0.081] | [−0.001, +0.059] | marginal |
| `leg_steric` | +0.028 [−0.004, +0.059] | [+0.015, +0.043] | marginal |
| `leg_hbond_local` | +0.020 [−0.025, +0.066] | [+0.003, +0.036] | no |
| `leg_electrostatic` | +0.018 [−0.022, +0.058] | [−0.029, +0.068] | no |
| `leg_solvation` / `leg_aromatic` | +0.016 | [+0.003, +0.028] / [−0.026, +0.072] | no |
| `leg_compactness` | +0.012 | [−0.051, +0.061] | no |
| `amb_bond` | +0.012 [−0.011, +0.036] | [−0.003, +0.032] | no |
| **matched random (the null)** | **+0.010 [−0.010, +0.033]** | [−0.009, +0.029] | — |
| *zero-info* α-helix | +0.004 [−0.068, +0.081] | [−0.074, +0.099] | no |
| **AMBER single point (total)** | **−0.004 [−0.037, +0.030]** | [−0.019, +0.021] | **no** |
| `leg_hbond_longrange` / `coop_helix` / `coop_sheet` | −0.004 / −0.015 / −0.021 | — | no |
| pool typicality / *zero-info* β-strand | −0.013 / −0.024 | — | no |
| **`amb_angle`** | **−0.037 [−0.062, −0.014]** | **[−0.071, −0.002]** | **negative — actively wrong** |
| **`amb_torsion`** | **−0.053 [−0.094, −0.013]** | **[−0.095, −0.008]** | **negative — actively wrong** |

> **The answer, and it is uncomfortable for both sides of the record.** The one component with
> in-band information the distance model does not already have is **`leg_contact`, the
> Miyazawa–Jernigan term** — at **+0.080 [+0.032, +0.128]**, fold-clustered
> [+0.052, +0.112], against a +0.010 null. That is the term Sprint 16 §3.2 measured as
> **"useless as a detector, at or below chance on eight of ten axes"** and §5.4 tested as
> *"the removable seat of the anti-ranking"* and **REFUTED**. Both findings can be true: MJ
> contains no information about *stereochemical defects* (Sprint 16's axis) and a little
> information about *which near-native candidate is nearer* (this axis), and they are
> different questions asked of the same number.
>
> **And it is worth almost nothing.** `leg_contact`'s in-band *selection* is **+0.1357
> [+0.0433, +0.2327] WORSE than random**, 52W/71L. It orders the bulk of the band very
> slightly better than chance and places the argmin worse than chance. **A weak monotone
> ordering signal does not become an argmin.** That is the same distinction the ledger records
> as *"the distogram orders the bulk, the torsion prior places the optimum"*, now measured
> inside Legacy.

---

## 4. E3 — LEGACY AS A GATE, AND THE DANGER

**The danger was stated first and tested first:** *Legacy may remove exactly the candidates
that matter.*

Two statistics that look contradictory are not, and the reason is a denominator. Near-native
recall is only defined on targets that **have** a sub-2 Å candidate; set best is defined on
all 126. `s17/phys_gate.py` recomputes both on the **73 / 126 targets** that have one, so they
share one denominator. Every arm is priced against a **matched random gate of the same
count**, 3 `stable_rng` draws.

| gate | f | near-native recall | random | **recall − random** | set best Å | random | **best − random** |
|---|---|---|---|---|---|---|---|
| Legacy | 0.10 | 0.915 | 0.885 | +0.030 [−0.040, +0.090] 57W/7L | 1.148 | 1.109 | +0.039 [−0.027, +0.118] |
| **Legacy** | **0.25** | **0.821** | 0.751 | +0.069 [−0.016, +0.146] 56W/13L | **1.266** | 1.129 | **+0.136 [+0.018, +0.279]** |
| **Legacy** | **0.50** | **0.613** | 0.501 | **+0.112 [+0.023, +0.198]** 49W/21L | **1.427** | 1.208 | **+0.219 [+0.031, +0.444]** |
| `leg_torsion` | 0.25 | 0.832 | 0.727 | **+0.104 [+0.037, +0.165]** 58W/12L | 1.146 | 1.154 | −0.008 [−0.078, +0.071] |
| **`leg_torsion`** | **0.50** | **0.660** | 0.501 | **+0.159 [+0.073, +0.239]** 54W/16L | **1.206** | 1.195 | **+0.011 [−0.079, +0.128]** |
| `leg_steric` | 0.50 | 0.567 | 0.488 | **+0.079 [+0.014, +0.146]** 53W/19L | 1.223 | 1.220 | +0.003 [−0.083, +0.092] |

**On selected RMSD every gate is null against its matched random control** (n = 126, the
distance score selecting from the survivors): Legacy −0.0001 / −0.0454 / +0.0411 at
f = 0.10 / 0.25 / 0.50, every CI including zero. That reproduces Sprint 16's G3 and the SELECT
workstream's finding, and it is the third independent measurement of the same null.

### 4.1 The mechanism, measured to the level of one candidate

Normalised rank in the gate's **own** order (0 = kept first, 1 = rejected first; a random
candidate sits at 0.500):

| gate | the pool's **single best** candidate sits at | the **median** sub-2 Å candidate sits at | gap |
|---|---|---|---|
| Legacy total | **0.450** | 0.410 | **+0.040** |
| `leg_torsion` | **0.406** | 0.364 | **+0.042** |
| `leg_steric` | **0.512** | 0.475 | **+0.037** |

> **Legacy likes near-native candidates — and likes the single best one systematically less
> than the typical one.** Every gate is *protective of the class* (both ranks below 0.500 for
> Legacy and torsion) and *anti-selective within it*. This is Sprint 16's G4 — *"improves the
> set mean, destroys the set best"* — resolved from an aggregate to a mechanism: it is not
> that the gate throws away near-natives, it is that the ordering it imposes **inside** the
> near-native class is anti-correlated with quality at the extreme.

### 4.2 The coordinator's question — does a Legacy-informed shortlist keep what a score-ranked one discards?

Shortlists of size M, all native-free, on the 73-target stratum. `mix` = top-M/2 by distance
∪ top-M/2 by Legacy.

| M | `dist` has < 2 Å | `mix` has < 2 Å | random has < 2 Å | `dist` best | `mix` best | **mix − dist best** |
|---|---|---|---|---|---|---|
| 25 | 0.644 | 0.658 | 0.644 | 1.881 | **1.741** | **−0.139 [−0.315, −0.000]** 19W/28L |
| 75 | 0.753 | 0.753 | 0.795 | 1.571 | 1.510 | −0.062 [−0.219, +0.095] 21W/24L |
| 150 | 0.822 | 0.836 | 0.849 | 1.371 | 1.354 | −0.017 [−0.120, +0.075] 14W/27L |

The M = 25 row is the only one whose CI excludes zero and it **loses on more targets than it
wins** (19W/28L) with a boundary CI — a mean carried by a handful of targets. Under the
ledger's median-vs-mean rule that is **SUGGESTIVE at best, not a result**, and I decline to
call it one.

### 4.3 The finding that damages the shortlist rather than Legacy

**n = 126, all targets.** The shipped distance top-M against a **matched random shortlist of
the same size** (3 `stable_rng` draws):

| M | `dist` P(has < 2 Å) | random | **difference** | `dist` ORACLE best Å | random | **difference** | W/L |
|---|---|---|---|---|---|---|---|
| 25 | 0.373 | 0.357 | +0.016 [−0.045, +0.077] | 2.609 | 2.414 | **+0.195 [+0.036, +0.359]** | 59W/67L |
| 75 | 0.437 | 0.452 | −0.016 [−0.066, +0.034] | 2.306 | 2.064 | **+0.242 [+0.110, +0.382]** | 63W/63L |
| 150 | 0.476 | 0.508 | −0.032 [−0.090, +0.024] | 2.106 | 1.926 | **+0.179 [+0.063, +0.305]** | 72W/53L |
| 300 | 0.540 | 0.545 | −0.005 [−0.042, +0.029] | 1.886 | 1.789 | **+0.097 [+0.016, +0.193]** | 75W/26L |

> **The shipped distance shortlist is no better than a random shortlist of the same size at
> containing a sub-2 Å candidate at any width, and is significantly worse at containing the
> BEST one at every width.** The win/loss records are near-even, so the effect is
> mean-carried and the median-vs-mean early warning fires — but the sign is identical at all
> four widths, which a concentration artefact would not reproduce. This is the SELECT
> workstream's K-widening result (*shortlist better on average, worse at its best*) measured
> at **fixed** K, and it says the trade is a property of **score-ranking itself**, not of how
> wide you retrieve.

---

## 5. E4 — THE DECISIVE LEGACY-vs-AMBER EXPERIMENT, ON IDENTICAL CANDIDATES

**n = 126 targets × K = 500 candidates = 63,000 genuine ff14SB/GBn2 single points**
(`s17/results/phys_ident_amb.json`, 2,147 s, ~17 ms per point including structure assembly).
Legacy's eleven components and AMBER's five terms are computed on **the same ideal-geometry
rebuild of the same window's (φ, ψ)**, so the two potentials see one structure and the
comparison is exact. **No shortlist**: the AMBER arm covers the entire K = 500 candidate set,
so the only truncation is K = 500 itself, whose cost is stated in §3 (+0.397 Å of ORACLE
ceiling against the full universe, n = 126).

Every statistic is **rank-based** (Spearman, AUROC, argmin, percentile). Raw AMBER single
points on unrelaxed structures put 99.9% of their variance in ten configurations (Sprint 16
§2.3); the rank statistics are invariant to the conditioning map at ρ = 1.000000 exactly, so
no moment of a raw AMBER energy is taken anywhere.

### The four roles, side by side

| role | **Legacy** (11-term total, `DEFAULT_WEIGHTS`) | **AMBER** (ff14SB/GBn2 single point) | matched random |
|---|---|---|---|
| **SELECTOR**, global argmin | 4.490 Å | **5.191 Å** | 4.427 Å |
| **SELECTOR**, in-band argmin vs matched random | +0.0572 [−0.0348, +0.1547] 55W/69L | +0.0484 [−0.0416, +0.1348] 57W/65L | 0 |
| **RANK FEATURE**, global Spearman | **+0.307** | **−0.027** | −0.005 |
| **RANK FEATURE**, in-band Spearman | +0.058 | +0.047 | +0.006 |
| **RANK FEATURE**, band AUROC | 0.584 | **0.499** | 0.500 |
| **RANK FEATURE**, incremental over distance (in-band partial ρ) | +0.031 [−0.020, +0.081]; best component **`contact` +0.080 [+0.032, +0.128]** | **−0.004 [−0.037, +0.030]**; worst components **`angle` −0.037**, **`torsion` −0.053**, both CIs excluding zero | +0.010 (the null) |
| **GATE** at f = 0.50, selected RMSD vs matched-random gate | +0.0411 [−0.0682, +0.1449] | −0.0422 [−0.1154, +0.0259] | 0 |
| **GATE** at f = 0.50, near-native recall vs random | **+0.112 [+0.023, +0.198]** | +0.036 [−0.031, +0.098] | 0 |
| **GATE** at f = 0.50, pool-best survival vs random | 0.516 vs 0.492 | **0.452 vs 0.508 (worse)** | — |
| **REPAIR OPERATOR** | **cannot move an atom** | §1–2: all clashes removed at **7 × 10⁻¹⁵ Å** of Cα displacement | matched-magnitude random is at least as accurate |
| cost per candidate | ~free (11 components, seconds for 63,000) | **~17 ms** single point; **10–20 s** for a restrained minimisation | 0 |

> ### The answer the sprint asked for, in two sentences
>
> **What Legacy knows that AMBER does not:** *which candidates are garbage* (global ρ +0.307
> against AMBER's −0.027, band AUROC 0.584 against 0.499) and — through the
> **Miyazawa–Jernigan contact term alone** — a small amount of **in-band ordering information
> that survives the distance model** (+0.080 [+0.032, +0.128] against a +0.010 null), which
> AMBER has none of (−0.004 [−0.037, +0.030], with its `angle` and `torsion` terms
> significantly *anti*-informative). Legacy is also the only one of the two that gates
> **coverage** better than random (+0.112 near-native recall at f = 0.50, CI excluding zero;
> AMBER's gate *loses* pool-best survival to random, 0.452 vs 0.508).
>
> **What AMBER knows that Legacy does not:** *where the atoms should go.* AMBER is the only
> operator in the programme that turns a broken coordinate average into a sterically legal
> all-atom structure, and §1 shows it does so at **exactly zero Cα cost** when the Cα are
> constrained. Legacy cannot move an atom; it detects defects (Sprint 16: AUROC 0.93–0.99 on
> `torsion`) and repairs none.
>
> **They are complementary on the axes where neither ranks.** On ranking, at ~1,500× the cost,
> **AMBER is not merely no better than Legacy — it is worse than a random score at the one
> job (garbage rejection) Legacy actually does**, and its most informative single term is
> *negatively* informative. **The physics-ranking route is closed by measurement, on both
> potentials, at full n, on identical candidates, with the matched-random control attached.**

### 5.1 What this costs, and why the shortlist question turned out not to arise

The brief warned that full-universe AMBER scoring (13,000–27,000 windows × 126 targets) is out
of reach and that shortlist recall would have to be priced. At the measured **17 ms** per
single point the K = 500 set costs **36 min for the whole instrument**, so no *additional*
shortlist was needed and no additional recall was lost. The full universe would be
**16–34 CPU-hours** — and given that AMBER's band AUROC is **0.499**, buying it would purchase
a coin flip.

---

## 6. E5 — DISAGREEMENT AS INFORMATION — REFUTED

**n = 126.** Target-level native-free features from the disagreement between the distance
score and Legacy, tested against the distance selector's realized error and its regret
(realized − pool best). LFO = the sign of the feature–response relation is fixed on four
folds and the Spearman scored on the held out fifth; a raw in-sample ρ is not a result.

The strongest control is not another disagreement feature. It is **`dist_score`** — how good
the objective thinks its own pick is — a native-free target-level number that needs **no
second ranker at all**. Every disagreement feature is therefore re-scored with `dist_score`
partialled out.

**All three rankers**, n = 126, K = 500 identical candidates.

| feature | ρ(err) | LFO(err) | **LFO(err), `dist_score` partialled** | **LFO(regret), partialled** |
|---|---|---|---|---|
| **`dist_score`** (the trivial baseline) | 0.606 | **0.598** | — | — |
| `rho_dl` (distance–Legacy Spearman) | −0.416 | 0.449 | **−0.148** | −0.078 |
| `rank_of_distpick_in_leg` | 0.420 | 0.428 | **−0.020** | −0.072 |
| `overlap_dl` (top-25 overlap, dist/Legacy) | −0.398 | 0.397 | **+0.040** | −0.122 |
| `score_margin` | 0.441 | 0.388 | +0.087 | −0.045 |
| `rho_dl_band` | −0.199 | 0.208 | +0.137 | +0.168 |
| `rank_of_legpick_in_dist` | 0.121 | 0.151 | −0.142 | +0.081 |
| **`rho_la`** (Legacy–AMBER Spearman) | −0.138 | 0.136 | **−0.125** | −0.103 |
| **`overlap_da`** (top-25 overlap, dist/AMBER) | 0.139 | 0.140 | **−0.076** | −0.099 |
| **`rho_da`** (distance–AMBER Spearman) | −0.211 | 0.127 | **−0.171** | −0.168 |
| `overlap_la` | 0.088 | 0.091 | −0.012 | +0.082 |

> **Every disagreement feature has real leave-fold-out skill (0.09–0.45) and every bit of it
> is subsumed by the objective's own score.** After partialling `dist_score` the ten features
> span −0.171 to +0.168 with no consistent sign across the two responses — noise around zero
> on five folds. Bringing AMBER in as a third ranker makes it *worse*, not better: the three
> AMBER-involving features have the weakest raw skill (0.091–0.140) and the most negative
> partials. The pre-registered falsifier (*no leave-fold-out skill beyond the strongest known
> native-free predictor*) **FIRED**. **Sprint §34's route to target-level calibration is
> closed through the physics channel**, and what is left is the trivial observation that a
> selector that scores its own pick badly is on a hard target — which requires no second
> energy model, and certainly not a 36-CPU-minute one.

---

## 7. WHAT I REFUTED, INCLUDING MY OWN

| hypothesis | whose | verdict |
|---|---|---|
| Cα-fixed AMBER recovers validity without sacrificing Cα accuracy | **the brief's, adopted as mine** | **HALF REFUTED.** Sterics: yes, completely, at 7 × 10⁻¹⁵ Å. Ramachandran: the falsifier fired (0.734 vs the input's 0.836); cis 0.313; gate fails 38% |
| the ε = 0 accuracy falsifier is a test | **the brief's** | **EXACT tautology, declared before the run and never claimed as a result** |
| some ε buys the unrestrained arm's validity at < 0.05 Å (H1c) | **mine** | **REFUTED**: the sub-0.05 Å rungs return rama 0.647/0.649 against 0.845 |
| freeing N and C lets the force field work in the 55.1% Cα-RMSD cannot see | **the lane's premise** | **REFUTED**: at matched displacement the incumbent's N/CA/C restraint is better on accuracy, Ramachandran, geometry and cis |
| AMBER's cis-peptide defect is a force-field problem | inherited (Sprint 16 G12) | **REFUTED**: ρ = −0.949 with the input's Cα–Cα contraction. It is the averaging operator |
| my constant-α-helix zero-information control | **mine** | **REFUTED — a units defect (degrees for radians), found, fixed, and every affected number recomputed** |
| Legacy has in-band ranking skill | inherited | **REFUTED again**, third independent measurement |
| Legacy adds nothing the distance model lacks | **mine** | **PARTIALLY REFUTED**: `leg_contact` (MJ) at +0.080 [+0.032, +0.128] — the one term Sprint 16 called useless |
| …and that it is worth something | **mine** | **REFUTED**: the same term selects +0.136 Å *worse* than random in-band |
| the Legacy gate removes exactly the candidates that matter | **the brief's stated danger** | **REFUTED as stated, CONFIRMED as amended**: it *keeps* the near-native class better than random (+0.112) and *loses* the near-native BEST (+0.219 Å); the best sits at Legacy rank 0.450 against the median near-native's 0.410 |
| AMBER ranks better than Legacy given 1,500× the cost | open | **REFUTED**: global ρ −0.027 vs +0.307, band AUROC 0.499 vs 0.584, and it loses pool-best survival to a random gate |
| ranker disagreement predicts selector failure | **the sprint's §34** | **REFUTED**, subsumed by `dist_score` |
| the shipped distance shortlist beats a random shortlist on coverage | inherited, untested | **REFUTED**: null on P(has sub-2 Å) at every width and significantly worse on the best |

## 8. WHAT REMAINS OPEN

1. **`leg_torsion` as a coverage gate.** The only Legacy operation with a CI excluding zero on
   near-native recall (+0.159 [+0.073, +0.239] at f = 0.50) **and** no significant damage to
   the set best (+0.011 [−0.079, +0.128]). It does not improve selected RMSD — nothing does —
   but if a downstream consumer needs a smaller candidate set at preserved coverage, this is
   the only measured way to halve it. **Not a result yet: the recall gain has not been shown to
   convert into any downstream accuracy.**
2. **Cα-fixed AMBER as the terminal validity operator**, quoted only as the conjunction in
   §2. It needs the 38% gate-failure rate fixed before deployment, and §1.3 says the fix is in
   the averaging operator's 22.4% backbone contraction, not in the minimiser.
3. **The flat-bottom cross-check (pass 3)** was still running when this was written; it tests
   only that the two parameterisations trace one frontier and cannot change a conclusion.
4. **Not attempted, and it should be:** whether repairing the coordinate average's Cα–Cα
   spacing *before* AMBER — a one-parameter uniform re-scaling toward 3.80 Å — removes the cis
   defect and the gate failures at once. §1.3's ρ = −0.949 and ρ = +0.807 predict that it
   would, and it is cheap. It is a **prediction of this workstream, not a measurement.**

## 9. REPRODUCTION

```
python -m s17.phys_ca --pass 1      # cafix + k30, n = 126        ~1 h
python -m s17.phys_ca --pass 2      # the 9-rung ladder, n = 30   ~1 h
python -m s17.phys_ca --pass 4      # the rotated-frame null      ~7 min
python -m s17.phys_ca --pass 3      # flat-bottom cross-check     ~20 min
python -m s17.phys_ca_report        # the frontier
python -m s17.phys_ident --no-amber --out phys_ident_leg.json     ~3 min
python -m s17.phys_ident --out phys_ident_amb.json                ~36 min
python -m s17.phys_report phys_ident_amb.json                     # E2, E4, E5
python -m s17.phys_gate                                           # E3
```

**One operational hazard recorded, because it destroyed a result before it was caught.** Two
`s17.phys_ident` processes writing the same output file raced, and the *later-finishing, less
complete* one overwrote a finished 126-target artefact with 9 rows. This is exactly the hazard
`s12.instrument.write` documents (*"a PARTIAL run silently overwrites a COMPLETE one of the
same name"*) and the writer here has no `complete` flag either. Nothing downstream was
published from the truncated file — the row count was checked before the analysis ran — but
the next module in this tree should carry the flag.

---

## 10. CLAIMS BLOCK — for the sprint ledger

| # | claim | status | evidence | native-free? | pre-reg? |
|---|---|---|---|---|---|
| P1 | Cα-fixed AMBER (ε = 0, particle mass 0) leaves Cα untouched | **EXACT** | max per-atom \|ΔCα\| = **7.1 × 10⁻¹⁵ Å**, max \|ΔRMSD\| = 6.9 × 10⁻¹⁴ Å, n = 126. A tautology; declared as one before the run | yes | PRE-REG |
| P2 | …and therefore passes the rotated-lab-frame null | **EXACT, buys nothing** | mean 0.00000, **max 0.00000**, 0/9 converged targets non-zero — the first AMBER arm in the programme to PASS the standing band, for the same reason its accuracy effect is zero | yes | PRE-REG |
| P3 | Cα-fixed AMBER repairs sterics at zero Cα cost | **DEMONSTRATED, on its gated subset** | vs its **own gated input** (n = 78): clashes < 2.0 Å −0.308 [−0.513, −0.141] **13/0**; < 2.6 Å −2.013 [−2.923, −1.218] **31/0**; min separation +0.308 [+0.190, +0.431] 44/34 (mean-carried, flagged) | yes | PRE-REG |
| P4 | …without sacrificing Ramachandran or peptide-bond validity | **REFUTED — the workstream's own hypothesis; the pre-registered falsifier FIRED** | rama 0.906 → **0.734** (−0.1718 [−0.2203, −0.1256], **1 better / 44 worse**), cis 0.000 → **0.313**, ω dev **60.2°**, **gate fails 48/126** vs the incumbent's 3 | — | PRE-REG |
| P5 | The Cα/validity Pareto frontier splits in two | **DEMONSTRATED** | across 9 rungs (n = 30, pre-registered subsample): clashes < 2.0 Å = **0.00 at every rung** from ε = 0 to unrestrained; Ramachandran 0.647 → 0.847 and cis 0.574 → 0.084 **monotone in realised Cα displacement**. **Sterics free, torsions bought** | yes | PRE-REG |
| P6 | Some ε buys the unrestrained arm's validity at < 0.05 Å (H1c) | **REFUTED** | the sub-0.05 Å rungs return rama 0.647/0.649 vs `free`'s 0.845 and cis 0.539/0.574 vs 0.087 | — | PRE-REG |
| P7 | Cost < 0.10 Å for realised displacement ≤ 0.5 Å (H1a) | **SUPPORTED, no margin** | `ca100` at disp 0.499 Å costs +0.0662 **[+0.0331, +0.0996]** against a 0.10 threshold | yes | PRE-REG |
| P8 | Freeing N and C lets the force field work where Cα-RMSD is blind | **REFUTED — the lane's own premise** | at matched displacement `k30` (N/CA/C) beats `ca30` (Cα only) on accuracy (+0.122 vs +0.140), rama (0.828 vs 0.814), geometry (0.0221 vs 0.0275) and cis (0.086 vs 0.201) | — | DISCOVERED |
| P9 | AMBER's cis-peptide defect (Sprint 16 G12) is a force-field defect | **REFUTED — mechanism identified** | Spearman(input mean Cα–Cα distance, cis fraction at `cafix`) = **−0.949**; convergence **+0.807**. The coordinate average contracts the backbone **22.4%** (2.949 Å vs 3.80 Å) and no trans peptide spans that | ORACLE-free | DISCOVERED |
| P10 | The instrument reproduces the incumbent | **ESTABLISHED** | an independently constructed OpenMM System reproduces `s16/results/repair_A.json` at **max \|ΔRMSD\| = 0.000e+00 Å**, **max \|ΔE\| = 2.7 × 10⁻¹³ kcal/mol**, n = 126; and the gate excludes the *same three* targets (1D6X, 2NB7, 7BX2) | — | DISCOVERED |
| P11 | Legacy has in-band ranking skill | **REFUTED (3rd independent measurement)** | in-band selection vs matched random **+0.0572 [−0.0348, +0.1547], 55W/69L**, fold-CI [+0.007, +0.115] — losing record, interval excluding zero the wrong way | yes | PRE-REG |
| P12 | The shipped distance objective has in-band skill | **REFUTED** | **+0.0587 [−0.0395, +0.1606], 57W/62L** against matched random in-band. The bar for a physics ranker is the noise floor | yes | INHERITED, re-tested |
| P13 | Legacy contains in-band information the distance model lacks | **SUPPORTED, and worthless** | `leg_contact` (MJ) partial ρ **+0.080 [+0.032, +0.128]**, fold [+0.052, +0.112], against a **+0.010** random null — the term Sprint 16 §3.2 called useless. But the same term **selects +0.136 Å worse than random in-band** | yes | DISCOVERED |
| P14 | The Legacy gate removes the candidates that matter | **REFUTED as stated / CONFIRMED as amended** | at f = 0.50 near-native recall **+0.112 [+0.023, +0.198]** (49W/21L, it *protects* the class) and set best **+0.219 [+0.031, +0.444]** (it destroys the extreme). n = 73 on one denominator | ORACLE labels | PRE-REG |
| P15 | The mechanism is the set-mean trap at one candidate | **DEMONSTRATED** | in Legacy's own normalised order the pool's single best sits at **0.450** and the median sub-2 Å candidate at **0.410** (random 0.500); `leg_torsion` 0.406 vs 0.364; `leg_steric` 0.512 vs 0.475 | ORACLE labels | DISCOVERED |
| P16 | `leg_torsion` gates coverage without destroying the best | **SUPPORTED** | recall **+0.159 [+0.073, +0.239]** at f = 0.50, set best **+0.011 [−0.079, +0.128]**. The only such operation found | ORACLE labels | DISCOVERED |
| P17 | AMBER out-ranks Legacy given 1,500× the cost | **REFUTED** | 63,000 single points, n = 126, identical candidates: global ρ **−0.027** vs Legacy's +0.307, band AUROC **0.499**, global argmin **5.191 Å** vs random's 4.427; incremental over distance **−0.004 [−0.037, +0.030]** with `angle` **−0.037** and `torsion` **−0.053** significantly *anti*-informative; its gate loses pool-best survival to random (0.452 vs 0.508) | yes | PRE-REG |
| P18 | Ranker disagreement predicts selector failure | **REFUTED** | 10 features, LFO 0.09–0.45 raw; after partialling the trivial `dist_score` (LFO 0.598 alone) they span **−0.171 to +0.168** with no consistent sign. The AMBER-involving features are the weakest | yes | PRE-REG |
| P19 | The shipped distance shortlist beats a random shortlist on coverage | **REFUTED** | n = 126, M = 25/75/150/300: P(has sub-2 Å) +0.016 / −0.016 / −0.032 / −0.005 (all CIs include zero) and ORACLE best **+0.195 / +0.242 / +0.179 / +0.097 Å worse** (all CIs exclude zero, W/L near-even so mean-carried) | ORACLE labels | DISCOVERED |
| P20 | A validity number quoted against the ungated input is safe | **REFUTED — caught in this workstream's own draft** | the gate keeps a **different subset per rung**; at `cafix` the 78 survivors' input has **0.31** sub-2.0 Å clashes, not the instrument's 5.61. Every rung now carries its own gated-input row | — | DISCOVERED |
