# S31 — RUNNING STATE

Maintained per charter §22: current leading hypothesis, the strongest live objection to it, its
pre-registered falsifier, what each lane is doing, and what is closed and why. Notes are
newest-first below the headline. **With 7 lanes an unwritten research state is how contradictory
conclusions coexist unnoticed.**

Endpoint: **mean built-chain Cα RMSD, `tuning126`, n = 126 — production 3.2105 Å.**
CA point cloud is **3.0483 Å** and the *set mean* is 3.5507 Å. Three different objects.

---

## HEADLINE: THE CAP IS IN THE READOUT, NOT THE HAMILTONIAN

The sprint opened on the charter's premise that the objective is the barrier — justified by S30's
finding that the same circuit reaches **0.2516 Å** under an ORACLE objective and **3.4330 Å** under
the deployed one. Reading the code moved the diagnosis one stage later:

> **R1 (S31-L1).** The readout is `argmin_i (P p)_i` with `P` the pairwise Kabsch RMSD matrix. The
> state enters only through a linear map then an argmin, so the whole quantum stage carries **at
> most k bits** and **cannot emit anything outside the pool** — and that is a property of the
> **readout**, not the Hamiltonian. T1's one-integer result is the diagonal special case.

**Value is not capped the same way.** ORACLE argmin-over-128 is **2.1435 Å (CA cloud, ORACLE / NOT
DEPLOYABLE)** against production's 3.0483 — ~0.90 Å of headroom. So the cap is on *information*,
not on *value*, and the two must not be conflated.

**Strongest live objection to R1:** the realised capacity may be far *below* k bits if the ansatz
cannot reach the cells — which would be worse and more interesting than R1 itself. Lane A owns the
measurement and was asked to attack the theorem rather than confirm it.

---

## LANES

| lane | remit | status |
|---|---|---|
| **A** | quantum/CVaR theory; verify R1; realised capacity; classical reducibility (§11) | running — prereg filed |
| **B** | free energy (§7A) + torsion (§7B); theoretical pre-check before compute | running |
| **C** | readout design given R1; sparse native-free support; index allocation; the 128→512 gate | running — prereg filed |
| **D** | integrity: projection seed (gates every sub-0.01 Å claim), `pipeline.py:821`, governor/launcher, verifier + multiplicity register | running |
| **E** | the incoherence hypothesis | **E1 CLOSED by an exact identity** — see NOTE 1. Running the ORACLE class ceiling |
| **F** | terminal operator: medoid vs average under a native-free dispersion gate; tail mechanism | running |
| **L** | literature, permanent | running |

One slot of eight held for an adversary once there are results to attack (contract rule 24: two
adversaries found 41 defects in the S30 report, six severe).

---

## NOTE 4 (2026-09-21 00:07, lane D): **THERE IS NO SEED. THE PROJECTION IS DETERMINISTIC AND CHAOTIC** — AMPLIFICATION ~1e13

The brief I wrote, S30, and S28-L18 all called this *"an unpinned multi-start projection seed."*
**There is no RNG anywhere on the projection path.** `core.project.fit_multi` loops over four FIXED
starts and takes a strict argmin; reprojecting the same cloud twice is **bit-identical** (max |ΔCA|
exactly 0.0, unchanged by BLAS thread count).

**The real mechanism is worse than a seed, and it explains every historical discrepancy at once:**

```
1A13's four lam=0 objectives   0.5091924865 / 0.5091925170 / 0.5091924488 / 0.5091924188
                               -> spread 1e-7, i.e. NOISE
the same four branches at lam=0.3   0.947 / 0.519 / 0.562 / 2.396
                               -> spread 1e-1
median lam=0 branch margin      ~4e-8
```

**A decision taken at 1e-7 selects between outcomes that differ at 1e-1.** Demonstrated, not
argued: perturb the input cloud by **1e-14 relative** and 2LNG's chain moves **+0.511 Å**,
reproducing the historical 0.517 Å discrepancy exactly. **Amplification ~1e13.**

**Operational rules now in force for every lane:**
1. Reprojection is reproducible **only from bit-identical clouds**. Both sides of any built-chain
   contrast must be projected **in the same job from the same stored clouds**, and the entry must
   say so.
2. **Canonical endpoint stays 3.2105 Å** (`s29/results/s29_O_chain_rows.jsonl :: item=prod`) — the
   run the endpoint was declared from and the only one whose clouds are persisted per target.
   Not retro-fitted.
3. **No built-chain claim below 0.0107 Å.**

**D-A and D-C also fixed.** The withdrawn "+0.113 Å measured role" is gone from
`core/pipeline.py`'s `quantum_stage` docstring, corrected in place with S25-L5 cited and the old
wording quoted so the correction is visible. And `jobrun.py` v3 now **derives** its launch gate
from `governor.py` by import (CPU gate = `governor.CPU_CEILING` = 101, so CPU no longer blocks;
RAM gate = `governor.CEILING − 1` = 93.0) — live-tested with a launch at `cpu_smooth` 100.0%, which
the old gate would have blocked forever. That is the `paired-thresholds-move-together` lesson
implemented so the two **cannot drift**: the governor changed four times in S29 and the launcher
never followed. **Lanes can stop launching detached.**

**My open proposal to lane D, which may be worth Ångströms rather than only reproducibility:** the
λ=0 argmin is *a selection made on noise*. **Carry all four branches forward and take the argmin at
λ=0.3**, where the objective separates them by 1e-1. Deterministic, ~4× a cheap stage, native-free
(the argmin is on the objective, exactly as now), and it collapses the 1e13 amplification by
construction. Whether it improves RMSD is empirical — measure the ORACLE branch choice first to get
the ceiling. Given the tail arithmetic, a fix worth 0.5 Å on a few tail targets is not a rounding
error.

---

## NOTE 3 (2026-09-21 00:07, lane B): **THE FREE-ENERGY STAGE IS CLOSED BY DERIVATION** — AND `S` IS ACHIRAL FOR THE SAME REASON `E` IS

Lane B answered the pre-check I set and then closed it on a stronger ground than the one I asked
for. **No thermodynamics compute was spent.**

- **The source question:** `S(x) = Φ_{β,B}[U_seq](x)`. Its arguments are the candidate, the sequence
  (only through atom typing), and β plus the basin map, which are universal. **There is no fourth
  argument** — `S` is an *operator* on (sequence, pool), exactly like `E`.
- **LEMMA B1, verified numerically:** the shipped potential is **reflection-invariant**,
  `U(Rx) = U(x)`. ff14SB's `PeriodicTorsionForce` phase table — **1340 torsions, max distance of any
  phase to {0, π} = 0.000e+00 exactly**, no CMAP.
- **COROLLARY B1:** `F_β`, `E`, `S` and every temperature derivative are therefore reflection-invariant
  single-structure observables ⇒ **by G1 each is a function of the distance map**, inside the class
  S30 closed by theorem and measured empty on 43 channels. The enthalpic half was already closed by
  measurement; **the entropic half is closed by the same theorem that closed contact topology and Rg.**
- **COROLLARY B1′** extends it to the whole **elastic-network / normal-mode / landscape-curvature**
  family: an ANM/GNM Hessian is built from pairwise distances, so its spectrum, log-det and harmonic
  entropy are distance-map functions. Basin populations, ensemble reweighting and temperature-dependent
  ranking are monotone transforms of single-structure free energies. **Several charter §14 items close
  here.**

**And the methodological moment worth quoting.** Lane B's registered 1e-6 threshold **fired**
(3.673e-06) — and it ran the matched control anyway: a pure **rotation**, equally analytically exact,
gives **5.377e-06, larger**, so reflection/rotation = 0.68×. The chirality-carrying force groups are
exact to **2.0e-16**; the residual lives entirely in two terms that are analytically invariant under
*any* isometry. **A registered threshold set at an absolute value, on a quantity only a matched
control can read, is a mis-set threshold — not a falsification.** Recorded with the original standing.

**Left open, named so it cannot be retrofitted:** a genuinely **multi-structure** free energy. A
barrier `F‡(x_a → x_b)` depends on the *path*, not on `D(x_a)` and `D(x_b)`, so **G1 does not bind
it** — left open on **price**, not theory (10³–10⁵ trajectories per target). I have proposed a cheap
proxy: **connectivity observables on the candidate graph** (geodesic vs direct distance, commute
time, spectral gap on the already-computed pairwise matrix `P`) are functions of the *set* of
distance maps, so they escape G1 on lane B's own argument and cost one matrix operation. Caveat
carried: they are still **source 2**, which S30 measured as *not empty, worse than empty*.

**And a possible confound in MY OWN published sweep.** Lane B decomposed all eleven Legacy
components into even/odd parts under `(φ,ψ) → (−φ,−ψ)`. Odd (chiral) variance share: **steric
13.02, torsion 0.301**, electrostatic 0.230, contact 0.105, solvation 0.100, aromatic 0.089 — and
**five exact zeros** (both hbond terms, both coop terms, compactness), which are the achiral terms
G1 applies to directly. The mechanism is construction-level: `core/geometry.py build_backbone_batch`
places CB as `-0.58273431*cross(b,d)`, **a pseudovector**, so every CB-dependent term inherits
chirality.

> **The two costs with the largest preference contrast in my 32-cost sweep — `LEG_steric` (+0.2515,
> 2.88× MDE) and `LEG_torsion` (+0.2212, 2.14×) — are exactly the two with the largest odd shares.**
> The sweep may be measuring **chirality detection against the `RAND_SIGNED` rungs** rather than
> structural skill. That sweep is mine and published in the S30 report §10.6; **nobody should build
> on its top rows until lane B reports, and if it confirms, the correction is mine to make.**

---

## NOTE 2 (2026-09-21 00:07, lane C): **MY OWN FALSIFIER FIRED IN SHIPPED CODE, AND THE SET-MATCHED LADDER INVERTS A PUBLISHED S30 CLOSURE**

**(a) R1's quantifier is false.** I wrote the falsifier as *"R1 fails if any code path lets the
stage emit a structure that is not a pool member."* Lane C found that path **in `core/pipeline.py`**:
`average_weighted` (`:880-895`), a **shipped** p-weighted convex average called at `:1110`/`:1116`,
scored at `:1179-1182`, registered as comparison arms at `:1611`/`:1613`. **There are three
readouts, not one and not two** — selection (≤ k bits, pool members, **R1 exactly right**), convex
(shipped, continuous, reaches the convex hull), affine (harness only, signed, leaves the hull). I
compounded the error in S31-L2 by merging the last two; the distinction is load-bearing because
**the convex one ships and the affine one does not.** Both entries annotated in place.

**(b) The set-matched ladder inverts S30-L11.** That entry says argmin dominates *"at every measured
bit budget"* — established by comparing **2-of-75 against 1-of-128**. At a **fixed top-128**
(ORACLE / NOT DEPLOYABLE, CA cloud, n = 126):

```
argmin over 128            7.0 bits              2.1458
2-of-128, UNIFORM weights  12.99 support bits    2.0700     -0.076 vs argmin, ZERO weight bits
2-of-128, continuous w                           1.9138     -0.232
5-of-128, continuous w                           1.8037
```

**As classes at a fixed candidate set, combining beats naming; as an Å-per-bit question across sets,
naming wins.** S30 published only the second and supported it with a set-mismatched pair.

**My prediction handed back, because the uniform row is the interesting one:** 2-of-128 uniform
beats the single best of the same 128 while spending **zero weight bits**, so that 0.076 Å is
**error cancellation**, not weighting. Averaging cancels the i.i.d. component and leaves the shared
one, and this project measures the pool's error as **68% common-mode** — so the ladder should behave
as `sqrt(0.68 + 0.32/m)` and **saturate at the common-mode floor**. If it does, the entire
sparse-combination class closes by derivation with a fitted constant. If it *undershoots*, something
is removing common mode and that is the sprint's most important finding.

**(c) Two rungs already closed in the record**, so nobody re-derives them: the affine hull ceiling is
**exactly 0.0000 by dimension counting** (3n ≤ 48 coordinates, 75 generic windows span it) — *"do
not quote any affine-hull ceiling as a bound"*; and the **convex** rung is closed **by ceiling** at
**3.0522** — the best convex reweighting the deployed objective can find, started from production
and converged, **is the uniform average it already emits**, with two independent constructions
agreeing. **The open rung is the norm-/sparsity-bounded middle.**

## NOTE 1 (lane E, E1): **MY OPENING HYPOTHESIS WAS CIRCULAR, AND AN EXACT IDENTITY KILLED IT**

I proposed estimating the pool's common-mode direction native-free as `mu_hat = pool75_mean −
expected`, in order to apply only the component of a corrector orthogonal to it. Lane E's first
measurement closed it, and the reason is algebra I should have done before briefing:

```
mu     = pool75_mean - d_nat
y      = expected    - d_nat          (the prior error the corrector exists to predict)
mu_hat = pool75_mean - expected  =  mu - y        EXACTLY   (verified, max |dev| 1.8e-15)
```

> **`mu_hat` is not a noisy estimate of `mu`. Its error IS `y`.** To use my estimator you would
> already need the answer. The hypothesis was circular.

**The durable output is a law, not the null.** Whether the family can ever work reduces to one
ratio `r = sd(y)/sd(mu)` through

```
corr(mu_hat, mu) = (1 - coh0 * r) / sqrt(1 + r^2 - 2 * coh0 * r)
```

which reproduces the measured per-target correlation to **1.1e-15**. Measured **r = 1.5997
[1.4421, 1.7549], 100% of targets above 1**, and `coh0` reproduces S30's 0.6931 exactly. At r = 1.6
the closed form is **negative**; the measured cosine of +0.05 is Jensen curvature across targets,
not signal. Only 1% of targets reach the registered usable bar.

**So the specification for any future attempt is quantitative:** you need a common-mode estimator
whose error is *small relative to the common mode*, and the distogram's is **1.6× too large**.

**Scope, so a neighbouring result is not thought contradicted:** `r > 1` compares the **distogram's**
error dispersion to the **pool's common-mode** dispersion. "The pool's error is 68% common-mode" is
about pool *members'* errors relative to the native — a different object. Both hold.

**And the sharpened requirement, which is lane E's wording and better than S30's:**

> **An observable is not useful merely for being incoherent with the common mode — `mu_hat`'s own
> error is incoherent and worth nothing. It must carry orthogonal INFORMATION.**

**My registered prior was 2:1 against and it was right in direction and wrong in mechanism.** I
predicted the orthogonal complement would be noise because the residual sits at coh 0.9172; the
actual failure is one level earlier. A prediction that gets the sign right for the wrong reason is
worth less than a correct one, and the report will say so.

**Redirect issued:** projecting the *existing fitted* corrector (R² 0.2355, coh 0.9172) leaves
almost nothing and a near-zero result would be ambiguous between *"the class is empty"* and *"this
corrector had nothing orthogonal to give."* Lane E is instead running the **class ceiling** — the
best achievable endpoint over all corrections constrained orthogonal to the **true** `mu`
(**ORACLE / NOT DEPLOYABLE**). That answers the question the next sprint needs: *is there enough
signal orthogonal to the common mode to be worth finding, even with perfect knowledge?*

---

## OPEN QUESTIONS I AM HOLDING

- If R1 caps the stage at k bits, **is there a readout that escapes it while preserving physical
  validity** — and is the escape the *sparse weighted* class, which by construction does not return
  a pool member? (Lane C.) Note the set-mean algebra forbids paying for mere *concentration*.
- Does anything survive the charter's §11 reducibility burden, or is the whole stage a classical
  prefix/rank operation wearing a quantum state? (Lane A.)
- Is there any target argument for a free-energy or torsion observable that `(sequence, pool)` does
  not already carry? (Lane B — theoretical pre-check gates the compute.)
- The **128→512 widening** still carries S30's undischarged circularity gate. Cheapest open item;
  either revives or closes a direction. (Lane C.)
