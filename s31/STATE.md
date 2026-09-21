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
