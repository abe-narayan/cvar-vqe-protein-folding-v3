# PREREG — S31 lane E: project the common mode out of the corrector

**Written 2026-09-20 23:45, branch `s26`, before the first number of this lane existed.**
Committed before `s31/s31_E_*.py` was run. Nothing below was chosen after seeing an S31 number.

What I *had* read when writing this: `s31/BRIEF.md` §3/§14, `s30/REPORT_S30.md` §0.2–§0.3 and §12,
`s30/s30_P_lr.py`, `s30/s30_P_prior.py`, `s30/s30_P_chain.py`, and the three S30 result JSONs
(`s30_P_lr.json`, `s30_P_prior.json`, `s30_P_chain.json`). Those are **S30** numbers and they are the
state of the art this lane is trying to beat; they are quoted below as the baseline.

---

## 0. The hypothesis in one line

S30 proved *applying* a fitted corrector is harmful and that its residual is **more** coherent with
the pool's common mode than the uncorrected error. It never tested applying **only the component
orthogonal to the common mode**. The common-mode direction is claimed native-free estimable as
`mu_hat = pool75_mean − expected`.

## 1. The algebra, derived before any measurement (this is not a result, it is arithmetic)

With `y = expected − d_nat` (the prior's error, ORACLE) and
`mu = pool75_mean − d_nat` (the pool's common-mode pair error, ORACLE; `s30_P_lr.py:202`), and
`mu_hat = pool75_mean − expected` (native-free — both terms are features, no `d_nat`):

```
mu_hat  =  (pool75_mean − d_nat) − (expected − d_nat)  =  mu − y        EXACTLY.
```

Three consequences I register *now*, because they decide how the results must be read:

1. **`mu_hat`'s error as an estimate of `mu` is exactly `−y`** — the very quantity the corrector is
   built to predict. The scheme is circular in a precise, quantifiable way: perfect knowledge of
   `mu` is equivalent to perfect knowledge of `y`.
2. **The E1 cosine is therefore a pure ratio statistic.** With `r = sd(y)/sd(mu)` and
   `coh0 = corr(y, mu) = 0.6931` (S30, within target), within-target *correlation* obeys
   `corr(mu_hat, mu) = (1 − coh0·r) / sqrt(1 + r² − 2·coh0·r)`. E1 is, in substance, a measurement
   of `r`. I will print the identity's prediction beside the direct measurement as a self-check; if
   they disagree, the code is wrong, not the theory.
3. **Projecting `yhat` orthogonal to `mu_hat` is projecting orthogonal to `mu − y`**, which contains
   `−y`. So the projection can destroy exactly the useful part of the corrector. This is a stronger
   reason for pessimism than the coordinator's, and it is registered here rather than discovered later.

## 2. E1 — bars, chosen before looking

**Statistic:** within-target mean of `cos(mu_hat, mu)` and of `corr(mu_hat, mu)` over the 126 dev
targets, all pairs at `min_sep = 2` (S30's `coh` space) as primary, `|i−j| ≥ 7` as secondary.
Fold-clustered 95% CI on the pinned folds. Matched control: a per-target random direction drawn with
the **same within-target marginal** as `mu_hat` (a random permutation of `mu_hat`'s own entries
within the target, so norm and marginal are matched exactly and only the pair-to-pair pairing is
destroyed), plus an isotropic Gaussian direction.

**Mechanistic bar.** Projection onto the complement of `mu_hat` removes the fraction `c²` of `mu`'s
energy that lies along `mu_hat`. So:

| within-target mean cosine `c` | verdict, registered now |
|---|---|
| `c ≥ 0.707` (`c² ≥ 1/2`) | **STRONG** — the estimate removes most of the common mode; proceed at full weight |
| `0.577 ≤ c < 0.707` (`c² ≥ 1/3`) | **USABLE** — proceed to E2/E3, flag the estimate as partial |
| `0.30 ≤ c < 0.577` | **WEAK** — proceed to E2/E3 only because they are cheap, and report the whole lane as a bounded test of a weak direction |
| `c < 0.30`, **or** not separated from the matched control by ≥ 1.0× its own MDE | **DEAD** — report immediately, run nothing below |

## 3. E2/E3 — falsifier in Å on the built chain

**Basis, stated every time:** built chain, production **3.2105 Å** (S29 `prod` row) with this lane's
own re-projection measured in the same run and reported beside it (lane P's re-projection was
3.2126 Å; the projection seed is unpinned and the range across instruments is 0.0107 Å). CA point
cloud **3.0483 Å**. All applied deltas are **paired against this lane's own PROD arm computed in the
same process**, so the projection spread cancels.

The arm to beat is S30's `N3_plus_pool` fitted corrector: out-of-fold R² 0.2355, applied
**+0.0554 Å (CA cloud, 0.68× MDE, NOT MEASURED)**, `coh` 0.9172 against uncorrected 0.6931.

Registered outcomes for **E2** (project after fitting) and **E3** (constrain during fitting),
each judged on its own MDE = 2.8016 × SE with fold-clustered CI:

* **RESULT (deployable candidate):** applied built-chain delta **negative** and `|delta| ≥ 1.0× MDE`,
  with ≥ 4/5 folds the same sign. Secondary requirement: `coh` falls **below 0.6931**.
* **NOT MEASURED:** `0.7× ≤ |delta| < 1.0× MDE`.
* **FALSIFIED — the route closes:** `|delta| < 0.7× MDE`, **or** delta positive at any size.

**Mandatory matched controls, registered now** (`control-must-match-the-operators-space` is this
project's most repeated error):

1. **Norm-matched shrinkage.** `yhat · (||yhat_perp|| / ||yhat||)` per target. The projection shrinks
   the corrector; if plain shrinkage buys the same thing, E2 measured shrinkage, not orthogonality.
   **A positive E2 that does not beat this control is not a result.**
2. **Random-direction projection.** Project out a per-target random direction with `mu_hat`'s own
   within-target marginal (same permutation control as E1). If projecting out *any* direction helps
   equally, the common mode is not the mechanism.
3. **PROD** (zero correction) and **S30's unprojected `N3_plus_pool`**, both recomputed in this run.

## 4. Native-free trace — to be verified in code, not asserted

S30 published its `coh` gate as native-free and it was ORACLE; the error was caught by an adversary,
not by the lane. I register the obligation: `s31_E_*.py` must contain an explicit assertion that the
deployable arms' inputs never touch `d_nat`, and the report must print the trace
(`mu_hat = pool75_mean − expected`; `pool75_mean` from `rec["sub"]` over the pool `W`; `expected`
from the leave-fold-out distogram; neither reads `u["nat_ca"]`). Any arm that uses `mu`, `y`, or the
i.i.d. construction is **ORACLE / NOT DEPLOYABLE** and must be labelled so in the same sentence as
its number. The `coh` statistic itself stays ORACLE and is reported as a development-time diagnostic
only.

## 5. Honest prior odds

The coordinator registered **2:1 against** E2/E3 producing a deployable gain. **I register 4:1
against for E2 and 3:1 against for E3**, and the extra pessimism is the §1.3 identity: the direction
being projected out is `mu − y`, so the projection removes corrector capacity that lies along `−y`,
i.e. along the signal. E3 is less bad than E2 only because the constrained fit can re-spend that
capacity on the remaining features rather than losing it.

Against that: `P(E1 usable)` I put at about **0.5**, driven entirely by whether `r = sd(y)/sd(mu)`
is small. If the pool's common-mode error is much larger than the distogram's own error, `c` is high
and the direction is usable; if they are comparable, `c` collapses. I genuinely do not know which,
and `coh0 = 0.6931` alone does not determine it.

**What a null buys, and why the lane is still worth running:** if E1 passes and E2/E3 still emit
nothing, that is the coordinator's clean result — the i.i.d. arm's −0.2466 Å comes from having
genuine orthogonal *signal*, not from being incoherent — and it converts S30's §12 requirement from
"find an incoherent observable" into the sharper "find an observable carrying orthogonal
information", which is a strictly harder and more honest specification. If E1 *fails*, the S30
admission test loses its only candidate native-free surrogate and the §12 requirement stands
unweakened, which is also worth knowing in one hour rather than one sprint.

## 6. What this lane will not do

No reopening of benchmark60, no regeneration of folds or clusters, no tuning of any deployable
parameter on native RMSD. The ridge path, the fold split and the feature list are S30's, unchanged.

---

# AMENDMENT 1 — written 2026-09-20 23:52, AFTER seeing E1 and labelled as such

**E1 came back DEAD** (within-target cosine +0.0495 against a 0.30 kill bar; not separated from the
matched permutation control, 0.51× MDE). Under §2 the registered consequence is "report immediately,
run nothing below". I am amending, and the amendment is recorded here with its timestamp and its
motive so nobody has to reconstruct it later.

**What changes.** §2's stop applies to the *deployable* claim, which is dead and stays dead. It does
not answer the coordinator's actual E2/E3 question, which is about the *operation* (projecting out
the common mode), not about this particular estimator of it. I therefore add an **ORACLE ceiling**:
run the identical projection and constrained-fit operators using the **true `mu`**, which is
**ORACLE / NOT DEPLOYABLE** and is labelled so in every sentence that carries a number from it.

**Why this is worth compute and is not a fishing expedition.** It is a *ceiling*, and a ceiling has a
pre-stated meaning in both directions:

* **If the ORACLE-direction projection emits nothing**, the projection idea closes *at its ceiling* —
  no better estimator of the common mode could rescue it — and S30's §12 requirement sharpens from
  "find an incoherent observable" to "find an observable carrying **orthogonal information**",
  because `mu_hat`'s own error is incoherent and is worth nothing.
* **If it emits a gain**, the next sprint has a named, sized target: estimate the common-mode
  direction. We already know `pool75_mean − expected` is not it, and why (`r = sd(y)/sd(mu) = 1.60`).

**Bars for the amendment, registered now, before the arms are run.** Same as §3 and applied to the
built chain against this lane's own PROD arm in the same process: RESULT at `|delta| ≥ 1.0× MDE` in
the better direction with ≥ 4/5 folds same sign; NOT MEASURED at 0.7–1.0×; **the projection route
closes at `|delta| < 0.7× MDE` or any positive delta**. The §3 matched controls are mandatory and
are now matched **per projection arm**: a norm-matched shrinkage control for each projection, and a
random-direction projection control.

**The registered native-free arms are run anyway**, despite the DEAD gate, because they cost one
extra vector each and a *measured* negative is worth more than an inferred one. They are reported as
registered arms whose gate had already failed, not as new hypotheses.

**Additional diagnostic registered now:** the per-target norm fraction `||P_along(yhat)|| / ||yhat||`
for both directions. The coordinator's stated reason for expecting a null is that the fitted
corrector is "almost entirely in the `mu` direction", so its orthogonal complement is nearly pure
noise. That fraction measures it directly and is reported whatever the endpoint does.

**Prior odds for the amendment, registered before running it:** ~3:1 against the ORACLE projection
producing a gain, for the §1.3 reason (the projection removes corrector capacity lying along `−y`)
which applies to the true `mu` just as it does to `mu_hat`. A null here is the more informative
outcome and is the one I expect.
