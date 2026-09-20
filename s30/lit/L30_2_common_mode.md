# S30 LANE L, TOPIC 2 -- THE COMMON MODE: WHAT IS AND IS NOT IDENTIFIABLE FROM A POOL

Brief priority 2: "methods that break common-mode error ... several of these literatures assume
ensemble errors are independent, which is false for this pool."

Status of the inputs, checked rather than assumed:

```
$ ls s23/results/errdecomp.json          -> present
$ git log --all --oneline -- s23/errdecomp.py
```
(both verified before anything below was written; numbers quoted from the project memory note
`pool-error-is-68-percent-common-mode`, whose BODY -- not its index line -- is the source. The
index line says "invisible from inside the pool"; the body says why, and the body is what this
file formalises.)

---

## 0. WHAT IS ALREADY IN THE RECORD, SO THAT I DO NOT CLAIM IT

The project's own memory note already states the invariance **in prose**:

> "a bias shared by every member moves `c` and every `w_k` together and leaves every within-pool
> statistic unchanged."

That sentence is correct and it is not mine. **This file does not discover the invariance.** What
it adds is four things the prose does not contain, and which change what follows from it:

1. the **exact quantifier** -- what the invariance is over, and, critically, what it is *not* over;
2. the identification of the **model class** under which it is a theorem rather than a tautology;
3. the **complete list of escapes** the theorem permits, which is three and not more; and
4. the observation that **the project has already measured all three**, so the accounting closes.

---

## 1. THE TAUTOLOGY TRAP, STATED FIRST SO THE REAL CLAIM IS NOT CONFUSED WITH IT

A wrong version of this argument is available and sounds identical. It goes: "for any candidate
native `t'`, the world `(t', {w_k - t'})` reproduces the data `{w_k}` exactly, so the native is
unidentifiable." **That is true and it is worthless.** It says only that a data set of structures
does not, by itself, name one of them as the answer. Every estimation problem has that property
before a model is imposed.

The content has to come from the model that restricts the errors, because that is what turns the
pool into evidence. So the theorem below is stated **relative to an error model**, and its force is
that the error model the whole aggregation literature uses is the one that fails.

---

## 2. THE THEOREM

**Setup.** One target. Native `t` in R^(3n). Pool members `w_1 ... w_K` in R^(3n), all in the
common frame the average is actually taken in (`s23/errdecomp.py`'s medoid frame). Member errors
`e_k = w_k - t`.

**Model class M(G, mu).** The members are generated as

        w_k = t + mu + d_k ,        d_1 ... d_K  i.i.d. ~ G,  E_G[d] = 0

with `t` in R^(3n) unknown, `mu` in R^(3n) an unknown **common bias**, and `G` an unknown
mean-zero distribution. This is the smallest class that contains both the i.i.d. model the
literature assumes (`mu = 0`) and the pool we actually have.

**Theorem (non-identifiability of the common bias).** Under M, the distribution of the observed
pool `{w_k}` depends on `(t, mu)` **only through the sum `t + mu`**. Hence `(t, mu)` is not
identifiable: for every `c` in R^(3n), the parameter values `(t, mu)` and `(t - c, mu + c)` induce
*exactly the same distribution* of the data, for every K.

*Proof.* `w_k = (t + mu) + d_k`, and `d_k ~ G` does not depend on `t` or `mu`. The likelihood is
`prod_k G(w_k - (t + mu))`. Substituting `(t - c, mu + c)` leaves every factor unchanged. The
invariance group is the full translation group of R^(3n), acting on `(t, mu)` and trivially on the
data. QED.

**Corollary 1 (what IS identifiable).** `t + mu` is identifiable and is consistently estimated by
the pool mean `c`, with error `O(K^-1/2)` under G. The idiosyncratic deviations `d_k = w_k - c` are
identifiable. **`mu` and `t` separately are not, at any K.**

**Corollary 2 (the exact deficiency).** The error of the pool average decomposes as
`ebar = mu + dbar`, and `dbar -> 0` as K grows. So the asymptotic error of the *best possible*
pool-only estimator under M is exactly `|mu|`, and no amount of pool is worth anything against it.
This is the S23 L9 identity read as an identification statement rather than a variance
decomposition:

        mean_k |e_k|^2  =  |ebar|^2  +  mean_k |d_k|^2
                           160.36        63.82
                           NOT IDENTIFIABLE   IDENTIFIABLE

**Corollary 3 (the effective sample size, and it is the number to quote).** Under the i.i.d. model
(`mu = 0`) the pool average's squared error would be `mean_k |e_k|^2 / K`. Observed, it is
`|ebar|^2`. The ratio is the pool's **effective number of independent members**:

        m_eff = mean_k |e_k|^2 / |ebar|^2 = 224.18 / 160.36 = 1.398

or, per-target-averaged from the note's own `f = 0.676`, `m_eff = 1/0.676 = 1.479`.

> **The deployed 75-member pool carries the statistical content of about 1.4 independent
> members.** Every aggregation, consensus, diversity, voting, bagging and ensemble-QA method in
> the literature prices its gain in K. Priced in `m_eff` instead, this pool is *K = 1.4*, and the
> honest conversion factor from any such paper's promised gain to ours is `(1.4/K_theirs)`.

That single number is the family-level rejection the brief asked for, in the form that transfers
to any future proposal without re-reading the paper.

---

## 3. THE QUANTIFIER, STATED EXACTLY, BECAUSE THE DIFFERENCE IS WHERE THE HOPE LIVES

The coordinator asked for this explicitly and it is the most important paragraph in the file.

**What the theorem says:** `mu` is not identifiable **from the pool, under model class M**. It is a
statement about one data matrix and one invariance group.

**What the theorem does NOT say:** it does not say `mu` is invisible to any method. The invariance
is broken by anything that makes the likelihood depend on `t` and `mu` separately. There are
exactly three ways to do that, and they are not a list I chose -- they are the only three places
the proof can fail:

| # | escape | what it supplies | breaks the proof at |
|---|---|---|---|
| E1 | a restriction on `mu` itself -- a prior on the bias's **form**, not its magnitude | outside knowledge of the shape of the error | the step "`mu` ranges over R^(3n)" |
| E2 | a constraint that `t` satisfies and `t + mu` does not | a property of the answer, not of the data | the step "`t` ranges over R^(3n)" |
| E3 | a second observation whose distribution depends on `(t, mu)` **differently** | information generated by a different process | the likelihood factorisation |

**Nothing else can work.** Any estimator built from `{w_k}` alone -- consensus, medoid, typicality,
dispersion, any weighting, any re-ranking, any selection, any re-averaging -- is a function of the
identifiable part and is therefore *constant in `mu`*. That is the formal content of the memory
note's prose, and it is why "every operator downstream of retrieval is capped".

---

## 4. THE ESCAPES ARE NOT HYPOTHETICAL -- THIS PROJECT HAS MEASURED ALL THREE

This is what makes the theorem worth having: the accounting closes on the record.

**E1, a prior on the bias's form -- CLOSED, and the literature says why it is dangerous.**
The natural form to posit is a scale error, because one is known to exist: coordinate averaging
**contracts the backbone 25.8%** (memory `averaging-space-beats-the-objective`), which is Jensen's
inequality and not a fittable artefact. But the memory note already derives the closure: the
optimal rescale is exactly `s* = <c,t>/|c|^2 = 1 - <c,ebar>/|c|^2`, **a function of the invisible
component and of nothing else** -- which is precisely Corollary 1 showing up as algebra, and which
*derives* the empirical placebo result that a mismatched same-length native serves the correction
just as well.

The literature that governs E1 is **Kennedy & O'Hagan (2001)** model calibration with a discrepancy
term, and its decisive follow-up **Brynjarsdottir & O'Hagan, Inverse Problems 30:114007 (2014),
"Learning about physical parameters: the importance of model discrepancy"**. Their result, in our
words: the calibration parameter and the discrepancy function are **not jointly identifiable** --
many configurations explain the data equally well -- and a discrepancy term improves the physical
parameter estimate **only in the presence of a strongly informative prior on the discrepancy's
shape**. This is E1 proved in a different field, and it carries a warning we already have in
another costume: *a wrong discrepancy prior is worse than none*, which is this project's
"confidently wrong costs 2-3x absent" (memory `torsion-restraints-reach-the-target`). Two
literatures, same finding, arrived at independently.

**E2, a constraint on the native -- OPEN, MEASURED, AND SMALL.**
A hard physical constraint satisfied by real backbones and violated by the contracted average is
exactly an E2. The project has run it: **AMBER-relaxing the average at k=30 is -0.022 A
[-0.036, -0.009] at n=126, valid geometry, 5/5 folds same sign** (memory
`averaging-space-beats-the-objective`). That is the theorem's E2 exit being taken, working, and
being worth 0.7% of baseline.

**This is the most useful single line in the file: E2 is the only escape this project has ever
gotten a signed, fold-consistent result out of, and its measured size is -0.022 A.** Anyone
proposing a geometry-repair, constraint-projection or physical-validity operator is proposing an
E2 and should be priced against that number, not against the 68%.

**E3, an outside source -- OPEN IN PRINCIPLE, AND THE RECORD SAYS WHY IT KEEPS FAILING.**
E3 requires a source whose bias is **not the same `mu`**. The project's S24 measurement is the
decisive one and it is not what one would guess:

> quality-matched **retrieval-free** provenance cosine **0.9432**, against a **0.9330**
> within-source control (memory `pool-error-is-68-percent-common-mode`, S24 corrections).

Read that carefully. Two candidates from *completely different sources* are **more aligned with
each other** than two draws from the *same* source. **`mu` is therefore not a property of where the
candidates come from.** It is a property of the prior they are scored against -- which is why the
same note records "selection aligns the output with WHATEVER prior it is scored against, including
a target-blind one", and why `prior-derivative-is-the-only-steep-lever` (-2.15 A per unit) is the
only steep direction anyone has found.

---

## 5. THE CONSEQUENCE FOR LANE X, STATED AS A CONDITIONAL

Lane X's premise -- "selection is capped by common-mode error, so consider generating instead" --
is **justified by the theorem**: selection is a function of the pool, hence constant in `mu`, hence
structurally capped. Section 4's E3 measurement then constrains the conclusion:

> **Changing the GENERATOR does not change `mu`.** A retrieval-free source shares the bias at
> cosine 0.9432, above the within-source control. The escape is not a different supplier of
> candidates; it is a different **prior**.

So the theorem justifies half of lane X's premise and refutes the obvious remedy. A generation arm
is worth running only if it can show its candidates are *scored against something else* -- not
merely *drawn from somewhere else*. That is a pre-registrable distinction and it should be lane X's
falsifier.

---

## 6. WHAT THIS FAMILY OF LITERATURE ASSUMES, AND THAT WE VIOLATE

The brief asked for findings of the form "this entire family assumes X, which our pool violates".
Here is the list, with the violated assumption named per family rather than per paper.

| family | the assumption | our violation |
|---|---|---|
| bagging / ensemble averaging / consensus QA (Pcons, ModFOLDclust, DAVIS-EMAconsensus) | member errors i.i.d. mean-zero; gain scales as 1/K | `m_eff = 1.4` out of 75; `mu != 0` |
| Krogh-Vedelsby ambiguity / Ueda-Nakano bias-variance-covariance (Brown, Wyatt & Tino 2005 eqs 9-10) | **nothing -- it is an identity and it is ours** | none; it is the correct framing and supplies the `(1 - 1/M)` covariance coefficient that does not decay |
| negative correlation learning (Liu & Yao 1999) | members are estimators **being trained** | ours are retrieved real windows |
| control variates (Glynn & Szechtman) | a control variate with **known mean** | the known mean is the native |
| multifidelity Monte Carlo (Peherstorfer et al. 2018) | `m0` samples of the **high-fidelity** quantity as an anchor | reduces variance around an anchor; never de-biases |
| Kennedy-O'Hagan calibration with discrepancy | field observations of the truth, **or** an informative prior on the discrepancy's form | E1 above; and B&O 2014 proves the requirement is not optional |
| factor models / PCA / ICA on the members | recovers the shared direction up to **sign and scale** (the standard gauge indeterminacy) | it recovers the *contrast* structure; `mu` is absorbed by the centring, so the offset is not in the data at all |
| blind source separation | a structural assumption (non-Gaussianity, sparsity, non-negativity) to fix the gauge | we have none that pins the offset |

The last two rows are worth one extra sentence, because they are the ones people re-propose. PCA on
the centred pool returns the common-mode **direction** but is *identically uninformative* about
`mu`'s component along it, because centring removes exactly that component before the
decomposition runs. **That is the reason lane O's ORACLE global `eta` for the PC1 family is
`+0.0000 exactly`** (S29-L47). An exact zero in a measurement is almost always something being
identically zero for a reason; this is the reason. Two of S29-L47's four scalars are exactly zero
and both are of this type.

---

## 7. WHAT I DID NOT PROVE, AND WHERE I COULD BE WRONG

- The theorem is **relative to model class M**. If the true error structure is not
  "common bias + i.i.d.", e.g. if the loadings vary in a way correlated with an observable, the
  invariance is weaker and something is recoverable. The measured `f` per target ranges 0.372 to
  0.961 (10th/90th), which is *consistent with* varying loadings and does not prove M.
- `m_eff = 1.398` is a ratio of two aggregate numbers from `errdecomp.json` and is arithmetic on an
  existing artefact, **not a new measurement**. The per-target version (1.479) uses the note's own
  per-target mean `f`. I have not recomputed either from the raw file.
- Corollary 2's "`O(K^-1/2)`" assumes finite second moments under G and independence of `d_k`,
  which is the very assumption the rest of the file says fails for `e_k`. It is fine here --
  `d_k` is defined as the residual after removing the common part -- but it is an assumption and
  the S24 correction ("`f` is NOT a screen") is the warning that this decomposition is descriptive.
- I have **not** shown that E3 is achievable. I have shown that the obvious E3 (a different
  candidate source) is measured to fail, and located where a working E3 would have to differ.
