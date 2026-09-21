# S32 CONTRACT — the rules every lane works under

Carried forward from S31's 32 rules (each earned by an incident) plus S32-specific additions.
**Read this before your first measurement. Violations invalidate a result regardless of its size.**

---

## 0. The endpoint

**The only score is mean built-chain Cα RMSD on `tuning126`, n = 126. Production = 3.2105 Å.**

Three different objects, never differenced against each other:

| object | value | what it is |
|---|---|---|
| **built chain** | **3.2105** | **THE ENDPOINT** |
| CA point cloud | 3.0483 | an intermediate, diagnostic only |
| set mean | 3.5507 | the mean over the top-75 members, a third object |

**Every number you report names its basis.** A number without a basis is not a result. S31 had
exactly one misstated basis and it was in the headline: it understated its own headroom by 0.16 Å.

## 1. Statistics

- `MDE = 2.8016 × SE`, **per comparison**, never a shared constant. Report SE beside every mean.
- `s24/stats_lib.compare(v, base, folds=..., names=...)` is **LOWER IS BETTER**. A positive effect
  is a regression. It refuses a verdict without `folds`.
- **< 0.7× MDE → NOT A RESULT.** **0.7–1.0× → NOT MEASURED** (say so; it is not a negative).
  **≥ 1.0× with fold CI excluding zero and ≥ 4/5 folds agreeing → RESULT.**
- A fold CI excluding zero while the effect is under MDE does **not** rescue it. The MDE gate binds.
- Report **median beside mean** and W/L. A near-even W/L with a CI excluding zero means
  concentration — compare against a uniform-effect null before claiming a broad effect.
- **Paired, in the same job**, wherever the chain is involved (see rule 3).

## 2. ORACLE discipline

Any arm that reads native coordinates — at any point, including feature choice, thresholds,
hyperparameters, architecture selection, or stopping — is **ORACLE / NOT DEPLOYABLE** and must be
labelled on **every** occurrence, in tables, prose and artefacts. An ORACLE result prices a ceiling.
It is never a breakthrough and never a deployment candidate.

## 3. The projection is deterministic but chaotic — and this now matters scientifically

There is **no RNG** on the projection path. The λ=0 multi-start argmin is decided at a **1e-7**
spread among branches **1e-1 apart in RMSD** — amplification ~1e13; **73/126** targets sit below a
1e-6 margin.

- **Project both sides of a chain contrast in the SAME job.** Cross-job chain comparison is legal
  only if you *check* that the shared `PROD` rows are bit-identical, and you print the check.
  (S31 verified 126/126 bit-identical across two jobs; the rule is "check", not "never".)
- **The per-target built-chain floor is 0.0134 mean / 0.0329 p90 / 0.2285 max**, measured between
  two implementations of the *same* operator. **0.0107 is a MEAN floor and does not bound
  per-target statements.**

## 4. A number carries its definition, not just its value

Every defect in S31's audit — all of them — was one shape: a quantity re-used across a boundary its
definition does not cross (an object, a search size, a bit budget, a basis, a data regime, a moment,
a key, a sign). **When quoting a number from another sprint, open the entry that produced it and
copy its definition into the sentence.**

The S31 instance that nearly reached the headline: `CTRL_SHRINK_ORACLE_Y` is shrunk to the *perp*
arm's norm; differencing the *along* arm against it reads an apparently controlled +2.06× positive
and is meaningless.

## 5. A verification must be able to fail

If the test data cannot exhibit the thing the caveat is about, the check is decoration. S31
"verified" a claim about ties on **random floats, which never tie**; on real data 125/126 targets
tie. Ask: *what input would make this fail?* If there is none, the test is wrong.
**Audits ship with a self-test on the defect that motivated them.**

## 6. Controls must match the operator's own space

The project's most repeated error. A control matched to a *different* arm's magnitude is not a
control. State what the control is matched to, in the same sentence as the number.

## 7. Cross-lane claims have no owner

**In S31 every single-lane result held and every cross-lane synthesis failed, four for four.** A
claim spanning two lanes is audited by nobody, because each lane sees half and assumes the other was
checked. **Any claim combining two lanes' numbers must be re-derived from both lanes' raw artefacts
by one person, in one script, before it is written down.**

## 8. Pre-registration

Before the first number exists, commit: HYPOTHESIS / MECHANISM / PREDICTION / FALSIFIER / CONTROL /
TEST / STATISTICAL RULE / DEPLOYMENT CONDITION. Commit it to git and cite the commit hash in the
ledger. **A falsifier decided after seeing the outcome is not a falsifier.**

## 9. Multiplicity

Write every emitted comparison to `s32/MULTIPLICITY.md` **as it is emitted**, not at the end. Record
registered vs exploratory. **A per-target minimum over K variants is mostly an order statistic** —
price best-of-K before calling it a lead, and prefer split-half transfer.

## 10. Draw distributions

A control with a random draw needs **its own distribution**, not its best draw. Report the draw mean
and the draw-to-draw sd; a single draw that clears MDE while another does not is NOT MEASURED.

## 11. Benchmark integrity — absolute

**Never** open or alter `benchmark60`. **Never** regenerate the frozen folds, the cluster
assignments, the sealed splits, or the target identities. **Never** tune a deployable parameter on
native RMSD.

## 12. FAIL18 is outcome-defined

It is a **circular** stratum: selected by the outcome being measured. Use it as a mechanism
diagnostic only, and prefer filter-independent strata. S31 measured the effect size rising
monotonically with a stratum's circularity (0.16× → 0.38× → 0.84× → 1.24×) — *that gradient is the
stratum's definition doing the work.*

## 13. Retractions are in place, not by deletion

Annotate the original wording with the retraction beside it. Never silently edit a claim that has
been quoted. Scientific history is not erased.

## 14. Never

Fabricate a result. Alter a result to fit a hypothesis. Tune deployable parameters on native RMSD.
Hide a failed run. Silently change an endpoint definition. Confuse cloud and built-chain metrics.
Call an ORACLE result deployable. Report "didn't work" without "because".

---

## S32-specific rules

## 15. An arm that CONSTRUCTS rather than SELECTS carries the geometry secondaries from its first row

S31's `AVG_SEP` re-embedded by MDS and was the only arm carrying none of the four registered
geometry checks; its own native-free `move` column already said it failed (0.792 mean / 1.746 max)
before any endpoint number existed. **Virtual-bond mean and sd, and the displacement from the
structure it claims to repair, ship with row one.**

## 16. The projection cost is a property of the OBJECT being projected, not a constant

Measured on the S29 O-ladder, built chain minus cloud:

```
real deposited member (a valid chain)      -0.0007 to -0.0030    ~FREE
sparse convex combination, K=500, s=10     +0.0002               ~FREE
sparse convex combination, K=500, s=20     -0.0052               ~FREE
dense prefix average, K=500                +0.1701
PRODUCTION (75-member uniform average)     +0.1622
```

**Do not quote "the projection costs 0.16 Å" as a constant.** But the mechanism is **NOT distance
from the manifold** — that framing is **corrected by lane R (S32-L9)**, which measured the actual
`d` per rung in one job:

```
rung                cloud   chain      d  |  obs price   orthogonal null  |   cos
prod               2.7337  2.8725  0.7054 |   +0.1387        +0.1303      |  -0.05
dense avg (bestm)  2.3832  2.5346  0.8001 |   +0.1514        +0.1579      |  +0.03
sparse K=500 s=10  1.1375  1.1446  0.8435 |   +0.0072        +0.2832      |  +0.40
```

**Production and the sparse combinations sit at essentially the SAME manifold distance (0.705 vs
0.844) and pay +0.139 against +0.007 — a factor of 20.** What separates them is **direction**:
`cos = (e² + d² − chain²)/(2ed)`, where 0 means the displacement is orthogonal to the native error
and +1 means it removes error one-for-one.

> **The projection does not tax you for leaving the valid-chain manifold. It taxes you for leaving it
> in a direction that has nothing to do with your error. It is not sparsity that is cheap — it is
> ALIGNMENT.**

> **AND ALIGNMENT IS NOT A LEVER — lane R's n=126 control (S32-L(R1)) closes it.** The *same*
> `s = 10` combination pays **+0.0002 when the native chose its ten members and +0.1988 when it did
> not**. Same sparsity, same averaging operator, same projection, same job, three pinned draws.
> **The cheap price was ORACLE-induced.** The deployable version — score-top-10 — is **+0.1371 at
> 1.46× MDE on the WRONG side of its own orthogonal null**, and production at **+0.1622 against a
> null of +0.1504 is 0.24×, NOT MEASURED**. ***Every object the native did not touch sits at or
> above its own orthogonal null.*** So `d` and `cos` are the right *diagnostics* and neither is an
> available *intervention*.

Any readout proposal must state its **`d` and its `price`** — not merely "sparse" or "dense".

> **But `cos` is an INTERPRETATION of the price, never corroboration of it.** Lane V verified (`s32/results/s32_V_cos_identity.json`, n = 79) that rebuilding the price from `(e, d, cos)` returns the observed price with **max |error| = 0.000e+00**. `cos = (e² + d² − chain²)/(2ed)` is the law of cosines inverted on three numbers the row already carries, so it is a **bijection** with the price given `(e, d)`. ***Any sentence quoting the price AND the cos as two pieces of support is double-counting one measurement.***
>
> **It does, however, mean what its name says** — which was not free. The three RMSDs are Kabsch distances in Kendall shape space, which is curved and in which the law of cosines is *not* an identity, so the inverted quantity could have been a triangle defect rather than an angle. Measured directly with all three objects in the cloud's own frame: `cos_algebraic` **−0.0519** against `cos_direct` **−0.0557**, mean difference **0.0038** (p90 0.0085), and `|cos| > 1` on **0 of 79**. *The triangle is near-Euclidean at this scale, so `cos` is a real alignment cosine to ±0.004.*
>
> **And one number to carry: production's mean `cos` is −0.052 — on the WRONG side of the orthogonal null.** The projection's displacement is not merely uncorrelated with the cloud's error; **it points slightly away from the native.** That sharpens S32-L9 rather than weakening it.

## 17. The pool is not the bottleneck; say what is, with the ladder in the sentence

Built chain, all ORACLE / NOT DEPLOYABLE:

```
best member in the K=500 pool                1.7078
best sparse convex combination, K=500, s=10  1.1139
PRODUCTION                                   3.2105
```

**There is 1.50 Å of headroom inside the existing pool's single members and 2.10 Å inside its sparse
convex combinations.** Any claim that candidate *generation* is the bottleneck must first explain
why 2.10 Å of already-present headroom is not the bottleneck.

## 18. Direction, not correlation

An observable that captures a real and large component of the error can still point the wrong way.
S31: the 94%-predictable component of the ideal correction is worth **+0.0747 Å — harmful**. Every
proposed channel reports **directional** utility, not just correlation or R².

## 19. Reopening a closed direction requires a named mechanism

Permitted (charter §13) and encouraged — but the burden is a concrete statement of *which*
assumption of the original closure fails, quoted from the original entry. "I dislike the negative
result" is not a mechanism.

## 20. Bit-identity, not value-identity, licenses a cross-job chain comparison

Measured by lane V on the canonical λ=0.3 arm (S32-L4): **a one-ULP change in the input cloud
(7.1e-15 Å) moves the built chain by 0.10–0.15 Å.** The stage is perfectly deterministic given
bit-identical input — three in-process repeats agree to 9 dp, thread count 1/2/4/8 is irrelevant —
and it is **discontinuous** in that input. 7e-15 × the ~1e13 amplification is 0.07–0.15 Å, which is
what is measured.

- **"The same cloud value" is not enough. It must be the same float64 bits.**
- Both sides of a chain contrast are projected **in one job**, from clouds that are bit-identical.
- **Unpaired cross-job chain claims below ~0.03 Å are not resolvable.** Paired comparisons are
  untouched: SE of the chain mean is 0.1543, so an unpaired MDE against production is **0.4324 Å**.
- An independent re-projection disagrees with the S29 canonical on **126/126** targets — mean
  +0.0021, p90 0.026, **max 0.5174**. That is not a bug to fix; it is the stage's character.

## 21. The endpoint is not a cached scalar — quote its definition, not just its value

**3.210534** = *the λ = 0.3 multi-start projection arm of `s12/instrument.project` applied to the
production top-75 coordinate average, CA-RMSD to native, averaged over `tuning126`* — concretely,
the 126 `item="prod"` rows of `s29/results/s29_O_chain_rows*.jsonl`.

Five distinct objects live near it and **they are not five estimates of one** (S32-L4):

```
3.048338  CA point cloud, unprojected                       rmsd_avg
3.204076  lambda = 0 arm, no Ramachandran penalty           rmsd_fit      (a DIAGNOSTIC)
3.214765  lambda = 0.3 arm -- what PRODUCTION ITSELF EMITS   rmsd_arm / ca
3.235460  the same chain after AMBER relaxation             rmsd_full     (AMBER costs +0.0207)
3.210534  THE CANONICAL ENDPOINT -- a RE-PROJECTION of the stored cloud, lambda = 0.3
```

**The canonical endpoint is 0.0043 Å better than the chain the production pipeline emits.** A reader
who assumes the endpoint is what the pipeline outputs is wrong by more than several historical
claims are large.

## 22. The Ramachandran prior is already active, and a penalty is not a selector

The λ-ladder is `(0.0 → 0.3)` and the **canonical arm is λ = 0.3**, so the prior is already choosing
among near-degenerate solutions — `lam_path`'s docstring states that as its purpose. **Do not
propose "use the Ramachandran prior to pick the branch" as a new lever; it ships.**

What is *not* shipped, and remains open, is **post-hoc ranking of converged branches** by
Ramachandran likelihood, which is a different operator from a penalty term inside the objective.
And `λ`, the number of starts, `maxiter`, the penalty and the tolerance are **science, not tuning
knobs** — varying them is a different pipeline and must be keyed as one; tuning them on native RMSD
violates rule 11.


## 23. `split_half_transfer` nulls the wrong question unless you give it production as the baseline

The contract and every lane brief point at `stats_lib.split_half_transfer` for any best-of-K arm, and
its docstring says it *"nulls itself and needs no null at all."* **That is true of the question it
answers — does *which* setting wins transfer? — and false of the question a lane usually wants — does
the winner beat PRODUCTION?**

`split_half_transfer(M - prod)` centres on `M.mean(1)`, **the mean over the grid's columns**. If any
column is implausible, the "transfer" is measuring *"the chosen arm beats the average arm, one of
which is catastrophically bad"* — and nobody would deploy the average of the grid.

**Lane V caught this in its own code before it reached the report.** On lane R's 16 criteria, one
column is `typicality` at **+0.3434**, and the wrong version read:

```
oracle-over-criteria -0.1168, SPLIT-HALF TRANSFER -0.0254, CI [-0.0356, -0.0153], 22% of the oracle
```

**A CI excluding zero and an effect nearly 3× the best single arm — and entirely an artefact of the
baseline.** Done correctly (choose `argmin` of the per-criterion effect **vs production** on one
half, evaluate that criterion **vs production** on the other, 400 repeats):

```
out-of-sample transfer vs production   +0.0004   CI [-0.0091, +0.0116]
```

***A CI centred on zero.*** **The rule:** when you use `split_half_transfer` to ask whether a
selected arm beats the incumbent, **the baseline must be the incumbent, not the grid mean** — and a
grid containing any arm nobody would deploy will inflate the number with a tight CI. This is contract
rule 6 (*a control matched to a different arm's magnitude is not a control*) in a new disguise, and
it is the same shape as S31's `CTRL_SHRINK_ORACLE_Y` in rule 4.
