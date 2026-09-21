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

**Do not quote "the projection costs 0.16 Å" as a constant.** It costs ~0 for anything near the
manifold of valid chains and ~0.16 for a dense average, which is not a valid chain. Any readout
proposal must state which regime it lands in.

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
