# S30 LANE X — FINDINGS

Standing format: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN, plus "what
damaged my own expectations" and "what I did not do and why". Pre-registration
`s30/PREREG_S30_X.md`, committed **d4305d17, 2026-09-20 12:46**, before any number below existed.

Lane X is the permanent divergent role. I was asked whether the quantum stage should GENERATE
rather than SELECT, and handed one gate to apply before spending compute: *"does this space
contain structures better than the pool's own ORACLE best? ... Measure the ceiling first,
always."* **This lane attacks that gate, and the gate loses.**

---

## 1. DEMONSTRATED — THE CEILING GATE IS AN ANTI-PREDICTOR: SIGN CORRECT ON 1 OF 10 CELLS

Two independent families in the record where a candidate SOURCE was changed and both its ORACLE
ceiling and the endpoint were measured. Scoring each on the only thing a gate is for — **does it
predict the DIRECTION the endpoint will move?**

| family | cells | ceiling gate sign correct | source law sign correct |
|---|---|---|---|
| S7-6 K-ladder (K = 25…2000 vs K = 500) | 6 | **1 / 6** | **5 / 6** |
| S24 generated-source unions (4 samplers) | 4 | **0 / 4** | n/a (set mean not stored) |
| **total** | **10** | **1 / 10** | 5 / 6 |

A coin is right 5 times in 10. `P(X ≤ 1 | n = 10, p = 0.5) = 11/1024 = 0.0107` — **indicative
only**: the six K-rungs share one reference rung and the four S24 arms share 126 targets and one
pool, so the cells are not independent. The size of the miss is the finding, not the p-value:
mean |residual| **0.333 Å for the ceiling gate against 0.0153 Å for the source law**, a factor of
**22**, on the same six cells.

Every row on both families is ORACLE (the ceiling is a per-target minimum against the native).

**Why a wrong-signed gate is worse than no gate:** it rejects exactly the spaces it should admit
and admits the ones it should reject. S29-L56 abandoned the chimera space on this gate.

## 2. DEMONSTRATED — THE SOURCE LAW. PRIMARY H-X1 PASSES ITS REGISTERED FALSIFIER.

Coefficients **frozen** at the Sprint 20 wide-set values (2,142 (arm,target) cells, 17 generator
arms at 8,192 draws; `operator-consumes-set-mean`). Nothing is fitted in the primary arm.

```
    d_out  =  0.803 * d_set_mean  +  0.298 * d_set_best
```

H-X1 PRIMARY, `docs/FINDINGS.md` §S7-6 / `s7/poolsize_kcurve.json`, n = 126 per rung, six rungs
predicted against the shipped K = 500 reference:

```
    mean |residual|  0.0153 A        registered bar 0.10 A        PASS
    max  |residual|  0.0263 A
    sign correct     5 / 6           registered bar < 1/3 wrong   PASS
    VERDICT: PASS
```

Matched controls in the operator's own space (contract rule 6), same six cells:

```
    two-term source law     mean |resid| 0.0153     sign 5/6
    mean-only               mean |resid| 0.1198     sign 5/6
    best-only (THE GATE)    mean |resid| 0.3330     sign 1/6
```

**DISCLOSURE (contract rule 14).** One cell — K = 25 → K = 2000, predicted +0.0843 against an
observed +0.0810 — I computed by hand from the published table *before* writing the
pre-registration. It motivated the hypothesis, it is **excluded from the bar above**, and it is
carried separately in the artefact as `disclosed_cell_NOT_out_of_sample`.

**IDENTIFIABILITY, registered in advance and confirmed.** Across the K-ladder, pool mean and pool
best correlate at **r = −0.9898**. This family **cannot** identify the two coefficients
separately, so none were fitted on it; only the frozen-coefficient prediction is identified.

H-X1 SECONDARY, S24 lane C (`s24/results/c_ladder.json`), a different sprint, corpus and sampler
family, **504 (arm, target) cells**, 4 arms × 126 targets:

```
    WITHIN-TARGET slope on the selected set's mean   0.9349
      target-clustered CI95  [0.8201, 1.0632]        R2 within 0.668
    pooled (cross-target, contaminated by difficulty)  1.0413   R2 0.929
```

The within-target slope is the operator's own response with each target's difficulty removed. It
is the identified version of the law's dominant term and its CI excludes zero by a wide margin.

## 3. DEMONSTRATED — THE ADMISSION CONDITION THAT REPLACES THE GATE

Combining §2's frozen `b = 0.298` with the measured within-target `a = 0.935` (and with the
frozen `a = 0.803` as the bracket):

```
    ADMIT a generative space G against pool P iff

        Delta(set mean)  <  -(b/a) * Delta(set best)          b/a in [0.319, 0.371]
```

> **An X Å improvement in a generative space's CEILING is worth having only if its TYPICAL member
> degrades by less than ~0.32–0.37 X. At fixed ceiling, improving the mean is worth ~2.7–3.1×
> what improving the best is worth.**

It is also the cheaper gate. A set MEAN is estimated from a handful of draws with small variance;
a set BEST over 10⁴ configurations is an order statistic needing `best_of_k_within` and a
split-half transfer arm before it can be read at all (`grid-oracles-are-order-statistics`).

**SCOPE, stated so it is not over-read.** Validated for changes to the candidate **SOURCE at a
fixed score**. It is **not** validated for changes to the **SCORE at a fixed source** — the
Sprint 18 amendment to `operator-consumes-set-mean` measured the law failing on exactly that
(Legacy gate: premise held, law predicted a gain, the output got +0.076 worse). Lane F's S30-L2
(the 500→75 filter is +1.767 Å worse than a random 75 on FAIL18) is a **gate** change and this
condition must NOT be composed with it.

### 3b. The design specification, stated so it can be falsified

Best ceiling available anywhere on this corpus is the whole-library best **1.313 Å**
(`docs/FINDINGS.md` §S8-4, `s8/triage.py`) against the shipped pool's **1.7108 Å** (ORACLE both).
Granting a generative space that ceiling *in full*, `Δbest = −0.398` contributes only
`0.298 × −0.398 = −0.119 Å`. To reach the sprint's −0.30 Å primary:

```
    Delta(set mean)  =  (-0.30 + 0.119) / a   =   -0.19 A  (a = 0.935)  to  -0.23 A  (a = 0.803)
```

> **A winning generative space must make its TYPICAL member ~0.2 Å better than the pool's typical
> member even after being granted a perfect ceiling. Generation must raise the space's FLOOR.**

Nothing in the record approaches this. S24's four samplers: standalone 3.789 / 3.243 / 3.207 /
3.175 against the 3.0483 incumbent point cloud — every one is worse in the mean while three of
them beat the pool at the ceiling.

### 3c. Why this inverts the case for a *quantum* generative space specifically

A `2**q` register over structural variables is by construction a **wide** space (S29-L56's chimera
was `8**S = 32,768`). Width is the encoding's selling point. But width at ρ ≈ 0 buys ceiling and
costs mean — it is the K-ladder move run further. **The property a quantum encoding is chosen FOR
is the property the admission condition penalises.**

The coefficients are not physics; they are a property of the **readout**. Under an argmin readout
they would be (0, 1) and the ceiling *would* be the endpoint. So:

> **Wide space + averaging readout consumes the wrong statistic. Wide space + argmin needs ρ,
> which is closed. Narrow, uniformly-good space + averaging readout needs no ρ at all — and is
> the only one of the three that is open.** Generation is closed *jointly with the readout*, and
> the joint closure is tighter than either half.

## 4. REFUTED (my own hypothesis, as I registered I expected) — THE NMR-ENSEMBLE FLOOR

`core/geometry.py:808` scores against **model 1** (`model_index=0`); `ca_rmsd_to_ensemble` and
`native_ensemble_from_pdb` exist in two modules and are **called by nothing**. 111 of 126 tuning
targets are multi-model depositions (mean 15.6 models, median 20). So the question was real.

The floor, ORACLE by construction: a method predicting the ensemble MEDOID perfectly still scores
`RMSD(medoid, model 1)` against the benchmark.

```
    floor mean 0.6136 A   median 0.3784   sd 0.7630   max 4.2943   21.4% above 1.0 A
    pairwise ensemble spread mean 0.9535 A
    endpoint (built chain) 3.2148     FAIL18 6.2872     other 108 2.7027
    corr(floor, endpoint) +0.1220     corr(spread, endpoint) +0.1009
```

**The tail enrichment has the OPPOSITE sign to the hypothesis.** FAIL18 minus the other 108:

```
    d endpoint   +3.5845 A
    d floor      -0.1776 A     share of the tail excess  -5.0%
    d spread     -0.4872 A     share                    -13.6%
```

**VERDICT: PRICED AND DEAD at this endpoint** (registered rule: dead if mean floor < 1.0 Å and it
explains < 25% of the FAIL18 excess). My registered prediction was ~0.6 Å mean and < 5% of the
tail; measured 0.614 Å and −5.0%.

**And it answers the coordinator's own open question.** *"Does the tail's difficulty have a cause
we can name, or is it just harder sequences?"* — **it is not conformational ambiguity in the
reference.** The hard 18 have *tighter* deposited ensembles than the easy 108. If anything the
tail is the RIGID end of the set. This is consistent with lane F's S30-L2 (the tail is not
pool-limited either) and narrows the cause further.

**Calibration, so the number is reusable:** the model-1 convention costs ≤ ~0.61 Å mean and does
not bind until the endpoint is below roughly 1.5 Å. Revisit it there, not before. Note
`corr(floor, pool_best_ORACLE) = +0.311` — a floppy peptide has a floppy retrieval neighbourhood,
which is the mechanism to expect if it is ever revisited.

## 5. A GAP IN A CLOSURE THAT HAS BEEN DRIVING LANE DESIGN FOR SIX SPRINTS

`prior-derivative-is-the-only-steep-lever` instructs "do not build another candidate generator —
closed on five instruments". Checking the artefacts: four of the five measured **only the endpoint
after selection** and carry no oracle arm at all (`s24/results/biasalign.json`, `qmatch.json` —
field lists verified). The fifth (S24 lane C) measured the generated source's oracle best in
isolation but **never formed `min(pool, generated)`**, the enlarged pool's actual ceiling.

Recomputed here from lane C's own stored per-target rows (`s30/s30_X_sourcelaw.py`), ORACLE:

| arm | targets (of 126) carrying a structure better than the WHOLE K=500 pool's best | union ceiling, matched 500 | Δ vs 1.7108 | union ceiling, all 2000 | Δ | endpoint Δ |
|---|---|---|---|---|---|---|
| T0_helix | 12 | 1.6969 | −0.0139 | 1.6814 | −0.0294 | **+0.0085** |
| T1_blind | 39 | 1.6088 | −0.1020 | 1.5297 | −0.1812 | **+0.0170** |
| T2_restype | **45** | 1.5946 | **−0.1162** | 1.5154 | −0.1955 | **+0.0157** |
| T3_pool | **47** | 1.6011 | −0.1098 | 1.5229 | −0.1879 | **+0.0509** |

An **untrained per-residue-type Ramachandran sampler** puts a structure better than the entire
retrieved pool's best member on **45 of 126 targets** and the endpoint still gets worse.

**This does not weaken S24's verdict — it is §1 a third time, and it is the strongest single
demonstration in the record that the ceiling is not the endpoint.** What it does change is the
wording: the closure should read *"no achievable source produces candidates the shipped score can
convert into Ångströms"*, not *"no alternative source contains better structures"*, which is false
on 45/126 targets. The stronger wording has been used to rule proposals out.

Inverting the law on these four arms, the Δ(set mean) each union must have suffered to produce its
observed endpoint is **+0.016 / +0.059 / +0.063 / +0.104 Å** — small, positive, and ordered exactly
like the generated share of the union's top-75 (0.097 / 0.235 / 0.256 / 0.523). Consistent, and
labelled an inversion rather than a prediction because the union's set mean is not stored.

## 6. THE DIVERGENT VERDICT, WHICH IS NOT THE ONE I WAS SENT TO FIND

I was asked to attack "retrieve, then select or average." The honest answer from the measurements
above is that **the frame is right for a reason nobody in this project has written down.**

At ρ ≈ 0 the endpoint is a near-unit-slope function of the candidate set's MEAN
(0.935 [0.820, 1.063] within target, 504 cells). The **only** stage in the pipeline that raises
that mean is retrieval — it turns a per-target universe of ~18,674 windows into a K = 500 set of
mean 4.4533 Å. No downstream operator can do it, and no generative sampler in the record does it
either. **Replacing retrieval with generation replaces the one stage that optimises the statistic
the terminal actually consumes.**

So my recommendation is the opposite of my remit: **do not build a generative structural space
this sprint.** The specification in §3b is a −0.2 Å improvement in the typical member, and the
gate that would have let a proposal through is the one shown in §1 to be wrong on the sign 9 times
in 10.

**Where the live freedom actually is:** in `a` and `b` themselves. They are readout properties.
Everything above is a statement about a coordinate average over a score-ordered top-m. It is the
same place S30's STATE already points ("a sparse weighted readout — the only ladder class not
closed by ceiling"), reached from the opposite direction and with a transfer function attached.

## 7. A FRAMING I DO NOT ACCEPT AS BINDING, AND WHY

My brief cites *"the pool's error is 68% common-mode — an exact identity, 50.7× the i.i.d.
prediction"* as the brutal cap that motivates generation. Reading the memory's **body** rather
than its index line (`read-the-memory-body-not-the-index-line`), it records two things that
weaken it as an argument for anything:

1. `|ebar|² = n·RMSD²` **exactly** — the "common-mode" numerator *is* the output RMSD. A quantity
   algebraically equal to the endpoint is not independent evidence about the endpoint.
2. Its own S24 correction: ***"f is NOT a screen and must never be used as one"*** — a blind
   library draw scores f = 0.4896 against the incumbent's 0.6758 while being 0.76 Å **worse**.

And the "50.7×" compares that identity against an i.i.d.-member model that retrieval exists to
violate: retrieval deliberately returns structures similar to one another. The durable content is
the weaker statement — **pool members resemble each other more than any resembles the native** —
which is true, useful, and is not a cap on generation. The same memory's own conclusion is that
*"the remaining leverage is in candidate GENERATION" is REFUTED*.

## 8. MULTIPLICITY, POWER, AND WHAT THESE NUMBERS ARE NOT

- **The primary is a 6-cell prediction test against a registered RESIDUAL bar, not a powered
  effect-size comparison.** The S7-6 rungs are published aggregates; no per-target rows survive,
  so no SE, MDE or fold CI is computable for them and none is quoted. Registered as a residual
  bar for exactly that reason.
- Comparisons made by this lane: 6 K-ladder cells × 3 predictors, 4 S24 union arms, 1 within-target
  slope, 1 pooled slope, 4 per-arm slopes, and the H-X2 battery. Only H-X1-primary and H-X2 were
  pre-registered with falsifiers and only those two are read as results.
- No per-target maximum over any family is taken anywhere in this lane, so no `best_of_k_within`
  arm is required.
- Every ceiling row is ORACLE and labelled ORACLE in the sentence that uses it. No native was used
  to set any parameter; the only fitted quantity is §2's within-target slope, which is a
  descriptive coefficient, not a deployable.

## 9. WHAT DAMAGED MY OWN EXPECTATIONS

1. **I expected to be able to state a "richer space" criterion and hand it over. I could not —
   the criterion I was given turned out to be wrong-signed, and the one that replaces it says
   don't build the space.** The divergent lane's output this round is a prohibition, not a design.
2. **I expected the NMR-ensemble floor to explain part of the tail.** It explains **−5%** — the
   hard targets have *tighter* ensembles. I registered that I expected to kill it and I did, but
   I did not expect the sign to reverse.
3. **I expected S24's five-instrument closure to be five instruments.** Four of the five carry no
   oracle arm at all, and the fifth never formed the union ceiling it was closest to measuring.
   The closure is still correct; its stated reason was not the one the artefacts support.
4. **I expected the ceiling gate to be merely insufficient.** 1 of 10 on the sign is not
   insufficient, it is inverted.

## 10. WHAT I DID NOT DO, AND WHY

- **Build a generative space.** §3b sets a −0.2 Å specification on the set mean that nothing in
  the record approaches, and building before deriving is what cost S29-L56 its lane. The
  derivation was the deliverable my remit asked for first.
- **Fit (a, b) on the K-ladder** — r = −0.9898 between the regressors, registered in advance as
  forbidding it.
- **Spend any VQE or pipeline compute.** Every number here comes from stored artefacts and 126
  PDB parses. Box impact ≈ nil, which matters under contract rule 21.
- **Compose §3's condition with lane F's S30-L2.** F's is a gate change; S18 measured this law
  failing on gate changes. Naming the non-composition is the point.
- **Quote the 12-target S29-L56 probe as instrument evidence** (contract rule 16). Its chimera
  numbers appear once, as a labelled approximate inversion.

## ARTEFACTS

`s30/PREREG_S30_X.md` (d4305d17, before every number); code `s30/s30_X_sourcelaw.py`,
`s30/s30_X_ensemble.py`; results `s30/results/s30_X_sourcelaw.json`,
`s30/results/s30_X_ensemble.json` (126 per-target rows). Sources read, all pre-existing:
`docs/FINDINGS.md` §S7-6 and §S8-4, `s24/results/c_ladder.json`, `s29/results/s29_X_probe_*.json`,
`bench_results/baseline_tuning126.json`, `results/monomer_manifest.json`, `pdbs/`, `pdbs_ext/`.
