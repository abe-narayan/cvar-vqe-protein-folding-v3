# PREREG S30 / LANE X — THE GATE I WAS HANDED IS THE ANTI-PREDICTOR

Written 2026-09-20, before any number produced by this lane exists. Appended, never edited,
after the first result (contract rule 12, 15). Contract: `s30/S30_CONTRACT.md`.

Lane X is the permanent divergent role. My remit was to ask whether the quantum stage should be
GENERATING rather than SELECTING, and I was given one test to apply before spending compute:

> *"Does this space contain structures better than the pool's own ORACLE best, and by how much?
> ... If your space's ceiling is worse than the pool's, no operator can recover the difference
> and you should abandon it and say so. Measure the ceiling first, always."*

**This pre-registration attacks that test.** Not the conclusion of S29-L56 — the test that reached
it. The record already contains a direct sweep of exactly this variable, and over that sweep the
ceiling and the endpoint move in OPPOSITE directions.

---

## 0. Contract rule 24 — which S29/record findings this attacks, which it accepts as binding

**ATTACKS**

- The S30 lane-X charter gate quoted above ("measure the ceiling first, always"), and S29-L56's
  reading of its own result. S29-L56 abandoned the chimera space because its ORACLE best was
  +0.6533 Å worse than the pool's. I claim that number is the *weaker* half of the explanation
  and the gate built on it selects against the wrong quantity.
- The framing of `pool-error-is-68-percent-common-mode` as an independent brutal cap. Its own body
  records `|ebar|^2 = n*RMSD^2` **exactly** and warns that `f` "is NOT a screen and must never be
  used as one". A quantity that is algebraically the output RMSD is not independent evidence about
  the output RMSD, and the "50.7× the i.i.d. prediction" compares an identity against a model
  (independent member errors) that retrieval is explicitly designed to violate. The durable content
  is the weaker and still useful statement: **pool members resemble each other more than any of
  them resembles the native.**

**ACCEPTS AS BINDING (my whole argument is downstream of these)**

- `objective-does-not-rank-the-native` / `nothing-ranks-within-the-pool` / S29's recognition
  closures: the native-free score's ordering skill is ~0. Everything below is a statement about
  what happens to a candidate set *when rho is approximately 0*, and is void if rho ever becomes
  large.
- `operator-consumes-set-mean`, **wide-set coefficients** (Sprint 20, 2,142 (arm,target) cells,
  17 generator arms at 8,192 draws each): `d_out = 0.803*d_set_mean + 0.298*d_set_best`, R^2 0.960.
  I use the WIDE-set coefficients throughout, never the narrow-set 1.16/0.04, because every
  generative space under discussion has 10^3–10^5 members.
- `grid-oracles-are-order-statistics`: every per-target minimum over a family is priced as one.

---

## 1. THE DERIVATION (written before the measurements in section 3)

### 1.1 The record already swept the ceiling, and the endpoint went the other way

`docs/FINDINGS.md` §S7-6, `s7/poolsize.py`, `s7/poolsize_kcurve.json`, same 126 targets
(verified in place by `sed -n '1750,1795p' docs/FINDINGS.md` before this file was written):

| K | selected | pool best (ORACLE) | pool mean |
|---|---|---|---|
| 25 | 3.439 | 2.350 | 4.177 |
| 50 | **3.399** | 2.177 | 4.232 |
| 100 | 3.425 | 1.970 | 4.316 |
| 250 | 3.461 | 1.808 | 4.393 |
| 500 | 3.454 | 1.711 | 4.453 |
| 1000 | 3.483 | 1.568 | 4.522 |
| 2000 | 3.520 | 1.504 | 4.596 |

Enlarging the candidate source improves the ORACLE ceiling **monotonically by 0.846 Å** and makes
the endpoint **monotonically worse by 0.081 Å**. S7 stated the direction ("selection tracks the pool
MEAN not its BEST ... *anti*-correlated with the pool BEST"). What it did not do is turn that into a
gate, and the project has since run at least two generative probes (S24 lane C, S29 lane X) whose
admission test was the ceiling.

**So the handed gate — "if your space's ceiling is worse than the pool's, abandon it" — is, on the
one family where the variable has actually been swept, an anti-predictor of the endpoint.**

### 1.2 What richness must mean instead: an admission condition with numbers in it

Take the wide-set operator law as the transfer function from a candidate SOURCE to the emitted
structure, holding the downstream (score, top-m, coordinate average) fixed:

```
    d_out  ~  a * mean(source)  +  b * best(source) + const,      (a, b) = (0.803, 0.298)
```

A generative space G is admissible against the shipped pool P iff `Δd_out < 0`, i.e.

```
    0.803 * [mean(G) - mean(P)]  +  0.298 * [best(G) - best(P)]  <  0

                    Δmean  <  -0.371 * Δbest                              (ADMISSION)
```

Read it in the direction that bites: **an X Å improvement in a generative space's ceiling is only
worth having if the space's TYPICAL member degrades by less than 0.371 X.** Equivalently, at fixed
ceiling, improving the mean is worth **2.7×** what improving the best is worth.

This is the gate I propose in place of the one I was given, and it is strictly cheaper: `mean` is
estimated from a handful of draws with tiny variance, while `best` over 10^4 configurations is an
order statistic needing `best_of_k_within` and a split-half transfer arm before it can be read at all.

### 1.3 Why this inverts the case for a quantum generative space specifically

A `2**q` register over structural variables is, by construction, a **wide** space: S29-L56's chimera
was `8**S = 32,768` configurations. Width is the encoding's entire selling point. But width at
rho ~ 0 buys ceiling and costs mean — it is exactly the K-ladder move, run further. **The property a
quantum encoding is chosen FOR is the property the ADMISSION condition penalises.**

The complementary statement, which is where the freedom actually is: the coefficients `(a, b)` are
not physics, they are a property of the **readout**. Under an argmin readout they would be `(0, 1)`
and the ceiling would be the endpoint. So:

> **Wide generative space + averaging readout consumes the wrong statistic. Wide space + argmin
> readout needs rho, which is closed. Narrow, uniformly-good space + averaging readout needs no rho
> at all — and is the only one of the three that is open.**

Generation is therefore closed **jointly with the readout**, not separately, and the joint closure
is tighter than either half. A generative stage is admissible only if it arrives with a readout
change, and a readout change is admissible only with discrimination.

### 1.4 The corollary that defends the frame I was told to attack

At rho ~ 0 the only stage in the pipeline that raises the **mean** of the candidate set is
retrieval: it converts a per-target universe of ~18,674 windows into a K=500 set of mean 4.453.
Nothing downstream can do that, and a generative space built from torsion priors does not do it
either (S24: T1_blind 3.2435, T2_restype 3.2065, T3_pool 3.1752 standalone, all worse than the
3.0483 incumbent while their ORACLE bests are better).

**So "retrieve then average" is not an accident of history to be escaped. It is the only
architecture in the record that optimises the one statistic the terminal consumes.** I was asked to
attack that frame; the honest divergent answer is that the frame is right for a reason nobody has
written down, and the leverage is not in replacing retrieval but in the readout's coefficients.

### 1.5 The design specification this yields, stated so it can be falsified

Absolute best ceiling available on this corpus is the whole-library best, **1.313 Å**
(`docs/FINDINGS.md` §S8-4, `s8/triage.py`), against the shipped pool's 1.711. Granting that in full,
`Δbest = -0.398` contributes `0.298 * -0.398 = -0.119 Å`. To reach the sprint's -0.30 Å primary:

```
    Δmean  =  (-0.30 + 0.119) / 0.803  =  -0.225 A
```

> **A winning generative space must make its TYPICAL member 0.225 Å better than the pool's typical
> member, even after being granted a perfect ceiling. Generation must raise the space's FLOOR, not
> its ceiling.**

### 1.6 A gap in the S24 closure, found while deriving the above

`prior-derivative-is-the-only-steep-lever` instructs "do not build another candidate generator —
closed on five instruments". Checking the artefacts behind that closure: four of the five measured
only the ENDPOINT after selection and contain no oracle arm at all. The fifth (S24 lane C,
`s24/c_ladder.py`, `s24/results/c_ladder.json`) did measure the generated source's oracle best in
isolation, but **never formed `min(pool, generated)` — the enlarged pool's actual ceiling.**
Recomputed from lane C's own stored per-target rows, an untrained per-residue-type Ramachandran
sampler (T2_restype) contains a structure better than the *entire* K=500 pool's best member on
**45 of 126 targets**, moving the union ceiling 1.7108 → 1.5946 (matched 500) / 1.5154 (all 2000).

The endpoint still did not move (+0.0157). **That is not a weakening of S24's verdict — it is
section 1.1 again, and the third independent family in which a better ceiling buys nothing.** It
does mean the closure should be stated as *"no achievable source produces candidates the shipped
score can convert into Ångströms"* rather than *"no alternative source contains better structures"*,
which is false on 45/126 targets. Recorded because the stronger wording has been driving lane
design for six sprints.

---

## 2. HYPOTHESES AND FALSIFIERS — registered before the numbers

### H-X1 — THE SOURCE LAW. Out-of-sample, unfitted, no free parameters.

The S20 wide-set coefficients, applied to a candidate SOURCE's `(mean, best)` with the shipped
downstream held fixed, predict the endpoint of every enlargement/generation family in the record.

- **Arms:** the six S7-6 rungs against K=500; S24's four generated-source unions; S29-L56's chimera.
- **Prediction:** `Δd_out_pred = 0.803*Δmean + 0.298*Δbest`, coefficients FROZEN at the S20 values.
  Nothing is fitted. This is a prediction test, not a regression.
- **FALSIFIER (registered):** H-X1 is REFUTED if the mean |residual| over the available cells
  exceeds **0.10 Å**, or if the predicted sign is wrong on **≥ 1/3** of cells.
- **DISCLOSURE, because it matters and rule 14 governs chronology.** Before writing this file I
  computed ONE cell of this ladder by hand from the published S7-6 table: K=25 → K=2000 predicts
  `0.803*(0.419) + 0.298*(-0.846) = +0.084` against an observed `+0.081`. That cell is therefore
  NOT out-of-sample and is excluded from the falsifier arithmetic; it is what motivated the
  hypothesis. Every other cell is registered here before I have computed it.
- **IDENTIFIABILITY, stated before the result so it cannot be quietly dropped.** Across the K-ladder
  `mean` and `best` are almost perfectly anti-correlated, so this family **cannot identify `a` and
  `b` separately** and I will not fit them on it. Only the frozen-coefficient PREDICTION is
  identified. I will report the mean/best correlation as the evidence for that limitation.
- **What a pass buys:** the ADMISSION condition of §1.2 becomes a measured transfer function rather
  than an argument, and every future generative proposal is gated on two cheap numbers.
- **What a failure buys:** the gate I am proposing is wrong too, and the ceiling test survives by
  default. Either is a result.

### H-X2 — THE DEPOSITED NATIVE IS ONE DRAW FROM AN ENSEMBLE (the information-ceiling probe)

`core/geometry.py:808` `parse_pdb(..., model_index=0)` — the benchmark's "native" is **model 1**.
`ca_rmsd_to_ensemble` and `native_ensemble_from_pdb` exist in both geometry modules and are called
by **nothing** in the pipeline (verified by grep before this file was written). Spot-checking
`pdbs/`, the files carry 20–55 deposited models; `results/monomer_manifest.json` confirms all 52
monomer targets are SOLUTION NMR. If the deposited natives are ensembles, RMSD-to-model-1 has an
irreducible floor equal to the ensemble's own spread, and no method with any information beats it.

- **Arms:** per-target ensemble spread on tuning126 (mean pairwise CA-RMSD over deposited models,
  and RMSD of model 1 to the ensemble medoid); correlation with production per-target RMSD;
  enrichment in the FAIL18 tail.
- **FALSIFIER / decision rule (registered):** the hypothesis is **PRICED AND DEAD at this endpoint**
  if mean spread < 1.0 Å AND spread accounts for < 25% of the FAIL18 tail's excess over the other
  108. It is **LIVE** otherwise.
- **MY OWN REGISTERED PREDICTION, so this cannot be read as a win either way.** I expect mean spread
  ~0.6 Å and max ~1.5 Å (that is the monomer52 distribution, which I looked at while assessing
  feasibility and disclose here), and I expect it to explain **< 5%** of a 6.28 Å tail. **I expect
  to kill my own hypothesis**, and the value of the arm is that it prices a plausible ceiling route
  in one cheap measurement instead of a sprint. If it is dead it should be revisited only once the
  endpoint is below ~1.5 Å.
- **Rule 11 note:** `results/monomer_manifest.json`'s `ensemble_spread` field is deposited metadata
  used for SET SELECTION, not a method result; reading it does not spend benchmark60.

### WHAT WOULD MAKE ME WITHDRAW THIS WHOLE LANE

If H-X1 fails its falsifier, §1.2–§1.5 are an argument with no transfer function behind them and the
ADMISSION condition must not be quoted. I will say so in the heading of the ledger entry.

---

## 3. WHAT I WILL NOT DO

- Build a generative space. §1.5 sets a specification (`Δmean ≤ -0.225 Å` even at a perfect ceiling)
  that nothing in the record approaches, and the derivation's whole point is that building before
  deriving is what cost S29-L56 its lane.
- Fit `(a, b)` on the K-ladder — §2 H-X1 registers the collinearity that forbids it.
- Quote the 12-target S29-L56 probe as instrument evidence (contract rule 16).
- Use the native to tune anything. Every ceiling row here is labelled ORACLE in the sentence that
  uses it.
