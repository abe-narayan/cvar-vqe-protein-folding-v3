# S31 — RUNNING STATE

Maintained per charter §22: current leading hypothesis, the strongest live objection to it, its
pre-registered falsifier, what each lane is doing, and what is closed and why. Notes are
newest-first below the headline. **With 7 lanes an unwritten research state is how contradictory
conclusions coexist unnoticed.**

**SPRINT CLOSED 2026-09-21 01:27. Report `s31/REPORT_S31.md` is FINAL (1,518 lines, 22 sections
+ 2 appendices); verifier `s31/s31_verify.py` at 271 matched / 0 mismatched / 0 not-found /
0 flagged. Published: <https://claude.ai/artifact/BJ6WTXG3cjtq5jQ95LqgSS> — rebuild with
`python s31/build_report_page.py`, which renders `REPORT_S31.md` client-side so the page and
the repository cannot drift. All nine lanes handed back. Do not respawn them.**

---

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

**Value is not capped the same way.** ORACLE argmin-over-128 is **2.1435 Å on the BUILT CHAIN**
(the CA-cloud value is 2.1458) against production's **3.2105 Å** — **1.067 Å of headroom**
(**ORACLE / NOT DEPLOYABLE**). So the cap is on *information*, not on *value*.

> **CORRECTED (2026-09-21 00:43, lane V's D3).** This sentence originally read *"2.1435 Å (CA cloud) against
> production's 3.0483 — ~0.90 Å of headroom."* **2.1435 is the built chain, not the cloud**, and it
> was differenced against the cloud production — so the headline **understated its own headroom by
> 0.16 Å** in the sentence carrying the framing claim. NOTE 6 has the same five numbers labelled
> correctly, so this was **synthesis error, not lane error**. And it exposes a gap in the
> instrument: `s31_verify.py` audits for an **unstated** basis and reports "none", but an audit for
> an unstated basis **does not catch a misstated one**, and the sprint had exactly one of the
> latter, in its headline. The basis audit should assert each number against the value on the basis
> it names.

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
| **F** | terminal operator; **AVG_SEP** (separation-profile correction, zero free parameters) is the live arm; now also the **per-target prefix length** | **CLOSED. `AVG_SEP` is REFUTED at 126/126: +0.4609 Å WORSE on the chain** (2.34× MDE, 5/5 folds, 41W/85L), +0.4483 on the cloud. *The “0.4–0.6 Å better” in this cell was a PARTIAL-ROWS reading with the SIGN INVERTED — `stats_lib.compare` is lower-is-better, and 0.4609 with fold CI [0.378, 0.565] is exactly that band. The lane published no aggregate; the inverted reading was mine.* |
| **L** | literature, permanent | **§L1 closed the deployed objective** (NOTE 5); continuing on L2, the new-observable question |
| **P** | the `p*` substitution — lane L's decisive experiment | running, eighth slot |

**Eight lanes, the charter's maximum.** The adversary slot is spent on lane P because its
experiment has no third outcome; the lanes have in practice been adversarial to each other —
lane C falsified my R1 quantifier **in shipped code**, lane B found a likely confound in my own
published 32-cost sweep, and lane E killed my opening hypothesis with an exact identity.

---

## NOTE 15 (2026-09-21 01:20, coordinator, S31-L22 controls): **ONLY THE *NEGATIVE* HALF OF THE ALONG/PERP RESULT IS CONTROLLED, AND THE LANE SAID SO AGAINST ITS OWN INTEREST**

Lane E's energy-matched direction control and shrink curve landed after it handed back. They were
built so a null would refute §20.1. **A null did not fire on the direction question** — at matched
energy an arbitrary direction's along-half is worth **+0.1317 more** (2.07× MDE, 5/5) and its
perp-half **−0.3797 less** (2.53×, 5/5). **`mu` is special on both sides.**

**But the two halves are not equally earned, and the report now says which is which.**
`y_along(mu)` is **84% of `y`'s energy**, the shrink curve **brackets** that norm
(0.75·3.70 = 2.78 < 3.12 < 3.33 = 0.90·3.70), and against both brackets the along arm is **NOT
MEASURED** — −0.0614 at **0.94×** and −0.0405 at **0.64× (65W/61L, a coin flip)**.

> ***"Correcting along `mu` is the whole prize" is, at matched magnitude, the same statement as
> "most of a perfect correction is the whole prize."*** The controlled half is the **negative** one:
> the orthogonal complement is *specifically* worthless (**+0.6169, 3.25×** against its own
> norm-matched shrink; **+0.3797, 2.53×** against energy-matched). **And the negative half is
> exactly what charter §14 and contract rule 29 proposed to build on.**

**Two things this changes.** §21 now states **how big the ask is** — *"supply the common mode" is
not a clever projection extracting a large gain from a small piece of information; it is a request
for most of a perfect prior* — and §0.2 no longer claims more than is controlled. What survives
undiminished is the part that is an **identity, not a fit**: the 94%-predictable orthogonal
component is worth **+0.0747 Å, worse than doing nothing**. ***The corrector is not weak; it is
aimed at the wrong component.***

**A trap recorded because I walked up to it.** `CTRL_SHRINK_ORACLE_Y` is shrunk to the **PERP**
arm's norm (1.99), **not** the along arm's (3.12). Differencing the along arm against it reads
**−0.2680 at 2.06×, 94W/32L — an apparently controlled positive**, and it is meaningless. One
boundary from being quoted as the headline; the sprint's signature defect (**a number carries its
definition**) in its purest form. `s31_verify.py` now asserts the trap's value *and* that the report
refuses it.

**Instrument fact.** These arms ran in a different process from the split they are differenced
against — legal because the two jobs' `PROD` rows are **bit-identical on 126/126 in both bases**.
**The projection is deterministic; lane D measured sensitivity to *differing* inputs, not
nondeterminism.** Verifier **230/230, 0 mismatched**.

---

## NOTE 14 (2026-09-21 00:47, lane P, S31-L20): **THE DECISIVE SUBSTITUTION IS A NULL AT 0.19x MDE — AND "THERE IS NO THIRD OUTCOME" WAS ITSELF THE HYPOTHESIS THAT DIED**

Lane L handed the sprint a clean dichotomy: replace `run_cvar_vqe`'s `p` with the closed-form
global minimiser `p*`, and *"endpoint improves -> the quantum layer is a softmax; endpoint
worsens -> the circuit's inability to optimise is the active ingredient. There is no third
outcome."* **Measured at n = 126 on the built chain, in one process from one cache, with arm A
re-projected to 3.2105 against the canonical 3.2105 at a worst per-target deviation of 0.0000 A:**

```
P1  E - D   p* vs VQE p, CONVEX readout (shipped)   -0.0112 A  SE 0.0207  MDE 0.0580  0.19x  NULL
P2  C - B   p* vs VQE p, SELECTION readout          -0.0280 A  SE 0.0348  MDE 0.0975  0.29x  NULL
```

**Neither branch fires.** The closed form is certified optimal on the real instrument (duality
gap 3.56e-09, KKT at `s*` 2.67e-15, mirror descent never below it, **circuit strictly worse
126/126**, F gap +0.1937, `p*` at 0.0012 s against the circuit's 0.0985 s) — and solving the
deployed objective *exactly* instead of *badly* is worth **nothing at the endpoint**.

**And the null is not "no effect."** Against the implementation-noise null measured in the same
job (S31-L18): `P1`'s per-target `|d|` is **0.1432 mean / 0.3706 p90 / 1.1457 max, above the
chain floor on 100 of 126 targets** — **10.7x the noise null's 0.0134**, with **74 of 126** above
the null's own p90. **The substitution changes the answer nearly everywhere and the changes
cancel.** `P2` is lane L's own caveat firing as written: the selection readout is
piecewise-constant, so `TV = 0.378` moves the medoid on **66 of 126** targets and is identically
zero on the other 60; on the 66 the effect is -0.0534 A at 0.29x, still a null.

**Why it had to be null, and this is S31-L17 arriving at the endpoint.** `E = _zrank(dis[o])` is
the standardised rank of 128 values, so it is the **same vector on every target** (max deviation
from the tie-free reference 4.07e-02). Both `p*` and `p_theta` are therefore **one fixed
weighting curve per alpha** — `H(p*) = 4.9137` bits at `alpha = 1` with **sd 3.1e-04** across
126 targets. **P1 compares two points inside a family none of whose members knows which target
it is being applied to.** An aggregate gain was never available from moving within it.

> **Consequence for the sprint, and it is a scope fix, not a result.** This experiment **does not
> test S20's law** ("optimising a bad objective harder makes things worse"). It cannot: both arms
> are target-independent rank curves, so neither carries the information S20 is about. **Do not
> record S31-L20 as evidence for or against S20.** What it does establish is that **the
> optimisation quality of the deployed CVaR-VQE is invisible at the endpoint**, which is what
> charter §11 needed.

**The row this project did not have.** `S3 = D - A`, the **built-chain** cost of the **shipped**
quantum stage: **+0.0175 A, SE 0.0180, MDE 0.0504, 0.35x — a NULL.** On the CA cloud +0.0178,
reproducing lane A's independent +0.0178 to four decimals on a different code path. It
decomposes exactly: `+0.0003` code path `+ 0.0307` prefix widening 75->128 `- 0.0134` weights.
**The largest component is the widening that exists only to feed the quantum register, and it is
a cost.** The **+0.2260 A parity deficit the sprint was carrying is withdrawn** — it was the
affine readout, not the shipped one.

**Two corrections to lane L's L1.1, both confirmed by the coordinator.** (1) Its 12 verification
cells ran at `T = 0.1, 0.05`; the **deployed `T` is 0.3 on every fold** (`core/pipeline.py:113`).
(2) Its *"the closed form is always more entropic than what the circuit finds"* **reverses at the
deployed `T` on 3 of 5 folds**: at `alpha = 1`, `H* = 4.9137` against the circuit's `5.6738`.
**That also falsified my own registered prediction that arm E would sit between D and F** — it
sits below both (E 3.2169, D 3.2281, F 3.2415) — and I registered 20-45 medoid disagreements
when the answer is 66. Both recorded in S31-L20 §6.

**Endpoint table (BUILT CHAIN, n = 126; CA cloud in brackets — a different object):**

```
A  quantum OFF (production)  3.2105 +- 0.1540  [3.0483]   D  ON, VQE p, CONVEX  3.2281 +- 0.1558  [3.0661]
B  ON, VQE p,  SELECTION     3.3117 +- 0.1648  [3.3135]   E  ON, p*,    CONVEX  3.2169 +- 0.1526  [3.0605]
C  ON, p*,     SELECTION     3.2838 +- 0.1586  [3.2852]   F  uniform,   CONVEX  3.2415 +- 0.1528  [3.0532]
A2 code-path control         3.2108 +- 0.1542  [3.0483]
```

The only measurably worse arm is the **selection readout, and only on the CA cloud**:
`B - A = +0.2652 at 2.15x` and `C - A = +0.2368 at 1.82x` (MEASURED, cloud); the same contrasts
on the built chain are 0.93x and 0.61x — **NOT MEASURED. State the basis.**

**One S32 candidate, gate first.** The decomposition isolates **`F - A2 = +0.0307 A` for widening
the retained prefix 75 -> 128** — paid **only** to fill the `2**7` register
(`core/pipeline.py:757`), and pointless when `quantum = False`, which is production.
**It is a NULL at 0.42x MDE and is NOT A RESULT.** Recorded because it is the only place this
sprint where turning something **off** has a measured sign. Its fold CI [+0.0028, +0.0557]
excludes zero while the effect sits at 0.42x — **the MDE gate binds and the CI does not rescue
it**; quoting the CI alone would repeat the sibling of the underpowered bug
(`s24/stats_lib.py:145`). **S32 should price the widening on its own, pre-registered.**

---

## NOTE 13 (2026-09-21 00:29, lane L, S31-L16): **AN AUDIT OF MY OWN SYNTHESES — AND THE STRUCTURAL REASON THEY ARE THE SPRINT'S WEAK POINT**

I asked lane L to audit my cross-lane claims and it ran the audit **ahead of the request**, on the
grounds that *a synthesis is cheapest to kill before it is written down.* It led with what holds —
lane C's S31-L5 is among the best-disciplined entries in the sprint — and then found four
**re-quotation hazards**, none of them a criticism of the lane whose numbers they use. **All four
are mine.**

### THE META-FINDING, WHICH IS THE DURABLE ONE

> **Every single-lane result this sprint HELD. Every cross-lane synthesis of mine FAILED — four for
> four.** The mechanism is **structural, not personal**: a single-lane claim is audited by the lane
> that owns the data; **a cross-lane claim is audited by nobody**, because each contributing lane
> sees only its own half and assumes the other half was checked. **The synthesis inherits both
> lanes' caveats and carries neither.**

**Rule adopted, and saved to memory:** *a claim combining two lanes' numbers needs a NAMED OWNER WHO
HOLDS BOTH, and must carry both caveats in the sentence that carries the number.* The lanes are
adversarial to each other by construction; **my cross-lane claims were the only ones in the sprint
with no reviewer.**

### (1) The lever comparison was made at UNEQUAL INFORMATION COST — the project's characteristic error, third sprint running

I wrote that *"the convex readout over the same 128 is 0.290 Å better than naming the best of them"*
and concluded the readout class is the larger lever. But:

```
best1_top128 = 2.1435 A   costs  7 BITS            (naming one of 128)
hull_top128  = 1.8538 A   costs  128 REAL NUMBERS  (S29 labels this rung [EXPRESSIVENESS], FREE WEIGHTS)
```

**So "larger lever" is not established — only "larger ceiling at unpriced cost."** S30 §9.4 already
named this exact comparison as the project's characteristic error (2 members with ORACLE *weights*
at 1.4315 against 75 with ORACLE *membership* at 2.3055, ~17.9 bits against 7), and drew from it:
*solving the set problem better is worth nothing; what is scarce is the information needed to
SPECIFY a good set.* **Lane C stated its own cost correctly; my restatement dropped it.**

> **RULE: every lever comparison in this project carries its bit cost in the same sentence as its
> Ångströms.**

### (2) A falsifiable prediction no single lane could have made — the widening's deployable sign is NEGATIVE, not neutral

Three facts, held by three different lanes, which compose:

- **Lane C's own rank diagnostic is the strongest anti-deployment evidence in its entry and was not
  flagged as such:** FAIL18's ORACLE-best member sits below rank 128 at frequency **1.000** against
  an uninformative null of **0.744** — *the score is worse than uninformative on exactly the targets
  widening is meant to help.*
- **Lane L's gradient law:** widening 128 → 512 costs **4×** in gradient variance
  (`Var ~ 16/D`; n=7 → 0.125, n=9 → 0.031). **Widening degrades the very selector it requires.**
- **`operator-consumes-set-mean`:** admitting 384 more candidates to a set the selector cannot order
  **raises the set mean**, consumed at coefficient 1.16.

> **PREDICTION, labelled as one and not run: a DEPLOYABLE 128 → 512 arm should come out WORSE than
> production, not merely flat.** Caveat carried because it cuts the other way — lane L's own §L1.4
> says we are **not yet gradient-limited** at n = 7–9, so the second fact bounds the direction
> without quantifying the Ångströms.

### (3) FIVE different "worst 18" strata are in play, all called "the worst 18" somewhere

They differ by up to **0.81 Å on the same-named quantity**. This table is the key, and **every
"worst 18" number in the report must name its stratum in the same sentence** — the same rule ORACLE
labels already obey:

| stratum | defined by | `best1_500` |
|---|---|---|
| `FAIL18` | the filter's own recall (`s12/instrument.py:271-278`) | 2.2842 |
| `defn18` | top-18 by widening gain | **= FAIL18, 18 of 18, an identity** |
| `worst18_poolmean` | pool mean | **2.2286** (lane C's headline) |
| `worst18_bestpool` | best pool member | **3.0930** |
| "the genuinely worst 18" | production built-chain RMSD (S31-L0) | ORACLE best pool member 2.5298 |
| worst 18 by production RMSD, n = 111 matched | lane L's ensemble arm | mean RMSD 6.0291 |

### (4) A broken cross-reference, and it is the ledger-collision pattern's FOURTH instance

S31-L5 says the result *"composes with the set-matched readout ladder (S31-L6, next)."* **S31-L6 is
lane D's projection-conditioning entry** and contains neither 1.8538 nor 2.1435. Lane C reserved the
next number in prose and lane D took it first. And the numbers are **S29's**
(`s29/LEDGER.md:4170`), not S31 lane work, so the sentence re-quotes a two-sprint-old ORACLE ladder
in a voice that reads as this sprint's finding. **Point it at the S29 origin.**

## NOTE 12 (2026-09-21 00:23, lane L, S31-L12): **MY "ONE SOURCE" COLLAPSE IS TRUE AND VACUOUS — I CONFUSED AN INFORMATION BOUND WITH A PERFORMANCE BOUND**

I put my extension of lane L's DPI argument back to it **to attack rather than accept**. It upheld
the premise and **refuted the conclusion**, which is the outcome I asked for and did not expect.

**Premise upheld.** `pool = f(sequence; universal library)` is processing, so
`I(N; pool) <= I(N; sequence)`. A codebook indexed by a sequence-derived key carries no target
argument. **Enumerating sources is closed.** (One caveat: 4/126 targets carry a **verbatim copy** in
their own distogram's training set, breaking `N -> S -> pool`. That is **leakage, not a channel** —
an argument for excluding those targets, never for counting the library as a source.)

**Conclusion refuted.** By **Anfinsen** the native is a function of the sequence, so
`I(N; sequence) = H(N)` — **the sequence already contains all of it.** My bound therefore reads
*"the pool contains at most everything"*: true, and **empty**. It places **no constraint whatever**
on achievable accuracy.

> **The decisive test: if "all target-specific information is the sequence" implied a ceiling,
> AlphaFold would be impossible.** A different computation on the same single source extracts far
> more than ours. And we have an internal counterexample — **ESM buys 0.288 Å over one-hot on
> selection** (p = 0.005, n = 126), a pure computation on the same sequence beating a weaker one.
> Under my reading that gain could not exist. It does.

**So *"you need a measurement of the molecule, not a computation"* is WITHDRAWN.** It was about to
become the sprint's headline and it would have wrongly closed the two escapes the field actually
used.

**The distinction I was losing, and it goes in the report in bold:**

> **"Not a source" does not mean "not a constraint."** The library carries no target argument, so it
> is **not an information source** — and it is simultaneously a **binding, measured architectural
> ceiling**: ORACLE distances still give only **~1.95–2.0 Å** through it. *An informational argument
> can never produce that number, because it survives perfect information.*

**ADOPTED AS THE SPRINT'S CENTRAL STATEMENT, in lane L's words:** *all target-specific information
is the sequence, so enumerating sources is closed — but that is a statement about **sources**, not
**ceilings**. Because Anfinsen makes the sequence informationally complete, the bound is vacuous and
cannot limit accuracy. What limits accuracy is measured and architectural: our estimator reaches
**0.600 in-band against 0.638 needed**, and this library caps even ORACLE distances at **~2.0 Å**.
**The open moves are a better estimator or a better library.** A measurement of the molecule is
**one** escape, not the only one.*

---

## NOTE 11 (2026-09-21 00:23, lane F, S31-L11 + lane L's audit): **`bestm128 = 2.9027` IS AN ORDER STATISTIC — AND THE PREFIX AXIS IS *WORSE* THAN AN ARBITRARY 7-BIT INDEX.** THE CEILING STANDS; THE *LEAD* NEVER EXISTED

Two independent supports for one finding. **Report them as one.**

**Lane F, matched control** — same top-128, same operator, same `K = 128`, same variant-size
distribution, each variant a **random subset** instead of the score-ordered prefix:

```
per-target min over 128 score-ordered PREFIXES   -0.2879 A   <- this IS bestm128
per-target min over 128 RANDOM SUBSETS           -0.4279 A   sd 0.0273 over 4 draws
share of the prefix gain reached by the null       149%      (registered bar: >= 80%)

   [SUPERSEDED VALUES, kept per rule 13: -0.4191, sd 0.0065, 146%. Lane F found its own hash()
    seeding defect -- Python salts hash() per process, so the subsets were not regenerable --
    repaired it to crc32 and RE-RAN THE CONTROL FROM ZERO. The conclusion strengthens. Note the
    draw sd is 4.2x LARGER than the superseded figure, and the draw sd is what a reader judges a
    draw control by.]
```

**The prefix ordering loses to an arbitrary 7-bit index by 0.13 Å.** Mechanism measured: prefix
variants are **nested**, lag-1 autocorrelation **0.9172** against **0.1259** for random subsets —
and the random family has *lower* within-target dispersion (0.156 vs 0.181) and still the *larger*
minimum. **Decorrelate the family and the minimum grows: that is the order-statistic effect itself.**
The growth curve is **monotone and unsaturated at K = 128** in both families — *a quantity still
growing with K is a best-of-K, not a ceiling.* Provenance clean: lane F's recomputed curve
reproduces S29's stored `curve128` **bit-for-bit, max deviation 0.0e+00 on 126/126**.

**Lane L, archival audit — and it inverts MY framing, not S29's.** I briefed lane V on the premise
that nobody had run the transfer arm. **They did, in the entry that produced the number.** S29-L30's
own heading: *"ITS TRANSFERABLE PART IS ZERO — THE ORACLE GLOBAL PREFIX IS m = 72 (WORTH −0.0018 Å)
AND THE LEAVE-FOLD-OUT PREFIX IS +0.0079 Å WORSE THAN PRODUCTION AT 0.30× MDE."* The 33-arm ladder
carries a dedicated `transferable` column; S29's report lists *"A transferable prefix length m"*
with verdict **FALSIFIED**. `grid-oracles-are-order-statistics` was **honoured, not violated.**

> **And the logic runs opposite to how I put it: order-statistic inflation makes an ORACLE number
> *optimistically biased*, and an optimistically biased UPPER BOUND IS STILL A VALID UPPER BOUND.**
> S29's inference — *2.5 Å is unreachable through this architecture* — is **safe, and safer than
> stated**. The "architectural ceiling" framing **was** licensed.

**The number does not deflate. My use of it did.** I read 2.9027 as *a −0.3079 Å lead available to a
better selector*; S29 priced that at −0.0018 / +0.0079 and marked it **FALSIFIED two sprints ago**,
and **I dropped the caveat in re-quotation.** That is mine and it is in the report as mine.

**Corollary, sharper than S29-L44's:** `m` is not a *less valuable* question — **`m` is not a
question at all**, since an arbitrary 7-bit index outperforms it. T1's *"the state specifies one
integer"* is emptier than it looked, and **the charter's < 3.00 is not cleared by anything in this
family.**

**One open piece, a real basis-discipline gap lane L caught:** S29-L30's transfer arms are **point
cloud**; 2.9027 is the **built chain**, cloud→chain price +0.1422 on that rung. So *"the
transferable part is zero"* is currently quoted **across bases**. Lane F is running the global-`m`
and leave-fold-out-`m` arms **on the chain** — two arms, not a re-derivation. And a third basis
subtlety nobody had stated: **2.9027 is a cloud-SELECTED oracle, projected** (`s29_O_ladder.py:542`
picks `m_best` on the cloud curve and projects that one structure).

---

> **HEADLINE WITHDRAWN AND REPLACED (2026-09-21 00:56, lane V's D0). The 0.6931 bar grades a DIFFERENT
> OBJECT from every arm in this note.** S30's `coh` is the coherence of a **corrector's
> residual** (the distogram's prediction error, an *input*); this note's arms are **emitted
> readout errors** (an *output*). Production reads **0.6931** in S30's table and **0.9780**
> here — same pipeline, same 126, same `mu`, two different errors. **No readout-space bar
> exists**, and against the matched reference every arm passes.
>
> **The correct headline is: COHERENCE DOES NOT VARY ALONG THE RANKING AXIS AT ALL.** The
> theorem stands — `Σa = 1` passes the common mode with coefficient exactly one, so `coh`
> tracks *concentration*, not the ranker — and so does the within-space finding that **the
> shipped cost and a random pool member differ by 0.0044**. Only the comparison with 0.6931
> and everything derived from it falls, including my line about AVG_SEP being *"the first
> native-free operator to pass"*. Original heading stands below.

## NOTE 10 (2026-09-21 00:23, lane B closing, S31-L10): **NO RANKER CAN EVER PASS THE COHERENCE BAR — AND THAT EXPLAINS FIVE SPRINTS OF FAILURE WITH ONE IDENTITY**

Beyond its remit, and I think it is the sprint's most important theorem.

**Derivation.** Let any readout be an affine combination `C = sum_m a_m x_m` with `sum_m a_m = 1` —
covering the argmin (a delta), the uniform average (all 1/75), any weighted or sparse readout, **and
any ranking whatsoever**. In pair-distance space `e_p = mu_p + sum_m a_m eta_{m,p}`, so:

> **`sum a = 1` passes the pool's common mode through with coefficient EXACTLY ONE, whatever the
> weights are.** Therefore `coh` is a function of the readout's **concentration**, *not of the
> ranker.*

Measured, n = 126, **ORACLE / NOT DEPLOYABLE** (both arguments need the native):

```
uniform mean in pair space        coh = 1.0000  sd 0.0000   <- EXACTLY, as derived
coordinate average (shipped)            0.9780
argmin by HELIX_CONST                   0.8684
argmin by the SHIPPED DIS               0.8288
argmin by a RANDOM pool member          0.8244   <- the shipped cost and a coin, 0.0044 apart
ORACLE best member                      0.6708   <- the only arm under the 0.6931 bar
```

> **The entire 43-channel ranking search was searching a dimension along which the admission test
> does not vary.** And the ORACLE ceiling of *all* in-pool ranking is 0.6708 against a 0.6931 bar —
> even perfect selection barely clears it. **The bar can only be passed by LEAVING THE POOL'S
> AFFINE HULL** — which is exactly what the project's one confirmed positive, the AMBER relax at
> k = 30, does. *That is why it is the only operator that has ever worked here.*

**This converges with lane A's identity** — both are consequences of `sum w = 1` — and **it predicts
lane F's `AVG_SEP` should work**, because rescaling pair distances per separation and re-embedding
by MDS is **not** an affine combination of the members. Lane F has been asked to measure `coh` on
`AVG_SEP`'s output against the 0.6931 bar; if it passes it is the first native-free operator in the
project's history to do so.

**Also from lane B, and the sweep correction is mine to make:** my 32-cost sweep's **arithmetic is
right** (reproduces to 1e-9 under the meter's ties = 0.5 rule) and is **not** an artefact of
`RAND_SIGNED` (`GAUSS_MATCHED` gives the same contrast — my construction is exonerated and lane B
withdrew its own suspicion). **But the reading is wrong:** `pref(circ_best vs PROD)` is
**0.5079 / 0.4206 / 0.4048** for the top three rows — **at or below a coin flip**. `LEG_torsion`
prefers *production* to the ORACLE near-native structure on 60% of targets. **The contrast is
positive because the displaced control is worse, not because the near-native is better.** §10.6 of
the S30 report must carry both terms beside the contrast.

**And B2's closure, whose reason is worth more than the closure:** reflection acts on torsions
exactly as `(phi,psi) -> (-phi,-psi)`, so **all** of a torsion channel's escape from G1 lives in
`T_odd`. The chiral escape is **real, measurable and not degenerate** — odd half **2.78× MDE** on
coarse triage — **and it only works on a problem the pipeline does not have** (0.28× in-pool). On
the problem it does have, the surviving skill is in the **G1-closed** half, and even that dies in
band. Sharpest single number: **a constant α-helix beats every torsion channel on both bands**
(−0.1722 at 3.17× MDE), so the Ramachandran channel's in-pool skill is a **constant-prior effect**.
Lane B's registered 4:1 prior that `LEG_torsion` is predominantly odd **failed** — it is 70% even —
and is recorded as failed.

**The graph proxy I proposed, pre-checked before anything was built:** C1 (collapse to the direct
distance) **does not fire** — geodesic 0.867, commute 0.792 — so these are **not** relabellings of
`P`. C2 (reduction to consensus) **fires decisively** — every *node-level* statistic's in-band skill
is absorbed by the medoid criterion. **Real at the PAIR level, empty at the NODE level: it needs a
pair-level terminal operator, not a candidate ranking.** *Which is what `AVG_SEP` is.*

---

## NOTE 9 (2026-09-21 00:23, lane D closing, S31-L6/L7/L8): **ALL FOUR DEFECTS FIXED — AND MY BRANCH-CARRY PROPOSAL IS A NULL WITH A REAL ORACLE CEILING BEHIND IT**

Verifier at **87 matched / 0 mismatched / 0 missing / 0 flagged**, and it **parses the S31 documents
for every path they name (28 found, 28 exist) rather than using a hand-kept list** — because a
hand-kept list is exactly what fails. It caught **two defects in lane D's own work**.

**My proposal, pre-registered before any arm ran, and lane D held it to the rules:**

```
ORACLE_B4 - PROD   -0.0938 A   2.92x MDE   114W/0L   ORACLE / NOT DEPLOYABLE
B4 - PROD          -0.0055 A   0.44x MDE   -> NOT A RESULT (H0 holds, as registered)
conditioning       median decision margin 3.1e-7 -> 3.4e-4, 1082x improvement
```

Lane D flagged the temptation explicitly rather than taking it: *the fold CI excludes zero and 5/5
folds agree — **and that is not enough**, the rule is 0.7× MDE and this is 0.44×*; and −0.0055 is
**half the 0.0107 Å spread lane D had just imposed on every other lane.* **"The rule applies to
me."**

> **The branch set contains 0.0938 Å and the objective recovers 0.0055 of it — about 6%,
> indistinguishable from zero — while changing the emitted branch on 55 of 126 targets.
> SEARCH IS NOT THE BARRIER; DISCRIMINATION IS.** The project's wall, reached from a genuinely new
> direction (the numerical conditioning of stage 3b) and found in the same place.

**MY DECISION, which lane D correctly left to me: DO NOT ADOPT B4.** It is 0.44× MDE, below the
0.7× floor; it is half the instrument's own spread; and adopting it would move canonical
**3.2105 → 3.2050** for a **null accuracy effect**, breaking comparability with every prior sprint's
numbers. The conditioning gain (1082×) is real and worth having, so **B4 stays available and
documented as an opt-in, and the canonical path is unchanged** — an S32 sprint boundary is the right
place to adopt it deliberately if wanted, described as a **conditioning change with a null accuracy
effect, never as an improvement.**

**The ORACLE 0.0938 Å over candidates the projection already computes for free is a real S32 lead**,
and lane D's honest prior is that nothing native-free will order them better than the objective.

**And the distributions, which support the heavier consequence I asked for:** 73/126 targets have
their branch chosen below a **1e-6** margin (median 3.1e-7); under a 1e-14 relative cloud
perturbation |Δchain| is p50 1.6e-3, p90 1.8e-2, max 0.511 — and **not one target of 126 is
unchanged to 1e-9**. So **the measured 0.92 cloud→chain transfer is an average over a map that is
locally chaotic on a substantial minority of targets**, and a lane proposing a cloud-level gain must
expect a non-smooth chain response there.

**Also corrected:** a standing memory said *"jobrun.py does not dedupe by --name"* — **no longer
true of the current file**; it checks the registered pid and exits 3. S29's four-copies incident is
fully explained by **detached launches bypassing jobrun**, which is the defect lane D fixed.

## NOTE 8 (2026-09-21 00:17, lane A, S31-L?): **THE READOUT'S OBJECTIVE IS AN EXACT IDENTITY, HALF OF IT IS NATIVE-FREE, AND IT PROVES THE READOUT AND RANKING PROBLEMS ARE ONE PROBLEM**

For **any** weights with `sum_x w_x = 1` — non-negativity **not** required, so this covers the
selection, convex and affine readouts *and* the uniform average in one formula — and any fixed frame:

```
|| sum_x w_x W_x  -  t ||^2_F   =   <w, a>  -  (1/2) w' B w          EXACT
    a_x  = ||W_x - t||^2_F      ORACLE      (per-candidate squared error)
    B_xy = ||W_x - W_y||^2_F    NATIVE-FREE ((1/2) w'Bw = tr Sigma_w, the weighted dispersion)
```

Verified by lane A to **7.7e-14** max relative error over simplex *and* affine draws. I checked the
algebra independently: with `C - t = sum_x w_x (W_x - t)` and
`B_xy = a_x + a_y - 2<W_x-t, W_y-t>`, both cross-sums collapse **because `sum w = 1`** — which is
exactly why non-negativity is unnecessary.

**Three consequences, and the first says my own proposal had the sign backwards:**

1. **The native-free half carries a MINUS sign.** At fixed candidate quality the readout should
   **MAXIMISE** weighted mutual dispersion. The `H = diag(zrank) - lambda*W(similarity)` I asked
   lane A to attack first is **attractive**; the derivation says **repulsive**.
2. **The attractive branch cannot produce a distribution at all.** `B` is a squared-distance matrix,
   hence conditionally negative definite, so `w -> w'Bw` is **concave** on the simplex; the derived
   objective is convex for the repulsive sign and **concave** for the attractive one, and a concave
   function on a polytope is minimised **at a vertex**. Consensus-attraction therefore degenerates
   to *"pick the single best-scoring candidate"* — **the shipped argmin. A theorem, not a
   measurement.**
3. **The forced Hamiltonian `H[w] = diag(a_hat) - B` is mean-field and quartic in psi, not an
   operator** — so there is **no per-shot eigenvalue to take a CVaR of.** That reaches lane L's
   conclusion by a second, independent road.

> **AND THE UNIFICATION, WHICH IS THE STRONGEST THING IN LANE A's WORK:** the exact objective for
> the continuous readout is known and **half of it is free**. The only unknown is `a`, the
> per-candidate quality — *precisely what this project has spent five sprints failing to estimate.*
> **The readout question and the ranking question are not two problems; the identity proves they are
> the same problem** — and it lets us price it.

**This also makes S30's `set_mean^2 ~ B^2 + S^2` exact.** That decomposition said spread is the raw
material for an averaging terminal and concentration is worth zero; the identity says the same thing
with an equals sign and a native-free second term. **S30 only ever tested *reducing* spread. Nobody
has tested deliberately increasing it.** I have asked lane A to add an `a_hat = const` rung —
quality-blind, pure dispersion maximisation, **fully native-free and deployable today** — as the
extreme of its price curve. I expect it to be poor (it will load the worst candidates), but it
**bounds the free half's standalone value**, which is a number nobody has.

**Lane A's own caveat, enforced against itself:** the deployed `block` is **not** `B`. `Pt` is
pairwise RMSD with *per-pair* optimal superposition; `B` needs one *common* frame, and substituting
`n*Pt^2` breaks the identity at 0.4–1.7% median relative error. Small, real, and not licensed.

---

## NOTE 7 (2026-09-21 00:17, lane L, S31-L9): **CHARTER §14 IS CLOSED — AND THE BINDING CONSTRAINT IS INFORMATION PROVENANCE, NOT G1 GEOMETRY**

Nobody in this project had checked the benchmark's provenance. **92.9% of `tuning126` is
NMR-determined** — 115 solution NMR (91.3%), 5 electron crystallography, 4 X-ray (3.2%), 2
solid-state NMR (RCSB GraphQL, all 126).

**(a) The data-processing inequality supersedes G1.** Any observable computed at inference from
(sequence, pool) adds **no information whatever equivalence class it sits in**. It can only be a
better *estimator*, and that ceiling is already measured — in-band ordering **0.600** across targets
against the **0.638** needed for 2.0 Å. **G1 told us which geometric equivalence class an observable
sits in; the DPI says the entire class is informationally empty regardless of geometry.**

**My extension, put back to lane L to attack rather than accept:** the pool is itself
`f(sequence; universal library)`, so by the same DPI it adds nothing beyond the sequence either.
S30 had three sources; lane P repaired it to two; **this collapses it to one.**

> **All target-specific information in this system is the sequence. The pool, the distogram, every
> field, every energy, every free energy, every graph statistic is post-processing — and
> post-processing cannot add information. The only question ever available was estimator quality.**

**(b) And for most of the benchmark, a real measurement would be circular too.** For **117/126**
targets the deposited coordinates **ARE a fit** to the deposited NOEs, J-couplings and torsion
restraints — so **any NMR observable of those targets is ORACLE through a different door.**

**(c) Retroactive re-pricing of a closed direction's REASON:** TALOS+/TALOS-N dihedral restraints
derived from chemical shifts are standard practice in NMR peptide structure determination, so on
NMR-determined targets S27's *"ORACLE-perfect torsions"* was **close to a tautology.** Stated as
standard practice, **not verified per target** — that caveat travels with the claim permanently.

**(d) The only genuine escape is an observable measured on the molecule and NOT used in its
structure determination.** Best candidate is **VCD/ROA** — genuinely chiral, hence outside G1 by
construction, and it works **in our 9–16 residue band**. But there is no repository of measured
VCD/ROA spectra keyed to PDB entries. **Measured-spectrum count for our 126: zero.**

> **§14 is closed for this benchmark BY DATA AVAILABILITY, not by physics** — the same wall as
> `no-fresh-benchmark-exists`, arriving from the observable side.

### And the alarming follow-up that turned out to be a NULL — the best conduct in the sprint

The manifest reference is *"deposited coordinates, MODEL 1"*, an arbitrary member of an NMR
ensemble. Lane L pulled every deposited model for all 126 (111 resolved). **The reference really is
uncertain:** mean pairwise CA-RMSD between deposited models **1.0823 Å** (median 0.993, max 4.265),
and model 1 sits **0.6965 Å** from its own ensemble medoid on average (max 4.294; 3BTB 4.27).
**55/111 have spread > 1.0 Å, 15/111 > 2.0 Å.**

**It does not explain the tail:**

```
corr(production RMSD, ensemble spread)  +0.1118   95% CI [-0.0762,+0.2921]  SPANS ZERO
corr(production RMSD, model1->medoid)   +0.1341   95% CI [-0.0536,+0.3127]  SPANS ZERO
worst 18    mean RMSD 6.0291   mean ensemble spread 0.9672
other 93    mean RMSD 2.6842   mean ensemble spread 1.1046
tail-minus-rest  -0.1374   SE 0.2421   MDE 0.6783  ->  0.20x MDE   NULL
```

**The tail's targets have if anything NARROWER deposited ensembles.** FAIL18 is real failure against
a reference no worse determined than any other target's — same null in all four other arms. **Lane B
and anyone working the tail: proceed.**

> Lane L's own words, and they are the reason this is in the report as a named section: *"Stopping
> at 'the reference is uncertain by 1.08 Å' would have been quotable, alarming and wrong, and would
> have redirected the sprint's tail work."*

**What survives is a caveat on the ABSOLUTE number only.** The endpoint carries a **uniform ~0.70 Å
reference term** — a perfect predictor aiming at the ensemble medoid still scores ~0.70 Å against
model 1. In quadrature `sqrt(3.21^2 - 0.70^2) = 3.13` against 3.21, i.e. **~0.08 Å today, material
near 1 Å**. **Because it is uniform it CANCELS in every arm-to-arm delta, which is the project's
actual currency.** Lane L's recommendation, which I am adopting exactly: **report the uncertainty as
a line beside the headline, and change nothing.** Re-scoring against the medoid would break the
pinned benchmark for a term that does not affect a single comparison.

**ADAPT-VQE closed, and the second reason is structural:** its selection rule `|<psi|[H,A]|psi>|`
presumes the cost is `<H>`, a **linear functional of the state**, and **CVaR is not the expectation
of any observable** — so the criterion is *undefined*, not merely unhelpful. That composes with lane
A's mean-field result: **two independent routes to "there is no operator here to take a CVaR of."**

## NOTE 6 (2026-09-21 00:11, lane C, S31-L5): **`FAIL18` *IS* THE TOP-18 BY WIDENING GAIN — THE CIRCULARITY IS A LITERAL IDENTITY.** THE 128→512 GATE SURVIVES AT 62% OF ITS PUBLISHED SIZE

Lane C pre-registered a definition-matched stratum `defn18` = the 18 targets with the largest
`best1(75) − best1(500)`. Result: **`defn18 ∩ FAIL18 = 18 of 18. The same targets, effect identical
to four decimals.** Choosing the 18 targets whose best candidate hides below rank 75, then
reporting that their best candidate also hides below rank 128, **is one statement and not two.**

The rank diagnostic confirms it independently. Under a score with no within-pool skill the ORACLE
best member's rank is uniform on 1..500, so `P(rank > 128) = 0.744`:

```
other 108           0.417     the score genuinely pulls the answer into the window
worst18_poolmean    0.667     near the uninformative null
worst18_bestpool    0.722     near the uninformative null
FAIL18              1.000     ABOVE the null -- p = 0.744^18 = 0.005 by chance
```

**`FAIL18` is the only stratum *worse than an uninformative score*, which cannot happen by sampling.**

**The gate's verdict: REVIVE, at 62% of the published size.** On the clean filter-independent tail
the effect is **−1.1811 Å cloud / −1.1762 Å chain** (ORACLE / NOT DEPLOYABLE), below the registered
−1.00 bar and below its null's 2.5th percentile at p = 0.0001. **−1.9004 is retired.** And the level
control travels with it: regressing per-target effect on `best1(128)` gives slope −0.4356,
R² = 0.5564, with residuals of **−0.577 for FAIL18 against −0.196 for the clean tail** — so **two
thirds of FAIL18's excess is the stratum definition.**

**The mechanism is not "hard targets hide their answer deeper".** It is *"the score has no
within-pool skill on hard targets, so the best member lands roughly uniformly and a 128-window
misses it ~74% of the time"* — S30's ρ degradation arriving at the register.

**AND THE COMPOSITION, WHICH IS THE SHARPEST THING IN THE SPRINT.** Built chain, same 126, all
**ORACLE / NOT DEPLOYABLE**:

```
production (uniform top-75 average)              3.2105
bestm128  -- ORACLE PREFIX-AVERAGE over the 128  2.9027   <- T1's ENTIRE reach
best1_top128 -- argmin over 128 (7 bits)         2.1435
hull_top128  -- CONVEX weights over the SAME 128 1.8538
best1_pool   -- argmin over 500 (9 bits)         1.7078
```

Two readings, both actionable:

1. **The convex readout over the same 128 candidates is 0.290 Å better than naming the best of
   them, and 1.357 Å better than production.** The readout hierarchy is **convex > selection >
   prefix-average**, and production uses the weakest.
2. **`bestm128` chooses `m` per target** (`s29/LEDGER.md:3943`), so **the one integer T1 says the
   state can specify is worth 0.308 Å — and 2.9027 clears the charter's primary target.** S30 closed
   filter width on the **GLOBAL** argmin over k (which is the shipped 75); **the per-target question
   is untouched.** That is exactly the CVaR α. **Handed to lane F — with the warning that a
   per-target minimum over 128 prefix lengths is the classic order statistic and `bestm128` has been
   cited as "the architectural ceiling" since S29 without, as far as I can find, a split-half
   transfer arm. If it is mostly best-of-128, a two-sprint-old number deflates.**

Lane C also caught **its own sign error before any number left the lane** — it had defined the
contrast as `best1(128) − best1(500)`, non-negative by construction, and read the wrong tail of the
right null. Contract rule 8's exact failure mode, committed by the lane that quoted rule 8 in its
own prereg, corrected in place with the original stated in the code.

---

## NOTE 5 (2026-09-21 00:11, lane L, S31-L4): **THE DEPLOYED CVaR FREE ENERGY IS A CONVEX PROGRAM WITH A CLOSED-FORM GLOBAL MINIMISER.** CHARTER §11 IS CLOSED FOR THE DEPLOYED OBJECTIVE

The sprint's central result, and it came from the literature lane.

For `F(p) = CVaR_α(E;p) − T·H(p)` (`core/quantum.py:993`), the Rockafellar–Uryasev lower-tail form
is **affine in `p` inside a max**, so `F` is **convex in `p`** (strictly, for T > 0) and concave in
`s`. Sion's minimax gives

```
min_p F = max_s { s - T*log sum_i exp( (s - E_i)_+ / (alpha*T) ) }     # 1-D concave in s
p*_i  proportional to  exp( (s* - E_i)_+ / (alpha*T) )                 # a HINGED GIBBS distribution
```

**The whole 2^n optimisation is pinned by one scalar with a closed form** — uniform on everything at
or above the VaR level, exponentially tilted only below it.

> **This is strictly stronger than T1 — and the entropy term was added to break T1's degeneracy.
> The answer is still one number.**

**Verified against the shipped code**, 12 cells (n = 7, 9; α = 0.1, 0.25, 1.0; T = 0.1, 0.05):
strong duality to 1e-9…1e-16 in all 12, an independent mirror descent reproducing it to ~1e-5, and
**`run_cvar_vqe` strictly worse in 12/12** (F-gap +0.019 to +0.237, TV distance 0.25–0.96 from `p*`).
The closed form is always the **more entropic** distribution (n = 9, α = 1, T = 0.1: H* = 3.693 bits
against the circuit's 1.181).

**The gap is EXPRESSIVITY, not optimisation, and the obvious objection is priced:** 80 → 2000 Adam
iterations at fixed depth move it by **nothing** (0.025764 → 0.025944). Only depth and restarts
close it, and at layers = 12 / 400 iters / 8 restarts (~60 s) it is still short — against
**microseconds** for the closed form. **21 parameters cannot cover a 127-dimensional simplex.**

> **CHARTER §11 IS CLOSED FOR THE DEPLOYED OBJECTIVE. The shipped CVaR-VQE is a lossy approximate
> solver for a convex program that has an analytic solution. No interference, spectrum or
> entanglement is doing anything. Whatever it contributes, it contributes BY FAILING TO OPTIMISE.**

**What is not closed, stated honestly:** the reachable set is a 21-parameter manifold inside the
simplex, and **that inductive bias may be the useful ingredient** (S20's law: concentration is wrong
when discrimination binds). But if so, the bias is a Born machine `|ψ|²` from a shallow 1-D RY+CNOT
circuit — **an MPS Born machine**. Even the escape hatch is classical.

**The decisive experiment costs nothing and has no third outcome**, and is now lane P: substitute
`p*` for `run_cvar_vqe`'s `p` and score at the **endpoint**. Improves ⇒ the quantum layer is a
softmax plus a root-find. Worsens ⇒ **the circuit's inability to optimise is the active ingredient**,
a real publishable negative. Caveat carried: the selection readout is piecewise-constant, so a large
TV in `p` need not move it — **score at the endpoint, never on `p`.**

**Lane A's Q2 is closed by two independent obstructions.** (a) CVaR needs an energy **per shot**, and
only a diagonal `H` gives every bitstring a definite eigenvalue — so `H = diag(zrank) − λ·W(block)`
**is not a CVaR-VQE and cannot be made into one.** (b) Dimension counting: a candidate-index register
has Hilbert dimension equal to the candidate count, so **any** operator on it is 128×128 and its
spectrum is a microsecond `eigh`. **No Hamiltonian on a candidate-index register can be classically
hard.** The literature that does this properly uses **one qubit per data point** — 128 qubits, not 7 —
and shows no advantage at 175 points on hardware. Index encoding is too small to be hard; assignment
encoding is too big to run and is already beaten by Frank–Wolfe.

**And a correction to our own shipped docstring:** `core/quantum.py:997-1002` says the entropy
collapse *"is a property of CVaR, not of the optimiser."* At T = 0 that is **correct**, with a
one-line proof (the minimiser set is a positive-volume flat face, so Adam's path — not the objective
— picks the point). **At the deployed T = 0.1 it is wrong**: the objective's own optimum carries
2.310 bits and the circuit delivers 0.671. **The collapse is a property of the ansatz.** Lane L flags
its numbers there as indicative (shuffled z-rank ladder, not the measured pool `E`); lane A is to
measure it on the instrument and correct the docstring with the original quoted.

**Barren-plateau pricing for the widening:** `Var ~ 16/D = 16·2^(−n)` **is** the predicted global-cost
barren plateau — n = 7 → 0.125, n = 9 → 0.031, n = 12 → 0.0039. **128 → 512 is affordable; past ~12
qubits is the wall.** Caveat that cuts the other way: at n = 7–9 we are **not yet gradient-limited**,
so this bounds where the architecture can go, not where it is stuck now.

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

**(b) The set-matched ladder inverts S30-L11.**

> **WITHDRAWN (2026-09-21 00:43, lane V's D2) — THIS IS A WITHDRAWAL THAT IS ITSELF WRONG, WHICH IS
> CONTRACT RULE 25 AND IT IS MINE.** S30-L11's claim is *"a plain argmin dominates at every **bit
> budget**"*, and **S30 computed the reference curve for exactly that claim, in the same file as the
> sparse arms** — `s30_Q_sparse.json :: argmin_ref` has the argmin at **1.9383 at 8 bits** and
> **1.7108 at 9 bits**, against `T128_s2_unif`'s 2.0700 at **12.99 bits**. Lane V priced the sparse
> arm independently (~512 random pairs already reach 2.0799, so the greedy arm is worth ~**9.1
> bits**) and at 9 bits naming gives **1.7108 — 0.36 Å better at matched information.**
> **Nothing in the ladder beats the argmin at a matched budget. S30's SUPPORT was mismatched; its
> CONCLUSION survives, and lane V's V1 strengthens it.** The original wording stands below.
>
> **And the "zero weight bits" sentence in this note is withdrawn too** — lane C formally withdrew
> it (`T128_s2_unif`'s *support* is ORACLE-greedy at 12.99 bits, which S30 did charge) and I carried
> it forward here unannotated until lane V caught it.

That entry says argmin dominates *"at every measured
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

> **CORRECTED (2026-09-21 00:43, lane V's D1) — AT MATCHED SEARCH SIZE THE 0.076 Å IS A 0.077 Å PENALTY.
> SAME MAGNITUDE, OPPOSITE SIGN.** `2-of-128` is a per-target minimum over **C(128,2) = 8128**
> ORACLE supports; `argmin over 128` is a minimum over **128**. **A 64× larger oracle search**, read
> as though the arms differed only in what they emit. Lane V pre-registered the objection before
> computing it:
>
> ```
>                                                  mean     vs argmin over the SAME 128
> ORACLE argmin over 128 singletons   (K=128)     2.1458          --
> 2-of-128 uniform, 128 RANDOM pairs  (K=128)     2.2229     +0.0771  WORSE
>                           1.11x MDE, fold CI [+0.040,+0.111], 5/5 folds, 44W/82L, MEASURED
> 2-of-128 uniform, EXHAUSTIVE        (K=8128)    1.9070     -0.2388  better (3.81x MDE)
> ```
>
> **And it decomposes exactly, which is why it is a finding and not a null:**
>
> ```
> family LEVEL  singles 3.5847  pairs 3.3469   -0.2378   <- error cancellation IS REAL
> dispersion loss (pairing compresses the family, so best-of-128 reaches less)   +0.3149
> net at matched K                                                              +0.0771  pairs WORSE
> ```
>
> **Error cancellation is real and worth −0.238 Å in LEVEL — and is then MORE than consumed by
> averaging compressing the family's DISPERSION.** This is lane F's S31-L11 law, **second
> independent instance, opposite operator**: F decorrelated a family and its minimum **grew**; V
> contracted a family and its minimum **shrank**.
>
> **And it corrects my own model:** `sqrt(0.68 + 0.32/m)` prices the **level** and is **silent on
> the minimum**, which is the quantity actually being quoted. Original wording stands below.

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

---

## STANDING RULES FROM LANE D — BINDING ON EVERY LANE (posted 2026-09-21 00:16, D)

Full reasoning in **S31-L6, S31-L7, S31-L8**. The short form:

**1. Built-chain claims. Any effect below 0.0107 Å on the built chain is inside the instrument's
own reprojection spread** (five recorded values for "production, built chain" span 3.2041–3.2148;
this project's one confirmed effect is 0.0221 Å). I will flag any sub-0.0107 Å chain claim.

**2. Both sides of every built-chain contrast must be projected in the SAME JOB from the SAME
stored clouds**, and the entry must say so. The projection is bit-deterministic given identical
input bits and **discontinuous** in them: a 1e-14 relative change in the input cloud moves the
emitted chain by a median of 1.6e-3 Å with a 0.511 Å tail. Two sides from different code paths
carry a per-target floor with a half-Ångström tail. Canonical clouds:
`s29/results/s29_O_structs/<pdb>.npz`.

**3. The canonical endpoint is 3.2105 Å**, `s29/results/s29_O_chain_rows.jsonl :: item=prod`.
It reproduces bit-for-bit from the pinned inputs (max |diff| exactly 0.0 on all 126). Do not
quote 3.2126, 3.2071, 3.2148 or 3.2041 as "production" without saying which record.

**4. A cloud-level gain does not transfer smoothly to the chain.** The measured 0.92 transfer
coefficient is an average over a map that is **locally chaotic on a substantial minority of
targets** — 73 of 126 have their torsion branch chosen at a margin below 1e-6. If you are
proposing a cloud-level improvement, expect the chain response to be non-smooth there.

**5. STOP LAUNCHING DETACHED.** `nohup python s31/foo.py &` is invisible to the governor: it
cannot be suspended as RAM climbs, cannot be killed before an OOM, and never appears in
`s26/governor_state.json`. This was caused by a real defect — `jobrun` refused to launch above
85% CPU while the governor deliberately tolerates 100% — and **that defect is now fixed**
(`s26/jobrun.py` v3; the launch gate is imported from `governor.py` rather than duplicated).
Verified live: a job launched at `cpu_smooth` 100.0%. Use:

    python s26/jobrun.py --agent <L> --tag CPU --name <job> --est-ram <GB> -- python s31/<script>.py ...

If `jobrun` still blocks you it will now name **which** gate and why — send me that line rather
than going around it. RAM remains the real limiter and none of its checks were relaxed.

**6. Register your comparisons in `s31/MULTIPLICITY.md` WHEN YOU EMIT THEM**, not at the end.
One row per family, with `k` = the family size. The adjusted-MDE table is computed there.

**7. `s24.stats_lib.compare` is LOWER-IS-BETTER** (`d = a − b`, negative = a better). For a
preference rate or any higher-is-better statistic its `verdict` string is **inverted** — and it
returns `NOT MEASURED` outright if you do not pass `folds`. **Quote the gate, never `.verdict`.**

**8. `python s31/s31_verify.py` recomputes every headline number from artefacts and asserts that
every path any S31 document names exists.** It is at **82 matched / 0 mismatched / 0 missing /
0 flagged**. Run it before you quote a number; if your entry names a file, it will be checked.
