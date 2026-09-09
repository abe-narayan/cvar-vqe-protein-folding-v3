# SPRINT 24 — WORKSTREAM C: THE FROM-SCRATCH GENERATOR

**Verdict: the from-scratch `p(φ,ψ | sequence, structural prior)` path is CLOSED. No network was
trained, and the reason is a diagnosis rather than a refusal: the conditioning channel is empty, the
chain correlation a sequence model would add is not present in the data, and the selector — not the
source — installs the common mode. The lane's compute goes elsewhere.**

Pre-registration: `s24/PREREG_C.md` (six forks, falsifiers, stopping rule; written before the
confirmatory run). **The lane falsifier I pre-registered fired.**

| artefact | status | n |
|---|---|---|
| `c_ladder.py` → `results/c_ladder.json` | CONFIRMATORY, complete, stamped | 126 |
| `c_final.py` → `results/c_final.json` | CONFIRMATORY, complete, stamped | 126 |
| `c_bias.py` → `results/c_bias.json` | CONFIRMATORY, complete, stamped | 126 |
| `c_markov.py` → `results/c_markov.json` | CONFIRMATORY, complete | 126 |
| `c_probe{,2,3,4}.py` | **EXPLORATORY**, superseded by the above | 26 |

All artefacts stamped via `stats_lib.save_atomic(..., module_file=...)` — module name, source
sha256, git commit, dirty flag, launch time, and `complete` set on the **full key set**.

**Scope.** This verdict covers the from-scratch, sequence-conditioned path only. Lane B's residual
architecture is a different object and nothing here is evidence against it.

---

## C1 — THE HEADLINE. THE COVERAGE CLAIM'S PREMISE IS TRUE AND ITS CONCLUSION IS FALSE.

Route (a) was the strongest available case for this lane: *the corpus is missing conformations the
score would rank highly if it ever saw them.* It needs no decorrelation and no trained model. I
measured it rather than betting on it. **n=126, matched count 500 generated against the shipped
K=500 pool:**

| arm | frac scoring better than pool's 75th | share of merged top-75 taken from generated | ORACLE best of generated |
|---|---|---|---|
| `T0_helix` | 0.031 | 0.097 | 2.942 |
| `T1_blind` | 0.070 | 0.235 | 2.166 |
| `T2_restype` | 0.078 | 0.256 | 2.118 |
| `T3_pool` | **0.249** | **0.523** | 2.020 |

*(pool ORACLE best **1.7108**; incumbent 3.0483)*

**The first half of route (a) holds, and holds strongly.** A quarter of `T3_pool`'s samples score
better than the pool's 75th-best retrieved window, and when a merged 1000 is scored by one
functional the score takes **52% of its final 75 from the generated half.** The corpus genuinely is
missing conformations the shipped functional ranks highly, and a trivial untrained sampler supplies
them in quantity.

**The second half fails.** Those chains are structurally *worse* than what they displace: generated
ORACLE-best 2.02–2.17 Å against the pool's **1.71 Å** at matched count, with ~2 samples in 500
beating the pool's best.

> **The premise of the coverage claim is true and its conclusion is false, because the score's
> preference is not aligned with structural quality.** This is
> `objective-does-not-rank-the-native` — the native sits at the 36.8th percentile of the distance
> objective — surfacing as the precise, quantified reason candidate generation cannot pay.

**The endpoint agrees.** Every matched-size mixture and every union is null or worse:

    T2_restype  m60_15  -0.0037  SE 0.0120  MDE 0.0336  fold[-0.012,+0.003]  75W/51L  NULL
                m50_25  +0.0038  SE 0.0213  MDE 0.0597  fold[-0.021,+0.038]  73W/53L  NULL
                UNION   +0.0157  SE 0.0201  MDE 0.0562  fold[-0.036,+0.064]  54W/72L  NULL
    T3_pool     UNION   +0.0509  SE 0.0239  MDE 0.0671  fold[+0.007,+0.094]  45W/81L  NULL   genfrac 0.523

`T3_pool` sits **on** the interpolation line between the mixture endpoints (departures −0.009 to
−0.011), meaning its bias is fully parallel to the incumbent's. `T1`/`T2` bow only −0.035 to −0.084,
against the blind library's −0.226 in L2(b). **A live generator supplying half the final set
reproduces L2(d)'s precise null.**

---

## C2 — THE φ RESULT, CONFRONTED. WHAT I EXPECTED TO DIFFER, AND WHAT ACTUALLY DID.

From `PREREG_C.md` §0, written before the run. `phi-carries-no-sequence-signal`: full sequence
context predicts φ at **36.1° against 36.4° sequence-blind** — 0.8%; the entire channel is **10.4° of
ψ**; direct sequence→torsion→build emits **4.151 Å**. I pre-registered three ways my object could
differ, and pre-committed to expecting the lane to fail anyway:

1. **ψ's error is genuinely bimodal**, so a mixture expresses what a point estimate cannot. **Real.**
2. **The readout averages 75 samples** rather than building one chain. **Real, and the larger of the
   two.** Direct build was 4.151 Å; the same information class through a top-75 average reaches
   **3.21–3.24 Å**.
3. **The distogram enters generation, not only scoring.** **This is the one that failed**, and for a
   reason I did not anticipate: conditioning on the distogram does not add information, it adds *the
   selector's own error* (C4).

> The two mechanisms I predicted would help did help, by roughly **0.9 Å**, and the lane still does
> not move the endpoint by a measurable amount. **The φ result was not overturned; it was routed
> around, and the thing on the other side of it was the same wall.**

---

## C3 — REPRESENTATION AND MODEL CLASS: WHY NO NETWORK WAS TRAINED.

Directive §38 asks which *component* failed. Here is the answer, with the evidence.

### The representation comparison (readout, functional, count, basis all fixed; n=126)

| representation | conditioning | standalone, scored top-75, POINT CLOUD | cos vs incumbent |
|---|---|---|---|
| REBUILD75 — the incumbent's own 75, rebuilt from their own torsions | — CALIBRATION | **3.0524** | 0.9912 |
| constant α-helix + 15° jitter | none — plausible zero-information control | 3.7892 | 0.7508 |
| resampled (φ,ψ) **pairs**, sequence-blind | corpus | 3.2435 | 0.8043 |
| per-residue-type Ramachandran | sequence, at its measured ceiling | 3.2065 | 0.8209 |
| **per-residue 2-component von Mises on (sin φ, cos φ, sin ψ, cos ψ)** | retrieved pool | **3.1752** | 0.9055 |
| the same von Mises + **first-order Markov basin chain** | + chain correlation | 3.1556 | 0.8994 |

*(incumbent 3.0483 point cloud; the built-chain figure is 3.2126, reported separately and never
compared across the basis)*

**Three findings decide the model class.**

1. **The torsion manifold is free.** Rebuilding the incumbent's own top-75 through ideal geometry
   and re-running the identical readout costs **+0.0041 Å**. Whatever closed this lane, it was not
   the representation — which matters for anyone working in the same space.
2. **Sequence conditioning bought nothing.** `T2_restype` and `T1_blind` differ by 0.037 Å. That is
   the φ memory reproduced in RMSD space rather than in degrees.
3. **Chain correlation — the one thing an autoregressive model, an HMM or s19's entangled latent
   adds over a product distribution — is not there to capture.** `c_markov.py` fits a first-order
   Markov chain to the pool's own basin-label sequence and compares it against a control with the
   transition matrix replaced by the product of its own marginals, so **the marginals are held
   exactly fixed and the only difference is correlation**:

       T4_markov - T5_shuffle   -0.0095  SE 0.0130  MDE 0.0363  fold[-0.026,+0.003]
                                63W/63L  worst +0.49  NULL  (effect/MDE 0.26x)

   And the mechanism: **adjacent-residue basin mutual information in the retrieved pool is 0.0499
   nats against a maximum of 0.6931 — 7.2%.** Neighbouring residues' conformational basins are
   nearly independent at this length. There is almost no chain correlation for a sequence model to
   model.

### The reference a network would have to beat

**`T3_pool`: three parameters per residue per component, zero training, 3.1752 Å standalone.** s14
measured this exact channel at **φ 33.6° / ψ 59.2°** — better than this project's own trained
leave-fold-out sequence predictor at 36.1° / 62.4°.

> **A cVAE, MDN, autoregressive torsion model or diffusion model would have to beat a free
> construction that already beats the best trained torsion predictor here, using a conditioning
> channel measured at approximately zero, to add a chain correlation that is measured at 7% of its
> maximum and worth 0.26× its own MDE.** The component that failed is the **conditioning channel**,
> not the model class, not the representation, not the optimiser, and not sample diversity. Per
> directive §12 the correct model is the smallest one showing genuine multimodality, and that model
> needs no gradient steps.

---

## C4 — THE MECHANISM: SELECTION INSTALLS THE COMMON MODE. THIRD INDEPENDENT ROUTE.

L3 predicted that quality obtained *via* the distogram score causes alignment. I had built this probe
before L3 arrived; it confirms L3 **with a zero-information control that neither L3 nor
`referent.py` has.** Cosines between **signed pair-distance error vectors**, n=126, MAP estimator:

| arm | scored top-75 | uniform 75 (no selection) | what selection adds |
|---|---|---|---|
| INCUMBENT | 0.6259 | — | — |
| `T0_helix` (knows nothing) | 0.5930 | 0.4558 | **+0.137** |
| `T1_blind` | 0.6279 | 0.4427 | **+0.185** |
| `T2_restype` | 0.6403 | 0.4467 | **+0.194** |
| `T3_pool` | 0.7024 | 0.5556 | **+0.147** |

> **Applying the selector raises alignment with the distogram's own prediction error by +0.166
> (MAP) / +0.177 (posterior mean) on average, in every arm, including one carrying zero information
> about the target.** A bias that survives replacing the entire source is not a property of the
> source.

*(Correction to my own exploratory report: I first quoted **+0.20** from the n=26 probe. At n=126 it
is **+0.166**. The confirmatory figure supersedes it.)*

This is the third independent route to the mechanism, alongside the coordinator's `referent.py`
(β = 0.520 against a 0.352 placebo floor) and L3's quality-matched retrieval-free source.

### RECIPE PINNING FOR LANE D — one of the five choices differed, and here are both numbers

| | Lane D (`referent.py::_stats`) | Lane C | |
|---|---|---|---|
| (a) `Dhat` | MAP, `grid[argmin(risk,1)]` | **`dg["expected"]`, the posterior mean** | **DIFFERED** |
| (b) `min_sep` | 2 | 2 (`I.pair_index` default) | match |
| (c) `Dc` | uniform coordinate average, own medoid frame, point cloud | same | match |
| (d) averaging | per target, then averaged | same | match |
| (e) centring | neither centred | same | match |

`c_final.py` recomputes everything under **both** estimators on all 126 so D can diff per target
without guessing. **The estimator choice is worth −0.007 to −0.017 in the cosine and does not change
any conclusion**; the headline moves from +0.177 to +0.166. Separately, my `_bias`/`_cos` were
verified **byte-identical** to `residlib.bias/cos` (max abs difference 0.0), which are in turn
identical to `qmatch._bias` — so every cosine this lane quotes is the same object as L2's and L3's
by construction.

---

## C5 — THE SIGNED-BIAS TENSION, RESOLVED AT n=126, AT 5× MDE.

The coordinator flagged an apparent contradiction: the distogram predicts distances **too long**
while the emitted cloud is **too short**, yet β says half the prior's error transfers. `c_bias.py`,
n=126, complete:

    the DISTOGRAM's own prediction              +0.4062  SE 0.1807  fold[+0.128,+0.658]
    the whole legal universe (real geometry)    +0.4419  SE 0.1887  fold[+0.089,+0.712]
    the 75 SELECTED members, before averaging   -0.0627  SE 0.1745  fold[-0.383,+0.297]
    the EMITTED coordinate average              -0.5852  SE 0.1729  fold[-0.943,-0.174]
      -> what AVERAGING alone contributes       -0.5226  SE 0.0373  MDE 0.1046  fold[-0.586,-0.467]

    projection onto the distogram's error direction   members +0.5665   average +0.5413
    cosine with the distogram's error direction       members +0.7079   average +0.6351

**There is no tension; two separate things were being added.**

1. **The distogram's signed offset is not a distogram defect.** The whole legal universe of *real*
   protein windows carries the same signed offset (+0.442 vs +0.406). The signed mean measures a
   corpus/native scale mismatch, not prediction error — a `shared-referent-floor` effect. Measure
   the floor before attributing the quantity.
2. **The contraction is the averaging operator**, at **−0.5226 Å, 5.0× its own MDE**, fold CI
   excluding zero — reproducing the known 25.8% backbone contraction.
3. **The two are orthogonal.** Averaging leaves the inherited-error component essentially untouched
   (0.5665 → 0.5413) while moving the scale mode by half an Ångström. **β prices the error
   *direction*; the signed mean prices the *scale* mode.**

> **This bounds what a readout change could buy.** The selected members already carry the inherited
> error at cosine 0.708 *before any averaging*. The common mode is installed by **selection**, so
> fixing the readout cannot remove it — consistent with s23 L9's finding that every within-pool
> operator works on the 32%.

---

## C6 — THE SPEC WAS RETIRED ON THIS LANE'S EVIDENCE.

The L2 spec's two halves were measured on different readouts: L2's 0.647 came from an **unscored**
library draw; the RMSD half and the deployed operator use the **scored** top-75. Both, n=126:

| arm | SCORED: RMSD / cos / within / ratio | UNIFORM: RMSD / cos / within / ratio |
|---|---|---|
| `T0_helix` | 3.789 / 0.751 / 0.988 / 0.760 | 3.997 / **0.706** / 0.978 / 0.722 |
| `T1_blind` | 3.244 / 0.804 / 0.949 / 0.848 | 3.806 / 0.577 / 0.933 / 0.618 |
| `T2_restype` | 3.207 / 0.821 / 0.952 / 0.863 | **3.793 / 0.577** / 0.938 / 0.616 |
| `T3_pool` | 3.175 / 0.906 / 0.987 / 0.917 | 3.465 / 0.900 / 0.972 / 0.925 |

Under the scored reading every arm fails the cosine half; under the unscored reading `T2_restype`
**passes both** (3.793 ≤ 3.9, 0.577 ≤ 0.65) — and its mixtures are still null. **The coordinator has
ruled that the cosine is defined against the scored readout and has withdrawn the spec as a
go/no-go on this lane's evidence.** The decisive number is `T0_helix`: **a zero-information constant
α-helix reaches cos 0.706 unscored, landing on L2's library ratio.** Decorrelation is free, so
passing the cosine bar carries no information about whether a source is useful — the same shape as
L1, where every source with a "better" `fcommon` was worse. The cosine survives as a diagnostic and
as the instrument that established C4; it is no longer a bar.

---

## C7 — AUDITS, VIA LANE B'S HARNESS. AND LANE B'S WARNING FIRED ON MY BEST ARM.

`residlib.collapse_audit` / `validity_audit` taken verbatim rather than rebuilt. Emitted sets are
**replayed from `c_ladder.py`'s exact seed stream**, so these audit the sets whose RMSD is quoted,
not a fresh draw. Collapse audit on a stated 250-sample subsample (it is O(B²)).

**Mode collapse — none.** Every arm: 250/250 unique, `dup_frac` 0.0000, ESS 250.0, max mode
occupancy 0.0040. No RMSD in this document comes from a sampler emitting copies of one structure.
Torsion entropy: T0 0.034, T1 0.472, T2 0.456, **T3 0.252 against the incumbent's own 0.261** — the
von Mises mixture reproduces the retrieved pool's torsion diversity almost exactly.

**Geometric validity.** Per Lane B, for torsion-built chains ω deviation, cis fraction, bond-length
and bond-angle deviation are **exactly zero by construction** and Cα–Cα is 3.80 Å by construction;
they are named and omitted rather than printed as a table of structural zeros. **The only axes
carrying information here are Ramachandran and clashes:**

| arm | rama_favoured | clash_frac | clash/struct |
|---|---|---|---|
| `T0_helix` | 0.9806 | 0.0087 | 0.464 |
| `T1_blind` | 0.9510 | 0.0153 | 0.900 |
| `T2_restype` | 0.9547 | 0.0151 | 0.886 |
| **`T3_pool`** | **0.8980** | 0.0077 | 0.464 |

> **Lane B's warning fired, and it fired on my best arm — not the one I expected.** I had flagged
> `T0_helix` as this lane's damped sampler; it is `T3_pool`. The von Mises mixture with capped
> concentration places density *between* the basins it was fitted to, and it is simultaneously the
> **best-RMSD arm (3.1752)** and the **worst-Ramachandran arm (0.898)** — 10% outliers, against
> 0.951–0.955 for arms that resample real torsion pairs directly and 0.955 for the incumbent's own
> sets. **The safe arm is the absurd one**, exactly as predicted. Anyone reusing `T3_pool` as a
> baseline should know that ~5–6% of its Ramachandran validity is bought by smoothing.

---

## C8 — A POST-HOC STRATIFICATION, EXPLICITLY NOT PROMOTED.

Requested by the coordinator *after* the run, on Lane A's L8 finding that the permitted corpus is
**H 0.375 / E 0.022 / C 0.604** — β-sheet essentially absent, and mechanically so. Native SS by a
**CA-only proxy** (i→i+3 and i→i+4 spans; DSSP is not computable, only `nat_ca` exists for the dev
natives). The predicted quiet failure is visible:

| class | n | incumbent | pool best | `T0_helix` | `T1_blind` | `T2_restype` | `T3_pool` |
|---|---|---|---|---|---|---|---|
| helix | 49 | 1.9145 | 0.9832 | +0.00 | +0.49 | +0.39 | +0.21 |
| extended | 34 | 3.8918 | 1.9249 | **+1.68** | +0.17 | +0.18 | +0.16 |
| coil | 43 | 3.6734 | 2.3706 | +0.84 | **−0.12** | **−0.12** | +0.01 |

The aggregate null is **not uniform**: `T2_restype` is +0.387 on helix (`worse`, 27W/33L) and −0.121
on coil (27W/16L). And `T0_helix` at +1.68 on extended targets validates the proxy — a pure helix
generator should fail on strands, and does.

> **This is NOT a result and I am not promoting it.** The stratification was not pre-registered, it
> was requested after the numbers existed, every cell that favours my arms sits in the **0.7–1.3×
> MDE Type-M zone** (coil `T2_restype` at 0.82× MDE), and directional subgroup structure inside a
> null aggregate is exactly what this project's own methodology warns about. It is recorded as a
> hypothesis for a future sprint — *the corpus's helix bias may make blind resampling relatively
> better on coil targets* — and nothing more.

---

## C9 — WHAT IS LIMITING, AND WHAT I RECOMMEND.

**Limiting, in order:**

1. **The selector, not the source.** Selection installs +0.166 of alignment with the distogram's own
   error in every arm including a zero-information one, and the selected members carry it at cosine
   0.708 before averaging. Every candidate source scored by this functional lands on the same bias.
2. **The score's preference is not aligned with structural quality.** Supplying candidates the score
   *prefers*, in quantity, up to 52% of the final set, does not move the answer.
3. **The conditioning channel is empty.** Sequence buys ~0 at this length; chain correlation is 7.2%
   of its maximum and worth 0.26× MDE; the distogram is the selector's own referent, so conditioning
   generation on it imports the common mode rather than new information.

**Not limiting, and measured so:** the torsion representation (+0.0041 Å), the model class, sample
diversity (zero collapse on every axis), and geometric validity.

**Recommendations.**

- **Do not train a from-scratch sequence-conditioned generator.** The evidence is C3.
- **The ideal-geometry torsion manifold is free (+0.0041 Å)** — a useful, transferable result for
  anyone working in torsion space, including Lane B's residuals.
- **Any lane whose samples pass through the shipped score should expect C4 to apply**, and should
  measure its own scored-vs-unscored cosine ladder *with a zero-information control* before quoting
  a decorrelation number. A raw cosine without that control is uninterpretable, and `T0_helix` shows
  why.
- **The open direction is the SELECTOR, not the generator.** The one measured quantity with headroom
  is that generated sets contain structures the score ranks highly and nativeness does not reward.
  That is directive §28 / Lane D's ground.
- **Register `T3_pool` as the baseline any future learned torsion model must beat** — with the
  Ramachandran caveat in C7 attached.

**Caveats that travel with all of the above.**

- The fork enumeration in `PREREG_C.md` was written by **this lane, which has a stake**. Per Rule 0
  that is a weakening and is declared as one; **Lane E should re-enumerate.**
- C8 is post-hoc and unpromoted.
- The `T0_helix` jitter width (15°) and the `dup_tol` (0.10 Å) were not varied; both are stated
  rather than tuned, and no conclusion rests on either.
- The CA-only SS classifier in C8 is a proxy, not DSSP, and is stated as one.
