# S30 LANE L, TOPIC 1 -- THE REFERENCE-STATE PROBLEM AT PEPTIDE LENGTH

Brief priority 1: native-free structure scoring at 9-16 aa, the recognition problem, "knowledge-
based potentials and their reference-state problem". Serves lane R.

S29 closed the *empirical* question: no published QA method has been trained or benchmarked below
40-50 residues (S29-L1), and pLDDT has no within-target skill on 588 peptides of 10-40 aa
(McDonald 2023). **This file answers the mechanism question S29 did not ask: why does the
length wall exist, and is it an artefact of the field's training sets or is it forced?**

Answer: it is forced, by the reference state, and the size of the effect is computable in closed
form.

---

## 1. THE ALGEBRA

A knowledge-based / statistical potential is a log-odds against a reference state:

        u(d) = -kT ln [ P_obs(d) / P_ref(d) ]

so a candidate's total score contains a **separable** term

        E(x) = -kT sum_pairs ln P_obs(d_ij)   +   kT sum_pairs ln P_ref(d_ij)
                                                  ^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                  depends on the candidate ONLY through
                                                  the distances, via a fixed function

For the ideal-gas / infinite-system reference `P_ref(d) ~ d^2`, that second term is
`2kT sum_pairs ln d_ij`. **Scale the whole structure by lambda and it changes by
`2 kT N_pairs ln lambda`, exactly.** It is a pure compactness term with no shape content.

For DOPE's reference state -- **Shen & Sali, Protein Sci 15:2507 (2006)**: non-interacting points
uniform in a ball of radius `a`, with (their words) *"the radius of gyration Rg as the sample
native structure; thus, a = sqrt(5/3) Rg"* -- the reference density is

        f(d; a) = (3 d^2 / a^3) ( 1 - (3/4)(d/a) + (1/16)(d/a)^3 ),     0 <= d <= 2a

(standard result for two uniform points in a ball; I verified it integrates to 1 over [0, 2a]).
This has the exact form `f(d;a) = (1/a) g(d/a)`, so

        ln f(d; a) = -ln a + ln g(d/a)

The `ln g(d/a)` part is **exactly scale-invariant**; the `-ln a` part is a per-candidate constant.
**Conclusion: if the reference state's size parameter tracks the candidate's own size, the
potential is scale-invariant by construction. If it does not, the potential carries a pure
compactness term whose size is computed below.**

## 2. THE ARITHMETIC (closed form; NOT a measurement)

Evaluated at this instrument's chain lengths, with `Rg = 2.2 n^0.38` -- the same folded-protein
scaling law `s27/ham_lib.py` uses for the RG_LAW channel. Script:
`s30/lit/s30_L_refstate.py` (arithmetic only, no project data read).

```
    n      Rg       a     2a = the reference state's ENTIRE SUPPORT
    9    5.07    6.55     13.09 A
   13    5.83    7.53     15.05 A
   16    6.31    8.15     16.29 A
   30    8.01   10.34     20.69 A
  150   14.77   19.07     38.13 A
```

**DOPE is tabulated to a 15 A cutoff.** For a 13-mer the reference state's support ends at
15.05 A -- the table and the reference run out together. For a 9-mer the support ends at
**13.09 A, so the top 13% of the tabulated range has ZERO reference density**: the potential is
being evaluated where its own reference says the configuration cannot occur. For a 150-residue
protein the support is 38 A and the 15 A cutoff sits comfortably inside.

Finite-size correction factor `R(d;a) = f_ball/f_ideal = 1 - (3/4)t + t^3/16`, `t = d/a`:

```
    n       4A       6A       8A      10A      12A      14A      15A
    9    0.556    0.361    0.197    0.077    0.010    0.000    0.000
   13    0.611    0.434    0.278    0.150    0.058    0.007    0.000
   16    0.639    0.473    0.323    0.195    0.095    0.028    0.009
  150    0.843    0.766    0.690    0.616    0.544    0.474    0.440
```

At 8 A -- near our CONTACT channel's 7.5 A cutoff -- the correction is **0.278 at n=13 against
0.690 at n=150, a factor of 2.5**, i.e. about 0.9 kT per pair of systematic mispricing for any
potential whose reference was fitted at protein size.

**The number that decides it.** Contract a 13-mer uniformly by 10% -- pure scale, zero shape
change -- and a FIXED-reference score moves by **-0.256 kT per pair**. Over the ~55 scored Ca
pairs of a 13-mer that is **about 14 kT of reward for being smaller, with no shape content
whatsoever**. A size-matched reference moves by **exactly 0.0000**, by scale invariance rather
than by fitting.

## 3. WHAT THIS EXPLAINS

**S29-L50's loadings become algebraic rather than empirical.** Lane T measured, native-free on
126 targets x 500 members, every fold CI excluding zero:

        LEG_compactness +0.956   RG_LAW +0.921   POOLGO +0.672
        LEG_solvation   +0.621   LEG    +0.606   DSSPHB +0.587

Those were reported as a measured pattern. Section 1 makes them **forced**: any distance-based
statistical potential with a fixed reference contains a separable pure-scale term, so a high
|rho(channel, Rg)| is the expected value, not a finding. S29's framing of its own result should be
corrected accordingly -- the surprise in S29-L50 was never the loadings, it was that the loadings
**fail to explain in-band skill** (F1a +0.083 against a +0.40 bar), and that half is untouched and
still stands.

**It explains the field's 40-50 residue wall mechanistically.** DOPE's authors write plainly that
*"DOPE, like other statistical potentials, is less accurate for smaller proteins"*. Section 2 says
why: at 9-16 residues the molecule's diameter is comparable to the potential's own interaction
cutoff, so the reference state is boundary-dominated and the long-range bins carry no reference
density at all. The wall S29-L1 found in the training sets has a reason underneath it.

**It is a published, author-validated comparison, not my inference.** Shen & Sali's own ablation:
a FIXED reference sphere (their DOPE-24) performs *substantially worse* than the size-adaptive
DOPE. And a separate source states the failure direction explicitly: *"a reference state that is
too small results in an erroneous preference for loosely packed structures"* -- i.e. the reference
state's size parameter **is the knob that sets the potential's compactness preference**. It is
tunable, and in our channels it is untuned.

## 4. WHAT OUR CHANNELS ACTUALLY DO -- CHECKED, NOT ASSUMED

From `s27/ham_lib.py` (docstring, lines 24-45):

| channel | reference state as built | size-corrected? |
|---|---|---|
| DISTPOT | *"Sippl / DFIRE-style ... reference = separation-only"* | **no** -- absorbs chain connectivity, not size |
| CONTACT | *"Miyazawa-Jernigan-style ... quasi-chemical"* | **no** -- absorbs composition, not size |
| ENV | `-log p(nbr bin | aa) / p(nbr bin)` | **no** -- neighbour counts are size-driven |
| CAGEO | universe `p(theta_i, tau_i)` | n/a -- angles are scale-free already |
| RAMA | per-fold `p(aa, phi, psi)` table | n/a -- scale-free |

So the two families that carry the artefact are exactly the pair-distance and burial families, and
the two torsion/angle families are already immune. **That immunity is very likely why LEG_torsion
(+0.181 partialled in-band, fold CI [+0.053, +0.236]) and CAGEO (+0.216) are among S29-L50's nine
survivors**, and it is a cleaner explanation than "torsion happens to be orthogonal".

One honest qualification: our channels are fitted per target on that target's own leakage-safe
window universe, whose members are of comparable length. So the reference is already at roughly the
right *pool-average* size. What is NOT absorbed is the **per-candidate** variation in Rg within a
pool -- and that residual is exactly the in-band compactness axis. The artefact is therefore
smaller than the n=13-vs-n=150 table suggests and is concentrated precisely where in-band ranking
happens.

## 5. THE CONSTRUCTIVE PROPOSAL, AND ITS PRICE

**Build a pair channel with a per-candidate reference `a_i = sqrt(5/3) Rg_i`.** The result is
provably scale-invariant, i.e. compactness-free **by construction** rather than by partialling.
Two uses, only one of which is worth anything:

1. **As a band statistic -- this is the real value.** S29-L19's urgent warning was that a
   compactness-loaded band statistic crushes `rho_SY.R` to zero BY CONSTRUCTION and yields a
   self-fulfilling null. A scale-invariant statistic has provably zero uniform-scale loading, so it
   is the correct band axis. This use does not pass through the terminal operator at all.
2. **As a deployable ranker -- worth approximately nothing, and I will not pretend otherwise.**
   It creates no new information; it is a re-parameterisation of an existing structure->score map.
   By S29 section 5.1 it reaches the endpoint only as a cosine, and by section 5.4 the averaging
   readout spends at most 0.04 of any ranking.

**The trade is explicit and is against us in one direction.** ANDIS (Yu et al., Bioinformatics
35:1499, 2019): *"native recognition and decoy discrimination cannot be optimized simultaneously
with the same parameter sets"*. Removing the size term removes real BETWEEN-band signal to buy
IN-BAND interpretability. That is a choice, not a free win, and any arm that makes it must report
the between-band loss as well as the in-band gain.

## 6. THE WARNING FOR LANE R, WHICH IS THE POINT OF THIS FILE

Lane R is testing whether nativeness is recognisable from a single structure's geometry. If the
scorer is a fixed-reference pair potential, **the quantity being measured is contaminated by a term
that reads size and nothing else**, at the magnitude in section 2. Then:

- a **positive** result is unfalsifiable as "recognition", because "the native is the right size"
  and "the native is native" are not separated; and
- a **null** is the S29-L19 self-fulfilling null in a new costume.

Minimum discipline: report `rho(score, Rg)` beside every recognition number. Better: build the
size-matched variant and report both.

## 7. WHAT I REJECTED IN THIS TOPIC, WITH REASONS

| source | why rejected |
|---|---|
| Fine-grained statistical torsion-angle potentials (sub-region Ramachandran, 137 regions x 20 aa) | the discrimination is **sequence-conditioned** Ramachandran, and that channel is measured dead here: memory `phi-carries-no-sequence-signal` -- full sequence context predicts phi at 36.1 deg vs 36.4 deg sequence-blind |
| the omega torsion as an extra discriminator | omega is ~180 deg by construction in the built chain; no information |
| neighbour-dependent Ramachandran (Ting & Dunbrack 2010) | same closure as row 1; the conditioner is sequence |
| Ca pseudo-torsion / pseudo-angle (theta, tau) potentials, reported to beat DFIRE/dDFIRE/RWPlus on model selection | **already built**: `s27/ham_lib.py` CAGEO is exactly this ("CA virtual-angle / virtual-torsion statistics ... Levitt-style CA-trace potential"), and S29-L50 already prices it at +0.216 partialled in-band. Checked the code before proposing it -- this would otherwise have been a re-import |
| DFIRE's `r^1.61` finite-size exponent (Zhou & Zhou 2002) | the exponent was **fitted** to protein-size spheres, not derived; at 9-16 aa the correct exponent differs and the published value carries the artefact rather than fixing it. NOTED as corroboration of section 1, rejected as an import |
| GOAP, and the general "use molecular volume instead of Rg" variant | same mechanism as DOPE, no additional information; and volume is a worse-conditioned estimate than Rg at this length |

## 8. WHERE I COULD BE WRONG

- The -0.256 kT/pair scale-sensitivity figure is computed for **DOPE's ball reference**. Our
  DISTPOT uses an empirical separation-conditioned reference, not a ball. The *mechanism* transfers
  exactly (any reference that does not track candidate size leaves a scale term); the *magnitude*
  is DOPE-specific and I have not computed the DISTPOT-specific one.
- Section 4's claim that the torsion channels' scale-freeness explains their survival in S29-L50 is
  an **inference, not a measurement**. It is consistent with the data and it predicts that any
  other scale-free channel should also survive partialling -- which is testable against the nine
  survivors and I have not run that test.
- `Rg = 2.2 n^0.38` is a folded-protein law applied to peptides. Real 9-16mers are less compact, so
  the true `a` is larger and the effects in section 2 are **overstated** in magnitude while correct
  in direction and ordering. Using a measured per-target Rg would sharpen this and I did not do it.
