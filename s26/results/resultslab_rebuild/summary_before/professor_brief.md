# PROFESSOR BRIEF — quantum-assisted peptide structure prediction

Prepared for a research meeting. Technically accurate, no overselling. Every number in this brief is
measured on the same instrument and states its basis. Commit `a15406c`.

---

## The problem

Predict the backbone of a short peptide — 9 to 16 residues — as a Cα trace, and measure it as
full-chain Cα-RMSD to the deposited structure after optimal superposition.

Short peptides are a genuinely hard regime and not a scaled-down version of protein folding. They are
often conformationally heterogeneous, their structures are frequently stabilised by context that a
9-mer window does not contain, and the sequence–structure channel that large folding models exploit is
very weak at this length: in this project, **the full sequence context predicts φ to 36.1° against
36.4° for a sequence-blind Ramachandran marginal.**

**Instrument:** 126 cluster-disjoint development targets, 5 pinned folds, a sealed 60-target benchmark
that was never opened for tuning.

---

## The method

```
sequence → ESM-2 embedding → distance posterior (distogram)
              +
        BLOSUM retrieval of 500 real protein windows
              ↓
        Bayes-risk score of each window against the posterior
              ↓
        CVaR-VQE selects a subset          ← the quantum component
              ↓
        uniform coordinate average
              ↓
        ideal-geometry projection → final structure
```

Two genuine physical Hamiltonians are available at the scoring stage and are independently evaluable:
**Legacy**, an 11-term empirical potential, and **AMBER** ff14SB/GBn2 through OpenMM. Neither is in the
production score, for reasons measured rather than assumed (below).

---

## The result

| | value | what it is |
|---|---|---|
| **Production** | **3.2148 Å** | the ideal-geometry chain the system emits — a real, valid, renderable structure |
| intermediate | 3.0483 Å | the raw coordinate average — **22.3% contracted, not a protein structure** |

The coordinate average has a mean virtual Cα–Cα bond of 2.9614 Å against a native 3.8122 Å, with a
worst case of 0.649 Å — shorter than a covalent C–C bond. It cannot be written as a valid PDB. The
project's own code labels that arm `"raw average (illegal)"`. **We report the built chain as the
result and keep the point cloud as a labelled intermediate.** That choice costs 0.166 Å against every
number in our earlier reports, and we made it anyway.

For scale: the best single window in the retrieval pool averages **1.71 Å** (an oracle, unreachable),
and a random 75-window subset averages **3.43 Å**.

---

## What is actually quantum

**Genuine, and verified independently this sprint:** an exact dense statevector simulation of a
21-parameter real-amplitude circuit (7 qubits, 3 layers), a Hamiltonian-derived objective
`F(p) = E_p[E] − T·H(p)` with a real CVaR of the state's own distribution, and Adam on the **exact
parameter-shift gradient** — cos 1.000000000 against finite differences, relative error 4.6e-10.

**And here is what it contributes, stated plainly.**

The CVaR tail's support is always a *subset of an initial prefix of the energy order*. The trained
state can **delete** a member of the classical top-*m*; it can never **add** one from outside it. That
is a theorem about the algorithm, not a property of our data.

The optimiser genuinely trains: it beats best-of-200 from the untrained circuit at every temperature,
closing **78–89%** of the available free-energy gap. But it lands **0.902 nats and 45% of its
probability mass** away from its own analytic optimum — which for this objective is exactly the Gibbs
state `exp(−E/T)/Z` — and **the endpoint cannot tell the difference** (0.24× the minimum detectable
effect).

> **The readout is insensitive to a distributional difference of nearly half the mass.**

That single fact explains why an exact classical Gibbs state, a trained quantum state 45% away from it,
and a flat average over the top 64 candidates all land within a few hundredths of an Ångström of each
other. The binding constraint is **readout slack** — and it would bind identically on a problem where
the Gibbs state were expensive to compute, so "the target happens to be classical and cheap" is *not*
the explanation.

**One more structural fact, and it explains five sprints of null results.** The Hamiltonian is built
from the *scores of the already-sorted top candidates*, so its spectrum is the standardised rank ladder
— identical across targets to within tie-averaging (worst deviation 1.18% of range). There are
effectively **two trained states in the entire deployment, not 126**, and all per-target information
enters through which candidate occupies which rank. **A more expressive ansatz has nothing to be
expressive about**, which is why deeper circuits and larger bond dimension were measured to order
nothing, five times, across five sprints.

**What we do not claim:** we previously reported "CVaR is worth +0.113 Å". That is withdrawn. It was
measured at a temperature the pipeline does not deploy, and as a paired contrast it was 0.51× its own
MDE with a confidence interval spanning zero — null by our own standing rule. It had been labelled off
marginal means where a paired statistic was required.

---

## What Legacy and AMBER contribute

Nothing positive, on this pool, and the measurement is unusually clean.

Used as the selector, with everything else held identical and genuine CVaR-VQE as the matched selector
for all seven configurations:

| configuration | mean Cα-RMSD |
|---|---|
| Distogram alone | **3.058** |
| AMBER + Distogram | 3.132 |
| Legacy + Distogram | 3.215 |
| Legacy + AMBER + Distogram | 3.253 |
| *random 75-subset* | *3.425* |
| Legacy + AMBER | 3.674 |
| Legacy alone | 3.755 |
| AMBER alone | 3.881 |

> **Both physics energies are measurably *worse* than a random subset of the same pool** — Legacy by
> +0.330 and AMBER by +0.455, each on 5 of 5 folds.
>
> **And permuting a physics channel — destroying its correspondence to candidates while preserving its
> marginal distribution exactly — *improves* the endpoint everywhere.**

This is not a claim that the force fields are wrong. It is a claim about what they rank *on this pool*:
these are real protein windows, all of them physically plausible, and within that band the energies
order candidates along an axis that is not nativeness. Legacy and AMBER respond to compactness with
**opposite signs** (their top-75 sets differ by −0.758 Å and +1.103 Å in radius of gyration), they
correlate with each other at only ρ = −0.09, and AMBER is **orthogonal to the distogram**
(ρ = −0.019, CI [−0.058, +0.021]).

We tested the physics channel as a filter, as a partition, as an additive score term, and as a
per-target audit supplying only the sign at inference. **All five forms are closed.** In the last of
them, the real signal and a version of it with all target-correspondence destroyed were
indistinguishable: +0.0049 against +0.0053.

---

## The main scientific finding

**The bottleneck is the accuracy of the distance prior, and it is not reachable by consuming that prior
differently.**

An oracle experiment that interpolates the posterior toward the truth — holding the pool, the score,
the selector and the readout fixed — moves the endpoint at **−2.15 Å per unit of progress**, and
2.2% of the way would reach 3.0 Å. That is a transfer function, not an achievement, and it comes with a
sharp caveat: it holds *along the native's own direction*. A real operator that travels 25% of the way
to truth at a 60° angle is worth **+0.024 Å**.

We then tested every way we could think of to consume the existing posterior better: recalibrating its
width (it is over-confident by 2×), shifting its location globally and per residue-separation,
mode-seeking on the 24% of pairs that are multimodal, changing the risk functional, projecting the
distance field onto the metrically realisable set, and removing a discretisation artefact we
discovered along the way — the score's effective per-pair target takes only 17 distinct values, with
gaps up to 4 Å.

**For three of four families, the best parameter chosen with complete leakage on the very targets it
scores is the shipped default.** Across 51 arms spanning a wide range of progress-toward-truth, the
correlation between that progress and the endpoint is **+0.054**.

> You cannot buy accuracy by re-reading the same distribution, however much closer to the truth you
> move its summary statistic. If sub-3 Å is reachable here, it is through a genuinely better distance
> predictor — an achievability question about a trained model, not a question about the pipeline that
> consumes it.

---

## What remains unresolved

1. **Whether the distance prior can actually be improved**, and by how much, from information already
   available. We measured what a better prior would *buy*, not whether one is obtainable.
2. **The 2/60 benchmark leak.** An identity-normalisation defect (`core/data.py:188`, normalising by
   the longer sequence — its own docstring calls it leaky) means 4 of 126 development targets and 2 of
   60 sealed-benchmark targets carry a verbatim self-copy in their own fold model's training set. On
   the development cases the leaked window is *worse* than the pool's own best on 3 of 4. **We did not
   quantify the benchmark effect, because doing so requires opening the sealed benchmark.** It is
   declared on every benchmark figure.
3. **Cis-peptides.** The ideal-geometry projection uses a constant 3.804 Å virtual bond and cannot
   represent a cis-peptide (cis-proline is ~2.9 Å).
4. **The dynamical Lie algebra of the ansatz** is not measured. We measured gradient variance directly
   instead (no exponential decay at any width tested; the CVaR non-linearity *flattens* the decay), and
   scoped it: with 21 parameters against dim so(2⁷) = 8128 this circuit is nowhere near a 2-design at
   any width we tested, so none of it may be read as a 2-design result.

---

## How to check us

Every number traces to a persisted artefact carrying a module hash, a git commit and a completion flag
gated on a full key set. Every exported structure reproduces its own reported RMSD through the
evaluation instrument to within the PDB coordinate quantisation (worst 2.6e-04 Å). The full decision
ledger, including every retraction, is in `s25/LEDGER.md`; the quantum specification with source
references is `s25/QUANTUM.md`; the frozen pipeline is `ARCHITECTURE.md`.

**Five claims made by this project were retracted during this sprint, all by internal audit.** They are
listed with their evidence in `s25/LEDGER.md` L5, L7, L9, L11 and L15. We would rather you read those
than take the rest on trust.
