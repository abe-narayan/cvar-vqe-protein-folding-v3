# SPRINT 14 — VQE-CENTRED FOLDING REDESIGN

**Research dossier.** Sections follow the sprint brief's section 59 (A–M).

Tiering, applied to every claim: **DEMONSTRATED** (legitimate inference-time information,
paired CI and null) · **ORACLE DIAGNOSTIC** (reads the native; prices a ceiling; never a
headline) · **LITERATURE-SUPPORTED** (with citation) · **HYPOTHESIS** · **REFUTED**.

Instrument verified unchanged at sprint start and end:
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`. Benchmark untouched; pinned artefacts
retain their pre-sprint timestamps; `git status` shows zero modified tracked files.

---

## EXECUTIVE SUMMARY

**The sprint did not reach 2.0 Å, and it establishes with a mechanism why no amount of
quantum or classical search in this representation will.**

The question was whether realistic structural information could make the torsion search
space informative enough that VQE/CVaR becomes a useful conformational optimiser. The answer
has eight parts. Items 7 and 8 are the load-bearing ones: the ordering axis that would solve
this problem is learnable but changes sign from target to target, and running the optimiser
is measurably worse than not running it.

**1. The requirement is now a number, and the shortfall is an order of magnitude.** Reaching
2.0 Å through a top-100 terminal operator requires in-band pairwise ordering accuracy of
**0.638**. Measured: Legacy **0.539**, a purpose-built learned objective **0.523**, the
1-local torsion prior **0.501**, AMBER **0.463** — below chance. Every negative in this sprint
reduces to that line. It is an information gap, not a tuning gap.

The same holds on the generative side. Reaching 2.0 Å needs per-torsion accuracy near
**σ = 15°** under the i.i.d. error model real predictors satisfy; the best channel obtainable
from sequence and retrieval achieves **φ 33.6°, ψ 59.2°**, roughly twice as coarse. And the
one experimental channel that could close it is unavailable: chemical shifts exist for
**54 of 126** targets, and oracle-perfect torsions on all of them still leave the instrument
at **2.021 Å** — 55 targets are needed.

**2. One thing did get repaired, and it should be said before the negative.** Sprint 13
established that the certified global optimum of the energy landscape was *worse than a
random draw*, and that had been treated as a property of the problem. It is a property of the
**physical energies**. Over the full 262,144-configuration enumeration, the native-free
structural objective built this sprint has a certified optimum **0.885 Å better than
random**, against Legacy's 0.139 Å *worse* — roughly six times the margin of the previous
best objective. Choosing a structural rather than an energetic objective fixed it, and the
result survives certification, which is the strongest available form of the claim.

**3. Better search still cannot substitute, because the bottleneck is discrimination, not
optimisation.** The certified optimum still sits **1.876 Å worse than the best structure
inside the objective's own lowest decile** — and nothing can beat a certified optimum.
Independently, sampling from 10 to 20,000 evaluations improves the objective monotonically
and moves the emitted structure by 0.17 Å, all of it arriving by evaluation 300, while the
best structure available improves from 2.907 to 1.532 Å and the **selection gap grows to
2.04 Å**. Two independent routes — exact enumeration on nine targets and sampling on 126 —
put the irreducible loss at 1.88 and 2.04 Å. No arm, classical or quantum, has a
selection-gap advantage, and CVaR's tail weighting buys no discrimination.

So the refined statement is: **the objective now points somewhere worth going, and still
cannot tell good structures from mediocre ones once it arrives.** One caution on the mean:
per-target certified optima span 0.675 to 5.256 Å and two of nine targets anti-rank, so the
average is a mixture of a few excellent cases and two failures.

**4. The mechanism has a number, and it belongs to the energy models rather than to the
problem.** Discrimination is a function of the quality gap between two structures, not their
separation. **No physical objective exceeds 0.511 pairwise accuracy when two structures differ
by less than 0.25 Å**; Legacy needs a 1.31 Å gap to reach 55%, AMBER 1.58 Å. The useful search
range is about 2.5 Å wide, so both energies resolve roughly one bit of it.

**But a learned native-free structural objective is not bound by that floor.** On the full
enumeration the pure distogram term reaches **0.654 at a 0.25 Å gap and 0.633 at zero gap** —
it separates structures of essentially equal quality. So the floor is a property of Legacy and
AMBER specifically, which is exactly how it was scoped.

**5. Selection is an exact additive identity, and that dissolves the paradox.**

```
sel  =  pool  +  FILTERING  +  ORDERING
```

At the lowest 1% of the space: Legacy has FILTERING **−0.124** and ORDERING **+0.263**
[+0.046, +0.484], netting **+0.139 Å** — *exactly* the Sprint 13 certified-optimum-versus-random
figure. A result that stood for a sprint as an empirical curiosity is now the sum of two
measured terms. Legacy's steric term has a filtering component and **no ordering component at
all** (−0.021, CI spanning zero); AMBER has the **best filtering term measured** (−0.141) and a
positive ordering term.

**The ORDERING term's sign *is* "optimise harder, get worse", measured directly.** It is
positive for both physical energies: they actively mis-order inside the region a search
occupies. This is the single number any future objective should be judged on.

So bulk discrimination and optimum quality are **different terms of one identity**, not
competing measures of one thing — the learned term having the best of one and the worst of the
other was never a contradiction. **Neither in-decile rank correlation nor bulk pairwise
accuracy predicts selection quality**, and the sprint used both as if they did. Of the four
proxies tried, only **tail-restricted pairwise accuracy** tracks the certified optimum.

**A caution that must travel with every tail number:** because gaps inside a tail are small by
construction, accuracy falls toward chance as the tail shrinks whether or not skill was lost.
The correct baseline is a **random-tail null at 0.524–0.527, not 0.500**. The result survives
it — Legacy scores 0.462 gap-matched inside its own lowest 1% against 0.523 for a random subset
of the same size. **An energy is blindest inside its own tail.**

**6. An inversion that decides the architecture.** Global rank correlation and
certified-optimum quality move in **opposite directions** as the objective's two terms are
reweighted: adding distogram weight improves global ordering monotonically while making the
certified optimum monotonically worse. **The distogram orders the bulk; the retrieval torsion
prior places the optimum.** Combined with the Sprint 12 law that the terminal operator
consumes the set *mean* rather than the set *best*, this means the objective's weighting must
be chosen jointly with the terminal operator, and a single scalar objective-quality number is
not decision-relevant.

**7. The deepest finding: the bottleneck is not representation, capacity, search or sample
size. It is that the ordering changes sign from target to target.** A purpose-built objective
reaches **0.986** in-band pairwise accuracy on a target it has seen — far above the 0.638 that
2.0 Å requires — and **held-out configurations of that same target score identically**
(0.986, converting to −0.868 Å). There is no overfitting to explain away. **Cross-target it
collapses to 0.600**, which converts to −0.047 Å. The gap between 0.986 and 0.600 is the whole
problem.

That closes the "build a richer model" branch with a positive measurement rather than a
failure: a linear pair potential already saturates the within-target problem, so no amount of
nonlinearity, capacity or equivariance can improve on 0.986, and none of it addresses transfer.
It converges with two independent diagnostics — a learning curve that *declines* (one or two
training targets beat twelve) and per-target skill correlating **+0.909** with the native's
radius of gyration. **The in-band discriminating axis is real, learnable, and per-target; its
correct sign is a property of the individual target that inference cannot observe.**

**8. The definitive negative: running the VQE is worse than not running it, 0/12.** Against the
only fair control — **the best of the same number of shots drawn from the untrained circuit** —
optimisation loses by **+0.65 to +1.32 Å** and wins 0 times out of 12 on four of six
initialisation arms. Best-of-N sampling is free, and the shipped pipeline already does it.

The mechanism is item 6 seen from the outcome side: optimisation concentrates the distribution
onto the objective's low-energy tail, every objective is at chance inside that tail, and the
returned structure is therefore a random draw from it. **The optimiser trades a broad sample
whose best member is good for a narrow one whose mean it samples blindly. Concentration is
exactly the wrong move when discrimination is the binding constraint.**

**And the warm-start confound is answered in the direction least convenient for VQE.** After
optimisation the result is 2.42–2.76 Å regardless of where it started, from initialisation
means spanning 2.966 to 4.131 Å; an oracle warm start at q=0.95 ends at 2.578 against a random
start's 2.606. **The optimiser erases its own initialisation.** The better the oracle start,
the *worse* the comparison — a better start has a better best-of-N while the VQE converges to
the same place anyway.

**What this means for the quantum question.** A variational optimiser is a better search.
Search saturates two orders of magnitude before the budget runs out, and the residual gap
survives certified global optimisation. There is therefore no route by which a better ansatz,
encoding, initialisation or CVaR schedule contributes through search quality in this
representation. That is a negative, but it is a *mechanistic* negative with a measured
threshold, and it is stronger than the Sprint 13 version because it was established on an
objective built to be good rather than on one already known to be bad.

**What was built that is worth keeping.** A native-free structural objective that orders the
space roughly fourfold better than the best physical energy without learning anything; an
exact price for physical validity (imposing ideal geometry costs 0.157 Å); a positional cost
law explaining a previously empirical correction; a coherence surface showing error *shape*
is worth as much as error *size*; and a leakage guard that mechanically rejects any emitter
which reads the native.

**Corrections carried.** Five claims were retracted during the sprint, four of them the
coordinator's own. They are listed in section L and kept beside what replaced them.

---

## A. FINAL ARCHITECTURE

### A.1 What was measured, as a data flow

```
sequence
  |
  +-- BLOSUM62 retrieval, K=500, leakage-safe library     [existing]
  |     |
  |     +-- window coordinates W        -> coordinate-space channel
  |     +-- window torsions PHI/PSI     -> NEW: position-specific torsion prior
  |
  +-- leave-fold-out distogram (ESM features)             [existing]
        |
        v
  NATIVE-FREE STRUCTURAL HAMILTONIAN            [new, s14/hamil.py]
        H(S) = (1-w) * z[E_prior(S)] + w * z[E_disto(S)]
        |
        v
  search over torsion configurations S in {0..k-1}^n
        random / greedy / annealing / GA / exact enumeration / VQE / CVaR-VQE
        |
        v
  low-energy SET (top m)                        [s14/consensus.py]
        |
        v
  coordinate-space consensus + projection to ideal geometry
        |
        v
  AMBER validation flag + single refinement pass          [validator, never a score]
```

### A.2 The architecture is reported, not recommended

Every stage above was measured and the composition **loses to the incumbent retrieval
pipeline** by +0.110 Å [+0.004, +0.214] at its best setting. It is documented because the
decomposition is informative, not because it should ship. The honest recommendation is in
section M.

### A.3 Component roles, as measured

| component | measured role | evidence |
|---|---|---|
| retrieval coordinate average | the strongest operator in the system | beats torsion averaging of identical windows by +1.024 Å |
| retrieval torsion prior | best torsion channel in the project; ordering skill conditional on the prior-typical region | φ 33.6°, ψ 59.2°; in-decile +0.046 uniform / +0.401 prior-sampled |
| distogram Bayes risk | orders the whole space; the durable half of the objective | global ρ +0.442 uniform, in-decile +0.161 |
| Legacy energy | **clash gate, not a ranking function** | 98.8% of variance is `steric`, constant across its own lowest decile |
| AMBER ff14SB/GBn2 | **physical validator; refinement safe but inert** | best clash detector measured; refinement +0.014 Å [−0.001,+0.028] |
| projection to ideal geometry | **costs 0.157 Å — the price of physical validity** | +0.157 Å [+0.124,+0.191], 18W/108L |
| coordinate consensus over a set | recovers 0.29 Å from members that do not improve | argmin 3.608 → m*=200 consensus 3.314 |

---

## B. MATHEMATICAL FORMULATION

### B.1 Representation

A configuration is an integer vector `S ∈ {0,…,k−1}^n`, one discrete torsion state per
residue, drawn from the sequence-conditioned leakage-safe library
`torsion_lib2.library_for(seq, k, exclude_seq=seq)`. Decoding is exact:
`(φ_i, ψ_i) = (PHI[i, S_i], PSI[i, S_i])`, and the backbone follows by ideal-geometry
construction.

**Live qubit count.** Nominal is `n·log₂k`. The true count is **`(n−2)·log₂k`**: four
torsions per chain — `φ₀, ψ₀, φ_{n−1}, ψ_{n−1}` — are inert for the CA trace, so **both
terminal residues are invisible to CA-RMSD**. This follows from the locality theorem rather
than from an implementation detail: `d_ij` depends on the residues *strictly between* i and
j, and residues 0 and n−1 are strictly between no pair. At n = 9, k = 4 the live count is
**14 of 18**; the k = 4 operating point across the instrument is **21.9 mean live qubits, not
25.9**. Found independently by two agents. Note `φ_{n−1}` still places C/CB/O, so Legacy and
AMBER depend on states the CA metric cannot see.

### B.2 The structural Hamiltonian

Both terms are diagonal in the computational basis and are standardised per target over the
same configuration sample, so `w` is interpretable and arms differ only in `w`:

```
H_w(S) = (1 − w) · z[ E_prior(S) ] + w · z[ E_disto(S) ]

E_prior(S) = − Σ_i log P_i(S_i)                     1-local, exactly diagonal
E_disto(S) = mean_{(i,j)} risk_ij( d_ij(S) )        many-body
```

`P_i` is the native-free position-specific state occupancy built from the retrieved windows;
`risk_ij` is the shipped leave-fold-out distogram's Bayes-risk function over CA–CA distances.

**Locality.** `E_prior` is exactly 1-local. `E_disto` is not 2-local and cannot be made so:
by the Sprint 13 locality theorem `d_ij` depends on exactly the `j−i−1` residues strictly
between i and j, so a separation-8 pair is a **14-qubit interaction at k = 4**. No 2-local
Ising form of a distance-based molecular objective exists in this encoding. This is a real
cost of the representation and is stated rather than hidden.

### B.3 CVaR

For an energy `E` over basis states with probabilities `p`, at level α, CVaR is the
conditional expectation over the α-quantile lower tail, with the partial bucket at the
quantile boundary weighted so total mass is exactly α. Implementation and gradient
correctness were re-derived and re-verified this sprint rather than assumed; see section I.

### B.4 Normalisation

Raw AMBER spans 16 orders of magnitude. Terms are put on comparable scales by gradient
standard deviation (`s14/cache/ener_norm.json`, key `grad_sd`) — **never by raw AMBER's own
standard deviation**. A monotone transformation leaves rank correlation unchanged while
collapsing dynamic range, which is the proof that the damage in the physical energies is in
the **ordering, not the scale**, and therefore that no rescaling repairs it.

---

## C. EXPERIMENTAL MATRIX (coordinator workstream)

Every arm below is scored on the same 126-target instrument, by the same metric, against the
same reference, through `s14/ladder.py`. Every native-free emitter passes `audit_emitter`,
which replaces the native trace with noise and requires bit-identical output.

| # | experiment | module | question | verdict |
|---|---|---|---|---|
| C1 | level-0/1 emission ladder | `ladder.py` | do the Sprint 13 anchors reproduce? | yes, to the third decimal |
| C2 | retrieval-conditioned torsion prior | `retprior.py` | is the retrieval pool a better torsion channel than a trained sequence model? | **yes**, and it still loses |
| C3 | MAE versus emitted accuracy | `retprior.py` | does angular error order emitted accuracy? | positively but imperfectly, ρ +0.467, with real inversions |
| C4 | error-coherence surface | `coherence.py` | at matched σ, does error *shape* matter? | **yes, up to 0.93 Å**; real predictors are near-i.i.d. |
| C5 | Sept 2026 preprint audit | `lit_FINDINGS.md` | does arXiv:2609.02113 change the direction? | no; it corroborates our central negative |
| C6 | positional cost of a torsion error | `position.py` | where along the chain does an error cost? | symmetric hump, mid-chain 4.0× the ends |
| C7 | native-free structural Hamiltonian | `hamil.py` | can a non-energy objective order the space? | **yes, ~4× better than any physical energy** |
| C7c | proposal + sequence-blind controls | `hamil_control.py` | is that ordering an artefact of the proposal? | **partly — headline retracted**; retrieval conditioning survives |
| C8 | budget curve on a good objective | `budgetcurve.py` | does searching a good objective harder help? | **no. saturates by evaluation 300 of 20,000** |
| C9 | set aggregation versus argmin | `consensus.py` | should a good objective select a set rather than a point? | yes, worth 0.29 Å; still loses by +0.110 |
| C10 | averaging space and projection cost | `avgspace.py` | does *where* you average matter, and what does validity cost? | +1.024 Å and **+0.157 Å** respectively |
| C11 | classical limit of a prior-reproducing VQE | `vqe_classical_limit.json` | what must a quantum claim be differenced against? | **3.485 Å**, not 4.072 |
| C12 | three-term objective at scale | `hamil3.py` | does `leg_torsion` replicate at n=126? | see section F |

## D. RMSD RESULTS — the complete native-free ladder

All values mean CA-RMSD in ångström **over the 126-target instrument**. ORACLE rows read the
native and are diagnostics, never results.

**Do not compare these against the certified numbers in sections E and I.** Those are measured
on the **nine fully enumerated targets**, which are all n=9 and systematically easier — the
certified optimum of 2.896 Å and the sampled argmin of 3.572 Å here are different quantities
on different populations, not a contradiction. Every certified figure in this dossier carries
its scope.

| arm | mean | median | <2 Å | FAIL18 | vs incumbent |
|---|---|---|---|---|---|
| ORACLE torsion-space ceiling, k=4, random start (S13) | 1.982 | — | — | — | — |
| ORACLE best of 4,000 sampled configurations | 1.647 | — | — | — | — |
| ORACLE best member of the top-200 selected set | 2.647 | — | — | — | — |
| coordinate average of top-75 windows, **not a valid backbone** | 3.048 | 2.837 | 0.29 | 5.832 | −0.157 |
| **incumbent retrieval pipeline** | **3.204** | 2.966 | 0.28 | 6.026 | reference |
| structural-Hamiltonian set consensus, m*=200 | 3.314 | — | 0.27 | 6.057 | +0.110 [+0.004,+0.214] |
| prior sampling + coordinate consensus, B=200 | 3.485 | 3.299 | 0.24 | **5.771** | +0.281 [+0.113,+0.461] |
| top-75 similarity-weighted torsion mean | 3.514 | 3.309 | 0.25 | 6.180 | +0.310 [+0.128,+0.492] |
| structural-Hamiltonian argmin, 20,000 evaluations | 3.572 | — | 0.25 | 6.396 | +0.368 [+0.257,+0.482] |
| top-75 state argmax, k=4 | 3.842 | 3.436 | 0.28 | 6.725 | +0.638 [+0.401,+0.884] |
| class prior argmax (1-local closed-form minimiser) | 3.969 | 3.948 | 0.24 | **5.670** | +0.765 [+0.514,+1.024] |
| pool-500 state argmax, k=4 | 4.004 | 4.147 | 0.29 | 5.837 | +0.800 [+0.525,+1.080] |
| **constant α-helix, zero information** | 4.065 | 4.253 | 0.29 | **5.887** | +0.861 [+0.574,+1.171] |
| top-75 torsion circular mean | 4.072 | 3.552 | 0.29 | 6.743 | +0.868 [+0.551,+1.202] |
| pool-500 torsion circular mean | 4.175 | 3.795 | 0.27 | 6.219 | +0.971 [+0.643,+1.332] |
| sample the class back-off prior | 4.851 | 4.542 | 0.00 | 5.856 | +1.647 [+1.380,+1.912] |
| uniform random k=4 library state | 5.033 | 4.867 | 0.00 | 6.146 | +1.829 [+1.559,+2.091] |

**Nothing native-free beats the incumbent.** The closest is set consensus at +0.110, a CI
that barely excludes zero. The only row below the incumbent is not a valid protein backbone.

**Four low-information arms beat the incumbent on the 18 hardest targets** (bold FAIL18):
the class prior argmax at 5.670, prior sampling with consensus at 5.771, the raw coordinate
average at 5.832, and a constant α-helix at 5.887, against the incumbent's 6.026. On the
failure class, this system's information is worse than none — Sprint 12's finding, reached
here by three independent routes.

## E. VQE CAUSALITY ANALYSIS

The decomposition, all on one instrument and one candidate space:

| stage | emitted | increment |
|---|---|---|
| torsion information, one committed vector | 4.072 | — |
| + sampled 200× and coordinate-averaged | 3.485 | **−0.587** |
| + ranked by the structural Hamiltonian, top 200 averaged | 3.314 | **−0.171** |
| incumbent | 3.204 | −0.110 |

**Aggregation is worth 3.4× more than the objective.** The entire contribution available to
any objective — and therefore to any optimiser that searches one, quantum or classical — is
**0.171 Å**, and that has already been banked classically by ranking 4,000 samples.

Two independent measurements say no optimiser can add to it.

**The selection gap reproduces exactly across two independent measurements.** The VQE
workstream recomputed the coordinator's statistic on their own enumeration — same objective,
same weight, nested prefixes, independent implementation, **disjoint target set**:

| budget | 10 | 300 | 3,000 | 20,000 |
|---|---|---|---|---|
| selection gap, 9 enumerated targets | 0.637 | 1.583 | 1.971 | 1.923 |
| selection gap, 126-target instrument | 0.835 | 1.577 | 1.935 | 2.040 |

Agreement is 0.006 Å at budget 300 and 0.036 Å at 3,000, within 0.12 Å everywhere.

**It survives certification.** On the enumerated nine, the certified global optimum sits
**1.876 Å worse than the best structure inside the objective's own lowest decile.** Nothing
can beat a certified optimum, so this cannot be attributed to insufficient search.

**The grid also reproduces Sprint 13 exactly**, which validates the instrument: at signal 0 the
certified optimum is **3.920 Å against a random draw of 3.781 Å**, +0.139 Å worse — the recorded
figure to three decimal places, from independent code on the same nine spaces, and the same
+0.139 the filtering/ordering identity nets to.

**The budget trap is a property of BAD OBJECTIVES, not of search.** At full n the non-monotone
curve **inverts by objective signal 0.50**, where budget 30,000 reaches **2.149 Å against a
certified optimum of 2.153 Å** and more search is strictly better. The complete statement,
superseding both the Sprint 13 warning and this dossier's earlier framing: *searching harder
makes structures worse on a bad objective, neither helps nor hurts on a mediocre one, and helps
monotonically on a good one.* It does not rescue VQE — at the signal levels reachable with real
native-free information the curve is in the flat-to-harmful regime, and VQE loses to greedy on
the objective axis regardless — but "optimise harder, get worse" must always carry its
objective-quality condition.

**And the saturation half is target-length dependent; the coordinator's first statement
of it was too general.** On the 126-target instrument the emitted structure is flat (3.742 →
3.572). On the enumerated nine it *falls* 0.38 Å monotonically (3.297 → 2.920). Those nine are
all n=9 and see 7.6% of their space at budget 20,000, where the full instrument sees a small
fraction of a much larger one — the Sprint 13 fraction-of-space law again. The correct joint
statement is: **the discrimination loss is the same everywhere and reproduces exactly; whether
search still buys anything on top of it depends on chain length.**

**No arm has a selection-gap advantage.** Across random, greedy, annealing, genetic search
and VQE at four CVaR levels, every arm returns roughly 1 Å worse than its own best evaluated
configuration, and a tenfold budget increase does not close it. **CVaR's tail weighting buys
no discrimination.**

**And on the objective axis, VQE is decisively the worst optimiser tested.** Fraction of
cells reaching the certified optimum at budget 30,000: greedy 1-opt **0.71**, annealing
**0.71**, VQE **0.00–0.12**. On the real Legacy objective, greedy and annealing reach it in
**100%** of cells and VQE in **0%**.

So the structural tie between VQE and the classical arms is not VQE holding its own — it is
**both methods hitting the same discrimination floor from different distances**, with the
floor close enough that a much worse optimiser also reaches it. Reporting only the structural
axis would have shown parity; reporting both shows there is none. This is precisely why the
brief required the two axes never be substituted for one another.

The conclusion is not that VQE performed badly. It is that **the quantity a variational
optimiser improves is not the quantity that limits this problem.**

---

## F. AMBER VERSUS LEGACY — the controlled comparison

Run with no prejudgement in either direction, on the same conformational ensembles, over the
full A–J arm matrix. The verdict is that **they have different legitimate roles and neither
is a ranking function.**

### F.1 Legacy is a clash gate

**98.8% of Legacy's variance is a single term, `steric`.** That term is exactly zero on
73.5% of the space and **constant across its entire lowest-energy decile** — a binary gate
carrying no ordering information in the region a search actually occupies. It places the
ORACLE-snap configuration at the **62nd percentile** of its own energy on a
Ramachandran-plausible population. Coordinate descent on it takes a 1.295 Å structure to
**3.684 Å**.

Its apparent −0.296 Å selection advantage has a top-2 concentration share of 1.02: drop two
targets and it becomes **+0.009**.

**Legitimate role:** a 10%-removal geometric veto, steric term only. Above 10% removal the
total is an unsafe filter — it costs +0.219 Å of the pool best by 25% removal and +0.385 Å
by 75%.

### F.2 AMBER is a physical validator and an inert refiner

AMBER is the **best clash detector measured**. The sampled ensemble is 13.07% clashed
(minimum separation 0.23 Å); AMBER's low decile is **0.09% clashed at 1.98 Å** against
Legacy's 0.56% at 1.62 Å. It is the **only physical objective with a positive excess decoy
AUC** against a matched random-anchor null (+0.056 to +0.113). It can discard **90% of a pool
and still hold the best structure**, where Legacy has already lost 0.380 Å.

**Refinement is safe but inert:** +0.014 Å [−0.001, +0.028], 13W/32L over 45 cells. From
1.295 Å it returns 1.321 Å; from 4.334 Å it returns 4.337 Å. It neither rescues bad
structures nor damages good ones. Descent on it takes 1.295 Å to 3.685 Å, the same failure
as Legacy.

**But AMBER is a legitimate replacement for the geometric projection, and a better-quality one.**
At n=126: relaxing a genuine all-atom coordinate average under restraint k=10 emits **3.207
against the projection's 3.205** — +0.001 with a CI **0.03 Å wide**. A genuine force field, a
strictly *weaker* constraint than the ideal-geometry manifold, recovers exactly the same
accuracy: none.

**Ideal geometry was never the binding constraint.** Averaging contracts the backbone by
**−25.8%** (N–CA 1.120 against 1.458; C–N 0.864 against 1.329), and contraction predicts
per-target cost at **ρ = −0.552** at full n. The 0.157 Å is the price of restoring *any* valid
bond length to a structure averaging shrank by a quarter. Free relaxation is the worst arm
(+0.272) and that is the signature: unrestrained, the field re-expands to its own preferred
geometry and drifts 1.526 Å away, discarding the ensemble information that made the average
good.

**And at the right restraint strength it buys accuracy too — the sprint's only positive
accuracy result.** The restraint frontier does not turn; tighter restraint keeps buying RMSD.
On 50 targets against arm A at 3.171:

| restraint | RMSD | vs A | CI | W/L | geom deviation | strain vs the force field's own equilibrium |
|---|---|---|---|---|---|---|
| k=10 | 3.169 | −0.001 | [−0.023,+0.019] | 17/33 | 0.0134 | tighter than equilibrium |
| **k=30** | 3.146 | **−0.025** | **[−0.046,−0.005]** | 26/24 | 0.0181 | **at equilibrium** |
| k=100 | 3.101 | −0.070 | [−0.098,−0.043] | 35/15 | 0.0374 | 2.2× strained |
| k=300 | 3.066 | −0.105 | [−0.140,−0.073] | 40/10 | 0.0665 | 3.9× strained |

The strain column uses a **measured** yardstick rather than an asserted one: free relaxation
is ff14SB's own equilibrium geometry, at 0.0170 deviation. Beyond k=30, RMSD is bought with
strain — exactly the trap the brief demanded be guarded against, caught by the workstream in
its own favourable result. So the honest form of the frontier is **"a better number, paid for
in geometric quality"**, not a free gain.

**And the one free setting survives, at n=126.** An intermediate verdict that it failed its
concentration check was **withdrawn**, and the reason is a correction to this project's own
standing rule: **a raw drop-top threshold is not a valid concentration test.** When an effect's
mean is small relative to per-target spread, discarding the ten most favourable targets removes
a large share of the total *even if every target carries an identical effect* — so at low
signal-to-noise the test must fail a uniform effect.

Compared against its own null (a uniform effect of the same mean and sd; mean/sd = 0.28,
4,000 simulations), every statistic sits at the **62nd percentile**: drop-top-10 −0.0080
against a null of −0.0094, drop-top-20 +0.0028 against +0.0006, top-10 share 0.668 against
0.698. **No evidence of concentration in any arm.**

The symmetric caution carries equal weight: *not concentrated* is not *demonstrated uniform* —
at this signal-to-noise the test has almost no power either way. What carries the result is the
CI and **the per-fold table**: all five pinned folds same sign, −0.014 to −0.036, a **2.6×
range** (it was 25× at n=50; that spread was noise narrowing).

**Final, n=126: k=30 emits 3.183 against the projection's 3.205 — −0.022 Å [−0.036, −0.009] —
at 1.13× the force field's own equilibrium strain, zero clashes, Ramachandran 0.874 against
0.734.** So "never a better number" is retracted, and now at *valid* geometry rather than only
at strained geometry.

**Three caveats travel with it, and the third was strengthened twice.** It is `tuning126` only
and has never been near the benchmark, where nothing has ever transferred. It is **0.7% of the
baseline**. And **the restraint constant was selected on dev-set RMSD — there is no native-free
rule in hand that picks it.**

That third caveat first read "chosen with the frontier visible", implying a principled rule
merely applied at the wrong time. **The rule is degenerate**, and the workstream retracted it
themselves: free relaxation sits at 1.00x equilibrium *by definition* — it *is* the equilibrium
— so "closest to equilibrium" trivially selects the **worst** arm (3.477 A). The non-circular
form, "tightest restraint not exceeding equilibrium", selects k=10, which gives **+0.001 A, no
gain at all**. Reaching k=30 needs a 1.15x threshold, and 1.15 was chosen with the frontier
visible.

**What the strain yardstick does and does not do.** It legitimately *rules out* k=100 and k=300
as strained — that is what stopped −0.105 A becoming the headline — but it does **not** select
k=30 from {k=2, k=10, k=30}, all of which sit at or near equilibrium. **A confirmation must
pre-register the restraint constant**; a pre-registered k could land anywhere from **+0.001
(k=10) to −0.022 (k=30)**.

**Two further checks, one favourable and one not.** The baseline is conservative: arm A uses
the *better* of the two projection arms (`fit_ca` at 3.2052 against `ca` at 3.2126), and
against the other the gain would be **larger** at −0.0295 A [−0.0465, −0.0122]. And a
superlative in an earlier draft — "the first native-free accuracy improvement this project has
measured" — is **false**: Sprint 12's score-filter plus consensus medoid measured
**−0.172 A [−0.316, −0.027]**, earlier and roughly eight times larger. The workstream declined
to carry that claim on the grounds they could not verify it from their own work, and they were
right to.

### The methodology this episode cost, and what it bought

This result was published, retracted, and un-retracted inside a few hours, and both errors
were about the *check* rather than the number.

**Error one — reading a win/loss as reassurance.** The workstream and the coordinator both
called 26W/24L with a CI excluding zero "a small consistent shift", and the coordinator
published it, before either had looked at the drop-top curve. The check had not been skipped:
`I.paired` returns `drop_top10_mean_diff` and `top10_share` **in the same dict** the W/L was
read from. *The failure mode is reading past a check, not forgetting to run one*, so the fix is
structural — the checker now emits one PASS/FAIL verdict block that cannot be read selectively.

**Error two — reading an uncalibrated drop-top curve as disqualifying.** That was equally
wrong: at mean/sd = 0.28 a failing drop-top curve is simply **what a uniform effect looks
like**. The mandated rule caught `leg_torsion` correctly and misfired here, and the difference
is signal-to-noise.

> **The concentration check is NECESSARY AND NOT SUFFICIENT. It needs its own null, and
> mean/sd must print beside it so a reader can see when it has no power.**

**What survives as a free diagnostic:** a near-even W/L *together with* a CI excluding zero is
suggestive of a **concentrated** effect and worth checking against a null — a uniform effect of
that size would give a lopsided W/L and a median near the mean. The **median-versus-mean
ratio** flags it at no cost, before any drop-top number is computed. Both fixes are now in the
tooling rather than left to judgement.

The stereochemistry gain stands regardless: **zero clashes against 45.2%**, and
**Ramachandran-allowed 0.860 against the projection's 0.734**, at about 12 s per target.

**Legitimate role:** post-hoc validation and a single refinement pass at the end, emitting a
flag and a structure — **never a score**.

### F.3 Complementarity: premise correct, route refuted

Truth-partialled error correlation between the two models is **+0.096** — they are genuinely
decorrelated, so hypothesis H7's premise holds. But it is not exploitable. Optimal linear
fusion is worth **+0.032 ρ under a per-target ORACLE weight**, and the identical machinery
**loses 0.240 ρ going from in-sample to leave-one-target-out**. In the arm matrix
`rank(Legacy) + rank(AMBER)` is **0.015 Å worse than Legacy alone**.

This reproduces the project's existing arithmetic: fusion gain goes as the *square* of the
weaker channel's skill, and both channels here are weak.

**Pareto is refuted too.** Mean |conflict| is 0.222, the frontier beats a same-size random
subset by only 0.048 Å, and it **discards 72% of the truly-best 1%** of structures.

### F.4 The discrimination threshold — the number that explains everything

The naive decoy experiment is confounded: discrimination is not a function of structural
*separation* (it is flat at fixed gap) but of the **quality gap** between the two structures.

| objective | quality gap needed for 55% pairwise accuracy |
|---|---|
| `leg_steric` | 0.58 Å |
| Legacy total | 1.31 Å |
| AMBER total | 1.58 Å (and 3.87 Å for 60%) |
| the 1-local prior, `leg_torsion`, Legacy-minus-steric | never reach 55% at any gap |

**No objective exceeds 0.511 pairwise accuracy when two structures differ by less than
0.25 Å.** The useful search range is about 2.5 Å wide, so both energies resolve roughly one
bit of it.

That single threshold explains the certified global optimum being +0.139 Å worse than random,
the 1.3 → 3.7 Å descent, the ~1 Å selection gaps in every search arm, and the budget
saturation. **It is the mechanism behind the whole sprint's negative.**

## G. TORSION-INFORMATION ANALYSIS

### G.1 The requirement

Under the i.i.d. error model that real predictors satisfy, reaching 2.0 Å requires
**σ ≈ 15.1°**. The incumbent pipeline is equivalent to **σ ≈ 27.1°** (refining Sprint 13's
~29°).

### G.2 What is available

| source | φ | ψ |
|---|---|---|
| sequence-blind corpus marginal (S13) | 36.4° | 72.8° |
| best learned full-context sequence predictor (S13) | 36.1° | 62.4° |
| **retrieval top-75 circular mean (new, S14)** | **33.6°** | **59.2°** |

The retrieval pool beats the trained leave-fold-out sequence model on both angles with no
training at all, and is the best torsion channel measured anywhere in this project. It is
still roughly **twice as coarse** as 2.0 Å requires.

### G.3 Error shape matters as much as error size

At matched marginal σ, along-chain error coherence moves emitted RMSD by up to **0.93 Å**.
The σ needed to reach 2.0 Å ranges from 10.6° (positively autocorrelated errors) to 20.2°
(anti-correlated).

**But all ten real emitters sit in the near-i.i.d. band** (lag-1 between −0.092 and +0.066).
So the Sprint 12 restraint surface, which used i.i.d. corruption, **stands and is mildly
conservative**. The coordinator's hypothesis that real predictors have coherent errors is
refuted.

### G.4 Where an error costs

The cost of a torsion error is a **symmetric hump peaked at mid-chain**: the middle 40% of
the chain carries 62.8% of the total single-residue cost, the outer 40% carries 15.7%, a
**4.0× ratio**, and the terminal tenth at each end is roughly **14× cheaper** than mid-chain.

The mechanism is geometric: a torsion at position p hinges two rigid segments of length p and
n−p, and under Kabsch superposition the displacement follows the lever-arm product p(n−p).
There is no privileged terminus (front/back ratio 0.886), which refutes the coordinator's
N-terminal-dominance prediction.

**This explains the Sprint 13 correction** that terminal dropout is 0.40–0.50 Å cheaper than
uniform, converting it from an empirical surprise into a mechanical consequence. It also
means **coverage must be position-weighted**: 75% coverage concentrated mid-chain beats 90%
concentrated at the termini.

## H. ENERGY-LANDSCAPE ANALYSIS

| objective | global ρ | in-decile ρ | argmin RMSD |
|---|---|---|---|
| raw AMBER | −0.036 | −0.088 | 5.202 |
| Legacy total | +0.198 | +0.043 | 4.120 |
| `leg_torsion` alone | +0.230 | **−0.073** | 4.933 |
| retrieval torsion prior | +0.073 | +0.049 | 4.725 |
| distogram Bayes risk | **+0.443** | +0.157 | 3.879 |
| **prior + distogram** | +0.388 | **+0.367** | **3.571** |

(Global and in-decile figures under a strictly uniform proposal except the last row, whose
in-decile is proposal-sensitive; see section L.)

**In-decile ρ is not monotone in objective quality.** On a signal-tunable family it peaks
near +0.446 at moderate signal and falls to +0.340 where the objective is genuinely much
better (global ρ +0.933), because a better objective's lowest decile is a narrower and more
homogeneous set with less orderable spread. Both axes must always be reported.

**The two terms are genuinely complementary.** Global ordering rises monotonically with
distogram weight while in-decile ordering peaks at a mixture and collapses at pure distogram.
The distogram orders the whole space; the torsion prior orders the good region.

---

## K. LITERATURE POSITION

Full ledger, 25+ papers with verification status, in `s14/lit_FINDINGS.md` (1,115 lines).
The primary text of the key paper was fetched and read in full; nothing below is secondary.

**The September 2026 preprint is arXiv:2609.02113** (Cumbo et al., Cleveland Clinic, method
QTF): a logarithmic-scale VQE for off-lattice protein structure prediction in continuous
torsional space, with chignolin and Trp-cage.

| what it claims | what it is |
|---|---|
| 0.623 Å | the **argmin over 2.68 M snapshots ranked by RMSD to the native** — an oracle selection |
| the predictive result | median final model **2.90 Å chignolin / 5.39 Å Trp-cage**, two targets, terminal residues excluded, RMSD never defined in Methods |
| logarithmic qubit scaling | a **reparameterisation, not a compression**: `P ≈ 2N+6n`, i.e. 132 classical parameters producing 43 torsions |
| a VQE | the energy is never a qubit observable; it is computed classically from the reconstructed Cartesian structure, with COBYLA/SLSQP |
| CVaR | **absent entirely** |

**They independently corroborate our central negative.** Over 7.3 M structures they find
all-atom potential energy *negatively* correlated with accuracy, and their own tables show
snapshot mining reaching sub-2 Å in ~99% of replicas while their energy functions select one
in 6%, 4% and 0.3% of cases. They print no rank correlation, no confidence interval and no
p-value in fifty pages. Our certified global optimum result is strictly stronger and remains
unpublished.

**An unnoticed defect in their decoder**, worth recording because we have measured the hazard
independently: their hardware decoder is `θᵢ = 2πCᵢ − π` with `Cᵢ` a cumulative probability,
so the decoded torsion vector is **always monotonically sorted with total variation ≤ 2π**.
That biases it toward near-constant torsion profiles — regular secondary structure — which
predicts their chignolin-works / Trp-cage-fails asymmetry from the decoder alone. It is the
same trap as our measurement that a constant α-helix beats uniform random sampling by
0.457 Å.

**Novelty audit.** Taken: off-lattice and torsion-space quantum folding; binned φ/ψ VQE at
peptide length; CVaR-VQE on lattices; continuous log-scaled torsion encoding with an
all-atom backend. **Surviving:** the exact locality theorem, the no-free-parameter
spectrum-to-gradient-variance chain, the certified global optimum on an enumerable space,
in-loop CVaR in torsion space, the causal VQE-versus-classical control at matched budget, and
the protected statistical instrument.

**Repositioning, not redirection.** Nothing invalidates the direction. The contribution is
not a new quantum encoding — that ground is taken and worth nothing — but **the first
properly controlled measurement of whether quantum conformational optimisation helps, on a
certified landscape.**

## L. CORRECTIONS

Kept beside what replaced them, per the project's standing rule. Four of the six retracted
claims were the coordinator's own; one killed the sprint's most attractive result.

| retracted | replaced by |
|---|---|
| *(coordinator)* Real torsion predictors have coherent errors, so the restraint surface is mispriced | They are near-i.i.d.; all ten emitters have lag-1 error autocorrelation in [−0.092, +0.066]. The surface stands and is mildly conservative. |
| *(coordinator)* Torsion-error cost is N-terminally dominated | No privileged terminus, front/back ratio 0.886. Optimal superposition rotates the whole structure, so a hinge argument replaces the lever argument. |
| *(coordinator)* The structural objective is eight times better ordered | Roughly fourfold, and distogram-driven. The proposal was drawn half from the prior it scores; under uniform sampling that term collapses +0.355 → +0.046. |
| *(coordinator)* An eightfold in-decile gain bought stability not accuracy | In-decile ρ is **not monotone** in objective quality — it peaks at moderate signal and falls where the objective is much better. Both axes now required. |
| *(coordinator)* AMBER cannot recover the 0.157 Å either | Invalid inference: the refinement experiment ran on structures already on the ideal manifold, where AMBER has nothing to fix. The off-manifold case is a different experiment. |
| *(coordinator)* The 0.157 Å is the price of ideal geometry | It is the price of **re-expanding a contracted structure**. Coordinate averaging shrinks the backbone 16.0% (N–CA 1.238 vs 1.458 Å) and contraction predicts cost at ρ = −0.829. Any operation restoring real bonds pays it. |
| `leg_torsion` selects 2.954 Å and is the best native-free selector | The **worst arm tested** at n=126 (argmin 4.933, in-decile −0.073), and it poisons every combination. Two of nine targets carried 68% — the concentration analysis called it in advance. |
| Three torsions per chain are inert; dead block `log₂k` | **Four** are inert; both terminal residues invisible; dead block `2·log₂k`. 14 of 18 live at n=9, k=4. Found independently by two agents. |
| Term reweighting of Legacy helps (in-sample +0.209 ρ) | Leave-one-target-out: **−0.031**. |
| Legacy's energy–structure decoupling rides on the invisible terminal DOFs | Only 3.8% of variance. |
| Legacy and AMBER are complementary **and exploitable** | Premise holds (partialled error correlation +0.096); route refuted. Fusion is +0.032 ρ under a per-target *oracle* weight and loses 0.240 ρ to LOTO. Pareto discards 72% of the truly-best 1%. |
| *(ENER's own, caught by a runtime assertion)* an 86% bond contraction in the averaged backbone | A frame bug: `build_backbone_batch` returns each window in the builder's frame, not the window's. Corrected figure is **16.0%**. |

**A data-integrity correction, not a claim.** The cached all-atom subset in the enumerated
files is **40% oracle-conditioned** — 0.401 Å better than the space it was drawn from, not a
uniform sample. Only `amber_kind == 0` rows are usable. Caught mid-sprint and propagated to
the agent training on those files before it could enter a learned objective.

## M. RECOMMENDATION

**Do not build the architecture in section A.** It is documented because its decomposition is
informative, not because it should ship. It loses to the incumbent at every setting.

**Stop funding search-side work in this representation.** Two independent measurements —
budget saturation at evaluation 300 of 20,000, and a 0.68 Å selection gap surviving certified
global optimisation — say that the quantity a variational optimiser improves is not the
quantity that limits this problem. That applies to VQE, CVaR, ansatz design, encoding choice
and initialisation alike.

**The one thing worth funding is anything that reveals the per-target sign — and the sprint's
last experiment shows that direction is NOT closed.** The in-band ordering axis is learnable to
**0.986** within a target and its direction flips between targets, so a predictor of that
direction would convert an existing capability into deployable accuracy.

The axis is compactness (per-target skill correlates **+0.909** with the native radius of
gyration), and that radius has a **native-free predictor already sitting in the pipeline**: the
shipped distogram predicts every pair distance, and `Rg^2 = (1/(2N^2)) sum_ij d_ij^2` follows in
closed form. Measured on all 126 targets, length-residualised so length alone cannot explain it:

| native-free proxy for the native Rg | length-residualised r | CI95 |
|---|---|---|
| retrieval pool mean Rg | 0.244 | [+0.031, +0.444] |
| distogram-predicted Rg | 0.313 | [+0.110, +0.518] |
| incumbent emitted Rg | **0.365** | [+0.155, +0.577] |

**All three intervals exclude zero — the signal is real and free to compute.** But whether it
*converts* is **unresolved, and the direct measurement is discouraging.**

Rather than let a composed estimate stand, the chain was measured end to end on the 12 targets
where both halves exist (`s14/signchain.py`; length-residualised, permutation p-values, which
are the right test at n=12):

| arm | Pearson | perm p |
|---|---|---|
| **ORACLE ceiling — native Rg → per-target skill** | **+0.416** | **0.177** |
| distogram-predicted Rg → per-target skill | +0.257 | 0.419 |
| incumbent emitted Rg → per-target skill | +0.196 | 0.537 |
| retrieval pool mean Rg → per-target skill | −0.044 | 0.888 |

**The ORACLE ceiling does not reproduce at +0.909 — it measures +0.416, not significant, with a
CI spanning zero.** The likely explanation is that the two measurements use *different skill
definitions*: the flip diagnostic correlated the headline arm's per-target skill against a
native-Rg z-score, while this correlates the ceiling experiment's `cross_rho` against a
length-residualised native Rg. Neither is necessarily in error. But the consequence stands —
the ceiling for *this* chain is +0.416, and **the "roughly 61% sign accuracy" figure should not
be quoted.**

Every native-free proxy has a permutation p above 0.4. The direction is not refuted — n=12 has
almost no power, and the 126-target proxy-to-compactness link stands on its own — but **it is
not supported by the only end-to-end measurement available.**

**So the next sprint's first job is not to build. It is to resolve that discrepancy** on a
common skill definition and a larger enumerated set. Only then is
proxy → per-target sign → emitted RMSD worth running.

**The most promising concrete instance:** 54 targets have deposited chemical shifts, and a
*forward* shift predictor scores candidate structures against them using **no native
coordinates at all**. That is a per-target signal available at inference, which is exactly the
missing ingredient — and it is a *discrimination* channel, unlike the generation channel that
was priced and closed. It was proposed unprompted by the SHIFT workstream and inherits none of
that route's negatives.

**Adopt the tail-restricted statistic as standing practice.** Of the four proxies this sprint
used for objective quality, only it tracks the optimum. Bulk pairwise accuracy below a 0.25 Å
gap *is* beatable — the distogram reaches 0.654 there against the energies' 0.511 — and it
predicts nothing. Read every tail number against a **random-tail null at 0.524–0.527, never
against 0.500**.

**And coordinate-space aggregation deserves study it has never had.** It is worth 3.4× the
objective, and a *perfect* ranker emits 2.474 Å through a decile against 1.219 Å through a
top-100 — a 1.26 Å swing from the terminal operator alone, on identical rankings.

**Publish the negative.** A certified selection gap with an independently corroborated
mechanism, on a protected instrument, is a stronger contribution than another encoding — and
the corroborating group published the correlation without ever printing a coefficient.

---

## I. TRAINABILITY, CVaR CORRECTNESS AND ENCODINGS

Re-derived from scratch rather than assumed, per the brief's instruction not to trust prior
tests. Unit tests in `s14/test_vqe.py`.

### I.1 The CVaR value estimator is correct. Three defects surround it.

The recorded `baseline="tail"` gradient defect is **real and reproduces bit-for-bit**, and its
bias now has a **closed form** rather than a cosine. Two further defects are new:

**A sampled-CVaR bias at non-integer α·N.** `cvar_from_samples` averages the lowest
`ceil(α·N)` samples; the correct estimator splits the boundary atom. They coincide exactly
when α·N is an integer and differ otherwise, always **upward**:

| N | α | α·N | bias (sd of sampled energy) |
|---|---|---|---|
| 13 | 0.10 | 1.3 | **+0.134** |
| 8 | 0.30 | 2.4 | **+0.112** |
| 10 | 0.25 | 2.5 | +0.081 |
| 1000 | 0.0333 | 33.3 | +0.008 |

It decays as 1/N and is **not a bug to fix** — it is the estimator, and it is consistent. But
at α = 0.05 with 256 shots the effective tail is 13 samples and the bias is **the same order
as the differences the optimiser is trying to resolve**. Any CVaR value quoted at small α and
modest shots is biased upward and must never be compared against an exact CVaR.

**An exactly-zero-gradient regime, with an iff:**

> `dCVaR/dp ≡ 0` **⟺** `p(x*) ≥ α`, where `x*` is the stable argmin of E.

Zero counterexamples in 3,000 random cases. This is a *hard* zero — not a barren plateau, not
a small gradient — and no estimator, shot count or baseline recovers it. It is **convergence,
not failure** (`CVaR_α ≥ min E` always, with equality exactly in this regime). But at
α = 0.01 roughly **one random initialisation in four to one in twenty starts dead**, and the
fraction grows as the distribution concentrates during training.

### I.2 α is not a learning rate — and at small α, CVaR stops being CVaR

Gradients at different α evaluated at the *same* point are **not parallel**, refuting the
workstream's own hypothesis. More consequentially, at small α the CVaR objective **collapses
to an argmin-finder**, proved as an identity rather than observed as a coincidence. That is
the mechanism behind the previously recorded trap in which two distinct AMBER variants became
literally the same objective at α ≤ 0.25.

**The α sweep** (18 qubits, budget 20,480, five seeds; initialisation RMSD 4.117 on every row
— the VQE starts at a random point, not at the answer). At objective signal 0.30, certified
optimum 1.929 Å:

| α | objective gap | RMSD returned | entropy (bits) | diversity | P(RMSD < 2.5) |
|---|---|---|---|---|---|
| **1.0** | 0.0119 | **2.266** | 4.07 | 0.178 | 0.277 |
| 0.25 | 0.0157 | 2.456 | 6.61 | 0.244 | **0.371** |
| 0.1 | 0.0125 | 2.300 | 7.51 | 0.289 | 0.247 |
| 0.01 | 0.0151 | 2.459 | 10.26 | 0.364 | 0.135 |

**Plain expectation-value VQE (α = 1.0) returns the best structure.** Lower α buys entropy and
diversity monotonically and does not buy accuracy. The one thing small α does buy is
probability mass on good states at α = 0.25 — which matters only if the terminal operator
consumes a *set*, and is the single place CVaR could still earn its keep.

### I.3 Encodings: binary, for four independent reasons

Binary attains the information-theoretic qubit bound exactly at power-of-two k; needs **510
Pauli terms where one-hot needs 454,463**; is surjective, so it requires no penalty and its
spectrum carries no penalty-dependent degree of freedom; and has the highest gradient
variance of any encoding tested. Gray coding buys nothing — a clean positive on one target
was **retracted after replication across 54 cells** showed a coin flip.

A methodological result worth keeping: **a non-surjective encoding's Pauli spectrum measures
its constraint, not its objective.** That invalidates naive spectrum comparisons across
encodings.

**No encoding change can rescue a bad objective**, which is the actual binding constraint.

### I.4 The QNG refutation holds only at depth 1 — a correction to the brief

Fubini–Study metric by central differences on the exact statevector, machinery wholly
independent of `core/quantum.py`, n = 8:

| entangler | depth | mean g_ii | max off-diagonal | condition number |
|---|---|---|---|---|
| every pattern | 1 | 0.250000 | 0.000000 | **1.000** |
| chain | 2 | 0.250000 | 0.142821 | 3.665 |
| ring | 2 | 0.250000 | 0.203157 | 9.937 |
| **block (torsion-aware)** | **2** | 0.250000 | 0.246720 | **155.2** |
| **block** | **3** | 0.250000 | 0.246720 | **478.9** |
| all-to-all | 2 | 0.250000 | 0.106178 | 2.478 |

The recorded claim — "QNG is refuted; the metric is exactly `I/4`" — is **true and complete at
depth 1, and true there for every entangler pattern. It does not extend to depth ≥ 2.**
`g_ii = 0.250000` exactly at every depth and pattern, but the off-diagonal structure appears
at depth 2 and the condition number reaches 155 for the torsion-aware block entangler and 479
at depth 3. **A torsion-aware ansatz at depth ≥ 2 does reshape the manifold**, so QNG is open
again in exactly the regime a problem-inspired ansatz would occupy. It remains a moot point
for accuracy while discrimination binds, but the blanket refutation was wrong.

## J. FAILURE-CLASS ANALYSIS

**The pattern that has now appeared four independent ways.** Four low-information arms beat
the incumbent on the 18 hardest targets: the class prior argmax (5.670), prior sampling with
coordinate consensus (5.771), the raw coordinate average (5.832) and a **zero-information
constant α-helix** (5.887), against the incumbent's 6.026. On the failure class, this system's
information is worse than none — Sprint 12's finding, reached here by three further routes.

**The failure class is not where the shift channel fails.** Availability of heteronuclear
chemical shifts is *higher* on the hard targets than elsewhere (fibril/lasso 0.500, FAIL18
0.444, other 0.398), because solid-state NMR of fibrils requires ¹³C/¹⁵N labelling. So the one
experimental channel that exists is preferentially available exactly where the system fails —
which is the opposite of the usual pattern and would matter if the arithmetic in C14 did not
close the route regardless.

**Per-target heterogeneity dominates every aggregate.** On the enumerated nine, the certified
optimum of the best objective spans **0.675 Å (2MK7) to 5.256 Å (1CS9, worse than its own
random draw)**, and two of nine targets **anti-rank**. Reweighting the objective is a
per-target coin flip rather than a smooth trade: w = 1.0 rescues 6EY3 from 3.820 to 1.948 and
destroys 2MK7 from 0.675 to 4.446. No aggregate in this dossier should be read without its
spread.

---

## STILL RUNNING AT WRITE-UP, AND WHAT THEY CAN AND CANNOT CHANGE

Two jobs were still in flight when this dossier was written. Both are **confirmations of
results already established at smaller n**, not open questions, and both write incrementally
so whatever they reach is readable. Stated here rather than omitted, per the sprint's own rule
that a partial result is reported with its n.

**AMBER relaxation versus geometric projection — complete at n=126 for k=10, and the restraint
frontier measured at n=50.** Reported in section F.2. The one live thread is **k=30 on all 126
targets**, which decides whether the sprint's only positive accuracy result survives at full
sample. Everything else is closed.

**Matched-budget arms on the structural objective against the certified optimum**
(`s14/results/vqe_signal_hamil.json`, at 6/9 targets; the blend grid at 5/9). This is the
formal closing experiment for the causality question. Its answer is already determined from
three directions that agree: the selection-gap reproduction across disjoint target sets, the
0/12 initialisation result, and the objective-axis result where VQE reaches the certified
optimum in 0–12% of cells against greedy's 71%. A contrary result here would be a genuine
surprise and should be treated as one.

**Not run, and honestly labelled.** Five gaps, none of them silent:

1. **The ansatz-family comparison** (`s14/vqe_ansatz.py`) — end-to-end comparison across
   entangler families, barren-plateau slopes, and a multi-seed replication of the Sprint 13
   expressivity table. **Stopped deliberately** under the CPU squeeze because it writes its
   JSON only on completion and would not have finished. Its two most valuable results were
   captured standalone and *are* in the findings: the depth-1-only scope of the QNG refutation,
   and the initialisation control that produced the 0/12 result. This is the one substantive
   experiment the sprint lost.
2. **Whether QNG *helps* at depth ≥ 2** is a HYPOTHESIS, not a finding — only the metric's
   condition number was measured, never an optimisation using it.
3. **Shots sensitivity** — untested.
4. **A QAOA arm** — absent.
5. **The out-of-distribution generator-shift test** (does a learned objective survive a VQE's
   own proposal distribution) — written, not run; it does not change the verdict.

And the **n=11 enumeration batch** was launched and abandoned under CPU contention; the
declining learning curve retrospectively justifies that, since more targets were making the
objective worse.

All five are recorded in the owning workstreams' own "what I did not establish" sections as
well as here, so no gap depends on this dossier to be visible.

## COMPUTE, HONESTLY REPORTED

The box was saturated for the whole sprint: CPU pinned at 100% with 4–10 concurrent Python
processes, RAM between 71% and 89% of 15.6 GB, against a brief targeting ~95% of both. RAM was
never the binding constraint; **cores were**. At the worst point individual agent processes were
measured at **3–5% of one core**, and one workstream reported ~0.07 cores per process.

Three consequences are recorded rather than hidden. The n=11 enumeration batch was abandoned.
One ansatz module was **deliberately stopped and replaced with a reduced form** because it wrote
its JSON only on completion and would have produced nothing. And the coordinator **killed its
own projected-consensus job at two hours** when a parallel workstream superseded it — holding a
core for a marginal duplicate while three agents owned the sprint's open questions was the wrong
trade.
