# WORKSTREAM C — STATUS  ***FINAL — every block landed at its full pre-registered configuration***

*(This file was written as an interim status on request. Sections 1-4 below are preserved as
written; the FINAL RESULT block at the top is appended after all blocks completed. Nothing in the
interim text has been softened.)*

## FINAL: ALL SIX BLOCKS COMPLETE, NOTHING PARTIAL, NOTHING NOT MEASURED FOR WANT OF COMPUTE

| block | n | status |
|---|---|---|
| gates GC21a/c/d | 3-9000 | PASS; **GC21b mis-specified by me and restated, firing count reported** |
| N normalisation | 30 | COMPLETE |
| C1 mechanism | 30 | COMPLETE |
| C2 lambda sweep (gradient / variance / tail / Hessian) | 30 | COMPLETE |
| C2t selection-side matrix, declared extension | **126** | COMPLETE |
| C3 + P staged schedules and preconditioners | 30 | COMPLETE |

**HEADLINE: F-C3 FIRES.** `LA_Nt - A_raw = +0.2331 [+0.0673, +0.3950]`, 12W/18L, and it **survives
the move-size correction** (`+0.2285 [+0.1140, +0.3647]`). Staged Legacy->AMBER continuation is
**worse** than direct AMBER, not merely null -- exactly as C1's mechanism result predicted.
**F-P fires on every preconditioning arm**: direct AMBER reaches the lowest and only physical final
energy. The one arm that beat the comparator on RMSD (trust region, radius/3, -0.2987 [-0.6307,
-0.0547], 20W/10L) is **MOVE-SIZE, NOT PHYSICS** -- against its own matched null it is
`-0.0207 [-0.3185, +0.1736]`, 15W/15L.

**And the selection side at n = 126: all 60 arms of the mandatory matrix's physics half are WORSE
than a matched-count random tail, 60/60 CIs excluding zero.** Sprint 20's argmin result **does**
transfer to the tail operator. Full detail in `s21/agentC_FINDINGS.md`.

---

# (interim status as originally written)

## 1. WHAT HAS LANDED, WHAT IS RUNNING, WHAT IS QUEUED

| deliverable | status | artefact | ETA |
|---|---|---|---|
| **Normalisation** (assigned axis 3) | **LANDED, COMPLETE 30/30** | `s21/results/c_norm.json`, `c_norm_report.txt` | — |
| **C1 — the continuation's MECHANISM** | **LANDED, COMPLETE 30/30** | `c_c1.json`, `c_c1_report.txt` | — |
| **C3 — staged schedules + Block P preconditioners** | **RUNNING** | `c_c3.json` | ~60 min |
| **C2 — the λ landscape sweep** | queued | `c_c2.json` | ~60 min after C3 |
| **C2t — selection-side λ at n = 126** | queued | `c_c2tail.json` | ~20 min after C2 |
| gates GC21a/c/d | **PASS**, firings reported | `c_gate.json` | — |

**Nothing will be reported COMPLETE that is not.** `_write` sets `complete` only when the row count
covers the **full** 30-target pre-registered subset **and** the subset actually requested is that
subset — a smoke run cannot write a COMPLETE flag. If C2 or C2t does not land I will name the
missing configuration and report those cells **NOT MEASURED**.

---

## 2. THE TWO RESULTS THAT ARE ALREADY DECIDED

### (a) `λ` in raw units is not a dial — the priority-one experiment's parameterisation is a units artefact

Measured at 150 starts (`c_norm.json`), the crossover at which AMBER's term first dominates the
mixed gradient,

    λ* = ‖∇E_Legacy‖ / (‖∇E_Legacy‖ + ‖∇E_AMBER‖)
    median 1.99e-05     min 1.51e-09     max 5.58e-02     spread 3.7e+07 x

`H(λ) = (1−λ)E_L + λE_A` in raw units is **already pure AMBER at λ = 10⁻⁴**, and the λ at which it
becomes so varies **seven orders of magnitude** across targets. A uniform raw λ grid is a **step
function whose step position is a property of the target's steric state**. Any "abrupt transition"
or "λ-dependent basin hopping" read off a raw grid would have been reported as physics; it is
kcal/mol. **Label: EXACT for the statement, ESTABLISHED for the magnitudes.** This is why Block N
ran before the sweep.

### (b) C1 — the directive's stated MECHANISM for the continuation is REFUTED, and the reason is Legacy's own character

n = 30, paired, fold-clustered CI, built-chain basis, AMBER = bare single point.

    log10 E_AMBER        θ_L − θ_0   +1.4947 [+1.1437, +2.0174]   worse on 25/30
    log10 ‖∇E_AMBER‖     θ_L − θ_0   +1.4550 [+1.0944, +1.9101]   worse on 24/30

> **A Legacy minimisation does not carry the state OUT of the steric singularity. It drives it
> ~31× DEEPER in energy and ~28× deeper in gradient norm**, with CIs excluding zero.
>
> **And the mechanism is Legacy's own established character, not an accident.** `s20` L8:
> Legacy-preferred candidates are **0.45 Å more compact, 124W/2L**. Compaction is precisely what
> closes heavy-atom contacts, and closed contacts are what the steric singularity *is*. **Legacy is
> the wrong preconditioner for AMBER for exactly the reason it is Legacy.**

**My own pre-registered primary endpoint was the wrong endpoint and its own matched control says
so.** I registered *participation ratio* as primary. It rises — `+0.0132 [+0.0019, +0.0226]`, up on
19/30 — so **F-C1 as literally written DOES NOT FIRE**. But against the matched-magnitude
realisable null it is `+0.0089 [−0.0039, +0.0197]`, **NOT MEASURED**: a generic move of the same
torus size does the same thing. I am recording this as *"the registered falsifier did not fire, and
the hypothesis is refuted anyway on the two direct energetic reads"* rather than claiming support
from the clause that survived. The pre-registration stays unedited.

---

## 3. ON YOUR THREE FINDINGS — how they change what I can claim

**1. Readout vs training, and the +0.056 [−0.020, +0.129] training-H null.** Accepted, and it
sharpens my lane rather than ending it. `H(λ)` is a **training-side** object, so the continuation
can only be claimed as an **optimisation-path** effect. My pre-registration already derives
(`PREREG_C.md` §4, EXACT) that a strictly increasing transform of a single Hamiltonian leaves
`argmin`, every CVaR tail **set**, and every critical point invariant — **so my Block-P
preconditioners cannot change any endpoint objective by construction, only the path.** I will quote
your training-H null beside every C3 row and will not claim an endpoint effect.
**And C1 already says the path goes the wrong way**, which is a stronger statement than the null.
Your requested control — readout held fixed while only the training H moves along λ — is what C3's
arms are: every arm's readout is the same (built chain from the minimiser, mean over starts).

**2. Readout must be fixed before an objective statement is well-posed.** Agreed, and the two
readouts in my lane are already separated and will be labelled on every row:
* C1/C3 — **built chain**, one structure per start, mean over starts (a *best*-flavoured readout);
* C2/C2t — **point cloud**, the coordinate average of the CVaR α-tail (a *set-mean* readout).
**These are not compared to each other anywhere**, and I will not merge them.

**3. Normalisation as a governed fork, `BRIEF` §7 rule 0.** This is the pattern I ran: `Nt` (asinh
on the pool median/MAD) declared **PRIMARY in a git-frozen pre-registration (`7d07d63`) before any
Sprint-21 RMSD existed**, with `raw`, `Nz` (robust z) and `Ng` (gradient-matched) computed
**alongside on every arm**, and `Nr` (tailprice's rank→normal) **declared and rejected for
continuation with its reason given in advance** — zero gradient a.e., undefined off-pool. The
pre-registration also carries the commitment you are asking for: *if the primary conclusion flips
between the normalisations, it is reported NORMALISATION-DEPENDENT and no arm is promoted*, and a
non-declared winner is recorded as **the declared choice being WRONG**, not swapped in.
**I will add the explicit rule-0 fork enumeration — functional / basis / readout / normalisation /
null, with the alternative not taken named — to the `c_norm.py` and `c_cont.py` docstrings.**

Worth flagging back: **an affine normalisation does not fix this problem.** On the worst real
structure in the set (`2MIG`, +5.458e+23 kcal/mol) the four forks map it to
**5.5e23 / 1.6e19 / 1.6e12 / 44.89** for `raw` / `Nz` / `Ng` / `Nt`. A z-score fixes the units and
leaves the tail nineteen orders wide. If any other lane is z-scoring AMBER before combining it,
that arm is a clash census.

---

## 4. YOUR OPPORTUNITY — I am taking the cheap half of it, and the continuation is NOT blocked

The continuation is running, not stalled, so per your instruction I am **not** substituting the
regime-separation question for it. But **the Legacy–AMBER disagreement signal is nearly free on my
instrument**: `c_c2tail` already evaluates genuine Legacy and genuine AMBER over each target's own
75-member pool at **n = 126**, so per-target `ρ(E_L, E_A)`, the Jaccard overlap of their α-tails and
their tail-set energy margins cost one extra pass over data I am computing anyway. I will add them
as a **declared extension** and report whether they separate per-target selector skill.

**One caveat you should price before counting on it:** my pool is the shipped **top-75 retrieved
windows**, not the objective's **top-512**. A disagreement signal that works on mine is a
*candidate* for yours, not an answer to it, and porting it is a real experiment rather than a
lookup. I will say so in the write-up rather than let the number travel.
