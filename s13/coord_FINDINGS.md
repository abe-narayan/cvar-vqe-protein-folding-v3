# Sprint 13 — COORDINATOR findings

My own experiments: the representation ceiling, the cost model, objective validity, the shape
of the AMBER objective, and the native-free search that everything else is upstream of.
Agent findings are in `s13/{tors,qarch,geo,walsh,lit}_FINDINGS.md`.

Every arm that reads a native quantity is labelled ORACLE or DIAGNOSTIC. benchmark60 was
never read; dev24 was never run; no tracked file was modified.

---

## C0. The representation is NOT the barrier — `s13/ceiling.py`

**Status: PROVEN (ORACLE DIAGNOSTIC — it reads native torsions to score).**

The proposed architecture optimises over a discrete, sequence-conditioned, leakage-safe
torsion library: `torsion_lib2.library_for(seq, k, exclude_seq=seq)` gives k states per
residue, so a peptide is `n_res · log2(k)` qubits. Before building anything on that space I
measured what it can express. Per residue the native torsions are snapped to their nearest
library state (`snap`, chain-blind), then coordinate descent on true CA-RMSD refines the
assignment (`descent`, chain-aware).

    ideal-geometry rebuild of the NATIVE continuous torsions: 0.347 A
    (the floor any torsion representation inherits)

| k | mean qubits | snap | **descent** | <2 Å | <1.5 Å | FAIL18 | other-108 |
|---|---|---|---|---|---|---|---|
| 4 | 25.9 | 2.860 | **1.594** | 0.73 | 0.48 | 2.078 | 1.513 |
| 8 | 38.9 | 2.218 | **1.184** | 0.91 | 0.70 | 1.539 | 1.125 |
| 16 | 51.8 | 1.844 | **0.876** | 0.95 | 0.89 | 1.296 | 0.806 |
| 32 | 64.8 | 1.375 | **0.634** | 1.00 | 0.97 | 0.783 | 0.609 |

**At k=4, in ~26 qubits, the space already contains a 1.594 Å answer** — against the shipped
retrieval pipeline's 3.213 Å and its 6.019 Å on the failure class. So the sprint's question
reduces cleanly to whether a *native-free objective* can find those configurations.

The `snap`→`descent` gap of 1.27 Å at k=4 is itself informative: the best per-residue choice
is not the best chain, because torsion error compounds. That is genuine interaction structure,
in contrast with Sprint 12's assembly Hamiltonian where the interaction was 2.4% of variance.

## C1. Cost model — and a correction to my own first measurement

**Status: PROVEN.** My first timing said AMBER single-point cost 6 ms. **That was wrong: it
re-evaluated the same state and hit `core.amber`'s result memo.** Measured on DISTINCT
configurations with one warm builder at threads=1:

| call | cost |
|---|---|
| `core.energy.components_batch` (Legacy) | **0.3 ms** |
| `core.amber.single_point`, distinct states | **28 ms** |
| `core.amber.single_point`, repeated state (memo hit) | 0.8 ms |
| `core.amber.refine_coords(k=10, steps=0)` | **9.1 s** |

**AMBER is ~90× Legacy per objective evaluation** — affordable in a VQE loop at ~9 minutes
per target for 19,200 evaluations, and the reason budget parity must be counted in
evaluations and never in wall time. The error is recorded because it was live in the shared
brief for an hour and agents were planning feasibility against it.

**An operational conflict worth recording:** `core.amber.builder_for` refuses to start an
OpenMM context above **92% physical memory**, while the sprint's compute target is ~95%.
Those cannot both be satisfied, so AMBER work has to be batched into windows where the rest
of the fleet leaves headroom. Catching the `MemoryError` and continuing would silently turn
a stalled column into "AMBER has no signal" — a fabricated negative.

## C2. Objective validity: Legacy separates the native from garbage; AMBER does not rank at all

**Status: STRONGLY SUPPORTED (DIAGNOSTIC).** `s13/coord_objval.py`, `s13/coord_ambershape.py`.
Per target, configurations drawn uniformly from the k=4 space, plus the native snapped to its
nearest states and the coordinate-descent optimum, scored by every candidate objective
against true CA-RMSD.

On 1A13 with 40 uniform samples:

| objective | ρ(E, CA-RMSD) | log₁₀ dynamic range | native percentile | descent percentile |
|---|---|---|---|---|
| Legacy (11-term) | **+0.398** | 2.38 | 0.000 | 0.000 |
| **AMBER raw single-point** | **−0.017** | **12.64** | 0.000 | **0.275** |
| AMBER, energy capped at p90 | −0.006 | 10.44 | 0.000 | 0.275 |
| AMBER, monotone log compression | −0.017 | **1.37** | 0.000 | 0.275 |
| AMBER, electrostatics + solvation only | **+0.245** | 2.42 | 0.050 | 0.150 |

Three readings, and the second is the one that matters:

1. **Raw AMBER single-point on an unrelaxed ideal-geometry build is a clash detector, not a
   ranker.** Twelve and a half orders of magnitude of dynamic range, ρ ≈ 0. It separates the
   native from random garbage trivially and then says nothing about which non-native
   configuration is better — and it scores the 1.888 Å descent structure worse than 27% of
   random ones.
2. **The damage is in the rank ORDER, not the scale.** The log compression is *monotone* by
   construction, so it cannot change the ordering, and indeed ρ is identical to raw
   (−0.017 both) while the dynamic range falls from 12.6 to 1.4 decades. **You cannot fix
   AMBER here by rescaling or by softening.** This control is the point of the arm.
3. **Deleting the steric term recovers ρ to +0.245.** The information in AMBER lives in the
   electrostatic and solvation terms and is masked by the `r⁻¹²` core.

### C2 COMPLETED, 20 targets x 120 configurations — and it CORRECTS the single-target reading above

`s13/results/coord_amber_shape_report.json`. `rho_low10` is the rank correlation INSIDE the
low-energy decile, which is where an optimiser actually lives and is the number that decides
usability.

| variant | rho over all | **rho in the low decile** | log10 range | native pct | descent pct | argmin RMSD |
|---|---|---|---|---|---|---|
| Legacy | +0.198 | **+0.043** | 2.14 | 0.319 | 0.308 | 4.120 |
| AMBER raw | -0.036 | -0.088 | **16.06** | 0.382 | 0.485 | 5.202 |
| AMBER capped at p90 | -0.036 | -0.088 | 10.68 | 0.375 | 0.480 | 5.202 |
| AMBER, monotone log compression | -0.036 | -0.088 | **1.81** | 0.382 | 0.485 | 5.202 |
| AMBER, electrostatics + solvation only | +0.015 | +0.027 | 2.42 | 0.399 | 0.483 | 4.781 |
| AMBER after restrained minimisation | +0.103 | — | 2.83 | — | — | 4.506 |

**My single-target reading was optimistic and is superseded.** On 1A13 I recorded the native
at percentile 0.000 under both energies and Legacy at rho +0.398. Across 20 targets the native
sits at the **32nd to 40th percentile** — roughly a third of random configurations score
better than the native — and **Legacy's low-decile rho is +0.043, i.e. nothing.** Legacy is a
garbage detector too, not only AMBER. The 1.6 A descent structure is scored worse than 31-48%
of random garbage by every variant.

What survives unchanged is the monotone control, and it is the point of the arm: the log
compression cannot alter the rank order, and indeed rho is identical to raw (-0.036, -0.088)
while the dynamic range falls from 16.1 to 1.8 decades. **The damage is in the ordering, not
the scale; no rescaling or softening fixes it.** Minimisation is the only variant that moves
rho at all (+0.103) and it costs 9.1 s per evaluation.

Two independent agents reproduced the native-percentile figure: the architecture agent
measured 32.5th percentile over 126 targets x 12,001 configurations, and on nine fully
enumerated targets found Legacy's **certified global optimum is +0.139 A worse than random
sampling**. Three routes, one conclusion: **neither energy ranks the native in torsion
space.**


---

## C3 (COMPLETED, all 126 targets). No native-free objective finds the answer the space contains

**Status: PROVEN.** `s13/coord_search.py` → `s13/results/coord_search_b5000.json`.
Simulated annealing over the k=4 torsion space (~26 qubits) at a **matched 5,000
objective-evaluation budget**, with uniform random sampling as the mandatory control.

| arm | mean | median | <2 Å | <1.5 Å | FAIL18 | other-108 |
|---|---|---|---|---|---|---|
| random sampling (argmin Legacy over 5,000 draws) | 4.522 | 4.414 | 0.01 | 0.00 | 5.690 | 4.328 |
| **SA on Legacy** | **4.624** | 4.433 | 0.02 | 0.00 | 5.835 | 4.422 |
| SA on the 1-local torsion prior | 3.969 | 3.948 | 0.24 | 0.15 | 5.670 | 3.686 |
| SA on prior + Legacy (standardised) | 3.980 | 4.030 | 0.25 | 0.17 | 5.763 | 3.683 |
| *ORACLE snap (native, per-residue nearest state)* | *2.860* | *2.297* | *0.42* | *0.32* | *4.262* | *2.627* |
| **ORACLE descent (the representation ceiling)** | **1.594** | 1.548 | 0.73 | 0.48 | 2.078 | 1.513 |
| *for reference: the shipped retrieval pipeline* | *3.213* | | | | *6.019* | *2.745* |

Paired against the random-sampling control at matched budget:

| arm | Δ vs random | W/L | drop-top-10 |
|---|---|---|---|
| SA on the prior | **−0.553 [−0.865, −0.238]** | 75/51 | −0.243 |
| SA on prior + Legacy | −0.542 [−0.862, −0.221] | 71/55 | −0.232 |
| **SA on Legacy** | **+0.102 [−0.054, +0.255]** | 58/67 | +0.243 |

**Two findings, and the second one kills the first.**

**(a) Optimising Legacy is worse than not optimising at all.** Simulated annealing on the
genuine 11-term energy returns structures 0.10 Å *worse* than picking the best of 5,000
uniform random draws, and it does so while driving the energy far below the native's (on
1A13: SA reaches −21.6 against the native's −16.0). The native is not the Legacy minimum and
the Legacy minimum is nowhere near the native. This is the project's oldest finding —
"Legacy in generation hurts" — reproduced on 126 targets in a new representation with the
mechanism in plain view, and it is exactly the pathology the objective-validity gate exists
to catch. **The gate worked: no VQE was built on this objective.**

**(b) The torsion prior's 0.553 Å advantage is entirely a helix artefact.** The prior beats
random, survives drop-top-10, and looks like the sprint's one positive — so I ran the
control. Native helix fraction against outcome (ORACLE diagnostic, native secondary structure
from `I.ss_of`):

| | n | SA on prior | random | ORACLE ceiling |
|---|---|---|---|---|
| helical targets (>50% helix) | 42 | **2.066** | 4.229 | 1.067 |
| non-helical targets (<10% helix) | 70 | **5.218** | 4.827 | 1.918 |

    rho(helix fraction, SA-prior RMSD)              = -0.744  (p < 1e-4)
    rho(helix fraction, prior's advantage over random) = -0.646  (p < 1e-4)

**On non-helical targets the prior is WORSE than random sampling.** The 1-local prior's modal
state is the α basin, so "optimise the prior" means "build an α-helix", and 42 of 126 targets
are substantially helical. It is a helix generator, not a predictor, and quoting its 3.969 Å
aggregate without this control would have been the same class of error the project's
corrections ledger already records four times.

**The complete ladder, and the gap it exposes.** The discrete torsion space at ~26 qubits
contains a 1.594 Å answer. The best native-free arm tested returns 3.969 Å and is an
artefact; the honest native-free arms return 4.5–4.6 Å; a leave-fold-out sequence-only
torsion predictor built by the torsion agent returns 4.151 Å on direct build. **Nothing
native-free reaches even the 3.213 Å of the retrieval pipeline it was meant to replace.**
The ~2.4 Å between the ceiling and the best honest arm is pure recognition failure — the same
wall Sprint 12 measured, now reproduced in a completely different representation, which makes
it a property of the problem rather than of the retrieval architecture.

## C5. What the sprint's two halves are now worth

The **folding half is a rigorous negative**: the torsion representation is excellent and
neither a learned sequence-only predictor nor any physical energy can locate good
configurations inside it. That is worth stating precisely because the representation was the
last untested structural hypothesis in the project.

The **trainability half stands on its own and is unaffected by that negative**, because it
does not depend on the folding working — it needs only that the two energy models are real,
that the encoding is real, and that the measurements are exact. On present evidence it
delivers: an exact locality theorem for distance terms in torsion space, the first Pauli-weight
spectrum measured for a genuine molecular force field, a closed no-free-parameter chain from
that spectrum to gradient variance, and an honest negative on whether cost-locality explains
trainability in the reachable regime.
