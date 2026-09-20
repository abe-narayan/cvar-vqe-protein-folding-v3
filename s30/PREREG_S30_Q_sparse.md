# PREREG S30-Q-1 — the sparse weighted readout: ceiling as a function of support size and of BITS

**Lane Q (quantum reformulation, leads L5 + L6). Registered before the first number exists.**
Written against `s29/s29_O_ladder.py` (verified present, 63,175 bytes) and the arms already on
disk in `s29/results/s29_O_cloud_rows.jsonl` (126 rows, supports and weights stored per arm) and
`s29/results/s29_O_chain_rows*.jsonl` (33 arms × 126 targets, built chain).

---

## 1. The question this answers, and why it comes before any Hamiltonian

S29 left exactly one ladder class unclosed by ceiling: a **sparse weighted readout**. Two members
with ORACLE weights emit **1.4315 Å** on the built chain where 75 members with ORACLE membership
emit **2.3055** and production emits **3.2105**.

The coordinator's framing, which I adopt: a low *parameter* count is not a low *information*
requirement. Choosing 2 of 500 is `log2 C(500,2) = 16.93` bits — **more** than the 7 bits the
deployed top-128 argmin needs, not less. And the sparse arms on disk fit **free continuous
weights per target against the native**, which is a second, unpriced information channel on top
of the support.

So the ORACLE number 1.4315 Å prices **support + continuous weights**, and nothing yet says which
of the two carries the gain, or what either costs in bits.

**This registration measures that decomposition before any circuit is designed.** If the gain lives
in the weights, a quantum *subset-selection* formulation cannot deliver it and the L5+L6
combination is dead on arrival. If it lives in the support, the bit price decides whether the
support can be chosen at all.

## 2. What I accept as binding, and what I attack (contract rule 24)

**Accepted as binding:**
- The terminal operator consumes the set mean (`operator-consumes-set-mean`, R² 0.89): a perfect
  rank-1 is worth ≈ −0.03 Å through the shipped m = 75 average.
- The deployed architecture's ORACLE ceiling is 2.9027 Å (S29 D1).
- Recognition is closed across bands, within bands and per target (S29).
- Per-target maxima over a family are order statistics (contract rule 8).

**Attacked:** S29's "a sparse weighted readout is the one class not closed by ceiling." I am
testing whether that class is closed by **information** rather than by ceiling — a different and,
if true, more final kind of closure.

## 3. The decomposition — a 2-way factorial on support and weights

Every arm is a convex combination `X = sum_j w_j * Wf[j]` over a support of size `s` drawn from a
candidate set `S`, in the S28 frame, scored by Cα RMSD to the native. ORACLE weights are solved by
`oracle_hull(..., support=...)` (the existing alternating convex-NNLS, unchanged).

| | ORACLE weights (continuous) | GRID weights (L levels) | UNIFORM weights (1/s) |
|---|---|---|---|
| **ORACLE support** (greedy, on disk) | A — *on disk* | Q(L) | B |
| **SCORE-PREFIX support** (first s of the DIS order) | C | — | D |
| **DIVERSITY support** (native-free farthest-point) | G | — | H |
| **RANDOM support** (8 seeded draws, MEAN not min) | E | — | F |

`s ∈ {2,3,5,10,20}`; `S ∈ {top75, pool(K=500)}`.

**To resolve the low-bit region** — where the comparison against the deployed argmin actually
happens — the same family is also run for `s ∈ {2,3,5}` over score-order prefixes
`T ∈ {4,8,16,32,64,128,256,500}` (`s < T`).

## 4. The bit ledger (the axis the curve is plotted against)

Both channels are counted the same way — as the log of the number of distinguishable choices:

- **support bits** = `log2 C(|S|, s)` for ORACLE and RANDOM support; **0** for SCORE-PREFIX and
  DIVERSITY support, which are deterministic functions of native-free quantities.
- **weight bits** = `log2 C(L+s-1, s-1)` for weights restricted to multiples of `1/L` summing to 1
  (the count of lattice points on the simplex); **0** for UNIFORM; **unbounded** for continuous
  ORACLE weights, which is why the continuous arms are reported as ceilings and never as routes.
- `L ∈ {1,2,3,4,8,16,32}`. `L = 1` puts all mass on one member, i.e. it degenerates to the argmin
  restricted to the support, and is carried as an internal consistency check.

**The reference curve at the same budget** is the deployed readout's own: the argmin over the top
`2^B` of the score order costs exactly `B` bits and delivers `best1_top(2^B)`. This is computed for
`B = 1..9` from the stored ORACLE per-member RMSDs.

## 5. Pre-registered falsifiers

**F1 — the support carries the gain.** At `s = 2` over the pool, ORACLE support with UNIFORM
weights (arm B) recovers **≥ 50%** of the gain `production − A`.
*Refuted if it recovers < 50%*, which would mean the sparse readout is not a subset-selection
problem and a subset Hamiltonian cannot deliver it.

**F2 — the weights are cheap in bits.** The grid arm Q(L) comes within **0.05 Å** of the continuous
ORACLE arm A at `L ≤ 4`.
*Refuted if it needs `L ≥ 16`*, which prices the weight channel at ≥ 4 bits per component and
pushes the total well past the argmin's budget.

**F3 — a native-free support rule has skill.** At `s = 2` over the pool, the better of SCORE-PREFIX
(C) and DIVERSITY (G) support, with ORACLE weights, beats the RANDOM-support ORACLE-weight **mean**
(E) by **≥ 0.20 Å**.
*Refuted otherwise.* Note the symmetric outcome is informative: if E is already close to A, the
support hardly matters and F1 is refuted by a second route.

**F4 — THE DECISIVE ONE. The sparse family beats the argmin at EQUAL bits.** For at least one
bit budget `B` in 4..12, some fully-priced sparse arm (support bits + weight bits ≤ B, no
continuous weights, no ORACLE support unless its bits are charged) achieves a **lower** mean RMSD
than `best1_top(2^B)`.
*Refuted if the argmin dominates the sparse family at every budget* — in which case the sparse
readout is not a route at any information price, the coordinator's suspicion is confirmed, and
**the L5+L6 combination is dead**. I will report that as the answer, not look for a fourth arm.

## 6. Order-statistic discipline (contract rule 8)

- The RANDOM-support arms are reported as the **mean over 8 seeded draws**, never the per-target
  minimum. The per-target minimum over draws is computed *only* to price the order statistic with
  `best_of_k_within` and is labelled as such.
- The ORACLE-support greedy path is itself a per-target selection over C(|S|,s); its bits are
  charged in the ledger, which is the correct pricing for this family.
- Any per-target maximum over `s`, `L` or `T` is reported with `best_of_k_within` and split-half
  transfer, or not reported.

## 7. Basis, labels, controls

- **Point cloud first** (fast, 126 targets, ~3 s each); the built chain is the endpoint
  (contract rule 1) and the decisive arms are chained afterwards. Every point-cloud figure says
  "point cloud" in the same sentence.
- The ORACLE/ORACLE corner (arm A) is **already chained** at every `s` for both sets, so the chain
  anchor of this family exists and is not re-derived.
- **Every arm in this registration except SCORE-PREFIX-UNIFORM and DIVERSITY-UNIFORM reads the
  native and is labelled ORACLE.** No deployable parameter is tuned here; this is a ceiling
  measurement and is reported as one.
- Greedy support selection is an approximation to the exhaustive `C(T,s)` optimum. One **exhaustive
  `s = 2` check at `T = 64`** (C(64,2) = 2016 pairs) bounds the greedy gap; if the gap exceeds
  0.05 Å the greedy arms are relabelled as lower bounds on the ceiling.

## 8. What would make me stop rather than continue

If **F4 is refuted**, I do not design a Hamiltonian for this family. The finding is reported as the
closure of the last open ladder class and the lead is returned. The charter's own instruction is
that we should both know this in an hour rather than a week.

---

*Registered by lane Q. No number in this document; none existed when it was committed.*
