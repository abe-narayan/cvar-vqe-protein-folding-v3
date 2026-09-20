# Sprint 29 — CVaR-VQE Protein Folding: the path below 2.5 Å

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await runs still in flight
(lane B's 126-target endpoint, lane M's F1 and F2, lane D's B3-at-n and the heavy test suite,
lane P's falsifier arms, lane X's 12-target configuration probes, lane T's 126-target compactness).
The charter requires the report only after everything is finished; this file is assembled as
results land so that nothing is reconstructed from memory at the end.

Branch `s26` · benchmark: 126 dev targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Ledger: `s29/LEDGER.md` (S29-L0 … S29-L46, no L43 — see the numbering note)
Running state: `s29/STATE.md` · Contract: `s29/S29_CONTRACT.md` · Charter: `s29/BRIEF.md`

---

## 0. The answer, up front

The charter asked for the path below 2.5 Å and named 3.0 Å as a meaningful intermediate.
**Neither is reachable through the current architecture, and this sprint established that as a
measured bound rather than as a failure to find a trick.** Production stands where it stood:

| | point cloud | **built chain** |
|---|---|---|
| production (deployable incumbent, n = 126) | 3.0483 | **3.2105** |

The sprint's contribution is not an Ångström. It is that the gap is now **located and priced**, in
three numbers that did not exist before, all on the charter's own endpoint:

1. **The achievable native-free bound (S29-L23, attacked and survived in S29-L35).**
   Every native-free operator this project can build is a displacement of the production cloud, and
   its entire value is one cosine: `RMSD = RMSD_prod · √(1 − ρ²)`. Reaching 3.00 Å needs ρ = 0.358;
   2.50 Å needs ρ = 0.628. Across **39 displacement fields** — 21 in the pre-registered attack plus
   18 more — **not one beats the random-shape reference of 0.1398**, and the best field, stepped by
   an amount chosen *with the native in hand*, is worth **0.019 Å**. The bound is ≈ 3.18–3.21 Å.

2. **The architectural ceiling (S29-L30, corrected to the endpoint in S29-L44).**
   The deployed quantum stage sees only the top-128 prefix and its tail *is* a prefix
   (Barkoutsos eq 12). Its ORACLE ceiling — the best it could emit **with the native in hand** —
   is **2.9027 Å built chain** (n = 126, vs production 3.2105, −0.3079, MDE 0.0982, **3.14× MDE**).
   The charter's 2.5 Å is therefore unreachable through this architecture *even with an oracle*.

3. **The readout, not the search and not the pool, is what binds (S29-L44).**
   Holding the candidate set fixed at the top-128 the quantum stage actually sees:

   | | built chain, n = 126 |
   |---|---|
   | prefix-average readout (what ships), ORACLE per-target *m* | **2.9027** |
   | single best member of the same 128 (order statistic, **7 bits**) | **2.1435** |
   | free convex combination of the same 128 (*expressiveness only*) | 1.8538 |
   | free convex combination of the full K=500 pool (*expressiveness only*) | 1.1235 |

   Choosing *m* and choosing *which member* cost the **same 7 bits** and are worth −0.3079 and
   −1.0657 against production. Same information budget, **3.5× the payoff**, differing only in what
   it is spent on. The architecture spends its information on the wrong question.

And the reason those 7 bits cannot simply be supplied: **the per-target sign is an incidental
parameter** in the sense of Neyman & Scott (1948) — not estimable from other targets' answers as a
matter of statistical theory (S29-L31). Recognition was then closed three independent ways
(§9). That is the sprint's unified finding, and it has a theorem rather than an anecdote.

> **A caveat that is load-bearing, stated here and repeated wherever the numbers appear.**
> The 0.7592 Å in row 2 of that table prices switching the terminal operator from an average to a
> selection **and** supplying the bits — *not* a ranking improvement fed to the operator that
> ships. The terminal operator consumes the set **mean**: `d_out = 1.16·d_set_mean + 0.04·d_set_best`
> (R² = 0.89). A perfect rank-1 signal is worth ≈ −1.74 Å through argmin and ≈ **−0.03 Å** through
> the m = 75 average. Lane M caught me stating this wrongly; the correction is S29-L44 addendum 2.
> The hull rows are **expressiveness, not ceilings** — they fit 128 or 500 free weights per target
> against the native, and grid oracles are order statistics.

---

## 1. What was tested

[PENDING — completed once the last lanes report. Structure: the eight lanes, their briefs, the
pre-registered falsifier for each experiment, and the multiplicity ledger.]

Eight lanes ran under `s29/S29_CONTRACT.md`, with three permanent roles the charter required:
adversary (D), divergent (X), literature (L).

| lane | brief | deliverables |
|---|---|---|
| **L** | literature, permanent | `s29/lit/L_1..L_8*.md`, `L_INDEX.md` (~2,600 lines), `s29_L_FINDINGS.md` (17 findings) |
| **M** | data path, convenience choices, harness audit | `DATAPATH.md`, `CONVENIENCE_CHOICES.md` (31 entries, 17 untested), `s29_M_harness_audit.md` |
| **T** | theory, permanent | `THEORY.md` (1,080 lines), `THEORY_SUMMARY.md` |
| **O** | ORACLE ceiling ladder | `s29_O_ladder_table.json`, `s29_O_chain_rows*.jsonl` |
| **D** | adversary, permanent; the cost–RMSD meter | `s29_D_cost_audit.py` and its six numbers |
| **X** | divergent, permanent: configuration-space encoding | `s29_X_probe*.json` (94 arms) |
| **P** | the projection stage's price | `s29_P_*.json`, `s29_P_rows_shard*.jsonl` |
| **B** | compatibility Hamiltonian, tail-then-aggregate | `s29_B_tta.py`, `s29_B_tta_*` rows |

---

## 2. What changed

**Nothing deployable changed, and that is a result rather than an absence of one.** No parameter,
no scorer, no stage of the shipped pipeline was altered, because no candidate change cleared its
pre-registered falsifier. The charter's rule — *below 0.7× MDE is not a result* — was applied to
this sprint's own positives as strictly as to anything else, and it disqualified several.

What changed is the **map**:

- The gap is no longer "somewhere in the pipeline." It is in the readout and in 7 bits per target,
  and both halves are now measured on the endpoint.
- Seven candidate routes were closed **without spending a full endpoint run**, five of them by
  derivation or by a measurement that already existed in the record (§14).
- One number that had been quoted throughout the sprint — the architectural ceiling — turned out to
  be on the wrong axis, and was corrected (§12, and the record of my own errors in §13).

---

## 3. What failed

[PENDING]

## 4. What succeeded

[PENDING]

## 5. Why, at the level of mechanism

[PENDING]

## 6. What the VQE contributed, with controls

The charter's hard constraint was that CVaR-VQE remain the central scientific object. It did, and
the sprint's clearest *positive* measurement is about it.

**The set-equality theorem fails on real pools (S29-L25, taken to the full instrument as S29-L45).**
CVaR is defined on *sorted* samples for diagonal Hamiltonians (Barkoutsos et al. eq 12), so the
deployed tail is provably a prefix of the energy order. Placing the objective on the tail's own
coordinate average, `V(S) = f(mean of S)`, breaks that. Exhaustively over all C(500,2) = 124,750
pairs of each of the **126** pools, with a tie rule corrected mid-sprint (below):

```
the f-optimal PAIR is not the energy prefix       114 / 126   (90.5%)
the f-optimal m = 5 SUBSET is not a prefix        124 / 126   (98.4%)
the per-state sort finds the exhaustive optimum    12 / 126   ( 9.5%)
mean objective gap over the prefix, m = 5             +0.1499
```

The optimum is reachable **neither by sorting E nor by sorting f**. This is the first measurement in
the project's record of *which set* as a genuine optimisation variable rather than the read-out of a
sort, and it answers charter §11 questions 5 and 7 affirmatively.

**And it buys nothing measurable.** The price, decomposed (ORACLE, point cloud, n = 126):

```
f-optimal m=5 subset − production   +0.2451   1.45× MDE   5/5 folds   power 0.98   WORSE
  of which m = 75 → m = 5 alone     +0.1647   1.36× MDE   5/5 folds   WORSE
  of which the NON-PREFIX CHOICE    +0.0804   0.73× MDE   4/5 folds   NOT MEASURED
```

Two thirds of the aggregate harm is simply that a 5-member average beats nothing — an m-ladder the
project has priced three times. **The part attributable to the hypothesis itself is inside its own
MDE.** Lane B reported it that way rather than banking the aggregate as a refutation, and read lane
D's field-survey caveat first: MSET_5 has signed cosine +0.063 with the direction to the native, so
a harmful result here is the *expected* value of an unsigned displacement.

This bounds the endpoint arm: a classical exhaustive/greedy search over the same objective is a
strict upper bound on what the CVaR-VQE could find under it, and that search does not beat
production.

**Controls.** [PENDING — lane X's 94-arm ladder, including the untrained-circuit and best-of-N
controls, and lane B's endpoint arm.]

**A defect found and fixed in our own code, mid-result.** `s29_B_tta.py::subset_target` carried the
comment *"average over the tied argmin set rather than reading array order"* and then took `tie[0]`
— array order, which on a DIS-sorted pool is the best-ranked member, biasing the rule **toward the
prefix**, the very hypothesis under test. This is the project's own named failure mode
(`tie-breaking-leaks-the-pool-order`, where an argmin on a tied signal once invented a 1.386 Å
winner). Fixed in `7b2e83e1`: one global argmin over all pairs, ties broken by a stable per-target
RNG, tie-set size recorded on every row, and a regression test whose name states the invariant. The
maximum tie set at n = 126 is 4, so the fix was not cosmetic — though it changed no number on the
original 12, which lane D verified independently (S29-L38). The pre-fix rows are kept, not deleted.

## 7. The final built-chain RMSD, with full statistics

[PENDING — lane B's 126-target endpoint is running.]

## 8. Uncertainty, MDE, fold CI, concentration

[PENDING]

## 9. Every control, and what it ruled out

[PENDING — partial: recognition closed three ways.]

## 10. The literature relied on, and rejected

[PENDING — lane L's 17 findings and the five-paper reading list for S30.]

## 11. The best architecture, and why

[PENDING]

## 12. What remains unresolved

[PENDING]

## 13. The next question, and the evidence that would settle it

[PENDING]

## 14. Every hypothesis entertained and killed, with the reason

[PENDING — the full table. Seven routes closed without an endpoint run; several more by measurement.]
