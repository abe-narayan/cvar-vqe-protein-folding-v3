# S29 LANE O -- FINDINGS

The ORACLE ceiling ladder by operator class, and the typicality-axis probe. Ledger entries
S29-L20 (rung 6), S29-L21 (rung 8), S29-L30 (rung 9), S29-L47 (the complete ladder). Pre-
registration `s29/PREREG_S29_O.md`. Code `s29/s29_O_ladder.py`, tests `tests/test_s29_O.py`
(13 pass). Basis: the BUILT CHAIN is primary; every point-cloud figure says so.
EVERY NUMBER IN THIS LANE IS ORACLE except the arms explicitly marked DEPLOYABLE.

## DEMONSTRATED

**D1. The deployed quantum architecture's ORACLE ceiling is 2.9027 A on the built chain,
against production's 3.2105 (-0.3079, 3.14x MDE, fold CI [-0.3625, -0.2499], 5/5, 117W/9L).**
The quantum stage sees the top-128 prefix (`core/pipeline.py:758`) and the set-equality theorem
makes the CVaR tail exactly a prefix of the energy order, consumed by the uniform average. The
best prefix-m average over the top-128 is therefore the best structure the architecture can
emit with the native in hand. **The charter's 2.5 A is unreachable through this architecture
even with perfect ORACLE selection.** `s29/results/s29_O_headline_contrasts.json`.

**D2. With the candidate set held fixed, switching the readout from the prefix average to one
member of the same top-128 is -0.7592 A (4.18x MDE, 5/5, 120W/6L).** It prices a READOUT SWITCH
plus 7 bits (log2 128), not a ranking improvement: the terminal operator consumes the set MEAN
(`operator-consumes-set-mean`, R^2 0.89), so a perfect rank-1 through the shipped m = 75 average
is worth about -0.03 A. Stated with that constraint in S29-L47 section 2.

**D3. Choosing m and choosing WHICH member cost the same 7 bits and are worth -0.3079 A and
-1.0670 A against production.** Same information budget, 3.5x the payoff. The architecture
spends its per-target information on the wrong question.

**D4. Which class could reach 2.5 A on the built chain, at all.** Ranking clears it at every
prefix (2.3055 / 2.1435 / 1.7078 for top-75 / top-128 / K=500); the weighted-average family
CONTAINS structures that clear it (2-member 2.1683 / 1.4315; hull 2.0646 / 1.8538 / 1.1235) --
expressiveness, not a ceiling, since those fit free weights per target against the native;
**the prefix-average class (the CVaR tail) and the basin-average class miss it at every
prefix** (2.9027 / 2.7763, and 2.6596 / 2.5620 at k = 8). The binding operation is WEIGHTING,
and the deployed readout has no weights.

**D5. The projection price is a monotone function of the cloud's own accuracy.** corr(cloud
RMSD, chain minus cloud) = +0.866 over 32 chained arms; the price is <= 0 on every arm under
2.31 A of cloud accuracy and above +0.10 on every arm above 2.41 A. S10-5 and S28-L26b's
pattern, now with a sample instead of two points, and it matches lane P's shape-distortion
curve: the arms that pay are the ones whose clouds carry the distortion.

**D6. Both of S10-5's published hull conventions reproduce.** `cf` (one shared transform solved
jointly with the weights, the emittable one) 1.9975 / 1.1167 against its 1.987 / 1.094; `oa`
(each candidate posed on the native individually, NOT emittable) 1.8016 / 0.9532 against its
1.802 / 0.953. The two are never mixed in the table.

## REFUTED (the lane's pre-registered kill)

**R1. H1, the typicality axis, is dead on both registered clauses (S29-L20).** The ORACLE
cosine between u = (shipped average minus sequence-blind average) and v = (native minus shipped
average), rigid-body removed, is **-0.0584** (SE 0.0324): indistinguishable from a random shape
field's SIGNED cosine (0.73x MDE) and WORSE than the random fields' |cos| 0.1443 at 2.19x MDE,
fold CI [-0.2355, -0.1663], 5/5 folds. The ORACLE best GLOBAL step is **exactly t = 0** -- one
scalar with the native in hand cannot beat doing nothing -- and the leave-fold-out step is
bit-identical to production on 126/126 targets on the built chain. The secondary blind
definition (S24's B') is the same picture at +0.0037 A on the chain (0.29x MDE).
Mechanism: rho(cos, t*) = +0.810, exactly as t* = (|v|/|u|) cos predicts, and the signs cancel
(54% want t < 0, 36% t > 0). On FAIL18 the axis points AWAY from the native three times harder
(-0.3171 vs -0.0153; random-18 null p = 0.0006), the OPPOSITE sign to the registered prior, and
the mechanism is in the record: on those targets the sequence-blind pipeline BEATS the shipped
one (`sequence-conditioning-hurts-the-failures`), so production is further from the native than
typical is and u points outward by construction.

**R2. Lane T's one-parameter PC1 family has an achievable ceiling of exactly zero (S29-L21).**
The ORACLE best GLOBAL eta along the pool's first shape mode is 0.0000 A; the ORACLE best eta is
positive on 52% of targets (the sign is a coin flip, as T derived); the two one-global-sign arms
are -0.2142 and -0.2455; the leave-fold-out eta is +0.0071 A WORSE than production (0.82x MDE).
T predicted "under 0.15 A" and the answer is the floor of that interval.

## THE UNIFYING MEASUREMENT

**U1. Every rung with a free scalar has an ORACLE GLOBAL value of 0.0 to 0.6% of its per-target
gain, and a leave-fold-out value on the wrong side of zero.**

| rung | free parameter | ORACLE per-target | ORACLE GLOBAL | LEAVE-FOLD-OUT |
|---|---|---|---|---|
| 2/9 | prefix size m inside the top-128 | -0.2879 | -0.0018 (0.6%) | +0.0079 |
| 2 | prefix size m over K=500 | -0.4421 | -0.0018 (0.4%) | +0.0079 |
| 6 | step t along the typicality axis | -0.3061 | **+0.0000 (exactly)** | +0.0000 |
| 8 | step eta along the pool's first shape mode | -0.4543 | **+0.0000 (exactly)** | +0.0071 |

With lane L's S29-L31 this is one fact, not four: the scalar is a Neyman-Scott INCIDENTAL
PARAMETER -- one nuisance parameter per target with a bounded number of observations per target
-- and it is not estimable from other targets' answers as a matter of statistical theory rather
than model capacity. The size of a rung's ORACLE gap therefore measures how much PER-TARGET
information it needs, not how much accuracy it offers. It also sharpens S22 L4: that entry
measured per-target m transferring 65% ACROSS POOL HALVES OF THE SAME TARGET (using that
target's own native); ACROSS TARGETS it transfers -3%.

## WHAT DAMAGED MY OWN EXPECTATIONS

1. **I expected the typicality axis to be weakly positive and merely unusable.** It is
   NEGATIVE, and the ORACLE global step is exactly zero. "Do nothing is the ORACLE optimum of a
   one-parameter family" is a far stronger statement than the falsifier I registered, and it
   happened twice (rungs 6 and 8) on two unrelated families.
2. **I expected the basin-average class to be the interesting one** (a CVaR-VQE that selects a
   coherent basin is the sprint's natural quantum story). It is the worst named class on the
   board: at k = 8 over K=500 it emits 2.5620 A, worse than two members with weights (1.4315)
   and worse than one member chosen well (1.7078).
3. **I expected the hull numbers to be the headline.** They are expressiveness, not ceilings --
   lane M's correction, adopted -- and the number that answers the charter is the dullest row in
   the table, the prefix average.
4. **I expected my convex solver to be the problem** when rung 5 came in 0.195 A above S10-5's
   published 1.802. It was a CONVENTION difference (S10-5 reports two and this lane initially
   compared against the wrong one); both anchors reproduce to 4 decimals once the `oa` arm was
   added. I nearly declared a solver defect that does not exist.

## WHAT I DID NOT DO AND WHY

- **No built chain for rungs 6's per-target step or rung 8** (the point cloud decides them: both
  achievable arms are at or on the wrong side of zero, and the coordinator's condition for
  chaining rung 8 was 0.15 A, which it does not approach).
- **No second seed**, because there is no positive to replicate: every deployable arm in this
  lane is NOT MEASURED and two are exact zeros.
- **No attempt to build a per-target predictor of m, t or eta.** Five router constructions
  (S22 L7, S23 L7) and three more in S28 already failed at it, and S29-L31 now says why.
- **No AMBER, no benchmark60, no pinned artefact touched.**

## OPEN

- The weighted-average family contains structures at 1.12 A on the built chain with 10 members.
  Nothing in this lane says an operator could find them, and lane M's expressiveness caveat
  applies; but the gap between "two members with ORACLE weights" (1.4315) and "75 members with
  ORACLE membership" (2.3055) is the cleanest statement of where the architecture's readout
  loses, and a readout that emits a sparse weighted combination is the one class this lane did
  not close by ceiling.
- The FAIL18 / 108 split is in every row of the ladder table and no rung behaves differently on
  the two strata except rung 6's cosine (p = 0.0006 against the random-18 null), which points
  the wrong way.
