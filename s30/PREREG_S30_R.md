# PREREG S30 lane R — IS NATIVENESS RECOGNISABLE FROM THE GEOMETRY OF A SINGLE STRUCTURE?

Charter lead L11. Written and committed BEFORE any aggregate over more than the one probe target
named in "what I have already seen". Contract rules 13 (0.7x MDE = NOT MEASURED), 17
(multiplicity), 20 (shrink signature), 27 (prereg precedes artefact) carry over from S29.

---

## 1. The question, and why the existing record does not answer it

S28-L48: **20 of 31 native-free scorers prefer the projected production average to a 0.25 A ORACLE
structure**, fold CI below 0.5. S29-L12 explains it: the perception-distortion theorem forbids any
realism measure from preferring the distortion-optimal answer ACROSS realism levels.

That explanation is also the gap. Every rung in that comparison differs in **kind** as well as in
nativeness — PROD is a contracted coordinate average, the ORACLE rungs are circuit outputs, the
controls are perturbations of PROD. So the measurement cannot separate

  (H-kind)      "scorers detect which construction produced the structure", from
  (H-native)    "scorers cannot see nativeness".

S29-L50 pushes the other way: **9 of 31 channels keep in-band ordering skill >= 0.15** after
partialling out both rank(Rg) and rank(|Rg - median Rg|) — LEG_torsion +0.181 among them, and it
is a function of the structure rather than of the distogram. So *something* orders structures
within a pool. What is not known is whether that something is **nativeness** or **pool typicality**,
and whether it survives into the near-native regime the pool never reaches (pool best ~1.8-2.4 A).

This lane builds the ladder that separates them.

## 2. The instrument

Every rung is an **ideal-geometry backbone built from (phi, psi)** via
`core.geometry.build_backbone_batch`. No projection, no coordinate averaging, no contraction:
bond lengths and bond angles are IDENTICAL at every rung by construction, so the 25.8%-contraction
confound and the projection step are both absent from the comparison.

Perturbed residues receive torsions drawn from the **fold's leakage-safe Ramachandran table**
`s8/generate_rama.npz[fold]` (the same table `ham_lib.h_rama` scores against, and the same
leave-fold-out rule the shipped distogram trains under). So a perturbed residue's local
conformation is a draw from the distribution the native's own torsions are typical draws of:
**local realism is matched by construction, not by assertion.** It is audited, not assumed
(section 5, audit A2).

Three families, per target, over all 126 dev targets:

| family | anchor | what varies | label |
|---|---|---|---|
| **A** | the native's own (phi, psi) | m in {1,2,3,4,6,8,12} residues resampled, R draws each | CA-RMSD to `cand.nat_ca` (ORACLE) |
| **B** | a random REAL pool member's (phi, psi) | same m grid, same R | CA-RMSD to the native AND to the anchor (both ORACLE) |
| **C** | fixed rungs | — | the rebuilt native; the 9 cached rungs of `s29/results/s29_D_ladder_structs/` (PROD, circ_best 0.29 A ORACLE, circ_s0, sub0, circ_opt, NATIVE, RAND_SIGNED[0], GAUSS_0.3[0], GAUSS_MATCHED[0]); the 500 real pool members |

**The key property, measured on the probe target 1A13 before this file was written** (this is the
whole of what I had seen): the rebuilt native sits **0.284 A** from `nat_ca` — that is the
ladder's floor and every sentence quoting it will say so — and **at m = 1, one single resampled
residue, CA-RMSD to the native spans 0.284 to 4.685 A**. At a FIXED perturbation budget, with
identical local-geometry statistics, two structures differ by one torsion and by 4.4 A of
nativeness. Ordering THAT is the experiment.

**Scorers.** All 18 `s27/ham_lib.py` CHANNELS + DIS + LEG + the 11 LEG_<term> components
(including LEG_torsion), ~31 in total — the same library S28-L48 and S29-L50 measured, so the
result is directly comparable to both.

**Pool-consistency reference.** CONS, DMAP_CONS, TORS_CONS and POOLGO need a reference set. It is
computed over the real 500-member pool UNION the ladder, so every rung — ladder, native,
production, pool member — is scored against an IDENTICAL reference. The ladder's contamination
share of that reference is reported, not hidden.

## 3. What is measured

Per target, then aggregated over the 126 with fold-clustered CIs through `s24.stats_lib.compare`
on the pinned folds. `partial_spearman(..., vshape=True)` (`s29/s29_T_compactness.py:68`) is used
for every partialled number: it removes rank(Rg) AND rank(|Rg - median Rg|), per S29-L50's repair.

1. **rho_A** = within-m Spearman(channel, RMSD-to-native) on ladder A, averaged over the m strata.
   Realism-matched and budget-matched ordering of nativeness.
2. **rho_B** = the same on ladder B, ordering RMSD-to-native. The anchor is not the native, so
   there is no anchor confound at all; this is the deployment-relevant regime.
3. **rho_ANCHOR** = within-m Spearman(channel, RMSD-to-**anchor**) on ladder B.
   **This is the control in the operator's own space.** A channel that orders distance from an
   arbitrary real anchor as well as it orders distance from the native is measuring PERTURBATION
   MAGNITUDE, not nativeness. Omitting this control is the project's most repeated error
   (`control-must-match-the-operators-space`, 3 instances in 2 sprints).
4. **PREF** = share of targets on which the channel scores the nearest-native ladder-A rung
   (<= 1 A) BELOW the projected production chain; ties 0.5. Beside it, the mandatory S28-L36
   **pool-member control** (share on which it scores a random real pool member below production)
   and their paired contrast through `ST.compare`.
5. Every one of 1-4 repeated **partialled** on rank(Rg) and rank(|Rg - median Rg|).
6. **LFO combination.** A ridge/logistic over the standardised channels fit on 4 folds and
   evaluated on the 5th, discriminating near-native (<= 1 A) from production. Priced against
   (a) the same fit with shuffled labels and (b) **the anchor-control labels** — because the S28
   nested logistic reached 0.960 held-out sign accuracy while its RAND_SIGNED control reached
   0.952 (S28-L48), i.e. an unmatched combination re-measures garbage rejection
   (`decoy-bank-not-a-pool-proxy`).

## 4. THE FALSIFIER — F-R1

> **Nativeness is recognisable from single-structure geometry on this instrument** if and only if
> at least one native-free channel, or the leave-fold-out combination, clears **BOTH** clauses:
>
> **(i) ORDERING.** partialled rho_A >= **+0.25**, fold CI excluding zero, **AND**
> (rho_A - rho_ANCHOR) >= **+0.10** with the paired fold CI excluding zero.
>
> **(ii) PREFERENCE.** PREF(<= 1 A rung vs production chain) >= **0.65** with the fold CI above
> 0.5, **AND** (PREF - pool-member control) >= **+0.10** with the paired fold CI excluding zero.
>
> **If nothing clears both, the answer to L11 is NO on this instrument**, and I will report it
> that way — plainly, with what it forecloses — rather than quoting the best cell.

Bars justified in advance, not after: +0.25 is below the best partialled in-band skill the record
owns (DMAP_CONS +0.316, S29-L50) and above the 9-channel survivor floor (+0.15), so it asks for a
real effect without asking for a new record. 0.65 is comfortably above the 0.5 coin toss and far
below the 0.93 at which DIS prefers production (S28-L49) — a channel that recognises nativeness
should beat a coin toss by a visible margin. The +0.10 margins are the smallest contrast the
per-comparison MDE is expected to resolve at n = 126.

**Registered prior (mine, before the run): F-R1 does not fire, at about 4 to 1.** I expect
rho_ANCHOR to track rho_A closely for the consistency and geometric families, and I expect
LEG_torsion and RAMA to carry a small genuine rho_A that fails clause (ii). I am recording this so
that a null cannot later be presented as a surprise, and so that a POSITIVE counts double against
my own expectation.

**Secondary questions, reported but NOT part of F-R1** (they cannot be promoted after the fact):
S1 does rho_A depend on m (does skill die as the perturbation grows)? S2 does any channel order
ladder B's RMSD-to-native (rho_B) better than ladder A's? S3 where does the native sit in the
percentile of its own ladder A?

## 4b. Three additions required by the coordinator, registered before any result

Added on the coordinator's instruction after the design was approved and BEFORE any aggregate
existed (the approval message is in this lane's transcript; this file is committed before the
first scoring pass). They are registered here so they cannot be read as chosen after the fact.

### D1 — the locality decomposition: does F-R1's negative half follow by construction?

If ONE resampled torsion swings global CA-RMSD by 4.4 A while every local geometric statistic is
still a draw from the same Ramachandran distribution, then **a purely local channel is blind to
that variation by construction**, and its ceiling on ladder A is set by how much of the
RMSD-to-native variance is reachable from local features at all.

Measured, within each m stratum, on ladder A:

- **R2_local** — the share of var(RMSD-to-native) explained by a least-squares fit on the full
  LOCAL feature block: the per-residue (phi, psi) sin/cos pairs, the per-residue Rama log-p, the
  identity of the resampled residues (a 0/1 indicator vector), and the signed circular torsion
  deltas from the anchor. This is an ORACLE regression fit IN-SAMPLE per target, so it is an
  **upper bound** on what any local channel could extract, and it is labelled as such everywhere.
- **R2_global** — the same for a global block: the candidate's full CA distance map, its Rg, and
  its contact map.
- The **gap** R2_global - R2_local is the size of the window in which a globally-reaching channel
  (DMAP_CONS, CONS, POOLGO, CONTACT, CONTACT_LL, DIS) could possibly operate where a local one
  (RAMA, LEG_torsion, DSSPHB, CAGEO) could not.

**If R2_local is near zero within m, then a null for the local channels is a property of the
instrument and not an empirical miss**, and the live question narrows to the globally-reaching
family. That converts the verdict into a mechanism, and it is reported whichever way it comes out.
Note the asymmetry that keeps this honest: a HIGH R2_local would make my registered prior *harder*
to defend, not easier.

### D2 — the resolution curve: at what delta-RMSD does ordering die?

F-R1 is pass/fail. It does not say at what separation discrimination exists, and a channel that
tells 4 A from 0.3 A but not 1.5 A from 1.0 A is coarse triage and nothing else. So, within each m
stratum, every pair of ladder-A rungs is binned by |delta RMSD| into
{0-0.25, 0.25-0.5, 0.5-1, 1-2, 2-4, >4} A, and for each channel and bin I report the **pairwise
concordance** (share of pairs the channel orders correctly, ties 0.5) with a fold-clustered CI and
the per-bin MDE. The reported number is **the smallest |delta RMSD| bin at which concordance stays
above 0.5 with the fold CI excluding it** — the channel's resolution. Chance is 0.5 exactly by
construction, so this needs no separate null.

### D3 — the anchor control's own confound, measured before the verdict is read

On ladder B the anchor is a real pool member, and pool members are themselves somewhat
native-like, so RMSD-to-anchor and RMSD-to-native are **positively correlated by construction**.
That correlation shrinks the rho_A - rho_ANCHOR contrast and biases F-R1(i) toward declaring a
channel a mere magnitude-detector.

So `corr_spearman(RMSD-to-anchor, RMSD-to-native)` within each m stratum is measured per target
and reported **before** the F-R1 verdict is read, with its distribution over the 126. **If its
median exceeds 0.5, clause (i)'s +0.10 margin is declared possibly unreachable for reasons
unrelated to the channels, and the clause is reported with that caveat attached rather than
scored as a clean fail.** The contrast is additionally reported in a partialled form —
rho(channel, RMSD-to-native) with rank(RMSD-to-anchor) removed — which is immune to this
confound and is the number to quote if the raw contrast is compressed.

## 5. Audits that must pass before any headline is read

- **A1 reproduction.** The 9 cached rungs' ORACLE RMSDs recomputed here must equal
  `s29/results/s29_D_ladder_structs/<pdb>.npz :: rmsd_*` to 1e-6. A discrepancy is the finding.
- **A2 realism flatness.** For ladder A, the per-rung RAMA score, the mean |Rg - median pool Rg|
  and the EXVOL clash term are regressed on RMSD-to-native WITHIN each m stratum. If realism is
  NOT flat in RMSD within m, the ladder is a garbage-rejection instrument and every positive on it
  is reported as such. **This audit can only weaken my own positives.**
- **A3 NaN poison.** The generator must never read `nat_ca` except to LABEL a rung. Labels are
  computed after all scores. Verified by running the scoring pass with `nat_ca` set to NaN and
  checking every channel value is bit-identical.
- **A4 multiplicity.** ~31 channels x 3 rho columns. The max-over-channels is priced against a
  per-target sign-flip null (500 draws), per S28-L34(e). A positive whose p_max exceeds 0.05 is
  reported as NOT SIGNIFICANT AFTER MULTIPLICITY.
- **A5 tie discipline.** No `np.argmin` on a possibly-tied signal; outcomes averaged over the tied
  argmin set (`ST.argmin_tied`), per `consensus-is-the-only-in-band-discriminator`.

## 6. What I had already seen when writing this

One probe target, 1A13: the rebuilt-native floor (0.284 A), the per-m mean/min/max RMSD table of
the generator, and the pool's RMSD quantiles. No channel value, no correlation, no preference,
no aggregate over any set of targets. The prototype that produced it is
`s30/s30_R_ladder.py`'s generator, run standalone.

## 7. Multiplicity budget for this lane

0 endpoint (deployable RMSD) comparisons. ~31 channels x 3 ordering columns x 2 partial settings
reported as ONE table, none selected; PREF reported for all ~31 with one max-over-channels null.
Nothing in this lane proposes an arm.
