# PREREG — LANE T — IS THE PHYSICS/REALISM CHANNEL FAMILY ON THE COMPACTNESS AXIS?

Written and committed **before** the 126-target result is read. The only compactness numbers I
have seen when writing this are the 4-target smoke test printed in S29-L41's launch note
(LEG_compactness +0.956, RG_LAW +0.904, LEG_solvation +0.766, POOLGO +0.701, DISTPOT +0.656,
ENV +0.649, DIS +0.438, AMB -0.252, CONTACT -0.287, and Rg's own in-band rho +0.190). Four
targets is a smoke test and no clause below is tuned to it. Provenance is checkable: this file's
commit precedes the commit that reads `s29/results/s29_T_compactness.json` at n = 126, and
contract rule 27 (chronology certified from git, not from anyone's word) applies to me.

## 0. Why this measurement exists

It is the substitute for the S8 free-energy diagnostic the coordinator assigned, which cannot be
run because that stage is absent from disk and from git history (S29-L41). Lane L's objection
(topic 7) is that what an entropy term tracks at peptide length is plausibly BASIN WIDTH, which
is plausibly compactness-like — in which case the free-energy class is not orthogonal to the
realism axis, and section 7's row 3 (a physics term on the emitted structure, in its free-energy
form) is not a live exit from assumption (B2) of the section 8 bound. That question is answerable
for every channel the project owns, without OpenMM.

## 1. The object

`s29/s29_T_compactness.py`, all 126 targets, 500 pool members each, 32 S27 channels
(`s27/cache/<pdb>.npz` via `s27.run_pool.channels_for`). Three quantities per (target, channel),
and they are never mixed:

| quantity | native-free? | what it says |
|---|---|---|
| `rho_rg = Spearman(channel, member Rg)` | **NATIVE-FREE** — the deployable measurement | whether the channel IS a compactness measure |
| `rho_inband = Spearman(channel, member CA-RMSD to the native)` on members with RMSD < 3 A | **ORACLE** diagnostic | whether the channel has in-band skill at all |
| `rho_inband_partial_rg` = the same, Rg partialled out of both ranks | **ORACLE** diagnostic | whether that skill survives removing compactness |

Member Rg is computed from the pool coordinates; the in-band subset is selected by the native and
is therefore ORACLE, labelled as such everywhere it appears. Nothing here selects a deployable
parameter, tunes anything, or touches the endpoint.

## 2. Statistics, fixed now

- Per channel, the fold-clustered CI on the mean of `rho_rg` over the 126 targets, via
  `s24.stats_lib.compare` against zero with the pinned folds, and `MDE = 2.8016 x SE`.
- **Contract rule 13: a channel whose |mean rho_rg| is below 0.7x its own MDE is NOT A RESULT and
  is reported as "not measured", not ranked.** The same rule applies to the two ORACLE columns.
- The 32 channels are reported together as one table. No channel is selected, so no max-over-32
  null is owed; if I later quote a single channel as the best of the 32 I will price it with one.
- Medians over targets are quoted beside means because the per-target rho distributions are
  bounded and skewed.

## 3. The falsifier, registered

> **F1 — LANE L'S OBJECTION IS CONFIRMED and section 7 row 3 CLOSES.** Both clauses must hold:
> (a) across the 32 channels, the rank correlation between `|mean rho_rg|` and `|mean
> rho_inband|` is **>= +0.4** — the channels with in-band skill are the compactness-loaded ones;
> and (b) for the channels with the most in-band skill, the median share of that skill removed by
> partialling Rg (`1 - |rho_partial| / |rho_inband|`) is **>= 0.4**.
>
> **F2 — THE OBJECTION FAILS and row 3 STAYS OPEN.** At least one channel has
> `|mean rho_inband_partial_rg| >= 0.15` with a fold CI excluding zero **and** `|mean rho_rg| <=
> 0.30`: a channel with real in-band skill that is not a compactness measure.
>
> **F1 and F2 are not exhaustive and I will not force the result into one of them.** If neither
> fires — for example if the skill columns are too weak to order at all, or if the compactness
> loading is high but the partialled skill survives on a channel that also has high `rho_rg` — I
> will report exactly that, and row 3 stays open with its status changed from "live" to
> "undecided by this instrument".

## 4. Prior, stated before the run

I expect **F1**, at about 3 to 1. The record's reasons, not the smoke test's: the in-band axis is
measured to be compactness (oracle Spearman +0.909 with the native's z-scored Rg, S29-L19, and
`in-band-ordering-is-per-target`); S7's Finding 10 already measured the depth half of exactly this
question on AMBER, where partialling Rg out cuts in-band rho to **0.159** against the distogram's
**0.479** (diff -0.320, CI [-0.425, -0.216], AMBER better on 12/70); and `docs/FINDINGS.md` states
the mechanism for the force field directly — "a force field whose minima are placed by generic
compaction rather than by sequence-specific structure will fail to refine and rank by compactness".
What would surprise me is a channel in the LEG or torsion family retaining partialled skill.

## 5. Two interpretive traps I am pre-committing against

1. **The shared-referent floor** (project memory). Rg itself has in-band skill, so a channel and
   the truth both correlating with Rg is not evidence of a mechanism. That is precisely why the
   partial is the decisive column and the raw `rho_inband` is not; I will not read a high
   `rho_inband` as skill if the partial kills it.
2. **High `rho_rg` does not by itself make a channel useless**, and I will not say it does. It
   makes the channel a compactness measure. The claim F1 licenses is narrower and is the only one
   I will make: *the in-band skill of this library lives on the compactness axis, so a new channel
   of the same family — including a free-energy channel whose entropy term tracks basin width —
   should be expected to land on that axis too and is not an independent exit from (B2).*

## 6. The caveat that travels with any verdict, whichever way it falls

Lane L's, restated: any free-energy arm would need its entropy term to rescue a channel measured
**+0.455 A WORSE than a random subset as a ranker, 5/5 folds** (S25 L16). And the gate lane L
proposed stands, amended by S29-L41: this diagnostic first; the stage would have to be **rebuilt**
from the `docs/FINDINGS.md` section B spec rather than resumed; an accuracy arm never in this
sprint.
