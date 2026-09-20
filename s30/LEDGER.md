# SPRINT 30 LEDGER

Entries are `## S30-L<n> -- TITLE (date time, lane)` with the verdict in the heading.
Run `date` in the same command as the append. Corrections are annotated in place with the
original wording left standing.

---

## S30-L0 -- THE CHARTER, TWO READING-LIST CORRECTIONS, AND THE ARITHMETIC THAT SHOULD GOVERN THE SPRINT: **FIXING TEN TARGETS BEATS IMPROVING ALL 126 BY 0.20 Å** (2026-09-20 12:35, coordinator)

### The charter
Saved verbatim as `s30/BRIEF.md`. One hard constraint (CVaR-VQE remains the spine and main
scientific object), one endpoint (mean built-chain Cα RMSD, 126 targets), current production
**3.2105 Å**, primary target < 3.0, ambitious < 2.5. Section 8's twelve leads are explicitly a
**leads register, not a task list**; the charter states that a sprint ignoring ten of them and
finding the mechanism is a success, and one executing all twelve and finding nothing is not.

### Two corrections to the reading list, made by checking rather than assuming
The charter names `s27/REPORT_S29.md` and `s27/LEDGER.md through the S29 close`. Neither is right:
the S29 report is at **`s29/REPORT_S29.md`** (1,619 lines) and `s27/LEDGER.md` is **S28's** ledger
(4,215 lines, ends S28-L50) while S29's is **`s29/LEDGER.md`** (5,823 lines, S29-L0..L57). Both
were read. There is no separate S29 retractions file; S29's 18 withdrawals live in its report
§14(e). Recorded because this project has four instances of prose naming a path that does not
exist, one of which cost a mis-planned lane-week.

### Section 7's meter already exists
`s29/s29_D_cost_audit.py` (767 lines) implements the cost-RMSD meter. S30 **verifies and extends**
rather than rebuilds; lane D owns it and copies it to `s30/s30_D_meter.py` so the S29 artefact
stays untouched.

### THE ARITHMETIC — computed at the open from lane O's S29 rows, before any new experiment

```
production, n = 126      mean 3.2105   median 2.9661
the worst 18 targets     mean 6.2758
the other 108            mean 2.6997

cap the worst 10 at 3.00 A  ->  mean 2.9074  (-0.3031)   BEATS the charter's primary target
cap the worst 18 at 3.00 A  ->  mean 2.7426  (-0.4680)
cap the worst 30 at 3.00 A  ->  mean 2.5778  (-0.6328)
worst 18 all the way to 2.50 ->  mean 2.1334  (-1.0772)

versus improving EVERY ONE of the 126 by 0.20 A  ->  mean 3.0105  (-0.2000)
```

**Fixing ten targets is worth more than improving all 126 by 0.20 Å.** For an endpoint that is a
mean over a distribution with median 2.9661 and a tail reaching 8.24 Å, a mechanism that works on
the typical target and leaves the tail alone is close to worthless, and a mechanism that only works
on hard targets can hit the primary target on its own.

**The caveat, stated with the number so it is not over-read:** capping is an ORACLE operation. It
bounds the prize; it does not deliver it. Two routes realise part of it — (a) a native-free regime
detector, or (b) a method simply better on hard pools *without needing to know they are hard*.
Route (b) requires no detection and is the stronger target.

**And a live, unexplained clue:** the record says the shipped pipeline is **worse than a
sequence-blind one on its 18 hardest targets** (blind 5.425, shipped 6.019) while sequence
conditioning is worth +0.776 Å overall. Something the pipeline does on hard targets is actively
harmful.

### Opening lane assignment (7 of a permitted 8; one slot held for what results demand)

| lane | remit |
|---|---|
| **R** | is nativeness recognizable from a single structure at all (L11) — judged to sit underneath the rest |
| **L** | literature, permanent role |
| **D** | adversary, permanent; owns and extends the cost/RMSD meter |
| **T** | theory: the bit accounting (L8), and when the CVaR tail stops being a prefix |
| **F** | the failure tail — highest leverage on the endpoint by the arithmetic above |
| **X** | divergent, permanent: should the quantum stage select at all, or generate? |
| **Q** | L5 + L6 **together** — a subset objective is pointless through an averaging readout, and a sparse readout is unusable without a way to choose its support |

Contract: `s30/S30_CONTRACT.md`, 28 rules, each annotated with the incident that paid for it.

## S30-L1 -- A CORRECTION TO THE PROJECT'S MOST-CITED NEGATIVE, POSTED BEFORE MY OWN NUMBERS EXIST: S28-L48's "20 OF 31 SCORERS PREFER PRODUCTION TO A 0.25 Å ORACLE STRUCTURE" IS A **CROSS-KIND** PREFERENCE, NOT EVIDENCE THAT NATIVENESS IS INVISIBLE -- EVERY RUNG IN THAT LADDER DIFFERS IN CONSTRUCTION AS WELL AS IN NATIVENESS, AND PERCEPTION-DISTORTION ALREADY PREDICTS THE CROSS-KIND HALF (2026-09-20 12:41, R)

Posted on the coordinator's instruction **before any lane R measurement exists**, so that the
correction cannot be read as motivated by whatever my ladder returns. My pre-registration
`s30/PREREG_S30_R.md` was committed at 12:40:35 (commit `7eabffee`) and contains no aggregate
over more than the single probe target it declares. This entry rests on reading the artefacts and
the definitions, not on new data.

### WHAT S28-L48 MEASURED, FROM ITS OWN DEFINITION

`s27/REPORT_S28.md` section 3.7, the lane C2 audit (entries S28-L35, L36, L37, L48, L49):

> BUILT CHAIN, 31 scorers (16 backbone + 15 CA on the projected chains), one max-over-31 null:
> **20 of 31 prefer the projected production average to a 0.25 Å ORACLE structure** with the fold
> CI below 0.5 (18 informative once the 5 tie-dominated are set aside; DIS on 93% of targets, LEG
> on 79%, RAMA on 59%); 6 coin tosses; CAGEO collapses to anti-recognition on ideal geometry
> (0.611 to 0.421) as registered; the one two-clause pass, CONTACT@chain 0.583, sits at the
> max-over-31 null's mean (0.579, p95 0.627, p_max 0.39).

The verdict attached to it, verbatim: *"the closure claim STANDS over the whole library. No
native-free scorer in the project recognises a near-native structure when it is offered one."*

### THE DEFECT: THE TWO STRUCTURES BEING COMPARED DIFFER IN **KIND**, NOT ONLY IN NATIVENESS

Read from `s29/s29_D_cost_audit.py:250 build_rungs`, which regenerates C2's rungs exactly and
asserts them against the S28 artefacts (`s27/results/s28_A_structs/`, 379 files, verified present):

| rung | how it is CONSTRUCTED | kind |
|---|---|---|
| PROD | `H.readout_uniform(cand, top)` -- the uniform coordinate average of the DIS top-75 | a coordinate average, contracted 25.8% (`averaging-space-beats-the-objective`), then projected |
| circ_best | `A.oracle_circuit_ceiling(...)` best-of-5 -- a quantum circuit's amplitude readout, fitted WITH the native | a circuit output |
| sub0 | `A.oracle_subspace_ls(...)` -- a least-squares fit in a random subspace, fitted WITH the native | a least-squares fit |
| RAND_SIGNED / GAUSS_* | `C2.rand_signed` / `C2.gauss_perturb` of PROD at an ORACLE-matched scale | perturbations of PROD |

So the comparison "does scorer f prefer circ_best to PROD?" varies **two things at once**:
nativeness (0.25 Å vs ~3.2 Å) and construction (a circuit amplitude readout vs a projected
coordinate average). A scorer that is entirely blind to nativeness but sensitive to construction
would produce exactly the reported table. The measurement therefore cannot separate

    (H-kind)    scorers detect WHICH PROCESS produced the structure, from
    (H-native)  scorers cannot see NATIVENESS.

### AND THE RECORD ALREADY CONTAINS THE MECHANISM FOR (H-kind)

S29-L12 imported Blau & Michaeli's perception-distortion theorem and the report states the
consequence itself (§9.1, §12.0): *for any distortion measure, the distortion-optimal estimator's
output distribution must diverge from the real one, so every realism-type scorer must disprefer
the RMSD-optimal answer.* S29-L12's own text names S28-L48 as "that statement measured".

That is the point. **If perception-distortion predicts the result from the construction alone,
then the result carries no additional information about whether nativeness is visible.** A
cross-kind preference is the theorem's content; it is not an independent test of recognition.
The two have been cited as if they were two facts. They are one fact and its explanation.

### WHAT S28-L48 DOES AND DOES NOT ESTABLISH

**It does establish**, and this survives intact and is deployably important: *no scorer in this
library will pick a near-native structure out of a set whose members were produced by different
processes* -- which is the real selection setting, and is why every native-free selector priced in
S12/S27/S28/S29 fails at the endpoint. It also stands as the measured instance of
perception-distortion on this instrument.

**It does not establish** the stronger sentence the record has been leaning on -- *"nativeness is
not recognisable"* -- because it never held construction fixed. The S29 report's own §12.0 points
the same way from the opposite side: **9 of 31 channels keep in-band ordering skill ≥ 0.15** after
partialling out both rank(Rg) and rank(|Rg − median Rg|) (S29-L50), which is not the behaviour of
a library that is blind to structural quality. §12.0 and §9.1 have been read as consistent because
one was about preference and the other about ordering; the kind confound is why they can both be
true.

### CONSEQUENCE, AND WHAT I AM DOING ABOUT IT

1. The sentence "no native-free scorer recognises a near-native structure when it is offered one"
   should be quoted as **"no native-free scorer prefers a near-native structure of a DIFFERENT
   CONSTRUCTION to the production average"**. The qualifier is not cosmetic; it is the whole
   content of the objection.
2. Lane R's ladder is built precisely to remove the confound: every rung is an ideal-geometry
   backbone built from (phi, psi) by the same function, perturbed residues draw their torsions
   from the fold's leakage-safe Ramachandran table, and the perturbation BUDGET is held fixed
   within each stratum. Kind, local realism and budget are matched; only nativeness varies.
3. The instrument's scope, stated now rather than later: its floor is the torsion rebuild at
   **0.284 Å from `nat_ca`** on the probe target, not 0 Å, and every sentence quoting it will say
   so.

**This entry changes no number.** It changes which sentence the numbers support. I am posting it
before my own results because the correction is worth exactly as much as its independence from
them, and because I was handed S28-L48 as settled and it is not.

Artefacts read: `s27/REPORT_S28.md` §3.7; `s27/LEDGER.md` (S28-L35/L36/L37/L48/L49 lines 3833,
3870, 4006, 4033, 4046); `s29/REPORT_S29.md` §9.1, §9.4, §12.0; `s29/LEDGER.md` S29-L12, S29-L50;
`s29/s29_D_cost_audit.py:250`; `s27/results/s28_A_structs/` (379 files, present);
`s29/results/s29_D_ladder_structs/` (126 files, present).
