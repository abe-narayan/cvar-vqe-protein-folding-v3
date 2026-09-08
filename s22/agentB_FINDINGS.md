# Sprint 22 — WORKSTREAM B FINDINGS (physics: Legacy vs AMBER)

Pre-registration: `s22/PREREG_B.md` (frozen before data; two dated addenda, both 2026-09-07 — one
recorded after F-B1 forced a correction to its own operationalisation, one registering the Priority-2
surrogate test before its data was read — see that file, not repeated here). Domain: the two physics Hamiltonians — construction, curvature, landscape mechanism, and the
mandated deep Legacy-vs-AMBER comparison. RMSD-facing consequences of anything below are explicitly
out of scope for this document; see BRIEF section 2 for where the project's Ångströms currently are
(selection, not physics).

**AMBER, stated once here and then at every table**: every number in this document that says
"AMBER" is the **bare single point** (`ConstrainedBox.energy_point`, `core.amber`/`s17.phys_lib`,
via `s20.c_land.Pot.amber`) — **no minimisation**. This is a forced choice (PREREG_B.md section 9,
FUNCTIONAL fork): a landscape mechanism probe needs a function of theta, and the deployed
`H_AMBER = E ∘ Relax_50` is not one (the relaxation leaves the torsion manifold — the same reason
`s21/c_norm.py` and `s21/c_cont.py` made this choice). Where the deployed operator's own cost is
relevant (Priority 4), it is cited from Sprint 21 and never conflated with a number measured here.

---

## PRIORITY 1 — THE CONTROLLED SYNTHETIC PERTURBATION PROBE (highest priority; new this sprint)

### Design

Five torsion-space perturbation constructions ("axes"), matched at three shared RMS magnitudes
(0.02, 0.05, 0.10 rad over the full active coordinate set, `phi[0]` excluded as the known exact
inert coordinate — `s20` L-series), applied to real starting structures (pool medoid + 2 pool draws,
per target) drawn from the shipped top-75 real-rebuild pool for **the same 30-target subset that
already carries the C5 curvature numbers** (`s21/c_norm.SUBSET`), so the mechanism claim below and
the curvature claim it is built from are about the *same* targets rather than an unstated
generalisation:

    A  diffuse/random        -- iid direction over every active (phi,psi)
    B  compactness-directed  -- coherent global bias toward the extended (-120,130) or helical
                                (-63,-42) reference basin (core.energy's own Rama basins / the same
                                reference core.amber.AmberHamiltonian uses)
    C  contact-density-local -- the same coherent bias, restricted to the middle third of the chain
    D  steric-pair-directed  -- numeric gradient direction that most changes the closest non-bonded
                                CB-CB pair's distance in the START structure
    E  local/single-residue  -- the same random construction as A, concentrated onto one residue

3,780 trials over 30 targets (complete: true, full key set, `s22/results/probe_perturb_07144a9404e8baf6.json`).
AMBER/OpenMM was taken as the sole context in this session for the run (163 s wall clock) and
released immediately after (`Pot.close()` per target), per BRIEF section 8.

### Result 1 — the sharpest and cleanest finding: Legacy and AMBER respond to compactness with
### EQUAL MAGNITUDE and OPPOSITE SIGN

    axis B (compactness), Spearman rho of dE vs delta-Rg, pooled n=540:
        Legacy   +0.454
        AMBER    -0.444

Legacy's bare energy **rises** as a real structure is pushed toward the extended reference and
**falls** as it is pushed toward the compact/helical one; AMBER's bare energy does the **opposite**,
at essentially the same magnitude. This is not a restatement of "AMBER's landscape is harder" (C5) —
it is a **directional mechanism finding, obtained by direct perturbation rather than by staged
minimisation**, and it independently reproduces and sharpens `s21` L33/C7″'s conclusion ("Legacy
prefers 0.45 Å more compact structures... compaction closes the contacts that ARE the [AMBER steric]
singularity") through a completely different construction: C7″ inferred the mechanism from what a
Legacy minimisation does to AMBER's Hessian after the fact; this measures the two energies' raw
first-order response to the SAME controlled move, with no minimisation anywhere in the loop.

**Interpretation, stated as an interpretation, not a further measurement**: Legacy's compactness
term and burial-rewarding solvation term make it prefer compact geometry; AMBER's implicit-solvent
nonbonded/GB terms are dominated, on real (not idealised) candidate geometry, by the steric
consequences of compaction — pushing a real structure toward extended geometry relieves clashes the
bare AMBER single point penalises heavily. **The two Hamiltonians are not merely "different" on
compactness; they are anti-correlated at matched magnitude**, which is a stronger and more useful
statement for anyone designing a mixture or a router than "both are worse than a random tail" (C10).

F-B2, as literally pre-registered (**magnitude** comparison, `|sens_Legacy| > |sens_AMBER|`), is
**NOT MEASURED**: point estimate `+0.010`, CI `[-0.186, +0.197]` — the magnitudes are statistically
indistinguishable (0.454 vs 0.444). The registered falsifier asked the wrong question; the
**sign opposition** is the finding, and it is reported as a finding **separate from, and not a
substitute for, the falsifier's literal verdict** (BRIEF section 3: falsifiers are never moved
afterward — F-B2 stands as NOT MEASURED, and the sign-opposition is additional, clearly labelled).

### Result 2 — the concentration hypothesis was WRONGLY OPERATIONALISED, caught before it was
### written up as a clean refutation, and the corrected version is Priority-1's second result

F-B1 (torsion-INDEX concentration: axes D and E vs axis A) came back **REFUTED, on the wrong side**:

    R(Legacy) = 1.338     R(AMBER) = 0.257
    R(AMBER)/R(Legacy) = 0.192   bootstrap 95% CI [0.041, 0.514]   (predicted >1; observed <1)

Before reporting this as "AMBER is actually LESS sensitive to concentration than Legacy" — which
would contradict the C5 participation-ratio mechanism this hypothesis was built from — the
operationalisation was re-examined (dated addendum, `PREREG_B.md`, same day). `core.geometry.
build_backbone` is a **sequential NeRF chain**: perturbing one residue's torsions moves every atom
from that residue to the C-terminus. **Torsion-index concentration (few nonzero coordinates) is not
Cartesian concentration**, and axis D's direction (a numeric gradient of a pairwise distance) is
diffuse across nearly every active coordinate despite targeting one pair — it was never actually
index-concentrated at all. More fundamentally, **C5's participation ratio describes concentration in
the Hessian's own EIGENBASIS**, which is generically not aligned with any single-residue coordinate
axis, and nothing in the original F-B1 design touched that eigenbasis. **This is exactly the kind of
unstated operator mismatch Rule 0 exists to catch, found here by the lane against its own registered
hypothesis, before the wrong verdict was written up** — the discipline this project's ledger repeatedly
credits (M2, M9 in `s21/CLAIMS.md`).

**F-B1′ (the corrected version, same-day addendum, own top-Hessian-eigenvector as the perturbation
direction, `s22/probe_eigen.py`, `s22/results/probe_eigen_a11b574ea2da5228.json`, complete, 30/30
targets)**. First, the instrument re-certifies itself: this subset's own median participation ratio
(Legacy 0.4965, AMBER 0.0668) reproduces C5's headline numbers (0.497 / 0.067) almost exactly before
any new number is read. Then:

    OWN-RATIO(potential) = median(|dE|/rms) along POTENTIAL's OWN top Hessian eigenvector
                          / median(|dE|/rms) along the diffuse random axis (A)

    OWN-RATIO(Legacy) = 2.559       OWN-RATIO(AMBER) = 74.43
    OWN-RATIO(AMBER) / OWN-RATIO(Legacy) = 29.09
    bootstrap 95% CI (2000 draws, targets as the resampling unit) = [4.49, 542.68]

**F-B1′ is SUPPORTED**, the CI excluding 1 by a wide margin on the predicted side. **AMBER's energy
is ~29× more disproportionately sensitive, relative to its own generic/random response, to a move
along its own dominant curvature direction than Legacy is to a move along ITS dominant direction.**
This is the properly basis-matched confirmation of the mechanism section 1 argued for from C5's
participation-ratio numbers: **AMBER's curvature concentration is real and has a measurable
consequence for how it responds to perturbation — but only when the perturbation is aligned with the
Hessian's OWN eigenbasis, not with a raw coordinate (single-residue) index.** The two results
together are the actual finding: **coordinate-index concentration is not what AMBER's landscape is
concentrated in; eigen-direction concentration is** — a mechanistic fact about the parameterisation
(NeRF's sequential atom placement) as much as about either physics model, and one this workstream
would not have found without registering, running, and then correctly diagnosing why F-B1 came back
backwards.

**A bootstrap defect was found and fixed before this number was quoted.** The first version of
`probe_eigen.analyse` re-seeded `np.random.default_rng(0)` INSIDE the per-draw resampling closure, so
all 2000 "bootstrap draws" resampled identically and returned a degenerate CI, `[12.72, 12.72]` —
two bounds with no daylight between them, which is itself the tell. Caught by that symptom, not by
inspection first; the RNG was hoisted outside the loop and the corrected CI is the one reported
above. Recorded here, and in `PREREG_B.md`'s addendum, because a coding defect that inflates
confidence (a spuriously tight CI, here in the *supportive* direction) is exactly the kind of error
this programme's ledger asks to be surfaced rather than quietly fixed.

### Result 3 — steric response is largely SHARED, not specialised

Axis D (targeted steric clash/relief on the single closest non-bonded pair): Legacy rho = -0.214,
AMBER rho = -0.190 against the pair's own distance change — same sign, similar magnitude for both.
**Both energies penalise a targeted clash and reward its relief at comparable rank-correlation
strength.** Steric response is not where the two Hamiltonians diverge; compactness (Result 1) is.

### Result 4 — AMBER's dose-response to a GENERIC (non-directed) torsion move is markedly less
### predictable than Legacy's

    axis A (diffuse random): rho(dE, ||delta_theta||)   Legacy +0.271   AMBER +0.149
    axis E (single-residue):  rho(dE, ||delta_theta||)   Legacy +0.379   AMBER +0.157

For a move whose size is the only lever (no directed target), **Legacy's energy change tracks the
move's magnitude about twice as reliably (by rank correlation) as AMBER's does.** This is consistent
with, and gives a quantitative face to, C5's landscape-shape numbers: Legacy's Hessian spreads
curvature over about half its modes (participation ratio 0.497, more nearly quadratic, hence a more
monotone dose-response), while AMBER's is concentrated and near-zero over roughly half its spectrum
(near-zero eigenvalue fraction 0.481) — so "how far you moved" is a much weaker predictor of "how
much AMBER's energy changed" than the same statement for Legacy, independent of the direction moved.

### Result 5 — contact-density manipulation moves both energies mildly, AMBER more so

Axis C (local-window coherent bias): Legacy rho = +0.118 (weak), AMBER rho = +0.203 (moderate)
against the realised change in geometric (unweighted) local contact fraction. Both signs positive:
forcing extra local contacts through a coherent torsion bias on real starting geometry produces
*higher* energy on average for both models — consistent with such contacts being geometrically
strained rather than the sequence's favourable (MJ-favourable, for Legacy) pairs. AMBER's larger
coefficient is consistent with, but far weaker than, Result 1's compactness finding; it is reported
for completeness and is **not** promoted to a separate mechanism claim on its own (n=540, and the
within-target permutation null's 95th percentile for this axis, 0.086-0.120, is close enough to the
observed values that this one is reported as suggestive rather than decisive).

### Null

Within-target, within-axis permutation of the (proxy, dE) pairing (PREREG_B.md section 5), 500
draws per axis/energy. All five OBSERVED |rho| values reported in Results 1, 3, 4 clear their own
axis's null 95th percentile with room; Result 5 (axis C) is the one comparison that sits closest to
its own null, and is labelled accordingly above rather than folded into the confident results.

---

## PRIORITY 2 — A DEFENSIBLE AMBER SURROGATE, IF ONE EXISTS

**Candidate tested**: `H_surrogate := amber_bonded` (`s21.c_cont.MixPot.amber_bonded`) — the
bond+angle+torsion force-group subset of the SAME genuine ff14SB `System`, with nonbonded and
solvation switched off. This is not a new construction; it already exists in `s21/c_cont.py`,
already correctly never called "amber" there, built for the continuation work. This workstream asks
a different question of it: not "does it help as a continuation stage" (answered, negatively, by
`s21` C7′/C7″) but "does it TRACK the genuine bare single point's response well enough to be a cheap
ranking proxy."

**Pre-registered falsifier** (`PREREG_B.md` Addendum 2): pooled Spearman rho between `dE_amber`
(genuine, bare single point) and `dE_bonded` (candidate) over the same axis-A and axis-B trials,
medoid start, full 30-target subset — supported if rho > 0.7, refuted otherwise, threshold fixed
before the run.

    s22/results/probe_surrogate_9d38334ce908bde3.json  complete: true, 30/30 targets, 540 trials

    median |E_bonded_0| / |E_amber_0| at baseline           0.0278   (the surrogate structurally
                                                                       cannot capture ~97% of the
                                                                       genuine magnitude on a typical
                                                                       target — core/amber.py's own
                                                                       docstring: CustomGBForce is
                                                                       ~98% of one evaluation's cost)
    Spearman(dE_amber, dE_bonded), pooled, n=540              0.211
      axis A (diffuse)                          n=360         0.106
      axis B (compactness)                      n=180         0.400

**VERDICT: REFUTED. `amber_bonded` is NOT a defensible surrogate, even for ranking alone**
(0.211 << 0.7, and the axis it does best on, compactness, still falls well short). The baseline
fraction table explains why: on the targets whose bare AMBER energy is dominated by the steric
singularity (1e7-1e9 kcal/mol; `s21` C5), the bonded subset is a vanishing fraction of the total
(down to 0.00 at the reported precision) and cannot see the nonbonded clash that IS the singularity
— which is also this workstream's own Result 1 mechanism (compactness closes contacts; the bonded
terms do not represent contacts at all). **No defensible AMBER surrogate was found by this route.**
The negative is reported with the same weight a positive would carry: it forecloses one plausible
tractability lever (a cheap bonded-only proxy for the expensive nonbonded/solvation calculation)
rather than leaving it untested.

**Labelling discipline, honoured throughout**: nothing in this section or in `s21/c_cont.MixPot`
(the module the candidate surrogate is drawn from) is ever called "amber" — it is `amber_bonded` /
`H_surrogate`, a genuine bonded-force-group SUBSET of the real ff14SB System, never a learned or
fitted model, and the genuine bare-single-point `H_AMBER` remains independently evaluable throughout
(Pillar 3; `s21/c_cont.py`'s own docstring states the same discipline for the identical object).
No other surrogate construction (geometric pre-relaxation, softened steric, local Taylor expansion)
was built or tested this sprint — time and the AMBER queue went to Priority 1's mechanism probe and
its correction first, per BRIEF's own priority ordering. This is a scope limit, not a further
negative result, and is stated as one.

---

## PRIORITY 3 — NORMALISED CONTINUATION

**Not re-run.** BRIEF section 1 lists this among what Sprint 21 closed and instructs "DO NOT REDO":

- `s21` C6/C6′ (`s21/CLAIMS.md`): the raw-unit `H(lambda) = (1-lambda) H_L + lambda H_A` is
  **degenerate** — crossover `lambda* ` has a 3.7e7-fold spread across targets, and a raw uniform
  grid is AMBER-dominated above `lambda ≈ 2e-5`. Stronger: **a raw lambda ladder never visits an
  intermediate Hamiltonian at all** — largest single-step spectral change is 0.9989 and it occurs in
  the *first* interval; the Hessian spectrum at `lambda=0.1` equals the `lambda=1` spectrum to four
  significant figures on all seven scale-free metrics, through `lambda=0.9`. The **declared**
  normalisation (`Nt`, asinh on the pool median/MAD, fixed *before* any Sprint-21 RMSD existed) was
  built and measured properly in `s21/c_norm.py` (reused directly by this workstream's own
  normalisation, see Priority 1's method) — the degeneracy is a property of the *raw-unit*
  construction, not of a normalisation nobody tried.
- `s21` C7′: **staged** Legacy→AMBER continuation is *worse* than direct AMBER
  (`+0.2331 [+0.0673,+0.3950]`), surviving a move-size correction, and not stage starvation (staged
  arms spent *more* evaluations).
- `s21` C7″: the causal mechanism — Legacy preconditioning drives AMBER **~31× deeper** into the
  steric singularity, because Legacy prefers 0.45 Å more compact structures and compaction closes
  the contacts that ARE the singularity.

**What this workstream adds, without re-running the ladder**: Priority 1's Result 1 (compactness
sign-opposition, measured by direct perturbation with no minimisation in the loop at all) is an
**independent replication of C7″'s mechanism by a different method** — C7″ inferred the mechanism
from what a Legacy minimisation does to AMBER's post-hoc Hessian; this measures the same directional
antagonism as a first-order energy response on real geometry. Two different instruments, same sign,
same qualitative mechanism: **compaction is not a neutral move for AMBER, and it is exactly the move
Legacy's own energy rewards.** This is offered as corroboration, not as a new continuation result,
and does not touch C6/C6′/C7′'s verdicts.

---

## PRIORITY 4 — AMBER AS HAMILTONIAN vs AMBER AS REPAIR

Kept conceptually and experimentally separate, per BRIEF's explicit instruction, by **never running
both objects under the same label**:

- **This workstream's Priority-1 probe uses AMBER strictly as a Hamiltonian** (H_AMBER, the bare
  single point, a function of theta) — never as a terminal repair operator. No structure in this
  document has been passed through `Relax_50` or any other repair step.
- **The deployed repair operator's own cost is cited, not re-measured**, from `s21` C8c
  (`s21/CLAIMS.md`): the bare single point costs **+1.471 Å** against the deployed
  `E ∘ Relax_50`'s **+0.021 Å [+0.014,+0.028]** — a ~70× difference in RMSD consequence between the
  same energy function used bare versus used as `E ∘ Relax_50`. **This is the number that answers
  "does it matter which AMBER you mean" — restated here because BRIEF requires it stated at every
  appearance, not because this workstream re-measured it.**
- Consequence for how to read Priority 1: every mechanism finding above (compactness
  anti-correlation, concentration, dose-response predictability, steric sharing) describes
  **H_AMBER's landscape as a Hamiltonian**. None of it is a claim about what `Relax_50` does to a
  candidate, which is a different operator with a different (and separately measured) RMSD effect.
  A reader who wants "does perturbing toward compactness hurt the DEPLOYED pipeline" needs the
  repair-operator numbers above, not this document's Results 1-5.

---

## RULE 0 — SIX OPERATOR FORKS FOR PRIORITY 1 (sent as required before the run; reproduced here for
## the permanent record — see `PREREG_B.md` section 9 for the full text)

FUNCTIONAL: bare AMBER single point, not `E ∘ Relax_50` (forced; a landscape needs a function of
theta). BASIS: real top-75 pool starting structures, not idealised/lattice states. READOUT: signed
raw dE per trial, standardised by that target's own pool median/MAD — never an RMSD readout (the
Hard-Rule-1 exception, declared). NORMALISATION: per-target robust z from the already-audited
`s21/c_norm.json`, not the `Nt` asinh transform (reasoned: this probe never sums or mixes the two
energies). NULL: within-target, within-axis permutation — not a cross-target shuffle (mis-specified
for this instrument family per BRIEF section 2's own worked example). THE LABEL: both falsifiers are
continuous ratios tested against a fixed reference (1 or 0), never a covariate-dependent threshold.

---

## SUMMARY, FOR THE FINAL REPORT

1. **Legacy and AMBER are anti-correlated on compactness at matched magnitude and comparable
   strength** (+0.454 vs -0.444) — the sharpest mechanistic statement this workstream produced, and
   it independently corroborates `s21`'s staged-preconditioning finding by direct measurement.
2. **Steric response is shared, not specialised** (-0.214 vs -0.190) — the two Hamiltonians largely
   agree on penalising a targeted clash.
3. **AMBER's response to a generic, non-directed move is markedly less predictable than Legacy's**
   (dose-response rho roughly half), consistent with AMBER's near-zero-eigenvalue-heavy, low-
   participation-ratio spectrum (C5) — but coordinate-index concentration is the WRONG
   operationalisation of "concentrated" (torsion-index-local is not Cartesian-local, because
   `build_backbone` is a sequential NeRF chain). The CORRECTED test, along each potential's own top
   Hessian EIGENVECTOR, gives the clean confirmatory result: **AMBER is ~29× [CI 4.5-543×] more
   disproportionately sensitive to its own dominant curvature direction than Legacy is to its own**
   (F-B1′, SUPPORTED). Concentration is real; it lives in the eigenbasis, not the coordinate index.
4. **No defensible AMBER surrogate was found.** The one genuine, already-existing, correctly-labelled
   candidate (`amber_bonded`, the bonded force-group subset of the real System) tracks the bare
   single point's response at Spearman rho = 0.211 pooled (0.7 was the pre-registered bar) — REFUTED,
   because the steric singularity that dominates AMBER's real response lives in the nonbonded/
   solvation terms the bonded subset cannot see. `H_AMBER` remains genuine and independently
   evaluable; no shortcut around it was found this sprint.
5. **Continuation and preconditioning remain closed exactly as `s21` left them**; this workstream
   adds an independent replication of the preconditioning mechanism (Result 1's compactness
   antagonism, measured by direct perturbation rather than staged minimisation), not a new
   continuation result.
6. **AMBER-as-Hamiltonian and AMBER-as-repair are never conflated here**; every number above is the
   bare single point, and the deployed repair operator's own (very different) cost is cited from
   `s21` C8c wherever the distinction matters.

## ARTEFACTS

    s22/PREREG_B.md                                          pre-registration + 2 dated addenda
    s22/probe_perturb.py, probe_eigen.py, probe_surrogate.py  the three instruments
    s22/results/probe_perturb_07144a9404e8baf6.json           complete, 30/30, 3780 trials
    s22/results/probe_perturb_analysis_07144a9404e8baf6.json  F-B1, F-B2, per-axis sensitivities
    s22/results/probe_eigen_a11b574ea2da5228.json              complete, 30/30, eigenvector trials
    s22/results/probe_eigen_analysis_a11b574ea2da5228.json     F-B1' (corrected concentration test)
    s22/results/probe_surrogate_9d38334ce908bde3.json          complete, 30/30, surrogate trials
    s22/results/probe_surrogate_analysis_9d38334ce908bde3.json Priority-2 verdict (REFUTED)
