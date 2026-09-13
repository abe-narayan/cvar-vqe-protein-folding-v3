# IDEA -- Where the steric singularity lives: builder sidechains or pool backbones? (lane PH; own idea 3)

## Hypothesis
The AMBER single point on an unrelaxed ideal-geometry rebuild measures its worst clash (S13
walsh headline; S20 L8/L12; S25 landscape: ten candidates carry 99.7% of the variance). The
rebuild places every sidechain with the deterministic builder (`core.amber.build_full_structure`
/ `sidechains.py`, fixed `CHI_ANGLES`, no rotamer scan: "unlike the legacy model this does not
scan chi1 rotamers", `core/amber.py` `energy_from_coords`). The hypothesis is that a large share
of the singularity is the BUILDER's, not the pool's: the clashes that put half of every pool
above 1e4 kcal/mol are sidechain-involving contacts created by a fixed chi1, and a per-residue
greedy chi1 relief (each residue independently takes the best of three staggered chi1 values,
holding the others fixed; at most 3n single points per candidate) removes most of the
catastrophe. If so, the relieved single point is a materially different observable from the raw
one: it measures backbone strain, and its ranking and reject behaviour must be re-measured
rather than inherited.

## Why the record does not already close it
- S20 L8 states "the collapse is a property of the ideal-geometry rebuild itself" but did not
  separate backbone-backbone from sidechain-involving contacts. S19 section 3.1 counted
  heavy-atom contacts on N/CA/C/O/CB only (2.66 per 75 below 2.0 A).
- S13 W2 capped the energy at its 99th percentile (a modified observable) and found the Pauli
  spectrum still above Legacy's. Capping is not the same as relieving the clash: it truncates the
  ordering at the top, it does not change which candidates are at the top.
- Rank standardisation (S24 D1-A, S25) is monotone in the raw energy and therefore cannot reach
  a different ordering; a chi1 relief can.
- `converged interaction-only AMBER` (memory) relaxes everything at k=10 and did not transfer;
  it moves the backbone, which is what this idea deliberately does not do.

## Exact falsifier
Part A (native-free census, the cheap half): on every top-75 rebuild of the 126 targets, the
closest heavy-atom contact classified bb-bb / bb-sc / sc-sc, on the all-atom rebuild, and the
share of the members above 1e4 whose closest contact involves a sidechain atom. (This half is
already in `s26/ph_reject.py census`, `singularity` block, and is the optional Part IV
measurement.) Part A closes the idea if fewer than half of the members above 1e4 have a
sidechain-involving closest contact: then the singularity is backbone-owned and no rotamer
relief can move it.
Part B (AMBER compute, only if A survives): the relieved single point E_relief for the 126 x 75
members; the fraction above 1e4 before and after; then the in-band rho against ORACLE RMSD with
an Rg control, and the same reject test as `s26/PREREG_amber_reject.md` on E_relief in place of
the raw energy, same controls, same falsifier.

## Expected effect against the MDE
Part A: I expect more than half of the catastrophes to be sidechain-involving (the builder places
Trp, Arg, Lys, Phe with one fixed chi1 into a backbone it did not see). Part B: the fraction above
1e4 falls, the ranking does not improve (in-band rho stays under 0.15), and the reject stays null.
Plausibility 0.6 for A, 0.15 for any accuracy effect in B.

## Cost
Part A: included in the reject census, ~3 min CPU. Part B: 126 targets x 75 members x up to 3n
single points at ~9 ms = 126 x 75 x 40 x 9 ms = 57 min, plus 126 builder constructions at 2.8 s;
~1.2 h under `--tag AMBER`, ~0.8 GB, checkpointed per target. Agent-hours: 3.

## Information value
Part A sharpens Part IV of the report whatever the answer (the singularity is named as the
builder's or the pool's, with a number). Part B is the only way to test whether "AMBER measures
its worst clash" is a property of the force field on this pool or of the sidechain placement,
which decides whether any future all-atom scoring of retrieved windows is worth attempting.

---
Part A measured (2026-09-13 08:41, ledger L23, `s26/results/ph_reject_census.json` `singularity`): of
the 5,057 top-75 members above 1e4 kcal/mol, the closest heavy-atom contact is bb-sc on 2,627
and sc-sc on 2,267 (96.8% side-chain-involving) and bb-bb on 163; rho(e, min heavy-atom
distance) within the top-75 is -0.74. Part A survives. Part B (the relief itself, ~1.2 h AMBER)
is offered to the tournament as written.
