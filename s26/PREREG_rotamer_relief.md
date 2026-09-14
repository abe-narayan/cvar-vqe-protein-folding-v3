# PREREG -- ROTAMER RELIEF OF THE AMBER SINGLE POINT (lane PH, S26; own idea, part B runs only if ranked)

Written 2026-09-13 09:20, before `s26/ph_relief.py` produced a number. Never edited after;
addenda appended. Idea: `s26/IDEA_rotamer_relief.md` (part A measured, L23, survives). Code:
`s26/ph_relief.py`. Results: `s26/results/ph_relief_{run,report}.json`, cells under
`s26/results/ph_relief_cells/`.

## 1. What is measured and why

L23: within a shipped top-75 the AMBER single point is the minimum heavy-atom distance of the
rebuild (rho -0.74), and 96.8% of the members above 1e4 kcal/mol owe their closest contact to a
side-chain atom placed by the deterministic builder (`core.amber.build_sidechain`: chi1 from
`CHI_ANGLES`, one fixed value per residue type, no scan; proline ignores the override; glycine
and alanine have no chi1). The single point therefore measures the builder's rotamer choice on
most candidates. Part B asks what the energy measures once that choice is relieved.

The relieved energy, native-free, per candidate: one greedy sweep in residue order over every
residue with a chi1 (all but G, A, P); for residue i the chi1 option set is {the builder's
default, 60, 180, 300} degrees (duplicates removed; the default is always included, so the
relieved energy can never exceed the raw one); the option with the lowest single-point energy is
fixed and the sweep continues. At most 4 x n_sc single points per candidate, each a genuine
ff14SB/GBn2 evaluation of the built structure with NO minimisation (`refine_coords(k_restraint=
0, steps=-1)`, the same call the cache was built with). Backbone atoms never move. `E_relief`
and the chosen chi1 vector are stored per member.

## 2. Gates, before any RMSD

G1 (bit-exactness): on four members of every target the single point with `chi1=None` must
equal the cached `s24/cache_amber` value with relative difference 0.0; a non-zero difference
stops the run (a gate that cannot fire is not evidence: the comparison count is printed).
G2 (monotone): `E_relief <= e_raw` on every member (the default is in the option set).
G3 (the probe): `python s26/jobrun.py --agent PH --tag AMBER --name ph_relief_probe --est-ram
0.8 -- python s26/ph_relief.py probe --pdb 1A13 --m 5`: five members relieved, single points
per member, wall per member and peak RSS from `s26/jobs_done/ph_relief_probe.json`; the
126-target cost is extrapolated from it and stated in an addendum before the run.

## 3. Part B(i), the census question (native-free, printed by the run itself)

The fraction of top-75 members above 1e4 kcal/mol before (0.535, L23) and after relief; the
fraction above 1e6; the change in the closest-contact class. FALSIFIER B(i): if the fraction
above 1e4 after relief is not below HALF the raw fraction (i.e. not below 0.27), the
singularity is not the builder's and the idea is closed at this step; B(ii) and B(iii) are still
reported but carry no claim.

## 4. Part B(ii), does the relieved energy rank? (gated, ORACLE evaluation)

Per target, on the top-75: Spearman(E_relief, ORACLE candidate RMSD `rr`) whole set and in-band
(members within 3 A of the best, as S25 defined it), with an Rg partial (Rg native-free from the
CA windows). The same for the raw energy, on the same members. Contrast: per-target rho
difference, relieved minus raw, `ST.compare` with folds. Claim only if the difference clears its
own MDE with the fold CI excluding zero. Expected: no measurable improvement (in-band rho stays
below 0.15).

## 5. Part B(iii), the reject on the relieved energy (gated)

Arm S of `s26/PREREG_amber_reject.md` (reject from the shipped top-75, no refill) on E_relief at
the same four thresholds, primary 1e4, with the identical controls (RANDS matched count, 16
draws; PERMS, 16 permutations; anchor). Arm R (refill from the ranked pool) needs the relieved
energy of candidates beyond the top-75, which costs 6.7x more; it is NOT run unless the
coordinator grants a top-150 budget (2x), in which case refill is limited to the relieved prefix
and shortfalls are counted. Endpoints: point cloud for everything; built chain (`I.project`) for
S at 1e4 and RANDS at 1e4 with 4 draws (about 9 projections per target, 1.4 h CPU). Falsifier as
in the reject prereg: S beats the anchor and RANDS past their own MDEs with fold CIs excluding
zero, 5/5 folds, on the built chain. Expected: null or harmful (S23 L5: removing members costs).

## 6. Cost, memory, checkpoints

Per member about 4 x 11 = 44 single points at ~9 ms plus the Python around each call (S13
measured 8.66 ms bare, 12.5 ms with components; `refine_coords` adds ~28 ms of Python per call,
so the fast path `AmberSP`-style, `_assemble` + `getState`, is used and gated bit-exact against
`refine_coords` on every target's first member). Estimate 0.5 s per member, 40 s per target,
126 x 40 s = 1.4 h plus 126 x 2.8 s builder construction, one process, est-ram 0.8 GB,
`--tag AMBER`, cells `s26/results/ph_relief_cells/<pdb>.json` after each target (75 energies, 75
chi1 vectors, single-point counts, walls), resumable. The probe fixes the real number before the
run. Report: 5 min CPU, plus 1.4 h CPU for the built-chain half of B(iii). Agent-hours: 3.

## 7. Replication

A positive B(ii) or B(iii) result is re-run with the residue sweep order reversed (a different
greedy path) and the fold order reversed, and must land inside its own fold CI.

---
## ADDENDUM 1 (2026-09-14 00:45) -- THE FIRST PROBE FIRED GATE G1, AND WHY

`ph_relief_probe` (jobs_done, exit 1, 10 s, peak 0.264 GB): on 1A13 member 387 the fast single
point equalled `refine_coords(k=0, steps=-1)` to 0.0 but differed from the cached value by a
relative 0.47. Cause, read from `s24/d2_amberscore.py:114`: the cache was built with
`s13.qarch_lib.Space(pdb, 4).rep` (k = 4 torsion states), while `ph_relief.py` built its
representation with k = 8. `core.amber.builder_for` keys its OpenMM context on `rep.n_states`
and calibrates the hydrogen frames on the representation's reference states, so the two
representations give two legitimate single points that differ by the hydrogen placement. The
fix is one constant (k = 4, the cache's instrument); the gate is unchanged and must now pass at
0.0 against both `refine_coords` and the cache, or the run stops. Nothing else changed.
