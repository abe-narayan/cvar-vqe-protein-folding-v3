"""s26/ph_append_census.py -- one-shot writer for lane PH's census ledger entries (2026-09-13).

Re-reads the ledger tail immediately before appending (L16b rule 4), takes the next three free
numbers, appends the entries, adds the STATUS lines, the three PREREG addenda and the IDEA
Part-A note, and replaces the placeholders LXX-CIS / LXX-REJ / LXX-C3 in the notes and findings
with the numbers actually taken.  Idempotent guard: refuses if the ledger already carries the
CIS CENSUS heading.
"""
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
NOW = time.strftime("%Y-%m-%d %H:%M")


def rd(p):
    with open(p, encoding="utf-8") as fh:
        return fh.read()


def wr(p, s):
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(s)


led = os.path.join(HERE, "LEDGER.md")
txt = rd(led)
if "CIS CENSUS: NO CIS PEPTIDE BOND" in txt:
    print("already appended; nothing done")
    sys.exit(0)
nums = [int(m) for m in re.findall(r"^## L(\d+)", txt, flags=re.M)]
n = max(nums) + 1
CIS, REJ, C3 = n, n + 1, n + 2

E_CIS = f"""## L{CIS} -- CIS CENSUS: NO CIS PEPTIDE BOND EXISTS ANYWHERE ON THE INSTRUMENT; THE DATABASE STEP GATE, NOT THE PROJECTION, SETS THE COST (2026-09-13, PH)

`s26/ph_cis.py census`, `s26/results/ph_cis_census.json` (complete 126/126), job
`s26/jobs_done/ph_cis_census.json` (exit 0, 90 s, peak RSS 0.038 GB). Pre-registered in
`s26/PREREG_cis.md` section 2 with the prediction "zero on model 1 and in the pool, because of
the step gate; some cis bonds in other ensemble models". Allowed before the gate by L5: ORACLE
DIAGNOSTIC on the natives (omega angles only, no RMSD), native-free on the pool.

    model-1 natives with a cis bond, |omega| < 30 deg        0 / 126
    the same by consecutive CA-CA < 3.3 A                    0 / 126   (the two criteria agree 126/126)
    deposited models with a cis bond, every ensemble         0 / 1,966
    universe windows with a CA-CA step < 3.3 A               0 / 2,352,893   minimum step 3.5045 A
    K=500 pools and production top-75 with such a window     0 and 0
    omega non-planarity |180 - |omega||, 1,507 bonds         mean 1.91 deg, median 0.34, p90 5.8, p99 17.4,
                                                             max 42.8 (9UV5, bonds 0 and 6: -137.2 and +144.8 deg);
                                                             0.53% of bonds beyond 20 deg, 0.13% beyond 30 deg

The universe minimum step IS the gate: `core/data.py:406-407` drops any peptide with a
consecutive CA-CA step below 3.5 A and `core/data.py:697-698` drops any fragment window with
one. The 126 are drawn from that database (`s7/debias.py:103-120`) and so is the sealed
benchmark, so neither can contain a cis target and no pool can contain a cis window. The
two-bond-length projection is therefore worth exactly 0.000 A on this instrument by
construction of the database, not by any property of the projection. The residual cost of the
constant omega here is the non-planarity tail (0.5% of bonds beyond 20 deg), which part 2
(gated) prices as the ideal-trans floor on the native's own torsions. The registered prediction
that other ensemble models carry cis bonds was WRONG: 0 of 1,966. The cis-peptide question is a
world-supply question (the 16 containment-fresh targets, 10 amyloid) and the design note in
`s26/agentPH_FINDINGS.md` section 1.4 is written for that supply.

---
"""

E_REJ = f"""## L{REJ} -- STERIC REJECT CENSUS: AT 1e4 kcal/mol THE REJECT REMOVES 40 OF 75, EMPTIES 11 TOP-75 SETS AND 8 WHOLE POOLS, AND 96.8% OF THE CATASTROPHES ARE SIDE-CHAIN CONTACTS (2026-09-13, PH)

`s26/ph_reject.py census`, `s26/results/ph_reject_census.json` (complete 126/126), job
`s26/jobs_done/ph_reject_census.json` (exit 0, 100 s, peak RSS 0.117 GB). Native-free: the
63,000 cached single points (`s24/cache_amber`, pool identity `universe_idx == I.pool_idx` and
production `sub` == score top-75 asserted on every target) and the ideal-geometry rebuilds; no
RMSD, no native. Required by `s26/briefs/PH.md` section 3.1 before any RMSD is read; the prior
was registered in `s26/PREREG_amber_reject.md` section 5.

    threshold     pool frac > T   top-75 rejected   zero-reject   all-75 rejected   no survivor in 500   refill depth (mean/median)   R overlap with anchor
    1e3             0.760           54.5 / 75          1              25                 20                270 / 234                    0.27
    1e4 PRIMARY     0.586           40.1 / 75          2              11                  8                201 / 147                    0.46
    1e5             0.446           30.3 / 75          6               3                  1                167 / 118                    0.60
    1e6             0.345           23.3 / 75          8               2                  0                135 / 103                    0.69

The 0.586 reproduces S25's 58.6% (`s25/results/phys_landscape.json`). At the primary threshold
the operator is not a small surgical reject: it removes more than half of every shipped set,
reaches rank 147 (median) of 500 to refill, empties the whole top-75 on 11 targets (1G89 1ID6
1LB7 2MAI 2NB7 2XL1 5MML 5Z5W 7BX2 8UN8 9S5G) and finds no survivor among all 500 candidates on
8 (1G89 1ID6 2MAI 2NB7 2XL1 5MML 7BX2 8UN8). Both empty cases fall back to the anchor, a rule
added in the prereg's addendum 1 before any RMSD is read. Native-free geometry of the retained
set at 1e4: R is +0.254 A more expanded in Rg than the anchor (SE 0.058) and +0.144 A in the
minimum |i-j| >= 3 CA-CA distance; S is +0.060 and +0.070. Same sign as S25's +1.10 A expansion
of AMBER's top-75, smaller because the distogram's order is kept among the survivors.

WHERE THE SINGULARITY LIVES (the optional Part IV measurement, `singularity` block; every top-75
member rebuilt with all heavy atoms through the reference builder `sidechains.py`, no OpenMM):

    closest heavy-atom contact, residue separation >= 2
      all 9,450 members:              bb-bb 1,054   bb-sc 5,300   sc-sc 3,096
      the 5,057 members above 1e4:    bb-bb   163   bb-sc 2,627   sc-sc 2,267    -> 96.8% side-chain-involving
    members per target with a heavy-atom pair closer than 2.0 A:   all-atom 40.7 of 75 (SE 1.9)
                                                                    backbone+CB 2.61 of 75 (SE 0.34)   [S19 section 3.1: 2.66, reproduced]
    Spearman(e_amber, minimum heavy-atom distance) within a top-75: mean -0.743, median -0.803
    minimum heavy-atom distance: members above 1e4, 1.48 A; members below 1e4, 2.31 A

The AMBER single point on the top-75 is, to rho -0.74, the minimum heavy-atom distance of the
rebuild, and 96.8% of the rebuilds it condemns are condemned by a contact involving a side chain
placed by the deterministic builder (fixed chi1, no rotamer scan). The physically impossible
class on the backbone is 2.6 per 75 (S19); the class the 1e4 reject removes is 40 per 75.
Stated before any RMSD is read: the steric reject at the primary threshold is a reject of the
builder's side-chain placement, not of the pool's backbones. This is Part A of
`s26/IDEA_rotamer_relief.md` and it survives.

---
"""

E_C3 = f"""## L{C3} -- C3 NATIVE-FREE PART: THE PRODUCTION RELAXATION MOVES THE CA TRACE 0.220 A RMS; 58.7% OF BUILT CHAINS START ABOVE 1e4 kcal/mol; 125 OF 126 CONVERGE (2026-09-13, PH)

`s26/ph_c3.py nativefree`, `s26/results/ph_c3_nativefree.json` (complete 126/126), job
`s26/jobs_done/ph_c3_nativefree.json` (exit 0, 5 s; the RSS sampler polls every 5 s and
under-reads a 5 s process, so no memory number is claimed). Read from
`bench_results/cache/1fc9f2dcf489e2fb` with every RMSD key stripped; nothing here reads a native.

    AMBER displacement of the built chain, per-atom RMS after superposition   0.220 A (SE 0.008, median 0.197, range 0.103 to 0.591)
    `amber_moved` (restraint RMSD on N/CA/C, unsuperposed)                     0.233 A
    the built chain's own AMBER energy e0 before relaxation                    median 8.6e4 kcal/mol, min -473, max 1.3e14; 58.7% above 1e4
    relaxed energy e1                                                          mean -560 (SE 30), max +1262 (9KAR)
    converged (e1 <= CONVERGE_MAX_KCAL = 1000)                                 125 / 126 (9KAR fails)
    bond + angle strain after                                                  mean 60 kcal/mol
    virtual CA-CA bond: built 3.80395 (exact, sd 1e-16) -> relaxed             mean 3.867, min 3.12 (1M02), max 5.38 (2BP4, e1 +845); 4.86 on 9KAR
    targets with a relaxed CA-CA outside [3.6, 4.0]                            6 / 126
    Rg change                                                                  +0.046 A

Two things the validity story must carry: the emission's own strain census (58.7% of built
chains sit above 1e4 kcal/mol before relaxation, the same fraction as the pool), and the two
targets where the relaxation itself breaks a virtual bond (2BP4, 9KAR). Derived prediction for
stage 1, registered in `s26/PREREG_c3_control.md` addendum 1 before the gate: an orthogonal move
of 0.220 A on a 3.21 A chain costs about m^2 / (2 RMSD) = 0.0075 A by the S16 identity; the
production step costs +0.0207; so AMBER is predicted WORSE than the matched random control by
about +0.013 A. Measured after the gate.

---
"""

wr(led, txt.rstrip("\n") + "\n\n" + E_CIS + "\n" + E_REJ + "\n" + E_C3)
print("ledger", CIS, REJ, C3)

# ---------------- STATUS
st = os.path.join(HERE, "STATUS.md")
lines = rd(st).split("\n")
for i, l in enumerate(lines):
    if l.startswith("## PH"):
        lines[i + 1:i + 1] = [
            f"- {NOW} running: nothing (gate closed). Done: reading list; PREREG_amber_reject/cis/c3_control; "
            f"five IDEA files; ph_lib/ph_cis/ph_reject/ph_c3 + synthetic tests; census jobs ph_cis_census "
            f"(0.038 GB), ph_reject_census (0.117 GB), ph_c3_nativefree; ledger L{CIS}-L{C3}; findings and Part IV notes; commit 5dc7a3a6.",
            f"- {NOW} next: READY TO LAUNCH the minute PHASE 0 SIGNED OFF lands: `ph_cis.py floor` (1 min), "
            f"`ph_c3.py stage1` (2 min), `ph_reject.py cloud` (5 min) then `ph_reject.py chain` (CPU, ~3 h, per-target cells); "
            f"C3 stage 2 when P delivers (AMBER, probe first). Awaiting tournament ranking for IDEA_branch_select / rotamer_relief part B.",
        ]
        break
wr(st, "\n".join(lines))

# ---------------- PREREG addenda and IDEA note
def append(name, text):
    with open(os.path.join(HERE, name), "a", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


append("PREREG_amber_reject.md", f"""
---
## ADDENDUM 1 ({NOW}) -- THE CENSUS, WRITTEN BEFORE ANY RMSD IS READ

`s26/results/ph_reject_census.json`, ledger L{REJ}. Measured what section 5 guessed: at 1e4 the
reject removes 40.1 of 75 (53.5%; prior "about half"), refill reaches rank 147 (median; prior
150 to 200), 11 targets lose the whole top-75 and 8 have no survivor in the whole pool (the prior
did not anticipate whole-pool emptiness). At 1e3 it is 54.5 of 75 with 25 empty sets and 20
empty pools; at 1e5 30.3 of 75; at 1e6 23.3 of 75.

Rules added before any endpoint is read: (1) an EMPTY retained set, S or R, falls back to the
anchor (do nothing), the only native-free deployable choice; the number of fallbacks per arm and
threshold, and the number of targets actually moved, are reported beside every MDE. (2) Every
contrast is shown on all 126 and on the moved subset. (3) The operator at the primary threshold
is now known to be a reject of the builder's side-chain placement (96.8% of condemned members
have a side-chain-involving closest contact; rho(e, min heavy-atom distance) = -0.74) and not of
the pool's backbones (2.6 of 75 have a backbone+CB contact below 2.0 A; S19 reproduced). That
does not change the falsifier; it changes what a positive result would mean, and it is stated
now so it cannot be discovered afterwards. Prior unchanged: null or harmful.

The census artefact was stamped by `ph_reject.py` at source sha ba028f7717d38725, before the
fallback lines (`R_empty`, `R_eff`) were added to `retained_sets`; the census computation does
not read those fields. The committed file (5dc7a3a6) is the one the endpoint runs will stamp.
""")

append("PREREG_cis.md", f"""
---
## ADDENDUM 1 ({NOW}) -- THE CENSUS RESULT AND THE RE-SCOPE OF PART 2

`s26/results/ph_cis_census.json`, ledger L{CIS}. The registered prediction held on model 1 and on
the pool (0/126, 0/2,352,893, universe minimum step 3.5045 A, which is the gate) and FAILED on
the ensembles: 0 of 1,966 deposited models carry a cis bond, where I predicted some. Written
before part 2 runs: the cis-vs-non-cis comparison is empty; part 2 measures the ideal-trans
representation floor (CA and N/CA/C bases) on all 126 as the residual cost of the constant
omega, its Spearman with the per-target maximum omega deviation, and reports 9UV5 (the one
target with a bond beyond 30 deg) by name, never as a subgroup. The falsifier for "the design is
worth writing on this instrument" cannot fire, so the design is written for the world supply
with the base rate stated as 0 of 1,507 bonds here and a literature value with its source.
""")

append("PREREG_c3_control.md", f"""
---
## ADDENDUM 1 ({NOW}) -- THE MEASURED MAGNITUDE AND THE DERIVED PREDICTION

`s26/results/ph_c3_nativefree.json`, ledger L{C3}. The production displacement is 0.220 A RMS
after superposition (median 0.197, range 0.103 to 0.591), not the ~0.15 A guessed from one
target in section 3. The S16 identity then predicts that a move of that size orthogonal to the
residual costs about 0.220^2 / (2 x 3.21) = 0.0075 A on the built chain; the production step
costs +0.0207. Registered prediction: AMBER minus the matched random control = about +0.013 A
(AMBER worse), ORACLE cosine of AMBER's displacement negative. 9KAR does not converge (e1 = 1262
> 1000); the convergence gate is REPORTED, never silently applied: every stage-1 contrast is
given on all 126 and on the 125 converged. 2BP4 (relaxed CA-CA 5.38 A) and 9KAR (4.86 A) are
the two emissions whose virtual bond the relaxation breaks; they stay in.
""")

append("IDEA_rotamer_relief.md", f"""
---
Part A measured ({NOW}, ledger L{REJ}, `s26/results/ph_reject_census.json` `singularity`): of
the 5,057 top-75 members above 1e4 kcal/mol, the closest heavy-atom contact is bb-sc on 2,627
and sc-sc on 2,267 (96.8% side-chain-involving) and bb-bb on 163; rho(e, min heavy-atom
distance) within the top-75 is -0.74. Part A survives. Part B (the relief itself, ~1.2 h AMBER)
is offered to the tournament as written.
""")

# ---------------- placeholders
for name in ("PH_PART_IV_NOTES.md", "agentPH_FINDINGS.md"):
    p = os.path.join(HERE, name)
    if os.path.exists(p):
        s = rd(p)
        s2 = s.replace("LXX-CIS", f"L{CIS}").replace("LXX-REJ", f"L{REJ}").replace("LXX-C3", f"L{C3}")
        if s2 != s:
            wr(p, s2)
            print("patched", name)
print("done")
