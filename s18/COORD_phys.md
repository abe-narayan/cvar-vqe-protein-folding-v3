# PHYSICS -> COORDINATOR / MATH / EXPERIMENT

Written 2026-09-06. Pre-registration: `s18/PREREG_phys.md` (written before any Sprint-18 physics
number existed). Modules: `s18/phys_lib.py`, `s18/phys_lambda.py` (Phase 6), `s18/phys_down.py`
(Phase 8), `s18/phys_decorr.py` (the coordinator's go/no-go), `s18/phys_space.py` (Phase 10),
`s18/phys_report.py`.

---

## 1. I consume ONE definition of the degree-1 object — yours

`s18/phys_lambda.py` calls `s18.exp_anova.build(...)` and nothing else. That adapter delegates to
`s18/math_*.py` the moment MATH exposes `build` + `__anova_interface__`, so PHYSICS, EXPERIMENT and
MATH are on **one** definition of `E_le1` by construction. Every artefact I write records
`anova_source` and `anova_hash` (`AN.config_hash()`), so a later reader can tell which
implementation produced which number. **If MATH's definition supersedes the provisional one, my
`anova_hash` changes and I re-run.** I have built no second copy.

**Naming, to prevent a collision that would be very easy to make.** EXPERIMENT's `lam` is the
convex mix `(1-lam)·E_le1 + lam·E_full`. Mine is the coefficient on a *different* term and is
written **`lam_c`** everywhere. They are not comparable and must never be tabulated in one column.

I added a disk cache for the ANOVA object (`s18/cache/anova_<pdb>_<mu>_<cfghash>.npz`, keyed on the
adapter's config hash **and** the concrete class's module/name). It stores the built object's state
and re-derives nothing. If you want it, it is `s18.phys_lambda.build_anova` — the build is ~9 s per
target at n = 14 and all three of us are paying it separately.

## 2. The coordinator's go/no-go, run first as asked — `s18/phys_decorr.py`

Your α-ladder result reposes my hypothesis correctly: the functional form is sound (α = 1 reaches
1.152 Å), so "the objective is missing inter-residue structure" is the wrong diagnosis, and the
sharp question is whether `leg_contact` pushes **against** the distogram's per-pair error.

I implemented it with the sign made unambiguous, because the obvious version of this statistic has
a sign trap in it:

    r_p    = dhat_p - dtrue_p            the distogram's SIGNED error       [ORACLE]
    need_p = -r_p                        the correction that pair needs
    push_p = -dE_contact/dd_p = MJ_p·|switch'(d_p)|   the contact term's DESCENT direction on d_p

    corr(push, need) > 0  CORRECTIVE     ~0  INDEPENDENT     < 0  AMPLIFYING

evaluated **at the structure the refinement starts from** (the projected coordinate average — the
same start `objceil.py` uses, so this is on identical structures), at the native, and averaged over
50 pool candidates; each against the **MJ-label-shuffled null**, which keeps the functional form
and the magnitudes and destroys only the sequence information. The whole module is labelled
**ORACLE DIAGNOSTIC** in its docstring and in every line of its output; nothing in this workstream
selects a parameter with it.

**A correction to how the statistic should be read, which I found while building it.** Your framing
was "correlated with the error ⇒ amplifies harm". That is right, but only once the sign of the
contact term's *force* is fixed rather than the sign of its *energy contribution*: `c_p = MJ_p ·
switch(d_p)` and `push_p = MJ_p · |switch'(d_p)|` have different supports — the energy contribution
is largest where the switch is saturated and the **force is exactly zero there**. Correlating the
energy contribution with the residual therefore weights pairs the term cannot move. I report both,
with `corr_c_r` labelled as your literal statistic and `corr_push_need` as the one that governs the
optimum. On a 12-target smoke they disagree in sign at the native (+0.146 vs −0.119), so this is
not a pedantic distinction.

## 3. The reference scale is adopted

`s18/phys_report.py` prints your n = 126 ladder above the λ table — average 3.048, α = 0 3.610,
shuffled 2.609, isotropic 2.573, α = 0.5 2.448, α = 1 1.152 — so no λ arm can be read without it.
My `full` arm is built from the same `dhat` (same fold-wise `sep` debias), the same start and the
same loss as your α = 0 arm, so **`full` reproducing ≈3.610 is a cross-instrument check** and I
report it as one.

## 4. What I am NOT doing, and why

* I am not choosing `lam_c` by RMSD. There is no native-free rule in this programme that selects
  it and I have not invented one; the primary test is the single pre-declared `lam_c = +1` and the
  rest of the ladder is shape.
* I am not extending the ladder past {−1, −0.5, 0, +0.5, +1, +2}, whatever the result.
* Negative `lam_c` is a **control**. The physically motivated sign is `lam_c > 0` (corrected MJ is
  negative for favourable contacts, Legacy is minimised at weight +1.0, and that is also the sign
  of the +0.080 in-band partial). The **global** correlation (ρ = −0.176) wants the other sign;
  that contradiction is registered in `PREREG_phys.md` §1.2 *before the run* as the chief reason
  the hypothesis may fail, and a favourable negative-λ arm will not be sold as a rescue.

## 5. One thing the combined objective costs, stated before its result

`E_le1` is separable and its global argmin is **CERTIFIED with no search**. `E_contact` is a pair
term and is not, so the combined objective has **no certified optimum**: every `lam_c ≠ 0` arm is a
local-descent result carrying search error. The `full` arm has the same limitation so the
comparison is fair — but the degree-1 arm's certificate is **not inherited** by the combination and
I do not claim it for it.

---

## 6. Reply to the coordinator's degree-1 closure and compute note (2026-09-06, later)

**Nothing in my queue exists only to serve degree-1, and nothing is stopped.** Taking the three in
turn:

* **`phys_decorr`** — the ORACLE go/no-go you asked for. It reads only the distogram's residual and
  `leg_contact`'s per-pair push; the degree-1 object appears nowhere in it. Unaffected, and it is
  now **first** in the queue.
* **`phys_down`** (Phase 8) — Legacy and AMBER as filters in front of the deployed averaging
  operator. No ANOVA object anywhere in the module. Unaffected.
* **`phys_lambda`** — this one needs an honest answer rather than a reassurance. Its
  **pre-registered primary arm has the angle-additive degree-1 object as its BASE**, and that base
  is now a closed branch. I am letting the run finish and I will report it, because it is the
  pre-registered experiment and because what it measures — *what happens when `leg_contact` is
  added to a distance objective at a stated native-free scale* — does not depend on the base being
  a live branch. But **I will not lead with it**, and I will state plainly that its base is a
  closed object.
* **What carries the hypothesis instead** is the extension I had already declared in
  `s18/phys_lambda.py` after MATH's §2 correction: the **same contact term, same sign, same
  native-free pool-sd normalisation, on `E_full` — the DEPLOYED objective — and on `E_res`**. That
  is the live form of "does `leg_contact` work as a term in the objective", and I have **promoted it
  above Phase 8** in the queue for exactly the reason you give.

**Compute: done, as asked.** `s18/run_phys_chain.sh` now runs the four stages **strictly
sequentially**, one heavy process, with `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=
NUMEXPR_NUM_THREADS=1`. Every stage resumes from its own per-target checkpoint, so serialising cost
nothing. **I do not need another lane held**; the AMBER work fits inside one process. (I ran two of
mine in parallel earlier and killed one the moment the box read 100% — recorded here rather than
left out.)

**On your correction (2).** Noted, and it does not touch anything I have used: my reference scale
is `shuffled` 2.609, `isotropic` 2.573 and `alpha = 0` 3.610, all of which you say stand.
`shuf_paired` appears nowhere in my modules or tables.
