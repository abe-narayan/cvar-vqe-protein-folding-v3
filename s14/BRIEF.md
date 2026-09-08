# SPRINT 14 — VQE-CENTRED FOLDING REDESIGN. Shared agent brief.

Read this completely before writing any code. It is the contract every agent works under.

> **CORRECTIONS BANNER, added 2026-09-05 late in the sprint.** Several statements in section
> 2 below were carried in from Sprints 12-13 and have since been **corrected or refuted by
> this sprint's own measurements**. They are left in place because this project never deletes
> a superseded claim. Before relying on any of them, check `s14/LEDGER.md` and section L of
> `s14/SPRINT14_DOSSIER.md`. The ones that moved:
>
> * "Coverage, not accuracy, is the make-or-break parameter" and the >= 90% coverage gate:
>   **REFUTED**. Measured gaps cost 0.014 A at sigma 12. The binding quantity is the fraction
>   of TARGETS with heteronuclear shifts (0.444), not residues within them (0.952).
> * The 1.486 A restraint figure: **superseded**. Under TALOS-N's real error mixture the
>   channel lands at 2.347 A, and the whole route is closed by arithmetic (54/126 targets;
>   ORACLE-perfect torsions still leave the instrument at 2.021 A, and 55 are needed).
> * "The incumbent is equivalent to sigma ~29 deg": refined to **27.1 deg** on an i.i.d.
>   error surface, and error COHERENCE is worth nearly a factor of two in required accuracy
>   (10.6 to 20.2 deg to reach 2.0 A) though all real predictors are near-i.i.d.
> * "Optimising harder makes structures worse" and the certified optimum being worse than
>   random: **both are properties of the PHYSICAL ENERGIES, not of the problem**. A
>   native-free structural objective has a certified optimum 0.885 A BETTER than random.
> * The three-inert-torsion / `log2 k` dead-qubit count and the blanket QNG refutation are
>   corrected in place further down this section.

---

## 0. THE QUESTION

> Can realistic, legitimately-available local structural information make the torsion
> search space informative enough that **VQE/CVaR itself becomes a useful conformational
> optimiser**, rather than a quantum decoration bolted onto a classical pipeline?

Target: **mean CA-RMSD < 2.0 A on the 126-target development instrument**, predictive,
with a demonstrable causal contribution from VQE/CVaR.

Four things are NOT disposable: **genuine VQE**, **genuine CVaR**, **a genuine
non-all-atom (Legacy) energy**, **a genuine all-atom AMBER energy**. Everything else
— the retrieval pipeline, the encoding, the ansatz, the objective — may be replaced.

---

## 1. HARD RULES. Violating any of these invalidates the work.

### 1.1 The benchmark is protected

`results/benchmark_manifest.json`, `peptide_folds.json`, `peptide_clusters.json` are
**pinned**. Do not read benchmark structures for design. Do not inspect benchmark RMSDs.
Do not tune on it, select on it, or regenerate any of it. It is evaluated once, by the
coordinator, after the architecture is frozen. If you are not sure whether something
touches it, ask the coordinator rather than guessing.

### 1.2 NO NATIVE INFORMATION AT INFERENCE. Absolute.

Native coordinates, native torsions, native distances, native RMSD, native energies,
native contact maps and native secondary structure may be used **only** for:

  * training supervised predictors under leave-fold-out discipline,
  * experiments explicitly labelled **ORACLE DIAGNOSTIC**,
  * post-hoc evaluation.

The Sprint 12 figure of ~1.486 A from torsion restraints is an **ORACLE feasibility
ceiling**, not a predictive result. Never quote it as performance. Never recreate that
leakage by accident — the classic route is a helper that silently calls
`Space.ORACLE_*` or `I.load_univ(pdb)["nat_ca"]` inside something that looks predictive.

Every experiment must state, in one line, exactly what information the model had.

### 1.3 Splits

`tuning126` is the development instrument. Folds and identity clusters are pinned; five
leave-fold-out folds. `dev24` is the cluster-disjoint dev split and must not be quietly
converted into a tuning target. Do not repeatedly optimise the same aggregate number
until it becomes a benchmark surrogate — reserve subgroups for confirmation.

### 1.4 No fake components

Do not call a classical search "VQE". Do not call enumeration "VQE". Do not call a
surrogate "AMBER". Do not label a classical result quantum. If you run a classical
control, say so in the arm name.

### 1.5 Preserve history

Never delete or silently overwrite an earlier finding. If you correct one, write the
correction next to it with the evidence. Corrections are the most valuable output this
project produces.

### 1.6 Report honestly

If VQE loses, say so and diagnose why. If CVaR does nothing, say so. If AMBER hurts,
say so. Do not manufacture a positive. A rigorous negative with a mechanism is a
success under this brief (section 57 of the sprint brief).

---

## 2. WHAT IS ALREADY ESTABLISHED. Do not re-run these; build on them.

These come from Sprints 5-13 and are evidence, not doctrine — but re-testing one needs
a specific new hypothesis, not a rerun.

**Accuracy state.** Shipped retrieval pipeline emits **3.213 A** on tuning126 (this is
`synthesis_fit` = 3.2040761603809194 through the shared instrument; 3.4540 is the older
`shipped` arm). No validated improvement has ever transferred to the held-out benchmark:
the last attempt was +0.0103 A, CI [-0.1596, +0.1803].

**The torsion space is not the barrier.** ORACLE coordinate descent on true CA-RMSD in
`torsion_lib2` k-state space:

| k | mean live qubits | ORACLE descent | frac < 2 A |
|---|---|---|---|
| 4 | 25.9 | 1.594 | 0.73 |
| 8 | 38.9 | 1.184 | 0.91 |
| 16 | 51.8 | 0.876 | 0.95 |

From a **random** start the same descent reaches only **1.982** — 0.388 A of that ceiling
is the privileged oracle start. Uniform-random torsion states reach 2.698, already better
than the 3.213 pipeline. Of the 1.104 A the sequence-conditioned library buys over an
information-free space, **88% is generic Ramachandran and 12% is sequence-related**; a
*wrong target's* library costs only +0.135 A.

**Neither energy ranks the native.** Legacy rank correlation with RMSD **inside the
low-energy decile is +0.043**; raw AMBER -0.088. The native sits at the **32nd-40th
percentile** of both. A monotone log compression leaves rho unchanged while collapsing
AMBER's range from 16.1 decades to 1.8 — **the damage is in the ordering, not the scale**,
so no rescaling or softening fixes it. On nine **fully enumerated** targets (262,144
configurations each) Legacy's **certified global optimum is +0.139 A worse than random
sampling**. Legacy's entire skill is clash rejection: remove its steric term and it
inverts to +0.395 A worse than random.

**Optimising harder makes structures worse, and the trap is shaped to fool you.**
RMSD vs evaluation budget on the enumerated targets: 3.764 at 10 evaluations, 3.667 at
300, **3.920 at the certified global optimum**. The turn is governed by the *fraction of
space explored*, not the evaluation count. **A VQE at n=12-16 with 1e4 evaluations sits
on the improving limb of a curve whose limit is known to be bad. It will look like it is
working.** Any VQE arm must be checked against this.

**Weak controls to avoid.** A zero-information **constant alpha-helix** (phi=-63 deg,
psi=-42 deg) emits 4.065 A and **beats uniform random sampling by 0.457 A**
[-0.822, -0.090]. "Beats random" proves nothing here. The constant helix is the baseline
to clear. Likewise the 1-local torsion prior's apparent advantage (3.969) is entirely a
helix artefact: rho(helix fraction, RMSD) = -0.744, and on non-helical targets its
advantage vanishes.

**Energy reached and structural accuracy are uncorrelated.** Within-target
rho(energy reached, RMSD) across six optimisers is **-0.006**. SPSA optimises the AMBER
objective best of five arms and returns the **worst** structure. Always report objective
quality and structural quality on separate axes.

**Sequence carries almost no phi.** Leave-fold-out mean absolute torsion error, degrees:
sequence-blind marginal 36.4 / 72.8; residue-class prior 36.9 / 68.9; full sequence
context and properties **36.1 / 62.4**. The entire measurable sequence channel is 10.4 deg
of psi. The best sequence-only builder achieves sigma 67.7 deg and emits 3.770 A; a
*perfect* confidence gate still loses at 3.478. **The incumbent 3.213 A pipeline is
equivalent to sigma ~29 deg** — that is the number any torsion channel must beat to tie.

**Torsion restraints (ORACLE) reach the target.** sigma 12 deg at full coverage builds
1.486 A directly; as a pool *filter* only 2.972. Coverage, not accuracy, is the
make-or-break parameter: at sigma 12 deg, coverage 1.00/0.90/0.75/0.50 gives
1.486/2.141/2.792/3.290. **CORRECTION carried forward:** terminal dropout — which is
where a shift-based predictor actually declines — is **0.40-0.50 A cheaper** than
uniform dropout, so uniform dropout *understates* the channel.

**Exact locality theorem.** Under an ideal-geometry backbone builder, CA-CA distance
`d_ij` depends on **exactly the `j-i-1` residues strictly between i and j** (agreement
1.0000, zero counterexamples at machine precision, holds at termini, on GLY/PRO, under
cis-omega). All-atom: CA `i<m<j`, N `i<=m<j`, C/O `i<m<=j`, CB `i<=m<=j`. Consequence: a
separation-8 pair is a **14-qubit interaction at k=4**, so **no 2-local Ising form of a
distance-based molecular objective exists in this encoding**. Both energies are
full-register; "AMBER is less local than Legacy" is a **category error** and is refuted.

**Spectra and trainability.** Both energies are diagonal, hence exact weighted sums of
Pauli-Z strings. After monotone rank-preserving conditioning, mean Pauli weight is
Legacy 2.236 vs AMBER 3.015 (AMBER higher on 79/79 cells). The chain spectrum -> gradient
variance closes with **no free parameter** (measured/predicted 1.006 and 1.001).
**Beware:** the *unconditioned* spectrum is a delta-spike artefact — raw AMBER's top-10
configurations of 4,096 carry a median 99.6% of its Walsh variance and its weight lands
on Binomial(m, 1/2) exactly. 99th-percentile winsorisation is NOT enough; condition
monotonically.

**Trainability negatives.** Cost-locality does *not* explain trainability at 6-18 qubits
(the measured ansatz kernel is flat in Pauli weight; what matters is n). The metric contains
**no Hamiltonian** — bit-identical across energy models at matched parameters (0.000e+00),
now a unit test. The energy model selects which region of the manifold is visited.

**QNG — CORRECTED 2026-09-05. The refutation holds ONLY at depth 1.** "The Fubini-Study
metric is full rank with `g_ii = 0.2500` exactly and is exactly `I/4`" is true and complete at
depth 1 for **every** entangler pattern, and **does not extend to depth >= 2**. `g_ii` stays
0.250000 at every depth, but off-diagonal structure appears at depth 2 and the condition
number reaches **155.2** for a torsion-aware block entangler at depth 2 and **478.9** at depth
3 (chain 3.665, ring 9.937, all-to-all 2.478). **A problem-inspired ansatz at depth >= 2 does
reshape the manifold**, so QNG is open again in exactly the regime such an ansatz occupies.

**CVaR — three defects, only one previously known.** The value estimator is correct.
(1) The recorded `baseline="tail"` gradient defect is real and now has a CLOSED-FORM bias.
(2) NEW: `cvar_from_samples` averages the lowest `ceil(alpha*N)` samples where the correct
estimator splits the boundary atom, giving an **upward bias at non-integer `alpha*N`**
(+0.134 sd at N=13, alpha=0.10; decays as 1/N). At alpha=0.05 with 256 shots that is the same
order as the differences an optimiser is resolving. Never compare a sampled CVaR at small
alpha against an exact one. (3) NEW: `dCVaR/dp` is **identically zero iff `p(x*) >= alpha`**
(0 counterexamples in 3,000 cases) — a hard zero no shot count or baseline recovers. It is
convergence rather than failure, but at alpha=0.01 one random initialisation in four to
twenty **starts dead**.

**alpha is NOT a learning rate** (gradients at different alpha at the same point are not
parallel), and at small alpha the CVaR objective **collapses to an argmin-finder**, proved as
an identity.

**SIGN CORRECTION 2026-09-06 to the long-recorded "two AMBER variants became the same objective
at alpha <= 0.25" trap: THE SIGN IS BACKWARDS.** Measured, **small alpha DISCRIMINATES BETTER
between objectives; it is CONCENTRATION that collapses them.** The collapse is real but was
mis-attributed to the alpha level rather than to the distribution narrowing during
optimisation. This corrected framing is consistent with the sprint's central finding that
concentration is the destructive operation.

**Also corrected: the `baseline="tail"` defect's "cosine -0.023" figure.** The bias has a
CLOSED FORM, `-c * grad P(E < q)` verified to < 1e-9 — it adds the gradient of the tail
PROBABILITY, not noise — and -0.023 is one draw from a spread measured at **+0.06 to +0.96**.
Quote the closed form, never the cosine. In the sweep, plain expectation-value VQE (alpha = 1.0)
returns the BEST structure; lower alpha buys entropy and diversity monotonically and does not
buy accuracy.

**Encoding: binary, for four independent reasons.** It attains the information-theoretic
qubit bound at power-of-two k; needs **510 Pauli terms where one-hot needs 454,463**; is
surjective so needs no penalty; and has the highest gradient variance tested. Gray buys
nothing (a clean single-target positive was retracted after replication over 54 cells).
Methodological trap: **a non-surjective encoding's Pauli spectrum measures its CONSTRAINT,
not its objective**, which invalidates naive cross-encoding spectrum comparisons.

**Known defect to respect — CORRECTED 2026-09-05; the earlier figure was a factor of two
low.** FOUR torsions per chain are inert for the CA trace: `phi[0]`, `psi[0]`, `phi[n-1]`,
`psi[n-1]`. **Both terminal residues are entirely invisible to CA-RMSD**, so the dead block
is **`2*log2(k)` qubits per chain**, not `log2(k)`. At n=9, k=4 the live count is **14 of
18**, and the k=4 operating point is **21.9 mean live qubits, not 25.9**. Every
`n_res*log2(k)` count in this project is an overcount by `2*log2(k)`.

This is a prediction of the locality theorem rather than a coding slip: `d_ij` depends on the
residues STRICTLY between i and j, and residues 0 and n-1 are strictly between no pair. It
was found independently by two agents from different directions. A 1-local prior therefore
spends **22%** of its per-residue decisions on states CA-RMSD cannot see. Note `phi[n-1]`
still places C/CB/O, so Legacy and AMBER *do* depend on states the CA metric cannot see.
Report live qubits, never nominal.

**DATA INTEGRITY — the cached AMBER subset is ORACLE-CONDITIONED.** In
`s13/results/qarch_enum_<PDB>.npz` the 2,955-configuration AMBER subset (`amber_total`,
`amb_*`) is **40% oracle-conditioned**: 0.401 A better in RMSD than the space it was drawn
from, and NOT a uniform sample. **Only rows with `amber_kind == 0` are usable as an unbiased
sample.** Any correlation, AUC, training label, feature or normalisation constant computed
over the full 2,955 rows is conditioned on a selection that partly encodes the native. The
`rmsd`, `legacy`, `prior` and `leg_*` columns are full-enumeration and unaffected.

**METRIC CAUTION — in-decile rank correlation is NOT monotone in objective quality.** On a
signal-tunable family it peaks near +0.446 at moderate signal and FALLS to +0.340 where the
objective is genuinely much better (global rho +0.933). As an objective improves, its lowest
decile becomes a narrower and more homogeneous set of good structures, so less RMSD spread
remains inside it to rank. A single in-decile number is ambiguous between objectives that
differ fourfold in real quality. **Always report global rho, in-decile rho, argmin RMSD and
decile-mean RMSD together.**

**Legacy is a clash gate, not a ranking function.** 98.8% of its variance is the `steric`
term, which is exactly zero on 73.5% of the space and **constant across its entire
lowest-energy decile** — a binary gate carrying no ordering information where a search lives.
Its legitimate role is a 10%-removal geometric veto, steric term only; the total is unsafe
as a filter above 10% removal.

**AMBER is a genuine physical validator and an inert refiner.** Best clash detector measured
(its low decile is 0.09% clashed at 1.98 A against Legacy's 0.56% at 1.62 A), and the only
physical objective with a positive excess decoy AUC. It can discard 90% of a pool and still
hold the best structure. But refinement is **safe and inert**: +0.014 A [-0.001,+0.028],
13W/32L — from 1.295 A it returns 1.321 A, from 4.334 A it returns 4.337 A. Use it once, at
the end, emitting a flag and a structure, **never as a score**.

**The decoy-discrimination floor.** Discrimination is a function of the QUALITY GAP, not
structural separation. **No objective exceeds 0.511 pairwise accuracy when two structures
differ by less than 0.25 A.** `leg_steric` needs 0.58 A to reach 55%, Legacy 1.31 A, AMBER
1.58 A (and 3.87 A for 60%). The useful search range is ~2.5 A wide, so both energies resolve
about one bit of it. Cached decoy sets with matched random-anchor nulls: `s14/ener_decoy.decoys()`.
Shared helpers: `s14/ener_lib.py` (tie-averaged `argmin_rmsd`, `decile_rho`, `pair_accuracy`),
normalisation constants `s14/cache/ener_norm.json` (**use `grad_sd`; never raw AMBER's sd**).

**Retrieval-side closures** (do not reopen without a new mechanism): in-band ranking over
retrieval pools is **signal-limited, not sample-limited** (a set-transformer over the full
signed deviation tensor has a FLAT learning curve, 3.043 at n=8 to 3.026 at n=75, while
the same harness with a leaked label reaches 2.534 at n=8); the terminal operator consumes
the set MEAN not the set BEST (`d_out = 1.16*d_set_mean + 0.04*d_set_best`, R^2 0.893);
post-hoc correction of the distance objective is dead twice over; growing the fragment
library buys ~0.43 A per decade.

**No fresh benchmark exists.** All 204 identity clusters of 9-16mers in the corpus are
spent. The containment-fresh world supply is 16 targets, 10 of them amyloid fibrils.

**Failure class.** 10 of the 18 catastrophic targets are fibril segments or lasso
peptides (10/18 vs 6/108, Fisher p = 1.2e-6). A linear-window architecture has no
mechanism for a covalent thread or a lattice contact.

---

## 3. THE MACHINE. Read this before launching anything.

15.6 GB RAM total, **and the box is shared** — at sprint start only ~2.7 GB was free
(browser + editor + agents). 8 logical cores, but they are **4 fast + 4 slow = 6.43
core-equivalents**; measure with the real workload before blaming contention.

Rules:

  * Set thread caps at import: `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=2`.
    `s13/qarch_lib.py` already does this; import it early.
  * Before any allocation over ~500 MB call `s13.qarch_lib.wait_for_memory(min_gb=...)`
    or `core.quantum.mem_gate(...)`. **Wait for headroom — never catch a MemoryError and
    return NaN.** A NaN column silently becomes "this model has no signal", which is a
    fabricated negative.
  * Do NOT load `esm_cache.npz` (1.5 GB). Use `s12/esm_bank.py`, the 100 MB compaction.
  * Two heavy jobs fit. Three do not. If your job is heavy, say so in your first message
    to the coordinator so the schedule can be arranged.
  * Cache expensive quantities to `s14/results/` or `s14/cache/`. AMBER single points cost
    **28 ms distinct** (not 6 ms — the original timing hit a result memo and was wrong for
    an hour); `refine_coords` costs 9.1 s; a Legacy evaluation costs 0.3 ms.

---

## 4. THE MACHINERY THAT ALREADY EXISTS. Use it; do not rebuild it.

**`s12/instrument.py`** — the shared measurement instrument. Every agent scores against
the same objects so the benchmark is never touched.

    from s12 import instrument as I
    I.targets()                 # the 126 dev targets: pdb, seq, n, fold
    I.load_univ(pdb)            # cached window universe: W/PHI/PSI/S/org/sim/order/rr/nat_ca
    I.build_ca(phi, psi)        # ideal-geometry CA trace; accepts (n,) or (B,n)
    I.ca_rmsd(a, b)             # Kabsch CA-RMSD
    I.kabsch_rmsd_batch(A, b)   # batched
    I.project(C, seq, fold)     # projection onto the ideal-geometry manifold
    I.paired(a, b, ...)         # paired bootstrap CI + W/L + drop-top-10/20
    I.write(name, obj)          # atomic, with n_expected/complete flags
    I.FAIL18                    # the 18 zero-recall targets

`python -m s12.instrument` must reproduce, exactly:
`shipped 3.4540004952559396, pool_best 1.7108244199364904, top75_best 2.3061526409453816,
synthesis_fit 3.2040761603809194, n_zero_recall 18`. **Run it at the start and end of your
work.** If it drifts, stop and tell the coordinator.

**`s13/qarch_lib.py`** — the torsion state space and its energies.

    from s13.qarch_lib import Space, empirical_prior, prior_energy, legacy_energy, \
                              legacy_components, amber_energies, wait_for_memory, \
                              spearman, percentile_of, write
    sp = Space("1CS9", k=4)     # sequence-conditioned, TARGET HELD OUT
    sp.PHI, sp.PSI              # (n, k) radian tables
    sp.uniform(B, rng)          # (B, n) random configurations
    sp.sample_prior(P, B, rng)
    sp.ca(S) / sp.coords(S)     # decode to structure
    sp.rmsd(S)                  # ORACLE post-hoc scoring only
    sp.ORACLE_snap()            # ORACLE
    empirical_prior(sp)         # leakage-safe (n, k) state occupancy

**`s13/results/qarch_enum_<PDB>.npz`** — nine targets **fully enumerated** at k=4, 262,144
configurations each, already scored. Keys: `rmsd`, `legacy`, `prior`, eleven `leg_*`
component columns, and `amber_total` + five `amb_*` columns on a 2,955-config subset.
Targets: 1CS9 2MK7 2P5H 6EY3 6F3V 6S0N 7N2I 8IS3 9UV5. **This is 2.36 M labelled
structures with true RMSD and is the single most valuable asset in the repository.** Use
it for exact analysis, for training a structural objective, and for exact VQE controls.

**`core/quantum.py`** — the quantum machinery.

    StatevectorCircuit(n, layers, ring)   # exact statevector, n <~ 20
    .state / .probs / .probs_batch / .n_params
    cvar_exact(E, p, alpha)               # exact CVaR from a full distribution
    grad_cvar_paramshift / grad_cvar_fd / grad_cvar_score
    run_cvar_vqe(E, alpha, T, n, ...) ; run_global_cvar_vqe
    MPSAnsatz, OneLayerAnsatz, Adam, all_bitstrings, BestSeenTracker
    cvar(energies, alpha) ; cvar_from_probs ; alpha_schedule
    mem_gate(tag, limit) ; limit_threads(n)

**KNOWN DEFECT:** `qansatz.cvar_gradient` with `baseline="tail"` centres the score-function
baseline on the tail only and is biased (cosine -0.023 with the true gradient).
`baseline="const"` is correct. **Verify the estimator yourself before trusting it** —
section 2.B of the sprint brief requires deriving it, checking alpha dependence, and
checking finite-difference consistency. Do not assume prior tests were adequate.

Other modules: `torsion_lib2.py` (libraries, `library_for(seq, k, exclude_seq)`,
`PerResidueTorsion`), `core/energy.py` (`components_batch`, the Legacy terms),
`core/amber.py` (`single_point`, `refine_coords` — real ff14SB/GBn2 via OpenMM),
`core/geometry.py` (`build_backbone_batch`), `core/project.py`, `budget.py`
(`BudgetedEnergyModel` for evaluation accounting), `peptide_db.py` (`holdout(seq)`),
`representations.py`, `distogram.py`, `pairnet.py`, `sidechains.py`.

---

## 5. STATISTICS AND REPORTING. Non-negotiable.

Every comparison reports: **paired mean difference, bootstrap 95% CI, median, mean/sd, W/L,
per-fold values, and a drop-top-10 / drop-top-20 concentration check.** Use `I.paired`.
A result carried by two targets is not a result; concentration analysis is mandatory, always.

**CORRECTION 2026-09-06, and it corrects this brief.** A raw drop-top threshold is **NOT a
valid concentration test on its own.** When an effect's mean is small relative to per-target
spread, discarding the ten most favourable targets removes a large share of the total **even if
every target carries an identical effect** — so at low signal-to-noise the check *must* fail a
uniform effect. It fired correctly on `leg_torsion` (a genuine two-of-nine artefact) and
**misfired** on a real result in the same sprint; the difference is nothing but signal-to-noise.

  * **Compare every concentration statistic to a simulated UNIFORM-EFFECT null** of the same
    mean and sd. A bare `top10_share > 1.0` is not a failure — that null's own mean was 1.410
    at n=50.
  * **Print `mean/sd` beside the check** so a reader can see when it has no power (0.28 in the
    case that misfired).
  * **State the sample size inline**: at n=50 drop-top-10 discards 20% of the sample, at n=126
    only 8%.
  * When the check has no power, **lean on the per-fold table**, which is the genuinely
    informative one.
  * The check is **NECESSARY AND NOT SUFFICIENT**, and *"not concentrated"* is **not**
    *"demonstrated uniform"*.

**A related reading hazard, learned the same day.** `I.paired` returns `drop_top10_mean_diff`
and `top10_share` **in the same dict** as the W/L. A claim was published on the reassuring
field while the disqualifying ones sat unread beside it. **The failure mode is reading PAST a
check, not forgetting to run one** — so emit these as one PASS/FAIL verdict block, never as
separate fields a reader can select from. And note: a near-even W/L *with* a CI excluding zero
is **suggestive** of concentration and a prompt to run the null-calibrated check — the
median-versus-mean gap flags it for free — but it is not proof of one.

Every attractive correlation gets a **null control** — a permutation, a label shuffle, a
matched random arm, or a sequence-blind twin. Report the null next to the effect.

Do not call 0.02 A an improvement without statistical support. Do not select seeds after
seeing results. Predefine selection logic. Do not cherry-pick targets.

**Tier every claim explicitly** in your findings file:

  * **DEMONSTRATED** — measured on legitimate inference-time information, with CI and null.
  * **ORACLE DIAGNOSTIC** — reads the native; prices a ceiling; never a headline.
  * **HYPOTHESIS** — proposed, not yet measured.
  * **LITERATURE-SUPPORTED** — from a citation, with the citation.
  * **REFUTED / INVALIDATED** — with the evidence that killed it, kept in place.

This tiering exists to prevent the project's most dangerous failure mode: an oracle
diagnostic quietly becoming a headline.

---

## 6. REPRODUCIBILITY

Every final experiment records: random seed, target ID, sequence, representation, qubit
count, ansatz, depth, Hamiltonian terms and coefficients, optimiser, learning rate, CVaR
alpha, iterations, initialisation, compute budget, energy, RMSD, post-processing,
refinement, final structure path, provenance. Write JSON to `s14/results/` via
`I.write`. Deterministic given the seed.

---

## 7. HOW TO WORK

Write reusable modules under `s14/`, not fragile one-off scripts. Prefix your files with
your workstream tag. Write your findings continuously to `s14/<TAG>_FINDINGS.md` — do not
save them all for the end, because the coordinator reads them to schedule the next wave.

Use the Write tool for file creation. Heredocs break on apostrophes and on `\n` inside
Python string literals; this has cost this project hours more than once.

Set `PYTHONIOENCODING=utf-8` and keep console output ASCII — the Windows console is
cp1252 and a stray arrow character will kill a long run at the last print.

When you finish, your final message must state: what you measured, the numbers with CIs,
what you refuted (including your own hypotheses), what is still open, and the exact
commands to reproduce.
