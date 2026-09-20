# AUDIT Z — adversarial audit of `s30/REPORT_S30.md` §0, §2, §5, §6, §9, §10, §11, §12, §13

Lane Z, second report adversary, contract rule 28. 2026-09-20.
Read-only. No project file was modified; this file is the only output.
Report state audited: **1313 lines** (the file grew from 1019 to 1313 during the audit — §6.4, §7
and §8 landed mid-pass; all quotations below were re-checked against the 1313-line state at
14:17:47).
Verifier state audited: `s30/s30_verify.py`, **run**, 36 matched / 0 mismatched / 0 missing.
Lane V's `AUDIT_V.md` read first, as instructed; nothing V confirmed is re-litigated here.

---

## VERDICT (read this first)

**§0 is NOT safe to publish as the sprint's headline as it stands.** Four defects change what it
says, and three of them are on the two sentences §0 exists to deliver:

- **§0.3 and §13.3 call the sprint's deliverable "native-free". It is ORACLE.** `coh` is
  `corr(r, μ)` where `μ = pool75_mean − d_nat` and `r = (expected − d_nat) − correction`
  (`s30/s30_P_lr.py:58,78,202,204`). **Both arguments require the native.** §12.2 prints the
  `d_nat` in the definition on the same page that §0.3 and §13.3 call it native-free. The
  LEDGER never claims native-free; the report added it. (D1)
- **§0.1 row 6's "−0.0406 Å" is not the AMBER relax gain.** It is the high-minus-low dispersion
  *contrast* within chain-length tertiles (`s30_G_disp2.json → length_matched_split.note`:
  *"high vs low dispersion WITHIN chain-length tertiles"*). The deployable gain is
  **−0.0221** (`endpoint_arithmetic.whole_sample_gain`), which is also where the *same cell's*
  "0.69% of baseline" comes from (`pct_of_baseline = 0.6875`). The cell states two numbers that
  cannot both describe one quantity, and §1.1 already prints the right one. (D2)
- **§0.1 row 6 calls +0.0554 Å "the best *measured*, applied prior correction". Its own artefact
  says `NOT MEASURED (|effect| 0.0554 <= its own MDE 0.0820, 0.68x)`** — below the 0.7× "not a
  result" floor that §0.1 row 8 states two lines later. The ratio appears **nowhere** in the
  1313-line report, in any of its five occurrences. (D3)
- **§0.3 and §13.3 strip the ORACLE / NOT-DEPLOYABLE label off −0.126 / −0.247.** Those are the
  two most quotable graphics in the document (the gate box and the architecture diagram), and
  they are the only two of six occurrences where the label is missing. This is the number you
  flagged as the sprint's most misquotable, and it is unlabelled exactly where it will be
  screenshotted. (D4)

On the counts (D6): I checked all seven. Two were wrong when I started and **you fixed both
mid-audit** ("thirty-odd, eleven of them mine" → the table's 26/10; "ten further defects" → 18).
**"Twelve lanes" is still wrong in two places** — §2.1's own table now lists thirteen — and the
"four of ten falsifiers" claim has the right number and the **wrong composition**: lane P's P1
*held*, by the ledger's own words, and the fourth is lane F's F2a/F2c, which Appendix A.5 already
records (D7).

**§5 is clean.** All four code citations hold, the four-row `circ_*` table reproduces to four
decimals from the cache, and the 9-qubit circuit genuinely ran in S30 — I traced it to
`ORACLE_ITERS = 300`, `N_QUBITS, LAYERS = 9, 3`, an exact parameter-shift Jacobian and `adam`.
One line-range citation is wrong (D13).

**§10 is clean** — all thirteen meter numbers reproduce exactly — except a stale `22/22` (D16).

**§9's citations all exist in the ledger and nine of ten are characterised faithfully.** Two
rows upgrade or merge what the ledger says (D10, D11), and D11 is a *third* instance of lane V's
D7 pattern: the report dropping a caveat the lane attached to its own number on purpose.

The surviving science is not damaged by any of this. Every defect below removes a label, a
stratum, a ratio or a count. **None removes a result** — but D1, D2 and D4 remove the three
things a reader would quote from §0 without reading further, and D1 changes what the deliverable
*is*.

---

## CONFIRMED

Recomputed from artefacts, code or the ledger. Nothing here is taken from the report's prose.

### §0.1 — the endpoint
`3.2105 Å, unchanged` ✅. Verifier reproduces it from `s29/results/s29_O_chain_rows.jsonl`
(`item=prod`), n = 126, and nothing was deployed. Rows 1–11 map one-to-one onto the charter's
questions 1–11 (`s30/BRIEF.md:558-569`). Row 8's MDE rule is the project's standing one. Row 9's
FAIL18 circularity matches S30-L23.

### §0.2 — the −0.3259 headline: **exact, and correctly attributed**
`s30/results/s30_P_chain.json → arms.ORACLEsign_LFOmag`:
```
chain_mean 2.886729   delta -0.325896   x_mde -2.0535   folds_same_sign 5   W/L 92/34   BETTER
PROD_chain_mean 3.212625      (cloud_delta -0.356687)
```
**−0.3259, 3.2126 → 2.8867, 2.05×, 5/5 all reproduce to four decimals.** 2.8867 < 3.00, so
"clears the charter's primary target" (`BRIEF.md:73`) is right. The arm is labelled ORACLE
throughout. The 3.2126-vs-3.2105 footnote is correct and matches §1.1 and the verifier
(chain disagreement 0.0021, cloud disagreement 0.0000).

### §0.3 — the `coh` numbers: **faithfully transcribed**
`s30/results/s30_P_lr.json → incoherence`:
```
UNCORRECTED            0.69310      N1 0.78574   N2 0.78273   N3 0.91717
IIDmatched R2 0.235    0.53652      IIDmatched R2 0.156  0.58676
```
0.6931 ✅; `0.783 / 0.786 / 0.917` ✅ (N2/N1/N3, order swapped, values exact);
`0.587 / 0.537` ✅. The 5.1% orthogonality discount and the 3.01× requirement are S30-L18
(`LEDGER.md:2561-2562`) ✅. §12.2's definition of `coh` matches `s30_P_lr.py:197-210` exactly,
including "within-target, mean over targets".

### §0.2's supporting measurements
- Sequence R² **0.83%** against the registered **1.96%** bar ✅ (`s30_P_lr.json → bars.B2 = 0.0196`;
  S30-L25 *"Registered bar 1.96%. Best arm 0.83%… P1 holds"*).
- Pool **+0.0758** ✅ = `nested.LONGRANGE_sep_ge7.N3_plus_pool.increment_over_N1 = 0.07583`.
- Uncollapsed posterior **−0.018** ✅ = the same block's N2 `increment_over_previous = −0.01840`.
- Coherence **0.6931 → 0.9172** ✅.
- **68% in |i−j| ≥ 7** ✅ and the five-signs reduction ✅ (S30-L25).
- `+0.0554 / −0.2466` ✅ as raw values (`applied.N3_plus_pool.delta`, `applied.IIDmatched_R2_0.235.delta`).

### §2 — the lane and entry arithmetic
Counted directly from `## S30-L…` headings: **28 entries, L0…L27, no gaps, no duplicates.**
Per-lane: D 7, L 5, T 3, F 3, X 2, R 2, Q 2, **P 2**, G 1, coordinator 1 = 28. The table's
current state (P = 2, coordinator row added) is **correct**; the earlier `P | 3` and the
`S30-L16` mis-citation for lane T's retraction were both fixed during this audit.
"Ten pre-registration files" ✅ — all ten exist, and the named set matches.

### §5 — audited against the code, not against lane W
Every line citation holds:

| claim | checked | result |
|---|---|---|
| `core/pipeline.py:179` `quantum: bool = False` | `grep -n` | **:179 exact** |
| `core/pipeline.py:241` `PROD = Config()` | `grep -n` | **:241 exact** |
| `core/pipeline.py:173-178` *"VQE/CVaR do not participate at all"* | read | **verbatim** |
| `core/pipeline.py:821` the withdrawn +0.113 Å | read | **:821 exact**, sentence spans 821-822 |
| `core/pipeline.py:838` = `H = diag(zrank(top-128 scores))` | read | **:838 = `E = _zrank(np.asarray(pool["sc"], float)[o])`** ✅ |
| `s30/THEORY.md:450-453` removes ansatz/optimiser/pool/search | read | **:450-452**, and also removes the encoding and the objective's functional form |

**"A 9-qubit circuit *did* run in S30" — CONFIRMED, from the code path, not from prose.**
`s29/s29_D_cost_audit.py:264-266` constructs `Q.StatevectorCircuit(A.N_QUBITS, A.LAYERS)` and calls
`A.oracle_circuit_ceiling(circ, frame, cand.nat_ca, starts=1)`. In `s27/s28_A_amp.py`:
`N_QUBITS, LAYERS = 9, 3` (`:54`), `ORACLE_STARTS, ORACLE_ITERS, ORACLE_LR = 5, 300, 0.05` (`:60`),
and the inner `fg` builds `J = 0.5 * circ.states_batch(th[None,:] + math.pi*np.eye(P))` — an
**exact parameter-shift Jacobian** — fed to `adam(fg, th0, iters, lr)`. So "9-qubit depth-3,
exact parameter-shift, 300 Adam iterations" is exact. `s30/results/s30_D_ladder_structs/` holds
**126** `.npz`, mtime 2026-09-20 12:44 (i.e. written during S30, not copied), and `1A13.npz`
carries `"dev_circ_s0": 0.0` against the 1e-6 assertion. **"No VQE" and "no quantum compute" are
indeed different claims and only the first is true.**

**The four-row table reproduces exactly** from `s30_D_meter_DIS_{chain,ca}.json → mean_rmsd`:
```
              chain      ca
circ_best     0.2516    0.2884
circ_s0       0.3175    0.3854
PROD          3.2071    3.0483
circ_opt      3.4330    3.3850
```
All eight values match the report to four decimals, and the `PROD`-is-the-meter's-re-projection
footnote is correct.

### §6 — the discipline claims that are code claims
- `s24.stats_lib._verdict` **does** refuse a verdict when `folds is None` rather than falling back
  to the IID CI — its own docstring says so and the branch is there ✅.
- `s30_D_meter.py:451` **is** `g["rand_signed"] = "PASS" if e["effect"] > 0 else "FAIL"` ✅, and
  the artefact's `contrast.verdict` really does read `WORSE` for the good-direction `+0.1716`
  while `verdict.gates.rand_signed` reads `PASS` — the verifier prints both ✅.
- Multiplicity **30 chain / 33 cloud** ✅ (`multiplicity.comparisons_emitted`), and the quoted
  sentence is verbatim from the artefact ✅.
- §6.3's median-vs-mean item ✅: `mean +0.0357 median +0.0000 drop-top10 +0.0603 … 51th pct of
  the null … flag False`.
- Lane R's max-over-43 sign-flip null at `p_max = 0.0` ✅ and lane Q's across-target null ✅.

### §9 — every citation exists in the ledger
Blau & Michaeli, Kennedy & O'Hagan, Brynjarsdottir & O'Hagan (*Inverse Problems* 30:114007, 2014),
Maehara (*ORL* 43:526, 2015), Wilder (AAAI 2018), Ohsaka & Yoshida (2017), Nemhauser–Wolsey–Fisher,
Das & Kempe, Goldberg (1984), Maurey, ANDIS/Yu et al. 2019, DOPE/Shen & Sali, Klenin & Langowski,
Abe et al. arXiv:2302.00704, PMC2662860 — **all fifteen present**, none invented.
Characterisations checked against the ledger body:

- **Maehara / Wilder / Ohsaka & Yoshida** — the report's rows and §9.3's "Mismatch stated rather
  than hidden" are **verbatim** from `LEDGER.md:799-807`. Nothing upgraded.
- **Nemhauser–Wolsey–Fisher / Das & Kempe** — "Monotonicity, not submodularity, is what we fail
  first" is the ledger's own sentence (`:771-772`) ✅.
- **Kennedy & O'Hagan / Brynjarsdottir** — verbatim from `:712-716`, and "closes E1" matches
  S30-L7 ✅.
- **ANDIS** — the quoted sentence is verbatim from `:639-640` ✅.
- **Abe et al.** — "already in the S29 index… used to reject without re-running" ✅ (`:774-775`).
- **Klenin & Langowski** ✅; the audit that caught lane G's own non-rotation-invariant first
  implementation is in the ledger and does not contradict the report's row.
- **§9.2 (PMC2662860)** — **all five numbers exact**: 63.0% → 1.09%, PULCHRA 3.64%, 3.36 vs 3.28
  (+0.08 worse), 2090 proteins under 200 residues. The ~2/n mechanism and the "large-`n` limit of
  ours" reading are S30-L13's own ✅. This is the best-argued row in §9.
- **§9.3's rejections** — DPPs (V constant on centroid-equivalence classes) ✅; submodularity ✅;
  QUBO/Ising "moot… Frank-Wolfe strictly upper-bounds any circuit" ✅ verbatim; Ramachandran
  (36.1° vs 36.4°) ✅.

### §10 — the meter, all of it
Recomputed from `s30_D_meter_DIS_{chain,ca}.json`:
```
ladder rho S28      chain -0.4023 [-0.477,-0.322] 5/5     ca -0.1818 [-0.308,-0.053] 4/5
ladder rho CHARTER  chain -0.1964                         ca +0.2603
ladder rho FULL     chain -0.3187                         ca +0.1179
cosine              ca -0.0339   z -2.4468
native_pctile       0.3676 [0.306,0.405]  (identical on both bases)
pref circ_best      chain 0.0714          ca 0.2063
pref pool member    chain 0.0201          ca 0.1265
production RMSD     chain 3.2071          ca 3.0483
GATE                chain BLOCK           ca BLOCK
```
**Every number in §10.1 and §10.2 matches.** The charter column matches `BRIEF.md:253-256` exactly.
§10.4's cosine nulls reproduce from `cosine_null`: `per_draw_abs_mean 0.13983`,
`target_mean_null.sd 0.01444`, `[p2.5, p97.5] = [-0.01689, +0.03075]`, `z = -2.4468`; the quoted
"a cosine of +0.05 is far below 0.140 and still several sigma above chance" is **verbatim from the
artefact's own `note`** ✅. §10.5's draw noise **59.4% / 27.1%** ✅ (verifier). The gate's rule is
quoted verbatim from `verdict.rule` ✅.

### §11 — the open items are real
- §11.1's "the named gap, stated by the lane that produced it" is at `STATE.md:443` ✅ and the
  report does not strengthen it.
- §11.3's launch-gate row: `s26/jobrun.py:39` is `CPU_START, MAX_CONCURRENT = 85.0, 4` ✅;
  the refusal at `:137` is `crowded = (not stale) and (ram > CEILING or cpu > CPU_START …)` ✅;
  `s26/governor.py:61` is `CPU_CEILING = 101.0` with the comment "CPU-triggered suspension is
  DISABLED" ✅. **All four code claims hold.**
- §11.4's "leads register, not a task list" is **verbatim** from `BRIEF.md:31` ✅.
- §11.2's chirality scoping ("the theorem is length-free, the emptiness is not") matches S30-L26 ✅.

### §12.4 — the three-way bit comparison
`0.376 Å/bit` recomputes from lane T's law: `0.2191 × (3.05 − 1.331) = 0.3766` ✅.
`0.0713 = 0.3567/5` (cloud) and `0.0652 = 0.3259/5` (chain) ✅; ratios 5.3× and 5.8× ✅; lane T's
independent 3× ✅ (V confirmed 0.132 vs 0.044). The "unregistered / not exchangeable" caveat is
attached ✅.

### §13.1 — the production diagram, checked against the code
`k = 500` (`core/pipeline.py:125`) ✅; `m = 75` (`:126`) ✅; BLOSUM62 (`:10`, `:673`) ✅;
ESM-2 650M (`core/data.py:886 ESM_MODEL = "esm2_t33_650M_UR50D"`) ✅; **17 bins**
(`core.predict.NBINS = 17`) ✅; **L1** Bayes risk (`core/predict.py:420`, `np.abs`) ✅;
coordinate average → CA cloud 3.0483 ✅; multi-start projection, seed not pinned ✅;
built chain 3.2105 ✅; quantum stage off the path with both line citations exact ✅.
**The two values on the arrows are right and their bases are right.**

### §13.2 — what was built
`s30/s30_verify.py` **runs, 36 matched / 0 mismatched / 0 missing** — I ran it. The
ledger-path assertion is real and lists nine files, all present ✅. The BRIEF-77-minutes claim
checks out: S30-L0 stamped 12:35, `s30/BRIEF.md` first committed 92b145c6 at **13:52:21** = 77 min ✅.
The meter's four extensions (`verify`, both bases, R = 8 draws, multiplicity counter) are all
in `s30_D_meter.py` ✅.

---

## NOT RECOMPUTABLE

| claim | why |
|---|---|
| §2.1 "never more than eight concurrently, against the charter's floor of four" | No concurrency log is on disk that I could read cheaply. Prose. |
| §5.1 "three closures and a theorem" | A framing, not a count with an artefact behind it. |
| §5.3 item 14's "the same ladder on every target to **1.18% of range**" | Lane W sources it to `s25/results/q_gibbs.json`; the report carries neither the file nor the key. Recomputable in principle, not from `s30/results/`. |
| §9's opening "three papers closed project directions that no experiment here had the power to close" | Which three is never stated. Either name them or drop the count. |
| §0.2's "physics an operator on either, not a third source (theorem G1, **repaired by lane P**)" | The repair exists only in `STATE.md:433-436` (*"Lane P's repair to lane G's enumeration, which is tighter than G1"*). **It never entered the LEDGER.** The report's sentence is faithful to STATE; a reader sent to the ledger will not find it. |
| §6.2's "8 feature arms × 5 strata (lane P)" | The coordinator's own arithmetic; no artefact enumerates it. The paragraph correctly labels the total a lower bound. |
| §0.4's "one withdrawal was itself withdrawn" | True per A.3/A.7, but "thirty-odd" is not recomputable — see D6 for the number that is. |

---

## DEFECTS

Ranked by severity.

### D1 (SEVERE) — §0.3 and §13.3: the sprint's deliverable is called **native-free** and it is an **ORACLE** quantity

**Report says (§0.3):** *"Not a number — **a gate**, cheap and **native-free**, that no measurement
in this project's history would have passed."*
**And (§13.3):** *"**This gate is the sprint's deliverable.** It is cheap, it is **native-free**…"*

**The code says** (`s30/s30_P_lr.py`):
```python
:58   d_nat = np.linalg.norm(nat[ii] - nat[jj], axis=1)            # ORACLE
:78   Y.append(exp - d_nat)                                        # ORACLE
:202  mu = X[:, FEATS.index("pool75_mean")] - X[:, FEATS.index("expected")] + y
:204  def coh(r):  ...  np.corrcoef(r[m], mu[m])
```
`y = expected − d_nat`, so `mu = pool75_mean − d_nat` and `r = y − yhat`. **Both arguments of the
correlation are functions of the native.** `coh` cannot be evaluated on a target whose structure
is unknown.

**The LEDGER never says native-free.** S30-L27 §4 calls it *"the specification for a fourth
source, as a testable number and not an adjective"* and writes the definition with `d_nat_{t,p}`
in it. **§12.2 of the report prints that same `d_nat`** — so the document contradicts itself two
pages apart, and the false half is the one in the headline and the architecture diagram.

This matters operationally, not just verbally. A native-free gate could screen a corrector *at
inference*, on any target. An ORACLE gate can only be run on a labelled benchmark — which is
still a real and cheap deliverable (it costs one correlation on 126 dev targets instead of a
pipeline run), but it is a **development-time** admission test, not a deployment-time one.

**Corrected sentences:**
> **§0.3:** Not a number — **a gate**: one correlation on the 126 labelled dev targets, evaluated
> before a corrector is ever run through the pipeline, that no measurement in this project's
> history would have passed. **It is an ORACLE statistic** — `coh` is computed against the native
> pair distances — so it screens candidates on the benchmark, not at inference. That is what makes
> it cheap: it replaces an endpoint run with a correlation, not a native with a surrogate.
>
> **§13.3:** **This gate is the sprint's deliverable.** It costs one correlation instead of a
> pipeline run, it is computed on the benchmark's own natives (ORACLE, and therefore a
> development-time admission test rather than an inference-time filter), and it discriminates a
> −0.25 Å corrector from a +0.06 Å one at identical out-of-fold accuracy.

---

### D2 (SEVERE) — §0.1 row 6 and §11.2: **−0.0406 Å is a between-halves contrast, not the AMBER relax gain**, and the same cell's "0.69% of baseline" is computed from the number it displaces

**Report says (§0.1 row 6):** *"The one confirmed **deployable** effect (AMBER relax at k = 30) is
**−0.0406 Å**, 3.56× MDE, 5/5 folds on a length-matched split — **0.69% of baseline**…"*
**And (§11.2):** *"The AMBER relax at k = 30 is confirmed (−0.0406, 3.56× MDE, length-matched,
5/5 folds)…"*

**`s30/results/s30_G_disp2.json` says:**
```
length_matched_split:  effect -0.040647   mean_hi -0.042720   mean_lo -0.002073
                       ratio -3.5597   5/5   note: "high vs low dispersion WITHIN chain-length tertiles"
endpoint_arithmetic:   whole_sample_gain          -0.022074
                       gain_if_applied_only_to_hi_half -0.021457
                       pct_of_baseline             0.68755     (= 0.022074 / 3.2105)
```
`effect = mean_hi − mean_lo`. **−0.0406 is the difference between two halves of the benchmark**;
it is the size of the *grading*, not of any operation you can apply. The deployable quantity is
**−0.0221** applied everywhere, or **−0.0215** applied only to the divergent half — and §A.7 states
both correctly. **0.69% is 0.0221/3.2105**, so row 6 pairs a number with a percentage that belongs
to a different number.

§1.1 gets this right in the same document: *"the sprint's one confirmed deployable effect — the
AMBER relax at k=30 — is **0.0221 Å**."* The verifier prints it too: *"E2 at 0.0221 is only ~2x the
full 0.0107 spread."*

**Second problem in the same cell:** row 6 calls it *"the one confirmed **deployable** effect"*
while §3's E2 row says **"NOT CLOSED"** and §11.2 says *"**it is not deployable**, independent of
its size."* Three statements, two of them in §0's own neighbourhood.

**Corrected sentence (§0.1 row 6):**
> **Zero — nothing was deployed.** The best *applied* prior correction emits **+0.0554 Å worse**,
> which is **0.68× its own MDE and therefore NOT MEASURED**. The one effect that is both confirmed
> and applicable — the AMBER relax at k = 30 — is **−0.0221 Å over all 126 (0.69% of baseline)**,
> graded by pool divergence (the high-minus-low-dispersion contrast is −0.0406 at 3.56× MDE, 5/5
> folds, length-matched), with **97.2% of the gain in the divergent half**. **It is not deployable:
> its restraint constant was picked on dev-set RMSD and has no native-free selection rule.**

§11.2 needs the same split: the *relax gain* is −0.0221; the *−0.0406* is what confirms the
grading.

---

### D3 (SEVERE) — **+0.0554 Å is NOT MEASURED (0.68× MDE) and is called "the best *measured*" correction; its ratio appears nowhere in 1313 lines**

**Report says (§0.1 row 6):** *"The best *measured, applied* prior correction emits **+0.0554 Å
WORSE**."*
**§0.2:** *"A corrector fitted on the distogram and pool emits **+0.0554 Å worse**…"* → and from it,
*"the pool… is **strictly harmful when applied**."*

**`s30_P_lr.json → applied.N3_plus_pool`:**
```
delta +0.055421   mde 0.081957   x_mde 0.67621   W/L 55/71
ci95_fold [-0.005995, +0.118924]          <- INCLUDES ZERO
verdict "NOT MEASURED (|effect| 0.0554 <= its own MDE 0.0820, 0.68x)"
```
**0.68× is below the 0.7× floor §0.1 row 8 states as "not a result", and the fold CI includes
zero.** The other two fitted arms are 0.91× and 0.95× — also in the NOT MEASURED band.

`grep` over the whole report: `0.0554` appears five times (§0.1, §0.2, §3, §7.1 #55, §8 L1) and
**not once with its MDE ratio or its verdict.** This is the same class lane V flagged as D5 and
D18 — a rule stated in §6 and broken in §0 — but here it is load-bearing for a headline
("strictly harmful when applied") rather than for a supporting clause.

The ledger is the origin: S30-L27's §3 table prints the fitted arms' deltas with `WORSE` and no
ratio, while giving `1.23x` and `1.77x` for the i.i.d. arms. The report inherited the asymmetry
and then added the word "measured".

**Corrected sentences:**
> **§0.1 row 6:** …The best *applied* prior correction emits **+0.0554 Å worse — 0.68× its own MDE,
> fold CI [−0.006, +0.119], NOT MEASURED.**
>
> **§0.2:** A corrector fitted on the distogram and pool emits **+0.0554 Å worse (0.68× MDE, NOT
> MEASURED — the claim is that it does not help, not that it measurably hurts)**; a synthetic
> i.i.d. one at matched out-of-fold R² emits **−0.2466 Å better (1.77× MDE, 97W/29L)**. *A 0.30 Å
> swing at matched accuracy, and the three fitted arms are 0.68×, 0.91× and 0.95× of their own
> MDEs while both i.i.d. arms clear theirs.* **The swing is the result; neither fitted arm is a
> measured harm on its own.**

Note this makes the finding *sharper*, not weaker: the point was never that +0.0554 is a big
positive number, it is that the fitted arms sit at zero-to-slightly-worse while the matched-R²
i.i.d. arms clear their MDE in the other direction.

---

### D4 (SEVERE) — §0.3 and §13.3 drop the **ORACLE / NOT DEPLOYABLE** label off −0.126 / −0.247, in the two most quotable graphics in the document

The i.i.d. prices appear six times. **Four carry the label; the two that do not are the gate box
in §0.3 and the gate diagram in §13.3:**

| line | context | ORACLE label? |
|---|---|---|
| 59-62 (§0.2) | *"The i.i.d. arm is **ORACLE-constructed and not deployable** — it is the price of a channel nobody has, never an achievement"* | ✅ |
| **76 (§0.3)** | `ADMIT iff coh < 0.6931.     R2 0.16 -> -0.126 A;  R2 0.24 -> -0.247 A` | ❌ **none** |
| 246 (§3) | "a synthetic i.i.d. one" | ✅ |
| 640 (§7.1) | "an ORACLE-constructed i.i.d. corrector" | ✅ |
| 810 (§7.4) | "the ORACLE arms that lower it" | ✅ |
| 1172-3 (§12.3) | `ORACLE-CONSTRUCTED, NOT DEPLOYABLE / a PRICE for a channel nobody has` | ✅ |
| **1284 (§13.3)** | `coh < 0.6931 -> admit and measure the endpoint.  R2 0.16 -> -0.126 A; R2 0.24 -> -0.247 A` | ❌ **none** |

In both unlabelled cases the numbers are framed as **what you get if you pass the gate** — the
strongest possible reading, in a code block and a diagram, i.e. exactly the two objects that get
lifted out of a report intact.

**Also missing everywhere in the report:** the two arms' MDE ratios (**1.23×** and **1.77×** —
both below the 2× line §6.2 names as the sprint's survival standard) and the artefact's own
`TYPE-M ZONE: magnitude inflated ~1.04x` flag on the −0.126 arm. `grep TYPE-M` over the report
returns one hit, and it is about lane F's retrieval row.

**Corrected blocks:**
```
    ADMIT iff coh < 0.6931.
    The only arms that clear it are ORACLE-CONSTRUCTED and NOT DEPLOYABLE:
    R2 0.16 -> -0.126 A (1.23x MDE, type-M inflated ~1.04x);  R2 0.24 -> -0.247 A (1.77x MDE)
    -- the PRICE of a channel nobody has, never an achievement.
```
and the same two lines inside §13.3's `coh < 0.6931` branch.

---

### D5 (MAJOR) — §0 and §12 quote a **CA point cloud** number under a **built-chain** heading, with no basis stated

`+0.0554`, `−0.2466`, `−0.126` and `−0.247` are all `mean_cloud` against `prod = 3.0483`
(`s30_P_lr.json → applied`). S30-L27 §3 states this explicitly and flags it under contract rule 1:
*"Point cloud, production 3.0483, an **INTERMEDIATE** (built chain 3.2041, contract rule 1)."*

**The report carries none of it.** §0.1's table is headed by *"Final mean built-chain RMSD 3.2105"*
two rows above the cell that quotes `+0.0554`; §0.2's "0.30 Å swing" and §0.3/§12.3's price ladder
give no basis at all. This is lane V's D10 pattern recurring in the section V never saw — and here
the same section *does* footnote the basis for the −0.3259 arm (correctly), which makes the silence
on the other four numbers look deliberate rather than absent.

The conversion is known and measured: lane P's own annotation gives the cloud→chain transfer for a
prior correction as **0.92** on two independent arms (*"the projection ABSORBS about 8% of the
gain"*), and explicitly warns that the 1.16 set-mean coefficient must not be used for it.

**Corrected sentence (§0.2, and mirror it in §0.3 and §12.3):**
> …emits **+0.0554 Å worse**; a synthetic i.i.d. one at matched out-of-fold R² emits **−0.2466 Å
> better**. *(Both on the **CA point cloud** against production's 3.0483 — an intermediate, per
> contract rule 1. The measured cloud→built-chain transfer for a prior correction is **0.92**, so
> −0.2466 cloud is ≈ −0.227 on the endpoint.)*

---

### D6 (MAJOR) — the counts in §0.4 and §2.1

I counted all seven. **Two were wrong when I started and you fixed both mid-audit**
("thirty-odd, eleven of them mine" → the table's **26, ten of them mine**; "ten further defects" →
**eighteen**). Recorded here so the sweep is auditable, with the one that is **still wrong**:

| claim | counted | source |
|---|---|---|
| 28 ledger entries | **28** ✅ | `grep -c '^## S30-L'`; L0…L27, no gaps, no duplicates |
| per-lane entries in §2.1 | **D7 L5 T3 F3 X2 R2 Q2 P2 G1 + coord 1 = 28** ✅ | the headings' own lane stamps (`P | 3` was wrong at 14:11 and is now `2`) |
| ten pre-registrations | **10** ✅ | the ten `PREREG_S30_*.md`, from **eight** lanes — see D18 |
| "Twenty-six claims withdrawn, ten of them mine" | **26 / 10** ✅ | A.1 10 + A.2 4 + A.3 4 + A.4 3 + A.5 5 |
| "eighteen further defects… two of which inverted a headline" | **18 / 2** ✅ | `AUDIT_V.md` ranks D1…D18; V's verdict names two inversions |
| **"Twelve lanes"** (§0.4 line 85 **and** §2.1 line 162) | **13** ❌ | §2.1's own table now lists D, L, T, F, X, R, Q, P, G, V, W, Y **and Z**, plus a separate coordinator row |
| "four of the ten registered falsifiers fired against the lane that wrote them" | right number, **wrong composition** ❌ | see D7 |

**Corrected sentence (both places):** *"**Thirteen lanes**, 28 ledger entries, ten
pre-registrations."* §2.1's prose should also say *"thirteen lanes ran"* now that lane Z is in its
own table.

---

### D7 (MAJOR) — §2.2 / §0.4: "four of ten falsifiers fired against their own lane" has the right number and the **wrong fourth lane**

**Report says:** three bullets (D, T, G), then *"Lane P registered its bars at 1.96% / 12.82% /
39.44% R² before the regression existed and came in at **0.83%**"*, then *"Four of ten fired
against the lane that wrote them."*

**Lane P's falsifier did not fire against lane P.** `PREREG_S30_P.md:138-139`:
> **P1 (headline):** out-of-fold excess R² over C-DIM is **< 1.96 %** for the full feature set.
> Falsified if the fold-clustered CI excludes 1.96 % from above.

Measured 0.83%, and S30-L25 states the outcome in its own words: ***"Registered bar 1.96%. Best arm
0.83%, excess 0.95%, 0.75x MDE. P1 holds."*** **Lane P's registered prediction was confirmed.** The
prediction that *was* falsified in that lane's work is **the coordinator's** 3:1 on long-range
R² ≈ 0 (+0.1959) — §7.1 row 53 labels it correctly as the coordinator's.

**The genuine fourth is lane F**, and Appendix A.5 already records it: *"Lane F's own F2a and F2c —
both refuted by its own measurements"* (S30-L16: *"F2a REFUTED"*). A fifth candidate is lane Q,
whose fourth registered falsifier was refuted (§8, L6). Lane R's F-R1 prior **held** (registered
"does not fire at 4 to 1"; it did not fire). Lane X's counting falsifier **tied** (24.80% vs 25%)
and the report elsewhere says it decided nothing.

**Corrected bullets (§2.2):**
> - **Lane D** — its own falsifier for the combination question was not met (S30-L21).
> - **Lane T** — retracted its own headline on a null lane Q told it to run (S30-L15 §3c).
> - **Lane F** — its own F2a and F2c were refuted by its own measurements (S30-L16).
> - **Lane G** — **both** registered priors were directionally wrong, and it said so first (S30-L26).
>
> Two registered nulls **held**, which is what makes the four meaningful: lane P's P1 (bar 1.96%,
> measured 0.83%) and lane R's F-R1 prior. The prediction falsified in lane P's work was **mine**,
> registered at 3:1, that the long-range R² would also be ≈ 0; it is +0.1959 raw.

---

### D8 (MAJOR) — §12.3: "accuracy ≈ 0.60 against a *free* baseline of 0.627" reads an **all-pairs** accuracy against a **long-range-only** baseline, and both matched comparisons point the other way

**Report says:** *"The deployable sign currently performs at accuracy ≈ 0.60 against a *free*
baseline of 0.627; ≈ 0.8 is needed."*

**`s30_P_sign.json`:** `LFO_LOGIT.sign_acc_all = 0.6000`; `0.627` is
`per_bin["4"].baseline_always_plus` — the always-positive baseline **in the |i−j| ≥ 7 bin only**
(the five-bin baseline row is `0.556 0.571 0.579 0.643 0.627`). Matched:

```
all pairs :  LFO_LOGIT 0.600   vs  always-plus 0.5562   ->  ABOVE by +0.044
bin 4     :  LFO_LOGIT 0.643   vs  always-plus 0.6270   ->  ABOVE by +0.016  (0.10x MDE)
```
On **both** matched framings the deployable channel is *above* the free baseline, not below it.
The sentence as written implies a 0.027 shortfall that does not exist in either space.

This is the "statistic read against another statistic's null" class — the fifth instance in a
sprint whose §6.3 names it as the methodological finding.

The lane's own conclusion is different and better: *"Not one channel clears 1× MDE over the free
baseline in any of the five bins… The inferred 0.6 is confirmed and refined: raw accuracy
**0.60–0.67**, and essentially all of it is free."*

**Corrected sentence:**
> **The prize is five bits per target.** The best native-free sign channel reaches raw accuracy
> **0.60 over all pairs and 0.64–0.67 in the long-range bin**, against an always-positive baseline
> of **0.556** and **0.627** in the same two spaces — i.e. **essentially all of the realised
> accuracy is free**, and no channel clears 1× MDE over the baseline in any of the five bins. The
> matched-accuracy i.i.d. ladder crosses zero at ~0.63–0.65 realised accuracy, exactly where the
> free baseline already sits. **≈ 0.8 is needed** (ORACLE sign, accuracy 1.0, is −0.3259 Å on the
> chain; the accuracy ladder gives −0.2744 Å at 0.8 on the cloud).

---

### D9 (MODERATE) — "Five ORACLE signs **on long-range pairs**" — they are one sign **per separation bin**, and the ~4× is a different arm's number

**Report says (§0.2 and §12.3, identically):** *"Five ORACLE signs on long-range pairs are worth
−0.3259 Å on the built chain… The prize is **five bits per target**, concentrated **~4×** on the
tail and **68%** in `|i−j| ≥ 7`."*

`s30/s30_P_prior.py:46`:
```python
SEP_EDGES = [(2, 2), (3, 3), (4, 4), (5, 6), (7, 99)]   # the five-number profile, fixed bins
```
**The five signs are one per separation bin, spanning |i−j| = 2 to 99. Only the fifth is
long-range.** The sentence's own second half says 68% of the value is in the ≥ 7 bin, which is only
a meaningful statement if the other four bins exist — so §0.2 contradicts itself inside one
paragraph.

Second: the **~4×** is `ORACLE_SEPPROF5`'s (S30-L25 §4: other108 −0.359, tailA −1.405, tailB
−1.164 → 3.91× and 3.24×). For the **five-sign arm** the chain annotation gives other108 −0.2470,
tailA −0.8018, tailB −0.4492 → **3.25× and 1.82×**.

**Corrected sentence:**
> **Five ORACLE signs — one per separation bin (|i−j| = 2, 3, 4, 5–6, ≥ 7) — are worth −0.3259 Å on
> the built chain…** The prize is **five bits per target**, with **68% of the underlying profile's
> value in `|i−j| ≥ 7`** and a tail concentration of **3.3× and 1.8×** on the two filter-independent
> tails (the full five-number profile concentrates 3.9× and 3.2×).

---

### D10 (MODERATE) — §9.1's DOPE row moves "the top 13% has zero reference density" from the **n = 9** row onto the **13-mer**

**Report says:** *"The reference state's entire support for a **13-mer** is [0, 15.05 Å] against a
15 Å table cutoff — **the top 13% has zero reference density**."*

**`LEDGER.md:601-603` (S30-L6's table):**
```
 n      Rg      a    2a = the reference state's ENTIRE support
 9    5.07   6.55    13.09 A   <- DOPE is tabulated to 15 A: the top 13% has ZERO reference density
13    5.83   7.53    15.05 A   <- table and reference run out together
```
The 13% is `(15 − 13.09)/15` on the **n = 9** row. The **13-mer's own annotation is the opposite**:
support and cutoff *"run out together"*, i.e. there is no dead range at n = 13. The report merged
two rows and attached the more dramatic consequence to the wrong one.

**Corrected sentence:**
> The reference state's entire support runs out where the table does: **13.09 Å at n = 9 against a
> 15 Å tabulation — the top 13% has zero reference density — and exactly 15.05 Å at n = 13**, where
> table and reference end together. At n = 150 the support is 38.13 Å and the 15 Å cutoff sits
> comfortably inside. The authors say plainly that DOPE *"is less accurate for smaller proteins"*;
> lane L turned that into the mechanism for the field's 40–50 residue wall.

---

### D11 (MODERATE) — §9.1's Maurey row drops **both** caveats lane L attached to its own number

**Report says:** *"Maurey bounds any hull point within `R/√m`, which with the project's own
dispersion `√63.82 = 7.99` gives **0.22 Å at m = 75** — i.e. the grid is not the constraint."*

**`LEDGER.md:831-838`, lane L's own "WHERE I COULD BE WRONG":**
> Maurey's bound is for **multisets** (weights on a `1/m` grid), not the uniform-weight **subsets**
> the readout uses — the subset family is strictly poorer and the `1.4315` vs `2.3055` gap is the
> size of that difference, so *"the relaxation is tight"* is a claim about the **weighted** problem
> and a reader could wrongly take it as one about subsets. `R` should be a **max** over members and
> I used the RMS dispersion, so **the table is optimistic by perhaps 2-3×**.

The report takes the number, states the conclusion the lane warned against, and carries neither
caveat. **This is lane V's D7 pattern — the report being less careful than its own ledger —
recurring in the section V did not audit.** It is the third instance the two audits have found
together.

**Corrected row:**
> **Goldberg (1984); Maurey's empirical method** | Fixed-`m` is the hard formulation
> (densest-`k`-subgraph); free-`m` is poly-time by max-flow. Maurey bounds any hull point within
> `R/√m`: with the project's RMS dispersion `√63.82 = 7.99` that is **0.22 Å at m = 75**, so the
> grid is not the constraint. **Lane L's own two caveats travel with it: the bound is for weighted
> multisets, not the uniform subsets the readout uses (the 1.4315 vs 2.3055 gap is the size of that
> difference), and `R` should be a max rather than an RMS, so the table is optimistic by perhaps
> 2–3×. The ordering across `m` is unaffected.**

---

### D12 (MODERATE) — §6.3 calls the `compare` sign-convention trap **"latent"**; it **fired, in lane P**, and the sprint's own state calls that the fifth instance

**Report says:** *"And a fourth, **latent**, caught by the verifier rather than by a reader:
`s24.stats_lib.compare` is lower-is-better… **no S30 document quotes the inverted field**."*

**`LEDGER.md:3354-3358` (S30-L25, lane P's own "Two things I got wrong or nearly wrong"):**
> The verdict strings in my first run were **inverted**: `ST.compare`'s convention is
> lower-is-better and explained variance is higher-is-better, so an arm that was 0.0337 WORSE than
> its control printed `BETTER`. I caught it on the first table and negated both arms. It is the
> same failure lane D found in S29-L23's printed verdict, one sprint later, **in a lane that had
> just read the entry about it.**

**`STATE.md:450`** calls it *"**Fifth** instance of the statistic/convention mismatch this
sprint."* So the trap is not latent — it produced a wrong table in a lane that had read the
warning, and the lane caught it itself. "No S30 *document* quotes the inverted field" remains
true and is the weaker half of the claim.

Losing "latent" makes §6.3's own generalisation stronger, not weaker: a trap that fires in the
lane that just read about it is better evidence that the failure is silent than a trap nobody hit.

**Corrected sentence:**
> And a fourth, which is the same defect in the library rather than in a null: **`s24.stats_lib.compare`
> is lower-is-better** (`d = a − b`) because its native statistic is RMSD. Fed a preference rate,
> its `verdict` string reads `WORSE` for the good direction. **It fired this sprint** — lane P's
> first run printed `BETTER` for an arm 0.0337 worse than its control, *in a lane that had just read
> the entry about it* (S30-L25), and lane P caught it on its own first table. Lane D's gate handles
> it correctly (`s30_D_meter.py:451`, `PASS if effect > 0`), no S30 document quotes the inverted
> field, and `s30/s30_verify.py` now asserts the trap's shape so it stays known. **Counting lane P's,
> that is five instances of one defect in one sprint.**

---

### D13 (MODERATE) — §5.4 and §11.3 cite `s25/LEDGER.md:233-244` for a figure that is at `:280`

**Report says:** *"**S25-L5 withdrew that number** (`s25/LEDGER.md:233-244`: *"never a measured
effect … by this project's own fixed rule that is a NULL"*) and **replaced it with −0.1405 Å at
0.68× MDE**."*

`s25/LEDGER.md:232-244` holds the withdrawal and the figure **−0.1126 Å at 0.51× MDE**. The
**−0.1405 / 0.68×** replacement is at **`s25/LEDGER.md:280`** (`VQE_LFO minus argmin −0.1405 SE
0.0732 MDE 0.2051 **0.68x MDE** 66W/48L 5/5`). **Lane W cited both correctly** (`QUANTUM_W.md:55`
gives `:280` for the replacement); the report merged two citations into one range.

Minor second point: the quoted phrase *"never a measured effect"* is from the **heading** at
`:232`, outside the stated range.

**Corrected sentence:**
> **S25-L5 withdrew that number** (`s25/LEDGER.md:232-244`: the +0.113 Å *"was never a measured
> effect"*, measured −0.1126 Å at **0.51× MDE** with a fold CI spanning zero — *"by this project's
> own fixed rule that is a NULL"*). The replacement headline is **−0.1405 Å at 0.68× MDE**
> (`s25/LEDGER.md:280`), **which is itself below the 0.7× floor and is not a result** — the lane
> says so in its own handover caveat (`s25/LEDGER.md:1196`).

---

### D14 (MODERATE) — §12.5's "Closed by theorem or by price" contradicts §0.1 row 1 and §7's own three-way taxonomy

**Report says (§12.5):** *"**Closed by theorem or by price**, with the closure named: …"* — nine
items.

**§0.1 row 1 says** *"closed — by measurement, **by theorem**, or by price"*, and **§7 makes the
distinction its most reusable output**, with a counted split (*"58 killed by measurement, 19 by
theorem, 8 by price"*) and an explicit rule that the three are not interchangeable: *"a hypothesis
killed by measurement **can be reopened by a better instrument**; one killed by theorem cannot."*

At least four of §12.5's nine were killed by **measurement**, not theorem or price:

| item | actual closure |
|---|---|
| the field combination (rank, not count) | **measurement** — ORACLE global ρ = 0.1693 measured on 21 fields (S30-L21) |
| the second-moment/quadric escape | **measurement** — +0.2059 (1.81×) and +0.086 (2.95×), two samplers (S30-L12/L15) |
| generative spaces | **measurement** + algebra — four arms, `set_mean² ≈ B² + S²` (S30-L20) |
| recognition from single-structure geometry | **measurement** — 43 channels on a matched ladder (S30-L19) |
| sparse weighted readouts | **price** ✅ |
| subset objectives (T1) | **theorem** ✅ |
| torsion encodings | **price/arithmetic** ✅ |
| common-mode correction | **theorem** ✅ |
| achiral single-structure channels (G1) | **theorem** ✅ |

Under §7's own rule, telling the next sprint that all nine are closed "by theorem or by price"
tells it that none can be reopened by a better instrument. Four of them can.

§12.5 also **omits filter width**, which §0.4 lists among the directions closed.

**Corrected opening:**
> Closed, with the closure named **and its kind**, because the kind is what says whether a better
> instrument could reopen it. **By theorem** (cannot be reopened within the model class): subset
> objectives through an averaging readout (T1); common-mode correction from pool data
> (non-identifiable at any K); achiral single-structure channels (G1). **By price** (reopenable only
> by changing the budget): sparse weighted readouts; torsion encodings. **By measurement**
> (reopenable by a better instrument, and say so): the field combination (rank, not count); the
> second-moment/quadric escape (twice, independently); generative spaces (closed *jointly* with the
> readout); recognition from single-structure geometry; filter width as a free lunch.

---

### D15 (MODERATE) — three unbounded generalisations from three arms of one regression

Three sentences in §0.3 and §13.3 quantify over the project or over all accuracies, from
`s30_P_lr.json`'s **three nested fitted correctors** (N1, N2, N3) and **two ORACLE arms**:

1. §0.3: *"**every corrector this project owns**"* — three nested blocks of one lane's
   regression, not an inventory.
2. §0.3 / §13.3: *"**no measurement in this project's history** would have passed it"* — `coh` was
   computed for five arms, all in S30-L27. No historical corrector was retrospectively scored.
3. §13.3: *"coh ≥ 0.6931 → REJECT. **It will emit WORSE at any accuracy.**"* — the evidence spans
   R² 0.152 to 0.235. "At any accuracy" is an extrapolation in both directions.

The ledger's own wording is narrower — *"Every channel this project owns raises it"* — and even
that is the lane's generalisation. The claim is interesting enough at its measured size.

**Corrected block (§0.3):**
```
coh = corr(a corrector's residual, the pool's common-mode pair error)      [ORACLE; see 12.2]

    uncorrected                                    0.6931
    all three fitted correctors measured (S30-L27)  0.783 / 0.786 / 0.917   <- ALL RAISE IT
    imposed-structure ORACLE arms only              0.587 / 0.537
```
and in §13.3: *"coh ≥ 0.6931 → REJECT. Every corrector measured so far that lands here emits
worse, at R² from 0.15 to 0.24."*

---

### D16 (MINOR) — §10.1's "22/22" is stale; the verifier runs **36**

§10.1: *"Verified by `s30/s30_verify.py`, **22/22 matched, 0 mismatched, 0 missing**."*
§7.1 row 17 repeats *"22/22 matched, 0 mismatched, 0 missing"*.
§13.2: *"36 headline numbers recomputed from artefacts, 0 mismatches."*

I ran it: **`MATCHED: 36    MISMATCHED: 0    KEYS/FILES NOT FOUND: 0`**. Update both, and keep lane
V's D16 caveat attached — the verifier's coverage is still concentrated on lane D's artefacts plus
the endpoint reproducibility and the ledger-path assertion, so *"verified 36/36"* should not be read
as coverage of §4's or §12's numbers.

---

### D17 (MINOR) — **"theorem G1" names two different theorems**, and §0 uses the label for the one §7 calls G1b

- §2.1, §3, §11.2, §12.5 use **G1** = the chirality dichotomy (*"every reflection-invariant
  single-structure channel is a distance-map reading"*, S30-L26).
- §0.2 and §12.1 use **G1** = the source enumeration (*"physics an operator, not a third source"*),
  which §7 correctly calls **Corollary G1b** (*"this project has exactly three references…"*).

Two different statements under one label, one of them in the most-read paragraph in the document.
Use **G1b** in §0.2 and §12.1. See also the NOT RECOMPUTABLE note: lane P's repair exists only in
`STATE.md:433-436` and never entered the ledger.

---

### D18 (MINOR) — §6.1's "every lane did" pre-register

The ten prereg files come from **eight** lanes (D, F, G, P, Q, R, T, X). Lane **L** — five ledger
entries, three theorems — registered nothing, and neither did V, W, Y, Z or the coordinator.
Corrected: *"every **measurement** lane did; four lanes' falsifiers then fired against them and
were reported as failures."* (Note "three" → "four" to match §2.2 after D7.)

---

### D19 (MINOR) — §10.1's "all four anchors recompute" needs §10.3's caveat inline

The native-percentile row reads *"charter ~36.9 | measured 0.3676"*. **0.3688 (= 36.88) is
`DIS_SURR`'s**; the shipped `DIS` is 0.3676 (= 36.76). §10.3 discloses exactly this
(*"I quoted the DIS_SURR value in lane D's brief"*), so the document is honest in aggregate — but
§10.1's row is the one a reader takes the number from, and as printed it shows the anchor missing
by 0.14 points under a heading that says all four recompute. Add *"(the charter's 36.9 is
`DIS_SURR`'s 0.3688 — §10.3)"*.

---

### D20 (MINOR) — "identical out-of-fold R² = 0.2355"

`N3_plus_pool.R2_oof = 0.235453`; `IIDmatched_R2_0.235.R2_oof = **0.234915**`. They differ by
0.00054. The ledger also says "identical"; **"matched to within 0.0005"** is the accurate word and
costs nothing. Likewise §12.3's `R2 = 0.16` is 0.16323 and `R2 = 0.24` is 0.23492.

---

### D21 (MINOR, outside my remit but it contradicts §A.7) — §7.1 row 46 reinstates the "third instance" that §A.7 corrects

§7.1 row 46: *"**Third instance of the cross-kind confound this sprint** (S30-L26)."*
§A.7: *"**Second instance of the cross-kind confound this sprint**, after S28-L48… *My first draft
said "third" and counted the widening result, which Appendix A correctly records as circular rather
than cross-kind — a different defect.*"*

This is lane V's **D13**, fixed in A.7 and then reintroduced by §7 (written after V's audit). The
document now states both numbers and explains why one of them is wrong. Change row 46 to "Second".

---

### D22 (MINOR) — §13.1's diagram quotes the **FAIL18** shape share on the **whole production path**

`[2] … its error is **83.1% shape**`. From `s30_F_score.json → predictors_of_rho_pool`:
`shape_err/err_mae` = **3.6304/4.3670 = 83.1% on FAIL18**, **2.1333/2.3386 = 91.2% over all 126**,
**1.8838/2.0005 = 94.2% on the other 108**. The diagram describes the pipeline on every target, so
the figure there should be **91.2% (83.1% on the tail)**. Lane V's D7 also requires the 3.23×
scale-error ratio to travel with any quotation of this row.

---

### D23 (MINOR) — §12.4's "three independent measurements"

The block lists four rows, two of which (`0.0713` cloud and `0.0652` chain) are **the same lane-P
measurement on two bases** — they are `0.3567/5` and `0.3259/5`. So it is **two** independent
measurements (lane P's prior-sign bit and lane T's 3× candidate-indexing result) set against lane
T's value-of-a-bit law. The direction holds; the word "independent" is doing more work than the
table supports.

---

## Things I checked and found clean

Recorded so the sweep is auditable rather than assumed.

- **The sub-MDE sweep over §0, §5, §10, §12.** Beyond D3 and D4, I checked every ratio-bearing
  number in my sections against its artefact's `x_mde`. The −0.3259 arm (2.05×), the −0.5335 arm
  (2.29×), the meter's ladder rows (all fold CIs excluding zero), the chain preference contrast
  (0.56×, correctly reported as below the floor in §10.5/§6.3) and the CA contrast (1.84×) are all
  quoted at their true ratios or correctly withheld.
- **The basis sweep over §5, §10 and §13.** §5.2 labels both columns; §10.2 labels both; §13.1's
  3.0483 and 3.2105 are the right pair from the right record. **No basis error outside D5.**
- **The `compare` sign-convention class in my sections.** Every preference rate, ρ, R², percentile
  and win rate in §5, §6, §10 and §12 is read in the direction its artefact defines.
  `pref(circ_best)` and `pref(pool member)` are read higher-is-better; `native_pctile` is read
  lower-is-better; `coh` is read higher-is-worse, which is the artefact's direction. **No inversion
  leaked into a claim in my sections.**
- **§0.1 against the charter.** Rows 1–11 are the charter's questions 1–11 in order, unedited.
  §0's preamble quote (*"a well-evidenced ceiling argument… is worth more than a fragile 2.98 Å"*)
  is **verbatim** from `BRIEF.md:81`.
- **§5's circuit-ran claim, attacked rather than accepted.** I tried to show the 9-qubit circuit was
  loaded from S28 rather than run. It was not: `circ_best` is loaded (`z["oracle_circ"]`) but
  `circ_s0` is **regenerated** through `oracle_circuit_ceiling(..., starts=1)`, which runs
  `adam(fg, th0, 300, 0.05)` with an exact parameter-shift Jacobian, and the assertion
  `dev_s0 < 1e-6` is what proves it reproduced. The cache's own `meta` carries `"dev_circ_s0": 0.0`.
  The report's sentence is exactly right, including the distinction it draws between "no VQE" and
  "no quantum compute".
- **§13.2's build list.** Every file named exists; the verifier's ledger-path assertion covers nine
  of them and passes.
- **§6.2's honest-limitation paragraph.** *"The sprint-wide total is a lower bound… nobody
  maintained a single register"* is the right admission and is not softened anywhere else.
- **§11.3's four open defects.** All four are real, all four are located, and all four code
  citations hold.

---

## VERDICT on §0

**Not yet.** §0's structure is right and its two headline measurements — `3.2105 Å unchanged` and
`−0.3259 Å from five ORACLE signs, 2.05× MDE, 5/5 folds` — are exact, correctly attributed and
correctly labelled ORACLE. **The section fails on its labels, not on its arithmetic**, and three of
the four failures are on the sentences it exists to deliver:

1. the deliverable is called **native-free** and is **ORACLE** (D1);
2. the one "deployable" effect is quoted at **1.84× its real size**, against a percentage computed
   from the real one, and is called deployable in a document that twice says it is not (D2);
3. the "**strictly harmful**" half of the central swing rests on a **0.68× MDE** number the report
   calls "measured" and never prints a ratio for (D3);
4. the **ORACLE label falls off the price pair in exactly the two graphics** a reader will lift out
   intact (D4).

Fix D1–D5 and §0 is publishable. **None of them removes a result.** The 0.30 Å swing at matched
accuracy survives all four — it gets *sharper* under D3's correction, because the honest statement
is that the fitted arms sit at zero-to-worse while the matched-R² i.i.d. arms clear their MDE in
the other direction. The `coh` gate survives D1 — it is still cheap, still discriminating, still
unprecedented in this project's record; it is a development-time admission test rather than an
inference-time filter, which is what §12.2 already implies and what the next sprint needs to know
before it tries to run it on an unsolved target.

D6 and D7 are not about the science but they are about §0's credibility as a summary. Two of §0.4's
counts were already corrected during this audit; **"twelve lanes" and the falsifier composition are
the two still standing**, and both are wrong in the direction that flatters the sprint.

**Report state at the close of this audit: 1405 lines, 14:25:25.** All twenty-three defects above
were re-checked against that state; D6's two fixed sub-items are recorded as fixed. §0 lines 29, 66
and 76 — the three that carry D1, D2, D3 and D4 — are unchanged from 14:11.
