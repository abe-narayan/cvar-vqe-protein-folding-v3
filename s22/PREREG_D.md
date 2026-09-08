# WORKSTREAM D — PRE-REGISTRATION, SPRINT 22

**Never edited after data. Corrections are dated addenda appended at the bottom of the relevant
entry, never edits to the original text.** Where an entry documents an AUDIT of data that already
existed before this document was written (the coordinator's L1-L3, computed and posted to
`s22/LEDGER.md` before Workstream D was briefed), it is labelled **AUDIT**, not **PREREG**, and the
falsifier is stated as the standard this audit was held to, not as a prediction made blind. Where an
entry is a genuine forward experiment, it is labelled **PREREG** and the hypothesis/falsifier was
written into the script's docstring (timestamped by the file's first commit-worthy version) before
the numbers in its output were read.

---

## D-AUDIT-1 — THE ROUTING CEILING (s22 LEDGER L1)

**What is under audit.** The coordinator's claim: over 13 deployable native-free arms at n=126, best
fixed arm 3.048, ORACLE per-target routing 2.566, headroom 0.482 Å, with stated concentration
diagnostics (50%/21 targets, 80%/50 targets, median 0.317, corr(headroom, incumbent)=+0.564,
12/126 incumbent-already-best).

**Standard held.** Recompute from the named sources (`s21/results/poolgap.json`,
`s21/results/latentsel.json`) independently. Check every arm's BASIS. Check the arm-set is
actually 13 and not ambiguous. Report agreement or disagreement with the exact figure, not just the
direction.

**Falsifier (of the audit, not of the coordinator's claim).** The audit fails to support the claim if
(a) no 12/13/14-arm reconstruction consistent with the coordinator's stated composition reproduces a
headroom within ~0.05 Å of 0.482, or (b) the concentration diagnostics do not replicate in shape.

**Result.** See `s22/agentD_FINDINGS.md` §1 and `s22/results/d_route_audit.py`. Headroom reproduces
at 0.505-0.506 Å across five plausible reconstructions (12, 13 and 14 arms), robust to the exact
arm-set ambiguity to within 0.002 Å — but this is 0.02-0.03 Å from the coordinator's 0.482, an
unresolved discrepancy. Concentration diagnostics replicate in shape (top21 51.1-51.2% vs 50%, top50
81.5-81.8% vs 80%, median 0.313 vs 0.317, corr +0.60 vs +0.564) but not to the exact digit. **Audit
verdict: the ceiling is REAL and ROBUST to the stated ambiguity in its own construction, but the
coordinator's exact figure is NOT independently reproduced and the precise arm list must be supplied
before 0.482 is quoted again as a specific number rather than "≈0.48-0.51."**

> **ADDENDUM, dated 2026-09-07, same session — RESOLVED.** The coordinator supplied the literal list
> and the source of the discrepancy: his first computation used K=13 (poolgap's 8 arms plus
> `lat_avg75, lat_avg75_rand, lat_argmin, lat_medoid` plus `ship_avg75`, giving 2.566/0.482), his
> second computation dropped `ship_avg75` (K=12, a different oracle), and he quoted the first number
> while running subsequent decompositions on the second — a bookkeeping error, not a different
> arm-selection philosophy, and not the duplicate-arm hazard raised above (`pool_argmin` and
> `avg75_medoid75` were never in his list). Independently reconstructing his exact first-computation
> set (found by chance in `s22/results/routerdata.json`'s `A_plus_latent` configuration, a Workstream
> C artefact that appeared mid-session) reproduces `headroom = 0.48195...` to five significant
> figures. **The audit is now fully closed: 0.482 is confirmed on its own stated construction; 0.505
> (this document's own reconstruction, using `lat_rand1` in place of `lat_avg75_rand` as the fourth
> latent arm) remains a legitimate but DIFFERENT 12-arm menu, and both numbers are correct for what
> they each measure.**

Two identity/near-identity hazards found in the source files themselves (not the coordinator's
routing construction, but load-bearing if his 13-arm list drew on them uncritically):
`poolgap.json`'s `avg75_medoid75 ≡ avg_75` and `pool_argmin ≡ avg_1` to floating-point precision —
both are the SAME quantity computed twice under two names, an instance of BRIEF §9's "an identity
counted as an empirical trial" hazard. Additionally, TWO physically different "pool argmin" /
"pool average" quantities exist across the project's own files: window-basis (`poolgap.py`, direct
retrieved coordinates) and rebuild-basis (`latentsel.py`, rebuilt onto the generative latent's own
torsion manifold). Aggregate basis price, measured independently here: pool_argmin rebuild−window
`+0.0528` [SE 0.0418, MDE 0.1172]; pool_avg75 rebuild−window `+0.0155` [SE 0.0102, MDE 0.0286] — the
second reproduces Sprint 21 L20's `+0.016` almost exactly, which is a genuine (if minor) independent
cross-check of L20's number.

---

## D-AUDIT-2 — THE CORRECT NULL FOR A PER-TARGET ROUTER OVER CORRELATED ARMS (s22 LEDGER L1's null)

**What is under audit.** The coordinator's within-arm permutation null returned 1.115, better than
the observed oracle (2.566), and he diagnosed this as mis-specified (destroying the shared
difficulty factor lets the min draw from easy targets), concluding "the ceiling is legitimate, the
estimate needs held-out folds."

**Standard held.** Derive rather than assert. State what null WOULD be well-specified, if any, and
say plainly whether the ceiling needs one at all.

**Derivation and result.** Full derivation in `s22/agentD_FINDINGS.md` §2 (Duty 3's mandated
derivation). Headline: **the ceiling requires no null.** `oracle_routing = E_t[min_a X_{a,t}]` is a
deterministic function of 13×126 real numbers, each itself a deterministic RMSD between a specific
candidate structure and a specific native structure — there is no sampling distribution under which
to ask "is this minimum surprising." A permutation null answers a different, well-defined question
("would this data, with the joint target-arm structure destroyed, produce a similarly low minimum")
and that question's answer (yes, even lower) is a fact about how much the observed HEADROOM depends
on cross-arm correlation, not a validity check on the ceiling. The coordinator's diagnosis of the
MECHANISM is correct: decorrelating arms lets the per-target minimum draw its low value from
whichever arm happens (post-shuffle) to be paired with an easy target, which is easier to achieve
than requiring the SAME 13-tuple to contain a good option on a genuinely hard target — hence the
lower (better) null value. His disposition ("legitimate ceiling; estimate needs held-out folds") is
CORRECT, and is sharpened here into a decision-theoretic argument that needs no permutation test at
all: with zero real per-target signal, the expected value of ANY function from (native-free features)
to (arm choice) is bounded below by the best FIXED arm's own value, by the standard argument that a
decision rule blind to the relevant conditioning variable cannot outperform the unconditional
optimum. That is, `3.048` is already the correct "no-skill" floor for an ACHIEVABLE router, exactly
as BRIEF states ("the correct baseline for a router is the best FIXED arm") — it does not need to be
manufactured by a permutation procedure, and no permutation of the routing table's rows or columns
recovers it as a null, because the ceiling and the floor are different objects answering different
questions.

**Falsifier of this derivation.** If a construction exists under which the ceiling ITSELF (not an
achieved router's estimate of it) is provably inflated by a describable statistical mechanism (e.g.
if the 13 arms were not fixed in advance but selected FROM a larger family using the same 126
targets — a genuine best-of-K-arms selection bias, BRIEF Rule 0 clause 3), then the ceiling DOES need
a correction, and "no null needed" is wrong. **This is checked and applies conditionally**: see
D-AUDIT-2a below.

### D-AUDIT-2a — an addendum found while deriving the above (dated 2026-09-07, same session)

If the 13-arm LIST was itself chosen after looking at which arms perform well on these 126 targets
(rather than being a fixed, pre-specified menu of "everything this project's deployable pipeline can
emit"), the ceiling DOES inherit a best-of-K selection bias at the level of arm choice, on top of the
per-target min. This is not tested here because the audit could not obtain the coordinator's exact
arm-selection procedure (D-AUDIT-1). **Recommendation, adopted as a standing rule for Workstream C**:
report which arms were excluded from the 13 and why, so a reader can judge whether the arm menu
itself was tuned.

---

## D-AUDIT-3 — THE RETRIEVAL/LATENT DIFFICULTY INTERACTION (s22 LEDGER L2/L3)

**What is under audit.** L2's claim that `lat_avg75 - ship_avg75` is monotone and sign-flipping
across quartiles of the incumbent's own (oracle) RMSD, read as "retrieval is better where retrieval
works; the generative latent is better exactly where retrieval fails" — a genuine difficulty-regime
complementarity. The coordinator's own stated objection: stratifying by the incumbent's own RMSD and
then measuring a difference that CONTAINS the incumbent's RMSD as a term could manufacture part of
the effect by regression to the mean.

**Standard held.** Quantify how much of the quartile pattern a null construction with NO genuine
regime-dependent information can produce, using (i) a closed-form decomposition into the mechanical
linear component implied by `corr(lat, inc)` and `beta = OLS slope of lat on inc`, and (ii) an
empirical placebo that destroys the true per-target pairing while preserving each variable's marginal
and its (weak) association with peptide length.

**Result — the objection is CONFIRMED, and it is large.** `corr(ship_avg75, lat_avg75) = +0.819`,
OLS slope `beta = 0.580` (R² = 0.671). For ANY two variables related this way, `d = lat - inc` is
ALGEBRAICALLY forced to trend downward in `inc` with slope exactly `beta - 1 = -0.420`
(`Cov(inc,d) = Cov(inc,lat) - Var(inc)`), independent of whether `lat` carries genuine
difficulty-regime information. Decomposing the four quartile means into this mechanical component and
the residual:

    Q1 (easiest): real +1.469   mechanical +1.193   residual +0.276 [+0.082,+0.464]
    Q2:           real +0.448   mechanical +0.679   residual -0.231 [-0.428,-0.016]
    Q3:           real +0.152   mechanical +0.307   residual -0.155 [-0.354,+0.066]
    Q4 (hardest): real -0.455   mechanical -0.552   residual +0.098 [-0.168,+0.370]

The mechanical component alone reproduces 81-202% of each quartile's real mean, overshooting on 3 of
4. The residual — the genuine information left after removing the single global linear relationship —
is NOT monotone, does not sign-flip in the reported direction, and two of its four CIs include zero.
A placebo (shuffle `lat_avg75` across targets within the same length-decile, breaking true pairing,
preserving marginals) reproduces an EVEN LARGER version of the reported sign-flip shape (+2.25 to
-1.70) from zero real information, confirming the shape is not, on its own, diagnostic of anything.

**Falsifier of THIS audit.** If the residual (after removing the mechanical linear component) had
itself been monotone and sign-flipping with CIs excluding zero at both ends, the objection would have
been REFUTED and L2's regime-complementarity claim would stand as stated. It is not: the falsifier
for the audit's own concern did not fire, i.e., the concern IS supported.

**Consequence, adopted.** L2 is downgraded from "a routing signal, cleaner than anything in Sprint
21" to: *the generative latent's error scales more weakly with difficulty than the incumbent's
(slope 0.58 vs 1 — a single global fact), and the residual beyond that single fact is small, not
monotone, and not established as a usable regime signal.* This also reframes L3: it is not that two
native-free proxies failed to find a real regime split at insufficient precision; there may not be
much regime-dependent structure left to find once the mechanical component is removed, and neither
proxy result should be read as "difficulty-based routing is close but not quite there."

---

## D-PREREG-1 — THE ÅNGSTRÖM CONVERSION OF THE COMPACTNESS CHANNEL (rg_z / rg_gap)

**Hypothesis.** BRIEF §2 names this the campaign's "central conversion": L27/L28 demonstrated
`rg_z`/`rg_gap` as a SIGNAL (partial ρ ≈ 0.34-0.38 against an ORACLE rank-percentile label, n=75,
n≤13 panel only) but explicitly left its Å value unmeasured. If the signal is real and usable, a
held-out-fold router that conditions the choice among ladder arms (top-75 vs a wider or narrower
average, or the medoid) on `rg_z`/`rg_gap` should beat the fixed incumbent (3.048) by a margin
exceeding its own MDE, on the FULL n=126 pool (not the n≤13 subset L27/L28 used).

**Mechanism.** L19's finding that the objective's top window is sometimes a "concentrated wrong
region," and L27/L28's finding that `rg_z`/`rg_gap` (predicted-vs-realised compactness disagreement)
correlates with that regime, jointly predict: on high-|rg_z| targets, trusting the score less (wider
average, or a set-typicality readout like the medoid) should do relatively better than the score's
own top-75.

**Falsifier.** If no combination of {rg_z, rg_gap, |rg_z|} × {2-arm and 8-arm candidate sets} × 3
quantile bins, evaluated by 5-fold held-out CV (bin edges and per-bin best-arm choice fit on training
folds only, applied to the held-out fold), beats the fixed incumbent by more than its own MDE, the
conversion FAILS on this construction and the Å value of the L27/L28 channel remains NOT MEASURED
(not "zero" — only this particular router construction is closed).

**Extension needed first, run and audited before this test.** L27/L28's `rg_disto`/`rg_pool_mean`
were computed only for the n≤13 panel (75 targets), inherited from `d_lrank.json`. Extended here to
the full n=126 pool (`s22/results/d_rg_full126.py` / `d_rg_full126.json`), using the identical
formula (`Rg² = (1/N²) Σ d̂²` from the shipped distogram's posterior-mean predicted distances, plus
the constant trans-virtual-bond term) — this is itself a new, audited artefact, not a re-quote.

**Result.** `s22/results/d_rg_route.py`. Every router construction tried is null:

    router[rg_z    -> avg_75/avg_150]     d = -0.021  SE 0.019  MDE 0.053   NOT MEASURED
    router[rg_z    -> full 8-arm ladder]  d = -0.002  SE 0.022  MDE 0.062   NOT MEASURED
    router[rg_gap  -> avg_75/avg_150]     d = -0.021  SE 0.019  MDE 0.053   NOT MEASURED
    router[|rg_z|  -> avg_75/avg_150]     d = -0.014  SE 0.017  MDE 0.046   NOT MEASURED
    router[* -> avg_75/avg_500]           d =  0.000  (avg_75 wins EVERY training fold, no crossover found)
    router[* -> avg_75/medoid_75]         d =  0.000  (same)

**Disposition: FALSIFIER DID NOT FIRE.** No construction beats the incumbent past its own MDE; two of
six literally never route away from `avg_75` in ANY training fold (no in-sample crossover exists for
those candidate pairs at all, let alone an out-of-sample one). **The Å value of the L27/L28 channel
is NOT MEASURED, and this specific, reasonably thorough attempt to measure it returns null.** This
does not close the channel (L27/L28's own ρ≈0.35 correlation, on its own label and panel, is
unaffected by this result — a coarse 3-bin threshold router is a weak instrument for a ρ≈0.35 signal,
and a continuous/regression-based router, or one using the channel jointly with other features, is
untested). It does mean: **nobody should currently plan a deployable pipeline change around this
channel without a stronger router than the one tested here**, and the "central conversion" BRIEF asks
for remains open, not merely unattempted.

---

## D-LITERATURE — SCOPE STATEMENT

Full write-up in `s22/agentD_FINDINGS.md` §4. Search priorities per BRIEF plus the coordinator's
request (learning-to-defer / selective prediction / mixture-of-experts routing / hybrid
retrieval-generation, as a literature for the L2 "aggregate hides a regime-dependent complementarity"
structure). Conducted via `WebSearch`/`WebFetch`, 2025-2026 emphasis, with explicit
known-method/known-combination/implementation-novelty/empirical-novelty/theoretical-novelty
classification for every claim in this campaign that could plausibly be called novel.

---

## D-MATH — SCOPE STATEMENT

Full derivations in `s22/agentD_FINDINGS.md` §2-3: (i) the router null (this document's D-AUDIT-2,
derived in full there), (ii) the CVaR face / degeneracy-breaking question, (iii) Type-M exaggeration
at this programme's typical sample sizes, tabulated as a function of observed-effect/MDE ratio rather
than asserted as a constant.
