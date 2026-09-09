# SPRINT 24 — WORKSTREAM B (RESIDUAL GENERATOR) — PRE-REGISTRATION

Basis discipline: **point cloud on every side of every contrast in this file.** Readout is the
shipped one — top-75 by the shipped Bayes-risk distogram score, superpose on the medoid, uniform
mean of exactly 75 members. Only the MEMBERS change between arms.

Seeds: `s15.seed.stable_rng`. Writes: tmp + `os.replace`. Results: `s24/results/`.

---

## B-0 — THE SUBSTRATE CHECK (no hypothesis, no direction; a machinery measurement)

The residual architecture requires rebuilding each retrieved window from its torsions before a
residual can be applied. An ideal-geometry rebuild of a real window is not that window. **If the
rebuild alone costs RMSD, that cost is the floor of the entire lane** and every residual arm must
buy it back before it buys anything.

    RUN, n=126, complete.   incumbent (real windows)   3.0483   <- reproduces the pinned constant
                            rebuild (torsion rebuild)  3.0524
                            rebuild - incumbent  +0.0041  SE 0.0047  MDE 0.0131  57W/69L   NULL

    mean |rebuilt window - its own real window| = 0.3730 A, and it cancels in the average.

**The substrate is free.** Recorded as a null, not promoted. This is not a directional experiment
and carries no fork list; it is reported because a non-null here would have closed the lane.

---

## B-1 — RESID-0: WHAT IS IN THE RESIDUAL CHANNEL AT ALL, AND WHAT KIND OF ACCURACY BUYS IT

**Why this runs before a model is trained, and before Lane A's corpus exists.** The lane's premise
is that a learned `(dphi, dpsi)` correction can move the emitted cloud where selection, aggregation,
routing, readout and repair have all failed to (s23 L5/L9/L11). That premise has an ORACLE ceiling
and a *shape*, and both are measurable today with no network, from the labels already in the
universe files. If the oracle residual's gain lives entirely in a component no native-free model can
supply, the lane is closed by arithmetic and no training run is needed to close it.

**The derivation this tests, stated before the run.** Write the retained member as `w_k = t + e_k`
and the residual as `d_k`. The readout is a uniform mean, so it sees `mean_k d_k` and nothing else
to first order. Therefore:

  * a residual whose per-member variation is **i.i.d.** contributes ~0 to the emitted cloud — it
    cancels exactly as s23 L9 says idiosyncratic member error cancels;
  * only the **common component** `dbar` moves the answer;
  * so a *stochastic* generator's stochasticity cannot itself move the mean, and the lane reduces to
    whether the CONDITIONAL MEAN residual points against the 68% shared bias.

That is a prediction, and B-1 measures it directly by decomposing the oracle residual into exactly
those two parts and pushing each through the deployed readout.

### Arms (all n=126, all point cloud, all the same 75 members, only the torsions change)

    P0    incumbent, real windows                                       [reference]
    P0R   rebuild, no residual                                          [THE BASELINE for B-1]
    ORAC(a)   member torsions + a * d*_k ,  d*_k = wrap(native - member)     ORACLE
    COMM(a)   member torsions + a * dbar  ,  dbar = circular mean of d*_k over the 75  ORACLE
    IDIO(a)   member torsions + a * (d*_k - dbar)                             ORACLE
    COH(s)    oracle direction corrupted by angular noise s, ONE draw shared by all 75 members
    IID(s)    oracle direction corrupted by angular noise s, drawn INDEPENDENTLY per member
    RAND(s)   zero-information residual of the same magnitude, no oracle direction   NULL

`a` in {0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0}; `s` in {10, 20, 30, 45, 60, 90} degrees.

`COH` / `IID` are the **matched-accuracy corruption null that project memory makes mandatory for
every learned corrector** (`error-coherence-decides-correctors`: at identical 0.688 sign accuracy,
coherent mistakes emitted +0.31 A and i.i.d. mistakes -0.14 A). They are the arms that price how
accurate the model must be, and in which direction its mistakes may point.

`RAND` is the plausible zero-information control in the operator's own space
(`zero-information-control-must-be-plausible`): a residual of matched magnitude and no direction,
NOT uniform-on-the-torus and NOT the absence of a residual.

### Rule 0 — SIX OPERATOR FORKS, each naming the alternative not taken

    functional     DECLARED the shipped Bayes-risk distogram score, unchanged, selecting the SAME
                   top-75 indices for every arm, fixed before any residual is applied.
                   NOT TAKEN re-scoring after the residual, which would confound the residual's
                   effect with a change of membership and reopen the flat lever of BRIEF SS1.3.
    basis          DECLARED point cloud on every side. Residual arms are compared against P0R
                   (rebuilt, no residual), never against P0, so the rebuild operator is held fixed.
                   NOT TAKEN comparing a torsion-rebuilt arm to the real-window incumbent, which
                   would charge the residual for the 0.156 A basis operator.
    readout        DECLARED uniform mean of exactly 75 members, medoid frame, as deployed.
                   NOT TAKEN re-optimising m for the residual arms, which would let the readout
                   absorb the effect.
    normalisation  DECLARED the residual applied in raw wrapped-angle space, with the sin/cos
                   parameterisation tested against it rather than assumed (the directive proposes
                   sin/cos; `s14.retprior.circ_mean` supplies a genuine circular mean for `dbar`).
                   NOT TAKEN assuming the sin/cos parameterisation is correct because the directive
                   names it.
    null           DECLARED RAND(s), a matched-magnitude directionless residual, PLUS the matched-
                   accuracy COH/IID pair.  NOT TAKEN "no residual" as the only null, which cannot
                   distinguish a correction from a perturbation of the same size.
    THE LABEL      DECLARED continuous Ca-RMSD of the emitted cloud, and every oracle arm carries
                   the word ORACLE at every appearance.  NOT TAKEN a torsion-space MAE or any
                   sigma/coverage summary, which project memory has three times shown does not
                   price emitted RMSD (`prior-mae-prices-selected-rmsd`).

### Hypotheses and their falsifiers

    H1  the residual channel has real headroom through the deployed readout: ORAC(a) falls steeply
        in a.
        Falsifier: ORAC(1) failing to approach the native-torsion rebuild floor (~0.35 A), which
        would mean the readout destroys even a perfect residual.

    H2  THE ONE THAT DECIDES THE LANE.  The gain is carried by the COMMON component.
        Falsifier: COMM(a) flat while IDIO(a) falls.  If IDIO carries the gain, then the lane needs
        per-member specificity, the i.i.d. part does NOT cancel as derived, and the derivation above
        is wrong -- which I would rather find here than after training.
        The reverse outcome is the one that damages my lane's stated architecture: if COMM carries
        it, a STOCHASTIC generator is decoration and the honest object is a single common
        correction, which is a deterministic corrector and is exactly what
        `error-coherence-decides-correctors` warns about.

    H3  the accuracy bar is reachable.  COH/IID stay below P0R out to some usable sigma.
        Falsifier: both arms crossing P0R by 20-30 degrees.  Project memory measures the achievable
        torsion channel at sigma = 69.7 degrees with phi essentially unpredictable from sequence
        (`phi-carries-no-sequence-signal`: 36.1 vs 36.4 degrees blind).  **If the crossing is below
        ~60-70 degrees the lane is closed by a number this project already owns.**

### The observation that must be recorded before the run, because it is the lane's real risk

`d*_k = psi_native - psi_k` and `psi_k` is FULLY OBSERVED at inference. So
`E[d*_k | features] = E[psi_native | features] - psi_k` **exactly**. A torsional residual model is a
native-torsion predictor with a known offset subtracted; the residual framing adds no information of
its own. What it does add is CONDITIONING the s13 torsion predictor never had -- the retrieved
template's own torsions and the distogram. B-1's accuracy ladder prices what that conditioning would
have to be worth. Written down now so it cannot be discovered after a result and presented as
insight.

---

## B-2 — RESID-1: THE DISTOGRAM CONDITIONING CHANNEL, COMPUTED RATHER THAN LEARNED

Registered now, run after B-1. Descend the shipped distogram objective in torsion space from each
retained member's own torsions for a fixed budget, average the 75 descended members through the
same readout. This is the **best case for the only target-specific conditioning a residual model
has**, obtained without training, and it is bounded above by what any distogram-conditioned learned
residual can do. Forks and falsifiers to be written before that run; the standing warning is
`objective-does-not-rank-the-native` (the native is the distance objective's argmin on 3/126
targets, 36.8th percentile), which predicts this arm turns over with budget.

---

## STATUS

Lane A's permitted corpus, exclusion list and hash **do not exist yet** (`s24/` contains only the
coordinator's two Tier-1 files). No training has begun and none will begin before that artefact
exists. B-0 and B-1 consume labels for EVALUATION ONLY and train nothing.
