# S29 THEORY (lane T)

Derivations for Sprint 29, one section per question of `s29/briefs/S29T.md`. Sections 2 and 3 are
written first because the coordinator gates every build on them; 1, 4, 5, 6, 7 follow in the
brief's order.

**Conventions used throughout.** A target has `N` residues; a structure is `C in R^{3N}`; the pair
index `alpha = (i,j)` runs over the `P` pairs with `|i-j| >= 2` that the shipped distogram scores
(`core/predict.py :: Distogram.__init__`). `D(C) in R^P` is the map `alpha -> d_alpha(C)`;
`Jc in R^{P x 3N}` is its Jacobian at `C`, with rows `grad d_alpha = (e_i - e_j) (x) u_alpha`,
`u_alpha = (C_i - C_j)/d_alpha` the pair's unit axis. `t` is the native, `c` the production point
(the DIS top-75 uniform coordinate average in its medoid frame, 3.0483 A point cloud,
`s27/results/vqe_rows.jsonl :: DIS`), `W_k` the pool windows, `p_alpha` the predicted per-pair
posterior with median `m_alpha`, `w_alpha` the shipped per-pair weight. Rigid-body components are
removed from every direction, as in S28-L23b. ORACLE quantities (anything computed from `t`) are
labelled in the sentence that uses them. No number below is quoted from prose: each carries its
ledger entry or its artefact.

---

## 2. THE LOCALLY INFORMATIVE CLASS

**The question.** Given only the per-pair marginals `p_alpha` (no joint), which objectives `f(C)`
have positive expected cosine between their steepest-descent direction at the production cloud and
the direction to the native? Measured for the shipped objective: `cos = -0.034` (SE 0.021,
56/126 positive), against a random-direction reference of `0.140`, with no S27 channel above the
reference (S28-L23b). This section says why that is a property of the whole class, not of the
distogram, and which assumption a new information source would have to break.

### 2.1 The class, and the general form of its gradient

Define **class M** (marginal objectives): `f` sees the structure only through its distance map and
sees the data only through the marginals,

    f(C) = Phi( D(C) ; {p_alpha} ),        g := grad_D Phi in R^P,      grad f(C) = Jc^T g.     (2.1)

This is wider than "separable": `Phi` may couple pairs (a learned function of the whole map, a
consensus term against the pool's own map, a triangle-repaired target). Every native-free channel
the project owns that has a gradient is in M: DIS, DIS_MEAN, DISTPOT, CONTACT, CONTACT_LL, CONS,
DMAP_CONS, RG_LAW, RG_UNIV, EXVOL (`s27/ham_lib.py`). The **separable** sub-class M0 has
`g_alpha = w_alpha phi'_alpha(d_alpha(C))` with `phi_alpha` a functional of `p_alpha` alone; the
shipped objective is in M0 with the L1 Bayes risk

    phi_alpha(d) = sum_b p_alpha(b) |d - c_b|,      phi'_alpha(d) = 2 F_alpha(d) - 1,           (2.2)

`F_alpha` the posterior CDF (`core/predict.py:418-422`). Two facts about (2.2) that matter later:
`phi'` **vanishes exactly at the posterior median** and **saturates at +-1**, so the shipped
gradient is a sum of bounded pair forces pulling each distance toward its own median, and a pair
that is 10 A wrong contributes no more force than a pair that is 1 A wrong. `phi'` is a step
function on the 17 bin centres (S25 L6), so those +-1 saturations are reached quickly.

### 2.2 The exact informativeness identity

Let `u = t - c` (ORACLE) and let

    r := Jc u in R^P,        r_alpha = d/ds|_0 d_alpha(c + s u)  ~=  d_alpha(t) - d_alpha(c),   (2.3)

the first-order per-pair signal (exact as a directional derivative; the approximation by the
distance difference is first order in `|u|`). Then, by the chain rule and nothing else,

    <-grad f(c), u> = - g^T Jc u = - <g, r>,                                                    (2.4)

    cos( -grad f(c), u )  =  - <g, r> / ( ||Jc^T g||  ||u|| ).                                  (2.5)

**(2.5) is an identity, not a model.** An objective in M is locally informative at production **iff
its per-pair coefficient vector `g` is negatively correlated (in the `r`-metric) with the
production structure's own signed per-pair error against the native.** Nothing about the
functional form, the normalisation, the temperature, the readout or the ansatz enters. This is the
whole of "the objective's gradient is blind" reduced to one scalar product in pair space.

Two immediate corollaries.

* **C1 (stationarity kills it).** If production were exactly the minimiser over structures of an
  objective in M, then `Jc^T g = 0` and the cosine is 0/0: the gradient carries no direction at
  all. The shipped case is near this: `|grad S| = 0.0488 A^-1` RMS per atom against `|u| = 3.048 A`
  (S28-L23b). The cosine is being read off a small residual, which is why its magnitude is a
  second-order quantity in 2.4 below.
* **C2 (the geometric compression is exact).** `grad f = Jc^T g` annihilates `ker(Jc^T)`, of
  dimension `P - (3N - 6)`. At `N = 12`, `P = 55` and `3N - 6 = 30`: **45% of pair space is
  invisible to any objective in M at the gradient level**, and two objectives whose coefficient
  vectors differ by anything in that kernel have *identical* gradients. Any proposal whose content
  is "make the target map self-consistent / embeddable / triangle-repaired" moves `g` partly inside
  that kernel and buys, to that extent, exactly zero at the gradient. (It can still move the
  objective's *minimiser*; C2 is a statement about the local direction, which is what the cost
  meter's second number reads.)

### 2.3 The typicality decomposition, and the theorem

The record's error model has one shape, stated three ways: every predictor emits the typical
peptide of that length (S18/S19; `sequence-conditioning-hurts-the-failures`), the pool's error is
68% common-mode (S23 L9, exact identity, 50.7x the i.i.d. prediction), and score-selected sources
have parallel biases (cos 0.943, S24 L3). Write, per pair, with `tau` the **typical** map (the
length-matched, target-blind expected distance; operationally the sequence-blind pool's mean map,
which lane O is computing for the H1 probe):

    d_alpha(c) = tau_alpha + a_alpha       a = production's deviation from typical     (small)
    m_alpha    = tau_alpha + b_alpha       b = the posterior's deviation from typical  (small)
    d_alpha(t) = tau_alpha + n_alpha       n = the NATIVE's deviation from typical     (O(1))   (2.6)

and adopt the model's assumptions explicitly:

> **(A1) Marginal class.** `f in M`: (2.1) holds.
> **(A2) Proper risk.** `g_alpha = w_alpha phi'_alpha(d_alpha(c))` with `phi'_alpha(m_alpha) = 0`
> and `phi'` increasing, so `g_alpha = w_alpha kappa_alpha (d_alpha(c) - m_alpha) + O((d-m)^2)`
> with `kappa_alpha = 2 p_alpha(m_alpha) > 0` (twice the posterior density at its own median).
> **(A3) Coherence.** `a` and `b` are jointly distributed with regression slope
> `beta_alpha = cov(a_alpha, b_alpha)/var(a_alpha)`; the record says they are nearly parallel.
> **(A4) No native channel.** Neither `a` nor `b` correlates with `n`: nothing in the system sees
> the native's deviation from typical. `E[a n] = E[b n] = 0`.

Substituting (2.6) and (A2) into (2.4), using `r = n - a` to first order,

    <g, r> = sum_alpha w kappa (a - b)(n - a),
    E[<g, r>] = - sum_alpha w_alpha kappa_alpha ( var(a_alpha) - cov(a_alpha, b_alpha) )
              = - sum_alpha w_alpha kappa_alpha var(a_alpha) ( 1 - beta_alpha ).                (2.7)

> **THEOREM 2 (marginal objectives are informative only to second order, and only about the pool).**
> Under (A1)-(A4), for every `f in M0`,
>
>     E[ cos( -grad f(c), u ) ] = ( sum w kappa var(a) (1 - beta) ) / ( ||Jc^T g|| ||u|| )  +  O(3),
>
> which contains **no term in `n`**. The expected cosine is (i) second order in the *small*
> quantities `a, b`, (ii) proportional to the failure of the posterior's deviation-from-typical to
> match production's, and (iii) **exactly zero when `beta = 1`** -- when the posterior's per-pair
> median map and the production structure deviate from typical in the same way and by the same
> amount. Its realised value fluctuates about that mean with the scale of a random direction, times
> `sqrt(P (1 - rho_ab) / 2) * sigma_a / sigma_n`.
>
> **Corollary 2a.** No function of the marginals alone -- separable or not -- can have an expected
> cosine whose *size* is set by the native's deviation from typical. The only quantity it can
> extract is the disagreement between two estimates of the typical map.
>
> **Corollary 2b (the sign).** `E[cos] < 0` iff `beta > 1`: the objective points away from the
> native exactly when the posterior's median map deviates from typical *more* than a real structure
> can. This is expected here, because `m` is an unconstrained per-pair object and is generally not
> a realisable distance matrix, while `c` is an average of real windows. The measured sign
> (`-0.034`, and `-0.143` on FAIL18 where the posterior is worst, S28-L23b) is the sign of
> `-(beta - 1)`.
>
> **Corollary 2c (non-separable objectives buy nothing by being non-separable).** For `f in M \ M0`
> the same computation runs with `g = grad_D Phi`; the only new freedom is that `g_alpha` may
> depend on the whole map. By (A4) that freedom cannot introduce a term in `n`, and by C2 the part
> of it living in `ker(Jc^T)` is annihilated. A joint, a triangle repair, a learned map-level head
> trained on native-free features -- all remain inside the theorem.

**Why (A4) is the load-bearing assumption, and not a rhetorical one.** It is the compact form of
four independent measurements: no native-free scorer in a 31-member library prefers a 0.25 A
structure to the contracted average once the pool-member control is applied (S28-L48, L49); the
native sits at the distogram's 36.9th percentile in its own pool (S28-L30); all 38 native-free
signals of S12 have no in-band skill except typicality (`consensus-is-the-only-in-band-discriminator`);
and the field has no method at this length either (S29-L1). If some channel *did* correlate with
`n`, (2.7) would acquire a first-order term `-sum w kappa cov(a - b, n)` and the cosine would move
to the `sigma_n` scale -- i.e. by an order of magnitude, not a factor.

### 2.4 Can anything in M have a positive cosine? Yes -- and that is a warning, not an opening

Set `beta < 1` deliberately: replace the target map `m` by the shrunk map
`m(s) = tau + s (m - tau)` with `s < 1`. Then `beta -> s beta` and (2.7) turns positive. **A
positive gradient cosine at production is purchasable with no new information at all, by shrinking
the objective's target toward the typical map.** What it buys structurally is motion toward `tau`,
i.e. *more* typicality and (since averaging already contracts 22%, S23 L1) plausibly a worse
structure. The sprint must therefore treat the meter's four numbers as jointly necessary and
individually insufficient:

> **Warning for rule 19 (lane D's cost-RMSD meter).** The gradient cosine alone can be moved into
> the positive by a one-parameter shrink of the objective's target toward typicality, with zero
> information added and with the emitted structure moving *away* from the native. Any candidate
> objective that gains cosine must be shown not to have gained it this way: report the shrink
> factor `s` implied by its target map (`beta` below), or pair the cosine with the native
> percentile, which the shrink moves the wrong way.

The bound on what genuine incoherence buys. With `sigma_a = sigma_b = sigma` and correlation
`rho_ab`, per-pair weights absorbed, the systematic term of (2.7) has size
`P sigma^2 (1 - rho)` and the fluctuation has sd `sqrt(P) sigma sqrt(2(1-rho)) sigma_n`, so

    E[cos] ~ cos_random * sqrt( P (1 - rho_ab) / 2 ) * (sigma_a / sigma_n),
    cos_random = sqrt( 2 / (pi (3N - 6)) ) = 0.140 at the instrument's mean length (S28-L23b).   (2.8)

At `P = 55`, `rho_ab = 0.94` (S24 L3's score-selected bias cosine, used as a proxy and labelled as
one) and `sigma_a/sigma_n = 0.2` (sequence conditioning is worth 0.776 A of a 3.989 A blind
pipeline, S12 `coord_null`), (2.8) gives `|E[cos]| ~ 0.036`, against the measured `0.034`. The
match of the *magnitude* is the model's one quantitative success and should not be over-read (two
of the three inputs are proxies); the *sign* is Corollary 2b's and is a separate statement.

### 2.5 What a new source would have to break, term by term

| candidate source | which assumption it attacks | what the theorem says it buys | record's price |
|---|---|---|---|
| a second, differently biased pool | (A3): lowers `rho_ab`, `beta != 1` | second order: `E[cos]` scales as `sqrt(1-rho) sigma_a/sigma_n`; at `rho = 0.65` (the *unselected* blind library) `|E[cos]| ~ 0.12`, still at the random reference | S24 L2/L3: independence 31% only when unselected, and then `q = 1.231` (0.76 A worse); score selection restores `rho = 0.943` |
| the pool's own dispersion | (A4)? No: dispersion is the 32% idiosyncratic part, a second moment | magnitude of the error, never its sign; `cos` needs the sign | S23 L9; `in-band-ordering-is-per-target` ("the only leverage supplies the per-target SIGN") |
| a joint `p(D)` replacing marginals | (A1) formally; (A4) not at all | realisability repair moves `g` partly into `ker(Jc^T)` (C2) and buys 0 there; it can move the minimiser (section 1) but not the local direction | S13 locality theorem: the joint is what a torsion parameterisation enforces exactly, and neither energy ranks (`torsion-space-neither-energy-ranks`) |
| an ESM-attention pairwise map | (A4) only if it carries `n` beyond the distogram | first order iff `cov(attention deviation, n) != 0` given DIS | S17: in-band content, length-gated, no shortlist value; the distogram already consumes ESM-2 650M PCA-32 (`core/data.py:886`) |
| a physics term on the emitted structure | (A4), in principle: it does not go through `p` | first order iff it correlates with `n`; it is not in M, so Theorem 2 does not apply to it | S25 L16: both physics energies are measurably worse than a random subset; S13: CVaR at small alpha is a steric-clash filter |
| a learned residual (predict `a - n`) | (A4) by construction, if trainable | first order by definition -- this is the only class that changes the order | S19/`error-coherence-decides-correctors`: a corrector trained on the predictor's own features inherits its error structure and emits +0.31 A at 0.688 sign accuracy |

**The one structural reading.** Everything in the table except the last two rows is an attack on
(A3), and (A3) is worth a *second-order* quantity. Only a channel that sees `n` -- the native's
deviation from the typical peptide -- changes the order of the answer. That is finding 8 in the
charter and the prior's derivative in S24 (`prior-derivative-is-the-only-steep-lever`, -2.15 A per
unit of prior improvement against a few hundredths for every downstream lever), derived here rather
than observed.

### 2.6 Prediction (checkable in minutes)

Three clauses, all on existing artefacts, all ORACLE diagnostics:

1. **The identity.** For every target, `<-grad S(c), u> = -<g, r>` with `g` the shipped per-pair
   coefficients `w_alpha (2 F_alpha(d_alpha(c)) - 1)/P` and `r = Jc u`, to floating-point. Lane D
   has both sides already: `s27/results/s28_A2_cosine_rows.jsonl` holds `grad S` and `u` per target,
   and `s12.instrument.distogram(pdb)` holds `risk`, `grid` and the pair list. Predicted agreement:
   relative error < 1e-6 on 126/126.
2. **The sign law (Corollary 2b).** Compute `beta` per target as the ordinary-least-squares slope
   of `b = m - tau` on `a = D(c) - tau`, with `tau` the sequence-blind pool's mean distance map
   (lane O's H1 artefact) or, failing that, the corpus mean map at that length. **Prediction:
   `beta > 1` on at least 2/3 of targets, and `sign(cos_DIS) = -sign(beta - 1)` on at least 70% of
   targets**, with the exceptions concentrated where `|beta - 1|` is smallest. Falsified if `beta`
   has median <= 1, or if the sign agreement is inside a coin-toss CI.
3. **The shrink experiment (the warning).** Re-score the shipped objective with the target map
   shrunk to `tau + s(m - tau)` for `s in {1.0, 0.75, 0.5}` and re-read the meter's cosine at
   production. **Prediction: the cosine rises monotonically with decreasing `s` and crosses zero
   near `s = 1/beta_median`, while the native percentile (meter number 3) gets *worse* over the
   same grid.** This is the cheapest available demonstration that meter number 2 is gameable, and
   it costs one re-score of an existing risk table.

Who checks: **lane D** (it owns `s29/s29_D_cost_audit.py` and the four numbers), using
`s27/results/s28_A2_cosine_rows.jsonl`, `s12/cache/disto_*.npz` and lane O's blind-pool map. Clause
1 is minutes; clauses 2 and 3 are one pass over 126 cached distograms.
