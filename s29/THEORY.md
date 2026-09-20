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

---

## 3. THE SPECTRAL CONDITION FOR A MEANINGFUL NON-DIAGONAL HAMILTONIAN

**The question.** S28 found the pool's Gaussian similarity graph near rank one (Perron vector 95
to 99% uniform) and the hopping term's gradient variance decaying at -1.7 to -1.8 per qubit
"because of the SPECTRUM, not the off-diagonality" (S28-L8b, L11). This section derives the
spectrum from the pool's geometry (3a), derives the gradient variance from the spectrum (3b),
measures the centered and agreement matrices on 12 real pools and answers whether the decay goes
flat (3c), and says what the ground state of `diag(E) - J A_c` is and whether it is the
consistency mechanism S27 section 6 measured as harmful (3d).

Measurements in this section are mine: `s29/s29_T_spectra.py` -> `s29/results/s29_T_spectra.json`
and `s29_T_spectra_rows.jsonl` (72 spectra cells) / `s29_T_grad_rows.jsonl` (216 gradient cells),
job `s26/jobs_done/s29T_spectra.json` (exit 0, 100.2 s, peak RSS 0.337 GB), 12 targets (S27's T11
set, `P.targets()[::11][:12]`), native-free: the module reads no native and no RMSD. The hop-only
gradient cells reproduce S28's `s28_B_train.json :: hop_only|J1` **exactly** at n = 4 to 8
(4.0612e-03, 7.2844e-04, 2.8167e-04, 5.8983e-05, 3.9483e-05) and to 7% at n = 9 (4.4466e-06 vs
4.1567e-06; S28's n = 9 register orders the 500 candidates differently under the same graph).

### 3a. Why `lambda_2/lambda_1` is 0.11 to 0.14, from the band

`A_ij = exp(-d_ij^2 / 2 sigma^2)` off the diagonal, `sigma` = median off-diagonal `d`, normalised
to unit spectral norm (`s27/s28_B_hop.py :: kernel_graph`). Write `a_ij = exp(-d_ij^2/2 sigma^2)`,
`abar` its mean, and expand about the mean exponent. For a pool whose members are a cloud
`x_1..x_M` in shape space (`d_ij = |x_i - x_j|` to the accuracy of the RMSD superposition),

    d_ij^2/2 sigma^2 = ( |x_i|^2 + |x_j|^2 - 2 <x_i, x_j> ) / 2 sigma^2,
    A = abar (1 1^T - I) + Delta,      Delta ~= -abar [ u 1^T + 1 u^T - (1/sigma^2) X X^T ],    (3.1)

with `u_i = (|x_i|^2 - mean)/2 sigma^2`. The rank-2 part mixes into the Perron direction; the
non-uniform structure is the pool's Gram matrix `X X^T`. Hence

    lambda_1(A) ~= abar M,   lambda_2(A) ~= (abar/sigma^2) lambda_1(X X^T) = (abar/sigma^2) M s_1^2,

    **lambda_2 / lambda_1  ~=  s_1^2 / sigma^2  =  f_1 * <d^2> / (2 sigma^2)**,                 (3.2)

where `s_1^2` is the largest eigenvalue of the pool's shape covariance, `f_1 = s_1^2 / sum_k s_k^2`
is its share of the pool's total shape variance, and `<d^2> = 2 sum_k s_k^2` is the mean squared
pairwise distance. **The ratio is the top shape mode's share of the pool's variance, discounted by
how much of the distance distribution the kernel width spans**, and nothing else to first order.

Measured on the 12 pools at the deployed register (n = 9, the 500-window pool): `f_1 = 0.324`,
effective shape dimension `(sum s^2)^2/sum s^4 = 6.11`, `<d^2>/2 sigma^2 = 0.534`, so (3.2)
predicts `0.173` against a measured `lambda_2/lambda_1 = 0.138` (median over 12, range [0.104,
0.211]) -- the S28 number reproduced (S28-L8b: 0.138 at n = 9) and predicted from the pool's own
geometry to 25%. The prediction over-estimates because `d_max/sigma = 2.56` is not small and the
exponential's curvature damps the far tail. (Convention note: I order eigenvalues by `|lambda|`,
so at n <= 6 "lambda_2" is sometimes the most negative eigenvalue and the ratio prints negative on
half the targets; the n >= 7 rows are unambiguous.)

> **The reading.** `lambda_2/lambda_1 ~ 0.13` is not a defect of the Gaussian kernel or of `sigma`.
> It says the pool has about 6 effective shape modes with the leading one carrying a third of the
> variance, and any similarity kernel that is a smooth decreasing function of `d` will have the
> same ratio to first order, because it is a property of the POOL, not of the kernel. Tuning
> `sigma` moves `<d^2>/2 sigma^2` and nothing else: to reach `lambda_2/lambda_1 = 0.5` one needs
> `sigma^2 ~ 0.65 <d^2>`, i.e. `sigma` about 1.15x the RMS pairwise distance, at which point every
> off-diagonal entry is within 20% of 1 and the graph is a constant plus noise. **The band is a
> trap: wide `sigma` flattens the spectrum by making the graph uninformative.**

### 3b. Gradient variance from the spectrum: the stable rank is the whole story

Let `F(theta) = <psi(theta)| A |psi(theta)>` on the real RY/CNOT ansatz. Each parameter enters
through `psi -> exp(theta_k Gt_k) psi` with `Gt_k` real antisymmetric and `Gt_k^2 = -I/4` (the
generator is `-i Y_q/2` conjugated by a Clifford; `g_kk = 1/4` exactly, S13 section 11). Then

    dF/dtheta_k = psi^T (A Gt_k + Gt_k^T A) psi = <psi| [A, Gt_k] |psi>,                        (3.3)

exact, with `[A, Gt_k]` symmetric. For `psi` distributed as a real 2-design on `S^{D-1}`,
`D = 2^n`, and any symmetric `M`, `E[psi^T M psi] = tr M / D` and
`Var[psi^T M psi] = 2 (D tr M^2 - (tr M)^2) / (D^2 (D+2))`. With `tr[A, Gt] = 0`,

    Var_theta[ dF/dtheta_k ] ~= 2 ||[A, Gt_k]||_F^2 / D^2,
    ||[A, Gt]||_F^2 = 2 tr((A Gt)^2) + tr(A^2)/2   (using Gt^2 = -I/4),
    => **Var[dF/dtheta] ~= tr(A_0^2) / D^2 = r_stable(A) / D^2 when ||A||_2 = 1**,               (3.4)

with `A_0 = A - (tr A/D) I` the traceless part (the identity has no gradient) and the **stable
rank** `r_stable = ||A||_F^2 / ||A||_2^2 = sum_a (lambda_a/lambda_1)^2`. The cross term
`2 tr((A Gt)^2)` is the only place the observable's BASIS enters; it is bounded by `tr(A^2)/2` and
averages small over the ansatz's generators, which is why (3.4) carries no diagonal-versus-
off-diagonal term at leading order.

**This derives and sharpens S28's mechanism.** S28 said the decay is "the near-rank-one spectrum's,
not off-diagonality's" and supported it with a diagonal same-spectrum control decaying at the same
rate (S28-L11 caveat (b)). (3.4) says why: the variance depends on `A` only through `||A||_F^2`,
which is basis-independent, so a diagonal matrix with the same spectrum must decay identically.
It also converts the observation into a design rule, because `r_stable` is a free design variable:

| observable, unit spectral norm | `r_stable` | predicted slope of `log2 Var` per qubit |
|---|---|---|
| near-rank-one dense kernel (S28's A) | `1 + (l2/l1)^2 + ... ~ 1.04` | -2 (measured -1.83) |
| a `k`-regular sparse graph, `A/k` | `M/k` | -1 (S28-B2 measured -1.00 / -1.03 at k = 5 / 10) |
| the rank-standardised diagonal energy `diag(E)/max E` | `D/3` (E has sd 1, range +-sqrt 3) | -1 (S25 measured -0.649 at depth 3) |
| a single Pauli string or a permutation | `D` | 0 |

Three checks with no free parameter. (i) **S28's Gaussian graph**: `r_stable = 1.036` at n = 9 gives
`Var = 1.036/512^2 = 3.95e-6` against the measured `4.16e-6` (S28) and `4.45e-6` (mine): **within
7%**. At n = 4, `1.117/256 = 4.36e-3` against `4.06e-3`: within 7%. (ii) **S28-B2's kNN graph**, for
which the entry recorded a mechanism-free "30 to 60x the variance at n = 9": a `k`-regular graph
has `||A||_F^2/lambda_1^2 = M/k`, so `Var = M/(k D^2) = 500/(10 * 512^2) = 1.9e-4`, i.e. **46x** the
Gaussian's -- the measured factor, predicted. (iii) The `J` at which the hopping term's gradient
variance equals the diagonal terms' is `J* = D sqrt(Var_diag / r_stable)`; with
`Var_diag = 3.05e-2` (S28-L8b) and `r_stable = 1.04`, `J* = 512 * 0.171 = 88`, against lane D's
independently computed `J* = 85.7` (S28-L11 item 1).

> **Design rule (the spectral condition).** An off-diagonal term is visible to the circuit's
> gradient at the deployed width if and only if its **stable rank grows with the register**:
> `Var_hop / Var_diag = r_stable / (D^2 Var_diag)`, so parity at `n = 9` needs
> `r_stable ~ 8000`, i.e. `r_stable ~ 16 D`. **No dense similarity kernel can do this, centered or
> not**; a `k`-regular sparse graph reaches `r_stable = M/k` (parity would need `k ~ 0.06`, below
> one neighbour: unreachable); and **any Gram matrix of structural deviations is capped by its
> rank, `r_stable <= rank <= 3 N_res - 6 <= 42` on this instrument**. The conclusion is not
> "off-diagonal Hamiltonians cannot be trained" -- it is that at `D = 512` **no pool-geometry
> operator can be gradient-dominant at `J = O(1)`**, and the honest options are (a) accept
> `J ~ 60 to 90`, at which the hopping term IS the energy and CVaR acts on a perturbation,
> (b) shrink the register (at `D = 64` parity needs `r_stable ~ 120`, above a Gram cap of 42 but
> within reach of a sparse graph), or (c) put the compatibility structure in the ENERGY (diagonal,
> `r_stable = D/3` by construction) rather than in a hopping term.

### 3c. The centered graph and the signed agreement matrix, measured on 12 real pools

Three matrices on the same sub-pool of the `2^n` best DIS candidates, each at unit spectral norm:
`A` (S28's Gaussian kernel), `A_c = H A H` with `H = I - 11^T/M` (double centering: the uniform,
typicality mode projected out exactly; measured uniform overlap of the top eigenvector `1.4e-31`),
and `G = Delta Delta^T / N_res`, the **signed agreement matrix** of the members' deviations from
the pool mean in the medoid frame.

    n = 9 (the deployed register), median over 12 pools    A        A_c       G
    lambda_2 / lambda_1                                  0.138    0.465    0.634
    lambda_3 / lambda_1                                  ~0.10    0.276    0.346
    r_stable = ||.||_F^2 / lambda_1^2                    1.036    1.591    1.675
    effective rank  (sum l^2)^2 / sum l^4                1.07     2.35     2.35
    participation ratio of the top eigenvector / dim     0.905    0.514    0.255
    top eigenvector's squared overlap with uniform       0.969    0.000    0.007

    hop-only gradient variance (n = 4..9, 120 theta draws, seed 1009, depth 3, the S28 law)
    matrix   n=4        n=5        n=6        n=7        n=8        n=9      slope   pred slope
    A      4.061e-03  7.284e-04  2.817e-04  5.898e-05  3.948e-05  4.447e-06  -1.830    -2.020
    A_c    9.590e-03  3.112e-03  6.011e-04  1.203e-04  2.196e-05  3.599e-06  -2.305    -2.512
    G      3.141e-03  1.375e-03  3.892e-04  9.907e-05  2.093e-05  5.051e-06  -1.900    -1.962
    (measured / predicted by (3.4), median over the 18 cells: 0.90; range 0.31 to 2.50, the
     extremes at n <= 6 for A_c, where D <= 64 and the 2-design approximation is worst)

> **THE ANSWER TO 3c, AND IT CORRECTS A READING ALREADY IN `s29/STATE.md`.** Centering the graph
> raises `lambda_2/lambda_1` by 3.4x (0.138 -> 0.465) and the signed agreement matrix by 4.6x
> (-> 0.634): in *that* sense the S28 degeneracy is removed, and the coordinator's integration
> note 2 is right that S28 closed one similarity measure and not the class. **But the
> gradient-variance decay does NOT go flat, and it does not even improve: the centered graph decays
> at -2.30 per qubit against the raw graph's -1.83, and the agreement matrix at -1.90.** At n = 9
> both sit within a factor of 1.5 of the uncentered Gaussian (3.60e-6 and 5.05e-6 against
> 4.45e-6). The mechanism is (3.4): the gradient sees the **stable rank**, which rises only from
> 1.04 to 1.59 / 1.68, because after centering the second mode is merely comparable to the first,
> not because the spectrum has become extensive. A matrix with `lambda_2/lambda_1 = 0.63` and a
> fast tail still has `||A||_F^2 = O(1)` at unit spectral norm.
>
> **Consequence for lane B, before it builds:** the centered or agreement Hamiltonian at `J = O(1)`
> is as gradient-invisible at `n = 9` as S28's was (`J* = 512 sqrt(0.0305/1.6) = 71`, against 88
> for the raw graph: a 19% reduction in the required `J`, not an order of magnitude), and the Gram
> cap `r_stable <= 3 N_res - 6` means no re-weighting of the agreement matrix fixes it. What is
> genuinely new in `A_c` and `G` is the **meaning of the ground state** (3d), not the trainability;
> a build justified by "the spectrum is no longer degenerate" is justified by the wrong number.

### 3d. What the ground state of `diag(E) - J A_c` is, and the symmetry that cancels it

For `G = Delta Delta^T`, the quadratic form is `<v|G|v> = |sum_i v_i delta_i|^2 / N_res`: the top
eigenvector is **exactly the signed combination of members whose deviations from the pool mean add
up to the largest displacement**, i.e. the pool's first principal shape mode, and
`sqrt(lambda_1(G))` is that displacement's length. For `A_c` the same holds to first order by
(3.1): the centered kernel is `-(abar/sigma^2) X X^T` with the uniform mode projected off, so its
top eigenvector is the same principal contrast up to sign conventions. So:

> At large `J`, the ground state of `diag(E) - J A_c` (or `-J G`) is the pool's **principal
> contrast**: amplitudes positive on the members that deviate one way along the pool's first shape
> mode and negative on those that deviate the other way. It is a coherent set only in the signed
> sense; in probability it is **bimodal over the two poles of that mode**.

**The symmetry that cancels it, and it is fatal for a `p`-based readout.** Every deployed readout
(R1 the CVaR tail, R2 the `p`-weighted average, R3 the `p`-top-m) is a function of `p_i = psi_i^2`
and is therefore blind to the amplitude's sign. Write the members as
`W_i = cbar + s_i |v_i| delta + eps_i` with `s_i = sign(v_i)` along the principal mode. The
`p`-weighted average is

    sum_i v_i^2 W_i  =  cbar + delta ( sum_{s_i>0} v_i^2 - sum_{s_i<0} v_i^2 ) + O(eps),        (3.5)

and the bracket is the **imbalance between the two poles**, which centering makes small by
construction (`v` is orthogonal to the uniform vector). **A centered off-diagonal Hamiltonian
consumed by any probability readout is self-cancelling to first order: it selects both poles of
the mode it discovered and averages them back to the pool mean.** The residual imbalance is
whatever `diag(E)` -- the distogram -- contributes, which is no new information (section 2).

Two corollaries, both cheap to check and both binding on lane B's design:

* **B must use a signed readout** (`C = sum_i psi_i W_i / sum_i psi_i`, S28 lane A's) or break the
  pole symmetry explicitly. S28 refuted that readout *for accuracy under the shipped objective*
  (+0.23 to +0.26 A, S28-L26b) -- but that was a refutation of the OBJECTIVE, not of the readout
  (S28's own reading: given sign freedom, the objective uses it to move away from the native).
  Pairing a centered Hamiltonian with a signed readout is legitimate and new; pairing it with a
  `p` readout is provably null by (3.5).
* **The whole family collapses to a one-parameter ORACLE question.** Under the signed readout the
  ground state emits `production +- eta * PC1(pool)` for a scalar `eta`. So the build's accuracy
  ceiling -- before any circuit, any `J`, any CVaR -- is
  `min_eta mean_targets RMSD(c + eta PC1, native)` with sign and step chosen ORACLE. **Lane O can
  measure this in minutes** from pools it already has, and it prices the whole build. My prior, and
  it is assumption (A3) of theorem 2 again: PC1 of the pool is a *within-pool* direction, i.e. it
  lives in the 32% idiosyncratic component that averaging already removes (S23 L9), and is by
  construction nearly orthogonal to the 68% common-mode error the answer needs.

**Is this the S27 section-6 mechanism?** Not the same, and the distinction matters. S27's
consistency channels (CONS, DMAP_CONS, POOLGO) select the pool's **mode** -- the lowest-dispersion
subset -- which by the set-mean law (`operator-consumes-set-mean`) moves the retained set's mean
toward the pool's centre, and was measured harmful (+0.286 A for CONS, S27 section 6). The centered
Hamiltonian selects the **extremes of a contrast**, which moves the set mean *away* from the pool
centre along `+-PC1`. It is a different operator and the S27 closure does not transfer. What does
transfer is the law behind both: the terminal operator consumes the retained set's mean, so the
only question is whether the displacement's SIGN is right -- and section 2 says the marginals do
not supply a sign. The predicted outcome is a wash (half the targets improve, half worsen, mean
inside MDE) unless an `n`-carrying channel picks the pole.

### 3e. Prediction (checkable in minutes; lanes B, O, D)

1. **The decay (lane B, the measurement it was spawned to make).** On the centered graph `A_c` and
   the agreement matrix `G`, the hop-only gradient variance at n = 4..9 has log2 slope
   **-2.3 +- 0.3 and -1.9 +- 0.3 respectively, NOT flat**, and at n = 9 both are within a factor of
   1.5 of the uncentered Gaussian's. Falsified if either slope is above -1.3, or if either n = 9
   variance exceeds 3e-5. Reproducible in 100 s: `python s29/s29_T_spectra.py --grad`.
2. **The stable-rank law (lane B or D).** For ANY unit-spectral-norm observable the project builds,
   `Var[dF/dtheta_0] = r_stable/D^2` within a factor of 2 at `n >= 7`, and the implied
   `J* = D sqrt(0.0305/r_stable)`. Predicted for `G` at n = 9: `Var = 6.4e-6`, `J* = 71`. Falsified
   by any observable at `n >= 7` whose measured variance departs from `r_stable/D^2` by over 3x.
3. **The pole symmetry (lane B, before any endpoint run).** With a `p`-based readout the emitted
   structure of the `J -> infinity` ground state of `diag(E) - J A_c` differs from the pool mean by
   less than the projection floor on at least 80% of targets. With the SIGNED readout it differs by
   `eta PC1`, `eta` the amplitude imbalance.
4. **The family's ORACLE ceiling (lane O, minutes).** `min_eta mean RMSD(c + eta PC1(pool), t)`
   with the best global `eta` chosen leave-fold-out: predicted **under 0.15 A** better than
   production, with the per-target sign an order statistic priced by `best_of_k_within`. If it is
   above 0.30 A my prediction fails and the B build is worth considerably more than I think.

---

## Q1. THE FLAT MINIMISER SET, AND WHY THE OPTIMISER REACHES THE OPTIMUM WITHOUT MOVING THE STRUCTURE

(The coordinator's question of 2026-09-20 00:2x, taken ahead of sections 1, 4, 5, 6, 7; it is
section 6's material, so section 6 is answered here and not repeated later.)

### Q1.1 The flat directions are exact, and they are most of the simplex

Barkoutsos et al.'s statement -- if `p(x*) >= alpha` for the optimal basis state `x*` then
`CVaR_alpha(E; p) = E(x*)`, its global minimum, so every such `p` is a global minimiser -- is the
special case at the optimum. The operative statement at the **production point** is stronger and
needs no overlap assumption. From the exact derivative (`core/quantum.py :: cvar_exact`, envelope
theorem),

    dCVaR_alpha/dp(x) = (E(x) - q)/alpha  on the strict tail,  0 for every x with E(x) > q,      (Q1.1)

`q` the VaR. **CVaR is exactly constant along every simplex direction supported above the VaR.**
At the deployed cell the realised tail is `m = 74.1` of `D = 512` (`s27/results/s28_B_rows.jsonl`,
seed 0, J = 0), so the CVaR term is blind to

    D - m - 1 = 437 of the 511 simplex directions  =  85.5%                                      (Q1.2)

at every point, not only at the minimiser. Barkoutsos's global-minimiser set
`{p : p(x*) >= alpha}` is, by comparison, exponentially rare in state space -- for a Haar-random
real state `P(p(x*) >= 0.18) ~ (1 - alpha)^((D-3)/2) = e^-46` at `D = 512` -- and the deployed
state is nowhere near it (`p_max ~ 1/PR = 1/428 = 0.0023 << 0.18`). **The flatness that matters
here is (Q1.2), not the overlap condition**: the objective is indifferent to how the mass above
the VaR is arranged, which is exactly the freedom an averaging readout would have to consume.

### Q1.2 What the entropy term selects inside the flat set: derived exactly

The project's CVaR is the LOWER tail, whose Rockafellar-Uryasev form is a maximum,
`CVaR_alpha(E; p) = max_t [ t - (1/alpha) sum_x p_x (t - E_x)_+ ]`. The bracket is linear in `p`,
concave in `t`, and the simplex is compact convex, so Sion's minimax theorem applies and

    F* = min_p max_t [...] = max_t [ t - T log sum_x exp( (t - E_x)_+ / (alpha T) ) ],
    p*(x)  proportional to  exp( (t* - E_x)_+ / (alpha T) ).                                     (Q1.3)

(The code asserts the alpha = 1 limit: (Q1.3) returns the Gibbs free energy `-T log sum e^{-E/T}`
and the Boltzmann law, reproducing S25 section 6.3 to 6 decimals,
`s29/results/s29_T_reach.json :: entropy`.)

Read (Q1.3): **`p*` is exactly uniform on every state above the VaR -- the 437 flat directions,
which the entropy term resolves by flattening them -- and rises exponentially below it with scale
`alpha T` in energy units.** On the deployed rank ladder the rank spacing is
`Delta E = 2 sqrt(3)/D = 0.00676`, so the enhancement per rank is `exp(Delta E/(alpha T)) = 1.078`
and the total enhancement across the prefix is 8.4x. Exact numbers at the deployed cell
(`n = 9`, `alpha = 0.18`, `T = 0.5`, rank ladder):

    state                          F         m (realised tail)   entropy    PR
    uniform over 512            -4.5394            92            9.000     512.0
    the EXACT simplex optimum   -4.7237            29            8.819     342.3      (t* = -1.534,
                                                                                       prefix 29)
    the deployed trained circuit -4.5610            74.1          8.836     407        (S28 rows,
                                                                                       seed 0, sd 0.022)

So the optimum is a **near-uniform state with a modest enhancement on a short prefix** -- the
brief's phrase is right in substance, with the correction that it is not "uniform ON a prefix" but
"uniform above the VaR, enhanced below it". The circuit sits between the uniform state and the
optimum, closing about a third of the `m` distance (92 -> 74 of 92 -> 29) and a fraction of the
free-energy gap; it is under-trained toward a target that is itself nearly uniform.

### Q1.3 Therefore: a target-independent rank-weight profile, and one scalar channel

`E = zrank(score[top[:D]])` is the standardised rank ladder on every target, identical up to
tie-averaging (worst deviation 1.18% of range on 8 targets, S25 L17). `p*` in (Q1.3) is a function
of `E`, `alpha` and `T` only. Hence:

> **The deployed quantum stage is, up to tie structure, a FIXED weight profile over RANKS,
> `w(rank; alpha, T)`, applied to each target's own sorted candidate list.** It has no other target
> dependence. The emitted structure is `sum_i w(rank_i) W_i`, whose only target-specific input is
> which candidate the distogram put at which rank -- i.e. no information the classical sort does
> not already have. S25 L17's "two trained states in the whole deployment" is this, derived.

Two measured corroborations, both native-free: the realised `m` over 126 targets at the deployed
cell has sd 6.74 and correlation **+0.026 with chain length** (seed 1: 71.2, sd 7.94, -0.104) (my computation over
`s27/results/s28_B_rows.jsonl`, seed 0, J = 0), against a between-cell spread of 29 to 92 as
`(alpha, T)` moves; and the tail is the classical top-`m` prefix to 1.1e-13 on 4,914 cells
(S28-L21).

**This closes the loop the sprint has been circling.** The optimiser reduces the objective on
126/126 (S28-L18b) and the structure does not move, because:
(i) 85.5% of the objective's directions are flat and the entropy term sets them to uniform;
(ii) the non-flat directions are the tail's internal weights, and the deployed readout R1 consumes
only the tail SET;
(iii) the tail set is a prefix of the energy order by the set-equality theorem, so it is fixed by
the ordering the classical sort already produced;
(iv) the single remaining scalar the optimisation can move is **where the prefix cuts**, `m`, and
`m` is a function of `(alpha, T)` alone.
**The deployed CVaR-VQE is therefore equivalent, at the endpoint, to choosing one number `m`** --
and the `m`-ladder has been priced three times (S22's ladder, S27 T5, S28's `s28_B_mladder.json`).

### Q1.4 Falsifiable clauses for lane D's meter

* **(M5) FLAT FRACTION.** For any proposed objective, report the fraction of readout-relevant
  directions along which it is exactly constant at the production point. For the deployed CVaR it
  is 437/511 = 85.5% in simplex space and, after the readout, everything except one scalar.
  **A candidate objective is worth a build only if its flat fraction is materially below that AND
  the non-flat directions are ones the readout consumes.** Cheap: the flat set of a CVaR-family
  objective is computable in closed form from the realised tail; for a general objective it is the
  rank of the Jacobian of (objective -> readout-visible coordinates).
* **(M6) THE FIXED-PROFILE CONTROL, which is stronger than "a classical equivalent".** Replace the
  whole quantum stage by the target-independent weight profile `p*(alpha, T)` of (Q1.3) applied to
  the target's own rank order -- no circuit, no optimiser, no per-target computation at all.
  **Prediction: the deployed arm's emitted structure equals this control to within the built-chain
  input floor (S28-L43) on at least 120 of 126 targets.** Any new formulation that claims the
  quantum stage contributes something must break this control; if it does not, the stage is a
  lookup table indexed by rank. Minutes to run, and it needs no VQE.

---

## Q2. THE UNTESTED NON-CLASSICAL CELL: A FREE ENERGY OVER A NON-COMMUTING HAMILTONIAN

Lane L's condition (S29-L13): a formulation is non-classical in the relevant sense only if
**(C1)** the Hamiltonian's terms do not commute, so its eigenbasis is not the computational basis,
**and (C2)** the prepared object is not an eigenvector, so an eigensolver is not its classical
counterpart either. The question: can that cell contain anything measurable *here*?

**The derived answer is NO for the candidate-index encoding with a pool-geometry coupling, and the
reason is a new second-order result, not an appeal to the earlier nulls.**

### Q2.1 A zero-diagonal coupling's thermal state is a classical reweighting at second order

Let `H = diag(E) - J A` with `A` a pool-similarity or agreement matrix, `A_xx = 0` (every graph in
the record has a zero diagonal by construction, `s27/s28_B_hop.py :: finish_graph`). The measured
distribution of the Gibbs state is `p_x = <x| e^{-H/T} |x> / Z`. Expand in `J` (Duhamel):

    <x| e^{-H/T} |x> = e^{-E_x/T} [ 1 + (J/T) A_xx + (J^2) sum_{y != x} A_xy^2 g(E_x, E_y; T) + ... ]

and **the first-order term vanishes identically because `A_xx = 0`**, so

    p_x  proportional to  e^{-E_x/T} [ 1 + J^2 sum_y A_xy^2 g(E_x, E_y; T) + O(J^3) ],           (Q2.1)

    g(E_x, E_y; T) = ( 1 - e^{-(E_y - E_x)/T} - (E_y - E_x)/T ) / (E_y - E_x)^2,  g -> 1/(2T^2) as E_y -> E_x.

> **The entire quantum content of a thermal state over a zero-diagonal coupling, at leading order,
> is a classical reweighting of the Boltzmann law by the candidate's own squared-similarity degree
> `sum_y A_xy^2` -- a native-free scalar feature computable without any circuit.** And that feature
> is precisely what S28's RAND control holds fixed: the degree-matched random graph (Sinkhorn-
> scaled to the same row sums) changed nothing at any `J`, readout or seed (F3/F4 silent on 72
> contrasts, S28-L41). The second-order cell has, in effect, already been controlled for.

### Q2.2 The three conditions the cell would have to satisfy, and which one is binding

(a) **Information.** The coupling must carry something outside the marginals -- it must break (A4)
of Theorem 2. Every pool-geometry operator (similarity, centering, agreement, disagreement between
the prior and the pool) is a function of the pool and the posterior, so it does not. **Binding.**
(b) **Trainability.** By section 3b the coupling is gradient-visible at `D = 512` only if its
stable rank grows with the register; a dense kernel has `r_stable ~ 1.6`, and any Gram of
structural deviations is capped at `3 N_res - 6 <= 42`. At `J ~ 70 to 90` the coupling IS the
energy and the free-energy term becomes a perturbation on an eigenvector problem -- which
violates (C2). Binding jointly with (a).
(c) **Readout.** By section 3d a `p`-based readout is blind to the sign structure that is the only
new content of a centered coupling, and cancels it by the pole symmetry. Binding unless the
readout changes.

**And a fourth, from the record rather than from theory:** even if the cell is entered, S25 L15
measured that this readout cannot resolve a distributional difference of 45% of the mass
(KL 0.902 nats, TV 0.453, endpoint 0.24x MDE). A thermal correction of order `J^2 A^2` is far
smaller than that. **Readout slack alone kills the measurability of the cell in the candidate-index
encoding.**

### Q2.3 What would have to be true instead, and the smallest formulation that qualifies

The obstruction in (b) is specific and it points at exactly one operator class. My variance law
says `Var[dF/dtheta] = r_stable/D^2` at unit spectral norm, so an off-diagonal term is trainable
only if `r_stable` grows with `D`. A dense kernel cannot; **a sum of local Pauli terms can, by
construction**: for a transverse field `M = sum_q X_q`, `||M||_F^2 = n D` and `||M||_2 = n`, so

    r_stable( sum_q X_q ) = D/n,   Var[dF/dtheta] ~= 1/(n D),   slope = -1 per qubit exactly,
    J* = D sqrt( Var_diag / r_stable ) = sqrt( n D Var_diag ) = 11.8 at n = 9.                   (Q2.2)

So the only non-commuting operators that are simultaneously (C1) non-commuting, gradient-visible at
the deployed width, and cheap are **local mixers**, not pool-geometry couplings -- and a local
mixer is only meaningful when the basis states have local structure, i.e. **in a configuration-
space encoding (lane X's), not in a candidate-index encoding** where a qubit is a bit of an
arbitrary label and `X_q` flips a candidate's index bit.

> **The smallest formulation that qualifies:** basis state = a per-residue configuration assignment
> (fragment or torsion bin); `H = H_diag(configuration posterior, 1- and 2-body) + Gamma sum_q X_q`;
> the prepared object a **free-energy / thermal state** (`F = CVaR_alpha - T H(p)` with `Gamma > 0`,
> or a variational Gibbs state), never an eigenvector; readout = the CVaR tail's coordinate
> average. It satisfies (C1) and (C2) by construction, satisfies (b) by (Q2.2), and satisfies (c)
> because the tail's ENSEMBLE, not a signed amplitude, is what is consumed.
>
> **Its cheapest falsifier, in this order.** (1) `Gamma = 0` must be a genuinely different answer
> from `Gamma > 0` on the SAMPLED distribution -- TV against the classical Boltzmann law of
> `H_diag` above 0.45, because below that the readout provably cannot resolve it (S25 L15). (2) The
> classical counterpart must be named correctly: for a local mixer it is **not** an eigensolver but
> a classical thermal sampler / simulated annealing over the same configuration space at matched
> evaluations. (3) The endpoint contrast against that sampler, built chain. If (1) fails, the cell
> is empty for this instrument and no build follows.
>
> **What it still cannot do:** by (a), nothing in this cell creates information about the native's
> deviation from typical. It can only change WHICH configurations the tail contains. So the
> honest upside is charter section 17's "a specific, measurable, classically irreproducible
> contribution at unchanged RMSD", not an accuracy result -- unless the configuration space's own
> posterior carries more than the pool's marginals, which is lane X's premise and is measured by
> lane D's meter, not by me.
