# S32 CAUSAL MAP — where RMSD is lost

Charter §55 calls this the backbone of the sprint. Built by the coordinator from artefacts, not from
prior reports' prose. **Every rung is on the BUILT CHAIN** unless it says otherwise; the CA cloud is
an intermediate and the set mean is a third object.

**Status: v1, hours into the sprint.** Rungs marked ⚠ are quoted from S29's O-ladder and are being
independently re-derived by lane V (contract rule 4 — a cross-sprint quotation is the class of
statement that has been wrong before). Rows marked **OPEN** are what the lanes are measuring.

---

## The ladder, as one object

```
                                                    chain     loss at this arrow
sequence
  │  BLOSUM62 retrieval, K = 500
  ▼
K=500 pool ......................... best member   1.7078
  │                            best sparse s=10    1.1139   ⚠  <- the pool ALREADY contains this
  │  DISTOGRAM SCORE prefix, 500 → 128  (128 = 2^7, the VQE register: "the quantum field of view")
  ▼                                                          +0.4350   57% order statistic,
                                                                       43% score WORSE than random
top-128 ............................ best member   2.1435
  │  prefix, 128 → 75
  ▼                                                          +0.1620
top-75 ............................. best member   2.3055
  │  UNIFORM COORDINATE AVERAGE  (the shipped readout)
  ▼                                                          +0.7378   (to cloud 3.0483)
average CA cloud ................................. 3.0483
  │  multi-start ideal-geometry projection, λ-ladder (0.0 → 0.3); the CANONICAL arm is λ = 0.3
  ▼                                                          +0.1622
BUILT CHAIN (PRODUCTION) ......................... 3.2105
```

**Total headroom already inside the existing pool: 2.10 Å.** Nothing below the first rung adds
information; every arrow after it only chooses or destroys.

---

## Arrow by arrow

### 1. sequence → K=500 pool  (retrieval)

| | |
|---|---|
| **transformation** | BLOSUM62 window matching against the fragment bank |
| **information in** | the target sequence only |
| **information out** | 500 real deposited structures, each a valid chain |
| **destroyed** | everything about the target not expressible as sequence similarity |
| **assumption** | that sequence similarity retrieves structural similarity |
| **known** | the best-matching window has **12% identity** — *structure and sequence are decoupled at this length*; yet BLOSUM beats random at **every** K, and a perfect scorer over the retrieved set caps at 1.838 Å |
| **could a better solution disappear here?** | Yes, and irreversibly — but **the pool it produces already supports 1.1139 Å**, so this is not where the *current* 3.21 Å is lost |

### 2. K=500 → top-128  (the DISTOGRAM SCORE PREFIX — *not* retrieval)  — +0.4350 Å, **57% of it an order statistic**

| | |
|---|---|
| **transformation** | keep the 128 best by the shipped L1 Bayes-risk **distogram score**. `s29_O_ladder.py` calls this *"rung 9: the top-128 prefix (**the quantum field of view**)"* — **128 = 2⁷, the VQE register width**. Retrieval is the earlier arrow (universe → K=500); calling this one "the retrieval filter", as I first did, attributes a score-stage loss to the retrieval stage |
| **destroyed** | 372 candidates, among them the ones that made 1.7078 reachable |
| **assumption** | that the score orders by structural quality |
| **known** | the score's **in-band** skill is ≈ 0: ρ_in-band **−0.0262** (DIS) against ρ_global +0.1176. Its global correlation is outlier rejection, not ranking |
| **MEASURED (lane V, 2000 draws, cloud, ORACLE)** | the +0.4350 is **+0.2477 pure set-size order statistic** (a *random* 128 of the 500 loses that much) **+ 0.1872 the score ordering performing WORSE THAN RANDOM** (1.06× MDE, 72W/54L, Type-M inflated ~1.10×). The 128→75 arrow splits the same way: +0.0907 size, **+0.0696 ordering, also wrong-signed** |
| **the mechanism, and it is clean** | over the pool the score lifts mean candidate quality 4.4533 → 3.5847 and **halves the spread** (sd 1.3206 → 0.6384), but the **5th percentile is unchanged** (2.6098 → 2.6184). ***It concentrates on the mode and buys nothing in the good tail***, so a random 128 samples the pool's fat left tail and the score's 128 does not |
| **why that cuts both ways at arrow 4** | by S32-L2 the readout rewards **low ⟨w,a⟩ AND high spread**. The score improves the first and **destroys the second**. Which dominates at the endpoint is **not** decidable from this measurement and is a live experiment (lane P) |
| **the survival statistic** | the K=500 best member survives into the top-128 on **63/126 targets — a coin flip** — at mean rank **170/500** |

### 3. top-128 → top-75  (the prefix)  — +0.1620 Å

| | |
|---|---|
| **known, and closed** | the per-target prefix length does not transfer: ORACLE global `m*`=72 is **−0.0044 (0.19× MDE)**, leave-fold-out `m` is **+0.0075, 63W/63L — a literal coin flip**. And a matched **random-subset** family reaches **141%** of the prefix family's gain on the chain: *the prefix ordering loses to an arbitrary 7-bit index*, because prefix variants are nested (lag-1 autocorr 0.917 vs 0.112) so a less-correlated family has a larger per-target minimum |
| **reading** | this arrow is an order statistic, not a lever. **Do not re-open the `m` axis.** |

### 4. top-75 → average CA cloud  (the readout)  — **+0.7378 Å to the cloud; the hull contains 1.11**

| | |
|---|---|
| **transformation** | uniform mean of 75 aligned structures |
| **exact objective** | for any `Σw = 1`: `‖Σ w_x W_x − t‖² = ⟨w,a⟩ − ½ w'Bw` (verified 1.66e-11), with `a_x = ‖W_x − t‖²` **ORACLE** and `B_xy = ‖W_x − W_y‖²` **native-free** |
| **what that says** | minimising wants **low ⟨w,a⟩** (good candidates) and **high w'Bw** (spread — the variance-cancellation term). `B` is conditionally negative semidefinite on `{Σv=0}`, so the program is **convex on the simplex** and tuning-free |
| **so the only unknown is `a`** | and the native-free channel captures **3.99%** of the ORACLE cross term. Setting `a` = const reduces it to "maximise spread", which is quality-blind dispersion: **+0.1436 Å at 1.22× MDE, WORSE** |
| **why averaging is not simply bad** | the pool's error is **68% common-mode** (exact identity, 50.7× the i.i.d. prediction). Averaging removes the other 32%; it cannot touch the common part. Set mean 3.5507 → average 3.2105 is that 32% being cashed |
| **OPEN (lane Q)** | **how precisely must `a` be known for `w*` to be useful?** `∂w*/∂a` is answerable analytically, with no circuit. If `w*` is insensitive along the directions a native-free estimator can resolve, the family is closed by derivation. If it is sensitive only on a low-dimensional subspace, *that subspace is the minimal information requirement* |

### 5. average cloud → built chain  (projection)  — **+0.1622 Å, and part of it is arithmetic**

| | |
|---|---|
| **transformation** | L-BFGS-B over (φ, ψ) minimising CA-RMSD to the cloud plus `λ·pen(φ,ψ)`, from four generic starts, down the λ-ladder **(0.0 → 0.3)**. **The canonical emitted chain is the λ = 0.3 arm**; `fit_ca` (λ=0) is a diagnostic. *An earlier version of this file said the final rung is λ=0 — wrong, and corrected here.* The Ramachandran prior is therefore **already active and already selecting among branches**: `lam_path`'s docstring — *"the prior's real job is choosing among near-degenerate solutions, and it can only do that if the optimiser can see them"* |
| **the cost is a property of the OBJECT, not a constant** | real member **−0.0007 to −0.0030**; sparse combination s=10 **+0.0002**; dense prefix average **+0.1701**; production **+0.1622**. A 75-structure mean has a contracted backbone and is **not a valid chain**; the projection must re-expand it and that repair is **5% of the endpoint** |
| **the degeneracy, from the module's own docstring** | *"A CA trace admits two ideal-geometry torsion solutions at near-equal objective distance, **one Ramachandran-plausible and one not**, and a warm-started optimiser cannot cross between them."* And: *"**THE REFERENCE DISAGREES WITH ITSELF, by up to 1.6 Å** … the structures this stage returns on those targets are **not determined by the objective; they are determined by the arithmetic**"* |
| **the degeneracy survives the prior** | lane V measured, **on the canonical λ=0.3 arm**, that a **one-ULP (7e-15 Å) change in the input cloud moves the built chain 0.10–0.15 Å**; an independent re-projection disagrees with the S29 canonical on **126/126** targets (mean +0.0021, p90 0.026, **max 0.5174**). The stage is perfectly deterministic given bit-identical input and **discontinuous in it** |
| **a penalty is not a selector** | the prior enters as a **term in the objective**, biasing the optimisation. **Post-hoc ranking of converged branches by Ramachandran likelihood is a different operator and is untested.** Hint, same job, paired: λ=0.3 is **+0.0107 Å worse on the chain than λ=0** (3.2148 vs 3.2041) — the prior costs about a hundredth while buying plausibility. Get the MDE before reading anything into it |
| **OPEN (lane R)** | ORACLE-best branch (priced as best-of-N, with split-half transfer); and whether a **native-free chiral** criterion picks it in band |

---

## The thesis this map suggests, stated so it can be killed

**Arrows 2–4 are all the same problem in different clothes: they need `a`, per-candidate quality, and
every native-free estimate of `a` has in-band skill that is zero or the wrong sign.** If that is
information-limited rather than algorithm-limited, no amount of optimisation at those arrows helps,
and S31's conclusion stands.

**Arrow 5 is a different kind of problem, and it is the one place the missing quantity may be
computable.** The branch decision is discrete, genuinely per-target, currently decided by rounding,
and the discriminating variable is **chirality of the torsion solution** — which is precisely what a
distance-map observable cannot see.

> **Every native-free ranking signal this project has ever tested is a function of the distance map,
> and is therefore ACHIRAL by theorem G1. A torsion-branch choice is exactly the degree of freedom an
> achiral observable is blind to. The Ramachandran surface is strongly chiral, native-free,
> target-specific through `res_classes(seq)`, and already implemented in `core/project.py`.**

**That is a mechanism, not a hope, and it is falsifiable three ways:** the ORACLE-best branch may be
close to production (no prize); the prize may vanish once best-of-N is priced (an order statistic);
or the Ramachandran score may have no in-band skill over branches (mechanism absent — which would
contradict the docstring and would itself be worth knowing).

**What this map does NOT claim.** That arrow 5 is the largest loss — it is not; arrow 4 is. It claims
arrow 5 is the largest loss whose missing information might be **native-free and computable today**.
Those are different statements and the report must not merge them.

---

## Rungs still to be filled

| rung | owner | why it is not yet here |
|---|---|---|
| the ladder re-derived from raw artefacts | **V** | contract rule 4; ⚠ rows are provisional until then |
| ORACLE-best branch, and its best-of-N price | **R** | the central unknown of the thesis above |
| in-band skill of a chiral criterion | **R**, **D** | R over branches, D over candidates — *different objects* |
| `∂w*/∂a` sensitivity | **Q** | closes or opens arrows 2–4 by derivation |
| pool quality vs diversity vs common-mode, separated | **P** | charter §10 demands they be measured independently |
| does this ladder have the same SHAPE at 40+ residues | **L** | several core findings are length-suspicious |
