# PREREG — S31 lane B (physical model): free energy (B1) and backbone torsion (B2)

**Written before the first number.** Committed at the timestamp in the commit that adds this file.
Lane B, charter §7A, §7B, §14. Contract rules 19, 28, 29, 30 govern this document.

**Endpoint basis for everything below.** B2's measurements are made in **pool-member CA-RMSD to
native** (`universe['rr']`, a CA point cloud, per candidate), which is *not* the endpoint. The
endpoint is the **built-chain mean, 3.2105 Å**; the CA point-cloud endpoint is **3.0483 Å**. I will
not quote an endpoint number unless I run the endpoint, and if I do I will say which basis.

---

## B1 — FREE ENERGY. The theoretical pre-check, registered before any thermodynamics is computed.

### B1.0 The question as the charter puts it

> Can a physically motivated free-energy or thermodynamic observable supply target-specific
> information not reducible to sequence plus the static pool?

### B1.1 What `S` is a function of — the question I was asked to answer plainly

Let `U` be the shipped potential (`core/amber.py`: `app.ForceField("amber14/protein.ff14SB.xml",
"implicit/gbn2.xml")`). For a candidate structure `x` and inverse temperature `β`, with `B(x)` the
basin containing `x`:

```
F_β(x) = -β⁻¹ ln ∫_{B(x)} e^{-βU(y)} dy
E(x)   = ⟨U⟩_{B(x)}
S(x)   = (E - F)/T = -∂F/∂T
```

**Arguments of `S`:** `x` — a member of the pool; the **sequence**, which enters only through the
atom typing that instantiates `U`; and `β` and the basin map `B`, which are universal constants and
a universal construction. **There is no fourth argument.**

So `S = Φ_{β,B}[U_seq](x)`: a universal functional evaluated at a pool member. By lane P's repair of
the source enumeration — *a universal function evaluated on a target-specific argument yields
target-specific output but has no target argument of its own* — **`S` is not a third source.** It is
an operator on (sequence, pool), exactly like `E`.

That answers the *source* question. It does **not** by itself answer the *usefulness* question,
because every cost in the pipeline is also a function of (sequence, pool). The usefulness question
is settled separately, below.

### B1.2 Lemma B1 (to be verified numerically, not assumed)

> **Lemma B1.** The shipped potential is reflection-invariant: `U(Rx) = U(x)` for `R` the point
> reflection `x ↦ −x`.

*Argument.* Bond, angle, Lennard-Jones and Coulomb terms are functions of interatomic distances,
which `R` preserves. Every `PeriodicTorsionForce` term has the form `(V_n/2)[1 + cos(nφ − γ)]`;
under `R` every dihedral `φ ↦ −φ`, and `cos(−nφ − γ) = cos(nφ − γ)` **iff** `γ ∈ {0, π}`. ff14SB
carries no CMAP term (CMAP is ff19SB), so there is no non-even torsion contribution. GBn2's Born
radii are functions of interatomic distances. ∎ *(subject to the phase table actually containing
only 0 and π — that is F2 below.)*

### B1.3 Corollary B1 — the closure

> **Corollary B1.** `F_β`, `E`, `S` and every temperature derivative of `F` are rotation-,
> translation- **and reflection-**invariant single-structure observables. By **G1** (S30, classical
> MDS) each is therefore a **function of the candidate's pairwise distance map**, and so lies inside
> the class S30 closed by theorem and measured empty on 43 channels (S30-R: ordering survives on
> 2/43; preference fails on **43/43**).

Novelty in `S` does not rescue the free energy, because **`S` is achiral for the same reason `E`
is.** The enthalpic half was already closed by measurement (AMBER at the 51st percentile on matched
pools, significantly worse than chance as an accuracy ranker); the entropic half is closed by the
same theorem that closed contact topology and Rg.

> **Corollary B1′.** The same argument closes the elastic-network / normal-mode / landscape-curvature
> family named in charter §14: an ANM or GNM Hessian is built from pairwise distances, so its
> spectrum, its log-determinant (the harmonic configurational entropy `S_harm = const − ½ ln det′ H`)
> and every spectral-graph observable derived from it are distance-map functions.

**Consequences registered here, so they cannot be retrofitted:**

- *Basin populations*, *ensemble reweighting*, *temperature-dependent ranking*: a softmax or a
  monotone reweighting of single-structure free energies is a transform of a distance-map function.
  Closed by Corollary B1.
- *Local configurational entropy*, *per-residue entropy*: closed **twice** — by Corollary B1 and by
  S30-R's locality result (local features ΔR² −0.089 against global +0.600, so a sum of per-residue
  terms cannot see a lever arm).
- *Entropy/enthalpy decomposition*: splitting one distance-map function into two distance-map
  functions.

### B1.4 What is NOT closed, named honestly

**(i) Genuinely multi-structure free energy.** G1 is a theorem about *single-structure* observables.
A barrier `F‡(x_a → x_b)` between two pool basins is not a function of `D(x_a)` and `D(x_b)` — it
depends on the path — and escapes G1. **This is left open on PRICE, not on theory:** path sampling
or umbrella sampling over pool pairs is 10³–10⁵ trajectories per target against a
one-AMBER-process-at-a-time rule on a 16.7 GB box. I am not closing it and I am not running it.

**(ii) The chiral part of a chiral force field.** If Lemma B1 fails — e.g. under ff19SB's CMAP —
then `F` has an odd part that escapes G1. **That odd part is carried by the backbone torsion term.**
So the free energy's only theoretical escape hatch *is* B2, which is why this lane spends its compute
there.

### B1.5 A registered directional prediction about restrained free energy

A free energy restrained toward the distogram prior is, **by construction**, coherent with the
prior's error, which S30 measured to be coherent with the pool's common mode. Contract rule 29: at
identical out-of-fold R² = 0.2355 a coherent corrector emits **+0.0554 Å worse** and an i.i.d. one
**−0.2466 Å better**. **Registered prediction: a restrained-vs-unrestrained free-energy ranker is
harmful, not neutral.** I am not running it; I am recording the sign so that a future lane that runs
it has a pre-registered expectation to falsify.

### B1.6 Falsifiers for B1, in the units I will measure

| id | falsifier | unit | closes/opens |
|---|---|---|---|
| **F1** | `max over sampled conformations of \|U(x) − U(Rx)\| / (\|U(x)\| + 1)` **> 1e-6** | relative energy, per term group and total | Lemma B1 false as implemented → the free-energy route escapes G1 through its chiral part → B1 reopens |
| **F2** | ff14SB's `PeriodicTorsionForce` phase table contains any `γ ∉ {0°, 180°}` (tolerance 1e-9 rad) | radians | Lemma B1 false algebraically → B1 reopens |
| **F3** | the same check on the **Legacy** potential's torsion term shows it is *even* in `(φ,ψ)` | relative score | would mean the Legacy torsion channel is achiral and **B2 closes immediately by G1** |

**F1 and F2 can only remove my own excuse** (contract rule 22): they are run before, not after, the
conclusion is written.

### B1.7 Prior odds, stated honestly before the measurement

- Lemma B1 holds numerically on the shipped code: **12 : 1 for**. The algebra is not in doubt; the
  implementation could carry a chirality restraint, a CB improper, or a guard that breaks it.
- A single-candidate free-energy observable is a channel materially new relative to the 43 already
  audited: **1 : 25 against**.
- B1 ends this sprint closed by derivation rather than by compute: **4 : 1 for**.

### B1.8 What I will therefore not spend

No thermodynamic integration. No quasi-harmonic entropy over 126 × 75. No basin-population ranking.
No temperature-dependent candidate ranking. Each is a distance-map function by Corollary B1 at a cost
of hours to days per arm. **The compute goes to B2.**

---

## B2 — BACKBONE TORSION

### B2.1 The decomposition — the mechanism, stated before the numbers

For an ideal-geometry backbone (fixed bond lengths and bond angles, as
`core/geometry.build_backbone_batch` builds), the point reflection `R` acts on torsions **exactly**
as `(φ, ψ) ↦ (−φ, −ψ)`, because `R` preserves all bond lengths and bond angles and negates every
dihedral. Hence, for any torsion observable `T`:

```
T_even(φ,ψ) = ½[ T(φ,ψ) + T(−φ,−ψ) ]        reflection-INVARIANT
T_odd (φ,ψ) = ½[ T(φ,ψ) − T(−φ,−ψ) ]        reflection-ANTI-invariant
T = T_even + T_odd
```

`T_even` is a rotation/translation/reflection-invariant single-structure observable, so **by G1 it is
a function of the distance map** and lies in the closed class.

> **All of a torsion channel's escape from G1 lives in `T_odd`.**

### B2.2 The second axis — separability — and the 2×2 that follows

The channel under test is
`LEG_torsion(x) = 0.15 · Σ_i rama_penalty(aa_i, φ_i, ψ_i)` (`core/energy.py:torsion_term`,
`DEFAULT_WEIGHTS["torsion"] = 0.15`) — a **sum of per-residue terms**, separable by construction and
documented as such in its own docstring. S30-R's locality result on this instrument says a sum of
per-residue terms cannot see a lever arm.

|                       | **separable** (Σ over residues) | **coupled** (non-separable) |
|---|---|---|
| **even** (achiral)    | closed twice: G1 **and** locality | closed by G1 |
| **odd** (chiral)      | closed by locality | **the only open cell** |

**A chiral, non-separable torsion functional is the one cell no existing theorem closes.** That is
what B2 is for.

### B2.3 Hypothesis H-B2, registered

`LEG_torsion`'s **+0.2212 preference contrast (2.14× MDE, 5/5 folds**, S31 meter sweep,
`s30/results/s30_D_meter_sw_LEG_torsion_chain.json`) is **coarse-triage skill only**, and:

1. it is carried predominantly by the **odd** part (the Ramachandran basins `_RAMA_BASINS` sit at
   negative `φ` with the α_L basin at depth 0.40 against α_R's 1.00 — a strongly chiral function);
2. its **in-pool** selection skill is ≈ 0, because the pool is assembled from real PDB fragments, so
   the Ramachandran channel is **saturated** inside the pool: its within-pool dispersion is small
   relative to the triage gap it was measured on.

If both hold, the meter sweep's positive on this channel is explained without any in-pool
information, and **B2 closes on the separable cell** — leaving only the coupled-chiral cell.

### B2.4 Why a torsion quantity could contain information the pool common mode does not

Required by contract rule 28 and by the charter's "do not invent physical justification after seeing
results". The mechanism I claim, in advance:

The pool's common mode is a **shared displacement of the candidates' distance maps** — it is an error
in *where the residues are*. A chiral, non-separable torsion functional measures *handedness of the
relation between two chain segments* (e.g. the signed crossing of two contacting segments, or the
sign pattern of coupled torsion transitions). Two facts make that plausibly incoherent with the
common mode:

- **Chirality is a sign, not a magnitude.** The common mode is identified in S30-L7 as a *mean
  displacement* `μ_{t,p}` of pair distances — an additive, achiral quantity. An observable that is
  odd under reflection is orthogonal to every even observable **on a reflection-closed ensemble**,
  and the common-mode displacement is even.
- **Handedness has a fixed biological ground truth.** Right-handed β-sheet twist and right-handed
  α-helices are sequence-independent, fold-independent facts. A functional that reads them is not
  reading the pool: it is reading a universal constraint the pool members may violate
  *independently* of their shared displacement.

**This is a mechanism, not a result.** Its falsifier is M5 below and I expect it to fail.

### B2.5 What will be measured — all on existing caches, no new physics compute

Data: `s8/generate_univ/*.npz` (126 targets, ≤13k candidates each with `W`, `PHI`, `PSI`, ORACLE
`rr`, `fold`), `s12.instrument.pool_idx` (the shipped K = 500), `shipped_record(pdb)["sub"]` (the
shipped top-75), `s27/cache/*.npz` (all 30 channels precomputed on the 500).

| id | measurement | control |
|---|---|---|
| **M0** | **mirror gap** — the distribution of `\|rmsd(x,nat) − rmsd(Rx,nat)\|` over the pool, against the pool's own RMSD IQR. Prices how much ranking signal sits in the chiral direction of configuration space. **Diagnostic, not a ceiling** — see the caveat in B2.7 | pool RMSD spread |
| **M1** | **chirality decomposition** of `LEG_torsion` and `RAMA`: variance share of `T_even` vs `T_odd` over the pool, and the in-pool skill of each part separately | — |
| **M2** | **saturation**: within-pool σ of the channel against its triage gap (pool members vs the meter's `RAND_SIGNED` rungs); and the within-pool share of total variance | — |
| **M3** | **in-pool skill**: per-target Spearman ρ(T, ORACLE `rr`) over the 500 and over the shipped top-75, fold-clustered CI, MDE | (a) `DIS`, the shipped cost; (b) a **constant α-helix** zero-information control (contract rule 9 — never uniform-on-the-torus); (c) a within-target label shuffle |
| **M4** | **partialling**: within-target partial Spearman ρ(T, `rr` \| `DIS`) — information not already in the shipped cost | same |
| **M5** | **coherence — the admission test**: `coh` = within-target correlation of the channel's implied residual with the pool common-mode pair error `μ_{t,p}` (S30 §12.2). **ADMIT iff `coh` < 0.6931** | uncorrected 0.6931 |
| **M6** | **the open cell**: build a **chiral, non-separable** torsion functional and run it through M1–M5 | all of the above |

Multiplicity is counted and reported; every ORACLE quantity (`rr`, `μ_{t,p}`, any percentile) is
labelled ORACLE / NOT DEPLOYABLE in the sentence carrying its number.

### B2.6 Falsifiers for B2, in the units I will measure

| id | falsifier | consequence |
|---|---|---|
| **G1-check** | `T_odd` variance share of `LEG_torsion` **< 10%** | the channel is essentially achiral → G1 applies directly → **B2 closes** |
| **G2** | in-pool ρ on the shipped top-75 band fails **1.0× MDE** with a fold-clustered CI containing zero | no in-pool selection skill; the +0.2212 is triage only. **Registered prediction: this fires.** |
| **G3** | partial ρ(T, `rr` \| `DIS`) **clears 1.0× MDE** with the fold CI excluding zero | the channel carries information the shipped cost does not → B2 survives to M6 **regardless of G2** |
| **G4** | M6's chiral non-separable functional clears **both** G3 **and** `coh < 0.6931` | the project's first non-empty escape from G1 → to the coordinator immediately |
| **G5** | the constant-α-helix control matches the informative arm on any of M3/M4 | that arm is measuring chain length or a build artefact, not torsion content |

### B2.7 A caveat I am registering now rather than being caught on later

**G1 is a statement about the novelty of a channel, not about information capacity on this pool.**
The distance map `D(x)` determines `x` up to reflection, and the actual pool contains no mirror
pairs (all L-peptides), so on the real pool an achiral observable is *not* information-limited —
being "a function of `D`" costs nothing there. What G1 forbids is a *new channel*, not a *good
function of the old one*. Accordingly **M0 is reported as a diagnostic of chiral dynamic range and
explicitly not as a ceiling**, and the closures above are closures on *novelty*, which is what the
charter asked about.

### B2.8 Prior odds, stated honestly before the measurement

- `LEG_torsion` is predominantly odd (chiral): **4 : 1 for**.
- `LEG_torsion` has in-pool selection skill on the shipped top-75 at ≥ 1.0× MDE: **1 : 6 against**.
- Any torsion functional I build clears the coherence admission test `coh < 0.6931`: **1 : 12
  against**.
- B2 ends this sprint **closed** rather than open: **3 : 1 for**.
- Lane B produces a deployable endpoint improvement: **1 : 30 against**. The value I expect to
  deliver is a derivation that redirects compute, not an Å.

---

## Multiplicity

Lane B pre-declares its comparison budget: **M0–M5 on two incumbent channels (`LEG_torsion`,
`RAMA`) and M6 on at most three constructed functionals**, each against three controls, on two
bands (500 and top-75). Every comparison emitted is counted in the artefact and reported to the
sprint register. No comparison outside this list is reported as a registered test; anything else is
labelled exploratory.
