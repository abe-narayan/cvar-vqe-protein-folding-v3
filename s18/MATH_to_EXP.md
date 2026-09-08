# MATH → EXPERIMENT : the object, the interface, and one correction that changes your labels

> **ADDENDUM, 2026-09-06, after the branch was closed.** The coordinator closed the degree-1
> branch (F2/F4/F5 fired) and stopped further compute in this lane. Consequences for you:
> the coefficient cache `s18/cache/anova_<pdb>_pool_512_16.npz` is **complete for 40 of 126
> targets, not 126** — `math_iface.build()` will rebuild any missing target on demand
> (≈40–90 s each) but nothing is pre-warmed beyond those 40. The interface, the verification
> and the naming correction below all stand unchanged and are correct; §2 is the reason the
> branch closed and is worth reading even though the experiment is over. Full write-up:
> `s18/math_FINDINGS.md`.

Written 2026-09-06 by the MATH workstream, in reply to `s18/COORD_exp_to_math.md`.
Everything below is native-free. Artefacts: `s18/results/math_lattice.json` (sealed,
`_config_hash 4cb0a8dfe9b5376b`), `s18/results/math_anova_selfcheck.json`, coefficient cache
`s18/cache/anova_<pdb>_pool_512_16.npz` (each carries `complete=True` and the config hash).

---

## 1. The module is live. Import it and delete the provisional tables.

```python
from s18 import math_iface as MI

o = MI.build("1A13")                  # cached, native-free, mu = "pool", S = 512, G = 16

f, g = o.E_full(phi, psi)             # deployed objective + exact torsion gradient (2n,)
f, g = o.E_res(phi, psi)              # RESIDUE-additive ANOVA  <-- the BRIEF's object
f, g = o.E_ang(phi, psi)              # ANGLE-additive ANOVA    <-- strictly coarser
f, g = o.E_ge2(phi, psi, base="res")  # the complement
f, g = o.E_lambda(phi, psi, lam, base="res")
o.E0
o.argmin_res(phi0, psi0)              # EXACT separable argmin, no search
o.argmin_ang(phi0, psi0)
o.with_dhat(d_new)                    # rebuild BOTH truncations for any dhat, no chain builds
```

Your names still work as aliases (`o.E_le1` = `E_ang`, `o.argmin_le1` = `argmin_ang`), but
please switch to the unambiguous ones — see §2, that aliasing is exactly the trap.

**Verified on 1A13** (`python -m s18.math_iface 1A13`):

| check | result |
|---|---|
| `E_full` vs `s17/refine.py::_obj` at the same torsions | **2.3e-13** (float noise) |
| gradients vs central differences, all three arms | rel err **≤ 2.9e-9** |
| `E_lambda(·,1) − E_full`, `E_lambda(·,0) − E_res` | **0.0 exactly** |
| CONTROL D — known additive synthetic recovered | **rel err 2.4e-15, PASS** |

---

## 2. THE CORRECTION. Your `E_le1` is not the continuous form of the 2.411 Å result.

Your provisional file labels the angle-additive object `E_le1` / "STRICT weight-≤1" and the
residue-additive object `E_res`. On the enumerated lattice those two roles are **measured**,
and the labels invite the wrong reading.

**(a) The BRIEF's bridge produces the RESIDUE-additive object.** A functional ANOVA is an
ANOVA under a *named* product factorisation. Take the factor to be the residue, μ = ⊗ᵣ μᵣ on
(T²)ⁿ, and you get `E_res`. **VERIFIED EXACT** on all 19 enumerated targets: the order-1
residue ANOVA under uniform μ equals Walsh weight-0 + weight-1 + *all intra-residue weight-2*
coefficients, worst discrepancy **8.9e-16** relative. The order-1 *qubit* ANOVA equals the
strict Walsh weight-≤1 projection, worst **8.9e-16**. Orthogonality and the variance budget
`Var(E) = Var(PE) + Var(E−PE)` hold to 1e-15. **The bridge is sound; F5 does not fire on the
construction.**

**(b) But the 2.411 Å number is the strict object, and the residue-additive object is null.**
Same 19 targets, certified argmins, ORACLE scoring:

| object | argmin RMSD | vs full, paired [95% CI] | W/L/ties | retained variance |
|---|---|---|---|---|
| full objective | 2.661 | — | — | 1.000 |
| **strict Walsh weight-≤1** | **2.411** | −0.249 [−0.650, +0.093] | 7/5/7 | 0.613 |
| **residue-additive ANOVA** | **2.657** | **−0.004 [−0.380, +0.307]** | **3/5/11** | 0.913 |

(Your sibling `s18/q_anova.py` computed the same two numbers independently and agrees to
15 digits.) So the effect Phase 0 reproduced belongs entirely to the object the bridge does
**not** produce.

**(c) There is no continuous image of the strict object.** The two lattice qubits of a residue
index k-means clusters of the *joint* (φ,ψ) library — for the GENERAL class, states
(58°,30°), (−102°,135°), (−68°,−38°), (−120°,−9°) in unsorted k-means order. The high bit
groups {L-helix, β} against {α, bridge}. They are not φ and ψ, and no relabelling makes them
so. The continuum has no qubits; the lattice has no φ/ψ factorisation. **`E_ang` is a
legitimate, strictly coarser truncation worth running on its own merits — it is not a mapping
of the lattice result, and must not be reported as one.**

**(d) And the strict object is not a property of the physics.** Permuting the four torsion
state labels per residue changes no structure, no energy and no ordering of configurations.
Over 24 stable-seeded relabellings the strict weight-≤1 argmin's 19-target mean is
**2.862 ± 0.135 Å** (range 2.600–3.221); **23 of 24 are worse than the full objective's
2.661**, and the shipped labelling (2.411) sits at the **0th percentile** of its own null.
The residue-additive projection is **exactly invariant** to the same relabelling on every
target. Matched-random control (a random Walsh subspace at the same retained variance) is
being run and will be reported in `math_FINDINGS.md`.

**What to do with this:** run `E_res` as the pre-registered degree-1 arm (it is what the
brief derives), run `E_ang` beside it as a second, coarser truncation, and label the second
one an analogy. Do **not** describe either as "the continuous version of 2.411 Å".

---

## 3. Four structural facts to build into the harness. All EXACT.

1. **The distance objective depends on the n−2 INTERIOR residues only.** `d_ab` depends on
   exactly residues {a+1,…,b−1}; ψ_a rotates about an axis through CA_a and φ_b only moves
   atoms placed after CA_b, so both drop out. Verified to 1e-9 at n = 9 and n = 12.
   **Consequence: fᵣ ≡ 0 for r = 0 and r = n−1 for *both* truncations *and* for the full
   objective.** Neither `argmin_res` nor a full refinement can place the terminal residues;
   the adapter leaves them at the start value. The frozen metric scores them, and for a
   9-mer that is 2/9 of the chain. Any terminal-residue choice is a modelling patch and must
   be reported as one — and it is identical across arms, so it does not confound the
   comparison.
2. **Both truncations are separable ⇒ their global argmin is closed form.** No search, no
   budget, no optimiser. `argmin_res` / `argmin_ang` return the certified optimum.
3. **λ = 1 reproduces the full objective bit-for-bit** (checked, 0.0 exactly), so your
   validity gate passes by construction.
4. **The tables are affine in `dhat`.** `E = Σ_p w_p(d_p² − 2 d̂_p d_p + d̂_p²)`, so
   tabulating the conditional moments `E_μ[d_p | θᵣ]` and `E_μ[d_p² | θᵣ]` (which is what the
   cache stores) rebuilds either truncation for **any** distance vector with no new chain
   builds. `o.with_dhat(d)` costs milliseconds. **That makes the coordinator's α-ladder free
   on the degree-1 arms, and it makes the argmin's sensitivity to d̂ exactly computable.**

---

## 4. The reference measure

Default `mu="pool"`: for each residue independently, resample (φᵣ,ψᵣ) from the target's own
shipped top-75 windows (`s14.retprior.windows(pdb,"top75")`). Native-free, target-conditioned,
and a genuine **product** measure — drawing whole windows would not be one, and the ANOVA is
undefined unless μ factorises. `mu="uniform"` (zero-information) and `mu="rama"` (shipped
leave-fold-out Ramachandran density of each residue's class) are built by the same call with
their own cache keys. Lattice sensitivity across the three is in `math_FINDINGS.md`.

**One thing μ cannot do:** `pool` and `rama` are products over residues but *not* over qubits —
a general law on 4 states is not a product of two Bernoullis, and neither factorises on any of
the 19 targets. A functional ANOVA needs the measure to factorise over the chosen factors, so
under either informative μ the qubit decomposition is not an ANOVA at all: its per-qubit terms
stop being mutually orthogonal, the variance budget stops closing, and "zero the coefficients
above weight 1" stops being the projection. The strict object exists as an ANOVA only under
the uniform lattice measure. That is a second, independent reason it cannot be carried
forward.

---

## 5. Quadrature and resolution

`E_le1` is defined against a **frozen** sample Θ^(1..S) ~ μ drawn from
`stable_rng(pdb, mu, S, grid, "mufreeze", salt="s18math")`, so it is deterministic and
bit-reproducible. Estimator (this is the part that must not be changed):

    h_r(t) = (1/S) Σ_s E(Θ_s with residue r set to t),   Ebar = (1/S) Σ_s E(Θ_s)
    f_r = h_r − Ebar,  E_0 = Ebar  ⇒  E_le1 = Σ_r h_r(θ_r) − (n−1)·Ebar

For an exactly additive E this is **identically** E, at any S — the Ebar terms telescope.
Re-centring each fᵣ on its own mesh mean breaks that and leaves an O(σ/√S) constant behind;
that is the bug Control D caught here. Production setting is S = 512, G = 16 (22.5°), with
common random numbers across the whole mesh so every pair that does not span r cancels to the
last bit.
