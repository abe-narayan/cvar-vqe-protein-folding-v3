# EXPERIMENT -> MATH : the interface I will consume

Written 2026-09-06 by the EXPERIMENT workstream, before any Sprint-18 run.

I own the 126-target continuous-torsion refinement test. The BRIEF tells me to consume MATH's
callable and **not** to build a second, subtly different definition of the degree-1 object. MATH
had not landed a module when I started, so I have built a **PROVISIONAL** implementation in
`s18/exp_anova.py` behind an adapter. **The moment `s18/math_*.py` exposes the interface below,
`exp_anova.py` will import it and the provisional code becomes dead.** If MATH's definition
differs from the provisional one, MATH's wins and I re-run.

## The interface

```python
build(dhat, sd, i, j, pool_phi, pool_psi, seed) -> obj
obj.E_full(phi, psi)  -> (f, grad_2n)      # the deployed distance objective, unchanged
obj.E_le1(phi, psi)   -> (f, grad_2n)      # first-order functional ANOVA under mu
obj.E_res(phi, psi)   -> (f, grad_2n)      # residue-additive (adds intra-residue weight-2)
obj.E_lambda(phi, psi, lam)                # = E_le1 + lam * (E_full - E_le1)
obj.E0                                     # the ANOVA mean
obj.argmin_le1()      -> (phi, psi)        # EXACT separable argmin, no search
```

## The provisional definitions (what MATH must confirm or overturn)

Reference measure `mu` is a **product measure over the 2n individual torsion coordinates**,
each factor being the **empirical marginal of that coordinate over the target's own shipped
top-75 retrieval pool** (`s14.avgspace.top75_windows`). This is native-free and
target-conditioned -- the BRIEF's "most informative legitimate choice". A uniform-torsion
sensitivity arm is also run.

    E0        = E_mu[E]
    f_i^phi   = E_mu[E | phi_i] - E0                      (n one-dim functions)
    f_i^psi   = E_mu[E | psi_i] - E0                      (n one-dim functions)
    E_le1     = E0 + sum_i f_i^phi(phi_i) + sum_i f_i^psi(psi_i)      STRICT weight-<=1
    g_i       = E_mu[E | phi_i, psi_i] - E0               (n two-dim functions)
    E_res     = E0 + sum_i g_i(phi_i, psi_i)              RESIDUE-ADDITIVE
    E_res - E_le1 = the intra-residue weight-2 mass       (BRIEF sec 4 point 1)

Conditionals are Monte-Carlo over `mu` on a regular periodic grid, then represented by exact
**trigonometric interpolation** (FFT), which gives an analytic gradient for L-BFGS.

## Two facts MATH should know because they change the experiment

1. **`E_le1` is separable, so its global argmin is EXACT and costs no search.** Coordinate-wise
   argmin of the 1-D tables. This is the direct analogue of the 19-target *certified* degree-1
   argmin, and it is the arm Phase 0's number should be compared against. I report it separately
   from the local-descent arm.
2. **The lambda ladder collapses to a convex mix**: `E_le1 + lam*(E_full - E_le1) =
   (1-lam) E_le1 + lam E_full`, so lambda=1 must reproduce s17's `refine_full` bit-for-bit.
   That is my internal validity check.

## What I need from MATH

- the lattice-equivalence check (uniform mu on the enumerated k=4 lattice == Walsh weight-<=1)
- confirmation or correction of the mu choice and of the coordinate-vs-residue split
- F5: is the continuous mapping valid at all
