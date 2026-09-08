"""SPRINT 14, ENER-4 -- how much of each energy rides on states CA-RMSD CANNOT SEE?

`core/project.py` leaves three torsions per chain inert for the CA trace: `phi[0]`,
`psi[n-1]` and `phi[n-1]`.  The k-state variable of residue n-1 selects BOTH `phi[n-1]` and
`psi[n-1]`, so **residue n-1's entire state is invisible to CA-RMSD** -- yet `phi[n-1]`
places C, CB and O, so Legacy and AMBER both depend on it.  Residue 0's state is half
invisible: `phi[0]` is inert, `psi[0]` is not.

That is a direct mechanism for energy-structure decoupling: any variance an objective puts
on residue n-1 is variance spent on a decision the metric is blind to, and an optimiser
will happily spend its budget there.  Nobody has measured it.

Measured exactly, with no sampling error, by the TOTAL SOBOL INDEX on the full 262,144-
configuration enumeration:

    T_i = 1 - Var_{s_-i}( E[ f | s_-i ] ) / Var(f)

which is the share of Var(f) that vanishes when variable i is averaged out -- main effect
plus every interaction involving i.  The enumeration is in odometer order so this is an
exact reshape-and-mean, not an estimate.

AMBER cannot be done this way (the cache holds 1.13% of the space), so it is measured by
direct construction: take M base configurations, sweep residue i through all k states, and
compare the within-sweep variance with the total.  That is the SAME quantity in
expectation and costs 28 ms per point.

    python -m s14.ener_invisible [--amber]
"""
from __future__ import annotations

import sys

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E

M_BASE = 120           # base configurations per target for the AMBER sweep


# ------------------------------------------------------------------ exact Sobol
def total_sobol(f, n, k):
    """Exact total Sobol index of every variable, for f indexed in odometer order."""
    f = np.asarray(f, float)
    v = f.var()
    if v <= 0:
        return np.full(n, np.nan)
    out = np.empty(n)
    for i in range(n):
        a = f.reshape(k ** i, k, k ** (n - 1 - i))
        cond = a.mean(axis=1)                      # E[f | s_-i]
        out[i] = 1.0 - cond.var() / v
    return out


def first_sobol(f, n, k):
    """Exact first-order (main-effect) Sobol index of every variable."""
    f = np.asarray(f, float)
    v = f.var()
    if v <= 0:
        return np.full(n, np.nan)
    out = np.empty(n)
    for i in range(n):
        a = f.reshape(k ** i, k, k ** (n - 1 - i))
        mu_i = a.mean(axis=(0, 2))                 # E[f | s_i]
        out[i] = mu_i.var() / v
    return out


def exact_target(pdb):
    z = E.enum(pdb)
    n, k = z.n, z.k
    rows = {}
    fields = {"rmsd": z.rmsd, "legacy": z.legacy, "prior": z.prior,
              "legacy_nosteric": z.obj("legacy_nosteric")}
    for t in E.LEG_TERMS:
        fields["leg_" + t] = z.leg[t]
    for nm, v in fields.items():
        T = total_sobol(v, n, k)
        S1 = first_sobol(v, n, k)
        rows[nm] = dict(total=T.tolist(), first=S1.tolist(),
                        T_last=float(T[n - 1]), T_first_res=float(T[0]),
                        T_interior_mean=float(np.nanmean(T[1:n - 1])),
                        # the operative number: share of the objective's variance that
                        # lives on the residue CA-RMSD cannot see, normalised by the
                        # average interior residue so chain-length effects cancel
                        T_last_over_interior=float(T[n - 1] / np.nanmean(T[1:n - 1]))
                        if np.nanmean(T[1:n - 1]) > 0 else float("nan"))
    return dict(pdb=pdb, n=n, k=k, per_field=rows)


# ------------------------------------------------------------------ AMBER by sweep
def amber_sweep(pdb, seed=0):
    """Variance of AMBER attributable to residue i, by explicit k-state sweeps.

    For a random base set, sweeping variable i through all k states and taking the
    within-sweep variance estimates E_s[ Var(f | s_-i) ] = T_i * Var(f) directly.
    """
    from s13 import qarch_lib as Q
    z = E.enum(pdb)
    sp = z.space()
    n, k = z.n, z.k
    rng = np.random.default_rng(seed)
    base = z.states(rng.choice(z.B, M_BASE, replace=False)).astype(int)
    # sweep the LAST residue, the FIRST residue, and one interior residue
    which = [n - 1, 0, n // 2]
    S = []
    for i in which:
        for b in range(M_BASE):
            for kk in range(k):
                s = base[b].copy(); s[i] = kk
                S.append(s)
    S = np.array(S, int)
    E.wait_for_memory(1.5, f"amber-invisible-{pdb}")
    print(f"  {pdb}: AMBER on {len(S)} configs (~{len(S)*0.028:.0f} s)", flush=True)
    vals = np.asarray(Q.amber_energies(sp, S, components=False, progress=500), float)
    vals = vals.reshape(len(which), M_BASE, k)
    # a monotone conditioning is mandatory: raw AMBER's variance is a delta spike on the
    # worst clash (project finding), so the raw variance share is an artefact
    soft = E.softcore(vals)
    tot_raw = vals.reshape(-1).var()
    tot_soft = soft.reshape(-1).var()
    out = {}
    for j, i in enumerate(which):
        out[str(i)] = dict(
            residue=int(i),
            within_var_raw=float(vals[j].var(axis=1).mean()),
            T_raw=float(vals[j].var(axis=1).mean() / tot_raw) if tot_raw > 0 else np.nan,
            within_var_soft=float(soft[j].var(axis=1).mean()),
            T_soft=float(soft[j].var(axis=1).mean() / tot_soft) if tot_soft > 0 else np.nan,
            # rank-scale: fraction of sweeps where the last-residue choice changes the
            # ARGMIN of the energy over that sweep
            frac_sweeps_nonconstant=float((soft[j].std(axis=1) > 1e-9).mean()))
    return dict(pdb=pdb, n=n, which=which, per_residue=out,
                note="residue n-1 is entirely invisible to CA-RMSD")


# ------------------------------------------------------------------ main
def main(with_amber=False):
    ex = [exact_target(p) for p in E.ENUM_TARGETS]
    print("=== EXACT TOTAL SOBOL INDEX (full enumeration, no sampling error) ===")
    print("T_i = share of the objective's variance that vanishes when residue i is "
          "averaged out.")
    print("Residue n-1 is INVISIBLE to CA-RMSD; residue 0 is half invisible "
          "(phi[0] inert).\n")
    print(f"{'field':22s} {'T[n-1]':>9s} {'T[0]':>9s} {'T interior':>11s} "
          f"{'T[n-1]/interior':>16s}")
    fields = list(ex[0]["per_field"])
    agg = {}
    for nm in fields:
        r = [x["per_field"][nm] for x in ex]
        agg[nm] = dict(
            T_last=float(np.nanmean([x["T_last"] for x in r])),
            T_first=float(np.nanmean([x["T_first_res"] for x in r])),
            T_int=float(np.nanmean([x["T_interior_mean"] for x in r])),
            ratio=float(np.nanmean([x["T_last_over_interior"] for x in r])))
        a = agg[nm]
        print(f"{nm:22s} {a['T_last']:9.5f} {a['T_first']:9.5f} {a['T_int']:11.5f} "
              f"{a['ratio']:16.3f}")
    # the headline: RMSD's T on the last residue must be EXACTLY zero
    tl = [x["per_field"]["rmsd"]["T_last"] for x in ex]
    print(f"\nCA-RMSD total Sobol index on residue n-1, per target: "
          f"{[f'{v:.3e}' for v in tl]}")
    print(f"CA-RMSD total Sobol index on residue 0, per target:   "
          f"{[f'{x['per_field']['rmsd']['first'][0]:.4f}' for x in ex]}")

    amb = None
    if with_amber:
        amb = [amber_sweep(p) for p in E.ENUM_TARGETS]
        print("\n=== AMBER, by explicit k-state sweep (28 ms/point) ===")
        print(f"{'residue role':22s} {'T_raw':>9s} {'T_softcore':>12s} "
              f"{'nonconstant sweeps':>19s}")
        for role, key in (("n-1 (INVISIBLE)", lambda r: str(r["n"] - 1)),
                          ("0 (half invisible)", lambda r: "0"),
                          ("interior n/2", lambda r: str(r["n"] // 2))):
            vs = [r["per_residue"][key(r)] for r in amb]
            print(f"{role:22s} {np.mean([v['T_raw'] for v in vs]):9.4f} "
                  f"{np.mean([v['T_soft'] for v in vs]):12.4f} "
                  f"{np.mean([v['frac_sweeps_nonconstant'] for v in vs]):19.3f}")

    E.write("ener_invisible", dict(
        what="variance of each objective carried by degrees of freedom CA-RMSD cannot see",
        exact=ex, exact_pooled=agg, amber_sweep=amb),
        n_expected=len(E.ENUM_TARGETS))
    return ex, agg, amb


if __name__ == "__main__":
    main(with_amber="--amber" in sys.argv)
