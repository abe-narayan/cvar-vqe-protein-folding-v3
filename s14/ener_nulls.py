"""SPRINT 14, ENER-8 -- the null control for the one attractive positive result.

`leg_torsion` was the ONLY objective in the E1 matrix whose selected-RMSD CI excluded zero
on the good side: 2.954 A against a random draw of 3.789, -0.835 A [-1.607, -0.105], 6W/3L.
The 1-local prior and `leg_torsion` were also the only two objectives with a large positive
excess AUC in the decoy test (+0.16 to +0.22).

The project has an explicit warning about exactly this shape of result: a zero-information
constant alpha-helix (phi -63, psi -42) beats uniform random sampling by 0.457 A, and the
1-local prior's apparent advantage is entirely a helix artefact with
rho(helix fraction, RMSD) = -0.744.  "Beats random" proves nothing here; the constant helix
is the baseline to clear.

So this module runs the controls:

    helix_frac        an objective that is PURELY the fraction of residues in the alpha
                      basin -- zero physics, zero sequence, zero energy
    constant_helix    the single configuration nearest (-63, -42) at every residue -- a
                      zero-information point estimate, not a ranking at all
    partialled        `leg_torsion` with helix fraction regressed out, to ask whether
                      anything survives
    per target        split by the target's own helicity, because the artefact should
                      vanish on non-helical targets

    python -m s14.ener_nulls
"""
from __future__ import annotations

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E
from s14.ener_matrix import ensemble


def helix_fraction(z, idx):
    S = z.states(np.asarray(idx)).astype(int)
    rows = np.arange(z.n)
    phi = np.degrees(z.PHI[rows[None, :], S]); psi = np.degrees(z.PSI[rows[None, :], S])
    return ((phi > -100) & (phi < -30) & (psi > -80) & (psi < -5)).mean(1)


def constant_helix_index(z, phi0=-63.0, psi0=-42.0):
    """ORACLE-FREE: the configuration nearest a constant alpha-helix at every residue."""
    d = (np.abs(_wrapdeg(np.degrees(z.PHI) - phi0))
         + np.abs(_wrapdeg(np.degrees(z.PSI) - psi0)))
    return int(z.index_of(np.argmin(d, axis=1)[None])[0])


def _wrapdeg(a):
    return (np.asarray(a, float) + 180.0) % 360.0 - 180.0


def _partial_out(x, c):
    x = E.rank_norm(x) - E.rank_norm(x).mean()
    c = E.rank_norm(c) - E.rank_norm(c).mean()
    return x - (x @ c) / (c @ c) * c if c @ c > 0 else x


def main():
    rows = []
    for pdb in E.ENUM_TARGETS:
        d = ensemble(pdb)
        z = d["_z"]
        R = d["rmsd"]
        hf = helix_fraction(z, d["idx"])
        ch = constant_helix_index(z)
        # native helicity of the target itself (ORACLE, for the split only)
        nat_hf = float(helix_fraction(z, np.array([z.snap_index]))[0])
        rows.append(dict(
            pdb=pdb, fold=int(z.fold), pool_mean=float(R.mean()), pool_best=float(R.min()),
            native_helix_frac=nat_hf,
            leg_torsion_sel=E.argmin_rmsd(d["leg_torsion"], R),
            prior_sel=E.argmin_rmsd(d["prior"], R),
            legacy_sel=E.argmin_rmsd(d["legacy"], R),
            amber_sel=E.argmin_rmsd(d["amber"], R),
            # NULL 1: pure helicity, no physics at all
            helix_sel=E.argmin_rmsd(-hf, R),
            # NULL 2: the zero-information constant alpha-helix point estimate
            const_helix_rmsd=float(z.rmsd[ch]),
            # NULL 3: leg_torsion with helicity removed
            torsion_partialled_sel=E.argmin_rmsd(_partial_out(d["leg_torsion"], hf), R),
            prior_partialled_sel=E.argmin_rmsd(_partial_out(d["prior"], hf), R),
            rho_torsion_helix=E.spearman(d["leg_torsion"], hf),
            rho_prior_helix=E.spearman(d["prior"], hf),
            rho_helix_rmsd=E.spearman(-hf, R),
            rho_legacy_helix=E.spearman(d["legacy"], hf),
            rho_amber_helix=E.spearman(d["amber"], hf)))
        print(f"{pdb}: nat_helix {nat_hf:.2f}  leg_torsion {rows[-1]['leg_torsion_sel']:.3f}"
              f"  pure-helix {rows[-1]['helix_sel']:.3f}"
              f"  const-helix {rows[-1]['const_helix_rmsd']:.3f}"
              f"  torsion|helix {rows[-1]['torsion_partialled_sel']:.3f}"
              f"  pool {R.mean():.3f}", flush=True)

    folds = [r["fold"] for r in rows]
    base = np.array([r["pool_mean"] for r in rows])
    print("\n=== NULL CONTROLS FOR THE ONE POSITIVE RESULT ===")
    print(f"{'arm':32s} {'sel':>7s} {'vs random':>10s} {'CI':>20s} {'W/L':>6s}")
    for nm, key in (("leg_torsion (the claim)", "leg_torsion_sel"),
                    ("NULL: pure helix fraction", "helix_sel"),
                    ("NULL: constant alpha-helix point", "const_helix_rmsd"),
                    ("leg_torsion | helix partialled", "torsion_partialled_sel"),
                    ("1-local prior", "prior_sel"),
                    ("prior | helix partialled", "prior_partialled_sel"),
                    ("Legacy total", "legacy_sel"),
                    ("AMBER total", "amber_sel")):
        v = np.array([r[key] for r in rows])
        p = I.paired(v, base, folds=folds)
        print(f"{nm:32s} {v.mean():7.3f} {p['mean_diff']:+10.3f} "
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}")
    # against the constant-helix baseline, which is the real baseline to clear
    ch = np.array([r["const_helix_rmsd"] for r in rows])
    print(f"\n{'arm':32s} {'vs CONSTANT-HELIX baseline':>28s}")
    for nm, key in (("leg_torsion", "leg_torsion_sel"),
                    ("pure helix fraction", "helix_sel"),
                    ("1-local prior", "prior_sel"),
                    ("Legacy total", "legacy_sel"),
                    ("AMBER total", "amber_sel")):
        v = np.array([r[key] for r in rows])
        p = I.paired(v, ch, folds=folds)
        print(f"{nm:32s} {p['mean_diff']:+10.3f} "
              f"[{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] {p['n_better']}/{p['n_worse']}")
    print("\ncorrelations with the helix fraction of the configuration:")
    for k in ("rho_torsion_helix", "rho_prior_helix", "rho_legacy_helix",
              "rho_amber_helix", "rho_helix_rmsd"):
        print(f"  {k:24s} {np.mean([r[k] for r in rows]):+.3f}")
    # split by the target's own helicity
    hi = [r for r in rows if r["native_helix_frac"] >= 0.5]
    lo = [r for r in rows if r["native_helix_frac"] < 0.5]
    print(f"\nsplit by the TARGET's native helicity "
          f"(helical n={len(hi)}, non-helical n={len(lo)}):")
    for nm, key in (("leg_torsion", "leg_torsion_sel"), ("pure helix", "helix_sel"),
                    ("prior", "prior_sel"), ("Legacy", "legacy_sel")):
        f = lambda g, k: (np.mean([r[k] for r in g]) - np.mean([r["pool_mean"] for r in g])
                          if g else float("nan"))
        print(f"  {nm:16s} helical {f(hi, key):+.3f}   non-helical {f(lo, key):+.3f}")
    E.write("ener_nulls", dict(
        what="null controls for leg_torsion and the 1-local prior",
        per_target=rows), n_expected=len(E.ENUM_TARGETS))
    return rows


if __name__ == "__main__":
    main()
