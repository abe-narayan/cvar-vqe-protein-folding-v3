"""SPRINT 15, ALIGN, TASK 1 -- CHARACTERISE THE LOW-RESPONSE TORSION SUBSPACE.

The INFO workstream measured that real torsion errors sit in the RMSD-quiet directions
(alignment 0.565 for the best fragment against a 0.945 random-direction null) but never
described the subspace itself.  This module does, on all 126 tuning targets:

  1. the singular-value spectrum of J = d(superposed CA)/d(phi, psi), normalised, with the
     participation ratio, the 90%-damage rank, and the dimension of the near-null space;
  2. where the four known-inert terminal torsions sit in that spectrum;
  3. how much ANGULAR ERROR BUDGET the quiet half of the spectrum absorbs at a fixed RMSD
     cost, relative to a random direction -- the "how much can be absorbed there" number;
  4. the decomposition of each real channel's measured alignment into WHICH directions its
     error occupies, i.e. the 0.565 taken apart;
  5. **the native-free question**: is the quiet subspace identifiable WITHOUT the native?
     Principal angles between the bottom-half subspace of J at the native torsions and the
     bottom-half subspace of J at the incumbent's emitted torsions.  If those subspaces do
     not overlap, every intervention in `s15/align_fit.py` is dead before it starts, and
     that is the cheapest possible way to learn it.

`J` at native torsions is ORACLE.  `J` at emitted torsions is NATIVE-FREE.  Both are computed
and both are labelled.

    python -m s15.align_jac
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy.linalg as la                    # noqa: E402
from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import info_lib as L                # noqa: E402
from s15.info_regime import channels         # noqa: E402
import peptide_db as pdb                     # noqa: E402
from s15 import seed as SD                   # noqa: E402

CH = ("retrieval top-75 circular mean", "top-75 medoid window",
      "constant alpha-helix (control)", "incumbent projected torsions",
      "ORACLE best pool window")


def principal_angles(A1, A2):
    """cos of the principal angles between two orthonormal column spaces."""
    return la.svd(A1.T @ A2, compute_uv=False)


def main():
    rows = []
    for t in I.targets():
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        n = len(nphi)
        J, _, _ = A.sup_jacobian(nphi, npsi)
        s, V = A.spectrum(J)
        s2 = s * s
        tot = s2.sum()
        m = 2 * n
        cum = np.cumsum(s2) / tot
        rank90 = int(np.searchsorted(cum, 0.90) + 1)
        pr = float(tot ** 2 / (s2 ** 2).sum())               # participation ratio
        near_null = int((s <= 1e-8 * s[0]).sum())
        quiet = int((s <= 0.10 * s[0]).sum())

        #: budget absorbed by the QUIET HALF vs a random direction.  An error placed
        #: isotropically in a set S of singular directions does RMS damage
        #: sqrt(mean_{i in S} s_i^2) per unit angle, so the tolerated angle scales inversely.
        h = m // 2
        lo = s[h:]                                           # the quiet half
        hi = s[:h]
        rms_all = float(np.sqrt(s2.mean()))
        rms_lo = float(np.sqrt((lo ** 2).mean()))
        rms_hi = float(np.sqrt((hi ** 2).mean()))

        rec = {"pdb": p, "n": n, "m": m, "s": s.tolist(), "s_max": float(s[0]),
               "rank90": rank90, "participation_ratio": pr, "near_null": near_null,
               "quiet_10pct": quiet,
               "budget_ratio_quiet_half": rms_all / max(rms_lo, 1e-30),
               "budget_ratio_loud_half": rms_all / max(rms_hi, 1e-30),
               "s_lo_over_s_hi": rms_lo / max(rms_hi, 1e-30)}

        #: the four inert torsions, as coordinate directions, against the spectrum
        inert = [0, n - 1, n, 2 * n - 1]                     # phi0, phi_{n-1}, psi0, psi_{n-1}
        cn = la.norm(J, axis=0)
        rec["inert_col_norms"] = [float(cn[k]) for k in inert]
        rec["max_col_norm"] = float(cn.max())
        rec["median_interior_col_norm"] = float(np.median(np.delete(cn, inert)))
        #: how much of each inert coordinate direction lies in the exact null space
        Vn = V[:, s <= 1e-8 * s[0]]
        rec["inert_in_nullspace"] = [float((Vn[k] ** 2).sum()) for k in inert] if \
            Vn.shape[1] else [0.0] * 4

        #: per-channel error decomposition
        rec["ch"] = {}
        chs = channels(t)
        for nm in CH:
            phi, psi = chs[nm]
            e = np.concatenate([A.wrap(np.asarray(phi, float) - nphi),
                                A.wrap(np.asarray(psi, float) - npsi)])
            if la.norm(e) < 1e-12:
                continue
            c = V.T @ e
            q = c * c
            q = q / q.sum()
            dmg = s2 * q
            dmg = dmg / dmg.sum()
            rec["ch"][nm] = {
                "align": A.alignment(J, e),
                "err_frac_quiet_half": float(q[h:].sum()),
                "dmg_frac_quiet_half": float(dmg[h:].sum()),
                "err_frac_bottom_quarter": float(q[3 * m // 4:].sum()),
                "err_frac_top_quarter": float(q[:m // 4].sum()),
                "err_frac_nullspace": float(q[s <= 1e-8 * s[0]].sum()),
                "rms_deg": float(np.degrees(np.sqrt((e ** 2).mean()))),
            }

        #: NATIVE-FREE identifiability of the subspace
        iphi, ipsi = chs["incumbent projected torsions"]
        Ji, _, _ = A.sup_jacobian(np.asarray(iphi, float), np.asarray(ipsi, float))
        si, Vi = A.spectrum(Ji)
        for frac, tag in ((0.5, "half"), (0.25, "quarter")):
            k = max(1, int(round(frac * m)))
            ca = principal_angles(V[:, m - k:], Vi[:, m - k:])
            rec[f"subspace_overlap_{tag}"] = float((ca ** 2).mean())
            rec[f"subspace_min_cos_{tag}"] = float(ca.min())
        #: and a random-rotation null for that overlap
        # BRIEF S3.7 fix (Sprint 16, REPAIR): was `np.random.default_rng(abs(hash(p)) % 7 + 5)`.
        # `hash()` is salted per process, so this null was not reproducible across
        # interpreters.  The reported value is unchanged -- s16/audit_FINDINGS.md S2.1 shows
        # the quantity is analytically k/m = 0.500 and both seedings return it -- and
        # `s16/repair_seedfix.py` re-measures it under both seedings to confirm.
        rng = SD.stable_rng("align_jac", "subspace_null", p)   # only a null; value irrelevant
        Q, _ = la.qr(rng.normal(size=(m, m)))
        k = m // 2
        rec["subspace_overlap_null"] = float((principal_angles(V[:, m - k:],
                                                              Q[:, :k]) ** 2).mean())
        #: spectral agreement of the two Jacobians
        rec["spec_corr_native_vs_emitted"] = L.pearson(s, si)
        rows.append(rec)

    # ------------------------------------------------------------------ report
    out = {"rows": rows}
    ns = np.asarray([r["n"] for r in rows])
    print("TASK 1 -- THE LOW-RESPONSE TORSION SUBSPACE, 126 targets, J at NATIVE torsions "
          "(ORACLE DIAGNOSTIC)\n")

    #: normalised spectrum on a common axis (fraction of the way through the spectrum)
    grid = np.linspace(0, 1, 21)
    curves = []
    for r in rows:
        s = np.asarray(r["s"], float)
        curves.append(np.interp(grid, np.linspace(0, 1, len(s)), s / s[0]))
    C = np.asarray(curves)
    print("normalised singular value s_i/s_1 against fractional index i/2n "
          "(mean over targets):")
    print("  i/2n : " + " ".join(f"{g:5.2f}" for g in grid[::2]))
    print("  s/s1 : " + " ".join(f"{v:5.3f}" for v in C.mean(0)[::2]))
    out["spectrum_curve"] = {"grid": grid.tolist(), "mean": C.mean(0).tolist(),
                             "sd": C.std(0).tolist()}

    def col(k, f=float):
        return np.asarray([f(r[k]) for r in rows])

    for k, lab in (("rank90", "directions carrying 90% of trace(J'J)"),
                   ("participation_ratio", "participation ratio (effective rank)"),
                   ("near_null", "exactly-null directions"),
                   ("quiet_10pct", "directions with s_i <= 0.10 s_1")):
        v = col(k)
        print(f"\n{lab:<48} mean {v.mean():7.3f}  median {np.median(v):7.3f}"
              f"  (2n mean {2 * ns.mean():.1f})")
        out[k] = {"mean": float(v.mean()), "median": float(np.median(v))}
    v = col("rank90") / (2 * ns)
    print(f"{'  ... as a fraction of 2n':<48} mean {v.mean():7.3f}")
    out["rank90_frac"] = float(v.mean())

    br = col("budget_ratio_quiet_half")
    lr = col("s_lo_over_s_hi")
    print(f"\nANGULAR BUDGET: an error placed isotropically in the QUIET HALF of the spectrum")
    print(f"  tolerates {br.mean():.2f}x (median {np.median(br):.2f}x) more RMS angle than a "
          f"random direction at the same RMSD cost")
    print(f"  the quiet half's RMS singular value is {lr.mean():.3f}x the loud half's")
    out["budget_ratio_quiet_half"] = {"mean": float(br.mean()),
                                      "median": float(np.median(br))}
    out["s_lo_over_s_hi"] = {"mean": float(lr.mean())}

    inert = np.asarray([r["inert_col_norms"] for r in rows])
    mx = col("max_col_norm")
    med = col("median_interior_col_norm")
    print(f"\nTHE FOUR INERT TORSIONS in the superposed Jacobian (column norms):")
    for k, nm in enumerate(("phi[0]", "phi[n-1]", "psi[0]", "psi[n-1]")):
        print(f"  {nm:<10} max over 126 targets {inert[:, k].max():.3e}")
    print(f"  interior median column norm {med.mean():.3f}   max column norm {mx.mean():.3f}"
          f"   ratio {inert.max() / med.mean():.1e}")
    out["inert_col_norm_max"] = float(inert.max())
    out["interior_median_col_norm"] = float(med.mean())

    print("\nPER-CHANNEL ERROR DECOMPOSITION (share of ||e||^2 by spectral region)\n")
    print(f"{'channel':<34}{'align':>8}{'top 1/4':>9}{'quiet 1/2':>11}"
          f"{'bottom 1/4':>12}{'null':>8}{'dmg quiet':>11}{'RMS deg':>9}")
    out["channels"] = {}
    for nm in CH:
        vv = [r["ch"][nm] for r in rows if nm in r["ch"]]
        if not vv:
            continue
        g = {k: float(np.mean([x[k] for x in vv])) for k in vv[0]}
        g["n"] = len(vv)
        out["channels"][nm] = g
        print(f"{nm:<34}{g['align']:>8.3f}{g['err_frac_top_quarter']:>9.3f}"
              f"{g['err_frac_quiet_half']:>11.3f}{g['err_frac_bottom_quarter']:>12.3f}"
              f"{g['err_frac_nullspace']:>8.3f}{g['dmg_frac_quiet_half']:>11.3f}"
              f"{g['rms_deg']:>9.1f}")
    print("\n  an ISOTROPIC error puts 0.250 / 0.500 / 0.250 in top-quarter / quiet-half / "
          "bottom-quarter\n  by construction, so those columns read directly against "
          "0.250 / 0.500 / 0.250.")

    print("\nNATIVE-FREE IDENTIFIABILITY of the quiet subspace")
    for tag in ("half", "quarter"):
        v = col(f"subspace_overlap_{tag}")
        mn = col(f"subspace_min_cos_{tag}")
        print(f"  bottom-{tag} subspace, J(native) vs J(incumbent emitted): "
              f"mean cos^2 {v.mean():.3f}  min cos {mn.mean():.3f}")
        out[f"subspace_overlap_{tag}"] = float(v.mean())
    nul = col("subspace_overlap_null")
    print(f"  random-subspace null (same dimensions): mean cos^2 {nul.mean():.3f}")
    out["subspace_overlap_null"] = float(nul.mean())
    sc = col("spec_corr_native_vs_emitted")
    print(f"  correlation of the two spectra: {sc.mean():+.3f}")
    out["spec_corr_native_vs_emitted"] = float(sc.mean())

    L.jwrite("align_jac", out)
    return out


if __name__ == "__main__":
    main()
