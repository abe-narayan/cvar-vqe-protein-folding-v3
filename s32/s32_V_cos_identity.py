"""S32 LANE V -- is lane R's `cos_align` a measurement, or an algebraic restatement?

`s32/s32_R_ladder_analyse.py:80`:

    cos_align = (e**2 + d**2 - chain**2) / (2*e*d)

That is the law of cosines INVERTED on three numbers the row already carries:
`e = RMSD(cloud, native)`, `d = RMSD(chain, cloud)`, `chain = RMSD(chain, native)`.

TWO SEPARATE QUESTIONS, and they have different answers:

 1. **Is `cos` independent evidence?**  No, by construction: it is a bijection with `price =
    chain - e` given (e, d).  Any sentence that quotes the price AND the cos as two pieces of
    support is double-counting one measurement.  Stated here and asserted numerically.

 2. **Is `cos` even the cosine it is named for?**  Each of e, d, chain is a *Kabsch* RMSD --
    a distance in Kendall SHAPE space, which is positively curved, with each pair superposed
    INDEPENDENTLY.  The law of cosines is a Euclidean identity.  So the inverted quantity is a
    triangle-defect statistic, not an angle, unless the three objects can be placed in ONE
    common frame without loss.  This measures the true cosine directly, in a common frame, and
    compares:

        superpose native onto the cloud  -> t_C ;  e_vec = t_C - C
        superpose chain  onto the cloud  -> X_C ;  d_vec = X_C - C
        cos_direct = <d_vec, e_vec> / (|d_vec| |e_vec|)
        chain_euclid = |X_C - t_C| / sqrt(n)          vs the Kabsch chain RMSD

    If `chain_euclid != chain_kabsch`, the triangle is not Euclidean and `cos_align` is biased.
    If `cos_direct != cos_align`, the number does not mean what its name says.

ORACLE / NOT DEPLOYABLE: the native is read.
"""
from __future__ import annotations
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
R_STRUCTS = os.path.join(RESULTS, "s32_R_structs")
S29_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(R_STRUCTS, "*.npz"))):
        pdb = os.path.basename(f)[:-4]
        with np.load(f, allow_pickle=True) as z:
            if "prod_ca" not in z.files:
                continue
            X = np.asarray(z["prod_ca"], float)
        cf = os.path.join(S29_STRUCTS, "%s.npz" % pdb)
        if not os.path.exists(cf):
            continue
        with np.load(cf) as z:
            C = np.asarray(z["prod"], float)
        nat = I.load_univ(pdb)["nat_ca"]
        n = len(C)

        e = I.ca_rmsd(C, nat)                    # cloud error       (Kabsch)
        d = I.ca_rmsd(X, C)                      # displacement      (Kabsch)
        ch = I.ca_rmsd(X, nat)                   # chain error       (Kabsch)
        den = 2.0 * e * d
        cos_alg = (e ** 2 + d ** 2 - ch ** 2) / den if den > 1e-12 else np.nan

        # --- the direct measurement, all three objects in ONE frame (the cloud's)
        tC = I.superpose_batch(nat[None], C)[0]
        XC = I.superpose_batch(X[None], C)[0]
        ev = (tC - C).ravel(); dv = (XC - C).ravel()
        ne, nd = np.linalg.norm(ev), np.linalg.norm(dv)
        cos_dir = float(ev @ dv / (ne * nd)) if ne > 1e-12 and nd > 1e-12 else np.nan
        chain_euclid = float(np.linalg.norm(XC - tC) / np.sqrt(n))

        rows.append(dict(pdb=pdb, n=n, e=e, d=d, chain=ch,
                         e_common=float(ne / np.sqrt(n)), d_common=float(nd / np.sqrt(n)),
                         cos_algebraic=float(cos_alg), cos_direct=cos_dir,
                         chain_kabsch=ch, chain_euclid=chain_euclid,
                         price=float(ch - e)))
    if not rows:
        print("no lane-R structs with prod_ca yet")
        return

    G = lambda k: np.array([r[k] for r in rows], float)
    ca, cd = G("cos_algebraic"), G("cos_direct")
    ck, ce = G("chain_kabsch"), G("chain_euclid")
    out = dict(n=len(rows), oracle="ORACLE / NOT DEPLOYABLE", per_target=rows)

    print("=" * 100)
    print("S32-V: lane R's cos_align -- identity check and direct measurement.  n=%d targets"
          % len(rows))
    print()
    print("1. IS IT INDEPENDENT EVIDENCE?")
    #: reconstruct the price from (e, d, cos) alone: if it returns `price` exactly, cos carries
    #: nothing the price does not already carry.
    back = np.sqrt(np.maximum(G("e") ** 2 + G("d") ** 2 - 2 * G("e") * G("d") * ca, 0)) - G("e")
    err = float(np.abs(back - G("price")).max())
    out["price_reconstructed_from_cos_max_abs_err"] = err
    print("   price rebuilt from (e, d, cos_align) alone: max |error| = %.3e over %d targets"
          % (err, len(rows)))
    print("   -> cos_align is a BIJECTION with the price given (e, d).  It is a re-parameterisation,")
    print("      not a second measurement; quoting both as support is double-counting one number.")
    print()
    print("2. IS IT THE COSINE IT IS NAMED FOR?")
    print("   cos_algebraic  mean %+.4f  median %+.4f  min %+.4f  max %+.4f  (|cos|>1 on %d)"
          % (np.nanmean(ca), np.nanmedian(ca), np.nanmin(ca), np.nanmax(ca), int((np.abs(ca) > 1).sum())))
    print("   cos_direct     mean %+.4f  median %+.4f  min %+.4f  max %+.4f"
          % (np.nanmean(cd), np.nanmedian(cd), np.nanmin(cd), np.nanmax(cd)))
    dif = ca - cd
    print("   difference     mean %+.4f  |mean| %.4f  p90 %.4f  max %.4f"
          % (np.nanmean(dif), np.nanmean(np.abs(dif)), np.nanpercentile(np.abs(dif), 90),
             np.nanmax(np.abs(dif))))
    print("   chain, Kabsch vs Euclidean-in-one-frame: mean %.4f vs %.4f, max |diff| %.4f"
          % (ck.mean(), ce.mean(), np.abs(ck - ce).max()))
    out["cos_algebraic"] = dict(mean=float(np.nanmean(ca)), median=float(np.nanmedian(ca)),
                                n_abs_gt_1=int((np.abs(ca) > 1).sum()))
    out["cos_direct"] = dict(mean=float(np.nanmean(cd)), median=float(np.nanmedian(cd)))
    out["cos_difference"] = dict(mean=float(np.nanmean(dif)), abs_mean=float(np.nanmean(np.abs(dif))),
                                 p90=float(np.nanpercentile(np.abs(dif), 90)),
                                 max=float(np.nanmax(np.abs(dif))))
    out["chain_kabsch_vs_euclid"] = dict(kabsch=float(ck.mean()), euclid=float(ce.mean()),
                                         max_abs_diff=float(np.abs(ck - ce).max()))
    agree = np.nanmean(np.abs(dif)) < 0.02
    out["verdict"] = (
        "cos_align is an exact re-parameterisation of the price and agrees with the directly "
        "measured cosine to %.4f -- it is not independent evidence, but it does mean what it says"
        % np.nanmean(np.abs(dif)) if agree else
        "cos_align is an exact re-parameterisation of the price AND differs from the directly "
        "measured cosine by %.4f on average -- the Kabsch triangle is not Euclidean, so the "
        "quantity is a triangle defect, not an angle" % np.nanmean(np.abs(dif)))
    print()
    print("   VERDICT: %s" % out["verdict"])
    print("=" * 100)
    with open(os.path.join(RESULTS, "s32_V_cos_identity.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out


if __name__ == "__main__":
    main()
