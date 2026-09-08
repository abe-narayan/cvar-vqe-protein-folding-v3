"""s16/verify_csteer.py -- VERIFY workstream, hostile re-reading of the PIVOT (s16/csteer.py).

  P1  IS `3.647 * (1 - s)` A MEASUREMENT OR A THEOREM?
      Claim (CSTEER_FINDINGS section 1): "Coordinate space is exactly linear, to three
      decimals at every rung, and the identity survives Kabsch superposition."

      THEOREM.  Let N be optimally superposed onto S, so the identity rotation is Kabsch-
      optimal for the pair; equivalently the cross-covariance H = Nc^T Sc is symmetric PSD
      and the centroids coincide.  For X_s = S + s (N - S) = (1-s) S + s N,

          Xc^T Nc = (1-s) Sc^T Nc + s Nc^T Nc = (1-s) H^T + s Nc^T Nc,

      a sum of two symmetric PSD matrices, hence symmetric PSD; so the identity is Kabsch-
      optimal for (X_s, N) as well and

          rmsd(X_s, N) = ||Xc - Nc|| / sqrt(n) = (1-s) ||Sc - Nc|| / sqrt(n)
                       = (1 - s) rmsd(S, N)                              EXACTLY, all s in [0,1].

      No protein, no pipeline, no fit enters.  This module demonstrates it on random point
      clouds that have nothing to do with the instrument.  If it holds there, section 1 is a
      tautology published as an exhibit.

  P2  THE ZERO-INFORMATION REFERENCE the AUDIT workstream made mandatory.
      CSTEER_FINDINGS section 2 claims that from the best start (`avg`, 3.048 A) every
      native-free direction has a cosine with the true residual indistinguishable from zero,
      with a random control at 0.005.  A matched RANDOM control is not a zero-information
      REFERENCE: random has no structure at all, whereas a constant ideal alpha-helix, an
      ideal beta-strand and an arbitrary retrieval window all carry real peptide structure
      and no target-specific information.  If one of THOSE has a cosine the sprint's channels
      do not, the section's framing is wrong.

      Also checked: the coordinator scales the random direction to another arm's norm.  Cosine
      is scale-invariant, so that matching cannot affect the cos column at all -- it only
      affects the LFO ladder.  Stated explicitly because the writeup presents it as the
      control's matching.

  P3  A ZERO-INFORMATION REFERENCE FOR THE LADDER ITSELF: how far apart are two arbitrary
      conformations of the same peptide?  This prices what "3.647 A" and "3.048 A" mean.

ORACLE.  The native is read only to form the residual `r` for the cosine, which is the
quantity under audit, and to score.  Every DIRECTION here is native-free or information-free.

Writes `s16/results/verify_csteer.json`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

HELIX = (np.deg2rad(-57.0), np.deg2rad(-47.0))
STRAND = (np.deg2rad(-139.0), np.deg2rad(135.0))
PPII = (np.deg2rad(-75.0), np.deg2rad(145.0))


def _sup(X, T):
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(T, float))[0]


def _cos(d, r):
    nd = float(np.linalg.norm(d)); nr = float(np.linalg.norm(r))
    if nd < 1e-12 or nr < 1e-12:
        return 0.0
    return float((d.ravel() @ r.ravel()) / (nd * nr))


# -------------------------------------------------------------------------- P1
def p1(out):
    print("=" * 96)
    print("P1  IS COORDINATE-SPACE LINEARITY A MEASUREMENT OR A THEOREM?")
    print("=" * 96)
    rng = SD.stable_rng("verify", "linearity")
    steps = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0])
    worst = 0.0
    for _ in range(400):
        n = int(rng.integers(5, 40))
        # two point clouds with NOTHING to do with peptides
        S = rng.standard_normal((n, 3)) * rng.uniform(0.5, 8.0)
        N = rng.standard_normal((n, 3)) * rng.uniform(0.5, 8.0) + rng.standard_normal(3) * 20
        Nsup = _sup(N, S)
        r = Nsup - S
        base = I.ca_rmsd(S, N)
        for s in steps:
            got = I.ca_rmsd(S + s * r, N)
            worst = max(worst, abs(got - (1.0 - s) * base))
    print(f"  400 random point-cloud pairs (n = 5..40), 9 rungs, NO protein involved.")
    print(f"  max | rmsd(S + s*r, N)  -  (1-s) * rmsd(S, N) |  =  {worst:.3e} A")
    print("  VERDICT: the linearity is an algebraic IDENTITY, not a property of this pipeline.")
    print("  It holds for arbitrary point clouds because a convex combination of two symmetric")
    print("  PSD cross-covariances is symmetric PSD, so the Kabsch rotation never moves off the")
    print("  identity along the segment.  `3.647*(1-s) to three decimals at every rung` is a")
    print("  check that floating point works.  See also `law_pred = base*sqrt(1-c^2)` in")
    print("  csteer.py, which is the same algebra written out in the module docstring.")
    out["p1"] = {"max_abs_deviation_A": float(worst), "n_trials": 400}


# -------------------------------------------------------------------------- P2 / P3
def p2(out):
    from s14 import avgspace as AV
    from s14 import retprior as R
    from s15 import distcal as C
    print("\n" + "=" * 96)
    print("P2  ZERO-INFORMATION REFERENCES FOR THE COSINES, from the best start (`avg`)")
    print("=" * 96)
    tg = I.targets()
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)

    names = ["med", "proj", "helix", "strand", "ppii", "poolwin", "poolwin_worst",
             "rand", "randrot"]
    cos = {k: [] for k in names}
    base_avg, base_helix, pool_pair, pool_to_avg = [], [], [], []

    for t in tg:
        p = t["pdb"]; d = data[p]
        n = int(d["n"]); seq = d["seq"]; fold = int(d["fold"])
        nat = np.asarray(d["nat"], float)
        W = np.asarray(AV.top75_windows(p)[0], float)
        P = I.pairwise_rmsd(W)
        avg_c, b_med = I.coordinate_average(W, P)
        avg = np.asarray(avg_c, float)
        med = np.asarray(W[b_med], float)
        proj = np.asarray(I.project(avg, seq, fold)["ca"], float)

        S = avg
        Nsup = _sup(nat, S)
        r = Nsup - S
        base_avg.append(float(I.ca_rmsd(S, nat)))

        # --- zero-information structures: constant ideal secondary structure
        hel = I.build_ca(np.full(n, HELIX[0]), np.full(n, HELIX[1]))
        strd = I.build_ca(np.full(n, STRAND[0]), np.full(n, STRAND[1]))
        pp2 = I.build_ca(np.full(n, PPII[0]), np.full(n, PPII[1]))
        base_helix.append(float(I.ca_rmsd(hel, nat)))

        rng = SD.stable_rng(p, "verify_csteer")
        k = int(rng.integers(0, len(W)))
        win = np.asarray(W[k], float)
        kw = int(np.argmax(P[b_med]))            # the window FURTHEST from the medoid
        winw = np.asarray(W[kw], float)

        dirs = {"med": med, "proj": proj, "helix": hel, "strand": strd, "ppii": pp2,
                "poolwin": win, "poolwin_worst": winw}
        for nm, T in dirs.items():
            cos[nm].append(_cos(_sup(T, S) - S, r))
        z = rng.standard_normal(S.shape); z -= z.mean(0)
        cos["rand"].append(_cos(z, r))
        # a random direction confined to the same 3n-6 subspace the real ones live in:
        # remove the three infinitesimal-rotation modes as well as translation.
        Sc = S - S.mean(0)
        G = np.stack([np.cross(np.array(e), Sc) for e in np.eye(3)])   # (3, n, 3)
        Z = z.copy()
        for g in G:
            gg = g / max(np.linalg.norm(g), 1e-12)
            Z -= (Z.ravel() @ gg.ravel()) * gg
        cos["randrot"].append(_cos(Z, r))

        pool_pair.append(float(P[np.triu_indices(len(W), 1)].mean()))
        pool_to_avg.append(float(np.mean([I.ca_rmsd(w, avg) for w in W])))

    print(f"  n = {len(tg)}, start = `avg` (top-75 coordinate average), mean "
          f"{np.mean(base_avg):.3f} A")
    print(f"  {'direction':<16}{'kind':<22}{'cos [95% CI fold-clustered]':>34}{'|cos| mean':>12}")
    kinds = {"med": "native-free channel", "proj": "native-free channel",
             "helix": "ZERO-INFORMATION", "strand": "ZERO-INFORMATION",
             "ppii": "ZERO-INFORMATION", "poolwin": "ZERO-INFO (arbitrary win)",
             "poolwin_worst": "native-free (outlier)",
             "rand": "random control", "randrot": "random, rot-modes out"}
    res = {}
    for nm in names:
        v = np.asarray(cos[nm], float)
        ci = I.paired(v, np.zeros(len(v)), folds=folds, names=pdbs)
        res[nm] = {"mean": float(v.mean()), "ci95": ci["ci95"],
                   "median": float(np.median(v)), "abs_mean": float(np.abs(v).mean()),
                   "wl": [int((v > 0).sum()), int((v < 0).sum())]}
        print(f"  {nm:<16}{kinds[nm]:<22}{v.mean():>+14.3f} "
              f"[{ci['ci95'][0]:+.3f},{ci['ci95'][1]:+.3f}]{np.abs(v).mean():>12.3f}")
    out["p2_cos_from_avg"] = res

    # the coordinator's own numbers, re-read from the artefact with the target as the unit
    print("\n  the artefact's own cosines from `avg`, re-read with fold-clustered CIs:")
    cj = json.load(open(os.path.join(RESULTS, "csteer.json")))
    crows = cj["rows"]
    cf = np.asarray([r["fold"] for r in crows], int)
    cp = [r["pdb"] for r in crows]
    rep = {}
    for dn in ("fit", "proj", "med", "ens", "rand0", "rand1", "rand2"):
        v = np.asarray([r["cells"]["avg"]["arms"][dn]["cos"] for r in crows], float)
        ci = I.paired(v, np.zeros(len(v)), folds=cf, names=cp)
        rep[dn] = {"mean": float(v.mean()), "ci95": ci["ci95"]}
        print(f"    {dn:<8}{v.mean():>+9.3f} [{ci['ci95'][0]:+.3f},{ci['ci95'][1]:+.3f}]")
    out["p2_artefact_reread"] = rep

    print("\n  NOTE ON THE RANDOM CONTROL'S MATCHING.  csteer.py scales the random direction to")
    print("  another arm's norm.  Cosine is scale-invariant, so the scaling has NO effect on the")
    print("  cos column that section 2 quotes; it matters only for the LFO ladder.  Verified:")
    v0 = np.asarray([r["cells"]["avg"]["arms"]["rand0"]["cos"] for r in crows], float)
    print(f"    rand0 cos mean {v0.mean():+.4f}; this module's unscaled Gaussian "
          f"{np.mean(cos['rand']):+.4f}; rot-modes removed {np.mean(cos['randrot']):+.4f}")
    print("    -- all three are the same null, as they must be.")

    print("\n" + "=" * 96)
    print("P3  WHAT DO THESE RMSDs MEAN?  zero-information reference points")
    print("=" * 96)
    print(f"  mean CA-RMSD, `avg` to native                          {np.mean(base_avg):.3f} A")
    print(f"  mean CA-RMSD, constant ideal alpha-helix to native     {np.mean(base_helix):.3f} A")
    print(f"  mean pairwise CA-RMSD among the top-75 windows         {np.mean(pool_pair):.3f} A")
    print(f"  mean CA-RMSD, a top-75 window to their own average     {np.mean(pool_to_avg):.3f} A")
    out["p3"] = {"avg_to_native": float(np.mean(base_avg)),
                 "helix_to_native": float(np.mean(base_helix)),
                 "pool_pairwise": float(np.mean(pool_pair)),
                 "pool_to_avg": float(np.mean(pool_to_avg))}


def main():
    out = {}
    p1(out)
    p2(out)
    with open(os.path.join(RESULTS, "verify_csteer.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote s16/results/verify_csteer.json")


if __name__ == "__main__":
    main()
