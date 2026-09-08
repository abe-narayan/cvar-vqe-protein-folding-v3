"""SPRINT 15, INFO, PART C.2 -- where the project actually sits on the phase diagram,
and whether the diagram predicts real emitters.

The phase surface is built from synthetic error models.  It is only useful if a REAL
channel's measured per-torsion error, read off the surface, predicts that channel's real
CA-RMSD.  This module

  1. measures every available torsion channel's per-torsion angular error directly
     (circular RMS and MAE of `wrap(pred - native)`, phi and psi separately),
  2. reads the primary surface at (coverage = 1, sigma = the channel's RMS error),
  3. compares that prediction with the channel's ACTUAL emitted CA-RMSD, and
  4. marks the regime on the plane.

Channels are native-free emitters; the angular errors and CA-RMSDs are post-hoc
evaluation, so this is an ORACLE-EVALUATED but not ORACLE-CONDITIONED measurement, except
for the two arms whose names say ORACLE.

    python -m s15.info_regime
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s15 import info_lib as L            # noqa: E402
from s14 import retprior as RP           # noqa: E402
import peptide_db as pdb                  # noqa: E402

ALPHA = (np.radians(-63.0), np.radians(-42.0))


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def channels(t):
    """dict name -> (phi, psi) for one target.  Native-free unless the name says ORACLE."""
    p = t["pdb"]
    u = I.load_univ(p)
    n = int(t["n"])
    out = {}
    P5, S5, sim5 = RP.windows(p, "pool")
    P7, S7, sim7 = RP.windows(p, "top75")
    out["retrieval pool-500 circular mean"] = (RP.circ_mean(P5, None, 0),
                                               RP.circ_mean(S5, None, 0))
    out["retrieval top-75 circular mean"] = (RP.circ_mean(P7, None, 0),
                                             RP.circ_mean(S7, None, 0))
    w = RP._softmax_w(sim7, 2.0)
    out["retrieval top-75 sim-weighted"] = (RP.circ_mean(P7, w, 0), RP.circ_mean(S7, w, 0))
    # a single random retrieved window, and the pool medoid (native-free)
    rng = np.random.default_rng(hash(p) % 2 ** 31)
    j = int(rng.integers(0, len(P7)))
    out["one random top-75 window"] = (P7[j], S7[j])
    W = np.asarray(u["W"], float)[I.pool_idx(u)]
    rec = I.shipped_record(p)
    sub = np.asarray(rec["sub"], int)
    Pm = I.pairwise_rmsd(W[sub])
    m = I.medoid(Pm)
    out["top-75 medoid window"] = (P7[m], S7[m])
    out["constant alpha-helix (control)"] = (np.full(n, ALPHA[0]), np.full(n, ALPHA[1]))
    out["incumbent projected torsions"] = (np.asarray(rec["phi"], float),
                                           np.asarray(rec["psi"], float))
    # ORACLE arms
    rr = np.asarray(u["rr"], float)[I.pool_idx(u)]
    out["ORACLE best pool window"] = (np.asarray(u["PHI"], float)[I.pool_idx(u)][int(rr.argmin())],
                                      np.asarray(u["PSI"], float)[I.pool_idx(u)][int(rr.argmin())])
    return out


def main():
    with open(os.path.join(L.RESULTS, "info_phase.json")) as fh:
        ph = json.load(fh)
    sig = np.asarray(ph["sigmas"], float)
    grid = np.asarray(ph["surfaces"]["primary_iid_uniform_pool"]["grid"])[:, -1]  # cov = 1

    rows = {}
    for t in I.targets():
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        nat = np.asarray(I.load_univ(p)["nat_ca"], float)
        for nm, (phi, psi) in channels(t).items():
            phi = np.asarray(phi, float); psi = np.asarray(psi, float)
            ep, es = wrap(phi - nphi), wrap(psi - npsi)
            r = rows.setdefault(nm, {"ep": [], "es": [], "rmsd": []})
            r["ep"].append(ep); r["es"].append(es)
            r["rmsd"].append(I.ca_rmsd(I.build_ca(phi, psi), nat))

    out = {"channels": {}}
    print("PART C.2 -- WHERE THE PROJECT SITS ON THE PHASE DIAGRAM\n")
    print(f"{'channel':<34}{'RMS phi':>9}{'RMS psi':>9}{'RMS':>8}{'MAE':>8}"
          f"{'predRMSD':>10}{'actual':>9}{'gap':>8}")
    for nm, r in rows.items():
        ep = np.degrees(np.concatenate(r["ep"]))
        es = np.degrees(np.concatenate(r["es"]))
        rp = float(np.sqrt((ep ** 2).mean())); rs = float(np.sqrt((es ** 2).mean()))
        rms = float(np.sqrt(((ep ** 2).mean() + (es ** 2).mean()) / 2))
        mae = float((np.abs(ep).mean() + np.abs(es).mean()) / 2)
        act = float(np.mean(r["rmsd"]))
        pred = float(np.interp(rms, sig, grid))
        # lag-1 autocorrelation of the per-residue error, the coherence check
        lag = []
        for e in r["ep"] + r["es"]:
            if len(e) > 2 and np.std(e) > 1e-9:
                lag.append(np.corrcoef(e[:-1], e[1:])[0, 1])
        out["channels"][nm] = {"rms_phi_deg": rp, "rms_psi_deg": rs, "rms_deg": rms,
                               "mae_deg": mae, "actual_rmsd": act, "predicted_rmsd": pred,
                               "gap": act - pred,
                               "lag1_mean": float(np.nanmean(lag)) if lag else None,
                               "n": len(r["rmsd"])}
        print(f"{nm:<34}{rp:>9.1f}{rs:>9.1f}{rms:>8.1f}{mae:>8.1f}"
              f"{pred:>10.3f}{act:>9.3f}{act - pred:>+8.3f}")

    print("\nlag-1 autocorrelation of the per-residue angular error (coherence check):")
    for nm, v in out["channels"].items():
        print(f"  {nm:<34}{v['lag1_mean']:+.3f}" if v["lag1_mean"] is not None else nm)

    inc = 3.2040761603809194
    k = np.searchsorted(-grid, -inc)
    out["incumbent_sigma_equivalent_deg"] = float(np.interp(-inc, -grid[::-1], sig[::-1]))
    print(f"\nsigma-equivalent of the incumbent pipeline (3.204 A through the surface at "
          f"coverage 1): {out['incumbent_sigma_equivalent_deg']:.1f} deg "
          f"(sprint 14 recorded 27.1 deg)")
    L.jwrite("info_regime", out)
    return out


if __name__ == "__main__":
    main()
