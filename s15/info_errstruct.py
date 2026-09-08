"""SPRINT 15, INFO, PART C.3 -- why the i.i.d. phase surface over-predicts every real
channel by 0.6-2.9 A, and what that does to the sprint-14 arithmetic.

`s15.info_regime` found that reading a real emitter's measured per-torsion RMS angular
error off the i.i.d. surface predicts a CA-RMSD that is 0.6 A (retrieval circular mean) to
2.9 A (the ORACLE best pool window: 64.4 deg RMS error, 1.77 A actual, 4.71 A predicted)
WORSE than what the emitter really produces.  Per-torsion angular error is therefore **not
a sufficient statistic for CA-RMSD**, and the sprint-14 line "per-torsion sigma required
for 2.0 A (i.i.d. error) = 15.1 deg" prices a channel that no real emitter resembles.

This module prices the structure directly by SURROGATE DESTRUCTION.  Take the real error
field `e = wrap(channel - native)` of an emitter, destroy one correlation at a time, add
the surrogate back to the native torsions, rebuild and re-measure:

    real            e itself                     -- identity check, must equal the emitter
    perm_res        residue order shuffled within the chain, (phi,psi) kept paired
    split_pair      phi and psi errors shuffled INDEPENDENTLY over residues
    gauss_matched   i.i.d. wrapped normal at the channel's own per-torsion sd
    boot_marginal   i.i.d. resample from the channel's pooled error distribution

The correlations themselves are reported too: same-residue r(dphi_i, dpsi_i), the
peptide-plane compensating pair r(dpsi_i, dphi_i+1), and the chain lag-1 terms.

    python -m s15.info_errstruct
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s15 import info_lib as L            # noqa: E402
from s15.info_regime import channels, wrap   # noqa: E402
import peptide_db as pdb                  # noqa: E402

REPS = 24


def surrogate(kind, ep, es, rng, sd_p, sd_s, pool_p, pool_s, R):
    n = len(ep)
    if kind == "real":
        return np.tile(ep, (R, 1)), np.tile(es, (R, 1))
    if kind == "perm_res":
        o = np.array([rng.permutation(n) for _ in range(R)])
        return ep[o], es[o]
    if kind == "split_pair":
        o1 = np.array([rng.permutation(n) for _ in range(R)])
        o2 = np.array([rng.permutation(n) for _ in range(R)])
        return ep[o1], es[o2]
    if kind == "gauss_matched":
        return rng.normal(0, sd_p, (R, n)), rng.normal(0, sd_s, (R, n))
    if kind == "boot_marginal":
        return (rng.choice(pool_p, (R, n)), rng.choice(pool_s, (R, n)))
    if kind == "rot_res":
        # cyclic rotation: keeps the whole adjacency structure of the error field and
        # only breaks its ALIGNMENT with the chain's position-dependent cost profile
        sh = rng.integers(1, max(2, n), R)
        idx = (np.arange(n)[None, :] + sh[:, None]) % n
        return ep[idx], es[idx]
    if kind == "sign_flip":
        # keeps |error| at every position EXACTLY, destroys only the directions
        return (ep[None, :] * rng.choice([-1.0, 1.0], (R, n)),
                es[None, :] * rng.choice([-1.0, 1.0], (R, n)))
    if kind == "sign_flip_paired":
        s = rng.choice([-1.0, 1.0], (R, n))
        return ep[None, :] * s, es[None, :] * s
    raise ValueError(kind)


KINDS = ["real", "rot_res", "sign_flip_paired", "sign_flip", "perm_res", "split_pair",
         "gauss_matched", "boot_marginal"]


def main():
    tg = I.targets()
    acc = {}
    corr = {}
    for t in tg:
        p = t["pdb"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        nat = np.asarray(I.load_univ(p)["nat_ca"], float)
        for nm, (phi, psi) in channels(t).items():
            ep = wrap(np.asarray(phi, float) - nphi)
            es = wrap(np.asarray(psi, float) - npsi)
            c = corr.setdefault(nm, {"same": [], "comp": [], "pp": [], "ss": []})
            c["same"].append((ep, es))
            acc.setdefault(nm, []).append((p, ep, es, nat, nphi, npsi))

    out = {"channels": {}}
    print("PART C.3 -- ERROR STRUCTURE, NOT ERROR MAGNITUDE\n")
    print("Correlation structure of the per-residue angular error, pooled over 126 "
          "targets:")
    print(f"{'channel':<34}{'r(dphi_i,dpsi_i)':>18}{'r(dpsi_i,dphi_i+1)':>20}"
          f"{'r(dphi lag1)':>14}{'r(dpsi lag1)':>14}")
    for nm, rows in acc.items():
        a = np.concatenate([r[1] for r in rows]); b = np.concatenate([r[2] for r in rows])
        same = L.pearson(a, b)
        cp, cn, lp, ls = [], [], [], []
        for _, ep, es, *_ in rows:
            if len(ep) > 2:
                cp.append(es[:-1]); cn.append(ep[1:])
                lp.append((ep[:-1], ep[1:])); ls.append((es[:-1], es[1:]))
        comp = L.pearson(np.concatenate(cp), np.concatenate(cn))
        lag_p = L.pearson(np.concatenate([x for x, _ in lp]),
                          np.concatenate([y for _, y in lp]))
        lag_s = L.pearson(np.concatenate([x for x, _ in ls]),
                          np.concatenate([y for _, y in ls]))
        out["channels"][nm] = {"r_same_residue": same, "r_compensating_pair": comp,
                               "r_lag1_phi": lag_p, "r_lag1_psi": lag_s}
        print(f"{nm:<34}{same:>18.3f}{comp:>20.3f}{lag_p:>14.3f}{ls and lag_s:>14.3f}")

    print(f"\nSURROGATE DESTRUCTION -- mean CA-RMSD over 126 targets, {REPS} reps\n")
    print(f"{'channel':<34}" + "".join(f"{k:>16}" for k in KINDS))
    for nm, rows in acc.items():
        allp = np.concatenate([r[1] for r in rows])
        alls = np.concatenate([r[2] for r in rows])
        sd_p, sd_s = float(allp.std()), float(alls.std())
        vals = {k: [] for k in KINDS}
        for a, (p, ep, es, nat, nphi, npsi) in enumerate(rows):
            rng = np.random.default_rng(97 * a + 13)
            for k in KINDS:
                R = 1 if k == "real" else REPS
                dp, ds = surrogate(k, ep, es, rng, sd_p, sd_s, allp, alls, R)
                ca = I.build_ca(nphi[None, :] + dp, npsi[None, :] + ds)
                vals[k].append(I.kabsch_rmsd_batch(ca, nat).mean())
        out["channels"][nm]["surrogates"] = {k: float(np.mean(v)) for k, v in vals.items()}
        out["channels"][nm]["per_target"] = {k: [float(x) for x in v]
                                             for k, v in vals.items()}
        print(f"{nm:<34}" + "".join(f"{np.mean(vals[k]):>16.3f}" for k in KINDS))

    print("\nPaired cost of destroying each correlation (vs the real error field):")
    for nm, v in out["channels"].items():
        pt = v["per_target"]
        real = np.asarray(pt["real"], float)
        line = [f"  {nm:<34}"]
        for k in KINDS[1:]:
            pr = L.report(np.asarray(pt[k], float), real)
            line.append(f"{k} {pr['mean_diff']:+.3f}[{pr['ci95'][0]:+.3f},"
                        f"{pr['ci95'][1]:+.3f}] ")
            v.setdefault("vs_real", {})[k] = pr
        print("".join(line))
    L.jwrite("info_errstruct", out)
    return out


if __name__ == "__main__":
    main()
