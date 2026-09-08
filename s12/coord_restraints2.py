"""COORDINATOR EXPERIMENT 5b -- the gap-filling arm the first pass got wrong.

`coord_restraints.py` measured the ceiling of a torsion-restraint (chemical-shift) channel
and found the headline: at sigma = 12 degrees on EVERY residue, simply building the chain
from the restrained torsions emits **1.475 A** with 86% of targets under 2.0 A, against the
shipped 3.213 -- and on the FAIL18 1.757 against 6.019.

But its coverage sweep is unfair and must not be quoted.  When a residue had no restraint
that pass set its torsions to ZERO, which is not "no information", it is a specific and
absurd conformation; torsion error compounds along the chain, so one nonsense residue in the
middle destroys everything downstream.  That is why `build` collapses from 1.475 (frac 1.0)
to 3.634 (frac 0.75).  It is an artefact of the fallback, not a property of partial coverage,
and coverage is exactly the axis that decides whether the modality is deployable -- TALOS-N
declines to predict for a real fraction of residues.

This module measures the honest version.  Unrestrained residues are filled from the
retrieval pool rather than from nothing, which is what a real hybrid system would do:

    fill_pool     pick the pool member whose torsions best agree with the RESTRAINED
                  residues, and take its torsions at the unrestrained ones.  Deployable in
                  shape: the restraints choose the donor, the pool supplies the gaps.
    fill_rama     fill unrestrained residues from the fold's Ramachandran mode (the most
                  common basin).  A pool-free control -- if this matches `fill_pool`, the
                  pool is contributing nothing and the restraints are doing all the work.
    fill_zero     the first pass's fallback, carried so the artefact is visible rather than
                  silently corrected.

ORACLE / DIAGNOSTIC throughout.  It reads native torsions to synthesise the restraints and
prices a FUTURE modality; it is not a method and is not deployable as written.  The leakage
caveat from `coord_restraints.py` stands: a shift-restrained predictor solves a different
problem from a sequence-only one and must be reported in its own column.

    python -m s12.coord_restraints2
"""
import os
import sys
import json
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s12.coord_restraints import _wrap, torsion_agreement, D2R   # noqa: E402
import peptide_db as pdb                   # noqa: E402

SIGMAS = (6.0, 12.0, 20.0)
FRACS = (1.0, 0.9, 0.75, 0.5, 0.25)
SEEDS = (0, 1, 2)


def run_target(t):
    u = I.load_univ(t["pdb"])
    n, seq, fold, nat = t["n"], t["seq"], t["fold"], u["nat_ca"]
    pool = I.pool_idx(u, I.K)
    W = u["W"][pool].astype(np.float32).astype(float)
    PHI = u["PHI"][pool]; PSI = u["PSI"][pool]
    dg = I.distogram(t["pdb"], seq, fold)
    i, j = I.pair_index(n)
    sc = I.shipped_score(dg, I.pair_dists(W, i, j))
    base_sub = np.argsort(sc, kind="stable")[:I.M]
    C, _ = I.coordinate_average(W[base_sub])
    base = I.ca_rmsd(I.project(C, seq, fold, lam=0.3)["ca"], nat)

    #: the pool's own Ramachandran mode, deployable and native-free
    mode_phi = float(np.median(PHI)); mode_psi = float(np.median(PSI))

    p = pdb.by_pdb(t["pdb"])
    phi0 = np.asarray(p.phi, float); psi0 = np.asarray(p.psi, float)   # ORACLE

    res = {"pdb": t["pdb"], "n": n, "fold": fold, "base": base, "arms": {}}
    for sig in SIGMAS:
        for frac in FRACS:
            acc = {"fill_pool": [], "fill_rama": [], "fill_zero": []}
            for s in SEEDS:
                rng = np.random.default_rng(int(1e6 * sig + 1e4 * frac * 100 + s
                                                + hash(t["pdb"]) % 9973))
                ph = _wrap(phi0 + rng.normal(0, sig * D2R, n))
                ps = _wrap(psi0 + rng.normal(0, sig * D2R, n))
                mask = rng.random(n) < frac
                if not mask.any():
                    mask[rng.integers(0, n)] = True

                # --- donor chosen by agreement on the RESTRAINED residues only
                ag = torsion_agreement(PHI, PSI, ph, ps, mask)
                d = int(np.argmin(ag))
                php = np.where(mask, ph, PHI[d]); psp = np.where(mask, ps, PSI[d])
                acc["fill_pool"].append(I.ca_rmsd(I.build_ca(php, psp), nat))

                phr = np.where(mask, ph, mode_phi); psr = np.where(mask, ps, mode_psi)
                acc["fill_rama"].append(I.ca_rmsd(I.build_ca(phr, psr), nat))

                phz = np.where(mask, ph, 0.0); psz = np.where(mask, ps, 0.0)
                acc["fill_zero"].append(I.ca_rmsd(I.build_ca(phz, psz), nat))
            res["arms"][f"s{sig:g}_f{frac:g}"] = {k: float(np.mean(v)) for k, v in acc.items()}
    return res


def main():
    tg = I.targets()
    path = os.path.join(I.RESULTS, "coord_restraints2.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t))
        if k % 10 == 0 or k == len(tg) - 1:
            I.write("coord_restraints2", {"what": "gap-filling for a torsion-restraint channel",
                                          "per_target": rows}, n_expected=len(tg))
            a = rows[-1]["arms"]["s12_f0.75"]
            print(f"  {k+1}/{len(tg)} {t['pdb']} base={rows[-1]['base']:.2f} | s12 f0.75: "
                  f"pool={a['fill_pool']:.2f} rama={a['fill_rama']:.2f} "
                  f"zero={a['fill_zero']:.2f}  [{time.time()-t0:.0f}s]", flush=True)
    I.write("coord_restraints2", {"what": "gap-filling for a torsion-restraint channel",
                                  "per_target": rows}, n_expected=len(tg))
    report(rows)


def report(rows):
    names = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    base = np.array([r["base"] for r in rows], float)
    out = {"n": len(rows), "base": float(base.mean()),
           "base_FAIL18": float(base[isf].mean()), "base_other108": float(base[~isf].mean()),
           "arms": {}}
    print(f"\nbase {base.mean():.4f}  FAIL18 {base[isf].mean():.4f}  "
          f"other108 {base[~isf].mean():.4f}\n")
    hdr = (f"{'cell':12s} {'fill_pool':>10s} {'fill_rama':>10s} {'fill_zero':>10s}  "
           f"{'pool: d vs base':>22s} {'<2A':>5s} {'F18':>6s}")
    print(hdr); print("-" * len(hdr))
    for key in rows[0]["arms"]:
        cell = {}
        for arm in ("fill_pool", "fill_rama", "fill_zero"):
            v = np.array([r["arms"][key][arm] for r in rows], float)
            cell[arm] = {"mean": float(v.mean()), "FAIL18": float(v[isf].mean()),
                         "other108": float(v[~isf].mean()),
                         "frac_under_2": float((v < 2.0).mean()),
                         "paired": I.paired(v, base, folds=folds, names=names)}
        out["arms"][key] = cell
        d = cell["fill_pool"]["paired"]
        print(f"{key:12s} {cell['fill_pool']['mean']:10.3f} {cell['fill_rama']['mean']:10.3f} "
              f"{cell['fill_zero']['mean']:10.3f}  {d['mean_diff']:+8.3f} "
              f"[{d['ci95'][0]:+.3f},{d['ci95'][1]:+.3f}] "
              f"{cell['fill_pool']['frac_under_2']:5.2f} {cell['fill_pool']['FAIL18']:6.2f}")
    I.write("coord_restraints2_report", out)


if __name__ == "__main__":
    main()
