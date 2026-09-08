"""COORDINATOR EXPERIMENT 5 -- pricing the sprint's top-ranked next experiment BEFORE it is built.

The dossier ranks "BMRB chemical shifts -> TALOS-N-style phi/psi restraints" as the single
most promising untested direction: 79 of the 126 tuning targets and 10 of the 18 failures
have deposited shifts, and TALOS-N converts them to backbone torsions at roughly 12 degrees
of accuracy.  Building that pipeline is days of work -- NMR-STAR parsing, residue mapping,
a torsion predictor -- and the project's own history says the expensive thing to do is build
a channel before measuring what it could possibly be worth.

So this module measures the CEILING of the whole modality without touching the BMRB at all.
It takes the NATIVE torsions, corrupts them to the accuracy a real shift-based predictor
delivers, and asks what the pipeline emits.  If native phi/psi known to +-12 degrees does
not get under 2.0 A, the direction is dead and nobody should build it.  If it does, the
prize is quantified and the build is justified.

EVERYTHING HERE IS ORACLE / DIAGNOSTIC.  It reads native torsions and is not deployable.
It prices a FUTURE modality; it is not a method and must never be quoted as one.

THE LEAKAGE POINT THAT MUST TRAVEL WITH ANY RESULT.  Chemical shifts are an independent
experimental observable -- they are measured, not computed from coordinates -- so using them
is not reading the answer.  But the deposited structure was DETERMINED using those shifts
together with NOEs, so a shift-restrained predictor is solving a different problem from a
sequence-only one.  Any future arm must be reported as NMR-restrained prediction, in its own
column, never as an improvement to the sequence-only system.

FOUR ARMS, each at several noise levels (0, 6, 12, 20, 30, 45 degrees):

    build     build the chain directly from the noisy torsions.  This is what a torsion
              predictor buys if you simply trust it.  Torsion error compounds along the
              chain, so this is expected to be poor at large n -- that expectation is
              exactly what needs measuring rather than assuming.
    filter    score the K=500 pool by agreement with the noisy torsions, take the top 75,
              then the UNCHANGED synthesis (coordinate average -> projection).  This is the
              deployable SHAPE of a shift channel: a better filter, not a better builder.
    fuse      the distogram's z-ranked score plus the torsion agreement's, equally weighted.
    project   project the noisy-torsion chain onto the pool's convex geometry by using it as
              the projection target instead of the coordinate average.

    python -m s12.coord_restraints
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
import peptide_db as pdb                   # noqa: E402

#: Per-torsion standard deviation, degrees.  12 is TALOS-N's published accuracy on
#: well-determined residues; 20-30 is what it delivers on flexible ones, and TALOS-N
#: declines to predict at all for a fraction of residues, which the `frac` sweep covers.
SIGMAS = (0.0, 6.0, 12.0, 20.0, 30.0, 45.0)
#: THE SWEEP IS SCORED ON THE RAW COORDINATE AVERAGE, NOT THROUGH THE PROJECTION.
#: The objective agent validated the raw top-m average as a proxy for the emitted structure
#: at r = 0.995 with a near-constant +0.16 A offset, and a full projection sweep here would
#: be ~324 L-BFGS solves per target.  The headline cells (sigma = 12, frac = 1.0) and the
#: baseline ARE run through the real multi-start projection, and both are reported.
#: Fraction of residues for which a restraint is available at all.  TALOS-N is typically
#: confident on ~70-90% of residues; the rest fall back to the pool's own information.
FRACS = (1.0, 0.75, 0.5)
SEEDS = (0, 1, 2)
D2R = np.pi / 180.0


def _wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def torsion_agreement(PHI, PSI, phi_t, psi_t, mask):
    """Mean absolute circular deviation of each pool member from the restraint. Lower better."""
    dphi = np.abs(_wrap(PHI - phi_t[None, :]))
    dpsi = np.abs(_wrap(PSI - psi_t[None, :]))
    w = mask.astype(float)
    if w.sum() == 0:
        return np.zeros(len(PHI))
    return ((dphi + dpsi) * w[None, :]).sum(1) / (2.0 * w.sum())


def _zrank(x):
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(r.std(), 1e-12)


def run_target(t):
    u = I.load_univ(t["pdb"])
    n, seq, fold, nat = t["n"], t["seq"], t["fold"], u["nat_ca"]
    pool = I.pool_idx(u, I.K)
    W = u["W"][pool].astype(np.float32).astype(float)
    PHI = u["PHI"][pool]; PSI = u["PSI"][pool]
    rr = u["rr"][pool]
    dg = I.distogram(t["pdb"], seq, fold)
    i, j = I.pair_index(n)
    sc = I.shipped_score(dg, I.pair_dists(W, i, j))
    zsc = _zrank(sc)

    p = pdb.by_pdb(t["pdb"])
    phi0 = np.asarray(p.phi, float); psi0 = np.asarray(p.psi, float)   # ORACLE

    def raw(idx):
        C, _ = I.coordinate_average(W[idx])
        return I.ca_rmsd(C, nat)

    def emit(idx):
        C, _ = I.coordinate_average(W[idx])
        o = I.project(C, seq, fold, lam=0.3)
        return I.ca_rmsd(o["ca"], nat)

    base_sub = np.argsort(sc, kind="stable")[:I.M]
    res = {"pdb": t["pdb"], "n": n, "fold": fold,
           "base_raw": raw(base_sub), "base_arm": emit(base_sub),
           "pool_best": float(rr.min()), "arms": {}}

    for sig in SIGMAS:
        for frac in FRACS:
            keyb, keyf, keyu, keyp = [], [], [], []
            for s in SEEDS:
                rng = np.random.default_rng(int(1e6 * sig + 1e3 * frac * 100 + s
                                                + hash(t["pdb"]) % 9973))
                ph = _wrap(phi0 + rng.normal(0, sig * D2R, n))
                ps = _wrap(psi0 + rng.normal(0, sig * D2R, n))
                mask = rng.random(n) < frac
                ph = np.where(mask, ph, phi0 * 0.0)      # unrestrained residues carry no info
                ps = np.where(mask, ps, psi0 * 0.0)
                # ---- build: trust the torsions directly
                ca = I.build_ca(np.where(mask, ph, 0.0), np.where(mask, ps, 0.0)) \
                    if mask.all() else I.build_ca(ph, ps)
                keyb.append(I.ca_rmsd(ca, nat))
                # ---- filter: torsion agreement as the pool score
                ag = torsion_agreement(PHI, PSI, ph, ps, mask)
                keyf.append(raw(np.argsort(ag, kind="stable")[:I.M]))
                # ---- fuse with the distogram
                fu = zsc + _zrank(ag)
                keyu.append(raw(np.argsort(fu, kind="stable")[:I.M]))
                # ---- project: the restrained chain IS already ideal geometry, so its
                # "projection" is itself; carried so the column has a defined meaning
                keyp.append(keyb[-1])
            res["arms"][f"s{sig:g}_f{frac:g}"] = {
                "build": float(np.mean(keyb)), "filter": float(np.mean(keyf)),
                "fuse": float(np.mean(keyu)), "project": float(np.mean(keyp))}

    #: THE HEADLINE CELL, through the REAL multi-start projection, so the sweep's proxy is
    #: anchored to an emitted number rather than trusted.
    rng = np.random.default_rng(4242 + hash(t["pdb"]) % 9973)
    ph = _wrap(phi0 + rng.normal(0, 12.0 * D2R, n))
    ps = _wrap(psi0 + rng.normal(0, 12.0 * D2R, n))
    mask = np.ones(n, bool)
    ag = torsion_agreement(PHI, PSI, ph, ps, mask)
    res["emitted_s12"] = {
        "filter": emit(np.argsort(ag, kind="stable")[:I.M]),
        "fuse": emit(np.argsort(zsc + _zrank(ag), kind="stable")[:I.M]),
        "build": I.ca_rmsd(I.build_ca(ph, ps), nat)}
    return res


def main():
    tg = I.targets()
    path = os.path.join(I.RESULTS, "coord_restraints.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t))
        if k % 5 == 0 or k == len(tg) - 1:
            json.dump({"what": "ORACLE ceiling of a torsion-restraint (chemical-shift) channel",
                       "per_target": rows}, open(path, "w"), indent=1)
            a = rows[-1]["arms"]["s12_f1"]; e = rows[-1]["emitted_s12"]
            print(f"  {k+1}/{len(tg)} {t['pdb']} base={rows[-1]['base_arm']:.2f} "
                  f"| s12 raw build={a['build']:.2f} filt={a['filter']:.2f} "
                  f"fuse={a['fuse']:.2f} | EMITTED filt={e['filter']:.2f} "
                  f"fuse={e['fuse']:.2f}  [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "ORACLE ceiling of a torsion-restraint (chemical-shift) channel",
               "per_target": rows}, open(path, "w"), indent=1)
    report(rows)


def report(rows):
    names = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    base = np.array([r["base_arm"] for r in rows], float)
    out = {"n": len(rows), "base": I.summary(base),
           "base_FAIL18": float(base[isf].mean()), "base_other108": float(base[~isf].mean()),
           "arms": {}}
    print(f"\nbase (shipped synthesis, lam=0.3): {base.mean():.4f}  "
          f"FAIL18 {base[isf].mean():.4f}  other108 {base[~isf].mean():.4f}\n")
    hdr = f"{'restraint':14s} {'build':>8s} {'filter':>8s} {'fuse':>8s} {'project':>8s}   {'best d vs base':>15s}"
    print(hdr); print("-" * len(hdr))
    for key in rows[0]["arms"]:
        cell = {}
        for arm in ("build", "filter", "fuse", "project"):
            v = np.array([r["arms"][key][arm] for r in rows], float)
            cell[arm] = {"mean": float(v.mean()), "FAIL18": float(v[isf].mean()),
                         "other108": float(v[~isf].mean()), "frac_under_2": float((v < 2.0).mean()),
                         "paired": I.paired(v, base, folds=folds, names=names)}
        out["arms"][key] = cell
        bestarm = min(cell, key=lambda a: cell[a]["mean"])
        d = cell[bestarm]["paired"]
        print(f"{key:14s} {cell['build']['mean']:8.3f} {cell['filter']['mean']:8.3f} "
              f"{cell['fuse']['mean']:8.3f} {cell['project']['mean']:8.3f}   "
              f"{bestarm}: {d['mean_diff']:+.3f} [{d['ci95'][0]:+.3f},{d['ci95'][1]:+.3f}]")
    I.write("coord_restraints_report", out)
    print("\nFAIL18 / other108 at sigma=12, frac=1.0:")
    c = out["arms"]["s12_f1"]
    for arm in ("build", "filter", "fuse", "project"):
        print(f"  {arm:8s} all {c[arm]['mean']:.3f}   FAIL18 {c[arm]['FAIL18']:.3f}   "
              f"other108 {c[arm]['other108']:.3f}   <2A {c[arm]['frac_under_2']:.3f}")


if __name__ == "__main__":
    main()
