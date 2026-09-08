"""SPRINT 14, coordinator -- WHERE a torsion error falls decides what it costs.

WHY.  Two things are now measured and neither explains the other.

  * At MATCHED marginal sigma, along-chain error coherence moves emitted RMSD by up to
    0.93 A (`s14/coherence.py`): at sigma 12 deg, anti-correlated 1.272, i.i.d. 1.617,
    positively autocorrelated 2.203.  Coherence is worth as much as accuracy.
  * But every REAL native-free emitter measured sits in the near-i.i.d. band, lag-1
    autocorrelation -0.09 to +0.07.  So coherence does NOT explain why the top-75
    similarity-weighted arm (phi 38.5 deg -> 3.514 A) beats the top-75 circular mean
    (phi 33.6 deg -> 4.072 A).  That anomaly is still unexplained.

THE REMAINING CANDIDATE.  Mean absolute error weights every residue equally.  A backbone
builder does not.  A torsion error at residue 2 of 13 rotates the entire downstream chain
about that vertex; the same error at residue 12 displaces one residue.  So the cost of an
error should fall steeply with sequence position, MAE should be a badly weighted summary of
what actually matters, and two predictors with the same MAE can differ by a lot depending on
WHERE their errors sit.

WHAT THIS DECIDES.  This is not bookkeeping.  It sets the value of partial coverage, which
is the make-or-break parameter for the chemical-shift route: if the N-terminal third
dominates, then coverage there is worth several times coverage elsewhere, and both the
coverage thresholds and the kill threshold for that route have to be re-derived
position-weighted.  It also predicts the sign of the Sprint 13 correction that terminal
dropout is 0.40-0.50 A CHEAPER than uniform -- C-terminal residues should be nearly free.

ORACLE DIAGNOSTIC.  Perturbs native torsions to price an error model.  Not a method.

Run:
    python -m s14.position
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s14 import coherence as CO            # noqa: E402
from s14 import ladder as L                # noqa: E402
from s14 import retprior as R              # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

wrap = CO.wrap


def single_residue_cost(sigma_deg=20.0, n_rep=24, targets=None):
    """ORACLE.  Perturb ONE residue's phi (and separately psi), everything else native.

    Reported on a fractional position axis so chains of different length are comparable.
    """
    tg = targets if targets is not None else I.targets()
    rows = []
    for t in tg:
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        nphi, npsi = CO.native_torsions(pdb)
        n = len(nphi)
        base = I.ca_rmsd(I.build_ca(nphi, npsi), nat)
        s = np.deg2rad(sigma_deg)
        for i in range(n):
            rng = np.random.default_rng(hash((pdb, i)) % (2 ** 32))
            dp, ds = [], []
            for _ in range(n_rep):
                phi = nphi.copy(); phi[i] = wrap(phi[i] + rng.normal() * s)
                dp.append(I.ca_rmsd(I.build_ca(phi, npsi), nat))
                psi = npsi.copy(); psi[i] = wrap(psi[i] + rng.normal() * s)
                ds.append(I.ca_rmsd(I.build_ca(nphi, psi), nat))
            rows.append({"pdb": pdb, "n": n, "i": i, "frac": i / max(n - 1, 1),
                         "base": float(base),
                         "d_phi": float(np.mean(dp) - base),
                         "d_psi": float(np.mean(ds) - base)})
    return rows


def bin_by_fraction(rows, nbin=10):
    f = np.asarray([r["frac"] for r in rows], float)
    dp = np.asarray([r["d_phi"] for r in rows], float)
    ds = np.asarray([r["d_psi"] for r in rows], float)
    edges = np.linspace(0, 1, nbin + 1)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (f >= a) & (f <= b if b == 1.0 else f < b)
        if m.sum():
            out.append({"lo": float(a), "hi": float(b), "n": int(m.sum()),
                        "d_phi": float(dp[m].mean()), "d_psi": float(ds[m].mean())})
    return out


def positional_weights(rows):
    """Normalised cost weight per fractional position -- what MAE SHOULD have used."""
    b = bin_by_fraction(rows, nbin=10)
    w = np.asarray([x["d_phi"] + x["d_psi"] for x in b], float)
    return (w / w.sum()).tolist(), b


def weighted_error_of_emitters(rows):
    """Re-score every real emitter's torsion error with the measured position weights.

    If position is the missing variable, a position-WEIGHTED error should order emitted
    RMSD where plain MAE does not.  This is the test; it can fail.
    """
    _, bins = positional_weights(rows)
    lo = np.asarray([x["lo"] for x in bins]); hi = np.asarray([x["hi"] for x in bins])
    wphi = np.asarray([x["d_phi"] for x in bins]); wpsi = np.asarray([x["d_psi"] for x in bins])

    emitters = dict(L.EMITTERS); emitters.update(R.EMITTERS)
    tg = I.targets()
    out = {}
    for name, (emit, _) in emitters.items():
        plain, weighted, rmsd = [], [], []
        for t in tg:
            pdb, n = t["pdb"], int(t["n"])
            u = I.load_univ(pdb)
            ephi, epsi = CO.error_field(emit, pdb, t["seq"], n, t["fold"])
            frac = np.arange(n) / max(n - 1, 1)
            b = np.clip(np.searchsorted(hi, frac, side="left"), 0, len(hi) - 1)
            ap = np.abs(np.rad2deg(ephi)); asx = np.abs(np.rad2deg(epsi))
            plain.append(float(np.mean(np.concatenate([ap[1:], asx[:-1]]))))
            num = float((ap * wphi[b]).sum() + (asx * wpsi[b]).sum())
            den = float(wphi[b].sum() + wpsi[b].sum())
            weighted.append(num / den)
            rng = np.random.default_rng(0)
            p2, s2 = emit(pdb, t["seq"], n, int(t["fold"]), rng)
            rmsd.append(I.ca_rmsd(I.build_ca(p2, s2), u["nat_ca"]))
        out[name] = {"plain_mae_deg": float(np.mean(plain)),
                     "position_weighted_deg": float(np.mean(weighted)),
                     "rmsd": float(np.mean(rmsd))}
    return out


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    return float((ra * rb).sum() / np.sqrt((ra * ra).sum() * (rb * rb).sum()))


def run():
    print("PART 1 -- cost of a single perturbed torsion, by fractional chain position")
    rows = single_residue_cost()
    bins = bin_by_fraction(rows)
    print(f"{'position':>12}{'d_phi (A)':>12}{'d_psi (A)':>12}{'n':>8}")
    for b in bins:
        print(f"{b['lo']:.1f}-{b['hi']:.1f}".rjust(12)
              + f"{b['d_phi']:>12.4f}{b['d_psi']:>12.4f}{b['n']:>8}")
    w, _ = positional_weights(rows)
    first, last = sum(w[:3]), sum(w[-3:])
    print(f"\nfirst 30% of the chain carries {first:.1%} of the total single-residue cost; "
          f"last 30% carries {last:.1%}   ratio {first / max(last, 1e-9):.2f}x")

    print("\nPART 2 -- does a position-weighted error order emitted RMSD where MAE fails?")
    we = weighted_error_of_emitters(rows)
    print(f"{'emitter':<28}{'RMSD':>8}{'plain MAE':>12}{'pos-weighted':>14}")
    for k, v in sorted(we.items(), key=lambda kv: kv[1]["rmsd"]):
        print(f"{k:<28}{v['rmsd']:>8.3f}{v['plain_mae_deg']:>12.1f}"
              f"{v['position_weighted_deg']:>14.1f}")
    names = list(we)
    rho_plain = _spearman([we[k]["plain_mae_deg"] for k in names],
                          [we[k]["rmsd"] for k in names])
    rho_w = _spearman([we[k]["position_weighted_deg"] for k in names],
                      [we[k]["rmsd"] for k in names])
    print(f"\nSpearman(error, emitted RMSD) across {len(names)} emitters:"
          f"   plain MAE {rho_plain:+.3f}    position-weighted {rho_w:+.3f}")

    out = {"per_residue": rows, "bins": bins, "weights": w,
           "front_back_ratio": float(first / max(last, 1e-9)),
           "emitters": we, "rho_plain_mae": rho_plain, "rho_position_weighted": rho_w}
    with open(os.path.join(RESULTS, "position.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_position", out, n_expected=len(I.targets()))
    return out


if __name__ == "__main__":
    run()
