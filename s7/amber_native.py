"""Does the native minimise the AMBER energy on real pools?

THE QUESTION. Finding 5 applied a native-percentile test to the distance objective and got a
devastating answer: on 126 tuning targets the native is the objective's argmin on only
**3/126** targets and sits at the **36.8th percentile** of the pool's score distribution. On
123 of 126 targets something scores better than the truth. Re-weighting cannot fix it (a full
sweep of the `1/(sd+0.5)^gamma` rule and per-shell weights moved selection by at most
-0.025 A, inside the 0.147 A SE), and a *perfect* distance oracle on these pools still only
selects 1.994 A. The distance channel is therefore exhausted as a selector.

All-atom Amber is a genuinely different channel, it is an architecturally mandatory component,
and the repo's own records disagree about it in a way this test can settle:

  FOR      `physics-ranks-real-geometry-not-lattice`: Amber in-band rho +0.320 against the
           distogram's +0.091 -- 3.5x the skill -- and "Sprint 4's negative was a build
           confound".
  AGAINST  `nothing-ranks-within-the-pool`: all-atom Amber reranking is worth +0.004 A. But
           that was the 24-target dev set at SE ~0.354, plausibly through a shortlist that
           `in-band-is-the-only-ranking-metric` says caps any reranker at 2.44 A.
  CAUTION  `amber-native-antiranking-distributed` and `amoeba-does-not-fix-antiranking` both
           report Amber failing to rank the native lowest -- but on DECOY BANKS, and
           `decoy-bank-not-a-pool-proxy` established that instrument does not transfer.

So: the anti-ranking claim has never been tested on REAL POOLS with adequate power. That is
what this does.

If Amber places the native near the 0th percentile where distances place it at the 37th,
Amber is a correctly-specified objective and the pipeline should pivot to it. If Amber ALSO
places the native mid-pool, then no available scoring channel ranks truth first, which is the
most important negative result of the program.

The confound this controls, which is the whole reason the old result is suspect. Pool
candidates are IDEAL-GEOMETRY backbones rebuilt from a real fragment's torsions via
`build_backbone_batch`. A native read from a PDB is an EXPERIMENTAL structure with its own
bond lengths and angles. Scoring one against the other compares geometry provenance, not
conformation -- exactly the "build confound" that is recorded as having invalidated Sprint 4's
negative. So the native is scored TWICE:

  native_raw       the experimental coordinates as deposited
  native_rebuilt   the native's OWN (phi, psi) pushed through the SAME builder as every
                   candidate, so it is geometrically commensurable with the pool

`native_rebuilt` is the honest comparator and its percentile is the headline number.
`native_raw` is reported alongside to quantify how large the provenance effect is.

ORACLE DIAGNOSTIC. Every native quantity here is oracle information. This module is a
diagnostic and is imported by nothing deployable. It reads natives on purpose, to answer
whether the objective is correctly specified; it never ranks a candidate using one.

Scoring recipe is the verified one from `s5/physics_probe.py`: `amber_refine.refine_coords`
(NOT `refine`, which would project continuous torsions onto the discrete state library),
`k_restraint=K_MODERATE`, `steps=0`, `components=True`, energy = nonbonded + solvation, with
the same bond+angle > 1000 strain rejection for numerical failures.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.stats import spearmanr

import distogram as dgm
import legacy_field as lfd
import peptide_db as db
import protein_geometry as geo
import torsion_lib2 as tl2
from s5.lib import B62, DEV, encode, kabsch_rmsd_batch, windows_full

K = int(os.environ.get("K", "100"))
K_AMBER = int(os.environ.get("K_AMBER", "40"))
BAND = 1.5
OUT = os.environ.get("OUT", "s7/amber_native.json")
_FRAG = {}


def pool_for(fold, target_seq):
    if fold not in _FRAG:
        _FRAG[fold] = list(dgm._fold_fragments(fold, 5))
    folds = db.folds(5)
    return ([q for q in db.load() if folds[q.seq] != fold and q.seq != target_seq]
            + _FRAG[fold])


def amber_energy(seq, rep, coords_single, ar):
    """nonbonded + solvation for one structure, or nan on numerical failure."""
    out = ar.refine_coords(seq, rep, coords_single, k_restraint=ar.K_MODERATE,
                           steps=0, tolerance=1.0, components=True)
    cm = out["components"]
    if cm["bond"] + cm["angle"] > 1000.0:
        return float("nan"), True
    return float(cm["nonbonded"] + cm["solvation"]), False


def measure(pid):
    import amber_refine as ar
    t0 = time.time()
    p = db.by_pdb(pid)
    fold = db.folds(5)[p.seq]
    CA, S, PHI, PSI, _ = windows_full(pool_for(fold, p.seq), p.n)

    # retrieval: BLOSUM62, sequence only -- no native involved
    sim = B62[S, encode(p.seq)[None, :]].sum(1)
    idx = np.argsort(-sim)[:K]
    phi, psi = PHI[idx].astype(float), PSI[idx].astype(float)
    coords = geo.build_backbone_batch(phi, psi)
    Cca = coords["CA"]

    tab = tl2.library_for(p.seq, 8, p.seq)          # excludes self
    rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)

    ka = min(K_AMBER, len(idx))
    amber = np.full(ka, np.nan)
    strain = 0
    for b in range(ka):
        try:
            e, bad = amber_energy(p.seq, rep, {k: v[b] for k, v in coords.items()}, ar)
            amber[b] = e
            strain += int(bad)
        except Exception:                                            # noqa: BLE001
            pass

    # ---- ORACLE DIAGNOSTIC from here on: natives are read deliberately ----
    # native_rebuilt: the native's own torsions through the SAME builder as the pool, so
    # geometry provenance is held constant. This is the honest comparator.
    nat_rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))
    try:
        e_rb, _ = amber_energy(p.seq, rep, nat_rb, ar)
    except Exception:                                                # noqa: BLE001
        e_rb = float("nan")
    # native_raw is NOT available: `peptide_db.Peptide` exposes only (ca, phi, psi, rebuild)
    # -- no experimental all-atom coordinates -- so the deposited structure cannot be scored
    # by a force field that needs a full topology. `native_rebuilt` was always the honest
    # comparator anyway (it holds geometry provenance constant against the pool); the size of
    # the provenance effect is quantified by `rebuild_gap` instead.
    e_raw = float("nan")
    ar.clear_cache()

    rms = kabsch_rmsd_batch(Cca, p.ca)[:ka]
    ok = np.isfinite(amber)
    rec = {"pdb": pid, "n": int(p.n), "fold": int(fold), "K": int(K), "ka": int(ka),
           "amber_ok": int(ok.sum()), "strain_rejects": int(strain),
           "pool_best": float(rms.min()), "pool_mean": float(rms.mean()),
           "e_native_rebuilt": e_rb, "e_native_raw": e_raw,
           "rebuild_gap": float(np.abs(kabsch_rmsd_batch(nat_rb["CA"][None], p.ca)[0])),
           "wall": round(time.time() - t0, 1)}

    # Compactness control. Amber's nonbonded + solvation strongly rewards compactness, and a
    # compact native would make energy correlate with RMSD through radius of gyration alone --
    # a geometry proxy, not physics discrimination. This is the same confound that reduced the
    # inter-generator agreement lead from +0.594 to -0.07 once controlled, so it is measured
    # here rather than assumed away. `rho_partial_rg` is the number to trust.
    rg_all = np.sqrt(((Cca - Cca.mean(1, keepdims=True)) ** 2).sum(2).mean(1))[:ka]

    if ok.sum() >= 5:
        a, r, g = amber[ok], rms[ok], rg_all[ok]
        rec["amber_sel"] = float(r[int(np.argmin(a))])
        rec["amber_rho"] = float(spearmanr(a, r).statistic)
        rec["rg_rho"] = float(spearmanr(g, r).statistic)
        rec["amber_rg_rho"] = float(spearmanr(a, g).statistic)
        # rank-linear partial correlation of (amber, rmsd) controlling for rg
        def _rk(x):
            s = spearmanr(x, x).statistic  # noqa: F841  (forces scipy import path)
            o = np.argsort(np.argsort(x)).astype(float)
            return (o - o.mean()) / max(o.std(), 1e-12)
        ra, rr_, rgz = _rk(a), _rk(r), _rk(g)
        p_ar, p_ag, p_rg = (float(np.mean(ra * rr_)), float(np.mean(ra * rgz)),
                            float(np.mean(rr_ * rgz)))
        den = np.sqrt(max((1 - p_ag ** 2) * (1 - p_rg ** 2), 1e-12))
        rec["rho_partial_rg"] = float((p_ar - p_ag * p_rg) / den)
        band = r <= r.min() + BAND
        rec["amber_rho_band"] = (float(spearmanr(a[band], r[band]).statistic)
                                 if band.sum() > 4 else None)
        for tag, e in (("rebuilt", e_rb), ("raw", e_raw)):
            if np.isfinite(e):
                rec[f"native_{tag}_pct"] = float(100.0 * (a < e).mean())
                rec[f"native_{tag}_is_argmin"] = bool(e < a.min())
    return rec


def main():
    pids = os.environ.get("PIDS", "").split(",") if os.environ.get("PIDS") else DEV
    rows = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["pdb"] for r in rows}
    for pid in pids:
        if pid in done:
            continue
        try:
            rec = measure(pid)
        except Exception as exc:                                     # noqa: BLE001
            rec = {"pdb": pid, "error": f"{type(exc).__name__}: {exc}"}
        rows.append(rec)
        tmp = OUT + ".tmp"
        json.dump(rows, open(tmp, "w"), indent=1)
        os.replace(tmp, OUT)
        if "error" in rec:
            print(f"{pid:6} ERROR {rec['error'][:70]}", flush=True)
        else:
            print(f"{pid:6} ok {rec['amber_ok']:3d}/{rec['ka']:3d}  "
                  f"pool {rec['pool_best']:5.2f}  "
                  f"amber_sel {rec.get('amber_sel', float('nan')):5.2f}  "
                  f"rho {rec.get('amber_rho', float('nan')):+.3f}  "
                  f"prg {rec.get('rho_partial_rg', float('nan')):+.3f}  "
                  f"band {rec.get('amber_rho_band') if rec.get('amber_rho_band') is None else round(rec['amber_rho_band'],3)}  "
                  f"nat_pct rb {rec.get('native_rebuilt_pct')} raw {rec.get('native_raw_pct')} "
                  f"({rec['wall']:.0f}s)", flush=True)

    v = [r for r in rows if "amber_sel" in r]
    if not v:
        print("\nno scored targets")
        return
    f = lambda k: float(np.mean([r[k] for r in v if r.get(k) is not None]))
    print(f"\n=== {len(v)} targets, K={K}, K_amber<={K_AMBER} ===")
    print(f"pool best {f('pool_best'):.3f}   pool mean {f('pool_mean'):.3f}   "
          f"amber selected {f('amber_sel'):.3f}")
    print(f"amber rho global {f('amber_rho'):+.3f}   in-band {f('amber_rho_band'):+.3f}")
    print(f"\nCOMPACTNESS CONTROL (the number to trust):")
    print(f"  rho(amber, rmsd)          {f('amber_rho'):+.3f}")
    print(f"  rho(rg, rmsd)             {f('rg_rho'):+.3f}   <- rg alone as a ranker")
    print(f"  rho(amber, rg)            {f('amber_rg_rho'):+.3f}")
    print(f"  PARTIAL, controlling rg   {f('rho_partial_rg'):+.3f}")
    for tag in ("rebuilt", "raw"):
        pct = [r[f"native_{tag}_pct"] for r in v if r.get(f"native_{tag}_pct") is not None]
        am = [r[f"native_{tag}_is_argmin"] for r in v
              if r.get(f"native_{tag}_is_argmin") is not None]
        if pct:
            print(f"\nNATIVE ({tag}) percentile in the Amber energy distribution:")
            print(f"  mean {np.mean(pct):.1f}   median {np.median(pct):.1f}   "
                  f"argmin on {sum(am)}/{len(am)} targets")
    print("\nthe distance objective, same test, 126 targets: argmin 3/126, "
          "mean percentile 36.8")
    print(f"mean rebuild gap (native torsions -> builder) {f('rebuild_gap'):.3f} A")


if __name__ == "__main__":
    main()
