"""SPRINT 16, ENERGY -- THE CAUSAL ABLATION.  What does Legacy contribute?  What does
AMBER contribute?  Do they interact?

ONE FIXED GENERATED ENSEMBLE.  Every arm consumes the SAME 75 windows per target -- the
shipped distogram filter's top-75 set, `s14.avgspace.top75_windows`, a pure cache read
with no RNG.  Nothing about the ensemble changes between arms; only the readout does.

THE DESIGN IS A 2 x 2 FACTORIAL WITH A MATCHED CLASSICAL CONTROL.

                            AMBER off                 AMBER on
    Legacy off       ctrl : avg -> projection    amb  : avg -> ff14SB/GBn2 k=30
    Legacy on        leg  : filter -> avg -> proj legamb : filter -> avg -> AMBER
    random control   rnd  : rndfilter -> avg -> proj rndamb: rndfilter -> avg -> AMBER

    LEGACY  = leg  - ctrl                 (main effect of Legacy, AMBER off)
    AMBER   = amb  - ctrl                 (main effect of AMBER,  Legacy off)
    INTER   = (legamb - amb) - (leg - ctrl)   (interaction)
    and every Legacy row is priced against `rnd`, the matched-count random filter, so a
    "Legacy helps" claim has to beat DROPPING THE SAME NUMBER OF WINDOWS AT RANDOM.

PRE-REGISTERED, WRITTEN BEFORE ANY FILTERED ARM'S RMSD EXISTED
--------------------------------------------------------------
* LEGACY-ON is `drop the worst quartile (ceil(0.25 * m) = 19 of 75) of windows by the
  genuine eleven-component Legacy total at DEFAULT_WEIGHTS`.  Native-free: it reads the
  target sequence and each window's torsions and nothing else.  The quartile is a fixed
  fraction, not a tuned threshold, and it was fixed before any RMSD was looked at.
* The choice of a RANK filter over Legacy's documented CLASH-GATE role is forced by a
  native-free measurement taken first and recorded in `s16/results/energy_ablate.json`
  as `steric_gate_stats`: the Legacy steric term is non-zero on only 0-12% of windows
  (mean ~4%), so a pure clash gate cannot move a 75-window coordinate average.  That
  measurement involves no native quantity.  The steric gate is measured anyway, as a
  DETECTOR, in `s16/energy_repair.py`, which is the axis where it has a defensible role.
* RANDOM control drops a uniformly random 19 of 75 (2 independent draws, `stable_rng`).
  Two draws because this arm DOES contain an RNG, so the sprint's multi-start
  false-positive floor CAN reach it -- unlike the AMBER k = 30 path, which has no RNG.
* SUCCESS for "Legacy contributes": mean CA-RMSD of `leg` beats `ctrl` by more than
  `rnd` beats `ctrl`, with a paired fold-aware bootstrap CI on (leg - rnd) excluding
  zero at n = 126.
* SUCCESS for "AMBER contributes to ACCURACY": (amb - ctrl) CI excludes zero on the
  CONVERGENCE-GATED set, with the excluded count printed.
* VALIDITY is scored on a separate axis and never substituted for accuracy: the
  stereochemistry panel (`energy_lib.panel`) on the SAME structure whose RMSD is quoted.

COST, at the right number.  `refine_coords` at k = 30 is a full restrained minimisation
measured at 12.57 s mean on this box -- ~1500x an AMBER single point (8-14 ms) and ~2.8x
the projection (4.46 s).  Four AMBER minimisations and four projections per target is
~68 s x 126 targets.

    python -m s16.energy_ablate --shard 0/2
    python -m s16.energy_ablate --shard 1/2
    python -m s16.energy_ablate --report
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import argparse
import json
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                       # noqa: E402
from s15.seed import stable_rng                       # noqa: E402
from s15.phys_repl import averaged_backbone_from, project_both   # noqa: E402
from s14.avgspace import top75_windows                # noqa: E402
from s16 import energy_lib as L                       # noqa: E402

K = 30.0
DROP_FRAC = 0.25                     # pre-registered
RND_DRAWS = (1, 2)                   # pre-registered
ARMS = ("ctrl", "amb", "leg", "legamb", "rnd1", "rnd1amb", "rnd2", "rnd2amb")


def keep_indices(m, drop, order):
    """Keep the m - drop windows that `order` puts first.  `order` is an argsort of the
    filter score (ascending = keep)."""
    return np.sort(np.asarray(order, int)[: m - drop])


def _panel_of_backbone(phi, psi, seq):
    from core import geometry as geo
    bb = geo.build_backbone(np.asarray(phi, float), np.asarray(psi, float))
    c = {a: np.asarray(bb[a], float) for a in L.ATOMS if a in bb}
    return L.panel(c, seq), c


def run_target(t, k=K):
    from core import amber as am
    import torsion_lib2 as tl2
    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    # One bounded yield per TARGET, not per AMBER call: the box is shared with sibling
    # agents, and a per-call hold cost ~6 minutes of pure waiting per target.
    L.cpu_hold(96.0, pdb, max_wait=30.0)
    W, PHI, PSI, u = top75_windows(pdb)
    nat = u["nat_ca"]
    m = len(W)
    drop = int(np.ceil(DROP_FRAC * m))

    # ---- the filters.  All native-free.
    comp = L.legacy_components_of_windows(seq, PHI, PSI)
    tot = L.legacy_total_from(comp)
    keep = {"ctrl": np.arange(m),
            "leg": keep_indices(m, drop, np.argsort(tot, kind="mergesort"))}
    for d in RND_DRAWS:
        rng = stable_rng("s16energy", "rndfilter", int(d), pdb)
        keep[f"rnd{d}"] = np.sort(rng.permutation(m)[: m - drop])

    rec = dict(pdb=pdb, fold=fold, n=len(seq), n_windows=int(m), n_drop=int(drop),
               ORACLE_pool_best_rmsd=float(np.min(I.kabsch_rmsd_batch(W, nat))),
               legacy_steric_nonzero=float((comp["steric"] > 1e-9).mean()),
               legacy_steric_max=float(comp["steric"].max()),
               legacy_total_med=float(np.median(tot)),
               ORACLE_window_rmsd_mean=float(np.mean(I.kabsch_rmsd_batch(W, nat))))
    # ORACLE DIAGNOSTIC: what the filter did to the ensemble's true quality.
    wr = I.kabsch_rmsd_batch(W, nat)
    for tag, ix in keep.items():
        rec[f"ORACLE_setmean_{tag}"] = float(wr[ix].mean())
        rec[f"ORACLE_setbest_{tag}"] = float(wr[ix].min())

    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)

    for tag, ix in keep.items():
        avg, C_ca, dev = averaged_backbone_from(np.asarray(W)[ix],
                                                np.asarray(PHI)[ix],
                                                np.asarray(PSI)[ix])
        rec[f"{tag}_setdiv"] = float(np.mean(I.pairwise_rmsd(np.asarray(W)[ix])))
        rec[f"{tag}_rawavg_rmsd"] = float(I.ca_rmsd(avg["CA"], nat))
        rec[f"{tag}_rawavg_panel"] = L.panel(avg, seq)

        # -- AMBER off: the ideal-geometry projection (the incumbent operator)
        t0 = time.time()
        pr = project_both(C_ca, seq, fold)
        pan, _ = _panel_of_backbone(pr["fit_phi"], pr["fit_psi"], seq)
        rec[tag] = dict(rmsd=float(I.ca_rmsd(pr["fit_ca"], nat)), panel=pan,
                        wall=round(time.time() - t0, 2))

        # -- AMBER on: restrained ff14SB/GBn2 relaxation of the same average
        t0 = time.time()
        L.mem_hold(1.6, pdb)
        try:
            r = am.refine_coords(seq, rep, avg, k_restraint=float(k), steps=0,
                                 tolerance=1.0, threads=1, memo=False)
        except NotImplementedError as exc:
            rec[tag + "amber"] = dict(error=str(exc)[:120])
            continue
        ca = np.asarray(r["ca"], float)
        bb = {a: np.asarray(r["backbone"][a], float) for a in L.ATOMS
              if a in r["backbone"]}
        rec[tag + "amber"] = dict(
            rmsd=float(I.ca_rmsd(ca, nat)), energy=float(r["energy"]),
            energy_initial=float(r["energy_initial"]),
            converged=bool(r["converged"]), converge_reason=r["converge_reason"],
            restraint_rmsd=float(r["restraint_rmsd"]),
            moved_ca=float(I.ca_rmsd(ca, avg["CA"])),
            panel=L.panel(bb, seq), wall=round(time.time() - t0, 2))
    return rec


# ------------------------------------------------------------------ driver
def _path(shard):
    return os.path.join(L.RESULTS, f"energy_ablate_s{shard}.json")


def _done_everywhere():
    """Every pdb already finished by ANY shard process.

    Re-read before each target so that extra worker processes can be added mid-run and
    cooperate instead of duplicating.  The race window is one target; duplicates are
    de-duplicated at read time by `load_all`.
    """
    done = set()
    for f in __import__("glob").glob(os.path.join(L.RESULTS, "energy_ablate_s*.json")):
        try:
            with open(f) as fh:
                done |= {r["pdb"] for r in json.load(fh).get("per_target", [])}
        except Exception:
            pass
    return done


def run(shard=0, nshard=1, k=K):
    tg = [t for i, t in enumerate(I.targets()) if i % nshard == shard]
    p = _path(shard)
    rows = []
    if os.path.exists(p):
        with open(p) as fh:
            rows = json.load(fh).get("per_target", [])
    t0 = time.time()
    for t in tg:
        if t["pdb"] in _done_everywhere():
            continue
        rows.append(run_target(t, k))
        r = rows[-1]
        print(f"[abl {shard} {len(rows)}/{len(tg)}] {r['pdb']} "
              f"ctrl {r['ctrl']['rmsd']:.3f} amb {r['ctrlamber'].get('rmsd', float('nan')):.3f} "
              f"leg {r['leg']['rmsd']:.3f} legamb {r['legamber'].get('rmsd', float('nan')):.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)
        _write(shard, rows, len(tg))      # checkpoint EVERY target: a mid-run
        #                                   MemoryError cost this run 8 targets once
    _write(shard, rows, len(tg))
    return rows


def _write(shard, rows, n_expected):
    with open(_path(shard), "w") as fh:
        json.dump(dict(what="S16 ENERGY causal ablation: Legacy x AMBER on one fixed ensemble",
                       k=K, drop_frac=DROP_FRAC, per_target=rows,
                       n_rows=len(rows), n_expected=int(n_expected)), fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def load_all():
    """Every shard's rows, DE-DUPLICATED by pdb (cooperating workers can race)."""
    seen, rows = {}, []
    for f in sorted(__import__("glob").glob(os.path.join(L.RESULTS, "energy_ablate_s*.json"))):
        with open(f) as fh:
            for r in json.load(fh)["per_target"]:
                seen[r["pdb"]] = r
    rows = sorted(seen.values(), key=lambda r: r["pdb"])
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", default="0/1")
    args = ap.parse_args()
    s, n = args.shard.split("/")
    run(int(s), int(n))
