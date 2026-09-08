"""s18/phys_down.py -- PHASE 8: THE CAUSAL DOWNSTREAM COMPARISON, ON IDENTICAL STRUCTURES.

Pre-registration: `s18/PREREG_phys.md` section 3, written before this module produced a number.

THE QUESTION.  Sprint 17 measured every physics filter as a RANKER or a GATE and scored it on
ranking statistics.  It never asked the only question a filter's user cares about:

    if you put this filter in front of the operator that actually emits the answer,
    does the ANSWER get better?

That is a causal question about the whole downstream chain, and the ledger says it need not agree
with the ranking statistics at all:

  * `d_out = 1.16 * d_set_mean + 0.04 * d_set_best`, R2 0.89 -- the terminal averaging operator
    consumes the set MEAN, not the set BEST.
  * Legacy's gate at f = 0.50 RAISES near-native recall (+0.112 [+0.023, +0.198]) and DESTROYS
    the set best (+0.219 [+0.031, +0.444]).

Those two together make a REGISTERED PREDICTION, written before the run: **through an averaging
operator, Legacy's gate should HELP**, because the quantity it improves (the class, hence the set
mean) is the quantity the operator reads, and the quantity it destroys (the single best member)
enters at a weight of 0.04.  If that prediction fails, the operator law is wrong or the gate's
recall gain does not reach the mean.  Either outcome is a result.

THE INSTRUMENT.  Every filter operates on the SAME 75 windows -- the shipped top-75 the deployed
operator averages (`s14.avgspace.top75_windows`) -- so no arm is handed a set built by its own
criterion, and the comparison across filters is exact.

ARMS (the brief's list, in its order):

    none            no filter                                the incumbent
    leg_torsion     the one Legacy gate that raises recall without damaging the set best
    leg_contact     the Miyazawa-Jernigan term this sprint is about
    legacy          the combined eleven-term total at DEFAULT_WEIGHTS (never fitted)
    amber_sp        genuine ff14SB/GBn2 single point, one per candidate
    rand            MATCHED-RANDOM gate of the same count, 3 `stable_rng` draws   [CONTROL]
    helix           ZERO-INFORMATION: the constant ideal alpha-helix               [CONTROL]

DOWNSTREAM OPERATORS applied to each filter's survivor set:

    avg    the raw all-atom coordinate average            (the 3.048 A structure)
    proj   its ideal-geometry projection                  (the 3.213 A incumbent)
    amb30  AMBER restrained repair, N/CA/C at k = 30      (the earned role)

so `legacy -> amb30` and `rand -> amb30` are the brief's `Legacy->AMBER` and
`matched-random gate->AMBER` arms, and every AMBER arm has **its own gated input** available by
construction because the input is recorded per arm.

WHAT IS NEVER DONE HERE.  AMBER's validity statistics are printed beside the RMSD table for the
SAME structures and neither is quoted without the other.  Legacy's energy correlations are
printed beside the selection delta against a matched-random gate of the same count.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s15.phys_repl import averaged_backbone_from, project_both    # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s16 import energy_lib as EL                                  # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402

ATOMS = ("N", "CA", "C", "O", "CB")
NEAR = 2.0                       # "near-native", angstrom
GATE_F = (0.50, 0.25)            # 0.50 PRIMARY (pre-registered), 0.25 as shape
N_RAND = 3                       # matched-random draws
K_RESTRAINT = 30.0               # the incumbent N/CA/C protocol
AMB_STEPS = 0
AMB_TOL = 1.0

#: Which filters get a full AMBER k = 30 repair.  Pre-registered in PREREG_phys.md section 3.
AMBER_ARMS = ("none", "legacy", "leg_torsion", "leg_contact", "rand")

FILTERS = ("none", "leg_torsion", "leg_contact", "legacy", "amber_sp", "rand", "helix")


def _keep(score, n_keep):
    """Indices of the `n_keep` LOWEST-energy candidates, ties broken deterministically by index."""
    s = np.asarray(score, float)
    s = np.where(np.isfinite(s), s, np.inf)
    return np.sort(np.lexsort((np.arange(len(s)), s))[:n_keep])


def _ensemble_stats(idx, d, W):
    """Everything about a SURVIVOR SET that does not need a downstream operator."""
    idx = np.asarray(idx, int)
    dd = np.asarray(d, float)[idx]
    return {"m": int(len(idx)),
            "set_mean": float(dd.mean()),
            "set_best": float(dd.min()),
            "set_median": float(np.median(dd)),
            "n_near": int((dd < NEAR).sum()),
            "diversity": PL.diversity(np.asarray(W, float)[idx])}


def _emit(W, PHI, PSI, idx, seq, fold, nat, do_proj=True):
    """The two cheap downstream operators on one survivor set: raw average and projection.

    The projection costs ~5 s (multi-start, maxiter 300) against the average's 0.03 s, so it is
    run on the PRIMARY f = 0.50 arms only; the f = 0.25 arms are shape and carry the average.
    """
    idx = np.asarray(idx, int)
    avg, C_ca, dev = averaged_backbone_from(W[idx], PHI[idx], PSI[idx])
    out = {"avg_rmsd": float(I.ca_rmsd(avg["CA"], nat)),
           "ca_identity_max_dev": float(dev),
           "avg_valid": EL.panel({a: np.asarray(avg[a], float) for a in ATOMS if a in avg}, seq),
           "avg_ca_spacing": float(np.linalg.norm(np.diff(np.asarray(avg["CA"], float), axis=0),
                                                  axis=1).mean())}
    if do_proj:
        from core import geometry as geo
        pr = project_both(C_ca, seq, fold)
        bb = geo.build_backbone(pr["phi"], pr["psi"])
        out["proj_rmsd"] = float(I.ca_rmsd(pr["ca"], nat))
        out["proj_valid"] = EL.panel({a: np.asarray(bb[a], float) for a in ATOMS if a in bb}, seq)
    return out, avg


# ------------------------------------------------------------------ one target
def run_target(t, want_amber=True):
    from core import amber as am
    import torsion_lib2 as tl2
    from s17 import phys_lib as P17

    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(W)

    #: ORACLE -- labels and ceilings only, never a predictor.
    d = np.asarray(I.kabsch_rmsd_batch(W, nat), float)

    #: NATIVE-FREE scores on the identical 75 structures.
    comp = EL.legacy_components_of_windows(seq, PHI, PSI)
    legacy = np.asarray(EL.legacy_total_from(comp), float)
    scores = {"leg_torsion": np.asarray(comp["torsion"], float),
              "leg_contact": np.asarray(comp["contact"], float),
              "legacy": legacy}

    #: ZERO-INFORMATION reference score: distance of each candidate to a constant alpha-helix.
    hb = P17.helix_backbone(n)
    scores["helix"] = np.asarray(I.kabsch_rmsd_batch(W, np.asarray(hb["CA"], float)), float)

    rec = {"pdb": pdb, "n": n, "fold": fold, "K": int(K),
           "pool_mean": float(d.mean()), "pool_best": float(d.min()),
           "n_near_pool": int((d < NEAR).sum()),
           "leg_terms": {k: np.asarray(v, float).tolist() for k, v in comp.items()},
           "d": d.tolist()}

    #: AMBER SINGLE POINT, one per candidate, on the same ideal-geometry rebuilds.
    box = None
    if want_amber:
        from core import geometry as geo
        tab = tl2.library_for(seq, 4, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        box = P17.ConstrainedBox(seq, rep)
        E = np.full(K, np.nan)
        t0 = time.time()
        for b in range(K):
            c = geo.build_backbone(PHI[b], PSI[b])
            E[b] = box.energy_point({a: np.asarray(v, float) for a, v in c.items()})["energy"]
        rec["amber_sp_wall"] = round(time.time() - t0, 2)
        scores["amber_sp"] = E
        rec["s_amber_sp"] = E.tolist()
        box.close()

    #: ---- the arms.  Every filter sees the SAME 75 structures; only the survivor set differs.
    rng = SD.stable_rng(pdb, "s18down")
    arms, arm_idx = {}, {}
    for f in GATE_F:
        n_keep = max(2, int(round(K * (1.0 - f))))
        for name in FILTERS:
            key = name if name == "none" else f"{name}@{f:.2f}"
            if key in arms:
                continue
            if name == "rand":
                idxs = [np.sort(rng.permutation(K)[:n_keep]) for _ in range(N_RAND)]
            elif name == "none":
                idxs = [np.arange(K)]           # the unfiltered arm does not depend on f
            elif name in scores:
                idxs = [_keep(scores[name], n_keep)]
            else:
                continue
            sub = []
            for r, idx in enumerate(idxs):
                e = _ensemble_stats(idx, d, W)
                #: PRIMARY f only, and for the matched-random arm only its first draw, get the
                #: 5 s projection; every draw gets every cheap statistic.
                em, _avg = _emit(W, PHI, PSI, idx, seq, fold, nat,
                                 do_proj=(f == GATE_F[0] and r == 0))
                e.update(em)
                sub.append(e)
            #: a matched-random arm is the MEAN over its draws; the draws are kept beside it.
            #: `proj_rmsd` lives on draw 0 only (see `_emit`), so it is carried, not averaged.
            agg = {}
            for k in sub[0]:
                vals = [s[k] for s in sub if k in s]
                agg[k] = (float(np.mean(vals)) if isinstance(sub[0][k], (int, float))
                          and len(vals) == len(sub) else sub[0][k])
            if len(sub) > 1:
                agg["draws"] = sub
            arms[key] = agg
            arm_idx[key] = [np.asarray(x, int) for x in idxs]

    #: ---- AMBER k = 30 restrained repair.  Each arm carries ITS OWN input, so a gated arm is
    #: never compared to the ungated instrument-wide input (the Sprint-17 near-miss).
    if want_amber:
        tab = tl2.library_for(seq, 4, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        for name in AMBER_ARMS:
            key = name if name == "none" else f"{name}@{GATE_F[0]:.2f}"
            if key not in arm_idx:
                continue
            outs = []
            #: draw 0 only for the matched-random arm -- an AMBER repair is ~8 s and the gate's
            #: matched-random control is already averaged over 3 draws on every cheap statistic.
            for idx in arm_idx[key][:1]:
                avg, _C, _dev = averaged_backbone_from(W[idx], PHI[idx], PSI[idx])
                coords = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}
                t0 = time.time()
                r = am.refine_coords(seq, rep, coords, k_restraint=K_RESTRAINT,
                                     steps=AMB_STEPS, tolerance=AMB_TOL, threads=1, memo=False)
                ca = np.asarray(r["ca"], float)
                bb = {a: np.asarray(r["backbone"][a], float) for a in ATOMS if a in r["backbone"]}
                outs.append({
                    "rmsd": float(I.ca_rmsd(ca, nat)),
                    "input_rmsd": float(I.ca_rmsd(coords["CA"], nat)),
                    "energy": float(r["energy"]),
                    "energy_initial": float(r["energy_initial"]),
                    "converged": bool(r["converged"]),
                    "ca_disp": float(I.ca_rmsd(ca, coords["CA"])),
                    "restraint_rmsd": float(r["restraint_rmsd"]),
                    "valid": EL.panel(bb, seq),
                    "input_valid": EL.panel(coords, seq),
                    "wall": round(time.time() - t0, 2)})
            arms[key]["amb30"] = outs[0] if len(outs) == 1 else {
                k: (float(np.mean([o[k] for o in outs])) if isinstance(outs[0][k], (int, float))
                    else outs[0][k]) for k in outs[0]}
            if len(outs) > 1:
                arms[key]["amb30"]["draws"] = outs
                #: a draw fails the gate independently; keep every flag, not the mean of them.
                arms[key]["amb30"]["converged"] = bool(all(o["converged"] for o in outs))
                arms[key]["amb30"]["energy"] = float(max(o["energy"] for o in outs))
    rec["arms"] = arms
    return rec


# ------------------------------------------------------------------ driver
def run(targets=None, out="down.json", want_amber=True, verbose=True):
    tg = targets if targets is not None else I.targets()
    path = os.path.join(PL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    cfg = {"K_RESTRAINT": K_RESTRAINT, "AMB_STEPS": AMB_STEPS, "GATE_F": list(GATE_F),
           "N_RAND": N_RAND, "NEAR": NEAR, "FILTERS": list(FILTERS),
           "AMBER_ARMS": list(AMBER_ARMS), "want_amber": bool(want_amber)}
    for c, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        PL.mem_ok(1.6)
        rows.append(run_target(t, want_amber=want_amber))
        if verbose:
            a = rows[-1]["arms"]
            print(f"  {len(rows)}/{len(tg)} {t['pdb']}  none={a['none']['avg_rmsd']:.3f} "
                  f"leg@.5={a.get('legacy@0.50',{}).get('avg_rmsd',float('nan')):.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        PL.write(out, {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg)},
                 n_expected=len(tg))
    PL.write(out, {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg)}, n_expected=len(tg))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-amber", action="store_true")
    ap.add_argument("--out", default="down.json")
    a = ap.parse_args()
    tg = I.targets()
    if a.limit:
        tg = tg[:a.limit]
    run(tg, out=a.out, want_amber=not a.no_amber)
