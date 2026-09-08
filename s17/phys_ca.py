"""s17/phys_ca.py -- Ca-PRESERVING AMBER REPAIR, run as a constrained optimisation.

Sprint 17 sections 29-31.  Pre-registration in `s17/PREREG_phys.md` (E1), written before
this module produced a number.

    minimise  E_ff14SB/GBn2(x)   subject to   ||CA_i(x) - CA_i(x0)|| <= eps

swept over a ladder of eps, from **Ca completely fixed** (eps = 0 exactly, by zero particle
mass) through extremely tight / tight / moderate / loose / unrestrained.  The ladder is swept
in the DUAL variable k (see `phys_lib.LADDER` for why), and the frontier is plotted against
the REALISED Ca displacement -- measured, never assumed.

THE INPUT is PHYS's input, unchanged so the numbers are comparable across three sprints:
the raw all-atom coordinate average of the shipped top-75 window set
(`s14.avgspace.top75_windows` -> `s15.phys_repl.averaged_backbone_from`).  Pure cache read,
no RNG.  Its do-nothing Ca-RMSD is 3.0498 A over the 126 targets and it is genuinely broken
(5.6 clashes below 2.0 A, 0.315 rms relative bond strain).

CONTROLS CARRIED AT EVERY RUNG, both mandatory:
  * ZERO-INFORMATION `none`   -- do nothing.  The operator's own input.  It won accuracy in
                                Sprint 16 against both repairers by 0.13-0.16 A on 99/123.
  * ZERO-INFORMATION `helix`  -- a constant ideal alpha-helix.  Ramachandran 1.000 and zero
                                clashes BY CONSTRUCTION.  No validity number is quotable
                                without it: it beat AMBER 0W/65L on Ramachandran.
  * MATCHED RANDOM            -- an isotropic random Ca displacement at each rung's own
                                realised RMS displacement, 3 `stable_rng` draws.  Sprint 16
                                found it at least as accurate as AMBER's displacement.

CONVERGENCE GATE reported with its exclusion count and the excluded ids at every rung.
ROTATED-LAB-FRAME NULL reported with its MAXIMUM, on a pre-registered subsample, against the
unchanged band |mean| <= 0.005 A AND max <= 0.05 A.

PASSES (fixed here before the run, and why):
  P1  `cafix` + `k30` at FULL n = 126.  cafix is the pre-registered hypothesis; k30 is the
      incumbent and its presence in the same process is the instrument's own reproduction
      check against `s16/results/repair_A.json`.
  P2  the whole harmonic ladder on the PRE-REGISTERED, NATIVE-FREE, fold x length stratified
      30-target subsample `s16.repair.shape_subsample` -- reused verbatim, so no new
      selection freedom is created.  This is the FRONTIER, and it is a shape measurement:
      no rung it favours is re-quoted at full n as an out-of-sample choice.
  P3  the flat-bottom cross-check (eps = 0.10 / 0.25 A at k_f = 1e5) on the first 8 of the
      subsample, to show the two parameterisations trace one frontier.
  P4  the rotated-frame null at `cafix` on the first 12 of the subsample.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                     # noqa: E402
from s15.seed import stable_rng                     # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s17 import phys_lib as P                       # noqa: E402

P1 = ("cafix", "k30")
P2 = ("cafix", "ca1000", "ca300", "ca100", "ca30", "ca10", "ca3", "free", "k30")
P3 = ("eps010", "eps025")
N_XCHECK = 8
N_FRAME = 12


def subsample():
    from s16.repair import shape_subsample
    return shape_subsample()


def _sup(X, T):
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(T, float))[0]


def run_target(t, rungs, frame_rungs=(), verbose=False):
    from s14.avgspace import top75_windows
    from s15.phys_repl import averaged_backbone_from
    import torsion_lib2 as tl2

    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    nat = np.asarray(u["nat_ca"], float)
    avg, C_ca, dev = averaged_backbone_from(W, PHI, PSI)
    ca_in = np.asarray(avg["CA"], float)

    rec = {"pdb": pdb, "n": n, "fold": fold,
           "input_ca": ca_in.tolist(), "nat_ca": nat.tolist(),
           "none": {"rmsd": float(I.ca_rmsd(ca_in, nat)),
                    "validity": P.validity(avg, seq),
                    "ca_disp_rms": 0.0, "ca_disp_max": 0.0}}

    #: ZERO-INFORMATION reference.  rama 1.000 / 0 clashes BY CONSTRUCTION.
    hb = P.helix_backbone(n)
    hca = _sup(hb["CA"], ca_in)
    rec["helix"] = {"rmsd": float(I.ca_rmsd(hb["CA"], nat)),
                    "validity": P.validity(hb, seq),
                    "ca_disp_rms": float(np.sqrt(((hca - ca_in) ** 2).sum(1).mean())),
                    "ca_disp_max": float(np.sqrt(((hca - ca_in) ** 2).sum(1)).max())}

    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    box = P.ConstrainedBox(seq, rep)
    try:
        for rung in rungs:
            r = box.relax(avg, rung)
            ca = np.asarray(r["ca"], float)
            v = P.validity(r["backbone"], seq)
            d = {"rmsd": float(I.ca_rmsd(ca, nat)),
                 "energy": float(r["energy"]),
                 "energy_initial": float(r["energy_initial"]),
                 "converged": bool(r["converged"]),
                 "converge_reason": str(r.get("converge_reason", "")),
                 "ca_disp_rms": float(r["ca_disp_rms"]),
                 "ca_disp_max": float(r["ca_disp_max"]),
                 "ca_violation": float(r["ca_violation"]),
                 "moved_ca_rmsd": float(I.ca_rmsd(ca, ca_in)),
                 "eps": float(r["eps"]), "wall": round(float(r["wall"]), 2),
                 "validity": v, "ca": ca.tolist()}
            #: MATCHED RANDOM: same realised RMS displacement, isotropic, 3 draws.
            rr = []
            for s in range(P_NRAND):
                rng = stable_rng("s17phys", "matchrand", rung, pdb, s)
                cr = P.matched_random_disp(ca_in, float(r["ca_disp_rms"]), rng)
                rr.append(float(I.ca_rmsd(cr, nat)))
            d["rand_rmsd"] = rr
            d["rand_rmsd_mean"] = float(np.mean(rr))
            rec[rung] = d
            if verbose:
                print(f"    {pdb} {rung:>8} E {r['energy']:>10.1f} conv "
                      f"{str(r['converged']):>5} disp {r['ca_disp_rms']:.3f} "
                      f"rmsd {d['rmsd']:.4f} rama {v['rama_favoured']:.3f} "
                      f"cl2 {v['n_clash_2A']:.0f} bond {v['bond_strain']:.4f} "
                      f"({r['wall']:.1f}s)", flush=True)

        #: ROTATED-LAB-FRAME NULL -- zero by construction (proper rotations only).
        for rung in frame_rungs:
            Q, tv = P.random_rigid(stable_rng("s17phys", "frame", 1, pdb))
            rot = {a: np.asarray(v, float) @ Q.T + tv for a, v in avg.items()}
            r2 = box.relax(rot, rung)
            rec.setdefault("frame", {})[rung] = {
                "rmsd": float(I.ca_rmsd(np.asarray(r2["ca"], float), nat)),
                "energy": float(r2["energy"]), "converged": bool(r2["converged"])}
    finally:
        box.close()
    return rec


P_NRAND = 3


def _load(path):
    if os.path.exists(path):
        try:
            return json.load(open(path))["rows"]
        except Exception:
            pass
    return []


def run(tag, tg, rungs, frame_rungs=(), verbose=True):
    path = os.path.join(RESULTS, f"phys_ca_{tag}.json")
    rows = _load(path)
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for c, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        EL.mem_hold(1.6, tag=f"phys_ca {tag}")
        #: `core.amber.memory_guard` raises MemoryError above a 92% physical-memory ceiling
        #: and killed a shard of the Sprint 16 ablation mid-run; it killed pass 2 here too.
        #: Yield and retry rather than losing the shard -- the per-target checkpoint means a
        #: retry costs one target, and a give-up costs the pass.
        rec = None
        for attempt in range(6):
            try:
                rec = run_target(t, rungs, frame_rungs, verbose=verbose)
                break
            except MemoryError as ex:
                print(f"  [{tag}] {t['pdb']} MemoryError ({ex}); yielding "
                      f"(attempt {attempt+1}/6)", flush=True)
                try:
                    from core import amber as _am
                    _am.clear_cache()
                except Exception:
                    pass
                EL.mem_hold(2.0, tag=f"phys_ca {tag} retry", max_wait=900.0)
        if rec is None:
            print(f"  [{tag}] {t['pdb']} SKIPPED after 6 MemoryError attempts", flush=True)
            continue
        rows.append(rec)
        json.dump({"rows": rows, "rungs": list(rungs), "tag": tag,
                   "converge_max_kcal": 1000.0,
                   "n_expected": len(tg)}, open(path, "w"))
        if verbose:
            print(f"  [{tag}] {len(rows)}/{len(tg)} {t['pdb']} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="ps", default="1")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    tg = I.targets()
    sub = set(subsample())
    subt = [t for t in tg if t["pdb"] in sub]
    print(f"PRE-REGISTERED SUBSAMPLE (s16.repair.shape_subsample, native-free, "
          f"fold x length stratified, n = {len(subt)}):")
    print("  " + " ".join(t["pdb"] for t in subt), flush=True)

    if a.ps == "1":
        run("p1", tg[:a.limit] if a.limit else tg, P1)
    elif a.ps == "2":
        run("p2", subt, P2)
    elif a.ps == "3":
        run("p3", subt[:N_XCHECK], P3)
    elif a.ps == "4":
        run("p4", subt[:N_FRAME], ("cafix",), frame_rungs=("cafix",))
