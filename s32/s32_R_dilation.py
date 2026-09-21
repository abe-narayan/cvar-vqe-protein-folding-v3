#!/usr/bin/env python
"""S32 LANE R -- the scalar-dilation arms at the CLOUD basis, and the separation-dependence
of the averaging contraction.

    python s32/s32_R_dilation.py

PRODUCER FOR:
    s32/results/s32_R_dilation_cloud.json      (R-14, the REGISTERED P1.3 arm)
    s32/results/s32_R_dilation_rg_cloud.json   (R-15, EXPLORATORY)

WHY THIS FILE EXISTS, STATED PLAINLY.  Both artefacts were first produced by an inline shell
heredoc and had **no committed producer**.  Lane V's AUDIT 11 caught it by asking, for every
`s32/results/*.json`, whether any `s32/*.py` writes that filename -- the sharper test,
because it catches a missing PRODUCER rather than a missing key.  Two RESULT-grade endpoint
numbers and a registered falsification of P1.3 were resting on artefacts nobody could
re-run.  This is the third instance of `findings-prose-is-not-evidence-of-code` in this
project and the first caught by machine.  The computation is reconstructed here EXACTLY as
it was run -- same clouds, same operator, same grid, same statistics -- and re-emitted, so
the artefacts have a producer and the numbers can be checked.

BASIS: **CA POINT CLOUD, n = 126, tuning126.** An INTERMEDIATE.  These numbers are never
differenced against a chain number.  The chain-basis versions of the same arms live in
`s32_R_repair.py` -> `s32_R_repair.json` (`SCALE_NF`, `SCALE_NF_MED`, `SCALE_GRID`).

THE TWO CALIBRATIONS, and why a scalar was tried at all.  The production average's virtual
CA-CA bond is far short of the builder's ideal, so the cloud is off the valid-chain manifold
partly by being CONTRACTED.  Two native-free scalars can repair that, and they disagree
about how much repair is needed, which is the point:

    s_nf  = ideal_virtual_bond / mean_virtual_bond(C)        -- fixes |i-j| = 1
    s_rg  = mean Rg of the 75 members / Rg(C)                -- fixes the global size

Both dilate about the cloud's own centroid, so each is a pure dilation of shape.  Neither
reads the native or any label.  The ORACLE grid is scored on native RMSD and is labelled
ORACLE / NOT DEPLOYABLE; it prices the ceiling of the scale axis and nothing else.
"""
from __future__ import annotations

import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s24 import stats_lib as ST                                          # noqa: E402
from s12 import instrument as I                                          # noqa: E402
from core import project as pj                                           # noqa: E402

RESULTS = os.path.join(ROOT, "s32", "results")
S29_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")

#: fixed here, in the committed source.  There is no RNG on this path -- the grid, the
#: clouds and the operator are all deterministic -- so the pinned quantity IS the grid.
GRID = (0.90, 0.94, 0.97, 1.00, 1.03, 1.06, 1.09, 1.12, 1.15, 1.20)
SEED = 32_140_003          # unused by the computation; recorded so the pin is explicit


def ideal_bond():
    """The builder's virtual CA-CA distance, MEASURED from the shipped builder rather than
    quoted.  It is a constant of the ideal geometry."""
    ca = np.asarray(pj.build_ca_exact(np.full((1, 8), -1.0), np.full((1, 8), 2.0)),
                    float)[0]
    return float(np.linalg.norm(np.diff(ca, axis=0), axis=1).mean())


def dilate(C, s):
    m = C.mean(0)
    return m + s * (C - m)


def top75(u, pdb):
    pool = I.pool_idx(u)
    sub = np.asarray(I.shipped_record(pdb)["sub"], int)
    W = np.asarray(u["W"], float)[pool[sub]].astype(np.float32).astype(float)
    return W


def main():
    ideal = ideal_bond()
    grid = np.array(GRID, float)
    rows = []
    for t in I.targets():
        p = t["pdb"]
        with np.load(os.path.join(S29_STRUCTS, f"{p}.npz"), allow_pickle=True) as z:
            C = np.asarray(z["prod"], float)
            fold = int(z["fold"])
        u = I.load_univ(p)
        nat = np.asarray(u["nat_ca"], float)              # ORACLE: scoring only
        W = top75(u, p)
        m = C.mean(0)
        vb = float(np.linalg.norm(np.diff(C, axis=0), axis=1).mean())
        vbn = float(np.linalg.norm(np.diff(nat, axis=0), axis=1).mean())
        rg_C = float(np.sqrt(((C - m) ** 2).sum(-1).mean()))
        rg_mem = float(np.mean([np.sqrt(((w - w.mean(0)) ** 2).sum(-1).mean()) for w in W]))
        rg_nat = float(np.sqrt(((nat - nat.mean(0)) ** 2).sum(-1).mean()))
        s_nf = ideal / vb
        s_rg = rg_mem / max(rg_C, 1e-9)
        eg = [float(I.ca_rmsd(dilate(C, s), nat)) for s in grid]
        rows.append(dict(pdb=p, fold=fold, n=int(t["n"]),
                         vb=vb, vb_nat=vbn, rg_C=rg_C, rg_mem=rg_mem, rg_nat=rg_nat,
                         s_nf=s_nf, s_rg=s_rg,
                         e0=float(I.ca_rmsd(C, nat)),
                         e_nf=float(I.ca_rmsd(dilate(C, s_nf), nat)),
                         e_rg=float(I.ca_rmsd(dilate(C, s_rg), nat)),
                         eg=eg, s_orc=float(grid[int(np.argmin(eg))]),
                         e_orc=float(min(eg))))

    g = lambda k: np.array([r[k] for r in rows], float)                  # noqa: E731
    F = np.array([r["fold"] for r in rows], int)
    names = [r["pdb"] for r in rows]
    E0, ENF, ERG, EO = g("e0"), g("e_nf"), g("e_rg"), g("e_orc")
    VB, RC, RM, RN = g("vb"), g("rg_C"), g("rg_mem"), g("rg_nat")

    common = dict(
        n=len(rows), ideal_bond=ideal, grid=list(GRID), seed=SEED,
        basis="CA POINT CLOUD, tuning126 -- an INTERMEDIATE; never differenced "
              "against a chain number",
        geometry=dict(
            cloud_vbond_mean=float(VB.mean()), cloud_vbond_sd=float(VB.std()),
            native_vbond_mean=float(g("vb_nat").mean()),
            contraction_at_sep1=float(1.0 - VB.mean() / ideal),
            cloud_rg_mean=float(RC.mean()), members_rg_mean=float(RM.mean()),
            native_rg_mean=float(RN.mean()),
            contraction_in_rg=float(1.0 - RC.mean() / RM.mean()),
            note="the contraction is SEPARATION-DEPENDENT: these two numbers are the same "
                 "distortion measured at |i-j|=1 and at the radius of gyration, and "
                 "neither may be substituted for the other"),
        cloud_mean=float(E0.mean()),
        oracle_grid_mean=float(EO.mean()),
        oracle_note="ORACLE / NOT DEPLOYABLE -- the per-target best factor is chosen on "
                    "native RMSD and is a best-of-%d order statistic" % len(GRID),
        rows=rows)

    nf = dict(common, arm="s_nf = ideal_virtual_bond / mean_virtual_bond(C)",
              registered="PREREG_S32_R.md P1.3 (R-14)",
              s_mean=float(g("s_nf").mean()), s_min=float(g("s_nf").min()),
              s_max=float(g("s_nf").max()), arm_mean=float(ENF.mean()),
              cmp=ST.compare(ENF, E0, F, names=names,
                             label="dilated-to-ideal-bond cloud - production cloud (CLOUD)"),
              provenance=ST.provenance(__file__))
    rg = dict(common, arm="s_rg = mean Rg of the 75 members / Rg(C)",
              registered="EXPLORATORY (R-15)",
              s_mean=float(g("s_rg").mean()), s_min=float(g("s_rg").min()),
              s_max=float(g("s_rg").max()), arm_mean=float(ERG.mean()),
              cmp=ST.compare(ERG, E0, F, names=names,
                             label="Rg-matched dilated cloud - production cloud (CLOUD)"),
              provenance=ST.provenance(__file__))

    for obj, path in ((nf, "s32_R_dilation_cloud.json"),
                      (rg, "s32_R_dilation_rg_cloud.json")):
        p = os.path.join(RESULTS, path)
        tmp = p + ".%d.tmp" % os.getpid()          # a SHARED .tmp is not atomic
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=1)
        os.replace(tmp, p)
        print("wrote", p)

    print("\nideal virtual CA-CA from the shipped builder: %.4f A" % ideal)
    print("production cloud vbond %.4f  -> contraction at |i-j|=1  %.2f%%"
          % (VB.mean(), 100 * common["geometry"]["contraction_at_sep1"]))
    print("production cloud Rg    %.4f  vs its own members %.4f  -> contraction in Rg %.2f%%"
          % (RC.mean(), RM.mean(), 100 * common["geometry"]["contraction_in_rg"]))
    print("native Rg %.4f, native vbond %.4f\n" % (RN.mean(), g("vb_nat").mean()))
    print("CLOUD BASIS (an intermediate, NOT the endpoint), n = %d:" % len(rows))
    print("  production cloud                       %.4f" % E0.mean())
    print("  R-14 dilated to ideal bond (s %.4f)   %.4f  eff %+.4f  %.2fxMDE  %d/5 folds  %s"
          % (nf["s_mean"], ENF.mean(), nf["cmp"]["effect"],
             abs(nf["cmp"]["effect_over_mde"]), nf["cmp"]["folds_same_sign"],
             nf["cmp"]["verdict"][:18]))
    print("  R-15 Rg-matched            (s %.4f)   %.4f  eff %+.4f  %.2fxMDE  %d/5 folds  %s"
          % (rg["s_mean"], ERG.mean(), rg["cmp"]["effect"],
             abs(rg["cmp"]["effect_over_mde"]), rg["cmp"]["folds_same_sign"],
             rg["cmp"]["verdict"][:18]))
    print("  ORACLE best-on-grid                    %.4f   ORACLE / NOT DEPLOYABLE"
          % EO.mean())


if __name__ == "__main__":
    main()
