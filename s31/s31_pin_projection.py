# -*- coding: utf-8 -*-
"""S31 defect D-B -- pin the cloud -> built-chain projection and prove it deterministic.

THE DEFECT AS REPORTED.  Identical CA point clouds were observed to give different built
chains across records:

    s29/results/s29_O_chain_rows.jsonl :: item=prod   chain 3.2105   cloud 3.0483  <- canonical
    s27/results/chain_rows.jsonl       :: config=DIS  chain 3.2126   cloud 3.0483
    independent re-projections                        chain 3.2071, 3.2126

and the cause was recorded as "an unpinned multi-start projection seed"
(BRIEF section 19/20B; S28-L18 at `s27/LEDGER.md:1025`, scope-corrected in S28-L27b at
`s27/LEDGER.md:1872`).

WHAT THIS SCRIPT ESTABLISHES.  There is NO RANDOM NUMBER GENERATOR anywhere on the
projection path.  `core.project.fit_multi` loops over the four FIXED generic starts in
`core.project.STARTS` and keeps the strict argmin; `fit_prior` is a deterministic
L-BFGS-B call; the builder and the Kabsch are deterministic.  So there is no seed to set.
The projection is a deterministic FUNCTION of the input cloud's bits -- and that is the
whole defect: it is a DISCONTINUOUS function of them.  The multi-start argmin selects
between near-degenerate torsion branches, so an input cloud that differs in the 14th
decimal (which is what "the same cloud" means across two code paths that build the
average differently) can select the other branch, and the two branches are 0.1-0.5 A
apart on a handful of targets.

The pin therefore has two halves and both are recorded here:

  1. THE OPERATOR PIN -- every constant that determines the output, digested, so that a
     future run can assert it is the same operator.  `core.project.PROJECTION_PIN`.
  2. THE INPUT PIN -- the sha256 of each canonical production cloud's float64 bytes, so
     that a future reprojection can assert it is being handed the same bits.  Without
     this half the operator pin is worthless, because the input is where the variation
     actually entered.

and the verification the brief asks for -- reproject the same cloud twice, get
bit-identical chains -- is `mode=determinism`, run at n=126.

It also emits the FLIP REGISTER: for every target, the margin by which the winning
branch won, at both rungs of the lambda path.  A small margin is a target whose chain
value is not stable under a 1e-14 input perturbation, and the register says how much of
the 126-target mean sits on such targets.  That is the instrument's own noise floor, and
it is what any future sub-0.01 A chain claim has to clear.

NOTHING HERE CHANGES ANY PRODUCED NUMBER.  No default is altered, no history is
retro-fitted.  The canonical endpoint remains 3.2105 A (see `canonical` in the output).

Usage:
    python s31/s31_pin_projection.py determinism [--limit N] [--pdbs A,B]
    python s31/s31_pin_projection.py report
"""
import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s24 import stats_lib as ST          # noqa: E402
from core import project as pj           # noqa: E402

RESULTS = os.path.join(ROOT, "s31", "results")
STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")
S29_ROWS = os.path.join(ROOT, "s29", "results", "s29_O_chain_rows.jsonl")
S27_ROWS = os.path.join(ROOT, "s27", "results", "chain_rows.jsonl")
OUT = os.path.join(RESULTS, "s31_D_projection_pin.json")
ROWS = os.path.join(RESULTS, "s31_D_projection_pin_rows.jsonl")

#: the production projection's configuration, as `s12.instrument.project` calls it.
LAM = 0.3
MAXITER = 300
GRAD = "exact"
PENALTY = "ramah"


# --------------------------------------------------------------------------- the pin
def operator_pin():
    """Every constant that determines the projection's output, plus a digest of them.

    If this digest changes, the operator changed, and no reprojection is comparable to a
    number produced before the change.  This is the replacement for the "seed" the brief
    asks to pin: there is no seed, there is this.
    """
    d = {
        "starts_deg": [list(s) for s in pj.STARTS],
        "n_starts": len(pj.STARTS),
        "grad_mode": GRAD,
        "grad_default": pj.GRAD,
        "fd_eps": float(pj.FD_EPS),
        "lam": LAM,
        "maxiter": MAXITER,
        "multi_start": True,
        "penalty": PENALTY,
        "tie_break": "strict argmin, first start wins an exact tie (fit_multi)",
        "rng": "NONE -- no random number generator is constructed on this path",
        "optimiser": "scipy.optimize.minimize L-BFGS-B, jac=True",
    }
    d["digest"] = hashlib.sha256(
        json.dumps(d, sort_keys=True).encode("utf-8")).hexdigest()[:16]
    return d


def cloud_sha(C):
    """sha256 of the cloud's float64 bytes, C-contiguous.  The INPUT half of the pin."""
    return hashlib.sha256(np.ascontiguousarray(np.asarray(C, float)).tobytes()).hexdigest()


def env_pin():
    import scipy
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "blas_threads_env": {k: os.environ.get(k)
                             for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                       "OPENBLAS_NUM_THREADS")},
    }


# ------------------------------------------------------------------- the flip register
def branch_margins(C, seq, fold):
    """How near-tied was the decision the multi-start argmin made, at both rungs?

    Rung lam=0: the four generic starts' objectives.  The margin is best minus runner-up.
    Rung lam=LAM: `lam_path(multi=True)` chooses between the CONTINUATION from the lam=0
    solution and a fresh multi-start.  Both are computed here and both margins reported.
    A margin at the scale of the input perturbation (1e-14) is a target whose emitted
    chain is a coin toss between branches that may be 0.5 A apart.
    """
    C = np.asarray(C, float)
    n = len(C)
    pen = pj.make_penalty(PENALTY, seq, int(fold))
    starts = [(np.full(n, math.radians(a)), np.full(n, math.radians(b)))
              for a, b in pj.STARTS]

    o0 = []
    best0 = None
    for ph, ps in starts:
        got = pj.fit_prior(C, ph, ps, pen=None, lam=0.0, maxiter=MAXITER, grad=GRAD)
        o0.append(float(got[3]))
        if best0 is None or got[3] < best0[3]:
            best0 = got
    s0 = sorted(o0)
    margin0 = float(s0[1] - s0[0])

    cont = pj.fit_prior(C, best0[1], best0[2], pen=pen, lam=LAM, maxiter=MAXITER, grad=GRAD)
    oL = []
    bestL = None
    for ph, ps in starts:
        got = pj.fit_prior(C, ph, ps, pen=pen, lam=LAM, maxiter=MAXITER, grad=GRAD)
        oL.append(float(got[3]))
        if bestL is None or got[3] < bestL[3]:
            bestL = got
    sL = sorted(oL)
    # what lam_path(multi=True) actually decides: continuation vs fresh multi-start
    chose_alt = bool(bestL[3] < cont[3])
    margin_final = float(abs(cont[3] - bestL[3]))
    return {
        "obj_lam0_starts": o0, "margin_lam0": margin0,
        "obj_lam_starts": oL, "margin_lam_starts": float(sL[1] - sL[0]),
        "obj_continuation": float(cont[3]), "obj_multistart": float(bestL[3]),
        "chose_multistart_over_continuation": chose_alt,
        "margin_final_decision": margin_final,
    }


# --------------------------------------------------------------------------- the run
def load_rows(path, key, val):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get(key) == val:
                out[r["pdb"]] = r
    return out


def one_target(pdb):
    f = os.path.join(STRUCTS, "%s.npz" % pdb)
    with np.load(f) as z:
        C = np.array(z["prod"], float)
        seq = str(z["seq"])
        fold = int(z["fold"])
    u = I.load_univ(pdb)
    nat = np.asarray(u["nat_ca"], float)

    t0 = time.time()
    a = I.project(C, seq, fold, lam=LAM, multi=True, maxiter=MAXITER)
    b = I.project(C, seq, fold, lam=LAM, multi=True, maxiter=MAXITER)
    ca_a, ca_b = np.asarray(a["ca"], float), np.asarray(b["ca"], float)

    bitwise = bool(ca_a.tobytes() == ca_b.tobytes())
    coord_maxdiff = float(np.abs(ca_a - ca_b).max())
    r_a = float(I.ca_rmsd(ca_a, nat))
    r_b = float(I.ca_rmsd(ca_b, nat))

    row = dict(pdb=pdb, n=int(len(C)), fold=fold,
               cloud_sha256=cloud_sha(C),
               rmsd_cloud=float(I.ca_rmsd(C, nat)),
               rmsd_chain_pass1=r_a, rmsd_chain_pass2=r_b,
               chain_bit_identical=bitwise,
               coord_maxdiff=coord_maxdiff,
               rmsd_diff=float(r_a - r_b),
               secs=float(time.time() - t0))
    row.update(branch_margins(C, seq, fold))
    return row


def cmd_determinism(args):
    os.makedirs(RESULTS, exist_ok=True)
    pdbs = ([p.strip() for p in args.pdbs.split(",")] if args.pdbs else
            sorted(x[:-4] for x in os.listdir(STRUCTS) if x.endswith(".npz")))
    if args.limit:
        pdbs = pdbs[:int(args.limit)]
    done = load_rows(ROWS, "kind", "pin")
    fh = open(ROWS, "a", encoding="utf-8")
    for i, pdb in enumerate(pdbs):
        if pdb in done and not args.force:
            continue
        row = one_target(pdb)
        row["kind"] = "pin"
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        print("[%3d/%3d] %-6s n=%3d  bit-identical %-5s  |dCA|max %.2e  chain %.6f  "
              "margin(final) %.3e  margin(lam0) %.3e  %.1fs"
              % (i + 1, len(pdbs), pdb, row["n"], row["chain_bit_identical"],
                 row["coord_maxdiff"], row["rmsd_chain_pass1"],
                 row["margin_final_decision"], row["margin_lam0"], row["secs"]),
              flush=True)
    fh.close()
    return cmd_report(args)


def cmd_report(args):
    rows = load_rows(ROWS, "kind", "pin")
    if not rows:
        print("no rows at %s" % ROWS)
        return
    pdbs = sorted(rows)
    s29 = load_rows(S29_ROWS, "item", "prod")
    s27 = load_rows(S27_ROWS, "config", "DIS")

    r1 = np.array([rows[p]["rmsd_chain_pass1"] for p in pdbs])
    r2 = np.array([rows[p]["rmsd_chain_pass2"] for p in pdbs])
    cl = np.array([rows[p]["rmsd_cloud"] for p in pdbs])
    allbit = all(rows[p]["chain_bit_identical"] for p in pdbs)
    maxcoord = max(rows[p]["coord_maxdiff"] for p in pdbs)

    out = {
        "what": "S31 D-B: the cloud->built-chain projection, pinned and proven deterministic",
        "provenance": ST.provenance(__file__),
        "operator_pin": operator_pin(),
        "environment": env_pin(),
        "n": len(pdbs),
        "determinism": {
            "claim": "reprojecting the SAME cloud twice gives BIT-IDENTICAL chains",
            "all_bit_identical": bool(allbit),
            "max_coord_abs_diff_over_targets": float(maxcoord),
            "max_rmsd_abs_diff_over_targets":
                float(np.abs(r1 - r2).max()),
            "n_targets": len(pdbs),
        },
        "reprojection": {
            "mean_chain_pass1": float(r1.mean()),
            "mean_chain_pass2": float(r2.mean()),
            "mean_cloud": float(cl.mean()),
        },
        "canonical": {
            "value_chain": 3.2105,
            "value_cloud": 3.0483,
            "artefact": "s29/results/s29_O_chain_rows.jsonl :: item=prod",
            "why": ("the s29 production ladder is the run the endpoint was declared from "
                    "and the only one whose clouds are persisted target-by-target "
                    "(s29/results/s29_O_structs/*.npz); every other chain mean in the "
                    "record is a reprojection of clouds rebuilt by a different code "
                    "path, agreeing to 1e-14 on the cloud and differing on the chain"),
        },
    }

    # ---- agreement with the canonical rows, and with S27's
    for tag, ref in (("vs_s29_prod_canonical", s29), ("vs_s27_DIS", s27)):
        common = [p for p in pdbs if p in ref]
        if not common:
            continue
        mine = np.array([rows[p]["rmsd_chain_pass1"] for p in common])
        theirs = np.array([ref[p]["rmsd_chain"] for p in common])
        d = mine - theirs
        cd = np.array([rows[p]["rmsd_cloud"] - ref[p].get("rmsd_cloud", np.nan)
                       for p in common])
        k = int(np.argmax(np.abs(d)))
        out[tag] = {
            "n": len(common),
            "mean_mine": float(mine.mean()), "mean_theirs": float(theirs.mean()),
            "mean_diff": float(d.mean()),
            "max_abs_diff": float(np.abs(d).max()),
            "worst_target": common[k], "worst_diff": float(d[k]),
            "n_above_1e-4": int((np.abs(d) > 1e-4).sum()),
            "n_above_0.02": int((np.abs(d) > 0.02).sum()),
            "n_above_0.1": int((np.abs(d) > 0.1).sum()),
            "cloud_max_abs_diff": float(np.nanmax(np.abs(cd))),
        }

    # ---- the flip register
    mf = np.array([rows[p]["margin_final_decision"] for p in pdbs])
    m0 = np.array([rows[p]["margin_lam0"] for p in pdbs])
    risk = [p for p in pdbs if min(rows[p]["margin_final_decision"],
                                   rows[p]["margin_lam0"]) < 1e-6]
    out["flip_register"] = {
        "what": ("the margin by which the multi-start argmin won.  A margin at the scale "
                 "of an input perturbation is a target whose chain value is not stable "
                 "under that perturbation"),
        "median_margin_final": float(np.median(mf)),
        "min_margin_final": float(mf.min()),
        "median_margin_lam0": float(np.median(m0)),
        "min_margin_lam0": float(m0.min()),
        "n_margin_below_1e-6": len(risk),
        "n_margin_below_1e-9": int(sum(1 for p in pdbs
                                       if min(rows[p]["margin_final_decision"],
                                              rows[p]["margin_lam0"]) < 1e-9)),
        "flip_risk_targets": risk[:40],
        "n_chose_multistart_over_continuation":
            int(sum(1 for p in pdbs if rows[p]["chose_multistart_over_continuation"])),
    }

    # ---- the instrument's own noise floor, measured against the confirmed effect
    known = [("s29_O_prod_canonical", 3.2105),
             ("s27_DIS", 3.2126),
             ("s28_A_reprojection", 3.2071),
             ("s31_D_reprojection", float(r1.mean()))]
    vals = [v for _, v in known]
    out["instrument_spread"] = {
        "records": dict(known),
        "spread_A": float(max(vals) - min(vals)),
        "confirmed_effect_for_scale": 0.0221,
        "note": ("any chain claim smaller than the spread is inside the instrument's own "
                 "reprojection noise; this is the number the brief's 'any future sub-0.01 A "
                 "chain claim is invalid' rests on"),
    }

    ST.save_atomic(OUT, out, module_file=__file__)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("provenance", "operator_pin", "environment")},
                     indent=1)[:4000])
    print("\nwrote %s" % OUT)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["determinism", "report"])
    ap.add_argument("--pdbs", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    {"determinism": cmd_determinism, "report": cmd_report}[a.mode](a)


if __name__ == "__main__":
    main()
