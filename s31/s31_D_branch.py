# -*- coding: utf-8 -*-
"""S31 lane D -- the lambda=0 branch selection is made on noise.  Is fixing it worth anything?

PRE-REGISTERED IN `s31/PREREG_S31_D_branch.md` BEFORE ANY ARM HERE WAS COMPUTED.  Read that
file first; this module only implements it.

Three arms, all projected IN THE SAME JOB from the SAME stored clouds
(`s29/results/s29_O_structs/<pdb>.npz`, key `prod`), per the D-B rule that both sides of a
built-chain contrast must come from one code path:

  PROD        `core.project.lam_path(pen, (0.0, 0.3), multi=True)` replicated exactly -- the
              incumbent.  Asserted per target against `s29_O_chain_rows.jsonl :: prod`.
  B4          continue EACH of the four lam=0 branch solutions to lam=0.3, union the four
              fresh lam=0.3 starts, argmin of the LAM=0.3 OBJECTIVE over all eight.
              Native-free, deterministic, DEPLOYABLE.
  ORACLE_B4   the same eight candidates, argmin of CA-RMSD TO THE NATIVE.
              **ORACLE / NOT DEPLOYABLE** -- a ceiling, never a method.

REGISTERED BEFORE THE RESULT (prereg section 3): B4's candidate set is a STRICT SUPERSET of
PROD's and both minimise the SAME objective, so B4's objective is <= PROD's on every target
BY CONSTRUCTION.  "B4 reaches a lower objective" is therefore not evidence of anything and
must never be reported as if it were.  The open question is only whether a lower objective is
a better structure -- and this project's standing finding is that the objective does not rank
the native, so the registered expectation is H0.

Usage:
    python s31/s31_D_branch.py run [--shard i --nshards k] [--pdbs A,B] [--limit N]
    python s31/s31_D_branch.py report
"""
import argparse
import glob
import json
import math
import os
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
OUT = os.path.join(RESULTS, "s31_D_branch.json")
ROWS = os.path.join(RESULTS, "s31_D_branch_rows.jsonl")

LAM, MAXITER, GRAD, PENALTY = 0.3, 300, "exact", "ramah"


def rows_path(shard=None, nshards=0):
    if shard is None:
        return ROWS
    return ROWS.replace(".jsonl", "_shard%dof%d.jsonl" % (int(shard), int(nshards)))


def generic_starts(n):
    return [(np.full(n, math.radians(a)), np.full(n, math.radians(b))) for a, b in pj.STARTS]


def one_target(pdb):
    with np.load(os.path.join(STRUCTS, "%s.npz" % pdb)) as z:
        C = np.array(z["prod"], float)
        seq, fold = str(z["seq"]), int(z["fold"])
    nat = np.asarray(I.load_univ(pdb)["nat_ca"], float)
    n = len(C)
    pen = pj.make_penalty(PENALTY, seq, int(fold))
    starts = generic_starts(n)
    t0 = time.time()

    # ---- rung lam = 0: the four generic starts.  This is where the noise-level argmin is.
    sol0 = [pj.fit_prior(C, ph, ps, pen=None, lam=0.0, maxiter=MAXITER, grad=GRAD)
            for ph, ps in starts]
    f0 = np.array([s[3] for s in sol0], float)
    i0 = int(np.argmin(f0))                       # strict argmin, first wins -- as fit_multi
    o0 = np.sort(f0)
    margin_lam0 = float(o0[1] - o0[0])

    # ---- rung lam = LAM.
    # the four FRESH multi-starts (production computes these too, as `alt`)
    fresh = [pj.fit_prior(C, ph, ps, pen=pen, lam=LAM, maxiter=MAXITER, grad=GRAD)
             for ph, ps in starts]
    # the continuation from EACH lam=0 branch (production computes only the one from i0)
    cont = [pj.fit_prior(C, s[1], s[2], pen=pen, lam=LAM, maxiter=MAXITER, grad=GRAD)
            for s in sol0]

    # ---- PROD, replicating lam_path(multi=True) exactly
    got = cont[i0]
    alt = min(fresh, key=lambda r: r[3])
    prod = alt if alt[3] < got[3] else got

    # ---- B4: argmin of the LAM objective over the union of all eight
    cands = cont + fresh
    fc = np.array([c[3] for c in cands], float)
    jb = int(np.argmin(fc))
    b4 = cands[jb]
    fs = np.sort(fc)
    margin_b4 = float(fs[1] - fs[0])

    # ---- ORACLE_B4: argmin of RMSD to native over the same eight.  NOT DEPLOYABLE.
    rr = np.array([float(I.ca_rmsd(np.asarray(c[0], float), nat)) for c in cands], float)
    jo = int(np.argmin(rr))

    r_prod = float(I.ca_rmsd(np.asarray(prod[0], float), nat))
    r_b4 = float(I.ca_rmsd(np.asarray(b4[0], float), nat))
    r_or = float(rr[jo])

    return dict(
        pdb=pdb, n=n, fold=fold, kind="branch",
        rmsd_cloud=float(I.ca_rmsd(C, nat)),
        rmsd_PROD=r_prod, rmsd_B4=r_b4, rmsd_ORACLE_B4=r_or,
        obj_PROD=float(prod[3]), obj_B4=float(b4[3]),
        # the conditioning quantities the whole exercise is about
        margin_lam0=margin_lam0, margin_B4=margin_b4,
        lam0_objs=[float(x) for x in f0], lam_objs=[float(x) for x in fc],
        lam0_argmin=i0, b4_argmin=jb, oracle_argmin=jo,
        b4_is_prod=bool(abs(r_b4 - r_prod) < 1e-9),
        oracle_is_prod=bool(abs(r_or - r_prod) < 1e-9),
        secs=float(time.time() - t0))


def load_rows(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except Exception:
            continue
        if r.get("kind") == "branch":
            out[r["pdb"]] = r
    return out


def all_rows():
    out = {}
    for f in sorted(glob.glob(ROWS.replace(".jsonl", "*.jsonl"))):
        for pdb, r in load_rows(f).items():
            if pdb in out and abs(out[pdb]["rmsd_PROD"] - r["rmsd_PROD"]) > 0.0:
                raise AssertionError("%s computed twice and disagrees" % pdb)
            out[pdb] = r
    return out


def cmd_run(a):
    os.makedirs(RESULTS, exist_ok=True)
    pdbs = ([p.strip() for p in a.pdbs.split(",")] if a.pdbs else
            sorted(x[:-4] for x in os.listdir(STRUCTS) if x.endswith(".npz")))
    if a.limit:
        pdbs = pdbs[:int(a.limit)]
    if a.nshards:
        pdbs = [p for k, p in enumerate(pdbs) if k % int(a.nshards) == int(a.shard)]
    done = all_rows()
    fh = open(rows_path(a.shard if a.nshards else None, a.nshards), "a", encoding="utf-8")
    for i, pdb in enumerate(pdbs):
        if pdb in done and not a.force:
            continue
        r = one_target(pdb)
        fh.write(json.dumps(r) + "\n")
        fh.flush()
        print("[%3d/%3d] %-6s n=%3d  PROD %.5f  B4 %.5f (%+.4f)  ORACLE %.5f (%+.4f)  "
              "margin lam0 %.2e -> B4 %.2e  %.1fs"
              % (i + 1, len(pdbs), pdb, r["n"], r["rmsd_PROD"], r["rmsd_B4"],
                 r["rmsd_B4"] - r["rmsd_PROD"], r["rmsd_ORACLE_B4"],
                 r["rmsd_ORACLE_B4"] - r["rmsd_PROD"], r["margin_lam0"], r["margin_B4"],
                 r["secs"]), flush=True)
    fh.close()
    if not a.nshards:
        cmd_report(a)


def cmd_report(a):
    rows = all_rows()
    if not rows:
        print("no rows")
        return
    pdbs = sorted(rows)
    folds = np.array([rows[p]["fold"] for p in pdbs])
    prod = np.array([rows[p]["rmsd_PROD"] for p in pdbs])
    b4 = np.array([rows[p]["rmsd_B4"] for p in pdbs])
    orc = np.array([rows[p]["rmsd_ORACLE_B4"] for p in pdbs])

    # ---- the identity check: PROD must reproduce the canonical rows target by target
    ref = {}
    for line in open(S29_ROWS, encoding="utf-8"):
        line = line.strip()
        if line:
            r = json.loads(line)
            if r.get("item") == "prod":
                ref[r["pdb"]] = r["rmsd_chain"]
    common = [p for p in pdbs if p in ref]
    dref = np.array([rows[p]["rmsd_PROD"] - ref[p] for p in common])

    out = {
        "what": ("S31-D: carrying all four lam=0 branches to lam=0.3 instead of selecting one "
                 "at 1e-7 noise.  Pre-registered in s31/PREREG_S31_D_branch.md."),
        "prereg": "s31/PREREG_S31_D_branch.md",
        "provenance": ST.provenance(__file__),
        "n": len(pdbs),
        "basis": "built chain, mean CA-RMSD to native, 126 dev targets",
        "clouds": "s29/results/s29_O_structs/<pdb>.npz :: prod (all arms, same job)",
        "identity_check_vs_s29_canonical": {
            "n": len(common),
            "max_abs_diff": float(np.abs(dref).max()) if len(common) else None,
            "mean_diff": float(dref.mean()) if len(common) else None,
            "PASS": bool(len(common) and np.abs(dref).max() < 1e-9),
            "note": ("PROD must reproduce s29_O_chain_rows.jsonl::prod bit-for-bit; if it "
                     "does not, the arms are not comparable to the canonical endpoint"),
        },
        "means": {"PROD": float(prod.mean()), "B4": float(b4.mean()),
                  "ORACLE_B4": float(orc.mean())},
        "objective_domination_check": {
            "n_B4_objective_le_PROD": int(sum(
                1 for p in pdbs if rows[p]["obj_B4"] <= rows[p]["obj_PROD"] + 1e-12)),
            "note": ("BY CONSTRUCTION this must be n.  It is registered as a structural fact, "
                     "NOT as evidence for the hypothesis -- see prereg section 3."),
        },
    }

    # ---- ORACLE CEILING FIRST (prereg section 4)
    out["ORACLE_B4_vs_PROD"] = ST.compare(
        orc, prod, folds, names=pdbs,
        label="ORACLE_B4 - PROD (built chain) -- ORACLE / NOT DEPLOYABLE, a ceiling",
        seed_parts=("s31D", "branch", "oracle"))
    # ---- the primary
    out["B4_vs_PROD"] = ST.compare(
        b4, prod, folds, names=pdbs,
        label="B4 - PROD (built chain) -- native-free, DEPLOYABLE, the PRIMARY",
        seed_parts=("s31D", "branch", "b4"))

    def gate(c):
        x = abs(c["effect"]) / c["mde"] if c["mde"] else float("nan")
        if x < 0.7:
            return "NOT A RESULT (<0.7x MDE)"
        if x < 1.0:
            return "NOT MEASURED (0.7-1.0x MDE)"
        return "BETTER (>=1.0x MDE)" if c["effect"] < 0 else "WORSE (>=1.0x MDE)"

    out["gates"] = {"ORACLE_B4_vs_PROD": gate(out["ORACLE_B4_vs_PROD"]),
                    "B4_vs_PROD": gate(out["B4_vs_PROD"])}
    out["gate_note"] = ("QUOTE THE GATE, NOT .verdict.  s24.stats_lib.compare is "
                        "LOWER-IS-BETTER (d = a - b, negative = a better); these arms are "
                        "RMSD so the convention is native here, but the rule stands.")

    # ---- conditioning: the whole point
    m0 = np.array([rows[p]["margin_lam0"] for p in pdbs])
    mb = np.array([rows[p]["margin_B4"] for p in pdbs])
    out["conditioning"] = {
        "what": ("the margin by which the winning candidate beat the runner-up.  PROD decides "
                 "at lam=0; B4 decides at lam=0.3."),
        "median_margin_PROD_lam0": float(np.median(m0)),
        "median_margin_B4_lam": float(np.median(mb)),
        "n_PROD_margin_below_1e-6": int((m0 < 1e-6).sum()),
        "n_B4_margin_below_1e-6": int((mb < 1e-6).sum()),
        "improvement_factor_median": float(np.median(mb) / max(float(np.median(m0)), 1e-300)),
    }

    # ---- the tail
    d = b4 - prod
    changed = [p for p in pdbs if not rows[p]["b4_is_prod"]]
    out["tail"] = {
        "n_targets_branch_changed": len(changed),
        "changed_targets": changed[:40],
        "median_effect": float(np.median(d)), "mean_effect": float(d.mean()),
        "max_gain": float(d.min()), "max_loss": float(d.max()),
        "effect_among_changed": float(d[[pdbs.index(p) for p in changed]].mean())
                                if changed else 0.0,
        "concentration": ST.concentration(d, seed_parts=("s31D", "branch", "conc")),
    }

    # ---- FAIL18 vs the 108
    try:
        f18 = set(I.FAIL18)
    except Exception:
        f18 = set()
    if f18:
        m = np.array([p in f18 for p in pdbs])
        if m.sum() and (~m).sum():
            out["FAIL18"] = ST.compare(b4[m], prod[m], folds[m],
                                       names=[p for p in pdbs if p in f18],
                                       label="B4 - PROD on FAIL18",
                                       seed_parts=("s31D", "branch", "f18"))
            out["OTHER108"] = ST.compare(b4[~m], prod[~m], folds[~m],
                                         names=[p for p in pdbs if p not in f18],
                                         label="B4 - PROD on the other 108",
                                         seed_parts=("s31D", "branch", "o108"))

    ST.save_atomic(OUT, out, module_file=__file__)
    for k in ("n", "basis", "means", "gates", "conditioning", "tail",
              "identity_check_vs_s29_canonical", "objective_domination_check"):
        print(k, "=", json.dumps(out[k], default=str)[:700])
    for k in ("ORACLE_B4_vs_PROD", "B4_vs_PROD"):
        c = out[k]
        print("%-22s effect %+.4f  SE %.4f  MDE %.4f  %.2fx  foldCI %s  %dW/%dL  median %+.4f  -> %s"
              % (k, c["effect"], c["se"], c["mde"], abs(c["effect"]) / c["mde"],
                 [round(x, 4) for x in (c.get("ci95_fold") or [])],
                 c["n_better"], c["n_worse"], c["median_effect"], out["gates"][k]))
    print("wrote", OUT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "report"])
    ap.add_argument("--pdbs", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[a.mode](a)


if __name__ == "__main__":
    main()
