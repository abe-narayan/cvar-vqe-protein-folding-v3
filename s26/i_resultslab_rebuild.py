#!/usr/bin/env python
"""s26/i_resultslab_rebuild.py -- the frozen results-lab rebuild, bracketed so it can be judged.

Lane I, Sprint 26 (extended window, coordinator item 3).  The build itself is the documented
command, run unchanged as a subprocess:

    python -m s25.resultslab.build --mode frozen --spec s25/results/real_pools/spec.json

`pre` snapshots the TRACKED `results/summary/*` (copies + sha256) and the sha256 of every
tracked `results/structures/**/*.pdb` (2,142 files) to `s26/results/resultslab_rebuild/`.
`build` runs the command (this is the governed job).  `post` compares the rebuilt tree with the
snapshot: every per-target RMSD on both bases and every leaderboard mean / gate / verdict,
the presence of the L7 WARN explanation in `pool_gate_rule`, the four gates, and the PDBs
(ATOM records compared exactly; REMARK provenance lines expected to differ) -- and writes a
verdict.  If any number moves, `post --restore` puts the tracked bytes back
(`git checkout -- results/summary results/structures`) so nothing unreproduced is left on disk.

    python s26/i_resultslab_rebuild.py pre
    python s26/jobrun.py --agent I --tag CPU --name resultslab_rebuild --est-ram 1.2 \\
        -- python s26/i_resultslab_rebuild.py build
    python s26/i_resultslab_rebuild.py post [--restore]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUMMARY = os.path.join(ROOT, "results", "summary")
STRUCT = os.path.join(ROOT, "results", "structures")
SNAP = os.path.join(ROOT, "s26", "results", "resultslab_rebuild")
SPEC = os.path.join(ROOT, "s25", "results", "real_pools", "spec.json")
PY = sys.executable
WARN_TEXT = "WARN = one or more individual targets below their own pool best"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atom_sha(path):
    """sha256 of the ATOM/HETATM/TER/END records only (provenance REMARKs excluded)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith((b"ATOM", b"HETATM", b"TER", b"END")):
                h.update(line)
    return h.hexdigest()


def pre():
    os.makedirs(os.path.join(SNAP, "summary_before"), exist_ok=True)
    snap = {"taken": time.strftime("%Y-%m-%dT%H:%M:%S"), "summary": {}, "structures": {}}
    for f in sorted(os.listdir(SUMMARY)):
        p = os.path.join(SUMMARY, f)
        shutil.copy2(p, os.path.join(SNAP, "summary_before", f))
        snap["summary"][f] = {"sha256": sha256(p), "bytes": os.path.getsize(p)}
    pdbs = sorted(glob.glob(os.path.join(STRUCT, "**", "*.pdb"), recursive=True))
    for p in pdbs:
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        snap["structures"][rel] = {"sha256": sha256(p), "atom_sha256": atom_sha(p)}
    snap["n_structures"] = len(pdbs)
    tracked = subprocess.run(["git", "ls-files", "results/structures", "results/summary"], cwd=ROOT,
                             capture_output=True, text=True).stdout.split()
    snap["n_tracked"] = len(tracked)
    with open(os.path.join(SNAP, "snapshot_before.json"), "w") as fh:
        json.dump(snap, fh, indent=1)
    print(f"snapshot: {len(snap['summary'])} summary files, {len(pdbs)} PDBs "
          f"({snap['n_tracked']} tracked paths) -> {SNAP}")
    return 0


def build():
    cmd = [PY, "-m", "s25.resultslab.build", "--mode", "frozen", "--spec", SPEC]
    print("running:", " ".join(cmd), flush=True)
    t0 = time.time()
    rc = subprocess.run(cmd, cwd=ROOT).returncode
    print(f"build exit {rc} after {time.time() - t0:.0f} s", flush=True)
    return rc


def _rows(results_json):
    with open(results_json, encoding="utf-8") as fh:
        d = json.load(fh)
    return d, {(r["configuration"], r["pdb_id"]): r for r in d["records"]}


def post(restore=False):
    with open(os.path.join(SNAP, "snapshot_before.json")) as fh:
        snap = json.load(fh)
    before_dir = os.path.join(SNAP, "summary_before")
    B, brows = _rows(os.path.join(before_dir, "results.json"))
    A, arows = _rows(os.path.join(SUMMARY, "results.json"))
    out = {"checked": time.strftime("%Y-%m-%dT%H:%M:%S"), "problems": [], "notes": []}

    # 1. per-target RMSDs, both bases, exact equality
    keys = sorted(set(brows) | set(arows))
    n_exact = n_close = n_diff = 0
    worst = 0.0
    diffs = []
    for k in keys:
        b, a = brows.get(k), arows.get(k)
        if b is None or a is None:
            out["problems"].append(f"record {k} present on one side only")
            continue
        for fld in ("rmsd", "rmsd_secondary"):
            x, y = b.get(fld), a.get(fld)
            if x is None and y is None:
                continue
            if x is None or y is None:
                out["problems"].append(f"{k}.{fld}: {x!r} vs {y!r}")
                continue
            if x == y:
                n_exact += 1
            elif abs(x - y) < 1e-9:
                n_close += 1
                worst = max(worst, abs(x - y))
            else:
                n_diff += 1
                worst = max(worst, abs(x - y))
                if len(diffs) < 30:
                    diffs.append({"key": list(k), "field": fld, "tracked": x, "rebuilt": y})
    out["per_target"] = {"n_records": len(keys), "n_values_exact": n_exact, "n_values_close_1e-9": n_close,
                         "n_values_different": n_diff, "worst_abs_diff": worst, "first_diffs": diffs}
    if n_diff:
        out["problems"].append(f"{n_diff} per-target RMSD values differ (worst {worst:.3e} A)")

    # 2. the leaderboard
    def lb(d):
        return {r["configuration"]: r for r in d["leaderboard"]}
    LB_B, LB_A = lb(B), lb(A)
    lbcmp = {}
    for cfg in LB_B:
        b, a = LB_B[cfg], LB_A.get(cfg)
        if a is None:
            out["problems"].append(f"leaderboard row {cfg} missing after rebuild")
            continue
        row = {}
        for fld in ("mean", "median", "mean_secondary", "pool_gate", "n_pool_violations",
                    "difficulty_gate", "corr_with_pool_best", "paired_effect", "mde", "verdict",
                    "ci95_fold_lo", "ci95_fold_hi", "wins", "losses"):
            x, y = b.get(fld), a.get(fld)
            same = (x == y) or (isinstance(x, float) and isinstance(y, float) and abs(x - y) < 1e-9)
            row[fld] = {"tracked": x, "rebuilt": y, "same": bool(same)}
            if not same:
                out["problems"].append(f"leaderboard {cfg}.{fld}: {x!r} -> {y!r}")
        lbcmp[cfg] = row
    out["leaderboard"] = lbcmp
    prod = LB_A.get("production", {})
    out["headline"] = {"built_chain_mean": prod.get("mean"), "point_cloud_mean": prod.get("mean_secondary"),
                       "pool_gate": prod.get("pool_gate"), "n_pool_violations": prod.get("n_pool_violations"),
                       "difficulty_gate": prod.get("difficulty_gate")}
    # the headline at the precision the record quotes (4 decimals), and exactly against the
    # tracked leaderboard (row comparison above). The lab's point cloud is the PDB-quantised
    # (%8.3f) cloud, 3.048329 in the tracked leaderboard against the cache's 3.048338: the two
    # agree at 4 decimals and differ at the 5th by PDB quantisation (S25 RESULTS 2, 2.4e-4 A worst).
    if not (prod.get("mean") is not None and round(prod["mean"], 4) == 3.2126):
        out["problems"].append(f"production built-chain mean {prod.get('mean')} != 3.2126 at 4 dp")
    if not (prod.get("mean_secondary") is not None and round(prod["mean_secondary"], 4) == 3.0483):
        out["problems"].append(f"production point-cloud mean {prod.get('mean_secondary')} != 3.0483 at 4 dp")

    # 3. the L7 WARN explanation and the gates
    out["pool_gate_rule_has_L7_text"] = WARN_TEXT in str(A.get("pool_gate_rule", ""))
    if not out["pool_gate_rule_has_L7_text"]:
        out["problems"].append("pool_gate_rule does not carry the L7 WARN explanation")
    gates = {}
    for cfg, r in LB_A.items():
        gates[cfg] = {"pool": r.get("pool_gate"), "difficulty": r.get("difficulty_gate")}
        if r.get("pool_gate") == "FAIL" or r.get("difficulty_gate") != "PASS":
            out["problems"].append(f"gate not passed for {cfg}: {gates[cfg]}")
    prov = {r.get("provenance") for r in A["records"]}
    gates["_provenance_values"] = sorted(str(p) for p in prov)
    if prov != {"genuine"}:
        out["problems"].append(f"provenance values after rebuild: {prov}")
    out["gates"] = gates
    out["status"] = A.get("status")
    out["contains_synthetic"] = A.get("contains_synthetic")

    # 4. the PDB files: ATOM records exact, headers expected to differ in provenance lines
    pdbs = sorted(glob.glob(os.path.join(STRUCT, "**", "*.pdb"), recursive=True))
    same_bytes = same_atoms = diff_atoms = new_files = 0
    atom_diffs = []
    for p in pdbs:
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        b = snap["structures"].get(rel)
        if b is None:
            new_files += 1
            continue
        if sha256(p) == b["sha256"]:
            same_bytes += 1
            same_atoms += 1
        elif atom_sha(p) == b["atom_sha256"]:
            same_atoms += 1
        else:
            diff_atoms += 1
            if len(atom_diffs) < 20:
                atom_diffs.append(rel)
    missing = [k for k in snap["structures"] if not os.path.exists(os.path.join(ROOT, k))]
    out["structures"] = {"n_before": snap["n_structures"], "n_after": len(pdbs),
                         "byte_identical": same_bytes, "atom_records_identical": same_atoms,
                         "atom_records_differ": diff_atoms, "new_files": new_files,
                         "missing_after": len(missing), "first_atom_diffs": atom_diffs}
    if diff_atoms or missing or new_files:
        out["problems"].append(f"structures: {diff_atoms} ATOM-record differences, {missing and len(missing)} "
                               f"missing, {new_files} new")

    # 5. byte-wise summary of results/summary
    bw = {}
    for f in sorted(os.listdir(SUMMARY)):
        p = os.path.join(SUMMARY, f)
        b = snap["summary"].get(f)
        bw[f] = {"bytes_before": b and b["bytes"], "bytes_after": os.path.getsize(p),
                 "byte_identical": bool(b and b["sha256"] == sha256(p))}
    out["summary_files"] = bw
    # what differs in results.json apart from numbers: provenance / timestamps / label
    for fld in ("provenance", "label", "pool_gate_rule", "status", "chain_lam"):
        if B.get(fld) != A.get(fld):
            out["notes"].append(f"results.json top-level {fld} differs (expected for provenance/rule)")
    ts_b = {r.get("timestamp") for r in B["records"]}
    ts_a = {r.get("timestamp") for r in A["records"]}
    out["notes"].append(f"record timestamps: tracked {sorted(ts_b)}, rebuilt {sorted(ts_a)}")
    gc_b = {r.get("git_commit") for r in B["records"]}
    gc_a = {r.get("git_commit") for r in A["records"]}
    out["notes"].append(f"record git_commit: tracked {sorted(gc_b)}, rebuilt {sorted(gc_a)}")
    mh_b = {r.get("module_hash") for r in B["records"]}
    mh_a = {r.get("module_hash") for r in A["records"]}
    out["notes"].append(f"record module_hash: tracked {sorted(mh_b)}, rebuilt {sorted(mh_a)}")

    out["verdict"] = "REPRODUCED" if not out["problems"] else "DOES NOT REPRODUCE"
    with open(os.path.join(SNAP, "post_verdict.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k not in ("leaderboard",)}, indent=1, default=str))
    print("VERDICT:", out["verdict"])
    if out["problems"] and restore:
        print("restoring the tracked results/summary and results/structures (git checkout)")
        subprocess.run(["git", "checkout", "--", "results/summary", "results/structures"], cwd=ROOT)
        subprocess.run(["git", "clean", "-fdq", "--", "results/structures"], cwd=ROOT)
    return 0 if not out["problems"] else 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("pre", "build", "post"))
    ap.add_argument("--restore", action="store_true")
    a = ap.parse_args(argv)
    return {"pre": pre, "build": build, "post": lambda: post(a.restore)}[a.cmd]()


if __name__ == "__main__":
    raise SystemExit(main())
