#!/usr/bin/env python
"""s26/i_verify_rerun.py -- re-run the standalone audits under verify/ WITHOUT touching their
tracked JSON, and diff the fresh output against the tracked one.

Lane I, Sprint 26 (extended window, coordinator item 2).  Every audit under `verify/` writes
its JSON beside itself with a hard-coded path (`os.path.join(here, "<name>.json")`, or
`core.project._outpath` for the projection audits).  Those files are tracked evidence and
must not be overwritten, so this runner executes each audit IN THIS PROCESS via `runpy` with
`builtins.open` and `os.replace` / `os.rename` wrapped: any WRITE to a path under `verify/`
(or under `bench_results/`, which `core.project harvest` would touch) is redirected to
`s26/results/verify/<same basename>`; reads are untouched.  The tracked file's sha256 is
recorded before and asserted unchanged after.  The fresh JSON is then compared with the
tracked one, leaf by leaf: counts of identical / numerically-close / different / added /
removed leaves, and the first differences with both values.  Keys that record the run
rather than the result (wall clocks, timestamps, hostnames, commit hashes) are listed
separately so a "differs" verdict is never made of a clock.

    python s26/i_verify_rerun.py list                      # the table below, with cost class
    python s26/i_verify_rerun.py run <audit>               # one audit, in-process, redirected
    python s26/i_verify_rerun.py diff <audit>              # compare fresh vs tracked only
    python s26/i_verify_rerun.py launch [--only a,b] [--skip a,b] [--include-heavy]
                                                           # each audit as its own governed job,
                                                           # sequentially (waits for each)
    python s26/i_verify_rerun.py report                    # s26/results/verify/REPORT.md + .json

Audits whose tracked JSON has no runner in `verify/` (`recon_containment_audit.json`, written
by a deleted script; `project_stability_partial77.json`, a partial run kept by name) are
listed as NOT RE-RUNNABLE.  Audits whose runner writes no tracked JSON (`determinism_audit`,
`hazard_audit`, `amber_affinity`, `amber_threads`) are run and recorded with "no tracked
JSON to diff".  The benchmark is never touched: `headline_audit` reads `s9/final_report.json`
(a cached aggregate, its own docstring says READ ONLY) and nothing here passes
`--i-am-spending-the-benchmark`.
"""
from __future__ import annotations

import argparse
import builtins
import re
import hashlib
import json
import math
import os
import runpy
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VERIFY = os.path.join(ROOT, "verify")
OUT = os.path.join(ROOT, "s26", "results", "verify")
PY = sys.executable

#: name -> (script relative to ROOT, argv, tag, est_ram_gb, tracked json basenames, cost, note)
#: cost: "light" (< 2 min), "medium" (minutes), "heavy" (tens of minutes; --include-heavy)
AUDITS = {
    "vqe_lfo_audit": ("verify/vqe_lfo_audit.py", [], "CPU", 0.2, ["vqe_lfo_audit.json"], "light",
                      "reads s8/integrate_vqe.json; no model, no pipeline"),
    "cvar_audit": ("verify/cvar_audit.py", [], "CPU", 0.6, ["cvar_audit.json"], "light",
                   "core.quantum on synthetic energies; pennylane import"),
    "ansatz_audit": ("verify/ansatz_audit.py", [], "CPU", 0.6, ["ansatz_audit.json"], "light",
                     "core.quantum circuit genuineness; pennylane import"),
    "headline_audit": ("verify/headline_audit.py", [], "CPU", 0.5, ["headline_audit.json"], "light",
                       "recorded constants vs s9/synth_cache and bench_results; READ ONLY on the "
                       "benchmark aggregate s9/final_report.json"),
    "legacy_audit": ("verify/legacy_audit.py", [], "CPU", 0.9, ["legacy_audit.json"], "medium",
                     "real retrieval pool + the 11 Legacy terms, core vs root"),
    "amber_audit": ("verify/amber_audit.py", [], "AMBER", 0.9, ["amber_audit.json"], "medium",
                    "OpenMM System on 1A13; the pinned interaction energy"),
    "amber_platform": ("verify/amber_platform.py", [], "AMBER", 1.0, ["amber_platform.json"], "medium",
                       "CPU vs OpenCL arms of the AMBER block; OpenCL may be unavailable"),
    "leak_audit": ("verify/leak_audit.py", [], "AMBER", 1.6, ["leak_audit.json"], "medium",
                   "NaN-poisons 3 targets through run_target (stage 4 included)"),
    "grad_key_collision": ("verify/grad_key_collision.py", [], "AMBER", 1.8, ["grad_key_collision.json"],
                           "medium", "2 targets x 2 child processes; NOTE: it sets PROJECT_GRAD in the "
                           "environment, which core/project.py no longer reads, so both children now "
                           "run the same mode -- a differing JSON here is the audit being stale, not "
                           "the code (see s26 L19)"),
    "equiv_compare": ("verify/equiv_compare.py", [], "CPU", 0.6, ["equivalence.json"], "light",
                      "per-target comparison of the two arms' cache directories (keys from core)"),
    "projection_divergence": ("verify/projection_divergence.py", ["464a0ddb5f283e04", "1fc9f2dcf489e2fb"],
                              "CPU", 0.6, ["projection_divergence.json"], "light",
                              "arms apart, from the baseline and production cache directories"),
    "project_arms": ("verify/project_arms.py", ["464a0ddb5f283e04", "1fc9f2dcf489e2fb", "1fc9f2dcf489e2fb"],
                     "CPU", 0.6, ["project_arms.json"], "light",
                     "three-arm compare over cache keys; the fd arm key is not on disk any more, "
                     "so the production key is passed twice and that arm reads as identical"),
    "determinism_audit": ("verify/determinism_audit.py",
                          ["--functions", "d3_key_types,d4_forged_collision,d5_config_key,d1_cross_process"],
                          "AMBER", 1.8, [], "medium",
                          "D1 cross-process (3 children on 1CS9, stage 4 included), D3/D4/D5 cache keys; "
                          "D2 (resume) is NOT run: as written it moves the last two records of the "
                          "PRODUCTION cache (bench_results/cache/1fc9f2dcf489e2fb) to a temp dir, "
                          "re-runs smoke8 (which does not contain them) and then rmtree's the backup, "
                          "i.e. it would delete two production records; writes determinism_audit.json, "
                          "which is NOT tracked"),
    "hazard_audit": ("verify/hazard_audit.py", [], "AMBER", 1.6, [], "medium",
                     "the known hazards as live experiments; hazard_audit.json is NOT tracked"),
    "project_selfcheck": ("core/project.py", ["selfcheck"], "CPU", 0.6, [], "light",
                          "builder vs reference, gradient vs finite differences; prints only"),
    "project_exactness": ("core/project.py", ["exactness"], "CPU", 0.8, ["project_exactness.json"],
                          "heavy", "126 targets, shipped mode against s8.project"),
    "project_equiv": ("core/project.py", ["equiv"], "CPU", 0.8, ["project_equiv.json"], "heavy",
                      "126 targets x four arms; the reference arm alone is ~8 s per target"),
    "project_degeneracy": ("core/project.py", ["degeneracy"], "CPU", 0.8, ["project_degeneracy.json"],
                           "heavy", "126 targets, objective gap between the best two starts"),
    "project_iters": ("core/project.py", ["iters"], "CPU", 0.8, ["project_iters.json"], "medium",
                      "24 targets by default (the tracked run's default)"),
    "project_stability": ("core/project.py", ["stability", "77"], "CPU", 0.8,
                          ["project_stability_partial77.json"], "heavy",
                          "the tracked file is a PARTIAL (77-target) run kept under its own name; "
                          "re-run with limit 77 so the comparison is like for like"),
}
NOT_RERUNNABLE = {
    "recon_containment_audit.json": "written by a deleted script (leak_audit_mine.json, moved in "
                                    "the 2026-09-04 consolidation); no runner in verify/",
    "project_inputs.json": "written by `core.project harvest` from the production cache; it is the "
                           "INPUT the projection audits read, re-harvesting it is covered by "
                           "project_exactness / project_equiv reading it",
}
#: leaves whose key names a property of the RUN, not the result
RUN_KEYS = ("wall", "elapsed", "time", "seconds", "_s", "ts", "date", "started", "host",
            "commit", "git", "pid", "speedup", "per_target_s", "cpu")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def _is_redirect_target(path):
    try:
        ap = os.path.abspath(str(path))
    except Exception:                                                    # noqa: BLE001
        return False
    return (ap.startswith(VERIFY + os.sep) or
            ap.startswith(os.path.join(ROOT, "bench_results") + os.sep))


def _redirected(path):
    ap = os.path.abspath(str(path))
    rel = os.path.relpath(ap, ROOT).replace(os.sep, "__")
    return os.path.join(OUT, rel)


class Redirect:
    """Wrap builtins.open / os.replace / os.rename so writes under verify/ or bench_results/
    land under s26/results/verify/ instead.  Reads pass through."""

    def __init__(self):
        self.writes = []

    def __enter__(self):
        import shutil
        self._shutil = shutil
        self._open, self._replace, self._rename = builtins.open, os.replace, os.rename
        self._remove, self._unlink = os.remove, os.unlink
        self._sh_move, self._sh_rmtree = shutil.move, shutil.rmtree
        os.makedirs(OUT, exist_ok=True)
        me = self

        def refuse(what, path):
            raise PermissionError(f"lane I guard: {what} of {path} refused -- it is under verify/ or "
                                  f"bench_results/ (tracked evidence or the production cache)")

        def remove_(path, *a, **k):
            if _is_redirect_target(path) and me._exists(path):
                refuse("delete", path)
            if _is_redirect_target(path):
                return me._remove(_redirected(path), *a, **k)
            return me._remove(path, *a, **k)

        def rmtree_(path, *a, **k):
            if _is_redirect_target(path):
                refuse("rmtree", path)
            return me._sh_rmtree(path, *a, **k)

        def move_(src, dst, *a, **k):
            if _is_redirect_target(src) and me._exists(src):
                refuse("move out", src)
            if _is_redirect_target(dst):
                target = _redirected(dst)
                me.writes.append((os.path.abspath(str(dst)), target))
                src2 = _redirected(src) if _is_redirect_target(src) else src
                return me._sh_move(src2, target, *a, **k)
            return me._sh_move(src, dst, *a, **k)

        def open_(file, mode="r", *a, **k):
            if isinstance(file, (str, os.PathLike)) and any(c in str(mode) for c in "wax+") \
                    and _is_redirect_target(file):
                target = _redirected(file)
                me.writes.append((os.path.abspath(str(file)), target))
                return me._open(target, mode, *a, **k)
            return me._open(file, mode, *a, **k)

        def replace_(src, dst, *a, **k):
            if _is_redirect_target(src) and me._exists(src):
                refuse("move out", src)
            if _is_redirect_target(dst):
                target = _redirected(dst)
                me.writes.append((os.path.abspath(str(dst)), target))
                src2 = _redirected(src) if _is_redirect_target(src) else src
                return me._replace(src2, target, *a, **k)
            return me._replace(src, dst, *a, **k)

        def rename_(src, dst, *a, **k):
            if _is_redirect_target(src) and me._exists(src):
                refuse("move out", src)
            if _is_redirect_target(dst):
                target = _redirected(dst)
                me.writes.append((os.path.abspath(str(dst)), target))
                src2 = _redirected(src) if _is_redirect_target(src) else src
                return me._rename(src2, target, *a, **k)
            return me._rename(src, dst, *a, **k)

        builtins.open, os.replace, os.rename = open_, replace_, rename_
        os.remove, os.unlink = remove_, remove_
        shutil.move, shutil.rmtree = move_, rmtree_
        return self

    def _exists(self, path):
        """True if the ORIGINAL path exists on disk (a redirected write never creates it)."""
        return os.path.exists(os.path.abspath(str(path)))

    def __exit__(self, *exc):
        builtins.open, os.replace, os.rename = self._open, self._replace, self._rename
        os.remove, os.unlink = self._remove, self._unlink
        self._shutil.move, self._shutil.rmtree = self._sh_move, self._sh_rmtree
        return False


def _leaves(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _leaves(v, f"{path}.{k}" if path else str(k))
    elif isinstance(obj, list):
        if obj and all(not isinstance(x, (dict, list)) for x in obj):
            yield path, obj
        else:
            for i, v in enumerate(obj):
                yield from _leaves(v, f"{path}[{i}]")
    else:
        yield path, obj


def _close(a, b, rtol=1e-9, atol=1e-12):
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if math.isnan(a) and math.isnan(b):
            return True
        return abs(a - b) <= atol + rtol * max(abs(a), abs(b))
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return all(_close(x, y, rtol, atol) for x, y in zip(a, b))
    return a == b


def diff_json(tracked_path, fresh_path):
    with open(tracked_path, encoding="utf-8") as fh:
        T = json.load(fh)
    with open(fresh_path, encoding="utf-8") as fh:
        F = json.load(fh)
    t = dict(_leaves(T))
    f = dict(_leaves(F))
    out = {"n_tracked_leaves": len(t), "n_fresh_leaves": len(f), "identical": 0, "close": 0,
           "different": 0, "run_property_different": 0, "added": [], "removed": [], "diffs": [],
           "run_property_diffs": []}
    for k in t:
        if k not in f:
            out["removed"].append(k)
            continue
        a, b = t[k], f[k]
        if a == b:
            out["identical"] += 1
        elif _close(a, b):
            out["close"] += 1
        elif any(s in k.lower() for s in RUN_KEYS):
            out["run_property_different"] += 1
            if len(out["run_property_diffs"]) < 20:
                out["run_property_diffs"].append({"key": k, "tracked": _short(a), "fresh": _short(b)})
        else:
            out["different"] += 1
            if len(out["diffs"]) < 40:
                out["diffs"].append({"key": k, "tracked": _short(a), "fresh": _short(b)})
    out["added"] = [k for k in f if k not in t][:40]
    out["n_added"] = sum(1 for k in f if k not in t)
    out["n_removed"] = len(out["removed"])
    out["removed"] = out["removed"][:40]
    out["verdict"] = ("IDENTICAL" if out["different"] == 0 and out["n_added"] == 0
                      and out["n_removed"] == 0 and out["close"] == 0 else
                      "IDENTICAL up to float noise" if out["different"] == 0 and out["n_added"] == 0
                      and out["n_removed"] == 0 else "DIFFERS")
    return out


def _short(v):
    s = json.dumps(v, default=str)
    return s if len(s) <= 160 else s[:157] + "..."


def run(name):
    script, argv, tag, est, tracked, cost, note = AUDITS[name]
    os.makedirs(OUT, exist_ok=True)
    before = {b: sha256(os.path.join(VERIFY, b)) for b in tracked if os.path.exists(os.path.join(VERIFY, b))}
    t0 = time.time()
    rec = {"audit": name, "script": script, "argv": argv, "tag": tag, "est_ram_gb": est, "cost": cost,
           "note": note, "tracked": tracked, "started": time.strftime("%Y-%m-%dT%H:%M:%S")}
    sys.argv = [os.path.join(ROOT, script)] + list(argv)
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    if VERIFY not in sys.path:
        sys.path.insert(0, VERIFY)
    err = None
    funcs = None
    if argv[:1] == ["--functions"]:
        funcs = argv[1].split(",")
        sys.argv = [os.path.join(ROOT, script)]
    with Redirect() as rd:
        try:
            if funcs:
                mod = runpy.run_path(os.path.join(ROOT, script), run_name="_verify_rerun_import")
                out = {}
                for fname in funcs:
                    try:
                        mod[fname](out)
                    except Exception as exc:                             # noqa: BLE001
                        out[fname + "_ERROR"] = f"{type(exc).__name__}: {str(exc)[:600]}"
                    print(f"  {fname} done", flush=True)
                out["_lane_I_note"] = ("run through s26/i_verify_rerun.py with functions "
                                       + ",".join(funcs) + "; see the AUDITS note for what was skipped")
                here = os.path.join(ROOT, os.path.dirname(script))
                base = os.path.basename(script).replace(".py", ".json")
                with open(os.path.join(here, base), "w") as fh:
                    json.dump(out, fh, indent=2, sort_keys=True, default=str)
            else:
                runpy.run_path(os.path.join(ROOT, script), run_name="__main__")
        except SystemExit as e:
            if e.code not in (None, 0):
                err = f"SystemExit({e.code})"
        except Exception as e:                                           # noqa: BLE001
            import traceback
            err = "".join(traceback.format_exception(type(e), e, e.__traceback__))[-4000:]
    rec["wall_s"] = round(time.time() - t0, 1)
    rec["error"] = err
    rec["writes"] = [{"intended": os.path.relpath(a, ROOT), "written": os.path.relpath(b, ROOT)}
                     for a, b in rd.writes]
    after = {b: sha256(os.path.join(VERIFY, b)) for b in before}
    rec["tracked_unchanged"] = before == after
    assert before == after, f"TRACKED FILE CHANGED during {name}: {before} -> {after}"
    rec["diffs"] = {}
    for b in tracked:
        fresh = _redirected(os.path.join(VERIFY, b))
        tp = os.path.join(VERIFY, b)
        if not os.path.exists(fresh):
            # the runner may name its output by limit / partial count (core.project._outpath):
            # take the freshest written JSON whose basename shares the tracked stem
            stem = re.sub(r"_(partial|limit)\d+$", "", b[:-5])
            cands = [w for i, w in rd.writes
                     if os.path.basename(i).startswith(stem) and i.endswith(".json")
                     and os.path.exists(w)]
            if cands:
                fresh = cands[-1]
        if os.path.exists(fresh) and os.path.exists(tp):
            rec["diffs"][b] = diff_json(tp, fresh)
            rec["diffs"][b]["fresh_path"] = os.path.relpath(fresh, ROOT)
        elif not os.path.exists(fresh):
            rec["diffs"][b] = {"verdict": "NO FRESH OUTPUT", "fresh_path": os.path.relpath(fresh, ROOT)}
    if not tracked:
        rec["diffs"]["_"] = {"verdict": "no tracked JSON to diff"}
    with open(os.path.join(OUT, f"{name}.rerun.json"), "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=1, default=str)
    print(f"\n==== {name}: {'ERROR ' + err.splitlines()[-1] if err else 'ok'} in {rec['wall_s']} s; "
          f"tracked unchanged: {rec['tracked_unchanged']}")
    for b, d in rec["diffs"].items():
        print(f"  {b}: {d.get('verdict')}"
              + (f" (identical {d['identical']}, close {d['close']}, different {d['different']}, "
                 f"run-property {d['run_property_different']}, added {d['n_added']}, removed {d['n_removed']})"
                 if "identical" in d else ""))
        for x in d.get("diffs", [])[:8]:
            print(f"     {x['key']}: tracked {x['tracked']}  fresh {x['fresh']}")
    return 0 if err is None else 1


def launch(only=None, skip=None, include_heavy=False, after=()):
    names = [n for n in AUDITS if (not only or n in only) and (not skip or n not in skip)
             and (include_heavy or AUDITS[n][5] != "heavy")]
    done_dir = os.path.join(ROOT, "s26", "jobs_done")
    for j in after:
        while not os.path.exists(os.path.join(done_dir, f"{j}.json")):
            print(f"launch: waiting for job {j} to finish before the verify chain starts", flush=True)
            time.sleep(30)
    for n in names:
        script, argv, tag, est, tracked, cost, note = AUDITS[n]
        cmd = [PY, os.path.join(ROOT, "s26", "jobrun.py"), "--agent", "I", "--tag", tag,
               "--name", f"verify_{n}", "--est-ram", str(est), "--", PY,
               os.path.join(ROOT, "s26", "i_verify_rerun.py"), "run", n]
        print(f"launching verify_{n} (tag {tag}, est {est} GB, {cost})", flush=True)
        rc = subprocess.run(cmd, cwd=ROOT).returncode
        print(f"verify_{n} -> exit {rc}", flush=True)
    return 0


def report():
    rows = []
    for n in AUDITS:
        p = os.path.join(OUT, f"{n}.rerun.json")
        if os.path.exists(p):
            with open(p, encoding="utf-8") as fh:
                rows.append(json.load(fh))
    L = ["# verify/ audits re-run under the S26 governor (lane I)", "",
         "Each audit was executed in-process by `s26/i_verify_rerun.py run <name>` with its writes "
         "redirected to `s26/results/verify/`; the tracked JSON's sha256 was asserted unchanged. "
         "Leaves compared: identical / close (float noise) / different / run-property (clocks, "
         "hosts, commits) / added / removed.", "",
         "| audit | tag | peak RSS | wall | tracked JSON | verdict | identical | close | different | run-prop | added | removed | error |",
         "|---|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        done = os.path.join(ROOT, "s26", "jobs_done", f"verify_{r['audit']}.json")
        peak = "?"
        if os.path.exists(done):
            with open(done, encoding="utf-8") as fh:
                peak = f"{json.load(fh).get('peak_rss_gb', '?')} GB"
        if not r["diffs"]:
            L.append(f"| {r['audit']} | {r['tag']} | {peak} | {r['wall_s']} s | - | - | | | | | | | {r['error'] or ''} |")
        for b, d in r["diffs"].items():
            L.append(f"| {r['audit']} | {r['tag']} | {peak} | {r['wall_s']} s | {b} | {d.get('verdict')} | "
                     f"{d.get('identical', '')} | {d.get('close', '')} | {d.get('different', '')} | "
                     f"{d.get('run_property_different', '')} | {d.get('n_added', '')} | {d.get('n_removed', '')} | "
                     f"{(r['error'] or '').splitlines()[-1][:80] if r['error'] else ''} |")
    L += ["", "## Not re-runnable", ""]
    for k, v in NOT_RERUNNABLE.items():
        L.append(f"- `verify/{k}`: {v}")
    L += ["", "## First differences per audit", ""]
    for r in rows:
        for b, d in r["diffs"].items():
            if d.get("diffs"):
                L.append(f"### {r['audit']} / {b}")
                for x in d["diffs"][:12]:
                    L.append(f"- `{x['key']}`: tracked `{x['tracked']}`, fresh `{x['fresh']}`")
                L.append("")
    with open(os.path.join(OUT, "REPORT.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    with open(os.path.join(OUT, "REPORT.json"), "w", encoding="utf-8") as fh:
        json.dump({"written": time.strftime("%Y-%m-%dT%H:%M:%S"), "audits": rows,
                   "not_rerunnable": NOT_RERUNNABLE}, fh, indent=1, default=str)
    print("\n".join(L))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("list", "run", "diff", "launch", "report"))
    ap.add_argument("audit", nargs="?")
    ap.add_argument("--only", default="")
    ap.add_argument("--skip", default="")
    ap.add_argument("--include-heavy", action="store_true")
    ap.add_argument("--after", default="", help="comma-separated job names to wait for first")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for n, (s, argv_, tag, est, tr, cost, note) in AUDITS.items():
            print(f"{n:24} {tag:5} est {est:.1f} GB  {cost:6}  {s} {' '.join(argv_)}  -> {tr or 'no tracked json'}")
        for k, v in NOT_RERUNNABLE.items():
            print(f"{'(not re-runnable)':24} {k}: {v}")
        return 0
    if a.cmd == "run":
        return run(a.audit)
    if a.cmd == "diff":
        script, argv_, tag, est, tracked, cost, note = AUDITS[a.audit]
        for b in tracked:
            d = diff_json(os.path.join(VERIFY, b), _redirected(os.path.join(VERIFY, b)))
            print(json.dumps(d, indent=1))
        return 0
    if a.cmd == "launch":
        return launch(only=[x for x in a.only.split(",") if x], skip=[x for x in a.skip.split(",") if x],
                      include_heavy=a.include_heavy, after=[x for x in a.after.split(",") if x])
    return report()


if __name__ == "__main__":
    raise SystemExit(main())
