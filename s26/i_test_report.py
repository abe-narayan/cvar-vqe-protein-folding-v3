#!/usr/bin/env python
"""s26/i_test_report.py -- per-file test counts, skip reasons and peak RSS for the S26 record.

Reads a pytest junit XML written by a job that ran under `s26/jobrun.py`, joins it with that
job's `s26/jobs_done/<name>.json` (exit code, wall, PEAK RSS), and merges the result into
`s26/results/test_run.json`; then rewrites `s26/TEST_RUN.md` from every job recorded so far.
A skip whose message names the memory ceiling is counted separately from a real skip.

    python s26/i_test_report.py record --xml s26/results/pytest_core.xml --job pytest_core \
        --commit a4db170c --note "non-AMBER files"
    python s26/i_test_report.py render
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s24 import stats_lib as ST                                          # noqa: E402

RES = os.path.join(ROOT, "s26", "results")
OUT_JSON = os.path.join(RES, "test_run.json")
OUT_MD = os.path.join(ROOT, "s26", "TEST_RUN.md")
DONE = os.path.join(ROOT, "s26", "jobs_done")


def is_memory_guard(msg: str) -> bool:
    m = (msg or "").lower()
    return "ceiling" in m and ("memory" in m or "physical" in m)


def parse_xml(path):
    tree = ET.parse(path)
    files = {}
    for tc in tree.iter("testcase"):
        cls = tc.get("classname", "")
        fname = cls.split(".")[-1] + ".py" if cls else tc.get("file", "?")
        f = files.setdefault(fname, {"passed": 0, "failed": 0, "errors": 0, "skipped": 0,
                                     "memory_guard_skips": 0, "skip_reasons": [],
                                     "failures": []})
        sk = tc.find("skipped")
        fl = tc.find("failure")
        er = tc.find("error")
        if er is not None:
            f["errors"] += 1
            f["failures"].append({"test": tc.get("name"), "kind": "error",
                                  "message": (er.get("message") or "")[:300]})
        elif fl is not None:
            f["failed"] += 1
            f["failures"].append({"test": tc.get("name"), "kind": "failure",
                                  "message": (fl.get("message") or "")[:300]})
        elif sk is not None:
            f["skipped"] += 1
            msg = sk.get("message") or ""
            f["skip_reasons"].append({"test": tc.get("name"), "reason": msg[:300],
                                      "memory_guard": is_memory_guard(msg)})
            f["memory_guard_skips"] += int(is_memory_guard(msg))
        else:
            f["passed"] += 1
    for f in files.values():
        f["total"] = f["passed"] + f["failed"] + f["errors"] + f["skipped"]
    return dict(sorted(files.items()))


def totals(files):
    t = {"passed": 0, "failed": 0, "errors": 0, "skipped": 0, "memory_guard_skips": 0, "total": 0}
    for f in files.values():
        for k in t:
            t[k] += f[k]
    return t


def load():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as fh:
            return json.load(fh)
    return {"label": "S26 test suite under the governor (lane I)", "jobs": {}}


def render(state):
    L = ["# S26 TEST RUN (under the governor)", "",
         "Every job below ran through `s26/jobrun.py` (registered with the governor, stdout in "
         "`s26/logs/<job>.log`, exit code / wall / PEAK RSS in `s26/jobs_done/<job>.json`). "
         "Counts come from the junit XML each job wrote to `s26/results/<job>.xml`; the machine-"
         "readable form is `s26/results/test_run.json`. A skip whose message names the memory "
         "ceiling is a memory-guard skip and is counted apart from real skips.", ""]
    grand = {"passed": 0, "failed": 0, "errors": 0, "skipped": 0, "memory_guard_skips": 0,
             "total": 0}
    for name, j in state["jobs"].items():
        t = j["totals"]
        L += [f"## job `{name}`", "",
              f"- command: `{' '.join(j.get('cmd', []))}`",
              f"- commit: `{j.get('commit', '?')}`  note: {j.get('note', '')}",
              f"- start {j.get('start', '?')}  end {j.get('end', '?')}  wall {j.get('wall_s', '?')} s  "
              f"exit {j.get('exit_code', '?')}  peak RSS {j.get('peak_rss_gb', '?')} GB  "
              f"tag {j.get('tag', '?')}  est {j.get('est_ram_gb', '?')} GB",
              f"- totals: {t['total']} tests, {t['passed']} passed, {t['failed']} failed, "
              f"{t['errors']} errors, {t['skipped']} skipped ({t['memory_guard_skips']} memory-guard)",
              "", "| file | tests | passed | failed | errors | skipped | memory-guard skips |",
              "|---|---:|---:|---:|---:|---:|---:|"]
        for fn, f in j["files"].items():
            L.append(f"| `{fn}` | {f['total']} | {f['passed']} | {f['failed']} | {f['errors']} | "
                     f"{f['skipped']} | {f['memory_guard_skips']} |")
        reasons = [(fn, r) for fn, f in j["files"].items() for r in f["skip_reasons"]]
        if reasons:
            L += ["", "skip reasons (`-rs`):", ""]
            for fn, r in reasons:
                tag = "MEMORY-GUARD" if r["memory_guard"] else "real"
                L.append(f"- `{fn}::{r['test']}` [{tag}]: {r['reason']}")
        fails = [(fn, x) for fn, f in j["files"].items() for x in f["failures"]]
        if fails:
            L += ["", "failures / errors:", ""]
            for fn, x in fails:
                L.append(f"- `{fn}::{x['test']}` ({x['kind']}): {x['message']}")
        L.append("")
        if j.get("counts_in_grand_total", True):
            for k in grand:
                grand[k] += t[k]
    L += ["## Combined (jobs marked as counting toward the suite total)", "",
          f"**{grand['total']} tests: {grand['passed']} passed, {grand['failed']} failed, "
          f"{grand['errors']} errors, {grand['skipped']} skipped "
          f"({grand['memory_guard_skips']} memory-guard skips).**", "",
          f"Rendered {time.strftime('%Y-%m-%d %H:%M')} by `s26/i_test_report.py`."]
    state["combined"] = grand
    with open(OUT_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("record", "render"))
    ap.add_argument("--xml")
    ap.add_argument("--job")
    ap.add_argument("--commit", default=None)
    ap.add_argument("--note", default="")
    ap.add_argument("--no-grand-total", action="store_true",
                    help="record the job but do not count it toward the combined total "
                         "(a re-run of a file already counted)")
    a = ap.parse_args(argv)
    state = load()
    if a.cmd == "record":
        files = parse_xml(a.xml)
        done = {}
        dp = os.path.join(DONE, f"{a.job}.json")
        if os.path.exists(dp):
            with open(dp) as fh:
                done = json.load(fh)
        commit = a.commit or subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                                            capture_output=True, text=True).stdout.strip()
        state["jobs"][a.job] = {
            "xml": os.path.relpath(a.xml, ROOT).replace(os.sep, "/"), "files": files,
            "totals": totals(files), "commit": commit, "note": a.note,
            "counts_in_grand_total": not a.no_grand_total,
            "cmd": done.get("cmd", []), "tag": done.get("tag"), "est_ram_gb": done.get("est_ram_gb"),
            "start": done.get("start"), "end": done.get("end"), "wall_s": done.get("wall_s"),
            "exit_code": done.get("exit_code"), "peak_rss_gb": done.get("peak_rss_gb"),
            "waited_s": done.get("waited_s"),
        }
    render(state)
    ST.save_atomic(OUT_JSON, state, module_file=__file__)
    print(f"wrote {OUT_JSON} and {OUT_MD}: jobs {list(state['jobs'])}; combined {state['combined']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
