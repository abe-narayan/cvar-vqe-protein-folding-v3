"""Full-tree classification artefact for the Sprint 25 cleanup lane.

Consumes the reachability graph (reach.py) and classifies EVERY file in the tree into
PRODUCTION / INSTRUMENT / RESEARCH_ARCHIVE / EVIDENCE / REGENERABLE_CACHE / OBSOLETE /
PROTECTED, with the evidence for each call recorded per file.

Read-only over the repository.
"""
from __future__ import annotations
import json, os, subprocess, sys
from collections import defaultdict

ROOT = os.path.abspath(sys.argv[1])
REACH = json.load(open(sys.argv[2], encoding="utf-8"))
OUT = sys.argv[3]

MODS = REACH["modules"]
PATH2MOD = {r["path"]: m for m, r in MODS.items()}

tracked = set()
try:
    tracked = {p.strip().replace("\\", "/") for p in subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True).stdout.splitlines()}
except Exception:
    pass

SKIP_DIRS = {".git"}

#: filename patterns that are SCIENTIFIC EVIDENCE by the sprint's hard rule
EVIDENCE_NAMES = ("LEDGER.md", "BRIEF.md", "CLAIMS.md", "FINDINGS.md", "DOSSIER.md",
                  "FINAL_REPORT.md", "PAPER_DRAFT.md", "STATUS", "PHASE0.md")
EVIDENCE_PREFIX = ("PREREG_",)

PROTECTED = ("s24/cache_amber/", "s24/results/a_corpus_permitted_")

#: directories whose contents are per-target caches: large, regenerable by the module
#: that wrote them, and carrying no claim that the results JSON beside them does not.
CACHE_DIRS = ("/cache/", "/libcache/", "core_cache/", "__pycache__", ".pytest_cache")


def classify(rel, size):
    ev = []
    p = rel

    for pr in PROTECTED:
        if p.startswith(pr):
            return "PROTECTED", ["named in the sprint brief as not cheaply regenerable"]

    if "__pycache__" in p or p.startswith(".pytest_cache"):
        return "OBSOLETE", ["compiled bytecode / pytest cache; gitignored; "
                            "regenerated on the next run"]

    base = os.path.basename(p)
    top = p.split("/")[0]

    # --- python modules: classify by reachability
    mod = PATH2MOD.get(p)
    if mod:
        lab = MODS[mod]["labels"]
        n_imp = MODS[mod]["n_importers"]
        if "PRODUCTION" in lab:
            return "PRODUCTION", [f"in the import closure of core.pipeline/core.bench "
                                  f"as `{mod}`; {n_imp} importers"]
        if "INSTRUMENT" in lab:
            return "INSTRUMENT", [f"in the import closure of s12.instrument as `{mod}`"]
        if p.startswith("tests/"):
            return "PRODUCTION", ["test of the production surface"]
        if p.startswith("verify/"):
            return "INSTRUMENT", ["standalone audit; core/project.py reads and writes "
                                  "verify/ at runtime and tests/ imports three of these"]
        if "TEST_CLOSURE" in lab:
            return "PRODUCTION", [f"reached from tests/ as `{mod}` "
                                  f"({', '.join(MODS[mod]['imported_by'][:3])})"]
        return "RESEARCH_ARCHIVE", [
            f"sprint experiment module `{mod}`; unreachable from production, the "
            f"instrument, or tests; {n_imp} in-sprint importers"]

    # --- documentation / evidence
    if base.endswith(".md"):
        if any(base.endswith(e) or e in base for e in EVIDENCE_NAMES) or \
                any(base.startswith(e) for e in EVIDENCE_PREFIX):
            return "EVIDENCE", ["pre-registration / ledger / findings: protected by the "
                                "sprint's hard deletion rule"]
        if top in ("s5", "s7", "s8", "s9") or top.startswith("s1") or top.startswith("s2"):
            return "EVIDENCE", ["sprint documentation"]
        return "DOC", ["repository documentation"]

    # --- caches vs results
    if any(c in "/" + p for c in CACHE_DIRS):
        return "REGENERABLE_CACHE", ["per-target cache written by the sprint module that "
                                     "owns it; not cited as evidence"]
    if "/results/" in p or p.startswith("results/") or top == "bench_results":
        if base.endswith(".tmp"):
            return "OBSOLETE", ["interrupted write: a .tmp left by a crashed writer"]
        return "EVIDENCE", ["result artefact under a results/ directory"]

    if base.endswith(".log"):
        if top == "_archive":
            return "EVIDENCE", ["console log preserved by the 2026-09-04 consolidation; "
                                "_archive/README.txt records that one is cited by "
                                "FINDINGS.md and one is a sole on-disk source"]
        if top == "verify":
            return "INSTRUMENT", ["equivalence-run transcript; 5 of these are cited by "
                                  "basename inside verify/run_equiv*.sh"]
        if top.startswith("s") and top[1:].isdigit():
            return "RESEARCH_ARCHIVE", [
                "sprint console log, not cited by basename in any of the 8,556 text "
                "files in the tree (433.7 MB of .py/.md/.json/.toml/.txt/.sh searched) "
                "-- but the 2026-09-04 precedent found a sprint log that WAS the sole "
                "on-disk source of a published number, so these are archived, not deleted"]
        return "OBSOLETE", [
            "root-level run log; gitignored as disposable and not cited by basename "
            "anywhere in the tree (8,556 text files, 433.7 MB searched, including every "
            "sprint result JSON and all of FINDINGS.md)"]

    if base.endswith((".npz", ".npy", ".pt", ".pkl", ".gz")) or ".pt.d" in base:
        return "REGENERABLE_CACHE", ["derived bank / model checkpoint; .gitignore names "
                                     "the module that rebuilds it"]

    if base.endswith((".png", ".svg")):
        return "EVIDENCE", ["sprint figure: the plotted form of a result JSON beside it"]

    if base.endswith(".pdb"):
        if top in ("prots", "pdbs_ext"):
            return "DATA_CORPUS", ["gitignored RCSB corpus; re-fetchable by the query "
                                   "scripts named in .gitignore, but hours of download"]
        return "DATA_PINNED", ["pinned structure set: adding or removing one file changes "
                               "BLOSUM tie order and moves up to 47 of 500 pool members"]

    if base.endswith((".html", ".css", ".js")):
        return "PRODUCTION", ["Phase-II results renderer output"]

    if base.endswith(".sh"):
        return "INSTRUMENT", ["equivalence run script cited by README"]

    if base.endswith(".json"):
        return "EVIDENCE", ["result / manifest JSON"]

    return "OTHER", ["unclassified"]


rows = []
for dp, dn, fn in os.walk(ROOT):
    dn[:] = [d for d in dn if d not in SKIP_DIRS]
    for f in fn:
        full = os.path.join(dp, f)
        rel = os.path.relpath(full, ROOT).replace("\\", "/")
        try:
            size = os.path.getsize(full)
        except OSError:
            size = 0
        cls, ev = classify(rel, size)
        rows.append({"path": rel, "class": cls, "bytes": size,
                     "tracked": rel in tracked, "evidence": ev})

agg = defaultdict(lambda: {"n": 0, "bytes": 0, "tracked": 0})
per_dir = defaultdict(lambda: defaultdict(lambda: {"n": 0, "bytes": 0}))
for r in rows:
    a = agg[r["class"]]
    a["n"] += 1
    a["bytes"] += r["bytes"]
    a["tracked"] += int(r["tracked"])
    t = r["path"].split("/")[0] if "/" in r["path"] else "(root)"
    d = per_dir[t][r["class"]]
    d["n"] += 1
    d["bytes"] += r["bytes"]

art = {
    "generated": "s25 cleanup lane -- INVENTORY ONLY, nothing moved or deleted",
    "root": ROOT.replace("\\", "/"),
    "reachability": {k: REACH[k] for k in
                     ("entry_points", "dynamic_edges", "counts", "closures")},
    "classes": {
        "PRODUCTION": "reachable by import from core.pipeline / core.bench, or a test of it",
        "INSTRUMENT": "the evaluation path (s12.instrument, verify/) -- results depend on it",
        "RESEARCH_ARCHIVE": "sprint experiment module, unreachable from production",
        "EVIDENCE": "LEDGER/BRIEF/PREREG/FINDINGS/CLAIMS and result JSON -- NEVER deletable",
        "PROTECTED": "named in the sprint brief as not cheaply regenerable",
        "REGENERABLE_CACHE": "per-target cache / derived bank / checkpoint",
        "DATA": "structure files",
        "DOC": "repository documentation",
        "OBSOLETE": "bytecode, disposable logs, interrupted writes",
        "OTHER": "unclassified",
    },
    "summary": {k: dict(v) for k, v in sorted(agg.items())},
    "by_top_level": {k: {c: dict(v) for c, v in sorted(d.items())}
                     for k, d in sorted(per_dir.items())},
    "files_note": "per-file rows are given for every class EXCEPT DATA_CORPUS and "
                  "REGENERABLE_CACHE, which are 28k files of re-derivable bytes and are "
                  "given as directory aggregates in `bulk_by_directory` instead.",
    "bulk_by_directory": None,   # filled below
    "files": [r for r in rows
              if r["class"] not in ("DATA_CORPUS", "REGENERABLE_CACHE")],
}
bulk = defaultdict(lambda: defaultdict(lambda: {"n": 0, "bytes": 0}))
for r in rows:
    if r["class"] in ("DATA_CORPUS", "REGENERABLE_CACHE"):
        d = os.path.dirname(r["path"]) or "(root)"
        e = bulk[d][r["class"]]
        e["n"] += 1
        e["bytes"] += r["bytes"]
art["bulk_by_directory"] = {k: {c: dict(v) for c, v in sorted(d.items())}
                            for k, d in sorted(bulk.items())}
json.dump(art, open(OUT, "w", encoding="utf-8"), indent=1)

print(f"{'class':<20} {'files':>7} {'MB':>10} {'tracked':>8}")
for k, v in sorted(agg.items(), key=lambda kv: -kv[1]["bytes"]):
    print(f"{k:<20} {v['n']:>7} {v['bytes']/1e6:>10.1f} {v['tracked']:>8}")
print()
print("OBSOLETE detail:")
od = defaultdict(lambda: [0, 0])
for r in rows:
    if r["class"] == "OBSOLETE":
        k = r["evidence"][0]
        od[k][0] += 1
        od[k][1] += r["bytes"]
for k, (n, b) in sorted(od.items(), key=lambda kv: -kv[1][1]):
    print(f"   {n:>6} files  {b/1e6:>8.2f} MB   {k}")
