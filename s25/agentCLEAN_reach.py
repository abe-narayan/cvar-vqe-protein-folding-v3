"""Import-graph reachability analysis over the whole repository.

Builds a static AST import graph of every .py file in the tree, adds the DYNAMIC
edges that `core/__init__.py:REPLACES` creates via importlib (a pure-AST walk
misses them), and computes reachability from four labelled entry-point sets:

    PRODUCTION   core.pipeline, core.bench, core.__init__  (+ every REPLACES target)
    INSTRUMENT   s12.instrument
    TESTS        tests/*.py
    VERIFY       verify/*.py

Writes a machine-readable artefact.  Read-only; touches nothing in the repo.
"""
from __future__ import annotations

import ast
import json
import os
import sys
from collections import deque

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
OUT = sys.argv[2] if len(sys.argv) > 2 else "reachability.json"

SKIP_DIRS = {"__pycache__", ".git", ".pytest_cache", ".vscode", ".venv", "venv"}


def walk_py():
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            if fn.endswith(".py"):
                yield os.path.join(dirpath, fn)


def to_modname(path):
    """Filesystem path -> dotted module name, relative to ROOT."""
    rel = os.path.relpath(path, ROOT).replace("\\", "/")
    if rel.endswith("/__init__.py"):
        rel = rel[: -len("/__init__.py")]
    elif rel.endswith(".py"):
        rel = rel[:-3]
    return rel.replace("/", ".")


FILES = sorted(walk_py())
MOD2PATH = {}
for p in FILES:
    MOD2PATH.setdefault(to_modname(p), p)
PATH2MOD = {v: k for k, v in MOD2PATH.items()}


def resolve(name, cur_mod):
    """Best-effort: dotted name -> a module in this repo, or None (stdlib/3rd-party)."""
    if name in MOD2PATH:
        return name
    # `from core.quantum import X` where X is a symbol, not a module
    parts = name.split(".")
    while parts:
        parts.pop()
        cand = ".".join(parts)
        if cand and cand in MOD2PATH:
            return cand
    return None


def imports_of(path, modname):
    try:
        src = open(path, "r", encoding="utf-8", errors="replace").read()
        tree = ast.parse(src, filename=path)
    except SyntaxError as e:
        return None, f"SyntaxError: {e}"
    out = set()
    pkg = modname.rsplit(".", 1)[0] if "." in modname else ""
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                out.add(a.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:  # relative import
                base = pkg
                for _ in range(node.level - 1):
                    base = base.rsplit(".", 1)[0] if "." in base else ""
                full = f"{base}.{node.module}" if node.module else base
                out.add(full)
                for a in node.names:
                    out.add(f"{full}.{a.name}")
            elif node.module:
                out.add(node.module)
                for a in node.names:
                    out.add(f"{node.module}.{a.name}")
        # importlib.import_module("literal")
        elif isinstance(node, ast.Call):
            f = node.func
            nm = None
            if isinstance(f, ast.Attribute) and f.attr == "import_module":
                nm = "import_module"
            elif isinstance(f, ast.Name) and f.id in ("import_module", "__import__"):
                nm = f.id
            if nm and node.args and isinstance(node.args[0], ast.Constant) \
                    and isinstance(node.args[0].value, str):
                out.add(node.args[0].value)
    return out, None


GRAPH = {}
ERRORS = {}
RAW = {}
for p in FILES:
    m = PATH2MOD[p]
    imps, err = imports_of(p, m)
    if err:
        ERRORS[m] = err
        GRAPH[m] = set()
        continue
    RAW[m] = sorted(imps)
    edges = set()
    for i in imps:
        r = resolve(i, m)
        if r and r != m:
            edges.add(r)
    GRAPH[m] = edges

# ---- dynamic edges from core/__init__.py:REPLACES -------------------------------
DYNAMIC = {}
init = os.path.join(ROOT, "core", "__init__.py")
if os.path.exists(init):
    tree = ast.parse(open(init, encoding="utf-8").read())
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "REPLACES" for t in node.targets):
            d = ast.literal_eval(node.value)
            for k, (opt, leg) in d.items():
                for tgt in (opt, leg):
                    r = resolve(tgt, "core")
                    if r:
                        DYNAMIC.setdefault("core", set()).add(r)
for k, v in DYNAMIC.items():
    GRAPH.setdefault(k, set()).update(v)


def closure(seeds):
    seen, q = set(), deque()
    for s in seeds:
        if s in GRAPH:
            seen.add(s)
            q.append(s)
    while q:
        m = q.popleft()
        for n in GRAPH.get(m, ()):
            if n not in seen:
                seen.add(n)
                q.append(n)
    return seen


PROD_SEEDS = ["core", "core.pipeline", "core.bench"]
INSTR_SEEDS = ["s12.instrument"]
TEST_SEEDS = [m for m in GRAPH if m.startswith("tests.") or m.startswith("tests")]
VERIFY_SEEDS = [m for m in GRAPH if m.startswith("verify.")]

prod = closure(PROD_SEEDS)
instr = closure(INSTR_SEEDS)
tests_c = closure(TEST_SEEDS)
verify_c = closure(VERIFY_SEEDS)

# reverse graph, for "who reaches X"
REV = {m: set() for m in GRAPH}
for m, es in GRAPH.items():
    for e in es:
        REV.setdefault(e, set()).add(m)


def size(m):
    try:
        return os.path.getsize(MOD2PATH[m])
    except OSError:
        return 0


def top(m):
    return m.split(".")[0]


records = {}
for m in sorted(GRAPH):
    labels = []
    if m in prod:
        labels.append("PRODUCTION")
    if m in instr:
        labels.append("INSTRUMENT")
    if m in tests_c:
        labels.append("TEST_CLOSURE")
    if m in verify_c:
        labels.append("VERIFY_CLOSURE")
    records[m] = {
        "path": os.path.relpath(MOD2PATH[m], ROOT).replace("\\", "/"),
        "top": top(m),
        "bytes": size(m),
        "labels": labels,
        "imports_repo": sorted(GRAPH[m]),
        "imported_by": sorted(REV.get(m, ())),
        "n_importers": len(REV.get(m, ())),
    }

unreached = [m for m in records if not records[m]["labels"]]
orphans = [m for m in unreached if records[m]["n_importers"] == 0]

by_top = {}
for m, r in records.items():
    t = r["top"]
    b = by_top.setdefault(t, {"n": 0, "bytes": 0, "prod": 0, "instr": 0,
                              "test": 0, "unreached": 0, "orphan": 0})
    b["n"] += 1
    b["bytes"] += r["bytes"]
    if "PRODUCTION" in r["labels"]:
        b["prod"] += 1
    if "INSTRUMENT" in r["labels"]:
        b["instr"] += 1
    if "TEST_CLOSURE" in r["labels"]:
        b["test"] += 1
    if not r["labels"]:
        b["unreached"] += 1
    if m in orphans:
        b["orphan"] += 1

art = {
    "root": ROOT.replace("\\", "/"),
    "generated_by": "s25 cleanup lane, static AST import graph + core.REPLACES dynamic edges",
    "n_py_files": len(FILES),
    "entry_points": {
        "production": PROD_SEEDS,
        "instrument": INSTR_SEEDS,
        "tests": sorted(TEST_SEEDS),
        "verify": sorted(VERIFY_SEEDS),
    },
    "dynamic_edges": {k: sorted(v) for k, v in DYNAMIC.items()},
    "parse_errors": ERRORS,
    "closures": {
        "production": sorted(prod),
        "instrument": sorted(instr),
        "tests": sorted(tests_c),
        "verify": sorted(verify_c),
    },
    "counts": {
        "production": len(prod), "instrument": len(instr),
        "tests": len(tests_c), "verify": len(verify_c),
        "unreached": len(unreached), "orphans": len(orphans),
    },
    "by_top_level_dir": dict(sorted(by_top.items())),
    "unreached_modules": sorted(unreached),
    "orphan_modules": sorted(orphans),
    "modules": records,
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(art, f, indent=1, sort_keys=False)

print(f"files={len(FILES)} parse_errors={len(ERRORS)}")
print(f"production closure = {len(prod)}")
print(f"instrument closure = {len(instr)}")
print(f"tests closure      = {len(tests_c)}")
print(f"verify closure     = {len(verify_c)}")
print(f"unreached          = {len(unreached)}  (orphans, nobody imports = {len(orphans)})")
print()
print("PRODUCTION closure:")
for m in sorted(prod):
    print("   ", m)
print()
print("INSTRUMENT closure (minus production):")
for m in sorted(instr - prod):
    print("   ", m)
print()
print("TEST closure (minus production+instrument):")
for m in sorted(tests_c - prod - instr):
    print("   ", m)
print()
print("by top-level dir:")
for t, b in sorted(by_top.items()):
    print(f"   {t:<14} py={b['n']:<4} {b['bytes']/1e6:6.2f}MB  prod={b['prod']:<3} "
          f"instr={b['instr']:<3} test={b['test']:<3} unreached={b['unreached']:<4} orphan={b['orphan']}")
