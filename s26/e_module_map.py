"""s26/e_module_map.py -- the repository import graph, by AST walk, for EXAMINATION.md.

Lane E (Examiner), Sprint 26.  Walks every ``.py`` under the repository root (skipping
``.git``, ``__pycache__``, ``_archive``, anything named ``QUARANTINE*`` and the quarantined
``distogram_models_large.STALE-*``), records for each module its line count, first docstring
line, every import statement (module level AND inside functions, because the production
modules import lazily), and resolves each import to a repository module where one exists.
The reverse graph (``imported_by``) is derived from that.  Also records both
``core.backend_report()`` resolutions (default and ``CORE_BACKENDS=legacy``) by importing
``core`` in two subprocesses, so the backend column of the module map is measured rather
than copied.

Output: ``s26/results/module_map.json``.  Pure stdlib plus one subprocess per backend
report; no database, model or ESM bank is loaded.

    python s26/e_module_map.py
"""
from __future__ import annotations

import ast
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "s26", "results", "module_map.json")
SKIP_DIRS = {".git", "__pycache__", "_archive", "node_modules", ".pytest_cache"}
SKIP_PREFIX = ("QUARANTINE", "distogram_models_large.STALE")

#: the module sets EXAMINATION.md tabulates explicitly
ROOT_REFERENCE = [
    "amber_hamiltonian", "amber_refine", "budget", "catrace", "distogram", "energy_terms",
    "esm_features", "floor", "foldvqe", "fragment_db", "hamiltonian", "legacy_field",
    "legacy_refine", "objective", "pairnet", "peptide_db", "priors", "protein_geometry",
    "qansatz", "refine2", "representations", "sidechains", "torsion_lib2", "vqe",
]


def py_files():
    for dp, dns, fns in os.walk(ROOT):
        rel = os.path.relpath(dp, ROOT)
        parts = [] if rel == "." else rel.split(os.sep)
        if any(p in SKIP_DIRS or p.startswith(SKIP_PREFIX) for p in parts):
            dns[:] = []
            continue
        dns[:] = sorted(d for d in dns
                        if d not in SKIP_DIRS and not d.startswith(SKIP_PREFIX))
        for fn in sorted(fns):
            if fn.endswith(".py"):
                yield os.path.join(dp, fn)


def modname(path):
    rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
    if rel.endswith("/__init__.py"):
        rel = rel[: -len("/__init__.py")]
    elif rel.endswith(".py"):
        rel = rel[:-3]
    return rel.replace("/", ".")


def first_doc_line(tree):
    doc = ast.get_docstring(tree, clean=True)
    if not doc:
        return ""
    return doc.strip().splitlines()[0][:200]


def collect_imports(tree, pkg):
    """Every (dotted name, is_relative, lineno, inside_function) in the module."""
    out = []
    func_depth = 0

    class V(ast.NodeVisitor):
        def visit_FunctionDef(self, node):
            nonlocal func_depth
            func_depth += 1
            self.generic_visit(node)
            func_depth -= 1
        visit_AsyncFunctionDef = visit_FunctionDef

        def visit_Import(self, node):
            for a in node.names:
                out.append((a.name, False, node.lineno, func_depth > 0))
            self.generic_visit(node)

        def visit_ImportFrom(self, node):
            base = node.module or ""
            if node.level:
                # relative import: resolve against the package of this module
                anchor = pkg.split(".")
                anchor = anchor[: len(anchor) - (node.level - 1)] if node.level > 1 else anchor
                base = ".".join([p for p in anchor if p] + ([base] if base else []))
            for a in node.names:
                full = f"{base}.{a.name}" if base else a.name
                out.append((full, bool(node.level), node.lineno, func_depth > 0))
            self.generic_visit(node)

    V().visit(tree)
    return out


def resolve(name, known):
    """Longest known-module prefix of a dotted import name, or None (third party)."""
    parts = name.split(".")
    for k in range(len(parts), 0, -1):
        cand = ".".join(parts[:k])
        if cand in known:
            return cand
    return None


def backend_reports():
    py = sys.executable
    code = "import core, json; print(json.dumps(core.backend_report()))"
    out = {}
    for label, env_extra in (("default", {}), ("legacy", {"CORE_BACKENDS": "legacy"})):
        env = dict(os.environ)
        env.pop("CORE_BACKENDS", None)
        env.update(env_extra)
        r = subprocess.run([py, "-c", code], cwd=ROOT, env=env, capture_output=True,
                           text=True, timeout=300)
        out[label] = json.loads(r.stdout.strip().splitlines()[-1]) if r.returncode == 0 \
            else {"error": r.stderr[-2000:]}
    return out


def main():
    files = list(py_files())
    mods = {}
    for p in files:
        name = modname(p)
        with open(p, "rb") as fh:
            raw = fh.read()
        try:
            tree = ast.parse(raw, filename=p)
        except SyntaxError as exc:
            mods[name] = {"path": os.path.relpath(p, ROOT).replace(os.sep, "/"),
                          "lines": raw.count(b"\n") + (0 if raw.endswith(b"\n") else 1),
                          "syntax_error": str(exc)}
            continue
        pkg = name if p.endswith("__init__.py") else name.rsplit(".", 1)[0] if "." in name else ""
        mods[name] = {
            "path": os.path.relpath(p, ROOT).replace(os.sep, "/"),
            "lines": raw.count(b"\n") + (0 if raw.endswith(b"\n") else 1),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "doc": first_doc_line(tree),
            "_raw_imports": collect_imports(tree, pkg),
        }
    known = set(mods)
    # packages (directories with __init__) are already modules; add bare package names
    # for directories without __init__ so `import s12.instrument` still resolves.
    for name in list(known):
        parts = name.split(".")
        for k in range(1, len(parts)):
            known.add(".".join(parts[:k]))
    for name, m in mods.items():
        if "_raw_imports" not in m:
            continue
        repo, third = [], []
        for full, rel, ln, infn in m.pop("_raw_imports"):
            r = resolve(full, known)
            if r is not None and r != name:
                repo.append({"module": r, "as_written": full, "line": ln,
                             "lazy": infn})
            elif r is None:
                third.append(full.split(".")[0])
        seen = set()
        m["imports_repo"] = [x for x in repo
                             if not (x["module"] in seen or seen.add(x["module"]))]
        m["imports_third_party"] = sorted(set(third))
        m["imported_by"] = []
    for name, m in mods.items():
        for x in m.get("imports_repo", []):
            tgt = x["module"]
            if tgt in mods:
                mods[tgt]["imported_by"].append(name)
            elif tgt + ".__init__" in mods:
                mods[tgt + ".__init__"]["imported_by"].append(name)
    for m in mods.values():
        m["imported_by"] = sorted(set(m.get("imported_by", [])))
        m["n_importers"] = len(m["imported_by"])
    reports = backend_reports()
    # which backend name resolves to which module, under each setting
    backend_of = {}
    for label, rep in reports.items():
        if "error" in rep:
            continue
        for bname, mod in rep.items():
            backend_of.setdefault(mod, {}).setdefault(label, []).append(bname)
    for name, m in mods.items():
        m["backend_names"] = backend_of.get(name, {})
    payload = {
        "root": ROOT,
        "n_modules": len(mods),
        "skipped_dirs": sorted(SKIP_DIRS),
        "skipped_prefixes": list(SKIP_PREFIX),
        "backend_report": reports,
        "root_reference_modules": ROOT_REFERENCE,
        "modules": dict(sorted(mods.items())),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1, sort_keys=False)
    print(f"wrote {OUT}: {len(mods)} modules")
    print("backend_report default:", json.dumps(reports.get("default")))
    print("backend_report legacy: ", json.dumps(reports.get("legacy")))
    # a short table for the terminal
    groups = {"core": [], "root": [], "s5/s7/s8/s9": [], "resultslab": []}
    for name, m in sorted(mods.items()):
        if name.startswith("core."):
            groups["core"].append(name)
        elif name in ROOT_REFERENCE:
            groups["root"].append(name)
        elif name.split(".")[0] in ("s5", "s7", "s8", "s9"):
            groups["s5/s7/s8/s9"].append(name)
        elif name.startswith("s25.resultslab"):
            groups["resultslab"].append(name)
    for g, names in groups.items():
        print(f"\n== {g} ({len(names)}) ==")
        for n in names:
            m = mods[n]
            print(f"{n:32} {m['lines']:5d}  importers={m['n_importers']:3d}  "
                  f"imports={[x['module'] for x in m.get('imports_repo', [])]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
