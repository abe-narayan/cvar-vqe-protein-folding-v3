"""Conservative AST hygiene pass: unused imports, dead private helpers, broken
literal paths, and duplicate top-level function names across the production surface.

Conservative by design -- anything that could be reached dynamically (getattr,
__all__, a string literal naming it, a star-import in the file) is NOT reported.
"""
from __future__ import annotations
import ast, json, os, re, sys

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
TARGETS = sys.argv[2:] or ["core", "s12/instrument.py", "tests", "verify"]


def files():
    for t in TARGETS:
        p = os.path.join(ROOT, t)
        if os.path.isfile(p):
            yield p
        else:
            for dp, dn, fn in os.walk(p):
                dn[:] = [d for d in dn if d != "__pycache__"]
                for f in fn:
                    if f.endswith(".py"):
                        yield os.path.join(dp, f)


UNUSED, DEAD, PATHS, STAR = [], [], [], []
NAMES = {}

for path in sorted(files()):
    rel = os.path.relpath(path, ROOT).replace("\\", "/")
    src = open(path, encoding="utf-8", errors="replace").read()
    try:
        tree = ast.parse(src, filename=path)
    except SyntaxError as e:
        print(f"SYNTAX {rel}: {e}")
        continue

    has_star = any(isinstance(n, ast.ImportFrom) and any(a.name == "*" for a in n.names)
                   for n in ast.walk(tree))
    if has_star:
        STAR.append(rel)

    # every Name/Attribute-root identifier actually used
    used = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Name):
            used.add(n.id)
        elif isinstance(n, ast.Attribute):
            b = n
            while isinstance(b, ast.Attribute):
                b = b.value
            if isinstance(b, ast.Name):
                used.add(b.id)
    # anything mentioned in a string (getattr / __all__ / doctest / f-string) is spared
    strings = " ".join(n.value for n in ast.walk(tree)
                       if isinstance(n, ast.Constant) and isinstance(n.value, str))

    # ---- unused imports (module level only; local imports are deliberate lazy loads)
    if not has_star:
        for n in tree.body:
            if isinstance(n, (ast.Import, ast.ImportFrom)):
                for a in n.names:
                    if a.name == "*":
                        continue
                    bind = a.asname or a.name.split(".")[0]
                    if bind in ("annotations",):
                        continue
                    if bind not in used and not re.search(r"\b%s\b" % re.escape(bind), strings):
                        UNUSED.append((rel, n.lineno, bind, a.name))

    # ---- dead private module-level helpers
    defs = {}
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            defs[n.name] = n.lineno
            NAMES.setdefault(n.name, []).append(rel)
    for name, ln in defs.items():
        if not name.startswith("_") or name.startswith("__"):
            continue
        # count Name loads of it excluding its own def
        n_uses = sum(1 for n in ast.walk(tree)
                     if isinstance(n, ast.Name) and n.id == name)
        if n_uses == 0 and not re.search(r"\b%s\b" % re.escape(name), strings):
            # also make sure nothing else in the repo imports it
            DEAD.append((rel, ln, name))

    # ---- literal paths that do not exist
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, str):
            v = n.value
            if ("/" in v or "\\" in v) and not v.startswith(("http", "{", "%")) \
                    and re.search(r"\.(json|npz|npy|md|pdb|csv|pt|txt|xml)$", v) \
                    and "*" not in v and "{" not in v and "%" not in v:
                cand = os.path.join(ROOT, v.replace("\\", "/"))
                if not os.path.exists(cand):
                    PATHS.append((rel, n.lineno, v))

# cross-file check for DEAD: is the private name imported anywhere in the repo?
allsrc = ""
for dp, dn, fn in os.walk(ROOT):
    dn[:] = [d for d in dn if d not in ("__pycache__", ".git", ".pytest_cache")]
    for f in fn:
        if f.endswith(".py"):
            try:
                allsrc += open(os.path.join(dp, f), encoding="utf-8", errors="replace").read()
            except OSError:
                pass
DEAD = [(r, l, n) for (r, l, n) in DEAD
        if len(re.findall(r"\b%s\b" % re.escape(n), allsrc)) <= 1]

print("=" * 78)
print(f"UNUSED MODULE-LEVEL IMPORTS: {len(UNUSED)}")
for r, l, b, full in UNUSED:
    print(f"  {r}:{l}  {b}   (import {full})")
print()
print(f"DEAD PRIVATE MODULE-LEVEL FUNCTIONS (0 references repo-wide): {len(DEAD)}")
for r, l, n in DEAD:
    print(f"  {r}:{l}  {n}")
print()
print(f"LITERAL PATHS THAT DO NOT RESOLVE FROM REPO ROOT: {len(PATHS)}")
for r, l, v in sorted(set(PATHS)):
    print(f"  {r}:{l}  {v}")
print()
print(f"FILES WITH STAR IMPORTS (unused-import check skipped there): {STAR}")
print()
dup = {k: v for k, v in NAMES.items() if len(set(v)) > 1 and not k.startswith("test_")}
print(f"TOP-LEVEL FUNCTION NAMES DEFINED IN >1 FILE OF THE SCANNED SET: {len(dup)}")
for k, v in sorted(dup.items()):
    if len(set(v)) >= 2:
        print(f"  {k:<34} {sorted(set(v))}")

json.dump({"unused_imports": UNUSED, "dead_private": DEAD,
           "broken_literal_paths": sorted(set(PATHS)), "star_imports": STAR,
           "duplicate_toplevel_names": {k: sorted(set(v)) for k, v in dup.items()}},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "hygiene.json"), "w"),
          indent=1)
