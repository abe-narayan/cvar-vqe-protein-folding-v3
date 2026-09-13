#!/usr/bin/env python
"""s26/i_ast_check.py -- is a production-path module AST-identical to its committed parent
once docstrings are stripped?  If not, print exactly what changed.

The S26 lane contract (section 5) requires every `.py` on the production path to stay
AST-identical to its parent modulo docstrings unless the change is deliberate, tested and in
the ledger.  This tool makes the check runnable and the deviation readable: it strips
docstrings from both trees, `ast.unparse`s them, and prints a unified diff.

    python s26/i_ast_check.py core/data.py core/amber.py ...        # against HEAD
    python s26/i_ast_check.py --ref ae86a124 core/data.py            # against a commit
"""
from __future__ import annotations

import argparse
import ast
import difflib
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def strip_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, "body", None)
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                node.body = body[1:] or [ast.Pass()]
    return tree


def canon(src, name):
    return ast.unparse(strip_docstrings(ast.parse(src, filename=name)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--ref", default="HEAD")
    a = ap.parse_args(argv)
    rc = 0
    for rel in a.paths:
        rel = rel.replace("\\", "/")
        old = subprocess.run(["git", "show", f"{a.ref}:{rel}"], cwd=ROOT, capture_output=True,
                             text=True, encoding="utf-8")
        if old.returncode != 0:
            print(f"{rel}: not in {a.ref} ({old.stderr.strip()[:80]})")
            rc = 1
            continue
        with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
            new = fh.read()
        co, cn = canon(old.stdout, rel), canon(new, rel)
        if co == cn:
            print(f"{rel}: AST-identical to {a.ref} modulo docstrings")
            continue
        rc = 1
        diff = list(difflib.unified_diff(co.splitlines(), cn.splitlines(),
                                         f"{a.ref}:{rel}", f"working:{rel}", lineterm="", n=1))
        print(f"{rel}: DIFFERS from {a.ref} modulo docstrings ({sum(1 for l in diff if l[:1] in '+-' and l[:3] not in ('+++', '---'))} changed lines)")
        for line in diff:
            print("    " + line)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
