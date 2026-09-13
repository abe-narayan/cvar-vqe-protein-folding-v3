#!/usr/bin/env python
"""s26/i_claim_check.py -- re-read every claimed number from the artefact it is claimed from.

Rule 6 of the S26 contract: never state a number without its artefact path.  This turns that
rule into a check.  `s26/results/claims.json` lists claims as {id, claim, value, artefact,
locator[, tol][, optional][, fallback_artefact]}; for each one the artefact is opened and the
locator applied:

    json_path   {"type": "json_path", "path": "science_delta.rmsd_arm"}     dotted; ints index lists
    json_row    {"type": "json_row", "rows": "rows", "match": {...}, "key": "mean"}
    glob_mean   {"type": "glob_mean", "key": "rmsd_arm", "n_expected": 126}   artefact is a glob
    regex       {"type": "regex", "pattern": "..."}       the pattern must match somewhere
    text        {"type": "text"}                          `value` must appear verbatim
    sha256      {"type": "sha256"}                        `value` is the hex digest (bytes only)

Numbers compare within `tol` (default 1e-9); strings compare exactly.  A missing artefact is
ABSENT (and the fallback artefact, if given, is tried); `optional: true` turns ABSENT into a
warning instead of a failure.  Output: a table, `s26/results/claim_check.json`, exit 1 on any
MISMATCH or non-optional ABSENT.

    python s26/i_claim_check.py [--claims s26/results/claims.json]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import re
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CLAIMS = os.path.join(ROOT, "s26", "results", "claims.json")
OUT = os.path.join(ROOT, "s26", "results", "claim_check.json")


def _walk(obj, path):
    for seg in [s for s in path.split(".") if s != ""]:
        if isinstance(obj, list):
            obj = obj[int(seg)]
        else:
            obj = obj[seg]
    return obj


def _read_json(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def _sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def locate(claim, artefact):
    """Return (found_value, note).  Raises FileNotFoundError when the artefact is absent."""
    loc = claim["locator"]
    t = loc["type"]
    if t == "glob_mean":
        files = sorted(glob.glob(os.path.join(ROOT, artefact)))
        if not files:
            raise FileNotFoundError(artefact)
        vals = [_read_json(f)[loc["key"]] for f in files]
        n_exp = loc.get("n_expected")
        note = f"n={len(vals)}" + ("" if n_exp is None or n_exp == len(vals)
                                   else f" (EXPECTED {n_exp})")
        if n_exp is not None and n_exp != len(vals):
            return None, note
        return float(sum(vals)) / len(vals), note
    p = os.path.join(ROOT, artefact)
    if not os.path.isfile(p):
        raise FileNotFoundError(artefact)
    if t == "sha256":
        return _sha256(p), f"{os.path.getsize(p)} bytes"
    if t == "json_path":
        return _walk(_read_json(p), loc["path"]), ""
    if t == "json_row":
        d = _read_json(p)
        rows = _walk(d, loc["rows"]) if loc.get("rows") else d
        for r in rows:
            if all(r.get(k) == v for k, v in loc["match"].items()):
                return r.get(loc["key"]), ""
        return None, f"no row matches {loc['match']}"
    with open(p, encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    if t == "regex":
        m = re.search(loc["pattern"], txt)
        return (m.group(0) if m else None), ("" if m else "pattern not found")
    if t == "text":
        return (claim["value"] if str(claim["value"]) in txt else None), ""
    raise ValueError(f"unknown locator type {t!r}")


def compare(expected, found, tol):
    if found is None:
        return False
    if isinstance(expected, bool) or isinstance(found, bool):
        return bool(expected) == bool(found)
    if isinstance(expected, (int, float)) and isinstance(found, (int, float)):
        return abs(float(expected) - float(found)) <= tol
    return str(expected) == str(found)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--claims", default=DEFAULT_CLAIMS)
    a = ap.parse_args(argv)
    spec = _read_json(a.claims)
    rows = []
    n_bad = 0
    print(f"{'id':44} {'status':9} {'expected':>22} {'found':>22}  artefact")
    for c in spec["claims"]:
        tol = float(c.get("tol", 1e-9))
        status, found, note, used = "ERROR", None, "", c["artefact"]
        try:
            found, note = locate(c, c["artefact"])
            status = "OK" if compare(c["value"], found, tol) else "MISMATCH"
        except FileNotFoundError:
            status = "ABSENT"
            fb = c.get("fallback_artefact")
            if fb:
                try:
                    found, note = locate(c, fb)
                    used = fb
                    status = ("OK(fallback)" if compare(c["value"], found, tol)
                              else "MISMATCH(fallback)")
                except FileNotFoundError:
                    status = "ABSENT"
        except Exception as exc:                                         # noqa: BLE001
            status, note = "ERROR", f"{type(exc).__name__}: {exc}"[:120]
        bad = status.startswith("MISMATCH") or status == "ERROR" or (
            status == "ABSENT" and not c.get("optional", False))
        n_bad += int(bad)
        exp_s, fnd_s = str(c["value"])[:22], ("" if found is None else str(found))[:22]
        print(f"{c['id'][:44]:44} {status:9} {exp_s:>22} {fnd_s:>22}  {used}"
              + (f"  [{note}]" if note else ""))
        rows.append({"id": c["id"], "claim": c.get("claim", ""), "status": status,
                     "expected": c["value"], "found": found, "artefact": used, "note": note,
                     "optional": bool(c.get("optional", False))})
    payload = {"checked": time.strftime("%Y-%m-%dT%H:%M:%S"), "claims_file":
               os.path.relpath(a.claims, ROOT).replace(os.sep, "/"), "n": len(rows),
               "n_failing": n_bad, "rows": rows}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    print(f"\n{len(rows)} claims, {n_bad} failing -> {os.path.relpath(OUT, ROOT)}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
