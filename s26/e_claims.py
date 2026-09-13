"""s26/e_claims.py -- locate every numerical claim of EXAMINATION.md section C in the repo.

Lane E (Examiner), Sprint 26.  Two passes over the repository, both read-only:

  1. LITERAL: grep every text file (md/py/csv/txt/json under a size cap) for the claim as it
     is quoted in the brief.  This finds where a number is CITED.
  2. STORED: walk every JSON artefact under the results directories and report every leaf
     whose value rounds to the claim at the quoted precision, with its key path and the value
     at STORED precision.  This finds where a number LIVES.

Never opens ``results/benchmark_manifest.json`` (Rule 1), ``pdbs/``, ``s26/`` (this sprint's
own outputs, to avoid citing ourselves), ``.git``, ``_archive``, or any binary.  Writes
``s26/results/claim_search.json`` and prints a compact table.

    python s26/e_claims.py
"""
from __future__ import annotations

import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "s26", "results", "claim_search.json")
SKIP_PARTS = {".git", "__pycache__", "_archive", "node_modules", ".pytest_cache", "s26",
              "pdbs", "distogram_models", "prots", "site", "structures"}
#: Rule 1: the sealed benchmark's manifest, its per-target cache and its final report hold
#: benchmark RMSDs and are never parsed here.  (The first run of this script did parse
#: `s9/final_report.json` before this exclusion existed and printed its published aggregate
#: `dist/shipped/mean`; that incident is recorded in s26/agentE_FINDINGS.md.)
SKIP_FILES = {"results/benchmark_manifest.json", "results/monomer_manifest.json",
              "s9/final_report.json", "s9/final_synth.json"}
SKIP_SUBSTR = ("benchmark", "s9/final_cache", "bench60", "benchmark60")
TEXT_EXT = {".md", ".py", ".csv", ".txt", ".json", ".yaml", ".yml", ".toml", ".cfg"}
MAX_TEXT_BYTES = 6_000_000
MAX_JSON_BYTES = 80_000_000
JSON_DIRS = ["bench_results", "results/summary", "s7", "s8", "s9", "s12", "s13", "s14",
             "s15", "s16", "s17", "s18", "s19", "s20", "s21", "s22", "s23", "s24", "s25",
             "verify"]

#: TARGETED lookups: where a document cites a number, look in THAT sprint's results
#: directory for a leaf that rounds to it, indexed paths included.  Each spec is
#: (claim id, directory, matcher, optional key-substring filter); a matcher is
#: {"val", "dec"} (round to `dec` places), {"val", "rel"} (relative tolerance) or
#: {"lo", "hi"} (a range).
TARGETED = [
    ("C10", "s13/results", {"val": 3.770, "dec": 3}, None),
    ("C11", "s24/results", {"val": 2.2261, "dec": 4}, None),
    ("C13", "s13/results", {"val": 1.594, "dec": 3}, None),
    ("C14", "s15/results", {"val": 0.611, "dec": 3}, None),
    ("C15", "s24/results", {"val": -2.1496, "dec": 4}, None),
    ("C16", "s24/results", {"val": 0.0225, "dec": 4}, None),
    ("C17", "s23/results", {"val": 0.676, "dec": 3}, None),
    ("C18", "s25/results", {"val": 0.902, "dec": 3}, None),
    ("C19", "s25/results", {"val": -0.0302, "dec": 4}, None),
    ("C20", "s25/results", {"val": 0.0406, "dec": 4}, None),
    ("C20", "s25/results", {"val": 0.0118, "dec": 4}, None),
    ("C21", "s25/results", {"val": 5.6e-17, "rel": 0.05}, None),
    ("C22", "s25/results", {"lo": 0.9999999995, "hi": 1.0000000005}, "cos"),
    ("C23", "s25/results", {"val": 4.663e-10, "rel": 0.05}, None),
    ("C24", "s25/results", {"val": 0.566586, "dec": 6}, None),
    ("C24", "s8", {"val": 0.524, "dec": 3}, "tail"),
    ("C25", "s25/results", {"val": -0.6492, "dec": 4}, None),
    ("C25", "s25/results", {"val": -0.2522, "dec": 4}, None),
    ("C25", "s25/results", {"val": -0.0472, "dec": 4}, None),
    ("C25", "s25/results", {"val": -0.3105, "dec": 4}, None),
    ("C26", "s13/results", {"val": 36.1, "dec": 1}, "phi"),
    ("C26", "s13/results", {"val": 36.4, "dec": 1}, "phi"),
    ("C27", "docs", {"val": 0.0004, "dec": 4}, "fit"),
    ("C30", "s25/results", {"val": 0.330, "dec": 3}, None),
    ("C30", "s25/results", {"val": 0.455, "dec": 3}, None),
    ("C31", "s25/results", {"val": -0.7423, "dec": 4}, None),
    ("C32", "s25/results", {"val": 2592, "dec": 0}, None),
    ("C33", "s25/results", {"lo": 0.78, "hi": 0.89}, "gap"),
    ("C33", "s25/results", {"lo": 78.0, "hi": 89.0}, "gap"),
    ("C34", "s25/results", {"val": 0.7529, "dec": 4}, None),
    ("C34", "s25/results", {"val": 0.8013, "dec": 4}, None),
    ("C34", "s25/results", {"val": 0.1127, "dec": 4}, None),
]

#: id, label, literal regexes, [(value, decimals), ...] for the stored-precision scan
CLAIMS = [
    ("C01", "3.2148 built-chain mean (rmsd_arm)", [r"3\.2148"], [(3.2148, 4)]),
    ("C02", "3.0483 point-cloud mean (rmsd_avg)", [r"3\.0483"], [(3.0483, 4)]),
    ("C03", "3.236 / 3.2355 repaired emission (rmsd_full)", [r"3\.2355", r"3\.236\b"],
     [(3.2355, 4)]),
    ("C04", "2.9610 benchmark full system", [r"2\.9610", r"2\.961\b"], [(2.961, 3)]),
    ("C05", "2.9507 benchmark shipped baseline", [r"2\.9507"], [(2.9507, 4)]),
    ("C06", "+0.0103 benchmark paired delta", [r"0\.0103"], [(0.0103, 4)]),
    ("C07", "1.711 / 1.7108 pool best", [r"1\.7108", r"1\.711\b"], [(1.7108, 4)]),
    ("C08", "3.425 / 3.4251", [r"3\.4251", r"3\.425\b"], [(3.4251, 4)]),
    ("C09", "4.065", [r"4\.065\b"], [(4.065, 3)]),
    ("C10", "3.770", [r"3\.770\b", r"3\.77\b"], [(3.770, 3)]),
    ("C11", "2.226 / 2.2261", [r"2\.2261", r"2\.226\b"], [(2.2261, 4)]),
    ("C12", "1.313", [r"1\.313\b"], [(1.313, 3)]),
    ("C13", "1.594", [r"1\.594\b"], [(1.594, 3)]),
    ("C14", "0.611", [r"0\.611\b"], [(0.611, 3)]),
    ("C15", "-2.15 A per gamma / -2.1496", [r"2\.1496", r"-2\.15\b"], [(-2.1496, 4)]),
    ("C16", "gamma = 0.0225", [r"0\.0225\b"], [(0.0225, 4)]),
    ("C17", "68% common-mode / 0.676", [r"0\.676\b", r"68%", r"68 per ?cent"], [(0.676, 3)]),
    ("C18", "0.902 nats", [r"0\.902\b"], [(0.902, 3)]),
    ("C19", "0.24x MDE", [r"0\.24 ?x\b", r"0\.24 MDE", r"0\.24 ?[x×] ?MDE"], []),
    ("C20", "1.18% worst deviation", [r"1\.18 ?%", r"1\.18 per ?cent", r"0\.0118\b"],
     [(0.0118, 4)]),
    ("C21", "5.6e-17", [r"5\.6\d*e-17"], [(5.6e-17, None)]),
    ("C22", "cos 1.000000000", [r"1\.000000000"], []),
    ("C23", "4.6e-10", [r"4\.6\d*e-10"], [(4.6e-10, None)]),
    ("C24", "0.656 / 0.655634 / 0.524 / 0.566586",
     [r"0\.655634", r"0\.656\b", r"0\.524\b", r"0\.566586"],
     [(0.655634, 6), (0.566586, 6), (0.524, 3)]),
    ("C25", "-0.649 / -0.252 / -0.047 / -0.311 log2 Var per qubit",
     [r"-0\.649\b", r"-0\.252\b", r"-0\.047\b", r"-0\.311\b"],
     [(-0.649, 3), (-0.252, 3), (-0.047, 3), (-0.311, 3)]),
    ("C26", "36.1 deg vs 36.4 deg", [r"36\.1\b", r"36\.4\b"], [(36.1, 1), (36.4, 1)]),
    ("C27", "+0.0004 / +0.0030 identity-leak prices",
     [r"\+0\.0004\b", r"\+0\.0030\b", r"0\.0004\b", r"0\.0030\b"],
     [(0.0004, 4), (0.0030, 4)]),
    ("C28", "2/60 and 4/126 self-copies", [r"2/60\b", r"4/126\b"], []),
    ("C29", "seven-configuration suite 3.058/3.132/3.215/3.253/3.674/3.755/3.881",
     [r"3\.058\b", r"3\.132\b", r"3\.215\b", r"3\.253\b", r"3\.674\b", r"3\.755\b",
      r"3\.881\b"],
     [(3.058, 3), (3.132, 3), (3.215, 3), (3.253, 3), (3.674, 3), (3.755, 3), (3.881, 3)]),
    ("C30", "+0.330 / +0.455", [r"\+0\.330\b", r"\+0\.455\b", r"0\.330\b", r"0\.455\b"],
     [(0.330, 3), (0.455, 3)]),
    ("C31", "rho -0.7423", [r"0\.7423"], [(-0.7423, 4)]),
    ("C32", "2,592 cells / 0 violations", [r"2,592", r"\b2592\b"], [(2592, 0)]),
    ("C33", "78-89% of the free-energy gap", [r"78-89", r"78 ?- ?89 ?%", r"78 to 89"], []),
    ("C34", "0.7529 / 0.8013 / 0.1127 mean |z|", [r"0\.7529", r"0\.8013", r"0\.1127"],
     [(0.7529, 4), (0.8013, 4), (0.1127, 4)]),
    ("C35", "355 passed / 13 skipped", [r"355 passed", r"13 skipped", r"\b355\b"], []),
]


def _skip(rel):
    parts = rel.split("/")
    if any(p in SKIP_PARTS for p in parts[:-1]):
        return True
    if rel in SKIP_FILES or any(s in rel for s in SKIP_SUBSTR):
        return True
    return False


def _match(x, m):
    if "dec" in m:
        return round(x, m["dec"]) == round(m["val"], m["dec"])
    if "rel" in m:
        return abs(x - m["val"]) <= m["rel"] * abs(m["val"])
    return m["lo"] <= x <= m["hi"]


def _walk_targeted(obj, path, m, keysub, out, cap=25):
    if isinstance(obj, dict):
        for k, v in obj.items():
            _walk_targeted(v, path + [str(k)], m, keysub, out, cap)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _walk_targeted(v, path + [str(i)], m, keysub, out, cap)
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)):
        if len(out) >= cap:
            return
        p = "/".join(path)
        if keysub and keysub.lower() not in p.lower():
            return
        try:
            if _match(float(obj), m):
                out.append({"path": p, "stored": obj})
        except (OverflowError, ValueError):
            pass


def targeted_pass():
    hits = []
    for cid, d, m, keysub in TARGETED:
        base = os.path.join(ROOT, d)
        found = []
        if os.path.isdir(base):
            for dp, dns, fns in os.walk(base):
                rel = os.path.relpath(dp, ROOT).replace(os.sep, "/")
                if any(p in SKIP_PARTS or p == "cache" for p in rel.split("/")):
                    dns[:] = []
                    continue
                for fn in sorted(fns):
                    if not fn.endswith(".json"):
                        continue
                    p = os.path.join(dp, fn)
                    r = os.path.relpath(p, ROOT).replace(os.sep, "/")
                    if _skip(r) or os.path.getsize(p) > MAX_JSON_BYTES:
                        continue
                    try:
                        with open(p, "r", encoding="utf-8") as fh:
                            obj = json.load(fh)
                    except Exception:                                    # noqa: BLE001
                        continue
                    loc = []
                    _walk_targeted(obj, [], m, keysub, loc)
                    for h in loc:
                        found.append({"file": r, **h})
        hits.append({"id": cid, "dir": d, "matcher": m, "key_filter": keysub,
                     "n_hits": len(found), "hits": found[:40]})
    return hits


def text_files():
    for dp, dns, fns in os.walk(ROOT):
        rel = os.path.relpath(dp, ROOT).replace(os.sep, "/")
        parts = [] if rel == "." else rel.split("/")
        if any(p in SKIP_PARTS for p in parts):
            dns[:] = []
            continue
        dns[:] = sorted(d for d in dns if d not in SKIP_PARTS)
        for fn in sorted(fns):
            ext = os.path.splitext(fn)[1].lower()
            if ext not in TEXT_EXT:
                continue
            p = os.path.join(dp, fn)
            r = os.path.relpath(p, ROOT).replace(os.sep, "/")
            if _skip(r):
                continue
            try:
                if os.path.getsize(p) > MAX_TEXT_BYTES:
                    continue
            except OSError:
                continue
            yield p, r


def literal_pass():
    pats = {cid: [re.compile(x) for x in lits] for cid, _, lits, _ in CLAIMS}
    hits = {cid: [] for cid, *_ in CLAIMS}
    for p, rel in text_files():
        try:
            with open(p, "r", encoding="utf-8", errors="replace") as fh:
                for ln, line in enumerate(fh, 1):
                    if len(line) > 4000:
                        line = line[:4000]
                    for cid, rs in pats.items():
                        for r in rs:
                            if r.search(line):
                                hits[cid].append({"file": rel, "line": ln,
                                                  "pattern": r.pattern,
                                                  "text": line.strip()[:220]})
                                break
        except OSError:
            continue
    return hits


def _walk(obj, path, out, targets, cap):
    if isinstance(obj, dict):
        for k, v in obj.items():
            _walk(v, path + [str(k)], out, targets, cap)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _walk(v, path + [str(i)], out, targets, cap)
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)):
        for cid, val, dec in targets:
            if len(out[cid]) >= cap:
                continue
            try:
                x = float(obj)
                if dec is None:                      # tiny magnitudes: 2% relative match
                    hit = abs(x - val) <= 0.02 * abs(val)
                else:
                    hit = round(x, dec) == round(val, dec)
                if hit:
                    out[cid].append({"path": "/".join(path), "stored": obj,
                                     "aggregate": is_aggregate(path)})
            except (OverflowError, ValueError):
                pass


def is_aggregate(path):
    """A leaf whose key path carries no array index: a summary value, not a per-row one."""
    return not any(p.isdigit() for p in path)


def stored_pass():
    targets = [(cid, v, d) for cid, _, _, nums in CLAIMS for v, d in nums]
    hits = {cid: [] for cid, *_ in CLAIMS}
    files = []
    for d in JSON_DIRS:
        base = os.path.join(ROOT, d)
        if not os.path.isdir(base):
            continue
        for dp, dns, fns in os.walk(base):
            rel = os.path.relpath(dp, ROOT).replace(os.sep, "/")
            if any(p in SKIP_PARTS or p == "cache" for p in rel.split("/")):
                dns[:] = []
                continue
            for fn in sorted(fns):
                if fn.endswith(".json"):
                    files.append(os.path.join(dp, fn))
    for p in sorted(set(files)):
        rel = os.path.relpath(p, ROOT).replace(os.sep, "/")
        if _skip(rel):
            continue
        try:
            if os.path.getsize(p) > MAX_JSON_BYTES:
                continue
            with open(p, "r", encoding="utf-8") as fh:
                obj = json.load(fh)
        except Exception:                                                # noqa: BLE001
            continue
        local = {cid: [] for cid, *_ in CLAIMS}
        _walk(obj, [], local, targets, cap=60)
        for cid, rows in local.items():
            for r in rows:
                hits[cid].append({"file": rel, **r})
    return hits


def main():
    try:                                   # Windows console defaults to cp1252
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                                    # noqa: BLE001
        pass
    lit = literal_pass()
    sto = stored_pass()
    payload = {"claims": []}
    lines = []
    for cid, label, lits, nums in CLAIMS:
        L, S = lit[cid], sto[cid]
        docs = [h for h in L if h["file"].endswith((".md", ".txt"))]
        arts = [h for h in L if not h["file"].endswith((".md", ".txt", ".py"))]
        code = [h for h in L if h["file"].endswith(".py")]
        agg = [h for h in S if h["aggregate"]]
        payload["claims"].append({
            "id": cid, "label": label, "literal_patterns": lits,
            "n_literal": len(L), "n_stored": len(S), "n_stored_aggregate": len(agg),
            "doc_hits": docs[:80], "artefact_hits": arts[:60], "code_hits": code[:30],
            "aggregate_hits": agg[:120], "stored_hits": S[:80]})
        lines.append(f"\n== {cid} {label}: literal {len(L)} (docs {len(docs)}, artefacts "
                     f"{len(arts)}, code {len(code)}), stored leaves {len(S)}, of which "
                     f"AGGREGATE (no array index in the key path) {len(agg)}")
        for h in docs[:12]:
            lines.append(f"   D {h['file']}:{h['line']}  {h['text'][:120]}")
        for h in agg[:14]:
            lines.append(f"   G {h['file']} :: {h['path']} = {h['stored']!r}")
        if not agg:
            for h in S[:4]:
                lines.append(f"   S {h['file']} :: {h['path']} = {h['stored']!r}")
    tg = targeted_pass()
    payload["targeted"] = tg
    lines.append("\n\n==== TARGETED LOOKUPS (the citing sprint's own results directory, "
                 "indexed paths included) ====")
    for t in tg:
        lines.append(f"\n-- {t['id']} in {t['dir']}  match {t['matcher']}  key~{t['key_filter']}"
                     f"  hits {t['n_hits']}")
        for h in t["hits"][:16]:
            lines.append(f"   T {h['file']} :: {h['path']} = {h['stored']!r}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1)
    txt = OUT[:-5] + ".txt"
    with open(txt, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {OUT} and {txt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
