#!/usr/bin/env python
"""s26/e_report_check.py -- every number in s26/REPORT.md Appendix B against its artefact.

Appendix B of the report is a table `| number | where used | artefact (file :: key, or
file:line) | stored value |`.  For every row this script

  1. finds every path in the artefact cell (backticked; `path :: key, key`, `path:line`,
     `path:a-b`, `tests/x.py::test_name`, brace and <placeholder> patterns) and checks that it
     exists on disk (a test name must exist as `def test_name` in its file);
  2. where the artefact is JSON, collects the leaves under the stated keys (all leaves when no
     key is stated; a list of numbers also contributes its mean, a container its length) and
     checks that every number quoted in the row is present within the rounding the quoted
     digits imply (half a unit in the last place; a percentage is tried as x and x/100);
  3. where the artefact is a markdown, text, log or source file, checks that the stated lines
     (or the stated ledger entries, or the whole file) contain the number;
  4. a stored-value cell of the form `derived: <expr> = <value>; ...` is evaluated (digits,
     + - * / ^ and parentheses only), every literal in the expression with four or more
     decimals must itself be found in the row's artefacts, and the results count as found.

Rule 1: the sealed benchmark record (`s9/final_report.json`, anything named `benchmark`) is
never opened; such a path is checked for existence only and the row says so.

Verdict per row: PASS, MISSING PATH, or NUMBER NOT FOUND (with the numbers listed).  `--style`
(default on) also runs the style check over the report's own prose (everything before Appendix C,
which reproduces the ledger verbatim): the six banned words, U+2014 and U+2013, and any sentence
that contrasts two RMSD-like numbers without a basis label on the sentence or its paragraph.  Output: a table on stdout and `s26/results/e_report_check.json` with provenance.

    python s26/e_report_check.py                # table + JSON, exit 1 on any failure
    python s26/e_report_check.py --quiet        # one line per failing row only
    python s26/e_report_check.py --no-style     # Appendix B only
"""
from __future__ import annotations

import argparse
import fnmatch
import glob
import hashlib
import io
import json
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

REPORT = os.path.join(HERE, "REPORT.md")
OUT = os.path.join(HERE, "results", "e_report_check.json")

NEVER_OPEN = ("benchmark", "bench60", "s9/final_report.json", "s9/final_synth.json",
              "s9/final_cache", "results/monomer_manifest.json")
TEXT_EXT = (".md", ".txt", ".py", ".log", ".csv", ".sh", ".bat", ".toml", ".cfg", ".ini")
BANNED = ("genuinely", "honestly", "leverage", "robust", "delve", "underscore")
BASIS = ("point cloud", "built chain", "relaxed chain", "single window", "same basis",
         "rmsd_avg", "rmsd_arm", "rmsd_full", "rmsd_fit", "selection basis", "s8 instrument",
         "lam = 0 chain", "basis not stated", "basis unstated", "both bases", "on both bases",
         "point-cloud", "built-chain", "relaxed-chain", "single-window")
CUES = (" against ", " minus ", " vs ", " versus ", " beats ", " worse ", " better ", " gap ",
        " difference ", " compared ", " than ", " costs ", " cost ", " buys ", " over ")

NUM_TOKEN = re.compile(
    r"(?<![\w.])([-+]?\d+(?:,\d{3})*(?:\.\d+)?(?:e[-+]?\d+)?)"
    r"(%|x|th|st|nd|rd|W|L|k|ms|s|GB|MB|deg|A)?(?![\w.]|\.\d)")
HEX_TOKEN = re.compile(r"(?<![\w])([0-9a-f]{6,64})(?![\w])")
LEDGER_TOKEN = re.compile(r"\bL(\d+)(?:-([A-Z]))?\b|\bL-B\b")


# ----------------------------------------------------------------------------- report parsing

def appendix_b_rows(text):
    """Return (row_number, line_number, [cells]) for every data row of Appendix B."""
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("## APPENDIX B"))
    end = next((i for i, l in enumerate(lines) if l.startswith("## APPENDIX C")), len(lines))
    rows = []
    k = 0
    for i in range(start, end):
        l = lines[i]
        if not l.startswith("| "):
            continue
        cells = [c.strip() for c in l.strip().strip("|").split(" | ")]
        if len(cells) < 4 or cells[0] == "number" or set(cells[0]) <= set("-"):
            continue
        k += 1
        rows.append((k, i + 1, cells))
    return rows


def number_tokens(cell):
    """Quoted numbers of a cell as (text, value, kind); kind in num / pct / hex."""
    s = cell
    s = s.replace("..", " to ").replace("/", " ")
    s = re.sub(r"\[|\]|\{|\}|\(|\)|:|;|,(?!\d{3}\b)", " ", s)
    out = []
    seen = set()
    for m in NUM_TOKEN.finditer(s):
        tok, suf = m.group(1), m.group(2)
        raw = tok.replace(",", "")
        try:
            v = float(raw)
        except ValueError:
            continue
        kind = "pct" if suf == "%" else "num"
        key = (raw, kind)
        if key in seen:
            continue
        seen.add(key)
        out.append({"text": raw, "value": v, "kind": kind, "tol": tolerance(raw)})
    for m in HEX_TOKEN.finditer(cell):
        h = m.group(1)
        if re.search(r"[a-f]", h) and re.search(r"\d", h):
            out.append({"text": h, "value": None, "kind": "hex", "tol": None})
    for m in re.finditer(r"([0-9a-f]{6,})\.\.\.([0-9a-f]{3,})", cell):
        out.append({"text": m.group(1), "value": None, "kind": "hex", "tol": None})
        out.append({"text": m.group(2), "value": None, "kind": "hex", "tol": None})
    return out


def tolerance(raw):
    m = re.match(r"[-+]?(\d+)(?:\.(\d+))?(?:e([-+]?\d+))?$", raw)
    ip, dp, ex = m.groups()
    d = len(dp or "")
    e = int(ex or 0)
    return 0.5 * 10.0 ** (e - d)


# ------------------------------------------------------------------------ citation parsing

def expand_braces(s):
    m = re.search(r"\{([^{}]*,[^{}]*)\}", s)
    if not m:
        return [s]
    out = []
    for alt in m.group(1).split(","):
        out.extend(expand_braces(s[:m.start()] + alt.strip() + s[m.end():]))
    return out


def parse_citations(cell):
    """Every backticked span of the artefact cell, classified."""
    cites = []
    last_dir = None
    for m in re.finditer(r"`([^`]+)`", cell):
        span = m.group(1).strip()
        keys = []
        test = None
        line_a = line_b = None
        path = span
        if " :: " in span:
            path, k = span.split(" :: ", 1)
            keys = [x.strip() for x in k.split(", ") if x.strip()]
        elif "::" in span and re.search(r"\.py::", span):
            path, test = span.split("::", 1)
        else:
            lm = re.match(r"^(.+?):(\d+)(?:-(\d+))?$", span)
            if lm:
                path, line_a, line_b = lm.group(1), int(lm.group(2)), int(lm.group(3) or lm.group(2))
        path = path.strip()
        if not (("/" in path) or re.search(r"\.(json|md|txt|py|log|csv|npz|pt|pdb|sh|bat|toml)$", path)):
            continue                                        # a bare name such as DEFAULT_WEIGHTS
        if "/" not in path:
            if os.path.exists(os.path.join(ROOT, path)):
                pass                                        # a root file such as README.md
            elif last_dir and os.path.exists(os.path.join(ROOT, last_dir, path)):
                path = last_dir + "/" + path                # `a/x.json`, `y.json` (same directory)
            else:
                hits = [h for h in glob.glob(os.path.join(ROOT, "**", path), recursive=True)
                        if "node_modules" not in h]
                if hits:
                    path = os.path.relpath(hits[0], ROOT).replace("\\", "/")
        last_dir = path.rsplit("/", 1)[0] if "/" in path else last_dir
        cites.append({"span": span, "path": path, "keys": keys, "test": test,
                      "lines": (line_a, line_b) if line_a else None, "pos": m.start()})
    return cites


def ledger_refs(cell, cites):
    """Map ledger paths to the entry tokens (L17, L10-A, L-B) named after them in the cell."""
    anchors = []
    for c in cites:
        if c["path"].endswith("LEDGER.md"):
            anchors.append((c["pos"], c["path"]))
    for m in re.finditer(r"\bS(\d\d) ledger\b", cell):
        anchors.append((m.start(), "s%s/LEDGER.md" % m.group(1)))
    anchors.sort()
    refs = {}
    for m in LEDGER_TOKEN.finditer(cell):
        tok = m.group(0)
        prior = [a for a in anchors if a[0] < m.start()]
        if not prior:
            continue
        refs.setdefault(prior[-1][1], []).append(tok)
    for _, path in anchors:
        refs.setdefault(path, [])
    return refs


# ------------------------------------------------------------------------- artefact access

def never_open(path):
    p = path.replace("\\", "/")
    return any(s in p for s in NEVER_OPEN)


def resolve_paths(path):
    """Expand braces and <placeholders>; return the list of existing relative paths."""
    out = []
    for cand in expand_braces(path):
        pat = re.sub(r"<[^>/]+>", "*", cand)
        if any(ch in pat for ch in "*?"):
            hits = sorted(glob.glob(os.path.join(ROOT, pat)))
            out.extend(os.path.relpath(h, ROOT).replace("\\", "/") for h in hits)
        elif os.path.exists(os.path.join(ROOT, cand)):
            out.append(cand)
    return out


_json_cache = {}


def load_json(rel):
    if rel not in _json_cache:
        with io.open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
            _json_cache[rel] = json.load(fh)
    return _json_cache[rel]


_text_cache = {}


def load_text(rel):
    if rel not in _text_cache:
        with io.open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace") as fh:
            _text_cache[rel] = fh.read()
    return _text_cache[rel]


def seg_match(name, key):
    if name == key:
        return True
    pat = re.sub(r"<[^>]+>", "*", name)
    return any(ch in pat for ch in "*?") and fnmatch.fnmatchcase(key, pat)


def resolve_key(obj, key):
    """Objects reached from `obj` by a key path such as rows[*]/MASS1.0, exact[*].terms, results/*."""
    cur = [obj]
    for seg in key.split("/"):
        nxt = []
        for o in cur:
            nxt.extend(_apply_segment(o, seg))
        cur = nxt
        if not cur:
            break
    return cur


def _apply_segment(o, seg):
    # 1. literal key (MASS1.0, n7_L2, ...)
    if isinstance(o, dict) and seg in o:
        return [o[seg]]
    if isinstance(o, list) and re.fullmatch(r"\d+", seg) and int(seg) < len(o):
        return [o[int(seg)]]
    # 2. name[idx][idx].attr.attr, where idx is *, an integer, or a field=value filter over a
    #    list of dicts (`concentration[model=legacy]`)
    m = re.fullmatch(r"([^\[\].]*)((?:\[(?:\*|\d+|[\w.]+=[^\]]+)\])*)(?:\.(.+))?", seg)
    if not m:
        return []
    name, idx, attrs = m.group(1), m.group(2), m.group(3)
    cur = [o]
    if name:
        nxt = []
        for c in cur:
            if isinstance(c, dict):
                nxt.extend(v for k, v in c.items() if seg_match(name, k))
            elif isinstance(c, list):
                for e in c:                                  # implicit map over a list
                    if isinstance(e, dict):
                        nxt.extend(v for k, v in e.items() if seg_match(name, k))
        cur = nxt
    for ix in re.findall(r"\[([^\]]+)\]", idx):
        nxt = []
        for c in cur:
            if "=" in ix and not ix.isdigit():
                field, val = ix.split("=", 1)
                if isinstance(c, list):
                    nxt.append([e for e in c if isinstance(e, dict) and str(e.get(field)) == val])
                elif isinstance(c, dict):
                    nxt.append([e for e in c.values() if isinstance(e, dict) and str(e.get(field)) == val])
            elif isinstance(c, list):
                nxt.extend(c if ix == "*" else ([c[int(ix)]] if int(ix) < len(c) else []))
            elif isinstance(c, dict):
                nxt.extend(c.values() if ix == "*" else [])
        cur = nxt
    if attrs:
        for a in attrs.split("."):
            nxt = []
            for c in cur:
                if isinstance(c, dict) and a in c:
                    nxt.append(c[a])
                elif isinstance(c, list):
                    nxt.extend(e[a] for e in c if isinstance(e, dict) and a in e)
            cur = nxt
    return cur


def collect(obj, label, nums, strs, cap=3_000_000):
    """Numeric leaves, string leaves, list means and container lengths under obj."""
    stack = [(obj, label)]
    while stack and len(nums) < cap:
        o, lab = stack.pop()
        if isinstance(o, bool):
            continue
        if isinstance(o, (int, float)):
            if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
                continue
            nums.append((float(o), lab))
        elif isinstance(o, str):
            strs.append((o, lab))
            for m in NUM_TOKEN.finditer(o):                  # "3.339e-01", "n=7 qubits, layers=3"
                try:
                    nums.append((float(m.group(1).replace(",", "")), lab + "#str"))
                except ValueError:
                    pass
        elif isinstance(o, dict):
            nums.append((float(len(o)), lab + "#len"))
            for k, v in o.items():
                stack.append((v, lab + "/" + str(k)))
        elif isinstance(o, list):
            nums.append((float(len(o)), lab + "#len"))
            vals = [x for x in o if isinstance(x, (int, float)) and not isinstance(x, bool)]
            if vals and len(vals) == len(o):
                nums.append((float(sum(vals)) / len(vals), lab + "#mean"))
                for i, v in enumerate(o):
                    if not (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
                        nums.append((float(v), lab + "[%d]" % i))
            else:
                for i, v in enumerate(o):
                    stack.append((v, lab + "[%d]" % i))


def norm_text(t):
    t = t.replace("−", "-").replace("–", "-").replace("—", "-")
    for _ in range(3):
        t = re.sub(r"(\d),(\d{3})(?!\d)", r"\1\2", t)
    return t


def ledger_entry_text(text, tok):
    """The text of ledger entry `tok` (L17, L10-A, L-B), or None."""
    lines = text.splitlines()
    pat = re.compile(r"^## %s(?![\w-])" % re.escape(tok))
    for i, l in enumerate(lines):
        if pat.match(l):
            j = i + 1
            while j < len(lines) and not lines[j].startswith("## "):
                j += 1
            return "\n".join(lines[i:j])
    return None


# ----------------------------------------------------------------------------- matching

def num_found(tok, nums):
    v, tol = abs(tok["value"]), tok["tol"]
    nums = [(abs(x), lab) for x, lab in nums]
    cands = [(v, tol)]
    if tok["kind"] == "pct":
        cands.append((v / 100.0, tol / 100.0))
    for target, t in cands:
        for x, lab in nums:
            if abs(x - target) <= t + 1e-12 * max(1.0, abs(target)):
                return lab
    return None


def text_found(tok, text):
    t = norm_text(text)
    if tok["kind"] == "hex":
        return tok["text"] in t
    raw = tok["text"].lstrip("+-")
    variants = {raw}
    if "e" in raw:
        mant, ex = raw.split("e")
        variants.add(mant + "e" + str(int(ex)))
        variants.add(mant + "e" + ("+%d" % int(ex) if int(ex) >= 0 else "%d" % int(ex)))
        variants.add(mant + "e" + ("%02d" % int(ex) if int(ex) >= 0 else "-%02d" % -int(ex)))
    if "." in raw:
        variants.add(raw.rstrip("0").rstrip(".") if raw.rstrip("0").endswith(".") is False else raw)
    for var in variants:
        if re.search(r"(?<![\d])" + re.escape(var) + r"(?![\d])", t):
            return True
    return False


_text_nums_cache = {}


def text_numbers(seg, label):
    """Every number written in a text segment, for tolerance matching (2.1496 matches -2.15)."""
    if label not in _text_nums_cache:
        out = []
        for m in NUM_TOKEN.finditer(norm_text(seg)):
            try:
                out.append((abs(float(m.group(1).replace(",", ""))), label + "#text"))
            except ValueError:
                pass
        _text_nums_cache[label] = out
    return _text_nums_cache[label]


SAFE_EXPR = re.compile(r"^[\d\s.eE+\-*/()^]+$")


def eval_derived(cell):
    """`derived: expr = value; expr = value` -> list of (expr, value, ok, literals)."""
    body = cell.split("derived:", 1)[1]
    out = []
    for part in body.split(";"):
        if "=" not in part:
            continue
        expr, val = part.rsplit("=", 1)
        expr, val = expr.strip(), val.strip()
        if not SAFE_EXPR.match(expr):
            out.append({"expr": expr, "value": val, "ok": False, "literals": []})
            continue
        try:
            got = eval(expr.replace("^", "**"), {"__builtins__": {}}, {})  # noqa: S307 (digits/operators only)
            want = float(val)
            ok = abs(got - want) <= tolerance(val) + 1e-12
        except Exception:
            got, ok = None, False
        lits = [x for x in re.findall(r"\d+\.\d{4,}", expr)]
        out.append({"expr": expr, "value": val, "got": got, "ok": bool(ok), "literals": lits})
    return out


# --------------------------------------------------------------------------- one row

def check_row(k, lineno, cells):
    num_cell, where, art_cell, stored_cell = cells[0], cells[1], cells[2], cells[3]
    toks = number_tokens(num_cell)
    if not stored_cell.lower().startswith(("as ", "derived")):
        for t in number_tokens(stored_cell):
            t["from_stored"] = True
            toks.append(t)
    cites = parse_citations(art_cell)
    lrefs = ledger_refs(art_cell, cites)
    notes, missing_paths = [], []
    nums, strs, texts = [], [], []      # candidates across every source of the row
    sources = []
    for c in cites:
        rels = resolve_paths(c["path"])
        if not rels:
            missing_paths.append(c["path"])
            continue
        if c["test"]:
            body = load_text(rels[0])
            if not re.search(r"^def %s\b" % re.escape(c["test"]), body, re.M):
                missing_paths.append(c["path"] + "::" + c["test"])
                continue
            sources.append(rels[0] + "::" + c["test"])
            continue
        for rel in rels[:12]:
            if never_open(rel):
                notes.append("not opened (Rule 1): " + rel)
                sources.append(rel + " [existence only]")
                continue
            ext = os.path.splitext(rel)[1].lower()
            if ext == ".json":
                try:
                    obj = load_json(rel)
                except Exception as e:        # noqa: BLE001
                    notes.append("unreadable json %s: %s" % (rel, e))
                    continue
                if c["keys"]:
                    for key in c["keys"]:
                        hits = resolve_key(obj, key)
                        if not hits:
                            notes.append("key not found: %s :: %s" % (rel, key))
                        for h in hits:
                            collect(h, rel + "::" + key, nums, strs)
                        vals = [float(h) for h in hits
                                if isinstance(h, (int, float)) and not isinstance(h, bool)]
                        if len(vals) > 1:                    # rows[*]/x: the mean and the count
                            nums.append((sum(vals) / len(vals), rel + "::" + key + "#mean"))
                            nums.append((float(len(vals)), rel + "::" + key + "#count"))
                    sources.append(rel + " :: " + ", ".join(c["keys"]))
                else:
                    collect(obj, rel, nums, strs)
                    sources.append(rel + " [all leaves]")
            elif ext in TEXT_EXT:
                body = load_text(rel)
                if c["lines"]:
                    a, b = c["lines"]
                    seg = "\n".join(body.splitlines()[a - 1:b])
                    texts.append((seg, "%s:%d-%d" % (rel, a, b)))
                    sources.append("%s:%d-%d" % (rel, a, b))
                elif rel in lrefs and lrefs[rel]:
                    for tok in lrefs[rel]:
                        ent = ledger_entry_text(body, tok)
                        if ent is None:
                            notes.append("ledger entry not found: %s %s" % (rel, tok))
                        else:
                            texts.append((ent, rel + " " + tok))
                    sources.append(rel + " " + " ".join(lrefs[rel]))
                else:
                    texts.append((body, rel))
                    sources.append(rel + " [whole file]")
            else:
                sources.append(rel + " [existence only]")
    # ledger prose references with no backticked path ("S26 ledger L14")
    for path, toklist in lrefs.items():
        if any(path == c["path"] for c in cites):
            continue
        if not os.path.exists(os.path.join(ROOT, path)):
            missing_paths.append(path)
            continue
        body = load_text(path)
        for tok in toklist:
            ent = ledger_entry_text(body, tok)
            if ent is None:
                notes.append("ledger entry not found: %s %s" % (path, tok))
            else:
                texts.append((ent, path + " " + tok))
        sources.append(path + " " + " ".join(toklist))
    # derived expressions
    derived_ok_values = []
    if stored_cell.lower().startswith("derived:"):
        for d in eval_derived(stored_cell):
            lit_missing = [x for x in d["literals"]
                           if num_found({"value": float(x), "tol": tolerance(x), "kind": "num"}, nums) is None
                           and not any(text_found({"text": x, "kind": "num"}, t) for t, _ in texts)]
            if d["ok"] and not lit_missing:
                derived_ok_values.append(d)
                notes.append("derived ok: %s = %s" % (d["expr"], d["value"]))
            else:
                notes.append("derived FAILED: %s = %s (literals missing: %s)"
                             % (d["expr"], d["value"], ", ".join(lit_missing) or "none"))
    # numbers
    found, missing = [], []
    for t in toks:
        if t["kind"] == "hex":
            hit = any(t["text"] in s for s, _ in strs) or any(text_found(t, x) for x, _ in texts)
            (found if hit else missing).append(t["text"])
            continue
        lab = num_found(t, nums)
        if lab is None:
            for seg, slab in texts:
                if text_found(t, seg) or num_found(t, text_numbers(seg, slab)) is not None:
                    lab = slab
                    break
        if lab is None:
            for d in derived_ok_values:
                targets = [(t["value"], t["tol"])]
                if t["kind"] == "pct":
                    targets.append((t["value"] / 100.0, t["tol"] / 100.0))
                if any(abs(float(d["value"]) - tv) <= tt + 1e-12 for tv, tt in targets):
                    lab = "derived"
                    break
        (found if lab else missing).append(t["text"] + ("%" if t["kind"] == "pct" else ""))
    if any("not opened (Rule 1)" in n for n in notes) and not sources.count(""):
        pass
    if missing_paths:
        status = "MISSING PATH"
    elif missing and not all("not opened (Rule 1)" in n for n in notes if n) or (missing and not notes):
        status = "NUMBER NOT FOUND"
    elif missing:
        status = "PASS"
        notes.append("numbers unverifiable behind Rule 1: " + ", ".join(missing))
    else:
        status = "PASS"
    return {"row": k, "line": lineno, "where": where, "numbers": num_cell[:160],
            "status": status, "missing_paths": missing_paths, "missing_numbers": missing,
            "found_numbers": found, "sources": sources, "notes": notes}


# ------------------------------------------------------------------------------ style

def style_check(text):
    lines = text.splitlines()
    banned = []
    for i, l in enumerate(lines, 1):
        for w in BANNED:
            for m in re.finditer(r"\b%s\w*" % w, l, re.I):
                banned.append({"line": i, "word": m.group(0)})
    dashes = [{"line": i, "char": "U+2014" if "—" in l else "U+2013"}
              for i, l in enumerate(lines, 1) if "—" in l or "–" in l]
    # contrasts: a sentence about RMSD (an Angstrom or RMSD marker) with two or more numbers of
    # RMSD size, a contrast cue, and no basis label; sentences about other quantities are skipped
    rmsd_num = re.compile(r"(?<![\w.])[-+]?\d\.\d{3,4}(?![\d])")
    rmsd_marker = re.compile(r"\bRMSD\b|\d A\b|\bA\)|\bAngstrom\b|\bemits?\b|\bendpoint\b", re.I)
    not_rmsd = ("radius of gyration", " rg ", "free energy", "kl ", "variance", "cosine", "pauli",
                "weight", "correlation", "corr(", "participation", "schmidt", "bond dimension",
                "gradient", "relative error", "entropy", "decay", "walsh", "log2", "kcal",
                "percentile of a uniform", "spearman", "rho(", "torsion error", "degrees",
                "theorem", "agreement")
    flagged = []
    para, start = [], 1
    for i, l in enumerate(lines + [""], 1):
        if l.strip():
            if not para:
                start = i
            para.append(l)
            continue
        if para:
            block = " ".join(x.strip() for x in para)
            if not block.startswith("|"):
                lowered = block.lower()
                para_basis = any(b in lowered for b in BASIS)
                for sent in re.split(r"(?<=[.;])\s+(?=[A-Z(`])", block):
                    nums = rmsd_num.findall(sent)
                    low = " " + sent.lower() + " "
                    if len(nums) >= 2 and any(c in low for c in CUES) and rmsd_marker.search(sent) \
                            and not any(w in low for w in not_rmsd) \
                            and not any(b in low for b in BASIS) and not para_basis:
                        flagged.append({"line": start, "sentence": sent[:220]})
            para = []
    return {"banned": banned, "dashes": dashes, "contrasts_without_basis": flagged}


# ------------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default=REPORT)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--no-style", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:   # noqa: BLE001
        pass
    t0 = time.time()
    text = io.open(a.report, encoding="utf-8").read()
    rows = appendix_b_rows(text)
    results = [check_row(k, ln, cells) for k, ln, cells in rows]
    counts = {s: sum(1 for r in results if r["status"] == s)
              for s in ("PASS", "MISSING PATH", "NUMBER NOT FOUND")}
    # the style rules apply to the report's own prose; Appendix C reproduces s26/LEDGER.md
    # verbatim and is exempt (its two banned-word occurrences are declared in its preamble)
    own_text = text.split("
## APPENDIX C")[0]
    style = None if a.no_style else style_check(own_text)
    if not a.quiet:
        print("%-4s %-17s %-6s %s" % ("row", "status", "line", "numbers / problem"))
    for r in results:
        if a.quiet and r["status"] == "PASS":
            continue
        prob = ""
        if r["missing_paths"]:
            prob = "missing: " + ", ".join(r["missing_paths"])
        elif r["missing_numbers"]:
            prob = "not found: " + ", ".join(r["missing_numbers"])
        print("%-4d %-17s %-6d %s" % (r["row"], r["status"], r["line"],
                                      (prob or r["numbers"][:90])))
        for n in r["notes"]:
            if not a.quiet or "FAILED" in n or "not found" in n:
                print("     . " + n)
    print("rows %d: PASS %d, MISSING PATH %d, NUMBER NOT FOUND %d"
          % (len(results), counts["PASS"], counts["MISSING PATH"], counts["NUMBER NOT FOUND"]))
    rc = int(counts["MISSING PATH"] + counts["NUMBER NOT FOUND"] > 0)
    if style is not None:
        print("style: banned words %d, en/em dashes %d, contrasts without a basis label %d"
              % (len(style["banned"]), len(style["dashes"]), len(style["contrasts_without_basis"])))
        for b in style["banned"]:
            print("     . banned word line %d: %s" % (b["line"], b["word"]))
        for d in style["dashes"]:
            print("     . dash line %d: %s" % (d["line"], d["char"]))
        for f in style["contrasts_without_basis"]:
            print("     . contrast line %d: %s" % (f["line"], f["sentence"]))
        rc |= int(len(style["banned"]) + len(style["dashes"]) > 0)
    out = {
        "kind": "report check: every Appendix B number of s26/REPORT.md against its artefact",
        "lane": "E", "sprint": 26,
        "report": os.path.relpath(a.report, ROOT).replace("\\", "/"),
        "report_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "counts": counts, "style": style, "wall_s": round(time.time() - t0, 3),
        "rows": results,
    }
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    try:
        from s24 import stats_lib
        stats_lib.save_atomic(a.out, out, complete_keys=("row", "status", "sources"),
                              rows=results, n_expected=len(rows), module_file=__file__)
    except Exception as e:      # noqa: BLE001
        out["provenance"] = {"error": "stats_lib.save_atomic unavailable: %s" % e}
        with io.open(a.out, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1)
    print("wrote", os.path.relpath(a.out, ROOT).replace("\\", "/"))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
