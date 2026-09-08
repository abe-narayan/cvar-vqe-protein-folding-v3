"""SHIFT agent -- BMRB access layer.

Fetches and caches BMRB depositions for the 126 development targets and measures what
backbone chemical-shift information is ACTUALLY available for each.

Nothing here reads native coordinates, native torsions or native RMSD.  Chemical shifts
are an experimental observable; see s14/shift_FINDINGS.md section 0 for the leakage
framing (the deposited PDB coordinates were solved USING these shifts plus NOEs, so any
arm driven by them is NMR-restrained structure determination in its own column).

    python -m s14.shift_bmrb          # build the cache + availability table
"""
from __future__ import annotations

import os
import sys
import json
import time
import urllib.request
import urllib.error

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

CACHE = os.path.join(ROOT, "s14", "cache", "bmrb")
RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

API = "https://api.bmrb.io/v2"
HDR = {"Application": "s14-shift-agent"}

# The nuclei a TALOS-class method consumes.  TALOS-N's input set is
# HN, N, HA, CA, CB, C' (Shen & Bax 2013).
BACKBONE = ("H", "N", "HA", "CA", "CB", "C")
# BMRB atom-name variants that map onto each canonical slot.
ATOM_ALIAS = {
    "H": "H", "HN": "H",
    "N": "N",
    "HA": "HA", "HA2": "HA", "HA3": "HA", "HA1": "HA",
    "CA": "CA",
    "CB": "CB",
    "C": "C", "C'": "C", "CO": "C",
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


# --------------------------------------------------------------------------- fetching
def _get(url, timeout=90, tries=3):
    last = None
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=HDR)
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.load(fh)
        except Exception as exc:                                  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (a + 1))
    raise RuntimeError(f"GET failed {url}: {last!r}")


def cached(name, fn):
    """Disk-cached JSON.  `name` is a bare file stem under s14/cache/bmrb/."""
    p = os.path.join(CACHE, name + ".json")
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    obj = fn()
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, p)
    return obj


def pdb_bmrb_map(match_type="exact"):
    d = cached(f"map_pdb_bmrb_{match_type}",
               lambda: _get(f"{API}/mappings/pdb/bmrb?match_type={match_type}"))
    return {r["pdb_id"].upper(): [str(x) for x in r["bmrb_ids"]] for r in d}


def fasta_search(seq):
    """BMRB's own BLAST over deposited polymer sequences."""
    return cached("fasta_" + seq,
                  lambda: _get(f"{API}/search/fasta/{seq}?type=polymer"))


def entry_entities(bid):
    d = cached(f"entity_{bid}",
               lambda: _get(f"{API}/entry/{bid}?saveframe_category=entity"))
    out = []
    for sf in d.get(str(bid), {}).get("entity", []) or []:
        tags = dict(sf.get("tags", []))
        out.append({
            "id": tags.get("ID"),
            "name": tags.get("Name"),
            "seq": (tags.get("Polymer_seq_one_letter_code") or "").replace("\n", "").strip(),
            "type": tags.get("Polymer_type"),
            "nmono": tags.get("Number_of_monomers"),
        })
    return out


def entry_shifts(bid):
    """-> list of rows {entity, seq_id, comp, atom, val, list_id}."""
    d = cached(f"cs_{bid}",
               lambda: _get(f"{API}/entry/{bid}?saveframe_category=assigned_chemical_shifts"))
    rows = []
    for sf in d.get(str(bid), {}).get("assigned_chemical_shifts", []) or []:
        for lp in sf.get("loops", []):
            if lp.get("category") != "_Atom_chem_shift":
                continue
            tg = {t: i for i, t in enumerate(lp["tags"])}
            for r in lp["data"]:
                def g(k):
                    i = tg.get(k)
                    return None if i is None else (None if r[i] in (".", "?", "") else r[i])
                try:
                    val = float(g("Val")) if g("Val") is not None else None
                except ValueError:
                    val = None
                sid = g("Seq_ID") or g("Comp_index_ID")
                if sid is None or val is None:
                    continue
                rows.append({
                    "entity": g("Entity_ID") or "1",
                    "seq_id": int(sid),
                    "comp": (g("Comp_ID") or "").upper(),
                    "atom": (g("Atom_ID") or "").upper(),
                    "val": val,
                    "list_id": g("Assigned_chem_shift_list_ID") or "1",
                })
    return rows


def entry_meta(bid):
    """Sample conditions / experiment metadata, for the solution-state caveat."""
    def _f():
        try:
            return _get(f"{API}/entry/{bid}?saveframe_category=sample_conditions")
        except Exception:                                          # noqa: BLE001
            return {}
    d = cached(f"meta_{bid}", _f)
    out = []
    for sf in d.get(str(bid), {}).get("sample_conditions", []) or []:
        for lp in sf.get("loops", []):
            if lp.get("category") != "_Sample_condition_variable":
                continue
            tg = {t: i for i, t in enumerate(lp["tags"])}
            cond = {}
            for r in lp["data"]:
                cond[r[tg["Type"]]] = r[tg["Val"]]
            out.append(cond)
    return out


# --------------------------------------------------------------------- sequence mapping
def locate(target_seq, entity_seq):
    """Map each target residue index -> entity Seq_ID (1-based), or None.

    Exact substring first (the overwhelmingly common case: our targets are windows of, or
    equal to, the deposited polymer).  Falls back to a simple ungapped best-offset scan.
    Returns (offset0, n_match, mapping list) or None if the best ungapped match has
    < 70 % identity.
    """
    t, e = target_seq.upper(), entity_seq.upper()
    if not e:
        return None
    j = e.find(t)
    if j >= 0:
        return {"off": j, "ident": 1.0, "map": [j + i + 1 for i in range(len(t))]}
    best, bo = -1, None
    for o in range(-len(t) + 1, len(e)):
        m = sum(1 for i in range(len(t))
                if 0 <= o + i < len(e) and e[o + i] == t[i])
        if m > best:
            best, bo = m, o
    ident = best / len(t)
    if ident < 0.70:
        return None
    return {"off": bo, "ident": ident,
            "map": [(bo + i + 1) if 0 <= bo + i < len(e) else None for i in range(len(t))]}


# ------------------------------------------------------------------------- per-target
def per_residue_nuclei(rows, entity_id, seq_ids):
    """seq_ids: list (len n) of BMRB Seq_ID or None.  -> list of sets of canonical nuclei."""
    have = {}
    for r in rows:
        if str(r["entity"]) != str(entity_id):
            continue
        a = ATOM_ALIAS.get(r["atom"])
        if a is None:
            continue
        have.setdefault(r["seq_id"], set()).add(a)
    return [set() if s is None else set(have.get(s, ())) for s in seq_ids]


def talos_predictable(nuclei, min_shifts=3, min_neigh=2):
    """TALOS-N shift-completeness gate (Shen & Bax 2013):
    'If at least two of the three residues [i-1, i, i+1] have at least three chemical
    shifts, the center residue is considered to be predictable.'
    """
    n = len(nuclei)
    ok = [len(s) >= min_shifts for s in nuclei]
    out = []
    for i in range(n):
        c = sum(1 for j in (i - 1, i, i + 1) if 0 <= j < n and ok[j])
        out.append(c >= min_neigh)
    return out


def main():
    from s12 import instrument as I

    tg = I.targets()
    fail18 = set(I.FAIL18)
    ex = pdb_bmrb_map("exact")
    au = pdb_bmrb_map("author")
    try:
        with open(os.path.join(ROOT, "s12", "results", "lit_bmrb_coverage.json")) as fh:
            s12map = json.load(fh)
    except Exception:                                              # noqa: BLE001
        s12map = {}

    out = {}
    for k, t in enumerate(tg):
        pdb, seq, n = t["pdb"], t["seq"], t["n"]
        cands, why = [], {}
        for src, ids in (("pdb_exact", ex.get(pdb, [])),
                         ("pdb_author", au.get(pdb, [])),
                         ("s12_rcsb", (s12map.get(pdb) or {}).get("bmrb", []))):
            for b in ids:
                b = str(b)
                if b not in cands:
                    cands.append(b)
                why.setdefault(b, []).append(src)
        # sequence search -- finds depositions of the same peptide not cross-linked to the PDB
        fas = []
        try:
            fas = fasta_search(seq)
        except Exception as exc:                                   # noqa: BLE001
            print(f"  [{pdb}] fasta failed {exc!r}")
        for h in fas:
            if float(h.get("percent_id", 0)) >= 99.0 and int(h.get("alignment_length", 0)) >= n:
                b = str(h["entry_id"])
                if b not in cands:
                    cands.append(b)
                why.setdefault(b, []).append("fasta100")

        recs = []
        for b in cands:
            try:
                ents = entry_entities(b)
                rows = entry_shifts(b)
            except Exception as exc:                               # noqa: BLE001
                print(f"  [{pdb}] bmrb {b} failed {exc!r}")
                continue
            for e in ents:
                if e["type"] and "polypeptide" not in e["type"]:
                    continue
                loc = locate(seq, e["seq"])
                if loc is None:
                    continue
                nuc = per_residue_nuclei(rows, e["id"], loc["map"])
                pred = talos_predictable(nuc)
                recs.append({
                    "bmrb": b, "entity": e["id"], "src": why.get(b, []),
                    "ident": loc["ident"], "off": loc["off"],
                    "entity_len": len(e["seq"]), "entity_name": e["name"],
                    "n_shift_rows": len(rows),
                    "nuclei": ["".join(sorted(s)) for s in nuc],
                    "per_nuc": {a: [int(a in s) for s in nuc] for a in BACKBONE},
                    "cov_any": sum(1 for s in nuc if s) / n,
                    "cov_ge3": sum(1 for s in nuc if len(s) >= 3) / n,
                    "cov_full6": sum(1 for s in nuc if len(s) == 6) / n,
                    "cov_talos": sum(pred) / n,
                    "talos_mask": [int(x) for x in pred],
                })
        # best record = highest TALOS-predictable coverage, tie-break on cov_ge3
        recs.sort(key=lambda r: (r["cov_talos"], r["cov_ge3"], r["ident"]), reverse=True)
        out[pdb] = {
            "n": n, "seq": seq, "fold": t["fold"], "fail18": pdb in fail18,
            "n_candidates": len(cands), "candidates": cands,
            "records": recs, "best": recs[0] if recs else None,
        }
        b = recs[0] if recs else None
        if b is None:
            tail = "-- no usable deposition"
        else:
            tail = ("bmrb {:>6s} cov_any {:.2f} ge3 {:.2f} full6 {:.2f} talos {:.2f}"
                    .format(b["bmrb"], b["cov_any"], b["cov_ge3"], b["cov_full6"],
                            b["cov_talos"]))
        print("[{:3d}/126] {} n={:2d} cand={:2d} {}".format(k + 1, pdb, n, len(cands), tail),
              flush=True)

    with open(os.path.join(RESULTS, "shift_availability.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote s14/results/shift_availability.json")


if __name__ == "__main__":
    main()
