"""SHIFT agent -- build a TRAINING corpus of (backbone chemical shifts -> phi/psi) pairs
from BMRB depositions with matched PDB coordinates.

This is the training set for a TALOS-class shift->torsion model.  It is built entirely
OUTSIDE the 126 development targets, and every candidate protein is rejected if it
contains any target's sequence as a substring or aligns to one above the project's identity
threshold -- the containment failure class Sprint 13 documented in `s13/tors_common.py`.

No development-target structure is read here.  Coordinates are read only for corpus
proteins, which are not in the instrument.

    python -m s14.shift_corpus            # incremental; safe to re-run
"""
from __future__ import annotations

import os
import sys
import json
import gzip
import time
import random
import urllib.request

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s14 import shift_bmrb as B                          # noqa: E402

PDBCACHE = os.path.join(ROOT, "s14", "cache", "pdb")
os.makedirs(PDBCACHE, exist_ok=True)
OUT = os.path.join(B.RESULTS, "shift_corpus.npz")

MIN_LEN, MAX_LEN = 30, 220
TARGET_N = 400          # proteins to accept
D2R = np.pi / 180.0


# ------------------------------------------------------------------------ coordinates
def fetch_pdb(pdb):
    p = os.path.join(PDBCACHE, pdb.upper() + ".pdb.gz")
    if not os.path.exists(p):
        url = f"https://files.rcsb.org/download/{pdb.upper()}.pdb.gz"
        req = urllib.request.Request(url, headers={"User-Agent": "s14-shift-agent"})
        with urllib.request.urlopen(req, timeout=120) as fh:
            raw = fh.read()
        tmp = p + ".tmp"
        with open(tmp, "wb") as fh:
            fh.write(raw)
        os.replace(tmp, p)
    with gzip.open(p, "rt", errors="ignore") as fh:
        return fh.read()


def backbone_from_pdb(text, want_chain=None):
    """Model 1 only.  -> {chain: {resseq: {'aa':X, 'N':xyz, 'CA':xyz, 'C':xyz}}}"""
    out = {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name not in ("N", "CA", "C"):
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        res = line[17:20].strip().upper()
        aa = B.AA3.get(res)
        if aa is None:
            continue
        ch = line[21]
        if want_chain is not None and ch != want_chain:
            continue
        try:
            rs = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        out.setdefault(ch, {}).setdefault(rs, {"aa": aa})[name] = xyz
    return out


def _dihedral(p0, p1, p2, p3):
    b0 = np.asarray(p0) - np.asarray(p1)
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    b1n = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.arctan2(y, x))


def torsions(chain):
    """chain: {resseq: {...}} -> {resseq: (phi, psi, aa)} for residues with both defined."""
    keys = sorted(chain)
    out = {}
    for i, r in enumerate(keys):
        c = chain[r]
        if not all(k in c for k in ("N", "CA", "C")):
            continue
        prev = chain.get(r - 1)
        nxt = chain.get(r + 1)
        if prev is None or nxt is None:
            continue
        if "C" not in prev or "N" not in nxt:
            continue
        # require genuine peptide bonds on both sides (numbering can lie)
        d1 = np.linalg.norm(np.array(prev["C"]) - np.array(c["N"]))
        d2 = np.linalg.norm(np.array(c["C"]) - np.array(nxt["N"]))
        if d1 > 2.0 or d2 > 2.0:
            continue
        phi = _dihedral(prev["C"], c["N"], c["CA"], c["C"])
        psi = _dihedral(c["N"], c["CA"], c["C"], nxt["N"])
        out[r] = (phi, psi, c["aa"])
    return out


# ------------------------------------------------------------------------- assembly
def seq_of(chain):
    keys = sorted(chain)
    return "".join(chain[k]["aa"] for k in keys), keys


def build(target_seqs, limit=TARGET_N, seed=0):
    ex = B.pdb_bmrb_map("exact")
    items = sorted(ex.items())
    random.Random(seed).shuffle(items)

    rows = []
    accepted, seen_seq = 0, set()
    t0 = time.time()
    for k, (pdb, bids) in enumerate(items):
        if accepted >= limit:
            break
        bid = str(bids[0])
        try:
            ents = B.entry_entities(bid)
        except Exception:                                          # noqa: BLE001
            continue
        ents = [e for e in ents
                if e["type"] and "polypeptide" in e["type"]
                and MIN_LEN <= len(e["seq"]) <= MAX_LEN]
        if not ents:
            continue
        # containment / identity guard against every development target
        bad = False
        for e in ents:
            s = e["seq"].upper()
            for ts in target_seqs:
                if ts in s or s in ts:
                    bad = True
                    break
            if bad:
                break
        if bad:
            continue
        try:
            shifts = B.entry_shifts(bid)
        except Exception:                                          # noqa: BLE001
            continue
        if not shifts:
            continue
        try:
            text = fetch_pdb(pdb)
        except Exception:                                          # noqa: BLE001
            continue
        chains = backbone_from_pdb(text)
        if not chains:
            continue

        got = 0
        for e in ents:
            eseq = e["seq"].upper()
            if eseq in seen_seq:
                continue
            # find the PDB chain whose sequence best contains the entity sequence
            best = None
            for ch, cd in chains.items():
                cseq, keys = seq_of(cd)
                j = cseq.find(eseq[:min(len(eseq), 30)])
                if j < 0:
                    continue
                if best is None or len(cseq) > len(best[1]):
                    best = (ch, cseq, keys, cd, j)
            if best is None:
                continue
            ch, cseq, keys, cd, j = best
            tors = torsions(cd)
            # per-residue shifts for this entity
            have = {}
            for r in shifts:
                if str(r["entity"]) != str(e["id"]):
                    continue
                a = B.ATOM_ALIAS.get(r["atom"])
                if a is None:
                    continue
                have.setdefault(r["seq_id"], {})[a] = r["val"]
            if len(have) < 10:
                continue
            for i in range(len(eseq)):
                sid = i + 1                       # entity Seq_ID is 1-based
                pk = keys[j + i] if 0 <= j + i < len(keys) else None
                if pk is None or pk not in tors:
                    continue
                phi, psi, aa = tors[pk]
                if aa != eseq[i]:
                    continue
                sh = have.get(sid)
                if not sh or len(sh) < 3:
                    continue
                rows.append({
                    "pdb": pdb, "bmrb": bid, "ent": e["id"], "aa": aa,
                    "prev": eseq[i - 1] if i > 0 else "-",
                    "next": eseq[i + 1] if i + 1 < len(eseq) else "-",
                    "phi": phi, "psi": psi,
                    "sh": sh,
                    "sh_prev": have.get(sid - 1, {}),
                    "sh_next": have.get(sid + 1, {}),
                })
                got += 1
            if got:
                seen_seq.add(eseq)
        if got:
            accepted += 1
            if accepted % 20 == 0:
                print("  accepted {:4d} proteins, {:6d} residues  [{:.0f}s, scanned {}]"
                      .format(accepted, len(rows), time.time() - t0, k + 1), flush=True)
    return rows


def main():
    from s12 import instrument as I
    tg = I.targets()
    tseqs = [t["seq"].upper() for t in tg]
    rows = build(tseqs)
    print("corpus: {} proteins-worth, {} residues".format(
        len({r["pdb"] for r in rows}), len(rows)))
    with open(os.path.join(B.RESULTS, "shift_corpus.json"), "w") as fh:
        json.dump(rows, fh)
    print("wrote s14/results/shift_corpus.json")


if __name__ == "__main__":
    main()
