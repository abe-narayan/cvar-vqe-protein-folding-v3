"""s26/e_hashes.py -- sha256 and size of every pinned artefact, so drift can be detected.

Lane E (Examiner), Sprint 26.  Hashes BYTES only.  It never parses, prints or otherwise
reads the contents of ``results/benchmark_manifest.json`` (Rule 1); the sealed benchmark's
hash is compared against the S20 record (a40581ad... 8002 bytes) and the verdict is printed
without opening the file as text.

    python s26/e_hashes.py            # write s26/results/pinned_hashes.json
    python s26/e_hashes.py --check    # re-hash and report any drift against that file

Nothing here loads ``esm_cache.npz`` (1.5 GB): it is not in the pinned list and its size is
recorded from ``os.stat`` only.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "s26", "results", "pinned_hashes.json")

PINNED_FILES = [
    "peptide_folds.json",
    "peptide_clusters.json",
    "catrace_prior.npz",
    "results/benchmark_manifest.json",
    "results/monomer_manifest.json",
    "esm_small.npz",
    "esm_pca.npz",
    "peptide_db.npz",
    "fragment_db.npz",
    "fragment_db_large.npz",
]
PINNED_GLOBS = [
    "distogram_models/*",
    "pdbs/*",
]
#: recorded in S20; the sealed benchmark manifest must still hash to this
S20_BENCHMARK_SHA256 = "a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d"
S20_BENCHMARK_BYTES = 8002
#: stat-only (never hashed, never opened): too large to hash casually
STAT_ONLY = ["esm_cache.npz"]


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def collect():
    entries = {}
    paths = list(PINNED_FILES)
    for g in PINNED_GLOBS:
        paths.extend(sorted(p.replace(os.sep, "/") for p in
                            glob.glob(os.path.join(ROOT, g))
                            if os.path.isfile(p)))
    for rel in paths:
        rel = os.path.relpath(rel, ROOT).replace(os.sep, "/") if os.path.isabs(rel) else rel
        ab = os.path.join(ROOT, rel)
        if not os.path.isfile(ab):
            entries[rel] = {"absent": True}
            continue
        st = os.stat(ab)
        entries[rel] = {"sha256": sha256_of(ab), "bytes": int(st.st_size),
                        "mtime": time.strftime("%Y-%m-%dT%H:%M:%S",
                                               time.localtime(st.st_mtime))}
    stat_only = {}
    for rel in STAT_ONLY:
        ab = os.path.join(ROOT, rel)
        stat_only[rel] = ({"bytes": int(os.stat(ab).st_size), "hashed": False}
                          if os.path.isfile(ab) else {"absent": True})
    return entries, stat_only


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="compare against the stored file and report drift")
    a = ap.parse_args(argv)
    entries, stat_only = collect()
    bm = entries.get("results/benchmark_manifest.json", {})
    bm_ok = (bm.get("sha256") == S20_BENCHMARK_SHA256
             and bm.get("bytes") == S20_BENCHMARK_BYTES)
    n_pdb = sum(1 for k in entries if k.startswith("pdbs/") and "sha256" in entries[k])
    n_dm = sum(1 for k in entries if k.startswith("distogram_models/")
               and "sha256" in entries[k])
    payload = {
        "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "root": ROOT,
        "benchmark_manifest_matches_S20_record": bool(bm_ok),
        "S20_record": {"sha256": S20_BENCHMARK_SHA256, "bytes": S20_BENCHMARK_BYTES},
        "n_pdbs_files": n_pdb,
        "n_distogram_model_files": n_dm,
        "stat_only": stat_only,
        "files": entries,
    }
    if a.check and os.path.exists(OUT):
        with open(OUT) as fh:
            old = json.load(fh)["files"]
        drift = []
        for k, v in entries.items():
            o = old.get(k)
            if o is None:
                drift.append((k, "NEW"))
            elif o.get("sha256") != v.get("sha256") or o.get("bytes") != v.get("bytes"):
                drift.append((k, "CHANGED"))
        for k in old:
            if k not in entries:
                drift.append((k, "MISSING"))
        print(f"benchmark manifest matches S20 record: {bm_ok}")
        if drift:
            print("DRIFT DETECTED:")
            for k, w in drift:
                print(f"  {w:8} {k}")
            return 1
        print(f"no drift: {len(entries)} entries identical to {OUT}")
        return 0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {OUT}")
    print(f"benchmark manifest matches S20 record (sha256 {S20_BENCHMARK_SHA256[:12]}..., "
          f"{S20_BENCHMARK_BYTES} bytes): {bm_ok}")
    print(f"pdbs/ files hashed: {n_pdb}   distogram_models/ files hashed: {n_dm}")
    for k, v in entries.items():
        if k.startswith("pdbs/"):
            continue
        if "sha256" in v:
            print(f"{v['sha256']}  {v['bytes']:>12d}  {k}")
        else:
            print(f"{'ABSENT':64}  {'':>12}  {k}")
    for k, v in stat_only.items():
        print(f"{'(stat only, not hashed)':64}  {v.get('bytes', 'ABSENT'):>12}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
