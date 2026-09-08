"""COORDINATOR AUDIT 2 -- can a genuinely fresh benchmark exist?

The record's three instruments (tuning126, dev24, benchmark60) between them spend ALL 204
identity clusters of 9-16-residue peptides in the frozen 787-peptide database.  There is
zero cluster-disjoint supply left inside the corpus.  So a fresh benchmark needs data the
corpus does not contain.

Two supplies were measured against RCSB (2026-09-05):

    date-fresh   9-16mer single-protein-entity entries released after 2026-01-01:  28
                 ... after 2026-06-01:                                             20
    containment-fresh   9-16mer single-protein-entity, no nucleic acid, ABSENT from
                 pdbs/ and pdbs_ext/ entirely:                                    146

The date-fresh supply cannot reach the pre-registered minimum of ~40 targets and is
reported as inadequate.  The containment-fresh supply is what this module builds and
freezes.  Its provenance claim is NOT "deposited after the models were trained" -- it is
the stronger and more directly relevant one:

    these structures were never in `pdbs/` or `pdbs_ext/`, so they were never parsed into
    `peptide_db`, never clustered, never assigned a fold, never in any distogram's training
    set, never in the fragment library, never in any retrieval pool, and never scored by
    any arm of this project.

Every target is additionally required to be below IDENTITY_THRESHOLD to every training
peptide and every library fragment, so no homolog of a benchmark target trained the model
that will score it.

DOWNLOADS GO TO `s12/newpdbs/` AND NOWHERE ELSE.  Adding one file to `pdbs/` or `pdbs_ext/`
re-permutes `_scan`, which re-draws dev_set(24), which changes the tuning126 manifest --
the most dangerous operation in the repository.  This module never writes outside `s12/`.

    python -m s12.coord_benchN fetch     # download the 146 candidates to s12/newpdbs/
    python -m s12.coord_benchN build     # gate, dedup, screen, cluster -> the frozen manifest
"""
import glob
import json
import os
import sys
import time
import urllib.request

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
import peptide_db as pdb                   # noqa: E402
import fragment_db as fdb                  # noqa: E402
import protein_geometry as geo             # noqa: E402

NEW = os.path.join(ROOT, "s12", "newpdbs")
SUPPLY = os.path.join(ROOT, "s12", "results", "coord_rcsb_supply.json")
MANIFEST = os.path.join(ROOT, "s12", "results", "coord_benchN_manifest.json")

MIN_LEN, MAX_LEN = 9, 16
REBUILD_TOL = 1.5                          # peptide_db.REBUILD_TOL
IDENTITY_THRESHOLD = 0.6                   # peptide_db.IDENTITY_THRESHOLD


def fetch():
    os.makedirs(NEW, exist_ok=True)
    want = json.load(open(SUPPLY))["missing_9_16"]
    print(f"{len(want)} candidate entries", flush=True)
    ok = 0
    for k, pid in enumerate(want):
        dst = os.path.join(NEW, f"{pid}.pdb")
        if os.path.exists(dst):
            ok += 1
            continue
        try:
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb",
                                        timeout=45) as fh:
                blob = fh.read()
            with open(dst, "wb") as out:
                out.write(blob)
            ok += 1
        except Exception as exc:                                       # noqa: BLE001
            print(f"  {pid}: {type(exc).__name__} {exc}", flush=True)
        if k % 25 == 0:
            print(f"  {k}/{len(want)} ({ok} on disk)", flush=True)
        time.sleep(0.1)
    print(f"{ok}/{len(want)} downloaded to {NEW}")


def _header(path):
    """EXPDTA / NUMMDL / TITLE / KEYWDS, for the experiment-class screen."""
    out = {"expdta": "", "nummdl": None, "title": "", "keywds": "", "resolution": None}
    with open(path, errors="ignore") as fh:
        for line in fh:
            rec = line[:6]
            if rec == "ATOM  " or rec == "MODEL ":
                break
            if rec == "EXPDTA":
                out["expdta"] += line[10:].strip() + " "
            elif rec == "NUMMDL":
                try:
                    out["nummdl"] = int(line[10:].strip())
                except ValueError:
                    pass
            elif rec == "TITLE ":
                out["title"] += line[10:].strip() + " "
            elif rec == "KEYWDS":
                out["keywds"] += line[10:].strip() + " "
            elif line.startswith("REMARK   2 RESOLUTION"):
                for tok in line.split():
                    try:
                        out["resolution"] = float(tok)
                        break
                    except ValueError:
                        continue
    for k in ("expdta", "title", "keywds"):
        out[k] = " ".join(out[k].split())
    return out


def build():
    paths = sorted(glob.glob(os.path.join(NEW, "*.pdb")))
    print(f"{len(paths)} downloaded files", flush=True)

    existing = list(pdb.load())
    ex_seqs = {p.seq for p in existing}
    frags = list(fdb.load(False))
    frag_seqs = [q.seq for q in frags]
    print(f"corpus: {len(existing)} peptides, {len(frags)} fragments", flush=True)

    gated, rejected = [], []
    for path in paths:
        pid = os.path.basename(path)[:-4].upper()
        head = _header(path)
        try:
            seq, coords, phi, psi = geo.native_coords_from_pdb(path)
        except Exception as exc:                                       # noqa: BLE001
            rejected.append({"pdb": pid, "why": f"parse: {type(exc).__name__}"}); continue
        n = len(seq)
        if not (MIN_LEN <= n <= MAX_LEN):
            rejected.append({"pdb": pid, "why": f"length {n}"}); continue
        if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(psi))):
            rejected.append({"pdb": pid, "why": "non-finite torsions"}); continue
        ca = np.asarray(coords["CA"], float)
        if len(ca) != n:
            rejected.append({"pdb": pid, "why": "CA count"}); continue
        step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        if step.min() < 3.5 or step.max() > 4.1:
            rejected.append({"pdb": pid, "why": f"CA step {step.min():.2f}-{step.max():.2f}"}); continue
        try:
            built = geo.build_backbone(phi, psi)
            keys = [k for k in ("N", "CA", "C") if k in coords and k in built]
            mob = np.vstack([built[k] for k in keys]); tgt = np.vstack([coords[k] for k in keys])
            rb = float(geo.rmsd(geo.kabsch_superpose(mob, tgt), tgt))
        except Exception as exc:                                       # noqa: BLE001
            rejected.append({"pdb": pid, "why": f"rebuild: {type(exc).__name__}"}); continue
        if not np.isfinite(rb) or rb > REBUILD_TOL:
            rejected.append({"pdb": pid, "why": f"rebuild {rb:.2f}"}); continue
        if "X" in seq:
            rejected.append({"pdb": pid, "why": "non-standard residue"}); continue
        if seq in ex_seqs:
            rejected.append({"pdb": pid, "why": "sequence already in corpus"}); continue
        gated.append({"pdb": pid, "seq": seq, "n": n, "rebuild": round(rb, 3),
                      "header": head, "ca": ca, "phi": phi, "psi": psi, "path": path})
        print(f"  KEEP {pid} n={n} rb={rb:.2f} {head['expdta'][:20]}", flush=True)

    print(f"\n{len(gated)} pass the structural gates, {len(rejected)} rejected", flush=True)

    # ---- deduplicate among themselves by exact sequence (best rebuild wins)
    best = {}
    for g in gated:
        if g["seq"] not in best or g["rebuild"] < best[g["seq"]]["rebuild"]:
            best[g["seq"]] = g
    uniq = sorted(best.values(), key=lambda r: r["pdb"])
    print(f"{len(uniq)} unique sequences", flush=True)

    # ---- leakage screen: identity AND containment to every training peptide and fragment
    kept, leaked = [], []
    from s12.coord_contain import containment
    for g in uniq:
        s = g["seq"]
        tk = {s[i:i + 3] for i in range(len(s) - 2)}
        max_id_pep = max_ct_pep = 0.0; who_pep = None
        for q in existing:
            if not (tk & {q.seq[i:i + 3] for i in range(len(q.seq) - 2)}):
                continue
            idv = pdb.identity(s, q.seq)
            if idv > max_id_pep:
                max_id_pep, who_pep = idv, q.pdb
            max_ct_pep = max(max_ct_pep, containment(s, q.seq))
        max_id_fr = max_ct_fr = 0.0
        for fs in frag_seqs:
            if not (tk & {fs[i:i + 3] for i in range(len(fs) - 2)}):
                continue
            max_id_fr = max(max_id_fr, pdb.identity(s, fs))
            max_ct_fr = max(max_ct_fr, containment(s, fs))
        verb_pep = any(s in q.seq for q in existing)
        verb_fr = any(s in fs for fs in frag_seqs)
        g["screen"] = {"max_identity_peptide": round(max_id_pep, 4),
                       "max_containment_peptide": round(max_ct_pep, 4),
                       "nearest_peptide": who_pep,
                       "max_identity_fragment": round(max_id_fr, 4),
                       "max_containment_fragment": round(max_ct_fr, 4),
                       "verbatim_in_peptides": bool(verb_pep),
                       "verbatim_in_fragments": bool(verb_fr)}
        #: THE SCREEN, AND WHY IT IS NOT A CONTAINMENT THRESHOLD.
        #: A containment cut against the 6,003-fragment bank is at the NULL at this
        #: length: 12 randomly generated sequences per length (composition-matched to the
        #: corpus) score mean max-containment 0.56-0.63 against the fragments and up to
        #: 0.78, with identity at 0.41-0.45 (s12/results/coord_containment_null.json).
        #: So `containment >= 0.6` rejects a random sequence about half the time and is
        #: evidence of nothing.  The three rules below each sit clear of that null:
        #:   * verbatim substring   -- p ~ 0 under the null at n >= 9
        #:   * identity >= 0.6      -- the production rule, ~0.15 above the null mean
        #:   * containment >= 0.85  -- ~0.25 above the null mean, so a real near-copy
        reasons = []
        if verb_pep or verb_fr:
            reasons.append("verbatim copy in the corpus")
        if max(max_id_pep, max_id_fr) >= IDENTITY_THRESHOLD:
            reasons.append(f"identity {max(max_id_pep, max_id_fr):.2f} >= {IDENTITY_THRESHOLD}")
        if max(max_ct_pep, max_ct_fr) >= 0.85:
            reasons.append(f"containment {max(max_ct_pep, max_ct_fr):.2f} >= 0.85")
        g["screen"]["reject_reasons"] = reasons
        (leaked if reasons else kept).append(g)
        print(f"  {g['pdb']} id_pep={max_id_pep:.2f} ct_pep={max_ct_pep:.2f} "
              f"id_frag={max_id_fr:.2f} ct_frag={max_ct_fr:.2f} "
              f"-> {('LEAK: ' + '; '.join(reasons)) if reasons else 'clean'}", flush=True)

    # ---- cluster the survivors among THEMSELVES, keep one per cluster
    reps, seen = [], []
    for g in sorted(kept, key=lambda r: r["pdb"]):
        if any(pdb.identity(g["seq"], h["seq"]) >= IDENTITY_THRESHOLD
               or g["seq"] in h["seq"] or h["seq"] in g["seq"] for h in seen):
            continue
        seen.append(g); reps.append(g)

    #: The EXPERIMENT-CLASS composition, reported rather than filtered on.  The 9-16mer
    #: single-chain supply that this corpus never fetched is dominated by cryo-EM amyloid
    #: fibril segments, whose deposited conformation is held by inter-chain packing -- the
    #: class `results/monomer_manifest.json` already excludes as unfoldable in isolation.
    import collections
    klass = collections.Counter()
    for g in reps:
        e = g["header"]["expdta"].upper(); t = (g["header"]["title"] + " " +
                                                g["header"]["keywds"]).upper()
        if "FIBRIL" in t or "AMYLOID" in t:
            klass["fibril/amyloid"] += 1
        elif "MICELLE" in t or "MEMBRANE" in t or "BICELLE" in t:
            klass["membrane/micelle"] += 1
        elif "NMR" in e:
            klass["solution NMR, other"] += 1
        elif "ELECTRON MICROSCOPY" in e:
            klass["cryo-EM, other"] += 1
        else:
            klass["X-ray/other"] += 1

    man = {
        "name": "benchN",
        "built": "2026-09-05",
        "provenance": "RCSB single-protein-entity, no nucleic acid, 9-16 residues, ABSENT "
                      "from pdbs/ and pdbs_ext/ at the time of the sprint-12 audit. Never "
                      "parsed into peptide_db, never clustered, never in a fold, never in "
                      "any distogram training set, never in the fragment library, never in "
                      "a retrieval pool, never scored by any arm of this project.",
        "gates": {"length": [MIN_LEN, MAX_LEN], "rebuild_tol": REBUILD_TOL,
                  "ca_step": [3.5, 4.1], "finite_torsions": True,
                  "dedup": "exact sequence, best rebuild",
                  "leakage": f"rejected if the sequence appears VERBATIM in any of the "
                             f"{len(existing)} corpus peptides or {len(frags)} library "
                             f"fragments, OR identity >= {IDENTITY_THRESHOLD} to any of "
                             f"them, OR containment >= 0.85 -- three rules each clear of "
                             f"the measured null",
                  "self_cluster": f"one representative per identity>={IDENTITY_THRESHOLD} "
                                  f"or substring cluster among the survivors"},
        "counts": {"downloaded": len(paths), "passed_gates": len(gated),
                   "unique_sequences": len(uniq), "clean_of_leakage": len(kept),
                   "leaked_and_dropped": len(leaked), "final": len(reps)},
        "experiment_classes": dict(klass),
        "null_for_the_screen": "s12/results/coord_containment_null.json -- random "
                               "composition-matched sequences score max-containment "
                               "0.56-0.63 against the fragment bank, so a containment "
                               "threshold at 0.6 is at the null and was not used",
        "reference": "deposited coordinates, model 1 (the project's declared endpoint)",
        "adequately_powered": len(reps) >= 40,
        "targets": [{"pdb": g["pdb"], "seq": g["seq"], "n": g["n"], "rebuild": g["rebuild"],
                     "expdta": g["header"]["expdta"], "nummdl": g["header"]["nummdl"],
                     "resolution": g["header"]["resolution"], "title": g["header"]["title"],
                     "screen": g["screen"]} for g in reps],
        "dropped_for_leakage": [{"pdb": g["pdb"], "seq": g["seq"], "screen": g["screen"]}
                                for g in leaked],
        "rejected": rejected,
    }
    with open(MANIFEST, "w") as fh:
        json.dump(man, fh, indent=1)
    print(f"\nbenchN: {len(reps)} targets  (adequately powered: {man['adequately_powered']})")
    print(f"wrote {MANIFEST}")
    return man


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "build"
    if cmd == "fetch":
        fetch()
    else:
        build()
